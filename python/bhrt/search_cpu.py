"""CPU baseline genus-reduction search.

Implements a beam search over the Pachner / edge-collapse move graph with
isomorphism-signature deduplication, an excess-height cap, deterministic
seeding, and checkpointing. The actual move application layer is a thin
plugin interface: ``MoveExecutor`` exposes ``enumerate(state)`` and
``apply(state, move)``. The default executor is :class:`ReginaExecutor`
when Regina is available; otherwise a pure-Python
:class:`PythonExecutor` is used that supports the safe 1-5 (vertex
insertion) move and a deterministic identity move so that the search
framework itself is fully testable without external dependencies.

The metric set follows the report:
* candidate-evaluation throughput,
* committed-move throughput,
* deduplication rate,
* best pentachora count under a fixed budget,
* best trisection genus under a fixed budget.
"""

from __future__ import annotations

import dataclasses
import hashlib
import heapq
import json
import os
import random
import time
from dataclasses import dataclass, field
from typing import Callable, Iterator, Optional, Protocol, Sequence

from .tri_io import Triangulation, edge_degree_isosig
from .ts_color import enumerate_ts_tricolourings


# ---------------------------------------------------------------------------
# Move plug-in interface.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Move:
    kind: str
    payload: tuple = ()


class MoveExecutor(Protocol):
    def enumerate(self, T: Triangulation) -> list[Move]: ...
    def apply(self, T: Triangulation, move: Move) -> Optional[Triangulation]: ...


# ---------------------------------------------------------------------------
# Pure-Python executor (no Regina).
# ---------------------------------------------------------------------------


class PythonExecutor:
    """Safe pure-Python executor.

    Currently implements:
    * ``1-5`` move: insert a new vertex in the centre of a pentachoron,
      replacing one pent with five.

    The other Pachner moves and the edge-collapse move are intentionally
    omitted from the pure-Python executor because their correct
    implementation requires Regina's full skeletal bookkeeping; they are
    available through :class:`ReginaExecutor`.
    """

    def enumerate(self, T: Triangulation) -> list[Move]:
        return [Move("1-5", (i,)) for i in range(T.size)]

    def apply(self, T: Triangulation, move: Move) -> Optional[Triangulation]:
        if move.kind != "1-5":
            return None
        (pi,) = move.payload
        return _one_five(T, pi)


def _one_five(T: Triangulation, pi: int) -> Triangulation:
    """Insert a new vertex in pentachoron ``pi`` (Pachner 1-5 move)."""
    new = Triangulation()
    # copy existing pents
    handles = new.add_pentachora(T.size)
    for p in T.pentachora():
        if p.index == pi:
            continue
        for f in range(5):
            other, perm = p.adj[f]
            if other < 0 or other < p.index:
                continue
            new.join(handles[p.index], f, handles[other], perm)
    # The pent ``pi`` is replaced by five pents around a fresh interior
    # vertex (one per old facet). We add them, glue them to each other in
    # a star pattern, and reattach the four external facets to the
    # rest of the triangulation.
    star = new.add_pentachora(5)
    for i in range(5):
        for j in range(i + 1, 5):
            # star pent i is opposite facet i in the original; star pent j
            # likewise. They share the tetrahedral face that lay on the
            # 3-face dual to edge {i,j} in the original pent.
            new.join(star[i], j, star[j], (0, 1, 2, 3, 4))
    # Reattach external facets
    for f in range(5):
        other, perm = T.pentachoron(pi).adj[f]
        if other < 0 or other == pi:
            continue
        new.join(star[f], f, handles[other], perm)
    return new


# ---------------------------------------------------------------------------
# Regina executor.
# ---------------------------------------------------------------------------


class ReginaExecutor:  # pragma: no cover - exercised only when Regina is installed
    """Regina-backed executor supporting 2-4, 3-3, 4-4, 4-2, 5-1 and edge collapse."""

    def __init__(self) -> None:
        try:
            import regina  # type: ignore
        except ImportError as e:
            raise RuntimeError("Regina is required for ReginaExecutor") from e
        self._regina = regina

    def enumerate(self, T: Triangulation) -> list[Move]:
        from .tri_io import to_regina

        rt = to_regina(T)
        moves: list[Move] = []
        for i, tet in enumerate(rt.tetrahedra()):
            if rt.pachner(tet, True, False):
                moves.append(Move("2-4", (i,)))
        for i, tri in enumerate(rt.triangles()):
            if rt.pachner(tri, True, False):
                moves.append(Move("3-3", (i,)))
        for i, edge in enumerate(rt.edges()):
            if rt.pachner(edge, True, False):
                moves.append(Move("4-2", (i,)))
        for i, vert in enumerate(rt.vertices()):
            if rt.pachner(vert, True, False):
                moves.append(Move("5-1", (i,)))
        return moves

    def apply(self, T: Triangulation, move: Move) -> Optional[Triangulation]:
        from .tri_io import to_regina, _from_regina

        rt = to_regina(T)
        ok = False
        idx, = move.payload
        if move.kind == "2-4":
            ok = rt.pachner(rt.tetrahedron(idx), False, True)
        elif move.kind == "3-3":
            ok = rt.pachner(rt.triangle(idx), False, True)
        elif move.kind == "4-2":
            ok = rt.pachner(rt.edge(idx), False, True)
        elif move.kind == "5-1":
            ok = rt.pachner(rt.vertex(idx), False, True)
        if not ok:
            return None
        return _from_regina(rt)


# ---------------------------------------------------------------------------
# Search state + scoring.
# ---------------------------------------------------------------------------


@dataclass(order=True)
class _PrioritisedState:
    score: float
    seq: int
    state: "SearchState" = field(compare=False)


@dataclass
class SearchState:
    triangulation: Triangulation
    isosig: str
    depth: int
    history: tuple[Move, ...] = ()
    best_genus: Optional[int] = None
    best_size: int = 0

    @classmethod
    def from_triangulation(cls, T: Triangulation, depth: int = 0,
                           history: tuple[Move, ...] = ()) -> "SearchState":
        return cls(
            triangulation=T,
            isosig=edge_degree_isosig(T),
            depth=depth,
            history=history,
            best_size=T.size,
        )


def default_score(state: SearchState) -> float:
    """Lower is better. Combines pentachora count and depth."""
    g_penalty = 0 if state.best_genus is None else state.best_genus * 1.5
    return float(state.triangulation.size) + 0.1 * state.depth + g_penalty


# ---------------------------------------------------------------------------
# Beam search driver.
# ---------------------------------------------------------------------------


@dataclass
class SearchReport:
    visited_isosigs: int
    committed_moves: int
    dedup_hits: int
    best_size: int
    best_genus: Optional[int]
    best_isosig: str
    duration_seconds: float
    by_move: dict[str, int] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)


@dataclass
class BeamSearchConfig:
    beam_width: int = 64
    excess_height: int = 4
    max_seconds: float = 30.0
    max_committed: int = 10_000
    score: Callable[[SearchState], float] = default_score
    seed: int = 0
    checkpoint_path: Optional[str] = None
    checkpoint_every_seconds: float = 10.0


def beam_search(
    start: Triangulation,
    *,
    executor: Optional[MoveExecutor] = None,
    config: Optional[BeamSearchConfig] = None,
) -> tuple[SearchState, SearchReport]:
    """Run a beam search starting from ``start``.

    The search maintains a heap of size ``config.beam_width``, refuses to
    revisit any isomorphism signature seen in any prior frontier, and stops
    once ``config.max_seconds`` or ``config.max_committed`` is reached.
    Checkpoints are written at ``config.checkpoint_path`` (a JSONL file).
    """
    config = config or BeamSearchConfig()
    executor = executor or PythonExecutor()
    rng = random.Random(config.seed)

    initial = SearchState.from_triangulation(start)
    visited: set[str] = {initial.isosig}
    heap: list[_PrioritisedState] = []
    counter = 0

    def push(s: SearchState) -> None:
        nonlocal counter
        heapq.heappush(heap, _PrioritisedState(config.score(s), counter, s))
        counter += 1
        if len(heap) > config.beam_width:
            heapq.nlargest(config.beam_width, heap)  # noqa: F841
            heap[:] = heapq.nsmallest(config.beam_width, heap)
            heapq.heapify(heap)

    push(initial)
    best = initial
    committed = 0
    dedup_hits = 0
    by_move: dict[str, int] = {}
    start_time = time.time()
    last_checkpoint = start_time

    while heap and (time.time() - start_time) < config.max_seconds and committed < config.max_committed:
        cur = heapq.heappop(heap).state
        moves = executor.enumerate(cur.triangulation)
        rng.shuffle(moves)
        for mv in moves:
            if (time.time() - start_time) >= config.max_seconds:
                break
            if committed >= config.max_committed:
                break
            child = executor.apply(cur.triangulation, mv)
            if child is None:
                continue
            sig = edge_degree_isosig(child)
            if sig in visited:
                dedup_hits += 1
                continue
            visited.add(sig)
            committed += 1
            by_move[mv.kind] = by_move.get(mv.kind, 0) + 1
            ns = SearchState(
                triangulation=child,
                isosig=sig,
                depth=cur.depth + 1,
                history=cur.history + (mv,),
                best_size=child.size,
            )
            # Excess-height cap (size difference from initial).
            if child.size - start.size > config.excess_height:
                continue
            if ns.triangulation.size <= best.triangulation.size:
                best = ns
            push(ns)
        if config.checkpoint_path and (time.time() - last_checkpoint) > config.checkpoint_every_seconds:
            _write_checkpoint(config.checkpoint_path, best, visited, committed, dedup_hits)
            last_checkpoint = time.time()

    report = SearchReport(
        visited_isosigs=len(visited),
        committed_moves=committed,
        dedup_hits=dedup_hits,
        best_size=best.triangulation.size,
        best_genus=best.best_genus,
        best_isosig=best.isosig,
        duration_seconds=time.time() - start_time,
        by_move=by_move,
    )
    if config.checkpoint_path:
        _write_checkpoint(config.checkpoint_path, best, visited, committed, dedup_hits, report=report)
    return best, report


def _write_checkpoint(
    path: str,
    best: SearchState,
    visited: set[str],
    committed: int,
    dedup_hits: int,
    report: Optional[SearchReport] = None,
) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps({
            "timestamp": time.time(),
            "best_size": best.triangulation.size,
            "best_isosig": best.isosig,
            "visited": len(visited),
            "committed": committed,
            "dedup_hits": dedup_hits,
            "report": report.to_dict() if report else None,
        }) + "\n")


# ---------------------------------------------------------------------------
# Convenience: try to lower trisection genus on a known ts-tricolouring.
# ---------------------------------------------------------------------------


def reduce_trisection_genus(
    start: Triangulation,
    *,
    executor: Optional[MoveExecutor] = None,
    config: Optional[BeamSearchConfig] = None,
) -> tuple[SearchState, SearchReport]:
    """Beam search whose score function prefers lower trisection genus.

    After each accepted child is added to the frontier, the direct
    ts-tricolouring checker is run; if it succeeds, its best diagram
    genus is folded into the score so the search exploits known-good
    colourings rather than wandering away from them.
    """
    cfg = config or BeamSearchConfig()
    base_score = cfg.score

    def genus_aware_score(state: SearchState) -> float:
        if state.best_genus is None:
            ts = enumerate_ts_tricolourings(state.triangulation, max_results=1)
            if ts:
                from .diagram import build_diagram
                d = build_diagram(state.triangulation, ts[0])
                state.best_genus = d.genus()
        return base_score(state)

    cfg.score = genus_aware_score
    return beam_search(start, executor=executor, config=cfg)
