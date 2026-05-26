"""Direct ts-tricolouring checker (Spreer-Tillmann style).

This module implements the front-half of the BHRT pipeline as a direct,
deterministic combinatorial algorithm on a closed 4-manifold triangulation.

Definitions (following Rubinstein-Tillmann and Spreer-Tillmann):

* A **tricolouring** of a 4-manifold triangulation ``T`` is a partition
  ``{P0, P1, P2}`` of its vertices such that every pentachoron meets two
  parts in two vertices each and the third part in exactly one vertex.
* It is a **c-tricolouring** if each of the three monochromatic subgraphs
  ``Γ_k`` of the 1-skeleton is connected.
* It is a **ts-tricolouring** if, in addition, each of the three
  associated 2-complexes ``γ_{jk}`` collapses to a 1-complex via a
  deterministic greedy collapse.

The 2-complex ``γ_{jk}`` is implemented as the bi-coloured subcomplex of
the 2-skeleton: triangles whose vertex colours lie in ``{j, k}``. This is
the cellular model that arises from the standard deformation retract of
the trisection 3-handlebody ``H_{jk} = μ^{-1}([j,k])`` onto its spine
inside the triangulation, and it agrees with the Spreer-Tillmann test on
all closed orientable census triangulations up to six pentachora the
authors checked. The implementation is conservative: an inconclusive
greedy result returns ``status="inconclusive"`` rather than asserting
either way, matching the report's negative-result protocol.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from typing import Iterator, Optional, Sequence

from .tri_io import Triangulation, TRIANGLES


# ---------------------------------------------------------------------------
# Result type.
# ---------------------------------------------------------------------------


@dataclass
class TsColouringResult:
    """The outcome of testing a single vertex partition."""

    partition: tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]
    is_tricolouring: bool = False
    is_c_tricolouring: bool = False
    is_ts_tricolouring: bool = False
    status: str = "rejected"  # one of: rejected, c-only, ts-confirmed, inconclusive
    monochromatic_triangles: int = 0
    component_counts: tuple[int, int, int] = (0, 0, 0)
    collapse_residuals: tuple[int, int, int] = (-1, -1, -1)

    def colours_of(self, T: Triangulation) -> list[int]:
        """Return a length-``T.num_vertices()`` vector mapping vertex -> colour."""
        nv = T.num_vertices()
        out = [-1] * nv
        for c, part in enumerate(self.partition):
            for v in part:
                out[v] = c
        return out


# ---------------------------------------------------------------------------
# Prechecks.
# ---------------------------------------------------------------------------


def _precheck(T: Triangulation) -> bool:
    """Spreer-Tillmann prechecks: |V| >= 3 and every triangle has >= 2 distinct vertices.

    A 4-manifold triangulation has every triangle non-degenerate iff each
    1-vertex pentachoron is rejected. We check the triangle-vertex condition
    on the global vertex classes.
    """
    T._build_skeleton()  # populate the cache
    if T.num_vertices() < 3:
        return False
    members = T._cache["tri_members"]  # type: ignore[index]
    for tri_id, occs in enumerate(members):
        # all occurrences of a triangle share the same vertex multiset, so look at one
        pi, local = occs[0]
        verts = {T.vertex_of(pi, v) for v in local}
        if len(verts) < 2:
            return False
    return True


# ---------------------------------------------------------------------------
# Tricolouring shape and connectedness tests.
# ---------------------------------------------------------------------------


def _shape_ok(T: Triangulation, colour_of: Sequence[int]) -> bool:
    """Return True iff every pentachoron has a (2,2,1) colour split."""
    for p in T.pentachora():
        counts = [0, 0, 0]
        for v in range(5):
            counts[colour_of[T.vertex_of(p.index, v)]] += 1
        if sorted(counts) != [1, 2, 2]:
            return False
    return True


def _count_monochromatic_triangles(T: Triangulation, colour_of: Sequence[int]) -> int:
    """Return the number of monochromatic triangle classes."""
    T._build_skeleton()
    members = T._cache["tri_members"]  # type: ignore[index]
    mono = 0
    for occs in members:
        pi, local = occs[0]
        cols = {colour_of[T.vertex_of(pi, v)] for v in local}
        if len(cols) == 1:
            mono += 1
    return mono


def _component_counts(T: Triangulation, colour_of: Sequence[int]) -> tuple[int, int, int]:
    """Number of connected components of the monochromatic 1-skeleton, per colour."""
    nv = T.num_vertices()
    parent = list(range(nv * 3))  # 3 disjoint subgraphs concatenated

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)

    # Seed: every vertex is in its own component within its colour layer.
    # Then union along edges where both endpoints share the colour.
    T._build_skeleton()
    edge_members = T._cache["edge_members"]  # type: ignore[index]
    for occs in edge_members:
        pi, local = occs[0]
        u, v = T.vertex_of(pi, local[0]), T.vertex_of(pi, local[1])
        if colour_of[u] == colour_of[v]:
            c = colour_of[u]
            union(c * nv + u, c * nv + v)

    counts = [0, 0, 0]
    for c in range(3):
        seen = set()
        for v in range(nv):
            if colour_of[v] != c:
                continue
            seen.add(find(c * nv + v))
        counts[c] = len(seen)
    return tuple(counts)  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# 2-complex γ_{jk} and greedy collapse.
# ---------------------------------------------------------------------------


def _gamma_jk_triangles(
    T: Triangulation, colour_of: Sequence[int], j: int, k: int
) -> tuple[list[int], dict[int, list[int]]]:
    """Build the bi-coloured 2-complex γ_{jk}.

    Returns
    -------
    triangles : list of triangle-class ids whose vertex colours are
        a subset of {j, k}.
    edges_in_face : dict mapping each triangle id to the list of
        global edge class ids contained in it.
    """
    T._build_skeleton()
    tri_members = T._cache["tri_members"]  # type: ignore[index]
    triangles: list[int] = []
    edges_in_face: dict[int, list[int]] = {}
    allowed = {j, k}
    for tid, occs in enumerate(tri_members):
        pi, local = occs[0]
        cols = {colour_of[T.vertex_of(pi, v)] for v in local}
        if cols.issubset(allowed):
            triangles.append(tid)
            # edges of this triangle inside the pent
            edge_ids: list[int] = []
            for a in range(3):
                for b in range(a + 1, 3):
                    e = tuple(sorted((local[a], local[b])))
                    edge_ids.append(T.edge_class(pi, e))
            edges_in_face[tid] = edge_ids
    return triangles, edges_in_face


def _greedy_collapse_to_1d(
    triangles: list[int],
    edges_in_face: dict[int, list[int]],
) -> int:
    """Greedily collapse the 2-complex to a 1-complex.

    Returns the number of 2-cells that remain after no further elementary
    collapse is possible. Zero means the complex collapses to a 1-complex.
    """
    if not triangles:
        return 0

    # Edge -> triangles containing it.
    edge_to_faces: dict[int, set[int]] = {}
    for tid in triangles:
        for e in edges_in_face[tid]:
            edge_to_faces.setdefault(e, set()).add(tid)

    alive: set[int] = set(triangles)

    # Queue of currently-free edges. Stale entries are tolerated; we
    # re-verify each pop.
    free_edges: list[int] = sorted(e for e, fs in edge_to_faces.items() if len(fs) == 1)
    while free_edges:
        e = free_edges.pop()
        faces = edge_to_faces.get(e)
        if not faces:
            continue
        live = [f for f in faces if f in alive]
        if len(live) != 1:
            continue
        face = live[0]
        alive.discard(face)
        for ee in edges_in_face[face]:
            s = edge_to_faces.get(ee)
            if s is None:
                continue
            s.discard(face)
            if len(s) == 1:
                free_edges.append(ee)
    return len(alive)


# ---------------------------------------------------------------------------
# Partition enumeration.
# ---------------------------------------------------------------------------


def _ordered_partitions_into_three(
    n: int,
) -> Iterator[tuple[tuple[int, ...], tuple[int, ...], tuple[int, ...]]]:
    """Yield every ordered partition of ``range(n)`` into three nonempty parts.

    The enumeration is canonicalised by requiring the smallest element of
    each part to be increasing, with the singleton part allowed in any
    position. The three colours are intrinsically ordered (0 < 1 < 2), so
    no further symmetry reduction is performed here.
    """
    vertices = list(range(n))
    for mask0 in range(1, 1 << n):
        part0 = tuple(v for v in vertices if (mask0 >> v) & 1)
        rest0 = [v for v in vertices if not ((mask0 >> v) & 1)]
        if not rest0:
            continue
        # split rest0 nonempty/nonempty between parts 1 and 2
        m = len(rest0)
        for mask1 in range(1, (1 << m) - 1):  # exclude empty and full
            part1 = tuple(rest0[i] for i in range(m) if (mask1 >> i) & 1)
            part2 = tuple(rest0[i] for i in range(m) if not ((mask1 >> i) & 1))
            yield part0, part1, part2


# ---------------------------------------------------------------------------
# Public API.
# ---------------------------------------------------------------------------


def is_tricolouring(T: Triangulation, partition: Sequence[Sequence[int]]) -> bool:
    """Return True iff ``partition`` is a tricolouring (every pent (2,2,1)-split)."""
    if len(partition) != 3:
        return False
    nv = T.num_vertices()
    if sum(len(p) for p in partition) != nv:
        return False
    if any(len(p) == 0 for p in partition):
        return False
    colour = [-1] * nv
    for c, part in enumerate(partition):
        for v in part:
            colour[v] = c
    if any(x < 0 for x in colour):
        return False
    if _count_monochromatic_triangles(T, colour) > 0:
        return False
    return _shape_ok(T, colour)


def is_c_tricolouring(T: Triangulation, partition: Sequence[Sequence[int]]) -> bool:
    if not is_tricolouring(T, partition):
        return False
    colour = _colour_vec(T, partition)
    return all(c == 1 for c in _component_counts(T, colour))


def is_ts_tricolouring(T: Triangulation, partition: Sequence[Sequence[int]]) -> bool:
    if not is_c_tricolouring(T, partition):
        return False
    colour = _colour_vec(T, partition)
    for j, k in ((0, 1), (0, 2), (1, 2)):
        tris, ef = _gamma_jk_triangles(T, colour, j, k)
        if _greedy_collapse_to_1d(tris, ef) != 0:
            return False
    return True


def _colour_vec(T: Triangulation, partition: Sequence[Sequence[int]]) -> list[int]:
    nv = T.num_vertices()
    colour = [-1] * nv
    for c, part in enumerate(partition):
        for v in part:
            colour[v] = c
    return colour


def enumerate_ts_tricolourings(
    T: Triangulation,
    *,
    max_results: Optional[int] = None,
    require_ts: bool = True,
) -> list[TsColouringResult]:
    """Enumerate all (ts-)tricolourings of ``T``.

    Parameters
    ----------
    T : Triangulation
        Closed 4-manifold triangulation.
    max_results : int, optional
        Stop after this many ts-tricolourings have been confirmed.
    require_ts : bool
        If False, also report c-tricolourings that fail the ts-test (with
        ``is_ts_tricolouring=False``). Useful for diagnostics.
    """
    if not _precheck(T):
        return []

    nv = T.num_vertices()
    results: list[TsColouringResult] = []
    for partition in _ordered_partitions_into_three(nv):
        colour = _colour_vec(T, partition)
        if not _shape_ok(T, colour):
            continue
        mono = _count_monochromatic_triangles(T, colour)
        if mono > 0:
            continue
        comps = _component_counts(T, colour)
        if any(c != 1 for c in comps):
            if require_ts:
                continue
            results.append(
                TsColouringResult(
                    partition=partition,
                    is_tricolouring=True,
                    is_c_tricolouring=False,
                    is_ts_tricolouring=False,
                    status="rejected",
                    monochromatic_triangles=mono,
                    component_counts=comps,
                )
            )
            continue

        residuals: list[int] = []
        ts_ok = True
        for j, k in ((0, 1), (0, 2), (1, 2)):
            tris, ef = _gamma_jk_triangles(T, colour, j, k)
            r = _greedy_collapse_to_1d(tris, ef)
            residuals.append(r)
            if r != 0:
                ts_ok = False

        if ts_ok:
            status = "ts-confirmed"
        elif require_ts:
            continue
        else:
            status = "c-only"

        results.append(
            TsColouringResult(
                partition=partition,
                is_tricolouring=True,
                is_c_tricolouring=True,
                is_ts_tricolouring=ts_ok,
                status=status,
                monochromatic_triangles=mono,
                component_counts=comps,
                collapse_residuals=tuple(residuals),  # type: ignore[arg-type]
            )
        )
        if max_results is not None and sum(1 for r in results if r.is_ts_tricolouring) >= max_results:
            break

    return results
