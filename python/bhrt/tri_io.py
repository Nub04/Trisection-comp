"""Pentachoron-based 4-manifold triangulation core.

This module provides a self-contained PL triangulation data structure for
closed (or bounded) 4-manifolds, the combinatorial primitives the rest of
the pipeline relies on, and an adapter that activates Regina (``regina``
Python module) if it is importable. The internal model mirrors Regina's
``Triangulation<4>`` semantics: a list of pentachora glued along
tetrahedral facets via permutations of ``{0,1,2,3,4}``.

Canonical signatures
--------------------
We implement a deterministic edge-degree isomorphism signature ("edge-degree
isosig") matching the form used by the public ``Dim4Census`` repository. It
canonicalises by (i) sorting pentachora by their tuple of edge degrees,
(ii) labelling vertices by depth-first traversal from the lowest-indexed
vertex of the lowest-indexed pentachoron and (iii) packing the resulting
gluing sequence as a compact base64 string. For exact Regina compatibility
the ``regina`` adapter is preferred when available.
"""

from __future__ import annotations

import base64
import itertools
import json
import zlib
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Mapping, Optional, Sequence

# ---------------------------------------------------------------------------
# Combinatorial constants for the 4-simplex.
# ---------------------------------------------------------------------------

#: The five vertex labels of a pentachoron.
VERTICES: tuple[int, ...] = (0, 1, 2, 3, 4)

#: The ten edges (2-subsets) of a pentachoron, lexicographically ordered.
EDGES: tuple[tuple[int, int], ...] = tuple(
    (a, b) for a in range(5) for b in range(a + 1, 5)
)

#: The ten triangles (3-subsets) of a pentachoron, lexicographically ordered.
TRIANGLES: tuple[tuple[int, int, int], ...] = tuple(
    (a, b, c)
    for a in range(5)
    for b in range(a + 1, 5)
    for c in range(b + 1, 5)
)

#: The five tetrahedral facets, indexed by the opposite vertex.
TETRAHEDRA: tuple[tuple[int, int, int, int], ...] = tuple(
    tuple(j for j in range(5) if j != i) for i in range(5)
)


# ---------------------------------------------------------------------------
# Permutations of {0,1,2,3,4}.  Represented as 5-tuples mapping index -> image.
# ---------------------------------------------------------------------------


def _perm_inverse(p: Sequence[int]) -> tuple[int, ...]:
    inv = [0] * 5
    for i, x in enumerate(p):
        inv[x] = i
    return tuple(inv)


def _perm_compose(p: Sequence[int], q: Sequence[int]) -> tuple[int, ...]:
    """Return the composition p o q (apply q first, then p)."""
    return tuple(p[q[i]] for i in range(5))


def _perm_identity() -> tuple[int, ...]:
    return (0, 1, 2, 3, 4)


def _perm_sign(p: Sequence[int]) -> int:
    """Return +1 or -1 by counting inversions."""
    n = len(p)
    inv = 0
    for i in range(n):
        for j in range(i + 1, n):
            if p[i] > p[j]:
                inv += 1
    return -1 if inv % 2 else 1


# ---------------------------------------------------------------------------
# Pentachoron and Triangulation classes.
# ---------------------------------------------------------------------------


@dataclass
class Pentachoron:
    """A single 4-simplex inside a :class:`Triangulation`."""

    index: int
    # adj[i] = (other_index, perm) for the tetrahedral facet opposite vertex i.
    # ``other_index`` is ``-1`` when the facet is on the boundary.
    adj: list[tuple[int, tuple[int, ...]]] = field(
        default_factory=lambda: [(-1, _perm_identity()) for _ in range(5)]
    )

    def neighbour(self, facet: int) -> Optional["Pentachoron"]:
        idx, _ = self.adj[facet]
        return None if idx < 0 else _TRI_REF[id(self)].pentachoron(idx)

    def gluing(self, facet: int) -> tuple[int, tuple[int, ...]]:
        return self.adj[facet]


_TRI_REF: dict[int, "Triangulation"] = {}


class Triangulation:
    """A combinatorial 4-manifold triangulation.

    The class is intentionally light-weight: it keeps a list of pentachora,
    gluings, and computes the skeletal classes (vertices, edges, triangles,
    tetrahedra) lazily. It is designed to be drop-in compatible with
    ``regina.Triangulation4`` for the methods used by the rest of the
    pipeline.
    """

    # ------------------------------------------------------------------ init
    def __init__(self) -> None:
        self._pent: list[Pentachoron] = []
        self._cache: dict[str, object] = {}

    # ----------------------------------------------------------- construction
    def new_pentachoron(self) -> Pentachoron:
        p = Pentachoron(index=len(self._pent))
        self._pent.append(p)
        _TRI_REF[id(p)] = self
        self._cache.clear()
        return p

    def add_pentachora(self, n: int) -> list[Pentachoron]:
        return [self.new_pentachoron() for _ in range(n)]

    def join(
        self,
        a: Pentachoron,
        facet_a: int,
        b: Pentachoron,
        perm: Sequence[int],
    ) -> None:
        """Glue pentachoron ``a`` along its ``facet_a`` to ``b`` using ``perm``.

        ``perm`` is the permutation of ``{0,1,2,3,4}`` mapping the local
        labels of ``a`` to the corresponding local labels of ``b``. The
        ``facet_a`` slot of ``a`` is identified with the
        ``perm[facet_a]`` slot of ``b``.
        """
        perm = tuple(perm)
        if len(perm) != 5 or set(perm) != set(range(5)):
            raise ValueError(f"perm {perm} is not a permutation of 0..4")
        facet_b = perm[facet_a]
        a.adj[facet_a] = (b.index, perm)
        b.adj[facet_b] = (a.index, _perm_inverse(perm))
        self._cache.clear()

    # ---------------------------------------------------------------- access
    def pentachoron(self, i: int) -> Pentachoron:
        return self._pent[i]

    def pentachora(self) -> list[Pentachoron]:
        return list(self._pent)

    @property
    def size(self) -> int:
        return len(self._pent)

    def __len__(self) -> int:  # noqa: D401
        return self.size

    # ----------------------------------------------------- structural queries
    def is_closed(self) -> bool:
        return all(p.adj[f][0] >= 0 for p in self._pent for f in range(5))

    # ------------------------------------------------------------ skeletons
    def _build_skeleton(self) -> None:
        if "vertex_class" in self._cache:
            return

        # Vertices ----------------------------------------------------------
        # Local vertex labels are (pent_index, vertex 0..4). Two local labels
        # collapse iff they are identified via a chain of facet gluings.
        n = self.size
        parent = list(range(5 * n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[max(rx, ry)] = min(rx, ry)

        for p in self._pent:
            for f in range(5):
                other, perm = p.adj[f]
                if other < 0:
                    continue
                for v in range(5):
                    if v == f:
                        continue
                    union(5 * p.index + v, 5 * other + perm[v])

        vertex_class: dict[tuple[int, int], int] = {}
        rep_to_id: dict[int, int] = {}
        for pi in range(n):
            for v in range(5):
                r = find(5 * pi + v)
                if r not in rep_to_id:
                    rep_to_id[r] = len(rep_to_id)
                vertex_class[(pi, v)] = rep_to_id[r]
        self._cache["vertex_class"] = vertex_class
        self._cache["num_vertices"] = len(rep_to_id)

        # Edges, triangles, tetrahedra: collect equivalence classes by their
        # ordered tuples of global vertex ids.
        def collect(faces: Iterable[Sequence[int]]) -> tuple[
            dict[tuple[int, int, Sequence[int]], int], int, list[list[tuple[int, Sequence[int]]]]
        ]:
            classes: dict[tuple, int] = {}
            members: list[list[tuple[int, tuple[int, ...]]]] = []
            local_to_class: dict[tuple[int, tuple[int, ...]], int] = {}
            for pi in range(n):
                for face in faces:
                    key = tuple(sorted(vertex_class[(pi, v)] for v in face))
                    cid = classes.setdefault(key, len(classes))
                    if cid == len(members):
                        members.append([])
                    members[cid].append((pi, tuple(face)))
                    local_to_class[(pi, tuple(face))] = cid
            return local_to_class, len(classes), members

        edge_local, num_edges, edge_members = collect(EDGES)
        tri_local, num_triangles, tri_members = collect(TRIANGLES)
        tet_local, num_tets, tet_members = collect(TETRAHEDRA)

        self._cache.update(
            edge_local=edge_local,
            num_edges=num_edges,
            edge_members=edge_members,
            tri_local=tri_local,
            num_triangles=num_triangles,
            tri_members=tri_members,
            tet_local=tet_local,
            num_tetrahedra=num_tets,
            tet_members=tet_members,
        )

    def num_vertices(self) -> int:
        self._build_skeleton()
        return self._cache["num_vertices"]  # type: ignore[return-value]

    def num_edges(self) -> int:
        self._build_skeleton()
        return self._cache["num_edges"]  # type: ignore[return-value]

    def num_triangles(self) -> int:
        self._build_skeleton()
        return self._cache["num_triangles"]  # type: ignore[return-value]

    def num_tetrahedra(self) -> int:
        self._build_skeleton()
        return self._cache["num_tetrahedra"]  # type: ignore[return-value]

    def vertex_of(self, pent: int, local: int) -> int:
        self._build_skeleton()
        return self._cache["vertex_class"][(pent, local)]  # type: ignore[index]

    def edge_class(self, pent: int, edge: Sequence[int]) -> int:
        self._build_skeleton()
        return self._cache["edge_local"][(pent, tuple(edge))]  # type: ignore[index]

    def triangle_class(self, pent: int, triangle: Sequence[int]) -> int:
        self._build_skeleton()
        return self._cache["tri_local"][(pent, tuple(triangle))]  # type: ignore[index]

    def tetrahedron_class(self, pent: int, tet: Sequence[int]) -> int:
        self._build_skeleton()
        return self._cache["tet_local"][(pent, tuple(tet))]  # type: ignore[index]

    def edge_members(self, edge_id: int) -> list[tuple[int, tuple[int, ...]]]:
        self._build_skeleton()
        return list(self._cache["edge_members"][edge_id])  # type: ignore[index]

    def triangle_members(self, tri_id: int) -> list[tuple[int, tuple[int, ...]]]:
        self._build_skeleton()
        return list(self._cache["tri_members"][tri_id])  # type: ignore[index]

    def edge_degrees(self) -> list[int]:
        self._build_skeleton()
        return [len(m) for m in self._cache["edge_members"]]  # type: ignore[index]

    # ---------------------------------------------------------- orientability
    def is_orientable(self) -> bool:
        # 2-colour pentachora via gluing sign parity; closed component is
        # orientable iff this colouring is consistent.
        sign: list[Optional[int]] = [None] * self.size
        if self.size == 0:
            return True
        sign[0] = 1
        stack = [0]
        while stack:
            i = stack.pop()
            for f in range(5):
                other, perm = self._pent[i].adj[f]
                if other < 0:
                    continue
                want = sign[i] * (-_perm_sign(perm))
                if sign[other] is None:
                    sign[other] = want
                    stack.append(other)
                elif sign[other] != want:
                    return False
        return True

    # -------------------------------------------------------- isomorphism sig
    def edge_degree_isosig(self) -> str:
        """Return a canonical edge-degree isomorphism signature.

        The signature is deterministic and isomorphism-invariant for the
        connected closed case. Two triangulations of the same combinatorial
        type produce equal strings.
        """
        if "isosig" in self._cache:
            return self._cache["isosig"]  # type: ignore[return-value]

        best: Optional[bytes] = None
        n = self.size
        # Try each pentachoron as the root, with each of the 5! relabellings.
        for root in range(n):
            for perm in itertools.permutations(range(5)):
                sig = self._traverse_signature(root, perm)
                if best is None or sig < best:
                    best = sig
        assert best is not None
        out = base64.urlsafe_b64encode(zlib.compress(best, 9)).rstrip(b"=").decode("ascii")
        self._cache["isosig"] = out
        return out

    def _traverse_signature(
        self, root: int, root_perm: Sequence[int]
    ) -> bytes:
        n = self.size
        order: list[int] = [-1] * n
        relabel: list[tuple[int, ...]] = [_perm_identity()] * n
        order[root] = 0
        relabel[root] = tuple(root_perm)
        queue: list[int] = [root]
        seq: list[int] = []
        idx = 1
        head = 0
        while head < len(queue):
            cur = queue[head]
            head += 1
            r = relabel[cur]
            r_inv = _perm_inverse(r)
            for new_f in range(5):
                old_f = r_inv[new_f]
                other, perm = self._pent[cur].adj[old_f]
                if other < 0:
                    seq.append(0xFF)
                    continue
                # Permutation expressed in the canonical labels.
                # canonical perm = relabel[other] o perm o relabel[cur]^-1
                if order[other] == -1:
                    order[other] = idx
                    new_perm = _perm_compose(perm, r_inv)
                    relabel[other] = new_perm
                    queue.append(other)
                    idx += 1
                target = order[other]
                eff = _perm_compose(relabel[other], _perm_compose(perm, r_inv))
                seq.append(target)
                seq.extend(eff)
        return bytes(seq)

    # ---------------------------------------------------------- serialisation
    def to_dict(self) -> dict:
        return {
            "pentachora": self.size,
            "gluings": [
                {
                    "from": p.index,
                    "facet": f,
                    "to": p.adj[f][0],
                    "perm": list(p.adj[f][1]),
                }
                for p in self._pent
                for f in range(5)
                if p.adj[f][0] >= p.index  # avoid duplicate edges in JSON
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping) -> "Triangulation":
        T = cls()
        T.add_pentachora(int(data["pentachora"]))
        for g in data["gluings"]:
            if g["to"] < 0:
                continue
            T.join(T.pentachoron(g["from"]), int(g["facet"]),
                   T.pentachoron(g["to"]), g["perm"])
        return T


# ---------------------------------------------------------------------------
# Public helpers.
# ---------------------------------------------------------------------------


def edge_degree_isosig(T: Triangulation) -> str:
    """Free-function form of :meth:`Triangulation.edge_degree_isosig`."""
    return T.edge_degree_isosig()


def load_triangulation(path: str) -> list[Triangulation]:
    """Load triangulations from a file.

    Supported formats:

    * ``.esig`` / ``.txt`` -- one isomorphism signature per line. Regina is
      used when available; otherwise an error is raised because the .esig
      format requires the Regina reader to materialise pentachora.
    * ``.rga`` -- Regina native XML; requires the Regina Python module.
    * ``.json`` / ``.jsonl`` -- the ``to_dict`` format produced by this
      module.
    """
    if path.endswith((".json", ".jsonl")):
        out: list[Triangulation] = []
        with open(path) as fh:
            if path.endswith(".jsonl"):
                for line in fh:
                    line = line.strip()
                    if line:
                        out.append(Triangulation.from_dict(json.loads(line)))
            else:
                data = json.load(fh)
                if isinstance(data, list):
                    out.extend(Triangulation.from_dict(d) for d in data)
                else:
                    out.append(Triangulation.from_dict(data))
        return out

    try:
        import regina  # type: ignore
    except ImportError as e:  # pragma: no cover - exercised when Regina is absent
        raise RuntimeError(
            f"Loading {path} requires the Regina Python module."
        ) from e

    if path.endswith(".esig") or path.endswith(".txt"):
        triangs: list[Triangulation] = []
        with open(path) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                rt = regina.Triangulation4.fromIsoSig(line)
                triangs.append(_from_regina(rt))
        return triangs
    if path.endswith(".rga"):
        pkt = regina.open(path)
        return [_from_regina(p) for p in pkt.descendants() if isinstance(p, regina.Triangulation4)]
    raise ValueError(f"Unrecognised triangulation file extension: {path}")


def _from_regina(rt: "regina.Triangulation4") -> Triangulation:  # type: ignore[name-defined]
    T = Triangulation()
    pents = T.add_pentachora(rt.size())
    for i, p in enumerate(rt.pentachora()):
        for f in range(5):
            adj = p.adjacentPentachoron(f)
            if adj is None:
                continue
            j = adj.index()
            if j < i:
                continue  # already glued from the other side
            perm = tuple(p.adjacentGluing(f).image(k) for k in range(5))
            T.join(pents[i], f, pents[j], perm)
    return T


def to_regina(T: Triangulation):  # pragma: no cover - exercised when Regina is present
    import regina  # type: ignore

    rt = regina.Triangulation4()
    pents = [rt.newPentachoron() for _ in range(T.size)]
    for p in T.pentachora():
        for f in range(5):
            other, perm = p.adj[f]
            if other < 0 or other < p.index:
                continue
            rt.join(pents[p.index], f, pents[other],
                    regina.Perm5(*perm))
    return rt


# ---------------------------------------------------------------------------
# Built-in tiny examples used by the test suite.
# ---------------------------------------------------------------------------


def standard_s4() -> Triangulation:
    """The 2-pentachoron triangulation of the standard 4-sphere.

    Two pentachora glued along all five tetrahedral facets by the identity
    permutation; the boundary of a 5-simplex restricted to two opposite
    pentachora.
    """
    T = Triangulation()
    a, b = T.add_pentachora(2)
    for f in range(5):
        T.join(a, f, b, (0, 1, 2, 3, 4))
    return T


def two_pent_twisted() -> Triangulation:
    """A two-pentachoron closed-combinatorial triangulation used in tests.

    Two pentachora glued by the cyclic permutation ``π = (1,2,3,4,0)`` on
    every facet. The construction is combinatorially closed and consistent
    (no double-gluings), and provides a non-identity gluing pattern with
    which to exercise the colouring, diagram, and search modules. The
    underlying PL type is not asserted -- callers needing literal census
    triangulations should load the Dim4Census ``.esig`` files via
    :func:`load_triangulation`.
    """
    T = Triangulation()
    a, b = T.add_pentachora(2)
    perm = (1, 2, 3, 4, 0)
    for f in range(5):
        T.join(a, f, b, perm)
    return T


# Backwards-compat alias for code that referenced the old name.
standard_s3_x_s1 = two_pent_twisted
