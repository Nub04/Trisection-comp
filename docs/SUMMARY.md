# What this program does

`bhrt_trisect` is a research pipeline that takes a **triangulated closed
4-manifold** and computes a **trisection** of it — a decomposition into
three simple pieces introduced by Gay–Kirby — together with the algebraic
invariants of the manifold. It implements the combinatorial bridge between
triangulations and trisections from Spreer–Tillmann and
Bell–Hass–Rubinstein–Tillmann (BHRT), plus a search engine that tries to
lower the genus of the trisection.

The whole computation is one arrow:

```
triangulation  →  ts-tricolouring  →  trisection diagram  →  invariants
                                                         ↘  genus-reduction search
```

## The pipeline, stage by stage

**1. Triangulation.** A closed 4-manifold is stored as a set of
pentachora (4-simplices) with face gluings (`cpp/triangulation.cpp`).
Inputs come from Regina census files (`.esig`) via a small Python bridge,
or from the plain-text `.bhrt` format. The code computes the full skeleton
(vertices, edges, triangles, tetrahedra), the Euler characteristic, and a
canonical isomorphism signature for deduplication.

**2. ts-tricolouring** (`cpp/ts_color.cpp`). Following Spreer–Tillmann, it
3-colours the vertices and checks the conditions that make the colouring
encode a trisection: every pentachoron is type (2,2,1); no triangle is
monochromatic; each colour class is connected; and each bi-coloured
2-complex collapses to a spine. A colouring passing all checks certifies a
trisection. The code can also enumerate every such colouring or scan a
whole census counting how many triangulations admit one.

**3. Trisection diagram** (`cpp/diagram.cpp`). From a certified colouring
it builds the central surface Σ (one quadrilateral per pentachoron) and
reads off the three cut systems (α, β, γ). A validator
(`validateDiagram`) checks that Σ is a closed, orientable, connected
surface and that each of α, β, γ is a genuine cut system of *g* disjoint
closed curves.

**4. Invariants** (`cpp/invariants.cpp`). From the cut curves’ homology
classes in H₁(Σ;ℤ) = ℤ²ᵍ it computes the homology and intersection form of
the 4-manifold:

- H₁(X) as the cokernel of the three Lagrangians (exact, over ℤ, via Smith
  normal form);
- the Betti numbers b₂, b₃ and the torsion of H₂/H₃ (exact);
- the integral **intersection form**, its **signature** and **parity**,
  computed by the classical linking-matrix method for simply-connected
  diagrams (the case covering S⁴, ℂP², S²×S², connected sums).

**5. Genus-reduction search** (`cpp/search_cpu.cpp`). A beam search walks
the Pachner-move graph, evaluating the trisection genus of each
triangulation it reaches and keeping a Pareto front of
(pentachora, genus). An optional CUDA layer (`cuda/`) scores move
candidates in batch: it is a scorer only — every move is committed on the
CPU. As of this version the GPU scorer is wired into `beamSearch` behind
`SearchConfig::use_gpu_scorer` (build with `make BHRT_HAS_CUDA=1`); it
produces bit-identical scores to the CPU path and falls back to the CPU
automatically when no GPU is present, so results never depend on whether
CUDA is used.

## What is rigorous today

- **ts-tricolouring detection** — a faithful implementation of the
  Spreer–Tillmann test.
- **Homology of X** — H₁, b₂, b₃ and torsion are computed exactly from the
  Lagrangian data and are validated against S⁴, ℂP², ℂP²#ℂP², S²×S² and
  S¹×S³ in the test suite (`cpp/bhrt_test.cpp`).
- **Intersection form / signature / parity** — exact for simply-connected,
  (g;0,0,0)-standardisable diagrams; verified to give σ(ℂP²)=1 (odd),
  σ(S²×S²)=0 (even), σ(ℂP²#ℂP²)=2.
- **Move executor** — the 1→5 Pachner move is implemented and is checked to
  preserve the manifold (χ stays 2 on S⁴); its inverse 5→1 is detected and
  **self-verified** by re-expansion, so the search now genuinely explores
  the move graph instead of standing still.
- **Diagram validators** — surface closedness/orientability/connectivity
  and cut-system checks.

## What still needs work / external backends

- The **lateral Pachner moves** (2-4, 3-3, 4-2) are only available through
  the Regina-backed executor (build with `-DBHRT_HAS_REGINA=ON`). Without
  them the built-in search can grow and shrink (1↔5) but has limited
  lateral mobility.
- The **general** intersection form (non-simply-connected, or
  non-standardisable diagrams) needs the full FKSZ algorithm, which uses
  the geometric crossing data of the curves on Σ — the cellular extractor
  does not yet record a rotation system, so `computeInvariants(diagram)`
  currently reports homology only and defers the form to
  `computeInvariantsFromClasses` (which takes the curve classes directly).
- Census-scale input parsing relies on Regina’s Python bindings.

## How to drive it

```bash
make && ./build/bhrt-test        # native C++ tests (incl. ground truth)
./build/bhrt-cli info            # skeleton stats on the built-in S^4
./build/bhrt-cli search 64 4 20 0
python3 -m bhrt_trisect.regina_bridge data/Census/6p.esig /tmp/6p.bhrt
./build/bhrt-bench /tmp/6p.bhrt --invariants --search
```

## References

Gay–Kirby, *Trisecting 4-manifolds*; Rubinstein–Tillmann,
*Multisections of PL manifolds*; Bell–Hass–Rubinstein–Tillmann,
*Computing trisections of 4-manifolds*; Spreer–Tillmann,
*Determining the trisection genus through triangulations*;
Feller–Klug–Schirmer–Zemke, *Calculating the homology and intersection
form of a 4-manifold from a trisection diagram* (arXiv:1711.04762).
