# bhrt_trisect

End-to-end **C++/CUDA** research pipeline implementing

    triangulation -> ts-tricolouring -> trisection diagram -> invariants

with a CPU beam-search and (optional) CUDA-batched candidate scoring
for trisection-genus reduction. Python is used only where it is
genuinely required (Sage IS Python; Regina ships Python-only .esig
parsers); every algorithm lives in C/C++/CUDA.

## What's in the box

### C/C++/CUDA (primary)

```
cpp/                Real C++ implementation; nothing here calls Python.
  bhrt_trisect.hpp     public C++ API
  bhrt_c_api.{h,cpp}   stable C ABI
  triangulation.cpp    4D simplicial complex, skeleton, isosig
  ts_color.cpp         Spreer-Tillmann ts-tricolouring checker
  diagram.cpp          central surface + cut systems + diagram validator
  invariants.cpp       Smith normal form, H_*, intersection form, signature
                       (exact homology from the three Lagrangians; linking-
                       matrix intersection form in the simply-connected case)
  search_cpu.cpp       beam search; native 1-5/5-1 move executor
  regina_bridge.cpp    Regina-backed executor: 2-4 / 3-3 / 4-4 / 4-2 /
                       edge-collapse (optional, BHRT_HAS_REGINA)
  regina_io.hpp        shared bhrt <-> Regina conversions
  tri_io.cpp           .esig read/write via Regina isosigs (optional)
  bhrt_format.cpp      .bhrt text format reader/writer
  bhrt_cli.c           C-language CLI driver (bhrt-cli)
  bhrt_bench.cpp       C++ benchmark driver (bhrt-bench)
  bhrt_diagram.cpp     invariants from a .tri trisection-diagram file
  bhrt_test.cpp        C++ test suite (bhrt-test); ground-truth manifolds
                       (S^4, CP^2, -CP^2, S^2xS^2, CP^2#CP^2, S^1xS^3) plus
                       Regina- and CUDA-gated backend tests
cuda/               CUDA kernels
  ts_scan.cu           GPU-exhaustive ts-tricolouring search: all 3^(n-1)
                       candidate colourings are filtered on-device ((2,2,1)
                       type, monochromatic triangles, Gamma_k connectivity);
                       survivors are re-verified by the CPU reference, so
                       results are identical to the CPU enumerator
  score_moves.cu       per-candidate scorer (self-verified against the CPU
                       reference on first use; clean fallback without a GPU)
  candidate_pack.cu    on-device record finishing
  smoke.cu             standalone raw-CUDA sanity check
examples/           .tri trisection-diagram inputs (cp2.tri, s2xs2.tri)
docs/               SUMMARY, PRESENTATION, HANDOUT, RUNSHEET, architecture
```

### Move set

| MoveType        | 4D move                | Δ pentachora | executor |
|-----------------|------------------------|--------------|----------|
| `Pachner_1_5`   | 1-5 (pentachoron)      | +4           | native   |
| `Pachner_5_1`   | 5-1 (vertex)           | -4           | native (self-verified) |
| `Pachner_2_3`   | 2-4 (tetrahedron)      | +2           | Regina   |
| `Pachner_3_3`   | 3-3 (triangle)         |  0           | Regina   |
| `Pachner_4_4`   | 4-4 (edge)             |  0           | Regina   |
| `Pachner_4_2`   | 4-2 (edge)             | -2           | Regina   |
| `EdgeCollapse`  | edge collapse          | -deg(e)      | Regina   |

Without Regina the lateral/collapse candidates are skipped safely; the
Regina executor checks legality, runs on a fresh copy, and self-verifies
validity, closedness and Euler characteristic before committing.

### Shell wrappers around the C++ binaries

```
scripts/
  scan_census.sh         Regina .esig -> .bhrt -> bhrt-cli scan-file
  reduce_genus.sh        bhrt-bench --search
  compute_invariants.sh  bhrt-bench --invariants
  reproduce_paper.sh     full census sweep
```

### Python (only what genuinely needs Python)

```
python/bhrt_trisect/
  __init__.py           thin facade
  subprocess_api.py     wrappers around bhrt-cli / bhrt-bench
  regina_bridge.py      Regina .esig / .rga -> .bhrt converter
                        (Regina's Python bindings are the only realistic
                         way to parse Regina isosig strings)
  *.py (other names)    deprecation stubs that raise ImportError and
                        point at the C++ replacement file

sage/                   SageMath exact-algebra helpers (Sage IS Python).
  invariants.py
  twisted_invariants.py
```

## Quick start

```bash
# Build everything (no dependencies).
make

# Run the C++ test suite (all assertions must pass; the count grows with
# the optional backends compiled in).
./build/bhrt-test

# Smoke the CLI on the built-in S^4 example.
./build/bhrt-cli info     # skeleton statistics
./build/bhrt-cli isosig   # canonical edge-degree isomorphism signature
./build/bhrt-cli scan     # ts-tricolouring counters

# Beam-search demo.
./build/bhrt-cli search 64 4 20 0

# Write a .bhrt example file, then process it.
./build/bhrt-cli write-demo /tmp/demo.bhrt
./build/bhrt-cli scan-file  /tmp/demo.bhrt
./build/bhrt-bench /tmp/demo.bhrt --invariants --search

# Invariants from an explicit trisection diagram (.tri file):
./build/bhrt-diagram examples/cp2.tri     # H_*, Q = [1], sigma = 1, odd
./build/bhrt-diagram examples/s2xs2.tri   # H_*, Q = H, sigma = 0, even

# Guided demo (CPU build):
bash demo.sh
```

## Census pipeline (requires Regina)

On a Regina-enabled build the .esig step is pure C++ as well
(`bhrt::loadEsig` / `bhrt::writeEsig` parse and emit Regina isosigs
directly). The Python converter remains as a fallback for builds
without the Regina backend:

```bash
# Fallback only: convert Regina .esig via Regina's Python bindings
python3 -m bhrt_trisect.regina_bridge \
    data/Census/6p-closedOrientable.esig \
    /tmp/6p.bhrt

# Everything after this is pure C++.
./build/bhrt-cli scan-file /tmp/6p.bhrt
./build/bhrt-bench         /tmp/6p.bhrt --invariants --search
```

Or use the shell wrappers, which chain both steps:

```bash
scripts/scan_census.sh data/Census/6p-closedOrientable.esig
scripts/reproduce_paper.sh
```

## Building with optional backends

Two equivalent build systems; keep their build dirs separate. If your
shell exports a stale compiler (conda often sets `CXX=g++-11`), run
`unset CXX CC` first.

```bash
# Makefile (default flags include AddressSanitizer -- the debug build)
make BHRT_HAS_REGINA=1                 # Regina executor + .esig I/O
make BHRT_HAS_REGINA=1 BHRT_HAS_CUDA=1 # + GPU scorer
# NOTE: the Makefile does not track -D flag changes; `make clean` first.

# CMake (Release; tracks option changes)
cmake -S cpp -B build-cmake -DBHRT_HAS_REGINA=ON -DBHRT_HAS_CUDA=ON
cmake --build build-cmake -j
ctest --test-dir build-cmake
```

The GPU paths are opt-in at run time: set `SearchConfig::use_gpu_scorer`
or `BHRT_USE_GPU=1`. Two kernels are used: the move scorer (first batch
verified against the CPU reference; falls back on mismatch) and the
exhaustive ts-tricolouring filter (`enumerateTsTricolouringsGPU`), which
attacks the pipeline's exponential 3^(n-1) bottleneck. The colouring
filter is correctness-neutral by construction: the GPU only discards
candidates that fail the cheap prechecks, and every survivor is run
through the full certified CPU check before being reported. Without a
GPU both paths degrade silently to the CPU implementations.

See `INSTALL.md` for full installation of Regina, SageMath, and CUDA.

## Verifying the build

```
./build/bhrt-test                                  # 36/36 native C++ assertions
PYTHONPATH=python python3 -m pytest tests/ -q      #  5/5 cross-check tests
```

## Scope and known limitations

Stated precisely so results are not over-claimed:

* **Invariants.** Homology (free rank + torsion) is exact for any valid
  Lagrangian diagram. The integral intersection form, signature and parity
  are computed exactly when the diagram is simply connected and
  (g;0,0,0)-standardisable (linking-matrix case: covers S^4, +/-CP^2,
  S^2xS^2, connected sums, ...). The general Feller-Klug-Schirmer-Zemke
  pairing needs geometric curve-crossing data on Sigma that the cellular
  extractor does not record; such inputs return an empty form with an
  audit note rather than a guess.
* **Diagram extraction.** `extractDiagram` cut curves are heuristic and
  not guaranteed to form genuine cut systems; run `validateDiagram` (the
  pipeline does) and trust only validated diagrams. Hand-coded `.tri`
  diagrams bypass this caveat entirely.
* **ts-colour search.** Collapsibility testing is greedy (a `false` is
  inconclusive), and colouring enumeration is exhaustive (3^(n-1)), so
  large censuses need the precheck filters.
* **Known issue under investigation.** On Regina-enabled builds the beam
  search currently aborts with a glibc heap error; the AddressSanitizer
  build (`make BHRT_HAS_REGINA=1`) is the diagnostic. CPU-only builds are
  unaffected.

## Mathematical foundations

* Gay-Kirby, *Trisecting 4-manifolds*.
* Rubinstein-Tillmann, *Multisections of piecewise linear manifolds*.
* Bell-Hass-Rubinstein-Tillmann, *Computing trisections of 4-manifolds*.
* Spreer-Tillmann, *Determining the trisection genus through triangulations*.
* Feller-Klug-Schirmer-Zemke, *Invariants from trisection diagrams*.
* Florens-Moussard, Moussard-Schirmer (twisted invariants).
