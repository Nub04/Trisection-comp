# Trisection-genus census

`bhrt-bench --census` runs a GPU-parallel census of the **trisection genus** of
triangulated 4-manifolds.

## What it computes

A *trisection* (Gay-Kirby) splits a closed 4-manifold into three 4-dimensional
handlebodies meeting in a central surface `Σ`; the **trisection genus** `g(M)`
is the minimal genus of `Σ` over all trisections. Bell-Hass-Rubinstein-Tillmann
(arXiv:1711.02763) and Spreer-Tillmann (arXiv:1712.04915) compute trisections
combinatorially from *ts-tricolourings* of a triangulation: a 3-colouring of the
vertices in which every pentachoron is bicoloured `(2,2,1)` and the induced
bi-coloured 2-complexes collapse to graphs. Each ts-tricolouring gives a
trisection whose central surface has genus `g(Σ) = (2 − χ(Σ))/2`; the trisection
genus of the triangulation is the minimum of `g(Σ)` over its ts-tricolourings.

The census enumerates ts-tricolourings on the GPU (the `3^(n−1)` colouring
space, with the cheap `(2,2,1)` / monochromatic / connectivity prechecks done on
device) and reports, per node, the minimum genus and the full genus multiset.

> Only the central surface and its genus are used. `g(Σ)` needs the colouring and
> `Σ` — both built correctly by `extractDiagram` — and does **not** depend on the
> cut-curve extraction, which is still degenerate (see
> `docs/SCOPE_AND_LIMITATIONS.md`). The census therefore reports genus, not the
> full FKSZ intersection form.

## Methodology

The design follows Burton's experimental-mathematics template for 3-manifolds
(arXiv:1110.6080, *Simplification paths in the Pachner graphs of closed
orientable 3-manifold triangulations*), lifted to 4-manifold trisections:

- **Nodes** are 4-manifold triangulations, deduplicated by the *complete*
  isomorphism signature `isoSig()` (Burton's `σ`). `edgeDegreeIsoSig` is only an
  edge-degree invariant and is deliberately not used for dedup.
- **`--census-generate N`** walks the Pachner graph outward from `S⁴` via the
  native 1-5/5-1 moves up to `N` pentachora (the full 2-4/3-3/4-4/4-2 catalogue
  is available in Regina builds), forming a finite level set.
- **`--census <file>`** runs over an explicit node list — `.bhrt` (native) or
  `.esig` (Regina isomorphism signatures).
- The **GPU** performs the per-node ts-tricolouring scan, the step that
  dominates cost and grows as `3^(n−1)`.

The new contribution is the GPU-parallel ts-tricolouring scan; Burton and
Spreer-Tillmann supply the method and the reference data this validates against.

## Validation data (paper-sourced)

`data/spreer_tillmann.esig` holds the verbatim Regina isomorphism signatures
from Spreer-Tillmann, SoCG 2018, Appendix A, for `S⁴`, `CP²` and `S²×S²`.
`data/known_trisection_genus.txt` holds the published genera.

| Manifold | Pentachora | Expected min genus `g(M)` | Source |
|----------|-----------:|--------------------------:|--------|
| S⁴       | 2 / —      | 0  | Gay-Kirby 2016 |
| CP²      | 8          | 1  | Spreer-Tillmann 2018 |
| S²×S²    | 14         | 2  | Spreer-Tillmann 2018 |
| K3       | 134        | 22 | Spreer-Tillmann 2018, Thm 1 |

Each simple-crystallisation triangulation admits 15 ts-tricolourings. The K3
signature (~700 chars) is not inlined; paste it from the paper's Appendix A as a
`K3|<sig>` line in the `.esig` to include it.

## Build and run

CPU (CMake), Regina enabled so `.esig` loads:

```bash
unset CXX CC
rm -rf build-cmake
cmake -S cpp -B build-cmake -DCMAKE_BUILD_TYPE=Release -DBHRT_HAS_REGINA=ON -DBHRT_HAS_CUDA=OFF
cmake --build build-cmake -j
LD_LIBRARY_PATH=/usr/local/lib ./build-cmake/bhrt-bench \
    --census data/spreer_tillmann.esig --known data/known_trisection_genus.txt
```

GPU (nvcc-direct) for the parallel scan, native `.bhrt` / generate input:

```bash
bash build_gpu.sh
./build-gpu/bhrt-bench --census-generate 6 --gpu     # machinery check (all S^4)
```

Expected: `CP2 → min_genus 1`, `S2xS2 → 2`, `S4 → 0`; the final
`{"summary":true,...}` line carries the genus histogram and, with `--known`,
`validation_pass` / `validation_fail`.

## Scope and honesty

- Genus only; the full FKSZ intersection form needs genuine cut curves (open
  problem — see `docs/SCOPE_AND_LIMITATIONS.md`).
- `--census-generate` reaches only `S⁴` nodes via the native moves, so it
  validates the machinery (dedup, GPU≡CPU, genus tabulation), not new genera;
  other manifolds' Pachner graphs need a Regina build.
- `.esig` input requires a Regina build; a GPU scan on `.esig` data means
  converting to `.bhrt` first, or scanning on the CPU Regina build.
- Census patterns are experimental data, not theorems — cf. Burton's own caveat
  that his "≤ 2 extra tetrahedra" regularity did not fully generalise.
