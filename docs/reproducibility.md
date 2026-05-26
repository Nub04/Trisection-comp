# Reproducibility

This document records the exact corpus, software versions, hardware,
and protocol needed to reproduce every table and figure in the paper.

## Software stack

| Component  | Version    | Purpose                                |
|------------|------------|----------------------------------------|
| Regina     | 7.4+       | 4D triangulations, Pachner moves       |
| SageMath   | 10.x       | exact integer linear algebra           |
| Python     | 3.10+      | orchestration, scripts                 |
| NumPy      | 1.26+      | arrays                                 |
| SymPy      | 1.12+      | exact algebra fallback                 |
| CUDA       | 12.0+      | optional GPU scorer                    |
| PyCUDA     | 2024.1+    | optional preferred GPU backend         |
| Numba      | 0.59+      | optional secondary GPU backend         |
| pytest     | 7.x+       | regression suite                       |
| CMake      | 3.20+      | C++ build                              |
| g++/clang  | C++20      | C++ core compiler                      |

Linux is the recommended development host because Regina's
`regina-python` shell command and CUDA toolchain are most fluent
there. macOS and Windows are usable but require manual workarounds
(see Regina's docs).

## Hardware baseline (paper-grade runs)

| Resource | Recommended                                |
|----------|--------------------------------------------|
| CPU      | 16+ cores                                  |
| RAM      | 64-128 GB (Pachner-graph search is memory-heavy) |
| GPU      | NVIDIA, 16-24 GB VRAM                      |
| Disk     | 200 GB free for census + checkpoints       |

## Corpus

All corpora are fetched from the public `Dim4Census` repository at
`https://github.com/regina-normal/Dim4Census`. The exact set used by
the paper is enumerated in `data/manifests/census.json`. Each entry
carries the file path, the documented size, and the expected
ts-tricolouring count that the regression suite validates against.

For full reproducibility, clone `Dim4Census`, place the .esig files
under `bhrt_trisect/data/Census/`, and run:

    python scripts/reproduce_paper.py \
        --out-dir results/paper \
        --manifest data/manifests/census.json

The paper-reproduction script writes JSON tables and the matplotlib
figures referenced from the manuscript.

## Negative-result protocol

Every search run is checkpointed by:

* exact move set used (Pachner 2-3 / 3-3 / 4-4 / 1-5 / 5-1 + edge
  collapse, plus any optional moves enabled by command-line flag);
* excess-height bound (`--height`);
* canonical edge-degree isomorphism signature (Regina canonical form);
* frontier size at each iteration;
* RNG seed list (`--seed`).

This provides the exact bounded-search statement required to publish
negative results in the style of Burke-Burton-Spreer (e.g., their
S^4_C bound of 12 pentachora).

## Determinism

All randomness in the search engine is seeded from `SearchConfig.seed`.
Re-running with the same seed and the same canonical isosig of the
input triangulation produces the same trace. The CUDA kernels are
deterministic given the same launch configuration; the only source of
non-determinism is the `top_k` reduction inside the shortlist kernel,
which is replaced by a deterministic sorted variant when
`BHRT_DETERMINISTIC_TOPK=1` is set in the environment.

## Smoke verification

After installation:

    pytest -q

should run in under a minute and exercise the pure-Python path. The
full census regression requires Regina and the `Dim4Census` data and
takes longer; see `scripts/reproduce_paper.py`.
