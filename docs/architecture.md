# Architecture

This document describes the data layout, module boundaries, and
execution model of the BHRT trisection pipeline + GPU genus-reduction
engine.

The high-level mantra: **Regina owns triangulations, Sage owns exact
algebra, CUDA scores move candidates, Python orchestrates.**

## Module map

```
bhrt_trisect/
  cpp/                    C++ core (Regina-backed)
    bhrt_trisect.hpp      public API
    tri_io.cpp            .esig / .rga loading
    ts_color.cpp          Spreer-Tillmann direct checker
    diagram.cpp           central surface + cut-system extractor
    search_cpu.cpp        beam-search baseline (CPU)
    bindings.cpp          pybind11 bindings -> bhrt_trisect_core
    CMakeLists.txt
  cuda/                   CUDA candidate-scoring kernels
    candidate_pack.cu     on-device record finishing
    score_moves.cu        scoring + top-K shortlist
  python/bhrt_trisect/    Python implementation (pure-Python fallback)
    __init__.py
    triangulation.py
    isosig.py
    tri_io.py
    ts_color.py
    bhrt.py
    diagram.py
    invariants.py
    moves.py
    search_cpu.py
    search_gpu.py
    bench.py
    cli.py
  sage/                   Sage helpers (exact algebra)
    invariants.py
    twisted_invariants.py
  data/
    manifests/census.json corpus manifest
  scripts/                user-facing entry points
    scan_census.py
    reduce_genus.py
    compute_invariants.py
    reproduce_paper.py
  tests/                  pytest regression suite
  docs/                   architecture + reproducibility
```

The Python implementation in `python/bhrt_trisect/` is functionally
complete; the C++/CUDA core in `cpp/` and `cuda/` is the performance
target referenced by the report. Both implementations expose the same
public surface, so a Python-only install runs the full mathematical
pipeline (just more slowly) and tests can validate either backend.

## Data flow

```
.esig / .rga  ->  Triangulation (Regina-backed)
                     |
                     v
              ts_color: direct ts-tricolouring
                     |
        success -----+----- failure
            |                  |
            v                  v
      (use directly)   bhrt: BHRT refinement fallback
            |                  |
            +-----+------------+
                  v
             diagram: central surface + cut systems
                  |
                  +--------------------+
                  v                    v
            invariants            search_cpu/gpu
            (Sage exact)          beam-search lowering
```

## Checkpoint record schema

Search-engine checkpoints are line-delimited JSON. Each record has the
following schema (see `python/bhrt_trisect/tri_io.py:CheckpointRecord`):

```json
{
  "id": "uuid4",
  "edge_degree_isosig": "<sorted_degrees>|<isosig>",
  "seed": 0,
  "move_history_hash": "blake2b-128 hex",
  "best_genus": 3,
  "pentachora": 6,
  "status": "frontier",
  "extra": {}
}
```

`status` is one of:

| value             | meaning                                            |
|-------------------|----------------------------------------------------|
| `open`            | newly enqueued, not yet expanded                   |
| `ts_supported`    | direct ts-checker certified a trisection           |
| `rejected`        | excess-height cap or colour-feasibility filter cut |
| `exhausted`       | move enumeration found no legal continuations      |
| `frontier`        | currently in the beam                              |

## GPU pipeline

The GPU layer is a **scorer**, not a topology engine. Every move is
committed on the CPU through Regina's safe API. Concretely:

1. CPU enumerates legal local moves from each frontier state.
2. CPU flattens each candidate into a fixed-width 5-column record:
   `[move_code, locator, delta_pent_estimate, current_pent_count, feasibility]`.
3. GPU kernel `bhrt_score_moves` computes a weighted heuristic score
   per record using `DEFAULT_WEIGHTS` (a 5-vector loaded into constant
   memory in the production build).
4. GPU kernel `bhrt_top_k` shortlists the top-K candidates per parent.
5. CPU applies only the shortlisted candidates through Regina,
   recomputes canonical isosigs, deduplicates, and re-runs the direct
   ts-checker on each child.
6. The Pareto front of `(pentachora, trisection_genus)` is updated.

The CUDA C++ kernel, the Numba CUDA kernel, and the NumPy fallback
produce identical scores for the same input by construction. The GPU
path is selected automatically; environment variable
`BHRT_GPU_BACKEND=cuda|numba|cpu` forces a specific backend.

## Why pure-Python at all

Regina's installation footprint is non-trivial and SageMath is even
heavier. We keep the entire mathematical pipeline runnable in pure
Python (with NumPy + SymPy) so that:

* casual users can experiment with small triangulations on a laptop
  without any C++ toolchain;
* the CI surface is small and reproducible;
* the Python code serves as the executable spec for the C++/CUDA core.

For census-scale runs and any production search the C++/CUDA stack
should be built.
