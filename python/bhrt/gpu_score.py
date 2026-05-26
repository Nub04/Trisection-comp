"""GPU candidate scorer for the beam search.

The CPU search emits *candidate move records* — flattened, fixed-width
numerical descriptors of a proposed local Pachner / collapse edit. The
GPU side ingests a batch of these records, evaluates a scoring function
on each one in parallel, and returns scores back to the CPU which then
commits only the top-K survivors via Regina's safe move API.

Record schema (one row = one candidate, ``dtype=int32``)::

    [ state_id,
      move_kind,
      local_vertex_degree_min,
      local_vertex_degree_max,
      local_pent_count,
      local_facet_external_count,
      delta_pent,            # expected pentachora change
      colour_singletons,     # number of colour-singletons in the local nbhd
      colour_doubletons,     # number of (2,2,1) splits in the local nbhd
      genus_lower_bound,     # cheap proxy for the new diagram genus
    ]

Scoring (weights configurable)::

    score = w_pent     * (-delta_pent)            # prefer pentachora reduction
          + w_genus    * (-genus_lower_bound)     # prefer lower genus
          + w_ts_feas  *  colour_doubletons       # prefer (2,2,1) preservation
          - w_singleton* colour_singletons        # mild penalty
          - w_external * local_facet_external_count  # discourage boundary blow-up

Two backends are provided:

* :class:`NumpyScorer` — always available; vectorised NumPy.
* :class:`NumbaCudaScorer` — Numba-CUDA kernel; activated when
  ``numba.cuda`` is importable and a device is present. Falls back to
  the NumPy scorer on import or device errors.

Both backends present the same interface and produce identical scores.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Record schema.
# ---------------------------------------------------------------------------


REC_FIELDS = (
    "state_id",
    "move_kind",
    "vdeg_min",
    "vdeg_max",
    "local_pent",
    "facet_ext",
    "delta_pent",
    "col_single",
    "col_double",
    "genus_lb",
)
NUM_FIELDS = len(REC_FIELDS)
DTYPE = np.int32


@dataclass
class ScoreWeights:
    w_pent: float = 1.0
    w_genus: float = 1.5
    w_ts_feas: float = 0.5
    w_singleton: float = 0.25
    w_external: float = 0.1


def empty_batch(n: int) -> np.ndarray:
    """Allocate an empty batch of ``n`` candidate records."""
    return np.zeros((n, NUM_FIELDS), dtype=DTYPE)


def pack_record(
    *,
    state_id: int,
    move_kind: int,
    vdeg_min: int,
    vdeg_max: int,
    local_pent: int,
    facet_ext: int,
    delta_pent: int,
    col_single: int,
    col_double: int,
    genus_lb: int,
) -> np.ndarray:
    return np.array(
        [state_id, move_kind, vdeg_min, vdeg_max, local_pent,
         facet_ext, delta_pent, col_single, col_double, genus_lb],
        dtype=DTYPE,
    )


# ---------------------------------------------------------------------------
# NumPy scorer (CPU fallback).
# ---------------------------------------------------------------------------


class NumpyScorer:
    """Vectorised CPU scorer; identical numerics to the GPU backend."""

    def __init__(self, weights: Optional[ScoreWeights] = None) -> None:
        self.weights = weights or ScoreWeights()

    def score(self, batch: np.ndarray) -> np.ndarray:
        w = self.weights
        delta_pent = batch[:, 6].astype(np.float32)
        genus_lb = batch[:, 9].astype(np.float32)
        col_single = batch[:, 7].astype(np.float32)
        col_double = batch[:, 8].astype(np.float32)
        facet_ext = batch[:, 5].astype(np.float32)
        return (
            -w.w_pent * delta_pent
            - w.w_genus * genus_lb
            + w.w_ts_feas * col_double
            - w.w_singleton * col_single
            - w.w_external * facet_ext
        )

    def topk(self, batch: np.ndarray, k: int) -> np.ndarray:
        s = self.score(batch)
        if k >= len(s):
            return np.argsort(-s)
        idx = np.argpartition(-s, k - 1)[:k]
        return idx[np.argsort(-s[idx])]


# ---------------------------------------------------------------------------
# Numba CUDA scorer.
# ---------------------------------------------------------------------------


class NumbaCudaScorer:
    """Numba-CUDA scorer. Falls back to NumPy if CUDA is unavailable."""

    def __init__(self, weights: Optional[ScoreWeights] = None) -> None:
        self.weights = weights or ScoreWeights()
        self._cuda_ok = False
        try:
            from numba import cuda  # type: ignore
            if cuda.is_available():  # pragma: no cover - depends on hardware
                self._cuda_ok = True
                self._cuda = cuda
                self._compile()
        except Exception:
            self._cuda_ok = False
        self._fallback = NumpyScorer(weights)

    def _compile(self) -> None:  # pragma: no cover - GPU-only path
        from numba import cuda

        w = self.weights

        @cuda.jit
        def _score_kernel(batch, scores,
                          w_pent, w_genus, w_ts, w_single, w_ext):
            i = cuda.grid(1)
            if i >= batch.shape[0]:
                return
            delta_pent = float(batch[i, 6])
            genus_lb = float(batch[i, 9])
            col_single = float(batch[i, 7])
            col_double = float(batch[i, 8])
            facet_ext = float(batch[i, 5])
            scores[i] = (
                -w_pent * delta_pent
                - w_genus * genus_lb
                + w_ts * col_double
                - w_single * col_single
                - w_ext * facet_ext
            )

        self._kernel = _score_kernel

    def score(self, batch: np.ndarray) -> np.ndarray:
        if not self._cuda_ok:
            return self._fallback.score(batch)
        # pragma: no cover - executed only on a GPU host
        cuda = self._cuda
        n = batch.shape[0]
        d_batch = cuda.to_device(batch.astype(np.int32))
        d_scores = cuda.device_array(n, dtype=np.float32)
        tpb = 128
        bpg = (n + tpb - 1) // tpb
        self._kernel[bpg, tpb](
            d_batch, d_scores,
            np.float32(self.weights.w_pent),
            np.float32(self.weights.w_genus),
            np.float32(self.weights.w_ts_feas),
            np.float32(self.weights.w_singleton),
            np.float32(self.weights.w_external),
        )
        return d_scores.copy_to_host()

    def topk(self, batch: np.ndarray, k: int) -> np.ndarray:
        s = self.score(batch)
        if k >= len(s):
            return np.argsort(-s)
        idx = np.argpartition(-s, k - 1)[:k]
        return idx[np.argsort(-s[idx])]


# ---------------------------------------------------------------------------
# Helpers for candidate-record construction from a Triangulation.
# ---------------------------------------------------------------------------


MOVE_KIND_ID = {
    "1-5": 1,
    "2-4": 2,
    "3-3": 3,
    "4-4": 4,
    "4-2": 5,
    "5-1": 6,
    "collapse": 7,
}


def build_candidate_batch(
    candidates: Iterable[tuple[int, str, dict]],
) -> np.ndarray:
    """Pack an iterable of ``(state_id, move_kind, fields)`` into a batch."""
    rows: list[np.ndarray] = []
    for state_id, kind, fields in candidates:
        rows.append(pack_record(
            state_id=state_id,
            move_kind=MOVE_KIND_ID.get(kind, 0),
            vdeg_min=fields.get("vdeg_min", 0),
            vdeg_max=fields.get("vdeg_max", 0),
            local_pent=fields.get("local_pent", 0),
            facet_ext=fields.get("facet_ext", 0),
            delta_pent=fields.get("delta_pent", 0),
            col_single=fields.get("col_single", 0),
            col_double=fields.get("col_double", 0),
            genus_lb=fields.get("genus_lb", 0),
        ))
    if not rows:
        return empty_batch(0)
    return np.stack(rows, axis=0)


def default_scorer(weights: Optional[ScoreWeights] = None):
    """Return the best available scorer (Numba CUDA if present, else NumPy)."""
    candidate = NumbaCudaScorer(weights)
    if candidate._cuda_ok:
        return candidate
    return NumpyScorer(weights)
