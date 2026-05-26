"""bhrt_trisect: triangulation -> ts-tricolouring -> trisection diagram -> invariants,
plus CPU/GPU genus reduction."""

__version__ = "0.1.0"

from .tri_io import (
    Triangulation,
    Pentachoron,
    load_triangulation,
    edge_degree_isosig,
)
from .ts_color import (
    TsColouringResult,
    enumerate_ts_tricolourings,
    is_tricolouring,
    is_c_tricolouring,
    is_ts_tricolouring,
)
from .diagram import TrisectionDiagram, build_diagram
from .invariants import compute_invariants, InvariantBundle
from .bhrt import bhrt_refine

__all__ = [
    "Triangulation",
    "Pentachoron",
    "load_triangulation",
    "edge_degree_isosig",
    "TsColouringResult",
    "enumerate_ts_tricolourings",
    "is_tricolouring",
    "is_c_tricolouring",
    "is_ts_tricolouring",
    "TrisectionDiagram",
    "build_diagram",
    "compute_invariants",
    "InvariantBundle",
    "bhrt_refine",
]
