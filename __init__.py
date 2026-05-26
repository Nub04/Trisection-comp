"""bhrt_trisect: thin Python surface over the C++ core.

All algorithmic logic now lives in C++ (cpp/) and CUDA (cuda/). This
Python package exposes only:

* :mod:`bhrt_trisect.subprocess_api` -- wrapper functions that call
  the compiled ``bhrt-cli`` / ``bhrt-bench`` binaries.
* :mod:`bhrt_trisect.regina_bridge`  -- .esig / .rga converter that
  uses Regina's Python bindings (the only realistic .esig parser).

For the exact-algebra helpers see ``sage/`` (Sage IS Python; no choice).

Old algorithmic modules (triangulation, ts_color, diagram, invariants,
search_cpu, ...) now raise ImportError pointing at the C++ replacement.
"""

from .subprocess_api import bench, info, isosig, scan_file, write_demo
from . import regina_bridge

__version__ = "0.2.0"
__all__ = ["bench", "info", "isosig", "scan_file", "write_demo", "regina_bridge"]
