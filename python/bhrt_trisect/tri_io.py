"""Deprecated: triangulation I/O moved to C++.

The .bhrt text format reader/writer is in ``cpp/bhrt_format.cpp``.
Regina .esig conversion is in :mod:`bhrt_trisect.regina_bridge`.
"""

from .regina_bridge import esig_to_bhrt, rga_to_bhrt  # re-export

__all__ = ["esig_to_bhrt", "rga_to_bhrt"]
