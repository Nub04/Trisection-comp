"""Deprecated: this module's logic now lives in C++.

The implementation moved to cpp/isosig.cpp (or the relevant C++
file). Drive it via the C-language CLI:

    ./build/bhrt-cli       # primary user-facing tool
    ./build/bhrt-bench     # bulk corpus driver
    ./build/bhrt-test      # native C++ test suite

The Python package now only exposes (1) a thin subprocess wrapper around
the C++ binaries (:mod:`bhrt_trisect.subprocess_api`), and (2) a Regina
.esig -> .bhrt converter (:mod:`bhrt_trisect.regina_bridge`) because
Regina ships Python bindings. Sage helpers live in ``sage/``.
"""

raise ImportError(
    "isosig has moved to C++ (cpp/isosig.cpp). "
    "Use ./build/bhrt-cli or bhrt_trisect.subprocess_api."
)
