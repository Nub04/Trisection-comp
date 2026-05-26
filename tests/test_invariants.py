"""Algorithmic tests moved to native C++ (cpp/bhrt_test.cpp).

Run them with:

    make test          # builds build/bhrt-test
    ./build/bhrt-test  # runs the 36 native assertions

The cross-validation harness that ties the C++ binary back to Python is
in test_cpp_binary.py.
"""
