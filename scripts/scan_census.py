#!/usr/bin/env python3
"""Deprecated. This script is now a shell wrapper around the C++ binary.

Use the corresponding .sh file instead:

    scan_census.py        -> scripts/scan_census.sh
    reduce_genus.py       -> scripts/reduce_genus.sh
    compute_invariants.py -> scripts/compute_invariants.sh
    reproduce_paper.py    -> scripts/reproduce_paper.sh
"""

import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
script = os.path.join(here, os.path.basename(__file__).replace(".py", ".sh"))
sys.stderr.write(
    f"deprecated: use {script} instead "
    "(the algorithmic logic is in the C++ binary build/bhrt-cli)\n")
sys.exit(2)
