#!/usr/bin/env bash
# compute_invariants.sh -- run the FKSZ invariants engine over a corpus.
#
# Output: one JSONL row per triangulation. When invoked under
#   sage -bash scripts/compute_invariants.sh ...
# the C++ engine still does the integer linear algebra; Sage is wired
# in via the optional sage/ helpers for twisted invariants.
#
# Usage:
#   scripts/compute_invariants.sh <input.bhrt>
set -euo pipefail
IN=${1:?usage: compute_invariants.sh <input.bhrt>}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
exec "$ROOT/build/bhrt-bench" "$IN" --invariants
