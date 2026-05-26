#!/usr/bin/env bash
# scan_census.sh -- ts-tricolouring scan of a Regina .esig census.
#
# Pipeline:
#   1. Convert the .esig to our portable .bhrt format using Regina's
#      Python bindings (this is the only step that needs Python).
#   2. Run the C++ bhrt-cli scan-file subcommand on the .bhrt result.
#
# Usage:
#   scripts/scan_census.sh <input.esig> [<intermediate.bhrt>]
set -euo pipefail
IN=${1:?usage: scan_census.sh <input.esig> [<intermediate.bhrt>]}
BHRT=${2:-/tmp/scan_census.bhrt}
ROOT=$(cd "$(dirname "$0")/.." && pwd)

echo "# stage 1: convert $IN -> $BHRT (via Regina)" >&2
python3 -m bhrt_trisect.regina_bridge "$IN" "$BHRT"

echo "# stage 2: scan $BHRT (C++)" >&2
exec "$ROOT/build/bhrt-cli" scan-file "$BHRT"
