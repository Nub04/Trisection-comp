#!/usr/bin/env bash
# reduce_genus.sh -- beam-search genus reduction over a corpus.
#
# Usage:
#   scripts/reduce_genus.sh <input.bhrt> [beam] [height] [iters] [seed]
set -euo pipefail
IN=${1:?usage: reduce_genus.sh <input.bhrt> [beam] [height] [iters] [seed]}
BEAM=${2:-128}
HEIGHT=${3:-4}
ITERS=${4:-50}
SEED=${5:-0}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
exec "$ROOT/build/bhrt-bench" "$IN" --search \
    --beam "$BEAM" --height "$HEIGHT" --iters "$ITERS" --seed "$SEED"
