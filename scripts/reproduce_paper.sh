#!/usr/bin/env bash
# reproduce_paper.sh -- rebuild every paper table/figure.
#
# Walks every corpus in data/manifests/census.json, converts .esig to
# .bhrt via Regina, then runs scan + invariants + bounded search on
# each. Output JSON goes to $OUT_DIR.
set -euo pipefail
OUT_DIR=${OUT_DIR:-results/paper}
ROOT=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$OUT_DIR"

for path in "$ROOT"/data/Census/2p-closedOrientable.esig \
            "$ROOT"/data/Census/4p-closedOrientable.esig \
            "$ROOT"/data/Census/6p-closedOrientable.esig; do
    [ -f "$path" ] || { echo "skip (missing): $path" >&2; continue; }
    name=$(basename "$path" .esig)
    bhrt="$OUT_DIR/$name.bhrt"
    python3 -m bhrt_trisect.regina_bridge "$path" "$bhrt"
    "$ROOT/build/bhrt-cli"   scan-file "$bhrt"   > "$OUT_DIR/${name}_scan.txt"
    "$ROOT/build/bhrt-bench" "$bhrt" --invariants > "$OUT_DIR/${name}_invariants.jsonl"
done
echo "# done. results in $OUT_DIR"
