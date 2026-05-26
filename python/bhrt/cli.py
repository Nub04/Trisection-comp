"""Command-line entry points for the bhrt_trisect toolchain.

Subcommands::

    bhrt scan-census    INPUT  --out FILE
    bhrt reduce-genus   --input FILE --beam 64 --height 4 [--device 0]
    bhrt invariants     --input FILE --out FILE
    bhrt benchmark      [--corpus FILE] [--seed N]
    bhrt smoke                       # built-in S^4 / S^3xS^1 end-to-end demo
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Sequence

from .tri_io import (
    Triangulation,
    load_triangulation,
    standard_s4,
    standard_s3_x_s1,
    edge_degree_isosig,
)
from .ts_color import enumerate_ts_tricolourings
from .diagram import build_diagram
from .invariants import compute_invariants
from .benchmark import BenchmarkConfig, run_benchmark
from .search_cpu import (
    BeamSearchConfig,
    PythonExecutor,
    reduce_trisection_genus,
)


# ---------------------------------------------------------------------------
# Helpers.
# ---------------------------------------------------------------------------


def _emit_scan_record(name: str, T: Triangulation) -> dict:
    ts = enumerate_ts_tricolourings(T)
    record = {
        "name": name,
        "isosig": edge_degree_isosig(T),
        "pentachora": T.size,
        "vertices": T.num_vertices(),
        "edges": T.num_edges(),
        "ts_count": len(ts),
        "diagrams": [],
    }
    for r in ts:
        diag = build_diagram(T, r)
        record["diagrams"].append({
            "partition": [list(p) for p in r.partition],
            "genus": diag.genus(),
            "vertices": diag.num_vertices(),
            "edges": diag.num_edges(),
            "faces": diag.num_faces(),
        })
    return record


# ---------------------------------------------------------------------------
# Subcommand implementations.
# ---------------------------------------------------------------------------


def cmd_scan_census(args: argparse.Namespace) -> int:
    if args.input == "builtin":
        triangs: list[tuple[str, Triangulation]] = [
            ("S^4", standard_s4()),
            ("S^3xS^1", standard_s3_x_s1()),
        ]
    else:
        loaded = load_triangulation(args.input)
        triangs = [(f"{args.input}#{i}", T) for i, T in enumerate(loaded)]
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    n_ts = 0
    with open(args.out, "w") as fh:
        for name, T in triangs:
            rec = _emit_scan_record(name, T)
            n_ts += rec["ts_count"]
            fh.write(json.dumps(rec) + "\n")
    print(f"scanned {len(triangs)} triangulations, wrote ts-tricolourings={n_ts} to {args.out}")
    return 0


def cmd_reduce_genus(args: argparse.Namespace) -> int:
    if args.input == "builtin":
        triangs: list[tuple[str, Triangulation]] = [
            ("S^4", standard_s4()),
            ("S^3xS^1", standard_s3_x_s1()),
        ]
    else:
        triangs = []
        with open(args.input) as fh:
            for line in fh:
                rec = json.loads(line)
                name = rec.get("name", "<anon>")
                if "triangulation" in rec:
                    triangs.append((name, Triangulation.from_dict(rec["triangulation"])))
                else:
                    print(f"skipping {name}: no embedded triangulation", file=sys.stderr)
    cfg = BeamSearchConfig(
        beam_width=args.beam,
        excess_height=args.height,
        max_seconds=args.max_seconds,
        max_committed=args.max_committed,
        seed=args.seed,
    )
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True) if args.out else None
    out_fh = open(args.out, "w") if args.out else sys.stdout
    try:
        for name, T in triangs:
            best, report = reduce_trisection_genus(T, executor=PythonExecutor(), config=cfg)
            row = {
                "name": name,
                "start_size": T.size,
                "best_size": best.triangulation.size,
                "best_isosig": best.isosig,
                "report": report.to_dict(),
            }
            out_fh.write(json.dumps(row) + "\n")
    finally:
        if args.out:
            out_fh.close()
    return 0


def cmd_invariants(args: argparse.Namespace) -> int:
    out_records: list[dict] = []
    with open(args.input) as fh:
        for line in fh:
            rec = json.loads(line)
            name = rec.get("name", "<anon>")
            if "triangulation" not in rec:
                continue
            T = Triangulation.from_dict(rec["triangulation"])
            ts = enumerate_ts_tricolourings(T, max_results=1)
            if not ts:
                out_records.append({"name": name, "status": "no ts-tricolouring"})
                continue
            d = build_diagram(T, ts[0])
            inv = compute_invariants(d)
            out_records.append({"name": name, "invariant": inv.to_dict()})
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True) if args.out else None
    if args.out:
        with open(args.out, "w") as fh:
            json.dump(out_records, fh, indent=2)
    else:
        json.dump(out_records, sys.stdout, indent=2)
    return 0


def cmd_benchmark(args: argparse.Namespace) -> int:
    rep = run_benchmark(BenchmarkConfig(
        corpus=args.corpus or "builtin",
        seed=args.seed,
        beam_width=args.beam,
        excess_height=args.height,
        max_seconds=args.max_seconds,
        max_committed=args.max_committed,
        include_search=not args.no_search,
        out_path=args.out,
    ))
    if not args.out:
        json.dump(rep, sys.stdout, indent=2)
    return 0


def cmd_smoke(args: argparse.Namespace) -> int:
    for name, T in [("S^4", standard_s4()), ("S^3xS^1", standard_s3_x_s1())]:
        print(f"=== {name} ===")
        print(f"  pentachora: {T.size}")
        print(f"  vertices:   {T.num_vertices()}")
        print(f"  closed:     {T.is_closed()}")
        print(f"  orientable: {T.is_orientable()}")
        print(f"  isosig:     {edge_degree_isosig(T)}")
        ts = enumerate_ts_tricolourings(T, max_results=1)
        print(f"  ts-tricolourings (>=1): {len(ts)}")
        if ts:
            d = build_diagram(T, ts[0])
            inv = compute_invariants(d)
            print(f"  diagram: V={d.num_vertices()} E={d.num_edges()} F={d.num_faces()} g={d.genus()}")
            print(f"  invariants: {json.dumps(inv.to_dict(), indent=2)}")
    return 0


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="bhrt", description="BHRT 4-manifold trisection toolchain")
    sub = parser.add_subparsers(dest="command", required=True)

    p_scan = sub.add_parser("scan-census", help="Detect ts-tricolourings and extract diagrams.")
    p_scan.add_argument("input", help="Path to census file or 'builtin'.")
    p_scan.add_argument("--out", required=True, help="Output JSONL path.")
    p_scan.set_defaults(func=cmd_scan_census)

    p_red = sub.add_parser("reduce-genus", help="Beam-search to lower triangulation size / trisection genus.")
    p_red.add_argument("--input", required=True)
    p_red.add_argument("--out", default=None)
    p_red.add_argument("--beam", type=int, default=64)
    p_red.add_argument("--height", type=int, default=4)
    p_red.add_argument("--max-seconds", type=float, default=30.0)
    p_red.add_argument("--max-committed", type=int, default=10_000)
    p_red.add_argument("--seed", type=int, default=0)
    p_red.add_argument("--device", type=int, default=0, help="GPU device id (ignored by CPU executor).")
    p_red.set_defaults(func=cmd_reduce_genus)

    p_inv = sub.add_parser("invariants", help="Compute invariant bundle for each triangulation.")
    p_inv.add_argument("--input", required=True)
    p_inv.add_argument("--out", default=None)
    p_inv.set_defaults(func=cmd_invariants)

    p_bench = sub.add_parser("benchmark", help="Run the reproducible benchmark suite.")
    p_bench.add_argument("--corpus", default=None)
    p_bench.add_argument("--seed", type=int, default=0)
    p_bench.add_argument("--beam", type=int, default=32)
    p_bench.add_argument("--height", type=int, default=4)
    p_bench.add_argument("--max-seconds", type=float, default=5.0)
    p_bench.add_argument("--max-committed", type=int, default=1000)
    p_bench.add_argument("--no-search", action="store_true")
    p_bench.add_argument("--out", default=None)
    p_bench.set_defaults(func=cmd_benchmark)

    p_smoke = sub.add_parser("smoke", help="Built-in end-to-end demo.")
    p_smoke.set_defaults(func=cmd_smoke)

    args = parser.parse_args(argv)
    return int(args.func(args) or 0)


if __name__ == "__main__":
    sys.exit(main())
