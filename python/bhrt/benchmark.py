"""Reproducible benchmark runner.

Emits a JSON report whose keys mirror the metric table in the architecture
note: ts-tricolouring counts, candidate throughput, dedup rate, best
pentachora count, best trisection genus, and timing.

Example::

    from bhrt.benchmark import BenchmarkConfig, run_benchmark
    rep = run_benchmark(BenchmarkConfig(corpus="builtin", seed=42))
"""

from __future__ import annotations

import json
import random
import time
from dataclasses import dataclass, field
from typing import Optional

from .tri_io import Triangulation, standard_s4, standard_s3_x_s1, load_triangulation
from .ts_color import enumerate_ts_tricolourings
from .diagram import build_diagram
from .invariants import compute_invariants
from .search_cpu import (
    BeamSearchConfig,
    PythonExecutor,
    reduce_trisection_genus,
)


@dataclass
class BenchmarkConfig:
    corpus: str = "builtin"            # "builtin" or path to a manifest file
    seed: int = 0
    beam_width: int = 32
    excess_height: int = 4
    max_seconds: float = 5.0
    max_committed: int = 1_000
    include_search: bool = True
    out_path: Optional[str] = None


@dataclass
class BenchmarkEntry:
    name: str
    pentachora: int
    ts_count: int
    diagram_genus: Optional[int]
    invariant: Optional[dict]
    ts_seconds: float
    diagram_seconds: float
    invariant_seconds: float
    search_seconds: Optional[float] = None
    search_best_size: Optional[int] = None
    search_visited: Optional[int] = None

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "pentachora": self.pentachora,
            "ts_count": self.ts_count,
            "diagram_genus": self.diagram_genus,
            "invariant": self.invariant,
            "ts_seconds": self.ts_seconds,
            "diagram_seconds": self.diagram_seconds,
            "invariant_seconds": self.invariant_seconds,
            "search_seconds": self.search_seconds,
            "search_best_size": self.search_best_size,
            "search_visited": self.search_visited,
        }


def _builtin_corpus() -> list[tuple[str, Triangulation]]:
    return [
        ("S^4 (2 pent)", standard_s4()),
        ("S^3 x S^1 (1 pent)", standard_s3_x_s1()),
    ]


def run_benchmark(config: BenchmarkConfig) -> dict:
    rng = random.Random(config.seed)
    if config.corpus == "builtin":
        corpus = _builtin_corpus()
    else:
        triangs = load_triangulation(config.corpus)
        corpus = [(f"{config.corpus}#{i}", T) for i, T in enumerate(triangs)]

    entries: list[BenchmarkEntry] = []
    for name, T in corpus:
        t0 = time.time()
        ts = enumerate_ts_tricolourings(T, max_results=1)
        ts_time = time.time() - t0

        diag_time = 0.0
        inv_time = 0.0
        diagram_genus: Optional[int] = None
        invariant_dict: Optional[dict] = None
        if ts:
            t1 = time.time()
            d = build_diagram(T, ts[0])
            diag_time = time.time() - t1
            diagram_genus = d.genus()
            t2 = time.time()
            inv = compute_invariants(d)
            inv_time = time.time() - t2
            invariant_dict = inv.to_dict()

        entry = BenchmarkEntry(
            name=name,
            pentachora=T.size,
            ts_count=len(ts),
            diagram_genus=diagram_genus,
            invariant=invariant_dict,
            ts_seconds=ts_time,
            diagram_seconds=diag_time,
            invariant_seconds=inv_time,
        )

        if config.include_search:
            t3 = time.time()
            best, report = reduce_trisection_genus(
                T,
                executor=PythonExecutor(),
                config=BeamSearchConfig(
                    beam_width=config.beam_width,
                    excess_height=config.excess_height,
                    max_seconds=config.max_seconds,
                    max_committed=config.max_committed,
                    seed=rng.randint(0, 2**31 - 1),
                ),
            )
            entry.search_seconds = time.time() - t3
            entry.search_best_size = best.triangulation.size
            entry.search_visited = report.visited_isosigs

        entries.append(entry)

    report = {
        "config": {
            "corpus": config.corpus,
            "seed": config.seed,
            "beam_width": config.beam_width,
            "excess_height": config.excess_height,
            "max_seconds": config.max_seconds,
            "max_committed": config.max_committed,
        },
        "entries": [e.to_dict() for e in entries],
    }
    if config.out_path:
        with open(config.out_path, "w") as fh:
            json.dump(report, fh, indent=2)
    return report
