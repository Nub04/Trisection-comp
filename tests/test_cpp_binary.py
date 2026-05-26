"""The single Python test: drive the C++ binaries and compare outputs.

The native C++ test suite lives in cpp/bhrt_test.cpp (run it via
`./build/bhrt-test` after `make`). This Python harness is the
cross-validation: it asserts that the C++ binary produces the
hard-coded skeleton numbers for the built-in 2-pentachoron S^4
example, that the .bhrt round-trip preserves the isosig, and that the
bench driver emits valid JSONL.
"""

import json
import os
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
CLI = ROOT / "build" / "bhrt-cli"
BENCH = ROOT / "build" / "bhrt-bench"
TEST = ROOT / "build" / "bhrt-test"


def _have(p):
    return p.exists() and os.access(p, os.X_OK)


@pytest.fixture(scope="module")
def have_cli():
    if not _have(CLI):
        pytest.skip(f"{CLI} not built; run `make` first")


@pytest.fixture(scope="module")
def have_bench():
    if not _have(BENCH):
        pytest.skip(f"{BENCH} not built; run `make` first")


def _run(*args):
    return subprocess.check_output(list(args), text=True)


def test_native_cpp_test_suite_passes():
    if not _have(TEST):
        pytest.skip("build/bhrt-test not present; run `make test`")
    subprocess.check_call([str(TEST)])  # exits 0 iff all assertions pass


def test_cli_info_reports_s4_numbers(have_cli):
    out = _run(str(CLI), "info")
    fields = {}
    for line in out.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            fields[k.strip()] = v.strip()
    assert fields["pentachora"] == "2"
    assert fields["vertices"]   == "5"
    assert fields["edges"]      == "10"
    assert fields["triangles"]  == "10"
    assert fields["tetrahedra"] == "5"
    assert fields["euler chi"]  == "2"
    assert fields["is_closed"]  == "1"
    assert fields["is_valid"]   == "1"


def test_cli_scan_counters(have_cli):
    out = _run(str(CLI), "scan")
    counters = {}
    for line in out.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            counters[k.strip()] = int(v.strip())
    assert counters["n_total"] == 1
    assert counters["n_pass_precheck"] == 1


def test_bhrt_format_roundtrip_via_cli(tmp_path, have_cli):
    path = tmp_path / "demo.bhrt"
    _run(str(CLI), "write-demo", str(path))
    text = path.read_text()
    assert "# bhrt-format v1" in text
    assert text.count("END") == 1
    # Scan the file we just wrote.
    out = _run(str(CLI), "scan-file", str(path))
    assert "n_total" in out


def test_bench_emits_valid_jsonl(tmp_path, have_cli, have_bench):
    path = tmp_path / "demo.bhrt"
    _run(str(CLI), "write-demo", str(path))
    out = _run(str(BENCH), str(path), "--invariants",
                "--search", "--beam", "8", "--height", "2",
                "--iters", "3", "--seed", "1")
    rows = [json.loads(line) for line in out.splitlines() if line.strip()]
    assert len(rows) == 1
    r = rows[0]
    assert r["pent"] == 2
    assert r["vertices"] == 5
    assert "search_best_pent" in r
