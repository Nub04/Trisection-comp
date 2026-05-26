"""Thin Python API over the C++ binaries.

The functions here exist for users who want to script the pipeline
from Python (for instance to chain it together with Sage). They are
all pure wrappers around ``./build/bhrt-cli`` and ``./build/bhrt-bench``
-- no algorithm logic lives in Python.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List


def _binary(name: str) -> str:
    here = Path(__file__).resolve().parent.parent.parent
    candidates = [here / "build" / name, Path.cwd() / "build" / name]
    for c in candidates:
        if c.exists() and os.access(c, os.X_OK):
            return str(c)
    on_path = shutil.which(name)
    if on_path:
        return on_path
    raise FileNotFoundError(
        f"could not find {name}; run `make` to build the C++ core first")


def info() -> dict:
    out = subprocess.check_output([_binary("bhrt-cli"), "info"], text=True)
    res = {}
    for line in out.splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            res[k.strip()] = v.strip()
    return res


def isosig() -> str:
    return subprocess.check_output(
        [_binary("bhrt-cli"), "isosig"], text=True).strip()


def scan_file(path: str) -> dict:
    out = subprocess.check_output(
        [_binary("bhrt-cli"), "scan-file", path], text=True)
    res = {}
    for line in out.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            res[k.strip()] = int(v.strip())
    return res


def bench(path: str, *,
          invariants: bool = False, search: bool = False,
          beam: int = 64, height: int = 4, iters: int = 50,
          seed: int = 0) -> List[dict]:
    cmd = [_binary("bhrt-bench"), path]
    if invariants: cmd.append("--invariants")
    if search:     cmd.append("--search")
    cmd += ["--beam", str(beam), "--height", str(height),
            "--iters", str(iters), "--seed", str(seed)]
    out = subprocess.check_output(cmd, text=True)
    return [json.loads(line) for line in out.splitlines() if line.strip()]


def write_demo(path: str) -> None:
    subprocess.check_call([_binary("bhrt-cli"), "write-demo", path])
