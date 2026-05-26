# Python in `bhrt_trisect` — what is live vs. deprecated

The mathematical pipeline runs in the C++ core (`cpp/`). Python is kept for
three jobs only. **Regina and Sage are required for the jobs marked LIVE.**

## LIVE — keep

| File / dir | Role | Needs |
|---|---|---|
| `bhrt_trisect/regina_bridge.py` | Convert Regina `.esig` / `.rga` census files to `.bhrt` (Regina's isosig parser lives in its Python bindings) | Regina 7.4+ Python bindings (`import regina`) |
| `bhrt_trisect/subprocess_api.py` | Thin Python wrappers that drive the compiled `bhrt-cli` / `bhrt-bench` binaries | the built C++ binaries |
| `../sage/invariants.py`, `../sage/twisted_invariants.py` | Exact-algebra helpers (twisted invariants, exact SNF) | SageMath |

Regina is **also** used at the C++ level by `cpp/tri_io.cpp` and the lateral
Pachner moves when you build with `-DBHRT_HAS_REGINA=ON` (CMake).

## DEPRECATED — safe to ignore or delete

These are `ImportError` stubs left from the C++ migration, plus an older
parallel pure-Python implementation that the C++ core supersedes. They are not
imported by anything LIVE above.

- `bhrt_trisect/ts_color.py`, `diagram.py`, `bhrt.py`, `invariants.py`,
  `moves.py`, `search_cpu.py`, `search_gpu.py`, `triangulation.py`,
  `tri_io.py`, `isosig.py`, `ts_color.py`, `bench.py`, `cli.py` — stubs that
  raise `ImportError` pointing at the C++ replacement.
- `python/bhrt/` — the earlier full pure-Python reference implementation
  (kept only as an executable spec; not needed to run or test the tool).
- root-level `*.py` in the project (e.g. `search_cpu.py`, `triangulation.py`)
  — older stubs.

To remove the dead weight without touching anything LIVE:

```bash
# from bhrt_trisect/
rm -rf python/bhrt
# delete the stub modules but KEEP regina_bridge.py and subprocess_api.py
```

(Leave `sage/` and `regina_bridge.py` in place.)
