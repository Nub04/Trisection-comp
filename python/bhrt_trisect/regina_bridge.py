"""Regina .esig / .rga -> .bhrt converter.

Regina ships first-party Python bindings (``import regina``). They are
the only realistic way to parse Regina isosig strings without
reimplementing the entire isosig format. So this tiny Python module is
the *only* code that needs Regina; everything downstream consumes the
plain-text .bhrt format that the C++ core understands.

Usage from a shell:

    python3 -m bhrt_trisect.regina_bridge \
        data/Census/6p-closedOrientable.esig \
        data/Census/6p-closedOrientable.bhrt

Usage from Python:

    from bhrt_trisect.regina_bridge import esig_to_bhrt
    esig_to_bhrt("Census/6p-closedOrientable.esig",
                 "Census/6p-closedOrientable.bhrt")
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


def _regina():
    try:
        import regina  # type: ignore
        return regina
    except Exception as e:  # pragma: no cover
        raise RuntimeError(
            "Regina is not importable. Install Regina 7.4+ and ensure its "
            "Python bindings are on sys.path. See INSTALL.md section 3."
        ) from e


def _emit_bhrt(R, label: str, fh) -> None:
    if label:
        fh.write(f"# label: {label}\n")
    fh.write(f"N {R.size()}\n")
    for p in range(R.size()):
        regP = R.pentachoron(p)
        for f in range(5):
            adj = regP.adjacentPentachoron(f)
            if adj is None:
                continue
            if adj.index() < p:
                continue
            perm = regP.adjacentGluing(f)
            vs = " ".join(str(perm[i]) for i in range(5))
            fh.write(f"G {p} {f} {adj.index()} {vs}\n")
    fh.write("END\n")


def esig_to_bhrt(in_path: str, out_path: str, limit: int | None = None) -> int:
    """Convert an .esig file to .bhrt. Returns the record count."""
    regina = _regina()
    opener = gzip.open if str(in_path).endswith(".gz") else open
    count = 0
    with opener(in_path, "rt") as src, open(out_path, "w") as dst:
        dst.write("# bhrt-format v1\n")
        for line in src:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            sig = line.split("|", 1)[1] if "|" in line else line
            R = regina.Triangulation4.fromIsoSig(sig)
            if R is None:
                raise ValueError(f"failed to parse isosig: {sig}")
            _emit_bhrt(R, f"esig:{count}", dst)
            count += 1
            if limit is not None and count >= limit:
                break
    return count


def rga_to_bhrt(in_path: str, out_path: str) -> int:
    """Convert an .rga (Regina XML) file to .bhrt. Returns the record count."""
    regina = _regina()
    container = regina.open(in_path)
    if container is None:
        raise ValueError(f"Regina could not open {in_path}")
    count = 0
    with open(out_path, "w") as dst:
        dst.write("# bhrt-format v1\n")
        for child in (container.descendants() if hasattr(container, "descendants") else []):
            if isinstance(child, regina.Triangulation4):
                _emit_bhrt(child, getattr(child, "label", lambda: "")() or f"rga:{count}", dst)
                count += 1
    return count


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input", help="Regina .esig or .rga file")
    p.add_argument("output", help="output .bhrt file")
    p.add_argument("--limit", type=int, default=None)
    args = p.parse_args(argv)
    src = Path(args.input)
    if src.suffix in (".rga",):
        n = rga_to_bhrt(args.input, args.output)
    else:
        n = esig_to_bhrt(args.input, args.output, limit=args.limit)
    print(f"wrote {n} triangulations to {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
