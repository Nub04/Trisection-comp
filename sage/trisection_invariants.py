"""Exact (SageMath) invariants of a custom trisection diagram.

Reads the same `.tri` format as the C++ tool `bhrt-diagram` and computes the
homology, intersection form, signature and parity of the closed 4-manifold X
using Sage's exact integer linear algebra (Smith normal form, exact matrix
inverse). This is the exact-arithmetic cross-check of
`cpp/invariants.cpp::computeInvariantsFromClasses`.

Run:

    sage -python sage/trisection_invariants.py examples/cp2.tri
    sage -python sage/trisection_invariants.py --demo s2xs2

The `.tri` format (see examples/):

    g 1
    a 1 0
    b 0 1
    c 1 1
"""

from __future__ import annotations

import sys

try:
    from sage.all import Matrix, ZZ, QQ, identity_matrix, vector  # type: ignore
except Exception:  # pragma: no cover
    print("SageMath is required: run with `sage -python ...`.")
    raise


DEMOS = {
    "s4":    (0, [], [], []),
    "cp2":   (1, [[1, 0]], [[0, 1]], [[1, 1]]),
    "s1xs3": (1, [[1, 0]], [[1, 0]], [[1, 0]]),
    "s2xs2": (2, [[1, 0, 0, 0], [0, 1, 0, 0]],
                 [[0, 0, 1, 0], [0, 0, 0, 1]],
                 [[0, 1, 1, 0], [1, 0, 0, 1]]),
}


def parse_tri(path):
    g = None
    a, b, c = [], [], []
    with open(path) as fh:
        for raw in fh:
            line = raw.split("#", 1)[0].split()
            if not line:
                continue
            tag, rest = line[0], [int(x) for x in line[1:]]
            if tag == "g":
                g = rest[0]
            elif tag == "a":
                a.append(rest)
            elif tag == "b":
                b.append(rest)
            elif tag == "c":
                c.append(rest)
    if g is None:
        raise ValueError("no genus 'g' line")
    return g, a, b, c


def rank(rows, ncols):
    if not rows:
        return 0
    return Matrix(ZZ, rows).rank()


def cokernel(rows, ambient):
    """Free rank and torsion of Z^ambient / rowspace(rows)."""
    if not rows:
        return ambient, []
    D = Matrix(ZZ, rows).smith_form()[0]
    divs = [int(D[i, i]) for i in range(min(D.nrows(), D.ncols())) if D[i, i] != 0]
    free = ambient - len(divs)
    tor = [d for d in divs if d != 1]
    return free, tor


def invariants(g, a, b, c):
    twog = 2 * g
    out = {"g": g}
    if g == 0:
        out.update(b1=0, b2=0, b3=0, h1_tor=[], signature=0, parity="even", Q=None)
        return out

    A, B, C = Matrix(ZZ, a), Matrix(ZZ, b), Matrix(ZZ, c)

    stacked = A.stack(B).stack(C)
    b1, h1_tor = cokernel(stacked.rows(), twog)

    k_ab = twog - A.stack(B).rank()
    k_bc = twog - B.stack(C).rank()
    k_ac = twog - A.stack(C).rank()
    chi = 2 + g - (k_ab + k_bc + k_ac)
    b2 = chi - 2 + 2 * b1

    out.update(b1=b1, b2=max(b2, 0), b3=b1, h1_tor=h1_tor)

    # Linking-matrix intersection form for the simply-connected (g;0,0,0) case.
    out["signature"] = 0
    out["parity"] = "unknown"
    out["Q"] = None
    if b1 == 0 and k_ab == 0:
        basis = A.stack(B)                    # 2g x 2g, unimodular here
        if basis.is_invertible():
            coeff = C * basis.inverse()       # g x 2g over QQ
            M = coeff[:, 0:g]
            N = coeff[:, g:twog]
            if N.is_invertible():
                Q = (M * N.inverse())
                Qz = Matrix(ZZ, [[int(round(Q[i, j])) for j in range(g)]
                                 for i in range(g)])
                # symmetrise defensively
                Qz = (Qz + Qz.transpose()) / 2
                Qz = Matrix(ZZ, Qz)
                eigs = Qz.change_ring(QQ).eigenvalues()
                pos = sum(1 for e in eigs if e > 0)
                neg = sum(1 for e in eigs if e < 0)
                out["signature"] = pos - neg
                out["parity"] = "even" if all(Qz[i, i] % 2 == 0 for i in range(g)) else "odd"
                out["Q"] = Qz
    return out


def show(inv):
    def grp(free, tor):
        s = "Z^%d" % free
        for t in tor:
            s += " (+) Z/%d" % t
        return s
    print("genus g          =", inv["g"])
    print("H_1              =", grp(inv["b1"], inv["h1_tor"]))
    print("H_2              =", grp(inv["b2"], inv["h1_tor"]))
    print("H_3              =", grp(inv["b3"], []))
    print("b_2              =", inv["b2"])
    print("signature sigma  =", inv["signature"])
    print("parity           =", inv["parity"])
    if inv["Q"] is not None:
        print("intersection form Q_X =")
        print(inv["Q"])


def main(argv):
    if len(argv) >= 2 and argv[1] == "--demo":
        g, a, b, c = DEMOS[argv[2]]
    elif len(argv) >= 2:
        g, a, b, c = parse_tri(argv[1])
    else:
        print("usage: sage -python sage/trisection_invariants.py <file.tri>")
        print("       sage -python sage/trisection_invariants.py --demo {s4|cp2|s2xs2|s1xs3}")
        return 2
    show(invariants(g, a, b, c))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
