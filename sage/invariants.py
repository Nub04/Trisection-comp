"""Sage-backed exact integer linear algebra for the invariants engine.

This module is imported by ``bhrt_trisect.invariants`` when the Sage
backend is selected. It uses Sage's Matrix(ZZ, ...) and Smith-normal-
form routines directly instead of the SymPy fallback.

Run inside Sage:

    sage -python -c "from bhrt_trisect.sage import invariants; invariants.smoke()"

Or inside a Sage notebook:

    load("bhrt_trisect/sage/invariants.py")
"""

from __future__ import annotations

from typing import List, Tuple

try:
    from sage.all import Matrix, ZZ, identity_matrix  # type: ignore
except Exception:  # pragma: no cover
    Matrix = None  # type: ignore
    ZZ = None  # type: ignore
    identity_matrix = None  # type: ignore


def smith_normal_form(M):
    """Return (D, U, V) such that U * M * V = D where D is in SNF.

    ``M`` is a Python list of lists of ints. The Sage call is exact.
    """
    if Matrix is None:
        raise RuntimeError("SageMath is not available")
    A = Matrix(ZZ, M)
    D, U, V = A.smith_form()
    return D, U, V


def symplectic_basis_for_intersection_form(J):
    """Return a basis change S so that S * J * S^T is symplectic.

    ``J`` is the integer-valued algebraic intersection matrix on
    H_1(Sigma; Z). Sage's symplectic_form_basis (when applicable)
    produces an exact integer basis change.
    """
    if Matrix is None:
        raise RuntimeError("SageMath is not available")
    A = Matrix(ZZ, J)
    return A.symplectic_form()


def smoke() -> None:
    """Quick check that the Sage backend imports correctly."""
    if Matrix is None:
        print("Sage not available; smoke test skipped.")
        return
    D, U, V = smith_normal_form([[2, 4, 4], [-6, 6, 12], [10, -4, -16]])
    print("Smith form diagonal:", [int(D[i, i]) for i in range(D.nrows())])
