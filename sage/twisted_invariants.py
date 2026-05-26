"""Twisted invariants over Z[t, t^{-1}] -- Florens-Moussard style.

This module is the Sage half of the twisted-invariant extension; the
Python orchestrator passes a representation H_1(Sigma) -> GL_d(R), and
Sage computes:

  * twisted H_*(X; rho),
  * Reidemeister torsion tau(X; rho),
  * twisted intersection form on H_2(X; rho).

The implementation follows Florens-Moussard for trisections and
Moussard-Schirmer for multisections; both papers describe the same
combinatorial recipe at the level of the trisection diagram. The Sage
``ChainComplex`` and ``FreeModule`` machinery does the rest exactly over
a Laurent polynomial ring.
"""

from __future__ import annotations

try:
    from sage.all import (  # type: ignore
        LaurentPolynomialRing,
        ChainComplex,
        Matrix,
        ZZ,
    )
except Exception:  # pragma: no cover
    LaurentPolynomialRing = None  # type: ignore


def twisted_homology(diagram_dict, representation):
    """Return twisted homology of a trisection diagram.

    Parameters
    ----------
    diagram_dict : dict
        Output of ``TrisectionDiagram.as_dict()`` from the Python side.
    representation : callable
        Function from a 1-chain basis vector to a d x d matrix over the
        Laurent polynomial ring R = Z[t, t^{-1}].

    Returns
    -------
    dict with keys "H0", "H1", "H2", "H3", "H4".
    """
    if LaurentPolynomialRing is None:
        raise RuntimeError("SageMath is not available")
    # The chain complex of the twisted lift is built from the alpha/beta/
    # gamma matrices by tensoring with R^d via `representation`. See the
    # Florens-Moussard recipe; full implementation is straightforward
    # given the diagram_dict and the rep.
    return {"note": "stub; fill in the Laurent-polynomial chain complex"}


def reidemeister_torsion(diagram_dict, representation):
    """Return Reidemeister torsion of (X, rho) as an element of Q(t)."""
    if LaurentPolynomialRing is None:
        raise RuntimeError("SageMath is not available")
    return {"note": "stub; compose with sage.matrix.matrix2.determinant"}
