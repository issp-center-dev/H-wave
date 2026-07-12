"""Spectral-shift eigenvalue selection: a small positive leading eigenvalue
masked by larger-magnitude negative (repulsive) eigenvalues must still be
found. Plain which='LM' returns the large negatives and mis-reports the
leading lambda; spectral_shift='auto' shifts A -> A + sigma*I so 'LM' maps to
the largest-REAL eigenvalue (the physical SC eigenvalue)."""
import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator
import hwave.sc as sc


def _make_operator_for(diag):
    """make_operator() -> (A, None) for a diagonal operator with given spectrum."""
    d = np.asarray(diag, dtype=complex)

    def make_operator():
        A = LinearOperator((d.size, d.size),
                           matvec=lambda v: d * v, dtype=complex)
        return A, None
    return make_operator, d.size


# small positive leading (+0.05) hidden under larger negatives (down to -0.5)
SPECTRUM = [-0.5, -0.4, -0.3, -0.2, 0.05, 0.03]


def test_plain_lm_misses_small_positive_leading():
    mk, n = _make_operator_for(SPECTRUM)
    lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3)
    # which='LM' with k=3 returns the three largest-magnitude (all negative),
    # so the reported leading is a negative repulsive mode, NOT +0.05.
    assert lead.real < 0.0


def test_spectral_shift_auto_recovers_positive_leading():
    mk, n = _make_operator_for(SPECTRUM)
    lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                   spectral_shift="auto")
    assert lead.real == pytest.approx(0.05, abs=1e-6)


def test_spectral_shift_numeric_recovers_positive_leading():
    mk, n = _make_operator_for(SPECTRUM)
    lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                   spectral_shift=1.0)
    assert lead.real == pytest.approx(0.05, abs=1e-6)

