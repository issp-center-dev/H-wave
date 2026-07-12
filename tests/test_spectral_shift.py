"""Spectral-shift eigenvalue selection: a small positive leading eigenvalue
masked by larger-magnitude negative (repulsive) eigenvalues must still be
found. Plain which='LM' returns the large negatives and mis-reports the
leading lambda; spectral_shift asks ARPACK for the largest-REAL eigenvalue
(which='LR') on the right-half-plane-shifted operator A + sigma*I, recovering
the physical SC eigenvalue even for complex (non-Hermitian) spectra."""
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


# largest REAL part is +0.05, but two complex modes 0.03 +/- 0.6j have larger
# |lambda + sigma| -- a which='LM' shift would pick one of those; which='LR'
# must still return +0.05 (non-Hermitian / complex-spectrum case).
SPECTRUM_CX = [-0.5, 0.05, 0.03 + 0.6j, 0.03 - 0.6j, -0.2, -0.1]


def test_spectral_shift_complex_spectrum_picks_largest_real():
    mk, n = _make_operator_for(SPECTRUM_CX)
    lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                   spectral_shift="auto")
    assert lead.real == pytest.approx(0.05, abs=1e-6)
    assert abs(lead.imag) < 1e-6


def _make_operator_dense(M):
    """make_operator() -> (A, None) wrapping a dense (possibly non-normal) matrix."""
    M = np.asarray(M, dtype=complex)

    def make_operator():
        A = LinearOperator(M.shape, matvec=lambda v: M @ v, dtype=complex)
        return A, None
    return make_operator, M.shape[0]


def test_spectral_shift_nonnormal_complex_picks_largest_real():
    # upper-triangular => eigenvalues are the diagonal (SPECTRUM_CX), but the
    # off-diagonal makes the operator genuinely NON-NORMAL. Largest |lambda| is
    # a complex mode; which='LR' must still return the largest-real +0.05.
    M = np.diag(np.asarray(SPECTRUM_CX, dtype=complex))
    M[0, 1] = 5.0
    mk, n = _make_operator_dense(M)
    lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                   spectral_shift="auto")
    assert lead.real == pytest.approx(0.05, abs=1e-4)
    assert abs(lead.imag) < 1e-4


def test_spectral_shift_rejects_nonpositive_numeric():
    mk, n = _make_operator_for(SPECTRUM)
    for bad in (0.0, -1.0, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="spectral_shift"):
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                              spectral_shift=bad)


def test_spectral_shift_rejects_bad_string():
    mk, n = _make_operator_for(SPECTRUM)
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                          spectral_shift="invalid")


def test_spectral_shift_rejected_for_non_arnoldi():
    mk, n = _make_operator_for(SPECTRUM)
    with pytest.raises(ValueError, match="arnoldi"):
        sc._solve_leading(mk, n, "shift-invert", num_eigenvalues=3,
                          spectral_shift="auto")


def test_spectral_shift_rejected_for_subspace():
    # method="subspace" bypasses _solve_leading; _solve_eigenvalue must reject
    # a spectral_shift itself rather than silently ignoring it.
    with pytest.raises(ValueError, match="arnoldi"):
        sc._solve_eigenvalue(None, None, 1, 1, 1, 1, num_eigenvalues=3,
                             method="subspace", spectral_shift="auto")


def test_negative_leading_warns_for_unseeded_arnoldi(caplog):
    mk, n = _make_operator_for(SPECTRUM)
    with caplog.at_level("WARNING"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3)
    assert any("negative" in r.message for r in caplog.records)


def test_negative_leading_no_warn_with_seed_vec(caplog):
    # A seed vector selects the continuation branch (here the -0.5 mode);
    # that branch is intentionally tracked, so no misleading tip must fire.
    mk, n = _make_operator_for(SPECTRUM)
    seed = np.zeros(n, dtype=complex)
    seed[0] = 1.0  # overlaps the -0.5 eigenvector
    with caplog.at_level("WARNING"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=5, seed_vec=seed)
    assert not any("negative" in r.message for r in caplog.records)


def test_negative_leading_no_warn_for_roundoff_scale(caplog):
    # An all-but-zero spectrum whose leading real part is only -machine-eps
    # must not trigger the recommendation.
    mk, n = _make_operator_for([-1.0e-14, -2.0e-14, 1.0e-15, 3.0e-15])
    with caplog.at_level("WARNING"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=2)
    assert not any("negative" in r.message for r in caplog.records)

