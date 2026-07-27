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


def test_spectral_shift_rejects_bool_true_on_both_validation_paths():
    """FIX 2 (review round 3): ``float(True) == 1.0`` succeeds, so a TOML
    ``spectral_shift = true`` used to be silently accepted as sigma=1.0
    instead of being rejected as a type error. Must be rejected for BOTH
    solver_mode values ``_validate_spectral_shift`` otherwise accepts a
    numeric shift for."""
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._validate_spectral_shift(True, "arnoldi")
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._validate_spectral_shift(False, "arnoldi")
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._validate_spectral_shift(True, "iteration")


def test_spectral_shift_rejects_bool_true_via_solve_leading_arnoldi():
    mk, n = _make_operator_for(SPECTRUM)
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                          spectral_shift=True)


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


# ---------------------------------------------------------------------------
# spectral_shift on the POWER-ITERATION path (solver_mode="iteration")
#
# A repulsive-dominant kernel (dominant eigenvalue negative) makes the plain
# power iteration flip the iterate's sign every step, so its convergence test
# ||sigma_new/||sigma_new|| - sigma_old|| sits at ~2 forever. Iterating on the
# shifted operator K + sigma*I (sigma above the spectral radius) makes the
# ALGEBRAICALLY LARGEST eigenvalue the dominant one -- and positive, so the
# iteration converges -- and sigma is subtracted back from the reported
# eigenvalue. The eigenvector is untouched by the shift.
# ---------------------------------------------------------------------------

_ITER_KW = dict(max_iter=4000, convergence_tol=1.0e-8, alpha=0.5)


def _iterate(diag, init=None, **kw):
    mk, n = _make_operator_for(diag)
    if init is None:
        # generic start: overlaps every eigenvector of a diagonal operator
        init = np.ones(n, dtype=complex) / np.sqrt(n)
    kwargs = dict(_ITER_KW)
    kwargs.update(kw)
    return sc._solve_leading(mk, n, "iteration", init_vec=init, **kwargs)


# dominant (largest-|lambda|) eigenvalue is NEGATIVE => plain power iteration
# cannot converge; the algebraically largest (physical) eigenvalue is +1.
NEG_DOMINANT = [-3.0, 1.0, 0.5]

# every eigenvalue negative: the algebraically largest is -3, and the reported
# lambda must carry that MINUS sign (not the unsigned iterate norm 8, and not
# the shifted value).
ALL_NEGATIVE = [-3.0, -8.0, -5.0]


def test_iteration_without_shift_does_not_converge_on_negative_dominant():
    lead, _, info = _iterate(NEG_DOMINANT, max_iter=200)
    assert info["converged"] is False


def test_iteration_numeric_shift_converges_and_unshifts():
    """K + 4I = diag(1, 5, 4.5) has a positive dominant eigenvalue, so the
    power iteration converges; subtracting the shift back must report the
    eigenvalue of the ORIGINAL operator (+1, the algebraically largest -- the
    physical Eliashberg mode, same selection as arnoldi's which='LR'), with
    the matching (un-shifted) eigenvector."""
    lead, vec, info = _iterate(NEG_DOMINANT, spectral_shift=4.0)
    assert info["converged"] is True
    assert lead == pytest.approx(1.0, abs=1e-6)
    # eigenvector of the lambda=+1 mode (index 1), unchanged by the shift
    assert abs(vec[1]) == pytest.approx(1.0, abs=1e-6)
    assert abs(vec[0]) < 1e-6 and abs(vec[2]) < 1e-6


def test_iteration_shift_reports_the_negative_lambda_with_its_sign():
    """All-negative spectrum: the reported value must be the SIGNED lambda
    (-3), not the shifted eigenvalue (+7) and not the unsigned iterate norm."""
    lead, vec, info = _iterate(ALL_NEGATIVE, spectral_shift=10.0)
    assert info["converged"] is True
    assert lead == pytest.approx(-3.0, abs=1e-6)
    assert abs(vec[0]) == pytest.approx(1.0, abs=1e-6)


def test_iteration_auto_shift_reports_the_negative_lambda():
    lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift="auto")
    assert info["converged"] is True
    assert lead == pytest.approx(-3.0, abs=1e-6)


def test_iteration_auto_shift_on_negative_dominant():
    lead, _, info = _iterate(NEG_DOMINANT, spectral_shift="auto")
    assert info["converged"] is True
    assert lead == pytest.approx(1.0, abs=1e-6)


def test_iteration_shift_is_a_noop_on_a_positive_dominant_operator():
    """Guards the subtract-back: on an operator the plain iteration already
    handles, the shift must not change the answer."""
    lead_plain, _, info_plain = _iterate([2.0, 1.0, 0.5])
    lead_shift, _, info_shift = _iterate([2.0, 1.0, 0.5], spectral_shift=5.0)
    lead_auto, _, info_auto = _iterate([2.0, 1.0, 0.5], spectral_shift="auto")
    assert info_plain["converged"] and info_shift["converged"] and \
        info_auto["converged"]
    assert lead_shift == pytest.approx(lead_plain, abs=1e-6)
    assert lead_auto == pytest.approx(lead_plain, abs=1e-6)
    assert lead_plain == pytest.approx(2.0, abs=1e-6)


def test_iteration_shift_respects_the_parity_projection():
    """The shift preserves the invariant subspaces, so a projected iteration
    must stay in its sector -- and the AUTO shift must be estimated inside that
    sector (component 1, the -8 mode, is projected out, so the sector's
    spectral radius is 5, not 8)."""
    mk, n = _make_operator_for(ALL_NEGATIVE)
    keep = np.array([1.0, 0.0, 1.0])

    def project_fn(v):
        return keep * v

    init = project_fn(np.ones(n, dtype=complex) / np.sqrt(n))
    lead, vec, info = sc._solve_leading(mk, n, "iteration", init_vec=init,
                                        project_fn=project_fn,
                                        spectral_shift="auto", **_ITER_KW)
    assert info["converged"] is True
    assert lead == pytest.approx(-3.0, abs=1e-6)
    assert abs(vec[1]) < 1e-12          # never leaves the sector
    # sector-restricted radius (5) => shift ~5.25, well below 1.05 * 8
    assert info["spectral_shift"] < 6.0


def test_iteration_shift_logs_sigma_and_the_unshift(caplog):
    with caplog.at_level("INFO", logger="hwave_sc"):
        _iterate(ALL_NEGATIVE, spectral_shift=10.0)
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "spectral shift" in msgs.lower()
    assert "10.0" in msgs or "10.000000" in msgs
    assert "un-shift" in msgs.lower() or "unshift" in msgs.lower()


def test_iteration_rejects_bad_spectral_shift_values():
    for bad in (0.0, -1.0, float("nan"), float("inf"), "invalid"):
        with pytest.raises(ValueError, match="spectral_shift"):
            _iterate(NEG_DOMINANT, spectral_shift=bad)


def test_iteration_rejects_bool_true_as_numeric_spectral_shift():
    """FIX 2 (review round 3), iteration path: same bool-as-1.0 hole must be
    closed on ``solver_mode='iteration'`` too, not just arnoldi."""
    with pytest.raises(ValueError, match="spectral_shift"):
        _iterate(NEG_DOMINANT, spectral_shift=True)


def test_spectral_shift_still_rejected_for_shift_invert():
    mk, n = _make_operator_for(SPECTRUM)
    with pytest.raises(ValueError, match="spectral_shift"):
        sc._solve_leading(mk, n, "shift-invert-gmres", num_eigenvalues=3,
                          spectral_shift="auto")


def test_negative_leading_no_warn_for_roundoff_scale(caplog):
    # An all-but-zero spectrum whose leading real part is only -machine-eps
    # must not trigger the recommendation.
    mk, n = _make_operator_for([-1.0e-14, -2.0e-14, 1.0e-15, 3.0e-15])
    with caplog.at_level("WARNING"):
        sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=2)
    assert not any("negative" in r.message for r in caplog.records)



# ---------------------------------------------------------------------------
# The shifted power iteration must not claim "eigenvalue" for a number that
# is not one.
#
# ``||(K + sigma I) v|| - sigma`` equals a signed eigenvalue of K only after
# the iteration has converged to a mode whose SHIFTED eigenvalue is positive
# real. An explicit spectral_shift is not required to be that large, and
# "auto" is only a finite-step estimate of the spectral radius -- e.g.
# diag(-3, -8, -5) with sigma = 1 drives the iterate toward the -8 mode, whose
# shifted eigenvalue is -7, and the raw quantity is |-7| - 1 = +6, which is
# not an eigenvalue of K at all.
#
# The solver therefore validates the returned vector with a signed Rayleigh
# quotient on the UNSHIFTED operator and only reports an eigenvalue when the
# residual says so; otherwise it labels the number as a shifted iterate-norm
# ESTIMATE.
# ---------------------------------------------------------------------------

def test_shifted_run_reports_a_validated_rayleigh_eigenvalue():
    """(b) A sufficient shift (sigma = 10 on ALL_NEGATIVE) still yields
    lambda = -3, now backed by a small Rayleigh residual and flagged as a
    genuine eigenvalue."""
    lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=10.0)
    assert info["converged"] is True
    assert lead == pytest.approx(-3.0, abs=1e-6)
    assert info["eigenvalue_validated"] is True
    assert info["eigenvalue_kind"] == "eigenvalue"
    assert info["rayleigh_eigenvalue"] == pytest.approx(-3.0, abs=1e-6)
    assert info["rayleigh_residual"] < 1.0e-6
    assert info["shift_sufficient"] is True


def test_insufficient_shift_never_reports_the_artefact_as_an_eigenvalue():
    """(a) diag(-3, -8, -5) with sigma = 1 and a short max_iter: the raw
    quantity is the ~6 artefact, which is not an eigenvalue of K. It must be
    labelled a shifted iterate-norm estimate, never an eigenvalue."""
    lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=3,
                             alpha=0.0)
    assert info["converged"] is False
    # the artefact itself (documented, not endorsed)
    assert 5.0 < lead < 6.5
    assert info["eigenvalue_validated"] is False
    assert info["eigenvalue_kind"] == "shifted-iterate-norm-estimate"
    assert info["rayleigh_eigenvalue"] is None
    assert info["rayleigh_residual"] > info["rayleigh_residual_tol"]
    # ... and no eigenvalue of K is anywhere near the reported number
    assert min(abs(lead - ev) for ev in ALL_NEGATIVE) > 1.0


def test_insufficient_shift_that_does_lock_on_reports_the_true_eigenvalue():
    """Run the same insufficient shift long enough for the iterate to lock
    onto the -8 mode: the reported number must then be that TRUE eigenvalue
    (-8), never the 6 artefact, and the run must warn that sigma was too
    small (the mode's shifted eigenvalue -7 is negative)."""
    lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=200,
                             alpha=0.0)
    assert lead == pytest.approx(-8.0, abs=1e-6)
    assert info["eigenvalue_validated"] is True
    assert info["shift_sufficient"] is False


def test_insufficient_shift_warns(caplog):
    with caplog.at_level("WARNING", logger="hwave_sc"):
        _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=200, alpha=0.0)
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "insufficient" in msgs.lower()


def test_unvalidated_shifted_run_warns_that_it_is_not_an_eigenvalue(caplog):
    with caplog.at_level("WARNING", logger="hwave_sc"):
        _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=3, alpha=0.0)
    msgs = " ".join(r.getMessage() for r in caplog.records)
    assert "not an eigenvalue" in msgs.lower()


def test_unshifted_iteration_info_is_unchanged():
    """The default (unshifted) path must be untouched: no Rayleigh matvec, no
    extra keys, exactly the historical info dict.

    The historical contract (commit 712b1a0, ``_solve_iteration``) is
    ``{"converged": bool, "n_iter": int}`` and NOTHING else -- not even
    ``spectral_shift: None``. Callers that compare or serialize this dict must
    see the same key set they always saw when no shift was requested.
    """
    _, _, info = _iterate([2.0, 1.0, 0.5])
    assert set(info) == {"converged", "n_iter"}
    assert "spectral_shift" not in info


# --- the note that reaches eigenvalue.dat -----------------------------------

def test_note_labels_a_validated_shifted_value_as_an_eigenvalue():
    info = {"eigenvalue_validated": True, "rayleigh_residual": 1.2e-12,
            "eigenvalue_kind": "eigenvalue"}
    note = sc._shifted_eigenvalue_note("iteration", "auto", True, 12, info)
    assert "SIGNED eigenvalue" in note
    assert "converged" in note


def test_note_labels_an_unvalidated_shifted_value_as_not_an_eigenvalue():
    info = {"eigenvalue_validated": False, "rayleigh_residual": 0.55,
            "eigenvalue_kind": "shifted-iterate-norm-estimate"}
    note = sc._shifted_eigenvalue_note("iteration", 1.0, False, 3, info)
    lower = note.lower()
    assert "not an eigenvalue" in lower
    assert "estimate" in lower
    # (c) the non-converged status is stated too
    assert "did not converge" in lower
