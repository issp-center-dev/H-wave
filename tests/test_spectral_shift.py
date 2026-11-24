"""Spectral-shift eigenvalue selection: a small positive leading eigenvalue
masked by larger-magnitude negative (repulsive) eigenvalues must still be
found. Plain which='LM' returns the large negatives and mis-reports the
leading lambda; spectral_shift asks ARPACK for the largest-REAL eigenvalue
(which='LR') on the right-half-plane-shifted operator A + sigma*I, recovering
the physical SC eigenvalue even for complex (non-Hermitian) spectra."""
import unittest

import numpy as np
from scipy.sparse.linalg import LinearOperator
import hwave.sc as sc
from tests.approx_util import _caplog


def _make_operator_for(diag):
    """make_operator() -> (A, None) for a diagonal operator with given spectrum."""
    d = np.asarray(diag, dtype=complex)

    def make_operator():
        A = LinearOperator((d.size, d.size),
                           matvec=lambda v: d * v, dtype=complex)
        return A, None
    return make_operator, d.size


# ``_caplog`` now lives in ``tests.approx_util`` (shared with the other
# bond/spectral test modules).

# small positive leading (+0.05) hidden under larger negatives (down to -0.5)
SPECTRUM = [-0.5, -0.4, -0.3, -0.2, 0.05, 0.03]

# largest REAL part is +0.05, but two complex modes 0.03 +/- 0.6j have larger
# |lambda + sigma| -- a which='LM' shift would pick one of those; which='LR'
# must still return +0.05 (non-Hermitian / complex-spectrum case).
SPECTRUM_CX = [-0.5, 0.05, 0.03 + 0.6j, 0.03 - 0.6j, -0.2, -0.1]


class TestArnoldiSpectralShiftSelection(unittest.TestCase):

    def test_plain_lm_set_without_the_leading_is_resolved(self):
        """which='LM' with k=3 returns the three largest-magnitude (all
        negative) eigenvalues; their largest real part (-0.3) is below their
        smallest magnitude (0.3), so the LM set cannot contain the largest
        real part and the driver re-solves for it: +0.05, not a repulsive
        mode."""
        mk, n = _make_operator_for(SPECTRUM)
        lead, _, info = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3)
        np.testing.assert_allclose(lead.real, 0.05, rtol=0, atol=1e-6)
        assert info["selection"] == "LR_retry"

    def test_lm_set_that_certifies_its_leading_is_kept(self):
        """A LM set whose largest real part is at least its smallest
        magnitude already holds the global largest real part: no re-solve."""
        mk, n = _make_operator_for([0.9, -0.5, 0.4, -0.1, 0.05, 0.03])
        lead, _, info = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3)
        np.testing.assert_allclose(lead.real, 0.9, rtol=0, atol=1e-8)
        assert info["selection"] == "LM"

    def test_spectral_shift_auto_recovers_positive_leading(self):
        mk, n = _make_operator_for(SPECTRUM)
        lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                       spectral_shift="auto")
        np.testing.assert_allclose(lead.real, 0.05, rtol=0, atol=1e-6)

    def test_spectral_shift_numeric_recovers_positive_leading(self):
        mk, n = _make_operator_for(SPECTRUM)
        lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                       spectral_shift=1.0)
        np.testing.assert_allclose(lead.real, 0.05, rtol=0, atol=1e-6)

    def test_spectral_shift_complex_spectrum_picks_largest_real(self):
        mk, n = _make_operator_for(SPECTRUM_CX)
        lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                       spectral_shift="auto")
        np.testing.assert_allclose(lead.real, 0.05, rtol=0, atol=1e-6)
        assert abs(lead.imag) < 1e-6

    def test_spectral_shift_nonnormal_complex_picks_largest_real(self):
        # upper-triangular => eigenvalues are the diagonal (SPECTRUM_CX), but the
        # off-diagonal makes the operator genuinely NON-NORMAL. Largest |lambda| is
        # a complex mode; which='LR' must still return the largest-real +0.05.
        M = np.diag(np.asarray(SPECTRUM_CX, dtype=complex))
        M[0, 1] = 5.0
        mk, n = self._make_operator_dense(M)
        lead, _, _ = sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                       spectral_shift="auto")
        np.testing.assert_allclose(lead.real, 0.05, rtol=0, atol=1e-4)
        assert abs(lead.imag) < 1e-4

    @staticmethod
    def _make_operator_dense(M):
        """make_operator() -> (A, None) wrapping a dense (possibly non-normal) matrix."""
        M = np.asarray(M, dtype=complex)

        def make_operator():
            A = LinearOperator(M.shape, matvec=lambda v: M @ v, dtype=complex)
            return A, None
        return make_operator, M.shape[0]

    def test_spectral_shift_rejects_nonpositive_numeric(self):
        mk, n = _make_operator_for(SPECTRUM)
        for bad in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "spectral_shift"):
                    sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                                      spectral_shift=bad)

    def test_spectral_shift_rejects_bad_string(self):
        mk, n = _make_operator_for(SPECTRUM)
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                              spectral_shift="invalid")

    def test_spectral_shift_rejects_bool_true_on_both_validation_paths(self):
        """FIX 2 (review round 3): ``float(True) == 1.0`` succeeds, so a TOML
        ``spectral_shift = true`` used to be silently accepted as sigma=1.0
        instead of being rejected as a type error. Must be rejected for BOTH
        solver_mode values ``_validate_spectral_shift`` otherwise accepts a
        numeric shift for."""
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._validate_spectral_shift(True, "arnoldi")
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._validate_spectral_shift(False, "arnoldi")
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._validate_spectral_shift(True, "iteration")

    def test_spectral_shift_rejects_bool_true_via_solve_leading_arnoldi(self):
        mk, n = _make_operator_for(SPECTRUM)
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3,
                              spectral_shift=True)

    def test_spectral_shift_rejected_for_non_arnoldi(self):
        mk, n = _make_operator_for(SPECTRUM)
        with self.assertRaisesRegex(ValueError, "arnoldi"):
            sc._solve_leading(mk, n, "shift-invert", num_eigenvalues=3,
                              spectral_shift="auto")

    def test_spectral_shift_rejected_for_subspace(self):
        # method="subspace" bypasses _solve_leading; _solve_eigenvalue must reject
        # a spectral_shift itself rather than silently ignoring it.
        with self.assertRaisesRegex(ValueError, "arnoldi"):
            sc._solve_eigenvalue(None, None, 1, 1, 1, 1, num_eigenvalues=3,
                                 method="subspace", spectral_shift="auto")

    def test_uncertified_lm_set_warns_for_unseeded_arnoldi(self):
        mk, n = _make_operator_for(SPECTRUM)
        with _caplog("WARNING") as records:
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=3)
        assert any("largest real part" in r.getMessage() for r in records)

    def test_negative_leading_no_warn_with_seed_vec(self):
        # A seed vector selects the continuation branch (here the -0.5 mode);
        # that branch is intentionally tracked, so no misleading tip must fire.
        mk, n = _make_operator_for(SPECTRUM)
        seed = np.zeros(n, dtype=complex)
        seed[0] = 1.0  # overlaps the -0.5 eigenvector
        with _caplog("WARNING") as records:
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=5, seed_vec=seed)
        assert not any("negative" in r.message for r in records)


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


class TestPowerIterationSpectralShift(unittest.TestCase):

    def test_iteration_without_shift_does_not_converge_on_negative_dominant(self):
        lead, _, info = _iterate(NEG_DOMINANT, max_iter=200)
        assert info["converged"] is False

    def test_iteration_numeric_shift_converges_and_unshifts(self):
        """K + 4I = diag(1, 5, 4.5) has a positive dominant eigenvalue, so the
        power iteration converges; subtracting the shift back must report the
        eigenvalue of the ORIGINAL operator (+1, the algebraically largest -- the
        physical Eliashberg mode, same selection as arnoldi's which='LR'), with
        the matching (un-shifted) eigenvector."""
        lead, vec, info = _iterate(NEG_DOMINANT, spectral_shift=4.0)
        assert info["converged"] is True
        np.testing.assert_allclose(lead, 1.0, rtol=0, atol=1e-6)
        # eigenvector of the lambda=+1 mode (index 1), unchanged by the shift
        np.testing.assert_allclose(abs(vec[1]), 1.0, rtol=0, atol=1e-6)
        assert abs(vec[0]) < 1e-6 and abs(vec[2]) < 1e-6

    def test_iteration_shift_reports_the_negative_lambda_with_its_sign(self):
        """All-negative spectrum: the reported value must be the SIGNED lambda
        (-3), not the shifted eigenvalue (+7) and not the unsigned iterate norm."""
        lead, vec, info = _iterate(ALL_NEGATIVE, spectral_shift=10.0)
        assert info["converged"] is True
        np.testing.assert_allclose(lead, -3.0, rtol=0, atol=1e-6)
        np.testing.assert_allclose(abs(vec[0]), 1.0, rtol=0, atol=1e-6)

    def test_iteration_auto_shift_reports_the_negative_lambda(self):
        lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift="auto")
        assert info["converged"] is True
        np.testing.assert_allclose(lead, -3.0, rtol=0, atol=1e-6)

    def test_iteration_auto_shift_on_negative_dominant(self):
        lead, _, info = _iterate(NEG_DOMINANT, spectral_shift="auto")
        assert info["converged"] is True
        np.testing.assert_allclose(lead, 1.0, rtol=0, atol=1e-6)

    def test_iteration_shift_is_a_noop_on_a_positive_dominant_operator(self):
        """Guards the subtract-back: on an operator the plain iteration already
        handles, the shift must not change the answer."""
        lead_plain, _, info_plain = _iterate([2.0, 1.0, 0.5])
        lead_shift, _, info_shift = _iterate([2.0, 1.0, 0.5], spectral_shift=5.0)
        lead_auto, _, info_auto = _iterate([2.0, 1.0, 0.5], spectral_shift="auto")
        assert info_plain["converged"] and info_shift["converged"] and \
            info_auto["converged"]
        np.testing.assert_allclose(lead_shift, lead_plain, rtol=0, atol=1e-6)
        np.testing.assert_allclose(lead_auto, lead_plain, rtol=0, atol=1e-6)
        np.testing.assert_allclose(lead_plain, 2.0, rtol=0, atol=1e-6)

    def test_iteration_shift_respects_the_parity_projection(self):
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
        np.testing.assert_allclose(lead, -3.0, rtol=0, atol=1e-6)
        assert abs(vec[1]) < 1e-12          # never leaves the sector
        # sector-restricted radius (5) => shift ~5.25, well below 1.05 * 8
        assert info["spectral_shift"] < 6.0

    def test_iteration_shift_logs_sigma_and_the_unshift(self):
        with _caplog("INFO", logger="hwave_sc") as records:
            _iterate(ALL_NEGATIVE, spectral_shift=10.0)
        msgs = " ".join(r.getMessage() for r in records)
        assert "spectral shift" in msgs.lower()
        assert "10.0" in msgs or "10.000000" in msgs
        assert "un-shift" in msgs.lower() or "unshift" in msgs.lower()

    def test_iteration_rejects_bad_spectral_shift_values(self):
        for bad in (0.0, -1.0, float("nan"), float("inf"), "invalid"):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "spectral_shift"):
                    _iterate(NEG_DOMINANT, spectral_shift=bad)

    def test_iteration_rejects_bool_true_as_numeric_spectral_shift(self):
        """FIX 2 (review round 3), iteration path: same bool-as-1.0 hole must be
        closed on ``solver_mode='iteration'`` too, not just arnoldi."""
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            _iterate(NEG_DOMINANT, spectral_shift=True)

    def test_spectral_shift_still_rejected_for_shift_invert(self):
        mk, n = _make_operator_for(SPECTRUM)
        with self.assertRaisesRegex(ValueError, "spectral_shift"):
            sc._solve_leading(mk, n, "shift-invert-gmres", num_eigenvalues=3,
                              spectral_shift="auto")

    def test_negative_leading_no_warn_for_roundoff_scale(self):
        # An all-but-zero spectrum whose leading real part is only -machine-eps
        # must not trigger the recommendation.
        mk, n = _make_operator_for([-1.0e-14, -2.0e-14, 1.0e-15, 3.0e-15])
        with _caplog("WARNING") as records:
            sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=2)
        assert not any("negative" in r.message for r in records)


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

class TestShiftedEigenvalueValidation(unittest.TestCase):

    def test_shifted_run_reports_a_validated_rayleigh_eigenvalue(self):
        """(b) A sufficient shift (sigma = 10 on ALL_NEGATIVE) still yields
        lambda = -3, now backed by a small Rayleigh residual and flagged as a
        genuine eigenvalue."""
        lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=10.0)
        assert info["converged"] is True
        np.testing.assert_allclose(lead, -3.0, rtol=0, atol=1e-6)
        assert info["eigenvalue_validated"] is True
        assert info["eigenvalue_kind"] == "eigenvalue"
        np.testing.assert_allclose(info["rayleigh_eigenvalue"], -3.0,
                                   rtol=0, atol=1e-6)
        assert info["rayleigh_residual"] < 1.0e-6
        assert info["shift_sufficient"] is True

    def test_insufficient_shift_never_reports_the_artefact_as_an_eigenvalue(self):
        """(a) diag(-3, -8, -5) with sigma = 1 and a short max_iter: the raw
        quantity ||(K + sigma I) v|| - sigma is the ~6 artefact, which is not
        an eigenvalue of K. A non-converged run reports NaN, never it."""
        lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=3,
                                 alpha=0.0)
        assert info["converged"] is False
        assert np.isnan(lead)
        assert info["eigenvalue_validated"] is False
        assert info["eigenvalue_kind"] == "not-converged"
        assert info["rayleigh_eigenvalue"] is None

    def test_insufficient_shift_that_locks_on_is_a_two_cycle(self):
        """Run the same insufficient shift long enough for the iterate to lock
        onto the -8 mode: its shifted eigenvalue -7 is negative, so the iterate
        alternates sign -- a two-cycle, recorded and reported as NaN at
        max_iter (no eigenvalue); the cycling mode's true eigenvalue (-8) is
        recorded."""
        lead, _, info = _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=200,
                                 alpha=0.0)
        assert np.isnan(lead)
        assert info["converged"] is False
        assert info["two_cycle"] is True
        # recorded, not stopped on: the run goes on to max_iter
        assert info["n_iter"] == 200
        np.testing.assert_allclose(info["cycle_eigenvalue"], -8.0, rtol=0,
                                   atol=1e-6)

    def test_insufficient_shift_warns(self):
        with _caplog("WARNING", logger="hwave_sc") as records:
            _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=200, alpha=0.0)
        msgs = " ".join(r.getMessage() for r in records)
        assert "insufficient" in msgs.lower()

    def test_unconverged_shifted_run_warns_that_no_eigenvalue_is_reported(self):
        with _caplog("WARNING", logger="hwave_sc") as records:
            _iterate(ALL_NEGATIVE, spectral_shift=1.0, max_iter=3, alpha=0.0)
        msgs = " ".join(r.getMessage() for r in records)
        assert "no eigenvalue is reported" in msgs.lower()

    def test_unshifted_iteration_reports_the_signed_rayleigh_eigenvalue(self):
        """The unshifted path (spectral_shift None / "none") also reports the
        SIGNED Rayleigh eigenvalue of its converged iterate, never the iterate
        norm, and carries no spectral_shift key."""
        lead, _, info = _iterate([2.0, 1.0, 0.5])
        assert info["converged"] is True
        np.testing.assert_allclose(lead, 2.0, rtol=0, atol=1e-8)
        np.testing.assert_allclose(info["rayleigh_eigenvalue"], 2.0, rtol=0,
                                   atol=1e-8)
        assert "spectral_shift" not in info
        assert info["two_cycle"] is False


# --- the note that reaches eigenvalue.dat -----------------------------------

class TestShiftedEigenvalueNote(unittest.TestCase):

    def test_note_labels_a_validated_shifted_value_as_an_eigenvalue(self):
        info = {"eigenvalue_validated": True, "rayleigh_residual": 1.2e-12,
                "eigenvalue_kind": "eigenvalue"}
        note = sc._shifted_eigenvalue_note("iteration", "auto", True, 12, info)
        assert "SIGNED eigenvalue" in note
        assert "converged" in note

    def test_note_labels_a_non_converged_value_as_not_an_eigenvalue(self):
        info = {"eigenvalue_validated": False, "rayleigh_residual": 0.55,
                "eigenvalue_kind": "not-converged"}
        note = sc._shifted_eigenvalue_note("iteration", 1.0, False, 3, info)
        lower = note.lower()
        assert "not an eigenvalue" in lower
        assert "nan" in lower
        assert "eigen_converged = false" in lower
        # (c) the non-converged status is stated too
        assert "did not converge" in lower

    def test_note_names_the_two_cycle(self):
        info = {"two_cycle": True, "cycle_eigenvalue": -1.93}
        note = sc._shifted_eigenvalue_note("iteration", None, False, 14, info)
        assert "two-cycle" in note
        assert "NEGATIVE" in note
        assert "no spectral_shift" in note


class TestShiftedIterationIsUnmixed(unittest.TestCase):
    """With a spectral shift (auto or explicit) the power iteration does not
    mix the old iterate in (``alpha`` would act as a second shift); without a
    shift ``alpha`` is applied as before."""

    @staticmethod
    def _op(seed=11, n=40):
        rng = np.random.default_rng(seed)
        Q, _ = np.linalg.qr(rng.standard_normal((n, n)))
        lam = np.concatenate([[0.9, 0.8], rng.uniform(-1.5, 0.6, n - 2)])
        M = (Q * lam) @ Q.T
        ev = np.linalg.eigvals(M)

        def mk():
            return LinearOperator((n, n), matvec=lambda v: M @ v,
                                  dtype=complex), n
        init = rng.standard_normal(n) + 0j
        return mk, n, init, float(np.max(ev.real))

    def _run(self, spectral_shift, alpha):
        mk, n, init, lead = self._op()
        lam, _vec, info = sc._solve_leading(
            mk, n, "iteration", init_vec=init, alpha=alpha, max_iter=5000,
            convergence_tol=1e-9, spectral_shift=spectral_shift)
        return lam, info, lead

    def test_auto_shift_ignores_alpha(self):
        lam_a, info_a, lead = self._run("auto", 0.5)
        lam_0, info_0, _ = self._run("auto", 0.0)
        self.assertTrue(info_a["converged"])
        self.assertAlmostEqual(lam_a, lead, delta=1e-7)
        self.assertEqual(info_a["n_iter"], info_0["n_iter"])
        self.assertEqual(info_a["alpha_applied"], 0.0)

    def test_explicit_shift_ignores_alpha(self):
        lam_a, info_a, lead = self._run(2.0, 0.5)
        _lam, info_0, _ = self._run(2.0, 0.0)
        self.assertTrue(info_a["converged"])
        self.assertAlmostEqual(lam_a, lead, delta=1e-7)
        self.assertEqual(info_a["n_iter"], info_0["n_iter"])

    def test_unshifted_iteration_keeps_alpha(self):
        mk, n, _init, _lead = self._op()
        d = np.array([1.0, 0.6, 0.3, 0.1] + [0.05] * 6, complex)

        def mkd():
            return LinearOperator((10, 10), matvec=lambda v: d * v,
                                  dtype=complex), 10
        init = np.ones(10, complex)
        runs = {}
        for alpha in (0.5, 0.0):
            lam, _v, info = sc._solve_leading(
                mkd, 10, "iteration", init_vec=init, alpha=alpha,
                max_iter=5000, convergence_tol=1e-9, spectral_shift=None)
            self.assertTrue(info["converged"])
            self.assertAlmostEqual(lam, 1.0, delta=1e-7)
            runs[alpha] = info
        self.assertEqual(runs[0.5]["alpha_applied"], 0.5)
        self.assertGreater(runs[0.5]["n_iter"], runs[0.0]["n_iter"])


if __name__ == "__main__":
    unittest.main()
