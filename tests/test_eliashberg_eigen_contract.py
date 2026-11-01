"""What the Eliashberg eigen-solvers may report as "the leading eigenvalue".

A repulsive-dominant pairing kernel has a NEGATIVE eigenvalue of the largest
magnitude (the on-site s-wave mode of the repulsive Hubbard model, about
-1.9 on an 8x8 U = 3 grid) next to the small positive eigenvalue of the
physical (d-wave) channel. The damped power iteration
``sigma <- (1 - alpha) A sigma / ||A sigma|| + alpha sigma`` then settles into
a sign-alternating two-cycle of amplitude ``s = (1 - alpha) / (1 + alpha)``
and never converges; its last iterate norm ``|lambda_-| s`` is not an
eigenvalue of anything.

The contract pinned here, on a synthetic kernel with exactly that spectrum
(fast, no FLEX run):

* a power iteration that does not converge -- by the iteration cap or by the
  detected two-cycle -- reports ``lambda = NaN`` and ``eigen_converged =
  False``, and the two-cycle is detected long before the cap;
* a converged iteration reports the SIGNED Rayleigh quotient of its iterate,
  equal to the dense/ARPACK eigenvalue;
* without an explicit ``spectral_shift`` the iteration runs with
  ``spectral_shift = "auto"`` and converges to the largest-REAL-part
  eigenvalue; ``spectral_shift = "none"`` keeps the unshifted iteration;
* every eigenvalue-family path reports the largest-real-part eigenvalue even
  when the requested largest-magnitude set does not contain it;
* the writers carry ``eigen_converged`` into ``eigenvalue.dat`` and the npz.

The end-to-end FLEX checks live in ``test_eliashberg_bond_repulsive_singlet``.
"""
import os
import shutil
import tempfile
import types
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import eliashberg_dynamic as ed


LAM_NEG = -1.93          # largest magnitude, repulsive (on-site s-like)
LAM_POS = 0.359          # largest real part, the physical channel eigenvalue
LAM_REST = 0.01          # every other mode

GAP_SHAPE = (1, 1, 4, 4, 1, 8)


def _channel_pair(shape=GAP_SHAPE, seed=7):
    """Two orthonormal states of the singlet even-frequency channel."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(2):
        g = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
        v = ed._project_parity_dynamic(g, "singlet").ravel()
        for u in out:
            v = v - np.vdot(u, v) * u
        out.append(v / np.linalg.norm(v))
    return out


def _repulsive_matvec(shape=GAP_SHAPE):
    """Kernel with spectrum {LAM_NEG, LAM_POS, LAM_REST...}; both special
    modes lie in the singlet channel, so channel projection cannot separate
    them -- only the eigenvalue selection can."""
    a, b = _channel_pair(shape)

    def matvec(x):
        x = np.asarray(x).ravel()
        ca, cb = np.vdot(a, x), np.vdot(b, x)
        rest = x - ca * a - cb * b
        return LAM_NEG * ca * a + LAM_POS * cb * b + LAM_REST * rest

    return matvec, a, b


def _seed(shape=GAP_SHAPE):
    """A channel seed with a component on both special modes."""
    a, b = _channel_pair(shape)
    rng = np.random.default_rng(3)
    g = rng.standard_normal(shape) + 0j
    phi = ed._project_parity_dynamic(g, "singlet").ravel()
    phi = phi + 0.3 * a + 0.3 * b
    return (phi / np.linalg.norm(phi)).reshape(shape)


def _run(eli, info=None, shape=GAP_SHAPE):
    matvec, _a, _b = _repulsive_matvec(shape)
    return ed.run_leading_eigenproblem(
        matvec, shape, eli, "singlet", phi0=_seed(shape), seed_vec=None,
        use_ir=False, axF=None, nmat=shape[-1], info_out=info)


class TestIterationNonConvergenceIsNotAnEigenvalue(unittest.TestCase):

    def test_unshifted_two_cycle_reports_nan(self):
        """The two-cycle is recorded for the diagnosis, but the run goes on
        to max_iter: no finite window can rule out a tiny positive dominant
        component that is still growing."""
        info = {}
        out = _run({"solver_mode": "iteration", "spectral_shift": "none",
                    "max_iter": 300, "alpha": 0.5,
                    "convergence_tol": 1.0e-8}, info)
        lam = out[0]
        self.assertTrue(np.isnan(lam), lam)
        self.assertIs(info["eigen_converged"], False)
        self.assertIs(info["two_cycle"], True)
        self.assertEqual(info["n_iter"], 300)

    def test_iteration_cap_reports_nan(self):
        info = {}
        out = _run({"solver_mode": "iteration", "max_iter": 3,
                    "convergence_tol": 1.0e-10}, info)
        self.assertTrue(np.isnan(out[0]))
        self.assertIs(info["eigen_converged"], False)
        self.assertIs(info["two_cycle"], False)

    def test_driver_two_cycle_with_other_alpha(self):
        """The two-cycle is detected for other dampings than alpha = 1/2.
        (For alpha >= ~0.7 the damped map has no stable two-cycle on this
        spectrum and the iteration converges to LAM_POS instead -- which the
        Rayleigh report then gets right.)"""
        matvec, a, b = _repulsive_matvec()
        from scipy.sparse.linalg import LinearOperator
        n = int(np.prod(GAP_SHAPE))
        A = LinearOperator((n, n), matvec=matvec, dtype=complex)
        for alpha in (0.0, 0.2, 0.6):
            with self.subTest(alpha=alpha):
                val, _vec, info = sc._solve_leading(
                    lambda: (A, n), n, "iteration", max_iter=5000,
                    convergence_tol=1.0e-8, alpha=alpha,
                    init_vec=_seed().ravel())
                self.assertTrue(np.isnan(val))
                self.assertFalse(info["converged"])
                self.assertTrue(info["two_cycle"])
                self.assertEqual(info["n_iter"], 5000)


class TestDefaultShiftSelectsLargestRealPart(unittest.TestCase):

    def test_default_iteration_converges_to_the_positive_eigenvalue(self):
        info = {}
        out = _run({"solver_mode": "iteration", "max_iter": 5000,
                    "convergence_tol": 1.0e-11}, info)
        self.assertIs(info["eigen_converged"], True)
        self.assertAlmostEqual(out[0], LAM_POS, delta=1.0e-8)

    def test_eigenvalue_mode_reports_the_positive_eigenvalue(self):
        info = {}
        out = _run({"solver_mode": "eigenvalue", "num_eigenvalues": 2}, info)
        self.assertIs(info["eigen_converged"], True)
        self.assertAlmostEqual(out[0], LAM_POS, delta=1.0e-8)

    def test_arpack_lm_set_without_the_lr_eigenvalue_is_resolved(self):
        """``which='LM'`` with one requested eigenvalue returns only LAM_NEG;
        the shared driver must still report the largest real part."""
        matvec, _a, _b = _repulsive_matvec()
        from scipy.sparse.linalg import LinearOperator
        n = int(np.prod(GAP_SHAPE))
        A = LinearOperator((n, n), matvec=matvec, dtype=complex)
        val, _vec, info = sc._solve_leading(lambda: (A, n), n, "arnoldi",
                                            num_eigenvalues=1)
        self.assertAlmostEqual(val.real, LAM_POS, delta=1.0e-8)
        self.assertEqual(info["selection"], "LR_retry")

    def test_spin_orbital_dense_branch_selects_largest_real_part(self):
        """The spin-orbital dense branch (num_eigenvalues >= dim - 1) used to
        keep the largest-MAGNITUDE eigenpairs, which drops LAM_POS when it is
        the smallest in magnitude."""
        shape = (2, 2, 2, 1, 1)
        basis = sc._SOAntisymmetricBasis(shape)
        dim = basis.dim
        rng = np.random.default_rng(5)
        q, _ = np.linalg.qr(rng.standard_normal((dim, dim)))
        diag = np.full(dim, -0.6)
        diag[0], diag[1] = LAM_NEG, LAM_POS
        red = q @ np.diag(diag) @ q.T

        from scipy.sparse.linalg import LinearOperator
        n = int(np.prod(shape))

        def mv(x):
            return basis.to_full(red @ basis.to_reduced(np.asarray(x).ravel()))

        A = LinearOperator((n, n), matvec=mv, dtype=complex)
        vals, _gaps = sc._so_solve_eigenvalue(A, basis, dim - 1, "arnoldi")
        self.assertAlmostEqual(vals[0].real, LAM_POS, delta=1.0e-10)


class TestConvergedIterationReportsSignedRayleigh(unittest.TestCase):

    def test_matches_dense_eigenvalue_shifted_and_unshifted(self):
        rng = np.random.default_rng(11)
        n = 40
        # non-Hermitian, real spectrum: dominant +0.8, next 0.5
        s = rng.standard_normal((n, n)) + 0.3 * np.eye(n)
        d = np.concatenate([[0.8, 0.5], 0.1 * rng.uniform(-1, 1, n - 2)])
        M = s @ np.diag(d) @ np.linalg.inv(s)
        from scipy.sparse.linalg import LinearOperator
        A = LinearOperator((n, n), matvec=lambda v: M @ v, dtype=complex)
        ref = np.max(np.linalg.eigvals(M).real)
        for shift in (None, "auto"):
            with self.subTest(shift=shift):
                val, _vec, info = sc._solve_leading(
                    lambda: (A, n), n, "iteration", max_iter=20000,
                    convergence_tol=1.0e-12, alpha=0.5,
                    init_vec=np.ones(n, dtype=complex), spectral_shift=shift)
                self.assertTrue(info["converged"])
                self.assertAlmostEqual(float(val), ref, delta=1.0e-8)
                self.assertAlmostEqual(info["rayleigh_eigenvalue"], ref,
                                       delta=1.0e-8)


def _dense_op(M):
    from scipy.sparse.linalg import aslinearoperator
    M = np.asarray(M, dtype=complex)
    n = M.shape[0]
    return (lambda: (aslinearoperator(M), n)), n


def _similar(diag, seed=5, strength=0.5):
    """A non-normal matrix with the given eigenvalues and its eigenvectors."""
    rng = np.random.default_rng(seed)
    n = len(diag)
    S = np.eye(n) + strength * rng.standard_normal((n, n)) / np.sqrt(n)
    return S @ np.diag(diag) @ np.linalg.inv(S), S


class TestSeedOnAnEigenvector(unittest.TestCase):
    """A seed that is (nearly) an eigenvector of K, or that K annihilates,
    cannot search its symmetry sector: never a converged value of that mode
    unless the iteration demonstrably left it."""

    def test_seed_on_a_negative_subdominant_mode(self):
        rng = np.random.default_rng(5)
        diag = [1.0, -0.9] + list(rng.uniform(-0.05, 0.05, 38))
        M, S = _similar(diag)
        mk, n = _dense_op(M)
        ref = max(np.linalg.eigvals(M).real)
        for eps in (0.0, 1e-8, 1e-6):
            v0 = (S[:, 1] + eps * S[:, 0]).astype(complex)
            for shift in ("none", "auto"):
                with self.subTest(eps=eps, shift=shift):
                    lam, _v, info = sc._solve_leading(
                        mk, n, "iteration", max_iter=5000,
                        convergence_tol=1e-8, init_vec=v0 / np.linalg.norm(v0),
                        alpha=0.5, spectral_shift=shift)
                    # used to be -0.9 "converged" at step 1
                    if info["converged"]:
                        self.assertAlmostEqual(lam, ref, delta=1e-6)
                    else:
                        self.assertTrue(np.isnan(lam))
                    if eps == 0.0:
                        self.assertEqual(info["eigenvalue_kind"],
                                         "seed-eigenvector")
                        self.assertEqual(info["n_iter"], 0)
                        self.assertAlmostEqual(info["seed_eigenvalue"], -0.9,
                                               delta=1e-10)

    def test_annihilated_seed_is_reported_not_converged(self):
        """K v0 = 0: under a shift v0 is an eigenvector of K + sigma*I; the
        iteration used to stop at step 1 with lambda ~ 0."""
        M = np.diag([0.0, 0.7, -0.9, 0.2]).astype(complex)
        mk, n = _dense_op(M)
        v0 = np.array([1.0, 0, 0, 0], dtype=complex)
        from tests.approx_util import _caplog
        with _caplog("WARNING", logger="hwave_sc") as records:
            lam, _v, info = sc._solve_leading(
                mk, n, "iteration", max_iter=5000, convergence_tol=1e-10,
                init_vec=v0, alpha=0.5, spectral_shift="auto")
        self.assertTrue(np.isnan(lam))
        self.assertFalse(info["converged"])
        self.assertEqual(info["eigenvalue_kind"], "seed-annihilated")
        msg = " ".join(r.message for r in records)
        self.assertIn("initial gap is annihilated", msg)
        self.assertIn('solver_mode="eigenvalue"', msg)

    def test_kernel_that_annihilates_everything_is_not_converged(self):
        mk, n = _dense_op(np.zeros((6, 6)))
        lam, _v, info = sc._solve_leading(
            mk, n, "iteration", max_iter=200, convergence_tol=1e-8,
            init_vec=np.ones(6, dtype=complex), alpha=0.5,
            spectral_shift="auto")
        self.assertTrue(np.isnan(lam))
        self.assertFalse(info["converged"])

    def test_exact_seed_eigenvector_is_not_perturbed_out_of_its_sector(self):
        """Two invariant sectors A (seed's, leading 0.5) and B (leading 0.8):
        the iteration must never report B's 0.8 -- the seed's sector is the
        documented search space. An exact eigenvector seed gives NaN; a seed
        close to the 0.5 eigenvector (the correct leading one of its sector)
        converges to 0.5."""
        rng = np.random.default_rng(1)
        dA = np.r_[0.5, 0.3, rng.uniform(-0.1, 0.1, 8)]
        dB = np.r_[0.8, rng.uniform(-0.1, 0.1, 9)]
        mk, n = _dense_op(np.diag(np.r_[dA, dB]))
        for w in (0.0, 5e-4, 5e-3):
            v0 = np.zeros(n, complex)
            v0[0], v0[1] = 1.0, w
            for shift in ("none", "auto"):
                with self.subTest(w=w, shift=shift):
                    lam, v, info = sc._solve_leading(
                        mk, n, "iteration", max_iter=3000,
                        convergence_tol=1e-5, init_vec=v0, alpha=0.5,
                        spectral_shift=shift)
                    self.assertEqual(np.linalg.norm(v[10:]), 0.0)
                    if w == 0.0:
                        self.assertTrue(np.isnan(lam))
                    else:
                        self.assertAlmostEqual(lam, 0.5, delta=1e-8)

    def test_large_shift_does_not_accept_a_slow_drift(self):
        """K = diag(0.3, 0.27, -40, ...), auto shift sigma ~ 42, seed on the
        0.27 mode plus 1e-6 of the 0.3 mode: the iterate drifts by tiny steps
        and used to be accepted as converged on 0.27 after 528 steps (here:
        6). It must never report 0.27."""
        rng = np.random.default_rng(0)
        d = rng.uniform(-0.05, 0.05, 20)
        d[0], d[1], d[2] = 0.3, 0.27, -40.0
        mk, n = _dense_op(np.diag(d))
        for amp in (1e-6, 2e-2):
            v0 = np.eye(n, dtype=complex)[1] + amp * np.eye(n)[0]
            with self.subTest(amp=amp):
                lam, _v, info = sc._solve_leading(
                    mk, n, "iteration", max_iter=1000, convergence_tol=1e-5,
                    init_vec=v0, alpha=0.5, spectral_shift="auto")
                if info["converged"]:
                    self.assertAlmostEqual(lam, 0.3, delta=1e-4)
                else:
                    self.assertTrue(np.isnan(lam))


class TestTwoCycleDetectorNeedsASustainedCycle(unittest.TestCase):

    def test_transient_anti_alignment_does_not_stop_the_run(self):
        """diag(lneg, 2, 0.01), plain power iteration (alpha = 0): the seed
        sits on the negative mode with a tiny component of the dominant +2
        mode (and a decaying 0.01 mode). The iterates anti-align for many
        steps while the +2 component grows; the run used to stop there as a
        "two-cycle" (after 5 steps for 1e-9, after 11 for 1e-16 at lneg =
        -1.9) and must converge to 2 when max_iter allows the component to
        grow -- and honestly report NaN when it does not."""
        for lneg in (-1.0, -1.9):
            for amp in (1e-9, 1e-16, 1e-20):
                with self.subTest(lneg=lneg, amp=amp):
                    M = np.diag([lneg, 2.0, 0.01]).astype(complex)
                    mk, n = _dense_op(M)
                    v0 = np.array([1.0, amp, 0.3], dtype=complex)
                    lam, _v, info = sc._solve_leading(
                        mk, n, "iteration", max_iter=3000,
                        convergence_tol=1e-6, init_vec=v0, alpha=0.0,
                        spectral_shift="none")
                    self.assertTrue(info["converged"])
                    self.assertAlmostEqual(lam, 2.0, delta=1e-8)
        # 1e-20 at lneg = -1.9 needs ~900 steps to grow (factor 2/1.9 per
        # step): 300 are not enough, and the run says so
        M = np.diag([-1.9, 2.0, 0.01]).astype(complex)
        mk, n = _dense_op(M)
        lam, _v, info = sc._solve_leading(
            mk, n, "iteration", max_iter=300, convergence_tol=1e-6,
            init_vec=np.array([1.0, 1e-20, 0.3], dtype=complex), alpha=0.0,
            spectral_shift="none")
        self.assertTrue(np.isnan(lam))
        self.assertFalse(info["converged"])

    def test_dominant_complex_pair_is_never_a_real_value(self):
        c, s_ = np.cos(1.0), np.sin(1.0)
        M = np.zeros((5, 5))
        M[:2, :2] = 2.0 * np.array([[c, -s_], [s_, c]])     # 2 exp(+-i)
        M[2, 2], M[3, 3], M[4, 4] = 0.3, -0.2, 0.1
        mk, n = _dense_op(M)
        from tests.approx_util import _caplog
        with _caplog("WARNING", logger="hwave_sc") as records:
            lam, _v, info = sc._solve_leading(
                mk, n, "iteration", max_iter=300, convergence_tol=1e-8,
                init_vec=np.ones(5, dtype=complex), alpha=0.5,
                spectral_shift="none")
        self.assertTrue(np.isnan(lam))
        self.assertFalse(info["converged"])
        self.assertFalse(info["two_cycle"])
        self.assertIn("complex", " ".join(r.message for r in records))

    def test_leading_eigenvalue_itself_complex_gives_nan_under_auto(self):
        M = np.diag([0.5 + 0.2j, -1.0, 0.05, -0.02])
        mk, n = _dense_op(M)
        lam, _v, info = sc._solve_leading(
            mk, n, "iteration", max_iter=400, convergence_tol=1e-8,
            init_vec=np.ones(4, dtype=complex), alpha=0.5,
            spectral_shift="auto")
        self.assertTrue(np.isnan(lam))
        self.assertFalse(info["converged"])


class TestMaxIterAdvice(unittest.TestCase):

    def test_small_gap_reports_the_steps_needed(self):
        """Gap 0.2 under a dominant -20: the auto shift (sigma ~ 21)
        compresses the gap and 1000 steps are not enough. The warning states
        the estimated count; that count does converge."""
        import re
        M = np.diag([0.5, 0.3, -20.0, 0.01, -0.03]).astype(complex)
        mk, n = _dense_op(M)
        kw = dict(convergence_tol=1e-8, init_vec=np.ones(5, dtype=complex),
                  alpha=0.5, spectral_shift="auto")
        from tests.approx_util import _caplog
        with _caplog("WARNING", logger="hwave_sc") as records:
            lam, _v, info = sc._solve_leading(mk, n, "iteration",
                                              max_iter=1000, **kw)
        self.assertTrue(np.isnan(lam))
        msg = " ".join(r.message for r in records)
        m = re.search(r"max_iter of about (\d+)", msg)
        self.assertIsNotNone(m, msg)
        need = int(m.group(1))
        lam, _v, info = sc._solve_leading(mk, n, "iteration",
                                          max_iter=int(1.2 * need), **kw)
        self.assertTrue(info["converged"])
        self.assertAlmostEqual(lam, 0.5, delta=1e-6)


class TestResidualCheckScale(unittest.TestCase):
    """The residual bound used to grow with the lambda being checked, so a
    strongly non-normal kernel's off-spectrum Rayleigh value passed at tol
    1e-3 (lambda = 12.34 with spectral radius 0.75, 42.64 with 7.12, 4.20
    with 0.776). Kernels and seeds saved from that search."""

    CASES = {
        # name: (spectral radius, largest real part)
        "1_109": (0.75, 0.3958),
        "1_231": (7.12, 0.6579),
        "2_61": (0.776, 0.7758),
    }

    def test_no_value_beyond_the_spectral_radius(self):
        data = np.load(os.path.join(os.path.dirname(__file__), "sc",
                                    "offspectrum_kernels.npz"))
        for name, (rho, lr) in self.CASES.items():
            M, v0, alpha = (data[name + "_M"], data[name + "_v0"],
                            float(data[name + "_alpha"]))
            mk, n = _dense_op(M)
            for tol in (1e-3, 1e-5):
                with self.subTest(case=name, tol=tol):
                    lam, _v, info = sc._solve_leading(
                        mk, n, "iteration", max_iter=3000,
                        convergence_tol=tol, init_vec=v0.copy(),
                        alpha=alpha, spectral_shift="auto")
                    if info["converged"]:
                        self.assertLessEqual(abs(lam), 1.01 * rho)
                        self.assertAlmostEqual(lam, lr, delta=0.02)
                    else:
                        self.assertTrue(np.isnan(lam))


class TestLargestRealSolveFailure(unittest.TestCase):
    """When ARPACK's largest-real-part solve does not converge, the
    largest-magnitude set must not be reported as the leading eigenvalue."""

    def _run(self, dense_max):
        from unittest import mock
        import scipy.sparse.linalg as spl
        real_eigs = sc.eigs
        M, _S = _similar([-3.0, 2.5, 0.2] + [-1.0 - 0.01 * k for k in range(37)])
        mk, n = _dense_op(M)

        def eigs(A, k=6, which="LM", **kw):
            if which == "LR":
                raise spl.ArpackNoConvergence("forced", np.array([]),
                                              np.array([]))
            return real_eigs(A, k=k, which=which, **kw)
        with mock.patch.object(sc, "eigs", eigs), \
                mock.patch.object(sc, "_DENSE_FALLBACK_MAX", dense_max):
            # one LM eigenvalue: {-3}, which cannot be the largest real part
            return sc._solve_leading(mk, n, "arnoldi", num_eigenvalues=1), M

    def test_small_operator_falls_back_to_dense(self):
        (val, _vec, info), M = self._run(2000)
        self.assertEqual(info["selection"], "dense-LR")
        self.assertAlmostEqual(val.real, max(np.linalg.eigvals(M).real),
                               delta=1e-9)

    def test_failed_auto_shift_probe_is_named_as_such(self):
        """With spectral_shift = "auto" the LM radius probe runs first; when
        it fails, the warning must name the probe, not the LR solve."""
        from unittest import mock
        import scipy.sparse.linalg as spl
        from tests.approx_util import _caplog
        M, _S = _similar([-3.0, 2.5, 0.2] + [-1.0 - 0.01 * k for k in range(37)])
        mk, n = _dense_op(M)
        real_eigs = sc.eigs

        def eigs(A, k=6, which="LM", **kw):
            if which == "LM":
                raise spl.ArpackNoConvergence("forced", np.array([]),
                                              np.array([]))
            return real_eigs(A, k=k, which=which, **kw)
        with mock.patch.object(sc, "eigs", eigs), \
                _caplog("WARNING", logger="hwave_sc") as records:
            val, _vec, info = sc._solve_leading(mk, n, "arnoldi",
                                                num_eigenvalues=2,
                                                spectral_shift="auto")
        msg = " ".join(r.message for r in records)
        self.assertIn("spectral-radius probe", msg)
        self.assertNotIn("largest-real-part solve (which='LR') did not", msg)
        self.assertEqual(info["selection"], "dense-LR")
        self.assertAlmostEqual(val.real, max(np.linalg.eigvals(M).real),
                               delta=1e-9)

    def test_large_operator_reports_nan(self):
        (val, _vec, info), _M = self._run(10)
        self.assertTrue(np.isnan(val))
        self.assertTrue(info["lr_retry_failed"])


class TestSpectralShiftOption(unittest.TestCase):

    def test_resolution(self):
        r = sc._resolve_spectral_shift
        self.assertEqual(r(None, "iteration"), "auto")
        self.assertIsNone(r("none", "iteration"))
        self.assertIsNone(r("None", "iteration"))
        self.assertEqual(r(2.5, "iteration"), 2.5)
        self.assertIsNone(r(None, "eigenvalue"))
        self.assertIsNone(r("none", "eigenvalue"))
        self.assertEqual(r("auto", "eigenvalue"), "auto")


class TestWritersCarryTheConvergenceFlag(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _write(self, lam, flag):
        gap = np.zeros(GAP_SHAPE, dtype=complex)
        gap[0, 0, 0, 0, 0, :] = 1.0
        path = os.path.join(self.tmp, "eigenvalue.dat")
        ed.write_eigenvalue_file(path, lam, None, None, None,
                                 eigen_converged=flag)
        ed.write_dynamic_outputs(
            self.tmp, gap, lam, 0.1, "singlet", np.zeros(4), np.zeros(4),
            np.zeros(1), 10.0, eigen_converged=flag)
        with open(path) as fh:
            text = fh.read()
        z = np.load(os.path.join(self.tmp, "gap_dynamic.npz"))
        return text, z

    def test_not_converged(self):
        text, z = self._write(float("nan"), False)
        self.assertIn("# eigen_converged=false", text)
        values = [l for l in text.splitlines() if not l.startswith("#")]
        self.assertTrue(np.isnan(float(values[0])))
        self.assertIs(bool(z["eigen_converged"]), False)
        self.assertTrue(np.isnan(float(z["eigenvalue"])))
        with open(os.path.join(self.tmp, "gap.dat")) as fh:
            self.assertIn("eigen_converged=false", fh.readline())

    def test_converged(self):
        text, z = self._write(0.5, True)
        self.assertIn("# eigen_converged=true", text)
        self.assertIs(bool(z["eigen_converged"]), True)


class TestCliExitStatus(unittest.TestCase):
    """``hwave_sc`` exits non-zero when the eigen-solver did not converge
    (after writing its outputs)."""

    def _main(self, result):
        from unittest import mock
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        path = os.path.join(tmp, "in.toml")
        with open(path, "w") as fh:
            fh.write("[mode]\nmode = \"SC\"\n")
        with mock.patch.object(sc, "calc_eliashberg", return_value=result), \
                mock.patch("sys.argv", ["hwave_sc", "-q", path]):
            try:
                sc.main()
            except SystemExit as exc:
                return exc.code
        return 0

    def test_nan_exits_nonzero(self):
        self.assertEqual(self._main(float("nan")), sc.EXIT_EIGEN_NOT_CONVERGED)
        self.assertNotEqual(sc.EXIT_EIGEN_NOT_CONVERGED, 0)

    def test_converged_exits_zero(self):
        self.assertEqual(self._main(0.4), 0)
        self.assertEqual(self._main(None), 0)


class TestHwaveScDefaultSolverMode(unittest.TestCase):

    def test_default_is_the_eigenvalue_solver(self):
        self.assertEqual(sc.DEFAULT_SOLVER_MODE, "eigenvalue")
        from hwave.solver.eliashberg_bond_io import PairingControls
        ctl = PairingControls.from_param({}, pairing_types=("singlet",))
        self.assertEqual(ctl.solver_mode, "eigenvalue")

    def test_dynamic_run_without_solver_mode_uses_arnoldi(self):
        matvec, _a, _b = _repulsive_matvec()
        info = {}
        out = ed.run_leading_eigenproblem(
            matvec, GAP_SHAPE, {"num_eigenvalues": 2}, "singlet",
            phi0=_seed(), seed_vec=None, use_ir=False, axF=None,
            nmat=GAP_SHAPE[-1], info_out=info)
        self.assertNotEqual(out[8], "iteration")
        self.assertAlmostEqual(out[0], LAM_POS, delta=1e-8)


class TestInProcessDefaultSolverMode(unittest.TestCase):

    def test_in_process_default_is_eigenvalue(self):
        from hwave.solver.flex import FLEX
        dummy = types.SimpleNamespace(longitudinal_bond_pairing=("singlet",))
        FLEX._install_pairing_controls(dummy, {})
        self.assertEqual(dummy._pairing_controls.solver_mode, "eigenvalue")
        FLEX._install_pairing_controls(dummy, {"solver_mode": "iteration"})
        self.assertEqual(dummy._pairing_controls.solver_mode, "iteration")



if __name__ == "__main__":
    unittest.main()
