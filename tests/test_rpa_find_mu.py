"""Direct unit tests for RPA._find_mu (spec
2026-08-28-mu-green-seam-160, Change 1) and the H0_k storage
(Change 1b). The mu tests drive _find_mu through a minimal duck-typed
host (only H0_eigenvalue and ene_cutoff are read), so no input files
are needed for the root-finder cases.
"""

import unittest
from contextlib import ExitStack

import numpy as np
from scipy import optimize

from hwave.solver.rpa import RPA
from hwave.solver.rpa import _masked_fermi_delta_n, _polish_mu_root


def _masked_fermi_count(w, T, mu, ene_cutoff):
    """Reference reimplementation of _find_mu's counter for asserting
    residuals (same masked arithmetic, independent code)."""
    x = (np.asarray(w).ravel() - mu) / T
    mask = x < ene_cutoff
    x1 = np.where(mask, x, 0.0)
    f = np.where(mask, 1.0 / (1.0 + np.exp(x1)), 0.0)
    return float(np.sum(f))


class _Host:
    ene_cutoff = 1.0e2

    def __init__(self, levels):
        self.H0_eigenvalue = np.asarray(levels, dtype=np.float64)


class TestH0kStorage(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from tests.equivalence_cells import GEEV_DIAGNOSTIC_FIXTURE
        from tests.test_rpa_flex_equivalence_table import build_solver
        cls._stack = ExitStack()
        cls.rpa, cls.green, _ = build_solver(
            GEEV_DIAGNOSTIC_FIXTURE, "rpa", cls._stack)
        cls.rpa._calc_epsilon_k(cls.green)

    @classmethod
    def tearDownClass(cls):
        cls._stack.close()

    def test_h0k_shape_dtype_and_spectral_consistency(self):
        H0 = self.rpa.H0_k
        ew = self.rpa.H0_eigenvalue
        ev = self.rpa.H0_eigenvector
        self.assertEqual(H0.dtype, np.complex128)
        self.assertEqual(H0.shape, ev.shape)
        rebuilt = np.matmul(ev * ew[:, :, np.newaxis, :],
                            np.conj(ev).swapaxes(-2, -1))
        self.assertLess(float(np.max(np.abs(rebuilt - H0))), 1e-12)

    def test_h0k_is_an_owning_non_writeable_copy(self):
        H0 = self.rpa.H0_k
        self.assertIsNone(H0.base)          # owning copy, no aliasing
        self.assertFalse(H0.flags.writeable)
        with self.assertRaises(ValueError):
            H0[0, 0, 0, 0] = 0.0


class TestMaskedFermiDeltaN(unittest.TestCase):

    def test_derivative_matches_finite_difference(self):
        rng = np.random.default_rng(7)
        w = rng.uniform(-3.0, 3.0, size=(1, 8, 3))
        T, target, cutoff, mu, h = 0.3, 10.0, 1.0e2, 0.4, 1.0e-6
        _n, dn = _masked_fermi_delta_n(w, T, mu, target, cutoff)
        n_plus, _ = _masked_fermi_delta_n(w, T, mu + h, target, cutoff)
        n_minus, _ = _masked_fermi_delta_n(w, T, mu - h, target, cutoff)
        fd = (n_plus - n_minus) / (2.0 * h)
        self.assertAlmostEqual(dn / fd, 1.0, delta=1e-5)


class TestPolishMuRoot(unittest.TestCase):

    def test_non_improving_candidate_is_rejected(self):
        # Constant residual: candidate cannot improve -> keep the root.
        calls = []
        def f(mu):
            calls.append(mu)
            return 1.0, 1.0
        self.assertEqual(_polish_mu_root(f, 0.0), 0.0)

    def test_non_finite_candidate_residual_is_rejected(self):
        def f(mu):
            if mu == 0.0:
                return 1.0, 1.0
            return float("nan"), 1.0
        self.assertEqual(_polish_mu_root(f, 0.0), 0.0)

    def test_zero_or_nonfinite_derivative_stops_the_polish(self):
        self.assertEqual(_polish_mu_root(lambda m: (1.0, 0.0), 0.5), 0.5)
        self.assertEqual(
            _polish_mu_root(lambda m: (1.0, float("inf")), 0.5), 0.5)

    def test_candidate_clipped_to_endpoint_still_improves(self):
        # delta_n = mu - 5 (true root 5, outside bracket (0, 2)):
        # cand 5 -> clipped to 2, residual 3 < 5 -> accepted; second
        # round clips to 2 again, no further improvement -> returns 2.
        f = lambda mu: (mu - 5.0, 1.0)
        self.assertEqual(_polish_mu_root(f, 0.0, bracket=(0.0, 2.0)), 2.0)

    def test_linear_function_converges_in_one_step(self):
        f = lambda mu: (2.0 * (mu - 1.25), 2.0)
        self.assertEqual(_polish_mu_root(f, 0.0), 1.25)


class TestFindMuResidual(unittest.TestCase):
    """_find_mu invoked through a duck-typed host: the acceptance
    criterion is the ABSOLUTE particle-number residual (spec:
    Correctness criterion), never a specific mu value."""

    def _residual(self, host, ncond, T, mu):
        return abs(_masked_fermi_count(
            host.H0_eigenvalue, T, mu, host.ene_cutoff) - ncond)

    def test_ordinary_monotone_root_reaches_roundoff(self):
        host = _Host([[[-2.0, -1.0, 1.0, 2.0]]])
        dist, mu = RPA._find_mu(host, 2.0, 0.5)
        self.assertLess(self._residual(host, 2.0, 0.5, mu), 1e-12)

    def test_degenerate_eigenvalues(self):
        host = _Host([[[-1.0, -1.0, -1.0, 1.0]]])
        _dist, mu = RPA._find_mu(host, 1.5, 0.2)
        self.assertLess(self._residual(host, 1.5, 0.2, mu), 1e-12)

    def test_large_energy_scale_where_rtol_dominates_xtol(self):
        host = _Host([[[-2.0e6, -1.0e6, 1.0e6, 2.0e6]]])
        _dist, mu = RPA._find_mu(host, 2.0, 0.5e6)
        self.assertLess(self._residual(host, 2.0, 0.5e6, mu), 1e-10)

    def test_no_sign_change_fallback_never_worsens(self):
        # target 0.1 is reachable only below ev[0]: delta_n has the
        # same sign at both ends of the eigenvalue span. Since #218 the
        # root is found by the widened bracket (newton is reached only
        # for targets with no root); on this fixture the legacy newton
        # happened to converge too, so the two agree to round-off.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        _dist, mu = RPA._find_mu(host, 0.1, 0.5)
        self.assertLess(self._residual(host, 0.1, 0.5, mu), 1e-10)
        w = host.H0_eigenvalue
        f = lambda m: _masked_fermi_delta_n(w, 0.5, m, 0.1, host.ene_cutoff)
        mu_newton, _r = optimize.newton(lambda m: f(m)[0], np.min(w),
                                        full_output=True)
        mu_legacy = _polish_mu_root(f, mu_newton, bracket=None)
        self.assertAlmostEqual(mu, mu_legacy, delta=1e-10)

    def test_charge_gap_plateau_accepts_any_residual_minimal_point(self):
        # T=0.01, gap [-1,-1 | 1,1], target exactly 2: any mu deep in
        # the gap has residual ~0; assert the residual, NOT the mu.
        host = _Host([[[-1.0, -1.0, 1.0, 1.0]]])
        _dist, mu = RPA._find_mu(host, 2.0, 0.01)
        self.assertLess(self._residual(host, 2.0, 0.01, mu), 1e-12)


def _cosine_band_8x8():
    """One-orbital -2(cos kx + cos ky) band on 8x8 (64 levels): the
    reproduction of issue #218."""
    k = 2.0 * np.pi * np.arange(8) / 8
    kx, ky = np.meshgrid(k, k, indexing="ij")
    return (-2.0 * (np.cos(kx) + np.cos(ky)))[:, :, np.newaxis]


class TestFindMuWidenedBracket(unittest.TestCase):
    """Issue #218: the root lies OUTSIDE the eigenvalue span (near-full /
    near-empty filling at finite T, large T, small ene_cutoff). The
    eigenvalue span is then not sign-bracketed and the unbracketed
    newton from ev[0] stalls; _find_mu now retries brentq on a widened,
    evaluated-sign bracket first. Residuals are measured with the file's
    independent counter."""

    def _residual(self, host, ncond, T, mu):
        return abs(_masked_fermi_count(
            host.H0_eigenvalue, T, mu, host.ene_cutoff) - ncond)

    # --- new behaviour --------------------------------------------------
    # Before #218 the near-full case and the mask-jump case raised
    # RuntimeError from the unbracketed newton; the near-empty, large-T
    # and small-cutoff cases happened to converge through newton because
    # their root lies below the ev[0] start point. All five now take the
    # widened bracket.

    def test_near_full_filling_root_above_band_top(self):
        host = _Host(_cosine_band_8x8())
        ncond = 0.999 * host.H0_eigenvalue.size
        _dist, mu = RPA._find_mu(host, ncond, 0.1)
        self.assertLess(self._residual(host, ncond, 0.1, mu), 1e-12)
        self.assertGreater(mu, np.max(host.H0_eigenvalue))

    def test_near_empty_filling_root_below_band_bottom(self):
        host = _Host(_cosine_band_8x8())
        ncond = 0.001 * host.H0_eigenvalue.size
        _dist, mu = RPA._find_mu(host, ncond, 0.1)
        self.assertLess(self._residual(host, ncond, 0.1, mu), 1e-12)
        self.assertLess(mu, np.min(host.H0_eigenvalue))

    def test_large_temperature_root_far_outside_a_fixed_margin(self):
        # T = 100: the root sits ~600 above the band, beyond any fixed
        # +-10 margin; the widened bracket scales with T.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        _dist, mu = RPA._find_mu(host, 3.99, 100.0)
        self.assertLess(self._residual(host, 3.99, 100.0, mu), 1e-10)
        self.assertGreater(mu, 2.0 + 100.0)

    def test_small_ene_cutoff_still_brackets_a_near_full_target(self):
        # With a margin of only T*ene_cutoff = 5 the upper endpoint gives
        # sum f = 3.999948 (the top level alone is 5e-5 short of 1), i.e.
        # delta_n < 0 for this target whose root sits at ~7.76 > 7; the
        # margin floor of 40 saturates every level exactly.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        host.ene_cutoff = 10.0
        w = host.H0_eigenvalue
        ncond = 3.99999
        f = lambda m: _masked_fermi_delta_n(w, 0.5, m, ncond, 10.0)[0]
        # not in-band (no strict sign change over the eigenvalue span)
        self.assertGreaterEqual(f(np.min(w)) * f(np.max(w)), 0.0)
        # ... and T*ene_cutoff alone would NOT bracket the root either
        self.assertLess(f(np.max(w) + 0.5 * 10.0), 0.0)
        _dist, mu = RPA._find_mu(host, ncond, 0.5)
        self.assertLess(self._residual(host, ncond, 0.5, mu), 1e-10)
        self.assertGreater(mu, np.max(w) + 0.5 * 10.0)

    def test_absorbed_seed_margin_is_replaced_by_the_adjacent_float(self):
        # ulp(1e20) = 16384: the seed T*40 = 4000 rounds back onto the
        # eigenvalue; the adjacent float 1e20 - 16384 (x = 163.8 >>
        # ene_cutoff = 10) is used instead and brackets the root at once
        # (a version that only doubled the margin raised RuntimeError
        # here). The exact root (1082 below 1e20) is NOT a representable
        # float, so the pin is the best REPRESENTABLE residual: no worse
        # than either float neighbour of the eigenvalue.
        host = _Host([[[1.0e20]]])
        host.ene_cutoff = 10.0
        ncond, T = 2.0e-5, 100.0
        self.assertEqual(1.0e20 - T * 40.0, 1.0e20)      # absorbed
        _dist, mu = RPA._find_mu(host, ncond, T)
        self.assertTrue(np.isfinite(mu))
        self.assertLess(mu, 1.0e20)                        # went outward
        best = min(self._residual(host, ncond, T, 1.0e20),
                   self._residual(host, ncond, T, np.nextafter(1.0e20, 0.0)))
        self.assertLessEqual(self._residual(host, ncond, T, mu), best)

    def test_target_inside_a_mask_jump_returns_a_bounded_residual(self):
        # ene_cutoff = 10: each level's mask boundary is a jump of
        # 1/(1+exp(10)) = 4.54e-5 in delta_n, so a target below that has
        # NO exact root (pre-existing property of the mask, shared by
        # the in-band branch). The widened search still returns a finite
        # mu whose residual is at most the jump, and the transactional
        # polish never worsens the brentq point.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        host.ene_cutoff = 10.0
        ncond, T = 2.0e-5, 0.5
        jump = 1.0 / (1.0 + np.exp(10.0))
        _dist, mu = RPA._find_mu(host, ncond, T)
        self.assertTrue(np.isfinite(mu))
        self.assertLessEqual(self._residual(host, ncond, T, mu), jump * 1.01)
        w = host.H0_eigenvalue
        f = lambda m: _masked_fermi_delta_n(w, T, m, ncond, host.ene_cutoff)[0]
        mu_b = optimize.brentq(f, np.min(w) - T * 40.0, np.max(w) + T * 40.0,
                               xtol=1e-14)
        self.assertLessEqual(self._residual(host, ncond, T, mu),
                             self._residual(host, ncond, T, mu_b) + 1e-15)

    # --- compatibility pins ------------------------------------------

    def test_target_outside_zero_to_count_still_raises(self):
        # No root anywhere: the legacy unbracketed newton is reached and
        # fails exactly as before.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        with self.assertRaises(RuntimeError):
            RPA._find_mu(host, 5.0, 0.5)
        with self.assertRaises(RuntimeError):
            RPA._find_mu(host, -1.0, 0.5)

    def test_absorbed_margin_uses_the_adjacent_float(self):
        # T*S is far below the eigenvalue ulp (1e20 - 0.04 == 1e20, and
        # no bounded number of doublings would change that): the
        # endpoints are replaced by the adjacent floats, which are >=
        # 2*T*S away and bracket the root. The exact root is not a
        # representable float here either, so the pin is the best
        # representable residual around the top level. (Raised
        # RuntimeError from the legacy path before.)
        host = _Host([[[1.0e20, 2.0e20]]])
        ncond, T = 1.999, 1.0e-3
        _dist, mu = RPA._find_mu(host, ncond, T)
        self.assertTrue(np.isfinite(mu))
        best = min(self._residual(host, ncond, T, 2.0e20),
                   self._residual(host, ncond, T, np.nextafter(2.0e20, np.inf)))
        self.assertLessEqual(self._residual(host, ncond, T, mu), best)

    def test_absorbed_margin_needing_many_doublings_still_brackets(self):
        # ene_cutoff = 100, T = 0.5: the seed displacement T*S = 50 and
        # its first eight doublings (up to 6400) are all absorbed by
        # ulp(1e20) = 16384; the adjacent float makes the very first
        # attempt bracket instead. Best-representable pin as above.
        host = _Host([[[1.0e20]]])
        ncond, T = 0.999, 0.5
        self.assertEqual(1.0e20 + T * 100.0 * 2.0 ** 7, 1.0e20)   # absorbed
        _dist, mu = RPA._find_mu(host, ncond, T)
        self.assertTrue(np.isfinite(mu))
        self.assertGreater(mu, 1.0e20)                            # went outward
        best = min(self._residual(host, ncond, T, 1.0e20),
                   self._residual(host, ncond, T, np.nextafter(1.0e20, np.inf)))
        self.assertLessEqual(self._residual(host, ncond, T, mu), best)

    def test_in_band_root_is_bit_identical_to_the_span_bracketed_path(self):
        # Local reconstruction of the pre-#218 in-band path (brentq on
        # the eigenvalue span, polish clipped to it). Each fixture first
        # proves it IS in-band (strict sign change on the span), so the
        # routing -- not the polisher -- is what is pinned.
        # (The degenerate fixture of test_degenerate_eigenvalues, target
        # 1.5, is NOT in-band: its root sits just below ev[0]; target 2.5
        # on the same levels is.)
        fixtures = [
            ([[[-2.0, -1.0, 1.0, 2.0]]], 2.0, 0.5),
            ([[[-1.0, -1.0, -1.0, 1.0]]], 2.5, 0.2),
            ([[[-2.0e6, -1.0e6, 1.0e6, 2.0e6]]], 2.0, 0.5e6),
        ]
        for levels, ncond, T in fixtures:
            with self.subTest(levels=levels, ncond=ncond, T=T):
                host = _Host(levels)
                w = host.H0_eigenvalue
                ev = np.sort(w.flatten())
                f = lambda m: _masked_fermi_delta_n(w, T, m, ncond,
                                                    host.ene_cutoff)
                self.assertLess(f(ev[0])[0] * f(ev[-1])[0], 0.0)
                mu_b, _r = optimize.brentq(lambda m: f(m)[0], ev[0], ev[-1],
                                           xtol=1e-14, full_output=True,
                                           disp=False)
                mu_ref = _polish_mu_root(f, mu_b, bracket=(ev[0], ev[-1]))
                _dist, mu = RPA._find_mu(host, ncond, T)
                self.assertEqual(mu, mu_ref)
