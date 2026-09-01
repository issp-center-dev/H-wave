"""Direct unit tests for RPA._find_mu (spec
2026-08-28-mu-green-seam-160, Change 1) and the H0_k storage
(Change 1b). The mu tests drive _find_mu through a minimal duck-typed
host (only H0_eigenvalue and ene_cutoff are read), so no input files
are needed for the root-finder cases.
"""

import unittest
from contextlib import ExitStack

import numpy as np

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
        # same sign at both bracket ends -> newton fallback + polish.
        host = _Host([[[-1.0, 0.0, 1.0, 2.0]]])
        _dist, mu = RPA._find_mu(host, 0.1, 0.5)
        self.assertLess(self._residual(host, 0.1, 0.5, mu), 1e-10)

    def test_charge_gap_plateau_accepts_any_residual_minimal_point(self):
        # T=0.01, gap [-1,-1 | 1,1], target exactly 2: any mu deep in
        # the gap has residual ~0; assert the residual, NOT the mu.
        host = _Host([[[-1.0, -1.0, 1.0, 1.0]]])
        _dist, mu = RPA._find_mu(host, 2.0, 0.01)
        self.assertLess(self._residual(host, 2.0, 0.01, mu), 1e-12)
