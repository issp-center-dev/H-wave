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
