"""Exact shapes of the vectorised Lehmann frequency sum
(``ed_oracle_util.lehmann_frequency_sum``) through ``TotalNED.green``."""

import unittest

import numpy as np

from tests import ed_oracle_util as U
from tests.test_ed_cache import _fx, _terms


class TestGreenEmptyShapes(unittest.TestCase):
    """``TotalNED.green`` / ``lehmann_frequency_sum`` return exact zero
    arrays of the requested shape when a dimension is empty, as the
    per-frequency loop they replaced did."""

    def setUp(self):
        self.fx = _fx()
        self.ed = U.TotalNED(self.fx, terms=_terms(self.fx, 0.2))
        self.iws = 1j * (2 * np.arange(-3, 3) + 1) * np.pi * self.fx.T

    def _check(self, out, shape):
        self.assertEqual(out.shape, shape)
        self.assertEqual(out.dtype, np.complex128)
        self.assertTrue(np.array_equal(out, np.zeros(shape, complex)))

    def test_empty_rows_cols_both_and_frequencies(self):
        n, nw = self.fx.nmode, len(self.iws)
        self._check(self.ed.green(self.iws, rows=[]), (nw, 0, n))
        self._check(self.ed.green(self.iws, cols=[]), (nw, n, 0))
        self._check(self.ed.green(self.iws, rows=[], cols=[]), (nw, 0, 0))
        self._check(self.ed.green(np.zeros(0, complex)), (0, n, n))
        self._check(self.ed.green(np.zeros(0, complex), rows=[1], cols=[]), (0, 1, 0))

    def test_frequency_sum_with_an_empty_operand(self):
        A = np.ones((2, 5), complex)
        w, dE = np.ones(5), np.arange(5.0)
        self._check(U.lehmann_frequency_sum(self.iws, A[:0], A, w, dE), (6, 0, 2))
        self._check(U.lehmann_frequency_sum(self.iws, A, A[:0], w, dE), (6, 2, 0))
        self._check(U.lehmann_frequency_sum(self.iws[:0], A, A, w, dE), (0, 2, 2))

    def test_single_state_sector(self):
        """L = 1, one orbital: the N = 0 and N = 2 sectors hold one state
        each; G equals the closed form of the Hubbard atom."""
        fx = U.EDFixture(L=1, norb=1, t={(0, 0): 0.0}, eps=(0.3,), T=0.5, mu=0.1)
        Uv = 0.7
        terms = U.canonical_density_terms(fx, [(0, 0, 0, Uv)])
        g = U.TotalNED(fx, terms=terms).green(self.iws, rows=[0], cols=[0, 1])
        self.assertEqual(g.shape, (len(self.iws), 1, 2))
        self.assertTrue(np.array_equal(g[:, 0, 1], np.zeros(len(self.iws), complex)))
        e = 0.3 - 0.1
        b = fx.beta
        z = 1 + 2 * np.exp(-b * e) + np.exp(-b * (2 * e + Uv))
        n1 = (np.exp(-b * e) + np.exp(-b * (2 * e + Uv))) / z    # <n_up>
        want = (1 - n1) / (self.iws - e) + n1 / (self.iws - e - Uv)
        np.testing.assert_allclose(g[:, 0, 0], want, rtol=1e-13, atol=0)


if __name__ == "__main__":
    unittest.main()
