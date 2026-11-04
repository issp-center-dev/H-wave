"""The run-to-run comparison helpers of :mod:`tests.so_run_compare`."""

import unittest

import numpy as np

from tests.so_run_compare import (degenerate_sets, eigenvalues_close, numeric_mismatches,
                                  same_eigenvectors, same_leading_vector)


def _unit(rng, n):
    v = rng.standard_normal(n) + 1j * rng.standard_normal(n)
    return v / np.linalg.norm(v)


class TestRunCompare(unittest.TestCase):

    def setUp(self):
        self.rng = np.random.default_rng(7)

    def test_eigenvalues(self):
        self.assertTrue(eigenvalues_close([0.3, 2.0], [0.3 + 9e-13, 2.0 + 1.9e-12]))
        self.assertFalse(eigenvalues_close([0.3], [0.3 + 2e-12]))
        self.assertFalse(eigenvalues_close([0.3], [0.3, 0.1]))

    def test_degenerate_sets(self):
        self.assertEqual(degenerate_sets([0.5, 0.5 + 1e-12, 0.2]), [[0, 1], [2]])

    def test_phase_and_degenerate_basis_ignored(self):
        v = np.array([_unit(self.rng, 6) for _ in range(3)])
        vals = np.array([0.5, 0.5, 0.2])
        rot = np.array([[np.cos(0.4), -np.sin(0.4)], [np.sin(0.4), np.cos(0.4)]])
        w = v.copy()
        w[:2] = rot @ v[:2]
        w[2] = np.exp(1j * 0.7) * v[2]
        self.assertTrue(same_eigenvectors(vals, v, w))
        # a non-degenerate pair may not mix
        self.assertFalse(same_eigenvectors(np.array([0.5, 0.4, 0.2]), v, w))

    def test_leading_vector(self):
        v = np.array([_unit(self.rng, 6) for _ in range(2)])
        self.assertTrue(same_leading_vector([0.5, 0.2], v[0], 1j * v[0]))
        self.assertFalse(same_leading_vector([0.5, 0.2], v[0], v[1]))
        self.assertFalse(same_leading_vector([0.5, 0.2], v[0], 2 * v[0]))
        # degenerate leading set: any vector of the other run's span
        mix = (v[0] + v[1]) / np.linalg.norm(v[0] + v[1])
        self.assertTrue(same_leading_vector([0.5, 0.5], v[0], mix, v, v))
        self.assertFalse(same_leading_vector([0.5, 0.5], v[0], _unit(self.rng, 6), v, v))

    def test_numeric_mismatches(self):
        v = np.array([_unit(self.rng, 4) for _ in range(2)])
        a = {("gap.npz", "eigenvalues"): np.array([0.5, 0.2]),
             ("gap.npz", "eigenvectors"): v, ("gap.npz", "gap"): v[0],
             ("chiq.npz", "chiq"): np.full(3, 2.0), ("chiq.npz", "n"): np.int64(3),
             ("return", "lambda"): np.float64(0.5)}
        b = dict(a)
        b[("gap.npz", "eigenvectors")] = -v
        b[("gap.npz", "gap")] = 1j * v[0]
        b[("chiq.npz", "chiq")] = np.full(3, 2.0 + 1e-13)
        b[("return", "lambda")] = np.float64(0.5 + 1e-13)
        self.assertEqual(numeric_mismatches(a, b), [])
        b[("chiq.npz", "chiq")] = np.full(3, 2.0 + 1e-12)
        b[("chiq.npz", "n")] = np.int64(4)
        b[("gap.npz", "gap")] = v[1]
        c = dict(b)
        del c[("return", "lambda")]
        self.assertEqual(numeric_mismatches(a, b),
                         [("chiq.npz", "chiq"), ("chiq.npz", "n"), ("gap.npz", "gap")])
        self.assertIn(("return", "lambda"), numeric_mismatches(a, c))


if __name__ == "__main__":
    unittest.main()
