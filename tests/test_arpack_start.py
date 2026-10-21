"""The Eliashberg leading-eigenvalue solve is reproducible (GitHub issue #255).

Without a start vector ARPACK draws a random one on every call, so the
converged eigenvalue of one input moved by ~1e-9 relative from run to run
(the in-process vs post-processing pairing comparison failed in 2 of 40
runs). Every unseeded call now starts from :func:`hwave.sc._arpack_start`.
"""
import unittest

import numpy as np
from scipy.sparse.linalg import LinearOperator

import hwave.sc as sc


def _random_operator(n=60, seed=3):
    rng = np.random.default_rng(seed)
    m = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    m = 0.5 * (m + m.conj().T) + 0.3j * np.diag(rng.normal(size=n))   # non-Hermitian
    return lambda: (LinearOperator((n, n), matvec=lambda v: m @ v, dtype=complex), n)


class TestArpackStart(unittest.TestCase):
    def test_start_vector_is_fixed_and_yields_to_a_seed(self):
        A, n = _random_operator()()
        v1, v2 = sc._arpack_start(A), sc._arpack_start(A)
        np.testing.assert_array_equal(v1, v2)
        self.assertAlmostEqual(float(np.linalg.norm(v1)), 1.0, places=12)
        seed = np.ones(n, dtype=complex)
        self.assertIs(sc._arpack_start(A, seed), seed)

    def test_start_vector_has_the_operator_precision(self):
        for dt in (np.float32, np.float64, np.complex64, np.complex128):
            with self.subTest(dtype=dt):
                A = LinearOperator((10, 10), matvec=lambda v: 2.0 * v, dtype=dt)
                self.assertEqual(sc._arpack_start(A).dtype, np.dtype(dt))
        # an integer operator keeps a float64 start (eigs promotes it to double)
        Ai = LinearOperator((10, 10), matvec=lambda v: 2 * v, dtype=np.int64)
        self.assertEqual(sc._arpack_start(Ai).dtype, np.dtype(np.float64))

    def test_odd_sector_leading_eigenvalue_is_found(self):
        """The reason the start is generic and not constant: on an operator
        that commutes with the reflection i -> n-1-i (a stand-in for a
        lattice symmetry), a constant start lies entirely in the EVEN
        sector, and so does the whole Krylov space it starts; the leading
        eigenvalue here sits in the ODD sector."""
        n = 40
        J = np.eye(n)[::-1]
        rng = np.random.default_rng(5)
        h = rng.normal(size=(n, n)); h = 0.5 * (h + h.T)
        P_even, P_odd = 0.5 * (np.eye(n) + J), 0.5 * (np.eye(n) - J)
        m = P_even @ h @ P_even + P_odd @ (h + 10.0 * np.eye(n)) @ P_odd
        make = lambda: (LinearOperator((n, n), matvec=lambda v: m @ v, dtype=complex), n)
        A, _ = make()
        v = sc._arpack_start(A)
        self.assertGreater(float(np.linalg.norm(P_odd @ v)), 0.1)
        self.assertGreater(float(np.linalg.norm(P_even @ v)), 0.1)
        lam = sc._solve_leading(make, n, "arnoldi", num_eigenvalues=4,
                                spectral_shift="auto")[0]
        want = float(np.max(np.linalg.eigvalsh(m)))
        self.assertAlmostEqual(float(np.real(lam)), want, places=8)
        const = np.ones(n) / np.sqrt(n)                     # the rejected alternative
        self.assertLess(float(np.linalg.norm(P_odd @ const)), 1e-12)

    def test_repeated_solves_are_bit_identical(self):
        make = _random_operator()
        cases = (("arnoldi", {}), ("arnoldi", {"spectral_shift": "auto"}),
                 ("shift-invert-bicgstab", {}), ("shift-invert-bicgstab", {"sigma_shift": 2.0}))
        for mode, kw in cases:
            with self.subTest(mode=mode, **kw):
                a = sc._solve_leading(make, 60, mode, num_eigenvalues=4, **kw)
                lam_a, vec_a = complex(a[0]), np.array(a[1], copy=True)
                b = sc._solve_leading(make, 60, mode, num_eigenvalues=4, **kw)
                self.assertEqual(lam_a, complex(b[0]))
                np.testing.assert_array_equal(vec_a, b[1])

    def test_a_seed_is_used_and_left_unchanged(self):
        make = _random_operator()
        A, n = make()
        seed = np.random.default_rng(9).normal(size=n) + 0j
        keep = seed.copy()
        sc._solve_leading(make, n, "arnoldi", num_eigenvalues=4, seed_vec=seed)
        np.testing.assert_array_equal(seed, keep)


if __name__ == "__main__":
    unittest.main()
