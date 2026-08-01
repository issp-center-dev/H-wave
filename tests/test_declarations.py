"""The single-sourced declaration symmetrisation (#108 increment 2).

The module owns the reduction in both representations its consumers use
(real-space dense for the ring assembly, momentum-space for the S/C
builders); these tests pin the algebra once -- most importantly the
FFT-commutation identity that makes the two forms the SAME reduction --
plus bit-parity of each delegation against the pre-refactor inline
expressions.
"""

import unittest

import numpy as np

from hwave.solver.declarations import symmetrise_dense, symmetrise_k
from hwave.solver.kgrid import reverse_fft_axes


def _rand(shape, seed):
    rng = np.random.default_rng(seed)
    return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)


def _fft_ab_q(arr):
    """(nx, ny, nz, a, b) real-space table -> (a, b, nx, ny, nz) with the
    e^{-iqR} phase convention (numpy's forward FFT)."""
    return np.fft.fftn(arr, axes=(0, 1, 2)).transpose(3, 4, 0, 1, 2)


class TestFftCommutation(unittest.TestCase):
    """fft(symmetrise_dense(T)) == symmetrise_k(fft(T)): the identity
    that makes the two representations one reduction. Odd and even
    lattice dimensions, all-complex tables, both rules."""

    def test_same_operator_rule_commutes(self):
        arr = _rand((4, 3, 5, 3, 3), 1)
        lhs = _fft_ab_q(symmetrise_dense(arr))
        rhs = symmetrise_k({"Hund": _fft_ab_q(arr)})["Hund"]
        np.testing.assert_allclose(lhs, rhs, atol=1e-12)

    def test_pairhop_hermitian_rule_commutes(self):
        arr = _rand((4, 3, 5, 3, 3), 2)
        lhs = _fft_ab_q(symmetrise_dense(arr, hermitian=True))
        rhs = symmetrise_k({"PairHop": _fft_ab_q(arr)})["PairHop"]
        np.testing.assert_allclose(lhs, rhs, atol=1e-12)


class TestDenseForm(unittest.TestCase):

    def test_matches_the_pre_refactor_inline_expression(self):
        """Bit-parity: the ring assembly's historical inline expression
        (reverse + orbital transpose (+ conjugate for PairHop), mean)."""
        for herm, seed in ((False, 3), (True, 4)):
            with self.subTest(hermitian=herm):
                arr = _rand((4, 3, 5, 2, 2), seed)
                rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                                   (0, 1, 2, 4, 3))
                if herm:
                    rev = np.conjugate(rev)
                old = 0.5 * (arr + rev)
                np.testing.assert_array_equal(
                    symmetrise_dense(arr, hermitian=herm), old)

    def test_one_sided_bond_splits_between_both_directions(self):
        """A single +x declaration means (T_ab + T_ba)/2: half the weight
        stays, half lands at the wrapped -x cell on the transposed
        orbital pair -- the reading the ring adjudicated (a one-sided
        off-site V has vertex v cos(qR), not v e^{-iqR})."""
        nx = 4
        arr = np.zeros((nx, 1, 1, 2, 2), dtype=complex)
        arr[1, 0, 0, 0, 1] = 1.0
        sym = symmetrise_dense(arr)
        self.assertAlmostEqual(sym[1, 0, 0, 0, 1].real, 0.5, places=15)
        self.assertAlmostEqual(sym[nx - 1, 0, 0, 1, 0].real, 0.5,
                               places=15)
        self.assertAlmostEqual(float(np.abs(sym).sum()), 1.0, places=15)

    def test_antisymmetric_declaration_vanishes(self):
        arr = np.zeros((1, 1, 1, 2, 2), dtype=complex)
        arr[0, 0, 0, 0, 1] = 1.0
        arr[0, 0, 0, 1, 0] = -1.0
        self.assertEqual(float(np.abs(symmetrise_dense(arr)).max()), 0.0)

    def test_hermitian_closed_pairhop_keeps_its_phase(self):
        p = 0.7 + 0.4j
        arr = np.zeros((1, 1, 1, 2, 2), dtype=complex)
        arr[0, 0, 0, 0, 1] = p
        arr[0, 0, 0, 1, 0] = np.conj(p)
        sym = symmetrise_dense(arr, hermitian=True)
        self.assertAlmostEqual(sym[0, 0, 0, 0, 1], p, places=15)


class TestKForm(unittest.TestCase):

    def test_matches_the_pre_refactor_inline_expression(self):
        """Bit-parity of the moved k-space body: fancy-gather reversal +
        orbital transpose mean; PairHop conjugate-transpose at fixed q."""
        M = _rand((2, 2, 4, 3, 5), 5)
        nx, ny, nz = 4, 3, 5
        Mrev = M[:, :, (-np.arange(nx)) % nx][:, :, :, (-np.arange(ny)) % ny][
            :, :, :, :, (-np.arange(nz)) % nz]
        old_plain = 0.5 * (M + Mrev.transpose(1, 0, 2, 3, 4))
        old_ph = 0.5 * (M + np.conj(M.transpose(1, 0, 2, 3, 4)))
        out = symmetrise_k({"Exchange": M, "PairHop": M})
        np.testing.assert_array_equal(out["Exchange"], old_plain)
        np.testing.assert_array_equal(out["PairHop"], old_ph)

    def test_idempotent(self):
        M = _rand((2, 2, 4, 3, 5), 6)
        once = symmetrise_k({"Hund": M, "PairHop": M})
        twice = symmetrise_k(once)
        for k in once:
            np.testing.assert_allclose(twice[k], once[k], atol=1e-15)


if __name__ == "__main__":
    unittest.main()
