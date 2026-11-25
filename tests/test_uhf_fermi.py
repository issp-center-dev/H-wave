"""The Fermi function of UHFk and UHFr is ``hartree_fock.masked_fermi``.

The ``_legacy_*`` functions are the bodies the two solvers carried before
(copied verbatim; do not edit them). The shared function must reproduce
their arithmetic bit for bit with the solver's ``ene_cutoff``.
"""
import unittest

import numpy as np

from hwave.solver.hartree_fock import masked_fermi


def _legacy_uhfk_fermi(t, mu, ev, ene_cutoff):
    w_ = (ev - mu) / t
    mask_ = w_ < ene_cutoff
    w1_ = np.where( mask_, w_, 0.0 )
    v1_ = 1.0 / (1.0 + np.exp(w1_))
    v_ = np.where( mask_, v1_, 0.0 )
    return v_


def _same_bits(a, b):
    return a.dtype == b.dtype and a.shape == b.shape and a.tobytes() == b.tobytes()


def _eigenvalues():
    rng = np.random.default_rng(20261125)
    return rng.uniform(-6.0, 6.0, size=(16, 8))           # (nvol, nd)


class TestUHFkFermi(unittest.TestCase):

    def test_equals_the_former_uhfk_closure_bit_for_bit(self):
        ev = _eigenvalues()
        for t in (0.3, 0.05, 0.005):
            for mu in (-1.3, 0.0, 0.7):
                for cutoff in (20.0, 1.0e2):
                    with self.subTest(t=t, mu=mu, cutoff=cutoff):
                        self.assertTrue(_same_bits(
                            masked_fermi(t, mu, ev, cutoff),
                            _legacy_uhfk_fermi(t, mu, ev, cutoff)))

    def test_zero_at_the_cutoff(self):
        # (ev - mu) / t = 20 exactly
        ev = np.array([19.5, 20.0, 20.5])
        np.testing.assert_array_equal(
            masked_fermi(1.0, 0.0, ev, 20.0),
            _legacy_uhfk_fermi(1.0, 0.0, ev, 20.0))
        self.assertEqual(masked_fermi(1.0, 0.0, ev, 20.0)[1], 0.0)


if __name__ == "__main__":
    unittest.main()
