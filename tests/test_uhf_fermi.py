"""The Fermi function of UHFk and UHFr is ``hartree_fock.masked_fermi``.

The ``_legacy_*`` functions are the bodies the two solvers carried before
(copied verbatim; do not edit them). The shared function must reproduce
their arithmetic bit for bit with the solver's ``ene_cutoff``.
"""
import unittest

import numpy as np

from hwave.solver.hartree_fock import masked_fermi
from hwave.solver.uhfr import UHFr


def _legacy_uhfk_fermi(t, mu, ev, ene_cutoff):
    w_ = (ev - mu) / t
    mask_ = w_ < ene_cutoff
    w1_ = np.where( mask_, w_, 0.0 )
    v1_ = 1.0 / (1.0 + np.exp(w1_))
    v_ = np.where( mask_, v1_, 0.0 )
    return v_


def _legacy_uhfr_fermi(self, mu, eigenvalue):
    x = (eigenvalue - mu) / self.T
    fermi = np.where(x > self.ene_cutoff, 0.0, 1.0 / (np.exp(x) + 1.0))
    return fermi


class _Host:
    """The two attributes ``UHFr._fermi`` reads."""

    def __init__(self, T, ene_cutoff=1.0e2):
        self.T = T
        self.ene_cutoff = ene_cutoff


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


class TestUHFrFermi(unittest.TestCase):
    """The former UHFr expression masked with ``x > ene_cutoff`` and took
    ``exp`` of every entry. The shared function masks with
    ``x >= ene_cutoff`` before ``exp``. The two agree bit for bit except
    for an entry with ``(eigenvalue - mu) / T == ene_cutoff`` exactly."""

    def test_equals_the_former_uhfr_expression_bit_for_bit(self):
        ev = _eigenvalues()
        for t in (0.3, 0.05, 0.005, 0.002):
            for mu in (-1.3, 0.0, 0.7):
                for cutoff in (20.0, 1.0e2):
                    host = _Host(t, cutoff)
                    with self.subTest(t=t, mu=mu, cutoff=cutoff):
                        with np.errstate(over="ignore"):
                            legacy = _legacy_uhfr_fermi(host, mu, ev)
                        self.assertTrue(_same_bits(
                            UHFr._fermi(host, mu, ev), legacy))

    def test_no_overflow_at_low_temperature(self):
        # (eigenvalue - mu) / T reaches 3000; exp overflows beyond 709.78
        host = _Host(0.002)
        ev = _eigenvalues()
        with self.assertRaises(FloatingPointError):
            with np.errstate(over="raise"):
                _legacy_uhfr_fermi(host, 0.0, ev)
        with np.errstate(over="raise"):
            f = UHFr._fermi(host, 0.0, ev)
        self.assertTrue(np.all((f >= 0.0) & (f <= 1.0)))

    def test_the_entry_exactly_at_the_cutoff(self):
        """The one difference: the former expression kept
        ``1 / (exp(ene_cutoff) + 1)`` (3.7e-44 for the default 100) for
        an entry exactly at the cutoff; it is now zero, as in UHFk."""
        host = _Host(1.0)
        ev = np.array([99.5, 100.0, 100.5])
        legacy = _legacy_uhfr_fermi(host, 0.0, ev)
        new = UHFr._fermi(host, 0.0, ev)
        self.assertEqual(new[0], legacy[0])
        self.assertEqual(legacy[1], 1.0 / (np.exp(100.0) + 1.0))
        self.assertLess(legacy[1], 1.0e-43)
        self.assertEqual(new[1], 0.0)
        self.assertEqual(new[2], 0.0)
        self.assertEqual(legacy[2], 0.0)


if __name__ == "__main__":
    unittest.main()
