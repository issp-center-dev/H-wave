"""The content-keyed ED result cache (``tests/ed_cache.py``) shares a
result only between bit-identical inputs, returns read-only arrays, and
serves exactly what the engine itself computes."""

import unittest

import numpy as np

from tests import ed_cache
from tests import ed_oracle_util as U


def _fx(mu=0.2):
    return U.EDFixture(L=2, norb=2, t={(0, 0): 0.7 + 0.3j, (1, 1): 0.5 - 0.2j,
                                       (0, 1): 0.25 + 0.15j, (1, 0): 0.25 + 0.15j},
                       eps=(0.10, -0.05), T=0.5, mu=mu)


def _terms(fx, v):
    return U.canonical_density_terms(fx, [(0, 0, 0, v), (0, 1, 1, 0.5 * v)])


CHANNELS = [(0, 0, 0), (1, 0, 1), (1, 1, 0)]


class TestEDCache(unittest.TestCase):

    def test_serves_the_engine_value_bit_for_bit(self):
        fx = _fx()
        terms = _terms(fx, 0.3)
        want = U.SectorED(fx, terms=terms).bond_correlator(CHANNELS)
        got = ed_cache.sector(fx, terms=terms).bond_correlator(CHANNELS)
        self.assertTrue(np.array_equal(got, want))
        iws = 1j * (2 * np.arange(-4, 4) + 1) * np.pi * fx.T
        self.assertTrue(np.array_equal(
            ed_cache.total(fx, terms=terms).green(iws),
            U.TotalNED(fx, terms=terms).green(iws)))

    def test_identical_content_is_shared_and_read_only(self):
        terms = _terms(_fx(), 0.2)
        a = ed_cache.sector(_fx(), terms=terms).pair_correlator(CHANNELS)
        before = ed_cache.stats().get("result_hit", 0)
        # a DIFFERENT but equal fixture object, a copied term list, an
        # explicit h1 equal to the default one: still the same content
        fx2 = _fx()
        b = ed_cache.sector(fx2, terms=list(terms), h1=fx2.build_h1().copy()
                            ).pair_correlator(list(CHANNELS))
        self.assertIs(a, b)
        self.assertEqual(ed_cache.stats()["result_hit"], before + 1)
        self.assertFalse(a.flags.writeable)
        with self.assertRaises(ValueError):
            a[0, 0, 0] = 0.0

    def test_any_difference_in_the_inputs_is_a_different_entry(self):
        fx = _fx()
        terms = _terms(fx, 0.2)
        base = ed_cache.sector(fx, terms=terms).bond_correlator(CHANNELS)
        h1 = fx.build_h1().copy()
        h1[0, 0] = np.nextafter(h1[0, 0].real, 1.0)                     # one ulp
        variants = [
            ed_cache.sector(_fx(mu=np.nextafter(0.2, 1.0)), terms=terms),  # mu
            ed_cache.sector(fx, terms=terms[::-1]),                       # term order
            ed_cache.sector(fx, terms=_terms(fx, np.nextafter(0.2, 1.0))),  # a coefficient
            ed_cache.sector(fx, terms=terms, h1=h1),                      # h1
            ed_cache.total(fx, terms=terms),                              # engine
        ]
        keys = {v.key for v in variants} | {ed_cache.sector(fx, terms=terms).key}
        self.assertEqual(len(keys), len(variants) + 1)
        before = ed_cache.stats().get("result_miss", 0)
        out = variants[1].bond_correlator(CHANNELS)
        self.assertEqual(ed_cache.stats()["result_miss"], before + 1)
        self.assertIsNot(out, base)

    def test_zero_coefficient_terms_leave_the_engine_bit_identical(self):
        """The one key normalisation (dropping exactly-zero terms) is exact:
        the engines with and without them agree bit for bit."""
        fx = _fx()
        zeros = _terms(fx, 0.0) + [(1, 1, 2, 2, -0.0j), (3, 3, 0, 0, 0.0)]
        mixed = zeros[:3] + _terms(fx, 0.3) + zeros[3:]
        for name, cls in (("SectorED", U.SectorED), ("TotalNED", U.TotalNED)):
            for with_zeros, without in ((zeros, []), (mixed, _terms(fx, 0.3))):
                a, b = cls(fx, terms=with_zeros), cls(fx, terms=without)
                for key in a._V:
                    self.assertEqual(a._V[key].tobytes(), b._V[key].tobytes(), name)
                    self.assertEqual(a._ev[key].tobytes(), b._ev[key].tobytes(), name)
                self.assertEqual(
                    a.bond_correlator_transverse(CHANNELS).tobytes(),
                    b.bond_correlator_transverse(CHANNELS).tobytes(), name)
                self.assertEqual(ed_cache.sector(fx, with_zeros).key,
                                 ed_cache.sector(fx, without).key)

    def _mutate_before_first_use(self, explicit):
        """Construct a stand-in, mutate the h1 it was given (explicitly, or
        the array EDFixture.build_h1 caches) BEFORE the lazy engine build,
        and check the entry still holds the original content."""
        fx = _fx(mu=0.2 + (1e-3 if explicit else 2e-3))   # a fresh key per case
        terms = _terms(fx, 0.25)
        original = fx.build_h1().copy()
        want = U.SectorED(fx, terms=terms, h1=original.copy()).bond_correlator(CHANNELS)
        if explicit:
            src = original.copy()
            cached = ed_cache.sector(fx, terms=terms, h1=src)
        else:
            src = fx.build_h1()                 # the fixture's own cached array
            cached = ed_cache.sector(fx, terms=terms)
        src[0, 2] += 0.5                       # mutate the source
        src[2, 0] += 0.5
        try:
            first = cached.bond_correlator(CHANNELS)
            self.assertTrue(np.array_equal(first, want))
            again = ed_cache.sector(fx, terms=terms, h1=original.copy()
                                    ).bond_correlator(CHANNELS)
            self.assertTrue(np.array_equal(again, want))
        finally:
            if not explicit:
                src[0, 2] -= 0.5
                src[2, 0] -= 0.5

    def test_mutating_an_explicit_h1_after_construction_is_harmless(self):
        self._mutate_before_first_use(explicit=True)

    def test_mutating_the_fixture_h1_after_construction_is_harmless(self):
        self._mutate_before_first_use(explicit=False)

    def test_only_public_observables_are_served(self):
        with self.assertRaises(AttributeError):
            ed_cache.sector(_fx())._build_operator


if __name__ == "__main__":
    unittest.main()
