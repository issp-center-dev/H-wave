"""Gate C1b-PD (issue #270): the spin-orbital bond FLEX in the S_z-packed
layout against the dense layout (``[mode.param] so_bond_layout = "dense"``)
on the S_z-conserving ``ZEEMAN_SZ_ALL`` (a Zeeman field, every interaction
type on- and off-site, a complex off-site PairHop).

* :class:`TestOneMap` -- one map at a fixed chemical potential: Sigma_HF,
  Sigma_fluct, chi0q and chiq of the two layouts agree within ``1e-12``;
  the two runs really use different layouts (the layout reason lines and
  the store types of ``so_bond.open_store``); the compared content is
  off-site and spin polarized; a broken packed dressing misses.
* :class:`TestConverged` -- the converged states (its fast smoke is
  :class:`TestOneMap`): Sigma, G, mu, NCond within ``max(10 eps, 1e-8)``.

Mutations patch module attributes looked up at call time; every reference
is computed before its patch is entered.
"""
import unittest
from unittest import mock

import numpy as np

from tests import so_bond_flex_fixtures as X
from tests.so_bond_oracles import rel

#: the array tolerance of the one-map gate
TOL = 1e-12
#: the miss floor of a deliberate break
MISS = 1e-3
#: the off-site content floor of Sigma_fluct
CONTENT = 0.1
#: the spin-polarization floor: |Sigma_upup - Sigma_dndn| >= SPIN_POL |Sigma|
SPIN_POL = 1e-2
#: the SCF threshold of the converged comparison: max(10 eps, 1e-8), EPS = 12
CONV_TOL = max(10 * 1e-12, 1e-8)

FX = X.ZEEMAN_SZ_ALL
DENSE = {"so_bond_layout": "dense"}
PACKED_REASON = ("so bond layout: sz_packed (S_z conserved: no spin-flip transfer, "
                 "no S_z-changing interaction)")
DENSE_REASON = 'so bond layout: dense (reason: forced by so_bond_layout = "dense")'
#: the payload entries compared by the one-map gate
ROWS = (("Sigma_HF", "static_new"), ("Sigma_fluct", "fluct_new"), ("chi0q", "chi0q"),
        ("chiq", "chiq_s"))


def _run(param=None):
    """The :class:`tests.so_bond_flex_fixtures.RunResult` of one map of
    ``FX`` (``param`` overrides ``[mode.param]``) and the types of the
    stores ``so_bond.open_store`` returned during it."""
    from hwave.solver import so_bond
    orig = so_bond.open_store
    types = []

    def spy(*args, **kwargs):
        store = orig(*args, **kwargs)
        types.append(type(store))
        return store
    par = dict(param or {}, IterationMax=1)
    with mock.patch.object(so_bond, "open_store", spy):
        res = X.run_flex(FX, param=par)
    if len(res.payloads) != 1:
        raise RuntimeError("_run: {} payloads".format(len(res.payloads)))
    return res, types


def _messages(records):
    return [r.getMessage() for r in records]


def _zero_minus_sector(orig):
    """A wrapper of ``so_bond.so_bond_dress`` that zeroes the ``"-1"``
    sector block of the packed result (a deliberate break); ``applied``
    counts the zeroed calls."""
    from hwave.solver import so_bond

    def wrapper(Xb, gamma_b, **kwargs):
        Y = orig(Xb, gamma_b, **kwargs)
        if isinstance(Y, so_bond.PackedArrays):
            Y.blocks[2][...] = 0.0
            wrapper.applied += 1
        return Y
    wrapper.applied = 0
    return wrapper


class TestOneMap(unittest.TestCase):
    """C1b-PD, one map: packed and dense layouts of ``ZEEMAN_SZ_ALL``."""

    @classmethod
    def setUpClass(cls):
        from hwave.solver import flex_bond, so_bond
        cls.packed, cls.packed_types = _run()
        cls.dense, cls.dense_types = _run(DENSE)
        cls.packed_store = so_bond.PackedBondBlockStore
        cls.dense_store = flex_bond.BondBlockStore

    def test_layouts_agree(self):
        p, d = self.packed.payloads[0], self.dense.payloads[0]
        self.assertEqual(self.packed.solver._so_admission.layout, "sz_packed")
        self.assertEqual(self.dense.solver._so_admission.layout, "dense")
        self.assertEqual(p["mu"], d["mu"])
        for label, key in ROWS:
            with self.subTest(row=label):
                self.assertEqual(p[key].shape, d[key].shape)
                self.assertTrue(np.all(np.isfinite(p[key])))
                self.assertTrue(np.all(np.isfinite(d[key])))
                self.assertGreater(float(np.abs(d[key]).max()), 0.0)
                dev = rel(p[key], d[key])
                print("C1b-PD one map %-11s %.3e" % (label, dev))
                self.assertLessEqual(dev, TOL)
        # the two runs really differ in layout: reason lines and store types
        pm, dm = _messages(self.packed.records), _messages(self.dense.records)
        self.assertEqual(pm.count(PACKED_REASON), 1)
        self.assertEqual(dm.count(DENSE_REASON), 1)
        self.assertNotIn(DENSE_REASON, pm)
        self.assertNotIn(PACKED_REASON, dm)
        self.assertEqual(self.packed_types, [self.packed_store])
        self.assertEqual(self.dense_types, [self.dense_store])
        self.assertIsNot(self.packed_store, self.dense_store)

    def test_non_vacuity(self):
        d = self.dense.payloads[0]
        n = FX.norb_phys
        fluct = d["fluct_new"]
        off_site = float(np.abs(fluct - fluct.mean(axis=-3, keepdims=True)).max()
                         / np.abs(fluct).max())
        print("C1b-PD off-site share of Sigma_fluct %.3e" % off_site)
        self.assertGreaterEqual(off_site, CONTENT)
        for label, key in (("Sigma_HF", "static_new"), ("Sigma_fluct", "fluct_new")):
            a = d[key]
            pol = float(np.abs(a[..., :n, :n] - a[..., n:, n:]).max() / np.abs(a).max())
            print("C1b-PD spin polarization of %s %.3e" % (label, pol))
            self.assertGreaterEqual(pol, SPIN_POL, label)

    def test_mutation(self):
        from hwave.solver import so_bond
        ref = self.dense.payloads[0]                     # the reference (dense run) first
        wrapper = _zero_minus_sector(so_bond.so_bond_dress)
        spy = mock.Mock(wraps=wrapper)
        self.assertEqual(spy.call_count, 0)
        with mock.patch.object(so_bond, "so_bond_dress", spy):
            res, types = _run()
        self.assertEqual(types, [self.packed_store])
        self.assertGreater(spy.call_count, 0)
        self.assertGreater(wrapper.applied, 0)           # the break was applied
        miss = rel(res.payloads[0]["fluct_new"], ref["fluct_new"])
        print("C1b-PD mutation Y@-1 = 0 miss %.3e" % miss)
        self.assertTrue(np.isfinite(miss))
        self.assertGreaterEqual(miss, MISS)


class TestConverged(unittest.TestCase):
    """C1b-PD, converged states (fast smoke :class:`TestOneMap`)."""

    def test_layouts_agree(self):
        p = X.converge(FX)
        d = X.converge(FX, param=DENSE)
        self.assertEqual(p.solver._so_admission.layout, "sz_packed")
        self.assertEqual(d.solver._so_admission.layout, "dense")
        a, b = p.green_info, d.green_info
        rows = {"Sigma": rel(a["sigma"], b["sigma"]),
                "G": rel(a["green"], b["green"]),
                "mu": abs(a["physics"]["mu"] - b["physics"]["mu"]),
                "NCond": abs(a["physics"]["NCond"] - b["physics"]["NCond"])}
        print("C1b-PD converged (%d / %d iterations) " % (
            p.solver.scf_iterations, d.solver.scf_iterations)
            + "  ".join("%s %.3e" % kv for kv in rows.items()))
        for k, v in rows.items():
            self.assertTrue(np.isfinite(v), k)
            self.assertLessEqual(v, CONV_TOL, k)


if __name__ == "__main__":
    unittest.main()
