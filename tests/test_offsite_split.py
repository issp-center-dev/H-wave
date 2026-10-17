"""The pre-fold locality split has ONE owner, ``hwave.solver.offsite``;
the FLEX general path and the RPA bond-resolved longitudinal gate both
consume it (GitHub issue #181, Tier 3 Phase A).

The split is DECLARATION-based: every declared off-site entry of every
two-body type is retained (zero-valued entries, Exchange, PairHop and
PairLift included). Policy -- what to reject, skip or warn about -- is
the caller's, so this module can be shared by callers with different
policies without either one silently inheriting the other's.
"""
import unittest

import numpy as np


class _FakeLattice:
    def __init__(self, subshape=(1, 1, 1)):
        self.subshape = tuple(subshape)


class _FakeHamInfo:
    def __init__(self, param_ham, param_ham_orig=None):
        self.param_ham = param_ham
        self.param_ham_orig = param_ham_orig

    def _reshape_interaction(self, tbl, enable_spin_orbital):
        # a tagging "fold" so the tests can see which table was folded
        return {k: ("folded", v) for k, v in tbl.items()}


class TestSplitLocality(unittest.TestCase):

    def test_routes_every_declared_type_by_irvec(self):
        from hwave.solver.offsite import split_locality
        ham = {"CoulombInter": {((0, 0, 0), (0, 1)): 0.5,
                                ((1, 0, 0), (0, 0)): 0.2},
               "Exchange": {((1, 0, 0), (0, 1)): 0.0},
               "PairLift": {((1, 0, 0), (0, 1)): 0.1},
               "Geometry": {"norb": 2}}
        s = split_locality(_FakeHamInfo(ham), _FakeLattice())
        self.assertFalse(s.has_fold)
        self.assertEqual(s.onsite_tbl,
                         {"CoulombInter": {((0, 0, 0), (0, 1)): 0.5}})
        # declaration-based: the zero-valued Exchange and the inert
        # PairLift are KEPT; the caller decides what to do with them
        self.assertEqual(s.offsite_types,
                         ("CoulombInter", "Exchange", "PairLift"))
        self.assertEqual(s.offsite_tbl["Exchange"],
                         {((1, 0, 0), (0, 1)): 0.0})
        self.assertEqual(s.offsite_tbl["CoulombInter"],
                         {((1, 0, 0), (0, 0)): 0.2})
        # without folding the pre-fold off-site table IS the off-site table
        self.assertEqual(s.offsite_prefold_tbl, s.offsite_tbl)
        # and the whole table is the reader's own object (no copy, no union)
        self.assertIs(s.whole_tbl, ham)

    def test_fold_without_prefold_table_fails_closed(self):
        from hwave.solver.offsite import split_locality
        with self.assertRaises(ValueError) as cm:
            split_locality(_FakeHamInfo({"Hund": {}}, None),
                           _FakeLattice((2, 1, 1)))
        self.assertIn("param_ham_orig", str(cm.exception))

    def test_folding_folds_each_part_and_keeps_the_prefold_offsite_table(self):
        from hwave.solver.offsite import split_locality
        orig = {"Hund": {((0, 0, 0), (0, 1)): 0.4, ((2, 0, 0), (0, 0)): 0.3}}
        folded = {"Hund": {((0, 0, 0), (0, 1)): 0.4, ((0, 0, 0), (0, 0)): 0.3}}
        s = split_locality(_FakeHamInfo(folded, orig), _FakeLattice((2, 1, 1)))
        self.assertTrue(s.has_fold)
        self.assertEqual(s.onsite_tbl,
                         {"Hund": {((0, 0, 0), (0, 1)): ("folded", 0.4)}})
        self.assertEqual(s.offsite_tbl,
                         {"Hund": {((2, 0, 0), (0, 0)): ("folded", 0.3)}})
        self.assertEqual(s.offsite_prefold_tbl,
                         {"Hund": {((2, 0, 0), (0, 0)): 0.3}})
        self.assertIs(s.whole_tbl, folded)     # the reader's folded table

    def test_aggregate_coulomb_is_split_before_routing(self):
        from hwave.solver.offsite import split_locality
        ham = {"Coulomb": {((0, 0, 0), (0, 0)): 3.0,
                           ((1, 0, 0), (0, 0)): 0.3}}
        s = split_locality(_FakeHamInfo(ham), _FakeLattice())
        self.assertEqual(s.onsite_tbl,
                         {"CoulombIntra": {((0, 0, 0), (0, 0)): 3.0}})
        self.assertEqual(s.offsite_tbl,
                         {"CoulombInter": {((1, 0, 0), (0, 0)): 0.3}})
        self.assertIn("CoulombInter", s.whole_tbl)
        self.assertNotIn("Coulomb", s.whole_tbl)

    def test_aggregate_coulomb_with_explicit_table_is_rejected(self):
        from hwave.solver.offsite import normalize_coulomb
        with self.assertRaises(ValueError):
            normalize_coulomb({"Coulomb": {}, "CoulombIntra": {}})

    def test_sc_matrices_from_split_filters_offsite_to_density_types(self):
        """The whole table keeps Exchange (its ON-site cross content stays
        in the S/C matrices); the off-site part handed to the builder is
        filtered to the density types, so an off-site Exchange entry
        contributes nothing and the result equals the Exchange-free case."""
        from hwave.solver.offsite import (
            split_locality, sc_matrices_from_split)
        base = {"CoulombInter": {((1, 0, 0), (0, 1)): 0.2,
                                 ((-1, 0, 0), (1, 0)): 0.2}}
        with_x = dict(base)
        with_x["Exchange"] = {((1, 0, 0), (0, 1)): 0.3,
                              ((-1, 0, 0), (1, 0)): 0.3}
        s0 = split_locality(_FakeHamInfo(base), _FakeLattice())
        s1 = split_locality(_FakeHamInfo(with_x), _FakeLattice())
        S0, C0 = sc_matrices_from_split(s0, ("CoulombInter",), 2, 4, 1, 1)
        S1, C1 = sc_matrices_from_split(s1, ("CoulombInter",), 2, 4, 1, 1)
        self.assertEqual(S0.shape, (4, 1, 1, 4, 4))
        self.assertTrue(np.array_equal(S0, S1))
        self.assertTrue(np.array_equal(C0, C1))
        # the off-site V enters the density slots of C (aa,bb) with a
        # q-dependence, and nothing enters S
        self.assertEqual(np.abs(S0).max(), 0.0)
        self.assertGreater(np.abs(C0[:, 0, 0, 0, 3]).max(), 0.0)


if __name__ == "__main__":
    unittest.main()
