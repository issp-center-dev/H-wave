"""Tests of ``hwave.initial.cells``: commensurate q, rationalisation,
representability on a SubShape and cell proposals (spec sections 5 step 2
and 5.4; validation cases V2a and V3(i)).

A q is three reduced fractions in [0, 1) of the reciprocal basis. UHFk
builds only DIAGONAL supercells (``SubShape = (Bx, By, Bz)``), so q is
representable on B iff each reduced denominator divides the matching B_i;
a smaller non-diagonal cell (the triangular sqrt3 x sqrt3 cell at K) is
named, never emitted.
"""
import unittest
from fractions import Fraction as F

from hwave.initial import cells
from hwave.initial.cells import (
    CellProposal, Rationalised, check_subshapes, commensurate_q,
    is_self_conjugate, order, propose, q_string, rationalise, representable,
    same_q, to_qfrac)

K = (F(2, 3), F(1, 3), F(0))


class TestQBasics(unittest.TestCase):

    def test_commensurate_2x2(self):
        qs = commensurate_q((2, 2, 1))
        self.assertEqual(qs, (
            (F(0), F(0), F(0)), (F(0), F(1, 2), F(0)),
            (F(1, 2), F(0), F(0)), (F(1, 2), F(1, 2), F(0))))
        for q in qs:
            self.assertTrue(all(isinstance(x, F) for x in q))
            self.assertTrue(representable(q, (2, 2, 1))[0])

    def test_commensurate_count_and_order(self):
        qs = commensurate_q([3, 2, 2])
        self.assertEqual(len(qs), 12)
        self.assertEqual(list(qs), sorted(qs))
        self.assertEqual(len(set(qs)), 12)
        self.assertEqual(commensurate_q((1, 1, 1)), ((F(0), F(0), F(0)),))

    def test_commensurate_refuses_bad_subshape(self):
        for bad in ((0, 1, 1), (2, 2), (2.0, 1, 1), (True, 1, 1), (-2, 1, 1)):
            with self.assertRaises(ValueError, msg=repr(bad)):
                commensurate_q(bad)

    def test_order_and_self_conjugate(self):
        q = (F(1, 2), F(1, 3), F(0))
        self.assertEqual(order(q), 6)
        self.assertFalse(is_self_conjugate(q))
        self.assertTrue(is_self_conjugate((F(1, 2), F(1, 2), F(0))))
        self.assertTrue(is_self_conjugate((0, 0, 0)))
        self.assertEqual(order((0, 0, 0)), 1)
        self.assertEqual(order(K), 3)
        self.assertFalse(is_self_conjugate(K))

    def test_to_qfrac_reduces_mod_one(self):
        self.assertEqual(to_qfrac(("-1/2", "4/3", 2)),
                         (F(1, 2), F(1, 3), F(0)))
        self.assertEqual(to_qfrac((F(-1, 3), 0, "0")), (F(2, 3), F(0), F(0)))
        with self.assertRaises(TypeError):
            to_qfrac((0.5, 0, 0))          # floats go through rationalise
        with self.assertRaises(ValueError):
            to_qfrac((F(1, 2), 0))
        with self.assertRaises(TypeError):
            to_qfrac((True, 0, 0))

    def test_same_q_mod_one(self):
        # PR 1a q_frac uses the representative in (-B/2, B/2]
        self.assertTrue(same_q((F(-1, 3), F(1, 3), F(0)), K))
        self.assertTrue(same_q(("1/2", "-1/2", "0"), ("-1/2", "1/2", "1")))
        self.assertFalse(same_q(K, (F(1, 3), F(2, 3), F(0))))

    def test_q_string(self):
        self.assertEqual(q_string(K), "2/3,1/3,0")
        self.assertEqual(q_string((F(-1, 2), 0, 0)), "1/2,0,0")


class TestRationalise(unittest.TestCase):

    def test_exact_within_tolerance(self):
        r = rationalise((2 / 3, 1 / 3 + 5e-7, 0.0))
        self.assertIsInstance(r, Rationalised)
        self.assertFalse(r.incommensurate)
        self.assertEqual(r.q, K)
        self.assertEqual(r.approximants, ())
        self.assertIsNone(r.reason)

    def test_mod_one_and_wraparound(self):
        r = rationalise((-0.5, 1.0 - 1e-9, -1e-9))
        self.assertEqual(r.q, (F(1, 2), F(0), F(0)))

    def test_v3i(self):
        # V3 (i): q = (7/16, 1/2) gives SubShape 16x2, which divides 64x32
        r = rationalise((7 / 16, 0.5, 0.0))
        self.assertFalse(r.incommensurate)
        self.assertEqual(r.q, (F(7, 16), F(1, 2), F(0)))
        prop = propose(r.q, (64, 32, 1))
        self.assertEqual(prop.proposed_subshape, (16, 2, 1))
        self.assertIsNone(prop.proposed_cellshape)
        check_subshapes([prop.proposed_subshape], (64, 32, 1))
        self.assertEqual(representable(r.q, (16, 2, 1)), (True, ""))
        # ... with max_denominator = 12 it is not representable
        r12 = rationalise((7 / 16, 0.5, 0.0), max_denominator=12)
        self.assertTrue(r12.incommensurate)
        self.assertIsNone(r12.q)
        self.assertEqual(r12.approximants, (
            (F(4, 9), F(1, 2), F(0)), (F(3, 7), F(1, 2), F(0))))
        self.assertIn("not representable within max_denominator = 12",
                      r12.reason)
        self.assertIn("incommensurate: approximated", r12.reason)
        self.assertIn("4/9,1/2,0", r12.reason)

    def test_default_max_denominator_is_16(self):
        self.assertTrue(rationalise((1 / 17, 0, 0)).incommensurate)
        self.assertFalse(rationalise((1 / 16, 0, 0)).incommensurate)

    def test_approximants_bracket_each_component(self):
        r = rationalise((0.3, 0.3, 0.0), max_denominator=2)
        # 0.3 lies between 0 and 1/2 on both axes; nearest first
        self.assertEqual(r.approximants[0], (F(1, 2), F(1, 2), F(0)))
        self.assertEqual(set(r.approximants), {
            (a, b, F(0)) for a in (F(0), F(1, 2)) for b in (F(0), F(1, 2))})
        # the upper neighbour of 0.9 at max_denominator 2 is 1 == 0 (mod 1)
        r = rationalise((0.9, 0, 0), max_denominator=2)
        self.assertEqual(r.approximants,
                         ((F(0), F(0), F(0)), (F(1, 2), F(0), F(0))))

    def test_rationalise_refuses_bad_input(self):
        with self.assertRaises(ValueError):
            rationalise((float("nan"), 0, 0))
        with self.assertRaises(ValueError):
            rationalise((0.5, 0))
        with self.assertRaises(ValueError):
            rationalise((0.5, 0, 0), max_denominator=0)
        with self.assertRaises(ValueError):
            rationalise((0.5, 0, 0), tol=-1)


class TestRepresentability(unittest.TestCase):

    def test_v2a_reason_and_proposal(self):
        ok, reason = representable(K, (2, 2, 1))
        self.assertFalse(ok)
        self.assertEqual(reason, "q_x denominator 3 does not divide B_x = 2")
        prop = propose(K, (12, 12, 1))
        self.assertIsInstance(prop, CellProposal)
        self.assertEqual(prop.proposed_subshape, (3, 3, 1))
        self.assertIsNone(prop.proposed_cellshape)
        self.assertIn("√3×√3", prop.rotated_note)
        self.assertIn("rebuilt unit cell", prop.rotated_note)
        self.assertEqual(representable(K, (3, 3, 1)), (True, ""))
        # PatternSpec strings are accepted as they are
        self.assertFalse(representable(("2/3", "1/3", "0"), [2, 2, 1])[0])

    def test_reason_names_first_failing_axis(self):
        self.assertEqual(
            representable((F(1, 2), F(1, 4), F(0)), (2, 2, 1)),
            (False, "q_y denominator 4 does not divide B_y = 2"))
        self.assertEqual(
            representable((0, 0, F(1, 2)), (2, 2, 3)),
            (False, "q_z denominator 2 does not divide B_z = 3"))

    def test_no_rotated_note_for_diagonal_minimum(self):
        # order(q) equals the diagonal volume: no smaller cell exists
        for q in ((F(1, 2), F(0), F(0)), (F(1, 3), 0, 0),
                  (F(1, 2), F(1, 3), F(0))):
            self.assertIsNone(propose(q, (48, 48, 1)).rotated_note, q)

    def test_rotated_note_general_index(self):
        # order 4 < 4*4: a smaller non-diagonal cell exists, but it is not
        # the triangular one
        prop = propose((F(1, 4), F(1, 4), F(0)), (8, 8, 1))
        self.assertEqual(prop.proposed_subshape, (4, 4, 1))
        self.assertIn("4 primitive cells", prop.rotated_note)
        self.assertIn("4x4x1", prop.rotated_note)
        self.assertNotIn("√3", prop.rotated_note)
        # Neel: the sqrt2 x sqrt2 cell has 2 primitive cells
        self.assertIn("2 primitive cells",
                      propose((F(1, 2), F(1, 2), F(0)), (8, 8, 1)).rotated_note)

    def test_lcm_growth(self):
        prop = propose((F(1, 16), F(0), F(0)), (12, 12, 1))
        self.assertEqual(prop.proposed_subshape, (16, 1, 1))
        self.assertEqual(prop.proposed_cellshape, (48, 12, 1))
        self.assertEqual(cells.lcm(12, 16), 48)
        prop = propose(K, (10, 12, 1))
        self.assertEqual(prop.proposed_cellshape, (30, 12, 1))

    def test_subshape_not_dividing_refused(self):
        with self.assertRaises(ValueError) as cm:
            check_subshapes([(2, 2, 1), (3, 3, 1)], (8, 8, 1))
        msg = str(cm.exception)
        self.assertIn("SubShape [3, 3, 1] is not compatible with "
                      "CellShape [8, 8, 1]", msg)
        self.assertIn("multiple of", msg)
        self.assertIn("CellShape_x = 8 is not a multiple of SubShape_x = 3",
                      msg)
        check_subshapes([(2, 2, 1), [4, 4, 1], (8, 8, 1)], (8, 8, 1))

    def test_check_subshapes_refuses_malformed(self):
        for bad in ([], [(2, 2)], [(0, 1, 1)], [(2.0, 1, 1)]):
            with self.assertRaises(ValueError, msg=repr(bad)):
                check_subshapes(bad, (8, 8, 1))
        with self.assertRaises(ValueError):
            check_subshapes([(1, 1, 1)], (8, 8))


if __name__ == "__main__":
    unittest.main()
