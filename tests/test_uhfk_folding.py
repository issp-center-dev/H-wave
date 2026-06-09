#!/usr/bin/env python3
"""Regression test for sublattice folding accumulation in UHFk.

When two distinct original (R, orbital) entries fold onto the same supercell
(ir, orbital) key -- which can happen when the hopping range reaches the cell
size (allowed under relax_checks) -- their amplitudes must be SUMMED, not
overwritten.
"""

import unittest

from hwave.solver.uhfk import UHFk


class TestReshapeInteractionAccumulation(unittest.TestCase):
    def _make_stub(self):
        solver = UHFk.__new__(UHFk)
        # A 1x1x1 block with no sublattice: every R-vector wraps to ir=(0,0,0),
        # so two different original R-vectors collide on the same folded key.
        solver.subshape = (1, 1, 1)
        solver.shape = (1, 1, 1)
        solver.norb_orig = 1
        solver.param_ham_orig = {"Geometry": {"norb": 1}}
        return solver

    def test_colliding_entries_are_summed(self):
        solver = self._make_stub()
        ham = {
            ((0, 0, 0), (0, 0)): 2.0,
            ((1, 0, 0), (0, 0)): 3.0,
        }
        out = solver._reshape_interaction(ham, False)
        # Both entries fold to the same key and must be accumulated: 2 + 3 = 5.
        self.assertEqual(len(out), 1)
        self.assertAlmostEqual(out[((0, 0, 0), (0, 0))], 5.0, places=12)


if __name__ == "__main__":
    unittest.main()
