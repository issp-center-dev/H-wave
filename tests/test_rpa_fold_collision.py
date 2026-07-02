"""Regression test: RPA sublattice fold must accumulate wrap-around collisions.

When the hopping range reaches the folded supercell size, distinct original
R-vectors (e.g. +r and -r) fold onto the same (R, orbital) key.  UHFk
canonicalizes the folded R to [0, n) and SUMS such contributions; RPA's
Interaction._reshape_interaction must do the same -- keeping negative R keys
means the consumer's numpy scatter (tab_r[ir] = v, negative index wrapping)
silently overwrites one contribution with the other.
"""
import unittest

import numpy as np

from hwave.solver.rpa import Interaction


class _FakeLattice:
    def __init__(self, shape, subshape):
        self.shape = shape
        self.subshape = subshape


class TestRPAFoldCollision(unittest.TestCase):
    def _make_interaction(self, geom_norb_orig, subshape, shape):
        obj = object.__new__(Interaction)
        obj.enable_spin_orbital = False
        obj.lattice = _FakeLattice(shape=shape, subshape=subshape)
        obj.param_ham_orig = {"Geometry": {"norb": geom_norb_orig}}
        return obj

    def test_wraparound_collisions_are_summed(self):
        # CellShape [4,1,1], SubShape [2,1,1] -> folded shape nx=2.
        # Hops r=+2 and r=-2 are distinct original R-vectors that fold onto
        # the same canonical (R=1, orbital) keys.
        obj = self._make_interaction(1, subshape=(2, 1, 1), shape=(2, 1, 1))
        ham = {
            ((2, 0, 0), (0, 0)): 1.0,
            ((-2, 0, 0), (0, 0)): 0.5,
        }
        out = obj._reshape_interaction(ham, enable_spin_orbital=False)

        # all folded R components canonical (non-negative, within shape)
        for (ir, _ov) in out.keys():
            self.assertTrue(0 <= ir[0] < 2 and ir[1] == 0 and ir[2] == 0,
                            "folded R must be canonicalized, got {}".format(ir))

        # collisions summed, not overwritten
        self.assertTrue(
            np.isclose(out[((1, 0, 0), (0, 0))], 1.5),
            "colliding +r/-r contributions must be summed: got {}".format(
                out.get(((1, 0, 0), (0, 0)))))

        # total weight preserved: each original entry appears subvol times
        total_in = sum(ham.values()) * 2
        total_out = sum(out.values())
        self.assertTrue(np.isclose(total_out, total_in),
                        "fold must preserve the total interaction weight")


if __name__ == "__main__":
    unittest.main()
