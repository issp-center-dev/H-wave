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
        solver.norb_phys_orig = 1
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


class TestReshapeInteractionStrideSO(unittest.TestCase):
    """Physical-orbital-indexed interactions (CoulombIntra/Inter/Coulomb, Hund,
    Ising, Exchange, PairHop) fold with the physical orbital count as stride,
    not the SO-count Geometry norb. Under enable_spin_orbital=True with
    SubShape > (1,1,1), a mis-strided fold would emit orbital indices out of
    range of the downstream uab_r/vab_r arrays (which are sized norb_phys *
    subvol) and _make_ham_inter would raise IndexError. Mirrors the rpa.py
    regression `TestReshapeInteractionStride` from tests/test_rpa_geom_norb.py.
    """

    def _make_stub(self, subshape, cellshape, norb_orig, enable_spin_orbital):
        solver = UHFk.__new__(UHFk)
        solver.subshape = tuple(subshape)
        solver.shape = tuple(c // s for c, s in zip(cellshape, subshape))
        solver.norb_orig = norb_orig
        solver.enable_spin_orbital = enable_spin_orbital
        solver.norb_phys_orig = norb_orig // 2 if enable_spin_orbital else norb_orig
        solver.param_ham_orig = {"Geometry": {"norb": norb_orig}}
        return solver

    def test_physical_indexed_fold_uses_norb_phys_stride_in_so_mode(self):
        # SO mode, geom norb (SO count) = 2 -> norb_phys = 1. A physical-
        # orbital-indexed entry (enable_spin_orbital=False call, as for
        # CoulombIntra) on physical orbital 0, folded over subshape (2,2,1),
        # must land inside [0, norb_phys * subvol) = [0, 4), not stride by
        # the SO count (which would produce indices 0, 2, 4, 6 and later
        # crash _make_ham_inter with an out-of-bounds write into uab_r).
        solver = self._make_stub(subshape=(2, 2, 1), cellshape=(4, 4, 1),
                                  norb_orig=2, enable_spin_orbital=True)
        ham = {((0, 0, 0), (0, 0)): 1.0}
        out = solver._reshape_interaction(ham, enable_spin_orbital=False)
        folded_indices = [i for (_ir, ov) in out.keys() for i in ov]
        self.assertTrue(all(0 <= i < 4 for i in folded_indices),
                        "physical-indexed fold must stride by norb_phys under "
                        "enable_spin_orbital=True; got {}".format(folded_indices))

    def test_transfer_fold_keeps_so_count_stride_in_so_mode(self):
        # Transfer (enable_spin_orbital=True call) preserves the SO-count
        # stride via _reshape_orbit_spin; identity subshape (1,1,1) is an
        # invariance point that would break under an accidental swap.
        solver = self._make_stub(subshape=(1, 1, 1), cellshape=(4, 4, 1),
                                  norb_orig=2, enable_spin_orbital=True)
        ham = {((0, 0, 0), (1, 1)): 1.0}  # spin-down (SOI index 1)
        out = solver._reshape_interaction(ham, enable_spin_orbital=True)
        self.assertIn(((0, 0, 0), (1, 1)), out)

    def test_non_so_mode_stride_unchanged(self):
        # Non-SO mode with SubShape=(2,1,1) and norb_orig=2 keeps the original
        # stride (=2) since norb_phys_orig == norb_orig; a physical-orbital=1
        # entry folded to sub cell (1,0,0) must land at 1 + 2*1 = 3.
        solver = self._make_stub(subshape=(2, 1, 1), cellshape=(2, 1, 1),
                                  norb_orig=2, enable_spin_orbital=False)
        ham = {((0, 0, 0), (1, 1)): 1.0}
        out = solver._reshape_interaction(ham, enable_spin_orbital=False)
        folded_indices = sorted({i for (_ir, ov) in out.keys() for i in ov})
        self.assertEqual(folded_indices, [1, 3])


if __name__ == "__main__":
    unittest.main()
