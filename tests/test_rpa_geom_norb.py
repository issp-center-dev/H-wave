#!/usr/bin/env python3
"""geom-norb SO-count convention: derivation, validation, fold stride, guards."""

import os
import unittest

import numpy as np

import hwave.solver.rpa as solver_rpa
from hwave.solver.rpa import _so_physical_norb


class TestSoPhysicalNorb(unittest.TestCase):
    def test_non_so_passthrough(self):
        self.assertEqual(_so_physical_norb(3, False), 3)

    def test_so_halves_even(self):
        self.assertEqual(_so_physical_norb(4, True), 2)

    def test_so_odd_raises(self):
        with self.assertRaises(ValueError):
            _so_physical_norb(3, True)

    def test_so_check_norb_targets_prefold(self):
        # divide the post-fold value (8) but validate evenness on pre-fold (3)
        with self.assertRaises(ValueError):
            _so_physical_norb(8, True, check_norb=3)

    def test_so_check_norb_even_but_different_returns_geom_half(self):
        # post-fold geom_norb=8 with even pre-fold check_norb=4 -> halve geom_norb
        self.assertEqual(_so_physical_norb(8, True, check_norb=4), 4)


class TestSoTransferIndexRange(unittest.TestCase):
    def _construct_with_transfer(self, transfer_name):
        import hwave.qlmsio.read_input_k as read_input_k
        info_mode = {
            "mode": "RPA",
            "param": {"T": 2.0, "filling": 0.5,
                      "CellShape": [8, 1, 1], "SubShape": [1, 1, 1], "Nmat": 32},
            "enable_spin_orbital": True,
            "calc_scheme": "general",
        }
        info_file = {"input": {"path_to_input": "tests/rpa/input",
                               "interaction": {"path_to_input": "tests/rpa/input",
                                               "Geometry": "geom_so_2orb.dat",
                                               "Transfer": transfer_name}},
                     "output": {"path_to_output": "tests/rpa/output"}}
        os.makedirs(info_file["input"]["interaction"]["path_to_input"], exist_ok=True)
        read_io = read_input_k.QLMSkInput(info_file["input"])
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_out_of_range_so_index_raises(self):
        # geom_so_2orb.dat declares nd=4 (indices 0..3); an index of 4 is illegal.
        with self.assertRaises(ValueError):
            self._construct_with_transfer("transfer_so_oob.dat")


class _FakeLattice:
    def __init__(self, shape, subshape):
        self.shape = shape
        self.subshape = subshape


class TestReshapeInteractionStride(unittest.TestCase):
    """P4: non-Transfer (physical-indexed) interactions fold with the physical
    stride, while Transfer (spin-orbital-indexed) folds with the SO-count stride."""

    def _make_interaction(self, geom_norb_orig, enable_spin_orbital, subshape):
        obj = object.__new__(solver_rpa.Interaction)
        obj.enable_spin_orbital = enable_spin_orbital
        obj.lattice = _FakeLattice(shape=(2, 1, 1), subshape=subshape)
        obj.param_ham_orig = {"Geometry": {"norb": geom_norb_orig}}
        return obj

    def test_two_body_stride_is_physical_in_so_mode(self):
        # SO mode, geom norb (SO count) = 4 -> physical = 2. A two-body term
        # (enable_spin_orbital=False call) on physical orbital index 1, folded
        # over subshape (2,1,1), must land within [0, norb_phys*subvol) = [0,4),
        # not stride by the SO count (which would reach 5).
        obj = self._make_interaction(4, True, subshape=(2, 1, 1))
        ham = {((0, 0, 0), (1, 1)): 1.0}
        out = obj._reshape_interaction(ham, enable_spin_orbital=False)
        folded_indices = [i for (_ir, ov) in out.keys() for i in ov]
        self.assertTrue(all(0 <= i < 4 for i in folded_indices),
                        "physical-indexed fold must stride by norb_phys; got "
                        "{}".format(folded_indices))

    def test_transfer_stride_is_so_count_in_so_mode(self):
        # Transfer (enable_spin_orbital=True call) keeps the SO-count stride;
        # subshape (1,1,1) => identity fold, index 3 preserved (< SO count 4).
        obj = self._make_interaction(4, True, subshape=(1, 1, 1))
        ham = {((0, 0, 0), (3, 3)): 1.0}
        out = obj._reshape_interaction(ham, enable_spin_orbital=True)
        self.assertIn(((0, 0, 0), (3, 3)), out)


class TestTransModSublatticeSupported(unittest.TestCase):
    """trans_mod + sublattice folding is now supported (was a NotImplemented
    guard over a pre-existing 2023 mis-wire; see test_rpa_trans_mod.py for the
    folded-vs-unfolded chiq invariance gate)."""

    def _solver(self, subshape):
        import hwave.qlmsio.read_input_k as read_input_k
        info_mode = {
            "mode": "RPA",
            "param": {"T": 2.0, "filling": 0.5,
                      "CellShape": [8, 1, 1], "SubShape": list(subshape), "Nmat": 32},
            "enable_spin_orbital": False,
            "calc_scheme": "general",
        }
        info_file = {"input": {"path_to_input": "tests/rpa/input",
                               "interaction": {"path_to_input": "tests/rpa/input",
                                               "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat"}},
                     "output": {"path_to_output": "tests/rpa/output"}}
        read_io = read_input_k.QLMSkInput(info_file["input"])
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_trans_mod_with_sublattice_folds(self):
        import tempfile
        solver = self._solver(subshape=(2, 1, 1))
        # geom.dat is norb=1 (non-SO): nd0 = ns*norb_orig = 2, cellvol = 8.
        ns, norb_orig, cellvol = 2, 1, 8
        nd0 = ns * norb_orig
        rng = np.random.default_rng(0)
        a = rng.standard_normal((cellvol, nd0, nd0))
        tab = (a + a.transpose(0, 2, 1)) * 0.5  # symmetric (real Hermitian)
        with tempfile.TemporaryDirectory() as d:
            fn = os.path.join(d, "trans_mod.npz")
            np.savez(fn, trans_mod=tab)
            tab_k = solver._read_trans_mod(fn)  # must NOT raise
        # folded supercell: nvol = 4, nd = norb*subvol*ns = 1*2*2 = 4
        self.assertEqual(tab_k.shape, (solver.lattice.nvol, solver.nd, solver.nd))


if __name__ == "__main__":
    unittest.main()
