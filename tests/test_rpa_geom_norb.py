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


if __name__ == "__main__":
    unittest.main()
