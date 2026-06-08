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


if __name__ == "__main__":
    unittest.main()
