#!/usr/bin/env python3
"""Round-trip and SO-safety tests for UHFk output writers."""

import os
import tempfile
import unittest

import numpy as np

from hwave.solver.uhfk import UHFk
from hwave.qlmsio import wan90


class TestExportInteractionRoundTrip(unittest.TestCase):
    def test_exported_orbital_indices_round_trip(self):
        # _export_interaction must write 1-based orbital indices so that
        # read_w90 (which subtracts 1) recovers the original 0-based indices.
        solver = UHFk.__new__(UHFk)
        solver.norb = 2
        original = {
            ((0, 0, 0), (0, 1)): 2.0 + 0.0j,
            ((1, 0, 0), (1, 1)): -1.0 + 0.0j,
        }
        solver.param_ham = {"CoulombInter": original}

        fd, path = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        try:
            solver._export_interaction("CoulombInter", path)
            recovered = wan90.read_w90(path)
        finally:
            os.remove(path)

        for key, val in original.items():
            self.assertIn(key, recovered)
            self.assertAlmostEqual(recovered[key].real, val.real, places=10)


class TestSaveGreenoneSpinOrbital(unittest.TestCase):
    def test_save_greenone_rejected_in_spin_orbital_mode(self):
        # In spin-orbital mode the Green function has ns=1 (spin folded into the
        # orbital index), so the (i,s,j,t) one-body output is not representable.
        # It must be rejected with a clear error rather than crashing.
        solver = UHFk.__new__(UHFk)
        solver.enable_spin_orbital = True
        solver.has_sublattice = False
        solver.ns = 1
        solver.norb = 2
        solver.nd = 2
        solver.cellshape = (1, 1, 1)
        solver.cellvol = 1
        solver.Green = np.zeros((1, 1, 2, 1, 2), dtype=complex)

        green_info = {
            "onebodyg_uhf": [(0, 1, 0, 1)],  # spin index s=t=1
            "geometry_uhf": {"site2vec": {0: (0, 0, 0, 0)}},
        }

        fd, path = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        try:
            # Must not raise (no IndexError); returns None to signal "not written".
            result = solver._save_greenone(path, green_info)
        finally:
            os.remove(path)
        self.assertIsNone(result)


if __name__ == "__main__":
    unittest.main()
