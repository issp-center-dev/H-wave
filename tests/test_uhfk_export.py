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
    def test_save_greenone_emits_packed_spin_orbital_indices(self):
        # In spin-orbital mode the Green function has ns=1 and the orbital
        # slot carries the packed 2*a_phys + s index. _save_greenone must
        # index into it as gr[idx, 0, 2*a_phys+s, 0, 2*b_phys+t] and emit a
        # row for every (i, s, j, t) requested in onebodyg_uhf (including the
        # spin-off-diagonal s != t combinations).
        #
        # norb_phys = 2 (so packed axis size = 4) so that some rows have
        # a_phys = 1 and b_phys = 1: this genuinely exercises the 2*a
        # coefficient in the packing formula. With norb_phys = 1 every
        # row would have a_phys = b_phys = 0 and the mock would fail to
        # distinguish 2*a+s from 2*a-s (both collapse to +/-s under
        # numpy's negative-index wrap on an axis of length 2), from
        # a_phys+s, or even from just s.
        solver = UHFk.__new__(UHFk)
        solver.enable_spin_orbital = True
        solver.has_sublattice = False
        solver.ns = 1
        solver.norb = 4       # 4 packed spin-orbitals (norb_phys=2 x 2 spins)
        solver.nd = 4
        solver.cellshape = (1, 1, 1)
        solver.cellvol = 1
        solver.boundary_periodic = True
        solver.boundary_theta = (0.0, 0.0, 0.0)
        # Green[idx=0, 0, packed_a, 0, packed_b] where packed = 2*a_phys+s.
        # Fill all 4 x 4 = 16 slots with distinct sentinels so any
        # misindexing (e.g. 2*a - s in place of 2*a + s, or a + s in
        # place of 2*a + s) yields a value that mismatches the row's
        # expected sentinel.
        green = np.zeros((1, 1, 4, 1, 4), dtype=complex)
        for pa in range(4):
            for pb in range(4):
                green[0, 0, pa, 0, pb] = 0.01 * (pa * 4 + pb + 1) + 0.0j
        solver.Green = green

        # site 0 -> intra-cell orbital 0, site 1 -> intra-cell orbital 1.
        # Both sit in cell (0, 0, 0), so kx = ky = kz = 0 and idx = 0
        # for every (i, j) pair. This keeps the packing formula the only
        # moving part of the test.
        green_info = {
            "onebodyg_uhf": [
                (i, s, j, t)
                for i in range(2) for s in range(2)
                for j in range(2) for t in range(2)
            ],
            "geometry_uhf": {"site2vec": {0: (0, 0, 0, 0),
                                          1: (0, 0, 0, 1)}},
        }

        fd, path = tempfile.mkstemp(suffix=".dat")
        os.close(fd)
        try:
            solver._save_greenone(path, green_info)
            with open(path) as fp:
                lines = [line.split() for line in fp if line.strip()]
        finally:
            os.remove(path)

        # Each row: i s j t re im. 2 sites * 2 spins * 2 sites * 2 spins.
        self.assertEqual(len(lines), 16)
        parsed = {
            (int(row[0]), int(row[1]), int(row[2]), int(row[3])): float(row[4])
            for row in lines
        }
        # Under our site2vec, site i maps to a_phys = i and site j maps
        # to b_phys = j, so row (i, s, j, t) must read
        # gr[0, 0, 2*i + s, 0, 2*j + t] whose sentinel value is
        # 0.01 * ((2*i + s) * 4 + (2*j + t) + 1).
        for i in range(2):
            for s in range(2):
                for j in range(2):
                    for t in range(2):
                        expected = 0.01 * ((2 * i + s) * 4 + (2 * j + t) + 1)
                        self.assertAlmostEqual(
                            parsed[(i, s, j, t)], expected, places=10,
                            msg=(
                                "row (i={}, s={}, j={}, t={}) misindexed"
                            ).format(i, s, j, t),
                        )


if __name__ == "__main__":
    unittest.main()
