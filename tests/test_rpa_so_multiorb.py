#!/usr/bin/env python3
"""Multi-orbital spin-orbital (SO) index-convention equivalence for RPA.

The SO transfer file uses the interleaved convention (file index = 2*orb+spin+1,
matching UHFk and the docs). RPA must interpret it that way. A spin-independent,
orbital-dependent 2-orbital system encoded in SO form must therefore:

  * be detected as spin-free (not spin-diag), and
  * yield the same chi0q as the equivalent non-SO (ns=2) calculation.

For norb_phys=1 the interleaved and spin-block conventions coincide, which is why
all pre-existing SO tests pass; this 2-orbital case is what exposes the bug.
"""

import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


def _run(transfer, spin_orbital, geom):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [8, 1, 1],
            "SubShape": [1, 1, 1],
            "Nmat": 32,
        },
        "enable_spin_orbital": spin_orbital,
        "calc_scheme": "general",
    }
    info_file = {
        "input": {
            "path_to_input": "tests/rpa/input",
            "interaction": {
                "path_to_input": "tests/rpa/input",
                "Geometry": geom,
                "Transfer": transfer,
            },
        },
        "output": {"path_to_output": "tests/rpa/output"},
    }
    os.makedirs(info_file["output"]["path_to_output"], exist_ok=True)
    read_io = read_input_k.QLMSkInput(info_file["input"])
    ham = read_io.get_param("ham")
    solver = solver_rpa.RPA(ham, {}, info_mode)
    green = read_io.get_param("green")
    solver.solve(green, info_file["output"]["path_to_output"])
    return solver, green


class TestRPAMultiOrbitalSO(unittest.TestCase):
    def test_spin_independent_2orbital_so_is_spin_free(self):
        solver, _ = _run("transfer_so_2orb.dat", True, "geom_2orb.dat")
        # A spin-independent system must be detected as spin-free, not spin-diag.
        self.assertEqual(solver.spin_mode, "spin-free")

    def test_so_chi0q_matches_nonso(self):
        s_nonso, g_nonso = _run("transfer_nonso_2orb.dat", False, "geom_2orb.dat")
        s_so, g_so = _run("transfer_so_2orb.dat", True, "geom_2orb.dat")
        self.assertEqual(s_nonso.spin_mode, "spin-free")
        self.assertEqual(s_so.spin_mode, "spin-free")
        self.assertTrue(
            np.allclose(g_nonso["chi0q"], g_so["chi0q"]),
            "SO chi0q must match the equivalent non-SO chi0q",
        )


if __name__ == "__main__":
    unittest.main()
