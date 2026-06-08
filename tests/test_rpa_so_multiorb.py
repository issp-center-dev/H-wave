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


def _run(transfer, spin_orbital, geom, subshape=(1, 1, 1)):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [8, 1, 1],
            "SubShape": list(subshape),
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


class TestRPAMultiOrbitalSOSublatticeFold(unittest.TestCase):
    """Folding (SubShape vol > 1) with multi-orbital SO is now supported."""

    def test_subcell_folding_runs(self):
        solver, _ = _run("transfer_so_2orb.dat", True, "geom_so_2orb.dat",
                         subshape=(2, 1, 1))
        self.assertEqual(solver.spin_mode, "spin-free")

    def test_folded_chi0q_matches_unfolded(self):
        # Physical observable must not depend on the (artificial) sublattice fold.
        _su, g_unfold = _run("transfer_so_2orb.dat", True,
                             "geom_so_2orb.dat", subshape=(1, 1, 1))
        _sf, g_fold = _run("transfer_so_2orb.dat", True,
                           "geom_so_2orb.dat", subshape=(2, 1, 1))
        #
        # The raw |chi0q|.sum() is NOT fold-invariant: folding by SubShape
        # [2,1,1] moves a factor of 2 from the q-index into the orbital index
        # (unfold shape (Nm,8,2,2,2,2) -> fold shape (Nm,4,4,4,4,4)), so the
        # array layout and BZ grid differ. The genuinely fold-invariant physical
        # quantity is the UNIFORM (q=0) static susceptibility per physical site:
        # the q=0 component is the same physical point (uniform response) in both
        # Brillouin zones, and the uniform response is independent of the choice
        # of unit cell. We sum over all orbital pairs chi_{a,a,b,b}(q=0) (over
        # both spin-orbitals and, in the folded supercell, both sublattice
        # sites -- capturing the full uniform response of the cell) and divide
        # by the number of physical sites in the cell (1 unfolded, 2 folded).
        def uniform_q0_per_site(chi0q, n_sites):
            no = chi0q.shape[2]
            s = 0j
            for a in range(no):
                for b in range(no):
                    s += chi0q[:, 0, a, a, b, b].sum()
            return s / n_sites

        u = uniform_q0_per_site(g_unfold["chi0q"], 1)
        f = uniform_q0_per_site(g_fold["chi0q"], 2)
        self.assertAlmostEqual(u.real, f.real, places=10)
        self.assertAlmostEqual(u.imag, f.imag, places=10)


class TestRPAMultiOrbitalSO(unittest.TestCase):
    def test_spin_independent_2orbital_so_is_spin_free(self):
        solver, _ = _run("transfer_so_2orb.dat", True, "geom_so_2orb.dat")
        # A spin-independent system must be detected as spin-free, not spin-diag.
        self.assertEqual(solver.spin_mode, "spin-free")

    def test_so_chi0q_matches_nonso(self):
        s_nonso, g_nonso = _run("transfer_nonso_2orb.dat", False, "geom_2orb.dat")
        s_so, g_so = _run("transfer_so_2orb.dat", True, "geom_so_2orb.dat")
        self.assertEqual(s_nonso.spin_mode, "spin-free")
        self.assertEqual(s_so.spin_mode, "spin-free")
        self.assertTrue(
            np.allclose(g_nonso["chi0q"], g_so["chi0q"]),
            "SO chi0q must match the equivalent non-SO chi0q",
        )


if __name__ == "__main__":
    unittest.main()
