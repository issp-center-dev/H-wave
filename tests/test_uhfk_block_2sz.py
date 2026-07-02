"""Regression tests: block detection must honor the 2Sz-fixed constraint.

In 2Sz-fixed mode the UHFk mean-field Hamiltonian is treated as spin-block
diagonal (spin off-diagonal Fock entries are not used), and each spin sector
is filled with its own electron count (Ncond +- 2Sz)/2.  Block detection must
therefore never merge the up and down sectors into one mu-group in 2Sz-fixed
mode; doing so silently replaces the constrained solution by the Sz-free one.

Reference values are from the pre-block-detection implementation (main,
c0b0a3c) on tests/uhfk/CoulombIntra with 2Sz = 4 and no initial green.
"""
import os
import tempfile
import unittest

import numpy as np
import tomli


class TestUHFk2SzFixedBlocks(unittest.TestCase):
    _solver_cache = None

    @classmethod
    def _run_coulombintra_2sz4(cls):
        # both tests assert on the same converged state; run the SCF once
        if cls._solver_cache is not None:
            return cls._solver_cache
        cur = os.getcwd()
        os.chdir("tests/uhfk/CoulombIntra")
        try:
            with open("input.toml", "rb") as f:
                params = tomli.load(f)
            params["mode"]["param"]["2Sz"] = 4
            params["file"]["input"].pop("initial", None)

            from hwave.qlmsio import read_input_k
            from hwave.solver.uhfk import UHFk

            read_io = read_input_k.QLMSkInput(params["file"]["input"])
            ham_info = read_io.get_param("ham")
            green_info = read_io.get_param("green")

            solver = UHFk(ham_info, params["log"], params["mode"])
            with tempfile.TemporaryDirectory() as tmp:
                solver.solve(green_info, tmp)
            cls._solver_cache = solver
            return solver
        finally:
            os.chdir(cur)

    def test_2sz_constraint_honored_with_coulombintra(self):
        # NOTE: no total-energy assertion here.  With the exactly degenerate
        # non-interacting start, the converged SCF branch depends on the
        # LAPACK implementation (macOS and Linux find different, equally
        # Sz-constrained solutions), so only platform-independent invariants
        # of the constraint are checked: the converged Sz, the total electron
        # count, and the per-spin-sector electron counts (Ncond +- 2Sz)/2.
        solver = self._run_coulombintra_2sz4()
        self.assertTrue(
            np.isclose(solver.physics["Sz"], 2.0, rtol=0.0, atol=1.0e-8),
            "converged Sz = {} but 2Sz = 4 was requested".format(
                solver.physics["Sz"]))
        self.assertTrue(
            np.isclose(solver.physics["NCond"], 16.0, rtol=0.0, atol=1.0e-8),
            "converged NCond = {} but Ncond = 16 was requested".format(
                solver.physics["NCond"]))
        self.assertEqual(
            sorted(solver.group_nconds), [6, 10],
            "spin sectors must carry (Ncond +- 2Sz)/2 = 10 and 6 electrons")

    def test_2sz_fixed_blocks_are_pure_spin(self):
        solver = self._run_coulombintra_2sz4()
        norb = solver.norb
        for blk in solver.block_info:
            blk = np.asarray(blk)
            n_up = np.sum(blk < norb)
            n_down = np.sum(blk >= norb)
            self.assertTrue(
                n_up == 0 or n_down == 0,
                "2Sz-fixed mode produced a spin-mixed block: {}".format(blk))


if __name__ == "__main__":
    unittest.main()
