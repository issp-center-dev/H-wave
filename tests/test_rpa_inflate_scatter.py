"""The spin inflation of chi0q in RPA.solve (issue #309).

The spin-free and spin-diagonal branches used to build the spin-orbital chi0q
with einsums against a 0/1 spin tensor; they now scatter chi0q into a zeroed
buffer, which does not depend on the input's memory layout. This is a pure
rearrangement, so the inflated chi0q handed to ``_solve_rpa`` must equal the
einsum expressions exactly (values, dtype, shape). The one representational
difference is the sign of the structural zeros: the einsum forms ``c * 0``
(-0.0 where ``c`` is negative), the scatter leaves +0.0, as flex.py's scatter
always has; the two compare equal and give the same sums. The reference
einsums are pasted here; the inflated array is captured through the public
``solve`` path.
"""
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as rpa_mod

_INPUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rpa", "input")


def _reference(chi0q_orig, scheme, spin_mode, norb):
    """The einsum inflation as it was written before #309."""
    ns = 2
    nd = ns * norb
    if spin_mode == "spin-free":
        nfreq, nvol = chi0q_orig.shape[:2]
        if scheme == "reduced":
            return np.einsum('lkab,st->lksatb', chi0q_orig.reshape(nfreq, nvol, norb, norb),
                             np.identity(ns)).reshape(nfreq, nvol, nd, nd)
        spin_tensor = np.zeros((2, 2, 2, 2), dtype=np.int32)
        spin_tensor[0, 0, 0, 0] = 1
        spin_tensor[1, 1, 1, 1] = 1
        return np.einsum('lkabcd,stuv->lksatbucvd',
                         chi0q_orig.reshape((nfreq, nvol) + (norb,) * 4),
                         spin_tensor).reshape(nfreq, nvol, nd, nd, nd, nd)
    nfreq, nvol = chi0q_orig.shape[1:3]
    if scheme == "reduced":
        return np.einsum('glkab,gh->lkgahb', chi0q_orig,
                         np.identity(2)).reshape(nfreq, nvol, nd, nd)
    spin_tensor = np.zeros((2, 2, 2, 2), dtype=np.int32)
    spin_tensor[0, 0, 0, 0] = 1
    spin_tensor[1, 1, 1, 1] = 1
    return np.einsum('glkabcd,gtuv->lkgatbucvd', chi0q_orig,
                     spin_tensor).reshape(nfreq, nvol, nd, nd, nd, nd)


class TestInflationEqualsEinsum(unittest.TestCase):

    def setUp(self):
        self.out = tempfile.mkdtemp(prefix="rpa_inflate_")

    def tearDown(self):
        shutil.rmtree(self.out, ignore_errors=True)

    def _run(self, scheme, two_orb, extern):
        files = {"path_to_input": _INPUT,
                 "Geometry": "geom_2orb.dat" if two_orb else "geom.dat",
                 "Transfer": "transfer_nonso_2orb.dat" if two_orb else "transfer.dat",
                 "CoulombIntra": "coulombintra_2orb.dat" if two_orb else "coulombintra.dat"}
        param = {"T": 0.5, "mu": 0.0, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
                 "Nmat": 16}
        if extern:
            files["Extern"] = "extern.dat"
            param["coeff_extern"] = 0.3
        io = read_input_k.QLMSkInput({"path_to_input": _INPUT, "interaction": files})
        solver = rpa_mod.RPA(io.get_param("ham"), {},
                             {"mode": "RPA", "param": param, "calc_scheme": scheme})
        green_info = io.get_param("green")
        seen = []
        real = solver._solve_rpa

        def spy(chi0q, ham):
            seen.append((np.array(chi0q, copy=True), np.array(ham, copy=True)))
            return real(chi0q, ham)

        with mock.patch.object(solver, "_solve_rpa", side_effect=spy):
            solver.solve(green_info, self.out)
        return solver, green_info, seen[0], real

    def test_spin_free_and_spin_diag_branches(self):
        cases = [(scheme, two_orb, extern)
                 for scheme in ("reduced", "general")
                 for two_orb in (False, True)
                 for extern in (False, True)]
        for scheme, two_orb, extern in cases:
            with self.subTest(scheme=scheme, two_orb=two_orb, extern=extern):
                solver, green_info, (got, ham), solve_rpa = self._run(scheme, two_orb, extern)
                self.assertEqual(solver.spin_mode, "spin-diag" if extern else "spin-free")
                ref = _reference(np.asarray(green_info["chi0q"]), scheme,
                                 solver.spin_mode, solver.norb)
                self.assertEqual(got.shape, ref.shape)
                self.assertEqual(got.dtype, ref.dtype)
                self.assertTrue(np.array_equal(got, ref))
                # nonzero entries carry the same bits, signs included
                nz = ref != 0
                self.assertTrue(np.array_equal(got[nz].view(np.float64),
                                               ref[nz].view(np.float64)))
                # the inflation is not vacuous: off-spin blocks are zero, the
                # equal-spin ones carry chi0q
                self.assertGreater(np.abs(got).max(), 0.0)
                # and the RPA solve on either gives the same chiq values
                # (a structural zero may differ in sign only)
                chi_got, chi_ref = solve_rpa(got, ham), solve_rpa(ref, ham)
                self.assertTrue(np.array_equal(chi_got, chi_ref))


if __name__ == "__main__":
    unittest.main()
