#!/usr/bin/env python3
"""Regression tests for sublattice-folded Green deflate/reshape in UHFk.

The k-space Green function ``self.Green`` produced by ``_green()`` stores the
density matrix in the orbital-transposed slot convention relative to the
Hamiltonian. As a consequence the unfolded ``green`` output must place the
cross-sublattice offset on the FIRST orbital slot (``aa = a + norb_orig*ir``),
opposite to the Hamiltonian/transfer deflate (which uses the SECOND slot). If
both use the same convention, the unfolded ``green`` of a complex
(cross-sublattice) Hamiltonian comes out as the Hermitian conjugate of the
reference ``SubShape=[1,1,1]`` result.

Minimal model: 1D complex Peierls chain.
  CellShape=[4,1,1], norb=1, Ncond=4, 2Sz=0, U=0
  T(R=+1) = -exp(+i*0.4), T(R=-1) = -exp(-i*0.4)
"""

import os
import tempfile
import shutil
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa
from hwave.solver.uhfk import UHFk


GEOM = """  1.000000000000   0.000000000000   0.000000000000
  0.000000000000   1.000000000000   0.000000000000
  0.000000000000   0.000000000000   1.000000000000
1
    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00
"""


def _transfer_text():
    tp = -np.exp(1j * 0.4)   # R = +1
    tm = -np.exp(-1j * 0.4)  # R = -1
    return (
        "Transfer Peierls 1D chain\n"
        "1\n"
        "2\n"
        " 1 1\n"
        "  1    0    0    1    1  {:.12f}  {:.12f}\n"
        " -1    0    0    1    1  {:.12f}  {:.12f}\n"
    ).format(tp.real, tp.imag, tm.real, tm.imag)


def _build_solver(case_dir, subshape):
    info_inputfile = {
        "path_to_input": "",
        "interaction": {
            "path_to_input": case_dir,
            "Geometry": "geom.dat",
            "Transfer": "transfer.dat",
        },
    }
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")
    info_log = {"print_level": 0, "print_step": 1}
    info_mode = {
        "mode": "UHFk",
        "param": {
            "Ncond": 4, "2Sz": 0, "IterationMax": 2000, "EPS": 12,
            "Mix": 0.5, "RndSeed": 123456789, "T": 0.0,
            "CellShape": [4, 1, 1], "SubShape": list(subshape),
        },
    }
    solver = UHFk(ham_info, info_log, info_mode)
    green_info = read_io.get_param("green")
    solver.solve(green_info, os.path.join(case_dir, "output"))
    return solver


# Spin-orbital fixtures: same physical Peierls chain, one physical orbital, with
# spin folded into the orbital index (interleaved 1=up, 2=down). The hopping is
# spin-diagonal, so the physics matches the non-SO chain.
GEOM_SO = ("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n2\n"
           "  0.0 0.0 0.0\n  0.0 0.0 0.0\n")


def _transfer_text_so():
    tp = -np.exp(1j * 0.4)
    tm = -np.exp(-1j * 0.4)

    def row(r, o, v):
        return " {:2d}    0    0    {}    {}  {:.12f}  {:.12f}\n".format(r, o, o, v.real, v.imag)

    return ("Transfer SO Peierls\n2\n2\n 1 1\n"
            + row(1, 1, tp) + row(1, 2, tp) + row(-1, 1, tm) + row(-1, 2, tm))


def _build_solver_so(case_dir, subshape):
    info_inputfile = {
        "path_to_input": "",
        "interaction": {
            "path_to_input": case_dir,
            "Geometry": "geom_so.dat",
            "Transfer": "transfer_so.dat",
        },
    }
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")
    info_log = {"print_level": 0, "print_step": 1}
    info_mode = {
        "mode": "UHFk",
        "enable_spin_orbital": True,
        "param": {
            "Ncond": 4, "IterationMax": 2000, "EPS": 12,
            "Mix": 0.5, "RndSeed": 123456789, "T": 0.0,
            "CellShape": [4, 1, 1], "SubShape": list(subshape),
        },
    }
    solver = UHFk(ham_info, info_log, info_mode)
    green_info = read_io.get_param("green")
    solver.solve(green_info, os.path.join(case_dir, "output"))
    return solver


class TestGreenDeflateSubShape(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp(prefix="hwave_peierls_")
        with open(os.path.join(cls.tmp, "geom.dat"), "w") as f:
            f.write(GEOM)
        with open(os.path.join(cls.tmp, "transfer.dat"), "w") as f:
            f.write(_transfer_text())
        with open(os.path.join(cls.tmp, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra Peierls\n1\n1\n 1\n"
                    "   0    0    0    1    1   2.000000000000   0.000000000000\n")
        with open(os.path.join(cls.tmp, "geom_so.dat"), "w") as f:
            f.write(GEOM_SO)
        with open(os.path.join(cls.tmp, "transfer_so.dat"), "w") as f:
            f.write(_transfer_text_so())
        cls.s111 = _build_solver(cls.tmp, (1, 1, 1))
        cls.s211 = _build_solver(cls.tmp, (2, 1, 1))

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_green_subshape_invariance(self):
        # unfolded green must match the SubShape=[1,1,1] reference
        g_ref = self.s111.Green                       # no sublattice -> already unfolded
        g_def = self.s211._deflate_green(self.s211.Green)
        self.assertLess(np.max(np.abs(g_def - g_ref)), 1.0e-9)

    def test_reshape_green_round_trip(self):
        # folding the unfolded reference must reproduce the internal folded Green
        g_ref = self.s111.Green
        g_folded = self.s211._reshape_green(g_ref)
        self.assertLess(np.max(np.abs(g_folded - self.s211.Green)), 1.0e-9)

    def test_trans_mod_unchanged(self):
        # the Hamiltonian/transfer deflate must remain SubShape-invariant
        s111, s211 = self.s111, self.s211

        def trans_r(solver):
            nx, ny, nz = solver.shape
            nvol = solver.nvol
            nd = solver.nd
            tab_k = solver.ham
            return np.fft.fftn(tab_k.reshape(nx, ny, nz, nd, nd),
                               axes=(0, 1, 2)).reshape(nvol, nd, nd) / nvol

        ns, norb_orig = s111.ns, s111.norb_orig
        lvol = s111.cellvol
        tr_ref = trans_r(s111).reshape(lvol, ns, norb_orig, ns, norb_orig)
        norb, nvol = s211.norb, s211.nvol
        tab_r = trans_r(s211).reshape(nvol, ns, norb, ns, norb)
        tr_def = s211._deflate_green(tab_r, hamiltonian=True)
        # full (spin, orbital, spin, orbital) tensor must match the reference
        self.assertEqual(tr_def.shape, tr_ref.shape)
        self.assertLess(np.max(np.abs(tr_def - tr_ref)), 1.0e-9)

    def test_green_subshape_invariance_spin_orbital(self):
        # The deflate/reshape path is shared by spin-orbital mode (ns=1, spin
        # folded into the orbital index). The same fix must keep the SO green
        # SubShape-invariant for complex cross-sublattice hopping.
        so111 = _build_solver_so(self.tmp, (1, 1, 1))
        so211 = _build_solver_so(self.tmp, (2, 1, 1))
        g_ref = so111.Green                            # no sublattice -> unfolded
        g_def = so211._deflate_green(so211.Green)
        self.assertEqual(g_def.shape, g_ref.shape)
        self.assertLess(np.max(np.abs(g_def - g_ref)), 1.0e-9)

    def test_rpa_green_init_fold_matches_uhfk(self):
        # RPA folds UHFk's saved (unfolded) green via RPA._reshape_green. It must
        # reconstruct UHFk's correct internal folded Green, i.e. RPA's green fold
        # must stay the inverse of UHFk._deflate_green's Green convention.
        info_mode = {
            "mode": "RPA",
            "param": {
                "T": 2.0, "filling": 0.5,
                "CellShape": [4, 1, 1], "SubShape": [2, 1, 1], "Nmat": 16,
            },
            "enable_spin_orbital": False,
            "calc_scheme": "general",
        }
        inter = {
            "path_to_input": self.tmp,
            "Geometry": "geom.dat",
            "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat",
        }
        read_io = read_input_k.QLMSkInput({"path_to_input": self.tmp, "interaction": inter})
        ham = read_io.get_param("ham")
        rpa = solver_rpa.RPA(ham, {}, info_mode)

        unfolded = self.s111.Green                    # SubShape-invariant reference
        folded_rpa = rpa._reshape_green(unfolded)
        self.assertEqual(folded_rpa.shape, self.s211.Green.shape)
        self.assertLess(np.max(np.abs(folded_rpa - self.s211.Green)), 1.0e-9)

    def test_rpa_read_green_file_path(self):
        # Exercise the real RPA._read_green plumbing (npz load, 5D collapse,
        # shape validation, sublattice fold) on a UHFk-saved green file.
        info_mode = {
            "mode": "RPA",
            "param": {
                "T": 2.0, "filling": 0.5,
                "CellShape": [4, 1, 1], "SubShape": [2, 1, 1], "Nmat": 16,
            },
            "enable_spin_orbital": False,
            "calc_scheme": "general",
        }
        inter = {
            "path_to_input": self.tmp,
            "Geometry": "geom.dat",
            "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat",
        }
        read_io = read_input_k.QLMSkInput({"path_to_input": self.tmp, "interaction": inter})
        ham = read_io.get_param("ham")
        rpa = solver_rpa.RPA(ham, {}, info_mode)

        # write a green file in UHFk's _save_green format (unfolded "green" key)
        npz_path = os.path.join(self.tmp, "green_init.npz")
        green_unfolded = self.s211._deflate_green(self.s211.Green)
        np.savez(npz_path, green=green_unfolded, green_sublattice=self.s211.Green, momentum_convention="e_plus_ikR")

        nvol, nd = rpa.lattice.nvol, rpa.nd
        g = rpa._read_green(npz_path)
        self.assertEqual(g.shape, (nvol, nd, nd))
        self.assertLess(
            np.max(np.abs(g - self.s211.Green.reshape(nvol, nd, nd))), 1.0e-9)

    def _rpa_for_read(self):
        info_mode = {
            "mode": "RPA",
            "param": {
                "T": 2.0, "filling": 0.5,
                "CellShape": [4, 1, 1], "SubShape": [2, 1, 1], "Nmat": 16,
            },
            "enable_spin_orbital": False,
            "calc_scheme": "general",
        }
        inter = {
            "path_to_input": self.tmp,
            "Geometry": "geom.dat",
            "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat",
        }
        read_io = read_input_k.QLMSkInput({"path_to_input": self.tmp, "interaction": inter})
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_rpa_read_old_sublattice_green_uses_green_sublattice(self):
        # Issue #36: a pre-#35 UHFk sublattice green.npz stored the "green" key
        # deflated with the OLD Hamiltonian/transfer convention, while
        # "green_sublattice" held the already-folded (still correct) Green.
        # RPA must NOT fold the ambiguous old "green" key (which would be
        # silently wrong); it must recover the correct Green from
        # "green_sublattice".
        rpa = self._rpa_for_read()
        nvol, nd = rpa.lattice.nvol, rpa.nd

        npz_path = os.path.join(self.tmp, "green_old.npz")
        old_green = self.s211._deflate_green(self.s211.Green, hamiltonian=True)
        np.savez(npz_path, green=old_green, green_sublattice=self.s211.Green, momentum_convention="e_plus_ikR")

        g = rpa._read_green(npz_path)
        self.assertEqual(g.shape, (nvol, nd, nd))
        self.assertLess(
            np.max(np.abs(g - self.s211.Green.reshape(nvol, nd, nd))), 1.0e-9)

    def test_rpa_read_ambiguous_sublattice_green_raises(self):
        # A sublattice green.npz with neither the new convention metadata nor a
        # "green_sublattice" array is ambiguous: the "green" key could be folded
        # with either convention. RPA must fail loudly rather than silently
        # produce a wrong result.
        rpa = self._rpa_for_read()

        npz_path = os.path.join(self.tmp, "green_ambiguous.npz")
        old_green = self.s211._deflate_green(self.s211.Green, hamiltonian=True)
        np.savez(npz_path, green=old_green, momentum_convention="e_plus_ikR")

        with self.assertRaises(SystemExit):
            rpa._read_green(npz_path)

    def test_uhfk_save_green_writes_convention_metadata(self):
        # New UHFk sublattice green.npz must be self-describing so RPA can fold
        # the "green" key with the correct convention.
        npz_path = os.path.join(self.tmp, "green_meta.npz")
        self.s211._save_green(npz_path)
        data = np.load(npz_path, allow_pickle=False)
        self.assertIn("green_convention", data.files)
        self.assertEqual(str(data["green_convention"]), "green_slot_first")

    def _rpa_so_for_read(self):
        info_mode = {
            "mode": "RPA",
            "param": {
                "T": 2.0, "filling": 0.5,
                "CellShape": [4, 1, 1], "SubShape": [2, 1, 1], "Nmat": 16,
            },
            "enable_spin_orbital": True,
            "calc_scheme": "general",
        }
        inter = {
            "path_to_input": self.tmp,
            "Geometry": "geom_so.dat",
            "Transfer": "transfer_so.dat",
        }
        read_io = read_input_k.QLMSkInput({"path_to_input": self.tmp, "interaction": inter})
        ham = read_io.get_param("ham")
        return solver_rpa.RPA(ham, {}, info_mode)

    def test_rpa_so_sublattice_green_sublattice_without_tag_raises(self):
        # In spin-orbital mode green_sublattice is stored INTERLEAVED while RPA
        # works in spin-block order, so the green_sublattice shortcut is not safe
        # in general. An SO sublattice file carrying green_sublattice but no
        # convention tag must fail-fast (not be silently misread) -- issue #36.
        so211 = _build_solver_so(self.tmp, (2, 1, 1))
        rpa = self._rpa_so_for_read()

        npz_path = os.path.join(self.tmp, "green_so_notag.npz")
        old_green = so211._deflate_green(so211.Green, hamiltonian=True)
        np.savez(npz_path, green=old_green, green_sublattice=so211.Green, momentum_convention="e_plus_ikR")

        with self.assertRaises(SystemExit):
            rpa._read_green(npz_path)

    def test_rpa_so_sublattice_green_with_tag_reads(self):
        # With the convention tag present, an SO sublattice green file is folded
        # through the tag-gated "green" path (which applies the interleaved->
        # spin-block remap). The SO Peierls chain (norb_phys=1, spin-diagonal)
        # is physically identical to the non-SO chain, so the loaded green_init
        # must reproduce the non-SO reference self.s211.Green to machine
        # precision -- a numerical check, not just a shape check.
        so211 = _build_solver_so(self.tmp, (2, 1, 1))
        rpa = self._rpa_so_for_read()
        nvol, nd = rpa.lattice.nvol, rpa.nd

        npz_path = os.path.join(self.tmp, "green_so_tag.npz")
        new_green = so211._deflate_green(so211.Green)
        np.savez(npz_path, green=new_green,
                 green_convention=np.array("green_slot_first"), momentum_convention="e_plus_ikR")

        g = rpa._read_green(npz_path)
        self.assertEqual(g.shape, (nvol, nd, nd))
        self.assertLess(
            np.max(np.abs(g - self.s211.Green.reshape(nvol, nd, nd))), 1.0e-9)

    def test_rpa_read_new_sublattice_green_with_metadata(self):
        # With the new-convention metadata present, folding the "green" key is
        # unambiguous and must reproduce the correct internal folded Green even
        # when "green_sublattice" is absent.
        rpa = self._rpa_for_read()
        nvol, nd = rpa.lattice.nvol, rpa.nd

        npz_path = os.path.join(self.tmp, "green_new_nogsub.npz")
        new_green = self.s211._deflate_green(self.s211.Green)
        np.savez(npz_path, green=new_green,
                 green_convention=np.array("green_slot_first"), momentum_convention="e_plus_ikR")

        g = rpa._read_green(npz_path)
        self.assertEqual(g.shape, (nvol, nd, nd))
        self.assertLess(
            np.max(np.abs(g - self.s211.Green.reshape(nvol, nd, nd))), 1.0e-9)


if __name__ == "__main__":
    unittest.main()
