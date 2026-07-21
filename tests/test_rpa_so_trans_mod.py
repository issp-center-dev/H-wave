#!/usr/bin/env python3
"""Spin-orbital (SO) multi-orbital ``trans_mod`` index-convention gate for RPA.

UHFk writes the SO ``trans_mod`` .npz with the orbital axis in INTERLEAVED order
(index = 2*orb + spin), the same convention as the SO Transfer file. RPA's
``_read_trans_mod`` / ``_reshape_green`` / ``_make_ham_trans`` consume H0 in
SPIN-BLOCK order (index = spin*norb_phys + orb). For >1 physical orbital these
differ, so the loaded ``trans_mod`` must be remapped interleaved->spin-block on
both orbital axes before reshape/fold, mirroring ``_make_ham_trans``'s remap of
the Transfer file.

Correctness gate (strongest feasible): build the interleaved real-space H0
directly from the SO Transfer fixture ``transfer_so_2orb.dat`` entries, feed it
as ``trans_mod`` to an SO RPA run, and require the resulting chi0q/chiq to match
the run that uses ``transfer_so_2orb.dat`` directly (which goes through the file
remap) to machine precision. Also a folded-vs-unfolded q=0 invariance check.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


INPUT_DIR = "tests/rpa/input"


def _build_interleaved_trans_mod_from_so_transfer(transfer_file, cellvol, nd0):
    """Build a real-space (cellvol, nd0, nd0) interleaved H0 from an SO Transfer.

    The SO Transfer fixture lists rows ``Rx Ry Rz i j Re Im`` with 1-based
    interleaved orbital indices (i = 2*orb + spin + 1). We place each entry at
    the matching real-space cell (R folded into [0, cellvol)) and 0-based
    interleaved index, so the returned array IS the transfer's H0 in interleaved
    layout -- exactly what UHFk would emit as trans_mod in SO mode.
    """
    tab = np.zeros((cellvol, nd0, nd0), dtype=np.complex128)
    with open(transfer_file) as f:
        lines = f.readlines()
    # header: comment, nlines/ncol?, ... mirror wan90 transfer header layout
    # Format used by these fixtures: line0 comment, line1 = number of orbitals,
    # line2 = number of cells(?), line3 = '1 1', then data rows.
    data_lines = []
    for ln in lines[4:]:
        s = ln.split()
        if len(s) >= 7:
            data_lines.append(s)
    for s in data_lines:
        rx = int(s[0])
        i = int(s[3]) - 1  # 0-based interleaved
        j = int(s[4]) - 1
        re = float(s[5])
        im = float(s[6])
        tab[rx % cellvol, i, j] += re + 1j * im
    return tab


def _run(info_mode_extra, inter, output_dir, trans_mod_file=None, subshape=(1, 1, 1)):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [8, 1, 1],
            "SubShape": list(subshape),
            "Nmat": 32,
        },
        "calc_scheme": "general",
    }
    info_mode.update(info_mode_extra)
    interaction = {"path_to_input": INPUT_DIR}
    interaction.update(inter)
    input_block = {
        "path_to_input": INPUT_DIR,
        "interaction": interaction,
    }
    if trans_mod_file is not None:
        input_block["trans_mod"] = trans_mod_file
    info_file = {
        "input": input_block,
        "output": {"path_to_output": output_dir},
    }
    os.makedirs(output_dir, exist_ok=True)
    read_io = read_input_k.QLMSkInput(info_file["input"])
    ham = read_io.get_param("ham")
    solver = solver_rpa.RPA(ham, {}, info_mode)
    green = read_io.get_param("green")
    info_in = solver.read_init(info_file["input"])
    for k, v in info_in.items():
        green[k] = v
    solver.solve(green, output_dir)
    return solver, green


class TestRPASOTransModMultiOrbital(unittest.TestCase):
    TM_NAME = "so_trans_mod_fixture.npz"

    def setUp(self):
        self._tmpdir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self._tmpdir, "output")
        os.makedirs(self.output_dir, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def _write_fixture(self):
        tab = _build_interleaved_trans_mod_from_so_transfer(
            os.path.join(INPUT_DIR, "transfer_so_2orb.dat"),
            cellvol=8, nd0=4)
        path = os.path.join(self._tmpdir, self.TM_NAME)
        np.savez(path, trans_mod=tab)
        return path

    def test_so_trans_mod_matches_so_transfer_file(self):
        """SO trans_mod (interleaved) must reproduce the SO transfer-file run."""
        tm_path = self._write_fixture()
        inter = {"Geometry": "geom_so_2orb.dat",
                 "Transfer": "transfer_so_2orb.dat",
                 "CoulombIntra": "coulombintra_2orb.dat"}
        # reference: SO run via the transfer-file remap (no trans_mod)
        s_ref, g_ref = _run({"enable_spin_orbital": True}, inter, self.output_dir)
        # under test: SO run that REPLACES H0 with the interleaved trans_mod
        s_tm, g_tm = _run({"enable_spin_orbital": True}, inter, self.output_dir,
                          trans_mod_file=tm_path)
        self.assertEqual(s_ref.spin_mode, "spin-free")
        self.assertEqual(s_tm.spin_mode, "spin-free")
        self.assertTrue(
            np.allclose(g_ref["chi0q"], g_tm["chi0q"]),
            "SO trans_mod chi0q must match the SO transfer-file chi0q",
        )
        self.assertTrue(
            np.allclose(g_ref["chiq"], g_tm["chiq"]),
            "SO trans_mod chiq must match the SO transfer-file chiq",
        )

    def test_so_trans_mod_matches_nonso(self):
        """SO trans_mod must also match the equivalent non-SO (ns=2) run."""
        tm_path = self._write_fixture()
        s_nonso, g_nonso = _run(
            {"enable_spin_orbital": False},
            {"Geometry": "geom_2orb.dat", "Transfer": "transfer_nonso_2orb.dat"},
            self.output_dir,
        )
        s_tm, g_tm = _run(
            {"enable_spin_orbital": True},
            {"Geometry": "geom_so_2orb.dat", "Transfer": "transfer_so_2orb.dat"},
            self.output_dir,
            trans_mod_file=tm_path,
        )
        self.assertEqual(s_nonso.spin_mode, "spin-free")
        self.assertEqual(s_tm.spin_mode, "spin-free")
        self.assertTrue(
            np.allclose(g_nonso["chi0q"], g_tm["chi0q"]),
            "SO trans_mod chi0q must match the equivalent non-SO chi0q",
        )


class TestRPASOTransModSublatticeFold(unittest.TestCase):
    TM_NAME = "so_trans_mod_fold_fixture.npz"

    def setUp(self):
        self._tmpdir = tempfile.mkdtemp()
        self.output_dir = os.path.join(self._tmpdir, "output")
        os.makedirs(self.output_dir, exist_ok=True)

    def tearDown(self):
        shutil.rmtree(self._tmpdir, ignore_errors=True)

    def _write_fixture(self):
        tab = _build_interleaved_trans_mod_from_so_transfer(
            os.path.join(INPUT_DIR, "transfer_so_2orb.dat"),
            cellvol=8, nd0=4)
        path = os.path.join(self._tmpdir, self.TM_NAME)
        np.savez(path, trans_mod=tab)
        return path

    @staticmethod
    def _uniform_q0_per_site(chiq, n_sites):
        no = chiq.shape[2]
        s = 0j
        for a in range(no):
            for b in range(no):
                s += chiq[:, 0, a, a, b, b].sum()
        return s / n_sites

    def test_folded_chiq_matches_unfolded(self):
        tm_path = self._write_fixture()
        inter = {"Geometry": "geom_so_2orb.dat",
                 "Transfer": "transfer_so_2orb.dat",
                 "CoulombIntra": "coulombintra_2orb.dat"}
        _su, g_unfold = _run({"enable_spin_orbital": True}, inter, self.output_dir,
                             trans_mod_file=tm_path, subshape=(1, 1, 1))
        _sf, g_fold = _run({"enable_spin_orbital": True}, inter, self.output_dir,
                           trans_mod_file=tm_path, subshape=(2, 1, 1))
        u = self._uniform_q0_per_site(g_unfold["chiq"], 1)
        f = self._uniform_q0_per_site(g_fold["chiq"], 2)
        self.assertAlmostEqual(u.real, f.real, places=10)
        self.assertAlmostEqual(u.imag, f.imag, places=10)


if __name__ == "__main__":
    unittest.main()
