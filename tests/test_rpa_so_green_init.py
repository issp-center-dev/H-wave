#!/usr/bin/env python3
"""Spin-orbital (SO) multi-orbital ``green_init`` index-convention gate for RPA.

UHFk writes the saved ``green`` .npz (consumed by RPA as ``green_init``) with the
orbital axis in INTERLEAVED order (index = 2*orb + spin), the same convention as
the SO Transfer file and the SO ``trans_mod``. RPA's ``_calc_trans_mod`` consumes
``green_init`` in SPIN-BLOCK order (index = spin*norb_phys + orb), using
``self.norb`` / ``self.nd`` and the spin-block transfer. For >1 physical orbital
these differ, so the loaded ``green_init`` must be remapped interleaved->spin-block
on BOTH (nd,nd) axes (and sublattice-folded) before use, mirroring
``_read_trans_mod``. For norb_phys=1 the remap is the identity.

Correctness gate (strongest feasible): build a spin-independent, translation-
invariant green in SPIN-BLOCK order, feed it to a non-SO (ns=2) run as
``green_init``; build the INTERLEAVED version of the SAME green, feed it to the
equivalent SO run as ``green_init``; require the resulting chi0q/chiq to match to
machine precision. Plus a folded-vs-unfolded q=0 invariance check for the SO path.
"""

import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


INPUT_DIR = "tests/rpa/input"
OUTPUT_DIR = "tests/rpa/output"

CELLVOL = 8
NORB_PHYS = 2


def _spin_block_green(cellvol, norb_phys, ns=2):
    """Build a spin-independent, translation-invariant green in SPIN-BLOCK order.

    Layout: index = spin*norb_phys + orb, shape (cellvol, ns*norb_phys, ns*norb_phys).
    The (nd,nd) block is Hermitian (real symmetric here) and identical for both
    spins; orbitals are made distinguishable (different diagonal weights and a
    nonzero 0-1 off-diagonal) so that an incorrect orbital remap is detectable.
    """
    nd = ns * norb_phys
    rng = np.random.default_rng(12345)
    tab = np.zeros((cellvol, nd, nd), dtype=np.complex128)
    # one spin-block (norb_phys x norb_phys), distinct per R, Hermitian
    for r in range(cellvol):
        blk = rng.standard_normal((norb_phys, norb_phys))
        blk = 0.5 * (blk + blk.T)  # symmetric -> Hermitian (real)
        # make orbital 0 and 1 clearly distinct so a swap is visible
        blk[0, 0] += 2.0
        blk[1, 1] -= 1.0
        for s in range(ns):
            off = s * norb_phys
            tab[r, off:off + norb_phys, off:off + norb_phys] = blk
    return tab


def _spin_block_to_interleaved(tab_sb, norb_phys, ns=2):
    """interleaved[2*orb+spin] <- spin_block[spin*norb_phys+orb] on both axes."""
    nd = tab_sb.shape[-1]
    # forward map P: interleaved index k = 2*orb + spin from spin-block j = s*norb+orb
    # we want tab_il[k] = tab_sb[j] where k corresponds to same (orb,spin).
    # spin-block j -> (s, orb); interleaved k = 2*orb + s. Build perm: il_from_sb[j]=k
    il_from_sb = np.empty(nd, dtype=int)
    for s in range(ns):
        for orb in range(norb_phys):
            j = s * norb_phys + orb
            k = 2 * orb + s
            il_from_sb[j] = k
    tab_il = np.zeros_like(tab_sb)
    tab_il[:, il_from_sb, :] = tab_sb
    tab_il[:, :, il_from_sb] = tab_il.copy()
    return tab_il


def _run(info_mode_extra, inter, green_init_file=None, subshape=(1, 1, 1)):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [CELLVOL, 1, 1],
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
    if green_init_file is not None:
        input_block["green_init"] = green_init_file
    info_file = {
        "input": input_block,
        "output": {"path_to_output": OUTPUT_DIR},
    }
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    read_io = read_input_k.QLMSkInput(info_file["input"])
    ham = read_io.get_param("ham")
    solver = solver_rpa.RPA(ham, {}, info_mode)
    green = read_io.get_param("green")
    info_in = solver.read_init(info_file["input"])
    for k, v in info_in.items():
        green[k] = v
    solver.solve(green, OUTPUT_DIR)
    return solver, green


class TestRPASOGreenInitMultiOrbital(unittest.TestCase):
    SB_NAME = "so_green_init_sb_fixture.npz"
    IL_NAME = "so_green_init_il_fixture.npz"

    def tearDown(self):
        for n in (self.SB_NAME, self.IL_NAME):
            p = os.path.join(INPUT_DIR, n)
            if os.path.exists(p):
                os.remove(p)

    def _write_fixtures(self):
        tab_sb = _spin_block_green(CELLVOL, NORB_PHYS)
        tab_il = _spin_block_to_interleaved(tab_sb, NORB_PHYS)
        np.savez(os.path.join(INPUT_DIR, self.SB_NAME), green=tab_sb)
        np.savez(os.path.join(INPUT_DIR, self.IL_NAME), green=tab_il)
        return tab_sb, tab_il

    def test_so_green_init_matches_nonso(self):
        """SO green_init (interleaved) must reproduce the non-SO (spin-block) run."""
        self._write_fixtures()
        # non-SO reference: spin-block green_init fed to a ns=2 run
        s_nonso, g_nonso = _run(
            {"enable_spin_orbital": False},
            {"Geometry": "geom_2orb.dat",
             "Transfer": "transfer_nonso_2orb.dat",
             "CoulombIntra": "coulombintra_2orb.dat"},
            green_init_file=self.SB_NAME,
        )
        # under test: interleaved green_init fed to the equivalent SO run
        s_so, g_so = _run(
            {"enable_spin_orbital": True},
            {"Geometry": "geom_so_2orb.dat",
             "Transfer": "transfer_so_2orb.dat",
             "CoulombIntra": "coulombintra_2orb.dat"},
            green_init_file=self.IL_NAME,
        )
        self.assertTrue(
            np.allclose(g_nonso["chi0q"], g_so["chi0q"], atol=1e-10),
            "SO green_init chi0q must match the equivalent non-SO chi0q",
        )
        self.assertTrue(
            np.allclose(g_nonso["chiq"], g_so["chiq"], atol=1e-10),
            "SO green_init chiq must match the equivalent non-SO chiq",
        )


class TestRPASOGreenInitSublatticeFold(unittest.TestCase):
    IL_NAME = "so_green_init_fold_fixture.npz"

    def tearDown(self):
        p = os.path.join(INPUT_DIR, self.IL_NAME)
        if os.path.exists(p):
            os.remove(p)

    def _write_fixture(self):
        tab_sb = _spin_block_green(CELLVOL, NORB_PHYS)
        tab_il = _spin_block_to_interleaved(tab_sb, NORB_PHYS)
        np.savez(os.path.join(INPUT_DIR, self.IL_NAME), green=tab_il)

    @staticmethod
    def _uniform_q0_per_site(chiq, n_sites):
        no = chiq.shape[2]
        s = 0j
        for a in range(no):
            for b in range(no):
                s += chiq[:, 0, a, a, b, b].sum()
        return s / n_sites

    def test_folded_chiq_matches_unfolded(self):
        self._write_fixture()
        inter = {"Geometry": "geom_so_2orb.dat",
                 "Transfer": "transfer_so_2orb.dat",
                 "CoulombIntra": "coulombintra_2orb.dat"}
        _su, g_unfold = _run({"enable_spin_orbital": True}, inter,
                             green_init_file=self.IL_NAME, subshape=(1, 1, 1))
        _sf, g_fold = _run({"enable_spin_orbital": True}, inter,
                           green_init_file=self.IL_NAME, subshape=(2, 1, 1))
        u = self._uniform_q0_per_site(g_unfold["chiq"], 1)
        f = self._uniform_q0_per_site(g_fold["chiq"], 2)
        self.assertAlmostEqual(u.real, f.real, places=10)
        self.assertAlmostEqual(u.imag, f.imag, places=10)


if __name__ == "__main__":
    unittest.main()
