#!/usr/bin/env python3
"""trans_mod + sublattice folding invariance for RPA.

RPA accepts an optional ``trans_mod`` input: a real-space, single-particle
Hamiltonian modification that REPLACES H0 in ``_calc_epsilon_k``. When the
lattice is folded (SubShape volume > 1) the read path folds the real-space
``trans_mod`` (shape ``(cellvol, nd0, nd0)``, spin-block ``s*norb_orig+a``)
into the supercell layout via ``RPA._reshape_green``.

Physics gate: a ``trans_mod`` representing the SAME translation-invariant
physical H0 must give a fold-invariant interacting susceptibility ``chiq``
(uniform q=0 response per physical site) whether run unfolded (SubShape
[1,1,1]) or folded (SubShape [2,1,1]). Both runs read the SAME .npz file;
the unfolded run uses it directly, the folded run folds it. A
translation-invariant H0 must therefore yield identical physics.

This exercises the (pre-existing, since 2023) sublattice trans_mod read path:
the mis-wired ``self.lattice._reshape_green`` -> ``self._reshape_green`` and the
previously unassigned ``self.norb_orig``.
"""

import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


INPUT_DIR = "tests/rpa/input"
OUTPUT_DIR = "tests/rpa/output"


def _make_trans_mod_npz(path, ncell=8, norb_phys=2, ns=2):
    """Translation-invariant, Hermitian, spin-independent real-space H0.

    Shape (ncell, nd0, nd0) with nd0 = ns*norb_phys, spin-block layout
    index = s*norb_phys + a. On-site orbital energies at R=0 plus a small
    real NN hopping at R=+/-1 (identical for every cell -> translation
    invariant). Spin-independent (same block for both spins).
    """
    nd0 = ns * norb_phys
    tab = np.zeros((ncell, nd0, nd0), dtype=np.complex128)

    # single-spin (norb_phys x norb_phys) real-space hopping, by displacement R
    h_onsite = np.array([[0.3, 0.1],
                         [0.1, -0.2]], dtype=np.complex128)   # R=0, Hermitian
    h_p1 = np.array([[-0.5, 0.07],
                     [0.04, -0.45]], dtype=np.complex128)     # R=+1

    def put(R, block):
        # embed orbital block into both spin diagonal blocks (spin-independent)
        for s in range(ns):
            sl = slice(s * norb_phys, (s + 1) * norb_phys)
            tab[R % ncell][sl, sl] += block

    put(0, h_onsite)
    put(1, h_p1)
    put(-1, h_p1.conj().T)  # Hermiticity: H(-R) = H(R)^dagger

    np.savez(path, trans_mod=tab)
    return path


def _run(subshape, trans_mod_file):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": [8, 1, 1],
            "SubShape": list(subshape),
            "Nmat": 32,
        },
        "enable_spin_orbital": False,
        "calc_scheme": "general",
    }
    inter = {
        "path_to_input": INPUT_DIR,
        "Geometry": "geom_2orb.dat",
        "Transfer": "transfer_nonso_2orb.dat",
        "CoulombIntra": "coulombintra_2orb.dat",
    }
    info_file = {
        "input": {
            "path_to_input": INPUT_DIR,
            "interaction": inter,
            "trans_mod": trans_mod_file,
        },
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


def _uniform_q0_per_site(chiq, n_sites):
    no = chiq.shape[2]
    s = 0j
    for a in range(no):
        for b in range(no):
            s += chiq[:, 0, a, a, b, b].sum()
    return s / n_sites


class TestRPATransModSublatticeFold(unittest.TestCase):
    TM_NAME = "trans_mod_fixture.npz"

    def tearDown(self):
        # The fixture must live under path_to_input (read_init joins it there);
        # remove the generated artifact so it is not left untracked.
        p = os.path.join(INPUT_DIR, self.TM_NAME)
        if os.path.exists(p):
            os.remove(p)

    def test_folded_chiq_with_trans_mod_matches_unfolded(self):
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        tm_name = self.TM_NAME
        tm_path = os.path.join(INPUT_DIR, tm_name)
        _make_trans_mod_npz(tm_path)

        _su, g_unfold = _run((1, 1, 1), tm_name)
        _sf, g_fold = _run((2, 1, 1), tm_name)

        u = _uniform_q0_per_site(g_unfold["chiq"], 1)
        f = _uniform_q0_per_site(g_fold["chiq"], 2)

        self.assertAlmostEqual(u.real, f.real, places=10)
        self.assertAlmostEqual(u.imag, f.imag, places=10)


if __name__ == "__main__":
    unittest.main()
