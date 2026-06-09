#!/usr/bin/env python3
"""2D sublattice-fold invariance for RPA._reshape_green (fix B2).

Companion to tests/test_rpa_trans_mod.py, which only exercises a 1D fold
(SubShape [2,1,1] on CellShape [8,1,1]). For a 1D chain the C-order and
Fortran-order flat-SITE index COINCIDE (only x varies), so a site-index
ordering bug in RPA._reshape_green is invisible there.

RPA._reshape_green folds a real-space (cellvol, nd0, nd0) H0/Green array
(from trans_mod / green_init input) into the supercell layout. The resulting
array is then reshape(nx,ny,nz)-FFT'd, which interprets the flat SITE index
in C-order (z fastest). _reshape_green must therefore emit its supercell SITE
index in C-order too -- matching UHFk._deflate_green's _pack_site
(iz + Nz*(iy + Ny*ix)). It previously used Fortran-order (x fastest), so each
supercell site got the wrong FFT phase whenever two fold axes were non-trivial
AND the super-grid was non-square.

Discriminating gate (this test): a NON-SQUARE super-grid is required. CellShape
[8,4,1] with SubShape [2,2,1] gives super-grid Nx=4, Ny=2 -- C-order and
Fortran-order genuinely differ (a square super-grid only transposes the q-grid,
leaving the H0(k) eigenvalue SET unchanged, which is why a 4x4->2x2 fold cannot
detect the bug). With an x/y-asymmetric translation-invariant trans_mod:

  - the folded H0(k) eigenvalue spectrum must equal the unfolded one (band
    folding is unitary), AND
  - the interacting chiq uniform q=0 response per physical site must be
    fold-invariant.

With the old Fortran-order code the folded H0(k) spectrum is wrong by O(0.6);
with the C-order fix both gates hold to machine precision.
"""

import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


INPUT_DIR = "tests/rpa/input"
OUTPUT_DIR = "tests/rpa/output"

# Non-square super-grid under SubShape [2,2,1]: Nx=4, Ny=2.
CELL = [8, 4, 1]
SUBVOL_2D = 4  # Bx*By*Bz for SubShape [2,2,1]


def _make_trans_mod_2d_npz(path, cell=CELL, norb_phys=2, ns=2):
    """Translation-invariant, Hermitian, spin-independent real-space H0 on a 2D
    lattice, hopping in BOTH x and y with x/y-ASYMMETRIC blocks.

    Shape (cellvol, nd0, nd0), nd0 = ns*norb_phys, spin-block layout
    index = s*norb_phys + a. The real-space cell axis is flattened in C-order
    (z fastest): flat = rz + Nz*(ry + Ny*rx), matching RPA's reshape(nx,ny,nz).
    The x and y hopping blocks are deliberately different so that any x<->y
    site mis-ordering changes the H0(k) spectrum.
    """
    Nx, Ny, Nz = cell
    cellvol = Nx * Ny * Nz
    nd0 = ns * norb_phys
    tab = np.zeros((cellvol, nd0, nd0), dtype=np.complex128)

    h_onsite = np.array([[0.3, 0.1],
                         [0.1, -0.2]], dtype=np.complex128)   # R=0, Hermitian
    h_px = np.array([[-0.5, 0.07],
                     [0.04, -0.45]], dtype=np.complex128)     # R=+x
    h_py = np.array([[0.9, -0.2],
                     [0.15, 0.8]], dtype=np.complex128)       # R=+y (asymmetric)

    def flat(rx, ry, rz):
        return (rz % Nz) + Nz * ((ry % Ny) + Ny * (rx % Nx))

    def put(R, block):
        rx, ry, rz = R
        idx = flat(rx, ry, rz)
        for s in range(ns):
            sl = slice(s * norb_phys, (s + 1) * norb_phys)
            tab[idx][sl, sl] += block

    put((0, 0, 0), h_onsite)
    put((1, 0, 0), h_px)
    put((-1, 0, 0), h_px.conj().T)   # Hermiticity H(-R)=H(R)^dagger
    put((0, 1, 0), h_py)
    put((0, -1, 0), h_py.conj().T)

    np.savez(path, trans_mod=tab)
    return path


def _build_solver(subshape, trans_mod_file):
    info_mode = {
        "mode": "RPA",
        "param": {
            "T": 2.0,
            "filling": 0.5,
            "CellShape": list(CELL),
            "SubShape": list(subshape),
            "Nmat": 16,
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
    return solver, green, info_file["input"], info_in


def _run(subshape, trans_mod_file):
    solver, green, _inp, info_in = _build_solver(subshape, trans_mod_file)
    for k, v in info_in.items():
        green[k] = v
    solver.solve(green, OUTPUT_DIR)
    return solver, green


def _h0k_spectrum(subshape, trans_mod_file):
    """Sorted eigenvalue spectrum of the folded q-space H0 read from trans_mod.

    This is the band structure RPA feeds into chi0q; it is invariant under any
    basis permutation but exposes a wrong supercell-site -> FFT-phase mapping.
    """
    _solver, _green, _inp, info_in = _build_solver(subshape, trans_mod_file)
    H0 = info_in["trans_mod"]  # q-space (nvol_super, nd, nd)
    return np.sort(np.linalg.eigvalsh(H0).ravel())


def _uniform_q0_per_site(chiq, n_sites):
    no = chiq.shape[2]
    s = 0j
    for a in range(no):
        for b in range(no):
            s += chiq[:, 0, a, a, b, b].sum()
    return s / n_sites


class TestRPATransMod2DSublatticeFold(unittest.TestCase):
    TM_NAME = "trans_mod_2d_fixture.npz"

    def tearDown(self):
        p = os.path.join(INPUT_DIR, self.TM_NAME)
        if os.path.exists(p):
            os.remove(p)

    def test_2d_folded_H0k_spectrum_matches_unfolded(self):
        """Band-folding gate: the folded H0(k) spectrum must equal the unfolded
        one. FAILS by O(0.6) with the old Fortran-order _reshape_green site
        index; PASSES to machine precision with the C-order fix."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        tm_path = os.path.join(INPUT_DIR, self.TM_NAME)
        _make_trans_mod_2d_npz(tm_path)

        eu = _h0k_spectrum((1, 1, 1), self.TM_NAME)
        ef = _h0k_spectrum((2, 2, 1), self.TM_NAME)

        self.assertEqual(eu.size, ef.size)
        np.testing.assert_allclose(eu, ef, atol=1.0e-10)

    def test_2d_folded_chiq_with_trans_mod_matches_unfolded(self):
        """Physics gate: interacting chiq uniform q=0 response per physical site
        is fold-invariant for the same translation-invariant H0."""
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        tm_path = os.path.join(INPUT_DIR, self.TM_NAME)
        _make_trans_mod_2d_npz(tm_path)

        _su, g_unfold = _run((1, 1, 1), self.TM_NAME)
        _sf, g_fold = _run((2, 2, 1), self.TM_NAME)

        # Per-physical-site normalization (mirrors test_rpa_trans_mod.py):
        # divide by the sublattice volume.
        u = _uniform_q0_per_site(g_unfold["chiq"], 1)
        f = _uniform_q0_per_site(g_fold["chiq"], SUBVOL_2D)

        self.assertAlmostEqual(u.real, f.real, places=10)
        self.assertAlmostEqual(u.imag, f.imag, places=10)


if __name__ == "__main__":
    unittest.main()
