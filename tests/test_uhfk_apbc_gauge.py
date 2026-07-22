"""Inverse gauge transform applied at the per-site Green output path.

Builds a synthetic internal Green tensor on a 1D L=4 chain, calls
``_save_greenone`` for both PBC and APBC variants, and verifies:

  G_ij (output) = inverse_gauge_phase(r_i, r_j, theta, L) * tilde G_ij

per (i, j) pair. Uses the same ``__new__`` stub pattern as
``tests/test_uhfk_export.py`` so no SCF or loader machinery is required.
"""
from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pytest

from hwave.solver.uhfk import UHFk
from hwave.solver._apbc_phase import inverse_gauge_phase


L = 4  # 1D chain length used throughout this test


def _make_stub_solver(boundary_theta):
    """Build a minimal UHFk instance suitable for _save_greenone."""
    s = UHFk.__new__(UHFk)
    s.enable_spin_orbital = False
    s.has_sublattice = False
    s.cellshape = (L, 1, 1)
    s.cellvol = L
    s.norb = 1
    s.ns = 2
    s.nd = s.ns * s.norb  # = 2
    # Internal "tilde" Green: deterministic but non-trivial.
    # Shape per uhfk convention: (lvol, ns, norb, ns, norb)
    rng = np.random.default_rng(20260626)
    s.Green = (
        rng.standard_normal((L, 2, 1, 2, 1))
        + 1j * rng.standard_normal((L, 2, 1, 2, 1))
    ).astype(np.complex128)
    s.boundary_theta = tuple(boundary_theta)
    s.boundary_periodic = all(t == 0.0 for t in boundary_theta)
    return s


def _greenone_pairs():
    """All (i, s, j, t) pairs we want written out."""
    return [(i, sp, j, t) for i in range(L) for j in range(L)
            for sp in range(2) for t in range(2)]


def _geom_info():
    # site2vec maps site id -> (rx, ry, rz, orb). For a 1D chain of L sites.
    return {"site2vec": {i: (i, 0, 0, 0) for i in range(L)}}


def _read_greenone(path: Path):
    out: dict[tuple[int, int, int, int], complex] = {}
    with open(path) as fr:
        for line in fr:
            parts = line.split()
            if len(parts) < 6:
                continue
            i, sp, j, t = (int(parts[k]) for k in range(4))
            re, im = float(parts[4]), float(parts[5])
            out[(i, sp, j, t)] = complex(re, im)
    return out


def test_pbc_path_is_byte_identical_to_pre_apbc_behavior():
    """Theta = 0 must produce the same file the un-modified _save_greenone would."""
    s = _make_stub_solver(boundary_theta=(0.0, 0.0, 0.0))
    with tempfile.TemporaryDirectory() as tmp:
        fname = Path(tmp) / "greenone.dat"
        s._save_greenone(
            str(fname),
            {"onebodyg_uhf": _greenone_pairs(), "geometry_uhf": _geom_info()},
        )
        out = _read_greenone(fname)
    # For PBC the inverse phase is 1 -> output equals direct gr[idx,s,a,t,b].
    gr = s.Green
    for (i, sp, j, t), v_out in out.items():
        idx = (j - i + L) % L
        expected = complex(gr[idx, sp, 0, t, 0])
        np.testing.assert_allclose(v_out, expected, atol=1e-12,
                                   err_msg=f"PBC mismatch at (i,j,s,t)=({i},{j},{sp},{t})")


def test_apbc_off_diagonal_carries_inverse_gauge_phase():
    theta = np.array([np.pi, 0.0, 0.0])
    L_vec = np.array([L, 1, 1], dtype=np.float64)
    s = _make_stub_solver(boundary_theta=tuple(theta))
    with tempfile.TemporaryDirectory() as tmp:
        fname = Path(tmp) / "greenone.dat"
        s._save_greenone(
            str(fname),
            {"onebodyg_uhf": _greenone_pairs(), "geometry_uhf": _geom_info()},
        )
        out = _read_greenone(fname)
    gr_tilde = s.Green
    for (i, sp, j, t), v_out in out.items():
        idx = (j - i + L) % L
        v_tilde = complex(gr_tilde[idx, sp, 0, t, 0])
        r_i = np.array([i, 0, 0], dtype=np.float64)
        r_j = np.array([j, 0, 0], dtype=np.float64)
        expected = v_tilde * inverse_gauge_phase(r_i, r_j, theta, L_vec)
        np.testing.assert_allclose(v_out, expected, atol=1e-12,
                                   err_msg=f"APBC mismatch at (i,j,s,t)=({i},{j},{sp},{t})")


def test_apbc_diagonal_pairs_unchanged():
    """When r_i == r_j the inverse phase is 1, so output equals tilde value."""
    theta = np.array([np.pi, np.pi, 0.0])
    s = _make_stub_solver(boundary_theta=tuple(theta))
    with tempfile.TemporaryDirectory() as tmp:
        fname = Path(tmp) / "greenone.dat"
        s._save_greenone(
            str(fname),
            {"onebodyg_uhf": _greenone_pairs(), "geometry_uhf": _geom_info()},
        )
        out = _read_greenone(fname)
    gr = s.Green
    for (i, sp, j, t), v_out in out.items():
        if i != j:
            continue
        idx = 0  # j - i = 0
        expected = complex(gr[idx, sp, 0, t, 0])
        np.testing.assert_allclose(v_out, expected, atol=1e-12,
                                   err_msg=f"diag mismatch at i={i}, s={sp}, t={t}")


def test_end_to_end_physical_green_matches_analytic_1d():
    """End-to-end gauge-convention check.

    1) Analytic physical Green for closed-shell half-filled L=4 1D free
       fermion APBC, computed directly from physical k = (2 pi n + pi) / L.
    2) H-wave production chain: _make_ham_trans applies the gauge phase;
       we diagonalize, identify occupied k_n, build tilde G as a plane-wave
       sum over occupied gauge-basis k_n; call _save_greenone (which
       applies the inverse gauge); read the output file.
    3) Compare output to analytic Green. The expected value never calls
       any APBC helper, so a wrong gauge sign or boundary-crossing
       displacement convention cannot hide.
    """
    L_chain = 4
    t = 1.0
    theta_x = np.pi  # APBC in x

    # ---- 1) Analytic physical Green ----
    # Allowed k_phys = (2 pi n + theta) / L; energy eps = -2t cos(k_phys).
    # Closed-shell half-filling: pick L/2 lowest-energy k_phys per spin.
    k_phys_all = np.array([(2 * np.pi * n + theta_x) / L_chain for n in range(L_chain)])
    eps_all = -2.0 * t * np.cos(k_phys_all)
    occ_phys = np.argsort(eps_all)[: L_chain // 2]
    k_phys_occ = k_phys_all[occ_phys]
    g_analytic = np.zeros((L_chain, L_chain), dtype=complex)
    for i in range(L_chain):
        for j in range(L_chain):
            g_analytic[i, j] = (1.0 / L_chain) * np.sum(
                np.exp(1j * k_phys_occ * (j - i))
            )

    # ---- 2) H-wave production chain ----
    s = UHFk.__new__(UHFk)
    s.shape = (L_chain, 1, 1)
    s.cellshape = (L_chain, 1, 1)
    s.cellvol = L_chain
    s.nvol = L_chain
    s.norb = 1
    s.ns = 2
    s.nd = 2
    s.enable_spin_orbital = False
    s.has_sublattice = False
    s.boundary_theta = (theta_x, 0.0, 0.0)
    s.boundary_periodic = False
    s.param_ham = {
        "Transfer": {((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t},
    }
    # Production wires the APBC gauge phase into Transfer in _init_interaction
    # (pre-fold, signed irvec), not in _make_ham_trans. This stub bypasses
    # __init__, so mimic that step explicitly before FFT-ing.
    from hwave.solver._apbc_phase import transfer_phase
    theta_arr = np.array(s.boundary_theta, dtype=np.float64)
    L_arr = np.array(s.cellshape, dtype=np.float64)
    s.param_ham["Transfer"] = {
        k: v * transfer_phase(np.asarray(k[0], dtype=np.float64), theta_arr, L_arr)
        for k, v in s.param_ham["Transfer"].items()
    }
    s._make_ham_trans()
    # ham_trans shape: (nvol, nd, nd) = (L, 2, 2). For our spin-symmetric
    # Hamiltonian, each H(k_n) is a scalar * eye(2) -> two degenerate eigs.
    energies_per_k = np.array(
        [np.linalg.eigvalsh(s.ham_trans[k_n])[0] for k_n in range(L_chain)]
    )
    # Per spin occupation: lowest L/2 gauge-basis k_n.
    occ_kn = np.argsort(energies_per_k)[: L_chain // 2]

    # Build tilde G from gauge-basis plane waves (a non-circular construction:
    # we use only the gauge-basis k-grid and the eigenvalue ordering from
    # ham_trans, NOT inverse_gauge_phase).
    #
    # H-wave uses numpy ifftn(norm='forward'): hat T(K_idx) = sum_R T(R)
    # exp(+2 pi i K_idx R / L). With H = -t (c^dag_{x+1} c_x + h.c.), a plane
    # wave e^{+i k x} has eigenvalue -2t cos(k); this means H-wave's K_idx
    # corresponds to physical momentum k = -2 pi K_idx / L (negative sign),
    # i.e., the gauge-basis stationary state at K_idx has wave function
    #   tilde phi(x) = e^{-i 2 pi K_idx x / L} / sqrt(L).
    # Hence
    #   tilde G[delta] = sum_K tilde phi*(0) tilde phi(delta)
    #                  = (1/L) sum_K exp(-i 2 pi K delta / L).
    tilde_g = np.zeros((L_chain, 2, 1, 2, 1), dtype=complex)
    for delta in range(L_chain):
        val = (1.0 / L_chain) * sum(
            np.exp(-1j * 2.0 * np.pi * int(k_n) * delta / L_chain) for k_n in occ_kn
        )
        for sp in range(2):
            tilde_g[delta, sp, 0, sp, 0] = val
    s.Green = tilde_g

    # ---- 3) Save and load via production path ----
    pairs = [
        (i, sp, j, sp) for i in range(L_chain) for j in range(L_chain) for sp in range(2)
    ]
    geom = {"site2vec": {i: (i, 0, 0, 0) for i in range(L_chain)}}
    with tempfile.TemporaryDirectory() as tmp:
        fname = Path(tmp) / "greenone.dat"
        s._save_greenone(str(fname), {"onebodyg_uhf": pairs, "geometry_uhf": geom})
        out = _read_greenone(fname)

    # ---- 4) Compare ----
    # Spin-diagonal output must match analytic G (independent of sigma).
    for i in range(L_chain):
        for j in range(L_chain):
            expected = g_analytic[i, j]
            for sp in range(2):
                got = out[(i, sp, j, sp)]
                np.testing.assert_allclose(
                    got, expected, atol=1e-10,
                    err_msg=f"end-to-end mismatch (i,j,s)=({i},{j},{sp}): "
                            f"got {got}, want {expected}",
                )
    # (Spin-off-diagonal (s != t) entries are not requested in `pairs` above.)
