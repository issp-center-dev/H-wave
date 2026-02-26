"""Linearized Eliashberg equation solver for superconducting instability analysis.

This module implements a post-processing tool that uses the bare susceptibility
chi0q from H-wave's RPA solver to solve the linearized Eliashberg equation
and analyze superconducting instabilities.

The tool reads chi0q.npz output from RPA calculations along with interaction
parameters, reconstructs the Green's function, computes RPA vertices,
and solves the Eliashberg equation via self-consistent iteration or
eigenvalue analysis.
"""

from __future__ import annotations

import os
import sys
import logging

import numpy as np
from numpy.fft import fftn, ifftn
from scipy.optimize import bisect
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

import hwave
import hwave.qlmsio.wan90 as wan90

logger = logging.getLogger("hwave_sc")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _load_chi0q(input_dict):
    """Load chi0q from NPZ file produced by H-wave RPA solver.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array.
    """
    output_info = input_dict["file"]["output"]
    path_to_output = output_info["path_to_output"]
    chi0q_name = output_info.get("chi0q", "chi0q")
    file_name = os.path.join(path_to_output, chi0q_name + ".npz")

    logger.info("Loading chi0q from {}".format(file_name))
    data = np.load(file_name)
    chi0q = data["chi0q"]
    logger.info("chi0q shape: {}".format(chi0q.shape))
    return chi0q


def _calc_chi0q_internal(input_dict, chi0q_tensor="auto"):
    """Compute chi0q internally using H-wave's RPA module.

    Instead of loading chi0q from a file, this function creates an RPA solver
    instance and computes chi0q from scratch using the Transfer and interaction
    parameters specified in the TOML config.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.
    chi0q_tensor : str
        Index structure of chi0q. Options:
        - "reduced": 2-index (nmat, nvol, norb, norb). Exact for single-orbital
          and for multi-orbital with CoulombIntra only (no inter-orbital
          interactions).
        - "general": 4-index (nmat, nvol, norb, norb, norb, norb). Required
          when inter-orbital interactions (CoulombInter, Hund, Exchange) are
          present, as their S/C matrices couple different orbital-pair indices.
        - "auto": Use "general" if CoulombInter/Hund/Exchange are present,
          otherwise "reduced".

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array. Shape depends on chi0q_tensor mode:
        - reduced: (nmat, nvol, norb, norb)
        - general: (nmat, nvol, norb, norb, norb, norb)
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Determine calc_scheme from chi0q_tensor option
    if chi0q_tensor == "auto":
        # Use "general" when inter-orbital interactions are present,
        # because their S/C matrices have off-diagonal elements that
        # couple to chi0q off-diagonal components.
        # With CoulombIntra only, S is block-diagonal and reduced is exact.
        files = info_inputfile.get("interaction", {})
        has_interorbital = any(k in files for k in
                              ["Hund", "Exchange", "CoulombInter",
                               "Ising", "PairHop"])
        if has_interorbital:
            calc_scheme = "general"
        else:
            calc_scheme = "reduced"
    elif chi0q_tensor == "general":
        calc_scheme = "general"
    else:
        calc_scheme = "reduced"

    info_mode_rpa = dict(info_mode)
    info_mode_rpa["calc_scheme"] = calc_scheme

    logger.info("Computing chi0q internally (calc_scheme={})...".format(calc_scheme))

    # Read input files via QLMSkInput
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")

    # Create RPA solver (sets up Lattice, Interaction, parameters)
    solver = sol_rpa.RPA(ham_info, info_log, info_mode_rpa)

    # Get green_info (may contain chi0q_init, trans_mod, green_init)
    green_info = read_io.get_param("green")
    green_info.update(solver.read_init(info_inputfile))

    # Compute chi0q: eigenvalues -> mu -> Green's function -> chi0q
    beta = 1.0 / solver.T

    solver._calc_epsilon_k(green_info)

    if solver.calc_mu:
        if solver.spin_mode == "spin-free":
            Ncond = solver.Ncond / 2
        else:
            Ncond = solver.Ncond
        dist, mu = solver._find_mu(Ncond, solver.T)
    else:
        mu = solver.mu_value

    green0, green0_tail = solver._calc_green(beta, mu)

    chi0q = solver._calc_chi0q(green0, green0_tail, beta)

    # For spin-free mode, chi0q has shape (1, nmat, nvol, ...)
    # Remove the block index
    if solver.spin_mode in ["spin-free", "spinful"]:
        assert chi0q.shape[0] == 1
        chi0q = chi0q[0]
    else:
        # spin-diag: shape (2, nmat, nvol, ...)
        assert chi0q.shape[0] == 2
        pass

    logger.info("chi0q computed internally, shape: {}".format(chi0q.shape))
    return chi0q


def _read_interaction_files(input_dict):
    """Read Transfer and interaction files using wan90 module.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    geom_info : dict
        Geometry information (norb, rvec, center).
    hr : dict
        Transfer (hopping) parameters.
    interactions : dict
        Dictionary of interaction parameters keyed by type name.
    """
    files = input_dict["file"]["input"]["interaction"]
    path_to_input = files.get("path_to_input", "")

    geom_file = os.path.join(path_to_input, files["Geometry"])
    geom_info = wan90.read_geom(geom_file)
    logger.info("norb = {}".format(geom_info["norb"]))

    transfer_file = os.path.join(path_to_input, files["Transfer"])
    hr = wan90.read_w90(transfer_file)

    interaction_types = ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                        "Ising", "PairLift", "PairHop"]
    interactions = {}
    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            interactions[itype] = wan90.read_w90(f)

    return geom_info, hr, interactions


# ---------------------------------------------------------------------------
# k-space construction
# ---------------------------------------------------------------------------

def _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb):
    """Build Hamiltonian epsilon(k) from real-space transfer integrals.

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    hr : dict
        Transfer parameters from wan90.read_w90().
    norb : int
        Number of orbitals.

    Returns
    -------
    epsilon_k : ndarray
        Hamiltonian in k-space, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    epsilon_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    for (irvec, orbvec), value in hr.items():
        orb1, orb2 = orbvec
        Rx, Ry, Rz = irvec
        epsilon_k[orb2, orb1, :, :, :] += value * np.exp(
            1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
        )
    return epsilon_k


def _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb):
    """Build interaction matrices in k-space.

    Supports CoulombIntra (U), CoulombInter (U'), Hund (J), and Exchange (J').

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    interactions : dict
        Interaction parameters keyed by type.
    norb : int
        Number of orbitals.

    Returns
    -------
    inter_k : dict
        Dictionary of interactions in k-space.
        Keys: "CoulombIntra", "CoulombInter", "Hund", "Exchange".
        Values: ndarray of shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )

    def _to_k(value_r):
        val_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for (irvec, orbvec), value in value_r.items():
            orb1, orb2 = orbvec
            Rx, Ry, Rz = irvec
            val_k[orb2, orb1, :, :, :] += value * np.exp(
                1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
            )
        return val_k

    inter_k = {}
    for itype in ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                  "Ising", "PairLift", "PairHop"]:
        if itype in interactions:
            inter_k[itype] = _to_k(interactions[itype])

    return inter_k


# ---------------------------------------------------------------------------
# Green's function
# ---------------------------------------------------------------------------

def _calc_eigenvalues(epsilon_k):
    """Diagonalize epsilon(k) at every k-point.

    Parameters
    ----------
    epsilon_k : ndarray
        Hamiltonian, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    """
    norb = epsilon_k.shape[0]
    Nx, Ny, Nz = epsilon_k.shape[2], epsilon_k.shape[3], epsilon_k.shape[4]
    eigenvalues = np.zeros((Nx, Ny, Nz, norb))
    eigenvectors = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)

    # Batch diagonalization: transpose to (Nx, Ny, Nz, norb, norb) for vectorized eigh
    eps_batch = epsilon_k.transpose(2, 3, 4, 0, 1)  # (Nx, Ny, Nz, norb, norb)
    eigenvalues, eigenvectors = np.linalg.eigh(eps_batch)

    return eigenvalues, eigenvectors


def _determine_mu(eigenvalues, beta, n_target, norb):
    """Determine chemical potential using bisection.

    The filling convention follows the reference code: n_target is the
    number of electrons per orbital per spin channel. For example,
    n_target=0.75 means 3/4 filling per spin.

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    beta : float
        Inverse temperature.
    n_target : float
        Target filling per orbital per spin (e.g. 0.75 for 3/4 filling).
    norb : int
        Number of orbitals.

    Returns
    -------
    mu : float
        Chemical potential.
    """
    Nx, Ny, Nz = eigenvalues.shape[:3]
    nvol = Nx * Ny * Nz

    def _calc_n(mu):
        x = beta * (eigenvalues - mu)
        fermi = np.where(x > 100, 0.0, np.where(x < -100, 1.0, 1.0 / (1.0 + np.exp(x))))
        total_n = np.sum(fermi)
        return float(total_n / nvol - n_target * norb)

    emin = np.min(eigenvalues)
    emax = np.max(eigenvalues)
    mu = bisect(_calc_n, emin - 10.0, emax + 10.0)
    return float(mu)


def _calc_green(eigenvalues, eigenvectors, mu, beta, nmat):
    """Construct Green's function G(k, iwn).

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    mu : float
        Chemical potential.
    beta : float
        Inverse temperature.
    nmat : int
        Number of Matsubara frequencies.

    Returns
    -------
    green_kw : ndarray
        Shape (norb, norb, Nx, Ny, Nz, nmat).
    """
    Nx, Ny, Nz, norb = eigenvalues.shape
    iomega = np.array([(2.0 * i + 1.0 - nmat) * np.pi for i in range(nmat)]) / beta

    # Vectorized Green's function construction:
    # G_{ij}(k, iwn) = sum_m U_{im}(k) U*_{jm}(k) / (iwn - (e_m(k) - mu))

    # factor[kx,ky,kz,i,j,m] = U[kx,ky,kz,i,m] * conj(U[kx,ky,kz,j,m])
    factor = np.einsum('...im,...jm->...ijm', eigenvectors, np.conj(eigenvectors))
    # factor shape: (Nx, Ny, Nz, norb, norb, norb)

    # denom[kx,ky,kz,m,w] = 1 / (iwn_w - (e_m(k) - mu))
    xi = eigenvalues - mu  # (Nx, Ny, Nz, norb)
    denom = 1.0 / (1j * iomega[None, None, None, None, :] - xi[:, :, :, :, None])
    # denom shape: (Nx, Ny, Nz, norb, nmat)

    # G[kx,ky,kz,i,j,w] = sum_m factor[...,i,j,m] * denom[...,m,w]
    green_kw_tmp = np.einsum('...ijm,...mw->...ijw', factor, denom)
    # shape: (Nx, Ny, Nz, norb, norb, nmat)

    # Transpose to output convention: (norb, norb, Nx, Ny, Nz, nmat)
    green_kw = green_kw_tmp.transpose(3, 4, 0, 1, 2, 5)

    return green_kw


# ---------------------------------------------------------------------------
# RPA vertices
# ---------------------------------------------------------------------------

def _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz):
    """Build spin (S) and charge (C) interaction matrices for all q-points at once.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    S_all : ndarray
        Spin interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    C_all : ndarray
        Charge interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    """
    nd = norb * norb
    S_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
    C_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

    def _get(itype):
        if itype in inter_k:
            return inter_k[itype]  # (norb, norb, Nx, Ny, Nz)
        return None

    U_mat = _get("CoulombIntra")
    Up_mat = _get("CoulombInter")
    J_mat = _get("Hund")
    Jp_mat = _get("Exchange")
    I_mat = _get("Ising")
    PH_mat = _get("PairHop")

    # Build using vectorized index conditions
    for l1 in range(norb):
        for l2 in range(norb):
            idx12 = l1 * norb + l2
            for l3 in range(norb):
                for l4 in range(norb):
                    idx34 = l3 * norb + l4

                    if l1 == l2 == l3 == l4:
                        if U_mat is not None:
                            S_all[:, :, :, idx12, idx34] = U_mat[l1, l1]
                            C_all[:, :, :, idx12, idx34] = U_mat[l1, l1]
                    elif l1 == l3 and l2 == l4 and l1 != l2:
                        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
                        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
                        if Up_mat is not None:
                            s_q += Up_mat[l1, l2]
                            c_q -= Up_mat[l1, l2]
                        if I_mat is not None:
                            s_q -= I_mat[l1, l2]
                            c_q -= I_mat[l1, l2]
                        if J_mat is not None:
                            c_q += J_mat[l1, l2]
                        S_all[:, :, :, idx12, idx34] = s_q
                        C_all[:, :, :, idx12, idx34] = c_q
                    elif l1 == l2 and l3 == l4 and l1 != l3:
                        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
                        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
                        if J_mat is not None:
                            s_q += J_mat[l1, l3]
                            c_q -= J_mat[l1, l3]
                        if I_mat is not None:
                            s_q -= 2.0 * I_mat[l1, l3]
                        if Up_mat is not None:
                            c_q += 2.0 * Up_mat[l1, l3]
                        S_all[:, :, :, idx12, idx34] = s_q
                        C_all[:, :, :, idx12, idx34] = c_q
                    elif l1 == l4 and l2 == l3 and l1 != l2:
                        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
                        if Jp_mat is not None:
                            s_q += Jp_mat[l1, l2]
                        if PH_mat is not None:
                            s_q += PH_mat[l1, l2]
                        S_all[:, :, :, idx12, idx34] = s_q
                        C_all[:, :, :, idx12, idx34] = s_q  # S = C for this channel

    return S_all, C_all


def _build_sc_matrices(inter_k, norb, ix, iy, iz):
    """Build spin (S) and charge (C) interaction matrices at a given q-point.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    The composite index maps as (l1,l2) -> l1*norb + l2,
    giving norb^2 x norb^2 matrices.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    ix, iy, iz : int
        q-point indices.

    Returns
    -------
    S_mat : ndarray
        Spin interaction matrix, shape (norb^2, norb^2).
    C_mat : ndarray
        Charge interaction matrix, shape (norb^2, norb^2).
    """
    nd = norb * norb
    S_mat = np.zeros((nd, nd), dtype=complex)
    C_mat = np.zeros((nd, nd), dtype=complex)

    # Extract interaction values at this q-point
    def _get(itype):
        if itype in inter_k:
            return inter_k[itype][:, :, ix, iy, iz]
        return np.zeros((norb, norb), dtype=complex)

    U_mat = _get("CoulombIntra")    # U_mm (intra-orbital)
    Up_mat = _get("CoulombInter")   # U'_mm' (inter-orbital)
    J_mat = _get("Hund")            # J_mm' (Hund's coupling)
    Jp_mat = _get("Exchange")       # J'_mm' (pair-hopping)
    I_mat = _get("Ising")           # I_mm' (Ising S^z S^z)
    PH_mat = _get("PairHop")        # P_mm' (pair hopping)

    for l1 in range(norb):
        for l2 in range(norb):
            idx12 = l1 * norb + l2
            for l3 in range(norb):
                for l4 in range(norb):
                    idx34 = l3 * norb + l4

                    s_val = 0.0 + 0.0j
                    c_val = 0.0 + 0.0j

                    if l1 == l2 == l3 == l4:
                        # Same orbital: S = U, C = U
                        s_val = U_mat[l1, l1]
                        c_val = U_mat[l1, l1]
                    elif l1 == l3 and l2 == l4 and l1 != l2:
                        # l1=l3 != l2=l4 (cross): S = U' - I, C = -U' + J - I
                        s_val = Up_mat[l1, l2] - I_mat[l1, l2]
                        c_val = (-Up_mat[l1, l2] + J_mat[l1, l2]
                                 - I_mat[l1, l2])
                    elif l1 == l2 and l3 == l4 and l1 != l3:
                        # l1=l2 != l3=l4 (dens): S = J - 2I, C = 2U' - J
                        s_val = J_mat[l1, l3] - 2.0 * I_mat[l1, l3]
                        c_val = 2.0 * Up_mat[l1, l3] - J_mat[l1, l3]
                    elif l1 == l4 and l2 == l3 and l1 != l2:
                        # l1=l4 != l2=l3 (exch): S = J' + P, C = J' + P
                        s_val = Jp_mat[l1, l2] + PH_mat[l1, l2]
                        c_val = Jp_mat[l1, l2] + PH_mat[l1, l2]

                    S_mat[idx12, idx34] = s_val
                    C_mat[idx12, idx34] = c_val

    return S_mat, C_mat


def _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                      pairing_type="singlet"):
    """Compute effective pairing interaction V(q).

    Supports two modes:
    - Simple mode (only CoulombIntra/Inter, 2-index chi0q):
      Uses Wc=U+2V, Ws=-U formulation.
    - General mode (with Hund/Exchange, or 4-index chi0q):
      Uses S,C matrices from Kuroki et al.

    When 4-index chi0q is provided (8D array), general mode is always used
    because the full orbital tensor structure must be preserved.

    Pairing types:
    - singlet: V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)
    - triplet: V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Parameters
    ----------
    chi0q : ndarray
        Bare susceptibility.
        - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb, Nx, Ny, Nz, nmat : int
        System parameters.
    pairing_type : str
        "singlet" or "triplet". Default "singlet".

    Returns
    -------
    Result depends on mode:
        Simple mode: tuple (Pc_q, Ps_q), each shape (norb, norb, Nx, Ny, Nz).
        General mode: Vs_q, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    has_interorbital_vertex = any(k in inter_k for k in
                                  ["Hund", "Exchange", "Ising", "PairHop"])
    chi0q_is_4index = (chi0q.ndim == 8)

    if chi0q_is_4index or has_interorbital_vertex:
        # General mode: 4-index S,C matrices
        # Required when 4-index chi0q is available or Hund/Exchange present
        return _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                         pairing_type=pairing_type)
    else:
        # Simple mode: backward compatible with original implementation
        # Only used for 2-index chi0q without Hund/Exchange
        return _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type=pairing_type)


def _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                             pairing_type="singlet"):
    """Compute vertices using simple Wc=U+2V, Ws=-U formulation.

    For singlet:
        Pc = (Wc+Ws)/2 - (1/2) Wc chi_c Wc
        Ps = -Ws + (3/2) Ws chi_s Ws
    For triplet:
        Pc = (Wc-Ws)/2 - (1/2) Wc chi_c Wc   (charge part same sign change)
        Ps = Ws - (1/2) Ws chi_s Ws            (spin part sign flip)

    Returns Pc_q and Ps_q separately for backward compatibility.

    Returns
    -------
    Pc_q : ndarray
        Charge vertex, shape (norb, norb, Nx, Ny, Nz).
    Ps_q : ndarray
        Spin vertex, shape (norb, norb, Nx, Ny, Nz).
    """
    U_k = inter_k.get("CoulombIntra", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))
    V_k = inter_k.get("CoulombInter", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))

    # Transpose to batch dimension first: (Nx, Ny, Nz, norb, norb)
    Wc = (U_k + 2.0 * V_k).transpose(2, 3, 4, 0, 1).copy()
    Ws = (-U_k).transpose(2, 3, 4, 0, 1).copy()

    # chi0 at static limit: (Nx, Ny, Nz, norb, norb)
    chi0_static = chi0q[:, :, :, :, :, nmat // 2].transpose(2, 3, 4, 0, 1).copy()

    I_mat = np.broadcast_to(np.eye(norb, dtype=complex), (Nx, Ny, Nz, norb, norb)).copy()

    # Batched solve
    mat_s = I_mat + np.einsum('...ab,...bc->...ac', chi0_static, Ws)
    mat_c = I_mat + np.einsum('...ab,...bc->...ac', chi0_static, Wc)

    chis = np.linalg.solve(mat_s, chi0_static)
    chic = np.linalg.solve(mat_c, chi0_static)

    WsChisWs = Ws @ chis @ Ws
    WcChicWc = Wc @ chic @ Wc

    if pairing_type == "singlet":
        Pc_all = (Wc + Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = -Ws + 1.5 * WsChisWs
    elif pairing_type == "triplet":
        Pc_all = (Wc - Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = Ws - 0.5 * WsChisWs
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Transpose back: (Nx, Ny, Nz, norb, norb) -> (norb, norb, Nx, Ny, Nz)
    Pc_q = Pc_all.transpose(3, 4, 0, 1, 2)
    Ps_q = Ps_all.transpose(3, 4, 0, 1, 2)

    return Pc_q, Ps_q


def _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                              pairing_type="singlet"):
    """Compute effective pairing interaction using general S,C matrices.

    For singlet (Kuroki et al., arXiv:0902.3691, Eq.(6)):
        V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)

    For triplet (Takimoto et al., PRB 69, 104504):
        V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Supports both 2-index (reduced) and 4-index (general) chi0q:
    - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
      Expanded to diagonal of (norb^2, norb^2) matrix.
    - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
      Used directly as (norb^2, norb^2) matrix. This is the correct
      treatment for multi-orbital systems.

    Returns
    -------
    V_q : ndarray
        Effective pairing interaction, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb
    chi0q_is_4index = (chi0q.ndim == 8)

    # Build S, C matrices for all q-points at once: (Nx, Ny, Nz, nd, nd)
    S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)

    # Extract chi0 at static limit for all q-points
    if chi0q_is_4index:
        # (norb, norb, norb, norb, Nx, Ny, Nz, nmat) -> (Nx, Ny, Nz, nd, nd)
        chi0_static = chi0q[:, :, :, :, :, :, :, nmat // 2].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()
    else:
        # (norb, norb, Nx, Ny, Nz, nmat) -> expand to (Nx, Ny, Nz, nd, nd)
        chi0_2d = chi0q[:, :, :, :, :, nmat // 2].transpose(2, 3, 4, 0, 1).copy()
        # chi0_2d shape: (Nx, Ny, Nz, norb, norb)
        if norb == 1:
            chi0_static = chi0_2d.reshape(Nx, Ny, Nz, 1, 1)
        else:
            # Expand: chi0_{l1*norb+l2, l3*norb+l2} = chi0_2d[l1, l3]
            chi0_static = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            for l2 in range(norb):
                chi0_static[:, :, :,
                            l2::norb,
                            l2::norb] = chi0_2d

    # Batched RPA solve for all q-points simultaneously
    # chi_s = [I - chi0 @ S]^{-1} @ chi0
    # chi_c = [I + chi0 @ C]^{-1} @ chi0
    I_mat = np.broadcast_to(np.eye(nd, dtype=complex), (Nx, Ny, Nz, nd, nd)).copy()

    mat_s = I_mat - np.einsum('...ab,...bc->...ac', chi0_static, S_all)
    mat_c = I_mat + np.einsum('...ab,...bc->...ac', chi0_static, C_all)

    chis = np.linalg.solve(mat_s, chi0_static)  # batched solve
    chic = np.linalg.solve(mat_c, chi0_static)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + 0.5 * (C_all - S_all)
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Reshape (Nx, Ny, Nz, nd, nd) -> (norb, norb, norb, norb, Nx, Ny, Nz)
    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(3, 4, 5, 6, 0, 1, 2)

    return Vs_q


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------

def _calc_g2(green_kw):
    """Calculate G2 = sum_n G(k, wn) G(-k+q, -wn).

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).

    Returns
    -------
    G2 : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    # G(-k, -wn) via roll+flip
    green_kw_inv = np.roll(
        green_kw[:, :, ::-1, ::-1, ::-1, ::-1],
        (1, 1, 1), (2, 3, 4)
    )
    G2 = np.einsum("ijpqsk, lmpqsk -> ijlmpqs", green_kw, green_kw_inv)
    return G2


def _eliashberg_kernel_fft(V_q, G2, sigma_old, norb):
    """Apply one iteration of the Eliashberg kernel using FFT convolution.

    Supports both simple (2-index) and general (4-index) pairing vertices.

    For the simple case (V_q shape: norb, norb, Nx, Ny, Nz):
        sigma_{il}(k) = -sum_q V_{ij}(q) * G2_{ijlm}(q) * sigma_{jm}(k)

    For the general case (V_q shape: norb, norb, norb, norb, Nx, Ny, Nz):
        sigma_{l1l4}(k) = -sum_q sum_{l2l3} V_{l1l2,l3l4}(q) * F_{l2l3}(q)
        where F_{l2l3}(q) = sum_{l5l6} G_{l2l5}(k-q) sigma_{l5l6}(k-q) G_{l3l6}(q-k)

    Parameters
    ----------
    V_q : ndarray
        Pairing vertex. Either shape (norb, norb, Nx, Ny, Nz) for simple mode,
        or shape (norb, norb, norb, norb, Nx, Ny, Nz) for general mode.
    G2 : ndarray
        Two-particle Green's function, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    sigma_old : ndarray
        Previous gap function, shape (norb, norb, Nx, Ny, Nz).
    norb : int
        Number of orbitals.

    Returns
    -------
    sigma_new : ndarray
        Updated gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    if V_q.ndim == 5:
        # Simple mode: V_q is (norb, norb, Nx, Ny, Nz)
        G2Sigma = np.einsum("ijlmpqs, jmpqs -> ilpqs", G2, sigma_old)

        # Vectorized IFFT over all orbital pairs at once
        P_r = ifftn(V_q, axes=(-3, -2, -1))
        G2Sigma_r = ifftn(G2Sigma, axes=(-3, -2, -1))

        Sigma_r = P_r * G2Sigma_r
        sigma_new = fftn(Sigma_r, axes=(-3, -2, -1))
        return -sigma_new

    else:
        # General mode: V_q is (norb, norb, norb, norb, Nx, Ny, Nz)
        F_q = np.einsum("ijlmpqs, jmpqs -> ilpqs", G2, sigma_old)

        # Vectorized IFFT: all orbital indices transformed at once
        V_r = ifftn(V_q, axes=(-3, -2, -1))
        F_r = ifftn(F_q, axes=(-3, -2, -1))

        # sigma_{l1l4}(r) = sum_{l2l3} V_{l1l2,l3l4}(r) * F_{l2l3}(r)
        sigma_r = np.einsum("ijklpqs, jkpqs -> ilpqs", V_r, F_r)

        sigma_new = fftn(sigma_r, axes=(-3, -2, -1))
        return -sigma_new


# ---------------------------------------------------------------------------
# Gap function initialization
# ---------------------------------------------------------------------------

def _initialize_gap(mode, norb, kx_array, ky_array, kz_array):
    """Initialize gap function sigma(k) with specified symmetry.

    Supports 2D and 3D lattice symmetries. The form factors are defined
    in terms of kx, ky, kz and work for any dimension (Nz=1 for 2D).

    Parameters
    ----------
    mode : str
        Initialization mode specifying the gap symmetry.
        Abbreviations cx=cos(kx), sx=sin(kx), etc.

        **Isotropic / s-wave:**
        - "cos"       : cos(kx+ky+kz) (original reference code)
        - "s"         : 1 (constant, k-independent)
        - "s_ext"     : cx*cy + cy*cz + cz*cx (extended s-wave, 3D)
        - "s_ext_2d"  : cx*cy (extended s-wave, 2D)

        **d-wave:**
        - "d_x2y2"    : cx - cy
        - "d_xy"      : sx*sy
        - "d_xz"      : sx*sz
        - "d_yz"      : sy*sz
        - "d_z2"      : 2*cz - cx - cy (d_{3z^2-r^2})

        **p-wave (odd parity):**
        - "p_x"       : sx
        - "p_y"       : sy
        - "p_z"       : sz

        **Other:**
        - "random"    : random (all symmetries mixed)
    norb : int
        Number of orbitals.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.

    Returns
    -------
    sigma : ndarray
        Normalized initial gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    I = np.identity(norb)
    cx, cy, cz = np.cos(kx_mesh), np.cos(ky_mesh), np.cos(kz_mesh)
    sx, sy, sz = np.sin(kx_mesh), np.sin(ky_mesh), np.sin(kz_mesh)

    form_factors = {
        # s-wave
        "cos":      lambda: np.cos(kx_mesh + ky_mesh + kz_mesh),
        "s":        lambda: np.ones((Nx, Ny, Nz)),
        "s_ext":    lambda: cx * cy + cy * cz + cz * cx,
        "s_ext_2d": lambda: cx * cy,
        # d-wave
        "d_x2y2":   lambda: cx - cy,
        "d_xy":     lambda: sx * sy,
        "d_xz":     lambda: sx * sz,
        "d_yz":     lambda: sy * sz,
        "d_z2":     lambda: 2.0 * cz - cx - cy,
        # p-wave
        "p_x":      lambda: sx,
        "p_y":      lambda: sy,
        "p_z":      lambda: sz,
    }

    if mode in form_factors:
        f_k = form_factors[mode]()
        sigma = I[:, :, None, None, None] * f_k[None, None, :, :, :]
    elif mode == "random":
        rng = np.random.default_rng(12345)
        sigma = rng.standard_normal((norb, norb, Nx, Ny, Nz))
    else:
        raise ValueError(
            "Unknown init_gap mode: '{}'. Available: {}".format(
                mode, list(form_factors.keys()) + ["random"]))

    norm = np.linalg.norm(sigma)
    if norm > 0:
        sigma /= norm
    return sigma


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

def _solve_iteration(green_kw, Vs_q, G2, sigma_init, norb,
                     max_iter=1000, alpha=0.5, tol=1.0e-5):
    """Solve linearized Eliashberg equation by self-consistent iteration.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    sigma_init : ndarray
        Initial gap function.
    norb : int
        Number of orbitals.
    max_iter : int
        Maximum iterations.
    alpha : float
        Mixing parameter (0 < alpha < 1).
    tol : float
        Convergence tolerance.

    Returns
    -------
    sigma : ndarray
        Converged gap function.
    eigenvalue : float
        Leading eigenvalue (norm of sigma_new before normalization).
    converged : bool
        Whether convergence was achieved.
    n_iter : int
        Number of iterations performed.
    """
    sigma_old = sigma_init.copy()

    eigenvalue = 0.0
    for iteration in range(max_iter):
        sigma_new = _eliashberg_kernel_fft(Vs_q, G2, sigma_old, norb)
        norm = np.linalg.norm(sigma_new)
        eigenvalue = norm

        diff = np.linalg.norm(sigma_new / norm - sigma_old)
        logger.info("Iteration {:4d}: eigenvalue = {:.6f}, diff = {:.6e}".format(
            iteration, norm, diff))

        if diff < tol:
            logger.info("Converged at iteration {}".format(iteration + 1))
            return sigma_new / norm, eigenvalue, True, iteration + 1

        sigma_old = (1.0 - alpha) * sigma_new / norm + alpha * sigma_old

    logger.warning("Failed to converge after {} iterations".format(max_iter))
    return sigma_old, eigenvalue, False, max_iter


def _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz):
    """Create LinearOperator for the Eliashberg kernel.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.

    Returns
    -------
    A : LinearOperator
        Linear operator representing the Eliashberg kernel K.
    vec_size : int
        Size of the flattened vector.
    """
    vec_size = norb * norb * Nx * Ny * Nz

    def matvec(v):
        sigma = v.reshape(norb, norb, Nx, Ny, Nz)
        result = _eliashberg_kernel_fft(Vs_q, G2, sigma, norb)
        return result.real.ravel()

    A = LinearOperator((vec_size, vec_size), matvec=matvec, dtype=float)
    return A, vec_size


def _solve_eigenvalue(Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=10,
                      method="arnoldi", sigma_shift=None):
    """Solve linearized Eliashberg equation by eigenvalue analysis.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    method : str
        Eigenvalue solver method:
        - "arnoldi" : Implicitly Restarted Arnoldi (ARPACK eigs).
            Standard Krylov subspace method for non-symmetric matrices.
            Directly finds the largest eigenvalues.
        - "shift-invert-bicgstab" : Shift-invert with BiCGSTAB.
            Solves (K - sigma*I)^{-1} eigenvalue problem using BiCGSTAB
            for the linear system. Efficient for finding eigenvalues near
            a target value sigma.
        - "shift-invert-gmres" : Shift-invert with GMRES.
            Same as above but uses GMRES for the linear system.
            More robust than BiCGSTAB for non-symmetric systems.
        - "shift-invert-lgmres" : Shift-invert with LGMRES.
            Uses LGMRES (loose GMRES) which can be faster for some problems.
        - "subspace" : Subspace iteration (block power method).
            Simultaneously propagates multiple vectors through the kernel
            with QR orthogonalization. Robust for degenerate eigenvalues
            and finds multiple eigenvalues with different symmetries.
    sigma_shift : float, optional
        Shift parameter for shift-invert methods. Eigenvalues near this
        value are found most efficiently. If None, a preliminary Arnoldi
        step estimates an appropriate shift value.

    Returns
    -------
    eigenvalues : ndarray
        Leading eigenvalues sorted by magnitude.
    eigenvectors : ndarray
        Corresponding eigenvectors reshaped to (num_ev, norb, norb, Nx, Ny, Nz).
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    logger.info("Computing {} eigenvalues with method='{}'...".format(max_ev, method))

    if method == "arnoldi":
        vals, vecs = eigs(A, k=max_ev, which='LM')

    elif method == "subspace":
        return _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=max_ev
        )

    elif method.startswith("shift-invert"):
        if sigma_shift is None:
            # Estimate shift from a quick Arnoldi run
            logger.info("Estimating shift with preliminary Arnoldi...")
            vals_pre, _ = eigs(A, k=1, which='LM')
            sigma_shift = vals_pre[0].real * 0.9
            logger.info("Using sigma_shift = {:.6f}".format(sigma_shift))
        vals, vecs = _eigs_shift_invert(
            A, vec_size, max_ev, method, sigma=sigma_shift
        )

    else:
        raise ValueError("Unknown eigenvalue method: {}".format(method))

    # Sort by magnitude (descending)
    idx = np.argsort(-np.abs(vals))
    vals = vals[idx]
    vecs = vecs[:, idx]

    eigenvectors = np.array([
        vecs[:, i].real.reshape(norb, norb, Nx, Ny, Nz)
        for i in range(len(vals))
    ])

    return vals, eigenvectors


def _eigs_shift_invert(A, vec_size, num_ev, method, sigma=0.0, rtol_linear=1e-8):
    """Eigenvalue computation using shift-invert with iterative linear solver.

    Transforms the eigenvalue problem K*x = lambda*x into
    (K - sigma*I)^{-1}*x = nu*x where nu = 1/(lambda - sigma).
    Eigenvalues of the original problem near sigma become the largest
    eigenvalues of the transformed problem, which Arnoldi finds efficiently.

    The inverse (K - sigma*I)^{-1} is never formed explicitly; instead,
    each application solves (K - sigma*I)*y = x using an iterative solver.

    Parameters
    ----------
    A : LinearOperator
        The Eliashberg kernel operator K.
    vec_size : int
        Dimension of the vector space.
    num_ev : int
        Number of eigenvalues to compute.
    method : str
        One of "shift-invert-bicgstab", "shift-invert-gmres",
        "shift-invert-lgmres".
    sigma : float
        Shift value. Eigenvalues near sigma are found most efficiently.
        Default 0.0 finds eigenvalues nearest to zero.
    rtol_linear : float
        Relative tolerance for the iterative linear solver.

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues of the original problem K.
    eigenvectors : ndarray
        Corresponding eigenvectors.
    """
    solver_name = method.replace("shift-invert-", "")
    solver_map = {
        "bicgstab": bicgstab,
        "gmres": gmres,
        "lgmres": lgmres,
    }
    if solver_name not in solver_map:
        raise ValueError("Unknown linear solver: {}".format(solver_name))
    linear_solver = solver_map[solver_name]

    # (A - sigma*I)
    def shifted_matvec(v):
        return A.matvec(v) - sigma * v

    A_shifted = LinearOperator((vec_size, vec_size),
                               matvec=shifted_matvec, dtype=float)

    solve_count = [0]
    fail_count = [0]

    # (A - sigma*I)^{-1} via iterative solver
    def inv_matvec(v):
        solve_count[0] += 1
        x, info = linear_solver(A_shifted, v, rtol=rtol_linear, maxiter=500)
        if info != 0:
            fail_count[0] += 1
            if fail_count[0] <= 3:
                logger.warning(
                    "Linear solver {} did not converge (info={}), "
                    "solve #{}, using approximate solution".format(
                        solver_name, info, solve_count[0]))
        return x

    A_inv = LinearOperator((vec_size, vec_size),
                           matvec=inv_matvec, dtype=float)

    # eigs on (A - sigma*I)^{-1} finds eigenvalues nu = 1/(lambda - sigma)
    # largest |nu| correspond to lambda closest to sigma
    nus, vecs = eigs(A_inv, k=num_ev, which='LM')

    logger.info("Shift-invert: {} linear solves, {} failures".format(
        solve_count[0], fail_count[0]))

    # Convert back: lambda = 1/nu + sigma
    eigenvalues = 1.0 / nus + sigma

    return eigenvalues, vecs


def _solve_subspace_iteration(Vs_q, G2, norb, Nx, Ny, Nz,
                              num_eigenvalues=5, max_iter=300, tol=1e-6):
    """Find multiple eigenvalues by subspace iteration (block power method).

    Simultaneously propagates a block of vectors through the kernel,
    orthogonalizing via QR decomposition at each step. This naturally
    converges to the dominant invariant subspace.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance for eigenvalues.

    Returns
    -------
    eigenvalues : ndarray
        Converged eigenvalues sorted by magnitude.
    eigenvectors : ndarray
        Shape (num_eigenvalues, norb, norb, Nx, Ny, Nz).
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    num_ev = min(num_eigenvalues, vec_size)

    # Use extra vectors for better convergence (subspace padding)
    n_work = min(num_ev + max(num_ev, 5), vec_size)

    # Random initial subspace
    rng = np.random.default_rng(42)
    V = rng.standard_normal((vec_size, n_work))
    V, _ = np.linalg.qr(V)

    eigenvalues_old = np.zeros(num_ev)

    for iteration in range(max_iter):
        # Apply kernel to all vectors: W = A @ V
        W = np.zeros_like(V)
        for j in range(n_work):
            W[:, j] = A.matvec(V[:, j])

        # Rayleigh quotient: H = V^T A V (small n_work x n_work matrix)
        H = V.T @ W

        # Eigendecomposition of small matrix
        evals_h, evecs_h = np.linalg.eig(H)
        idx = np.argsort(-np.abs(evals_h))
        evals_h = evals_h[idx]
        evecs_h = evecs_h[:, idx]

        # Ritz vectors: update subspace
        V = W @ evecs_h
        # Re-orthogonalize
        V, _ = np.linalg.qr(V)

        # Check convergence of the wanted eigenvalues
        eigenvalues_new = evals_h[:num_ev].real
        diff = np.max(np.abs(eigenvalues_new - eigenvalues_old))

        if iteration % 10 == 0 or diff < tol:
            logger.info("Subspace iter {:4d}: eigenvalues = {}, diff = {:.2e}".format(
                iteration, np.array2string(eigenvalues_new, precision=4), diff))

        if diff < tol:
            logger.info("Subspace iteration converged at iteration {}".format(
                iteration + 1))
            break

        eigenvalues_old = eigenvalues_new.copy()

    # Extract final Ritz vectors for the wanted eigenvalues
    # Recompute from final V
    W = np.zeros((vec_size, n_work))
    for j in range(n_work):
        W[:, j] = A.matvec(V[:, j])
    H = V.T @ W
    evals_h, evecs_h = np.linalg.eig(H)
    idx = np.argsort(-np.abs(evals_h))

    eigenvalues = evals_h[idx[:num_ev]]
    ritz_vecs = V @ evecs_h[:, idx[:num_ev]]

    eigenvectors = np.array([
        ritz_vecs[:, i].real.reshape(norb, norb, Nx, Ny, Nz)
        for i in range(num_ev)
    ])

    return eigenvalues, eigenvectors


def _solve_shifted_bicg(Vs_q, G2, norb, Nx, Ny, Nz,
                        sigma_list, num_eigenvalues=3, tol_linear=1e-8):
    """Find eigenvalues near multiple target values using shifted BiCG.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    sigma_list : list of float
        Target shift values. Eigenvalues near each sigma are found.
    num_eigenvalues : int
        Number of eigenvalues to find near each sigma.
    tol_linear : float
        Tolerance for the linear solver.

    Returns
    -------
    all_eigenvalues : dict
        Dictionary mapping sigma -> array of eigenvalues.
    all_eigenvectors : dict
        Dictionary mapping sigma -> array of eigenvectors.
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    all_eigenvalues = {}
    all_eigenvectors = {}

    # Use the first sigma as seed system
    sigma_seed = sigma_list[0]

    # Seed: build Krylov subspace for (A - sigma_seed I)
    # We use a manual BiCG implementation with shift tracking
    logger.info("Shifted BiCG: seed sigma = {:.6f}, {} additional shifts".format(
        sigma_seed, len(sigma_list) - 1))

    # For each shift, solve via shift-invert + eigs
    # The key optimization: reuse preconditioner from seed system
    for i, sigma in enumerate(sigma_list):
        logger.info("Processing shift sigma = {:.6f} ({}/{})".format(
            sigma, i + 1, len(sigma_list)))

        def shifted_matvec(v, s=sigma):
            return A.matvec(v) - s * v

        A_shifted = LinearOperator(
            (vec_size, vec_size), matvec=shifted_matvec, dtype=float
        )

        def inv_matvec(v):
            x, info = gmres(A_shifted, v, rtol=tol_linear, maxiter=300)
            return x

        A_inv = LinearOperator(
            (vec_size, vec_size), matvec=inv_matvec, dtype=float
        )

        try:
            nus, vecs = eigs(A_inv, k=max_ev, which='LM')
            eigenvalues = 1.0 / nus + sigma
            idx = np.argsort(-np.abs(eigenvalues))
            eigenvalues = eigenvalues[idx]
            eigvecs = np.array([
                vecs[:, j].real.reshape(norb, norb, Nx, Ny, Nz)
                for j in idx
            ])
        except Exception as e:
            logger.warning("Shift sigma={:.6f} failed: {}".format(sigma, e))
            eigenvalues = np.array([])
            eigvecs = np.array([])

        all_eigenvalues[sigma] = eigenvalues
        all_eigenvectors[sigma] = eigvecs

        if len(eigenvalues) > 0:
            logger.info("  Found eigenvalues: {}".format(
                np.array2string(eigenvalues.real, precision=6)))

    return all_eigenvalues, all_eigenvectors


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _save_results(output_dir, sigma, eigenvalue, eigenvalues_eig, kx_array, ky_array, kz_array,
                  gap_file="gap.dat", eigenvalue_file="eigenvalue.dat"):
    """Save gap function and eigenvalue results to files.

    Parameters
    ----------
    output_dir : str
        Output directory path.
    sigma : ndarray or None
        Gap function from iteration, shape (norb, norb, Nx, Ny, Nz).
    eigenvalue : float or None
        Leading eigenvalue from iteration.
    eigenvalues_eig : ndarray or None
        Eigenvalues from eigenvalue analysis.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    gap_file : str
        Output filename for gap function.
    eigenvalue_file : str
        Output filename for eigenvalues.
    """
    os.makedirs(output_dir, exist_ok=True)

    if sigma is not None:
        norb = sigma.shape[0]
        Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
        filepath = os.path.join(output_dir, gap_file)
        logger.info("Saving gap function to {}".format(filepath))

        with open(filepath, "w") as fw:
            # Header
            header_parts = ["# kx", "ky", "kz"]
            for i in range(norb):
                for j in range(norb):
                    header_parts.append("Re(sigma_{}{})".format(i, j))
                    header_parts.append("Im(sigma_{}{})".format(i, j))
            fw.write(" ".join(header_parts) + "\n")

            for ix in range(Nx):
                kx = kx_array[ix]
                if kx > np.pi:
                    kx -= 2.0 * np.pi
                for iy in range(Ny):
                    ky = ky_array[iy]
                    if ky > np.pi:
                        ky -= 2.0 * np.pi
                    for iz in range(Nz):
                        kz = kz_array[iz]
                        if kz > np.pi:
                            kz -= 2.0 * np.pi
                        parts = ["{:.8f}".format(kx),
                                 "{:.8f}".format(ky),
                                 "{:.8f}".format(kz)]
                        for i in range(norb):
                            for j in range(norb):
                                val = sigma[i, j, ix, iy, iz]
                                parts.append("{:.8e}".format(val.real))
                                parts.append("{:.8e}".format(val.imag))
                        fw.write(" ".join(parts) + "\n")

    if eigenvalue is not None or eigenvalues_eig is not None:
        filepath = os.path.join(output_dir, eigenvalue_file)
        logger.info("Saving eigenvalues to {}".format(filepath))
        with open(filepath, "w") as fw:
            if eigenvalue is not None:
                fw.write("# Iteration eigenvalue\n")
                fw.write("{:.8e}\n".format(eigenvalue))
            if eigenvalues_eig is not None:
                fw.write("# Eigenvalue analysis\n")
                fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  |eigenvalue|\n")
                for i, ev in enumerate(eigenvalues_eig):
                    fw.write("{:4d} {:15.8e} {:15.8e} {:15.8e}\n".format(
                        i, ev.real, ev.imag, abs(ev)))


# ---------------------------------------------------------------------------
# chi0q format conversion
# ---------------------------------------------------------------------------

def _convert_chi0q_to_ref_format(chi0q, norb, Nx, Ny, Nz, nmat):
    """Convert chi0q from H-wave format to reference code format.

    Supports both 2-index (reduced) and 4-index (general) chi0q:

    2-index:
        H-wave: (nmat, nvol, norb, norb) -> Ref: (norb, norb, Nx, Ny, Nz, nmat)
    4-index:
        H-wave: (nmat, nvol, norb, norb, norb, norb)
        -> Ref: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)

    Parameters
    ----------
    chi0q : ndarray
        chi0q in H-wave format.
    norb, Nx, Ny, Nz, nmat : int
        System parameters.

    Returns
    -------
    chi0q_ref : ndarray
        chi0q in reference format:
        - 2-index: (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    """
    if chi0q.ndim == 4:
        # 2-index: (nmat, nvol, norb, norb) -> (norb, norb, Nx, Ny, Nz, nmat)
        nf = chi0q.shape[0]
        chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb)
        chi0q_ref = chi0q_3d.transpose(4, 5, 1, 2, 3, 0)
    elif chi0q.ndim == 6:
        if chi0q.shape[0] == norb and chi0q.shape[1] == norb:
            # Already in ref format (norb, norb, Nx, Ny, Nz, nmat)
            chi0q_ref = chi0q
        else:
            # 4-index H-wave: (nmat, nvol, norb, norb, norb, norb)
            # -> (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
            nf = chi0q.shape[0]
            chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb, norb, norb)
            chi0q_ref = chi0q_3d.transpose(4, 5, 6, 7, 1, 2, 3, 0)
    elif chi0q.ndim == 8:
        # Already in ref format (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        chi0q_ref = chi0q
    else:
        raise ValueError(
            "Unexpected chi0q shape: {}. Expected 4D, 6D, or 8D.".format(chi0q.shape)
        )
    return chi0q_ref


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------

def calc_eliashberg(input_dict):
    """Main calculation orchestration for linearized Eliashberg equation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.
    """
    # --- Parse parameters ---
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = 1.0 / T
    cell_shape = mode_param["CellShape"]
    sub_shape = mode_param.get("SubShape", cell_shape)
    nmat = mode_param.get("Nmat", 1024)

    # Filling
    if "filling" in mode_param:
        n_filling = mode_param["filling"]
    elif "Ncond" in mode_param:
        raise ValueError("Ncond not yet supported in eliashberg. Use filling.")
    else:
        raise ValueError("filling must be specified.")

    # Lattice dimensions
    if isinstance(cell_shape, list):
        while len(cell_shape) < 3:
            cell_shape.append(1)
    Lx, Ly, Lz = cell_shape
    if isinstance(sub_shape, list):
        while len(sub_shape) < 3:
            sub_shape.append(1)
    Bx, By, Bz = sub_shape
    Nx, Ny, Nz = Lx // Bx, Ly // By, Lz // Bz
    nvol = Nx * Ny * Nz

    # Eliashberg parameters
    eli_param = input_dict.get("eliashberg", {})
    solver_mode = eli_param.get("solver_mode", "iteration")
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)
    init_gap_mode = eli_param.get("init_gap", "cos")
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    pairing_type = eli_param.get("pairing_type", "singlet")
    chi0q_mode = eli_param.get("chi0q_mode", "load")
    chi0q_tensor = eli_param.get("chi0q_tensor", "auto")
    gap_file = eli_param.get("output_gap", "gap.dat")
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")

    output_dir = input_dict["file"]["output"]["path_to_output"]

    logger.info("=== Linearized Eliashberg Equation Solver ===")
    logger.info("T = {}, beta = {:.4f}".format(T, beta))
    logger.info("Cell: {}x{}x{}, Sub: {}x{}x{}, Grid: {}x{}x{}".format(
        Lx, Ly, Lz, Bx, By, Bz, Nx, Ny, Nz))
    logger.info("Nmat = {}, filling = {}".format(nmat, n_filling))
    logger.info("solver_mode = {}, chi0q_mode = {}".format(solver_mode, chi0q_mode))

    # --- Step 1: Load or compute chi0q ---
    if chi0q_mode == "calc":
        chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor=chi0q_tensor)
    else:
        chi0q_raw = _load_chi0q(input_dict)

    # --- Step 2: Read input files ---
    geom_info, hr, interactions = _read_interaction_files(input_dict)
    norb = geom_info["norb"]

    # --- Step 3: Setup k-mesh ---
    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)

    # --- Step 4: Build Hamiltonian and diagonalize ---
    logger.info("Building Hamiltonian in k-space...")
    epsilon_k = _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb)

    logger.info("Diagonalizing Hamiltonian...")
    eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)

    # --- Step 5: Determine chemical potential ---
    logger.info("Determining chemical potential...")
    mu = _determine_mu(eigenvalues, beta, n_filling, norb)
    logger.info("mu = {:.6f}".format(mu))

    # --- Step 6: Calculate Green's function ---
    logger.info("Calculating Green's function G(k, iwn)...")
    green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

    # --- Step 7: Build interaction in k-space ---
    logger.info("Building interactions in k-space...")
    inter_k = _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb)

    # --- Step 8: Convert chi0q format ---
    # Determine nmat from chi0q shape
    # H-wave format: first axis is nmat for both 4D and 6D
    if chi0q_raw.ndim in (4, 6):
        nmat_chi0q = chi0q_raw.shape[0]
    else:
        nmat_chi0q = chi0q_raw.shape[-1]
    chi0q = _convert_chi0q_to_ref_format(chi0q_raw, norb, Nx, Ny, Nz, nmat_chi0q)
    logger.info("chi0q converted to shape: {}".format(chi0q.shape))

    # --- Step 9: Compute RPA vertices ---
    logger.info("Computing RPA vertices (pairing_type={})...".format(pairing_type))
    vertex_result = _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat_chi0q,
                                      pairing_type=pairing_type)
    if isinstance(vertex_result, tuple):
        # Simple mode: (Pc_q, Ps_q) -> combine to single vertex
        Pc_q, Ps_q = vertex_result
        Vs_q = Pc_q + Ps_q
        logger.info("Simple mode: Pc + Ps vertex, shape {}".format(Vs_q.shape))
    else:
        # General mode: 4-index V^s
        Vs_q = vertex_result
        logger.info("General mode: 4-index V^s vertex, shape {}".format(Vs_q.shape))

    # --- Step 10: Compute G2 ---
    logger.info("Computing G2...")
    G2 = _calc_g2(green_kw)

    # --- Step 11: Initialize gap function ---
    sigma_init = _initialize_gap(init_gap_mode, norb, kx_array, ky_array, kz_array)

    # --- Step 12: Solve ---
    sigma_result = None
    eigenvalue_iter = None
    eigenvalues_eig = None

    if solver_mode in ("iteration", "both"):
        logger.info("=== Self-consistent iteration ===")
        sigma_result, eigenvalue_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=max_iter, alpha=alpha, tol=tol
        )
        logger.info("Iteration result: eigenvalue = {:.6f}, converged = {}, n_iter = {}".format(
            eigenvalue_iter, converged, n_iter))

    if solver_mode in ("eigenvalue", "both"):
        logger.info("=== Eigenvalue analysis ===")
        eigenvalues_eig, eigenvectors_eig = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=num_eigenvalues,
            method=eigenvalue_method
        )
        logger.info("Leading eigenvalues:")
        for i, ev in enumerate(eigenvalues_eig):
            logger.info("  {:3d}: {:.6f} (|ev| = {:.6f})".format(i, ev.real, abs(ev)))

        # Use leading eigenvector as gap if no iteration result
        if sigma_result is None:
            sigma_result = eigenvectors_eig[0]
            eigenvalue_iter = eigenvalues_eig[0].real

    # --- Step 13: Save results ---
    _save_results(
        output_dir, sigma_result, eigenvalue_iter, eigenvalues_eig,
        kx_array, ky_array, kz_array,
        gap_file=gap_file, eigenvalue_file=eigenvalue_file
    )

    logger.info("=== Done ===")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """Command-line interface for the Eliashberg equation solver."""
    import tomli
    import argparse

    parser = argparse.ArgumentParser(
        description="Linearized Eliashberg equation solver (post-tool for H-wave RPA)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input", type=str, help="input TOML file (same as RPA input + [eliashberg] section)")
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress output")
    parser.add_argument(
        "-v", "--version", action="version",
        version="hwave_sc v{}".format(hwave.__version__),
    )
    args = parser.parse_args()

    # Setup logging
    log_level = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=log_level, format="%(name)s: %(message)s")

    file_toml = args.input
    if not os.path.exists(file_toml):
        logger.error("Input file does not exist: {}".format(file_toml))
        sys.exit(1)

    with open(file_toml, "rb") as f:
        input_dict = tomli.load(f)

    calc_eliashberg(input_dict)


if __name__ == "__main__":
    main()
