"""MYO-convention spin/charge interaction matrices for full-vertex FLEX.

This module provides :func:`build_sc_matrices_myo`, a sibling of
``hwave.sc._build_sc_matrices_all_q`` that builds the spin (S) and charge (C)
interaction matrices in the **MYO (Mochizuki-Yanase-Ogata)** convention,
following cond-mat/0407094 Eq.(6), rather than the Kuroki convention used in
``hwave.sc`` (arXiv:0902.3691 Eq.(5)).

The two conventions differ in EXACTLY ONE matrix element: the charge (ab,ab)
element (Case 2 below). MYO uses ``-U'+2J`` while Kuroki (``sc.py``) uses
``-U'+J``. All spin elements and all other charge elements are identical.

The structure of this builder is copied verbatim from
``hwave.sc._build_sc_matrices_all_q`` so that the two stay element-for-element
comparable; only the single Case 2 charge term is modified.
"""

import numpy as np


def build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz):
    """Build MYO-convention spin (S) and charge (C) matrices for all q-points.

    Follows MYO (Mochizuki-Yanase-Ogata), cond-mat/0407094 Eq.(6):
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k. Each entry is an
        ndarray of shape (norb, norb, Nx, Ny, Nz).
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

    # Build using precomputed index arrays to avoid Python loops
    # for small norb (1-3), the loop overhead is negligible;
    # for larger norb the vectorized approach helps
    l1_arr, l2_arr, l3_arr, l4_arr = np.meshgrid(
        np.arange(norb), np.arange(norb), np.arange(norb), np.arange(norb),
        indexing='ij')
    l1f, l2f, l3f, l4f = l1_arr.ravel(), l2_arr.ravel(), l3_arr.ravel(), l4_arr.ravel()
    idx12 = l1f * norb + l2f
    idx34 = l3f * norb + l4f

    # Case 1: l1 == l2 == l3 == l4
    # Accumulate rather than assign: this element is also reached by Case 3
    # below (which now includes l1 == l3), where an inter-site same-orbital
    # CoulombInter contributes 2 V_aa(q) to the charge channel.
    mask1 = (l1f == l2f) & (l2f == l3f) & (l3f == l4f)
    if U_mat is not None and np.any(mask1):
        for i in np.where(mask1)[0]:
            _l = l1f[i]
            S_all[:, :, :, idx12[i], idx34[i]] += U_mat[_l, _l]
            C_all[:, :, :, idx12[i], idx34[i]] += U_mat[_l, _l]

    # Case 2: l1==l3, l2==l4, l1!=l2
    mask2 = (l1f == l3f) & (l2f == l4f) & (l1f != l2f)
    for i in np.where(mask2)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        if Up_mat is not None:
            s_q += Up_mat[_l1, _l2]
            c_q -= Up_mat[_l1, _l2]
        if I_mat is not None:
            s_q -= I_mat[_l1, _l2]
            c_q -= I_mat[_l1, _l2]
        if J_mat is not None:
            # MYO cond-mat/0407094 Eq.(6): charge (ab,ab) = -U'+2J
            # (Kuroki sc.py uses -U'+J, i.e. coefficient 1 here).
            c_q += 2.0 * J_mat[_l1, _l2]
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 3: l1==l2, l3==l4 -- INCLUDING l1 == l3, but for CoulombInter ONLY.
    # The l1 != l3 exclusion dropped the inter-site same-orbital CoulombInter
    # from the charge diagonal C[(a,a),(a,a)], which must be U_a + 2 V_aa(q)
    # (issue #95); the simple two-index formulation used by chi0q_mode="load"
    # builds exactly that (`Wc = U_k + 2 V_k`, _compute_vertices_simple).
    # Case 1 above writes U_a into the same element, so both accumulate.
    # Hund and Ising stay restricted to l1 != l3: an orbital has no Hund or
    # Ising coupling with itself, and letting a stray diagonal entry through
    # here would silently move S as well.
    mask3 = (l1f == l2f) & (l3f == l4f)
    for i in np.where(mask3)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l3 = l1f[i], l3f[i]
        if _l1 != _l3:
            if J_mat is not None:
                s_q += J_mat[_l1, _l3]
                c_q -= J_mat[_l1, _l3]
            if I_mat is not None:
                s_q -= 2.0 * I_mat[_l1, _l3]
        if Up_mat is not None:
            c_q += 2.0 * Up_mat[_l1, _l3]
        S_all[:, :, :, idx12[i], idx34[i]] += s_q
        C_all[:, :, :, idx12[i], idx34[i]] += c_q

    # Case 4: l1==l4, l2==l3, l1!=l2
    mask4 = (l1f == l4f) & (l2f == l3f) & (l1f != l2f)
    for i in np.where(mask4)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        if Jp_mat is not None:
            s_q += Jp_mat[_l1, _l2]
        if PH_mat is not None:
            s_q += PH_mat[_l1, _l2]
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = s_q  # S = C for this channel

    return S_all, C_all
