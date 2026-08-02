"""(k, -k) time-reversal pair Fij builder.

See docs/en/source/algorithm/uhfk_to_mvmc.rst for the Bloch convention,
sublattice unfolding, and canonical-pair construction.

Physical-basis amplitude convention (positive-Bloch, matches H-wave's
folded-eigenvector convention: negative-gauge APBC transform
``tilde_c_r = exp(-i theta r / L_phys) c_r`` combined with positive-
Bloch tilde-annihilation Fourier ``c_R = (1/sqrt N) sum_k c_k exp(+i
k R)`` — the latter follows from H-wave's ``ifftn(..., norm='forward')``
on the Hamiltonian). The physical Bloch amplitude for the folded band
``(k_folded, l)`` at site ``r_i`` is

    psi_phys(r_i, k_folded, l) = (1/sqrt(N_folded))
                              * v[k_folded, sub_offset(r_i), l]
                              * exp(+i k_folded . folded_cell(r_i))
                              * exp(+i theta . r_i / L_phys).

For the (k, -k) time-reversal partner the down side reads the same
plane-wave envelope AT THE PARTNER row and the partner's eigenvector
column:

    A_up_{i, alpha(k, n)}   = plane_wave[k, i]         v[k,       row_up_i, n_up]
    A_down_{j, alpha(k, n)} = plane_wave[k_p, j]       v[k_p,   row_down_j, n_down]

with plane_wave[k, r] = (1/sqrt(N_folded)) exp(+i k folded_cell(r))
exp(+i theta r/L_phys) and ``k_p`` = partner row (``-k`` residue
under APBC twist). Under this convention F_phys[i, j] = sum_alpha
A_up[i, alpha] A_down[j, alpha] reproduces H-wave's ``greenone.dat``
element-wise, both spins, for arbitrary SubShape (verified to 1e-14
for the SubShape=[2,1,1] APBC L=8 fixture).

"""
from __future__ import annotations

import numpy as np

from .partner_index import find_partner_rows


def build_amplitudes(
    wavevector_index,
    eigenvector,
    stepped_occupation,
    column_spin,
    column_mu_group,
    site_positions,
    norb_orig,
    theta,
    L,
    *,
    cell_shape=None,
    subshape=None,
):
    """Build the physical-basis A^up, A^down occupied-orbital matrices.

    Returns
    -------
    A_up, A_down : (Ns_phys, N_occ) complex arrays
        Each column is one occupied pair alpha(k, n) built from the
        physical-basis eigenstates. ``A_up @ A_down.T`` yields F_phys.

    The down orbital column index aligns with the up column index via the
    (k, -k) partner lookup; for self-pair k the partner equals the up's
    own row.

    Sublattice support
    ------------------
    When ``cell_shape`` and ``subshape`` are supplied and ``SubShape >
    [1, 1, 1]``, H-wave's eigenvectors live on the folded BZ of size
    ``nvol_folded = prod(cell_shape // subshape)``. Each physical site
    maps to a ``(folded_cell, sub_offset)`` pair via
    ``sublattice_unfold.decode_physical_site``, and the up / down row
    within the folded ``nd = 2 * norb_folded`` block depends on the
    folded orbital. Passing ``cell_shape=None`` selects
    ``SubShape = [1, 1, 1]``.

    Raises
    ------
    ValueError
        If the partner-mapped (k, -k) occupation set is not consistent
        with the actual occupied down set (in magnetic
        UHF or spin-dependent fillings, the up-and-down occupied sets
        need not coincide through the time-reversal partner map).
    ValueError
        If ``subshape`` does not divide ``cell_shape`` in every direction
        before amplitude construction.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    wavevector_index = np.asarray(wavevector_index, dtype=np.int64)
    eigenvector = np.asarray(eigenvector, dtype=np.complex128)
    stepped_occupation = np.asarray(stepped_occupation, dtype=np.float64)
    column_spin = np.asarray(column_spin, dtype=np.int64)
    column_mu_group = np.asarray(column_mu_group, dtype=np.int64)
    site_positions = np.asarray(site_positions, dtype=np.float64)
    theta = np.asarray(theta, dtype=np.float64)
    L = np.asarray(L, dtype=np.int64)

    if cell_shape is None:
        # Single-orbital fallback: SubShape=[1,1,1], L is already CellShape.
        cell_shape = np.asarray(L, dtype=np.int64)
        subshape = np.array([1, 1, 1], dtype=np.int64)
    else:
        cell_shape = np.asarray(cell_shape, dtype=np.int64)
        subshape = np.asarray(subshape, dtype=np.int64)
    if not np.all(cell_shape % subshape == 0):
        raise ValueError(
            f"SubShape {subshape.tolist()} does not divide "
            f"CellShape {cell_shape.tolist()} in every direction"
        )
    subvol = int(np.prod(subshape))
    norb_folded = int(norb_orig) * subvol

    nvol, nd = stepped_occupation.shape
    Ns_phys = site_positions.shape[0]
    if norb_orig != 1:
        raise NotImplementedError(
            "This conversion supports only norb_orig == 1; got "
            f"{norb_orig}. See "
            "docs/en/source/uhfk/tools/uhfk_to_mvmc.rst."
        )

    partner_rows, is_self_pair = find_partner_rows(wavevector_index, theta, L)

    from .sublattice_unfold import unfold_amplitude_columns

    plane_wave_up, plane_wave_down, row_up_per_site, row_down_per_site = (
        unfold_amplitude_columns(
            folded_wavevector_index=wavevector_index,
            cell_shape=cell_shape,
            subshape=subshape,
            site_positions=site_positions.astype(np.int64),
            norb_orig=norb_orig,
            theta=theta,
        )
    )

    up_cols = np.where(column_spin == 0)[0]
    down_cols = np.where(column_spin == 1)[0]
    if len(up_cols) == 0 or len(down_cols) == 0:
        raise ValueError(
            "Sz-fixed mode requires both up-only and down-only column "
            "blocks; got column_spin = " + str(column_spin.tolist())
        )

    # Enforce pair closure over (k_row, local_band); see
    # docs/en/source/algorithm/uhfk_to_mvmc.rst. local_band is the
    # position of a column within its spin block. In single-orbital cases
    # norb_folded == 1 so local_band is always 0
    # and the tuple check reduced to the row-only check.
    up_cols_list = list(up_cols)
    down_cols_list = list(down_cols)
    up_col_to_local = {int(c): idx for idx, c in enumerate(up_cols_list)}
    down_col_to_local = {int(c): idx for idx, c in enumerate(down_cols_list)}

    occ_up_kb = set()
    occ_down_kb = set()
    for n_row in range(nvol):
        for col in up_cols_list:
            if stepped_occupation[n_row, col] >= 0.5:
                occ_up_kb.add((n_row, up_col_to_local[int(col)]))
        for col in down_cols_list:
            if stepped_occupation[n_row, col] >= 0.5:
                occ_down_kb.add((n_row, down_col_to_local[int(col)]))
    partner_of_occ_up = {
        (int(partner_rows[n_row]), local_band)
        for (n_row, local_band) in occ_up_kb
    }
    if partner_of_occ_up != occ_down_kb:
        missing = sorted(partner_of_occ_up - occ_down_kb)[:10]
        extra = sorted(occ_down_kb - partner_of_occ_up)[:10]
        raise ValueError(
            "(k, -k) pair-closure violated over (k_row, local_band): "
            f"partner(occ_up) has {len(partner_of_occ_up)} entries, "
            f"occ_down has {len(occ_down_kb)}; missing in down: {missing}; "
            f"extra in down: {extra}. Sublattice + Sz-fixed requires the "
            "time-reversal partner of every occupied up local band to "
            "coincide with the same local band on the down side."
        )

    A_up_list = []
    A_down_list = []
    for n_row in range(nvol):
        for col_up_idx, col_up in enumerate(up_cols_list):
            if stepped_occupation[n_row, col_up] < 0.5:
                continue
            partner_n = int(partner_rows[n_row])
            # Pair up_col_local with down_col at the SAME local_band index;
            # see docs/en/source/algorithm/uhfk_to_mvmc.rst.
            col_down = down_cols_list[col_up_idx]
            a_up = np.empty(Ns_phys, dtype=np.complex128)
            a_down = np.empty(Ns_phys, dtype=np.complex128)
            # A_down uses plane_wave AT THE PARTNER row (not ``n_row``).
            # Under the documented positive-Bloch convention,
            # ``plane_wave_up[k, i] = v[k, s(i), l] * exp(+i k_folded
            # R(i)) * exp(+i theta r/L) / sqrt(N_folded)`` is exactly
            # the physical Bloch amplitude at row ``k``, so the down
            # partner at ``(-k)`` uses ``plane_wave_up[partner_n, i]``
            # combined with the partner-row eigenvector column.
            for i in range(Ns_phys):
                a_up[i] = plane_wave_up[n_row, i] * eigenvector[
                    n_row, row_up_per_site[i], col_up
                ]
                a_down[i] = plane_wave_down[partner_n, i] * eigenvector[
                    partner_n, row_down_per_site[i], col_down
                ]
            A_up_list.append(a_up)
            A_down_list.append(a_down)

    if not A_up_list:
        raise ValueError("no occupied up-spin states found")

    A_up = np.stack(A_up_list, axis=1)
    A_down = np.stack(A_down_list, axis=1)
    return A_up, A_down


def build_fij_phys(A_up, A_down):
    """Return F^phys_{ij} = (A_up @ A_down.T)_{ij}, shape (Ns, Ns) complex.

    F is built from ``c^dag_i↑ c^dag_j↓`` pair coefficients under the
    positive-Bloch construction (module docstring). For a
    translation-invariant Slater state on the unfolded lattice the
    result depends only on the physical displacement ``r_j - r_i``.
    See docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    A_up = np.asarray(A_up, dtype=np.complex128)
    A_down = np.asarray(A_down, dtype=np.complex128)
    return A_up @ A_down.T
