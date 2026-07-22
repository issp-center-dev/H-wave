"""Sublattice unfold helpers for SubShape > [1, 1, 1].

Physical site coordinates come from ``geometry_uhf.dat``. The decode
functions consume ``(cx, cy, cz)`` integer lattice coords directly,
never a row-major site ordinal.

Encoding matches ``_apbc_phase.sublattice_offset`` inverse:
    folded_orb = orig_orb + norb_orig * (sx + Sx * (sy + Sy * sz))

For the down-spin row offset in H-wave's non-spin-orbital
mode, ``nd = 2 * norb_folded`` and down rows begin at ``norb_folded``.
``folded_row_indices`` exposes both row indices so the bridge never
addresses the up block when it means the down block. See
docs/en/source/algorithm/uhfk_to_mvmc.rst.
"""
from __future__ import annotations

import numpy as np


def decode_physical_site(
    site_coord: np.ndarray, subshape: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return (folded_cell, sub_offset) as (3,) int arrays.

    ``site_coord`` is a (3,) integer lattice-unit coordinate from
    ``geometry_uhf.dat``. Never derive coordinates from a site ordinal.
    """
    site_coord = np.asarray(site_coord, dtype=np.int64)
    subshape = np.asarray(subshape, dtype=np.int64)
    folded_cell = site_coord // subshape
    sub_offset = site_coord % subshape
    return folded_cell, sub_offset


def encode_folded_orbital(
    orig_orb: int,
    sub_offset: np.ndarray,
    norb_orig: int,
    subshape: np.ndarray,
) -> int:
    """Encode ``(orig_orb, sub_offset)`` as ``folded_orb``.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    sub_offset = np.asarray(sub_offset, dtype=np.int64)
    subshape = np.asarray(subshape, dtype=np.int64)
    sx, sy, sz = int(sub_offset[0]), int(sub_offset[1]), int(sub_offset[2])
    Sx, Sy, _ = int(subshape[0]), int(subshape[1]), int(subshape[2])
    return int(orig_orb) + int(norb_orig) * (sx + Sx * (sy + Sy * sz))


def folded_row_indices(folded_orb: int, norb_folded: int) -> tuple[int, int]:
    """Return (row_up, row_down) in the folded nd index."""
    return int(folded_orb), int(norb_folded) + int(folded_orb)


def unfold_amplitude_columns(
    folded_wavevector_index: np.ndarray,   # (nvol_folded, 3)
    cell_shape: np.ndarray,
    subshape: np.ndarray,
    site_positions: np.ndarray,            # (Ns_phys, 3) integer coords
    norb_orig: int,
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (plane_wave_up, plane_wave_down, row_up_per_site,
    row_down_per_site).

    plane_wave_up[k, i]   = exp(+i k_folded · folded_cell_i) * exp(+i theta · r_i / L_phys) / sqrt(nvol_folded)
    plane_wave_down[k, i] = plane_wave_up[k, i]  (partner-row lookup lives in the caller)
    row_up_per_site[i]   = folded_orb_i
    row_down_per_site[i] = norb_folded + folded_orb_i

    Sign convention
    ---------------
    H-wave's folded eigenvector ``v[k_folded, aa, l]`` is diagonalized
    under the negative-gauge APBC transformation ``tilde_c_r =
    exp(-i theta r / L_phys) c_r`` (see :file:`_apbc_phase.py`), and its
    tilde-side annihilation Fourier convention is POSITIVE Bloch
    (``c_R = (1/sqrt N) sum_k c_k exp(+i k R)``, matching H-wave's
    ``ifftn(..., norm='forward')`` on the Hamiltonian). The physical
    Bloch amplitude at site ``r`` for the folded band ``(k_folded, l)``
    is therefore

        psi_phys(r) = (1/sqrt(N_folded))
                    * v[k_folded, sub_offset(r), l]
                    * exp(+i k_folded . folded_cell(r))
                    * exp(+i theta . r / L_phys).

    For ``SubShape = [1, 1, 1]`` the (k, -k) time-reversal pair sum
    symmetrises the plane-wave factor, so both signs produce the same
    density. See docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    folded_wavevector_index = np.asarray(folded_wavevector_index, dtype=np.int64)
    cell_shape = np.asarray(cell_shape, dtype=np.int64)
    subshape = np.asarray(subshape, dtype=np.int64)
    site_positions = np.asarray(site_positions, dtype=np.int64)
    theta = np.asarray(theta, dtype=np.float64)

    nvol_folded = folded_wavevector_index.shape[0]
    Ns_phys = site_positions.shape[0]
    subvol = int(np.prod(subshape))
    norb_folded = int(norb_orig) * subvol

    L_folded = cell_shape // subshape
    L_phys = cell_shape.astype(np.float64)

    row_up = np.empty(Ns_phys, dtype=np.int64)
    row_down = np.empty(Ns_phys, dtype=np.int64)
    folded_cell = np.empty((Ns_phys, 3), dtype=np.int64)
    for i in range(Ns_phys):
        fc, so = decode_physical_site(site_positions[i], subshape)
        folded_cell[i] = fc
        folded_orb = encode_folded_orbital(0, so, norb_orig, subshape)
        row_up[i], row_down[i] = folded_row_indices(folded_orb, norb_folded)

    # k_folded_phys per row: 2 pi * n / L_folded, shape (nvol_folded, 3)
    k_folded = 2.0 * np.pi * folded_wavevector_index.astype(np.float64) / L_folded

    # plane_wave_up[k, i] = exp(+i k_folded · folded_cell_i) / sqrt(nvol_folded)
    #   * exp(+i theta · r_i / L_phys)
    inv_sqrt = 1.0 / np.sqrt(float(nvol_folded))
    kf_dot_fc = np.einsum("kd,id->ki", k_folded, folded_cell.astype(np.float64))
    plane_wave_up_folded = np.exp(+1j * kf_dot_fc) * inv_sqrt

    theta_over_L = theta / L_phys
    phys_arg = np.einsum("d,id->i", theta_over_L, site_positions.astype(np.float64))
    phys_up = np.exp(+1j * phys_arg)   # exp(+i theta r / L)

    # Down uses the SAME plane-wave envelope as up. The (k, -k)
    # time-reversal partnership is resolved via the caller's
    # ``eigenvector[partner_n, row_down_per_site[i], col_down]`` lookup
    # (fij_builder.build_amplitudes), NOT via a separately-conjugated
    # plane_wave_down. This unified form reproduces H-wave's
    # ``greenone.dat`` on both spins; see
    # docs/en/source/algorithm/uhfk_to_mvmc.rst.
    plane_wave_up = plane_wave_up_folded * phys_up[np.newaxis, :]
    plane_wave_down = plane_wave_up_folded * phys_up[np.newaxis, :]

    return plane_wave_up, plane_wave_down, row_up, row_down
