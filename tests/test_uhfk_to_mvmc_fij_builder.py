"""Tests for ``fij_builder``.

See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

from tools._uhfk_to_mvmc.fij_builder import (
    build_amplitudes,
    build_fij_phys,
)


def _klist(n):
    return np.roll(np.arange(n) - (n // 2), -(n // 2))


def test_build_amplitudes_pbc_l4_single_orbital_freeparticle():
    """Free particles, PBC L=4, single orbital, half filling (Ne_up = 2).
    Expect A_up to be 2 orthonormal plane waves at occupied k's."""
    L = 4
    norb = 1
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L)], dtype=np.float64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)
    k_vals = 2.0 * np.pi * wv[:, 0] / L
    ws = -2.0 * np.cos(k_vals)
    eigenvalue = np.stack([ws, ws], axis=1)
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0

    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)
    pytest.skip(
        "L=4 half filling is degenerate; covered by quarter-filling test "
        "in Step 5"
    )


def test_build_amplitudes_pbc_l8_quarter_filling_freeparticle():
    """Free particles, PBC L=8, single orbital, Ne=1 (only k=0 occupied)."""
    L = 8
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L)], dtype=np.float64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)

    k_vals = 2.0 * np.pi * wv[:, 0] / L
    ws = -2.0 * np.cos(k_vals)
    eigenvalue = np.stack([ws, ws], axis=1)
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)
    stepped_occupation = np.zeros((L, 2), dtype=np.float64)
    n_k0 = list(wv[:, 0]).index(0)
    stepped_occupation[n_k0, 0] = 1.0
    stepped_occupation[n_k0, 1] = 1.0

    A_up, A_down = build_amplitudes(
        wavevector_index=wv,
        eigenvector=eigenvector,
        stepped_occupation=stepped_occupation,
        column_spin=column_spin,
        column_mu_group=column_mu_group,
        site_positions=site_positions,
        norb_orig=1,
        theta=theta,
        L=L_vec,
    )

    assert A_up.shape == (L, 1)
    assert A_down.shape == (L, 1)
    expected = np.full(L, 1.0 / np.sqrt(L), dtype=np.complex128)
    np.testing.assert_allclose(A_up[:, 0], expected, atol=1e-12)
    np.testing.assert_allclose(A_down[:, 0], expected, atol=1e-12)


def test_build_fij_phys_pbc_l8_k0_only():
    """For Ne=1 (k=0 only) PBC, F_{ij} = 1/L for all (i, j)."""
    L = 8
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L)], dtype=np.float64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
    stepped_occupation = np.zeros((L, 2), dtype=np.float64)
    n_k0 = list(wv[:, 0]).index(0)
    stepped_occupation[n_k0, 0] = 1.0
    stepped_occupation[n_k0, 1] = 1.0
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)

    A_up, A_down = build_amplitudes(
        wavevector_index=wv, eigenvector=eigenvector,
        stepped_occupation=stepped_occupation, column_spin=column_spin,
        column_mu_group=column_mu_group,
        site_positions=site_positions, norb_orig=1, theta=theta, L=L_vec,
    )
    F = build_fij_phys(A_up, A_down)

    np.testing.assert_allclose(
        F, np.full((L, L), 1.0 / L, dtype=np.complex128), atol=1e-12
    )


def test_build_amplitudes_rejects_open_pair_closure():
    """If the up occupied set's (k, -k) partner does not equal the down
    occupied set, build_amplitudes must raise because
    magnetic / spin-dependent occupations cannot be silently encoded by
    the (k, -k) construction)."""
    import pytest

    L = 8
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L)], dtype=np.float64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)

    # Up occupies tilde_k = 0 (row 0); partner under PBC is row 0 itself
    # (self-pair). Down occupies tilde_k = 1 (row 1) instead of the
    # self-pair: this breaks (k, -k) closure.
    stepped_occupation = np.zeros((L, 2), dtype=np.float64)
    n_k0 = list(wv[:, 0]).index(0)
    n_k1 = list(wv[:, 0]).index(1)
    stepped_occupation[n_k0, 0] = 1.0
    stepped_occupation[n_k1, 1] = 1.0

    with pytest.raises(ValueError, match="pair-closure"):
        build_amplitudes(
            wavevector_index=wv, eigenvector=eigenvector,
            stepped_occupation=stepped_occupation, column_spin=column_spin,
            column_mu_group=column_mu_group,
            site_positions=site_positions, norb_orig=1,
            theta=theta, L=L_vec,
        )


def test_build_amplitudes_subshape_2_free_particle():
    """L=4 SubShape=[2,1,1] PBC free particle (v=1). Occupy k=0 only in
    folded BZ (nvol_folded=2). Verify A_up shape = (Ns_phys=4, 1)."""
    L_phys = 4
    subshape = np.array([2, 1, 1], dtype=np.int64)
    cell_shape = np.array([L_phys, 1, 1], dtype=np.int64)
    L_folded = cell_shape // subshape
    wv = np.array([[v, 0, 0] for v in _klist(L_folded[0])], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L_phys)], dtype=np.int64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    # nd = 2 * norb_folded = 2 * 2 = 4. Up rows [0,1], down rows [2,3].
    # Single k=0 occupied on each spin block, band 0.
    eigenvector = np.zeros((2, 4, 4), dtype=np.complex128)
    for k in range(2):
        # Diagonal identity per block (v=1 for the single occupied band)
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
        eigenvector[k, 2, 2] = 1.0
        eigenvector[k, 3, 3] = 1.0
    stepped_occupation = np.zeros((2, 4), dtype=np.float64)
    n_k0 = list(wv[:, 0]).index(0)
    stepped_occupation[n_k0, 0] = 1.0  # occupied up band 0
    stepped_occupation[n_k0, 2] = 1.0  # occupied down band 0
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    column_mu_group = np.array([0, 0, 1, 1], dtype=np.int64)

    A_up, A_down = build_amplitudes(
        wavevector_index=wv,
        eigenvector=eigenvector,
        stepped_occupation=stepped_occupation,
        column_spin=column_spin,
        column_mu_group=column_mu_group,
        site_positions=site_positions,
        norb_orig=1,
        theta=theta,
        L=L_folded,
        cell_shape=cell_shape,
        subshape=subshape,
    )

    # 1 up occupied, 1 down partner -> single alpha column
    assert A_up.shape == (L_phys, 1)
    assert A_down.shape == (L_phys, 1)
    # site 0: folded_cell=0, sub_offset=0, folded_orb=0 -> row_up=0
    #   A_up[0, 0] = 1/sqrt(nvol_folded=2) * eigenvector[0, 0, 0] = 1/sqrt(2)
    assert A_up[0, 0] == pytest.approx(1.0 / np.sqrt(2), abs=1e-12)
    # site 1: folded_cell=0, sub_offset=1, folded_orb=1 -> row_up=1
    #   A_up[1, 0] = 1/sqrt(2) * eigenvector[0, 1, 0]. eigenvector[0, 1, 0]=0
    #   (column 0 is up band 0, which lives on row 0). So A_up[1, 0] = 0.
    assert A_up[1, 0] == pytest.approx(0.0 + 0.0j, abs=1e-12)


def test_pair_closure_rejects_local_band_mismatch():
    """SubShape=[2,1,1] with 2 folded orbitals per spin block: if the
    up occupied set on band 0 has partner rows also occupied on band 0,
    but band 1 down side has a different occupation, row-only pair
    closure would pass. Must fail-fast on (k_row, local_band) tuple."""
    L_phys = 4
    subshape = np.array([2, 1, 1], dtype=np.int64)
    cell_shape = np.array([L_phys, 1, 1], dtype=np.int64)
    L_folded = cell_shape // subshape
    wv = np.array([[v, 0, 0] for v in _klist(L_folded[0])], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L_phys)], dtype=np.int64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    # nd=4 (2 up bands, 2 down bands); occupy up band 0 at k=0 AND k=1,
    # but only down band 1 at k=0 AND k=1. Rows overlap; local bands
    # don't.
    eigenvector = np.zeros((2, 4, 4), dtype=np.complex128)
    for k in range(2):
        eigenvector[k, 0, 0] = eigenvector[k, 1, 1] = 1.0
        eigenvector[k, 2, 2] = eigenvector[k, 3, 3] = 1.0
    stepped_occupation = np.zeros((2, 4), dtype=np.float64)
    for n_row in range(2):
        stepped_occupation[n_row, 0] = 1.0  # up band 0
        stepped_occupation[n_row, 3] = 1.0  # down band 1
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    column_mu_group = np.array([0, 0, 1, 1], dtype=np.int64)

    with pytest.raises(ValueError, match="pair-closure"):
        build_amplitudes(
            wavevector_index=wv,
            eigenvector=eigenvector,
            stepped_occupation=stepped_occupation,
            column_spin=column_spin,
            column_mu_group=column_mu_group,
            site_positions=site_positions,
            norb_orig=1,
            theta=theta,
            L=L_folded,
            cell_shape=cell_shape,
            subshape=subshape,
        )


def test_a_down_uses_partner_row_eigenvector_v2_1():
    """`A_down` reads BOTH the plane_wave AND the eigenvector at the
    partner row, not the up's row. The two disagree whenever the
    partner-row eigenvector column carries a non-trivial complex
    entry (as in the SubShape=[2,1,1] APBC L=8 fixture at
    tilde_k_folded=-pi/2, where ``v[3, 3, col_down] = -0.5 + 0.5j``).

    Synthetic setup: SubShape=[2,1,1] APBC L=4 with an eigenvector
    that has DIFFERENT complex phases at (row, partner) pairs so that
    "use partner row" vs "use n_row" produce distinguishable A_down.
    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    L_phys = 4
    subshape = np.array([2, 1, 1], dtype=np.int64)
    cell_shape = np.array([L_phys, 1, 1], dtype=np.int64)
    L_folded = cell_shape // subshape  # (2, 1, 1) -> nvol_folded = 2
    wv = np.array([[v, 0, 0] for v in _klist(L_folded[0])], dtype=np.int64)
    site_positions = np.array(
        [[i, 0, 0] for i in range(L_phys)], dtype=np.int64
    )
    theta = np.array([np.pi, 0.0, 0.0], dtype=np.float64)  # APBC

    # nd = 4 (2 up rows, 2 down rows). Set complex partner-distinguishing
    # values on the down block. Use both rows at each k to test summation.
    # up at k=0 (row 0, col 0) with v = (1, 0)^T on up rows.
    # up at k=partner (row 1, col 0) with same trivial v so pair closure
    # holds. Down block: use v[row 0, col 2] = (1, 0), v[row 1, col 2] =
    # (1j, 1). The 1j entry on partner row's down orbital 0 is the
    # discriminator: the correct A_down[..., alpha=(n_row=0)] reads
    # eigenvector[partner=row 1, down_row, col=2], picking up the 1j.
    n_k0 = list(wv[:, 0]).index(0)
    n_partner0 = list(wv[:, 0]).index(-1)  # APBC partner of row 0
    eigenvector = np.zeros((2, 4, 4), dtype=np.complex128)
    # up block on both rows: v[k, 0, 0] = 1 (up band 0 lives on aa=0).
    for k in range(2):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0  # up band 1 unused
    # down block: at n_k0 row, v[..., 2, 2] = 1 (real). At partner row,
    # v[..., 2, 2] = 1j (complex — this is what the partner-row lookup
    # must find).
    eigenvector[n_k0, 2, 2] = 1.0
    eigenvector[n_k0, 3, 3] = 1.0  # unused
    eigenvector[n_partner0, 2, 2] = 1j
    eigenvector[n_partner0, 3, 3] = 1j  # unused

    stepped_occupation = np.zeros((2, 4), dtype=np.float64)
    stepped_occupation[n_k0, 0] = 1.0
    stepped_occupation[n_k0, 2] = 1.0
    stepped_occupation[n_partner0, 0] = 1.0
    stepped_occupation[n_partner0, 2] = 1.0
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    column_mu_group = np.array([0, 0, 1, 1], dtype=np.int64)

    A_up, A_down = build_amplitudes(
        wavevector_index=wv,
        eigenvector=eigenvector,
        stepped_occupation=stepped_occupation,
        column_spin=column_spin,
        column_mu_group=column_mu_group,
        site_positions=site_positions,
        norb_orig=1,
        theta=theta,
        L=L_folded,
        cell_shape=cell_shape,
        subshape=subshape,
    )

    # 2 columns (one per occupied up state). A_down[:, alpha=(n_k0)]
    # picks up the 1j from partner row's down eigenvector.
    assert A_up.shape == (L_phys, 2)
    assert A_down.shape == (L_phys, 2)
    # A_down[:, alpha=(n_k0)] non-zero site: r=0 (folded_orb=0). Value:
    #   plane_wave_down[partner=n_partner0, r=0] * v[partner, 2, 2]
    #   = plane_wave_up[partner, 0] * 1j
    # plane_wave_up[partner=n_partner0, r=0] = exp(+i * kf * 0) * exp(0) * 1/sqrt(2) = 0.707
    # So A_down[0, alpha=(n_k0)] = 0.707 * 1j = 0.707j.
    #
    # If the implementation regressed to plane_wave_down[n_k0, r=0] *
    # v[partner, 2, 2] the plane wave value would still be 0.707 (same
    # at r=0 for both rows in the trivial folded_cell=0 case), so a
    # regression at r=0 would be masked. Test at r=1 (folded_cell=0,
    # sub_offset=1, row_down=3) which is unoccupied on the down side of
    # eigenvector[partner, 3, 2] = 0. Then A_down[1, alpha=(n_k0)] = 0.
    #
    # At r=2 (folded_cell=1, sub_offset=0, row_down=2): plane_wave_up[
    # partner=n_partner0, 2] = exp(+i * kf_partner * 1) * exp(+i*pi*2/4)
    # * 1/sqrt(2). kf_partner = 2*pi*(-1)/2 = -pi. exp(-i*pi) = -1.
    # exp(+i*pi/2) = +1j. So plane_wave = (-1)*(+1j)*(1/sqrt(2)) = -1j/sqrt(2).
    # A_down[2, alpha=(n_k0)] = -1j/sqrt(2) * 1j = 1/sqrt(2).
    #
    # If the impl used plane_wave_down[n_k0, 2] instead:
    # kf_n_k0 = 0, exp(0) * exp(+i pi/2) * (1/sqrt(2)) = +1j/sqrt(2).
    # A_down[2, alpha=(n_k0)] = +1j/sqrt(2) * 1j = -1/sqrt(2). DIFFERENT.
    # The r=2 assertion below distinguishes the two implementations.
    partner_kf_r2_pw = (
        np.exp(+1j * 2 * np.pi * (-1) / 2 * 1)
        * np.exp(+1j * np.pi * 2 / L_phys)
        / np.sqrt(2)
    )
    expected_A_down_r2 = partner_kf_r2_pw * 1j
    assert A_down[2, 0] == pytest.approx(
        complex(expected_A_down_r2), abs=1e-12
    )
