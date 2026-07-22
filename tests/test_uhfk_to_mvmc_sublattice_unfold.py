"""Sublattice unfold helpers with 3D generic coverage.

Adversarial cases:
- Down-row block offset (folded_row_indices)
- Site decoding driven by coordinates, not by row-major ordinal
- SubShape[d] must divide CellShape[d]
"""
from __future__ import annotations

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

from tools._uhfk_to_mvmc.sublattice_unfold import (
    decode_physical_site,
    encode_folded_orbital,
    folded_row_indices,
)


def test_decode_encode_round_trip_1d():
    """1D SubShape=[2,1,1]: sub_offset should recover from any coord."""
    subshape = np.array([2, 1, 1], dtype=np.int64)
    # site at cx=3 (unfolded), Sx=2 -> folded_cell=1, sub_offset=(1,0,0)
    coord = np.array([3, 0, 0], dtype=np.int64)
    folded_cell, sub_offset = decode_physical_site(coord, subshape)
    np.testing.assert_array_equal(folded_cell, [1, 0, 0])
    np.testing.assert_array_equal(sub_offset, [1, 0, 0])

    # encode round-trip: (0, (1,0,0)) -> 0 + 1*(1 + 2*0 + ...) = 1
    idx = encode_folded_orbital(0, sub_offset, norb_orig=1, subshape=subshape)
    assert idx == 1


def test_decode_encode_round_trip_3d():
    """3D SubShape=[2,2,2]: (sx,sy,sz)=(1,1,1) -> encoded index = 1+2+4=7."""
    subshape = np.array([2, 2, 2], dtype=np.int64)
    coord = np.array([3, 5, 7], dtype=np.int64)
    folded_cell, sub_offset = decode_physical_site(coord, subshape)
    np.testing.assert_array_equal(folded_cell, [1, 2, 3])
    np.testing.assert_array_equal(sub_offset, [1, 1, 1])

    idx = encode_folded_orbital(0, sub_offset, norb_orig=1, subshape=subshape)
    # 0 + 1 * (1 + 2*(1 + 2*1)) = 0 + 1*(1 + 6) = 7
    assert idx == 7


def test_folded_row_indices_down_has_offset():
    """Down row must include norb_folded offset, not just folded_orb.

    With norb_folded=2, including or omitting the block offset gives
    distinguishable rows. See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    # norb_folded=2 (single orbital × subvol=2)
    row_up, row_down = folded_row_indices(folded_orb=0, norb_folded=2)
    assert row_up == 0
    assert row_down == 2, "down row must have +norb_folded offset"

    row_up, row_down = folded_row_indices(folded_orb=1, norb_folded=2)
    assert row_up == 1
    assert row_down == 3


def test_decode_non_row_major_coordinates():
    """Site 5 in a lexicographic ordering has coord (2, 1, 0) for
    Lx=3, Ly=2; but geometry_uhf.dat may list sites in a different
    order. decode_physical_site must use the coordinate directly, not
    infer it from the row-major ordinal.
    """
    subshape = np.array([1, 2, 1], dtype=np.int64)
    coord_a = np.array([2, 1, 0], dtype=np.int64)
    fa, sa = decode_physical_site(coord_a, subshape)
    np.testing.assert_array_equal(fa, [2, 0, 0])
    np.testing.assert_array_equal(sa, [0, 1, 0])

    # Different site with a permuted layout hits the same sub_offset
    # only if the coord itself matches. The helper does NOT depend on
    # any global site index.
    coord_b = np.array([0, 1, 0], dtype=np.int64)
    fb, sb = decode_physical_site(coord_b, subshape)
    np.testing.assert_array_equal(fb, [0, 0, 0])
    np.testing.assert_array_equal(sb, [0, 1, 0])


# ---------------------------------------------------------------------
# unfold_amplitude_columns
# ---------------------------------------------------------------------

from tools._uhfk_to_mvmc.sublattice_unfold import unfold_amplitude_columns


def _klist(n):
    return np.roll(np.arange(n) - (n // 2), -(n // 2))


def test_unfold_amplitude_columns_1d_pbc_l4():
    """L=4 SubShape=[1,1,1] PBC gives the positive-Bloch plane wave.

    ``plane_wave_up[k, i] = exp(+i k r_i) / sqrt(L)``. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    L = 4
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array([[i, 0, 0] for i in range(L)], dtype=np.int64)
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    cell_shape = np.array([L, 1, 1], dtype=np.int64)
    subshape = np.array([1, 1, 1], dtype=np.int64)

    pw_up, pw_dn, row_up, row_dn = unfold_amplitude_columns(
        wv, cell_shape, subshape, site_positions, norb_orig=1, theta=theta,
    )

    assert pw_up.shape == (L, L)
    np.testing.assert_array_equal(row_up, np.zeros(L, dtype=np.int64))
    # norb_folded = 1 -> down row index = 1
    np.testing.assert_array_equal(row_dn, np.ones(L, dtype=np.int64))

    # k=0, i=1: exp(+i*0*1)/sqrt(4) = 0.5
    assert pw_up[list(wv[:, 0]).index(0), 1] == pytest.approx(0.5 + 0.0j)
    # k=1 (k_folded = 2pi/4), i=1: exp(+i*pi/2)/2 = +0.5i.
    n_k1 = list(wv[:, 0]).index(1)
    assert pw_up[n_k1, 1] == pytest.approx(0.0 + 0.5j, abs=1e-12)


def test_unfold_amplitude_columns_1d_subshape_2():
    """L=4 SubShape=[2,1,1]: folded lattice has nvol_folded=2. Each
    physical site maps to (folded_cell, sub_offset). Row indices track
    up/down block offsets."""
    L = 4
    subshape = np.array([2, 1, 1], dtype=np.int64)
    cell_shape = np.array([L, 1, 1], dtype=np.int64)
    L_folded = cell_shape // subshape          # (2, 1, 1)
    wv = np.array([[v, 0, 0] for v in _klist(L_folded[0])], dtype=np.int64)
    site_positions = np.array([[i, 0, 0] for i in range(L)], dtype=np.int64)
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)

    pw_up, pw_dn, row_up, row_dn = unfold_amplitude_columns(
        wv, cell_shape, subshape, site_positions, norb_orig=1, theta=theta,
    )

    # norb_folded = subvol = 2, so up rows are {0, 1}, down rows are {2, 3}.
    assert pw_up.shape == (2, L)
    # site 0: coord (0,0,0) -> folded_cell (0,0,0), sub (0,0,0) -> folded_orb 0
    assert row_up[0] == 0 and row_dn[0] == 2
    # site 1: coord (1,0,0) -> folded_cell (0,0,0), sub (1,0,0) -> folded_orb 1
    assert row_up[1] == 1 and row_dn[1] == 3
    # site 2: coord (2,0,0) -> folded_cell (1,0,0), sub (0,0,0) -> folded_orb 0
    assert row_up[2] == 0 and row_dn[2] == 2
    # site 3: coord (3,0,0) -> folded_cell (1,0,0), sub (1,0,0) -> folded_orb 1
    assert row_up[3] == 1 and row_dn[3] == 3

    # plane_wave_up[k=0, i=2] = exp(+i*0*1)/sqrt(2) = 1/sqrt(2)
    n_k0 = list(wv[:, 0]).index(0)
    assert pw_up[n_k0, 2] == pytest.approx(1.0 / np.sqrt(2), abs=1e-12)


def test_unfold_plane_wave_down_equals_up_v2_1():
    """``plane_wave_down`` shares the SAME positive-Bloch
    envelope as ``plane_wave_up``. The (k, -k) time-reversal partnership
    is resolved via the caller's ``eigenvector[partner_n, ...]`` lookup,
    not via a conjugated plane wave. This prevents a regression that
    reintroduces
    reintroducing ``conj(plane_wave_up_folded) * exp(-i theta r / L)``
    and silently breaks SubShape > 1 APBC density. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    # 1D APBC SubShape=[2,1,1] — the exact fixture whose density check
    # historically failed under the negative-Bloch down envelope.
    L = 8
    subshape = np.array([2, 1, 1], dtype=np.int64)
    cell_shape = np.array([L, 1, 1], dtype=np.int64)
    L_folded = cell_shape // subshape  # (4, 1, 1)
    wv = np.array([[v, 0, 0] for v in _klist(L_folded[0])], dtype=np.int64)
    site_positions = np.array([[i, 0, 0] for i in range(L)], dtype=np.int64)
    theta = np.array([np.pi, 0.0, 0.0], dtype=np.float64)  # APBC

    pw_up, pw_dn, _row_up, _row_dn = unfold_amplitude_columns(
        wv, cell_shape, subshape, site_positions, norb_orig=1, theta=theta,
    )

    # The down envelope is the up envelope element-wise.
    np.testing.assert_array_almost_equal(pw_dn, pw_up, decimal=14)
