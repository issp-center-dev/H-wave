"""Tests for (k, -k) partner index lookup.

See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import numpy as np
import pytest

# Make tools/ importable. Tests are run from the repo root.
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.partner_index import (
    find_partner_rows,
)


def _klist(n):
    """Mirror UHFk._init_wavevec._klist."""
    return np.roll(np.arange(n) - (n // 2), -(n // 2))


def test_partner_pbc_1d_l4():
    """PBC 1D L=4 ⇒ wavevector_index = [0, 1, -2, -1]; partners:
    n=0 self-paired; n=1 ↔ n=3 (-1); n=2 self-paired (L/2 = 2 → -2)."""
    wv = np.array([[v, 0, 0] for v in _klist(4)], dtype=np.int64)
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)  # PBC
    L = np.array([4, 1, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)

    # n=0 (wv=0): partner is itself
    assert partners[0] == 0 and is_self[0]
    # n=2 (wv=-2 ≡ 2 mod 4): partner satisfies -(-2) ≡ 2 mod 4 → self
    assert partners[2] == 2 and is_self[2]
    # n=1 (wv=1): partner satisfies -(1) ≡ -1 ≡ 3 mod 4; wv=-1 is at row 3
    assert partners[1] == 3 and not is_self[1]
    # n=3 (wv=-1): partner satisfies -(-1) ≡ 1; wv=1 is at row 1
    assert partners[3] == 1 and not is_self[3]


def test_partner_apbc_1d_l4_no_self_pair():
    """APBC 1D L=4 ⇒ no self-pair; partner ≡ -n - 1 mod 4."""
    wv = np.array([[v, 0, 0] for v in _klist(4)], dtype=np.int64)
    theta = np.array([np.pi, 0.0, 0.0], dtype=np.float64)
    L = np.array([4, 1, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)

    assert not np.any(is_self), "APBC even L has no self-pair"
    # Verify partner-of-partner is identity
    for n_row in range(4):
        assert partners[partners[n_row]] == n_row


def test_partner_apbc_1d_l5_has_self_pair():
    """APBC odd L=5 ⇒ self-pair at residue class n = (L-1)/2 = 2."""
    wv = np.array([[v, 0, 0] for v in _klist(5)], dtype=np.int64)
    theta = np.array([np.pi, 0.0, 0.0], dtype=np.float64)
    L = np.array([5, 1, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)

    # _klist(5) = [0, 1, 2, -2, -1]; residues mod 5: [0, 1, 2, 3, 4]
    # For APBC partner residue is (-n - 1) mod 5
    # n=2 (residue 2): partner residue = (-2 - 1) mod 5 = 2 → self
    self_idx = list(wv[:, 0]).index(2)
    assert is_self[self_idx]


def test_partner_pbc_2d_4x4():
    """PBC 2D L=4x4: self-pair locus = {(0,0), (0,2), (2,0), (2,2)} residues."""
    wv = np.array(
        [[ix, iy, 0] for ix in _klist(4) for iy in _klist(4)], dtype=np.int64
    )
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L = np.array([4, 4, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)

    self_residues = {
        tuple(int(v[d] % L[d]) for d in range(3))
        for v in wv[is_self]
    }
    expected = {(0, 0, 0), (0, 2, 0), (2, 0, 0), (2, 2, 0)}
    assert self_residues == expected


def test_partner_mixed_pbc_apbc_2d():
    """Mixed PBC-x + APBC-y, L=4: x has self-pair at {0, 2}; y has none.
    Total k self-pair requires both directions self-paired → none."""
    wv = np.array(
        [[ix, iy, 0] for ix in _klist(4) for iy in _klist(4)], dtype=np.int64
    )
    theta = np.array([0.0, np.pi, 0.0], dtype=np.float64)
    L = np.array([4, 4, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)
    assert not np.any(is_self)
    # partner-of-partner identity
    for n_row in range(len(wv)):
        assert partners[partners[n_row]] == n_row


def test_partner_apbc_folded_l4_from_l8_subshape_2():
    """L=8 CellShape reduced by SubShape=[2,1,1] to folded L_folded=4.
    APBC partner should follow (-n - 1) mod 4."""
    L_folded = 4
    wv = np.array([[v, 0, 0] for v in _klist(L_folded)], dtype=np.int64)
    theta = np.array([np.pi, 0.0, 0.0], dtype=np.float64)
    L = np.array([L_folded, 1, 1], dtype=np.int64)

    partners, is_self = find_partner_rows(wv, theta, L)
    assert not np.any(is_self)
    for n_row in range(L_folded):
        assert partners[partners[n_row]] == n_row
