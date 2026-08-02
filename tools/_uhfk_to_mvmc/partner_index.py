"""(k, -k) partner index lookup against H-wave UHFk wavevector_index.

H-wave emits wavevector_index in signed integer form
(_klist convention: np.roll(np.arange(n) - n//2, -(n//2))). The partner
of row n is the row whose emitted index m satisfies
    m ≡ -n - 2 * twist_offset_d (mod L_d)
per direction d. Self-pair: m_d ≡ n_d for all d. See
docs/en/source/algorithm/uhfk_to_mvmc.rst.
"""
from __future__ import annotations

import numpy as np


def find_partner_rows(
    wavevector_index: np.ndarray,
    theta: np.ndarray,
    L: np.ndarray,
):
    """Return ``(partner_rows, is_self_pair)``.

    Parameters
    ----------
    wavevector_index : (nvol, 3) int
        Signed wavevector indices as emitted by ``UHFk._init_wavevec``.
    theta : (3,) float
        Twist angle per direction (0 for PBC, pi for APBC).
    L : (3,) int
        ``CellShape``.

    Returns
    -------
    partner_rows : (nvol,) int
        For each row n, the row whose wavevector_index matches the
        ``(k, -k)`` partner residue class.
    is_self_pair : (nvol,) bool
        True if ``partner_rows[n] == n``.
    """
    wavevector_index = np.asarray(wavevector_index, dtype=np.int64)
    theta = np.asarray(theta, dtype=np.float64)
    L = np.asarray(L, dtype=np.int64)

    # twist_offset[d] in {0, 1/2}: 0 for PBC, 1/2 for APBC
    twist_offset = theta / (2.0 * np.pi)
    # 2 * twist_offset is 0 or 1
    shift = np.rint(2.0 * twist_offset).astype(np.int64)

    nvol = wavevector_index.shape[0]
    partner_rows = np.full(nvol, -1, dtype=np.int64)

    # Build a lookup: (residue tuple) → row index
    res_lookup = {}
    for n_row, wv in enumerate(wavevector_index):
        key = tuple(int(v % L[d]) for d, v in enumerate(wv))
        res_lookup[key] = n_row

    for n_row, wv in enumerate(wavevector_index):
        m_residue = tuple(
            int((-int(wv[d]) - int(shift[d])) % int(L[d])) for d in range(3)
        )
        if m_residue not in res_lookup:
            raise ValueError(
                f"partner residue {m_residue} not present in wavevector_index "
                f"(row {n_row}, wv={wv.tolist()}, theta={theta.tolist()}, "
                f"L={L.tolist()}); check that the (theta, L) combination "
                f"matches the emitted grid"
            )
        partner_rows[n_row] = res_lookup[m_residue]

    is_self_pair = partner_rows == np.arange(nvol)
    return partner_rows, is_self_pair
