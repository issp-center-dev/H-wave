"""Allowlist predicate shared by the CLI and static coverage checker.

Any change to the allowed shape classes MUST land here; both callers
(tools/uhfk_to_mvmc.py CLI + tools/_uhfk_to_mvmc/allowlist_coverage.py)
share this predicate so they cannot drift.
"""
from __future__ import annotations

import numpy as np


_V37_ALLOWED_APBC_MASKS = frozenset({
    (1, 1, 0),  # xy: case_soc_rashba_3d_sub_apbc_xy
    (1, 0, 1),  # xz: case_soc_rashba_3d_sub_apbc_xz
    (0, 1, 1),  # yz: case_soc_rashba_3d_sub_apbc_yz
    (1, 1, 1),  # xyz: case_soc_rashba_3d_sub_apbc_xyz
})
_V36_ALLOWED_APBC_MASKS = frozenset({
    (1, 0, 0), (0, 1, 0), (0, 0, 1),  # single-direction APBC
})

_V36_LATTICE = ((2, 2, 1), (6, 4, 1))  # (sub_shape, cell_shape)
_V37_LATTICE = ((2, 2, 2), (4, 4, 4))


def apbc_mask_of(theta):
    """Return the 0/1 mask indicating which theta components are ±pi."""
    theta = np.asarray(theta, dtype=np.float64)
    return tuple(int(np.abs(np.abs(t) - np.pi) < 1e-12) for t in theta)


def is_supported_triple(theta, sub_shape, cell_shape, is_soc_mode):
    """Return True iff the (theta, sub_shape, cell_shape, is_soc_mode)
    triple is supported. False -> CLI reject.

    Early-return branches for non-SOC / SubShape=[1,1,1] / SOC+PBC are
    preserved unchanged; only SOC + APBC + SubShape > 1 hits the
    allowlist frozensets.
    """
    sub_shape = tuple(int(s) for s in sub_shape)
    cell_shape = tuple(int(c) for c in cell_shape)
    if not is_soc_mode:
        return True
    if not any(s != 1 for s in sub_shape):
        return True
    apbc_mask = apbc_mask_of(theta)
    if apbc_mask == (0, 0, 0):
        return True
    if (apbc_mask in _V36_ALLOWED_APBC_MASKS
            and sub_shape == _V36_LATTICE[0]
            and cell_shape == _V36_LATTICE[1]):
        return True
    if (apbc_mask in _V37_ALLOWED_APBC_MASKS
            and sub_shape == _V37_LATTICE[0]
            and cell_shape == _V37_LATTICE[1]):
        return True
    return False


REJECT_MESSAGE = (
    "ERROR: unsupported SOC + APBC + SubShape combination. Supported "
    "combinations are single-direction (x, y, or z) APBC with "
    "CellShape=[6,4,1]/SubShape=[2,2,1], or multi-direction "
    "(xy, xz, yz, or xyz) APBC with "
    "CellShape=[4,4,4]/SubShape=[2,2,2]. Add a fixture and gate "
    "validation before expanding the allowlist."
)
