"""Bridge-side BoundaryCondition parsing + eigen twist consistency.

Delegates alias/case handling to hwave.solver._apbc_phase.normalize_boundary_condition
(the single H-wave source of truth). Adds only:

  - absent-key default (matches uhfk.py:243-245)
  - container-shape / element-type checks that H-wave's helper does not perform
  - eigen.npz twist_offset consistency guard
"""
from __future__ import annotations

import numpy as np

from hwave.solver._apbc_phase import (
    normalize_boundary_condition,
    twist_offset,
)


class BoundaryInputError(ValueError):
    """Raised by bridge-side pre-dispatch boundary validation."""


_TWIST_TOL = 1e-12


def normalize_boundary_condition_list(entries):
    """Return the (theta_x, theta_y, theta_z) tuple for ``entries``.

    - ``None`` (absent key) -> all-PBC default (matches UHFk).
    - Container must be a list/tuple; every element must be a string.
    - Delegates alias handling + length check to H-wave's helper.
    """
    if entries is None:
        entries = ["periodic", "periodic", "periodic"]
    if not isinstance(entries, (list, tuple)):
        raise BoundaryInputError(
            f"BoundaryCondition must be a list or tuple of strings; "
            f"got {type(entries).__name__}"
        )
    for i, e in enumerate(entries):
        if not isinstance(e, str):
            raise BoundaryInputError(
                f"BoundaryCondition entry at index {i} must be a str; "
                f"got {type(e).__name__} ({e!r})"
            )
    return normalize_boundary_condition(list(entries))


def check_eigen_twist_consistency(theta_from_toml, twist_offset_from_eigen):
    """Reject if eigen.npz's ``twist_offset`` does not match ``theta_from_toml``.

    Parameters
    ----------
    theta_from_toml : tuple[float, float, float]
        Output of :func:`normalize_boundary_condition_list`.
    twist_offset_from_eigen : numpy.ndarray | None
        Value stored in ``eigen.npz["twist_offset"]``, or ``None`` if the
        field is absent (legacy H-wave runs).
    """
    theta_arr = np.asarray(theta_from_toml, dtype=np.float64)
    expected_twist = twist_offset(theta_arr)
    if twist_offset_from_eigen is None:
        # Legacy eigen.npz. Accept only when input toml claims all-PBC
        # (theta == 0); otherwise refuse because we cannot verify
        # consistency.
        if not np.allclose(theta_arr, 0.0, atol=_TWIST_TOL):
            raise BoundaryInputError(
                "eigen.npz does not contain twist_offset (legacy H-wave "
                f"output); input.toml BoundaryCondition implies "
                f"theta = {tuple(theta_arr)}, which requires eigen.npz "
                "produced with twist_offset present"
            )
        return
    twist_saved = np.asarray(twist_offset_from_eigen, dtype=np.float64)
    if twist_saved.shape != expected_twist.shape:
        raise BoundaryInputError(
            f"eigen.npz twist_offset shape {twist_saved.shape} does not "
            f"match expected shape {expected_twist.shape}"
        )
    if not np.allclose(twist_saved, expected_twist, atol=_TWIST_TOL):
        raise BoundaryInputError(
            f"eigen.npz was produced under twist_offset = "
            f"{tuple(twist_saved)}, but input.toml BoundaryCondition "
            f"canonicalizes to twist_offset = {tuple(expected_twist)}. "
            "Regenerate eigen.npz with matching BoundaryCondition, or "
            "update input.toml to match the saved eigen."
        )
