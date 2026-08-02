"""T=0 occupation step rounding with finite-T guards.

The bridge consumes occupation.npz, sorts each mu-group's
eigenvalues, and projects to a T=0 step function. fail-fast on:
  - any state with f in (FRACTIONAL_LOW, FRACTIONAL_HIGH) — there are
    more than one such state per mu-group → SCF was not close to a Slater
    determinant
  - Fermi-level degeneracy (ncond-th vs (ncond+1)-th eigenvalue equal)
  - Ne per group disagreeing with input toml's Ncond/2Sz derivation
"""
from __future__ import annotations

import numpy as np

FRACTIONAL_LOW = 1e-3
FRACTIONAL_HIGH = 1.0 - FRACTIONAL_LOW
DEGEN_TOL = 1e-10


class OccupationGuardError(RuntimeError):
    """Raised when finite-T residual or degeneracy makes T=0 projection
    unsafe."""


def step_occupation(
    occupation: np.ndarray,
    eigenvalue: np.ndarray,
    column_spin: np.ndarray,
    column_mu_group: np.ndarray,
    T: float,
    ncond_per_group,
    is_soc_mode: bool = False,
):
    """Return ``(stepped_occupation, summary)``.

    Parameters
    ----------
    occupation : (nvol, nd) float
        From occupation.npz.
    eigenvalue : (nvol, nd) float
        From eigen.npz.
    column_spin : (nd,) int
        From occupation.npz. -1 (Sz-free / mixed) is rejected here for
        the non-SOC path; accepted when ``is_soc_mode=True``.
    column_mu_group : (nd,) int
        Mu-group index per column.
    T : float
        SCF temperature; used only for the residual check threshold.
    ncond_per_group : list[int]
        Expected occupied count per mu-group (from input toml's Ncond/2Sz).
    is_soc_mode : bool, default False
        SOC dispatch flag threaded through by the CLI. When True,
        ``column_spin = -1`` (single mixed-spin block with a single
        global mu-group) is accepted and the lowest ``Ncond``
        eigenvalues are occupied regardless of spin label.

    Returns
    -------
    stepped_occupation : (nvol, nd) float
        0 or 1, with the lowest ``ncond_per_group[g]`` states (sorted
        across all (k, n) in group g) marked occupied.
    summary : dict with keys ``ne_per_group``, ``fractional_count``,
        ``max_residual``.
    """
    occupation = np.asarray(occupation, dtype=np.float64)
    eigenvalue = np.asarray(eigenvalue, dtype=np.float64)
    column_spin = np.asarray(column_spin, dtype=np.int64)
    column_mu_group = np.asarray(column_mu_group, dtype=np.int64)

    if np.any(column_spin < 0) and not is_soc_mode:
        raise OccupationGuardError(
            "occupation_step encountered column_spin = -1 (Sz-free / "
            "mixed block); this requires is_soc_mode=True. Enable the "
            "SOC path or reject the input before calling occupation_step"
        )

    n_groups = len(ncond_per_group)
    stepped = np.zeros_like(occupation)
    ne_per_group = []
    fractional_count = 0
    max_residual = 0.0

    for g in range(n_groups):
        cols = np.where(column_mu_group == g)[0]
        if len(cols) == 0:
            raise OccupationGuardError(
                f"mu-group {g} has no eigenvector columns"
            )
        ws_block = eigenvalue[:, cols]            # (nvol, blk_size)
        fs_block = occupation[:, cols]            # (nvol, blk_size)
        flat_w = ws_block.reshape(-1)
        flat_f = fs_block.reshape(-1)
        order = np.argsort(flat_w, kind="stable")

        ncond_g = int(ncond_per_group[g])
        if ncond_g < 0 or ncond_g > flat_w.size:
            raise OccupationGuardError(
                f"mu-group {g} ncond={ncond_g} out of range [0, "
                f"{flat_w.size}] derived from input toml"
            )

        # Fermi-level degeneracy: ncond-th vs (ncond+1)-th
        if 0 < ncond_g < flat_w.size:
            w_homo = flat_w[order[ncond_g - 1]]
            w_lumo = flat_w[order[ncond_g]]
            if abs(w_homo - w_lumo) < DEGEN_TOL:
                raise OccupationGuardError(
                    f"mu-group {g} Fermi-level degeneracy: "
                    f"w[{ncond_g-1}]={w_homo:.6e}, "
                    f"w[{ncond_g}]={w_lumo:.6e} differ by < {DEGEN_TOL}"
                )

        # Fractional residual: count states where SCF occupation is not
        # clearly 0 or 1.
        frac_mask = (flat_f > FRACTIONAL_LOW) & (flat_f < FRACTIONAL_HIGH)
        fractional_count += int(frac_mask.sum())
        max_residual = max(
            max_residual,
            float(np.max(np.minimum(flat_f, 1.0 - flat_f)))
        )
        if frac_mask.sum() >= 2:
            raise OccupationGuardError(
                f"mu-group {g} has {int(frac_mask.sum())} fractional "
                f"occupations (T={T}); at most one is supported. Rerun "
                "SCF with smaller T or use the BCS-like path"
            )

        # Step rounding: lowest ncond_g positions occupied
        mask_flat = np.zeros_like(flat_w, dtype=np.float64)
        mask_flat[order[:ncond_g]] = 1.0
        stepped_block = mask_flat.reshape(ws_block.shape)

        # Place back into the full (nvol, nd) array
        for ci, col_n in enumerate(cols):
            stepped[:, col_n] = stepped_block[:, ci]
        ne_per_group.append(ncond_g)

    summary = {
        "ne_per_group": ne_per_group,
        "fractional_count": fractional_count,
        "max_residual": max_residual,
    }
    return stepped, summary
