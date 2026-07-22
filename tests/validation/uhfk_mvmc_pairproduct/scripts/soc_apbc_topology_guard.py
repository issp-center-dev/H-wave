"""SOC + APBC + SubShape topology and mutation guard for the G4 gate.

Semantic contract:

  1. Load ``composite_element.json`` via ``--composite-manifest``.
  2. Load CURRENT-run ``${WORK_DIR}/hwave/green.npz`` (NEVER trust
     manifest values alone).
  3. Assert the manifest's ``(i_c, s_c, j_c, t_c)`` element exists on
     the current run with ``abs(G_current) >= 0.8 * G_c_abs``.
  4. Assert ``s_c != t_c``, ``sub_offset_d(i_c) != sub_offset_d(j_c)``
     for every active APBC direction ``d`` (``theta_d != 0``), and
     ``abs(G_current) >= 1e-3``.
  5. Run all mutations on the CURRENT run at the composite element;
     assert each ``delta_M >= T_M`` per the manifest's
     ``T_M_per_mutation``. The manifest schema is auto-detected (see
     ``_iter_mutation_ids``): the 10-entry whole-vector mutators
     (``M-gauge-1..5`` + ``M-ship-1..5``) or the 30-entry
     per-direction mutators (``M-gauge-1-x`` .. ``M-ship-5-z``).
  6. On success emit the anchored PASS record. On any
     failure exit code 2 + empty stdout.

See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np


# Single source of truth for G4 metadata.
G4_MODE = "g4"
G4_ARTIFACT_SOURCE = "hwave+bridge+composite-manifest"
G4_HELPER = "soc_apbc_topology_guard"


# Whole-vector mutation ids (10 entries). Each mutator
# below is applied to the FULL theta/L/dr vector at once.
_V36_MUTATION_IDS = tuple(
    f"M-gauge-{n}" for n in range(1, 6)
) + tuple(
    f"M-ship-{n}" for n in range(1, 6)
)

# Per-direction mutation ids (30 entries). The
# 10 mutation kinds are split into x/y/z variants that mutate only the
# named axis component, leaving the other two axes baseline. M-4 has
# intentionally different schema-specific semantics: the
# whole-vector id flips sub_offset. Per-direction M-gauge-4 omits sub_offset
# on its named axis, while M-ship-4 halves that contribution so it is
# distinct from M-ship-5, which omits it. See the mutator branches below.
_V37_MUTATION_IDS = tuple(
    f"M-gauge-{n}-{d}"
    for n in range(1, 6) for d in ("x", "y", "z")
) + tuple(
    f"M-ship-{n}-{d}"
    for n in range(1, 6) for d in ("x", "y", "z")
)


def _iter_mutation_ids(manifest):
    """Return the mutation-id tuple appropriate for the manifest's
    ``T_M_per_mutation`` schema. Requires an exact-key match against
    30-entry or 10-entry schema; otherwise it raises because the
    manifest matches neither known schema."""
    t_m = manifest.get("T_M_per_mutation")
    if not isinstance(t_m, dict):
        raise ValueError(
            f"manifest 'T_M_per_mutation' must be a dict; got "
            f"{type(t_m).__name__}"
        )
    keys = set(t_m.keys())
    v37 = set(_V37_MUTATION_IDS)
    v36 = set(_V36_MUTATION_IDS)
    if keys == v37:
        return _V37_MUTATION_IDS
    if keys == v36:
        return _V36_MUTATION_IDS
    # Fall through: hybrid or partial schema.
    if keys < v37 and keys > v36:
        raise ValueError(
            f"manifest T_M_per_mutation looks partially upgraded "
            f"from 10 to 30 entries ({len(keys)} keys, expected 10 or 30). "
            f"Missing per-direction keys: {sorted(v37 - keys)}. "
            f"Extra: {sorted(keys - v37)}."
        )
    raise ValueError(
        f"manifest T_M_per_mutation matches neither the 10-key nor "
        f"30-key schema; got {len(keys)} keys: {sorted(keys)}"
    )


def _split_mutator_id(mutator):
    """Split a mutator id into ``(base_id, axis)``.

    ``"baseline"`` and whole-vector ids (``M-gauge-1`` ..
    ``M-ship-5``) return ``axis=None``: the whole-vector mutation (or
    no mutation at all) applies.

    Per-direction ids (e.g. ``M-gauge-1-x``) return
    ``base_id="M-gauge-1"`` and ``axis=0``: the mutation applies to
    that single vector component only. Mutations 1/2/3/5
    splice the whole-vector formula into an otherwise-baseline vector; M-4
    instead changes the named-axis sub_offset coefficient because the
    sign flip is structurally invisible on the shipped
    ``L_folded=2`` lattice. M-gauge-4 omits the contribution;
    M-ship-4 halves it, keeping it distinct from M-ship-5's omission.
    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    if mutator == "baseline":
        return mutator, None
    parts = mutator.rsplit("-", 1)
    if len(parts) == 2 and parts[1] in ("x", "y", "z"):
        return parts[0], "xyz".index(parts[1])
    return mutator, None


def _require_finite_float(value, field_name):
    """Return ``value`` as float, rejecting non-numeric and non-finite data."""
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"manifest {field_name} must be numeric; got {value!r}"
        ) from exc
    if not np.isfinite(number):
        raise ValueError(
            f"manifest {field_name} must be finite; got {number!r}"
        )
    return number


def _require_finite_int(value, field_name):
    """Return an exactly integral, finite manifest scalar as ``int``."""
    number = _require_finite_float(value, field_name)
    if not number.is_integer():
        raise ValueError(
            f"manifest {field_name} must be an integer; got {number!r}"
        )
    return int(number)


def _require_finite_vector(manifest, field_name, *, integral):
    """Validate and return one finite length-three manifest vector."""
    try:
        vector = np.asarray(manifest[field_name], dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"manifest {field_name} must be a numeric length-3 vector"
        ) from exc
    if vector.shape != (3,):
        raise ValueError(
            f"manifest {field_name} must have shape (3,); got "
            f"{vector.shape}"
        )
    if not np.all(np.isfinite(vector)):
        raise ValueError(
            f"manifest {field_name} must contain only finite values; "
            f"got {vector.tolist()}"
        )
    if integral and not np.all(vector == np.floor(vector)):
        raise ValueError(
            f"manifest {field_name} must contain only integers; "
            f"got {vector.tolist()}"
        )
    return vector.astype(np.int64) if integral else vector


def _validate_manifest_numeric_policy(manifest):
    """Validate every numeric manifest field consumed by G4.

    Thresholds are derived data, not manifest-controlled policy. Recompute
    ``max(1e-5, 0.10 * G_c_abs)`` from the manifest's finite ``G_c_abs``.
    Every whole-vector entry and every active-axis per-direction entry must
    equal that value exactly; inactive-axis entries must be exactly zero.
    """
    mutation_ids = _iter_mutation_ids(manifest)
    G_c_abs = _require_finite_float(manifest["G_c_abs"], "G_c_abs")
    theta = _require_finite_vector(
        manifest, "theta_radians", integral=False,
    )
    cell_shape = _require_finite_vector(
        manifest, "cell_shape", integral=True,
    )
    subshape = _require_finite_vector(
        manifest, "sub_shape", integral=True,
    )
    if np.any(cell_shape <= 0) or np.any(subshape <= 0):
        raise ValueError(
            "manifest cell_shape and sub_shape entries must be positive"
        )
    for field_name in ("i_c", "s_c", "j_c", "t_c", "N_pairs"):
        _require_finite_int(manifest[field_name], field_name)
    if "ncond" in manifest:
        _require_finite_int(manifest["ncond"], "ncond")
    _require_finite_vector(manifest, "sub_offset_i", integral=True)
    _require_finite_vector(manifest, "sub_offset_j", integral=True)

    expected_active = max(1e-5, 0.10 * G_c_abs)
    thresholds = manifest["T_M_per_mutation"]
    for mutator_id in mutation_ids:
        threshold = _require_finite_float(
            thresholds[mutator_id],
            f"threshold T_M_per_mutation[{mutator_id!r}]",
        )
        _, axis = _split_mutator_id(mutator_id)
        is_active = axis is None or abs(float(theta[axis])) > 1e-12
        expected = expected_active if is_active else 0.0
        if threshold != expected:
            axis_detail = "whole-vector"
            if axis is not None:
                axis_detail = (
                    f"{'active' if is_active else 'inactive'} axis "
                    f"{'xyz'[axis]}"
                )
            raise ValueError(
                f"manifest threshold policy mismatch for {mutator_id} "
                f"({axis_detail}): got {threshold!r}, expected "
                f"{expected!r} = "
                f"{'max(1e-5, 0.10 * G_c_abs)' if is_active else '0.0'}"
            )
    return mutation_ids


def _load_manifest(path):
    with open(path) as fp:
        return json.load(fp)


def _folded_cell(r_phys, subshape):
    return r_phys // subshape


def _gauge_lift_element(green_sublattice, i, s, j, t, subshape, cell_shape,
                        site_positions, boundary_theta, mutator="baseline"):
    """Element-level gauge_lift with optional M-gauge mutation.

    See scripts/phase2_produce_manifest.py for the derivation of the
    per-mutation phase surgery; this guard reproduces those formulas so
    the workspace-run deltas can be compared directly with the manifest
    T_M values.
    """
    gs_soc = green_sublattice[:, 0, :, 0, :]
    L_folded = cell_shape // subshape
    LFX, LFY, LFZ = int(L_folded[0]), int(L_folded[1]), int(L_folded[2])
    gs = gs_soc.reshape(LFX, LFY, LFZ, gs_soc.shape[1], gs_soc.shape[2])
    G_k = np.fft.ifftn(gs, axes=(0, 1, 2), norm="forward")

    r_phys_i = site_positions[i].astype(np.int64)
    r_phys_j = site_positions[j].astype(np.int64)
    fc_i = _folded_cell(r_phys_i, subshape)
    fc_j = _folded_cell(r_phys_j, subshape)
    so_i = r_phys_i - fc_i * subshape
    so_j = r_phys_j - fc_j * subshape
    folded_orb_i = int(
        so_i[0] + subshape[0] * (so_i[1] + subshape[1] * so_i[2])
    )
    folded_orb_j = int(
        so_j[0] + subshape[0] * (so_j[1] + subshape[1] * so_j[2])
    )
    aa = 2 * folded_orb_i + int(s)
    bb = 2 * folded_orb_j + int(t)

    dr_folded_base = (fc_j.astype(np.float64) + so_j.astype(np.float64)) \
                     - (fc_i.astype(np.float64) + so_i.astype(np.float64))
    dr_full = r_phys_j.astype(np.float64) - r_phys_i.astype(np.float64)
    L_full = (subshape * L_folded).astype(np.float64)

    dr_folded = dr_folded_base.copy()
    # NOTE: must be a defensive copy (not a view) -- np.asarray would
    # alias the caller's persistent theta array when it is already
    # float64, and the per-axis branch below mutates ``theta`` in
    # place. The caller reuses the same array across every mutator
    # call in soc_apbc_topology_guard's loop, so an in-place mutation
    # on an aliased view would silently corrupt subsequent mutations.
    theta = np.array(boundary_theta, dtype=np.float64, copy=True)
    twist_L = L_full.copy()
    twist_dr = dr_full.copy()
    twist_sign = -1.0

    base_id, axis = _split_mutator_id(mutator)
    if base_id not in (
        "baseline", "M-gauge-1", "M-gauge-2", "M-gauge-3", "M-gauge-4",
        "M-gauge-5",
    ):
        raise ValueError(f"_gauge_lift_element: unknown mutator {mutator!r}")

    if axis is None:
        # Whole-vector mutators and baseline use the full-vector formula.
        if base_id == "M-gauge-1":
            twist_sign = +1.0
        elif base_id == "M-gauge-2":
            theta = theta / (2.0 * np.pi)
        elif base_id == "M-gauge-3":
            twist_L = L_folded.astype(np.float64)
        elif base_id == "M-gauge-4":
            so_diff = so_j.astype(np.float64) - so_i.astype(np.float64)
            dr_folded = (
                fc_j.astype(np.float64) - fc_i.astype(np.float64) - so_diff
            )
        elif base_id == "M-gauge-5":
            twist_dr = dr_folded_base
    else:
        # Per-direction mutators splice only component
        # ``axis`` into an otherwise-baseline vector; the other two axes
        # keep their baseline (unmutated) value. M-gauge-4-{x,y,z}
        # deliberately omits sub_offset on the named axis. This differs
        # from the whole-vector M-gauge-4 sign flip above: those are
        # distinct complete mutation ids with schema-frozen semantics.
        if base_id == "M-gauge-1":
            theta[axis] = -theta[axis]
        elif base_id == "M-gauge-2":
            theta[axis] = theta[axis] / (2.0 * np.pi)
        elif base_id == "M-gauge-3":
            twist_L[axis] = L_folded.astype(np.float64)[axis]
        elif base_id == "M-gauge-4":
            dr_folded[axis] = (
                fc_j.astype(np.float64)[axis]
                - fc_i.astype(np.float64)[axis]
            )
        elif base_id == "M-gauge-5":
            twist_dr[axis] = dr_folded_base[axis]

    phase_twist = np.exp(twist_sign * 1j * np.dot(theta, twist_dr / twist_L))
    accum = 0.0j
    for kx in range(LFX):
        for ky in range(LFY):
            for kz in range(LFZ):
                k_vec = 2.0 * np.pi * np.array(
                    [kx, ky, kz], dtype=np.float64
                ) / L_folded
                phase_folded = np.exp(-1j * np.dot(k_vec, dr_folded))
                accum += G_k[kx, ky, kz, aa, bb] * phase_folded
    accum *= phase_twist
    return accum / float(LFX * LFY * LFZ)


def _build_A_ship_mutated(mutator, eigen, occ, site_positions, cell_shape,
                          subshape, ncond, boundary_theta):
    """Real M-ship mutation: rebuild the shipping A matrix with the
    mutation applied to the SOC branch of build_slater_orbitals so the
    G4 gate exercises the mutation formulas, not a
    closed-form surrogate that can disagree on fixture-specific phase
    cancellations."""
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows
    from tools._uhfk_to_mvmc.sublattice_unfold import decode_physical_site

    L_folded = (cell_shape // subshape).astype(np.float64)
    L_phys = cell_shape.astype(np.float64)
    theta = np.asarray(boundary_theta, dtype=np.float64)
    Ns = site_positions.shape[0]

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[int(ncond)], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )
    wvi = eigen["wavevector_index"].astype(np.int64)
    ev = eigen["eigenvector"].astype(np.complex128)
    nvol_folded = wvi.shape[0]
    inv_sqrt = 1.0 / np.sqrt(float(nvol_folded))

    folded_cell = np.empty((Ns, 3), dtype=np.int64)
    sub_offset = np.empty((Ns, 3), dtype=np.int64)
    folded_orb = np.empty(Ns, dtype=np.int64)
    for i in range(Ns):
        fc, so = decode_physical_site(site_positions[i], subshape)
        folded_cell[i] = fc
        sub_offset[i] = so
        folded_orb[i] = (
            so[0] + subshape[0] * (so[1] + subshape[1] * so[2])
        )

    base_id, axis = _split_mutator_id(mutator)
    if base_id not in (
        "baseline", "M-ship-1", "M-ship-2", "M-ship-3", "M-ship-4",
        "M-ship-5",
    ):
        raise ValueError(
            f"_build_A_ship_mutated: unknown mutator {mutator!r}"
        )

    if axis is None:
        # Whole-vector mutators, baseline, and M-ship-4/5 reuse phys_dn.
        if base_id == "M-ship-1":
            phys_arg = np.einsum(
                "d,id->i", theta / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(+1j * phys_arg)  # sign flipped
        elif base_id == "M-ship-2":
            theta_wrong = theta / (2.0 * np.pi)
            phys_arg = np.einsum(
                "d,id->i", theta_wrong / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        elif base_id == "M-ship-3":
            phys_arg = np.einsum(
                "d,id->i", theta / L_folded,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
        else:
            phys_arg = np.einsum(
                "d,id->i", theta / L_phys,
                site_positions.astype(np.float64),
            )
            phys_dn = np.exp(-1j * phys_arg)
    else:
        # Per-direction mutators apply the same formula restricted to
        # component ``axis``; the other
        # two axes use the baseline theta/L_phys value. M-ship-4/-5
        # leave phys_dn at its baseline value (they only touch
        # kf_dot_r below), matching the whole-vector else-branch.
        theta_v37 = theta.copy()
        L_phys_v37 = L_phys.copy()
        if base_id == "M-ship-1":
            theta_v37[axis] = -theta_v37[axis]
        elif base_id == "M-ship-2":
            theta_v37[axis] = theta_v37[axis] / (2.0 * np.pi)
        elif base_id == "M-ship-3":
            L_phys_v37[axis] = L_folded[axis]
        phys_arg = np.einsum(
            "d,id->i", theta_v37 / L_phys_v37,
            site_positions.astype(np.float64),
        )
        phys_dn = np.exp(-1j * phys_arg)

    k_folded_all = 2.0 * np.pi * wvi.astype(np.float64) / L_folded
    if axis is None:
        # Whole-vector mutators, baseline, and M-ship-1/2/3 reuse kf_dot_r.
        if base_id == "M-ship-4":
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all,
                (folded_cell - sub_offset).astype(np.float64),
            )
        elif base_id == "M-ship-5":
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all, folded_cell.astype(np.float64),
            )
        else:
            kf_dot_r = np.einsum(
                "kd,id->ki", k_folded_all,
                (folded_cell + sub_offset).astype(np.float64),
            )
    else:
        # Per-direction mutators keep the other two axes
        # keep the baseline (folded_cell + sub_offset) value.
        # M-ship-4-{x,y,z} deliberately halves sub_offset on the named
        # axis, while M-ship-5 omits it. Coefficients 1 (baseline), 1/2
        # (M-4), and 0 (M-5) are distinct by construction. This differs
        # from the whole-vector M-ship-4 sign flip above: those are
        # distinct complete mutation ids with schema-frozen semantics.
        # M-ship-1/-2/-3 leave kf_dot_r at its baseline value because
        # they only touch phys_dn above.
        pos_v37 = (folded_cell + sub_offset).astype(np.float64)
        if base_id == "M-ship-4":
            pos_v37[:, axis] = (
                folded_cell[:, axis].astype(np.float64)
                + 0.5 * sub_offset[:, axis].astype(np.float64)
            )
        elif base_id == "M-ship-5":
            pos_v37[:, axis] = folded_cell[:, axis].astype(np.float64)
        kf_dot_r = np.einsum("kd,id->ki", k_folded_all, pos_v37)

    plane_wave = (
        np.exp(-1j * kf_dot_r) * inv_sqrt * phys_dn[np.newaxis, :]
    )
    A = np.zeros(
        (2 * Ns, 2 * len(pair_list)), dtype=np.complex128,
    )
    for p, pair in enumerate(pair_list):
        for col_idx, member in enumerate(pair):
            k_row, alpha = int(member[0]), int(member[1])
            slater_col = 2 * p + col_idx
            pw_m = plane_wave[k_row]
            for spin in (0, 1):
                hw = (2 * folded_orb + spin).astype(np.int64)
                v = np.array(
                    [ev[k_row, int(hw[i]), alpha] for i in range(Ns)],
                    dtype=np.complex128,
                )
                A[spin * Ns:(spin + 1) * Ns, slater_col] = pw_m * v
    return A


def _resolve_workspace(workspace):
    """Canonicalize the workspace path and every consumed
    artifact; reject any resolved artifact under ``tests/data``. The
    shared guard covers the workspace root, subdir-symlink cases, and
    single-artifact-symlink cases."""
    from pathlib import Path

    _scripts_dir = os.path.dirname(os.path.abspath(__file__))
    if _scripts_dir not in sys.path:
        sys.path.insert(0, _scripts_dir)
    from _snapshot_guard import (
        reject_snapshot_workspace,
        SnapshotWorkspaceRejected,
    )

    try:
        reject_snapshot_workspace(
            workspace, helper_name="soc_apbc_topology_guard",
        )
    except SnapshotWorkspaceRejected as exc:
        raise ValueError(str(exc)) from exc
    ws = Path(workspace).resolve(strict=True)
    return str(ws)


def soc_apbc_topology_guard(workspace, composite_manifest_path):
    """Run the semantic composite and mutation check on the current-run
    workspace. Returns a dict with ``max_abs_delta`` and ``tol`` fields
    for the PASS record; raises on any semantic failure so the CLI
    entry point can convert to exit code 2 and empty stdout."""
    _resolve_workspace(workspace)
    manifest = _load_manifest(composite_manifest_path)
    mutation_ids = _validate_manifest_numeric_policy(manifest)
    i_c = _require_finite_int(manifest["i_c"], "i_c")
    s_c = _require_finite_int(manifest["s_c"], "s_c")
    j_c = _require_finite_int(manifest["j_c"], "j_c")
    t_c = _require_finite_int(manifest["t_c"], "t_c")
    G_c_abs_manifest = _require_finite_float(
        manifest["G_c_abs"], "G_c_abs",
    )
    T_M_manifest = manifest["T_M_per_mutation"]

    # Cross-spin invariant (from manifest — defensively re-checked).
    if s_c == t_c:
        raise ValueError(
            f"composite manifest requires cross-spin; got s_c=t_c={s_c}"
        )

    cell_shape = _require_finite_vector(
        manifest, "cell_shape", integral=True,
    )
    subshape = _require_finite_vector(
        manifest, "sub_shape", integral=True,
    )
    theta = _require_finite_vector(
        manifest, "theta_radians", integral=False,
    )
    Nsite = int(np.prod(cell_shape))

    # sub_offset differs along every ACTIVE APBC direction (from
    # manifest — defensively re-checked). Every active direction must have
    # differing offsets. Active axes are derived from theta directly so the
    # check also covers manifests without active_apbc_axes. See
    # ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    so_i_vec = [int(x) for x in manifest["sub_offset_i"]]
    so_j_vec = [int(x) for x in manifest["sub_offset_j"]]
    for _d in range(3):
        if abs(float(theta[_d])) <= 1e-12:
            continue  # Inactive axes are not gated.
        if so_i_vec[_d] == so_j_vec[_d]:
            _axis_char = "xyz"[_d]
            raise ValueError(
                f"composite manifest requires sub_offset_{_axis_char} "
                f"differs (active APBC axis {_axis_char}, "
                f"theta_{_axis_char}={float(theta[_d]):.6f}); got "
                f"so_i_{_axis_char}=so_j_{_axis_char}={so_i_vec[_d]}"
            )

    # Magnitude floor from manifest.
    if G_c_abs_manifest < 1e-3:
        raise ValueError(
            f"composite manifest G_c_abs = {G_c_abs_manifest:.3e} is below "
            "the magnitude floor of 1e-3."
        )

    # Load CURRENT-run green_sublattice — NEVER trust manifest values.
    green_npz_path = os.path.join(workspace, "hwave", "green.npz")
    if not os.path.isfile(green_npz_path):
        # Try the SCF standard layout (case_dir/output/green.npz).
        alt = os.path.join(workspace, "output", "green.npz")
        if not os.path.isfile(alt):
            raise FileNotFoundError(
                f"soc_apbc_topology_guard: missing {green_npz_path} "
                f"and {alt}"
            )
        green_npz_path = alt
    green_sublattice = np.load(green_npz_path)["green_sublattice"]

    site_positions = np.array(
        [
            (ix, iy, iz)
            for iz in range(int(cell_shape[2]))
            for iy in range(int(cell_shape[1]))
            for ix in range(int(cell_shape[0]))
        ],
        dtype=np.int64,
    )

    G_current = _gauge_lift_element(
        green_sublattice, i_c, s_c, j_c, t_c,
        subshape=subshape, cell_shape=cell_shape,
        site_positions=site_positions, boundary_theta=theta,
        mutator="baseline",
    )
    G_current_abs = _require_finite_float(
        abs(G_current), "current-workspace G_current_abs",
    )
    magnitude_ratio = G_current_abs / max(G_c_abs_manifest, 1e-16)
    if G_current_abs < 0.8 * G_c_abs_manifest:
        raise ValueError(
            f"composite element on current workspace: "
            f"|G_current[i_c,s_c,j_c,t_c]| = {G_current_abs:.3e} < "
            f"0.8 * G_c_abs_manifest = {0.8 * G_c_abs_manifest:.3e} "
            f"(magnitude_ratio = {magnitude_ratio:.3f}); composite "
            "did not survive; regenerate manifest with a fresh SCF run."
        )

    # M-gauge mutations run at element level on the current-run
    # green_sublattice.
    max_shortfall = 0.0
    for m_gauge in (m for m in mutation_ids if m.startswith("M-gauge-")):
        G_mut = _gauge_lift_element(
            green_sublattice, i_c, s_c, j_c, t_c,
            subshape=subshape, cell_shape=cell_shape,
            site_positions=site_positions, boundary_theta=theta,
            mutator=m_gauge,
        )
        delta = _require_finite_float(
            abs(G_mut - G_current),
            f"current-workspace delta[{m_gauge!r}]",
        )
        T = _require_finite_float(
            T_M_manifest[m_gauge],
            f"T_M_per_mutation[{m_gauge!r}]",
        )
        if delta < T:
            raise ValueError(
                f"{m_gauge}: delta = {delta:.3e} < T_M = {T:.3e} "
                "on current workspace; composite lost mutation "
                "sensitivity vs manifest."
            )
        shortfall = T / max(delta, 1e-16)
        max_shortfall = max(max_shortfall, shortfall)

    # M-ship mutations require the shipping A on the current run. Rebuild
    # from eigen/occupation NPZ files under the workspace (SCF standard
    # layout ``output/`` or bridge-linked ``hwave/``); the baseline A and
    # 5 mutated A matrices are constructed once and the composite element
    # extracted from each.
    eigen_path = None
    occ_path = None
    for base in ("hwave", "output"):
        e = os.path.join(workspace, base, "eigen.npz")
        o = os.path.join(workspace, base, "occupation.npz")
        if os.path.isfile(e) and os.path.isfile(o):
            eigen_path, occ_path = e, o
            break
    if eigen_path is None:
        raise FileNotFoundError(
            "soc_apbc_topology_guard: missing eigen.npz + occupation.npz "
            "under workspace's hwave/ or output/ directory; cannot run "
            "M-ship mutations."
        )
    n_pairs = _require_finite_int(manifest["N_pairs"], "N_pairs")
    ncond = _require_finite_int(
        manifest.get("ncond", 2 * n_pairs), "ncond",
    )
    eigen = np.load(eigen_path)
    occ = np.load(occ_path)
    A_base = _build_A_ship_mutated(
        "baseline", eigen, occ, site_positions, cell_shape, subshape,
        ncond, theta,
    )

    # Verify that the shadow-copy baseline equals the shipping A produced
    # by build_slater_orbitals before running the mutation matrix. Otherwise,
    # drift between the shadow copy and canonical kernel would be invisible.
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, build_slater_orbitals, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows
    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[ncond], is_soc_mode=True,
    )
    partner_rows_ship, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows_ship, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows_ship,
        is_soc_mode=True,
    )
    A_ship_canonical = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions,
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    shadow_drift = _require_finite_float(
        np.max(np.abs(A_base - A_ship_canonical)),
        "current-workspace shadow_drift",
    )
    if shadow_drift > 1e-10:
        raise ValueError(
            f"G4 shadow-copy baseline drifted from shipping "
            f"build_slater_orbitals: |A_shadow - A_ship|_max = "
            f"{shadow_drift:.3e} > 1e-10. Either the shadow copy or the "
            "shipping kernel has changed; re-align "
            "_build_A_ship_mutated with the shipping SOC branch or "
            "extract a shared kernel."
        )
    G_ship_base = np.conj(A_base) @ A_base.T
    for m_ship in (m for m in mutation_ids if m.startswith("M-ship-")):
        A_m = _build_A_ship_mutated(
            m_ship, eigen, occ, site_positions, cell_shape, subshape,
            ncond, theta,
        )
        G_ship_m = np.conj(A_m) @ A_m.T
        all_i = i_c + s_c * Nsite
        all_j = j_c + t_c * Nsite
        delta = _require_finite_float(
            abs(G_ship_m[all_i, all_j] - G_ship_base[all_i, all_j]),
            f"current-workspace delta[{m_ship!r}]",
        )
        T = _require_finite_float(
            T_M_manifest[m_ship],
            f"T_M_per_mutation[{m_ship!r}]",
        )
        if delta < T:
            raise ValueError(
                f"{m_ship}: delta = {delta:.3e} < T_M = {T:.3e} "
                "on current workspace."
            )

    return {
        "max_abs_delta": max(0.0, 1.0 - magnitude_ratio),
        "tol": 0.2,  # 20% magnitude drift allowed
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "G4 topology and mutation guard. The semantic guard verifies the pinned "
            "composite element on the current-run workspace and runs "
            "the manifest's mutation matrix (10 whole-vector "
            "mutations, or 30 per-direction mutations; schema "
            "auto-detected from the manifest)."
        )
    )
    parser.add_argument("--workspace", required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument(
        "--composite-manifest", dest="composite_manifest",
        required=True,
    )
    args = parser.parse_args(argv)

    if args.mode != G4_MODE:
        print(
            f"soc_apbc_topology_guard: unsupported mode={args.mode}; "
            f"only {G4_MODE!r} is defined.",
            file=sys.stderr,
        )
        return 2

    try:
        result = soc_apbc_topology_guard(
            args.workspace, args.composite_manifest,
        )
    except (FileNotFoundError, ValueError, KeyError) as e:
        print(f"soc_apbc_topology_guard: {e}", file=sys.stderr)
        return 2

    max_abs_delta = float(result.get("max_abs_delta", 0.0))
    tol = float(result.get("tol", 0.2))
    print(
        f"G4 PASS mode={G4_MODE} artifact_source={G4_ARTIFACT_SOURCE} "
        f"helper={G4_HELPER} max_abs_delta={max_abs_delta:.6e} "
        f"tol={tol:.6e}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
