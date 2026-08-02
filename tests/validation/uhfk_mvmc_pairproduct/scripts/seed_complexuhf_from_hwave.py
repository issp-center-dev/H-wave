"""Seed ComplexUHF UHF with H-wave's converged density.

ComplexUHF's random initialization + independent SCF can land on a
different broken-symmetry minimum than H-wave (both are valid mean-field
minima; U > 0 Hubbard + SOC + APBC has multiple degenerate ones). This
divergence makes element-level G2a / G2b density comparisons meaningless
even when both solvers converged to correct physics.

The Initial keyword in ComplexUHF's namelist.def (mVMC-1.4.0 src/
ComplexUHF/readdef.c KWInitial + initial.c) supplies the full 2Ns x 2Ns
one-body Green as the SCF starting point:

    ====================
    NInitial <N>
    ====================
    ====================
    site_0 spin_0 site_1 spin_1 re im
    ...

Given the same Hamiltonian + a starting density near a specific minimum,
ComplexUHF converges to that minimum rather than a random one; the
resulting zvo_UHF_cisajs.dat then matches H-wave's density at the SCF
convergence tolerance rather than at "which basin did we land in"
tolerance.

This helper takes an H-wave workspace + ComplexUHF workspace, builds the
shipping A density via build_slater_orbitals, writes initial.def, and
patches namelist.def to reference it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

import numpy as np


def _find_workspace_config(hwave_workspace):
    """Load geometry and SOC configuration from a validation workspace."""
    _here = os.path.dirname(os.path.abspath(__file__))
    repo = os.path.abspath(os.path.join(_here, "..", "..", "..", ".."))
    if repo not in sys.path:
        sys.path.insert(0, repo)
    src = os.path.join(repo, "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from tools._uhfk_to_mvmc.boundary_input import (
        normalize_boundary_condition_list,
    )
    from tools._uhfk_to_mvmc.input_loader import (
        load_geometry_uhf, load_input_toml,
    )

    input_toml = os.path.join(hwave_workspace, "input.toml")
    geometry = os.path.join(hwave_workspace, "geometry_uhf.dat")
    toml_param = load_input_toml(input_toml)
    cell_shape = np.asarray(toml_param["CellShape"], dtype=np.int64)
    sub_shape = np.asarray(
        toml_param.get("SubShape", list(cell_shape)), dtype=np.int64,
    )
    boundary_theta = np.asarray(
        normalize_boundary_condition_list(toml_param["BoundaryCondition"]),
        dtype=np.float64,
    )
    Ncond = int(toml_param["Ncond"])
    is_soc_mode = bool(toml_param.get("enable_spin_orbital", False))
    unit_vec, site_R_int, norb = load_geometry_uhf(geometry)
    site_positions = site_R_int.astype(np.int64)
    Ns = site_positions.shape[0]
    return {
        "cell_shape": cell_shape,
        "sub_shape": sub_shape,
        "boundary_theta": boundary_theta,
        "Ncond": Ncond,
        "is_soc_mode": is_soc_mode,
        "site_positions": site_positions,
        "Ns": Ns,
    }


def _build_shipping_G(hwave_workspace, cfg):
    """Build G = conj(A_ship) @ A_ship.T on the physical basis."""
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, build_slater_orbitals, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    eigen = np.load(os.path.join(hwave_workspace, "output", "eigen.npz"))
    occ = np.load(
        os.path.join(hwave_workspace, "output", "occupation.npz"),
    )
    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[cfg["Ncond"]], is_soc_mode=cfg["is_soc_mode"],
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], cfg["boundary_theta"],
        cfg["cell_shape"] // cfg["sub_shape"],
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=cfg["is_soc_mode"],
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=cfg["site_positions"],
        cell_shape=cfg["cell_shape"],
        subshape=cfg["sub_shape"],
        theta=cfg["boundary_theta"],
        pair_list=pair_list,
        is_soc_mode=cfg["is_soc_mode"],
    )
    return np.conj(A_ship) @ A_ship.T


def _write_initial_def(
    G_full, out_path, mag_threshold=1e-12,
    perturb_scale=0.0, perturb_seed=0,
):
    """Emit initial.def with every non-negligible G[i, s, j, t] element.

    Skips machine-precision-zero elements to keep the file compact; the
    ComplexUHF loader defaults missing entries to 0 anyway.

    When ``perturb_scale > 0``, add a small random Hermitian perturbation
    of magnitude
    ``perturb_scale * |G_full|_max`` to the initial density before
    writing. Otherwise ComplexUHF is seeded at the exact SCF fixed
    point of the same Hamiltonian, terminates at 0 steps, and the
    resulting cross-solver comparison degenerates to reading-back the
    supplied seed. With a non-zero perturbation the SCF must actually
    iterate; convergence back to the same basin at 1e-6 tolerance
    against H-wave's density then provides genuine independent
    verification of the shared Hamiltonian + SCF-update path."""
    Ns = G_full.shape[0] // 2
    G_seed = np.array(G_full, copy=True)
    if perturb_scale > 0:
        # Deterministic reproducible perturbation. Hermitian noise so
        # ComplexUHF's density-matrix invariants (real diagonal, G = G^H)
        # remain plausible starting points.
        rng = np.random.default_rng(perturb_seed)
        scale = perturb_scale * float(np.max(np.abs(G_full)))
        noise = (
            rng.standard_normal(G_seed.shape)
            + 1j * rng.standard_normal(G_seed.shape)
        ) * scale
        noise = 0.5 * (noise + noise.conj().T)
        G_seed = G_seed + noise
    rows = []
    for a in range(2 * Ns):
        i = a % Ns
        s = a // Ns
        for b in range(2 * Ns):
            j = b % Ns
            t = b // Ns
            v = G_seed[a, b]
            if abs(v) < mag_threshold:
                continue
            rows.append((i, s, j, t, float(v.real), float(v.imag)))
    with open(out_path, "w") as fp:
        # ComplexUHF ReadDefFile skips IgnoreLinesInDef (=5) header
        # lines on the data-read pass; write 5 header lines exactly.
        fp.write("=============================================\n")
        fp.write(f"NInitial          {len(rows)}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")
        for (i, s, j, t, re, im) in rows:
            fp.write(f"{i:6d} {s:6d} {j:6d} {t:6d} "
                     f"{re: 24.16e} {im: 24.16e}\n")


def _write_provenance(initial_path, hwave_workspace, perturb_scale):
    """Write the canonical-JSON provenance sidecar."""
    initial_path = Path(initial_path)
    hwave_output = Path(hwave_workspace) / "output"

    def sha256_file(path):
        digest = hashlib.sha256()
        with open(path, "rb") as fp:
            for chunk in iter(lambda: fp.read(65536), b""):
                digest.update(chunk)
        return digest.hexdigest()

    energy_total = None
    energy_path = hwave_output / "energy.dat"
    with open(energy_path) as fp:
        for line in fp:
            if re.match(r"\s*Energy_Total\s*=", line):
                energy_total = float(line.split("=", 1)[1].strip())
                break
    if energy_total is None:
        raise ValueError(f"{energy_path}: Energy_Total entry missing")

    payload = {
        "sha256_initial_def": sha256_file(initial_path),
        "sha256_hwave_green": sha256_file(hwave_output / "green.npz"),
        "sha256_hwave_eigen": sha256_file(hwave_output / "eigen.npz"),
        "energy_total": energy_total,
        "perturb_scale": perturb_scale,
    }
    sidecar_path = initial_path.with_name("initial.def.provenance")
    sidecar_path.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")),
        encoding="utf-8",
    )
    return sidecar_path


def _patch_namelist(namelist_path, initial_def_relpath):
    """Append 'Initial <path>' entry to namelist.def if not present."""
    with open(namelist_path) as fp:
        content = fp.read()
    if "Initial " in content:
        return
    lines = content.rstrip().splitlines()
    lines.append(f"Initial {initial_def_relpath}")
    with open(namelist_path, "w") as fp:
        fp.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--hwave-workspace", required=True,
                    help="Path to ${WORK}/hwave/ with output/{eigen,"
                    "occupation}.npz, input.toml, geometry_uhf.dat")
    ap.add_argument("--complexuhf-workspace", required=True,
                    help="Path to ${WORK}/complexuhf/ with namelist.def")
    ap.add_argument("--initial-def", default="initial.def",
                    help="Name of the emitted initial.def under the "
                    "ComplexUHF workspace (default: initial.def)")
    ap.add_argument("--perturb-scale", type=float, default=1e-3,
                    help="Hermitian random perturbation applied to the "
                    "H-wave density before writing. A zero perturbation makes "
                    "ComplexUHF terminate at 0 SCF steps (circular "
                    "verification). Default 1e-3 forces real iteration "
                    "while remaining inside H-wave's basin. Pass 0 to "
                    "disable and reproduce the pre-hardening seed.")
    ap.add_argument("--perturb-seed", type=int, default=42,
                    help="RNG seed for the perturbation (default 42).")
    args = ap.parse_args()

    cfg = _find_workspace_config(args.hwave_workspace)
    G_ship = _build_shipping_G(args.hwave_workspace, cfg)

    initial_path = os.path.join(args.complexuhf_workspace, args.initial_def)
    _write_initial_def(
        G_ship, initial_path,
        perturb_scale=args.perturb_scale,
        perturb_seed=args.perturb_seed,
    )
    provenance_path = _write_provenance(
        initial_path, args.hwave_workspace, args.perturb_scale,
    )
    namelist_path = os.path.join(args.complexuhf_workspace, "namelist.def")
    _patch_namelist(namelist_path, args.initial_def)
    print(f"seed_complexuhf_from_hwave: wrote {initial_path} "
          f"(perturb_scale={args.perturb_scale})")
    print(f"seed_complexuhf_from_hwave: wrote {provenance_path}")
    print(f"seed_complexuhf_from_hwave: patched {namelist_path}")


if __name__ == "__main__":
    main()
