"""Seven-gate dispatcher and legacy energy-only compatibility.

Two invocation modes:

1. **Legacy**: three positional args
   ``compare.py <hwave_energy_dat> <mvmc_zvo_out_dat> <case_name>``.
   Reports H-wave-vs-mVMC energy delta with the same output the existing
   PBC / APBC fixtures rely on.

2. **Seven-gate dispatcher**:
   ``compare.py --workspace ${WORK_DIR} --mode {g0-writer-check|g1|
     g2a-emitted-F|g2a-in-memory|g2b|g3}``.
   Emits ONE anchored PASS record using the canonical helper registered
   in :data:`EXPECTED_MODE_DISPATCH`. See
   ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.

The dispatcher enforces:
  - EXPECTED_MODE_DISPATCH is a deep-frozen (MappingProxyType +
    tuple-of-tuples) reference table.
  - MODE_DISPATCH is a mutable copy compared against EXPECTED_MODE_DISPATCH
    before every mode's execution to catch registry corruption.
  - Every canonical helper resolves via ``importlib.import_module`` plus
    ``getattr``; the resolved object's ``__module__`` and ``__qualname__``
    MUST match the table (blocks a same-qualname-different-module
    hostile substitution).
  - Any dispatch integrity failure exits code 2 BEFORE any PASS line is
    printed.
"""
from __future__ import annotations

import argparse
import importlib
import os
import re
import sys
from statistics import mean, stdev
from types import MappingProxyType
from typing import List


# ---------------------------------------------------------------------
# Legacy energy-only helpers.
# ---------------------------------------------------------------------


def parse_hwave_energy(path: str) -> float:
    """Parse H-wave UHFk's ``energy.dat`` and return Energy_Total."""
    with open(path) as fp:
        for line in fp:
            m = re.match(
                r"\s*Energy_Total\s*=\s*(-?\d+\.\d+(?:[eE][+-]?\d+)?)",
                line,
            )
            if m:
                return float(m.group(1))
    raise RuntimeError(f"could not find Energy_Total in {path}")


def parse_mvmc_zvo_out(path: str) -> List[float]:
    """Return the list of per-bin <H> values from a zvo_out_*.dat file."""
    samples: List[float] = []
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if not toks:
                continue
            try:
                samples.append(float(toks[0]))
            except ValueError:
                continue
    if not samples:
        raise RuntimeError(f"no numeric rows parsed from {path}")
    return samples


# ---------------------------------------------------------------------
# Single source of truth for gate metadata.
# ---------------------------------------------------------------------


# Deep-freeze via tuples + MappingProxyType per entry so no in-place
# mutation is possible via slicing / append or nested dict assignment.
EXPECTED_MODE_DISPATCH = MappingProxyType({
    "g0-writer-check": MappingProxyType({
        "helpers": (
            ("tools._uhfk_to_mvmc.pair_product_density", "parse_emitted_F"),
        ),
        "artifact_source": "bridge_zeronoise",
    }),
    "g1": MappingProxyType({
        "helpers": (
            ("tools._uhfk_to_mvmc.general_fij_builder", "build_slater_orbitals"),
            ("tools._uhfk_to_mvmc.density_check",       "gauge_lift"),
        ),
        "artifact_source": "bridge+hwave",
    }),
    "g2a-emitted-F": MappingProxyType({
        "helpers": (
            ("tools._uhfk_to_mvmc.pair_product_density", "parse_emitted_F"),
            ("tools._uhfk_to_mvmc.pair_product_density",
             "pair_product_density_from_F"),
        ),
        "artifact_source": "bridge+complexuhf",
    }),
    "g2a-in-memory": MappingProxyType({
        "helpers": (
            ("tools._uhfk_to_mvmc.general_fij_builder", "build_slater_orbitals"),
        ),
        "artifact_source": "hwave+complexuhf",
    }),
    "g2b": MappingProxyType({
        "helpers": (("tools._uhfk_to_mvmc.density_check", "gauge_lift"),),
        "artifact_source": "hwave+complexuhf",
    }),
    "g3": MappingProxyType({
        "helpers": (
            ("tools._uhfk_to_mvmc.energy_compare", "energy_relative_delta"),
        ),
        "artifact_source": "hwave+mvmc",
    }),
})


# Runtime dispatch table; wiring tests corrupt this to verify
# _validate_dispatch catches registry drift.
MODE_DISPATCH = {
    mode: {
        "helpers": list(entry["helpers"]),
        "artifact_source": entry["artifact_source"],
    }
    for mode, entry in EXPECTED_MODE_DISPATCH.items()
}


class DispatchError(Exception):
    """Raised when the mode/helper registry fails integrity checks."""


def _validate_dispatch(requested_mode: str):
    """Validate MODE_DISPATCH against EXPECTED_MODE_DISPATCH and resolve
    the canonical helpers for ``requested_mode``. Exits code 2 on any
    integrity failure, BEFORE any PASS line can be printed."""
    if requested_mode not in EXPECTED_MODE_DISPATCH:
        print(
            f"MODE_DISPATCH: unknown requested_mode={requested_mode}",
            file=sys.stderr,
        )
        sys.exit(2)
    # Structural comparison against the deep-frozen reference table.
    for expected_mode, expected_entry in EXPECTED_MODE_DISPATCH.items():
        runtime_entry = MODE_DISPATCH.get(expected_mode)
        if runtime_entry is None:
            print(
                f"MODE_DISPATCH corrupted: missing mode={expected_mode}",
                file=sys.stderr,
            )
            sys.exit(2)
        if tuple(runtime_entry["helpers"]) != expected_entry["helpers"]:
            print(
                f"MODE_DISPATCH corrupted: helpers mismatch for "
                f"mode={expected_mode}",
                file=sys.stderr,
            )
            sys.exit(2)
        if runtime_entry["artifact_source"] != expected_entry["artifact_source"]:
            print(
                f"MODE_DISPATCH corrupted: artifact_source mismatch for "
                f"mode={expected_mode}",
                file=sys.stderr,
            )
            sys.exit(2)
    entry = MODE_DISPATCH[requested_mode]
    resolved = []
    for module_dotted, qualname in entry["helpers"]:
        try:
            mod = importlib.import_module(module_dotted)
        except ImportError as e:
            print(
                f"MODE_DISPATCH: cannot import {module_dotted} for "
                f"mode={requested_mode}: {e}",
                file=sys.stderr,
            )
            sys.exit(2)
        obj = getattr(mod, qualname, None)
        if obj is None:
            print(
                f"MODE_DISPATCH: {module_dotted} has no attribute "
                f"{qualname} for mode={requested_mode}",
                file=sys.stderr,
            )
            sys.exit(2)
        if getattr(obj, "__module__", None) != module_dotted:
            print(
                f"MODE_DISPATCH: {qualname} resolved via {module_dotted} "
                f"but obj.__module__="
                f"{getattr(obj, '__module__', None)}; wrong-module dispatch",
                file=sys.stderr,
            )
            sys.exit(2)
        if obj.__qualname__ != qualname:
            print(
                f"MODE_DISPATCH: expected qualname={qualname}, got "
                f"{obj.__qualname__}",
                file=sys.stderr,
            )
            sys.exit(2)
        resolved.append(obj)
    return entry, resolved


def _emit_pass(
    mode: str,
    max_abs_delta: float,
    tol: float,
    *,
    g2_evidence=None,
):
    """Print the anchored PASS record. All metadata fields come
    from EXPECTED_MODE_DISPATCH (never MODE_DISPATCH) so a runtime table
    mutation cannot poison the emitted PASS."""
    expected = EXPECTED_MODE_DISPATCH[mode]
    helpers = "+".join(qn for _, qn in expected["helpers"])
    gate_name = _mode_to_gate_name(mode)
    record = (
        f"{gate_name} PASS mode={mode} "
        f"artifact_source={expected['artifact_source']} "
        f"helper={helpers} max_abs_delta={max_abs_delta:.6e} "
        f"tol={tol:.6e}"
    )
    if g2_evidence is not None:
        initial_delta, final_delta, contraction_ratio = g2_evidence
        record += (
            f" initial_delta={initial_delta:.6e} "
            f"final_delta={final_delta:.6e} "
            f"contraction_ratio={contraction_ratio:.6e}"
        )
    print(record)


def _mode_to_gate_name(mode: str) -> str:
    return {
        "g0-writer-check": "G0-writer-check",
        "g1": "G1",
        "g2a-emitted-F": "G2a-emitted-F",
        "g2a-in-memory": "G2a-in-memory-A",
        "g2b": "G2b",
        "g3": "G3",
    }[mode]


# ---------------------------------------------------------------------
# Workspace config loader (shared by G1 / G2a-* / G2b).
# ---------------------------------------------------------------------


def _load_workspace_config(workspace):
    """Return the geometry / SOC config the G1/G2* dispatchers need.

    Reads workspace/hwave/input.toml + workspace/hwave/geometry_uhf.dat.
    Under the seven-gate producer contract run.sh copies both artefacts
    into ${WORK_DIR}/hwave/; missing files raise FileNotFoundError with
    the canonical path so a hostile stripped-down workspace cannot
    silently trigger a False-neutral PASS.
    """
    # Lazy import to keep compare.py importable without the tools tree
    # on sys.path (used by non-dispatch unit tests).
    from tools._uhfk_to_mvmc.boundary_input import (
        normalize_boundary_condition_list,
    )
    from tools._uhfk_to_mvmc.input_loader import (
        load_geometry_uhf, load_input_toml,
    )

    import numpy as np

    input_toml = os.path.join(workspace, "hwave", "input.toml")
    geometry = os.path.join(workspace, "hwave", "geometry_uhf.dat")
    if not os.path.isfile(input_toml):
        raise FileNotFoundError(input_toml)
    if not os.path.isfile(geometry):
        raise FileNotFoundError(geometry)
    toml_param = load_input_toml(input_toml)
    cell_shape = np.asarray(toml_param["CellShape"], dtype=np.int64)
    sub_shape = np.asarray(
        toml_param.get("SubShape", list(cell_shape)), dtype=np.int64,
    )
    boundary_theta = np.asarray(
        normalize_boundary_condition_list(
            toml_param["BoundaryCondition"]
        ),
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


def _build_shipping_A(workspace, cfg):
    """Build the shipping A matrix from workspace/hwave/{eigen,occupation}.npz
    using the canonical build_slater_orbitals; returns (A, pair_list)."""
    import numpy as np
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_pair_list, build_slater_orbitals, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    eigen = np.load(os.path.join(workspace, "hwave", "eigen.npz"))
    occ = np.load(os.path.join(workspace, "hwave", "occupation.npz"))
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
    return A_ship, pair_list


def _build_G_from_gauge_lift(workspace, cfg, gauge_lift_callable):
    """Build the full 2Ns x 2Ns G matrix via the resolved gauge_lift.

    Uses the boundary_theta from workspace config, NOT (0, 0, 0). See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    import numpy as np

    green_path = os.path.join(workspace, "hwave", "green.npz")
    green_sub = np.load(green_path)["green_sublattice"]
    Ns = cfg["Ns"]
    G = np.zeros((2 * Ns, 2 * Ns), dtype=np.complex128)
    folded_cell_of = lambda r_phys: r_phys // cfg["sub_shape"]
    for i in range(Ns):
        for s in (0, 1):
            for j in range(Ns):
                for t in (0, 1):
                    G[i + s * Ns, j + t * Ns] = gauge_lift_callable(
                        green_sub, i, s, j, t,
                        subshape=cfg["sub_shape"],
                        cell_shape=cfg["cell_shape"],
                        site_positions=cfg["site_positions"],
                        folded_cell_of=folded_cell_of,
                        boundary_theta=cfg["boundary_theta"],
                    )
    return G


class ComplexUHFParseError(RuntimeError):
    """Raised when zvo_UHF_cisajs.dat fails a strict validation check.

    Short or unparsable lines and empty files must not silently produce
    zeros: ``G_cuhf = 0`` and ``G_bridge = 0`` would then agree trivially
    on an invalid ComplexUHF artifact.
    """


def _load_complexuhf_G(workspace, Ns):
    """Parse ComplexUHF's zvo_UHF_cisajs.dat into a full 2Ns x 2Ns G
    matrix (spin-block ordering: all_i = i + s * Ns).

    Missing file: FileNotFoundError. Empty / short / partial /
    malformed / duplicate / out-of-range / non-finite: raises
    ComplexUHFParseError. Every (i, s, j, t) with 0 <= i, j < Ns and
    s, t in {0, 1} MUST appear exactly once and MUST parse to a
    finite complex number. Missing coverage is a wrong artifact and
    fails closed."""
    import numpy as np

    path = os.path.join(workspace, "complexuhf", "zvo_UHF_cisajs.dat")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    G = np.zeros((2 * Ns, 2 * Ns), dtype=np.complex128)
    seen = set()
    with open(path) as fp:
        for lineno, ln in enumerate(fp, start=1):
            stripped = ln.strip()
            if not stripped:
                continue
            toks = stripped.split()
            if len(toks) != 6:
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: expected exactly 6 tokens "
                    f"(i s j t re im); got {len(toks)}: {stripped!r}"
                )
            try:
                i = int(toks[0]); s = int(toks[1])
                j = int(toks[2]); t = int(toks[3])
                re = float(toks[4]); im = float(toks[5])
            except ValueError as exc:
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: token parse error ({exc}): "
                    f"{stripped!r}"
                )
            if not (0 <= i < Ns) or not (0 <= j < Ns):
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: site index out of range "
                    f"[0, {Ns}): i={i}, j={j}"
                )
            if s not in (0, 1) or t not in (0, 1):
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: spin out of range {{0, 1}}: "
                    f"s={s}, t={t}"
                )
            key = (i, s, j, t)
            if key in seen:
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: duplicate (i, s, j, t) = {key}"
                )
            seen.add(key)
            val = complex(re, im)
            if not (np.isfinite(re) and np.isfinite(im)):
                raise ComplexUHFParseError(
                    f"{path}:{lineno}: non-finite value {val!r}"
                )
            G[i + s * Ns, j + t * Ns] = val
    expected = (2 * Ns) ** 2
    if len(seen) != expected:
        raise ComplexUHFParseError(
            f"{path}: expected {expected} (i, s, j, t) rows; got "
            f"{len(seen)}. Truncated or partial ComplexUHF artifact."
        )
    return G


class ComplexUHFInitialParseError(RuntimeError):
    """Raised when ComplexUHF's sparse ``initial.def`` is malformed."""


def _load_complexuhf_initial_G(workspace, Ns):
    """Parse ``initial.def`` into ComplexUHF's site-spin matrix order.

    ComplexUHF indexes the rows as ``2*i+s``. The five-line header is
    followed by sparse ``i s j t re im`` rows; omitted entries are zero.
    """
    import numpy as np

    path = os.path.join(workspace, "complexuhf", "initial.def")
    if not os.path.isfile(path):
        raise FileNotFoundError(path)
    with open(path) as fp:
        header = [fp.readline() for _ in range(5)]
        if any(line == "" for line in header):
            raise ComplexUHFInitialParseError(
                f"{path}: expected a five-line header"
            )
        count_match = re.fullmatch(r"\s*NInitial\s+(\d+)\s*", header[1])
        if count_match is None:
            raise ComplexUHFInitialParseError(
                f"{path}:2: expected 'NInitial <count>'"
            )
        expected_rows = int(count_match.group(1))
        G = np.zeros((2 * Ns, 2 * Ns), dtype=np.complex128)
        seen = set()
        for lineno, line in enumerate(fp, start=6):
            stripped = line.strip()
            if not stripped:
                continue
            tokens = stripped.split()
            if len(tokens) != 6:
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: expected exactly 6 tokens "
                    f"(i s j t re im); got {len(tokens)}"
                )
            try:
                i, s, j, t = (int(token) for token in tokens[:4])
                real, imag = (float(token) for token in tokens[4:])
            except ValueError as exc:
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: token parse error ({exc})"
                ) from exc
            if not (0 <= i < Ns) or not (0 <= j < Ns):
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: site index out of range "
                    f"[0, {Ns}): i={i}, j={j}"
                )
            if s not in (0, 1) or t not in (0, 1):
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: spin out of range {{0, 1}}: "
                    f"s={s}, t={t}"
                )
            key = (i, s, j, t)
            if key in seen:
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: duplicate (i, s, j, t) = {key}"
                )
            if not (np.isfinite(real) and np.isfinite(imag)):
                raise ComplexUHFInitialParseError(
                    f"{path}:{lineno}: non-finite value "
                    f"{complex(real, imag)!r}"
                )
            seen.add(key)
            G[2 * i + s, 2 * j + t] = complex(real, imag)
    if len(seen) != expected_rows:
        raise ComplexUHFInitialParseError(
            f"{path}: NInitial declares {expected_rows} rows; got "
            f"{len(seen)}"
        )
    return G


def _check_g2_contraction(mode, workspace, G_reference, G_final, tol):
    """Require a G2 trajectory from clearly outside to inside ``tol``."""
    import numpy as np

    Ns = G_reference.shape[0] // 2
    initial_site_spin = _load_complexuhf_initial_G(workspace, Ns)
    gate_name = _mode_to_gate_name(mode)
    finite_inputs = (
        ("tol", np.asarray(tol)),
        ("reference density", np.asarray(G_reference)),
        ("final density", np.asarray(G_final)),
        ("initial density", np.asarray(initial_site_spin)),
    )
    for label, values in finite_inputs:
        if not np.all(np.isfinite(values)):
            print(f"{gate_name}: {label} contains non-finite values", file=sys.stderr)
            return 1

    spin_block_indices = np.asarray(
        [i + s * Ns for i in range(Ns) for s in (0, 1)],
        dtype=np.int64,
    )
    reference_site_spin = G_reference[
        np.ix_(spin_block_indices, spin_block_indices)
    ]
    initial_delta = float(
        np.max(np.abs(initial_site_spin - reference_site_spin))
    )
    final_delta = float(np.max(np.abs(G_final - G_reference)))
    if not (np.isfinite(initial_delta) and np.isfinite(final_delta)):
        print(
            f"{gate_name}: non-finite contraction delta "
            f"(initial={initial_delta!r}, final={final_delta!r})",
            file=sys.stderr,
        )
        return 1

    with np.errstate(divide="ignore", invalid="ignore"):
        contraction_ratio = float(np.divide(final_delta, initial_delta))
    if not np.isfinite(contraction_ratio):
        print(
            f"{gate_name}: contraction ratio is non-finite: "
            f"{contraction_ratio!r}",
            file=sys.stderr,
        )
        return 1

    minimum_initial_delta = 10.0 * tol
    if initial_delta < minimum_initial_delta:
        print(
            f"{gate_name}: initial density delta = {initial_delta:.3e} "
            f"< required 10*tol = {minimum_initial_delta:.3e}",
            file=sys.stderr,
        )
        return 1
    if final_delta >= tol:
        print(
            f"{gate_name}: final density delta = {final_delta:.3e} "
            f">= tol = {tol:.3e}",
            file=sys.stderr,
        )
        return 1
    _emit_pass(
        mode,
        final_delta,
        tol,
        g2_evidence=(initial_delta, final_delta, contraction_ratio),
    )
    return 0


# ---------------------------------------------------------------------
# Mode dispatchers (real comparisons).
# ---------------------------------------------------------------------


def _validated_max_abs_delta(
    gate_name, lhs_name, lhs, rhs_name, rhs, gtol
):
    """Return a finite max delta for valid threshold-comparison inputs."""
    import numpy as np

    if not (np.isfinite(gtol) and gtol > 0.0):
        print(
            f"{gate_name}: tolerance must be finite and positive; "
            f"got {gtol!r}",
            file=sys.stderr,
        )
        return None
    for label, values in ((lhs_name, lhs), (rhs_name, rhs)):
        if not np.all(np.isfinite(values)):
            print(
                f"{gate_name}: {label} contains non-finite values",
                file=sys.stderr,
            )
            return None
    max_abs_delta = float(np.max(np.abs(lhs - rhs)))
    if not np.isfinite(max_abs_delta):
        print(
            f"{gate_name}: max absolute delta is non-finite "
            f"({max_abs_delta!r})",
            file=sys.stderr,
        )
        return None
    return max_abs_delta


def _dispatch_g0_writer_check(workspace: str, resolved, gtol: float) -> int:
    parse_emitted_F = resolved[0]
    subdir = os.path.join(workspace, "bridge_zeronoise")
    F_emitted = parse_emitted_F(subdir)
    import numpy as np
    F_pre = np.load(os.path.join(subdir, "F_pre_noise.npz"))["F"]
    max_abs_delta = _validated_max_abs_delta(
        "G0-writer-check",
        "emitted pair matrix",
        F_emitted,
        "pre-noise pair matrix",
        F_pre,
        gtol,
    )
    if max_abs_delta is None:
        return 1
    if max_abs_delta > gtol:
        return 1
    _emit_pass("g0-writer-check", max_abs_delta, gtol)
    return 0


def _dispatch_g1(workspace: str, resolved, gtol: float) -> int:
    """G1: shipping A density vs gauge_lift-lifted green_sublattice
    at ``gtol`` (1e-10 default). Both callables resolve via
    EXPECTED_MODE_DISPATCH; the composite helper string
    ``build_slater_orbitals+gauge_lift`` is bound to the PASS record. See
    ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    """
    import numpy as np
    build_slater_orbitals, gauge_lift = resolved

    cfg = _load_workspace_config(workspace)
    A_ship, _ = _build_shipping_A(workspace, cfg)
    G_ship = np.conj(A_ship) @ A_ship.T
    G_lift = _build_G_from_gauge_lift(workspace, cfg, gauge_lift)
    max_abs_delta = _validated_max_abs_delta(
        "G1",
        "shipping density",
        G_ship,
        "gauge-lifted density",
        G_lift,
        gtol,
    )
    if max_abs_delta is None:
        return 1
    if max_abs_delta > gtol:
        print(
            f"G1: |G_ship - G_lift|_max = {max_abs_delta:.3e} > "
            f"tol = {gtol:.3e}",
            file=sys.stderr,
        )
        return 1
    _emit_pass("g1", max_abs_delta, gtol)
    return 0


def _dispatch_g2a_emitted_F(workspace: str, resolved, gtol: float) -> int:
    """G2a-emitted-F: pair_product_density_from_F(parse_emitted_F(
    workspace/bridge)) vs ComplexUHF one-body Green at ``gtol`` (1e-6
    default). Fails closed on missing ComplexUHF file."""
    parse_emitted_F, pair_product_density_from_F = resolved

    F = parse_emitted_F(os.path.join(workspace, "bridge"))
    cfg = _load_workspace_config(workspace)
    n_pairs = cfg["Ncond"] // 2
    G_bridge = pair_product_density_from_F(F, n_pairs, rank_tol=1e-6)
    G_cuhf = _load_complexuhf_G(workspace, cfg["Ns"])
    return _check_g2_contraction(
        "g2a-emitted-F", workspace, G_bridge, G_cuhf, gtol
    )


def _dispatch_g2a_in_memory(workspace: str, resolved, gtol: float) -> int:
    """G2a-in-memory-A: conj(A_ship) @ A_ship.T vs ComplexUHF
    one-body Green at ``gtol`` (1e-6). A_ship is built via the canonical
    build_slater_orbitals; no shadow copy."""
    import numpy as np
    (build_slater_orbitals,) = resolved

    cfg = _load_workspace_config(workspace)
    A_ship, _ = _build_shipping_A(workspace, cfg)
    G_ship = np.conj(A_ship) @ A_ship.T
    G_cuhf = _load_complexuhf_G(workspace, cfg["Ns"])
    return _check_g2_contraction(
        "g2a-in-memory", workspace, G_ship, G_cuhf, gtol
    )


def _dispatch_g2b(workspace: str, resolved, gtol: float) -> int:
    """G2b: gauge_lift-lifted green_sublattice vs ComplexUHF
    one-body Green at ``gtol`` (1e-6). Uses boundary_theta from the
    workspace's input.toml. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    (gauge_lift,) = resolved

    cfg = _load_workspace_config(workspace)
    G_lift = _build_G_from_gauge_lift(workspace, cfg, gauge_lift)
    G_cuhf = _load_complexuhf_G(workspace, cfg["Ns"])
    return _check_g2_contraction(
        "g2b", workspace, G_lift, G_cuhf, gtol
    )


def _dispatch_g3(workspace: str, resolved, tol: float) -> int:
    (energy_relative_delta,) = resolved
    # Reuse the canonical G3 resolver so the artifact source contract
    # is honored.
    from tools._uhfk_to_mvmc.energy_compare import (
        EnergyCompareError,
        _resolve_g3_paths,
    )
    hwave_energy, mvmc_selected = _resolve_g3_paths(workspace)
    try:
        e_hwave, e_mvmc, delta_rel = energy_relative_delta(
            hwave_energy, mvmc_selected,
        )
    except EnergyCompareError as exc:
        print(f"G3: {exc}", file=sys.stderr)
        return 1
    # Same fail-open the G2 gates had: `nan > tol` is False, so a NaN
    # anywhere in the mVMC samples would fall straight through to a PASS
    # and silently disable the only H-wave/mVMC energy-agreement gate.
    import math

    for label, value in (
        ("tolerance", tol),
        ("H-wave energy", e_hwave),
        ("mVMC energy", e_mvmc),
        ("relative delta", delta_rel),
    ):
        if not math.isfinite(value):
            print(
                f"G3: {label} is non-finite ({value!r})",
                file=sys.stderr,
            )
            return 1
    if delta_rel > tol:
        return 1
    _emit_pass("g3", delta_rel, tol)
    return 0


_DISPATCHERS = {
    "g0-writer-check": _dispatch_g0_writer_check,
    "g1": _dispatch_g1,
    "g2a-emitted-F": _dispatch_g2a_emitted_F,
    "g2a-in-memory": _dispatch_g2a_in_memory,
    "g2b": _dispatch_g2b,
    "g3": _dispatch_g3,
}


# ---------------------------------------------------------------------
# CLI plumbing.
# ---------------------------------------------------------------------


def _run_workspace_mode(args) -> int:
    # Reject snapshot workspaces BEFORE any dispatcher can
    # print a PASS record. Import lazily so the legacy energy-only path
    # does not pay the price and so this module remains importable when
    # the scripts/ helper is not on sys.path.
    _scripts_dir = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "scripts")
    )
    if _scripts_dir not in sys.path:
        sys.path.insert(0, _scripts_dir)
    from _snapshot_guard import (
        reject_snapshot_workspace, SnapshotWorkspaceRejected,
    )
    try:
        reject_snapshot_workspace(args.workspace, helper_name="compare.py")
    except SnapshotWorkspaceRejected as exc:
        print(str(exc), file=sys.stderr)
        return 2

    _, resolved = _validate_dispatch(args.mode)
    dispatcher = _DISPATCHERS[args.mode]
    tol = _pick_tol_for_mode(args)
    try:
        return dispatcher(args.workspace, resolved, tol)
    except FileNotFoundError as e:
        print(f"compare.py --mode {args.mode}: {e}", file=sys.stderr)
        return 1


def _pick_tol_for_mode(args) -> float:
    m = args.mode
    if m == "g0-writer-check":
        return args.gtol_writer
    if m == "g1":
        return args.gtol_g1
    if m in ("g2a-emitted-F", "g2a-in-memory"):
        return args.gtol_ship
    if m == "g2b":
        return args.gtol_gauge
    if m == "g3":
        return args.etol
    return 0.0


def _legacy_energy_main(argv: List[str]) -> int:
    if len(argv) != 4:
        print(__doc__, file=sys.stderr)
        return 2
    hwave_path, mvmc_path, case = argv[1], argv[2], argv[3]
    e_uhf = parse_hwave_energy(hwave_path)
    bins = parse_mvmc_zvo_out(mvmc_path)
    e_mvmc = mean(bins)
    se_mvmc = stdev(bins) / (len(bins) ** 0.5) if len(bins) > 1 else 0.0
    delta = e_mvmc - e_uhf
    rel = abs(delta) / max(abs(e_uhf), 1e-12)
    if len(bins) >= 2 and se_mvmc > 0:
        sigma_sigmas = abs(delta) / se_mvmc
        print(f"[{case}] H-wave UHFk total = {e_uhf:.10e}")
        print(
            f"[{case}] mVMC initial mean = {e_mvmc:.10e}  "
            f"stderr = {se_mvmc:.2e}  (n_bins={len(bins)})"
        )
        print(
            f"[{case}] delta = {delta:+.4e}  ({sigma_sigmas:.2f} sigma, "
            f"{rel * 100:.2f}% rel)"
        )
        if sigma_sigmas <= 3.0 and rel < 0.05:
            print(f"[{case}] OK  (within 3 sigma and 5% rel)")
            return 0
        print(f"[{case}] FAIL  (delta exceeds 3 sigma or 5% rel)")
        return 1
    print(f"[{case}] H-wave UHFk total = {e_uhf:.10e}")
    print(
        f"[{case}] mVMC initial mean = {e_mvmc:.10e}  "
        f"(n_bins={len(bins)}; stderr undefined)"
    )
    print(f"[{case}] delta = {delta:+.4e}  ({rel * 100:.2f}% rel)")
    if rel < 0.01:
        print(f"[{case}] OK  (within 1% rel; single-bin VMC sample)")
        return 0
    print(f"[{case}] FAIL  (delta exceeds 1% rel)")
    return 1


def main(argv=None) -> int:
    argv = sys.argv if argv is None else argv
    # Legacy: exactly three positional args and none of the workspace flags.
    legacy_flag = ("--workspace", "--mode")
    if len(argv) == 4 and not any(a.startswith(legacy_flag) for a in argv[1:]):
        return _legacy_energy_main(argv)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", required=True)
    parser.add_argument(
        "--mode", required=True,
        choices=sorted(EXPECTED_MODE_DISPATCH.keys()),
    )
    parser.add_argument("--gtol_writer", type=float, default=1e-10)
    parser.add_argument("--gtol_g1", type=float, default=1e-10)
    parser.add_argument("--gtol_ship", type=float, default=1e-6)
    parser.add_argument("--gtol_gauge", type=float, default=1e-6)
    parser.add_argument("--etol", type=float, default=0.01)
    args = parser.parse_args(argv[1:])
    return _run_workspace_mode(args)


if __name__ == "__main__":
    sys.exit(main())
