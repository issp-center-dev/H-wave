"""Static-gate orchestrator and live UHF smoke launcher.

``trans.def`` substitution and marker writes use a temporary file, fsync,
and ``os.replace``, with temporary-file cleanup on failure. The mandatory
canonical-JSON marker binds the ``trans.def`` and seed-provenance hashes,
workspace SCF fingerprints, expected phase mask, ``namelist.def``, and every
resolved namelist target. The live smoke re-resolves and re-hashes the exact
bundle, rejects a missing or corrupt marker, and launches UHF only after the
static checks succeed. A non-zero UHF exit status is an error.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import subprocess
from pathlib import Path


_MARKER_FIELDS = (
    "expected_phase_mask",
    "namelist_sha256",
    "trans_def_sha256",
    "namelist_targets",
    "initial_def_provenance_sha256",
    "hwave_green_sha256",
    "hwave_eigen_sha256",
    "hwave_energy_total",
)
_NAMELIST_KEYS = (
    "ModPara",
    "LocSpin",
    "Trans",
    "CoulombIntra",
    "Orbital",
    "OrbitalParallel",
    "Initial",
)
_PROVENANCE_FIELDS = {
    "energy_total",
    "perturb_scale",
    "sha256_hwave_eigen",
    "sha256_hwave_green",
    "sha256_initial_def",
}
_ENERGY_ABS_TOL = 1e-12
# Lower bound, not an exact match. The seed perturbation exists so the
# G2 cross-solver density gates start OUTSIDE their own tolerance and
# have to contract into it; a perturbation at or below that tolerance
# lets a no-op solver pass by handing the seed straight back (the
# defect that shipped when this was pinned at 1e-6 == the G2 tolerance).
# The seeded density lands at roughly 0.5 * perturb_scale from the
# H-wave reference, and compare.py requires the initial delta to exceed
# 10 * tol with tol = 1e-6, so the seeding is only sound for
# perturb_scale >= 2e-5. The floor below keeps a safety factor on top of
# that while leaving the scale free to be tuned upward.
_MIN_PERTURB_SCALE = 1e-4


class Phase5StaticFailure(RuntimeError):
    """Static tier failure. Marker NOT written; live smoke MUST refuse."""


class Phase5LiveSmokeGuardError(RuntimeError):
    """Live smoke launcher refused because marker gate failed."""


def _sha256_bytes(data: bytes) -> str:
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest()


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fp:
        for chunk in iter(lambda: fp.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _atomic_write_bytes(dst: Path, data: bytes) -> None:
    """Atomic write via tmp-file + fsync + os.replace with cleanup on
    failure. On failure, dst is unchanged and tmp is cleaned up.

    Shared by ``atomic_substitute_trans_def`` and the
    ``phase5_static_ok.marker`` writer: the marker write must be
    atomic for the same reason trans.def substitution must be atomic
    -- a crash or concurrent write must never leave a truncated /
    corrupt file for a downstream reader to trip over.
    """
    dst = Path(dst)
    tmp = dst.with_name(dst.name + ".tmp")
    try:
        with open(tmp, "wb") as fp:
            fp.write(data)
            fp.flush()
            os.fsync(fp.fileno())
        os.replace(tmp, dst)
    except Exception:
        tmp.unlink(missing_ok=True)
        raise


def atomic_substitute_trans_def(src: bytes, dst: Path) -> None:
    """Atomically substitute dst with src bytes via tmp-file +
    fsync + os.replace. On failure, dst is unchanged and tmp is
    cleaned up.
    """
    _atomic_write_bytes(Path(dst), src)


def _raise(failure_type: type[RuntimeError], message: str) -> None:
    raise failure_type(message)


def _require_nonempty_file(
    path: Path, failure_type: type[RuntimeError],
) -> None:
    if not path.is_file():
        _raise(failure_type, f"{path} missing")
    if path.stat().st_size == 0:
        _raise(failure_type, f"{path} is empty")


def _resolve_workspace_file(
    declared_path: Path,
    workspace: Path,
    failure_type: type[RuntimeError],
) -> Path:
    """Resolve an existing non-empty file without permitting workspace escape."""
    workspace_root = Path(workspace).resolve(strict=True)
    try:
        resolved = Path(declared_path).resolve(strict=True)
    except (FileNotFoundError, OSError) as exc:
        _raise(failure_type, f"{declared_path} missing: {exc}")
    try:
        resolved.relative_to(workspace_root)
    except ValueError:
        _raise(
            failure_type,
            f"{declared_path} resolves outside workspace {workspace_root}: "
            f"{resolved}",
        )
    _require_nonempty_file(resolved, failure_type)
    return resolved


def _parse_namelist_bundle(
    workspace: Path,
    failure_type: type[RuntimeError],
) -> tuple[Path, dict[str, Path]]:
    """Strictly parse the exact ComplexUHF input bundle used at launch."""
    workspace = Path(workspace)
    launch_dir = (workspace / "complexuhf").resolve(strict=True)
    namelist_declared = launch_dir / "namelist.def"
    namelist_path = _resolve_workspace_file(
        namelist_declared, workspace, failure_type,
    )
    values = {key: [] for key in _NAMELIST_KEYS}
    for line_number, line in enumerate(namelist_path.read_text().splitlines(), 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        tokens = stripped.split()
        if len(tokens) != 2:
            _raise(
                failure_type,
                f"{namelist_path}:{line_number}: expected exactly "
                "'<key> <path>'",
            )
        key, value = tokens
        if key not in values:
            _raise(
                failure_type,
                f"{namelist_path}:{line_number}: unexpected key {key!r}",
            )
        values[key].append(value)

    for key in _NAMELIST_KEYS:
        if len(values[key]) != 1:
            _raise(
                failure_type,
                f"{namelist_path}: {key} must appear exactly once; "
                f"found {len(values[key])}",
            )

    resolved_targets = {}
    for key in _NAMELIST_KEYS:
        declared = Path(values[key][0])
        if not declared.is_absolute():
            declared = launch_dir / declared
        resolved_targets[key] = _resolve_workspace_file(
            declared, workspace, failure_type,
        )
    return namelist_path, resolved_targets


def _parse_expected_phase_mask(
    value, failure_type: type[RuntimeError], field: str = "expected_phase_mask",
) -> tuple[int, int, int]:
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        _raise(failure_type, f"{field} must contain exactly three entries")
    if any(type(entry) is not int or entry not in (0, 1) for entry in value):
        _raise(failure_type, f"{field} entries must be integers in {{0, 1}}")
    return tuple(value)


def _phase_values_to_mask(
    values, path: Path, failure_type: type[RuntimeError],
) -> tuple[int, int, int]:
    mask = []
    for axis, value in enumerate(values):
        if math.isclose(value, 0.0, rel_tol=0.0, abs_tol=1e-12):
            mask.append(0)
        elif math.isclose(abs(value), 180.0, rel_tol=0.0, abs_tol=1e-12):
            mask.append(1)
        else:
            _raise(
                failure_type,
                f"{path}: phase{axis} must be 0 or 180 degrees, got {value!r}",
            )
    return tuple(mask)


def _parse_stan_in(
    path: Path, failure_type: type[RuntimeError],
) -> tuple[int, tuple[int, int, int]]:
    _require_nonempty_file(path, failure_type)
    assignments = {}
    pattern = re.compile(r"^\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*([^\s/]+)")
    for line in path.read_text().splitlines():
        match = pattern.match(line)
        if match:
            assignments[match.group(1)] = match.group(2).strip('"')
    required = ("W", "L", "Height", "phase0", "phase1", "phase2")
    missing = [key for key in required if key not in assignments]
    if missing:
        _raise(failure_type, f"{path}: missing assignment(s) {missing}")
    try:
        dimensions = tuple(int(assignments[key]) for key in required[:3])
        phases = tuple(float(assignments[key]) for key in required[3:])
    except ValueError as exc:
        _raise(failure_type, f"{path}: malformed dimension or phase value: {exc}")
    if any(dimension <= 0 for dimension in dimensions):
        _raise(failure_type, f"{path}: dimensions must be positive: {dimensions}")
    nsite = math.prod(dimensions)
    return nsite, _phase_values_to_mask(phases, path, failure_type)


def _parse_modpara_nsite(
    path: Path, failure_type: type[RuntimeError],
) -> int:
    _require_nonempty_file(path, failure_type)
    pattern = re.compile(r"^\s*Nsite\s+([0-9]+)\s*$")
    values = [
        int(match.group(1))
        for line in path.read_text().splitlines()
        if (match := pattern.match(line))
    ]
    if len(values) != 1:
        _raise(
            failure_type,
            f"{path}: expected exactly one anchored 'Nsite <integer>' line, "
            f"found {len(values)}",
        )
    return values[0]


def _parse_geometry(
    path: Path, nsite: int, failure_type: type[RuntimeError],
) -> tuple[int, int, int]:
    _require_nonempty_file(path, failure_type)
    lines = [line.split() for line in path.read_text().splitlines() if line.strip()]
    expected_rows = nsite + 7
    if len(lines) != expected_rows:
        _raise(
            failure_type,
            f"{path}: geometry.dat has {len(lines) - 7} site rows; "
            f"expected {nsite}",
        )
    try:
        for tokens in lines[:4]:
            if len(tokens) != 3:
                raise ValueError("basis/phase rows require three columns")
            tuple(float(token) for token in tokens)
        for tokens in lines[4:7]:
            if len(tokens) != 3:
                raise ValueError("cell rows require three columns")
            tuple(int(token) for token in tokens)
        for tokens in lines[7:]:
            if len(tokens) != 4:
                raise ValueError("site rows require four columns")
            tuple(int(token) for token in tokens)
    except ValueError as exc:
        _raise(failure_type, f"{path}: malformed geometry.dat row: {exc}")
    phases = tuple(float(token) for token in lines[3])
    return _phase_values_to_mask(phases, path, failure_type)


def _numeric_rows(
    path: Path, columns: int, failure_type: type[RuntimeError],
) -> list[list[str]]:
    _require_nonempty_file(path, failure_type)
    rows = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("="):
            continue
        tokens = stripped.split()
        try:
            int(tokens[0])
        except (ValueError, IndexError):
            continue
        if len(tokens) != columns:
            _raise(
                failure_type,
                f"{path}: numeric row has {len(tokens)} columns; expected {columns}",
            )
        rows.append(tokens)
    return rows


def _validate_indexed_two_column_file(
    path: Path, nsite: int, failure_type: type[RuntimeError],
) -> list[list[str]]:
    rows = _numeric_rows(path, 2, failure_type)
    if len(rows) != nsite:
        _raise(
            failure_type,
            f"{path}: got {len(rows)} data rows; expected {nsite}",
        )
    try:
        indices = [int(row[0]) for row in rows]
    except ValueError as exc:
        _raise(failure_type, f"{path}: malformed row index: {exc}")
    if sorted(indices) != list(range(nsite)):
        _raise(failure_type, f"{path}: row indices must be exactly 0..{nsite - 1}")
    return rows


def _parse_unique_header_int(
    path: Path, key: str, failure_type: type[RuntimeError],
) -> int:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s+([0-9]+)\s*$")
    values = [
        int(match.group(1))
        for line in path.read_text().splitlines()
        if (match := pattern.match(line))
    ]
    if len(values) != 1:
        _raise(failure_type, f"{path}: expected exactly one {key} declaration")
    return values[0]


def _validate_coulombintra(
    path: Path, nsite: int, failure_type: type[RuntimeError],
) -> None:
    rows = _validate_indexed_two_column_file(path, nsite, failure_type)
    declared = _parse_unique_header_int(path, "NCoulombIntra", failure_type)
    if declared != len(rows):
        _raise(
            failure_type,
            f"{path}: NCoulombIntra declares {declared}, got {len(rows)} rows",
        )


def _validate_orbitalidx(
    path: Path, nsite: int, parallel: bool, failure_type: type[RuntimeError],
) -> None:
    _require_nonempty_file(path, failure_type)
    declared = _parse_unique_header_int(path, "NOrbitalIdx", failure_type)
    mapping = []
    optimize = []
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("="):
            continue
        tokens = stripped.split()
        if tokens[0] in ("NOrbitalIdx", "ComplexType"):
            continue
        try:
            ints = tuple(int(token) for token in tokens)
        except ValueError as exc:
            _raise(failure_type, f"{path}: malformed integer row: {exc}")
        if len(ints) == 4:
            mapping.append(ints)
        elif len(ints) == 2:
            optimize.append(ints)
        else:
            _raise(failure_type, f"{path}: expected a two- or four-column row")

    expected_mapping = nsite * (nsite - 1) // 2 if parallel else nsite * nsite
    if len(mapping) != expected_mapping:
        _raise(
            failure_type,
            f"{path}: got {len(mapping)} mapping rows; expected {expected_mapping}",
        )
    if len(optimize) != declared:
        _raise(
            failure_type,
            f"{path}: got {len(optimize)} optimize rows; "
            f"NOrbitalIdx declares {declared}",
        )
    pairs = [(row[0], row[1]) for row in mapping]
    if len(set(pairs)) != len(pairs):
        _raise(failure_type, f"{path}: duplicate site-pair mapping row")
    if any(not (0 <= i < nsite and 0 <= j < nsite) for i, j in pairs):
        _raise(failure_type, f"{path}: mapping references site outside 0..{nsite - 1}")
    expected_pairs = (
        {(i, j) for i in range(nsite) for j in range(i + 1, nsite)}
        if parallel
        else {(i, j) for i in range(nsite) for j in range(nsite)}
    )
    if set(pairs) != expected_pairs:
        _raise(failure_type, f"{path}: site-pair mapping is incomplete")
    optimize_indices = [row[0] for row in optimize]
    if sorted(optimize_indices) != list(range(declared)):
        _raise(failure_type, f"{path}: optimize indices must be exactly 0..{declared - 1}")


def _parse_trans_def(
    path: Path, failure_type: type[RuntimeError],
) -> dict[tuple[int, int, int, int], complex]:
    _require_nonempty_file(path, failure_type)
    declared = _parse_unique_header_int(path, "NTransfer", failure_type)
    rows = {}
    actual_rows = 0
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("=") or stripped.startswith("NTransfer"):
            continue
        tokens = stripped.split()
        if len(tokens) != 6:
            _raise(failure_type, f"{path}: trans.def row must have six columns")
        try:
            key = tuple(int(token) for token in tokens[:4])
            value = complex(float(tokens[4]), float(tokens[5]))
        except ValueError as exc:
            _raise(failure_type, f"{path}: malformed trans.def row: {exc}")
        if not math.isfinite(value.real) or not math.isfinite(value.imag):
            _raise(failure_type, f"{path}: non-finite transfer value for {key}")
        if key in rows:
            _raise(failure_type, f"{path}: duplicate transfer row {key}")
        rows[key] = value
        actual_rows += 1
    if declared != actual_rows:
        _raise(
            failure_type,
            f"{path}: NTransfer declares {declared}, got {actual_rows} rows",
        )
    for (i, spin_i, j, spin_j), value in rows.items():
        partner_key = (j, spin_j, i, spin_i)
        partner = rows.get(partner_key)
        if partner is None or partner != value.conjugate():
            _raise(
                failure_type,
                f"{path}: Hermiticity failure for {(i, spin_i, j, spin_j)}; "
                f"partner {partner_key} is {partner!r}, expected "
                f"{value.conjugate()!r}",
            )
    return rows


def _validate_runtime_bundle(
    workspace: Path,
    expected_phase_mask: tuple[int, int, int],
    failure_type: type[RuntimeError],
    namelist_targets: dict[str, Path],
    fixture_root: Path | None = None,
) -> None:
    complexuhf = workspace / "complexuhf"
    workspace_nsite, workspace_mask = _parse_stan_in(
        complexuhf / "stan.in", failure_type,
    )
    if fixture_root is not None:
        fixture_nsite, fixture_mask = _parse_stan_in(
            Path(fixture_root) / "stan.in", failure_type,
        )
        if fixture_mask != expected_phase_mask:
            _raise(
                failure_type,
                f"{fixture_root / 'stan.in'}: phase mask {fixture_mask} does "
                f"not match expected phase mask {expected_phase_mask}",
            )
        if workspace_nsite != fixture_nsite:
            _raise(
                failure_type,
                f"{complexuhf / 'stan.in'}: Nsite {workspace_nsite} does not "
                f"match fixture-derived Nsite {fixture_nsite}",
            )
        nsite = fixture_nsite
    else:
        nsite = workspace_nsite
    if workspace_mask != expected_phase_mask:
        _raise(
            failure_type,
            f"{complexuhf / 'stan.in'}: phase mask {workspace_mask} does not "
            f"match expected phase mask {expected_phase_mask}",
        )
    modpara_path = namelist_targets["ModPara"]
    modpara_nsite = _parse_modpara_nsite(modpara_path, failure_type)
    if modpara_nsite != nsite:
        _raise(
            failure_type,
            f"{modpara_path}: Nsite is {modpara_nsite}; "
            f"fixture-derived expectation is {nsite}",
        )
    geometry_mask = _parse_geometry(
        complexuhf / "geometry.dat", nsite, failure_type,
    )
    if geometry_mask != expected_phase_mask:
        _raise(
            failure_type,
            f"{complexuhf / 'geometry.dat'}: phase mask {geometry_mask} does "
            f"not match expected phase mask {expected_phase_mask}",
        )
    _validate_indexed_two_column_file(
        namelist_targets["LocSpin"], nsite, failure_type,
    )
    _validate_coulombintra(
        namelist_targets["CoulombIntra"], nsite, failure_type,
    )
    _validate_orbitalidx(
        namelist_targets["Orbital"], nsite, False, failure_type,
    )
    _validate_orbitalidx(
        namelist_targets["OrbitalParallel"], nsite, True, failure_type,
    )
    trans_def = namelist_targets["Trans"]
    emitted_trans_def = workspace / "mvmc" / "trans.def"
    _require_nonempty_file(trans_def, failure_type)
    _require_nonempty_file(emitted_trans_def, failure_type)
    if trans_def.read_bytes() != emitted_trans_def.read_bytes():
        _raise(
            failure_type,
            f"{trans_def} is not byte-equal to bridge-emitted {emitted_trans_def}",
        )
    _parse_trans_def(trans_def, failure_type)


def _read_energy_total(
    path: Path, failure_type: type[RuntimeError],
) -> tuple[str, float]:
    _require_nonempty_file(path, failure_type)
    pattern = re.compile(r"^\s*Energy_Total\s*=\s*(\S+)\s*$")
    values = [
        match.group(1)
        for line in path.read_text().splitlines()
        if (match := pattern.match(line))
    ]
    if len(values) != 1:
        _raise(failure_type, f"{path}: expected exactly one Energy_Total line")
    try:
        numeric = float(values[0])
    except ValueError as exc:
        _raise(failure_type, f"{path}: malformed Energy_Total: {exc}")
    if not math.isfinite(numeric):
        _raise(failure_type, f"{path}: Energy_Total must be finite")
    return values[0], numeric


def _collect_scf_fingerprint(
    workspace: Path, failure_type: type[RuntimeError],
) -> tuple[str, str, str, float]:
    output = workspace / "hwave" / "output"
    green = output / "green.npz"
    eigen = output / "eigen.npz"
    energy_dat = output / "energy.dat"
    for path in (green, eigen, energy_dat):
        _require_nonempty_file(path, failure_type)
    energy_text, energy_value = _read_energy_total(energy_dat, failure_type)
    return (
        _sha256_file(green),
        _sha256_file(eigen),
        energy_text,
        energy_value,
    )


def _validate_provenance(
    initial_def: Path,
    green_sha: str,
    eigen_sha: str,
    energy_value: float,
    failure_type: type[RuntimeError],
) -> Path:
    provenance_path = initial_def.with_name(initial_def.name + ".provenance")
    _require_nonempty_file(initial_def, failure_type)
    _require_nonempty_file(provenance_path, failure_type)
    try:
        provenance = json.loads(provenance_path.read_text())
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        _raise(failure_type, f"{provenance_path}: malformed JSON: {exc}")
    if not isinstance(provenance, dict):
        _raise(failure_type, f"{provenance_path}: JSON root must be an object")
    missing = sorted(_PROVENANCE_FIELDS - set(provenance))
    unexpected = sorted(set(provenance) - _PROVENANCE_FIELDS)
    if missing:
        _raise(failure_type, f"{provenance_path}: missing field(s) {missing}")
    if unexpected:
        _raise(failure_type, f"{provenance_path}: unexpected field(s) {unexpected}")

    expected_hashes = {
        "sha256_initial_def": _sha256_file(initial_def),
        "sha256_hwave_green": green_sha,
        "sha256_hwave_eigen": eigen_sha,
    }
    sha_pattern = re.compile(r"^[0-9a-f]{64}$")
    for field, expected in expected_hashes.items():
        value = provenance[field]
        if not isinstance(value, str) or not sha_pattern.fullmatch(value):
            _raise(failure_type, f"{provenance_path}: malformed {field}")
        if value != expected:
            _raise(
                failure_type,
                f"{provenance_path}: {field} mismatch: records {value!r}, "
                f"current artifact is {expected!r}",
            )
    provenance_energy = provenance["energy_total"]
    if type(provenance_energy) not in (int, float) or not math.isfinite(
        provenance_energy,
    ):
        _raise(failure_type, f"{provenance_path}: malformed energy_total")
    if not math.isclose(
        float(provenance_energy), energy_value,
        rel_tol=0.0, abs_tol=_ENERGY_ABS_TOL,
    ):
        _raise(
            failure_type,
            f"{provenance_path}: energy_total mismatch: records "
            f"{provenance_energy!r}, energy.dat is {energy_value!r}",
        )
    perturb_scale = provenance["perturb_scale"]
    if type(perturb_scale) not in (int, float) or not math.isfinite(perturb_scale):
        _raise(failure_type, f"{provenance_path}: malformed perturb_scale")
    if float(perturb_scale) < _MIN_PERTURB_SCALE:
        _raise(
            failure_type,
            f"{provenance_path}: perturb_scale is {perturb_scale!r}, "
            f"below the {_MIN_PERTURB_SCALE!r} floor required to keep "
            "the G2 density gates non-vacuous; a seed at or inside the "
            "G2 tolerance lets a no-op solver pass by returning it "
            "unchanged",
        )
    return provenance_path


def _current_marker_payload(
    workspace: Path,
    expected_phase_mask: tuple[int, int, int],
    failure_type: type[RuntimeError],
    namelist_path: Path,
    namelist_targets: dict[str, Path],
) -> dict:
    green_sha, eigen_sha, energy_text, energy_value = _collect_scf_fingerprint(
        workspace, failure_type,
    )
    provenance_path = _validate_provenance(
        namelist_targets["Initial"],
        green_sha,
        eigen_sha,
        energy_value,
        failure_type,
    )
    trans_def = namelist_targets["Trans"]
    _require_nonempty_file(trans_def, failure_type)
    workspace_root = Path(workspace).resolve(strict=True)
    target_records = {
        key: {
            "path": path.relative_to(workspace_root).as_posix(),
            "sha256": _sha256_file(path),
        }
        for key, path in namelist_targets.items()
    }
    return {
        "expected_phase_mask": list(expected_phase_mask),
        "namelist_sha256": _sha256_file(namelist_path),
        "namelist_targets": target_records,
        "trans_def_sha256": _sha256_file(trans_def),
        "initial_def_provenance_sha256": _sha256_file(provenance_path),
        "hwave_green_sha256": green_sha,
        "hwave_eigen_sha256": eigen_sha,
        "hwave_energy_total": energy_text,
    }


def run_static_bundle_validator(
    fixture_root: Path,
    workspace: Path,
    expected_phase_mask: tuple,
    _skip_bundle_files_for_test: bool = False,
) -> Path:
    """Run all static checks and write a canonical-JSON marker.

    Returns Path to the marker; raises Phase5StaticFailure on any
    check failure (no marker written on failure).

    ``_skip_bundle_files_for_test`` skips only fixture/bundle structure
    validation. SCF artifacts, trans.def, initial.def, and parsed seed
    provenance remain mandatory and verified. The sanctioned real-run
    entry point, :func:`run_phase5_gate`, does not expose this test-only
    switch, and the live launcher always validates the complete bundle.
    """
    workspace = Path(workspace)
    expected_phase_mask = _parse_expected_phase_mask(
        expected_phase_mask, Phase5StaticFailure,
    )
    namelist_path, namelist_targets = _parse_namelist_bundle(
        workspace, Phase5StaticFailure,
    )
    if not _skip_bundle_files_for_test:
        _validate_runtime_bundle(
            workspace,
            expected_phase_mask,
            Phase5StaticFailure,
            namelist_targets,
            fixture_root=Path(fixture_root),
        )

    marker_payload = _current_marker_payload(
        workspace,
        expected_phase_mask,
        Phase5StaticFailure,
        namelist_path,
        namelist_targets,
    )
    for field in _MARKER_FIELDS:
        if marker_payload[field] in (None, "", []):
            raise Phase5StaticFailure(
                f"refusing to write marker with empty {field}"
            )
    for key, record in marker_payload["namelist_targets"].items():
        if record["path"] in (None, "") or record["sha256"] in (None, ""):
            raise Phase5StaticFailure(
                f"refusing to write marker with empty namelist target {key}"
            )
    marker_path = workspace / "phase5_static_ok.marker"
    marker_bytes = json.dumps(
        marker_payload, sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    _atomic_write_bytes(marker_path, marker_bytes)
    return marker_path


def run_live_uhf_smoke(
    workspace: Path, uhf_binary: str,
) -> subprocess.CompletedProcess:
    """Live UHF smoke launcher. Requires marker; recomputes every
    marker-referenced hash from the CURRENT workspace and refuses
    to invoke UHF if anything mismatches.

    There is no `require_marker` opt-out.
    """
    workspace = Path(workspace)
    marker_path = workspace / "phase5_static_ok.marker"
    if not marker_path.is_file():
        raise Phase5LiveSmokeGuardError(
            f"{marker_path} missing; static gate has not run or failed."
        )
    try:
        marker_raw = marker_path.read_text()
        marker = json.loads(marker_raw)
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise Phase5LiveSmokeGuardError(
            f"{marker_path} is corrupt / not valid JSON: {exc}"
        ) from exc
    if not isinstance(marker, dict):
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: marker JSON root must be an object"
        )
    missing = sorted(set(_MARKER_FIELDS) - set(marker))
    unexpected = sorted(set(marker) - set(_MARKER_FIELDS))
    if missing:
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: missing marker field(s) {missing}"
        )
    if unexpected:
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: unexpected marker field(s) {unexpected}"
        )
    marker_phase_mask = _parse_expected_phase_mask(
        marker["expected_phase_mask"],
        Phase5LiveSmokeGuardError,
    )
    for field in (
        "namelist_sha256",
        "trans_def_sha256",
        "initial_def_provenance_sha256",
        "hwave_green_sha256",
        "hwave_eigen_sha256",
        "hwave_energy_total",
    ):
        if not isinstance(marker[field], str) or not marker[field]:
            raise Phase5LiveSmokeGuardError(
                f"{marker_path}: {field} must be a non-empty string"
            )
    marker_targets = marker["namelist_targets"]
    if not isinstance(marker_targets, dict):
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: namelist_targets must be an object"
        )
    missing_targets = sorted(set(_NAMELIST_KEYS) - set(marker_targets))
    unexpected_targets = sorted(set(marker_targets) - set(_NAMELIST_KEYS))
    if missing_targets or unexpected_targets:
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: namelist_targets keys mismatch; missing "
            f"{missing_targets}, unexpected {unexpected_targets}"
        )
    sha_pattern = re.compile(r"^[0-9a-f]{64}$")
    for key in _NAMELIST_KEYS:
        record = marker_targets[key]
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            raise Phase5LiveSmokeGuardError(
                f"{marker_path}: namelist_targets.{key} must contain exactly "
                "path and sha256"
            )
        if not isinstance(record["path"], str) or not record["path"]:
            raise Phase5LiveSmokeGuardError(
                f"{marker_path}: namelist_targets.{key}.path must be non-empty"
            )
        if not isinstance(record["sha256"], str) or not sha_pattern.fullmatch(
            record["sha256"],
        ):
            raise Phase5LiveSmokeGuardError(
                f"{marker_path}: namelist_targets.{key}.sha256 is malformed"
            )
    canonical_marker = json.dumps(
        marker, sort_keys=True, separators=(",", ":"),
    )
    if marker_raw != canonical_marker:
        raise Phase5LiveSmokeGuardError(
            f"{marker_path}: marker is not canonical JSON"
        )

    _, workspace_phase_mask = _parse_stan_in(
        workspace / "complexuhf" / "stan.in",
        Phase5LiveSmokeGuardError,
    )
    namelist_path, namelist_targets = _parse_namelist_bundle(
        workspace, Phase5LiveSmokeGuardError,
    )
    current = _current_marker_payload(
        workspace,
        workspace_phase_mask,
        Phase5LiveSmokeGuardError,
        namelist_path,
        namelist_targets,
    )
    for field in _MARKER_FIELDS:
        if marker[field] != current[field]:
            raise Phase5LiveSmokeGuardError(
                f"{field} mismatch: marker records {marker[field]!r}, "
                f"current workspace is {current[field]!r}"
            )

    # The test-only static skip cannot weaken a real run: live smoke
    # independently re-runs every runtime bundle check.
    _validate_runtime_bundle(
        workspace,
        workspace_phase_mask,
        Phase5LiveSmokeGuardError,
        namelist_targets,
    )

    # Re-resolve and re-hash the launch manifest and every selected target
    # immediately before exec. This catches a swap after bundle validation.
    final_namelist_path, final_namelist_targets = _parse_namelist_bundle(
        workspace, Phase5LiveSmokeGuardError,
    )
    final = _current_marker_payload(
        workspace,
        workspace_phase_mask,
        Phase5LiveSmokeGuardError,
        final_namelist_path,
        final_namelist_targets,
    )
    for field in _MARKER_FIELDS:
        if marker[field] != final[field]:
            raise Phase5LiveSmokeGuardError(
                f"pre-exec {field} mismatch: marker records "
                f"{marker[field]!r}, current workspace is {final[field]!r}"
            )

    res = subprocess.run(
        [uhf_binary, "namelist.def"],
        cwd=str(workspace / "complexuhf"),
        capture_output=True, text=True,
    )
    if res.returncode != 0:
        raise Phase5LiveSmokeGuardError(
            f"UHF exited with returncode={res.returncode}. "
            f"stderr={res.stderr!r} stdout={res.stdout!r}"
        )
    return res


def run_phase5_gate(
    fixture_root: Path,
    workspace: Path,
    expected_phase_mask: tuple,
    uhf_binary: str,
) -> "subprocess.CompletedProcess":
    """Sole entry point for the static and live gates.

    Enforces static-first ordering BY CONSTRUCTION: this function is
    the only sanctioned way to run both tiers, and it always runs
    the static bundle validator first. If the static tier fails,
    Phase5StaticFailure propagates and the live UHF smoke is
    NEVER invoked. Callers MUST use this entry point rather than
    calling the two sub-functions independently — that keeps the
    ordering guarantee documented in the module docstring safely
    inside the process.

    Returns the CompletedProcess from the successful live smoke;
    raises Phase5StaticFailure or Phase5LiveSmokeGuardError on any
    failure.
    """
    run_static_bundle_validator(
        fixture_root=fixture_root,
        workspace=workspace,
        expected_phase_mask=expected_phase_mask,
    )
    return run_live_uhf_smoke(
        workspace=workspace,
        uhf_binary=uhf_binary,
    )
