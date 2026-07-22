"""phase5_gate.py orchestrator contract tests.

Includes tests for residual risks: marker gate cannot be
bypassed by direct live-smoke call, static failure prevents UHF exec,
stale marker rejected, atomic trans.def substitution on failure.
"""
from __future__ import annotations

import json
import hashlib
import importlib.util
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc import phase5_gate as phase5_gate_module
from tools._uhfk_to_mvmc.phase5_gate import (
    Phase5StaticFailure,
    Phase5LiveSmokeGuardError,
    run_static_bundle_validator,
    run_live_uhf_smoke,
    atomic_substitute_trans_def,
)
from tools._uhfk_to_mvmc.trans_emit import emit_trans_def


_REPO_ROOT = Path(__file__).resolve().parents[1]
_REAL_CASE_NAME = "case_soc_rashba_3d_sub_apbc_xy"
_REAL_CASE = (
    _REPO_ROOT / "tests/validation/uhfk_mvmc_pairproduct" / _REAL_CASE_NAME
)
_REAL_FIXTURE = (
    _REPO_ROOT / "tests/validation/apbc_complexuhf" / _REAL_CASE_NAME
)
_REAL_PHASE_MASK = (1, 1, 0)
_REAL_SNAPSHOT_PREFIX = "v37_case_soc_rashba_3d_sub_apbc_xy"
_SEED_SCRIPT = (
    _REPO_ROOT
    / "tests/validation/uhfk_mvmc_pairproduct/scripts/seed_complexuhf_from_hwave.py"
)
_REAL_WORKSPACE_TEMPLATE: Path | None = None
# Must match what run.sh passes to seed_complexuhf_from_hwave.py and stay
# at or above phase5_gate._MIN_PERTURB_SCALE. Seeding below that floor is
# what made the G2 density gates passable by a no-op solver.
_SHIPPED_PERTURB_SCALE = 1e-3
_NAMELIST_TARGETS = {
    "ModPara": "modpara.def",
    "LocSpin": "locspn.def",
    "Trans": "trans.def",
    "CoulombIntra": "coulombintra.def",
    "Orbital": "orbitalidx.def",
    "OrbitalParallel": "orbitalidxpara.def",
    "Initial": "initial.def",
}


def _write_namelist(complexuhf: Path, targets=None) -> None:
    selected = _NAMELIST_TARGETS if targets is None else targets
    (complexuhf / "namelist.def").write_text(
        "".join(f"{key:>16}  {value}\n" for key, value in selected.items())
    )


@pytest.fixture(scope="module", autouse=True)
def _real_workspace_template(tmp_path_factory):
    global _REAL_WORKSPACE_TEMPLATE
    template = tmp_path_factory.mktemp("phase5-real-workspace") / "work"
    hwave = template / "hwave"
    complexuhf = template / "complexuhf"
    mvmc = template / "mvmc"
    output = hwave / "output"
    output.mkdir(parents=True)
    # Build from tracked snapshots so this module runs on a clean checkout.
    # trans.def is deliberately generated below through the production bridge
    # emitter instead of adding a second golden copy that could drift.
    for name in ("green", "eigen", "occupation"):
        shutil.copy2(
            _REPO_ROOT / "tests/data" / f"{_REAL_SNAPSHOT_PREFIX}_{name}.npz",
            output / f"{name}.npz",
        )
    (output / "energy.dat").write_text(
        "Energy_Total = -25.371716580385547\n"
    )
    for name in ("input.toml", "geometry_uhf.dat"):
        shutil.copy2(_REAL_CASE / name, hwave / name)
    shutil.copytree(_REAL_FIXTURE, complexuhf)
    mvmc.mkdir()
    emit_trans_def(
        _REAL_CASE / "Transfer.dat",
        [4, 4, 4],
        mvmc / "trans.def",
        boundary_theta=[math.pi, math.pi, 0.0],
    )
    shutil.copy2(mvmc / "trans.def", complexuhf / "trans.def")

    spec = importlib.util.spec_from_file_location("phase5_seed", _SEED_SCRIPT)
    seed = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(seed)
    cfg = seed._find_workspace_config(hwave)
    seed._write_initial_def(
        seed._build_shipping_G(hwave, cfg),
        complexuhf / "initial.def",
        perturb_scale=_SHIPPED_PERTURB_SCALE,
        perturb_seed=42,
    )
    seed._write_provenance(
        complexuhf / "initial.def", hwave, perturb_scale=_SHIPPED_PERTURB_SCALE,
    )
    _REAL_WORKSPACE_TEMPLATE = template
    yield
    _REAL_WORKSPACE_TEMPLATE = None


def _copy_real_workspace(tmp_path: Path) -> tuple[Path, Path]:
    """Copy a shipped workspace at its pre-static-check layout."""
    assert _REAL_WORKSPACE_TEMPLATE is not None
    workspace = tmp_path / "work"
    fixture_root = tmp_path / "fixture"
    shutil.copytree(_REAL_WORKSPACE_TEMPLATE, workspace)
    shutil.copytree(_REAL_FIXTURE, fixture_root)
    return fixture_root, workspace


def _write_valid_provenance(workspace: Path) -> None:
    hwave_output = workspace / "hwave" / "output"
    complexuhf = workspace / "complexuhf"
    energy = -25.371716580385547
    payload = {
        "energy_total": energy,
        "perturb_scale": _SHIPPED_PERTURB_SCALE,
        "sha256_hwave_eigen": hashlib.sha256(
            (hwave_output / "eigen.npz").read_bytes(),
        ).hexdigest(),
        "sha256_hwave_green": hashlib.sha256(
            (hwave_output / "green.npz").read_bytes(),
        ).hexdigest(),
        "sha256_initial_def": hashlib.sha256(
            (complexuhf / "initial.def").read_bytes(),
        ).hexdigest(),
    }
    (complexuhf / "initial.def.provenance").write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")),
    )


def _make_workspace_and_bundle(tmp_path, phase_mask=(0, 0, 0)):
    """Build a minimal workspace + bundle that passes the static tier
    for these unit tests. Live UHF binary itself is mocked out."""
    (tmp_path / "hwave" / "output").mkdir(parents=True)
    (tmp_path / "complexuhf").mkdir()
    (tmp_path / "mvmc").mkdir()
    # Placeholder hwave outputs (real content is not exercised in unit
    # tests; the fingerprint just needs to be reproducible).
    for fname in ("green.npz", "eigen.npz"):
        (tmp_path / "hwave" / "output" / fname).write_bytes(b"placeholder")
    (tmp_path / "hwave" / "output" / "energy.dat").write_text(
        "Energy_Total = -25.371716580385547\n"
    )
    # Minimal internally consistent ComplexUHF bundle.
    (tmp_path / "complexuhf" / "modpara.def").write_text(
        "-----\nNsite   1\n-----\n"
    )
    (tmp_path / "complexuhf" / "stan.in").write_text(
        "W = 1\nL = 1\nHeight = 1\n"
        "phase0 = 0.0\nphase1 = 0.0\nphase2 = 0.0\n"
    )
    (tmp_path / "complexuhf" / "geometry.dat").write_text(
        "1 0 0\n0 1 0\n0 0 1\n0 0 0\n"
        "1 0 0\n0 1 0\n0 0 1\n0 0 0 0\n"
    )
    (tmp_path / "complexuhf" / "locspn.def").write_text(
        "NlocalSpin 0\n0 0\n"
    )
    (tmp_path / "complexuhf" / "coulombintra.def").write_text(
        "NCoulombIntra 1\n0 2.0\n"
    )
    (tmp_path / "complexuhf" / "orbitalidx.def").write_text(
        "NOrbitalIdx 1\nComplexType 1\n0 0 0 1\n0 1\n"
    )
    (tmp_path / "complexuhf" / "orbitalidxpara.def").write_text(
        "NOrbitalIdx 0\nComplexType 1\n"
    )
    trans = "NTransfer 1\n0 0 0 0 1.0 0.0\n"
    (tmp_path / "complexuhf" / "trans.def").write_text(trans)
    (tmp_path / "mvmc" / "trans.def").write_text(trans)
    (tmp_path / "complexuhf" / "initial.def").write_bytes(b"initial-def\n")
    _write_namelist(tmp_path / "complexuhf")
    _write_valid_provenance(tmp_path)
    return tmp_path


def test_static_success_produces_marker(tmp_path):
    """Static tier success writes phase5_static_ok.marker with a
    canonical JSON digest binding all post-substitution artifacts."""
    ws = _make_workspace_and_bundle(tmp_path)
    marker = run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    assert marker.exists()
    payload = json.loads(marker.read_text())
    assert "trans_def_sha256" in payload
    assert "initial_def_provenance_sha256" in payload
    assert "hwave_green_sha256" in payload
    assert payload["namelist_sha256"]
    assert set(payload["namelist_targets"]) == set(_NAMELIST_TARGETS)
    for key, target in payload["namelist_targets"].items():
        assert target["path"] == f"complexuhf/{_NAMELIST_TARGETS[key]}"
        assert target["sha256"]


def test_static_rejects_namelist_initial_outside_workspace(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    outside = tmp_path.parent / f"{tmp_path.name}-outside-initial.def"
    outside.write_bytes(b"outside seed\n")
    targets = dict(_NAMELIST_TARGETS, Initial=str(outside))
    _write_namelist(ws / "complexuhf", targets)

    with pytest.raises(Phase5StaticFailure, match="outside workspace"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_rejects_duplicate_namelist_key(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    namelist = ws / "complexuhf" / "namelist.def"
    namelist.write_text(namelist.read_text() + "Initial initial.def\n")

    with pytest.raises(Phase5StaticFailure, match="Initial.*exactly once"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_rejects_missing_namelist_key(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    targets = dict(_NAMELIST_TARGETS)
    del targets["Initial"]
    _write_namelist(ws / "complexuhf", targets)

    with pytest.raises(Phase5StaticFailure, match="Initial.*exactly once"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_rejects_namelist_symlink_escaping_workspace(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    outside = tmp_path.parent / f"{tmp_path.name}-symlink-target.def"
    outside.write_bytes(b"outside seed\n")
    (ws / "complexuhf" / "escaping-initial.def").symlink_to(outside)
    targets = dict(_NAMELIST_TARGETS, Initial="escaping-initial.def")
    _write_namelist(ws / "complexuhf", targets)

    with pytest.raises(Phase5StaticFailure, match="outside workspace"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_rejects_missing_namelist_target(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    targets = dict(_NAMELIST_TARGETS, Initial="missing-initial.def")
    _write_namelist(ws / "complexuhf", targets)

    with pytest.raises(Phase5StaticFailure, match="missing-initial.def.*missing"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_validates_namelist_selected_initial(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    alternate = ws / "complexuhf" / "alternate-initial.def"
    alternate.write_bytes(b"different seed\n")
    targets = dict(_NAMELIST_TARGETS, Initial=alternate.name)
    _write_namelist(ws / "complexuhf", targets)

    with pytest.raises(Phase5StaticFailure, match="provenance.*missing"):
        run_static_bundle_validator(
            fixture_root=tmp_path / "fake_fixture",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=True,
        )


def test_static_failure_raises_without_marker(tmp_path):
    """If any static check fails, Phase5StaticFailure raises and NO
    marker is written."""
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "modpara.def").write_text("Nsite   32\n")
    with pytest.raises(Phase5StaticFailure, match="Nsite"):
        run_static_bundle_validator(
            fixture_root=ws / "complexuhf",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=False,
        )
    assert not (ws / "phase5_static_ok.marker").exists()


def test_live_smoke_refuses_without_static_marker(tmp_path):
    """live_smoke MUST refuse to spawn UHF if
    the marker is missing. No require_marker opt-out parameter exists."""
    ws = _make_workspace_and_bundle(tmp_path)
    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError, match="marker"):
            run_live_uhf_smoke(
                workspace=ws, uhf_binary="/nonexistent/UHF",
            )
        assert not run_mock.called, (
            "subprocess.run was invoked despite missing marker; "
            "gate bypass."
        )


def test_live_smoke_refuses_stale_marker(tmp_path):
    """Mutating an artifact hash referenced by
    the marker AFTER static success MUST cause live_smoke to refuse
    UHF exec."""
    ws = _make_workspace_and_bundle(tmp_path)
    # Provide a trans.def whose SHA the marker will record.
    (ws / "complexuhf" / "trans.def").write_bytes(b"initial-trans-def")
    marker = run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    # Attacker mutates the post-substitution trans.def AFTER marker
    # was written.
    (ws / "complexuhf" / "trans.def").write_bytes(b"hostile-trans-def")
    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError, match="sha256 mismatch"):
            run_live_uhf_smoke(
                workspace=ws, uhf_binary="/nonexistent/UHF",
            )
        assert not run_mock.called


def test_atomic_trans_def_writeout_on_failure(tmp_path):
    """If os.replace fails after
    writing to trans.def.tmp, the target trans.def is unchanged and
    the temp file is cleaned up."""
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "trans.def").write_bytes(b"existing")
    # Simulate failure at the atomic step.
    with patch("os.replace", side_effect=OSError("simulated")):
        with pytest.raises(OSError):
            atomic_substitute_trans_def(
                src=b"replacement-bytes",
                dst=ws / "complexuhf" / "trans.def",
            )
    # Existing trans.def is preserved.
    assert (ws / "complexuhf" / "trans.def").read_bytes() == b"existing"
    # tmp file is cleaned up.
    assert not (ws / "complexuhf" / "trans.def.tmp").exists()


def test_run_phase5_gate_static_failure_aborts_before_live_smoke(tmp_path):
    """If the static tier fails, the
    live smoke launcher is NEVER reached — the safety property (static
    must precede live) holds by construction."""
    from tools._uhfk_to_mvmc.phase5_gate import run_phase5_gate
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "modpara.def").write_text("Nsite   32\n")
    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5StaticFailure):
            run_phase5_gate(
                fixture_root=tmp_path / "fake_fixture",
                workspace=ws,
                expected_phase_mask=(0, 0, 0),
                uhf_binary="/nonexistent/UHF",
            )
        assert not run_mock.called, (
            "subprocess.run was invoked despite static failure; the "
            "wrapper does not enforce static-first ordering."
        )


def test_live_smoke_refuses_deleted_trans_def(tmp_path):
    """If trans.def is DELETED (not mutated) after the marker records
    its hash, the current live_smoke MUST refuse rather than silently
    fall through to subprocess.run."""
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "trans.def").write_bytes(b"initial-trans-def")
    run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    (ws / "complexuhf" / "trans.def").unlink()  # attacker deletes it
    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError, match="missing"):
            run_live_uhf_smoke(
                workspace=ws, uhf_binary="/nonexistent/UHF",
            )
        assert not run_mock.called


def test_static_validator_fails_on_missing_trans_def(tmp_path):
    """Critical #1 closure: workspace without trans.def at static-
    check time raises Phase5StaticFailure, not silent empty-SHA."""
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "trans.def").unlink()
    with pytest.raises(Phase5StaticFailure, match="trans.def missing"):
        run_static_bundle_validator(
            fixture_root=ws / "complexuhf",
            workspace=ws,
            expected_phase_mask=(0, 0, 0),
            _skip_bundle_files_for_test=False,  # runtime mode
        )


def test_live_smoke_refuses_mutated_green_npz(tmp_path):
    """Important #2 closure: marker records hwave_green_sha256; if
    green.npz is mutated after the marker was written, live smoke
    MUST refuse."""
    ws = _make_workspace_and_bundle(tmp_path)
    run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    (ws / "hwave" / "output" / "green.npz").write_bytes(
        b"mutated-green-npz",
    )
    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError,
                           match="sha256_hwave_green"):
            run_live_uhf_smoke(
                workspace=ws, uhf_binary="/nonexistent/UHF",
            )
        assert not run_mock.called


def test_marker_write_is_atomic(tmp_path):
    """Important #3 closure: marker uses tmp+fsync+replace, same as
    trans.def substitution. On write failure the target marker is
    unchanged."""
    ws = _make_workspace_and_bundle(tmp_path)
    (ws / "complexuhf" / "trans.def").write_bytes(b"stub")
    marker = run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    original_bytes = marker.read_bytes()
    with patch("os.replace", side_effect=OSError("simulated")):
        with pytest.raises(OSError):
            run_static_bundle_validator(
                fixture_root=tmp_path / "other_fixture",
                workspace=ws,
                expected_phase_mask=(1, 0, 0),
                _skip_bundle_files_for_test=True,
            )
    assert marker.read_bytes() == original_bytes, (
        "marker was clobbered when atomic write failed"
    )


def test_live_smoke_raises_on_uhf_nonzero_exit(tmp_path):
    """Important #4 closure: UHF exiting non-zero raises
    Phase5LiveSmokeGuardError instead of returning silently."""
    import subprocess as _sp
    ws = _make_workspace_and_bundle(tmp_path)
    run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    fake_result = _sp.CompletedProcess(
        args=[], returncode=1, stdout="", stderr="uhf-failed"
    )
    with patch("subprocess.run", return_value=fake_result) as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError,
                           match="returncode=1"):
            run_live_uhf_smoke(
                workspace=ws, uhf_binary="/nonexistent/UHF",
            )
        assert run_mock.called  # UHF was actually launched


def test_live_smoke_refuses_initial_swapped_after_bundle_validation(tmp_path):
    ws = _make_workspace_and_bundle(tmp_path)
    run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    real_validate = phase5_gate_module._validate_runtime_bundle
    swapped = False

    def validate_then_swap(*args, **kwargs):
        nonlocal swapped
        result = real_validate(*args, **kwargs)
        if not swapped:
            (ws / "complexuhf" / "initial.def").write_bytes(b"swapped seed\n")
            swapped = True
        return result

    with patch.object(
        phase5_gate_module,
        "_validate_runtime_bundle",
        side_effect=validate_then_swap,
    ):
        with patch("subprocess.run") as run_mock:
            with pytest.raises(Phase5LiveSmokeGuardError, match="initial"):
                run_live_uhf_smoke(ws, "/nonexistent/UHF")
            assert not run_mock.called


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ("namelist", "namelist_sha256"),
        ("target", "namelist_targets"),
    ],
)
def test_live_smoke_recomputes_new_namelist_marker_fields(
    tmp_path, mutation, match,
):
    ws = _make_workspace_and_bundle(tmp_path)
    run_static_bundle_validator(
        fixture_root=tmp_path / "fake_fixture",
        workspace=ws,
        expected_phase_mask=(0, 0, 0),
        _skip_bundle_files_for_test=True,
    )
    if mutation == "namelist":
        namelist = ws / "complexuhf" / "namelist.def"
        namelist.write_bytes(namelist.read_bytes() + b"\n")
    else:
        locspn = ws / "complexuhf" / "locspn.def"
        locspn.write_bytes(locspn.read_bytes() + b"\n")

    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError, match=match):
            run_live_uhf_smoke(ws, "/nonexistent/UHF")
        assert not run_mock.called


def test_real_output_only_scf_layout_mints_nonempty_marker(tmp_path):
    """Regression for D1: seed-time SCF files live under hwave/output/."""
    fixture_root, workspace = _copy_real_workspace(tmp_path)

    marker = run_static_bundle_validator(
        fixture_root=fixture_root,
        workspace=workspace,
        expected_phase_mask=_REAL_PHASE_MASK,
    )

    raw = marker.read_text()
    payload = json.loads(raw)
    assert raw == json.dumps(payload, sort_keys=True, separators=(",", ":"))
    assert payload["hwave_green_sha256"]
    assert payload["hwave_eigen_sha256"]
    assert payload["hwave_energy_total"]
    assert "fixture_root" not in payload


@pytest.mark.parametrize("name", ["green.npz", "eigen.npz", "energy.dat"])
@pytest.mark.parametrize("failure", ["missing", "empty"])
def test_real_workspace_rejects_missing_or_empty_scf_artifact(
    tmp_path, name, failure,
):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    artifact = workspace / "hwave" / "output" / name
    if failure == "missing":
        artifact.unlink()
    else:
        artifact.write_bytes(b"")

    with pytest.raises(Phase5StaticFailure, match=name):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )
    assert not (workspace / "phase5_static_ok.marker").exists()


def test_real_workspace_rejects_unanchored_nsite_substrings(tmp_path):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    modpara = workspace / "complexuhf" / "modpara.def"
    text = modpara.read_text().replace("Nsite 64", "Nsite 8\nUnrelated 64")
    modpara.write_text(text)

    with pytest.raises(Phase5StaticFailure, match="Nsite"):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


def test_real_workspace_rejects_phase_mask_not_matching_fixture(tmp_path):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    stan = fixture_root / "stan.in"
    stan.write_text(stan.read_text().replace("phase1  = 180.0", "phase1  = 0.0"))

    with pytest.raises(Phase5StaticFailure, match="phase mask"):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


def _remove_one_data_row(path: Path) -> None:
    lines = path.read_text().splitlines(keepends=True)
    if path.name == "geometry.dat":
        del lines[-1]
    else:
        expected_columns = {
            "locspn.def": 2,
            "coulombintra.def": 2,
            "orbitalidx.def": 4,
            "orbitalidxpara.def": 4,
            "trans.def": 6,
        }[path.name]
        for index, line in enumerate(lines):
            tokens = line.split()
            if len(tokens) != expected_columns:
                continue
            try:
                int(tokens[0])
            except ValueError:
                continue
            del lines[index]
            break
        else:  # pragma: no cover - shipped fixture contract
            raise AssertionError(f"no data row found in {path}")
    path.write_text("".join(lines))


@pytest.mark.parametrize(
    "name",
    [
        "geometry.dat",
        "locspn.def",
        "coulombintra.def",
        "orbitalidx.def",
        "orbitalidxpara.def",
        "trans.def",
    ],
)
def test_real_workspace_rejects_each_wrong_bundle_row_count(tmp_path, name):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    bundle_file = workspace / "complexuhf" / name
    _remove_one_data_row(bundle_file)
    if name == "trans.def":
        shutil.copy2(bundle_file, workspace / "mvmc" / "trans.def")

    with pytest.raises(Phase5StaticFailure, match=name):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


def test_real_workspace_rejects_nonhermitian_trans_def(tmp_path):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    trans = workspace / "complexuhf" / "trans.def"
    lines = trans.read_text().splitlines(keepends=True)
    for index, line in enumerate(lines):
        tokens = line.split()
        if len(tokens) == 6:
            tokens[4] = str(float(tokens[4]) + 0.125)
            lines[index] = " ".join(tokens) + "\n"
            break
    altered = "".join(lines)
    trans.write_text(altered)
    # Preserve the authoritative-copy equality so this exercises Hermiticity.
    (workspace / "mvmc" / "trans.def").write_text(altered)

    with pytest.raises(Phase5StaticFailure, match="Hermiticity"):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


def test_real_workspace_rejects_trans_def_different_from_mvmc(tmp_path):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    trans = workspace / "complexuhf" / "trans.def"
    trans.write_bytes(trans.read_bytes() + b"\n")

    with pytest.raises(Phase5StaticFailure, match="byte-equal"):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


@pytest.mark.parametrize("too_small", [1e-6, 1e-5, 0.0])
def test_real_workspace_rejects_perturb_scale_below_g2_floor(
    tmp_path, too_small,
):
    """A seed perturbation at or below the G2 tolerance makes the
    cross-solver density gates passable by a solver that returns the
    seed unchanged. 1e-6 is the value that actually shipped, and it
    equals the G2 tolerance exactly, so it must now be refused.
    """
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    provenance = workspace / "complexuhf" / "initial.def.provenance"
    payload = json.loads(provenance.read_text())
    payload["perturb_scale"] = too_small
    provenance.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")),
    )

    with pytest.raises(Phase5StaticFailure, match="perturb_scale"):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


@pytest.mark.parametrize(
    ("mutation", "match"),
    [
        ("initial_sha", "sha256_initial_def"),
        ("energy", "energy_total"),
        ("missing", "missing field"),
    ],
)
def test_real_workspace_rejects_tampered_provenance(
    tmp_path, mutation, match,
):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    provenance = workspace / "complexuhf" / "initial.def.provenance"
    payload = json.loads(provenance.read_text())
    if mutation == "initial_sha":
        payload["sha256_initial_def"] = "0" * 64
    elif mutation == "energy":
        payload["energy_total"] += 1.0
    else:
        del payload["sha256_hwave_green"]
    provenance.write_text(
        json.dumps(payload, sort_keys=True, separators=(",", ":")),
    )

    with pytest.raises(Phase5StaticFailure, match=match):
        run_static_bundle_validator(
            fixture_root=fixture_root,
            workspace=workspace,
            expected_phase_mask=_REAL_PHASE_MASK,
        )


@pytest.mark.parametrize(
    "field",
    [
        "expected_phase_mask",
        "namelist_sha256",
        "namelist_targets",
        "trans_def_sha256",
        "initial_def_provenance_sha256",
        "hwave_green_sha256",
        "hwave_eigen_sha256",
        "hwave_energy_total",
    ],
)
def test_live_smoke_refuses_each_blanked_marker_field(tmp_path, field):
    fixture_root, workspace = _copy_real_workspace(tmp_path)
    marker = run_static_bundle_validator(
        fixture_root=fixture_root,
        workspace=workspace,
        expected_phase_mask=_REAL_PHASE_MASK,
    )
    payload = json.loads(marker.read_text())
    payload[field] = ""
    marker.write_text(json.dumps(payload, sort_keys=True, separators=(",", ":")))

    with patch("subprocess.run") as run_mock:
        with pytest.raises(Phase5LiveSmokeGuardError, match=field):
            run_live_uhf_smoke(workspace, "/nonexistent/UHF")
        assert not run_mock.called
