"""Mandatory negative tests for the semantic G4 topology guard.

Each test constructs a small stub manifest / workspace and
asserts the guard rejects with exit code 2 + empty stdout on the
targeted invariant violation.
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest


_GUARD_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "validation",
        "uhfk_mvmc_pairproduct",
        "scripts",
        "soc_apbc_topology_guard.py",
    )
)
_REPO_ROOT = Path(__file__).resolve().parents[1]
_DATA_DIR = _REPO_ROOT / "tests" / "data"
_CASE_ROOT = (
    _REPO_ROOT / "tests" / "validation" / "uhfk_mvmc_pairproduct"
)
_REAL_CASE_SOURCE = _CASE_ROOT / "case_soc_rashba_2d_sub_apbc"
_REAL_CASE_DIR = str(_REAL_CASE_SOURCE)
# Fixture whose active APBC axes are (y, z) -- NO x. Used to
# regression-test that the sub_offset-differs check gates on the
# fixture's actual active axes rather than a hardcoded x. See
# ``docs/en/source/algorithm/uhfk_to_mvmc.rst`` for gauge composition.
_REAL_CASE_SOURCE_YZ = _CASE_ROOT / "case_soc_rashba_3d_sub_apbc_yz"
_REAL_CASE_DIR_YZ = str(_REAL_CASE_SOURCE_YZ)
_REAL_CASE_SOURCE_XY = _CASE_ROOT / "case_soc_rashba_3d_sub_apbc_xy"
_REAL_CASE_DIR_XY = str(_REAL_CASE_SOURCE_XY)


@pytest.fixture(scope="module", autouse=True)
def _real_case_workspaces(tmp_path_factory):
    """Build clean-checkout G4 workspaces from tracked snapshots."""
    global _REAL_CASE_DIR, _REAL_CASE_DIR_YZ, _REAL_CASE_DIR_XY

    workspace_root = tmp_path_factory.mktemp("g4-real-cases")
    case_specs = (
        (
            _REAL_CASE_SOURCE,
            "v36_case_soc_rashba_2d_sub_apbc",
            "case_soc_rashba_2d_sub_apbc",
        ),
        (
            _REAL_CASE_SOURCE_YZ,
            "v37_case_soc_rashba_3d_sub_apbc_yz",
            "case_soc_rashba_3d_sub_apbc_yz",
        ),
        (
            _REAL_CASE_SOURCE_XY,
            "v37_case_soc_rashba_3d_sub_apbc_xy",
            "case_soc_rashba_3d_sub_apbc_xy",
        ),
    )
    workspaces = []
    for source, snapshot_prefix, case_name in case_specs:
        workspace = workspace_root / case_name
        output = workspace / "output"
        output.mkdir(parents=True)
        shutil.copy2(source / "composite_element.json", workspace)
        for name in ("green", "eigen", "occupation"):
            shutil.copy2(
                _DATA_DIR / f"{snapshot_prefix}_{name}.npz",
                output / f"{name}.npz",
            )
        workspaces.append(str(workspace))

    original_dirs = (_REAL_CASE_DIR, _REAL_CASE_DIR_YZ, _REAL_CASE_DIR_XY)
    _REAL_CASE_DIR, _REAL_CASE_DIR_YZ, _REAL_CASE_DIR_XY = workspaces
    yield
    _REAL_CASE_DIR, _REAL_CASE_DIR_YZ, _REAL_CASE_DIR_XY = original_dirs


def _load_guard_module():
    spec = importlib.util.spec_from_file_location(
        "_v36_topology_guard_under_test", _GUARD_PATH,
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def guard_mod():
    return _load_guard_module()


def _load_real_manifest():
    manifest_path = os.path.join(_REAL_CASE_DIR, "composite_element.json")
    with open(manifest_path) as fp:
        return json.load(fp)


def _run_cli(workspace, manifest_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = "src:."
    return subprocess.run(
        [
            sys.executable, _GUARD_PATH,
            "--workspace", workspace,
            "--mode", "g4",
            "--composite-manifest", manifest_path,
        ],
        capture_output=True, text=True, env=env,
    )


def _write_manifest(path, manifest):
    with open(path, "w") as fp:
        json.dump(manifest, fp, indent=2, sort_keys=True)
        fp.write("\n")


def _make_stub_workspace(tmp_path, green_source):
    """Create a workspace with only hwave/green.npz linked from a real
    SCF output. Callers can then override the manifest or delete files
    to trigger specific negative paths."""
    hwave = tmp_path / "hwave"
    hwave.mkdir(parents=True, exist_ok=True)
    dst = hwave / "green.npz"
    dst.write_bytes(green_source.read_bytes())
    return str(tmp_path)


# ---------------------------------------------------------------------
# 1. Missing composite (magnitude falls below 0.8 * G_c_abs)
# ---------------------------------------------------------------------


def test_g4_rejects_missing_composite_element(tmp_path):
    """Manifest points at a composite element whose magnitude on the
    fresh workspace falls below 0.8 * G_c_abs (in this test, we bump
    G_c_abs 10x higher than the real value so no current-run element
    meets the floor)."""
    ws = _make_stub_workspace(
        tmp_path, green_source=(
            _load_real_manifest_dir_green_file()
        ),
    )
    manifest = _load_real_manifest()
    manifest["G_c_abs"] = 10.0 * manifest["G_c_abs"]
    expected_threshold = max(1e-5, 0.10 * manifest["G_c_abs"])
    for mutator_id in manifest["T_M_per_mutation"]:
        manifest["T_M_per_mutation"][mutator_id] = expected_threshold
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "composite element on current workspace" in res.stderr


# ---------------------------------------------------------------------
# 2. Same-spin composite in manifest.
# ---------------------------------------------------------------------


def test_g4_rejects_same_spin_composite(tmp_path):
    ws = _make_stub_workspace(
        tmp_path, green_source=_load_real_manifest_dir_green_file(),
    )
    manifest = _load_real_manifest()
    manifest["s_c"] = manifest["t_c"] = 0
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "cross-spin" in res.stderr


# ---------------------------------------------------------------------
# 3. Manifest sub_offset_x does not differ between i and j.
# ---------------------------------------------------------------------


def test_g4_rejects_no_sub_offset_diff(tmp_path):
    ws = _make_stub_workspace(
        tmp_path, green_source=_load_real_manifest_dir_green_file(),
    )
    manifest = _load_real_manifest()
    manifest["sub_offset_i"] = [0, 0, 0]
    manifest["sub_offset_j"] = [0, 1, 0]  # same x, different y
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "sub_offset_x differs" in res.stderr


# ---------------------------------------------------------------------
# 3b. Regression: the sub_offset-differs check must gate on the
# fixture's ACTUAL active APBC axis, not a hardcoded x. A manifest
# whose only active axis is y (theta = (0, pi, 0)) with sub_offset
# differing in x but NOT in y must still be rejected. This pins the
# requirement against a hardcoded axis-0 check. See
# ``docs/en/source/algorithm/uhfk_to_mvmc.rst`` for gauge composition.
# ---------------------------------------------------------------------


def test_g4_rejects_no_sub_offset_diff_on_non_x_active_axis(tmp_path):
    ws = _make_stub_workspace(
        tmp_path, green_source=_load_real_manifest_dir_green_file(),
    )
    manifest = _load_real_manifest()
    manifest["theta_radians"] = [0.0, 3.141592653589793, 0.0]  # y-only APBC
    manifest["sub_offset_i"] = [0, 0, 0]
    manifest["sub_offset_j"] = [1, 0, 0]  # x differs, active-axis y does not
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "sub_offset_y differs" in res.stderr


# ---------------------------------------------------------------------
# 4. Manifest G_c_abs below the 1e-3 magnitude floor.
# ---------------------------------------------------------------------


def test_g4_rejects_low_magnitude_composite(tmp_path):
    ws = _make_stub_workspace(
        tmp_path, green_source=_load_real_manifest_dir_green_file(),
    )
    manifest = _load_real_manifest()
    manifest["G_c_abs"] = 5e-4  # below 1e-3
    expected_threshold = max(1e-5, 0.10 * manifest["G_c_abs"])
    for mutator_id in manifest["T_M_per_mutation"]:
        manifest["T_M_per_mutation"][mutator_id] = expected_threshold
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "magnitude floor" in res.stderr


# ---------------------------------------------------------------------
# 5. Mutation delta below the policy-derived T_M.
# ---------------------------------------------------------------------


def test_g4_rejects_mutation_below_T_M(tmp_path):
    """Keep the manifest policy valid, but construct a current-run green
    tensor whose target element has only folded k=0 support. M-gauge-4
    changes only the folded phase, so its delta is exactly zero while the
    composite magnitude and the preceding gauge mutations remain strong."""
    ws = _make_stub_workspace(
        tmp_path, green_source=_load_real_manifest_dir_green_file(),
    )
    manifest = _load_real_manifest()
    with np.load(_load_real_manifest_dir_green_file()) as green_npz:
        synthetic = np.zeros_like(green_npz["green_sublattice"])
    folded_orb_i = (
        manifest["sub_offset_i"][0]
        + manifest["sub_shape"][0] * (
            manifest["sub_offset_i"][1]
            + manifest["sub_shape"][1] * manifest["sub_offset_i"][2]
        )
    )
    folded_orb_j = (
        manifest["sub_offset_j"][0]
        + manifest["sub_shape"][0] * (
            manifest["sub_offset_j"][1]
            + manifest["sub_shape"][1] * manifest["sub_offset_j"][2]
        )
    )
    aa = 2 * folded_orb_i + manifest["s_c"]
    bb = 2 * folded_orb_j + manifest["t_c"]
    synthetic[:, 0, aa, 0, bb] = manifest["G_c_abs"]
    np.savez(
        os.path.join(ws, "hwave", "green.npz"),
        green_sublattice=synthetic,
    )
    m_path = os.path.join(ws, "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(ws, m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "M-gauge-4" in res.stderr
    assert "T_M" in res.stderr


@pytest.mark.parametrize(
    ("mutator_id", "threshold", "error_fragment"),
    [
        ("M-gauge-1-x", float("nan"), "finite"),
        ("M-gauge-1-x", 1e-300, "threshold policy"),
        ("M-gauge-1-x", 0.0, "threshold policy"),
        ("M-gauge-1-z", 1e-5, "inactive"),
    ],
)
def test_g4_rejects_manifest_threshold_policy_drift(
    tmp_path, mutator_id, threshold, error_fragment,
):
    manifest_path = os.path.join(_REAL_CASE_DIR_XY, "composite_element.json")
    with open(manifest_path) as fp:
        manifest = json.load(fp)
    manifest["T_M_per_mutation"][mutator_id] = threshold
    drifted_path = tmp_path / "manifest.json"
    _write_manifest(drifted_path, manifest)

    res = _run_cli(_REAL_CASE_DIR_XY, str(drifted_path))

    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert error_fragment in res.stderr.lower()


# ---------------------------------------------------------------------
# 6. Guard reads CURRENT workspace, NOT manifest-cached values.
# ---------------------------------------------------------------------


def test_g4_reads_current_workspace_not_manifest_alone(tmp_path):
    """Swap the workspace's green.npz for a random-noise green_sublattice
    of the same shape. The guard MUST reject because the recomputed
    G_current at the manifest coordinates no longer matches the
    manifest's cached G_c_abs; this proves the guard is reading the
    CURRENT workspace and not trusting the manifest values alone."""
    manifest = _load_real_manifest()
    real_green_path = os.path.join(_REAL_CASE_DIR, "output", "green.npz")
    with np.load(real_green_path) as gz:
        real_gs = gz["green_sublattice"]
    # Random-noise green_sublattice of the same shape.
    rng = np.random.default_rng(31337)
    noise = (
        rng.standard_normal(real_gs.shape).astype(np.complex128)
        + 1j * rng.standard_normal(real_gs.shape).astype(np.complex128)
    ) * 1e-4  # small so the composite magnitude collapses
    hwave = tmp_path / "hwave"
    hwave.mkdir(parents=True, exist_ok=True)
    np.savez(hwave / "green.npz", green_sublattice=noise)
    m_path = str(tmp_path / "manifest.json")
    _write_manifest(m_path, manifest)
    res = _run_cli(str(tmp_path), m_path)
    assert res.returncode == 2, res.stderr
    assert res.stdout == "", res.stdout
    assert "composite element on current workspace" in res.stderr


# ---------------------------------------------------------------------
# Positive path (baseline): guard passes on the real fixture.
# ---------------------------------------------------------------------


def test_g4_passes_on_real_case_soc_rashba_2d_sub_apbc(tmp_path):
    """Sanity: the guard MUST PASS on the committed fixture with its own
    committed composite_element.json. If this fails, all six negative
    tests above are trivially satisfied for the wrong reason."""
    manifest_path = os.path.join(_REAL_CASE_DIR, "composite_element.json")
    res = _run_cli(_REAL_CASE_DIR, manifest_path)
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("G4 PASS mode=g4 "), res.stdout
    assert "artifact_source=hwave+bridge+composite-manifest" in res.stdout
    assert "helper=soc_apbc_topology_guard" in res.stdout


def test_g4_passes_on_real_case_soc_rashba_3d_sub_apbc_yz(tmp_path):
    """The guard MUST PASS on the committed yz fixture, whose active
    APBC axes are (y, z) with
    NO x. This is the concrete case that exposed the hardcoded-x
    sub_offset bug: the composite this fixture's own producer selects
    has sub_offset differing in y and z but NOT in x (x is not an
    APBC direction here, so the active-axis condition places no requirement
    on it), which a hardcoded ``sub_offset_x`` check rejects even
    though the fixture is valid. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    manifest_path = os.path.join(_REAL_CASE_DIR_YZ, "composite_element.json")
    res = _run_cli(_REAL_CASE_DIR_YZ, manifest_path)
    assert res.returncode == 0, res.stderr
    assert res.stdout.startswith("G4 PASS mode=g4 "), res.stdout
    assert "artifact_source=hwave+bridge+composite-manifest" in res.stdout
    assert "helper=soc_apbc_topology_guard" in res.stdout


# ---------------------------------------------------------------------
# Helper: load the real fixture's green.npz for cloning.
# ---------------------------------------------------------------------


def _load_real_manifest_dir_green_file():
    return Path(_REAL_CASE_DIR) / "output" / "green.npz"
