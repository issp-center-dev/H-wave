"""Snapshot-workspace rejection guard for compare.py and the G4 guard.

Both ``compare.py --workspace <ws> --mode <mode>`` and
``soc_apbc_topology_guard.py --workspace <ws> --mode g4 --composite-manifest ...``
MUST refuse to consume any workspace whose canonical path (or the
canonical path of a shipped-in artifact under it) lives under
``tests/data``. Otherwise a hostile / careless invocation could point
either helper at the unit-test snapshot fixtures and mint a green
anchored PASS record without exercising the producer chain.

The test matrix covers:

1. workspace-root literal ``tests/data``  → reject
2. workspace-root literal ``tests/data/<file>`` (nonsense subdir under
   the snapshot root)                                         → reject
3. workspace-root symlink to ``tests/data``                    → reject
4. artifact-inside-workspace symlink to ``tests/data/<npz>``   → reject

Each case runs against both compare.py (mode=g0-writer-check) and the
topology_guard (mode=g4). Total: 4 * 2 = 8 negative tests.

See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst`` for the validation
contract.
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


_REPO_ROOT = Path(__file__).resolve().parents[1]
_TESTS_DATA_ROOT = _REPO_ROOT / "tests" / "data"
_COMPARE_PATH = str(
    _REPO_ROOT / "tests" / "validation" / "uhfk_mvmc_pairproduct" / "compare.py"
)
_GUARD_PATH = str(
    _REPO_ROOT / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "scripts" / "soc_apbc_topology_guard.py"
)
_REAL_CASE_SOURCE = (
    _REPO_ROOT / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "case_soc_rashba_2d_sub_apbc"
)
_REAL_CASE_DIR = str(_REAL_CASE_SOURCE)
_REAL_MANIFEST = os.path.join(_REAL_CASE_DIR, "composite_element.json")


@pytest.fixture(scope="module", autouse=True)
def _real_case_workspace(tmp_path_factory):
    """Build the positive G4 workspace from tracked snapshots."""
    global _REAL_CASE_DIR, _REAL_MANIFEST

    workspace = (
        tmp_path_factory.mktemp("snapshot-guard-real-case")
        / "case_soc_rashba_2d_sub_apbc"
    )
    output = workspace / "output"
    output.mkdir(parents=True)
    shutil.copy2(_REAL_CASE_SOURCE / "composite_element.json", workspace)
    for name in ("green", "eigen", "occupation"):
        shutil.copy2(
            _TESTS_DATA_ROOT
            / f"v36_case_soc_rashba_2d_sub_apbc_{name}.npz",
            output / f"{name}.npz",
        )

    original_case_dir = _REAL_CASE_DIR
    original_manifest = _REAL_MANIFEST
    _REAL_CASE_DIR = str(workspace)
    _REAL_MANIFEST = str(workspace / "composite_element.json")
    yield
    _REAL_CASE_DIR = original_case_dir
    _REAL_MANIFEST = original_manifest


def _run(cmd):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(_REPO_ROOT / "src"), str(_REPO_ROOT)]
    )
    return subprocess.run(
        cmd, capture_output=True, text=True, env=env, cwd=str(_REPO_ROOT),
    )


def _run_compare(workspace, mode="g0-writer-check"):
    return _run([
        sys.executable, _COMPARE_PATH,
        "--workspace", str(workspace), "--mode", mode,
    ])


def _run_guard(workspace, manifest=None):
    if manifest is None:
        manifest = _REAL_MANIFEST
    return _run([
        sys.executable, _GUARD_PATH,
        "--workspace", str(workspace),
        "--mode", "g4",
        "--composite-manifest", manifest,
    ])


def _assert_rejected(res, needle):
    assert res.returncode == 2, (
        f"expected exit code 2; got {res.returncode}. "
        f"stdout={res.stdout!r} stderr={res.stderr!r}"
    )
    assert res.stdout == "", (
        f"expected empty stdout (no PASS record); got stdout={res.stdout!r}"
    )
    assert needle in res.stderr, (
        f"expected {needle!r} substring in stderr; got {res.stderr!r}"
    )


# ---------------------------------------------------------------------
# compare.py rejections
# ---------------------------------------------------------------------


def test_compare_rejects_workspace_literal_tests_data():
    """compare.py --workspace tests/data --mode g0-writer-check must
    reject before running any dispatcher."""
    res = _run_compare(_TESTS_DATA_ROOT)
    _assert_rejected(res, "resolves under tests/data")


def test_compare_rejects_workspace_literal_tests_data_subdir(tmp_path):
    """Literal subdirectory under tests/data (created for the test)
    must be rejected as well."""
    subdir = _TESTS_DATA_ROOT / "_phase3e_probe_compare"
    subdir.mkdir(exist_ok=True)
    try:
        res = _run_compare(subdir)
        _assert_rejected(res, "resolves under tests/data")
    finally:
        subdir.rmdir()


def test_compare_rejects_workspace_symlinked_to_tests_data(tmp_path):
    """When --workspace is a symlink whose target lives under tests/data,
    canonicalization must catch it."""
    link = tmp_path / "ws_symlink"
    link.symlink_to(_TESTS_DATA_ROOT, target_is_directory=True)
    res = _run_compare(link)
    _assert_rejected(res, "resolves under tests/data")


def test_compare_rejects_artifact_symlinked_to_tests_data(tmp_path):
    """Workspace root is a benign tmp dir but ``hwave/green.npz`` is a
    symlink into tests/data; the per-artifact resolve MUST catch this."""
    (tmp_path / "hwave").mkdir()
    real_snap = (
        _TESTS_DATA_ROOT / "v36_case_soc_rashba_2d_sub_apbc_green.npz"
    )
    (tmp_path / "hwave" / "green.npz").symlink_to(real_snap)
    res = _run_compare(tmp_path)
    _assert_rejected(res, "hwave/green.npz")


# ---------------------------------------------------------------------
# topology_guard.py rejections
# ---------------------------------------------------------------------


def test_guard_rejects_workspace_literal_tests_data():
    res = _run_guard(_TESTS_DATA_ROOT)
    _assert_rejected(res, "resolves under tests/data")


def test_guard_rejects_workspace_literal_tests_data_subdir():
    subdir = _TESTS_DATA_ROOT / "_phase3e_probe_guard"
    subdir.mkdir(exist_ok=True)
    try:
        res = _run_guard(subdir)
        _assert_rejected(res, "resolves under tests/data")
    finally:
        subdir.rmdir()


def test_guard_rejects_workspace_symlinked_to_tests_data(tmp_path):
    link = tmp_path / "ws_symlink"
    link.symlink_to(_TESTS_DATA_ROOT, target_is_directory=True)
    res = _run_guard(link)
    _assert_rejected(res, "resolves under tests/data")


def test_guard_rejects_artifact_symlinked_to_tests_data(tmp_path):
    (tmp_path / "hwave").mkdir()
    real_snap = (
        _TESTS_DATA_ROOT / "v36_case_soc_rashba_2d_sub_apbc_green.npz"
    )
    (tmp_path / "hwave" / "green.npz").symlink_to(real_snap)
    # Copy manifest into the workspace so argparse succeeds; the guard
    # must still reject on the artifact resolution.
    manifest_dst = tmp_path / "composite_element.json"
    with open(_REAL_MANIFEST) as fp:
        manifest = json.load(fp)
    with open(manifest_dst, "w") as fp:
        json.dump(manifest, fp)
    res = _run_guard(tmp_path, manifest=str(manifest_dst))
    _assert_rejected(res, "hwave/green.npz")


# ---------------------------------------------------------------------
# Positive baseline: the real fixture is NOT under tests/data and MUST
# still pass through the guard (topology_guard) after this hardening.
# ---------------------------------------------------------------------


def test_guard_passes_real_fixture_after_snapshot_guard_landed():
    """The positive baseline still works. If snapshot rejection
    is over-broad and starts refusing legitimate case dirs, this trips."""
    res = _run_guard(_REAL_CASE_DIR)
    assert res.returncode == 0, (
        f"unexpected reject on real fixture; stderr={res.stderr!r}"
    )
    assert res.stdout.startswith("G4 PASS mode=g4 "), res.stdout


# ---------------------------------------------------------------------
# Non-repository CWD regression: both helpers must locate ``tests/data``
# relative to their modules, even when invoked from elsewhere.
# ---------------------------------------------------------------------


def _run_from_non_repo_cwd(cmd):
    """Invoke the CLI from /tmp so its CWD is not the repository."""
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(_REPO_ROOT / "src"), str(_REPO_ROOT)]
    )
    return subprocess.run(
        cmd, capture_output=True, text=True, env=env, cwd="/tmp",
    )


def test_compare_rejects_tests_data_when_run_from_non_repo_cwd():
    """An absolute ``tests/data`` workspace must be rejected from /tmp."""
    res = _run_from_non_repo_cwd([
        sys.executable, _COMPARE_PATH,
        "--workspace", str(_TESTS_DATA_ROOT),
        "--mode", "g2b",
    ])
    _assert_rejected(res, "resolves under tests/data")


def test_guard_rejects_tests_data_when_run_from_non_repo_cwd():
    """Same for the G4 topology guard."""
    res = _run_from_non_repo_cwd([
        sys.executable, _GUARD_PATH,
        "--workspace", str(_TESTS_DATA_ROOT),
        "--mode", "g4",
        "--composite-manifest", _REAL_MANIFEST,
    ])
    _assert_rejected(res, "resolves under tests/data")


def test_compare_rejects_artifact_symlink_from_non_repo_cwd(tmp_path):
    """The per-artifact resolve must also survive a non-repo CWD."""
    (tmp_path / "hwave").mkdir()
    real_snap = (
        _TESTS_DATA_ROOT / "v36_case_soc_rashba_2d_sub_apbc_green.npz"
    )
    (tmp_path / "hwave" / "green.npz").symlink_to(real_snap)
    res = _run_from_non_repo_cwd([
        sys.executable, _COMPARE_PATH,
        "--workspace", str(tmp_path), "--mode", "g0-writer-check",
    ])
    _assert_rejected(res, "hwave/green.npz")
