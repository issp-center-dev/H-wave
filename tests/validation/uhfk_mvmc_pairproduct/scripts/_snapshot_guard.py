"""Shared snapshot-workspace rejection guard.

Both ``compare.py`` (seven-gate dispatcher) and
``soc_apbc_topology_guard.py`` (G4 helper) refuse to consume any
workspace whose canonical path (or the canonical path of a shipped-in
artifact under it) lives under ``tests/data/``. That directory holds
unit-test snapshots, not live producer output; running a real E2E gate
against them would silently give a green record without exercising the
producer chain.

Contract:

- Resolve the workspace with ``Path.resolve(strict=True)``. If it lives
  under any of the reachable ``tests/data`` roots, raise
  :class:`SnapshotWorkspaceRejected`.
- For each expected artifact under the workspace, resolve the file's
  own path. If it lives under ``tests/data`` (via a symlink chain the
  workspace root itself would not catch), raise the same exception.
- The caller converts the exception into exit code 2 + empty stdout so
  the anchored PASS record cannot be printed.

See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import os
from pathlib import Path
from typing import Iterable


class SnapshotWorkspaceRejected(RuntimeError):
    """Raised when a workspace or one of its artifacts canonicalizes
    under ``tests/data/``. Callers translate this into exit code 2."""


# Resolve the ``tests/data`` root from this module's location, not the
# process CWD, so invocation from elsewhere cannot bypass the snapshot check.
_SCRIPTS_DIR = Path(__file__).resolve().parent
_REPO_ROOT_FROM_SCRIPT = _SCRIPTS_DIR.parents[3]  # scripts/../../../..
_TESTS_DATA_ROOTS_ABS = (_REPO_ROOT_FROM_SCRIPT / "tests" / "data",)


# Artifacts the gate helpers read, relative to the workspace. If
# any of these paths (after full realpath) canonicalize under
# ``tests/data``, the helper must refuse the workspace even if the
# workspace root itself is a benign tmp directory.
_DEFAULT_ARTIFACTS = (
    "hwave/green.npz",
    "hwave/eigen.npz",
    "hwave/occupation.npz",
    "output/green.npz",
    "output/eigen.npz",
    "output/occupation.npz",
    "bridge/F_pre_noise.npz",
    "bridge/F_post_aggregate.npz",
    "bridge_zeronoise/F_pre_noise.npz",
    "bridge_zeronoise/F_post_aggregate.npz",
    "complexuhf/zvo_UHF_cisajs.dat",
    "mvmc/zvo_out_selected.dat",
    "composite_element.json",
)


def _resolved_tests_data_roots():
    """Resolve every absolute tests/data root anchored at THIS module's
    location. Independent of the caller's CWD.

    Fail-closed contract: if the anchored root does not exist on disk
    (e.g. the module was shipped without the tests tree — extremely
    unusual), the guard MUST still refuse the workspace rather than
    silently allowing bypass. Callers translate `SnapshotWorkspaceRejected`
    into exit code 2.
    """
    roots = []
    for anchored in _TESTS_DATA_ROOTS_ABS:
        try:
            roots.append(anchored.resolve(strict=True))
        except FileNotFoundError:
            raise SnapshotWorkspaceRejected(
                f"snapshot_guard: anchored tests/data root {anchored} "
                "does not exist; refusing to run without a reachable "
                "reference root."
            )
    return roots


def _is_under(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def reject_snapshot_workspace(
    workspace,
    artifacts: Iterable[str] = _DEFAULT_ARTIFACTS,
    helper_name: str = "snapshot_guard",
):
    """Canonicalize ``workspace`` and every existing artifact under it;
    raise :class:`SnapshotWorkspaceRejected` if any lands under a
    ``tests/data`` root reachable from CWD.

    Parameters
    ----------
    workspace : str | os.PathLike
        Directory the caller wants to consume as a live workspace.
    artifacts : iterable of str
        Workspace-relative artifact paths to canonicalize individually.
        A per-artifact check catches the symlink-inside-benign-workspace
        pattern (e.g. ``/tmp/ws/hwave/green.npz -> tests/data/...``)
        that the workspace-root check alone would miss.
    helper_name : str
        Included in the exception message so callers can trace which
        gate helper refused which workspace.
    """
    roots = _resolved_tests_data_roots()

    ws_path = Path(workspace)
    try:
        ws_resolved = ws_path.resolve(strict=True)
    except FileNotFoundError as exc:
        raise SnapshotWorkspaceRejected(
            f"{helper_name}: workspace {workspace} does not exist"
        ) from exc
    for root in roots:
        if _is_under(ws_resolved, root):
            raise SnapshotWorkspaceRejected(
                f"{helper_name}: workspace {ws_resolved} resolves under "
                f"tests/data root {root}; refusing to consume unit "
                "snapshots as live gate input."
            )

    for artifact_rel in artifacts:
        artifact_path = ws_path / artifact_rel
        # os.path.realpath follows symlinks even when the target does
        # not exist under strict=True, so we use it for the artifact
        # check.
        if not artifact_path.exists() and not artifact_path.is_symlink():
            continue
        try:
            resolved = Path(os.path.realpath(str(artifact_path)))
        except OSError:
            continue
        for root in roots:
            if _is_under(resolved, root):
                raise SnapshotWorkspaceRejected(
                    f"{helper_name}: artifact {artifact_rel} resolves "
                    f"to {resolved} under tests/data root {root}; "
                    "refusing to consume unit snapshots as live gate "
                    "input."
                )
