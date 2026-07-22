"""initial.def provenance producer contract tests."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess

import pytest


_REPO_ROOT = Path(__file__).resolve().parents[1]
_SEED_SCRIPT = (
    _REPO_ROOT
    / "tests/validation/uhfk_mvmc_pairproduct/scripts/seed_complexuhf_from_hwave.py"
)
_RUN_SH = (
    _REPO_ROOT / "tests/validation/uhfk_mvmc_pairproduct/run.sh"
)


def _load_seed_module():
    spec = importlib.util.spec_from_file_location(
        "seed_complexuhf_from_hwave", _SEED_SCRIPT,
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _phase5_mask_dispatch_block() -> str:
    text = _RUN_SH.read_text()
    start = text.index("  # Static/live gates apply to the 3D APBC cases.")
    end = text.index("  # normalize_mvmc_output", start)
    return text[start:end]


def _run_phase5_mask_dispatch(tmp_path, case_name):
    harness = tmp_path / "phase5_mask_dispatch.sh"
    harness.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        f"{_phase5_mask_dispatch_block()}\n"
        'printf "MASK=%s\\n" "${PHASE5_EXPECTED_MASK}"\n'
    )
    env = os.environ.copy()
    env["CASE"] = case_name
    return subprocess.run(
        ["bash", str(harness)],
        capture_output=True,
        text=True,
        env=env,
        cwd=tmp_path,
    )


def test_phase5_mask_dispatch_rejects_unknown_v37_case(tmp_path):
    result = _run_phase5_mask_dispatch(
        tmp_path, "case_soc_rashba_3d_sub_apbc_future",
    )

    assert result.returncode == 1
    assert "unknown case name" in result.stderr
    assert "define PHASE5_EXPECTED_MASK" in result.stderr


@pytest.mark.parametrize(
    ("case_name", "expected_mask"),
    [
        ("case_soc_rashba_3d_sub_apbc_xy", "1, 1, 0"),
        ("case_soc_rashba_3d_sub_apbc_xz", "1, 0, 1"),
        ("case_soc_rashba_3d_sub_apbc_yz", "0, 1, 1"),
        ("case_soc_rashba_3d_sub_apbc_xyz", "1, 1, 1"),
        ("case_soc_rashba_2d_sub_apbc", ""),
    ],
)
def test_phase5_mask_dispatch_preserves_supported_cases(
    tmp_path, case_name, expected_mask,
):
    result = _run_phase5_mask_dispatch(tmp_path, case_name)

    assert result.returncode == 0, result.stderr
    assert result.stdout == f"MASK={expected_mask}\n"


def test_phase5_provenance_is_canonical_and_uses_anchored_energy(tmp_path):
    seed = _load_seed_module()
    hwave_output = tmp_path / "hwave" / "output"
    complexuhf = tmp_path / "complexuhf"
    hwave_output.mkdir(parents=True)
    complexuhf.mkdir()

    initial_def = complexuhf / "initial.def"
    green = hwave_output / "green.npz"
    eigen = hwave_output / "eigen.npz"
    initial_def.write_bytes(b"initial-def\n")
    green.write_bytes(b"green-npz")
    eigen.write_bytes(b"eigen-npz")
    (hwave_output / "energy.dat").write_text(
        "FreeEnergy_Total = -9.5\n"
        "  Energy_Total = -1.25e+00\n"
    )

    sidecar = seed._write_provenance(
        initial_def, tmp_path / "hwave", perturb_scale=1e-6,
    )

    raw = sidecar.read_text()
    payload = json.loads(raw)
    assert raw == json.dumps(payload, sort_keys=True, separators=(",", ":"))
    assert set(payload) == {
        "sha256_initial_def",
        "sha256_hwave_green",
        "sha256_hwave_eigen",
        "energy_total",
        "perturb_scale",
    }
    assert payload["sha256_initial_def"] == hashlib.sha256(
        initial_def.read_bytes(),
    ).hexdigest()
    assert payload["sha256_hwave_green"] == hashlib.sha256(
        green.read_bytes(),
    ).hexdigest()
    assert payload["sha256_hwave_eigen"] == hashlib.sha256(
        eigen.read_bytes(),
    ).hexdigest()
    assert payload["energy_total"] == -1.25
    assert payload["perturb_scale"] == 1e-6


def _phase5_launch_dispatch_block() -> str:
    text = _RUN_SH.read_text()
    seed = text.index(
        'python3 "${HERE}/scripts/seed_complexuhf_from_hwave.py"'
    )
    start = text.index(
        '    if [[ -n "${PHASE5_EXPECTED_MASK}" ]]; then', seed,
    )
    end = text.index(
        "    # Residual SCF-log sanity check.", start,
    )
    return text[start:end]


@pytest.mark.parametrize(
    ("phase5_mask", "expected_message", "excluded_message"),
    [
        (
            "1, 1, 0",
            "Gated ComplexUHF run failed (gate check or solver); tail uhf.log:",
            "ComplexUHF UHF failed; tail uhf.log:",
        ),
        (
            "",
            "ComplexUHF UHF failed; tail uhf.log:",
            "Gated ComplexUHF run failed (gate check or solver); tail uhf.log:",
        ),
    ],
)
def test_phase5_launch_dispatch_distinguishes_gate_and_solver_failures(
    tmp_path, phase5_mask, expected_message, excluded_message,
):
    harness = tmp_path / "phase5_launch_dispatch.sh"
    harness.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "python3() { return 9; }\n"
        "run_native() { return 9; }\n"
        'ROOT="unused"\n'
        'COMPLEXUHF_CASE="unused"\n'
        'WORK="unused"\n'
        'UHF="unused"\n'
        f'PHASE5_EXPECTED_MASK="{phase5_mask}"\n'
        f"{_phase5_launch_dispatch_block()}\n"
    )

    result = subprocess.run(
        ["bash", str(harness)],
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )

    assert result.returncode == 1
    assert expected_message in result.stderr
    assert excluded_message not in result.stderr


def test_phase5_run_sh_uses_sole_entry_point_and_preserves_v36_launch():
    text = _RUN_SH.read_text()
    phase5_block = text[
        text.index('python3 "${HERE}/scripts/seed_complexuhf_from_hwave.py"'):
        text.index("scf_steps=$(awk")
    ]

    assert "run_phase5_gate" in phase5_block
    assert "run_static_bundle_validator" not in phase5_block
    assert "run_live_uhf_smoke" not in phase5_block
    assert phase5_block.count("from tools._uhfk_to_mvmc.phase5_gate import") == 1
    assert (
        'echo "Gated ComplexUHF run failed (gate check or solver); tail uhf.log:" >&2'
        in phase5_block
    )
    assert (
        '    else\n'
        '      # Preserve the established launch path byte-for-byte.\n'
        '      run_native "${UHF}" namelist.def > uhf.log 2>&1 || {\n'
        '        echo "ComplexUHF UHF failed; tail uhf.log:" >&2\n'
        '        tail -50 uhf.log >&2\n'
        '        exit 1\n'
        '      }\n'
        '    fi\n'
    ) in phase5_block


def _complexuhf_preparation_block() -> str:
    text = _RUN_SH.read_text()
    staging = text.index("# BEGIN task6a2-authoritative-bundle-staging")
    start = text.rfind(
        '  (\n    cd "${COMPLEXUHF_WORK}"',
        0,
        staging,
    )
    end = text.index(
        'python3 "${HERE}/scripts/seed_complexuhf_from_hwave.py"', staging,
    )
    return text[start:end]


def _run_complexuhf_bundle_staging(tmp_path, workspace_file):
    block = _complexuhf_preparation_block()
    fixture = tmp_path / "fixture"
    workspace = tmp_path / "workspace"
    fixture.mkdir()
    workspace.mkdir()
    for bundle_file in (
        "namelist.def",
        "modpara.def",
        "locspn.def",
        "geometry.dat",
        "coulombintra.def",
        "orbitalidx.def",
        "orbitalidxpara.def",
    ):
        (fixture / bundle_file).write_text(f"{bundle_file}\n")
    (workspace / workspace_file).write_text("generated by vmcdry\n")

    staging_start = block.index("# BEGIN task6a2-authoritative-bundle-staging")
    staging_end = block.index("# END task6a2-authoritative-bundle-staging")
    staging_block = block[staging_start:staging_end]
    harness = tmp_path / "complexuhf_bundle_staging.sh"
    harness.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        'CASE="case_soc_rashba_3d_sub_apbc_xy"\n'
        f'COMPLEXUHF_CASE="{fixture}"\n'
        f'complexuhf_vmcdry_files=("{workspace_file}")\n'
        f"{staging_block}\n"
    )
    result = subprocess.run(
        ["bash", str(harness)],
        capture_output=True,
        text=True,
        cwd=workspace,
    )
    return result, workspace


def test_phase5_run_sh_stages_v37_bundle_between_vmcdry_and_trans():
    block = _complexuhf_preparation_block()
    expected_bundle = (
        "namelist.def modpara.def locspn.def geometry.dat "
        "coulombintra.def orbitalidx.def orbitalidxpara.def"
    )

    staging_start = block.index(
        'case_soc_rashba_3d_sub_apbc_*)\n'
        '        complexuhf_bundle_files=('
    )
    trans_substitution = block.index(
        'cp "${MVMC_WORK}/trans.def" "./trans.def"'
    )
    pre_vmcdry_snapshot = block.index(
        "complexuhf_files_before_vmcdry"
    )
    vmcdry_run = block.index(
        'run_native "${VMCDRY}" stan.in </dev/null > vmcdry.log 2>&1'
    )
    post_vmcdry_snapshot = block.index(
        "complexuhf_files_after_vmcdry"
    )

    assert pre_vmcdry_snapshot < vmcdry_run < post_vmcdry_snapshot
    assert staging_start < trans_substitution
    assert expected_bundle in " ".join(block.split())
    assert (
        'cp "${COMPLEXUHF_CASE}/${bundle_file}" "./${bundle_file}"'
        in block
    )
    assert (
        'cmp -s "${COMPLEXUHF_CASE}/${bundle_file}" "./${bundle_file}"'
        in block
    )


def test_phase5_run_sh_v37_staging_deletes_all_new_vmcdry_files(tmp_path):
    block = _complexuhf_preparation_block()
    expected_allowlist = (
        "namelist.def|modpara.def|locspn.def|geometry.dat|"
        "coulombintra.def|orbitalidx.def|orbitalidxpara.def|"
        "trans.def|initial.def|vmcdry.log)"
    )

    assert 'find . -mindepth 1 -maxdepth 1 -type f -printf \'%f\\0\'' in block
    assert "comm -z -13" in block
    assert block.count(
        'for workspace_file in "${complexuhf_vmcdry_files[@]}"; do'
    ) == 2
    assert 'rm -f -- "./${workspace_file}"' in block
    assert expected_allowlist in block
    assert 'unexpected vmcdry file survived staging' in block
    assert '${COMPLEXUHF_CASE}' in block
    result, workspace = _run_complexuhf_bundle_staging(tmp_path, "lattice.xsf")
    assert result.returncode == 0, result.stderr
    assert not (workspace / "lattice.xsf").exists()
    assert "./*.def" not in block


def test_phase5_run_sh_guards_modpara_override_with_explicit_case_split():
    block = _complexuhf_preparation_block()
    override = 'bash "${COMPLEXUHF_CASE}/complexuhf_modpara_override.txt"'

    assert block.count(override) == 1
    assert (
        'case_soc_rashba_2d_sub_apbc)\n'
        f'        {override}\n'
        '        ;;'
    ) in block
    assert (
        'case_soc_rashba_3d_sub_apbc_*)\n'
        '        # The committed modpara.def is authoritative.\n'
        '        ;;'
    ) in block
    assert '[[ -f "${COMPLEXUHF_CASE}/complexuhf_modpara_override.txt" ]]' not in block
