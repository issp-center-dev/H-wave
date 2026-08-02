"""Regression tests for non-vacuous G2 cross-solver density gates."""
from __future__ import annotations

import importlib.util
import re
from pathlib import Path

import numpy as np
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
VALIDATION_ROOT = (
    REPO_ROOT / "tests" / "validation" / "uhfk_mvmc_pairproduct"
)


def _load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_minimal_g2b_workspace(
    root: Path,
    reference: np.ndarray,
    *,
    perturb_scale: float | None = None,
    return_seed_unchanged: bool,
    flag_fock: bool = True,
) -> None:
    ns = reference.shape[0] // 2
    hwave = root / "hwave"
    complexuhf = root / "complexuhf"
    hwave.mkdir(parents=True)
    complexuhf.mkdir(parents=True)

    # Every shipping fixture that runs G2 shares ComplexUHF's Fock
    # functional, and the contraction precondition is unconditional.
    # The false variant is retained only for the adversarial regression
    # proving that input metadata cannot disable that precondition.
    (hwave / "input.toml").write_text(
        "[mode]\nmode = \"UHFk\"\nenable_spin_orbital = true\n"
        f"flag_fock = {str(flag_fock).lower()}\n"
        "[mode.param]\nNcond = 2\nCellShape = [2, 1, 1]\n"
        "SubShape = [1, 1, 1]\n"
        "BoundaryCondition = [\"periodic\", \"periodic\", \"periodic\"]\n"
        "T = 0.0\n",
        encoding="utf-8",
    )
    (hwave / "geometry_uhf.dat").write_text(
        "1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n"
        "0.0 0.0 0.0\n2 0 0\n0 1 0\n0 0 1\n"
        "0 0 0 0\n1 0 0 0\n",
        encoding="utf-8",
    )
    np.savez(
        hwave / "green.npz",
        green_sublattice=np.zeros((2, 1, 4, 1, 4), dtype=np.complex128),
    )

    if perturb_scale is None:
        # Track the single shipped perturbation literal so this test cannot
        # drift from the value used by run.sh.
        run_text = (VALIDATION_ROOT / "run.sh").read_text(encoding="utf-8")
        match = re.search(
            r'^\s*PERTURB_SCALE="([^"]+)"\s*$',
            run_text,
            flags=re.MULTILINE,
        )
        assert match is not None, "run.sh no longer pins one G2 perturb-scale"
        perturb_scale = float(match.group(1))

    seed_mod = _load_module(
        "_seed_complexuhf_for_g2_test",
        VALIDATION_ROOT / "scripts" / "seed_complexuhf_from_hwave.py",
    )
    initial_path = complexuhf / "initial.def"
    seed_mod._write_initial_def(
        reference,
        initial_path,
        perturb_scale=perturb_scale,
        perturb_seed=42,
    )

    seeded_values = {}
    for line in initial_path.read_text(encoding="utf-8").splitlines()[5:]:
        i, s, j, t, real, imag = line.split()
        seeded_values[(int(i), int(s), int(j), int(t))] = (
            float(real),
            float(imag),
        )
    with (complexuhf / "zvo_UHF_cisajs.dat").open("w") as fp:
        for i in range(ns):
            for s in (0, 1):
                for j in range(ns):
                    for t in (0, 1):
                        if return_seed_unchanged:
                            real, imag = seeded_values.get(
                                (i, s, j, t), (0.0, 0.0)
                            )
                        else:
                            value = reference[i + ns * s, j + ns * t]
                            real, imag = float(value.real), float(value.imag)
                        fp.write(f"{i} {s} {j} {t} {real:.16e} {imag:.16e}\n")


def _run_g2b(compare_mod, root, monkeypatch, reference, *, tol="1e-6"):
    from tools._uhfk_to_mvmc import density_check

    ns = reference.shape[0] // 2
    original = density_check.gauge_lift

    def _reference_gauge_lift(_green, i, s, j, t, **_kwargs):
        return reference[i + ns * s, j + ns * t]

    _reference_gauge_lift.__module__ = original.__module__
    _reference_gauge_lift.__qualname__ = original.__qualname__
    monkeypatch.setattr(density_check, "gauge_lift", _reference_gauge_lift)

    return compare_mod.main(
        [
            "compare.py",
            "--workspace",
            str(root),
            "--mode",
            "g2b",
            f"--gtol_gauge={tol}",
        ]
    )


def _replace_real_token(path: Path, line_index: int, token: str) -> None:
    lines = path.read_text(encoding="utf-8").splitlines()
    fields = lines[line_index].split()
    fields[4] = token
    lines[line_index] = " ".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_shipped_seed_readback_fails_g2b(tmp_path, monkeypatch, capsys):
    """Returning the shipped seed unchanged must never satisfy G2b."""
    compare_mod = _load_module(
        "_compare_for_g2_contraction_test",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path, reference, return_seed_unchanged=True
    )

    from tools._uhfk_to_mvmc import density_check

    original = density_check.gauge_lift

    def _reference_gauge_lift(_green, i, s, j, t, **_kwargs):
        return reference[i + 2 * s, j + 2 * t]

    _reference_gauge_lift.__module__ = original.__module__
    _reference_gauge_lift.__qualname__ = original.__qualname__
    monkeypatch.setattr(density_check, "gauge_lift", _reference_gauge_lift)

    exit_code = compare_mod.main(
        [
            "compare.py",
            "--workspace",
            str(tmp_path),
            "--mode",
            "g2b",
            "--gtol_gauge",
            "1e-6",
        ]
    )

    assert exit_code != 0, capsys.readouterr()


def test_g2b_rejects_converged_output_when_seed_started_inside_margin(
    tmp_path, monkeypatch, capsys
):
    """A good final density cannot rescue an initially vacuous G2 setup."""
    compare_mod = _load_module(
        "_compare_for_g2_inside_margin_test",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=1e-6,
        return_seed_unchanged=False,
    )

    from tools._uhfk_to_mvmc import density_check

    original = density_check.gauge_lift

    def _reference_gauge_lift(_green, i, s, j, t, **_kwargs):
        return reference[i + 2 * s, j + 2 * t]

    _reference_gauge_lift.__module__ = original.__module__
    _reference_gauge_lift.__qualname__ = original.__qualname__
    monkeypatch.setattr(density_check, "gauge_lift", _reference_gauge_lift)

    exit_code = compare_mod.main(
        [
            "compare.py",
            "--workspace",
            str(tmp_path),
            "--mode",
            "g2b",
            "--gtol_gauge",
            "1e-6",
        ]
    )

    assert exit_code != 0, capsys.readouterr()


def test_g2b_pass_record_reports_contraction_evidence(
    tmp_path, monkeypatch, capsys
):
    compare_mod = _load_module(
        "_compare_for_g2_evidence_test",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=1e-3,
        return_seed_unchanged=False,
    )

    from tools._uhfk_to_mvmc import density_check

    original = density_check.gauge_lift

    def _reference_gauge_lift(_green, i, s, j, t, **_kwargs):
        return reference[i + 2 * s, j + 2 * t]

    _reference_gauge_lift.__module__ = original.__module__
    _reference_gauge_lift.__qualname__ = original.__qualname__
    monkeypatch.setattr(density_check, "gauge_lift", _reference_gauge_lift)

    exit_code = compare_mod.main(
        [
            "compare.py",
            "--workspace",
            str(tmp_path),
            "--mode",
            "g2b",
            "--gtol_gauge",
            "1e-6",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 0, captured
    assert re.search(
        r"^G2b PASS mode=g2b .* initial_delta=[0-9.e+-]+ "
        r"final_delta=[0-9.e+-]+ contraction_ratio=[0-9.e+-]+$",
        captured.out,
    )


def test_flag_fock_false_cannot_skip_g2b_contraction_precondition(
    tmp_path, monkeypatch, capsys
):
    compare_mod = _load_module(
        "_compare_for_g2_no_exemption_test",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=0.0,
        return_seed_unchanged=False,
        flag_fock=False,
    )

    exit_code = _run_g2b(compare_mod, tmp_path, monkeypatch, reference)

    captured = capsys.readouterr()
    assert exit_code != 0, captured
    assert " PASS " not in captured.out
    assert not hasattr(compare_mod, "_fixture_is_fock_consistent")


@pytest.mark.parametrize(
    "token,value",
    [("nan", np.nan), ("inf", np.inf), ("-inf", -np.inf)],
)
def test_g2b_rejects_non_finite_reference_density(
    tmp_path, monkeypatch, capsys, token, value
):
    compare_mod = _load_module(
        f"_compare_for_g2_nonfinite_reference_{token}",
        VALIDATION_ROOT / "compare.py",
    )
    clean_reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        clean_reference,
        perturb_scale=1e-3,
        return_seed_unchanged=False,
    )
    bad_reference = clean_reference.copy()
    bad_reference[0, 0] = value

    exit_code = _run_g2b(compare_mod, tmp_path, monkeypatch, bad_reference)

    captured = capsys.readouterr()
    assert exit_code != 0, captured
    assert " PASS " not in captured.out


@pytest.mark.parametrize("token", ["nan", "inf", "-inf"])
def test_g2b_rejects_non_finite_final_density(
    tmp_path, monkeypatch, capsys, token
):
    compare_mod = _load_module(
        f"_compare_for_g2_nonfinite_final_{token}",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=1e-3,
        return_seed_unchanged=False,
    )
    _replace_real_token(
        tmp_path / "complexuhf" / "zvo_UHF_cisajs.dat", 0, token
    )

    with pytest.raises(compare_mod.ComplexUHFParseError, match="non-finite"):
        _run_g2b(compare_mod, tmp_path, monkeypatch, reference)

    assert " PASS " not in capsys.readouterr().out


@pytest.mark.parametrize("token", ["nan", "inf", "-inf"])
def test_g2b_rejects_non_finite_initial_density(
    tmp_path, monkeypatch, capsys, token
):
    compare_mod = _load_module(
        f"_compare_for_g2_nonfinite_initial_{token}",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=1e-3,
        return_seed_unchanged=False,
    )
    _replace_real_token(tmp_path / "complexuhf" / "initial.def", 5, token)

    with pytest.raises(
        compare_mod.ComplexUHFInitialParseError, match="non-finite"
    ):
        _run_g2b(compare_mod, tmp_path, monkeypatch, reference)

    assert " PASS " not in capsys.readouterr().out


@pytest.mark.parametrize("tol", ["nan", "inf", "-inf"])
def test_g2b_rejects_non_finite_tolerance(
    tmp_path, monkeypatch, capsys, tol
):
    compare_mod = _load_module(
        f"_compare_for_g2_nonfinite_tol_{tol}",
        VALIDATION_ROOT / "compare.py",
    )
    reference = 0.25 * np.eye(4, dtype=np.complex128)
    _write_minimal_g2b_workspace(
        tmp_path,
        reference,
        perturb_scale=1e-3,
        return_seed_unchanged=False,
    )

    exit_code = _run_g2b(
        compare_mod, tmp_path, monkeypatch, reference, tol=tol
    )

    captured = capsys.readouterr()
    assert exit_code != 0, captured
    assert " PASS " not in captured.out
