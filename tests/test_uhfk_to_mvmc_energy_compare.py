"""Pinning tests for energy_relative_delta and the G3 workspace resolver.

See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst`` for the validation
contract.
"""
from __future__ import annotations

import pytest

from tools._uhfk_to_mvmc.energy_compare import (
    EnergyCompareError,
    _parse_mvmc_zvo_out,
    _resolve_g3_paths,
    energy_relative_delta,
)


def _write_hwave_energy(tmp_path, energy_total: float, extra_lines=None):
    path = tmp_path / "energy.dat"
    lines = [f"Energy_Total = {energy_total:.15e}\n"]
    if extra_lines:
        lines.extend(extra_lines)
    path.write_text("".join(lines))
    return str(path)


def _write_mvmc_zvo_out(tmp_path, samples):
    path = tmp_path / "zvo_out_001.dat"
    lines = []
    for h in samples:
        lines.append(f"{h:.15e} {h*h:.15e} 0.0 0.0 0.0 0.0\n")
    path.write_text("".join(lines))
    return str(path)


def test_energy_relative_delta_hand(tmp_path):
    """Hand-written E_h=1.0, E_m=1.001, delta_rel=1e-3."""
    hwave = _write_hwave_energy(tmp_path, 1.0)
    mvmc = _write_mvmc_zvo_out(tmp_path, [1.001])
    e_h, e_m, delta = energy_relative_delta(hwave, mvmc)
    assert e_h == pytest.approx(1.0, abs=1e-15)
    assert e_m == pytest.approx(1.001, abs=1e-15)
    assert delta == pytest.approx(1e-3, abs=1e-15)


def test_energy_relative_delta_hand_negative_hwave(tmp_path):
    """Verify the |E_hwave| in the denominator handles negatives."""
    hwave = _write_hwave_energy(tmp_path, -3.75)
    mvmc = _write_mvmc_zvo_out(tmp_path, [-3.74625])
    e_h, e_m, delta = energy_relative_delta(hwave, mvmc)
    assert delta == pytest.approx(1e-3, abs=1e-15)


def test_energy_relative_delta_averages_mvmc_bins(tmp_path):
    """mVMC returns the mean of the per-bin <H> samples (column 0)."""
    hwave = _write_hwave_energy(tmp_path, -2.0)
    mvmc = _write_mvmc_zvo_out(tmp_path, [-2.0, -2.0, -2.0])
    _, e_m, delta = energy_relative_delta(hwave, mvmc)
    assert e_m == pytest.approx(-2.0, abs=1e-15)
    assert delta == 0.0

    sub = tmp_path / "sub"
    sub.mkdir()
    hwave2 = _write_hwave_energy(sub, -2.0)
    mvmc2 = _write_mvmc_zvo_out(sub, [-1.999, -2.001, -2.003])
    _, e_m2, _ = energy_relative_delta(hwave2, mvmc2)
    assert e_m2 == pytest.approx(-2.001, abs=1e-14)


def test_energy_relative_delta_missing_hwave_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="hwave energy file not found"):
        energy_relative_delta(
            str(tmp_path / "does_not_exist.dat"),
            _write_mvmc_zvo_out(tmp_path, [1.0]),
        )


def test_energy_relative_delta_missing_mvmc_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="mvmc zvo_out file not found"):
        energy_relative_delta(
            _write_hwave_energy(tmp_path, 1.0),
            str(tmp_path / "does_not_exist.dat"),
        )


def test_energy_relative_delta_hwave_no_energy_total_raises(tmp_path):
    path = tmp_path / "energy.dat"
    path.write_text("Energy_Band = -1.0\n")
    with pytest.raises(EnergyCompareError, match="Energy_Total not found"):
        energy_relative_delta(
            str(path), _write_mvmc_zvo_out(tmp_path, [1.0])
        )


def test_energy_relative_delta_mvmc_empty_raises(tmp_path):
    path = tmp_path / "zvo_out_001.dat"
    path.write_text("# only a comment\n\n")
    with pytest.raises(EnergyCompareError, match="no numeric rows"):
        energy_relative_delta(
            _write_hwave_energy(tmp_path, 1.0), str(path)
        )


def test_parse_mvmc_zvo_out_rejects_truncated_row(tmp_path):
    path = tmp_path / "zvo_out_001.dat"
    path.write_text("-1.25\n")
    with pytest.raises(EnergyCompareError, match="expected exactly 6 columns"):
        _parse_mvmc_zvo_out(str(path))


def test_parse_mvmc_zvo_out_rejects_non_finite_non_energy_column(tmp_path):
    path = tmp_path / "zvo_out_001.dat"
    path.write_text("-1.25 nan 0.0 0.0 0.0 0.0\n")
    with pytest.raises(EnergyCompareError, match="non-finite"):
        _parse_mvmc_zvo_out(str(path))


def test_parse_mvmc_zvo_out_rejects_malformed_row(tmp_path):
    path = tmp_path / "zvo_out_001.dat"
    path.write_text("-1.25 invalid 0.0 0.0 0.0 0.0\n")
    with pytest.raises(EnergyCompareError, match="numeric"):
        _parse_mvmc_zvo_out(str(path))


@pytest.mark.skip(
    reason="A regenerated case_pbc energy reference is not available"
)
def test_energy_relative_delta_v1_case_pbc():
    """Deferred case_pbc energy comparison.

    On the case_pbc fixture output, assert delta_rel matches the existing
    compare.py energy report within 5e-7. This requires a regenerated
    fixture reference.
    """
    pass


def test_resolve_g3_paths_returns_canonical_paths(tmp_path):
    """Happy path: both files present."""
    hwave_dir = tmp_path / "hwave"
    hwave_dir.mkdir()
    (hwave_dir / "energy.dat").write_text("Energy_Total = -1.0\n")
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    (mvmc_dir / "zvo_out_selected.dat").write_text("-1.0 1.0 0 0 0 0\n")
    h_path, m_path = _resolve_g3_paths(str(tmp_path))
    assert h_path == str(hwave_dir / "energy.dat")
    assert m_path == str(mvmc_dir / "zvo_out_selected.dat")


def test_resolve_g3_paths_missing_hwave_energy(tmp_path):
    """Missing hwave/energy.dat raises with the documented message."""
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    (mvmc_dir / "zvo_out_selected.dat").write_text("-1.0 1.0 0 0 0 0\n")
    with pytest.raises(FileNotFoundError, match="G3: missing .*energy.dat"):
        _resolve_g3_paths(str(tmp_path))


def test_resolve_g3_paths_missing_zvo_out_selected(tmp_path):
    """Missing mvmc/zvo_out_selected.dat raises with the
    documented run.sh normalization message. This test explicitly does NOT
    scan raw zvo_out_*.dat files — that is `normalize_mvmc_output`'s job in
    run.sh."""
    hwave_dir = tmp_path / "hwave"
    hwave_dir.mkdir()
    (hwave_dir / "energy.dat").write_text("Energy_Total = -1.0\n")
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    # Include a raw zvo_out_001.dat — resolver MUST NOT pick it up.
    (mvmc_dir / "zvo_out_001.dat").write_text("-1.0 1.0 0 0 0 0\n")
    with pytest.raises(
        FileNotFoundError, match="G3: missing.*zvo_out_selected.dat.*normalize"
    ):
        _resolve_g3_paths(str(tmp_path))
