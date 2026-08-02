"""Tests for the InOrbitalGeneral writer and aggregator."""
from __future__ import annotations

import json
import os
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.general_output_writer import (
    ClassInconsistencyError,
    aggregate_general_orbital_params,
    write_zqp_orbital_general,
)
from tools._uhfk_to_mvmc.density_check import (
    DensityMismatchError,
    compare_against_onebodyg_uhf_general,
)
from tools._uhfk_to_mvmc.pair_product_density import parse_emitted_F


def _antisym_F(nsite):
    """Return a random-but-reproducible antisymmetric (2Ns, 2Ns) F."""
    rng = np.random.default_rng(1234)
    n = 2 * nsite
    R = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (R - R.T) / 2  # antisymmetric


def test_aggregate_general_orbital_params_consistent_class_averages_signed_value():
    """Two rows mapped to the same idx with consistent signed values →
    params[idx] = the shared value (average of equal numbers)."""
    F = _antisym_F(nsite=2)
    # Ns = 2, so 2Ns = 4, upper triangle has 2*Ns*(2Ns-1)/2 = 6 rows.
    # Put (all_i=0, all_j=1) and (all_i=2, all_j=3) into the same class 0,
    # both with sign=+1. F[0,1] and F[2,3] must be equal for the class-
    # consistency check to pass; force them equal here.
    F[0, 1] = 0.42 + 0.13j
    F[1, 0] = -F[0, 1]
    F[2, 3] = 0.42 + 0.13j
    F[3, 2] = -F[2, 3]
    mapping = {
        (0, 1): (0, 1),
        (2, 3): (0, 1),
        (0, 2): (1, 1),
        (0, 3): (2, 1),
        (1, 2): (3, 1),
        (1, 3): (4, 1),
    }
    params = aggregate_general_orbital_params(
        F, mapping=mapping, n_orbital_idx=5,
        epsilon_noise=0.0,
    )
    assert params[0] == pytest.approx(0.42 + 0.13j, abs=1e-14)


def test_aggregate_general_raises_on_inconsistent_class_entries():
    """Two rows mapped to the same idx with unequal signed values →
    ClassInconsistencyError including idx and max residual."""
    F = _antisym_F(nsite=2)
    F[0, 1] = 0.5
    F[1, 0] = -F[0, 1]
    F[2, 3] = 0.5 + 0.1j  # differs from F[0, 1] by 0.1j
    F[3, 2] = -F[2, 3]
    mapping = {
        (0, 1): (0, 1),
        (2, 3): (0, 1),
        (0, 2): (1, 1),
        (0, 3): (2, 1),
        (1, 2): (3, 1),
        (1, 3): (4, 1),
    }
    with pytest.raises(ClassInconsistencyError) as excinfo:
        aggregate_general_orbital_params(
            F, mapping=mapping, n_orbital_idx=5,
            epsilon_noise=0.0,
        )
    msg = str(excinfo.value)
    assert "idx" in msg and "0" in msg
    assert "residual" in msg


def test_write_zqp_orbital_general_roundtrip_matches_input_params():
    """Round-trip: write params to file then re-parse; expect exact match
    on real/imag parts."""
    params = np.array([0.1 + 0.0j, -0.2 + 0.5j, 0.0 + 0.0j], dtype=np.complex128)
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as tmp:
        write_zqp_orbital_general(tmp.name, params)
        with open(tmp.name) as fp:
            body = fp.read()
    # Header assertions
    assert "NOrbitalIdx" in body
    # Row body: "<idx> <real> <imag>" per line
    data_lines = [
        ln for ln in body.splitlines()
        if ln.strip() and not ln.startswith("==") and not ln.startswith("N")
        and not ln.startswith("i") and not ln.startswith("=")
    ]
    parsed = []
    for ln in data_lines:
        toks = ln.strip().split()
        if len(toks) == 3:
            parsed.append(complex(float(toks[1]), float(toks[2])))
    assert len(parsed) == len(params)
    for p_got, p_want in zip(parsed, params):
        assert p_got == pytest.approx(p_want, abs=1e-15)


def test_compare_against_onebodyg_uhf_general_accepts_matching_g():
    """Synthetic G matching greenone.dat entries → no exception."""
    G_all = np.zeros((4, 4), dtype=np.complex128)
    G_all[0, 1] = 0.25
    G_all[2, 3] = 0.10
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as tmp:
        # greenone.dat format: i s j t re im
        tmp.write("0 0 1 0  0.25 0.0\n")
        tmp.write("0 1 1 1  0.10 0.0\n")  # spin-down entry: G[0+Ns, 1+Ns]=G[2, 3]=0.10
        tmp.write("0 0 0 1  0.0 0.0\n")  # spin off-diagonal: expected 0
        tmp.name
        path = tmp.name
    compare_against_onebodyg_uhf_general(G_all, path, tol=1e-10)


def test_compare_against_onebodyg_uhf_general_flags_mismatch():
    G_all = np.zeros((4, 4), dtype=np.complex128)
    G_all[0, 1] = 0.25
    with tempfile.NamedTemporaryFile("w", suffix=".dat", delete=False) as tmp:
        tmp.write("0 0 1 0  0.10 0.0\n")  # want 0.10 but bridge gives 0.25
        path = tmp.name
    with pytest.raises(DensityMismatchError, match="differ"):
        compare_against_onebodyg_uhf_general(G_all, path, tol=1e-10)


def _write_golden_bridge_workspace(workspace_dir, golden):
    """Write the minimal namelist.def, modpara.def, orbitalidx_general.def and
    zqp_orbital_uhfk.dat files that ``parse_emitted_F`` consumes.

    Uses the mVMC namelist -> ModPara -> Nsite chain. Only
    the fields exercised by parse_emitted_F are populated; the other
    namelist entries mirror the real bridge output for realism but are not
    read by the parser.
    """
    Nsite = golden["Nsite"]
    orbitalidx_rows = golden["orbitalidx"]
    zqp_rows = golden["zqp"]

    namelist_path = os.path.join(workspace_dir, "namelist.def")
    modpara_path = os.path.join(workspace_dir, "modpara.def")
    orbitalidx_path = os.path.join(
        workspace_dir, "orbitalidx_general.def"
    )
    zqp_path = os.path.join(workspace_dir, "zqp_orbital_uhfk.dat")

    with open(namelist_path, "w") as fp:
        fp.write(f"         ModPara  {os.path.basename(modpara_path)}\n")
        fp.write("  OrbitalGeneral  orbitalidx_general.def\n")
        fp.write(
            f"InOrbitalGeneral  {os.path.basename(zqp_path)}\n"
        )
    with open(modpara_path, "w") as fp:
        fp.write("--------------------\n")
        fp.write("Model_Parameters   0\n")
        fp.write("--------------------\n")
        fp.write(f"Nsite          {Nsite}\n")
    # orbitalidx_general.def: 5 header lines + 15 mapping rows for Ns=3 +
    # NOrbitalIdx optimize-flag rows.
    n_orbital_idx = max(int(r["class_idx"]) for r in orbitalidx_rows) + 1
    with open(orbitalidx_path, "w") as fp:
        fp.write("=============================================\n")
        fp.write(f"NOrbitalIdx {n_orbital_idx:10d}\n")
        fp.write(f"ComplexType {1:10d}\n")
        fp.write("=============================================\n")
        fp.write("=============================================\n")
        for r in orbitalidx_rows:
            fp.write(
                f"    {int(r['i']):3d}  {int(r['spn_i']):d}"
                f"     {int(r['j']):3d}  {int(r['spn_j']):d}"
                f"     {int(r['class_idx']):6d}  {int(r['sign']):2d}\n"
            )
        for k in range(n_orbital_idx):
            fp.write(f"    {k:3d}      1\n")
    with open(zqp_path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx  {n_orbital_idx}\n")
        fp.write("======================\n")
        fp.write("== i_j_OrbitalIdx ===\n")
        fp.write("======================\n")
        # Sort by class_idx to write in ascending index order.
        for r in sorted(zqp_rows, key=lambda x: int(x["class_idx"])):
            fp.write(
                "{:d} {: .18e} {: .18e}\n".format(
                    int(r["class_idx"]), float(r["re"]), float(r["im"])
                )
            )


def test_parse_emitted_F_matches_golden(tmp_path):
    """Load ``tests/data/v36_writer_check_golden.json``, materialize a minimal
    but real bridge workspace, call ``parse_emitted_F`` and assert the
    reconstructed F matches the golden expected F at atol=1e-15.

    The adversarial fixture has Nsite=3, non-sequential class ids,
    sign=-1 row on a cross-spin upper-triangle pair, spin-block-
    disambiguating row (i=1, spn_i=0, j=0, spn_j=1) -> all_i=1, all_j=3.
    A parser using site-major ``all = spn + i * 2`` would place the
    disambiguating row on ``(all_i=2, all_j=0)`` (lower triangle) and either
    fail the antisymmetry pin or mis-locate the value; the golden expected
    F traps that.
    """
    golden_path = os.path.join(
        os.path.dirname(__file__), "data", "v36_writer_check_golden.json"
    )
    with open(golden_path) as fp:
        golden = json.load(fp)

    workspace = tmp_path / "bridge"
    workspace.mkdir()
    _write_golden_bridge_workspace(str(workspace), golden)

    F = parse_emitted_F(str(workspace))
    F_expected = np.array(golden["expected_F"]["real"], dtype=np.float64) + \
        1j * np.array(golden["expected_F"]["imag"], dtype=np.float64)
    assert F.shape == F_expected.shape == (
        2 * golden["Nsite"], 2 * golden["Nsite"]
    )
    delta = float(np.max(np.abs(F - F_expected)))
    assert delta < 1e-15, (
        f"parse_emitted_F does not reproduce golden F at 1e-15; "
        f"max abs delta = {delta}\n"
        f"F[0,4] = {F[0,4]}, expected {F_expected[0,4]}\n"
        f"F[1,3] = {F[1,3]}, expected {F_expected[1,3]}\n"
    )
