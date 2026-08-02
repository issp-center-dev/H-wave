"""Tests for tools/_uhfk_to_mvmc/transfer_hermiticity.py.

See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest

from tools._uhfk_to_mvmc.trans_emit import parse_hwave_transfer, TransEmitError
from tools._uhfk_to_mvmc.transfer_hermiticity import (
    check_transfer_dat_hermiticity,
    TransferHermiticityError,
)


_REPO = Path(__file__).resolve().parents[1]
_V36_TRANSFER = str(
    _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "case_soc_rashba_2d_sub_apbc" / "Transfer.dat"
)
_V37_TRANSFER = (
    _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "case_soc_rashba_3d_sub_apbc_xy" / "Transfer.dat"
)
_MANIFEST_PRODUCER = (
    _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "scripts" / "phase2_produce_manifest.py"
)


def test_accepts_v36_shipping_fixture():
    """case_soc_rashba_2d_sub_apbc Transfer.dat MUST pass."""
    assert len(parse_hwave_transfer(_V36_TRANSFER)) == 16
    check_transfer_dat_hermiticity(_V36_TRANSFER)


def test_accepts_v37_shipping_fixture():
    """case_soc_rashba_3d_sub_apbc_xy Transfer.dat MUST pass."""
    assert len(parse_hwave_transfer(_V37_TRANSFER)) == 24
    check_transfer_dat_hermiticity(_V37_TRANSFER)


def test_rejects_header_only_transfer_with_zero_declared_rows(tmp_path):
    """A parsed Transfer.dat with no data rows is not a valid gate input."""
    path = tmp_path / "Transfer.dat"
    path.write_text("Transfer\n1\n0\n")

    with pytest.raises(TransferHermiticityError, match="no data rows"):
        check_transfer_dat_hermiticity(str(path))


def test_rejects_fewer_rows_than_declared_inventory(tmp_path):
    """A complete Hermitian subset cannot satisfy a larger declared inventory."""
    path = tmp_path / "Transfer.dat"
    path.write_text(
        "Transfer\n"
        "1\n"
        "4\n"
        "1 1 1 1\n"
        "     1    0    0    1    2  0.100  0.200\n"
        "    -1    0    0    2    1  0.100 -0.200\n"
    )

    with pytest.raises(TransferHermiticityError, match="declared nr=4.*found 2"):
        check_transfer_dat_hermiticity(str(path))


@pytest.mark.parametrize("coefficient", ["nan", "inf", "-inf"])
@pytest.mark.parametrize("component", ["real", "imaginary"])
def test_rejects_nonfinite_hermitian_pair(tmp_path, coefficient, component):
    """Non-finite coefficients fail before Hermitian partner comparison."""
    path = tmp_path / "Transfer.dat"
    if component == "real":
        first_re, first_im = coefficient, "0.200"
        partner_re, partner_im = coefficient, "-0.200"
    else:
        conjugate = {"nan": "nan", "inf": "-inf", "-inf": "inf"}
        first_re, first_im = "0.100", coefficient
        partner_re, partner_im = "0.100", conjugate[coefficient]
    path.write_text(
        "Transfer\n"
        "1\n"
        "2\n"
        "1 1\n"
        f"     1    0    0    1    2  {first_re} {first_im}\n"
        f"    -1    0    0    2    1  {partner_re} {partner_im}\n"
    )

    with pytest.raises(
        TransferHermiticityError,
        match=r"non-finite coefficient.*row 1",
    ):
        check_transfer_dat_hermiticity(str(path))


def test_rejects_raw_hermitian_rows_made_asymmetric_by_ndegen(tmp_path):
    """Hermiticity is checked after Wannier90 ndegen division."""
    path = tmp_path / "Transfer.dat"
    path.write_text(
        "Transfer\n"
        "1\n"
        "2\n"
        "1 2\n"
        "     1    0    0    1    2  1.000  0.000\n"
        "    -1    0    0    2    1  1.000  0.000\n"
    )

    with pytest.raises(TransferHermiticityError, match="wrong magnitude"):
        check_transfer_dat_hermiticity(str(path))


def test_accepts_post_division_hermitian_rows_with_nonunit_ndegen(tmp_path):
    """Different raw values pass when ndegen division makes them partners."""
    path = tmp_path / "Transfer.dat"
    path.write_text(
        "Transfer\n"
        "1\n"
        "2\n"
        "2 4\n"
        "     1    0    0    1    2  1.000  0.500\n"
        "    -1    0    0    2    1  2.000 -1.000\n"
    )

    check_transfer_dat_hermiticity(str(path))


def test_rejects_missing_partner(tmp_path):
    """Row with dR=(1,0,0), s1=1, s2=2 without matching partner
    (-1,0,0), s1=2, s2=1 with the conjugate value fails."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer\n1\n1\n1\n"
        "     1    0    0    1    2  0.100  0.200\n"
    )
    with pytest.raises(TransferHermiticityError, match="missing partner"):
        check_transfer_dat_hermiticity(str(path))


def test_rejects_missing_real_same_spin_offsite_z_partner(tmp_path):
    """Deleting one real z-hopping partner from the fixture fails."""
    lines = _V37_TRANSFER.read_text().splitlines()
    missing_partner = (
        "     0    0   -1    1    1  -1.000000000000  0.000000000000"
    )
    assert lines.count(missing_partner) == 1
    lines.remove(missing_partner)

    path = tmp_path / "Transfer.dat"
    path.write_text("\n".join(lines) + "\n")

    with pytest.raises(TransferHermiticityError, match="missing partner"):
        check_transfer_dat_hermiticity(str(path))


def test_rejects_wrong_conj_sign(tmp_path):
    """Partner exists but with wrong sign on imag part."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer\n1\n2\n1 1\n"
        "     1    0    0    1    2  0.100  0.200\n"
        "    -1    0    0    2    1  0.100  0.200\n"  # should be -0.200
    )
    with pytest.raises(TransferHermiticityError, match="wrong conj sign"):
        check_transfer_dat_hermiticity(str(path))


def test_rejects_wrong_magnitude(tmp_path):
    """Partner exists with matching sign but different magnitude."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer\n1\n2\n1 1\n"
        "     1    0    0    1    2  0.100  0.200\n"
        "    -1    0    0    2    1  0.500 -0.200\n"  # should be 0.100
    )
    with pytest.raises(TransferHermiticityError, match="wrong magnitude"):
        check_transfer_dat_hermiticity(str(path))


def test_diagonal_rows_pass(tmp_path):
    """Spin-diagonal rows (s1 == s2) with real values pass trivially."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer\n1\n2\n1 1\n"
        "     1    0    0    1    1 -1.000  0.000\n"
        "    -1    0    0    1    1 -1.000  0.000\n"
    )
    check_transfer_dat_hermiticity(str(path))


def test_real_same_spin_onsite_row_passes(tmp_path):
    """A real same-spin R=(0,0,0) row is self-conjugate."""
    path = tmp_path / "Transfer.dat"
    path.write_text(
        "Transfer\n1\n1\n1\n"
        "     0    0    0    1    1  0.100  0.000\n"
    )
    check_transfer_dat_hermiticity(str(path))


def test_manifest_producer_fails_closed_on_nonhermitian_transfer(tmp_path):
    """Producer startup reaches the checker before fixture data loading."""
    case = tmp_path / "case"
    case.mkdir()
    (case / "Transfer.dat").write_text(
        "Transfer\n1\n1\n1\n"
        "     0    0    1    1    1  -1.000  0.000\n"
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        (str(_REPO / "src"), str(_REPO))
    )

    result = subprocess.run(
        [
            sys.executable,
            str(_MANIFEST_PRODUCER),
            "--case", str(case),
            "--out", str(tmp_path / "manifest.json"),
        ],
        cwd=_REPO,
        env=env,
        capture_output=True,
        text=True,
    )

    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "Transfer.dat Hermiticity check failed" in output
    assert "missing partner" in output


def test_rejects_nonunit_ndegen_with_fewer_r_points_than_nr(tmp_path):
    """Non-unit ndegen requires the declared dense R-point listing."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer comment\n"
        "1\n"
        "7\n"
        "1 1 1 1 1 1 2\n"
        "     1    0    0    1    2  0.100  0.200\n"
        "    -1    0    0    2    1  0.100 -0.200\n"
    )
    with pytest.raises(TransferHermiticityError, match="dense file listing"):
        check_transfer_dat_hermiticity(str(path))


def test_preserves_canonical_parser_error_message(tmp_path):
    """Canonical parser failures are retyped without changing the message."""
    path = tmp_path / "trans.dat"
    path.write_text("")

    with pytest.raises(TransEmitError) as canonical_error:
        parse_hwave_transfer(str(path))
    with pytest.raises(TransferHermiticityError) as checker_error:
        check_transfer_dat_hermiticity(str(path))

    assert str(checker_error.value) == str(canonical_error.value)


def test_missing_file_still_raises_file_not_found(tmp_path):
    """A missing path keeps the checker's existing exception contract."""
    path = tmp_path / "missing" / "Transfer.dat"

    with pytest.raises(FileNotFoundError) as error:
        check_transfer_dat_hermiticity(str(path))

    assert str(error.value) == f"Transfer.dat not found: {path}"


def test_rejects_onsite_same_spin_with_nonzero_imag(tmp_path):
    """A same-spin on-site row with im != 0 violates the invariant
    H[s,s](0) = conj(H[s,s](0)) = H[s,s](0)* which implies im must be
    0. Currently the trivial-skip branch also catches this via the
    partner-key mechanism (self-referential lookup fails the sign
    check), but this test pins the intended rejection."""
    path = tmp_path / "trans.dat"
    path.write_text(
        "Transfer comment\n1\n1\n1\n"
        "     0    0    0    1    1  0.100  0.500\n"
    )
    with pytest.raises(TransferHermiticityError):
        check_transfer_dat_hermiticity(str(path))
