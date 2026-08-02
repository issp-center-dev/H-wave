"""Allowlist coverage checker tests.

Pins the static checker that iterates every case_*/input.toml
under a fixture root, feeds each through the CLI's own predicate
(no drift), and enforces expected-inventory contract (all required
present, no duplicates on allowed triples).
"""
from __future__ import annotations

from pathlib import Path

import pytest

from tools._uhfk_to_mvmc.allowlist_coverage import (
    check_fixture_allowlist_coverage,
    AllowlistCoverageError,
)


_REPO = Path(__file__).resolve().parents[1]


def test_allowlist_accepts_existing_non_soc_fixture_via_cli_predicate(
    tmp_path,
):
    """Non-SOC fixture (single-orbital 1D Hubbard style) MUST pass via
    the CLI predicate's early-return branch. The checker must reuse
    is_supported_triple so non-SOC /
    SubShape=[1,1,1] fixtures are admissible without listing them in
    expected_map."""
    fixture = tmp_path / "case_stub_non_soc"
    fixture.mkdir()
    (fixture / "input.toml").write_text(
        """
[mode]
mode = "UHFk"
[mode.param]
Ncond = 4
CellShape = [8, 1, 1]
SubShape = [1, 1, 1]
BoundaryCondition = ["antiperiodic", "periodic", "periodic"]
"""
    )
    check_fixture_allowlist_coverage(
        tmp_path, expected_map={}, require_all_expected=False,
    )


def test_allowlist_rejects_unlisted_soc_apbc_subshape_shape(tmp_path):
    """An unlisted SOC + APBC + SubShape > 1 shape fails the CLI predicate."""
    fixture = tmp_path / "case_stub_novel_shape"
    fixture.mkdir()
    (fixture / "input.toml").write_text(
        """
[mode]
mode = "UHFk"
enable_spin_orbital = true
[mode.param]
Ncond = 4
CellShape = [8, 8, 8]
SubShape = [2, 2, 2]
BoundaryCondition = ["antiperiodic", "antiperiodic", "periodic"]
"""
    )
    with pytest.raises(AllowlistCoverageError,
                       match="not in the supported allowlist"):
        check_fixture_allowlist_coverage(
            tmp_path, expected_map={}, require_all_expected=False,
        )


def test_allowlist_rejects_missing_expected_fixture(tmp_path):
    """require_all_expected=True and the expected fixture is missing
    from the discovered set → fail with 'missing required fixture'."""
    expected_map = {
        "case_stub_missing/case_soc_rashba_3d_sub_apbc_xy":
            ((1, 1, 0), (2, 2, 2), (4, 4, 4)),
    }
    with pytest.raises(AllowlistCoverageError,
                       match="missing required fixture"):
        check_fixture_allowlist_coverage(
            tmp_path, expected_map=expected_map,
            require_all_expected=True,
        )


def test_allowlist_rejects_duplicate_triple(tmp_path):
    """Two fixture directories both claiming the SAME allowed triple
    are rejected as duplicates (unless both entries are in expected_map).
    """
    for name in ("case_stub_dup_a", "case_stub_dup_b"):
        (tmp_path / name).mkdir()
        (tmp_path / name / "input.toml").write_text(
            """
[mode]
mode = "UHFk"
enable_spin_orbital = true
[mode.param]
Ncond = 4
CellShape = [4, 4, 4]
SubShape = [2, 2, 2]
BoundaryCondition = ["antiperiodic", "antiperiodic", "periodic"]
"""
        )
    with pytest.raises(AllowlistCoverageError, match="duplicate"):
        check_fixture_allowlist_coverage(
            tmp_path, expected_map={}, require_all_expected=False,
        )


def test_allowlist_accepts_mirror_pair_across_two_family_roots(tmp_path):
    """Mirror-fixture convention: same basename `case_soc_stub_xy` under
    two different parent dirs (mimicking uhfk_mvmc_pairproduct/ +
    apbc_complexuhf/ layout) with the same SOC+APBC+SubShape triple MUST
    pass, since each shipping fixture is present under BOTH roots by
    design."""
    for root in ("uhfk_mvmc_pairproduct", "apbc_complexuhf"):
        (tmp_path / root / "case_soc_stub_xy").mkdir(parents=True)
        (tmp_path / root / "case_soc_stub_xy" / "input.toml").write_text(
            """
[mode]
mode = "UHFk"
enable_spin_orbital = true
[mode.param]
Ncond = 4
CellShape = [4, 4, 4]
SubShape = [2, 2, 2]
BoundaryCondition = ["antiperiodic", "antiperiodic", "periodic"]
"""
        )
    check_fixture_allowlist_coverage(
        tmp_path, expected_map={}, require_all_expected=False,
    )


def test_allowlist_rejects_three_mirror_basenames(tmp_path):
    """Bounded mirror tolerance: MORE than 2 same-basename fixtures
    sharing a gated triple should still be rejected -- mirror
    convention is strictly a 2-fixture pattern."""
    for root in ("uhfk_mvmc_pairproduct", "apbc_complexuhf", "third_root"):
        (tmp_path / root / "case_soc_stub_xy").mkdir(parents=True)
        (tmp_path / root / "case_soc_stub_xy" / "input.toml").write_text(
            """
[mode]
mode = "UHFk"
enable_spin_orbital = true
[mode.param]
Ncond = 4
CellShape = [4, 4, 4]
SubShape = [2, 2, 2]
BoundaryCondition = ["antiperiodic", "antiperiodic", "periodic"]
"""
        )
    with pytest.raises(AllowlistCoverageError, match="duplicate"):
        check_fixture_allowlist_coverage(
            tmp_path, expected_map={}, require_all_expected=False,
        )


def test_allowlist_covers_real_shipping_fixtures():
    """Iterating tests/validation/ MUST accept every fixture on disk."""
    fixture_root = _REPO / "tests" / "validation"
    check_fixture_allowlist_coverage(
        fixture_root, expected_map={}, require_all_expected=False,
    )
