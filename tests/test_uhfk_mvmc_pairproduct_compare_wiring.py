"""Wiring tests for the seven-gate dispatcher in
``tests/validation/uhfk_mvmc_pairproduct/compare.py``.

The compare.py dispatcher MUST enforce:
  - EXPECTED_MODE_DISPATCH is a deep-frozen reference table.
  - MODE_DISPATCH is a mutable copy compared structurally against
    EXPECTED_MODE_DISPATCH before every mode.
  - Every canonical helper resolves via importlib + getattr, and the
    resolved object's __module__ / __qualname__ MUST match the table.
  - Any integrity failure exits code 2 BEFORE any PASS line.

These tests target the dispatch framework itself. See
``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst`` for the validation contract.
"""
from __future__ import annotations

import importlib
import os
import sys
import types

import numpy as np
import pytest


# The compare.py under test is not a package member; load it by path.
_COMPARE_PATH = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "tests",
        "validation",
        "uhfk_mvmc_pairproduct",
        "compare.py",
    )
)
if not os.path.isfile(_COMPARE_PATH):
    _COMPARE_PATH = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__),
            "validation",
            "uhfk_mvmc_pairproduct",
            "compare.py",
        )
    )


def _load_compare_module():
    """Import compare.py as a fresh module every call so tests can
    mutate MODE_DISPATCH without leaking to other tests."""
    spec = importlib.util.spec_from_file_location(
        "_v36_compare_under_test", _COMPARE_PATH,
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def compare_mod():
    return _load_compare_module()


# ---------------------------------------------------------------------
# Deep immutability of EXPECTED_MODE_DISPATCH.
# ---------------------------------------------------------------------


def test_expected_mode_dispatch_top_level_frozen(compare_mod):
    with pytest.raises(TypeError):
        compare_mod.EXPECTED_MODE_DISPATCH["g1"] = "poisoned"  # type: ignore


def test_expected_mode_dispatch_nested_helpers_frozen(compare_mod):
    with pytest.raises(TypeError):
        compare_mod.EXPECTED_MODE_DISPATCH["g1"]["helpers"] = ()  # type: ignore


def test_expected_mode_dispatch_helpers_tuple_frozen(compare_mod):
    helpers = compare_mod.EXPECTED_MODE_DISPATCH["g1"]["helpers"]
    with pytest.raises(TypeError):
        helpers[0] = ("evil.module", "evil_qualname")  # type: ignore


def test_expected_mode_dispatch_nested_artifact_source_frozen(compare_mod):
    with pytest.raises(TypeError):
        compare_mod.EXPECTED_MODE_DISPATCH["g1"]["artifact_source"] = "evil"  # type: ignore


# ---------------------------------------------------------------------
# _validate_dispatch returns the requested_mode entry (no shadowing).
# ---------------------------------------------------------------------


def test_validate_dispatch_returns_requested_mode_entry(compare_mod):
    for mode in compare_mod.EXPECTED_MODE_DISPATCH:
        entry, resolved = compare_mod._validate_dispatch(mode)
        expected = compare_mod.EXPECTED_MODE_DISPATCH[mode]
        assert tuple(entry["helpers"]) == expected["helpers"]
        assert entry["artifact_source"] == expected["artifact_source"]
        # Every canonical helper resolved to its own module.
        for (module_dotted, qualname), obj in zip(expected["helpers"], resolved):
            assert obj.__module__ == module_dotted
            assert obj.__qualname__ == qualname


def test_validate_dispatch_unknown_mode_exits_2(compare_mod, capsys):
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g_bogus")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "unknown requested_mode" in err


# ---------------------------------------------------------------------
# MODE_DISPATCH corruption catches.
# ---------------------------------------------------------------------


def test_validate_dispatch_corrupt_helpers_exits_2(compare_mod, capsys):
    compare_mod.MODE_DISPATCH["g1"]["helpers"] = [
        ("tools._uhfk_to_mvmc.general_fij_builder", "build_slater_orbitals"),
    ]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g1")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "MODE_DISPATCH corrupted" in err
    assert "helpers mismatch" in err


def test_validate_dispatch_corrupt_artifact_source_exits_2(compare_mod, capsys):
    compare_mod.MODE_DISPATCH["g3"]["artifact_source"] = "hwave+evil"
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g3")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "MODE_DISPATCH corrupted" in err
    assert "artifact_source mismatch" in err


def test_validate_dispatch_missing_mode_exits_2(compare_mod, capsys):
    del compare_mod.MODE_DISPATCH["g2b"]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g2b")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "unknown requested_mode" in err or "missing mode" in err


def test_validate_dispatch_missing_other_mode_still_exits_2(compare_mod, capsys):
    """Corrupt a mode OTHER than the requested one; the structural check
    still fires because it iterates the full expected table."""
    del compare_mod.MODE_DISPATCH["g3"]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g1")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "missing mode=g3" in err


# ---------------------------------------------------------------------
# Wrong-module / missing-qualname / same-qualname-wrong-module.
# ---------------------------------------------------------------------


def test_validate_dispatch_wrong_module_exits_2(compare_mod, capsys):
    """Corrupt MODE_DISPATCH so the helper resolves via an IMPORTABLE
    module that does NOT expose the qualname; the getattr(...) returns
    None branch fires exit 2 with the ``has no attribute`` message."""
    compare_mod.MODE_DISPATCH["g1"]["helpers"] = [
        ("tools._uhfk_to_mvmc.energy_compare", "build_slater_orbitals"),
        ("tools._uhfk_to_mvmc.density_check",       "gauge_lift"),
    ]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g1")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "MODE_DISPATCH corrupted" in err or "has no attribute" in err


def test_validate_dispatch_helper_missing_exits_2(compare_mod, capsys):
    compare_mod.MODE_DISPATCH["g2b"]["helpers"] = [
        ("tools._uhfk_to_mvmc.density_check", "does_not_exist_qualname"),
    ]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g2b")
    assert ex.value.code == 2
    err = capsys.readouterr().err
    assert "MODE_DISPATCH corrupted" in err or "has no attribute" in err


def test_validate_dispatch_same_qualname_wrong_module_exits_2(
    compare_mod, capsys
):
    """A stub module exposes the same
    qualname the table expects but from a hostile module path. The
    obj.__module__ check MUST fire exit 2."""
    stub = types.ModuleType("_wiring_stub_wrong_module_")

    def gauge_lift(*args, **kwargs):
        return 0.0j

    stub.gauge_lift = gauge_lift
    sys.modules["_wiring_stub_wrong_module_"] = stub
    try:
        compare_mod.MODE_DISPATCH["g2b"]["helpers"] = [
            ("_wiring_stub_wrong_module_", "gauge_lift"),
        ]
        with pytest.raises(SystemExit) as ex:
            compare_mod._validate_dispatch("g2b")
        assert ex.value.code == 2
        err = capsys.readouterr().err
        assert "corrupted" in err or "wrong-module dispatch" in err
    finally:
        del sys.modules["_wiring_stub_wrong_module_"]


# ---------------------------------------------------------------------
# Composite-gate: missing-helper corruption catches for G1 and G2a-emitted-F.
# ---------------------------------------------------------------------


def test_compare_composite_gate_g1_missing_helper_exits_2(
    compare_mod, capsys
):
    # Drop the second canonical helper (gauge_lift) from g1's tuple.
    compare_mod.MODE_DISPATCH["g1"]["helpers"] = [
        ("tools._uhfk_to_mvmc.general_fij_builder", "build_slater_orbitals"),
    ]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g1")
    assert ex.value.code == 2


def test_compare_composite_gate_g2a_emitted_F_missing_helper_exits_2(
    compare_mod, capsys
):
    # Drop parse_emitted_F, leave only pair_product_density_from_F.
    compare_mod.MODE_DISPATCH["g2a-emitted-F"]["helpers"] = [
        ("tools._uhfk_to_mvmc.pair_product_density",
         "pair_product_density_from_F"),
    ]
    with pytest.raises(SystemExit) as ex:
        compare_mod._validate_dispatch("g2a-emitted-F")
    assert ex.value.code == 2


# ---------------------------------------------------------------------
# Numerical gate failures must fail closed on invalid comparison inputs.
# ---------------------------------------------------------------------


def _run_g0_writer_check(compare_mod, tmp_path, emitted, reference, gtol):
    bridge = tmp_path / "bridge_zeronoise"
    bridge.mkdir(exist_ok=True)
    np.savez(bridge / "F_pre_noise.npz", F=np.asarray(reference))

    def parse_emitted_F(_path):
        return np.asarray(emitted)

    return compare_mod._dispatch_g0_writer_check(
        str(tmp_path), (parse_emitted_F,), gtol
    )


@pytest.mark.parametrize(
    ("emitted", "reference"),
    [([np.nan], [0.0]), ([0.0], [np.nan])],
)
def test_g0_writer_check_rejects_non_finite_array(
    compare_mod, tmp_path, capsys, emitted, reference
):
    rc = _run_g0_writer_check(
        compare_mod, tmp_path, emitted, reference, 1.0e-10
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G0-writer-check PASS" not in captured.out
    assert "non-finite" in captured.err


@pytest.mark.parametrize("gtol", [np.nan, np.inf, -np.inf, 0.0, -1.0])
def test_g0_writer_check_rejects_invalid_gtol(
    compare_mod, tmp_path, capsys, gtol
):
    rc = _run_g0_writer_check(
        compare_mod, tmp_path, [1.0], [0.0], gtol
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G0-writer-check PASS" not in captured.out
    assert "tolerance" in captured.err


def test_g0_writer_check_normal_match_passes(compare_mod, tmp_path, capsys):
    rc = _run_g0_writer_check(
        compare_mod, tmp_path, [0.0], [0.0], 1.0e-10
    )
    captured = capsys.readouterr()
    assert rc == 0
    assert "G0-writer-check PASS" in captured.out


def test_g0_writer_check_normal_mismatch_fails(compare_mod, tmp_path, capsys):
    rc = _run_g0_writer_check(
        compare_mod, tmp_path, [1.0], [0.0], 1.0e-10
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G0-writer-check PASS" not in captured.out


def test_g0_writer_check_rejects_non_finite_computed_delta(
    compare_mod, tmp_path, capsys
):
    largest = np.finfo(np.float64).max
    with np.errstate(over="ignore"):
        rc = _run_g0_writer_check(
            compare_mod, tmp_path, [largest], [-largest], 1.0e-10
        )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G0-writer-check PASS" not in captured.out
    assert "max absolute delta is non-finite" in captured.err


def _run_g1(
    compare_mod, tmp_path, monkeypatch, shipping_orbitals, lifted_density, gtol
):
    monkeypatch.setattr(compare_mod, "_load_workspace_config", lambda _path: {})
    monkeypatch.setattr(
        compare_mod,
        "_build_shipping_A",
        lambda _workspace, _cfg: (np.asarray(shipping_orbitals), None),
    )
    monkeypatch.setattr(
        compare_mod,
        "_build_G_from_gauge_lift",
        lambda _workspace, _cfg, _gauge_lift: np.asarray(lifted_density),
    )
    return compare_mod._dispatch_g1(
        str(tmp_path), (object(), object()), gtol
    )


@pytest.mark.parametrize(
    ("shipping_orbitals", "lifted_density"),
    [([[np.nan]], [[0.0]]), ([[0.0]], [[np.nan]])],
)
def test_g1_rejects_non_finite_array(
    compare_mod,
    tmp_path,
    monkeypatch,
    capsys,
    shipping_orbitals,
    lifted_density,
):
    rc = _run_g1(
        compare_mod,
        tmp_path,
        monkeypatch,
        shipping_orbitals,
        lifted_density,
        1.0e-10,
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G1 PASS" not in captured.out
    assert "non-finite" in captured.err


@pytest.mark.parametrize("gtol", [np.nan, np.inf, -np.inf, 0.0, -1.0])
def test_g1_rejects_invalid_gtol(
    compare_mod, tmp_path, monkeypatch, capsys, gtol
):
    rc = _run_g1(
        compare_mod, tmp_path, monkeypatch, [[1.0]], [[0.0]], gtol
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G1 PASS" not in captured.out
    assert "tolerance" in captured.err


def test_g1_normal_match_passes(compare_mod, tmp_path, monkeypatch, capsys):
    rc = _run_g1(
        compare_mod, tmp_path, monkeypatch, [[1.0]], [[1.0]], 1.0e-10
    )
    captured = capsys.readouterr()
    assert rc == 0
    assert "G1 PASS" in captured.out


def test_g1_normal_mismatch_fails(compare_mod, tmp_path, monkeypatch, capsys):
    rc = _run_g1(
        compare_mod, tmp_path, monkeypatch, [[1.0]], [[0.0]], 1.0e-10
    )
    captured = capsys.readouterr()
    assert rc != 0
    assert "G1 PASS" not in captured.out


# ---------------------------------------------------------------------
# metadata-preserving monkey-patch verifies dispatch actually calls the
# canonical helper.
# ---------------------------------------------------------------------


def _write_minimal_g2b_workspace(root, Ns=2):
    """Build a workspace bare enough for compare.py --mode g2b to
    complete once monkey-patches force gauge_lift's output to zero.

    Uses a tiny 2-site PBC SOC fixture (CellShape=[2,1,1], SubShape=[1,1,1])
    so parse_input.toml / load_geometry_uhf / _resolve_g3_paths style
    loaders all succeed but the numerical work is trivial. The sentinel
    is expected to return 0.0 for every (i, s, j, t) so G_lift == 0 and
    the ComplexUHF file we write is all zeros, giving max_abs_delta = 0."""
    import numpy as np

    hwave = root / "hwave"
    hwave.mkdir(parents=True, exist_ok=True)
    complexuhf = root / "complexuhf"
    complexuhf.mkdir(parents=True, exist_ok=True)

    # Minimal input.toml.
    (hwave / "input.toml").write_text(
        "[log]\nprint_level = 0\n\n[mode]\nmode = \"UHFk\"\nenable_spin_orbital = true\n"
        "[mode.param]\nNcond = 2\nCellShape = [2, 1, 1]\nSubShape = [1, 1, 1]\n"
        "BoundaryCondition = [\"periodic\", \"periodic\", \"periodic\"]\n"
        "IterationMax = 1\nEPS = 8\nMix = 0.5\nRndSeed = 12345\nT = 0.0\n"
        "[file]\n[file.output]\npath_to_output = \"output\"\n"
    )
    # Minimal geometry_uhf.dat (2 sites, 1 orbital).
    (hwave / "geometry_uhf.dat").write_text(
        "  1.000000000000   0.000000000000   0.000000000000\n"
        "  0.000000000000   1.000000000000   0.000000000000\n"
        "  0.000000000000   0.000000000000   1.000000000000\n"
        "  0.000000000000   0.000000000000   0.000000000000\n"
        "  2 0 0\n  0 1 0\n  0 0 1\n"
        "  0 0 0 0\n  1 0 0 0\n"
    )
    # Minimal green.npz. gauge_lift is monkey-patched so the numerical
    # content does not matter beyond having the right shape and dtype.
    green_sub = np.zeros((2, 1, 4, 1, 4), dtype=np.complex128)
    np.savez(hwave / "green.npz", green_sublattice=green_sub)
    # G2 contraction precondition: start clearly outside the default
    # 1e-6 tolerance, then use the all-zero converged artifact below.
    (complexuhf / "initial.def").write_text(
        "=============================================\n"
        "NInitial          1\n"
        "=============================================\n"
        "=============================================\n"
        "=============================================\n"
        "0 0 0 0 1.0e-3 0.0\n"
    )
    # ComplexUHF density (all zeros, matching the sentinel's output).
    # The strict parser requires exactly (2*Ns)^2 = 16 (i, s, j, t)
    # rows for Ns = 2.
    with open(complexuhf / "zvo_UHF_cisajs.dat", "w") as fp:
        for i in range(Ns):
            for s in (0, 1):
                for j in range(Ns):
                    for t in (0, 1):
                        fp.write(
                            f"    {i}    {s}    {j}    {t}     "
                            "0.0000000000     0.0000000000\n"
                        )


def test_compare_g2b_dispatch_invokes_canonical_gauge_lift(
    compare_mod, tmp_path, monkeypatch, capsys
):
    """Wrap the canonical density_check.gauge_lift with a
    metadata-preserving sentinel; run compare.py --mode g2b; assert the
    sentinel was invoked (via a side-effect counter) AND the emitted
    PASS record still carries the canonical helper= from the table.

    The real G2b dispatcher performs actual artefact I/O and numeric
    comparison, so the wiring test provides a minimal synthetic workspace whose
    gauge-lifted density is defined by the sentinel and whose
    ComplexUHF reference is zero. The test therefore still isolates
    the dispatch metadata contract from the numerical physics."""
    from tools._uhfk_to_mvmc import density_check

    _write_minimal_g2b_workspace(tmp_path)

    original = density_check.gauge_lift
    calls = {"n": 0}

    def _sentinel(*args, **kwargs):
        calls["n"] += 1
        return 0.0 + 0.0j

    _sentinel.__module__ = original.__module__
    _sentinel.__qualname__ = original.__qualname__
    _sentinel.__name__ = original.__name__

    monkeypatch.setattr(density_check, "gauge_lift", _sentinel)
    exit_code = compare_mod.main(
        ["compare.py", "--workspace", str(tmp_path), "--mode", "g2b"],
    )
    assert exit_code == 0, capsys.readouterr()
    out = capsys.readouterr().out
    assert (
        "G2b PASS mode=g2b artifact_source=hwave+complexuhf "
        "helper=gauge_lift " in out
    )
    # Real G2b dispatcher sweeps every (i, s, j, t) — the sentinel MUST
    # have been called at least once (was 0 under the pre-hardening
    # stub). This confirms dispatch actually invokes the canonical
    # helper on the workspace data.
    assert calls["n"] >= 1


def test_compare_g3_invokes_canonical_energy_relative_delta(
    compare_mod, tmp_path, monkeypatch, capsys
):
    # Build a minimal G3 workspace.
    hwave_dir = tmp_path / "hwave"
    hwave_dir.mkdir()
    (hwave_dir / "energy.dat").write_text("Energy_Total = -1.0\n")
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    (mvmc_dir / "zvo_out_selected.dat").write_text(
        "-1.0 1.0 0 0 0 0\n"
    )

    from tools._uhfk_to_mvmc import energy_compare

    original = energy_compare.energy_relative_delta
    calls = {"n": 0, "last_args": None}

    def _sentinel(hwave_energy_dat, mvmc_zvo_out_dat):
        calls["n"] += 1
        calls["last_args"] = (hwave_energy_dat, mvmc_zvo_out_dat)
        return original(hwave_energy_dat, mvmc_zvo_out_dat)

    _sentinel.__module__ = original.__module__
    _sentinel.__qualname__ = original.__qualname__
    _sentinel.__name__ = original.__name__

    monkeypatch.setattr(energy_compare, "energy_relative_delta", _sentinel)
    exit_code = compare_mod.main(
        ["compare.py", "--workspace", str(tmp_path), "--mode", "g3"],
    )
    assert exit_code == 0
    assert calls["n"] >= 1
    out = capsys.readouterr().out
    assert (
        "G3 PASS mode=g3 artifact_source=hwave+mvmc "
        "helper=energy_relative_delta " in out
    )


@pytest.mark.parametrize(
    ("zvo_row", "diagnostic"),
    [
        ("-1.25\n", "expected exactly 6 columns"),
        ("-1.25 nan 0 0 0 0\n", "non-finite"),
    ],
)
def test_compare_g3_rejects_invalid_zvo_row(
    compare_mod, tmp_path, capsys, zvo_row, diagnostic
):
    hwave_dir = tmp_path / "hwave"
    hwave_dir.mkdir()
    (hwave_dir / "energy.dat").write_text("Energy_Total = -1.25\n")
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    (mvmc_dir / "zvo_out_selected.dat").write_text(zvo_row)

    exit_code = compare_mod.main(
        ["compare.py", "--workspace", str(tmp_path), "--mode", "g3"],
    )
    captured = capsys.readouterr()
    assert exit_code != 0
    assert "G3 PASS" not in captured.out
    assert diagnostic in captured.err


@pytest.mark.parametrize("bad_sample", ["nan", "inf", "-inf"])
def test_compare_g3_rejects_non_finite_mvmc_energy(
    compare_mod, tmp_path, capsys, bad_sample
):
    """A non-finite mVMC sample must fail G3, not pass it.

    `nan > tol` is False, so before the finiteness guard a NaN sailed
    through to the anchored PASS record and silently disabled the only
    H-wave/mVMC energy-agreement gate. run.sh accepts that record by
    prefix alone, so nothing downstream would have caught it.
    """
    hwave_dir = tmp_path / "hwave"
    hwave_dir.mkdir()
    (hwave_dir / "energy.dat").write_text("Energy_Total = -1.0\n")
    mvmc_dir = tmp_path / "mvmc"
    mvmc_dir.mkdir()
    (mvmc_dir / "zvo_out_selected.dat").write_text(
        f"{bad_sample} 1.0 0 0 0 0\n"
    )

    exit_code = compare_mod.main(
        ["compare.py", "--workspace", str(tmp_path), "--mode", "g3"],
    )
    captured = capsys.readouterr()
    assert exit_code != 0
    assert "G3 PASS" not in captured.out
    assert "non-finite" in captured.err


# ---------------------------------------------------------------------
# CLI edge cases.
# ---------------------------------------------------------------------


def test_legacy_positional_args_still_work(compare_mod, tmp_path, capsys):
    """The 3-positional-arg invocation MUST continue to work for
    existing PBC / APBC fixtures (backward compat)."""
    hwave = tmp_path / "energy.dat"
    hwave.write_text("Energy_Total = -1.0\n")
    mvmc = tmp_path / "zvo_out_001.dat"
    mvmc.write_text("-1.0 1.0 0 0 0 0\n")
    rc = compare_mod.main(
        ["compare.py", str(hwave), str(mvmc), "unit_test_case"],
    )
    assert rc == 0
    out = capsys.readouterr().out
    assert "unit_test_case" in out


# ---------------------------------------------------------------------
# The ComplexUHF parser must fail closed on empty, partial, duplicate,
# out-of-range, malformed, or non-finite input. Returning zeros for an
# empty file would give spurious G2a/G2b PASS records.
# ---------------------------------------------------------------------


def test_load_complexuhf_G_rejects_empty_file(compare_mod, tmp_path):
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    (complexuhf / "zvo_UHF_cisajs.dat").write_text("")
    with pytest.raises(compare_mod.ComplexUHFParseError, match="expected"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)


def test_load_complexuhf_G_rejects_truncated_file(compare_mod, tmp_path):
    """Two entries out of the required 16 for Ns=2."""
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    (complexuhf / "zvo_UHF_cisajs.dat").write_text(
        "0 0 0 0 0.0 0.0\n0 1 0 0 0.0 0.0\n"
    )
    with pytest.raises(compare_mod.ComplexUHFParseError, match="expected"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)


def test_load_complexuhf_G_rejects_duplicate_entry(compare_mod, tmp_path):
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    # 16 rows for Ns=2 with (0, 0, 0, 0) duplicated.
    rows = ["0 0 0 0 0.0 0.0\n"] * 2
    for i in range(2):
        for s in (0, 1):
            for j in range(2):
                for t in (0, 1):
                    if (i, s, j, t) == (0, 0, 0, 0):
                        continue
                    rows.append(f"{i} {s} {j} {t} 0.0 0.0\n")
    (complexuhf / "zvo_UHF_cisajs.dat").write_text("".join(rows))
    with pytest.raises(compare_mod.ComplexUHFParseError, match="duplicate"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)


def test_load_complexuhf_G_rejects_out_of_range_site(compare_mod, tmp_path):
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    # Ns=2 → sites {0, 1}; site=5 is out of range.
    (complexuhf / "zvo_UHF_cisajs.dat").write_text(
        "5 0 0 0 0.0 0.0\n"
    )
    with pytest.raises(compare_mod.ComplexUHFParseError, match="out of range"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)


def test_load_complexuhf_G_rejects_malformed_line(compare_mod, tmp_path):
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    (complexuhf / "zvo_UHF_cisajs.dat").write_text(
        "not_a_line\n"
    )
    with pytest.raises(compare_mod.ComplexUHFParseError, match="6 tokens"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)


def test_load_complexuhf_G_rejects_non_finite(compare_mod, tmp_path):
    complexuhf = tmp_path / "complexuhf"
    complexuhf.mkdir()
    (complexuhf / "zvo_UHF_cisajs.dat").write_text(
        "0 0 0 0 nan 0.0\n"
    )
    with pytest.raises(compare_mod.ComplexUHFParseError, match="non-finite"):
        compare_mod._load_complexuhf_G(str(tmp_path), Ns=2)
