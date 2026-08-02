"""Mutation schema-detection tests for
tests/validation/uhfk_mvmc_pairproduct/scripts/soc_apbc_topology_guard.py.
"""
import copy
import importlib.util
import json
import os
from pathlib import Path

import numpy as np
import pytest


_REPO = Path(__file__).resolve().parents[1]
_GUARD_PATH = str(
    _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "scripts" / "soc_apbc_topology_guard.py"
)
_PRODUCER_PATH = str(
    _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
    / "scripts" / "phase2_produce_manifest.py"
)
_DATA_DIR = _REPO / "tests" / "data"


@pytest.fixture(scope="module")
def guard_mod():
    spec = importlib.util.spec_from_file_location(
        "_v37_topology_guard_schema_under_test", _GUARD_PATH,
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def producer_mod():
    spec = importlib.util.spec_from_file_location(
        "_v37_phase2_manifest_producer_under_test", _PRODUCER_PATH,
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _site_positions(cell_shape):
    return np.array(
        [
            [x, y, z]
            for z in range(int(cell_shape[2]))
            for y in range(int(cell_shape[1]))
            for x in range(int(cell_shape[0]))
        ],
        dtype=np.int64,
    )


def _manifest_and_geometry(prefix):
    with open(_DATA_DIR / f"{prefix}_composite_element.json") as fp:
        manifest = json.load(fp)
    cell_shape = np.asarray(manifest["cell_shape"], dtype=np.int64)
    subshape = np.asarray(manifest["sub_shape"], dtype=np.int64)
    theta = np.asarray(manifest["theta_radians"], dtype=np.float64)
    return manifest, cell_shape, subshape, theta, _site_positions(cell_shape)


def test_iter_mutation_ids_returns_v37_on_30_key_manifest(guard_mod):
    manifest = {
        "T_M_per_mutation": {mid: 0.5 for mid in guard_mod._V37_MUTATION_IDS},
    }
    ids = guard_mod._iter_mutation_ids(manifest)
    assert ids == guard_mod._V37_MUTATION_IDS


def test_iter_mutation_ids_returns_v36_on_10_key_manifest(guard_mod):
    manifest = {
        "T_M_per_mutation": {mid: 0.5 for mid in guard_mod._V36_MUTATION_IDS},
    }
    ids = guard_mod._iter_mutation_ids(manifest)
    assert ids == guard_mod._V36_MUTATION_IDS


def test_iter_mutation_ids_rejects_hybrid_manifest(guard_mod):
    """A hybrid manifest with 15 keys MUST raise, not match either schema."""
    partial_v37 = list(guard_mod._V37_MUTATION_IDS)[:5]
    manifest = {"T_M_per_mutation": {}}
    for mid in guard_mod._V36_MUTATION_IDS:
        manifest["T_M_per_mutation"][mid] = 0.5
    for mid in partial_v37:
        manifest["T_M_per_mutation"][mid] = 0.5
    with pytest.raises(ValueError):
        guard_mod._iter_mutation_ids(manifest)


def test_iter_mutation_ids_rejects_non_dict_T_M(guard_mod):
    """Non-dict T_M_per_mutation MUST raise ValueError, not AttributeError."""
    for bad_value in (None, [], "not a dict", 42):
        manifest = {"T_M_per_mutation": bad_value}
        with pytest.raises(ValueError, match="dict"):
            guard_mod._iter_mutation_ids(manifest)


def test_iter_mutation_ids_rejects_missing_T_M(guard_mod):
    """Manifest without T_M_per_mutation MUST raise ValueError."""
    manifest = {"i_c": 0, "s_c": 0, "j_c": 0, "t_c": 0, "G_c_abs": 1.0}
    with pytest.raises(ValueError, match="dict"):
        guard_mod._iter_mutation_ids(manifest)


def test_split_mutator_id_v36_forms(guard_mod):
    """Whole-vector ids return axis=None."""
    base, axis = guard_mod._split_mutator_id("M-gauge-1")
    assert base == "M-gauge-1"
    assert axis is None
    base, axis = guard_mod._split_mutator_id("M-ship-5")
    assert base == "M-ship-5"
    assert axis is None


def test_split_mutator_id_v37_forms(guard_mod):
    """Per-direction ids return the correct base and axis."""
    for axis_idx, axis_char in enumerate("xyz"):
        base, axis = guard_mod._split_mutator_id(f"M-gauge-1-{axis_char}")
        assert base == "M-gauge-1"
        assert axis == axis_idx
        base, axis = guard_mod._split_mutator_id(f"M-ship-5-{axis_char}")
        assert base == "M-ship-5"
        assert axis == axis_idx


def test_split_mutator_id_baseline_form(guard_mod):
    """The 'baseline' pseudo-id (used by non-mutation guard paths) has
    axis=None."""
    base, axis = guard_mod._split_mutator_id("baseline")
    assert base == "baseline"
    assert axis is None


def test_split_mutator_id_invalid_ids_are_not_whitelisted(guard_mod):
    """Invalid mutator ids (out-of-range gauge/ship index, or a bad
    axis suffix) MUST split into a base_id that is NOT one of the
    known M-gauge-1..5 / M-ship-1..5 base ids. _split_mutator_id itself
    does not validate; it is the caller's whitelist check (the
    ``base_id not in (...)`` guard in _gauge_lift_element and
    _build_A_ship_mutated) that rejects these downstream, so this test
    pins the split output that whitelist check depends on."""
    whitelist = set(guard_mod._V36_MUTATION_IDS)
    # Out-of-range index (only 1..5 are valid): axis suffix still
    # parses, so the base id alone must fail the whitelist.
    base, axis = guard_mod._split_mutator_id("M-gauge-6-x")
    assert base == "M-gauge-6"
    assert axis == 0
    assert base not in whitelist
    # Malformed axis suffix ("w" is not x/y/z): falls through to the
    # whole-string-as-base branch, which also fails the whitelist.
    base, axis = guard_mod._split_mutator_id("M-gauge-1-w")
    assert base == "M-gauge-1-w"
    assert axis is None
    assert base not in whitelist


@pytest.mark.parametrize("base_id", ["M-gauge-4", "M-ship-4"])
def test_v37_m4_is_not_structurally_degenerate_on_l_folded_two(
    producer_mod, base_id,
):
    assert not producer_mod._is_structurally_degenerate_mutator(
        base_id,
        axis=0,
        L_folded=np.array([2, 2, 2], dtype=np.int64),
        sub=np.array([2, 2, 2], dtype=np.int64),
    )


@pytest.mark.parametrize("base_id", ["M-gauge-4", "M-ship-4"])
def test_v37_m4_remains_degenerate_without_folded_or_sublattice_extent(
    producer_mod, base_id,
):
    assert producer_mod._is_structurally_degenerate_mutator(
        base_id,
        axis=0,
        L_folded=np.array([1, 2, 2], dtype=np.int64),
        sub=np.array([2, 2, 2], dtype=np.int64),
    )
    assert producer_mod._is_structurally_degenerate_mutator(
        base_id,
        axis=0,
        L_folded=np.array([2, 2, 2], dtype=np.int64),
        sub=np.array([1, 2, 2], dtype=np.int64),
    )


def test_manifest_producer_rejects_structurally_degenerate_active_axis(
    producer_mod,
):
    with pytest.raises(RuntimeError, match="inadequate.*mutation.*lattice"):
        producer_mod._mutation_threshold(
            "M-gauge-4-x",
            axis_char="x",
            is_active=True,
            is_degenerate=True,
            threshold_floor=0.01,
        )


def test_manifest_producer_rejects_zero_threshold_on_active_axis(producer_mod):
    with pytest.raises(RuntimeError, match="thresholds must be positive"):
        producer_mod._mutation_threshold(
            "M-gauge-1-x",
            axis_char="x",
            is_active=True,
            is_degenerate=False,
            threshold_floor=0.0,
        )


def test_manifest_producer_rejects_nan_threshold_on_active_axis(producer_mod):
    with pytest.raises(RuntimeError, match="finite"):
        producer_mod._mutation_threshold(
            "M-gauge-1-x",
            axis_char="x",
            is_active=True,
            is_degenerate=False,
            threshold_floor=float("nan"),
        )


def test_manifest_producer_keeps_zero_threshold_on_inactive_axis(producer_mod):
    assert producer_mod._mutation_threshold(
        "M-gauge-4-z",
        axis_char="z",
        is_active=False,
        is_degenerate=True,
        threshold_floor=0.01,
    ) == 0.0


@pytest.mark.parametrize("suffix", ["xy", "xz", "yz", "xyz"])
def test_shipped_v37_manifests_have_positive_active_axis_thresholds(suffix):
    case_dir = (
        _REPO / "tests" / "validation" / "uhfk_mvmc_pairproduct"
        / f"case_soc_rashba_3d_sub_apbc_{suffix}"
    )
    with open(case_dir / "composite_element.json") as fp:
        manifest = json.load(fp)
    active_axes = set(manifest["active_apbc_axes"])
    for mutator_id, threshold in manifest["T_M_per_mutation"].items():
        if mutator_id[-1] in active_axes:
            assert threshold > 0.0, mutator_id


def test_v37_gauge_m4_has_nonzero_delta_on_l_folded_two(guard_mod):
    prefix = "v37_case_soc_rashba_3d_sub_apbc_xy"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    green = np.load(_DATA_DIR / f"{prefix}_green.npz")["green_sublattice"]
    args = (
        green,
        int(manifest["i_c"]),
        int(manifest["s_c"]),
        int(manifest["j_c"]),
        int(manifest["t_c"]),
    )
    kwargs = {
        "subshape": sub,
        "cell_shape": cell,
        "site_positions": positions,
        "boundary_theta": theta,
    }
    baseline = guard_mod._gauge_lift_element(
        *args, **kwargs, mutator="baseline",
    )
    mutated = guard_mod._gauge_lift_element(
        *args, **kwargs, mutator="M-gauge-4-x",
    )
    assert abs(mutated - baseline) > 1e-12


def test_v37_shipping_m4_has_nonzero_delta_on_l_folded_two(guard_mod):
    prefix = "v37_case_soc_rashba_3d_sub_apbc_xy"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    eigen = np.load(_DATA_DIR / f"{prefix}_eigen.npz")
    occ = np.load(_DATA_DIR / f"{prefix}_occupation.npz")
    build_args = (eigen, occ, positions, cell, sub, int(manifest["ncond"]), theta)
    baseline = guard_mod._build_A_ship_mutated("baseline", *build_args)
    mutated = guard_mod._build_A_ship_mutated("M-ship-4-x", *build_args)
    nsite = int(np.prod(cell))
    all_i = int(manifest["i_c"]) + int(manifest["s_c"]) * nsite
    all_j = int(manifest["j_c"]) + int(manifest["t_c"]) * nsite
    baseline_g = (np.conj(baseline) @ baseline.T)[all_i, all_j]
    mutated_g = (np.conj(mutated) @ mutated.T)[all_i, all_j]
    assert abs(mutated_g - baseline_g) > 1e-12


@pytest.mark.parametrize("suffix", ["xy", "xz", "yz", "xyz"])
def test_v37_shipping_m4_differs_from_m5_on_every_active_axis(
    guard_mod, suffix,
):
    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{suffix}"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    eigen = np.load(_DATA_DIR / f"{prefix}_eigen.npz")
    occ = np.load(_DATA_DIR / f"{prefix}_occupation.npz")
    build_args = (
        eigen, occ, positions, cell, sub, int(manifest["ncond"]), theta,
    )

    for axis, axis_char in enumerate("xyz"):
        if abs(float(theta[axis])) <= 1e-12:
            continue
        m4 = guard_mod._build_A_ship_mutated(
            f"M-ship-4-{axis_char}", *build_args,
        )
        m5 = guard_mod._build_A_ship_mutated(
            f"M-ship-5-{axis_char}", *build_args,
        )
        max_abs_diff = float(np.max(np.abs(m4 - m5)))
        assert max_abs_diff > 0.0, (
            f"{suffix} active axis {axis_char}: M-ship-4 and "
            "M-ship-5 produced identical A matrices"
        )


@pytest.mark.parametrize("suffix", ["xy", "xz", "yz", "xyz"])
def test_v37_gauge_m4_differs_from_m5_on_every_active_axis(
    guard_mod, suffix,
):
    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{suffix}"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    green = np.load(_DATA_DIR / f"{prefix}_green.npz")["green_sublattice"]
    args = (
        green,
        int(manifest["i_c"]),
        int(manifest["s_c"]),
        int(manifest["j_c"]),
        int(manifest["t_c"]),
    )
    kwargs = {
        "subshape": sub,
        "cell_shape": cell,
        "site_positions": positions,
        "boundary_theta": theta,
    }

    for axis, axis_char in enumerate("xyz"):
        if abs(float(theta[axis])) <= 1e-12:
            continue
        m4 = guard_mod._gauge_lift_element(
            *args, **kwargs, mutator=f"M-gauge-4-{axis_char}",
        )
        m5 = guard_mod._gauge_lift_element(
            *args, **kwargs, mutator=f"M-gauge-5-{axis_char}",
        )
        assert abs(m4 - m5) > 0.0, (
            f"{suffix} active axis {axis_char}: M-gauge-4 and "
            "M-gauge-5 produced the same composite value"
        )


@pytest.mark.parametrize(
    ("mutator_id", "threshold"),
    [
        ("M-gauge-1-x", float("nan")),
        ("M-gauge-1-x", 1e-300),
        ("M-gauge-1-x", 0.0),
        ("M-gauge-1-z", 1e-5),
    ],
)
def test_manifest_producer_rejects_threshold_policy_drift_before_writing(
    producer_mod, tmp_path, capsys, mutator_id, threshold,
):
    manifest, _, _, _, _ = _manifest_and_geometry(
        "v37_case_soc_rashba_3d_sub_apbc_xy",
    )
    drifted = copy.deepcopy(manifest)
    drifted["T_M_per_mutation"][mutator_id] = threshold
    output = tmp_path / "composite_element.json"

    exit_code = producer_mod._write_manifest_after_self_check(
        drifted, output,
    )

    assert exit_code != 0
    assert not output.exists()
    assert "threshold" in capsys.readouterr().err.lower()


def test_manifest_producer_self_check_failure_returns_nonzero(
    producer_mod, tmp_path, capsys,
):
    manifest, _, _, _, _ = _manifest_and_geometry(
        "v37_case_soc_rashba_3d_sub_apbc_xy",
    )
    subthreshold = copy.deepcopy(manifest)
    subthreshold["delta_M_per_mutation_at_producer_time"][
        "M-gauge-1-x"
    ] = 0.0
    output = tmp_path / "composite_element.json"

    exit_code = producer_mod._write_manifest_after_self_check(
        subthreshold, output,
    )

    assert exit_code != 0
    assert not output.exists()
    assert "self-check" in capsys.readouterr().err.lower()


def test_v36_whole_vector_gauge_m4_matches_frozen_manifest(guard_mod):
    prefix = "v36_case_soc_rashba_2d_sub_apbc"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    green = np.load(_DATA_DIR / f"{prefix}_green.npz")["green_sublattice"]
    args = (
        green,
        int(manifest["i_c"]),
        int(manifest["s_c"]),
        int(manifest["j_c"]),
        int(manifest["t_c"]),
    )
    kwargs = {
        "subshape": sub,
        "cell_shape": cell,
        "site_positions": positions,
        "boundary_theta": theta,
    }
    baseline = guard_mod._gauge_lift_element(
        *args, **kwargs, mutator="baseline",
    )
    mutated = guard_mod._gauge_lift_element(
        *args, **kwargs, mutator="M-gauge-4",
    )
    expected = manifest["delta_M_per_mutation_at_producer_time"]["M-gauge-4"]
    assert abs(mutated - baseline) == pytest.approx(expected, abs=1e-12)


def test_v36_whole_vector_shipping_m4_matches_frozen_manifest(guard_mod):
    prefix = "v36_case_soc_rashba_2d_sub_apbc"
    manifest, cell, sub, theta, positions = _manifest_and_geometry(prefix)
    eigen = np.load(_DATA_DIR / f"{prefix}_eigen.npz")
    occ = np.load(_DATA_DIR / f"{prefix}_occupation.npz")
    build_args = (eigen, occ, positions, cell, sub, int(manifest["ncond"]), theta)
    baseline = guard_mod._build_A_ship_mutated("baseline", *build_args)
    mutated = guard_mod._build_A_ship_mutated("M-ship-4", *build_args)
    nsite = int(np.prod(cell))
    all_i = int(manifest["i_c"]) + int(manifest["s_c"]) * nsite
    all_j = int(manifest["j_c"]) + int(manifest["t_c"]) * nsite
    baseline_g = (np.conj(baseline) @ baseline.T)[all_i, all_j]
    mutated_g = (np.conj(mutated) @ mutated.T)[all_i, all_j]
    expected = manifest["delta_M_per_mutation_at_producer_time"]["M-ship-4"]
    assert abs(mutated_g - baseline_g) == pytest.approx(expected, abs=1e-12)
