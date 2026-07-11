import logging
import pytest
import hwave.tsweep as ts


def test_explicit_temperatures_pass_through():
    assert ts.build_ladder({"temperatures": [0.01, 0.008, 0.006]}) == [0.01, 0.008, 0.006]


def test_generated_linear_ladder():
    out = ts.build_ladder({"T_start": 0.01, "T_stop": 0.002, "num": 5, "spacing": "linear"})
    assert out[0] == pytest.approx(0.01)
    assert out[-1] == pytest.approx(0.002)
    assert len(out) == 5
    assert all(out[i] > out[i + 1] for i in range(len(out) - 1))


def test_generated_log_ladder():
    out = ts.build_ladder({"T_start": 0.01, "T_stop": 0.001, "num": 4, "spacing": "log"})
    assert out[0] == pytest.approx(0.01)
    assert out[-1] == pytest.approx(0.001)
    assert len(out) == 4


def test_neither_temperatures_nor_generator_raises():
    with pytest.raises(ValueError):
        ts.build_ladder({})


def test_non_descending_warns_but_returns(caplog):
    with caplog.at_level(logging.WARNING):
        out = ts.build_ladder({"temperatures": [0.002, 0.004], "warm_start": True})
    assert out == [0.002, 0.004]
    assert any("descend" in r.message.lower() for r in caplog.records)


def test_resolve_sigma_name_default_and_custom():
    assert ts.resolve_sigma_name({"file": {"output": {}}}) == "sigma.npz"
    assert ts.resolve_sigma_name({"file": {"output": {"sigma": "sig"}}}) == "sig.npz"
    assert ts.resolve_sigma_name({"file": {"output": {"sigma": "s.npz"}}}) == "s.npz"


def test_resolve_gap_name_dynamic_vs_static():
    assert ts.resolve_gap_name({"eliashberg": {"frequency": "dynamic"}}) == "gap_dynamic.npz"
    assert ts.resolve_gap_name({"eliashberg": {}}) == "gap.npz"
    assert ts.resolve_gap_name({}) == "gap.npz"


def test_rung_dir_format():
    assert ts.rung_dir("/out", 0, 0.01).endswith("000_T0.01")
    assert ts.rung_dir("/out", 12, 0.002).endswith("012_T0.002")


def _valid_base():
    return {"mode": {"param": {"CellShape": [2, 2, 1], "Nmat": 16, "filling": 0.5}},
            "file": {"input": {}, "output": {}},
            "eliashberg": {"frequency": "dynamic"}}


def test_preflight_ok():
    ts.preflight(_valid_base(), {"run_eliashberg": True})


@pytest.mark.parametrize("drop", ["CellShape", "Nmat"])
def test_preflight_missing_shape_field(drop):
    base = _valid_base()
    del base["mode"]["param"][drop]
    with pytest.raises(ValueError, match=drop):
        ts.preflight(base, {"run_eliashberg": True})


def test_preflight_missing_filling_and_ncond():
    base = _valid_base()
    del base["mode"]["param"]["filling"]
    with pytest.raises(ValueError, match="filling"):
        ts.preflight(base, {"run_eliashberg": True})


def test_preflight_run_eliashberg_without_section():
    base = _valid_base()
    del base["eliashberg"]
    with pytest.raises(ValueError, match="eliashberg"):
        ts.preflight(base, {"run_eliashberg": True})
    # FLEX-only is fine without [eliashberg]:
    ts.preflight(base, {"run_eliashberg": False})
