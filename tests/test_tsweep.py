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


def test_make_rung_dicts_first_rung_no_seeds():
    base = _valid_base()
    base["file"]["input"]["sigma_init"] = "stale.npz"      # must be stripped
    base["eliashberg"]["seed_eigenvector"] = "stale_gap.npz"
    base["continuation"] = {"temperatures": [0.01]}         # must be dropped
    flex, eli = ts.make_rung_dicts(base, 0.01, "/o/000_T0.01", run_eliashberg=True)
    assert flex["mode"]["param"]["T"] == 0.01
    assert flex["file"]["output"]["path_to_output"] == "/o/000_T0.01"
    assert flex["file"]["input"]["path_to_flex_output"] == "/o/000_T0.01"
    assert "sigma_init" not in flex["file"]["input"]          # stripped, no seed
    assert "seed_eigenvector" not in eli["eliashberg"]
    assert "continuation" not in flex and "continuation" not in eli


def test_make_rung_dicts_seeded_rung():
    base = _valid_base()
    flex, eli = ts.make_rung_dicts(base, 0.008, "/o/001_T0.008", run_eliashberg=True,
                                   sigma_init="/o/000_T0.01/output/sigma.npz",
                                   seed_gap="/o/000_T0.01/output/gap_dynamic.npz")
    assert flex["file"]["input"]["sigma_init"] == "/o/000_T0.01/output/sigma.npz"
    assert eli["eliashberg"]["seed_eigenvector"] == "/o/000_T0.01/output/gap_dynamic.npz"


def test_make_rung_dicts_flex_only():
    base = _valid_base()
    del base["eliashberg"]
    flex, eli = ts.make_rung_dicts(base, 0.01, "/o/000_T0.01", run_eliashberg=False)
    assert eli is None
    assert "path_to_flex_output" not in flex["file"]["input"]


def test_make_rung_dicts_invariance_except_seed_keys():
    import copy
    base = _valid_base()
    flex, eli = ts.make_rung_dicts(base, 0.008, "/o/r", run_eliashberg=True,
                                   sigma_init="S", seed_gap="G")
    # canonical = flex without its sigma_init:
    canon = copy.deepcopy(flex)
    del canon["file"]["input"]["sigma_init"]
    eli_canon = copy.deepcopy(eli)
    del eli_canon["eliashberg"]["seed_eigenvector"]
    assert canon == eli_canon      # differ only by their own seed keys


def test_parse_leading_eig_with_match(tmp_path):
    p = tmp_path / "eigenvalue.dat"
    p.write_text("# Eigenvalue analysis\n"
                 "# index  Re  Im  |ev|  match\n"
                 "   0  6.63000000e-01  0.00000000e+00  6.63e-01 1\n"
                 "   1  3.20000000e-01  0.00000000e+00  3.20e-01 0\n")
    re, im, match = ts.parse_leading_eig(str(p))
    assert re == pytest.approx(0.663)
    assert im == pytest.approx(0.0)
    assert match == 1


def test_parse_leading_eig_no_match_column(tmp_path):
    p = tmp_path / "eigenvalue.dat"
    p.write_text("# index Re Im |ev|\n   0  4.00e-01  1.00e-02  4.0e-01\n")
    re, im, match = ts.parse_leading_eig(str(p))
    assert re == pytest.approx(0.4) and im == pytest.approx(0.01) and match == -1


def test_parse_leading_eig_missing_raises(tmp_path):
    with pytest.raises(ValueError):
        ts.parse_leading_eig(str(tmp_path / "nope.dat"))


def test_write_summary_schema(tmp_path):
    p = tmp_path / "lambda_vs_T.dat"
    ts.write_summary(str(p), [
        {"idx": 0, "T": 0.01, "status": "ok", "error_stage": "none",
         "re": 0.361, "im": 0.0, "match": 1, "flex_converged": 1, "flex_iter": 14},
        {"idx": 1, "T": 0.006, "status": "error", "error_stage": "eliashberg",
         "re": float("nan"), "im": float("nan"), "match": -1,
         "flex_converged": 1, "flex_iter": 12},
    ])
    lines = p.read_text().splitlines()
    assert lines[0].startswith("#")
    assert "0.361000" in lines[1] and "ok" in lines[1] and "none" in lines[1]
    assert "nan" in lines[2] and "error" in lines[2] and "eliashberg" in lines[2]
