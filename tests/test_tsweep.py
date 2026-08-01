import copy
import json
import logging
import os
import numpy as np
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


def test_generated_single_point_ladder():
    # num == 1 is a distinct branch: yields [T_start], ignoring T_stop/spacing.
    assert ts.build_ladder(
        {"T_start": 0.01, "T_stop": 0.002, "num": 1, "spacing": "log"}) == [0.01]


def test_invalid_spacing_raises():
    with pytest.raises(ValueError):
        ts.build_ladder({"T_start": 0.01, "T_stop": 0.002, "num": 3, "spacing": "cubic"})


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
            "eliashberg": {"frequency": "dynamic", "solver_mode": "eigenvalue"}}


def test_preflight_ok():
    ts.preflight(_valid_base(), {"run_eliashberg": True})


def test_preflight_requires_eigenvalue_solver_mode():
    # default solver_mode "iteration" writes no parity/index-0 block -> reject
    base = _valid_base()
    base["eliashberg"]["solver_mode"] = "iteration"
    with pytest.raises(ValueError, match="solver_mode"):
        ts.preflight(base, {"run_eliashberg": True})
    del base["eliashberg"]["solver_mode"]            # absent -> defaults iteration
    with pytest.raises(ValueError, match="solver_mode"):
        ts.preflight(base, {"run_eliashberg": True})
    # FLEX-only ladders are unaffected
    ts.preflight(base, {"run_eliashberg": False})


def test_documented_example_passes_preflight():
    # Mirror the complete [eliashberg]/[continuation] example in
    # docs/{en,ja}/source/rpa/tutorial/sc-index.rst so a copy-paste of the
    # tutorial does not fail hwave_tsweep preflight (issue #61).
    base = {
        "mode": {"param": {"T": 0.02, "CellShape": [32, 32, 1],
                           "Nmat": 512, "filling": 0.75}},
        "file": {"input": {}, "output": {"path_to_output": "output"}},
        "eliashberg": {"frequency": "dynamic", "chi0q_mode": "flex",
                       "pairing_type": "singlet",
                       "solver_mode": "eigenvalue"},
    }
    cont = {"T_start": 0.02, "T_stop": 0.005, "num": 6, "spacing": "log",
            "run_eliashberg": True, "warm_start": True, "seed_gap": True}
    ts.preflight(base, cont)


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


def test_make_rung_dicts_does_not_mutate_base():
    base = _valid_base()
    base["file"]["input"]["sigma_init"] = "stale.npz"
    base["eliashberg"]["seed_eigenvector"] = "stale_gap.npz"
    base["continuation"] = {"temperatures": [0.01]}
    base["mode"]["param"]["T"] = 0.5
    ts.make_rung_dicts(base, 0.008, "/o/r", run_eliashberg=True,
                       sigma_init="S", seed_gap="G")
    # caller's base must be untouched (function operates on deep copies only)
    assert base["file"]["input"]["sigma_init"] == "stale.npz"
    assert base["eliashberg"]["seed_eigenvector"] == "stale_gap.npz"
    assert base["continuation"] == {"temperatures": [0.01]}
    assert base["mode"]["param"]["T"] == 0.5


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


def _cont_base(tmp_path, temps, run_eli=True):
    base = _valid_base()
    base["continuation"] = {"temperatures": temps, "output_dir": "sweep",
                            "run_eliashberg": run_eli, "warm_start": True,
                            "seed_gap": True, "summary_file": "lambda_vs_T.dat"}
    return base


def _install_fake_solvers(monkeypatch, tmp_path, fail_at=None):
    """Fake qlms.run: writes a sigma.npz, returns converged dict, records calls.
    Fake sc.calc_eliashberg: writes eigenvalue.dat + gap_dynamic.npz."""
    calls = {"flex": [], "eli": []}

    def fake_flex(*, input_dict):
        calls["flex"].append(copy.deepcopy(input_dict))
        out = input_dict["file"]["output"]["path_to_output"]
        os.makedirs(out, exist_ok=True)
        np.savez(os.path.join(out, "sigma.npz"), sigma=np.ones((1,)))
        return {"scf_converged": True, "scf_iterations": 10}

    def fake_eli(input_dict):
        calls["eli"].append(copy.deepcopy(input_dict))
        T = input_dict["mode"]["param"]["T"]
        if fail_at is not None and abs(T - fail_at) < 1e-12:
            raise RuntimeError("boom")
        out = input_dict["file"]["output"]["path_to_output"]
        os.makedirs(out, exist_ok=True)
        np.savez(os.path.join(out, "gap_dynamic.npz"), gap=np.ones((1,)))
        with open(os.path.join(out, "eigenvalue.dat"), "w") as fw:
            fw.write("# index Re Im |ev| match\n0 0.5 0.0 0.5 1\n")

    monkeypatch.setattr("hwave.qlms.run", fake_flex)
    monkeypatch.setattr("hwave.sc.calc_eliashberg", fake_eli)
    return calls


def test_run_chains_seeds(monkeypatch, tmp_path):
    calls = _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    rows = ts.run(base, base_dir=str(tmp_path))
    assert [r["status"] for r in rows] == ["ok", "ok", "ok"]
    # rung 0 has no seeds; rung 1,2 point at previous rung outputs, absolute
    assert "sigma_init" not in calls["flex"][0]["file"]["input"]
    si1 = calls["flex"][1]["file"]["input"]["sigma_init"]
    assert os.path.isabs(si1) and "000_T0.01" in si1 and si1.endswith("sigma.npz")
    seed1 = calls["eli"][1]["eliashberg"]["seed_eigenvector"]
    assert os.path.isabs(seed1) and "000_T0.01" in seed1 and seed1.endswith("gap_dynamic.npz")
    # summary file written under absolute output_dir
    assert os.path.exists(os.path.join(tmp_path, "sweep", "lambda_vs_T.dat"))


def test_run_error_stops_chain(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path, fail_at=0.008)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    rows = ts.run(base, base_dir=str(tmp_path))
    assert [r["status"] for r in rows] == ["ok", "error"]           # stopped
    assert rows[1]["error_stage"] == "eliashberg"


def test_run_keep_going_reseeds(monkeypatch, tmp_path):
    calls = _install_fake_solvers(monkeypatch, tmp_path, fail_at=0.008)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    rows = ts.run(base, base_dir=str(tmp_path), keep_going=True)
    assert [r["status"] for r in rows] == ["ok", "error", "ok"]
    # rung 2 (0.006) cold-started (no seed: prev rung failed)
    assert "sigma_init" not in calls["flex"][2]["file"]["input"]


def test_run_dry_run_invokes_no_solver(monkeypatch, tmp_path):
    calls = _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    rows = ts.run(base, base_dir=str(tmp_path), dry_run=True)
    assert calls["flex"] == [] and calls["eli"] == []
    assert all(r["status"] == "dry" for r in rows)


def test_run_injects_sigma_output_for_warm_start(monkeypatch, tmp_path):
    calls = _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])       # _valid_base has no file.output.sigma
    assert "sigma" not in base["file"]["output"]
    ts.run(base, base_dir=str(tmp_path))
    # every FLEX dict must carry a sigma output name so the self-energy is written
    assert all(f["file"]["output"].get("sigma") == "sigma" for f in calls["flex"])


def test_run_respects_custom_sigma_output(monkeypatch, tmp_path):
    calls = _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    base["file"]["output"]["sigma"] = "selfE"        # user-set name is preserved
    ts.run(base, base_dir=str(tmp_path))
    assert all(f["file"]["output"]["sigma"] == "selfE" for f in calls["flex"])


def test_run_does_not_mutate_caller_input(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    base["file"]["input"]["path_to_input"] = "rel/inputs"   # relative -> resolved inside
    snapshot = copy.deepcopy(base)
    ts.run(base, base_dir=str(tmp_path))
    assert base == snapshot        # run() operates on its own deep copy


def test_main_dry_run(monkeypatch, tmp_path, capsys):
    # tomli_w is not a project dependency; write the TOML as plain text
    # instead of relying on it (see plan Task 8 fallback).
    toml_path = tmp_path / "in.toml"
    toml_path.write_text(
        '[mode.param]\nCellShape=[2,2,1]\nNmat=16\nfilling=0.5\n'
        '[file.input]\n[file.output]\n'
        '[continuation]\ntemperatures=[0.01,0.008]\noutput_dir="sweep"\n'
        'run_eliashberg=false\n')
    monkeypatch.setattr("sys.argv", ["hwave_tsweep", str(toml_path), "--dry-run"])
    ts.main()
    # dry-run writes a summary of planned rungs, invokes no solver
    assert (tmp_path / "sweep" / "lambda_vs_T.dat").exists()


# --- resume / restart (issue #64) -------------------------------------------

def _rout(tmp_path, idx, T):
    return os.path.join(tmp_path, "sweep", ts.rung_dir("", idx, T).lstrip("/"),
                        "output")


def test_resume_after_eliashberg_complete_sweep_is_noop(monkeypatch, tmp_path):
    """A fully completed sweep resumed reruns NO solver and reproduces every
    rung from the recorded state."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)     # fresh call log
    rows = ts.run(base, base_dir=str(tmp_path), resume=True)
    assert calls2["flex"] == [] and calls2["eli"] == []
    assert [r["status"] for r in rows] == ["ok", "ok", "ok"]


def test_resume_after_flex_interruption_reseeds(monkeypatch, tmp_path):
    """Interrupted mid-sweep (rung 1's Eliashberg failed): resume reuses the
    completed rung 0 and restarts at rung 1, seeded from rung 0."""
    _install_fake_solvers(monkeypatch, tmp_path, fail_at=0.008)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    rows1 = ts.run(base, base_dir=str(tmp_path))
    assert [r["status"] for r in rows1] == ["ok", "error"]

    calls2 = _install_fake_solvers(monkeypatch, tmp_path)     # healthy now
    rows2 = ts.run(base, base_dir=str(tmp_path), resume=True)
    # rung 0 reused (not recomputed): first solver call is rung 1
    assert calls2["flex"], "resume ran no rungs"
    assert calls2["flex"][0]["mode"]["param"]["T"] == 0.008
    # ... and it is seeded from rung 0's outputs
    si = calls2["flex"][0]["file"]["input"].get("sigma_init")
    assert si and "000_T0.01" in si and si.endswith("sigma.npz")
    seed = calls2["eli"][0]["eliashberg"].get("seed_eigenvector")
    assert seed and "000_T0.01" in seed
    assert [r["status"] for r in rows2] == ["ok", "ok", "ok"]


def test_resume_detects_corrupt_eigenvalue(monkeypatch, tmp_path):
    """A rung whose eigenvalue.dat is unparseable is NOT reused even if the
    summary marked it ok; it (and everything after) is recomputed."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    # corrupt rung 1's eigenvalue.dat
    eig = os.path.join(_rout(tmp_path, 1, 0.008), "eigenvalue.dat")
    with open(eig, "w") as fw:
        fw.write("# garbage\nnot a number here\n")
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    ts.run(base, base_dir=str(tmp_path), resume=True)
    assert [c["mode"]["param"]["T"] for c in calls2["flex"]] == [0.008, 0.006]


def test_resume_detects_missing_output(monkeypatch, tmp_path):
    """A rung missing its eigenvalue.dat on disk is recomputed (file existence
    is checked, not only the summary)."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    os.remove(os.path.join(_rout(tmp_path, 1, 0.008), "eigenvalue.dat"))
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    ts.run(base, base_dir=str(tmp_path), resume=True)
    assert [c["mode"]["param"]["T"] for c in calls2["flex"]] == [0.008, 0.006]


@pytest.mark.parametrize("name", ["sigma.npz", "gap_dynamic.npz"])
def test_resume_detects_corrupt_seed(monkeypatch, tmp_path, name):
    """A present but unreadable seed forces its producing rung to rerun."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    with open(os.path.join(_rout(tmp_path, 1, 0.008), name), "wb") as fw:
        fw.write(b"not an npz")
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    ts.run(base, base_dir=str(tmp_path), resume=True)
    assert [c["mode"]["param"]["T"] for c in calls2["flex"]] == [0.008, 0.006]


def test_resume_config_mismatch_fails_fast(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    ts.run(base, base_dir=str(tmp_path))
    base2 = copy.deepcopy(base)
    base2["mode"]["param"]["Nmat"] = 999                     # shape change
    with pytest.raises(ValueError, match="configuration mismatch"):
        ts.run(base2, base_dir=str(tmp_path), resume=True)


def test_resume_ladder_mismatch_fails_fast(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    ts.run(base, base_dir=str(tmp_path))
    base2 = _cont_base(tmp_path, [0.01, 0.008, 0.006])       # different ladder
    with pytest.raises(ValueError, match="ladder mismatch"):
        ts.run(base2, base_dir=str(tmp_path), resume=True)


def test_resume_without_manifest_fails_fast(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    with pytest.raises(ValueError, match="manifest"):
        ts.run(base, base_dir=str(tmp_path), resume=True)


def test_resume_via_continuation_flag(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    base["continuation"]["resume"] = True                    # config-driven resume
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    rows = ts.run(base, base_dir=str(tmp_path))
    assert calls2["flex"] == [] and calls2["eli"] == []
    assert [r["status"] for r in rows] == ["ok", "ok", "ok"]


def test_fresh_run_writes_manifest(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    ts.run(base, base_dir=str(tmp_path))
    assert os.path.exists(os.path.join(tmp_path, "sweep", ts.MANIFEST_NAME))


def test_fresh_run_invalidates_old_summary_before_new_manifest(monkeypatch,
                                                               tmp_path):
    """An interrupted fresh run cannot pair old rows with its new manifest."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    ts.run(base, base_dir=str(tmp_path))

    original = ts.write_manifest

    def interrupt_after_checkpoint_reset(*args, **kwargs):
        original(*args, **kwargs)
        raise KeyboardInterrupt

    monkeypatch.setattr(ts, "write_manifest", interrupt_after_checkpoint_reset)
    with pytest.raises(KeyboardInterrupt):
        ts.run(base, base_dir=str(tmp_path))
    assert ts.read_summary_rows(
        os.path.join(tmp_path, "sweep", "lambda_vs_T.dat")) == []


def test_resume_rejects_unknown_manifest_version(monkeypatch, tmp_path):
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008])
    ts.run(base, base_dir=str(tmp_path))
    manifest = os.path.join(tmp_path, "sweep", ts.MANIFEST_NAME)
    with open(manifest) as f:
        data = json.load(f)
    data["version"] += 1
    with open(manifest, "w") as f:
        json.dump(data, f)
    with pytest.raises(ValueError, match="version mismatch"):
        ts.run(base, base_dir=str(tmp_path), resume=True)


def test_summary_roundtrip_and_atomic(tmp_path):
    p = str(tmp_path / "lambda_vs_T.dat")
    rows = [
        {"idx": 0, "T": 0.01, "status": "ok", "error_stage": "none",
         "re": 0.36, "im": 0.0, "match": 1, "flex_converged": 1, "flex_iter": 14},
        {"idx": 1, "T": 0.008, "status": "not_converged", "error_stage": "none",
         "re": 0.41, "im": 0.0, "match": 1, "flex_converged": 0, "flex_iter": 50},
    ]
    ts.write_summary(p, rows)
    # atomic: no staging file (`.tmp.<pid>`) left behind after a clean write
    import glob
    assert glob.glob(p + ".tmp*") == []
    back = ts.read_summary_rows(p)
    assert [r["idx"] for r in back] == [0, 1]
    assert back[0]["status"] == "ok" and back[1]["status"] == "not_converged"
    assert abs(back[0]["re"] - 0.36) < 1e-9
    assert back[1]["flex_converged"] == 0 and back[1]["flex_iter"] == 50


def test_fingerprint_ignores_temperature_but_tracks_shape():
    base = _valid_base()
    base["mode"]["param"]["T"] = 0.01
    fp_a = ts.config_fingerprint(base, run_eli=True)
    base["mode"]["param"]["T"] = 0.5                         # T must not matter
    fp_b = ts.config_fingerprint(base, run_eli=True)
    assert fp_a == fp_b
    base["mode"]["param"]["Nmat"] = 999                      # shape must matter
    assert ts.config_fingerprint(base, run_eli=True) != fp_a


# --- review follow-ups: fingerprint completeness, final-rung validation,
#     NaN ladder, strict summary parsing (Codex Phase-1 on #64) --------------

def test_atomic_write_failed_replace_preserves_old(tmp_path, monkeypatch):
    """If os.replace fails mid-update, the previous checkpoint must survive
    intact and no staging file should be observable to a reader."""
    import glob
    p = str(tmp_path / "s.dat")
    with open(p, "w") as fw:
        fw.write("OLD")
    monkeypatch.setattr(ts.os, "replace",
                        lambda *a, **k: (_ for _ in ()).throw(OSError("boom")))
    with pytest.raises(OSError):
        ts._atomic_write_text(p, "NEW")
    assert open(p).read() == "OLD"                          # old content intact
    # (the leaked .tmp.<pid> is the caller's to clean up; the target is safe)


@pytest.mark.parametrize("key,val", [("Coulomb", "coulomb.dat"),
                                     ("Extern", "extern.dat")])
def test_config_fingerprint_tracks_all_supported_interaction_files(
        tmp_path, key, val):
    """`Coulomb` (combined) and `Extern` are accepted by the k-space reader and
    must be part of the fingerprint (filename + content)."""
    (tmp_path / val).write_text("1 0 0\n")
    base = {"mode": {"mode": "FLEX",
                     "param": {"CellShape": [2, 2, 1], "Nmat": 16, "filling": 0.5}},
            "file": {"input": {"interaction": {"path_to_input": ".", key: val}}},
            "eliashberg": {"frequency": "dynamic", "solver_mode": "eigenvalue"}}
    fp1 = ts.config_fingerprint(base, run_eli=True, base_dir=str(tmp_path))
    (tmp_path / val).write_text("2 0 0\n")                  # in-place edit
    assert ts.config_fingerprint(base, run_eli=True,
                                 base_dir=str(tmp_path)) != fp1


def test_config_fingerprint_tracks_interaction_content(tmp_path):
    """Editing an interaction file in place (same filename) must change the
    fingerprint, so resume cannot silently reuse rungs from a different system."""
    (tmp_path / "geom.dat").write_text("1 0 0\n")
    base = {"mode": {"mode": "FLEX",
                     "param": {"CellShape": [2, 2, 1], "Nmat": 16, "filling": 0.5}},
            "file": {"input": {"interaction": {"path_to_input": ".",
                                               "Geometry": "geom.dat"}}},
            "eliashberg": {"frequency": "dynamic", "solver_mode": "eigenvalue"}}
    fp1 = ts.config_fingerprint(base, run_eli=True, base_dir=str(tmp_path))
    (tmp_path / "geom.dat").write_text("2 0 0\n")            # in-place edit
    fp2 = ts.config_fingerprint(base, run_eli=True, base_dir=str(tmp_path))
    assert fp1 != fp2


def test_config_fingerprint_tracks_all_solver_params(tmp_path):
    """The fingerprint must cover the full mode.param and eliashberg blocks,
    not just a hand-picked subset."""
    base = _valid_base()
    fp = ts.config_fingerprint(base, run_eli=True, base_dir=str(tmp_path))
    b2 = copy.deepcopy(base)
    b2["eliashberg"]["pairing_type"] = "triplet"            # physics param
    assert ts.config_fingerprint(b2, run_eli=True, base_dir=str(tmp_path)) != fp
    b3 = copy.deepcopy(base)
    b3["mode"]["param"]["Mix"] = 0.9                        # any mode param
    assert ts.config_fingerprint(b3, run_eli=True, base_dir=str(tmp_path)) != fp


def test_resume_final_rung_corrupt_output_recomputed(monkeypatch, tmp_path):
    """A corrupt output on the FINAL ladder rung must be detected and that rung
    recomputed -- the final rung is no longer exempt from output validation."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006])
    ts.run(base, base_dir=str(tmp_path))
    with open(os.path.join(_rout(tmp_path, 2, 0.006), "sigma.npz"), "w") as fw:
        fw.write("garbage not an npz")
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    ts.run(base, base_dir=str(tmp_path), resume=True)
    assert [c["mode"]["param"]["T"] for c in calls2["flex"]] == [0.006]


def test_resume_flex_only_with_dynamic_block_is_noop(monkeypatch, tmp_path):
    """A FLEX-only sweep (run_eliashberg=false) that still carries a dynamic
    [eliashberg] block must not demand a gap file: it seeds only sigma, and a
    completed sweep resumes with no solver calls."""
    _install_fake_solvers(monkeypatch, tmp_path)
    base = _cont_base(tmp_path, [0.01, 0.008, 0.006], run_eli=False)
    base["eliashberg"] = {"frequency": "dynamic", "solver_mode": "eigenvalue"}
    ts.run(base, base_dir=str(tmp_path))
    calls2 = _install_fake_solvers(monkeypatch, tmp_path)
    rows = ts.run(base, base_dir=str(tmp_path), resume=True)
    assert calls2["flex"] == [] and calls2["eli"] == []
    assert [r["status"] for r in rows] == ["ok", "ok", "ok"]


def test_validate_resume_rejects_nan_ladder():
    """A NaN manifest temperature must fail fast, not slip through the tolerance
    comparison (NaN comparisons are always False)."""
    manifest = {"version": ts._MANIFEST_VERSION,
                "ladder": [0.01, float("nan")], "fingerprint": "x"}
    with pytest.raises(ValueError, match="ladder mismatch"):
        ts._validate_resume(manifest, [0.01, 0.008], "x")


def test_read_summary_rows_truncates_on_out_of_order(tmp_path):
    """A duplicate/out-of-order (or malformed) row truncates the parsed prefix,
    so a damaged checkpoint shortens the reusable prefix instead of being
    silently patched over."""
    p = str(tmp_path / "s.dat")
    with open(p, "w") as fw:
        fw.write(ts._SUMMARY_HEADER)
        fw.write("0 0.01 ok none 0.5 0.0 1 1 10\n")
        fw.write("2 0.006 ok none 0.5 0.0 1 1 10\n")        # expected idx 1
    assert [r["idx"] for r in ts.read_summary_rows(p)] == [0]


def test_read_summary_rows_truncates_on_malformed(tmp_path):
    p = str(tmp_path / "s.dat")
    with open(p, "w") as fw:
        fw.write(ts._SUMMARY_HEADER)
        fw.write("0 0.01 ok none 0.5 0.0 1 1 10\n")
        fw.write("1 notanumber ok none 0.5 0.0 1 1 10\n")   # malformed T
    assert [r["idx"] for r in ts.read_summary_rows(p)] == [0]


import tempfile
import unittest


class TestFingerprintCaseInsensitivity(unittest.TestCase):
    """Third surfacing of the case-sensitivity defect class (PR #128
    round 6): a lowercase interaction key EXECUTED correctly but was
    omitted from the resume fingerprint, so a resume could silently
    reuse results computed with a different interaction file. unittest
    TestCase on purpose -- the module's pytest-style functions do not
    gate in CI."""

    def _base(self, d, key, fname):
        return {
            "mode": {"mode": "RPA", "param": {"T": 1.0,
                                              "CellShape": [2, 2, 1]}},
            "file": {"input": {"interaction": {"path_to_input": d,
                                               key: fname}}},
        }

    def test_every_interaction_key_participates_lowercase(self):
        import hwave.tsweep as ts

        for key in ts._INTERACTION_KEYS:
            with self.subTest(key=key):
                d = tempfile.mkdtemp()
                fn = os.path.join(d, "x.dat")
                with open(fn, "w") as f:
                    f.write("one\n")
                fp_lower = ts.config_fingerprint(
                    self._base(d, key.lower(), "x.dat"), False, d)
                fp_canon = ts.config_fingerprint(
                    self._base(d, key, "x.dat"), False, d)
                # canonical and lowercase configurations are the same run
                self.assertEqual(fp_lower, fp_canon)
                # an in-place content edit must invalidate the resume,
                # regardless of the key's case
                with open(fn, "w") as f:
                    f.write("two\n")
                fp_edited = ts.config_fingerprint(
                    self._base(d, key.lower(), "x.dat"), False, d)
                self.assertNotEqual(fp_lower, fp_edited)


class TestReaderKeyCaseParity(unittest.TestCase):
    """QLMSkInput accepted 'geometry' in validation but dispatched on the
    exact-case key, feeding the geometry file to the wannier90
    interaction parser (loud ValueError). Dispatch is normalized now."""

    def test_lowercase_geometry_and_transfer_load(self):
        import hwave.qlmsio.read_input_k as read_input_k

        d = tempfile.mkdtemp()
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n"
                    "1\n 0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n1\n1\n 1\n"
                    "   1    0    0    1    1   1.000000   0.0\n")
        canon = read_input_k.QLMSkInput(
            {"path_to_input": d,
             "interaction": {"path_to_input": d, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat"}}).get_param("ham")
        lower = read_input_k.QLMSkInput(
            {"path_to_input": d,
             "interaction": {"path_to_input": d, "geometry": "geom.dat",
                             "transfer": "transfer.dat"}}).get_param("ham")
        self.assertEqual(lower["Geometry"]["norb"],
                         canon["Geometry"]["norb"])
        self.assertEqual(sorted(lower["Transfer"].keys()),
                         sorted(canon["Transfer"].keys()))
