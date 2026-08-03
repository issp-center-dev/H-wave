"""Fail-fast guards for IR-native FLEX files (Stage 3, design
docs/design/ir-matsubara-stage3.md Sec. 4/4.3).

Uniform-grid readers must reject sparse-node files BEFORE any permissive
legacy fallback (e.g. the missing-freq_index center-row path) can silently
misread them.

Import-safe under both pytest and unittest discovery. Tests run from the
repository root.
"""
import os

import numpy as np
import pytest

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

pytestmark = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                                reason="sparse-ir not installed")


def _make_solver(nmat=64, extra_param=None, iteration_max=3):
    info_mode = {'mode': 'FLEX', 'param': {
        'T': 0.5, 'mu': 0.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
        'Nmat': nmat, 'IterationMax': iteration_max, 'Mix': 0.3, 'EPS': 8}}
    info_mode['calc_scheme'] = 'reduced'
    if extra_param:
        info_mode['param'].update(extra_param)
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'}}
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    rio = read_input_k.QLMSkInput(info_input)
    ham = rio.get_param("ham")
    gi = rio.get_param("green")
    ham["CoulombIntra"] = {k: complex(2.5) for k in ham["CoulombIntra"]}
    return solver_flex.FLEX(ham, {}, info_mode), gi


def _write_run(outdir, native):
    extra = {'matsubara_basis': 'ir'}
    if native:
        extra['write_densified'] = False
    solver, gi = _make_solver(extra_param=extra)
    solver.solve(gi, outdir)
    solver.save_results({"path_to_output": outdir, "chi0q": "chi0q",
                         "sigma": "sigma", "green": "green"}, gi)
    return solver, gi


def test_static_chi0q_loader_rejects_ir_native(tmp_path):
    import hwave.sc as sc
    out = str(tmp_path)
    _write_run(out, native=True)
    input_dict = {"file": {"output": {"path_to_output": out,
                                      "chi0q": "chi0q"}},
                  "mode": {}}
    with pytest.raises(ValueError, match="write_densified"):
        sc._load_chi0q(input_dict)


def test_static_flex_chi_loader_rejects_ir_native(tmp_path):
    import hwave.sc as sc
    out = str(tmp_path)
    _write_run(out, native=True)
    input_dict = {"file": {"input": {"path_to_flex_output": out}},
                  "eliashberg": {}}
    with pytest.raises(ValueError, match="write_densified"):
        sc._read_flex_chi_raw(input_dict)


def test_rpa_read_chi0q_rejects_ir_native(tmp_path):
    out = str(tmp_path)
    solver, _ = _write_run(out, native=True)
    with pytest.raises(ValueError, match="write_densified"):
        solver._read_chi0q(os.path.join(out, "chi0q.npz"))


def test_ir_native_sigma_init_into_uniform_run_rejected(tmp_path):
    out = str(tmp_path)
    _write_run(out, native=True)
    uniform_solver, _ = _make_solver()
    with pytest.raises(ValueError, match="matsubara_basis"):
        uniform_solver.read_init({"path_to_input": out,
                                  "sigma_init": "sigma.npz"})


def test_densified_files_still_load_everywhere(tmp_path):
    """Migration pin: densified IR files (and thus pre-Stage-3 files, which
    share the schema) pass every guard unchanged."""
    import hwave.sc as sc
    out = str(tmp_path)
    _write_run(out, native=False)
    input_dict = {"file": {"output": {"path_to_output": out,
                                      "chi0q": "chi0q"}},
                  "mode": {}}
    chi0q, static_index = sc._load_chi0q(input_dict)
    assert static_index == 64 // 2
    input_dict2 = {"file": {"input": {"path_to_flex_output": out}},
                   "eliashberg": {}}
    chi_s, chi_c, conv, _tags = sc._read_flex_chi_raw(input_dict2)
    assert conv == "kuroki"
    assert chi_s.shape[0] == 64


# ---------------------------------------------------------------------------
# Dynamic Eliashberg ingestion of IR-native files (design Sec. 4.1)
# ---------------------------------------------------------------------------

def _dyn_input(outdir, nmat=64, extra=None):
    eli = {"chi0q_mode": "flex", "frequency": "dynamic",
           "solver_mode": "iteration", "max_iter": 200,
           "convergence_tol": 1e-7, "matsubara_basis": "ir"}
    if extra:
        eli.update(extra)
    return {
        "mode": {"param": {"T": 0.5, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": nmat,
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": "tests/rpa/input",
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": outdir}},
        "eliashberg": eli,
    }


def _write_converged_run(outdir, native):
    extra = {'matsubara_basis': 'ir'}
    if native:
        extra['write_densified'] = False
    solver, gi = _make_solver(extra_param=extra, iteration_max=60)
    solver.solve(gi, outdir)
    assert solver.scf_converged
    solver.save_results({"path_to_output": outdir, "chi0q": "chi0q",
                         "sigma": "sigma", "green": "green",
                         "chiq_s": "chiq_s", "chiq_c": "chiq_c"}, gi)
    return solver, gi


def test_dynamic_uniform_mode_rejects_ir_native(tmp_path):
    from hwave.solver import eliashberg_dynamic as ed
    out = str(tmp_path)
    _write_converged_run(out, native=True)
    inp = _dyn_input(out, extra={"matsubara_basis": "uniform"})
    with pytest.raises(ValueError, match="matsubara_basis"):
        ed.solve_dynamic(inp)


def test_dynamic_mixed_encoding_rejected(tmp_path):
    import shutil
    from hwave.solver import eliashberg_dynamic as ed
    a = tmp_path / "native"
    b = tmp_path / "dens"
    a.mkdir(); b.mkdir()
    _write_converged_run(str(a), native=True)
    _write_converged_run(str(b), native=False)
    shutil.copy(str(b / "chiq_c.npz"), str(a / "chiq_c.npz"))
    with pytest.raises(ValueError, match="one encoding"):
        ed.solve_dynamic(_dyn_input(str(a)))


def test_dynamic_ir_native_chain_matches_densified_chain(tmp_path):
    """FLEX(ir) -> dynamic(ir) lambda: the IR-native file chain must agree
    with the densified chain built from THE SAME converged run (offline
    densify), at eps scale (the refit replaces the full-grid compress)."""
    from hwave.solver import eliashberg_dynamic as ed
    nat = tmp_path / "nat"
    den = tmp_path / "den"
    nat.mkdir(); den.mkdir()
    solver, gi = _write_converged_run(str(nat), native=True)
    axF, axB = solver._ir_axF, solver._ir_axB
    nmat = solver.nmat
    # offline densify of the SAME arrays -> the reference densified files
    gi_d = dict(gi)
    gi_d["sigma"] = np.moveaxis(axF.eval_to_uniform(
        axF.fit_from_freq(np.moveaxis(gi["sigma"], 1, -1)), nmat), -1, 1)
    gi_d["green"] = np.moveaxis(axF.eval_to_uniform(
        axF.fit_from_freq(np.moveaxis(gi["green"], 1, -1)), nmat), -1, 1)
    for key in ("chi0q", "chiq_s", "chiq_c"):
        gi_d[key] = np.moveaxis(axB.eval_to_uniform(
            axB.fit_from_freq(np.moveaxis(gi[key], 0, -1)), nmat), -1, 0)
    solver.write_densified = True
    solver.save_results({"path_to_output": str(den), "chi0q": "chi0q",
                         "sigma": "sigma", "green": "green",
                         "chiq_s": "chiq_s", "chiq_c": "chiq_c"}, gi_d)

    lam_n = ed.solve_dynamic(_dyn_input(str(nat)))
    lam_d = ed.solve_dynamic(_dyn_input(str(den)))
    assert abs(lam_n - lam_d) < 1e-4 * abs(lam_d), (lam_n, lam_d)


def test_refit_nodes_pass_through_and_beta_guard():
    """Unit contract of the ingestion helper: exact node-set equality is a
    pass-through (values untouched); a beta mismatch on chi/green input is
    a hard error printing both values."""
    from hwave.solver import eliashberg_dynamic as ed
    from hwave.solver.ir_axis import IRAxis
    beta = 2.0
    ax = IRAxis(beta=beta, wmax=12.0, eps=1e-8, statistics="B")
    vals = (1.0 / (1.0 + (ax.freq_n * np.pi / beta) ** 2)).astype(complex)
    vals = vals[None, :]
    meta = {"freq_n": ax.freq_n, "beta": beta, "wmax": ax.wmax,
            "tol": ax.eps, "L": ax.L, "statistics": "B"}
    out = ed._ir_refit_nodes(vals, meta, ax, "chiq_s", beta)
    np.testing.assert_array_equal(out, vals)

    bad = dict(meta, beta=beta * 1.1)
    with pytest.raises(ValueError, match="beta"):
        ed._ir_refit_nodes(vals, bad, ax, "chiq_s", beta)


def test_ir_native_discriminator_rejects_incomplete_native_markers():
    """A partial native schema must not fall through to a permissive
    uniform-grid reader; only IR provenance by itself denotes densified data.
    """
    from hwave.solver.ir_axis import is_ir_native
    data = {"matsubara_basis": "ir",
            "frequency_grid": "sparse_ir_nodes",
            "ir_freq_n": np.array([-2, 0, 2], dtype=np.int64)}
    assert is_ir_native(data)
    assert not is_ir_native({"matsubara_basis": "ir"})
    with pytest.raises(ValueError, match="frequency_grid"):
        is_ir_native({k: v for k, v in data.items()
                      if k != "frequency_grid"})
    wrong = dict(data, frequency_grid="uniform")
    with pytest.raises(ValueError, match="frequency_grid"):
        is_ir_native(wrong)


@pytest.mark.parametrize("freq_n", [
    np.array([-2.0, 0.5, 2.0]),
    np.array([-2, 2, 0], dtype=np.int64),
])
def test_ir_native_meta_rejects_malformed_nodes(freq_n):
    """Node metadata must not be normalized by dtype coercion or sorting;
    malformed serialized order/type is rejected at the schema boundary.
    """
    from hwave.solver.ir_axis import ir_native_meta
    data = {"ir_freq_n": freq_n, "ir_beta": 2.0, "ir_wmax": 12.0,
            "ir_tol": 1e-8, "ir_L": 3, "ir_statistics": "B"}
    with pytest.raises(ValueError, match="ir_freq_n"):
        ir_native_meta(data)


def test_refit_nodes_rejects_wrong_sector_and_node_count():
    """Even the exact-node pass-through validates the file's statistics and
    stored array length before accepting native values.
    """
    from hwave.solver import eliashberg_dynamic as ed
    from hwave.solver.ir_axis import IRAxis
    beta = 2.0
    ax = IRAxis(beta=beta, wmax=12.0, eps=1e-8, statistics="B")
    vals = np.ones((1, ax.n_freq), dtype=complex)
    meta = {"freq_n": ax.freq_n, "beta": beta, "wmax": ax.wmax,
            "tol": ax.eps, "L": ax.L, "statistics": "F"}
    with pytest.raises(ValueError, match="statistics"):
        ed._ir_refit_nodes(vals, meta, ax, "chiq_s", beta)

    meta["statistics"] = "B"
    with pytest.raises(ValueError, match="frequency-axis length"):
        ed._ir_refit_nodes(vals[..., :-1], meta, ax, "chiq_s", beta)


def test_dynamic_native_memory_guard_uses_node_count(tmp_path, monkeypatch):
    from hwave.solver import eliashberg_dynamic as ed
    out = str(tmp_path)
    solver, _ = _write_converged_run(out, native=True)
    captured = {}
    real = ed.check_memory

    def spy(norb, Nk, nmat, mem_limit_gb=None):
        captured["nmat"] = nmat
        return real(norb, Nk, nmat, mem_limit_gb)

    monkeypatch.setattr(ed, "check_memory", spy)
    ed.solve_dynamic(_dyn_input(out))
    assert captured["nmat"] == solver._ir_axB.n_freq


def test_ir_native_sigma_init_seeds_ir_run(tmp_path, caplog):
    """Sweep chain: run A's native sigma seeds run B (same T and cross-T).
    Same-T seeding reaches the same fixed point as the unseeded run; the
    cross-T seed logs both betas and still converges (design OQ-S3-G)."""
    import logging
    out = str(tmp_path)
    _write_converged_run(out, native=True)

    # same-T seed
    sB, giB = _make_solver(extra_param={'matsubara_basis': 'ir'},
                           iteration_max=60)
    info = sB.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})
    assert info.get("sigma_init_ir") is not None
    giB.update({k: info[k] for k in ("sigma_init", "sigma_init_ir")})
    sB.solve(giB, out)
    assert sB.scf_converged
    s0, gi0 = _make_solver(extra_param={'matsubara_basis': 'ir'},
                           iteration_max=60)
    s0.solve(gi0, out)
    ref = gi0["sigma"]
    np.testing.assert_allclose(giB["sigma"], ref,
                               atol=1e-5 * np.abs(ref).max())

    # cross-T seed (T=0.45): INFO log carries both betas, run converges
    sC, giC = _make_solver(extra_param={'matsubara_basis': 'ir', 'T': 0.45},
                           iteration_max=60)
    infoC = sC.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})
    giC.update({k: infoC[k] for k in ("sigma_init", "sigma_init_ir")})
    with caplog.at_level(logging.INFO, logger="hwave.solver.flex"):
        sC.solve(giC, out)
    assert sC.scf_converged
    assert any("beta" in r.getMessage().lower()
               and "sigma_init" in r.getMessage()
               for r in caplog.records)
