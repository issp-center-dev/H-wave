"""Task 8: wiring the dynamic bond kernel into ``solve_dynamic``.

State machine (spec ``2026-07-27-dynamic-bond-channels-design.md`` S5.2), the
chi-file bypass, the green-consistency guards, provenance/diagnostics, and the
invariance of the DEFAULT (bond off) dynamic path.
"""

import logging
import os

import numpy as np
import pytest


# ---------------------------------------------------------------------------
# fixtures: a 4x4 single-band square lattice with U = 4 and NN V = 1
# ---------------------------------------------------------------------------

def _write_model(input_dir, U=4.0, V=1.0):
    """geometry + NN transfer + CoulombIntra(U) + CoulombInter(V), the same
    file formats ``tests/test_eliashberg_dynamic.py`` uses."""
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
        f.write("1\n0.0 0.0 0.0\n")
    with open(os.path.join(input_dir, "transfer.dat"), "w") as f:
        f.write("Transfer\n1\n4\n 1 1 1 1\n")
        for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
            f.write("{:4d} {:4d} {:4d}  1  1  1.0 0.0\n".format(*r))
    with open(os.path.join(input_dir, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n1\n1\n 1\n")
        f.write("  0  0  0  1  1  {:.12f} 0.0\n".format(U))
    with open(os.path.join(input_dir, "coulombinter.dat"), "w") as f:
        f.write("CoulombInter\n1\n4\n 1 1 1 1\n")
        for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
            f.write("{:4d} {:4d} {:4d}  1  1  {:.12f} 0.0\n".format(*r, V))
    return input_dir


def _bond_input(input_dir, output_dir, nmat=16, cell=(4, 4, 1), T=0.5,
                with_inter=True, **eli):
    interaction = {"path_to_input": input_dir,
                   "Geometry": "geom.dat",
                   "Transfer": "transfer.dat",
                   "CoulombIntra": "coulombintra.dat"}
    if with_inter:
        interaction["CoulombInter"] = "coulombinter.dat"
    eli_param = {"frequency": "dynamic", "bond_channels": True,
                 "solver_mode": "eigenvalue", "eigenvalue_method": "arnoldi",
                 "num_eigenvalues": 4}
    eli_param.update(eli)
    return {
        "mode": {"param": {"T": T, "CellShape": list(cell),
                           "SubShape": [1, 1, 1], "Nmat": nmat,
                           "filling": 0.5}},
        "file": {"input": {"interaction": interaction},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": eli_param,
    }


def _write_bare_green_npz(path, input_dict, nmat=None):
    """Write a ``green.npz`` in the FLEX layout holding the BARE green of the
    fixture model -- a legitimate external ``bond_green`` for the file path
    (and, being the same object the fallback builds, it lets a test compare
    the two sources directly)."""
    import hwave.sc as sc

    mode_param = input_dict["mode"]["param"]
    beta = 1.0 / mode_param["T"]
    Lx, Ly, Lz = mode_param["CellShape"]
    nmat = int(nmat if nmat is not None else mode_param["Nmat"])
    geom_info, hr, _ = sc._read_interaction_files(input_dict)
    norb = geom_info["norb"]
    kx = np.linspace(0, 2.0 * np.pi, Lx, endpoint=False)
    ky = np.linspace(0, 2.0 * np.pi, Ly, endpoint=False)
    kz = np.linspace(0, 2.0 * np.pi, Lz, endpoint=False)
    eps = sc._build_hamiltonian_k(kx, ky, kz, hr, norb)
    evals, evecs = sc._calc_eigenvalues(eps)
    mu = sc._determine_mu(evals, beta, mode_param["filling"], norb)
    green = sc._calc_green(evals, evecs, mu, beta, nmat)
    # sc layout (norb, norb, Nx, Ny, Nz, nmat) -> file layout
    # (nblock, nfreq, nvol, norb, norb)
    g = green.transpose(5, 2, 3, 4, 0, 1).reshape(
        1, nmat, Lx * Ly * Lz, norb, norb)
    np.savez(path, green=g)
    return path


# ---------------------------------------------------------------------------
# S5.2 state machine
# ---------------------------------------------------------------------------

def test_dynamic_bond_requires_no_chi_files(tmp_path, caplog):
    """Row 3: bond_channels=true + chi0q_mode set -> the chi files are NOT
    read (none exist here) and the ignored setting is named in a warning."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir, chi0q_mode="flex")
    # No chiq_s.npz / chiq_c.npz / green.npz anywhere: reaching the loader at
    # all would raise.
    assert not os.path.exists(os.path.join(output_dir, "chiq_s.npz"))

    def _must_not_run(*a, **k):
        raise AssertionError("load_flex_chi_dynamic must not run on the "
                             "bond dynamic path")

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(ed, "load_flex_chi_dynamic", _must_not_run)
        with caplog.at_level(logging.WARNING, logger="qlms"):
            lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    assert any("chi0q_mode" in rec.message and "IGNORED" in rec.message
               for rec in caplog.records)


def test_dynamic_bond_bare_green_fallback_warns(tmp_path, caplog):
    """Row 4: no bond_green -> the bare H0 green is used, with a prominent
    warning, and provenance records bond_green_source='bare'."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    assert any("BARE Green function" in rec.message for rec in caplog.records)
    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert str(meta["bond_green_source"]) == "bare"
    assert bool(meta["bond_channels"])


def test_dynamic_bond_file_green_matches_bare_equivalent(tmp_path):
    """Row 2: with bond_green set the file is used (provenance says so) and,
    for a file holding the very green the fallback would build, the two give
    the same lambda -- i.e. the file path is a faithful substitution, not a
    different convention."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    out_bare = str(tmp_path / "out_bare")
    out_file = str(tmp_path / "out_file")
    inp_bare = _bond_input(input_dir, out_bare)
    lam_bare = ed.solve_dynamic(inp_bare)

    green_path = str(tmp_path / "green.npz")
    _write_bare_green_npz(green_path, inp_bare)
    lam_file = ed.solve_dynamic(
        _bond_input(input_dir, out_file, bond_green=green_path))
    assert np.isclose(lam_file, lam_bare, rtol=1e-10, atol=1e-12)
    meta = np.load(os.path.join(out_file, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert str(meta["bond_green_source"]) == "file"
    assert str(meta["bond_green"]) == green_path


@pytest.mark.parametrize("mismatch", ["nmat", "cellshape"])
def test_dynamic_bond_green_mismatch_errors(tmp_path, mismatch):
    """A bond_green on a different Matsubara grid or k-grid than the run
    config is an ERROR (spec S3.4: no silent mixed-convention runs)."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir, nmat=16)
    green_path = str(tmp_path / "green.npz")
    if mismatch == "nmat":
        _write_bare_green_npz(green_path, inp, nmat=8)
        pattern = "Matsubara frequencies"
    else:
        other = _bond_input(input_dir, output_dir, nmat=16, cell=(2, 2, 1))
        _write_bare_green_npz(green_path, other)
        pattern = "does not match the model"
    inp["eliashberg"]["bond_green"] = green_path
    with pytest.raises(ValueError, match=pattern):
        ed.solve_dynamic(inp)


def test_dynamic_bond_requires_coulomb_inter(tmp_path):
    """S6: the bond guards of the static path apply unchanged -- a run with
    no CoulombInter term declared is refused."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "output"), with_inter=False)
    with pytest.raises(ValueError, match="(?i)CoulombInter"):
        ed.solve_dynamic(inp)


def test_dynamic_bond_prereqs_no_longer_require_flex_chi0q_mode():
    """The chi0q_mode='flex' requirement exists because the SCALAR dynamic
    vertex ingests FLEX chi files; the bond branch reads none, so it must not
    inherit the requirement (while the scalar path still does)."""
    import hwave.sc as sc

    base = {"mode": {"param": {"Nmat": 8}},
            "eliashberg": {"frequency": "dynamic"}}
    with pytest.raises(ValueError, match="flex"):
        sc._validate_dynamic_prereqs(base)
    bond = {"mode": {"param": {"Nmat": 8}},
            "eliashberg": {"frequency": "dynamic", "bond_channels": True}}
    sc._validate_dynamic_prereqs(bond)          # must NOT raise
    odd = {"mode": {"param": {"Nmat": 7}},
           "eliashberg": {"frequency": "dynamic", "bond_channels": True}}
    with pytest.raises(ValueError, match="even"):
        sc._validate_dynamic_prereqs(odd)


# ---------------------------------------------------------------------------
# default-path invariance
# ---------------------------------------------------------------------------

def test_dynamic_default_path_unchanged(tmp_path, monkeypatch):
    """With bond_channels absent, NONE of the new code paths may be entered:
    every dynamic bond entry point is monkeypatched to raise, and the
    existing scalar dynamic run must still produce its usual outputs.

    (The complementary half of the invariance evidence is the unchanged
    tests/test_eliashberg_dynamic.py and tests/test_eliashberg_ir.py suites,
    which pin the scalar path's numbers themselves.)"""
    from hwave.solver import bond_channels as bc
    from hwave.solver import eliashberg_dynamic as ed
    from tests.test_eliashberg_dynamic import (_write_flex_fixture,
                                               _write_geom_transfer_coulomb,
                                               _dynamic_input_dict)

    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1,
                            Nx=2, Ny=2, Nz=1)

    def _boom(*a, **k):
        raise AssertionError("bond code path entered with bond_channels off")

    for name in ("bond_bubble_dynamic", "dress_bond_dynamic",
                 "make_bond_kernel_dynamic", "make_bond_kernel_dynamic_ir"):
        monkeypatch.setattr(bc, name, _boom)
    monkeypatch.setattr(ed, "_build_bond_dynamic_context", _boom)

    lam = ed.solve_dynamic(_dynamic_input_dict(input_dir, output_dir,
                                               m["nmat"]))
    assert np.isfinite(lam)
    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    # the default output key set is untouched: no bond provenance leaks in
    assert "bond_channels" not in meta.files
    assert "gap" in meta.files and "eigenvalue" in meta.files


# ---------------------------------------------------------------------------
# end-to-end smoke (both parities, both Matsubara bases)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
@pytest.mark.parametrize("basis", ["uniform", "ir"])
def test_dynamic_bond_e2e_smoke(tmp_path, pairing, basis):
    """4x4, Nmat=16, U=4, V=1: solve_dynamic returns a finite real lambda and
    writes the provenance + diagnostics the spec requires."""
    if basis == "ir":
        pytest.importorskip("sparse_ir")
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir, nmat=16,
                      pairing_type=pairing, matsubara_basis=basis,
                      bond_diagnostics=True)
    lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    assert isinstance(lam, float)

    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert bool(meta["bond_channels"])
    assert int(meta["bond_n_channels"]) == 5          # Delta r = 0, +-x, +-y
    assert str(meta["bond_green_source"]) == "bare"
    assert str(meta["bond_frequency_grid"]) == basis
    assert float(meta["bond_cond_min_spin"]) > 0.0
    assert float(meta["bond_cond_min_charge"]) > 0.0
    assert np.isfinite(float(meta["bond_lambda_imag"]))
    # spec S2 (amended): the Hermitian-metric feature set is uniform-only
    assert str(meta["hermitian_metric_features"]) == (
        "uniform-only" if basis == "ir" else "uniform")
    # the gap itself is always written back on the uniform output grid
    assert meta["gap"].shape == (1, 1, 4, 4, 1, 16)

    diag = np.load(os.path.join(output_dir, "bond_diagnostics.npz"),
                   allow_pickle=False)
    assert diag["cond_min_spin"].shape == (4, 4, 1, 16)
    assert diag["cond_min_charge"].shape == (4, 4, 1, 16)
    assert "bond_diag_axes" in diag.files
    assert np.isfinite(float(diag["tail_est"]))
    assert str(diag["tail_est_branch"]) == pairing
    assert bool(diag["tail_est_unreliable"]) in (True, False)


def test_dynamic_bond_ir_retains_the_static_chi_component(tmp_path, caplog):
    """The dressed BOND susceptibility is static-dominated, so its
    frequency-independent component is retained on the IR path.

    At Nmat=32 on this fixture the DROP policy trips the issue-#57 guard
    ("the discarded frequency-independent component exceeds the data
    scale"), so a run that completes here is direct evidence the retain
    branch is the one taken -- and an explicit ir_keep_static_chi=false must
    be refused loudly rather than silently overridden."""
    pytest.importorskip("sparse_ir")
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "out"), nmat=32,
                      matsubara_basis="ir", ir_keep_static_chi=False)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    assert any("ir_keep_static_chi=false is NOT honoured" in rec.message
               for rec in caplog.records)
    meta = np.load(os.path.join(str(tmp_path / "out"), "gap_dynamic.npz"),
                   allow_pickle=False)
    assert str(meta["bond_ir_static_chi"]) == "retained"


def test_dynamic_bond_ir_matches_the_uniform_operator(tmp_path):
    """Cross-basis agreement (spec S3.5): with the same physical input, the
    IR-axis bond kernel and the oracle-verified uniform one must give the
    same leading lambda up to the instantaneous-term convention difference,
    which shrinks with the window. At Nmat=64 that is ~3e-3 relative."""
    pytest.importorskip("sparse_ir")
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    lam_u = ed.solve_dynamic(_bond_input(input_dir, str(tmp_path / "u"),
                                         nmat=64))
    lam_ir = ed.solve_dynamic(_bond_input(input_dir, str(tmp_path / "i"),
                                          nmat=64, matsubara_basis="ir"))
    assert abs(lam_ir - lam_u) <= 5.0e-3 * abs(lam_u), (lam_u, lam_ir)


def test_dynamic_bond_diagnostics_are_opt_in(tmp_path):
    """No bond_diagnostics key -> no diagnostics npz (and no provenance flag
    claiming one was written)."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    ed.solve_dynamic(_bond_input(input_dir, output_dir))
    assert not os.path.exists(
        os.path.join(output_dir, "bond_diagnostics.npz"))
    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert "bond_diagnostics" not in meta.files


def test_dynamic_bond_rejects_asymmetric_green(tmp_path):
    """Spec S2: a green outside the v1 symmetry class is REJECTED, not
    silently whitened -- the real-spectrum solver contract depends on it."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir, nmat=16)
    green_path = str(tmp_path / "green.npz")
    _write_bare_green_npz(green_path, inp)
    data = dict(np.load(green_path))
    rng = np.random.default_rng(0)
    data["green"] = data["green"] + 0.1 * (
        rng.standard_normal(data["green"].shape)
        + 1j * rng.standard_normal(data["green"].shape))
    np.savez(green_path, **data)
    inp["eliashberg"]["bond_green"] = green_path
    with pytest.raises(ValueError, match="symmetry class"):
        ed.solve_dynamic(inp)


# ---------------------------------------------------------------------------
# preflight accounting (hand-off note: the dynamic budget must not undercount)
# ---------------------------------------------------------------------------

def test_dynamic_preflight_counts_every_frequency_resolved_buffer():
    """The dynamic estimate must scale EVERY frequency-resolved family --
    dressed chi_s/chi_c, the kernel's hoisted vertex, the matvec buffers and
    G2p_w -- not only chi_bar. A preflight that undercounts waves through a
    run that then OOMs anyway."""
    import hwave.sc as sc
    from hwave.solver import bond_channels as bc

    inter = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
             ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
    bset = bc.resolve_interactions(inter, np.eye(3), 1)
    assert bset.n_channels == 5
    kw = dict(norb=1, bond_set=bset, Nx=8, Ny=8, Nz=1, nmat=32)
    est = sc._bond_memory_estimate(dynamic_nmat=32, **kw)

    Nq, ND, nd = 64, 5, 1
    unit_dyn = Nq * 32 * ND * ND * 16
    unit_mv = Nq * 32 * ND * 16
    assert est["g2_bytes"] == nd * nd * Nq * 32 * 16
    # dressing: chi_bar + I + mat_s + mat_c + chi_s + chi_c
    assert est["dress_bytes"] == sc._BOND_DRESS_DYN_ARRAYS * unit_dyn
    # kernel build: caller-held chi_s/chi_c + the hoist high-water + G2p_w
    assert est["kernel_build_bytes"] == (
        (2 + bc.BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS) * unit_dyn
        + est["g2_bytes"])
    # matvec: chi_s/chi_c + persistent F_rt + matvec-class buffers + G2p_w
    assert est["kernel_matvec_bytes"] == (
        (2 + bc.BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT) * unit_dyn
        + bc.BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS * unit_mv + est["g2_bytes"])
    # the peak is a phase MAXIMUM plus what lives across every phase
    assert est["peak"] == (est["green_bytes"] + est["q_bytes"]
                           + est["vertex_bytes"]
                           + max(est["bubble_bytes"] + unit_dyn,
                                 est["dress_bytes"],
                                 est["kernel_build_bytes"],
                                 est["kernel_matvec_bytes"]))
    # and it must dominate the static one at the same nmat
    assert est["peak"] > sc._bond_memory_estimate(**kw)["peak"] * 10


def test_dynamic_preflight_refuses_over_cap(tmp_path):
    """bond_memory_cap_gb applies to the dynamic path with the same
    fail-before-allocating semantics."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "output"),
                      bond_memory_cap_gb=1.0e-7)
    with pytest.raises(ValueError, match="bond_memory_cap_gb"):
        ed.solve_dynamic(inp)
