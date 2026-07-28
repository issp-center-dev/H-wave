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
    # dressing: chi_bar (caller-held) + mat_s + mat_c + chi_s + chi_c + the
    # matmul temporary -- a MEASURED constant, pinned by
    # test_preflight_dress_budget_covers_the_measured_dress_peak
    assert est["dress_bytes"] == sc._BOND_DRESS_DYN_ARRAYS * unit_dyn
    # kernel build: caller-held chi_s/chi_c + the hoist high-water + G2p_w
    assert est["kernel_build_bytes"] == (
        (2 + bc.BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS) * unit_dyn
        + est["g2_bytes"])
    # matvec: persistent F_rt + matvec-class buffers + G2p_w + the ARPACK
    # basis. chi_s/chi_c are NOT here -- the wiring releases both before any
    # matvec runs.
    assert est["kernel_matvec_bytes"] == (
        bc.BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT * unit_dyn
        + bc.BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS * unit_mv + est["g2_bytes"]
        + est["arpack_bytes"])
    # the peak is a phase MAXIMUM plus what lives across every phase
    assert est["peak"] == (est["green_bytes"] + est["q_bytes"]
                           + est["vertex_bytes"]
                           + max(est["bubble_bytes"] + unit_dyn,
                                 est["dress_bytes"],
                                 est["kernel_build_bytes"],
                                 est["kernel_matvec_bytes"]))
    # and it must dominate the static one at the same nmat
    assert est["peak"] > sc._bond_memory_estimate(**kw)["peak"] * 10


def test_preflight_dress_budget_covers_the_measured_dress_peak():
    """The dressing budget must not UNDERCOUNT ``dress_bond_dynamic``'s real
    high-water mark -- the frequency-resolved sibling of
    ``test_preflight_bubble_budget_covers_the_measured_bond_bubble_peak``,
    and the test class that would have caught the (now fixed)
    ``S_full``/``C_full`` reshape copies: a broadcast view cannot be reshaped
    to a flat batch without materialising it, so the function really needed
    8 units where the budget granted 6."""
    from hwave.solver import bond_channels as bc
    from tests.test_sc_bond import _measure_peak_bytes
    import hwave.sc as sc

    norb, Nx, Ny, Nz, nmat = 1, 4, 4, 1, 16
    inter = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
             ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
    bond_set = bc.resolve_interactions(inter, np.eye(3), norb)
    ND = norb * norb * bond_set.n_channels

    rng = np.random.default_rng(0)

    def rc(shape):
        return 0.05 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))

    chi_bar = rc((Nx, Ny, Nz, nmat, ND, ND))
    S_bond = 6.0 * rc((Nx, Ny, Nz, ND, ND))
    C_bond = 6.0 * rc((Nx, Ny, Nz, ND, ND))

    bc.dress_bond_dynamic(chi_bar, S_bond, C_bond, cond_tol=None)  # warm up
    peak = _measure_peak_bytes(
        lambda: bc.dress_bond_dynamic(chi_bar, S_bond, C_bond, cond_tol=None))
    if peak is None:
        pytest.skip("tracemalloc does not track numpy array data here")

    est = sc._bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat,
                                   dynamic_nmat=nmat)
    # chi_bar is allocated outside the measured region but IS inside the
    # dressing budget (the caller holds it across the call), so the budget
    # the measured peak must fit into is the dressing line minus that array.
    unit_dyn = Nx * Ny * Nz * nmat * ND * ND * 16
    budget = est["dress_bytes"] - unit_dyn
    assert budget >= peak, (
        "the preflight budgets {:.3f} MB for dress_bond_dynamic's own "
        "buffers but the measured peak is {:.3f} MB -- "
        "_BOND_DRESS_DYN_ARRAYS has desynced from dress_bond_dynamic".format(
            budget / 1e6, peak / 1e6))
    assert budget <= 3.0 * peak


def test_dress_bond_dynamic_matches_the_materialized_vertex_form():
    """The broadcast-matmul memory fix must be BIT-identical to the explicit
    ``S_full``/``C_full`` form it replaced (reproduced here, so the equality
    is against the algebra rather than against the same code)."""
    from hwave.solver import bond_channels as bc

    rng = np.random.default_rng(3)
    Nx, Ny, Nz, nmat, ND = 3, 2, 1, 8, 4

    def rc(shape):
        return 0.05 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))

    chi_bar = rc((Nx, Ny, Nz, nmat, ND, ND))
    S_bond = 5.0 * rc((Nx, Ny, Nz, ND, ND))
    C_bond = 5.0 * rc((Nx, Ny, Nz, ND, ND))

    N = Nx * Ny * Nz * nmat
    chi_flat = chi_bar.reshape(N, ND, ND)
    S_flat = np.broadcast_to(S_bond[:, :, :, None, :, :],
                             chi_bar.shape).reshape(N, ND, ND)
    C_flat = np.broadcast_to(C_bond[:, :, :, None, :, :],
                             chi_bar.shape).reshape(N, ND, ND)
    I_mat = np.broadcast_to(np.eye(ND, dtype=complex), (N, ND, ND)).copy()
    ref_s = np.linalg.solve(I_mat - chi_flat @ S_flat,
                            chi_flat).reshape(chi_bar.shape)
    ref_c = np.linalg.solve(I_mat + chi_flat @ C_flat,
                            chi_flat).reshape(chi_bar.shape)

    chi_s, chi_c, _, _ = bc.dress_bond_dynamic(chi_bar, S_bond, C_bond,
                                               cond_tol=None)
    assert np.abs(chi_s - ref_s).max() == 0.0
    assert np.abs(chi_c - ref_c).max() == 0.0


def test_dynamic_preflight_counts_the_eigensolver_basis():
    """The ARPACK basis is its own budget line (spec S7): it scales with
    num_eigenvalues and, at B = 1, is comparable to the whole matvec family
    -- it was previously uncounted."""
    import hwave.sc as sc
    from hwave.solver import bond_channels as bc

    inter = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0}
    bset = bc.resolve_interactions(inter, np.eye(3), 1)
    kw = dict(norb=1, bond_set=bset, Nx=8, Ny=8, Nz=1, nmat=32,
              dynamic_nmat=32)
    vec_size_bytes = 1 * 64 * 32 * 16
    est4 = sc._bond_memory_estimate(num_eigenvalues=4, **kw)
    est16 = sc._bond_memory_estimate(num_eigenvalues=16, **kw)
    assert est4["arpack_bytes"] == 9 * vec_size_bytes
    assert est16["arpack_bytes"] == 33 * vec_size_bytes
    assert (est16["kernel_matvec_bytes"] - est4["kernel_matvec_bytes"]
            == est16["arpack_bytes"] - est4["arpack_bytes"])
    assert sc._bond_memory_estimate(**kw)["arpack_bytes"] == (
        (2 * sc._BOND_DEFAULT_NUM_EIGENVALUES + 1) * vec_size_bytes)


def test_bond_memory_estimate_refuses_mixed_grids():
    """``nmat`` and ``dynamic_nmat`` describe ONE grid on this path; two
    different values would budget a run that does not exist."""
    import hwave.sc as sc
    from hwave.solver import bond_channels as bc

    bset = bc.resolve_interactions({((1, 0, 0), (0, 0)): 1.0,
                                    ((-1, 0, 0), (0, 0)): 1.0}, np.eye(3), 1)
    with pytest.raises(ValueError, match="must equal nmat"):
        sc._bond_memory_estimate(1, bset, 4, 4, 1, 16, dynamic_nmat=32)


def test_dynamic_bond_rejects_ir_native_green(tmp_path):
    """An IR-NATIVE (sparse-node) bond_green must be diagnosed as such, not
    as a Matsubara-count mismatch: the header peek sees its NODE COUNT, so
    before the fix the user was told to 'regenerate at Nmat = <node count>'."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "output"), nmat=16)
    green_path = str(tmp_path / "green_ir.npz")
    _write_bare_green_npz(green_path, inp, nmat=16)
    data = dict(np.load(green_path))
    # the three markers ir_axis.is_ir_native requires, on a node count that
    # is NOT the run's Nmat -- exactly the case that used to mis-report
    data["green"] = data["green"][:, :7]
    data["matsubara_basis"] = "ir"
    data["frequency_grid"] = "sparse_ir_nodes"
    data["ir_freq_n"] = np.arange(7)
    np.savez(green_path, **data)
    inp["eliashberg"]["bond_green"] = green_path
    with pytest.raises(ValueError) as exc:
        ed.solve_dynamic(inp)
    msg = str(exc.value)
    assert "sparse_ir_nodes" in msg and "write_densified" in msg
    assert "Nmat = 7" not in msg


@pytest.mark.parametrize("case", ["match", "mismatch"])
def test_bond_green_metadata_temperature_check(tmp_path, case):
    """``_bond_green_metadata_consistency`` compares the file's own beta/nmat
    metadata with the run: equal passes silently, different is an ERROR (a
    green at another temperature lives on another Matsubara grid entirely)."""
    from hwave.solver import eliashberg_dynamic as ed

    path = str(tmp_path / "g.npz")
    beta_run = 2.0
    np.savez(path, green=np.zeros((1, 16, 4, 1, 1), dtype=complex),
             beta=(beta_run if case == "match" else 4.0), nmat=16)
    if case == "match":
        ed._bond_green_metadata_consistency(path, 16, beta_run)   # no raise
    else:
        with pytest.raises(ValueError, match="SAME temperature"):
            ed._bond_green_metadata_consistency(path, 16, beta_run)

    # the nmat metadata is checked independently of the temperature
    np.savez(path, green=np.zeros((1, 16, 4, 1, 1), dtype=complex),
             beta=beta_run, nmat=32)
    with pytest.raises(ValueError, match="records nmat = 32"):
        ed._bond_green_metadata_consistency(path, 16, beta_run)

    # a file with no metadata at all is REPORTED, not refused
    np.savez(path, green=np.zeros((1, 16, 4, 1, 1), dtype=complex))
    ed._bond_green_metadata_consistency(path, 16, beta_run)


def test_dynamic_bond_iteration_mode_warns_about_the_repulsive_kernel(
        tmp_path, caplog):
    """solver_mode='iteration' without spectral_shift is a silent-wrong-result
    trap on a repulsion-dominated bond kernel; it must be named."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "out"), nmat=16,
                      solver_mode="iteration", max_iter=20)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        lam_iter = ed.solve_dynamic(inp)
    assert any("repulsion-" in rec.message and "spectral_shift" in rec.message
               for rec in caplog.records)
    # and the trap is real: the ARPACK leading value differs
    lam_eig = ed.solve_dynamic(_bond_input(input_dir, str(tmp_path / "out2"),
                                           nmat=16))
    assert not np.isclose(lam_iter, lam_eig, rtol=1e-3)

    # no warning once a shift is supplied
    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="qlms"):
        ed.solve_dynamic(_bond_input(input_dir, str(tmp_path / "out3"),
                                     nmat=16, solver_mode="iteration",
                                     max_iter=20, spectral_shift="auto"))
    assert not any("repulsion-" in rec.message for rec in caplog.records)


def test_dynamic_bond_does_not_claim_the_zero_chi_diagnostic(tmp_path, caplog):
    """zero_chi_s/zero_chi_c are ignored on the bond branch, so neither
    eigenvalue.dat nor the gap npz may record them: an output naming a
    diagnostic that never ran is worse than no output."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    output_dir = str(tmp_path / "output")
    inp = _bond_input(input_dir, output_dir, nmat=16, zero_chi_s=True,
                      zero_chi_c=True)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        ed.solve_dynamic(inp)
    assert any("zero_chi_s/zero_chi_c are IGNORED" in rec.message
               for rec in caplog.records)
    with open(os.path.join(output_dir, "eigenvalue.dat")) as f:
        assert "zero_chi_s" not in f.read()
    meta = np.load(os.path.join(output_dir, "gap_dynamic.npz"),
                   allow_pickle=False)
    assert "zero_chi_s" not in meta.files
    assert "zero_chi_c" not in meta.files


def test_dynamic_preflight_refuses_over_cap(tmp_path):
    """bond_memory_cap_gb applies to the dynamic path with the same
    fail-before-allocating semantics."""
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = _write_model(str(tmp_path / "input"))
    inp = _bond_input(input_dir, str(tmp_path / "output"),
                      bond_memory_cap_gb=1.0e-7)
    with pytest.raises(ValueError, match="bond_memory_cap_gb"):
        ed.solve_dynamic(inp)
