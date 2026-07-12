"""IR-basis (sparse-ir) FLEX for the general (full-vertex / MYO) scheme.

Mirrors tests/test_flex_ir.py but for calc_scheme='general' on a norb=2
ON-SITE fixture. Equivalence is convergence of the uniform result toward the
IR result as Nmat grows (the uniform path carries the O(beta/Nmat) artifact).
Run from the repository root.
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

# Reuse the on-site 2-orbital fixture writer from the general-scheme tests.
from tests.test_flex_general import _write_2d_2orb_onsite_fixture


def _make_general_solver(nmat, matsubara_basis="uniform", T=2.0,
                         iteration_max=60, mix=0.3):
    """A spin-free general FLEX solver on the on-site 2-orbital fixture."""
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    dirpath = os.path.join(tempfile.gettempdir(), "hwave_flex_2d_2orb_onsite")
    _write_2d_2orb_onsite_fixture(dirpath)
    info_input = {'path_to_input': dirpath, 'interaction': {
        'path_to_input': dirpath, 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat',
        'CoulombInter': 'coulombinter.dat'}}
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    param = {'T': T, 'mu': 0.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
             'Nmat': nmat, 'IterationMax': iteration_max, 'Mix': mix, 'EPS': 8,
             'matsubara_basis': matsubara_basis}
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': 'general'}
    gi = read_input_k.QLMSkInput(info_input).get_param("green")
    solver = solver_flex.FLEX(ham, {}, info_mode)
    return solver, gi


def test_grev_index_map_tau_flip_only_spatial_roll_flip():
    """§6.10: the reverse+roll used by the general IR bubble maps tau by
    flip-only (j -> nt-1-j) and each spatial axis by roll(-1)+flip
    ((n-x) mod n); block/orbital axes untouched."""
    xp = np
    g, nt, nx, ny, nz, no2 = 1, 5, 4, 4, 1, 4
    rng = np.random.default_rng(0)
    a = (rng.standard_normal((g, nt, nx, ny, nz, no2))
         + 1j * rng.standard_normal((g, nt, nx, ny, nz, no2)))
    ref = -xp.flip(xp.roll(a, -1, axis=(2, 3, 4)), axis=(1, 2, 3, 4))
    # explicit elementwise map
    out = np.empty_like(a)
    for j in range(nt):
        jm = nt - 1 - j
        for x in range(nx):
            xm = (nx - x) % nx
            for y in range(ny):
                ym = (ny - y) % ny
                for z in range(nz):
                    zm = (nz - z) % nz
                    out[:, j, x, y, z, :] = -a[:, jm, xm, ym, zm, :]
    np.testing.assert_allclose(ref, out, atol=1e-12)


def test_chi0_general_gate_ir_matches_uniform_large_nmat():
    """GATE: the general IR-native chi0 (rank-6) must agree with the uniform
    general chi0 at large Nmat, compressed onto the same bosonic nodes."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 2.0
    beta = 1.0 / T
    diffs = []
    for nmat in (256, 1024):
        s, gi = _make_general_solver(nmat, "ir", T=T)
        s._calc_epsilon_k(gi)
        s._ir_setup(beta)
        axB = s._ir_axB
        g_nodes, _ = s._calc_green_ir(beta, 0.0)
        chi0_ir = s._calc_chi0q_general_ir(g_nodes, beta)[0]      # strip block

        su, giu = _make_general_solver(nmat, "uniform", T=T)
        su._calc_epsilon_k(giu)
        gu, gtail = su._calc_green(beta, 0.0)
        chi0_u = su._calc_chi0q(gu, gtail, beta)[0]               # rank-6 (l,r,a,c,b,d)
        chi0_u_nodes = ed._ir_compress(
            np.moveaxis(chi0_u, 0, -1), axB, nmat, "chi0_u_gen",
            drop_constant=True)
        chi0_u_nodes = np.moveaxis(chi0_u_nodes, -1, 0)

        scale = np.abs(chi0_u_nodes).max()
        diffs.append(np.abs(chi0_ir - chi0_u_nodes).max() / scale)
    assert diffs[-1] < 1e-2, "general chi0 gate failed: {}".format(diffs)
    assert diffs[-1] < 0.6 * diffs[0], \
        "no Nmat convergence toward IR chi0: {}".format(diffs)


def test_chi0_general_ir_orbital_covariance():
    """A norb-permutation P of the bare G permutes the general IR chi0 output
    orbital axes correspondingly (guards swapped MYO pair indices)."""
    T = 2.0
    beta = 1.0 / T
    s, gi = _make_general_solver(256, "ir", T=T)
    s._calc_epsilon_k(gi)
    s._ir_setup(beta)
    g_nodes, _ = s._calc_green_ir(beta, 0.0)          # (g, nF, r, no, no)
    chi0 = s._calc_chi0q_general_ir(g_nodes, beta)     # (g, nB, r, a,c,b,d)
    P = np.array([1, 0])                               # swap the two orbitals
    g_perm = g_nodes[:, :, :, P][:, :, :, :, P]
    chi0_perm = s._calc_chi0q_general_ir(g_perm, beta)
    chi0_from_chi0 = chi0[:, :, :, P][:, :, :, :, P][:, :, :, :, :, P] \
        [:, :, :, :, :, :, P]
    np.testing.assert_allclose(chi0_perm, chi0_from_chi0, atol=1e-10)


def test_sigma_general_ir_direct_matches_uniform():
    """Direct method call: general IR self-energy vs uniform general
    self-energy on the same dressed G and v_eff, compressed onto nodes."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 2.0
    beta = 1.0 / T
    nmat = 1024
    su, giu = _make_general_solver(nmat, "uniform", T=T)
    su._calc_epsilon_k(giu)
    gu, gtail = su._calc_green(beta, 0.0)
    chi0_u = su._calc_chi0q(gu, gtail, beta)[0]
    _, v_eff_u, _, _ = su._flex_compute_veff_general(
        chi0_u, su.ham_info.ham_inter_q)
    sig_u = su._calc_self_energy_general(gu, v_eff_u, beta)   # (g,nmat,r,no,no)

    s, gi = _make_general_solver(nmat, "ir", T=T)
    s._calc_epsilon_k(gi)
    s._ir_setup(beta)
    axF, axB = s._ir_axF, s._ir_axB
    g_nodes, _ = s._calc_green_ir(beta, 0.0)
    chi0_ir = s._calc_chi0q_general_ir(g_nodes, beta)[0]
    _, v_eff_ir, _, _ = s._flex_compute_veff_general(
        chi0_ir, s.ham_info.ham_inter_q)
    sig_ir = s._calc_self_energy_general_ir(g_nodes, v_eff_ir, beta)

    sig_u_nodes = np.moveaxis(ed._ir_compress(
        np.moveaxis(sig_u[0], 0, -1), axF, nmat, "sig_u"), -1, 0)
    scale = np.abs(sig_u_nodes).max()
    assert np.abs(sig_ir[0] - sig_u_nodes).max() / scale < 2e-2


def test_sigma_general_gate_one_iteration():
    """GATE: one general SCF step (sigma=0 -> sigma_new) on the IR path must
    match the uniform general step at large Nmat, compressed onto the
    fermionic nodes."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 2.0
    beta = 1.0 / T
    diffs = []
    for nmat in (256, 1024):
        s, gi = _make_general_solver(nmat, "ir", T=T, iteration_max=1)
        os.makedirs('tests/flex/output', exist_ok=True)
        s.solve(gi, 'tests/flex/output')
        sig_ir = s.sigma                                   # densified (uniform grid)

        su, giu = _make_general_solver(nmat, "uniform", T=T, iteration_max=1)
        su.solve(giu, 'tests/flex/output')
        sig_u = su.sigma
        axF = s._ir_axF
        # compress both onto fermionic nodes for a node-space comparison
        a = np.moveaxis(ed._ir_compress(
            np.moveaxis(sig_ir, 1, -1), axF, nmat, "sig_ir_gen"), -1, 1)
        b = np.moveaxis(ed._ir_compress(
            np.moveaxis(sig_u, 1, -1), axF, nmat, "sig_u_gen"), -1, 1)
        scale = np.abs(b).max()
        diffs.append(np.abs(a - b).max() / scale)
    assert diffs[-1] < 2e-2, "general sigma gate failed: {}".format(diffs)
    assert diffs[-1] < 0.7 * diffs[0], \
        "no Nmat convergence of sigma: {}".format(diffs)


def test_e2e_general_flex_ir_vs_uniform():
    """Converged general FLEX on IR nodes reproduces the uniform general run
    (densified outputs), within the input-Nmat-limited IR tolerance."""
    T = 2.0
    s_ir, gi_ir = _make_general_solver(1024, "ir", T=T, iteration_max=60)
    os.makedirs('tests/flex/output', exist_ok=True)
    s_ir.solve(gi_ir, 'tests/flex/output')
    s_u, gi_u = _make_general_solver(1024, "uniform", T=T, iteration_max=60)
    s_u.solve(gi_u, 'tests/flex/output')
    for key in ("sigma", "chi_s", "chi_c"):
        a = getattr(s_ir, key)
        b = getattr(s_u, key)
        scale = np.abs(b).max()
        assert np.abs(a - b).max() / scale < 3e-2, key


def test_dispatch_routes_general_uniform_and_ir(monkeypatch):
    """The SCF calls exactly _calc_chi0q_general_ir + _calc_self_energy_general_ir
    for general+IR, and the uniform general methods for general+uniform."""
    import hwave.solver.flex as flex
    for basis, chi_name, sig_name in (
            ("ir", "_calc_chi0q_general_ir", "_calc_self_energy_general_ir"),
            ("uniform", "_calc_chi0q", "_calc_self_energy_general")):
        s, gi = _make_general_solver(64, basis, T=2.0, iteration_max=1)
        called = {"chi": 0, "sig": 0}
        orig_chi = getattr(flex.FLEX, chi_name)
        orig_sig = getattr(flex.FLEX, sig_name)

        def wrap_chi(self, *a, _o=orig_chi, **k):
            called["chi"] += 1
            return _o(self, *a, **k)

        def wrap_sig(self, *a, _o=orig_sig, **k):
            called["sig"] += 1
            return _o(self, *a, **k)

        monkeypatch.setattr(flex.FLEX, chi_name, wrap_chi)
        monkeypatch.setattr(flex.FLEX, sig_name, wrap_sig)
        os.makedirs('tests/flex/output', exist_ok=True)
        s.solve(gi, 'tests/flex/output')
        assert called["chi"] >= 1 and called["sig"] >= 1, basis
        monkeypatch.undo()


# ---------------------------------------------------------------------------
# Task 4: dynamic Eliashberg end-to-end (general + IR, densified + native)
# ---------------------------------------------------------------------------

def _fixture_dirpath():
    """Same fixed location `_make_general_solver` writes the on-site
    2-orbital fixture to, so `_run_dynamic_lambda` can point the dynamic
    Eliashberg solver's [file.input.interaction] at the identical geometry
    without re-deriving the path."""
    import tempfile
    return os.path.join(tempfile.gettempdir(), "hwave_flex_2d_2orb_onsite")


def _outdict(path):
    """FLEX output filenames: solve() only populates in-memory arrays, so
    save_results must be called explicitly (mirrors tests/test_flex_ir.py)
    for the dynamic Eliashberg loader to have files to read."""
    return {"path_to_output": path, "chi0q": "chi0q", "sigma": "sigma",
            "green": "green", "chiq_s": "chiq_s", "chiq_c": "chiq_c"}


def _run_dynamic_lambda(flex_out_dir, T, matsubara_basis="uniform"):
    """Run the dynamic Eliashberg solver consuming FLEX outputs in flex_out_dir."""
    from hwave.solver import eliashberg_dynamic as ed
    dirpath = _fixture_dirpath()
    _write_2d_2orb_onsite_fixture(dirpath)
    eli = {"chi0q_mode": "flex", "frequency": "dynamic",
           "pairing_type": "singlet", "solver_mode": "eigenvalue",
           "num_eigenvalues": 4}
    if matsubara_basis != "uniform":
        eli["matsubara_basis"] = matsubara_basis
    inp = {"mode": {"param": {"T": T, "CellShape": [4, 4, 1],
                              "SubShape": [1, 1, 1], "Nmat": 256,
                              "filling": 0.5}},
           "file": {"input": {"path_to_flex_output": flex_out_dir,
                              "interaction": {
                                  "path_to_input": dirpath,
                                  "Geometry": "geom.dat",
                                  "Transfer": "transfer.dat",
                                  "CoulombIntra": "coulombintra.dat",
                                  "CoulombInter": "coulombinter.dat"}},
                    "output": {"path_to_output": flex_out_dir}},
           "eliashberg": eli}
    return ed.solve_dynamic(inp)


def test_dynamic_consumes_general_ir_densified(tmp_path):
    """End-to-end: dynamic Eliashberg lambda from general+IR (densified) FLEX
    outputs matches general+uniform FLEX outputs within tolerance."""
    T = 2.0
    out_ir = str(tmp_path / "gen_ir")
    out_u = str(tmp_path / "gen_u")
    os.makedirs(out_ir, exist_ok=True)
    os.makedirs(out_u, exist_ok=True)
    s_ir, gi_ir = _make_general_solver(256, "ir", T=T, iteration_max=40)
    s_ir.solve(gi_ir, out_ir)
    s_ir.save_results(_outdict(out_ir), gi_ir)
    s_u, gi_u = _make_general_solver(256, "uniform", T=T, iteration_max=40)
    s_u.solve(gi_u, out_u)
    s_u.save_results(_outdict(out_u), gi_u)
    lam_ir = _run_dynamic_lambda(out_ir, T)
    lam_u = _run_dynamic_lambda(out_u, T)
    assert abs(lam_ir - lam_u) / abs(lam_u) < 5e-2, \
        "lam_ir={} lam_u={}".format(lam_ir, lam_u)


def test_dynamic_consumes_general_ir_native(tmp_path):
    """End-to-end with IR-native output (write_densified=false): the dynamic
    solver (matsubara_basis='ir') consumes MYO sparse-node susceptibilities,
    and the resulting lambda matches a general+uniform run within tolerance
    (guards against a finite-but-wrong-parity result on the IR-native
    sparse-node path -- measured lam_ir=~lam_u to well within 5e-2)."""
    T = 2.0
    out_ir_native = str(tmp_path / "gen_ir_native")
    out_u = str(tmp_path / "gen_u")
    os.makedirs(out_ir_native, exist_ok=True)
    os.makedirs(out_u, exist_ok=True)
    s, gi = _make_general_solver(256, "ir", T=T, iteration_max=40)
    s.write_densified = False
    s.solve(gi, out_ir_native)
    s.save_results(_outdict(out_ir_native), gi)
    s_u, gi_u = _make_general_solver(256, "uniform", T=T, iteration_max=40)
    s_u.solve(gi_u, out_u)
    s_u.save_results(_outdict(out_u), gi_u)
    lam = _run_dynamic_lambda(out_ir_native, T, matsubara_basis="ir")
    lam_u = _run_dynamic_lambda(out_u, T)
    assert np.isfinite(lam)
    assert abs(lam - lam_u) / abs(lam_u) < 5e-2, \
        "lam_ir_native={} lam_u={}".format(lam, lam_u)


def test_dynamic_rejects_ir_native_missing_convention(tmp_path):
    """An IR-native chi npz pair (chiq_s AND chiq_c) missing chi_convention
    must be rejected via the missing-convention fail-fast path, not silently
    read as the default (kuroki) convention. (Stripping only one of the two
    files would instead trip the s/c convention-mismatch guard, which
    doesn't exercise the missing-convention check on its own.)"""
    import glob
    out = str(tmp_path / "bad")
    os.makedirs(out, exist_ok=True)
    s, gi = _make_general_solver(64, "ir", T=2.0, iteration_max=1)
    s.write_densified = False
    s.solve(gi, out)
    s.save_results(_outdict(out), gi)
    # strip chi_convention from both chiq_s.npz and chiq_c.npz
    for pattern in ("chiq_s*.npz", "chiq_c*.npz"):
        path = glob.glob(os.path.join(out, pattern))[0]
        d = dict(np.load(path, allow_pickle=True))
        d.pop("chi_convention", None)
        np.savez(path, **d)
    with pytest.raises((ValueError, KeyError)):
        _run_dynamic_lambda(out, 2.0, matsubara_basis="ir")


# ---------------------------------------------------------------------------
# Task 5: CPU/GPU parity, dtype, and static-solver graceful handling
# ---------------------------------------------------------------------------

def _has_cupy_device():
    try:
        import cupy
        cupy.zeros(1)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not _has_cupy_device(), reason="no CuPy device")
def test_general_ir_gpu_matches_cpu():
    """general+IR SCF on GPU reproduces the CPU run (host outputs)."""
    T = 2.0
    s_c, gi_c = _make_general_solver(256, "ir", T=T, iteration_max=10)
    s_c.solve(gi_c, 'tests/flex/output')
    s_g, gi_g = _make_general_solver(256, "ir", T=T, iteration_max=10)
    s_g.use_gpu = True
    s_g.solve(gi_g, 'tests/flex/output')
    for key in ("sigma", "chi_s", "chi_c"):
        a, b = getattr(s_c, key), getattr(s_g, key)
        assert isinstance(b, np.ndarray)               # host-backed
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-8)


def test_general_ir_chi0_dtype_is_complex128():
    T = 2.0
    beta = 1.0 / T
    s, gi = _make_general_solver(64, "ir", T=T)
    s._calc_epsilon_k(gi)
    s._ir_setup(beta)
    g_nodes, _ = s._calc_green_ir(beta, 0.0)
    chi0 = s._calc_chi0q_general_ir(g_nodes, beta)
    assert chi0.dtype == np.complex128


def test_static_chi0q_loader_rejects_general_ir_native(tmp_path):
    """§6.9: the static Eliashberg solver (hwave_sc, [eliashberg] frequency
    left at its "static" default) already fails fast on IR-native chi0q
    files via sc._reject_ir_native/is_ir_native, which inspects only the
    npz's frequency_grid tag and is agnostic to calc_scheme. This asserts
    that existing guard also covers calc_scheme='general' (previously only
    exercised for calc_scheme='reduced' in test_ir_native_readers.py) --
    the static path does not otherwise know about matsubara_basis and never
    silently misreads sparse IR nodes as a uniform grid."""
    import hwave.sc as sc
    out = str(tmp_path)
    s, gi = _make_general_solver(64, "ir", T=2.0, iteration_max=1)
    s.write_densified = False
    s.solve(gi, out)
    s.save_results(_outdict(out), gi)
    input_dict = {"file": {"output": {"path_to_output": out,
                                      "chi0q": "chi0q"}},
                  "mode": {}}
    with pytest.raises(ValueError, match="write_densified"):
        sc._load_chi0q(input_dict)
