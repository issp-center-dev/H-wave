"""IR-basis (sparse-ir) path of the dynamic Eliashberg solver (Stage 1).

Design: docs/design/ir-matsubara.md. The fixtures run a REAL small FLEX
calculation first so every frequency object is physically smooth (pole
structure) and therefore IR-representable -- random-tensor fixtures are
white noise on the frequency axis and must not be used here.

Tests must run from the repository root.
"""
import os
import unittest

import numpy as np
import pytest

# Import-safe under BOTH pytest and unittest discovery (the CI runs
# `unittest discover`, where pytest.importorskip at module level would turn
# a missing optional dependency into an import ERROR instead of a skip).
try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

pytestmark = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                                reason="sparse-ir not installed")


BETA_T = 0.5          # temperature of the mini pipeline
NMAT = 128
LX = LY = 4


def _run_flex_to(outdir, gpu=False, nmat=NMAT):
    """Run the small 1-orbital Hubbard FLEX and write full-grid npz outputs."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    info_mode = {'mode': 'FLEX', 'param': {
        'T': BETA_T, 'mu': 0.0, 'CellShape': [LX, LY, 1],
        'SubShape': [1, 1, 1], 'Nmat': nmat,
        'IterationMax': 30, 'Mix': 0.3, 'EPS': 8, 'gpu': gpu},
        'calc_scheme': 'reduced'}
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'}}
    rio = read_input_k.QLMSkInput(info_input)
    ham = rio.get_param("ham")
    gi = rio.get_param("green")
    s = solver_flex.FLEX(ham, {}, info_mode)
    os.makedirs(outdir, exist_ok=True)
    s.solve(gi, outdir)
    s.save_results({'path_to_output': outdir, 'chiq_s': 'chiq_s.npz',
                    'chiq_c': 'chiq_c.npz', 'green': 'green.npz',
                    'sigma': 'sigma.npz'}, gi)
    return gi


def _eliashberg_input(outdir, extra=None):
    eli = {"chi0q_mode": "flex", "frequency": "dynamic",
           "solver_mode": "iteration", "max_iter": 200,
           "convergence_tol": 1e-7}
    if extra:
        eli.update(extra)
    return {
        "mode": {"param": {"T": BETA_T, "CellShape": [LX, LY, 1],
                           "SubShape": [1, 1, 1], "Nmat": NMAT,
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": "tests/rpa/input",
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": outdir}},
        "eliashberg": eli,
    }


@pytest.fixture(scope="module")
def flex_outdir(tmp_path_factory):
    outdir = str(tmp_path_factory.mktemp("flexout"))
    _run_flex_to(outdir)
    return outdir


def test_ir_matvec_matches_uniform_kernel(flex_outdir):
    """Acceptance gate (design Sec. 2): the IR kernel applied to an
    IR-representable probe must reproduce the full-grid kernel after
    densification, BEFORE any eigensolver runs."""
    from hwave.solver import eliashberg_dynamic as ed
    import hwave.sc as sc

    inp = _eliashberg_input(flex_outdir)
    norb, Nx, Ny, Nz = 1, LX, LY, 1
    beta = 1.0 / BETA_T

    # full-grid inputs
    chis_w, chic_w, green_w, conv = ed.load_flex_chi_dynamic(
        inp, norb, Nx, Ny, Nz)
    kx = np.linspace(0, 2*np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2*np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2*np.pi, Nz, endpoint=False)
    geom, hr, inter = sc._read_interaction_files(inp)
    inter_k = sc._build_interaction_k(kx, ky, kz, inter, norb)

    # IR-path objects from the SAME npz
    axF, axB = ed._ir_axes_for_run(inp["eliashberg"], beta, hr, inter_k, norb)
    chis_n = ed._ir_compress(chis_w, axB, NMAT, "chiq_s", drop_constant=True)
    chic_n = ed._ir_compress(chic_w, axB, NMAT, "chiq_c", drop_constant=True)
    green_n = ed._ir_compress(green_w, axF, NMAT, "green")
    V_n = ed.compute_vertices_flex_dynamic(chis_n, chic_n, inter_k, norb,
                                           Nx, Ny, Nz,
                                           pairing_type="singlet",
                                           convention=conv)
    G2_n = ed.calc_g2_dynamic(green_n, beta)
    V_rt_tau = ed._ir_vertex_to_rtau(V_n, axB, axF)

    # KERNEL-ALGEBRA gate: feed the uniform kernel the DENSIFIED
    # (IR-representable) chi/G so both kernels see identical data; the raw
    # uniform-FFT chi additionally carries finite-Nmat artifacts (the
    # delta(tau) constant and aliasing images) that no basis of bandwidth
    # wmax can or should represent -- data quality is covered by the
    # end-to-end lambda test below, not by this gate.
    chis_d = axB.eval_to_uniform(axB.fit_from_freq(chis_n), NMAT)
    chic_d = axB.eval_to_uniform(axB.fit_from_freq(chic_n), NMAT)
    green_d = axF.eval_to_uniform(axF.fit_from_freq(green_n), NMAT)
    V_w = ed.compute_vertices_flex_dynamic(chis_d, chic_d, inter_k, norb,
                                           Nx, Ny, Nz,
                                           pairing_type="singlet",
                                           convention=conv)
    G2_w = ed.calc_g2_dynamic(green_d, beta)

    # IR-representable probe: the frequency-flat d-wave-like seed
    gap_static = sc._initialize_gap("cos", norb, kx, ky, kz)
    phi_u = np.broadcast_to(gap_static[..., None],
                            gap_static.shape + (NMAT,)).astype(complex)
    phi_n = np.broadcast_to(gap_static[..., None],
                            gap_static.shape + (axF.n_freq,)).astype(complex)

    out_u = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi_u.copy(), norb, beta)
    out_n = ed.eliashberg_kernel_ir(V_rt_tau, G2_n, phi_n.copy(), axF, beta)
    out_n_dense = axF.eval_to_uniform(axF.fit_from_freq(out_n), NMAT)

    scale = np.abs(out_u).max()
    diff_default = np.abs(out_n_dense - out_u).max() / scale

    # The two kernels are different discretizations of the same operator.
    # On REAL FLEX data the max-norm difference is bounded by the input's
    # out-of-basis content (uniform-FFT artifacts propagated into the
    # vertex): it concentrates on the HIGH-|n| nodes (measured 1e-4 at
    # |n| <= 5 vs O(0.1) at the tail on this fixture) and is NON-monotonic
    # in wmax, while the low-frequency sector that determines lambda
    # agrees -- which is what the end-to-end lambda test pins. Assert the
    # loose whole-axis bound here...
    assert diff_default < 5e-1
    eli2 = dict(inp["eliashberg"]); eli2["ir_wmax"] = 200.0
    axF2, axB2 = ed._ir_axes_for_run(eli2, beta, hr, inter_k, norb)
    # This is a kernel-ALGEBRA gate (operator equivalence on identical
    # densified data), not a data-fidelity check: at this deliberately large
    # wmax the delta(tau) constant is big (ill-conditioned drop), which the
    # production solve would reject -- opt out of that guard here.
    chis2 = ed._ir_compress(chis_w, axB2, NMAT, "chiq_s", drop_constant=True,
                            error_on_large_constant=False)
    chic2 = ed._ir_compress(chic_w, axB2, NMAT, "chiq_c", drop_constant=True,
                            error_on_large_constant=False)
    green2 = ed._ir_compress(green_w, axF2, NMAT, "green")
    V2 = ed.compute_vertices_flex_dynamic(chis2, chic2, inter_k, norb,
                                          Nx, Ny, Nz, pairing_type="singlet",
                                          convention=conv)
    G22 = ed.calc_g2_dynamic(green2, beta)
    Vrt2 = ed._ir_vertex_to_rtau(V2, axB2, axF2)
    chis_d2 = axB2.eval_to_uniform(axB2.fit_from_freq(chis2), NMAT)
    chic_d2 = axB2.eval_to_uniform(axB2.fit_from_freq(chic2), NMAT)
    green_d2 = axF2.eval_to_uniform(axF2.fit_from_freq(green2), NMAT)
    Vw2 = ed.compute_vertices_flex_dynamic(chis_d2, chic_d2, inter_k, norb,
                                           Nx, Ny, Nz, pairing_type="singlet",
                                           convention=conv)
    G2w2 = ed.calc_g2_dynamic(green_d2, beta)
    phi_n2 = np.broadcast_to(gap_static[..., None],
                             gap_static.shape + (axF2.n_freq,)).astype(complex)
    out_u2 = ed.eliashberg_kernel_dynamic(Vw2, G2w2, phi_u.copy(), norb, beta)
    out_n2 = ed.eliashberg_kernel_ir(Vrt2, G22, phi_n2.copy(), axF2, beta)
    out_nd2 = axF2.eval_to_uniform(axF2.fit_from_freq(out_n2), NMAT)
    diff_wide = np.abs(out_nd2 - out_u2).max() / np.abs(out_u2).max()
    assert diff_wide < 1e-2, "wide-basis kernels disagree: {}".format(diff_wide)
    assert diff_wide < 0.2 * diff_default, (
        "no systematic wmax convergence: {} vs {}".format(diff_wide,
                                                          diff_default))
    # ...and prove the KERNEL ALGEBRA itself is exact with a fully
    # IR-representable synthetic vertex (no chi artifacts, no fit escape
    # hatches): a smooth Lorentzian-profile V and a bare-G2, where both
    # kernels must agree to numerical precision (issue #57).
    assert _smooth_vertex_gate(2, 4, 4, 1, beta=2.0, nmat=256,
                               wmax=20.0) < 1e-6
    assert _smooth_vertex_gate(2, 1, 8, 8, beta=200.0, nmat=1024,
                               wmax=5.0) < 1e-4


def test_solve_dynamic_ir_matches_uniform_lambda(flex_outdir, tmp_path):
    """End-to-end: the IR and uniform leading eigenvalues are two
    discretizations of the same continuum answer, and their difference is
    bounded by the input-data quality (the finite-Nmat artifacts of the
    uniform-FFT FLEX output): percent-level at Nmat=128 and shrinking
    systematically as the input Nmat grows (measured 1.5e-2/0.395 at 128
    with the issue-#57 dispersion-based auto wmax, previously 1.0e-2 with
    the inflated pre-#57 estimate whose extra bandwidth over-resolved the
    input artifacts). The IR run writes the gap on the SAME uniform
    grid/metadata, with IR provenance recorded."""
    from hwave.solver import eliashberg_dynamic as ed

    lam_u = ed.solve_dynamic(_eliashberg_input(flex_outdir))
    lam_ir = ed.solve_dynamic(_eliashberg_input(
        flex_outdir, extra={"matsubara_basis": "ir"}))
    diff_128 = abs(lam_ir - lam_u)
    assert diff_128 < 5e-2 * abs(lam_u)

    d = np.load(os.path.join(flex_outdir, "gap_dynamic.npz"))
    assert d["gap"].shape[-1] == NMAT
    assert str(d["matsubara_basis"]) == "ir"
    assert d["iomega"].shape == (NMAT,)

    # input-Nmat convergence: regenerate the pipeline at 4x Nmat
    out512 = str(tmp_path / "flex512")
    _run_flex_to(out512, nmat=4 * NMAT)
    inp_u = _eliashberg_input(out512)
    inp_u["mode"]["param"]["Nmat"] = 4 * NMAT
    inp_i = _eliashberg_input(out512, extra={"matsubara_basis": "ir"})
    inp_i["mode"]["param"]["Nmat"] = 4 * NMAT
    lam_u5 = ed.solve_dynamic(inp_u)
    lam_i5 = ed.solve_dynamic(inp_i)
    diff_512 = abs(lam_i5 - lam_u5)
    # measured 1.7e-3/0.398 with the dispersion-based auto wmax (5.4e-4
    # under the inflated pre-#57 estimate); the systematic-shrink assertion
    # below is the strong claim.
    assert diff_512 < 5e-3 * abs(lam_u5)
    assert diff_512 < 0.5 * diff_128, (
        "no Nmat convergence: {} vs {}".format(diff_512, diff_128))


def test_solve_dynamic_ir_requires_sparse_ir(flex_outdir, monkeypatch):
    import hwave.solver.ir_axis as ir_axis
    from hwave.solver import eliashberg_dynamic as ed
    monkeypatch.setattr(ir_axis, "_import_sparse_ir",
                        lambda: (_ for _ in ()).throw(ImportError("nope")))
    with pytest.raises(ImportError, match="sparse-ir"):
        ed.solve_dynamic(_eliashberg_input(
            flex_outdir, extra={"matsubara_basis": "ir"}))


def test_solve_dynamic_rejects_unknown_basis(flex_outdir):
    from hwave.solver import eliashberg_dynamic as ed
    with pytest.raises(ValueError, match="matsubara_basis"):
        ed.solve_dynamic(_eliashberg_input(
            flex_outdir, extra={"matsubara_basis": "nonsense"}))


def test_solve_dynamic_ir_gpu_matches_cpu(flex_outdir):
    """The IR kernel is xp-dispatched (transform matmuls + einsums); gpu=true
    must reproduce the CPU IR lambda."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    from hwave.solver import eliashberg_dynamic as ed
    lam_cpu = ed.solve_dynamic(_eliashberg_input(
        flex_outdir, extra={"matsubara_basis": "ir"}))
    lam_gpu = ed.solve_dynamic(_eliashberg_input(
        flex_outdir, extra={"matsubara_basis": "ir", "gpu": True}))
    assert abs(lam_gpu - lam_cpu) < 1e-8 * max(1.0, abs(lam_cpu))


# ---------------------------------------------------------------------------
# Issue #57: instantaneous (frequency-independent) vertex part, auto-wmax,
# and the drop_constant guard on static-dominated data
# ---------------------------------------------------------------------------

def _const_vertex_gate(norb, Nx, Ny, Nz, beta, nmat, wmax):
    """Kernel gate with a PURE frequency-independent vertex (chi terms = 0):
    the bare 0.5*(S+C) piece is a delta(tau) in imaginary time, which the
    bosonic IR basis cannot represent -- it must be handled analytically."""
    from hwave.solver import eliashberg_dynamic as ed
    from hwave.solver.ir_axis import IRAxis

    rng = np.random.default_rng(7)
    Vc = rng.standard_normal((norb,)*4 + (Nx, Ny, Nz)) * 0.3
    V_w = np.broadcast_to(Vc[..., None], Vc.shape + (nmat,)).astype(complex)
    axF = IRAxis(beta=beta, wmax=wmax, eps=1e-8, statistics="F")
    axB = IRAxis(beta=beta, wmax=wmax, eps=1e-8, statistics="B")

    kx = np.linspace(0, 2*np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2*np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2*np.pi, Nz, endpoint=False)
    ek = -0.3*2*(np.cos(kx)[:, None, None] + np.cos(ky)[None, :, None]
                 + np.cos(kz)[None, None, :])
    wF_u = (2*np.arange(nmat) + 1 - nmat) * np.pi / beta
    g_u = 1.0/(1j*wF_u[None, None, None, :] - ek[..., None])
    wF_n = axF.freq_n * np.pi / beta
    g_n = 1.0/(1j*wF_n[None, None, None, :] - ek[..., None])
    G2_u = np.zeros((norb,)*4 + (Nx, Ny, Nz, nmat), dtype=complex)
    G2_n = np.zeros((norb,)*4 + (Nx, Ny, Nz, axF.n_freq), dtype=complex)
    for a in range(norb):
        for b in range(norb):
            G2_u[a, b, a, b] = -g_u * g_u[..., ::-1] / beta
            G2_n[a, b, a, b] = -g_n * g_n[..., ::-1] / beta
    phi = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex)
    phi_u = np.broadcast_to(phi[..., None],
                            phi.shape + (nmat,)).astype(complex)
    phi_n = np.broadcast_to(phi[..., None],
                            phi.shape + (axF.n_freq,)).astype(complex)

    out_u = ed.eliashberg_kernel_dynamic(V_w, G2_u, phi_u.copy(), norb, beta)
    # IR path: the dynamic part of this vertex is ZERO; everything is
    # instantaneous and must flow through the analytic term.
    V_dyn_nodes = np.zeros(Vc.shape + (axB.n_freq,), dtype=complex)
    V_rt_tau = ed._ir_vertex_to_rtau(V_dyn_nodes, axB, axF)
    V_inst_rt = ed._spatial_ifftn(Vc.astype(complex), axes=(4, 5, 6))
    out_n = ed.eliashberg_kernel_ir(V_rt_tau, G2_n, phi_n.copy(), axF, beta,
                                    V_inst_rt=V_inst_rt)
    # The instantaneous term's output is nu-INDEPENDENT -- deliberately out
    # of the fermionic basis -- so compare on the NODES (where the solver
    # lives) by subsampling the uniform result at the node frequencies,
    # NOT via a fit/densify round trip (which would alias the constant and
    # test the comparison, not the kernel).
    k_idx = (axF.freq_n - 1 + nmat) // 2
    inside = (k_idx >= 0) & (k_idx < nmat)   # IR nodes can exceed the grid
    out_u_nodes = out_u[..., k_idx[inside]]
    diff = np.abs(out_n[..., inside] - out_u_nodes)
    return diff.max() / np.abs(out_u_nodes).max()


def test_ir_kernel_instantaneous_vertex_matches_uniform():
    # small-beta (fixture-like) and large-beta (production-like) regimes
    assert _const_vertex_gate(2, 4, 4, 1, beta=2.0, nmat=256, wmax=5.0) < 2e-2
    assert _const_vertex_gate(2, 1, 8, 8, beta=200.0, nmat=1024,
                              wmax=5.0) < 2e-2


@pytest.fixture(scope="module")
def flex_outdir_offsite(tmp_path_factory):
    """FLEX outputs for a model WITH off-site CoulombInter: there S+C != 0
    (unlike pure CoulombIntra, where the Kuroki S+C cancels), so the
    pairing vertex carries a genuine instantaneous part -- the model class
    issue #57 hit on beta'-(ET)2ICl2 (UVg)."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    outdir = str(tmp_path_factory.mktemp("flexout_uv"))
    info_mode = {'mode': 'FLEX', 'param': {
        'T': BETA_T, 'mu': 0.0, 'CellShape': [LX, LY, 1],
        'SubShape': [1, 1, 1], 'Nmat': NMAT,
        'IterationMax': 60, 'Mix': 0.3, 'EPS': 8},
        'calc_scheme': 'reduced'}
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat',
        'CoulombInter': 'coulombinter.dat'}}
    rio = read_input_k.QLMSkInput(info_input)
    ham = rio.get_param("ham")
    gi = rio.get_param("green")
    s = solver_flex.FLEX(ham, {}, info_mode)
    s.solve(gi, outdir)
    s.save_results({'path_to_output': outdir, 'chiq_s': 'chiq_s',
                    'chiq_c': 'chiq_c', 'green': 'green'}, gi)
    return outdir


def _offsite_input(outdir, extra=None):
    inp = _eliashberg_input(outdir, extra=extra)
    inp["file"]["input"]["interaction"]["CoulombInter"] = "coulombinter.dat"
    return inp


def test_ir_matvec_gate_offsite_interaction(flex_outdir_offsite):
    """End-to-end kernel gate on the off-site model, BOTH pairing channels
    (the triplet instantaneous part is 0.5*(C-S); a sign regression there
    must not hide behind the singlet-only gate)."""
    from hwave.solver import eliashberg_dynamic as ed
    import hwave.sc as sc

    inp = _offsite_input(flex_outdir_offsite)
    norb, Nx, Ny, Nz = 1, LX, LY, 1
    beta = 1.0 / BETA_T
    chis_w, chic_w, green_w, conv = ed.load_flex_chi_dynamic(
        inp, norb, Nx, Ny, Nz)
    kx = np.linspace(0, 2*np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2*np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2*np.pi, Nz, endpoint=False)
    geom, hr, inter = sc._read_interaction_files(inp)
    inter_k = sc._build_interaction_k(kx, ky, kz, inter, norb)

    axF, axB = ed._ir_axes_for_run(inp["eliashberg"], beta, hr, inter_k,
                                   norb)
    chis_n = ed._ir_compress(chis_w, axB, NMAT, "chiq_s", drop_constant=True)
    chic_n = ed._ir_compress(chic_w, axB, NMAT, "chiq_c", drop_constant=True)
    green_n = ed._ir_compress(green_w, axF, NMAT, "green")
    G2_n = ed.calc_g2_dynamic(green_n, beta)
    chis_d = axB.eval_to_uniform(axB.fit_from_freq(chis_n), NMAT)
    chic_d = axB.eval_to_uniform(axB.fit_from_freq(chic_n), NMAT)
    green_d = axF.eval_to_uniform(axF.fit_from_freq(green_n), NMAT)
    G2_w = ed.calc_g2_dynamic(green_d, beta)
    gap_static = sc._initialize_gap("cos", norb, kx, ky, kz)
    phi_u = np.broadcast_to(gap_static[..., None],
                            gap_static.shape + (NMAT,)).astype(complex)
    phi_n = np.broadcast_to(gap_static[..., None],
                            gap_static.shape + (axF.n_freq,)).astype(complex)

    for pairing in ("singlet", "triplet"):
        V_inst = ed._instantaneous_vertex(inter_k, norb, Nx, Ny, Nz,
                                          pairing_type=pairing,
                                          convention=conv)
        if pairing == "singlet":
            # 0.5*(S+C) != 0 with off-site V -- the issue-#57 case
            assert np.abs(V_inst).max() > 1e-10
        else:
            # 0.5*(C-S) on this 1-orbital fixture is exactly V(q): the on-site
            # U cancels between the two channels (Pauli -- two same-spin
            # electrons cannot share a site, so U cannot drive triplet
            # pairing), while the inter-site V survives, since it enters the
            # charge channel only.  This used to assert the vertex was
            # identically zero, which was the missing 2 V_aa(q) on the charge
            # diagonal (issue #95) showing up as "an extended Hubbard model
            # has no triplet pairing interaction".
            V_q = np.squeeze(inter_k["CoulombInter"][0, 0])
            assert np.abs(V_inst).max() > 1e-10
            np.testing.assert_allclose(
                np.squeeze(V_inst), V_q, rtol=0.0, atol=1e-12,
                err_msg="triplet instantaneous vertex must equal V(q)")
        V_n = ed.compute_vertices_flex_dynamic(chis_n, chic_n, inter_k,
                                               norb, Nx, Ny, Nz,
                                               pairing_type=pairing,
                                               convention=conv)
        V_rt_tau = ed._ir_vertex_to_rtau(V_n - V_inst[..., None], axB, axF)
        V_inst_rt = ed._spatial_ifftn(V_inst.astype(complex),
                                      axes=(4, 5, 6))
        V_w = ed.compute_vertices_flex_dynamic(chis_d, chic_d, inter_k,
                                               norb, Nx, Ny, Nz,
                                               pairing_type=pairing,
                                               convention=conv)
        out_u = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi_u.copy(),
                                             norb, beta)
        out_n = ed.eliashberg_kernel_ir(V_rt_tau, G2_n, phi_n.copy(), axF,
                                        beta, V_inst_rt=V_inst_rt)
        # node-subsampled comparison: the instantaneous part's output is
        # nu-independent (out of the fermionic basis on purpose), so a
        # fit/densify round trip would alias it -- compare where the
        # solver lives, on the nodes.
        k_idx = (axF.freq_n - 1 + NMAT) // 2
        inside = (k_idx >= 0) & (k_idx < NMAT)
        out_u_nodes = out_u[..., k_idx[inside]]
        diff = (np.abs(out_n[..., inside] - out_u_nodes).max()
                / np.abs(out_u_nodes).max())
        assert diff < 2e-2, (pairing, diff)


def test_solve_dynamic_ir_gpu_matches_cpu_offsite(flex_outdir_offsite):
    """gpu=true with a NONZERO instantaneous vertex: exercises the device
    handling of the V_inst split (the pure-CoulombIntra GPU test cannot --
    there V_inst cancels exactly)."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    from hwave.solver import eliashberg_dynamic as ed
    lam_cpu = ed.solve_dynamic(_offsite_input(
        flex_outdir_offsite, extra={"matsubara_basis": "ir"}))
    lam_gpu = ed.solve_dynamic(_offsite_input(
        flex_outdir_offsite, extra={"matsubara_basis": "ir", "gpu": True}))
    assert abs(lam_gpu - lam_cpu) < 1e-8 * max(1.0, abs(lam_cpu))





def _reverse_q_and_orbitals(V):
    """``V_{abcd}(q) -> V_{dcba}(-q)`` on the FFT grid.

    A pairing vertex INVARIANT under this map is exactly what the dynamic
    kernel needs in order to commute with the combined parity operator
    ``Delta_{ab}(k, iw) -> Delta_{ba}(-k, -iw)`` (the vertices the solver
    actually builds from a symmetric susceptibility have it by
    construction; a random tensor does not).
    """
    from hwave.solver.kgrid import reverse_fft_axes
    return reverse_fft_axes(V, (4, 5, 6)).transpose(3, 2, 1, 0, 4, 5, 6)


def _smooth_vertex_pieces(norb, Nx, Ny, Nz, beta, nmat, wmax, g=1.3,
                          symmetrize=False):
    """Ingredients of the kernel-algebra gate below: a fully IR-representable
    synthetic vertex (Lorentzian bosonic profile) on both frequency grids, the
    bare pair bubble on both, and a flat trial gap.

    ``symmetrize`` additionally projects the random momentum profile onto the
    parity-invariant vertices of :func:`_reverse_q_and_orbitals`, which the
    parity-commutation cases below need and the plain agreement gate does not.
    """
    from hwave.solver.ir_axis import IRAxis

    rng = np.random.default_rng(11)
    Vc = rng.standard_normal((norb,)*4 + (Nx, Ny, Nz)) * 0.4
    if symmetrize:
        Vc = 0.5 * (Vc + _reverse_q_and_orbitals(Vc))
    nuB_u = (2*np.arange(nmat) - nmat) * np.pi / beta
    V_w = (Vc[..., None] * (g*g/(nuB_u**2 + g*g))).astype(complex)
    axF = IRAxis(beta=beta, wmax=wmax, eps=1e-8, statistics="F")
    axB = IRAxis(beta=beta, wmax=wmax, eps=1e-8, statistics="B")
    nuB_n = axB.freq_n * np.pi / beta
    V_n = (Vc[..., None] * (g*g/(nuB_n**2 + g*g))).astype(complex)
    kx = np.linspace(0, 2*np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2*np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2*np.pi, Nz, endpoint=False)
    ek = -0.3*2*(np.cos(kx)[:, None, None] + np.cos(ky)[None, :, None]
                 + np.cos(kz)[None, None, :])
    wF_u = (2*np.arange(nmat) + 1 - nmat) * np.pi / beta
    g_u = 1.0/(1j*wF_u[None, None, None, :] - ek[..., None])
    wF_n = axF.freq_n * np.pi / beta
    g_n = 1.0/(1j*wF_n[None, None, None, :] - ek[..., None])
    G2_u = np.zeros((norb,)*4 + (Nx, Ny, Nz, nmat), dtype=complex)
    G2_n = np.zeros((norb,)*4 + (Nx, Ny, Nz, axF.n_freq), dtype=complex)
    for a in range(norb):
        for b in range(norb):
            G2_u[a, b, a, b] = -g_u * g_u[..., ::-1] / beta
            G2_n[a, b, a, b] = -g_n * g_n[..., ::-1] / beta
    phi = np.ones((norb, norb, Nx, Ny, Nz), dtype=complex)
    phi_u = np.broadcast_to(phi[..., None],
                            phi.shape + (nmat,)).astype(complex)
    phi_n = np.broadcast_to(phi[..., None],
                            phi.shape + (axF.n_freq,)).astype(complex)
    return dict(rng=rng, axF=axF, axB=axB, V_w=V_w, V_n=V_n, G2_u=G2_u,
                G2_n=G2_n, phi_u=phi_u, phi_n=phi_n, norb=norb, beta=beta,
                nmat=nmat, shape=(Nx, Ny, Nz))


def _smooth_vertex_gate(norb, Nx, Ny, Nz, beta, nmat, wmax, g=1.3):
    """Kernel-algebra exactness gate: a fully IR-representable synthetic
    vertex (Lorentzian bosonic profile) and bare pair bubble -- the uniform
    and IR kernels must agree to numerical precision (no chi-artifact
    contamination, unlike the FLEX-data gate)."""
    from hwave.solver import eliashberg_dynamic as ed

    p = _smooth_vertex_pieces(norb, Nx, Ny, Nz, beta, nmat, wmax, g=g)
    axF = p["axF"]
    out_u = ed.eliashberg_kernel_dynamic(p["V_w"], p["G2_u"],
                                         p["phi_u"].copy(), norb, beta)
    Vrt = ed._ir_vertex_to_rtau(p["V_n"], p["axB"], axF)
    out_n = ed.eliashberg_kernel_ir(Vrt, p["G2_n"], p["phi_n"].copy(), axF, beta)
    k_idx = (axF.freq_n - 1 + nmat) // 2
    inside = (k_idx >= 0) & (k_idx < nmat)
    return (np.abs(out_n[..., inside] - out_u[..., k_idx[inside]]).max()
            / np.abs(out_u).max())


def test_channel_decomposition_ir_offsite_bare_only_keeps_instantaneous_term(
    flex_outdir_offsite, monkeypatch
):
    """Both zero flags on the IR/off-site path remove fluctuations while the
    nonzero singlet bare vertex is subtracted and reinserted analytically."""
    from hwave.solver import eliashberg_dynamic as ed

    captured = {}
    real_kernel = ed.eliashberg_kernel_ir
    real_compress = ed._ir_compress

    def capture(vertex_rt, *args, **kwargs):
        captured["dynamic_max"] = float(np.max(np.abs(vertex_rt)))
        captured["instantaneous_max"] = float(np.max(np.abs(kwargs["V_inst_rt"])))
        return real_kernel(vertex_rt, *args, **kwargs)

    def reject_unused_channel_compression(arr, axis, nmat, label, **kwargs):
        if label in ("chiq_s", "chiq_c"):
            raise AssertionError("zeroed channel must bypass IR compression")
        return real_compress(arr, axis, nmat, label, **kwargs)

    monkeypatch.setattr(ed, "eliashberg_kernel_ir", capture)
    monkeypatch.setattr(ed, "_ir_compress", reject_unused_channel_compression)
    inp = _offsite_input(
        flex_outdir_offsite,
        extra={
            "matsubara_basis": "ir",
            "pairing_type": "singlet",
            "zero_chi_s": True,
            "zero_chi_c": True,
        },
    )
    ed.solve_dynamic(inp)
    assert captured["dynamic_max"] < 1e-10
    assert captured["instantaneous_max"] > 1e-8


@pytest.mark.parametrize("pairing_type", ["singlet", "triplet"])
@pytest.mark.parametrize(
    "flags, zeroed",
    [
        ({"zero_chi_s": True}, {"chiq_s"}),
        ({"zero_chi_c": True}, {"chiq_c"}),
        ({"zero_chi_s": True, "zero_chi_c": True}, {"chiq_s", "chiq_c"}),
    ],
)
def test_channel_decomposition_ir_bypasses_only_zeroed_channel(
    flex_outdir_offsite, monkeypatch, pairing_type, flags, zeroed
):
    """On the IR/off-site path a zeroed susceptibility channel (a) bypasses IR
    compression entirely and (b) reaches the vertex builder identically zero,
    while every retained channel is still compressed and nonzero -- for both
    pairing parities and for single- as well as both-channel zeroing (the
    both-zeroed case is covered above; here the single-channel routing is the
    point). (a) guards against a regression that densifies a discarded channel;
    (b) is the scientifically decisive check that the fluctuations are actually
    removed from the vertex, not merely routed differently."""
    from hwave.solver import eliashberg_dynamic as ed

    compressed = []
    captured = {}
    real_compress = ed._ir_compress
    real_vertex = ed.compute_vertices_flex_dynamic

    def record_compress(arr, axis, nmat, label, **kwargs):
        if label in ("chiq_s", "chiq_c"):
            compressed.append(label)
        return real_compress(arr, axis, nmat, label, **kwargs)

    def capture_vertex(chis_w, chic_w, *args, **kwargs):
        captured["chis_max"] = float(np.max(np.abs(chis_w)))
        captured["chic_max"] = float(np.max(np.abs(chic_w)))
        return real_vertex(chis_w, chic_w, *args, **kwargs)

    monkeypatch.setattr(ed, "_ir_compress", record_compress)
    monkeypatch.setattr(ed, "compute_vertices_flex_dynamic", capture_vertex)
    inp = _offsite_input(
        flex_outdir_offsite,
        extra=dict(flags, matsubara_basis="ir", pairing_type=pairing_type),
    )
    lam = ed.solve_dynamic(inp)

    # (a) discarded channels are never densified/compressed
    seen = set(compressed)
    assert seen.isdisjoint(zeroed), (
        "zeroed channel(s) {} must bypass IR compression, saw {}".format(
            zeroed, seen))
    retained = {"chiq_s", "chiq_c"} - zeroed
    assert retained <= seen, (
        "retained channel(s) {} must still be IR-compressed, saw {}".format(
            retained, seen))
    # (b) the vertex builder sees an exactly-zero discarded channel and a
    #     nonzero retained channel
    if "chiq_s" in zeroed:
        assert captured["chis_max"] == 0.0
    else:
        assert captured["chis_max"] > 0.0
    if "chiq_c" in zeroed:
        assert captured["chic_max"] == 0.0
    else:
        assert captured["chic_max"] > 0.0
    assert np.isfinite(lam)


class TestIRInstantaneousVertexTerm(unittest.TestCase):
    """The frequency-flat (instantaneous) vertex term of the IR kernel.

    With a frequency-INDEPENDENT vertex the kernel's frequency convolution
    collapses to ``-V_inst * (1/beta) sum_n F(i w_n)``. That Matsubara sum is
    UNREGULARIZED (no ``e^{i w 0^+}`` factor), so for an ``F`` with a jump at
    ``tau = 0`` it converges to the jump MIDPOINT ``0.5 (F(0^+) + F(0^-))``,
    not to the one-sided ``F(0^+)``. The surplus half-jump is the one piece of
    ``F`` that ANTI-commutes with the frequency reversal ``i w -> -i w``, so
    using ``F(0^+)`` (a) breaks the kernel's commutation with the combined
    parity operator, and (b) does not converge to the uniform-grid kernel's
    value at any ``Nmat``. Both are asserted below.

    These cases need a frequency-ODD component in the probe: the shipped
    agreement gates drive the kernel with a gap that is flat in frequency, and
    for a frequency-EVEN ``F`` the even-``l`` coefficients vanish, the two
    functionals coincide identically, and the defect is invisible.
    """

    def setUp(self):
        if not _HAVE_SPARSE_IR:
            raise unittest.SkipTest("sparse-ir not installed")

    def test_matsubara_sum_weights_are_the_truncated_uniform_sum_limit(self):
        """``IRAxis.u_matsubara_sum`` at the basis level, with no kernel in the
        way: for an IR-representable function the plain sum over a centered
        uniform fermionic grid converges to ``coeffs @ u_matsubara_sum`` as
        O(beta / Nmat) -- each doubling halves the gap -- and stays O(1) away
        from ``coeffs @ u_zero_plus`` forever. That is what identifies the
        midpoint, and not the one-sided limit, as the functional the uniform
        kernel's single tau bin is a truncation of."""
        from hwave.solver.ir_axis import IRAxis

        beta = 2.0
        ax = IRAxis(beta=beta, wmax=20.0, eps=1e-8, statistics="F")
        rng = np.random.default_rng(3)
        c = rng.standard_normal(ax.L) + 1j * rng.standard_normal(ax.L)
        mid = c @ ax.u_matsubara_sum
        one_sided = c @ ax.u_zero_plus
        self.assertGreater(abs(mid - one_sided), 1e-3 * abs(mid),
                           "the fixture has no tau jump: the two functionals "
                           "coincide and the case proves nothing")
        gaps = {}
        for nmat in (128, 256, 512, 1024):
            n = 2 * np.arange(nmat) + 1 - nmat
            uhat = np.asarray(ax._basis.uhat(n)).reshape(ax.L, nmat)
            total = (c @ uhat).sum() / beta
            gaps[nmat] = (abs(total - mid), abs(total - one_sided))
        for lo, hi in ((128, 256), (256, 512), (512, 1024)):
            # measured: 1.24 -> 0.655 -> 0.333 -> 0.168 against the midpoint,
            # and a flat 1.03e+01 against the one-sided value
            self.assertLess(gaps[hi][0], 0.6 * gaps[lo][0],
                            "no O(beta/Nmat) convergence to u_matsubara_sum "
                            "at Nmat {} -> {}".format(lo, hi))
            self.assertGreater(gaps[hi][1], 0.9 * gaps[lo][1],
                               "u_zero_plus must NOT be the limit of the "
                               "truncated sum (Nmat {} -> {})".format(lo, hi))

    def test_matsubara_sum_weights_follow_the_statistics(self):
        """The bosonic half of the test above. The extension is PERIODIC, not
        anti-periodic, so the tau = 0 midpoint is ``0.5 (f(0^+) + f(beta^-))``
        and the bosonic axis needs the PLUS sign the fermionic one does not.
        Against the fermionic combination the truncated bosonic sum does not
        converge at all."""
        from hwave.solver.ir_axis import IRAxis

        beta = 2.0
        ax = IRAxis(beta=beta, wmax=20.0, eps=1e-8, statistics="B")
        self.assertEqual(ax.statistics, "B")
        rng = np.random.default_rng(3)
        c = rng.standard_normal(ax.L) + 1j * rng.standard_normal(ax.L)
        mid = c @ ax.u_matsubara_sum
        fermionic = c @ (0.5 * (np.asarray(ax.u_zero_plus)
                                - np.asarray(ax.u_beta_minus)))
        self.assertGreater(abs(mid - fermionic), 1e-3 * abs(mid),
                           "the fixture has no tau jump: the two combinations "
                           "coincide and the case proves nothing")
        gaps = {}
        for nmat in (128, 256, 512, 1024):
            n = 2 * np.arange(nmat) - nmat              # a CENTRED bosonic grid
            uhat = np.asarray(ax._basis.uhat(n)).reshape(ax.L, nmat)
            total = (c @ uhat).sum() / beta
            gaps[nmat] = (abs(total - mid), abs(total - fermionic))
        for lo, hi in ((128, 256), (256, 512), (512, 1024)):
            # measured: 8.65 -> 4.54 -> 2.30 -> 1.15 against the bosonic
            # midpoint, and a nearly flat 1.20e+01 -> 9.77 against the
            # fermionic combination
            self.assertLess(gaps[hi][0], 0.6 * gaps[lo][0],
                            "no O(beta/Nmat) convergence to the bosonic "
                            "u_matsubara_sum at Nmat {} -> {}".format(lo, hi))
            self.assertGreater(gaps[hi][1], 0.7 * gaps[lo][1],
                               "the FERMIONIC combination must NOT be the limit "
                               "of the bosonic truncated sum (Nmat {} -> {})"
                               .format(lo, hi))

    @staticmethod
    def _instantaneous_vertex(p, scale=0.6, seed=29):
        """A parity-invariant, frequency-flat vertex on the ``(q)`` grid, plus
        its spatial transform (what the IR kernel consumes)."""
        from hwave.solver import eliashberg_dynamic as ed
        norb = p["norb"]
        rng = np.random.default_rng(seed)
        V = rng.standard_normal((norb,)*4 + p["shape"]) * scale
        V = 0.5 * (V + _reverse_q_and_orbitals(V))
        return V, ed._spatial_ifftn(V.astype(complex), axes=(4, 5, 6))

    def test_ir_kernel_with_instantaneous_vertex_commutes_with_parity(self):
        """Measured before the fix: 8.05e-02, against 8e-16 for the same IR
        kernel without the instantaneous term and 9e-16 for the uniform kernel
        with it. After: 2.0e-14."""
        from scipy.sparse.linalg import LinearOperator
        from hwave.solver import eliashberg_dynamic as ed

        p = _smooth_vertex_pieces(2, 4, 4, 1, beta=2.0, nmat=256, wmax=20.0,
                                  symmetrize=True)
        axF, beta, norb = p["axF"], p["beta"], p["norb"]
        V_inst, V_inst_rt = self._instantaneous_vertex(p)
        self.assertGreater(float(np.abs(V_inst_rt).max()), 0.0)
        Vrt = ed._ir_vertex_to_rtau(p["V_n"], p["axB"], axF)

        def operator(shape, matvec):
            n = int(np.prod(shape))
            return LinearOperator((n, n), matvec=matvec, dtype=complex)

        # control: the SAME vertex on the uniform kernel (where the flat term
        # is just one more bosonic frequency) commutes -- so the fixture is
        # admissible and any IR leakage is the IR flat term, not the vertex
        gap_u = (norb, norb) + p["shape"] + (p["nmat"],)
        V_w_tot = p["V_w"] + V_inst[..., np.newaxis].astype(complex)
        A_u = operator(gap_u, lambda v: ed.eliashberg_kernel_dynamic(
            V_w_tot, p["G2_u"], v.reshape(gap_u), norb, beta).ravel())
        gap_n = (norb, norb) + p["shape"] + (axF.n_freq,)
        A_n = operator(gap_n, lambda v: ed.eliashberg_kernel_ir(
            Vrt, p["G2_n"], v.reshape(gap_n), axF, beta,
            V_inst_rt=V_inst_rt).ravel())
        for eta in ("singlet", "triplet"):
            self.assertLessEqual(ed._parity_leakage(A_u, gap_u, eta), 1e-10,
                                 "uniform control, " + eta)
            self.assertLessEqual(ed._parity_leakage(A_n, gap_n, eta), 1e-10,
                                 "IR kernel, " + eta)

    def _flat_term_vs_uniform(self, nmat, norb=2, shape=(4, 4, 1), beta=2.0,
                              wmax=20.0, seed=17):
        """Max relative IR-vs-uniform difference of the FLAT term ALONE (the
        dynamic vertex is identically zero), on a probe whose ``F`` really
        carries a ``1/(i w)`` tail: a frequency-FLAT pair bubble times a
        single-pole gap. The physical bubble decays as ``1/w^2`` and would
        hide the tau jump the two functionals differ by."""
        from hwave.solver import eliashberg_dynamic as ed
        from hwave.solver.ir_axis import IRAxis

        Nx, Ny, Nz = shape
        rng = np.random.default_rng(seed)
        axF = IRAxis(beta=beta, wmax=wmax, eps=1e-8, statistics="F")
        V_inst = rng.standard_normal((norb,)*4 + shape) * 0.6
        V_inst_rt = ed._spatial_ifftn(V_inst.astype(complex), axes=(4, 5, 6))
        G2_k = rng.standard_normal((norb,)*4 + shape).astype(complex)
        ek = -0.6 * (np.cos(2*np.pi*np.arange(Nx)/Nx)[:, None, None]
                     + np.cos(2*np.pi*np.arange(Ny)/Ny)[None, :, None]
                     + np.cos(2*np.pi*np.arange(Nz)/Nz)[None, None, :])
        base = (rng.standard_normal((norb, norb) + shape)
                + 1j * rng.standard_normal((norb, norb) + shape))
        w_u = (2*np.arange(nmat) + 1 - nmat) * np.pi / beta
        w_n = axF.freq_n * np.pi / beta

        def gap(w):
            return (base[..., None]
                    / (1j*w - ek[None, None, ..., None])).astype(complex)

        G2_u = np.broadcast_to(G2_k[..., None],
                               G2_k.shape + (nmat,)).astype(complex)
        G2_n = np.broadcast_to(G2_k[..., None],
                               G2_k.shape + (axF.n_freq,)).astype(complex)
        V_u = np.broadcast_to(V_inst[..., None].astype(complex),
                              G2_k.shape + (nmat,)).copy()
        out_u = ed.eliashberg_kernel_dynamic(V_u, G2_u, gap(w_u), norb, beta)
        zero_rt = np.zeros((norb,)*4 + shape + (axF.n_tau,), dtype=complex)
        out_n = ed.eliashberg_kernel_ir(zero_rt, G2_n, gap(w_n), axF, beta,
                                        V_inst_rt=V_inst_rt)
        idx = (axF.freq_n - 1 + nmat) // 2
        inside = (idx >= 0) & (idx < nmat)
        return float(np.abs(out_n[..., inside] - out_u[..., idx[inside]]).max()
                     / np.abs(out_u).max())

    def test_ir_instantaneous_term_converges_to_uniform(self):
        """The uniform grid represents the flat term as a single tau bin, i.e.
        the Matsubara sum TRUNCATED at Nmat, so the two must differ by
        O(beta / Nmat): doubling Nmat halves the difference. Measured before
        the fix: 1.846e+00 -> 1.843e+00 (ratio 0.998, O(1) and non-convergent);
        after: 4.379e-03 -> 2.185e-03 (ratio 0.499)."""
        d128 = self._flat_term_vs_uniform(128)
        d256 = self._flat_term_vs_uniform(256)
        self.assertLess(d256, 0.6 * d128,
                        "no Nmat convergence of the flat term: "
                        "{:.3e} -> {:.3e}".format(d128, d256))
        self.assertGreater(d256, 0.4 * d128,
                           "the flat-term difference must HALVE (O(beta/Nmat)), "
                           "not collapse: {:.3e} -> {:.3e}".format(d128, d256))
        self.assertLess(d128, 1e-1,
                        "flat term off by O(1): {:.3e}".format(d128))


if __name__ == "__main__":
    unittest.main()
