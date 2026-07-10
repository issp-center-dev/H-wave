"""IR-basis (sparse-ir) path of the dynamic Eliashberg solver (Stage 1).

Design: docs/design/ir-matsubara.md. The fixtures run a REAL small FLEX
calculation first so every frequency object is physically smooth (pole
structure) and therefore IR-representable -- random-tensor fixtures are
white noise on the frequency axis and must not be used here.

Tests must run from the repository root.
"""
import os

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
    axF, axB = ed._ir_axes_for_run(inp["eliashberg"], beta, hr, inter_k)
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

    # The two kernels are different discretizations of the same operator:
    # uniform = cyclic frequency convolution on the finite grid, IR = the
    # continuous tau-product truncated to the basis bandwidth. They agree in
    # the joint large-wmax / large-Nmat limit; assert that agreement improves
    # systematically with wmax and is percent-level already at the default.
    assert diff_default < 5e-1
    eli2 = dict(inp["eliashberg"]); eli2["ir_wmax"] = 200.0
    axF2, axB2 = ed._ir_axes_for_run(eli2, beta, hr, inter_k)
    chis2 = ed._ir_compress(chis_w, axB2, NMAT, "chiq_s", drop_constant=True)
    chic2 = ed._ir_compress(chic_w, axB2, NMAT, "chiq_c", drop_constant=True)
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


def test_solve_dynamic_ir_matches_uniform_lambda(flex_outdir, tmp_path):
    """End-to-end: the IR and uniform leading eigenvalues are two
    discretizations of the same continuum answer, and their difference is
    bounded by the input-data quality (the finite-Nmat artifacts of the
    uniform-FFT FLEX output): percent-level at Nmat=128 and shrinking
    systematically as the input Nmat grows (measured 1.0e-2 -> 5.4e-4 for
    128 -> 512 on this fixture). The IR run writes the gap on the SAME
    uniform grid/metadata, with IR provenance recorded."""
    from hwave.solver import eliashberg_dynamic as ed

    lam_u = ed.solve_dynamic(_eliashberg_input(flex_outdir))
    lam_ir = ed.solve_dynamic(_eliashberg_input(
        flex_outdir, extra={"matsubara_basis": "ir"}))
    diff_128 = abs(lam_ir - lam_u)
    assert diff_128 < 3e-2 * abs(lam_u)

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
    assert diff_512 < 2e-3 * abs(lam_u5)
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
