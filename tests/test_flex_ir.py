"""IR-basis (sparse-ir) FLEX SCF (Stage 2 of docs/design/ir-matsubara.md).

Unlike Stage 1 (which consumes finite-Nmat uniform-FFT FLEX outputs), the IR
FLEX computes chi0 and Sigma natively on sparse nodes, so its agreement with
the uniform path is limited by the UNIFORM path's own O(beta/Nmat)
discretization artifacts: the equivalence tests therefore assert convergence
of the uniform result TOWARD the IR result as Nmat grows.

Import-safe under both pytest and unittest discovery (CI lesson from #54).
Tests must run from the repository root.
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


def _make_solver(nmat, extra_param=None, U=2.5, T=0.5, Lx=4, Ly=4,
                 iteration_max=60, mix=0.3):
    info_mode = {'mode': 'FLEX', 'param': {
        'T': T, 'mu': 0.0, 'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
        'Nmat': nmat, 'IterationMax': iteration_max, 'Mix': mix, 'EPS': 8}}
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
    if "CoulombIntra" in ham:
        ham["CoulombIntra"] = {k: complex(U) for k in ham["CoulombIntra"]}
    solver = solver_flex.FLEX(ham, {}, info_mode)
    return solver, gi


def _run(nmat, extra_param=None, **kw):
    solver, gi = _make_solver(nmat, extra_param, **kw)
    os.makedirs('tests/flex/output', exist_ok=True)
    solver.solve(gi, 'tests/flex/output')
    return solver, gi


def test_chi0_gate_ir_matches_uniform_large_nmat():
    """GATE (design Sec. 2 analogue): the IR-native chi0 on bosonic nodes
    must agree with the uniform-FFT chi0 computed at LARGE Nmat and
    compressed (constant dropped) onto the same nodes -- the uniform result
    converges to the IR one as its Nmat artifacts vanish."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 0.5
    beta = 1.0 / T

    diffs = []
    for nmat in (256, 1024):
        solver, gi = _make_solver(nmat, {'matsubara_basis': 'ir'})
        solver._calc_epsilon_k(gi)
        solver._ir_setup(beta)
        axF, axB = solver._ir_axF, solver._ir_axB

        # IR-native chi0 from the bare G on nodes
        g_nodes, _ = solver._calc_green_ir(beta, 0.0)
        chi0_ir = solver._calc_chi0q_ir(g_nodes, beta)[0]   # strip block

        # uniform chi0 from the same solver machinery at this Nmat
        solver_u, gi_u = _make_solver(nmat)
        solver_u._calc_epsilon_k(gi_u)
        g_u, gtail_u = solver_u._calc_green(beta, 0.0)
        chi0_u = solver_u._calc_chi0q(g_u, gtail_u, beta)[0]
        # (nmat, nvol, nd, nd) -> nodes, dropping the delta(tau) constant
        chi0_u_nodes = ed._ir_compress(
            np.moveaxis(chi0_u, 0, -1), axB, nmat, "chi0_u",
            drop_constant=True)
        chi0_u_nodes = np.moveaxis(chi0_u_nodes, -1, 0)

        scale = np.abs(chi0_u_nodes).max()
        diffs.append(np.abs(chi0_ir - chi0_u_nodes).max() / scale)

    # measured: 1.89e-2 (Nmat=256) -> 3.1e-3 (1024) -> 1.5e-3 (2048),
    # i.e. ~1/Nmat -- the uniform path's own artifact scale beyond the
    # removed delta(tau) constant (aliasing images).
    assert diffs[-1] < 5e-3, "chi0 gate failed: {}".format(diffs)
    assert diffs[-1] < 0.5 * diffs[0], \
        "no Nmat convergence toward IR chi0: {}".format(diffs)


def test_sigma_gate_one_iteration():
    """GATE: one SCF step (sigma=0 -> sigma_new) on the IR path must match
    the uniform step at large Nmat, compressed onto the fermionic nodes."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 0.5
    beta = 1.0 / T

    diffs = []
    for nmat in (256, 1024):
        s_ir, gi_ir = _run(nmat, {'matsubara_basis': 'ir'}, iteration_max=1,
                           mix=1.0)
        s_u, gi_u = _run(nmat, None, iteration_max=1, mix=1.0)
        axF = s_ir._ir_axF
        sig_u = np.moveaxis(gi_u["sigma"], 1, -1)       # freq axis last
        sig_u_nodes = np.moveaxis(
            ed._ir_compress(sig_u, axF, nmat, "sigma_u"), -1, 1)
        # IR run stores the DENSIFIED sigma; re-compress to nodes for the
        # node-space comparison (round trip is eps-exact)
        sig_ir = np.moveaxis(gi_ir["sigma"], 1, -1)
        sig_ir_nodes = np.moveaxis(
            ed._ir_compress(sig_ir, axF, nmat, "sigma_ir"), -1, 1)
        scale = np.abs(sig_u_nodes).max()
        diffs.append(np.abs(sig_ir_nodes - sig_u_nodes).max() / scale)

    assert diffs[-1] < 1e-2, "sigma gate failed: {}".format(diffs)
    assert diffs[-1] < 0.5 * diffs[0], \
        "no Nmat convergence toward IR sigma: {}".format(diffs)


def test_e2e_converged_flex_ir_vs_uniform():
    """Full SCF: the uniform observables converge toward the IR ones as
    Nmat grows (chi_s peak, sigma, and the fixed-mu particle number)."""
    s_ir, gi_ir = _run(1024, {'matsubara_basis': 'ir'})
    assert s_ir.scf_converged

    gaps = []
    for nmat in (256, 1024):
        s_u, gi_u = _run(nmat)
        assert s_u.scf_converged
        # chi_s static (nu=0 slice at the bosonic center) peak value
        chis_u = gi_u["chiq_s"][nmat // 2].real.max()
        chis_ir = gi_ir["chiq_s"][1024 // 2].real.max()
        n_u = gi_u["physics"]["NCond"]
        n_ir = gi_ir["physics"]["NCond"]
        gaps.append((abs(chis_u - chis_ir) / abs(chis_ir),
                     abs(n_u - n_ir) / max(abs(n_ir), 1e-30)))

    assert gaps[-1][0] < 2e-3, "chi_s peak mismatch: {}".format(gaps)
    assert gaps[-1][1] < 2e-3, "NCond mismatch: {}".format(gaps)
    assert gaps[-1][0] <= gaps[0][0] * 1.05, \
        "chi_s not converging toward IR: {}".format(gaps)


def test_ir_filling_driven_mu_conserves_ncond():
    """The IR dressed mu search (n from G(beta^-) through the basis) must
    deliver the requested particle number after the SCF."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    Ncond = 8.0     # Nstate = 4*4*2 = 32 -> doped
    info_mode = {'mode': 'FLEX', 'param': {
        'T': 0.5, 'Ncond': Ncond, 'CellShape': [4, 4, 1],
        'SubShape': [1, 1, 1], 'Nmat': 256, 'IterationMax': 60,
        'Mix': 0.3, 'EPS': 8, 'matsubara_basis': 'ir'},
        'calc_scheme': 'reduced'}
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'}}
    rio = read_input_k.QLMSkInput(info_input)
    ham = rio.get_param("ham")
    gi = rio.get_param("green")
    ham["CoulombIntra"] = {k: complex(2.5) for k in ham["CoulombIntra"]}
    s = solver_flex.FLEX(ham, {}, info_mode)
    os.makedirs('tests/flex/output', exist_ok=True)
    s.solve(gi, 'tests/flex/output')
    assert s.scf_converged
    assert abs(gi["physics"]["NCond"] - Ncond) < 1e-6


def test_ir_outputs_are_densified_uniform_grid():
    """IR FLEX outputs must land on the run's uniform Nmat grid so
    downstream consumers (incl. the Stage-1 dynamic Eliashberg loader)
    work unchanged."""
    nmat = 256
    s, gi = _run(nmat, {'matsubara_basis': 'ir'}, iteration_max=5)
    assert gi["sigma"].shape[1] == nmat
    assert gi["green"].shape[1] == nmat
    assert gi["chiq_s"].shape[0] == nmat
    assert gi["chi0q"].shape[0] == nmat
    for key in ("sigma", "green", "chiq_s", "chiq_c", "chi0q"):
        assert isinstance(gi[key], np.ndarray)


def test_ir_rejects_general_scheme():
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    info_mode = {'mode': 'FLEX', 'param': {
        'T': 0.5, 'mu': 0.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
        'Nmat': 64, 'matsubara_basis': 'ir'}, 'calc_scheme': 'general'}
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'}}
    rio = read_input_k.QLMSkInput(info_input)
    with pytest.raises(ValueError, match="matsubara_basis"):
        solver_flex.FLEX(rio.get_param("ham"), {}, info_mode)


def test_ir_rejects_unknown_basis():
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    info_mode = {'mode': 'FLEX', 'param': {
        'T': 0.5, 'mu': 0.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
        'Nmat': 64, 'matsubara_basis': 'nonsense'}, 'calc_scheme': 'reduced'}
    info_input = {'path_to_input': 'tests/rpa/input', 'interaction': {
        'path_to_input': 'tests/rpa/input', 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'}}
    rio = read_input_k.QLMSkInput(info_input)
    with pytest.raises(ValueError, match="matsubara_basis"):
        solver_flex.FLEX(rio.get_param("ham"), {}, info_mode)


def test_ir_with_anderson_mixing():
    """IR path composes with Anderson acceleration and reaches the same
    fixed point as IR + linear mixing."""
    s_lin, gi_lin = _run(256, {'matsubara_basis': 'ir'})
    s_and, gi_and = _run(256, {'matsubara_basis': 'ir',
                               'mixing_scheme': 'anderson'})
    assert s_and.scf_converged
    ref = gi_lin["sigma"]
    np.testing.assert_allclose(gi_and["sigma"], ref,
                               atol=1e-5 * np.abs(ref).max())
    assert s_and.scf_iterations < s_lin.scf_iterations


def test_ir_gpu_matches_cpu():
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    s_c, gi_c = _run(256, {'matsubara_basis': 'ir'}, iteration_max=10)
    s_g, gi_g = _run(256, {'matsubara_basis': 'ir', 'gpu': True},
                     iteration_max=10)
    np.testing.assert_allclose(gi_g["sigma"], gi_c["sigma"], atol=1e-9)
    assert isinstance(gi_g["sigma"], np.ndarray)
