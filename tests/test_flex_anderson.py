"""Tests for Anderson-accelerated self-energy mixing in FLEX.

The default SCF update is linear mixing
``sigma <- (1-Mix)*sigma + Mix*sigma_new``. ``mixing_scheme = "anderson"``
turns on Anderson acceleration (Pulay/DIIS-type extrapolation over a short
history of iterates and residuals), which reaches the same fixed point in
fewer iterations. Tests must run from the repository root.
"""
import os

import numpy as np
import pytest


def _run_flex(mixing_scheme=None, depth=None, gpu=False, U=2.5, T=1.0,
              Lx=8, Ly=8, Nmat=64, mix=0.2, iteration_max=200):
    info_log = {}
    param = {'T': T, 'mu': 0.0, 'CellShape': [Lx, Ly, 1],
             'SubShape': [1, 1, 1], 'Nmat': Nmat,
             'IterationMax': iteration_max, 'Mix': mix, 'EPS': 8,
             'gpu': gpu}
    if mixing_scheme is not None:
        param['mixing_scheme'] = mixing_scheme
    if depth is not None:
        param['anderson_depth'] = depth
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': 'reduced'}
    info_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {'path_to_input': 'tests/rpa/input',
                        'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat'},
    }
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    read_io = read_input_k.QLMSkInput(info_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    if "CoulombIntra" in ham_info:
        ham_info["CoulombIntra"] = {k: complex(U)
                                    for k in ham_info["CoulombIntra"]}
    os.makedirs('tests/flex/output', exist_ok=True)
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    solver.solve(green_info, 'tests/flex/output')
    return solver, green_info


def test_linear_run_reports_scf_iterations():
    """The solver must expose how the SCF ended (count + converged flag) so
    acceleration schemes can be compared and sweeps can inspect runs."""
    solver, _ = _run_flex()
    assert solver.scf_converged is True
    assert 1 <= solver.scf_iterations <= 200


def test_anderson_reaches_same_fixed_point_in_fewer_iterations():
    """Anderson mixing must converge to the SAME self-energy as linear mixing
    (the fixed point is a property of the equations, not the mixer) and take
    fewer SCF iterations at the same Mix."""
    sol_lin, gi_lin = _run_flex(mixing_scheme="linear")
    sol_and, gi_and = _run_flex(mixing_scheme="anderson")

    assert sol_lin.scf_converged and sol_and.scf_converged
    ref = gi_lin["sigma"]
    np.testing.assert_allclose(gi_and["sigma"], ref,
                               atol=1e-5 * np.abs(ref).max())
    assert sol_and.scf_iterations < sol_lin.scf_iterations, (
        "anderson ({}) not faster than linear ({})".format(
            sol_and.scf_iterations, sol_lin.scf_iterations))


def test_anderson_depth_one_matches_linear():
    """With no usable history (depth=1) the Anderson update degenerates to
    exactly the linear mixing step, so the two runs must agree closely."""
    sol_lin, gi_lin = _run_flex(mixing_scheme="linear", iteration_max=6)
    sol_and, gi_and = _run_flex(mixing_scheme="anderson", depth=1,
                                iteration_max=6)
    np.testing.assert_allclose(gi_and["sigma"], gi_lin["sigma"],
                               rtol=1e-10, atol=1e-12)


def test_unknown_mixing_scheme_raises():
    with pytest.raises(ValueError, match="mixing_scheme"):
        _run_flex(mixing_scheme="nonsense", iteration_max=2)


def test_anderson_gpu_matches_cpu():
    """Anderson mixing under gpu=true must reproduce the CPU fixed point."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    sol_c, gi_c = _run_flex(mixing_scheme="anderson", Lx=4, Ly=4, Nmat=32)
    sol_g, gi_g = _run_flex(mixing_scheme="anderson", Lx=4, Ly=4, Nmat=32,
                            gpu=True)
    assert sol_g.scf_converged
    np.testing.assert_allclose(gi_g["sigma"], gi_c["sigma"], atol=1e-8)
