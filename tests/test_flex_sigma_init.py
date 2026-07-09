"""FLEX self-energy warm-start (``sigma_init``).

FLEX starts its SCF loop from Sigma = 0 by default. Passing
``[file.input] sigma_init = "sigma.npz"`` seeds the loop from a previously
saved self-energy instead -- decisive for convergence near a magnetic
instability, where the Sigma=0 transient makes the SCF oscillate. These tests
pin: read_init loads it, solve warm-starts from it (a converged Sigma is a
fixed point that one iteration reproduces), and a grid mismatch fails fast.
"""
import os
import numpy as np
import pytest


def _make_solver(fixing, Lx=4, Ly=4, Nmat=32, T=4.0, U=1.0,
                 iteration_max=50, mix=0.2):
    """Weak-coupling FLEX solver (converges) in fixed-mu or fixed-N mode."""
    info_log = {}
    param = {'T': T, 'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
             'Nmat': Nmat, 'IterationMax': iteration_max, 'Mix': mix, 'EPS': 6}
    param.update(fixing)
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': 'reduced'}
    info_file_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {'path_to_input': 'tests/rpa/input',
                        'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat'},
    }
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    if "CoulombIntra" in ham_info:
        ham_info["CoulombIntra"] = {k: complex(U) for k in ham_info["CoulombIntra"]}
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    return solver, green_info


def _converged_sigma():
    solver, gi = _make_solver({'mu': 0.0}, iteration_max=60)
    solver.solve(gi, 'tests/flex/output')
    sigma = gi["sigma"]
    assert sigma is not None
    return sigma


def test_read_init_loads_sigma_init(tmp_path):
    sigma = _converged_sigma()
    np.savez(tmp_path / "sigma.npz", sigma=sigma)
    solver, _ = _make_solver({'mu': 0.0}, iteration_max=1)
    info = solver.read_init({'path_to_input': str(tmp_path),
                             'sigma_init': 'sigma.npz'})
    assert "sigma_init" in info
    assert info["sigma_init"].shape == sigma.shape


def test_warm_start_from_converged_is_a_fixed_point():
    """A converged self-energy is a fixed point: warm-starting a fresh run from
    it and taking ONE iteration reproduces it, whereas the Sigma=0 start after
    one iteration is far from it (so the check is non-trivial)."""
    sigma_star = _converged_sigma()

    # warm-started single iteration -> should stay at the fixed point
    solver_w, gi_w = _make_solver({'mu': 0.0}, iteration_max=1)
    gi_w["sigma_init"] = sigma_star
    solver_w.solve(gi_w, 'tests/flex/output')
    sig_w = gi_w["sigma"]
    rel_w = np.linalg.norm(sig_w - sigma_star) / np.linalg.norm(sigma_star)

    # zero-started single iteration -> far from the fixed point
    solver_z, gi_z = _make_solver({'mu': 0.0}, iteration_max=1)
    solver_z.solve(gi_z, 'tests/flex/output')
    sig_z = gi_z["sigma"]
    rel_z = np.linalg.norm(sig_z - sigma_star) / np.linalg.norm(sigma_star)

    assert rel_w < 1.0e-2, "warm-start did not reproduce the fixed point (rel={:.3e})".format(rel_w)
    assert rel_z > 10 * rel_w, "test trivial: zero-start too close (rel_z={:.3e}, rel_w={:.3e})".format(rel_z, rel_w)


def test_sigma_init_shape_mismatch_raises():
    solver, gi = _make_solver({'mu': 0.0}, iteration_max=1)
    gi["sigma_init"] = np.zeros((1, 8, 16, 2, 2), dtype=complex)  # wrong Nmat
    with pytest.raises(ValueError, match="sigma_init shape"):
        solver.solve(gi, 'tests/flex/output')
