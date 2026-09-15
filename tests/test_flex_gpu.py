"""GPU (CuPy) execution and spatial-FFT parallelism tests for the FLEX solver.

The CuPy equivalence test skips on machines without a usable CUDA device; the
fallback and fft_workers tests run everywhere. Tests must run from the
repository root (they use the shared fixtures under tests/rpa/input).
"""
import logging
import os

import numpy as np
import pytest


def _run_flex(gpu=False, fft_workers=1, Lx=4, Ly=4, Nmat=32, T=2.0, mu=0.0,
              filling=None, return_solver=False):
    """Run a small 1-orbital Hubbard FLEX solve; return green_info."""
    info_log = {}
    param = {
        'T': T,
        'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
        'Nmat': Nmat,
        'gpu': gpu,
        'fft_workers': fft_workers,
    }
    if filling is not None:
        param['filling'] = filling
    else:
        param['mu'] = mu
    info_mode = {
        'mode': 'FLEX',
        'param': param,
        'calc_scheme': 'reduced',
    }
    info_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {
            'path_to_input': 'tests/rpa/input',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }
    os.makedirs('tests/flex/output', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    read_io = read_input_k.QLMSkInput(info_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    solver.solve(green_info, 'tests/flex/output')
    if return_solver:
        return green_info, solver
    return green_info


def _assert_results_close(a, b, atol):
    for key in ("sigma", "green", "chiq_s", "chiq_c", "chi0q"):
        np.testing.assert_allclose(
            a[key], b[key], atol=atol,
            err_msg="FLEX output '{}' mismatch".format(key))
    assert np.isclose(a["physics"]["NCond"], b["physics"]["NCond"],
                      atol=atol)
    assert np.isclose(a["physics"]["mu"], b["physics"]["mu"], atol=atol)


def test_flex_gpu_falls_back_without_cupy(monkeypatch, caplog):
    """[mode.param] gpu=true without CuPy must warn and produce the identical
    (numpy-path) result."""
    from hwave.solver import backend

    ref = _run_flex(gpu=False)

    def _no_cupy():
        raise ImportError("No module named 'cupy'")

    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        out = _run_flex(gpu=True)
    assert any("cupy" in rec.message.lower() for rec in caplog.records)
    _assert_results_close(out, ref, atol=1e-12)


def test_flex_gpu_matches_cpu():
    """gpu=true on a real CUDA device must reproduce the CPU FLEX results to
    fp64 round-off, including the outputs stored in green_info (which must be
    host numpy arrays)."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")

    ref = _run_flex(gpu=False)
    out = _run_flex(gpu=True)
    for key in ("sigma", "green", "chiq_s", "chiq_c", "chi0q"):
        assert isinstance(out[key], np.ndarray), \
            "green_info['{}'] must be a host numpy array".format(key)
    _assert_results_close(out, ref, atol=1e-10)


def test_flex_gpu_calc_mu_matches_cpu():
    """gpu=true with a filling target (calc_mu) must reproduce the CPU
    results: this exercises the host/device boundary of the dressed
    chemical-potential search (_matsubara_number_operator brings its
    non-Hermitian operator to the host). Also checks that no CuPy arrays
    leak into public solver attributes or ham_info after the solve."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")

    ref, _ = _run_flex(gpu=False, filling=0.4, return_solver=True)
    out, solver = _run_flex(gpu=True, filling=0.4, return_solver=True)
    _assert_results_close(out, ref, atol=1e-9)
    assert np.isclose(out["physics"]["Sz"], ref["physics"]["Sz"], atol=1e-9)
    for name in ("H0_eigenvalue", "H0_eigenvector", "green0", "green0_tail",
                 "sigma", "green_kw", "chi_s", "chi_c"):
        assert isinstance(getattr(solver, name), np.ndarray), \
            "solver.{} must be a host numpy array after solve".format(name)
    assert isinstance(solver.ham_info.ham_inter_q, np.ndarray), \
        "ham_info.ham_inter_q must not be mutated to a device array"


def test_flex_fft_workers_matches_serial():
    """The scipy-parallel spatial-FFT path must match the serial numpy path
    to machine precision through a full FLEX solve."""
    from hwave.solver import backend
    if backend._SFFT is None:
        pytest.skip("scipy.fft unavailable")
    serial = _run_flex(fft_workers=1)
    par = _run_flex(fft_workers=-1)
    _assert_results_close(par, serial, atol=1e-11)


# --- #181 follow-up: the local second-order kernel on the GPU -------------

def _run_flex_second_order_local(gpu=False, Nmat=32, iteration_max=3):
    """A general-scheme FLEX solve of the 2-orbital ON-SITE (U, U', Hund) plus
    OFF-SITE (V) input under ``flex_second_order = "local"``; returns
    ``(green_info, solver)``. The off-site rows are what make ``factors.vpair``
    non-None, so the device check below covers the off-site operand too."""
    import shutil
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    from tests.test_flex_second_order_scf import _write_inputs
    d = tempfile.mkdtemp(prefix="hwave_gpu_so_")
    try:
        idict = _write_inputs(d)
        read_io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        param = {'T': 1.0, 'filling': 0.5, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                 'Nmat': Nmat, 'IterationMax': iteration_max, 'Mix': 0.5, 'EPS': 8,
                 'flex_second_order': 'local', 'gpu': gpu}
        info_mode = {'mode': 'FLEX', 'param': param, 'enable_spin_orbital': False,
                     'calc_scheme': 'general'}
        solver = solver_flex.FLEX(read_io.get_param("ham"), {}, info_mode)
        green_info = read_io.get_param("green")
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return green_info, solver


def _operand_checking_accumulate_batch(module_prefix, seen):
    """A drop-in for :func:`hwave.solver.second_order.accumulate_batch` that
    asserts EVERY array operand lives on the expected backend before doing the
    real work: the output view, the bubble batch, both factor matrices of every
    spin triple, and the off-site ``vpair``.

    ``module_prefix`` is ``"cupy"`` for a device run and ``"numpy"`` for a host
    run -- the host spelling is what makes this wrapper itself testable without
    a CUDA device (see the CPU twin below), since the only difference between
    the two is the string.
    """
    import hwave.solver.second_order as so_mod
    orig = so_mod.accumulate_batch          # captured NOW: build every wrapper
                                            # before the first monkeypatch, or a
                                            # second one would wrap the first.

    def _check(name, x):
        assert type(x).__module__.split(".")[0] == module_prefix, \
            "{} is a {} array, expected {}".format(name, type(x).__module__, module_prefix)

    def wrapper(out_b, chibar_b, l0, factors, work=None):
        _check("out_b", out_b)
        _check("chibar_b", chibar_b)
        for i, (A, B) in enumerate(zip(factors.A_on, factors.B_on)):
            _check("factors.A_on[{}]".format(i), A)
            _check("factors.B_on[{}]".format(i), B)
        assert factors.vpair is not None, "the off-site fixture must produce a vpair"
        _check("factors.vpair", factors.vpair)
        seen["calls"] += 1
        return orig(out_b, chibar_b, l0, factors, work=work)

    return wrapper


def test_flex_second_order_local_operands_follow_the_backend_cpu(monkeypatch):
    """CPU twin of the GPU case below: the very same operand-checking wrapper,
    with the expected backend spelled ``"numpy"``, must see every operand of
    every ``accumulate_batch`` call on a host run. This is what makes the
    device assertion above meaningful -- it pins that the wrapper reaches the
    production call site, reads the right operands and really would fire, on a
    machine with no CUDA device."""
    import hwave.solver.second_order as so_mod
    seen, seen_device = {"calls": 0}, {"calls": 0}
    host_wrapper = _operand_checking_accumulate_batch("numpy", seen)
    device_wrapper = _operand_checking_accumulate_batch("cupy", seen_device)
    monkeypatch.setattr(so_mod, "accumulate_batch", host_wrapper)
    _run_flex_second_order_local(gpu=False)
    assert seen["calls"] > 0, "accumulate_batch was never called"
    # and the wrapper is not vacuous: asking for the device backend on this
    # host run must fail on the first operand.
    monkeypatch.setattr(so_mod, "accumulate_batch", device_wrapper)
    with pytest.raises(AssertionError, match="expected cupy"):
        _run_flex_second_order_local(gpu=False)
    assert seen_device["calls"] == 0


def test_flex_gpu_second_order_local_matches_cpu(monkeypatch):
    """``flex_second_order = "local"`` on a real CUDA device must reproduce the
    CPU self-energy to fp64 round-off, and every operand of the second-order
    kernel must be a device array (no silent host round trip of the factor
    pack, the bubble batch or the off-site vpair)."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")

    ref, _ = _run_flex_second_order_local(gpu=False)
    import hwave.solver.second_order as so_mod
    seen = {"calls": 0}
    monkeypatch.setattr(so_mod, "accumulate_batch",
                        _operand_checking_accumulate_batch("cupy", seen))
    out, solver = _run_flex_second_order_local(gpu=True)
    assert seen["calls"] > 0, "accumulate_batch was never called on the device run"
    assert solver.flex_second_order == "local"
    assert isinstance(out["sigma"], np.ndarray), "green_info['sigma'] must be a host array"
    scale = np.abs(ref["sigma"]).max()
    assert scale > 1e-6                                   # anti-vacuity
    np.testing.assert_allclose(out["sigma"], ref["sigma"], rtol=0, atol=1e-10 * scale,
                               err_msg="gpu/cpu sigma mismatch under the local kernel")
