"""GPU (CuPy) execution and spatial-FFT parallelism tests for the FLEX solver.

The CuPy equivalence test skips on machines without a usable CUDA device; the
fallback and fft_workers tests run everywhere. Tests must run from the
repository root (they use the shared fixtures under tests/rpa/input).
"""
import logging
import os

import numpy as np
import pytest


def _run_flex(gpu=False, fft_workers=-1, Lx=4, Ly=4, Nmat=32, T=2.0, mu=0.0):
    """Run a small 1-orbital Hubbard FLEX solve; return green_info."""
    info_log = {}
    info_mode = {
        'mode': 'FLEX',
        'param': {
            'T': T, 'mu': mu,
            'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'gpu': gpu,
            'fft_workers': fft_workers,
        },
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


def test_flex_fft_workers_matches_serial():
    """The scipy-parallel spatial-FFT path must match the serial numpy path
    to machine precision through a full FLEX solve."""
    from hwave.solver import backend
    if backend._SFFT is None:
        pytest.skip("scipy.fft unavailable")
    serial = _run_flex(fft_workers=1)
    par = _run_flex(fft_workers=-1)
    _assert_results_close(par, serial, atol=1e-11)
