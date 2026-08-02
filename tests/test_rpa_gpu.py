"""GPU (CuPy) execution tests for the RPA solver.

The CuPy equivalence tests skip on machines without a usable CUDA device; the
fallback test runs everywhere. Tests must run from the repository root (they
use the shared fixtures under tests/rpa/input).
"""
import logging
import os

import numpy as np
import pytest


def _run_rpa(gpu=False, calc_type="ring", Lx=4, Ly=4, Nmat=32, T=2.0, mu=0.0,
             coeff_tail=0.0):
    """Run a small 1-orbital Hubbard RPA solve; return green_info."""
    info_log = {}
    info_mode = {
        'mode': 'RPA',
        'param': {
            'T': T, 'mu': mu,
            'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1],
            'Nmat': Nmat,
            'coeff_tail': coeff_tail,
            'gpu': gpu,
        },
        'calc_scheme': 'general' if calc_type == "ring+ladder" else 'reduced',
        'calc_type': calc_type,
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
    os.makedirs('tests/rpa/output', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as solver_rpa
    read_io = read_input_k.QLMSkInput(info_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    solver = solver_rpa.RPA(ham_info, info_log, info_mode)
    solver.solve(green_info, 'tests/rpa/output')
    return green_info


def _skip_without_device():
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")


def test_rpa_gpu_falls_back_without_cupy(monkeypatch, caplog):
    """[mode.param] gpu=true without CuPy must warn and produce the identical
    (numpy-path) result."""
    from hwave.solver import backend

    ref = _run_rpa(gpu=False)

    def _no_cupy():
        raise ImportError("No module named 'cupy'")

    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        out = _run_rpa(gpu=True)
    assert any("cupy" in rec.message.lower() for rec in caplog.records)
    np.testing.assert_allclose(out["chi0q"], ref["chi0q"], atol=1e-12)
    np.testing.assert_allclose(out["chiq"], ref["chiq"], atol=1e-12)


@pytest.mark.parametrize("coeff_tail", [0.0, 1.0])
def test_rpa_gpu_matches_cpu(coeff_tail):
    """gpu=true on a real CUDA device must reproduce the CPU chi0q/chiq to
    fp64 round-off, stored as host numpy arrays. coeff_tail=0 covers the
    inactive-tail branch (the prior baseline); coeff_tail=1 exercises the
    active-tail endpoint branch (issue #134) on the device path."""
    _skip_without_device()
    ref = _run_rpa(gpu=False, coeff_tail=coeff_tail)
    out = _run_rpa(gpu=True, coeff_tail=coeff_tail)
    for key in ("chi0q", "chiq"):
        assert isinstance(out[key], np.ndarray), \
            "green_info['{}'] must be a host numpy array".format(key)
        np.testing.assert_allclose(out[key], ref[key], atol=1e-10,
                                   err_msg="RPA '{}' mismatch".format(key))


@pytest.mark.parametrize("coeff_tail", [0.0, 1.0])
def test_rpa_gpu_ladder_matches_cpu(coeff_tail):
    """The transverse (ring+ladder) channel must also match on the GPU,
    with the endpoint-correction branch both inactive and active."""
    _skip_without_device()
    ref = _run_rpa(gpu=False, calc_type="ring+ladder", coeff_tail=coeff_tail)
    out = _run_rpa(gpu=True, calc_type="ring+ladder", coeff_tail=coeff_tail)
    np.testing.assert_allclose(out["chiq"], ref["chiq"], atol=1e-10)
    np.testing.assert_allclose(out["chiq_pm"], ref["chiq_pm"], atol=1e-10)
