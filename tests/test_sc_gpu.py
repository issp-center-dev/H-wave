"""GPU (CuPy) support for the static Eliashberg solver.

The real-CuPy parity tests skip on machines without a usable CUDA device; the
routing, fallback, and orchestration tests run everywhere (they use a numpy-
delegating spy backend or a monkeypatched missing-CuPy import). Tests must run
from the repository root.
"""
import logging
import types

import numpy as np
import pytest

import hwave.sc as sc
from hwave.solver import backend


def _spy_backend(record):
    """A numpy-delegating stand-in for a device backend module. Every op runs on
    numpy (results stay correct), but asarray/asnumpy calls are recorded so a
    test can assert the operator moved its inputs onto the backend and its
    outputs back to the host. `array_module_of` is monkeypatched to return this,
    so backend.to_host/spatial_* route through it too."""
    def asarray(a, *args, **kwargs):
        record.setdefault("asarray_shapes", []).append(np.asarray(a).shape)
        return np.asarray(a, *args, **kwargs)

    def asnumpy(a):
        record.setdefault("asnumpy_calls", 0)
        record["asnumpy_calls"] += 1
        return np.asarray(a)

    fft = types.SimpleNamespace(ifftn=np.fft.ifftn, fftn=np.fft.fftn)
    return types.SimpleNamespace(
        asarray=asarray, asnumpy=asnumpy, einsum=np.einsum,
        moveaxis=np.moveaxis, fft=fft, newaxis=np.newaxis)


def _small_simple_operator_inputs(seed=0):
    """Tiny simple-mode (5-D vertex) inputs for the kernel operator."""
    norb, Nx, Ny, Nz = 1, 2, 2, 1
    rng = np.random.default_rng(seed)
    def cplx(shape):
        return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
    Vs_q = cplx((norb, norb, Nx, Ny, Nz))
    G2 = cplx((norb, norb, norb, norb, Nx, Ny, Nz))
    return Vs_q, G2, norb, Nx, Ny, Nz


def test_kernel_operator_cpu_returns_host_complex128():
    """On the plain CPU path matvec/matmat return real numpy complex128 arrays
    of the documented shapes (guards against silent CuPy leakage)."""
    Vs_q, G2, norb, Nx, Ny, Nz = _small_simple_operator_inputs()
    A, n = sc._make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    out = A.matvec(np.ones(n, dtype=complex))
    assert isinstance(out, np.ndarray) and out.shape == (n,)
    assert out.dtype == np.complex128
    B = np.ones((n, 3), dtype=complex)
    outm = A.matmat(B)
    assert isinstance(outm, np.ndarray) and outm.shape == (n, 3)
    assert outm.dtype == np.complex128


def test_kernel_operator_routes_input_and_output_through_backend(monkeypatch):
    """matvec/matmat must (a) move their host input onto the invariants'
    backend (xp.asarray) BEFORE any kernel op, and (b) return host arrays via
    backend.to_host. Catches the review must-fix: a missing host->device input
    transfer would break real CuPy but not a spy that returns numpy."""
    record = {}
    spy = _spy_backend(record)
    monkeypatch.setattr(backend, "array_module_of", lambda a: spy)

    Vs_q, G2, norb, Nx, Ny, Nz = _small_simple_operator_inputs()
    # reference (unpatched) result for correctness
    A_ref, n = sc._make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    ref = A_ref.matvec(np.ones(n, dtype=complex))

    A, n = sc._make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    v = np.ones(n, dtype=complex)
    out = A.matvec(v)

    # (a) the input vector (shape (n,)) was passed through xp.asarray
    assert (n,) in record["asarray_shapes"]
    # (b) an output was host-ized
    assert record.get("asnumpy_calls", 0) >= 1
    # correctness preserved
    assert isinstance(out, np.ndarray)
    np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-12)


def test_kernel_operator_asserts_matched_backends(monkeypatch):
    """Vs_q and G2 must live on the same backend."""
    Vs_q, G2, norb, Nx, Ny, Nz = _small_simple_operator_inputs()
    calls = {"n": 0}
    real = backend.array_module_of

    def alternating(a):
        calls["n"] += 1
        return np if calls["n"] == 1 else types.SimpleNamespace()

    monkeypatch.setattr(backend, "array_module_of", alternating)
    with pytest.raises(ValueError, match="same backend"):
        sc._make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)


# --- calc_eliashberg static-path GPU orchestration ------------------------

from tests.test_eliashberg_ir import flex_outdir  # noqa: E402  (module fixture)


def _static_flex_input(outdir):
    """A small static (frequency omitted -> 'static') Eliashberg input that
    reads the FLEX fixture's susceptibilities. Nmat matches the fixture."""
    from tests.test_eliashberg_ir import BETA_T, LX, LY, NMAT
    return {
        "mode": {"param": {"T": BETA_T, "CellShape": [LX, LY, 1],
                           "SubShape": [1, 1, 1], "Nmat": NMAT,
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": "tests/rpa/input",
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": outdir}},
        "eliashberg": {"chi0q_mode": "flex", "solver_mode": "iteration",
                       "max_iter": 50, "convergence_tol": 1e-7},
    }


def _no_cupy(*args, **kwargs):
    raise ImportError("No module named 'cupy'")


def test_static_gpu_true_no_longer_raises_and_falls_back(flex_outdir, tmp_path,
                                                         monkeypatch, caplog):
    """gpu=true on the static path no longer raises; without CuPy it warns and
    completes on numpy."""
    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    inp = _static_flex_input(flex_outdir)
    inp["eliashberg"]["gpu"] = True
    with caplog.at_level(logging.WARNING):
        sc.calc_eliashberg(inp)   # must not raise
    assert any("cupy" in r.getMessage().lower() for r in caplog.records)
    import os
    assert os.path.exists(os.path.join(inp["file"]["output"]["path_to_output"],
                                       "eigenvalue.dat"))


def test_static_gpu_required_fails_fast_without_cupy(flex_outdir, tmp_path,
                                                     monkeypatch):
    """gpu_required=true must raise (not silently fall back) when CuPy is
    missing, on the static path too (issue #63 contract)."""
    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    inp = _static_flex_input(flex_outdir)
    inp["eliashberg"]["gpu"] = True
    inp["eliashberg"]["gpu_required"] = True
    with pytest.raises(RuntimeError, match="gpu_required"):
        sc.calc_eliashberg(inp)


def test_static_gpu_false_default_unchanged(flex_outdir, tmp_path):
    """The default (gpu=false) static solve still runs and writes results."""
    import os
    inp = _static_flex_input(flex_outdir)
    sc.calc_eliashberg(inp)
    assert os.path.exists(os.path.join(inp["file"]["output"]["path_to_output"],
                                       "eigenvalue.dat"))
