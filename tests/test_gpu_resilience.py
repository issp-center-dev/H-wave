"""GPU execution resilience tests (issue #63): strict fallback, VRAM preflight,
and public-state restoration after a GPU-path exception.

These run without a real CUDA device: a tiny fake CuPy backend routes
``array_module_of`` through a wrapper type whose module is reported as
``cupy``, so ``to_host`` exercises the host-restore path. Tests must run from
the repository root (they use the shared fixtures under tests/rpa/input).
"""
import logging
import os
import types

import numpy as np
import pytest


# --- fake device backend ----------------------------------------------------

class _DevArr:
    """Stand-in for a cupy ndarray. ``__module__`` is forced to "cupy" so
    ``backend.array_module_of`` routes it through ``_import_cupy`` (our fake),
    and it is deliberately NOT an ``np.ndarray`` so a still-on-device attribute
    is distinguishable from a host-restored one."""
    __module__ = "cupy"

    def __init__(self, a):
        self.a = np.asarray(a)

    @property
    def shape(self):
        return self.a.shape


def _fake_cupy(free_bytes=None, total_bytes=None):
    def _mem_get_info():
        if free_bytes is None:
            raise RuntimeError("no meminfo in fake backend")
        return (free_bytes, total_bytes)

    return types.SimpleNamespace(
        asarray=lambda a: a if isinstance(a, _DevArr) else _DevArr(a),
        asnumpy=lambda d: d.a if isinstance(d, _DevArr) else np.asarray(d),
        cuda=types.SimpleNamespace(
            runtime=types.SimpleNamespace(
                getDeviceCount=lambda: 1,
                memGetInfo=_mem_get_info)),
    )


def _install_fake_backend(monkeypatch, fake):
    """Route the solver's backend module at the fake device so gpu_active=True
    and to_host uses the fake asnumpy."""
    from hwave.solver import backend
    monkeypatch.setattr(backend, "_import_cupy", lambda: fake)
    monkeypatch.setattr(backend, "get_backend",
                        lambda use_gpu, logger=None, required=False: (fake, True))


# --- unsolved solver factories ----------------------------------------------

_INPUT = {
    'path_to_input': 'tests/rpa/input',
    'interaction': {
        'path_to_input': 'tests/rpa/input',
        'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat',
        'CoulombIntra': 'coulombintra.dat',
    },
}


def _make_solver(kind, extra_param=None):
    import hwave.qlmsio.read_input_k as read_input_k
    param = {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
             'SubShape': [1, 1, 1], 'Nmat': 16, 'gpu': True}
    if extra_param:
        param.update(extra_param)
    info_mode = {'mode': kind, 'param': param, 'calc_scheme': 'reduced'}
    if kind == 'RPA':
        info_mode['calc_type'] = 'ring'
    read_io = read_input_k.QLMSkInput(_INPUT)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    if kind == 'FLEX':
        import hwave.solver.flex as m
        solver = m.FLEX(ham_info, {}, info_mode)
    else:
        import hwave.solver.rpa as m
        solver = m.RPA(ham_info, {}, info_mode)
    os.makedirs('tests/rpa/output', exist_ok=True)
    return solver, green_info


# --- strict mode (opt-in fail-fast) -----------------------------------------

@pytest.mark.parametrize("kind", ["RPA", "FLEX"])
def test_gpu_required_fails_fast_without_cupy(kind, monkeypatch):
    """gpu_required=true must raise (not silently fall back to CPU) when CuPy
    is unavailable, so a large scheduler job fails fast."""
    from hwave.solver import backend
    monkeypatch.setattr(
        backend, "_import_cupy",
        lambda: (_ for _ in ()).throw(ImportError("no cupy")))
    solver, green_info = _make_solver(kind, {'gpu_required': True})
    assert solver.gpu_required is True
    with pytest.raises(RuntimeError, match="gpu_required"):
        solver.solve(green_info, 'tests/rpa/output')


@pytest.mark.parametrize("kind", ["RPA", "FLEX"])
def test_gpu_required_string_false_disables_strict(kind, monkeypatch):
    """A programmatic string "false" must NOT enable strict mode (as_bool
    coercion): the run falls back to CPU and completes."""
    from hwave.solver import backend
    monkeypatch.setattr(
        backend, "_import_cupy",
        lambda: (_ for _ in ()).throw(ImportError("no cupy")))
    solver, green_info = _make_solver(kind, {'gpu_required': "false"})
    assert solver.gpu_required is False
    solver.solve(green_info, 'tests/rpa/output')   # must not raise


# --- public-state restoration after a GPU-path exception --------------------

@pytest.mark.parametrize("kind", ["RPA", "FLEX"])
def test_public_state_restored_after_gpu_exception(kind, monkeypatch):
    """If a GPU-path operation raises mid-solve, the solver's public H0
    eigenpairs must be host (numpy) arrays afterward, not device arrays."""
    _install_fake_backend(monkeypatch, _fake_cupy())
    solver, green_info = _make_solver(kind)

    def _boom(*a, **k):
        raise RuntimeError("injected GPU-path failure")

    # _calc_green is the first heavy op after H0 is moved to the device.
    monkeypatch.setattr(solver, "_calc_green", _boom)

    with pytest.raises(RuntimeError, match="injected"):
        solver.solve(green_info, 'tests/rpa/output')

    # H0 was converted to the fake device type before the failure; the finally
    # cleanup must have restored it to a host numpy array.
    assert isinstance(solver.H0_eigenvalue, np.ndarray)
    assert isinstance(solver.H0_eigenvector, np.ndarray)
    assert not isinstance(solver.H0_eigenvalue, _DevArr)


# --- VRAM preflight (advisory) ----------------------------------------------

@pytest.mark.parametrize("kind", ["RPA", "FLEX"])
def test_vram_preflight_warns_when_short(kind, monkeypatch, caplog):
    """FLEX/RPA must run a VRAM preflight that warns (advisory) when the
    estimated resident tensors exceed free device memory, naming the solver."""
    _install_fake_backend(monkeypatch, _fake_cupy(free_bytes=1, total_bytes=2))
    solver, green_info = _make_solver(kind)

    # Stop right after the preflight so the test does not run the full solve.
    def _boom(*a, **k):
        raise RuntimeError("stop after preflight")

    monkeypatch.setattr(solver, "_calc_green", _boom)

    with caplog.at_level(logging.WARNING, logger="qlms"):
        with pytest.raises(RuntimeError):
            solver.solve(green_info, 'tests/rpa/output')

    assert any("GPU memory" in rec.message or "free" in rec.message.lower()
               for rec in caplog.records), \
        "no VRAM preflight warning emitted"
