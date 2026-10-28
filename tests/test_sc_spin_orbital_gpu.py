"""``[eliashberg] gpu = true`` in spin-orbital mode (issue #274): the static
(#269) and the dynamic (#284) spin-orbital Eliashberg solvers give the same
eigenvalues and gaps on the device backend as on the host.

With ``gpu = true`` both solvers park their two large invariants (the
pairing vertex and the pair bubble G2) on the device; the kernel operator
then moves each host vector (``matvec``) or host block of columns
(``matmat``, which the antisymmetric basis ``sc._SOAntisymmetricBasis``
forwards on the dense full-spectrum branch) to the device and returns host
arrays to the host-side eigensolvers.

* ``TestFakeDeviceStatic`` / ``TestFakeDeviceDynamic`` run everywhere. A
  numpy-backed fake CuPy (:class:`_FakeCupy`) is installed behind
  ``backend._import_cupy``, so ``get_backend`` selects it and
  ``backend.array_module_of`` routes its arrays to it -- the production
  device branch runs unmodified. The fake is STRICT like CuPy: its array
  type is not an ``np.ndarray``, implicit host conversion
  (``np.asarray(device_array)``) raises ``TypeError``, and an operation
  mixing a host ``ndarray`` with a device array raises ``TypeError``. Every
  host->device transfer and device->host copy is recorded, so each test
  asserts that the invariants were parked and that the operator's
  ``matvec`` / ``matmat`` inputs crossed to the device.
* ``TestRealDeviceStatic`` / ``TestRealDeviceDynamic`` run the same
  comparisons on an actual CUDA device and are skipped (``setUpClass``,
  ``unittest.SkipTest``: a clean skip under both pytest and unittest) when
  CuPy or a device is unavailable. ``gpu_required = true`` makes a silent
  CPU fallback an error there.

Tests must run from the repository root.
"""
import contextlib
import logging
import shutil
import tempfile
import types
import unittest
from unittest import mock

import numpy as np
from numpy.lib.mixins import NDArrayOperatorsMixin

import hwave.sc as sc
from hwave.solver import backend
import tests.test_sc_spin_orbital as S

#: Nmat of the dynamic runs (as the end-to-end tests of
#: tests/test_eliashberg_dynamic_so.py)
_NMAT_DYN = 8


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


# ---------------------------------------------------------------------------
# a strict, numpy-backed fake CuPy
# ---------------------------------------------------------------------------

class _FakeDeviceArray(NDArrayOperatorsMixin):
    """A device array stand-in. ``__module__`` starts with ``cupy`` so
    ``backend.array_module_of`` sends it to ``backend._import_cupy`` (the
    fake module). It is NOT an ``np.ndarray``: ``__array__`` raises as
    CuPy's does, and ufuncs / numpy functions refuse host ``ndarray``
    operands (0-d arrays and scalars are accepted, as by CuPy)."""

    def __init__(self, a):
        self._a = a

    def __array__(self, *args, **kwargs):
        raise TypeError("Implicit conversion to a NumPy array is not allowed "
                        "(fake device array)")

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if "out" in kwargs:
            kwargs["out"] = tuple(_dev_unwrap(o) for o in kwargs["out"])
        return _dev_wrap(getattr(ufunc, method)(
            *[_dev_unwrap(x) for x in inputs], **kwargs))

    def __array_function__(self, func, types_, args, kwargs):
        return _dev_wrap(func(*_dev_unwrap_tree(args),
                              **_dev_unwrap_tree(kwargs)))

    def __getattr__(self, name):
        if name.startswith("__"):
            # never expose __array_interface__ / __array_struct__ etc.: numpy
            # would read the host buffer through them
            raise AttributeError(name)
        attr = getattr(self._a, name)
        if callable(attr):
            def method(*args, **kwargs):
                return _dev_wrap(attr(*_dev_unwrap_tree(args),
                                      **_dev_unwrap_tree(kwargs)))
            return method
        return _dev_wrap(attr)

    def __getitem__(self, key):
        return _dev_wrap(self._a[_dev_unwrap_tree(key)])

    def __setitem__(self, key, value):
        self._a[_dev_unwrap_tree(key)] = _dev_unwrap(value)

    def __len__(self):
        return len(self._a)

    def __iter__(self):
        return (_dev_wrap(x) for x in self._a)

    def __complex__(self):
        return complex(self._a)

    def __float__(self):
        return float(self._a)

    def get(self):
        return self._a.copy()


_FakeDeviceArray.__module__ = "cupy._core.core"


def _dev_unwrap(x):
    if isinstance(x, _FakeDeviceArray):
        return x._a
    if isinstance(x, np.ndarray) and x.ndim > 0:
        raise TypeError("host numpy.ndarray operand mixed with a (fake) "
                        "device array, shape {}".format(x.shape))
    return x


def _dev_unwrap_tree(x):
    if isinstance(x, tuple):
        return tuple(_dev_unwrap_tree(y) for y in x)
    if isinstance(x, list):
        return [_dev_unwrap_tree(y) for y in x]
    if isinstance(x, dict):
        return {k: _dev_unwrap_tree(v) for k, v in x.items()}
    return _dev_unwrap(x)


def _dev_wrap(x):
    if isinstance(x, np.ndarray):
        return _FakeDeviceArray(x)
    if isinstance(x, tuple):
        return tuple(_dev_wrap(y) for y in x)
    return x


class _FakeCupy:
    """The module side: ``asarray`` (host->device, recorded), ``asnumpy``
    (device->host, recorded), ``fft``, ``cuda.runtime``, a memory pool, and
    every other ``numpy`` function applied to unwrapped device operands
    (host ``ndarray`` operands refused) with device results. ``ops`` counts
    the device-side calls by name."""

    ndarray = _FakeDeviceArray
    newaxis = None

    def __init__(self):
        self.to_device = []        # shapes of host arrays moved to the device
        self.to_host = 0           # device->host copies
        self.ops = {}
        self.cuda = types.SimpleNamespace(
            runtime=types.SimpleNamespace(
                getDeviceCount=lambda: 1,
                memGetInfo=lambda: (64 * 2 ** 30, 80 * 2 ** 30)),
            memory=types.SimpleNamespace(OutOfMemoryError=MemoryError))
        pool = types.SimpleNamespace(used_bytes=lambda: 0, total_bytes=lambda: 0,
                                     free_all_blocks=lambda: None)
        self.get_default_memory_pool = lambda: pool
        self.fft = types.SimpleNamespace(**{
            n: self._device_fn("fft." + n, getattr(np.fft, n))
            for n in ("fft", "ifft", "fftn", "ifftn")})

    def _device_fn(self, name, fn):
        def call(*args, **kwargs):
            self.ops[name] = self.ops.get(name, 0) + 1
            return _dev_wrap(fn(*_dev_unwrap_tree(args),
                                **_dev_unwrap_tree(kwargs)))
        return call

    def __getattr__(self, name):
        if name.startswith("__"):
            raise AttributeError(name)
        return self._device_fn(name, getattr(np, name))

    def asarray(self, a, dtype=None):
        if isinstance(a, _FakeDeviceArray):
            return a if dtype is None else _FakeDeviceArray(a._a.astype(dtype))
        host = np.array(a, dtype=dtype)          # a copy, as a transfer is
        self.to_device.append(host.shape)
        return _FakeDeviceArray(host)

    def asnumpy(self, a):
        self.to_host += 1
        if isinstance(a, _FakeDeviceArray):
            return a._a.copy()
        return np.asarray(a)


class TestFakeCupyIsStrict(unittest.TestCase):
    """The fake refuses what CuPy refuses, so a missing host->device transfer
    in the solver fails here instead of passing silently."""

    def test_refusals(self):
        xp = _FakeCupy()
        d = xp.asarray(np.ones((2, 2)))
        self.assertEqual(type(d).__module__.split(".")[0], "cupy")
        with mock.patch.object(backend, "_import_cupy", lambda: xp):
            self.assertIs(backend.array_module_of(d), xp)
            np.testing.assert_array_equal(backend.to_host(d * 2), 2 * np.ones((2, 2)))
        with self.assertRaises(TypeError):
            np.asarray(d)
        with self.assertRaises(TypeError):
            d * np.ones((2, 2))
        with self.assertRaises(TypeError):
            xp.einsum("ij,jk->ik", d, np.ones((2, 2)))
        with self.assertRaises(TypeError):
            d[0] = np.ones(2)
        self.assertEqual(xp.to_device, [(2, 2)])


# ---------------------------------------------------------------------------
# comparison helpers
# ---------------------------------------------------------------------------

def _match_multiset(test, want, got, rtol):
    """One-to-one greedy matching of ``got`` onto ``want``; asserts every
    pair agrees to ``rtol`` of the spectral scale."""
    want, got = np.asarray(want), list(np.asarray(got))
    test.assertEqual(len(want), len(got))
    scale = max(1.0, float(np.abs(want).max())) if want.size else 1.0
    worst = 0.0
    for w in sorted(want, key=lambda z: (-z.real, z.imag)):
        j = int(np.argmin([abs(w - g) for g in got]))
        worst = max(worst, abs(w - got[j]))
        got.pop(j)
    test.assertLessEqual(worst, rtol * scale)
    return worst


def _assert_same_eigenspaces(test, vals_a, vecs_a, vals_b, vecs_b,
                             cluster_rtol, tol, truncated):
    """Each eigenvalue cluster (equal to ``cluster_rtol`` of the spectral
    scale) spans the same subspace in both runs: the sine of the largest
    principal angle between the two spans is at most ``tol`` (phase- and
    degeneracy-invariant). ``vecs_*`` are sequences of eigenvectors, any
    shape each. ``truncated``: the sets are a selection, so the cluster at
    the selection boundary (smallest magnitude / real part) may be cut
    differently in the two runs -- a degenerate partner returned by one run
    and not the other -- and is skipped. Returns the number of clusters
    compared."""
    va = np.asarray(vals_a)
    vb = np.asarray(vals_b)
    A = np.stack([np.asarray(v).ravel() for v in vecs_a], axis=1)
    B = np.stack([np.asarray(v).ravel() for v in vecs_b], axis=1)
    atol = cluster_rtol * max(1.0, float(np.abs(va).max()))
    done = np.zeros(len(va), dtype=bool)
    n_checked = 0
    for i in range(len(va)):
        if done[i]:
            continue
        ia = np.abs(va - va[i]) <= atol
        ib = np.abs(vb - va[i]) <= atol
        done |= ia
        if truncated and (
                abs(abs(va[i]) - min(np.abs(va).min(), np.abs(vb).min())) <= atol
                or abs(va[i].real - min(va.real.min(), vb.real.min())) <= atol):
            continue
        with test.subTest(eigenvalue=va[i]):
            test.assertEqual(int(ia.sum()), int(ib.sum()))
            qa, _ = np.linalg.qr(A[:, ia])
            qb, _ = np.linalg.qr(B[:, ib])
            test.assertLessEqual(
                np.linalg.norm(qa - qb @ (qb.conj().T @ qa), 2), tol)
        n_checked += 1
    return n_checked


def _phase_aligned_distance(a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    ov = np.vdot(a, b)
    b = b * np.exp(-1j * np.angle(ov))
    return float(np.linalg.norm(a - b))


# ---------------------------------------------------------------------------
# the cases (shared by the fake and the real device)
# ---------------------------------------------------------------------------

def _dyn_run(inp):
    """``calc_eliashberg`` on the dynamic spin-orbital path with the solver
    results captured: eigenvalue leg ``vals`` / ``vecs`` (reduced
    coordinates, columns), iteration leg ``gap`` / ``lam``, and the basis."""
    cap = {}
    real_eig = sc._so_solve_eigenvalue
    real_it = sc._so_solve_iteration

    def eig(A, basis, *a, **kw):
        out = real_eig(A, basis, *a, **kw)
        cap["basis"] = basis
        cap["vals"], cap["vecs"] = out
        return out

    def it(A, basis, *a, **kw):
        out = real_it(A, basis, *a, **kw)
        cap["basis"] = basis
        cap["gap"], cap["it_lam"] = out[0], out[1]
        return out

    with _quiet(), mock.patch.object(sc, "_so_solve_eigenvalue", side_effect=eig), \
            mock.patch.object(sc, "_so_solve_iteration", side_effect=it):
        cap["lam"] = sc.calc_eliashberg(inp)
    return cap


class _Cases:
    """The CPU-vs-device harness. A subclass supplies ``_device()`` (a
    context manager for the device run, yielding a recorder or None) and the
    tolerances: ``EV_RTOL`` (eigenvalues, relative to the spectral scale),
    ``SPAN_TOL`` (sine of the largest principal angle between the two
    runs' eigenspaces), ``GAP_TOL`` (phase-aligned
    leading gap, unit norm), ``ITER_RTOL`` / ``ITER_GAP_TOL`` (power
    iteration)."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="sc_so_gpu_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def _both(self, write, make_input, run):
        """``(cpu, device, recorder)`` captures of the same input."""
        files = write(self.d)
        cpu = run(make_input(files, dict(gpu=False), "out_cpu"))
        with self._device() as rec:
            dev = run(make_input(files, dict(gpu=True, gpu_required=True),
                                 "out_gpu"))
        return cpu, dev, rec

    def _static(self, so_kwargs, eli):
        def write(d):
            return S._write_model(d, so=True, **so_kwargs)

        def make_input(files, gpu, out):
            return S._input(self.d, files, True, dict(eli, **gpu), out=out)
        return self._both(write, make_input, S._run)

    def _dynamic(self, eli, cell=S.CELL):
        def write(d):
            return S._write_model(d, so=True, soc=0.3)

        def make_input(files, gpu, out):
            e = {"frequency": "dynamic", "chi0q_mode": "calc"}
            e.update(eli)
            e.update(gpu)
            return S._input(self.d, files, True, e, nmat=_NMAT_DYN, out=out,
                            CellShape=list(cell))
        return self._both(write, make_input, _dyn_run)


class _StaticCases(_Cases):
    """The static spin-orbital solver (#269)."""

    def test_static_su2_full_spectrum(self):
        """SU(2) model, all 120 antisymmetric eigenpairs (the dense branch:
        one ``matmat`` of the whole reduced identity through the kernel)."""
        cpu, dev, rec = self._static(
            {}, {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                 "num_eigenvalues": 254})
        self.assertEqual(len(cpu["evs"]), 120)
        _match_multiset(self, cpu["evs"], dev["evs"], self.EV_RTOL)
        n = _assert_same_eigenspaces(self, cpu["evs"], cpu["vecs"], dev["evs"],
                                     dev["vecs"], 1e-8, self.SPAN_TOL,
                                     truncated=False)
        self.assertGreater(n, 20)      # many distinct (degenerate) clusters
        np.testing.assert_array_equal(dev["match"], cpu["match"])
        if rec is not None:
            n_full = S.ND * S.ND * int(np.prod(S.CELL))
            # the matmat block (all 120 reduced columns) crossed to the device
            self.assertIn((n_full, 120), rec.to_device)

    def test_static_complex_soc_arnoldi(self):
        """Complex spin-orbit hopping, ARPACK (``matvec``) leg; the leading
        gap is non-degenerate and compared directly (phase aligned)."""
        cpu, dev, rec = self._static(
            {"soc": 0.3}, {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                           "num_eigenvalues": 8})
        np.testing.assert_allclose(dev["evs"], cpu["evs"], rtol=0,
                                   atol=self.EV_RTOL * np.abs(cpu["evs"]).max())
        self.assertGreater(cpu["evs"][0].real - cpu["evs"][1].real, 0.05)
        self.assertLessEqual(_phase_aligned_distance(cpu["gap"], dev["gap"]),
                             self.GAP_TOL)
        self.assertGreater(_assert_same_eigenspaces(
            self, cpu["evs"], cpu["vecs"], dev["evs"], dev["vecs"], 1e-8,
            self.SPAN_TOL, truncated=True), 3)
        self.assertIsInstance(dev["gap"], np.ndarray)
        if rec is not None:
            self.assertIn((S.ND * S.ND * int(np.prod(S.CELL)),), rec.to_device)

    def test_static_complex_soc_iteration(self):
        """The shifted power iteration on the device operator."""
        cpu, dev, rec = self._static(
            {"soc": 0.3}, {"chi0q_mode": "calc", "solver_mode": "iteration",
                           "convergence_tol": 1e-10, "max_iter": 20000,
                           "spectral_shift": "auto"})
        self.assertLessEqual(abs(dev["ev"] - cpu["ev"]),
                             self.ITER_RTOL * abs(cpu["ev"]))
        self.assertLessEqual(_phase_aligned_distance(cpu["gap"], dev["gap"]),
                             self.ITER_GAP_TOL)


class _DynamicCases(_Cases):
    """The dynamic spin-orbital solver (#284)."""

    def test_dynamic_arnoldi(self):
        cpu, dev, rec = self._dynamic({"solver_mode": "eigenvalue",
                                       "num_eigenvalues": 6})
        _match_multiset(self, cpu["vals"], dev["vals"], self.EV_RTOL)
        self.assertLessEqual(abs(dev["lam"] - cpu["lam"]),
                             self.EV_RTOL * abs(cpu["lam"]))
        self.assertGreater(_assert_same_eigenspaces(
            self, cpu["vals"], cpu["vecs"].T, dev["vals"], dev["vecs"].T, 1e-8,
            self.SPAN_TOL, truncated=True), 2)
        if rec is not None:
            self.assertIn((S.ND, S.ND) + tuple(S.CELL) + (_NMAT_DYN,),
                          rec.to_device)

    def test_dynamic_full_spectrum(self):
        """Every antisymmetric eigenpair (dense branch: ``matmat`` of the
        reduced identity, applied column by column by the dynamic operator),
        on a 2 x 2 cell to keep the dense problem small."""
        cell = (2, 2, 1)
        dim = S.ND * S.ND * int(np.prod(cell)) * _NMAT_DYN // 2
        cpu, dev, rec = self._dynamic({"solver_mode": "eigenvalue",
                                       "num_eigenvalues": dim}, cell=cell)
        self.assertEqual(len(cpu["vals"]), dim)
        _match_multiset(self, cpu["vals"], dev["vals"], self.EV_RTOL)
        self.assertGreater(_assert_same_eigenspaces(
            self, cpu["vals"], cpu["vecs"].T, dev["vals"], dev["vecs"].T, 1e-8,
            self.SPAN_TOL, truncated=False), 50)

    def test_dynamic_iteration(self):
        cpu, dev, rec = self._dynamic({"solver_mode": "iteration",
                                       "spectral_shift": "auto",
                                       "max_iter": 3000,
                                       "convergence_tol": 1e-10})
        self.assertLessEqual(abs(dev["lam"] - cpu["lam"]),
                             self.ITER_RTOL * abs(cpu["lam"]))
        self.assertLessEqual(_phase_aligned_distance(cpu["gap"], dev["gap"]),
                             self.ITER_GAP_TOL)


class _FakeDevice:
    """CPU vs the fake device: the same numpy arithmetic in the same order,
    so the agreement is to round-off."""

    EV_RTOL = 1e-12
    SPAN_TOL = 1e-10
    GAP_TOL = 1e-10
    ITER_RTOL = 1e-12
    ITER_GAP_TOL = 1e-10

    @contextlib.contextmanager
    def _device(self):
        rec = _FakeCupy()
        with mock.patch.object(backend, "_import_cupy", lambda: rec):
            yield rec
        # the device branch ran: the invariants were parked, the kernel ran
        # its FFTs and contractions on the device, and results came back
        self.assertGreaterEqual(len(rec.to_device), 3)
        self.assertGreater(rec.to_host, 0)
        self.assertGreater(sum(v for k, v in rec.ops.items()
                               if k.startswith("fft.")), 0)


class TestFakeDeviceStatic(_FakeDevice, _StaticCases, unittest.TestCase):

    def test_vertex_and_g2_are_parked(self):
        """The static entry moves the (nd, nd, nd, nd, Nx, Ny, Nz) vertex and
        pair bubble to the device once each: no four-leg tensor crosses
        again per kernel application."""
        _, _, rec = self._static(
            {"soc": 0.3}, {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                           "num_eigenvalues": 2})
        four_leg = (S.ND,) * 4 + tuple(S.CELL)
        self.assertEqual(rec.to_device.count(four_leg), 2)


class TestFakeDeviceDynamic(_FakeDevice, _DynamicCases, unittest.TestCase):

    def test_vertex_and_g2_are_parked(self):
        """The dynamic entry moves the (nd, nd, nd, nd, Nx, Ny, Nz, Nmat)
        vertex and pair bubble to the device once each (the vertex before its
        (r, tau) transform): no four-leg tensor crosses per matvec."""
        _, _, rec = self._dynamic({"solver_mode": "eigenvalue",
                                   "num_eigenvalues": 2})
        four_leg = (S.ND,) * 4 + tuple(S.CELL) + (_NMAT_DYN,)
        self.assertEqual(rec.to_device.count(four_leg), 2)


# ---------------------------------------------------------------------------
# the real device (skipped without CuPy / CUDA)
# ---------------------------------------------------------------------------

class _RealDevice:
    """CPU vs an actual CUDA device: cuFFT / cuBLAS against pocketfft /
    host BLAS, so the agreement is to round-off amplified by the eigensolver;
    the power iteration agrees to its convergence tolerance."""

    EV_RTOL = 1e-10
    SPAN_TOL = 1e-7
    GAP_TOL = 1e-8
    ITER_RTOL = 1e-8
    ITER_GAP_TOL = 1e-4

    @classmethod
    def setUpClass(cls):
        if not backend.gpu_available():
            raise unittest.SkipTest("no usable CUDA device / cupy")

    @contextlib.contextmanager
    def _device(self):
        seen = []
        real = backend.get_backend

        def spy(*a, **kw):
            out = real(*a, **kw)
            seen.append(out[1])
            return out

        with mock.patch.object(backend, "get_backend", side_effect=spy):
            yield None
        self.assertTrue(seen and all(seen))     # the device backend was active


class TestRealDeviceStatic(_RealDevice, _StaticCases, unittest.TestCase):
    pass


class TestRealDeviceDynamic(_RealDevice, _DynamicCases, unittest.TestCase):
    pass


if __name__ == "__main__":
    unittest.main()
