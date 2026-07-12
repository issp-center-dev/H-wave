"""Array-backend selection (numpy / cupy) for GPU-accelerated solvers.

The dynamic Eliashberg solver (and, prospectively, FLEX) works on large dense
complex arrays through FFTs, einsums and elementwise products -- operations
with 1:1 CuPy equivalents. This module centralizes the choice of array module
so solver code can stay backend-agnostic: pick the module once with
``get_backend`` and dispatch per-array with ``array_module_of``.

GPU use is strictly opt-in (``gpu = true`` in the ``[eliashberg]`` section);
when CuPy or a CUDA device is unavailable the run falls back to numpy with a
warning -- the result is identical, only slower.
"""
import logging

import numpy as np

try:
    import scipy.fft as _SFFT
except ImportError:                       # pragma: no cover - scipy is a dep
    _SFFT = None


def as_bool(value):
    """Coerce a config value to bool. TOML booleans arrive as real bools;
    programmatic inputs sometimes carry strings, where plain truthiness
    would read "false"/"off"/"0" as True."""
    if isinstance(value, str):
        return value.strip().lower() in ("true", "yes", "on", "1")
    return bool(value)


def _import_cupy():
    """Import and return cupy (separated out so tests can monkeypatch a
    missing installation)."""
    import cupy
    return cupy


def get_backend(use_gpu, logger=None, required=False):
    """Select the array module.

    Parameters
    ----------
    use_gpu : bool
        Whether GPU execution was requested.
    logger : logging.Logger, optional
        Logger for the fallback warning.
    required : bool, optional
        Strict mode. When True and GPU execution was requested but no usable
        CuPy/CUDA backend exists, raise ``RuntimeError`` instead of silently
        falling back to the (much slower) numpy path -- so a large scheduler
        job fails fast rather than turning a short GPU run into a very long CPU
        run. Opt-in via ``gpu_required`` so the default stays compatible.

    Returns
    -------
    xp : module
        ``cupy`` when GPU execution is requested and usable, else ``numpy``.
    gpu_active : bool
        True only when ``xp`` is cupy with at least one CUDA device.

    Raises
    ------
    RuntimeError
        When ``required`` and ``use_gpu`` but the GPU backend is unusable.
    """
    if not use_gpu:
        return np, False
    try:
        cupy = _import_cupy()
    except ImportError:
        msg = ("CuPy is not installed. Install the precompiled wheel matching "
               "your CUDA version (e.g. 'pip install cupy-cuda12x' for CUDA "
               "12.x); see https://docs.cupy.dev/en/stable/install.html")
        if required:
            raise RuntimeError(
                "gpu_required=true but " + msg[0].lower() + msg[1:])
        if logger is not None:
            logger.warning(
                "gpu=true requested but %s falling back to the numpy (CPU) "
                "backend.", msg)
        return np, False
    try:
        ndev = cupy.cuda.runtime.getDeviceCount()
        if ndev < 1:
            raise RuntimeError("no CUDA device")
    except Exception as exc:
        if required:
            raise RuntimeError(
                "gpu_required=true but no usable CUDA device was found "
                "({}). Check the driver/runtime, or unset gpu_required to "
                "allow the CPU fallback.".format(exc))
        if logger is not None:
            logger.warning(
                "gpu=true requested but no usable CUDA device was found "
                "(%s); falling back to the numpy (CPU) backend.", exc)
        return np, False
    return cupy, True


def array_module_of(arr):
    """Return the array module (numpy or cupy) that owns ``arr``."""
    if type(arr).__module__.split(".")[0] == "cupy":
        return _import_cupy()
    return np


def warn_if_device_memory_short(required_bytes, logger, label=""):
    """Warn when ``required_bytes`` exceeds the free memory of the current
    CUDA device. Advisory only (CuPy itself raises a clear OutOfMemoryError
    on allocation failure); a pre-warning lets long pipelines flag doomed
    configurations before spending time on the transfer."""
    try:
        cupy = _import_cupy()
        free_b, total_b = cupy.cuda.runtime.memGetInfo()
    except Exception:
        return
    if required_bytes > free_b and logger is not None:
        logger.warning(
            "%s requires an estimated %.2f GB of GPU memory but only "
            "%.2f GB of %.2f GB is free; expect a CuPy OutOfMemoryError. "
            "Reduce Nmat / the k-mesh, or run with gpu=false.",
            label or "the GPU transfer", required_bytes / 1e9,
            free_b / 1e9, total_b / 1e9)


def restore_host_attrs(obj, names):
    """Host-restore the named array attributes of ``obj`` in place.

    For each name present and non-None, replace the attribute with
    ``to_host(value)``. Identity for numpy/None attributes, so it is safe to
    call unconditionally (including on the CPU path). Used in the solvers'
    ``finally`` cleanup so public state (H0 eigenpairs, stored Green functions)
    is NumPy-backed after both normal completion and a GPU-path exception.

    Restoration is best-effort and exception-safe: because this runs inside a
    ``finally``, a failing device->host copy (e.g. a dead CUDA context after
    the original error) must NOT mask the original solver exception nor stop
    the remaining attributes from being restored. Such a failure is logged and
    the offending attribute is left as-is.
    """
    for name in names:
        try:
            val = getattr(obj, name, None)
            if val is None:
                continue
            setattr(obj, name, to_host(val))
        except Exception:                       # noqa: BLE001 - never mask caller's error
            logging.getLogger("qlms").warning(
                "restore_host_attrs: could not host-restore %r; leaving it "
                "backend-local.", name, exc_info=True)


def to_host(arr):
    """Return ``arr`` as a numpy array (device-to-host copy for cupy input,
    identity for numpy input)."""
    xp = array_module_of(arr)
    if xp is np:
        return arr
    return xp.asnumpy(arr)


def spatial_ifftn(a, axes, workers=1):
    """Spatial inverse FFT shared by the k-space solvers (RPA/FLEX/dynamic
    Eliashberg). On the numpy backend this is parallelized via scipy.fft when
    ``workers`` asks for it (``!= 1``) and scipy is available (scipy's result
    matches numpy to machine precision); on the cupy backend the FFT already
    runs on the GPU and ``workers`` is ignored."""
    xp = array_module_of(a)
    if xp is np and _SFFT is not None and workers not in (None, 0, 1):
        return _SFFT.ifftn(a, axes=axes, workers=workers)
    return xp.fft.ifftn(a, axes=axes)


def spatial_fftn(a, axes, workers=1):
    """Spatial forward FFT; see :func:`spatial_ifftn`."""
    xp = array_module_of(a)
    if xp is np and _SFFT is not None and workers not in (None, 0, 1):
        return _SFFT.fftn(a, axes=axes, workers=workers)
    return xp.fft.fftn(a, axes=axes)
