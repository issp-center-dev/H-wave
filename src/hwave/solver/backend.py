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
import numpy as np

try:
    import scipy.fft as _SFFT
except ImportError:                       # pragma: no cover - scipy is a dep
    _SFFT = None


def _import_cupy():
    """Import and return cupy (separated out so tests can monkeypatch a
    missing installation)."""
    import cupy
    return cupy


def get_backend(use_gpu, logger=None):
    """Select the array module.

    Parameters
    ----------
    use_gpu : bool
        Whether GPU execution was requested.
    logger : logging.Logger, optional
        Logger for the fallback warning.

    Returns
    -------
    xp : module
        ``cupy`` when GPU execution is requested and usable, else ``numpy``.
    gpu_active : bool
        True only when ``xp`` is cupy with at least one CUDA device.
    """
    if not use_gpu:
        return np, False
    try:
        cupy = _import_cupy()
    except ImportError:
        if logger is not None:
            logger.warning(
                "gpu=true requested but CuPy is not installed; "
                "falling back to the numpy (CPU) backend. Install the "
                "precompiled wheel matching your CUDA version (e.g. "
                "'pip install cupy-cuda12x' for CUDA 12.x); see "
                "https://docs.cupy.dev/en/stable/install.html")
        return np, False
    try:
        ndev = cupy.cuda.runtime.getDeviceCount()
        if ndev < 1:
            raise RuntimeError("no CUDA device")
    except Exception as exc:
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
