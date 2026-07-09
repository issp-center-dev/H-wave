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
                "falling back to the numpy (CPU) backend.")
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


def to_host(arr):
    """Return ``arr`` as a numpy array (device-to-host copy for cupy input,
    identity for numpy input)."""
    xp = array_module_of(arr)
    if xp is np:
        return arr
    return xp.asnumpy(arr)
