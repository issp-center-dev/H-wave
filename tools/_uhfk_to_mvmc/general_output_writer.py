"""InOrbitalGeneral aggregator and writer with class-consistency checks.

See docs/en/source/algorithm/uhfk_to_mvmc.rst for F construction and
docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for the consistency check.
"""
from __future__ import annotations

import numpy as np


class ClassInconsistencyError(RuntimeError):
    """Raised when F entries mapped to the same orbitalidx class differ
    beyond class_consistency_tol; indicates a mismatch between the
    input Slater's invariances and the orbitalidx_general.def class
    structure."""


def aggregate_general_orbital_params(
    F: np.ndarray,
    mapping: dict,
    n_orbital_idx: int,
    *,
    epsilon_noise: float = 1.0e-8,
    complex_type: int = 1,
    rng: np.random.Generator = None,
    class_consistency_tol: float = 1.0e-8,
) -> np.ndarray:
    """Aggregate F[all_i, all_j] * sign into per-class averages, then add
    small uniform rank-lift noise.

    A fail-fast consistency check runs BEFORE averaging: for every idx,
    the maximum residual |value_i - value_j| across all signed F entries
    assigned to that idx must be <= class_consistency_tol. Otherwise
    raise ClassInconsistencyError with the offending idx and observed
    max residual.

    Parameters
    ----------
    F : (2Ns, 2Ns) complex antisymmetric ndarray
    mapping : dict[(all_i, all_j)] = (idx, sign)
        From orbitalidx_general_reader.parse_orbitalidx_general_def.
    n_orbital_idx : int
        Total parameter count (from orbitalidx_general.def header).
    epsilon_noise : float, optional
        Rank-lift noise amplitude (default 1e-8).
    complex_type : int, optional
        1 → noise both real and imag parts; 0 → real only.
    rng : numpy.random.Generator, optional
        Defaults to np.random.default_rng(7919) for reproducibility.
    class_consistency_tol : float, optional
        Max allowed residual per class before averaging (default 1e-8).

    Returns
    -------
    params : (n_orbital_idx,) complex128
    """
    F = np.asarray(F, dtype=np.complex128)
    values_per_idx = {}
    for (all_i, all_j), (idx, sign) in mapping.items():
        if idx < 0:
            continue
        values_per_idx.setdefault(idx, []).append(sign * F[all_i, all_j])

    # Class-consistency check (fail-fast BEFORE averaging).
    for idx, vals in values_per_idx.items():
        if len(vals) <= 1:
            continue
        arr = np.asarray(vals, dtype=np.complex128)
        residual = float(np.max(np.abs(arr - arr.mean())))
        if residual > class_consistency_tol:
            raise ClassInconsistencyError(
                f"orbitalidx class idx={idx} has inconsistent signed F "
                f"entries: max residual {residual:.3e} exceeds "
                f"class_consistency_tol {class_consistency_tol:.3e}. "
                f"Values (first 5): {arr[:5].tolist()!r}. "
                "This indicates orbitalidx_general.def encodes a symmetry "
                "the Slater state does not respect; regenerate it with "
                "classes matching the state."
            )

    params = np.zeros(n_orbital_idx, dtype=np.complex128)
    for idx, vals in values_per_idx.items():
        arr = np.asarray(vals, dtype=np.complex128)
        params[idx] = arr.mean()

    if epsilon_noise > 0:
        if rng is None:
            rng = np.random.default_rng(7919)
        noise_re = rng.uniform(-0.5, 0.5, size=n_orbital_idx) * epsilon_noise
        if complex_type == 1:
            noise_im = rng.uniform(-0.5, 0.5, size=n_orbital_idx) * epsilon_noise
            params += noise_re + 1j * noise_im
        else:
            params += noise_re
    return params


def write_zqp_orbital_general(out_path: str, params: np.ndarray) -> None:
    """Write mVMC-compatible zqp_orbital_uhfk.dat (same 5-line header +
    ``<idx> <real> <imag>`` body as ``output_writer.write_zqp_orbital``).
    mVMC's namelist.def entry becomes ``InOrbitalGeneral <path>``.
    """
    params = np.asarray(params, dtype=np.complex128)
    n = int(len(params))
    with open(out_path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx  {n}\n")
        fp.write("======================\n")
        fp.write("== i_j_OrbitalIdx ===\n")
        fp.write("======================\n")
        for i in range(n):
            fp.write(
                "{:d} {: .18e} {: .18e}\n".format(
                    i, float(params[i].real), float(params[i].imag)
                )
            )
