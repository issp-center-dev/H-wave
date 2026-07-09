"""Dynamic Eliashberg equation solver (frequency-resolved).

This module implements the full-frequency Eliashberg equation solver for
analyzing superconducting instabilities with frequency-dependent kernels.
"""

import os
import logging
import numpy as np

logger = logging.getLogger("qlms").getChild("eliashberg_dynamic")


def estimate_memory_bytes(norb, Nk, nmat):
    v = 16 * norb**4 * Nk * nmat          # vertex
    g2 = v                                # G2, same shape
    vec = 16 * norb**2 * Nk * nmat        # one gap/Arnoldi vector
    return int((v + g2 + 4 * vec) * 1.5)  # + ~4 work vectors, 1.5x FFT headroom


def _available_ram_bytes():
    try:
        import psutil
        return int(psutil.virtual_memory().available)
    except Exception:
        try:
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) * 1024
        except Exception:
            pass
    return 8 * 2**30  # fallback: 8 GiB


def check_memory(norb, Nk, nmat, mem_limit_gb=None):
    need = estimate_memory_bytes(norb, Nk, nmat)
    if mem_limit_gb == 0:
        logger.info("dynamic Eliashberg est. peak %.2f GiB (guard disabled)",
                    need / 2**30)
        return
    limit = ((0.8 * _available_ram_bytes()) if mem_limit_gb is None
             else mem_limit_gb * 2**30)
    logger.info("dynamic Eliashberg est. peak %.2f GiB (limit %.2f GiB)",
                need / 2**30, limit / 2**30)
    if need > limit:
        raise MemoryError(
            "dynamic Eliashberg estimated peak {:.1f} GiB exceeds the limit "
            "{:.1f} GiB; reduce Nmat/Nk or set [eliashberg] mem_limit_gb "
            "(0 disables).".format(need / 2**30, limit / 2**30))


# NOTE: do NOT write a private chi loader here — path/name/convention/
# spin-orbital-expansion logic lives ONLY in sc._load_flex_susceptibilities_full
# (Task 1c). load_flex_chi_dynamic delegates to it to avoid contract drift.
# The module also needs the FFT alias for the Task-6 kernel:
#   from hwave.solver.perf import FFT

def load_flex_chi_dynamic(input_dict, norb, Nx, Ny, Nz):
    import hwave.sc as sc
    cfg_nmat = int(input_dict["mode"]["param"].get("Nmat", 1024))
    if cfg_nmat % 2 != 0:
        raise ValueError("dynamic Eliashberg requires even Nmat; got {}".format(cfg_nmat))
    # reuse the exact static path/name/convention/spin-orbital-expansion logic
    chis_w, chic_w, green_w, chi_convention = \
        sc._load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz)
    if not (chis_w.shape[-1] == chic_w.shape[-1] == green_w.shape[-1] == cfg_nmat):
        raise ValueError(
            "dynamic Eliashberg grid mismatch: nmat differs — chis={}, chic={}, "
            "green={}, config Nmat={}".format(
                chis_w.shape[-1], chic_w.shape[-1], green_w.shape[-1], cfg_nmat))
    check_memory(norb, Nx*Ny*Nz, cfg_nmat,
                 input_dict["eliashberg"].get("mem_limit_gb"))
    return chis_w, chic_w, green_w, chi_convention


def compute_vertices_flex_dynamic(chis_w, chic_w, inter_k, norb,
                                  Nx, Ny, Nz, pairing_type, convention):
    """Full-frequency pairing vertex: apply sc._compute_vertices_flex per
    bosonic Matsubara frequency and stack along the trailing axis.

    Parameters
    ----------
    chis_w, chic_w : ndarray
        Spin/charge susceptibilities, shape (Nx, Ny, Nz, nd, nd, nmat)
        with nd = norb**2.
    inter_k : dict
        Interactions in k-space from sc._build_interaction_k.
    norb, Nx, Ny, Nz : int
        Orbital count and grid dimensions.
    pairing_type : str
        "singlet" or "triplet".
    convention : str
        "kuroki" or "myo" — MUST match the orbital convention chis_w/chic_w
        were produced in (chi_convention from the FLEX loader); forwarded
        unchanged to sc._compute_vertices_flex so the matching S/C matrices
        are used at every frequency.

    Returns
    -------
    Vs_q_w : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat).
    """
    import hwave.sc as sc
    nmat = chis_w.shape[-1]

    def _one(l):
        return sc._compute_vertices_flex(
            chis_w[..., l], chic_w[..., l], inter_k, norb, Nx, Ny, Nz,
            pairing_type=pairing_type, convention=convention)

    v0 = _one(0)
    out = np.empty(v0.shape + (nmat,), dtype=v0.dtype)
    out[..., 0] = v0
    for l in range(1, nmat):
        out[..., l] = _one(l)
    return out


def calc_g2_dynamic(green_kw, beta):
    """Frequency-resolved pair bubble: identical to sc._calc_g2 except the
    fermionic Matsubara sum is NOT taken (the frequency axis is kept).

    sc._calc_g2 computes G2[i,j,l,m,x,y,z] = (1/beta) * sum_n
    green_kw[i,j,x,y,z,n] * green_kw_inv[l,m,x,y,z,n], where green_kw_inv is
    G(-k,-wn) built via roll+flip. This function drops the sum over n and
    returns the per-frequency summand, so calc_g2_dynamic(...).sum(axis=-1)
    reproduces sc._calc_g2(...) to machine precision (see
    tests/test_eliashberg_dynamic.py::test_g2_dynamic_sums_to_static).

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    beta : float
        Inverse temperature.

    Returns
    -------
    G2_w : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat).
    """
    norb = green_kw.shape[0]
    Nx, Ny, Nz, nmat = green_kw.shape[2], green_kw.shape[3], green_kw.shape[4], green_kw.shape[5]
    nvol = Nx * Ny * Nz

    # G(-k, -wn) via roll+flip -- SAME construction as sc._calc_g2.
    green_kw_inv = np.roll(
        green_kw[:, :, ::-1, ::-1, ::-1, ::-1],
        (1, 1, 1), (2, 3, 4)
    )
    # Same reshape/index layout as sc._calc_g2's A/B (ij, site, n) and
    # (lm, site, n), but keep the per-frequency product instead of summing
    # (matmul-)contracting over n.
    A = green_kw.reshape(norb * norb, nvol, nmat)      # (ij, site, n)
    B = green_kw_inv.reshape(norb * norb, nvol, nmat)  # (lm, site, n)
    # G2[ij, lm, site, n] = A[ij, site, n] * B[lm, site, n]  (no sum over n)
    G2 = A[:, np.newaxis, :, :] * B[np.newaxis, :, :, :]
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    return G2 / beta


def solve_dynamic(input_dict):
    """Solve the dynamic (frequency-resolved) Eliashberg equation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Raises
    ------
    NotImplementedError
        This is a stub; implementation is in later tasks.
    """
    raise NotImplementedError("dynamic Eliashberg solver: implemented in later tasks")
