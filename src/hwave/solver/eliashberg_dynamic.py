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
