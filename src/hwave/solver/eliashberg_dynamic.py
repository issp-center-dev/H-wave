"""Dynamic Eliashberg equation solver (frequency-resolved).

This module implements the full-frequency Eliashberg equation solver for
analyzing superconducting instabilities with frequency-dependent kernels.
"""

import os
import logging
import numpy as np

from hwave.solver import matsubara as ms
from hwave.solver.perf import FFT

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
    # green_w may be None (missing green.npz); the caller (solve_dynamic) is the
    # fail-fast site for that with a user-facing message. Only cross-check the
    # green frequency axis when it was actually loaded.
    green_nmat = green_w.shape[-1] if green_w is not None else cfg_nmat
    if not (chis_w.shape[-1] == chic_w.shape[-1] == green_nmat == cfg_nmat):
        raise ValueError(
            "dynamic Eliashberg grid mismatch: nmat differs — chis={}, chic={}, "
            "green={}, config Nmat={}".format(
                chis_w.shape[-1], chic_w.shape[-1], green_nmat, cfg_nmat))
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


def eliashberg_kernel_dynamic(Vs_q_w, G2_w, phi_w, norb, beta):
    r"""Apply the frequency-resolved (tau-product) Eliashberg kernel.

    Implements one action of the linearized Eliashberg operator on a trial
    gap ``phi``, keeping the full fermionic Matsubara axis. The structure
    mirrors the static kernel (``sc._make_kernel_operator``) and the FLEX
    self-energy (``flex._calc_self_energy``): an orbital contraction with the
    pair bubble G2, then an imaginary-time PRODUCT with the pairing vertex
    (NOT a circular frequency convolution) done via matsubara transforms.

    Orbital contraction (matches the static kernel):
        F_{l2,l3}(k,m)     = sum_{l5,l6} G2_{l2,l5,l3,l6}(k,m) phi_{l5,l6}(k,m)
        phi_out_{l1,l4}    = - sum_{l2,l3} V_{l1,l2,l3,l4}(r,tau) F_{l2,l3}(r,tau)

    The temperature factor 1/beta lives inside G2 (``calc_g2_dynamic``); the
    kernel carries only the spatial -(1/N) via the ifftn (numpy) convention.
    ``beta`` is accepted for signature symmetry with the static path but is
    not applied here.

    Parameters
    ----------
    Vs_q_w : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        on the bosonic Matsubara axis.
    G2_w : ndarray
        Pair bubble, shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        on the fermionic Matsubara axis (already divided by beta).
    phi_w : ndarray
        Trial gap, shape (norb, norb, Nx, Ny, Nz, nmat).
    norb : int
        Number of orbitals (for signature symmetry; inferred from shapes).
    beta : float
        Inverse temperature (unused here; T is inside G2).

    Returns
    -------
    phi_out_w : ndarray
        Same shape as ``phi_w``.
    """
    # F_{l2,l3}(k, m) = sum_{l5,l6} G2_{l2,l5,l3,l6}(k,m) phi_{l5,l6}(k,m)
    F = np.einsum('iljmxyzn,lmxyzn->ijxyzn', G2_w, phi_w)
    # spatial k->r on F (per orbital pair, per fermionic freq); freq fermion->tau
    F_rt = FFT.ifftn(ms.fermion_to_tau(F, axis=-1), axes=(2, 3, 4))
    # V(q, iv_l) -> (r, tau): spatial q->r, freq boson->tau
    V_rt = FFT.ifftn(ms.boson_to_tau(Vs_q_w, axis=-1), axes=(4, 5, 6))
    # phi_out_{l1,l4}(r,tau) = - sum_{l2,l3} V_{l1,l2,l3,l4}(r,tau) F_{l2,l3}(r,tau)
    prod = -np.einsum('abcdxyzt,bcxyzt->adxyzt', V_rt, F_rt)
    # back: spatial r->k (fftn), freq tau->fermion. The single spatial fold's
    # -(1/N) is already carried by the ifftn above (numpy divides by N on the
    # inverse transform; the r->k fftn does not re-multiply the pair bubble).
    phi_out = ms.tau_to_fermion(FFT.fftn(prod, axes=(2, 3, 4)), axis=-1)
    return phi_out


def frequency_inner(a, b):
    """Inner product over ALL components (orbital, k, and frequency axes)."""
    return np.vdot(a, b)


def _fix_gauge(phi_w):
    r"""Fix the eigenvector gauge deterministically.

    An eigenvector of the (complex, non-Hermitian) Eliashberg kernel is only
    defined up to an overall complex scale. To make the reported gap
    reproducible across runs / linear-algebra backends we pin two freedoms:

    1. **Magnitude** — L2-normalize over *all* components (orbital, k, and the
       full Matsubara axis) so ``||phi|| = 1``.
    2. **Global phase** — multiply by the single phase that makes the
       largest-``|magnitude|`` component real and positive. Ties in the
       magnitude are broken by the first component in lexicographic
       ``(orb1, orb2, kx, ky, kz, iomega)`` order, which is exactly C-order of
       the ``(norb, norb, Nx, Ny, Nz, nmat)`` array (so ``np.argmax`` on the
       raveled magnitudes -- which returns the first maximizing index -- gives
       the tie-break for free).

    The magnitudes are invariant under a global phase, so applying an arbitrary
    phase before ``_fix_gauge`` reproduces the identical array (see
    ``tests/test_eliashberg_dynamic.py::test_gauge_deterministic``).

    Parameters
    ----------
    phi_w : ndarray
        Gap / eigenvector, shape (norb, norb, Nx, Ny, Nz, nmat).

    Returns
    -------
    ndarray
        Same shape as ``phi_w``, L2-normalized and phase-fixed.
    """
    phi = np.asarray(phi_w).astype(complex, copy=True)
    nrm = np.linalg.norm(phi)
    if nrm > 0:
        phi /= nrm
    flat = phi.ravel()  # C-order == lexicographic (orb1,orb2,kx,ky,kz,iomega)
    pivot = flat[int(np.argmax(np.abs(flat)))]
    if pivot != 0:
        phi /= (pivot / abs(pivot))  # rotate so the pivot becomes real-positive
    return phi


def write_dynamic_outputs(output_dir, gap_w, eigenvalue, T, pairing_type,
                          kx_array, ky_array, kz_array, beta,
                          gap_file="gap.dat", npz_file="gap_dynamic.npz"):
    r"""Write the dynamic-Eliashberg gap outputs.

    Produces two files under ``output_dir``:

    * ``gap_dynamic.npz`` -- the full frequency-resolved gap plus metadata.
      Keys: ``gap`` (norb, norb, Nx, Ny, Nz, nmat), ``iomega`` (the centered
      fermionic Matsubara frequencies :math:`\omega_n=(2n+1-N_{mat})\pi T`),
      ``T``, ``pairing_type``, ``frequency`` (== "dynamic"), ``eigenvalue``,
      ``axis_order``, and ``normalization`` (documenting the ``_fix_gauge``
      convention).
    * ``gap.dat`` -- a plain-text slice at a single Matsubara frequency (the
      smallest positive :math:`\omega_n`, index ``nmat//2``), mirroring the
      static ``sc._save_results`` column layout (kx ky kz then Re/Im per orbital
      pair). Its FIRST line is a ``#``-prefixed header that carries
      ``frequency=dynamic`` together with the slice index and its
      :math:`\omega_n` value.

    Parameters
    ----------
    output_dir : str
        Destination directory (created if absent).
    gap_w : ndarray
        Gauge-fixed gap, shape (norb, norb, Nx, Ny, Nz, nmat).
    eigenvalue : float
        Leading eigenvalue lambda.
    T : float
        Temperature.
    pairing_type : str
        "singlet" or "triplet".
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    beta : float
        Inverse temperature (accepted for signature symmetry; iomega uses T).
    gap_file, npz_file : str
        Output filenames.
    """
    os.makedirs(output_dir, exist_ok=True)
    norb = gap_w.shape[0]
    Nx, Ny, Nz, nmat = gap_w.shape[2:6]

    # centered fermionic Matsubara frequencies: iw_n = (2n + 1 - Nmat) pi T
    n_idx = np.arange(nmat)
    iomega = (2.0 * n_idx + 1.0 - nmat) * np.pi * T

    axis_order = "(orb1, orb2, kx, ky, kz, iomega)"
    normalization = ("L2-normalized over all (orb,k,iomega) components; global "
                     "phase fixes the largest-|magnitude| component real-"
                     "positive (lexicographic (orb1,orb2,kx,ky,kz,iomega) "
                     "tie-break)")

    np.savez(
        os.path.join(output_dir, npz_file),
        gap=gap_w,
        iomega=iomega,
        T=T,
        pairing_type=pairing_type,
        frequency="dynamic",
        eigenvalue=eigenvalue,
        axis_order=axis_order,
        normalization=normalization,
    )

    # gap.dat: the fermionic slice nearest omega = 0^+ (smallest positive w_n).
    n0 = nmat // 2
    logger.info("Saving dynamic gap slice (index %d, iw_n=%.6e) to %s",
                n0, iomega[n0], os.path.join(output_dir, gap_file))
    with open(os.path.join(output_dir, gap_file), "w") as fw:
        header = ["# frequency=dynamic",
                  "index={}".format(n0),
                  "iomega_n={:.8e}".format(iomega[n0]),
                  "pairing_type={}".format(pairing_type),
                  "eigenvalue={:.8e}".format(eigenvalue),
                  "T={:.8e}".format(T)]
        fw.write("  ".join(header) + "\n")
        cols = ["# kx", "ky", "kz"]
        for i in range(norb):
            for j in range(norb):
                cols.append("Re(sigma_{}{})".format(i, j))
                cols.append("Im(sigma_{}{})".format(i, j))
        fw.write(" ".join(cols) + "\n")
        for ix in range(Nx):
            kx = kx_array[ix]
            if kx > np.pi:
                kx -= 2.0 * np.pi
            for iy in range(Ny):
                ky = ky_array[iy]
                if ky > np.pi:
                    ky -= 2.0 * np.pi
                for iz in range(Nz):
                    kz = kz_array[iz]
                    if kz > np.pi:
                        kz -= 2.0 * np.pi
                    parts = ["{:.8f}".format(kx), "{:.8f}".format(ky),
                             "{:.8f}".format(kz)]
                    for i in range(norb):
                        for j in range(norb):
                            val = gap_w[i, j, ix, iy, iz, n0]
                            parts.append("{:.8e}".format(val.real))
                            parts.append("{:.8e}".format(val.imag))
                    fw.write(" ".join(parts) + "\n")


def solve_dynamic(input_dict):
    """Solve the dynamic (frequency-resolved) Eliashberg equation.

    Reads the FLEX outputs (full-frequency chi_s/chi_c and the dressed
    Green's function), builds the frequency-resolved pairing vertex and pair
    bubble, and finds the leading eigenpair of the tau-product Eliashberg
    kernel via the shared driver ``sc._solve_leading`` (the same eigenvalue
    ordering / shift-invert / iteration machinery as the static path).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Returns
    -------
    float
        The leading (largest real part) Eliashberg eigenvalue lambda.
    """
    import hwave.sc as sc
    from scipy.sparse.linalg import LinearOperator

    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = 1.0 / T
    cell_shape = mode_param["CellShape"]
    sub_shape = mode_param.get("SubShape", cell_shape)
    if isinstance(cell_shape, list):
        cell_shape = list(cell_shape)
        while len(cell_shape) < 3:
            cell_shape.append(1)
    Lx, Ly, Lz = cell_shape
    if isinstance(sub_shape, list):
        sub_shape = list(sub_shape)
        while len(sub_shape) < 3:
            sub_shape.append(1)
    Bx, By, Bz = sub_shape
    Nx, Ny, Nz = Lx // Bx, Ly // By, Lz // Bz
    Nk = Nx * Ny * Nz

    eli_param = input_dict.get("eliashberg", {})
    pairing_type = eli_param.get("pairing_type", "singlet")
    # Mirror calc_eliashberg's config -> _solve_leading string mapping.
    solver_mode = eli_param.get("solver_mode", "iteration")
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)
    init_gap_mode = sc._resolve_init_gap(eli_param.get("init_gap"), pairing_type)

    # --- Geometry / interactions (norb from the geometry file) ---
    geom_info, hr, interactions = sc._read_interaction_files(input_dict)
    norb = geom_info["norb"]

    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)
    inter_k = sc._build_interaction_k(kx_array, ky_array, kz_array,
                                      interactions, norb)

    # --- FLEX inputs (full frequency) ---
    chis_w, chic_w, green_w, chi_convention = load_flex_chi_dynamic(
        input_dict, norb, Nx, Ny, Nz)
    if green_w is None:
        raise ValueError(
            "dynamic Eliashberg requires the dressed green.npz from the FLEX "
            "run (the pair bubble G2 is built from it); none was found. Check "
            "[file.input] path_to_flex_output / [eliashberg] flex_green.")
    nmat = chis_w.shape[-1]

    # --- Vertex and pair bubble on the full frequency grid ---
    logger.info("Computing dynamic FLEX pairing vertex (pairing_type=%s, "
                "convention=%s)...", pairing_type, chi_convention)
    Vs_q_w = compute_vertices_flex_dynamic(
        chis_w, chic_w, inter_k, norb, Nx, Ny, Nz,
        pairing_type=pairing_type, convention=chi_convention)
    logger.info("Computing frequency-resolved pair bubble G2...")
    G2_w = calc_g2_dynamic(green_w, beta)

    # --- Seed: static init_gap form factor, broadcast flat across omega ---
    gap_shape = (norb, norb, Nx, Ny, Nz, nmat)
    sigma_static = sc._initialize_gap(init_gap_mode, norb,
                                      kx_array, ky_array, kz_array)
    phi0 = np.broadcast_to(sigma_static[..., np.newaxis], gap_shape).copy()
    phi0 = phi0.astype(complex)
    n0 = np.linalg.norm(phi0)
    if n0 > 0:
        phi0 /= n0

    vec_size = norb * norb * Nk * nmat
    assert phi0.size == vec_size

    def _matvec(x):
        return eliashberg_kernel_dynamic(
            Vs_q_w, G2_w, x.reshape(gap_shape), norb, beta).ravel()

    def make_operator():
        op = LinearOperator((vec_size, vec_size), matvec=_matvec,
                            dtype=complex)
        return op, vec_size

    # Map [eliashberg] controls to the _solve_leading solver_mode string,
    # exactly as calc_eliashberg does for the static path.
    if solver_mode == "iteration":
        eigenvalue, sigma_flat, info = sc._solve_leading(
            make_operator, vec_size, "iteration",
            max_iter=max_iter, convergence_tol=tol, alpha=alpha,
            init_vec=phi0.ravel())
        eigenvalues_all = None
    else:
        # "eigenvalue" / "both": use the ARPACK/shift-invert eigen family.
        eigenvalue, sigma_flat, info = sc._solve_leading(
            make_operator, vec_size, eigenvalue_method,
            num_eigenvalues=num_eigenvalues)
        eigenvalues_all = info.get("eigenvalues")

    lam = float(np.real(eigenvalue))
    logger.info("Dynamic Eliashberg leading eigenvalue lambda = %.6f", lam)

    # --- Outputs ---
    # Gauge-fix the eigenvector (deterministic phase/normalization) so the
    # written gap is reproducible across runs and linear-algebra backends.
    gap_w = _fix_gauge(np.asarray(sigma_flat).reshape(gap_shape))
    output_dir = input_dict["file"]["output"]["path_to_output"]
    os.makedirs(output_dir, exist_ok=True)
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")
    with open(os.path.join(output_dir, eigenvalue_file), "w") as fw:
        fw.write("# Dynamic Eliashberg leading eigenvalue\n")
        fw.write("{:.8e}\n".format(lam))
        if eigenvalues_all is not None:
            fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  |eigenvalue|\n")
            for i, ev in enumerate(eigenvalues_all):
                fw.write("{:4d} {:15.8e} {:15.8e} {:15.8e}\n".format(
                    i, ev.real, ev.imag, abs(ev)))

    gap_file = eli_param.get("output_gap", "gap.dat")
    write_dynamic_outputs(output_dir, gap_w, lam, T, pairing_type,
                          kx_array, ky_array, kz_array, beta,
                          gap_file=gap_file)

    return lam
