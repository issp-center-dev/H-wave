"""The dynamic (frequency-resolved) Eliashberg solver in spin-orbital mode
(``[mode] enable_spin_orbital = true`` with ``[eliashberg] frequency =
"dynamic"``, issue #284).

The gap is one spin-orbital matrix per momentum and fermionic Matsubara
frequency, ``Delta_AB(k, i w_n)`` over ``nd = 2 N`` spin-orbitals in
spin-block order (``A = s * N + a``), on the uniform centered grid
``w_n = (2 n + 1 - Nmat) pi T``. Fermion antisymmetry requires

    Delta_AB(k, i w_n) = -Delta_BA(-k, -i w_n),

where ``-i w_n`` is the index ``Nmat - 1 - n`` (the map
``eliashberg_dynamic._reverse_kw_and_orbital``). A fermionic frequency is
never zero, so no element is paired with itself and the antisymmetric space
has exactly half the dimension of the gap.

Pipeline (every piece below is shared with an existing, validated path):

* the pairing vertex at every bosonic frequency ``nu_m = (2 m - Nmat) pi T``,
  ``V(q, i nu_m) = 1/2 Gamma + (-Gamma chi(q, i nu_m) Gamma)``
  (:func:`hwave.solver.so_pairing.so_pairing_vertex_dynamic`), with ``Gamma``
  the antisymmetrized spinful vertex of the interaction files and ``chi`` the
  spinful susceptibility of a spin-orbital FLEX run (``chi0q_mode = "flex"``)
  or ``[1 + chi0 Gamma]^-1 chi0`` of the bare spinful bubble (``"calc"``);
* the frequency-resolved pair bubble ``calc_g2_dynamic`` and the tau-product
  kernel ``eliashberg_kernel_dynamic`` of the non-spin-orbital dynamic solver,
  at ``nd`` spin-orbitals (the kernel is generic in the orbital count; with
  this vertex it equals, to round-off, the ``BondPairKernel`` application
  that ``tests/test_so_pairing_ed.py`` validates against exact
  diagonalization with spin-orbit coupling);
* the eigenproblem on the fermion-antisymmetric coordinates ``E^dagger K E``
  (``sc._SOAntisymmetricBasis`` with the frequency-reversing map) through the
  drivers of the static spin-orbital solver (``sc._so_solve_iteration`` /
  ``sc._so_solve_eigenvalue``, including the dense full-spectrum branch).

Frequency parity ``Delta_AB(k, i w_n) -> Delta_AB(k, -i w_n)`` is not a
symmetry of the spin-orbital kernel in general, so the whole antisymmetric
space (even- and odd-frequency gaps together) is solved; the even/odd-
frequency content of each gap is reported, not selected.

Scope: the uniform grid, ``chi0q_mode = "flex"`` or ``"calc"``, on-site
interactions, no bond channels (the refusals are in
:func:`hwave.solver.so_pairing.check_so_config`).
"""
import logging
import os
import time

import numpy as np
from requests.structures import CaseInsensitiveDict
from scipy.sparse.linalg import LinearOperator

from hwave.solver import backend
from hwave.solver.eliashberg_dynamic import (
    _available_ram_bytes, _fix_gauge, _npz_freq_size, _reverse_kw_and_orbital,
    _reverse_w_dynamic, calc_g2_dynamic, eliashberg_kernel_dynamic)

logger = logging.getLogger("qlms").getChild("eliashberg_dynamic_so")

#: default ``init_gap`` of the spin-orbital dynamic solver: a fixed-seed
#: random gap over (A, B, k, w_n). A frequency-flat seed is purely even in
#: frequency, and frequency parity is not a symmetry of the spin-orbital
#: kernel in general, so a power iteration started there could stay in the
#: even-frequency sector of a kernel that happens to conserve it.
DEFAULT_INIT_GAP = "random"
_RANDOM_GAP_SEED = 12345

#: Above this, ``pauli_leakage`` (the symmetric part of ``K`` applied to
#: antisymmetric probes) is reported with a warning. On the uniform grid the
#: kernel commutes with the Pauli map to round-off (~1e-15); the solve itself
#: is restricted to the antisymmetric space either way.
PAULI_LEAKAGE_WARN = 1.0e-8

#: labels of the two frequency-parity sectors, in output order
FREQUENCY_PARITY_LABELS = ("even_w", "odd_w")

_OUTPUT_NPZ = "gap_dynamic.npz"


# ---------------------------------------------------------------------------
# Antisymmetry and the weights of a dynamic spin-orbital gap
# ---------------------------------------------------------------------------

def antisymmetrize(gap_w):
    """``(Delta - R Delta) / 2`` with ``R Delta_AB(k, n) = Delta_BA(-k,
    nmat - 1 - n)``, on ``(nd, nd, Nx, Ny, Nz, nmat)``."""
    gap_w = np.asarray(gap_w)
    return 0.5 * (gap_w - _reverse_kw_and_orbital(gap_w))


def antisymmetry_residual(gap_w):
    """``||Delta + R Delta|| / (2 ||Delta||)`` (the symmetric part relative to
    the gap); ``inf`` for a zero or non-finite gap. Compare with
    ``sc._SO_ANTISYM_TOL``."""
    gap_w = np.asarray(gap_w)
    n = float(np.linalg.norm(gap_w))
    if n == 0.0 or not np.isfinite(n):
        return float("inf")
    return float(np.linalg.norm(0.5 * (gap_w + _reverse_kw_and_orbital(gap_w)))) / n


def gap_weights(gap_w, norb_phys):
    """Spin and frequency-parity content of a dynamic spin-orbital gap.

    With ``Delta_ab(k, i w_n) = (psi_ab + d_ab . sigma) i sigma_y`` per
    orbital pair (``so_pairing.so_spin_components``) and the frequency-parity
    halves ``Delta^(+/-) = (Delta(i w_n) +/- Delta(-i w_n)) / 2``, returns a
    dict of

    * ``spin`` ``(4,)``: ``w_X = sum_{n,k,a,b} |X_ab(k, n)|^2 / W``, X = psi,
      dx, dy, dz, ``W`` the sum over all four (the static #269 weights, now
      summed over frequency too);
    * ``spin_w`` ``(nmat, 4)``: the same per frequency, normalized by the
      same ``W`` (the whole array sums to 1);
    * ``parity`` ``(2,)``: ``|Delta^(+)|^2 / |Delta|^2`` and
      ``|Delta^(-)|^2 / |Delta|^2`` (even and odd frequency);
    * ``spin_parity`` ``(4, 2)``: ``w_X`` of the even and odd halves (sums
      to 1).

    On the antisymmetric space the momentum parity is minus the frequency
    parity, so these two sectors are the only ones. In the SU(2) limit
    ``(psi, even)`` is the even-frequency singlet, ``(d, even)`` the
    even-frequency triplet, ``(psi, odd)`` and ``(d, odd)`` the odd-frequency
    triplet-parity singlet and singlet-parity triplet. NaN for a zero gap."""
    from hwave.solver.so_pairing import so_spin_components
    gap_w = np.asarray(backend.to_host(gap_w))
    nmat = gap_w.shape[-1]
    comps = so_spin_components(gap_w, norb_phys)
    per = np.array([[float(np.sum(np.abs(c[..., n]) ** 2)) for c in comps]
                    for n in range(nmat)])                      # (nmat, 4)
    tot = float(per.sum())
    rev = _reverse_w_dynamic(gap_w)
    halves = (0.5 * (gap_w + rev), 0.5 * (gap_w - rev))
    sp = np.array([[float(np.sum(np.abs(c) ** 2))
                    for c in so_spin_components(h, norb_phys)] for h in halves]).T
    norm2 = float(np.linalg.norm(gap_w) ** 2)
    if not (tot > 0.0 and norm2 > 0.0):
        nan = float("nan")
        return {"spin": np.full(4, nan), "spin_w": np.full((nmat, 4), nan),
                "parity": np.full(2, nan), "spin_parity": np.full((4, 2), nan)}
    return {"spin": per.sum(axis=0) / tot,
            "spin_w": per / tot,
            "parity": np.array([float(np.linalg.norm(h) ** 2) / norm2
                                for h in halves]),
            "spin_parity": sp / float(sp.sum())}


# ---------------------------------------------------------------------------
# Memory
# ---------------------------------------------------------------------------

def estimate_memory_bytes(nd, Nk, nmat, solver_mode="iteration",
                          num_eigenvalues=10, chi0q_mode="flex"):
    """Estimated peak host memory (bytes) of a spin-orbital dynamic solve.

    With ``T = 16 nd^4 Nk nmat`` (one four-leg tensor: chi, chi0, the vertex,
    its (r, tau) form, G2), ``g = 16 nd^2 Nk nmat`` (the Green function, one
    gap vector), ``n = nd^2 Nk nmat`` and ``dim = n / 2``, the stages are
    (each measured with ``tracemalloc``, issue #284)

    * input: chi + chi0 + green, ``2 T + g`` (``"calc"``: the RPA bubble
      routine peaks at ``4.5 T``);
    * vertex: chi + V, ``2 T + g``;
    * vertex transform (slab by slab, :func:`vertex_to_rt`): V + V(r, tau)
      + slab buffers, ``2 T + 3 T / nd^2 + g``;
    * pair bubble: V(r, tau) + G2 (scaled in place) + G(-k), ``2 T + 2 g``;
    * solve: V(r, tau) + G2 + the eigensolver workspace (ARPACK: ``(ncv + 4)
      dim`` complex vectors, ``ncv = max(2 k + 1, 20)``; the dense branch for
      ``k >= dim - 1``: ``3 dim^2``; the power iteration: a few vectors) + the
      basis index arrays (``24 n`` bytes) + six gap-sized kernel temporaries.

    The estimate is the largest stage plus 10 % for allocator slack."""
    nd, Nk, nmat = int(nd), int(Nk), int(nmat)
    T = 16 * nd ** 4 * Nk * nmat
    g = 16 * nd ** 2 * Nk * nmat
    n = nd ** 2 * Nk * nmat
    dim = n // 2
    k = int(num_eigenvalues)
    if solver_mode in ("eigenvalue", "both"):
        if k >= dim - 1:
            work = 3 * 16 * dim * dim
        else:
            ncv = max(2 * k + 1, 20)
            work = 16 * (ncv + 4) * dim + 3 * 16 * ncv * ncv
    else:
        work = 16 * 4 * dim
    stages = (
        (4.5 * T if chi0q_mode == "calc" else 2 * T) + g,
        2 * T + g,
        2 * T + 3 * T / nd ** 2 + g,
        2 * T + 2 * g,
        2 * T + work + 24 * n + 6 * g,
    )
    return int(1.1 * max(stages))


def check_memory(need, mem_limit_gb):
    """Refuse a run whose estimate exceeds ``[eliashberg] mem_limit_gb``
    (default: 80 % of the available RAM; 0 disables the guard)."""
    if mem_limit_gb == 0:
        logger.info("spin-orbital dynamic Eliashberg est. peak %.2f GiB (guard "
                    "disabled)", need / 2 ** 30)
        return
    limit = ((0.8 * _available_ram_bytes()) if mem_limit_gb is None
             else float(mem_limit_gb) * 2 ** 30)
    logger.info("spin-orbital dynamic Eliashberg est. peak %.2f GiB (limit %.2f "
                "GiB)", need / 2 ** 30, limit / 2 ** 30)
    if need > limit:
        raise MemoryError(
            "[mode] enable_spin_orbital = true with frequency = \"dynamic\": the "
            "estimated peak memory {:.1f} GiB exceeds the limit {:.1f} GiB "
            "([eliashberg] mem_limit_gb; default 80% of the available RAM). Every "
            "four-leg tensor scales as nd^4 Nk Nmat with nd = 2 x the physical "
            "orbitals: reduce Nmat or the k grid, or raise mem_limit_gb (0 "
            "disables the guard).".format(need / 2 ** 30, limit / 2 ** 30))


# ---------------------------------------------------------------------------
# The operator
# ---------------------------------------------------------------------------

def vertex_to_rt(V_w, workers=1):
    """``vertex_qw_to_rt`` (bosonic frequency -> tau, then q -> r) of the
    vertex ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)``, one ``(a, b)`` slab at a
    time: the transforms act on the last four axes only, so this is the same
    map with a working set of one slab (``1 / nd^2`` of the tensor) instead of
    two to three full-size FFT buffers. Same backend as ``V_w``."""
    from hwave.solver import matsubara as ms
    from hwave.solver.backend import spatial_ifftn
    xp = backend.array_module_of(V_w)
    out = xp.empty(V_w.shape, dtype=complex)
    nd = V_w.shape[0]
    for a in range(nd):
        for b in range(V_w.shape[1]):
            out[a, b] = spatial_ifftn(ms.boson_to_tau(V_w[a, b], axis=-1),
                                      axes=(2, 3, 4), workers=workers)
    return out


def make_operator(Vs_rt, G2_w, beta, workers=1):
    """``(A, basis)``: the dynamic kernel as a flat ``LinearOperator`` on
    ``(nd, nd, Nx, Ny, Nz, nmat)`` gaps (``eliashberg_kernel_dynamic`` with the
    hoisted vertex ``Vs_rt``; a cupy ``G2_w`` / ``Vs_rt`` runs it on the
    device, the vector crosses as a host array) and the fermion-antisymmetric
    basis of that space with the frequency-reversing Pauli map."""
    import hwave.sc as sc
    nd = int(G2_w.shape[0])
    gap_shape = (nd, nd) + tuple(int(x) for x in G2_w.shape[4:8])
    n = int(np.prod(gap_shape))

    def _matvec(x):
        out = eliashberg_kernel_dynamic(
            None, G2_w, np.asarray(x).reshape(gap_shape), nd, beta, Vs_rt=Vs_rt,
            workers=workers)
        return np.asarray(backend.to_host(out)).ravel()

    A = LinearOperator((n, n), matvec=_matvec, dtype=complex)
    basis = sc._SOAntisymmetricBasis(gap_shape, reverse=_reverse_kw_and_orbital)
    return A, basis


def build_operator(gam, chi_w, green_w, beta, xp=np, workers=1):
    """The production operator from its inputs (the steps of
    :func:`solve_dynamic_so`, without the memory ordering): the vertex
    (:func:`hwave.solver.so_pairing.so_pairing_vertex_dynamic`) of ``gam``
    ``(nvol, nd^2, nd^2)`` and ``chi_w`` ``(nmat, nvol, nd^2, nd^2)``, its
    (r, tau) form, the pair bubble of ``green_w`` ``(nd, nd, Nx, Ny, Nz,
    nmat)``, then :func:`make_operator`. ``xp = cupy`` keeps the vertex and
    the pair bubble on the device."""
    from hwave.solver.so_pairing import so_pairing_vertex_dynamic
    nd = int(green_w.shape[0])
    V_w = so_pairing_vertex_dynamic(gam, chi_w, nd, tuple(green_w.shape[2:5]))
    Vs_rt = vertex_to_rt(xp.asarray(V_w), workers=workers)
    del V_w
    G2_w = xp.asarray(calc_g2_dynamic(np.asarray(green_w), beta))
    return make_operator(Vs_rt, G2_w, beta, workers=workers)


def pauli_leakage(A, basis, n_probe=3, seed=0):
    """``max ||sym(K E c)|| / ||K E c||`` over ``n_probe`` fixed-seed random
    reduced vectors ``c``: how far the kernel maps antisymmetric gaps out of
    the antisymmetric space (zero when ``K`` commutes with the Pauli map; the
    solve is restricted to the antisymmetric space regardless)."""
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(int(n_probe)):
        c = rng.standard_normal(basis.dim) + 1j * rng.standard_normal(basis.dim)
        y = np.asarray(A.matvec(basis.to_full(c))).reshape(basis.shape)
        ny = float(np.linalg.norm(y))
        if ny == 0.0:
            continue
        sym = 0.5 * (y + _reverse_kw_and_orbital(y))
        worst = max(worst, float(np.linalg.norm(sym)) / ny)
    return worst


# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------

def initial_gap(init_gap, norb_phys, kx, ky, kz, gap_shape):
    """The (normalized, antisymmetric) initial gap of the power iteration.

    ``"random"`` (the default, :data:`DEFAULT_INIT_GAP`): a fixed-seed random
    gap over every (A, B, k, w_n), projected onto the antisymmetric space --
    both frequency parities. A named form factor: the static spin-orbital seed
    (``so_pairing.so_initial_gap``) broadcast flat over frequency (an
    even-frequency gap)."""
    from hwave.solver.so_pairing import so_initial_gap
    if init_gap == "random":
        rng = np.random.default_rng(_RANDOM_GAP_SEED)
        g = rng.standard_normal(gap_shape) + 1j * rng.standard_normal(gap_shape)
    else:
        g0 = so_initial_gap(init_gap, norb_phys, kx, ky, kz)
        g = np.broadcast_to(g0[..., None], gap_shape).astype(complex)
    g = antisymmetrize(g)
    nrm = float(np.linalg.norm(g))
    if nrm == 0.0:
        raise ValueError(
            "[eliashberg] init_gap = {!r} has no fermion-antisymmetric component "
            "on this grid; choose another init_gap (the spin-orbital dynamic "
            "default is \"random\").".format(init_gap))
    return g / nrm


def load_seed(path, gap_shape):
    """The flat, antisymmetrized eigenvector-continuation seed of
    ``[eliashberg] seed_eigenvector``: a ``gap_dynamic.npz`` of a
    spin-orbital dynamic run on the same grid. A gap without the
    ``spin_orbital = True`` marker is refused (a non-spin-orbital gap over 2N
    orbitals has the same shape)."""
    from hwave.solver import spin_orbital
    from hwave.solver.rpa import (check_momentum_marker,
                                  validate_momentum_convention)
    with np.load(path) as data:
        check_momentum_marker(data, path)
        mark = np.asarray(data["spin_orbital"]) if "spin_orbital" in data else None
        if mark is None or mark.ndim != 0 or mark.dtype != np.bool_ or not bool(mark):
            raise ValueError(
                "seed_eigenvector '{}' is not the gap of a spin-orbital run (no "
                "spin_orbital = True marker); a spin-orbital dynamic run is seeded "
                "only by the gap_dynamic.npz of another one.".format(path))
        if ("index_convention" in data and str(data["index_convention"])
                != spin_orbital.INDEX_CONVENTION):
            raise ValueError(
                "seed_eigenvector '{}' carries index_convention = {!r}, expected "
                "{!r}.".format(path, str(data["index_convention"]),
                               spin_orbital.INDEX_CONVENTION))
        gap = np.asarray(data["gap"])
        if gap.ndim != 6:
            raise ValueError(
                "seed_eigenvector '{}' has ndim = {}; a dynamic gap is (nd, nd, Nx, "
                "Ny, Nz, Nmat).".format(path, gap.ndim))
        validate_momentum_convention(data, path, gap, (2, 3, 4),
                                     tuple(gap.shape[2:5]))
    if gap.shape != tuple(gap_shape):
        raise ValueError(
            "seed_eigenvector '{}' has shape {}, this run {}; eigenvector "
            "continuation needs the same orbitals, k grid and Nmat.".format(
                path, gap.shape, tuple(gap_shape)))
    return np.ascontiguousarray(antisymmetrize(gap.astype(complex))).ravel()


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def _grid(mode_param):
    cell = mode_param["CellShape"]
    cs = list(cell) if isinstance(cell, (list, tuple)) else [cell]
    while len(cs) < 3:
        cs.append(1)
    ss = list(mode_param.get("SubShape", cs))
    while len(ss) < 3:
        ss.append(1)
    if ss != [1, 1, 1]:
        raise ValueError(
            "SubShape (sublattice folding) is not supported by the "
            "Eliashberg module: fold the model into the unit cell "
            "yourself, or set SubShape = [1, 1, 1] explicitly. (Note: "
            "omitting SubShape defaults it to CellShape -- the whole "
            "cell as one supercell -- per the package convention, so it "
            "must be set explicitly here.)")
    return tuple(int(x) for x in cs)


def _check_vertex_every_frequency(gam, chi_w, chi0_w, chi_path, T):
    """``sc._check_so_flex_vertex`` at every bosonic frequency: the stored
    chi must solve ``(1 + chi0 Gamma) chi = chi0`` with the stored chi0 and
    the Gamma of this run's interaction files, frequency by frequency."""
    import hwave.sc as sc
    nmat = chi_w.shape[0]
    errs = np.array([sc._so_flex_backward_error(gam, chi_w[m], chi0_w[m])
                     for m in range(nmat)])
    bad = np.where(~(errs <= sc._SO_FLEX_GAMMA_TOL))[0]
    worst = float(np.nanmax(errs)) if np.any(np.isfinite(errs)) else float("nan")
    logger.info("spin-orbital FLEX: vertex consistency backward error, max over "
                "%d bosonic frequencies %.3e", nmat, worst)
    if bad.size:
        m = int(bad[np.argmax(np.nan_to_num(errs[bad], nan=np.inf))])
        raise ValueError(
            "the interaction files of this run do not reproduce the vertex of "
            "the FLEX run that wrote '{}' (backward error {:.3e} > {:.1e} of "
            "(1 + chi0 Gamma) chi = chi0 at bosonic index {}, nu = {:.6g}; {} of "
            "{} bosonic frequencies fail). Use the same interaction files and "
            "[mode.param] as the FLEX run.".format(
                chi_path, errs[m], sc._SO_FLEX_GAMMA_TOL, m,
                (2 * m - nmat) * np.pi * T, bad.size, nmat))
    return worst


# ---------------------------------------------------------------------------
# Outputs
# ---------------------------------------------------------------------------

def _gap_layout(norb_phys):
    return ("Delta_AB(k, iw_n), A = spin * {} + orb (spin 0 = up, spin 1 = down), "
            "iw_n = (2n + 1 - Nmat) pi T; fermion-antisymmetric, Delta_AB(k, iw_n) "
            "= -Delta_BA(-k, -iw_n) (-iw_n is index Nmat - 1 - n)".format(norb_phys))


def _fmt_weights(labels, w):
    return " ".join("{}={:.8e}".format(k, float(v)) for k, v in zip(labels, w))


def write_outputs(output_dir, gap_file, eigenvalue_file, *, gap_w, lam,
                  eigenvalues, eigen_weights, match, norb_phys, T, kx, ky, kz,
                  eigenvalue_selection, leakage, note=None, extra_prov=None):
    """Write ``eigenvalue.dat``, ``gap.dat`` and ``gap_dynamic.npz`` in the
    spin-orbital schema of the static solver (#269) plus the frequency
    content (see the user guide)."""
    from hwave.solver import so_pairing, spin_orbital
    os.makedirs(output_dir, exist_ok=True)
    comps = so_pairing.SPIN_COMPONENTS
    nd = gap_w.shape[0]
    Nx, Ny, Nz, nmat = gap_w.shape[2:6]
    iomega = (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi * T
    w = gap_weights(gap_w, norb_phys)
    layout = _gap_layout(norb_phys)

    prov = {
        "spin_orbital": "true",
        "index_convention": spin_orbital.INDEX_CONVENTION,
        "gap_layout": layout,
        "gap_npz": _OUTPUT_NPZ,
        "match_column": "1 = fermion-antisymmetric eigenvector",
        "spin_decomposition": (
            "Delta_ab(k, iw_n) = (psi_ab + d_ab . sigma) i sigma_y per orbital "
            "pair (a, b); w_X = sum_{n,k,a,b} |X_ab(k, iw_n)|^2 / sum_{n,k,a,b} "
            "(|psi_ab|^2 + |d_ab|^2), X = psi, dx, dy, dz"),
        "frequency_parity": (
            "w_even_w, w_odd_w = |Delta(iw_n) +/- Delta(-iw_n)|^2 / (4 |Delta|^2); "
            "frequency parity is not a symmetry of the spin-orbital kernel in "
            "general, so both sectors are solved together"),
        "weight_columns": ("w_psi w_dx w_dy w_dz w_even_w w_odd_w (after the "
                           "match column)"),
        "gap_spin_weights": _fmt_weights(comps, w["spin"]),
        "gap_frequency_parity_weights": _fmt_weights(FREQUENCY_PARITY_LABELS,
                                                     w["parity"]),
        "gap_spin_frequency_weights": " ".join(
            "{}_{}={:.8e}".format(c, p, w["spin_parity"][i, j])
            for i, c in enumerate(comps)
            for j, p in enumerate(FREQUENCY_PARITY_LABELS)),
        "pauli_leakage": "{:.3e}".format(leakage),
    }
    prov.update(extra_prov or {})
    with open(os.path.join(output_dir, eigenvalue_file), "w") as fw:
        fw.write("# Dynamic Eliashberg leading eigenvalue\n")
        for key in sorted(prov):
            fw.write("# {} = {}\n".format(key, prov[key]))
        fw.write("# eigenvalue_selection: {}\n".format(eigenvalue_selection))
        if note:
            for line in str(note).splitlines():
                fw.write("# {}\n".format(line))
        fw.write("{:.8e}\n".format(lam))
        if eigenvalues is not None:
            fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  |eigenvalue|  "
                     "match(1=fermion-antisymmetric)  w_psi  w_dx  w_dy  w_dz  "
                     "w_even_w  w_odd_w\n")
            for i, ev in enumerate(eigenvalues):
                row = "{:4d} {:15.8e} {:15.8e} {:15.8e} {:d}".format(
                    i, ev.real, ev.imag, abs(ev), int(bool(match[i])))
                row += "".join(" {:15.8e}".format(x) for x in eigen_weights[i])
                fw.write(row + "\n")

    record = {
        "gap": gap_w, "iomega": iomega, "T": T, "frequency": "dynamic",
        "eigenvalue": lam,
        "axis_order": "(A, B, kx, ky, kz, iomega), A = spin * norb_phys + orb",
        "normalization": ("L2-normalized over all (A,B,k,iomega) components; "
                          "global phase fixes the largest-|magnitude| component "
                          "real-positive (lexicographic tie-break)"),
        "momentum_convention": "e_plus_ikR",
        "index_convention": np.array(spin_orbital.INDEX_CONVENTION),
        "spin_orbital": np.array(True),
        "norb_phys": np.array(int(norb_phys)),
        "kx": np.asarray(kx), "ky": np.asarray(ky), "kz": np.asarray(kz),
        "spin_components": np.array(comps),
        "gap_spin_weights": w["spin"],
        "gap_spin_weights_w": w["spin_w"],
        "gap_frequency_parity_labels": np.array(FREQUENCY_PARITY_LABELS),
        "gap_frequency_parity_weights": w["parity"],
        "gap_spin_frequency_weights": w["spin_parity"],
        "eigenvalue_selection": str(eigenvalue_selection),
        "pauli_leakage": float(leakage),
    }
    if eigenvalues is not None:
        ew = np.asarray(eigen_weights)
        record.update({
            "eigenvalues": np.asarray(eigenvalues),
            "eigenvector_spin_weights": ew[:, :4],
            "eigenvector_frequency_parity_weights": ew[:, 4:6],
            "eigenvalue_match": np.asarray(match, dtype=bool)})
    np.savez(os.path.join(output_dir, _OUTPUT_NPZ), **record)

    n0 = nmat // 2
    with open(os.path.join(output_dir, gap_file), "w") as fw:
        fw.write("  ".join(["# frequency=dynamic", "index={}".format(n0),
                            "iomega_n={:.8e}".format(iomega[n0]),
                            "eigenvalue={:.8e}".format(lam),
                            "T={:.8e}".format(T)]) + "\n")
        fw.write("# spin_orbital = true\n")
        fw.write("# index_convention = {}\n".format(spin_orbital.INDEX_CONVENTION))
        fw.write("# gap_layout = {}; sigma_ij below is Delta_ij at iw_n\n".format(
            layout))
        cols = ["# kx", "ky", "kz"]
        for i in range(nd):
            for j in range(nd):
                cols += ["Re(sigma_{}{})".format(i, j), "Im(sigma_{}{})".format(i, j)]
        fw.write(" ".join(cols) + "\n")
        for ix in range(Nx):
            kxv = kx[ix] - 2.0 * np.pi if kx[ix] > np.pi else kx[ix]
            for iy in range(Ny):
                kyv = ky[iy] - 2.0 * np.pi if ky[iy] > np.pi else ky[iy]
                for iz in range(Nz):
                    kzv = kz[iz] - 2.0 * np.pi if kz[iz] > np.pi else kz[iz]
                    parts = ["{:.8f}".format(v) for v in (kxv, kyv, kzv)]
                    for i in range(nd):
                        for j in range(nd):
                            val = gap_w[i, j, ix, iy, iz, n0]
                            parts += ["{:.8e}".format(val.real),
                                      "{:.8e}".format(val.imag)]
                    fw.write(" ".join(parts) + "\n")
    return w


# ---------------------------------------------------------------------------
# The entry
# ---------------------------------------------------------------------------

def solve_dynamic_so(input_dict):
    """Solve the spin-orbital dynamic Eliashberg equation (issue #284).

    Reached from ``sc.calc_eliashberg`` (and from a direct
    ``eliashberg_dynamic.solve_dynamic`` call) when ``[mode]
    enable_spin_orbital`` is set and ``[eliashberg] frequency = "dynamic"``.
    Returns the leading eigenvalue (the real part of the eigenvalue analysis'
    first eigenvalue, or the power iteration's value)."""
    import hwave.sc as sc
    from hwave.solver import so_pairing, spin_orbital

    t_start = time.time()
    # this IS the dynamic entry (also when eliashberg_dynamic.solve_dynamic is
    # called directly, without a frequency key): check the configuration as a
    # dynamic one, before anything is read
    eli_raw = {k: v for k, v in (input_dict.get("eliashberg", {}) or {}).items()
               if str(k).lower() != "frequency"}
    eli_raw["frequency"] = "dynamic"
    so_pairing.check_so_config(dict(input_dict, eliashberg=eli_raw))
    eli = CaseInsensitiveDict(eli_raw)
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = sc._coerce_run_beta(T)
    Nx, Ny, Nz = _grid(mode_param)
    Nk = Nx * Ny * Nz
    nmat = int(mode_param.get("Nmat", sc._DEFAULT_NMAT))
    if nmat <= 0 or nmat % 2 != 0:
        raise ValueError(
            "[mode] enable_spin_orbital = true with frequency = \"dynamic\" "
            "requires an even, positive [mode.param] Nmat; got {}".format(nmat))
    chi0q_mode = eli.get("chi0q_mode")
    solver_mode = eli.get("solver_mode", "iteration")
    num_eigenvalues = int(eli.get("num_eigenvalues", 10))
    eigenvalue_method = eli.get("eigenvalue_method", "arnoldi")
    spectral_shift = eli.get("spectral_shift")
    init_gap_mode = eli.get("init_gap", DEFAULT_INIT_GAP)
    fft_workers = eli.get("fft_workers", 1)
    for key in ("ir_wmax", "ir_tol", "ir_keep_static_chi"):
        if key in eli:
            logger.warning("[eliashberg] %s is not used on the uniform grid "
                           "(matsubara_basis = \"uniform\"); it is ignored.", key)
    gap_file = eli.get("output_gap", "gap.dat")
    eigenvalue_file = eli.get("output_eigenvalue", "eigenvalue.dat")
    for name, what in ((gap_file, "output_gap"), (eigenvalue_file, "output_eigenvalue")):
        if os.path.normpath(name) == os.path.normpath(_OUTPUT_NPZ):
            raise ValueError(
                "[eliashberg] {} = {!r} would overwrite the spin-orbital dynamic "
                "gap archive '{}'; choose another name.".format(what, name,
                                                                _OUTPUT_NPZ))
    output_dir = input_dict["file"]["output"]["path_to_output"]

    # --- model (no heavy work yet) ---
    geom_info, hr, interactions = sc._read_interaction_files(input_dict)
    nd = geom_info["norb"]
    norb_phys = spin_orbital.physical_norb(
        nd, True, source="the Geometry file ({})".format(
            input_dict["file"]["input"]["interaction"].get("Geometry", "")))
    so_pairing.check_so_supported(
        input_dict, interactions, sc._nonbond_site(input_dict, norb_phys),
        norb_phys=norb_phys)
    logger.info("=== Spin-orbital dynamic Eliashberg solver (issue #284) ===")
    logger.info("T = %s, beta = %.6g, grid %dx%dx%d, Nmat = %d, nd = %d "
                "spin-orbitals (%d physical orbitals), chi0q_mode = %s, "
                "solver_mode = %s", T, beta, Nx, Ny, Nz, nmat, nd, norb_phys,
                chi0q_mode, solver_mode)

    # --- memory guard, from the configuration and the archive headers ---
    nmat_est = nmat
    if chi0q_mode == "flex":
        flex_dir = sc._resolve_flex_dir(input_dict)
        for fname, keys, axis in (
                (eli.get("flex_chiq", "chiq.npz"), ("chiq",), 0),
                (eli.get("flex_chi0q", "chi0q.npz"), ("chi0q",), 0),
                (eli.get("flex_green", "green.npz"), ("green",), 1)):
            nf = _npz_freq_size(os.path.join(flex_dir, fname), keys, axis)
            if nf is not None and nf != nmat:
                raise ValueError(
                    "FLEX archive '{}' holds {} Matsubara frequencies, but "
                    "[mode.param] Nmat = {}: frequency = \"dynamic\" needs the "
                    "FLEX outputs on the run's grid. Re-run FLEX with the same "
                    "Nmat.".format(os.path.join(flex_dir, fname), nf, nmat))
    need = estimate_memory_bytes(nd, Nk, nmat_est, solver_mode, num_eigenvalues,
                                 chi0q_mode)
    check_memory(need, eli.get("mem_limit_gb"))

    use_gpu = backend.as_bool(eli.get("gpu", False))
    xp, gpu_active = backend.get_backend(
        use_gpu, logger=logger,
        required=backend.as_bool(eli.get("gpu_required", False)))

    kx = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)

    # --- inputs: chi at every bosonic frequency, and G ---
    gam = sc._so_rpa_vertex(input_dict)
    if chi0q_mode == "flex":
        chi_w, green_w, chi0_w, chi_path = sc._load_flex_so(
            input_dict, nd, Nx, Ny, Nz, interactions=interactions, full=True)
        logger.info("Using the spin-orbital FLEX chiq (all %d bosonic "
                    "frequencies) and dressed Green's function", nmat)
        consistency = _check_vertex_every_frequency(gam, chi_w, chi0_w,
                                                    chi_path, T)
        del chi0_w
    else:
        epsilon_k = sc._build_hamiltonian_k(kx, ky, kz, hr, nd)
        epsilon_k = spin_orbital.interleaved_to_spinblock(
            epsilon_k, norb_phys, axes=(0, 1))
        eigenvalues, eigenvectors = sc._calc_eigenvalues(epsilon_k)
        n_filling = mode_param.get("filling")
        if n_filling is None:
            raise ValueError("filling must be specified.")
        mu = sc._determine_mu(eigenvalues, beta, n_filling, nd)
        logger.info("mu = %.6f", mu)
        green_w = sc._build_bond_green(
            eigenvalues, eigenvectors, mu, beta, nmat, 0.0, None).full_sc
        chi0q_raw = sc._calc_chi0q_internal(input_dict, chi0q_tensor="general",
                                            precomputed_mu=mu)
        nvol, nd2 = Nk, nd * nd
        if np.shape(chi0q_raw) != (nmat, nvol) + (nd,) * 4:
            raise ValueError(
                "internal error: the spinful chi0q has shape {}, expected "
                "{}".format(np.shape(chi0q_raw), (nmat, nvol) + (nd,) * 4))
        chi_w = np.empty((nmat, nvol, nd2, nd2), dtype=complex)
        eye = np.eye(nd2, dtype=complex)
        for m in range(nmat):
            c0 = np.asarray(chi0q_raw[m], dtype=complex).reshape(nvol, nd2, nd2)
            chi_w[m] = np.linalg.solve(eye[None] + c0 @ gam, c0)
        del chi0q_raw, c0       # c0 is a view: it would keep chi0q alive
        consistency = None

    # --- vertex (q, i nu) -> (r, tau), then the pair bubble ---
    t0 = time.time()
    V_w = so_pairing.so_pairing_vertex_dynamic(gam, chi_w, nd, (Nx, Ny, Nz))
    del chi_w, gam
    if gpu_active:
        backend.warn_if_device_memory_short(
            3 * V_w.nbytes, logger, label="the spin-orbital dynamic Eliashberg "
            "kernel")
        logger.info("GPU backend active (CuPy): moving the pairing vertex and G2 "
                    "to the device (%.2f GB each).", V_w.nbytes / 1e9)
        V_w = xp.asarray(V_w)
    Vs_rt = vertex_to_rt(V_w, workers=fft_workers)
    del V_w
    G2_w = calc_g2_dynamic(green_w, beta)
    del green_w
    if gpu_active:
        G2_w = xp.asarray(G2_w)
    A, basis = make_operator(Vs_rt, G2_w, beta, workers=fft_workers)
    logger.info("vertex, pair bubble and basis built in %.2f s; antisymmetric "
                "space dimension %d of %d (self-paired %d)", time.time() - t0,
                basis.dim, basis.n, basis.n_self)
    if basis.n_self != 0:
        raise ValueError("internal error: the dynamic Pauli map has {} "
                         "self-paired elements".format(basis.n_self))

    leakage = pauli_leakage(A, basis)
    logger.info("Pauli leakage of the kernel (antisymmetric probes): %.3e", leakage)
    if not leakage <= PAULI_LEAKAGE_WARN:
        logger.warning(
            "the spin-orbital dynamic kernel maps antisymmetric gaps out of the "
            "antisymmetric space by %.3e (> %.0e); the solve is restricted to the "
            "antisymmetric space, but the kernel does not commute with the Pauli "
            "map to round-off.", leakage, PAULI_LEAKAGE_WARN)

    gap_shape = basis.shape
    seed_path = eli.get("seed_eigenvector")
    seed_full = load_seed(seed_path, gap_shape) if seed_path else None

    # --- solve ---
    eigenvalues_eig = eig_weights = match = None
    note = None
    t0 = time.time()
    if solver_mode == "iteration":
        sigma_init = initial_gap(init_gap_mode, norb_phys, kx, ky, kz, gap_shape)
        info = {}
        gap, lam, converged, n_iter = sc._so_solve_iteration(
            A, basis, sigma_init, int(eli.get("max_iter", 1000)),
            eli.get("alpha", 0.5), eli.get("convergence_tol", 1.0e-5),
            spectral_shift=spectral_shift, info_out=info)
        logger.info("Iteration result: eigenvalue = %.6f, converged = %s, n_iter = "
                    "%d", lam, converged, n_iter)
        selection = "iteration"
        if spectral_shift is not None:
            note = sc._shifted_eigenvalue_note(
                "iteration", spectral_shift, converged, n_iter, info,
                kernel_label="dynamic kernel")
        gap_w = _fix_gauge(gap)
        if not antisymmetry_residual(gap_w) <= sc._SO_ANTISYM_TOL:
            raise ValueError(
                "[mode] enable_spin_orbital = true: the power iteration returned "
                "a gap that is not fermion-antisymmetric (Delta_AB(k, iw_n) = "
                "-Delta_BA(-k, -iw_n)) to round-off; it is not written.")
        lam = float(np.real(lam))
    else:
        if solver_mode == "both":
            logger.warning(
                "Dynamic solver_mode='both' runs the eigenvalue leg only; the "
                "power-iteration cross-check is skipped.")
        seed_red = basis.to_reduced(seed_full) if seed_full is not None else None
        vals, vecs = sc._so_solve_eigenvalue(
            A, basis, num_eigenvalues, eigenvalue_method,
            sigma_shift=eli.get("sigma_shift"), spectral_shift=spectral_shift,
            seed_vec=seed_red, map_back=False)
        selection = "LR" if spectral_shift is not None else "LM"
        if seed_red is not None:
            selection += "_seed_overlap"
        # As in the non-spin-orbital dynamic solver: a plain which='LM' set
        # (no spectral_shift, no seed, not the dense full-spectrum branch)
        # whose leading eigenvalue is negative may have lost the physical
        # leading eigenvalue to large repulsive ones, which grow with Nmat;
        # re-solve ONCE with spectral_shift = "auto" (largest real part) on
        # the same reduced operator.
        if (selection == "LM" and eigenvalue_method == "arnoldi"
                and num_eigenvalues < basis.dim - 1):
            vals_arr = np.asarray(vals)
            scale = float(np.max(np.abs(vals_arr))) if vals_arr.size else 0.0
            neg_tol = 1.0e-8 * max(1.0, scale)
            leading_re = float(np.real(vals_arr[0])) if vals_arr.size else float("nan")
            if not vals_arr.size or leading_re < -neg_tol:
                logger.warning(
                    "%s: arnoldi (which='LM', num_eigenvalues=%d) returned "
                    "no positive '%s' eigenvalue (leading Re(lambda) = "
                    "%.4g); re-solving with spectral_shift = \"auto\" "
                    "(largest real part) -- set [eliashberg] spectral_shift "
                    "= \"auto\" or use solver_mode = \"iteration\" to avoid "
                    "the extra solve",
                    "Spin-orbital dynamic Eliashberg", num_eigenvalues,
                    "fermion-antisymmetric", leading_re)
                vals, vecs = sc._so_solve_eigenvalue(
                    A, basis, num_eigenvalues, eigenvalue_method,
                    sigma_shift=eli.get("sigma_shift"), spectral_shift="auto",
                    seed_vec=None, map_back=False)
                selection = "LR_retry"
                note = (
                    "eigenvalue selection: which='LM' with "
                    "num_eigenvalues={} returned no positive '{}' "
                    "eigenvalue (Re(lambda) = {:.4g}); the values below "
                    "come from a second solve with spectral_shift='auto' "
                    "(largest real part)".format(
                        num_eigenvalues, "fermion-antisymmetric", leading_re))
        rows, match = [], []
        gap_w = None
        for i in range(len(vals)):
            g = _fix_gauge(basis.to_full(vecs[:, i]).reshape(gap_shape))
            match.append(antisymmetry_residual(g) <= sc._SO_ANTISYM_TOL)
            wi = gap_weights(g, norb_phys)
            rows.append(np.concatenate([wi["spin"], wi["parity"]]))
            if i == 0:
                gap_w = g
        match = np.array(match, dtype=bool)
        sc._so_check_antisymmetric_eigenpairs(vals, match)
        eigenvalues_eig, eig_weights = np.asarray(vals), np.array(rows)
        lam = float(np.real(vals[0]))
        for i, ev in enumerate(vals):
            logger.info("  %3d: %.6f%+.6fi (|ev| = %.6f)", i, ev.real, ev.imag,
                        abs(ev))
    logger.info("Spin-orbital dynamic leading eigenvalue lambda = %.6f (solve "
                "%.2f s)", lam, time.time() - t0)

    extra = {"chi0q_mode": str(chi0q_mode)}
    if consistency is not None:
        extra["flex_vertex_backward_error_max"] = "{:.3e}".format(consistency)
    if solver_mode == "iteration":
        extra["init_gap"] = str(init_gap_mode)
    w = write_outputs(
        output_dir, gap_file, eigenvalue_file, gap_w=gap_w, lam=lam,
        eigenvalues=eigenvalues_eig, eigen_weights=eig_weights, match=match,
        norb_phys=norb_phys, T=T, kx=kx, ky=ky, kz=kz,
        eigenvalue_selection=selection, leakage=leakage, note=note,
        extra_prov=extra)
    logger.info("leading gap: spin weights %s; frequency parity %s",
                _fmt_weights(so_pairing.SPIN_COMPONENTS, w["spin"]),
                _fmt_weights(FREQUENCY_PARITY_LABELS, w["parity"]))
    logger.info("=== Done (%.2f s) ===", time.time() - t_start)
    return lam
