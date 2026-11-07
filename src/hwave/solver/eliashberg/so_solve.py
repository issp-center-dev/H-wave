"""The spin-orbital static solve of hwave_sc: state and provenance records, the static vertex, the sector solve and its report.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 4); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import NamedTuple
import numpy as np

from hwave.solver import backend
from hwave.solver.eliashberg.config import _validate_spectral_shift
from hwave.solver.eliashberg.band import _build_hamiltonian_k
from hwave.solver.eliashberg.kernel import _green_sc_to_canonical
from hwave.solver.eliashberg.eigen import _operator_is_exactly_zero, _order_by_seed_overlap, _order_eigenpairs, _reverse_k_and_orbital, _normalize_spectral_shift, _solve_leading

logger = logging.getLogger("hwave_sc")


#: pp_sector of a result that is not sector-resolved (dense layout, or none)
_SO_PP_DENSE = 2


#: pair-S_z sectors in their tie-break order
_SO_PP_ORDER = (1, 0, -1)


#: per-pair residual bound ||K v - lambda v|| <= tol * n_K * ||v||
_SO_SECTOR_RESIDUAL_TOL = 1.0e-10


#: a sector part of the seed at most this fraction of the seed (2-norm) is
#: round-off: the sector is solved unseeded
_SO_SEED_PART_RTOL = 1.0e-8


#: the seeded union warns when no selected pair overlaps the seed this much
_SO_SEED_OVERLAP_WARN = 0.5


#: eigenvalues closer than this (relative, floor 1) form one tie cluster
_SO_SECTOR_TIE_TOL = 1.0e-12


#: steps and relative-change stop of the operator-norm power estimate
_SO_NORM_STEPS = 30


_SO_NORM_RTOL = 1.0e-3


#: Bound on the relative antisymmetry residual ``||v - P v|| / ||v||`` of a
#: spin-orbital gap (:func:`_so_antisymmetry_residual`). The solvers return
#: gaps that are antisymmetric by construction (residual exactly 0, see
#: :class:`_SOAntisymmetricBasis`); anything above round-off is refused.
_SO_ANTISYM_TOL = 1.0e-10


#: Backward-error bound of the Gamma consistency check of
#: :func:`_check_so_flex_vertex` (round-off of the FLEX solve is ~1e-15; an
#: interaction that differs from the FLEX run's shows up at O(dU/U)).
_SO_FLEX_GAMMA_TOL = 1.0e-8


def _so_coeff_extern(input_dict):
    """``[mode.param] coeff_extern`` (case-insensitive, default 0.0): a
    finite real number."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    raw = CaseInsensitiveDict(input_dict["mode"].get("param", {}) or {}).get(
        "coeff_extern", 0.0)
    return spin_orbital.validate_coeff_extern(raw)


class _SORun(NamedTuple):
    """This spin-orbital hwave_sc run as the archive checks see it (#88):
    its band (:func:`hwave.solver.fingerprint.band_inputs`), its state
    (``beta``, ``nmat``, ``calc_mu = True``, ``target`` in electrons per cell;
    ``mu``, ``density`` and ``mu_tol`` are NaN, in memory only, never
    compared and never serialized) and the entry context
    (:class:`hwave.solver.spin_orbital.SOEntryContext`: the interaction
    compiled once, the ``orbital_site`` WARNING state) that the nested RPA
    runs reuse."""
    band: object
    state: object
    context: object = None


def _so_run_inputs(input_dict, hr, extern_table, norb_phys, beta, nmat,
                   n_filling, context=None):
    """The :class:`_SORun` of a spin-orbital hwave_sc run, built from its own
    input before any archive is read. hwave_sc always searches mu
    (``filling`` is mandatory), so the target is ``filling * nd`` electrons
    per cell over all ``nd = 2 * norb_phys`` spin-orbitals. The band is the
    Transfer table as read (interleaved), the physical-orbital Extern table
    (or None), and the metadata of the spin-orbital RPA/FLEX producers
    (``coeff_extern``, the index convention, the orbital counts, CellShape
    and the resolved SubShape). ``context`` is stored as
    :attr:`_SORun.context`."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import fingerprint, spin_orbital
    # case-insensitive, like every [mode.param] key of the producers
    param = CaseInsensitiveDict(input_dict["mode"].get("param", {}) or {})
    cs = param["CellShape"]
    cs = list(cs) if isinstance(cs, (list, tuple)) else [cs]
    cs = cs + [1] * (3 - len(cs))
    ss = param.get("SubShape", cs)
    ss = list(ss) if isinstance(ss, (list, tuple)) else [ss]
    ss = ss + [1] * (3 - len(ss))
    meta = fingerprint.BandMeta(
        coeff_extern=_so_coeff_extern(input_dict), enable_spin_orbital=True,
        index_convention=spin_orbital.INDEX_CONVENTION,
        norb=2 * int(norb_phys), norb_phys=int(norb_phys),
        cell_shape=tuple(int(x) for x in cs), sub_shape=tuple(int(x) for x in ss))
    nan = float("nan")
    state = fingerprint.SOState(float(beta), int(nmat), nan, nan, True,
                                float(n_filling) * 2 * int(norb_phys), nan)
    return _SORun(fingerprint.band_inputs(hr, extern_table, meta), state, context)


def _so_own_state(eigenvalues, beta, mu, nmat, target):
    """The state of a kernel built on this run's own band and mu: the masked
    Fermi count of the eigenvalues at ``mu`` (``rpa._masked_fermi_delta_n``,
    the count ``_determine_mu`` roots) per cell, the ``calc_mu = True``
    branch."""
    from hwave.solver import fingerprint
    from hwave.solver import rpa as _rpa
    w = np.asarray(backend.to_host(eigenvalues))
    nvol = int(np.prod(w.shape[:3]))
    count = _rpa._masked_fermi_delta_n(w, 1.0 / beta, mu, 0.0, 100.0)[0]
    return fingerprint.so_state(beta=beta, nmat=nmat, mu=mu,
                                density=float(count) / nvol, calc_mu=True,
                                target=target)


def _so_flex_consumed_state(green_sc, epsilon_k, sigma_static, mu, beta,
                            target):
    """The state of a borrowed spin-orbital FLEX Green function as consumed:
    the closure density FLEX itself uses
    (``hartree_fock.equal_time_density`` with the reference
    ``heff_eigenpairs(H0_k, sigma_static)``, flex.py
    ``_phase_b_density_and_hf``) per cell, at the mu the kernel applied.
    ``green_sc`` is ``(nd, nd, Nx, Ny, Nz, nmat)``, ``epsilon_k`` this run's
    spin-block band ``(nd, nd, Nx, Ny, Nz)`` and ``sigma_static`` the FLEX
    run's static self-energy ``(nvol, nd, nd)``."""
    from hwave.solver import fingerprint, hartree_fock as _hf
    eps = np.asarray(backend.to_host(epsilon_k))
    nd = eps.shape[0]
    shape = tuple(int(x) for x in eps.shape[2:5])
    nvol = int(np.prod(shape))
    H0 = np.moveaxis(eps.reshape(nd, nd, nvol), -1, 0)
    g = _green_sc_to_canonical(np.asarray(backend.to_host(green_sc)))
    dens = _hf.equal_time_density(g, _hf.heff_eigenpairs(H0, sigma_static), mu,
                                  beta, shape, sym_tol=None)
    return fingerprint.so_state(beta=beta, nmat=int(g.shape[1]), mu=mu,
                                density=float(np.real(dens.n_per_spin)) / nvol,
                                calc_mu=True, target=target)


def _so_writer_members(stage, label, so_run, fp_green, state):
    """The fingerprint members of a hwave_sc gap archive: this run's band
    and interaction fingerprints, and the Green function and state its
    kernel used (never the producer's); ``fp_green`` is
    :func:`_so_green_fp` of that Green function."""
    from hwave.solver import fingerprint
    return fingerprint.so_provenance_members(
        producer_stage=stage, producer_label=label,
        transfer=so_run.band.transfer, extern=so_run.band.extern,
        meta=so_run.band.meta, fp_green=fp_green, state=state,
        interaction=so_run.context.interaction)


def _so_epsilon_k(input_dict, kx_array, ky_array, kz_array, hr, extern_table,
                  norb_phys):
    """The spin-orbital one-body band ``(nd, nd, Nx, Ny, Nz)`` in spin-block
    order (``nd = 2 * norb_phys``), as the spin-orbital RPA assembles it::

        interleaved_to_spinblock(T(k)) + coeff_extern * H1(k) (x) sigma^z

    ``T(k)`` from the interleaved Transfer table ``hr``, ``H1(k)`` from the
    physical-orbital Extern table (``None``: no field). The Zeeman sign is
    that of the Extern definition: ``+h`` shifts spin up by ``+h`` and spin
    down by ``-h``. A Zeeman field written into the Transfer file as well
    (Route 1) adds to it, with a WARNING."""
    from hwave.solver import spin_orbital
    nd = 2 * norb_phys
    epsilon_k = spin_orbital.interleaved_to_spinblock(
        _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, nd),
        norb_phys, axes=(0, 1))
    if extern_table is None:
        return epsilon_k
    coeff = _so_coeff_extern(input_dict)
    if coeff == 0.0:
        return epsilon_k
    up, dn = 0, 1
    for a in range(norb_phys):
        t_up = hr.get(((0, 0, 0), (spin_orbital.interleaved_index(a, up),
                                   spin_orbital.interleaved_index(a, up))), 0.0)
        t_dn = hr.get(((0, 0, 0), (spin_orbital.interleaved_index(a, dn),
                                   spin_orbital.interleaved_index(a, dn))), 0.0)
        if t_up != t_dn:
            logger.warning(
                "Extern adds to spin-dependent on-site terms already present "
                "in the Transfer file (orbital {}); the two fields add".format(
                    a + 1))
            break
    h1 = coeff * _build_hamiltonian_k(kx_array, ky_array, kz_array,
                                      extern_table, norb_phys)
    if np.any(h1 != 0):
        logger.info("Zeeman sign: +h raises spin-up (favours spin-down "
                    "occupation)")
    epsilon_k[:norb_phys, :norb_phys] += h1
    epsilon_k[norb_phys:, norb_phys:] -= h1
    return epsilon_k


def _so_rpa_vertex(input_dict, so_context=None):
    """The antisymmetrized spinful vertex ``Gamma`` ``(nvol, nd^2, nd^2)``
    (:func:`hwave.solver.rpa.spinful_bubble_vertex`) of a spin-orbital
    config, built from the interaction objects of the spin-orbital RPA
    (``ham_inter_q`` and ``spinful_exchange()``, the tensors the RPA and FLEX
    SO resum; interaction files in physical orbital indices). With
    ``so_context`` (the entry context of the hwave_sc run) the RPA reuses it:
    no second compile of the interaction, no second ``orbital_site``
    WARNING."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_inputfile = dict(input_dict.get("file", {}).get("input", {}))
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_mode_rpa = dict(input_dict.get("mode", {}))
    info_mode_rpa["calc_scheme"] = "general"
    info_mode_rpa["_orbital_site_reported"] = True
    info_mode_rpa["enable_spin_orbital"] = True
    if so_context is not None:
        from hwave.solver import spin_orbital
        spin_orbital.context_for_nested(info_mode_rpa, so_context)
    read_io = read_input_k.QLMSkInput(info_inputfile)
    solver = sol_rpa.RPA(read_io.get_param("ham"),
                         {"print_level": 1, "print_step": 1}, info_mode_rpa)
    nvol, nd = solver.lattice.nvol, 2 * solver.norb
    ham = np.asarray(backend.to_host(solver.ham_info.ham_inter_q)).reshape(
        nvol, nd, nd, nd, nd)
    return sol_rpa.spinful_bubble_vertex(ham, solver.ham_info.spinful_exchange())


def _so_gap_npz_name(gap_file):
    """The spin-orbital gap archive: the gap file name with its extension
    replaced by ``.npz`` (``gap.dat`` -> ``gap.npz``); a gap file that is
    itself named ``*.npz`` gets ``.npz`` appended instead, so the two never
    collide."""
    name = os.path.splitext(gap_file)[0] + ".npz"
    return gap_file + ".npz" if name == gap_file else name


def _so_flex_writer_members(stage, label, so_run, borrowed, green_sc, epsilon_k,
                            beta):
    """The writer members of a run whose kernel used a borrowed spin-orbital
    FLEX Green function: its digest as consumed (the digest of the archive
    check), the mu the kernel applied (the validated archive value) and the
    closure density of that G (:func:`_so_flex_consumed_state`)."""
    state = _so_flex_consumed_state(
        green_sc, epsilon_k, borrowed.sigma_static, borrowed.info.producer_state.mu,
        beta, so_run.state.target)
    return _so_writer_members(stage, label, so_run, borrowed.fp_green, state)


def _so_prov_line(v):
    """The ``# key = value`` text of a scalar provenance member."""
    a = np.asarray(v)
    if a.dtype.kind == "b":
        return "true" if bool(a) else "false"
    if a.dtype.kind == "f":
        return repr(float(a))
    if a.dtype.kind in "iu":
        return str(int(a))
    return str(a[()])


def _so_provenance_lines(members):
    """The scalar members of ``members`` as ``{key: text}`` (the arrays stay
    in the npz only)."""
    return {k: _so_prov_line(v) for k, v in (members or {}).items()
            if np.ndim(v) == 0}


def _so_output_record(output_dir, gap_file, norb_phys, gap, eigenvalues,
                      eigenvectors, match, kx_array, ky_array, kz_array,
                      sz_residual=None, provenance_members=None,
                      eigen_converged=None):
    """The spin-orbital output schema (issue #83).

    ``eigen_converged`` (bool, optional) is stored as the archive key of
    that name: False when the gap is the last iterate of a power iteration
    that did not converge (not an eigenvector).

    Writes the gap archive (:func:`_so_gap_npz_name`, default ``gap.npz``)
    and returns ``(provenance, weights, gap_comments)`` for
    :func:`_save_results`:

    * ``provenance`` -- the ``# key = value`` lines of ``eigenvalue.dat``
      (``spin_orbital = true``, ``index_convention = spin_block``, the gap
      layout, the spin decomposition and the weights of the written gap);
    * ``weights`` -- ``(n_eig, 4)`` spin weights ``(psi, dx, dy, dz)`` of the
      eigenvectors (:func:`hwave.solver.so_pairing.so_spin_decomposition`),
      or None without an eigenvalue analysis;
    * ``gap_comments`` -- the comment lines of ``gap.dat``.

    The archive holds ``gap`` ``(nd, nd, Nx, Ny, Nz)`` (the gap of
    ``gap.dat``), ``index_convention``, ``spin_orbital``, ``norb_phys``,
    ``kx``/``ky``/``kz``, ``spin_components`` and ``gap_spin_weights``; with
    an eigenvalue analysis also ``eigenvalues``, ``eigenvectors``
    ``(n_eig, nd, nd, Nx, Ny, Nz)``, ``eigenvector_spin_weights`` and
    ``eigenvalue_match``.

    With ``provenance_members`` (#88: the fingerprint members of
    :func:`hwave.solver.fingerprint.so_provenance_members` for this run's band
    and the Green function and state its kernel used) the archive gets every
    member and ``provenance`` gets the scalar ones.

    With ``sz_residual`` (the S_z-conservation residual of the run,
    :func:`hwave.solver.so_pairing.sz_conservation_residual`) the S_z report
    of :func:`hwave.solver.so_pairing.sz_report` is logged and its members
    (``sz_weight_upup``, ``sz_weight_dndn``, ``sz_weight_singlet0``,
    ``sz_weight_triplet0`` ``(n_ev,)``, ``sz_conservation_residual``,
    ``sz_conserved``) are added to the archive; ``n_ev`` counts the
    eigenvectors, or is 1 (the written gap) without an eigenvalue analysis."""
    from hwave.solver import so_pairing, spin_orbital
    comps = so_pairing.SPIN_COMPONENTS
    nd = 2 * int(norb_phys)
    layout = ("Delta_AB(k), A = spin * {} + orb (spin 0 = up, spin 1 = down); "
              "fermion-antisymmetric, Delta_AB(k) = -Delta_BA(-k)"
              .format(norb_phys))
    npz_name = _so_gap_npz_name(gap_file)
    prov = {
        "spin_orbital": "true",
        "index_convention": spin_orbital.INDEX_CONVENTION,
        "gap_layout": layout,
        "gap_npz": npz_name,
        "match_column": "1 = fermion-antisymmetric eigenvector",
        "spin_decomposition": (
            "Delta_ab(k) = (psi_ab(k) + d_ab(k) . sigma) i sigma_y per orbital "
            "pair (a, b); w_X = sum_{k,a,b} |X_ab(k)|^2 / sum_{k,a,b} "
            "(|psi_ab(k)|^2 + |d_ab(k)|^2), X = psi, dx, dy, dz"),
        "weight_columns": "w_psi w_dx w_dy w_dz (after the match column)",
    }
    record = {"index_convention": np.array(spin_orbital.INDEX_CONVENTION),
              "spin_orbital": np.array(True),
              "norb_phys": np.array(int(norb_phys)),
              "kx": np.asarray(kx_array), "ky": np.asarray(ky_array),
              "kz": np.asarray(kz_array),
              "spin_components": np.array(comps)}
    if gap is not None:
        gap = np.asarray(backend.to_host(gap))
        w = so_pairing.so_spin_decomposition(gap, norb_phys)
        prov["gap_spin_weights"] = " ".join(
            "{}={:.8e}".format(c, w[c]) for c in comps)
        record["gap"] = gap.reshape((nd, nd) + gap.shape[2:])
        record["gap_spin_weights"] = np.array([w[c] for c in comps])
    if eigen_converged is not None:
        record["eigen_converged"] = np.array(bool(eigen_converged))
    weights = None
    if eigenvalues is not None:
        vecs = np.stack([np.asarray(backend.to_host(v)) for v in eigenvectors])
        weights = np.array([[so_pairing.so_spin_decomposition(v, norb_phys)[c]
                             for c in comps] for v in vecs])
        record["eigenvalues"] = np.asarray(eigenvalues)
        record["eigenvectors"] = vecs
        record["eigenvector_spin_weights"] = weights
        if match is not None:
            record["eigenvalue_match"] = np.asarray(match, dtype=bool)
    if sz_residual is not None:
        if eigenvalues is not None:
            sz_gaps = vecs
        elif gap is not None:
            sz_gaps = [gap]
        else:
            sz_gaps = []
        if len(sz_gaps):
            record.update(so_pairing.sz_report(
                [so_pairing.so_sz_decomposition(v, norb_phys) for v in sz_gaps],
                sz_residual, logger))
    if provenance_members:
        # #88: the fingerprint members of this run (its band, the G and state its
        # kernel used): the npz holds all of them, eigenvalue.dat the scalars
        record.update(provenance_members)
        prov.update(_so_provenance_lines(provenance_members))
    gap_comments = [
        "spin_orbital = true",
        "index_convention = {}".format(spin_orbital.INDEX_CONVENTION),
        "gap_layout = {}; sigma_ij below is Delta_ij".format(layout),
    ]
    if gap is not None:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, npz_name)
        logger.info("Saving the spin-orbital gap to {}".format(path))
        np.savez(path, **record)
    return prov, weights, gap_comments


def _so_static_vertex_from_chi(gam, chi_static, nd, spatial_shape):
    """The static spin-orbital pairing vertex ``Vs_q``
    ``(nd, nd, nd, nd, Nx, Ny, Nz)`` from the spinful vertex ``gam`` and a
    static spinful susceptibility ``chi_static`` (both ``(nvol, nd^2,
    nd^2)``, bubble ``[a, c, b, d]`` layout, spin-block order):
    :func:`hwave.solver.so_pairing.so_pairing_vertex`. Pure; the susceptibility
    may come from the RPA (:func:`_so_static_vertex`) or from FLEX."""
    from hwave.solver.so_pairing import so_pairing_vertex
    return so_pairing_vertex(gam, chi_static, nd, spatial_shape)


def _so_static_chi0(chi0q_raw, nd, Nx, Ny, Nz, static_index):
    """The static spinful bubble ``(nvol, nd^2, nd^2)`` (bubble ``[a, c, b,
    d]`` layout) of a general spinful chi0q -- raw ``(nmat, nvol, nd, nd, nd,
    nd)`` or reference ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)`` -- at the zero
    bosonic frequency ``static_index`` (``nmat // 2`` when ``None``). Any other
    layout (reduced 2-index, spin-free / spin-diagonal) is refused."""
    nvol = Nx * Ny * Nz
    chi0q_raw = np.asarray(chi0q_raw)
    if chi0q_raw.ndim == 6 and chi0q_raw.shape[1:] == (nvol,) + (nd,) * 4:
        nmat = chi0q_raw.shape[0]
        si = nmat // 2 if static_index is None else static_index
        chi0 = chi0q_raw[si].reshape(nvol, nd * nd, nd * nd)
    elif chi0q_raw.ndim == 8 and chi0q_raw.shape[:7] == (nd,) * 4 + (Nx, Ny, Nz):
        nmat = chi0q_raw.shape[-1]
        si = nmat // 2 if static_index is None else static_index
        chi0 = chi0q_raw[..., si].reshape(nd * nd, nd * nd, nvol).transpose(2, 0, 1)
    else:
        raise ValueError(
            "[mode] enable_spin_orbital = true: the chi0q has shape {}, not the "
            "general (4-index) spinful layout over nd = {} spin-orbitals "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}). A reduced (2-index) or a "
            "spin-free / spin-diagonal chi0q cannot carry the spin-flip "
            "particle-hole pairs; use chi0q_mode = \"calc\", or load the "
            "chi0q of a spin-orbital RPA run with calc_scheme = \"general\" and a "
            "spin-mixing transfer.".format(chi0q_raw.shape, nd, nvol, nd=nd))
    return chi0


def _so_static_vertex(input_dict, chi0q_raw, norb_phys, Nx, Ny, Nz,
                      static_index, so_context=None):
    """The static spin-orbital pairing vertex of the RPA paths
    (``chi0q_mode = "calc"`` / ``"load"``): ``Gamma`` from the spin-orbital RPA
    objects (:func:`_so_rpa_vertex`), ``chi = [1 + chi0 Gamma]^-1 chi0`` at the
    zero bosonic frequency ``static_index`` (``nmat // 2`` when ``None``) of
    the general spinful ``chi0q_raw`` -- raw ``(nmat, nvol, nd, nd, nd, nd)``
    or reference ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)`` -- then
    :func:`_so_static_vertex_from_chi`. ``so_context`` is passed to
    :func:`_so_rpa_vertex`."""
    nd = 2 * int(norb_phys)
    chi0 = _so_static_chi0(chi0q_raw, nd, Nx, Ny, Nz, static_index)
    chi0 = np.ascontiguousarray(chi0, dtype=complex)
    gam = _so_rpa_vertex(input_dict, so_context=so_context)
    eye = np.eye(nd * nd, dtype=complex)
    chi = np.linalg.solve(eye[None] + chi0 @ gam, chi0)
    return _so_static_vertex_from_chi(gam, chi, nd, (Nx, Ny, Nz))


def _so_flex_backward_error(gam, chi, chi0):
    """The backward error of :func:`_check_so_flex_vertex`: the maximum over
    the q blocks of ``|(1 + chi0 gam) chi - chi0| / (|chi0| + |chi0 gam|
    |chi|)`` (Frobenius per block; a zero-scale block counts 0 if its residual
    vanishes and infinity otherwise; NaN propagates). ``chi`` and ``chi0`` are
    ``(..., nd^2, nd^2)`` with any leading axes."""
    c0g = chi0 @ gam
    res = np.linalg.norm(chi + c0g @ chi - chi0, axis=(-2, -1))
    scale = (np.linalg.norm(chi0, axis=(-2, -1))
             + np.linalg.norm(c0g, axis=(-2, -1))
             * np.linalg.norm(chi, axis=(-2, -1)))
    pos = scale > 0
    per_q = np.where(res == 0, 0.0, np.inf)
    per_q = np.where(pos, res / np.where(pos, scale, 1.0), per_q)
    # NaN in any block propagates (and is refused by the callers)
    per_q = np.where(np.isnan(res) | np.isnan(scale), np.nan, per_q)
    return float(np.max(per_q)) if per_q.size else 0.0


def _check_so_flex_vertex(gam, chi, chi0, chi_path):
    """Refuse when ``gam`` (built here from the interaction files) is not the
    vertex the FLEX run resummed: the stored ``chi`` must solve
    ``(1 + chi0 gam) chi = chi0`` with the stored ``chi0`` of the same map.
    Measured as a backward error, ``|(1 + chi0 gam) chi - chi0| /
    (|chi0| + |chi0 gam| |chi|)`` (Frobenius) PER q block, and the maximum
    over q is compared with the tolerance: the bound does not depend on how
    close the run is to an instability, nor on the number of q points (a
    global norm would dilute one bad block by ``1/sqrt(nvol)``). A block of
    zero scale (``chi0_q = 0`` and ``chi0_q gam chi_q = 0``) has error 0 if
    its residual vanishes and infinity otherwise."""
    err = _so_flex_backward_error(gam, chi, chi0)
    logger.info("spin-orbital FLEX: vertex consistency backward error %.3e", err)
    if not err <= _SO_FLEX_GAMMA_TOL:
        raise ValueError(
            "the interaction files of this run do not reproduce the vertex of "
            "the FLEX run that wrote '{}' (backward error {:.3e} > {:.1e} of "
            "(1 + chi0 Gamma) chi = chi0 at the zero bosonic frequency). Use the "
            "same interaction files and [mode.param] as the FLEX run.".format(
                chi_path, err, _SO_FLEX_GAMMA_TOL))
    return err


# --- the per-sector solvers of the S_z-packed spin-orbital kernel -------------------
#
# With S_z conservation the pairing kernel is block-diagonal in the pair-S_z
# sector p = S_z(B) + S_z(C) of a gap component Delta_BC (p in +1, 0, -1). Each
# sector is a SectorOp: an operator on the sector vector and the subset basis
# of its antisymmetric coordinates. The solvers below work sector by sector and
# merge the results; the dense layout keeps _so_solve_eigenvalue /
# _so_solve_iteration, and both layouts meet in the one report
# (_so_bond_report).


class SectorOp(NamedTuple):
    """One pair-S_z sector of the spin-orbital pairing kernel.

    ``operator`` acts on the sector vector of ``basis`` (a
    :meth:`_SOAntisymmetricBasis.from_subset` basis); ``has_adjoint = True``
    promises a working ``operator.rmatvec``, and a solver then reduces it
    with ``basis.reduce(operator, adjoint=True)``."""
    pp_sector: int
    operator: object
    basis: object
    has_adjoint: bool


def _so_antisymmetry_residual(gap):
    """``||v - P v|| / ||v||`` of a spin-orbital gap ``(nd, nd, Nx, Ny, Nz)``
    with ``P`` the projection onto fermion-antisymmetric gaps; ``inf`` for a
    zero or non-finite gap. Unlike :func:`_is_gap_parity` (a dominant-sector
    test used by the non-spin-orbital channels) this is a round-off-scale
    criterion: compare with :data:`_SO_ANTISYM_TOL`."""
    gap = np.asarray(gap)
    n = float(np.linalg.norm(gap))
    if n == 0.0 or not np.isfinite(n):
        return float("inf")
    sym = 0.5 * (gap + _reverse_k_and_orbital(gap))
    return float(np.linalg.norm(sym)) / n


def _so_is_antisymmetric(gap):
    return _so_antisymmetry_residual(gap) <= _SO_ANTISYM_TOL


def _so_solve_iteration(A, basis, sigma_init, max_iter, alpha, tol,
                        spectral_shift=None, info_out=None):
    """The power iteration of the spin-orbital kernel on the antisymmetric
    coordinates of ``basis`` (the shared ``_solve_leading`` driver on
    ``E^dagger K E``); returns ``(gap, eigenvalue, converged, n_iter)`` with
    the gap mapped back by ``E``."""
    red = basis.reduce(A)
    c0 = basis.to_reduced(np.asarray(sigma_init).ravel())
    eigenvalue, c, info = _solve_leading(
        lambda: (red, basis.dim), basis.dim, "iteration",
        max_iter=max_iter, convergence_tol=tol, init_vec=c0, alpha=alpha,
        project_fn=None, spectral_shift=spectral_shift)
    if info_out is not None:
        info_out.update(info)
    gap = basis.to_full(c).reshape(basis.shape)
    return gap, eigenvalue, info["converged"], info["n_iter"]


def _so_dense_reduced_eigen(red, dim, num_eigenvalues, seed_vec=None):
    """The dense diagonalization of the reduced operator ``red`` (``dim x
    dim``): ``min(max(1, num_eigenvalues), dim)`` eigenpairs selected by the
    largest real part, ordered by descending real part (by the overlap with
    ``seed_vec`` when given). Returns ``(vals, vecs, M)`` with the
    eigenvectors as columns and ``M`` the dense reduced matrix. Logs
    nothing (the dense branch of :func:`_so_solve_eigenvalue`)."""
    dense = np.asarray(red.matmat(np.eye(dim, dtype=complex)))
    vals, vecs = np.linalg.eig(dense)
    n_keep = min(max(1, int(num_eigenvalues)), dim)
    # the dense spectrum is complete, so select by the largest REAL part
    # (the physical leading eigenvalue) whether or not a spectral_shift
    # was given; a largest-magnitude cut could drop it behind repulsive
    # modes
    sel = np.argsort(-vals.real, kind="stable")[:n_keep]
    vals, vecs = _order_eigenpairs(vals[sel], vecs[:, sel])
    if seed_vec is not None:
        vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
    return vals, vecs, dense


def _so_solve_eigenvalue(A, basis, num_eigenvalues, method, sigma_shift=None,
                         spectral_shift=None, seed_vec=None, map_back=True,
                         info_out=None):
    """The eigenvalue analysis of the spin-orbital kernel on the
    antisymmetric coordinates of ``basis``; returns ``(vals, gaps)`` with
    the eigenvectors mapped back by ``E``, shape ``(n_ev,) + shape``,
    ordered by descending real part.

    ``seed_vec`` (reduced coordinates, optional): an eigenvector-continuation
    seed, used as the ARPACK start vector and to order the eigenpairs by
    their overlap with it (``_solve_leading``; on the dense branch the
    selected eigenpairs are ordered the same way). ``map_back = False``
    returns the reduced eigenvectors as columns ``(dim, n_ev)`` instead of
    the gaps (the dynamic solver maps them back one at a time: each gap is
    ``nd^2 Nk nmat`` long). The defaults are the static solver's call.

    ARPACK needs ``k < dim - 1``; when ``num_eigenvalues >= dim - 1`` the
    reduced operator (then small) is diagonalized densely instead and
    ``min(num_eigenvalues, dim)`` eigenpairs are returned, selected by the
    largest real part (the dense spectrum is complete), so no physical
    eigenpair is dropped at the full-spectrum boundary. Only ``method = "arnoldi"`` is supported.

    ``info_out`` (dict, optional) receives ``selection``: the criterion that
    produced the eigenpairs (``"dense-LR"`` on the dense branch, else
    ``_solve_leading``'s ``"LM"`` / ``"LR"`` / ``"LR_retry"``)."""
    if method != "arnoldi":
        raise ValueError(
            "[mode] enable_spin_orbital = true supports only "
            "eigenvalue_method = \"arnoldi\", got {!r}.".format(method))
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, method)
    red = basis.reduce(A)
    dim = basis.dim
    if num_eigenvalues >= dim - 1:
        if sigma_shift is not None:
            # the same advisory the ARPACK path (_solve_leading) gives
            logger.warning(
                "[eliashberg] sigma_shift is ignored for eigenvalue_method="
                "'arnoldi' (it targets the shift-invert methods); use "
                "spectral_shift to bias the arnoldi selection.")
        vals, vecs, _ = _so_dense_reduced_eigen(red, dim, num_eigenvalues,
                                                seed_vec)
        n_keep = len(vals)
        logger.info("spin-orbital: num_eigenvalues = %d >= dim - 1 = %d of the "
                    "antisymmetric space; dense diagonalization, %d eigenpairs.",
                    num_eigenvalues, dim - 1, n_keep)
        selection = "dense-LR"
    else:
        _, _, eig = _solve_leading(
            lambda: (red, dim), dim, method,
            num_eigenvalues=num_eigenvalues, sigma_shift=sigma_shift,
            spectral_shift=spectral_shift, seed_vec=seed_vec)
        vals, vecs = eig["eigenvalues"], eig["eigenvectors"]
        selection = eig.get("selection")
    if info_out is not None:
        info_out["selection"] = selection
    if not map_back:
        return vals, vecs
    gaps = np.array([basis.to_full(vecs[:, i]).reshape(basis.shape)
                     for i in range(len(vals))])
    return vals, gaps


def _pp_label(p):
    return "+1" if p == 1 else str(int(p))


class SectorLeader(NamedTuple):
    """The eigen solve of one pair-S_z sector.

    ``value``: the largest real part among the returned eigenvalues (``0``
    for a null sector, NaN when the sector is empty or did not converge).
    ``selection``: ``"LM"``, ``"LR_retry"``, ``"LR"`` (a user
    ``spectral_shift``), ``"dense-LR"``, ``"null-operator"``,
    ``"not-converged"``, ``"residual-failed"`` or ``"empty"``.
    ``residual``: the largest ``||K v - lambda v|| / (n_K ||v||)`` over the
    returned pairs (certified when ``<= 1e-10``). ``norm_estimate`` /
    ``norm_kind``: the operator-norm scale of ``n_K`` -- ``"exact"`` (the
    2-norm of the dense matrix), ``"adjoint"`` (power iteration on ``K^dagger
    K``) or ``"forward"`` (power iteration on ``K``, a lower bound).
    ``n_returned``: the number of pairs the sector solve returned."""
    pp_sector: int
    dim: int
    value: complex
    selection: str
    certified: bool
    residual: float
    norm_estimate: float
    norm_kind: str
    n_returned: int


@dataclass(frozen=True)
class SectorSolve:
    """The union of the sector eigen solves.

    ``vals (n,)``, ``gaps (n,) + gap_shape`` and ``pp_sector (n,) int8``:
    the certified pairs of the sectors that did not fail, at most ``k`` of
    them, ordered by descending real part (ties: see
    :func:`_so_solve_sectors`; by the seed overlap in a seeded solve).
    ``leaders``: one :class:`SectorLeader` per given sector, in the order
    ``(+1, 0, -1)``. ``failed_sectors``: the
    sectors whose solve did not converge or returned a pair that failed the
    residual check. ``converged``: no failed sector and at least one pair;
    ``leading`` is then ``vals[0]``, else NaN. ``selection``: the per-sector
    selections, e.g. ``"+1:dense-LR;0:LR;-1:null-operator"``."""
    vals: np.ndarray
    gaps: np.ndarray
    pp_sector: np.ndarray
    leaders: tuple
    failed_sectors: tuple
    converged: bool
    leading: complex
    selection: str


class IterationLeader(NamedTuple):
    """The power iteration of one pair-S_z sector.

    ===============  ===================  ==========================  =========
    status           value                gap                         n_iter
    ===============  ===================  ==========================  =========
    ``"empty"``      NaN                  ``None``                    0
    ``"null"``       0.0                  the first lifted unit       0
                                          vector ``e_0``
    ``"converged"``  the signed Rayleigh  the converged iterate       the steps
                     eigenvalue
    ``"failed"``     NaN                  the last iterate            the steps
    ===============  ===================  ==========================  =========

    A one-dimensional sector is solved exactly (``"converged"``, its single
    matrix element, ``n_iter = 0``; ``"failed"`` when that element is not
    finite). A run that reports convergence with a non-finite value is
    ``"failed"`` too. ``gap`` is lifted to the full gap shape."""
    pp_sector: int
    status: str
    value: object
    gap: object
    n_iter: int


@dataclass(frozen=True)
class SectorIteration:
    """The reduction of the sector power iterations.

    Without a failed sector the leader is the largest real part over the
    non-empty sectors (null sectors at 0; ties as in
    :func:`_so_solve_sectors`) and ``value``, ``gap``, ``pp_sector`` are its.
    With a failed sector ``value`` is NaN, ``converged`` false, and ``gap`` /
    ``pp_sector`` are the last iterate of the first failed sector in the
    order ``(+1, 0, -1)``. With every sector empty ``value`` is NaN, ``gap``
    ``None`` and ``pp_sector`` 2. ``n_iter`` is the largest sector step
    count; ``n_iter_per_sector`` lists them in sector order (0 for empty,
    null and exactly solved sectors). ``selection`` is ``"null-operator"``
    when the leader is a null sector, else ``"iteration"``."""
    value: object
    gap: object
    pp_sector: int
    leaders: tuple
    failed_sectors: tuple
    converged: bool
    n_iter: int
    n_iter_per_sector: tuple
    selection: str


class EigenLeg(NamedTuple):
    """The eigenvalue leg of a spin-orbital run, in either layout."""
    vals: np.ndarray
    gaps: np.ndarray
    pp_sector: np.ndarray
    converged: bool
    leading: complex
    selection: str
    failed_sectors: tuple

    @classmethod
    def from_dense(cls, vals, gaps, *, selection):
        """The result of :func:`_so_solve_eigenvalue` (the arrays are passed
        through, not copied); converged when the first eigenvalue is finite."""
        converged = bool(len(vals) > 0 and np.isfinite(vals[0]))
        leading = complex(vals[0]) if len(vals) else complex(np.nan)
        return cls(vals, gaps, np.full(len(vals), _SO_PP_DENSE, dtype=np.int8),
                   converged, leading, selection, ())

    @classmethod
    def from_sectors(cls, solve):
        return cls(solve.vals, solve.gaps, solve.pp_sector, solve.converged,
                   solve.leading, solve.selection, solve.failed_sectors)


class IterationLeg(NamedTuple):
    """The power-iteration leg of a spin-orbital run, in either layout."""
    value: object
    gap: object
    pp_sector: int
    converged: bool
    n_iter: int
    failed_sectors: tuple
    selection: str

    @classmethod
    def from_dense(cls, gap, value, converged, n_iter, info):
        """The result of :func:`_so_solve_iteration` (``info`` is its
        ``info_out`` record, accepted for the call's symmetry)."""
        return cls(value, gap, _SO_PP_DENSE, bool(converged), int(n_iter), (),
                   "iteration")

    @classmethod
    def from_sectors(cls, iteration):
        return cls(iteration.value, iteration.gap, iteration.pp_sector,
                   iteration.converged, iteration.n_iter,
                   iteration.failed_sectors, iteration.selection)


class SectorReport(NamedTuple):
    """What a spin-orbital run reports, one schema for both layouts.

    ``reported``: the reported eigenvalue (NaN when not certified).
    ``gap`` (``gap_shape``) with ``gap_status``: ``"certified"`` (the gap of
    the certified leader; exactly when ``reported`` is finite),
    ``"uncertified"`` (a real candidate -- a listed pair or the last iterate
    -- written while the result is not certified) or ``"none"`` (no candidate
    exists; a zero array). ``pp_sector``: the sector of the written gap (2
    when not sector-resolved or none). ``eigen_converged``;
    ``iteration_converged`` (only in ``"both"``); ``selection``; the listing
    ``vals``, ``gaps``, ``pp_sectors`` (empty in ``"iteration"``);
    ``failed_sectors`` (of the authoritative leg) and
    ``iteration_failed_sectors`` (the iteration leg's, whenever it ran)."""
    reported: float
    gap: np.ndarray
    gap_status: str
    pp_sector: int
    eigen_converged: bool
    iteration_converged: object
    selection: str
    vals: np.ndarray
    gaps: np.ndarray
    pp_sectors: np.ndarray
    failed_sectors: tuple
    iteration_failed_sectors: object


def _so_canonical_sector_ops(sector_ops, solver):
    """``sector_ops`` validated (each pp_sector +1, 0 or -1, none twice) and
    put in the order ``(+1, 0, -1)``."""
    seen = {}
    for op in sector_ops:
        p = op.pp_sector
        if (isinstance(p, (bool, np.bool_))
                or not isinstance(p, (int, np.integer)) or p not in _SO_PP_ORDER):
            raise ValueError("{}: pp_sector must be +1, 0 or -1, got {}"
                             .format(solver, p))
        if int(p) in seen:
            raise ValueError("{}: pp_sector {} appears twice".format(solver, p))
        seen[int(p)] = op
    return tuple(seen[p] for p in _SO_PP_ORDER if p in seen)


def _so_check_sector_op(op, gap_shape, solver):
    basis = op.basis
    if tuple(int(x) for x in basis.shape) != tuple(gap_shape):
        raise ValueError("{}: the basis of sector {} lives on a gap of shape {}, "
                         "not {}".format(solver, _pp_label(op.pp_sector),
                                         tuple(basis.shape), tuple(gap_shape)))
    if tuple(op.operator.shape) != (basis.n, basis.n):
        raise ValueError("{}: the operator of sector {} has shape {}, its basis "
                         "acts on vectors of length {}".format(
                             solver, _pp_label(op.pp_sector),
                             tuple(op.operator.shape), basis.n))


def _so_lifted_unit_gaps(basis, count, gap_shape):
    """The first ``count`` reduced unit vectors lifted to full gaps of
    ``gap_shape``, each normalized to unit 2-norm."""
    out = np.zeros((count,) + tuple(gap_shape), dtype=complex)
    for i in range(count):
        e = np.zeros(basis.dim, dtype=complex)
        e[i] = 1.0
        g = basis.scatter(basis.to_full(e))
        g /= np.linalg.norm(g)
        out[i] = g.reshape(gap_shape)
    return out


def _so_sector_norm_estimate(red, dim, adjoint):
    """Fixed-seed power estimate of ``||red||_2``.

    ``adjoint``: ``x_0`` from ``default_rng(0)``, normalized; ``y = K x_j``,
    ``z = K^dagger y``, ``mu_j = ||z||``, ``x_{j+1} = z / ||z||``; stop after
    30 steps or when ``|mu_j - mu_{j-1}| < 1e-3 mu_j``; the estimate is
    ``sqrt(mu_last)``. Without an adjoint: ``r_j = ||K x_j||`` over the
    normalized steps ``x_{j+1} = K x_j / ||K x_j||`` with the same stop, the
    estimate ``max_j r_j`` (a lower bound). A zero image at step 0 redraws
    ``x_0`` once; a second zero gives 0. Returns ``(estimate, steps,
    relative change, converged)``."""
    rng = np.random.default_rng(0)

    def draw():
        x = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        return x / np.linalg.norm(x)

    def image(x):
        y = np.asarray(red.matvec(x)).ravel()
        if adjoint:
            y = np.asarray(red.rmatvec(y)).ravel()
        return y

    x = draw()
    z = image(x)
    if not np.any(z):
        x = draw()
        z = image(x)
        if not np.any(z):
            return 0.0, 1, 0.0, True
    mu_prev = None
    best = 0.0
    change = float("inf")
    steps = 0
    converged = False
    for j in range(_SO_NORM_STEPS):
        if j > 0:
            z = image(x)
        mu = float(np.linalg.norm(z))
        steps = j + 1
        best = max(best, mu)
        if mu == 0.0:
            break
        if mu_prev is not None:
            change = abs(mu - mu_prev) / mu
            if change < _SO_NORM_RTOL:
                converged = True
                break
        mu_prev = mu
        x = z / mu
    if adjoint:
        est = float(np.sqrt(mu))
    else:
        est = float(best)
    return est, steps, change, converged


def _so_union(cands, num_eigenvalues):
    """Order the candidate pairs ``(value, rank, in-sector index, gap,
    pp_sector)``: tie clusters by single linkage (``|l_i - l_j| <= 1e-12
    max(1, |l_i|, |l_j|)``) ordered by descending largest real part, then
    descending largest imaginary part, then the smallest sector rank of a
    member; inside a cluster by sector rank, then in-sector index. Returns
    the first ``num_eigenvalues``."""
    n = len(cands)
    parent = list(range(n))

    def root(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    vals = [complex(c[0]) for c in cands]
    for i in range(n):
        for j in range(i + 1, n):
            tol = _SO_SECTOR_TIE_TOL * max(1.0, abs(vals[i]), abs(vals[j]))
            if abs(vals[i] - vals[j]) <= tol:
                parent[root(i)] = root(j)
    clusters = {}
    for i in range(n):
        clusters.setdefault(root(i), []).append(i)
    groups = sorted(
        clusters.values(),
        key=lambda m: (-max(vals[i].real for i in m),
                       -max(vals[i].imag for i in m),
                       min(cands[i][1] for i in m)))
    order = [i for m in groups
             for i in sorted(m, key=lambda i: (cands[i][1], cands[i][2]))]
    return [cands[i] for i in order[:num_eigenvalues]]


def _so_order_by_seed(kept, seed_unit):
    """Order the kept union pairs ``(value, rank, index, gap, pp_sector)``
    by descending normalized overlap ``|<seed, gap>| / ||gap||`` with the
    unit full seed, then descending real part, then sector rank (the
    sector form of :func:`_order_by_seed_overlap`). Logs a WARNING when
    the best overlap is below 0.5: the seed's branch is then not among the
    selected pairs."""
    ov = [abs(np.vdot(seed_unit, np.ravel(c[3])))
          / (float(np.linalg.norm(c[3])) + 1.0e-300) for c in kept]
    order = sorted(range(len(kept)),
                   key=lambda i: (-ov[i], -complex(kept[i][0]).real, kept[i][1]))
    best = ov[order[0]]
    if best < _SO_SEED_OVERLAP_WARN:
        logger.warning("spin-orbital: the largest overlap of a selected eigenvector with "
                       "the seed is %.3f (< %.1f); the seed's branch is not among the "
                       "selected eigenpairs.", best, _SO_SEED_OVERLAP_WARN)
    return [kept[i] for i in order]


def _so_solve_sectors(sector_ops, num_eigenvalues, method, *, gap_shape,
                      spectral_shift=None, seed_vec=None):
    """The eigen solve of the spin-orbital kernel sector by sector.

    ``sector_ops``: :class:`SectorOp` per pair-S_z sector (any subset of
    ``+1, 0, -1``, each at most once; taken in that order). Per sector, with
    ``k = num_eigenvalues`` and ``dim`` its antisymmetric dimension:

    * ``dim = 0``: ``"empty"``, no pairs, certified;
    * an exactly zero reduced operator (:func:`_operator_is_exactly_zero`):
      ``"null-operator"``, ``min(k, dim)`` eigenvalues 0 with the first
      reduced unit vectors as gaps (lifted, unit 2-norm), certified;
    * ``dim <= max(2k + 1, 20)``: dense diagonalization (``"dense-LR"``);
    * else :func:`_so_solve_eigenvalue` on the sector (ARPACK, ``"LM"`` /
      ``"LR_retry"`` / ``"LR"``).

    The seed ``seed_vec`` (a full gap, optional) is split per sector into
    reduced coordinates; a part with ``||seed_p||_2 <= 1e-8 ||seed||_2``
    means no seed for that sector (an INFO line when it is not exactly
    zero). With a nonzero seed the kept pairs are then ordered by their
    normalized overlap with the seed, then real part, then sector rank, as
    the dense layout orders a seeded solve (:func:`_so_order_by_seed`); the
    selection of the ``k`` pairs is unchanged. Every returned pair must
    satisfy ``||K v - lambda v|| <= 1e-10 n_K ||v||`` with ``n_K = max(1,
    |lambda|, est)``, ``est`` the operator-norm estimate
    (:func:`_so_sector_norm_estimate`; exact on the dense branch), else the
    whole sector is ``"residual-failed"``. A ``"not-converged"`` or
    ``"residual-failed"`` sector is failed: its pairs are not listed and the
    result is not converged (``leading = NaN``). The certified pairs are
    merged by :func:`_so_union` and the first ``k`` kept, their gaps placed
    in full gaps of ``gap_shape``. Only ``method = "arnoldi"`` is
    supported."""
    solver = "_so_solve_sectors"
    if method != "arnoldi":
        raise ValueError(
            "{}: the per-sector solve supports eigenvalue_method = \"arnoldi\" "
            "only, got {}".format(solver, method))
    ops = _so_canonical_sector_ops(sector_ops, solver)
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, method)
    if (isinstance(num_eigenvalues, (bool, np.bool_))
            or not isinstance(num_eigenvalues, (int, np.integer))
            or num_eigenvalues < 1):
        raise ValueError("{}: num_eigenvalues must be a positive integer, got {}"
                         .format(solver, num_eigenvalues))
    k = int(num_eigenvalues)
    gap_shape = tuple(int(x) for x in gap_shape)
    for op in ops:
        _so_check_sector_op(op, gap_shape, solver)
    seed_flat = None if seed_vec is None else np.asarray(seed_vec).ravel()
    if seed_flat is not None and seed_flat.size != int(np.prod(gap_shape)):
        raise ValueError("{}: seed_vec has {} entries, the gap {} has {}".format(
            solver, seed_flat.size, gap_shape, int(np.prod(gap_shape))))
    if seed_flat is not None and not np.all(np.isfinite(seed_flat)):
        raise ValueError("{}: seed_vec has non-finite entries".format(solver))
    seed_norm = 0.0 if seed_flat is None else float(np.linalg.norm(seed_flat))
    if seed_norm == 0.0:
        seed_flat = None                    # a zero seed is no seed

    leaders = []
    cands = []          # (value, rank, index, gap, pp_sector)
    failed = []
    for op in ops:
        p = int(op.pp_sector)
        rank = _SO_PP_ORDER.index(p)
        basis = op.basis
        dim = int(basis.dim)
        if dim == 0:
            leaders.append(SectorLeader(p, 0, complex(np.nan), "empty", True,
                                        0.0, 0.0, "exact", 0))
            continue
        red = basis.reduce(op.operator, adjoint=bool(op.has_adjoint))
        n_ret = min(k, dim)
        if _operator_is_exactly_zero(red, dim):
            gaps = _so_lifted_unit_gaps(basis, n_ret, gap_shape)
            for i in range(n_ret):
                cands.append((0j, rank, i, gaps[i], p))
            leaders.append(SectorLeader(p, dim, 0j, "null-operator", True, 0.0,
                                        0.0, "exact", n_ret))
            continue
        seed = None
        if seed_flat is not None:
            part = basis.to_reduced(seed_flat[basis.indices]
                                    if basis.indices is not None else seed_flat)
            part_norm = float(np.linalg.norm(part))
            if part_norm > _SO_SEED_PART_RTOL * seed_norm:
                seed = part
            elif part_norm > 0.0:
                logger.info("spin-orbital sector %s: the seed part (norm %.3e, %.3e of "
                            "the seed) is round-off; the sector is solved unseeded.",
                            _pp_label(p), part_norm, part_norm / seed_norm)
        M = None
        if dim <= max(2 * k + 1, 20):
            try:
                vals, vecs, M = _so_dense_reduced_eigen(red, dim, k, seed)
                selection = "dense-LR"
            except np.linalg.LinAlgError as exc:
                logger.warning("spin-orbital sector %s: the dense "
                               "diagonalization failed (%s).", _pp_label(p), exc)
                vals = np.full(n_ret, np.nan + 0j)
                vecs = np.zeros((dim, n_ret), dtype=complex)
                selection = "not-converged"
        else:
            info = {}
            vals, vecs = _so_solve_eigenvalue(
                op.operator, basis, k, method, spectral_shift=spectral_shift,
                seed_vec=seed, map_back=False, info_out=info)
            selection = info.get("selection")
            if selection == "not-converged" or not np.all(np.isfinite(vals)):
                selection = "not-converged"
        vals = np.asarray(vals)
        vecs = np.asarray(vecs)
        kind = ("exact" if dim <= max(2 * k + 1, 20)
                else "adjoint" if op.has_adjoint else "forward")
        if selection == "not-converged":
            logger.warning("spin-orbital sector %s: the eigen solve did not "
                           "converge; the sector is failed.", _pp_label(p))
            leaders.append(SectorLeader(p, dim, complex(np.nan), selection,
                                        False, float("inf"), float("nan"),
                                        kind, len(vals)))
            failed.append(p)
            continue
        if kind == "exact":
            est = float(np.linalg.norm(M, 2))
            images = M @ vecs
        else:
            est, steps, change, conv = _so_sector_norm_estimate(
                red, dim, op.has_adjoint)
            images = np.asarray(red.matmat(vecs))
        scale = max(1.0, float(np.max(np.abs(vals))), est)
        if kind == "forward" or (kind == "adjoint" and not conv):
            logger.info(
                "spin-orbital sector %s: operator-norm estimate %s after %d "
                "steps (relative change %.2e); residual scale %.6e",
                _pp_label(p), kind, steps, change, scale)
        resid = 0.0
        for i in range(len(vals)):
            v = vecs[:, i]
            n_k = max(1.0, abs(vals[i]), est)
            r = (float(np.linalg.norm(images[:, i] - vals[i] * v))
                 / (n_k * float(np.linalg.norm(v))))
            resid = max(resid, r) if np.isfinite(r) else float("inf")
        certified = bool(resid <= _SO_SECTOR_RESIDUAL_TOL)
        lead_value = complex(vals[np.argmax(vals.real)]) if len(vals) else complex(np.nan)
        if not certified:
            logger.warning(
                "spin-orbital sector %s: a returned eigenpair fails the "
                "residual check (%.3e > %.0e relative to the scale %.6e); "
                "the sector is failed.", _pp_label(p), resid,
                _SO_SECTOR_RESIDUAL_TOL, scale)
            leaders.append(SectorLeader(p, dim, lead_value, "residual-failed",
                                        False, resid, est, kind, len(vals)))
            failed.append(p)
            continue
        leaders.append(SectorLeader(p, dim, lead_value, selection, True, resid,
                                    est, kind, len(vals)))
        for i in range(len(vals)):
            g = basis.scatter(basis.to_full(vecs[:, i])).reshape(gap_shape)
            cands.append((complex(vals[i]), rank, i, g, p))

    if all(lead.dim == 0 for lead in leaders):
        logger.warning("spin-orbital: the antisymmetric gap space is empty on "
                       "this grid; no eigenvalue is reported")
    kept = _so_union(cands, k)
    if seed_flat is not None and kept:
        kept = _so_order_by_seed(kept, seed_flat / seed_norm)
    vals = np.array([c[0] for c in kept], dtype=complex)
    gaps = (np.array([c[3] for c in kept], dtype=complex) if kept
            else np.zeros((0,) + gap_shape, dtype=complex))
    pp = np.array([c[4] for c in kept], dtype=np.int8)
    converged = bool(not failed and len(kept) > 0)
    leading = complex(vals[0]) if converged else complex(np.nan)
    selection = ";".join("{}:{}".format(_pp_label(lead.pp_sector), lead.selection)
                         for lead in leaders)
    return SectorSolve(vals, gaps, pp, tuple(leaders), tuple(failed), converged,
                       leading, selection)


def _so_iterate_sectors(sector_ops, sigma_init, *, gap_shape, max_iter, alpha,
                        tol, spectral_shift=None):
    """The power iteration of the spin-orbital kernel sector by sector.

    Per sector (order ``+1, 0, -1``): ``dim = 0`` is skipped (``"empty"``);
    an exactly zero reduced operator is ``"null"`` (``lambda = 0``, the first
    lifted unit vector, not iterated); ``dim = 1`` is solved exactly;
    otherwise ``_solve_leading(..., "iteration")`` runs on the reduced sector
    operator from the sector's part of ``sigma_init`` (a zero part starts
    from the complex ``default_rng(0)`` vector of the reduced space), and a
    sector whose run did not converge is failed. The reduction is
    :class:`SectorIteration`'s."""
    solver = "_so_iterate_sectors"
    ops = _so_canonical_sector_ops(sector_ops, solver)
    spectral_shift = _normalize_spectral_shift(spectral_shift)
    _validate_spectral_shift(spectral_shift, "iteration")
    gap_shape = tuple(int(x) for x in gap_shape)
    for op in ops:
        _so_check_sector_op(op, gap_shape, solver)
    init_flat = np.asarray(sigma_init).ravel()
    if init_flat.size != int(np.prod(gap_shape)):
        raise ValueError("{}: sigma_init has {} entries, the gap {} has {}".format(
            solver, init_flat.size, gap_shape, int(np.prod(gap_shape))))

    leaders = []
    for op in ops:
        p = int(op.pp_sector)
        basis = op.basis
        dim = int(basis.dim)
        if dim == 0:
            leaders.append(IterationLeader(p, "empty", float("nan"), None, 0))
            continue
        red = basis.reduce(op.operator)
        if _operator_is_exactly_zero(red, dim):
            gap = _so_lifted_unit_gaps(basis, 1, gap_shape)[0]
            leaders.append(IterationLeader(p, "null", 0.0, gap, 0))
            continue
        if dim == 1:
            # a one-dimensional sector: its single element is the eigenvalue
            lam = complex(np.asarray(red.matvec(np.ones(1, dtype=complex))).ravel()[0])
            gap = _so_lifted_unit_gaps(basis, 1, gap_shape)[0]
            if np.isfinite(lam):
                leaders.append(IterationLeader(p, "converged", lam, gap, 0))
            else:
                logger.warning("spin-orbital sector %s: the matrix element of the "
                               "one-dimensional sector is not finite (%s); the sector "
                               "is failed.", _pp_label(p), lam)
                leaders.append(IterationLeader(p, "failed", float("nan"), gap, 0))
            continue
        part = basis.to_reduced(init_flat[basis.indices]
                                if basis.indices is not None else init_flat)
        if not np.any(part != 0):
            rng = np.random.default_rng(0)
            part = rng.standard_normal(dim) + 1j * rng.standard_normal(dim)
        value, vec, info = _solve_leading(
            lambda: (red, dim), dim, "iteration", max_iter=max_iter,
            convergence_tol=tol, init_vec=part, alpha=alpha,
            spectral_shift=spectral_shift)
        gap = basis.scatter(basis.to_full(np.asarray(vec).ravel())).reshape(gap_shape)
        n_iter = int(info.get("n_iter", 0))
        if info.get("converged") and np.isfinite(value):
            leaders.append(IterationLeader(p, "converged", value, gap, n_iter))
        else:
            leaders.append(IterationLeader(p, "failed", float("nan"), gap, n_iter))

    leaders = tuple(leaders)
    per_sector = tuple(lead.n_iter for lead in leaders)
    n_iter = max(per_sector) if per_sector else 0
    failed = tuple(lead.pp_sector for lead in leaders if lead.status == "failed")
    if failed:
        first = next(lead for lead in leaders if lead.status == "failed")
        return SectorIteration(float("nan"), first.gap, first.pp_sector, leaders,
                               failed, False, n_iter, per_sector, "iteration")
    live = [lead for lead in leaders if lead.status != "empty"]
    if not live:
        # the empty-space WARNING is the eigen leg's (one per run)
        return SectorIteration(float("nan"), None, _SO_PP_DENSE, leaders, (),
                               False, n_iter, per_sector, "iteration")
    best = _so_union([(complex(lead.value), _SO_PP_ORDER.index(lead.pp_sector), 0,
                       None, lead.pp_sector) for lead in live], 1)[0]
    lead = next(x for x in live if x.pp_sector == best[4])
    selection = "null-operator" if lead.status == "null" else "iteration"
    return SectorIteration(lead.value, lead.gap, lead.pp_sector, leaders, (),
                           True, n_iter, per_sector, selection)


def _so_bond_report(solver_mode, eigen, iteration, *, gap_shape):
    """The one report of a spin-orbital run (:class:`SectorReport`) from its
    legs (:class:`EigenLeg` for ``"eigenvalue"`` and ``"both"``,
    :class:`IterationLeg` for ``"iteration"`` and ``"both"``; the other
    ``None``).

    * ``"eigenvalue"``: the eigen leg's leader (``reported = Re leading``,
      NaN if not converged); the gap of the leader, else of the first listed
      pair (``"uncertified"``), else zeros (``"none"``).
    * ``"iteration"``: the iteration value and gap (``eigen_converged`` is the
      iteration's convergence, which needs a finite value); zeros and
      ``"none"`` when it has no gap.
    * ``"both"``: the eigen leg decides ``reported``, ``eigen_converged``,
      the listing and ``failed_sectors``; the gap is the eigen leader's when
      the eigen leg converged, else the iteration's (``"uncertified"``),
      else zeros; ``iteration_converged`` is the iteration's convergence
      (again with a finite value).

    A listing whose first value is NaN (the failed-solve sentinel of the
    dense layout) holds no candidate: in ``"eigenvalue"`` mode its gap is
    zeros with ``"none"``.

    ``gap_status`` is ``"certified"`` exactly when ``reported`` is finite."""
    if solver_mode not in ("eigenvalue", "iteration", "both"):
        raise ValueError("_so_bond_report: solver_mode must be \"eigenvalue\", "
                         "\"iteration\" or \"both\", got {!r}".format(solver_mode))
    if solver_mode in ("eigenvalue", "both") and eigen is None:
        raise ValueError("_so_bond_report: solver_mode {!r} needs the eigenvalue "
                         "leg".format(solver_mode))
    if solver_mode in ("iteration", "both") and iteration is None:
        raise ValueError("_so_bond_report: solver_mode {!r} needs the iteration "
                         "leg".format(solver_mode))
    gap_shape = tuple(int(x) for x in gap_shape)

    def zeros():
        return np.zeros(gap_shape, dtype=complex)

    def shaped(g):
        g = np.asarray(g)
        return g if g.shape == gap_shape else g.reshape(gap_shape)

    nan = float("nan")
    iteration_converged = None
    iteration_failed = None
    def iteration_ok():
        # convergence of the iteration leg needs a finite value
        return bool(iteration.converged) and bool(np.isfinite(iteration.value))

    if solver_mode == "iteration":
        converged = iteration_ok()
        reported = float(np.real(iteration.value)) if converged else nan
        if iteration.gap is None:
            gap, pp = zeros(), _SO_PP_DENSE
        else:
            gap, pp = shaped(iteration.gap), int(iteration.pp_sector)
        selection = iteration.selection
        vals = np.zeros(0, dtype=complex)
        gaps = np.zeros((0,) + gap_shape, dtype=complex)
        pps = np.zeros(0, dtype=np.int8)
        failed = tuple(iteration.failed_sectors)
        iteration_failed = failed
        has_candidate = iteration.gap is not None
    else:
        converged = bool(eigen.converged)
        reported = float(np.real(eigen.leading)) if converged else nan
        if not np.isfinite(reported):
            converged = False
        n = len(eigen.vals)
        if converged:
            gap, pp = shaped(eigen.gaps[0]), int(eigen.pp_sector[0])
            has_candidate = True
        elif solver_mode == "both":
            has_candidate = iteration.gap is not None
            gap, pp = ((shaped(iteration.gap), int(iteration.pp_sector))
                       if has_candidate else (zeros(), _SO_PP_DENSE))
        else:
            # a NaN-valued listing (the failed-solve sentinel of the dense
            # layout) holds placeholders, not candidates
            has_candidate = n > 0 and bool(np.isfinite(eigen.vals[0]))
            gap, pp = ((shaped(eigen.gaps[0]), int(eigen.pp_sector[0]))
                       if has_candidate else (zeros(), _SO_PP_DENSE))
        if solver_mode == "both":
            iteration_converged = iteration_ok()
            iteration_failed = tuple(iteration.failed_sectors)
        selection = eigen.selection
        vals, gaps, pps = eigen.vals, eigen.gaps, eigen.pp_sector
        failed = tuple(eigen.failed_sectors)
    if np.isfinite(reported):
        status = "certified"
    elif has_candidate:
        status = "uncertified"
    else:
        status = "none"
    return SectorReport(reported, gap, status, pp,
                        bool(converged), iteration_converged, selection, vals,
                        gaps, pps, failed, iteration_failed)


def _so_check_antisymmetric_eigenpairs(vals, match):
    """Defensive check of the spin-orbital eigenpairs (``match`` from
    :func:`_so_is_antisymmetric`): refuse when the LEADING eigenvector --
    the gap written under the fermion-antisymmetric header -- is not
    antisymmetric, and warn about any further pair that is not (match
    column 0). The solvers work on antisymmetric coordinates, so this fires
    only if that construction is broken."""
    if np.size(vals) and not np.all(np.isfinite(np.asarray(vals))):
        # the eigenvalue solver could not determine the largest real part
        # (NaN sentinel, _solve_leading): the vectors are placeholders, not
        # eigenvectors, and the run is reported with eigen_converged = false
        logger.warning(
            "spin-orbital: no eigenvalue was determined (NaN); the written gap "
            "is a placeholder, not an eigenvector.")
        return
    match = np.asarray(match, dtype=bool)
    if match.size and not match[0]:
        raise ValueError(
            "[mode] enable_spin_orbital = true: the leading eigenvector is not "
            "fermion-antisymmetric (Delta_AB(k) = -Delta_BA(-k)) to round-off; "
            "it is Pauli-forbidden and is not written (eigenvalues {})."
            .format(np.asarray(vals)))
    bad = np.where(~match)[0]
    if bad.size:
        logger.warning(
            "spin-orbital eigenpairs %s are not fermion-antisymmetric to "
            "round-off (eigenvalues %s); they are written with match = 0.",
            bad.tolist(), np.asarray(vals)[bad])
