"""Linearized Eliashberg equation solver for superconducting instability analysis.

This module implements a post-processing tool that uses the bare susceptibility
chi0q from H-wave's RPA solver to solve the linearized Eliashberg equation
and analyze superconducting instabilities.

The tool reads chi0q.npz output from RPA calculations along with interaction
parameters, reconstructs the Green's function, computes RPA vertices,
and solves the Eliashberg equation via self-consistent iteration or
eigenvalue analysis.
"""

from __future__ import annotations

import os
import sys
import logging

import numpy as np

from hwave.solver.vertex_table import sc_coefficients
from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.declarations import symmetrise_k
from numpy.fft import fftn, ifftn
from scipy.optimize import bisect
from scipy.sparse.linalg import LinearOperator, eigs, bicgstab, gmres, lgmres

import hwave
import hwave.qlmsio.wan90 as wan90
from hwave.solver.rpa import validate_chi0q_index_convention
from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver import backend

logger = logging.getLogger("hwave_sc")

# Number of independent random probes used by _solve_iteration to test whether
# the Eliashberg kernel commutes with parity before enabling gap-parity
# projection. More than one guards against a single probe accidentally lying
# near the nullspace of the commutator and missing a non-centrosymmetric kernel.
_PARITY_GUARD_PROBES = 3

# Default Matsubara grid size (mode.param.Nmat). Single definition so the
# legacy-file fallback in _static_freq_position and the main solver loop can
# never silently disagree on the grid.
_DEFAULT_NMAT = 1024


# ---------------------------------------------------------------------------
# Eliashberg frequency mode dispatch
# ---------------------------------------------------------------------------

def _eliashberg_frequency(input_dict):
    """Validate and return the eliashberg frequency mode.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Returns
    -------
    str
        Either "static" (default) or "dynamic".

    Raises
    ------
    ValueError
        If the frequency mode is not in the allowed set.
    """
    freq = input_dict.get("eliashberg", {}).get("frequency", "static")
    if freq not in ("static", "dynamic"):
        raise ValueError(
            "eliashberg.frequency must be 'static' or 'dynamic', got '{}'"
            .format(freq))
    return freq


def _validate_dynamic_prereqs(input_dict):
    """Validate prerequisites for dynamic Eliashberg calculation.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.

    Raises
    ------
    ValueError
        If chi0q_mode is not "flex" or if Nmat is odd.
    """
    eli = input_dict.get("eliashberg", {})
    if eli.get("chi0q_mode") != "flex":
        raise ValueError(
            "eliashberg.frequency='dynamic' requires chi0q_mode='flex' "
            "(full-frequency chiq_s/chiq_c and a dressed green are only "
            "produced by the FLEX path); got chi0q_mode='{}'"
            .format(eli.get("chi0q_mode")))
    nmat = int(input_dict["mode"]["param"].get("Nmat", 1024))
    if nmat % 2 != 0:
        raise ValueError(
            "eliashberg.frequency='dynamic' requires an even Nmat "
            "(centered Matsubara grid); got Nmat={}".format(nmat))


def _warn_if_static_ignores_channel_flags(eli_param):
    """Warn that the channel-decomposition flags are inert on the static path.

    ``zero_chi_s``/``zero_chi_c`` split the DYNAMIC pairing vertex into its
    spin-/charge-fluctuation and instantaneous-bare pieces (see
    ``eliashberg_dynamic``). The static solver builds its kernel from the RPA/
    FLEX static susceptibility and never reads these flags, so a user who sets
    one with ``frequency='static'`` gets no decomposition. Warn loudly instead
    of silently ignoring the request.
    """
    from hwave.solver.backend import as_bool

    ignored = [name for name in ("zero_chi_s", "zero_chi_c")
               if as_bool(eli_param.get(name, False))]
    if ignored:
        logger.warning(
            "[eliashberg] %s set but frequency='static'; the channel-"
            "decomposition flags only apply to frequency='dynamic' and are "
            "ignored here.", ", ".join(ignored))


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _static_freq_position(freq_index, nfreq, config_nmat, file_name,
                          file_nmat=None):
    """Locate the zero bosonic frequency along the chi0q frequency axis.

    RPA writes ``freq_index`` (the ORIGINAL Matsubara indices, zero frequency
    at original index Nmat//2) into chi0q.npz; the matsubara_frequency option
    can restrict the stored axis, so the static slice is generally NOT at
    nfreq//2 of the stored array.  Newer files also record ``nmat`` (the full
    grid size of the producing run), which resolves the position without any
    assumption about the hwave_sc configuration.

    Parameters
    ----------
    freq_index : array-like or None
        The freq_index metadata from the npz file (None for legacy files).
    nfreq : int
        Length of the stored frequency axis.
    config_nmat : int
        Nmat from mode.param of the hwave_sc run (fallback only).
    file_name : str
        For error messages.
    file_nmat : int or None
        The nmat metadata from the npz file, if present.

    Returns
    -------
    int or None
        Position of the zero bosonic frequency along the stored axis, or
        None when the file carries no usable metadata -- the caller must
        then slice the CENTER of the frequency axis it actually uses (only
        the caller can identify that axis reliably, e.g. for the 6D
        reference format).
    """
    if freq_index is None:
        logger.warning(
            "chi0q file '{}' has no freq_index metadata; using the center "
            "of the stored frequency axis as the static slice.".format(
                file_name))
        return None

    freq_index = np.asarray(freq_index).ravel()
    if freq_index.size == 0:
        raise ValueError(
            "chi0q file '{}' records no frequencies (matsubara_frequency="
            "\"none\"?); it cannot provide the static susceptibility."
            .format(file_name))
    if freq_index.size != nfreq:
        # Legacy FLEX files store the FULL grid but a restricted freq_index
        # (FLEX never applied the matsubara_frequency filter): the DATA axis
        # is authoritative, so fall back to the pre-metadata behavior.
        logger.warning(
            "chi0q file '{}': freq_index length {} does not match the "
            "frequency axis length {}; ignoring the metadata and using the "
            "center of the stored frequency axis as the static slice."
            .format(file_name, freq_index.size, nfreq))
        return None

    def _find(nmat_orig):
        # the zero bosonic frequency has ORIGINAL index nmat_orig//2
        pos = np.where(freq_index == nmat_orig // 2)[0]
        if pos.size == 0:
            raise ValueError(
                "chi0q file '{}' does not contain the zero bosonic frequency "
                "(original index Nmat//2 = {}; stored freq_index = {}..{}). "
                "Regenerate chi0q with a matsubara_frequency range covering "
                "Nmat//2.".format(file_name, nmat_orig // 2,
                                  freq_index[0], freq_index[-1]))
        return int(pos[0])

    if file_nmat is not None:
        # authoritative: the producing run recorded its own grid size
        return _find(int(file_nmat))

    if np.array_equal(freq_index, np.arange(nfreq)):
        # Legacy file (no nmat metadata) with a contiguous 0-based
        # freq_index: EITHER a full grid of its own (zero at nfreq//2) OR a
        # matsubara_frequency=[0,K] restriction of SOME run's grid (config
        # or larger).  Resolve only when the readings coincide; silently
        # picking either would return a finite frequency whenever the other
        # reading is the true one.
        if nfreq == config_nmat:
            return nfreq // 2
        raise ValueError(
            "chi0q file '{}' is ambiguous: freq_index = 0..{} with "
            "mode.param.Nmat = {} could be either a full {}-point grid "
            "(zero frequency at {}) or a [0,{}] restriction of a run with "
            "a different Nmat (zero frequency elsewhere or absent). If the "
            "file holds a full grid, set mode.param.Nmat = {}; otherwise "
            "regenerate it with a newer version (which records nmat in "
            "the file)."
            .format(file_name, nfreq - 1, config_nmat, nfreq, nfreq // 2,
                    nfreq - 1, nfreq))

    # Legacy restricted range: fall back to the configured Nmat, which must
    # be the producing run's grid size (the standard workflow shares one
    # TOML).  Indices reaching past the config grid disprove that reading,
    # so refuse to guess (looking up config_nmat//2 could land on a finite
    # frequency of the true, larger grid).
    if int(freq_index.max()) >= config_nmat:
        raise ValueError(
            "chi0q file '{}': freq_index reaches {} >= mode.param.Nmat = "
            "{}, so the file was produced by a run with a different Nmat "
            "and the zero-frequency position cannot be determined. Set "
            "mode.param.Nmat to the producing run's value, or regenerate "
            "the file with a newer version (which records nmat)."
            .format(file_name, int(freq_index.max()), config_nmat))
    return _find(config_nmat)


def _read_freq_meta(data):
    """Decode the frequency metadata of a chi0q/chiq npz file.

    Returns
    -------
    (freq_index or None, nmat or None)
        Single definition shared by all loaders so chi0q and chiq_s/chiq_c
        can never interpret the same file's frequency axis differently.
    """
    freq_index = data["freq_index"] if "freq_index" in data else None
    file_nmat = int(data["nmat"]) if "nmat" in data else None
    return freq_index, file_nmat


def _reject_ir_native(data, file_name, hint):
    """Fail fast when a uniform-grid reader is handed sparse-node data
    (design ir-matsubara-stage3.md Sec. 4). MUST run before any freq_index
    metadata interpretation: the legacy fallbacks (missing freq_index ->
    center row) would otherwise silently read a sparse node as the static
    slice."""
    if is_ir_native(data):
        raise ValueError(
            "file '{}' holds sparse-IR node data "
            "(frequency_grid=sparse_ir_nodes); {}".format(file_name, hint))


def _load_chi0q(input_dict):
    """Load chi0q from NPZ file produced by H-wave RPA solver.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array.
    static_index : int
        Position of the zero bosonic frequency along the frequency axis
        (determined from the freq_index metadata).
    """
    output_info = input_dict["file"]["output"]
    path_to_output = output_info["path_to_output"]
    chi0q_name = output_info.get("chi0q", "chi0q")
    file_name = os.path.join(path_to_output, chi0q_name + ".npz")

    logger.info("Loading chi0q from {}".format(file_name))
    data = np.load(file_name)
    _reject_ir_native(
        data, file_name,
        "the static Eliashberg solver requires uniform-grid files. Re-run "
        "FLEX with [mode.param] write_densified = true, or switch to "
        "frequency = \"dynamic\" with [eliashberg] matsubara_basis = "
        "\"ir\".")
    enable_spin_orbital = input_dict.get("mode", {}).get(
        "enable_spin_orbital", False)
    validate_chi0q_index_convention(data, enable_spin_orbital, file_name)
    chi0q = data["chi0q"]
    logger.info("chi0q shape: {}".format(chi0q.shape))

    # coeff_tail provenance check (issue #80): the Matsubara tail correction
    # changes chi0q at O(1) (lambda differed by 50% in the reproduction), so
    # a config whose effective coeff_tail differs from the file's would
    # recompute DIFFERENT physics under chi0q_mode="calc". Warn so load-vs-
    # calc comparisons are not silently inconsistent. Files without the key
    # (older builds) load silently as before.
    if "coeff_tail" in data:
        file_tail = float(data["coeff_tail"])
        config_tail = float(input_dict.get("mode", {}).get("param", {})
                            .get("coeff_tail", 0.0))
        if file_tail != config_tail:
            logger.warning(
                "chi0q file '{}' was produced with coeff_tail = {} but this "
                "config's effective value is coeff_tail = {}; results are "
                "NOT comparable with a chi0q_mode=\"calc\" run under this "
                "config. Set [mode.param] coeff_tail = {} to match the "
                "file.".format(file_name, file_tail, config_tail, file_tail))

    freq_index, file_nmat = _read_freq_meta(data)
    # Identify the frequency axis from the array LAYOUT, never from the
    # freq_index length (a restricted freq_index can coincidentally match
    # an orbital axis).  4D raw: axis 0.  8D ref: last axis.  6D is either
    # raw (nmat, nvol, norb^4; the last four axes are equal) or ref
    # (norb, norb, Nx, Ny, Nz, nmat; the first two axes are equal) --
    # disambiguate structurally, with the freq_index length only as the
    # tiebreaker for degenerate shapes.  Without metadata
    # _static_freq_position returns None and the caller slices the center
    # of the axis it actually uses.
    if freq_index is not None:
        nfi = np.asarray(freq_index).size
        if chi0q.ndim == 4:
            nfreq = chi0q.shape[0]
        elif chi0q.ndim == 8:
            nfreq = chi0q.shape[-1]
        elif chi0q.ndim == 6:
            raw_like = len(set(chi0q.shape[2:])) == 1
            ref_like = chi0q.shape[0] == chi0q.shape[1]
            if raw_like and not ref_like:
                nfreq = chi0q.shape[0]
            elif ref_like and not raw_like:
                nfreq = chi0q.shape[-1]
            elif nfi == chi0q.shape[0]:
                nfreq = chi0q.shape[0]
            else:
                nfreq = chi0q.shape[-1]
        else:
            nfreq = chi0q.shape[-1]
    else:
        nfreq = 0  # unused: no metadata -> _static_freq_position gives None
    config_nmat = input_dict.get("mode", {}).get("param", {}).get("Nmat",
                                                                _DEFAULT_NMAT)
    static_index = _static_freq_position(freq_index, nfreq, config_nmat,
                                         file_name, file_nmat=file_nmat)
    if static_index is not None and static_index != nfreq // 2:
        logger.info("static (zero bosonic frequency) slice at index {} "
                    "of {}".format(static_index, nfreq))
    return chi0q, static_index


def _calc_chi0q_internal(input_dict, chi0q_tensor="auto",
                         precomputed_mu=None):
    """Compute chi0q internally using H-wave's RPA module.

    Instead of loading chi0q from a file, this function creates an RPA solver
    instance and computes chi0q from scratch using the Transfer and interaction
    parameters specified in the TOML config.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.
    chi0q_tensor : str
        Index structure of chi0q. Options:
        - "reduced": 2-index (nmat, nvol, norb, norb). Exact for single-orbital
          and for multi-orbital with CoulombIntra only (no inter-orbital
          interactions).
        - "general": 4-index (nmat, nvol, norb, norb, norb, norb). Required
          when inter-orbital interactions (CoulombInter, Hund, Exchange) are
          present, as their S/C matrices couple different orbital-pair indices.
        - "auto": Use "general" if CoulombInter/Hund/Exchange are present,
          otherwise "reduced".
    precomputed_mu : float, optional
        If provided, skip chemical potential determination and use this value
        directly. This avoids redundant eigenvalue decomposition and mu search
        when the caller has already computed mu (e.g. calc_eliashberg).

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array. Shape depends on chi0q_tensor mode:
        - reduced: (nmat, nvol, norb, norb)
        - general: (nmat, nvol, norb, norb, norb, norb)
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Determine calc_scheme from chi0q_tensor option
    if chi0q_tensor == "auto":
        # Use "general" when inter-orbital interactions are present,
        # because their S/C matrices have off-diagonal elements that
        # couple to chi0q off-diagonal components.
        # With CoulombIntra only, S is block-diagonal and reduced is exact.
        # CaseInsensitiveDict, like _read_interaction_files (round-5
        # review): an exact-case check classified a 'coulombinter' run as
        # reduced, silently omitting the very components the comment
        # above says inter-orbital interactions require. 'Coulomb' counts
        # too -- its split can produce a CoulombInter part.
        from requests.structures import CaseInsensitiveDict
        files = CaseInsensitiveDict(info_inputfile.get("interaction", {}))
        has_interorbital = any(k in files for k in
                              ["Hund", "Exchange", "CoulombInter",
                               "Ising", "PairHop", "Coulomb"])
        if has_interorbital:
            calc_scheme = "general"
        else:
            calc_scheme = "reduced"
    elif chi0q_tensor == "general":
        calc_scheme = "general"
    else:
        calc_scheme = "reduced"

    info_mode_rpa = dict(info_mode)
    info_mode_rpa["calc_scheme"] = calc_scheme

    logger.info("Computing chi0q internally (calc_scheme={})...".format(calc_scheme))

    # Read input files via QLMSkInput
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")

    # Create RPA solver (sets up Lattice, Interaction, parameters)
    solver = sol_rpa.RPA(ham_info, info_log, info_mode_rpa)

    # Get green_info (may contain chi0q_init, trans_mod, green_init)
    green_info = read_io.get_param("green")
    green_info.update(solver.read_init(info_inputfile))

    # Compute chi0q: eigenvalues -> mu -> Green's function -> chi0q
    beta = 1.0 / solver.T

    solver._calc_epsilon_k(green_info)

    if precomputed_mu is not None:
        # Reuse mu from caller to skip redundant bisection
        mu = precomputed_mu
        logger.info("Using precomputed mu = {}".format(mu))
    elif solver.calc_mu:
        if solver.spin_mode == "spin-free":
            Ncond = solver.Ncond / 2
        else:
            Ncond = solver.Ncond
        dist, mu = solver._find_mu(Ncond, solver.T)
    else:
        mu = solver.mu_value

    green0, green0_tail = solver._calc_green(beta, mu)

    chi0q = solver._calc_chi0q(green0, green0_tail, beta)

    # For spin-free mode, chi0q has shape (1, nmat, nvol, ...)
    # Remove the block index
    if solver.spin_mode in ["spin-free", "spinful"]:
        assert chi0q.shape[0] == 1
        chi0q = chi0q[0]
    else:
        # spin-diag: shape (2, nmat, nvol, ...)
        assert chi0q.shape[0] == 2
        pass

    logger.info("chi0q computed internally, shape: {}".format(chi0q.shape))
    return chi0q


def _read_interaction_files(input_dict):
    """Read Transfer and interaction files using wan90 module.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.

    Returns
    -------
    geom_info : dict
        Geometry information (norb, rvec, center).
    hr : dict
        Transfer (hopping) parameters.
    interactions : dict
        Dictionary of interaction parameters keyed by type name.
    """
    # Case-insensitive on purpose (round-4 review): the FLEX reader stores
    # these keys in a CaseInsensitiveDict, so a run configured with e.g.
    # 'coulombinter' PRODUCES a susceptibility with that interaction while
    # an exact-case lookup here would silently read an EMPTY interaction
    # set -- the compatibility gate and the pairing vertex would both be
    # built from nothing.
    from requests.structures import CaseInsensitiveDict
    files = CaseInsensitiveDict(input_dict["file"]["input"]["interaction"])
    path_to_input = files.get("path_to_input", "")

    geom_file = os.path.join(path_to_input, files["Geometry"])
    geom_info = wan90.read_geom(geom_file)
    logger.info("norb = {}".format(geom_info["norb"]))

    transfer_file = os.path.join(path_to_input, files["Transfer"])
    hr = wan90.read_w90(transfer_file)

    interaction_types = ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                        "Ising", "PairLift", "PairHop"]
    interactions = {}
    from hwave.solver.declarations import validate_hermitian_closure

    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            tbl = wan90.read_w90(f)
            # issue #93: fail fast on non-Hermitian-closed declarations,
            # with the same rule and tolerance as the k-space reader
            validate_hermitian_closure(itype, tbl, source=f)
            interactions[itype] = tbl

    # combined 'Coulomb' input: same decomposition as UHFk/RPA
    # (wan90.split_coulomb; r=0 diagonal -> CoulombIntra, rest -> CoulombInter)
    if "Coulomb" in files:
        if "CoulombIntra" in interactions or "CoulombInter" in interactions:
            raise ValueError(
                "Coulomb cannot be specified together with "
                "CoulombIntra or CoulombInter")
        f = os.path.join(path_to_input, files["Coulomb"])
        logger.info("Reading Coulomb from {}".format(f))
        _coulomb_tbl = wan90.read_w90(f)
        validate_hermitian_closure("Coulomb", _coulomb_tbl, source=f)
        coulomb_intra, coulomb_inter = wan90.split_coulomb(_coulomb_tbl)
        if coulomb_intra:
            interactions["CoulombIntra"] = coulomb_intra
        if coulomb_inter:
            interactions["CoulombInter"] = coulomb_inter

    return geom_info, hr, interactions


# ---------------------------------------------------------------------------
# k-space construction
# ---------------------------------------------------------------------------

def _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb):
    """Build Hamiltonian epsilon(k) from real-space transfer integrals.

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    hr : dict
        Transfer parameters from wan90.read_w90().
    norb : int
        Number of orbitals.

    Returns
    -------
    epsilon_k : ndarray
        Hamiltonian in k-space, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    epsilon_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    # Solver-core convention (rpa.py _make_ham_trans: tab_r[R,orb1,orb2] +
    # fftn == e^{-ikR}): epsilon[a,b](k) = sum_R t_R[a,b] e^{-ikR}. This keeps
    # sc-built quantities element-wise consistent with arrays loaded from
    # FLEX/RPA output files. (The previous [orb2,orb1] + e^{+ikR} form is the
    # orbital transpose at -k; identical for real hoppings, different for
    # complex Hermitian ones.)
    for (irvec, orbvec), value in hr.items():
        orb1, orb2 = orbvec
        Rx, Ry, Rz = irvec
        epsilon_k[orb1, orb2, :, :, :] += value * np.exp(
            -1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
        )
    return epsilon_k


def _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb):
    """Build interaction matrices in k-space.

    Supports CoulombIntra (U), CoulombInter (U'), Hund (J), and Exchange (J').

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    interactions : dict
        Interaction parameters keyed by type.
    norb : int
        Number of orbitals.

    Returns
    -------
    inter_k : dict
        Dictionary of interactions in k-space.
        Keys: "CoulombIntra", "CoulombInter", "Hund", "Exchange".
        Values: ndarray of shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )

    def _to_k(value_r, transpose=True):
        # Same Fourier phase as the solver core, e^{-iqR}, but the ORBITAL
        # PAIR is stored transposed: an entry (R, (a, b)) lands at [b, a].
        #
        # The interaction is a four-index object, and its MATRIX form carries a
        # pair-index transpose that the one-body Hamiltonian does not. H-wave's
        # own paper (arXiv:2308.00324) makes this explicit: Eq.(12) defines the
        # tensor W_ij^{aa'bb'} with the first pair at site i and the second at
        # site j, Eq.(16) puts the first pair at momentum +q, and Eq.(21) then
        # defines the matrix used in the RPA equation as
        #     [W(q)]^{ab} = W_q^{ba}
        # so for a density-density term W^{aabb} = V_ab(R) the matrix element is
        # [W(q)]^{(aa),(bb)} = V_ba(q), i.e. the matrix is V(q)^T.
        #
        # `rpa.py::_make_ham_inter` has always built that (it stores the entry's
        # first orbital in the SECOND pair slot); this builder did not, so the
        # susceptibility was solved with one orientation and the pairing vertex
        # assembled from the other (issue #96).
        #
        # Both consumers need the transpose, confirmed by exact
        # diagonalization on a bond set with V_ab(R) != V_ba(R):
        #   * the RPA ladder [I + chi0 W]^-1 chi0 -- solved here too, in
        #     _compute_vertices_simple -- selects V^T with a residual that
        #     scales linearly in V (pure O(V^2) truncation), against 98% for V;
        #   * the bare pair-scattering amplitude
        #     <k' a up, -k' b down| H_int |k a up, -k b down>, which is exact at
        #     first order, matches V^T to 6e-16, against 99% for V.
        #
        # NOTE: `_build_hamiltonian_k` is a one-body object and correctly keeps
        # the plain [a, b] placement. Do not "align" this builder to it.
        #
        # PairHop is EXCLUDED (transpose=False below), and unlike the rest that
        # is now a MEASURED result rather than a gap. It is the one type
        # `rpa.py` does not place through `_append_inter`: `_append_pairhop`
        # uses the slots (b, a, a, b) rather than the density-density
        # (b, b, a, a), and in the S/C matrices it lands in the pair-ANTIdiagonal
        # Case 4 rather than Case 1/3, so the density-density determination
        # above says nothing about it.
        #
        # Issue #100 settled it end to end (TestPairHopEndToEnd in
        # tests/test_vertex_orientation.py), with no index reasoning anywhere
        # between the input file and the answer: feed this builder an on-site
        # PairHop matrix, build S from it, form the solver's own linear-order
        # response -chi0 . S . chi0 with chi0 from `rpa.py`'s kernel and an
        # EXACT Green function, and compare against the cross-spin response of
        # the same Hamiltonian by exact diagonalization. The untransposed
        # placement -- i.e. S[(a,b),(b,a)] = P[a,b] -- reproduces it, and the
        # test locates the residual by checking chi0 against the exact bubble
        # in the same run, so what is left is the imaginary-time grid. Refining
        # that grid off-line, 192 -> 384 -> 768 points, takes the residual
        # 2.3e-3 -> 1.1e-3 -> 5.5e-4 (the test itself runs only the first);
        # transposing PairHop sits at 0.55 and does not move.
        #
        # Beware the plausible-looking derivation that says otherwise. The
        # vertex sits BETWEEN two chi0 factors, so its row index lives in
        # chi0's COLUMN convention and its column index in chi0's ROW
        # convention. chi0's row pair is reversed relative to its bilinear and
        # its column pair is not, so the reversal lands on the vertex's COLUMN
        # side. Applying it to the row side instead gives the transpose, and is
        # wrong.
        val_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for (irvec, orbvec), value in value_r.items():
            orb1, orb2 = orbvec
            Rx, Ry, Rz = irvec
            i1, i2 = (orb2, orb1) if transpose else (orb1, orb2)
            val_k[i1, i2, :, :, :] += value * np.exp(
                -1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
            )
        return val_k

    inter_k = {}
    for itype in ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                  "Ising", "PairLift", "PairHop"]:
        if itype in interactions:
            inter_k[itype] = _to_k(interactions[itype],
                                   transpose=(itype != "PairHop"))

    return inter_k


# ---------------------------------------------------------------------------
# Green's function
# ---------------------------------------------------------------------------

def _calc_eigenvalues(epsilon_k):
    """Diagonalize epsilon(k) at every k-point.

    Parameters
    ----------
    epsilon_k : ndarray
        Hamiltonian, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    """
    norb = epsilon_k.shape[0]
    Nx, Ny, Nz = epsilon_k.shape[2], epsilon_k.shape[3], epsilon_k.shape[4]
    eigenvalues = np.zeros((Nx, Ny, Nz, norb))
    eigenvectors = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)

    # Batch diagonalization: transpose to (Nx, Ny, Nz, norb, norb) for vectorized eigh
    eps_batch = epsilon_k.transpose(2, 3, 4, 0, 1)  # (Nx, Ny, Nz, norb, norb)
    eigenvalues, eigenvectors = np.linalg.eigh(eps_batch)

    return eigenvalues, eigenvectors


def _determine_mu(eigenvalues, beta, n_target, norb):
    """Determine chemical potential using bisection.

    The filling convention follows the reference code: n_target is the
    number of electrons per orbital per spin channel. For example,
    n_target=0.75 means 3/4 filling per spin.

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    beta : float
        Inverse temperature.
    n_target : float
        Target filling per orbital per spin (e.g. 0.75 for 3/4 filling).
    norb : int
        Number of orbitals.

    Returns
    -------
    mu : float
        Chemical potential.
    """
    Nx, Ny, Nz = eigenvalues.shape[:3]
    nvol = Nx * Ny * Nz

    def _calc_n(mu):
        x = beta * (eigenvalues - mu)
        fermi = np.where(x > 100, 0.0, np.where(x < -100, 1.0, 1.0 / (1.0 + np.exp(x))))
        total_n = np.sum(fermi)
        return float(total_n / nvol - n_target * norb)

    emin = np.min(eigenvalues)
    emax = np.max(eigenvalues)
    mu = bisect(_calc_n, emin - 10.0, emax + 10.0)
    return float(mu)


def _calc_green(eigenvalues, eigenvectors, mu, beta, nmat):
    """Construct Green's function G(k, iwn).

    Parameters
    ----------
    eigenvalues : ndarray
        Shape (Nx, Ny, Nz, norb).
    eigenvectors : ndarray
        Shape (Nx, Ny, Nz, norb, norb).
    mu : float
        Chemical potential.
    beta : float
        Inverse temperature.
    nmat : int
        Number of Matsubara frequencies.

    Returns
    -------
    green_kw : ndarray
        Shape (norb, norb, Nx, Ny, Nz, nmat).
    """
    Nx, Ny, Nz, norb = eigenvalues.shape
    iomega = np.array([(2.0 * i + 1.0 - nmat) * np.pi for i in range(nmat)]) / beta

    # Vectorized Green's function construction:
    # G_{ij}(k, iwn) = sum_m U_{im}(k) U*_{jm}(k) / (iwn - (e_m(k) - mu))

    # factor[kx,ky,kz,i,j,m] = U[kx,ky,kz,i,m] * conj(U[kx,ky,kz,j,m])
    factor = np.einsum('...im,...jm->...ijm', eigenvectors, np.conj(eigenvectors))
    # factor shape: (Nx, Ny, Nz, norb, norb, norb)

    # denom[kx,ky,kz,m,w] = 1 / (iwn_w - (e_m(k) - mu))
    xi = eigenvalues - mu  # (Nx, Ny, Nz, norb)
    denom = 1.0 / (1j * iomega[None, None, None, None, :] - xi[:, :, :, :, None])
    # denom shape: (Nx, Ny, Nz, norb, nmat)

    # G[kx,ky,kz,i,j,w] = sum_m factor[...,i,j,m] * denom[...,m,w]
    # Batched GEMM over the flattened spatial axes (C-order, reshape-only):
    #   (nv, ij, m) @ (nv, m, w) -> (nv, ij, w). numpy.einsum does NOT lower
    #   the '...ijm,...mw->...ijw' form to BLAS GEMM, so reshape to matmul.
    nv = Nx * Ny * Nz
    G = factor.reshape(nv, norb * norb, norb) @ denom.reshape(nv, norb, nmat)
    green_kw_tmp = G.reshape(Nx, Ny, Nz, norb, norb, nmat)
    # shape: (Nx, Ny, Nz, norb, norb, nmat)

    # Transpose to output convention: (norb, norb, Nx, Ny, Nz, nmat)
    green_kw = green_kw_tmp.transpose(3, 4, 0, 1, 2, 5)

    return green_kw


# ---------------------------------------------------------------------------
# RPA vertices
# ---------------------------------------------------------------------------

def _symmetrise_interactions_k(inter_k):
    """Reduce each interaction to its physical symmetric coefficient.

    Thin delegation: the reduction and its full derivation -- the
    momentum-space form of the reversed-bond partner, the PairHop
    Hermitian rule, and the UHFk-relation adjudication note -- are
    single-sourced in :func:`hwave.solver.declarations.symmetrise_k`
    (#108). Idempotent, so it is safe that both the all-q and per-q
    builders apply it.
    """
    return symmetrise_k(inter_k)


def _accumulate_coeff(dst, coeff, value):
    """dst += coeff * value with IEEE-preserving forms.

    NumPy evaluates ``1.0 * (Inf+0j)`` as ``Inf+NaNj`` (a full complex
    multiply), while a direct add preserves the zero imaginary part, so
    the +-1 coefficients of the vertex table must use direct add/subtract
    to keep the pre-table behavior for non-finite couplings. Zero
    coefficients suppress the contribution entirely (0.0 * Inf is NaN).
    ``dst`` must be an ndarray or writeable view; the update is in place.
    """
    if coeff == 0.0:
        return
    if coeff == 1.0:
        dst += value
    elif coeff == -1.0:
        dst -= value
    else:
        dst += coeff * value


def _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz,
                             _presymmetrised=False):
    """Build spin (S) and charge (C) interaction matrices for all q-points at once.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    S_all : ndarray
        Spin interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    C_all : ndarray
        Charge interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    """
    nd = norb * norb
    S_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
    C_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

    # _presymmetrised is set by the per-q wrapper, which has already
    # symmetrised on the FULL grid: the -q partner of an off-site entry is
    # unreachable from a single-point slice, so re-symmetrising here would
    # average with the wrong (same-q) partner and corrupt off-site input.
    if not _presymmetrised:
        inter_k = _symmetrise_interactions_k(inter_k)

    def _get(itype):
        if itype in inter_k:
            return inter_k[itype]  # (norb, norb, Nx, Ny, Nz)
        return None

    U_mat = _get("CoulombIntra")
    Up_mat = _get("CoulombInter")
    J_mat = _get("Hund")
    Jp_mat = _get("Exchange")
    I_mat = _get("Ising")
    PH_mat = _get("PairHop")

    # Build using precomputed index arrays to avoid Python loops
    # for small norb (1-3), the loop overhead is negligible;
    # for larger norb the vectorized approach helps
    l1_arr, l2_arr, l3_arr, l4_arr = np.meshgrid(
        np.arange(norb), np.arange(norb), np.arange(norb), np.arange(norb),
        indexing='ij')
    l1f, l2f, l3f, l4f = l1_arr.ravel(), l2_arr.ravel(), l3_arr.ravel(), l4_arr.ravel()
    idx12 = l1f * norb + l2f
    idx34 = l3f * norb + l4f

    # Case 1: l1 == l2 == l3 == l4
    # Accumulate rather than assign: this element is also reached by Case 3
    # below (which now includes l1 == l3), where an inter-site same-orbital
    # CoulombInter contributes 2 V_aa(q) to the charge channel.
    mask1 = (l1f == l2f) & (l2f == l3f) & (l3f == l4f)
    if U_mat is not None and np.any(mask1):
        sU, cU = sc_coefficients("CoulombIntra", "diag")
        for i in np.where(mask1)[0]:
            _l = l1f[i]
            _accumulate_coeff(S_all[:, :, :, idx12[i], idx34[i]],
                              sU, U_mat[_l, _l])
            _accumulate_coeff(C_all[:, :, :, idx12[i], idx34[i]],
                              cU, U_mat[_l, _l])

    # Case 2: l1==l3, l2==l4, l1!=l2
    # Coefficients come from the single adjudicated table
    # (hwave.solver.vertex_table). Signs, slots and per-type factors were
    # established against exact diagonalization in #113 (Ising sign, the
    # previously missing Hund S term, Exchange moved here from the
    # pair-antidiagonal Case 4); the SU(2) Kanamori consistency check --
    # Hund + Exchange at equal J giving S(ab,ab) without J and
    # C(ab,ab) = -U' + 2J -- is pinned in the tests.
    mask2 = (l1f == l3f) & (l2f == l4f) & (l1f != l2f)
    cross_terms = [(Up_mat, "CoulombInter"), (I_mat, "Ising"),
                   (J_mat, "Hund"), (Jp_mat, "Exchange")]
    for i in np.where(mask2)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        for mat, itype in cross_terms:
            if mat is not None:
                sco, cco = sc_coefficients(itype, "cross")
                _accumulate_coeff(s_q, sco, mat[_l1, _l2])
                _accumulate_coeff(c_q, cco, mat[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 3: l1==l2, l3==l4 -- INCLUDING l1 == l3, but for CoulombInter ONLY.
    # The l1 != l3 exclusion dropped the inter-site same-orbital CoulombInter
    # from the charge diagonal C[(a,a),(a,a)], which must be U_a + 2 V_aa(q)
    # (issue #95); the simple two-index formulation used by chi0q_mode="load"
    # builds exactly that (`Wc = U_k + 2 V_k`, _compute_vertices_simple).
    # Case 1 above writes U_a into the same element, so both accumulate.
    # Hund and Ising stay restricted to l1 != l3: an orbital has no Hund or
    # Ising coupling with itself, and letting a stray diagonal entry through
    # here would silently move S as well.
    mask3 = (l1f == l2f) & (l3f == l4f)
    for i in np.where(mask3)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l3 = l1f[i], l3f[i]
        if _l1 != _l3:
            for mat, itype in ((J_mat, "Hund"), (I_mat, "Ising")):
                if mat is not None:
                    sco, cco = sc_coefficients(itype, "density")
                    _accumulate_coeff(s_q, sco, mat[_l1, _l3])
                    _accumulate_coeff(c_q, cco, mat[_l1, _l3])
        if Up_mat is not None:
            sco, cco = sc_coefficients("CoulombInter", "density")
            _accumulate_coeff(s_q, sco, Up_mat[_l1, _l3])
            _accumulate_coeff(c_q, cco, Up_mat[_l1, _l3])
        S_all[:, :, :, idx12[i], idx34[i]] += s_q
        C_all[:, :, :, idx12[i], idx34[i]] += c_q

    # Case 4: l1==l4, l2==l3, l1!=l2
    mask4 = (l1f == l4f) & (l2f == l3f) & (l1f != l2f)
    for i in np.where(mask4)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        # Exchange used to sit here; exact diagonalization (issue #113) puts
        # its vertex on the pair-DIAGONAL slot family (Case 2), and end to end
        # the antidiagonal placement produced the right magnitude at the wrong
        # slots in both channels. Only PairHop belongs here (#100/#102).
        if PH_mat is not None:
            sco, cco = sc_coefficients("PairHop", "antidiag")
            _accumulate_coeff(s_q, sco, PH_mat[_l1, _l2])
            _accumulate_coeff(c_q, cco, PH_mat[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    return S_all, C_all


def _build_sc_matrices(inter_k, norb, ix, iy, iz):
    """Spin (S) and charge (C) matrices at a single q-point.

    Delegates to :func:`_build_sc_matrices_all_q` and slices, so there is ONE
    implementation of the S/C content. The previous hand-maintained copy had
    already drifted from the all-q builder once; after the issue #113 vertex
    corrections a second parallel copy would be a liability.
    """
    # Symmetrise on the FULL grid first -- the same-operator partner of an
    # off-site entry lives at -q, so slicing before symmetrising would discard
    # it -- then slice the (small, norb^2 per point) interaction arrays down to
    # the requested q and delegate on a 1x1x1 grid, so the (nd^2 per point) S/C
    # matrices are never built for the whole grid. Indexing goes through
    # np.arange so numpy negative indices keep their usual meaning and
    # out-of-range indices raise instead of returning empty slices.
    inter_sym = _symmetrise_interactions_k(inter_k)
    inter_1 = {}
    for t, M in inter_sym.items():
        jx = np.arange(M.shape[2])[ix]
        jy = np.arange(M.shape[3])[iy]
        jz = np.arange(M.shape[4])[iz]
        inter_1[t] = np.ascontiguousarray(
            M[:, :, jx:jx+1, jy:jy+1, jz:jz+1])
    S_all, C_all = _build_sc_matrices_all_q(inter_1, norb, 1, 1, 1,
                                            _presymmetrised=True)
    return S_all[0, 0, 0], C_all[0, 0, 0]

def _declarations_partner_closed(interactions, Nx, Ny, Nz, norb):
    """True iff every CoulombIntra/CoulombInter declaration table is
    EXACTLY closed under the reversed-bond partner (R,a,b) <-> (-R,b,a).

    Decided on the RAW tables, before any exponential transform, with
    bitwise comparison: a both-ends declaration carries identical
    literals, and 0.5*(x + x) == x exactly in IEEE, so closure detection
    is exact -- while the k-space arrays of a closed table are only
    ALGEBRAICALLY symmetric (the +R and -R exponentials are evaluated
    independently and differ by roundoff), which is why the decision
    cannot be made after the transform (PR #129 round 5: re-averaging a
    closed configuration changed saved artifacts at the 1e-18 level).
    Non-finite input reads as not closed (explicitly, before any
    comparison), which fails toward symmetrising. The closure test
    compares the table DIRECTLY with its reversed/transposed partner --
    no mean is formed, so a closed coefficient above float64.max / 2
    cannot overflow into a false 'open' (PR #129 round 6 reproduced
    exactly that with the earlier mean-based comparison)."""
    for itype in ("CoulombIntra", "CoulombInter"):
        tbl = interactions.get(itype) if interactions else None
        if not tbl:
            continue
        arr = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)
        for (irvec, orbvec), v in tbl.items():
            arr[(irvec[0] % Nx, irvec[1] % Ny, irvec[2] % Nz,
                 *orbvec)] += v
        if not np.all(np.isfinite(arr)):
            return False
        rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                           (0, 1, 2, 4, 3))
        if not (np.array_equal(arr.real, rev.real)
                and np.array_equal(arr.imag, rev.imag)):
            return False
    return True


def _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                      pairing_type="singlet", static_index=None,
                      declarations_closed=False):
    """Compute effective pairing interaction V(q).

    Supports two modes:
    - Simple mode (only CoulombIntra/Inter, 2-index chi0q):
      Uses Wc=U+2V, Ws=-U formulation.
    - General mode (with Hund/Exchange, or 4-index chi0q):
      Uses S,C matrices from Kuroki et al.

    When 4-index chi0q is provided (8D array), general mode is always used
    because the full orbital tensor structure must be preserved.

    Pairing types:
    - singlet: V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)
    - triplet: V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Parameters
    ----------
    chi0q : ndarray
        Bare susceptibility.
        - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb, Nx, Ny, Nz, nmat : int
        System parameters.
    pairing_type : str
        "singlet" or "triplet". Default "singlet".

    Returns
    -------
    Result depends on mode:
        Simple mode: tuple (Pc_q, Ps_q), each shape (norb, norb, Nx, Ny, Nz).
        General mode: Vs_q, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    has_interorbital_vertex = any(k in inter_k for k in
                                  ["Hund", "Exchange", "Ising", "PairHop"])
    chi0q_is_4index = (chi0q.ndim == 8)

    # PairLift does not enter the spin/charge (particle-hole) pairing vertex in
    # either mode (its S=C=0 contribution, verified against the full 4-index
    # RPA), so a configured PairLift term is silently inert. Warn the user.
    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    if chi0q_is_4index or has_interorbital_vertex:
        # General mode: 4-index S,C matrices
        # Required when 4-index chi0q is available or Hund/Exchange present
        return _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                         pairing_type=pairing_type,
                                         static_index=static_index)
    else:
        # Simple mode: backward compatible with original implementation
        # Only used for 2-index chi0q without Hund/Exchange
        if "CoulombInter" in inter_k and norb > 1:
            logger.warning(
                "Multi-orbital CoulombInter in the simple (2-index) vertex mode "
                "uses the reduced density-density approximation: inter-orbital "
                "cross-channel contributions are dropped. Provide a 4-index "
                "chi0q (general mode) for the full Kuroki S/C treatment.")
        return _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type=pairing_type,
                                        static_index=static_index,
                                        declarations_closed=declarations_closed)


def _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                             pairing_type="singlet", static_index=None,
                             declarations_closed=False):
    """Compute vertices using simple Wc=U+2V, Ws=-U formulation.

    For singlet:
        Pc = (Wc+Ws)/2 - (1/2) Wc chi_c Wc
        Ps = -Ws + (3/2) Ws chi_s Ws
    For triplet:
        Pc = (Wc-Ws)/2 - (1/2) Wc chi_c Wc   (charge part same sign change)
        Ps = Ws - (1/2) Ws chi_s Ws            (spin part sign flip)

    Returns Pc_q and Ps_q separately for backward compatibility.

    Returns
    -------
    Pc_q : ndarray
        Charge vertex, shape (norb, norb, Nx, Ny, Nz).
    Ps_q : ndarray
        Spin vertex, shape (norb, norb, Nx, Ny, Nz).
    """
    # Symmetrise FIRST (PR #129 round 3): this path read the raw tables,
    # so a one-sided off-site declaration entered as v e^{-iqR} while the
    # ring and the general S/C route read the same Hamiltonian as
    # v cos(qR) -- measured drift 0.7 at q = pi/2 for V(R=+x) = 0.7.
    # SKIPPED when the caller proved the raw declarations partner-closed
    # (round 5): the reduction is only ALGEBRAICALLY idempotent after the
    # exponential transform -- re-averaging a closed configuration mixed
    # independently rounded +R/-R exponentials and changed saved
    # artifacts at the 1e-18 level, violating byte parity for
    # previously-working symmetric runs.
    if not declarations_closed:
        inter_k = _symmetrise_interactions_k(
            {k: inter_k[k] for k in ("CoulombIntra", "CoulombInter")
             if k in inter_k})
    U_k = inter_k.get("CoulombIntra", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))
    V_k = inter_k.get("CoulombInter", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))

    # Transpose to batch dimension first: (Nx, Ny, Nz, norb, norb)
    Wc = (U_k + 2.0 * V_k).transpose(2, 3, 4, 0, 1).copy()
    Ws = (-U_k).transpose(2, 3, 4, 0, 1).copy()

    # chi0 at static limit: (Nx, Ny, Nz, norb, norb)
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    chi0_static = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()

    I_mat = np.broadcast_to(np.eye(norb, dtype=complex), (Nx, Ny, Nz, norb, norb)).copy()

    # Batched solve
    mat_s = I_mat + chi0_static @ Ws
    mat_c = I_mat + chi0_static @ Wc

    chis = np.linalg.solve(mat_s, chi0_static)
    chic = np.linalg.solve(mat_c, chi0_static)

    WsChisWs = Ws @ chis @ Ws
    WcChicWc = Wc @ chic @ Wc

    if pairing_type == "singlet":
        Pc_all = (Wc + Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = -Ws + 1.5 * WsChisWs
    elif pairing_type == "triplet":
        Pc_all = (Wc - Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = Ws - 0.5 * WsChisWs
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Transpose back: (Nx, Ny, Nz, norb, norb) -> (norb, norb, Nx, Ny, Nz)
    Pc_q = Pc_all.transpose(3, 4, 0, 1, 2)
    Ps_q = Ps_all.transpose(3, 4, 0, 1, 2)

    return Pc_q, Ps_q


def _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                              pairing_type="singlet", static_index=None):
    """Compute effective pairing interaction using general S,C matrices.

    For singlet (Kuroki et al., arXiv:0902.3691, Eq.(6)):
        V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)

    For triplet (Takimoto et al., PRB 69, 104504):
        V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Supports both 2-index (reduced) and 4-index (general) chi0q:
    - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
      Expanded to diagonal of (norb^2, norb^2) matrix.
    - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
      Used directly as (norb^2, norb^2) matrix. This is the correct
      treatment for multi-orbital systems.

    Returns
    -------
    V_q : ndarray
        Effective pairing interaction, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb
    chi0q_is_4index = (chi0q.ndim == 8)

    # Build S, C matrices for all q-points at once: (Nx, Ny, Nz, nd, nd)
    S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)

    # Extract chi0 at static limit for all q-points
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    if chi0q_is_4index:
        # (norb, norb, norb, norb, Nx, Ny, Nz, nmat) -> (Nx, Ny, Nz, nd, nd)
        chi0_static = chi0q[:, :, :, :, :, :, :, si].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()
    else:
        # (norb, norb, Nx, Ny, Nz, nmat) -> expand to (Nx, Ny, Nz, nd, nd)
        chi0_2d = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()
        # chi0_2d shape: (Nx, Ny, Nz, norb, norb)
        # The guard runs for EVERY 2-index input, including norb == 1: the
        # rejection of Exchange/PairHop applies there too (the one-orbital
        # builder gives them zero weight, so accepting them silently omits
        # the interaction -- the norb shortcut below bypassed the helper
        # entirely, round 8). Only the partial-dressing warning is
        # norb-gated, inside the helper.
        _warn_reduced_flex_missing_components(
            inter_k, norb, Nx, Ny, Nz, source="a reduced 2-index chi0q",
            sc_matrices=(S_all, C_all))
        if norb == 1:
            chi0_static = chi0_2d.reshape(Nx, Ny, Nz, 1, 1)
        else:
            # A 2-index (reduced/squashed) chi0q is the density-density diagonal
            # of the susceptibility: chi0_2d[a, b] IS chi0_{(a,a),(b,b)} (the
            # matching interaction reduction is einsum('kaabb->kab', ...) in
            # RPA._inflate_chi0q_and_ham).  With the orbital-pair flat index
            # (l1,l2) -> l1*norb + l2 used by _build_sc_matrices_all_q, it
            # therefore belongs at the density-pair positions:
            #
            #     chi0_{(a,a),(b,b)} = chi0_2d[a, b],  everything else zero.
            #
            # The historical placement chi0_{(l1,l2),(l3,l2)} = chi0_2d[l1,l3]
            # (a delta_{l2,l4} scatter, i.e. kron(chi0_2d, I_norb)) read the
            # density-pair index as an ordinary orbital index: it dropped the
            # inter-orbital density coupling chi0_{(0,0),(1,1)} and scattered
            # chi0_2d onto pair indices the reduced scheme never computed.  This
            # is the same defect that _expand_flex_chi carried on the
            # chi0q_mode="flex" route.
            #
            # Off-density rows/columns stay exactly zero.  For the interaction
            # terms that put S/C weight there (CoulombInter, Hund, Ising,
            # Exchange, PairHop) those channels are then undressed -- an honest
            # reflection of what a reduced chi0q contains, rather than a
            # fabricated dressing.
            chi0_static = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            dens = np.arange(norb) * norb + np.arange(norb)
            chi0_static[..., dens[:, None], dens[None, :]] = chi0_2d

    # Batched RPA solve for all q-points simultaneously
    # chi_s = [I - chi0 @ S]^{-1} @ chi0
    # chi_c = [I + chi0 @ C]^{-1} @ chi0
    I_mat = np.broadcast_to(np.eye(nd, dtype=complex), (Nx, Ny, Nz, nd, nd)).copy()

    mat_s = I_mat - chi0_static @ S_all
    mat_c = I_mat + chi0_static @ C_all

    chis = np.linalg.solve(mat_s, chi0_static)  # batched solve
    chic = np.linalg.solve(mat_c, chi0_static)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + 0.5 * (C_all - S_all)
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Reshape (Nx, Ny, Nz, nd, nd) -> (norb, norb, norb, norb, Nx, Ny, Nz)
    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(3, 4, 5, 6, 0, 1, 2)

    return Vs_q


#: Interaction terms whose Kuroki S/C matrices have entries OUTSIDE the
#: density-pair block, keyed by the _build_sc_matrices_all_q case that puts
#: them there.  Case 2 is S/C[(a,b),(a,b)] and case 4 is S/C[(a,b),(b,a)], both
#: with a != b; a reduced/squashed run never computes chi on those pair indices,
#: so the corresponding fluctuation dressing is missing from the vertex.
#: Types whose vertex is PARTIALLY represented by a reduced chi: their
#: density-slot content is dressed, their cross-slot content is not.
_REDUCED_FLEX_PARTIAL = ("CoulombInter", "Hund", "Ising")
#: Types with NO density-diagonal vertex content at all
#: (hwave.solver.vertex_table): with a reduced chi nothing of them is
#: dressed, and since the #120 policy no reduced/squashed FLEX or RPA run
#: can even be produced with them -- such input is stale or mismatched,
#: and is REJECTED rather than warned about.
_REDUCED_FLEX_REJECTED = ("Exchange", "PairHop")



def _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz, sc_matrices=None):
    """Largest |S| or |C| on the blocks a reduced chi cannot dress.

    The reduced susceptibility is zero on every pair index (a,b) with a != b, so
    the vertex is undressed exactly where S or C is nonzero there --
    S/C[(a,b),(a,b)] (case 2) and S/C[(a,b),(b,a)] (case 4) of
    :func:`_build_sc_matrices_all_q`.

    This inspects the COMBINED matrices rather than testing each configured term
    on its own. Those blocks are sums: under the adjudicated slot table
    (#113) case 2 mixes CoulombInter, Ising, Hund AND Exchange, and case 4
    carries PairHop alone. Terms can cancel there (equal CoulombInter and
    Hund cancel case 2 exactly in both channels), and a per-term test would
    then announce missing dressing that does not exist.
    """
    # Reuse the caller's matrices when it already has them: at the full grid
    # these are O(Nq * norb^4) each, and building a second pair while the first
    # is live doubles the allocation for what is only a diagnostic.
    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    nd = norb * norb
    density = np.zeros(nd, dtype=bool)
    density[np.arange(norb) * norb + np.arange(norb)] = True
    off = ~density
    if not off.any():
        return 0.0
    weight = 0.0
    for M in (S_all, C_all):
        # anything outside the density x density sub-block, over all q
        weight = max(weight,
                     float(np.max(np.abs(M[..., off, :]))),
                     float(np.max(np.abs(M[..., :, off]))))
    return weight


def _build_vertex_sc_matrices(convention, inter_k, norb, Nx, Ny, Nz):
    """Build the (S, C) interaction matrices for the given orbital convention.

    Single dispatcher shared by :func:`_compute_vertices_flex` and its
    callers, so a dynamic run can build the frequency-INDEPENDENT matrices
    once and pass them to every per-frequency contraction instead of
    rebuilding ~nmat identical full-grid pairs.
    """
    conv = str(convention).lower()
    if conv == "myo":
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        return build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz)
    if conv == "kuroki":
        return _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    raise ValueError(
        "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(convention))


def _reject_reduced_flex_unsupported(inter_k, convention="kuroki",
                                     source=None):
    """Reject Exchange/PairHop wherever a Kuroki-convention vertex is built.

    This is the cheap half of :func:`_warn_reduced_flex_missing_components`:
    a dict scan with no S/C construction, so it can (and must) run at EVERY
    Kuroki vertex boundary -- :func:`_compute_vertices_flex` itself, hence
    also ``eliashberg_dynamic._instantaneous_vertex`` and every
    per-frequency dynamic contraction -- making the rejection an enforced
    invariant rather than an ordering artifact of which builder the
    orchestrator happens to call first (round-9 review: the IR
    instantaneous route silently accepted both terms and returned a zero
    vertex when called directly).

    Exchange and PairHop carry NO density-diagonal vertex content
    (hwave.solver.vertex_table), so a Kuroki (reduced-origin) chi dresses
    none of it and the result would silently omit the interaction; the
    general (myo) convention represents them fully and is not restricted.
    """
    if str(convention).lower() != "kuroki":
        return
    rejected = [k for k in _REDUCED_FLEX_REJECTED
                if k in inter_k and np.any(np.asarray(inter_k[k]) != 0)]
    if rejected:
        raise ValueError(
            "the Eliashberg vertex cannot be built from {} together with "
            "{}: those interactions have no density-diagonal vertex "
            "content at all (hwave.solver.vertex_table), so reduced "
            "(density-only) data dresses none of it and the result would "
            "silently omit the interaction. (H-wave's own reduced/squashed "
            "runs cannot even be produced with these terms since the "
            "unified scheme policy.) Provide a general (four-index) "
            "susceptibility instead -- for a FLEX source, re-run with "
            "calc_scheme='general'.".format(
                source or ("a REDUCED (calc_scheme='reduced' or "
                           "'squashed') FLEX susceptibility"),
                ", ".join(rejected)))


def _warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz,
                                          convention="kuroki",
                                          source=None, sc_matrices=None):
    """Warn when a reduced (kuroki) FLEX chi cannot support the interaction.

    Call this ONCE per run, from the place that is about to build the pairing
    vertex -- NOT from inside ``_compute_vertices_flex``. That function is
    invoked once per bosonic Matsubara frequency by the dynamic kernel (so the
    warning would repeat ``Nmat`` times, ~1000 in production runs), and again by
    ``eliashberg_dynamic._instantaneous_vertex`` with chi = 0, where the
    message would be doubly misleading. The Exchange/PairHop REJECTION, by
    contrast, is enforced inside ``_compute_vertices_flex`` on every call
    (:func:`_reject_reduced_flex_unsupported`) -- it is a cheap scan and
    must hold on every route, not just the orchestrated one.

    A ``calc_scheme="reduced"``/``"squashed"`` FLEX run stores only the
    density-density diagonal chi_{(a,a),(b,b)} of the susceptibility.  The
    Kuroki S/C matrices built from CoulombIntra alone live entirely on that
    density-pair block, so the reduced treatment is exact.  CoulombInter,
    Hund and Ising also populate the off-density block S/C[(a,b),(a,b)]
    with a != b -- and there chi is identically zero simply because the
    reduced run never computed it.  Those channels then keep only the bare
    0.5*(S+C) term with no fluctuation dressing: a silent approximation
    rather than a solver error, hence the WARNING.  Exchange and PairHop
    carry NO density-diagonal vertex content at all (vertex_table), so a
    reduced chi dresses nothing of them; combined with the #120 scheme
    policy (no reduced/squashed run can be produced with them) that
    combination is REJECTED rather than warned about.

    This is a genuine limitation of the stored data, not of the loader: it
    cannot be repaired on the Eliashberg side.  The universal remedy is a
    general (four-index) susceptibility; for a FLEX source that means
    re-running with ``calc_scheme="general"`` (which stores the full
    orbital-pair chi).
    """
    if str(convention).lower() != "kuroki":
        # Only the reduced/squashed route stores a density-only chi; the
        # general (myo) path carries the full orbital-pair susceptibility.
        return
    # The rejection must run BEFORE the norb == 1 shortcut: a one-orbital
    # Exchange/PairHop has zero S/C weight in this builder too, so the
    # interaction would be silently omitted rather than represented --
    # exactly what the rejection exists to prevent (round-7 review; the
    # bypass returned zero vertices for both terms on both routes).
    _reject_reduced_flex_unsupported(inter_k, convention, source)
    if norb <= 1:
        # norb == 1 has no off-density pair index, so no PARTIAL dressing
        # can be missing (the rejection above still applies).
        return
    configured = [k for k in _REDUCED_FLEX_PARTIAL if k in inter_k]
    if not configured:
        return
    # Ask the assembled S/C matrices, not the individual terms: the off-density
    # blocks are sums and the configured terms can cancel there (e.g. equal
    # CoulombInter and Hund cancel case 2 in both channels under the
    # adjudicated table). Only a nonzero combined block means dressing is
    # actually missing.
    if _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz,
                              sc_matrices) == 0.0:
        return
    # Decision above is on the assembled matrices; attribution here is per term,
    # so an inert term that happens to be configured is not named as a cause.
    # Attribute through the same assembled-and-symmetrised path as the
    # decision: a raw declaration that symmetrises to zero (#112) must not
    # be named as a cause.
    missing = [k for k in configured
               if _off_density_sc_weight({k: inter_k[k]}, norb,
                                         Nx, Ny, Nz) != 0.0] or configured
    logger.warning(
        "The Eliashberg vertex is being built from %s, which carries only the "
        "density-density components chi_{(a,a),(b,b)}, together with "
        "inter-orbital interaction(s) %s. Together those terms give the S/C "
        "matrices nonzero "
        "off-density components S/C[(a,b),(a,b)] and S/C[(a,b),(b,a)] (a != b) "
        "for which the reduced run computed no susceptibility, so those "
        "channels enter the pairing vertex undressed (bare 0.5*(S+C) only) and "
        "the resulting lambda is an approximation. calc_scheme='general' stores "
        "the full orbital-pair chi and removes this gap; note its off-site "
        "support is limited to same-orbital CoulombInter without sublattice "
        "folding -- for a model with off-site INTER-orbital interactions "
        "there is currently no FLEX-dressed vertex without this "
        "approximation.",
        source or ("a REDUCED (calc_scheme='reduced' or 'squashed') FLEX "
                   "susceptibility"),
        ", ".join(missing))


def _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                           pairing_type="singlet", convention="kuroki",
                           sc_matrices=None):
    """Compute pairing vertex from pre-computed FLEX susceptibilities.

    Uses the same formula as _compute_vertices_general but takes
    chi_s and chi_c directly instead of computing them from chi0q via RPA.
    This allows using dressed susceptibilities from the FLEX solver.

    ``convention`` records which FLEX path produced the susceptibilities:
    "kuroki" (default; reduced path) or "myo" (general full-vertex path). The
    two S/C builders historically differed in the charge ``C(ab,ab)`` element
    (``-U'+2J`` vs ``-U'+J``); the exact-diagonalization adjudication of the
    per-type vertex content (issue #113: Hund ``+J`` and Exchange ``+J'``
    there) made them IDENTICAL, so the flag no longer selects different
    matrices -- it remains a provenance / layout discriminator, and legacy
    files are guarded by ``sc_vertex_version`` instead. ``chis``/``chic`` are
    expected in the native [a,c,b,d] orbital-pair order regardless
    (issue #78).

    Parameters
    ----------
    chis : ndarray
        Spin susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd)
        where nd = norb^2.
    chic : ndarray
        Charge susceptibility at static limit, same shape.
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.
    pairing_type : str
        "singlet" or "triplet".
    sc_matrices : tuple of ndarray, optional
        Precomputed ``(S_all, C_all)`` for this convention
        (:func:`_build_vertex_sc_matrices`). The matrices are
        frequency-independent, so a dynamic run passes one pair to every
        per-frequency call instead of rebuilding ~nmat identical full-grid
        pairs (round-9 review).

    Returns
    -------
    Vs_q : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb

    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    # Validate the convention BEFORE the cache branch: with precomputed
    # sc_matrices the dispatcher below is skipped, and an unknown tag would
    # otherwise be silently treated as "not kuroki" = unrestricted
    # (round-10 review reproduced convention="invalid" returning a vertex).
    conv = str(convention).lower()
    if conv not in ("kuroki", "myo"):
        raise ValueError(
            "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(
                convention))

    # Enforced at THIS boundary (not only in the orchestrator): every Kuroki
    # vertex construction -- including chi = 0 via _instantaneous_vertex --
    # must reject Exchange/PairHop, or the interaction is silently omitted.
    _reject_reduced_flex_unsupported(inter_k, conv)

    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_vertex_sc_matrices(conv, inter_k,
                                                 norb, Nx, Ny, Nz)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + 0.5 * (C_all - S_all)
    else:
        raise ValueError("Unknown pairing_type: '{}'".format(pairing_type))

    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
        3, 4, 5, 6, 0, 1, 2)

    return Vs_q


def _resolve_flex_paths(input_dict):
    """Resolve the FLEX output directory and chi_s/chi_c/green file paths.

    Shared by `_load_flex_susceptibilities_full` (which reads the files) and
    `_load_flex_susceptibilities` (which re-derives the static frequency
    index from the same files' metadata), so the two never disagree about
    which files are in play.
    """
    file_input = input_dict.get("file", {}).get("input", {})
    flex_dir = file_input.get("path_to_flex_output",
                              input_dict.get("file", {}).get("output", {}).get(
                                  "path_to_output", "output"))

    eli_param = input_dict.get("eliashberg", {})
    chi_s_file = eli_param.get("flex_chi_s", "chiq_s.npz")
    chi_c_file = eli_param.get("flex_chi_c", "chiq_c.npz")
    green_file = eli_param.get("flex_green", "green.npz")

    chi_s_path = os.path.join(flex_dir, chi_s_file)
    chi_c_path = os.path.join(flex_dir, chi_c_file)
    green_path = os.path.join(flex_dir, green_file)
    return chi_s_path, chi_c_path, green_path


def _load_flex_susceptibilities_full(input_dict, norb, Nx, Ny, Nz,
                                     allow_ir=False, interactions=None):
    """Load FLEX-computed susceptibilities from NPZ files (full frequency axis).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    chis_w : ndarray
        Spin susceptibility, shape (Nx, Ny, Nz, nd, nd, nmat) with nd = norb^2
        (after the spin-orbital block-extraction/expansion, if needed) and
        the full bosonic Matsubara axis as the last dimension.
    chic_w : ndarray
        Charge susceptibility, same shape/convention as `chis_w`.
    green_w : ndarray or None
        Dressed Green's function if available, shape
        (norb, norb, Nx, Ny, Nz, nmat) with the full fermionic axis.
    chi_convention : str
        Provenance/S-C-convention tag of chis_w/chic_w: "myo" (general
        full-vertex FLEX) or "kuroki" (reduced FLEX / legacy files). It selects
        which S/C interaction matrices _compute_vertices_flex pairs the
        susceptibilities with (one implementation since the #113
        adjudication; the tag is a layout/provenance discriminator)
        and which orbital-pair shape family the file uses (orbital-pair
        nd=norb^2 for "myo", spin-orbital nd=norb*ns for "kuroki"). It is NOT an
        index-layout flag: chis_w/chic_w are returned in the public RPA
        [a,c,b,d] orbital-pair order for both (issue #78).
    ir_meta : dict or None
        Only when ``allow_ir=True``: ``None`` for uniform files, else the
        per-file IR node metadata ``{"chis":..., "chic":..., "green":...}``
        ("green" absent when there is no green file). Mixed encodings raise.
    """
    if not allow_ir:
        (chi_s_raw, chi_c_raw, chi_convention,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, interactions=interactions)
        ir_meta = None
        green_w = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    else:
        (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, allow_ir=True, interactions=interactions)
        green_w, green_meta = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                               allow_ir=True)
        if green_w is not None and (green_meta is None) != (ir_meta is None):
            raise ValueError(
                "all dynamic-Eliashberg inputs must share one encoding; "
                "re-run FLEX so chiq_s/chiq_c/green are all densified or "
                "all IR-native (chi: {}, green: {}).".format(
                    "IR-native" if ir_meta is not None else "densified",
                    "IR-native" if green_meta is not None else "densified"))
        if ir_meta is not None and green_meta is not None:
            ir_meta = dict(ir_meta, green=green_meta)

    # Expand the FULL frequency axis (the static slice is selected by the
    # caller). The frequency axis is moved from leading to trailing position.
    legacy_s, legacy_c = legacy_tags
    if _accept_up_block_only(input_dict):
        legacy_s = legacy_c = "config_override"
    _check_spin_block_discarded(chi_s_raw, norb, chi_convention, "chi_s",
                                legacy_s)
    _check_spin_block_discarded(chi_c_raw, norb, chi_convention, "chi_c",
                                legacy_c)
    chis_w = np.moveaxis(
        _expand_flex_chi(chi_s_raw, norb, Nx, Ny, Nz, chi_convention), 0, -1)
    chic_w = np.moveaxis(
        _expand_flex_chi(chi_c_raw, norb, Nx, Ny, Nz, chi_convention), 0, -1)

    if not allow_ir:
        return chis_w, chic_w, green_w, chi_convention
    return chis_w, chic_w, green_w, chi_convention, ir_meta


_STATIC_IR_HINT = (
    "the static Eliashberg solver requires uniform-grid files. Re-run FLEX "
    "with [mode.param] write_densified = true, or switch to frequency = "
    "\"dynamic\" with [eliashberg] matsubara_basis = \"ir\".")


def _accept_up_block_only(input_dict):
    """Whether ``[eliashberg] accept_up_block_only`` asserts the layout.

    The file tag is the better signal, but files written before this check
    existed cannot carry one -- so without a configuration route the escape
    hatch would be unreachable for exactly the files that need it. The user
    takes responsibility here rather than the loader guessing from values it
    cannot disambiguate.
    """
    eli = input_dict.get("eliashberg", {})
    value = eli.get("accept_up_block_only", False)
    if isinstance(value, bool):
        return value
    # Not bool(value): a TOML typo like the string "false" is truthy in Python,
    # and silently enabling an override that relaxes a correctness guard is the
    # worst possible reading of a malformed value.
    raise ValueError(
        "[eliashberg] accept_up_block_only must be a boolean (true/false), "
        "got {!r}.".format(value))


def _legacy_up_block_only(data):
    """Whether an npz declares the legacy single-populated-block layout.

    Opt-in only: the values cannot distinguish "the other blocks were never
    filled" from "this run is fully polarized", so the file has to say so.
    """
    if "chi_spin_blocks" not in data:
        return False
    return str(data["chi_spin_blocks"]).strip().lower() == "up_only"


def _onsite_transpose_asymmetry(table):
    """Max |X[0,(a,b)] - X[0,(b,a)]| over the on-site entries of a raw
    declaration table. Off-site entries are ignored on purpose: their
    same-operator partner lives at -R (and possibly under a wrapped
    canonical key), and the general FLEX path rejects off-site two-body
    terms anyway -- the reachable orientation-sensitive case (#101) is
    an asymmetric ON-SITE inter-orbital entry."""
    worst = 0.0
    for (irvec, orbvec), v in table.items():
        if tuple(irvec) != (0, 0, 0):
            continue
        a, b = orbvec
        if a == b:
            continue
        partner = table.get((tuple(irvec), (b, a)), 0.0)
        diff = abs(complex(v) - complex(partner))
        if not np.isfinite(diff):
            # fail CLOSED: max() drops NaN (nan > x is False), so a
            # non-finite declaration would otherwise read as symmetric
            return float("inf")
        worst = max(worst, diff)
    return worst


def _onsite_hermitian_pair_mismatch(table):
    """Max |X[0,(a,b)] - conj(X[0,(b,a)])| over on-site entries.

    PairHop's two declarations of one pair are HERMITIAN partners, so its
    orientation never depended on the #99 transpose -- but the version-2
    stamp introduced the conjugated-mean symmetrisation of those
    declarations, and a pair that is NOT Hermitian-closed was read
    differently by older builds (round-4 review measured |dS| = 0.15 in
    the myo matrices for a real asymmetric declaration). Fails closed on
    non-finite values like the transpose helper."""
    worst = 0.0
    for (irvec, orbvec), v in table.items():
        if tuple(irvec) != (0, 0, 0):
            continue
        a, b = orbvec
        if a == b:
            continue
        partner = table.get((tuple(irvec), (b, a)), 0.0)
        diff = abs(complex(v) - np.conj(complex(partner)))
        if not np.isfinite(diff):
            return float("inf")
        worst = max(worst, diff)
    return worst


def _read_flex_chi_raw(input_dict, allow_ir=False, interactions=None):
    """Read the raw FLEX chi_s / chi_c NPZ arrays and their orbital convention.

    Returns ``(chi_s_raw, chi_c_raw, chi_convention, legacy_tags)`` in the
    H-wave layout. ``legacy_tags`` is the pair of ``chi_spin_blocks`` flags read
    from the SAME open handles as the arrays: reopening the files afterwards
    would let a tagged replacement authorize the untagged -- possibly polarized
    -- array already in memory. No values from the replacement would be used,
    but the wrong ones would be accepted in silence.
    ``(nfreq, nvol, nd, nd)`` -- no reshape/expansion, so callers that only
    need one static frequency can slice before expanding.

    The tuple is ``(chi_s_raw, chi_c_raw, chi_convention, legacy_tags)``; with
    ``allow_ir=True`` it is ``(chi_s_raw, chi_c_raw, chi_convention, ir_meta,
    legacy_tags)``. In that case the return value
    gains a fourth element ``ir_meta``: ``None`` for uniform files, or
    ``{"chis": meta, "chic": meta}`` for IR-native ones (each file carries
    its own node set; the caller refits both independently). The arity
    changes ONLY under ``allow_ir=True``, so every existing call site keeps
    its unpacking. Mixed encodings (one native, one densified) are rejected
    -- refitting one of the pair would silently pair susceptibilities of
    different provenance."""
    chi_s_path, chi_c_path, _ = _resolve_flex_paths(input_dict)

    logger.info("Loading FLEX chi_s from: {}".format(chi_s_path))
    data_s = np.load(chi_s_path)
    if not allow_ir:
        _reject_ir_native(data_s, chi_s_path, _STATIC_IR_HINT)
    chi_s_raw = data_s["chiq_s"] if "chiq_s" in data_s else data_s["chiq"]
    # Orbital convention tag (general FLEX writes "myo", reduced "kuroki").
    # The tag and the general-path MYO consumption ship together, so any
    # untagged file from a released build is necessarily a reduced/Kuroki
    # output -> default to "kuroki". (A pre-tag general output could only exist
    # as a transient artifact of an unreleased dev build; the s/c-agreement
    # check below still guards against accidentally mixing conventions.)
    # This default is legacy-only: an IR-native file missing the tag is
    # rejected below instead of silently defaulting (an IR-native norb=2 MYO
    # file would otherwise be mis-read as spin-orbital "kuroki").
    chi_convention_present_s = "chi_convention" in data_s
    chi_convention = (str(data_s["chi_convention"])
                      if chi_convention_present_s else "kuroki")

    logger.info("Loading FLEX chi_c from: {}".format(chi_c_path))
    data_c = np.load(chi_c_path)
    if not allow_ir:
        _reject_ir_native(data_c, chi_c_path, _STATIC_IR_HINT)
    chi_c_raw = data_c["chiq_c"] if "chiq_c" in data_c else data_c["chiq"]
    # The spin and charge files must share one convention; combining e.g. an MYO
    # chi_s with a Kuroki chi_c would build a meaningless pairing vertex.
    chi_convention_present_c = "chi_convention" in data_c
    chi_convention_c = (str(data_c["chi_convention"])
                        if chi_convention_present_c else "kuroki")
    if chi_convention_c != chi_convention:
        raise ValueError(
            "FLEX chi_s and chi_c have different conventions ('{}' vs '{}'); "
            "they must come from the same run. Check flex_chi_s/flex_chi_c.".format(
                chi_convention, chi_convention_c))

    # S/C vertex-content versioning (#113). The Hund / Exchange / Ising vertex
    # entries were corrected against exact diagonalization; susceptibilities
    # computed with the OLD entries must not be paired with the corrected
    # matrices when the interaction set contains an affected type -- the
    # kernel would silently mix two different interactions. U/V-only inputs
    # are provably unchanged, so legacy files stay usable there.
    if interactions is not None:
        # activation requires actual CONTENT, not key presence: an explicitly
        # configured but empty Hund/Exchange/Ising file contributes nothing
        affected = [t for t in ("Hund", "Exchange", "Ising")
                    if interactions.get(t)]
        # #101: a SECOND independent reason a no-version file is unusable.
        # PR #99 changed the orbital orientation _build_interaction_k stores
        # (which the general/myo S/C matrices consume), and version-2 files
        # are necessarily post-#99 (the version tag came later), so for a
        # myo file the version requirement doubles as the orientation
        # marker. Orientation only matters when an interaction is not
        # invariant under the orbital transpose. PairHop is excluded (its
        # declaration partner is the HERMITIAN entry, a different pairing);
        # PairLift is excluded because its particle-hole S/C contribution
        # is exactly zero on both the producer and consumer sides -- an
        # asymmetric PairLift cannot change either the stored chi or the
        # vertex, so rejecting on it would be a pure false positive.
        orientation = []
        resym = []
        if str(chi_convention).lower() == "myo":
            orientation = [
                t for t in ("CoulombInter", "Hund", "Ising", "Exchange")
                if interactions.get(t)
                and _onsite_transpose_asymmetry(interactions[t]) > 0.0]
            # PairHop is NOT exempt wholesale: its orientation never
            # changed, but the conjugated-mean reading of its Hermitian-
            # partner declarations arrived with the version stamp, so a
            # non-Hermitian-closed declaration has unverifiable
            # symmetrisation semantics in an unversioned file.
            if (interactions.get("PairHop")
                    and _onsite_hermitian_pair_mismatch(
                        interactions["PairHop"]) > 0.0):
                resym = ["PairHop"]
        if affected or orientation or resym:
            versions = {}
            for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
                if "sc_vertex_version" not in data:
                    reasons = []
                    if affected:
                        reasons.append(
                            "the interaction set contains {} whose vertex "
                            "content changed in the #113 corrections".format(
                                ", ".join(affected)))
                    if resym:
                        reasons.append(
                            "PairHop declares on-site entries that are "
                            "not Hermitian partners, whose symmetrised "
                            "reading cannot be verified for a file this "
                            "old: the conjugated-mean reading of PairHop "
                            "declarations arrived with the version stamp "
                            "(#101)")
                    if orientation:
                        reasons.append(
                            "{} declare(s) an asymmetric on-site "
                            "inter-orbital coupling, whose orientation and "
                            "symmetrisation semantics cannot be verified "
                            "for a file this old: the interaction "
                            "orientation changed in #99 and the "
                            "declaration symmetrisation arrived with the "
                            "version-2 stamp, so an unversioned file may "
                            "pair either or both differently from the "
                            "current vertex (#101)".format(
                                ", ".join(orientation)))
                    raise ValueError(
                        "FLEX susceptibility file '{}' predates the current "
                        "vertex conventions (no sc_vertex_version field), "
                        "and {}. Pairing the old susceptibilities with the "
                        "current vertices would silently mix two different "
                        "interactions -- regenerate the susceptibilities "
                        "with the current code.".format(
                            path, "; and ".join(reasons)))
                # strict decode: the tag must be a single integral scalar.
                # A plain int() would silently truncate (2.9 -> 2, accepted)
                # and non-finite values would escape as OverflowError.
                try:
                    arr = np.asarray(data["sc_vertex_version"])
                    if arr.size != 1:
                        raise ValueError("not a scalar")
                    val = complex(arr.reshape(())[()])
                    if val.imag != 0.0 or not (
                            np.isfinite(val.real)
                            and float(val.real).is_integer()):
                        raise ValueError("not an integer")
                    versions[path] = int(val.real)
                except (TypeError, ValueError, OverflowError):
                    raise ValueError(
                        "FLEX susceptibility file '{}' carries a malformed "
                        "sc_vertex_version field ({!r}).".format(
                            path, data["sc_vertex_version"]))
            if len(set(versions.values())) != 1:
                raise ValueError(
                    "FLEX chi_s and chi_c carry different sc_vertex_version "
                    "values ({}); they must come from the same run.".format(
                        versions))
            ver = next(iter(versions.values()))
            if ver != 2:
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = {} "
                    "but this code supports version 2 (the #113 exact-"
                    "diagonalization vertex content); regenerate the "
                    "susceptibilities with the current code.".format(ver))
    # Self-describing index-order marker (issue #78 follow-up). Current files
    # always store [a,c,b,d]; the marker exists so a pre-#78 dev output (the
    # general path stored MYO-transposed arrays under the SAME "myo" tag,
    # indistinguishable by tag alone) is rejected instead of silently building
    # a transposed pairing vertex, and so any future layout change fails fast.
    # Marker-less "kuroki"/untagged files stay readable: the reduced-path
    # layout never changed.
    for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
        if "chi_orbital_layout" in data:
            layout = str(data["chi_orbital_layout"])
            if layout != "acbd":
                raise ValueError(
                    "FLEX chi file '{}' declares chi_orbital_layout='{}' but "
                    "this reader only supports 'acbd'.".format(path, layout))
        elif chi_convention == "myo":
            raise ValueError(
                "FLEX chi file '{}' is tagged chi_convention='myo' but lacks "
                "the chi_orbital_layout marker: it was produced by a pre-fix "
                "development build that stored the arrays orbital-pair-"
                "TRANSPOSED (issue #78) and would yield a wrong pairing "
                "vertex. Re-run FLEX with the current build to regenerate "
                "it.".format(path))
    if not allow_ir:
        return (chi_s_raw, chi_c_raw, chi_convention,
                (_legacy_up_block_only(data_s),
                 _legacy_up_block_only(data_c)))

    native_s, native_c = is_ir_native(data_s), is_ir_native(data_c)
    if native_s != native_c:
        raise ValueError(
            "all dynamic-Eliashberg inputs must share one encoding; re-run "
            "FLEX so chiq_s/chiq_c/green are all densified or all IR-native "
            "(chiq_s: {}, chiq_c: {}).".format(
                "IR-native" if native_s else "densified",
                "IR-native" if native_c else "densified"))
    if native_s and not (chi_convention_present_s and chi_convention_present_c):
        raise ValueError(
            "FLEX chi_s/chi_c are IR-native but missing 'chi_convention'; "
            "the file appears to be an MYO/IR-native FLEX output and its "
            "orbital layout (spin-orbital 'kuroki' vs orbital-pair 'myo') "
            "cannot be safely defaulted -- re-run FLEX with a build that "
            "tags chi_convention, or re-tag the npz explicitly.")
    ir_meta = ({"chis": ir_native_meta(data_s),
                "chic": ir_native_meta(data_c)} if native_s else None)
    return (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
            (_legacy_up_block_only(data_s),
             _legacy_up_block_only(data_c)))


#: Relative size, against the kept spin-up block, below which discarded spin
#: content is attributed to round-off rather than physics.
#:
#: A few hundred ulp of double precision. This is a NARROW MARGIN, not a
#: validated error bound: the error a linear solve can accumulate depends on
#: dimension, conditioning, backend and operation order, and near a
#: susceptibility pole a legitimate paramagnetic producer could exceed it. The
#: claim it rests on is only the measured one -- every producer exercised here
#: (CPU, uniform and IR axes, static and dynamic, plus production multi-orbital
#: output) is bit-exact, so the margin is never used in practice and exists so
#: that a backend whose solve is not bit-symmetric degrades to a warning rather
#: than aborting.
#:
#: It is NOT a claim about how weak a physical field can be. An earlier draft
#: used 1e-8, justified by a measured Zeeman case at ratio ~1.6 -- which sets no
#: lower bound at all, and would have relabelled a genuinely weak field as
#: round-off. If a supported backend is ever found to exceed this margin
#: legitimately, widen it deliberately with that measurement in hand rather than
#: treating the constant as already covering it.
_SPIN_DISCARD_ROUNDOFF_RATIO = 256 * np.finfo(float).eps


def _check_spin_block_discarded(chi_raw, norb, convention, label="chi",
                                legacy_up_block_only=False):
    """Reject, or at minimum report, spin content the embedding would drop.

    The reduced spin-orbital index is spin-block ordered ``s*norb + a``, and the
    embedding below keeps ONLY the spin-up block ``[:norb, :norb]``. Everything
    else -- the down block and both cross-spin blocks -- is dropped, because
    nothing downstream carries spin: the Kuroki S/C matrices and the
    singlet/triplet vertex formulas are norb^2-sized and paramagnetic.

    Dropping it is lossless exactly when it is redundant: the down block equals
    the up block and the cross blocks are zero. That holds bit-for-bit for a
    paramagnetic run -- the inflation writes the same array into both spin
    blocks and never touches the cross blocks
    (``FLEX._inflate_chi0q_and_ham``), and the channel vertices are spin-block
    diagonal, so the RPA solve preserves both properties.

    When it does not hold the discarded data is real and the eigenvalue is not a
    controlled approximation of anything, so this RAISES rather than warning:
    a spin-polarized model is outside what this formulation can express, and
    returning a number for it is worse than refusing.

    The test is on the DISCARDED DATA, not on the run's spin mode. Inferring
    "spin-polarized" from unequal diagonal blocks would be wrong in both
    directions: a spinful/spin-orbit run can have unequal blocks while still
    being time-reversal symmetric, and it can equally have equal diagonal blocks
    with nonzero cross-spin blocks -- discarded just the same. The stored npz
    records no ``spin_mode``, so the data is the only signal available.

    Two severities, because a false positive costs very different amounts:
    anything nonzero is reported, but only content above
    ``_SPIN_DISCARD_ROUNDOFF_RATIO`` -- a few hundred ulp of double precision --
    aborts. Exact bit-equality was confirmed on the uniform and IR axes, static
    and dynamic, and on production multi-orbital output, so the window is there
    purely so a backend whose solve is not bit-symmetric (a GPU batched solve,
    say) degrades to a warning instead of killing a legitimate paramagnetic run.
    It is deliberately anchored to machine precision and not to any assumption
    about how weak a physical field can be.
    """
    if str(convention).lower() != "kuroki":
        # Only the reduced/squashed layout has spin blocks to compare.
        return
    ns = 2
    chi = np.asarray(chi_raw)
    if norb <= 0:
        return
    if chi.ndim != 4 or chi.shape[-1] != chi.shape[-2]:
        # Rank AND squareness, before any dimension-based early return. Checking
        # only the last axis let a rectangular array such as (nfreq, nvol, 2, 8)
        # slip past -- its last dimension is not norb*ns, so the check returned,
        # and _expand_flex_chi then flattened 2*8 = 16, took sqrt, and
        # reinterpreted the rectangle as a 4x4 susceptibility.
        raise ValueError(
            "spin-block check expects a (nfreq, nvol, nd, nd) susceptibility "
            "with a square trailing matrix, got shape {}.".format(chi.shape))
    if chi.shape[-1] != norb * ns:
        # Not the reduced spin-orbital layout; _expand_flex_chi gives a better
        # diagnostic for this (it can recognise the spin-orbital-mode case).
        return
    up = chi[..., :norb, :norb]
    if up.size == 0:
        return
    down = chi[..., norb:, norb:]
    cross_ud = chi[..., :norb, norb:]
    cross_du = chi[..., norb:, :norb]

    # Accumulate per frequency rather than over the whole axis. Keeping the
    # whole-axis POLICY (see the loaders) does not require whole-axis
    # TEMPORARIES: abs(up), up-down and abs(up-down) at full size are each of
    # the order of the stored array, which on a production file would add
    # multi-GB allocations and could fail a valid file on memory alone.
    nfreq = chi.shape[0]

    # Per-frequency RATIOS, not a global difference over a global scale. The
    # ratio decides acceptance, and a single global pair lets a large redundant
    # frequency mask a small completely non-redundant one: with frequency 0
    # redundant at scale 1e14 and frequency 1 holding up=1, down=0, the true
    # mismatch is 100% while the global ratio is 1e-14 -- under the threshold.
    #
    # Non-finite values are rejected outright here rather than folded into the
    # maxima: max(0.0, nan) is 0.0 in Python, so a NaN in a discarded block
    # would leave the running maximum untouched and the array would be accepted
    # as exactly redundant.
    scale = d_down = d_cross = 0.0
    worst_ratio = 0.0
    worst_at = (0, 0.0, 0.0, 0.0)
    any_down = any_cross = False
    for w in range(nfreq):
        u, dn = up[w], down[w]
        cu, cd = cross_ud[w], cross_du[w]
        for name, blk in (("up", u), ("down", dn),
                          ("up-down cross", cu), ("down-up cross", cd)):
            if blk.size and not np.all(np.isfinite(blk)):
                raise ValueError(
                    "The reduced FLEX susceptibility {} has non-finite values "
                    "(NaN or Inf) in its {} block at frequency index {}. NaN "
                    "compares false against everything, so without this check "
                    "it would leave the redundancy maxima untouched and the "
                    "array would be accepted as exactly "
                    "redundant.".format(label, name, w))
        w_scale = float(np.max(np.abs(u)))
        w_down = float(np.max(np.abs(u - dn)))
        w_cross = max(float(np.max(np.abs(cu))), float(np.max(np.abs(cd))))
        scale = max(scale, w_scale)
        d_down = max(d_down, w_down)
        d_cross = max(d_cross, w_cross)
        w_ratio = max(w_down, w_cross) / max(w_scale, np.finfo(float).tiny)
        if w_ratio > worst_ratio:
            worst_ratio = w_ratio
            worst_at = (w, w_scale, w_down, w_cross)
        any_down = any_down or bool(np.any(dn))
        any_cross = any_cross or bool(np.any(cu)) or bool(np.any(cd))

    if d_down == 0.0 and d_cross == 0.0:
        return

    # An all-zero down/cross block is NOT a reliable marker of the legacy
    # single-populated representation: a saturated or projected spin sector, or
    # externally produced spin-resolved data, can look exactly the same. Absent
    # data and fully polarized data are indistinguishable from the values alone,
    # so the file has to SAY which it is. Accept it only when it does.
    if scale > 0.0 and not any_down and not any_cross:
        if legacy_up_block_only:
            # Name the actual source of the authorization: a file tag and a user
            # override carry different weight, and reporting one as the other
            # would misrepresent who vouched for the data.
            via = ("the [eliashberg] accept_up_block_only setting"
                   if legacy_up_block_only == "config_override"
                   else "the file's chi_spin_blocks='up_only' tag")
            logger.warning(
                "The reduced FLEX susceptibility %s has an all-zero down-spin "
                "block and all-zero cross-spin blocks. Per %s this is read as "
                "the layout in which only the block the Eliashberg step "
                "consumes was ever populated -- not as a spin-polarized run -- "
                "and the up-spin block is used.", label, via)
            return
        raise ValueError(
            "The reduced FLEX susceptibility {} has a nonzero up-spin block "
            "with an identically zero down-spin block and zero cross-spin "
            "blocks. That is either a file in which only the consumed block was "
            "ever populated, or a fully spin-polarized/projected susceptibility "
            "-- the values alone cannot tell them apart, and the second is not "
            "something this paramagnetic formulation can use. If you know it is "
            "the former, say so: set [eliashberg] accept_up_block_only = true, "
            "or tag the file itself with chi_spin_blocks='up_only'. Files "
            "written before this check existed carry no tag, so the "
            "configuration flag is the route for them.".format(label))

    # Describe the frequency that decided it: acceptance is a per-frequency
    # ratio, so quoting global maxima could pair a small tail's mismatch with a
    # huge unrelated scale from elsewhere on the axis.
    w_at, scale, d_down, d_cross = worst_at
    parts = []
    if d_down != 0.0:
        parts.append("the down-spin block differs from the up-spin block "
                     "by {:.3e}".format(d_down))
    if d_cross != 0.0:
        parts.append("the cross-spin blocks are nonzero "
                     "(max {:.3e})".format(d_cross))
    detail = " and ".join(parts)
    ratio = worst_ratio

    if ratio <= _SPIN_DISCARD_ROUNDOFF_RATIO:
        logger.warning(
            "The reduced FLEX susceptibility %s is not exactly spin-redundant: "
            "%s, against an up-block scale of %.3e at frequency index %d. "
            "That is %.1e of the kept block there -- within round-off of "
            "double precision, so it is treated as noise and the calculation "
            "continues using the up-spin block. Please report it: every tested "
            "path produces "
            "bit-identical spin blocks.", label, detail, scale, w_at, ratio)
        return

    raise ValueError(
        "The reduced FLEX susceptibility {} is SPIN-RESOLVED: {}, against an "
        "up-block scale of {:.3e} at frequency index {}. The Eliashberg "
        "pairing vertex here is "
        "paramagnetic -- the Kuroki S/C matrices carry no spin index and the "
        "singlet/triplet decomposition assumes spin-rotational symmetry -- so "
        "only the up-spin block could be used and the rest would be discarded. "
        "That is not a controlled approximation of the spin-resolved problem, "
        "so it is refused rather than silently returning a number. If you meant "
        "to run a paramagnetic calculation and a field such as coeff_extern is "
        "the only source of the splitting, drop it and re-run FLEX. If the spin "
        "structure is intrinsic to your model, this solver cannot describe it: "
        "a spin-resolved treatment needs an S_z-resolved vertex and the "
        "transverse susceptibility chi^(+-), which FLEX does not compute at all "
        "(calc_type='ring+ladder' is rejected on every FLEX scheme).".format(
            label, detail, scale, w_at))


def _expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention):
    """Reshape H-wave chi ``(nfreq, nvol, nd, nd)`` to
    ``(nfreq, Nx, Ny, Nz, nd, nd)`` in the ``nd = norb^2`` Eliashberg space,
    resolving the spin-orbital-vs-orbital-pair layout from ``convention``.

    The two FLEX conventions have DIFFERENT physical layouts that shape alone
    cannot always tell apart:

    - ``"myo"`` (general full-vertex FLEX) is already in orbital-pair space
      ``nd_chi = norb^2``; passed through unchanged.
    - ``"kuroki"`` (reduced / squashed FLEX) is in spin-orbital reduced space
      ``nd_chi = norb*ns`` (spin-block ordered ``s*norb + a``), and its matrix
      index is a DENSITY PAIR: ``X[a, b]`` is ``chi_{(a,a),(b,b)}``.  The
      spin-up block ``[:norb, :norb]`` is extracted and embedded at the
      density-pair positions ``[a*norb + a, b*norb + b]`` of the
      ``norb^2 x norb^2`` orbital-pair space, with every off-density component
      left exactly zero (see the body for why any other placement is wrong).

    Whether the spin-orbital block must be extracted is decided by ``nd_chi``:
    ``nd_chi == norb*ns`` means spin-orbital (extract), ``nd_chi == norb^2``
    means orbital-pair (pass through). These two collide ONLY for ``norb == 2``
    (``norb^2 == norb*ns == 4``); there the shape is ambiguous and the layout is
    resolved from ``convention`` (``"kuroki"`` -> spin-orbital extract,
    ``"myo"`` -> orbital-pair). The previous shape-only heuristic silently
    treated a norb=2 kuroki spin-orbital chi as orbital-pair, skipping the
    extraction and building a wrong pairing vertex.

    For ``norb != 2`` the shape is unambiguous, but ``convention`` is still
    REQUIRED to agree with the shape-inferred layout (and to be a recognized
    value): the caller forwards ``convention`` unchanged to
    ``_compute_vertices_flex``, which uses it (not the shape) to pick the
    MYO vs Kuroki S/C matrices, so a shape/tag mismatch would silently build
    the wrong pairing vertex just as in the norb=2 case.

    The mapping is elementwise in frequency, so it may be applied to a single
    static slice or the full axis identically.
    """
    ns = 2
    nd = norb * norb          # orbital-pair dimension
    nd_so = norb * ns         # spin-orbital reduced dimension
    nfreq = chi_raw.shape[0]
    # Derive nd_chi from the trailing AXES, not from sqrt of the flattened
    # element count: the latter accepts any rectangle whose product happens to
    # be a perfect square -- (nfreq, nvol, 2, 8) has 16 trailing elements and
    # was silently reshaped into a 4x4 susceptibility.
    #
    # Two layouts arrive here: the reduced/kuroki one is (nfreq, nvol, nd, nd),
    # the general/myo one is (nfreq, nvol, norb, norb, norb, norb) -- an
    # orbital-pair matrix stored with its pair indices unflattened.
    tail = chi_raw.shape[2:]
    if len(tail) == 2 and tail[0] == tail[1]:
        nd_chi = tail[0]
    elif len(tail) == 4 and len(set(tail)) == 1:
        nd_chi = tail[0] * tail[1]
    else:
        raise ValueError(
            "FLEX chi must be (nfreq, nvol, nd, nd) or "
            "(nfreq, nvol, norb, norb, norb, norb); got shape {}.".format(
                chi_raw.shape))
    chi_full = chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_chi, nd_chi)

    if nd_chi == nd and nd_chi == nd_so:
        # ambiguous (norb == 2): the convention tag is the only disambiguator,
        # and choosing the wrong branch silently corrupts the pairing vertex, so
        # require an explicitly known tag rather than defaulting on mismatch.
        if convention == "kuroki":
            is_spin_orbital = True
        elif convention == "myo":
            is_spin_orbital = False
        else:
            raise ValueError(
                "norb=2 FLEX chi has the shape-ambiguous dimension nd_chi={} "
                "(norb^2 == norb*ns); a known chi_convention ('myo' or "
                "'kuroki') is required to resolve the layout, got '{}'.".format(
                    nd_chi, convention))
    elif nd_chi == nd_so and nd_chi != nd:
        # Unambiguous spin-orbital shape: the convention tag must still agree,
        # otherwise a file shaped spin-orbital but mistagged "myo" (or with an
        # unrecognized tag) would be extracted correctly here yet the wrong
        # tag would later select the MYO S/C matrices downstream, silently
        # building the wrong pairing vertex.
        if convention not in ("myo", "kuroki"):
            raise ValueError(
                "FLEX chi has unrecognized chi_convention='{}' (expected "
                "'myo' or 'kuroki').".format(convention))
        if convention != "kuroki":
            raise ValueError(
                "FLEX chi dimension nd_chi={} (norb={}) is unambiguously "
                "spin-orbital (nd_so=norb*ns={}, nd=norb^2={}) but is tagged "
                "chi_convention='{}'; expected 'kuroki'.".format(
                    nd_chi, norb, nd_so, nd, convention))
        is_spin_orbital = True
    elif nd_chi == nd and nd_chi != nd_so:
        # Unambiguous orbital-pair shape: same agreement check, mirrored.
        if convention not in ("myo", "kuroki"):
            raise ValueError(
                "FLEX chi has unrecognized chi_convention='{}' (expected "
                "'myo' or 'kuroki').".format(convention))
        if convention != "myo":
            raise ValueError(
                "FLEX chi dimension nd_chi={} (norb={}) is unambiguously "
                "orbital-pair (nd=norb^2={}, nd_so=norb*ns={}) but is tagged "
                "chi_convention='{}'; expected 'myo'.".format(
                    nd_chi, norb, nd, nd_so, convention))
        is_spin_orbital = False
    else:
        hint = ""
        if nd_chi == norb and norb % 2 == 0:
            # The spin-orbital case. FLEX writes chi in its reduced spin-orbital
            # space, nd = norb_phys * ns; with geom.dat's norb being the
            # spin-orbital count (= 2 * norb_phys) that lands on nd_chi == norb.
            # Say so -- the bare dimensions look like a malformed file.
            hint = (" nd_chi == norb, which is what a "
                    "spin-orbital FLEX run (mode.enable_spin_orbital = true) "
                    "produces: FLEX writes chi with the PHYSICAL orbital count "
                    "while hwave_sc reads norb from geom.dat, where "
                    "spin-orbital mode stores the spin-orbital count "
                    "(= 2 x physical). The Eliashberg step does not support "
                    "spin-orbital / spin-mixing models -- its pairing vertex is "
                    "paramagnetic -- so this combination is not usable rather "
                    "than merely misconfigured.")
        raise ValueError(
            "FLEX chi dimension nd_chi={} matches neither the orbital-pair "
            "size norb^2={} nor the spin-orbital size norb*ns={} (norb={})."
            "{}".format(nd_chi, nd, nd_so, norb, hint))

    if not is_spin_orbital:
        return chi_full

    # Spin-orbital reduced (spin-block ordered s*norb+a): extract the spin-up
    # block and embed it at the DENSITY-PAIR positions of the orbital-pair
    # space.
    #
    # The reduced/squashed scheme keeps only the density-density diagonal of the
    # susceptibility: the stored X[a, b] IS chi_{(a,a),(b,b)} (its companion
    # interaction reduction in FLEX._inflate_chi0q_and_ham is
    # einsum('ksasatbtb->ksatb', ...), the density-density diagonal of the
    # vertex).  With the orbital-pair flat index (l1,l2) -> l1*norb + l2 used by
    # the S/C builders, the one faithful embedding is therefore
    #
    #     out[(a,a), (b,b)] = X[a, b],   every other component zero.
    #
    # The historical placement out[(l1,l2), (l3,l2)] = X[l1, l3] (a delta_{l2,l4}
    # "spectator" scatter, i.e. kron(X, I_norb)) instead read the density-pair
    # index as an ordinary orbital index.  That dropped the genuine
    # inter-orbital density coupling chi_{(0,0),(1,1)} -- which S @ chi @ S does
    # reference -- and scattered X into pair positions the reduced scheme never
    # computed, so chi0q_mode="flex" on a reduced run disagreed with the
    # equivalent "load" run for identical Sigma=0 physics.
    #
    # Off-density rows/columns stay exactly zero: a reduced run carries no
    # information about pair indices (a,b) with a != b, and fabricating them
    # from the density block is what caused the mismatch.  See
    # _compute_vertices_flex, which warns when the interaction actually needs
    # those missing components.
    chi_orb = chi_full[:, :, :, :, :norb, :norb]
    out = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
    dens = np.arange(norb) * norb + np.arange(norb)   # flat index of (a,a)
    out[..., dens[:, None], dens[None, :]] = chi_orb
    return out


def _load_flex_green(input_dict, norb, Nx, Ny, Nz, allow_ir=False):
    """Load the FLEX dressed Green's function, or ``None`` if absent.

    Returns shape ``(norb, norb, Nx, Ny, Nz, nfreq)`` (full fermionic axis;
    ``nfreq`` = Nmat for uniform files, the node count for IR-native ones).
    With ``allow_ir=True`` returns ``(green, ir_meta)`` instead (``ir_meta``
    is ``None`` for uniform files) -- arity changes only under that flag.
    """
    _, _, green_path = _resolve_flex_paths(input_dict)
    if not os.path.exists(green_path):
        return (None, None) if allow_ir else None
    logger.info("Loading FLEX dressed Green from: {}".format(green_path))
    data_g = np.load(green_path)
    if not allow_ir:
        _reject_ir_native(data_g, green_path, _STATIC_IR_HINT)
    green_raw = data_g["green"]
    # H-wave format: (nblock, nfreq, nvol, norb, norb)
    if green_raw.ndim != 5:
        raise ValueError(
            "green file '{}': expected a rank-5 array (nblock, nfreq, "
            "nvol, norb, norb), got shape {}.".format(
                green_path, green_raw.shape))
    nblock, nmat_g, nvol, norb1, norb2 = green_raw.shape
    # Validate every dimension BEFORE any use: an element-count-compatible
    # but misshapen file would otherwise be silently reinterpreted by the
    # reshape below (measured: (1, 3, 1, 4, 4) on a 4-site 2-orbital
    # system was accepted and scrambled), an empty frequency axis would
    # contract to an all-zero kernel and return a bogus eigenvalue, and a
    # zero-block file died on an incidental IndexError.
    expected_nvol = Nx * Ny * Nz
    if nblock < 1 or nmat_g < 1:
        raise ValueError(
            "green file '{}': needs at least one spin block and one "
            "frequency, got shape {}.".format(green_path, green_raw.shape))
    if nvol != expected_nvol or norb1 != norb or norb2 != norb:
        raise ValueError(
            "green file '{}': shape {} does not match this run "
            "(nvol = Nx*Ny*Nz = {}, norb = {}).".format(
                green_path, green_raw.shape, expected_nvol, norb))
    # A spin-diag FLEX run writes TWO blocks, G_up and G_down. Taking block 0
    # would discard the down-spin propagator exactly as the reduced chi loader
    # used to discard the down-spin susceptibility -- and the pair bubble built
    # from it feeds the Eliashberg kernel directly. The susceptibility guard
    # normally rejects such a run first, but relying on that is fragile: it is
    # a different file, and a run whose chi happens to be redundant while its
    # Green functions are not would slip straight through. (The dynamic loader
    # even reads the Green function BEFORE the susceptibility check.) So the
    # same policy applies here: discarding the other blocks is lossless
    # exactly when they are redundant; when they are not, the pair bubble
    # cannot represent the stored physics and the loader RAISES unless the
    # user takes responsibility via [eliashberg] accept_up_block_only.
    # Non-finite content is rejected up front (below); the block
    # comparison then runs on finite data with plain component equality.
    # Non-finite content is rejected at the file boundary, as the
    # susceptibility loader already does: a NaN/Inf dressed Green function
    # reaches the pair bubble G2 directly, and the iteration solver can
    # turn the resulting non-finite norm into a SAVED zero eigenvalue with
    # converged=False (round-7 review). Per-frequency scan, bounded
    # temporaries.
    for b in range(nblock):
        for f in range(nmat_g):
            if not np.all(np.isfinite(green_raw[b, f])):
                raise ValueError(
                    "green file '{}' contains non-finite values (block "
                    "{}, frequency index {}); a NaN/Inf dressed Green "
                    "function cannot feed the pair bubble. Re-run the "
                    "producing FLEX calculation.".format(green_path, b, f))

    if nblock > 1:
        def _differs(b):
            # Per-frequency, per-component comparison (finiteness is
            # already established above, so plain equality suffices);
            # frequency-sliced scanning bounds the comparison temporaries
            # to one (nvol, norb, norb) slice. .real/.imag are views.
            for f in range(green_raw.shape[1]):
                gb = green_raw[b, f]
                g0 = green_raw[0, f]
                if not (np.array_equal(gb.real, g0.real)
                        and np.array_equal(gb.imag, g0.imag)):
                    return True
            return False

        for b in range(1, nblock):
            if _differs(b):
                if _accept_up_block_only(input_dict):
                    logger.warning(
                        "green file '%s' carries %d spin blocks with "
                        "DIFFERING content; [eliashberg] "
                        "accept_up_block_only is set, so the up block "
                        "alone feeds the pair bubble and the discarded "
                        "spin content is the user's responsibility.",
                        green_path, nblock)
                    break
                raise ValueError(
                    "green file '{}' carries {} spin blocks whose contents "
                    "differ (block {} != block 0). The Eliashberg pair "
                    "bubble is built from the up block only, so a "
                    "spin-resolved dressed Green function cannot be "
                    "represented -- the discarded block is real physics, "
                    "not redundancy. Re-run FLEX paramagnetically, or set "
                    "[eliashberg] accept_up_block_only = true to take "
                    "responsibility for the approximation.".format(
                        green_path, nblock, b))
    # Grid/temperature provenance (round-3 review). The consumer rebuilds
    # the Matsubara frequencies from ITS OWN beta and assumes the file holds
    # the full centered fermionic grid, so a Green function produced at a
    # different temperature -- or a densified/restricted frequency axis --
    # would be consumed silently with a wrong grid, and the tail correction
    # would then apply a wrong analytic complement. Uniform files written
    # since this change carry beta; older files get a warning, not an error.
    if not is_ir_native(data_g):
        if "freq_index" in data_g:
            fi = np.asarray(data_g["freq_index"]).ravel()
            if fi.size != nmat_g or not np.array_equal(
                    fi, np.arange(nmat_g)):
                raise ValueError(
                    "green file '{}': freq_index does not describe the "
                    "full centered frequency grid of the stored axis "
                    "({} entries for {} frequencies); the pair bubble "
                    "needs the full fermionic grid. Re-run the producing "
                    "FLEX calculation with the full grid.".format(
                        green_path, fi.size, nmat_g))
        run_T = input_dict.get("mode", {}).get("param", {}).get("T")
        if "beta" in data_g:
            beta_arr = np.asarray(data_g["beta"]).ravel()
            # One finite positive scalar, checked BEFORE any comparison:
            # NaN/Inf make every tolerance comparison False (silently
            # "matching"), an empty array would die on an incidental
            # IndexError, and a vector would silently use its first element.
            if (beta_arr.size != 1 or not np.isfinite(beta_arr[0])
                    or not beta_arr[0] > 0):
                raise ValueError(
                    "green file '{}': beta metadata must be a single "
                    "finite positive scalar, got {!r}. Regenerate the "
                    "file.".format(green_path,
                                   np.asarray(data_g["beta"])))
            file_beta = float(beta_arr[0])
            # Symmetric relative tolerance, no absolute floor: an absolute
            # term would wave through large relative mismatches at small
            # beta (high temperature), where the grid differs the most.
            if run_T is not None and not np.isclose(
                    file_beta, 1.0 / float(run_T), rtol=1.0e-8, atol=0.0):
                raise ValueError(
                    "green file '{}' was produced at beta = {} but this "
                    "run uses beta = {}; the Matsubara grid (and the tail "
                    "correction built on it) would be inconsistent. Use "
                    "the producing run's temperature or regenerate the "
                    "file.".format(green_path, file_beta,
                                   1.0 / float(run_T)))
        elif run_T is not None:
            logger.warning(
                "green file '%s' carries no beta metadata (written before "
                "this field existed); the run's beta = %g is assumed to "
                "match the producing FLEX run. Verify the temperatures "
                "agree -- a mismatch silently corrupts the pair bubble.",
                green_path, 1.0 / float(run_T))
    # Convert to sc.py format: (norb, norb, Nx, Ny, Nz, nfreq)
    green = green_raw[0].reshape(
        nmat_g, Nx, Ny, Nz, norb, norb
    ).transpose(4, 5, 1, 2, 3, 0).copy()
    if not allow_ir:
        return green
    meta = ir_native_meta(data_g) if is_ir_native(data_g) else None
    return green, meta


def _load_flex_susceptibilities(input_dict, norb, Nx, Ny, Nz,
                                interactions=None):
    """Load FLEX-computed susceptibilities at the static (zero bosonic
    frequency) limit from NPZ files.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.

    Returns
    -------
    chis : ndarray
        Spin susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd).
    chic : ndarray
        Charge susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd).
    green_dressed : ndarray or None
        Dressed Green's function if available, shape (norb, norb, Nx, Ny, Nz, nmat).
    chi_convention : str
        Orbital convention of chis/chic: "myo" (general full-vertex FLEX) or
        "kuroki" (reduced FLEX / legacy files). Pass this to
        _compute_vertices_flex so the matching S/C matrices are used.
    """
    # Read the raw chi (H-wave layout, frequency = axis 0) WITHOUT expanding,
    # so the static slice is taken before the spin-orbital expansion -- the full
    # loader would otherwise allocate the whole Nmat-long expanded array only to
    # keep one frequency (a memory regression proportional to Nmat).
    (chi_s_raw, chi_c_raw, chi_convention,
     legacy_tags) = _read_flex_chi_raw(
        input_dict, interactions=interactions)

    # The zero bosonic frequency is located via the freq_index/nmat metadata
    # (RPA chiq files can carry a restricted matsubara_frequency axis whose
    # center is NOT the static limit); FLEX files always hold the full grid.
    # Re-derive the same paths the raw reader used so the metadata read here
    # is guaranteed to describe the same files.
    chi_s_path, chi_c_path, _ = _resolve_flex_paths(input_dict)
    config_nmat = input_dict.get("mode", {}).get("param", {}).get(
        "Nmat", _DEFAULT_NMAT)

    def _static_center(path, nfreq):
        data = np.load(path)
        fi, fn = _read_freq_meta(data)
        pos = _static_freq_position(fi, nfreq, config_nmat, path,
                                    file_nmat=fn)
        # no usable metadata: center of the actual data axis (axis 0 in the
        # H-wave chiq layout)
        return nfreq // 2 if pos is None else pos

    center_s = _static_center(chi_s_path, chi_s_raw.shape[0])
    center_c = _static_center(chi_c_path, chi_c_raw.shape[0])

    # Slice the static frequency FIRST, then expand only that single slice.
    # Check the WHOLE stored axis, not just the slice consumed here.
    #
    # Both choices have a failure mode and neither case is common: validating
    # only the static slice lets through a run that is redundant at omega=0 and
    # polarized elsewhere, which then returns a paramagnetic eigenvalue in
    # silence; validating everything rejects a file whose unused frequencies are
    # malformed. The first failure is silent and the second is loud and
    # diagnosable, so the axis is checked in full. Paramagnetism is a property
    # of the producing run rather than of one frequency, which is the same
    # reason.
    legacy_s, legacy_c = legacy_tags
    if _accept_up_block_only(input_dict):
        legacy_s = legacy_c = "config_override"
    _check_spin_block_discarded(chi_s_raw, norb, chi_convention, "chi_s",
                                legacy_s)
    _check_spin_block_discarded(chi_c_raw, norb, chi_convention, "chi_c",
                                legacy_c)
    chis = _expand_flex_chi(chi_s_raw[center_s:center_s + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]
    chic = _expand_flex_chi(chi_c_raw[center_c:center_c + 1],
                            norb, Nx, Ny, Nz, chi_convention)[0]

    green_dressed = _load_flex_green(input_dict, norb, Nx, Ny, Nz)
    return chis, chic, green_dressed, chi_convention


# ---------------------------------------------------------------------------
# G2 and Eliashberg kernel
# ---------------------------------------------------------------------------

def _coerce_g2_tail(value):
    """Strictly parse the [eliashberg] g2_tail switch.

    backend.as_bool reads any unrecognized string ("ture", "garbage") as
    False and any nonzero integer as True -- for a switch that changes the
    physics of every reported eigenvalue, a spelling error must fail loudly
    instead of silently flipping the result (round-3 review). Accepted:
    real booleans, integers 0/1, and the strings true/false/yes/no/on/off/0/1
    (case- and whitespace-insensitive).
    """
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, np.integer)) and value in (0, 1):
        return bool(value)
    if isinstance(value, str):
        v = value.strip().lower()
        if v in ("true", "yes", "on", "1"):
            return True
        if v in ("false", "no", "off", "0"):
            return False
    raise ValueError(
        "[eliashberg] g2_tail = {!r} is not a recognized boolean; use "
        "true/false (or yes/no, on/off, 0/1).".format(value))


# Threshold for the asymptotic-regime diagnostic below: at the window edge
# |i wn G(k, i wn) - I| measures how far the Green function still is from
# its leading 1/(i wn) tail. 0.5 flags the regime where the neglected
# scales are at least comparable to the largest retained frequency.
_G2_TAIL_EDGE_DEV_THRESHOLD = 0.5


def _warn_if_g2_tail_outside_asymptotic_regime(green_kw, beta):
    """Warn when the tail correction's asymptotic premise looks violated.

    The correction subtracts the model I/(i wn) x -I/(i wn) inside the
    window and adds its exact full sum, which IMPROVES the result only when
    the largest retained frequency exceeds the relevant energy/self-energy
    scales -- otherwise the model does not describe the summand anywhere in
    the window and the added beta/4-complement can overshoot (measured:
    H = diag(1, -1), beta = 10, Nmat = 2 makes the corrected cross-orbital
    channel WORSE than the bare sum while staying PSD, so the positivity
    diagnostic alone cannot catch it). The same check flags a loaded Green
    function whose leading coefficient is not the identity (e.g. a scaled
    or truncation-damaged file): for those the fixed unit-coefficient
    correction is wrong by construction.

    Returns the measured edge deviation max_k |i wn G(k, i wn) - I| over
    both window endpoints.
    """
    norb = green_kw.shape[0]
    nmat = green_kw.shape[-1]
    if nmat < 1:
        # Nothing to measure on an empty axis; return neutrally so the
        # public ordering (this diagnostic, then _calc_g2) surfaces the
        # actionable even-grid ValueError instead of an IndexError here.
        return 0.0
    eye = np.eye(norb).reshape(norb, norb, 1, 1, 1)
    dev = 0.0
    for idx in (0, nmat - 1):
        wn = (2.0 * idx + 1.0 - nmat) * np.pi / beta
        dev = max(dev, float(np.abs(1j * wn * green_kw[..., idx]
                                    - eye).max()))
    if dev > _G2_TAIL_EDGE_DEV_THRESHOLD:
        logger.warning(
            "g2_tail: at the largest retained Matsubara frequency the Green "
            "function deviates from its asymptotic tail I/(i wn) by %.2f "
            "(threshold %.2f). The tail correction is asymptotic and can "
            "OVERSHOOT in this regime; the result may be less accurate than "
            "the bare sum, and passing the positivity check does not certify "
            "accuracy. Increase Nmat until the window edge exceeds the "
            "relevant energy scales, or check the loaded Green function's "
            "normalization/provenance.",
            dev, _G2_TAIL_EDGE_DEV_THRESHOLD)
    return dev


def _calc_g2(green_kw, beta, tail=True):
    """Calculate G2 = T * sum_n G(k, wn) G(-k+q, -wn).

    The temperature factor T = 1/beta is included so that the
    Eliashberg kernel correctly computes:
        lambda * sigma(k) = -(T/N_L) sum_{k',n'} P(k-k') G(k') G(-k') sigma(k')

    With ``tail=True`` (the default) the truncated Matsubara sum gets the
    analytic high-frequency tail correction (issue #86). The summand's exact
    leading tail is delta_ij delta_lm / wn^2 -- the 1/(i wn) coefficient of
    G is the identity by completeness of the eigenbasis, and a self-energy
    only enters at the next order, so this holds for dressed (FLEX) Green
    functions too, with no free coefficient (unlike chi0q's user-tunable
    ``coeff_tail``). Subtracting that model inside the window and adding its
    exact full sum (1/beta) sum_{n in Z} 1/wn^2 = beta/4 amounts to adding

        c = beta/4 - (1/beta) sum_{n in window} 1/wn^2   (> 0)

    to G2[i,i,l,l], i.e. c times the identity on the (il) gap space. This
    reduces the truncation error from O(1/Nmat) to the next order and in
    practice restores the positive semi-definiteness of G2 that makes the
    static Eliashberg kernel similar to a Hermitian matrix (real spectrum);
    the higher-order remainder can still leave a small exact eigenvalue
    negative, which _warn_if_g2_indefinite reports. The bare truncated sum
    can be slightly indefinite, which then injects spurious imaginary parts
    into the reported eigenvalues at small Nmat.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    beta : float
        Inverse temperature.

    Returns
    -------
    G2 : ndarray
        Shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    norb = green_kw.shape[0]
    Nx, Ny, Nz, nmat = green_kw.shape[2], green_kw.shape[3], green_kw.shape[4], green_kw.shape[5]
    nvol = Nx * Ny * Nz

    # G(-k, -wn): centered-Matsubara flip on the frequency axis, then the
    # shared FFT-grid map k -> -k on the spatial axes.
    green_kw_inv = reverse_fft_axes(green_kw[..., ::-1], (2, 3, 4))
    # einsum("ijpqsk, lmpqsk -> ijlmpqs") sums over k (nmat dimension)
    # Reshape to use tensordot for BLAS: contract last axis (nmat) after
    # merging spatial dims
    # A: (norb, norb, nvol, nmat) -> (norb^2, nvol*nmat) -- but need per-site sum
    # Better: reshape to (norb*norb, nvol, nmat) then per-site outer product
    A = green_kw.reshape(norb * norb, nvol, nmat)       # (ij, site, k)
    B = green_kw_inv.reshape(norb * norb, nvol, nmat)   # (lm, site, k)
    # G2[ij, lm, site] = sum_k A[ij, site, k] * B[lm, site, k]
    # Per-site batched GEMM over the Matsubara axis n. numpy.einsum does NOT
    # lower 'isn,jsn->ijs' to BLAS GEMM; move site to the batch axis and matmul:
    #   As[s] @ Bs[s].T -> [i, j] per site; moveaxis(...,0,2) -> (i, j, s).
    As = np.moveaxis(A, 1, 0)             # (nvol, norb^2, nmat)
    Bs = np.moveaxis(B, 1, 0)             # (nvol, norb^2, nmat)
    G2 = np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)  # (norb^2, norb^2, nvol)
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz)
    G2 = G2 / beta
    if tail:
        # The grid below must match _calc_green's centered Matsubara grid
        # wn = (2n + 1 - nmat) * pi / beta, n = 0..nmat-1. That grid is only
        # fermionic for EVEN nmat: an odd nmat puts wn = 0 in the window
        # (divide by zero below, and a bosonic frequency in a fermionic sum),
        # and nmat <= 0 would return the bare analytic shift beta/4 with no
        # Green-function samples at all. Loaded FLEX Green functions reach
        # this point without any earlier grid validation, so guard here.
        if nmat <= 0 or nmat % 2 != 0:
            raise ValueError(
                "the Matsubara tail correction (g2_tail) requires an even, "
                "positive number of frequencies on the centered fermionic "
                "grid; the Green function has nmat = {}. Fix the frequency "
                "axis of the input Green function (or Nmat), or set "
                "[eliashberg] g2_tail = false.".format(nmat))
        wn = (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi / beta
        coeff = beta / 4.0 - np.sum(1.0 / wn**2) / beta
        di = np.arange(norb)
        # G2[i, i, l, l] += coeff for every (i, l), all k
        G2[di[:, None], di[:, None], di[None, :], di[None, :]] += coeff
    return G2


# _warn_if_g2_indefinite skip thresholds (module-level so tests can patch
# them): eigvalsh work scales as (norb^2)^3 per k-point, i.e. norb^6 * nvol
# flops up to a constant. The peak memory is ~three full complex copies of
# the (nvol, norb^2, norb^2) view (the transpose/reshape materializes one,
# the residual |M - M^dag| evaluation another, the Hermitized matrix a
# third; tracemalloc-measured peak ~1.5x the two-copy figure), each
# 16 * norb^4 * nvol bytes.
_G2_CHECK_MAX_WORK = 20_000_000_000
_G2_CHECK_MAX_BYTES = 2_000_000_000


def _warn_if_g2_indefinite(G2, norb, tail_enabled):
    """Warn when the pair bubble G2 has a significantly negative eigenvalue.

    G2 viewed as a matrix on the (i,l) gap space, M[(i,l),(j,m)] = G2[i,j,l,m],
    is Hermitian positive semi-definite in the exact Matsubara sum; that is
    what makes the static kernel K = -Gamma W similar to a Hermitian matrix
    and its spectrum real. Truncation breaks positivity (issue #86), and
    users then see complex eigenvalues that look like a broken symmetry.
    Point them at Nmat / g2_tail instead of leaving that misread silent.

    A significantly non-Hermitian G2 gets its own warning first (a loaded
    Green function can be malformed, non-causal, or on the wrong grid;
    silently Hermitizing would hide exactly that defect), and the PSD test
    then runs on the Hermitian part. The check diagonalizes an
    (norb^2 x norb^2) block per k-point, so it is skipped (with a log line)
    when the eigensolver work estimate norb^6 * nvol or the temporaries'
    size would rival the solve itself.
    """
    nvol = int(np.prod(G2.shape[4:]))
    work = norb**6 * nvol
    nbytes = 3 * 16 * norb**4 * nvol
    if work > _G2_CHECK_MAX_WORK or nbytes > _G2_CHECK_MAX_BYTES:
        logger.info("G2 positivity check skipped (work estimate norb^6 * "
                    "nvol = %d, temporaries = %d bytes).", work, nbytes)
        return None
    M = G2.reshape(norb, norb, norb, norb, nvol)
    M = M.transpose(4, 0, 2, 1, 3).reshape(nvol, norb * norb, norb * norb)
    scale = float(np.abs(M).max()) if M.size else 0.0
    herm_residual = (float(np.abs(M - np.conj(M.transpose(0, 2, 1))).max())
                     if M.size else 0.0)
    if herm_residual > 1.0e-8 * max(scale, 1.0e-300):
        logger.warning(
            "G2 (pair bubble) is significantly non-Hermitian on the gap "
            "space: max |M - M^dag| = %.3e against max |M| = %.3e. The "
            "exact pair bubble is Hermitian, so check the input Green "
            "function's provenance (grid, temperature, causality); the "
            "positivity diagnostic below uses only the Hermitian part.",
            herm_residual, scale)
    M = 0.5 * (M + np.conj(M.transpose(0, 2, 1)))
    ev = np.linalg.eigvalsh(M)
    min_ev, max_ev = ev.min(), ev.max()
    if min_ev < -1.0e-6 * max(max_ev, 1.0e-300):
        hint = ("increase Nmat" if tail_enabled
                else "increase Nmat or re-enable the tail correction "
                     "([eliashberg] g2_tail = true)")
        logger.warning(
            "G2 (pair bubble) is not positive semi-definite: min eigenvalue "
            "= %.3e (max = %.3e). The most likely cause is Matsubara "
            "truncation error rather than a broken symmetry, and complex "
            "Eliashberg eigenvalues reported below are then spurious at "
            "this level -- %s. If the Green function was loaded from file, "
            "also check that its grid and temperature match this run.",
            min_ev, max_ev, hint)
    return min_ev


def _eliashberg_kernel_fft(V_q, G2, sigma_old, norb):
    """Apply one iteration of the Eliashberg kernel using FFT convolution.

    Supports both simple (2-index) and general (4-index) pairing vertices.

    For the simple case (V_q shape: norb, norb, Nx, Ny, Nz):
        sigma_{il}(k) = -sum_q V_{ij}(q) * G2_{ijlm}(q) * sigma_{jm}(k)

    For the general case (V_q shape: norb, norb, norb, norb, Nx, Ny, Nz):
        sigma_{l1l4}(k) = -sum_q sum_{l2l3} V_{l1l2,l3l4}(q) * F_{l2l3}(q)
        where F_{l2l3}(q) = sum_{l5l6} G_{l2l5}(k-q) sigma_{l5l6}(k-q) G_{l3l6}(q-k)

    Parameters
    ----------
    V_q : ndarray
        Pairing vertex. Either shape (norb, norb, Nx, Ny, Nz) for simple mode,
        or shape (norb, norb, norb, norb, Nx, Ny, Nz) for general mode.
    G2 : ndarray
        Two-particle Green's function, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    sigma_old : ndarray
        Previous gap function, shape (norb, norb, Nx, Ny, Nz).
    norb : int
        Number of orbitals.

    Returns
    -------
    sigma_new : ndarray
        Updated gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = sigma_old.shape[-3], sigma_old.shape[-2], sigma_old.shape[-1]
    nvol = Nx * Ny * Nz

    # G2Sigma contraction: "ijlmpqs, jmpqs -> ilpqs"
    # = matrix-vector product treating (j,m) as contracted index
    # G2: (norb,norb,norb,norb,Nx,Ny,Nz) -> (norb, norb*norb, norb, nvol)
    # sigma: (norb,norb,Nx,Ny,Nz) -> (norb*norb, nvol)
    sigma_flat = sigma_old.reshape(norb * norb, nvol)  # (jm, site)
    G2_flat = G2.reshape(norb, norb * norb, norb, nvol)  # (i, jm, l, site) -- wrong
    # Need (i, l, jm, site) for matmul on jm
    # Actually einsum is: G2[i,j,l,m,s] * sigma[j,m,s] -> result[i,l,s]
    # Reshape G2 to (norb, norb, norb*norb, nvol): G2[i, l, jm, s]
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_r = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)
    # G2_r[i, l, jm, s], sigma_flat[jm, s]
    # result[i,l,s] = sum_jm G2_r[i,l,jm,s] * sigma_flat[jm,s]
    G2Sigma = np.einsum('iljs,js->ils', G2_r, sigma_flat).reshape(norb, norb, Nx, Ny, Nz)

    if V_q.ndim == 5:
        # Simple mode: V_q is (norb, norb, Nx, Ny, Nz)
        P_r = ifftn(V_q, axes=(-3, -2, -1))
        G2Sigma_r = ifftn(G2Sigma, axes=(-3, -2, -1))

        Sigma_r = P_r * G2Sigma_r
        sigma_new = fftn(Sigma_r, axes=(-3, -2, -1))
        return -sigma_new

    else:
        # General mode: V_q is (norb, norb, norb, norb, Nx, Ny, Nz)
        F_q = G2Sigma  # (norb, norb, Nx, Ny, Nz) = (l2, l3, ...)

        V_r = ifftn(V_q, axes=(-3, -2, -1))
        F_r = ifftn(F_q, axes=(-3, -2, -1))

        # sigma[i,l,s] = sum_{j,k} V_r[i,j,k,l,s] * F_r[j,k,s]
        # Reshape for matmul: V_r -> (norb, norb^2, norb, nvol), F_r -> (norb^2, nvol)
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        F_r_flat = F_r.reshape(norb * norb, nvol)
        sigma_r = np.einsum('ijls,js->ils', V_r_flat, F_r_flat).reshape(
            norb, norb, Nx, Ny, Nz)

        sigma_new = fftn(sigma_r, axes=(-3, -2, -1))
        return -sigma_new


# ---------------------------------------------------------------------------
# Gap function initialization
# ---------------------------------------------------------------------------

def _initialize_gap(mode, norb, kx_array, ky_array, kz_array):
    """Initialize gap function sigma(k) with specified symmetry.

    Supports 2D and 3D lattice symmetries. The form factors are defined
    in terms of kx, ky, kz and work for any dimension (Nz=1 for 2D).

    Parameters
    ----------
    mode : str
        Initialization mode specifying the gap symmetry.
        Abbreviations cx=cos(kx), sx=sin(kx), etc.

        **Isotropic / s-wave:**
        - "cos"       : cos(kx+ky+kz) (original reference code)
        - "s"         : 1 (constant, k-independent)
        - "s_ext"     : cx*cy + cy*cz + cz*cx (extended s-wave, 3D)
        - "s_ext_2d"  : cx*cy (extended s-wave, 2D)

        **d-wave:**
        - "d_x2y2"    : cx - cy
        - "d_y2z2"    : cy - cz (in-plane d_{y^2-z^2}, even parity / singlet;
                        opposite-sign anti-nodes at (pi,0) and (0,pi))
        - "d_xy"      : sx*sy
        - "d_xz"      : sx*sz
        - "d_yz"      : sy*sz
        - "d_z2"      : 2*cz - cx - cy (d_{3z^2-r^2})

        **p-wave (odd parity):**
        - "p_x"       : sx
        - "p_y"       : sy
        - "p_z"       : sz

        **Other:**
        - "random"    : random (all symmetries mixed)
    norb : int
        Number of orbitals.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.

    Returns
    -------
    sigma : ndarray
        Normalized initial gap function, shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )
    I = np.identity(norb)
    cx, cy, cz = np.cos(kx_mesh), np.cos(ky_mesh), np.cos(kz_mesh)
    sx, sy, sz = np.sin(kx_mesh), np.sin(ky_mesh), np.sin(kz_mesh)

    form_factors = {
        # s-wave
        "cos":      lambda: np.cos(kx_mesh + ky_mesh + kz_mesh),
        "s":        lambda: np.ones((Nx, Ny, Nz)),
        "s_ext":    lambda: cx * cy + cy * cz + cz * cx,
        "s_ext_2d": lambda: cx * cy,
        # d-wave
        "d_x2y2":   lambda: cx - cy,
        "d_y2z2":   lambda: cy - cz,
        "d_xy":     lambda: sx * sy,
        "d_xz":     lambda: sx * sz,
        "d_yz":     lambda: sy * sz,
        "d_z2":     lambda: 2.0 * cz - cx - cy,
        # p-wave
        "p_x":      lambda: sx,
        "p_y":      lambda: sy,
        "p_z":      lambda: sz,
    }

    if mode in form_factors:
        f_k = form_factors[mode]()
        sigma = I[:, :, None, None, None] * f_k[None, None, :, :, :]
    elif mode == "random":
        rng = np.random.default_rng(12345)
        sigma = rng.standard_normal((norb, norb, Nx, Ny, Nz))
    else:
        raise ValueError(
            "Unknown init_gap mode: '{}'. Available: {}".format(
                mode, list(form_factors.keys()) + ["random"]))

    norm = np.linalg.norm(sigma)
    if norm > 0:
        sigma /= norm
    elif mode != "random":
        logger.warning(
            "init_gap='%s' evaluates to zero on the current k-grid: a form "
            "factor built from sin(k_a) vanishes when the a-axis is squashed "
            "(e.g. p_x/d_xy at Nx=1, or p_z/d_xz at Nz=1). The seed is all-"
            "zero -- choose a symmetry that is non-vanishing on this grid "
            "(e.g. p_y/p_z, d_yz or d_y2z2 for an in-plane [1, Ny, Nz] cell).",
            mode)
    return sigma


def _resolve_init_gap(init_gap, pairing_type):
    """Resolve the initial-gap symmetry, defaulting to the channel's parity.

    The Eliashberg kernel preserves parity (it commutes with k -> -k), so an
    even seed cannot reach the odd triplet solution and vice versa. When the
    user does not specify ``init_gap``, default to an even s-wave ('cos') seed
    for singlet and an odd p-wave ('p_x') seed for triplet.

    Parameters
    ----------
    init_gap : str or None
        User-specified gap symmetry, or None to use the default.
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    str
        The gap symmetry mode to pass to _initialize_gap.
    """
    if init_gap is not None:
        return init_gap
    return "p_x" if pairing_type == "triplet" else "cos"


# ---------------------------------------------------------------------------
# Solvers
# ---------------------------------------------------------------------------

def _solve_iteration(green_kw, Vs_q, G2, sigma_init, norb,
                     max_iter=1000, alpha=0.5, tol=1.0e-5, pairing_type=None):
    """Solve linearized Eliashberg equation by self-consistent iteration.

    Parameters
    ----------
    green_kw : ndarray
        Green's function, shape (norb, norb, Nx, Ny, Nz, nmat).
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    sigma_init : ndarray
        Initial gap function.
    norb : int
        Number of orbitals.
    max_iter : int
        Maximum iterations.
    alpha : float
        Mixing parameter (0 < alpha < 1).
    tol : float
        Convergence tolerance.
    pairing_type : str, optional
        "singlet" or "triplet". When given, every iterate is projected onto the
        channel's parity sector (even for singlet, odd for triplet). For a
        centrosymmetric system the kernel commutes with parity, so in exact
        arithmetic the iteration would stay in the seed's sector; the projection
        removes the wrong-parity component that numerical noise (FFT/GEMM
        round-off) reintroduces each step and that the power iteration would
        otherwise amplify toward the dominant eigenpair of the *other* channel.
        The projection is only legitimate when the kernel commutes with parity:
        a cross-sector leakage test (several random probes) is run up front and,
        for a non-centrosymmetric kernel (complex/SOC hopping, vertex odd in k,
        chiral gaps), projection is automatically disabled with a warning and the
        un-projected iteration is used. When None, no projection is applied (the
        historical behavior).

    Returns
    -------
    sigma : ndarray
        Converged gap function.
    eigenvalue : float
        Leading eigenvalue (norm of sigma_new before normalization).
    converged : bool
        Whether convergence was achieved.
    n_iter : int
        Number of iterations performed.
    """
    # Precompute the Eliashberg kernel operator once. The vertex IFFT
    # (V_r = ifftn(Vs_q)) and the G2 reshape/transpose are invariant across
    # iterations (only sigma_old changes); _make_kernel_operator hoists them out
    # of the per-matvec closure, saving one full vertex IFFT (and the G2
    # preprocessing) on every one of the up-to-max_iter iterations. The matvec
    # is numerically identical to _eliashberg_kernel_fft(Vs_q, G2, sigma, norb).
    Nx, Ny, Nz = sigma_init.shape[-3], sigma_init.shape[-2], sigma_init.shape[-1]
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    shape = (norb, norb, Nx, Ny, Nz)

    # Parity projection is only legitimate when the kernel commutes with the
    # parity operator P (Delta_{ab}(k) -> Delta_{ba}(-k), i.e.
    # _reverse_k_and_orbital), which holds for centrosymmetric systems. For a
    # non-centrosymmetric kernel (complex/SOC hopping, a vertex odd in k, chiral
    # gaps) [A, P] != 0; projecting each iterate would then discard physical
    # components and steer the power iteration to the wrong eigenpair. We
    # measure this directly as the cross-sector leakage of the kernel: apply A
    # to even- and odd-projected random probes and ask how much of the result
    # lands in the *opposite* parity sector. Zero leakage <=> [A, P] = 0 <=>
    # projection is valid. Several independent probes are used (one probe could
    # accidentally lie near the nullspace of [A, P] and miss a non-commuting
    # kernel). If any probe leaks above threshold, warn and fall back to the
    # historical un-projected iteration.
    do_project = pairing_type is not None
    if do_project and pairing_type not in ("singlet", "triplet"):
        # Validate up front so an invalid channel always raises, even if the
        # centrosymmetry guard below later disables projection (in which case
        # _project_gap_parity -- which also validates -- is never reached).
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    if do_project:
        rng = np.random.default_rng(0)
        leakage = 0.0
        for _ in range(_PARITY_GUARD_PROBES):
            x = (rng.standard_normal(vec_size)
                 + 1j * rng.standard_normal(vec_size)).reshape(shape)
            for sign in (1.0, -1.0):  # even (+) and odd (-) projected probes
                xp = 0.5 * (x + sign * _reverse_k_and_orbital(x))
                axp = A.matvec(xp.ravel()).reshape(shape)
                # component of A xp in the OPPOSITE parity sector
                opp = 0.5 * (axp - sign * _reverse_k_and_orbital(axp))
                denom = np.linalg.norm(axp) + 1.0e-300
                leakage = max(leakage, np.linalg.norm(opp) / denom)
        if leakage > 1.0e-8:
            logger.warning(
                "Eliashberg kernel does not commute with parity (cross-sector "
                "leakage = {:.2e}); the system appears non-centrosymmetric, so "
                "parity projection for the '{}' channel is disabled and the "
                "un-projected iteration is used.".format(leakage, pairing_type))
            do_project = False

    def _project(sigma):
        return (_project_gap_parity(sigma, pairing_type)
                if do_project else sigma)

    sigma_old = _project(sigma_init.copy())
    if do_project and np.linalg.norm(sigma_old) < 1.0e-12 * (
            np.linalg.norm(sigma_init) + 1.0e-300):
        raise ValueError(
            "Initial gap has no component in the '{}' parity sector; "
            "choose an init_gap of the matching parity (even for singlet, "
            "odd for triplet).".format(pairing_type))

    # The power-iteration loop itself (matvec, projection, normalize,
    # convergence check, mixing) is factored into _solve_leading so a future
    # dynamic-frequency kernel can reuse the identical driver. `A` was already
    # built above (for the leakage probe), so `make_operator` just hands back
    # that same operator rather than rebuilding it. The gap-parity projector
    # operates on the (norb, norb, Nx, Ny, Nz) tensor, so it is wrapped as a
    # flat-vector -> flat-vector `project_fn` for the generic driver; when no
    # projection applies (`do_project` False) `project_fn` is None, matching
    # the historical (unprojected) iteration exactly.
    make_operator = lambda: (A, vec_size)
    project_fn = ((lambda flat: _project(flat.reshape(shape)).ravel())
                 if do_project else None)

    eigenvalue, sigma_flat, info = _solve_leading(
        make_operator, vec_size, "iteration",
        max_iter=max_iter, convergence_tol=tol, init_vec=sigma_old.ravel(),
        alpha=alpha, project_fn=project_fn,
    )

    return sigma_flat.reshape(shape), eigenvalue, info["converged"], info["n_iter"]


def _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz):
    """Create LinearOperator for the Eliashberg kernel.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.

    Returns
    -------
    A : LinearOperator
        Linear operator representing the Eliashberg kernel K.
    vec_size : int
        Size of the flattened vector.
    """
    vec_size = norb * norb * Nx * Ny * Nz
    nvol = Nx * Ny * Nz

    # Backend of the resident invariants (numpy on CPU, cupy when calc_eliashberg
    # parked them on the device). Everything the kernel builds/derives must live
    # on this same backend; assert Vs_q and G2 agree.
    xp = backend.array_module_of(Vs_q)
    if backend.array_module_of(G2) is not xp:
        raise ValueError(
            "_make_kernel_operator: Vs_q and G2 must be on the same backend "
            "(both host or both device).")

    # Precompute invariants that don't change per matvec call:
    # 1. V_r = IFFT(Vs_q) -- used every matvec
    V_r = backend.spatial_ifftn(Vs_q, axes=(-3, -2, -1))

    # 2. Precompute G2 reshaped for contraction
    G2_r = G2.reshape(norb, norb, norb, norb, nvol)
    G2_pre = G2_r.transpose(0, 2, 1, 3, 4).reshape(norb, norb, norb * norb, nvol)

    is_simple = (Vs_q.ndim == 5)

    # General (4-index) mode: the two contractions
    #   G2Sigma[i,l,s] = sum_j G2_pre[i,l,j,s] * sigma_flat[j,s]
    #   sigma_r[i,l,s] = sum_j V_r_flat[i,j,l,s] * F_r_flat[j,s]
    # are, per spatial point s, a small dense matvec over j (size norb*norb).
    # numpy.einsum does NOT lower these to batched BLAS GEMM, so we precompute
    # the operator tensors ONCE in GEMM-friendly (nvol, M=il, K=j) layout and
    # apply them with matmul (s, and the column axis z for matmat, as leading
    # batch axes). This is NOT bit-identical to the einsum (GEMM changes the
    # reduction order) but matches it to ~1e-13.
    if not is_simple:
        # G2_pre is (i, l, j, s); we need (s, (i,l), j).
        G2_gemm = G2_pre.transpose(3, 0, 1, 2).reshape(
            nvol, norb * norb, norb * norb)
        # V_r_flat is (i, j, l, s); the GEMM M-axis (i,l) needs i (axis0) and
        # l (axis2) adjacent, so transpose to (s, i, l, j) -> (s, (i,l), j).
        V_r_flat = V_r.reshape(norb, norb * norb, norb, nvol)
        V_r_gemm = V_r_flat.transpose(3, 0, 2, 1).reshape(
            nvol, norb * norb, norb * norb)

    def matvec(v):
        # scipy hands us a HOST numpy vector; move it onto the invariants'
        # backend before any device op (mixed numpy/cupy operands raise on real
        # CuPy). On the numpy backend this is a no-op.
        v = xp.asarray(v)
        sigma = v.reshape(norb, norb, Nx, Ny, Nz)
        sigma_flat = sigma.reshape(norb * norb, nvol)

        if is_simple:
            G2Sigma = xp.einsum('iljs,js->ils', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            Sigma_r = V_r * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(-3, -2, -1))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)[..., np.newaxis]
            G2Sigma = (G2_gemm @ sig)[..., 0]            # (s, (i,l))
            G2Sigma = G2Sigma.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(-3, -2, -1))
            F_r_flat = F_r.reshape(norb * norb, nvol)
            f = xp.moveaxis(F_r_flat, 1, 0)[..., np.newaxis]
            sigma_r = (V_r_gemm @ f)[..., 0]             # (s, (i,l))
            sigma_r = sigma_r.reshape(nvol, norb, norb).transpose(
                1, 2, 0).reshape(norb, norb, Nx, Ny, Nz)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(-3, -2, -1))

        # Keep complex: for complex hopping (e.g. spin-orbit coupling),
        # non-centrosymmetric models, or chiral gaps the kernel is genuinely
        # complex; projecting to real would discard physical components. Return
        # a HOST array so scipy's ARPACK/power/subspace driver stays on numpy.
        return backend.to_host(-sigma_new).ravel()

    def matmat(B):
        # Batched form of matvec: apply the SAME kernel to every column of the
        # (vec_size, k) block B in a single FFT. The batch/column axis is the
        # LAST axis so the spatial FFT axes keep absolute positions (2, 3, 4)
        # and V_r / G2_pre (no column axis) broadcast over it.
        B = xp.asarray(B)
        if B.ndim == 1:
            return matvec(B)
        k = B.shape[1]
        sigma = B.reshape(norb, norb, Nx, Ny, Nz, k)
        sigma_flat = sigma.reshape(norb * norb, nvol, k)

        if is_simple:
            G2Sigma = xp.einsum('iljs,jsz->ilsz', G2_pre, sigma_flat).reshape(
                norb, norb, Nx, Ny, Nz, k)
            G2Sigma_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            Sigma_r = V_r[..., np.newaxis] * G2Sigma_r
            sigma_new = backend.spatial_fftn(Sigma_r, axes=(2, 3, 4))
        else:
            sig = xp.moveaxis(sigma_flat, 1, 0)          # (nvol, norb*norb, k)
            G2Sigma = G2_gemm @ sig                      # (s, (i,l), k)
            G2Sigma = G2Sigma.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)

            F_r = backend.spatial_ifftn(G2Sigma, axes=(2, 3, 4))
            F_r_flat = F_r.reshape(norb * norb, nvol, k)
            f = xp.moveaxis(F_r_flat, 1, 0)              # (nvol, norb*norb, k)
            sigma_r = V_r_gemm @ f                       # (s, (i,l), k)
            sigma_r = sigma_r.reshape(nvol, norb, norb, k).transpose(
                1, 2, 0, 3).reshape(norb, norb, Nx, Ny, Nz, k)
            sigma_new = backend.spatial_fftn(sigma_r, axes=(2, 3, 4))

        return backend.to_host(-sigma_new).reshape(vec_size, k)

    A = LinearOperator((vec_size, vec_size), matvec=matvec, matmat=matmat,
                       dtype=complex)
    return A, vec_size


def _order_by_seed_overlap(vals, vecs, seed_vec):
    """Order eigenpairs by descending overlap ``|<seed, vec>|`` with a seed
    eigenvector (columns already L2-normalized by ARPACK).

    Used for eigenvector continuation: when a converged eigenvector from a
    neighbouring parameter (e.g. the next temperature) is supplied as a seed,
    the physical branch is the eigenpair whose eigenvector maximally overlaps
    it -- NOT the algebraically largest one (which, near an exceptional point
    of the non-Hermitian kernel, can jump to a different branch). Ties and a
    zero seed fall back to real-part ordering.
    """
    s = np.asarray(seed_vec).ravel()
    ns = np.linalg.norm(s)
    if ns == 0:
        return _order_eigenpairs(vals, vecs)
    s = s / ns
    ov = np.abs(vecs.conj().T @ s) / (
        np.linalg.norm(vecs, axis=0) + 1.0e-300)
    # Primary key: descending overlap. Secondary key: descending real part, so
    # an exact overlap TIE deterministically falls back to the physical
    # real-part ordering (as the docstring promises) instead of the arbitrary
    # order np.linalg.eig / ARPACK happened to return. np.lexsort takes the
    # primary key last.
    idx = np.lexsort((-vals.real, -ov))
    return vals[idx], vecs[:, idx]


def _order_eigenpairs(vals, vecs):
    """Order eigenpairs by descending real part (largest first).

    The superconducting eigenvalue is the algebraically largest one
    (lambda -> 1 at Tc); only positive eigenvalues are physically relevant.
    The leading eigenpair must therefore be the largest by real part, NOT the
    largest by magnitude: a large-magnitude *negative* eigenvalue is an
    unphysical repulsive mode and must not be reported first (nor mask a
    smaller positive eigenvalue).

    Note: ARPACK still finds eigenvalues by magnitude (which='LM'), so a small
    positive eigenvalue masked by much larger negative ones may not be among
    the returned set; request more eigenvalues, or use a shift-invert method
    targeting a positive shift, to resolve such cases.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues.
    vecs : ndarray
        Eigenvectors as columns, shape (vec_size, n).

    Returns
    -------
    vals, vecs ordered so that vals[0] has the largest real part.
    """
    idx = np.argsort(-vals.real)
    return vals[idx], vecs[:, idx]


def _reverse_k_and_orbital(gap):
    """Return the gap evaluated at (-k) with orbital indices transposed.

    For the pairing gap Delta_{ab}(k), the singlet/triplet parity condition is
    Delta_{ab}(k) = +/- Delta_{ba}(-k). This helper builds Delta_{ba}(-k):
    it reverses each spatial axis (k -> -k on the FFT grid, index i -> (N-i)%N)
    and swaps the two orbital indices.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).

    Returns
    -------
    ndarray
        Delta_{ba}(-k), same shape.
    """
    rev = reverse_fft_axes(gap, (2, 3, 4))  # k -> -k on the FFT grid
    rev = np.swapaxes(rev, 0, 1)            # orbital transpose a <-> b
    return rev


def _project_gap_parity(gap, pairing_type):
    """Project a gap onto the parity sector of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). With P the parity operator
    (P: Delta_{ab}(k) -> Delta_{ba}(-k), built by _reverse_k_and_orbital), the
    even/odd projectors are (1 +/- P)/2. The Eliashberg kernel commutes with P,
    so these projectors commute with it and the power iteration can be confined
    to the physical sector -- removing wrong-parity components that numerical
    noise would otherwise let the iteration amplify (see _solve_iteration).

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" (even projector) or "triplet" (odd projector).

    Returns
    -------
    ndarray
        The component of ``gap`` in the channel's parity sector, same shape.
    """
    if pairing_type == "singlet":
        sign = 1.0
    elif pairing_type == "triplet":
        sign = -1.0
    else:
        raise ValueError(
            "Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))
    return 0.5 * (gap + sign * _reverse_k_and_orbital(gap))


def _is_gap_parity(gap, pairing_type, tol=0.9):
    """Test whether a gap has the parity of the pairing channel.

    Singlet gaps are even (Delta_{ab}(k) = +Delta_{ba}(-k)); triplet gaps are
    odd (Delta_{ab}(k) = -Delta_{ba}(-k)). Returns True when the projection of
    the gap onto the requested parity sector retains at least ``tol`` of its
    norm.

    Parameters
    ----------
    gap : ndarray
        Gap function, shape (norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".
    tol : float
        Minimum fraction of the norm that must survive the parity projection.

    Returns
    -------
    bool
    """
    proj = _project_gap_parity(gap, pairing_type)
    n = np.linalg.norm(gap)
    if n == 0:
        return False
    return np.linalg.norm(proj) / n >= tol


def _reorder_eigenpairs_by_parity(vals, gaps, pairing_type):
    """Stable-reorder eigenpairs so those matching the channel parity come first.

    The Eliashberg kernel preserves parity, so its eigenvectors split into even
    (singlet) and odd (triplet) sectors. When solving for a given channel, the
    physical solution is the leading eigenpair of the matching parity, which is
    not necessarily the globally leading one. This keeps every eigenpair (so the
    returned count is unchanged) but promotes the requested-parity ones,
    preserving their existing order.

    Parameters
    ----------
    vals : ndarray
        Eigenvalues, shape (n,).
    gaps : ndarray
        Eigenvectors as gaps, shape (n, norb, norb, Nx, Ny, Nz).
    pairing_type : str
        "singlet" or "triplet".

    Returns
    -------
    vals, gaps reordered with matching-parity eigenpairs first.
    """
    match = np.array([_is_gap_parity(g, pairing_type) for g in gaps])
    if not np.any(match):
        logger.warning(
            "No computed eigenpair matches the requested '%s' parity; the "
            "reported leading gap belongs to the other channel. Increase "
            "num_eigenvalues or check the pairing_type.", pairing_type)
    idx = np.concatenate([np.where(match)[0], np.where(~match)[0]])
    return vals[idx], gaps[idx]


def _shift_from_eigenvalues(vals, factor=0.9):
    """Estimate a shift-invert target near the physical (largest-real) eigenvalue.

    The superconducting eigenvalue is the algebraically largest one, so when a
    positive eigenvalue is present among the sampled values, aim the shift at
    the largest real part (NOT the largest magnitude, which could be a large
    negative repulsive mode that masks it). If no positive eigenvalue was
    sampled (no SC instability in range), fall back to the largest-magnitude
    eigenvalue so shift-invert still tracks the dominant mode.

    Parameters
    ----------
    vals : ndarray
        Sampled eigenvalues (e.g. a few largest-magnitude ARPACK values).
    factor : float
        Fraction of the target eigenvalue to use as the shift.

    Returns
    -------
    float
    """
    real = vals.real
    scale = float(np.max(np.abs(vals)))
    # "positive" means significantly positive relative to the spectrum scale,
    # so numerical-noise eigenvalues (~1e-16) do not pull the shift to zero.
    if scale > 0 and np.any(real > 1e-8 * scale):
        return float(np.max(real)) * factor
    # No significant positive eigenvalue: track the dominant magnitude. Note a
    # purely imaginary / all-zero sample yields a 0.0 shift, which sits at a
    # possible zero mode; physical FLEX inputs have a real dominant eigenvalue.
    return float(vals[np.argmax(np.abs(vals))].real) * factor


def _solve_leading(make_operator, vec_size, solver_mode, num_eigenvalues=10,
                   max_iter=1000, convergence_tol=1.0e-5, init_vec=None,
                   sigma_shift=None, alpha=0.5, project_fn=None, seed_vec=None,
                   spectral_shift=None):
    """Shared leading-eigenpair driver behind the static Eliashberg solvers.

    This holds the ARPACK/shift-invert eigen-selection-and-ordering body of
    ``_solve_eigenvalue`` and the power-iteration loop of ``_solve_iteration``,
    generalized to act on any ``make_operator``/``vec_size`` pair (not just the
    static ``_make_kernel_operator`` result) so a future dynamic-frequency
    kernel can reuse the exact same driver.

    Parameters
    ----------
    make_operator : callable
        Zero-argument callable returning ``(A, op_vec_size)``, i.e. the same
        convention as ``_make_kernel_operator`` (a
        ``scipy.sparse.linalg.LinearOperator`` acting on flat vectors of
        length ``vec_size``, plus that size). Called exactly once.
    vec_size : int
        Size of the flattened vector space that ``A`` acts on.
    solver_mode : str
        "arnoldi", "shift-invert-bicgstab", "shift-invert-gmres", or
        "shift-invert-lgmres" select the ARPACK/shift-invert eigenanalysis
        (the former body of ``_solve_eigenvalue``, excluding its "subspace"
        branch which stays in that wrapper). "iteration" selects the power
        loop (the former body of ``_solve_iteration``).
    num_eigenvalues : int
        Number of eigenvalues to request from ARPACK. Only used by the
        eigenvalue-family modes.
    max_iter : int
        Maximum number of power-iteration steps. Only used by "iteration".
    convergence_tol : float
        Power-iteration convergence tolerance on the normalized-iterate
        difference. Only used by "iteration".
    init_vec : ndarray, optional
        Flat initial vector, shape ``(vec_size,)``, for "iteration" mode
        (already projected onto the desired sector by the caller, if
        applicable). Required for "iteration" mode.
    sigma_shift : float, optional
        Shift-invert target for the "shift-invert-*" modes. If None, it is
        estimated from a preliminary Arnoldi pass, exactly as
        ``_solve_eigenvalue`` does.
    alpha : float
        Power-iteration mixing parameter. Only used by "iteration".
    project_fn : callable, optional
        Flat-vector -> flat-vector projector applied to every power-iteration
        iterate (e.g. the gap-parity projector in ``_solve_iteration``). Only
        used by "iteration"; when None, no projection is applied.

    Returns
    -------
    leading_eigenvalue : complex or float
        The dominant eigenvalue (eigenvalue-family: largest real part, per
        ``_order_eigenpairs``; iteration: the converged/last iterate norm).
    leading_eigenvector : ndarray
        Flat, shape ``(vec_size,)``.
    eig_analysis : dict
        Eigenvalue-family modes: ``{"eigenvalues": vals, "eigenvectors": vecs,
        "sigma_shift": sigma_shift}`` with ``vals`` ordered by descending real
        part (``_order_eigenpairs``) and ``vecs`` the matching eigenvectors as
        columns. "iteration" mode: ``{"converged": bool, "n_iter": int}``.
    """
    if solver_mode == "iteration":
        if init_vec is None:
            raise ValueError("init_vec is required for solver_mode='iteration'")
        A, _ = make_operator()
        sigma_old = init_vec
        eigenvalue = 0.0

        for iteration in range(max_iter):
            sigma_new = A.matvec(sigma_old)
            if project_fn is not None:
                sigma_new = project_fn(sigma_new)
            norm = np.linalg.norm(sigma_new)
            eigenvalue = norm

            # A (possibly projected) iterate can collapse to the zero vector
            # when the kernel annihilates the requested sector; normalizing by
            # `norm` would then produce NaN. `np.linalg.norm` squares its
            # inputs, so it returns *exactly* 0.0 once the iterate underflows;
            # guard precisely on that, plus any non-finite norm, and report
            # non-convergence with the last finite iterate instead.
            if norm == 0.0 or not np.isfinite(norm):
                logger.warning(
                    "Eliashberg iterate collapsed to zero norm at iteration "
                    "{}; the kernel annihilates the requested sector, so the "
                    "eigenvalue cannot be normalized. Returning the previous "
                    "iterate as non-converged.".format(iteration))
                return 0.0, sigma_old, {"converged": False,
                                        "n_iter": iteration + 1}

            diff = np.linalg.norm(sigma_new / norm - sigma_old)
            logger.info("Iteration {:4d}: eigenvalue = {:.6f}, diff = {:.6e}".format(
                iteration, norm, diff))

            if diff < convergence_tol:
                logger.info("Converged at iteration {}".format(iteration + 1))
                return eigenvalue, sigma_new / norm, {"converged": True,
                                                       "n_iter": iteration + 1}

            sigma_old = (1.0 - alpha) * sigma_new / norm + alpha * sigma_old

        logger.warning("Failed to converge after {} iterations".format(max_iter))
        return eigenvalue, sigma_old, {"converged": False, "n_iter": max_iter}

    # Eigenvalue family: ARPACK Arnoldi or shift-invert.
    A, _ = make_operator()

    if not (solver_mode == "arnoldi" or solver_mode.startswith("shift-invert")):
        raise ValueError("Unknown eigenvalue method: {}".format(solver_mode))

    # Validate spectral_shift up front (before any branch, incl. the small-dense
    # early return) so invalid values / incompatible modes fail fast everywhere.
    if spectral_shift is not None:
        if isinstance(spectral_shift, str):
            if spectral_shift != "auto":
                raise ValueError(
                    "[eliashberg] spectral_shift string must be \"auto\", got "
                    "{!r}".format(spectral_shift))
        else:
            try:
                sv = float(spectral_shift)
            except (TypeError, ValueError, OverflowError):
                sv = float("nan")
            if not np.isfinite(sv) or sv <= 0.0:
                raise ValueError(
                    "[eliashberg] spectral_shift must be a positive finite "
                    "number or \"auto\", got {!r}".format(spectral_shift))
        if solver_mode != "arnoldi":
            raise ValueError(
                "[eliashberg] spectral_shift is only supported for "
                "eigenvalue_method='arnoldi', not {!r}".format(solver_mode))

    # sigma_shift is the shift-invert target; it has no effect on the plain
    # arnoldi path (which uses spectral_shift instead). Warn rather than fail so
    # existing configs keep working.
    if sigma_shift is not None and solver_mode == "arnoldi":
        logger.warning(
            "[eliashberg] sigma_shift is ignored for eigenvalue_method="
            "'arnoldi' (it targets the shift-invert methods); use "
            "spectral_shift to bias the arnoldi selection.")

    if vec_size < 1:
        raise ValueError("Eliashberg operator has empty vector space")

    if vec_size <= 2:
        # scipy.sparse.linalg.eigs requires k < N - 1 for LinearOperator input,
        # so the smallest valid dynamic grid (e.g. norb=1, Nk=1, Nmat=2) cannot
        # go through ARPACK. Reconstruct the tiny dense operator directly.
        dense = np.empty((vec_size, vec_size), dtype=complex)
        basis = np.eye(vec_size, dtype=complex)
        for j in range(vec_size):
            dense[:, j] = A.matvec(basis[:, j])
        vals, vecs = np.linalg.eig(dense)
        # Match the ARPACK path below: with a seed eigenvector, track the branch
        # that overlaps it (eigenvector continuation); otherwise order by
        # largest real part (the physical SC eigenvalue), not magnitude.
        if seed_vec is not None:
            vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
        else:
            vals, vecs = _order_eigenpairs(vals, vecs)
        n_keep = min(max(1, num_eigenvalues), vec_size)
        vals = vals[:n_keep]
        vecs = vecs[:, :n_keep]
        return vals[0], vecs[:, 0], {"eigenvalues": vals,
                                     "eigenvectors": vecs,
                                     "sigma_shift": sigma_shift}

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    logger.info("Computing {} eigenvalues with method='{}'...".format(
        max_ev, solver_mode))

    if solver_mode == "arnoldi":
        if spectral_shift is not None:
            # Select the LARGEST-REAL eigenvalue (the physical SC eigenvalue,
            # lambda -> 1 at Tc), which plain which='LM' misses when a small
            # positive lambda is masked by larger repulsive (negative)
            # eigenvalues. We ask ARPACK for the largest real part directly
            # (which='LR') -- unlike which='LM', this is the correct criterion
            # even for the non-Hermitian kernel's complex eigenvalues, where
            # |lambda+sigma| would otherwise be dominated by a large-|Im| mode.
            # The real spectral shift A -> A + sigma*I preserves the real-part
            # ordering (Re(lambda+sigma) = Re(lambda) + sigma) but moves the
            # spectrum into the right half-plane, which conditions ARPACK's LR
            # iteration; we subtract sigma afterwards.
            from scipy.sparse.linalg import LinearOperator as _LinOp
            # spectral_shift is validated (positive-finite / "auto", arnoldi)
            # near the top of _solve_leading.
            if isinstance(spectral_shift, str):  # "auto"
                k_pre = min(6, vec_size - 2)
                if k_pre >= 1:
                    vals_pre, _ = eigs(A, k=k_pre, which='LM')
                    sig = float(np.max(np.abs(vals_pre))) * 1.5 + 1.0e-6
                else:
                    sig = 1.0
                if not np.isfinite(sig):
                    raise ValueError(
                        "auto spectral_shift overflowed to a non-finite value "
                        "(spectral radius too large); pass an explicit shift.")
            else:
                sig = float(spectral_shift)
            logger.info("Spectral shift sigma={:.6f}: eigs(A+sigma*I, LR)".format(sig))
            A_sh = _LinOp(A.shape, matvec=lambda v: A.matvec(v) + sig * v,
                          dtype=A.dtype)
            vals, vecs = eigs(A_sh, k=max_ev, which='LR', v0=seed_vec)
            vals = vals - sig
        else:
            vals, vecs = eigs(A, k=max_ev, which='LM', v0=seed_vec)

    elif solver_mode.startswith("shift-invert"):
        if sigma_shift is None:
            # Estimate shift from a quick Arnoldi run. Sample a few
            # largest-magnitude eigenvalues and aim at the largest *real* part
            # (the physical SC eigenvalue), not the largest magnitude (which
            # can be a large negative repulsive mode).
            k_pre = min(6, vec_size - 2)
            if k_pre < 1:
                # Operator too small for a preliminary ARPACK pass (ARPACK
                # needs k < N-1); fall back to a neutral shift.
                sigma_shift = 0.0
            else:
                logger.info("Estimating shift with preliminary Arnoldi...")
                vals_pre, _ = eigs(A, k=k_pre, which='LM')
                sigma_shift = _shift_from_eigenvalues(vals_pre)
            logger.info("Using sigma_shift = {:.6f}".format(sigma_shift))
        vals, vecs = _eigs_shift_invert(
            A, vec_size, max_ev, solver_mode, sigma=sigma_shift,
            seed_vec=seed_vec
        )

    # With a seed eigenvector, track the branch that overlaps it (eigenvector
    # continuation); otherwise order by largest real part (the physical SC
    # eigenvalue), not magnitude.
    if seed_vec is not None:
        vals, vecs = _order_by_seed_overlap(vals, vecs, seed_vec)
    else:
        vals, vecs = _order_eigenpairs(vals, vecs)

    # A negative leading eigenvalue from plain which='LM' (no spectral_shift) is
    # often an artifact: a small positive lambda masked by larger repulsive
    # modes. Tip the user toward spectral_shift='auto'. Skip this when a
    # seed_vec is given -- there vals[0] is the seed-overlap continuation
    # branch, which may be intentionally negative, not a masked leading mode.
    # Require a meaningfully negative value (relative to the spectral scale)
    # so roundoff-scale negatives near a numerically-zero leading eigenvalue
    # do not trigger a misleading recommendation.
    if (solver_mode == "arnoldi" and spectral_shift is None
            and seed_vec is None and len(vals)):
        scale = float(np.max(np.abs(vals))) if len(vals) else 0.0
        neg_tol = 1.0e-8 * max(scale, 1.0)
        if vals[0].real < -neg_tol:
            logger.warning(
                "Leading eigenvalue Re(lambda)=%.4g is negative; if a positive "
                "(attractive) mode is expected, set [eliashberg] spectral_shift="
                "\"auto\" so the largest-REAL eigenvalue is selected instead of "
                "the largest-magnitude one.", vals[0].real)

    return vals[0], vecs[:, 0], {"eigenvalues": vals, "eigenvectors": vecs,
                                 "sigma_shift": sigma_shift}


def _solve_eigenvalue(Vs_q, G2, norb, Nx, Ny, Nz, num_eigenvalues=10,
                      method="arnoldi", sigma_shift=None, spectral_shift=None):
    """Solve linearized Eliashberg equation by eigenvalue analysis.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    method : str
        Eigenvalue solver method:
        - "arnoldi" : Implicitly Restarted Arnoldi (ARPACK eigs).
            Standard Krylov subspace method for non-symmetric matrices.
            Directly finds the largest eigenvalues.
        - "shift-invert-bicgstab" : Shift-invert with BiCGSTAB.
            Solves (K - sigma*I)^{-1} eigenvalue problem using BiCGSTAB
            for the linear system. Efficient for finding eigenvalues near
            a target value sigma.
        - "shift-invert-gmres" : Shift-invert with GMRES.
            Same as above but uses GMRES for the linear system.
            More robust than BiCGSTAB for non-symmetric systems.
        - "shift-invert-lgmres" : Shift-invert with LGMRES.
            Uses LGMRES (loose GMRES) which can be faster for some problems.
        - "subspace" : Subspace iteration (block power method).
            Simultaneously propagates multiple vectors through the kernel
            with QR orthogonalization. Robust for degenerate eigenvalues
            and finds multiple eigenvalues with different symmetries.
    sigma_shift : float, optional
        Shift parameter for shift-invert methods. Eigenvalues near this
        value are found most efficiently. If None, a preliminary Arnoldi
        step estimates an appropriate shift value.

    Returns
    -------
    eigenvalues : ndarray
        Leading eigenvalues ordered by descending real part (largest first).
    eigenvectors : ndarray
        Corresponding eigenvectors reshaped to (num_ev, norb, norb, Nx, Ny, Nz).
    """
    vec_size = norb * norb * Nx * Ny * Nz

    if method == "subspace":
        # spectral_shift routes through the ARPACK largest-real path in
        # _solve_leading; the subspace driver never reaches it, so reject a
        # non-None spectral_shift here (same arnoldi-only contract) rather
        # than silently ignoring a misconfigured static input.
        if spectral_shift is not None:
            raise ValueError(
                "[eliashberg] spectral_shift is only supported for "
                "eigenvalue_method='arnoldi', not '{}'".format(method))
        # Subspace (block power) iteration has its own dedicated driver
        # (magnitude-based Ritz selection, not the ARPACK/shift-invert path),
        # so it is not routed through _solve_leading; call it directly, exactly
        # as before.
        max_ev = min(num_eigenvalues, vec_size - 2)
        if max_ev < 1:
            max_ev = 1
        logger.info("Computing {} eigenvalues with method='{}'...".format(
            max_ev, method))
        return _solve_subspace_iteration(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=max_ev
        )

    # ARPACK Arnoldi / shift-invert: delegate the eigen-selection, shift
    # estimation, and descending-real-part ordering to the shared driver.
    make_operator = lambda: _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    _, _, eig_analysis = _solve_leading(
        make_operator, vec_size, method,
        num_eigenvalues=num_eigenvalues, sigma_shift=sigma_shift,
        spectral_shift=spectral_shift,
    )
    vals = eig_analysis["eigenvalues"]
    vecs = eig_analysis["eigenvectors"]

    eigenvectors = np.array([
        vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(len(vals))
    ])

    return vals, eigenvectors


def _eigs_shift_invert(A, vec_size, num_ev, method, sigma=0.0, rtol_linear=1e-8,
                       seed_vec=None):
    """Eigenvalue computation using shift-invert with iterative linear solver.

    Transforms the eigenvalue problem K*x = lambda*x into
    (K - sigma*I)^{-1}*x = nu*x where nu = 1/(lambda - sigma).
    Eigenvalues of the original problem near sigma become the largest
    eigenvalues of the transformed problem, which Arnoldi finds efficiently.

    The inverse (K - sigma*I)^{-1} is never formed explicitly; instead,
    each application solves (K - sigma*I)*y = x using an iterative solver.

    Parameters
    ----------
    A : LinearOperator
        The Eliashberg kernel operator K.
    vec_size : int
        Dimension of the vector space.
    num_ev : int
        Number of eigenvalues to compute.
    method : str
        One of "shift-invert-bicgstab", "shift-invert-gmres",
        "shift-invert-lgmres".
    sigma : float
        Shift value. Eigenvalues near sigma are found most efficiently.
        Default 0.0 finds eigenvalues nearest to zero.
    rtol_linear : float
        Relative tolerance for the iterative linear solver.

    Returns
    -------
    eigenvalues : ndarray
        Eigenvalues of the original problem K.
    eigenvectors : ndarray
        Corresponding eigenvectors.
    """
    solver_name = method.replace("shift-invert-", "")
    solver_map = {
        "bicgstab": bicgstab,
        "gmres": gmres,
        "lgmres": lgmres,
    }
    if solver_name not in solver_map:
        raise ValueError("Unknown linear solver: {}".format(solver_name))
    linear_solver = solver_map[solver_name]

    # (A - sigma*I)
    def shifted_matvec(v):
        return A.matvec(v) - sigma * v

    A_shifted = LinearOperator((vec_size, vec_size),
                               matvec=shifted_matvec, dtype=complex)

    solve_count = [0]
    fail_count = [0]

    # (A - sigma*I)^{-1} via iterative solver
    def inv_matvec(v):
        solve_count[0] += 1
        x, info = linear_solver(A_shifted, v, rtol=rtol_linear, maxiter=500)
        if info != 0:
            fail_count[0] += 1
            if fail_count[0] <= 3:
                logger.warning(
                    "Linear solver {} did not converge (info={}), "
                    "solve #{}, using approximate solution".format(
                        solver_name, info, solve_count[0]))
        return x

    A_inv = LinearOperator((vec_size, vec_size),
                           matvec=inv_matvec, dtype=complex)

    # eigs on (A - sigma*I)^{-1} finds eigenvalues nu = 1/(lambda - sigma)
    # largest |nu| correspond to lambda closest to sigma. A seed vector (v0)
    # biases the Arnoldi start toward the physical branch being tracked.
    nus, vecs = eigs(A_inv, k=num_ev, which='LM', v0=seed_vec)

    logger.info("Shift-invert: {} linear solves, {} failures".format(
        solve_count[0], fail_count[0]))

    # Convert back: lambda = 1/nu + sigma
    eigenvalues = 1.0 / nus + sigma

    return eigenvalues, vecs


def _solve_subspace_iteration(Vs_q, G2, norb, Nx, Ny, Nz,
                              num_eigenvalues=5, max_iter=300, tol=1e-6):
    """Find multiple eigenvalues by subspace iteration (block power method).

    Simultaneously propagates a block of vectors through the kernel,
    orthogonalizing via QR decomposition at each step. This naturally
    converges to the dominant invariant subspace.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    num_eigenvalues : int
        Number of eigenvalues to compute.
    max_iter : int
        Maximum iterations.
    tol : float
        Convergence tolerance for eigenvalues.

    Returns
    -------
    eigenvalues : ndarray
        Converged eigenvalues ordered by descending magnitude.
    eigenvectors : ndarray
        Shape (num_eigenvalues, norb, norb, Nx, Ny, Nz).
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)
    num_ev = min(num_eigenvalues, vec_size)

    # Use extra vectors for better convergence (subspace padding)
    n_work = min(num_ev + max(num_ev, 5), vec_size)

    # Random initial subspace
    rng = np.random.default_rng(42)
    V = rng.standard_normal((vec_size, n_work))
    V, _ = np.linalg.qr(V)

    eigenvalues_old = np.zeros(num_ev, dtype=complex)

    for iteration in range(max_iter):
        # Apply kernel to all vectors at once: W = A @ V (kernel may be
        # complex). matmat batches all columns into a single FFT (numerically
        # identical to the per-column matvec).
        W = A.matmat(V)

        # Rayleigh quotient: H = V^dagger A V (conjugate transpose for complex)
        H = V.conj().T @ W

        # Eigendecomposition of small matrix
        evals_h, evecs_h = np.linalg.eig(H)
        idx = np.argsort(-np.abs(evals_h))
        evals_h = evals_h[idx]
        evecs_h = evecs_h[:, idx]

        # Ritz vectors: update subspace
        V = W @ evecs_h
        # Re-orthogonalize
        V, _ = np.linalg.qr(V)

        # Check convergence of the wanted eigenvalues (track the full complex
        # value so imaginary motion is not ignored for a complex kernel).
        eigenvalues_new = evals_h[:num_ev]
        diff = np.max(np.abs(eigenvalues_new - eigenvalues_old))

        if iteration % 10 == 0 or diff < tol:
            logger.info("Subspace iter {:4d}: eigenvalues = {}, diff = {:.2e}".format(
                iteration, np.array2string(eigenvalues_new, precision=4), diff))

        if diff < tol:
            logger.info("Subspace iteration converged at iteration {}".format(
                iteration + 1))
            break

        eigenvalues_old = eigenvalues_new.copy()

    # Extract final Ritz vectors for the wanted eigenvalues
    # Recompute from final V (batched matmat = per-column matvec)
    W = A.matmat(V)
    H = V.conj().T @ W
    evals_h, evecs_h = np.linalg.eig(H)
    # Subspace (block power) iteration is magnitude-based: it converges to the
    # dominant-magnitude invariant subspace, so report its modes by magnitude.
    # The physical SC selection (largest real) is applied in the default
    # arnoldi path; use that to read off the leading SC eigenvalue.
    idx = np.argsort(-np.abs(evals_h))

    eigenvalues = evals_h[idx[:num_ev]]
    ritz_vecs = V @ evecs_h[:, idx[:num_ev]]

    eigenvectors = np.array([
        ritz_vecs[:, i].reshape(norb, norb, Nx, Ny, Nz)
        for i in range(num_ev)
    ])

    return eigenvalues, eigenvectors


def _solve_shifted_bicg(Vs_q, G2, norb, Nx, Ny, Nz,
                        sigma_list, num_eigenvalues=3, tol_linear=1e-8):
    """Find eigenvalues near multiple target values using shifted BiCG.

    Parameters
    ----------
    Vs_q : ndarray
        Effective pairing interaction. Either shape (norb, norb, Nx, Ny, Nz)
        for simple mode or (norb, norb, norb, norb, Nx, Ny, Nz) for general.
    G2 : ndarray
        Two-particle Green's function.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        System dimensions.
    sigma_list : list of float
        Target shift values. Eigenvalues near each sigma are found.
    num_eigenvalues : int
        Number of eigenvalues to find near each sigma.
    tol_linear : float
        Tolerance for the linear solver.

    Returns
    -------
    all_eigenvalues : dict
        Dictionary mapping sigma -> array of eigenvalues.
    all_eigenvectors : dict
        Dictionary mapping sigma -> array of eigenvectors.
    """
    A, vec_size = _make_kernel_operator(Vs_q, G2, norb, Nx, Ny, Nz)

    max_ev = min(num_eigenvalues, vec_size - 2)
    if max_ev < 1:
        max_ev = 1

    all_eigenvalues = {}
    all_eigenvectors = {}

    # Use the first sigma as seed system
    sigma_seed = sigma_list[0]

    # Seed: build Krylov subspace for (A - sigma_seed I)
    # We use a manual BiCG implementation with shift tracking
    logger.info("Shifted BiCG: seed sigma = {:.6f}, {} additional shifts".format(
        sigma_seed, len(sigma_list) - 1))

    # For each shift, solve via shift-invert + eigs
    # The key optimization: reuse preconditioner from seed system
    for i, sigma in enumerate(sigma_list):
        logger.info("Processing shift sigma = {:.6f} ({}/{})".format(
            sigma, i + 1, len(sigma_list)))

        def shifted_matvec(v, s=sigma):
            return A.matvec(v) - s * v

        A_shifted = LinearOperator(
            (vec_size, vec_size), matvec=shifted_matvec, dtype=complex
        )

        def inv_matvec(v):
            x, info = gmres(A_shifted, v, rtol=tol_linear, maxiter=300)
            return x

        A_inv = LinearOperator(
            (vec_size, vec_size), matvec=inv_matvec, dtype=complex
        )

        try:
            nus, vecs = eigs(A_inv, k=max_ev, which='LM')
            eigenvalues = 1.0 / nus + sigma
            # Eigenvalues found near this shift; report by magnitude.
            idx = np.argsort(-np.abs(eigenvalues))
            eigenvalues = eigenvalues[idx]
            eigvecs = np.array([
                vecs[:, j].reshape(norb, norb, Nx, Ny, Nz)
                for j in idx
            ])
        except Exception as e:
            logger.warning("Shift sigma={:.6f} failed: {}".format(sigma, e))
            eigenvalues = np.array([])
            eigvecs = np.array([])

        all_eigenvalues[sigma] = eigenvalues
        all_eigenvectors[sigma] = eigvecs

        if len(eigenvalues) > 0:
            logger.info("  Found eigenvalues: {}".format(
                np.array2string(eigenvalues.real, precision=6)))

    return all_eigenvalues, all_eigenvectors


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def _save_results(output_dir, sigma, eigenvalue, eigenvalues_eig, kx_array, ky_array, kz_array,
                  gap_file="gap.dat", eigenvalue_file="eigenvalue.dat",
                  eigenvalue_match=None):
    """Save gap function and eigenvalue results to files.

    Parameters
    ----------
    output_dir : str
        Output directory path.
    sigma : ndarray or None
        Gap function from iteration, shape (norb, norb, Nx, Ny, Nz).
    eigenvalue : float or None
        Leading eigenvalue from iteration.
    eigenvalues_eig : ndarray or None
        Eigenvalues from eigenvalue analysis.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    gap_file : str
        Output filename for gap function.
    eigenvalue_file : str
        Output filename for eigenvalues.
    """
    os.makedirs(output_dir, exist_ok=True)

    if sigma is not None:
        norb = sigma.shape[0]
        Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
        filepath = os.path.join(output_dir, gap_file)
        logger.info("Saving gap function to {}".format(filepath))

        with open(filepath, "w") as fw:
            # Header
            header_parts = ["# kx", "ky", "kz"]
            for i in range(norb):
                for j in range(norb):
                    header_parts.append("Re(sigma_{}{})".format(i, j))
                    header_parts.append("Im(sigma_{}{})".format(i, j))
            fw.write(" ".join(header_parts) + "\n")

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
                        parts = ["{:.8f}".format(kx),
                                 "{:.8f}".format(ky),
                                 "{:.8f}".format(kz)]
                        for i in range(norb):
                            for j in range(norb):
                                val = sigma[i, j, ix, iy, iz]
                                parts.append("{:.8e}".format(val.real))
                                parts.append("{:.8e}".format(val.imag))
                        fw.write(" ".join(parts) + "\n")

    if eigenvalue is not None or eigenvalues_eig is not None:
        filepath = os.path.join(output_dir, eigenvalue_file)
        logger.info("Saving eigenvalues to {}".format(filepath))
        with open(filepath, "w") as fw:
            if eigenvalue is not None:
                fw.write("# Iteration eigenvalue\n")
                fw.write("{:.8e}\n".format(eigenvalue))
            if eigenvalues_eig is not None:
                fw.write("# Eigenvalue analysis\n")
                # The optional trailing `match` column is 1 when the eigenvector
                # lies (dominantly) in the channel's parity sector (even for
                # singlet, odd for triplet) and 0 when it lies in the opposite
                # sector. The leading index/Re/Im/|ev| columns are unchanged, so
                # position-based parsers (e.g. reading column 1 as Re) are
                # unaffected; only the column count grows from 4 to 5.
                if eigenvalue_match is not None:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|  match(1=channel-parity)\n")
                else:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|\n")
                for i, ev in enumerate(eigenvalues_eig):
                    row = "{:4d} {:15.8e} {:15.8e} {:15.8e}".format(
                        i, ev.real, ev.imag, abs(ev))
                    if eigenvalue_match is not None:
                        row += " {:d}".format(int(bool(eigenvalue_match[i])))
                    fw.write(row + "\n")


# ---------------------------------------------------------------------------
# chi0q format conversion
# ---------------------------------------------------------------------------

def _convert_chi0q_to_ref_format(chi0q, norb, Nx, Ny, Nz):
    """Convert chi0q from H-wave format to reference code format.

    Supports both 2-index (reduced) and 4-index (general) chi0q:

    2-index:
        H-wave: (nmat, nvol, norb, norb) -> Ref: (norb, norb, Nx, Ny, Nz, nmat)
    4-index:
        H-wave: (nmat, nvol, norb, norb, norb, norb)
        -> Ref: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)

    Parameters
    ----------
    chi0q : ndarray
        chi0q in H-wave format.
    norb, Nx, Ny, Nz : int
        System parameters.

    Returns
    -------
    chi0q_ref : ndarray
        chi0q in reference format:
        - 2-index: (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    """
    if chi0q.ndim == 4:
        # 2-index: (nmat, nvol, norb, norb) -> (norb, norb, Nx, Ny, Nz, nmat)
        nf = chi0q.shape[0]
        chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb)
        chi0q_ref = chi0q_3d.transpose(4, 5, 1, 2, 3, 0)
    elif chi0q.ndim == 6:
        if chi0q.shape[0] == norb and chi0q.shape[1] == norb:
            # Already in ref format (norb, norb, Nx, Ny, Nz, nmat)
            chi0q_ref = chi0q
        else:
            # 4-index H-wave: (nmat, nvol, norb, norb, norb, norb)
            # -> (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
            nf = chi0q.shape[0]
            chi0q_3d = chi0q.reshape(nf, Nx, Ny, Nz, norb, norb, norb, norb)
            chi0q_ref = chi0q_3d.transpose(4, 5, 6, 7, 1, 2, 3, 0)
    elif chi0q.ndim == 8:
        # Already in ref format (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
        chi0q_ref = chi0q
    else:
        raise ValueError(
            "Unexpected chi0q shape: {}. Expected 4D, 6D, or 8D.".format(chi0q.shape)
        )
    return chi0q_ref


# ---------------------------------------------------------------------------
# Main calculation
# ---------------------------------------------------------------------------

def calc_eliashberg(input_dict):
    """Main calculation orchestration for linearized Eliashberg equation.

    Supports ``[eliashberg] gpu = true`` for both ``frequency = "static"`` and
    ``"dynamic"``: the kernel matvec/matmat runs on the GPU (CuPy) while the
    host eigensolvers (ARPACK/power/subspace/shift-invert) consume host arrays.
    Falls back to CPU with a warning when CuPy/CUDA is unusable, unless
    ``gpu_required = true`` (then it fails fast).

    Parameters
    ----------
    input_dict : dict
        Parsed TOML configuration dictionary.
    """
    # --- Parse parameters ---
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = 1.0 / T
    cell_shape = mode_param["CellShape"]
    # Resolve SubShape by the PACKAGE convention (documented default:
    # CellShape, i.e. the whole cell as one supercell) and guard the
    # RESOLVED value: a guard keyed on the explicit key alone let an
    # omitted SubShape default to a fully folded configuration and reach
    # the file reader with the very mismatch the guard exists to stop
    # (round-7 review).
    if isinstance(cell_shape, (list, tuple)):
        _cs = list(cell_shape)
    else:
        _cs = [cell_shape]
    while len(_cs) < 3:
        _cs.append(1)
    _ss = list(mode_param.get("SubShape", _cs))
    while len(_ss) < 3:
        _ss.append(1)
    if _ss != [1, 1, 1]:
        # supported nowhere in this module: the geometry and
        # interactions are consumed UNFOLDED here, so a folded
        # susceptibility mismatches the expected orbital count and an
        # off-site bond would fold onto an on-site supercell entry
        # (round-4 review); failing late produced an unhelpful shape
        # error instead of this actionable one
        raise ValueError(
            "SubShape (sublattice folding) is not supported by the "
            "Eliashberg module: fold the model into the unit cell "
            "yourself, or set SubShape = [1, 1, 1] explicitly. (Note: "
            "omitting SubShape defaults it to CellShape -- the whole "
            "cell as one supercell -- per the package convention, so it "
            "must be set explicitly here.)")
    sub_shape = _ss
    nmat = mode_param.get("Nmat", _DEFAULT_NMAT)

    # Filling
    if "filling" in mode_param:
        n_filling = mode_param["filling"]
    elif "Ncond" in mode_param:
        raise ValueError("Ncond not yet supported in eliashberg. Use filling.")
    else:
        raise ValueError("filling must be specified.")

    # Lattice dimensions
    if isinstance(cell_shape, list):
        while len(cell_shape) < 3:
            cell_shape.append(1)
    Lx, Ly, Lz = cell_shape
    if isinstance(sub_shape, list):
        while len(sub_shape) < 3:
            sub_shape.append(1)
    Bx, By, Bz = sub_shape
    Nx, Ny, Nz = Lx // Bx, Ly // By, Lz // Bz
    nvol = Nx * Ny * Nz

    # Eliashberg parameters
    eli_param = input_dict.get("eliashberg", {})

    # Dispatch to dynamic Eliashberg if requested
    if _eliashberg_frequency(input_dict) == "dynamic":
        _validate_dynamic_prereqs(input_dict)
        from hwave.solver import eliashberg_dynamic
        return eliashberg_dynamic.solve_dynamic(input_dict)

    # Static path from here on: the channel-decomposition flags only affect the
    # dynamic vertex, so warn rather than silently ignore them.
    _warn_if_static_ignores_channel_flags(eli_param)

    # GPU backend for the static kernel (the dynamic branch returned above and
    # handles its own gpu flag). get_backend falls back to numpy with a warning
    # when CuPy/CUDA is unusable, unless gpu_required=true (then it raises).
    use_gpu = backend.as_bool(eli_param.get("gpu", False))
    gpu_required = backend.as_bool(eli_param.get("gpu_required", False))
    xp, gpu_active = backend.get_backend(use_gpu, logger, required=gpu_required)

    solver_mode = eli_param.get("solver_mode", "iteration")
    max_iter = eli_param.get("max_iter", 1000)
    alpha = eli_param.get("alpha", 0.5)
    tol = eli_param.get("convergence_tol", 1.0e-5)
    num_eigenvalues = eli_param.get("num_eigenvalues", 10)
    eigenvalue_method = eli_param.get("eigenvalue_method", "arnoldi")
    pairing_type = eli_param.get("pairing_type", "singlet")
    # Default the initial gap to the channel's parity (even for singlet, odd
    # for triplet) when the user did not specify one. For a centrosymmetric
    # system the kernel preserves parity in exact arithmetic, but numerical
    # noise leaks the wrong parity; _solve_iteration projects every iterate back
    # onto the channel's sector (pairing_type below), so the seed stays in -- and
    # converges within -- the physical sector. (For a non-centrosymmetric kernel
    # _solve_iteration detects [A,P] != 0 and disables the projection.)
    init_gap_mode = _resolve_init_gap(eli_param.get("init_gap"), pairing_type)
    chi0q_mode = eli_param.get("chi0q_mode", "load")
    chi0q_tensor = eli_param.get("chi0q_tensor", "auto")
    gap_file = eli_param.get("output_gap", "gap.dat")
    eigenvalue_file = eli_param.get("output_eigenvalue", "eigenvalue.dat")

    output_dir = input_dict["file"]["output"]["path_to_output"]

    logger.info("=== Linearized Eliashberg Equation Solver ===")
    logger.info("T = {}, beta = {:.4f}".format(T, beta))
    logger.info("Cell: {}x{}x{}, Sub: {}x{}x{}, Grid: {}x{}x{}".format(
        Lx, Ly, Lz, Bx, By, Bz, Nx, Ny, Nz))
    logger.info("Nmat = {}, filling = {}".format(nmat, n_filling))
    logger.info("solver_mode = {}, chi0q_mode = {}".format(solver_mode, chi0q_mode))

    # --- Step 2: Read input files ---
    geom_info, hr, interactions = _read_interaction_files(input_dict)
    norb = geom_info["norb"]

    # --- Step 3: Setup k-mesh ---
    kx_array = np.linspace(0, 2.0 * np.pi, Nx, endpoint=False)
    ky_array = np.linspace(0, 2.0 * np.pi, Ny, endpoint=False)
    kz_array = np.linspace(0, 2.0 * np.pi, Nz, endpoint=False)

    # --- Step 4: Build Hamiltonian and diagonalize ---
    logger.info("Building Hamiltonian in k-space...")
    epsilon_k = _build_hamiltonian_k(kx_array, ky_array, kz_array, hr, norb)

    logger.info("Diagonalizing Hamiltonian...")
    eigenvalues, eigenvectors = _calc_eigenvalues(epsilon_k)

    # --- Step 5: Determine chemical potential ---
    logger.info("Determining chemical potential...")
    mu = _determine_mu(eigenvalues, beta, n_filling, norb)
    logger.info("mu = {:.6f}".format(mu))

    # --- Step 7: Build interaction in k-space ---
    logger.info("Building interactions in k-space...")
    inter_k = _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb)

    if chi0q_mode == "flex":
        # --- FLEX mode: use pre-computed dressed susceptibilities ---
        logger.info("FLEX mode: loading dressed susceptibilities and Green's function")

        chis, chic, green_dressed, chi_convention = _load_flex_susceptibilities(
            input_dict, norb, Nx, Ny, Nz, interactions=interactions)
        logger.info("FLEX susceptibility convention: {}".format(chi_convention))

        # Use dressed Green's function if available, otherwise use bare
        if green_dressed is not None:
            green_kw = green_dressed
            logger.info("Using FLEX dressed Green's function")
        else:
            logger.info("Calculating bare Green's function G(k, iwn)...")
            green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        # Compute pairing vertex from FLEX susceptibilities
        logger.info("Computing FLEX vertices (pairing_type={})...".format(pairing_type))
        # Reject BEFORE the S/C build: the pair is O(Nq * norb^4) and an
        # unsupported calculation must not allocate heavily on its way to
        # the validation error (round-10 review).
        _reject_reduced_flex_unsupported(inter_k, chi_convention)
        # One S/C build for the run: the pair feeds the missing-component
        # diagnostic AND the vertex contraction (round-9 review; previously
        # each built its own full-grid pair).
        sc_mats = _build_vertex_sc_matrices(chi_convention, inter_k,
                                            norb, Nx, Ny, Nz)
        _warn_reduced_flex_missing_components(
            inter_k, norb, Nx, Ny, Nz, chi_convention,
            sc_matrices=(sc_mats if str(chi_convention).lower() == "kuroki"
                         else None))
        Vs_q = _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                                      pairing_type=pairing_type,
                                      convention=chi_convention,
                                      sc_matrices=sc_mats)
        logger.info("FLEX vertex shape: {}".format(Vs_q.shape))

    else:
        # --- Standard RPA mode ---

        # Step 6: Calculate bare Green's function
        logger.info("Calculating Green's function G(k, iwn)...")
        green_kw = _calc_green(eigenvalues, eigenvectors, mu, beta, nmat)

        # Step 1: Load or compute chi0q
        if chi0q_mode == "calc":
            chi0q_raw = _calc_chi0q_internal(input_dict, chi0q_tensor=chi0q_tensor,
                                                precomputed_mu=mu)
            # internally computed chi0q always carries the full frequency grid
            static_index = None
        else:
            chi0q_raw, static_index = _load_chi0q(input_dict)

        # Step 8: Convert chi0q format; the frequency axis is last in the
        # reference format, so read its length after the conversion (a 6D
        # input can be either raw H-wave (nmat first) or ref (nmat last))
        chi0q = _convert_chi0q_to_ref_format(chi0q_raw, norb, Nx, Ny, Nz)
        nmat_chi0q = chi0q.shape[-1]
        logger.info("chi0q converted to shape: {}".format(chi0q.shape))

        # Step 9: Compute RPA vertices
        logger.info("Computing RPA vertices (pairing_type={})...".format(pairing_type))
        vertex_result = _compute_vertices(
            chi0q, inter_k, norb, Nx, Ny, Nz, nmat_chi0q,
            pairing_type=pairing_type, static_index=static_index,
            declarations_closed=_declarations_partner_closed(
                interactions, Nx, Ny, Nz, norb))
        if isinstance(vertex_result, tuple):
            Pc_q, Ps_q = vertex_result
            Vs_q = Pc_q + Ps_q
            logger.info("Simple mode: Pc + Ps vertex, shape {}".format(Vs_q.shape))
        else:
            Vs_q = vertex_result
            logger.info("General mode: 4-index V^s vertex, shape {}".format(Vs_q.shape))

    # --- Step 10: Compute G2 ---
    logger.info("Computing G2...")
    # New config keys are read case-insensitively (the config layers disagree
    # on case handling; see the PR #128 sweep) and through the strict parser:
    # a misspelled value must fail, not silently flip the physics.
    from requests.structures import CaseInsensitiveDict
    g2_tail = _coerce_g2_tail(
        CaseInsensitiveDict(eli_param).get("g2_tail", True))
    logger.info("g2_tail = %s (Matsubara tail correction for the pair "
                "bubble); Green frequency axis = %d points", g2_tail,
                green_kw.shape[-1])
    if g2_tail:
        _warn_if_g2_tail_outside_asymptotic_regime(green_kw, beta)
    G2 = _calc_g2(green_kw, beta, tail=g2_tail)
    _warn_if_g2_indefinite(G2, norb, g2_tail)

    # --- Step 11: Initialize gap function ---
    sigma_init = _initialize_gap(init_gap_mode, norb, kx_array, ky_array, kz_array)

    # GPU: park the two large invariants (pairing vertex and pair bubble) on the
    # device once; each matvec then only moves the gap vector across PCIe. The
    # solver entry points below pass Vs_q/G2 straight to _make_kernel_operator,
    # which derives its backend from them, so no other change is needed here.
    if gpu_active:
        # Resident device tensors are more than the two inputs: the operator
        # also holds V_r (~Vs_q) and G2_pre (~G2), and in general mode
        # G2_gemm (~G2) + V_r_gemm (~Vs_q) -- i.e. up to ~3x the inputs. The
        # per-matmat block workspace (gap-sized x num columns) adds on top; a
        # subspace/eigenvalue run with a wide block can still exceed this. It is
        # only a pre-warning -- CuPy raises a clear OutOfMemoryError regardless.
        est_bytes = 3 * (Vs_q.nbytes + G2.nbytes)
        backend.warn_if_device_memory_short(
            est_bytes, logger, label="the static Eliashberg kernel")
        logger.info("GPU backend active (CuPy): moving the pairing vertex and "
                    "G2 to the device (%.2f GB).",
                    (Vs_q.nbytes + G2.nbytes) / 1e9)
        Vs_q = xp.asarray(Vs_q)
        G2 = xp.asarray(G2)

    # --- Step 12: Solve ---
    sigma_result = None
    eigenvalue_iter = None
    eigenvalues_eig = None
    eigenvalue_match = None

    if solver_mode in ("iteration", "both"):
        logger.info("=== Self-consistent iteration ===")
        sigma_result, eigenvalue_iter, converged, n_iter = _solve_iteration(
            green_kw, Vs_q, G2, sigma_init, norb,
            max_iter=max_iter, alpha=alpha, tol=tol, pairing_type=pairing_type
        )
        logger.info("Iteration result: eigenvalue = {:.6f}, converged = {}, n_iter = {}".format(
            eigenvalue_iter, converged, n_iter))

    if solver_mode in ("eigenvalue", "both"):
        logger.info("=== Eigenvalue analysis ===")
        eigenvalues_eig, eigenvectors_eig = _solve_eigenvalue(
            Vs_q, G2, norb, Nx, Ny, Nz,
            num_eigenvalues=num_eigenvalues,
            method=eigenvalue_method,
            sigma_shift=eli_param.get("sigma_shift"),
            spectral_shift=eli_param.get("spectral_shift"),
        )
        # The kernel preserves parity; promote the eigenpairs whose gap has the
        # requested channel parity (singlet even / triplet odd) so the reported
        # leading eigenpair is the physical solution for this channel.
        eigenvalues_eig, eigenvectors_eig = _reorder_eigenpairs_by_parity(
            eigenvalues_eig, eigenvectors_eig, pairing_type)
        # Flag whether each eigenvector lies (dominantly) in the channel's
        # parity sector (even for singlet, odd for triplet). _is_gap_parity is a
        # threshold test, so this is "dominantly the channel parity" rather than
        # an exact label; for a centrosymmetric model parity is an exact kernel
        # symmetry and the opposite-parity modes are unphysical for the channel
        # (e.g. an even mode in the triplet kernel is Pauli-forbidden).
        eigenvalue_match = np.array(
            [_is_gap_parity(g, pairing_type) for g in eigenvectors_eig])
        logger.info("Leading eigenvalues:")
        for i, ev in enumerate(eigenvalues_eig):
            tag = "" if eigenvalue_match[i] else "  [opposite-parity sector]"
            logger.info("  {:3d}: {:.6f} (|ev| = {:.6f}){}".format(
                i, ev.real, abs(ev), tag))

        # Use leading eigenvector as gap if no iteration result
        if sigma_result is None:
            sigma_result = eigenvectors_eig[0]
            eigenvalue_iter = eigenvalues_eig[0].real

    # --- Step 13: Save results ---
    _save_results(
        output_dir, sigma_result, eigenvalue_iter, eigenvalues_eig,
        kx_array, ky_array, kz_array,
        gap_file=gap_file, eigenvalue_file=eigenvalue_file,
        eigenvalue_match=eigenvalue_match
    )

    logger.info("=== Done ===")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    """Command-line interface for the Eliashberg equation solver."""
    import tomli
    import argparse

    parser = argparse.ArgumentParser(
        description="Linearized Eliashberg equation solver (post-tool for H-wave RPA)",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument("input", type=str, help="input TOML file (same as RPA input + [eliashberg] section)")
    parser.add_argument("-q", "--quiet", action="store_true", help="suppress output")
    parser.add_argument(
        "-v", "--version", action="version",
        version="hwave_sc v{}".format(hwave.__version__),
    )
    args = parser.parse_args()

    # Setup logging
    log_level = logging.WARNING if args.quiet else logging.INFO
    logging.basicConfig(level=log_level, format="%(name)s: %(message)s")

    file_toml = args.input
    if not os.path.exists(file_toml):
        logger.error("Input file does not exist: {}".format(file_toml))
        sys.exit(1)

    with open(file_toml, "rb") as f:
        input_dict = tomli.load(f)

    calc_eliashberg(input_dict)


if __name__ == "__main__":
    main()
