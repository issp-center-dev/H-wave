"""Reading and checking the archives of hwave_sc: chi0q, FLEX susceptibilities and Green functions, interaction files, and the chi0q internal computation.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 3); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import os
import zipfile
from typing import NamedTuple
import numpy as np

from hwave.solver import npy_header as _npy_header
import hwave.qlmsio.wan90 as wan90
from hwave.solver.rpa import validate_chi0q_index_convention, TAIL_ENDPOINT_CONVENTION
from hwave.solver.ir_axis import is_ir_native, ir_native_meta
from hwave.solver.eliashberg.vertex import placement_sensitive_types, check_placement_stamp
from hwave.solver.eliashberg.config import _reject_ir_native, _check_spin_block_discarded, _finite_positive_float64, _coerce_run_beta, _resolve_spin_orbital_flag
from hwave.solver.eliashberg.band import _MeanFieldBand, _checked_band, _check_green_init_spin_free, SPIN_FREE_ATOL, _ONE_BODY_FILE_KEYS

logger = logging.getLogger("hwave_sc")


#: IR-native refusal hint of the spin-orbital dynamic reader (issue #284)
_SO_DYNAMIC_IR_HINT = (
    "the spin-orbital dynamic Eliashberg solver reads uniform-grid files only "
    "(the IR basis is not supported in spin-orbital mode yet). Re-run FLEX "
    "with [mode.param] write_densified = true.")


_STATIC_IR_HINT = (
    "the static Eliashberg solver requires uniform-grid files. Re-run FLEX "
    "with [mode.param] write_densified = true, or switch to frequency = "
    "\"dynamic\" with [eliashberg] matsubara_basis = \"ir\".")


#: ``chi_convention`` of the single antisymmetrized spinful channel a
#: spin-orbital FLEX run stores as ``chiq`` (flex.py save_results).
_SO_FLEX_CHI_CONVENTION = "spin_orbital"


_SO_FLEX_PRODUCER = "the spin-orbital FLEX run (hwave FLEX, enable_spin_orbital = true)"


# Default Matsubara grid size (mode.param.Nmat). Single definition so the
# legacy-file fallback in _static_freq_position and the main solver loop can
# never silently disagree on the grid.
_DEFAULT_NMAT = 1024


#: the producers a spin-orbital hwave_sc run names when an archive it reads
#: must be regenerated: the RPA run of a chi0q load, the FLEX run of a flex
#: entry (static and dynamic)
_SO_RPA_PRODUCER = "the spin-orbital RPA run (hwave RPA, enable_spin_orbital = true)"


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
        The nmat metadata from the npz file, if present.  When present it
        is the authoritative full-grid reference everywhere in this
        function -- including the two metadata-poor fallbacks, where a
        stored axis equal to the producer's own nmat is provably a full
        grid even if the consuming run's Nmat differs.

    Returns
    -------
    int or None
        Position of the zero bosonic frequency along the stored axis, or
        None when the file carries no usable metadata -- the caller must
        then slice the CENTER of the frequency axis it actually uses (only
        the caller can identify that axis reliably, e.g. for the 6D
        reference format).
    """
    # The full-grid reference: the PRODUCER's own recorded grid size when the
    # file carries one (authoritative provenance -- a stored axis equal to it
    # is provably a full grid whatever the consuming run's Nmat is), else the
    # configured Nmat.
    ref_nmat = int(file_nmat) if file_nmat is not None else config_nmat
    ref_desc = ("the file's recorded nmat = {}".format(ref_nmat)
                if file_nmat is not None
                else "mode.param.Nmat = {}".format(ref_nmat))
    # the same reference phrased as a requirement, for the message that
    # states what nfreq would have to equal
    ref_req = ("the file's recorded nmat ({})".format(ref_nmat)
               if file_nmat is not None
               else "mode.param.Nmat ({})".format(ref_nmat))
    if file_nmat is not None:
        ref_advice = (
            "The file records nmat = {}, so its stored axis is a "
            "restriction of that grid: regenerate it holding the full "
            "grid, or with a matsubara_frequency range covering nmat//2 "
            "(which records freq_index).".format(ref_nmat))
    else:
        ref_advice = (
            "If the file holds a full grid, set mode.param.Nmat = {}; "
            "otherwise regenerate it with a newer version (which records "
            "nmat in the file).".format(nfreq))

    if freq_index is None:
        # Pre-provenance file: with no metadata at all, the stored axis can
        # only be centered safely when the choice of slice is unambiguous.
        # A single-frequency axis (nfreq == 1) is such a case: index 0 is the
        # ONLY slice, so centering cannot pick a different (finite) frequency
        # -- this is the common static-only chi0q reduced to its zero
        # bosonic-frequency component.  A full grid of the reference size
        # (nfreq == ref_nmat) is the other.  Any other length could be a
        # matsubara_frequency restriction of some run's grid, so nfreq//2
        # might point at a finite frequency -- refuse to guess.
        if nfreq == 1:
            # Ruling (#186 review): a static-only chi0q holds exactly the
            # zero bosonic frequency by construction, so keep accepting the
            # singleton -- but state the assumption instead of reporting a
            # generic centering.
            logger.warning(
                "chi0q file '{}' has no freq_index metadata and a single "
                "stored frequency slice; assuming it is the zero bosonic "
                "frequency (a static-only file) -- regenerate the file "
                "with a newer version to record freq_index/nmat if it is "
                "not.".format(file_name))
            return None
        if nfreq == ref_nmat:
            logger.warning(
                "chi0q file '{}' has no freq_index metadata; using the "
                "center of the stored frequency axis as the static slice."
                .format(file_name))
            return None
        raise ValueError(
            "chi0q file '{}' has no freq_index metadata and its stored "
            "frequency axis length {} does not match {}, so the "
            "zero-frequency position cannot be determined: centering the "
            "axis could pick a finite frequency of a restricted grid. {}"
            .format(file_name, nfreq, ref_desc, ref_advice))

    freq_index = np.asarray(freq_index).ravel()
    if freq_index.size == 0:
        raise ValueError(
            "chi0q file '{}' records no frequencies (matsubara_frequency="
            "\"none\"?); it cannot provide the static susceptibility."
            .format(file_name))
    if freq_index.size != nfreq:
        # Legacy FLEX files store the FULL grid but a restricted freq_index
        # (FLEX never applied the matsubara_frequency filter): the DATA axis
        # is authoritative, so fall back to the pre-metadata behavior -- but
        # ONLY when the data axis is provably a full grid of the reference
        # size (nfreq == ref_nmat, the file's own recorded nmat when it has
        # one, else the configured Nmat) AND freq_index is a 0-based subset
        # of it (set(freq_index) <= set(range(nfreq))).  Otherwise the stored
        # axis could be a restriction of some other grid and centering could
        # pick a finite frequency, so refuse to guess.
        if nfreq == ref_nmat and set(freq_index.tolist()) <= set(
                range(nfreq)):
            logger.warning(
                "chi0q file '{}': freq_index length {} does not match the "
                "frequency axis length {}; ignoring the metadata and using "
                "the center of the stored frequency axis as the static "
                "slice.".format(file_name, freq_index.size, nfreq))
            return None
        raise ValueError(
            "chi0q file '{}': freq_index length {} does not match the "
            "frequency axis length {}, and the axis is not provably a full "
            "grid (nfreq must equal {}, with a 0-based freq_index), so the "
            "zero-frequency position cannot be determined: centering the "
            "axis could pick a finite frequency of a restricted grid. {}"
            .format(file_name, freq_index.size, nfreq, ref_req, ref_advice))

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


def _require_so_chi0q_marker(data, file_name, hint=None, producer="RPA"):
    """Refuse a chi0q that does not record ``spin_orbital = True``.

    Every RPA run writes ``index_convention = "spin_block"``, so that tag does
    not identify a spin-orbital file, and the shape cannot either: a
    non-spin-orbital chi0q over 2N orbitals has the shape of a spin-orbital
    one over N physical orbitals. The RPA therefore writes an explicit
    ``spin_orbital = True`` marker, ONLY in a spin-orbital run (a
    non-spin-orbital archive carries no such key); spin-orbital mode accepts
    only a 0-d boolean ``True``. Checked on the open npz before the chi0q
    array itself is read. Also used for the chi0q companion of a
    spin-orbital FLEX run (``producer = "FLEX"``)."""
    if hint is None:
        hint = ("Regenerate the file with a spin-orbital RPA run "
                "([mode] enable_spin_orbital = true, calc_scheme = \"general\") "
                "that computes chi0q itself or reads a chi0q_init marked "
                "spin_orbital = True.")
    if "spin_orbital" not in data:
        raise ValueError(
            "chi0q file '{}' has no spin_orbital marker, so it cannot be told "
            "from a non-spin-orbital chi0q of the same shape. {}".format(
                file_name, hint))
    val = np.asarray(data["spin_orbital"])
    if val.ndim != 0 or val.dtype != np.bool_:
        raise ValueError(
            "chi0q file '{}': malformed spin_orbital marker ({!r}); expected a "
            "boolean. {}".format(file_name, data["spin_orbital"], hint))
    if not bool(val[()]):
        raise ValueError(
            "chi0q file '{}' records spin_orbital = False: it was not written "
            "by a spin-orbital {} run and cannot be used with [mode] "
            "enable_spin_orbital = true. {}".format(file_name, producer, hint))


def _refuse_so_archive_in_non_so_mode(data, file_name,
                                      check_index_convention=False):
    """Refuse a spin-orbital susceptibility archive in a non-spin-orbital run.

    A spin-orbital archive over N physical orbitals has nd = 2N and so the
    shape of a non-spin-orbital general archive over 2N orbitals; read as the
    latter it would be solved with the wrong index interpretation and no
    error. Spin-orbital RPA/FLEX producers write the scalar bool
    ``spin_orbital = True`` (and a FLEX chiq the tag ``chi_convention =
    "spin_orbital"``); non-spin-orbital archives carry neither, so they load
    unchanged. A marker that is not a scalar bool is refused too (the
    spin-orbital side does the same). Checked on the open npz, before any
    array is read.

    ``check_index_convention`` (Green-function files only) also refuses
    ``index_convention == "spin_block"``: a spin-orbital FLEX green.npz can
    carry that tag alone, whereas the tag sits on EVERY chi0q/chiq archive
    (non-spin-orbital runs included), so the susceptibility readers leave
    it off."""
    hint = ("Enable [mode] enable_spin_orbital = true, or regenerate the "
            "archive with a non-spin-orbital run.")
    if "spin_orbital" in data:
        # a spin-orbital FLEX green.npz carries the marker too; it is named
        # by its payload member (every other archive keeps this wording)
        noun = "Green-function file" if "green" in data else "susceptibility file"
        val = np.asarray(data["spin_orbital"])
        if val.ndim != 0 or val.dtype != np.bool_:
            raise ValueError(
                "{} '{}': malformed spin_orbital marker "
                "({!r}); expected a boolean. {}".format(
                    noun, file_name, data["spin_orbital"], hint))
        if bool(val[()]):
            raise ValueError(
                "{} '{}' was written by a spin-orbital run "
                "(spin_orbital = True) and cannot be used without [mode] "
                "enable_spin_orbital = true: its shape would be misread as a "
                "non-spin-orbital file over twice as many orbitals. {}".format(
                    noun, file_name, hint))
    if ("chi_convention" in data
            and str(data["chi_convention"]) == _SO_FLEX_CHI_CONVENTION):
        raise ValueError(
            "FLEX file '{}' is tagged chi_convention = {!r} (a spin-orbital "
            "FLEX output) and cannot be used without [mode] "
            "enable_spin_orbital = true. {}".format(
                file_name, _SO_FLEX_CHI_CONVENTION, hint))
    from hwave.solver import spin_orbital as _spin_orbital
    if (check_index_convention and "index_convention" in data
            and str(data["index_convention"]) == _spin_orbital.INDEX_CONVENTION):
        raise ValueError(
            "Green-function file '{}' is tagged index_convention = {!r} (a "
            "spin-orbital run's output, spin_orbital = True) and cannot be "
            "used without [mode] enable_spin_orbital = true: its orbital "
            "axis is the spin-block spin-orbital index, not a physical "
            "orbital. {}".format(file_name, _spin_orbital.INDEX_CONVENTION,
                                 hint))


def _refuse_so_green_file(green_path, label="bond_green"):
    """Open a Green-function npz lazily and refuse a spin-orbital one.

    O(1): only the scalar marker keys are read, never the ``green`` array.
    A missing / unreadable file is left to :func:`_load_green_npz`, whose
    error is the more specific one."""
    if not os.path.isfile(green_path):
        return
    try:
        data = np.load(green_path)
    except (OSError, ValueError, EOFError, zipfile.BadZipFile,
            AttributeError):
        return
    try:
        _refuse_so_archive_in_non_so_mode(data, green_path,
                                          check_index_convention=True)
    finally:
        close = getattr(data, "close", None)
        if close is not None:
            close()


def _preflight_refuse_so_archives(input_dict, chi0q_mode):
    """Metadata-only form of :func:`_refuse_so_archive_in_non_so_mode`.

    Run by the public entries before the Hamiltonian is diagonalized or any
    Green function / interaction table is built, so a refusal costs nothing.
    Opens each archive that the non-spin-orbital run will read (the chi0q for
    ``chi0q_mode = "load"``, the FLEX chi_s/chi_c for ``"flex"``) and reads
    only the scalar markers. A missing file is left to the loader, whose
    error is the more specific one. The loader checks stay as defence in
    depth."""
    if chi0q_mode == "load":
        output_info = input_dict["file"]["output"]
        paths = [os.path.join(output_info["path_to_output"],
                              output_info.get("chi0q", "chi0q") + ".npz")]
    elif chi0q_mode == "flex":
        paths = list(_resolve_flex_paths(input_dict)[:2])
    else:
        return
    for path in paths:
        if not os.path.isfile(path):
            continue
        try:
            data = np.load(path)
        except (OSError, ValueError, EOFError, zipfile.BadZipFile,
                AttributeError):
            continue    # unreadable / not an npz: the loader reports it
        if not hasattr(data, "files"):
            continue    # a bare .npy under an .npz name: the loader reports it
        with data:
            _refuse_so_archive_in_non_so_mode(data, path)


def _load_chi0q(input_dict, norb=None, so_run=None):
    """Load chi0q from NPZ file produced by H-wave RPA solver.

    Parameters
    ----------
    input_dict : dict
        Parsed TOML input dictionary.
    so_run : _SORun, optional
        Spin-orbital mode (#88): the archive's provenance is
        checked against this run (:func:`_so_check_archive`, the load-mode
        nmat rule) after the marker checks and before the chi0q array is
        read, and the result is returned as a third item.

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array.
    static_index : int
        Position of the zero bosonic frequency along the frequency axis
        (determined from the freq_index metadata).
    so_info : SOArchiveInfo
        Only with ``so_run``.
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
    enable_spin_orbital = _resolve_spin_orbital_flag(input_dict)
    validate_chi0q_index_convention(data, enable_spin_orbital, file_name)
    so_info = None
    if enable_spin_orbital:
        _require_so_chi0q_marker(data, file_name)
        if so_run is not None:
            so_info = _so_check_archive(data, file_name, so_run, "load", "chi0q",
                                        _SO_RPA_PRODUCER)
    else:
        _refuse_so_archive_in_non_so_mode(data, file_name)
    chi0q = data["chi0q"]
    logger.info("chi0q shape: {}".format(chi0q.shape))

    # coeff_tail provenance check (issue #80): the Matsubara tail correction
    # changes chi0q at O(1) (lambda differed by 50% in the reproduction), so
    # a config whose effective coeff_tail differs from the file's would
    # recompute DIFFERENT physics under chi0q_mode="calc". Warn so load-vs-
    # calc comparisons are not silently inconsistent. Files without the key
    # (older builds) load silently as before.
    if "coeff_tail" in data:
        # type-strict provenance normalization, mirroring
        # rpa._validate_chi0q_provenance (round-8 review)
        import numbers
        _val = np.asarray(data["coeff_tail"])
        _item = _val[()] if _val.ndim == 0 else None
        if (_item is None or isinstance(_item, (bool, np.bool_))
                or not isinstance(_item, numbers.Real)
                or not np.isfinite(float(_item))):
            raise ValueError(
                "chi0q file '{}': malformed coeff_tail ({!r})".format(
                    file_name, data["coeff_tail"]))
        file_tail = float(_item)
        config_tail = float(input_dict.get("mode", {}).get("param", {})
                            .get("coeff_tail", 0.0))
        if file_tail != config_tail:
            logger.warning(
                "chi0q file '{}' was produced with coeff_tail = {} but this "
                "config's effective value is coeff_tail = {}; results are "
                "NOT comparable with a chi0q_mode=\"calc\" run under this "
                "config. Set [mode.param] coeff_tail = {} to match the "
                "file.".format(file_name, file_tail, config_tail, file_tail))
        # Endpoint-convention gate (issue #134): a nonzero-tail chi0q
        # produced before the branch-mean endpoint fix carries a pre-fix
        # O(1/Nmat) error indistinguishable from the array itself; refuse
        # it rather than feed it to the Eliashberg solver. Zero-tail files
        # and legacy files without the coeff_tail key are exempt (same
        # policy as the RPA chi0q_init reader).
        if file_tail != 0.0:
            from hwave.solver.rpa import TAIL_ENDPOINT_CONVENTION
            te = None
            if "tail_endpoint" in data:
                val = np.asarray(data["tail_endpoint"])
                te = str(val[()]) if val.ndim == 0 else repr(val)
            if te is None:
                raise ValueError(
                    "chi0q file '{}' records coeff_tail = {} but no "
                    "tail_endpoint marker: it was produced before the "
                    "equal-time endpoint fix (issue #134) and its "
                    "tail-corrected values carry the pre-fix O(1/Nmat) "
                    "endpoint error. Recompute the bubble with this "
                    "version.".format(file_name, file_tail))
            if te != TAIL_ENDPOINT_CONVENTION:
                raise ValueError(
                    "chi0q file '{}' carries unrecognized tail_endpoint = "
                    "{!r} (this build implements {!r}); refusing to use a "
                    "bubble whose endpoint treatment is unknown.".format(
                        file_name, te, TAIL_ENDPOINT_CONVENTION))

    # Unconditional marker validation (round-10 review): a PRESENT
    # marker must be checked even when the layout/grid cannot be
    # identified (absent CellShape, unknown shape).
    from hwave.solver.rpa import check_momentum_marker
    check_momentum_marker(data, file_name)
    # Fourier-sign provenance gate (issue #133): chi0q is q-labeled; a
    # pre-#133 file carries flipped labels. Legacy files are accepted only
    # when the payload is elementwise q-even (see the validator). The q
    # axes are chosen STRUCTURALLY from the accepted layouts (round-2
    # review: a size search probed the frequency axis when nfreq == nvol,
    # and skipped the already-expanded reference layouts entirely):
    #   4D/6D raw  (nfreq, nvol, ...)            -> axis 1
    #   5D/7D spin-diag (2, nfreq, nvol, ...)    -> axis 2
    #   6D ref  (no, no, Nx, Ny, Nz, nfreq)      -> axes (2, 3, 4)
    #   8D ref  (no, no, no, no, Nx, Ny, Nz, nfreq) -> axes (4, 5, 6)
    from hwave.solver.rpa import validate_momentum_convention
    _cs = input_dict.get("mode", {}).get("param", {}).get("CellShape")
    if _cs is None:
        # no lattice configured (unit-level callers): the grid, and with
        # it the momentum axes, cannot be identified -- production
        # configs always carry CellShape, so no production bypass
        _grid = None
    else:
        _cs = list(_cs)
        while len(_cs) < 3:
            _cs.append(1)
        _grid = tuple(int(x) for x in _cs)
    _nvol = int(np.prod(_grid)) if _grid is not None else -1
    # With the orbital count known (the production entry always passes
    # it), validate the EXACT supported layouts up front (round-6 review:
    # without norb a malformed 8D file passed both marked and unmarked;
    # with norb the (8,8,2,2,2,2)-style raw/ref ambiguity also vanishes,
    # since norb cannot be two values at once).
    _qax = None
    _resolved = False
    if norb is not None and _grid is not None:
        nv, no = _nvol, int(norb)
        # Full trailing shapes for every layout (round-7 review: partial
        # slices accepted malformed spin-diag arrays), and the layout
        # matched HERE directly selects the q axes -- re-running the
        # norb-blind structural routing below raised 'matches BOTH' for
        # shapes this gate had already disambiguated.
        if chi0q.ndim == 4 and tuple(chi0q.shape[1:]) == (nv, no, no):
            _qax = 1
        elif (chi0q.ndim == 6
                and tuple(chi0q.shape[1:]) == (nv, no, no, no, no)):
            _qax = 1
        elif (chi0q.ndim == 6
                and tuple(chi0q.shape[:5]) == (no, no) + _grid):
            _qax = (2, 3, 4)
        elif (chi0q.ndim == 8
                and tuple(chi0q.shape[:7]) == (no,) * 4 + _grid):
            _qax = (4, 5, 6)
        elif (chi0q.ndim == 5 and chi0q.shape[0] == 2
                and tuple(chi0q.shape[2:]) == (nv, no, no)):
            _qax = 2
        elif (chi0q.ndim == 7 and chi0q.shape[0] == 2
                and tuple(chi0q.shape[2:]) == (nv, no, no, no, no)):
            _qax = 2
        else:
            raise ValueError(
                "chi0q file '{}': shape {} matches no supported layout "
                "for norb = {} and CellShape {}; the file is malformed "
                "or from a different system. Regenerate it with the "
                "current version.".format(file_name, chi0q.shape, no,
                                          list(_grid)))
        _resolved = True
    if _resolved or _grid is None:
        pass
    elif chi0q.ndim in (5, 7) and chi0q.shape[0] == 2:
        _qax = 2
    elif chi0q.ndim == 8:
        _qax = (4, 5, 6)
    elif chi0q.ndim == 6:
        # raw (nfreq, nvol, norb x4) vs ref (norb, norb, Nx, Ny, Nz,
        # nfreq): decide by the FULL structure of both patterns
        # (round-3 review: testing shape[1] != nvol misrouted a ref file
        # with norb == nvol onto the orbital axis). Truly ambiguous
        # shapes fail closed unless the marker decides.
        _is_ref = (chi0q.shape[0] == chi0q.shape[1]
                   and tuple(chi0q.shape[2:5]) == _grid)
        _is_raw = (chi0q.shape[1] == _nvol
                   and len(set(chi0q.shape[2:6])) == 1)
        if _is_ref and _is_raw:
            if "momentum_convention" in getattr(data, "files", []):
                _qax = (2, 3, 4)   # any axis: the validator accepts on tag
            else:
                raise ValueError(
                    "chi0q file '{}': shape {} matches BOTH the raw and "
                    "the reference 6D layout and carries no "
                    "momentum_convention marker (issue #133) -- the "
                    "momentum axes cannot be identified safely. "
                    "Regenerate the file with the current version, which "
                    "stamps the marker.".format(file_name, chi0q.shape))
        elif _is_ref:
            _qax = (2, 3, 4)
        elif _is_raw:
            _qax = 1
        else:
            # neither complete pattern matches: fail closed HERE,
            # independent of provenance (round-5 review: routing through
            # the validator let a matching marker return early, and the
            # downstream converter then silently RESHAPED the unknown
            # layout, reinterpreting orbital axes as data -- a marker can
            # establish the Fourier sign, never the array layout)
            raise ValueError(
                "chi0q file '{}': shape {} matches neither the raw "
                "(nfreq, nvol, norb x4) nor the reference "
                "(norb, norb, Nx, Ny, Nz, nfreq) 6D layout for CellShape "
                "{}; the file is malformed or from a different lattice. "
                "Regenerate it with the current version.".format(
                    file_name, chi0q.shape, list(_grid)))
    elif chi0q.ndim == 4:
        _qax = 1
    if _qax is not None:
        validate_momentum_convention(data, file_name, chi0q, _qax, _grid)

    freq_index, file_nmat = _read_freq_meta(data)
    # Identify the frequency axis from the array LAYOUT, never from the
    # freq_index length (a restricted freq_index can coincidentally match
    # an orbital axis).  4D raw: axis 0.  8D ref: last axis.  5D/7D
    # spin-diagonal (leading 2): axis 1.  6D is either
    # raw (nmat, nvol, norb^4; the last four axes are equal) or ref
    # (norb, norb, Nx, Ny, Nz, nmat; the first two axes are equal) --
    # disambiguate structurally, with the freq_index length only as the
    # tiebreaker for degenerate shapes.  The length is computed even for
    # metadata-less files: _static_freq_position now needs it to decide
    # whether the stored axis IS the configured full grid (nfreq ==
    # config_nmat, the only case where centering is unambiguous).
    nfi = np.asarray(freq_index).size if freq_index is not None else None
    if _qax == 1:
        # raw layout (nfreq, nvol, norb...): frequency axis is axis 0
        nfreq = chi0q.shape[0]
    elif _qax == 2:
        # spin-diagonal block-leading layout (2, nfreq, nvol, norb...)
        nfreq = chi0q.shape[1]
    elif _qax in ((2, 3, 4), (4, 5, 6)):
        # reference layout (norb..., Nx, Ny, Nz, nfreq): axis is last
        nfreq = chi0q.shape[-1]
    elif chi0q.ndim in (5, 7) and chi0q.shape[0] == 2:
        # spin-diagonal block-leading layout, identifiable from the array
        # structure alone: the routing above leaves _qax unset when the run
        # carries no CellShape, and the last axis is an ORBITAL axis here,
        # so reading it as the frequency length would feed the full-grid
        # check the wrong number.
        nfreq = chi0q.shape[1]
    elif chi0q.ndim == 4:
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
        elif nfi is not None and nfi == chi0q.shape[0]:
            nfreq = chi0q.shape[0]
        else:
            nfreq = chi0q.shape[-1]
    else:
        # layout not resolved by any of the above: keep the historical
        # last-axis reading, but say so -- for an unsupported shape that
        # axis is a guess, and the static-slice guard downstream judges a
        # length this function could not identify.
        logger.debug(
            "chi0q file '{}': shape {} matches no known layout; taking the "
            "last axis as the frequency axis".format(file_name,
                                                     chi0q.shape))
        nfreq = chi0q.shape[-1]
    config_nmat = input_dict.get("mode", {}).get("param", {}).get("Nmat",
                                                                _DEFAULT_NMAT)
    static_index = _static_freq_position(freq_index, nfreq, config_nmat,
                                         file_name, file_nmat=file_nmat)
    if static_index is not None and static_index != nfreq // 2:
        logger.info("static (zero bosonic frequency) slice at index {} "
                    "of {}".format(static_index, nfreq))
    if so_run is not None:
        return chi0q, static_index, so_info
    return chi0q, static_index


def _auto_chi0q_tensor_scheme(files):
    """hwave_sc's ``chi0q_tensor = "auto"`` rule, behaviour-frozen at its
    1.0.x semantics via ``scheme.CAPABILITIES[...].sc_legacy_forcing``
    (True for CoulombInter/Hund/Ising/Exchange/PairHop/Coulomb; False for
    CoulombIntra/PairLift). Case-insensitive; an unknown key is IGNORED
    (fail-closed discovery is deliberately NOT adopted here -- recorded
    follow-up). Pinned by tests/test_sc_auto_scheme_freeze.py."""
    from hwave.solver.scheme import CAPABILITIES, canonical_name
    for k in (files or {}):
        canon = canonical_name(k)
        if canon is not None and CAPABILITIES[canon].sc_legacy_forcing:
            return "general"
    return "reduced"


def _calc_chi0q_internal(input_dict, chi0q_tensor="auto",
                         precomputed_mu=None, prepared=None,
                         static_only=False, so_context=None):
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
    prepared : _MeanFieldBand, optional
        A band already assembled by :func:`_mean_field_band` (issue #273).
        Its RPA solver -- input read once, ``_calc_epsilon_k`` already run --
        is reused, so chi0q is built on exactly the H0(k) the caller used
        for mu and the kernel Green function. ``chi0q_tensor`` is then
        ignored (the solver was constructed with the matching scheme).
    static_only : bool, optional
        If true, compute only the ``nmat // 2`` bosonic slice (the zero
        bosonic frequency the static Eliashberg solve reads) and return it as
        a one-frequency axis (``nmat`` -> 1 below; read it with static index
        0). It equals the full-grid slice to round-off; memory drops from
        ``O(nmat * nvol * nd**4)`` to ``O(nvol * nd**4)`` (issue #272).
        Applies to the ``prepared`` band path too. Dynamic callers must leave
        it false: they need every bosonic frequency.
    so_context : spin_orbital.SOEntryContext, optional
        The entry context of a spin-orbital hwave_sc run: the nested RPA
        reuses it (no second compile of the interaction, no second
        ``orbital_site`` WARNING).

    Returns
    -------
    chi0q : ndarray
        Bare susceptibility array. Shape depends on chi0q_tensor mode:
        - reduced: (nmat, nvol, norb, norb)
        - general: (nmat, nvol, norb, norb, norb, norb)
    """
    if prepared is not None:
        solver = prepared.solver
        logger.info("Computing chi0q internally on the mean-field one-body "
                    "band ({}; calc_scheme={})...".format(
                        ", ".join(prepared.sources), solver.calc_scheme))
        return _chi0q_from_prepared_solver(solver, precomputed_mu,
                                           static_only=static_only)

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
    calc_scheme = _chi0q_tensor_scheme(chi0q_tensor,
                                       info_inputfile.get("interaction", {}))

    info_mode_rpa = dict(info_mode)
    info_mode_rpa["calc_scheme"] = calc_scheme
    # hwave_sc has already reported an undeclared orbital_site (#234)
    info_mode_rpa["_orbital_site_reported"] = True
    # Forward the RESOLVED boolean, not the raw value: RPA mixes truthiness
    # and == True checks on this flag, so a string "false" would diverge
    # internally (round-1 review of #83's guard).
    info_mode_rpa["enable_spin_orbital"] = _resolve_spin_orbital_flag(
        input_dict)
    if so_context is not None:
        from hwave.solver import spin_orbital
        spin_orbital.context_for_nested(info_mode_rpa, so_context)

    logger.info("Computing chi0q internally (calc_scheme={})...".format(calc_scheme))

    # Read input files via QLMSkInput
    read_io = read_input_k.QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")

    # Create RPA solver (sets up Lattice, Interaction, parameters)
    solver = sol_rpa.RPA(ham_info, info_log, info_mode_rpa)
    if info_mode_rpa["enable_spin_orbital"]:
        # the spin-orbital Eliashberg vertex needs the full spinful chi0
        # (nd = 2 * norb spin-orbitals) even when the transfer happens to be
        # spin-diagonal, as the spin-orbital FLEX does (#83)
        solver._keep_spinful = True

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

    chi0q = solver._calc_chi0q(green0, green0_tail, beta,
                               static_only=static_only)
    del green0, green0_tail

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


def _chi0q_tensor_scheme(chi0q_tensor, interaction_files):
    """RPA ``calc_scheme`` for an ``[eliashberg] chi0q_tensor`` value:
    ``"auto"`` follows :func:`_auto_chi0q_tensor_scheme`, ``"general"`` is
    general, anything else is reduced (the historical rule)."""
    if chi0q_tensor == "auto":
        return _auto_chi0q_tensor_scheme(interaction_files)
    if chi0q_tensor == "general":
        return "general"
    return "reduced"


def _chi0q_from_prepared_solver(solver, precomputed_mu, static_only=False):
    """chi0q of an RPA solver whose ``_calc_epsilon_k`` has already run.

    Same arithmetic as the tail of :func:`_calc_chi0q_internal`; only the
    spin-free layout reaches here (:func:`_mean_field_band` refuses the
    others), so the block axis is always length 1.
    """
    beta = 1.0 / solver.T
    if precomputed_mu is not None:
        mu = precomputed_mu
        logger.info("Using precomputed mu = {}".format(mu))
    elif solver.calc_mu:
        dist, mu = solver._find_mu(solver.Ncond / 2, solver.T)
    else:
        mu = solver.mu_value
    green0, green0_tail = solver._calc_green(beta, mu)
    chi0q = solver._calc_chi0q(green0, green0_tail, beta,
                               static_only=static_only)
    del green0, green0_tail
    assert chi0q.shape[0] == 1
    chi0q = chi0q[0]
    logger.info("chi0q computed internally, shape: {}".format(chi0q.shape))
    return chi0q


def _read_spin_free_tol(eli_param):
    """``[eliashberg] spin_free_tol`` (case-insensitive): finite float >= 0,
    default :data:`SPIN_FREE_ATOL`."""
    import numbers
    from requests.structures import CaseInsensitiveDict
    raw = CaseInsensitiveDict(eli_param or {}).get("spin_free_tol")
    if raw is None:
        return SPIN_FREE_ATOL
    # type-strict (the _resolve_bond_coeff_tail pattern): a boolean would
    # become 1.0 and a quoted number would be coerced -- both silently
    # change the tolerance, so both are refused
    if isinstance(raw, (bool, np.bool_)) or not isinstance(raw, numbers.Real):
        raise ValueError(
            "[eliashberg] spin_free_tol must be a real number, got "
            "{!r}.".format(raw))
    tol = float(raw)
    if not np.isfinite(tol) or tol < 0.0:
        raise ValueError(
            "[eliashberg] spin_free_tol must be a finite number >= 0; got "
            "{!r}.".format(raw))
    return float(tol)


def _load_band_into_solver(solver, up):
    """Make ``solver`` diagonalize exactly ``up`` (the checked band) for
    chi0q: it is handed to RPA's own ``_calc_epsilon_k`` as an explicitly
    spin-free trans_mod (both spin blocks equal to ``up``, no mixing), with
    the Extern term switched off for that call because it is already inside
    ``up``. RPA's H0_k must then equal ``up`` exactly."""
    nvol, norb, _ = up.shape
    tm = np.zeros((nvol, 2 * norb, 2 * norb), dtype=np.complex128)
    tm[:, :norb, :norb] = up
    tm[:, norb:, norb:] = up
    ham = solver.ham_info
    saved_extern = ham.ham_extern_q
    ham.ham_extern_q = None
    try:
        solver._calc_epsilon_k({"trans_mod": tm})
    finally:
        ham.ham_extern_q = saved_extern
    H0 = np.asarray(solver.H0_k)
    if (solver.spin_mode != "spin-free" or H0.shape != (1,) + up.shape
            or not np.array_equal(H0[0], up)):
        raise RuntimeError(
            "hwave_sc: the RPA solver did not reproduce the checked "
            "one-body band (spin_mode={}, shape {}); refusing to compute "
            "chi0q on a different band.".format(solver.spin_mode, H0.shape))


def _mean_field_band(input_dict, sources, calc_scheme, spatial_shape, norb,
                     spin_free_tol=SPIN_FREE_ATOL, allow_average=True,
                     allow_hermitize=True):
    """Assemble the one-body band through the RPA machinery (issue #273).

    The band is built exactly as the internal chi0q builds it --
    ``QLMSkInput`` + ``RPA.read_init`` + ``RPA._calc_epsilon_k`` -- so that
    trans_mod (T + Sigma_HF written by UHFk), T + Sigma_HF[green_init] (UHFk's
    Hartree + Fock mean field of G, issue #278) and the Extern field define
    the SAME H0(k) for chi0q, the chemical potential
    and the kernel Green function.

    Only a spin-independent band is accepted: the singlet/triplet kernel
    carries no spin index, so a spin-diagonal or spin-mixing mean field
    cannot be represented. The band is assembled ONCE by RPA
    (:func:`_checked_band`: finiteness, then the absolute tolerance
    :data:`SPIN_FREE_ATOL`), after the input files are read by the RPA
    readers and before H0 is diagonalized -- so before any chi0q, pair
    bubble or kernel is built, and before a chi0q file is loaded. The
    checked band (its up-up block, Extern included) is the band used for
    epsilon(k), mu, the kernel Green function and -- through
    :func:`_load_band_into_solver`, which verifies RPA reproduces it
    exactly -- chi0q.
    """
    from requests.structures import CaseInsensitiveDict
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    fin = (input_dict.get("file", {}) or {}).get("input", {}) or {}
    for key in _ONE_BODY_FILE_KEYS:
        spellings = [k for k in fin if str(k).lower() == key]
        if len(spellings) > 1:
            raise ValueError(
                "[file.input] {} is given more than once with different "
                "spellings ({}); keys are case-insensitive.".format(
                    key, ", ".join(repr(k) for k in spellings)))
    ci = CaseInsensitiveDict(fin)
    # Only the keys the band needs: chi0q_init / UHF initial files are not
    # read here (an RPA read_init would otherwise load a chi0q_init array).
    # The lower-case spelling is what RPA.read_init looks up.
    info_inputfile = {"path_to_input": ci.get("path_to_input", ""),
                      "interaction": ci.get("interaction", {})}
    for key in _ONE_BODY_FILE_KEYS:
        if key in ci:
            info_inputfile[key] = ci[key]

    info_log = dict(input_dict.get("log", {}) or {})
    info_log.setdefault("print_level", 1)
    info_log.setdefault("print_step", 1)
    import copy
    # deep copy: the RPA constructor normalizes [mode.param] in place, and
    # the caller's dict has to stay as the user wrote it
    info_mode_rpa = copy.deepcopy(dict(input_dict.get("mode", {}) or {}))
    info_mode_rpa["calc_scheme"] = calc_scheme
    info_mode_rpa["_orbital_site_reported"] = True
    info_mode_rpa["enable_spin_orbital"] = _resolve_spin_orbital_flag(
        input_dict)

    logger.info("One-body band from {} (RPA H0 assembly; used for chi0q, "
                "mu and the kernel Green function)".format(", ".join(sources)))
    read_io = read_input_k.QLMSkInput(info_inputfile)
    solver = sol_rpa.RPA(read_io.get_param("ham"), info_log, info_mode_rpa)
    green_info = read_io.get_param("green")
    green_info.update(solver.read_init(info_inputfile))
    if "green_init" in green_info and "trans_mod" not in green_info:
        # G-level refusal, before the band is assembled
        _check_green_init_spin_free(green_info["green_init"], solver.norb,
                                    spin_free_tol)

    # Assemble the band ONCE, check it (spin independence, finiteness),
    # before any diagonalization and before any chi0q / pair-bubble /
    # kernel work; the checked band is the band that is used.
    up = _checked_band(solver, green_info, sources, spin_free_tol,
                       allow_average=allow_average,
                       allow_hermitize=allow_hermitize)

    Nx, Ny, Nz = spatial_shape
    if up.shape != (Nx * Ny * Nz, norb, norb):
        raise ValueError(
            "hwave_sc: the one-body band assembled from {} has shape {}, "
            "expected {} for the {}x{}x{} grid with norb={}.".format(
                " / ".join(sources), up.shape,
                (Nx * Ny * Nz, norb, norb), Nx, Ny, Nz, norb))
    # chi0q is computed by this solver on exactly the checked band
    _load_band_into_solver(solver, up)
    H0 = up[np.newaxis]
    # RPA (1, nvol, a, b) with nvol in C order over (Nx, Ny, Nz) -- the same
    # flattening as _green_sc_to_canonical -- to sc.py's (a, b, Nx, Ny, Nz).
    epsilon_k = np.ascontiguousarray(
        H0[0].reshape(Nx, Ny, Nz, norb, norb).transpose(3, 4, 0, 1, 2),
        dtype=np.complex128)
    bare_q = getattr(solver.ham_info, "ham_trans_q", None)
    same_as_bare = (bare_q is not None
                    and np.shape(bare_q) == H0[0].shape
                    and np.array_equal(H0[0], bare_q))
    return _MeanFieldBand(epsilon_k=epsilon_k, solver=solver,
                          sources=tuple(sources), same_as_bare=same_as_bare)


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
    from hwave.solver.declarations import (validate_hermitian_closure,
                                           validate_no_onsite_same_orbital_inter)

    for itype in interaction_types:
        if itype in files:
            f = os.path.join(path_to_input, files[itype])
            logger.info("Reading {} from {}".format(itype, f))
            tbl = wan90.read_w90(f)
            # issue #93: fail fast on non-Hermitian-closed declarations,
            # with the same rule and tolerance as the k-space reader
            validate_hermitian_closure(itype, tbl, source=f)
            validate_no_onsite_same_orbital_inter(itype, tbl, source=f)
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


def _read_so_extern(input_dict, norb_phys):
    """The spin-orbital ``Extern`` table, read and checked (issue #88).

    In spin-orbital mode the Extern file holds ``H_ab(R)`` in PHYSICAL
    orbital indices ``1..norb_phys``; the operator is ``coeff_extern *
    sum H_ab(R) sigma^z_{ss'} c+_{as} c_{bs'}``. The header (line 2) must
    declare ``norb_phys`` and every index must be in range
    (:func:`hwave.solver.spin_orbital.check_extern_indices`); a malformed
    header is refused by :func:`hwave.qlmsio.wan90.read_w90_header_dim`.
    Run in every ``chi0q_mode`` (calc, load, flex), before any heavy work.

    Returns
    -------
    (extern_table, file_name) : tuple
        ``(None, None)`` when ``[file.input.interaction]`` names no Extern
        file (the key is matched case-insensitively).
    """
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    files = CaseInsensitiveDict(input_dict["file"]["input"]["interaction"])
    if "Extern" not in files:
        return None, None
    file_name = os.path.join(files.get("path_to_input", ""), files["Extern"])
    dim = wan90.read_w90_header_dim(file_name)
    table = wan90.read_w90(file_name)
    spin_orbital.check_extern_indices(table, norb_phys, file_name,
                                      declared_norb=dim)
    logger.info("Reading Extern from {} ({} physical orbitals)".format(
        file_name, norb_phys))
    return table, file_name


def _so_check_archive(data, path, so_run, nmat_rule, kind, producer, **kw):
    """:func:`hwave.solver.fingerprint.check_so_archive` of the ``kind``
    archive (``"chi0q"``, ``"chiq"``, ``"green"`` or ``"sigma"``) against
    ``so_run``: its band, state and the interaction fingerprints of its
    entry context; this run consumes on-site (``"so_onsite"``) artifacts.
    ``producer`` (:data:`_SO_RPA_PRODUCER` or :data:`_SO_FLEX_PRODUCER`)
    is the run that writes the archive."""
    from hwave.solver import fingerprint
    return fingerprint.check_so_archive(
        data, path, kind=kind, run_band=so_run.band, run_state=so_run.state,
        nmat_rule=nmat_rule, run_interaction=so_run.context.interaction,
        need="so_onsite", expected_producer=producer, **kw)


def _so_green_fp(green_kw):
    """``fp_green`` of a Green function ``(nd, nd, Nx, Ny, Nz, nmat)`` as a
    hwave_sc gap archive records it."""
    from hwave.solver import fingerprint
    return fingerprint.green_digest(green_kw, "sc_o_o_xyz_f")


def _resolve_flex_dir(input_dict):
    """The directory the FLEX outputs are read from.

    ``[file.input] path_to_flex_output``, falling back to ``[file.output]
    path_to_output`` and finally to ``"output"``. Factored out so every
    consumer of the FLEX outputs (the chi/green loaders here and the bond
    archive of the dynamic bond path) resolves the SAME directory.
    """
    file_input = input_dict.get("file", {}).get("input", {})
    return file_input.get("path_to_flex_output",
                          input_dict.get("file", {}).get("output", {}).get(
                              "path_to_output", "output"))


def _resolve_flex_paths(input_dict):
    """Resolve the FLEX output directory and chi_s/chi_c/green file paths.

    Shared by `_load_flex_susceptibilities_full` (which reads the files) and
    `_load_flex_susceptibilities` (which re-derives the static frequency
    index from the same files' metadata), so the two never disagree about
    which files are in play.
    """
    flex_dir = _resolve_flex_dir(input_dict)

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
        green_w = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                   interactions=interactions)
    else:
        (chi_s_raw, chi_c_raw, chi_convention, ir_meta,
         legacy_tags) = _read_flex_chi_raw(
            input_dict, allow_ir=True, interactions=interactions)
        green_w, green_meta = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                               allow_ir=True,
                                               interactions=interactions)
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
    _non_so = not _resolve_spin_orbital_flag(input_dict)
    if _non_so:
        _refuse_so_archive_in_non_so_mode(data_s, chi_s_path)
    if not allow_ir:
        _reject_ir_native(data_s, chi_s_path, _STATIC_IR_HINT)
    chi_s_raw = data_s["chiq_s"] if "chiq_s" in data_s else data_s["chiq"]
    from hwave.solver.rpa import (validate_momentum_convention,
                                  check_momentum_marker)
    # unconditional (round-10): a present marker is validated even when
    # the layout-dependent evenness gate below cannot run
    check_momentum_marker(data_s, chi_s_path)
    # Fourier-sign provenance (issue #133): the H-wave chi layout --
    # uniform AND IR-native alike -- has the flattened q volume on axis 1
    # (round-4 review: IR-native files were wrongly exempted; no IR
    # consumer gates them).
    _cs = list(input_dict.get("mode", {}).get("param", {}).get(
        "CellShape", [1, 1, 1]))
    while len(_cs) < 3:
        _cs.append(1)
    _grid = tuple(int(x) for x in _cs)
    if chi_s_raw.ndim >= 2 and chi_s_raw.shape[1] == int(np.prod(_grid)):
        validate_momentum_convention(data_s, chi_s_path, chi_s_raw, 1,
                                     _grid)
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
    if _non_so:
        _refuse_so_archive_in_non_so_mode(data_c, chi_c_path)
    if not allow_ir:
        _reject_ir_native(data_c, chi_c_path, _STATIC_IR_HINT)
    chi_c_raw = data_c["chiq_c"] if "chiq_c" in data_c else data_c["chiq"]
    check_momentum_marker(data_c, chi_c_path)
    if chi_c_raw.ndim >= 2 and chi_c_raw.shape[1] == int(np.prod(_grid)):
        validate_momentum_convention(data_c, chi_c_path, chi_c_raw, 1,
                                     _grid)
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
        # #232: version 3 puts an off-site row (R, a, b) on the bond of a
        # Transfer row (R, a, b); every earlier file was computed with it on
        # the opposite bond. That only matters for a table that is not
        # symmetric under R -> -R.
        # PairLift is exempt for the same reason as in the orientation rule
        # above: its particle-hole content is exactly zero on both sides.
        placement = placement_sensitive_types(
            interactions,
            grid=(_grid if "CellShape" in input_dict.get("mode", {}).get(
                "param", {}) else None))
        if affected or orientation or resym or placement:
            versions = {}
            for data, path in ((data_s, chi_s_path), (data_c, chi_c_path)):
                if "sc_vertex_version" not in data:
                    reasons = []
                    if placement:
                        reasons.append(
                            "{} declare(s) off-site rows that are not "
                            "symmetric under R -> -R, which a file this old "
                            "placed on the bond opposite their Transfer rows "
                            "(#232)".format(", ".join(placement)))
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
                shown = "<unreadable>"
                try:
                    raw = data["sc_vertex_version"]
                    shown = repr(raw)
                    arr = np.asarray(raw)
                    if arr.dtype.kind not in "iufc":    # no bool, str or object stamps
                        raise ValueError("not numeric")
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
                        "sc_vertex_version field ({}).".format(path, shown))
            if len(set(versions.values())) != 1:
                raise ValueError(
                    "FLEX chi_s and chi_c carry different sc_vertex_version "
                    "values ({}); they must come from the same run.".format(
                        versions))
            ver = next(iter(versions.values()))
            if ver not in (2, 3):
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = {} "
                    "but this code supports version 2 or version 3 (the "
                    "#113 exact-diagonalization vertex content; version 3 "
                    "also places off-site rows on the bond of their Transfer "
                    "rows, #232); "
                    "regenerate the susceptibilities with the current "
                    "code.".format(ver))
            if ver == 2 and placement:
                raise ValueError(
                    "FLEX susceptibility files carry sc_vertex_version = 2, "
                    "and {} declare(s) off-site rows that are not symmetric "
                    "under R -> -R. Version-2 runs placed such a row on the "
                    "bond opposite its Transfer row (#232), so pairing those "
                    "susceptibilities with the current vertices would "
                    "silently mix two different interactions -- regenerate "
                    "the susceptibilities with the current code.".format(
                        ", ".join(placement)))
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


def _expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention):
    """Reshape H-wave chi ``(nfreq, nvol, nd, nd)`` to
    ``(nfreq, Nx, Ny, Nz, nd, nd)`` in the ``nd = norb^2`` Eliashberg space,
    resolving the spin-orbital-vs-orbital-pair layout from ``convention``.

    The two FLEX conventions have DIFFERENT physical layouts that shape alone
    cannot always tell apart:

    - ``"myo"`` (general full-vertex FLEX) is already in orbital-pair space
      ``nd_chi = norb^2``; passed through unchanged.
    - ``"kuroki"`` (reduced FLEX) is in spin-orbital reduced space
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
    # The reduced scheme keeps only the density-density diagonal of the
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


def _load_green_npz(green_path, norb, Nx, Ny, Nz, label="Green"):
    """Load a Green-function npz written in the H-wave layout, path-explicit.

    Shared by ``_load_flex_green`` (which resolves the path from the FLEX
    output directory) and by the ``[eliashberg] bond_green`` key (an
    externally supplied Green function, spec Goal). The file must hold a
    ``green`` array of shape ``(nblock, nfreq, nvol, norb, norb)`` -- what
    ``flex.save_results`` writes -- with ``nvol = Nx*Ny*Nz`` and the same
    ``norb`` as the model; anything else raises rather than being reshaped
    into silent nonsense.

    Returns ``(green, nfreq)`` with ``green`` in the sc.py layout
    ``(norb, norb, Nx, Ny, Nz, nfreq)``.
    """
    if not os.path.exists(green_path):
        raise FileNotFoundError(
            "{} file not found: {}".format(label, green_path))
    logger.info("Loading {} from: {}".format(label, green_path))
    data_g = np.load(green_path)
    _reject_ir_native(data_g, green_path, _STATIC_IR_HINT)
    if "green" not in data_g.files:
        raise ValueError(
            "{} file {} has no 'green' array (keys: {}); the expected format "
            "is the H-wave green.npz written by the FLEX/UHF solvers."
            .format(label, green_path, list(data_g.files)))
    green_raw = data_g["green"]
    if green_raw.ndim != 5:
        raise ValueError(
            "{} file {}: 'green' must have shape (nblock, nfreq, nvol, norb, "
            "norb), got {}".format(label, green_path, green_raw.shape))
    nblock_g, nfreq, nvol_g, norb1, norb2 = green_raw.shape
    if nblock_g != 1:
        raise ValueError(
            "{} file {}: 'green' carries {} spin blocks, but this loader "
            "unconditionally consumes block 0 -- the bond/FLEX Green "
            "consumers here are single-block. A file with nblock=2 (e.g. a "
            "legitimate spin-diagonal H-wave Green with separate G_up/G_down "
            "blocks) would otherwise be silently halved to its up block "
            "while being recorded as a trusted full Green. Re-run the "
            "producing calculation paramagnetically (nblock=1), or extract "
            "and save the single block you intend to use."
            .format(label, green_path, nblock_g))
    if (nvol_g != Nx * Ny * Nz) or (norb1 != norb) or (norb2 != norb):
        raise ValueError(
            "{} file {} does not match the model: it holds nvol={}, norb={}x{} "
            "while this run has nvol={} ({}x{}x{}) and norb={}."
            .format(label, green_path, nvol_g, norb1, norb2,
                    Nx * Ny * Nz, Nx, Ny, Nz, norb))
    green = green_raw[0].reshape(
        nfreq, Nx, Ny, Nz, norb, norb).transpose(4, 5, 1, 2, 3, 0).copy()
    return green, nfreq


class _UnsupportedNpyHeaderVersion(Exception):
    """Internal to :func:`_peek_green_npz_nfreq`: raised when a ``.npy``
    header version is one neither numpy's public ``read_array_header_X_Y``
    wrappers nor its internal version-generic dispatcher can parse.

    This is deliberately NOT a plain ``ValueError`` catch target inside
    :func:`_peek_green_npz_nfreq`'s own except clause: a genuinely unsupported
    version must surface as a clear error (review fix), not be swallowed into
    a silent ``None`` that reopens the double-preflight path the header peek
    exists to close.
    """


def _read_npy_header_shape(fh):
    """Read only the ``shape`` from a ``.npy`` header, for any format version
    ``np.load`` itself accepts (1.0, 2.0, 3.0, and any future version numpy
    grows a dedicated ``read_array_header_X_Y`` wrapper for).

    Dispatches on the version tuple from ``np.lib.format.read_magic`` to the
    matching public ``read_array_header_X_Y`` function when one exists (1.0,
    2.0). Format 3.0 (and any future version numpy's public API has not grown
    a dedicated wrapper for) has no such function, so this falls back to
    ``numpy.lib.format``'s own internal version-generic dispatcher -- exactly
    what ``np.load`` uses under the hood to parse ANY version it accepts, so
    the fallback keeps every such version working here too instead of
    silently giving up.

    Raises :class:`_UnsupportedNpyHeaderVersion` if the version is genuinely
    unknown to :mod:`hwave.solver.npy_header`.
    """
    try:
        return _npy_header.read_npy_header_shape(fh)
    except _npy_header.UnsupportedNpyHeaderVersion as exc:
        raise _UnsupportedNpyHeaderVersion(str(exc))


def _peek_green_npz_nfreq(green_path, label="Green"):
    """Return the ``nfreq`` axis of a Green npz WITHOUT materialising it.

    The bond path's resource preflight must run with the frequency count the
    run will really use -- and for an external ``bond_green`` that is the
    FILE's, not the configured ``Nmat`` (review fix: the preflight used to run
    with the configured value first, so an over-large ``Nmat`` could refuse a
    run that fits). Reading the npz member's array header alone keeps the
    preflight strictly ahead of the (possibly large) allocation.

    Supports every ``.npy`` header version ``np.load`` accepts (see
    :func:`_read_npy_header_shape`); a genuinely unsupported version raises
    rather than silently returning ``None``, since falling back to the
    configured Nmat for the first preflight is exactly the double-preflight
    behaviour this function exists to eliminate.

    Returns ``None`` when the frequency count cannot be determined from the
    header alone for an ordinary reason (missing file, no ``green`` member,
    unexpected rank); the caller then falls back to :func:`_load_green_npz`,
    which raises the informative error.
    """
    try:
        if not zipfile.is_zipfile(green_path):
            return None
        with zipfile.ZipFile(green_path) as zf:
            member = None
            for name in zf.namelist():
                if name in ("green.npy", "green"):
                    member = name
                    break
            if member is None:
                return None
            with zf.open(member) as fh:
                shape = _read_npy_header_shape(fh)
    except _UnsupportedNpyHeaderVersion as exc:
        raise ValueError(
            "{} file {}: {}".format(label, green_path, exc)) from exc
    except (OSError, ValueError):
        # unreadable/corrupt: let the real loader produce the diagnostic
        return None
    if len(shape) != 5:
        return None
    logger.info("%s %s carries %d Matsubara frequencies (read from the npz "
                "header)", label, green_path, int(shape[1]))
    return int(shape[1])


def _load_flex_green(input_dict, norb, Nx, Ny, Nz, allow_ir=False,
                     interactions=None):
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
    if interactions is not None:
        _grid_cs = input_dict.get("mode", {}).get("param", {}).get("CellShape")
        check_placement_stamp(
            data_g, green_path, interactions,
            grid=(tuple(int(x) for x in list(_grid_cs) + [1] * (3 - len(_grid_cs)))
                  if _grid_cs is not None else None),
            what="FLEX Green file")
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
        # An odd uniform axis cannot be a centered fermionic grid (it
        # contains wn = 0); reject at the file boundary, before any vertex
        # or pair-bubble construction (round-6 review).
        if nmat_g % 2 != 0:
            raise ValueError(
                "green file '{}': the frequency axis has {} points, which "
                "cannot be a centered fermionic Matsubara grid (even count "
                "required). Re-run the producing FLEX calculation.".format(
                    green_path, nmat_g))
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
        run_beta = None if run_T is None else _coerce_run_beta(run_T)
        if "beta" in data_g:
            # One REAL-NUMERIC finite positive scalar, validated on the
            # NARROWED float64 (see _finite_positive_float64): NaN/Inf make
            # every tolerance comparison False (silently "matching"), an
            # empty array would die on an incidental IndexError, a vector
            # would silently use its first element, float() would silently
            # discard a complex value's imaginary part, booleans would read
            # as 0/1, a string would raise an incidental TypeError, and a
            # wider-than-binary64 longdouble could pass elementwise checks
            # yet narrow to inf.
            file_beta = _finite_positive_float64(data_g["beta"])
            if file_beta is None:
                raise ValueError(
                    "green file '{}': beta metadata must be a single "
                    "finite positive real scalar, got {!r}. Regenerate "
                    "the file.".format(green_path,
                                       np.asarray(data_g["beta"])))
            # Explicitly symmetric relative tolerance (max of both
            # magnitudes), no absolute floor: an absolute term would wave
            # through large relative mismatches at small beta (high
            # temperature), where the grid differs the most.
            if run_beta is not None and abs(file_beta - run_beta) \
                    > 1.0e-8 * max(abs(file_beta), abs(run_beta)):
                raise ValueError(
                    "green file '{}' was produced at beta = {} but this "
                    "run uses beta = {}; the Matsubara grid (and the tail "
                    "correction built on it) would be inconsistent. Use "
                    "the producing run's temperature or regenerate the "
                    "file.".format(green_path, file_beta, run_beta))
        elif run_beta is not None:
            logger.warning(
                "green file '%s' carries no beta metadata (written before "
                "this field existed); the run's beta = %g is assumed to "
                "match the producing FLEX run. Verify the temperatures "
                "agree -- a mismatch silently corrupts the pair bubble.",
                green_path, run_beta)
    from hwave.solver.rpa import validate_momentum_convention
    # green layout (nblock, nfreq_or_nodes, nvol, norb, norb): q axis 2,
    # for uniform AND IR-native files alike (round-4 review)
    validate_momentum_convention(data_g, green_path, green_raw, 2,
                                 (Nx, Ny, Nz))
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

    green_dressed = _load_flex_green(input_dict, norb, Nx, Ny, Nz,
                                     interactions=interactions)
    return chis, chic, green_dressed, chi_convention


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


def _so_flex_static_slice(data, arr, path, nd, nvol, config_nmat, what):
    """The zero-bosonic-frequency slice ``(nvol, nd^2, nd^2)`` of a raw
    spinful FLEX array ``arr`` ``(nmat, nvol, nd, nd, nd, nd)`` of the open
    npz ``data`` read from ``path``, located through the file's
    ``freq_index`` / ``nmat`` metadata."""
    if arr.ndim != 6 or arr.shape[1:] != (nvol,) + (nd,) * 4:
        raise ValueError(
            "FLEX {} in '{}' has shape {}, not the general spinful layout "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}) over nd = {nd} spin-orbitals of "
            "this run. Re-run FLEX on the same model and grid.".format(
                what, path, arr.shape, nvol, nd=nd))
    fi, fn = _read_freq_meta(data)
    pos = _static_freq_position(fi, arr.shape[0], config_nmat, path,
                                file_nmat=fn)
    si = arr.shape[0] // 2 if pos is None else pos
    out = np.ascontiguousarray(arr[si].reshape(nvol, nd * nd, nd * nd),
                               dtype=complex)
    if not np.all(np.isfinite(out)):
        raise ValueError(
            "FLEX {} in '{}' contains non-finite values at the zero bosonic "
            "frequency; re-run the producing FLEX calculation.".format(what, path))
    return out


def _so_flex_full_axis(data, arr, path, nd, nvol, config_nmat, what):
    """The whole bosonic axis ``(nmat, nvol, nd^2, nd^2)`` of a raw spinful
    FLEX array ``arr`` ``(nmat, nvol, nd, nd, nd, nd)`` for the dynamic
    solver (issue #284): the stored axis must be the full uniform grid of the
    configured ``Nmat`` (``freq_index``, when recorded, ``0 .. Nmat - 1``;
    ``nmat``, when recorded, ``Nmat``), with bosonic index ``m`` at
    ``nu_m = (2 m - Nmat) pi T``."""
    if arr.ndim != 6 or arr.shape[1:] != (nvol,) + (nd,) * 4:
        raise ValueError(
            "FLEX {} in '{}' has shape {}, not the general spinful layout "
            "(nmat, {}, {nd}, {nd}, {nd}, {nd}) over nd = {nd} spin-orbitals of "
            "this run. Re-run FLEX on the same model and grid.".format(
                what, path, arr.shape, nvol, nd=nd))
    fi, fn = _read_freq_meta(data)
    nmat = int(config_nmat)
    full = (arr.shape[0] == nmat and (fn is None or int(fn) == nmat)
            and (fi is None or np.array_equal(np.asarray(fi).ravel(),
                                              np.arange(nmat))))
    if not full:
        raise ValueError(
            "FLEX {} in '{}' does not hold the full bosonic grid of [mode.param] "
            "Nmat = {} (stored axis {}, recorded nmat {}): frequency = "
            "\"dynamic\" needs the susceptibility at every bosonic frequency of "
            "the run's grid. Re-run FLEX with the same Nmat.".format(
                what, path, nmat, arr.shape[0], fn))
    out = np.ascontiguousarray(arr.reshape(nmat, nvol, nd * nd, nd * nd),
                               dtype=complex)
    for m in range(nmat):
        if not np.all(np.isfinite(out[m])):
            raise ValueError(
                "FLEX {} in '{}' contains non-finite values at bosonic index {}; "
                "re-run the producing FLEX calculation.".format(what, path, m))
    return out


class _SOFlexBorrowed(NamedTuple):
    """What a spin-orbital FLEX consumer borrowed (#88): ``info``, the
    :class:`hwave.solver.fingerprint.SOArchiveInfo` of the chiq archive;
    ``fp_green``, the digest of the Green function as consumed (computed once,
    in the verification against the archive); ``sigma_static``, the FLEX
    run's static self-energy ``(nvol, nd, nd)`` (for the closure density of
    that Green function)."""
    info: object
    fp_green: object
    sigma_static: object


def _so_flex_sigma_path(eli, flex_dir):
    """The sigma archive of a spin-orbital FLEX run: ``[eliashberg]
    flex_sigma`` (default ``"sigma.npz"``) in ``flex_dir``."""
    return os.path.join(flex_dir, eli.get("flex_sigma", "sigma.npz"))


def _load_flex_so_sigma(eli, flex_dir, nd, nvol, so_run, fp_green, green_path,
                        hint):
    """The static self-energy ``(nvol, nd, nd)`` of the spin-orbital FLEX run
    whose Green function digest is ``fp_green``, from ``[eliashberg]
    flex_sigma`` (default ``"sigma.npz"``): the archive is validated
    (:func:`_so_check_archive`) and must record the same ``fp_green``."""
    from hwave.solver import spin_orbital
    sigma_path = _so_flex_sigma_path(eli, flex_dir)
    if not os.path.exists(sigma_path):
        raise ValueError(
            "FLEX self-energy archive '{}' not found. {}, and also its sigma: "
            "hwave_sc records the closure density of the borrowed Green "
            "function, which needs the static self-energy. Set [file.output] "
            "sigma in the FLEX run and [eliashberg] flex_sigma to that "
            "file.".format(sigma_path, hint))
    with np.load(sigma_path) as ds:
        sinfo = _so_check_archive(ds, sigma_path, so_run, "exact", "sigma",
                                  _SO_FLEX_PRODUCER)
        if sinfo.producer_fp_green != fp_green:
            raise ValueError(
                "{} and {} come from different runs: '{}' records fp_green {} "
                "but the Green function '{}' has digest {}. Use the sigma and "
                "green archives of one spin-orbital FLEX run.".format(
                    os.path.basename(sigma_path), os.path.basename(green_path),
                    sigma_path, (sinfo.producer_fp_green or "none")[:12],
                    green_path, fp_green[:12]))
        conv = str(ds["sigma_convention"]) if "sigma_convention" in ds else None
        tag = str(ds["index_convention"]) if "index_convention" in ds else None
        if conv != "split" or tag != spin_orbital.INDEX_CONVENTION:
            raise ValueError(
                "FLEX self-energy archive '{}' records sigma_convention = {!r} "
                "and index_convention = {!r}; expected 'split' and {!r} (a "
                "spin-orbital Hartree-Fock FLEX run writes them).".format(
                    sigma_path, conv, tag, spin_orbital.INDEX_CONVENTION))
        st = np.asarray(ds["sigma_static"]) if "sigma_static" in ds else None
        want = (1, 1, nvol, nd, nd)
        # complex128 in either byte order (the archive check accepts both);
        # returned in native order
        if (st is None or st.dtype.kind != "c" or st.dtype.itemsize != 16
                or st.shape != want):
            raise ValueError(
                "FLEX self-energy archive '{}': sigma_static must be complex128 "
                "of shape {}, got {}. Regenerate it with a spin-orbital "
                "Hartree-Fock FLEX run.".format(
                    sigma_path, want,
                    None if st is None else "{} {}".format(st.dtype, st.shape)))
    return np.ascontiguousarray(st[0, 0], dtype=np.complex128)


def _load_flex_so(input_dict, nd, Nx, Ny, Nz, interactions=None, full=False,
                  so_run=None):
    """The outputs of a spin-orbital FLEX run (``calc_scheme = "general"``,
    ``[mode] enable_spin_orbital = true``) for ``chi0q_mode = "flex"``.

    ``full = True`` (the dynamic solver, issue #284) returns ``chi`` and
    ``chi0`` on the whole bosonic axis, ``(nmat, nvol, nd^2, nd^2)``
    (:func:`_so_flex_full_axis`), and requires the Green function on the
    configured Nmat; the tag checks are the same.

    Returns ``(chi, green_kw, chi0, chi_path, borrowed)``:

    * ``chi`` ``(nvol, nd^2, nd^2)``: the stored ``chiq`` -- the single
      antisymmetrized spinful channel ``[1 + chi0 Gamma]^-1 chi0`` in the
      bubble ``[a, c, b, d]`` pair layout over spin-block spin-orbitals
      (``Flex._flex_compute_veff_spinful``) -- at the zero bosonic frequency;
    * ``green_kw`` ``(nd, nd, Nx, Ny, Nz, nmat)``: the stored dressed Green
      function (spin-block order, as stored; no reorder);
    * ``chi0`` ``(nvol, nd^2, nd^2)``: the stored static bubble of the same
      map (for :func:`_check_so_flex_vertex`). The chi0q archive is REQUIRED:
      ``chiq`` carries no fingerprint of the vertex it was resummed with, so
      without it interaction files that differ from the FLEX run's would give
      silently wrong eigenvalues;
    * ``chi_path``: the chiq archive read;
    * ``borrowed`` (:class:`_SOFlexBorrowed`), or None without ``so_run``.

    With ``so_run`` (:class:`_SORun`, #88) the provenance of the chiq and
    chi0q archives is checked against this run (:func:`_so_check_archive`,
    exact nmat) after their tag checks, and so are the green and sigma
    archives of the same run, all before any payload is read; the chiq,
    chi0q and green archives are checked again with the digest of the
    loaded Green function, and the sigma archive is then read
    (:func:`_load_flex_so_sigma`, which binds it to the same digest).

    The archive is identified by its TAGS, never by its shape (a one-orbital
    spin-orbital file has nd = 2, the shape of a two-orbital non-spin-orbital
    one): ``chi_convention = "spin_orbital"`` and ``index_convention`` equal
    to :data:`hwave.solver.spin_orbital.INDEX_CONVENTION`. A non-spin-orbital
    ("myo"), reduced ("kuroki") or untagged archive is refused before any
    heavy work, and so is a Green function without the spin-block tag and a
    chi0q without the ``spin_orbital = True`` marker (written only by a
    spin-orbital run)."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver import spin_orbital
    from hwave.solver.rpa import (validate_momentum_convention,
                                  check_momentum_marker)
    eli = CaseInsensitiveDict(input_dict.get("eliashberg", {}) or {})
    flex_dir = _resolve_flex_dir(input_dict)
    chi_path = os.path.join(flex_dir, eli.get("flex_chiq", "chiq.npz"))
    chi0_path = os.path.join(flex_dir, eli.get("flex_chi0q", "chi0q.npz"))
    # the Green file name is read case-insensitively like the keys above; the
    # shared loader below is handed a view of the config with that name
    green_name = eli.get("flex_green", "green.npz")
    green_path = os.path.join(flex_dir, green_name)
    green_input = dict(input_dict)
    green_input["eliashberg"] = dict(
        {k: v for k, v in (input_dict.get("eliashberg", {}) or {}).items()
         if str(k).lower() != "flex_green"}, flex_green=green_name)
    nvol = Nx * Ny * Nz
    config_nmat = input_dict.get("mode", {}).get("param", {}).get(
        "Nmat", _DEFAULT_NMAT)
    hint = ("[mode] enable_spin_orbital = true with chi0q_mode = \"flex\" reads "
            "the outputs of a spin-orbital FLEX run (calc_scheme = \"general\", "
            "[mode] enable_spin_orbital = true, [file.output] chiq, chi0q, green "
            "and sigma)")

    if not os.path.exists(chi_path):
        raise ValueError(
            "FLEX susceptibility archive '{}' not found. {}; set [eliashberg] "
            "flex_chiq to the file the FLEX run wrote as [file.output] "
            "chiq.".format(chi_path, hint))
    if not os.path.exists(chi0_path):
        # the stored chiq carries no fingerprint of the vertex it was resummed
        # with; the run's chi0q is what pins Gamma (_check_so_flex_vertex)
        raise ValueError(
            "FLEX chi0q archive '{}' not found. {}, and also its chi0q: the "
            "check that this run's interaction files reproduce the FLEX run's "
            "vertex needs it. Set [file.output] chi0q in the FLEX run and "
            "[eliashberg] flex_chi0q to that file.".format(chi0_path, hint))
    ir_hint = _SO_DYNAMIC_IR_HINT if full else _STATIC_IR_HINT
    logger.info("Loading spin-orbital FLEX chiq from: {}".format(chi_path))
    with np.load(chi_path) as data:
        _reject_ir_native(data, chi_path, ir_hint)
        if "chi_convention" not in data:
            raise ValueError(
                "FLEX archive '{}' has no chi_convention tag: an untagged file "
                "cannot be told from a non-spin-orbital one by its shape (one "
                "spin-orbital orbital has the shape of two orbitals). {}.".format(
                    chi_path, hint))
        conv = str(data["chi_convention"])
        if conv != _SO_FLEX_CHI_CONVENTION:
            what = {"myo": "a non-spin-orbital general FLEX output (separate "
                           "spin and charge channels)",
                    "kuroki": "a reduced FLEX output (separate spin and "
                              "charge channels)"}.get(
                        conv.lower(), "not a spin-orbital FLEX output")
            raise ValueError(
                "FLEX archive '{}' carries chi_convention = {!r}: {}, which "
                "cannot carry the spin-flip pairs of a spin-orbital model. "
                "{}.".format(chi_path, conv, what, hint))
        if "index_convention" not in data:
            raise ValueError(
                "FLEX archive '{}' has no index_convention tag; the spin-orbital "
                "order of its axes cannot be verified. {}.".format(chi_path, hint))
        if str(data["index_convention"]) != spin_orbital.INDEX_CONVENTION:
            raise ValueError(
                "FLEX archive '{}' carries index_convention = {!r}, expected {!r}. "
                "{}.".format(chi_path, str(data["index_convention"]),
                             spin_orbital.INDEX_CONVENTION, hint))
        if "chiq" not in data:
            raise ValueError(
                "FLEX archive '{}' has no 'chiq' array (the spinful channel). "
                "{}.".format(chi_path, hint))
        if ("chi_orbital_layout" in data
                and str(data["chi_orbital_layout"]) != "acbd"):
            raise ValueError(
                "FLEX archive '{}' declares chi_orbital_layout = {!r} but this "
                "reader only supports 'acbd'.".format(
                    chi_path, str(data["chi_orbital_layout"])))
        if "scf_converged" in data and not bool(np.asarray(data["scf_converged"])):
            logger.warning(
                "FLEX archive '%s' records scf_converged = false: the FLEX run "
                "did not converge, and its chiq and green may not belong to one "
                "self-consistent state.", chi_path)
        check_momentum_marker(data, chi_path)
        # the chi0q companion must be a spin-orbital FLEX output: its marker
        # is checked after the lightweight chiq metadata above, but before
        # any array (chiq or chi0q) is materialized and before the Green
        # function or any vertex is built (_check_so_flex_vertex still
        # follows: the marker does not establish equal interactions)
        with np.load(chi0_path) as d0:
            _require_so_chi0q_marker(
                d0, chi0_path, producer="FLEX",
                hint="{}; and its chi0q ([file.output] chi0q).".format(hint))
        if so_run is not None:
            _so_check_archive(data, chi_path, so_run, "exact", "chiq",
                              _SO_FLEX_PRODUCER)
            # the chi0q companion carries the same members (FLEX writes one
            # build into every archive); checked before any array is read
            with np.load(chi0_path) as d0:
                _so_check_archive(d0, chi0_path, so_run, "exact", "chi0q",
                                  _SO_FLEX_PRODUCER)
            # the green and sigma archives of the same run against this run
            # (structure, band, interaction, state), also before any payload
            # is read; their Green digest is compared once the Green function
            # is loaded (a missing file gets its own message below)
            sigma_path = _so_flex_sigma_path(eli, flex_dir)
            for path, kind in ((green_path, "green"), (sigma_path, "sigma")):
                if os.path.exists(path):
                    with np.load(path) as other:
                        _so_check_archive(other, path, so_run, "exact", kind,
                                          _SO_FLEX_PRODUCER)
        chiq = data["chiq"]
        if chiq.ndim >= 2 and chiq.shape[1] == nvol:
            validate_momentum_convention(data, chi_path, chiq, 1, (Nx, Ny, Nz))
        if full:
            chi = _so_flex_full_axis(data, chiq, chi_path, nd, nvol, config_nmat,
                                     "chiq")
        else:
            chi = _so_flex_static_slice(data, chiq, chi_path, nd, nvol,
                                        config_nmat, "chiq")
        del chiq

    # the Green function: required (the stored chi belongs to the dressed G)
    if not os.path.exists(green_path):
        raise ValueError(
            "FLEX Green archive '{}' not found. {}; set [eliashberg] flex_green "
            "to the file the FLEX run wrote as [file.output] green.".format(
                green_path, hint))
    with np.load(green_path) as gdata:
        gtag = (str(gdata["index_convention"]) if "index_convention" in gdata
                else None)
    if gtag != spin_orbital.INDEX_CONVENTION:
        raise ValueError(
            "FLEX Green archive '{}' carries index_convention = {!r}, expected "
            "{!r} (a spin-orbital FLEX run tags it). {}.".format(
                green_path, gtag, spin_orbital.INDEX_CONVENTION, hint))
    green_kw = _load_flex_green(green_input, nd, Nx, Ny, Nz,
                                interactions=interactions)
    if green_kw is None or _resolve_flex_paths(green_input)[2] != green_path:
        raise ValueError(
            "internal error: the FLEX Green archive '{}' resolved to a different "
            "path in the shared loader.".format(green_path))
    if full and green_kw.shape[-1] != int(config_nmat):
        raise ValueError(
            "FLEX Green archive '{}' holds {} fermionic frequencies, but "
            "[mode.param] Nmat = {}: frequency = \"dynamic\" needs the Green "
            "function on the run's grid. Re-run FLEX with the same Nmat.".format(
                green_path, green_kw.shape[-1], int(config_nmat)))

    with np.load(chi0_path) as d0:
        _reject_ir_native(d0, chi0_path, ir_hint)
        validate_chi0q_index_convention(d0, True, chi0_path)
        check_momentum_marker(d0, chi0_path)
        if full:
            chi0 = _so_flex_full_axis(d0, d0["chi0q"], chi0_path, nd, nvol,
                                      config_nmat, "chi0q")
        else:
            chi0 = _so_flex_static_slice(d0, d0["chi0q"], chi0_path, nd, nvol,
                                         config_nmat, "chi0q")
    borrowed = None
    if so_run is not None:
        consumed_fp = _so_green_fp(green_kw)
        with np.load(chi_path) as data:
            so_info = _so_check_archive(data, chi_path, so_run, "exact", "chiq",
                                        _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        with np.load(chi0_path) as d0:
            _so_check_archive(d0, chi0_path, so_run, "exact", "chi0q",
                              _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        with np.load(green_path) as gd:
            _so_check_archive(gd, green_path, so_run, "exact", "green",
                              _SO_FLEX_PRODUCER, green_fp=consumed_fp)
        sigma_static = _load_flex_so_sigma(eli, flex_dir, nd, nvol, so_run,
                                           consumed_fp, green_path, hint)
        borrowed = _SOFlexBorrowed(so_info, consumed_fp, sigma_static)
    return chi, green_kw, chi0, chi_path, borrowed
