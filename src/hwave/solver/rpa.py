from __future__ import annotations
from typing import Optional

import sys, os
import numpy as np
import numpy.fft as FFT
import itertools
import copy
from requests.structures import CaseInsensitiveDict

try:
    from .perf import do_profile
except ImportError:
    from functools import wraps
    def do_profile(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)
        return wrapper

import logging
logger = logging.getLogger(__name__)

#import read_input_k
import hwave.qlmsio.read_input_k as read_input_k
import hwave.qlmsio.wan90 as wan90
from hwave.solver.vertex_table import fierz_coefficients, ring_spin_table
from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.declarations import symmetrise_dense
from hwave.solver.density_projection import project_density_pairs
from . import backend as _bk
from . import bond_channels
from . import bubble
from . import fold
from . import green as green_mod


def validate_chi0q_index_convention(data, enable_spin_orbital, file_name=""):
    """Reject a stored chi0q whose spin-orbital index convention is unknown.

    RPA writes ``index_convention="spin_block"`` (spin*norb+orb) into its
    chi0q.npz / chiq.npz. UHFk instead uses the interleaved (2*orb+spin)
    ordering, and chi0q files produced before the SO convention fix
    (commit 9dd9a21) carry no ``index_convention`` key at all. In spin-orbital
    mode the two orderings differ, so silently accepting such a file would mix
    conventions and corrupt the result. Outside spin-orbital mode the
    convention is irrelevant and this is a no-op.

    Parameters
    ----------
    data : Mapping
        Loaded npz (e.g. ``numpy.load(...)``); must support ``in`` and ``[]``.
    enable_spin_orbital : bool
        Whether the consuming calculation runs in spin-orbital mode.
    file_name : str, optional
        Source path, used only for the error message.

    Raises
    ------
    ValueError
        If spin-orbital mode is on and the file lacks a ``spin_block``
        ``index_convention`` marker.
    """
    if not enable_spin_orbital:
        return
    if "index_convention" not in data:
        raise ValueError(
            "chi0q file '{}' lacks an 'index_convention' marker. It predates "
            "the spin-orbital index convention fix and may use the interleaved "
            "(2*orb+spin) ordering; regenerate it with the current RPA solver "
            "(which writes index_convention='spin_block').".format(file_name)
        )
    conv = str(data["index_convention"])
    if conv != "spin_block":
        raise ValueError(
            "chi0q file '{}' has index_convention='{}', but spin-orbital mode "
            "requires 'spin_block' (spin*norb+orb). Regenerate it with the "
            "current RPA solver.".format(file_name, conv)
        )

MOMENTUM_CONVENTION = "e_plus_ikR"

# Default [mode.param] transverse_bond_memory_cap_gb (Phase W experimental
# gate, spec 2026-08-15-bond-transverse-design.md "Production surface"):
# same default sc.py's [eliashberg] bond_memory_cap_gb uses for its
# longitudinal sibling (_BOND_MEMORY_CAP_GB), so the two bond-channel
# gates share one documented default rather than drifting apart.
TRANSVERSE_BOND_MEMORY_CAP_GB_DEFAULT = 8.0


def check_momentum_marker(data, file_name):
    """Strictly validate a PRESENT momentum_convention marker.

    Returns True when a valid marker is present, False when the file has
    none. A malformed marker (not exactly one value) or a mismatching tag
    raises. Consumers whose legacy files are known to already use the
    documented sign (e.g. UHFk green outputs) call this alone, accepting
    unmarked files without a content check.
    """
    files = getattr(data, "files", data)
    if "momentum_convention" not in files:
        return False
    tag_arr = np.asarray(data["momentum_convention"]).ravel()
    # exactly one value (round-5 review): a multi-element marker
    # previously authorized via its first element and an empty one
    # died on an incidental IndexError
    if tag_arr.size != 1:
        raise ValueError(
            "file '{}': momentum_convention must be a single value, "
            "got {!r}; regenerate the file.".format(
                file_name, np.asarray(data["momentum_convention"])))
    tag = str(tag_arr[0])
    if tag != MOMENTUM_CONVENTION:
        raise ValueError(
            "file '{}' records momentum_convention = '{}' but this "
            "build uses '{}'; regenerate the file.".format(
                file_name, tag, MOMENTUM_CONVENTION))
    return True


# Chunk bound for the validator's scans: temporaries per processed block
# stay at a few times this many elements, independent of the payload size
# (round-9 review: the previous whole-tensor pipeline allocated ~5.5x the
# payload -- ~200 GiB for a large general-scheme chi file).
_MOMENTUM_SCAN_BLOCK = 1 << 18


def validate_momentum_convention(data, file_name, payload, q_axis,
                                 lattice_shape):
    """Fail closed on momentum-convention mismatches (issue #133).

    Every k/q-space NPZ written since #133 carries
    momentum_convention = "e_plus_ikR" (the documented Wannier90-style
    sign). Files written before the fix carry q labels of the OPPOSITE
    sign; silently combining one with current-convention quantities would
    mix q and -q. Legacy files without the marker are accepted ONLY when
    the stored payload is elementwise even under q -> -q on the FFT grid
    (then the two conventions coincide bit-for-bit -- true for every
    centrosymmetric fixture); otherwise they are rejected with a
    regeneration hint, following the sc_vertex_version fail-closed
    precedent of gating on content, not age.

    Both the finiteness scan (which runs for marked files too -- a valid
    marker must not authorize NaN/Inf content) and the q/-q comparison
    are CHUNKED: no full-size reversed, normalized, or boolean tensor is
    materialized (round-9 review).

    Parameters: payload is the stored array, q_axis its flattened
    (nx, ny, nz) momentum axis (or a 3-tuple of explicit axes),
    lattice_shape that grid.
    """
    marked = check_momentum_marker(data, file_name)
    arr = np.asarray(payload)
    if arr.size == 0:
        return
    nx, ny, nz = lattice_shape
    nvol = nx * ny * nz
    from hwave.solver.kgrid import reverse_fft_axes

    def _num(block):
        # signed-integer payloads overflow np.abs at INT_MIN (round-10
        # review); cast PER BLOCK to float64 -- bounded, and a no-op for
        # the float/complex arrays production writes
        b = np.asarray(block)
        if b.dtype.kind not in "fc":
            b = b.astype(np.float64)
        return b

    # Two-level chunking (round-10 review): chunk the q axis AND, when a
    # single q plane alone exceeds the element bound (large frequency/
    # orbital axes, tiny grids), also the leading non-q axis. Extraction
    # always indexes the ORIGINAL array (slices / fancy indexing copy
    # only the selected block; operating on a moveaxis view made numpy
    # consolidate the whole tensor). Reversed q indices are computed per
    # chunk -- no full-size index grids.
    def _rev_flat(idx):
        ix = idx // (ny * nz)
        iy = (idx // nz) % ny
        iz = idx % nz
        return (((-ix) % nx) * ny + ((-iy) % ny)) * nz + ((-iz) % nz)

    if isinstance(q_axis, (tuple, list)):
        axes = [a % arr.ndim for a in q_axis]
        dims = [nx, ny, nz]
        o = int(np.argmax(dims))
        outer_ax, n_outer = axes[o], dims[o]
        inner_block_axes = tuple(a - (a > outer_ax)
                                 for i, a in enumerate(axes) if i != o)
        q_set = set(axes)
    else:
        outer_ax = q_axis % arr.ndim
        n_outer = arr.shape[outer_ax]
        inner_block_axes = None
        q_set = {outer_ax}
    # lead-chunk over an axis OUTSIDE the momentum set (round-11: slicing
    # an inner q axis would break the reversal pairing); if every axis is
    # a momentum axis, lead chunking is disabled
    chunk_ax = next((i for i in range(arr.ndim) if i not in q_set), None)
    lead = arr.shape[chunk_ax] if chunk_ax is not None else 1
    plane = max(1, arr.size // n_outer)

    if plane <= _MOMENTUM_SCAN_BLOCK:
        q_step = max(1, _MOMENTUM_SCAN_BLOCK // plane)
        lead_step = lead
    else:
        q_step = 1
        per_lead = max(1, plane // lead)
        lead_step = max(1, _MOMENTUM_SCAN_BLOCK // per_lead)

    def _lead_slices():
        if chunk_ax is None:
            yield slice(None)
        else:
            for l0 in range(0, lead, lead_step):
                yield slice(l0, l0 + lead_step)

    def _index(lsl, qsel):
        idx = [slice(None)] * arr.ndim
        if chunk_ax is not None:
            idx[chunk_ax] = lsl
        idx[outer_ax] = qsel
        return arr[tuple(idx)]

    def _finite_blocks():
        for lsl in _lead_slices():
            for q0 in range(0, n_outer, q_step):
                yield _index(lsl, slice(q0, q0 + q_step))

    def _pairs():
        for lsl in _lead_slices():
            if inner_block_axes is None:
                for q0 in range(0, n_outer, q_step):
                    idx = np.arange(q0, min(q0 + q_step, n_outer))
                    yield (_index(lsl, idx),
                           _index(lsl, _rev_flat(idx)))
            else:
                for ix in range(n_outer):
                    a = np.asarray(_index(lsl, ix))
                    b = np.asarray(_index(lsl, (-ix) % n_outer))
                    yield a, reverse_fft_axes(b, inner_block_axes)

    # chunked finiteness + component-scale pass; runs for MARKED files
    # too (a marker must not authorize NaN/Inf content, round-9 review)
    scale = 0.0
    for block in _finite_blocks():
        b = _num(block)
        if not np.all(np.isfinite(b.real)) or (
                np.iscomplexobj(b) and not np.all(np.isfinite(b.imag))):
            raise ValueError(
                "file '{}' contains non-finite values; regenerate the "
                "file.".format(file_name))
        if np.iscomplexobj(b):
            m = max(float(np.abs(b.real).max()),
                    float(np.abs(b.imag).max()))
        else:
            m = float(np.abs(b).max())
        scale = max(scale, m)
    if marked or scale == 0.0:
        return
    norm = max(scale, 1.0e-300)

    for a, b in _pairs():
        an = _num(a) / norm
        bn = _num(b) / norm
        diff = np.abs(an - bn)
        tol = 1.0e-8 * (np.abs(an) + np.abs(bn)) + 1.0e-15
        if bool(np.any(diff > tol)):
            asym_rel = float(diff.max())
            raise ValueError(
                "file '{}' carries no momentum_convention marker and its "
                "content is not elementwise even under q -> -q (max pair "
                "deviation {:.3e} relative to the global scale {:.3e}): "
                "it was written before the #133 Fourier-sign alignment "
                "and its momentum labels are flipped relative to this "
                "build. Regenerate the file with the current version. "
                "(Files whose content is q-even are accepted: for them "
                "the two conventions coincide.)".format(
                    file_name, asym_rel, scale))



def _so_physical_norb(geom_norb, enable_spin_orbital, *, check_norb=None,
                      source="geom.dat"):
    """Physical orbital count from a geometry ``norb``.

    In spin-orbital mode ``geom.dat``'s ``norb`` is the spin-orbital count
    (= 2 * physical orbitals = Wannier90 num_wann), matching UHFk, so halve it.
    Evenness is validated on ``check_norb`` (the *original*, pre-sublattice-fold
    value) so the error names the user's actual ``geom.dat`` entry, while the
    returned count is derived from ``geom_norb`` (which may be the post-fold
    value ``orig * subvol``).
    """
    if not enable_spin_orbital:
        return geom_norb
    cn = check_norb if check_norb is not None else geom_norb
    if cn % 2 != 0:
        raise ValueError(
            "spin-orbital mode requires an even Geometry norb (the spin-orbital "
            "count = 2 * physical orbitals); got {} in {}".format(cn, source))
    return geom_norb // 2


def _chi0q_fingerprint(arr):
    """Content binding for chi0q provenance.

    Provenance metadata travels in a caller-owned dict, so a caller can
    replace the array while keeping the metadata; shape checks alone
    cannot see a same-shaped replacement. The FULL canonical byte content
    is hashed -- a strided sample was shown to miss single-element edits
    of unsampled positions (round-6 review). blake2b throughput makes the
    cost linear but small next to the solve itself; a non-contiguous view
    or a device-backed array is converted first, which may copy.
    Legitimate copies hash equal.
    """
    import hashlib

    from . import backend as _bk_local

    a = np.ascontiguousarray(_bk_local.to_host(arr))
    h = hashlib.blake2b(digest_size=16)
    h.update(str(a.shape).encode())
    h.update(str(a.dtype).encode())
    h.update(memoryview(a).cast("B"))
    return h.hexdigest()


def _validate_chi0q_provenance(meta, nfreq, source):
    """Validate and NORMALIZE a chi0q provenance mapping.

    Shared by the in-memory reuse route and the chi0q_init file route so
    neither can smuggle malformed provenance past the other (round-5
    review: nmat was validated only when freq_index was present, a
    fractional nmat was written back unnormalized, and a non-mapping
    value raised AttributeError instead of the promised ValueError).

    Returns {"freq_index": int ndarray or None, "nmat": int or None,
    "coeff_tail": float or None, "tail_endpoint": str or None}. Raises
    ValueError on any defect.
    """
    import operator
    from collections.abc import Mapping

    if not isinstance(meta, Mapping):
        raise ValueError(
            "chi0q provenance from {} must be a mapping, got {}".format(
                source, type(meta).__name__))
    out = {"freq_index": None, "nmat": None, "coeff_tail": None,
           "tail_endpoint": None}
    m_nmat = meta.get("nmat")
    if m_nmat is not None:
        # strict integral scalar: operator.index accepts int and NumPy
        # integer scalars (including the 0-d arrays npz files store) and
        # rejects floats, strings, complex, containers and booleans --
        # every one of which slipped past size-1 coercion (round 6)
        try:
            val = np.asarray(m_nmat)
            if val.ndim != 0:
                raise TypeError("not a scalar")
            item = val[()]
            if isinstance(item, (bool, np.bool_)):
                raise TypeError("boolean")
            n = operator.index(item)
        except TypeError:
            raise ValueError(
                "chi0q provenance from {}: nmat must be an integral "
                "scalar, got {!r}".format(source, m_nmat))
        if n <= 0:
            raise ValueError(
                "chi0q provenance from {}: nmat must be positive, got "
                "{}".format(source, n))
        out["nmat"] = n
    fi = meta.get("freq_index")
    if fi is not None:
        fi = np.asarray(fi)
        if fi.ndim != 1 or not np.issubdtype(fi.dtype, np.integer):
            raise ValueError(
                "chi0q provenance from {}: freq_index must be a "
                "one-dimensional integer array, got dtype {} with shape "
                "{}".format(source, fi.dtype, fi.shape))
        if len(fi) != nfreq:
            raise ValueError(
                "chi0q provenance from {} does not describe the supplied "
                "chi0q (freq_index of length {} for a {}-frequency "
                "array); the metadata is stale -- recompute or drop "
                "it.".format(source, len(fi), nfreq))
        if len(np.unique(fi)) != len(fi):
            raise ValueError(
                "chi0q provenance from {}: freq_index carries duplicate "
                "entries".format(source))
        if (out["nmat"] is not None and len(fi) > 0
                and (fi.min() < 0 or fi.max() >= out["nmat"])):
            raise ValueError(
                "chi0q provenance from {}: freq_index range [{}, {}] "
                "exceeds nmat = {}".format(
                    source, int(fi.min()), int(fi.max()), out["nmat"]))
        out["freq_index"] = fi
    ct = meta.get("coeff_tail")
    if ct is not None:
        # type-strict like the config validation (round-8 review):
        # float() would also normalize booleans, numeric strings and
        # non-finite values, letting malformed foreign provenance pass
        # the endpoint gate and be re-saved. Unwrap the 0-d scalar an
        # npz file stores, then require a finite non-boolean Real.
        import numbers
        val = np.asarray(ct)
        item = val[()] if val.ndim == 0 else None
        if (item is None or isinstance(item, (bool, np.bool_))
                or not isinstance(item, numbers.Real)
                or not np.isfinite(float(item))):
            raise ValueError(
                "chi0q provenance from {}: malformed coeff_tail "
                "({!r})".format(source, ct))
        out["coeff_tail"] = float(item)
    te = meta.get("tail_endpoint")
    if te is not None:
        # npz stores strings as 0-d unicode arrays; unwrap and require a
        # plain string (np.str_ passes the isinstance check)
        val = np.asarray(te)
        item = val[()] if val.ndim == 0 else None
        if not isinstance(item, str):
            raise ValueError(
                "chi0q provenance from {}: malformed tail_endpoint "
                "({!r})".format(source, te))
        out["tail_endpoint"] = str(item)
    return out


# Equal-time endpoint convention of the coeff_tail machinery (issue #134).
# chi0q computed with a NONZERO coeff_tail before the branch-mean endpoint
# fix carries an O(1/Nmat) error that a current recomputation at the same
# Nmat does not; files stamp this marker so loaders can tell the two
# apart. Bump the value if the endpoint treatment ever changes again.
TAIL_ENDPOINT_CONVENTION = "branch_mean_v1"


def _to_bubble_pair_convention(ham):
    """Intra-pair transpose taking the interaction tensor from the
    Hamiltonian slot convention to the bubble's (issue #139).

    ``_make_ham_inter`` stores ``W[..., b, b', a, a']`` as the
    coefficient of ``c^+_a c_a' c^+_b' c_b``; the ring equation
    contracts it against ``chi0`` whose pair slot ``(b, b')`` carries
    ``c^+_b c_b'``. Swapping the two indices of BOTH pairs converts
    between them. Density (pair-diagonal) content and real
    Hermitian-closed declarations are fixed points of this map.

    Applied when the LONGITUDINAL vertex is assembled, i.e. on the
    rank-4 orbital tensor before any density projection: the reduced
    projection keeps only pair-diagonal slots, where the
    map is the identity, so it is unaffected either way, while the
    general scheme needs it. The transverse (ladder) assembly is NOT
    routed through this helper: it re-pairs the tensor itself and
    therefore consumes the Hamiltonian convention directly.
    """
    nlead = ham.ndim - 4
    if nlead < 0:
        raise ValueError(
            "interaction tensor must carry four orbital axes, got shape "
            "{}".format(ham.shape))
    lead = tuple(range(nlead))
    return ham.transpose(*lead, nlead + 1, nlead, nlead + 3, nlead + 2)


def _enforce_tail_endpoint(meta, source):
    """Endpoint-convention gate on NORMALIZED chi0q provenance.

    Fail-closed for any nonzero-coeff_tail bubble whose endpoint
    treatment is absent (pre-#134) or unrecognized -- nothing in the
    array itself can reveal which treatment produced it. Zero-tail
    bubbles and unknown-tail (no coeff_tail record) bubbles pass:
    the endpoint is vacuous respectively unknowable there. One gate is
    shared by the chi0q_init file reader and the in-memory reuse route
    (round-7 review: the in-memory route bypassed the file gate).
    """
    tail = meta.get("coeff_tail")
    if tail is None or tail == 0.0:
        return
    te = meta.get("tail_endpoint")
    if te is None:
        raise ValueError(
            "chi0q provenance from {} records coeff_tail = {} but no "
            "tail_endpoint marker: the bubble was produced before the "
            "equal-time endpoint fix (issue #134) and its tail-corrected "
            "values carry the pre-fix O(1/Nmat) endpoint error. Recompute "
            "the bubble with this version instead of reusing it.".format(
                source, tail))
    if te != TAIL_ENDPOINT_CONVENTION:
        raise ValueError(
            "chi0q provenance from {} carries unrecognized tail_endpoint "
            "= {!r} (this build implements {!r}); refusing to reuse a "
            "bubble whose endpoint treatment is unknown.".format(
                source, te, TAIL_ENDPOINT_CONVENTION))


class Lattice:
    """
    Lattice parameters:

      Lattice.cellshape = (Lx, Ly, Lz)
      Lattice.cellvol
        original box shape and its volume

      Lattice.subshape = (Bx, By, Bz)
      Lattice.subvol
        shape of supercell and its volume

      Lattice.shape = (Nx, Ny, Nz)
      Lattice.nvol
        box shape in units of supercell and its volume

      Lattice.has_sublattice = True | False
        whether sublattice (except equiv to original lattice) is defined

    Constructor:
      Lattice(param) where param contains
        CellShape (mandatory)
        SubShape  (optional: default to same as CellShape)

    Note:
      Bx, By, Bz must be chosen so that Lx, Ly, Lz are divisable 
      by Bx, By, Bz, respectively. otherwise, initialization fails.

    """
    def __init__(self, params):
        logger.debug(">>> Lattice.__init__")

        self._init_lattice(params)
        self._show_params()

    def _init_lattice(self, params):
        logger.debug(">>> Lattice._init_lattice")

        if "CellShape" not in params:
            logger.error("Lattice initialization failed: 'CellShape' not found.")
            sys.exit(1)

        cell = params.get("CellShape")
        if type(cell) is not list:
            cell = [ cell ]
        cell_len = len(cell)
        if cell_len < 1 or cell_len > 3:
            logger.error("dimension of CellShape must be one, two, or three.")
            sys.exit(1)
        if cell_len < 3:
            cell.extend([1] * (3 - cell_len))

        Lx,Ly,Lz = cell
        self.cellshape = (Lx,Ly,Lz)
        self.cellvol = Lx * Ly * Lz
        self.celldim = cell_len

        if self.cellvol == 0:
            logger.error("invalid CellShape.")
            sys.exit(1)

        subcell = params.get("SubShape", [Lx,Ly,Lz])
        if type(subcell) is not list:
            subcell = [ subcell ]
        if len(subcell) != cell_len:
            logger.error("dimension of SubShape does not match with that of CellShape.")
            sys.exit(1)
        if len(subcell) < 3:
            subcell.extend([1] * (3 - len(subcell)))

        Bx,By,Bz = subcell
        self.subshape = (Bx,By,Bz)
        self.subvol = Bx * By * Bz

        if self.subvol == 0:
            logger.error("invalid SubShape.")
            sys.exit(1)

        self.has_sublattice = (self.subvol > 1)

        # check consistency
        # XXX use reciprocal lattice
        if not all([ self.cellshape[i] % self.subshape[i] == 0 for i in range(3) ]):
            logger.error("SubShape is not compatible with CellShape.")
            sys.exit(1)

        # replace by lattice of supercells
        nx, ny, nz = Lx//Bx, Ly//By, Lz//Bz
        nvol = nx * ny * nz

        self.shape = (nx, ny, nz)
        self.nvol = nvol

    def _show_params(self):
        logger.info("Lattice parameters:")
        logger.info("    CellShape       = {}".format(self.cellshape))
        logger.info("    cell volume     = {}".format(self.cellvol))
        logger.info("    cell dimension  = {}".format(self.celldim))
        logger.info("    SubShape        = {}".format(self.subshape))
        logger.info("    subshape volume = {}".format(self.subvol))
        logger.info("    Shape           = {}".format(self.shape))
        logger.info("    shape volume    = {}".format(self.nvol))
        logger.info("    has_sublattice  = {}".format(self.has_sublattice))

class Interaction:
    """
    Construct Hamiltonian from input
    """
    def __init__(self, lattice, param_ham, info_mode):
        logger.debug(">>> Interaction.__init__")

        self.lattice = lattice
        self.param_ham = param_ham

        self._has_interaction = False
        self._has_interaction_coulomb = False
        self._has_interaction_exchange = False
        self._has_interaction_pairhop = False


        # mode options
        self.enable_spin_orbital = info_mode.get("enable_spin_orbital", False)

        # initialize, and reshape if use sublattice
        self._init_interaction()

        # geom norb is the spin-orbital count in SO mode (UHFk/W90 convention).
        # _init_interaction may have folded the geometry (norb -> norb*subvol),
        # so read the POST-fold value here; validate evenness on the pre-fold
        # original (param_ham_orig exists only when has_sublattice).
        post_fold_norb = param_ham["Geometry"]["norb"]
        if self.lattice.has_sublattice:
            orig_norb = self.param_ham_orig["Geometry"]["norb"]
        else:
            orig_norb = post_fold_norb
        self.norb = _so_physical_norb(post_fold_norb, self.enable_spin_orbital,
                                      check_norb=orig_norb, source="Geometry")
        # Pre-fold physical orbital count (per ORIGINAL cell). Equals self.norb
        # without sublattice; under folding self.norb = norb_orig * subvol.
        # Reused by RPA._reshape_green to decode the folded orbital index.
        self.norb_orig = _so_physical_norb(orig_norb, self.enable_spin_orbital)

        # create hamiltonian
        self._make_ham_trans()
        self._make_ham_inter()

        pass

    def has_interaction(self):
        return self._has_interaction

    def has_interaction_coulomb(self):
        return self._has_interaction_coulomb

    def has_interaction_exchange(self):
        return self._has_interaction_exchange

    def has_interaction_pairhop(self):
        return self._has_interaction_pairhop

    def _init_interaction(self):
        logger.debug(">>> Interaction._init_interaction")

        # reinterpret interaction coefficient on sublattice
        if self.lattice.has_sublattice:
            # backup
            self.param_ham_orig = copy.deepcopy(self.param_ham)

            # replace by sublatticed versions
            # dispatch on the NORMALIZED name (PR #128 round 7, fifth
            # surfacing of the case-defect class): param_ham is a
            # CaseInsensitiveDict that PRESERVES the declared case, so a
            # lowercase 'geometry' key fell into _reshape_interaction and
            # crashed unpacking the geometry table
            for type in self.param_ham.keys():
                if type.lower() == "initial":
                    pass
                elif type.lower() == "geometry":
                    tbl = self._reshape_geometry(self.param_ham[type])
                    self.param_ham[type] = tbl
                elif type.lower() == "transfer":
                    tbl = self._reshape_interaction(self.param_ham[type], self.enable_spin_orbital)
                    self.param_ham[type] = tbl
                else:
                    tbl = self._reshape_interaction(self.param_ham[type], False)
                    self.param_ham[type] = tbl
        pass

    def _reshape_geometry(self, geom):
        logger.debug(">>> Interaction._reshape_geometry")
        return fold.reshape_geometry(geom, self.lattice.subshape)

    def _reshape_interaction(self, ham, enable_spin_orbital):
        logger.debug(">>> Interaction._reshape_interaction")

        # In SO mode, geom norb is the spin-orbital count; two-body
        # interactions use physical orbital indices, so the stride for the
        # non-SO fold is the physical count (see fold.reshape_interaction,
        # the single definition shared with UHFk).
        geom_norb_orig = self.param_ham_orig["Geometry"]["norb"]
        norb_phys_orig = _so_physical_norb(geom_norb_orig, self.enable_spin_orbital)

        return fold.reshape_interaction(
            ham, self.lattice.subshape, self.lattice.shape,
            norb_so_orig=geom_norb_orig,
            norb_phys_orig=norb_phys_orig,
            enable_spin_orbital=enable_spin_orbital)

    def _export_interaction(self, type, file_name):
        logger.debug(">>> Interaction._export_interaction")

        intr = self.param_ham[type]

        min_r = [0,0,0]
        max_r = [0,0,0]
        for (irvec,orbvec), v in self.param_ham[type].items():
            for k in range(3):
                min_r[k] = irvec[k] if irvec[k] < min_r[k] else min_r[k]
                max_r[k] = irvec[k] if irvec[k] > max_r[k] else max_r[k]
        rshape = [ max_r[i]-min_r[i]+1 for i in range(3) ]
        rsize = rshape[0] * rshape[1] * rshape[2]

        with open(file_name, "w") as fw:
            # write header
            fw.write("{} with sublattice for uhfk\n".format(type))
            # write number of orbitals
            fw.write("{}\n".format(self.norb))
            # write number of points of box enclosing transport vectors
            fw.write("{}\n".format(rsize))
            # write multiplicity factors (nominal)
            for i in range(rsize):
                if i > 0 and i % 15 == 0:
                    fw.write("\n")
                fw.write(" 1")
            fw.write("\n")
            # write index and elements
            for (irvec,orbvec), v in self.param_ham[type].items():
                if (abs(v) > 1.0e-12):
                    fw.write("{:3} {:3} {:3} {:3} {:3}  {:.12f} {:.12f}\n".format(
                        *irvec, orbvec[0]+1, orbvec[1]+1, v.real, v.imag
                    ))

    def _make_ham_trans(self):
        logger.debug(">>> Interaction._make_ham_trans")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        if 'Transfer' not in self.param_ham.keys():
            logger.warning("Transfer not found")
            self.ham_trans_r = None
            self.ham_trans_q = None
            self.ham_extern_r = None
            self.ham_extern_q = None
            return

        if self.enable_spin_orbital == True:
            # The SO transfer file uses the interleaved convention
            # (index = 2*orb + spin, matching UHFk and the docs), while RPA
            # works internally in spin-block order (index = spin*norb + orb).
            # Remap each orbital index P(i) = (i % 2) * norb + i // 2 on both
            # the row and column so the (spin, orbital) reshapes downstream are
            # correct. For norb_phys = 1 this is the identity.
            def _so_interleaved_to_spinblock(i):
                return (i % 2) * norb + i // 2

            tab_r = np.zeros((nx,ny,nz,nd,nd), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                if not (0 <= orbvec[0] < nd and 0 <= orbvec[1] < nd):
                    raise ValueError(
                        "spin-orbital Transfer index {} out of range [0,{}); "
                        "geom norb (SO count) must cover all transfer indices"
                        .format(orbvec, nd))
                a = _so_interleaved_to_spinblock(orbvec[0])
                b = _so_interleaved_to_spinblock(orbvec[1])
                # accumulate: R and R+-L are distinct bonds that wrap onto
                # the same array slot when the hopping range reaches the
                # lattice size (numpy negative indexing)
                tab_r[(*irvec, a, b)] += v

            # Fourier transform
            tab_q = FFT.ifftn(tab_r, axes=(0,1,2)) * nvol

            self.ham_trans_r = tab_r.reshape(nvol,nd,nd)
            self.ham_trans_q = tab_q.reshape(nvol,nd,nd)

        else:
            tab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Transfer"].items():
                if orbvec[0] < norb and orbvec[1] < norb:
                    # accumulate: wrap-around R vectors share the array slot
                    tab_r[(*irvec,*orbvec)] += v
                else:
                    pass  # skip spin dependence

            # Fourier transform: the documented convention (issue #133)
            # is eps(k) = sum_R e^{+ikR} t(R) -- ifftn * nvol, matching
            # UHFk and the spin-orbital branch above. The previous fftn
            # computed eps(-k), i.e. every k/q label of this mode's output
            # was globally negated for non-centrosymmetric models.
            tab_q = FFT.ifftn(tab_r, axes=(0,1,2)) * nvol

            # N.B. spin degree of freedom not included
            self.ham_trans_r = tab_r.reshape(nvol,norb,norb)
            self.ham_trans_q = tab_q.reshape(nvol,norb,norb)

        logger.debug("ham_trans_r shape={}, size={}, nonzero_count={}".format(
            self.ham_trans_r.shape,
            self.ham_trans_r.size,
            self.ham_trans_r[abs(self.ham_trans_r) > 1.0e-8].size,
        ))
        logger.debug("ham_trans_q shape={}, size={}, nonzero_count={}".format(
            self.ham_trans_q.shape,
            self.ham_trans_q.size,
            self.ham_trans_q[abs(self.ham_trans_q) > 1.0e-8].size,
        ))

        if 'Extern' in self.param_ham.keys():
            logger.info("read External field")

            hab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Extern"].items():
                if orbvec[0] < norb and orbvec[1] < norb:
                    # accumulate: wrap-around R vectors share the array slot
                    hab_r[(*irvec,*orbvec)] += v
                else:
                    pass  # skip spin dependence

            # Fourier transform: e^{+ikR}, same convention as the
            # transfer term (issue #133)
            hab_q = FFT.ifftn(hab_r, axes=(0,1,2)) * nvol

            # N.B. spin degree of freedom not included
            self.ham_extern_r = hab_r.reshape(nvol,norb,norb)
            self.ham_extern_q = hab_q.reshape(nvol,norb,norb)
        else:
            self.ham_extern_r = None
            self.ham_extern_q = None


    def _make_ham_inter(self):
        logger.debug(">>> Interaction._make_ham_inter")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        norb = self.norb
        ns = 2
        nd = norb * ns

        # Interaction Hamiltonian W[r,b,bp,a,ap]
        #   H = W(r)^{\beta\beta^\prime\alpha\alpah^\prime}
        #        * c_{i\alpha}^\dagger c_{i\alpha^\prime} c_{j\beta^\prime}^\dagger c_{j\beta}
        ham_r = np.zeros((nx,ny,nz,*(ns,norb)*4), dtype=np.complex128)
        # Longitudinal Fierz (cross-slot) corrections live in a SEPARATE
        # tensor: the ring's longitudinal solve consumes ham_inter_q +
        # ham_fierz_q, while the transverse (ladder) assembly reads
        # ham_inter_q alone. The two channels were adjudicated against exact
        # diagonalization independently (#105 transverse, #113 longitudinal),
        # and the cross-spin Exchange correction would otherwise land in the
        # very tensor block the transverse assembly reads via the crossing
        # relation, doubling the adjudicated ladder vertex (measured: the
        # extracted Exchange coupling went -0.7 -> -1.4 when the correction
        # shared the tensor).
        fierz_r = np.zeros_like(ham_r)

        # spin(a,ap,bp,b)  0: up, 1: down -- derived per type from the
        # adjudicated vertex table (density slots via the channel
        # decomposition, spin-flip slots from the #105 transverse data)
        # instead of hand-coded here; hwave.solver.vertex_table is the one
        # source of this content.
        spin_table = {t: ring_spin_table(t)
                      for t in ('CoulombIntra', 'CoulombInter', 'Hund',
                                'Ising', 'PairLift', 'Exchange', 'PairHop')}

        # coulomb-type interactions
        def _symmetrised(type, tbl):
            # Reduce each declaration to its physical symmetric coefficient,
            # the same reading as uhfk.py and hwave.sc (#106/#113/#114): the
            # two declarations of one bond are (R, a, b) and (-R, b, a), and
            # for every type except PairHop they multiply the SAME operator,
            # so the table entering the vertex is the mean
            #     T~[r, a, b] = (T[r, a, b] + T[-r, b, a]) / 2.
            # Without this the ring read a one-sided off-site declaration as
            # v e^{+iqR} (the documented sign, #133) where the operator it
            # declares -- v n_a(i) n_a(i+R),
            # even in R by the site sum -- has the exact vertex v cos(qR)
            # (measured: chiq differed by 1.2e-2 from the symmetric reading
            # of the same Hamiltonian). PairHop's partner is the HERMITIAN
            # entry, so its mean conjugates and its complex phase survives.
            #
            # The reversal is done on a dense (nx, ny, nz) array with the
            # same shared FFT-grid reversal uhfk.py uses (kgrid, index
            # i -> (-i) mod n), NOT by a
            # sign-flipped dictionary-key lookup: table keys may sit in a
            # wrapped canonical form ((n-1, 0, 0) for a -x bond, and folded
            # tables in particular store canonicalized displacements), where
            # a (-R) key lookup would silently miss the partner and halve
            # the coefficient.
            arr = np.zeros((nx, ny, nz, norb, norb), dtype=np.complex128)
            for (irvec, orbvec), v in tbl.items():
                arr[(irvec[0] % nx, irvec[1] % ny, irvec[2] % nz,
                     *orbvec)] += v
            sym = symmetrise_dense(arr, hermitian=(type == "PairHop"))
            out = {}
            for ix, iy, iz, a, b in zip(*np.nonzero(sym)):
                out[((int(ix), int(iy), int(iz)), (int(a), int(b)))] = \
                    sym[ix, iy, iz, a, b]
            return out

        def _append_inter(type, tbl=None):
            logger.debug("_append_inter {}".format(type))
            spins = spin_table[type]
            if tbl is None:
                tbl = self.param_ham[type]
            tbl = _symmetrised(type, tbl)
            for (irvec,orbvec), v in tbl.items():
                a, b = orbvec
                for spinvec, w in spins.items():
                    s1,s2,s3,s4 = spinvec
                    # beta beta' alpha alpha'
                    orb = (s4, b, s3, b, s1, a, s2, a)
                    ham_r[(*irvec, *orb)] += v * w

        # pairhop-type interaction
        #   H^PH = P^{\alpha\alpha^\prime}
        #        * c_{i\alpha\up}^\dagger c_{j\alpha^\prime\up}
        #            c_{i\alpha\down}^\dagger c_{j\alpha^\prime\down}
        #        + (up <-> down)
        def _append_pairhop(type):
            spins = spin_table[type]

            # Off-site PairHop physics is not implemented (spec
            # "Deferred (recorded)": "Off-site PairHop physics (warning +
            # tracking issue only)."). Only the on-site (irvec=(0,0,0))
            # part is read below; warn LOUDLY instead of silently
            # discarding the rest, naming the ORIGINAL (pre-fold) user
            # declarations. Locality is judged on the PRE-fold table
            # (self.param_ham_orig when the lattice has a sublattice,
            # else self.param_ham) -- the same pre-fold locality rule
            # _append_onsite_direct/_append_inter_cross above apply, for
            # the same reason: reading it off the FOLDED self.param_ham
            # table could let an off-site bond folded onto r=(0,0,0)
            # between supercell orbitals escape detection here.
            has_sub = getattr(self.lattice, "has_sublattice", False)
            orig_tbl = (self.param_ham_orig.get(type, {}) if has_sub
                        else self.param_ham.get(type, {}))
            offsite = [(irvec, orbvec, v)
                       for (irvec, orbvec), v in orig_tbl.items()
                       if tuple(int(x) for x in irvec) != (0, 0, 0)]
            if offsite:
                shown = offsite[:8]
                logger.warning(
                    "PairHop declares %d off-site term(s) (irvec != "
                    "(0,0,0)) that this solver silently discards: only "
                    "the on-site part of PairHop is represented (off-site "
                    "PairHop physics is not implemented; see issue #157). "
                    "Restrict PairHop to on-site declarations "
                    "if this is unintended. Declarations dropped: %s%s",
                    len(offsite),
                    "; ".join(
                        "irvec={} orb=({},{}) v={}".format(
                            tuple(int(x) for x in irv), ov[0], ov[1], v)
                        for irv, ov, v in shown),
                    " ... ({} more)".format(len(offsite) - 8)
                    if len(offsite) > 8 else "")

            tbl = _symmetrised(type, self.param_ham[type])
            for (irvec,orbvec), v in tbl.items():
                # take account of same-site interaction only
                if (irvec == (0,0,0)):
                    a, b = orbvec
                    for spinvec, w in spins.items():
                        s1,s2,s3,s4 = spinvec
                        # beta beta' alpha alpha'
                        orb = (s4, b, s3, a, s1, a, s2, b)
                        ham_r[(*irvec, *orb)] += v * w

        # Pre-fold on-site direct tensor feeding the spinful exchange
        # crossing (issue #137). Locality is judged on the PRE-fold
        # declarations, exactly like the Fierz builder below: sublattice
        # folding maps off-site bonds to r = 0 between supercell
        # orbitals, and reading ham_r[0, 0, 0] after folding would cross
        # off-site content that the ring-form resummation cannot
        # represent (round-1 review; the same trap the Fierz builder
        # documents). Only built for spin-orbital runs -- no other mode
        # consumes it.
        # (getattr: tests drive _make_ham_inter on __new__-built stubs
        # without __init__, the same pattern save_results documents)
        onsite_r = (np.zeros((*(ns, norb) * 4,), dtype=np.complex128)
                    if getattr(self, "enable_spin_orbital", False)
                    else None)

        def _append_onsite_direct(type, tbl=None, pairhop=False):
            if onsite_r is None:
                return
            has_sub = getattr(self.lattice, "has_sublattice", False)
            if tbl is None:
                if has_sub:
                    tbl = self.param_ham_orig.get(type, {})
                else:
                    tbl = self.param_ham.get(type, {})
            filtered = {}
            for (irvec, orbvec), v in tbl.items():
                if tuple(irvec) == (0, 0, 0):
                    filtered[(irvec, orbvec)] = v
            if not filtered:
                return
            if has_sub:
                filtered = self._reshape_interaction(filtered, False)
            filtered = _symmetrised(type, filtered)
            spins = spin_table[type]
            for (irvec, orbvec), v in filtered.items():
                if tuple(irvec) != (0, 0, 0):
                    continue
                a, b = orbvec
                for spinvec, w in spins.items():
                    s1, s2, s3, s4 = spinvec
                    # same slot layouts as the ham_r builders above
                    orb = ((s4, b, s3, a, s1, a, s2, b) if pairhop
                           else (s4, b, s3, b, s1, a, s2, a))
                    onsite_r[orb] += v * w

        def _append_inter_cross(type, tbl=None):
            # Longitudinal cross-slot (Fierz) content, adjudicated against
            # exact diagonalization in #113: the ring's W carried the
            # density (aa,bb) slots of each type but not the (ab,ab)
            # pair-cross slots, which is why its effective spin vertex was
            # diag(U0, 0, 0, U1) where the adjudicated table (and the
            # spin/charge builders in hwave.sc) have diag(U0, U', U', U1)
            # (#104; measured 1e-2 in chiq). Restricted to ON-SITE entries
            # with a != b: that is the adjudicated domain -- off-site
            # cross-orbital content is not representable by a q-only vertex
            # (#105 measurement) and stays as it was.
            # `tbl`, when given, must be a PRE-fold table. Locality is
            # judged BEFORE sublattice folding: folding maps an off-site
            # bond to (0,0,0) between supercell orbitals, which would
            # smuggle unadjudicated off-site content into the correction
            # (measured: a one-orbital +-x bond folded with SubShape=[2,1,1]
            # produced a spurious Fierz tensor of magnitude V). The pre-fold
            # table is filtered to its on-site a != b entries, and only the
            # filtered table is folded.
            has_sub = getattr(self.lattice, "has_sublattice", False)
            if tbl is None:
                if has_sub:
                    tbl = self.param_ham_orig.get(type, {})
                else:
                    tbl = self.param_ham.get(type, {})
            filtered = {}
            for (irvec, orbvec), v in tbl.items():
                if tuple(irvec) == (0, 0, 0) and orbvec[0] != orbvec[1]:
                    filtered[(irvec, orbvec)] = v
            if not filtered:
                return
            if has_sub:
                filtered = self._reshape_interaction(filtered, False)
            filtered = _symmetrised(type, filtered)
            # spin-resolved coefficients DERIVED from the adjudicated S/C
            # table via the channel decomposition W_same = (C - S)/2,
            # W_cross = (C + S)/2 (hwave.solver.vertex_table) -- the values
            # previously hand-coded here per type now have one source
            w_same, w_cross = fierz_coefficients(type)
            spins = {}
            if w_same != 0.0:
                spins[(0, 0, 0, 0)] = w_same
                spins[(1, 1, 1, 1)] = w_same
            if w_cross != 0.0:
                spins[(0, 0, 1, 1)] = w_cross
                spins[(1, 1, 0, 0)] = w_cross
            if not spins:
                return
            for (irvec, orbvec), v in filtered.items():
                if tuple(irvec) != (0, 0, 0):
                    continue
                a, b = orbvec
                if a == b:
                    continue
                for spinvec, w in spins.items():
                    s1, s2, s3, s4 = spinvec
                    # pair-diagonal placement -- bilinear (a,b) coupled to
                    # bilinear (a,b) -- selected by measurement: it makes
                    # the ring's chiq equal the FLEX channel solve to 3e-16
                    # in every adjudicated cell, while the pair-antidiagonal
                    # placement reproduces the old 1e-2 deficiency
                    orb = (s4, a, s3, b, s1, a, s2, b)
                    fierz_r[(*irvec, *orb)] += v * w

        if 'Coulomb' in self.param_ham.keys():
            # The aggregate 'Coulomb' input provides both the intra and inter
            # parts (shared decomposition, see wan90.split_coulomb).
            # Combining it with explicit CoulombIntra/CoulombInter is
            # ambiguous.
            if ('CoulombIntra' in self.param_ham.keys()
                    or 'CoulombInter' in self.param_ham.keys()):
                logger.error(
                    "Coulomb cannot be specified together with "
                    "CoulombIntra or CoulombInter")
                sys.exit(1)

            coulomb_intra, coulomb_inter = wan90.split_coulomb(
                self.param_ham['Coulomb'])

            _append_inter('CoulombIntra', coulomb_intra)
            _append_inter('CoulombInter', coulomb_inter)
            if getattr(self.lattice, "has_sublattice", False):
                _pre_intra, _pre_inter = wan90.split_coulomb(
                    self.param_ham_orig['Coulomb'])
            else:
                _pre_intra, _pre_inter = coulomb_intra, coulomb_inter
            _append_inter_cross('CoulombInter', tbl=_pre_inter)
            _append_onsite_direct('CoulombIntra', tbl=_pre_intra)
            _append_onsite_direct('CoulombInter', tbl=_pre_inter)
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'CoulombIntra' in self.param_ham.keys():
            _append_inter('CoulombIntra')
            _append_onsite_direct('CoulombIntra')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'CoulombInter' in self.param_ham.keys():
            _append_inter('CoulombInter')
            _append_inter_cross('CoulombInter')
            _append_onsite_direct('CoulombInter')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'Hund' in self.param_ham.keys():
            _append_inter('Hund')
            _append_inter_cross('Hund')
            _append_onsite_direct('Hund')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'Ising' in self.param_ham.keys():
            _append_inter('Ising')
            _append_inter_cross('Ising')
            _append_onsite_direct('Ising')
            self._has_interaction = True
            self._has_interaction_coulomb = True

        if 'PairLift' in self.param_ham.keys():
            _append_inter('PairLift')
            _append_onsite_direct('PairLift')
            self._has_interaction = True
            self._has_interaction_exchange = True

        if 'Exchange' in self.param_ham.keys():
            _append_inter('Exchange')
            _append_inter_cross('Exchange')
            _append_onsite_direct('Exchange')
            self._has_interaction = True
            self._has_interaction_exchange = True

        if 'PairHop' in self.param_ham.keys():
            _append_pairhop('PairHop')
            _append_onsite_direct('PairHop', pairhop=True)
            self._has_interaction = True
            self._has_interaction_pairhop = True

        # reshape to W(r)^{bb'aa'}, a,b=(spin,alpha)
        ham_r = ham_r.reshape(nx,ny,nz,*(nd,)*4)

        logger.debug("ham_inter_r shape={}, size={}".format(ham_r.shape, ham_r.size))
        logger.debug("ham_inter_r nonzero count={}".format(ham_r[abs(ham_r) > 1.0e-8].size))

        # Fourier transform W(q)^{bb'aa'}
        # e^{+iqR} (issue #133): the spin-orbital transfer already used
        # the documented sign while this shared FFT used the opposite one,
        # so spin-orbital chiq combined chi0(q) with W(-q) -- invisible for
        # R-symmetric real interactions, wrong for Hermitian-closed complex
        # off-site declarations.
        ham_q = FFT.ifftn(ham_r, axes=(0,1,2)) * nvol

        logger.debug("ham_inter_q shape={}, size={}".format(ham_q.shape, ham_q.size))
        logger.debug("ham_inter_q nonzero count={}".format(ham_q[abs(ham_q) > 1.0e-8].size))

        self.ham_inter_r = ham_r
        self.ham_inter_q = ham_q

        fierz_r = fierz_r.reshape(nx,ny,nz,*(nd,)*4)
        if np.any(fierz_r):
            self.ham_fierz_q = FFT.ifftn(fierz_r, axes=(0,1,2)) * nvol
        else:
            self.ham_fierz_q = None

        # Spinful exchange content (issue #137). The spinful general solve
        # resums chi = [1 - chi0 W]^-1 chi0 with ONE vertex tensor; the
        # physically complete first order is the ANTISYMMETRIZED bare
        # particle-hole vertex Gamma = D + X, where X is the crossed
        # (exchange) wiring of the same interaction. Re-pairing the on-site
        # two-body term
        #     W^{b b' a a'} c^+_a c_a' c^+_b' c_b
        #       = - W^{b b' a a'} c^+_a c_b c^+_b' c_a' + (one-body)
        # shows the crossed wiring is the direct wiring of the coefficient
        # tensor with the b and a' slots swapped and negated:
        #     X[b, b', a, a'] = - D[a', b', a, b]   (on-site block only).
        # Off-site exchange depends on both fermionic momenta and cannot be
        # written as W(q); it stays outside the ring-form resummation (the
        # same limitation the non-spin-orbital ladder has). In the
        # spin-conserving limit X reproduces the adjudicated transverse
        # (ring+ladder) vertex; ED confirmed Gamma = D + X for CoulombIntra
        # (issue #137 reproduction).
        if onsite_r is not None:
            exch_onsite = -np.transpose(
                onsite_r.reshape(*(nd,) * 4), (3, 1, 2, 0))
            self.ham_spinful_exchange = (exch_onsite
                                         if np.any(exch_onsite) else None)
        else:
            self.ham_spinful_exchange = None

class RPA:
    """
    RPA calculation
    """
    @do_profile
    def __init__(self, param_ham, info_log, info_mode):
        logger.debug(">>> RPA.__init__")

        self.param_ham = param_ham
        self.info_log = info_log
        self.param_mod = CaseInsensitiveDict(info_mode.get("param", {}))

        if str(info_mode.get("calc_scheme", "auto")).lower() == "squashed":
            # Removed in 2.0 (issue #144): squashed computed the same
            # susceptibility as reduced at several times the cost, and the
            # spin-off-diagonal slots of its 8-axis output were structurally
            # zero. Only the CONFIGURATION is rejected -- 'squashed' recorded
            # in the in-memory provenance metadata (chi0q_freq_meta) of a
            # reused chi0q is still accepted by the reuse route below (the
            # two share one bubble representation). Checked before
            # Lattice/Interaction construction so that no other validation
            # error can mask this message.
            raise ValueError(
                "calc_scheme='squashed' was removed in H-wave 2.0: it "
                "computed the same susceptibility as calc_scheme='reduced' "
                "at several times the cost, and the extra spin-resolved "
                "slots of its output were structurally zero (issue #144). "
                "Use calc_scheme='reduced'. The physics is identical; the "
                "output layout changes from (l,q,s1,s2,a,s3,s4,b) to "
                "(l,q,a,b) over generalized indices a = s*norb + orb. See "
                "the migration note in the RPA output-file documentation.")

        self.lattice = Lattice(self.param_mod)
        self.ham_info = Interaction(self.lattice, param_ham, info_mode)

        self._set_scheme(info_mode)

        self._init_param()
        self._show_params()

        pass

    def _set_scheme(self, info_mode):
        # handle calc_scheme: must be called after setting up interactions

        self.calc_scheme = info_mode.get("calc_scheme", "auto")

        # calc_type: "ring" (default) or "ring+ladder"
        self.calc_type = info_mode.get("calc_type", "ring")
        if self.calc_type not in ["ring", "ring+ladder"]:
            logger.error("calc_type must be 'ring' or 'ring+ladder', got '{}'".format(self.calc_type))
            sys.exit(1)

        # auto choose
        if self.calc_scheme == "auto":
            if not self.ham_info.has_interaction():
                logger.error("calc_scheme must be specified for chi0q-only mode.")
                sys.exit(1)
            else:
                if self.calc_type == "ring+ladder":
                    # ladder diagrams require general scheme (full rank-4 tensor)
                    self.calc_scheme = "general"
                    logger.info("auto mode for calc_scheme: set to general (ring+ladder)")
                elif any(self.param_ham.get(t)
                         for t in ("Exchange", "PairHop")):
                    # Exchange and PairHop carry NO density-diagonal vertex
                    # content: only the general scheme represents them. The
                    # historical auto choice silently dropped them
                    # entirely (#107).
                    self.calc_scheme = "general"
                    logger.info(
                        "auto mode for calc_scheme: set to general "
                        "(Exchange/PairHop have no density-diagonal vertex)")
                else:
                    # PairLift alone does not need the general scheme: its
                    # particle-hole vertex is exactly zero everywhere.
                    self.calc_scheme = "reduced"
                    logger.info("auto mode for calc_scheme: set to reduced")

        # consistency check. The reduced scheme solves with the
        # density-diagonal part of the interaction, and the adjudicated
        # particle-hole vertex of Exchange and PairHop has NO
        # density-diagonal content (hwave.solver.vertex_table: their
        # entries live on the cross and antidiagonal slot families). The
        # projection therefore does not approximate them -- it drops them
        # entirely, with exactly zero effect on the result. That was two
        # different policies before #107: reduced rejected exchange-type
        # input, FLEX warned of an 'approximation'. One policy now, for
        # both solvers (FLEX inherits this check): reject with a pointer
        # to the scheme that keeps them.
        if self.calc_scheme == "reduced":
            dropped = [t for t in ("Exchange", "PairHop")
                       if self.param_ham.get(t)]
            if dropped:
                raise ValueError(
                    "calc_scheme='{}' solves with the density-diagonal "
                    "part of the interaction, and the particle-hole vertex "
                    "of {} has no density-diagonal content: the interaction "
                    "would have exactly zero effect (not an approximation). "
                    "Use calc_scheme='general', which carries the full "
                    "vertex. (Note: FLEX's general scheme is spin-free "
                    "only; with enable_spin_orbital or a spin-polarized "
                    "setup no current FLEX scheme supports these "
                    "interactions.)".format(
                        self.calc_scheme, ", ".join(dropped)))
            if self.param_ham.get("PairLift"):
                logger.warning(
                    "PairLift's particle-hole vertex is exactly zero "
                    "(adjudicated against exact diagonalization), so it "
                    "has no effect on the susceptibility channels in any "
                    "scheme.")
        if self.calc_type == "ring+ladder" and self.calc_scheme != "general":
            logger.error("calc_type='ring+ladder' requires calc_scheme='general' or 'auto'.")
            sys.exit(1)

        # calc chiq if interaction term exists; otherwise chi0q-only mode
        self.calc_chiq = self.ham_info.has_interaction()

        # chi0q in reduced mode if calc_scheme is reduced
        self.enable_reduced = self.calc_scheme.lower() == "reduced"
        
    def _init_param(self):
        logger.debug(">>> RPA._init_param")

        self.T = self.param_mod.get("T", 0.0)
        self.ene_cutoff = self.param_mod.get("ene_cutoff", 1.0e+2)

        self.nmat = self.param_mod.get("Nmat", 1024)

        # Stay consistent with the Interaction's physical-orbital count
        # (already SO-halved and validated); avoids re-deriving / re-checking.
        self.norb = self.ham_info.norb
        self.norb_orig = self.ham_info.norb_orig
        self.ns = 2  # spin dof
        self.nd = self.norb * self.ns

        # finite-real validation (issue #134 deep review): NaN/inf were
        # accepted and silently produced non-finite susceptibilities; large
        # values cause catastrophic cancellation. Booleans are rejected --
        # float(True) == 1.0 would silently enable the correction.
        # Spinful antisymmetrized vertex (issue #137): default ON; False
        # reproduces the pre-#137 ring-only spinful numbers. Strict bool
        # (the g2_tail/#83 lesson: coercion turns "false" into True).
        _sve = self.param_mod.get("spinful_vertex_exchange", True)
        if not isinstance(_sve, (bool, np.bool_)):
            raise ValueError(
                "[mode.param] spinful_vertex_exchange must be a boolean, "
                "got {!r}".format(_sve))
        self.spinful_vertex_exchange = bool(_sve)

        import numbers
        _ct = self.param_mod.get("coeff_tail", 0.0)
        # type-strict, not coercive (round-7 review): float() would also
        # accept numeric strings (a quoted TOML number), booleans and 0-d
        # arrays. numbers.Real keeps Python and NumPy real scalars;
        # booleans are excluded explicitly -- float(True) == 1.0 would
        # silently enable the correction.
        if (isinstance(_ct, (bool, np.bool_))
                or not isinstance(_ct, numbers.Real)):
            raise ValueError(
                "[mode.param] coeff_tail must be a real number, got "
                "{!r}".format(_ct))
        _ct = float(_ct)
        if not np.isfinite(_ct):
            raise ValueError(
                "[mode.param] coeff_tail must be finite, got {}".format(_ct))
        self.coeff_tail = _ct
        self.ext = self.param_mod.get("coeff_extern", 0.0)

        # GPU (CuPy) execution and CPU spatial-FFT parallelism. use_gpu is the
        # user request; consumers resolve it via backend.get_backend (warn and
        # fall back to numpy when CuPy / a CUDA device is unavailable).
        # fft_workers defaults to 1 (the serial numpy path, bit-compatible
        # with previous releases); scipy-parallel FFTs are opt-in.
        self.use_gpu = _bk.as_bool(self.param_mod.get("gpu", False))
        # Strict GPU mode: when gpu_required=true, get_backend raises instead of
        # falling back to the (much slower) CPU path if CuPy/CUDA is unusable,
        # so a large scheduler job fails fast. Default false keeps the existing
        # warn-and-fall-back behavior. Inherited by the FLEX solver (subclass).
        self.gpu_required = _bk.as_bool(self.param_mod.get("gpu_required",
                                                           False))
        self.fft_workers = int(self.param_mod.get("fft_workers", 1))

        # exclusive options: mu and Ncond/filling
        have_mu = "mu" in self.param_mod.keys()
        have_Ncond = "Ncond" in self.param_mod.keys() or "Nelec" in self.param_mod.keys()
        have_filling = "filling" in self.param_mod.keys()

        if have_mu and (have_Ncond or have_filling):
            # conflicting parameters
            logger.error("both mu and Ncond or filling are specified")
            sys.exit(1)
        elif not (have_mu or have_Ncond or have_filling):
            # none specified
            logger.error("none of mu, Ncond, nor filling is specified")
            sys.exit(1)
        elif have_Ncond and have_filling:
            # both Ncond or filling
            logger.error("both Ncond and filling are specified")
            sys.exit(1)

        self.Nstate = self.lattice.nvol * self.nd

        if have_Ncond or have_filling:
            self.calc_mu = True
            if "Ncond" in self.param_mod:
                self.Ncond = self.param_mod["Ncond"]
                self.filling = 1.0 * self.Ncond / self.Nstate
            elif "Nelec" in self.param_mod:
                self.Ncond = self.param_mod["Nelec"]
                self.filling = 1.0 * self.Ncond / self.Nstate
            elif "filling" in self.param_mod:
                self.filling = self.param_mod['filling']
                self.Ncond = self.Nstate * self.filling
                # force Ncond to be integer when T=0
                if self.T == 0.0:
                    round_mode = self.param_mod.get("Ncond_round_mode", "strict")
                    self.Ncond = self._round_to_int(self.Ncond, round_mode)

        if have_mu:
            self.calc_mu = False
            self.mu_value = self.param_mod.get("mu", None)

        # range of matsubara frequency in matrix calculation and output
        self.freq_range = self.param_mod.get("matsubara_frequency", "all")
        self.freq_index = self._find_index_range(self.freq_range)
        logger.debug("freq_index = {}".format(self.freq_index))

        # check parameters
        err = 0
        # Finite-temperature Matsubara formalism: T must be strictly positive
        # (beta = 1/T is used directly).
        if self.T <= 0.0:
            logger.error("T must be greater than zero: T={}".format(self.T))
            err += 1
        # The chemical-potential search needs a partially-filled band; full
        # (Ncond == Nstate) or empty filling leaves mu unbracketed.
        if self.calc_mu and not (0 < self.Ncond < self.Nstate):
            logger.error("Ncond must satisfy 0 < Ncond < Nstate ({}): Ncond={}".format(
                self.Nstate, self.Ncond))
            err += 1
        if self.nmat <= 0:
            logger.error("Nmat must be greater than zero: Nmat={}".format(self.nmat))
            err += 1
        # Fermionic Matsubara grid iomega = (2n+1-Nmat)*pi/beta is symmetric and
        # never zero only when Nmat is even; an odd Nmat injects omega=0.
        if self.nmat % 2 != 0:
            logger.error("Nmat must be even: Nmat={}".format(self.nmat))
            err += 1
        if err > 0:
            sys.exit(1)

        self._init_transverse_bond_config()

        pass

    def _init_transverse_bond_config(self):
        """Parse ``[mode.param] transverse_bond_channels`` and its
        companion options (Phase W experimental gate, spec
        2026-08-15-bond-transverse-design.md "Production surface").

        Mirrors sc.py's ``_read_bond_config`` for the (eliashberg)
        longitudinal bond gate: the switch defaults to ``False``, is
        parsed ONLY under ``calc_type='ring+ladder'`` (the transverse
        ladder channel itself is unreachable otherwise), and every
        companion option set while the switch is stale (calc_type='ring',
        or transverse_bond_channels=false) is IGNORED WITH A WARNING
        rather than silently doing nothing or failing an otherwise valid
        run. ``self.param_mod`` is a ``CaseInsensitiveDict``, so every
        lookup here is case-robust by construction.
        """
        import numbers

        _tbc_keys = ("transverse_bond_channels",
                     "transverse_bond_max_shells",
                     "transverse_bond_memory_cap_gb")
        present = [k for k in _tbc_keys if k in self.param_mod]

        if self.calc_type != "ring+ladder":
            if present:
                logger.warning(
                    "[mode.param] %s set but calc_type='%s'; the "
                    "bond-resolved transverse gate only applies to "
                    "calc_type='ring+ladder' and %s ignored here.",
                    ", ".join(present), self.calc_type,
                    "is" if len(present) == 1 else "are")
            self.transverse_bond_channels = False
            self.transverse_bond_max_shells = None
            self.transverse_bond_memory_cap_gb = \
                TRANSVERSE_BOND_MEMORY_CAP_GB_DEFAULT
            return

        _flag = self.param_mod.get("transverse_bond_channels", False)
        if not isinstance(_flag, (bool, np.bool_)):
            raise ValueError(
                "[mode.param] transverse_bond_channels must be a boolean, "
                "got {!r}".format(_flag))
        self.transverse_bond_channels = bool(_flag)

        if not self.transverse_bond_channels:
            stale = [k for k in _tbc_keys
                     if k != "transverse_bond_channels" and k in self.param_mod]
            if stale:
                logger.warning(
                    "[mode.param] %s set but transverse_bond_channels="
                    "false; these options only apply to "
                    "transverse_bond_channels=true and are ignored here.",
                    ", ".join(stale))
            self.transverse_bond_max_shells = None
            self.transverse_bond_memory_cap_gb = \
                TRANSVERSE_BOND_MEMORY_CAP_GB_DEFAULT
            return

        # transverse_bond_max_shells: int >= 0, or absent/None (keep every
        # declared off-site shell) -- same semantics
        # resolve_transverse_topology documents for its own max_shells.
        _max_shells = self.param_mod.get("transverse_bond_max_shells", None)
        if _max_shells is not None:
            if (isinstance(_max_shells, (bool, np.bool_))
                    or not isinstance(_max_shells, numbers.Integral)):
                raise ValueError(
                    "[mode.param] transverse_bond_max_shells must be a "
                    "non-negative integer, got {!r}".format(_max_shells))
            _max_shells = int(_max_shells)
            if _max_shells < 0:
                raise ValueError(
                    "[mode.param] transverse_bond_max_shells must be >= 0 "
                    "(shell 0 = the on-site Delta r = 0 point), got "
                    "{}".format(_max_shells))
        self.transverse_bond_max_shells = _max_shells

        _cap_gb = self.param_mod.get("transverse_bond_memory_cap_gb",
                                      TRANSVERSE_BOND_MEMORY_CAP_GB_DEFAULT)
        if (isinstance(_cap_gb, (bool, np.bool_))
                or not isinstance(_cap_gb, numbers.Real)):
            raise ValueError(
                "[mode.param] transverse_bond_memory_cap_gb must be a "
                "real number, got {!r}".format(_cap_gb))
        _cap_gb = float(_cap_gb)
        if not np.isfinite(_cap_gb) or _cap_gb <= 0.0:
            raise ValueError(
                "[mode.param] transverse_bond_memory_cap_gb must be a "
                "positive finite number, got {}".format(_cap_gb))
        self.transverse_bond_memory_cap_gb = _cap_gb

    def _round_to_int(self, val, mode):
        import math
        mode = mode.lower()  # case-insensitive
        if   mode == "as-is":
            ret = val  # not rounding to int
        elif mode == "round":
            ret = round(val)
        elif mode == "round-up":
            ret = math.ceil(val)
        elif mode == "round-down":
            ret = math.floor(val)
        elif mode == "round-to-zero":
            ret = int(val)
        elif mode == "round-off":
            nn = math.floor(val)
            rr = val - nn
            ret = nn if rr < 0.5 else nn+1
        elif mode == "strict":
            if val != round(val):
                logger.error("value not integer")
                sys.exit(1)
            ret = round(val)
        elif mode == "exact":  # "round" with warning
            if val != round(val):
                logger.warning("value not integer")
            ret = round(val)
        else:
            logger.error("round_to_int: unknown mode {}".format(mode))
            sys.exit(1)
        return ret

    def _show_params(self):
        logger.debug(">>> RPA._show_params")
        logger.info("RPA parameters:")
        logger.info("    norbit          = {}".format(self.norb))
        logger.info("    nspin           = {}".format(self.ns))
        logger.info("    nd              = {}".format(self.nd))
        logger.info("    Nmat            = {}".format(self.nmat))
        logger.info("    coeff_tail      = {}".format(self.coeff_tail))

        if self.calc_mu == True:
            logger.info("    Ncond           = {}".format(self.Ncond))
            logger.info("    filling         = {:.3f}".format(self.filling))
        else:
            logger.info("    mu              = {}".format(self.mu_value))

        logger.info("    T               = {}".format(self.T))
        logger.info("    E_cutoff        = {:e}".format(self.ene_cutoff))
        logger.info("    coeff_extern    = {}".format(self.ext))

        # logger.info("    RndSeed         = {}".format(self.param_mod["RndSeed"]))
        # logger.info("    strict_hermite  = {}".format(self.strict_hermite))
        # logger.info("    hermite_tol     = {}".format(self.hermite_tolerance))
        logger.info("    freq_range      = {}".format(self.freq_range))
        logger.info("    calc_chiq       = {}".format(self.calc_chiq))
        logger.info("    spin_orbital    = {}".format(self.ham_info.enable_spin_orbital))
        logger.info("    calc_scheme     = {}".format(self.calc_scheme))
        logger.info("    calc_type       = {}".format(self.calc_type))
        logger.info("    transverse_bond_channels = {}".format(
            self.transverse_bond_channels))
        pass

    @do_profile
    def solve(self, green_info, path_to_output):
        """Solve the RPA equation, restoring host-backed public state.

        Thin wrapper around :meth:`_solve_impl` that guarantees the solver's
        public array attributes (``H0_eigenvalue``/``H0_eigenvector`` and the
        stored ``green0``/``green0_tail``) are NumPy-backed after the call --
        on normal completion AND after a GPU-path exception. Under GPU
        execution ``_solve_impl`` converts these to CuPy in place; without the
        ``finally`` a mid-solve error would leave a reused or inspected solver
        object holding device arrays (issue #63).
        """
        try:
            return self._solve_impl(green_info, path_to_output)
        finally:
            _bk.restore_host_attrs(
                self, ("H0_eigenvalue", "H0_eigenvector",
                       "green0", "green0_tail"))

    def _solve_impl(self, green_info, path_to_output):
        """Solve the RPA equation to calculate susceptibility.

        This is the main method that performs RPA calculations. It either calculates
        or loads chi0q, transforms interaction Hamiltonians based on spin state,
        and solves the RPA equation.

        Parameters
        ----------
        green_info : dict
            Dictionary containing Green's function information and calculation parameters.
            May include pre-calculated chi0q.
        path_to_output : str
            Path to directory where output files will be saved.

        Notes
        -----
        The calculation flow is:
        1. Calculate/load chi0q
        2. Transform interactions based on spin state (spin-free/spinful/spin-diag)
        3. Transform tensors based on calculation scheme (reduced/general)
        4. Solve RPA equation to get chiq
        """
        logger.info("Start RPA calculations")

        # A reused green_info must not carry results from a previous solve:
        # a stale chiq / chiq_pm would survive into save_results and be
        # labelled as part of this result. Dropped at ENTRY -- before any
        # validation and regardless of calc_chiq -- so a failed validation
        # or a chi0q-only run cannot leave them behind either (adversarial
        # review, round 3).
        green_info.pop("chiq", None)
        green_info.pop("chiq_pm", None)
        green_info.pop("chiq_pm_bond_static", None)
        green_info.pop("chiq_pm_static", None)

        beta = 1.0/self.T

        # GPU (CuPy) execution: resolve the backend once. The heavy work --
        # the bare Green's function, the chi0q FFT pair bubble, the spin
        # inflation einsums, and the batched chiq solve -- all dispatch on
        # their input arrays, so moving the inputs to the device here runs the
        # whole solve on the GPU; outputs are stored back as host arrays.
        xp, gpu_active = _bk.get_backend(self.use_gpu, logger=logger,
                                         required=self.gpu_required)

        if "chi0q" in green_info and green_info["chi0q"] is not None:
            # use chi0q input; a green_info stored by a previous solve
            # arrives here, so establish spin_mode from the shape exactly
            # as the file-based chi0q_init route does (issue #109)
            chi0q = green_info["chi0q"]
            if not gpu_active:
                # normalize a device array to the selected (host) backend
                # FIRST: real device arrays forbid the implicit conversion
                # a NumPy-based dtype probe would attempt
                chi0q = _bk.to_host(chi0q)
                green_info["chi0q"] = chi0q
            if not np.issubdtype(chi0q.dtype, np.number):
                raise ValueError(
                    "chi0q from green_info has non-numeric dtype {}".format(
                        chi0q.dtype))
            self._validate_chi0q_shape(chi0q, source="green_info")
            # spin-diag arrays carry the spin-block axis first; the
            # frequency axis is shape[1] there (shape[0] == 2 == nblock,
            # which used to be misreported as a partial frequency range)
            nfreq = (chi0q.shape[1] if self.spin_mode == "spin-diag"
                     else chi0q.shape[0])
            # Frequency provenance for save_results: inherit the axis of
            # the run that PRODUCED this chi0q. Priority: metadata stored
            # by a previous in-memory solve; else metadata already set by
            # the chi0q_init file reader on this instance; else -- an
            # untagged partial array -- fall back to the ambiguous 0..n-1
            # labeling rather than fabricating a full-axis claim.
            mem_meta = green_info.get("chi0q_freq_meta")
            if mem_meta is not None:
                # Validate the provenance against the array and this solver
                # (shared validator; round-3/5 reviews): a caller can
                # replace chi0q while leaving the previous solve's metadata
                # behind, and trusting it silently would mislabel saved
                # output or solve a semantically different tensor.
                norm = _validate_chi0q_provenance(
                    mem_meta, nfreq, source="green_info")
                # the validator has already guaranteed a Mapping, so read
                # the fingerprint unconditionally -- a dict-only check let
                # a UserDict's stale fingerprint pass unverified
                fp = mem_meta.get("fingerprint")
                if fp is not None and fp != _chi0q_fingerprint(chi0q):
                    raise ValueError(
                        "chi0q_freq_meta does not belong to the supplied "
                        "chi0q (content fingerprint mismatch): the array "
                        "was replaced while the previous solve's metadata "
                        "was kept. The metadata is stale -- recompute or "
                        "drop the key.")
                m_norb = mem_meta.get("norb")
                if m_norb is not None and int(m_norb) != int(self.norb):
                    raise ValueError(
                        "chi0q was produced with norb = {} but this solver "
                        "has norb = {}; the bubble does not describe this "
                        "system.".format(m_norb, self.norb))
                m_scheme = mem_meta.get("calc_scheme")
                if m_scheme is not None:
                    # reduced and squashed share one bubble representation
                    # (the reduced 'kab' layout; measured: a reduced-produced
                    # bubble solved under squashed matches an internal
                    # squashed recomputation exactly), so only the
                    # REPRESENTATION class must agree
                    rep = {"reduced": "reduced", "squashed": "reduced"}
                    if (rep.get(m_scheme, m_scheme)
                            != rep.get(self.calc_scheme, self.calc_scheme)):
                        raise ValueError(
                            "chi0q was produced under calc_scheme = '{}' "
                            "whose bubble representation is incompatible "
                            "with this solver's '{}'.".format(
                                m_scheme, self.calc_scheme))
                m_spin = mem_meta.get("spin_mode")
                if m_spin is not None and m_spin != self.spin_mode:
                    # the producer knows its spin structure; shape inference
                    # cannot separate spin-free norb-orbital from spinful
                    # norb/2-orbital input, so a declared mismatch is fatal
                    raise ValueError(
                        "chi0q was produced in spin mode '{}' but its shape "
                        "reads as '{}' for this solver's configuration; "
                        "refusing to solve a semantically different "
                        "tensor.".format(m_spin, self.spin_mode))
                # endpoint gate on the SAME normalized form the file
                # route checks (round-7 review: this route accepted a
                # nonzero-tail bubble without or with an unknown marker)
                _enforce_tail_endpoint(norm, "green_info chi0q_freq_meta")
                self._chi0q_init_meta = dict(norm)
            elif getattr(self, "_chi0q_init_meta", None) is not None:
                # metadata set by the chi0q_init file reader on this
                # instance: verify it still describes THIS array before
                # trusting it (the caller may have replaced info["chi0q"]
                # between read_init and solve)
                file_fp = self._chi0q_init_meta.get("fingerprint")
                if (file_fp is not None
                        and file_fp != _chi0q_fingerprint(chi0q)):
                    raise ValueError(
                        "the chi0q_init metadata on this solver does not "
                        "belong to the supplied chi0q (content fingerprint "
                        "mismatch): the array was replaced after the file "
                        "was read. Recompute or re-read the file.")
            else:
                # untagged external tensor -- full-frequency included
                # (round-7 review: a full-Nmat external array previously
                # fell through with no metadata, and save_results then
                # stamped the CURRENT run's coeff_tail and endpoint
                # marker onto data of unknown provenance)
                self._chi0q_init_meta = {
                    "freq_index": None, "nmat": None, "coeff_tail": None,
                    "tail_endpoint": None,
                }
            # Write the normalized provenance back so a CHAIN of reuses --
            # including a bubble that entered through the chi0q_init file
            # route -- keeps its producing axis (round-3 review).
            eff = getattr(self, "_chi0q_init_meta", None)
            if eff is not None:
                fi_eff = eff.get("freq_index")
                green_info["chi0q_freq_meta"] = {
                    "freq_index": (None if fi_eff is None
                                   else list(np.asarray(fi_eff))),
                    "nmat": eff.get("nmat"),
                    "coeff_tail": eff.get("coeff_tail"),
                    "tail_endpoint": eff.get("tail_endpoint"),
                    "spin_mode": self.spin_mode,
                    "norb": int(self.norb),
                    "calc_scheme": self.calc_scheme,
                    "fingerprint": _chi0q_fingerprint(chi0q),
                }
            # must_fix 6 guard state: the spin-diag transverse channel needs
            # Green functions bound to THIS bubble; an externally supplied
            # chi0q has none (a same-instance green0 may belong to an older
            # bubble)
            self._chi0q_external = True
            if (getattr(self, "calc_type", "ring") == "ring+ladder"
                    and self.spin_mode == "spin-diag"):
                # fail HERE, before the longitudinal solve populates chiq:
                # rejecting inside the transverse build would leave a fresh
                # chiq in green_info from a run that then failed
                raise ValueError(
                    "spin-diag transverse (ladder) channel cannot be "
                    "combined with an externally supplied chi0q: the "
                    "channel needs the Green's functions that produced "
                    "this exact bubble. Recompute chi0q internally.")
            if nfreq != self.nmat:
                logger.info("partial range in matsubara frequency: {} in {}".format(nfreq, self.nmat))
                #self.nmat = chi0q.shape[0]
            if gpu_active:
                # VRAM preflight for the externally-supplied chi0q path: the
                # transfer below plus the same-sized chiq solve workspace.
                # Advisory only (CuPy raises OutOfMemoryError on the actual
                # allocation).
                _bk.warn_if_device_memory_short(
                    2 * chi0q.nbytes, logger,
                    label="the RPA chiq solve (supplied chi0q)")
                chi0q = xp.asarray(chi0q)
        else:
            self._calc_epsilon_k(green_info)

            if self.calc_mu:
                if self.spin_mode == "spin-free":
                    Ncond = self.Ncond/2
                else:
                    Ncond = self.Ncond
                dist, mu = self._find_mu(Ncond, self.T)
            else:
                mu = self.mu_value

            if gpu_active:
                logger.info("RPA: GPU backend active (CuPy); moving H0 "
                            "eigenpairs to the device.")
                # VRAM preflight: the largest resident device tensor is the
                # inflated chi0q / chiq, ~ Nmat*Nvol*nd^4 complex128 in the full
                # (rank-4 orbital) channel; the chiq solve holds a same-sized
                # workspace. This nd^4 figure is an upper bound for the reduced
                # (rank-2) scheme and a rough order-of-magnitude
                # estimate otherwise -- advisory only (CuPy raises
                # OutOfMemoryError on the actual allocation).
                # H0_eigenvector shape = (nblock, Nvol, nd, nd).
                nd0 = self.H0_eigenvector.shape[-1]
                chi_bytes = self.nmat * self.lattice.nvol * (nd0 ** 4) * 16
                _bk.warn_if_device_memory_short(
                    2 * chi_bytes, logger, label="the RPA chi0q/chiq solve")
                self.H0_eigenvalue = xp.asarray(self.H0_eigenvalue)
                self.H0_eigenvector = xp.asarray(self.H0_eigenvector)

            green0, green0_tail = self._calc_green(beta, mu)
            #XXX
            self.green0 = green0
            self.green0_tail = green0_tail
            self._chi0q_external = False
            # a previous chi0q_init on this instance must not relabel the
            # axis of a bubble THIS run computes
            self._chi0q_init_meta = None

            chi0q = self._calc_chi0q(green0, green0_tail, beta)

            # filter by matsubara freq range
            if len(self.freq_index) < self.nmat:
                chi0q = chi0q[:,self.freq_index]
                logger.info("filter range in matsubara frequency: {} in {}".format(chi0q.shape[0], self.nmat))

            # nblock: defense in depth behind the kernel's structural
            # check -- under python -O the former asserts vanished and a
            # wrong block count was silently truncated (round-5 review)
            if self.spin_mode in [ "spin-free", "spinful" ]:
                if chi0q.shape[0] != 1:
                    raise ValueError(
                        "chi0q block count {} does not match spin_mode "
                        "'{}' (expected 1)".format(
                            chi0q.shape[0], self.spin_mode))
                chi0q = chi0q[0]
            else:
                if chi0q.shape[0] != 2:
                    raise ValueError(
                        "chi0q block count {} does not match spin_mode "
                        "'{}' (expected 2)".format(
                            chi0q.shape[0], self.spin_mode))

            green_info["chi0q"] = _bk.to_host(chi0q)
            # Record the producing run's frequency metadata alongside: a
            # later solver reusing this green_info (issue #109) must label
            # its outputs with THIS axis, exactly as the file route does via
            # _chi0q_init_meta -- stamping the reusing run's own
            # freq_index/nmat mislabeled a partial-range bubble as a full
            # axis (measured: 9 frequencies saved with freq_index 0..31).
            green_info["chi0q_freq_meta"] = {
                "freq_index": list(self.freq_index),
                "nmat": int(self.nmat),
                "coeff_tail": float(getattr(self, "coeff_tail", 0.0)),
                "tail_endpoint": TAIL_ENDPOINT_CONVENTION,
                # producer identity: shape-only spin detection cannot
                # distinguish a spin-free two-orbital bubble from a spinful
                # one-orbital one, so the consumer validates against these
                "spin_mode": self.spin_mode,
                "norb": int(self.norb),
                "calc_scheme": self.calc_scheme,
                # content binding: a same-shaped replacement of the array
                # must not inherit this metadata (round-5 review)
                "fingerprint": _chi0q_fingerprint(green_info["chi0q"]),
            }

        if self.calc_chiq:
            # ham_inter_q is built on the host at init; mirror it to chi0q's
            # backend so the inflation einsums and the solve stay on one device.
            ham_inter_q = self.ham_info.ham_inter_q
            # The longitudinal channels solve with the Fierz-corrected
            # tensor; ham_inter_q itself stays as the transverse assembly's
            # input (see _make_ham_inter). The correction's (a,b,a,b) slots
            # sit off the 'kaabb' diagonal, so the reduced reads
            # below are unaffected by construction.
            fierz_q = getattr(self.ham_info, "ham_fierz_q", None)
            if gpu_active:
                ham_inter_q = xp.asarray(ham_inter_q)

            def _fierz_long():
                # assembled lazily: the spinful antisymmetrized branch
                # replaces the Fierz-corrected tensor entirely, and
                # building + device-transferring it there would allocate
                # a large temporary that is immediately discarded
                # (round-1 review). The result is handed to the ring
                # solve in the BUBBLE's pair convention (issue #139).
                if fierz_q is None:
                    return _to_bubble_pair_convention(ham_inter_q)
                fz = xp.asarray(fierz_q) if gpu_active else fierz_q
                return _to_bubble_pair_convention(ham_inter_q + fz)

            if self.spin_mode == "spinful":
                chi0q_orig = chi0q
                # The transverse (ladder) assembly re-pairs the tensor
                # itself, so it must receive the HAMILTONIAN convention,
                # NOT the bubble one the ring solve needs (issue #139):
                # converting here as well conjugated a complex PairHop's
                # ladder vertex (exact diagonalization on a three-site
                # ring: relative error 1.2 against 2e-7 unconverted).
                ham_orig = ham_inter_q
                ham_long = None
                # Antisymmetrized vertex (issue #137): resum with
                # Gamma = D + crossed(D)|on-site so the spin-flip pair
                # slots are corrected (ring-only left them at the bare
                # bubble). The longitudinal Fierz tensor is the DENSITY-
                # slot projection of the same crossing (built for the
                # non-spin-orbital solver, which has no spin-flip slots),
                # so it must NOT be added on top -- that double-counts
                # every on-site cross-orbital slot. Opt-out via
                # [mode.param] spinful_vertex_exchange = false reproduces
                # the pre-#137 ring-only numbers (ham_inter + fierz).
                exch = getattr(self.ham_info, "ham_spinful_exchange", None)
                if exch is not None and self.spinful_vertex_exchange:
                    exch = xp.asarray(exch) if gpu_active else exch
                    ham_long = _to_bubble_pair_convention(
                        ham_inter_q + exch[None, ...])
                else:
                    ham_long = _fierz_long()

                if self.calc_scheme == "reduced":
                    # Treat combined spin-orbital indices as general orbitals.
                    # Block structure is exploited by _find_block_diagonal
                    # inside _solve_rpa.
                    nvol = self.lattice.nvol
                    nd = self.nd
                    ham = project_density_pairs(ham_long, nvol, nd, xp)
                else:
                    ham = ham_long

            elif self.spin_mode == "spin-diag":
                chi0q_orig = chi0q
                ham_orig = ham_inter_q
                ham_long = _fierz_long()

                if self.calc_scheme == "reduced":
                    nblock,nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = xp.identity(2)
                    chi0q = xp.einsum('glkab,gh->lkgahb',
                                      chi0q_orig,
                                      spin_tensor).reshape(nfreq,nvol,nd,nd)

                    ham = project_density_pairs(ham_long, nvol, nd, xp)

                else:
                    nblock,nfreq,nvol,norb1,norb2,norb3,norb4 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = xp.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = xp.einsum('glkabcd,gtuv->lkgatbucvd',
                                      chi0q_orig,
                                      spin_tensor).reshape(nfreq,nvol,nd,nd,nd,nd)
                    ham = ham_long

            elif self.spin_mode == "spin-free":
                # introduce spin degree of freedom
                chi0q_orig = chi0q
                ham_orig = ham_inter_q
                ham_long = _fierz_long()

                if self.calc_scheme == "reduced":
                    # alpha=alpha', beta=beta' case
                    nfreq,nvol,norb1,norb2 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = xp.identity(ns)
                    chi0q = xp.einsum('lkab,st->lksatb',
                                      chi0q_orig.reshape(nfreq,nvol,norb,norb),
                                      spin_tensor).reshape(nfreq,nvol,nd,nd)

                    # same spin-major extraction as the combined-index
                    # pair diagonal (verified bitwise); shared projection
                    ham = project_density_pairs(ham_long, nvol, nd, xp)

                else:
                    # general nd**4 case
                    nfreq,nvol,norb1,norb2,norb3,norb4 = chi0q_orig.shape

                    norb = self.norb
                    ns = self.ns
                    nd = norb * ns

                    spin_tensor = xp.zeros((2,2,2,2), dtype=np.int32)
                    spin_tensor[0,0,0,0] = 1
                    spin_tensor[1,1,1,1] = 1

                    chi0q = xp.einsum('lkabcd,stuv->lksatbucvd',
                                      chi0q_orig.reshape(nfreq,nvol,norb,norb,norb,norb),
                                      spin_tensor).reshape(nfreq,nvol,nd,nd,nd,nd)
                    ham = ham_long

            # For ring+ladder, validate the transverse channel BEFORE the
            # longitudinal solve: an unrepresentable input should fail here,
            # not after the expensive solve has already populated chiq.
            # transverse_bond_channels=true (Phase W experimental gate)
            # represents off-site cross-spin/spin-flip content that
            # _check_transverse_representable would otherwise reject, so
            # the representability check is bypassed ONLY on that path --
            # with the gate off this branch is unchanged, so the bypass is
            # unreachable there.
            if self.calc_type == "ring+ladder":
                if self.transverse_bond_channels:
                    self._transverse_bond_topo = \
                        self._validate_transverse_bond_prereqs()
                    # Resource preflight immediately AFTER the (cheap)
                    # prereq validation but BEFORE ANY expensive solve --
                    # including the longitudinal (ring) solve just below,
                    # not only the bond-resolved solve -- so a cap
                    # rejection fires before either expensive step runs.
                    self._transverse_bond_resource_preflight(
                        self._transverse_bond_topo)
                else:
                    self._check_transverse_representable(ham_orig)

            # solve longitudinal (ring) RPA
            sol = self._solve_rpa(chi0q, ham)

            # adhoc store (as a host array: the writers are numpy)
            green_info["chiq"] = _bk.to_host(sol)

            # Solve transverse (ladder) RPA if requested
            if self.calc_type == "ring+ladder":
                if self.transverse_bond_channels:
                    chi_pm_bond_static, chiq_pm_static = \
                        self._run_transverse_bond_pipeline(
                            self.green0, self.green0_tail, beta,
                            self._transverse_bond_topo)
                    green_info["chiq_pm_bond_static"] = \
                        _bk.to_host(chi_pm_bond_static)
                    green_info["chiq_pm_static"] = \
                        _bk.to_host(chiq_pm_static)
                    # Gated-run output ownership (spec): the plain ladder
                    # is NOT additionally executed, so green_info["chiq_pm"]
                    # (the legacy key) is intentionally never populated here.
                elif self.spin_mode == "spinful":
                    # Phase S (issue #110): the spinful plain-ladder path
                    # no longer slices chi0 and re-solves it (that dropped
                    # every cross-spin/spin-mixing vertex contribution the
                    # RPA equation never got to resum). Instead it slices
                    # the ALREADY-DRESSED longitudinal `sol` (computed
                    # just above, reused verbatim -- no second solve) --
                    # slice-AFTER-solve, not solve-after-slice. See
                    # `_extract_transverse_from_dressed`'s docstring for
                    # the equation and its Gate-S0 ED adjudication.
                    green_info["chiq_pm"] = _bk.to_host(
                        self._extract_transverse_from_dressed(sol))
                else:
                    chi0q_pm, ham_pm = self._build_transverse_channel(
                        chi0q_orig, ham_orig)
                    sol_pm = self._solve_rpa(chi0q_pm, ham_pm)
                    green_info["chiq_pm"] = _bk.to_host(sol_pm)

        # Restore the solver's public attributes to host arrays so the
        # post-solve object state is backend-independent for downstream
        # consumers (save_results, tests).
        if gpu_active:
            if getattr(self, "green0", None) is not None:
                self.green0 = _bk.to_host(self.green0)
                self.green0_tail = _bk.to_host(self.green0_tail)
            if getattr(self, "H0_eigenvalue", None) is not None:
                self.H0_eigenvalue = _bk.to_host(self.H0_eigenvalue)
                self.H0_eigenvector = _bk.to_host(self.H0_eigenvector)

        logger.info("End RPA calculations")
        pass

    @do_profile
    def save_results(self, info_outputfile, green_info):
        """Save calculation results to files.

        Parameters
        ----------
        info_outputfile : dict
            Dictionary containing output file configuration.
        green_info : dict
            Dictionary containing calculation information and results.

        Notes
        -----
        Saves the following information:
        - Calculated correlation functions (chi0q/chiq)
        - Matsubara frequency indices
        - Wave vector unit vectors
        - Wave vector indices
        """
        logger.info("Save RPA results")
        path_to_output = info_outputfile["path_to_output"]

        self._init_wavevec()

        # Frequency metadata (freq_index + nmat, the full grid size of the
        # producing run) lets consumers locate the zero bosonic frequency
        # (original index nmat//2).  A chi0q loaded via chi0q_init passes
        # through solve() untouched AND chiq is computed from it, so both
        # outputs inherit the input file's frequency axis: re-save them with
        # the metadata of the file that produced the chi0q -- stamping the
        # CURRENT run's values would mislabel the axis.
        def _freq_meta_kwargs(arr):
            # the frequency axis of a spin-diag bare bubble is axis 1
            # (axis 0 is the two-spin block); every other saved array
            # leads with frequency
            if arr.ndim in (5, 7) and arr.shape[0] == 2:
                nfreq_arr = arr.shape[1]
            else:
                nfreq_arr = arr.shape[0]
            init_meta = getattr(self, "_chi0q_init_meta", None)
            if init_meta is None:
                # coeff_tail provenance (issue #80): the tail correction
                # changes chi0q at O(1), so record which value produced the
                # file. Lets loaders warn when a chi0q_mode="calc" config
                # would recompute with a different setting. (getattr: tests
                # drive save_results on __new__-built stubs without __init__.)
                return {"freq_index": self.freq_index, "nmat": self.nmat,
                        "coeff_tail": getattr(self, "coeff_tail", 0.0),
                        # endpoint convention of THIS build (issue #134)
                        "tail_endpoint": TAIL_ENDPOINT_CONVENTION}
            kwargs = {}
            if init_meta["freq_index"] is not None:
                kwargs["freq_index"] = init_meta["freq_index"]
            else:
                # The documented format guarantees a freq_index key.  For a
                # legacy chi0q_init file without one, describe the stored
                # axis (0..n-1) without fabricating an nmat claim; a
                # downstream loader then treats the file explicitly as
                # ambiguous unless its config Nmat matches.
                kwargs["freq_index"] = np.arange(nfreq_arr)
            if init_meta["nmat"] is not None:
                kwargs["nmat"] = init_meta["nmat"]
            # pass-through: re-save the INPUT file's coeff_tail (stamping the
            # current run's value would mislabel data this run did not
            # compute); omit the key when the input did not record one.
            if init_meta.get("coeff_tail") is not None:
                kwargs["coeff_tail"] = init_meta["coeff_tail"]
            # same pass-through for the endpoint marker: never fabricate
            # one for data this run did not compute
            if init_meta.get("tail_endpoint") is not None:
                kwargs["tail_endpoint"] = init_meta["tail_endpoint"]
            return kwargs

        if "chiq" in info_outputfile.keys():
            if self.calc_chiq == True:
                file_name = os.path.join(path_to_output, info_outputfile["chiq"])
                save_kwargs = dict(
                    chiq = green_info["chiq"],
                    wavevector_unit = self.kvec,
                    wavevector_index = self.wavenum_table,
                    # RPA orders spin-orbital axes as spin-block (spin*norb+orb),
                    # unlike UHFk's interleaved (2*orb+spin) output; record it so
                    # consumers do not silently mix the two conventions.
                    index_convention = "spin_block",
                    momentum_convention = MOMENTUM_CONVENTION,
                    **_freq_meta_kwargs(green_info["chiq"]),
                )
                # transverse channel chi_+-(q), present for calc_type ring+ladder
                if green_info.get("chiq_pm") is not None:
                    save_kwargs["chiq_pm"] = green_info["chiq_pm"]
                # Bond-resolved transverse channel (Phase W experimental
                # gate, transverse_bond_channels=true): gate-owned output,
                # mutually exclusive with the legacy chiq_pm key above --
                # solve() only ever populates one of the two. Schema keys
                # verbatim per spec ("Production surface").
                if green_info.get("chiq_pm_bond_static") is not None:
                    topo = self._transverse_bond_topo
                    save_kwargs["chiq_pm_bond_static"] = \
                        green_info["chiq_pm_bond_static"]
                    save_kwargs["chiq_pm_static"] = \
                        green_info["chiq_pm_static"]
                    save_kwargs["transverse_bond_schema_version"] = 1
                    save_kwargs["transverse_output_kind"] = "bond_static"
                    save_kwargs["transverse_bond_delta_r"] = \
                        np.asarray(topo.delta_r)
                    save_kwargs["transverse_bond_reverse"] = \
                        np.asarray(topo.reverse)
                    save_kwargs["transverse_bond_index_order"] = \
                        "m*norb**2 + a*norb + b"
                    save_kwargs["transverse_bond_max_shells"] = (
                        -1 if self.transverse_bond_max_shells is None
                        else int(self.transverse_bond_max_shells))
                    save_kwargs["transverse_spatial_shape"] = \
                        np.array(self.lattice.shape, dtype=np.int64)
                    save_kwargs["transverse_q_convention"] = \
                        "q_d = 2*pi*n_d/N_d, C-order flattening"
                    save_kwargs["transverse_spin_mode"] = self.spin_mode
                    save_kwargs["transverse_normalization"] = \
                        "per-site, 1/sqrt(Nvol) bilinears"
                np.savez(file_name, **save_kwargs)
                logger.info("save_results: save chiq in file {}".format(file_name))
            else:
                logger.info("save_results: chiq not calculated. skip")

        if "chi0q" in info_outputfile.keys():
            file_name = os.path.join(path_to_output, info_outputfile["chi0q"])
            save_kwargs = dict(
                chi0q = green_info["chi0q"],
                wavevector_unit = self.kvec,
                wavevector_index = self.wavenum_table,
                # spin-orbital axes are spin-block ordered (spin*norb+orb)
                index_convention = "spin_block",
                momentum_convention = MOMENTUM_CONVENTION,
                **_freq_meta_kwargs(green_info["chi0q"]),
            )
            np.savez(file_name, **save_kwargs)
            logger.info("save_results: save chi0q in file {}".format(file_name))

        pass

    def _init_wavevec(self):
        # wave vectors on sublatticed geometry
        def _klist(n):
            return np.roll( (np.arange(n)-(n//2)), -(n//2) )

        geom = self.param_ham["Geometry"]
        rvec = geom["rvec"]
        omg = np.dot(rvec[0], np.cross(rvec[1], rvec[2]))
        kvec = np.array([
            np.cross(rvec[(i+1)%3], rvec[(i+2)%3])/omg * 2*np.pi/self.lattice.shape[i]
            for i in range(3) ])

        self.kvec = kvec  # store reciprocal lattice vectors

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol

        self.wavenum_table = np.array([(i,j,k) for i in _klist(nx) for j in _klist(ny) for k in _klist(nz)])

        kx = _klist(nx)
        ky = _klist(ny)
        kz = _klist(nz)
        # Build (nx,ny,nz) grids for each k-component
        kx_g, ky_g, kz_g = np.meshgrid(kx, ky, kz, indexing='ij')
        # wtable[ix,iy,iz,:] = kvec[0]*kx + kvec[1]*ky + kvec[2]*kz
        wtable = (kx_g[..., np.newaxis] * kvec[0]
                + ky_g[..., np.newaxis] * kvec[1]
                + kz_g[..., np.newaxis] * kvec[2])
        self.wave_table = wtable.reshape(nvol, 3)

    def _find_index_range(self, freq_range):
        # decode matsubara frequency index list
        #   1. single index
        #   2. min, max, step
        #   3. keyword
        # note index n in [0 .. Nmat-1] corresponds to
        #   w_n = (2*n-Nmat) * pi / beta

        nmat = self.nmat

        if type(freq_range) == int:
            # e.g. freq_range = 0
            freq_index = [ freq_range ]
        elif type(freq_range) == list:
            if len(freq_range) == 1:
                # e.g. freq_range = [index]
                freq_index = [ i for i in freq_range ]
            elif len(freq_range) == 2:
                # e.g. freq_range = [min,max]
                freq_index = [ i for i in range(freq_range[0], freq_range[1]+1) ]
            elif len(freq_range) >= 3:
                # e.g. freq_range = [min,max,step]
                freq_index = [ i for i in range(freq_range[0], freq_range[1]+1, freq_range[2]) ]
            else:
                raise ValueError("invalid value for matsubara_frequency")
        elif type(freq_range) == str:
            freq_range = freq_range.lower()
            if freq_range == "all":
                freq_index = [ i for i in range(nmat) ]
            elif freq_range == "center" or freq_range == "zero":
                freq_index = [ nmat//2 ]
            elif freq_range == "none":
                freq_index = []
            else:
                raise ValueError("invalid value for matsubara_frequency")
        else:
            raise ValueError("invalid value type for matsubara_frequency")

        return freq_index

    @do_profile
    def read_init(self, info_inputfile):
        logger.info("RPA read initial configs")
        info = {}

        path_to_input = info_inputfile.get("path_to_input", "")

        if "chi0q_init" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["chi0q_init"])
            info["chi0q"] = self._read_chi0q(file_name)

        if "trans_mod" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["trans_mod"])
            info["trans_mod"] = self._read_trans_mod(file_name)

        if "green_init" in info_inputfile.keys():
            file_name = os.path.join(path_to_input, info_inputfile["green_init"])
            info["green_init"] = self._read_green(file_name)

        return info

    def _validate_chi0q_shape(self, chi0q, source):
        """Validate a supplied chi0q against the lattice/scheme and set
        spin_mode from its shape.

        Shared by every route that accepts a chi0q the solver did not
        compute itself: the file-based ``chi0q_init`` reader and the
        in-memory reuse of a ``green_info`` stored by a previous solve.
        The latter used to run NO validation, so a fresh solver reached
        the spin_mode dispatch unset and crashed with AttributeError
        (issue #109).
        """
        if not np.issubdtype(chi0q.dtype, np.number):
            raise ValueError(
                "chi0q from {}: non-numeric dtype {}".format(
                    source, chi0q.dtype))
        if self.calc_scheme == "general":
            if len(chi0q.shape) == 6:
                # spin-free or spinful
                #   shape = (nmat,nvol,nd,nd,nd,nd) where nd = norb or norb*nspin
                cs = chi0q.shape
                if not (cs[1] == self.lattice.nvol):
                    raise ValueError(
                        "chi0q from {}: lattice volume (shape {})".format(
                            source, chi0q.shape))
                nd = cs[2]
                if not (nd == self.nd or nd == self.norb):
                    raise ValueError(
                        "chi0q from {}: shape[2] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[3] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[3] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[4] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[4] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[5] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[5] (shape {})".format(
                            source, chi0q.shape))

                if nd == self.nd:
                    self.spin_mode = "spinful"
                else:
                    self.spin_mode = "spin-free"
            elif len(chi0q.shape) == 7:
                # spin-diagonal
                #   shape = (nblock,nmat,nvol,norb,norb,norb,norb)
                cs = chi0q.shape
                if not (cs[0] == 2):
                    raise ValueError(
                        "chi0q from {}: spin block (shape {})".format(
                            source, chi0q.shape))
                if not (cs[2] == self.lattice.nvol):
                    raise ValueError(
                        "chi0q from {}: lattice volume (shape {})".format(
                            source, chi0q.shape))
                nd = cs[3]
                if not (nd == self.norb):
                    raise ValueError(
                        "chi0q from {}: shape[3] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[4] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[4] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[5] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[5] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[6] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[6] (shape {})".format(
                            source, chi0q.shape))

                self.spin_mode = "spin-diag"
            else:
                raise ValueError(
                    "chi0q from {}: unexpected shape for general scheme: "
                    "{}".format(source, chi0q.shape))

        elif self.calc_scheme == "reduced":
            # reduced: shape = (nmat,nvol,nd,nd) where nd = norb or norb*nspin
            if len(chi0q.shape) == 4:
                # spin-free or spinful
                #   shape = (nmat,nvol,nd,nd) where nd = norb or norb*nspin
                cs = chi0q.shape
                if not (cs[1] == self.lattice.nvol):
                    raise ValueError(
                        "chi0q from {}: lattice volume (shape {})".format(
                            source, chi0q.shape))
                nd = cs[2]
                if not (nd == self.nd or nd == self.norb):
                    raise ValueError(
                        "chi0q from {}: shape[2] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[3] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[3] (shape {})".format(
                            source, chi0q.shape))

                if nd == self.nd:
                    self.spin_mode = "spinful"
                else:
                    self.spin_mode = "spin-free"
            elif len(chi0q.shape) == 5:
                # spin-diagonal
                #   shape = (nblock,nmat,nvol,norb,norb)
                cs = chi0q.shape
                if not (cs[0] == 2):
                    raise ValueError(
                        "chi0q from {}: spin block (shape {})".format(
                            source, chi0q.shape))
                if not (cs[2] == self.lattice.nvol):
                    raise ValueError(
                        "chi0q from {}: lattice volume (shape {})".format(
                            source, chi0q.shape))
                nd = cs[3]
                if not (nd == self.norb):
                    raise ValueError(
                        "chi0q from {}: shape[3] (shape {})".format(
                            source, chi0q.shape))
                if not (cs[4] == nd):
                    raise ValueError(
                        "chi0q from {}: shape[4] (shape {})".format(
                            source, chi0q.shape))

                self.spin_mode = "spin-diag"
            else:
                raise ValueError(
                    "chi0q from {}: unexpected shape for reduced scheme: "
                    "{}".format(source, chi0q.shape))
        else:
            raise ValueError(
                "unknown scheme: {}".format(self.calc_scheme))

        logger.debug("chi0q from {}: shape={}, spin_mode={}".format(
            source, chi0q.shape, self.spin_mode))

    def _read_chi0q(self, file_name):
        """Read chi0q from file and validate its shape.

        Parameters
        ----------
        file_name : str
            Path to the file containing chi0q data.

        Returns
        -------
        ndarray
            The loaded chi0q array.

        Raises
        ------
        ValueError
            If the loaded chi0q's SHAPE doesn't match the lattice and
            calculation scheme. (Standardized by the #109 validation
            rework: this previously surfaced as a bare AssertionError,
            which python -O silently disabled.) File-open and missing-key
            failures are reported separately and terminate the run.

        Notes
        -----
        Performs checks for:
        - Lattice volume consistency
        - Number of orbitals consistency 
        - Spin freedom consistency
        - Calculation scheme compatibility
        """
        logger.debug(">>> RPA._read_chi0q")

        try:
            logger.debug("read chi0q from {}".format(file_name))
            data = np.load(file_name)
            chi0q = data["chi0q"]
            logger.debug("chi0q: shape={}".format(chi0q.shape))
        except Exception as e:
            logger.error("read_chi0q failed: {}".format(e))
            sys.exit(1)

        # Stage 3 (design ir-matsubara-stage3.md Sec. 4): chi0q_init assumes
        # a uniform frequency grid; sparse-node files must not fall through
        # to the positional metadata handling below. (Outside the try block
        # on purpose -- the except above would turn this into sys.exit.)
        from hwave.solver.ir_axis import is_ir_native
        if is_ir_native(data):
            raise ValueError(
                "file '{}' holds sparse-IR node data "
                "(frequency_grid=sparse_ir_nodes); chi0q_init requires a "
                "uniform-grid file. Re-run FLEX with [mode.param] "
                "write_densified = true.".format(file_name))

        validate_chi0q_index_convention(
            data, self.ham_info.enable_spin_orbital, file_name)

        # Keep the frequency metadata of the file that PRODUCED this chi0q:
        # solve() passes an input chi0q through untouched, so save_results
        # must not relabel its axis with the current run's freq_index/nmat.
        # coeff_tail rides along for the same reason (issue #80).
        nfreq_file = (chi0q.shape[1]
                      if (chi0q.ndim in (5, 7) and chi0q.shape[0] == 2)
                      else chi0q.shape[0])
        self._chi0q_init_meta = _validate_chi0q_provenance(
            {"freq_index": (data["freq_index"]
                            if "freq_index" in data else None),
             "nmat": data["nmat"] if "nmat" in data else None,
             "coeff_tail": (data["coeff_tail"]
                            if "coeff_tail" in data else None),
             "tail_endpoint": (data["tail_endpoint"]
                               if "tail_endpoint" in data else None)},
            nfreq_file, source=file_name)
        # bind the metadata to the LOADED tensor: if the caller replaces
        # info["chi0q"] with a same-shaped array before solve(), the file's
        # axis must not be inherited by the replacement (round-6 review)
        self._chi0q_init_meta["fingerprint"] = _chi0q_fingerprint(chi0q)
        # The tail correction changes chi0q at O(1); a config whose
        # coeff_tail differs from the file's describes different physics
        # than this chi0q_init actually contains (issue #80).
        file_tail = self._chi0q_init_meta["coeff_tail"]
        if file_tail is not None and file_tail != self.coeff_tail:
            logger.warning(
                "chi0q_init file '{}' was produced with coeff_tail = {} but "
                "the current config uses coeff_tail = {}; the loaded chi0q "
                "keeps the file's tail treatment, which is NOT comparable "
                "with a recomputation under this config.".format(
                    file_name, file_tail, self.coeff_tail))
        # Endpoint-convention gate (issue #134, fail-closed like the
        # momentum-convention gate); shared with the in-memory reuse
        # route -- see _enforce_tail_endpoint.
        _enforce_tail_endpoint(self._chi0q_init_meta,
                               "file '{}'".format(file_name))

        self._validate_chi0q_shape(chi0q, source=file_name)

        # Fourier-sign provenance gate (issue #133), AFTER the IR rejection
        # and shape validation so a malformed file gets its own diagnosis
        # first. The q axis is chosen STRUCTURALLY from the now-validated
        # layout, never by dimension-size search (round-2 review: a
        # spin-diag (2, nfreq, nvol, ...) file with nfreq == nvol would
        # otherwise be probed on its frequency axis and slip through):
        # raw layouts put the flattened volume on axis 1, spin-diag
        # layouts (ndim 5/7, leading 2) on axis 2.
        _qax = 2 if (chi0q.ndim in (5, 7) and chi0q.shape[0] == 2) else 1
        validate_momentum_convention(data, file_name, chi0q, _qax,
                                     self.lattice.shape)

        logger.info("read_chi0q: shape={}, spin_mode={}".format(chi0q.shape, self.spin_mode))

        return chi0q

    def _read_trans_mod(self, file_name):
        """Read transfer integral modifications from file.

        Parameters
        ----------
        file_name : str
            Path to the file containing transfer integral modifications.

        Returns
        -------
        ndarray
            The transfer integrals in k-space.

        Raises
        ------
        RuntimeError
            If file reading fails.

        Notes
        -----
        - Converts layout if sublattice exists
        - Performs Fourier transform to k-space
        """
        logger.debug(">>> RPA._read_trans_mod")

        try:
            logger.info("read trans_mod from {}".format(file_name))
            data = np.load(file_name)
            tab_r = data["trans_mod"]
            logger.debug("read_trans_mod: shape={}".format(tab_r.shape))
        except Exception as e:
            logger.error("read_trans_mod failed: {}".format(e))
            sys.exit(1)

        expected = (self.lattice.cellvol, self.ns * self.norb_orig, self.ns * self.norb_orig)
        if tab_r.shape != expected:
            raise ValueError(
                "trans_mod array shape {} does not match expected {} "
                "(cellvol, ns*norb_orig, ns*norb_orig)".format(tab_r.shape, expected))

        if self.ham_info.enable_spin_orbital:
            # UHFk writes the SO trans_mod with the orbital axis in INTERLEAVED
            # order (index = 2*orb + spin), matching the SO Transfer file, but
            # _reshape_green / the (ns, norb) reshapes downstream consume H0 in
            # SPIN-BLOCK order (index = spin*norb_phys + orb). Reorder both
            # orbital axes interleaved->spin-block, mirroring _make_ham_trans's
            # remap of the Transfer file. For norb_phys=1 this is the identity,
            # so non-multi-orbital SO and non-SO paths are unaffected.
            norb_phys = self.norb_orig            # pre-fold physical orbital count
            nd0 = tab_r.shape[-1]                 # = 2 * norb_phys
            inv = [2 * (j % norb_phys) + (j // norb_phys) for j in range(nd0)]
            tab_r = tab_r[:, inv, :][:, :, inv]   # interleaved -> spin-block

        if self.lattice.has_sublattice:
            # use reshape green to convert layout (Hamiltonian/transfer convention)
            tab_r = self._reshape_green(tab_r, hamiltonian=True)

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd

        # e^{+ikR} (issue #133), consistent with _make_ham_trans
        tab_k = (np.fft.ifftn(tab_r.reshape(nx,ny,nz,nd,nd), axes=(0,1,2))
                 * nvol).reshape(nvol,nd,nd)

        return tab_k

    def _read_green(self, file_name):
        logger.debug(">>> RPA._read_green")
        try:
            logger.info("read green from {}".format(file_name))
            data = np.load(file_name)
        except Exception as e:
            logger.error("read_green failed: {}".format(e))
            sys.exit(1)

        # Fourier-sign provenance (issue #133, round-9 review): validate a
        # PRESENT marker strictly, on every branch including the
        # green_sublattice early return below. Unmarked files stay
        # accepted without a content check: UHFk has always written its
        # green with the documented e^{+ikR} sign, so legacy files carry
        # correct labels -- but a file explicitly recording a DIFFERENT
        # convention must not be consumed.
        check_momentum_marker(data, file_name)

        nvol = self.lattice.nvol
        nd = self.nd

        # Sublattice green_init needs folding from the original (deflated) basis
        # into the supercell basis. The "green" key carries a fold-sign
        # convention: since PR #35 UHFk writes it with the Green convention
        # (within-cell offset on the FIRST orbital slot) and tags the file with
        # green_convention="green_slot_first". Files written before #35 stored
        # "green" with the opposite (Hamiltonian) convention and carry no tag,
        # so folding their "green" key here would be SILENTLY WRONG (issue #36).
        #
        # The "green_sublattice" array, when present, is the already-folded
        # internal Green and is correct regardless of the deflate convention, so
        # prefer it. Fall back to folding "green" only when the convention is
        # unambiguous (new-style tag present); otherwise fail loudly.
        #
        # NOTE: green_sublattice is UHFk's self.Green, stored with the orbital
        # axis in INTERLEAVED order in spin-orbital mode, whereas RPA works in
        # SPIN-BLOCK order. For norb_phys>1 the two differ by a permutation that
        # the "green" path applies (interleaved->spin-block) but a direct read of
        # green_sublattice would not. So the green_sublattice shortcut is taken
        # only in non-SO mode; SO files go through the tag-gated "green" path,
        # which performs the remap.
        if self.lattice.has_sublattice:
            if (not self.ham_info.enable_spin_orbital
                    and "green_sublattice" in data.files):
                logger.debug("read_green: use green_sublattice (fold-convention independent)")
                gsub = data["green_sublattice"]
                if gsub.ndim == 5:
                    lvol, s1, o1, s2, o2 = gsub.shape
                    gsub = gsub.reshape(lvol, s1 * o1, s2 * o2)
                expected = (nvol, nd, nd)
                if gsub.shape != expected:
                    raise ValueError(
                        "green_sublattice array shape {} does not match expected "
                        "{} (nvol, nd, nd)".format(gsub.shape, expected))
                return gsub.reshape(nvol, nd, nd)

            tag = str(data["green_convention"]) if "green_convention" in data.files else None
            if tag != "green_slot_first":
                if "green_sublattice" in data.files:
                    # SO mode reached here: green_sublattice exists but is stored
                    # interleaved (RPA is spin-block) and folded in UHFk's
                    # orbital order, so it cannot be consumed directly. The "green"
                    # key needs the tag to be folded with the correct sign.
                    reason = ("green_sublattice is present but cannot be used "
                              "directly in spin-orbital mode (it is interleaved "
                              "and folded in UHFk order), and the 'green' key "
                              "carries no green_convention tag")
                else:
                    reason = ("no green_convention tag and no green_sublattice "
                              "array")
                logger.error(
                    "read_green: {} is a sublattice green file with an "
                    "ambiguous/old fold convention (green_convention={}; {}). "
                    "Folding its 'green' key could be silently wrong (see issue "
                    "#36). Regenerate green_init with the current UHFk.".format(
                        file_name, tag, reason))
                sys.exit(1)

        try:
            green = data["green"]
            logger.debug("read_green: shape={}".format(green.shape))
        except Exception as e:
            logger.error("read_green failed: {}".format(e))
            sys.exit(1)

        # UHFk saves green as 5D (Lvol, ns, norb_orig, ns, norb_orig); collapse
        # the (ns, norb) pairs into single orbital axes to get the (Lvol, nd0, nd0)
        # layout the rest of this method expects. trans_mod is already saved 3D.
        if green.ndim == 5:
            lvol, s1, o1, s2, o2 = green.shape
            green = green.reshape(lvol, s1 * o1, s2 * o2)

        # green_init is produced by UHFk's _save_green (saves self.Green), the
        # DEFLATED pre-fold green for sublattice -- the same (cellvol, nd0, nd0)
        # layout as trans_mod, and consumed in spin-block order. It does NOT
        # share trans_mod's fold sign: a Green function carries the within-cell
        # offset on the first orbital slot, so the sublattice fold below uses
        # the Green convention (default), the inverse of UHFk._deflate_green.
        expected = (self.lattice.cellvol, self.ns * self.norb_orig, self.ns * self.norb_orig)
        if green.shape != expected:
            raise ValueError(
                "green array shape {} does not match expected {} "
                "(cellvol, ns*norb_orig, ns*norb_orig)".format(green.shape, expected))

        if self.ham_info.enable_spin_orbital:
            # UHFk writes the SO green with the orbital axis in INTERLEAVED order
            # (index = 2*orb + spin); _calc_trans_mod consumes green_init in
            # SPIN-BLOCK order (index = spin*norb_phys + orb), using self.norb /
            # self.nd. Reorder both orbital axes interleaved->spin-block, mirroring
            # _read_trans_mod (and _make_ham_trans's remap of the Transfer file).
            # For norb_phys=1 this is the identity.
            norb_phys = self.norb_orig            # pre-fold physical orbital count
            nd0 = green.shape[-1]                 # = 2 * norb_phys
            inv = [2 * (j % norb_phys) + (j // norb_phys) for j in range(nd0)]
            green = green[:, inv, :][:, :, inv]   # interleaved -> spin-block

        if self.lattice.has_sublattice:
            # use reshape green to convert layout
            green = self._reshape_green(green)

        nvol = self.lattice.nvol
        nd = self.nd
        return green.reshape(nvol,nd,nd)

    def _reshape_green(self, green_, hamiltonian=False):
        # convert green function into sublattice
        #
        # ``hamiltonian`` selects the cross-sublattice slot convention, mirroring
        # UHFk._deflate_green. A Green function saved by UHFk._save_green carries
        # the within-cell offset on the FIRST orbital slot, so folding it back
        # uses dr = ri - rj (the default), the inverse of that deflate.
        # Hamiltonian/transfer quantities (trans_mod) carry the offset on the
        # SECOND slot and fold with dr = rj - ri; pass hamiltonian=True for those.

        Lx,Ly,Lz = self.lattice.cellshape
        Lvol = Lx * Ly * Lz
        Bx,By,Bz = self.lattice.subshape
        Bvol = Bx * By * Bz
        Nx,Ny,Nz = self.lattice.shape
        Nvol = Nx * Ny * Nz

        norb_orig = self.norb_orig
        norb = self.norb
        ns = self.ns

        # check array size
        #assert(green.shape == (Lvol,ns,norb_orig,ns,norb_orig))
        green = green_.reshape(Lvol,ns,norb_orig,ns,norb_orig)

        def _pack_index(x, n):
            _ix, _iy, _iz = x
            _nx, _ny, _nz = n
            return _ix + _nx * (_iy + _ny * (_iz))

        def _unpack_index(x, n):
            _nx, _ny, _nz = n
            _ix = x % _nx
            _iy = (x // _nx) % _ny
            _iz = (x // (_nx * _ny)) % _nz
            return (_ix, _iy, _iz)

        # Build index mapping tables (vectorized)
        # Supercell site indices.
        # Flat SITE indices are C-order (z fastest): isite = iz + Nz*(iy + Ny*ix),
        # matching the data layout used everywhere else -- UHFk._deflate_green's
        # _pack_site/_unpack_site and RPA's reshape(nx,ny,nz) lattice flattening.
        # (1D folds coincide with Fortran order, which masked this for years.)
        isite_arr = np.arange(Nvol)
        izz = isite_arr % Nz
        iyy = (isite_arr // Nz) % Ny
        ixx = (isite_arr // (Nz * Ny)) % Nx
        ix0 = ixx * Bx  # (Nvol,)
        iy0 = iyy * By
        iz0 = izz * Bz

        # Orbital decomposition
        orb_arr = np.arange(norb)
        a_arr = orb_arr % norb_orig          # original orbital index
        ri_arr = orb_arr // norb_orig         # sublattice index
        rix = ri_arr % Bx
        riy = (ri_arr // Bx) % By
        riz = (ri_arr // (Bx * By)) % Bz

        # Compute jsite for all (isite, aa, bb) combinations (rix indexes the
        # within-cell coordinate of each orbital slot).
        # Hamiltonian/transfer convention: drx[aa,bb] = rix[bb] - rix[aa].
        # Green convention (default): drx[aa,bb] = rix[aa] - rix[bb], i.e. the
        # offset sits on the first slot -- the inverse of UHFk._deflate_green.
        if hamiltonian:
            drx = rix[np.newaxis, :] - rix[:, np.newaxis]  # (norb, norb)
            dry = riy[np.newaxis, :] - riy[:, np.newaxis]
            drz = riz[np.newaxis, :] - riz[:, np.newaxis]
        else:
            drx = rix[:, np.newaxis] - rix[np.newaxis, :]  # (norb, norb)
            dry = riy[:, np.newaxis] - riy[np.newaxis, :]
            drz = riz[:, np.newaxis] - riz[np.newaxis, :]

        # ix[isite, aa, bb] = (ix0[isite] + drx[aa, bb]) % Lx
        jx = (ix0[:, np.newaxis, np.newaxis] + drx[np.newaxis, :, :]) % Lx  # (Nvol, norb, norb)
        jy = (iy0[:, np.newaxis, np.newaxis] + dry[np.newaxis, :, :]) % Ly
        jz = (iz0[:, np.newaxis, np.newaxis] + drz[np.newaxis, :, :]) % Lz
        # Pack target SITE in C-order (z fastest), consistent with the unpack
        # above and UHFk._deflate_green's _pack_site.
        jsite_map = jz + Lz * (jy + Ly * jx)  # (Nvol, norb, norb)

        # Source orbital indices: a_arr[aa], a_arr[bb]
        a_src = a_arr  # (norb,) - maps aa -> original orbital a
        b_src = a_arr  # (norb,) - maps bb -> original orbital b

        # Gather using advanced indexing
        # green[jsite, s, a, t, b] -> green_sub[isite, s, aa, t, bb]
        green_sub = green[
            jsite_map[:, :, :, np.newaxis, np.newaxis],   # (Nvol, norb, norb, 1, 1)
            np.arange(ns)[np.newaxis, np.newaxis, np.newaxis, :, np.newaxis],  # s
            a_src[np.newaxis, :, np.newaxis, np.newaxis, np.newaxis],          # a
            np.arange(ns)[np.newaxis, np.newaxis, np.newaxis, np.newaxis, :],  # t
            b_src[np.newaxis, np.newaxis, :, np.newaxis, np.newaxis],          # b
        ]  # shape: (Nvol, norb, norb, ns, ns)

        # Transpose to match (Nvol, ns, norb, ns, norb)
        green_sub = green_sub.transpose(0, 3, 1, 4, 2)

        return green_sub

    def _calc_trans_mod(self, g0):
        logger.debug(">>> RPA._calc_trans_mod")

        nx,ny,nz = self.lattice.shape
        nvol = self.lattice.nvol
        nd = self.nd
        norb = self.norb

        gg = g0[0]
        ww = self.ham_info.ham_inter_r.reshape(nvol,nd,nd,nd,nd)

        hh1 = np.einsum('rbacd,cd->rab', ww, gg)
        hh2 = np.einsum('rcdab,dc->rab', ww, gg)
        hh3 = np.sum(hh1+hh2, axis=0)/2

        if self.ham_info.enable_spin_orbital:
            H0r = self.ham_info.ham_trans_r.reshape(nvol,nd,nd)
        else:
            H0r = np.einsum('kab,st->ksatb',
                            self.ham_info.ham_trans_r.reshape(nvol,norb,norb),
                            np.eye(2)).reshape(nvol,nd,nd)
        H0r[0] += hh3

        # e^{+ikR} (issue #133), consistent with _make_ham_trans
        H0 = (np.fft.ifftn(H0r.reshape(nx,ny,nz,nd,nd), axes=(0,1,2))
              * nvol).reshape(nvol,nd,nd)
        return H0

    @do_profile
    def _calc_epsilon_k(self, green_info):
        logger.debug(">>> RPA._calc_epsilon_k")

        nvol = self.lattice.nvol

        # Find transfer term H0(k)
        if "trans_mod" in green_info:
            logger.debug("calc_epsilon_k: use trans_mod")
            H0 = green_info['trans_mod']
            do_spin_orbital = True

        elif "green_init" in green_info:
            logger.debug("calc_epsilon_k: use initial green")
            H0 = self._calc_trans_mod(green_info["green_init"])
            do_spin_orbital = True

        else:
            H0 = self.ham_info.ham_trans_q
            do_spin_orbital = self.ham_info.enable_spin_orbital

        if do_spin_orbital:
            # H0(k) = T_{a~b~}(k) + h sigma_z_{ss'} H_{ab}(k)
            #   T_{a~b~} with extended orbital index
            #   H_{ab} with bare orbital index
            #   sigma_z = diag(1,-1) Pauli matrix
            #   h coefficient

            if self.ham_info.ham_extern_q is not None:
                H1 = self.ham_info.ham_extern_q * self.ext
                Sz = np.diag((1,-1))
                H0 += np.einsum('kab,st->ksatb', H1, Sz).reshape(H0.shape)

            # check if diagonal
            ns = self.ns
            norb = self.norb

            Htmp = H0.reshape(nvol,ns,norb,ns,norb)
            if np.allclose(Htmp[:,0,:,1,:], 0) and np.allclose(Htmp[:,1,:,0,:], 0):
                if np.allclose(Htmp[:,0,:,0,:], Htmp[:,1,:,1,:]):
                    logger.info("H is spin-free")
                    H0 = Htmp[:,0,:,0,:].reshape(1,nvol,norb,norb)
                    self.spin_mode = "spin-free"
                else:
                    logger.info("H is spin-diagnoal")
                    Hnew = np.zeros((2,nvol,norb,norb), dtype=np.complex128)
                    Hnew[0] = Htmp[:,0,:,0,:]
                    Hnew[1] = Htmp[:,1,:,1,:]
                    H0 = Hnew
                    self.spin_mode = "spin-diag"
            else:
                logger.debug("H is spinful")
                H0 = H0.reshape(1,*H0.shape)
                self.spin_mode = "spinful"

        else:
            # H0(k) = T_{ab}(k) x 1_{ss'} + h Sz_{ss'} H_{ab}(k)
            #   T_{a~b~} with bare orbital index

            if self.ham_info.ham_extern_q is not None:
                H1 = self.ham_info.ham_extern_q * self.ext

                ns = self.ns
                norb = self.norb

                Hnew = np.zeros((2,nvol,norb,norb), dtype=np.complex128)
                Hnew[0] = H0 + H1
                Hnew[1] = H0 - H1
                H0 = Hnew
                self.spin_mode = "spin-diag"
                logger.info("H is spin-diagnoal")
            else:
                logger.debug("H is spin-free")
                H0 = H0.reshape(1,*H0.shape)
                self.spin_mode = "spin-free"

        # diagonalize H0(k) with optional block decomposition
        nblock_spin = H0.shape[0]
        nd_block = H0.shape[-1]
        blocks = self._find_block_diagonal(H0.reshape(nblock_spin * nvol, nd_block, nd_block))

        if blocks is not None and len(blocks) > 1:
            logger.info("_calc_epsilon_k: orbital block structure detected, "
                        "nd={} -> {} blocks of sizes {}".format(
                            nd_block, len(blocks), [len(b) for b in blocks]))
            w = np.zeros((nblock_spin, nvol, nd_block), dtype=np.float64)
            v = np.zeros((nblock_spin, nvol, nd_block, nd_block), dtype=np.complex128)
            col_offset = 0
            for blk_idx in blocks:
                idx = np.array(blk_idx)
                ix = np.ix_(idx, idx)
                H0_blk = H0[:, :, ix[0], ix[1]]
                wb, vb = np.linalg.eigh(H0_blk)
                nb = len(idx)
                w[:, :, col_offset:col_offset + nb] = wb
                v[np.ix_(np.arange(nblock_spin), np.arange(nvol), idx,
                  np.arange(col_offset, col_offset + nb))] = vb
                col_offset += nb
        else:
            w,v = np.linalg.eigh(H0)

        self.H0_eigenvalue = w
        self.H0_eigenvector = v

    @do_profile
    def _find_mu(self, Ncond, T):
        logger.debug(">>> RPA._find_mu")

        from scipy import optimize

        # load eigenvalues (eigenvectors not needed thanks to unitarity)
        w = self.H0_eigenvalue
        # fetch parameters
        ene_cutoff = self.ene_cutoff

        ev = np.sort(w.flatten())
        occupied_number = Ncond

        def _fermi(t, mu, ev):
            w_ = (ev - mu) / t
            mask_ = w_ < ene_cutoff
            w1_ = np.where( mask_, w_, 0.0 )
            v1_ = 1.0 / (1.0 + np.exp(w1_))
            v_ = np.where( mask_, v1_, 0.0 )
            return v_

        # Exploit unitarity of eigenvectors:
        # Tr[V† diag(f) V] = sum_l f(ε_l)
        # This eliminates the O(nd²) einsum per iteration.
        def _calc_delta_n(mu):
            ff = _fermi(T, mu, w)
            return np.sum(ff).real - occupied_number

        # find mu s.t. <n>(mu) = N0
        is_converged = False
        if (_calc_delta_n(ev[0]) * _calc_delta_n(ev[-1])) < 0.0:
            logger.debug("RPA._find_mu: try bisection")
            mu, r = optimize.bisect(_calc_delta_n, ev[0], ev[-1], full_output=True, disp=False)
            is_converged = r.converged
        if not is_converged:
            logger.debug("RPA._find_mu: try newton")
            mu, r = optimize.newton(_calc_delta_n, ev[0], full_output=True)
            is_converged = r.converged
        if not is_converged:
            logger.error("RPA._find_mu: not converged. abort")
            sys.exit(1)

        logger.info("RPA._find_mu: mu = {}".format(mu))

        dist = _fermi(T, mu, w)

        return dist, mu

    @do_profile
    def _calc_green(self, beta, mu):
        """
        Thin wrapper around ``green.build_green`` (the shared Green's
        function builder used by both RPA and sc.py's bond carrier).

        ew: eigenvalues  ew[g,k,i]    i-th eigenvalue of wavenumber k in block g
        ev: eigenvectors ev[g,k,a,i]  i-th eigenvector corresponding to ew[g,k,i]
        beta: inverse temperature
        mu: chemical potential

        Returns
        -------
        (green, green_tail) : the deflated Green's function (with the
        coeff_tail/(iw_n) term subtracted for faster tau-space decay) and
        its analytic high-frequency tail. This preserves the legacy
        contract exactly: when ``coeff_tail == 0`` ``build_green`` returns
        ``green0_tail=None``, but callers such as ``_calc_chi0q_transverse``
        access ``green0_tail.shape`` unconditionally, so this wrapper
        materializes a shape-matched all-zero tail in that case.

        ``want_full=False`` is passed to ``build_green``: this wrapper only
        ever returns the deflated Green/tail pair, so the full-size
        ``full_kw`` reconstruction (one whole extra Green-sized allocation
        when ``coeff_tail != 0``) is never even built here.
        """
        logger.debug(">>> RPA._calc_green")

        # load eigenvalues and eigenvectors
        ew = self.H0_eigenvalue
        ev = self.H0_eigenvector

        nblock, nvol, nd = ew.shape
        assert nvol == self.lattice.nvol

        _full, green, green_tail = green_mod.build_green(
            ew, ev, mu, beta, self.nmat, self.coeff_tail, want_full=False)

        if green_tail is None:
            xp = _bk.array_module_of(ew)
            green_tail = xp.zeros_like(green)

        return green, green_tail

    @do_profile
    def _calc_chi0q(self, green_kw, green0_tail, beta):
        """Calculate the bare susceptibility chi0q.

        Parameters
        ----------
        green_kw : ndarray
            Green's function in k-space and Matsubara frequency.
            Shape: [g,l,k,a,b] where:
            - g: block index
            - l: Matsubara frequency index
            - k: wave number index
            - a,b: orbital and spin indices
        green0_tail : ndarray
            High-frequency tail correction for Green's function.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            The calculated chi0q.

        Notes
        -----
        Validates the input shapes (see the inline guards below), then
        delegates the dense-grid bubble calculation to
        ``bubble.dense_bubble``. The guards here duplicate some of the
        kernel's own validation deliberately -- their messages are the
        user-facing diagnostics for this solver and are kept even though
        the kernel would also catch the malformed input.
        """
        logger.debug(">>> RPA._calc_chi0q")

        workers = getattr(self, "fft_workers", 1)

        nx,ny,nz = self.lattice.shape
        #nvol = self.lattice.nvol

        nblock,nmat,nvol,nd,nd2 = green_kw.shape
        # ValueError, not assert (issue #125, the longitudinal analogue of
        # the transverse block-count fix): these validate INPUT array data,
        # and a bare assert disappears under python -O. A wrong frequency
        # count would then proceed into the kernel and produce
        # plausible-looking wrong output; a wrong volume fails later at
        # the lattice reshape, but with a diagnostic that names neither
        # the axis nor the expectation.
        if nvol != self.lattice.nvol:
            raise ValueError(
                "chi0q kernel: Green's function volume axis ({}) does not "
                "match the lattice ({})".format(nvol, self.lattice.nvol))
        if nmat != self.nmat:
            raise ValueError(
                "chi0q kernel: Green's function frequency axis ({}) does "
                "not match Nmat ({})".format(nmat, self.nmat))
        if nblock not in (1, 2):
            # 1 = spin-free/spinful, 2 = spin-diag. Anything else is
            # malformed input: nblock=0 returned a finite EMPTY result
            # and nblock=3 a finite three-block one, which the caller's
            # block handling would then silently truncate (round-5
            # review reproduced both).
            raise ValueError(
                "chi0q kernel: Green's function block axis ({}) must be "
                "1 (spin-free/spinful) or 2 (spin-diag)".format(nblock))
        if nd != nd2 or nd < 1:
            raise ValueError(
                "chi0q kernel: orbital axes must be square and nonempty, "
                "got ({}, {})".format(nd, nd2))
        if green0_tail.shape != green_kw.shape:
            # STRUCTURAL pairing check only: shape equality cannot prove
            # the tail came from the same _calc_green call (a stale
            # same-shaped tail passes and shifts chi0q). Provenance
            # machinery is deliberately out of scope here -- unlike
            # chi0q reuse (#116), this pair never crosses the
            # green_info/file boundary: every caller passes the two
            # arrays one _calc_green call returned together. The guard
            # exists because a same-SIZE tail of a different shape would
            # be silently reshaped below and corrupt chi0q with finite,
            # plausible-looking values (round-3 review reproduced it).
            raise ValueError(
                "chi0q kernel: green0_tail shape {} does not match the "
                "Green's function {} -- the tail must be the paired "
                "array from the same _calc_green call".format(
                    green0_tail.shape, green_kw.shape))

        # Delegate to the shared bubble kernel (spec: "Module layout").
        # green0_tail is forwarded AS-IS -- an all-zero tail is not
        # special-cased here; the kernel's own data-driven tail gate
        # (mirroring the legacy _tail_on predicate) handles it.
        return bubble.dense_bubble(
            green_kw, green0_tail, beta,
            spatial_shape=(nx, ny, nz),
            scheme="reduced" if self.enable_reduced else "general",
            workers=workers)

    def _calc_chi0q_transverse(self, green_kw, green0_tail, beta):
        """Calculate the transverse bare susceptibility chi0_+-(q,iω).

        chi0_+-[a,c,b,d](r,τ) = -G_↑[a,b](r,τ) * G_↓[d,c](-r,-τ)

        This crosses spin-up and spin-down Green's functions, unlike the
        longitudinal chi0 which uses same-spin products.

        Parameters
        ----------
        green_kw : ndarray, shape (2, nmat, nvol, norb, norb)
            Green's function with block 0=↑, block 1=↓.
        green0_tail : ndarray
            High-frequency tail correction.
        beta : float
            Inverse temperature.

        Returns
        -------
        ndarray
            Transverse chi0_+- with block dimension removed (shape depends
            on enable_reduced).

        Notes
        -----
        Validates the input shapes (see the inline guards below), then
        delegates the cross-block bubble calculation to
        ``bubble.transverse_bubble`` -- mirroring ``_calc_chi0q``'s
        wrapper (spec: "Module layout"), the guards here duplicate the
        kernel's own validation deliberately, keeping this solver's
        user-facing diagnostics even though the kernel would also catch
        the malformed input. Two of them are worth calling out
        specifically: the ``nblock == 2`` check keeps its own
        transverse-specific message (distinct from the kernel's own,
        looser ``nblock in {1, 2}`` guard), and the ``green0_tail`` shape
        guard stays UNCONDITIONAL -- this wrapper's contract still
        requires a real tail array; the kernel's ``green0_tail=None``
        convenience is NOT exposed through it.
        """
        logger.debug(">>> RPA._calc_chi0q_transverse")

        workers = getattr(self, "fft_workers", 1)

        nx, ny, nz = self.lattice.shape
        nblock, nmat, nvol, nd, nd2 = green_kw.shape
        if nblock != 2:
            # ValueError, not assert: an assert disappears under python -O,
            # and a wrong block count here would silently ignore the extra
            # blocks (measured: nblock=3 was accepted with the third block
            # dropped) -- plausible-looking wrong output.
            raise ValueError(
                "transverse chi0 requires a spin-diag Green's function with "
                "exactly 2 spin blocks (G_up, G_down), got nblock={}".format(
                    nblock))
        if nvol != self.lattice.nvol:
            raise ValueError(
                "transverse chi0: Green's function volume axis ({}) does "
                "not match the lattice ({})".format(nvol, self.lattice.nvol))
        if nmat != self.nmat:
            raise ValueError(
                "transverse chi0: Green's function frequency axis ({}) "
                "does not match Nmat ({})".format(nmat, self.nmat))
        if nd != nd2 or nd < 1:
            raise ValueError(
                "transverse chi0: orbital axes must be square and "
                "nonempty, got ({}, {})".format(nd, nd2))
        if green0_tail.shape != green_kw.shape:
            # same paired-tail invariant as the longitudinal kernel
            raise ValueError(
                "transverse chi0: green0_tail shape {} does not match the "
                "Green's function {} -- the tail must be the paired "
                "array from the same _calc_green call".format(
                    green0_tail.shape, green_kw.shape))

        # Delegate to the shared bubble kernel (spec: "Module layout").
        return bubble.transverse_bubble(
            green_kw, green0_tail, beta,
            spatial_shape=(nx, ny, nz),
            scheme="reduced" if self.enable_reduced else "general",
            workers=workers)

    def _assemble_transverse_vertex(self, ham_orig):
        """Build the transverse vertex ham_pm from the interaction tensor.

        ``ham_orig`` is in the HAMILTONIAN convention (the tensor as
        ``_make_ham_inter`` builds it), NOT the bubble-pair convention
        the ring solve consumes: the block reads and the crossing done
        below already re-pair it (issue #139).

        The vertex draws on exactly two spin blocks of the four-index tensor.
        Measured block occupancy (on-site fixture):

            CoulombIntra   uudd, dduu        Hund       uuuu, dddd
            CoulombInter   uuuu/uudd/dduu/dddd
            Exchange       udud, dudu        PairLift   uddu, duud
            Ising          uuuu/uudd/dduu/dddd          PairHop    uudd, dduu

        The same-spin blocks must NOT enter: a same-spin interaction cannot
        connect the up and down propagators of the transverse loop, so it
        contributes self-energy but no vertex (measured: 1.4e-5 against 0.63
        for the cross-spin part). PairLift sits outside both used blocks, but
        harmlessly -- its transverse vertex vanishes identically.

        Each block is averaged over the two redundant declarations of the same
        operator, with the MEAN, which is what an interaction file already
        means in this codebase (`uhfk.py` builds every two-body table as
        `(jab_r + jba)/2`; the two solvers agree exactly on symmetric input
        and diverge only on asymmetric declarations, #106). Without it, an
        on-site CoulombInter with V_01 = +1, V_10 = -1 -- an identically ZERO
        Hamiltonian, since n_a n_b = n_b n_a -- produced a vertex of +-1.

        The partner in the mean is decided PER SLOT FAMILY, not per block,
        because the redundancy being quotiented is a property of the slots:

        * PALINDROMIC pairs, ((a,a),(b,b)): the two declarations multiply the
          SAME operator (n_a n_b = n_b n_a for the density types; X_ab = X_ba
          for Exchange, since different-spin bilinears commute). The physical
          coefficient is the sum, real for any Hermitian declaration, and the
          declared imaginary parts multiply the null operator -- so the mean
          must NOT conjugate, which is exactly what drops them. A conjugated
          mean wrongly kept an imaginary part both for complex Hermitian-closed
          Exchange (spin-flip block) and for complex Hermitian-closed
          CoulombInter / Ising (cross block): measured -0.7 -/+ 0.4i in each
          case where the physical coupling is -0.7 (or +0.7 for Ising).

        * HETEROGENEOUS pairs, ((a,b),(b,a)) with a != b: the two declarations
          are coefficients of HERMITIAN-PARTNER operators (PairHop:
          Y_ab^dagger = Y_ba), the complex phase is physical, and the mean must
          conjugate -- that keeps a complex Hermitian-closed PairHop
          (P, conj(P)) at its full complex value (issue #100); a plain mean
          collapsed it to Re(P).

        A per-block rule looked sufficient at first because PairHop is the only
        heterogeneous occupant, but density-density types live on the SAME
        cross block in palindromic slots, so the block is not the right
        granularity.

        Resulting vertices, verified against exact diagonalization end to end
        with one common scale across all seven types (residuals at the
        imaginary-time discretisation floor):

            CoulombIntra  -U            Ising     +I
            CoulombInter  -U'           PairLift   0
            Hund           0            PairHop   -P
            Exchange      -(J + J^T)/2

        A non-Hermitian-closed TABLE reaching this builder directly (a
        reader-bypassing internal path) is projected onto its Hermitian
        part; FILE input cannot reach here unclosed since the #93
        read-time validation.

        The orbital ordering is one of two the measurements cannot separate;
        they agree exactly on every physically valid on-site input and differ
        only for off-site terms (rejected) and non-Hermitian declarations.
        """
        xp = _bk.array_module_of(ham_orig)
        norb = self.norb
        nd = self.norb * self.ns
        nvol = self.lattice.nvol
        ham_4d = ham_orig.reshape(nvol, nd, nd, nd, nd)
        cross_block = ham_4d[:, norb:, norb:, :norb, :norb]
        spin_flip_block = ham_4d[:, :norb, norb:, :norb, norb:]
        pair_swap = (0, 3, 4, 1, 2)

        # slot-family mask over (i, j, k, l): True where both pairs are
        # palindromic (i == j and k == l)
        oi = xp.arange(norb)
        palin = ((oi[:, None, None, None] == oi[None, :, None, None])
                 & (oi[None, None, :, None] == oi[None, None, None, :]))[None]

        # Swapping the two bilinears of an OFF-SITE term also reverses the
        # displacement, r -> -r, so the same-operator partner lives at -q. At
        # fixed q the plain swap manufactured false cancellations: the
        # off-site declaration V_ab(+r) = +1, V_ba(+r) = -1 -- a genuinely
        # NONZERO Hamiltonian, since its redundancy partner would be
        # V_ba(-r) -- assembled to an identically zero vertex, slipped
        # through the guard, and its physics was silently dropped, while the
        # truly redundant pair (+r, a, b) / (-r, b, a) was rejected. The
        # heterogeneous (Hermitian-partner) family needs no q reversal:
        # conj at fixed q is the Fourier image of {conjugate, r -> -r}, which
        # is exactly the h.c. redundancy. On-site input is q-independent and
        # unaffected either way.
        nx, ny, nz = self.lattice.shape

        def _mean(block):
            swapped = block.transpose(*pair_swap)
            # q -> -q on the flattened q axis: unflatten, apply the shared
            # FFT-grid map, flatten back (same gather as the historical
            # coordinate-wise index computation).
            swapped_rev = reverse_fft_axes(
                swapped.reshape(nx, ny, nz, *swapped.shape[1:]),
                (0, 1, 2)).reshape(swapped.shape)
            return 0.5 * (block
                          + xp.where(palin, swapped_rev, xp.conj(swapped)))

        cross_sym = _mean(cross_block)
        flip_sym = _mean(spin_flip_block)
        return -cross_sym.transpose(0, 2, 4, 1, 3) + flip_sym

    def _check_transverse_representable(self, ham_orig):
        """Reject input whose transverse vertex is not a function of q alone.

        ``ham_orig`` is in the HAMILTONIAN convention, as built by
        ``_make_ham_inter`` (issue #139).

        Called BEFORE the longitudinal solve, so invalid input fails without
        burning the full solve or leaving a partially populated green_info.

        The criterion is the q-dependence of the ASSEMBLED vertex -- the same
        `ham_pm` the channel will use -- not of the raw blocks. An off-site
        term makes the transverse pair c+_(i a up) c_(j b down) non-local, so
        its vertex cannot be a q-only matrix: measured, the extracted vertex
        for off-site CoulombInter / Ising / Exchange is not proportional to
        V(q) at all (residuals 1.0 / 1.0 / 0.33), while off-site Hund and
        PairLift give no vertex and are therefore fine. Checking the assembled
        vertex also means an internal TABLE whose off-site parts cancel in the
        symmetrised sum is accepted, which is correct: what the channel uses
        is then well-defined. (File input never reaches here unclosed: since
        #93 the readers reject declarations that are not Hermitian-closed.)
        (`_append_pairhop` discards off-site PairHop before this point,
        with a ``logger.warning`` naming the dropped declarations; see
        also the documentation warning.)
        """
        nvol = self.lattice.nvol
        if nvol <= 1:
            return
        xp = _bk.array_module_of(ham_orig)
        ham_pm = self._assemble_transverse_vertex(ham_orig)
        spread = float(_bk.to_host(xp.max(xp.abs(ham_pm - ham_pm[0:1]))))
        scale = float(_bk.to_host(xp.max(xp.abs(ham_pm))))
        # RELATIVE tolerance, with no absolute floor: a floor of max(scale, 1)
        # let a purely q-dependent vertex of small amplitude (< 1e-10) through,
        # and it would then be used by the unsupported calculation. Roundoff
        # from the Fourier transform sits at ~1e-16 relative, so 1e-10 relative
        # separates it cleanly from genuine q-dependence at ANY amplitude.
        # scale == 0 (no interaction reaching these blocks) implies spread == 0
        # and passes.
        if spread > 1e-10 * scale:
            logger.error(
                "calc_type='ring+ladder' requires the transverse vertex to be "
                "independent of q, which fails when a two-body interaction "
                "with a nonzero cross-spin or spin-flip part is off-site: the "
                "transverse pair c+_(i a up) c_(j b down) is then non-local "
                "and its vertex is not a function of q alone. Found "
                "q-dependence {:.3e} on a scale of {:.3e}. Use "
                "calc_type='ring', or restrict those two-body terms to "
                "on-site. (Off-site Hund and PairLift are accepted: their "
                "transverse vertex vanishes.)".format(spread, scale))
            raise ValueError(
                "calc_type='ring+ladder' does not support off-site "
                "cross-spin/spin-flip two-body interactions")

    def _build_transverse_channel(self, chi0q_orig, ham_orig):
        """Build transverse (ladder) channel for RPA.

        The transverse susceptibility chi_+-(q) describes spin-flip
        correlations <S^+(q) S^-(-q)>.

        For paramagnetic systems, the transverse bare susceptibility
        has the same numerical values as the longitudinal chi0:
            chi0_+-[a,c,b,d] = -G_↑[a,b] * G_↓[d,c] = chi0_orb[a,c,b,d]

        The transverse vertex W_+- is obtained by crossing the Hartree
        vertex (Fock exchange):
            W_+-[a,c,b,d] = Gamma_H[(↑,d),(↓,c),(↑,a),(↓,b)]

        Parameters
        ----------
        chi0q_orig : ndarray
            Original bare susceptibility (before spin inflation).
        ham_orig : ndarray
            Interaction tensor in spin-orbital space, in the
            HAMILTONIAN convention as built by _make_ham_inter: the
            transverse assembly re-pairs it itself, so it must NOT be
            pre-converted to the bubble-pair convention (issue #139).

        Returns
        -------
        chi0q_pm : ndarray
            Transverse bare susceptibility, shape matches general scheme.
        ham_pm : ndarray
            Transverse vertex, shape (nvol, norb, norb, norb, norb).
        """
        xp = _bk.array_module_of(chi0q_orig)
        norb = self.norb
        ns = self.ns
        nd = norb * ns
        nvol = self.lattice.nvol

        # --- Build chi0_+- ---
        # For paramagnetic and spin-diagonal cases, chi0_+- = chi0_orb
        # (same numerical values, just interpreted as transverse bubble)
        if self.spin_mode == "spin-free":
            # chi0q_orig shape: (nfreq, nvol, norb, norb, norb, norb).
            # The ladder channel is reachable only with calc_scheme="general"
            # (_set_scheme sys.exit()s otherwise), so enable_reduced is always
            # False here and chi0q is always the full rank-4 orbital tensor.
            if chi0q_orig.ndim != 6:
                # ValueError, consistently with the spinful malformed-rank
                # guard below: both protect the same normally unreachable
                # boundary (a round-3 review found the earlier
                # AssertionError here rested on no real distinction; the
                # spin-diag checks are different -- they reject reachable
                # runtime/provenance conditions, not malformed ranks).
                raise ValueError(
                    "transverse channel expects a general (rank-4 orbital) "
                    "chi0q, got ndim={}. calc_type='ring+ladder' requires "
                    "calc_scheme='general'; if that constraint is ever relaxed, "
                    "a reduced chi0q must be embedded at its DENSITY-PAIR "
                    "positions chi0_+-[a,a,b,b] = chi0_red[a,b] -- NOT scattered "
                    "as delta_{{l2,l4}}, which drops the very components the "
                    "vertex reads (see sc._expand_flex_chi).".format(
                        chi0q_orig.ndim))
            chi0q_pm = chi0q_orig.copy()

        elif self.spin_mode == "spin-diag":
            # Compute exact chi0_+- from G_↑ and G_↓ Green's functions.
            # chi0_+-[a,c,b,d](r,τ) = -G_↑[a,b](r,τ) * G_↓[d,c](-r,-τ)
            if getattr(self, "_chi0q_external", False):
                # A same-instance green0 may belong to an OLDER bubble than
                # the externally supplied chi0q; silently pairing them gives
                # a chiq_pm from different physics (round-3 review).
                raise ValueError(
                    "spin-diag transverse (ladder) channel cannot be "
                    "combined with an externally supplied chi0q: the "
                    "channel needs the Green's functions that produced "
                    "this exact bubble. Recompute chi0q internally.")
            if hasattr(self, 'green0') and self.green0 is not None:
                chi0q_pm_full = self._calc_chi0q_transverse(
                    self.green0, self.green0_tail, 1.0 / self.T)
                # Filter by freq_index if needed
                if len(self.freq_index) < self.nmat:
                    chi0q_pm = chi0q_pm_full[self.freq_index]
                else:
                    chi0q_pm = chi0q_pm_full
                # No reduced-to-general expansion here either: the ladder
                # channel forces calc_scheme="general", so _calc_chi0q_transverse
                # returns the full rank-4 orbital tensor (see the spin-free
                # branch above for what a reduced chi0q would require).
            else:
                # The transverse bubble chi0_+- = -G_up * G_down cannot be
                # reconstructed from the longitudinal chi0q alone (which carries
                # G_up*G_up and G_down*G_down). When the Green's functions are
                # not available -- e.g. chi0q was supplied externally
                # (chi0q_init) on a spin-split system -- using the up-up block
                # would silently give a wrong chi_+-. Fail loudly instead.
                logger.error(
                    "spin-diag transverse (ladder) channel requires the Green's "
                    "functions to build chi0_+- = G_up*G_down, but they are not "
                    "available (chi0q was supplied externally). Recompute chi0q "
                    "internally (do not use chi0q_init) for the ladder channel.")
                raise ValueError(
                    "spin-diag transverse (ladder) channel cannot be computed "
                    "from an externally-supplied chi0q; recompute chi0q "
                    "internally.")

        elif self.spin_mode == "spinful":
            # Phase S (issue #110): the spinful plain-ladder path no
            # longer routes through here at all. It used to slice chi0
            # (the BARE bubble) at the Sz-conserving block and re-solve
            # that reduced slice -- which structurally could not include
            # any cross-spin/spin-mixing vertex content for a genuinely
            # spin-orbit-coupled H0 (the "cross terms are not included"
            # warning this branch used to log). The fix is slice-AFTER-
            # solve, not solve-after-slice: `_solve_impl`'s plain-
            # transverse branch now calls
            # `RPA._extract_transverse_from_dressed` on the ALREADY-
            # DRESSED longitudinal `sol` instead of calling this method
            # for `spin_mode == "spinful"`. Reaching this branch means
            # that routing regressed (a caller invoking
            # `_build_transverse_channel` directly on a spinful solver) --
            # fail loudly rather than silently reproducing the #110 bug.
            raise RuntimeError(
                "_build_transverse_channel: unreachable for "
                "spin_mode='spinful' (issue #110) -- the spinful plain "
                "transverse channel is built by "
                "RPA._extract_transverse_from_dressed, slicing the "
                "DRESSED longitudinal solve, not by slicing chi0 and "
                "re-solving it here. This method must not be called "
                "directly for a spinful solver.")

        # --- Build W_+- (transverse vertex) ---
        # The transverse vertex comes from the CROSS-SPIN block alone. The
        # same-spin block must not appear: a same-spin interaction cannot
        # connect the up and down propagators of the transverse loop, so it
        # produces only self-energy and no vertex at all. That is measured, not
        # argued -- diagonalizing an explicit model with the same-spin part of
        # an interaction switched on, and removing the O(lambda) self-energy
        # exactly (it is the linear response to the Hartree-Fock one-body
        # potential built from the non-interacting density matrix), leaves a
        # vertex of 0 to 1e-6, while the cross-spin part leaves the full value.
        #
        # This used to be `up_block - cross_block.transpose(0, 4, 3, 2, 1)`,
        # which is wrong twice over (issue #90): the `up_block` term should not
        # be there, and the reordering of the cross block was not the right one.
        # For CoulombInter the two errors combined into the antisymmetric
        # remainder V(q) - V(q)^T, which vanishes only when V(q) happens to be
        # symmetric under the orbital-pair transpose -- so one-orbital fixtures
        # never saw it.
        #
        # Measured vertices, by exact diagonalization, all interactions on-site:
        #
        #   CoulombIntra U : -U at (a,a)x(a,a)      table said -U   correct
        #   CoulombInter U': -U' at (a,b)x(b,a)     table said 0    WRONG
        #   Hund J         :  0                     table said -J   WRONG
        #   Ising J        : +J                     table said 2J   WRONG
        #   Exchange J     : -J at (a,a)x(b,b)      table said 0    WRONG
        #   PairLift J     :  0                     table said 0    correct
        #   PairHop J      : -J at (a,b)x(a,b)      table had no entry
        #
        # Every case is q-independent to better than 1e-6, which is the
        # diagnostic that the extraction is clean: an on-site interaction
        # cannot produce a q-dependent vertex. The CoulombInter row is
        # additionally confirmed end to end in chi0's own labels, with a
        # residual that tracks the imaginary-time discretisation (8.9e-3 ->
        # 4.4e-3 -> 2.2e-3 over grids of 192 -> 384 -> 768) while the previous
        # behaviour sits at 1.00.
        #
        # Vertex assembly (and all measurement provenance) lives in
        # _assemble_transverse_vertex; representability was already validated
        # before the longitudinal solve.
        ham_pm = self._assemble_transverse_vertex(ham_orig)

        logger.info("ring+ladder: built transverse channel "
                    "(chi0_pm shape={}, ham_pm shape={})".format(
                        chi0q_pm.shape, ham_pm.shape))

        return chi0q_pm, ham_pm

    def _extract_transverse_from_dressed(self, sol):
        """Extract the physical spinful transverse response by slicing
        the DRESSED (RPA-solved) general spinful tensor -- the Phase-S
        fix for issue #110.

        Background: for a genuinely spin-orbit-coupled H0
        (``spin_mode == "spinful"``), the pre-Phase-S plain-ladder path
        sliced the Sz-conserving block off the BARE bubble ``chi0`` and
        then solved RPA on that reduced slice. That is backwards: the
        RPA equation resums the vertex over whatever chi0 it is given,
        so slicing chi0 FIRST discards every cross-spin/spin-mixing
        vertex contribution before the resummation ever sees it (the
        deleted "cross terms are not included" warning was this defect
        made explicit). The fix is to solve RPA over the FULL
        ``(2*norb)**2`` spin-orbital space first (the antisymmetrized
        general solve `_solve_impl` already performs for the
        longitudinal channel -- `sol = self._solve_rpa(chi0q, ham)`, a
        few lines above this method's call site), and only THEN slice
        out the physical transverse block -- slice-AFTER-solve, not
        solve-after-slice. No second solve is needed: this method is
        pure index selection on the tensor the longitudinal solve
        already produced.

        Equation (spin-block generalized-index convention
        ``g = s*norb + o``, spin-major: ``g < norb`` is spin-up,
        ``g >= norb`` is spin-down -- see `Interaction._make_ham_trans`'s
        own comment, "RPA works internally in spin-block order (index =
        spin*norb + orb)", and this convention's independent
        confirmation in Gate S0's "Operator-phase / index-convention
        check", tests/test_spinful_transverse_ed.py):

            chiq_pm[freq, q, a, c, b, d]
                = sol[freq, q, a, norb+c, b, norb+d]

        i.e. the creation-leg orbital ``a`` reads off the spin-up block,
        the annihilation-leg orbital ``c`` off the spin-down block, and
        likewise ``(b, d)`` for the second (conjugate) leg pair -- the
        SAME index selection the removed pre-Phase-S spinful branch of
        `_build_transverse_channel` applied to chi0 (the bare bubble),
        applied here to ``sol`` (the dressed tensor) instead.

        Pre-adjudicated ahead of this production change (Gate S0, Phase
        S Task 2,
        tests/test_spinful_transverse_ed.py::TestGateS0ExtractionAdjudication)
        against exact diagonalization (`TotalNED`, a total-N-sector ED
        oracle built specifically to host Sz-breaking one-body terms) on
        a genuinely Sz-breaking L=3/norb=2 fixture: PASS at a measured
        ~34,000x tolerance margin (delta_rich=3.30e-07, tol=3.30e-06,
        max_signal=1.11e-01); a deliberate g2/g4 leg-order mutation of
        the same equation FAILs the identical granule at 98% of the
        signal scale, confirming the PASS is not vacuous. This method's
        slice is byte-for-byte the same equation as the test module's
        `extract_pm_from_dressed` reference function (cross-pinned bit
        for bit in Gate S2, Phase S Task 3,
        tests/test_spinful_transverse_ed.py::TestGateS2ProductionRouting).

        Pure index selection: works unchanged on either backend (no
        `float()`/host-only conversion here -- the caller converts to
        host, via `_bk.to_host`, only after this method returns).

        Parameters
        ----------
        sol : ndarray, shape (nfreq, nvol, nd, nd, nd, nd), nd = 2*norb
            The RPA-DRESSED general spinful tensor -- this solver's own
            longitudinal ``sol = self._solve_rpa(chi0q, ham)`` output,
            REUSED verbatim (never recomputed).

        Returns
        -------
        ndarray, shape (nfreq, nvol, norb, norb, norb, norb)
            The full-frequency (dynamic) physical transverse chiq_pm.
        """
        norb = self.norb
        nd = 2 * norb
        if sol.ndim != 6:
            raise ValueError(
                "_extract_transverse_from_dressed: sol must be the "
                "6-dim (nfreq, nvol, nd, nd, nd, nd) dressed general "
                "spinful tensor, got ndim={}".format(sol.ndim))
        if sol.shape[2:] != (nd, nd, nd, nd):
            raise ValueError(
                "_extract_transverse_from_dressed: sol's trailing "
                "orbital axes {} do not match nd=2*norb={} (norb={}, "
                "from the solver) in every slot".format(
                    sol.shape[2:], nd, norb))
        return sol[..., 0:norb, norb:2 * norb, 0:norb, norb:2 * norb]

    @do_profile
    def _solve_rpa(self, chi0q, ham):
        """Solve the RPA equation.

        Parameters
        ----------
        chi0q : ndarray
            Bare susceptibility.
        ham : ndarray
            Interaction Hamiltonian.

        Returns
        -------
        ndarray
            The RPA susceptibility chiq.

        Notes
        -----
        Solves the equation:
        chiq = [1 + chi0q * W]^(-1) * chi0q
        where W is the interaction vertex.

        When the matrices have block-diagonal structure (e.g. spin-diagonal case),
        the solver automatically detects and exploits this to reduce problem size.
        Block structure is cached per ham shape+content to avoid re-detection.

        Frequency parallelization: since ham is frequency-independent, each
        frequency slice can be solved independently. When nmat is large enough,
        the solve is distributed across threads using concurrent.futures.
        """
        logger.debug(">>> RPA._solve_rpa")

        xp = _bk.array_module_of(chi0q)

        nvol = self.lattice.nvol
        nmat = chi0q.shape[0]
        chi_shape = chi0q.shape  # [nmat,nvol,(spin_orbital structure)]
        ndx = int(np.prod(chi_shape[2:2+(len(chi_shape)-2)//2]))

        chi0q_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
        ham_2d = ham.reshape(nvol, ndx, ndx)

        # Detect block structure from the COMBINED sparsity of chi0q AND ham.
        # Block-solving chi = [1 + chi0 ham]^{-1} chi0 is only valid when
        # neither chi0q nor ham couples indices across blocks: if ham is
        # block-diagonal but chi0q has off-block entries (e.g. spin-mixing
        # bands with a spin-diagonal interaction), 1 + chi0 ham acquires
        # off-block entries and the per-block solve is wrong. Detecting from
        # ham alone would miss this, so include chi0q's connectivity.
        # (The sum over (nmat, nvol) is O(nmat*nvol*ndx^2), cheaper than the
        #  O(nmat*nvol*ndx^3) solve. We reduce chi0q to its (ndx, ndx)
        #  connectivity first to avoid materializing a large concatenated
        #  array.)
        # Block detection is pure-python graph analysis on a small (ndx, ndx)
        # connectivity matrix, so bring the reduced connectivity to the host.
        conn_stack = _bk.to_host(xp.stack([
            xp.sum(xp.abs(ham_2d), axis=0),
            xp.sum(xp.abs(chi0q_2d), axis=(0, 1)),
        ]))  # (2, ndx, ndx); _find_block_diagonal sums |.| over axis 0
        blocks = self._find_block_diagonal(conn_stack)

        # Determine thread-parallel chunking for frequency axis
        # LAPACK releases the GIL, so threading gives real parallelism.
        # Only parallelize when there is enough work per thread.
        # LAPACK's batched zgesv already uses internal threads for large batches,
        # so explicit threading only helps for very large problems where the
        # per-chunk overhead is negligible relative to compute.
        # Users can force threading via HWAVE_RPA_THREADS env var.
        import os as _os
        n_workers = int(_os.environ.get("HWAVE_RPA_THREADS", "0"))
        if n_workers == 0:
            import multiprocessing
            n_workers = min(multiprocessing.cpu_count(), 4)
        # Heuristic: only parallelize for large multi-orbital problems
        # where ndx >= 16 (8+ orbitals spinful) and enough frequency points.
        # For small ndx, batched LAPACK is faster than thread pool overhead.
        # cuSOLVER's batched solve already saturates the GPU, so the CPU
        # thread-pool path is numpy-only.
        use_parallel = (xp is np and n_workers > 1 and nmat >= 4 * n_workers
                        and ndx >= 16 and nmat * nvol * ndx * ndx >= 1000000)

        if blocks is not None and len(blocks) > 1:
            logger.info("_solve_rpa: block-diagonal structure detected, "
                        "ndx={} -> {} blocks of sizes {}".format(
                            ndx, len(blocks), [len(b) for b in blocks]))
            sol = xp.zeros_like(chi0q_2d)
            for block_idx in blocks:
                idx = xp.array(block_idx)
                ix = xp.ix_(idx, idx)
                chi0q_blk = chi0q_2d[:, :, ix[0], ix[1]]
                ham_blk = ham_2d[:, ix[0], ix[1]]
                nb = len(block_idx)

                if use_parallel:
                    sol[:, :, ix[0], ix[1]] = self._solve_rpa_parallel(
                        chi0q_blk, ham_blk, nb, n_workers)
                else:
                    mat_blk = (xp.eye(nb, dtype=np.complex128)
                               + (chi0q_blk @ ham_blk[np.newaxis, :, :, :]))
                    sol[:, :, ix[0], ix[1]] = xp.linalg.solve(mat_blk, chi0q_blk)
        else:
            if use_parallel:
                sol = self._solve_rpa_parallel(chi0q_2d, ham_2d, ndx, n_workers)
            else:
                mat = (xp.eye(ndx, dtype=np.complex128)
                       + (chi0q_2d @ ham_2d[np.newaxis, :, :, :]))
                sol = xp.linalg.solve(mat, chi0q_2d)

        return sol.reshape(chi_shape)

    @staticmethod
    def _solve_rpa_parallel(chi0q_2d, ham_2d, ndx, n_workers):
        """Solve RPA equation with frequency-axis thread parallelism.

        Parameters
        ----------
        chi0q_2d : ndarray, shape (nmat, nvol, ndx, ndx)
            Bare susceptibility in 2D matrix form.
        ham_2d : ndarray, shape (nvol, ndx, ndx)
            Frequency-independent interaction Hamiltonian.
        ndx : int
            Matrix dimension.
        n_workers : int
            Number of threads.

        Returns
        -------
        sol : ndarray, shape (nmat, nvol, ndx, ndx)
            RPA susceptibility.
        """
        from concurrent.futures import ThreadPoolExecutor
        nmat = chi0q_2d.shape[0]
        sol = np.empty_like(chi0q_2d)
        eye = np.eye(ndx, dtype=np.complex128)
        ham_bc = ham_2d[np.newaxis, :, :, :]  # broadcast-ready

        # Split frequency axis into contiguous chunks
        chunk_size = (nmat + n_workers - 1) // n_workers
        slices = [slice(i, min(i + chunk_size, nmat))
                  for i in range(0, nmat, chunk_size)]

        def _solve_chunk(sl):
            chi0q_chunk = chi0q_2d[sl]
            mat = eye + (chi0q_chunk @ ham_bc)
            sol[sl] = np.linalg.solve(mat, chi0q_chunk)

        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            list(pool.map(_solve_chunk, slices))

        return sol

    def _find_block_diagonal(self, ham_2d):
        """Detect block-diagonal structure from the interaction Hamiltonian.

        Parameters
        ----------
        ham_2d : ndarray, shape (nvol, ndx, ndx)
            Interaction Hamiltonian reshaped to 2D matrices.

        Returns
        -------
        list of list of int, or None
            List of index groups forming independent blocks.
            Returns None if no block structure is found (single block).
        """
        ndx = ham_2d.shape[-1]
        if ndx <= 1:
            return None

        # Sum absolute values over nvol to get connectivity pattern
        connectivity = np.sum(np.abs(ham_2d), axis=0)

        # Build adjacency from non-zero off-diagonal entries
        threshold = 1.0e-12
        adj = (np.abs(connectivity) > threshold) | (np.abs(connectivity.T) > threshold)

        # Union-Find with path halving for connected components
        parent = list(range(ndx))

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]  # path halving
                x = parent[x]
            return x

        rows, cols = np.where(np.triu(adj, k=1))
        for i, j in zip(rows, cols):
            ri, rj = find(i), find(j)
            if ri != rj:
                if ri > rj:
                    ri, rj = rj, ri
                parent[rj] = ri

        # Collect components
        components = {}
        for i in range(ndx):
            r = find(i)
            if r not in components:
                components[r] = []
            components[r].append(i)

        if len(components) <= 1:
            return None

        return list(components.values())

    # -------------------------------------------------------------------
    # Phase W (bond-resolved transverse channel): the internal pipeline
    # entry (spec 2026-08-15-bond-transverse-design.md, "Phase W --
    # implementation + ED campaign"; plan Task 5). Appended for this task
    # only -- the methods above are unmodified.
    # -------------------------------------------------------------------

    def _build_ham_pm_onsite(self):
        """Build the on-site-only transverse vertex that feeds channel 0
        of ``W_pm_bond`` (spec "The dressing (static)", step 1).

        Locality is judged on the ORIGINAL PRE-FOLD declarations
        (``self.ham_info.param_ham_orig`` when the lattice has a
        sublattice, else ``self.ham_info.param_ham``), filtered to
        ``irvec == (0, 0, 0)`` -- the SAME pattern ``_make_ham_inter``'s
        ``onsite_r``/``_append_onsite_direct`` closures use for the
        spinful Fierz-exchange correction (issue #137), generalized here
        across every interaction type and independent of
        ``enable_spin_orbital``. Reading locality off the FOLDED table
        (``self.ham_info.param_ham`` under sublattice folding) would fold
        an off-site bond onto ``r=(0,0,0)`` between supercell orbitals
        and silently smuggle off-site content into the on-site channel --
        the pre-fold locality trap ``_append_inter_cross``/
        ``_append_onsite_direct`` document and PR history has hit
        repeatedly.

        Each filtered table is symmetrised with the SAME per-slot-family
        rule ``_make_ham_inter``'s ``_symmetrised`` closure applies
        (plain mean for the density/flip types, Hermitian mean for
        PairHop): restricted to a single (on-site) cell, the two
        closures are algebraically identical, because
        ``reverse_fft_axes`` "leaves index 0 in place" for any mesh
        shape -- the reversal partner of an ``r=(0,0,0)`` entry is always
        the SAME cell with orbitals swapped.

        The resulting on-site tensor is q-INDEPENDENT (an on-site
        interaction cannot produce a q-dependent vertex --
        ``_assemble_transverse_vertex``'s own docstring measurement), so
        it is broadcast identically across ``nvol`` before assembly.

        Returns
        -------
        ndarray, complex128, shape (nvol, norb, norb, norb, norb)
            ``_assemble_transverse_vertex`` applied to the on-site-only
            interaction tensor.
        """
        ham_info = self.ham_info
        lattice = self.lattice
        nvol = lattice.nvol
        norb = self.norb
        ns = self.ns
        nd = norb * ns
        has_sub = getattr(lattice, "has_sublattice", False)
        source = ham_info.param_ham_orig if has_sub else ham_info.param_ham

        spin_table_cache = {
            t: ring_spin_table(t)
            for t in ('CoulombIntra', 'CoulombInter', 'Hund', 'Ising',
                      'PairLift', 'Exchange', 'PairHop')
        }

        def _onsite_table(tbl):
            # Filter a PRE-FOLD table to irvec == (0, 0, 0), fold
            # orbital indices under sublattice (never spatial content --
            # an on-site declaration stays on-site under folding, but the
            # re-filter below stays defensive, mirroring
            # _append_onsite_direct exactly).
            filtered = {}
            for (irvec, orbvec), v in (tbl or {}).items():
                if tuple(int(x) for x in irvec) == (0, 0, 0):
                    filtered[((0, 0, 0),
                              (int(orbvec[0]), int(orbvec[1])))] = v
            if not filtered:
                return {}
            if has_sub:
                filtered = ham_info._reshape_interaction(filtered, False)
                filtered = {
                    (irvec, orbvec): v
                    for (irvec, orbvec), v in filtered.items()
                    if tuple(irvec) == (0, 0, 0)
                }
            return filtered

        def _symmetrised_onsite(type_name, tbl):
            # Algebraic equivalent of _make_ham_inter's _symmetrised
            # closure, restricted to the on-site cell (see docstring).
            arr = np.zeros((norb, norb), dtype=np.complex128)
            for (_irvec, orbvec), v in tbl.items():
                arr[orbvec] += v
            if type_name == "PairHop":
                sym = 0.5 * (arr + np.conjugate(arr.T))
            else:
                sym = 0.5 * (arr + arr.T)
            out = {}
            for a, b in zip(*np.nonzero(sym)):
                out[((0, 0, 0), (int(a), int(b)))] = sym[a, b]
            return out

        ham_r8 = np.zeros((*(ns, norb) * 4,), dtype=np.complex128)

        def _append_density_like(type_name, tbl=None):
            src_tbl = tbl if tbl is not None else source.get(type_name)
            sym = _symmetrised_onsite(type_name, _onsite_table(src_tbl))
            for (_irvec, orbvec), v in sym.items():
                a, b = orbvec
                for spinvec, w in spin_table_cache[type_name].items():
                    s1, s2, s3, s4 = spinvec
                    # beta beta' alpha alpha' -- same slot layout as
                    # _make_ham_inter's _append_inter.
                    orb = (s4, b, s3, b, s1, a, s2, a)
                    ham_r8[orb] += v * w

        def _append_pairhop_like(type_name, tbl=None):
            src_tbl = tbl if tbl is not None else source.get(type_name)
            sym = _symmetrised_onsite(type_name, _onsite_table(src_tbl))
            for (_irvec, orbvec), v in sym.items():
                a, b = orbvec
                for spinvec, w in spin_table_cache[type_name].items():
                    s1, s2, s3, s4 = spinvec
                    # same slot layout as _make_ham_inter's
                    # _append_pairhop.
                    orb = (s4, b, s3, a, s1, a, s2, b)
                    ham_r8[orb] += v * w

        if 'Coulomb' in source:
            # Aggregate 'Coulomb' input splits into intra/inter parts on
            # the PRE-FOLD table, mirroring _make_ham_inter's own
            # _append_onsite_direct call sites for the 'Coulomb' key.
            coulomb_intra_pre, coulomb_inter_pre = wan90.split_coulomb(
                source['Coulomb'])
            _append_density_like('CoulombIntra', tbl=coulomb_intra_pre)
            _append_density_like('CoulombInter', tbl=coulomb_inter_pre)
        if 'CoulombIntra' in source:
            _append_density_like('CoulombIntra')
        if 'CoulombInter' in source:
            _append_density_like('CoulombInter')
        if 'Hund' in source:
            _append_density_like('Hund')
        if 'Ising' in source:
            _append_density_like('Ising')
        if 'PairLift' in source:
            _append_density_like('PairLift')
        if 'Exchange' in source:
            _append_density_like('Exchange')
        if 'PairHop' in source:
            _append_pairhop_like('PairHop')

        ham_r4 = ham_r8.reshape(nd, nd, nd, nd)
        ham_onsite_nvol = np.broadcast_to(
            ham_r4[None, :, :, :, :], (nvol, nd, nd, nd, nd)).copy()
        return self._assemble_transverse_vertex(ham_onsite_nvol)

    def _run_transverse_bond_pipeline(self, green_kw, green0_tail, beta,
                                       topo):
        """The bond-resolved transverse (Phase W) internal pipeline
        entry: on-site vertex assembly -> ``W_pm_bond`` ->
        ``transverse_bond_bubble_static`` -> dressing (the EXISTING
        ``_solve_rpa``'s ``I + chi0.ham`` convention) -> collapse to the
        plain ``(nvol, norb, norb, norb, norb)`` shape (spec "The
        dressing (static)" / "Collapse rule (exact)").

        Test-invocable, NO prereq fallback (anti-vacuity, spec "Phase W
        -- implementation + ED campaign"): every granule and gate W1 call
        this entry directly; the production gate's own fallback behavior
        (Task 7) is validated separately by its own config matrix.

        Parameters
        ----------
        green_kw : ndarray, complex, shape (2, nmat, nvol, norb, norb)
            Canonical two-block Green's function (block 0 = up, block 1
            = down) -- same contract as ``_calc_chi0q_transverse``.
        green0_tail : ndarray or None
            Paired tau-space tail add-back, forwarded to
            ``bubble.transverse_bond_bubble_static`` (``None`` disables
            the tail correction).
        beta : float
            Inverse temperature.
        topo : bond_channels.TransverseTopology
            The master transverse bond topology (e.g.
            ``bond_channels.resolve_transverse_topology`` or an
            equivalent hand-built, invariant-satisfying instance); its
            off-site channels feed ``W_pm_bond``'s cross/flip families
            and ``chi0``'s bond-pair structure. On-site (channel 0)
            content is NOT read from ``topo`` -- it is built
            independently by ``_build_ham_pm_onsite`` from the solver's
            OWN pre-fold declarations (spec step 1); ``W_pm_bond`` places
            it verbatim into channel 0.

        Returns
        -------
        chi_pm_bond_static : ndarray, complex128, shape (nvol, ND, ND)
            The dressed bond-transverse static susceptibility, ``ND =
            B * norb**2``.
        chiq_pm_static : ndarray, complex128, shape
            (nvol, norb, norb, norb, norb)
            The collapsed (m=0, m'=0) sub-block (spec "Collapse rule
            (exact)": ``chiq_pm_static[q, a, c, b, d] =
            chi_pm_bond[q, (0, a, c), (0, b, d)]``) -- the same object
            gate W1 compares against today's plain ``chiq_pm``.
        """
        spatial_shape = tuple(int(x) for x in self.lattice.shape)
        nvol = self.lattice.nvol
        workers = getattr(self, "fft_workers", 1)
        norb = self.norb

        ham_pm_onsite = self._build_ham_pm_onsite()

        W = bond_channels.W_pm_bond(
            topo, ham_pm_onsite, spatial_shape=spatial_shape)
        chi0 = bubble.transverse_bond_bubble_static(
            green_kw, green0_tail, beta, topo,
            spatial_shape=spatial_shape, workers=workers)

        ND = W.shape[-1]
        # W_pm_bond assembles on the host (numpy) by design; under GPU
        # execution chi0 lives on the device, and _solve_rpa dispatches on
        # chi0's backend -- move W there first (numpy.asarray is a no-op
        # on the CPU path).
        W = _bk.array_module_of(chi0).asarray(W)
        # The dressing (spec "The dressing (static)"): the EXISTING
        # _solve_rpa's I + chi0 @ ham convention, with the leading
        # length-1 axis satisfying its frequency-axis interface. No sign
        # flip anywhere -- W's channel-0 block IS ham_pm_onsite verbatim.
        sol = self._solve_rpa(chi0.reshape(1, nvol, ND, ND), W)
        chi_pm_bond_static = sol[0]

        # Collapse (spec "Collapse rule (exact)"): the elementwise
        # (m=0, m'=0) sub-block. Channel 0 occupies indices [0, nd) on
        # both axes (W_pm_bond step 1), and the local orbital-pair index
        # within that block is (a, c) -> a*norb+c (the SAME bond-major/
        # orbital-minor convention transverse_bond_bubble_static and
        # W_pm_bond's channel-zero embedding use) -- so a bare slice +
        # C-order reshape recovers (q, a, c, b, d) directly, with no
        # further transpose.
        nd = norb * norb
        chiq_pm_static = chi_pm_bond_static[:, :nd, :nd].reshape(
            nvol, norb, norb, norb, norb)

        return chi_pm_bond_static, chiq_pm_static

    # -------------------------------------------------------------------
    # Phase W Task 7: the experimental gate's own prereq validation and
    # resource preflight. Both are called from solve() ONLY when
    # transverse_bond_channels=true (calc_type='ring+ladder' already
    # required to reach either), mirroring where sc.py's
    # _validate_bond_prereqs / _bond_resource_preflight run for the
    # (eliashberg) longitudinal bond gate: after the prerequisite runtime
    # state (here, self.spin_mode) is known, and BEFORE the expensive
    # solve, so an invalid/oversized request fails fast.
    # -------------------------------------------------------------------

    def _validate_transverse_bond_prereqs(self):
        """Top-level guards for ``transverse_bond_channels=true`` (spec
        "Production surface"): no active off-site transverse-resolved
        declaration, ``spin_mode='spinful'``, and an externally supplied
        ``chi0q_init`` are each a REJECT, naming the gate as inapplicable
        rather than silently falling back to the plain path -- every
        gate-on run that starts is bond-owned (spec: "the rev-3
        warn-and-fall-back is RETRACTED").

        Must be called AFTER ``self.spin_mode`` is known (i.e. from
        inside ``solve()``, at the same point ``_check_transverse_
        representable`` would otherwise run) and BEFORE the expensive
        solve.

        Returns
        -------
        bond_channels.TransverseTopology
            The resolved master transverse topology, computed once here
            so the activity predicate and the actual pipeline call see
            the identical object (never re-resolved, and therefore never
            able to silently disagree with what was just validated).
        """
        if self.spin_mode == "spinful":
            raise ValueError(
                "[mode.param] transverse_bond_channels=true cannot be "
                "combined with spin_mode='spinful' (a genuinely "
                "spin-orbit-coupled H0): the bond-resolved dressing is "
                "implemented for the spin-diagonal transverse channel "
                "only. The PLAIN (non-bond) spinful transverse channel "
                "is fully supported -- set transverse_bond_channels=false "
                "to use it. Only the composition of the bond-resolved "
                "path with the full spinful space is deferred and not "
                "yet implemented. Use a spin-diagonal system (no "
                "spin-mixing transfer/enable_spin_orbital term) if "
                "transverse_bond_channels=true is required.")

        if getattr(self, "_chi0q_external", False):
            raise ValueError(
                "[mode.param] transverse_bond_channels=true cannot be "
                "combined with an externally supplied chi0q (chi0q_init): "
                "the bond path needs the Green's functions that produced "
                "the bubble (chi0_pm_bond is built from G_up/G_down "
                "directly), which a file-based chi0q does not carry. "
                "Recompute chi0q internally (omit chi0q_init) or set "
                "transverse_bond_channels=false.")

        if getattr(self.lattice, "has_sublattice", False):
            raise ValueError(
                "[mode.param] transverse_bond_channels=true is not "
                "supported with a sublattice (SubShape < CellShape): the "
                "bond-resolved transverse channel does not yet implement "
                "the sublattice folding map for its off-site channels. "
                "Run without a sublattice (SubShape == CellShape) or set "
                "transverse_bond_channels=false.")

        has_sub = getattr(self.lattice, "has_sublattice", False)
        interactions = (self.ham_info.param_ham_orig if has_sub
                        else self.ham_info.param_ham)

        # The aggregate 'Coulomb' table (wan90.split_coulomb's shared
        # decomposition -- the SAME split _make_ham_inter and
        # _build_ham_pm_onsite each apply to feed the Hamiltonian and the
        # on-site vertex) is invisible to resolve_transverse_topology,
        # which only ever reads the 'CoulombInter'/'Ising'/'Exchange'
        # keys: without this, an off-site inter-orbital term declared
        # ONLY via aggregate 'Coulomb' would be seen everywhere else in
        # the pipeline but not by this gate's topology (wrongly reporting
        # "no active declaration", or silently computing without that
        # channel). Merge its inter-orbital off-site part into the
        # 'CoulombInter' table fed to the resolver. The ambiguity guard
        # ('Coulomb' cannot coexist with an explicit 'CoulombIntra'/
        # 'CoulombInter' declaration) already ran at construction time
        # (Interaction.__init__ -> _make_ham_inter, long before solve()
        # -- and therefore this method -- is ever reached), so
        # 'CoulombInter' is always empty here whenever 'Coulomb' is
        # present; the merge below is still written generally (update,
        # not overwrite) so it stays correct even if that guard's scope
        # ever changes.
        #
        # ``interactions`` (``self.ham_info.param_ham``/``param_ham_orig``)
        # is a ``CaseInsensitiveDict`` (see ``Interaction._init_interaction``'s
        # own comment on this exact defect class): callers may declare
        # 'ising'/'exchange'/'coulomb' in any case and ``_make_ham_inter``
        # still honors it. A plain ``dict(interactions)`` copy would
        # DOWNGRADE that -- it preserves each key's stored case as an
        # ordinary (case-SENSITIVE) dict key, so a lowercase 'ising' would
        # silently stop matching ``resolve_transverse_topology``'s
        # canonical-cased ``interactions.get('Ising', {})`` lookup and
        # vanish from the topology. Copy into a fresh
        # ``CaseInsensitiveDict`` instead, which preserves case-insensitive
        # lookup for every key untouched by this merge.
        if 'Coulomb' in interactions:
            _coulomb_intra, coulomb_inter_agg = wan90.split_coulomb(
                interactions['Coulomb'])
            interactions = CaseInsensitiveDict(interactions)
            merged_inter = dict(interactions.get('CoulombInter', {}) or {})
            merged_inter.update(coulomb_inter_agg)
            interactions['CoulombInter'] = merged_inter

        topo = bond_channels.resolve_transverse_topology(
            interactions, np.eye(3), self.norb,
            max_shells=self.transverse_bond_max_shells)

        # "Active" is defined ONCE, operationally, by the shared
        # predicate (spec): the same function the topology resolver's
        # own truncation guard implicitly relies on. A topology that
        # resolves to nothing off-site (nothing declared, or everything
        # truncated away -- resolve_transverse_topology itself already
        # refuses a truncation that would drop declared nonzero content,
        # so "truncated to nothing" here only ever means "nothing was
        # declared to begin with") has no transverse bond physics to
        # represent.
        if not bond_channels.transverse_effective_activity(topo):
            raise ValueError(
                "[mode.param] transverse_bond_channels=true requires an "
                "active off-site transverse-resolved declaration "
                "(CoulombInter, Ising or Exchange -- or the off-site part "
                "of an aggregate Coulomb declaration -- with a nonzero "
                "off-site coefficient after duplicate summation, Hermitian "
                "projection and transverse_bond_max_shells truncation); "
                "none is present, so the bond-resolved gate has nothing "
                "to represent. Declare an off-site term, relax "
                "transverse_bond_max_shells, or set "
                "transverse_bond_channels=false.")

        return topo

    def _transverse_bond_resource_preflight(self, topo):
        """Estimate the peak memory/compute cost of the bond-resolved
        transverse gate and refuse to exceed
        ``transverse_bond_memory_cap_gb`` (spec "Phase W -- Budget
        (stated, not deferred)").

        The pipeline (``_run_transverse_bond_pipeline``, rpa.py) has two
        phases that do NOT overlap in time -- the bubble phase
        (``bubble.transverse_bond_bubble_static``, called first) returns
        its result before the dressing phase (``_solve_rpa``, called
        second) allocates anything, so the two phases' peaks are taken as
        a ``max(...)``, not summed: whichever phase is more expensive for
        the given shapes bounds the run.

        Byte estimate, phase by phase (``itemsize = 16`` for complex128,
        ``P = norb**2`` the single-channel orbital-pair count, ``ND = B *
        P`` the bond-enlarged dimension, ``Nq = self.lattice.nvol``,
        ``Nmat = self.nmat`` -- available here because ``_init_param``
        (which sets ``self.nmat``, see ``__init__``'s ordered
        ``_init_mode -> _init_param -> ... -> _init_interaction``
        sequence) has already run by the time ``solve()`` reaches this
        preflight call, well before the (still-unrun) longitudinal
        solve):

        - **Dressing (solve) phase**: persistent storage for
          ``chi0_pm_bond`` + ``W_pm_bond`` + the dressed output + the
          ``_solve_rpa`` solve workspace, ``solve_bytes = (3 + K_solve) *
          Nq * ND**2 * itemsize`` with ``K_solve = 3`` -- a documented
          CONSERVATIVE POLICY ALLOWANCE for the numpy backend (LU factor
          copy + pivots + solve output; the ``_solve_rpa`` block-detection
          temporary is counted inside it), NOT a guaranteed LAPACK bound.
          This term has NO ``Nmat`` dependence: the solve only ever sees
          the static (``Omega=0``) slice.

        - **Bubble phase** (previously OMITTED entirely -- this is the
          fix; a first round undercounted its pair-processing sub-phase
          and never bounded its preparation sub-phase -- this is the
          correction). ``_run_transverse_bond_pipeline`` builds ``W``
          (via ``bond_channels.W_pm_bond``, shape ``(Nq, ND, ND)``)
          BEFORE calling ``bubble.transverse_bond_bubble_static``, and
          that call holds ``W`` resident throughout (it is consumed only
          afterwards, by the dressing phase). The bubble call itself has
          two sub-phases with DIFFERENT resident sets -- ``W`` is the
          only ``ND``-scale tensor alive during preparation (``chi_bar``
          below does not exist yet), so the two sub-phases' peaks are
          bounded SEPARATELY and combined with ``max(...)`` (never
          summed, for the same non-overlap reason as the two top-level
          phases):

          *Preparation sub-phase* (``_prepare_dense``, called once on
          the full ``nblock=2`` up/down Green tensor, BEFORE the
          per-block split): let ``U = 2 * Nmat * Nq * P * itemsize`` be
          the byte size of one full two-block ``(2, Nmat, Nq, norb,
          norb)`` tensor. Every internal step that is not a pure
          ``reshape`` (view) allocates a NEW ``U``-sized buffer while its
          predecessor is still referenced -- ``matsubara.fermion_to_tau``
          internally computes ``xp.fft.fft(arr, axis=1) *
          _bcast(omg, ...)``, an ``fft`` output buffer times a
          RESHAPED (not tiled) phase, so exactly one extra ``U``-buffer
          coexists with the fft output during that multiply;
          ``backend.spatial_ifftn`` (new buffer) runs while its input is
          still bound; ``kgrid.reverse_fft_axes`` (its own docstring:
          "New array") runs while ``green_rt`` is still bound because it
          is returned alongside ``green_rev`` in ``prepped``. Each of
          these is a 2-buffer transient, so the tail-OFF peak is
          ``2*U = 4*Nmat*Nq*P*itemsize``. Tail-ON additionally builds
          ``jump_f``/``jump_r_rev``/``fwd0_p``/``rev0_p``, four buffers
          of size ``S = 2*Nq*P*itemsize`` each (a single-frequency slice
          of the two-block tensor -- no ``Nmat`` factor) that all persist
          to the return, on top of the still-resident ``green_rt``/
          ``green_rev`` (``2*U``): tail-ON peak = ``2*U + 4*S =
          4*Nmat*Nq*P*itemsize + 8*Nq*P*itemsize``. Tail-ON is the
          uniformly larger (and therefore used) bound -- whether the
          tail is actually active depends on runtime Green-function
          content, not on anything known at preflight time:
          ``prep_bytes = itemsize * Nq * (ND**2 + 4*P*(Nmat + 2))``
          (the ``ND**2`` term is ``W`` alone).

          *Pair-processing sub-phase* (the ``B x B`` channel-pair loop):
          ``_iter_transverse_bond_channel_pairs`` extracts and holds, for
          the ENTIRE loop, two full-frequency Green carriers
          ``green_fwd_sgn`` and ``green_rev``, each shape ``(Nmat, Nq,
          norb, norb)`` (one block out of the ``nblock=2`` pair, the
          other half released immediately per the Task-4 per-block
          ``del`` discipline) -- contributing ``2 * Nmat * Nq * P *
          itemsize`` (``P = norb**2``) alongside the already-resident
          ``W`` AND, from this point on, ``chi_bar`` (the ``(Nq, ND,
          ND)`` static accumulator, allocated right before the loop
          starts): ``2 * Nq * ND**2 * itemsize``. Per pair,
          ``_bond_pair_full_block`` materializes the FULL-frequency block
          ``(Nmat, Nq, norb, norb, norb, norb)`` == ``(Nmat, Nq, P, P)``
          (``contract_general``'s outer-product buffer, immediately
          reinterpreted as the ``(npair, npair)`` pair block, ``npair =
          P``). Because that buffer is a non-contiguous ``swapaxes``
          view, the subsequent ``.reshape`` forces a fresh contiguous
          copy, then ``spatial_fftn`` allocates its own output -- each of
          these is an ordinary 2-buffer (old + new) transient. The THIRD
          step, ``matsubara.tau_to_boson``, is NOT: it computes
          ``xp.fft.ifft(arr * _bcast(omg_inv, ...), axis=0)`` where
          ``arr`` is the caller's ``chi0_qt`` (still bound in
          ``_bond_pair_full_block`` -- not ``del``'d until AFTER
          ``tau_to_boson`` returns). ``arr * _bcast(...)`` allocates ONE
          new buffer (``_bcast`` reshapes the phase vector, it does not
          tile it, so the multiply is an ordinary broadcast producing a
          single output), and ``xp.fft.ifft`` on THAT allocates a SECOND
          new buffer while the first is still on the evaluation stack --
          so at the peak, THREE full ``(Nmat, Nq, P, P)`` buffers coexist
          simultaneously: ``arr`` (== ``chi0_qt``, the caller's still-live
          reference), the multiply's output, and the ``ifft`` output (a
          first review round only counted two, missing that ``chi0_qt``
          itself remains alive across the call). The tail-correction
          block (``chi0_rt[0] = ...``, when active) runs BEFORE this
          chain, while only ONE full ``(Nmat, Nq, P, P)`` buffer
          (``chi0_rt`` itself) is resident, and adds a bounded ``O(P**2)``
          (no ``Nmat`` factor -- its own operands are single-frequency
          slices): its own peak, ``(Nmat + 3) * P**2``, is algebraically
          ``<= 3 * Nmat * P**2`` for every valid (even, positive)
          ``Nmat >= 2``, so it never exceeds the THREE-buffer peak.

          RETAINED TAIL CARRIERS (a second review round's fix -- these
          were bounded inside ``prep_bytes`` but omitted from
          ``pair_bytes``): ``prepped`` (the ``_DensePrepared`` instance
          ``_prepare_dense`` returns) is held alive by
          ``transverse_bond_bubble_static``'s own local variable for the
          ENTIRE function, including the whole pair-processing loop --
          not just during preparation. ``_iter_transverse_bond_channel_
          pairs`` only clears ``prepped.green_rt``/``prepped.green_rev``
          (setting them ``None`` right after extracting the single-block
          carriers above); it never touches ``prepped.fwd0_p``/
          ``rev0_p``/``jump_f``/``jump_r_rev``. When tail-on, those are
          the SAME four full two-block ``S = 2*Nq*P*itemsize`` buffers
          ``prep_bytes`` already counts -- ``tail_pack`` (built once,
          also for the whole loop) holds only single-block VIEWS into
          them, so the underlying two-block buffers are what actually
          stay resident, not merely their views. This is genuinely
          additional storage throughout the pair loop, on top of the
          THREE-buffer per-pair transient: ``+ 4*S = 8*Nq*P*itemsize``.
          Whether these four arrays exist at all is decided by
          ``_prepare_dense``'s own ``tail_on = bool(xp.any(green0_tail
          ...[:, 0] != 0))`` check on actual Green-function content --
          ``self.green0_tail`` is available at preflight time (set by
          ``solve()`` before this method runs), but RPA's own
          ``_calc_green`` NEVER returns ``None`` for it (its docstring:
          "materializes a shape-matched all-zero tail" when
          ``coeff_tail == 0``), so an ``is not None`` presence check
          would always be ``True`` here and provide no tightening --
          replicating ``xp.any(...)`` itself would duplicate a
          host/device-synchronizing reduction ``_prepare_dense`` already
          performs once, which is not "cheap" to redo speculatively at
          preflight time. So, like ``prep_bytes``, this term is counted
          UNCONDITIONALLY (the conservative bound, per this function's
          existing policy):
          ``pair_bytes = itemsize * Nq * (2*ND**2 + 3*Nmat*P**2 +
          2*Nmat*P + 8*P)``.

          CARRIER-EXTRACTION TRANSITION (verified, not merely assumed):
          at the instant ``_iter_transverse_bond_channel_pairs`` begins
          -- ``prepped.green_rt`` (2-block, ``2*Nmat*P``) and
          ``prepped.green_rev`` (2-block, ``2*Nmat*P``) both still whole,
          plus the newly-forming single-block ``green_fwd_sgn``
          (``Nmat*P``), plus the same retained tail carriers (``8*P``),
          plus ``W``+``chi_bar`` (``2*ND**2``) -- the total is
          ``2*ND**2 + 5*Nmat*P + 8*P``. Comparing to the steady-state
          ``pair_bytes`` above (``2*ND**2 + 3*Nmat*P**2 + 2*Nmat*P +
          8*P``), the difference is ``3*Nmat*P**2 - 3*Nmat*P =
          3*Nmat*P*(P - 1)``, which is exactly ``0`` at ``norb=1``
          (``P=1``) and strictly positive for ``norb>1``: the
          steady-state ``pair_bytes`` formula bounds the extraction
          transition too, with equality (not slack) at the ``norb=1``
          corner -- so no separate term is needed, but the bound is
          tight there, not merely "large enough by coincidence".

          ``bubble_bytes = max(prep_bytes, pair_bytes)``. With the
          retained-tail-carrier term now in BOTH, ``pair_bytes -
          prep_bytes = ND**2 + Nmat*P*(3*P - 2)``, which is strictly
          positive for every valid ``P >= 1``, ``Nmat >= 2``, ``ND >= 1``
          -- so ``pair_bytes`` now provably dominates ``prep_bytes``
          everywhere (unlike before this fix, where ``prep_bytes`` won
          at the smallest corner, ``norb=1``/``B=1``/``Nmat=2``: ``17``
          vs the OLD ``pair_bytes``'s ``12``, in units of ``itemsize *
          Nq`` -- with the fix, the same corner gives ``17`` vs the NEW
          ``pair_bytes``'s ``20``). The ``max`` is still computed
          explicitly in code rather than simplified away, for robustness
          against future changes to either sub-phase.

        ``peak_bytes = max(bubble_bytes, solve_bytes)``. For small ``B``
        (few bond channels, so ``ND`` and therefore ``solve_bytes`` stay
        small) and large ``Nmat``, ``bubble_bytes`` -- which has no
        counterpart bounding it in the old solve-only estimate -- can
        exceed ``solve_bytes`` arbitrarily; this is exactly the case the
        old estimate missed and could OOM on after passing the cap. Named
        an ESTIMATE in every user-facing message, per spec.

        Also emits a WARNING (never a refusal) when the estimated solve
        op-count ``Nq * ND**3`` exceeds ``1e12``: a memory-fitting run can
        still be compute-prohibitive.
        """
        B = int(len(topo.delta_r))
        norb = self.norb
        P = norb ** 2
        ND = B * P
        Nq = int(self.lattice.nvol)
        Nmat = int(self.nmat)
        K_solve = 3
        itemsize = 16  # complex128

        solve_bytes = (3 + K_solve) * Nq * (ND ** 2) * itemsize
        prep_bytes = itemsize * Nq * ((ND ** 2) + 4 * P * (Nmat + 2))
        pair_bytes = itemsize * Nq * (
            2 * (ND ** 2) + 3 * Nmat * (P ** 2) + 2 * Nmat * P + 8 * P)
        bubble_bytes = max(prep_bytes, pair_bytes)
        peak_bytes = max(bubble_bytes, solve_bytes)
        op_count = Nq * (ND ** 3)

        logger.info(
            "Transverse bond-channel preflight (ESTIMATE): B = %d "
            "channels, ND = B*norb**2 = %d, N_q = %d, Nmat = %d, "
            "K_solve = %d (numpy-backend conservative allowance), "
            "estimated bubble-prep-phase peak = %.3f GB, estimated "
            "bubble-pair-phase peak = %.3f GB, estimated solve-phase "
            "peak = %.3f GB, estimated overall peak memory = %.3f GB "
            "(cap %.3f GB), estimated solve op-count Nq*ND**3 = %.3e",
            B, ND, Nq, Nmat, K_solve, prep_bytes / 1.0e9,
            pair_bytes / 1.0e9, solve_bytes / 1.0e9, peak_bytes / 1.0e9,
            self.transverse_bond_memory_cap_gb, op_count)

        if peak_bytes > self.transverse_bond_memory_cap_gb * 1.0e9:
            raise ValueError(
                "[mode.param] transverse_bond_channels: ESTIMATED peak "
                "memory {:.3f} GB exceeds transverse_bond_memory_cap_gb = "
                "{:.3f} GB (bubble-phase estimate {:.3f} GB [prep {:.3f} "
                "GB, pair-loop {:.3f} GB, peak = max of the two "
                "sub-phases], solve-phase estimate {:.3f} GB; overall "
                "peak = max of the bubble and solve phases, they do not "
                "overlap in time). The solve-phase estimate is "
                "(3 + K_solve) * Nq * ND**2 * 16 bytes with ND = "
                "B*norb**2 = {}, Nq = {}, K_solve = {}. The bubble "
                "prep-phase estimate is 16 * Nq * (ND**2 + 4*P*(Nmat+2)) "
                "bytes and the pair-loop estimate is 16 * Nq * (2*ND**2 "
                "+ 3*Nmat*P**2 + 2*Nmat*P + 8*P) bytes, both with P = "
                "norb**2 and Nmat = {} (all are numpy-backend "
                "conservative allowances, not guaranteed LAPACK/BLAS/FFT "
                "bounds). Restrict or remove some off-site "
                "CoulombInter/Ising/Exchange declarations to shrink the "
                "channel set B, use a coarser k-mesh, reduce Nmat (a "
                "lever for BOTH bubble sub-phase estimates, since the "
                "solve-phase estimate does not depend on it), or "
                "raise transverse_bond_memory_cap_gb. "
                "transverse_bond_max_shells only helps if the outer "
                "shells you would drop are declared-zero -- "
                "resolve_transverse_topology refuses any truncation that "
                "would discard a declared nonzero off-site "
                "coefficient.".format(
                    peak_bytes / 1.0e9, self.transverse_bond_memory_cap_gb,
                    bubble_bytes / 1.0e9, prep_bytes / 1.0e9,
                    pair_bytes / 1.0e9, solve_bytes / 1.0e9,
                    ND, Nq, K_solve, Nmat))

        if op_count > 1.0e12:
            logger.warning(
                "[mode.param] transverse_bond_channels: estimated solve "
                "op-count Nq*ND**3 = %.3e exceeds 1e12 (ND = %d, Nq = "
                "%d); this run fits the memory cap but may be "
                "compute-prohibitive.", op_count, ND, Nq)


def run(*, input_dict: Optional[dict] = None, input_file: Optional[str] = None):
    """Main entry point for running RPA calculations.

    Parameters
    ----------
    input_dict : dict, optional
        Dictionary containing input parameters. Must be provided if input_file is None.
    input_file : str, optional
        Path to input file. Not used if input_dict is provided.

    Raises
    ------
    RuntimeError
        If input_dict is None.

    Notes
    -----
    The input dictionary should contain the following sections:
    - log: Logging configuration
        - print_level: Verbosity level (default: 1)
        - print_step: Print frequency (default: 1)
    - mode: Calculation mode configuration
        - mode: Must be "RPA" for RPA calculations
    - file: File paths configuration
        - input: Input file paths
            - path_to_input: Input directory (default: "")
        - output: Output file paths
            - path_to_output: Output directory (default: "output")

    The function performs the following steps:
    1. Initializes logging and file paths
    2. Reads interaction definitions
    3. Creates RPA solver instance
    4. Performs RPA calculations
    5. Saves results
    """
    logger = logging.getLogger("run")

    if input_dict is None:
        raise RuntimeError("input_dict not passed")

    # Initialize information about log
    info_log = input_dict.get("log", {})
    info_log["print_level"] = info_log.get("print_level", 1)
    info_log["print_step"] = info_log.get("print_step", 1)

    # Initialize information about mode
    info_mode = input_dict.get("mode", {})
    info_file = input_dict.get("file", {"input": {}, "output": {}})
    # Initialize information about input files
    info_inputfile = info_file.get("input", {})
    info_inputfile["path_to_input"] = info_inputfile.get("path_to_input", "")

    # Initialize information about output files
    info_outputfile = info_file.get("output", {})
    info_outputfile["path_to_output"] = info_outputfile.get("path_to_output", "output")
    path_to_output = info_outputfile["path_to_output"]
    os.makedirs(path_to_output, exist_ok=True)

    if "mode" not in info_mode:
        logger.error("mode is not defined in [mode].")
        sys.exit(1)

    mode = info_mode["mode"]

    if mode == "RPA":
        logger.info("RPA mode")

        logger.info("Read interaction definitions from files")
        read_io = read_input_k.QLMSkInput(info_inputfile)
        ham_info = read_io.get_param("ham")

        solver = RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        green_info.update( solver.read_init(info_inputfile) )

        logger.info("Start RPA calculation")
        solver.solve(green_info, path_to_output)

        logger.info("Save calculation results")
        solver.save_results(info_outputfile, green_info)

        logger.info("all done.")

    else:
        logger.warning("mode is incorrect: mode={}.".format(mode))
        sys.exit(0)

    pass

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="increase verbosity", action="count", default=0)
    parser.add_argument("-q", "--quiet", help="decrease verbosity", action="count", default=0)
    parser.add_argument("input_file", nargs='?', default="input.toml", help="parameter file in TOML format")
    args = parser.parse_args()

    fmt = "%(asctime)s %(levelname)s %(name)s: %(message)s"
    log_level = 20 - (args.verbose - args.quiet) * 10
    logging.basicConfig(level=log_level, format=fmt)

    try:
        import toml
        params = toml.load(args.input_file)
    except Exception as e:
        print(e)
        sys.exit(1)

    run(input_dict=params)
