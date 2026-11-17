"""The schema-4 spin-orbital bond archive (design section 6.2).

The spin-orbital FLEX with bond channels (``[mode] enable_spin_orbital =
true``, ``[mode.param] longitudinal_bond_channels = true`` and
``longitudinal_bond_output_full = true``) writes the bond-resolved
susceptibility ``chi_w`` of its final Green function into one archive: an
index ``.npz`` holding every small member and ``Gamma_bond``, and the slot
blocks ``chi_w`` (and the optional ``chibar_w``) either inline (layout
``"npz"``) or as ``.npy`` sidecar files next to the index (layout
``"sidecar"``), named ``<block>.<archive_nonce>.npy``. The companion Green
function the susceptibility was built from is named by ``green_file`` (a
path relative to the index's directory).

This module holds the archive's fixed values and message texts, the layout
choice (:func:`archive_layout`), the byte counts the memory plan uses
(:func:`slot_member_bytes`, :func:`inline_max_bytes`), the writer
(:func:`write_so_bond_archive`) and the validating reader
(:func:`load_so_bond_archive`, returning a :class:`SOBondArchive` whose
:meth:`SOBondArchive.read_slot` reads the slot payloads). The rules the writer and the reader share
-- the ``green_file`` rule (:func:`green_file_ok`), the pole group
(:func:`pole_group_problem`) and the reading of an index's sidecar names
(:func:`index_sidecars`) -- are defined once, here.

The writer never copies a slot array: a sidecar is written straight from
the array, an inline member costs numpy's zip writer one chunk of at most
16 MiB. The index is written last, to a temporary name, and renamed over
the previous one; a failure before that rename leaves the directory as it
was. After the rename the sidecars of the replaced archive are removed,
except those another index in the directory may still name. Concurrent
writers to one directory are not supported.

This module does not import the FLEX solver.
"""
import math
import os
import re
import secrets
import types
from pathlib import PurePath
from typing import NamedTuple, Optional

import numpy as np

from . import fingerprint as _fp
from . import so_bond
from .eliashberg_bond_io import _tmp_name

# ---------------------------------------------------------------- fixed values
SCHEMA = 4
INDEX_ORDER = "bond_major_pair_code_slot"
FREQ_AXIS = "uniform_boson_centered"
GREEN_LAYOUT = _fp.GREEN_LAYOUT
INDEX_CONVENTION = "spin_block"
PRODUCER_STAGE = "flex_so_bond"
PRODUCER_LABEL = "hwave FLEX (spin-orbital, bond channels)"
#: the storage of the slot blocks: one npz, or an index with ``.npy`` sidecars
OUTPUT_LAYOUTS = ("npz", "sidecar")
#: the slot layout of the blocks (``so_bond.LAYOUTS``)
BLOCK_LAYOUTS = ("dense", "sz_packed")
#: the slot payloads; ``Gamma_bond`` is always stored inline
SLOT_PAYLOADS = ("Gamma_bond", "chi_w", "chibar_w")
SIDECAR_PAYLOADS = ("chi_w", "chibar_w")
#: the pole members (the final susceptibility past its pole), all or none
POLE_MEMBERS = ("rpa_past_pole", "rpa_pole_indicator_spinful", "rpa_pole_q_spinful")
#: the operational scalars the caller gives (the bond topology members
#: ``delta_r``, ``reverse``, ``slot_convention``, ``max_shells`` come from
#: the provenance only)
SCALAR_KEYS = ("norb_phys", "nd", "nmat", "B", "cell_shape", "beta", "types")
MOMENTUM_KEYS = ("momentum_convention", "wavevector_unit", "wavevector_index")
TYPES_DTYPE = "<U12"

_ENTRY = 16                       # bytes per complex128 entry
_ZIP_CHUNK = 16 * 2 ** 20         # numpy's zip writer copies one chunk of this size

# ---------------------------------------------------------------- message texts
SO_BOND_FLEX_HEAD = ("[mode] enable_spin_orbital = true with [mode.param] "
                     "longitudinal_bond_channels = true (FLEX)")
_H = SO_BOND_FLEX_HEAD

T_NOT_SO_ARCHIVE = (
    "'{path}' is not a spin-orbital bond archive (bond_archive_schema {schema}, "
    "spin_orbital {so}); write it with FLEX, [mode] enable_spin_orbital = true, "
    "[mode.param] longitudinal_bond_channels = true and longitudinal_bond_output_full = "
    "true.")
T_COMPANION_MISSING = (
    "companion green file '{green}' named by '{archive}' is missing (keep the archive, its "
    "sidecars and the green file together, or re-run FLEX)")
T_NEED_GREEN = (
    _H + ": longitudinal_bond_output_full = true needs [file.output] green (the bond archive "
    "names the Green function its susceptibility was built from as its companion "
    "green_file; [file.output] keys are case-sensitive); add green to [file.output].")
T_GREEN_DIR = (
    _H + ": the bond archive '{archive}' and its companion green file '{green}' must be in "
    "the same directory or the green file below it (green_file is stored as a relative "
    "path without '..').")
T_NPZ_TOO_LARGE = (
    _H + ": longitudinal_bond_output_layout = \"npz\" cannot hold the bond archive member "
    "{member} ({size}, above the small-array bound {bound} = min(1 GiB, cap / 8)); set "
    "longitudinal_bond_output_layout = \"sidecar\" or remove it (the sidecar layout is "
    "then chosen).")
T_FORCED = (
    "spin-orbital bond archive: sidecar layout (the member {member} of {size} is above the "
    "small-array bound {bound}); the susceptibility is written as .npy files next to the "
    "index.")
T_PASS = (
    "spin-orbital bond archive: chi_w is chi_B of the final Green function (one extra "
    "susceptibility pass after the self-consistency).")
T_PASS_REFUSAL_NOTE = (
    " This happened in the final susceptibility pass of the bond archive; set "
    "longitudinal_bond_output_full = false to skip that pass.")
T_OLD_SIDECARS = (
    "spin-orbital bond archive '{path}': removed {n} sidecar file(s) of the archive it "
    "replaced: {removed}; kept {k} (named by another index, a malformed candidate or this "
    "one): {kept}; kept by: {kept_by}{nonso}")
#: :data:`T_OLD_SIDECARS` when no index was replaced (only non-spin-orbital
#: sidecars at the stem to report)
T_OLD_SIDECARS_NO_INDEX = T_OLD_SIDECARS.replace(" of the archive it replaced", "")
#: the ``{nonso}`` field of :data:`T_OLD_SIDECARS` when such files exist
T_NONSO_TAIL = "; the non-spin-orbital sidecars {names} at this stem are left in place"
T_OLD_SIDECARS_FAILED = (
    "spin-orbital bond archive '{path}': could not remove the old sidecar file '{name}' "
    "({error}); remove it by hand if no other archive uses it.")
#: the removal of the replaced archive's sidecars failed as a whole (the new
#: archive is complete)
T_OLD_SIDECARS_SCAN_FAILED = (
    "spin-orbital bond archive '{path}': the removal of the replaced archive's sidecar "
    "files failed ({error}); the new archive is complete; remove stale sidecar files by "
    "hand if no other archive uses them.")
T_ARCHIVE_KEY_IGNORED = (
    "[file.output] longitudinal_bond is written only with [mode.param] "
    "longitudinal_bond_output_full = true; no bond archive is written by this run.")
T_NO_ARCHIVE_ITER0 = (
    "spin-orbital bond archive: IterationMax = 0 runs no map, so no bond archive is "
    "written.")
T_REMEDY_ARCHIVE = (
    "longitudinal_bond_output_full = true adds the archive phases (archive_bubble, "
    "archive_dressing, archive_write); setting it to false removes them.")
T_OLD_INDEX_REMOVED = (
    "spin-orbital bond archive '{path}': writing the new index failed after the companion "
    "green file '{green}' was replaced; the previous index was removed so that no archive "
    "names a different Green function ({error}).")
T_OLD_INDEX_KEPT = (
    "spin-orbital bond archive '{path}': writing the new index failed after the companion "
    "green file '{green}' was replaced, and the previous index could not be removed "
    "({error}); it names a different Green function, so load_so_bond_archive refuses it; "
    "remove it or re-run FLEX.")
T_CLEANUP_FAILED = (
    "spin-orbital bond archive: could not remove '{name}' while cleaning up after a failed "
    "write ({error}); remove it by hand.")

_W = "write_so_bond_archive: "


def format_bytes(n) -> str:
    """``"{:.2f} MiB"`` of ``n / 2**20`` below 1 GiB, else ``"{:.2f} GiB"``."""
    if n < 2 ** 30:
        return "{:.2f} MiB".format(n / 2 ** 20)
    return "{:.2f} GiB".format(n / 2 ** 30)


# ------------------------------------------------------------ sizes and layout
def _positive(name, v, who):
    if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or int(v) < 1:
        raise ValueError("{}: {} must be an integer >= 1, got {!r}".format(who, name, v))
    return int(v)


def _sector_sizes(B, nd, layout):
    """The trailing slot sizes: ``(ND,)`` dense, ``(n_0, n_+1, n_-1)`` packed."""
    ND = B * nd * nd
    if layout == "dense":
        return (ND,)
    # n_0 = 2 B norb_phys^2, n_+1 = n_-1 = B norb_phys^2 (nd = 2 norb_phys)
    return (ND // 2, ND // 4, ND // 4)


def _block_names(p, layout):
    if layout == "dense":
        return [p]
    return [so_bond._packed_key(p, s) for s in so_bond.SECTOR_ORDER]


def slot_member_bytes(*, nmat, nvol, B, nd, layout) -> dict:
    """The bytes of each ``chi_w`` block: ``{"chi_w": 16 nmat nvol ND^2}``
    (``layout = "dense"``, ``ND = B nd^2``) or ``{"chi_w@0": ...,
    "chi_w@+1": ..., "chi_w@-1": ...}`` with the sector sizes ``n_0 = B
    nd^2 / 2``, ``n_+1 = n_-1 = B nd^2 / 4`` (``"sz_packed"``)."""
    who = "slot_member_bytes"
    nmat, nvol, B = (_positive(k, v, who) for k, v in (("nmat", nmat), ("nvol", nvol),
                                                        ("B", B)))
    nd = _positive("nd", nd, who)
    if layout not in BLOCK_LAYOUTS:
        raise ValueError("slot_member_bytes: layout must be 'dense' or 'sz_packed', got "
                         "{!r}".format(layout))
    return {name: _ENTRY * nmat * nvol * n * n
            for name, n in zip(_block_names("chi_w", layout), _sector_sizes(B, nd, layout))}


def inline_max_bytes(*, nmat, nvol, B, nd, block_layout, output_layout, n_transfer=0,
                     n_extern=0, n_rec=0, n_topology=0) -> int:
    """The largest member the archive stores inline, from the shapes: the
    ``Gamma_bond`` blocks (``16 nvol n^2``), the ``chi_w`` blocks in the npz
    layout (:func:`slot_member_bytes`), ``wavevector_index`` (int64
    ``(nvol, 3)``) and the provenance tables (``prov_transfer_keys`` int64
    ``(n_transfer, 5)`` and its values, the Extern tables, the interaction
    records -- ``prov_rec_cells`` int64 ``(n_rec, 4, 3)`` the largest per row,
    ``prov_rec_types`` bounded by ``<U12`` -- and ``prov_topology`` int64
    ``(n_topology, 3)``), the sector index and the bond tables."""
    who = "inline_max_bytes"
    nmat, nvol, B, nd = (_positive(k, v, who) for k, v in
                         (("nmat", nmat), ("nvol", nvol), ("B", B), ("nd", nd)))
    if block_layout not in BLOCK_LAYOUTS:
        raise ValueError("inline_max_bytes: block_layout must be 'dense' or 'sz_packed', got "
                         "{!r}".format(block_layout))
    if output_layout not in OUTPUT_LAYOUTS:
        raise ValueError("inline_max_bytes: output_layout must be 'npz' or 'sidecar', got "
                         "{!r}".format(output_layout))
    counts = []
    for name, v in (("n_transfer", n_transfer), ("n_extern", n_extern), ("n_rec", n_rec),
                    ("n_topology", n_topology)):
        if isinstance(v, bool) or not isinstance(v, (int, np.integer)) or int(v) < 0:
            raise ValueError("inline_max_bytes: {} must be an integer >= 0, got {!r}"
                             .format(name, v))
        counts.append(int(v))
    n_transfer, n_extern, n_rec, n_topology = counts
    sizes = [_ENTRY * nvol * n * n for n in _sector_sizes(B, nd, block_layout)]
    if output_layout == "npz":
        sizes += list(slot_member_bytes(nmat=nmat, nvol=nvol, B=B, nd=nd,
                                        layout=block_layout).values())
    ND = B * nd * nd
    sizes += [8 * 3 * nvol,                                   # wavevector_index
              8 * 5 * n_transfer, _ENTRY * n_transfer,        # prov_transfer_*
              8 * 5 * n_extern, _ENTRY * n_extern,            # prov_extern_*
              8 * 4 * 3 * n_rec, 8 * 4 * n_rec, _ENTRY * n_rec,
              np.dtype(TYPES_DTYPE).itemsize * n_rec,         # prov_rec_*
              8 * 3 * n_topology,                             # prov_topology
              8 * ND if block_layout == "sz_packed" else 0,   # sz_sector_index
              8 * 3 * B]                                      # delta_r
    return int(max(sizes))


class ArchiveLayout(NamedTuple):
    """The output layout chosen by :func:`archive_layout`: ``layout``
    (``"npz"`` or ``"sidecar"``), whether it was ``forced`` to the sidecar
    layout (no request, a member above the bound), the largest ``member``
    and its ``member_bytes``, and the small-array ``bound``."""
    layout: str
    forced: bool
    member: str
    member_bytes: int
    bound: float


def archive_layout(*, member_bytes, small_array_bytes, request) -> ArchiveLayout:
    """The output layout of the archive: the single-file npz unless a slot
    block is larger than the small-array bound, and never against an
    explicit choice.

    ``member_bytes`` maps the slot blocks to their bytes
    (:func:`slot_member_bytes`; ``Gamma_bond`` is always inline and not
    listed), ``small_array_bytes`` is the plan's small-array bound
    (``so_bond.small_array_bound``) and ``request`` the configured
    ``longitudinal_bond_output_layout`` (``None`` when unset; compared
    after ``strip().lower()``). Unset: ``"npz"`` unless the largest member
    is above the bound, then ``"sidecar"`` with ``forced`` (the caller logs
    :data:`T_FORCED`). An explicit ``"npz"`` above the bound is refused
    (``ValueError`` :data:`T_NPZ_TOO_LARGE`); an explicit choice is never
    overridden."""
    if not member_bytes:
        raise ValueError("archive_layout: member_bytes is empty")
    if request is None:
        req = None
    else:
        req = request.strip().lower() if isinstance(request, str) else None
        if req not in OUTPUT_LAYOUTS:
            raise ValueError("archive_layout: request must be None, 'npz' or 'sidecar', got "
                             "{!r}".format(request))
    member, size = max(sorted(member_bytes.items()), key=lambda kv: kv[1])
    above = size > small_array_bytes
    if req is None:
        return ArchiveLayout("sidecar" if above else "npz", above, member, size,
                             small_array_bytes)
    if req == "npz" and above:
        raise ValueError(T_NPZ_TOO_LARGE.format(member=member, size=format_bytes(size),
                                                bound=format_bytes(small_array_bytes)))
    return ArchiveLayout(req, False, member, size, small_array_bytes)


# ------------------------------------------------------------ shared rules
def green_file_ok(green_file, archive_dir) -> bool:
    """Whether ``green_file`` names a companion file below ``archive_dir``:
    a non-empty string, not absolute, without a ``..`` part, not ``.``,
    without a trailing separator and not an existing directory there."""
    if not isinstance(green_file, str) or not green_file:
        return False
    seps = tuple(s for s in ("/", os.sep, os.altsep) if s)
    if os.path.isabs(green_file) or PurePath(green_file).is_absolute() \
            or green_file.startswith(seps):
        return False
    if ".." in PurePath(green_file).parts or green_file.endswith(seps):
        return False
    if os.path.normpath(green_file) == ".":
        return False
    return not os.path.isdir(os.path.join(archive_dir, green_file))


def pole_group_problem(members) -> Optional[tuple]:
    """``None`` when ``members`` (the archive's members named ``rpa_*``) is
    empty or exactly the pole group :data:`POLE_MEMBERS` --
    ``rpa_past_pole`` a 0-d bool True, ``rpa_pole_indicator_spinful`` a 0-d
    finite float64, ``rpa_pole_q_spinful`` a finite float64 ``(3,)`` --
    else the first problem as ``(check, stored, expected)``."""
    names = sorted(members)
    stray = [n for n in names if n not in POLE_MEMBERS]
    if stray:
        return ("member {} is not a pole member".format(stray[0]), "present", "absent")
    if not names:
        return None
    for n in POLE_MEMBERS:
        if n not in members:
            return ("pole member {} is missing".format(n), "absent", "present")
    past = np.asarray(members["rpa_past_pole"])
    if past.ndim != 0 or past.dtype.kind != "b" or not bool(past[()]):
        return ("rpa_past_pole is not true", repr(past.tolist()), "True")
    ind = np.asarray(members["rpa_pole_indicator_spinful"])
    what = "rpa_pole_indicator_spinful is not a finite float64 scalar"
    if ind.ndim != 0 or not _fp._schema_dtype_ok(ind.dtype, "float64"):
        return (what, "{}/{}d".format(ind.dtype, ind.ndim), "float64/0d")
    if not np.isfinite(ind[()]):
        return (what, repr(float(ind[()])), "finite")
    q = np.asarray(members["rpa_pole_q_spinful"])
    what = "rpa_pole_q_spinful is not a finite float64 (3,) array"
    if q.shape != (3,) or not _fp._schema_dtype_ok(q.dtype, "float64"):
        return (what, "{}/{}".format(q.dtype, q.shape), "float64/(3,)")
    if not np.all(np.isfinite(q)):
        return (what, repr(float(q[~np.isfinite(q)][0])), "finite")
    return None


_BLOCK = r"(?:Gamma_bond|chi_w|chibar_w)(?:@(?:0|\+1|-1))?"
_BLOCK_RE = re.compile(r"^{}$".format(_BLOCK))
_SIDECAR_RE = re.compile(r"^({})\.([0-9a-f]{{32}})\.npy$".format(_BLOCK))


def sidecar_name(block, nonce) -> str:
    """The base name ``<block>.<nonce>.npy`` of a sidecar file."""
    return "{}.{}.npy".format(block, nonce)


def parse_sidecar_name(name):
    """``(block, nonce)`` of a sidecar base name ``<block>.<nonce>.npy``
    (``block`` a slot block, ``nonce`` 32 lowercase hex), else ``None``."""
    m = _SIDECAR_RE.match(name) if isinstance(name, str) else None
    return (m.group(1), m.group(2)) if m else None


def _schema_four(data, files):
    if "bond_archive_schema" not in files:
        return False
    a = np.asarray(data["bond_archive_schema"])
    return a.ndim == 0 and a.dtype.kind in "iu" and int(a[()]) == SCHEMA


def _zero_d_string(data, name, files):
    if name not in files:
        return None
    a = np.asarray(data[name])
    return str(a[()]) if a.ndim == 0 and a.dtype.kind == "U" else None


class IndexSidecars(NamedTuple):
    """What :func:`index_sidecars` reads from an npz: whether it is a
    schema-4 index (``bond_archive_schema`` a 0-d integer equal to 4), its
    ``bond_archive_layout`` and ``archive_nonce`` (0-d strings, else
    ``None``), the file names its ``<block>_file`` members carry, and
    whether it is ``malformed`` (any ``*_file`` member is not a 0-d
    string)."""
    index: bool
    layout: Optional[str]
    nonce: Optional[str]
    names: frozenset
    malformed: bool

    def own_sidecars(self):
        """The names this index's sidecar layout writes: plain base names
        ``<block>.<archive_nonce>.npy``; empty unless it is a sidecar
        index."""
        if not self.index or self.layout != "sidecar" or self.nonce is None:
            return frozenset()
        return frozenset(n for n in self.names
                         if (parse_sidecar_name(n) or (None, None))[1] == self.nonce)


def index_sidecars(data) -> IndexSidecars:
    """The sidecar names an archive ``data`` (an open npz) carries; a file
    that is not a schema-4 index carries none."""
    files = set(getattr(data, "files", data))
    if not _schema_four(data, files):
        return IndexSidecars(False, None, None, frozenset(), False)
    names, malformed = set(), False
    for f in sorted(files):
        if not (isinstance(f, str) and f.endswith("_file")):
            continue
        v = _zero_d_string(data, f, files)
        if v is None:
            # any *_file member that is not a 0-d string makes the index
            # malformed (it keeps every old name)
            malformed = True
        elif _BLOCK_RE.match(f[:-5]):
            names.add(v)
    return IndexSidecars(True, _zero_d_string(data, "bond_archive_layout", files),
                         _zero_d_string(data, "archive_nonce", files), frozenset(names),
                         malformed)


def _open_index(path) -> IndexSidecars:
    """:func:`index_sidecars` of the npz at ``path`` (raises when it cannot
    be read as an npz)."""
    with np.load(path, allow_pickle=False) as d:
        if not hasattr(d, "files"):  # pragma: no cover - np.load of npy in a context
            raise ValueError("not an npz archive")
        return index_sidecars(d)


def _read_old(path):
    """``(replaced, old sidecar names)`` of the file at ``path`` before the
    write: ``replaced`` when it is a readable npz with a
    ``bond_archive_schema`` member; the names are those of a schema-4
    sidecar index (an unreadable file yields none)."""
    if not os.path.isfile(path):
        return False, frozenset()
    try:
        d = np.load(path, allow_pickle=False)
    except Exception:
        return False, frozenset()
    if not hasattr(d, "files"):
        return False, frozenset()
    try:
        with d:
            replaced = "bond_archive_schema" in d.files
            return replaced, index_sidecars(d).own_sidecars()
    except Exception:
        return False, frozenset()


# ------------------------------------------------------------------ the writer
class SOBondArchivePayload(NamedTuple):
    """What the final susceptibility pass hands to the writer: ``chi_w``
    and ``gamma_b`` (dense arrays or ``so_bond.PackedArrays``) and the pole
    members (``bond_channels.RpaPoleReport.npz_members()``, empty below the
    pole)."""
    chi_w: object
    gamma_b: object
    pole_members: dict


def _new_nonce():
    """A fresh 128-bit ``archive_nonce`` (32 lowercase hex)."""
    return secrets.token_hex(16)


def _is_int(v):
    return not isinstance(v, (bool, np.bool_)) and isinstance(v, (int, np.integer))


def _check_scalars(scalars, provenance):
    if not isinstance(scalars, dict):
        raise ValueError(_W + "scalars must be a dict, got {}".format(type(scalars).__name__))
    for key in sorted(scalars):
        if key in provenance or key in _fp.SO_BOND_MEMBERS:
            raise ValueError(_W + "scalars must not carry {} (it comes from provenance "
                             "only)".format(key))
    if set(scalars) != set(SCALAR_KEYS):
        raise ValueError(_W + "scalars must have exactly the keys {}, got {}".format(
            ", ".join(SCALAR_KEYS), ", ".join(sorted(scalars))))
    rules = (("norb_phys", ">= 1", lambda v: v >= 1),
             ("nd", "equal to 2 * norb_phys", None),
             ("nmat", "> 0 and even", lambda v: v > 0 and v % 2 == 0),
             ("B", ">= 1", lambda v: v >= 1))
    out = {}
    for key, rule, ok in rules:
        v = scalars[key]
        good = _is_int(v) and (ok(int(v)) if ok else int(v) == 2 * out["norb_phys"])
        if not good:
            raise ValueError(_W + "scalars[{!r}] must be an integer {}, got {!r}".format(
                key, rule, v))
        out[key] = int(v)
    beta = scalars["beta"]
    if not isinstance(beta, (float, np.floating)) or not math.isfinite(beta) or not beta > 0:
        raise ValueError(_W + "scalars['beta'] must be a finite float > 0, got {!r}"
                         .format(beta))
    out["beta"] = float(beta)
    cs = np.asarray(scalars["cell_shape"])
    if cs.shape != (3,) or cs.dtype.kind not in "iu" or not np.all(cs >= 1):
        raise ValueError(_W + "scalars['cell_shape'] must be an integer array of shape (3,) "
                         "with entries >= 1, got {!r}".format(cs.tolist()))
    out["cell_shape"] = np.array(cs, dtype=np.int64)
    ty = scalars["types"]
    good = (isinstance(ty, np.ndarray) and ty.ndim == 1 and ty.dtype == np.dtype(TYPES_DTYPE)
            and list(ty) == sorted(set(ty.tolist())))
    if not good:
        raise ValueError(_W + "scalars['types'] must be a sorted 1-d <U12 array of distinct "
                         "names, got {!r}".format(np.asarray(ty).tolist()))
    out["types"] = np.array(ty, dtype=TYPES_DTYPE)
    return out


def _check_provenance(provenance, sc):
    if not isinstance(provenance, dict):
        raise ValueError(_W + "provenance must be a dict, got {}".format(
            type(provenance).__name__))
    need = set(_fp.SO_PROVENANCE_MEMBERS + _fp.SO_BOND_MEMBERS)
    missing = sorted(need - set(provenance))
    extra = sorted(set(provenance) - need - set(_fp.OPTIONAL_MEMBERS))
    if missing or extra:
        raise ValueError(_W + "provenance must carry the spin-orbital and bond provenance "
                         "members; missing {}; unexpected {}".format(
                             ", ".join(missing) or "none", ", ".join(extra) or "none"))

    def text(k):
        a = np.asarray(provenance[k])
        return str(a[()]) if a.ndim == 0 else repr(a.tolist())

    if text("producer_stage") != PRODUCER_STAGE:
        raise ValueError(_W + "provenance producer_stage must be {!r}, got {!r}".format(
            PRODUCER_STAGE, text("producer_stage")))
    if text("prov_index_convention") != INDEX_CONVENTION:
        raise ValueError(_W + "provenance prov_index_convention must be {!r}, got {!r}"
                         .format(INDEX_CONVENTION, text("prov_index_convention")))
    dr = np.asarray(provenance["delta_r"])
    pairs = (("B", sc["B"], "len(delta_r)", int(dr.shape[0]) if dr.ndim else None),
             ("norb_phys", sc["norb_phys"], "prov_norb_phys",
              np.asarray(provenance["prov_norb_phys"])[()]),
             ("nmat", sc["nmat"], "state_nmat", np.asarray(provenance["state_nmat"])[()]),
             ("beta", sc["beta"], "state_beta", np.asarray(provenance["state_beta"])[()]))
    for name, v, pname, pv in pairs:
        if pv is not None and np.ndim(pv) != 0:
            raise ValueError(_W + "provenance {} must be a 0-d scalar, got shape {}".format(
                pname, np.shape(pv)))
        if pv is None or v != pv:
            raise ValueError(_W + "scalars disagree with provenance: {} = {!r} but {} = {!r}"
                             .format(name, v, pname, pv.item() if hasattr(pv, "item")
                                     else pv))
    pcs = np.asarray(provenance["prov_cell_shape"])
    if pcs.shape != (3,) or not np.array_equal(pcs, sc["cell_shape"]):
        raise ValueError(_W + "scalars disagree with provenance: cell_shape = {} but "
                         "prov_cell_shape = {}".format(sc["cell_shape"].tolist(),
                                                       pcs.tolist()))
    return np.asarray(provenance["reverse"])


def _check_block(name, a, shape):
    if a.dtype != np.complex128:
        raise ValueError(_W + "{} must be complex128, got {}".format(name, a.dtype))
    if not a.flags.c_contiguous:
        raise ValueError(_W + "{} must be C-contiguous".format(name))
    if tuple(a.shape) != tuple(shape):
        raise ValueError(_W + "{} must have shape {}, got {}".format(name, tuple(shape),
                                                                     tuple(a.shape)))


def _slot_blocks(name, arr, block_layout, sector_map, lead, ND, payload):
    """``{block name: array}`` of a slot array after its checks."""
    if block_layout == "dense":
        if not isinstance(arr, np.ndarray):
            raise ValueError(_W + "{} of a dense archive must be a numpy array, got {}"
                             .format(name, type(arr).__name__))
        _check_block(name, arr, lead + (ND, ND))
        return {payload: arr}
    if not isinstance(arr, so_bond.PackedArrays):
        raise ValueError(_W + "{} of a sz_packed archive must be a PackedArrays, got {}"
                         .format(name, type(arr).__name__))
    if not so_bond._same_sector_map(arr.sector_map, sector_map):
        raise ValueError(_W + "{} is packed on another sector map".format(name))
    out = {}
    for s, b, n in zip(so_bond.SECTOR_ORDER, arr.blocks, sector_map.sizes):
        _check_block("{}@{}".format(name, s), b, lead + (n, n))
        out[so_bond._packed_key(payload, s)] = b
    return out


def _check_momentum(momentum, nvol):
    from .rpa import MOMENTUM_CONVENTION
    if not isinstance(momentum, dict) or set(momentum) != set(MOMENTUM_KEYS):
        got = ", ".join(sorted(momentum)) if isinstance(momentum, dict) else \
            type(momentum).__name__
        raise ValueError(_W + "momentum must have exactly the keys {}, got {}".format(
            ", ".join(MOMENTUM_KEYS), got))
    mc = momentum["momentum_convention"]
    if not isinstance(mc, str) or mc != MOMENTUM_CONVENTION:
        raise ValueError(_W + "momentum_convention must be {!r}, got {!r}".format(
            MOMENTUM_CONVENTION, mc))
    u = np.asarray(momentum["wavevector_unit"])
    if u.shape != (3, 3) or u.dtype.kind not in "fiu" or not np.all(np.isfinite(u)):
        raise ValueError(_W + "wavevector_unit must be a finite real array of shape (3, 3), "
                         "got {} {}".format(u.dtype, u.shape))
    w = np.asarray(momentum["wavevector_index"])
    if w.shape != (nvol, 3) or w.dtype.kind not in "iu":
        raise ValueError(_W + "wavevector_index must be an integer array of shape ({}, 3), "
                         "got {} {}".format(nvol, w.dtype, w.shape))
    return {"momentum_convention": np.str_(mc),
            "wavevector_unit": np.array(u, dtype=np.float64),
            "wavevector_index": np.array(w, dtype=np.int64)}


def _remove_quietly(full, log):
    """Remove ``full``; a missing file is fine; another error is one
    :data:`T_CLEANUP_FAILED` WARNING (never raised)."""
    try:
        os.remove(full)
    except FileNotFoundError:
        pass
    except BaseException as exc:  # noqa: B036 - never mask the original exception
        log.warning(T_CLEANUP_FAILED.format(name=full, error=exc))


def write_so_bond_archive(path, *, chi_w, gamma_b, block_layout, output_layout, sector_map,
                          scalars, provenance, green_file, momentum, scf_converged,
                          chibar_w=None, extra_members=None, commit_companion=None,
                          log) -> dict:
    """Write the schema-4 spin-orbital bond archive to ``path`` (ending in
    ``.npz``).

    ``chi_w`` (``(nmat, nvol, ., .)``), ``gamma_b`` (``(nvol, ., .)``) and
    the optional ``chibar_w`` are dense complex128 C-order arrays
    (``block_layout = "dense"``, trailing ``(ND, ND)``, ``ND = B nd^2``,
    ``sector_map = None``) or ``so_bond.PackedArrays`` on ``sector_map``,
    the sector map of the provenance's ``(reverse, nd)``
    (``"sz_packed"``). ``output_layout`` ``"npz"`` stores every block
    inline; ``"sidecar"`` writes the ``chi_w``/``chibar_w`` blocks as
    ``<block>.<archive_nonce>.npy`` next to the index (``<block>_file``,
    ``<block>_nbytes`` in the index); ``Gamma_bond`` is always inline.

    ``scalars`` holds exactly ``norb_phys, nd, nmat, B`` (integers),
    ``cell_shape`` (integer ``(3,)``), ``beta`` (a finite float > 0) and
    ``types`` (the sorted distinct interaction types, ``<U12``);
    ``provenance`` the spin-orbital and bond provenance members of stage
    ``flex_so_bond`` (``fingerprint.so_provenance_members``; ``delta_r``,
    ``reverse``, ``slot_convention`` and ``max_shells`` come from it only),
    in agreement with the scalars. ``green_file`` names the companion Green
    archive relative to the index's directory (:func:`green_file_ok`).
    ``momentum`` holds ``momentum_convention``, ``wavevector_unit`` and
    ``wavevector_index``; ``scf_converged`` a bool. ``extra_members`` is
    empty or exactly the pole group (:func:`pole_group_problem`).

    Every precondition is a ``ValueError`` before anything is written.
    Order: (i) the sidecar names of the index at ``path`` (if any) are read;
    (ii) the new sidecars, each registered before it is written; (iii) the
    index to a temporary name next to ``path``, then ``commit_companion()``
    when given; (iv) the rename over ``path``. Any exception in (ii)-(iv)
    removes the new sidecars and the temporary index (a removal error is
    one :data:`T_CLEANUP_FAILED` WARNING) and propagates; when the companion
    was already committed the previous index at ``path`` is removed
    (:data:`T_OLD_INDEX_REMOVED`, or :data:`T_OLD_INDEX_KEPT` when that
    fails). (v) After the rename, best effort and never raising: the
    replaced index's sidecars are removed except the names this index or
    another schema-4 index in the directory carries (an unreadable or
    malformed ``*.npz`` there keeps them all), one INFO
    :data:`T_OLD_SIDECARS`; the non-spin-orbital sidecars at the stem are
    listed, never removed.

    Returns ``{"path", "archive_nonce", "bond_archive_layout",
    "block_layout", "sidecars", "removed", "kept", "kept_by", "nonso"}``
    (the lists are sorted base names)."""
    # ---- preconditions: nothing is written before all of them hold
    if not isinstance(path, (str, os.PathLike)) or not os.fspath(path).endswith(".npz"):
        raise ValueError(_W + "path must end in .npz, got {!r}".format(path))
    path = os.fspath(path)
    if block_layout not in BLOCK_LAYOUTS:
        raise ValueError(_W + "block_layout must be 'dense' or 'sz_packed', got {!r}"
                         .format(block_layout))
    if output_layout not in OUTPUT_LAYOUTS:
        raise ValueError(_W + "output_layout must be 'npz' or 'sidecar', got {!r}"
                         .format(output_layout))
    out_dir = os.path.dirname(path)
    if not green_file_ok(green_file, out_dir or "."):
        raise ValueError(_W + "green_file must be a relative file path without '..', got "
                         "{!r}".format(green_file))
    sc = _check_scalars(scalars, provenance if isinstance(provenance, dict) else {})
    reverse = _check_provenance(provenance, sc)
    nd, B, nmat = sc["nd"], sc["B"], sc["nmat"]
    nvol = int(np.prod(sc["cell_shape"]))
    ND = B * nd * nd
    packed = block_layout == "sz_packed"
    if packed:
        want = so_bond.sz_sector_map(types.SimpleNamespace(reverse=reverse), nd)
        if not isinstance(sector_map, so_bond.SectorMap):
            raise ValueError(_W + "a sz_packed archive needs the sector map of the "
                             "provenance's (reverse, nd), got {}".format(
                                 type(sector_map).__name__))
        if not so_bond._same_sector_map(sector_map, want):
            raise ValueError(_W + "a sz_packed archive needs the sector map of the "
                             "provenance's (reverse, nd), got another SectorMap")
    elif sector_map is not None:
        raise ValueError(_W + "a dense archive takes sector_map = None, got {}".format(
            type(sector_map).__name__))
    blocks = {}
    blocks.update(_slot_blocks("gamma_b", gamma_b, block_layout, sector_map, (nvol,), ND,
                               "Gamma_bond"))
    slot = {}
    slot.update(_slot_blocks("chi_w", chi_w, block_layout, sector_map, (nmat, nvol), ND,
                             "chi_w"))
    if chibar_w is not None:
        slot.update(_slot_blocks("chibar_w", chibar_w, block_layout, sector_map,
                                 (nmat, nvol), ND, "chibar_w"))
    mom = _check_momentum(momentum, nvol)
    if not isinstance(scf_converged, (bool, np.bool_)):
        raise ValueError(_W + "scf_converged must be a bool, got {!r}".format(scf_converged))
    if commit_companion is not None and not callable(commit_companion):
        raise ValueError(_W + "commit_companion must be callable or None, got {}".format(
            type(commit_companion).__name__))
    extra = dict(extra_members or {})

    nonce = _new_nonce()
    head = {
        "bond_archive_schema": np.int64(SCHEMA),
        "spin_orbital": np.bool_(True),
        "index_convention": np.array(INDEX_CONVENTION, dtype="<U16"),
        "index_order": np.array(INDEX_ORDER, dtype="<U32"),
    }
    head.update(so_bond.sector_map_members(sector_map if packed else None))
    head.update({
        "freq_axis": np.array(FREQ_AXIS, dtype="<U32"),
        "bond_archive_layout": np.array(output_layout, dtype="<U8"),
        "norb_phys": np.int64(sc["norb_phys"]), "nd": np.int64(nd), "nmat": np.int64(nmat),
        "B": np.int64(B), "cell_shape": sc["cell_shape"], "beta": np.float64(sc["beta"]),
        "types": sc["types"],
    })
    tail = {"green_file": np.str_(green_file),
            "archive_nonce": np.array(nonce, dtype="<U32"),
            "green_layout": np.array(GREEN_LAYOUT, dtype="<U32")}
    tail.update(provenance)
    tail.update(mom)
    tail["scf_converged"] = np.bool_(scf_converged)
    taken = set(head) | set(blocks) | set(tail)
    for k in slot:
        taken |= {k, k + "_file", k + "_nbytes"}
    for k in sorted(extra):
        if k in taken:
            raise ValueError(_W + "extra_members key {!r} collides with an archive member"
                             .format(k))
    problem = pole_group_problem(extra)
    if problem is not None:
        raise ValueError(_W + "extra_members is not the pole group: {} ({} vs {})"
                         .format(*problem))
    # the archive as the reader sees it (the slot blocks inline): never
    # write one its structural check refuses
    try:
        _fp.validate_structure({**head, **blocks, **slot, **tail, **extra}, path, kind="bond")
    except ValueError as exc:
        raise ValueError(_W + "the archive would not pass the structural check: {}"
                         .format(exc)) from None

    # ---- (i) the replaced index's sidecar names
    replaced, old = _read_old(path)

    # ---- (ii)-(iv): the failure window
    side_dir = out_dir or "."
    registered = []
    tmp = None
    committed = False
    sidecar_members = {}
    try:
        if output_layout == "sidecar":
            for name, arr in slot.items():
                base = sidecar_name(name, nonce)
                full = os.path.join(side_dir, base)
                registered.append(full)
                np.save(full, arr, allow_pickle=False)
                sidecar_members[name + "_file"] = np.str_(base)
                sidecar_members[name + "_nbytes"] = np.int64(os.path.getsize(full))
            inline = {**head, **blocks, **sidecar_members, **tail, **extra}
        else:
            inline = {**head, **blocks, **slot, **tail, **extra}
        tmp = _tmp_name(path)
        with open(tmp, "wb") as fh:
            np.savez(fh, **inline)
        if commit_companion is not None:
            commit_companion()
            committed = True
        os.replace(tmp, path)
    except BaseException as exc:
        for full in registered:
            _remove_quietly(full, log)
        if tmp is not None:
            _remove_quietly(tmp, log)
        if committed:
            green = os.path.normpath(os.path.join(out_dir, green_file))
            try:
                os.remove(path)
            except FileNotFoundError:
                pass
            except BaseException as exc_rm:  # noqa: B036 - logged, never masks exc
                log.warning(T_OLD_INDEX_KEPT.format(path=path, green=green, error=exc_rm))
            else:
                log.warning(T_OLD_INDEX_REMOVED.format(path=path, green=green, error=exc))
        raise

    new_names = sorted(os.path.basename(f) for f in registered)
    result = {"path": path, "archive_nonce": nonce, "bond_archive_layout": output_layout,
              "block_layout": block_layout, "sidecars": new_names, "removed": [], "kept": [],
              "kept_by": [], "nonso": []}
    # ---- (v) best effort, never raises
    try:
        _clean_old_sidecars(path, side_dir, replaced, old, frozenset(new_names), result, log)
    except Exception as exc:
        log.warning(T_OLD_SIDECARS_SCAN_FAILED.format(path=path, error=exc))
    return result


def _scan_candidates(side_dir, path):
    """``(other names, keep all, kept_by candidates)`` of the other ``*.npz``
    files in ``side_dir`` (names starting with ``.`` and ``path`` itself
    skipped): the sidecar names each schema-4 index carries, and the
    unreadable or malformed candidates (each keeps every old name)."""
    me = os.path.basename(path)
    carried = {}
    bad = []
    for n in sorted(os.listdir(side_dir)):
        if not n.endswith(".npz") or n.startswith(".") or n == me:
            continue
        full = os.path.join(side_dir, n)
        if not os.path.isfile(full):
            continue
        try:
            info = _open_index(full)
        except Exception:
            bad.append(n)
            continue
        if not info.index:
            continue
        if info.malformed:
            bad.append(n)
        else:
            carried[n] = info.names
    return carried, bad


def _clean_old_sidecars(path, side_dir, replaced, old, new_names, result, log):
    removed, kept, kept_by = [], [], set()
    if old:
        carried, bad = _scan_candidates(side_dir, path)
        if bad:
            keep = set(old)
            kept_by.update(bad)
        else:
            keep = set(old) & set(new_names)
            for n, names in carried.items():
                if set(old) & names:
                    kept_by.add(n)
                    keep |= set(old) & names
        kept = sorted(keep)
        for base in sorted(set(old) - keep):
            full = os.path.join(side_dir, base)
            try:
                os.remove(full)
            except FileNotFoundError:
                continue
            except OSError as exc:
                log.warning(T_OLD_SIDECARS_FAILED.format(path=path, name=full, error=exc))
                continue
            removed.append(base)
    base = os.path.basename(path)
    stem = base[:-4]
    nonso = sorted(n for n in ("{}_chi_s_w.npy".format(stem), "{}_chi_c_w.npy".format(stem))
                   if os.path.isfile(os.path.join(side_dir, n)))
    result.update(removed=removed, kept=kept, kept_by=sorted(kept_by), nonso=nonso)
    if removed or kept or nonso:
        join = lambda v: ", ".join(sorted(v)) if v else "none"  # noqa: E731
        text = T_OLD_SIDECARS if replaced else T_OLD_SIDECARS_NO_INDEX
        log.info(text.format(
            path=path, n=len(removed), removed=join(removed), k=len(kept), kept=join(kept),
            kept_by=join(kept_by),
            nonso=T_NONSO_TAIL.format(names=", ".join(nonso)) if nonso else ""))


# ------------------------------------------------------------------ the reader
#: the longest stored value a refusal quotes (longer ones end in ``...``)
_REPR_BOUND = 80
_NONCE_RE = re.compile(r"^[0-9a-f]{32}$")
_GREEN_FILE_RULE = "a file name or a relative path below the archive directory"


def _bounded(text) -> str:
    """``text`` cut to at most :data:`_REPR_BOUND` characters."""
    text = str(text)
    if len(text) <= _REPR_BOUND:
        return text
    return text[:_REPR_BOUND - 3] + "..."


def _value_text(a) -> str:
    return _bounded(repr(a.tolist()))


def _is_so_archive(data, names, path):
    """(1b) The markers: ``bond_archive_schema`` a 0-d integer equal to 4
    and ``spin_orbital`` a 0-d bool True, else :data:`T_NOT_SO_ARCHIVE`."""
    shown = {}
    good = True
    for key, kinds, want in (("bond_archive_schema", "iu", SCHEMA), ("spin_orbital", "b", True)):
        if key not in names:
            shown[key] = "absent"
            good = False
            continue
        a = np.asarray(data[key])
        shown[key] = _value_text(a)
        if not (a.ndim == 0 and a.dtype.kind in kinds and a[()] == want):
            good = False
    if not good:
        raise _fp.SOArchiveError(T_NOT_SO_ARCHIVE.format(
            path=path, schema=shown["bond_archive_schema"], so=shown["spin_orbital"]))


#: (1b) the 0-d scalars in their order: (name, kind, range text, range rule)
_SCALAR_RULES = (
    ("nmat", "int", "> 0 and even", lambda v: v > 0 and v % 2 == 0),
    ("B", "int", ">= 1", lambda v: v >= 1),
    ("nd", "int", ">= 2 and even", lambda v: v >= 2 and v % 2 == 0),
    ("norb_phys", "int", ">= 1", lambda v: v >= 1),
    ("state_nmat", "int", None, None),
    ("prov_norb_phys", "int", None, None),
    ("beta", "float", "finite and > 0", lambda v: math.isfinite(v) and v > 0),
    ("state_beta", "float", "finite and > 0", lambda v: math.isfinite(v) and v > 0),
)


def _peek_scalars(data, names, malformed) -> dict:
    """(1b) The operational scalars, the cell and bond tables and their
    relations, from the small members only (no payload header is read)."""
    def member(k):
        if k not in names:
            malformed("member {} is missing".format(k), "absent", "present")
        return np.asarray(data[k])

    v = {}
    for k, kind, rule_text, rule in _SCALAR_RULES:
        a = member(k)
        kinds = "iu" if kind == "int" else "f"
        if a.ndim != 0 or a.dtype.kind not in kinds:
            malformed("member {} is not a 0-d {}".format(k, "integer" if kind == "int"
                                                         else "float"),
                      "{}/{}d".format(a.dtype, a.ndim), "{}/0d".format(kind))
        x = int(a[()]) if kind == "int" else float(a[()])
        if rule is not None and not rule(x):
            malformed("member {} is out of range".format(k), repr(x), rule_text)
        v[k] = x
    for k in ("cell_shape", "prov_cell_shape"):
        a = member(k)
        if a.shape != (3,) or a.dtype.kind not in "iu":
            malformed("member {} is not an integer (3,) array".format(k),
                      "{}/{}".format(a.dtype, a.shape), "int/(3,)")
        if not np.all(a >= 1):
            malformed("member {} is out of range".format(k), _value_text(a), "entries >= 1")
        v[k] = a.astype(np.int64)
    dr = member("delta_r")
    if dr.ndim != 2 or dr.shape[1] != 3 or dr.dtype.kind not in "iu":
        malformed("member delta_r is not an integer (B', 3) array",
                  "{}/{}".format(dr.dtype, dr.shape), "int/(B', 3)")

    # the relations, in this order; the copies of a recorded state field
    # use the structure check's own rule and wording
    if v["nd"] != 2 * v["norb_phys"]:
        malformed("nd differs from 2 * norb_phys", v["nd"], 2 * v["norb_phys"])
    _fp.check_state_relations(data, malformed.path, malformed.label, members=("norb_phys",))
    if v["B"] != dr.shape[0]:
        malformed("B differs from len(delta_r)", v["B"], dr.shape[0])
    _fp.check_state_relations(data, malformed.path, malformed.label, members=("nmat", "beta"))
    if not np.array_equal(v["cell_shape"], v["prov_cell_shape"]):
        malformed("cell_shape differs from prov_cell_shape", _value_text(v["cell_shape"]),
                  _value_text(v["prov_cell_shape"]))
    return v


class _Malformed:
    """The malformed refusal of :mod:`fingerprint` for ``path`` naming
    ``label``; keeps both for the shared relation check."""

    def __init__(self, path, label):
        self.path, self.label = path, label
        self._raise = _fp._malformed_raiser(path, label)

    def __call__(self, check, stored, expected):
        self._raise(check, stored, expected)


def _str_member(data, names, k, malformed) -> str:
    if k not in names:
        malformed("member {} is missing".format(k), "absent", "present")
    a = np.asarray(data[k])
    if a.ndim != 0 or a.dtype.kind != "U":
        malformed("member {} is not a 0-d string".format(k), "{}/{}d".format(a.dtype, a.ndim),
                  "str/0d")
    return str(a[()])


def _check_extras(data, names, path, sv, sources, malformed):
    """(3) The schema-4 members the structure check does not know; returns
    ``(bond_archive_layout, green_file, nvol)``."""
    for k, want in (("index_order", INDEX_ORDER), ("freq_axis", FREQ_AXIS)):
        got = _str_member(data, names, k, malformed)
        if got != want:
            malformed("{} is not this version's".format(k), _bounded(repr(got)), repr(want))
    out_layout = _str_member(data, names, "bond_archive_layout", malformed)
    if out_layout not in OUTPUT_LAYOUTS:
        malformed("bond_archive_layout is not a known layout", _bounded(repr(out_layout)),
                  "'npz' or 'sidecar'")
    blocks = [b for p in SLOT_PAYLOADS for b in _block_names(p, sv["layout"]) if b in sources]
    for b in blocks:
        side = sources[b][0] != "inline"
        want_side = out_layout == "sidecar" and not b.startswith("Gamma_bond")
        if side != want_side:
            malformed("payload member {} does not follow bond_archive_layout {}".format(
                b, out_layout), "sidecar" if side else "inline",
                "sidecar" if want_side else "inline")
    gl = _str_member(data, names, "green_layout", malformed)
    if gl != GREEN_LAYOUT:
        malformed("green_layout is not the canonical layout", _bounded(repr(gl)),
                  repr(GREEN_LAYOUT))
    if "scf_converged" not in names:
        malformed("member scf_converged is missing", "absent", "present")
    a = np.asarray(data["scf_converged"])
    if a.ndim != 0 or a.dtype.kind != "b":
        malformed("member scf_converged is not a 0-d bool", "{}/{}d".format(a.dtype, a.ndim),
                  "bool/0d")
    if "types" not in names:
        malformed("member types is missing", "absent", "present")
    a = np.asarray(data["types"])
    if a.ndim != 1 or a.dtype.kind != "U":
        malformed("member types is not a 1-d string array", "{}/{}d".format(a.dtype, a.ndim),
                  "str/1d")
    nvol = int(np.prod(sv["cell_shape"]))
    for b in blocks:
        lead = tuple(sources[b][1][:-2])
        want = (nvol,) if b.startswith("Gamma_bond") else (sv["nmat"], nvol)
        if lead != want:
            malformed("payload member {} has the wrong leading shape".format(b), str(lead),
                      str(want))
    nonce = _str_member(data, names, "archive_nonce", malformed)
    if not _NONCE_RE.match(nonce):
        malformed("archive_nonce is not 32 lowercase hex characters", _bounded(repr(nonce)),
                  "32 hex characters")
    for b in blocks:
        if sources[b][0] == "inline":
            continue
        base = os.path.basename(sources[b][0])
        if parse_sidecar_name(base) != (b, nonce):
            malformed("sidecar file name of {} does not carry the archive_nonce".format(b),
                      _bounded(repr(base)), repr(sidecar_name(b, nonce)))
    green_file = _str_member(data, names, "green_file", malformed)
    if not green_file_ok(green_file, os.path.dirname(path) or "."):
        malformed("green_file is not a relative file path", _bounded(repr(green_file)),
                  _GREEN_FILE_RULE)
    if "momentum_convention" not in names:
        malformed("member momentum_convention is missing", "absent", "present")
    from .rpa import check_momentum_marker
    check_momentum_marker(data, path)
    problem = pole_group_problem({k: data[k] for k in names
                                  if isinstance(k, str) and k.startswith("rpa_")})
    if problem is not None:
        check, stored, expected = problem
        malformed(check, _bounded(stored), expected)
    return out_layout, green_file, nvol


_SLOT_SET = "Gamma_bond, chi_w, chibar_w"
_WANTS = (None, "dense", "sz_packed")


def _is_index(v):
    return not isinstance(v, (bool, np.bool_)) and isinstance(v, (int, np.integer))


class SOBondArchive(NamedTuple):
    """A schema-4 spin-orbital bond archive accepted by
    :func:`load_so_bond_archive`.

    ``path`` (as given), ``info`` (the ``fingerprint.SOArchiveInfo`` of the
    compatibility check), ``layout`` (the stored slot layout, ``"dense"`` or
    ``"sz_packed"``), ``bond_archive_layout`` (``"npz"`` or ``"sidecar"``),
    the dimensions ``B``, ``nd``, ``nmat``, ``nvol``, ``sector_map`` (the
    archive's S_z sector map when it is packed, else ``None``),
    ``green_path`` (the companion Green file) and ``sources`` (block ->
    ``("inline", shape)`` or ``(<sidecar path>, shape)``).

    The payload is read only by :meth:`read_slot`, which opens the files
    on every call."""
    path: str
    info: object
    layout: str
    bond_archive_layout: str
    B: int
    nd: int
    nmat: int
    nvol: int
    sector_map: object
    green_path: str
    sources: object

    def _blocks(self, who, name):
        if not (isinstance(name, str) and name in SLOT_PAYLOADS):
            raise ValueError("{}: unknown member {!r} ({})".format(who, name, _SLOT_SET))
        return _block_names(name, self.layout)

    def has_slot(self, name) -> bool:
        """Whether the archive stores the slot payload ``name``
        (``"Gamma_bond"``, ``"chi_w"`` or ``"chibar_w"``)."""
        return all(b in self.sources for b in self._blocks("has_slot", name))

    def read_slot(self, name, *, want=None, l0=None, l1=None):
        """The slot payload ``name`` (``"Gamma_bond"``, ``"chi_w"`` or
        ``"chibar_w"``), frequencies ``[l0, l1)`` (``l0`` None -> 0, ``l1``
        None -> ``nmat``; ``Gamma_bond`` has no frequency axis), in the
        layout ``want``: ``None`` (the stored layout), ``"dense"`` (an
        ndarray ``(..., ND, ND)``) or ``"sz_packed"`` (a
        ``so_bond.PackedArrays``).

        Every call opens the files: an inline member is read whole (it is
        below the small-array bound), a sidecar through ``np.load(...,
        mmap_mode="r")``; the slice is copied. A packed archive read dense
        is expanded with exact zeros off the sectors. A dense archive read
        packed uses the S_z sector map of its bond set and is refused when
        the requested slice has a nonzero off-sector entry, with the
        measured leak of that slice."""
        blocks = self._blocks("read_slot", name)
        if not (want is None or (isinstance(want, str) and want in _WANTS)):
            raise ValueError("read_slot: want must be None, 'dense' or 'sz_packed', got {!r}"
                             .format(want))
        if not all(b in self.sources for b in blocks):
            raise ValueError("'{}' has no {} (this producer does not write it)".format(
                self.path, name))
        if name == "Gamma_bond":
            if l0 is not None or l1 is not None:
                raise ValueError("read_slot: Gamma_bond has no frequency axis")
            sl = None
        else:
            s0 = 0 if l0 is None else l0
            s1 = self.nmat if l1 is None else l1
            if not (_is_index(s0) and _is_index(s1) and 0 <= s0 < s1 <= self.nmat):
                raise ValueError("read_slot: frequency range [{}, {}) outside [0, {})".format(
                    s0, s1, self.nmat))
            sl = slice(int(s0), int(s1))
        arrays = [self._read_block(b, sl) for b in blocks]
        target = self.layout if want is None else want
        if self.layout == "dense":
            if target == "dense":
                return arrays[0]
            return so_bond.PackedArrays.from_dense(arrays[0], self._bond_sector_map(),
                                                   source=self.path)
        packed = so_bond.PackedArrays(self.sector_map, tuple(arrays))
        return packed if target == "sz_packed" else packed.to_dense()

    def _read_block(self, block, sl):
        src = self.sources[block][0]
        if src == "inline":
            with np.load(self.path, allow_pickle=False) as z:
                a = z[block]
            if sl is None:
                return a.astype(np.complex128, copy=False)       # native byte order
            return np.array(a[sl], dtype=np.complex128)
        m = np.load(src, mmap_mode="r", allow_pickle=False)
        try:
            return np.array(m if sl is None else m[sl], dtype=np.complex128)
        finally:
            del m

    def _bond_sector_map(self):
        """The S_z sector map of the archive's bond set (its ``reverse`` and
        ``nd``), for a dense archive read packed."""
        with np.load(self.path, allow_pickle=False) as z:
            rev = np.asarray(z["reverse"])
        return so_bond.sz_sector_map(types.SimpleNamespace(reverse=rev), self.nd)


def load_so_bond_archive(path, *, run_band, run_state, run_interaction, run_bond,
                         expected_producer=PRODUCER_LABEL) -> SOBondArchive:
    """Load the schema-4 spin-orbital bond archive at ``path`` for a run of
    band ``run_band``, state ``run_state``, interaction ``run_interaction``
    and bond set ``run_bond`` (the ``fingerprint`` inputs).

    The steps, each refusal a single-line ``ValueError``; no slot payload is
    read (only npy headers, after every scalar check):

    1. ``np.load(path, allow_pickle=False)``.
    1b. Before any header: the markers ``bond_archive_schema = 4`` and
        ``spin_orbital = True`` (:data:`T_NOT_SO_ARCHIVE`); the 0-d scalars
        ``nmat`` (even, > 0), ``B`` (>= 1), ``nd`` (even, >= 2),
        ``norb_phys`` (>= 1), ``state_nmat``, ``prov_norb_phys``
        (integers), ``beta``, ``state_beta`` (finite floats > 0); the
        tables ``cell_shape``, ``prov_cell_shape`` (integer ``(3,)``, >= 1)
        and ``delta_r`` (integer ``(B', 3)``); then the relations ``nd = 2
        norb_phys``, ``norb_phys = prov_norb_phys``, ``B = len(delta_r)``,
        ``nmat = state_nmat``, ``beta = state_beta``, ``cell_shape =
        prov_cell_shape``. These refusals name ``expected_producer``.
    2. One structural pass, ``fingerprint.validate_structure`` of kind
       ``"bond"`` with the sidecar directory (the archive's).
    3. The schema-4 members: ``index_order``, ``freq_axis``,
       ``bond_archive_layout`` (``Gamma_bond`` always inline; the
       ``chi_w``/``chibar_w`` blocks sidecars exactly in the sidecar
       layout), ``green_layout``, ``scf_converged``, ``types``, the
       leading shapes ``(nvol,)`` and ``(nmat, nvol)`` with ``nvol =
       prod(cell_shape)``, ``archive_nonce`` and the sidecar names
       ``<block>.<archive_nonce>.npy``, ``green_file``
       (:func:`green_file_ok`), ``momentum_convention`` and the pole members
       (all or none, :func:`pole_group_problem`; not exposed).
    4. The companion green file must exist (:data:`T_COMPANION_MISSING`).
    5. The companion: its structure (kind ``"green"``) and its compatibility
       with the run (its ``fp_green`` against its Green function).
    6. The archive's compatibility with the run, its ``fp_green`` against
       the digest of the companion's Green function (an archive whose
       companion was replaced by another Green function is refused here).
    """
    path = os.fspath(path)
    side_dir = os.path.dirname(path) or "."
    with np.load(path, allow_pickle=False) as data:
        names = set(data.files)
        # (1b) the peek, before any header
        _is_so_archive(data, names, path)
        sv = _peek_scalars(data, names, _Malformed(path, expected_producer))
        # (2) one structural pass
        sources = {}
        kind, stage, version = _fp.validate_structure(
            data, path, kind="bond", expected_producer=expected_producer,
            sidecar_dir=side_dir, sources_out=sources)
        label = _fp._readable_label(data, names)
        if label == _fp.UNKNOWN_PRODUCER and expected_producer:
            label = expected_producer
        malformed = _fp._malformed_raiser(path, label)
        sv["layout"] = str(np.asarray(data["block_layout"])[()])
        # (3) the schema-4 members
        out_layout, green_file, nvol = _check_extras(data, names, path, sv, sources,
                                                     malformed)
        # (4) the companion exists
        green_path = os.path.normpath(os.path.join(os.path.dirname(path), green_file))
        if not os.path.isfile(green_path):
            raise _fp.SOArchiveError(T_COMPANION_MISSING.format(green=green_path,
                                                                archive=path))
        run = dict(run_band=run_band, run_state=run_state, nmat_rule="exact",
                   run_interaction=run_interaction, need="so_bond",
                   expected_producer=expected_producer, run_bond=run_bond)
        # (5) the companion's structure and compatibility
        with np.load(green_path, allow_pickle=False) as g:
            g_names = set(g.files)
            g_valid = _fp.validate_structure(g, green_path, kind="green",
                                             expected_producer=expected_producer)
            g_label = _fp._readable_label(g, g_names)
            green_fp = _fp.green_digest(g["green"], GREEN_LAYOUT)
            _fp.check_so_compatibility(g, green_path, validated=g_valid + (g_label,),
                                       green_fp=green_fp, **run)
        # (6) the archive's compatibility, bound to the companion's G
        info = _fp.check_so_compatibility(data, path, validated=(kind, stage, version, label),
                                          green_fp=green_fp, **run)
        sector_map = None
        if sv["layout"] == "sz_packed":
            sector_map = so_bond.sz_sector_map(
                types.SimpleNamespace(reverse=np.asarray(data["reverse"])), sv["nd"])
    return SOBondArchive(path=path, info=info, layout=sv["layout"],
                         bond_archive_layout=out_layout, B=sv["B"], nd=sv["nd"],
                         nmat=sv["nmat"], nvol=nvol, sector_map=sector_map,
                         green_path=green_path, sources=types.MappingProxyType(sources))
