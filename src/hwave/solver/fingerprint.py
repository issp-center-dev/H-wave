"""Canonical compatibility digests for spin-orbital archives.

Each digest is the lowercase-hex SHA-256 of a canonical little-endian byte
stream (header with domain separation, then the component's content).  The
values are hashed as their exact bytes, with no rounding, after ``-0.0`` is
mapped to ``+0.0``; NaN and infinities are refused.  The digests are
run-versus-file compatibility keys, not an integrity layer.
"""
import functools
import hashlib
import logging
import os
import struct
import subprocess
from typing import Iterable, NamedTuple, Optional

import numpy as np

from . import backend

FINGERPRINT_VERSION = 1
CHUNK_BYTES = 64 * 2**20
_MAX_FREQ_GATHER = 16

PRODUCER_STAGES = (
    "rpa_so",
    "flex_so_onsite",
    "flex_so_bond",
    "sc_so_static",
    "sc_so_dynamic",
    "sc_so_bond_static",
    "sc_so_bond_dynamic",
)

logger = logging.getLogger(__name__)

_MAGIC = "hwave-fp"
_C16_BYTES = 16


# ---------------------------------------------------------------- primitives
def encode_int(v) -> bytes:
    return struct.pack("<q", int(v))


def _real(v, what="value") -> float:
    v = float(v)
    if not np.isfinite(v):
        raise ValueError(f"non-finite {what}: {v!r}")
    return v + 0.0  # -0.0 -> +0.0


def encode_real(v) -> bytes:
    return struct.pack("<d", _real(v))


def encode_complex(v) -> bytes:
    v = complex(v)
    return struct.pack("<dd", _real(v.real, "real part"), _real(v.imag, "imaginary part"))


def encode_bool(v) -> bytes:
    return b"\x01" if v else b"\x00"


def encode_string(s: str) -> bytes:
    b = s.encode("utf-8")
    return struct.pack("<I", len(b)) + b


def _header(component: str) -> bytes:
    return encode_string(_MAGIC) + encode_int(FINGERPRINT_VERSION) + encode_string(component)


def digest(component: str, chunks: Iterable[bytes]) -> str:
    """Hex SHA-256 of the domain-separated header followed by ``chunks``."""
    h = hashlib.sha256()
    h.update(_header(component))
    for c in chunks:
        h.update(c)
    return h.hexdigest()


# -------------------------------------------------------------------- tables
def _table_chunks(component, table):
    if table is None:
        yield encode_bool(False)
        return
    yield encode_bool(True)
    n = len(table)
    yield encode_int(n)
    if n == 0:
        return
    keys = np.array([[*k[0], *k[1]] for k in table], dtype="<i8").reshape(n, 5)
    vals = np.array(list(table.values()), dtype="<c16")
    bad = ~np.isfinite(vals)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        raise ValueError(
            f"non-finite value in {component} entry {tuple(int(x) for x in keys[i])}"
        )
    order = np.lexsort(keys.T[::-1])  # ascending by (Rx, Ry, Rz, a, b)
    rec = np.empty(n, dtype=[("k", "<i8", (5,)), ("v", "<f8", (2,))])
    rec["k"] = keys[order]
    rec["v"] = vals[order].view("<f8").reshape(n, 2) + 0.0
    rows = max(1, CHUNK_BYTES // rec.dtype.itemsize)
    for a in range(0, n, rows):
        yield np.ascontiguousarray(rec[a:a + rows]).view(np.uint8).data


def table_digest(component: str, table: Optional[dict]) -> str:
    """Digest of a ``wan90.read_w90``-style table ``{(irvec, orbvec): complex}``.

    ``None`` (absent) encodes only ``False``; a present table, even an empty
    one, encodes ``True``, the count and the sorted entries.
    """
    return digest(component, _table_chunks(component, table))


# ---------------------------------------------------------------------- band
class BandMeta(NamedTuple):
    coeff_extern: float
    enable_spin_orbital: bool
    index_convention: str
    norb: int
    norb_phys: int
    cell_shape: tuple
    sub_shape: tuple


class BandInputs(NamedTuple):
    digests: dict
    transfer: dict
    extern: Optional[dict]
    meta: BandMeta


def _meta_chunks(meta: BandMeta):
    yield encode_real(meta.coeff_extern)
    yield encode_bool(meta.enable_spin_orbital)
    yield encode_string(meta.index_convention)
    yield encode_int(meta.norb)
    yield encode_int(meta.norb_phys)
    for v in meta.cell_shape:
        yield encode_int(v)
    for v in meta.sub_shape:
        yield encode_int(v)


def band_digests(transfer: dict, extern: Optional[dict], meta: BandMeta) -> dict:
    fp_transfer = table_digest("transfer", transfer)
    fp_extern = table_digest("extern", extern)
    fp_meta = digest("meta", _meta_chunks(meta))
    fp_band = digest("band", (encode_string(s) for s in (fp_transfer, fp_extern, fp_meta)))
    return {"fp_transfer": fp_transfer, "fp_extern": fp_extern,
            "fp_meta": fp_meta, "fp_band": fp_band}


def band_inputs(transfer, extern, meta) -> BandInputs:
    return BandInputs(band_digests(transfer, extern, meta), transfer, extern, meta)


# --------------------------------------------------------------------- green
def _canonical_bytes(block, what):
    f = np.ascontiguousarray(block).astype("<c16", copy=True).view("<f8")
    f += 0.0  # -0.0 -> +0.0
    if not np.isfinite(f).all():
        raise ValueError(f"non-finite value in {what}")
    return memoryview(f)


def _green_chunks(green, layout):
    shape = tuple(green.shape)
    if layout == "b_f_k_o_o_v1":
        if len(shape) != 5 or shape[0] != 1:
            raise ValueError(f"layout {layout!r} needs shape (1, Nmat, nvol, nd, nd), got {shape}")
        nmat, nvol, nd = shape[1], shape[2], shape[3]
        view = green[0]
        frequency_last = False
    elif layout == "sc_o_o_xyz_f":
        if len(shape) != 6:
            raise ValueError(f"layout {layout!r} needs shape (nd, nd, Nx, Ny, Nz, Nmat), got {shape}")
        nd, nmat = shape[0], shape[5]
        nvol = shape[2] * shape[3] * shape[4]
        view = green
        frequency_last = True
    else:
        raise ValueError(f"unknown green layout {layout!r}")

    def block(f0, f1, k0, k1):
        if frequency_last:  # (nd, nd, Nx, Ny, Nz, f) -> (f, k, nd, nd)
            if k0 == 0 and k1 == nvol:
                b = backend.to_host(view[..., f0:f1])
            else:  # one frequency, a k sub-range gathered by index
                xp = backend.array_module_of(view)
                ix = np.unravel_index(np.arange(k0, k1), shape[2:5])
                idx = tuple(xp.asarray(i) for i in ix)
                b = backend.to_host(view[(slice(None), slice(None)) + idx + (slice(f0, f1),)])
            b = b.reshape(nd, nd, -1, f1 - f0)
            return np.ascontiguousarray(b.transpose(3, 2, 0, 1))
        return backend.to_host(view[f0:f1, k0:k1])

    return _block_chunks(nmat, nvol, nd, block,
                         _MAX_FREQ_GATHER if frequency_last else None)


def _block_chunks(nmat, nvol, nd, block, max_freq=None):
    """The canonical green stream (header, then the ``(Nmat, nvol, nd, nd)``
    values in C order) drawn from ``block(f0, f1, k0, k1)``, which returns
    the ``[f0:f1, k0:k1]`` block as a host array of shape
    ``(f1 - f0, k1 - k0, nd, nd)``. Each requested block holds at most
    ``max(CHUNK_BYTES, one k-row of nd*nd complex)``."""
    yield encode_string("c16")
    yield encode_int(4)
    for n in (nmat, nvol, nd, nd):
        yield encode_int(n)

    slice_bytes = nvol * nd * nd * _C16_BYTES
    if slice_bytes <= CHUNK_BYTES:
        nf = max(1, CHUNK_BYTES // slice_bytes)
        if max_freq is not None:
            nf = min(nf, max_freq)
        for f0 in range(0, nmat, nf):
            yield _canonical_bytes(block(f0, min(f0 + nf, nmat), 0, nvol), "green")
    else:
        nk = max(1, CHUNK_BYTES // (nd * nd * _C16_BYTES))
        for f0 in range(nmat):
            for k0 in range(0, nvol, nk):
                yield _canonical_bytes(block(f0, f0 + 1, k0, min(k0 + nk, nvol)), "green")


def green_digest(green, layout: str) -> str:
    """Digest of the canonical ``(Nmat, nvol, nd, nd)`` complex Green array.

    The array is hashed block by block (a CuPy array is moved to the host one
    block at a time); the caller's array is never modified.
    """
    return digest("green", _green_chunks(green, layout))


def green_digest_blocks(nmat: int, nvol: int, nd: int, block) -> str:
    """:func:`green_digest` (layout ``"b_f_k_o_o_v1"``) of a canonical
    ``(1, Nmat, nvol, nd, nd)`` Green array that is never materialized:
    ``block(f0, f1, k0, k1)`` returns its ``[0, f0:f1, k0:k1]`` block
    (shape ``(f1 - f0, k1 - k0, nd, nd)``, host or device). Blocks are
    requested in order and each holds at most ``max(CHUNK_BYTES, one k-row
    of nd*nd complex)``, so a caller that builds the array lazily (a sum, a
    layout expansion) never holds more than one block of it. The digest is
    byte-identical to ``green_digest`` of the materialized array."""
    def host_block(f0, f1, k0, k1):
        b = backend.to_host(block(f0, f1, k0, k1))
        if tuple(b.shape) != (f1 - f0, k1 - k0, nd, nd):
            raise ValueError("green block [{}:{}, {}:{}] has shape {}, expected {}".format(
                f0, f1, k0, k1, tuple(b.shape), (f1 - f0, k1 - k0, nd, nd)))
        return b
    return digest("green", _block_chunks(int(nmat), int(nvol), int(nd), host_block))


# ---------------------------------------------------------------- provenance
class SOState(NamedTuple):
    """The realized normal state of a spin-orbital producer.

    A tagged union on ``calc_mu``:

    * ``calc_mu = True``: ``target`` is the electron density per cell summed
      over all spin-orbitals, and ``mu_tol = |density - target|`` (the
      achieved residual).
    * ``calc_mu = False``: ``target`` is the configured mu (energy units),
      ``mu_tol = 0.0``, and ``density`` (electrons per cell over all
      spin-orbitals) is a diagnostic only.
    """
    beta: float
    nmat: int
    mu: float
    density: float
    calc_mu: bool
    target: float
    mu_tol: float


def so_state(*, beta, nmat, mu, density, calc_mu, target) -> SOState:
    """An :class:`SOState` with ``mu_tol`` set by the ``calc_mu`` tag."""
    calc_mu = bool(calc_mu)
    density = float(density)
    target = float(target)
    mu_tol = abs(density - target) if calc_mu else 0.0
    return SOState(float(beta), int(nmat), float(mu), density, calc_mu, target, mu_tol)


_VERSION_MAX = 256

#: The version-1 members ``so_provenance_members`` writes (the schema table
#: of the design's archive section), in a fixed order. Consumers validate an
#: archive by THIS list, never by a name prefix such as ``state_``: FLEX
#: archives already carry an unrelated ``state_iteration`` member.
SO_PROVENANCE_MEMBERS = (
    "fingerprint_version", "producer_stage", "producer_label", "producer_version",
    "fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green",
    "prov_transfer_keys", "prov_transfer_vals", "prov_extern_present",
    "prov_extern_keys", "prov_extern_vals",
    "prov_coeff_extern", "prov_enable_spin_orbital", "prov_index_convention",
    "prov_norb", "prov_norb_phys", "prov_cell_shape", "prov_sub_shape",
    "state_beta", "state_nmat", "state_mu", "state_density", "state_calc_mu",
    "state_target", "state_mu_tol",
)


@functools.lru_cache(maxsize=None)
def _git_describe() -> Optional[str]:
    """``git describe --always --dirty`` of the repository the package lives
    in (2 s limit per git call), or None. Accepted only when the package
    really is that repository's ``src/hwave`` (the top level holds a
    ``pyproject.toml``), so an installed copy inside some unrelated
    repository never borrows its description. Cached per process; never
    raises."""
    try:
        pkg = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        r = subprocess.run(["git", "rev-parse", "--show-toplevel"], cwd=pkg,
                           capture_output=True, text=True, timeout=2.0)
        top = r.stdout.strip()
        if r.returncode != 0 or not top:
            return None
        if not (os.path.isfile(os.path.join(top, "pyproject.toml"))
                and os.path.realpath(os.path.join(top, "src", "hwave"))
                == os.path.realpath(pkg)):
            return None
        r = subprocess.run(["git", "describe", "--always", "--dirty"], cwd=top,
                           capture_output=True, text=True, timeout=2.0)
        out = r.stdout.strip()
        return out if r.returncode == 0 and out else None
    except Exception:
        return None


def _producer_version() -> str:
    """``hwave.__version__`` (read at call time) plus `` (<git describe>)``
    when available; longer than 256 characters it is cut to 255 plus
    ``…``. Never fails."""
    try:
        import hwave
        v = str(getattr(hwave, "__version__", "unknown"))
    except Exception:
        v = "unknown"
    g = _git_describe()
    if g:
        v = "{} ({})".format(v, g)
    if len(v) > _VERSION_MAX:
        v = v[:_VERSION_MAX - 1] + "…"
        logger.info("producer_version truncated to %d characters", _VERSION_MAX)
    return v


def table_arrays(table: Optional[dict]):
    """``(keys int64 (n, 5), vals complex128 (n,))`` of a ``read_w90``-style
    table, sorted by the signed ``(Rx, Ry, Rz, a, b)`` key (the digest
    order); empty arrays for ``None``."""
    if not table:
        return np.zeros((0, 5), dtype=np.int64), np.zeros((0,), dtype=np.complex128)
    keys = np.array([[*k[0], *k[1]] for k in table], dtype=np.int64).reshape(len(table), 5)
    vals = np.array(list(table.values()), dtype=np.complex128)
    order = np.lexsort(keys.T[::-1])
    return np.ascontiguousarray(keys[order]), np.ascontiguousarray(vals[order])


def so_provenance_members(*, producer_stage, producer_label, transfer, extern,
                          meta: BandMeta, fp_green: str, state: SOState) -> dict:
    """The version-1 fingerprint and state members of a spin-orbital archive.

    Build it once per run and merge the same dict into every archive of
    that run. Refuses (ValueError naming the member) any non-finite float
    or complex member.
    """
    if producer_stage not in PRODUCER_STAGES:
        raise ValueError("unknown producer_stage {!r}".format(producer_stage))
    if not (isinstance(fp_green, str) and len(fp_green) == 64
            and all(c in "0123456789abcdef" for c in fp_green)):
        raise ValueError("fp_green must be a 64-character lowercase hex digest, got {!r}"
                         .format(fp_green))
    tk, tv = table_arrays(transfer)
    ek, ev = table_arrays(extern)
    m = {
        "fingerprint_version": np.int64(FINGERPRINT_VERSION),
        "producer_stage": np.str_(producer_stage),
        "producer_label": np.str_(producer_label),
        "producer_version": np.str_(_producer_version()),
        "fp_green": np.str_(fp_green),
        "prov_transfer_keys": tk,
        "prov_transfer_vals": tv,
        "prov_extern_present": np.bool_(extern is not None),
        "prov_extern_keys": ek,
        "prov_extern_vals": ev,
        "prov_coeff_extern": np.float64(meta.coeff_extern),
        "prov_enable_spin_orbital": np.bool_(meta.enable_spin_orbital),
        "prov_index_convention": np.str_(meta.index_convention),
        "prov_norb": np.int64(meta.norb),
        "prov_norb_phys": np.int64(meta.norb_phys),
        "prov_cell_shape": np.array(meta.cell_shape, dtype=np.int64).reshape(3),
        "prov_sub_shape": np.array(meta.sub_shape, dtype=np.int64).reshape(3),
        "state_beta": np.float64(state.beta),
        "state_nmat": np.int64(state.nmat),
        "state_mu": np.float64(state.mu),
        "state_density": np.float64(state.density),
        "state_calc_mu": np.bool_(state.calc_mu),
        "state_target": np.float64(state.target),
        "state_mu_tol": np.float64(state.mu_tol),
    }
    for k, v in m.items():
        if v.dtype.kind in "fc" and not np.all(np.isfinite(v)):
            raise ValueError("non-finite spin-orbital provenance member {}: {!r}".format(k, v))
    if not state.beta > 0:
        raise ValueError("state_beta must be > 0, got {!r}".format(state.beta))
    if not state.calc_mu and state.mu_tol != 0.0:
        raise ValueError("state_mu_tol must be 0.0 for a fixed-mu state, got {!r}"
                         .format(state.mu_tol))
    fps = band_digests(transfer, extern, meta)
    for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band"):
        m[k] = np.str_(fps[k])
    if set(m) != set(SO_PROVENANCE_MEMBERS):  # pragma: no cover - schema drift guard
        raise RuntimeError("spin-orbital provenance members drifted from the schema")
    return {k: m[k] for k in SO_PROVENANCE_MEMBERS}


# ---------------------------------------------------------- consumer checks
#: dtype kind and ndim of every version-1 member (the schema table of the
#: design's archive section); validated by :func:`check_so_archive`.
SO_MEMBER_SCHEMA = {
    "fingerprint_version": ("i", 0),
    "producer_stage": ("U", 0), "producer_label": ("U", 0), "producer_version": ("U", 0),
    "fp_transfer": ("U", 0), "fp_extern": ("U", 0), "fp_meta": ("U", 0),
    "fp_band": ("U", 0), "fp_green": ("U", 0),
    "prov_transfer_keys": ("i", 2), "prov_transfer_vals": ("c", 1),
    "prov_extern_present": ("b", 0),
    "prov_extern_keys": ("i", 2), "prov_extern_vals": ("c", 1),
    "prov_coeff_extern": ("f", 0), "prov_enable_spin_orbital": ("b", 0),
    "prov_index_convention": ("U", 0),
    "prov_norb": ("i", 0), "prov_norb_phys": ("i", 0),
    "prov_cell_shape": ("i", 1), "prov_sub_shape": ("i", 1),
    "state_beta": ("f", 0), "state_nmat": ("i", 0), "state_mu": ("f", 0),
    "state_density": ("f", 0), "state_calc_mu": ("b", 0),
    "state_target": ("f", 0), "state_mu_tol": ("f", 0),
}

UNKNOWN_PRODUCER = "<unknown producer>"

#: the consumer-check messages go to the ``qlms`` hierarchy, beside the
#: consumers' own log lines
check_logger = logging.getLogger("qlms").getChild("solver").getChild("fingerprint")

#: the relative tolerances of the compatibility checks
BETA_RTOL = 1e-12
TARGET_RTOL = 1e-12
MU_RTOL = 1e-12
DENSITY_RTOL = 1e-8
#: the cap on a producer's density-closure residual, relative to max(1, |t|)
RESIDUAL_CAP = 1e-8
#: |state_mu_tol - |state_density - state_target|| allowed, relative
MU_TOL_CONSISTENCY = 1e-12

NMAT_RULES = ("exact", "load")


class SOArchiveError(ValueError):
    """A consumed spin-orbital archive refused by :func:`check_so_archive`
    (a ValueError; callers that rewrap layout errors let it through)."""


class SOArchiveInfo(NamedTuple):
    """What :func:`check_so_archive` learned about a consumed archive.

    Validation and message data only: the hwave_sc writers record their own
    state and never copy these fields. A legacy (unversioned) archive gives
    ``SOArchiveInfo(None, None, None, legacy=True)``."""
    producer_state: Optional[SOState]
    producer_fp_green: Optional[str]
    producer_label: Optional[str]
    legacy: bool


def _members(data):
    return set(data.files) if hasattr(data, "files") else set(data)


def _readable_label(data, names):
    if "producer_label" not in names:
        return UNKNOWN_PRODUCER
    try:
        v = np.asarray(data["producer_label"])
        if v.ndim == 0 and v.dtype.kind == "U" and str(v[()]).strip():
            return str(v[()])
    except Exception:
        pass
    return UNKNOWN_PRODUCER


def _table_from_arrays(keys, vals):
    return {(tuple(int(x) for x in k[:3]), tuple(int(x) for x in k[3:])): complex(v)
            for k, v in zip(keys, vals)}


def _sorted_key(k):
    return tuple(k[0]) + tuple(k[1])


def _fmt_value(v):
    return "absent" if v is None else repr(v)


def _table_diff(name, file_tbl, run_tbl):
    """``<Table>: file <n> rows, run <m> rows; first differing row (R=..,
    a=.., b=..): file <v>, run <w>`` (1-based orbitals)."""
    file_tbl = file_tbl or {}
    run_tbl = run_tbl or {}
    head = "{}: file {} rows, run {} rows".format(name, len(file_tbl), len(run_tbl))
    for k in sorted(set(file_tbl) | set(run_tbl), key=_sorted_key):
        a, b = file_tbl.get(k), run_tbl.get(k)
        if a is None or b is None or complex(a) != complex(b):
            return "{}; first differing row (R={}, a={}, b={}): file {}, run {}".format(
                head, tuple(int(x) for x in k[0]), int(k[1][0]) + 1, int(k[1][1]) + 1,
                _fmt_value(a), _fmt_value(b))
    return head + " (the rows agree; the digests differ)"


def _scalar(data, k):
    return np.asarray(data[k])[()]


def check_so_archive(data, path: str, *, run_band: BandInputs, run_state: SOState,
                     run_has_extern: bool, nmat_rule: str, green=None,
                     green_layout: Optional[str] = None,
                     green_fp: Optional[str] = None) -> SOArchiveInfo:
    """Validate the version-1 provenance of a spin-orbital archive against
    this run (band, state and Green function only; no interaction or
    topology digest is recorded in this version).

    ``data`` is the open npz (or a mapping of its members), ``path`` its
    file name for the messages. ``run_band`` is this run's band
    (:func:`band_inputs`), ``run_state`` its state: ``beta``, ``nmat``,
    ``calc_mu`` and ``target`` (and ``mu`` for a fixed-mu run) are compared;
    ``density`` and ``mu_tol`` are not read. ``nmat_rule`` is ``"exact"``
    (a borrowed G) or ``"load"`` (load mode: the existing stored-axis rule
    applies instead, Nmat is not compared here). With ``green`` (in
    ``green_layout``) or ``green_fp`` (its digest, already computed by the
    caller) the archive's ``fp_green`` is verified; otherwise it is a
    producer assertion, logged at INFO.

    Order: (0) an archive without ``fingerprint_version`` takes the legacy
    policy -- refused when the run has an Extern field, a WARNING otherwise;
    (1) structural validation; (2) band; (3) state; (4) green. Every refusal
    is a single-line ValueError."""
    if nmat_rule not in NMAT_RULES:
        raise ValueError("nmat_rule must be one of {}, got {!r}".format(NMAT_RULES, nmat_rule))
    names = _members(data)
    label = _readable_label(data, names)

    # (0) legacy first: never reported as malformed
    if "fingerprint_version" not in names:
        why = ("written before this version, or by a run whose one-body band or "
               "state the version-1 fingerprints do not cover: trans_mod/"
               "green_init, or a supplied chi0q")
        if run_has_extern:
            raise SOArchiveError(
                "'{}' has no fingerprint_version ({}); the run has an Extern "
                "field, so its band cannot be verified. Regenerate it.".format(path, why))
        check_logger.warning(
            "'%s' has no fingerprint_version (%s); band and state compatibility "
            "are not checked. Regenerate it to enable the checks.", path, why)
        return SOArchiveInfo(None, None, None, legacy=True)

    def malformed(check, stored, expected):
        raise SOArchiveError(
            "'{}' is malformed or truncated: {} ({} vs {}). Regenerate it with "
            "{}.".format(path, check, stored, expected, label))

    # (1) structural
    ver = np.asarray(data["fingerprint_version"])
    if ver.ndim != 0 or ver.dtype.kind not in "iu":
        malformed("fingerprint_version is not an integer scalar",
                  "{}/{}d".format(ver.dtype, ver.ndim), "int64/0d")
    if int(ver[()]) != FINGERPRINT_VERSION:
        malformed("unsupported fingerprint_version", int(ver[()]), FINGERPRINT_VERSION)
    for k in SO_PROVENANCE_MEMBERS:
        if k not in names:
            malformed("member {} is missing".format(k), "absent", "present")
        kind, ndim = SO_MEMBER_SCHEMA[k]
        a = np.asarray(data[k])
        if a.dtype.kind != kind or a.ndim != ndim:
            malformed("member {} has the wrong dtype or rank".format(k),
                      "{}/{}d".format(a.dtype, a.ndim),
                      "{}/{}d".format({"i": "int64", "U": "str", "c": "complex128",
                                       "f": "float64", "b": "bool"}[kind], ndim))
        if kind in "fc" and not np.all(np.isfinite(a)):
            malformed("member {} is not finite".format(k), repr(a.tolist()), "finite")
    stage = str(_scalar(data, "producer_stage"))
    if stage not in PRODUCER_STAGES:
        malformed("producer_stage is not a known stage", repr(stage),
                  "one of {}".format(", ".join(PRODUCER_STAGES)))
    for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green"):
        s = str(_scalar(data, k))
        if len(s) != 64 or any(c not in "0123456789abcdef" for c in s):
            malformed("{} is not a 64-character lowercase hex digest".format(k),
                      repr(s), "64 hex characters")
    tk, tv = np.asarray(data["prov_transfer_keys"]), np.asarray(data["prov_transfer_vals"])
    ek, ev = np.asarray(data["prov_extern_keys"]), np.asarray(data["prov_extern_vals"])
    for nm, keys, vals in (("transfer", tk, tv), ("extern", ek, ev)):
        if keys.shape[1:] != (5,) or vals.shape != (keys.shape[0],):
            malformed("prov_{}_keys/vals shapes".format(nm),
                      "{}/{}".format(keys.shape, vals.shape), "(n, 5)/(n,)")
    for k in ("prov_cell_shape", "prov_sub_shape"):
        if np.asarray(data[k]).shape != (3,):
            malformed("{} shape".format(k), np.asarray(data[k]).shape, (3,))
    ext_present = bool(_scalar(data, "prov_extern_present"))
    if not ext_present and ek.shape[0] != 0:
        malformed("prov_extern_keys of an absent Extern", "{} rows".format(ek.shape[0]),
                  "0 rows")
    beta_f = float(_scalar(data, "state_beta"))
    if not beta_f > 0:
        malformed("state_beta must be > 0", repr(beta_f), "> 0")
    if "beta" in names:
        b = np.asarray(data["beta"])
        if b.ndim != 0 or b.dtype.kind not in "fi" or float(b[()]) != beta_f:
            malformed("beta differs from state_beta", repr(b.tolist()), repr(beta_f))
    nmat_f = int(_scalar(data, "state_nmat"))
    if "nmat" in names:
        n = np.asarray(data["nmat"])
        if n.ndim != 0 or n.dtype.kind not in "iu" or int(n[()]) != nmat_f:
            malformed("nmat differs from state_nmat", repr(n.tolist()), nmat_f)
    # the diff source must reproduce the stored digests
    file_transfer = _table_from_arrays(tk, tv)
    file_extern = _table_from_arrays(ek, ev) if ext_present else None
    file_meta = BandMeta(
        float(_scalar(data, "prov_coeff_extern")),
        bool(_scalar(data, "prov_enable_spin_orbital")),
        str(_scalar(data, "prov_index_convention")),
        int(_scalar(data, "prov_norb")), int(_scalar(data, "prov_norb_phys")),
        tuple(int(x) for x in np.asarray(data["prov_cell_shape"])),
        tuple(int(x) for x in np.asarray(data["prov_sub_shape"])))
    file_fps = {k: str(_scalar(data, k))
                for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band")}
    recomputed = band_digests(file_transfer, file_extern, file_meta)
    for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band"):
        if recomputed[k] != file_fps[k]:
            malformed("{} does not match the digest recomputed from its prov_* "
                      "fields".format(k), file_fps[k][:12] + "...",
                      recomputed[k][:12] + "...")
    state = SOState(beta_f, nmat_f, float(_scalar(data, "state_mu")),
                    float(_scalar(data, "state_density")),
                    bool(_scalar(data, "state_calc_mu")),
                    float(_scalar(data, "state_target")),
                    float(_scalar(data, "state_mu_tol")))
    t = state.target
    scale = max(1.0, abs(t))
    if state.calc_mu:
        resid = abs(state.density - t)
        if abs(state.mu_tol - resid) > MU_TOL_CONSISTENCY * scale:
            malformed("state_mu_tol is not |state_density - state_target|",
                      repr(state.mu_tol), repr(resid))
        if state.mu_tol > RESIDUAL_CAP * scale:
            raise SOArchiveError(
                "'{}' cannot be reused: the producer did not reach its filling "
                "(residual {:.3e} electrons per cell > 1e-8·max(1,|t|) = {:.3e}). "
                "Regenerate it with {} so that its density closes to that "
                "tolerance.".format(path, state.mu_tol, RESIDUAL_CAP * scale, label))
    elif state.mu_tol != 0.0:
        malformed("state_mu_tol of a fixed-mu state", repr(state.mu_tol), "0.0")

    # (2) band
    run_fps = run_band.digests
    diffs = []
    if file_fps["fp_meta"] != run_fps["fp_meta"]:
        for field in BandMeta._fields:
            a, b = getattr(file_meta, field), getattr(run_band.meta, field)
            if isinstance(b, (list, tuple)):
                b = tuple(int(x) for x in b)
            if a != b:
                diffs.append("{}: file {!r}, run {!r}".format(field, a, b))
        if not diffs:
            diffs.append("meta: file {}, run {}".format(file_fps["fp_meta"][:12],
                                                       run_fps["fp_meta"][:12]))
    if file_fps["fp_transfer"] != run_fps["fp_transfer"]:
        diffs.append(_table_diff("Transfer", file_transfer, run_band.transfer))
    if file_fps["fp_extern"] != run_fps["fp_extern"]:
        run_present = run_band.extern is not None
        if ext_present != run_present:
            diffs.append("Extern: file {}, run {}".format(
                "present" if ext_present else "absent",
                "present" if run_present else "absent"))
        else:
            diffs.append(_table_diff("Extern", file_extern, run_band.extern))
    if diffs:
        raise SOArchiveError(
            "'{}' (written by {}) does not match this run: band differs: {}. "
            "Archives are reused only at the identical state point (not a warm "
            "start across a sweep); re-run {} with this input (the same "
            "Transfer/Extern/coeff_extern).".format(path, label, "; ".join(diffs), label))

    # (3) state
    diffs = []
    if state.calc_mu != bool(run_state.calc_mu):
        diffs.append("calc_mu: file {}, run {}".format(state.calc_mu,
                                                       bool(run_state.calc_mu)))
    elif state.calc_mu:
        tr = float(run_state.target)
        sr = max(1.0, abs(tr))
        if abs(state.target - tr) > TARGET_RTOL * sr:
            diffs.append("target {!r} (file) vs {!r} (run) electrons per cell".format(
                state.target, tr))
        tol = max(DENSITY_RTOL * sr, 10.0 * state.mu_tol)
        if abs(state.density - tr) > tol:
            diffs.append("density {!r} (file) vs target {!r} (run, tol {:.3g}) "
                         "electrons per cell".format(state.density, tr, tol))
    else:
        mr = float(run_state.mu)
        if not abs(state.mu - mr) <= MU_RTOL * max(1.0, abs(mr)):
            diffs.append("fixed mu {!r} (file) vs {!r} (run) energy units".format(
                state.mu, mr))
    br = float(run_state.beta)
    if not abs(state.beta - br) <= BETA_RTOL * br:
        diffs.append("beta {!r} (file) vs {!r} (run)".format(state.beta, br))
    if nmat_rule == "exact" and state.nmat != int(run_state.nmat):
        diffs.append("nmat {} (file) vs {} (run)".format(state.nmat, int(run_state.nmat)))
    if diffs:
        if len(diffs) == 1 and diffs[0].startswith("nmat "):
            fix = "set [mode.param] Nmat = {} (the producer's grid)".format(state.nmat)
        else:
            fix = "re-run {} at this filling and temperature".format(label)
        raise SOArchiveError(
            "'{}' (written by {}) does not match this run: state differs: {}. "
            "Archives are reused only at the identical state point; {}.".format(
                path, label, "; ".join(diffs), fix))

    # (4) green
    file_green = str(_scalar(data, "fp_green"))
    if green_fp is None and green is not None:
        green_fp = green_digest(green, green_layout)
    if green_fp is not None:
        if green_fp != file_green:
            raise SOArchiveError(
                "'{}' (written by {}): fp_green mismatch: the Green function this "
                "run loaded has digest {}..., the archive records {}...; the Green "
                "function and the susceptibility do not come from the same run. "
                "Regenerate both with {}.".format(path, label, green_fp[:12],
                                                 file_green[:12], label))
    else:
        check_logger.info("fp_green: producer assertion (%s)", path)
    return SOArchiveInfo(state, file_green, label, legacy=False)
