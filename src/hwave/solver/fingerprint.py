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

    yield encode_string("c16")
    yield encode_int(4)
    for n in (nmat, nvol, nd, nd):
        yield encode_int(n)

    slice_bytes = nvol * nd * nd * _C16_BYTES
    if slice_bytes <= CHUNK_BYTES:
        nf = max(1, CHUNK_BYTES // slice_bytes)
        if frequency_last:
            nf = min(nf, _MAX_FREQ_GATHER)
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
