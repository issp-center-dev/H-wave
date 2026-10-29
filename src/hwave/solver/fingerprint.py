"""Canonical compatibility digests for spin-orbital archives.

Each digest is the lowercase-hex SHA-256 of a canonical little-endian byte
stream (header with domain separation, then the component's content).  The
values are hashed as their exact bytes, with no rounding, after ``-0.0`` is
mapped to ``+0.0``; NaN and infinities are refused.  The digests are
run-versus-file compatibility keys, not an integrity layer.
"""
import hashlib
import struct
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
