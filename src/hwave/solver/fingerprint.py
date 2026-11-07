"""Canonical compatibility digests for spin-orbital archives.

Each digest is the lowercase-hex SHA-256 of a canonical little-endian byte
stream (header with domain separation, then the component's content).  The
values are hashed as their exact bytes, with no rounding, after ``-0.0`` is
mapped to ``+0.0``; NaN and infinities are refused.  The digests are
run-versus-file compatibility keys, not an integrity layer.

Every spin-orbital archive is written at :data:`FINGERPRINT_VERSION` (2):
the band, Green and state fingerprints plus the interaction fingerprints,
the capability and the resolved bond topology. A consumer validates an
archive in three tables before any payload is read: the provenance schema
by archive kind and producer stage (:data:`PROVENANCE_SCHEMA`,
:func:`required_members`), the payload names by kind
(:data:`PAYLOAD_MEMBERS`) and the payload metadata read from the npy headers
alone (:data:`PAYLOAD_RANKS`). A version-1 or unversioned archive is
refused. See :func:`validate_structure` and :func:`check_so_archive`.
"""
import functools
import hashlib
import logging
import os
import struct
import subprocess
import types
import zipfile
from typing import Iterable, NamedTuple, Optional

import numpy as np

from . import backend

#: the fingerprint version every spin-orbital archive is written at, and
#: the only one a consumer reads
FINGERPRINT_VERSION = 2
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

#: ``fp_capability`` values: "so_onsite" = not bond-resolved, "so_bond" =
#: bond-resolved (literal strings, not digests)
CAPABILITIES = ("so_onsite", "so_bond")
#: the literal ``fp_resolved_topology`` of a stage that is not bond-resolved
ONSITE_RESOLVED_SENTINEL = "none:onsite"
#: the capability each producer stage writes
STAGE_CAPABILITY = {
    "rpa_so": "so_onsite",
    "flex_so_onsite": "so_onsite",
    "flex_so_bond": "so_bond",
    "sc_so_static": "so_onsite",
    "sc_so_dynamic": "so_onsite",
    "sc_so_bond_static": "so_bond",
    "sc_so_bond_dynamic": "so_bond",
}

#: the archive kinds and the producer stages this version reads for each
ARTIFACT_KINDS = ("chi0q", "chiq", "green", "sigma", "gap")
KIND_STAGES = {
    "chi0q": ("rpa_so", "flex_so_onsite"),
    "chiq": ("rpa_so", "flex_so_onsite"),
    "green": ("flex_so_onsite",),
    "sigma": ("flex_so_onsite",),
    "gap": ("sc_so_static", "sc_so_dynamic", "sc_so_bond_static"),
}
#: stages that are recognised but not read yet, with the path that brings them
DEFERRED_STAGES = {
    "sc_so_bond_dynamic": "the dynamic spin-orbital bond path",
    "flex_so_bond": "the spin-orbital bond FLEX",
}
#: the payload members each kind must carry
PAYLOAD_MEMBERS = {
    "chi0q": ("chi0q",),
    "chiq": ("chiq",),
    "green": ("green", "green_layout"),
    "sigma": ("sigma_static",),
    "gap": ("gap",),
}
#: the rank(s) of the complex payload of each (kind, stage), read from the
#: npy header. A spin-orbital RPA writes its chi0q in every layout of the
#: band it reduces to -- the spinful general (rank 6) and reduced (rank 4)
#: layouts, and the physical-orbital spin-diagonal ones with a leading spin
#: axis (ranks 7 and 5) -- and its chiq at rank 6 or 4; the layout refusal of
#: the hwave_sc loader stays with the loader. A dynamic gap carries the
#: frequency axis.
PAYLOAD_RANKS = {
    ("chi0q", "rpa_so"): (4, 5, 6, 7), ("chi0q", "flex_so_onsite"): (4, 6),
    ("chiq", "rpa_so"): (4, 6), ("chiq", "flex_so_onsite"): (6,),
    ("green", "flex_so_onsite"): (5,),
    ("sigma", "flex_so_onsite"): (5,),
    ("gap", "sc_so_static"): (5,), ("gap", "sc_so_dynamic"): (6,),
    ("gap", "sc_so_bond_static"): (5,),
}
#: the layout a green archive's ``green_layout`` names
GREEN_LAYOUT = "b_f_k_o_o_v1"
#: payload members that are 0-d strings, not complex arrays: checked by
#: value, never by their npy header
_STRING_PAYLOADS = frozenset({"green_layout"})

#: the layout members of a bond-resolved archive: ``block_layout`` (the slot
#: layout, ``"dense"`` or ``"sz_packed"``) and the three sector-map members a
#: packed slot archive carries (``so_bond._MAP_MEMBERS``)
LAYOUT_MEMBERS = ("block_layout", "sz_sector_order", "sz_sector_offsets", "sz_sector_index")
#: the reserved member-name prefix: a member whose name starts with it is
#: refused wherever the sector-map members are not allowed, and a packed slot
#: archive admits only the three sector-map members under it. Every other
#: member not named by a rule is accepted and ignored (forward compatible).
RESERVED_PREFIX = "sz_sector_"
_SLOT_LAYOUTS = ("dense", "sz_packed")


class SlotRules(NamedTuple):
    """The archive kinds whose payloads are code-slot arrays in a
    ``block_layout`` (dense ``<p>``, or packed ``<p>@0``, ``<p>@+1``,
    ``<p>@-1``): ``slot_kinds``; per kind the required ``payloads`` and the
    ``optional`` payloads; ``ranks[(kind, stage)][payload]`` the admitted
    ranks of each payload block."""
    slot_kinds: frozenset
    payloads: dict
    optional: dict
    ranks: dict


#: the slot kinds this version reads: none (the bond-resolved slot archives
#: arrive with the stages of :data:`DEFERRED_STAGES`)
PRODUCTION_SLOT_RULES = SlotRules(frozenset(), {}, {}, {})

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


def _header(component: str, version: int = FINGERPRINT_VERSION) -> bytes:
    return encode_string(_MAGIC) + encode_int(version) + encode_string(component)


def digest(component: str, chunks: Iterable[bytes], *, version: int = FINGERPRINT_VERSION) -> str:
    """Hex SHA-256 of the domain-separated header followed by ``chunks``."""
    h = hashlib.sha256()
    h.update(_header(component, version))
    for c in chunks:
        h.update(c)
    return h.hexdigest()


# ------------------------------------------------- interaction / topology
def _record_chunks(records):
    n = len(records)
    yield encode_int(n)
    if n == 0:
        return
    types = np.asarray(records["type"])
    tbytes = [str(t).encode("utf-8") for t in types]
    modes = np.asarray(records["modes"], dtype="<i8").reshape(n, 4)
    cells = np.asarray(records["cells"], dtype="<i8").reshape(n, 12)
    v = np.asarray(records["v"], dtype="<c16")
    bad = ~np.isfinite(v)
    if bad.any():
        i = int(np.flatnonzero(bad)[0])
        key = (str(types[i]), tuple(int(x) for x in modes[i]),
               tuple(tuple(int(x) for x in row) for row in cells[i].reshape(4, 3)))
        raise ValueError(f"non-finite value in interaction entry {key}")
    ranks = {b: i for i, b in enumerate(sorted(set(tbytes)))}
    trank = np.array([ranks[b] for b in tbytes], dtype=np.int64)
    keys = np.concatenate([trank[:, None], modes, cells], axis=1)
    order = np.lexsort(keys.T[::-1])
    for start in range(0, n, 4096):
        out = []
        for i in order[start:start + 4096]:
            out.append(encode_string(str(types[i])))
            out.append(modes[i].tobytes())
            out.append(cells[i].tobytes())
            out.append(encode_complex(v[i]))
        yield b"".join(out)


def interaction_digest(records, *, version: int = FINGERPRINT_VERSION) -> str:
    """Digest of compiled interaction records, independent of their order.

    Pass ``CompiledInteraction.records`` (``so_records.compile_interaction``),
    the records with exact-zero coefficients excluded; its
    ``records_all`` (zeros kept) gives a different digest for the same
    interaction."""
    return digest("interaction", _record_chunks(records), version=version)


def _topology_chunks(rows):
    rows = np.asarray(rows, dtype="<i8").reshape(-1, 3)
    n = len(rows)
    yield encode_int(n)
    if n:
        rows = rows[np.lexsort(rows.T[::-1])]
        yield np.ascontiguousarray(rows).tobytes()


def topology_digest(declared_topology, *, version: int = FINGERPRINT_VERSION) -> str:
    """Digest of the declared displacement set, independent of row order."""
    return digest("topology", _topology_chunks(declared_topology), version=version)


def _i8_array_chunk(a) -> bytes:
    a = np.ascontiguousarray(a, dtype="<i8")
    return (encode_string("i8") + encode_int(a.ndim)
            + b"".join(encode_int(n) for n in a.shape) + a.tobytes())


def _slot_string(slot_convention) -> str:
    to_string = getattr(slot_convention, "to_string", None)
    return to_string() if callable(to_string) else str(slot_convention)


def resolved_topology_digest(delta_r, reverse, slot_convention, *,
                             version: int = FINGERPRINT_VERSION) -> str:
    """Digest of a resolved bond topology: the ``delta_r`` ``(B, 3)`` and
    ``reverse`` ``(B,)`` arrays in slot order (``"i8"`` arrays) and the slot
    convention string (``SlotConvention.to_string()``; a ``SlotConvention``
    is converted). The shell cutoff is not hashed: two cutoffs that keep the
    same channels are compatible."""
    dr = np.asarray(delta_r)
    rv = np.asarray(reverse)
    if dr.ndim != 2 or dr.shape[1] != 3:
        raise ValueError("delta_r must have shape (B, 3), got {}".format(dr.shape))
    if rv.shape != (dr.shape[0],):
        raise ValueError("reverse must have shape ({},), got {}".format(dr.shape[0], rv.shape))
    return digest("resolved_topology",
                  (_i8_array_chunk(dr), _i8_array_chunk(rv),
                   encode_string(_slot_string(slot_convention))), version=version)


class InteractionInputs(NamedTuple):
    """The compiled interaction of a run with its fingerprints."""
    records: np.ndarray
    declared_topology: np.ndarray
    fp_interaction: str
    fp_topology: str


def interaction_inputs(compiled, *, version: int = FINGERPRINT_VERSION
                       ) -> InteractionInputs:
    """:class:`InteractionInputs` of a ``so_records.CompiledInteraction``:
    ``fp_interaction`` covers ``compiled.records`` (exact zeros excluded),
    ``fp_topology`` covers ``compiled.declared_topology``."""
    return InteractionInputs(compiled.records, compiled.declared_topology,
                             interaction_digest(compiled.records, version=version),
                             topology_digest(compiled.declared_topology, version=version))


class BondInputs(NamedTuple):
    """A resolved bond topology as a bond-resolved producer records it:
    ``delta_r`` ``(B, 3)``, ``reverse`` ``(B,)``, the slot convention (string
    or ``SlotConvention``) and the shell cutoff ``max_shells``."""
    delta_r: np.ndarray
    reverse: np.ndarray
    slot_convention: object
    max_shells: int


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


def table_digest(component: str, table: Optional[dict], *,
                 version: int = FINGERPRINT_VERSION) -> str:
    """Digest of a ``wan90.read_w90``-style table ``{(irvec, orbvec): complex}``.

    ``None`` (absent) encodes only ``False``; a present table, even an empty
    one, encodes ``True``, the count and the sorted entries.
    """
    return digest(component, _table_chunks(component, table), version=version)


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


def band_digests(transfer: dict, extern: Optional[dict], meta: BandMeta, *,
                 version: int = FINGERPRINT_VERSION) -> dict:
    fp_transfer = table_digest("transfer", transfer, version=version)
    fp_extern = table_digest("extern", extern, version=version)
    fp_meta = digest("meta", _meta_chunks(meta), version=version)
    fp_band = digest("band", (encode_string(s) for s in (fp_transfer, fp_extern, fp_meta)),
                     version=version)
    return {"fp_transfer": fp_transfer, "fp_extern": fp_extern,
            "fp_meta": fp_meta, "fp_band": fp_band}


def band_inputs(transfer, extern, meta, *, version: int = FINGERPRINT_VERSION) -> BandInputs:
    return BandInputs(band_digests(transfer, extern, meta, version=version),
                      transfer, extern, meta)


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
            raise ValueError(
                f"layout {layout!r} needs shape (nd, nd, Nx, Ny, Nz, Nmat), got {shape}")
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


def green_digest(green, layout: str, *, version: int = FINGERPRINT_VERSION) -> str:
    """Digest of the canonical ``(Nmat, nvol, nd, nd)`` complex Green array.

    The array is hashed block by block (a CuPy array is moved to the host one
    block at a time); the caller's array is never modified.
    """
    return digest("green", _green_chunks(green, layout), version=version)


def green_digest_blocks(nmat: int, nvol: int, nd: int, block, *,
                        version: int = FINGERPRINT_VERSION) -> str:
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
    return digest("green", _block_chunks(int(nmat), int(nvol), int(nd), host_block),
                  version=version)


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

#: The members every spin-orbital archive carries, in a fixed order: the
#: band, Green and state fingerprints and their provenance tables, then the
#: interaction fingerprints, the capability, the resolved topology (or its
#: on-site sentinel) and their provenance tables. Consumers validate an
#: archive by THIS list, never by a name prefix such as ``state_``: FLEX
#: archives also carry an unrelated ``state_iteration`` member.
SO_PROVENANCE_MEMBERS = (
    "fingerprint_version", "producer_stage", "producer_label", "producer_version",
    "fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green",
    "prov_transfer_keys", "prov_transfer_vals", "prov_extern_present",
    "prov_extern_keys", "prov_extern_vals",
    "prov_coeff_extern", "prov_enable_spin_orbital", "prov_index_convention",
    "prov_norb", "prov_norb_phys", "prov_cell_shape", "prov_sub_shape",
    "state_beta", "state_nmat", "state_mu", "state_density", "state_calc_mu",
    "state_target", "state_mu_tol",
    "fp_interaction", "fp_topology", "fp_resolved_topology", "fp_capability",
    "prov_rec_types", "prov_rec_modes", "prov_rec_cells", "prov_rec_v", "prov_topology",
)
#: the members a bond-resolved ("so_bond") stage adds
SO_BOND_MEMBERS = ("delta_r", "reverse", "slot_convention", "max_shells")
#: optional in every archive: the digest of the Green function a producer
#: consumed, when it differs from the one it records in fp_green
OPTIONAL_MEMBERS = ("fp_green_source",)

#: The provenance schema: member -> (dtype, shape rule).
#: "str" is any unicode width; numeric dtypes match in either byte order. A
#: shape rule is a tuple of ints and dimension names; one name has one value
#: across the members of an archive.
PROVENANCE_SCHEMA = {
    "fingerprint_version": ("int64", ()),
    "producer_stage": ("str", ()), "producer_label": ("str", ()),
    "producer_version": ("str", ()),
    "fp_transfer": ("str", ()), "fp_extern": ("str", ()), "fp_meta": ("str", ()),
    "fp_band": ("str", ()), "fp_green": ("str", ()),
    "prov_transfer_keys": ("int64", ("n_transfer", 5)),
    "prov_transfer_vals": ("complex128", ("n_transfer",)),
    "prov_extern_present": ("bool", ()),
    "prov_extern_keys": ("int64", ("n_extern", 5)),
    "prov_extern_vals": ("complex128", ("n_extern",)),
    "prov_coeff_extern": ("float64", ()), "prov_enable_spin_orbital": ("bool", ()),
    "prov_index_convention": ("str", ()),
    "prov_norb": ("int64", ()), "prov_norb_phys": ("int64", ()),
    "prov_cell_shape": ("int64", (3,)), "prov_sub_shape": ("int64", (3,)),
    "state_beta": ("float64", ()), "state_nmat": ("int64", ()), "state_mu": ("float64", ()),
    "state_density": ("float64", ()), "state_calc_mu": ("bool", ()),
    "state_target": ("float64", ()), "state_mu_tol": ("float64", ()),
    "fp_interaction": ("str", ()), "fp_topology": ("str", ()),
    "fp_resolved_topology": ("str", ()), "fp_capability": ("str", ()),
    "prov_rec_types": ("str", ("n_rec",)),
    "prov_rec_modes": ("int64", ("n_rec", 4)),
    "prov_rec_cells": ("int64", ("n_rec", 4, 3)),
    "prov_rec_v": ("complex128", ("n_rec",)),
    "prov_topology": ("int64", ("n_topology", 3)),
    "fp_green_source": ("str", ()),
    "delta_r": ("int64", ("n_bond", 3)),
    "reverse": ("int64", ("n_bond",)),
    "slot_convention": ("str", ()),
    "max_shells": ("int64", ()),
}


def required_members(kind: str, stage: str) -> frozenset:
    """The provenance members a ``kind`` archive written by ``stage`` must
    carry (``fp_green_source`` is optional and not listed)."""
    if kind not in KIND_STAGES:
        raise ValueError("unknown artifact kind {!r}".format(kind))
    if stage not in KIND_STAGES[kind]:
        raise ValueError("artifact kind {!r} has no producer_stage {!r} (its stages: {})"
                         .format(kind, stage, ", ".join(KIND_STAGES[kind])))
    bond = SO_BOND_MEMBERS if STAGE_CAPABILITY[stage] == "so_bond" else ()
    return frozenset(SO_PROVENANCE_MEMBERS + bond)


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


def _is_hex64(s) -> bool:
    return (isinstance(s, str) and len(s) == 64
            and all(c in "0123456789abcdef" for c in s))


def so_provenance_members(*, producer_stage, producer_label, transfer, extern,
                          meta: BandMeta, fp_green: str, state: SOState,
                          interaction: InteractionInputs,
                          bond: Optional[BondInputs] = None,
                          fp_green_source: Optional[str] = None) -> dict:
    """The fingerprint and state members of a spin-orbital archive
    (:data:`SO_PROVENANCE_MEMBERS`, at :data:`FINGERPRINT_VERSION`).

    Build it once per run and merge the same dict into every archive of
    that run. Refuses (ValueError naming the member) any non-finite float
    or complex member.

    ``interaction`` is the run's :func:`interaction_inputs`; ``bond``
    (:class:`BondInputs`) is given exactly when the stage is bond-resolved;
    ``fp_green`` is a Green digest at :data:`FINGERPRINT_VERSION`.
    ``fp_green_source``, when given, is recorded and must differ from
    ``fp_green``.
    """
    if producer_stage not in PRODUCER_STAGES:
        raise ValueError("unknown producer_stage {!r}".format(producer_stage))
    if not _is_hex64(fp_green):
        raise ValueError("fp_green must be a 64-character lowercase hex digest, got {!r}"
                         .format(fp_green))
    is_bond = STAGE_CAPABILITY[producer_stage] == "so_bond"
    if interaction is None:
        raise ValueError("the spin-orbital provenance needs the interaction fingerprints "
                         "(interaction=)")
    if is_bond and bond is None:
        raise ValueError("producer_stage '{}' needs the resolved bond topology"
                         .format(producer_stage))
    if not is_bond and bond is not None:
        raise ValueError("producer_stage '{}' is not bond-resolved; it takes no "
                         "resolved bond topology".format(producer_stage))
    if fp_green_source is not None:
        if not _is_hex64(fp_green_source):
            raise ValueError("fp_green_source must be a 64-character lowercase hex "
                             "digest, got {!r}".format(fp_green_source))
        if fp_green_source == fp_green:
            raise ValueError("fp_green_source must differ from fp_green")
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
    m.update(_interaction_members(interaction))
    m["fp_capability"] = np.str_(STAGE_CAPABILITY[producer_stage])
    order = SO_PROVENANCE_MEMBERS
    if is_bond:
        m.update(_bond_members(bond))
        order = order + SO_BOND_MEMBERS
    else:
        m["fp_resolved_topology"] = np.str_(ONSITE_RESOLVED_SENTINEL)
    if fp_green_source is not None:
        m["fp_green_source"] = np.str_(fp_green_source)
        order = order + OPTIONAL_MEMBERS
    if set(m) != set(order):  # pragma: no cover - schema drift guard
        raise RuntimeError("spin-orbital provenance members drifted from the schema")
    return {k: m[k] for k in order}


def _interaction_members(interaction: InteractionInputs) -> dict:
    rec = interaction.records
    n = len(rec)
    for k in ("fp_interaction", "fp_topology"):
        if not _is_hex64(getattr(interaction, k)):
            raise ValueError("{} must be a 64-character lowercase hex digest, got {!r}"
                             .format(k, getattr(interaction, k)))
    v = np.ascontiguousarray(np.asarray(rec["v"], dtype=np.complex128).reshape(n))
    if not np.all(np.isfinite(v)):
        raise ValueError("non-finite spin-orbital provenance member prov_rec_v")
    return {
        "fp_interaction": np.str_(interaction.fp_interaction),
        "fp_topology": np.str_(interaction.fp_topology),
        "prov_rec_types": np.array([str(t) for t in np.asarray(rec["type"]).reshape(n)],
                                   dtype=np.str_).reshape(n),
        "prov_rec_modes": np.ascontiguousarray(np.asarray(rec["modes"], dtype=np.int64)
                                               .reshape(n, 4)),
        "prov_rec_cells": np.ascontiguousarray(np.asarray(rec["cells"], dtype=np.int64)
                                               .reshape(n, 4, 3)),
        "prov_rec_v": v,
        "prov_topology": np.ascontiguousarray(
            np.asarray(interaction.declared_topology, dtype=np.int64).reshape(-1, 3)),
    }


def _bond_members(bond: BondInputs) -> dict:
    slot = _slot_string(bond.slot_convention)
    dr = np.ascontiguousarray(np.asarray(bond.delta_r, dtype=np.int64))
    rv = np.ascontiguousarray(np.asarray(bond.reverse, dtype=np.int64))
    return {
        "fp_resolved_topology": np.str_(resolved_topology_digest(dr, rv, slot)),
        "delta_r": dr,
        "reverse": rv,
        "slot_convention": np.str_(slot),
        "max_shells": np.int64(bond.max_shells),
    }


# ---------------------------------------------------------- consumer checks
def _schema_dtype_ok(dtype, name):
    """``dtype`` is the schema dtype ``name`` (any byte order; "str": any
    unicode width)."""
    if name == "str":
        return dtype.kind == "U"
    return dtype.newbyteorder("<") == np.dtype(name).newbyteorder("<")


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
    state and never copy these fields."""
    producer_state: SOState
    producer_fp_green: str
    producer_label: str
    version: int


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


def _malformed_raiser(path, label):
    def malformed(check, stored, expected):
        raise SOArchiveError(
            "'{}' is malformed or truncated: {} ({} vs {}). Regenerate it with "
            "{}.".format(path, check, stored, expected, label))
    return malformed


def _point(nmat_rule):
    # load mode follows the existing stored-axis rule for Nmat (see the state check)
    return ("the identical state point (Nmat excepted in chi0q_mode = \"load\")"
            if nmat_rule == "load" else "the identical state point")


def _file_band(data):
    """The band tables and metadata an archive's ``prov_*`` members record
    (for the diffs) and whether it records an Extern table."""
    tk, tv = np.asarray(data["prov_transfer_keys"]), np.asarray(data["prov_transfer_vals"])
    ek, ev = np.asarray(data["prov_extern_keys"]), np.asarray(data["prov_extern_vals"])
    ext_present = bool(_scalar(data, "prov_extern_present"))
    file_meta = BandMeta(
        float(_scalar(data, "prov_coeff_extern")),
        bool(_scalar(data, "prov_enable_spin_orbital")),
        str(_scalar(data, "prov_index_convention")),
        int(_scalar(data, "prov_norb")), int(_scalar(data, "prov_norb_phys")),
        tuple(int(x) for x in np.asarray(data["prov_cell_shape"])),
        tuple(int(x) for x in np.asarray(data["prov_sub_shape"])))
    return (_table_from_arrays(tk, tv),
            _table_from_arrays(ek, ev) if ext_present else None, file_meta, ext_present)


def _file_state(data) -> SOState:
    return SOState(float(_scalar(data, "state_beta")), int(_scalar(data, "state_nmat")),
                   float(_scalar(data, "state_mu")), float(_scalar(data, "state_density")),
                   bool(_scalar(data, "state_calc_mu")), float(_scalar(data, "state_target")),
                   float(_scalar(data, "state_mu_tol")))


def _check_mu_tol(state, malformed):
    """The ``calc_mu`` tagged union of a recorded state (structural)."""
    t = state.target
    scale = max(1.0, abs(t))
    if state.calc_mu:
        resid = abs(state.density - t)
        if abs(state.mu_tol - resid) > MU_TOL_CONSISTENCY * scale:
            malformed("state_mu_tol is not |state_density - state_target|",
                      repr(state.mu_tol), repr(resid))
    elif state.mu_tol != 0.0:
        malformed("state_mu_tol of a fixed-mu state", repr(state.mu_tol), "0.0")


def _check_residual(path, label, state):
    """A searched-mu producer must have reached its filling."""
    scale = max(1.0, abs(state.target))
    if state.calc_mu and state.mu_tol > RESIDUAL_CAP * scale:
        raise SOArchiveError(
            "'{}' cannot be reused: the producer did not reach its filling "
            "(residual {:.3e} electrons per cell > 1e-8·max(1,|t|) = {:.3e}). "
            "Regenerate it with {} so that its density closes to that "
            "tolerance.".format(path, state.mu_tol, RESIDUAL_CAP * scale, label))


def _check_band(path, label, file_fps, run_fps, file_band, run_band, point):
    """Refuse when a band digest differs, naming the differing fields and
    rows from the archive's ``prov_*`` tables (diff only)."""
    file_transfer, file_extern, file_meta, ext_present = file_band
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
    if not diffs and file_fps["fp_band"] != run_fps["fp_band"]:
        # fail closed: the combined digest differs although its parts agree
        diffs.append("fp_band: file {}, run {} (fp_transfer, fp_extern and fp_meta agree; "
                     "the band digest is inconsistent)".format(file_fps["fp_band"][:12],
                                                               run_fps["fp_band"][:12]))
    if diffs:
        raise SOArchiveError(
            "'{}' (written by {}) does not match this run: band differs: {}. "
            "Archives are reused only at {} (not a warm start across a sweep); "
            "re-run {} with this input (the same Transfer/Extern/coeff_extern)."
            .format(path, label, "; ".join(diffs), point, label))


def _check_state(path, label, state, run_state, nmat_rule, point):
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
            "Archives are reused only at {}; {}.".format(
                path, label, "; ".join(diffs), point, fix))
    if nmat_rule == "load" and state.nmat != int(run_state.nmat):
        check_logger.warning(
            "'%s' (written by %s) was summed on another Matsubara grid: "
            "state_nmat %d (file) vs Nmat %d (run); the frequency truncation of "
            "the loaded chi0 differs from this run's. Set [mode.param] Nmat = %d "
            "to reproduce the producer.", path, label, state.nmat,
            int(run_state.nmat), state.nmat)


def _check_green(path, label, data, green, green_layout, green_fp):
    file_green = str(_scalar(data, "fp_green"))
    if green_fp is None and green is not None:
        green_fp = green_digest(green, green_layout)
    if green_fp is not None:
        if green_fp != file_green:
            at = " (both at fingerprint_version {})".format(FINGERPRINT_VERSION)
            raise SOArchiveError(
                "'{}' (written by {}): fp_green mismatch: the Green function this "
                "run loaded has digest {}..., the archive records {}...{}; the Green "
                "function and the susceptibility do not come from the same run. "
                "Regenerate both with {}.".format(path, label, green_fp[:12],
                                                 file_green[:12], at, label))
    else:
        check_logger.info("fp_green: producer assertion (%s)", path)
    return file_green


# --------------------------------------------------------- structure (step 1)
def _zip_of(data, archive):
    """``(ZipFile or None, close)`` holding the payload members' headers."""
    if archive is None:
        z = getattr(data, "zip", None)
        return (z, False) if isinstance(z, zipfile.ZipFile) else (None, False)
    if isinstance(archive, zipfile.ZipFile):
        return archive, False
    z = getattr(archive, "zip", None)
    if isinstance(z, zipfile.ZipFile):
        return z, False
    if isinstance(archive, (str, os.PathLike)):
        return zipfile.ZipFile(archive), True
    raise ValueError("archive must be a zip path, a ZipFile or an NpzFile, got {!r}"
                     .format(type(archive).__name__))


def _read_header(zf, member, p, malformed):
    """``(npy header, data offset)`` of the zip member ``member`` (payload
    ``p`` in the messages); the array is never read."""
    from . import npy_header
    try:
        with zf.open(member) as fh:
            h = npy_header.read_npy_header(fh)
            return h, fh.tell()
    except (ValueError, OSError, zipfile.BadZipFile,
            npy_header.UnsupportedNpyHeaderVersion) as exc:
        malformed("payload member {} has an unreadable npy header".format(p),
                  repr(str(exc)), "a valid npy header")


def _check_payload_headers(zf, kind, stage, malformed):
    """The payload header check: dtype complex128, C order and the recorded rank of each
    payload member, from its npy header, and its stored size (header plus the declared
    data), from the zip directory; the array is never read."""
    names = set(zf.namelist())
    for p in PAYLOAD_MEMBERS[kind]:
        if p in _STRING_PAYLOADS:
            continue
        member = p + ".npy" if p + ".npy" in names else p
        if member not in names:
            malformed("payload member {} of a {} archive is missing".format(p, kind),
                      "absent", "present")
        h, data_offset = _read_header(zf, member, p, malformed)
        if not _schema_dtype_ok(np.dtype(h.dtype), "complex128"):
            malformed("payload member {} has the wrong dtype".format(p),
                      str(h.dtype), "complex128")
        if h.fortran_order:
            malformed("payload member {} is in Fortran order".format(p), "Fortran", "C")
        ranks = PAYLOAD_RANKS[(kind, stage)]
        if len(h.shape) not in ranks:
            malformed("payload member {} has the wrong rank".format(p), str(len(h.shape)),
                      " or ".join(str(r) for r in ranks))
        # the stored size is the header plus the declared data, read from the
        # zip directory (a truncated member would otherwise pass this step)
        size = zf.getinfo(member).file_size
        want = data_offset + int(np.prod(h.shape, dtype=np.int64)) * np.dtype(h.dtype).itemsize
        if size != want:
            malformed("payload member {} has the wrong size".format(p), str(size), str(want))


# ------------------------------------------------- the slot layout (step 1)
def _check_rules(rules, kind, stage):
    """A :class:`SlotRules` table is complete for the archive's ``kind`` and
    ``stage``: every slot kind has its ``payloads`` entry, and a slot kind
    has a ``ranks`` entry for ``stage`` naming every required and optional
    payload. An incomplete table is a programming error: a ``ValueError``
    before any archive member is read."""
    def missing(k, what):
        raise ValueError("SlotRules: slot kind {} has no {} entry".format(k, what))

    for k in sorted(rules.slot_kinds):
        if k not in rules.payloads:
            missing(k, "payloads")
    if kind in rules.slot_kinds:
        ranks = rules.ranks.get((kind, stage))
        if ranks is None:
            missing(kind, "ranks for stage {}".format(stage))
        for p in tuple(rules.payloads[kind]) + tuple(rules.optional.get(kind, ())):
            if p not in ranks:
                missing(kind, "rank for payload {} on stage {}".format(p, stage))


def _is_reserved(name):
    return isinstance(name, str) and name.startswith(RESERVED_PREFIX)


def check_layout_members(data, names, *, kind, stage, B, nd, malformed, path, label,
                         rules=PRODUCTION_SLOT_RULES) -> Optional[str]:
    """The layout members of an archive of ``kind`` written by ``stage``;
    returns its ``block_layout`` (``"dense"`` or ``"sz_packed"``), or
    ``None`` for a stage that is not bond-resolved.

    ``names`` are the archive's member names, ``B`` its channel count
    (``len(delta_r)``) and ``nd = 2 norb_phys``; ``malformed`` raises the
    single-line refusal. The rules, each naming the first offending member
    in sorted order:

    * a ``so_onsite`` stage carries neither ``block_layout`` nor any name
      with the reserved prefix :data:`RESERVED_PREFIX`;
    * a ``so_bond`` stage carries ``block_layout``, a 0-d string
      ``"dense"`` or ``"sz_packed"``;
    * a gap (always stored full) carries no reserved name;
    * a slot kind of ``rules``: a dense archive carries no reserved name; a
      packed one carries only the three sector-map members under the
      prefix, and they equal the map of ``(B, nd)``
      (``so_bond.check_sector_map_members``).

    Any other member is accepted (forward compatible)."""
    _check_rules(rules, kind, stage)
    reserved = sorted(n for n in names if _is_reserved(n))
    if STAGE_CAPABILITY[stage] != "so_bond":
        bad = sorted((["block_layout"] if "block_layout" in names else []) + reserved)
        if bad:
            malformed("member {} is not written by a so_onsite stage".format(bad[0]),
                      "present", "absent")
        return None
    if "block_layout" not in names:
        malformed("member block_layout is missing", "absent", "present")
    v = np.asarray(data["block_layout"])
    if v.ndim != 0 or v.dtype.kind != "U" or str(v[()]) not in _SLOT_LAYOUTS:
        malformed("block_layout is not a known layout", repr(v.tolist()),
                  "'dense' or 'sz_packed'")
    layout = str(v[()])
    if kind == "gap":
        if reserved:
            malformed("member {} is not written by a gap archive".format(reserved[0]),
                      "present", "absent")
        return layout
    if kind in rules.slot_kinds:
        if layout == "dense":
            if reserved:
                malformed("member {} is not written by a dense slot archive"
                          .format(reserved[0]), "present", "absent")
            return layout
        stray = [n for n in reserved if n not in LAYOUT_MEMBERS[1:]]
        if stray:
            malformed("member {} is not a sector-map member (the prefix {} is reserved)"
                      .format(stray[0], RESERVED_PREFIX), "present", "absent")
        from . import so_bond
        try:
            so_bond.check_sector_map_members(data, B, nd, file=path, producer_label=label)
        except SOArchiveError:
            raise
        except ValueError as exc:
            raise SOArchiveError(str(exc)) from None
    return layout


def check_slot_payloads(data, names, zf, *, kind, stage, layout, B, nd, sector_dims,
                        malformed, rules) -> None:
    """The payload members of a slot archive of ``kind`` (a slot kind of
    ``rules``) in ``layout`` (from :func:`check_layout_members`).

    Names: every required payload ``p``, and every optional one present
    under any of its names, is stored as exactly ``<p>`` (dense) or exactly
    ``<p>@0``, ``<p>@+1``, ``<p>@-1`` (packed). Then each block's metadata:
    complex128 (either byte order), C order, a rank of
    ``rules.ranks[(kind, stage)]``, the stored size (zip only) and the
    trailing pair ``(ND, ND)`` (``ND = B nd**2``; dense) or ``(n_s, n_s)``
    (``n_s = sector_dims[i]`` in sector order ``0, +1, -1``; packed). With
    a zip ``zf`` the metadata comes from the npy headers and no payload is
    read; for a plain mapping (``zf is None``) from its arrays, so the shape
    check is never skipped."""
    _check_rules(rules, kind, stage)
    if layout not in _SLOT_LAYOUTS:
        raise ValueError("check_slot_payloads: layout must be 'dense' or 'sz_packed', got "
                         "{!r}".format(layout))
    packed = layout == "sz_packed"
    if packed:
        from . import so_bond
        order = so_bond.SECTOR_ORDER
        if sector_dims is None or len(sector_dims) != len(order):
            raise ValueError("check_slot_payloads: a packed archive needs the {} sector "
                             "dimensions, got {!r}".format(len(order), sector_dims))
        pairs = [(int(n), int(n)) for n in sector_dims]
    ND = int(B) * int(nd) ** 2
    ranks = rules.ranks[(kind, stage)]

    def stored_as(p):
        return sorted(n for n in names
                      if isinstance(n, str) and (n == p or n.startswith(p + "@")))

    present_optional = tuple(p for p in rules.optional.get(kind, ()) if stored_as(p))
    blocks = []
    for p in tuple(rules.payloads[kind]) + present_optional:
        if packed:
            expected = [so_bond._packed_key(p, s) for s in order]
            blocks.extend((name, p, pair) for name, pair in zip(expected, pairs))
        else:
            expected = [p]
            blocks.append((p, p, (ND, ND)))
        found = stored_as(p)
        if found != sorted(expected):
            malformed("payload member {} does not follow block_layout {}".format(p, layout),
                      str(found), str(expected))

    zip_names = set(zf.namelist()) if zf is not None else None
    for name, p, pair in blocks:
        if zf is not None:
            member = name + ".npy" if name + ".npy" in zip_names else name
            if member not in zip_names:
                malformed("payload member {} of a {} archive is missing".format(name, kind),
                          "absent", "present")
            h, data_offset = _read_header(zf, member, name, malformed)
            dtype, shape, order_c = np.dtype(h.dtype), tuple(h.shape), not h.fortran_order
            stored = "Fortran"
        else:
            a = np.asarray(data[name])
            dtype, shape, order_c = a.dtype, tuple(a.shape), bool(a.flags.c_contiguous)
            stored = "Fortran" if a.flags.f_contiguous else "not contiguous"
        if not _schema_dtype_ok(dtype, "complex128"):
            malformed("payload member {} has the wrong dtype".format(name), str(dtype),
                      "complex128")
        if not order_c:
            malformed("payload member {} is in Fortran order".format(name), stored, "C")
        if len(shape) not in ranks[p]:
            malformed("payload member {} has the wrong rank".format(name), str(len(shape)),
                      " or ".join(str(r) for r in ranks[p]))
        if zf is not None:
            size = zf.getinfo(member).file_size
            want = data_offset + int(np.prod(shape, dtype=np.int64)) * dtype.itemsize
            if size != want:
                malformed("payload member {} has the wrong size".format(name), str(size),
                          str(want))
        if shape[-2:] != pair:
            malformed("payload member {} has the wrong slot shape".format(name),
                      str(shape[-2:]), str(pair))


def _first_non_finite(a):
    if a.ndim == 0:
        return repr(a.tolist())
    i = int(np.flatnonzero(~np.isfinite(a.reshape(-1)))[0])
    return "entry {}: {!r}".format(i, a.reshape(-1)[i].item())


#: the refusal of a version-1 archive
_VERSION_1_TEXT = (
    "'{path}' carries fingerprint_version 1 (written before interaction fingerprints "
    "existed); this version requires fingerprint_version 2. Regenerate it by re-running "
    "{producer} with this version (the run and the [file.output] key that write each "
    "archive are listed under 'Archive checks' in the hwave_sc documentation).")
#: the refusal of an unversioned archive
_UNVERSIONED_TEXT = (
    "'{path}' has no fingerprint_version (written before this version, or by a run whose "
    "band the fingerprints do not cover: trans_mod, green_init or a supplied chi0q; such "
    "outputs cannot be read by a spin-orbital consumer); this version requires "
    "fingerprint_version 2. Regenerate it by re-running {producer} with this version, "
    "from the Transfer/Extern input (the run and the [file.output] key that write each "
    "archive are listed under 'Archive checks' in the hwave_sc documentation).")


def _validate(data, path, kind, archive, label, expected_producer):
    """:func:`validate_structure` with the producer named ``label`` (the
    readable ``producer_label``, else ``expected_producer``)."""
    if kind not in KIND_STAGES:
        raise ValueError("unknown artifact kind {!r}".format(kind))
    names = _members(data)
    malformed = _malformed_raiser(path, label)

    # the version first: an archive of another version is never reported
    # as malformed
    if "fingerprint_version" not in names:
        raise SOArchiveError(_UNVERSIONED_TEXT.format(
            path=path, producer=expected_producer or label))
    ver = np.asarray(data["fingerprint_version"])
    if ver.ndim != 0 or ver.dtype.kind not in "iu":
        malformed("fingerprint_version is not an integer scalar",
                  "{}/{}d".format(ver.dtype, ver.ndim), "int64/0d")
    version = int(ver[()])
    if version == 1:
        raise SOArchiveError(_VERSION_1_TEXT.format(path=path, producer=label))
    if version != FINGERPRINT_VERSION:
        malformed("unsupported fingerprint_version", version, FINGERPRINT_VERSION)

    if "producer_stage" not in names:
        malformed("member producer_stage is missing", "absent", "present")
    st = np.asarray(data["producer_stage"])
    if st.ndim != 0 or st.dtype.kind != "U":
        malformed("member producer_stage has the wrong dtype or rank",
                  "{}/{}d".format(st.dtype, st.ndim), "str/0d")
    stage = str(st[()])
    if stage in DEFERRED_STAGES:
        raise SOArchiveError(
            "'{}' was written by producer_stage '{}', which this version does not read "
            "yet ({} arrives in a later version).".format(path, stage, DEFERRED_STAGES[stage]))
    if stage not in PRODUCER_STAGES:
        malformed("producer_stage is not a known stage", repr(stage),
                  "one of {}".format(", ".join(PRODUCER_STAGES)))
    if stage not in KIND_STAGES[kind]:
        malformed("producer_stage does not write a {} archive".format(kind), repr(stage),
                  "one of {}".format(", ".join(KIND_STAGES[kind])))
    is_bond = STAGE_CAPABILITY[stage] == "so_bond"

    # (a) the provenance schema
    listed = SO_PROVENANCE_MEMBERS + (SO_BOND_MEMBERS if is_bond else ())
    listed = listed + tuple(k for k in OPTIONAL_MEMBERS if k in names)
    dims = {}
    m = {}
    for k in listed:
        if k not in names:
            malformed("member {} is missing".format(k), "absent", "present")
        dname, rule = PROVENANCE_SCHEMA[k]
        a = np.asarray(data[k])
        if not _schema_dtype_ok(a.dtype, dname) or a.ndim != len(rule):
            malformed("member {} has the wrong dtype or rank".format(k),
                      "{}/{}d".format(a.dtype, a.ndim), "{}/{}d".format(dname, len(rule)))
        want = tuple(dims.setdefault(r, n) if isinstance(r, str) else r
                     for n, r in zip(a.shape, rule))
        if tuple(a.shape) != want:
            malformed("member {} has the wrong shape".format(k), str(tuple(a.shape)),
                      str(want))
        if a.dtype.kind in "fc" and not np.all(np.isfinite(a)):
            malformed("member {} is not finite".format(k), _first_non_finite(a), "finite")
        m[k] = a

    def text(k):
        return str(m[k][()])

    for k in ("fp_transfer", "fp_extern", "fp_meta", "fp_band", "fp_green",
              "fp_interaction", "fp_topology"):
        if not _is_hex64(text(k)):
            malformed("{} is not a 64-character lowercase hex digest".format(k),
                      repr(text(k)), "64 hex characters")
    cap = text("fp_capability")
    if cap not in CAPABILITIES:
        malformed("fp_capability is not a known capability", repr(cap),
                  "one of {}".format(", ".join(CAPABILITIES)))
    if cap != STAGE_CAPABILITY[stage]:
        malformed("fp_capability does not match producer_stage '{}'".format(stage),
                  repr(cap), repr(STAGE_CAPABILITY[stage]))
    rt = text("fp_resolved_topology")
    if is_bond and not _is_hex64(rt):
        malformed("fp_resolved_topology is not a 64-character lowercase hex digest",
                  repr(rt), "64 hex characters")
    if not is_bond:
        if rt != ONSITE_RESOLVED_SENTINEL:
            malformed("fp_resolved_topology of a so_onsite stage", repr(rt),
                      repr(ONSITE_RESOLVED_SENTINEL))
        for k in SO_BOND_MEMBERS:
            if k in names:
                malformed("member {} is not written by a so_onsite stage".format(k),
                          "present", "absent")
    if "fp_green_source" in m:
        src = text("fp_green_source")
        if not _is_hex64(src):
            malformed("fp_green_source is not a 64-character lowercase hex digest",
                      repr(src), "64 hex characters")
        if src == text("fp_green"):
            malformed("fp_green_source equals fp_green", repr(src), "a different digest")
    if is_bond:
        from . import so_bond
        slot = text("slot_convention")
        try:
            parsed = so_bond.SlotConvention.from_string(slot)
        except ValueError:
            malformed("slot_convention does not parse", repr(slot),
                      "the slot-convention grammar")
        if parsed != so_bond.SLOT_CONVENTION:
            malformed("slot_convention is not this version's", repr(slot),
                      repr(so_bond.SLOT_CONVENTION.to_string()))

    # the recorded band and state
    if not bool(m["prov_extern_present"][()]) and m["prov_extern_keys"].shape[0] != 0:
        malformed("prov_extern_keys of an absent Extern",
                  "{} rows".format(m["prov_extern_keys"].shape[0]), "0 rows")
    beta_f = float(m["state_beta"][()])
    if not beta_f > 0:
        malformed("state_beta must be > 0", repr(beta_f), "> 0")
    _check_mu_tol(_file_state(data), malformed)

    # duplicated operational fields agree (a folded cell_shape is not compared)
    if "beta" in names:
        b = np.asarray(data["beta"])
        if b.ndim != 0 or b.dtype.kind not in "fi" or float(b[()]) != beta_f:
            malformed("beta differs from state_beta", repr(b.tolist()), repr(beta_f))
    nmat_f = int(m["state_nmat"][()])
    if "nmat" in names:
        n = np.asarray(data["nmat"])
        if n.ndim != 0 or n.dtype.kind not in "iu" or int(n[()]) != nmat_f:
            malformed("nmat differs from state_nmat", repr(n.tolist()), nmat_f)
    norb_phys = int(m["prov_norb_phys"][()])
    if "norb_phys" in names:
        n = np.asarray(data["norb_phys"])
        if n.ndim != 0 or n.dtype.kind not in "iu" or int(n[()]) != norb_phys:
            malformed("norb_phys differs from prov_norb_phys", repr(n.tolist()), norb_phys)

    # the layout members: block_layout of a bond-resolved stage, read first;
    # the reserved sz_sector_ names
    nd = 2 * norb_phys
    n_bonds = int(m["delta_r"].shape[0]) if is_bond else None
    layout = check_layout_members(data, names, kind=kind, stage=stage, B=n_bonds, nd=nd,
                                  malformed=malformed, path=path, label=label)

    # (b) the payload names and markers
    for p in PAYLOAD_MEMBERS.get(kind, ()):
        if p not in names:
            malformed("payload member {} of a {} archive is missing".format(p, kind),
                      "absent", "present")
    if "index_convention" not in names:
        malformed("marker index_convention is missing", "absent", "present")
    ic = np.asarray(data["index_convention"])
    if ic.ndim != 0 or ic.dtype.kind != "U":
        malformed("member index_convention has the wrong dtype or rank",
                  "{}/{}d".format(ic.dtype, ic.ndim), "str/0d")
    if str(ic[()]) != text("prov_index_convention"):
        malformed("index_convention differs from prov_index_convention", repr(str(ic[()])),
                  repr(text("prov_index_convention")))
    if "spin_orbital" not in names:
        malformed("marker spin_orbital is missing", "absent", "present")
    so = np.asarray(data["spin_orbital"])
    if so.ndim != 0 or so.dtype.kind != "b" or not bool(so[()]):
        malformed("marker spin_orbital is not true", repr(so.tolist()), "True")
    if kind == "green":
        gl = np.asarray(data["green_layout"])
        if gl.ndim != 0 or gl.dtype.kind != "U" or str(gl[()]) != GREEN_LAYOUT:
            malformed("green_layout is not the canonical layout", repr(gl.tolist()),
                      repr(GREEN_LAYOUT))

    # (c) the payload headers
    zf, close = _zip_of(data, archive)
    try:
        if kind in PRODUCTION_SLOT_RULES.slot_kinds:
            # a slot archive is checked from its arrays when there is no zip
            sector_dims = None
            if layout == "sz_packed":
                from . import so_bond
                sector_dims = so_bond.sz_sector_map(types.SimpleNamespace(
                    reverse=m["reverse"]), nd).sizes
            check_slot_payloads(data, names, zf, kind=kind, stage=stage, layout=layout,
                                B=n_bonds, nd=nd, sector_dims=sector_dims,
                                malformed=malformed, rules=PRODUCTION_SLOT_RULES)
        elif zf is not None:
            _check_payload_headers(zf, kind, stage, malformed)
    finally:
        if close:
            zf.close()
    return kind, stage, version


def validate_structure(data, path: str, *, kind: str, archive=None,
                       expected_producer: Optional[str] = None):
    """Step 1 of the archive check: the structure of a spin-orbital archive,
    before any payload is read. Returns ``(kind, stage, version)``.

    ``data`` is the open npz (or a mapping of its members), ``path`` its file
    name for the messages, ``kind`` one of :data:`ARTIFACT_KINDS`.
    ``expected_producer`` names the run that writes such an archive (for an
    unversioned archive, or one whose ``producer_label`` is unreadable). The
    checks, each a single-line :class:`SOArchiveError`:

    * ``fingerprint_version`` is :data:`FINGERPRINT_VERSION`; a version-1 or
      unversioned archive gets its own message naming the producer to
      re-run;
    * ``producer_stage`` is a stage the kind is written by
      (:data:`KIND_STAGES`; a stage of :data:`DEFERRED_STAGES` gets its own
      message);
    * (a) every member of the provenance schema is present with its dtype
      and shape (:data:`PROVENANCE_SCHEMA`), the digests are 64 lowercase
      hex, the capability and resolved topology agree with the stage, the
      recorded state is consistent and the duplicated fields (``beta``,
      ``nmat``, ``norb_phys``, ``index_convention``) agree with it;
    * the layout members (:func:`check_layout_members`): a bond-resolved
      stage records ``block_layout`` (``"dense"`` or ``"sz_packed"``), a
      gap carries no member with the reserved prefix ``sz_sector_``, and an
      on-site stage carries none of these members;
    * (b) the kind's payload names (:data:`PAYLOAD_MEMBERS`) and the
      ``index_convention`` and ``spin_orbital`` markers are present,
      ``spin_orbital`` is true and a green archive's ``green_layout`` is the
      canonical layout;
    * (c) when ``archive`` is a zip path, a ``ZipFile`` or an ``NpzFile``
      (or ``data`` is an ``NpzFile``), each complex payload member's npy
      header says complex128, C order and the recorded rank
      (:data:`PAYLOAD_RANKS`); the payload itself is never read. For a
      plain mapping the headers of the on-site kinds are not checked; a slot
      kind (:func:`check_slot_payloads`) is checked from the mapping's
      arrays, so its shapes and dtypes are never skipped.
    """
    label = _readable_label(data, _members(data))
    if label == UNKNOWN_PRODUCER and expected_producer:
        label = expected_producer
    return _validate(data, path, kind, archive, label, expected_producer)


def _record_table(types, modes, cells, vals):
    return {(str(t), tuple(int(x) for x in mo),
             tuple(tuple(int(x) for x in row) for row in np.asarray(c).reshape(4, 3))):
            complex(v) for t, mo, c, v in zip(types, modes, cells, vals)}


def _interaction_diff(data, run_records):
    """``file <n> records, run <m> records; first differing record <type>
    modes (...) cells (...): file <v>, run <w>`` from the archive's
    ``prov_rec_*`` members (diff only) and the run's records."""
    ft = np.asarray(data["prov_rec_types"])
    file_tbl = _record_table(ft, np.asarray(data["prov_rec_modes"]),
                             np.asarray(data["prov_rec_cells"]), np.asarray(data["prov_rec_v"]))
    n = len(run_records)
    run_tbl = _record_table(np.asarray(run_records["type"]).reshape(n),
                            np.asarray(run_records["modes"]).reshape(n, 4),
                            np.asarray(run_records["cells"]).reshape(n, 4, 3),
                            np.asarray(run_records["v"]).reshape(n))
    head = "file {} records, run {} records".format(len(ft), n)
    for k in sorted(set(file_tbl) | set(run_tbl)):
        a, b = file_tbl.get(k), run_tbl.get(k)
        if a is None or b is None or a != b:
            return "{}; first differing record {} modes {} cells {}: file {}, run {}".format(
                head, k[0], k[1], k[2], _fmt_value(a), _fmt_value(b))
    return head + " (the records agree; the digests differ)"


def check_so_archive(data, path: str, *, kind: str, run_band: BandInputs,
                     run_state: SOState, nmat_rule: str, run_interaction: InteractionInputs,
                     need: str, expected_producer: Optional[str] = None,
                     run_bond: Optional[BondInputs] = None, green=None,
                     green_layout: Optional[str] = None,
                     green_fp: Optional[str] = None) -> SOArchiveInfo:
    """Validate the provenance of a spin-orbital archive against this run.

    ``data`` is the open npz (or a mapping of its members), ``path`` its
    file name for the messages, ``kind`` its archive kind (one of
    :data:`ARTIFACT_KINDS`). ``run_band`` is this run's band
    (:func:`band_inputs`), ``run_state`` its state: ``beta``, ``nmat``,
    ``calc_mu`` and ``target`` (and ``mu`` for a fixed-mu run) are compared;
    ``density`` and ``mu_tol`` are not read. ``nmat_rule`` is ``"exact"``
    (a borrowed G) or ``"load"`` (load mode: the existing stored-axis rule
    applies instead; a different Nmat is not refused here, only reported
    with a WARNING). ``run_interaction`` is this run's
    :func:`interaction_inputs` and ``need`` the capability it consumes
    (``"so_onsite"`` or ``"so_bond"``). ``expected_producer`` names the run
    that writes such an archive (for an unversioned archive, or one whose
    label is unreadable). With ``green`` (in ``green_layout``) or
    ``green_fp`` (its digest, already computed by the caller) the archive's
    ``fp_green`` is verified; otherwise it is a producer assertion, logged
    at INFO.

    The checks: (1) :func:`validate_structure`; (2) band (this run's digests
    recomputed from its own tables; no digest is recomputed from the
    archive's ``prov_*`` members, which only feed the diffs), interaction,
    capability, declared and resolved topology (``need = "so_bond"`` only,
    against ``run_bond``), state, green.

    Every refusal is a single-line ValueError (:class:`SOArchiveError`)."""
    if nmat_rule not in NMAT_RULES:
        raise ValueError("nmat_rule must be one of {}, got {!r}".format(NMAT_RULES, nmat_rule))
    if run_interaction is None:
        raise ValueError("check_so_archive needs run_interaction")
    if need not in CAPABILITIES:
        raise ValueError("need must be one of {}, got {!r}".format(CAPABILITIES, need))
    label = _readable_label(data, _members(data))
    if label == UNKNOWN_PRODUCER and expected_producer:
        label = expected_producer
    kind, stage, version = _validate(data, path, kind, None, label, expected_producer)
    state = _file_state(data)
    _check_residual(path, label, state)
    point = _point(nmat_rule)

    # band: this run's digests, computed from its own tables
    keys = ("fp_transfer", "fp_extern", "fp_meta", "fp_band")
    file_fps = {k: str(_scalar(data, k)) for k in keys}
    run_fps = run_band.digests
    if any(file_fps[k] != run_fps[k] for k in keys):
        _check_band(path, label, file_fps, run_fps, _file_band(data), run_band, point)

    # interaction
    if str(_scalar(data, "fp_interaction")) != run_interaction.fp_interaction:
        raise SOArchiveError(
            "'{}' (written by {}) does not match this run: interaction differs: {}. "
            "Archives are reused only at {} (not a warm start across a sweep); re-run {} "
            "with this input (the same interaction definitions).".format(
                path, label, _interaction_diff(data, run_interaction.records), point, label))

    # capability
    cap = str(_scalar(data, "fp_capability"))
    if cap != need:
        raise SOArchiveError("{} was produced by {} ({}), but this run needs a {} artifact"
                             .format(path, label, cap, need))

    # topology: bond consumers only
    if need == "so_bond":
        if run_bond is None:
            raise ValueError("a so_bond consumer needs run_bond")
        if str(_scalar(data, "fp_topology")) != run_interaction.fp_topology:
            raise SOArchiveError(
                "'{}' (written by {}) does not match this run: declared topology differs: "
                "file {} displacements, run {} displacements. Archives are reused only at {} "
                "(not a warm start across a sweep); re-run {} with this input (the same "
                "interaction definitions).".format(
                    path, label, len(np.asarray(data["prov_topology"])),
                    len(np.asarray(run_interaction.declared_topology).reshape(-1, 3)),
                    point, label))
        run_rt = resolved_topology_digest(run_bond.delta_r, run_bond.reverse,
                                          run_bond.slot_convention)
        if str(_scalar(data, "fp_resolved_topology")) != run_rt:
            raise SOArchiveError(
                "'{}' (written by {}) does not match this run: resolved bond topology "
                "differs: file delta_r {} reverse {} (max_shells {}), run delta_r {} reverse "
                "{} (max_shells {}). Archives are reused only at {} (not a warm start across "
                "a sweep); re-run {} with this input (the same bond settings).".format(
                    path, label, np.asarray(data["delta_r"]).tolist(),
                    np.asarray(data["reverse"]).tolist(), int(_scalar(data, "max_shells")),
                    np.asarray(run_bond.delta_r).tolist(),
                    np.asarray(run_bond.reverse).tolist(), int(run_bond.max_shells),
                    point, label))

    _check_state(path, label, state, run_state, nmat_rule, point)
    file_green = _check_green(path, label, data, green, green_layout, green_fp)
    return SOArchiveInfo(state, file_green, label, version)
