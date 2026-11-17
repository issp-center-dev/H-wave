"""Synthetic inputs of the schema-4 spin-orbital bond archive (test only).

:func:`archive_inputs` builds the keyword arguments of
``so_bond_archive.write_so_bond_archive`` for a small synthetic case: ``B``
bond channels along x (``0, +1, -1, ...``), one physical orbital by
default, a ``cell`` mesh and ``nmat`` bosonic frequencies. The slot arrays
are random and block diagonal in the S_z sectors, so the same case can be
written dense or S_z-packed. The provenance is the real
``fingerprint.so_provenance_members`` of stage ``flex_so_bond`` whose
``fp_green`` is the digest of the case's companion Green function;
:func:`write_companion_green` writes that Green function as a ``green``
archive the fingerprint check accepts.

The writer, the reader and the cross-refusal tests of the archive share
these inputs.
"""
import types
from typing import NamedTuple

import numpy as np

from hwave.solver import fingerprint as fp
from hwave.solver import so_bond
from hwave.solver.so_records import RECORD_DTYPE

#: the producer label the synthetic provenance records
LABEL = "hwave FLEX (spin-orbital, bond channels)"


class ArchiveInputs(NamedTuple):
    """``kwargs`` of the writer (without ``path``, ``output_layout`` and
    ``log``), the companion ``green`` ``(1, nmat, nvol, nd, nd)``, the
    dense ``chi_w`` / ``gamma_b`` / ``chibar_w`` arrays (``chibar_w`` None
    when not requested), the case's sector map and ``run``, the run-side
    keyword arguments of ``so_bond_archive.load_so_bond_archive``
    (``run_band``, ``run_state``, ``run_interaction``, ``run_bond``) of the
    same case."""
    kwargs: dict
    green: np.ndarray
    chi_w: np.ndarray
    gamma_b: np.ndarray
    chibar_w: object
    sector_map: object
    run: dict = None


def _interaction(v=1.0):
    r = np.zeros(1, dtype=RECORD_DTYPE)
    r[0]["type"] = "CoulombIntra"
    r[0]["modes"] = (0, 1, 1, 0)
    r[0]["v"] = v
    return fp.interaction_inputs(types.SimpleNamespace(
        records=r, declared_topology=np.zeros((1, 3), np.int64)))


def bond_inputs(B):
    """``B`` (odd) channels along x: 0, +1, -1, +2, -2, ... and their reverse."""
    rows = [0] + [s * k for k in range(1, B // 2 + 1) for s in (1, -1)]
    dr = np.zeros((B, 3), np.int64)
    dr[:, 0] = rows
    rev = np.array([0] + [i + 1 if i % 2 else i - 1 for i in range(1, B)], np.int64)
    return fp.BondInputs(dr, rev, so_bond.SLOT_CONVENTION.to_string(), 1)


def block_diagonal(sm, lead, rng):
    """A random code-slot array ``lead + (ND, ND)`` without off-sector entries."""
    M = np.zeros(tuple(lead) + (sm.ND, sm.ND), np.complex128)
    for s in range(3):
        idx = sm.index[sm.offsets[s]:sm.offsets[s + 1]]
        n = len(idx)
        M[..., idx[:, None], idx] = (rng.normal(size=tuple(lead) + (n, n))
                                     + 1j * rng.normal(size=tuple(lead) + (n, n)))
    return M


def wavevector_index(cell):
    nx, ny, nz = cell
    return np.array([(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)],
                    dtype=np.int64).reshape(nx * ny * nz, 3)


def archive_inputs(*, block_layout="dense", B=3, norb_phys=1, cell=(4, 1, 1), nmat=4,
                   seed=0, chibar=False, green_file="green.npz", beta=2.0,
                   scf_converged=True, chi_fill=None, prov_overrides=None):
    """The writer inputs of a synthetic case (see the module docstring).
    ``chi_fill`` (a number, dense layout only) fills ``chi_w`` with that
    value instead of random entries (a large array built fast).
    ``prov_overrides`` replaces keyword arguments of the archive's
    ``fingerprint.so_provenance_members`` call only (``transfer``,
    ``interaction``, ``bond``, ...); ``run`` stays the case's."""
    rng = np.random.default_rng(seed)
    nd = 2 * norb_phys
    cell = tuple(int(n) for n in cell)
    nvol = int(np.prod(cell))
    bond = bond_inputs(B)
    sm = so_bond.sz_sector_map(types.SimpleNamespace(reverse=bond.reverse), nd)
    green = (rng.normal(size=(1, nmat, nvol, nd, nd))
             + 1j * rng.normal(size=(1, nmat, nvol, nd, nd)))
    meta = fp.BandMeta(0.0, True, "spin_block", nd, norb_phys, cell, cell)
    state = fp.so_state(beta=beta, nmat=nmat, mu=-0.3, density=0.8, calc_mu=True,
                        target=0.8)
    transfer = {((0, 0, 0), (0, 0)): -1.0 + 0j}
    prov_kw = dict(
        producer_stage="flex_so_bond", producer_label=LABEL,
        transfer=transfer, extern=None, meta=meta,
        fp_green=fp.green_digest(green, fp.GREEN_LAYOUT), state=state,
        interaction=_interaction(), bond=bond)
    prov_kw.update(prov_overrides or {})
    prov = fp.so_provenance_members(**prov_kw)
    run = dict(run_band=fp.band_inputs(transfer, None, meta), run_state=state,
               run_interaction=_interaction(), run_bond=bond)
    if chi_fill is None:
        chi = block_diagonal(sm, (nmat, nvol), rng)
    else:
        chi = np.full((nmat, nvol, sm.ND, sm.ND), chi_fill, dtype=np.complex128)
    gam = block_diagonal(sm, (nvol,), rng)
    chibar_w = block_diagonal(sm, (nmat, nvol), rng) if chibar else None
    packed = block_layout == "sz_packed"

    def store(a):
        if a is None or not packed:
            return a
        return so_bond.PackedArrays.from_dense(a, sm, source="fixture")

    kwargs = dict(
        chi_w=store(chi), gamma_b=store(gam), block_layout=block_layout,
        sector_map=sm if packed else None,
        scalars=dict(norb_phys=norb_phys, nd=nd, nmat=nmat, B=B,
                     cell_shape=np.array(cell, np.int64), beta=float(beta),
                     types=np.array(["CoulombIntra"], dtype="<U12")),
        provenance=prov, green_file=green_file,
        momentum=dict(momentum_convention="e_plus_ikR",
                      wavevector_unit=2 * np.pi * np.eye(3),
                      wavevector_index=wavevector_index(cell)),
        scf_converged=scf_converged)
    if chibar:
        kwargs["chibar_w"] = store(chibar_w)
    return ArchiveInputs(kwargs, green, chi, gam, chibar_w, sm, run)


def companion_green_members(inputs, green=None):
    """The members of the ``green`` archive of ``inputs`` (its provenance,
    the markers, ``block_layout`` and the Green function ``green``, by
    default the case's own)."""
    kw = inputs.kwargs
    d = dict(kw["provenance"])
    d["spin_orbital"] = np.bool_(True)
    d["index_convention"] = np.str_("spin_block")
    d["block_layout"] = np.array(kw["block_layout"], dtype="<U16")
    d["green"] = inputs.green if green is None else green
    d["green_layout"] = np.str_(fp.GREEN_LAYOUT)
    return d


def write_companion_green(path, inputs, green=None):
    """Write :func:`companion_green_members` to ``path`` (``np.savez``)."""
    np.savez(path, **companion_green_members(inputs, green))
    return path
