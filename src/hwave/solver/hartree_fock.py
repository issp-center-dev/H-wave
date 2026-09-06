"""The Hartree-Fock (mean-field) kernel shared by UHFk and FLEX
(GitHub issue #181, Tier 3 Phase B; spec
docs/superpowers/specs/2026-09-06-flex-bond-sigma-phase-b-181-design.md,
section 2.2).

Two pure functions carry the interaction part of ``UHFk._make_ham_inter``
and ``UHFk._make_ham`` (NORMAL spin-major mode; UHFk's ``enable_spin_orbital``
branch is not shared in this version):

* :func:`build_interaction_tables` -- the raw pre-fold interaction dict ->
  per-type reversal-symmetrised real-space coefficient tables and the
  ``(2, 2, 2, 2)`` spin tables, plus a ``discarded`` report of every entry
  the per-type semantics drop (UHFk ignores it, FLEX refuses on it);
* :func:`accumulate_hf` -- the Hartree (+ Fock) contributions of every
  type added IN PLACE, in the original order and with the original
  expressions, onto the caller's ``(nvol, 2*norb, 2*norb)`` accumulator
  (UHFk starts it from the transfer Hamiltonian, FLEX from zero), so
  UHFk's per-iteration Hamiltonian is bit-identical to the pre-refactor
  code (tests/test_hartree_fock_kernel.py carries that code verbatim).

The FLEX-side density evaluator of the same spec lives here too
(:func:`equal_time_density`, section 2.1).
"""
from collections import namedtuple
from dataclasses import dataclass

import numpy as np

from .kgrid import reverse_fft_axes
from ..qlmsio import wan90

InteractionTables = namedtuple("InteractionTables",
                               ["inter_table", "spin_table", "discarded"])

# The fixed type order of UHFk's interaction loop (accumulate_hf) and of
# the discarded report.
HF_TYPE_ORDER = ("CoulombIntra", "CoulombInter", "Hund", "Ising",
                 "PairLift", "Exchange", "PairHop")


class NonFiniteError(FloatingPointError):
    """A numerical checkpoint found a non-finite value (spec section 6)."""


def is_hermitian_batch(a, rel=1e-10):
    """``(ok, err)`` for a batch ``(..., n, n)``: ``err = max|a - a^dag|``
    relative to ``max(1, max|a|)``; ``ok`` iff ``err <= rel``."""
    a = np.asarray(a)
    err = float(np.max(np.abs(a - np.conj(np.swapaxes(a, -1, -2)))))
    scale = max(1.0, float(np.max(np.abs(a)))) if a.size else 1.0
    err = err / scale
    return (err <= rel), err


def _reverse_closed(tab_r):
    """``(V + V_ba(-r)^*) / 2`` -- UHFk's reversal/Hermitian symmetrisation."""
    tba = np.conjugate(
        np.transpose(
            reverse_fft_axes(tab_r, (0, 1, 2)),
            (0, 1, 2, 4, 3)
        )
    )
    return (tab_r + tba) / 2


def build_interaction_tables(param_ham, norb, shape):
    """UHFk's ``_make_ham_inter`` (normal mode) as a pure function.

    Returns :class:`InteractionTables`; ``inter_table[type]`` is ``None``
    for an absent type. The ambiguous aggregate ``Coulomb`` +
    ``CoulombIntra``/``CoulombInter`` case raises ``ValueError`` (UHFk's
    wrapper turns it into its legacy ``logger.error`` + ``SystemExit``).
    ``discarded`` lists, in :data:`HF_TYPE_ORDER` and insertion order,
    every declared entry the per-type semantics drop -- today the
    ``CoulombIntra`` entries that are not ``r = 0, a == b``.
    """
    nx, ny, nz = shape
    inter_table = {}
    spin_table = {}
    discarded = []

    def _zeros():
        return np.zeros((nx, ny, nz, norb, norb), dtype=np.complex128)

    if 'Coulomb' in param_ham.keys():
        if ('CoulombIntra' in param_ham.keys()
                or 'CoulombInter' in param_ham.keys()):
            raise ValueError(
                "Coulomb cannot be specified together with "
                "CoulombIntra or CoulombInter")
        coulomb_intra, coulomb_inter = wan90.split_coulomb(
            param_ham["Coulomb"])
        uab_r = _zeros()
        vab_r = _zeros()
        for (irvec, orbvec), v in coulomb_intra.items():
            uab_r[(*irvec, *orbvec)] += v
        for (irvec, orbvec), v in coulomb_inter.items():
            vab_r[(*irvec, *orbvec)] += v
        inter_table["CoulombIntra"] = uab_r  # r=0 component
        spin_table["CoulombIntra"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["CoulombIntra"][0, 1, 1, 0] = 1
        spin_table["CoulombIntra"][1, 0, 0, 1] = 1
        inter_table["CoulombInter"] = _reverse_closed(vab_r)
        spin_table["CoulombInter"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["CoulombInter"][0, 0, 0, 0] = 1
        spin_table["CoulombInter"][1, 1, 1, 1] = 1
        spin_table["CoulombInter"][0, 1, 1, 0] = 1
        spin_table["CoulombInter"][1, 0, 0, 1] = 1
    else:
        inter_table["CoulombIntra"] = None
        inter_table["CoulombInter"] = None
        if 'CoulombIntra' in param_ham.keys():
            uab_r = _zeros()
            for (irvec, orbvec), v in param_ham["CoulombIntra"].items():
                alpha, beta = orbvec
                if irvec == (0, 0, 0) and alpha == beta:
                    uab_r[(*irvec, *orbvec)] += v
                else:
                    discarded.append(("CoulombIntra", tuple(irvec),
                                      tuple(orbvec), v))
            inter_table["CoulombIntra"] = uab_r
            spin_table["CoulombIntra"] = np.zeros((2, 2, 2, 2), dtype=int)
            spin_table["CoulombIntra"][0, 1, 1, 0] = 1
            spin_table["CoulombIntra"][1, 0, 0, 1] = 1
        if 'CoulombInter' in param_ham.keys():
            vab_r = _zeros()
            for (irvec, orbvec), v in param_ham["CoulombInter"].items():
                vab_r[(*irvec, *orbvec)] += v
            inter_table["CoulombInter"] = _reverse_closed(vab_r)
            spin_table["CoulombInter"] = np.zeros((2, 2, 2, 2), dtype=int)
            spin_table["CoulombInter"][0, 0, 0, 0] = 1
            spin_table["CoulombInter"][1, 1, 1, 1] = 1
            spin_table["CoulombInter"][0, 1, 1, 0] = 1
            spin_table["CoulombInter"][1, 0, 0, 1] = 1
    if 'Hund' in param_ham.keys():
        jab_r = _zeros()
        for (irvec, orbvec), v in param_ham["Hund"].items():
            jab_r[(*irvec, *orbvec)] += v
        inter_table["Hund"] = -_reverse_closed(jab_r)
        spin_table["Hund"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["Hund"][0, 0, 0, 0] = 1
        spin_table["Hund"][1, 1, 1, 1] = 1
    else:
        inter_table["Hund"] = None
    if 'Ising' in param_ham.keys():
        jab_r = _zeros()
        for (irvec, orbvec), v in param_ham["Ising"].items():
            jab_r[(*irvec, *orbvec)] += v
        inter_table["Ising"] = _reverse_closed(jab_r)
        spin_table["Ising"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["Ising"][0, 0, 0, 0] = 1
        spin_table["Ising"][1, 1, 1, 1] = 1
        spin_table["Ising"][0, 1, 1, 0] = -1
        spin_table["Ising"][1, 0, 0, 1] = -1
    else:
        inter_table["Ising"] = None
    if 'PairLift' in param_ham.keys():
        jab_r = _zeros()
        for (irvec, orbvec), v in param_ham["PairLift"].items():
            jab_r[(*irvec, *orbvec)] += v
        inter_table["PairLift"] = _reverse_closed(jab_r)
        spin_table["PairLift"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["PairLift"][0, 0, 1, 1] = 1
        spin_table["PairLift"][1, 1, 0, 0] = 1
    else:
        inter_table["PairLift"] = None
    if 'Exchange' in param_ham.keys():
        jab_r = _zeros()
        for (irvec, orbvec), v in param_ham["Exchange"].items():
            jab_r[(*irvec, *orbvec)] += v
        inter_table["Exchange"] = -_reverse_closed(jab_r)
        spin_table["Exchange"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["Exchange"][0, 1, 0, 1] = 1
        spin_table["Exchange"][1, 0, 1, 0] = 1
    else:
        inter_table["Exchange"] = None
    if 'PairHop' in param_ham.keys():
        jab_r = _zeros()
        for (irvec, orbvec), v in param_ham["PairHop"].items():
            jab_r[(*irvec, *orbvec)] += v
        inter_table["PairHop"] = _reverse_closed(jab_r)
        spin_table["PairHop"] = np.zeros((2, 2, 2, 2), dtype=int)
        spin_table["PairHop"][0, 1, 1, 0] = 1
        spin_table["PairHop"][1, 0, 0, 1] = 1
    else:
        inter_table["PairHop"] = None
    return InteractionTables(inter_table, spin_table, tuple(discarded))


def accumulate_hf(out, rho_so_r, inter_table, spin_table, shape, *,
                  include_fock):
    """Add the Hartree (+ Fock) mean-field terms of every type onto ``out``
    IN PLACE (UHFk's ``_make_ham`` interaction loop, normal mode, verbatim
    expressions and order). ``rho_so_r`` is the real-space equal-time
    density ``(nvol, 2, norb, 2, norb)``, ``rho[r, s, a, t, b] =
    <c^dag_{s a}(0) c_{t b}(r)>`` (UHFk's ``Green``); ``out`` is
    ``(nvol, 2*norb, 2*norb)`` complex128. Returns ``out``."""
    nx, ny, nz = shape
    nvol = int(nx * ny * nz)
    gab_r = rho_so_r
    norb_inter = gab_r.shape[2]
    nd = 2 * norb_inter
    nd_virt = nd
    gbb = np.diagonal(gab_r, axis1=2, axis2=4)[0, :, :, :]
    ham = out
    for type in ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising', 'PairLift', 'Exchange']:
        if inter_table[type] is not None:
            jab_r = inter_table[type].reshape(nvol, norb_inter, norb_inter)
            spin = spin_table[type]
            hh0 = np.einsum('uvb, suvt -> stb', gbb, spin)
            hh1 = np.einsum('rab, stb -> rsta', jab_r, hh0)
            hh2 = np.einsum('rsta, ab -> rsatb', hh1, np.eye(norb_inter, norb_inter))
            hh3 = np.sum(hh2, axis=0)  # shape: (2, norb_inter, 2, norb_inter)
            hh3_nd = np.broadcast_to(hh3.reshape(nd, nd), (nvol, nd, nd))
            ham += hh3_nd
            if include_fock:
                hh4 = np.einsum('rab, rubva, sutv -> rsatb', jab_r, gab_r, spin, optimize=True)
                hh5 = np.fft.ifftn(hh4.reshape(nx, ny, nz, nd_virt, nd_virt), axes=(0, 1, 2), norm='forward')
                ham -= hh5.reshape(nvol, nd, nd)
    for type in ['PairHop']:
        if inter_table[type] is not None:
            jab_r = inter_table[type].reshape(nvol, norb_inter, norb_inter)
            spin = spin_table[type]
            if include_fock:
                hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                hh2 = np.einsum('rvbua, sutv -> rsbta', np.conjugate(gab_r), spin)
                hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, (hh1 - hh2))
                hh4 = np.fft.ifftn(hh3.reshape(nx, ny, nz, nd_virt, nd_virt), axes=(0, 1, 2), norm='forward')
            else:
                hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, hh1)
                hh4 = np.fft.ifftn(hh3.reshape(nx, ny, nz, nd_virt, nd_virt), axes=(0, 1, 2), norm='forward')
            ham += hh4.reshape(nvol, nd, nd)
    return out
