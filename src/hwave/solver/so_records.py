"""Canonical quartic records of spin-orbital interaction declarations.

:func:`compile_interaction` turns the interaction tables of the spin-orbital
path into canonical monomials ``v c+_1 c_2 c+_3 c_4``. The modes are in the
SPIN-BLOCK index ``X = s * norb_phys + a`` (``s = 0`` up, ``1`` down,
``spin_orbital.INDEX_CONVENTION``); every monomial is translated so that
operator 1 sits in cell 0, and the cells are integer lattice vectors (no
mesh, no wrapping).

Reading of the tables
---------------------
Keys are ``((Rx, Ry, Rz), (a, b))`` with 0-based physical orbitals; the type
keys are case-insensitive. The entry ``(R, a, b)`` and its partner
``(-R, b, a)`` declare one bond. The table is read symmetrised on its
literal keys::

    T~[(R, a, b)] = (T[(R, a, b)] + P(T[(-R, b, a)])) / 2

a missing entry counting zero, ``P = conj`` for PairHop and the identity
otherwise (the content of :func:`declarations.symmetrise_dense`, without its
dense grid, on which a wrapped displacement would merge with another one).
Every entry of ``T~`` (both members of a bond) is then expanded with the
documented operator, its 1/2 and its "+ h.c." as written, orbital ``a`` in
cell 0 and ``b`` in cell ``R`` (``n_{as}[c] = c+_{as}[c] c_{as}[c]``,
``sigma = +1`` up, ``-1`` down)::

    CoulombIntra  U c+_{a up} c_{a up} c+_{a dn} c_{a dn}              (R = 0, a = b)
    CoulombInter  V/2 sum_{s,s'} n_{as}[0] n_{bs'}[R]
    Hund          -J/2 sum_s n_{as}[0] n_{bs}[R]
    Ising         J/2 sum_{s,s'} sigma_s sigma_s' n_{as}[0] n_{bs'}[R]
    Exchange      J c+_{a up}[0] c_{b up}[R] c+_{b dn}[R] c_{a dn}[0]
    PairHop       J/2 c+_{a up}[0] c_{b up}[R] c+_{a dn}[0] c_{b dn}[R]
                  + J*/2 c+_{b dn}[R] c_{a dn}[0] c+_{b up}[R] c_{a up}[0]
    PairLift      J/2 c+_{a up}[0] c_{a dn}[0] c+_{b up}[R] c_{b dn}[R]
                  + J*/2 c+_{b dn}[R] c_{b up}[R] c+_{a dn}[0] c_{a up}[0]

so the sum of the records is the interaction Hamiltonian as documented.

Remainders
----------
The particle-hole re-pairing ``c+1 c2 c+3 c4 = -c+1 c4 c+3 c2 + d23 c+1 c4
+ d34 c+1 c2`` (``d23``: modes 2 and 3 coincide in the same cell, likewise
``d34``) leaves one-body terms. They are listed in ``remainders`` and are
never part of a vertex.
"""

import logging
import math
import numbers
from dataclasses import dataclass

import numpy as np

from hwave.qlmsio.wan90 import split_coulomb
from hwave.solver.declarations import (REVERSE_ATOL, REVERSE_RTOL,
                                       validate_hermitian_closure,
                                       validate_no_onsite_same_orbital_inter)

logger = logging.getLogger(__name__)

#: one record: the canonical type name, the four spin-block modes of
#: ``c+_1 c_2 c+_3 c_4``, their integer cells (operator 1 in cell 0) and the
#: coefficient
RECORD_DTYPE = np.dtype([("type", "<U12"), ("modes", "<i8", (4,)), ("cells", "<i8", (4, 3)),
                         ("v", "<c16")])
#: one one-body remainder ``v c+_1 c_2`` (operator 1 in cell 0)
REMAINDER_DTYPE = np.dtype([("modes", "<i8", (2,)), ("cells", "<i8", (2, 3)), ("v", "<c16")])
#: a bond shell whose records all have ``|v| <= INSIGNIFICANT_RTOL * max|v|``
#: is reported at INFO level (diagnostic only)
INSIGNIFICANT_RTOL = 1e-14

#: the interaction types, in the canonical order of the expansion (the
#: documented declaration types, ``declarations.KNOWN_TYPES``)
TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange", "PairLift", "PairHop")
_CANONICAL = {t.lower(): t for t in TYPES}
_IGNORED = frozenset(("geometry", "transfer", "extern"))
#: types whose imaginary part is refused (CoulombIntra is required real by
#: the Hermitian-closure validation)
_REAL_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift")

_SAME_OPERATOR_MESSAGE = (
    "{t} coefficient at R={R} orbitals ({a}, {b}) has an imaginary part ({im:.6g}). "
    "In H-wave's declaration convention the entry and its partner (−R,b,a) multiply the "
    "same operator, so the imaginary part has no operator content; only PairHop carries a "
    "physical complex phase (complex couplings beyond PairHop: issue #302). Remove the "
    "imaginary part.")
_EXCHANGE_MESSAGE = (
    "Exchange coefficient at R={R} orbitals ({a}, {b}) has an imaginary part ({im:.6g}). "
    "The partner entry (−R,b,a) multiplies the Hermitian conjugate, so an imaginary part "
    "would be a Dzyaloshinskii-Moriya-type term; H-wave's declaration convention folds it, "
    "and that extension is tracked in issue #302. Remove the imaginary part.")


@dataclass(frozen=True)
class CompiledInteraction:
    """The compiled interaction (every array read-only).

    Attributes
    ----------
    norb_phys : int
        Physical orbitals per cell (``nd = 2 * norb_phys`` spin-block modes).
    records_all : ndarray of RECORD_DTYPE
        One record per distinct ``(type, modes, cells)``, its coefficient the
        correctly rounded sum (``math.fsum``, real and imaginary parts
        separately) of every contribution; sorted by ``(type, modes, cells)``.
        Exact zeros (declared-zero rows) are kept.
    records : ndarray of RECORD_DTYPE
        ``records_all`` without the records whose coefficient is exactly
        ``0 + 0j``.
    declared_topology : ndarray of int64, shape (n_d, 3)
        Every particle-hole displacement of every declared monomial
        (zero-valued included): for cells ``(0, c2, c3, c4)`` the direct
        pairs give ``c2`` and ``c4 - c3``, the crossed pairs ``c4`` and
        ``c2 - c3``. Closed under negation, contains 0, sorted ascending.
    remainders : ndarray of REMAINDER_DTYPE
        The one-body remainders of the crossed re-pairing of every record of
        ``records_all``, sorted by ``(modes, cells)``, equal keys summed with
        ``math.fsum``; exact zeros kept.
    """

    norb_phys: int
    records_all: np.ndarray
    records: np.ndarray
    declared_topology: np.ndarray
    remainders: np.ndarray


def _fmt_key(R, a, b):
    """``R=(Rx, Ry, Rz) orbitals (a, b)`` with 1-based orbitals."""
    return "R={} orbitals ({}, {})".format(tuple(R), a + 1, b + 1)


def _canonical_tables(tables):
    """``{canonical type: raw table}``: case-insensitive type keys, the
    aggregate ``Coulomb`` table split, Geometry/Transfer/Extern ignored."""
    by_lower = {}
    for key in tables.keys():
        by_lower.setdefault(str(key).lower(), []).append(key)
    for low, keys in by_lower.items():
        if len(keys) > 1:
            names = sorted(str(k) for k in keys)
            raise ValueError("interaction type given twice with different spelling: "
                             "{!r} and {!r}".format(names[0], names[1]))
    out = {}
    for low, (key,) in by_lower.items():
        if low in _IGNORED or low == "coulomb":
            continue
        name = _CANONICAL.get(low)
        if name is None:
            raise ValueError("unknown interaction type {!r}: the spin-orbital interaction "
                             "accepts {}".format(str(key), ", ".join(TYPES + ("Coulomb",))))
        out[name] = tables[key]
    if "coulomb" in by_lower:
        if "CoulombIntra" in out or "CoulombInter" in out:
            raise ValueError("Coulomb cannot be specified together with CoulombIntra or "
                             "CoulombInter")
        raw = _normalised("Coulomb", tables[by_lower["coulomb"][0]], None)
        intra, inter = split_coulomb(raw)
        out["CoulombIntra"], out["CoulombInter"] = intra, inter
    return out


def _normalised(name, table, norb_phys):
    """The table with integer keys ``((Rx, Ry, Rz), (a, b))``; orbitals are
    range-checked when ``norb_phys`` is given (first offending entry in
    ascending key order)."""
    out = {}
    for raw_key, v in dict(table).items():
        key = _entry_key(name, raw_key)
        if key in out:
            raise ValueError("{} declares {} twice".format(name, _fmt_key(key[0], *key[1])))
        out[key] = v
    if norb_phys is not None:
        for (R, (a, b)) in sorted(out):
            if not (0 <= a < norb_phys and 0 <= b < norb_phys):
                raise ValueError("{} coefficient at {} names an orbital outside 1..{}".format(
                    name, _fmt_key(R, a, b), norb_phys))
    return out


def _integral(x):
    """``int(x)`` for an exact integer (an integral float included); None
    for anything else (a bool, a fraction, a string, ...)."""
    if isinstance(x, (bool, np.bool_)):
        return None
    if isinstance(x, numbers.Integral):
        return int(x)
    if isinstance(x, numbers.Real) and math.isfinite(x) and float(x).is_integer():
        return int(x)
    return None


def _entry_key(name, key):
    """``((Rx, Ry, Rz), (a, b))`` as ints; a malformed key or a component
    that is not an exact integer (a fraction, a bool, a string) is refused,
    never truncated. Keys are checked in the table's own order: a malformed
    key cannot be placed in the ascending key order."""
    try:
        R, ab = key
        R, ab = tuple(R), tuple(ab)
    except (TypeError, ValueError):
        R, ab = None, None
    if R is None or len(R) != 3 or len(ab) != 2:
        raise ValueError("{}: an entry key must be ((Rx, Ry, Rz), (a, b)), got "
                         "{!r}".format(name, key))
    ints = [_integral(x) for x in R + ab]
    if any(i is None for i in ints):
        raise ValueError("{}: entry key {!r} has a component that is not an exact integer; "
                         "displacements and orbitals must be integers".format(name, key))
    return (tuple(ints[:3]), (ints[3], ints[4]))


def _signed_key(key):
    (R, (a, b)) = key
    return (*R, a, b)


def symmetrised_table(itype, table):
    """The symmetrised reading ``T~`` of one table on its literal keys.

    ``T~[(R, a, b)] = (T[(R, a, b)] + P(T[(-R, b, a)])) / 2`` over the keys
    and their partners, a missing entry counting zero, ``P = conj`` for
    PairHop and the identity otherwise; no displacement is wrapped. Keys are
    ``((Rx, Ry, Rz), (a, b))``; values are complex; the result is ordered by
    ascending key.
    """
    hermitian = itype == "PairHop"
    tbl = {}
    for raw_key, v in dict(table).items():
        key = _entry_key(itype, raw_key)
        if key in tbl:
            raise ValueError("{} declares {} twice".format(itype, _fmt_key(key[0], *key[1])))
        tbl[key] = complex(v)
    keys = set(tbl) | {(tuple(-x for x in R), (b, a)) for (R, (a, b)) in tbl}
    out = {}
    for (R, (a, b)) in sorted(keys, key=_signed_key):
        p = tbl.get((tuple(-x for x in R), (b, a)), 0j)
        out[(R, (a, b))] = 0.5 * (tbl.get((R, (a, b)), 0j) + (p.conjugate() if hermitian else p))
    return out


def _expand(name, R, a, b, w, norb):
    """The monomials ``(coefficient, modes, cells)`` of one symmetrised entry,
    operator 1 in cell 0, in the order of the operator table."""
    def X(s, orb):
        return s * norb + orb
    o = (0, 0, 0)
    mR = tuple(-x for x in R)
    if name == "CoulombIntra":
        return [(w, (X(0, a), X(0, a), X(1, a), X(1, a)), (o, o, o, o))]
    if name in ("CoulombInter", "Ising"):
        out = []
        for s, ss in ((0, 1.0), (1, -1.0)):
            for s2, ss2 in ((0, 1.0), (1, -1.0)):
                c = 0.5 * w * (ss * ss2 if name == "Ising" else 1.0)
                out.append((c, (X(s, a), X(s, a), X(s2, b), X(s2, b)), (o, o, R, R)))
        return out
    if name == "Hund":
        return [(-0.5 * w, (X(s, a), X(s, a), X(s, b), X(s, b)), (o, o, R, R)) for s in (0, 1)]
    if name == "Exchange":
        return [(w, (X(0, a), X(0, b), X(1, b), X(1, a)), (o, R, R, o))]
    if name == "PairHop":
        return [(0.5 * w, (X(0, a), X(0, b), X(1, a), X(1, b)), (o, R, o, R)),
                (0.5 * w.conjugate(), (X(1, b), X(1, a), X(0, b), X(0, a)), (o, mR, o, mR))]
    # PairLift
    return [(0.5 * w, (X(0, a), X(1, a), X(0, b), X(1, b)), (o, o, R, R)),
            (0.5 * w.conjugate(), (X(1, b), X(0, b), X(1, a), X(0, a)), (o, o, mR, mR))]


def _fsum_complex(values):
    return complex(math.fsum(v.real for v in values), math.fsum(v.imag for v in values))


def _sub(u, v):
    return tuple(x - y for x, y in zip(u, v))


def _displacements(cells):
    """The particle-hole displacements of a monomial (operator 1 in cell 0):
    direct ``c2``, ``c4 - c3``; crossed ``c4``, ``c2 - c3``."""
    c1, c2, c3, c4 = cells
    return (_sub(c2, c1), _sub(c4, c3), _sub(c4, c1), _sub(c2, c3))


def _readonly(arr):
    arr.flags.writeable = False
    return arr


def compile_interaction(tables, norb_phys: int) -> CompiledInteraction:
    """Compile spin-orbital interaction declarations into canonical records.

    Parameters
    ----------
    tables : mapping
        ``{type: {((Rx, Ry, Rz), (a, b)): value}}`` with 0-based physical
        orbitals; type keys are case-insensitive (a plain mapping must not
        carry one type under two spellings). ``Coulomb`` is split into its
        on-site same-orbital (CoulombIntra) and other (CoulombInter) entries
        and is refused together with either. ``Geometry``, ``Transfer`` and
        ``Extern`` are ignored; any other unknown type is refused.
    norb_phys : int
        Physical orbitals per cell.

    Returns
    -------
    CompiledInteraction

    Raises
    ------
    ValueError
        In this order: an unknown or doubly spelled type; a malformed key or
        a key component that is not an exact integer (named in the table's
        own order); an orbital out of range; a non-finite coefficient; a
        table that is not Hermitian-closed or an on-site same-orbital
        CoulombInter row; an imaginary part
        ``|Im X| > REVERSE_ATOL + REVERSE_RTOL * |X|`` on a CoulombInter,
        Hund, Ising, Exchange or PairLift entry (checked on the entries as
        declared). The orbital, non-finite and imaginary-part refusals name
        the first offending entry in ascending ``(Rx, Ry, Rz, a, b)`` order.
    """
    if isinstance(norb_phys, bool) or int(norb_phys) != norb_phys or int(norb_phys) < 1:
        raise ValueError("norb_phys must be a positive integer, got {!r}".format(norb_phys))
    norb = int(norb_phys)
    raw = _canonical_tables(tables)
    tabs = {name: _normalised(name, raw[name], norb) for name in TYPES if name in raw}

    # (1) non-finite coefficients, before anything is sorted or summed
    for name, tbl in tabs.items():
        for key in sorted(tbl, key=_signed_key):
            v = tbl[key]
            if not np.isfinite(complex(v)):
                raise ValueError("{} coefficient at {} is not finite ({})".format(
                    name, _fmt_key(key[0], *key[1]), v))
    # (2) Hermitian closure and the on-site same-orbital CoulombInter row
    for name, tbl in tabs.items():
        validate_hermitian_closure(name, tbl)
        validate_no_onsite_same_orbital_inter(name, tbl)
    # (3) the complex contract, on the entries as declared
    for name, tbl in tabs.items():
        if name not in _REAL_TYPES:
            continue
        for key in sorted(tbl, key=_signed_key):
            v = complex(tbl[key])
            if abs(v.imag) > REVERSE_ATOL + REVERSE_RTOL * abs(v):
                (R, (a, b)) = key
                msg = _EXCHANGE_MESSAGE if name == "Exchange" else _SAME_OPERATOR_MESSAGE
                raise ValueError(msg.format(t=name, R=R, a=a + 1, b=b + 1, im=v.imag))

    # (4) symmetrised reading, (5) expansion and canonical reduction
    contributions = {}
    for name, tbl in tabs.items():
        for (R, (a, b)), w in symmetrised_table(name, tbl).items():
            for idx, (c, modes, cells) in enumerate(_expand(name, R, a, b, w, norb)):
                key = (name, modes, cells)
                contributions.setdefault(key, []).append(((*R, a, b), idx, complex(c)))
    keys = sorted(contributions)
    records_all = np.zeros(len(keys), dtype=RECORD_DTYPE)
    topology = {(0, 0, 0)}
    for i, key in enumerate(keys):
        name, modes, cells = key
        contrib = sorted(contributions[key], key=lambda t: (t[0], t[1]))
        records_all[i] = (name, modes, cells, _fsum_complex([t[2] for t in contrib]))
        for d in _displacements(cells):
            topology.add(d)
            topology.add(tuple(-x for x in d))
    records = records_all[records_all["v"] != 0]
    declared_topology = np.array(sorted(topology), dtype=np.int64).reshape(-1, 3)

    rem = {}
    for r in records_all:
        modes = tuple(int(m) for m in r["modes"])
        cells = tuple(tuple(int(x) for x in c) for c in r["cells"])
        v = complex(r["v"])
        if modes[1] == modes[2] and cells[1] == cells[2]:
            rem.setdefault(((modes[0], modes[3]), (cells[0], cells[3])), []).append(v)
        if modes[2] == modes[3] and cells[2] == cells[3]:
            rem.setdefault(((modes[0], modes[1]), (cells[0], cells[1])), []).append(v)
    rkeys = sorted(rem)
    remainders = np.zeros(len(rkeys), dtype=REMAINDER_DTYPE)
    for i, key in enumerate(rkeys):
        remainders[i] = (key[0], key[1], _fsum_complex(sorted(rem[key], key=lambda z: (z.real, z.imag))))

    _report_insignificant_shells(records, declared_topology)
    return CompiledInteraction(norb_phys=norb,
                               records_all=_readonly(records_all),
                               records=_readonly(records),
                               declared_topology=_readonly(declared_topology),
                               remainders=_readonly(remainders))


def _report_insignificant_shells(records, declared_topology):
    """One INFO line naming the declared off-site shells ``{-d, d}`` whose
    records all have ``|v| <= INSIGNIFICANT_RTOL * max|v|`` (every shell when
    there is no nonzero record). Nothing else depends on it."""
    shells = sorted({tuple(sorted({tuple(d), tuple(-x for x in d)}))
                     for d in declared_topology.tolist() if any(d)})
    if not shells:
        return
    vmax = float(np.abs(records["v"]).max()) if len(records) else 0.0
    largest = {}
    for r in records:
        mag = abs(complex(r["v"]))
        for d in _displacements(tuple(tuple(int(x) for x in c) for c in r["cells"])):
            if any(d):
                s = tuple(sorted({d, tuple(-x for x in d)}))
                largest[s] = max(largest.get(s, 0.0), mag)
    quiet = [s for s in shells if largest.get(s, 0.0) <= INSIGNIFICANT_RTOL * vmax]
    if quiet:
        logger.info("so_records: bond shell(s) {} carry only couplings with |v| <= "
                    "{:g} * max|v|".format(quiet, INSIGNIFICANT_RTOL))
