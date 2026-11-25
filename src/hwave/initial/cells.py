"""Cells and wave vectors: commensurate q, rationalisation, representability
on a SubShape and cell proposals.

A wave vector ``q`` is a :data:`QFrac`: three ``Fraction`` components of the
reciprocal basis, reduced into ``[0, 1)``; ``q.R = sum_i q_i R_i`` for a
primitive-cell translation ``R`` in lattice coordinates. The q string of
this package uses that representative (``"2/3,1/3,0"``); a q with another
representative (for example a ``q_frac`` in ``(-1/2, 1/2]``) is compared
with :func:`same_q` (equal modulo 1).

UHFk builds a supercell from ``SubShape = (Bx, By, Bz)`` only, a DIAGONAL
cell, and refuses a ``SubShape`` that does not divide ``CellShape``. So

* the q commensurate with ``B`` are ``(m1/Bx, m2/By, m3/Bz)``;
* ``q`` is representable on ``B`` iff each reduced denominator of ``q_i``
  divides ``B_i``;
* the smallest diagonal cell for ``q`` is ``B_i = denominator(q_i)``. A
  non-diagonal cell can be smaller (``order(q)`` primitive cells, for
  example the triangular sqrt3 x sqrt3 cell at K against the diagonal
  3 x 3 cell); it needs a rebuilt unit cell and is only named
  (:class:`CellProposal` ``rotated_note``), never proposed.

No orbital index appears in this module.
"""
import math
import numbers
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Tuple

__all__ = [
    "MAX_DENOMINATOR", "RATIONALISE_TOL", "QFrac", "Rationalised",
    "CellProposal", "to_qfrac", "commensurate_q", "order",
    "is_self_conjugate", "same_q", "q_string", "rationalise",
    "representable", "propose", "check_subshapes", "lcm",
]

#: the largest denominator :func:`rationalise` accepts
MAX_DENOMINATOR = 16

#: the per-component tolerance of :func:`rationalise`
RATIONALISE_TOL = 1e-6

#: a wave vector: three reduced ``Fraction`` components in ``[0, 1)``
QFrac = Tuple[Fraction, Fraction, Fraction]

_AXES = ("x", "y", "z")

#: K and K' of a triangular lattice in the 60-degree basis a1 = (1, 0, 0),
#: a2 = (1/2, sqrt3/2, 0)
_K_60 = ((Fraction(2, 3), Fraction(1, 3), Fraction(0)),
         (Fraction(1, 3), Fraction(2, 3), Fraction(0)))

#: how many approximants :attr:`Rationalised.reason` names
_REASON_APPROXIMANTS = 2


def lcm(a: int, b: int) -> int:
    """The least common multiple of two positive integers."""
    return a * b // math.gcd(a, b)


def _component(x, path):
    """One exact rational q component, reduced into ``[0, 1)``."""
    if isinstance(x, bool):
        raise TypeError("{}: a q component must be a rational number or a "
                        "rational string, not {!r}".format(path, x))
    if isinstance(x, numbers.Rational):
        f = Fraction(x)
    elif isinstance(x, str):
        try:
            f = Fraction(x)
        except (ValueError, ZeroDivisionError):
            raise ValueError("{}: {!r} is not a rational number".format(
                path, x)) from None
        if "." in x or "e" in x.lower():
            raise ValueError("{}: {!r} must be written as an integer or a "
                             "fraction p/q".format(path, x))
    else:
        raise TypeError(
            "{}: a q component must be exact (int, Fraction or a string "
            "such as '2/3'), not {!r}; rationalise a floating-point q with "
            "rationalise()".format(path, x))
    return f % 1


def to_qfrac(q) -> QFrac:
    """``q`` (three ints, ``Fraction``s or rational strings) as a
    :data:`QFrac`: each component reduced into ``[0, 1)``. A float is
    refused (``TypeError``); use :func:`rationalise`."""
    if isinstance(q, (str, bytes)) or not hasattr(q, "__len__") \
            or len(q) != 3:
        raise ValueError("q must have three components, got {!r}".format(q))
    return tuple(_component(x, "q_" + ax) for x, ax in zip(q, _AXES))


def _shape(shape, name) -> Tuple[int, int, int]:
    """A cell shape as three positive Python ints."""
    if isinstance(shape, (str, bytes)) or not hasattr(shape, "__len__") \
            or len(shape) != 3:
        raise ValueError("{} must have three components, got {!r}".format(
            name, shape))
    out = []
    for v in shape:
        if isinstance(v, bool) or not isinstance(v, numbers.Integral) \
                or int(v) < 1:
            raise ValueError("{} must be three positive integers, got "
                             "{!r}".format(name, list(shape)))
        out.append(int(v))
    return tuple(out)


def commensurate_q(subshape) -> Tuple[QFrac, ...]:
    """Every q commensurate with the diagonal cell ``subshape``:
    ``(m1/Bx, m2/By, m3/Bz)`` with ``0 <= m_i < B_i``, in lexicographic
    order (``Bx*By*Bz`` vectors, Gamma first)."""
    bx, by, bz = _shape(subshape, "SubShape")
    return tuple((Fraction(m1, bx), Fraction(m2, by), Fraction(m3, bz))
                 for m1 in range(bx) for m2 in range(by) for m3 in range(bz))


def _denominators(q) -> Tuple[int, int, int]:
    return tuple(x.denominator for x in to_qfrac(q))


def order(q) -> int:
    """The order ``n`` of q: the least ``n >= 1`` with ``n q = 0`` modulo
    1, the lcm of the reduced denominators. It is also the number of
    primitive cells of the smallest (possibly non-diagonal) cell on which
    q is commensurate."""
    n = 1
    for d in _denominators(q):
        n = lcm(n, d)
    return n


def is_self_conjugate(q) -> bool:
    """Whether ``q = -q`` modulo 1 (``2q = 0``): every component 0 or 1/2."""
    return all(d in (1, 2) for d in _denominators(q))


def same_q(q1, q2) -> bool:
    """Whether two q are equal modulo 1 (any representatives)."""
    return to_qfrac(q1) == to_qfrac(q2)


def q_string(q) -> str:
    """The q string with the ``[0, 1)`` representative, e.g. ``"2/3,1/3,0"``."""
    return ",".join(str(x) for x in to_qfrac(q))


@dataclass(frozen=True)
class Rationalised:
    """The result of :func:`rationalise`.

    ``q`` is the exact :data:`QFrac` when every component lies within the
    tolerance of a fraction with denominator ``<= max_denominator``, else
    ``None``; then ``incommensurate`` is true and ``approximants`` holds the
    nearest admissible q (every combination of the two neighbouring
    fractions of each unmatched component, nearest first by the largest
    component deviation, then lexicographic). ``max_denominator`` and
    ``original`` record the request."""
    q: Optional[QFrac]
    approximants: Tuple[QFrac, ...]
    incommensurate: bool
    original: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    max_denominator: int = MAX_DENOMINATOR

    @property
    def reason(self) -> Optional[str]:
        """``None`` when commensurate; otherwise the message of the
        not-representable row, naming the two nearest approximants (the
        full list is :attr:`approximants`)."""
        if not self.incommensurate:
            return None
        shown = self.approximants[:_REASON_APPROXIMANTS]
        rest = len(self.approximants) - len(shown)
        more = " ({} more)".format(rest) if rest else ""
        return ("q = ({}) is not representable within max_denominator = {}"
                "; incommensurate: approximated by {}{}".format(
                    ", ".join(repr(x) for x in self.original),
                    self.max_denominator,
                    " or ".join(q_string(a) for a in shown), more))


def _neighbours(x, max_denominator):
    """The largest fraction ``<= x`` and the smallest ``>= x`` with a
    denominator ``<= max_denominator`` (``0 <= x < 1``; the upper one may
    be 1)."""
    lo, hi = Fraction(0), Fraction(1)
    xf = Fraction(x)
    for d in range(1, max_denominator + 1):
        p = math.floor(xf * d)
        a, b = Fraction(p, d), Fraction(p + 1, d)
        if a > lo:
            lo = a
        if xf == a:
            b = a
        if b < hi:
            hi = b
    return lo, hi


def rationalise(q, *, max_denominator: int = MAX_DENOMINATOR,
                tol: float = RATIONALISE_TOL) -> Rationalised:
    """Rationalise a floating-point q (three finite components, any
    representative) into a :data:`QFrac`.

    Each component is reduced into ``[0, 1)`` and replaced by the nearest
    fraction with denominator ``<= max_denominator`` when it lies within
    ``tol`` of it (absolute, per component). If any component does not, q
    is incommensurate within ``max_denominator`` and the nearest
    approximants are returned instead (:class:`Rationalised`)."""
    if isinstance(max_denominator, bool) \
            or not isinstance(max_denominator, numbers.Integral) \
            or max_denominator < 1:
        raise ValueError("max_denominator must be a positive integer, got "
                         "{!r}".format(max_denominator))
    if isinstance(tol, bool) or not isinstance(tol, numbers.Real) \
            or not math.isfinite(tol) or tol < 0:
        raise ValueError("tol must be a finite number >= 0, got "
                         "{!r}".format(tol))
    if isinstance(q, (str, bytes)) or not hasattr(q, "__len__") \
            or len(q) != 3:
        raise ValueError("q must have three components, got {!r}".format(q))
    original = []
    for x, ax in zip(q, _AXES):
        if isinstance(x, bool) or not isinstance(x, numbers.Real) \
                or not math.isfinite(float(x)):
            raise ValueError("q_{} must be a finite real number, got "
                             "{!r}".format(ax, x))
        original.append(float(x))
    max_denominator = int(max_denominator)

    exact, options = [], []
    for x in original:
        x = x % 1.0
        f = Fraction(x).limit_denominator(max_denominator)
        if abs(float(f) - x) <= tol:
            exact.append(f % 1)
            options.append(((f % 1, 0.0),))
        else:
            exact.append(None)
            lo, hi = _neighbours(x, max_denominator)
            opts = sorted({(lo % 1, abs(float(lo) - x)),
                           (hi % 1, abs(float(hi) - x))},
                          key=lambda t: (t[1], t[0]))
            options.append(tuple(opts))
    if all(f is not None for f in exact):
        return Rationalised(q=tuple(exact), approximants=(),
                            incommensurate=False, original=tuple(original),
                            max_denominator=max_denominator)
    combos = {}
    for a in options[0]:
        for b in options[1]:
            for c in options[2]:
                qq = (a[0], b[0], c[0])
                dev = max(a[1], b[1], c[1])
                combos[qq] = min(dev, combos.get(qq, dev))
    approx = tuple(sorted(combos, key=lambda qq: (combos[qq], qq)))
    return Rationalised(q=None, approximants=approx, incommensurate=True,
                        original=tuple(original),
                        max_denominator=max_denominator)


def representable(q, subshape) -> Tuple[bool, str]:
    """Whether q is commensurate with the diagonal cell ``subshape``.

    Returns ``(True, "")`` or ``(False, reason)``, the reason naming the
    first failing axis exactly: ``"q_x denominator 3 does not divide B_x =
    2"``. Other axes may fail too; the complete remedy is the cell of
    :func:`propose`, not the reason."""
    dens = _denominators(q)
    shape = _shape(subshape, "SubShape")
    for d, b, ax in zip(dens, shape, _AXES):
        if b % d != 0:
            return False, "q_{} denominator {} does not divide B_{} = {}".format(
                ax, d, ax, b)
    return True, ""


@dataclass(frozen=True)
class CellProposal:
    """The cell :func:`propose` suggests for a q.

    ``proposed_subshape`` is the smallest diagonal cell ``B_i =
    denominator(q_i)``. ``proposed_cellshape`` is ``lcm(CellShape_i, B_i)``
    per axis when ``B`` does not divide ``CellShape``, else ``None``.
    ``rotated_note`` is set when a non-diagonal cell of ``order(q)``
    primitive cells is smaller than the diagonal one; that cell needs a
    rebuilt unit cell and is not proposed."""
    proposed_subshape: Tuple[int, int, int]
    proposed_cellshape: Optional[Tuple[int, int, int]]
    rotated_note: Optional[str]


def _shape_str(shape):
    return "x".join(str(v) for v in shape)


def propose(q, cellshape) -> CellProposal:
    """The smallest diagonal cell for q and, when it does not divide
    ``cellshape``, the enlarged ``CellShape' = lcm(CellShape, B)``."""
    dens = _denominators(q)
    cell = _shape(cellshape, "CellShape")
    sub = tuple(dens)
    new_cell = None
    if any(c % b for c, b in zip(cell, sub)):
        new_cell = tuple(lcm(c, b) for c, b in zip(cell, sub))
    n = order(q)
    vol = sub[0] * sub[1] * sub[2]
    note = None
    if n < vol:
        tri = ""
        if to_qfrac(q) in _K_60:
            # the K points of a triangular lattice in the 60-degree basis
            # (a1 = (1, 0, 0), a2 = (1/2, sqrt3/2, 0)); the kernel of R -> q.R is
            # then the sqrt3 x sqrt3 sublattice. In another basis the same
            # fractions are not K, so the shape is only named conditionally.
            tri = (" (√3×√3 if the lattice is triangular with a 60° angle "
                   "between a1 and a2, where this q is a K point)")
        note = ("a non-diagonal supercell of {} primitive cells{} is smaller "
                "than the diagonal {} cell ({} primitive cells); SubShape "
                "builds only diagonal cells, so it needs a rebuilt unit cell "
                "(the hwave.model input builder or a user file) and is not "
                "proposed".format(n, tri, _shape_str(sub), vol))
    return CellProposal(proposed_subshape=sub, proposed_cellshape=new_cell,
                        rotated_note=note)


def check_subshapes(subshapes, cellshape) -> None:
    """Refuse a ``SubShape`` list the solver would refuse.

    UHFk requires every ``CellShape`` component to be a multiple of the
    matching ``SubShape`` component and otherwise aborts
    (``UHFk._init_lattice``, src/hwave/solver/uhfk.py:262-268). Raises
    ``ValueError`` quoting that rule; also for an empty list or a
    malformed shape."""
    cell = _shape(cellshape, "CellShape")
    if isinstance(subshapes, (str, bytes)) or not hasattr(subshapes, "__iter__"):
        raise ValueError("subshapes must be a list of SubShape triples, got "
                         "{!r}".format(subshapes))
    subshapes = list(subshapes)
    if not subshapes:
        raise ValueError("subshapes is empty: give at least one SubShape")
    for raw in subshapes:
        sub = _shape(raw, "SubShape")
        bad = [(ax, c, b) for ax, c, b in zip(_AXES, cell, sub) if c % b]
        if bad:
            raise ValueError(
                "SubShape {} is not compatible with CellShape {}: every "
                "CellShape component must be a multiple of the SubShape "
                "component, and the solver aborts otherwise ({})".format(
                    list(sub), list(cell), "; ".join(
                        "CellShape_{0} = {1} is not a multiple of "
                        "SubShape_{0} = {2}".format(ax, c, b)
                        for ax, c, b in bad)))
