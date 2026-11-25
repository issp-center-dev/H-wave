"""PatternSpec v1: the canonical description of a seed pattern.

A pattern is a list of bond-slot entries in the sigma (x) orbital-pair
basis, together with a wave vector q, a phase, a sign and a modulation.

* An entry ``(orb_a, orb_b, r)`` is the bond slot that multiplies
  ``c^dag_{orb_a}(R + r) c_{orb_b}(R)`` (R is the cell of the annihilation
  operator). ``orb_a`` and ``orb_b`` are global PHYSICAL orbital indices of
  the primitive cell, 0-based, as in the Wannier90/Geometry files. In
  spin-orbital mode they are still physical: the spin lives in the
  coefficients, never in an interleaved ``2*orb + spin`` index.
* The entries form an amplitude A, and the field is the anchored
  modulation ``O(R) = (e^{i(q.R+phi)} A + e^{-i(q.R+phi)} A^dag)/2``, so
  Hermiticity comes from the A^dag term. A pattern may give both ``+r``
  and ``-r`` entries; they are two entries of A, and nothing is added or
  merged. Two entries with the same ``(orb_a, orb_b, r)`` are refused.
* ``modulation = "cos"``: each entry carries the complex coefficients
  ``c0, cx, cy, cz`` (absent = 0) of the 2x2 spin block
  ``1/2 sum_mu c^mu sigma^mu``.
  ``modulation = "spiral"``: each entry carries a scalar shape ``s``, and
  the pattern a pair ``axes = (mu1, mu2)``:
  ``O(R) = cos(q.R + phi) s sigma^mu1 + sin(q.R + phi) s sigma^mu2``.
  ``modulation = "none"`` is the para pattern (no entries, q = 0,
  phase 0, sign 1). ``modulation = "random"`` has no entries but a
  ``seed``; q, phase and sign are null and the cell-resolved random
  pattern is regenerated from the seed, the ``subshape`` and the model.

The JSON form writes a complex number as ``[re, im]`` and q as three
rational strings (``"2/3"``, ``"0"``, ``"-1/2"``). The JSON schema is
published next to this module (:data:`SCHEMA_PATH`); the checks of
:func:`validate_patternspec` are the authority and include the rules a
schema cannot state (duplicate entries, finite numbers, the optional
physical-orbital bound).
"""
import cmath
import json
import math
import numbers
import operator
import os
import re
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional, Sequence, Tuple

import numpy as np

__all__ = [
    "PATTERNSPEC_VERSION", "MODULATIONS", "SPIN_AXES", "SCHEMA_PATH",
    "Entry", "PatternSpec", "validate_patternspec", "canonical_json",
    "q_fractions",
    "PHI0_TOL", "SIGN_TIE_RTOL", "NEG_TRANSLATED", "PhaseRow", "PatternRow",
    "realise", "trAA", "phi0", "phase_rows", "fix_sign", "canonical_basis",
    "cos_spec", "spiral_spec", "pattern_rows",
]

#: the PatternSpec version this module reads and writes
PATTERNSPEC_VERSION = 1

#: the allowed values of ``modulation``
MODULATIONS = ("cos", "spiral", "none", "random")

#: the allowed members of a spiral's ``axes`` pair
SPIN_AXES = ("x", "y", "z")

#: the published JSON schema of PatternSpec v1
SCHEMA_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "patternspec.schema.json")

_COEFFS = ("c0", "cx", "cy", "cz")
_ENTRY_KEYS = ("orb_a", "orb_b", "r") + _COEFFS + ("s",)
_TOP_KEYS = ("patternspec_version", "q", "subshape", "phase", "sign",
             "modulation", "axes", "entries", "seed")
_OPTIONAL_TOP_KEYS = ("axes", "seed")
_Q_RE = re.compile(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?")
_PARA_Q = ("0", "0", "0")


def canonical_json(obj) -> str:
    """The canonical JSON text of ``obj``: sorted keys, no whitespace, ASCII
    only, no NaN or infinity. Equal JSON values give equal bytes."""
    return json.dumps(obj, sort_keys=True, separators=(",", ":"),
                      ensure_ascii=True, allow_nan=False)


def q_fractions(q) -> Tuple[Fraction, Fraction, Fraction]:
    """The three components of a PatternSpec q as ``Fraction``s."""
    return tuple(Fraction(x) for x in q)


# ---------------------------------------------------------------------------
# scalar checks shared by the validator and the dataclass constructors

def _is_int(v):
    return isinstance(v, numbers.Integral) and not isinstance(v, bool)


def _is_real(v):
    return isinstance(v, numbers.Real) and not isinstance(v, bool)


def _as_int(v, path):
    if not _is_int(v):
        raise ValueError(f"{path}: must be an integer, got {v!r}")
    return int(operator.index(v))


def _as_float(v, path):
    if not _is_real(v):
        raise ValueError(f"{path}: must be a number, got {v!r}")
    out = float(v)
    if not math.isfinite(out):
        raise ValueError(f"{path}: must be finite, got {v!r}")
    return out


def _as_complex(v, path):
    if isinstance(v, bool) or not isinstance(v, numbers.Number):
        raise ValueError(f"{path}: must be a number, got {v!r}")
    out = complex(v)
    if not (math.isfinite(out.real) and math.isfinite(out.imag)):
        raise ValueError(f"{path}: must be finite, got {v!r}")
    return out


def _int_triple(v, path, *, positive=False):
    if not isinstance(v, (list, tuple)) or len(v) != 3:
        raise ValueError(f"{path}: must be a list of three integers, "
                         f"got {v!r}")
    out = tuple(_as_int(x, f"{path}[{i}]") for i, x in enumerate(v))
    for i, x in enumerate(out):
        if positive and x < 1:
            raise ValueError(f"{path}[{i}]: must be a positive integer, "
                             f"got {x}")
    return out


def _complex_pair(v, path):
    """Read a JSON complex number ``[re, im]``."""
    if not isinstance(v, (list, tuple)) or len(v) != 2:
        raise ValueError(f"{path}: a complex number must be written "
                         f"[re, im], got {v!r}")
    re_, im_ = (_as_float(x, f"{path}[{i}]") for i, x in enumerate(v))
    return complex(re_, im_)


def _complex_json(z):
    # ``+ 0.0`` maps -0.0 to 0.0: equal values give equal JSON bytes
    return [float(z.real) + 0.0, float(z.imag) + 0.0]


def _q_string(v, path):
    if not isinstance(v, str) or _Q_RE.fullmatch(v) is None:
        raise ValueError(f"{path}: must be an integer or a rational "
                         f"'p/q' string, got {v!r}")
    Fraction(v)
    return v


def _reduced_q(v, path):
    """A checked q string in its reduced form (``"2/4"`` -> ``"1/2"``,
    ``"-0"`` -> ``"0"``)."""
    return str(Fraction(_q_string(v, path)))


# ---------------------------------------------------------------------------
# the validator

def validate_patternspec(obj, *, norb_phys: Optional[int] = None) -> None:
    """Check a PatternSpec v1 JSON document; raise ``ValueError`` naming
    the JSON path of the first violation (``$`` is the document root).

    ``norb_phys`` (optional): the number of physical orbitals of the
    primitive cell. When given, ``orb_a``/``orb_b`` must be below it; this
    is how a caller with a spin-orbital model catches an interleaved
    ``2*orb + spin`` index given where a physical index is required.
    """
    if not isinstance(obj, dict):
        raise ValueError(f"$: a PatternSpec must be a JSON object, "
                         f"got {type(obj).__name__}")
    for key in obj:
        if key not in _TOP_KEYS:
            raise ValueError(f"$.{key}: unknown key (allowed: "
                             f"{', '.join(_TOP_KEYS)})")
    for key in _TOP_KEYS:
        if key not in obj and key not in _OPTIONAL_TOP_KEYS:
            raise ValueError(f"$.{key}: required key missing")
    if norb_phys is not None:
        norb_phys = _as_int(norb_phys, "norb_phys")

    version = obj["patternspec_version"]
    if not _is_int(version) or version != PATTERNSPEC_VERSION:
        raise ValueError(f"$.patternspec_version: must be "
                         f"{PATTERNSPEC_VERSION}, got {version!r}")
    modulation = obj["modulation"]
    if modulation not in MODULATIONS:
        raise ValueError(f"$.modulation: must be one of "
                         f"{', '.join(MODULATIONS)}, got {modulation!r}")
    _int_triple(obj["subshape"], "$.subshape", positive=True)

    q = obj["q"]
    if q is not None:
        if not isinstance(q, (list, tuple)) or len(q) != 3:
            raise ValueError(f"$.q: must be a list of three rational "
                             f"strings or null, got {q!r}")
        for i, x in enumerate(q):
            _q_string(x, f"$.q[{i}]")
    phase = obj["phase"]
    if phase is not None:
        _as_float(phase, "$.phase")
    sign = obj["sign"]
    if sign is not None and (not _is_int(sign) or sign not in (1, -1)):
        raise ValueError(f"$.sign: must be 1 or -1 (or null), got {sign!r}")
    axes = obj.get("axes")
    if axes is not None:
        if (not isinstance(axes, (list, tuple)) or len(axes) != 2
                or any(a not in SPIN_AXES for a in axes)
                or axes[0] == axes[1]):
            raise ValueError(f"$.axes: must be a pair of distinct members "
                             f"of x, y, z (or null), got {axes!r}")
    seed = obj.get("seed")
    if seed is not None and (not _is_int(seed) or seed < 0):
        raise ValueError(f"$.seed: must be a non-negative integer "
                         f"(or null), got {seed!r}")

    entries = obj["entries"]
    if not isinstance(entries, (list, tuple)):
        raise ValueError(f"$.entries: must be a list, got {entries!r}")
    seen = {}
    for i, e in enumerate(entries):
        path = f"$.entries[{i}]"
        _check_entry(e, path, norb_phys)
        key = (e["orb_a"], e["orb_b"], tuple(e["r"]))
        if key in seen:
            raise ValueError(
                f"{path}: duplicate of $.entries[{seen[key]}] "
                f"(orb_a, orb_b, r) = ({key[0]}, {key[1]}, {list(key[2])}); "
                f"give each bond slot once")
        seen[key] = i

    if modulation in ("cos", "spiral"):
        for key, value in (("q", q), ("phase", phase), ("sign", sign)):
            if value is None:
                raise ValueError(f"$.{key}: required (not null) for "
                                 f"modulation {modulation!r}")
        if seed is not None:
            raise ValueError("$.seed: must be null unless modulation is "
                             "'random'")
        if len(entries) == 0:
            raise ValueError(f"$.entries: must not be empty for modulation "
                             f"{modulation!r} (the para pattern is "
                             f"modulation 'none')")
    if modulation == "cos":
        if axes is not None:
            raise ValueError("$.axes: only allowed with modulation 'spiral'")
        for i, e in enumerate(entries):
            if "s" in e:
                raise ValueError(f"$.entries[{i}].s: a 'cos' entry carries "
                                 f"c0, cx, cy, cz only")
        if not any(_complex_pair(e[c], "$") != 0
                   for e in entries for c in _COEFFS if c in e):
            raise ValueError("$.entries: every coefficient is zero, so the "
                             "pattern has no amplitude (the para pattern "
                             "is modulation 'none')")
    elif modulation == "spiral":
        if axes is None:
            raise ValueError("$.axes: required for modulation 'spiral'")
        for i, e in enumerate(entries):
            if "s" not in e:
                raise ValueError(f"$.entries[{i}].s: required in a 'spiral' "
                                 f"entry")
            for c in _COEFFS:
                if c in e:
                    raise ValueError(f"$.entries[{i}].{c}: a 'spiral' entry "
                                     f"carries s only")
        if not any(_complex_pair(e["s"], "$") != 0 for e in entries):
            raise ValueError("$.entries: every s is zero, so the pattern "
                             "has no amplitude")
    elif modulation == "none":
        if len(entries) != 0:
            raise ValueError("$.entries: must be empty for modulation "
                             "'none'")
        if q is None or tuple(q) != _PARA_Q:
            raise ValueError(f"$.q: must be [\"0\", \"0\", \"0\"] for "
                             f"modulation 'none', got {q!r}")
        if phase is None or phase != 0:
            raise ValueError(f"$.phase: must be 0 for modulation 'none', "
                             f"got {phase!r}")
        if sign != 1:
            raise ValueError(f"$.sign: must be 1 for modulation 'none', "
                             f"got {sign!r}")
        if axes is not None:
            raise ValueError("$.axes: only allowed with modulation 'spiral'")
        if seed is not None:
            raise ValueError("$.seed: must be null unless modulation is "
                             "'random'")
    else:  # random
        if len(entries) != 0:
            raise ValueError("$.entries: must be empty for modulation "
                             "'random'")
        if seed is None:
            raise ValueError("$.seed: required for modulation 'random'")
        for key, value in (("q", q), ("phase", phase), ("sign", sign)):
            if value is not None:
                raise ValueError(f"$.{key}: must be null for modulation "
                                 f"'random', got {value!r}")
        if axes is not None:
            raise ValueError("$.axes: only allowed with modulation 'spiral'")


def _check_entry(e, path, norb_phys):
    if not isinstance(e, dict):
        raise ValueError(f"{path}: an entry must be a JSON object, "
                         f"got {e!r}")
    for key in e:
        if key not in _ENTRY_KEYS:
            raise ValueError(f"{path}.{key}: unknown key (allowed: "
                             f"{', '.join(_ENTRY_KEYS)})")
    for key in ("orb_a", "orb_b", "r"):
        if key not in e:
            raise ValueError(f"{path}.{key}: required key missing")
    for key in ("orb_a", "orb_b"):
        v = _as_int(e[key], f"{path}.{key}")
        if v < 0:
            raise ValueError(f"{path}.{key}: must be a 0-based physical "
                             f"orbital index (>= 0), got {v}")
        if norb_phys is not None and v >= norb_phys:
            raise ValueError(
                f"{path}.{key}: {v} is not a physical orbital index of a "
                f"cell with {norb_phys} physical orbitals (indices are "
                f"physical, 0..{norb_phys - 1}, also in spin-orbital mode; "
                f"the spin is in the coefficients)")
    _int_triple(e["r"], f"{path}.r")
    for key in _COEFFS + ("s",):
        if key in e:
            _complex_pair(e[key], f"{path}.{key}")


# ---------------------------------------------------------------------------
# the dataclasses

@dataclass(frozen=True)
class Entry:
    """One bond slot of the amplitude A: it multiplies
    ``c^dag_{orb_a}(R + r) c_{orb_b}(R)``. ``orb_a``/``orb_b`` are physical
    orbital indices (0-based). ``c0, cx, cy, cz`` are the coefficients of a
    'cos' pattern (0 = absent); ``s`` is the scalar shape of a 'spiral'
    pattern (``None`` otherwise)."""
    orb_a: int
    orb_b: int
    r: Tuple[int, int, int]
    c0: complex = 0j
    cx: complex = 0j
    cy: complex = 0j
    cz: complex = 0j
    s: Optional[complex] = None

    def __post_init__(self):
        for key in ("orb_a", "orb_b"):
            v = _as_int(getattr(self, key), f"Entry.{key}")
            if v < 0:
                raise ValueError(f"Entry.{key}: must be a 0-based physical "
                                 f"orbital index (>= 0), got {v}")
            object.__setattr__(self, key, v)
        object.__setattr__(self, "r", _int_triple(self.r, "Entry.r"))
        for key in _COEFFS:
            object.__setattr__(self, key,
                               _as_complex(getattr(self, key), f"Entry.{key}"))
        if self.s is not None:
            object.__setattr__(self, "s", _as_complex(self.s, "Entry.s"))

    @property
    def key(self):
        """``(orb_a, orb_b, r)``, the bond slot."""
        return (self.orb_a, self.orb_b, self.r)

    def to_json(self) -> dict:
        """The JSON entry: nonzero ``c*`` only (absent = 0), ``s`` when
        set."""
        out = {"orb_a": self.orb_a, "orb_b": self.orb_b, "r": list(self.r)}
        for key in _COEFFS:
            z = getattr(self, key)
            if z != 0:
                out[key] = _complex_json(z)
        if self.s is not None:
            out["s"] = _complex_json(self.s)
        return out

    @classmethod
    def from_json(cls, obj, *, path="$") -> "Entry":
        """Read a JSON entry (checked as by :func:`validate_patternspec`)."""
        _check_entry(obj, path, None)
        kw = {key: _complex_pair(obj[key], f"{path}.{key}")
              for key in _COEFFS + ("s",) if key in obj}
        return cls(obj["orb_a"], obj["orb_b"], tuple(obj["r"]), **kw)


@dataclass(frozen=True, kw_only=True)
class PatternSpec:
    """A PatternSpec v1 pattern (see the module docstring). Constructing
    one checks every v1 rule; entries are stored sorted by
    ``(orb_a, orb_b, r)``, so equal patterns compare and serialise equal
    whatever order the entries were given in."""
    patternspec_version: int = PATTERNSPEC_VERSION
    q: Optional[Tuple[str, str, str]]
    subshape: Tuple[int, int, int]
    phase: Optional[float]
    sign: Optional[int]
    modulation: str
    axes: Optional[Tuple[str, str]] = None
    entries: Tuple[Entry, ...] = field(default=())
    seed: Optional[int] = None

    def __post_init__(self):
        set_ = object.__setattr__
        set_(self, "patternspec_version",
             _as_int(self.patternspec_version, "$.patternspec_version"))
        if self.q is not None:
            if not isinstance(self.q, (list, tuple)) or len(self.q) != 3:
                raise ValueError(f"$.q: must be three rational strings or "
                                 f"None, got {self.q!r}")
            set_(self, "q", tuple(_reduced_q(x, f"$.q[{i}]")
                                  for i, x in enumerate(self.q)))
        set_(self, "subshape",
             _int_triple(self.subshape, "$.subshape", positive=True))
        if self.phase is not None:
            set_(self, "phase", _as_float(self.phase, "$.phase") + 0.0)
        if self.sign is not None:
            set_(self, "sign", _as_int(self.sign, "$.sign"))
        if not isinstance(self.modulation, str):
            raise ValueError(f"$.modulation: must be a string, got "
                             f"{self.modulation!r}")
        if self.axes is not None:
            if not isinstance(self.axes, (list, tuple)):
                raise ValueError(f"$.axes: must be a pair, got "
                                 f"{self.axes!r}")
            set_(self, "axes", tuple(self.axes))
        if isinstance(self.entries, (str, bytes)) or not hasattr(
                self.entries, "__iter__"):
            raise ValueError(f"$.entries: must be a sequence of Entry, got "
                             f"{self.entries!r}")
        entries = tuple(self.entries)
        for i, e in enumerate(entries):
            if not isinstance(e, Entry):
                raise ValueError(f"$.entries[{i}]: must be an Entry, got "
                                 f"{e!r}")
        if self.seed is not None:
            set_(self, "seed", _as_int(self.seed, "$.seed"))
        # duplicates are checked on the given order, so the message names
        # the position the caller used; then store the canonical order
        validate_patternspec(self._json(entries))
        set_(self, "entries", tuple(sorted(entries, key=lambda e: e.key)))

    def _json(self, entries):
        return {
            "patternspec_version": self.patternspec_version,
            "q": None if self.q is None else list(self.q),
            "subshape": list(self.subshape),
            "phase": self.phase,
            "sign": self.sign,
            "modulation": self.modulation,
            "axes": None if self.axes is None else list(self.axes),
            "entries": [e.to_json() for e in entries],
            "seed": self.seed,
        }

    def to_json(self) -> dict:
        """The canonical JSON document (every key written, plain Python
        values only, so :func:`canonical_json` of it is byte-stable)."""
        return self._json(self.entries)

    @classmethod
    def from_json(cls, obj, *, norb_phys: Optional[int] = None
                  ) -> "PatternSpec":
        """Read a JSON document; ``ValueError`` (naming the JSON path) on
        any v1 violation. ``norb_phys`` as in
        :func:`validate_patternspec`."""
        validate_patternspec(obj, norb_phys=norb_phys)
        entries = tuple(Entry.from_json(e, path=f"$.entries[{i}]")
                        for i, e in enumerate(obj["entries"]))
        return cls(
            patternspec_version=obj["patternspec_version"],
            q=None if obj["q"] is None else tuple(obj["q"]),
            subshape=tuple(obj["subshape"]),
            phase=obj["phase"],
            sign=obj["sign"],
            modulation=obj["modulation"],
            axes=None if obj.get("axes") is None else tuple(obj["axes"]),
            entries=entries,
            seed=obj.get("seed"),
        )


# ---------------------------------------------------------------------------
# realisation on a folded cell (design section 6.3)
#
# An amplitude A is a sequence of 'cos' entries (:class:`Entry` with the
# coefficients c0, cx, cy, cz): the bond slot ``(a, b, r)`` carries the 2x2
# spin block ``X_ab(r) = 1/2 sum_mu c^mu sigma^mu`` (rows: the spin of the
# creation operator). As a matrix on spin-orbitals ``(a s), (b t)`` this is
# ``A_{(as),(bt)}(r) = X_ab(r)[s, t]``. The anchored adjoint is
# ``(A^dag)_ab(r) = e^{-iq.r} conj(A_ba(-r))``; ``q.r`` is in units of
# 2 pi (``e^{iq.r} = exp(2 pi i sum_i q_i r_i)``).

#: ``|tr AA| <= PHI0_TOL * ||A||^2`` sets ``phi0 = 0`` (and, at q = -q,
#: adds the row ``phi0 + pi/2``)
PHI0_TOL = 1e-8

#: moduli are rounded to this fraction of the largest before a tie is
#: broken by the lowest index (sign rule, pivot rule)
SIGN_TIE_RTOL = 1e-10

#: the ``equivalent_to`` of the second row at odd n without a sign map
NEG_TRANSLATED = "-O(phi0), translated"

_PAULI = {
    "0": np.eye(2, dtype=np.complex128),
    "x": np.array([[0, 1], [1, 0]], dtype=np.complex128),
    "y": np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
    "z": np.array([[1, 0], [0, -1]], dtype=np.complex128),
}
_MU_KEYS = (("c0", "0"), ("cx", "x"), ("cy", "y"), ("cz", "z"))


def _q(q):
    from hwave.initial.cells import to_qfrac
    return to_qfrac(q)


def _qdot(q, r):
    """``q.r`` modulo 1, exact (``q`` a QFrac)."""
    return sum((x * int(y) for x, y in zip(q, r)), Fraction(0)) % 1


def _eiqr(q, r, sign=1):
    return cmath.exp(sign * 2j * math.pi * float(_qdot(q, r)))


def _amplitude(A, what="A"):
    """The 2x2 spin blocks ``{(a, b, r): X}`` of a 'cos' amplitude."""
    out = {}
    for i, e in enumerate(A):
        if not isinstance(e, Entry):
            raise ValueError(f"{what}[{i}]: must be an Entry, got {e!r}")
        if e.s is not None:
            raise ValueError(f"{what}[{i}]: an amplitude carries the 'cos' "
                             f"coefficients c0, cx, cy, cz, not s")
        if e.key in out:
            raise ValueError(f"{what}[{i}]: duplicate bond slot {e.key}")
        out[e.key] = 0.5 * sum(getattr(e, c) * _PAULI[m] for c, m in _MU_KEYS)
    return out


def _adjoint(blocks, q):
    """The anchored adjoint: ``(A^dag)_ba(-r) = e^{iq.r} conj(A_ab(r))``
    (the spin block conjugate-transposed)."""
    out = {}
    for (a, b, r), X in blocks.items():
        key = (b, a, tuple(-x for x in r))
        out[key] = out.get(key, 0) + _eiqr(q, r) * np.conj(X).T
    return out


def _coefficients(X):
    """``c^mu = tr(sigma^mu X)`` of a 2x2 block, in the order 0, x, y, z."""
    return [complex(np.trace(_PAULI[m] @ X)) for _, m in _MU_KEYS]


def trAA(A, q) -> complex:
    """The bond-slot bilinear ``sum_{ab,r} A_ab(r) A_ba(-r) e^{iq.r}``
    (spin indices summed: ``sum tr(X_ab(r) X_ba(-r)) e^{iq.r}``). It is
    the cross term of the anchored adjoint: ``sum |e^{i phi} A + e^{-i phi}
    A^dag|^2`` over the slots of one cell is ``2 ||A||^2 + 2 Re(e^{2 i
    phi} tr AA)``."""
    q = _q(q)
    blocks = _amplitude(A)
    total = 0j
    for (a, b, r), X in blocks.items():
        Y = blocks.get((b, a, tuple(-x for x in r)))
        if Y is not None:
            total += complex(np.trace(X @ Y)) * _eiqr(q, r)
    return total


def _norm2(blocks):
    return float(sum(np.sum(np.abs(X) ** 2) for X in blocks.values()))


def _degenerate(A, q):
    return abs(trAA(A, q)) <= PHI0_TOL * _norm2(_amplitude(A))


def phi0(A, q) -> float:
    """``phi0 = -1/2 arg tr AA`` (in ``[-pi/2, pi/2)``), the phase that
    maximises the Hermitian norm at R = 0; 0 when ``|tr AA| <= PHI0_TOL *
    ||A||^2``."""
    t = trAA(A, q)
    if abs(t) <= PHI0_TOL * _norm2(_amplitude(A)):
        return 0.0
    return -0.5 * cmath.phase(t) + 0.0


@dataclass(frozen=True)
class PhaseRow:
    """One emitted phase of an eigenvector: ``phase`` (phi), ``neg``
    (the -O row) and ``equivalent_to`` (``"-O(phi0), translated"`` for the
    second row at odd n without a verified sign map, else ``None``)."""
    phase: float
    neg: bool
    equivalent_to: Optional[str] = None


def phase_rows(A, q, *, sign_map: bool) -> Tuple[PhaseRow, ...]:
    """The rows of one eigenvector (design section 5.1 "Rows per
    eigenvector", section 6.3).

    ``sign_map``: whether a verified map takes O to -O. At q = -q (q = 0
    included) pass whether any verified map does (a translation included,
    :func:`hwave.initial.symmetry.find_map`); at q not equal to -q pass
    :func:`hwave.initial.symmetry.sign_map_verified` (time reversal or a
    spin rotation; a translation alone is a phase shift there).

    * q = -q: phi0, plus phi0 + pi/2 when ``|tr AA| <= PHI0_TOL ||A||^2``;
      each followed by its -O row (``neg``) unless ``sign_map``.
    * q not equal to -q, order n: exactly two rows, no -O row. Even n:
      phi0 and phi0 + pi/n. Odd n without a map: phi0 and phi0 + pi/n, the
      second being the translated -O(phi0) (``equivalent_to``). Odd n with
      a map: phi0 and phi0 + pi/(2n).
    """
    from hwave.initial.cells import is_self_conjugate, order
    qf = _q(q)
    phi = phi0(A, qf)
    if is_self_conjugate(qf):
        phases = [phi]
        if _degenerate(A, qf):
            phases.append(phi + math.pi / 2)
        rows = []
        for p in phases:
            rows.append(PhaseRow(p, False))
            if not sign_map:
                rows.append(PhaseRow(p, True))
        return tuple(rows)
    n = order(qf)
    if n % 2 == 0:
        return (PhaseRow(phi, False), PhaseRow(phi + math.pi / n, False))
    if not sign_map:
        return (PhaseRow(phi, False),
                PhaseRow(phi + math.pi / n, False, NEG_TRANSLATED))
    return (PhaseRow(phi, False), PhaseRow(phi + math.pi / (2 * n), False))


def _largest_index(moduli):
    """The index of the largest modulus after rounding to
    :data:`SIGN_TIE_RTOL` of the largest (ties: the lowest index)."""
    moduli = np.asarray(moduli, dtype=float)
    top = float(moduli.max()) if moduli.size else 0.0
    if top == 0.0:
        return None
    digits = -int(round(math.log10(SIGN_TIE_RTOL)))
    rounded = np.round(moduli / top, digits)
    return int(np.flatnonzero(rounded == rounded.max())[0])


def _negate(e):
    return Entry(e.orb_a, e.orb_b, e.r, c0=-e.c0, cx=-e.cx, cy=-e.cy, cz=-e.cz)


def fix_sign(A, q) -> Tuple[Entry, ...]:
    """``A`` or ``-A``, so that the largest-modulus coefficient of the phi0
    row at R = 0, ``O(0) = (e^{i phi0} A + e^{-i phi0} A^dag)/2``, read in
    the label basis (``c^mu_ab(r)``, slots in ``(orb_a, orb_b, r)`` order,
    mu in the order 0, x, y, z), is positive: its real part, or its
    imaginary part when the real part is below :data:`SIGN_TIE_RTOL` of
    its modulus. Ties: the lowest index after rounding the moduli to
    :data:`SIGN_TIE_RTOL` relative. This chooses between phi0 and
    phi0 + pi; the second phase row of the same eigenvector uses the same
    A and is not re-signed.

    ``q`` is needed because ``A^dag`` carries ``e^{-iq.r}`` on bond
    slots."""
    A = tuple(A)
    qf = _q(q)
    blocks = _amplitude(A)
    phi = phi0(A, qf)
    adj = _adjoint(blocks, qf)
    coeffs = []
    for key in sorted(set(blocks) | set(adj)):
        o0 = 0.5 * (cmath.exp(1j * phi) * blocks.get(key, 0)
                    + cmath.exp(-1j * phi) * adj.get(key, 0))
        coeffs.extend(_coefficients(np.asarray(o0, dtype=np.complex128)))
    i = _largest_index([abs(z) for z in coeffs])
    if i is None:
        raise ValueError("fix_sign: the amplitude is zero")
    z = coeffs[i]
    flip = z.real < 0 if abs(z.real) > SIGN_TIE_RTOL * abs(z) else z.imag < 0
    return tuple(_negate(e) for e in A) if flip else A


def canonical_basis(V) -> np.ndarray:
    """A deterministic orthonormal basis of the column span of ``V`` ``(n,
    g)`` (a degenerate eigenspace in the fixed sigma (x) orbital-pair
    basis), independent of the basis ``V`` is given in (up to rounding).

    Pivots by a hand-written greedy rule (not LAPACK ``geqp3``): ``g``
    times, the coordinate whose column of the deflated ``V^H`` has the
    largest norm (ties: the lowest index after rounding to
    :data:`SIGN_TIE_RTOL` relative), then that column is projected out.
    ``W = V (V[piv, :])^{-1}``, Gram-Schmidt in pivot order, and each
    column's phase set so that its pivot entry is real positive. The
    pattern-level phase and sign (:func:`phi0`, :func:`fix_sign`) are
    applied to each column afterwards."""
    V = np.asarray(V, dtype=np.complex128)
    if V.ndim != 2 or V.shape[1] == 0 or V.shape[1] > V.shape[0]:
        raise ValueError(f"canonical_basis: V must be (n, g) with 1 <= g <= n, "
                         f"got shape {V.shape}")
    n, g = V.shape
    R = V.conj().T.copy()
    piv = []
    for _ in range(g):
        norms = np.linalg.norm(R, axis=0)
        norms[piv] = -1.0
        i = _largest_index(np.maximum(norms, 0.0))
        if i is None or norms[i] <= 0:
            raise ValueError("canonical_basis: the columns of V are linearly "
                             "dependent")
        piv.append(i)
        u = R[:, i] / norms[i]
        R = R - np.outer(u, u.conj() @ R)
    W = V @ np.linalg.inv(V[piv, :])
    out = np.empty_like(W)
    for j in range(g):
        w = W[:, j].copy()
        for k in range(j):
            w = w - out[:, k] * np.vdot(out[:, k], w)
        w = w / np.linalg.norm(w)
        p = w[piv[j]]
        if abs(p) > 0:
            w = w * (abs(p) / p)
        out[:, j] = w
    return out


def _default_subshape(q):
    return tuple(x.denominator for x in _q(q))


def _q_strings(q):
    return tuple(str(x) for x in _q(q))


def cos_spec(A, q, *, phase, sign=1, subshape=None) -> "PatternSpec":
    """The 'cos' PatternSpec of amplitude ``A`` at ``q`` with ``phase`` and
    ``sign``; ``subshape`` defaults to the minimal diagonal cell (the
    denominators of q)."""
    return PatternSpec(q=_q_strings(q), subshape=(subshape if subshape is not None
                                                  else _default_subshape(q)),
                       phase=phase, sign=sign, modulation="cos", entries=tuple(A))


def spiral_spec(A, q, phi0, *, subshape=None) -> "PatternSpec":
    """The coplanar spiral of a spin-channel amplitude ``A = s (x) sigma^z``
    (the S^z sector; ``s_ab(r) = c^z_ab(r)/2``), taken with its raw phase:
    ``O(R) = M_s(R; phi0) sigma^x + M_s(R; phi0 - pi/2) sigma^y`` with the
    anchored modulation ``M_s``, stored as ``modulation = "spiral"``,
    ``axes = ("x", "y")`` and ``phase = phi0``. For ``A = e^{-i phi0} s``
    with a real Hermitian ``s`` this is ``cos(q.R) s sigma^x + sin(q.R) s
    sigma^y``. ``ValueError`` when A has a charge or transverse part."""
    A = tuple(A)
    _amplitude(A)                       # entries checked
    scale = max((abs(e.cz) for e in A), default=0.0)
    if scale == 0.0 or any(max(abs(e.c0), abs(e.cx), abs(e.cy)) > SIGN_TIE_RTOL * scale
                           for e in A):
        raise ValueError("spiral_spec: the amplitude must be in the S^z sector "
                         "(A = s (x) sigma^z: only cz coefficients)")
    entries = tuple(Entry(e.orb_a, e.orb_b, e.r, s=e.cz / 2) for e in A)
    return PatternSpec(q=_q_strings(q), subshape=(subshape if subshape is not None
                                                  else _default_subshape(q)),
                       phase=phi0, sign=1, modulation="spiral", axes=("x", "y"),
                       entries=entries)


@dataclass(frozen=True)
class PatternRow:
    """One pattern of an eigenvector: ``spec``, its ``phase_row`` (``None``
    for the spiral) and ``spiral``."""
    spec: "PatternSpec"
    phase_row: Optional[PhaseRow]
    spiral: bool = False


def pattern_rows(A, q, *, sign_map: bool, spiral: bool = False,
                 subshape=None) -> Tuple[PatternRow, ...]:
    """The patterns of one eigenvector amplitude ``A`` at ``q``: the sign
    fixed once (:func:`fix_sign`), the rows of :func:`phase_rows` as 'cos'
    specs (``sign = -1`` for a ``neg`` row), and with ``spiral`` (a spin
    eigenvector under a passed SU(2) check, q not equal to -q) one spiral
    at phi0 (:func:`spiral_spec`), so 3 rows in that case."""
    from hwave.initial.cells import is_self_conjugate
    qf = _q(q)
    if spiral and is_self_conjugate(qf):
        raise ValueError("pattern_rows: a spiral needs q not equal to -q, got q = "
                         f"{_q_strings(q)}")
    A = fix_sign(A, qf)
    out = [PatternRow(cos_spec(A, qf, phase=row.phase, sign=-1 if row.neg else 1,
                               subshape=subshape), row)
           for row in phase_rows(A, qf, sign_map=sign_map)]
    if spiral:
        out.append(PatternRow(spiral_spec(A, qf, phi0(A, qf), subshape=subshape), None,
                              True))
    return tuple(out)


def _cell_dims(facts, subshape):
    from hwave.initial.cells import check_subshapes
    B = _int_triple(subshape, "subshape", positive=True)
    L = tuple(int(x) for x in facts.cellshape)
    check_subshapes([B], L)
    N = tuple(l // b for l, b in zip(L, B))
    return B, N


class _Layout:
    """The folded index of ``(spin, physical orbital, in-cell offset)``:
    ``s*norb_B + a + norb_orig*ir`` in normal mode, ``interleaved_index(a,
    s) + norb_orig*ir`` in spin-orbital mode (``ir = ix + Bx*(iy +
    By*iz)``)."""

    def __init__(self, facts, B):
        from hwave.solver import initial_facade
        self.B = B
        self.bvol = B[0] * B[1] * B[2]
        self.so = bool(facts.spin_orbital)
        self.nphys = int(facts.norb_phys)
        self.norb_orig = int(facts.norb)
        self.nd = 2 * self.nphys * self.bvol
        self._il = initial_facade.interleaved_index

    def ir(self, off):
        B = self.B
        return off[0] + B[0] * (off[1] + B[1] * off[2])

    def index(self, s, a, ir):
        if self.so:
            return int(self._il(a, s)) + self.norb_orig * ir
        return s * (self.norb_orig * self.bvol) + a + self.norb_orig * ir


def _accumulate(o_h, blocks, layout, N, q, phi, sign):
    """Add ``1/2 e^{sign i(q.R + phi)} X_ab(r)`` for every slot: the term
    ``c^dag_a(R + r) c_b(R)``, R running over the in-cell offsets of the
    annihilation operator (``e^{iq.R}`` is periodic in the cell). In the
    Hamiltonian convention ``o_h[D][aa, bb]`` couples ``aa`` in supercell 0
    with ``bb`` in supercell D, so the creation cell ``c_sup`` gives
    ``D = -c_sup``."""
    B = layout.B
    offsets = [(x, y, z) for z in range(B[2]) for y in range(B[1]) for x in range(B[0])]
    for (a, b, r), X in blocks.items():
        for rb in offsets:
            ph = 0.5 * cmath.exp(sign * 1j * (2 * math.pi * float(_qdot(q, rb)) + phi))
            pos = [rb[i] + r[i] for i in range(3)]
            c_sup = [pos[i] // B[i] for i in range(3)]
            ra = tuple(pos[i] % B[i] for i in range(3))
            D = tuple((-c_sup[i]) % N[i] for i in range(3))
            ia, ib = layout.ir(ra), layout.ir(rb)
            for s in range(2):
                i = layout.index(s, a, ia)
                for t in range(2):
                    if X[s, t] != 0:
                        o_h[D][i, layout.index(t, b, ib)] += ph * X[s, t]


def _finish(o_h, layout, normalise):
    from hwave.solver import initial_facade
    if normalise:
        n2 = float(np.sum(np.abs(o_h) ** 2)) / layout.bvol
        if n2 == 0.0:
            raise ValueError("realise: the pattern is zero on this cell")
        o_h = o_h / math.sqrt(n2)
    return initial_facade.to_k(o_h)


def _check_orbitals(spec, nphys):
    for i, e in enumerate(spec.entries):
        for key in ("orb_a", "orb_b"):
            if getattr(e, key) >= nphys:
                raise ValueError(
                    f"$.entries[{i}].{key}: {getattr(e, key)} is not a physical "
                    f"orbital index of a cell with {nphys} physical orbitals")


def realise(spec: "PatternSpec", facts, subshape, *, normalise=True) -> np.ndarray:
    """The Hamiltonian-convention ``O(k)`` ``(nvol_B, nd_B, nd_B)`` of
    ``spec`` on the folded cell ``subshape`` of the model ``facts``
    (:class:`~hwave.initial.facts.ModelFacts`).

    ``O(k) = sum_D e^{ikD} O_H(D)`` (the solver's ``to_k``), where
    ``O_H(D)[aa, bb]`` is the coefficient of ``c^dag_aa`` in supercell 0
    and ``c_bb`` in supercell D, in the solver's index order (normal mode
    ``s*norb_B + a + norb_orig*ir``; spin-orbital mode the interleaved
    ``2*a + s`` through ``spin_orbital.interleaved_index``, plus
    ``norb_orig*ir``); the expectation of the pattern is ``sum_ab O_ab
    G_ab`` over k and orbitals.

    * ``cos``: the anchored modulation ``O(R) = (e^{i(q.R+phi)} A +
      e^{-i(q.R+phi)} A^dag)/2`` with ``(A^dag)_ab(r) = e^{-iq.r}
      conj(A_ba(-r))``, times ``sign``; an explicit ``-r`` entry is one
      more entry of A (nothing is added or merged).
    * ``spiral``: ``M_s(R; phi) sigma^mu1 + M_s(R; phi - pi/2) sigma^mu2``,
      ``M_s`` the anchored modulation of the scalar shape ``s``, times
      ``sign``.
    * ``random``: a Hermitian ``O`` on ``R = 0`` of the folded cell,
      regenerated from ``numpy.random.default_rng(seed)`` on the
      interacting orbitals ``I_c`` (spin-orbitals ``2a, 2a + 1`` of each
      interacting physical orbital in spin-orbital mode) within the solver
      blocks of ``facts.problem`` (spin-diagonal in 2Sz mode); ``facts``
      must be loaded on ``subshape`` and ``spec.subshape`` must be it.
    * ``none`` (para): zeros.

    With ``normalise`` (default) ``||O||_F = 1`` per primitive cell
    (``sum_D ||O_H(D)||_F^2 / Bvol``). ``ValueError`` when q is not
    commensurate with ``subshape``, ``subshape`` does not divide
    ``CellShape``, an orbital index is out of range, or the pattern is zero.
    """
    from hwave.initial.cells import representable
    if not isinstance(spec, PatternSpec):
        raise ValueError(f"realise: spec must be a PatternSpec, got {spec!r}")
    B, N = _cell_dims(facts, subshape)
    layout = _Layout(facts, B)
    nvol = N[0] * N[1] * N[2]
    if spec.modulation == "none":
        return np.zeros((nvol, layout.nd, layout.nd), dtype=np.complex128)
    if spec.modulation == "random":
        return _realise_random(spec, facts, B, N, layout, normalise)
    q = _q(spec.q)
    ok, reason = representable(q, B)
    if not ok:
        raise ValueError(f"realise: the pattern is not representable on SubShape "
                         f"{list(B)}: {reason}")
    _check_orbitals(spec, layout.nphys)
    o_h = np.zeros(N + (layout.nd, layout.nd), dtype=np.complex128)
    if spec.modulation == "cos":
        terms = [(_amplitude(spec.entries, "$.entries"), spec.phase)]
    else:
        terms = []
        for ax, shift in zip(spec.axes, (0.0, -math.pi / 2)):
            terms.append(({e.key: e.s * _PAULI[ax] for e in spec.entries},
                          spec.phase + shift))
    for blocks, phi in terms:
        _accumulate(o_h, blocks, layout, N, q, phi, +1)
        _accumulate(o_h, _adjoint(blocks, q), layout, N, q, phi, -1)
    if spec.sign != 1:
        o_h = o_h * spec.sign
    return _finish(o_h, layout, normalise)


def _realise_random(spec, facts, B, N, layout, normalise):
    problem = facts.problem
    if tuple(spec.subshape) != tuple(B) or tuple(problem.subshape) != tuple(B):
        raise ValueError(
            f"realise: a random pattern is cell-resolved and runs only on its own "
            f"SubShape {list(spec.subshape)}; got subshape {list(B)} and a model "
            f"loaded on {list(problem.subshape)}")
    if problem.nd != layout.nd:
        raise ValueError("realise: the model's folded cell does not match subshape")
    idx = sorted(layout.index(s, a, ir) for ir in range(layout.bvol)
                 for a in facts.interacting_orbitals for s in range(2))
    nd = layout.nd
    block_of = np.full(nd, -1)
    for k, blk in enumerate(problem.blocks.block_info):
        block_of[np.asarray(blk, dtype=int)] = k
    rng = np.random.default_rng(spec.seed)
    n = len(idx)
    x = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    h = (x + x.conj().T) / 2
    m = np.zeros((nd, nd), dtype=np.complex128)
    if n:
        m[np.ix_(idx, idx)] = h
    m[block_of[:, None] != block_of[None, :]] = 0
    if problem.two_sz is not None and not layout.so:
        spin = np.arange(nd) // (nd // 2)
        m[spin[:, None] != spin[None, :]] = 0
    o_h = np.zeros(N + (nd, nd), dtype=np.complex128)
    o_h[0, 0, 0] = m
    return _finish(o_h, layout, normalise)
