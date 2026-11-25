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
import json
import math
import numbers
import operator
import os
import re
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Optional, Tuple

__all__ = [
    "PATTERNSPEC_VERSION", "MODULATIONS", "SPIN_AXES", "SCHEMA_PATH",
    "Entry", "PatternSpec", "validate_patternspec", "canonical_json",
    "q_fractions",
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
