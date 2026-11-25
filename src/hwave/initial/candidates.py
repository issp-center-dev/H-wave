"""Candidate initial states: ``Candidate``, ``CandidateRun``, ``CandidateSet``
and their deterministic identifiers.

* A :class:`Candidate` names a **pattern** (a :class:`PatternSpec`). Its q,
  subshape, phase, sign and seed are read-only properties of the pattern,
  never stored twice.
* A :class:`CandidateRun` is the executable unit: one candidate evaluated
  on one ``SubShape``. Its ``run_id`` is derived,
  ``<candidate_id>-b<Bx>x<By>x<Bz>``, so it cannot disagree with the run.
* A :class:`CandidateSet` holds the model facts, the tier, the cells, the
  candidates, the runs, the report rows and the limits, with an exact JSON
  round trip (it is handed between processes as a JSON file).

Identifiers: ``candidate_id = <family>-q<num_den...>-<label slug>-r<rank>-
phi<k>[-neg]`` (random rows ``random-s<seed>``); a collision gets ``-2``,
``-3`` in call order (:func:`make_candidate_id`).
"""
import dataclasses
import math
import numbers
import operator
import re
import unicodedata
from dataclasses import dataclass
from types import MappingProxyType
from typing import Mapping, Optional, Tuple

from hwave.initial.patterns import (
    PatternSpec, _int_triple, _is_int, q_fractions)

__all__ = [
    "FAMILIES", "TIERS", "DEFAULT_AMPLITUDE", "CANDIDATE_SET_VERSION",
    "Candidate", "CandidateRun", "CandidateSet", "make_candidate_id",
    "label_slug",
]

#: the candidate families
FAMILIES = ("para", "ferro", "random", "instability", "baseline", "orbital",
            "user", "literature")

#: the tiers a CandidateSet can carry ("auto" is resolved before)
TIERS = ("large", "medium", "small")

#: the default requested amplitude (fraction of the saturation response)
DEFAULT_AMPLITUDE = 0.5

#: the version of the CandidateSet JSON document
CANDIDATE_SET_VERSION = 1

_CANDIDATE_KEYS = ("candidate_id", "family", "sources", "pattern",
                   "amplitude_requested", "citation", "label")
_SET_KEYS = ("candidate_set_version", "model", "tier", "cells", "candidates",
             "runs", "rows", "limits")


# ---------------------------------------------------------------------------
# identifiers

def label_slug(label: Optional[str]) -> str:
    """The id slug of a label: accents folded to ASCII, lower case, every
    other run of characters a single ``-``. ``None`` or an empty result
    gives ``"pattern"``."""
    if label is None:
        return "pattern"
    if not isinstance(label, str):
        raise ValueError(f"label: must be a string or None, got {label!r}")
    text = unicodedata.normalize("NFKD", label)
    text = text.encode("ascii", "ignore").decode("ascii").lower()
    slug = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    return slug or "pattern"


def _q_token(q):
    parts = []
    for f in q_fractions(q):
        n = ("m" if f.numerator < 0 else "") + str(abs(f.numerator))
        parts.append(n if f.denominator == 1 else f"{n}_{f.denominator}")
    return "q" + "_".join(parts)


def make_candidate_id(family, q, label, rank, phase_index, neg, taken,
                      *, seed=None) -> str:
    """A deterministic candidate id.

    ``<family>-q<num_den...>-<label slug>-r<rank>-phi<k>[-neg]``: q is
    written component by component as reduced ``num_den`` (``1/2`` ->
    ``1_2``, ``0`` -> ``0``, a minus sign as ``m``), joined by ``_``; the
    label slug is :func:`label_slug`. For ``family = "random"`` the id is
    ``random-s<seed>`` (``seed`` required; q, label, rank, phase_index and
    neg are not used). If the id is already in ``taken``, ``-2``, ``-3``,
    ... is appended, in call order; the returned id is added to ``taken``
    (pass ``None`` to skip the bookkeeping).
    """
    if family not in FAMILIES:
        raise ValueError(f"family: must be one of {', '.join(FAMILIES)}, "
                         f"got {family!r}")
    if family == "random":
        if not _is_int(seed) or seed < 0:
            raise ValueError(f"seed: a random candidate id needs a "
                             f"non-negative integer seed, got {seed!r}")
        base = f"random-s{int(seed)}"
    else:
        if q is None or len(q) != 3:
            raise ValueError(f"q: three rational strings required, "
                             f"got {q!r}")
        for name, v in (("rank", rank), ("phase_index", phase_index)):
            if not _is_int(v) or v < 0:
                raise ValueError(f"{name}: must be a non-negative integer, "
                                 f"got {v!r}")
        if not isinstance(neg, bool):
            raise ValueError(f"neg: must be a bool, got {neg!r}")
        base = (f"{family}-{_q_token(q)}-{label_slug(label)}"
                f"-r{int(rank)}-phi{int(phase_index)}")
        if neg:
            base += "-neg"
    if taken is None:
        return base
    cid, n = base, 2
    while cid in taken:
        cid = f"{base}-{n}"
        n += 1
    taken.add(cid)
    return cid


# ---------------------------------------------------------------------------
# JSON values held by a CandidateSet

def _frozen_json(v, path):
    """A deep read-only copy of a JSON value: objects become read-only
    mappings with string keys (in sorted order), arrays tuples; numbers
    plain ``int``/``float`` (finite). Anything else is refused."""
    if v is None or isinstance(v, (bool, str)):
        return v
    if _is_int(v):
        return int(operator.index(v))
    if isinstance(v, numbers.Real):
        out = float(v)
        if not math.isfinite(out):
            raise ValueError(f"{path}: not a JSON value (non-finite "
                             f"number {v!r})")
        return out
    if isinstance(v, Mapping):
        for key in v:
            if not isinstance(key, str):
                raise ValueError(f"{path}: not a JSON object (key {key!r} "
                                 f"is not a string)")
        # sorted keys: the JSON text does not depend on insertion order
        return MappingProxyType({key: _frozen_json(v[key], f"{path}.{key}")
                                 for key in sorted(v)})
    if isinstance(v, (list, tuple)):
        return tuple(_frozen_json(x, f"{path}[{i}]")
                     for i, x in enumerate(v))
    raise ValueError(f"{path}: not a JSON value ({type(v).__name__})")


def _thaw_json(v):
    """The plain dict/list JSON value of :func:`_frozen_json` output."""
    if isinstance(v, Mapping):
        return {key: _thaw_json(item) for key, item in v.items()}
    if isinstance(v, tuple):
        return [_thaw_json(x) for x in v]
    return v


def _check_keys(obj, keys, path):
    if not isinstance(obj, dict):
        raise ValueError(f"{path}: must be a JSON object, got {obj!r}")
    for key in obj:
        if key not in keys:
            raise ValueError(f"{path}.{key}: unknown key")
    for key in keys:
        if key not in obj:
            raise ValueError(f"{path}.{key}: required key missing")


# ---------------------------------------------------------------------------
# Candidate

@dataclass(frozen=True)
class Candidate:
    """One candidate pattern. ``citation`` is required (non-empty) for the
    literature family and must be ``None`` otherwise. ``family = "para"``
    goes with ``modulation = "none"`` and ``family = "random"`` with
    ``modulation = "random"``, both ways."""
    candidate_id: str
    family: str
    sources: Tuple[str, ...]
    pattern: PatternSpec
    amplitude_requested: float = DEFAULT_AMPLITUDE
    citation: Optional[str] = None
    label: Optional[str] = None

    def __post_init__(self):
        if not isinstance(self.candidate_id, str) or not self.candidate_id:
            raise ValueError(f"candidate_id: must be a non-empty string, "
                             f"got {self.candidate_id!r}")
        if self.family not in FAMILIES:
            raise ValueError(f"family: must be one of {', '.join(FAMILIES)}"
                             f", got {self.family!r}")
        if isinstance(self.sources, str) or not isinstance(
                self.sources, (list, tuple)):
            raise ValueError(f"sources: must be a sequence of strings, "
                             f"got {self.sources!r}")
        sources = tuple(self.sources)
        if not all(isinstance(s, str) for s in sources):
            raise ValueError(f"sources: must be strings, got {sources!r}")
        object.__setattr__(self, "sources", sources)
        if not isinstance(self.pattern, PatternSpec):
            raise ValueError(f"pattern: must be a PatternSpec, got "
                             f"{type(self.pattern).__name__}")
        a = self.amplitude_requested
        if (isinstance(a, bool) or not isinstance(a, numbers.Real)
                or not (0.0 < float(a) < 1.0)):
            raise ValueError(f"amplitude_requested: must lie in (0, 1), "
                             f"got {a!r}")
        object.__setattr__(self, "amplitude_requested", float(a))
        if self.family == "literature":
            if not isinstance(self.citation, str) or not self.citation.strip():
                raise ValueError("citation: a literature candidate needs a "
                                 "non-empty citation")
        elif self.citation is not None:
            raise ValueError(f"citation: only literature candidates carry "
                             f"a citation (family {self.family!r})")
        if self.label is not None and not isinstance(self.label, str):
            raise ValueError(f"label: must be a string or None, got "
                             f"{self.label!r}")
        mod = self.pattern.modulation
        if (self.family == "para") != (mod == "none"):
            raise ValueError(f"pattern: family 'para' and modulation 'none' "
                             f"go together (family {self.family!r}, "
                             f"modulation {mod!r})")
        if (self.family == "random") != (mod == "random"):
            raise ValueError(f"pattern: family 'random' and modulation "
                             f"'random' go together (family "
                             f"{self.family!r}, modulation {mod!r})")

    # derived, read-only
    @property
    def q(self):
        return self.pattern.q

    @property
    def subshape(self):
        return self.pattern.subshape

    @property
    def phase(self):
        return self.pattern.phase

    @property
    def sign(self):
        return self.pattern.sign

    @property
    def seed(self):
        return self.pattern.seed

    @classmethod
    def _external(cls, family, pattern, citation, label, amplitude_requested,
                  taken):
        # check everything before an id is taken from ``taken``
        cand = cls(candidate_id="pending", family=family, sources=(family,),
                   pattern=pattern, amplitude_requested=amplitude_requested,
                   citation=citation, label=label)
        cid = make_candidate_id(family, pattern.q, label, 0, 0, False, taken)
        return dataclasses.replace(cand, candidate_id=cid)

    @classmethod
    def literature(cls, *, pattern, citation, label=None,
                   amplitude_requested=DEFAULT_AMPLITUDE,
                   taken=None) -> "Candidate":
        """A literature candidate (``citation`` required). ``taken``: the
        ids already used, as in :func:`make_candidate_id`."""
        return cls._external("literature", pattern, citation, label,
                             amplitude_requested, taken)

    @classmethod
    def user(cls, *, pattern, label=None,
             amplitude_requested=DEFAULT_AMPLITUDE,
             taken=None) -> "Candidate":
        """A user candidate. ``taken`` as in :func:`make_candidate_id`."""
        return cls._external("user", pattern, None, label,
                             amplitude_requested, taken)

    def to_json(self) -> dict:
        return {
            "candidate_id": self.candidate_id,
            "family": self.family,
            "sources": list(self.sources),
            "pattern": self.pattern.to_json(),
            "amplitude_requested": self.amplitude_requested,
            "citation": self.citation,
            "label": self.label,
        }

    @classmethod
    def from_json(cls, obj, *, path="$") -> "Candidate":
        _check_keys(obj, _CANDIDATE_KEYS, path)
        try:
            pattern = PatternSpec.from_json(obj["pattern"])
        except ValueError as exc:
            raise ValueError(f"{path}.pattern: {exc}") from None
        return cls(candidate_id=obj["candidate_id"], family=obj["family"],
                   sources=obj["sources"], pattern=pattern,
                   amplitude_requested=obj["amplitude_requested"],
                   citation=obj["citation"], label=obj["label"])


# ---------------------------------------------------------------------------
# CandidateRun

@dataclass(frozen=True)
class CandidateRun:
    """One candidate evaluated on one ``SubShape``: the executable unit.
    ``run_id`` is derived, never stored."""
    candidate: Candidate
    subshape: Tuple[int, int, int]

    def __post_init__(self):
        if not isinstance(self.candidate, Candidate):
            raise ValueError(f"candidate: must be a Candidate, got "
                             f"{type(self.candidate).__name__}")
        object.__setattr__(self, "subshape", _int_triple(
            self.subshape, "subshape", positive=True))

    @property
    def run_id(self) -> str:
        bx, by, bz = self.subshape
        return f"{self.candidate.candidate_id}-b{bx}x{by}x{bz}"

    def to_json(self, *, inline: bool = True) -> dict:
        """``run_id``, ``candidate_id``, ``subshape`` and (``inline``) the
        candidate itself; a CandidateSet writes its runs with
        ``inline=False`` and refers to its candidates by id."""
        out = {"run_id": self.run_id,
               "candidate_id": self.candidate.candidate_id,
               "subshape": list(self.subshape)}
        if inline:
            out["candidate"] = self.candidate.to_json()
        return out

    @classmethod
    def from_json(cls, obj, *, candidates=None, path="$") -> "CandidateRun":
        """Read a run; without an inline ``candidate`` its ``candidate_id``
        is looked up in ``candidates`` (a mapping id -> Candidate). A
        stored ``run_id`` that differs from the derived one is refused."""
        if not isinstance(obj, dict):
            raise ValueError(f"{path}: must be a JSON object, got {obj!r}")
        allowed = ("run_id", "candidate_id", "subshape", "candidate")
        for key in obj:
            if key not in allowed:
                raise ValueError(f"{path}.{key}: unknown key")
        for key in ("candidate_id", "subshape"):
            if key not in obj:
                raise ValueError(f"{path}.{key}: required key missing")
        if "candidate" in obj:
            cand = Candidate.from_json(obj["candidate"],
                                       path=f"{path}.candidate")
            if cand.candidate_id != obj["candidate_id"]:
                raise ValueError(f"{path}.candidate_id: "
                                 f"{obj['candidate_id']!r} differs from the "
                                 f"inline candidate {cand.candidate_id!r}")
        else:
            if candidates is None or obj["candidate_id"] not in candidates:
                raise ValueError(f"{path}.candidate_id: unknown candidate "
                                 f"{obj['candidate_id']!r}")
            cand = candidates[obj["candidate_id"]]
        run = cls(cand, tuple(obj["subshape"]))
        if "run_id" in obj and obj["run_id"] != run.run_id:
            raise ValueError(f"{path}.run_id: stored {obj['run_id']!r} "
                             f"differs from the derived {run.run_id!r}")
        return run


# ---------------------------------------------------------------------------
# CandidateSet

@dataclass(frozen=True)
class CandidateSet:
    """Everything one enumeration produced. ``model`` and ``limits`` are
    JSON objects, ``cells`` and ``rows`` sequences of JSON values; they are
    stored as deep read-only copies (objects as read-only mappings, arrays
    as tuples), so the JSON round trip is exact. Every run refers to a
    candidate of the set; candidate and run ids are unique."""
    model: Mapping
    tier: str
    cells: tuple
    candidates: Tuple[Candidate, ...]
    runs: Tuple[CandidateRun, ...]
    rows: tuple
    limits: Mapping

    def __post_init__(self):
        set_ = object.__setattr__
        if self.tier not in TIERS:
            raise ValueError(f"tier: must be one of {', '.join(TIERS)}, "
                             f"got {self.tier!r}")
        for name in ("model", "limits"):
            value = getattr(self, name)
            if not isinstance(value, Mapping):
                raise ValueError(f"{name}: must be a JSON object, got "
                                 f"{type(value).__name__}")
            set_(self, name, _frozen_json(value, name))
        for name in ("cells", "rows"):
            value = getattr(self, name)
            if isinstance(value, (str, Mapping)) or not isinstance(
                    value, (list, tuple)):
                raise ValueError(f"{name}: must be a sequence, got "
                                 f"{type(value).__name__}")
            set_(self, name, _frozen_json(value, name))
        cands = tuple(self.candidates)
        by_id = {}
        for i, c in enumerate(cands):
            if not isinstance(c, Candidate):
                raise ValueError(f"candidates[{i}]: must be a Candidate")
            if c.candidate_id in by_id:
                raise ValueError(f"candidates[{i}]: duplicate candidate_id "
                                 f"{c.candidate_id!r}")
            by_id[c.candidate_id] = c
        runs = tuple(self.runs)
        run_ids = set()
        for i, r in enumerate(runs):
            if not isinstance(r, CandidateRun):
                raise ValueError(f"runs[{i}]: must be a CandidateRun")
            if by_id.get(r.candidate.candidate_id) != r.candidate:
                raise ValueError(f"runs[{i}]: the run's candidate "
                                 f"{r.candidate.candidate_id!r} is not in "
                                 f"candidates")
            if r.run_id in run_ids:
                raise ValueError(f"runs[{i}]: duplicate run_id "
                                 f"{r.run_id!r}")
            run_ids.add(r.run_id)
        set_(self, "candidates", cands)
        set_(self, "runs", runs)

    def to_json(self) -> dict:
        return {
            "candidate_set_version": CANDIDATE_SET_VERSION,
            "model": _thaw_json(self.model),
            "tier": self.tier,
            "cells": _thaw_json(self.cells),
            "candidates": [c.to_json() for c in self.candidates],
            "runs": [r.to_json(inline=False) for r in self.runs],
            "rows": _thaw_json(self.rows),
            "limits": _thaw_json(self.limits),
        }

    @classmethod
    def from_json(cls, obj) -> "CandidateSet":
        _check_keys(obj, _SET_KEYS, "$")
        version = obj["candidate_set_version"]
        if not _is_int(version) or version != CANDIDATE_SET_VERSION:
            raise ValueError(f"$.candidate_set_version: must be "
                             f"{CANDIDATE_SET_VERSION}, got {version!r}")
        if not isinstance(obj["candidates"], list):
            raise ValueError("$.candidates: must be a list")
        if not isinstance(obj["runs"], list):
            raise ValueError("$.runs: must be a list")
        cands = [Candidate.from_json(c, path=f"$.candidates[{i}]")
                 for i, c in enumerate(obj["candidates"])]
        by_id = {c.candidate_id: c for c in cands}
        runs = [CandidateRun.from_json(r, candidates=by_id,
                                       path=f"$.runs[{i}]")
                for i, r in enumerate(obj["runs"])]
        return cls(model=obj["model"], tier=obj["tier"], cells=obj["cells"],
                   candidates=cands, runs=runs, rows=obj["rows"],
                   limits=obj["limits"])
