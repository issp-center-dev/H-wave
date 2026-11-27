"""The candidate enumerator: reserved and baseline candidates, the runs
they give on each ``SubShape``,
deduplication under verified maps, and the per-``SubShape`` quotas.

:func:`candidate_states` is the entry point. Its steps, each a public
function:

1. :func:`reserved_candidates`: para and ferro, pattern-level (PatternSpec
   ``subshape = [1, 1, 1]``, q = 0), one candidate each (every ferro
   direction times every sign no verified map merges), and ``n_random``
   random patterns assigned round-robin over the ``SubShape`` list
   (``seed = base_seed + i``, regenerated as ``default_rng(seed)``).
2. :func:`baseline_candidates`: per ``SubShape`` B, the site-summed spin and
   charge patterns at every q commensurate with B (see **Order** below).
   The large tier has these only; the small and medium tiers add the
   orbital patterns ``n_a - n_b`` per site (family ``orbital``, baseline
   quota).
3. :func:`expand_runs`: which candidate runs on which ``SubShape``. User
   and literature candidates are evaluated on every B, para and ferro run
   once on every B, a random or generated candidate only on the B it was
   written for (its PatternSpec ``subshape``). A pair whose q is not
   commensurate with B is a not-tried row with its cell proposal.
4. The response checks (only para when the model has no
   interacting orbital or every mu-group is empty or full; the block check;
   the constraint modes), before any seeding.
5. :func:`deduplicate`: per evaluated B, under the maps of
   :func:`hwave.initial.symmetry.find_map` (a translation of that cell,
   time reversal, a spin rotation), in run order with the reserved runs
   first; a duplicate row names the kept run on the same B.
6. :func:`apply_quotas`: per ``SubShape``, on the generated rows only.

**Order.** Runs and rows are listed per ``SubShape`` (in the order given),
and within one in the order para, ferro, random, user and literature
(as given), baseline. The baseline order is
round-robin over q by sign-pattern rank: for rank r = 0, 1, ... (the r-th
representative eta of :func:`sign_patterns` in each (q, channel)), every
commensurate q in q-order (by increasing order n, so Gamma first and the
q = -q points next, then lexicographically), within one q the spin
directions (z, x, y, or z alone), then charge, then the orbital pairs
(small and medium tiers); within one sign pattern
its phase rows, the -O row of each phase, and the spiral. Every q thus
gets its low-rank patterns before any q gets a higher-rank one, so a
quota is spent evenly over q (lattice-equivalent q such as (0, 1/2) and
(1/2, 0) are treated alike). Order decides only the run order, the
display and which rows a quota cuts.

**The uniform charge at q = 0.** The all-+ sign pattern of the q = 0
charge channel is a uniform density change, fixed by N; it is skipped
before the sign-pattern cap (it takes no slot) and listed in the
``skipped_patterns`` of the SubShape's ``cells`` record with that reason
(:data:`UNIFORM_CHARGE_REASON`).

**Sign patterns and +-O.** The site sign patterns eta run over
the interacting sites of the primitive cell (site groups with at least one
interacting orbital; a site without one carries no pattern). At q = -q (q
= 0 included) the -O row of every pattern is listed. Whether a map merges
eta and -eta is decided on the first non-zero pattern of the (q, channel)
from :func:`hwave.initial.symmetry.find_map` (O to -O, translations
included): for a site-summed pattern the answer
does not depend on eta (a translation multiplies the whole pattern by
``e^{iq.R}``, time reversal and a spin rotation act on its sigma part).
When merged, the -O rows are still listed and :func:`deduplicate` turns
them into duplicate rows; they do not count toward ``max_site_patterns``.
When not merged, eta and -eta are a consecutive pair of
:func:`sign_patterns` (``both_signs=True``) and both count. At q not
equal to -q the second phase row replaces the -O row; its sign map is
:func:`hwave.initial.symmetry.sign_map_verified`.
No map is assumed: SU(2) and time reversal are the numeric checks of
:mod:`hwave.initial.symmetry`.

**Names.** Candidates are named by their seeding pattern, not by the state
it induces:
``Candidate.label`` is :func:`hwave.initial.labels.pattern_label` and the
id slug is built from its ASCII form. The sign pattern is passed as a
tuple, so a lattice name such as "Néel" is never derived.

**Cell notes.** The
``rotated_note`` of :func:`hwave.initial.cells.propose` (a smaller
non-diagonal cell exists) is shown only on a not-representable row, as
part of its proposal (``proposal_note``). A representable row, Néel on a
2x2 cell included, carries no note: the run is possible as it is.

**Indices.** Every public orbital index is a physical orbital index (the
convention of :mod:`hwave.results`); in
spin-orbital mode the patterns are realised on the interleaved layout by
:func:`hwave.initial.patterns.realise`.

**Tiers and hints.** :func:`choose_tier` picks the tier before any para
step. In the small tier the para Hartree-Fock reference and the hint
eigenpairs of :mod:`hwave.initial.hints` give the instability rows
(:func:`instability_candidates`, family ``instability``): the eigenvalues
order them (best first, :func:`hint_pairs`) and fill the instability
quota of each ``SubShape``; a pair that is not realised is still listed
with its eigenvalue and the reason. The eigenvalues never remove a row.

**T_hint** enters only the hints: of the candidate set it changes the
instability rows and nothing else (the reserved, baseline, orbital, user
and literature rows and the cells are the same at any T_hint). The seeds
of the small tier are built on the para reference at T_hint (the writer's
``para_hf`` base).
"""
import dataclasses
import itertools
import math
import numbers
import os
from dataclasses import dataclass
from fractions import Fraction
from types import MappingProxyType
from typing import Iterator, List, Mapping, Tuple

import numpy as np

from hwave.initial import cells, labels, symmetry
from hwave.initial.candidates import (
    Candidate, CandidateRun, CandidateSet, DEFAULT_AMPLITUDE, make_candidate_id)
from hwave.initial.facts import NO_SUBSHAPE_WARNING, ModelFacts, model_facts
from hwave.initial import hints as hints_mod
from hwave.initial.hints import hint_record, t_hint_bare
from hwave.initial.patterns import (
    NEG_TRANSLATED, Entry, PatternSpec, fix_sign, pattern_rows, phi0, realise)
from hwave.initial.seed import (
    BLOCK_REASON, CONSTRAINT_REASON, DEGENERATE_REASON, block_check, constraint_mode)

__all__ = [
    "Limits", "RESERVED_FAMILIES", "EXEMPT_FAMILIES", "GENERATED_FAMILIES", "PRECEDENCE",
    "NC_ZERO_NOTE", "EMPTY_FULL_NOTE", "LARGE_TIER_SKIP", "UNIFORM_CHARGE_REASON",
    "sign_patterns", "skipped_patterns",
    "interacting_sites", "baseline_candidates", "reserved_candidates", "expand_runs",
    "deduplicate", "apply_quotas", "candidate_states", "orbital_pairs",
    "HintPair", "hint_pairs", "InstabilityCandidates", "instability_candidates",
    "QUOTA_INSTABILITY", "TierChoice", "choose_tier", "small_tier_bytes",
    "SMALL_TIER_MAX_PAIRS", "SMALL_TIER_MAX_Q", "SMALL_TIER_MAX_OPERATIONS",
    "SMALL_TIER_LIMIT_WARNING",
]

#: families outside every quota, never cut, kept first by deduplication
RESERVED_FAMILIES = ("para", "ferro", "random")
#: families that are never cut and never counted
EXEMPT_FAMILIES = ("user", "literature")
#: families counted by the quotas
GENERATED_FAMILIES = ("baseline", "orbital", "instability")
#: the precedence of the families
PRECEDENCE = ("user", "literature", "para", "ferro", "random", "baseline", "orbital",
              "instability")

#: the note of a model without interacting orbitals
NC_ZERO_NOTE = ("n_c = 0: no orbital takes part in an interaction, so UHF is "
                "single-shot; only para is written")
#: the note of an empty or full model
EMPTY_FULL_NOTE = ("every mu-group is empty or full (Ncond = 0 or Ncond = Nstate): "
                   "only para is written; every pattern row is degenerate with para")
#: why the all-+ charge pattern at q = 0 is not emitted
UNIFORM_CHARGE_REASON = "uniform density change, fixed by N"
#: the hint skip reason of the large tier
LARGE_TIER_SKIP = "large tier: no instability analysis"

#: the reason of a pair whose q is not commensurate with the SubShape
_UNREPRESENTABLE = "unrepresentable"
#: the reason of a pattern row when the model has no interacting orbital
_NO_RESPONSE = "no response"

_COEFF = {"0": "c0", "x": "cx", "y": "cy", "z": "cz"}
_Q0 = ("0", "0", "0")
_UNIT = (1, 1, 1)

#: entries below this are dropped from a site-summed amplitude
_ZERO = 1e-14


# --------------------------------------------------------------------------
# limits
# --------------------------------------------------------------------------

def _check_int(name, v, low):
    if isinstance(v, bool) or not isinstance(v, numbers.Integral) or int(v) < low:
        raise ValueError("{}: must be an integer >= {}, got {!r}".format(name, low, v))
    return int(v)


@dataclass(frozen=True)
class Limits:
    """The caps of one enumeration, with their defaults. The ``**limits``
    keywords of :func:`candidate_states` are
    exactly these field names (:meth:`from_keywords`).

    ``top_k`` (eigenvectors per q, hint tiers), ``max_candidates`` (the
    generated rows of one ``SubShape``: baseline, orbital and instability,
    sign twins included), ``quota_instability`` (instability rows of one
    ``SubShape``), ``max_site_patterns`` (sign patterns per (q, channel),
    an unmerged -eta counted), ``n_random`` (random rows in total),
    ``base_seed`` (random row i uses ``seed = base_seed + i``),
    ``amplitude`` (the requested seed amplitude, in (0, 1)),
    ``max_denominator`` (the largest q denominator a cell is proposed
    for)."""
    top_k: int = 3
    max_candidates: int = 32
    quota_instability: int = 12
    max_site_patterns: int = 16
    n_random: int = 4
    base_seed: int = 0
    amplitude: float = DEFAULT_AMPLITUDE
    max_denominator: int = cells.MAX_DENOMINATOR

    def __post_init__(self):
        set_ = object.__setattr__
        for name, low in (("top_k", 1), ("max_candidates", 0), ("quota_instability", 0),
                          ("max_site_patterns", 1), ("n_random", 0), ("base_seed", 0),
                          ("max_denominator", 1)):
            set_(self, name, _check_int(name, getattr(self, name), low))
        a = self.amplitude
        if isinstance(a, bool) or not isinstance(a, numbers.Real) \
                or not math.isfinite(float(a)) or not 0.0 < float(a) < 1.0:
            raise ValueError("amplitude: must lie in (0, 1), got {!r}".format(a))
        set_(self, "amplitude", float(a))

    @classmethod
    def from_keywords(cls, **kw) -> "Limits":
        """The limits of the given keywords; ``TypeError`` naming every
        keyword that is not a field."""
        names = [f.name for f in dataclasses.fields(cls)]
        unknown = sorted(k for k in kw if k not in names)
        if unknown:
            raise TypeError("unknown limit keyword(s) {}; the accepted keywords are "
                            "{}".format(", ".join(unknown), ", ".join(names)))
        return cls(**kw)

    def to_json(self) -> dict:
        """The ``limits`` record: every field, the family precedence and
        the quota rules."""
        out = dataclasses.asdict(self)
        out["precedence"] = list(PRECEDENCE)
        out["reserved"] = list(RESERVED_FAMILIES)
        out["exempt"] = list(EXEMPT_FAMILIES)
        out["quota_scope"] = ("per SubShape, generated rows only (baseline, orbital, "
                              "instability; sign twins included); instability is cut first")
        return out


# --------------------------------------------------------------------------
# sign patterns
# --------------------------------------------------------------------------

def sign_patterns(n_sites: int, *, cap: int = 16,
                  both_signs: bool = False) -> Iterator[Tuple[int, ...]]:
    """The site sign patterns eta over ``n_sites`` sites, lazily, in a fixed
    order.

    Every class ``{eta, -eta}`` is represented by the pattern with
    ``eta[0] = +1``; the representatives are complete (all ``2^(n_sites -
    1)``) and ordered: all ``+1`` first, then by the number k of ``-1``
    entries (k = 1 ... n_sites - 1, the ``-1`` positions drawn from sites
    1 ... n_sites - 1 only), within one k lexicographically on those
    positions. With ``both_signs`` (no map merges eta and -eta, decided by
    the caller) each representative is followed at once by ``-eta``.
    Every yielded pattern counts toward ``cap``, so a cap can split a
    pair."""
    n = _check_int("n_sites", n_sites, 1)
    cap = _check_int("cap", cap, 1)
    both = bool(both_signs)

    def gen():
        count = 0
        for k in range(n):
            for pos in itertools.combinations(range(1, n), k):
                eta = [1] * n
                for p in pos:
                    eta[p] = -1
                eta = tuple(eta)
                for e in ((eta, tuple(-x for x in eta)) if both else (eta,)):
                    if count >= cap:
                        return
                    yield e
                    count += 1
    return gen()


# --------------------------------------------------------------------------
# patterns
# --------------------------------------------------------------------------

def interacting_sites(facts) -> Tuple[Tuple[int, ...], ...]:
    """The interacting orbitals of each site group that has one, in site
    order (physical orbital indices)."""
    ic = set(int(a) for a in facts.interacting_orbitals)
    out = []
    for group in facts.site_groups.groups:
        members = tuple(int(a) for a in group if int(a) in ic)
        if members:
            out.append(members)
    return tuple(out)


def _site_summed(sites, eta, mu, *, traceless):
    """``sum_site eta_site sigma^mu`` on the interacting orbitals, as on-site
    entries (``traceless``: the orbital mean is removed, the charge pattern
    at q = 0)."""
    vals = [(a, float(e)) for site, e in zip(sites, eta) for a in site]
    if traceless:
        mean = sum(v for _, v in vals) / len(vals)
        vals = [(a, v - mean) for a, v in vals]
    return tuple(Entry(a, a, (0, 0, 0), **{_COEFF[mu]: v}) for a, v in vals
                 if abs(v) > _ZERO)


def orbital_pairs(facts) -> Tuple[Tuple[Tuple[int, int], ...], ...]:
    """The orbital pairs of the orbital family: for pair rank p = 0, 1, ...
    the pair ``(a_p, a_(p+1))`` of consecutive interacting orbitals of every
    site that has at least p + 2 of them (physical orbital indices, sites
    in site order). The patterns ``n_a - n_b`` of the consecutive pairs of
    a site span its traceless orbital-diagonal densities."""
    sites = [s for s in interacting_sites(facts) if len(s) >= 2]
    out = []
    for p in range(max((len(s) for s in sites), default=1) - 1):
        out.append(tuple((s[p], s[p + 1]) for s in sites if len(s) >= p + 2))
    return tuple(out)


def _orbital_pattern(pairs, eta):
    """``sum_site eta_site (n_a - n_b)`` over the orbital pairs of the sites,
    as on-site entries."""
    out = []
    for (a, b), e in zip(pairs, eta):
        out.append(Entry(a, a, (0, 0, 0), c0=float(e)))
        out.append(Entry(b, b, (0, 0, 0), c0=-float(e)))
    return tuple(sorted(out, key=lambda x: x.key))


def _spin_axes(su2, problem):
    """z, x, y; z alone when the SU(2) check passed or 2Sz is fixed.

    z (the quantisation axis of the input's spinor basis) comes first, so a
    quota keeps the z moment of an Ising-like model before x and y. The
    per-axis residuals of the SU(2) check are not used to order the axes:
    they measure which rotations the HF map commutes with, not which
    moment direction it favours (an easy-plane exchange and an Ising term
    both commute with rotations about z)."""
    if (su2 is not None and su2.status == "passed") or not problem.sz_free:
        return ("z",)
    return ("z", "x", "y")


def _named(family, spec, facts, subshape, *, eta, spiral, rank, phase_index, neg, taken,
           amplitude, sources, patterns=None):
    """The candidate of ``spec``, named by its seeding pattern, and its
    realised ``O(k)`` on ``subshape``."""
    O = realise(spec, facts, subshape)
    lab = labels.pattern_label(O, facts, q=spec.q, eta=eta, phase=spec.phase,
                               sign=spec.sign, spiral=spiral)
    cid = make_candidate_id(family, spec.q, lab.ascii_name, rank, phase_index, neg, taken)
    cand = Candidate(candidate_id=cid, family=family, sources=tuple(sources), pattern=spec,
                     amplitude_requested=amplitude, label=lab.name)
    if patterns is not None:
        patterns[CandidateRun(cand, subshape).run_id] = O
    return cand, O


def _neg(spec):
    return dataclasses.replace(spec, sign=-spec.sign)


def _etas(n_s, cap, both, skip_uniform):
    """The sign patterns of one (q, channel): :func:`sign_patterns` with the
    all-+ class removed before the cap when ``skip_uniform``."""
    if not skip_uniform:
        return list(sign_patterns(n_s, cap=cap, both_signs=both))
    plus = (1,) * n_s
    minus = (-1,) * n_s
    it = (e for e in sign_patterns(n_s, cap=1 << n_s, both_signs=both)
          if e != plus and e != minus)
    return list(itertools.islice(it, cap))


def _channel_units(facts, B, q, kind, mu, sites, *, su2, tr, limits, taken, patterns):
    """The baseline candidates of one (q, channel) on the cell B, one list
    per sign-pattern representative in rank order. ``kind``: ``"spin"``,
    ``"charge"`` or ``"orbital"`` (then ``mu`` is the orbital-pair rank
    and ``sites`` the orbital pairs ``(a, b)`` of the sites that have
    one, family ``orbital``)."""
    p = facts.problem
    self_conj = cells.is_self_conjugate(q)
    traceless = kind == "charge" and all(x == 0 for x in q)
    n_s = len(sites)
    family = "orbital" if kind == "orbital" else "baseline"
    if kind == "orbital":
        def amp(eta):
            return _orbital_pattern(sites, eta)
    else:
        def amp(eta):
            return _site_summed(sites, eta, mu, traceless=traceless)
    # the merge of eta and -eta, decided on the first non-zero pattern
    first = None
    for eta in sign_patterns(n_s, cap=1 << (n_s - 1)):
        if traceless and eta == (1,) * n_s:
            continue
        A = amp(eta)
        if A:
            first = A
            break
    if first is None:
        return []
    merged = None
    if self_conj:
        spec0 = pattern_rows(first, q, sign_map=True, subshape=B)[0].spec
        O0 = realise(spec0, facts, B)
        merged = symmetry.find_map(O0, -O0, problem=p, su2=su2, tr=tr)
    both = self_conj and merged is None
    etas = _etas(n_s, limits.max_site_patterns, both, traceless)
    yielded = set(etas)
    spiral_ok = (kind == "spin" and not self_conj and su2 is not None
                 and su2.status == "passed")
    units = []
    rank = -1
    for eta in etas:
        if eta[0] != 1:
            continue                 # a -eta twin: the -O rows of its representative
        rank += 1
        A = amp(eta)
        if not A:
            continue
        out = []
        neg_listed = self_conj and (merged is not None or tuple(-x for x in eta) in yielded)
        flipped = fix_sign(A, q) != tuple(A)
        eta_name = tuple(-x for x in eta) if flipped else eta
        if self_conj:
            rows = pattern_rows(A, q, sign_map=True, subshape=B)
        else:
            O0 = realise(pattern_rows(A, q, sign_map=True, subshape=B)[0].spec, facts, B)
            smap = symmetry.sign_map_verified(O0, problem=p, su2=su2, tr=tr)
            rows = pattern_rows(A, q, sign_map=smap, spiral=spiral_ok, subshape=B)
        k = 0
        for row in rows:
            spec = row.spec
            index = 0 if row.spiral else k
            cand, _ = _named(family, spec, facts, B, eta=eta_name, spiral=row.spiral,
                             rank=rank, phase_index=index, neg=False, taken=taken,
                             amplitude=limits.amplitude, sources=(family,),
                             patterns=patterns)
            out.append(cand)
            if neg_listed and not row.spiral:
                twin, _ = _named(family, _neg(spec), facts, B, eta=eta_name,
                                 spiral=False, rank=rank, phase_index=index, neg=True,
                                 taken=taken, amplitude=limits.amplitude,
                                 sources=(family,), patterns=patterns)
                out.append(twin)
            if not row.spiral:
                k += 1
        units.append(out)
    return units


def _q_order(B):
    return sorted(cells.commensurate_q(B), key=lambda q: (cells.order(q), q))


def skipped_patterns(facts, subshape) -> list:
    """The patterns :func:`baseline_candidates` leaves out on purpose, as
    JSON records ``{"q", "channel", "eta", "reason"}``: the all-+ charge
    pattern at q = 0 (:data:`UNIFORM_CHARGE_REASON`)."""
    sites = interacting_sites(facts)
    if not sites:
        return []
    return [{"q": "0,0,0", "channel": "charge", "eta": [1] * len(sites),
             "reason": UNIFORM_CHARGE_REASON}]


def baseline_candidates(facts, subshape, *, su2, tr, limits=None, tier="large",
                        taken=None, patterns=None) -> List[Candidate]:
    """The baseline candidates of the cell ``subshape``: for every q
    commensurate with it, the site-summed spin patterns
    (z, x, y; z alone under a passed SU(2) check or with 2Sz fixed; plus
    one coplanar spiral per sign pattern at q not equal to -q under a
    passed SU(2) check), the site-summed charge patterns (traceless
    over the cell) and, in the small and medium tiers, the orbital
    patterns ``n_a - n_b`` (family ``orbital``, :func:`orbital_pairs`:
    uniform and with the (q, eta) signs over the sites that hold the
    pair), each in the order of the module docstring.

    ``facts`` must be loaded on ``subshape`` (``model_facts(...,
    subshape=subshape)``); the PatternSpec ``subshape`` of every candidate
    is ``subshape`` (a generated row runs only on its own cell). ``su2``,
    ``tr``: the checks of :mod:`hwave.initial.symmetry`. ``limits``:
    :class:`Limits` (``max_site_patterns``, ``amplitude``). ``tier``: the
    large tier has spin and charge only; the small and medium tiers add
    the orbital patterns. ``taken``: the ids already used
    (:func:`make_candidate_id`). ``patterns``: an optional dict that
    receives each candidate's realised ``O(k)`` keyed by its run id on
    ``subshape`` (the ``cache`` of :func:`deduplicate`)."""
    if tier not in ("large", "medium", "small"):
        raise ValueError("baseline_candidates: tier must be 'large', 'medium' or 'small', "
                         "got {!r}".format(tier))
    limits = Limits() if limits is None else limits
    B = tuple(int(x) for x in subshape)
    if tuple(facts.problem.subshape) != B:
        raise ValueError("baseline_candidates: facts are loaded on SubShape {}, not {}"
                         .format(list(facts.problem.subshape), list(B)))
    taken = set() if taken is None else taken
    sites = interacting_sites(facts)
    if not sites:
        return []
    channels = [("spin", mu, sites) for mu in _spin_axes(su2, facts.problem)]
    channels.append(("charge", "0", sites))
    if tier != "large":
        channels.extend(("orbital", p, pairs)
                        for p, pairs in enumerate(orbital_pairs(facts)))
    units = [_channel_units(facts, B, q, kind, mu, members, su2=su2, tr=tr, limits=limits,
                            taken=taken, patterns=patterns)
             for q in _q_order(B) for kind, mu, members in channels]
    out = []
    for r in range(max((len(u) for u in units), default=0)):
        for u in units:
            if r < len(u):
                out.extend(u[r])
    return out


def _para(taken):
    spec = PatternSpec(q=_Q0, subshape=_UNIT, phase=0.0, sign=1, modulation="none")
    cid = make_candidate_id("para", _Q0, "para", 0, 0, False, taken)
    return Candidate(candidate_id=cid, family="para", sources=("para",), pattern=spec,
                     label="para")


def reserved_candidates(facts, subshapes, *, su2, tr, n_random=4, base_seed=0,
                        amplitude=DEFAULT_AMPLITUDE, taken=None) -> List[Candidate]:
    """The reserved candidates: para, ferro and random.

    * para and ferro are pattern-level: PatternSpec ``subshape = [1, 1, 1]``
      and q = 0, one :class:`Candidate` each, run on every ``SubShape``
      (:func:`expand_runs`). Ferro is the uniform ``sigma^mu`` on the
      interacting orbitals, one per direction (z alone under a passed SU(2)
      check or with 2Sz fixed, else z, x, y), and its -O as a second
      candidate unless a verified map (:func:`symmetry.find_map` of O to
      -O on the cell of ``facts``) merges them; the merge is recorded in
      the candidate's ``sources`` (``"-O merged (<map>)"``).
    * ``n_random`` random patterns: row i (``seed = base_seed + i``,
      ``id random-s<seed>``) is written for ``subshapes[i % len(subshapes)]``
      (round-robin in list order) and runs only there.

    Ferro and random live on the interacting orbitals: a model without
    any (``n_c = 0``) gets para alone. ``facts``: any loaded cell (the
    ferro sign map is checked on it)."""
    taken = set() if taken is None else taken
    subshapes = [tuple(int(x) for x in B) for B in subshapes]
    if not subshapes:
        raise ValueError("subshapes is empty: give at least one SubShape")
    n_random = _check_int("n_random", n_random, 0)
    base_seed = _check_int("base_seed", base_seed, 0)
    out = [_para(taken)]
    sites = interacting_sites(facts)
    if not sites:
        return out
    B0 = tuple(int(x) for x in facts.problem.subshape)
    for mu in _spin_axes(su2, facts.problem):
        A = _site_summed(sites, (1,) * len(sites), mu, traceless=False)
        spec = PatternSpec(q=_Q0, subshape=_UNIT, phase=0.0, sign=1, modulation="cos",
                           entries=A)
        O = realise(spec, facts, B0)
        merged = symmetry.find_map(O, -O, problem=facts.problem, su2=su2, tr=tr)
        sources = ("ferro",) + (("-O merged ({})".format(merged),) if merged else ())
        cand, _ = _named("ferro", spec, facts, B0, eta=None, spiral=False, rank=0,
                         phase_index=0, neg=False, taken=taken, amplitude=amplitude,
                         sources=sources)
        out.append(cand)
        if merged is None:
            twin, _ = _named("ferro", _neg(spec), facts, B0, eta=None, spiral=False,
                             rank=0, phase_index=0, neg=True, taken=taken,
                             amplitude=amplitude, sources=("ferro",))
            out.append(twin)
    for i in range(n_random):
        seed = base_seed + i
        spec = PatternSpec(q=None, subshape=subshapes[i % len(subshapes)], phase=None,
                           sign=None, modulation="random", seed=seed)
        cid = make_candidate_id("random", None, None, 0, 0, False, taken, seed=seed)
        out.append(Candidate(candidate_id=cid, family="random", sources=("random",),
                             pattern=spec, amplitude_requested=amplitude, label="random"))
    return out


# --------------------------------------------------------------------------
# rows and runs
# --------------------------------------------------------------------------

def _channel(spec):
    if spec.modulation == "none":
        return "para"
    if spec.modulation == "random":
        return "random"
    if spec.modulation == "spiral":
        return "spin-spiral"
    used = {key for e in spec.entries for key in ("c0", "cx", "cy", "cz")
            if getattr(e, key) != 0}
    names = {"c0": "charge", "cx": "spin-x", "cy": "spin-y", "cz": "spin-z"}
    return names[used.pop()] if len(used) == 1 else "mixed"


def _equivalent_to(cand):
    """``"-O(phi0), translated"`` for the second phase row of a generated
    pattern at an odd-order q not equal to -q without a sign map (its
    phase is phi0 + pi/n), else ``None``."""
    spec = cand.pattern
    if cand.family not in ("ferro",) + GENERATED_FAMILIES or spec.modulation != "cos":
        return None
    q = cells.to_qfrac(spec.q)
    if cells.is_self_conjugate(q):
        return None
    n = cells.order(q)
    if n % 2 == 0:
        return None
    base = phi0(spec.entries, q)
    return NEG_TRANSLATED if abs(spec.phase - base - math.pi / n) <= 1e-12 else None


def _row(run, cellshape, *, representable=True, reason=None):
    c = run.candidate
    spec = c.pattern
    return {
        "run_id": run.run_id, "candidate_id": c.candidate_id, "family": c.family,
        "sources": list(c.sources),
        "q": None if spec.q is None else list(spec.q),
        "q_string": None if spec.q is None else cells.q_string(spec.q),
        "phase": spec.phase, "sign": spec.sign, "channel": _channel(spec),
        "pattern_label": c.label, "label": None, "rank": None, "lambda": None,
        "lambda_imag": None, "subshape": list(run.subshape), "cellshape": list(cellshape),
        "representable": bool(representable), "reason": reason, "proposed_subshape": None,
        "proposed_cellshape": None, "proposal_note": None,
        "equivalent_to": _equivalent_to(c), "duplicate_of": None,
        "citation": c.citation, "seed": spec.seed,
        "amplitude_requested": c.amplitude_requested,
        "tried": bool(representable),
        "not_tried_reason": None if representable else _UNREPRESENTABLE,
        "phase_marker": "pre_run",
    }


def _divides(sub, cellshape):
    return all(int(c) % int(b) == 0 for c, b in zip(cellshape, sub))


def _proposal(spec, cellshape, max_denominator):
    """``(proposed_subshape, proposed_cellshape, note)`` of a q that a
    ``SubShape`` cannot hold: the PatternSpec ``subshape`` when it divides
    ``CellShape`` and holds q, else :func:`cells.propose`."""
    q = cells.to_qfrac(spec.q)
    prop = cells.propose(q, cellshape)
    if max(x.denominator for x in q) > max_denominator:
        return None, None, ("q has a denominator above max_denominator = {}: no cell is "
                            "proposed".format(max_denominator))
    own = tuple(spec.subshape)
    if _divides(own, cellshape) and cells.representable(q, own)[0]:
        return list(own), None, prop.rotated_note
    cell = None if prop.proposed_cellshape is None else list(prop.proposed_cellshape)
    return list(prop.proposed_subshape), cell, prop.rotated_note


def _own_cell_only(cand):
    return cand.family == "random" or cand.family in GENERATED_FAMILIES


def expand_runs(cands, subshapes, *, cellshape,
                max_denominator: int = cells.MAX_DENOMINATOR):
    """The runs of ``cands`` on ``subshapes`` and one row per evaluated
    pair, in the order ``SubShape`` by ``SubShape``,
    candidates in the given order.

    User and literature candidates are evaluated on every B; para and
    ferro get one run per B; a random or generated candidate runs only on
    the B it was written for (its PatternSpec ``subshape``). A pair whose
    q is not commensurate with B is a not-tried row (``representable =
    False``, ``not_tried_reason = "unrepresentable"``) with the reason of
    :func:`cells.representable` and the proposal: the PatternSpec
    ``subshape`` when it divides ``cellshape`` (V2a), else the cell of
    :func:`cells.propose` (with ``proposed_cellshape`` when that one does
    not divide ``cellshape``); ``proposal_note`` holds the rotated-cell
    note. Runs precede deduplication, so every target has a ``run_id``.

    Returns ``(runs, rows)``; ``rows`` are JSON dicts, ``tried = True`` for
    every run until a later step says otherwise.
    """
    cellshape = tuple(int(x) for x in cellshape)
    subshapes = [tuple(int(x) for x in B) for B in subshapes]
    runs, rows = [], []
    for B in subshapes:
        for c in cands:
            if not isinstance(c, Candidate):
                raise ValueError("expand_runs: every candidate must be a Candidate, got "
                                 "{!r}".format(c))
            if _own_cell_only(c) and tuple(c.pattern.subshape) != B:
                continue
            run = CandidateRun(c, B)
            if c.pattern.q is None or c.pattern.modulation == "none":
                ok, reason = True, ""
            else:
                ok, reason = cells.representable(c.pattern.q, B)
            if ok:
                runs.append(run)
                rows.append(_row(run, cellshape))
            else:
                row = _row(run, cellshape, representable=False, reason=reason)
                sub, cell, note = _proposal(c.pattern, cellshape, max_denominator)
                row.update(proposed_subshape=sub, proposed_cellshape=cell,
                           proposal_note=note)
                rows.append(row)
    return runs, rows


def _realised(run, facts, cache):
    key = run.run_id
    if cache is not None and key in cache:
        return cache[key]
    f = facts[tuple(run.subshape)]
    O = realise(run.candidate.pattern, f, run.subshape)
    if cache is not None:
        cache[key] = O
    return O


def _not_tried(row, reason):
    row = dict(row)
    row["tried"] = False
    row["not_tried_reason"] = reason
    return row


def _only_para(problem):
    """Whether every mu-group of ``problem`` is empty or full."""
    blocks = problem.blocks
    sizes = [0] * len(blocks.group_nconds)
    for blk, g in zip(blocks.block_info, blocks.block_to_group):
        sizes[int(g)] += len(blk) * int(problem.nvol)
    for n, size in zip(blocks.group_nconds, sizes):
        n = float(n)
        if 0 < n < size:
            return False
    return True


def _response_checks(runs, rows, facts, *, nc_zero, only_para, cache):
    """The response checks, before any seeding: with ``nc_zero`` or
    ``only_para`` every non-para run is a no-response row; otherwise a
    pattern outside the solver blocks is refused and a constraint mode is
    a no-response row."""
    by_id = {r["run_id"]: i for i, r in enumerate(rows)}
    rows = list(rows)
    kept = []
    for run in runs:
        fam = run.candidate.family
        reason = None
        if fam != "para":
            if nc_zero:
                reason = _NO_RESPONSE
            elif only_para:
                reason = DEGENERATE_REASON
            else:
                f = facts[tuple(run.subshape)]
                p = f.problem
                try:
                    O = _realised(run, facts, cache)
                except ValueError as exc:
                    if fam not in EXEMPT_FAMILIES:
                        raise
                    reason = "refused: pattern: {}".format(exc)
                else:
                    msgs = block_check(O, p.blocks, ns=p.ns, norb=p.norb)
                    if msgs:
                        reason = "{}: {}".format(BLOCK_REASON, "; ".join(msgs))
                    elif constraint_mode(O, f):
                        reason = CONSTRAINT_REASON
        if reason is None:
            kept.append(run)
        else:
            i = by_id[run.run_id]
            rows[i] = _not_tried(rows[i], reason)
    return kept, rows


#: two pattern invariants closer than this (relative) may be related by a
#: map; farther apart the pair is not compared (a pre-filter only)
_INVARIANT_RTOL = 1e-6


def _same_invariants(a, b):
    for x, y in zip(a, b):
        if x.shape != y.shape:
            return False
        scale = max(1.0, float(np.max(np.abs(x), initial=0.0)),
                    float(np.max(np.abs(y), initial=0.0)))
        if float(np.max(np.abs(x - y), initial=0.0)) > _INVARIANT_RTOL * scale:
            return False
    return True


def _candidates_by_invariants(kept, inv, i):
    """The kept runs whose invariants match those of run ``i``
    (:func:`_same_invariants`, vectorised over ``kept``), in order."""
    if not kept:
        return []
    ok = np.ones(len(kept), dtype=bool)
    for a, x in enumerate(inv[i]):
        ys = [inv[j][a] for j in kept]
        same = np.array([y.shape == x.shape for y in ys])
        ok &= same
        if not ok.any():
            return []
        sel = np.flatnonzero(ok)
        Y = np.stack([ys[k] for k in sel])
        scale = np.maximum(np.maximum(1.0, np.max(np.abs(x), initial=0.0)),
                           np.max(np.abs(Y), axis=1, initial=0.0))
        diff = np.max(np.abs(Y - x[None, :]), axis=1, initial=0.0)
        ok[sel] = diff <= _INVARIANT_RTOL * scale
    return [kept[k] for k in np.flatnonzero(ok)]


def deduplicate(runs, *, facts, su2, tr, cache=None):
    """Merge runs that a verified map takes into one another, per evaluated
    ``SubShape``.

    ``facts``: a mapping ``SubShape -> ModelFacts`` loaded on that cell
    (the map's problem: a translation is a translation of that cell).
    ``su2``, ``tr``: the checks of :mod:`hwave.initial.symmetry`. Within one
    B the runs are taken in the given order with the reserved runs (para,
    ferro, random) first; each is compared with the runs kept so far by
    :func:`symmetry.find_map`. A match is a duplicate row,
    ``"duplicate of <run_id> (<map>)"`` naming the kept run on the same B,
    and the kept run's row ``sources`` gain the duplicate's sources and
    ``"merged <run_id> (<map>)"`` (the frozen ``Candidate.sources`` is not
    changed). Reserved runs are never removed, and neither are user and
    literature runs: such a run that
    matches a kept one stays a run, and its row is only marked,
    ``duplicate_of = "<run_id> (<map>)"`` (``None`` otherwise), the kept
    row's ``sources`` gaining ``"equivalent <run_id> (<map>)"``. ``cache``:
    an optional dict ``run_id -> O(k)`` shared with other steps.

    Returns ``(kept runs, rows)``: one row per input run, in input order
    (rows of kept runs ``tried = True``)."""
    runs = list(runs)
    order = {}
    for i, run in enumerate(runs):
        order.setdefault(tuple(run.subshape), []).append(i)
    rows = [None] * len(runs)
    removed = set()
    for B, idx in order.items():
        f = facts[B]
        idx = sorted(idx, key=lambda i: (runs[i].candidate.family not in RESERVED_FAMILIES, i))
        kept = []
        inv = {}
        # O(R) of each run, transformed once: the invariants and every map
        # search of the run read the same array
        o_r = {}
        for i in idx:
            run = runs[i]
            O = _realised(run, facts, cache)
            o_r[i] = symmetry._to_r(O, f.problem)
            inv[i] = symmetry._map_invariants_r(O, o_r[i], f.problem, su2)
            rows[i] = _row(run, f.cellshape)
            target = None
            for j in _candidates_by_invariants(kept, inv, i):
                m = symmetry._find_map_r(o_r[j], o_r[i], f.problem, su2, tr)
                if m is not None:
                    target = (j, m)
                    break
            if target is None:
                kept.append(i)
                continue
            j, m = target
            if run.candidate.family in RESERVED_FAMILIES + EXEMPT_FAMILIES:
                rows[i]["duplicate_of"] = "{} ({})".format(runs[j].run_id, m)
                rows[j]["sources"].append("equivalent {} ({})".format(run.run_id, m))
                kept.append(i)
                continue
            rows[i] = _not_tried(rows[i], "duplicate of {} ({})".format(runs[j].run_id, m))
            src = rows[j]["sources"]
            for s in run.candidate.sources:
                if s not in src:
                    src.append(s)
            src.append("merged {} ({})".format(run.run_id, m))
            removed.add(i)
    return [r for i, r in enumerate(runs) if i not in removed], rows


def apply_quotas(runs, rows, limits):
    """The per-``SubShape`` quotas.

    Only generated runs are counted: baseline and orbital together (the
    baseline quota) and instability. Per B, instability is first capped at
    ``quota_instability``; then the baseline rows keep up to
    ``max_candidates`` and instability the rest of ``max_candidates``
    (instability is cut first). Reserved, user and literature runs are
    never cut. A cut run leaves ``runs`` and its row reads ``"quota:
    <family>"`` (``tried = False``); rows of not-tried pairs consume no
    quota. Runs keep their order; the later rows of a family are cut.

    A duplicate row (``"duplicate of <X> (<map>)"``) whose referent X is cut
    is promoted to X's place in the quota order and carries the note
    ``"referent not run"``. A duplicate always follows its referent;
    reserved, user and literature runs are never cut, so a cut referent is
    generated, and the duplicate is generated too: of the same family, or
    an instability row duplicating a baseline one (when a baseline row is
    cut the baseline quota is full and instability gets nothing). X's place
    therefore lies past the quota and the promoted row is cut with it. It
    reads ``"quota: <its family> (duplicate of <X>, a <X's family> row cut
    by its quota; referent not run)"``, so no row claims a basin that was
    not run.

    Returns ``(runs, rows)``; ``rows`` is a new list, one entry per input
    row."""
    by_B = {}
    for i, run in enumerate(runs):
        by_B.setdefault(tuple(run.subshape), []).append(i)
    cut = set()
    for B, idx in by_B.items():
        base = [i for i in idx if runs[i].candidate.family in ("baseline", "orbital")]
        inst = [i for i in idx if runs[i].candidate.family == "instability"]
        n_base = min(len(base), limits.max_candidates)
        n_inst = min(len(inst), limits.quota_instability, limits.max_candidates - n_base)
        cut.update(base[n_base:])
        cut.update(inst[n_inst:])
    reasons = {runs[i].run_id: "quota: {}".format(runs[i].candidate.family) for i in cut}
    families = {runs[i].run_id: runs[i].candidate.family for i in cut}
    out = []
    for r in rows:
        if r.get("run_id") in reasons:
            out.append(_not_tried(r, reasons[r["run_id"]]))
            continue
        why = r.get("not_tried_reason") or ""
        if why.startswith("duplicate of "):
            target = why[len("duplicate of "):].split(" ", 1)[0]
            if target in reasons:
                out.append(_not_tried(r, "quota: {} (duplicate of {}, a {} row cut by "
                                         "its quota; referent not run)".format(
                                             r.get("family"), target, families[target])))
                continue
        out.append(dict(r))
    return [r for i, r in enumerate(runs) if i not in cut], out


# --------------------------------------------------------------------------
# the instability rows
# --------------------------------------------------------------------------

#: the reason of an instability row that its SubShape could run but the
#: instability quota of that SubShape leaves out (it is not realised)
QUOTA_INSTABILITY = "quota: instability"


@dataclass(frozen=True, eq=False)
class HintPair:
    """One eigenpair of the hints, as the instability rows use it.

    ``q`` (reduced Fractions), ``index`` (its place among the pairs of its
    q, merged over the sectors: the ``r<rank>`` of the candidate id),
    ``rank`` (its place among every pair of every q: by descending
    degenerate set of ``lam`` (:func:`hwave.initial.hints.degenerate_clusters`),
    inside a set by q order, increasing order n then lexicographic, then
    ``index``; recorded as the row's ``rank``), ``lam`` and ``lam_imag``
    (the eigenvalue), ``sector``, ``group`` (its degenerate set within
    q), ``amplitude`` (the conjugate field ``O = Gamma x`` as PatternSpec
    entries) and ``general_solve`` (the vertex was not self-adjoint at
    this q, so the pair comes from the general solve and its rows are
    marked)."""
    q: Tuple[Fraction, Fraction, Fraction]
    index: int
    rank: int
    lam: float
    lam_imag: float
    sector: str
    group: int
    amplitude: Tuple[Entry, ...]
    general_solve: bool


def hint_pairs(eigenpairs) -> Tuple[HintPair, ...]:
    """The eigenpairs of :func:`hwave.initial.hints.hint_eigenpairs` as
    :class:`HintPair`, in rank order (the order of the instability
    rows)."""
    flat = []
    for q, res in eigenpairs.items():
        qf = cells.to_qfrac(q)
        for i, ((lam, _x), lim, sec, grp, amp) in enumerate(zip(
                res.pairs, res.lambda_imag, res.sectors, res.groups, res.amplitudes)):
            flat.append((qf, i, float(lam), float(lim), sec, int(grp), tuple(amp),
                         not res.gamma_self_adjoint))
    # degenerate sets across q first (lambda(q) = lambda(-q) and lattice
    # images are equal by symmetry), then q order and the index at q, so
    # round-off never decides the rank
    cluster = hints_mod.degenerate_clusters([t[2] for t in flat])
    key = {id(t): c for t, c in zip(flat, cluster)}
    flat.sort(key=lambda t: (key[id(t)], cells.order(t[0]), t[0], t[1]))
    return tuple(HintPair(q=t[0], index=t[1], rank=k, lam=t[2], lam_imag=t[3],
                          sector=t[4], group=t[5], amplitude=t[6], general_solve=t[7])
                 for k, t in enumerate(flat))


@dataclass(frozen=True)
class InstabilityCandidates:
    """The instability rows of :func:`instability_candidates`.

    ``runnable``: the realised candidates, per ``SubShape`` in the given
    order, in rank order of their pairs (each pair's phase rows, the -O
    row of each at q = -q, and the spiral). ``listed``: ``(candidate,
    subshape, reason)`` of every pair that is not realised, named from its
    amplitude: ``"quota: instability"`` on a ``SubShape`` that holds its q
    (the quota of that cell leaves it out), or ``"unrepresentable"`` (no
    ``SubShape`` holds its q; listed once, on the first ``SubShape``, with
    the PatternSpec ``subshape`` its minimal cell). ``hint_of``: the
    :class:`HintPair` of every candidate id. ``realised``: per
    ``SubShape``, the number of pairs realised (or tried) there; ``held``:
    per ``SubShape``, the number of pairs whose q it holds."""
    runnable: Tuple[Candidate, ...]
    listed: Tuple[Tuple[Candidate, Tuple[int, int, int], str], ...]
    hint_of: Mapping
    # default_factory: a mapping proxy is unhashable before Python 3.12, and
    # dataclasses refuse an unhashable default there
    realised: Mapping = dataclasses.field(default_factory=lambda: MappingProxyType({}))
    held: Mapping = dataclasses.field(default_factory=lambda: MappingProxyType({}))


def _site_of(facts):
    out = {}
    for k, group in enumerate(facts.site_groups.groups):
        for a in group:
            out[int(a)] = k
    return out


def _amplitude_name(A, facts, q, phase):
    """The name of a listed (not realised) instability row from its
    amplitude: the dominant channels (weight at least
    :data:`hwave.initial.labels.DOMINANT_WEIGHT`) among ``charge``,
    ``spin-<axis>`` (on-site, orbital-diagonal), ``orbital-hybridization``,
    ``inter-orbital-spin-<axis>`` (on-site, two orbitals) and
    ``bond-charge``, ``bond-spin-<axis>`` (two sites), then q and the
    phase."""
    site = _site_of(facts)
    weights = {}
    for e in A:
        same = e.r == (0, 0, 0) and site.get(e.orb_a) == site.get(e.orb_b)
        for key, mu in (("c0", None), ("cx", "x"), ("cy", "y"), ("cz", "z")):
            w = abs(getattr(e, key)) ** 2
            if w == 0:
                continue
            if not same:
                name = "bond-charge" if mu is None else "bond-spin-" + mu
            elif e.orb_a == e.orb_b:
                name = "charge" if mu is None else "spin-" + mu
            else:
                name = ("orbital-hybridization" if mu is None
                        else "inter-orbital-spin-" + mu)
            weights[name] = weights.get(name, 0.0) + w
    total = sum(weights.values())
    ranked = sorted(weights, key=lambda k: (-weights[k], k))
    names = [k for k in ranked if weights[k] >= labels.DOMINANT_WEIGHT * total] or ranked[:1]
    return " ".join([" + ".join(names)] + labels._name_suffix(q, None, phase, None, False))


def _instability_rows(facts, B, h, *, su2, tr, limits, taken, patterns):
    """The realised candidates of one pair on the cell B."""
    p = facts.problem
    q = h.q
    A = h.amplitude
    sources = ("instability",) + ((hints_mod.GENERAL_SOLVE_SOURCE,)
                                  if h.general_solve else ())
    if cells.is_self_conjugate(q):
        rows = pattern_rows(A, q, sign_map=True, subshape=B)
        twins = True
    else:
        O0 = realise(pattern_rows(A, q, sign_map=True, subshape=B)[0].spec, facts, B)
        smap = symmetry.sign_map_verified(O0, problem=p, su2=su2, tr=tr)
        spiral = (h.sector == hints_mod.SECTOR_SZ and su2 is not None
                  and su2.status == "passed")
        rows = pattern_rows(A, q, sign_map=smap, spiral=spiral, subshape=B)
        twins = False
    out, k = [], 0
    for row in rows:
        index = 0 if row.spiral else k
        cand, _ = _named("instability", row.spec, facts, B, eta=None, spiral=row.spiral,
                         rank=h.index, phase_index=index, neg=False, taken=taken,
                         amplitude=limits.amplitude, sources=sources, patterns=patterns)
        out.append(cand)
        if twins and not row.spiral:
            twin, _ = _named("instability", _neg(row.spec), facts, B, eta=None,
                             spiral=False, rank=h.index, phase_index=index, neg=True,
                             taken=taken, amplitude=limits.amplitude, sources=sources,
                             patterns=patterns)
            out.append(twin)
        if not row.spiral:
            k += 1
    return out


def _listed(facts, h, subshape, taken, limits):
    A = fix_sign(h.amplitude, h.q)
    phase = phi0(A, h.q)
    spec = PatternSpec(q=tuple(str(x) for x in h.q), subshape=tuple(subshape), phase=phase,
                       sign=1, modulation="cos", entries=A)
    name = _amplitude_name(A, facts, h.q, phase)
    sources = ("instability",) + ((hints_mod.GENERAL_SOLVE_SOURCE,)
                                  if h.general_solve else ())
    cid = make_candidate_id("instability", spec.q, labels.ascii_name(name), h.index, 0,
                            False, taken)
    return Candidate(candidate_id=cid, family="instability", sources=sources, pattern=spec,
                     amplitude_requested=limits.amplitude, label=name)


def instability_candidates(facts, eigenpairs, subshapes, limits, *, su2, tr, taken=None,
                           patterns=None, budget=None, memo=None) -> InstabilityCandidates:
    """The instability rows of the hint eigenpairs on the ``SubShape``
    list.

    ``facts``: a mapping ``SubShape -> ModelFacts`` loaded on that cell;
    ``eigenpairs``: :func:`hwave.initial.hints.hint_eigenpairs` (or the
    :class:`HintPair` tuple of :func:`hint_pairs`); ``limits``:
    :class:`Limits`; ``su2``, ``tr``: the checks of
    :mod:`hwave.initial.symmetry`.

    Per ``SubShape`` B, the pairs whose q B holds are taken in rank order
    (best lambda first). The first of them are realised: the conjugate
    field ``O = Gamma x`` of the pair becomes the anchored modulation of
    :mod:`hwave.initial.patterns` (sign fixed once on the phi0 row,
    :func:`~hwave.initial.patterns.fix_sign`; the phase rows of
    :func:`~hwave.initial.patterns.phase_rows`, with the sign map of
    :func:`hwave.initial.symmetry.sign_map_verified` at q not equal to -q;
    at q = -q each phase row followed by its -O row, which deduplication
    merges when a verified map exists; a spiral for an ``S^z`` pair at q
    not equal to -q under a passed spin-rotation check), named by its
    pattern (family ``instability``, id rank = the pair's index at its q).
    ``budget`` maps B to the number of pairs to realise; by default pairs
    are realised until their rows reach ``limits.quota_instability``,
    the rest of a degenerate set included. Every other pair is listed,
    not realised (:class:`InstabilityCandidates`); a pair no B holds is
    listed as unrepresentable. Eigenvalues only order the rows: no row is
    removed for its eigenvalue. ``memo``: an optional dict kept between
    calls with a larger ``budget`` and the same ``taken`` start: per
    ``SubShape`` the ids taken at the start of that cell and its realised
    candidates by pair rank, so a pair is realised once; the memo of a
    cell is dropped when the ids taken before it differ (an earlier cell
    realised more), since its ids would no longer be unique."""
    pairs = eigenpairs if isinstance(eigenpairs, tuple) else hint_pairs(eigenpairs)
    subshapes = [tuple(int(x) for x in B) for B in subshapes]
    taken = set() if taken is None else taken
    runnable, listed, hint_of = [], [], {}
    held, realised, n_held = set(), {}, {}
    for B in subshapes:
        f = facts[B]
        mine = [h for h in pairs if cells.representable(h.q, B)[0]]
        held.update(id(h) for h in mine)
        n_held[B] = len(mine)
        n = None if budget is None or B not in budget else int(budget[B])
        count = 0
        last = None
        realised[B] = 0
        # the memo of this cell holds the ids taken from the set as it was at
        # the start of the cell; when another cell has taken other ids since
        # (its budget grew), every realisation of this cell is redone
        if memo is not None:
            start = frozenset(taken)
            mine_memo = memo.get(B)
            if mine_memo is None or mine_memo[0] != start:
                mine_memo = (start, {})
                memo[B] = mine_memo
            cell_memo = mine_memo[1]
        else:
            cell_memo = None
        for i, h in enumerate(mine):
            take = (count < limits.quota_instability if n is None else i < n)
            if not take and last is not None and (h.q, h.group) == last:
                take = True                   # the rest of a degenerate set
            reason = QUOTA_INSTABILITY
            if take:
                realised[B] += 1
                key = h.rank
                if cell_memo is not None and key in cell_memo:
                    rows = cell_memo[key]
                    if rows is not None:
                        taken.update(c.candidate_id for c in rows)
                else:
                    try:
                        rows = tuple(_instability_rows(f, B, h, su2=su2, tr=tr,
                                                       limits=limits, taken=taken,
                                                       patterns=patterns))
                    except ValueError:
                        # the field vanishes on this cell (realise refuses a
                        # zero pattern): listed as giving no response
                        rows = None
                    if cell_memo is not None:
                        cell_memo[key] = rows
                if rows is None:
                    reason = _NO_RESPONSE
                else:
                    for c in rows:
                        hint_of[c.candidate_id] = h
                    runnable.extend(rows)
                    count += len(rows)
                    last = (h.q, h.group)
                    continue
            c = _listed(f, h, B, taken, limits)
            hint_of[c.candidate_id] = h
            listed.append((c, B, reason))
    first = facts[subshapes[0]]
    for h in pairs:
        if id(h) in held:
            continue
        c = _listed(first, h, tuple(x.denominator for x in h.q), taken, limits)
        hint_of[c.candidate_id] = h
        listed.append((c, subshapes[0], _UNREPRESENTABLE))
    return InstabilityCandidates(runnable=tuple(runnable), listed=tuple(listed),
                                 hint_of=MappingProxyType(hint_of),
                                 realised=MappingProxyType(realised),
                                 held=MappingProxyType(n_held))


# --------------------------------------------------------------------------
# the entry point
# --------------------------------------------------------------------------

#: the interim size limits of the small tier (deterministic: the same input
#: gets the same tier on every host): the pair space, the q points (the
#: k-mesh of the primitive cell) and ``N_q |P|^3`` (the dense eigen-solves)
SMALL_TIER_MAX_PAIRS = 400
SMALL_TIER_MAX_Q = 4096
SMALL_TIER_MAX_OPERATIONS = 10 ** 11
#: the warning when the interim rule does not choose the small tier
SMALL_TIER_LIMIT_WARNING = ("the small tier exceeds the preview limits; the medium tier "
                            "arrives in a later release")
#: the warning when a small tier given explicitly is over the size limits
SMALL_TIER_OVER_BUDGET_WARNING = (
    "the small tier given explicitly exceeds the preview limits ({}); it runs anyway")
_TIERS_IN = ("auto", "large", "medium", "small")


def small_tier_bytes(n_pair, n_q, n_k, n_bands, top_k) -> int:
    """The small-tier memory estimate ``16 (|P|^2 + N_q top_k |P| + N_k nb^2
    + N_k nb |P|)`` bytes."""
    P, nq, nk, nb, k = (int(v) for v in (n_pair, n_q, n_k, n_bands, top_k))
    return 16 * (P * P + nq * k * P + nk * nb * nb + nk * nb * P)


@dataclass(frozen=True)
class TierChoice:
    """The tier of an enumeration: ``name`` (``"large"`` or ``"small"``),
    ``chosen_by`` (``"auto (interim rule)"`` or ``"user"``), ``record`` (the
    estimate, budget, constants, rule and preflight of the ``tier`` record
    of ``coverage.json``) and ``warnings``."""
    name: str
    chosen_by: str
    record: Mapping
    warnings: Tuple[str, ...] = ()


def choose_tier(base, tier="auto", *, su2=None, limits=None) -> TierChoice:
    """The tier of the model of ``base`` (primitive-cell facts), before any
    para step.

    The interim rule of this release for ``"auto"`` depends on the model
    only, so an input gets the same tier on every host: the small tier when
    ``|P| <=`` :data:`SMALL_TIER_MAX_PAIRS` (400), ``N_q <=``
    :data:`SMALL_TIER_MAX_Q` (4096) and ``N_q |P|^3 <=``
    :data:`SMALL_TIER_MAX_OPERATIONS` (1e11), with ``|P|`` from
    :func:`hwave.initial.hints.pair_space` and ``N_q = prod(CellShape)``;
    otherwise the large tier with the warning
    :data:`SMALL_TIER_LIMIT_WARNING`. ``"small"`` given explicitly runs
    the same size check (a warning when over a limit). The host memory
    (``initial_facade.available_host_bytes()``) never chooses: it only
    refuses a small tier whose memory estimate (:func:`small_tier_bytes`)
    exceeds the available memory (``ValueError``); unknown memory refuses
    nothing. ``"large"`` needs no estimate; ``"medium"`` raises
    ``NotImplementedError`` (a later release)."""
    from hwave.initial.report import BUDGET_FRACTION, INTERIM_TIER_RULE
    from hwave.solver import initial_facade
    if tier not in _TIERS_IN:
        raise ValueError("tier: must be 'auto', 'large', 'medium' or 'small', got {!r}"
                         .format(tier))
    if tier == "medium":
        raise NotImplementedError("tier 'medium' arrives in a later release; this release "
                                  "has the large and small tiers")
    limits = Limits() if limits is None else limits
    chosen_by = "auto (interim rule)" if tier == "auto" else "user"
    record = {
        "estimate": {"bytes": None, "seconds": None, "class": None},
        "budget": {"budget_fraction": BUDGET_FRACTION, "bytes": None,
                   "available_bytes": None, "source": None},
        "constants": {"small_max_pairs": SMALL_TIER_MAX_PAIRS,
                      "small_max_q": SMALL_TIER_MAX_Q,
                      "small_max_operations": SMALL_TIER_MAX_OPERATIONS},
        "rule": INTERIM_TIER_RULE if tier == "auto" else None,
        "preflight": None,
    }
    if tier == "large":
        return TierChoice("large", chosen_by, record)
    p = base.problem
    n_pair = hints_mod.pair_space(base, su2=su2).size
    nk = int(np.prod(base.cellshape))
    nb = int(p.nd)
    nbytes = small_tier_bytes(n_pair, nk, nk, nb, limits.top_k)
    ops = nk * n_pair ** 3
    over = []
    if n_pair > SMALL_TIER_MAX_PAIRS:
        over.append("|P| = {} > {}".format(n_pair, SMALL_TIER_MAX_PAIRS))
    if nk > SMALL_TIER_MAX_Q:
        over.append("N_q = {} > {}".format(nk, SMALL_TIER_MAX_Q))
    if ops > SMALL_TIER_MAX_OPERATIONS:
        over.append("N_q |P|^3 = {} > {:.0e}".format(ops, SMALL_TIER_MAX_OPERATIONS))
    record["preflight"] = {"pair_space": n_pair, "n_q": nk, "n_k": nk, "bands": nb,
                           "top_k": limits.top_k, "small_bytes": nbytes,
                           "small_operations": ops, "over_limits": over}
    warnings = []
    if tier == "auto" and over:
        return TierChoice("large", chosen_by, record, (SMALL_TIER_LIMIT_WARNING,))
    if tier == "small" and over:
        warnings = [SMALL_TIER_OVER_BUDGET_WARNING.format("; ".join(over))]
    mem = initial_facade.available_host_bytes()
    record["budget"].update(
        bytes=None if mem.nbytes is None else int(BUDGET_FRACTION * mem.nbytes),
        available_bytes=mem.nbytes, source=mem.source)
    if mem.nbytes is not None and nbytes > mem.nbytes:
        raise ValueError("tier small{}: the memory estimate {} bytes exceeds the available "
                         "host memory {} bytes ({}); give the large tier".format(
                             " (chosen by the interim rule)" if tier == "auto" else "",
                             nbytes, mem.nbytes, mem.source))
    record["estimate"]["bytes"] = nbytes
    return TierChoice("small", chosen_by, record, tuple(warnings))


def _input_path(obj):
    if isinstance(obj, ModelFacts):
        raise TypeError("candidate_states needs the input path (a TOML or a directory): "
                        "it loads the model once per SubShape")
    if not isinstance(obj, (str, os.PathLike)):
        raise TypeError("input_dir_or_model: a TOML file or an input directory, got {!r}"
                        .format(obj))
    return obj


def _subshape_list(subshapes):
    if isinstance(subshapes, (str, bytes)) or not hasattr(subshapes, "__iter__"):
        raise ValueError("subshapes: a list of SubShape triples, got {!r}".format(subshapes))
    out = []
    for raw in subshapes:
        if isinstance(raw, (str, bytes)) or not hasattr(raw, "__len__") or len(raw) != 3:
            raise ValueError("subshapes: each SubShape is three positive integers, got "
                             "{!r}".format(raw))
        B = tuple(_check_int("SubShape", x, 1) for x in raw)
        if B in out:
            raise ValueError("subshapes: SubShape {} is given twice".format(list(B)))
        out.append(B)
    return out


def _extras(extra):
    out = []
    for i, c in enumerate(extra or ()):
        if not isinstance(c, Candidate) or c.family not in EXEMPT_FAMILIES:
            raise ValueError("extra[{}]: a user or literature Candidate is required "
                             "(Candidate.user / Candidate.literature), got {!r}".format(i, c))
        out.append(c)
    seen = {}
    for i, c in enumerate(out):
        if c.candidate_id in seen:
            raise ValueError(
                "extra[{}]: candidate_id {!r} is also the id of extra[{}]; build the "
                "candidates with one shared taken= set (Candidate.user(..., taken=ids)) so "
                "that a collision gets a -2 suffix".format(i, c.candidate_id,
                                                           seen[c.candidate_id]))
        seen[c.candidate_id] = i
    return out


def _json_num(v):
    if v is None or isinstance(v, (bool, str)):
        return v
    if isinstance(v, Fraction):
        return str(v)
    if isinstance(v, numbers.Integral):
        return int(v)
    return float(v)


#: the statuses of a spin-rotation check that ran
_SU2_APPLIED = ("passed", "failed")


def _model_record(base, su2, tr, hint_rec, warnings, tier_record):
    sg = base.site_groups
    return {
        "spin_orbital": bool(base.spin_orbital), "norb": int(base.norb),
        "norb_phys": int(base.norb_phys), "iflag_fock": bool(base.iflag_fock),
        "interacting_orbitals": [int(a) for a in base.interacting_orbitals],
        "n_c": int(base.n_c),
        "site_groups": {"source": sg.source, "groups": [list(g) for g in sg.groups],
                        "ambiguous": list(sg.ambiguous)},
        "Ncond": _json_num(base.ncond), "Ncond_per_cell": _json_num(base.ncond_per_cell),
        "T": float(base.T), "2Sz": base.two_sz, "cellshape": list(base.cellshape),
        # the report's su2 record: max_residual and method are null when the
        # check is not applicable (the S_z-only check of 2Sz mode records its
        # result in sz_residual)
        "su2": {"status": su2.status,
                "max_residual": su2.max_residual if su2.status in _SU2_APPLIED else None,
                "method": su2.method if su2.status in _SU2_APPLIED else None,
                "sz_residual": su2.sz_residual},
        "tr": {"invariant": bool(tr.invariant), "h0_residual": float(tr.h0_residual),
               "hf_residual": float(tr.hf_residual)},
        "warnings": list(warnings),
        "input_sha256": dict(base.input_sha256),
        "hint": hint_rec,
        "tier_record": tier_record,
    }


#: the warning when a para check is asked for but the tier has no para step
PARA_CHECK_LARGE_WARNING = ("the para check is not run: the large tier has no "
                            "paramagnetic Hartree-Fock step")


def _stage(cands, subshapes, base, facts, lim, cache):
    """Runs and rows of ``cands`` after :func:`expand_runs` and the response
    checks."""
    runs, rows = expand_runs(cands, subshapes, cellshape=base.cellshape,
                             max_denominator=lim.max_denominator)
    return _response_checks(runs, rows, facts, nc_zero=base.n_c == 0,
                            only_para=_only_para(base.problem), cache=cache)


def _merge_by_subshape(first, second, subshapes):
    """Runs (or rows) of two stages in the order of :func:`expand_runs`:
    SubShape by SubShape, the first stage before the second."""
    def B_of(x):
        return tuple(x.subshape) if isinstance(x, CandidateRun) else tuple(x["subshape"])
    out = []
    for B in subshapes:
        out.extend(x for x in first if B_of(x) == B)
        out.extend(x for x in second if B_of(x) == B)
    return out


def _finish_stage(runs, rows, facts, lim, cache, *, su2, tr):
    """:func:`deduplicate` and :func:`apply_quotas` of staged runs."""
    kept, dup_rows = deduplicate(runs, facts=facts, su2=su2, tr=tr, cache=cache)
    index = {r["run_id"]: i for i, r in enumerate(rows)}
    rows = list(rows)
    for r in dup_rows:
        rows[index[r["run_id"]]] = r
    return apply_quotas(kept, rows, lim)


def _saturated(rows, B):
    """Whether the instability quota of B cut a realised row."""
    return any(r["family"] == "instability" and tuple(r["subshape"]) == B
               and (r.get("not_tried_reason") or "").startswith("quota: instability")
               for r in rows)


def _listed_row(cand, B, reason, cellshape, lim):
    run = CandidateRun(cand, B)
    if reason == _UNREPRESENTABLE:
        ok, why = cells.representable(cand.pattern.q, B)
        row = _row(run, cellshape, representable=False, reason=why)
        sub, cell, note = _proposal(cand.pattern, cellshape, lim.max_denominator)
        row.update(proposed_subshape=sub, proposed_cellshape=cell, proposal_note=note)
        return row
    return _not_tried(_row(run, cellshape), reason)


def candidate_states(input_dir_or_model, *, subshapes=None, tier="auto", t_hint=None,
                     extra=(), orbital_site=None, para=None, para_check=False,
                     para_dir=None, **limits) -> CandidateSet:
    """The candidate set of a UHFk input.

    Parameters
    ----------
    input_dir_or_model : path-like
        A UHFk TOML, or a directory holding exactly one (the ``Model`` of
        the input builder is not supported yet).
    subshapes : list of three-int sequences, optional
        The candidate ``SubShape`` list; default the input's ``SubShape``.
        Without either, ``ValueError`` with the warning
        :data:`hwave.initial.facts.NO_SUBSHAPE_WARNING`. Each must divide
        ``CellShape``; a duplicate is refused.
    tier : ``"auto"``, ``"large"``, ``"medium"`` or ``"small"``
        Chosen first, before any para step (:func:`choose_tier`; ``"auto"``
        follows the interim rule there). ``"medium"`` raises
        ``NotImplementedError`` (a later release).
    t_hint : float, optional
        The user T_hint (the hints only; never the other rows).
    extra : sequence of Candidate
        User and literature candidates (evaluated on every SubShape).
    orbital_site : list of int, optional
        A site label per physical orbital (overrides the input and the
        clustering).
    para : ParaResult, optional
        The small tier: a para reference computed before
        (:func:`hwave.initial.hints.para_hf` or
        :func:`hwave.initial.hints.para_from_work`) for this input and
        T_hint; ``None`` runs the para step here, after the tier is chosen.
    para_check : bool
        The small tier: the para step also runs from the uniform orbital
        occupation (:func:`hwave.initial.hints.para_hf`). With the large
        tier given explicitly ``ValueError``; when ``"auto"`` chose the
        large tier, a warning.
    para_dir : path-like, optional
        The work directory of the para step (its runs, ``green.npz`` and
        ``para_hf.json`` are kept there).
    **limits
        The fields of :class:`Limits`; any other keyword raises
        ``TypeError``.

    In the small tier the hints run on the converged para reference
    (:func:`hwave.initial.hints.hint_eigenpairs` over the whole k-mesh of
    the primitive cell); their pairs give the instability rows
    (:func:`instability_candidates`): the pairs a ``SubShape`` holds are
    realised in rank order until its instability quota is reached (after
    deduplication), every other pair is listed with its eigenvalue and the
    reason it is not tried. The hints change only the instability rows;
    every other row is the same at any T_hint. A para reference that is
    not converged skips the hints with its reason.
    """
    lim = Limits.from_keywords(**limits)
    if tier not in _TIERS_IN:
        raise ValueError("tier: must be 'auto', 'large', 'medium' or 'small', got {!r}"
                         .format(tier))
    if tier == "medium":
        choose_tier(None, tier)                 # NotImplementedError
    if para_check and tier == "large":
        raise ValueError("para_check: the para step runs only in the small tier (the large "
                         "tier has no hints)")
    if para is not None and tier == "large":
        raise ValueError("para: the large tier has no hints, so a para reference is not "
                         "used; give tier small")
    path = _input_path(input_dir_or_model)
    extras = _extras(extra)
    base = model_facts(path, subshape=_UNIT, orbital_site=orbital_site)
    if subshapes is None:
        if NO_SUBSHAPE_WARNING in base.warnings:
            raise ValueError(NO_SUBSHAPE_WARNING)
        subshapes = [model_facts(path, orbital_site=orbital_site).problem.subshape]
    subshapes = _subshape_list(subshapes)
    cells.check_subshapes(subshapes, base.cellshape)

    su2 = symmetry.su2_check(base.problem)
    tr = symmetry.tr_check(base.problem)
    choice = choose_tier(base, tier, su2=su2, limits=lim)
    tier_name = choice.name
    facts = {B: (base if B == _UNIT else model_facts(path, subshape=B,
                                                     orbital_site=orbital_site))
             for B in subshapes}
    hint = t_hint_bare(base.problem, t_user=t_hint)
    warnings_extra = list(choice.warnings)

    pairs = ()
    if tier_name == "small":
        if para is None:
            para = hints_mod.para_hf(base, hint, para_check=para_check, work_dir=para_dir)
        else:
            from hwave.initial.writer import _check_para
            _check_para(para, base, t_hint)
        warnings_extra.extend(para.warnings)
        if para.status == "converged":
            P = hints_mod.pair_space(base, su2=su2)
            eigen = hints_mod.hint_eigenpairs(base, para, P,
                                              cells.commensurate_q(base.cellshape),
                                              top_k=lim.top_k)
            hint_rec = hint_record(hint, para, vertex=hints_mod.vertex_fields(eigen))
            pairs = hint_pairs(eigen)
            for res in eigen.values():
                for w in res.warnings:
                    if w not in warnings_extra:
                        warnings_extra.append(w)
        else:
            hint_rec = hint_record(hint, para)
    else:
        if para_check:
            warnings_extra.append(PARA_CHECK_LARGE_WARNING)
        hint_rec = hint_record(hint, skip_reason=LARGE_TIER_SKIP)

    taken = {c.candidate_id for c in extras}
    reserved = reserved_candidates(base, subshapes, su2=su2, tr=tr, n_random=lim.n_random,
                                   base_seed=lim.base_seed, amplitude=lim.amplitude,
                                   taken=taken)
    cache = {}
    baseline = []
    for B in subshapes:
        baseline.extend(baseline_candidates(facts[B], B, su2=su2, tr=tr, limits=lim,
                                            tier=tier_name, taken=taken, patterns=cache))
    cands = reserved + extras + baseline

    inst = InstabilityCandidates((), (), MappingProxyType({}))
    budget = None
    # the other rows are staged once; a larger budget realises the new pairs
    # only (memo) and repeats deduplication and the quotas
    base_runs, base_rows = _stage(cands, subshapes, base, facts, lim, cache)
    memo = {}
    while True:
        if pairs:
            inst = instability_candidates(facts, pairs, subshapes, lim, su2=su2, tr=tr,
                                          taken=set(taken), patterns=cache, budget=budget,
                                          memo=memo)
        i_runs, i_rows = _stage(list(inst.runnable), subshapes, base, facts, lim, cache)
        runs, rows = _finish_stage(
            _merge_by_subshape(base_runs, i_runs, subshapes),
            _merge_by_subshape(base_rows, i_rows, subshapes), facts, lim, cache,
            su2=su2, tr=tr)
        grow = {}
        for B in subshapes:
            done = inst.realised.get(B, 0)
            if done < inst.held.get(B, 0) and not _saturated(rows, B):
                grow[B] = done + max(1, done)
        if not grow:
            break
        budget = dict(inst.realised)
        budget.update(grow)
    taken.update(c.candidate_id for c in inst.runnable)
    taken.update(c.candidate_id for c, _B, _r in inst.listed)
    rows = list(rows)
    for c, B, reason in inst.listed:
        rows.append(_listed_row(c, B, reason, base.cellshape, lim))
    for r in rows:
        h = inst.hint_of.get(r["candidate_id"])
        if h is not None:
            r.update(rank=h.rank, **{"lambda": h.lam, "lambda_imag": h.lam_imag})

    warnings = list(base.warnings) + list(su2.warnings) + list(hint.warnings)
    warnings.extend(w for w in warnings_extra if w not in warnings)
    if base.n_c == 0:
        warnings.append(NC_ZERO_NOTE)
    elif _only_para(base.problem):
        warnings.append(EMPTY_FULL_NOTE)
    all_cands = cands + list(inst.runnable) + [c for c, _B, _r in inst.listed]
    with_rows = {r["candidate_id"] for r in rows}
    cell_records = [{"subshape": list(B),
                     "commensurate_q": [cells.q_string(q) for q in cells.commensurate_q(B)],
                     "cellshape": list(base.cellshape), "cellshape_compatible": True,
                     "skipped_patterns": skipped_patterns(facts[B], B)}
                    for B in subshapes]
    record = dict(choice.record, chosen_by=choice.chosen_by)
    return CandidateSet(model=_model_record(base, su2, tr, hint_rec, warnings, record),
                        tier=tier_name, cells=cell_records,
                        candidates=[c for c in all_cands if c.candidate_id in with_rows],
                        runs=runs, rows=rows, limits=lim.to_json())
