"""Labels and frames of candidate patterns (design section 6.5).

Each local pattern on a site (or between two sites) is decomposed in the
sigma (x) orbital-pair basis

    O_{(as),(bt)} = 1/2 sum_mu c^mu_ab sigma^mu_st,   mu = 0, x, y, z,

so ``c^mu_ab = tr_s(sigma^mu O_ab)``: mu = 0 is charge/orbital, mu = x, y,
z spin. **The spin axes are the Cartesian frame of the input's spinor
basis** (sigma^z is the spin index of the input files); no lattice-plane
wording is used. For a density matrix ``c^0`` is the density and ``c^i``
twice the spin moment.

* :func:`decompose` / :func:`compose`: the change of basis on a ``(2n,
  2n)`` block, in the interleaved ``(a, s)`` order (``2*a + s``, the
  spin-orbital layout and the literal ``c (x) sigma``) or the spin-block
  ``(s, a)`` order (``s*n + a``, the normal-mode layout).
* :func:`tr_parity`: time reversal in a real orbital gauge maps ``c^0 ->
  conj(c^0)`` and ``c^i -> -conj(c^i)``. TR-even is real (on a Hermitian
  block: real-symmetric) ``c^0`` and imaginary (imaginary-antisymmetric)
  ``c^i``; TR-odd the converse. The parity is ``"even"`` or ``"odd"``
  when the other parity carries less than :data:`TR_MINOR_WEIGHT` of the
  weight, ``"mixed"`` otherwise (owner decision O8).
* :func:`gauge_dependent`: whether the spin-independent part of the
  Transfer table is complex anywhere (then the TR parities, read in the
  orbital gauge of the input, depend on that gauge).
* :func:`pattern_label`: the name of a candidate, read from its seeding
  pattern ``O`` itself (owner decision O8 (a); this name gives the id
  slug). :func:`label`: the name of an induced or converged state, read
  from its ``delta rho`` (O8 (b)). Both name a state by its dominant
  components only (weight at least :data:`DOMINANT_WEIGHT`), on-site
  channels before bond channels.

**Names.** ``Label.name`` is display text and may hold ``φ``, ``π``,
``η`` and a caller's ``eta`` name verbatim (``"Néel"``).
``Label.ascii_name`` is the same name in ASCII (``phi``, ``pi``, ``eta``,
accents folded, other non-ASCII characters dropped); id slugs are built
from it.

**Convention.** :func:`pattern_label` takes the Hamiltonian-convention
``O(k)`` of the folded cell that ``realise`` returns (as
:mod:`hwave.initial.symmetry` does). :func:`label` takes ``delta rho`` in the rho convention,
``rho = G^T`` per displacement (``rho(r)_{ij} = <c^dag_j(r) c_i(0)>`` with
``G[r, s, a, t, b] = <c^dag_{sa}(r) c_{tb}(0)>`` the stored slot); the
callers transpose a stored Green function before calling it. Orbital
indices are physical (ruling R8): in spin-orbital mode the interleaved
``2*orb + spin`` index of the solver is split into the physical orbital
and the spin here, internally.
"""
import math
import numbers
import unicodedata
import weakref
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, Literal, Mapping

import numpy as np

from hwave.initial.patterns import q_fractions
from hwave.initial.symmetry import hermitian_green, _to_r

__all__ = [
    "CHANNELS", "TR_KEYS", "DOMINANT_WEIGHT", "TR_MINOR_WEIGHT", "PARA_TOL", "GAUGE_TOL",
    "HERMITIAN_TOL", "THRESHOLDS", "LAYOUTS", "Label", "decompose", "compose", "tr_parity",
    "gauge_dependent", "label", "pattern_label", "ascii_name",
]

#: the channels a label is read from, in their fixed (tie-break) order.
#: On-site entries (one site, ``r = 0``): ``charge`` (the per-site mean of
#: the diagonal ``c^0``), ``orbital-polarization`` (the rest of the
#: diagonal ``c^0``), ``orbital-hybridization`` / ``orbital-current``
#: (real / imaginary off-diagonal ``c^0``), ``spin-<axis>`` (diagonal
#: ``c^i``), ``inter-orbital-spin-<axis>`` / ``spin-orbital-even`` (real /
#: imaginary off-diagonal ``c^i``, the latter summed over the axes). Bond
#: entries (two sites, or ``r != 0``): ``bond-charge`` / ``bond-current``
#: (real / imaginary ``c^0``), ``bond-spin-<axis>`` /
#: ``bond-spin-orbital-even`` (real / imaginary ``c^i``).
CHANNELS = (
    "charge", "orbital-polarization", "orbital-hybridization", "orbital-current",
    "spin-x", "spin-y", "spin-z",
    "inter-orbital-spin-x", "inter-orbital-spin-y", "inter-orbital-spin-z",
    "spin-orbital-even",
    "bond-charge", "bond-current", "bond-spin-x", "bond-spin-y", "bond-spin-z",
    "bond-spin-orbital-even",
)

#: the TR-even channels (the others are TR-odd)
_EVEN = frozenset(("charge", "orbital-polarization", "orbital-hybridization",
                   "spin-orbital-even", "bond-charge", "bond-spin-orbital-even"))

#: the two time-reversal weights listed with the channel weights
TR_KEYS = ("TR-even", "TR-odd")

#: a channel names the label when its weight is at least this (O8)
DOMINANT_WEIGHT = 0.1

#: the parity is even (odd) when the TR-odd (TR-even) weight is below this;
#: otherwise it is mixed (O8)
TR_MINOR_WEIGHT = 0.1

#: a ``delta rho`` (or pattern) whose largest entry is at most this is
#: ``para`` (a converged para run leaves rounding noise, ~1e-12)
PARA_TOL = 1e-8

#: the gauge warning fires when the spin-independent part of a Transfer
#: entry, per spin, has ``|Im| >`` this (input energy units); see
#: :func:`gauge_dependent` for the normalisation
GAUGE_TOL = 1e-8

#: the thresholds a label was read with (recorded in every :class:`Label`
#: so that the report can show them)
THRESHOLDS = {"dominant_weight": DOMINANT_WEIGHT, "tr_minor_weight": TR_MINOR_WEIGHT,
              "para_tol": PARA_TOL, "gauge_tol": GAUGE_TOL}

#: ``delta rho`` must be Hermitian to this fraction of its norm
HERMITIAN_TOL = 1e-8

#: diagonal ``c^i`` entries below this fraction of the largest one do not
#: decide between ferro- and antiferro-orbital
_QUALIFIER_FLOOR = 0.1

#: the block orders :func:`decompose` reads
LAYOUTS = ("interleaved", "spin_block")

_MU = ("c0", "cx", "cy", "cz")
_SIGMA = np.array([np.eye(2),
                   [[0, 1], [1, 0]],
                   [[0, -1j], [1j, 0]],
                   [[1, 0], [0, -1]]], dtype=np.complex128)


@dataclass(frozen=True)
class Label:
    """A label (design section 6.5).

    ``name``: the readable display name (may contain non-ASCII text).
    ``weights``: every channel of :data:`CHANNELS` (they sum to 1; all 0
    for ``para``) and the two :data:`TR_KEYS` (which sum to 1 as well).
    ``tr``: ``"even"``, ``"odd"`` or ``"mixed"`` (the :data:`TR_MINOR_WEIGHT`
    rule). ``gauge_dependent``: the Transfer table has a complex
    spin-independent part (:func:`gauge_dependent`), so ``tr`` and the
    current/hybridization split depend on the orbital gauge.
    ``ascii_name``: :func:`ascii_name` of ``name`` (for id slugs).
    ``source``: ``"pattern"`` (:func:`pattern_label`) or ``"state"``
    (:func:`label`). ``thresholds``: :data:`THRESHOLDS`, as used.
    """
    name: str
    weights: Dict[str, float]
    tr: str
    gauge_dependent: bool
    ascii_name: str = ""
    source: str = "state"
    thresholds: Mapping[str, float] = None


# --------------------------------------------------------------------------
# the sigma (x) orbital-pair basis
# --------------------------------------------------------------------------

def _check_layout(layout):
    if layout not in LAYOUTS:
        raise ValueError("layout: must be one of {}, got {!r}".format(LAYOUTS, layout))


def _as_asbt(block, layout):
    """``block`` ``(..., 2n, 2n)`` as ``O[..., a, s, b, t]``."""
    _check_layout(layout)
    b = np.asarray(block, dtype=np.complex128)
    if b.ndim < 2 or b.shape[-1] != b.shape[-2] or b.shape[-1] % 2:
        raise ValueError("a block must be square with an even size (2n, 2n), "
                         "got shape {}".format(b.shape))
    n = b.shape[-1] // 2
    lead = b.shape[:-2]
    if layout == "interleaved":
        return b.reshape(lead + (n, 2, n, 2))
    return np.moveaxis(b.reshape(lead + (2, n, 2, n)), (-4, -3, -2, -1), (-3, -4, -1, -2))


def decompose(block, *, layout="interleaved"):
    """The coefficients ``c^mu`` of ``block = 1/2 sum_mu c^mu (x) sigma^mu``.

    Parameters
    ----------
    block : array ``(..., 2n, 2n)``
        Leading axes are kept (e.g. one block per displacement).
    layout : ``"interleaved"`` or ``"spin_block"``
        The index order of ``block``: ``2*a + s`` (spin-orbital mode; the
        literal ``numpy.kron(c, sigma)``) or ``s*n + a`` (normal mode).

    Returns
    -------
    dict
        ``{"c0", "cx", "cy", "cz"}``, each ``(..., n, n)``, with ``c^mu_ab
        = sum_st O_{(as),(bt)} sigma^mu_ts``.
    """
    o = _as_asbt(block, layout)
    return {k: np.einsum("...asbt,ts->...ab", o, _SIGMA[m]) for m, k in enumerate(_MU)}


def compose(c, *, layout="interleaved"):
    """The block ``1/2 sum_mu c^mu (x) sigma^mu`` (the inverse of
    :func:`decompose`; missing keys are zero)."""
    _check_layout(layout)
    arrs = [np.asarray(c[k], dtype=np.complex128) for k in _MU if k in c]
    if not arrs:
        raise ValueError("compose: no coefficient given (keys c0, cx, cy, cz)")
    shape = arrs[0].shape
    o = 0.5 * sum(np.einsum("...ab,st->...asbt", np.asarray(c[k], dtype=np.complex128),
                            _SIGMA[m]) for m, k in enumerate(_MU) if k in c)
    n = shape[-1]
    lead = shape[:-2]
    if layout == "interleaved":
        return o.reshape(lead + (2 * n, 2 * n))
    return np.moveaxis(o, (-3, -4, -1, -2), (-4, -3, -2, -1)).reshape(lead + (2 * n, 2 * n))


def _tr_parts(c):
    """``(even, odd)`` squared norms of ``c`` under ``c^0 -> conj(c^0)``,
    ``c^i -> -conj(c^i)``."""
    even = odd = 0.0
    for k in _MU:
        if k not in c:
            continue
        v = np.asarray(c[k], dtype=np.complex128)
        re, im = float(np.sum(v.real ** 2)), float(np.sum(v.imag ** 2))
        if k == "c0":
            even, odd = even + re, odd + im
        else:
            even, odd = even + im, odd + re
    return even, odd


def _parity(even, odd):
    """The O8 rule on the squared norms of the two parities."""
    total = even + odd
    if total == 0.0 or odd < TR_MINOR_WEIGHT * total:
        return "even"
    if even < TR_MINOR_WEIGHT * total:
        return "odd"
    return "mixed"


def tr_parity(c) -> Literal["even", "odd", "mixed"]:
    """The time-reversal parity of the coefficients ``c`` (a
    :func:`decompose` dict; missing keys are zero) in a real orbital gauge.

    ``Theta = i sigma^y K`` maps ``c^0 -> conj(c^0)`` and ``c^i ->
    -conj(c^i)``. ``"even"``: the TR-odd part (imaginary ``c^0``, real
    ``c^i``; real-symmetric ``c^0`` and imaginary-antisymmetric ``c^i`` are
    TR-even on a Hermitian block) carries less than
    :data:`TR_MINOR_WEIGHT` of the squared norm; ``"odd"``: the converse;
    ``"mixed"`` otherwise (owner decision O8). Zero coefficients are
    ``"even"``.
    """
    return _parity(*_tr_parts(c))


# --------------------------------------------------------------------------
# the gauge warning
# --------------------------------------------------------------------------

def _problem(obj):
    return getattr(obj, "problem", obj)


def _layout(problem):
    return "interleaved" if problem.enable_spin_orbital else "spin_block"


def gauge_dependent(problem) -> bool:
    """Whether the labels of ``problem`` are gauge-dependent: the
    spin-independent part of some entry of the whole Transfer table (every
    R) has ``|Im| > GAUGE_TOL``.

    The table is ``problem.ham_trans_r`` (the folded Transfer, dense, in
    the solver's index order), decomposed per entry by :func:`decompose`
    (ruling R5, for 9-column and spin-orbital tables). **Normalisation.**
    With ``O = 1/2 sum c^mu sigma^mu`` the coefficient is ``c^0 = T_upup +
    T_dndn``, twice the spin-independent hopping. The threshold applies to
    the hopping per spin, ``c^0 / 2 = (T_upup + T_dndn) / 2`` (R5), so a
    spin-independent ``t`` warns when ``|Im t| > GAUGE_TOL``; in terms of
    ``c^0`` the test is ``|Im c^0| > 2 GAUGE_TOL``. Spin parts (e.g. the
    SOC ``i lambda sigma^z``, an imaginary ``c^z``) are exempt.

    ``problem``: a :class:`~hwave.solver.initial_facade.UHFkProblem` or a
    :class:`~hwave.initial.facts.ModelFacts`. The result is cached per
    problem (its arrays are read-only).
    """
    p = _problem(problem)
    try:
        return _GAUGE_CACHE[p]
    except (KeyError, TypeError):
        pass
    c0 = decompose(p.ham_trans_r, layout=_layout(p))["c0"]
    out = bool(c0.size and float(np.max(np.abs(c0.imag))) / 2 > GAUGE_TOL)
    try:
        _GAUGE_CACHE[p] = out
    except TypeError:   # not weak-referenceable: no cache
        pass
    return out


#: :func:`gauge_dependent` per problem (weak keys: a problem is never kept)
_GAUGE_CACHE = weakref.WeakKeyDictionary()


# --------------------------------------------------------------------------
# the label
# --------------------------------------------------------------------------

def _delta_rho(rho, problem, what="delta rho"):
    """``rho`` as ``(nvol, nd, nd)``, checked to be Hermitian."""
    r = np.asarray(rho, dtype=np.complex128)
    nvol, nd = problem.nvol, problem.nd
    shapes = ((nvol, nd, nd), (nvol, problem.ns, problem.norb, problem.ns, problem.norb))
    if r.shape not in shapes:
        raise ValueError("{} of this problem has shape {} or {}, got {}".format(
            what, shapes[0], shapes[1], r.shape))
    m = r.reshape(nvol, nd, nd)
    # m - (m + partner)/2 = (m - partner)/2
    g5 = m.reshape(shapes[1])
    residual = 2 * float(np.linalg.norm(g5 - hermitian_green(g5, problem.shape)))
    if residual > HERMITIAN_TOL * float(np.linalg.norm(m)):
        raise ValueError("{} is not Hermitian (X(r)_ij != conj(X(-r)_ji)); pass a "
                         "density matrix in the rho = G^T convention, or a Hermitian "
                         "pattern".format(what))
    return m


def _folded_sites(facts, n):
    """The site of each folded physical orbital ``A = a + norb_phys * ir``:
    ``(site of a, ir)``."""
    nphys = int(facts.norb_phys)
    if n % nphys:
        raise ValueError("the folded cell has {} physical orbitals, not a multiple of "
                         "norb_phys = {}".format(n, nphys))
    site_of = np.empty(nphys, dtype=int)
    for g, members in enumerate(facts.site_groups.groups):
        for a in members:
            site_of[a] = g
    A = np.arange(n)
    return site_of[A % nphys] + len(facts.site_groups.groups) * (A // nphys)


def _channel_norms(c, sites):
    """The squared norm of every channel of :data:`CHANNELS`."""
    n = sites.size
    onsite = np.zeros(c["c0"].shape, dtype=bool)
    onsite[0] = sites[:, None] == sites[None, :]
    diag = np.zeros_like(onsite)
    diag[0] = np.eye(n, dtype=bool)
    off = onsite & ~diag
    bond = ~onsite
    w = dict.fromkeys(CHANNELS, 0.0)

    def sq(x):
        return float(np.sum(np.abs(x) ** 2))

    d0 = np.diagonal(c["c0"][0]).real
    for s in np.unique(sites):
        v = d0[sites == s]
        w["charge"] += v.size * float(np.mean(v)) ** 2
        w["orbital-polarization"] += sq(v - np.mean(v))
    # imaginary diagonal parts vanish for a Hermitian delta rho
    w["orbital-hybridization"] = sq(c["c0"].real[off])
    w["orbital-current"] = sq(c["c0"].imag[off])
    w["bond-charge"] = sq(c["c0"].real[bond])
    w["bond-current"] = sq(c["c0"].imag[bond])
    for ax in "xyz":
        ci = c["c" + ax]
        w["spin-" + ax] = sq(ci.real[diag])
        w["inter-orbital-spin-" + ax] = sq(ci.real[off])
        w["spin-orbital-even"] += sq(ci.imag[off])
        w["bond-spin-" + ax] = sq(ci.real[bond])
        w["bond-spin-orbital-even"] += sq(ci.imag[bond])
    return w


def _orbital_qualifier(ci, sites):
    """``"ferro-orbital"`` when the diagonal on-site ``c^i`` has one sign on
    every site (entries below ``_QUALIFIER_FLOOR`` of the largest ignored),
    ``"antiferro-orbital"`` otherwise."""
    d = np.diagonal(ci[0]).real
    big = np.abs(d) > _QUALIFIER_FLOOR * float(np.max(np.abs(d)))
    for s in np.unique(sites):
        v = d[(sites == s) & big]
        if v.size and (v.min() < 0 < v.max()):
            return "antiferro-orbital"
    return "ferro-orbital"


def _q_text(q):
    if len(q) != 3:
        raise ValueError("q: must have three components, got {!r}".format(q))
    for x in q:
        if isinstance(x, bool) or not isinstance(x, (str, numbers.Rational)):
            raise ValueError("q: components must be rational (a 'p/q' string, an "
                             "integer or a Fraction), got {!r}".format(x))
    return "q=({})".format(",".join(str(f) for f in q_fractions(q)))


def _eta_text(eta):
    if isinstance(eta, str):
        return eta
    eta = tuple(eta)
    if not eta or any(isinstance(e, bool) or e not in (1, -1) for e in eta):
        raise ValueError("eta: must be a name or a sequence of +1/-1, got {!r}".format(eta))
    if eta == (1,):
        return None
    return "η=({})".format(",".join("+" if e == 1 else "-" for e in eta))


def _phase_text(phase):
    if isinstance(phase, bool) or not isinstance(phase, numbers.Real) \
            or not math.isfinite(phase):
        raise ValueError("phase: must be a finite number, got {!r}".format(phase))
    f = Fraction(phase / math.pi).limit_denominator(48)
    if abs(phase - float(f) * math.pi) <= 1e-9:
        if f == 0:
            return "φ=0"
        num = {1: "", -1: "-"}.get(f.numerator, str(f.numerator))
        return "φ={}π{}".format(num, "" if f.denominator == 1 else "/{}".format(f.denominator))
    return "φ={:.6g}".format(phase)


_ASCII = {"φ": "phi", "π": "pi", "η": "eta"}


def ascii_name(name: str) -> str:
    """``name`` in ASCII: ``φ``, ``π``, ``η`` spelled ``phi``, ``pi``,
    ``eta``; accents folded (``Néel`` -> ``Neel``); any other non-ASCII
    character dropped. Id slugs are built from this form."""
    text = "".join(_ASCII.get(ch, ch) for ch in name)
    text = unicodedata.normalize("NFKD", text)
    return text.encode("ascii", "ignore").decode("ascii")


def _name_suffix(q, eta, phase, sign, spiral):
    parts = []
    if q is not None:
        parts.append(_q_text(q))
    if eta is not None:
        text = _eta_text(eta)
        if text:
            parts.append(text)
    if phase is not None:
        parts.append(_phase_text(phase))
    if sign is not None:
        if isinstance(sign, bool) or sign not in (1, -1):
            raise ValueError("sign: must be 1, -1 or None, got {!r}".format(sign))
        if sign == -1:
            parts.append("sign=-1")
    if spiral:
        parts.append("spiral")
    return parts


def _dominant_tokens(weights, c, sites):
    """The dominant channels (weight >= DOMINANT_WEIGHT): on-site channels
    first, then bond channels, each heaviest first (ties in CHANNELS
    order); the heaviest channel alone when none is dominant."""
    order = sorted(range(len(CHANNELS)), key=lambda i: (-round(weights[CHANNELS[i]], 12), i))
    dominant = [CHANNELS[i] for i in order if weights[CHANNELS[i]] >= DOMINANT_WEIGHT]
    dominant = ([ch for ch in dominant if not ch.startswith("bond-")]
                + [ch for ch in dominant if ch.startswith("bond-")])
    tokens = []
    for ch in dominant or [CHANNELS[order[0]]]:
        if ch in ("spin-x", "spin-y", "spin-z"):
            ch = "{} {}".format(ch, _orbital_qualifier(c["c" + ch[-1]], sites))
        tokens.append(ch)
    return " + ".join(tokens)


def _label(m, facts, source, q, eta, phase, sign, spiral):
    problem = facts.problem
    suffix = _name_suffix(q, eta, phase, sign, spiral)
    if float(np.max(np.abs(m), initial=0.0)) <= PARA_TOL:
        weights = dict.fromkeys(CHANNELS + TR_KEYS, 0.0)
        head, tr = "para", "even"
    else:
        c = decompose(m, layout=_layout(problem))
        sites = _folded_sites(facts, c["c0"].shape[-1])
        norms = _channel_norms(c, sites)
        total = sum(norms.values())
        even = sum(v for k, v in norms.items() if k in _EVEN)
        odd = total - even
        weights = {k: v / total for k, v in norms.items()}
        weights["TR-even"], weights["TR-odd"] = even / total, odd / total
        head, tr = _dominant_tokens(weights, c, sites), _parity(even, odd)
    name = " ".join([head] + suffix)
    return Label(name=name, weights=weights, tr=tr, gauge_dependent=gauge_dependent(problem),
                 ascii_name=ascii_name(name), source=source, thresholds=dict(THRESHOLDS))


_NAME_DOC = """
    q : three rationals or None
        Written ``q=(1/2,1/2,0)``; omitted when None.
    eta : str, sequence of +1/-1, or None
        The site sign pattern: a name (``"Néel"``) is written as given, a
        sequence as ``η=(+,-)`` (a single ``+1`` is omitted).
    phase : float or None
        Written ``φ=0``, ``φ=π/4`` (a multiple of pi with a denominator up
        to 48) or ``φ=0.123457``; omitted when None.
    sign : +1, -1 or None
        ``-1`` is written ``sign=-1``; ``+1`` and None are omitted.
    spiral : bool
        ``True`` appends ``spiral``; a collinear pattern is the unmarked
        default.

    The name lists the dominant channels, those with a weight of at least
    :data:`DOMINANT_WEIGHT` (the heaviest one if none is): on-site channels
    first, then bond channels, each group heaviest first, joined by
    ``" + "``. A ``spin-<axis>`` channel is followed by ``ferro-orbital``
    or ``antiferro-orbital``. ``tr`` follows the :data:`TR_MINOR_WEIGHT`
    rule. An input whose largest entry is at most :data:`PARA_TOL` is
    named ``para``. The thresholds are recorded in ``Label.thresholds``.
"""


def pattern_label(pattern_k, facts, *, q, eta, phase, sign, spiral: bool) -> Label:
    """The label of a seeding pattern ``O`` (owner decision O8 (a): a
    candidate and its id slug are named by its pattern, not by the induced
    ``delta rho``).

    Parameters
    ----------
    pattern_k : array ``(nvol, nd, nd)``
        The Hamiltonian-convention ``O(k)`` of the folded cell of
        ``facts.problem`` (as ``realise`` returns it), in the solver's
        index order; it is read in real space, ``O(k) = sum_R e^{ikR}
        O(R)``. ``O(R)`` must be Hermitian (``O(R)_ab = conj(O(-R)_ba)``).
    facts : :class:`~hwave.initial.facts.ModelFacts`
        Of the cell the pattern lives on.
    """
    problem = facts.problem
    m = _delta_rho(_to_r(pattern_k, problem).reshape(problem.nvol, problem.nd, problem.nd),
                   problem, what="the pattern")
    return _label(m, facts, "pattern", q, eta, phase, sign, spiral)


def label(rho_delta_cell, facts, *, q, eta, phase, sign, spiral: bool) -> Label:
    """The label of an induced or converged state, from its ``delta rho``
    (design section 6.5; owner decision O8 (b): named by its dominant
    components only).

    Parameters
    ----------
    rho_delta_cell : array
        ``delta rho`` of the folded cell of ``facts.problem`` in the rho
        convention (``rho = G^T`` per displacement; the callers transpose
        a stored Green function), ``(nvol, nd, nd)`` or ``(nvol, ns, norb,
        ns, norb)`` in the solver's index order. It must be Hermitian
        (``rho(r)_ij = conj(rho(-r)_ji)``).
    facts : :class:`~hwave.initial.facts.ModelFacts`
        Of the cell ``delta rho`` lives on (its sites decide what is
        on-site and what is a bond).
    """
    problem = facts.problem
    return _label(_delta_rho(rho_delta_cell, problem), facts, "state", q, eta, phase, sign,
                  spiral)


if pattern_label.__doc__ and label.__doc__:   # absent under python -OO
    pattern_label.__doc__ += _NAME_DOC
    label.__doc__ += _NAME_DOC
