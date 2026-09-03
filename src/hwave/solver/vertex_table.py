"""The adjudicated spin/charge vertex content, single-sourced.

Every value here was established by exact diagonalization of H-wave's
documented file operators (issues #105/#113; methodology in #107): the
O(lambda) self-energy is removed exactly via the Hartree-Fock response,
the extraction is validated by q-independence for on-site input, and the
per-channel scale is calibrated on the CoulombIntra control. This module
is the ONE place that content lives; the two consumers -- the pair-space
S/C builders in :mod:`hwave.sc` (FLEX general and Eliashberg) and the
Fierz CROSS-SLOT correction of the ring vertex in
:mod:`hwave.solver.rpa` (#104) -- derive their coefficients from it, so
those pieces cannot drift apart again (they did, repeatedly: #96, #104,
#113). The ring's BASE content is covered too
(:func:`ring_spin_table`): its density-slot weights are the channel
decomposition of the same (S, C) entries, executed in code, and the
transverse spin-flip weights (adjudicated in #105, underivable from the
longitudinal channels) are recorded here as data.

Slot families of the pair-space matrices (slots are (l1*norb+l2,
l3*norb+l4)):

    ``diag``      (aa, aa)   same orbital on both pairs
    ``cross``     (ab, ab)   pair-diagonal, orbital off-diagonal (Fierz)
    ``density``   (aa, bb)   two densities
    ``antidiag``  (ab, ba)   pair-antidiagonal

``a != b`` throughout, with ONE exception: CoulombInter's ``density``
entry also applies at a == b, where an inter-site same-orbital V
contributes 2 V_aa(q) to the charge diagonal (issue #95); Hund and
Ising stay restricted to a != b there for ON-site input (an orbital has
no on-site Hund or Ising coupling with itself), while an inter-site
same-orbital Hund/Ising is a physical density-density term and takes
the ``density`` entry at a == b as well (#181 Tier 1, the locality
split in :mod:`hwave.sc`).

OFF-SITE content (R != 0): only the ``density`` family is
q-representable (Hartree vertex V_ab(q)); the ``cross``/``antidiag``
families describe a non-local pair off-site and are on-site only. In
particular, off-site Exchange and PairLift have NO longitudinal entry:
adjudicated by exact diagonalization (#181 Tier 2,
tests/test_offsite_exchange_ed_longitudinal.py) -- Exchange's density
slots measure 1e-3 against the controls' 2.0 / 1.0 per unit coupling and
a purely imaginary J gives an identically vanishing longitudinal
response; PairLift's longitudinal response vanishes off-site as it does
on-site. Their content is transverse (RING_SPIN_FLIP below), which is
why the table carries no off-site entry for either.

``S`` enters the spin channel as ``[1 - chi0 S]^-1 chi0`` and ``C`` the
charge channel as ``[1 + chi0 C]^-1 chi0``. Entries are per unit
coupling of the interaction file.

The ring's rank-4 vertex is spin-resolved instead of channel-resolved;
for a slot with content (S, C) the same-spin / cross-spin coefficients
in the solver's repulsion-positive convention are the channel
decomposition

    W_same  = (C - S) / 2
    W_cross = (C + S) / 2

anchored on the two families the ring always had right (CoulombIntra
diag and CoulombInter density; see the #104 fix). Executing this
decomposition in code -- rather than recording it in a commit message --
is the point of this module.
"""

from types import MappingProxyType

# (S, C) per unit coupling, per slot family. Types absent from a family
# contribute nothing there. PairLift has no particle-hole S/C content
# (adjudicated zero). Read-only: a consumer mutating the table would
# silently poison BOTH builders at once, which is worse than the drift
# this module exists to prevent.
_ADJUDICATED_SC_RAW = {
    "CoulombIntra": {"diag": (+1.0, +1.0)},
    "CoulombInter": {"cross": (+1.0, -1.0), "density": (0.0, +2.0)},
    "Hund":         {"cross": (-1.0, +1.0), "density": (+1.0, -1.0)},
    "Exchange":     {"cross": (+1.0, +1.0)},
    "Ising":        {"cross": (+1.0, -1.0), "density": (-2.0, 0.0)},
    "PairLift":     {},
    "PairHop":      {"antidiag": (+1.0, +1.0)},
}

ADJUDICATED_SC = MappingProxyType(
    {k: MappingProxyType(v) for k, v in _ADJUDICATED_SC_RAW.items()})
# no mutable backing name is retained: the proxies above hold the only
# references to the inner dicts, and the module namespace exposes only
# the read-only view
del _ADJUDICATED_SC_RAW


def sc_coefficients(itype, family):
    """(S, C) coefficients of ``itype`` at ``family``, or (0, 0)."""
    return ADJUDICATED_SC.get(itype, {}).get(family, (0.0, 0.0))


def fierz_coefficients(itype):
    """(W_same, W_cross) of the ring's cross-slot correction for ``itype``.

    The channel decomposition of the ``cross`` family content; zero for
    types without cross-slot content.
    """
    s, c = sc_coefficients(itype, "cross")
    return ((c - s) / 2.0, (c + s) / 2.0)


# Transverse spin-flip weights of the ring's base tensor, adjudicated in
# issue #105 (the transverse ladder channel was measured against exact
# diagonalization independently of the longitudinal channels). These are
# NOT derivable from the (S, C) table above: the longitudinal channels see
# only spin-conserving bilinears, while these entries flip spin on each
# bilinear -- Exchange couples the spin-raising density (a,a) to the
# spin-lowering (b,b), PairLift raises on both. Keys are the ring's
# ``spin(a, ap, bp, b)`` convention (0: up, 1: down); read-only for the
# same reason as ADJUDICATED_SC.
_RING_SPIN_FLIP_RAW = {
    "Exchange": {(0, 1, 1, 0): -1.0, (1, 0, 0, 1): -1.0},
    "PairLift": {(0, 1, 0, 1): +1.0, (1, 0, 1, 0): +1.0},
}

RING_SPIN_FLIP = MappingProxyType(
    {k: MappingProxyType(v) for k, v in _RING_SPIN_FLIP_RAW.items()})
del _RING_SPIN_FLIP_RAW


def ring_spin_table(itype):
    """Spin-resolved base weights of the ring vertex for ``itype``.

    Returns ``{(s1, s2, s3, s4): w}`` in the ring's ``spin(a, ap, bp, b)``
    convention. The spin-conserving part is the channel decomposition
    ``W_same = (C - S)/2``, ``W_cross = (C + S)/2`` of the type's
    adjudicated (S, C) content at the slot family its base placement
    occupies -- ``diag`` for CoulombIntra (same orbital on both
    bilinears), ``antidiag`` for PairHop (its ring placement couples
    (a,b) to (b,a)), ``density`` for everything else. The spin-flip part
    is the #105 transverse content, recorded data (see RING_SPIN_FLIP).

    This reproduces the hand-coded per-type spin tables the ring carried
    since its adjudication, now derived from the one table instead of
    duplicated next to it.
    """
    if itype == "CoulombIntra":
        family = "diag"
    elif itype == "PairHop":
        family = "antidiag"
    else:
        family = "density"
    s, c = sc_coefficients(itype, family)
    w_same = (c - s) / 2.0
    w_cross = (c + s) / 2.0
    table = {}
    if w_same != 0.0:
        table[(0, 0, 0, 0)] = w_same
        table[(1, 1, 1, 1)] = w_same
    if w_cross != 0.0:
        table[(0, 0, 1, 1)] = w_cross
        table[(1, 1, 0, 0)] = w_cross
    table.update(RING_SPIN_FLIP.get(itype, {}))
    return table
