"""The adjudicated spin/charge vertex content, single-sourced.

Every value here was established by exact diagonalization of H-wave's
documented file operators (issues #105/#113; methodology in #107): the
O(lambda) self-energy is removed exactly via the Hartree-Fock response,
the extraction is validated by q-independence for on-site input, and the
per-channel scale is calibrated on the CoulombIntra control. This module
is the ONE place that content lives; the two consumers -- the pair-space
S/C builders in :mod:`hwave.sc` (FLEX general and Eliashberg) and the
spin-resolved rank-4 ring vertex in :mod:`hwave.solver.rpa` (#104) --
derive their coefficients from it, so they cannot drift apart again
(they did, repeatedly: #96, #104, #113).

Slot families of the pair-space matrices (slots are (l1*norb+l2,
l3*norb+l4)):

    ``diag``      (aa, aa)   same orbital on both pairs
    ``cross``     (ab, ab)   pair-diagonal, orbital off-diagonal (Fierz)
    ``density``   (aa, bb)   two densities
    ``antidiag``  (ab, ba)   pair-antidiagonal

``a != b`` throughout, with ONE exception: CoulombInter's ``density``
entry also applies at a == b, where an inter-site same-orbital V
contributes 2 V_aa(q) to the charge diagonal (issue #95); Hund and
Ising stay restricted to a != b there (an orbital has no Hund or Ising
coupling with itself).

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

# (S, C) per unit coupling, per slot family. Types absent from a family
# contribute nothing there. PairLift has no particle-hole S/C content
# (adjudicated zero).
ADJUDICATED_SC = {
    "CoulombIntra": {"diag": (+1.0, +1.0)},
    "CoulombInter": {"cross": (+1.0, -1.0), "density": (0.0, +2.0)},
    "Hund":         {"cross": (-1.0, +1.0), "density": (+1.0, -1.0)},
    "Exchange":     {"cross": (+1.0, +1.0)},
    "Ising":        {"cross": (+1.0, -1.0), "density": (-2.0, 0.0)},
    "PairLift":     {},
    "PairHop":      {"antidiag": (+1.0, +1.0)},
}


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
