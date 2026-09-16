#!/usr/bin/env python3

"""G4: exact diagonalisation of the OFF-SITE second order on a chain.

Task 10 (``tests/test_flex_second_order_ed_onsite.py``, gate G3) settled the
ON-SITE second order on a single site, where every interaction entry is
local and the locality split of ``flex_second_order = "local"`` has nothing
to drop. This module is the off-site half: on periodic chains
(``CellShape = [L, 1, 1]``) with a NEIGHBOUR bond declared, the O(v^2)
coefficient of the exact fluctuation remainder

    Sigma_ED - Sigma_HF[rho_ED]

must equal the production ``Sigma_fluct`` of the general path PLUS the
off/off uncrossed exchange class that the local kernel deliberately drops
(``tests/test_second_order_oracle.py``'s ``"dropped"`` weighting); and with
the bond gate on -- which resums that class -- production alone must equal
the ED remainder.

Two fixtures, and what each adjudicates
--------------------------------------
* ``_fx_chain()`` -- ``L = 4``, one orbital, ``t = -1``, an on-site ``U``
  and a nearest-neighbour ``V``. Asserted: the O(V^2), O(U V) and O(U^2)
  coefficients of ``local + dropped`` and, for ``V``, the gate-on
  coefficient. Off-site Hund and Ising on the same chain are RECORDED
  (spec G2 (c) / G4), not asserted.
* ``_fx_orbital()`` -- ``L = 3``, TWO orbitals, and an ORBITAL-ASYMMETRIC
  inter-orbital bond (``v_01(+x) = v``, ``v_10(+x) = 0.6 v``). This is the
  fixture that adjudicates the two open ledger items of the campaign: the
  orientation of the off-site density vertex (the documented reading, the
  orientation fix of issue #193) and the bond gate's deviation on
  inter-orbital bonds (Task 9, issue #192). Its verdicts are in
  :class:`TestG4`'s methods and in the module's measured table below. The
  same fixture also carries the MIXED on-site/off-site entry (on-site
  ``Hund`` against the off-site ``V`` rows,
  :meth:`TestG4.test_hund_times_offsite_v_chain`).

Which first order is subtracted
-------------------------------
The remainder the gates compare is built with the TEST SIDE's own
first-order functional (:func:`_hf_sigma`, the Wick derivative of the ED
term list, shared with Task 10), never with production's -- so it is not
built with any part of the object under test. UHFk's own mean-field kernel
appears only as a CONVENTION PIN, in
:meth:`TestChainHamiltonian.test_ed_hamiltonian_is_the_production_mean_field`
(which is what fixes the orbital placement and the mirrored-row weight),
:meth:`TestChainHamiltonian.test_independent_functional_equals_production_on_random_densities`
(the two functionals agree away from the free density too) and
:meth:`TestFirstOrder.test_first_order_remainder_vanishes` (which subtracts
the production kernel precisely so that it can be two-sided in the weight;
the test-side functional tracks whichever Hamiltonian it is handed and
cannot discriminate one).

Which Hamiltonian the chain ED must be
--------------------------------------
The declared row ``(R, a, b, v)`` of an off-site two-body type places
orbital ``a`` in the ORIGINAL cell and orbital ``b`` in the cell displaced
by ``R``,

    v n_{j, a} n_{j+R, b}

-- the DOCUMENTED reading (``docs/en`` UHFk interaction file page:
"``[alpha]`` corresponds to the orbital alpha in the original cell, and
``[beta]`` corresponds to the orbital beta in the cell displaced by r"),
which the mean-field kernel implements since the orientation fix of issue
#192. Note that this is the OPPOSITE placement from the one a
``transfer.dat`` row gets in this code base -- ``EDFixture.build_h1`` puts
``t[(a, b)]`` on ``(j + R, a) <- (j, b)``, which
:meth:`TestChainHamiltonian.test_bare_green_matches_the_production_solver`
pins against the solver's own band. The one-body and two-body readings
therefore do not agree, and that asymmetry is what the fix makes explicit.

This is NOT assumed here. :meth:`TestChain
Hamiltonian.test_ed_hamiltonian_is_the_production_mean_field` compares the
first-order functional of the ED term list against UHFk's OWN mean-field
kernel (``hwave.solver.hartree_fock.accumulate_hf``, Hartree AND Fock, at
the same density) and fails with the opposite (orbital-swapped) placement
by 5.4e-2 / 1.4e-2 of the mean field on the two inter-orbital fixtures,
while the documented placement matches at 1e-16. On the single-orbital
chain the two placements are the same Hamiltonian, which is why only the
two-orbital fixture can decide it.
:meth:`TestChainHamiltonian.test_first_order_gate_every_type_documented_reading`
extends that adjudication to EVERY off-site type the kernel accepts
(``CoulombInter``, ``Hund``, ``Ising``, ``Exchange``, ``PairLift``,
``PairHop``) and to both settings of the kernel's ``include_fock`` switch,
on spin-coherent densities -- which the spin-flip types need in order to
have a mean field at all.

SCOPE OF THAT CLAIM, and of every verdict below it. What is adjudicated
here is the second-order vertex against UHFk's MEAN-FIELD READING of a
declared row: ``accumulate_hf`` is what defines the orbital placement the
chain ED implements, and the gate then asks whether the second order is
consistent with THAT Hamiltonian. Since the orientation fix that reading
IS the documented one, so the pin is now anchored to the file-format
definition and not only to internal consistency; what the module still
cannot decide on its own is whether the DOCUMENTATION is the intended
convention. The k-space bridge the comparison rides on (exponent sign,
orbital index order) is pinned separately, on a band that can see it, by
:meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`
-- without that pin a flipped bridge would silently select the opposite
placement and invert every verdict below.

The mirrored-row HALF weight (each of the two declared orientations of a
bond carries half of it, ``_MIRRORED_TYPES``) is the convention Task 10
established on the on-site types; here it is re-pinned two-sidedly by
:meth:`TestFirstOrder.test_first_order_remainder_vanishes` -- vanishing at
the half weight for every fixture, and 0.48 ... 0.53 of the self-energy's
own linear size at the full weight.

Stencils, working point, tolerance
----------------------------------
The coefficient extraction is Task 10's, imported rather than recopied:
the two-point stencils ``_coeff2`` / ``_coeff11`` with the O(x) refinement
``2 c(x/2) - c(x)``. Working point: ``nmat = 4096``, ``x = 3.125e-3``, and
the innermost eighth of the Matsubara grid -- Task 10's window, every one
of its 512 frequencies on both sides. The two error sources are the same
ones Task 10 measured: the finite ``x`` truncation (which the refinement
leaves at O(x^2)) and the production kernel's finite Matsubara window
(which is a fixed fraction of the coefficient and shrinks as ``1/nmat``:
9.2e-4 at nmat = 1024, 4.3e-4 at 2048, 2.3e-4 at 4096 on the V^2 entry).

Measured at that working point (relative to the ED coefficient's own
maximum over the window; tolerance 2e-3):

    fixture     entry              local + dropped     gate on
    L=4 norb=1  V^2                2.3e-4              2.4e-4
    L=4 norb=1  U V                4.8e-4              --
    L=4 norb=1  U^2                3.2e-4              --
    L=4 norb=1  Hund   (recorded)  2.6e-4              2.6e-4
    L=4 norb=1  Ising  (recorded)  2.2e-4              2.2e-4
    L=3 norb=2  v^2 (asymmetric)   3.5e-4              4.0e-2  <-- recorded
    L=3 norb=2  J V (mixed)        4.5e-4              --

The v^2 row is the campaign's open item: the STANDALONE local path plus
the dropped class reproduces the exact inter-orbital off-site second order
(so the DOCUMENTED orientation of the off-site density vertex -- orbital
``a`` in the original cell -- is confirmed by exact diagonalisation), while
the BOND GATE misses it by 4.0e-2 -- twenty
times the tolerance, and nearly twice the whole off/off uncrossed class it
is supposed to be resumming (2.2e-2 of the coefficient). That deviation is
a Phase B defect, tracked as issue #192; it is not a second-order kernel
finding and this module does not try to fix it.

It is RECORDED rather than asserted as an agreement -- but the record is
PINNED (:data:`_PINNED`): every printed quantity is held to its measured
value within :data:`_PINNED_BAND` and checked for finiteness, so it cannot
drift in either direction unnoticed. A value that moves OUT of the band is
a finding whichever way it moves.
"""
import os
import tempfile
import unittest

import numpy as np

from tests import ed_oracle_util as edu
from tests.heavy_tests import heavy
from tests.test_flex_second_order_ed_onsite import (_coeff2, _coeff11, _quiet, _refine,
                                                    wick_hf_sigma)
from tests.test_second_order_factors import _write_wan
from tests.test_second_order_oracle import (_MIRRORED_ROW_TYPES as _MIRRORED_TYPES,
                                             oracle_records, oracle_sigma2)

UP, DN = 0, 1

# ``_MIRRORED_TYPES`` above is the oracle's own list, imported rather than
# recopied: the types whose two declared orientations name the SAME bond, so
# each declared row is half of it. Pinned two-sidedly here by
# :meth:`TestFirstOrder.test_first_order_remainder_vanishes`.

#: The off-site rows of the per-type first-order gate
#: (:meth:`TestChainHamiltonian.test_first_order_gate_every_type_documented_reading`):
#: ONE inter-orbital nearest-neighbour bond per type, declared in both
#: orientations as a reader delivers it. Orbital 1 sits in the original
#: cell and orbital 2 in the cell displaced by ``+x``, so the orbital swap
#: maps each declaration onto a DIFFERENT Hamiltonian. ``PairHop`` carries a
#: complex amplitude (its mirrored row is the conjugate), which is the only
#: type whose orientation the kernel can get wrong on a REAL declaration.
_GATE_ROWS = {
    "CoulombInter": [(1, 0, 0, 1, 2, 0.4, 0.0), (-1, 0, 0, 2, 1, 0.4, 0.0)],
    "Hund": [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)],
    "Ising": [(1, 0, 0, 1, 2, 0.25, 0.0), (-1, 0, 0, 2, 1, 0.25, 0.0)],
    "Exchange": [(1, 0, 0, 1, 2, 0.2, 0.0), (-1, 0, 0, 2, 1, 0.2, 0.0)],
    "PairLift": [(1, 0, 0, 1, 2, 0.15, 0.0), (-1, 0, 0, 2, 1, 0.15, 0.0)],
    "PairHop": [(1, 0, 0, 1, 2, 0.1, 0.04), (-1, 0, 0, 2, 1, 0.1, -0.04)],
}

#: The on-site row the gate always carries alongside, so that every case
#: also exercises the mixed on-site/off-site accumulation.
_GATE_ONSITE = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.5, 0.0)]}

#: Types whose HARTREE-ONLY mean field cannot see the orientation of a
#: declared row at all. For each of them UHFk's direct term is built from
#: the EQUAL-SITE density ``gbb`` alone (``accumulate_hf``'s ``hh0``), so
#: on a translation-invariant density it is independent of the
#: displacement ``r`` -- and a Hermitian-closed declaration is invariant
#: under the orbital swap once ``r`` is dropped. The gate therefore asserts
#: EQUALITY of the two placements there instead of a miss; with the Fock
#: term on, every type including these misses by 9.4e-2 ... 4.5e-1.
#: ``PairHop`` is absent because its direct term reads the INTER-SITE
#: density, so it sees the orientation in both Fock settings.
_HARTREE_ORIENTATION_BLIND = ("CoulombInter", "Hund", "Ising", "Exchange",
                              "PairLift")

#: Working point of the heavy comparisons (see the module docstring).
_NMAT = 4096
_X = 3.125e-3
_WINDOW = _NMAT // 16          # half-width of the comparison window
_TOL = 2.0e-3

#: Anti-vacuity floor on the ED coefficient of a compared entry: a relative
#: comparison against a numerically zero reference would pass for any
#: production value at all. Every entry compared here measures
#: 5.8e-2 ... 3.6e-1.
_FLOOR = 1.0e-6

#: How large the dropped (off/off uncrossed) class must be, relative to the
#: coefficient, for the ``local + dropped`` assertion to be distinguishing
#: at all. Measured: 1.2e-1 for V^2 and 2.3e-1 for Hund on the
#: single-orbital chain, 2.2e-2 on the inter-orbital one.
_DROPPED_FLOOR = 1.0e-3

#: RECORDED quantities, pinned. These are outcomes this gate does not
#: ASSERT as agreements -- the bond gate's deviation on an inter-orbital
#: off-site bond (issue #192) and the off-site Hund/Ising results, which
#: belong to a follow-up class of issue #181 -- but a recorded number that
#: nothing checks is a number that can drift in either direction unnoticed.
#: Each entry is the measured value at this module's working point; the
#: band below is what a re-measurement may move it by.
#:
#: A change OUTSIDE the band is a finding either way: a larger deviation
#: means a regression, a smaller one means the defect has been partly fixed
#: (in which case the expectation, and issue #192, want updating).
_PINNED = {
    # inter-orbital off-site bond, L = 3 two-orbital chain
    "interorbital gate vs ED": 4.020e-2,          # issue #192
    "interorbital dropped class": 2.222e-2,
    # off-site Hund / Ising, L = 4 single-orbital chain
    "Hund local+dropped vs ED": 2.560e-4,
    "Hund gate vs ED": 2.560e-4,
    "Hund dropped class": 2.311e-1,
    "Ising local+dropped vs ED": 2.170e-4,
    "Ising gate vs ED": 2.170e-4,
    "Ising dropped class": 1.156e-1,
}

#: Relative band of :data:`_PINNED`. The pinned numbers are set by stencil
#: truncation and the finite Matsubara window, both deterministic, so the
#: band only has to absorb a different BLAS's round-off; it is deliberately
#: far tighter than the factor 20 that separates the issue-#192 deviation
#: from this module's tolerance.
_PINNED_BAND = 0.20


def _inner():
    """Indices of the comparison window on the production Matsubara grid --
    the innermost eighth, every frequency of it (Task 10's window)."""
    return np.arange(_NMAT // 2 - _WINDOW, _NMAT // 2 + _WINDOW)


def _fx_chain():
    """``L = 4``, one orbital, ``t = -1``: the brief's single-orbital chain."""
    return edu.EDFixture(L=4, norb=1, t={(0, 0): -1.0}, eps=(0.0,), T=0.5, mu=0.2)


def _fx_orbital():
    """``L = 3``, two orbitals, asymmetric intra/inter-orbital hopping and
    distinct on-site energies: the fixture that can see an orbital index
    order (``tests/test_bond_vs_ed_oracle.py``'s case-M shape, at the same
    ``T`` and ``mu``). ``nmode = 12``, Fock dimension 4096, largest
    ``(N_up, N_dn)`` block 400 -- comfortably inside ``SectorED``."""
    return edu.EDFixture(L=3, norb=2,
                         t={(0, 0): -1.0, (1, 1): -0.7, (0, 1): -0.2, (1, 0): -0.2},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


def _fx_complex():
    """``L = 3``, two orbitals, COMPLEX hopping in every orbital channel and
    ``t_01 != conj(t_10)``: a band that is neither inversion-symmetric nor
    transpose-symmetric.

    The two fixtures the gate itself runs on both have real hopping, so
    ``eps(k) = eps(-k)`` and ``eps(k) = eps(k)^T`` hold to round-off there
    and NEITHER the Fourier sign nor the orbital index order of
    :func:`_to_k` can be seen. This fixture was introduced to see them
    (:meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`).

    It has a second consumer:
    :meth:`TestChainHamiltonian.test_first_order_gate_every_type_documented_reading`
    runs the per-type first-order gate on it with the off-site rows of
    :data:`_GATE_ROWS`, because the same asymmetry that makes the Fourier
    sign visible also makes ``+x`` differ from its reverse, so the two
    orbital placements of a declared row are genuinely different
    Hamiltonians here. Neither consumer diagonalises it: the gate needs
    only the mean-field functional at a supplied density, so this fixture
    is still never handed to ``SectorED``.

    ``EDFixture.build_h1`` places ``t[(a, b)]`` on ``(j+1, a) <- (j, b)``
    and ``conj(t[(b, a)])`` on the return hop, so ``h1`` is Hermitian for
    ANY pair ``t_01``, ``t_10`` -- the asymmetry below is legal, not a
    broken Hamiltonian, and ``_write_chain_inputs`` writes exactly those
    two amplitudes into ``transfer.dat``."""
    return edu.EDFixture(L=3, norb=2,
                         t={(0, 0): -1.0 + 0.3j, (1, 1): -0.7 - 0.2j,
                            (0, 1): -0.2 + 0.15j, (1, 0): -0.25 - 0.1j},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


def _fx_dense():
    """``L = 2``, two orbitals, COMPLEX hopping in every orbital channel:
    ``nmode = 8``, Fock dimension 256 -- small enough to diagonalise DENSELY
    over the whole Fock space, which is what makes it a reference for the
    sector engine.

    Complex, transpose-asymmetric hopping is deliberate: on a real symmetric
    band a sector-pairing or conjugation error in the Lehmann contraction can
    hide."""
    return edu.EDFixture(L=2, norb=2,
                         t={(0, 0): -1.0 + 0.3j, (1, 1): -0.7 - 0.2j,
                            (0, 1): -0.2 + 0.15j, (1, 0): -0.25 - 0.1j},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


def _dense_spectrum(fx, terms=()):
    """``(E, W, w, Z, C, CD, K)`` of the grand-canonical Hamiltonian
    ``K = H1 + H_int - mu N`` diagonalised on the FULL Fock space.

    Nothing here knows about particle-number sectors, which is the point:
    ``SectorED`` splits the spectrum by ``(N_up, N_dn)`` and pairs the
    blocks by hand, and this is the reference that pairing is checked
    against."""
    C = fx.annihilators()
    CD = [c.conj().T for c in C]
    h1 = fx.build_h1()
    K = edu.h_int_from_terms(fx, terms) if terms else np.zeros((fx.dim, fx.dim), dtype=complex)
    for p in range(fx.nmode):
        for q in range(fx.nmode):
            if h1[p, q] != 0:
                K = K + h1[p, q] * (CD[p] @ C[q])
    K = K - fx.mu * sum(CD[p] @ C[p] for p in range(fx.nmode))
    ev, W = np.linalg.eigh(K)
    E = ev - ev.min()
    w = np.exp(-fx.beta * E)
    return E, W, w, w.sum(), C, CD, K


def _dense_green_and_density(fx, terms, iws):
    """``(G, rho)`` from the dense Lehmann sums over the full Fock space:

        G_{pq}(i w) = (1/Z) sum_{m n} (w_m + w_n) <m|c_p|n> conj(<m|c_q|n>)
                                      / (i w - (E_n - E_m))
        rho_{pq}    = (1/Z) sum_n w_n <n| c^dag_p c_q |n>
    """
    E, W, w, Z, C, CD, _K = _dense_spectrum(fx, terms)
    cm = [W.conj().T @ C[p] @ W for p in range(fx.nmode)]
    boltz = w[:, None] + w[None, :]
    dE = E[None, :] - E[:, None]
    G = np.zeros((len(iws), fx.nmode, fx.nmode), dtype=complex)
    for p in range(fx.nmode):
        for q in range(fx.nmode):
            num = cm[p] * np.conj(cm[q]) * boltz
            for i, iw in enumerate(iws):
                G[i, p, q] = (num / (iw - dE)).sum()
    rho = np.zeros((fx.nmode, fx.nmode), dtype=complex)
    for p in range(fx.nmode):
        for q in range(fx.nmode):
            rho[p, q] = (np.diag(W.conj().T @ (CD[p] @ C[q]) @ W) * w).sum()
    return G / Z, rho / Z


def _dense_first_moment(fx, terms=()):
    """``M1_{pq} = <{[c_p, K], c^dag_q}>``, the first high-frequency moment of
    ``G``: ``G(i w) = I / (i w) + M1 / (i w)^2 + O(w^-3)``.

    Built from the operators themselves, not from any expansion. The sign
    convention (``[c_p, K]``, not ``[K, c_p]``) is not assumed either: at
    zero interaction this must be ``h1 - mu``, which
    :meth:`TestGreenFunction.test_high_frequency_moments` asserts."""
    E, W, w, Z, C, CD, K = _dense_spectrum(fx, terms)
    M1 = np.zeros((fx.nmode, fx.nmode), dtype=complex)
    for p in range(fx.nmode):
        comm = C[p] @ K - K @ C[p]
        for q in range(fx.nmode):
            A = comm @ CD[q] + CD[q] @ comm
            M1[p, q] = (np.diag(W.conj().T @ A @ W) * w).sum()
    return M1 / Z


# --------------------------------------------------------------------------
# Declared rows -> the exact Hamiltonian
# --------------------------------------------------------------------------

def _rows_u_v(u, v):
    """On-site ``U`` and nearest-neighbour ``V`` of the single-orbital chain,
    as the reader delivers them (both orientations of the bond)."""
    rows = {}
    if u != 0.0:
        rows["CoulombIntra"] = [(0, 0, 0, 1, 1, u, 0.0)]
    if v != 0.0:
        rows["CoulombInter"] = [(1, 0, 0, 1, 1, v, 0.0), (-1, 0, 0, 1, 1, v, 0.0)]
    return rows


def _rows_offsite(itype, v):
    """One off-site single-orbital bond of ``itype``, both orientations."""
    return {itype: [(1, 0, 0, 1, 1, v, 0.0), (-1, 0, 0, 1, 1, v, 0.0)]}


def _rows_interorbital(v, ratio=0.6):
    """The ORBITAL-ASYMMETRIC inter-orbital bond of the two-orbital chain:
    ``v_01(+x) = v`` and ``v_10(+x) = ratio * v`` with ``ratio != 1``, each
    with its reversed partner. Both orbital assignments of the SAME
    displacement are declared with DIFFERENT amplitudes, so nothing here is
    invariant under the orbital swap the ED adjudicates."""
    return {"CoulombInter": [(1, 0, 0, 1, 2, v, 0.0), (-1, 0, 0, 2, 1, v, 0.0),
                             (1, 0, 0, 2, 1, ratio * v, 0.0), (-1, 0, 0, 1, 2, ratio * v, 0.0)]}


def _rows_interorbital_one(v):
    """ONE inter-orbital bond of the two-orbital chain, both orientations:
    ``v_01(+x) = v`` alone. The orbital swap maps this declaration onto a
    DIFFERENT Hamiltonian (unlike :func:`_rows_interorbital` at ratio 1,
    where the two declared bonds exchange and the swap is a symmetry), so
    this is the shape the orientation check needs."""
    return {"CoulombInter": [(1, 0, 0, 1, 2, v, 0.0), (-1, 0, 0, 2, 1, v, 0.0)]}


def _rows_hund_offsite_v(j, v):
    """On-site ``Hund`` BETWEEN the two orbitals plus the asymmetric
    inter-orbital off-site ``V`` rows of :func:`_rows_interorbital`, both
    orientations of each.

    This is the fixture of the mixed on-site/off-site class: the kernel's
    ``A_on chibar B_v + A_v chibar B_on`` products, with an on-site vertex
    that is NOT a plain density (``Hund`` is same-spin only), against the
    off-site density vertex. Every such diagram is local by the weight rule
    of spec 2.4 -- which
    :meth:`TestG4.test_hund_times_offsite_v_chain` asserts directly."""
    rows = {}
    if j != 0.0:
        rows["Hund"] = [(0, 0, 0, 1, 2, j, 0.0), (0, 0, 0, 2, 1, j, 0.0)]
    if v != 0.0:
        rows.update(_rows_interorbital(v))
    return rows


def _ed_terms(fx, rows_by_type, mirrored_weight=0.5, swap_orbitals=False):
    """The quartic term list (``(p, q, r, s, coeff)`` for
    ``coeff c^dag_p c_q c^dag_r c_s``, the form ``SectorED`` consumes) of
    the Hamiltonian the DECLARED rows describe.

    A row ``(R, a, b, v)`` places orbital ``a`` in the ORIGINAL cell ``j``
    and orbital ``b`` in the cell displaced by ``R`` -- the DOCUMENTED
    reading of an interaction row (``docs/en`` UHFk interaction file page:
    "``[alpha]`` corresponds to the orbital alpha in the original cell, and
    ``[beta]`` ... in the cell displaced by r"), which the mean-field
    kernel implements since the orientation fix of issue #192. ``w`` is the
    mirrored-row weight.

    Writing ``j' = j + R``, the monomials are the ``docs/en`` definitions
    with ``i -> j`` (orbital ``a``) and ``j -> j'`` (orbital ``b``):

    ``CoulombIntra``
        ``v n_{j,a,up} n_{j,a,dn}`` (on-site, same orbital).
    density types
        ``w v sum_j sum_{s1 s2} sigma(s1, s2) n_{j,a,s1} n_{j',b,s2}`` with
        ``sigma`` the type's spin structure: all four pairings for
        ``CoulombInter``, only the same-spin ones with a MINUS for
        ``Hund``, the density-DIFFERENCE signs for ``Ising``.
    ``Exchange``
        ``w v sum_j sum_s c^dag_{j,a,s} c_{j',b,s} c^dag_{j',b,-s} c_{j,a,-s}``
        -- the doc's ``up/dn`` monomial plus its ``dn/up`` partner --
        written here in the ANTICOMMUTED form
        ``- w v sum_j sum_s (c^dag_{j,a,s} c_{j,a,-s})(c^dag_{j',b,-s} c_{j',b,s})``,
        i.e. as a product of the two ON-SITE spin-flip bilinears. The four
        modes are distinct (the two sites differ, and so do the spins), so
        the reordering costs two transpositions and one overall sign and is
        the SAME operator --
        :meth:`TestChainHamiltonian.test_the_exchange_monomial_is_the_documented_one`
        multiplies both writings out and compares them at zero tolerance.
        The grouping matters only for the SPLIT: the Wick functional's
        direct lines are the ones that pair the operators as they are
        written, and UHFk's ``flag_fock = false`` keeps exactly the
        on-site-grouped pair for this type (its ``hh0``/``hh3`` block reads
        the equal-site density ``gbb``). Writing the monomial in UHFk's own
        grouping is what makes :func:`_hf_sigma`'s ``include_fock`` switch
        mean the same thing as the kernel's, for every type at once.
    ``PairLift``
        ``w v sum_j c^dag_{j,a,up} c_{j,a,dn} c^dag_{j',b,up} c_{j',b,dn}
        + h.c.`` -- the ``+ h.c.`` of the definition, per declared row.
    ``PairHop``
        ``v sum_j c^dag_{j,a,up} c_{j',b,up} c^dag_{j,a,dn} c_{j',b,dn}``
        with NO ``+ h.c.`` and no mirrored halving: ``PairHop`` is absent
        from :data:`_MIRRORED_TYPES` because its mirrored row
        ``(-R, b, a, conj v)`` IS the Hermitian conjugate of this one, so a
        Hermitian-closed declaration already carries the ``h.c.``.

    ``swap_orbitals`` builds the OPPOSITE placement (``a`` in the displaced
    cell) and exists only so the orientation can be checked two-sidedly.
    """
    terms = []
    for itype, rows in rows_by_type.items():
        w = mirrored_weight if itype in _MIRRORED_TYPES else 1.0
        for (rx, ry, rz, a1, b1, vr, vi) in rows:
            v = w * complex(vr, vi)
            if v == 0.0:
                continue
            a, b, R = a1 - 1, b1 - 1, rx
            if swap_orbitals:
                a, b = b, a
            if itype == "CoulombIntra":
                if (rx, ry, rz) != (0, 0, 0) or a != b:
                    raise ValueError("CoulombIntra must be an on-site same-orbital row")
                for j in range(fx.L):
                    terms.append((fx.mode(j, a, UP), fx.mode(j, a, UP),
                                  fx.mode(j, a, DN), fx.mode(j, a, DN), v))
                continue
            sigma = None
            if itype == "CoulombInter":
                def sigma(s1, s2, v=v):
                    return v
            elif itype == "Hund":
                def sigma(s1, s2, v=v):
                    return -v if s1 == s2 else 0.0
            elif itype == "Ising":
                def sigma(s1, s2, v=v):
                    return v if s1 == s2 else -v
            elif itype not in ("Exchange", "PairLift", "PairHop"):
                raise ValueError("chain ED: unsupported interaction type {!r}".format(itype))
            for j in range(fx.L):
                jp = (j + R) % fx.L
                ma = [fx.mode(j, a, s) for s in (UP, DN)]     # orbital a, ORIGINAL cell
                mb = [fx.mode(jp, b, s) for s in (UP, DN)]    # orbital b, cell displaced by R
                if sigma is not None:
                    for s1 in range(2):
                        for s2 in range(2):
                            c = sigma(s1, s2)
                            if c == 0.0:
                                continue
                            terms.append((ma[s1], ma[s1], mb[s2], mb[s2], c))
                elif itype == "Exchange":
                    for s in (UP, DN):
                        # the documented monomial, anticommuted into the
                        # product of the two ON-SITE spin-flip bilinears
                        terms.append((ma[s], ma[1 - s], mb[1 - s], mb[s], -v))
                elif itype == "PairLift":
                    terms.append((ma[UP], ma[DN], mb[UP], mb[DN], v))
                    terms.append((mb[DN], mb[UP], ma[DN], ma[UP], np.conj(v)))
                else:                                          # PairHop
                    terms.append((ma[UP], mb[UP], ma[DN], mb[DN], v))
    return terms


def _hf_sigma(fx, terms, rho, include_fock=True):
    """First-order (Hartree-Fock) self-energy of ``terms`` at the density
    ``rho[p, q] = <c^dag_p c_q>``, in mode space.

    The Wick engine is ``ed_oracle_util.hf_h1_from_terms``'s ``add()``
    (``E = c [rho_pq rho_rs + rho_ps (delta_qr - rho_rq)]`` and
    ``H_MF = sum dE/d rho_xy c^dag_x c_y``), evaluated at an ARBITRARY
    density rather than the free one -- which is the only difference, and
    the reason the gates need it: the subtraction has to be made at the
    EXACT density. The engine itself is Task 10's
    :func:`tests.test_flex_second_order_ed_onsite.wick_hf_sigma`, shared
    rather than copied, so both exact-diagonalization gates subtract the
    same functional.

    ``include_fock=False`` keeps only the Hartree (direct) contractions, so
    the mean-field gate can confront the production kernel in BOTH of its
    Fock settings.

    This IS the functional the gates subtract; UHFk's own mean-field kernel
    is compared against it separately, as a convention pin, by
    :meth:`TestChainHamiltonian.test_ed_hamiltonian_is_the_production_mean_field`."""
    return wick_hf_sigma(fx.nmode, terms, rho, include_fock=include_fock)


def _production_hf_sigma_k(fx, rows_by_type, rho, include_fock=True):
    """UHFk's OWN mean-field Hamiltonian for the declared rows at the
    density ``rho``, as ``Sigma_HF(k)[(s, a), (t, b)]``.

    ``hwave.solver.hartree_fock.accumulate_hf`` is ``UHFk._make_ham``'s
    interaction loop as a pure function -- the definition of the
    Hamiltonian the solver implements, Hartree AND Fock, for on-site and
    off-site rows alike. It is what makes the first-order gate below
    two-sided in both the mirrored-row weight and the orbital placement,
    exactly as Task 10's gate is two-sided through ``compile_onsite``.
    ``include_fock`` is forwarded to the kernel, so the gate can pin the
    Hartree-only setting (UHFk's ``flag_fock = false``) as well."""
    from hwave.solver.hartree_fock import accumulate_hf, build_interaction_tables
    norb, L = fx.norb, fx.L
    shape = (L, 1, 1)
    param_ham = {}
    for itype, rows in rows_by_type.items():
        tbl = {}
        for (rx, ry, rz, a1, b1, vr, vi) in rows:
            key = ((rx % L, ry, rz), (a1 - 1, b1 - 1))
            tbl[key] = tbl.get(key, 0.0) + complex(vr, vi)
        param_ham[itype] = tbl
    tabs = build_interaction_tables(param_ham, norb, shape)
    rho_so_r = np.zeros((L, 2, norb, 2, norb), dtype=complex)
    for r in range(L):
        for s in range(2):
            for a in range(norb):
                for t in range(2):
                    for b in range(norb):
                        rho_so_r[r, s, a, t, b] = rho[fx.mode(0, a, s), fx.mode(r, b, t)]
    out = np.zeros((L, 2 * norb, 2 * norb), dtype=complex)
    accumulate_hf(out, rho_so_r, tabs.inter_table, tabs.spin_table, shape,
                  include_fock=include_fock)
    return out


# --------------------------------------------------------------------------
# Mode space -> the solver's (k, orbital) frame
# --------------------------------------------------------------------------

def _k_phases(fx, sign=1.0):
    """``phase[kx, R] = e^{sign i k R}`` on ``k = 2 pi kx / L``. The ONE
    place the Fourier sign of this module's real-space-to-k bridge is
    written; ``sign = -1`` builds the reversed-k candidate that
    :meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`
    has to reject."""
    R = np.arange(fx.L)
    return np.exp(sign * 2j * np.pi * np.outer(R, R) / fx.L)


def _to_k(fx, M, sign=1.0):
    """The UP-spin block of a translation-invariant mode-space matrix, as
    ``A(k)_{ab} = sum_R A[(R, a), (0, b)] e^{i k R}``, ``k = 2 pi kx / L``.

    This is the convention ``EDFixture.build_h1`` and the solver's
    ``epsilon_k[a,b] += t_ab e^{+ikR}`` share (a ``transfer.dat`` row at
    ``R`` carries ``(j + R, a) <- (j, b)``), and
    :meth:`TestChainHamiltonian.test_bare_green_matches_the_production_solver`
    pins the whole chain of it end to end -- the SIGN of the exponent and
    the orbital index order specifically by
    :meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`,
    which needs a band that is neither inversion- nor transpose-symmetric to
    see either."""
    norb, L = fx.norb, fx.L
    phase = _k_phases(fx, sign)
    out = np.zeros((L, norb, norb), dtype=complex)
    for kx in range(L):
        for R in range(L):
            for a in range(norb):
                for b in range(norb):
                    out[kx, a, b] += phase[kx, R] * M[fx.mode(R, a, UP), fx.mode(0, b, UP)]
    return out


def _to_k_so(fx, M):
    """:func:`_to_k` on the FULL spin-orbital block: ``A(k)[(s a), (t b)]``
    in the ``s * norb + a`` order of ``hartree_fock.accumulate_hf``'s
    output, so a mode-space matrix can be compared with UHFk's own
    mean-field kernel on every spin block rather than only the up one."""
    norb, L = fx.norb, fx.L
    phase = _k_phases(fx)
    out = np.zeros((L, 2 * norb, 2 * norb), dtype=complex)
    for kx in range(L):
        for R in range(L):
            for s in range(2):
                for a in range(norb):
                    for t in range(2):
                        for b in range(norb):
                            out[kx, s * norb + a, t * norb + b] += (
                                phase[kx, R] * M[fx.mode(R, a, s), fx.mode(0, b, t)])
    return out


def _random_translation_invariant_density(fx, n=3, seed=20260915,
                                          spin_coherent=False):
    """``n`` deterministic random Hermitian, translation-invariant densities
    in mode space; spin-block-DIAGONAL by default.

    Translation invariance is required by UHFk's kernel (it reads one
    reference site's row), Hermiticity by the functional; the spin blocks
    are independent, so the densities are magnetic as well as
    orbital-coherent -- everything the free density of the pin below is
    not.

    ``spin_coherent=True`` populates the off-diagonal SPIN blocks
    ``<c^dag_up c_dn>`` too. That is not a refinement but a requirement for
    three of the types: ``PairLift``'s whole mean field, and the spin-flip
    contractions of ``Exchange`` and ``PairHop``, are proportional to a
    spin coherence and vanish identically on any spin-block-diagonal
    density -- so a gate run on such a density would be comparing zero with
    zero for them. The spin-diagonal default is kept so that the existing
    callers, which pin the density types, are unchanged."""
    rng = np.random.default_rng(seed)
    norb, L = fx.norb, fx.L
    nd = 2 * norb
    out = []
    for _ in range(n):
        A = rng.normal(size=(L, nd, nd)) + 1j * rng.normal(size=(L, nd, nd))
        if not spin_coherent:
            blk = A.copy()
            A = np.zeros_like(A)
            A[:, :norb, :norb] = blk[:, :norb, :norb]
            A[:, norb:, norb:] = blk[:, norb:, norb:]
        g = np.zeros_like(A)
        for R in range(L):
            # g[R] = conj(g[-R]^T) makes rho Hermitian
            g[R] = 0.5 * (A[R] + np.conj(A[(-R) % L].T))
        g *= 0.15
        g[0] += 0.5 * np.eye(nd)
        rho = np.zeros((fx.nmode, fx.nmode), dtype=complex)
        for j in range(L):
            for j2 in range(L):
                for s in range(2):
                    for a in range(norb):
                        for t in range(2):
                            for b in range(norb):
                                rho[fx.mode(j, a, s), fx.mode(j2, b, t)] = (
                                    g[(j - j2) % L, s * norb + a, t * norb + b])
        out.append(rho)
    return out


def _g0_k(fx, iws, sign=1.0):
    """``G0(k, i w) = [(i w + mu) - eps(k)]^{-1}`` of the free chain."""
    ek = _to_k(fx, fx.build_h1(), sign)
    eye = np.eye(fx.norb)
    return np.array([[np.linalg.inv((iw + fx.mu) * eye - ek[kx]) for kx in range(fx.L)]
                     for iw in iws])


def _green_ed_k(ed, fx, iws, sign=1.0):
    """The interacting ``G(k, i w)`` of the up spin, from ``SectorED.green``.

    Only the ``norb`` columns at site 0 are asked for: the fixture is
    translation invariant, so ``G[(R, a), (0, b)]`` already carries every
    displacement, and the restriction turns an ``nmode^2`` Lehmann
    contraction into an ``nmode * norb`` one."""
    L, norb = fx.L, fx.norb
    rows = [fx.mode(j, a, UP) for j in range(L) for a in range(norb)]
    cols = [fx.mode(0, b, UP) for b in range(norb)]
    gm = ed.green(iws, rows=rows, cols=cols)
    phase = _k_phases(fx, sign)
    gk = np.zeros((len(iws), L, norb, norb), dtype=complex)
    for kx in range(L):
        for R in range(L):
            gk[:, kx] += phase[kx, R] * gm[:, R * norb:(R + 1) * norb, :]
    return gk


def _sigma_ed(fx, rows_by_type, iws, mirrored_weight=0.5, hf_source="wick"):
    """``(Sigma_ED, Sigma_HF)`` on ``iws``: the exact self-energy of the
    chain from the Dyson inversion ``G0^{-1} - G^{-1}``, and a first-order
    functional evaluated with the EXACT density.

    ``hf_source`` selects WHOSE first order is subtracted:

    ``"wick"`` (the default, and what the second-order gates use)
        the TEST SIDE's own :func:`_hf_sigma` on the ED term list, carried
        to ``(k, orbital)`` by :func:`_to_k`. The remainder the gates then
        compare is production-free.
    ``"production"``
        UHFk's own mean-field kernel (:func:`_production_hf_sigma_k`). Used
        only by :meth:`TestFirstOrder.test_first_order_remainder_vanishes`,
        where confronting the production kernel with this Hamiltonian is
        the point -- and the reason that gate can be two-sided in the
        mirrored-row weight at all (the ``"wick"`` functional tracks
        whatever weight it is handed, so it cannot discriminate one).
    """
    ed = edu.SectorED(fx, terms=_ed_terms(fx, rows_by_type, mirrored_weight))
    sigma = np.linalg.inv(_g0_k(fx, iws)) - np.linalg.inv(_green_ed_k(ed, fx, iws))
    rho = ed.density_matrix()
    if hf_source == "production":
        return sigma, _production_hf_sigma_k(fx, rows_by_type, rho)[:, :fx.norb, :fx.norb]
    if hf_source != "wick":
        raise ValueError("hf_source must be 'wick' or 'production'")
    terms = _ed_terms(fx, rows_by_type, mirrored_weight)
    return sigma, _to_k(fx, _hf_sigma(fx, terms, rho))


def _sigma_fluct_ed(fx, rows_by_type, iws, mirrored_weight=0.5, hf_source="wick"):
    """``Sigma_ED - Sigma_HF[rho_ED]`` -- the exact fluctuation remainder,
    whose O(v^2) coefficient is the bare second-order skeleton. See
    :func:`_sigma_ed` for ``hf_source``."""
    sigma, hf = _sigma_ed(fx, rows_by_type, iws, mirrored_weight, hf_source)
    return sigma - hf[None]


# --------------------------------------------------------------------------
# The production side
# --------------------------------------------------------------------------

def _write_chain_inputs(d, fx, rows_by_type):
    """``geom.dat`` / ``transfer.dat`` / one file per interaction type
    describing ``fx``'s chain on ``CellShape = [L, 1, 1]``: a one-cell
    geometry with ``norb`` orbitals, the hopping as ``(+1, a, b, t_ab)``
    and its Hermitian partner ``(-1, a, b, conj(t_ba))``, and the on-site
    energies at ``R = 0``."""
    norb = fx.norb
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n{}\n".format(norb))
        for _ in range(norb):
            f.write("0.0 0.0 0.0\n")
    rows = []
    for (a, b), tv in fx.t.items():
        back = np.conj(fx.t[(b, a)])
        rows.append((1, 0, 0, a + 1, b + 1, tv.real, tv.imag))
        rows.append((-1, 0, 0, a + 1, b + 1, back.real, back.imag))
    for o, e in enumerate(fx.eps):
        rows.append((0, 0, 0, o + 1, o + 1, float(np.real(e)), 0.0))
    _write_wan(os.path.join(d, "transfer.dat"), "Transfer", norb, rows)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for itype, irows in rows_by_type.items():
        _write_wan(os.path.join(d, itype.lower() + ".dat"), itype, norb, irows)
        idict[itype] = itype.lower() + ".dat"
    return idict


def _chain_solver(d, fx, rows_by_type, gate, nmat):
    """A FLEX general solver on ``fx``'s chain with ``flex_second_order =
    "local"``, at the FIXED ``mu = fx.mu`` (no SCF, no density update), and
    with the bond gate switched on through its own two parameters when
    ``gate`` -- so the Phase B preflight and vertex preparation run exactly
    as in a real solve."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    idict = _write_chain_inputs(d, fx, rows_by_type)
    reader = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": fx.T, "mu": fx.mu, "CellShape": [fx.L, 1, 1], "SubShape": [1, 1, 1],
           "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
           "flex_second_order": "local"}
    if gate:
        par["flex_hartree_fock"] = True
        par["longitudinal_bond_channels"] = True
    return flex_mod.FLEX(reader.get_param("ham"), {},
                         {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
                          "calc_scheme": "general"})


def _production_sigma(fx, rows_by_type, gate, nmat):
    """Production ``Sigma_fluct`` of ONE map on the BARE Green function at
    the fixed ``mu``, and that Green function.

    ``gate=False`` is the standalone general path
    (``_flex_compute_veff_general`` -> ``_calc_self_energy_general``);
    ``gate=True`` is Phase B's bond-resolved path
    (``dress_and_build_w`` -> ``calc_self_energy_bond``), driven exactly as
    ``tests/test_flex_second_order_sopt.py::_one_map_sigma`` drives it."""
    from hwave.solver import flex_bond
    beta, L, norb = fx.beta, fx.L, fx.norb
    shape = (L, 1, 1)
    with tempfile.TemporaryDirectory() as d, _quiet():
        s = _chain_solver(d, fx, rows_by_type, gate, nmat)
        if gate:
            green_info = {}
            s._phase_b_reset(green_info)
            s._phase_b_preflight(green_info)
        s._calc_epsilon_k({})
        G = s._calc_dressed_green(beta, fx.mu,
                                  np.zeros((1, nmat, L, norb, norb), complex))
        if not gate:
            chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), beta)[0]
            _, v_eff, _, _ = s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q)
            return s._calc_self_energy_general(G, v_eff, beta)[0], G
        nd = norb * norb
        nb = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, L, nb * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, np.zeros_like(G), beta, s._bond_view, shape, 1)
            flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on,
                                        C_on=s._bond_C_on, nb=nmat, output_full=False,
                                        nmat=nmat, nvol=L, nd=nd, spatial_shape=shape,
                                        factors=s._second_order_factors, second_order="local")
            return flex_bond.calc_self_energy_bond(store, G, beta, s._bond_view,
                                                   shape, norb, 1)[0], G


class _Maps(object):
    """Memoised ED / production maps of one fixture: the pure-type points of
    a cross stencil are the same maps the pure-type stencils already asked
    for, and the gate-on and gate-off runs of a type share nothing but the
    ED side."""

    def __init__(self, fx, rows_of, nmat=_NMAT):
        self.fx, self.rows_of, self.nmat = fx, rows_of, nmat
        self.idx = _inner()
        self.iws = 1j * (2 * self.idx + 1 - nmat) * np.pi / fx.beta
        self._ed, self._pr = {}, {}
        self._bare = None

    def _zero(self):
        return np.zeros((len(self.idx), self.fx.L, self.fx.norb, self.fx.norb), complex)

    def ed(self, *args):
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        key = args
        if key not in self._ed:
            self._ed[key] = _sigma_fluct_ed(self.fx, rows, self.iws)
        return self._ed[key]

    def prod(self, *args, **kw):
        gate = kw.pop("gate", False)
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        key = (gate,) + args
        if key not in self._pr:
            self._pr[key] = _production_sigma(self.fx, rows, gate, self.nmat)[0][self.idx]
        return self._pr[key]

    def bare_green(self, rows):
        """The bare Green function the oracle's dropped class is evaluated
        on -- the SAME propagator both sides expand around. It is
        coupling-independent (no SCF, fixed mu), so it is built once from
        whichever map asks for it first and reused."""
        if self._bare is None:
            self._bare = _production_sigma(self.fx, rows, False, self.nmat)[1]
        return self._bare

    def dropped(self, *args):
        """The oracle's off/off UNCROSSED class at the given couplings. The
        oracle's skeleton is EXACTLY quadratic in the couplings, so a
        coefficient needs no stencil: the value at unit coupling IS it."""
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        recs = oracle_records(rows, self.fx.norb, (self.fx.L, 1, 1))
        return oracle_sigma2(self.bare_green(rows), self.fx.beta, recs, self.fx.norb,
                             "dropped", (self.fx.L, 1, 1))[0][self.idx]


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

class TestGreenFunction(unittest.TestCase):
    """``SectorED.green`` / ``density_matrix``: the single-particle
    observables the chain gate is built on, pinned where they are known in
    closed form or by an independent dense evaluation."""

    def test_free_chain_green_equals_the_one_body_resolvent(self):
        """With no interaction the cross-sector Lehmann sum must collapse to
        ``[i w - (h1 - mu)]^{-1}`` on every mode pair -- which pins the
        Boltzmann weighting, the grand-canonical denominators (``mu``
        enters them, the sectors differ by one particle) and the global
        partition function all at once."""
        fx = _fx_chain()
        iws = 1j * (2 * np.arange(8) + 1 - 8) * np.pi / fx.beta
        G = edu.green_function(fx, iws)
        h1 = fx.build_h1()
        eye = np.eye(fx.nmode)
        for i, iw in enumerate(iws):
            np.testing.assert_allclose(
                G[i], np.linalg.inv(iw * eye - (h1 - fx.mu * eye)), atol=1e-12)

    def test_single_site_interacting_green_equals_the_dense_lehmann_sum(self):
        """One site, one orbital, ``U n_up n_dn``: the sector engine against
        the DENSE Lehmann sum over the full Fock space
        (``EDFixture.annihilators``), which knows nothing about sectors.
        The free pin above cannot see a sector-pairing error that the
        interaction makes matter; this one can."""
        fx = edu.EDFixture(L=1, norb=1, t={(0, 0): 0.0}, eps=(0.3,), T=0.5, mu=0.2)
        terms = edu.canonical_density_terms(fx, [(0, 0, 0, 0.7)])
        iws = 1j * (2 * np.arange(6) + 1 - 6) * np.pi / fx.beta
        G = edu.green_function(fx, iws, terms)

        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        h1 = fx.build_h1()
        H = edu.h_int_from_terms(fx, terms)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                if h1[p, q] != 0:
                    H = H + h1[p, q] * (CD[p] @ C[q])
        N = sum(CD[p] @ C[p] for p in range(fx.nmode))
        E, V = np.linalg.eigh(H - fx.mu * N)
        E = E - E.min()
        w = np.exp(-fx.beta * E)
        cm = [V.conj().T @ C[p] @ V for p in range(fx.nmode)]
        boltz = w[:, None] + w[None, :]
        dE = E[None, :] - E[:, None]
        ref = np.zeros_like(G)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                num = cm[p] * np.conj(cm[q]) * boltz
                for i, iw in enumerate(iws):
                    ref[i, p, q] = (num / (iw - dE)).sum()
        ref /= w.sum()
        self.assertGreater(np.abs(ref).max(), 1e-3)            # anti-vacuity
        np.testing.assert_allclose(G, ref, atol=1e-12)

    def test_multiorbital_interacting_green_and_density_equal_the_dense_sum(self):
        """``SectorED.green`` and ``SectorED.density_matrix`` on an
        INTERACTING MULTI-ORBITAL fixture, against a dense Lehmann
        evaluation over the whole 256-dimensional Fock space.

        The two pins above cover the free multi-orbital case and the
        interacting SINGLE-orbital one; neither can see a sector pairing that
        goes wrong only when the interaction mixes orbitals -- which is the
        regime the chain gates actually run in. This fixture has complex,
        transpose-asymmetric hopping in every orbital channel and both
        on-site and inter-site inter-ORBITAL density terms."""
        fx = _fx_dense()
        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.9, 0.0), (0, 0, 0, 2, 2, 0.5, 0.0)],
                "CoulombInter": [(0, 0, 0, 1, 2, 0.4, 0.0), (0, 0, 0, 2, 1, 0.4, 0.0),
                                 (1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]}
        terms = _ed_terms(fx, rows)
        self.assertEqual(fx.dim, 256)
        iws = 1j * (2 * np.arange(8) + 1 - 8) * np.pi / fx.beta
        G_ref, rho_ref = _dense_green_and_density(fx, terms, iws)
        ed = edu.SectorED(fx, terms=terms)
        self.assertGreater(np.abs(G_ref).max(), 1e-3)            # anti-vacuity
        self.assertGreater(np.abs(rho_ref).max(), 1e-3)
        # the fixture really is orbital-mixing and complex, so the
        # comparison is not blind to either
        self.assertGreater(np.abs(np.asarray(G_ref).imag).max(), 1e-3)
        off = [(p, q) for p in range(fx.nmode) for q in range(fx.nmode)
               if p % 2 == q % 2 and (p // 2) != (q // 2)]
        self.assertGreater(max(abs(rho_ref[p, q]) for p, q in off), 1e-3)
        np.testing.assert_allclose(ed.green(iws), G_ref, rtol=0, atol=1e-12)
        np.testing.assert_allclose(ed.density_matrix(), rho_ref, rtol=0, atol=1e-12)

    def test_high_frequency_moments(self):
        """``i w G(i w) -> I`` and the first moment
        ``(i w)^2 (G - I / (i w)) -> M1 = <{[c_p, K], c^dag_q}>``.

        These are the two statements the chain gate's Dyson inversion
        ``G0^{-1} - G^{-1}`` silently relies on: if ``G`` did not carry the
        free ``1/(i w)`` head with the right normalisation, the subtraction
        would not remove the free part, and every coefficient downstream
        would be measured against the wrong reference. Both are checked as
        LIMITS -- the deviation has to shrink as the grid is extended --
        rather than at one frequency, where a wrong constant could hide.

        Measured on this fixture: ``max|i w G - I|`` is 4.1e-2, 1.0e-2 and
        2.5e-3 at the largest frequency of the nmat = 32, 128 and 512 grids
        (falling as ``|M1| / |w|``, ``|M1| = 2.0``), and the first-moment
        residual is 6.2e-3 -- 0.3% of ``|M1|`` -- at the largest frequency of
        the 512 grid, against 1.36 at the middle of it."""
        fx = _fx_dense()
        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.9, 0.0)],
                "CoulombInter": [(0, 0, 0, 1, 2, 0.4, 0.0), (0, 0, 0, 2, 1, 0.4, 0.0)]}
        terms = _ed_terms(fx, rows)
        ed = edu.SectorED(fx, terms=terms)
        eye = np.eye(fx.nmode)

        # the sign convention of M1, pinned where it is known in closed form
        free = _dense_first_moment(fx)
        np.testing.assert_allclose(free, fx.build_h1() - fx.mu * eye, rtol=0, atol=1e-10)

        head = []
        for nmat in (32, 128, 512):
            iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
            g = ed.green(iws)
            dev = np.abs(iws[:, None, None] * g - eye[None]).max(axis=(1, 2))
            head.append((abs(iws[-1]), dev[-1], g, iws))
        # the head is a LIMIT: extending the grid must shrink the deviation
        for (w0, d0, _g0, _i0), (w1, d1, _g1, _i1) in zip(head, head[1:]):
            self.assertGreater(w1, w0)
            self.assertLess(d1, d0)
        wmax, dmax, g, iws = head[-1]
        self.assertLess(dmax, 5.0e-3,
                        "i w G(i w) does not approach the identity: {:.3e} at |w| = {:.1f}"
                        .format(dmax, wmax))

        M1 = _dense_first_moment(fx, terms)
        scale = np.abs(M1).max()
        self.assertGreater(scale, 1e-3)                          # anti-vacuity
        est = (iws[:, None, None] ** 2) * (g - eye[None] / iws[:, None, None])
        far = np.abs(est[-1] - M1).max()
        near = np.abs(est[len(iws) // 2] - M1).max()
        self.assertLess(far, 1.0e-2 * scale,
                        "the first moment does not converge: {:.3e} of |M1| = {:.3e}"
                        .format(far / scale, scale))
        self.assertLess(far, 0.05 * near)       # and it really is a limit

    def test_free_density_matrix_equals_the_fermi_occupation(self):
        """``SectorED.density_matrix`` at zero coupling is the free
        ``<c^dag_p c_q> = [W f(ev) W^dagger]^T`` of ``h1 - mu``, the same
        expression ``hf_h1_from_terms`` Wick-contracts against."""
        fx = _fx_orbital()
        h1 = fx.build_h1()
        ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
        f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
        np.testing.assert_allclose(edu.SectorED(fx).density_matrix(),
                                   ((W * f) @ W.conj().T).T, atol=1e-12)


class TestChainHamiltonian(unittest.TestCase):
    """What the chain ED must be for the comparison to mean anything: the
    solver's own band structure, and the solver's own Hamiltonian."""

    def test_bare_green_matches_the_production_solver(self):
        """The bare ``G(k, i w)`` the solver builds from the chain's
        ``transfer.dat`` equals ``_g0_k``'s ``[(i w + mu) - eps(k)]^{-1}``,
        on BOTH fixtures.

        On these two fixtures the band is inversion- and transpose-symmetric
        (real hopping), so what this pin actually fixes is the
        ``(j + R, a) <- (j, b)`` reading of a transfer row, the ``kx``
        ordering and the site-major ``j * norb + a`` row layout -- checked
        against production rather than re-derived. The exponent SIGN and
        the orbital index order cannot be seen here and are pinned
        separately, on a complex-hopping band, by
        :meth:`test_fourier_sign_and_orbital_order_are_pinned`."""
        nmat = 32
        for fx in (_fx_chain(), _fx_orbital()):
            rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 1e-9, 0.0)]}
            _sig, G = _production_sigma(fx, rows, False, nmat)
            iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
            g0 = _g0_k(fx, iws)
            with self.subTest(L=fx.L, norb=fx.norb):
                self.assertGreater(np.abs(g0).max(), 1e-3)      # anti-vacuity
                np.testing.assert_allclose(G[0], g0, atol=1e-12)
                # and the ED engine's own free Green function is that same
                # object, so the Dyson inversion below subtracts like with like
                np.testing.assert_allclose(_green_ed_k(edu.SectorED(fx), fx, iws),
                                           g0, atol=1e-12)

    def test_fourier_sign_and_orbital_order_are_pinned(self):
        """On a band that is neither inversion- nor transpose-symmetric, the
        solver's own bare ``G(k, i w)`` selects ``_to_k``'s exponent SIGN
        and its orbital index order -- both of which the gate's two
        fixtures are blind to.

        Why this is load-bearing rather than decorative: ``_fx_chain`` and
        ``_fx_orbital`` have real hopping, so ``eps(k) = eps(-k)`` to
        round-off; with the exponent sign flipped every comparison in this
        module would still pass, and
        :meth:`test_ed_hamiltonian_is_the_production_mean_field` would then
        select the OPPOSITE orbital placement for a declared off-site row,
        inverting every downstream verdict with no pin going red. The same
        goes for the orbital index order on a real-symmetric ``eps(k)``.

        Three-sided: the ``e^{+ikR}`` convention matches production at
        1e-12; the reversed-k candidate (``sign = -1``, i.e. the array read
        at ``-kx``) and the orbital-transposed one both deviate by more
        than 1e-2. The two asymmetries are asserted to be O(1) first, so
        neither rejection can be vacuous."""
        nmat = 32
        fx = _fx_complex()
        iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
        eps = _to_k(fx, fx.build_h1())
        rev_idx = (-np.arange(fx.L)) % fx.L
        self.assertGreater(np.abs(eps - eps[rev_idx]).max(), 0.1,
                           "the fixture's band is inversion-symmetric: the Fourier "
                           "sign cannot be seen on it")
        self.assertGreater(np.abs(eps - np.swapaxes(eps, -1, -2)).max(), 0.1,
                           "the fixture's band is transpose-symmetric: the orbital "
                           "index order cannot be seen on it")

        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 1e-9, 0.0)]}
        _sig, G = _production_sigma(fx, rows, False, nmat)
        good = _g0_k(fx, iws)
        self.assertGreater(np.abs(good).max(), 1e-3)             # anti-vacuity
        np.testing.assert_allclose(G[0], good, atol=1e-12)
        # and the ED engine's own free Green function, through the same
        # transform, is that same object
        np.testing.assert_allclose(_green_ed_k(edu.SectorED(fx), fx, iws),
                                   good, atol=1e-12)
        reversed_k = _g0_k(fx, iws, sign=-1.0)
        self.assertGreater(np.abs(reversed_k - G[0]).max(), 1e-2,
                           "the reversed-k candidate also matches production: the "
                           "Fourier sign of _to_k is not pinned")
        transposed = np.swapaxes(good, -1, -2)
        self.assertGreater(np.abs(transposed - G[0]).max(), 1e-2,
                           "the orbital-transposed candidate also matches production: "
                           "the orbital index order of _to_k is not pinned")

    def test_ed_hamiltonian_is_the_production_mean_field(self):
        """The first-order functional of :func:`_ed_terms` equals UHFk's own
        mean-field kernel (``hartree_fock.accumulate_hf``, Hartree AND
        Fock) at the same density, for every fixture's declared rows -- and
        the two variants that would be wrong do NOT.

        This is what fixes the chain Hamiltonian, and it is two-sided in
        both conventions at once:

        * the ORBITAL PLACEMENT of an off-site row -- now the DOCUMENTED
          one (orbital ``a`` in the original cell). Measured deviation from
          the production mean field with the opposite placement: 0 on the
          single-orbital chain (where the two are the same Hamiltonian),
          5.434e-2 of the mean field on the single inter-orbital bond and
          1.359e-2 on the asymmetric pair, while the documented placement
          matches at 2.1e-16 / 3.9e-16 (every case here is at or below
          2.8e-15). This is the convention the campaign's second-order
          comparison rests on, and it is measured here rather than assumed.
        * the MIRRORED-ROW weight: at the full weight the off-site part of
          the mean field is doubled -- 1.000 of the mean field's own size
          on each of the four off-site cases.

        The free density used here is paramagnetic and spin-block-diagonal,
        which is enough for the density types but leaves ``Exchange``,
        ``PairLift`` and ``PairHop`` with no mean field at all;
        :meth:`test_first_order_gate_every_type_documented_reading` carries
        those, on spin-coherent densities.
        """
        cases = [("chain U", _fx_chain(), _rows_u_v(0.3, 0.0), False),
                 ("chain V", _fx_chain(), _rows_u_v(0.0, 0.3), False),
                 ("chain J", _fx_chain(), _rows_offsite("Hund", 0.3), False),
                 ("chain I", _fx_chain(), _rows_offsite("Ising", 0.3), False),
                 ("orbital single bond", _fx_orbital(), _rows_interorbital_one(0.3), True),
                 ("orbital asymmetric", _fx_orbital(), _rows_interorbital(0.3), True)]
        for (name, fx, rows, orbital_sensitive) in cases:
            h1 = fx.build_h1()
            ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
            f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
            rho = ((W * f) @ W.conj().T).T
            ref = _production_hf_sigma_k(fx, rows, rho)
            scale = np.abs(ref).max()
            with self.subTest(case=name):
                self.assertGreater(scale, 1e-3)                  # anti-vacuity
                # the production mean field is spin-diagonal and paramagnetic
                # on this density, so its up block is the whole story
                norb = fx.norb
                np.testing.assert_allclose(ref[:, :norb, norb:], 0.0, atol=1e-12)
                np.testing.assert_allclose(ref[:, norb:, :norb], 0.0, atol=1e-12)
                np.testing.assert_allclose(ref[:, norb:, norb:], ref[:, :norb, :norb],
                                           atol=1e-12)
                up = ref[:, :norb, :norb]
                mine = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows), rho))
                np.testing.assert_allclose(mine, up, atol=1e-11 * scale)
                swapped = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows, swap_orbitals=True), rho))
                if orbital_sensitive:
                    self.assertGreater(np.abs(swapped - up).max(), 1e-3 * scale,
                                       "{}: the orbital placement of the off-site row is not "
                                       "pinned -- both readings match the production mean "
                                       "field".format(name))
                else:
                    np.testing.assert_allclose(swapped, up, atol=1e-11 * scale)
                full = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows, mirrored_weight=1.0), rho))
                if any(t in _MIRRORED_TYPES for t in rows):
                    self.assertGreater(np.abs(full - up).max(), 1e-3 * scale,
                                       "{}: the FULL-weight Hamiltonian also matches the "
                                       "production mean field -- the mirrored-row weight is "
                                       "not pinned".format(name))

    def test_the_exchange_monomial_is_the_documented_one(self):
        """:func:`_ed_terms` writes ``Exchange`` in UHFk's own grouping --
        the product of the two ON-SITE spin-flip bilinears -- and that is
        the SAME operator as the documentation's

            J c^dag_{i a up} c_{j b up} c^dag_{j b dn} c_{i a dn}

        (``docs/en`` UHFk interaction file page), not a different vertex
        that happens to agree on the mean field.

        Multiplied out on the dense ``L = 2`` fixture and compared at ZERO
        tolerance, because the regrouping is an exact anticommutation (two
        transpositions and one sign) and nothing else. The grouping is what
        makes :func:`_hf_sigma`'s ``include_fock`` switch select the same
        pair of contractions as the kernel's, so it has to be checked
        rather than asserted in a comment."""
        fx = _fx_dense()
        v, R, a, b = 0.2, 1, 0, 1
        documented, grouped = [], []
        for j in range(fx.L):
            jp = (j + R) % fx.L
            ma = [fx.mode(j, a, s) for s in (UP, DN)]
            mb = [fx.mode(jp, b, s) for s in (UP, DN)]
            for s in (UP, DN):
                documented.append((ma[s], mb[s], mb[1 - s], ma[1 - s], v))
                grouped.append((ma[s], ma[1 - s], mb[1 - s], mb[s], -v))
        hd = edu.h_int_from_terms(fx, documented)
        hg = edu.h_int_from_terms(fx, grouped)
        self.assertGreater(np.abs(hd).max(), 1e-3)               # anti-vacuity
        self.assertTrue(np.array_equal(hd, hg),
                        "the grouped Exchange monomial is not the documented one")
        # and it is what _ed_terms actually emits
        rows = {"Exchange": [(R, 0, 0, a + 1, b + 1, v, 0.0)]}
        emitted = _ed_terms(fx, rows, mirrored_weight=1.0)

        def by_modes(ts):
            """``{(p, q, r, s): coeff}`` -- the term list as a multiset
            keyed on the mode indices alone. Sorting the raw tuples would
            compare their COMPLEX last element, which raises; and the
            mapping is only faithful if no two monomials share the four
            indices, which is asserted rather than assumed so that a future
            fixture fails here and not silently."""
            out = {}
            for (p, q, r, s, c_) in ts:
                self.assertNotIn((p, q, r, s), out,
                                 "two monomials share the mode indices "
                                 "{}".format((p, q, r, s)))
                out[(p, q, r, s)] = complex(c_)
            return out

        self.assertEqual(by_modes(emitted), by_modes(grouped))

    def test_first_order_gate_every_type_documented_reading(self):
        """G-HF, every type: for each off-site interaction type the kernel
        accepts, the production Hartree-Fock self-energy on a spin-coherent
        density equals the Wick derivative of the DOCUMENTED-reading ED
        Hamiltonian -- in BOTH of the kernel's Fock settings, on the FULL
        spin-orbital matrix -- and the opposite orbital placement does not.

        Why this fixture and this density. ``_fx_complex`` is the
        inversion- and transpose-asymmetric two-orbital ``L = 3`` band, so
        ``+x`` is not its own reverse and the two placements are genuinely
        different Hamiltonians. The densities are deterministic random
        Hermitian translation-invariant ones WITH SPIN COHERENCES: without
        them the whole ``PairLift`` mean field, and the spin-flip
        contractions of ``Exchange`` and ``PairHop``, vanish identically
        and the gate would be comparing zero with zero. Each case also
        carries an on-site ``CoulombIntra`` row, so the mixed on-site /
        off-site accumulation is exercised too -- and the type's OWN share
        of the mean field is measured (against the same case with the
        off-site rows removed) and required to be non-negligible, so the
        on-site row cannot make the comparison vacuous.

        Measured, relative to the mean field's own maximum, and in each
        column the value that is WORST for the assertion it supports: the
        LARGEST match residual over the two densities, and the SMALLEST
        swapped miss:

            type           flag_fock=true            flag_fock=false
                           match     miss (min)      match     miss (min)
            CoulombInter   2.0e-16   1.4e-1          1.4e-16   blind
            Hund           2.0e-16   4.5e-1          1.3e-16   blind
            Ising          2.3e-16   2.1e-1          1.7e-16   blind
            Exchange       8.0e-17   1.8e-1          0.0       blind
            PairLift       1.7e-16   1.3e-1          0.0       blind
            PairHop        1.3e-16   9.4e-2          4.1e-17   9.0e-2

        The tightest cell of the whole table is therefore ``PairHop``'s
        9.0e-2, ninety times the 1e-3 the gate demands; the largest match
        residual anywhere is 2.3e-16, four orders inside the 1e-12 it
        demands.

        "blind" is not a weaker check but a different one, asserted as an
        EQUALITY: see :data:`_HARTREE_ORIENTATION_BLIND` -- with the Fock
        term off, those types' direct contribution is built from the
        equal-site density alone and cannot depend on the displacement, so
        the two placements are the same mean field and the gate pins that
        they are. ``PairHop`` is the one type whose direct term reads the
        inter-site density, so it discriminates in both settings."""
        fx = _fx_complex()
        rhos = _random_translation_invariant_density(fx, n=2, seed=20260916,
                                                     spin_coherent=True)
        for itype, rows in _GATE_ROWS.items():
            rows_by_type = dict(_GATE_ONSITE)
            rows_by_type[itype] = rows
            terms = _ed_terms(fx, rows_by_type)
            swap_terms = _ed_terms(fx, rows_by_type, swap_orbitals=True)
            onsite_terms = _ed_terms(fx, _GATE_ONSITE)
            for fock in (True, False):
                for k, rho in enumerate(rhos):
                    with self.subTest(type=itype, include_fock=fock, density=k):
                        self.assertLess(np.abs(rho - rho.conj().T).max(), 1e-14)
                        ref = _to_k_so(fx, _hf_sigma(fx, terms, rho, include_fock=fock))
                        prod = _production_hf_sigma_k(fx, rows_by_type, rho,
                                                      include_fock=fock)
                        scale = max(np.abs(ref).max(), 1e-300)
                        self.assertGreater(scale, 1e-3)          # anti-vacuity
                        # ... and the type under test really carries part of it
                        onsite = _to_k_so(fx, _hf_sigma(fx, onsite_terms, rho,
                                                        include_fock=fock))
                        self.assertGreater(np.abs(ref - onsite).max() / scale, 1e-2,
                                           "{}: the off-site rows contribute nothing to "
                                           "the mean field on this density -- the "
                                           "comparison is vacuous".format(itype))
                        self.assertLess(np.abs(prod - ref).max() / scale, 1e-12,
                                        "{}: the production mean field is not the Wick "
                                        "derivative of the documented-reading "
                                        "Hamiltonian".format(itype))
                        swapped = _to_k_so(fx, _hf_sigma(fx, swap_terms, rho,
                                                         include_fock=fock))
                        if fock or itype not in _HARTREE_ORIENTATION_BLIND:
                            self.assertGreater(np.abs(prod - swapped).max() / scale, 1e-3,
                                               "{}: the opposite orbital placement also "
                                               "matches the production mean field -- the "
                                               "orientation is not pinned".format(itype))
                        else:
                            self.assertLess(np.abs(swapped - ref).max() / scale, 1e-12,
                                            "{}: the Hartree-only mean field DOES see the "
                                            "orientation -- _HARTREE_ORIENTATION_BLIND is "
                                            "wrong and this case must assert a miss"
                                            .format(itype))

    def test_independent_functional_equals_production_on_random_densities(self):
        """The test side's own first-order functional (:func:`_hf_sigma`,
        which is what the gates SUBTRACT) equals UHFk's mean-field kernel on
        several deterministic random Hermitian densities -- not only on the
        free one the placement pin above uses.

        The free density is paramagnetic, real and highly symmetric; an
        agreement there leaves room for a functional that differs on the
        orbital coherences, on a magnetic density, or on a complex one. The
        densities here are translation-invariant (which UHFk's kernel needs)
        and spin-block-diagonal, but otherwise generic: independent spin
        blocks, complex orbital coherences, every displacement populated."""
        cases = [("chain U+V", _fx_chain(), _rows_u_v(0.3, 0.2)),
                 ("chain Hund", _fx_chain(), _rows_offsite("Hund", 0.3)),
                 ("chain Ising", _fx_chain(), _rows_offsite("Ising", 0.3)),
                 ("orbital asymmetric", _fx_orbital(), _rows_interorbital(0.3)),
                 ("orbital Hund + V", _fx_orbital(), _rows_hund_offsite_v(0.25, 0.3))]
        for (name, fx, rows) in cases:
            terms = _ed_terms(fx, rows)
            for k, rho in enumerate(_random_translation_invariant_density(fx)):
                self.assertLess(np.abs(rho - rho.conj().T).max(), 1e-14)   # Hermitian
                ref = _production_hf_sigma_k(fx, rows, rho)
                mine = _to_k_so(fx, _hf_sigma(fx, terms, rho))
                scale = np.abs(ref).max()
                with self.subTest(case=name, density=k):
                    self.assertGreater(scale, 1e-3)                        # anti-vacuity
                    self.assertLess(np.abs(mine - ref).max(), 1e-11 * scale)

    def test_onsite_block_matches_the_second_order_compiler(self):
        """For the on-site rows alone, UHFk's mean field is
        ``hf_first_order(compile_onsite(...))`` -- the object Task 10's gate
        is built on. Recorded here so the two ED gates of the campaign are
        known to subtract the SAME first order, not two independent ones."""
        from hwave.solver.second_order import compile_onsite, hf_first_order
        fx = _fx_orbital()
        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.3, 0.0), (0, 0, 0, 2, 2, 0.2, 0.0)],
                "CoulombInter": [(0, 0, 0, 1, 2, 0.15, 0.0), (0, 0, 0, 2, 1, 0.15, 0.0)]}
        h1 = fx.build_h1()
        ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
        f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
        rho = ((W * f) @ W.conj().T).T
        ref = _production_hf_sigma_k(fx, rows, rho)
        # k-independent, as an on-site mean field must be
        np.testing.assert_allclose(ref, np.broadcast_to(ref[0], ref.shape), atol=1e-12)
        norb = fx.norb
        tbl = {t: {((0, 0, 0), (r[3] - 1, r[4] - 1)): complex(r[5], r[6]) for r in rr}
               for t, rr in rows.items()}
        rho_gen = np.zeros((2 * norb, 2 * norb), dtype=complex)
        for p in range(2 * norb):
            for q in range(2 * norb):
                sp, ap = divmod(p, norb)
                sq, aq = divmod(q, norb)
                rho_gen[p, q] = rho[fx.mode(0, ap, sp), fx.mode(0, aq, sq)]
        sig1 = hf_first_order(compile_onsite(tbl, norb), rho_gen)
        self.assertGreater(np.abs(sig1).max(), 1e-3)             # anti-vacuity
        np.testing.assert_allclose(ref[0], sig1, atol=1e-11 * np.abs(sig1).max())


class TestFirstOrder(unittest.TestCase):
    """The gate that licenses the second-order comparison: what is left
    after the Hartree-Fock subtraction really does start at O(v^2)."""

    def test_first_order_remainder_vanishes(self):
        """``Sigma_ED - Sigma_HF[rho_ED]`` has NO first order -- with the
        mirrored-row half weight, and only with it.

        Measured (O(v) remainder over the self-energy's own linear size,
        half weight / full weight): chain U 4.0e-5 (no mirrored row),
        chain V 3.3e-4 / 5.4e-1, chain U+V 2.9e-4 / 5.0e-1, chain Hund
        1.3e-4 / 4.8e-1, chain Ising 6.7e-4 / 4.9e-1, orbital asymmetric
        1.9e-4 / 5.3e-1."""
        cases = [("chain U", _fx_chain(), lambda v: _rows_u_v(v, 0.0)),
                 ("chain V", _fx_chain(), lambda v: _rows_u_v(0.0, v)),
                 ("chain U+V", _fx_chain(), lambda v: _rows_u_v(v, v)),
                 ("chain Hund", _fx_chain(), lambda v: _rows_offsite("Hund", v)),
                 ("chain Ising", _fx_chain(), lambda v: _rows_offsite("Ising", v)),
                 ("orbital asymmetric", _fx_orbital(), _rows_interorbital)]
        nmat, x = 32, 0.05
        for (name, fx, rows_of) in cases:
            iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
            mirrored = any(t in _MIRRORED_TYPES for t in rows_of(1.0))
            for weight, must_vanish in ((0.5, True), (1.0, False)):
                if weight == 1.0 and not mirrored:
                    continue      # the two conventions build the same Hamiltonian
                lin = _refine(lambda h: _sigma_fluct_ed(fx, rows_of(h), iws, weight,
                                                       "production") / h, x)
                sig, _hf = _sigma_ed(fx, rows_of(x), iws, weight, "production")
                scale = np.abs(sig).max() / x      # the self-energy's own linear size
                with self.subTest(case=name, mirrored_weight=weight):
                    self.assertGreater(scale, 1e-3)              # anti-vacuity
                    rel = np.abs(lin).max() / scale
                    if must_vanish:
                        self.assertLess(rel, _TOL,
                                        "{}: the first-order remainder does not vanish with the "
                                        "half weight ({:.2e} of the self-energy)".format(name, rel))
                    else:
                        self.assertGreater(rel, 0.1,
                                           "{}: the FULL-weight Hamiltonian also passes the "
                                           "first-order gate -- the weight is not pinned"
                                           .format(name))


class TestG4(unittest.TestCase):
    """Exact diagonalisation of the off-site second order on a chain
    (design 2026-09-15, gate G4)."""

    def _compare(self, label, ed, prod, dropped, gate=None):
        scale = np.abs(ed).max()
        self.assertGreater(scale, _FLOOR, "{}: the ED coefficient is at the floor".format(label))
        # max-norm RELATIVE error against the coefficient's own size: an
        # rtol+atol pair would let a large entry hide a proportionally large
        # error in a small one
        self.assertLess(np.abs(prod + dropped - ed).max() / scale, _TOL,
                        "{}: local + dropped vs exact diagonalisation".format(label))
        if gate is not None:
            self.assertLess(np.abs(gate - ed).max() / scale, _TOL,
                            "{}: bond gate vs exact diagonalisation".format(label))
        return scale

    def _record(self, key, value):
        """Print a RECORDED quantity and hold it to its pinned value.

        The value is not an agreement this gate claims -- it is an outcome
        it measures and routes elsewhere -- but pinning it is what turns the
        record into something a regression can fail."""
        expect = _PINNED[key]
        print("RECORDED {}: {:.3e} (pinned {:.3e} +- {:.0%})"
              .format(key, value, expect, _PINNED_BAND))
        self.assertTrue(np.isfinite(value), "{}: not finite".format(key))
        self.assertGreater(value, expect * (1.0 - _PINNED_BAND),
                           "{}: measured {:.3e}, pinned {:.3e}. A SMALLER value means the "
                           "recorded behaviour has changed for the better -- update the "
                           "expectation (and issue #192) rather than widening the band."
                           .format(key, value, expect))
        self.assertLess(value, expect * (1.0 + _PINNED_BAND),
                        "{}: measured {:.3e}, pinned {:.3e}. A LARGER value is a "
                        "regression in the recorded behaviour (issue #192 for the "
                        "inter-orbital bond gate).".format(key, value, expect))

    @heavy
    def test_coulombinter_chain(self):
        """``L = 4`` chain, on-site ``U`` and nearest-neighbour ``V``: the
        O(V^2), O(U V) and O(U^2) coefficients of the exact remainder equal
        production's ``local`` plus the oracle's dropped class, and the
        gate-on production equals the exact remainder for ``V``.

        The dropped class is load-bearing here (1.2e-1 of the O(V^2)
        coefficient), so the first assertion really does distinguish the
        local weighting from the exact one."""
        maps = _Maps(_fx_chain(), _rows_u_v)
        with _quiet():
            ed_v = _refine(lambda h: _coeff2(lambda v: maps.ed(0.0, v), h), _X)
            pr_v = _refine(lambda h: _coeff2(lambda v: maps.prod(0.0, v), h), _X)
            dr_v = maps.dropped(0.0, 1.0)
            gate_v = _refine(lambda h: _coeff2(lambda v: maps.prod(0.0, v, gate=True), h), _X)
            with self.subTest(entry="V^2"):
                scale = self._compare("V^2", ed_v, pr_v, dr_v, gate=gate_v)
                self.assertGreater(np.abs(dr_v).max(), _DROPPED_FLOOR * scale,
                                   "the dropped class is negligible on this fixture: the "
                                   "assertion would pass against the exact oracle too")

            ed_u = _refine(lambda h: _coeff2(lambda u: maps.ed(u, 0.0), h), _X)
            pr_u = _refine(lambda h: _coeff2(lambda u: maps.prod(u, 0.0), h), _X)
            with self.subTest(entry="U^2"):
                # on-site only: nothing is off-site, so nothing is dropped
                self.assertEqual(np.abs(maps.dropped(1.0, 0.0)).max(), 0.0)
                self._compare("U^2", ed_u, pr_u, 0.0)

            ed_x = _refine(lambda h: _coeff11(lambda u, v: maps.ed(u, v), h, h), _X)
            pr_x = _refine(lambda h: _coeff11(lambda u, v: maps.prod(u, v), h, h), _X)
            dr_x = maps.dropped(1.0, 1.0) - dr_v - maps.dropped(1.0, 0.0)
            with self.subTest(entry="U V"):
                # mixed on/off pairs are never dropped (one entry is on-site,
                # so the pair is always "crossed"); measured 1.7e-18
                self.assertLess(np.abs(dr_x).max(), _FLOOR * np.abs(ed_x).max())
                self._compare("U V", ed_x, pr_x, dr_x)

    @heavy
    def test_interorbital_bond_chain(self):
        """``L = 3``, two orbitals, ORBITAL-ASYMMETRIC inter-orbital bond:
        the adjudication of the campaign's two open items.

        ASSERTED -- the standalone local path plus the dropped class
        reproduces the exact O(v^2) remainder (3.5e-4 of it). Since the ED
        Hamiltonian is pinned to UHFk's own mean field
        (:meth:`TestChainHamiltonian.test_ed_hamiltonian_is_the_production_
        mean_field`), this is an independent confirmation of the ORIENTATION
        of the off-site density vertex -- the DOCUMENTED reading, orbital
        ``a`` in the original cell: the OPPOSITE orbital placement of
        the same declaration, run through the identical recipe, misses
        ``local + dropped`` by 2.2e-2 -- eleven times the tolerance, and the
        assertion below does fail loudly on it (that is exactly the value it
        failed by while the local kernel still read the reversed
        orientation, before the fix of issue #193).

        RECORDED -- the bond gate's deviation. Phase B's bond-resolved path
        is supposed to resum exactly the class the local kernel drops
        (2.2e-2 of the coefficient here), and on the single-orbital chain it
        does (2.4e-4). On this inter-orbital bond it misses by 4.0e-2 --
        twenty times the tolerance, and nearly twice the class it is
        recovering, so the deviation is not that class being mishandled but
        a bond-gate defect on inter-orbital off-site bonds. It is printed
        and routed to a Phase B follow-up of issue #181; it is NOT a
        second-order kernel finding, and it is not asserted here."""
        maps = _Maps(_fx_orbital(), _rows_interorbital)
        with _quiet():
            ed = _refine(lambda h: _coeff2(lambda v: maps.ed(v), h), _X)
            pr = _refine(lambda h: _coeff2(lambda v: maps.prod(v), h), _X)
            dr = maps.dropped(1.0)
            gate = _refine(lambda h: _coeff2(lambda v: maps.prod(v, gate=True), h), _X)
            scale = self._compare("inter-orbital v^2", ed, pr, dr)
            self.assertGreater(np.abs(dr).max(), _DROPPED_FLOOR * scale,
                               "the dropped class is negligible on this fixture")
            for arr, what in ((ed, "ED"), (pr, "local"), (dr, "dropped"), (gate, "gate")):
                self.assertTrue(np.all(np.isfinite(arr)),
                                "the {} coefficient is not finite".format(what))
            # Phase B follow-up of issue #181, tracked as issue #192
            self._record("interorbital gate vs ED", np.abs(gate - ed).max() / scale)
            self._record("interorbital dropped class", np.abs(dr).max() / scale)

    @heavy
    def test_hund_times_offsite_v_chain(self):
        """``L = 3``, two orbitals: the MIXED O(J V) coefficient of on-site
        ``Hund`` against the asymmetric inter-orbital off-site ``V``.

        The pure entries of this campaign exercise the on-site kernel and
        the off-site one; this is the mixed class -- the cross products
        ``A_on chibar B_v + A_v chibar B_on`` -- with an on-site vertex that
        is not a plain density (on-site ``Hund`` carries same-spin monomials
        only), which is the shape the weight rule of spec 2.4 says is local
        with weight 1 regardless of the off-site vertex's placement.

        Two statements:

        * the mixed coefficient of ``local + dropped`` equals the exact
          remainder's, at this module's tolerance and above its floor;
        * the mixed part of the DROPPED class is identically zero -- every
          mixed on/off diagram has a representable copy, so the local
          kernel drops none of them. This is asserted rather than assumed:
          were it nonzero, the first statement would be comparing against a
          correction the kernel is not supposed to need.

        Measured at this module's working point: the ED coefficient is
        1.07e-1, ``local + dropped`` reproduces it to 4.5e-4 of it (a factor
        4.4 inside the tolerance) and the mixed dropped part is 4.6e-18 of
        it -- zero to round-off.
        """
        maps = _Maps(_fx_orbital(), _rows_hund_offsite_v)
        with _quiet():
            ed = _refine(lambda h: _coeff11(lambda j, v: maps.ed(j, v), h, h), _X)
            pr = _refine(lambda h: _coeff11(lambda j, v: maps.prod(j, v), h, h), _X)
            # the oracle's skeleton is exactly quadratic, so the mixed part of
            # the dropped class is the unit-coupling value minus the two pure ones
            dr = maps.dropped(1.0, 1.0) - maps.dropped(1.0, 0.0) - maps.dropped(0.0, 1.0)
        scale = np.abs(ed).max()
        self.assertGreater(scale, _FLOOR,                        # anti-vacuity
                           "the mixed Hund x V coefficient is at the floor")
        for arr, what in ((ed, "ED"), (pr, "local"), (dr, "dropped")):
            self.assertTrue(np.all(np.isfinite(arr)),
                            "the mixed {} coefficient is not finite".format(what))
        # every mixed on/off diagram is local: nothing is dropped
        self.assertLess(np.abs(dr).max(), _FLOOR * scale,
                        "the mixed on/off class has a dropped part: the weight rule of "
                        "spec 2.4 says every mixed diagram is local")
        self.assertLess(np.abs(pr + dr - ed).max() / scale, _TOL,
                        "Hund x off-site V: local + dropped vs exact diagonalisation")

    @heavy
    def test_hund_ising_chain_recorded(self):
        """Off-site Hund and Ising on the ``L = 4`` chain, RECORDED (spec
        G2 (c) / G4): the same recipe as the CoulombInter chain, printed
        rather than asserted, because the off-site Hund/Ising second order
        is the follow-up class of issue #181's Tier 3 and not part of this
        gate's asserted surface. Measured: 2.5e-4 (Hund) and 2.2e-4
        (Ising) for both the standalone and the gate-on path."""
        for itype in ("Hund", "Ising"):
            maps = _Maps(_fx_chain(), lambda v, itype=itype: _rows_offsite(itype, v))
            with _quiet():
                ed = _refine(lambda h: _coeff2(lambda v: maps.ed(v), h), _X)
                pr = _refine(lambda h: _coeff2(lambda v: maps.prod(v), h), _X)
                gate = _refine(lambda h: _coeff2(lambda v: maps.prod(v, gate=True), h), _X)
                dr = maps.dropped(1.0)
            scale = np.abs(ed).max()
            self.assertGreater(scale, _FLOOR)                    # anti-vacuity
            for arr, what in ((ed, "ED"), (pr, "local"), (dr, "dropped"), (gate, "gate")):
                self.assertTrue(np.all(np.isfinite(arr)),
                                "the {} coefficient of {} is not finite".format(what, itype))
            self._record("{} local+dropped vs ED".format(itype),
                         np.abs(pr + dr - ed).max() / scale)
            self._record("{} gate vs ED".format(itype), np.abs(gate - ed).max() / scale)
            self._record("{} dropped class".format(itype), np.abs(dr).max() / scale)


if __name__ == "__main__":
    unittest.main()
