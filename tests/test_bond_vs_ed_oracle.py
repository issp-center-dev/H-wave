#!/usr/bin/env python3

"""ED adjudication of the bond-resolved ph/pp channels (#151) -- the
campaign's main comparison module.

See ``docs/superpowers/specs/2026-08-13-bond-ed-adjudication-design.md``
for the full design and ``docs/superpowers/plans/2026-08-13-bond-ed-
adjudication.md`` for the task-by-task build order. This module owns the
shared fixture factories (lazy, cached -- constructing an ``EDFixture``
costs nothing, but the first Lehmann diagonalization on an ``L=3, norb=2``
system is exponential-in-Hilbert-space and must happen at most once per
fixture) plus the module-wide constants every later task in the campaign
reuses.

Two interaction representations recur throughout the campaign and must
never be confused (design doc, "Declaration counting"): the ED side
(Hamiltonian and Hartree-Fock) takes a CANONICAL list -- one entry per
unordered displaced-density-pair class, via
``ed_oracle_util.canonical_density_terms`` -- while the SOLVER side takes
the CLOSED TWO-SIDED declaration dict production's own
``resolve_interactions`` expects (both ``V_ab(R)`` and ``V_ba(-R)``
declared explicitly). Feeding the closed two-sided list, undeduplicated,
to the ED side double-counts the Hamiltonian; this module's first test
pins that factor of 2 so the convention cannot silently drift (see the
class docstring below).

The second and third tests pin the design's claimed "null direction": for
a real density-density interaction, ``V_ab(R) = v + i*eps`` declared
alongside its Hermitian-closed partner ``V_ba(-R) = v - i*eps`` is, after
site relabeling, one operator weighted by ``(v+i*eps)+(v-i*eps) = 2v`` --
independent of ``eps``. The design claims this holds "ED and solver
alike". The ED half holds at machine precision (verified below). The
solver half does NOT survive contact with ``bare_bond_vertices`` -- see
``TestNullDirectionSolverSide`` and this task's report for the measured
finding.
"""

import collections
import functools
import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa
from hwave.sc import _build_bond_m0_blocks
from hwave.solver.bond_channels import (
    bare_bond_vertices, bond_bubble, dress_bond, resolve_interactions,
)
from tests import ed_oracle_util
from tests.approx_util import ApproxTestCase, assert_approx_array


V1 = 0.02      # base coupling used throughout the campaign for calibration
               # pins and static (non-Richardson) identity checks -- the
               # LEGACY value; NOT the ph-adjudication Richardson step (see
               # CAMPAIGN_V1 below). Left untouched: every remaining use of
               # V1 in this module is a representative fixed coupling
               # amplitude (canonical-list counting, the null-direction
               # pin, Pin 2a, the two-shell dW joint test, etc.), never a
               # finite-difference step size.
CAMPAIGN_V1 = 0.005
               # The campaign's ph-adjudication Richardson step (design doc,
               # "Fixtures", amended 2026-08-14 from the draft's 0.02): Task
               # 6 MEASURED that at V1 = 0.02 the off-site (g1/g2) unit
               # directions' Richardson convergence was too loose to clear
               # the granule power check (delta_rich ~1e-3, dominated by a
               # near-null cross-term cell and ordinary O(V1^2) curvature on
               # the off-site directions) even though every cell already
               # agreed with the solver prediction -- a finer scan confirmed
               # clean convergence toward the solver values and delta_rich
               # scaling ~ V1^2, so this halves the design's original step
               # to 0.005 (V2 = 2*CAMPAIGN_V1, halved pair CAMPAIGN_V1/2 =
               # 0.0025) to restore the power margin. The 2-site legacy ED
               # oracle (tests/test_rpa_vs_ed_oracle.py) keeps its own 0.02
               # and is NOT touched by this constant.
NMAT = 1024    # solver-side Matsubara count for the calibration/adjudication
               # pins (Tasks 5+); defined here as the shared module constant.


@functools.lru_cache(maxsize=1)
def fx5():
    """5-site, single-orbital ring: the single-orbital ED fixture.

    ``t`` is a dict keyed by orbital pair, per ``EDFixture``'s actual
    contract (``build_h1`` reads ``t.items()``) -- for norb=1 this is the
    single entry ``t[(0, 0)]``.
    """
    t = {(0, 0): 0.7 * np.exp(0.3j)}
    return ed_oracle_util.EDFixture(L=5, norb=1, t=t, eps=(0.0,), T=0.5, mu=0.2)


@functools.lru_cache(maxsize=1)
def fx3():
    """3-site, two-orbital ring: the case-M ED fixture (complex intra- and
    inter-orbital hopping, distinct per-orbital on-site energies).

    ``t`` entries are the four (a, b) elements of the campaign's hopping
    matrix ``[[0.7e^{0.3i}, 0.25+0.15i], [0.25+0.15i, 0.5e^{-0.2i}]]``.
    """
    t = {
        (0, 0): 0.7 * np.exp(0.3j),
        (1, 1): 0.5 * np.exp(-0.2j),
        (0, 1): 0.25 + 0.15j,
        (1, 0): 0.25 + 0.15j,
    }
    return ed_oracle_util.EDFixture(
        L=3, norb=2, t=t, eps=(0.10, -0.05), T=0.5, mu=0.2)


@functools.lru_cache(maxsize=1)
def fx2():
    """2-site, two-orbital ring: the extraction-gate fixture -- matches
    ``tests/test_rpa_vs_ed_oracle.py``'s ``FX2``/``HOP`` exactly."""
    t = {
        (0, 0): 0.7 + 0.3j,
        (1, 1): 0.5 - 0.2j,
        (0, 1): 0.25 + 0.15j,
        (1, 0): 0.25 + 0.15j,
    }
    return ed_oracle_util.EDFixture(
        L=2, norb=2, t=t, eps=(0.10, -0.05), T=0.5, mu=0.2)


# ---------------------------------------------------------------------------
# Canonical-list pins
# ---------------------------------------------------------------------------

class TestCanonicalListCounting(ApproxTestCase):
    """Pins the canonical-list convention (design doc: "Declaration
    counting, fixed"): a displaced density pair ``(a, b, R)`` and its
    Hermitian-closed partner ``(b, a, -R)`` are the SAME physical operator
    after site relabeling. ``canonical_density_terms`` must be fed one
    entry per unordered class -- feeding both is a silent factor-2 hazard
    the design flags explicitly, and could otherwise cancel between the ED
    Hamiltonian and the HF subtraction while leaving a production
    comparison wrong by 2. Uses ``fx5`` (single-orbital, off-site same-
    orbital bond): the duplication/nullity argument only needs one bond
    class and is exercised in full generality (all four spin pairings) by
    the ``otherwise`` branch of ``canonical_density_terms`` regardless of
    ``norb``.
    """

    def test_duplicated_closed_list_is_exactly_twice_canonical(self):
        fx = fx5()
        canonical = ed_oracle_util.canonical_density_terms(fx, [(0, 0, 1, V1)])
        duplicated = ed_oracle_util.canonical_density_terms(
            fx, [(0, 0, 1, V1), (0, 0, -1, V1)])
        h_canon = ed_oracle_util.h_int_from_terms(fx, canonical)
        h_dup = ed_oracle_util.h_int_from_terms(fx, duplicated)
        assert_approx_array(h_dup, 2.0 * h_canon, rel=0, abs=1e-14)

    def test_imaginary_closed_pair_direction_is_null(self):
        # The design's "null direction": declaring V_R(+1) = V1 + i*eps
        # alongside its Hermitian-closed partner V_R(-1) = V1 - i*eps sums,
        # after site relabeling, to (V1+i*eps) + (V1-i*eps) = 2*V1 on the
        # SAME operator -- independent of eps. Perturbing eps must move the
        # ED Hamiltonian by EXACTLY zero, not just approximately.
        fx = fx5()
        eps = 1e-3
        base = ed_oracle_util.canonical_density_terms(
            fx, [(0, 0, 1, V1), (0, 0, -1, V1)])
        pert = ed_oracle_util.canonical_density_terms(
            fx, [(0, 0, 1, V1 + 1j * eps), (0, 0, -1, V1 - 1j * eps)])
        h_base = ed_oracle_util.h_int_from_terms(fx, base)
        h_pert = ed_oracle_util.h_int_from_terms(fx, pert)
        assert_approx_array(h_pert, h_base, rel=0, abs=1e-16)


# ---------------------------------------------------------------------------
# Solver-side null-direction pin
# ---------------------------------------------------------------------------

_Fx3Resolved = collections.namedtuple(
    "_Fx3Resolved", ["bond_set", "S_bond", "C_bond", "Vpp_s", "Vpp_t"])


def _two_sided_decl(entries):
    """The CLOSED TWO-SIDED ``CoulombInter`` declaration dict
    ``resolve_interactions`` expects (Global Constraints: solver
    representation), built from a list of ONE-SIDED physical directions
    ``[(R, a, b, value), ...]`` -- ``R`` a single-axis ring displacement
    (the fixtures in this module are 1-D rings, so only the x-component of
    the Bravais-cell displacement is ever nonzero), ``a``, ``b`` PHYSICAL
    orbital indices. Declares both ``V_ab(R) = value`` and its Hermitian-
    closed partner ``V_ba(-R) = conj(value)`` explicitly -- the only pairing
    ``resolve_interactions`` accepts without raising, since it enforces
    ``V_ab(R) == conj(V_ba(-R))``. Shared by every closed two-sided
    declaration this module builds (Task 3's ``_resolve_fx3`` and Task 5's
    ``dW_matrices``/frame-map fixtures), so the convention cannot drift
    between them.
    """
    decl = {}
    for (R, a, b, v) in entries:
        decl[((R, 0, 0), (a, b))] = v
        decl[((-R, 0, 0), (b, a))] = np.conj(v)
    return decl


def _resolve_fx3(perturbation):
    """Resolve fx3's ring interaction through the SOLVER's own builders
    (``resolve_interactions`` then ``bare_bond_vertices``), for the
    (a=0, b=1, R=+1) bond direction the design's null-direction example
    uses (spec, "direction g+/g-").

    ``perturbation`` is ``{}`` for the base declaration
    ``{V_01(+1): V1, V_10(-1): V1}``, or ``{"im_eps": eps}`` to add the
    Hermitian-closed null perturbation ``V_01(+1) += +i*eps``,
    ``V_10(-1) += -i*eps`` -- the SAME sign convention pinned by
    ``test_imaginary_closed_pair_direction_is_null`` above (opposite signs
    are the only choice ``resolve_interactions`` accepts without raising:
    same-sign eps on both entries would violate the V_ab(R) ==
    conj(V_ba(-R)) Hermiticity check it enforces).

    ``lattice_vectors=np.eye(3)`` matches every ``resolve_interactions``
    caller in ``tests/test_bond_channels.py`` -- it is used only for
    Euclidean shell-length sorting, not for the fixture's actual ring
    geometry, and no ``bond_max_shells`` cutoff is in play here.

    ``bare_bond_vertices`` assumes ``C0_q`` already carries the FULL q=0
    Hartree ``2*V_ab(q=0)`` -- its Case-2-correction step subtracts the
    R!=0 part back out to isolate the local Cooper block (see its
    docstring / spec S4.3 star, and
    ``tests/test_bond_channels.py``'s ``_build_complex_two_orbital_bond``,
    which documents and exercises the identical precondition). Passing a
    bare zero ``C0_q`` would violate that precondition and inject a
    spurious residual into the local block unrelated to the null-
    direction question under test, so ``C0_q``'s (aa,bb) Hartree entries
    are filled from ``bond_set.v_bond`` exactly as that helper does:
    ``2 * sum_{m>=1} V_ab(R_m)``.

    Returns the resolved topology plus the ``bare_bond_vertices`` outputs.
    """
    eps = perturbation.get("im_eps", 0.0)
    coulomb_inter = _two_sided_decl([(1, 0, 1, V1 + 1j * eps)])
    norb = fx3().norb
    bond_set = resolve_interactions(coulomb_inter, np.eye(3), norb=norb)
    nd = norb * norb
    S0_q = np.zeros((1, 1, 1, nd, nd), dtype=complex)
    C0_q = np.zeros((1, 1, 1, nd, nd), dtype=complex)
    for a in range(norb):
        for b in range(norb):
            total_ab = sum(
                bond_set.v_bond[m][a, b] for m in range(1, bond_set.n_channels))
            C0_q[0, 0, 0, a * norb + a, b * norb + b] = 2.0 * total_ab
    S_bond, C_bond, Vpp_s, Vpp_t = bare_bond_vertices(bond_set, S0_q, C0_q, norb)
    return _Fx3Resolved(bond_set, S_bond, C_bond, Vpp_s, Vpp_t)


class TestNullDirectionSolverSide(ApproxTestCase):
    """Solver-side half of the null-direction pin (review round 2): the ph
    matrices ``S_bond``/``C_bond`` AND the pp vertices ``Vpp_s``/``Vpp_t``
    from ``bare_bond_vertices`` must be unchanged (abs=1e-15) under the
    imaginary null perturbation -- mirroring
    ``test_imaginary_closed_pair_direction_is_null``'s ED-side result.
    Each of the four matrices is checked in its own ``subTest`` so a
    failure on one does not hide the others' pass/fail status.

    MEASURED FINDING (discrepancy protocol -- do not touch
    ``bond_channels.py``): the PP sector holds -- ``Vpp_s``/``Vpp_t`` are
    exactly null-invariant (diff 0.0), because the ``D + D^dagger``
    construction that feeds the bond Cooper block takes ``2*Re(.)`` of
    each bond-diagonal entry (killing the imaginary offset before any
    sandwich), and the local block's R!=0 Hartree subtraction cancels
    against a correctly-filled ``C0_q`` exactly, independent of ``eps``.

    The PH sector does NOT hold. ``S_bond``/``C_bond``'s m!=0 bond-
    diagonal Fock element is set directly from ``v_bond[m][l1, l2]``
    (S) / ``-v_bond[m][l1, l2]`` (C) with no Hermitizing combination at
    all -- ``V_01(+1)`` and ``V_10(-1)`` live at DIFFERENT enlarged
    indices ``(m, l1, l2)`` and are never summed onto one physical slot
    the way the ED canonical list sums them onto one operator -- so the
    imaginary offset survives linearly: ``S_bond`` responds at magnitude
    ``eps`` (the bond-diagonal element alone), ``C_bond`` at magnitude
    ``2*eps`` (that same bond-diagonal element PLUS the m=0 Hartree
    sub-block, which is ``2*eps``-responsive through
    ``V_01(q=0) = V_01(+1) + V_01(-1)``). This is an INTERMEDIATE-VERTEX
    finding at the ``bare_bond_vertices`` output -- whether it survives
    into a fully-dressed physical observable (``dress_bond``,
    ``make_bond_kernel``, the eventual chi/gap) is exactly what Tasks 4-9
    adjudicate; that projection step is out of this task's scope. This
    test is left asserting the null (and failing on ``S_bond``/``C_bond``)
    rather than adjusted to match production, per the discrepancy
    protocol; see this task's report for the measured values. No custom
    ``msg`` is passed to ``assert_approx_array`` below -- its default
    failure message already carries the mismatch count, first bad index,
    actual/expected values, diff and tolerance, which the discrepancy
    protocol's "record the measured values" step needs; a custom message
    would replace (not augment) that detail.
    """

    def test_null_direction_moves_nothing_solver_side(self):
        base = _resolve_fx3({})
        pert = _resolve_fx3({"im_eps": 1e-3})
        for name in ("S_bond", "C_bond", "Vpp_s", "Vpp_t"):
            with self.subTest(matrix=name):
                assert_approx_array(
                    getattr(pert, name), getattr(base, name), rel=0, abs=1e-15)


# ---------------------------------------------------------------------------
# Sector-block engine vs the dense Lehmann path (#151, Task 4)
# ---------------------------------------------------------------------------

class TestSectorBlockEngine(ApproxTestCase):
    """``SectorED`` must reproduce the dense ``chi_connected`` path exactly
    (round-off), including a genuine interaction, and must produce the
    right-shaped tensor on the case-M fixture the dense path cannot afford
    (fx3: nmode=12, dim=4096 -- 12 dense (4096, 4096) annihilators alone
    would be multiple GB; SectorED never builds those)."""

    def test_block_matches_dense_with_interaction(self):
        fx = fx2()
        # Mixes both canonical_density_terms branches: an on-site U
        # (a==b, R==0), an on-site inter-orbital V (a!=b, R==0), and an
        # off-site same-orbital V (a==b, R!=0) -- exercises the general
        # nonlocal-Fock HF path is NOT what's under test here (h1 stays
        # the bare one-body matrix); this is purely the interacting-H
        # sector-vs-dense comparison.
        terms = ed_oracle_util.canonical_density_terms(
            fx, [(0, 0, 0, V1), (0, 1, 0, V1), (0, 0, 1, V1)])
        hint = ed_oracle_util.h_int_from_terms(fx, terms)
        dense = ed_oracle_util.chi_connected(fx, hint=hint)
        block = ed_oracle_util.SectorED(fx, terms=terms).chi_connected()
        assert_approx_array(block, dense, rel=0, abs=1e-11)

    def test_chi_connected_shape_case_m(self):
        fx = fx3()
        out = ed_oracle_util.SectorED(fx).chi_connected()
        self.assertEqual(out.shape, (3, 4, 4, 4, 4))


class TestBondPairCorrelatorApi(ApproxTestCase):
    """``bond_correlator``/``pair_correlator`` shape and the internal-
    consistency pin against the already-verified ``chi_connected`` (#151,
    Task 4)."""

    def test_bond_correlator_shape(self):
        fx = fx5()
        sector = ed_oracle_util.SectorED(fx)
        channels = [(0,), (1,), (-1,)]
        out = sector.bond_correlator(channels)
        self.assertEqual(out.shape, (5, 2, 2, 3, 3))

    def test_bond_correlator_diagonal_matches_chi_connected(self):
        # Index correspondence (derived, not scanned): for the R=0, a=b
        # bond channel at spin sigma, B_{0,a,a,sigma}(q) = sum_j e^{iqj}
        # c^dag_{j,a,sigma} c_{j,a,sigma} is LITERALLY chi_connected's
        # density operator O(q, x, x) with generalized index
        # x = sigma*norb + a (the a = s*norb+o convention pinned in the
        # plan's Global Constraints).
        #
        # SAME-SPIN (sigma == sigma'): bond_correlator's <B; B^dagger> and
        # chi_connected's <(c^dag_a c_c)(q); (c^dag_d c_b)(-q)> both reduce,
        # at a=c=b=d=x, to sum_{m,n} K[m,n] O(q,x,x)[m,n] conj(O(q,x,x)[m,n])
        # minus the SAME disconnected piece beta*|<O(q,x,x)>|^2 -- because
        # O(-q,x,x)[n,m] = conj(O(q,x,x)[m,n]) exactly for any Hermitian
        # single-site density operator's Fourier transform (a q -> -q
        # conjugation identity, independent of the fixture).
        #
        # CROSS-SPIN (sigma != sigma', review-derived): the SAME
        # Hermiticity identity applies independently to EACH leg, since it
        # only needs a single-index density operator, not a==c between
        # the two legs. With A = B_{0,a,a,sigma}(q) = O(q, x_sigma, x_sigma)
        # and B = B_{0,a,a,sigma'}(q) = O(q, x_sigma', x_sigma'):
        #   B^dagger(q) = O(-q, x_sigma', x_sigma')
        #   B[m,n]^dagger's matrix element conj(B(q)[m,n])
        #     = O(-q, x_sigma', x_sigma')[n,m]   (same identity, x = x_sigma')
        # so Xph[q,sigma,sigma',I,I]
        #   = sum K[m,n] O(q,x_sigma,x_sigma)[m,n] O(-q,x_sigma',x_sigma')[n,m]
        #   = chi[q, x_sigma, x_sigma, x_sigma', x_sigma']
        # (matching chi_connected's out[qi, a, c, b, d] slot with
        # a=c=x_sigma, b=d=x_sigma'), and the disconnected pieces coincide
        # by the same per-leg conjugation identity applied to the average.
        # Both use the SAME 1/L normalization, so the slots match at
        # round-off for ANY interaction, not just V=0.
        fx = fx3()
        terms = ed_oracle_util.canonical_density_terms(fx, [(0, 1, 1, V1)])
        sector = ed_oracle_util.SectorED(fx, terms=terms)
        channels = [(0, 0, 0), (0, 1, 1)]
        xph = sector.bond_correlator(channels)
        chi = sector.chi_connected()
        for i, (r, a, b) in enumerate(channels):
            for sigma in range(2):
                for sigmap in range(2):
                    x_sigma = sigma * fx.norb + a
                    x_sigmap = sigmap * fx.norb + a
                    assert_approx_array(
                        xph[:, sigma, sigmap, i, i],
                        chi[:, x_sigma, x_sigma, x_sigmap, x_sigmap],
                        rel=0, abs=1e-11)

    def test_pair_correlator_smoke(self):
        # channels includes an R != 0 entry (fx2 is L=2, so R=1 is the one
        # nontrivial displacement) so _apply_pair's (j+R) % L offset path
        # actually executes at least once -- with only R=0 channels the
        # offset was dead code (review finding, #151 Task 4 fix loop).
        fx = fx2()
        sector = ed_oracle_util.SectorED(fx)
        channels = [(0, 0, 0), (0, 1, 1), (0, 0, 1), (1, 0, 0)]
        xpp = sector.pair_correlator(channels)
        self.assertEqual(xpp.shape, (fx.L, 4, 4))
        self.assertTrue(np.all(np.isfinite(xpp)))
        self.assertTrue(np.any(xpp != 0))
        # Xpp[q] is Hermitian AT FIXED q (not a q <-> -q relation): both
        # Xpp[q,i,j] = <Delta_i(q); Delta_j(q)^dagger> and
        # Xpp[q,j,i] = <Delta_j(q); Delta_i(q)^dagger> sum over the SAME
        # (m, n) sector-pair domain (the pair operator's delta is always
        # (-1, -1)) against the SAME real, symmetric kernel K[m,n] = K[n,m]
        # (the static Lehmann kernel is real and symmetric by construction
        # -- see _static_kernel), so term-by-term
        # Xpp[q,j,i] = conj(Xpp[q,i,j]).
        #
        # CAUTION (review finding, #151 Task 4 fix loop): this identity is
        # ALGEBRAICALLY VACUOUS as a correctness check on Delta itself --
        # it holds for ANY family of operators {A_i} sharing one sector
        # shift, against ANY real-symmetric kernel, regardless of whether
        # Delta's site offset, up/down leg assignment or annihilation-
        # order sign are right. It is a (cheap, worth keeping) structural
        # sanity check on _lehmann_dagger/_pair_kernel, not a numerical
        # pin on the pair path. The pair path's actual correctness --
        # including the mu-sensitive cross-sector Lehmann denominator this
        # (N_up-1, N_dn-1) operator is the only in-repo path to exercise
        # (chi_connected's cross-sector slots are all Delta-N = 0, where mu
        # cancels between the legs) -- is deferred to Task 7's pin 3b
        # (eigenbasis_pair_bubble vs pair_correlator at V=0), which is
        # therefore LOAD-BEARING and must not be weakened to a diagnostic.
        for qi in range(fx.L):
            assert_approx_array(xpp[qi], xpp[qi].conj().T, rel=0, abs=1e-11)


# ---------------------------------------------------------------------------
# Task 5: solver-side sandwich, frame map and calibration pins (#151)
# ---------------------------------------------------------------------------
#
# Five interfaces:
#   free_green            -- the free k-space Matsubara Green function, in
#                             sc.py's (norb, norb, Nx, Ny, Nz, nmat) layout.
#   bond_setup             -- (bond_set, chi_bar) from a closed two-sided
#                             declaration, via resolve_interactions +
#                             free_green + bond_bubble.
#   dW_matrices            -- (dS, dC), the bare-vertex response to a closed
#                             two-sided declaration (bare_bond_vertices at
#                             the declaration minus at zero).
#   ed_to_solver_bond_map   -- the DERIVED index/conjugation map between
#                             SectorED.bond_correlator's (q, I, J) and
#                             bond_bubble's (q, m*nd+idx, m'*nd+idx').
#   pred_first_order        -- the exact-in-coupling solver-side sandwich
#                             prediction (+chibar@dS@chibar, -chibar@dC@chibar).
#
# Three calibration pins (both fx5 and fx3):
#   Pin 1  ED chi_connected(V=0)      vs solver's TAIL-CORRECTED chi0q
#   Pin 2a bond_bubble's Delta r=0    vs solver's RAW chi0 (same Nmat)
#   Pin 2b ed_to_solver_bond_map,     full (m, m') matrix, zero coupling


def free_green(fx, nmat):
    """The free k-space Matsubara Green function, sc.py layout
    ``(norb, norb, Nx, Ny, Nz, nmat)`` with ``Nx=fx.L, Ny=Nz=1``.

    ``eps(k) = t e^{+ik} + t^dagger e^{-ik} + diag(eps)`` on
    ``k = 2*pi*kx/L`` (``t`` assembled from ``fx.t`` -- ``(t^dagger)_{ab} =
    conj(t_{ba})``, matching ``EDFixture.build_h1``'s own hopping
    convention exactly, so this function and the ED Hamiltonian can never
    silently disagree on which orbital pair carries which direction);
    ``iomega_n = (2n+1-nmat)*pi/beta``; ``G = [(iomega + mu) I -
    eps(k)]^{-1}``. This is ALSO the ``epsilon_k[a,b] += t_ab *
    e^{+ikR}`` convention ``hwave.sc._build_hamiltonian_k`` uses for a
    Wannier90-style ``transfer.dat`` declaration at ``R=+1``, orbital pair
    ``(a,b)`` -- verified end to end by Pin 1/2a below (the solver's OWN
    Green function, built from a ``transfer.dat`` written with ``fx.t``,
    reproduces this function's output through ``bond_bubble``/``chi0q`` at
    round-off).

    ``mu``/``beta``/``eps`` are read from ``fx`` (``fx.mu``, ``fx.beta``,
    ``fx.eps``); orbital indices index ``fx.t``, a dict keyed by physical
    orbital pair (EDFixture's own contract).
    """
    beta = fx.beta
    L = fx.L
    norb = fx.norb
    tmat = np.zeros((norb, norb), dtype=complex)
    for (a, b), v in fx.t.items():
        tmat[a, b] = v
    iomega = (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi / beta
    kx = 2.0 * np.pi * np.arange(L) / L
    phase = np.exp(1j * kx)                              # (L,)
    eps_k = (tmat[None, :, :] * phase[:, None, None]
             + tmat.conj().T[None, :, :] * np.conj(phase)[:, None, None])
    eps_k = eps_k + np.diag(np.asarray(fx.eps, dtype=complex))[None, :, :]
    ident = np.eye(norb, dtype=complex)
    mat = ((1j * iomega + fx.mu)[None, :, None, None]
           * ident[None, None, :, :]
           - eps_k[:, None, :, :])                        # (L, nmat, norb, norb)
    green_ln = np.linalg.inv(mat)                          # (L, nmat, norb, norb)
    green = np.moveaxis(green_ln, (2, 3), (0, 1))           # (norb, norb, L, nmat)
    return green[:, :, :, np.newaxis, np.newaxis, :]         # (norb, norb, L, 1, 1, nmat)


def bond_setup(fx, closed_decls, **kw):
    """(bond_set, chi_bar) from a closed two-sided production declaration.

    ``closed_decls``: the TWO-SIDED production declaration dict
    (``_two_sided_decl``'s output) -- the SOLVER representation, per the
    Global Constraints (never fed a canonical/ED-side list here).
    ``**kw`` splits into ``nmat`` (default the module ``NMAT`` constant,
    the solver-side Matsubara count) and any remaining keyword arguments
    forwarded verbatim to ``resolve_interactions`` (e.g. ``bond_max_shells``).
    """
    nmat = kw.pop("nmat", NMAT)
    bond_set = resolve_interactions(closed_decls, np.eye(3), norb=fx.norb, **kw)
    green = free_green(fx, nmat)
    chi_bar = bond_bubble(green, bond_set, beta=fx.beta)
    return bond_set, chi_bar


def _bare_vertices_at(fx, closed_decls):
    """``(S_bond, C_bond)`` at ``q=(kx, 0, 0)`` for every ``kx`` on
    ``fx``'s ring, from ``closed_decls`` through the SAME q-dependent
    Case-2-corrected ``m=0`` block production itself uses
    (``hwave.sc._build_bond_m0_blocks``) -- not a local reimplementation,
    so ``dW_matrices`` can never silently drift from the production
    Hartree convention (``V_ab(q) = V_ab(0) + sum_{m>=1} V_ab(Delta r_m)
    e^{-i q.Delta r_m}``, spec S4.3/S4.3-star).
    """
    bond_set = resolve_interactions(closed_decls, np.eye(3), norb=fx.norb)
    kx = 2.0 * np.pi * np.arange(fx.L) / fx.L
    ky = np.array([0.0])
    kz = np.array([0.0])
    S0_q, C0_q = _build_bond_m0_blocks(bond_set, {}, {}, fx.norb, kx, ky, kz)
    S_bond, C_bond, _Vpp_s, _Vpp_t = bare_bond_vertices(bond_set, S0_q, C0_q, fx.norb)
    return S_bond[:, 0, 0], C_bond[:, 0, 0]


def dW_matrices(fx, closed_decls):
    """``(dS, dC)`` per physical direction: ``bare_bond_vertices`` at the
    (two-sided, closed) ``closed_decls`` declaration minus at the SAME
    declared topology scaled to zero (so the declared-but-zero shells stay
    part of ``ND`` on both sides -- a same-shape subtraction, not a
    smaller-vs-larger one).

    For L=5, one shell declared two-sided at amplitude g (single orbital):
    the C0 m=0 Hartree slot (aa),(bb) carries ``dC0/dg = 2*(e^{+iq} +
    e^{-iq}) = 4*cos(q)``; the m=+1 and m=-1 diagonal blocks carry
    ``dS/dg = +1`` and ``dC/dg = -1`` each; all other entries zero (signs
    and the ``2*V(q)`` convention from ``bare_bond_vertices``' documented
    C0 Hartree ``+2 V_ab(q)``). Two-shell direction g2 (m=+-2) carries the
    same pattern with ``C0 = 4*cos(2q)``. Case-M's g+/g- (inter-orbital,
    a!=b) directions carry the SAME m!=0 diagonal pattern at the
    (a,b)/(b,a) orbital-pair slot, plus a q-dependent CROSS Hartree entry
    at the (aa),(bb) slot with a!=b (``_build_bond_m0_blocks`` sums over
    ALL orbital pairs, not just the diagonal ones) -- these are the
    "(m, orbital-pair) blocks" later tasks consume.

    Linearity, asserted by the caller (not here): ``W(0.5)`` elementwise
    equals ``0.5 * W(1.0)`` (``assert_approx_array(rel=0, abs=1e-14)``),
    AND the boolean supports of ``W(0.5)`` and ``W(1.0)`` are identical --
    both pinned by evaluating ``dW_matrices`` on a SEPARATE ``closed_decls``
    already scaled to the desired amplitude (this function itself takes no
    scale argument).
    """
    zero_decls = {k: 0.0 * v for k, v in closed_decls.items()}
    S1, C1 = _bare_vertices_at(fx, closed_decls)
    S0, C0 = _bare_vertices_at(fx, zero_decls)
    return S1 - S0, C1 - C0


def ed_to_solver_bond_map(channels, bond_set):
    """DERIVED index/conjugation map between ``SectorED.bond_correlator``'s
    ``(q, I, J)`` and ``bond_bubble``'s ``(q, m*nd+idx, m'*nd+idx')``.

    Derivation. ``bond_bubble``'s validated closed form (Onari Eq.14, see
    ``tests/test_bond_channels.py``'s independent ``_direct_bond_bubble``
    reference)::

        chi_bar[m,(l1,l2); m',(l3,l4)](q)
            = -(T/N) sum_k e^{ik.(dr_m - dr_m')} G_{l1,l3}(k+q) G_{l4,l2}(k)

    is, by direct momentum-space Fourier transform of the bond bilinear
    ``B_{R,a,b}(q) = sum_j e^{iqj} c^+_{j+R,a} c_{j,b} = sum_p e^{-ipR}
    c^+_{p,a} c_{p-q,b}`` and the standard free-fermion Wick contraction of
    ``<B_{R,a,b}(q); B_{R',c,d}(q)^dagger>`` (same q on both legs, matching
    ``SectorED.bond_correlator``'s own ``<B;B^dagger>`` definition), the
    SAME-SPIN-DIAGONAL entry at orbital assignment ``a=l2, b=l1, c=l4,
    d=l3`` -- an ORBITAL-INDEX SWAP relative to ``bond_bubble``'s own
    ``idx = l1*norb+l2`` labeling, invisible at ``norb=1`` (fx5, where the
    swap is a no-op) and LOAD-BEARING at ``norb=2`` (fx3, where an
    unswapped map leaves an O(1) residual, not a finite-Nmat one -- see
    the module docstring's "orientation-sensitive pin" note and Task 4's
    report). The bond-channel index ``m`` itself is NOT reversed or
    swapped, no additional complex conjugation is needed, and cross-spin
    (sigma != sigma') blocks are exactly null (verified alongside the pin
    below, not asserted here). This is stated ONCE, from the two
    definitions; the pin below VERIFIES it numerically (a finite-Nmat
    residual, never zero, since ``bond_bubble`` sums a truncated Matsubara
    grid while ``bond_correlator`` is the exact Lehmann value) -- it does
    not search over candidate maps.

    Parameters
    ----------
    channels : list
        The EXACT ordered channel list to be passed to
        ``SectorED.bond_correlator`` -- ``[(R, a, b), ...]`` (``[(R,)]`` at
        norb=1). Must contain, for every ``bond_set`` channel ``m`` and
        every physical orbital pair ``(l1, l2)``, the swapped entry
        ``(bond_set.delta_r[m][0], l2, l1)`` (typically the caller
        enumerates ALL ``(m, l1, l2)`` combinations, which contains both
        orbital orders automatically).
    bond_set : ResolvedInteractionSet

    Returns
    -------
    slot_to_channel : ndarray of int, shape (ND,)
        ``slot_to_channel[m*nd + l1*norb + l2]`` is the index into
        ``channels`` (hence into ``bond_correlator``'s I/J axis) of the
        matching ED channel.
    """
    norb = bond_set.v_bond[0].shape[0]
    nd = norb * norb
    B = int(bond_set.n_channels)
    index_of = {ch: i for i, ch in enumerate(channels)}
    slot_to_channel = np.full(B * nd, -1, dtype=int)
    for m in range(B):
        if bond_set.delta_r[m][1] != 0 or bond_set.delta_r[m][2] != 0:
            # This map (and every fixture in this module) is derived for
            # the 1-D ring convention (channels indexed by the single ring
            # displacement R = delta_r[m][0]) -- a nonzero y/z component
            # would silently alias onto the wrong channel key below rather
            # than failing loudly (rider M-3, coordinator review).
            raise ValueError(
                "ed_to_solver_bond_map: bond_set.delta_r[{}] = {} has a "
                "nonzero y/z displacement; this map only supports the "
                "1-D ring fixtures used in this module (R = delta_r[m][0])"
                .format(m, bond_set.delta_r[m]))
        R = bond_set.delta_r[m][0]
        for l1 in range(norb):
            for l2 in range(norb):
                idx = l1 * norb + l2
                key = (R,) if norb == 1 else (R, l2, l1)
                if key not in index_of:
                    raise ValueError(
                        "ed_to_solver_bond_map: channels is missing {} "
                        "(needed for bond_bubble slot m={}, (l1,l2)=({},{}))"
                        .format(key, m, l1, l2))
                slot_to_channel[m * nd + idx] = index_of[key]
    return slot_to_channel


def pred_first_order(chibar, dS, dC):
    """``(+chibar @ dS @ chibar, -chibar @ dC @ chibar)`` per q -- the
    first-order sandwich prediction, EXACT in the coupling (linear, since
    ``bare_bond_vertices`` is linear in the declared interaction); no
    Richardson stencil on this side (design doc, "the particle-particle
    identity", and ``ed_oracle_util.richardson``'s own docstring: "the
    solver-side prediction is exactly linear ... and has no Richardson
    stencil"). Signs from ``dress_bond``'s ``chi_s = solve(I - chibar S,
    chibar)`` / ``chi_c = solve(I + chibar C, chibar)`` expanded to first
    order in the coupling.
    """
    return chibar @ dS @ chibar, -(chibar @ dC @ chibar)


def _richardson_nmat(fn, n1, order):
    """Two-point Richardson extrapolation of ``fn(nmat)`` as
    ``nmat -> infinity``, assuming the leading finite-Nmat error is
    ``O(1/nmat**order)``: combining ``fn(n1)`` and ``fn(2*n1)`` cancels
    that leading term exactly, leaving the next order. ``order=1`` for a
    RAW (no tail correction) Matsubara sum -- ``bond_bubble``, which
    documents that it lacks the tail correction -- and ``order=2`` for a
    TAIL-CORRECTED one (``coeff_tail=1.0``, the same acceleration
    ``hwave.sc._calc_green``'s docstring describes). Measured empirically
    on fx5 before being used as a pin here: at ``order=1`` the raw
    ``bond_bubble`` diagonal (m=m') error is C/Nmat to 4 significant
    figures over a 16x Nmat range; at ``order=2`` the tail-corrected
    chi0q error is C/Nmat**2 over the same range -- both clean single-term
    laws, so a single two-point combination is sufficient (no 3-point
    stencil needed).
    """
    f1 = fn(n1)
    f2 = fn(2 * n1)
    p = 2 ** order
    return (p * f2 - f1) / (p - 1)


def _write_free_inputs(d, fx):
    """Write ``geom.dat``/``transfer.dat``/``coulombintra.dat`` describing
    ``fx``'s free ring Hamiltonian (``fx.t``, ``fx.eps``) as a solver input
    directory, plus a NEGLIGIBLE (``1e-9``) ``CoulombIntra`` declaration on
    orbital 0 so the solver has a nonzero interaction file to resolve
    (mirroring ``tests/test_rpa_vs_ed_oracle.py``'s ``_write_inputs``
    pattern) -- irrelevant to chi0q, which is coupling-independent by
    construction.
    """
    # %.17g round-trips a float64 exactly (17 significant digits, per
    # IEEE-754 double precision) -- review finding: %.12f/%.12e truncation
    # would inject a deterministic ~5e-13 input-rounding error ahead of
    # Pin 2a's abs=1e-12 comparison, making that pin fragile/platform-
    # sensitive rather than genuinely round-off.
    norb = fx.norb
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % norb)
        for _ in range(norb):
            f.write("0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n3\n1 1 1\n" % norb)
        for (a, b), tv in fx.t.items():
            f.write(" 1 0 0 %d %d %.17g %.17g\n" % (a + 1, b + 1, tv.real, tv.imag))
            back = np.conj(fx.t[(b, a)])
            f.write("-1 0 0 %d %d %.17g %.17g\n"
                    % (a + 1, b + 1, back.real, back.imag))
        for o, e in enumerate(fx.eps):
            f.write(" 0 0 0 %d %d %.17g 0.0\n" % (o + 1, o + 1, e))
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("hdr\n%d\n1\n1\n" % norb)
        f.write(" 0 0 0 1 1 %.17g 0.0\n" % 1e-9)
    return {"path_to_input": d, "Geometry": "geom.dat",
            "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}


def _run_chi0q(fx, nmat, coeff_tail):
    """The RPA solver's own ``chi0q`` (NOT the interaction-dressed
    ``chiq``) for ``fx``'s free ring, static (zero bosonic frequency)
    slice, shape ``(L, norb, norb, norb, norb)``. ``calc_scheme="general"``
    (the full ``(a,c,b,d)`` index -- the only scheme carrying the shape
    ``bond_bubble``'s m=0 block can be compared against).
    """
    tmp = tempfile.TemporaryDirectory()
    try:
        inter = _write_free_inputs(tmp.name, fx)
        info_mode = {"mode": "RPA",
                     "param": {"T": fx.T, "mu": fx.mu,
                               "CellShape": [fx.L, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": nmat,
                               "coeff_tail": coeff_tail},
                     "calc_scheme": "general"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": tmp.name, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        chi0 = np.asarray(gi["chi0q"])
    finally:
        tmp.cleanup()
    nf = chi0.shape[0]
    no = fx.norb
    return chi0.reshape(nf, fx.L, no, no, no, no)[nf // 2]


class TestFreeGreen(ApproxTestCase):
    """``free_green`` solves its own documented definition: ``[(iomega +
    mu) I - eps(k)] @ G = I`` at every (k, iomega_n), both fixtures. A
    direct algebraic check independent of the RPA pipeline Pin 1/2a also
    exercise ``free_green`` through."""

    def _check(self, fx, nmat=32):
        green = free_green(fx, nmat)
        self.assertEqual(green.shape, (fx.norb, fx.norb, fx.L, 1, 1, nmat))
        beta = fx.beta
        iomega = (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi / beta
        tmat = np.zeros((fx.norb, fx.norb), dtype=complex)
        for (a, b), v in fx.t.items():
            tmat[a, b] = v
        for kx in range(fx.L):
            k = 2.0 * np.pi * kx / fx.L
            eps_k = (tmat * np.exp(1j * k) + tmat.conj().T * np.exp(-1j * k)
                      + np.diag(np.asarray(fx.eps, dtype=complex)))
            for n in range(nmat):
                mat = (1j * iomega[n] + fx.mu) * np.eye(fx.norb) - eps_k
                resid = mat @ green[:, :, kx, 0, 0, n] - np.eye(fx.norb)
                assert_approx_array(resid, np.zeros_like(resid), rel=0, abs=1e-10)

    def test_fx5(self):
        self._check(fx5())

    def test_fx3(self):
        self._check(fx3())


class TestPin1TailCorrectedChi0q(ApproxTestCase):
    """Pin 1: ED ``chi_connected`` at V=0 vs the RPA solver's TAIL-
    CORRECTED ``chi0q`` (``coeff_tail=1.0``, chi0q-only mode -- NOT an
    interaction-dressed ``chiq``), ``abs=1e-10``, both fixtures. The raw
    tail-corrected error is O(1/Nmat**2) (measured, see
    ``_richardson_nmat``'s docstring) -- at the module's ``NMAT`` alone
    that is ~2e-7 on fx5, well above 1e-10, so a two-point Nmat-Richardson
    extrapolation (``NMAT``, ``2*NMAT``) is used to reach the pin's
    tolerance cheaply.

    Spin-block combos (design doc, "external-leg indexing"): the bare
    bubble is diagonal in spin for a spin-independent free Hamiltonian --
    a single closed fermion loop cannot connect an up- to a down-
    propagator -- so BOTH same-spin blocks (0,0) and (1,1) of
    ``to_solver_slots(chi_connected)`` must equal chi0q, and BOTH cross-
    spin blocks (0,1) and (1,0) must be exactly null (also asserted here,
    at the same tolerance, so a regression that starts leaking spin
    off-diagonal content would be caught).
    """

    def _check(self, fx):
        # SectorED, NOT the dense ed_oracle_util.chi_connected -- fx3 (case
        # M, dim=4096) is exponential-in-Hilbert-space for the dense path
        # (Task 4's report: 12 dense (4096, 4096) annihilators alone would
        # be multiple GB); SectorED reproduces it exactly (Task 4's own
        # round-off pin) at a fraction of the cost.
        ed = ed_oracle_util.to_solver_slots(
            ed_oracle_util.SectorED(fx).chi_connected())
        no = fx.norb
        chi0 = _richardson_nmat(
            lambda n: _run_chi0q(fx, n, coeff_tail=1.0), NMAT, order=2)
        for (s1, s2) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            with self.subTest(spins=(s1, s2)):
                block = ed[:, s1 * no:(s1 + 1) * no, s1 * no:(s1 + 1) * no,
                            s2 * no:(s2 + 1) * no, s2 * no:(s2 + 1) * no]
                expected = chi0 if s1 == s2 else np.zeros_like(chi0)
                assert_approx_array(block, expected, rel=0, abs=1e-10)

    def test_fx5(self):
        self._check(fx5())

    def test_fx3(self):
        self._check(fx3())


class TestPin2aRawChi0(ApproxTestCase):
    """Pin 2a (spec's pin 2, solver-internal, both fixtures):
    ``bond_bubble``'s Delta r=0 block vs the RPA solver's RAW chi0
    (``coeff_tail=0.0``, SAME Nmat as ``bond_bubble``, no Richardson) --
    round-off, ``abs=1e-12``: both sides are the identical finite (no
    tail correction) Matsubara sum, assembled by ``bond_bubble`` and
    ``RPA._calc_chi0q`` respectively (rider M-4, coordinator review: the
    two share the underlying ``matsubara.fermion_to_tau``/spatial-FFT
    backends and are not fully independent implementations -- what THIS
    pin actually covers is the bond-specific machinery layered on top,
    the ``m=m'=0`` phase/roll/block-extraction bookkeeping, plus
    ``free_green`` against the eigenbasis Green function the solver
    builds from ``transfer.dat``), so this pin does not replace the ED
    frame-map comparison (Pin 2b) -- it is solver-internal self-
    consistency only.
    """

    def _check(self, fx):
        # ANY nonzero declared shell works (bond_bubble's m=0 block does
        # not depend on the declared V value); reuse fx's own hopping
        # shell (R=1) purely as a nonempty two-sided declaration.
        decls = _two_sided_decl([(1, 0, 0, V1)]) if fx.norb == 1 else \
            _two_sided_decl([(1, 0, 1, V1)])
        bond_set, chi_bar = bond_setup(fx, decls, nmat=NMAT)
        nd = fx.norb * fx.norb
        m0_block = chi_bar[:, 0, 0, 0:nd, 0:nd]
        chi0_raw = _run_chi0q(fx, NMAT, coeff_tail=0.0)
        no = fx.norb
        assert_approx_array(
            m0_block.reshape(fx.L, no, no, no, no), chi0_raw, rel=0, abs=1e-12)

    def test_fx5(self):
        self._check(fx5())

    def test_fx3(self):
        self._check(fx3())


class TestEdToSolverBondMapGuard(unittest.TestCase):
    """Rider M-3 (coordinator review): exercises the nonzero-y/z-
    displacement guard directly -- every fixture in this module is a 1-D
    ring, so nothing else in the suite reaches this branch, and an
    untested rejection path could silently regress into a wrong-but-
    finite alias instead of a loud failure."""

    def test_nonzero_z_displacement_raises(self):
        decls = {((1, 0, 1), (0, 0)): 0.1, ((-1, 0, -1), (0, 0)): 0.1}
        bond_set = resolve_interactions(decls, np.eye(3), norb=1)
        with self.assertRaises(ValueError):
            ed_to_solver_bond_map([(0,)], bond_set)

    def test_nonzero_y_displacement_raises(self):
        decls = {((0, 1, 0), (0, 0)): 0.1, ((0, -1, 0), (0, 0)): 0.1}
        bond_set = resolve_interactions(decls, np.eye(3), norb=1)
        with self.assertRaises(ValueError):
            ed_to_solver_bond_map([(0,)], bond_set)


class TestPin2bFrameMap(ApproxTestCase):
    """Pin 2b (LOAD-BEARING, Task 4 review): ``ed_to_solver_bond_map``,
    zero coupling, FULL (m, m') matrix, both fixtures. ``bond_bubble``
    (RAW, no tail correction, per its own documented limitation) is
    compared against ``SectorED.bond_correlator`` (exact Lehmann) through
    the mapped index.

    TWO distances are measured, printed and recorded, and must not be
    confused (coordinator review, round 2):

    - ``eps2_rich``: the Nmat-Richardson-extrapolated distance (``order=1``,
      matching ``bond_bubble``'s raw O(1/Nmat) convergence -- see
      ``_richardson_nmat``), computed at the (``n1``, ``2*n1``) base pair
      stated per fixture below. This is the number that actually adjudicates
      the frame map's ORIENTATION (whether ``ed_to_solver_bond_map`` is
      right): an unswapped/mis-oriented candidate fails it by O(1), not by
      a finite-Nmat amount, so driving it comfortably under the fixed
      ``< 1e-6`` safety ceiling is what makes the pin discriminating. It is
      a RECORDED DIAGNOSTIC with that fixed ceiling, never an input to any
      tolerance formula, and it is NOT representative of the raw production
      path's accuracy (see below).
    - ``eps2_raw``: the distance at the single RAW ``NMAT`` (the module
      constant, 1024) with NO Richardson extrapolation -- i.e. exactly the
      quantity Tasks 6-8 consume when they call ``bond_bubble`` directly at
      production Nmat. This is ~0.2/Nmat on the m==m' (delta-r = delta-r')
      diagonal slots (the raw O(1/Nmat) law measured elsewhere in this
      module), so at NMAT=1024 it sits around 2e-4 -- FIVE TO SIX ORDERS
      larger than ``eps2_rich``. A caller sizing a later tolerance from
      ``eps2_rich`` would be far too tight; ``eps2_raw`` is the number that
      bounds the raw production path. Asserted only against a loose sanity
      ceiling (``< 1e-3``, catching a gross regression, not a precision
      target) since it is not what this pin exists to adjudicate.

    Both same-spin diagonal blocks (sigma=sigma'=0 and sigma=sigma'=1) are
    checked against the SAME mapped matrix (spin-independent free
    Hamiltonian), and the cross-spin blocks are checked to be exactly
    null -- the ORIENTATION STRUCTURE (which off-diagonal (m, m') block is
    the conjugate of which) is pinned by construction: the full ND x ND
    matrix is compared elementwise, not reduced to a scalar or a diagonal
    slice, so a swapped/transposed/conjugated off-diagonal block cannot
    hide behind an accidental cancellation.
    """

    def _check(self, fx, decls, n1, raw_nmat=NMAT):
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        norb = fx.norb
        nd = norb * norb
        channels = []
        for m in range(bond_set.n_channels):
            R = bond_set.delta_r[m][0]
            for l1 in range(norb):
                for l2 in range(norb):
                    channels.append((R,) if norb == 1 else (R, l1, l2))
        smap = ed_to_solver_bond_map(channels, bond_set)

        # memoized so a raw_nmat matching one of the Richardson samples
        # (e.g. the default raw_nmat=NMAT=2*n1) is not recomputed (review
        # round 3, optional-turned-should_fix).
        @functools.lru_cache(maxsize=None)
        def _chibar(nmat):
            green = free_green(fx, nmat)
            return bond_bubble(green, bond_set, beta=fx.beta)[:, 0, 0]

        chibar_rich = _richardson_nmat(_chibar, n1, order=1)
        chibar_raw = _chibar(raw_nmat)

        sector = ed_oracle_util.SectorED(fx)
        xph = sector.bond_correlator(channels)   # (L, 2, 2, n_i, n_i)

        def _measure(chibar, label):
            worst = 0.0
            for (s1, s2) in ((0, 0), (0, 1), (1, 0), (1, 1)):
                with self.subTest(distance=label, spins=(s1, s2)):
                    xs = xph[:, s1, s2][:, smap][:, :, smap]  # (L, ND, ND)
                    expected = chibar if s1 == s2 else np.zeros_like(chibar)
                    diff = float(np.abs(xs - expected).max())
                    # review finding: max(worst, diff) silently discards a
                    # NaN diff (NaN compares False to everything, so
                    # Python's max keeps the OTHER, possibly-zero, operand)
                    # -- which would let this load-bearing pin pass on a
                    # broken (NaN) result. Fail loudly instead.
                    self.assertTrue(
                        np.isfinite(diff),
                        "Pin2b: non-finite {} discrepancy (spins={}, "
                        "diff={})".format(label, (s1, s2), diff))
                    worst = max(worst, diff)
            return worst

        eps2_rich = _measure(chibar_rich, "rich")
        eps2_raw = _measure(chibar_raw, "raw")
        print("Pin2b (fixture L={}, norb={}): eps2_rich (Nmat base pair "
              "({}, {})) = {:.3e}; eps2_raw (single Nmat={}, the raw "
              "production-path distance) = {:.3e}"
              .format(fx.L, fx.norb, n1, 2 * n1, eps2_rich, raw_nmat, eps2_raw))
        self.assertLess(eps2_rich, 1e-6)
        self.assertLess(eps2_raw, 1e-3)

    def test_fx5(self):
        decls = _two_sided_decl([(1, 0, 0, V1), (2, 0, 0, 0.5 * V1)])
        self._check(fx5(), decls, n1=512)

    def test_fx3(self):
        decls = _two_sided_decl([(1, 0, 0, V1), (1, 0, 1, V1)])
        self._check(fx3(), decls, n1=512)


class TestRichardsonNmatConvergenceOrder(ApproxTestCase):
    """Round-2 review finding: the O(1/Nmat)/O(1/Nmat**2) leading-error
    laws ``_richardson_nmat`` relies on (docstring) were empirically
    checked against ED only on fx5 during development; Pin 1/2b apply the
    SAME assumption on fx3 (multi-orbital) too. Pin the assumption
    directly there by comparing the extrapolated value at (n1, 2*n1)
    against (2*n1, 4*n1) -- if the assumed order were wrong for fx3, these
    would disagree by O(the leading un-cancelled term), not agree tightly.
    """

    def test_pin2b_raw_bond_bubble_order1_stable_on_fx3(self):
        fx = fx3()
        decls = _two_sided_decl([(1, 0, 0, V1), (1, 0, 1, V1)])
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)

        def _chibar(nmat):
            green = free_green(fx, nmat)
            return bond_bubble(green, bond_set, beta=fx.beta)[:, 0, 0]

        r1 = _richardson_nmat(_chibar, 256, order=1)
        r2 = _richardson_nmat(_chibar, 512, order=1)
        assert_approx_array(r1, r2, rel=0, abs=1e-7)

    def test_pin1_tail_corrected_chi0q_order2_stable_on_fx3(self):
        fx = fx3()
        r1 = _richardson_nmat(
            lambda n: _run_chi0q(fx, n, coeff_tail=1.0), 256, order=2)
        r2 = _richardson_nmat(
            lambda n: _run_chi0q(fx, n, coeff_tail=1.0), 512, order=2)
        assert_approx_array(r1, r2, rel=0, abs=1e-9)


class TestPredFirstOrder(ApproxTestCase):
    """``pred_first_order`` exercised end to end (review finding: it was
    otherwise built but never called by any test): its output must match a
    finite-coupling Richardson derivative (in the coupling, NOT Nmat) of
    ``dress_bond``'s own ``chi_s``/``chi_c`` at the SAME zero-coupling
    ``chi_bar`` -- i.e. ``pred_first_order`` reproduces the exact g -> 0
    derivative of the production dressing path itself, not merely of its
    ingredients in isolation. Uses the production ``dress_bond`` directly
    (not a hand-rolled solve), so a sign/matrix-order regression in either
    ``pred_first_order`` or the ``dW_matrices`` convention it consumes
    would show up here.
    """

    def _check(self, fx, decls_at):
        bond_set0, chi_bar = bond_setup(fx, decls_at(0.0), nmat=NMAT)
        dS, dC = dW_matrices(fx, decls_at(1.0))
        # dW_matrices/pred_first_order both work in the (L, ND, ND)
        # convenience shape (Ny=Nz=1 dropped); dress_bond needs chi_bar's
        # full (Nx, Ny, Nz, ND, ND) shape, so that one keeps it.
        chi_bar_flat = chi_bar[:, 0, 0]
        pred_s, pred_c = pred_first_order(chi_bar_flat, dS, dC)

        def _chi(g):
            S_bond, C_bond = _bare_vertices_at(fx, decls_at(g))
            return dress_bond(chi_bar, S_bond[:, np.newaxis, np.newaxis],
                              C_bond[:, np.newaxis, np.newaxis])

        chi_s0, chi_c0 = _chi(0.0)
        chi_s0, chi_c0 = chi_s0[:, 0, 0], chi_c0[:, 0, 0]
        g1 = 1e-4
        chi_s1, chi_c1 = _chi(g1)
        chi_s1, chi_c1 = chi_s1[:, 0, 0], chi_c1[:, 0, 0]
        chi_s2, chi_c2 = _chi(2 * g1)
        chi_s2, chi_c2 = chi_s2[:, 0, 0], chi_c2[:, 0, 0]
        num_s = 2 * (chi_s1 - chi_s0) / g1 - (chi_s2 - chi_s0) / (2 * g1)
        num_c = 2 * (chi_c1 - chi_c0) / g1 - (chi_c2 - chi_c0) / (2 * g1)
        # measured finite-g1 residual ~1e-9/1e-10 on both fixtures; 1e-6 is
        # a safe ceiling well above that and well below the ~0.1-0.25
        # signal scale (a sign/transpose bug would fail by O(1), not O(g1)).
        assert_approx_array(pred_s, num_s, rel=0, abs=1e-6)
        assert_approx_array(pred_c, num_c, rel=0, abs=1e-6)

    def test_fx5(self):
        self._check(fx5(), lambda g: _two_sided_decl([(1, 0, 0, g)]))

    def test_fx3(self):
        # inter-orbital direction: exercises a non-diagonal dS/dC block,
        # not just the single-orbital Hartree-diagonal case above.
        self._check(fx3(), lambda g: _two_sided_decl([(1, 0, 1, g)]))


class TestDwMatrices(ApproxTestCase):
    """``dW_matrices`` element-by-element pins (brief, "dW assertions")."""

    def test_single_shell_fx5(self):
        fx = fx5()
        decls = _two_sided_decl([(1, 0, 0, 1.0)])
        dS, dC = dW_matrices(fx, decls)
        bond_set = resolve_interactions(decls, np.eye(3), norb=1)
        m_plus1 = bond_set.delta_r.index((1, 0, 0))
        m_minus1 = bond_set.delta_r.index((-1, 0, 0))
        for q in range(fx.L):
            qv = 2.0 * np.pi * q / fx.L
            expected_S = np.zeros((3, 3), dtype=complex)
            expected_S[m_plus1, m_plus1] = 1.0
            expected_S[m_minus1, m_minus1] = 1.0
            expected_C = np.zeros((3, 3), dtype=complex)
            expected_C[0, 0] = 4.0 * np.cos(qv)
            expected_C[m_plus1, m_plus1] = -1.0
            expected_C[m_minus1, m_minus1] = -1.0
            with self.subTest(q=q):
                assert_approx_array(dS[q], expected_S, rel=0, abs=1e-10)
                assert_approx_array(dC[q], expected_C, rel=0, abs=1e-10)

    def test_linearity_pin(self):
        fx = fx5()
        decls_unit = _two_sided_decl([(1, 0, 0, 1.0)])
        decls_half = _two_sided_decl([(1, 0, 0, 0.5)])
        dS1, dC1 = dW_matrices(fx, decls_unit)
        dS_h, dC_h = dW_matrices(fx, decls_half)
        assert_approx_array(dS_h, 0.5 * dS1, rel=0, abs=1e-14)
        assert_approx_array(dC_h, 0.5 * dC1, rel=0, abs=1e-14)
        self.assertTrue(np.array_equal(dS_h != 0, dS1 != 0))
        self.assertTrue(np.array_equal(dC_h != 0, dC1 != 0))

    def test_two_shell_direction_g2(self):
        fx = fx5()
        decls = _two_sided_decl([(2, 0, 0, 1.0)])
        dS, dC = dW_matrices(fx, decls)
        bond_set = resolve_interactions(decls, np.eye(3), norb=1)
        m_plus2 = bond_set.delta_r.index((2, 0, 0))
        m_minus2 = bond_set.delta_r.index((-2, 0, 0))
        for q in range(fx.L):
            qv = 2.0 * np.pi * q / fx.L
            expected_S = np.zeros((3, 3), dtype=complex)
            expected_S[m_plus2, m_plus2] = 1.0
            expected_S[m_minus2, m_minus2] = 1.0
            expected_C = np.zeros((3, 3), dtype=complex)
            expected_C[0, 0] = 4.0 * np.cos(2.0 * qv)
            expected_C[m_plus2, m_plus2] = -1.0
            expected_C[m_minus2, m_minus2] = -1.0
            with self.subTest(q=q):
                assert_approx_array(dS[q], expected_S, rel=0, abs=1e-10)
                assert_approx_array(dC[q], expected_C, rel=0, abs=1e-10)

    def test_two_shell_joint_ray_dw(self):
        """Coordinator review (Important 2): no other dW test declares TWO
        shells at once, so ``v_bond[m]``-to-``delta_r[m]`` pairing across
        ``B=5`` (m=0, +-1, +-2 all present together) is unexercised -- equal
        amplitudes or single-shell declarations cannot see a channel-
        ordering bug (e.g. the g1/g2 shells silently swapped). Declares
        BOTH shells at UNEQUAL amplitude (g1=1.0, g2=0.6) simultaneously and
        pins the exact elements, plus the superposition identity
        ``dW(joint) == dW(g1) + 0.6*dW(g2)`` elementwise (comparing full
        ``B=5`` matrices throughout: the g1-only/g2-only declarations below
        also declare the OTHER shell at amplitude 0 so all three share the
        SAME B=5 topology/index ordering -- "declared topology, not
        magnitude" per ``resolve_interactions``, matching the linearity
        pin's convention above)."""
        fx = fx5()
        decls_joint = _two_sided_decl([(1, 0, 0, 1.0), (2, 0, 0, 0.6)])
        dS, dC = dW_matrices(fx, decls_joint)
        bond_set = resolve_interactions(decls_joint, np.eye(3), norb=1)
        m_plus1 = bond_set.delta_r.index((1, 0, 0))
        m_minus1 = bond_set.delta_r.index((-1, 0, 0))
        m_plus2 = bond_set.delta_r.index((2, 0, 0))
        m_minus2 = bond_set.delta_r.index((-2, 0, 0))
        for q in range(fx.L):
            qv = 2.0 * np.pi * q / fx.L
            expected_S = np.zeros((5, 5), dtype=complex)
            expected_S[m_plus1, m_plus1] = 1.0
            expected_S[m_minus1, m_minus1] = 1.0
            expected_S[m_plus2, m_plus2] = 0.6
            expected_S[m_minus2, m_minus2] = 0.6
            expected_C = np.zeros((5, 5), dtype=complex)
            expected_C[0, 0] = 4.0 * np.cos(qv) + 2.4 * np.cos(2.0 * qv)
            expected_C[m_plus1, m_plus1] = -1.0
            expected_C[m_minus1, m_minus1] = -1.0
            expected_C[m_plus2, m_plus2] = -0.6
            expected_C[m_minus2, m_minus2] = -0.6
            with self.subTest(q=q):
                assert_approx_array(dS[q], expected_S, rel=0, abs=1e-10)
                assert_approx_array(dC[q], expected_C, rel=0, abs=1e-10)

        # superposition: g1-only/g2-only, each declared over the SAME B=5
        # topology (the other shell present at amplitude 0) so the three
        # results share one index ordering and can be added elementwise.
        decls_g1_full = _two_sided_decl([(1, 0, 0, 1.0), (2, 0, 0, 0.0)])
        decls_g2_full = _two_sided_decl([(1, 0, 0, 0.0), (2, 0, 0, 1.0)])
        dS1, dC1 = dW_matrices(fx, decls_g1_full)
        dS2, dC2 = dW_matrices(fx, decls_g2_full)
        assert_approx_array(dS, dS1 + 0.6 * dS2, rel=0, abs=1e-14)
        assert_approx_array(dC, dC1 + 0.6 * dC2, rel=0, abs=1e-14)

    def _case_m_direction(self, R_decl, a, b):
        """Shared support-pattern check for case-M's g+/g- directions:
        the bond-diagonal Fock slot at (R_decl, a, b) and its Hermitian-
        closed partner (-R_decl, b, a) each carry dS=+1/dC=-1, and the
        (aa),(bb) Hartree cross slot carries dC0 = 2*e^{-i*q*R_decl}
        (``_build_bond_m0_blocks``' ``V_ab(q) = sum_{m>=1} V_ab(dr_m)
        e^{-i q.dr_m}``, here only the ``R_decl`` shell nonzero at (a,b))."""
        fx = fx3()
        decls = _two_sided_decl([(R_decl, a, b, 1.0)])
        dS, dC = dW_matrices(fx, decls)
        bond_set = resolve_interactions(decls, np.eye(3), norb=2)
        m_fwd = bond_set.delta_r.index((R_decl, 0, 0))
        m_bwd = bond_set.delta_r.index((-R_decl, 0, 0))
        nd = 4
        I_ab = m_fwd * nd + a * 2 + b
        I_ba = m_bwd * nd + b * 2 + a
        I_aa = a * 2 + a
        I_bb = b * 2 + b
        for q in range(fx.L):
            qv = 2.0 * np.pi * q / fx.L
            expected_S = np.zeros((3 * nd, 3 * nd), dtype=complex)
            expected_S[I_ab, I_ab] = 1.0
            expected_S[I_ba, I_ba] = 1.0
            expected_C = np.zeros((3 * nd, 3 * nd), dtype=complex)
            expected_C[I_ab, I_ab] = -1.0
            expected_C[I_ba, I_ba] = -1.0
            # C0 is Hermitian (bond_set's reversal closure gives
            # V_q[b,a] = conj(V_q[a,b])), so the (bb,aa) cross slot carries
            # the conjugate of (aa,bb), not zero.
            expected_C[I_aa, I_bb] = 2.0 * np.exp(-1j * qv * R_decl)
            expected_C[I_bb, I_aa] = 2.0 * np.exp(1j * qv * R_decl)
            with self.subTest(q=q):
                assert_approx_array(dS[q], expected_S, rel=0, abs=1e-10)
                assert_approx_array(dC[q], expected_C, rel=0, abs=1e-10)

    def test_case_m_g_plus(self):
        # g+ = [V_01(+1) = V_10(-1)]
        self._case_m_direction(1, 0, 1)

    def test_case_m_g_minus(self):
        # g- = [V_01(-1) = V_10(+1)]
        self._case_m_direction(-1, 0, 1)

    def test_fx5_onsite_U_direction(self):
        """Task 6's U (on-site) direction, pinned at the ELEMENT level
        (review should_fix): unlike g1/g2 above, which already have
        single-shell pins, the on-site U direction's dW had no element-
        level pin on the ACTUAL topology Task 6's granules use (fx5's
        shared B=5 grid, with both V shells declared at zero -- see that
        task's module note). Confirms the module note's code-derivation
        (S0 = U, C0 = U, touching ONLY the m=0 sub-block) against the real
        ``bare_bond_vertices`` output directly, not just the reading of
        ``_build_bond_m0_blocks`` the note describes. ``_decls_U`` is
        defined later in this module (Task 6's section) but resolved at
        call time, after the whole module has been imported."""
        fx = fx5()
        decls = _decls_U(fx, 1.0)
        dS, dC = dW_matrices(fx, decls)
        expected_S = np.zeros((5, 5), dtype=complex)
        expected_S[0, 0] = 1.0
        expected_C = np.zeros((5, 5), dtype=complex)
        expected_C[0, 0] = 1.0
        for q in range(fx.L):
            with self.subTest(q=q):
                assert_approx_array(dS[q], expected_S, rel=0, abs=1e-10)
                assert_approx_array(dC[q], expected_C, rel=0, abs=1e-10)


# ---------------------------------------------------------------------------
# Task 6: ph adjudication -- verdict machinery (Step 1: adjudicate_granule
# unit tests on synthetic arrays)
# ---------------------------------------------------------------------------

class TestAdjudicateGranule(ApproxTestCase):
    """Unit tests for ``ed_oracle_util.adjudicate_granule`` on synthetic
    (non-physical) arrays -- shape (L=1, ND=2, ND=2), matching the shape the
    real granules below actually use (q, external leg I, external leg J).
    """

    def _synthetic(self, bearing_value=1.0, zero_mask=None):
        if zero_mask is None:
            zero_mask = np.array([[True, False], [False, False]])
        base = np.zeros((1, 2, 2), dtype=complex)
        base[0, 1, 1] = bearing_value
        return zero_mask, base

    def test_pass(self):
        zero_mask, D = self._synthetic(bearing_value=1.0)
        rec = ed_oracle_util.adjudicate_granule(D, D, D, D, zero_mask, "t/pass")
        self.assertEqual(rec["status"], "PASS")
        self.assertEqual(rec["failures"], [])
        # review finding (should_fix): assertAlmostEqual's default
        # decimal-place semantics would also accept a spuriously-zero tol
        # here (7 places rounds both 1e-12 and 0.0 to 0.0000000) -- pin the
        # campaign's own approx_util rule (rel=0, an EXACT abs floor)
        # instead, matching this module's convention everywhere else.
        self.assertApprox(rec["tol"], 1e-12, rel=0, abs=1e-20)
        self.assertGreaterEqual(rec["max_signal"], 100.0 * rec["tol"])

    def test_fail_bearing_mismatch(self):
        zero_mask, D_ed = self._synthetic(bearing_value=1.0)
        D_pred = D_ed.copy()
        D_pred[0, 1, 1] = 1.5   # off by 0.5 -- far above the round-off tol
        rec = ed_oracle_util.adjudicate_granule(
            D_ed, D_ed, D_pred, D_pred, zero_mask, "t/fail")
        self.assertEqual(rec["status"], "FAIL")
        kinds = {f[0] for f in rec["failures"]}
        self.assertEqual(kinds, {"bearing"})
        self.assertEqual(rec["failures"][0][1], (0, 1, 1))

    def test_zero_cell_violation(self):
        zero_mask, D_ed = self._synthetic(bearing_value=1.0)
        D_pred = D_ed.copy()
        D_pred[0, 0, 0] = 5.0   # zero_mask[0,0] is True: this must be ~0
        rec = ed_oracle_util.adjudicate_granule(
            D_ed, D_ed, D_pred, D_pred, zero_mask, "t/zero-violation")
        self.assertEqual(rec["status"], "FAIL")
        kinds = {f[0] for f in rec["failures"]}
        self.assertIn("zero", kinds)
        zero_failure = [f for f in rec["failures"] if f[0] == "zero"][0]
        self.assertEqual(zero_failure[1], (0, 0, 0))

    def test_inconclusive_small_signal(self):
        zero_mask, D_ed_vhalf = self._synthetic(bearing_value=0.5)
        D_ed_v1 = D_ed_vhalf.copy()
        D_ed_v1[0, 1, 1] += 1e-2   # delta_rich = 1e-2 -> tol = 1e-1
        D_pred = D_ed_vhalf.copy()  # solver matches the finer ED value exactly
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred, D_pred, zero_mask, "t/inconclusive")
        self.assertApprox(rec["tol"], 1e-1, rel=0, abs=1e-15)
        self.assertEqual(rec["failures"], [])
        self.assertLess(rec["max_signal"], 100.0 * rec["tol"])
        self.assertEqual(rec["status"], "INCONCLUSIVE")

    def test_1e13_floor(self):
        zero_mask, D = self._synthetic(bearing_value=1.0)
        rec = ed_oracle_util.adjudicate_granule(D, D, D, D, zero_mask, "t/floor")
        # delta_rich = delta_nmat = 0.0 exactly -> the absolute floor alone
        # sets tol, not zero.
        self.assertEqual(rec["delta_rich"], 0.0)
        self.assertEqual(rec["delta_nmat"], 0.0)
        self.assertEqual(rec["tol"], 10.0 * 1e-13)

    def test_pass_zero_granule(self):
        zero_mask = np.ones((2, 2), dtype=bool)
        D = np.zeros((1, 2, 2), dtype=complex)
        rec = ed_oracle_util.adjudicate_granule(D, D, D, D, zero_mask, "t/zero-only")
        self.assertEqual(rec["status"], "PASS-ZERO")
        self.assertEqual(rec["max_signal"], 0.0)
        self.assertEqual(rec["failures"], [])

    def test_never_mutates_globals(self):
        zero_mask, D = self._synthetic(bearing_value=1.0)
        D_ed_v1 = D.copy()
        D_ed_vhalf = D.copy()
        D_pred_nmat = D.copy()
        D_pred_2nmat = D.copy()
        snapshots = [a.copy() for a in
                     (D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask)]
        ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask, "t/pure")
        for arr, snap in zip(
                (D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask),
                snapshots):
            self.assertTrue(np.array_equal(arr, snap))
        # calling twice with the same inputs gives the same record (no
        # hidden state accumulated across calls).
        rec1 = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask, "t/pure")
        rec2 = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask, "t/pure")
        self.assertEqual(rec1["status"], rec2["status"])
        self.assertEqual(rec1["tol"], rec2["tol"])

    def test_nan_input_raises_rather_than_false_pass_zero(self):
        # Review finding (must_fix): NaN compares False against every
        # threshold in adjudicate_granule (|x| > tol, |x| <= tol alike), so
        # a broken upstream computation could otherwise silently produce a
        # PASS-ZERO verdict with an empty failures list. Confirmed to
        # reproduce before the np.isfinite guard was added; now must raise.
        zero_mask = np.ones((2, 2), dtype=bool)
        D = np.full((1, 2, 2), np.nan, dtype=complex)
        zeros = np.zeros((1, 2, 2), dtype=complex)
        with self.assertRaises(ValueError):
            ed_oracle_util.adjudicate_granule(D, zeros, zeros, zeros, zero_mask, "t/nan")

    def test_inf_input_raises(self):
        zero_mask = np.array([[True, False], [False, False]])
        D_ed = np.zeros((1, 2, 2), dtype=complex)
        D_pred = D_ed.copy()
        D_pred[0, 1, 1] = np.inf
        with self.assertRaises(ValueError):
            ed_oracle_util.adjudicate_granule(
                D_ed, D_ed, D_pred, D_pred, zero_mask, "t/inf")


class TestProjectedStructuralZeroMask(unittest.TestCase):
    """``projected_structural_zero_mask`` on small synthetic supports:
    dense ``free_support`` collapses the propagated support to all-True
    whenever the effective vertex has ANY nonzero entry (matching this
    task's ph granules below, all of which are dense-bubble fixtures), and
    to all-False (zero-only) only when the effective vertex is entirely
    empty."""

    def test_dense_free_any_vertex_entry_gives_all_bearing(self):
        free = np.ones((3, 3), dtype=bool)
        eff = np.zeros((3, 3), dtype=bool)
        eff[1, 1] = True
        zero_mask = ed_oracle_util.projected_structural_zero_mask(free, eff)
        self.assertFalse(zero_mask.any())

    def test_empty_vertex_gives_zero_only(self):
        free = np.ones((3, 3), dtype=bool)
        eff = np.zeros((3, 3), dtype=bool)
        zero_mask = ed_oracle_util.projected_structural_zero_mask(free, eff)
        self.assertTrue(zero_mask.all())

    def test_shape_mismatch_raises(self):
        with self.assertRaises(ValueError):
            ed_oracle_util.projected_structural_zero_mask(
                np.ones((2, 2), dtype=bool), np.ones((3, 3), dtype=bool))

    def test_numeric_input_raises_rather_than_silently_casting(self):
        # Review finding (should_fix): dtype=bool on a NUMERIC array would
        # silently treat every nonzero float as True -- accidentally
        # "working" but defeating the documented "never thresholded
        # numeric matrix" contract. A caller that passes a raw dW matrix by
        # mistake (instead of its boolean support) must be rejected, not
        # silently coerced.
        free = np.ones((3, 3), dtype=bool)
        numeric_eff = np.zeros((3, 3))
        numeric_eff[1, 1] = 0.02
        with self.assertRaises(ValueError):
            ed_oracle_util.projected_structural_zero_mask(free, numeric_eff)

    def test_non_square_input_raises(self):
        with self.assertRaises(ValueError):
            ed_oracle_util.projected_structural_zero_mask(
                np.ones((2, 3), dtype=bool), np.ones((2, 3), dtype=bool))


# ---------------------------------------------------------------------------
# Task 6: ph adjudication -- single-orbital unit directions (U, g1, g2) on
# fx5, THE CAMPAIGN'S FIRST ADJUDICATION CELLS (#151). Each direction is an
# INDEPENDENT physical coupling (design doc: "Parameter basis, defined on
# the physical quotient"), adjudicated as its own (S, C) granule pair, plus
# two joint-ray superposition checks (U+V and the two-shell ray) that are
# NOT granules and NOT sensitivity inputs.
#
# U control note (brief): U is ON-SITE (Delta r = 0); its bond topology
# needs a declared V topology to exist inside the bond machinery (an
# onsite-only ``coulomb_inter`` dict resolves fine on its own -- B=1 -- but
# would leave U's dW living in a DIFFERENT, smaller ND index space than
# g1/g2's B=5 fx5 topology). So every direction below declares the SAME
# two-shell (Delta r = +-1, +-2) topology, with the inactive shell(s) at
# ZERO amplitude for U (and one of the two shells at zero for g1/g2,
# matching Task 5's ``test_two_shell_joint_ray_dw`` "declared topology, not
# magnitude" pattern) -- so dS/dC/chi_bar share one (L, 5, 5) grid across
# all three directions, which is what later alignment (Task 9's
# sensitivity_rank) needs.
#
# Verified directly against ``_build_bond_m0_blocks``/``bare_bond_vertices``
# (not guessed): for norb=1, an on-site CoulombInter declaration V_00(0)=U
# routes through ``bond_set.v_onsite`` into ``_build_bond_m0_blocks``'s
# ``V_r0``, giving Hartree ``C0 += 2*V_r0`` and Fock ``S0 += V_r0``,
# ``C0 -= V_r0`` on the SAME (aa,bb)=(0,0) slot -- i.e. S0 = U, C0 = 2U-U
# = U (the Kuroki S/C matrices), and ``bare_bond_vertices`` places this
# ENTIRELY in the m=0 sub-block (``S_bond[...,0:nd,0:nd] = S0_q``); no other
# ND slot is touched, so U's dW support is EXACTLY (m=0, m=0)-only for both
# S and C.
# ---------------------------------------------------------------------------

def _decls_U(fx, v):
    """U (on-site) direction at amplitude ``v``, carried on fx5's shared
    B=5 topology (the other two shells declared at zero amplitude -- see
    the module note above)."""
    decl = _two_sided_decl([(1, 0, 0, 0.0), (2, 0, 0, 0.0)])
    decl[((0, 0, 0), (0, 0))] = v
    return decl


def _terms_U(fx, v):
    return ed_oracle_util.canonical_density_terms(fx, [(0, 0, 0, v)])


def _decls_g1(fx, v):
    """g1 (Delta r = +-1) direction at amplitude ``v``, on the shared B=5
    topology (Delta r = +-2 declared at zero)."""
    return _two_sided_decl([(1, 0, 0, v), (2, 0, 0, 0.0)])


def _terms_g1(fx, v):
    return ed_oracle_util.canonical_density_terms(fx, [(0, 0, 1, v)])


def _decls_g2(fx, v):
    """g2 (Delta r = +-2) direction at amplitude ``v``, on the shared B=5
    topology (Delta r = +-1 declared at zero)."""
    return _two_sided_decl([(1, 0, 0, 0.0), (2, 0, 0, v)])


def _terms_g2(fx, v):
    return ed_oracle_util.canonical_density_terms(fx, [(0, 0, 2, v)])


_UNIT_DIRECTIONS = {
    "U": (_decls_U, _terms_U),
    "g1": (_decls_g1, _terms_g1),
    "g2": (_decls_g2, _terms_g2),
}


@functools.lru_cache(maxsize=None)
def _bond_topology_and_map(fx):
    """(bond_set, channels, smap) shared by every unit direction on ``fx``:
    all of them declare the SAME two-shell topology (see the module note
    above), so the bond-channel enumeration and the derived
    ``ed_to_solver_bond_map`` are identical across U/g1/g2 and can be built
    once."""
    decls = _two_sided_decl([(1, 0, 0, 1.0), (2, 0, 0, 1.0)])
    bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
    norb = fx.norb
    channels = []
    for m in range(bond_set.n_channels):
        R = bond_set.delta_r[m][0]
        for l1 in range(norb):
            for l2 in range(norb):
                channels.append((R,) if norb == 1 else (R, l1, l2))
    smap = ed_to_solver_bond_map(channels, bond_set)
    return bond_set, channels, smap


@functools.lru_cache(maxsize=None)
def _shared_chibar(fx, nmat):
    """The bare bond bubble at ``nmat``, shared across every direction on
    ``fx`` (chi_bar depends only on the topology and the free Green
    function, never on a direction's declared amplitude)."""
    bond_set, _channels, _smap = _bond_topology_and_map(fx)
    green = free_green(fx, nmat)
    return bond_bubble(green, bond_set, beta=fx.beta)[:, 0, 0]


def _dw_support_masks(direction, bond_set):
    """The ANALYTIC (boolean, not thresholded-numeric) support of dS/dC
    for one fx5 unit direction, derived from the same element pattern
    Task 5's ``dW_matrices`` element pins established (single-orbital: nd=1,
    so the enlarged index I equals the bond-channel index m directly).
    U touches ONLY the (m=0, m=0) slot (see the module note above, verified
    against the builder). g1/g2 touch the m=+-R bond-diagonal Fock slots
    (S and C) PLUS the (m=0, m=0) Hartree slot (C only -- S carries no
    Hartree term, per ``bare_bond_vertices``' documented m=0 sub-block =
    S0_q/C0_q, m!=0 diagonal = +-V_ab(R))."""
    B = bond_set.n_channels
    S_supp = np.zeros((B, B), dtype=bool)
    C_supp = np.zeros((B, B), dtype=bool)
    C_supp[0, 0] = True   # every direction here carries SOME (0,0) charge content
    if direction == "U":
        S_supp[0, 0] = True
        return S_supp, C_supp
    R = {"g1": 1, "g2": 2}[direction]
    m_plus = bond_set.delta_r.index((R, 0, 0))
    m_minus = bond_set.delta_r.index((-R, 0, 0))
    for m in (m_plus, m_minus):
        S_supp[m, m] = True
        C_supp[m, m] = True
    return S_supp, C_supp


class TestDwSupportMasksMatchNumericDw(unittest.TestCase):
    """Review should-fix: ``_dw_support_masks``' hand-derived ANALYTIC
    pattern is cross-checked here against the boolean nonzero pattern of
    the ACTUAL numeric ``dS``/``dC`` (``dW_matrices``, ANY q, OR-reduced --
    a q-dependent Hartree coefficient like ``4*cos(qR)`` is structurally
    nonzero even where it happens to vanish at one q). This does not
    replace the analytic derivation (the whole point of ``_dw_support_masks``
    is to be boolean-and-not-thresholded, never fed a numeric epsilon), but
    it does catch a mismatch between the hand-built pattern and what the
    production builder actually emits for each of the three unit
    directions this task adjudicates."""

    def _check(self, direction):
        fx = fx5()
        bond_set, _channels, _smap = _bond_topology_and_map(fx)
        decls_at, _terms_at = _UNIT_DIRECTIONS[direction]
        dS, dC = dW_matrices(fx, decls_at(fx, 1.0))
        numeric_S = np.any(np.abs(dS) > 0, axis=0)
        numeric_C = np.any(np.abs(dC) > 0, axis=0)
        S_supp, C_supp = _dw_support_masks(direction, bond_set)
        self.assertTrue(np.array_equal(numeric_S, S_supp),
                         "direction={} S support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_S, S_supp))
        self.assertTrue(np.array_equal(numeric_C, C_supp),
                         "direction={} C support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_C, C_supp))

    def test_U(self):
        self._check("U")

    def test_g1(self):
        self._check("g1")

    def test_g2(self):
        self._check("g2")


@functools.lru_cache(maxsize=None)
def _unit_direction_raw(direction):
    """Raw (un-adjudicated) ED- and solver-side first-order response
    arrays for one fx5 unit direction, S and C channels: the ED side is
    HF-SUBTRACTED (design doc: "the same HF subtraction (the derivative of
    the free bond-bilinear correlator with the HF one-body matrix inserted
    on both legs)") -- ``full(v) = SectorED(fx, terms=terms(v))
    .bond_correlator(channels)`` minus ``hf_only(v) = SectorED(fx,
    terms=(), h1=hf_h1_from_terms(fx, terms(v))).bond_correlator(channels)``
    -- so the compared quantity is the genuine VERTEX correction with the
    trivial mean-field propagator renormalization removed, matching
    ``pred_first_order``'s fixed-chi_bar sandwich exactly (both sides are
    zero at v=0 by construction: at v=0, terms=() so full == hf_only).
    Spin combos X_S = X^uu - X^ud, X_C = X^uu + X^ud are formed AT THE CALL
    SITE from the sigma-resolved HF-subtracted blocks, mapped into the
    solver's (q, ND, ND) frame via the Task-5 ``ed_to_solver_bond_map``.

    Solver side: ``pred_first_order`` (exact in the coupling) at Nmat and
    2*Nmat, using dS/dC evaluated at UNIT declared amplitude (the SLOPE,
    not a finite response) so its units match the ED-side Richardson
    derivative directly.

    Shared by ``case_fx5_unit_direction`` (wraps this with
    ``adjudicate_granule``) and the joint-ray superposition checks (which
    reuse a unit direction's OWN authoritative D_ed = D_ed_vhalf without
    recomputation, per the brief). lru_cache'd: Task 9 aggregates without
    recomputation, no test-order dependence.
    """
    fx = fx5()
    decls_at, terms_at = _UNIT_DIRECTIONS[direction]
    bond_set, channels, smap = _bond_topology_and_map(fx)

    # -- solver side --------------------------------------------------
    dS, dC = dW_matrices(fx, decls_at(fx, 1.0))
    chibar_n = _shared_chibar(fx, NMAT)
    chibar_2n = _shared_chibar(fx, 2 * NMAT)
    D_pred_nmat_S, D_pred_nmat_C = pred_first_order(chibar_n, dS, dC)
    D_pred_2nmat_S, D_pred_2nmat_C = pred_first_order(chibar_2n, dS, dC)

    # -- ED side --------------------------------------------------------
    @functools.lru_cache(maxsize=None)
    def _chi_hf_sub(v):
        terms = terms_at(fx, v)
        full = ed_oracle_util.SectorED(fx, terms=terms).bond_correlator(channels)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).bond_correlator(channels)
        return full - hf_only

    def _mapped(v, sigma, sigmap):
        chi = _chi_hf_sub(v)
        return chi[:, sigma, sigmap][:, smap][:, :, smap]

    def X_S_of_v(v):
        return _mapped(v, 0, 0) - _mapped(v, 0, 1)

    def X_C_of_v(v):
        return _mapped(v, 0, 0) + _mapped(v, 0, 1)

    D_ed_v1_S = ed_oracle_util.richardson(X_S_of_v, CAMPAIGN_V1)
    D_ed_vhalf_S = ed_oracle_util.richardson(X_S_of_v, CAMPAIGN_V1 / 2)
    D_ed_v1_C = ed_oracle_util.richardson(X_C_of_v, CAMPAIGN_V1)
    D_ed_vhalf_C = ed_oracle_util.richardson(X_C_of_v, CAMPAIGN_V1 / 2)

    # -- structural zero mask --------------------------------------------
    S_supp, C_supp = _dw_support_masks(direction, bond_set)
    free_support = np.ones_like(S_supp)   # dense free-fermion bubble (stated)
    zero_mask_S = ed_oracle_util.projected_structural_zero_mask(
        free_support, S_supp)
    zero_mask_C = ed_oracle_util.projected_structural_zero_mask(
        free_support, C_supp)

    return dict(
        D_ed_v1_S=D_ed_v1_S, D_ed_vhalf_S=D_ed_vhalf_S,
        D_pred_nmat_S=D_pred_nmat_S, D_pred_2nmat_S=D_pred_2nmat_S,
        D_ed_v1_C=D_ed_v1_C, D_ed_vhalf_C=D_ed_vhalf_C,
        D_pred_nmat_C=D_pred_nmat_C, D_pred_2nmat_C=D_pred_2nmat_C,
        zero_mask_S=zero_mask_S, zero_mask_C=zero_mask_C,
    )


@functools.lru_cache(maxsize=None)
def case_fx5_unit_direction(direction):
    """ph adjudication (S and C granules) for one independent unit coupling
    direction on fx5 (brief Step 2: ``U`` on-site, ``g1`` = Delta r = +-1,
    ``g2`` = Delta r = +-2). Returns ``dict(S=<record>, C=<record>)``, each
    the canonical full-grid record from ``adjudicate_granule``.
    ``functools.lru_cache``'d module function: Task 9 aggregates the SAME
    records without recomputation or test-order dependence (brief's
    explicit instruction)."""
    raw = _unit_direction_raw(direction)
    rec_S = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_S"], raw["D_ed_vhalf_S"],
        raw["D_pred_nmat_S"], raw["D_pred_2nmat_S"],
        raw["zero_mask_S"], "fx5/{}/S".format(direction))
    rec_C = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_C"], raw["D_ed_vhalf_C"],
        raw["D_pred_nmat_C"], raw["D_pred_2nmat_C"],
        raw["zero_mask_C"], "fx5/{}/C".format(direction))
    return {"S": rec_S, "C": rec_C}


class TestPhUnitDirectionsFx5(unittest.TestCase):
    """Step 2: the campaign's first actual verdict cells (#151 Task 6) --
    independent unit directions U / g1 / g2 on fx5, ph channels S and C.
    THE VERDICT MATTERS MORE THAN GREEN (brief): a FAIL is kept asserting
    (not weakened), recording a genuine production finding; see this
    task's report for the measured audit tuples."""

    def _check(self, direction):
        rec = case_fx5_unit_direction(direction)
        for channel in ("S", "C"):
            with self.subTest(direction=direction, channel=channel):
                r = rec[channel]
                self.assertIn(
                    r["status"], ("PASS", "PASS-ZERO"),
                    "granule {} status={} (delta_rich={:.3e} delta_nmat={:.3e} "
                    "tol={:.3e} max_signal={:.3e} first_failures={})".format(
                        r["label"], r["status"], r["delta_rich"],
                        r["delta_nmat"], r["tol"], r["max_signal"],
                        r["failures"][:5]))

    def test_U(self):
        self._check("U")

    def test_g1(self):
        self._check("g1")

    def test_g2(self):
        self._check("g2")


# ---------------------------------------------------------------------------
# Task 6: joint-ray superposition checks (U+V at (1,1), two-shell at
# (1, 0.6)) -- NOT granules, NOT sensitivity inputs (brief). Pure ED-side
# linearity checks: the joint ray's OWN Richardson derivative must equal
# the weighted SUM of the already-measured unit directions' own D_ed.
# ---------------------------------------------------------------------------

def _terms_ray_U_g1(fx, t):
    return ed_oracle_util.canonical_density_terms(
        fx, [(0, 0, 0, t), (0, 0, 1, t)])


def _terms_ray_g1_g2(fx, t):
    return ed_oracle_util.canonical_density_terms(
        fx, [(0, 0, 1, t), (0, 0, 2, 0.6 * t)])


def _joint_ray_superposition_check(testcase, directions_alphas, terms_at_ray,
                                    ray_label):
    fx = fx5()
    _bond_set, channels, smap = _bond_topology_and_map(fx)

    @functools.lru_cache(maxsize=None)
    def _chi_hf_sub(t):
        terms = terms_at_ray(fx, t)
        full = ed_oracle_util.SectorED(fx, terms=terms).bond_correlator(channels)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).bond_correlator(channels)
        return full - hf_only

    def _mapped(t, sigma, sigmap):
        chi = _chi_hf_sub(t)
        return chi[:, sigma, sigmap][:, smap][:, :, smap]

    def X_S_of_t(t):
        return _mapped(t, 0, 0) - _mapped(t, 0, 1)

    def X_C_of_t(t):
        return _mapped(t, 0, 0) + _mapped(t, 0, 1)

    for chan_name, X_of_t in (("S", X_S_of_t), ("C", X_C_of_t)):
        D_joint_v1 = ed_oracle_util.richardson(X_of_t, CAMPAIGN_V1)
        D_joint_vhalf = ed_oracle_util.richardson(X_of_t, CAMPAIGN_V1 / 2)
        delta_rich_joint = float(np.max(np.abs(D_joint_v1 - D_joint_vhalf)))

        unit_tols = [case_fx5_unit_direction(d)[chan_name]["tol"]
                     for d, _alpha in directions_alphas]
        tol_joint = 10.0 * max(delta_rich_joint, max(unit_tols))

        D_sum = sum(
            alpha * _unit_direction_raw(d)["D_ed_vhalf_" + chan_name]
            for d, alpha in directions_alphas)

        diff = float(np.max(np.abs(D_joint_vhalf - D_sum)))
        print("joint[{}/{}]: delta_rich_joint={:.6e} tol_joint={:.6e} "
              "diff={:.6e}".format(ray_label, chan_name, delta_rich_joint,
                                    tol_joint, diff))
        with testcase.subTest(ray=ray_label, channel=chan_name):
            testcase.assertLessEqual(diff, tol_joint)


class TestJointRaySuperposition(unittest.TestCase):
    """Step 2: the two joint-ray superposition checks (design doc: "at
    first order d/dV along the joint ray is the SUM of the separate
    derivatives -- no cross-term exists at this order")."""

    def test_U_plus_V_at_1_1(self):
        _joint_ray_superposition_check(
            self, [("U", 1.0), ("g1", 1.0)], _terms_ray_U_g1, "U+V(1,1)")

    def test_two_shell_ray_at_1_0p6(self):
        _joint_ray_superposition_check(
            self, [("g1", 1.0), ("g2", 0.6)], _terms_ray_g1_g2,
            "g1+0.6*g2(1,0.6)")


if __name__ == "__main__":
    unittest.main()
