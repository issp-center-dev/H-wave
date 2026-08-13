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
    _g2_from_green, bare_bond_vertices, bond_bubble, dress_bond,
    resolve_interactions,
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


# ---------------------------------------------------------------------------
# Task 7: pp reference path and adjudication (#151)
# ---------------------------------------------------------------------------
#
# The acceptance identity (design doc, "The particle-particle identity is
# DERIVED, not assumed"):
#
#     dChi_pp/dV|0 = -1/2 * X0_pp * (dV^pp/dV) * X0_pp
#
# in the ORDERED-PAIR matrix basis {Delta_i} = {Delta_{R,a,b}(q=0)} (the
# total Cooper-pair momentum is fixed at zero -- Vpp is q-independent, and
# the identity carries no q at all). production Vpp from
# bare_bond_vertices; X0_pp = production_pair_bubble; the ED side is
# richardson() of SectorED.pair_correlator with the pp HF insertion (free
# pair matrix derivative with hf_h1_from_terms on both legs), projected on
# BOTH sides with exchange_projector's Q (never bond_channels' internal
# P/Q_s/Q_t) before comparison, per channel s/t.
#
# Derivation (verified numerically against SectorED.pair_correlator at V=0,
# not guessed -- see the task report for the two independent cross-checks):
# for channels = [(R, a, b), ...] (the full ordered-pair triple, always --
# never the norb=1 (R,) shorthand),
#
#     X0_pp[(R,a,b),(R',c,d)] = (1/L) sum_k e^{ik(R'-R)} G2[a,c,b,d](k)
#
# with G2 = _g2_from_green(free_green(fx, nmat), fx.beta, tail=False) --
# mirroring bare_bond_vertices' A_coef/Bcoef contraction structure inside
# make_bond_kernel's instantaneous part (spec S4.5): A_coef[m,cd] =
# sum_{m',ef} X0_pp[(m,cd),(m',ef)] phi^{(m')}_{ef} is exactly the same
# object as bond_channels._bond_kernel_operators builds by hand for a
# specific phi. Both this formula and an INDEPENDENTLY derived one (direct
# Wick contraction of Delta_{R,a,b}(0) against Delta_{R',c,d}(0)^dagger,
# Fourier-transformed with free_green's OWN k-convention) agree exactly,
# and both match SectorED.pair_correlator at V=0 to the expected O(1/Nmat)
# finite-window residual (pin 3b).


def production_pair_bubble(fx, channels, nmat):
    """X0_pp assembled per the PRODUCTION conventions: from ``free_green``
    via the SAME finite-Nmat pp bubble ``_g2_from_green`` builds for
    ``make_bond_kernel``'s instantaneous part (``g2_tail=False``, matching
    the raw/no-tail convention ``bond_bubble`` itself uses for the ph
    channel), reconstructed into the ordered-pair basis via the phase
    ``e^{ik(R'-R)}`` (see the module note above for the derivation and its
    numeric cross-check). This is the matrix that multiplies ``dVpp`` in the
    acceptance identity; ``channels`` is the full ``(R, a, b)`` ordered-pair
    triple list (see ``exchange_projector``), never the norb=1 ``(R,)``
    shorthand.
    """
    # Review finding (must_fix): matches eigenbasis_pair_bubble's own guard
    # -- an odd/zero/negative nmat is not a fermionic Matsubara grid, and
    # pin 3a (round-off agreement with eigenbasis_pair_bubble) could not
    # catch either side silently computing the wrong quantity from the same
    # broken grid.
    # Round-2 review finding (should_fix): matches eigenbasis_pair_bubble's
    # own fix -- normalize through float() so NaN/inf/None/non-numeric
    # values raise the intended ValueError rather than a TypeError/
    # OverflowError from a bare int() call.
    try:
        nmat_f = float(nmat)
    except (TypeError, ValueError):
        raise ValueError(
            "production_pair_bubble: nmat must be a positive even integer "
            "(the centered fermionic Matsubara grid) -- got {!r}".format(nmat))
    if (not np.isfinite(nmat_f) or nmat_f != int(nmat_f) or nmat_f <= 0
            or int(nmat_f) % 2 != 0):
        raise ValueError(
            "production_pair_bubble: nmat must be a positive even integer "
            "(the centered fermionic Matsubara grid) -- got {!r}".format(nmat))
    green = free_green(fx, nmat)
    G2 = _g2_from_green(green, fx.beta, tail=False)[:, :, :, :, :, 0, 0]
    L = fx.L
    kx = 2.0 * np.pi * np.arange(L) / L
    channels = list(channels)
    n = len(channels)
    X = np.zeros((n, n), dtype=complex)
    for i, (R, a, b) in enumerate(channels):
        for j, (Rp, c, d) in enumerate(channels):
            phase = np.exp(1j * kx * (Rp - R))
            X[i, j] = np.sum(phase * G2[a, c, b, d, :]) / L
    return X


def _ordered_pair_channels(bond_set, norb):
    """The full ``(R, l1, l2)`` ordered-pair channel list -- ALWAYS the
    3-tuple, for every ``(bond channel, orbital pair)`` combination of
    ``bond_set`` -- what ``exchange_projector``/``production_pair_bubble``/
    ``eigenbasis_pair_bubble`` all need (unlike ``_bond_topology_and_map``'s
    ph channel list, which collapses to ``(R,)`` at norb=1)."""
    channels = []
    for m in range(bond_set.n_channels):
        R = bond_set.delta_r[m][0]
        for l1 in range(norb):
            for l2 in range(norb):
                channels.append((R, l1, l2))
    return channels


def _to_sector_channels(fx, channels):
    """``channels`` (full ``(R,a,b)`` triples) converted to
    ``SectorED.pair_correlator``'s own convention -- ``(R,)`` at norb=1,
    the same triple otherwise."""
    if fx.norb == 1:
        return [(R,) for (R, a, b) in channels]
    return list(channels)


def dVpp_matrices(fx, closed_decls, channels):
    """``(dVpp_s, dVpp_t)``: ``bare_bond_vertices``' bare Cooper vertices at
    the (two-sided, closed) ``closed_decls`` declaration minus at the SAME
    declared topology scaled to zero (mirrors ``dW_matrices``' same-shape
    subtraction pattern, but for ``Vpp_s``/``Vpp_t`` -- both are
    q-INDEPENDENT, built from the ``q=0`` slice only, per
    ``bare_bond_vertices``' own local-block construction, so ``S0_q``/``C0_q``
    here are evaluated at a single ``q=(0,0,0)`` point via
    ``_build_bond_m0_blocks``, not the full ring grid), REINDEXED into the
    DIRECT ordered-pair basis ``channels`` -- the SAME literal ``(R,a,b)``
    basis ``production_pair_bubble``/``exchange_projector``/ED all share
    (design doc: "``dV^pp/dV`` is the production tensor's derivative in
    that SAME basis" as ``Chi0_pp``).

    Why the reindex is needed (round-2 finding, verified independently of
    this task's own ``dVpp_matrices``/``production_pair_bubble`` code):
    ``bare_bond_vertices``' enlarged slot ``I = m*nd + l1*norb+l2`` reuses
    ``S_bond``/``C_bond``'s SAME enlarged-index convention, and Task 5's
    ``ed_to_solver_bond_map`` (reviewed, load-bearing) established that
    THAT slot corresponds, when compared against an ED channel list, to the
    orbital-SWAPPED channel ``(delta_r[m], l2, l1)`` -- not
    ``(delta_r[m], l1, l2)`` directly. Confirmed HERE for ``Vpp``
    specifically by an independent numerical cross-check on the
    ALREADY-validated ph identity (``dW_matrices``/``pred_first_order``/
    ``ed_to_solver_bond_map``, untouched by this task): feeding
    ``canonical_density_terms`` the UNSWAPPED ``(a, b, R)`` argument
    reproduces the ph identity's dW-sandwiched prediction to the expected
    finite-Nmat/Richardson residual (~1.4e-4), while the SWAPPED argument
    fails by O(1) (~0.14, matching the signal) -- i.e. the ED-side
    convention is the SAME UNSWAPPED one used throughout this whole
    campaign (design doc: "Sigma_j n_{j,a} n_{j+R,b}"), and it is
    ``bare_bond_vertices``' Vpp slot -- reusing the ph enlarged-index
    convention for a DIFFERENT physical operator (a pair, not a bond
    bilinear) -- that needs reindexing to align with it, not the ED term.
    (An earlier draft of this task instead swapped the ED-side term to
    match the RAW, unreindexed ``dVpp`` -- numerically equivalent on the
    single direction it was tried on, but a round-2 review correctly
    flagged that as attributing the fix to the wrong object; seeing it
    fail to justify itself against the independently-validated ph
    identity is what surfaced this cleaner, correct account.)

    A NO-OP at norb=1 (``a=b`` always, so the involution's every channel is
    a fixed point): every fx5 unit direction this task adjudicates (U/g1/
    g2) is numerically IDENTICAL with or without this reindex -- it matters
    only for a genuinely inter-orbital direction (case M, Task 8), where it
    is REQUIRED.
    """
    channels = list(channels)
    bond_set = resolve_interactions(closed_decls, np.eye(3), norb=fx.norb)
    # Review round 3 (should_fix): dVpp's raw slot I = m*nd + l1*norb+l2 is
    # POSITIONAL, built from bond_set's own (m, l1, l2) nested-loop order --
    # a caller-supplied channels list in ANY other order would silently
    # misalign the reindex below with no error. Require channels to equal
    # the canonical positional order exactly (_ordered_pair_channels on
    # THIS bond_set), not merely a permutation of it.
    canonical_channels = _ordered_pair_channels(bond_set, fx.norb)
    if channels != canonical_channels:
        raise ValueError(
            "dVpp_matrices: channels must be exactly bond_set's own "
            "positional (R, l1, l2) order (_ordered_pair_channels(bond_set, "
            "norb)) -- got a different order/content, which would silently "
            "misalign the orbital-swap reindex below")

    index_of = {ch: i for i, ch in enumerate(channels)}
    # perm[I] = index of the orbital-swapped channel (R, b, a) for
    # channels[I] = (R, a, b) -- an involution (perm[perm[I]] == I), the
    # SAME (R,a,b)->(R,b,a) relabeling ed_to_solver_bond_map established
    # for the ph enlarged index, NOT exchange_projector's (R,a,b)->(-R,b,a).
    perm = np.array([index_of[(R, b, a)] for (R, a, b) in channels])

    zero_decls = {k: 0.0 * v for k, v in closed_decls.items()}

    def _vpp_at(decls):
        bs = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        kx = np.array([0.0])
        ky = np.array([0.0])
        kz = np.array([0.0])
        S0_q, C0_q = _build_bond_m0_blocks(bs, {}, {}, fx.norb, kx, ky, kz)
        _S, _C, Vpp_s, Vpp_t = bare_bond_vertices(bs, S0_q, C0_q, fx.norb)
        return Vpp_s, Vpp_t

    Vpp_s1, Vpp_t1 = _vpp_at(closed_decls)
    Vpp_s0, Vpp_t0 = _vpp_at(zero_decls)
    dVpp_s = (Vpp_s1 - Vpp_s0)[np.ix_(perm, perm)]
    dVpp_t = (Vpp_t1 - Vpp_t0)[np.ix_(perm, perm)]
    return dVpp_s, dVpp_t


# ---------------------------------------------------------------------------
# Task 7, Step 1: exchange_projector unit tests; pin 3a (production vs
# eigenbasis pair bubble, abs=1e-12); pin 3b (eigenbasis pair bubble vs
# SectorED.pair_correlator at V=0, eps3 diagnostic < 1e-5).
# ---------------------------------------------------------------------------

class TestExchangeProjector(ApproxTestCase):
    """``P @ P = I``, ``Q_s + Q_t = I``, ``Q^2 = Q`` (both projectors,
    idempotent and complementary), and the ``ValueError`` on a channel list
    not closed under the involution ``(R,a,b) -> (-R,b,a)`` (round-4 fix:
    the invariant was implicit)."""

    def _check(self, channels):
        n = len(channels)
        P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)
        Id = np.eye(n, dtype=complex)
        assert_approx_array(P @ P, Id, rel=0, abs=1e-14)
        assert_approx_array(Q_s + Q_t, Id, rel=0, abs=1e-14)
        assert_approx_array(Q_s @ Q_s, Q_s, rel=0, abs=1e-14)
        assert_approx_array(Q_t @ Q_t, Q_t, rel=0, abs=1e-14)
        assert_approx_array(Q_s @ Q_t, np.zeros((n, n)), rel=0, abs=1e-14)

    def test_fx5_pp_channels(self):
        bond_set = resolve_interactions(
            _two_sided_decl([(1, 0, 0, 1.0), (2, 0, 0, 1.0)]), np.eye(3), norb=1)
        self._check(_ordered_pair_channels(bond_set, 1))

    def test_fx3_pp_channels(self):
        bond_set = resolve_interactions(
            _two_sided_decl([(1, 0, 0, V1), (1, 0, 1, V1)]), np.eye(3), norb=2)
        self._check(_ordered_pair_channels(bond_set, 2))

    def test_fixed_point_channel_is_singlet_only(self):
        # An on-site same-orbital channel (R=0, a=b) is its own partner
        # under the involution -- Q_s must be nonzero and Q_t exactly zero
        # there (the triplet's fixed-point cancellation the brief and
        # projected_structural_zero_mask's docstring both name explicitly).
        channels = [(0, 0, 0), (1, 0, 0), (-1, 0, 0)]
        _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)
        self.assertApprox(Q_s[0, 0], 1.0, rel=0, abs=1e-14)
        self.assertApprox(Q_t[0, 0], 0.0, rel=0, abs=1e-14)

    def test_unclosed_list_raises(self):
        with self.assertRaises(ValueError):
            ed_oracle_util.exchange_projector([(1, 0, 0)])   # missing (-1,0,0)

    def test_duplicate_channel_raises(self):
        with self.assertRaises(ValueError):
            ed_oracle_util.exchange_projector([(0, 0, 0), (0, 0, 0)])


class TestPin3aProductionVsEigenbasisPairBubble(ApproxTestCase):
    """Pin 3a: ``production_pair_bubble`` vs ``eigenbasis_pair_bubble`` --
    two DIFFERENT computational routes (momentum-space ``_g2_from_green``
    plus a phase reconstruction, vs a direct real-space sum over the
    one-body eigenbasis) evaluating the SAME finite-Nmat window, so they
    must agree at round-off (``abs=1e-12``), both fixtures."""

    def _check(self, fx, decls):
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        channels = _ordered_pair_channels(bond_set, fx.norb)
        Xp = production_pair_bubble(fx, channels, NMAT)
        Xe = ed_oracle_util.eigenbasis_pair_bubble(fx, channels, NMAT)
        diff = float(np.max(np.abs(Xp - Xe)))
        print("Pin3a (fixture L={}, norb={}): diff={:.3e}".format(
            fx.L, fx.norb, diff))
        assert_approx_array(Xp, Xe, rel=0, abs=1e-12)

    def test_fx5(self):
        self._check(fx5(), _two_sided_decl([(1, 0, 0, V1), (2, 0, 0, 0.5 * V1)]))

    def test_fx3(self):
        self._check(fx3(), _two_sided_decl([(1, 0, 0, V1), (1, 0, 1, V1)]))


class TestPin3bEigenbasisPairBubbleVsEd(ApproxTestCase):
    """Pin 3b: ``eigenbasis_pair_bubble`` (Nmat-Richardson extrapolated,
    ``order=1`` matching its raw/no-tail O(1/Nmat) convergence -- same
    stencil as Pin 2b's ``eps2_rich``) vs the EXACT
    ``SectorED.pair_correlator`` at V=0, q=0. ``eps3`` is measured, printed
    and asserted ``< 1e-5`` -- a CALIBRATION DIAGNOSTIC ONLY (round-3 fix:
    NOT coupled into any pp granule's tolerance; each granule's own
    ``delta_nmat`` covers the derivative-level Nmat guard). A raw single-Nmat
    comparison would NOT clear this ceiling (measured ~2e-4 at NMAT=1024,
    ~4x the design's O(1/Nmat) law) -- the Richardson extrapolation is what
    makes this pin discriminating without an enormous Nmat, mirroring Pin
    2b's ``eps2_rich`` exactly.
    """

    def _check(self, fx, decls, n1):
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        channels = _ordered_pair_channels(bond_set, fx.norb)
        channels_ed = _to_sector_channels(fx, channels)
        exact = ed_oracle_util.SectorED(fx).pair_correlator(channels_ed)[0]
        rich = _richardson_nmat(
            lambda n: ed_oracle_util.eigenbasis_pair_bubble(fx, channels, n),
            n1, order=1)
        eps3 = float(np.max(np.abs(rich - exact)))
        print("Pin3b (fixture L={}, norb={}): eps3={:.3e}".format(
            fx.L, fx.norb, eps3))
        self.assertLess(eps3, 1e-5)

    def test_fx5(self):
        self._check(fx5(), _two_sided_decl([(1, 0, 0, V1), (2, 0, 0, 0.5 * V1)]),
                    n1=512)

    def test_fx3(self):
        self._check(fx3(), _two_sided_decl([(1, 0, 0, V1), (1, 0, 1, V1)]),
                    n1=512)


class TestPairBubbleNmatGuard(unittest.TestCase):
    """Review must_fix: an odd/zero/negative ``nmat`` is not a fermionic
    Matsubara grid (the centered-grid antisymmetry
    ``iwn[nmat-1-n] == -iwn[n]`` that ``eigenbasis_pair_bubble``'s
    ``G(-iwn)`` shortcut relies on only holds for even ``nmat``) -- both
    ``production_pair_bubble`` and ``eigenbasis_pair_bubble`` must reject
    it rather than silently building a broken grid (pin 3a could not have
    caught the two sides silently agreeing on the wrong quantity)."""

    def _channels(self):
        bond_set = resolve_interactions(
            _two_sided_decl([(1, 0, 0, 1.0)]), np.eye(3), norb=1)
        return _ordered_pair_channels(bond_set, 1)

    def test_odd_nmat_raises(self):
        channels = self._channels()
        with self.assertRaises(ValueError):
            production_pair_bubble(fx5(), channels, 33)
        with self.assertRaises(ValueError):
            ed_oracle_util.eigenbasis_pair_bubble(fx5(), channels, 33)

    def test_zero_nmat_raises(self):
        channels = self._channels()
        with self.assertRaises(ValueError):
            production_pair_bubble(fx5(), channels, 0)
        with self.assertRaises(ValueError):
            ed_oracle_util.eigenbasis_pair_bubble(fx5(), channels, 0)

    def test_negative_nmat_raises(self):
        channels = self._channels()
        with self.assertRaises(ValueError):
            production_pair_bubble(fx5(), channels, -4)
        with self.assertRaises(ValueError):
            ed_oracle_util.eigenbasis_pair_bubble(fx5(), channels, -4)

    def test_nonintegral_nmat_raises(self):
        # Review round 2: int(nmat) alone would truncate (not reject) a
        # non-integral value -- must be rejected explicitly.
        channels = self._channels()
        with self.assertRaises(ValueError):
            production_pair_bubble(fx5(), channels, 32.5)
        with self.assertRaises(ValueError):
            ed_oracle_util.eigenbasis_pair_bubble(fx5(), channels, 32.5)

    def test_nonfinite_nmat_raises(self):
        # Review round 2: bare int(nan)/int(inf) raise ValueError/
        # OverflowError, not the intended, consistent ValueError this
        # module's callers expect -- pinned directly.
        channels = self._channels()
        for bad in (float("nan"), float("inf"), None, "not-a-number"):
            with self.subTest(nmat=bad):
                with self.assertRaises(ValueError):
                    production_pair_bubble(fx5(), channels, bad)
                with self.assertRaises(ValueError):
                    ed_oracle_util.eigenbasis_pair_bubble(fx5(), channels, bad)


class TestPpInterOrbitalSanityFx3(ApproxTestCase):
    """Should-fix (review round 1): every pp granule adjudicated by this
    task is fx5 (norb=1, U/g1/g2 all have a=b=0), so orbital-index
    bookkeeping -- where a genuine bug class lives (design doc: "orbital
    transposition is this codebase's most recurrent defect class") -- has
    NO coverage anywhere else in this task. This is a LIGHTWEIGHT
    two-direction sanity check on case M (fx3), explicitly NOT a full
    ``adjudicate_granule`` campaign cell (that is Task 8's scope, which
    owns case M's formal pp granules) -- it exists only to catch what the
    fx5 granules structurally cannot see.

    FINDING, and how a round-2 review correctly caught the first draft's
    mistake (receiving-code-review discipline: verify, don't just
    implement, and don't just defend a first answer either): the first
    draft of this test made the pp identity agree by feeding
    ``canonical_density_terms`` a SWAPPED ``(b, a, R)`` orbital argument.
    Round 2 objected that this "can conceal an orbital-index error ...
    invalidating this test as an independent oracle" and asked for the
    UNSWAPPED ``(a, b, R)`` term instead. That objection was investigated
    on the ALREADY-validated ph identity first (``dW_matrices``/
    ``pred_first_order``/``ed_to_solver_bond_map``, untouched by this
    task): feeding ``canonical_density_terms`` the UNSWAPPED argument
    reproduces the ph prediction at the expected residual (~1.4e-4), the
    SWAPPED one fails by O(1) (~0.14) -- i.e. round 2 was RIGHT that the ED
    side should stay unswapped (matching the design doc's own definition,
    "Sigma_j n_{j,a} n_{j+R,b}"). The real defect was in
    ``dVpp_matrices``: ``bare_bond_vertices``' enlarged Vpp slot reuses
    ``S_bond``/``C_bond``'s SAME index convention, and Task 5's
    ``ed_to_solver_bond_map`` already established (reviewed, load-bearing)
    that THAT slot needs an orbital-swap reindex to align with a direct
    ``(R,a,b)`` channel list -- ``dVpp_matrices`` now performs that
    reindex explicitly (see its docstring) instead of the ED side being
    swapped to compensate. Both give numerically IDENTICAL results on the
    direction originally tried (confirmed directly, not assumed) --
    round 2's objection was about attribution/generality, not about a
    wrong number, and attribution is exactly what matters for a
    reindex that must also work on case M's other directions.

    Both g+ (``V_01(+1)=V_10(-1)``) and g- (``V_01(-1)=V_10(+1)``, Task 5's
    naming) are checked here, with the UNSWAPPED ED canonical term and the
    reindexed ``dVpp_matrices``: both converge to the expected finite-
    Nmat/Richardson residual scale (measured: diff_s=8.09e-5,
    diff_t=7.83e-5 for EITHER direction, despite g+/g- filling DISJOINT
    ``Vpp_s``/``Vpp_t`` slots on this fixture's B=3 topology -- e.g.
    ``{6,9}`` vs ``{5,10}``, confirmed by direct inspection). This gap
    (dVpp's reindex requirement) was invisible to every existing test:
    Task 5's ``TestPredFirstOrder.test_fx3`` (the only prior test
    exercising an inter-orbital declaration through the solver) never
    compares against ED, and every ED-side fx5 comparison across Tasks 3,
    5, 6 and this task's own granules uses a=b=0, where the reindex is a
    no-op. LOAD-BEARING for Task 8: its pp case-M granules must build
    ``dVpp`` via ``dVpp_matrices`` (which now carries the reindex), not a
    raw ``bare_bond_vertices`` call, for any inter-orbital direction.
    """

    def _check(self, R_decl):
        fx = fx3()
        # solver: V_01(R_decl) = V_10(-R_decl) (Task 5's g+/g- naming)
        decls = _two_sided_decl([(R_decl, 0, 1, 1.0)])
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        channels = _ordered_pair_channels(bond_set, fx.norb)
        _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)

        dVpp_s, dVpp_t = dVpp_matrices(fx, decls, channels)
        X0 = production_pair_bubble(fx, channels, NMAT)

        def _pred(dVpp, Q):
            D = -0.5 * (X0 @ dVpp @ X0)
            return Q @ D @ Q

        D_pred_s = _pred(dVpp_s, Q_s)
        D_pred_t = _pred(dVpp_t, Q_t)

        # NEGATIVE CONTROL (review round 3, should_fix -- make this an
        # executable assertion, not prose only): the RAW (unreindexed)
        # dVpp -- bare_bond_vertices' Vpp taken at its own positional
        # (m, l1, l2) slot with NO orbital-swap correction -- must FAIL
        # against the same ED target by a large margin, confirming the
        # reindex in dVpp_matrices is load-bearing here and not a no-op
        # dressed up as a fix.
        def _raw_dvpp():
            zero_decls = {k: 0.0 * v for k, v in decls.items()}

            def _vpp_at(dd):
                bs = resolve_interactions(dd, np.eye(3), norb=fx.norb)
                kx = np.array([0.0]); ky = np.array([0.0]); kz = np.array([0.0])
                s0, c0 = _build_bond_m0_blocks(bs, {}, {}, fx.norb, kx, ky, kz)
                _s, _c, vs, vt = bare_bond_vertices(bs, s0, c0, fx.norb)
                return vs, vt

            vs1, vt1 = _vpp_at(decls)
            vs0, vt0 = _vpp_at(zero_decls)
            return vs1 - vs0, vt1 - vt0

        raw_s, raw_t = _raw_dvpp()
        D_pred_s_raw = _pred(raw_s, Q_s)
        D_pred_t_raw = _pred(raw_t, Q_t)

        # ED canonical term: UNSWAPPED (a, b, R_decl) -- the design doc's
        # own convention, independently confirmed on the ph identity (see
        # class docstring); dVpp_matrices carries the reindex instead.
        terms_at = lambda v: ed_oracle_util.canonical_density_terms(
            fx, [(0, 1, R_decl, v)])

        @functools.lru_cache(maxsize=None)
        def _xpp_hf_sub(v):
            terms = terms_at(v)
            full = ed_oracle_util.SectorED(fx, terms=terms).pair_correlator(channels)
            hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
            hf_only = ed_oracle_util.SectorED(
                fx, terms=(), h1=hf_h1).pair_correlator(channels)
            return (full - hf_only)[0]

        D_ed = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1)
        D_ed_s = Q_s @ D_ed @ Q_s
        D_ed_t = Q_t @ D_ed @ Q_t

        diff_s = float(np.max(np.abs(D_pred_s - D_ed_s)))
        diff_t = float(np.max(np.abs(D_pred_t - D_ed_t)))
        diff_s_raw = float(np.max(np.abs(D_pred_s_raw - D_ed_s)))
        diff_t_raw = float(np.max(np.abs(D_pred_t_raw - D_ed_t)))
        max_signal_s = float(np.max(np.abs(D_pred_s)))
        max_signal_t = float(np.max(np.abs(D_pred_t)))
        print("TestPpInterOrbitalSanityFx3[R_decl={:+d}]: diff_s={:.3e} "
              "diff_t={:.3e} (raw/unreindexed: diff_s={:.3e} diff_t={:.3e})"
              .format(R_decl, diff_s, diff_t, diff_s_raw, diff_t_raw))
        # A loose sanity ceiling (measured residual ~8e-5; this is a
        # single-(V, Nmat) smoke check, not an adjudicated tolerance).
        assert_approx_array(D_pred_s, D_ed_s, rel=0, abs=1e-3)
        assert_approx_array(D_pred_t, D_ed_t, rel=0, abs=1e-3)
        # NEGATIVE CONTROL, asserted (review round 3): the raw/unreindexed
        # candidate must fail by O(1) -- at least half EACH channel's OWN
        # signal scale (review round 4: the singlet scale alone could
        # falsely pass/fail the triplet check if the two magnitudes
        # differ) -- confirming the reindex is load-bearing here, not a
        # documented-but-untested claim.
        self.assertGreater(diff_s_raw, 0.5 * max_signal_s)
        self.assertGreater(diff_t_raw, 0.5 * max_signal_t)

    def test_g_plus(self):
        self._check(1)

    def test_g_minus(self):
        self._check(-1)


# ---------------------------------------------------------------------------
# Task 7, Step 2: the U-anchor -- a shared prerequisite (review, must-fix
# 13). Runs the U-only pp control for BOTH channels; dependent V-case
# granules (Step 3) require it to hold, with NO skip path.
# ---------------------------------------------------------------------------

def _dvpp_support_masks(direction, channels):
    """The ANALYTIC (boolean, symbolic -- never thresholded-numeric) support
    of ``Q_eta (dVpp_eta) Q_eta`` for one fx5 pp unit direction, carrying the
    triplet's FIXED-POINT CANCELLATION exactly (``projected_structural_zero_
    mask``'s docstring: boolean ``I OR P`` propagation cannot represent
    this -- the caller must supply the analytic support directly).

    Verified against the builder (not guessed): ``bare_bond_vertices``' RAW
    (pre-projection) ``dVpp`` is DIAGONAL in the channel-orbital-pair basis
    for every fx5 unit direction (norb=1, so the local crossing of a scalar
    is trivially diagonal; the bond Cooper block ``D`` is diagonal by
    construction, and ``D + D^dagger`` stays diagonal). U (on-site, R=0)
    touches only the SELF-PAIRED channel ``(0,0,0)`` -- a FIXED POINT of the
    exchange involution, where ``Q_s`` is nonzero (``=1``) but ``Q_t`` is
    EXACTLY zero (the brief's "U direction: pp-t support all-False").
    g1/g2 (off-site, R=+-1/+-2) touch the NON-fixed-point pair
    ``{(+R,0,0),(-R,0,0)}``; ``Q_eta (D+D^dagger) Q_eta`` spreads this into
    the FULL 2x2 block spanning both channels for BOTH S and T (no
    cancellation -- confirmed numerically against the real
    ``bare_bond_vertices``/``dVpp_matrices`` output during this task's
    derivation).
    """
    channels = list(channels)
    n = len(channels)
    index_of = {ch: i for i, ch in enumerate(channels)}
    raw_channels = {
        "U": [(0, 0, 0)],
        "g1": [(1, 0, 0), (-1, 0, 0)],
        "g2": [(2, 0, 0), (-2, 0, 0)],
    }[direction]
    S_supp = np.zeros((n, n), dtype=bool)
    T_supp = np.zeros((n, n), dtype=bool)
    for ch in raw_channels:
        i = index_of[ch]
        R, a, b = ch
        j = index_of[(-R, b, a)]
        fixed = (i == j)
        for p in (i, j):
            for q in (i, j):
                S_supp[p, q] = True
                if not fixed:
                    T_supp[p, q] = True
    return S_supp, T_supp


class TestDvppMatricesChannelOrderGuard(unittest.TestCase):
    """Review round 3 (should_fix): ``dVpp_matrices``' orbital-swap reindex
    is positional -- it assumes ``channels[I]`` is bond_set's own
    ``(m, l1, l2)`` nested-loop order -- so a caller-supplied ``channels``
    in a DIFFERENT order would silently misalign the reindex with no
    error. Pinned directly: a reordered (but otherwise valid/complete)
    channel list must raise, not silently mis-permute."""

    def test_reordered_channels_raises(self):
        fx = fx3()
        decls = _two_sided_decl([(1, 0, 1, 1.0)])
        bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
        canonical = _ordered_pair_channels(bond_set, fx.norb)
        reordered = canonical[::-1]
        self.assertNotEqual(canonical, reordered)
        with self.assertRaises(ValueError):
            dVpp_matrices(fx, decls, reordered)


class TestDvppSupportMasksMatchNumericDvpp(unittest.TestCase):
    """Review round 2 should-fix: ``_dvpp_support_masks``' hand-derived
    ANALYTIC pattern, cross-checked against the boolean nonzero pattern of
    the ACTUAL numeric ``Q_eta @ dVpp_eta @ Q_eta`` (``dVpp_matrices`` +
    ``exchange_projector``, the real production/independent-projector
    output) -- mirrors ``TestDwSupportMasksMatchNumericDw``'s ph-channel
    pattern. Does not replace the analytic derivation (the whole point of
    ``_dvpp_support_masks`` is to be boolean-and-not-thresholded, never fed
    a numeric epsilon) but catches a mismatch between the hand-built
    pattern and what the real builder emits for each of the three unit
    directions this task adjudicates."""

    def _check(self, direction):
        fx = fx5()
        bond_set, _channels_ph, _smap = _bond_topology_and_map(fx)
        channels = _pp_channels(fx)
        _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)
        decls_at, _terms_at = _UNIT_DIRECTIONS[direction]
        dVpp_s, dVpp_t = dVpp_matrices(fx, decls_at(fx, 1.0), channels)
        numeric_S = np.abs(Q_s @ dVpp_s @ Q_s) > 1e-12
        numeric_T = np.abs(Q_t @ dVpp_t @ Q_t) > 1e-12
        S_supp, T_supp = _dvpp_support_masks(direction, channels)
        self.assertTrue(np.array_equal(numeric_S, S_supp),
                         "direction={} pp_s support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_S, S_supp))
        self.assertTrue(np.array_equal(numeric_T, T_supp),
                         "direction={} pp_t support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_T, T_supp))

    def test_U(self):
        self._check("U")

    def test_g1(self):
        self._check("g1")

    def test_g2(self):
        self._check("g2")


@functools.lru_cache(maxsize=None)
def _pp_channels(fx):
    """The shared fx5 B=5 topology's full ``(R,a,b)`` ordered-pair channel
    list, reusing ``_bond_topology_and_map``'s bond_set (same shared
    topology U/g1/g2 all declare, per Task 6's module note) so the pp
    granules align on the same grid as the ph ones."""
    bond_set, _channels_ph, _smap = _bond_topology_and_map(fx)
    return _ordered_pair_channels(bond_set, fx.norb)


@functools.lru_cache(maxsize=None)
def _pp_unit_direction_raw(direction):
    """Raw (un-adjudicated) ED- and solver-side first-order pp response
    arrays for one fx5 unit direction, pp_s and pp_t channels. Mirrors
    ``_unit_direction_raw``'s ph structure: solver side is exact-in-coupling
    (``dVpp`` at unit declared amplitude, evaluated at NMAT and 2*NMAT); ED
    side is ``richardson()`` of the HF-subtracted ``pair_correlator`` at
    ``q=0`` (the pp identity carries no q -- see the module note above),
    with the pp HF insertion (design doc: "the derivative of the free pair
    correlator after inserting the same normal-state HF one-body matrix on
    both legs"). Both sides are projected with the SAME (independently
    built) ``exchange_projector`` Q before comparison. lru_cache'd: Task 9
    aggregates without recomputation.
    """
    fx = fx5()
    decls_at, terms_at = _UNIT_DIRECTIONS[direction]
    channels = _pp_channels(fx)
    channels_ed = _to_sector_channels(fx, channels)
    _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)

    # -- solver side --------------------------------------------------
    dVpp_s, dVpp_t = dVpp_matrices(fx, decls_at(fx, 1.0), channels)
    X0_n = production_pair_bubble(fx, channels, NMAT)
    X0_2n = production_pair_bubble(fx, channels, 2 * NMAT)

    def _pred(X0, dVpp, Q):
        D = -0.5 * (X0 @ dVpp @ X0)
        return Q @ D @ Q

    D_pred_nmat_S = _pred(X0_n, dVpp_s, Q_s)
    D_pred_2nmat_S = _pred(X0_2n, dVpp_s, Q_s)
    D_pred_nmat_T = _pred(X0_n, dVpp_t, Q_t)
    D_pred_2nmat_T = _pred(X0_2n, dVpp_t, Q_t)

    # -- ED side ----------------------------------------------------------
    @functools.lru_cache(maxsize=None)
    def _xpp_hf_sub(v):
        terms = terms_at(fx, v)
        full = ed_oracle_util.SectorED(fx, terms=terms).pair_correlator(channels_ed)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).pair_correlator(channels_ed)
        return (full - hf_only)[0]   # q=0 slice -- the pp identity has no q

    D_ed_v1 = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1)
    D_ed_vhalf = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1 / 2)
    D_ed_v1_S = Q_s @ D_ed_v1 @ Q_s
    D_ed_vhalf_S = Q_s @ D_ed_vhalf @ Q_s
    D_ed_v1_T = Q_t @ D_ed_v1 @ Q_t
    D_ed_vhalf_T = Q_t @ D_ed_vhalf @ Q_t

    # -- structural zero mask ----------------------------------------------
    S_supp, T_supp = _dvpp_support_masks(direction, channels)
    free_support = np.ones_like(S_supp)   # dense free-fermion bubble (stated)
    zero_mask_S = ed_oracle_util.projected_structural_zero_mask(
        free_support, S_supp)
    zero_mask_T = ed_oracle_util.projected_structural_zero_mask(
        free_support, T_supp)

    return dict(
        D_ed_v1_S=D_ed_v1_S, D_ed_vhalf_S=D_ed_vhalf_S,
        D_pred_nmat_S=D_pred_nmat_S, D_pred_2nmat_S=D_pred_2nmat_S,
        D_ed_v1_T=D_ed_v1_T, D_ed_vhalf_T=D_ed_vhalf_T,
        D_pred_nmat_T=D_pred_nmat_T, D_pred_2nmat_T=D_pred_2nmat_T,
        zero_mask_S=zero_mask_S, zero_mask_T=zero_mask_T,
    )


@functools.lru_cache(maxsize=None)
def case_fx5_pp_unit_direction(direction):
    """pp adjudication (pp_s and pp_t granules) for one independent unit
    coupling direction on fx5. Returns ``dict(pp_s=<record>,
    pp_t=<record>)`` -- Task 9's aggregation interface, matching
    ``case_fx5_unit_direction``'s ph shape."""
    raw = _pp_unit_direction_raw(direction)
    rec_s = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_S"], raw["D_ed_vhalf_S"],
        raw["D_pred_nmat_S"], raw["D_pred_2nmat_S"],
        raw["zero_mask_S"], "fx5/{}/pp_s".format(direction))
    rec_t = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_T"], raw["D_ed_vhalf_T"],
        raw["D_pred_nmat_T"], raw["D_pred_2nmat_T"],
        raw["zero_mask_T"], "fx5/{}/pp_t".format(direction))
    return {"pp_s": rec_s, "pp_t": rec_t}


@functools.lru_cache(maxsize=1)
def _u_anchor():
    """Runs the U-only pp control for BOTH channels; returns
    ``{"pp_s": record, "pp_t": record}``. Dependent V-case tests require
    ``pp_s["status"] == "PASS"`` AND ``pp_t["status"] == "PASS-ZERO"`` (U is
    structurally null in the triplet channel -- a pp_t "PASS" would itself
    be a FINDING). An anchor violation fails LOUDLY, naming the anchor --
    NO skip path (raises ``AssertionError`` rather than returning a record a
    caller might not check).
    """
    rec = case_fx5_pp_unit_direction("U")
    if rec["pp_s"]["status"] != "PASS":
        raise AssertionError(
            "pp U-ANCHOR VIOLATION: pp_s status={} (required PASS) -- the "
            "U-only pp control anchors the -1/2 * X0_pp * dVpp * X0_pp "
            "identity's sign/normalization (design doc); every dependent "
            "V-case pp granule is unreliable until this is resolved. "
            "record={}".format(rec["pp_s"]["status"], rec["pp_s"]))
    if rec["pp_t"]["status"] != "PASS-ZERO":
        raise AssertionError(
            "pp U-ANCHOR VIOLATION: pp_t status={} (required PASS-ZERO -- U "
            "is structurally null in the triplet channel; a pp_t 'PASS' "
            "would itself be a FINDING, not a pass). record={}".format(
                rec["pp_t"]["status"], rec["pp_t"]))
    return rec


class TestUAnchor(unittest.TestCase):
    """Step 2: the U-only pp control, run for BOTH pp_s and pp_t. The U
    control anchors the acceptance identity's constants: production's
    on-site ``L_s = 2U`` under the ``-1/2`` convention must reproduce the
    physical ``-U`` pair kernel."""

    def test_pp_s_pass_pp_t_pass_zero(self):
        rec = _u_anchor()
        self.assertEqual(rec["pp_s"]["status"], "PASS")
        self.assertEqual(rec["pp_t"]["status"], "PASS-ZERO")


# ---------------------------------------------------------------------------
# Task 7, Step 3: pp cells for V one-shell (g1) and two-shell (g2) through
# adjudicate_granule. Dependent on the U-anchor holding (Step 2) -- calling
# _u_anchor() first fails loudly, naming the anchor, per its own "no skip
# path" contract.
# ---------------------------------------------------------------------------

class TestPpUnitDirectionsFx5(unittest.TestCase):
    """Step 3: the campaign's pp verdict cells (#151 Task 7) -- independent
    unit directions g1 (one shell, Delta r = +-1) and g2 (two shell,
    Delta r = +-2) on fx5, pp channels s (singlet) and t (triplet).
    THE VERDICT MATTERS MORE THAN GREEN: a FAIL is kept asserting (not
    weakened), recording a genuine production finding on #82's Cooper
    vertex; see this task's report for the measured audit tuples."""

    def _check(self, direction):
        _u_anchor()   # NO skip path: raises loudly if the anchor is violated
        rec = case_fx5_pp_unit_direction(direction)
        for channel in ("pp_s", "pp_t"):
            with self.subTest(direction=direction, channel=channel):
                r = rec[channel]
                self.assertIn(
                    r["status"], ("PASS", "PASS-ZERO"),
                    "granule {} status={} (delta_rich={:.3e} delta_nmat={:.3e} "
                    "tol={:.3e} max_signal={:.3e} first_failures={})".format(
                        r["label"], r["status"], r["delta_rich"],
                        r["delta_nmat"], r["tol"], r["max_signal"],
                        r["failures"][:5]))

    def test_g1(self):
        self._check("g1")

    def test_g2(self):
        self._check("g2")


# ---------------------------------------------------------------------------
# Task 8: Case M -- multi-orbital ph + pp adjudication on fx3 (#151)
# ---------------------------------------------------------------------------
#
# Case M exercises the orbital-pair bookkeeping invisible at norb=1: the ph
# and pp channels both PASSED on fx5 (single orbital, Tasks 6-7), but every
# genuine orbital-index defect in this codebase (design doc: "orbital
# transposition is this codebase's most recurrent defect class") is a no-op
# there. This section repeats BOTH channel families on fx3 (case M, norb=2),
# directions g+ (Delta r=+1, (a,b)=(0,1)) and g- (Delta r=-1, (a,b)=(0,1)) --
# Task 5's naming, matching TestPpInterOrbitalSanityFx3.
#
# Term lists (brief, review round 2 -- exact dict literals):
#   g+: ED canonical_density_terms(fx3, [(0, 1, +1, g)]);
#       solver closed {((+1,0,0),(0,1)): g, ((-1,0,0),(1,0)): g}
#         == _two_sided_decl([(1, 0, 1, g)]) for real g (pinned below).
#   g-: ED canonical_density_terms(fx3, [(0, 1, -1, g)]);
#       solver closed {((-1,0,0),(0,1)): g, ((+1,0,0),(1,0)): g}
#         == _two_sided_decl([(-1, 0, 1, g)]) for real g (pinned below).
# NEVER feed the closed two-sided dict to the ED side (Task 3's factor-of-2
# double count); the pp side MUST build dVpp via Task 7's dVpp_matrices --
# the orbital-swap reindex is LOAD-BEARING here (Task 7's finding), not a
# no-op the way it is on every fx5 unit direction.

def _decls_case_m_gplus(fx, v):
    return _two_sided_decl([(1, 0, 1, v)])


def _terms_case_m_gplus(fx, v):
    return ed_oracle_util.canonical_density_terms(fx, [(0, 1, 1, v)])


def _decls_case_m_gminus(fx, v):
    return _two_sided_decl([(-1, 0, 1, v)])


def _terms_case_m_gminus(fx, v):
    return ed_oracle_util.canonical_density_terms(fx, [(0, 1, -1, v)])


_CASE_M_DIRECTIONS = {
    "g+": (_decls_case_m_gplus, _terms_case_m_gplus),
    "g-": (_decls_case_m_gminus, _terms_case_m_gminus),
}


class TestCaseMDeclLiterals(unittest.TestCase):
    """Pins the brief's exact closed-declaration dict literals against
    ``_two_sided_decl``'s own output, so a typo in the convention above
    cannot silently drift from the reviewed spec."""

    def test_g_plus_matches_brief_literal(self):
        g = 0.37
        got = _decls_case_m_gplus(fx3(), g)
        want = {((1, 0, 0), (0, 1)): g, ((-1, 0, 0), (1, 0)): g}
        self.assertEqual(got, want)

    def test_g_minus_matches_brief_literal(self):
        g = 0.37
        got = _decls_case_m_gminus(fx3(), g)
        want = {((-1, 0, 0), (0, 1)): g, ((1, 0, 0), (1, 0)): g}
        self.assertEqual(got, want)


@functools.lru_cache(maxsize=1)
def _case_m_topology_and_map():
    """(bond_set, channels, smap) for case M's single-shell (Delta r=+-1)
    inter-orbital topology, shared by g+ and g- (both channels m=+-1 exist
    in the SAME bond_set regardless of which one carries the coupling --
    only the declared VALUE differs between directions, never the shell
    set; mirrors ``_bond_topology_and_map``'s fx5 pattern, Task 6). Channels
    are always the full ``(R, l1, l2)`` triple (never the norb=1 ``(R,)``
    shorthand) -- the orbital-resolved rows the brief requires."""
    fx = fx3()
    decls = _two_sided_decl([(1, 0, 1, 1.0)])
    bond_set = resolve_interactions(decls, np.eye(3), norb=fx.norb)
    norb = fx.norb
    channels = []
    for m in range(bond_set.n_channels):
        R = bond_set.delta_r[m][0]
        for l1 in range(norb):
            for l2 in range(norb):
                channels.append((R, l1, l2))
    smap = ed_to_solver_bond_map(channels, bond_set)
    return bond_set, channels, smap


@functools.lru_cache(maxsize=None)
def _case_m_chibar(nmat):
    """The bare bond bubble at ``nmat`` on case M's shared topology (shared
    by g+ and g-, mirrors ``_shared_chibar``)."""
    bond_set, _channels, _smap = _case_m_topology_and_map()
    fx = fx3()
    green = free_green(fx, nmat)
    return bond_bubble(green, bond_set, beta=fx.beta)[:, 0, 0]


def _dw_support_masks_case_m(direction, bond_set):
    """ANALYTIC support of dS/dC for one case-M direction, derived directly
    from ``bare_bond_vertices``' documented element equations (``dW_matrices``'
    module note: "Case-M's g+/g- ... carry the SAME m!=0 diagonal pattern at
    the (a,b)/(b,a) orbital-pair slot, plus a q-dependent CROSS Hartree entry
    at the (aa),(bb) slot with a!=b").

    g+ declares V_01(+1)=g, V_10(-1)=g (no R=0 component): the m!=0 Fock
    diagonal touches I=(m_{+1}, idx=(0,1)) and I=(m_{-1}, idx=(1,0)) in BOTH
    S and C (``bare_bond_vertices``: ``S_bond[I,I]=+V``, ``C_bond[I,I]=-V``);
    since ``v_onsite`` is untouched (zero) here, the m=0 Fock sub-block is
    exactly zero in S -- S carries NO Hartree term, matching the fx5
    pattern. The R!=0 Hartree DOES reach the m=0 charge block, but at the
    CROSS (aa,bb)=(0,0),(1,1) slot (idx 0, 3) -- not the diagonal (0,0)/
    (3,3) slot fx5's a=b case used -- since ``_build_bond_m0_blocks``'s
    ``C0[...,a*norb+a,b*norb+b] += 2*V_ab(q)`` is keyed by (a,b)=(0,1) here,
    off-diagonal in the ``idx=a*norb+a`` labeling. g- carries the SAME
    pattern with the m_{+1}/m_{-1} roles of the two orbital pairs swapped.
    """
    norb = bond_set.v_bond[0].shape[0]
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    S_supp = np.zeros((ND, ND), dtype=bool)
    C_supp = np.zeros((ND, ND), dtype=bool)
    m_plus = bond_set.delta_r.index((1, 0, 0))
    m_minus = bond_set.delta_r.index((-1, 0, 0))
    if direction == "g+":
        touched = [(m_plus, 0, 1), (m_minus, 1, 0)]
    else:
        touched = [(m_minus, 0, 1), (m_plus, 1, 0)]
    for (m, l1, l2) in touched:
        I = m * nd + l1 * norb + l2
        S_supp[I, I] = True
        C_supp[I, I] = True
    C_supp[0, 3] = True   # (aa,bb)=(0,0),(1,1) cross Hartree, C only
    C_supp[3, 0] = True
    return S_supp, C_supp


class TestDwSupportMasksMatchNumericDwCaseM(unittest.TestCase):
    """Cross-checks ``_dw_support_masks_case_m`` against the actual numeric
    ``dW_matrices`` output (mirrors ``TestDwSupportMasksMatchNumericDw``,
    Task 6's should-fix precedent) -- catches a mismatch between the
    hand-derived case-M pattern above and what ``bare_bond_vertices``
    actually emits."""

    def _check(self, direction):
        bond_set, _channels, _smap = _case_m_topology_and_map()
        decls_at, _terms_at = _CASE_M_DIRECTIONS[direction]
        dS, dC = dW_matrices(fx3(), decls_at(fx3(), 1.0))
        numeric_S = np.any(np.abs(dS) > 0, axis=0)
        numeric_C = np.any(np.abs(dC) > 0, axis=0)
        S_supp, C_supp = _dw_support_masks_case_m(direction, bond_set)
        self.assertTrue(np.array_equal(numeric_S, S_supp),
                         "direction={} S support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_S, S_supp))
        self.assertTrue(np.array_equal(numeric_C, C_supp),
                         "direction={} C support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_C, C_supp))

    def test_g_plus(self):
        self._check("g+")

    def test_g_minus(self):
        self._check("g-")


@functools.lru_cache(maxsize=None)
def _case_m_unit_direction_raw(direction):
    """Raw ph ED/solver first-order response arrays for one case-M
    direction (g+ or g-), S and C channels. Mirrors ``_unit_direction_raw``
    exactly (Task 6), on fx3 instead of fx5."""
    fx = fx3()
    decls_at, terms_at = _CASE_M_DIRECTIONS[direction]
    bond_set, channels, smap = _case_m_topology_and_map()

    # -- solver side --------------------------------------------------
    dS, dC = dW_matrices(fx, decls_at(fx, 1.0))
    chibar_n = _case_m_chibar(NMAT)
    chibar_2n = _case_m_chibar(2 * NMAT)
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
    S_supp, C_supp = _dw_support_masks_case_m(direction, bond_set)
    free_support = np.ones_like(S_supp)   # dense fx3 bubble (t_01 != 0)
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
def case_fx3_unit_direction(direction):
    """ph adjudication (S and C granules) for one case-M direction (g+/g-)
    on fx3. Returns ``dict(S=<record>, C=<record>)``, matching
    ``case_fx5_unit_direction``'s shape -- Task 9's aggregation interface.
    ``functools.lru_cache``'d: no recomputation, no test-order dependence.
    """
    raw = _case_m_unit_direction_raw(direction)
    rec_S = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_S"], raw["D_ed_vhalf_S"],
        raw["D_pred_nmat_S"], raw["D_pred_2nmat_S"],
        raw["zero_mask_S"], "fx3/{}/S".format(direction))
    rec_C = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_C"], raw["D_ed_vhalf_C"],
        raw["D_pred_nmat_C"], raw["D_pred_2nmat_C"],
        raw["zero_mask_C"], "fx3/{}/C".format(direction))
    return {"S": rec_S, "C": rec_C}


@functools.lru_cache(maxsize=None)
def _case_m_pp_channels():
    """The case-M shared topology's full ``(R,a,b)`` ordered-pair pp
    channel list (mirrors ``_pp_channels``, Task 7), reusing
    ``_case_m_topology_and_map``'s bond_set so the pp granules align on the
    same grid as the ph ones."""
    bond_set, _channels_ph, _smap = _case_m_topology_and_map()
    return _ordered_pair_channels(bond_set, fx3().norb)


def _dvpp_support_masks_case_m(direction, channels):
    """ANALYTIC support of ``Q_eta (dVpp_eta) Q_eta`` for one case-M pp unit
    direction (mirrors ``_dvpp_support_masks``, Task 7). g+/g- both touch
    the NON-fixed-point pair ``{(-1,0,1),(1,1,0)}`` / ``{(1,0,1),(-1,1,0)}``
    -- an inter-orbital bond channel is NEVER its own exchange partner
    (a!=b always here, unlike fx5's on-site same-orbital U), so BOTH pp_s
    and pp_t get the SAME full 2x2 block (no triplet fixed-point
    cancellation).

    NOTE the channel pairing here is the MIRROR IMAGE of the raw
    (pre-reindex) ``bare_bond_vertices`` diagonal touched by each direction
    (g+'s raw D touches ``(1,0,1)``/``(-1,1,0)``; g-'s touches
    ``(-1,0,1)``/``(1,1,0)``) -- ``dVpp_matrices``' own orbital-swap
    reindex permutes ``(R,a,b) -> (R,b,a)`` on BOTH matrix axes before this
    module ever sees it (Task 7's finding, load-bearing here), which swaps
    which of the two channels in each pp pair ends up carrying the nonzero
    diagonal entry in the REINDEXED basis this function's caller compares
    against. Verified against the actual numeric ``dVpp_matrices`` output
    below (``TestDvppSupportMasksMatchNumericDvppCaseM``), not guessed --
    an initial draft assigned the raw (unswapped) pairing here and failed
    that cross-check by construction, which is what surfaced this note."""
    channels = list(channels)
    n = len(channels)
    index_of = {ch: i for i, ch in enumerate(channels)}
    raw_channels = {
        "g+": [(-1, 0, 1), (1, 1, 0)],
        "g-": [(1, 0, 1), (-1, 1, 0)],
    }[direction]
    S_supp = np.zeros((n, n), dtype=bool)
    T_supp = np.zeros((n, n), dtype=bool)
    for ch in raw_channels:
        i = index_of[ch]
        R, a, b = ch
        j = index_of[(-R, b, a)]
        fixed = (i == j)
        for p in (i, j):
            for q in (i, j):
                S_supp[p, q] = True
                if not fixed:
                    T_supp[p, q] = True
    return S_supp, T_supp


class TestDvppSupportMasksMatchNumericDvppCaseM(unittest.TestCase):
    """Cross-checks ``_dvpp_support_masks_case_m`` against the actual
    numeric ``Q_eta @ dVpp_eta @ Q_eta`` output (mirrors
    ``TestDvppSupportMasksMatchNumericDvpp``, Task 7's should-fix
    precedent)."""

    def _check(self, direction):
        fx = fx3()
        channels = _case_m_pp_channels()
        _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)
        decls_at, _terms_at = _CASE_M_DIRECTIONS[direction]
        dVpp_s, dVpp_t = dVpp_matrices(fx, decls_at(fx, 1.0), channels)
        numeric_S = np.abs(Q_s @ dVpp_s @ Q_s) > 1e-12
        numeric_T = np.abs(Q_t @ dVpp_t @ Q_t) > 1e-12
        S_supp, T_supp = _dvpp_support_masks_case_m(direction, channels)
        self.assertTrue(np.array_equal(numeric_S, S_supp),
                         "direction={} pp_s support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_S, S_supp))
        self.assertTrue(np.array_equal(numeric_T, T_supp),
                         "direction={} pp_t support mismatch: numeric={} "
                         "analytic={}".format(direction, numeric_T, T_supp))

    def test_g_plus(self):
        self._check("g+")

    def test_g_minus(self):
        self._check("g-")


@functools.lru_cache(maxsize=None)
def _case_m_pp_unit_direction_raw(direction):
    """Raw pp ED/solver first-order response arrays for one case-M
    direction, pp_s/pp_t. Mirrors ``_pp_unit_direction_raw`` (Task 7), on
    fx3. ``dVpp`` is built via Task 7's ``dVpp_matrices`` -- the
    orbital-swap reindex is LOAD-BEARING here (Task 7's finding), not a raw
    ``bare_bond_vertices`` call."""
    fx = fx3()
    decls_at, terms_at = _CASE_M_DIRECTIONS[direction]
    channels = _case_m_pp_channels()
    channels_ed = _to_sector_channels(fx, channels)
    _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)

    # -- solver side --------------------------------------------------
    dVpp_s, dVpp_t = dVpp_matrices(fx, decls_at(fx, 1.0), channels)
    X0_n = production_pair_bubble(fx, channels, NMAT)
    X0_2n = production_pair_bubble(fx, channels, 2 * NMAT)

    def _pred(X0, dVpp, Q):
        D = -0.5 * (X0 @ dVpp @ X0)
        return Q @ D @ Q

    D_pred_nmat_S = _pred(X0_n, dVpp_s, Q_s)
    D_pred_2nmat_S = _pred(X0_2n, dVpp_s, Q_s)
    D_pred_nmat_T = _pred(X0_n, dVpp_t, Q_t)
    D_pred_2nmat_T = _pred(X0_2n, dVpp_t, Q_t)

    # -- ED side ----------------------------------------------------------
    @functools.lru_cache(maxsize=None)
    def _xpp_hf_sub(v):
        terms = terms_at(fx, v)
        full = ed_oracle_util.SectorED(fx, terms=terms).pair_correlator(channels_ed)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).pair_correlator(channels_ed)
        return (full - hf_only)[0]   # q=0 slice -- the pp identity has no q

    D_ed_v1 = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1)
    D_ed_vhalf = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1 / 2)
    D_ed_v1_S = Q_s @ D_ed_v1 @ Q_s
    D_ed_vhalf_S = Q_s @ D_ed_vhalf @ Q_s
    D_ed_v1_T = Q_t @ D_ed_v1 @ Q_t
    D_ed_vhalf_T = Q_t @ D_ed_vhalf @ Q_t

    # -- structural zero mask ----------------------------------------------
    S_supp, T_supp = _dvpp_support_masks_case_m(direction, channels)
    free_support = np.ones_like(S_supp)   # dense fx3 bubble (stated)
    zero_mask_S = ed_oracle_util.projected_structural_zero_mask(
        free_support, S_supp)
    zero_mask_T = ed_oracle_util.projected_structural_zero_mask(
        free_support, T_supp)

    return dict(
        D_ed_v1_S=D_ed_v1_S, D_ed_vhalf_S=D_ed_vhalf_S,
        D_pred_nmat_S=D_pred_nmat_S, D_pred_2nmat_S=D_pred_2nmat_S,
        D_ed_v1_T=D_ed_v1_T, D_ed_vhalf_T=D_ed_vhalf_T,
        D_pred_nmat_T=D_pred_nmat_T, D_pred_2nmat_T=D_pred_2nmat_T,
        zero_mask_S=zero_mask_S, zero_mask_T=zero_mask_T,
    )


@functools.lru_cache(maxsize=None)
def case_fx3_pp_unit_direction(direction):
    """pp adjudication (pp_s and pp_t granules) for one case-M direction on
    fx3. Returns ``dict(pp_s=<record>, pp_t=<record>)``, matching
    ``case_fx5_pp_unit_direction``'s shape. Case M's g+/g- have no U-anchor
    precondition (unlike fx5's g1/g2, which required Task 7's ``_u_anchor``)
    -- the acceptance identity's constants were already anchored by Task 7's
    fx5 U control; this direction pair only needs the ORBITAL bookkeeping
    (the reindex) to hold, which is exactly what ``dVpp_matrices`` supplies.
    """
    raw = _case_m_pp_unit_direction_raw(direction)
    rec_s = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_S"], raw["D_ed_vhalf_S"],
        raw["D_pred_nmat_S"], raw["D_pred_2nmat_S"],
        raw["zero_mask_S"], "fx3/{}/pp_s".format(direction))
    rec_t = ed_oracle_util.adjudicate_granule(
        raw["D_ed_v1_T"], raw["D_ed_vhalf_T"],
        raw["D_pred_nmat_T"], raw["D_pred_2nmat_T"],
        raw["zero_mask_T"], "fx3/{}/pp_t".format(direction))
    return {"pp_s": rec_s, "pp_t": rec_t}


# --- joint ray (1,1) = g+ + g-, both ph and pp (brief's superposition check) ---

def _terms_ray_case_m(fx, t):
    return ed_oracle_util.canonical_density_terms(
        fx, [(0, 1, 1, t), (0, 1, -1, t)])


def _joint_ray_superposition_check_case_m(testcase, directions_alphas,
                                           terms_at_ray, ray_label):
    """ph joint-ray superposition check on case M (mirrors
    ``_joint_ray_superposition_check``, Task 6, generalized off its fx5/
    ``_bond_topology_and_map``/``case_fx5_unit_direction`` hardcoding).
    Pure ED-side linearity check: the joint ray's OWN Richardson derivative
    must equal the weighted sum of the already-measured unit directions'
    own ``D_ed`` -- NOT a granule, NOT a sensitivity input (Task 6's
    distinction, carried forward unchanged)."""
    fx = fx3()
    _bond_set, channels, smap = _case_m_topology_and_map()

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

        unit_tols = [case_fx3_unit_direction(d)[chan_name]["tol"]
                     for d, _alpha in directions_alphas]
        tol_joint = 10.0 * max(delta_rich_joint, max(unit_tols))

        D_sum = sum(
            alpha * _case_m_unit_direction_raw(d)["D_ed_vhalf_" + chan_name]
            for d, alpha in directions_alphas)

        diff = float(np.max(np.abs(D_joint_vhalf - D_sum)))
        print("joint[{}/{}]: delta_rich_joint={:.6e} tol_joint={:.6e} "
              "diff={:.6e}".format(ray_label, chan_name, delta_rich_joint,
                                    tol_joint, diff))
        with testcase.subTest(ray=ray_label, channel=chan_name):
            testcase.assertLessEqual(diff, tol_joint)


def _pp_joint_ray_superposition_check_case_m(testcase, directions_alphas,
                                              terms_at_ray, ray_label):
    """pp joint-ray superposition check, case M -- the FIRST pp joint-ray
    superposition check in the campaign (Task 6 built only the ph version
    on fx5; this generalizes the same idea to the pp channel, projected
    with the SAME (independently built) ``exchange_projector`` Q used by
    every other pp comparison in this module)."""
    fx = fx3()
    channels = _case_m_pp_channels()
    channels_ed = _to_sector_channels(fx, channels)
    _P, Q_s, Q_t = ed_oracle_util.exchange_projector(channels)

    @functools.lru_cache(maxsize=None)
    def _xpp_hf_sub(t):
        terms = terms_at_ray(fx, t)
        full = ed_oracle_util.SectorED(fx, terms=terms).pair_correlator(channels_ed)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).pair_correlator(channels_ed)
        return (full - hf_only)[0]

    D_joint_v1 = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1)
    D_joint_vhalf = ed_oracle_util.richardson(_xpp_hf_sub, CAMPAIGN_V1 / 2)

    for chan_name, Q, raw_key in (("pp_s", Q_s, "D_ed_vhalf_S"),
                                   ("pp_t", Q_t, "D_ed_vhalf_T")):
        Dj_v1 = Q @ D_joint_v1 @ Q
        Dj_vhalf = Q @ D_joint_vhalf @ Q
        delta_rich_joint = float(np.max(np.abs(Dj_v1 - Dj_vhalf)))

        unit_tols = [case_fx3_pp_unit_direction(d)[chan_name]["tol"]
                     for d, _alpha in directions_alphas]
        tol_joint = 10.0 * max(delta_rich_joint, max(unit_tols))

        D_sum = sum(
            alpha * _case_m_pp_unit_direction_raw(d)[raw_key]
            for d, alpha in directions_alphas)

        diff = float(np.max(np.abs(Dj_vhalf - D_sum)))
        print("joint[{}/{}]: delta_rich_joint={:.6e} tol_joint={:.6e} "
              "diff={:.6e}".format(ray_label, chan_name, delta_rich_joint,
                                    tol_joint, diff))
        with testcase.subTest(ray=ray_label, channel=chan_name):
            testcase.assertLessEqual(diff, tol_joint)


class TestCaseM(unittest.TestCase):
    """Step 1 (#151 Task 8): the campaign's case-M multi-orbital verdict
    cells -- directions g+ (Delta r=+1, (a,b)=(0,1)) / g- (Delta r=-1,
    (a,b)=(0,1)) on fx3, ph channels S/C and pp channels s/t (8 granules:
    (fx3,g+,S), (fx3,g+,C), (fx3,g-,S), (fx3,g-,C), (fx3,g+,pp-s),
    (fx3,g+,pp-t), (fx3,g-,pp-s), (fx3,g-,pp-t)), plus the joint ray
    (1,1) = g+ + g- superposition check on BOTH channel families (not
    granules -- pure ED-side linearity checks, per Task 6's distinction).
    THE VERDICT MATTERS MORE THAN GREEN (brief): a FAIL or INCONCLUSIVE is
    kept asserting (a production finding on #82's orbital-pair bookkeeping,
    the codebase's most recurrent defect class), never adjusted, never
    loosened; see this task's report for the measured audit tuples.

    Reproducible timing (brief Step 2, informational, NOT CI-enforced):
        /usr/bin/time -p python -m unittest tests.test_bond_vs_ed_oracle.TestCaseM
    """

    def _check_ph(self, direction):
        rec = case_fx3_unit_direction(direction)
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

    def _check_pp(self, direction):
        rec = case_fx3_pp_unit_direction(direction)
        for channel in ("pp_s", "pp_t"):
            with self.subTest(direction=direction, channel=channel):
                r = rec[channel]
                self.assertIn(
                    r["status"], ("PASS", "PASS-ZERO"),
                    "granule {} status={} (delta_rich={:.3e} delta_nmat={:.3e} "
                    "tol={:.3e} max_signal={:.3e} first_failures={})".format(
                        r["label"], r["status"], r["delta_rich"],
                        r["delta_nmat"], r["tol"], r["max_signal"],
                        r["failures"][:5]))

    def test_ph_g_plus(self):
        self._check_ph("g+")

    def test_ph_g_minus(self):
        self._check_ph("g-")

    def test_pp_g_plus(self):
        self._check_pp("g+")

    def test_pp_g_minus(self):
        self._check_pp("g-")

    def test_joint_ray_ph_1_1(self):
        _joint_ray_superposition_check_case_m(
            self, [("g+", 1.0), ("g-", 1.0)], _terms_ray_case_m,
            "caseM-ph-g+g-(1,1)")

    def test_joint_ray_pp_1_1(self):
        _pp_joint_ray_superposition_check_case_m(
            self, [("g+", 1.0), ("g-", 1.0)], _terms_ray_case_m,
            "caseM-pp-g+g-(1,1)")


def _fake_record(pred_full, bearing_mask, status="PASS", delta_rich=1e-6,
                  label="synthetic"):
    """A synthetic ``adjudicate_granule``-shaped record for
    ``TestSensitivityRank``'s unit tests below -- fast, no ED/solver cost,
    mirrors exactly the fields ``sensitivity_rank`` reads
    (``status``, ``delta_rich``, ``pred_full``, ``bearing_mask``)."""
    pred_full = np.asarray(pred_full, dtype=float)
    bearing_mask = np.asarray(bearing_mask, dtype=bool)
    return dict(label=label, delta_rich=delta_rich, delta_nmat=delta_rich,
                tol=max(10.0 * delta_rich, 1e-13),
                max_signal=float(np.max(np.abs(pred_full[bearing_mask])))
                if bearing_mask.any() else 0.0,
                status=status, failures=[], pred_full=pred_full,
                bearing_mask=bearing_mask)


class TestSensitivityExpectedTablesInvariant(unittest.TestCase):
    """Optional finding (review, round 1): nothing previously asserted
    that ``SENSITIVITY_EXPECTED_ACTIVE`` and ``SENSITIVITY_EXPECTED_NULL``
    share the same label set with DISJOINT active/null direction sets per
    label -- both hold today, but a future table edit could silently
    break the partition reasoning ``sensitivity_rank`` relies on."""

    def test_same_label_set(self):
        self.assertEqual(
            set(ed_oracle_util.SENSITIVITY_EXPECTED_ACTIVE),
            set(ed_oracle_util.SENSITIVITY_EXPECTED_NULL))

    def test_active_and_null_disjoint_per_label(self):
        for label, active in ed_oracle_util.SENSITIVITY_EXPECTED_ACTIVE.items():
            null = ed_oracle_util.SENSITIVITY_EXPECTED_NULL[label]
            with self.subTest(label=label):
                self.assertEqual(set(active) & set(null), set())

    def test_no_duplicate_directions_within_a_table_entry(self):
        # optional finding (review, round 2): a duplicated name within one
        # label's active/null tuple would be silently collapsed by every
        # SET-based completeness check above yet still appear TWICE in
        # ``active_names`` (a list, built by iterating the tuple in
        # order), producing a duplicate SVD column and a misleading rank
        # failure that the set-equality checks alone cannot catch.
        for label, active in ed_oracle_util.SENSITIVITY_EXPECTED_ACTIVE.items():
            null = ed_oracle_util.SENSITIVITY_EXPECTED_NULL[label]
            with self.subTest(label=label, table="active"):
                self.assertEqual(len(active), len(set(active)))
            with self.subTest(label=label, table="null"):
                self.assertEqual(len(null), len(set(null)))


class TestSensitivityRank(unittest.TestCase):
    """Step 1 unit tests for ``ed_oracle_util.sensitivity_rank`` on
    synthetic records (fast, no ED/solver cost) -- added during this
    task's own review round (codex_review_diff via the ai-review-cycle
    MCP) to close two must_fix gaps and exercise two should_fix guards
    the campaign's own (all-clean) data never happens to trigger:
    (1) an INCOMPLETE ``granule_records`` input (missing a null-anchor
    direction) must raise, not silently pass with only the present
    directions checked; (2) an UNEXPECTED extra direction must also
    raise; (3) an underdetermined matrix (fewer bearing rows than active
    columns) must be marked INCONCLUSIVE without ever calling ``svd`` on
    a degenerate shape; (4) genuinely parallel/rank-deficient columns (a
    LOW ``sv_ratio``, rows >= columns) must be marked INCONCLUSIVE via the
    ``SENS_SV_FLOOR`` gate specifically -- the campaign's real fx3 g+/g-
    columns turned out NOT to hit this path (Task 9 report, "rank
    analysis"), so this is the only place that gate's actual behavior is
    exercised at all.
    """

    def test_missing_null_anchor_direction_raises(self):
        # fx5/pp_t's a-priori set is active=(g1, g2), null=(U,) -- an
        # input missing U entirely (only g1/g2 present) must raise, not
        # silently read as "active == {g1, g2}, matches expectation".
        recs = {
            "g1": _fake_record([1.0, 0.0, 0.0, 1.0], [True, False, False, True]),
            "g2": _fake_record([0.0, 1.0, 1.0, 0.0], [False, True, True, False]),
        }
        with self.assertRaises(AssertionError):
            ed_oracle_util.sensitivity_rank(recs, "fx5/pp_t")

    def test_unexpected_extra_direction_raises(self):
        recs = {
            "U": _fake_record([1.0, 0.0], [True, False]),
            "g1": _fake_record([0.0, 1.0], [False, True]),
            "g2": _fake_record([1.0, 1.0], [True, True]),
            "g3": _fake_record([2.0, 2.0], [True, True]),   # not expected
        }
        with self.assertRaises(AssertionError):
            ed_oracle_util.sensitivity_rank(recs, "fx5/S")

    def test_active_direction_measured_null_raises(self):
        # fx5/pp_t expects active=(g1, g2), null=(U,). Supply the FULL
        # expected key set (U present and correctly PASS-ZERO, satisfying
        # the completeness check) but make g2 come back PASS-ZERO too --
        # a column active in expectation but structurally null in fact,
        # which must be caught by the active/null PARTITION assertion
        # specifically (review should_fix, round 2: the original version
        # of this test omitted U entirely and only ever exercised the
        # earlier completeness check, never this one).
        recs = {
            "U": _fake_record([0.0, 0.0], [False, False], status="PASS-ZERO"),
            "g1": _fake_record([1.0, 0.0], [True, False]),
            "g2": _fake_record([0.0, 0.0], [False, False], status="PASS-ZERO"),
        }
        with self.assertRaisesRegex(AssertionError, "active"):
            ed_oracle_util.sensitivity_rank(recs, "fx5/pp_t")

    def test_underdetermined_matrix_marks_inconclusive(self):
        # fx3/pp_t expects active=(g+, g-). Both directions bear on the
        # SAME single canonical row (index 0) and are structurally zero
        # everywhere else -- note this is NOT the same as two directions
        # bearing on DIFFERENT rows: two disjoint-support columns would
        # give a (2, 2) diagonal matrix (well-conditioned, not
        # underdetermined -- "structurally zero on a union row
        # contributes exact 0.0" there, not a missing row). Here the
        # union bearing mask has exactly ONE True row, so the restricted
        # matrix is (1, 2): fewer rows than active columns, structurally
        # unable to have full column rank regardless of the raw values.
        recs = {
            "g+": _fake_record([5.0, 0.0, 0.0, 0.0],
                                [True, False, False, False]),
            "g-": _fake_record([3.0, 0.0, 0.0, 0.0],
                                [True, False, False, False]),
        }
        result = ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")
        self.assertEqual(result["status"], "INCONCLUSIVE")
        self.assertEqual(result["n_rows"], 1)

    def test_parallel_columns_low_sv_ratio_marks_inconclusive(self):
        # rows (4) >= columns (2), so this goes through the actual SVD
        # gate (not the underdetermined shortcut) -- exactly parallel
        # columns give sigma_min == 0 exactly, sv_ratio == 0.0 < floor.
        recs = {
            "g+": _fake_record([1.0, 2.0, 3.0, 4.0],
                                [True, True, True, True], delta_rich=1e-8),
            "g-": _fake_record([2.0, 4.0, 6.0, 8.0],
                                [True, True, True, True], delta_rich=1e-8),
        }
        result = ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")
        self.assertEqual(result["status"], "INCONCLUSIVE")
        self.assertLess(result["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)

    def test_sigma_min_below_delta_rich_floor_marks_inconclusive(self):
        # well-separated (orthogonal) columns -- sv_ratio == 1.0, clears
        # SENS_SV_FLOOR easily -- but delta_rich_max is deliberately huge
        # relative to sigma_min, so the OTHER gate (100*delta_rich) fails.
        recs = {
            "g+": _fake_record([1.0, 0.0], [True, True], delta_rich=1.0),
            "g-": _fake_record([0.0, 1.0], [True, True], delta_rich=1e-8),
        }
        result = ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")
        self.assertEqual(result["status"], "INCONCLUSIVE")
        self.assertGreaterEqual(result["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)
        self.assertLess(result["sigma_min"], 100.0 * result["delta_rich_max"])

    def test_well_conditioned_active_matrix_passes(self):
        recs = {
            "g+": _fake_record([1.0, 0.0, 0.5, 0.0],
                                [True, True, True, True], delta_rich=1e-8),
            "g-": _fake_record([0.0, 1.0, 0.0, 0.5],
                                [True, True, True, True], delta_rich=1e-8),
        }
        result = ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")
        self.assertEqual(result["status"], "PASS")
        self.assertEqual(result["active"], ["g+", "g-"])
        self.assertEqual(result["null"], [])

    def test_mismatched_grid_length_across_directions_raises(self):
        recs = {
            "g+": _fake_record([1.0, 0.0], [True, True]),
            "g-": _fake_record([1.0, 0.0, 0.0], [True, True, True]),
        }
        with self.assertRaises(ValueError):
            ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")

    def test_mismatched_pred_full_bearing_mask_shape_within_record_raises(self):
        # a single direction whose OWN pred_full/bearing_mask disagree in
        # length -- the within-record branch (distinct from the
        # cross-direction grid-length mismatch above).
        bad = _fake_record([1.0, 0.0], [True, True])
        bad["pred_full"] = np.array([1.0, 0.0, 3.0])   # now len 3 vs mask's 2
        recs = {
            "g+": bad,
            "g-": _fake_record([1.0, 0.0], [True, True]),
        }
        with self.assertRaises(ValueError):
            ed_oracle_util.sensitivity_rank(recs, "fx3/pp_t")

    def test_null_direction_grid_length_mismatch_raises(self):
        # review should_fix (round 2): the shape/length validation must
        # cover EVERY supplied direction, not just the active columns --
        # here fx5/pp_t's null U anchor has a DIFFERENT canonical-grid
        # length than g1/g2, which must be caught even though U itself
        # never enters the matrix.
        recs = {
            "U": _fake_record([0.0, 0.0, 0.0], [False, False, False],
                               status="PASS-ZERO"),
            "g1": _fake_record([1.0, 0.0], [True, False]),
            "g2": _fake_record([0.0, 1.0], [False, True]),
        }
        with self.assertRaises(ValueError):
            ed_oracle_util.sensitivity_rank(recs, "fx5/pp_t")


class TestCampaignAggregationRobustness(unittest.TestCase):
    """Regression coverage (review should_fix, two rounds) for
    ``_campaign_verdict``'s completeness/non-vacuous-``all()`` guards,
    using a synthetic ``state`` dict directly (``_campaign_verdict`` is
    pure logic over its ``state`` argument -- exercising it does not need
    to run any ED/solver machinery). Round 2 of review found that the
    original version of this class built states under FABRICATED labels
    (``"fx5/0"`` etc.), which only exercised the per-fixture non-emptiness
    guard, never ``_campaign_verdict``'s OWN completeness assertion added
    in that same round (a partial-but-nonempty state under the REAL
    expected labels would previously have slipped past the old
    non-emptiness-only check). States here are now built over the exact
    same ``_CAMPAIGN_EXPECTED_GRANULES``/``_CAMPAIGN_EXPECTED_SENSITIVITY``
    label sets ``_campaign_state`` itself would produce.
    """

    def _synthetic_state(self, overrides=None, omit_fixtures=(),
                          omit_labels=()):
        """A synthetic, all-clean-by-default ``state`` over the REAL
        expected granule/sensitivity label sets. ``overrides``: dict
        mapping a specific label to a specific status. ``omit_fixtures``:
        drop every label for the given fixture name(s) entirely
        (simulates a fixture's evidence never aggregating).
        ``omit_labels``: drop specific individual labels (simulates a
        partial, non-fixture-aligned aggregation gap)."""
        overrides = overrides or {}

        def _build(expected_labels):
            out = {}
            for lbl in expected_labels:
                fixture = lbl.split("/", 1)[0]
                if fixture in omit_fixtures or lbl in omit_labels:
                    continue
                out[lbl] = dict(status=overrides.get(lbl, "PASS"))
            return out

        return dict(granules=_build(_CAMPAIGN_EXPECTED_GRANULES),
                    sensitivity=_build(_CAMPAIGN_EXPECTED_SENSITIVITY))

    def test_missing_fixture_evidence_raises(self):
        # fx3 has NO granules/sensitivity entries at all (simulates a
        # broken aggregation loop) -- must raise via the completeness
        # assertion, never read as a vacuous PASS.
        state = self._synthetic_state(omit_fixtures=("fx3",))
        with self.assertRaises(AssertionError):
            _campaign_verdict(state)

    def test_partial_incomplete_granules_state_raises(self):
        # review should_fix (round 2): a NONEMPTY-per-fixture but
        # incomplete state (missing just ONE granule label, keeping every
        # other real label -- including every sensitivity label -- intact)
        # must also raise -- the old per-fixture-non-emptiness-only guard
        # would have missed this.
        one_granule_label = next(iter(_CAMPAIGN_EXPECTED_GRANULES))
        state = self._synthetic_state(omit_labels=(one_granule_label,))
        with self.assertRaises(AssertionError):
            _campaign_verdict(state)

    def test_partial_incomplete_sensitivity_state_raises(self):
        # review should_fix (round 3): the round-2 tests above only ever
        # exercised the GRANULE completeness assertion (the sensitivity
        # side was always either fully present or fully absent together
        # with the granules) -- nothing independently proved the
        # SENSITIVITY key-set assertion actually fires. Here the granule
        # set stays FULLY complete and only one sensitivity label is
        # dropped.
        one_sensitivity_label = next(iter(_CAMPAIGN_EXPECTED_SENSITIVITY))
        state = self._synthetic_state(omit_labels=(one_sensitivity_label,))
        self.assertEqual(set(state["granules"]), _CAMPAIGN_EXPECTED_GRANULES)
        with self.assertRaises(AssertionError):
            _campaign_verdict(state)

    def test_fixture_status_fail_on_a_failing_granule(self):
        one_fx5_granule = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_GRANULES
            if lbl.startswith("fx5/"))
        state = self._synthetic_state(overrides={one_fx5_granule: "FAIL"})
        verdict = _campaign_verdict(state)
        self.assertEqual(verdict["fixtures"]["fx5"], "FAIL")
        self.assertEqual(verdict["fixtures"]["fx3"], "PASS")
        self.assertEqual(verdict["campaign"], "FAIL")

    def test_fixture_status_inconclusive_on_a_granule_power_shortfall(self):
        # review must_fix (round 3): a granule-level INCONCLUSIVE (a
        # power/resolution shortfall, per the design doc -- NOT a
        # demonstrated numeric mismatch) must report as fixture
        # "INCONCLUSIVE", distinct from "FAIL".
        one_fx5_granule = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_GRANULES
            if lbl.startswith("fx5/"))
        state = self._synthetic_state(
            overrides={one_fx5_granule: "INCONCLUSIVE"})
        verdict = _campaign_verdict(state)
        self.assertEqual(verdict["fixtures"]["fx5"], "INCONCLUSIVE")
        self.assertEqual(verdict["campaign"], "INCONCLUSIVE")

    def test_fixture_status_inconclusive_on_a_rank_gate_failure(self):
        # the design doc's OWN named case: "a rank failure marks the
        # fixture INCONCLUSIVE" -- distinct from FAIL.
        one_fx5_sensitivity = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_SENSITIVITY
            if lbl.startswith("fx5/"))
        state = self._synthetic_state(
            overrides={one_fx5_sensitivity: "INCONCLUSIVE"})
        verdict = _campaign_verdict(state)
        self.assertEqual(verdict["fixtures"]["fx5"], "INCONCLUSIVE")
        self.assertEqual(verdict["campaign"], "INCONCLUSIVE")

    def test_fail_takes_precedence_over_inconclusive(self):
        granule_lbl = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_GRANULES
            if lbl.startswith("fx5/"))
        sensitivity_lbl = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_SENSITIVITY
            if lbl.startswith("fx5/"))
        state = self._synthetic_state(overrides={
            granule_lbl: "FAIL", sensitivity_lbl: "INCONCLUSIVE"})
        verdict = _campaign_verdict(state)
        self.assertEqual(verdict["fixtures"]["fx5"], "FAIL")
        self.assertEqual(verdict["campaign"], "FAIL")

    def test_skipped_granule_failure_sensitivity_placeholder_status(self):
        # review should_fix (round 4): the tests above inject FAIL/
        # INCONCLUSIVE only at the granule level (leaving the synthetic
        # sensitivity side "PASS"), never exercising the REAL
        # ``_campaign_state``-produced "SKIPPED-GRANULE-FAILURE"
        # sensitivity placeholder a channel gets when its OWN granules
        # already have a FAIL/INCONCLUSIVE (``sensitivity_rank`` refuses
        # non-PASS input by design, so ``_campaign_state`` never calls it
        # in that case -- see its docstring). Mirror that exact
        # combination directly: a FAIL granule paired with its channel's
        # sensitivity result reading "SKIPPED-GRANULE-FAILURE" (not
        # "PASS"), and separately an INCONCLUSIVE granule paired with the
        # same placeholder -- both must still resolve to the correct
        # fixture/campaign status via the precedence rule.
        granule_lbl = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_GRANULES
            if lbl.startswith("fx5/"))
        sensitivity_lbl = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_SENSITIVITY
            if lbl.startswith("fx5/"))

        state_fail = self._synthetic_state(overrides={
            granule_lbl: "FAIL",
            sensitivity_lbl: "SKIPPED-GRANULE-FAILURE"})
        verdict_fail = _campaign_verdict(state_fail)
        self.assertEqual(verdict_fail["fixtures"]["fx5"], "FAIL")
        self.assertEqual(verdict_fail["campaign"], "FAIL")

        state_inconclusive = self._synthetic_state(overrides={
            granule_lbl: "INCONCLUSIVE",
            sensitivity_lbl: "SKIPPED-GRANULE-FAILURE"})
        verdict_inconclusive = _campaign_verdict(state_inconclusive)
        self.assertEqual(verdict_inconclusive["fixtures"]["fx5"],
                          "INCONCLUSIVE")
        self.assertEqual(verdict_inconclusive["campaign"], "INCONCLUSIVE")

    def test_all_clean_synthetic_state_passes(self):
        state = self._synthetic_state()
        verdict = _campaign_verdict(state)
        self.assertEqual(verdict["fixtures"]["fx5"], "PASS")
        self.assertEqual(verdict["fixtures"]["fx3"], "PASS")
        self.assertEqual(verdict["campaign"], "PASS")

    def test_unrecognized_status_value_raises(self):
        one_fx5_granule = next(
            lbl for lbl in _CAMPAIGN_EXPECTED_GRANULES
            if lbl.startswith("fx5/"))
        state = self._synthetic_state(
            overrides={one_fx5_granule: "NOT-A-REAL-STATUS"})
        with self.assertRaises(AssertionError):
            _campaign_verdict(state)


# ---------------------------------------------------------------------------
# Task 9: sensitivity diagnostic and campaign verdict (#151) -- aggregates
# EVERY granule record and EVERY (fixture, channel) sensitivity_rank result
# produced by Tasks 6-8, by CALLING their already-``functools.lru_cache``'d
# case helpers (never recomputing anything, never depending on which test
# ran first). ``_campaign_state`` is itself ``lru_cache``'d so every test
# method in ``TestCampaignVerdict`` (and the verdict-writing step outside
# this module) sees the SAME aggregation.
# ---------------------------------------------------------------------------

_CAMPAIGN_FX5_PH_DIRECTIONS = ("U", "g1", "g2")
_CAMPAIGN_FX5_PP_DIRECTIONS = ("U", "g1", "g2")
_CAMPAIGN_FX3_PH_DIRECTIONS = ("g+", "g-")
_CAMPAIGN_FX3_PP_DIRECTIONS = ("g+", "g-")

# The a-priori COMPLETE set of campaign granule labels
# (``"fixture/direction/channel"``) and (fixture, channel) sensitivity
# labels -- review finding (must_fix): ``_campaign_verdict``'s per-fixture
# ``all(...)`` over a label-prefix filter is VACUOUSLY True on an EMPTY
# filtered collection, so a fixture whose evidence silently failed to
# aggregate (a bug in ``_collect``'s loop, a renamed channel key, ...)
# would read as PASS with zero granules checked. Declaring the expected
# COMPLETE key set up front and asserting it against what
# ``_campaign_state`` actually built closes that gap: an incomplete
# aggregation now raises loudly here, before ``_campaign_verdict`` ever
# runs its (now guaranteed non-vacuous) per-fixture ``all(...)``.
_CAMPAIGN_EXPECTED_GRANULES = frozenset(
    "{}/{}/{}".format(fixture, direction, channel)
    for fixture, directions, channels in (
        ("fx5", _CAMPAIGN_FX5_PH_DIRECTIONS, ("S", "C")),
        ("fx5", _CAMPAIGN_FX5_PP_DIRECTIONS, ("pp_s", "pp_t")),
        ("fx3", _CAMPAIGN_FX3_PH_DIRECTIONS, ("S", "C")),
        ("fx3", _CAMPAIGN_FX3_PP_DIRECTIONS, ("pp_s", "pp_t")),
    )
    for direction in directions
    for channel in channels)
# The expected (fixture, channel) sensitivity labels are exactly the keys
# ``sensitivity_rank``'s own a-priori table declares (single source of
# truth -- never duplicated by hand here).
_CAMPAIGN_EXPECTED_SENSITIVITY = frozenset(
    ed_oracle_util.SENSITIVITY_EXPECTED_ACTIVE)


@functools.lru_cache(maxsize=1)
def _campaign_state():
    """Deterministic aggregation of every campaign granule record plus
    every (fixture, channel) ``sensitivity_rank`` result. Builds
    ``by_channel`` (``"fx5/S" -> {direction: record}``, the exact input
    shape ``ed_oracle_util.sensitivity_rank`` consumes) alongside the flat
    ``granules`` dict (``"fx5/U/S" -> record``) used for the granule rule.
    ``sensitivity_rank`` is only invoked for a channel whose OWN granules
    are all PASS/PASS-ZERO already (it refuses non-PASS input by design --
    see its docstring); a channel with a FAIL/INCONCLUSIVE granule instead
    gets a ``status="SKIPPED-GRANULE-FAILURE"`` sensitivity placeholder so
    the granule rule alone (not a raised exception) is what fails that
    fixture, keeping the two rules independently diagnosable.

    Asserts the resulting ``granules``/``sensitivity`` key sets exactly
    match ``_CAMPAIGN_EXPECTED_GRANULES``/``_CAMPAIGN_EXPECTED_SENSITIVITY``
    (review must_fix) -- an incomplete aggregation raises here rather than
    silently letting a downstream ``all(...)`` read an empty fixture slice
    as PASS.
    """
    granules = {}
    by_channel = {}

    def _collect(fixture, direction, chan_map):
        for chan, rec in chan_map.items():
            key = "{}/{}/{}".format(fixture, direction, chan)
            # Review finding (should_fix, round 3): the key is
            # SYNTHESIZED from the loop's own (fixture, direction, chan)
            # arguments, never cross-checked against the record's OWN
            # ``label`` field (every ``adjudicate_granule`` record carries
            # one, e.g. ``"fx5/{direction}/S"``) -- a swapped/mislabeled
            # record from a future case-helper edit could otherwise be
            # silently adjudicated under the WRONG fixture/direction/
            # channel label while still satisfying every completeness
            # check downstream (those only check KEY sets, never that a
            # key's VALUE actually is what it claims to be).
            if rec["label"] != key:
                raise AssertionError(
                    "_campaign_state._collect: record's own label {!r} != "
                    "the synthesized aggregation key {!r} -- a mislabeled "
                    "or swapped record must never be silently aggregated "
                    "under the wrong key".format(rec["label"], key))
            granules[key] = rec
            by_channel.setdefault(
                "{}/{}".format(fixture, chan), {})[direction] = rec

    for d in _CAMPAIGN_FX5_PH_DIRECTIONS:
        _collect("fx5", d, case_fx5_unit_direction(d))
    for d in _CAMPAIGN_FX5_PP_DIRECTIONS:
        _collect("fx5", d, case_fx5_pp_unit_direction(d))
    for d in _CAMPAIGN_FX3_PH_DIRECTIONS:
        _collect("fx3", d, case_fx3_unit_direction(d))
    for d in _CAMPAIGN_FX3_PP_DIRECTIONS:
        _collect("fx3", d, case_fx3_pp_unit_direction(d))

    if set(granules) != _CAMPAIGN_EXPECTED_GRANULES:
        raise AssertionError(
            "_campaign_state: incomplete/unexpected granule aggregation -- "
            "expected {} but got {} (missing {}, unexpected {})".format(
                sorted(_CAMPAIGN_EXPECTED_GRANULES), sorted(granules),
                sorted(_CAMPAIGN_EXPECTED_GRANULES - set(granules)),
                sorted(set(granules) - _CAMPAIGN_EXPECTED_GRANULES)))
    if set(by_channel) != _CAMPAIGN_EXPECTED_SENSITIVITY:
        raise AssertionError(
            "_campaign_state: incomplete/unexpected (fixture, channel) "
            "aggregation -- expected {} but got {}".format(
                sorted(_CAMPAIGN_EXPECTED_SENSITIVITY), sorted(by_channel)))

    sensitivity = {}
    for label, recs in by_channel.items():
        if all(r["status"] in ("PASS", "PASS-ZERO") for r in recs.values()):
            sensitivity[label] = ed_oracle_util.sensitivity_rank(recs, label)
        else:
            sensitivity[label] = dict(
                label=label, active=[], null=[], n_rows=0,
                singular_values=np.array([]), sv_ratio=float("nan"),
                sigma_min=float("nan"), sigma_max=float("nan"),
                delta_rich_max=float("nan"),
                status="SKIPPED-GRANULE-FAILURE")

    return dict(granules=granules, by_channel=by_channel,
                sensitivity=sensitivity)


def _campaign_verdict(state):
    """Fixture rule (brief): a fixture is PASS iff none of its granules is
    FAIL/INCONCLUSIVE AND every one of its (fixture, channel)
    sensitivity_rank gates is PASS. Campaign rule: PASS iff every fixture
    PASSes. The Task-3 S_bond/C_bond null-direction finding is OUTSIDE
    this verdict by design (a separately recorded finding, not a
    granule) -- see TestNullDirectionSolverSide and the verdict document.

    Each fixture (and the campaign as a whole) gets a TRI-STATE status
    string -- ``"PASS"``, ``"INCONCLUSIVE"``, or ``"FAIL"`` -- not a
    collapsed boolean (review must_fix: the design doc explicitly names
    INCONCLUSIVE as its OWN distinct outcome for a rank-gate failure, "a
    rank failure marks the fixture INCONCLUSIVE", which is a genuinely
    different scientific finding from a granule FAIL -- a demonstrated
    numeric mismatch vs. insufficient power/rank to see one either way --
    and collapsing both to a bare boolean would erase that distinction
    from the printed verdict). Precedence per fixture: any granule
    ``FAIL`` -> ``"FAIL"``; else any granule ``INCONCLUSIVE`` or any
    sensitivity-rank ``INCONCLUSIVE`` -> ``"INCONCLUSIVE"``; else
    (everything ``PASS``/``PASS-ZERO``/``PASS``) -> ``"PASS"``. Campaign:
    the same precedence over the fixtures' own statuses (any fixture
    ``FAIL`` -> campaign ``"FAIL"``; else any fixture ``INCONCLUSIVE`` ->
    campaign ``"INCONCLUSIVE"``; else ``"PASS"``) -- matching the design
    doc's "the campaign PASSES ... iff every fixture passes; any other
    outcome keeps the gate" (both FAIL and INCONCLUSIVE keep the gate;
    only the PRINTED label distinguishes which).

    ``state``'s ``granules``/``sensitivity`` key sets are asserted here
    directly against ``_CAMPAIGN_EXPECTED_GRANULES``/
    ``_CAMPAIGN_EXPECTED_SENSITIVITY`` (review should_fix, round 2:
    ``_campaign_state`` already enforces this on its OWN return value,
    but ``_campaign_verdict`` is a separate, independently callable/
    testable function -- relying only on its caller's discipline would
    let a partial-but-nonempty ``state`` built any other way slip past
    the per-fixture non-emptiness check below and still read as PASS).
    """
    if set(state["granules"]) != _CAMPAIGN_EXPECTED_GRANULES:
        raise AssertionError(
            "_campaign_verdict: incomplete/unexpected granule state -- "
            "expected {} but got {}".format(
                sorted(_CAMPAIGN_EXPECTED_GRANULES),
                sorted(state["granules"])))
    if set(state["sensitivity"]) != _CAMPAIGN_EXPECTED_SENSITIVITY:
        raise AssertionError(
            "_campaign_verdict: incomplete/unexpected sensitivity state -- "
            "expected {} but got {}".format(
                sorted(_CAMPAIGN_EXPECTED_SENSITIVITY),
                sorted(state["sensitivity"])))

    fixtures = {}
    for fixture in ("fx5", "fx3"):
        granule_items = [r for label, r in state["granules"].items()
                          if label.startswith(fixture + "/")]
        sensitivity_items = [s for label, s in state["sensitivity"].items()
                              if label.startswith(fixture + "/")]
        if not granule_items or not sensitivity_items:
            raise AssertionError(
                "_campaign_verdict: fixture {!r} has NO granules ({}) or NO "
                "sensitivity results ({}) -- an empty slice must never be "
                "read as a vacuous PASS".format(
                    fixture, len(granule_items), len(sensitivity_items)))
        # Defensive: the precedence logic below treats "not FAIL, not
        # INCONCLUSIVE" as PASS-worthy -- an unrecognized status string
        # (contract drift in adjudicate_granule/sensitivity_rank) must
        # not silently fall through that else branch as a PASS.
        bad_granule = [r["status"] for r in granule_items
                       if r["status"] not in
                       ("PASS", "PASS-ZERO", "FAIL", "INCONCLUSIVE")]
        bad_sensitivity = [s["status"] for s in sensitivity_items
                            if s["status"] not in
                            ("PASS", "INCONCLUSIVE", "SKIPPED-GRANULE-FAILURE")]
        if bad_granule or bad_sensitivity:
            raise AssertionError(
                "_campaign_verdict: fixture {!r} has unrecognized status "
                "value(s) -- granule {} / sensitivity {}".format(
                    fixture, bad_granule, bad_sensitivity))
        if any(r["status"] == "FAIL" for r in granule_items):
            status = "FAIL"
        elif (any(r["status"] == "INCONCLUSIVE" for r in granule_items)
                or any(s["status"] != "PASS" for s in sensitivity_items)):
            # A sensitivity_rank result that is anything other than
            # "PASS" (INCONCLUSIVE from a failed rank gate, or the
            # "SKIPPED-GRANULE-FAILURE" placeholder _campaign_state emits
            # when a channel's own granules already have a FAIL/
            # INCONCLUSIVE -- the granule branch above already covers
            # that FAIL case, so what reaches here is genuinely a rank/
            # power finding) is, per the design doc, a rank/power
            # finding, not a demonstrated numeric mismatch.
            status = "INCONCLUSIVE"
        else:
            status = "PASS"
        fixtures[fixture] = status

    if any(v == "FAIL" for v in fixtures.values()):
        campaign = "FAIL"
    elif any(v == "INCONCLUSIVE" for v in fixtures.values()):
        campaign = "INCONCLUSIVE"
    else:
        campaign = "PASS"
    return dict(fixtures=fixtures, campaign=campaign)


def _print_campaign_verdict_table(state, verdict):
    print("\n=== CAMPAIGN VERDICT TABLE (#151 Task 9) ===")
    print("-- granules --")
    for label in sorted(state["granules"]):
        r = state["granules"][label]
        print("  {:14s} status={:10s} delta_rich={:.3e} delta_nmat={:.3e} "
              "tol={:.3e} max_signal={:.3e} n_failures={}".format(
                  label, r["status"], r["delta_rich"], r["delta_nmat"],
                  r["tol"], r["max_signal"], len(r["failures"])))
    print("-- sensitivity_rank (fixture/channel) --")
    for label in sorted(state["sensitivity"]):
        s = state["sensitivity"][label]
        print("  {:11s} active={} null={} rows={} sv_ratio={:.3e} "
              "sigma_min={:.3e} 100*delta_rich_max={:.3e} status={}".format(
                  label, s["active"], s["null"], s["n_rows"], s["sv_ratio"],
                  s["sigma_min"], 100.0 * s["delta_rich_max"], s["status"]))
    print("-- fixtures --")
    for fixture in ("fx5", "fx3"):
        print("  {}: {}".format(fixture, verdict["fixtures"][fixture]))
    print("CAMPAIGN VERDICT: {}".format(verdict["campaign"]))
    print("(Task-3 S_bond/C_bond null-direction finding: OUTSIDE this "
          "verdict, recorded separately -- see TestNullDirectionSolverSide "
          "and the verdict document.)\n")


class TestCampaignVerdict(unittest.TestCase):
    """Task 9 Step 1: the campaign-level verdict, aggregating every
    granule record (Tasks 6-8) and every (fixture, channel)
    ``sensitivity_rank`` result via ``_campaign_state`` (this module's own
    ``lru_cache``'d wrapper around the case helpers -- no recomputation, no
    test-order dependence). THE VERDICT MATTERS MORE THAN GREEN: as of this
    task, every granule and every sensitivity gate PASSES (verified
    directly, not assumed -- fx3's g+/g- columns turned out NOT to be
    parallel on the full canonical grid despite their audit TUPLES being
    bit-identical by inversion symmetry: their per-cell ``pred_full``
    values differ by up to ~0.43, giving sv_ratio in [0.87, 0.98] on all
    four fx3 channels, comfortably above SENS_SV_FLOOR), so this class is
    expected green; a future FAIL/INCONCLUSIVE here is a real campaign
    finding and must be reported, not adjusted away.

    Deferred (review should_fix, investigated and NOT implemented): the
    ``lru_cache``'d case helpers return their record dicts/arrays BY
    REFERENCE, so ``_campaign_state`` does not itself guarantee immutable
    read-only results -- a consumer that mutated a returned record in
    place could change what a LATER call sees. Checked directly: no test
    anywhere in this module mutates a record (every consumption reads a
    field, e.g. ``r["status"]``, ``r["failures"][:5]``); adding a
    subprocess-based cross-order regression would cost ~4x this module's
    already-3-minute fx3 wall time for a hazard nothing here currently
    exercises. Left as a documented invariant (records are read-only by
    convention, not by enforcement) rather than adding that machinery,
    matching this campaign's established practice of deferring optional/
    expensive robustness items with an explicit reason (see e.g. Task 4's
    and Task 5's reports).
    """

    def test_no_granule_fails_or_is_inconclusive(self):
        state = _campaign_state()
        bad = [(label, r["status"]) for label, r in state["granules"].items()
               if r["status"] not in ("PASS", "PASS-ZERO")]
        self.assertEqual(
            bad, [], "granule(s) not PASS/PASS-ZERO: {}".format(bad))

    def test_sensitivity_rank_gates(self):
        state = _campaign_state()
        for label, s in state["sensitivity"].items():
            with self.subTest(label=label):
                self.assertEqual(
                    s["status"], "PASS",
                    "sensitivity_rank[{}] status={} sv_ratio={:.3e} "
                    "sigma_min={:.3e} 100*delta_rich_max={:.3e} active={} "
                    "null={}".format(
                        label, s["status"], s["sv_ratio"], s["sigma_min"],
                        100.0 * s["delta_rich_max"], s["active"], s["null"]))

    def test_campaign_verdict(self):
        state = _campaign_state()
        verdict = _campaign_verdict(state)
        _print_campaign_verdict_table(state, verdict)
        for fixture, status in verdict["fixtures"].items():
            with self.subTest(fixture=fixture):
                self.assertEqual(
                    status, "PASS",
                    "fixture {} verdict status={} (not PASS) -- see the "
                    "printed table above".format(fixture, status))
        self.assertEqual(
            verdict["campaign"], "PASS",
            "CAMPAIGN VERDICT: {} -- see the printed table above".format(
                verdict["campaign"]))


if __name__ == "__main__":
    unittest.main()
