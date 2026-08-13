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


V1 = 0.02      # base coupling used throughout the campaign
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
    tail correction) Matsubara sum, computed by two independent code
    paths (``bond_bubble`` vs ``RPA._calc_chi0q``), so this pin does not
    replace the ED frame-map comparison (Pin 2b) -- it is solver-internal
    self-consistency only.
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


class TestPin2bFrameMap(ApproxTestCase):
    """Pin 2b (LOAD-BEARING, Task 4 review): ``ed_to_solver_bond_map``,
    zero coupling, FULL (m, m') matrix, both fixtures. ``bond_bubble``
    (RAW, no tail correction, per its own documented limitation) is
    compared against ``SectorED.bond_correlator`` (exact Lehmann) through
    the mapped index; the measured finite-Nmat distance eps2 is a
    RECORDED DIAGNOSTIC with a fixed safety ceiling (``< 1e-6``), never an
    input to any tolerance formula. A two-point Nmat-Richardson
    (``order=1``, matching ``bond_bubble``'s raw O(1/Nmat) convergence --
    see ``_richardson_nmat``) brings eps2 comfortably below the ceiling at
    a small Nmat base (512), rather than requiring an impractically large
    single Nmat (the raw O(1/Nmat) law measured a coefficient that would
    need Nmat ~ 2e5 for 1e-6 without extrapolation).

    Both same-spin diagonal blocks (sigma=sigma'=0 and sigma=sigma'=1) are
    checked against the SAME mapped matrix (spin-independent free
    Hamiltonian), and the cross-spin blocks are checked to be exactly
    null -- the ORIENTATION STRUCTURE (which off-diagonal (m, m') block is
    the conjugate of which) is pinned by construction: the full ND x ND
    matrix is compared elementwise, not reduced to a scalar or a diagonal
    slice, so a swapped/transposed/conjugated off-diagonal block cannot
    hide behind an accidental cancellation.
    """

    def _check(self, fx, decls, n1):
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

        def _chibar(nmat):
            green = free_green(fx, nmat)
            return bond_bubble(green, bond_set, beta=fx.beta)[:, 0, 0]

        chibar_rich = _richardson_nmat(_chibar, n1, order=1)

        sector = ed_oracle_util.SectorED(fx)
        xph = sector.bond_correlator(channels)   # (L, 2, 2, n_i, n_i)

        eps2 = 0.0
        for (s1, s2) in ((0, 0), (0, 1), (1, 0), (1, 1)):
            with self.subTest(spins=(s1, s2)):
                xs = xph[:, s1, s2][:, smap][:, :, smap]  # (L, ND, ND)
                expected = chibar_rich if s1 == s2 else np.zeros_like(chibar_rich)
                diff = float(np.abs(xs - expected).max())
                # review finding: max(eps2, diff) silently discards a NaN
                # diff (NaN compares False to everything, so Python's max
                # keeps the OTHER, possibly-zero, operand) -- which would
                # let this load-bearing pin pass on a broken (NaN) result.
                # Fail loudly instead.
                self.assertTrue(
                    np.isfinite(diff),
                    "Pin2b: non-finite discrepancy (spins={}, diff={})"
                    .format((s1, s2), diff))
                eps2 = max(eps2, diff)
                print("Pin2b eps2 (fixture L={}, norb={}, spins={}): {:.3e}"
                      .format(fx.L, fx.norb, (s1, s2), diff))
        self.assertLess(eps2, 1e-6)

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


if __name__ == "__main__":
    unittest.main()
