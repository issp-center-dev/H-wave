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
import unittest

import numpy as np

from hwave.solver.bond_channels import bare_bond_vertices, resolve_interactions
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
    coulomb_inter = {
        ((1, 0, 0), (0, 1)): V1 + 1j * eps,
        ((-1, 0, 0), (1, 0)): V1 - 1j * eps,
    }
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


if __name__ == "__main__":
    unittest.main()
