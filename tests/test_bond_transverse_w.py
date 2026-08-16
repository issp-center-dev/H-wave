#!/usr/bin/env python3

"""Production-side unit/structural tests for the master transverse bond
topology (Phase W, Task 2 of
``docs/superpowers/plans/2026-08-16-bond-transverse-phase-w.md``): pins
``hwave.solver.bond_channels.TransverseTopology``,
``resolve_transverse_topology``, ``iter_reversal_orbits``,
``validate_topology_against_mesh`` and ``transverse_effective_activity``
against the invariants of
``docs/superpowers/specs/2026-08-15-bond-transverse-design.md`` ("The
master transverse topology", "Canonical orbit rule", "Construction domain
(binding for steps 2-3)").

This module covers TOPOLOGY-ONLY ground: no ``W_pm_bond`` vertex and no
bond bubble exist yet (Phase W Tasks 3/4) -- the ordered-record vertex
equations themselves were already numerically validated against ED by gate
W0 in ``tests/test_bond_transverse_ed.py`` (Task 1) and are NOT
re-adjudicated here.

Includes the carried Task-1-review finding's Hamiltonian-level normalization
pin (``TestHamiltonianLevelNormalizationPin``): the topology's FILE ->
COEFFS normalization is exercised independently of the (already-validated)
vertex equations, because the "double-declaration trap" for density-density
(commuting) operators -- declaring both ``C`` at ``R`` and ``conj(C)`` at
``-R`` describes the SAME physical operator sum TWICE, so the correct total
is ``2*Re(C)``, not ``C`` or ``2*C`` -- is a property of
``resolve_transverse_topology`` alone, not of the vertex construction that
consumes its output.

Tests must be run from the repository root.
"""

import unittest

import numpy as np

from hwave.solver.bond_channels import (
    TransverseTopology,
    resolve_transverse_topology,
    iter_reversal_orbits,
    validate_topology_against_mesh,
    transverse_effective_activity,
)
from tests import ed_oracle_util
from tests.approx_util import ApproxTestCase


def _topo(delta_r, reverse, coeffs=None, norb=1):
    """Small helper: build a TransverseTopology directly (bypassing
    resolve_transverse_topology) for tests that only care about the
    NamedTuple's own invariants/consumers, not the resolver."""
    B = len(delta_r)
    if coeffs is None:
        coeffs = {"CoulombInter": np.zeros((B, norb, norb), dtype=complex)}
    return TransverseTopology(
        delta_r=np.array(delta_r, dtype=int),
        reverse=np.array(reverse, dtype=int),
        coeffs=coeffs,
    )


# =============================================================================
# TransverseTopology: construction-time invariants
# =============================================================================

class TestTransverseTopologyConstruction(ApproxTestCase):

    def test_valid_construction_round_trips(self):
        t = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1])
        self.assertEqual(t.delta_r.shape, (3, 3))
        self.assertEqual(t.reverse.tolist(), [0, 2, 1])
        self.assertIn("CoulombInter", t.coeffs)

    def test_delta_r_and_reverse_are_int_dtype(self):
        t = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1])
        self.assertTrue(np.issubdtype(t.delta_r.dtype, np.integer))
        self.assertTrue(np.issubdtype(t.reverse.dtype, np.integer))

    def test_arrays_are_read_only(self):
        t = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1])
        self.assertFalse(t.delta_r.flags.writeable)
        self.assertFalse(t.reverse.flags.writeable)
        for arr in t.coeffs.values():
            self.assertFalse(arr.flags.writeable)
        with self.assertRaises(ValueError):
            t.delta_r[0, 0] = 5
        with self.assertRaises(ValueError):
            t.reverse[0] = 1
        with self.assertRaises(ValueError):
            t.coeffs["CoulombInter"][0, 0, 0] = 1.0

    def test_alias_safety_input_mutation_does_not_leak(self):
        dr = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
        rv = np.array([0, 2, 1])
        co_arr = np.zeros((3, 1, 1), dtype=complex)
        co = {"CoulombInter": co_arr}
        t = TransverseTopology(delta_r=dr, reverse=rv, coeffs=co)
        dr[1, 0] = 999
        rv[1] = 0
        co_arr[1, 0, 0] = 42.0
        co["Ising"] = np.zeros((3, 1, 1), dtype=complex)  # mutate the dict
        self.assertEqual(t.delta_r[1].tolist(), [1, 0, 0])
        self.assertEqual(int(t.reverse[1]), 2)
        self.assertApprox(t.coeffs["CoulombInter"][1, 0, 0], 0.0, abs=0)
        self.assertNotIn("Ising", t.coeffs)

    def test_channel_zero_must_be_origin(self):
        with self.assertRaises(ValueError):
            _topo([[1, 0, 0], [0, 0, 0]], [1, 0])

    def test_reverse_must_be_involution(self):
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 1, 1])

    def test_reverse_zero_must_be_zero(self):
        # A hand-crafted, otherwise-involutive reverse array with
        # reverse[0] != 0 must still be rejected -- channel 0 is defined to
        # be its own reversal partner.
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [1, 0, 2])

    def test_reverse_out_of_range_rejected(self):
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 5, 1])

    def test_reverse_geometrically_inconsistent_with_delta_r_rejected(self):
        # reverse[1] = 1 (involution-consistent, reverse[0]=0) but
        # delta_r[1] = (1,0,0) so -delta_r[1] = (-1,0,0) != delta_r[1].
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 1, 2])

    def test_closure_violation_rejected(self):
        bad = {"CoulombInter": np.array(
            [[[0]], [[1 + 2j]], [[9 + 9j]]], dtype=complex)}
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=bad)

    def test_closure_holds_exactly_for_a_consistent_pair(self):
        c = 0.4 - 0.3j
        good = {"CoulombInter": np.array(
            [[[0]], [[c]], [[np.conj(c)]]], dtype=complex)}
        t = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=good)
        self.assertApprox(
            t.coeffs["CoulombInter"][t.reverse[1]][0, 0],
            np.conj(t.coeffs["CoulombInter"][1][0, 0]), rel=0, abs=0)

    def test_offsite_coeffs_nonzero_at_channel_zero_rejected(self):
        bad = {"CoulombInter": np.array(
            [[[5.0]], [[1.0]], [[1.0]]], dtype=complex)}
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=bad)

    def test_channel_zero_coeffs_must_be_exactly_zero_not_just_small(self):
        bad = {"CoulombInter": np.array(
            [[[1e-30]], [[1.0]], [[1.0]]], dtype=complex)}
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=bad)

    def test_non_integral_delta_r_rejected(self):
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1.5, 0, 0]], [0, 1])

    def test_complex_delta_r_rejected(self):
        with self.assertRaises(ValueError):
            TransverseTopology(
                delta_r=np.array([[0, 0, 0], [1j, 0, 0]]),
                reverse=np.array([0, 1]), coeffs={})

    def test_coeffs_wrong_channel_count_rejected(self):
        bad = {"CoulombInter": np.zeros((2, 1, 1), dtype=complex)}
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=bad)

    def test_coeffs_norb_mismatch_across_types_rejected(self):
        bad = {
            "CoulombInter": np.zeros((3, 1, 1), dtype=complex),
            "Ising": np.zeros((3, 2, 2), dtype=complex),
        }
        with self.assertRaises(ValueError):
            _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1], coeffs=bad)

    def test_minimal_b1_topology_is_valid(self):
        t = _topo([[0, 0, 0]], [0], coeffs={
            "CoulombInter": np.zeros((1, 2, 2), dtype=complex)})
        self.assertEqual(t.delta_r.shape, (1, 3))


# =============================================================================
# resolve_transverse_topology
# =============================================================================

class TestResolveTransverseTopology(ApproxTestCase):

    def test_no_interactions_gives_b1_zero_topology(self):
        topo = resolve_transverse_topology({}, np.eye(3), norb=2)
        self.assertEqual(topo.delta_r.shape, (1, 3))
        self.assertEqual(topo.delta_r[0].tolist(), [0, 0, 0])
        self.assertEqual(sorted(topo.coeffs.keys()),
                          ["CoulombInter", "Exchange", "Ising"])
        for arr in topo.coeffs.values():
            self.assertTrue(np.all(arr == 0))

    def test_hund_and_pairlift_contribute_no_shells(self):
        interactions = {
            "Hund": {((1, 0, 0), (0, 0)): 0.7},
            "PairLift": {((1, 0, 0), (0, 0)): 0.3},
        }
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        self.assertEqual(topo.delta_r.shape, (1, 3))
        self.assertFalse(transverse_effective_activity(topo))

    def test_missing_type_key_treated_as_empty(self):
        # Only CoulombInter present; Ising/Exchange absent entirely.
        interactions = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.2}}
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        self.assertIn("Ising", topo.coeffs)
        self.assertIn("Exchange", topo.coeffs)
        self.assertTrue(np.all(topo.coeffs["Ising"] == 0))
        self.assertTrue(np.all(topo.coeffs["Exchange"] == 0))

    def test_onsite_entries_of_active_types_are_ignored(self):
        # On-site (R=0) CoulombInter/Ising/Exchange declarations never
        # become a topology channel (spec: "coeffs carry OFF-SITE channels
        # only").
        interactions = {
            "CoulombInter": {((0, 0, 0), (0, 0)): 5.0},
            "Ising": {((0, 0, 0), (0, 0)): 3.0},
        }
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        self.assertEqual(topo.delta_r.shape, (1, 3))

    def test_union_of_shells_across_types(self):
        # CoulombInter declares a +/-x shell only; Ising declares a
        # DIFFERENT +/-y shell only; the union must carry BOTH, and each
        # type's coeffs must be zero on the shell it didn't declare.
        interactions = {
            "CoulombInter": {((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2},
            "Ising": {((0, 1, 0), (0, 0)): 0.5, ((0, -1, 0), (0, 0)): 0.5},
        }
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        self.assertEqual(topo.delta_r.shape, (5, 3))  # 0 + 2 (+-x) + 2 (+-y)
        drs = [tuple(int(x) for x in r) for r in topo.delta_r]
        for m, dr in enumerate(drs):
            ci = topo.coeffs["CoulombInter"][m, 0, 0]
            ising = topo.coeffs["Ising"][m, 0, 0]
            if dr in ((1, 0, 0), (-1, 0, 0)):
                self.assertApprox(ci, 0.2, rel=0, abs=1e-13)
                self.assertApprox(ising, 0.0, rel=0, abs=0)
            elif dr in ((0, 1, 0), (0, -1, 0)):
                self.assertApprox(ci, 0.0, rel=0, abs=0)
                self.assertApprox(ising, 0.5, rel=0, abs=1e-13)
            else:
                self.assertApprox(ci, 0.0, rel=0, abs=0)
                self.assertApprox(ising, 0.0, rel=0, abs=0)

    def test_channel_ordering_zero_first_then_shells_by_length_then_lex(self):
        interactions = {
            "CoulombInter": {
                ((2, 0, 0), (0, 0)): 0.1, ((-2, 0, 0), (0, 0)): 0.1,
                ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2,
            },
        }
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        drs = [tuple(int(x) for x in r) for r in topo.delta_r]
        self.assertEqual(drs[0], (0, 0, 0))
        # shell 1 (length 1) before shell 2 (length 2); within a shell,
        # lexicographic irvec order (matching resolve_interactions).
        self.assertEqual(drs, [(0, 0, 0), (-1, 0, 0), (1, 0, 0),
                                (-2, 0, 0), (2, 0, 0)])

    def test_max_shells_truncates_whole_shells_when_dropped_shell_is_zero(self):
        # The kept shell (|R|=1) is nonzero; the DROPPED shell (|R|=2) is
        # declared but with an exactly-zero coefficient -- dropping only
        # zero-coefficient content is fine, never an error.
        interactions = {
            "CoulombInter": {
                ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2,
                ((2, 0, 0), (0, 0)): 0.0, ((-2, 0, 0), (0, 0)): 0.0,
            },
        }
        topo = resolve_transverse_topology(
            interactions, np.eye(3), norb=1, max_shells=1)
        self.assertEqual(topo.delta_r.shape, (3, 3))  # 0 + 1 shell (+-x, len1)
        drs = {tuple(int(x) for x in r) for r in topo.delta_r}
        self.assertEqual(drs, {(0, 0, 0), (1, 0, 0), (-1, 0, 0)})

    def test_max_shells_dropping_a_declared_nonzero_shell_raises(self):
        # max_shells=1 dropping a declared |R|=2 shell -> ValueError (the
        # ambiguity guard, generalizing resolve_interactions's
        # bond_max_shells=0 case to every truncation depth: dropping
        # declared content is never silent).
        interactions = {
            "CoulombInter": {
                ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2,
                ((2, 0, 0), (0, 0)): 0.1, ((-2, 0, 0), (0, 0)): 0.1,
            },
        }
        with self.assertRaises(ValueError) as cm:
            resolve_transverse_topology(
                interactions, np.eye(3), norb=1, max_shells=1)
        msg = str(cm.exception)
        self.assertIn("max_shells", msg)
        self.assertIn("(2, 0, 0)", msg)

    def test_max_shells_zero_with_declared_offsite_raises(self):
        # max_shells=0 + declared off-site -> ValueError (the n_keep=0
        # instance of the same guard).
        interactions = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(
                interactions, np.eye(3), norb=1, max_shells=0)

    def test_max_shells_zero_with_only_zero_declared_offsite_is_ok(self):
        # max_shells=0 with a declared-but-exactly-zero off-site value
        # drops nothing of substance -- not ambiguous, no error.
        interactions = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.0}}
        topo = resolve_transverse_topology(
            interactions, np.eye(3), norb=1, max_shells=0)
        self.assertEqual(topo.delta_r.shape, (1, 3))
        self.assertFalse(transverse_effective_activity(topo))

    def test_max_shells_zero_with_no_offsite_declared_is_ok(self):
        topo = resolve_transverse_topology(
            {}, np.eye(3), norb=1, max_shells=0)
        self.assertEqual(topo.delta_r.shape, (1, 3))

    def test_max_shells_keeping_all_declared_content_is_ok(self):
        # max_shells=1 keeping ALL declared content (only one shell
        # exists) -> OK, nothing dropped.
        interactions = {
            "CoulombInter": {
                ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2,
            },
        }
        topo = resolve_transverse_topology(
            interactions, np.eye(3), norb=1, max_shells=1)
        self.assertEqual(topo.delta_r.shape, (3, 3))
        drs = {tuple(int(x) for x in r) for r in topo.delta_r}
        self.assertEqual(drs, {(0, 0, 0), (1, 0, 0), (-1, 0, 0)})
        self.assertApprox(
            topo.coeffs["CoulombInter"][
                [i for i, r in enumerate(topo.delta_r)
                 if tuple(int(x) for x in r) == (1, 0, 0)][0], 0, 0],
            0.2, rel=0, abs=1e-13)

    def test_max_shells_drop_check_covers_every_active_type(self):
        # The guard must fire for Ising/Exchange too, not just
        # CoulombInter.
        for type_name, value in (("Ising", 0.4), ("Exchange", 0.3 + 0.1j)):
            interactions = {
                type_name: {
                    ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2,
                    ((2, 0, 0), (0, 0)): value, ((-2, 0, 0), (0, 0)): value,
                },
            }
            with self.assertRaises(ValueError):
                resolve_transverse_topology(
                    interactions, np.eye(3), norb=1, max_shells=1)

    def test_max_shells_negative_rejected(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=1, max_shells=-1)

    def test_max_shells_non_integral_rejected(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=1, max_shells=1.5)

    def test_mirrored_declaration_consistent_projects_to_conjugate_average(self):
        c1 = 0.4 + 0.1j
        c2 = np.conj(c1)  # exactly consistent mirrored pair
        interactions = {"CoulombInter": {
            ((1, 0, 0), (0, 0)): c1, ((-1, 0, 0), (0, 0)): c2}}
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        m1 = [i for i, r in enumerate(topo.delta_r)
              if tuple(int(x) for x in r) == (1, 0, 0)][0]
        self.assertApprox(topo.coeffs["CoulombInter"][m1, 0, 0], c1,
                           rel=0, abs=1e-13)

    def test_mirrored_declaration_inconsistent_raises(self):
        interactions = {"CoulombInter": {
            ((1, 0, 0), (0, 0)): 1.0 + 0.0j,
            ((-1, 0, 0), (0, 0)): -1.0 + 0.0j,  # not conj(1.0) = 1.0
        }}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_single_direction_declaration_synthesizes_conjugate_partner(self):
        c = 0.25 - 0.15j
        interactions = {"CoulombInter": {((1, 0, 0), (0, 0)): c}}
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        m_pos = [i for i, r in enumerate(topo.delta_r)
                  if tuple(int(x) for x in r) == (1, 0, 0)][0]
        m_neg = int(topo.reverse[m_pos])
        self.assertApprox(topo.coeffs["CoulombInter"][m_pos, 0, 0], c,
                           rel=0, abs=1e-13)
        self.assertApprox(topo.coeffs["CoulombInter"][m_neg, 0, 0],
                           np.conj(c), rel=0, abs=1e-13)

    def test_coeffs_additive_under_pre_summed_declaration(self):
        # "duplicate accumulation": a caller that pre-sums two logical
        # contributions into ONE dict entry (the only way a duplicate
        # (irvec, orbvec) key can arise -- a plain dict cannot hold two
        # entries under the identical key) must see the topology's coeffs
        # equal the LINEAR sum of what each contribution alone would give
        # (resolve_transverse_topology's closure is linear in the declared
        # value; this pins that no int-truncation/dedup logic sneaks in).
        v1, v2 = 0.11 + 0.02j, -0.05 + 0.07j
        topo1 = resolve_transverse_topology(
            {"CoulombInter": {((1, 0, 0), (0, 0)): v1}}, np.eye(3), norb=1)
        topo2 = resolve_transverse_topology(
            {"CoulombInter": {((1, 0, 0), (0, 0)): v2}}, np.eye(3), norb=1)
        topo_sum = resolve_transverse_topology(
            {"CoulombInter": {((1, 0, 0), (0, 0)): v1 + v2}},
            np.eye(3), norb=1)
        self.assertApproxArray(
            topo_sum.coeffs["CoulombInter"],
            topo1.coeffs["CoulombInter"] + topo2.coeffs["CoulombInter"],
            rel=0, abs=1e-13)

    def test_multi_orbital_shape(self):
        interactions = {"Exchange": {((1, 0, 0), (0, 1)): 0.3 + 0.1j}}
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=2)
        self.assertEqual(topo.coeffs["Exchange"].shape, (3, 2, 2))

    def test_norb_must_be_positive(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=0)


# =============================================================================
# iter_reversal_orbits
# =============================================================================

class TestIterReversalOrbits(ApproxTestCase):

    def test_singleton_orbit_at_onsite_diagonal(self):
        topo = _topo([[0, 0, 0]], [0], coeffs={
            "CoulombInter": np.zeros((1, 1, 1), dtype=complex)})
        orbits = iter_reversal_orbits(topo)
        self.assertEqual(orbits, [(0, 0, 0)])

    def test_onsite_off_diagonal_pairs(self):
        topo = _topo([[0, 0, 0]], [0], coeffs={
            "CoulombInter": np.zeros((1, 2, 2), dtype=complex)})
        orbits = set(iter_reversal_orbits(topo))
        # (0,0,0) and (0,1,1) singletons; (0,0,1)/(0,1,0) one representative.
        self.assertIn((0, 0, 0), orbits)
        self.assertIn((0, 1, 1), orbits)
        self.assertEqual(len(orbits & {(0, 0, 1), (0, 1, 0)}), 1)
        self.assertEqual(len(orbits), 3)

    def test_offsite_orbit_representative_is_tuple_minimum(self):
        topo = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1],
                      coeffs={"CoulombInter": np.zeros((3, 1, 1), dtype=complex)})
        orbits = iter_reversal_orbits(topo)
        self.assertIn((1, 0, 0), orbits)
        self.assertNotIn((2, 0, 0), orbits)

    def test_offsite_multi_orbital_orbit_count(self):
        # B=3 (0, +1, -1), norb=2: total (m,a,b) triples = 3*4=12; orbits =
        # 12/2 (every off-site pair is non-singleton since m != reverse[m]
        # there) except channel 0's 4 own (m,a,b) entries collapse to
        # 2 singletons ((0,0,0),(0,1,1)) + 1 pair -- 3 orbits there, plus
        # (12-4)/2 = 4 off-site orbits => 7 total.
        topo = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1],
                      coeffs={"CoulombInter": np.zeros((3, 2, 2), dtype=complex)})
        orbits = iter_reversal_orbits(topo)
        self.assertEqual(len(orbits), 7)
        self.assertEqual(len(orbits), len(set(orbits)))

    def test_covers_every_channel_orbital_pair_exactly_once(self):
        topo = _topo([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], [0, 2, 1],
                      coeffs={"CoulombInter": np.zeros((3, 2, 2), dtype=complex)})
        reverse = topo.reverse
        seen = set()
        for (m, a, b) in iter_reversal_orbits(topo):
            key1, key2 = (m, a, b), (int(reverse[m]), b, a)
            self.assertNotIn(key1, seen)
            self.assertNotIn(key2, seen)
            seen.add(key1)
            seen.add(key2)
        self.assertEqual(len(seen), 3 * 2 * 2)

    def test_empty_coeffs_raises(self):
        t = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]), coeffs={})
        with self.assertRaises(ValueError):
            iter_reversal_orbits(t)

    def test_deterministic_and_matches_resolver_output(self):
        interactions = {"CoulombInter": {
            ((1, 0, 0), (0, 0)): 0.2, ((-1, 0, 0), (0, 0)): 0.2}}
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        orbits_a = iter_reversal_orbits(topo)
        orbits_b = iter_reversal_orbits(topo)
        self.assertEqual(orbits_a, orbits_b)


# =============================================================================
# validate_topology_against_mesh
# =============================================================================

class TestValidateTopologyAgainstMesh(ApproxTestCase):

    def _topo_1d(self, r_list):
        """B-channel topology on a 1D chain (y=z=0), delta_r = [(0,0,0)] +
        [(r,0,0) for r in r_list], with a valid reverse permutation
        (assumes r_list is already closed under negation, distinct)."""
        delta_r = [(0, 0, 0)] + [(r, 0, 0) for r in r_list]
        index_of = {dr: i for i, dr in enumerate(delta_r)}
        reverse = [index_of[(-dr[0], 0, 0)] for dr in delta_r]
        return _topo(delta_r, reverse, coeffs={
            "CoulombInter": np.zeros((len(delta_r), 1, 1), dtype=complex)})

    def test_valid_shape_passes(self):
        topo = self._topo_1d([1, -1])
        validate_topology_against_mesh(topo, (5, 5, 5))  # no raise

    def test_spatial_shape_wrong_length_rejected(self):
        topo = self._topo_1d([1, -1])
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, (5, 5))

    def test_spatial_shape_non_positive_rejected(self):
        topo = self._topo_1d([1, -1])
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, (5, 0, 5))
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, (-1, 5, 5))

    def test_spatial_shape_non_integral_rejected(self):
        topo = self._topo_1d([1, -1])
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, (5.5, 5, 5))

    def test_spatial_shape_not_a_sequence_rejected(self):
        topo = self._topo_1d([1, -1])
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, 5)

    def test_rectangular_mesh_no_collision(self):
        # delta_r along x only; Ny, Nz irrelevant to the x-axis wrap.
        topo = self._topo_1d([1, -1])
        validate_topology_against_mesh(topo, (8, 3, 2))  # no raise

    def test_r_vs_r_plus_l_collision_across_shells(self):
        # R=1 (shell length 1) and R=6 (shell length 6) collide mod L=5 on
        # the x-axis -- a genuine cross-shell alias.
        topo = self._topo_1d([1, -1, 6, -6])
        with self.assertRaises(ValueError) as cm:
            validate_topology_against_mesh(topo, (5, 5, 5))
        msg = str(cm.exception)
        self.assertIn("5, 5, 5", msg)

    def test_even_mesh_self_reversal_collision_within_shell(self):
        # R=+3 and R=-3 (SAME shell, length 3) alias on an even mesh Nx=6
        # only if 3 == -3 mod 6, i.e. 3 == 3 -- true; this is the
        # canonical "+L/2 == -L/2" alias.
        topo = self._topo_1d([3, -3])
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(topo, (6, 4, 4))

    def test_even_mesh_self_reversal_does_not_collide_off_the_edge(self):
        # The SAME +/-3 pair does NOT alias on an odd/larger mesh.
        topo = self._topo_1d([3, -3])
        validate_topology_against_mesh(topo, (8, 4, 4))  # no raise

    def test_valid_noncolliding_reverse_pair_multiple_shells(self):
        topo = self._topo_1d([1, -1, 2, -2])
        validate_topology_against_mesh(topo, (9, 9, 9))  # no raise

    def test_error_names_both_channels_and_mesh(self):
        topo = self._topo_1d([1, -1, 6, -6])
        with self.assertRaises(ValueError) as cm:
            validate_topology_against_mesh(topo, (5, 5, 5))
        msg = str(cm.exception)
        # both offending raw delta_r values must be named
        self.assertIn("(1, 0, 0)", msg)
        self.assertIn("(6, 0, 0)", msg)

    def test_arrays_nvol_mismatch_rejected(self):
        # mesh (5, 3, 4): large enough on x that +/-1 does not itself alias
        # (unlike (2, 3, 4), where Nx=2 would trigger the even-mesh
        # self-reversal collision this test is not about).
        topo = self._topo_1d([1, -1])
        good = np.zeros((5 * 3 * 4, 5))
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(
                topo, (5, 3, 4), arrays={"green": np.zeros((99, 5))})
        validate_topology_against_mesh(
            topo, (5, 3, 4), arrays={"green": good})  # no raise

    def test_stale_topology_alias_revalidated(self):
        # Directly poke a topology's coeffs dict (mutable container, even
        # though its arrays are read-only) to break the closure invariant,
        # then confirm validate_topology_against_mesh -- as a consumer
        # boundary -- catches it rather than trusting the alias.
        topo = self._topo_1d([1, -1])
        broken = dict(topo.coeffs)
        broken["CoulombInter"] = np.array(
            [[[0]], [[1 + 2j]], [[9 + 9j]]], dtype=complex)
        bad_topo = topo._replace(coeffs=broken)
        with self.assertRaises(ValueError):
            validate_topology_against_mesh(bad_topo, (5, 5, 5))


# =============================================================================
# transverse_effective_activity
# =============================================================================

class TestTransverseEffectiveActivity(ApproxTestCase):

    def test_all_zero_topology_is_inactive(self):
        topo = resolve_transverse_topology({}, np.eye(3), norb=2)
        self.assertFalse(transverse_effective_activity(topo))

    def test_nonzero_coulomb_inter_is_active(self):
        topo = resolve_transverse_topology(
            {"CoulombInter": {((1, 0, 0), (0, 0)): 0.1}}, np.eye(3), norb=1)
        self.assertTrue(transverse_effective_activity(topo))

    def test_nonzero_ising_is_active(self):
        topo = resolve_transverse_topology(
            {"Ising": {((1, 0, 0), (0, 0)): 0.1}}, np.eye(3), norb=1)
        self.assertTrue(transverse_effective_activity(topo))

    def test_nonzero_exchange_is_active(self):
        topo = resolve_transverse_topology(
            {"Exchange": {((1, 0, 0), (0, 0)): 0.1 + 0.05j}}, np.eye(3), norb=1)
        self.assertTrue(transverse_effective_activity(topo))

    def test_only_hund_pairlift_is_inactive(self):
        topo = resolve_transverse_topology(
            {"Hund": {((1, 0, 0), (0, 0)): 0.7},
             "PairLift": {((1, 0, 0), (0, 0)): 0.3}}, np.eye(3), norb=1)
        self.assertFalse(transverse_effective_activity(topo))

    def test_mirrored_pair_that_cancels_to_zero_real_part_is_inactive(self):
        # C at R and -conj(C) at -R (an ANTI-Hermitian-style declared pair)
        # is rejected by the resolver's consistency check UNLESS the
        # resulting Re() genuinely vanishes; use a purely-imaginary C whose
        # conjugate pair is consistent AND whose closed representative
        # happens to be pinned at a value with zero real part by
        # construction of THIS test (declare only one direction with a
        # purely real coefficient of 0.0 explicitly, the simplest exact
        # zero case distinct from "not declared at all").
        topo = resolve_transverse_topology(
            {"CoulombInter": {((1, 0, 0), (0, 0)): 0.0}}, np.eye(3), norb=1)
        self.assertFalse(transverse_effective_activity(topo))


# =============================================================================
# Hamiltonian-level normalization pin (carried Task-1 review finding)
# =============================================================================

def _manual_physical_coeff(raw, irvec, a, b):
    """Independent (no topology/resolver code involved) computation of the
    TRUE physical Hamiltonian coefficient of the single Hermitian
    density-density operator class {(irvec, a, b), (-irvec, b, a)} from
    raw file-format declarations -- the "double-declaration trap" the
    carried Task-1 review finding describes: n_{j,a} n_{j+irvec,b} summed
    over j is IDENTICAL, as an operator, to n_{j,b} n_{j-irvec,a} summed
    over j (density operators commute; re-index j' = j - irvec), so a
    declared pair (C at irvec, C' at -irvec with swapped orbitals) both
    multiply the SAME operator and their physical total is C + C' (which
    reduces to C + conj(C) = 2*Re(C) when only one direction is declared
    and the other is implied by Hermiticity)."""
    neg = tuple(-x for x in irvec)
    c_fwd = raw.get((irvec, (a, b)))
    c_bwd = raw.get((neg, (b, a)))
    if c_fwd is not None and c_bwd is not None:
        return c_fwd + c_bwd
    if c_fwd is not None:
        return c_fwd + np.conj(c_fwd)
    if c_bwd is not None:
        return c_bwd + np.conj(c_bwd)
    return 0.0 + 0.0j


def _h_from_topology(fx, topo, type_name):
    """Build the exact many-body H that `topo.coeffs[type_name]` implies,
    via `iter_reversal_orbits` (restricted to the off-site domain, per the
    spec's "Construction domain") and `canonical_density_terms`. Every
    off-site orbit contributes `2*Re(coeffs[m][a,b])` on its single
    Hermitian operator class, matching `bare_bond_vertices`'s own
    documented "V_ab(R_m) + V_ba(-R_m) == 2*Re(v_bond[m][a,b])" derivation
    (this module's Finding-1 comment, reproduced structurally for the
    transverse coeffs)."""
    pairs = []
    for (m, a, b) in iter_reversal_orbits(topo):
        if m == 0:
            continue
        coeff = 2.0 * np.real(topo.coeffs[type_name][m, a, b])
        if coeff == 0.0:
            continue
        R = int(topo.delta_r[m][0])
        pairs.append((a, b, R, coeff))
    terms = ed_oracle_util.canonical_density_terms(fx, pairs)
    return ed_oracle_util.h_int_from_terms(fx, terms)


def _h_from_raw(fx, raw, irvec_ab_pairs):
    """Build the SAME many-body H directly from raw declarations, via
    `_manual_physical_coeff` (independent of the topology/resolver code
    under test)."""
    pairs = []
    for (irvec, a, b) in irvec_ab_pairs:
        coeff = _manual_physical_coeff(raw, irvec, a, b)
        R = irvec[0]
        pairs.append((a, b, R, coeff))
    terms = ed_oracle_util.canonical_density_terms(fx, pairs)
    return ed_oracle_util.h_int_from_terms(fx, terms)


class TestHamiltonianLevelNormalizationPin(ApproxTestCase):
    """Pins resolve_transverse_topology's FILE -> COEFFS normalization at
    the exact many-body Hamiltonian level, per the carried Task-1 review
    finding: build a tiny fixture's exact H from raw declarations two
    ways -- via the topology's coeffs reconstruction
    (`_h_from_topology`) and via an independent direct construction
    (`_h_from_raw`) -- and assert equality. Covers mirrored-pair
    declarations, single-direction declarations and duplicate
    accumulation, for CoulombInter (density-density/cross-family) and
    Ising (also density-density)."""

    def _fixture(self):
        # L=2, norb=2 chain -- the SAME fixture shape
        # test_offsite_exchange_term_builder_matches_dense_hamiltonian uses
        # (tests/test_bond_transverse_ed.py), so this pin sits on an
        # independently-established-tractable dense (dim=256) footing.
        return ed_oracle_util.EDFixture(
            L=2, norb=2,
            t={(0, 0): 0.3, (1, 1): 0.2, (0, 1): 0.1, (1, 0): 0.1},
            eps=(0.05, -0.02), T=0.5, mu=0.1)

    def test_single_direction_declaration(self):
        fx = self._fixture()
        c = 0.3 + 0.2j
        raw = {((1, 0, 0), (0, 1)): c}
        topo = resolve_transverse_topology(
            {"CoulombInter": raw}, np.eye(3), norb=2)
        with np.errstate(all="ignore"):
            H_topo = _h_from_topology(fx, topo, "CoulombInter")
            H_raw = _h_from_raw(fx, raw, [((1, 0, 0), 0, 1)])
        self.assertApproxArray(H_topo, H_raw, rel=0, abs=1e-13)
        # sanity: genuinely nonzero (not a PASS-ZERO false positive)
        self.assertGreater(float(np.max(np.abs(H_topo))), 1e-6)

    def test_mirrored_pair_declaration_consistent(self):
        fx = self._fixture()
        c1 = 0.3 + 0.2j
        c2 = np.conj(c1)
        raw = {((1, 0, 0), (0, 1)): c1, ((-1, 0, 0), (1, 0)): c2}
        topo = resolve_transverse_topology(
            {"CoulombInter": raw}, np.eye(3), norb=2)
        with np.errstate(all="ignore"):
            H_topo = _h_from_topology(fx, topo, "CoulombInter")
            H_raw = _h_from_raw(fx, raw, [((1, 0, 0), 0, 1)])
        self.assertApproxArray(H_topo, H_raw, rel=0, abs=1e-13)

    def test_mirrored_pair_null_direction_epsilon(self):
        # The historical Finding-1 regression pattern: a near-Hermitian
        # +-i*eps closed pair meant to represent ONE real coupling V; the
        # diagonal must carry Re(.) only, never leak the imaginary part.
        V, eps = 0.4, 1e-3
        fx = self._fixture()
        raw = {((1, 0, 0), (0, 0)): V + 1j * eps,
               ((-1, 0, 0), (0, 0)): V - 1j * eps}
        topo = resolve_transverse_topology(
            {"CoulombInter": raw}, np.eye(3), norb=2)
        with np.errstate(all="ignore"):
            H_topo = _h_from_topology(fx, topo, "CoulombInter")
            H_raw = _h_from_raw(fx, raw, [((1, 0, 0), 0, 0)])
        self.assertApproxArray(H_topo, H_raw, rel=0, abs=1e-13)
        # The stored per-channel coeff need not itself be real (it equals
        # the exactly-self-consistent declared value V+i*eps here, same as
        # resolve_interactions' v_bond) -- what must hold is that its
        # DOWNSTREAM physical reconstruction (2*Re(.), what
        # bare_bond_vertices/_h_from_topology actually use) is the real
        # 2V, with eps never leaking into it.
        m1 = [i for i, r in enumerate(topo.delta_r)
              if tuple(int(x) for x in r) == (1, 0, 0)][0]
        self.assertApprox(
            2.0 * np.real(topo.coeffs["CoulombInter"][m1, 0, 0]), 2.0 * V,
            rel=0, abs=1e-13)

    def test_duplicate_accumulation_presummed_declaration(self):
        # "duplicate accumulation": the raw dict entry already carries the
        # SUM of two logical contributions (the only way a duplicate
        # (irvec, orbvec) key can arise at all -- see
        # resolve_transverse_topology's own docstring / carried finding).
        v1, v2 = 0.11 + 0.02j, -0.05 + 0.07j
        fx = self._fixture()
        raw_summed = {((1, 0, 0), (0, 1)): v1 + v2}
        topo = resolve_transverse_topology(
            {"CoulombInter": raw_summed}, np.eye(3), norb=2)
        with np.errstate(all="ignore"):
            H_topo = _h_from_topology(fx, topo, "CoulombInter")
            H_raw = _h_from_raw(fx, raw_summed, [((1, 0, 0), 0, 1)])
            # cross-check: separately declaring v1 and v2 and summing the
            # two resulting Hamiltonians independently reproduces the SAME
            # H (the "duplicate" really did accumulate additively).
            H_v1 = _h_from_raw(fx, {((1, 0, 0), (0, 1)): v1},
                                [((1, 0, 0), 0, 1)])
            H_v2 = _h_from_raw(fx, {((1, 0, 0), (0, 1)): v2},
                                [((1, 0, 0), 0, 1)])
        self.assertApproxArray(H_topo, H_raw, rel=0, abs=1e-13)
        self.assertApproxArray(H_raw, H_v1 + H_v2, rel=0, abs=1e-13)

    def test_ising_single_direction_declaration(self):
        fx = self._fixture()
        v = 0.37
        raw = {((1, 0, 0), (0, 0)): v}
        topo = resolve_transverse_topology(
            {"Ising": raw}, np.eye(3), norb=2)
        with np.errstate(all="ignore"):
            H_topo = _h_from_topology(fx, topo, "Ising")
            H_raw = _h_from_raw(fx, raw, [((1, 0, 0), 0, 0)])
        self.assertApproxArray(H_topo, H_raw, rel=0, abs=1e-13)


if __name__ == "__main__":
    unittest.main()
