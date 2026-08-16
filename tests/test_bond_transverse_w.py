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

Also (Phase W, Task 3, appended below the Task-2 topology classes):
production tests for ``hwave.solver.bond_channels.W_pm_bond``, the
bond-resolved transverse vertex built from ``TransverseTopology`` per the
spec's "The vertex -- element equations" (steps 1-3, the AMENDED
2026-08-16 flip-family assignment). The ordered-record equations
themselves were already numerically validated against ED by Gate W0 in
``tests/test_bond_transverse_ed.py`` (Task 1) -- NOT re-adjudicated here;
this module instead CROSS-PINS the production ``W_pm_bond`` against Gate
W0's test-local oracle reference, ``w_expected_from_records`` (imported
module-qualified from that module), on the same three W0 ED granule
fixtures, plus structural (Hermiticity, per-entry, B=1 reduction, q-phase
convention) and validation tests.

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

import logging
import os
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
from hwave import sc as _sc
from hwave.solver import bubble as bubble_mod
from hwave.solver import rpa as rpa_mod
from hwave.solver.bond_channels import (
    TransverseTopology,
    resolve_transverse_topology,
    resolve_interactions,
    iter_reversal_orbits,
    validate_topology_against_mesh,
    transverse_effective_activity,
    W_pm_bond,
)
from tests import ed_oracle_util
from tests import test_bond_transverse_ed as _ted
from tests.approx_util import ApproxTestCase, assert_approx_array


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

    # -----------------------------------------------------------------
    # Strict input validation (review fix): norb, irvec components and
    # orbital indices used to pass silently through int(), so norb=1.5
    # became 1 and a fractional irvec became the wrong channel; an
    # out-of-range orbital index used to leak a bare IndexError instead
    # of an actionable ValueError.
    # -----------------------------------------------------------------

    def test_norb_non_integral_rejected(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=1.5)

    def test_norb_bool_rejected(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=True)

    def test_norb_nan_rejected(self):
        with self.assertRaises(ValueError):
            resolve_transverse_topology({}, np.eye(3), norb=float("nan"))

    def test_irvec_non_integral_component_rejected(self):
        interactions = {"CoulombInter": {((1.5, 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_irvec_wrong_length_rejected(self):
        interactions = {"CoulombInter": {((1, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_orbvec_non_integral_component_rejected(self):
        interactions = {"CoulombInter": {((1, 0, 0), (0.5, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_orbvec_wrong_length_rejected(self):
        interactions = {"CoulombInter": {((1, 0, 0), (0, 0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_orbital_index_out_of_range_rejected(self):
        # norb=1 -> only orbital index 0 is valid; index 1 used to leak
        # a bare IndexError from the (norb, norb) array assignment
        # instead of an actionable ValueError.
        interactions = {"CoulombInter": {((1, 0, 0), (1, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_orbital_index_negative_rejected(self):
        interactions = {"CoulombInter": {((1, 0, 0), (-1, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_numpy_bool_component_rejected(self):
        # np.bool_ is not a Python bool subclass; the guard must reject
        # it explicitly rather than let it pass as the integer 1.
        interactions = {"CoulombInter": {((np.bool_(True), 0, 0),
                                          (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_numeric_string_component_rejected(self):
        # "1" must not coerce (no float()/int() round-trip on strings).
        interactions = {"CoulombInter": {(("1", 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_non_scalar_component_rejected(self):
        # A hashable non-scalar (a nested tuple) must be rejected with
        # the guard's ValueError, not fall through to int() semantics.
        interactions = {"CoulombInter": {(((1, 2), 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    def test_large_exact_int_preserved_not_float_rounded(self):
        # 2**53 + 1 is not representable in float64; an int()-via-float
        # conversion would silently round it to 2**53. The exact-int
        # path must preserve it bit for bit (it then fails no validation
        # here: it is a legitimate -- if absurd -- displacement).
        big = 2**53 + 1
        interactions = {"CoulombInter": {((big, 0, 0), (0, 0)): 0.2},
                        }
        topo = resolve_transverse_topology(interactions, np.eye(3), norb=1)
        self.assertIn(big, set(int(r[0]) for r in topo.delta_r))

    def test_huge_float_component_rejected(self):
        # A float above the 2**53 contiguous-integer bound is ambiguous
        # as an integer and must be rejected, not silently converted.
        interactions = {"CoulombInter": {((2.0**60, 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    @unittest.skipIf(np.finfo(np.longdouble).nmant
                     <= np.finfo(np.float64).nmant,
                     "longdouble is float64 on this platform")
    def test_longdouble_fractional_rejected_at_native_precision(self):
        # A fractional longdouble whose value would round to an INTEGRAL
        # float64 must still be rejected: the check runs at native
        # precision, never through a float() narrowing.
        frac = np.longdouble(1) + np.finfo(np.longdouble).eps
        self.assertEqual(float(frac), 1.0)  # narrows to integral float64
        interactions = {"CoulombInter": {((frac, 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)

    @unittest.skipIf(np.finfo(np.longdouble).nmant
                     <= np.finfo(np.float64).nmant,
                     "longdouble is float64 on this platform")
    def test_longdouble_just_above_2p53_rejected(self):
        # 2**53 + 1 is exactly representable in extended precision but
        # narrows onto the accepted 2**53 boundary in float64 -- it must
        # be rejected at native precision.
        val = np.longdouble(2) ** 53 + np.longdouble(1)
        self.assertNotEqual(val, np.longdouble(2) ** 53)
        interactions = {"CoulombInter": {((val, 0, 0), (0, 0)): 0.2}}
        with self.assertRaises(ValueError):
            resolve_transverse_topology(interactions, np.eye(3), norb=1)


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


# =============================================================================
# W_pm_bond (Phase W, Task 3)
# =============================================================================
#
# Production tests for hwave.solver.bond_channels.W_pm_bond -- the
# bond-resolved transverse vertex built from a TransverseTopology per the
# spec's "The vertex -- element equations" (steps 1-3, the AMENDED
# 2026-08-16 flip-family assignment). The ordered-record equations
# themselves were already numerically validated against ED by Gate W0
# (tests/test_bond_transverse_ed.py, Task 1, imported module-qualified
# below as ``_ted``); this section cross-pins the PRODUCTION function
# against that gate's oracle, plus structural/validation coverage.

def _closed_coeffs(B, norb, reverse, declared):
    """Build a Hermitian-closed ``(B, norb, norb)`` complex coeffs array
    from a SINGLE-DIRECTION declaration dict ``{m: (norb, norb) array}``
    -- exactly the shape gate W0's granule fixtures declare (tests/
    test_bond_transverse_ed.py's ``TestW0Granules``/
    ``w_expected_from_records``): every declared channel's mirror partner
    ``reverse[m]`` is synthesized as ``conj(block).T``, the same
    single-direction synthesis ``resolve_transverse_topology`` performs
    (``_close_offsite_hermitian``) -- reproduced directly here (not via
    that resolver) so a ``TransverseTopology`` can be hand-built on the
    EXACT ``_TopoLike`` channel layout gate W0 used, with ``coeffs``
    satisfying the closure invariant :class:`TransverseTopology`'s own
    constructor enforces. Missing channels stay exactly zero (the
    ``TransverseTopology`` on-site-channel-0-is-zero invariant is
    automatic here since gate W0's declarations dicts never key on
    ``m == 0`` -- "coeffs carry OFF-SITE channels only", spec)."""
    arr = np.zeros((B, norb, norb), dtype=complex)
    for m, block in declared.items():
        block = np.asarray(block, dtype=complex)
        arr[m] = block
        mr = int(reverse[m])
        if mr != m:
            arr[mr] = np.conj(block).T
    return arr


def _topo_from_topolike(topo_like, norb, declarations):
    """Build a production :class:`TransverseTopology` on the EXACT channel
    layout of a gate-W0 ``_TopoLike`` (tests/test_bond_transverse_ed.py),
    with ``coeffs`` closed from the SAME single-direction declarations
    dict a W0 granule used (``declarations["CoulombInter"]``/``["Ising"]``/
    ``["Exchange"]``, each ``{m: (norb, norb) array}``; a missing type key
    is treated as "declares nothing off-site", matching
    ``w_expected_from_records``'s own ``declarations.get(kind, {})``
    contract)."""
    delta_r = np.asarray(topo_like.delta_r)
    reverse = np.asarray(topo_like.reverse)
    B = delta_r.shape[0]
    coeffs = {}
    for type_name in ("CoulombInter", "Ising", "Exchange"):
        declared = declarations.get(type_name, {})
        coeffs[type_name] = _closed_coeffs(B, norb, reverse, declared)
    return TransverseTopology(delta_r=delta_r, reverse=reverse, coeffs=coeffs)


class TestWPmBondCrossPinAgainstW0Oracle(ApproxTestCase):
    """THE TASK'S CENTERPIECE: production ``W_pm_bond`` matches Gate W0's
    ED-validated test-local oracle, ``w_expected_from_records`` (tests/
    test_bond_transverse_ed.py, Task 1, imported module-qualified as
    ``_ted``), at ``rel=0, abs=1e-13``, on ALL THREE W0 ED granule
    fixtures -- (a) multi-orbital off-site Exchange (complex J,
    non-self-inverse q, both orientations), (b) complex off-site
    CoulombInter (g1/g2 shells), (c) off-site Ising (g1/g2 shells).

    The topology layout (``delta_r``/``reverse``), declared coefficients
    and ``q_mesh`` are copied VERBATIM from ``_ted.TestW0Granules``'s
    three granule bodies -- the pure numeric inputs those granules PASSED
    Gate W0's ED adjudication with (the ED machinery itself -- ``fx``,
    Richardson, ``SectorED`` -- is NOT re-run here; this class exercises
    only the production ``W_pm_bond`` vs. the oracle formula, both fed the
    identical topology/coefficients/mesh)."""

    def _run_cross_pin(self, topo_like, norb, declarations, q_mesh):
        topo = _topo_from_topolike(topo_like, norb, declarations)
        onsite = np.asarray(declarations["onsite_ham_pm"], dtype=complex)
        nvol = q_mesh[0] * q_mesh[1] * q_mesh[2]
        onsite_nvol = np.broadcast_to(
            onsite[None, ...], (nvol,) + onsite.shape).copy()
        got = W_pm_bond(topo, onsite_nvol, spatial_shape=q_mesh)
        expected = _ted.w_expected_from_records(
            topo_like, declarations, q_mesh)
        assert_approx_array(got, expected, rel=0, abs=1e-13)

    def test_granule_a_multiorbital_offsite_exchange(self):
        # Verbatim data from
        # _ted.TestW0Granules.test_granule_a_multiorbital_offsite_exchange:
        # L=3, norb=2, off-site Exchange at R=1 with a=0 != b=1 and a
        # genuinely complex J; two directions ("real"/"imag" phase).
        norb = 2
        topo_like = _ted._TopoLike(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]),
            reverse=np.array([0, 2, 1]))
        q_mesh = (3, 1, 1)
        onsite_ham_pm = np.zeros((2, 2, 2, 2), dtype=complex)
        a, b = 0, 1
        for phase in (1.0 + 0.0j, 0.0 + 1.0j):
            block = np.zeros((2, 2), dtype=complex)
            block[a, b] = phase
            declarations = {"onsite_ham_pm": onsite_ham_pm,
                             "Exchange": {1: block}}
            self._run_cross_pin(topo_like, norb, declarations, q_mesh)

    def test_granule_b_complex_offsite_coulomb_inter(self):
        # Verbatim data from
        # _ted.TestW0Granules.test_granule_b_complex_offsite_coulomb_inter:
        # L=5, norb=1, off-site CoulombInter, g1 (R=+1) and g2 (R=+2).
        norb = 1
        topo_like = _ted._TopoLike(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                               [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1, 4, 3]))
        q_mesh = (5, 1, 1)
        onsite_ham_pm = np.zeros((1, 1, 1, 1), dtype=complex)
        C = 0.6 + 0.0j
        for m_pos in (1, 3):
            declarations = {
                "onsite_ham_pm": onsite_ham_pm,
                "CoulombInter": {m_pos: np.array([[C]], dtype=complex)},
            }
            self._run_cross_pin(topo_like, norb, declarations, q_mesh)

    def test_granule_c_offsite_ising(self):
        # Verbatim data from _ted.TestW0Granules.test_granule_c_offsite_ising:
        # L=5, norb=1, off-site Ising, g1 (R=+1) and g2 (R=+2).
        norb = 1
        topo_like = _ted._TopoLike(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                               [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1, 4, 3]))
        q_mesh = (5, 1, 1)
        onsite_ham_pm = np.zeros((1, 1, 1, 1), dtype=complex)
        for m_pos in (1, 3):
            declarations = {
                "onsite_ham_pm": onsite_ham_pm,
                "Ising": {m_pos: np.array([[1.0 + 0j]], dtype=complex)},
            }
            self._run_cross_pin(topo_like, norb, declarations, q_mesh)


class TestWPmBondStructural(ApproxTestCase):
    """Structural pins independent of the ED cross-pin above: Hermiticity
    at every q (including a Hermitian-closed COMPLEX CoulombInter +-i*eps
    pair -- the W0 blind spot the Task-1 review flagged: a bug that stored
    the raw complex coefficient instead of ``Re(.)`` on the bond-diagonal
    would put a genuinely complex value on a MATRIX DIAGONAL, breaking
    Hermiticity), the two Exchange records pinned individually at a
    non-self-inverse q with a genuinely complex J, and the B=1 on-site
    reduction (spec step 1's algebraic basis for gate W1)."""

    def _offsite_topo(self, norb, coeffs):
        delta_r = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
        reverse = np.array([0, 2, 1])
        full = {}
        for type_name in ("CoulombInter", "Ising", "Exchange"):
            full[type_name] = coeffs.get(
                type_name, np.zeros((3, norb, norb), dtype=complex))
        return TransverseTopology(
            delta_r=delta_r, reverse=reverse, coeffs=full)

    def test_hermitian_at_every_q_including_complex_coulomb_inter(self):
        norb = 1
        V, eps = 0.4, 1.0e-3
        ci = np.zeros((3, 1, 1), dtype=complex)
        ci[1, 0, 0] = V + 1j * eps    # channel m=1, R=(1,0,0), declared
        ci[2, 0, 0] = V - 1j * eps    # channel m=2 = reverse[1], R=(-1,0,0)
        ex = np.zeros((3, 1, 1), dtype=complex)
        ex[1, 0, 0] = 0.5 - 0.2j
        ex[2, 0, 0] = np.conj(ex[1, 0, 0])
        ising = np.zeros((3, 1, 1), dtype=complex)
        ising[1, 0, 0] = 0.2
        ising[2, 0, 0] = 0.2
        topo = self._offsite_topo(
            norb, {"CoulombInter": ci, "Exchange": ex, "Ising": ising})
        Nx = 5
        onsite = np.zeros((Nx, 1, 1, 1, 1), dtype=complex)
        W = W_pm_bond(topo, onsite, spatial_shape=(Nx, 1, 1))
        got_dagger = np.conj(np.transpose(W, (0, 2, 1)))
        self.assertApproxArray(got_dagger, W, rel=0, abs=1e-12)

    def test_exchange_entries_pinned_individually(self):
        norb = 2
        J = 0.6 - 0.35j
        ex = np.zeros((3, 2, 2), dtype=complex)
        ex[1, 0, 1] = J
        ex[2, 1, 0] = np.conj(J)
        topo = self._offsite_topo(norb, {"Exchange": ex})
        Nx = 5
        onsite = np.zeros((Nx, 2, 2, 2, 2), dtype=complex)
        W = W_pm_bond(topo, onsite, spatial_shape=(Nx, 1, 1))
        kx = 2.0 * np.pi * np.arange(Nx) / Nx
        q_idx = 1                      # q = 2*pi/5, non-self-inverse (L odd)
        q = kx[q_idx]
        phase = np.exp(1j * q * 1.0)   # R = delta_r[1] = (1, 0, 0)
        idx_aa, idx_bb = 0, 3           # channel 0, (0,0) and (1,1)
        self.assertApproxArray(
            np.array([W[q_idx, idx_aa, idx_bb]]),
            np.array([-np.conj(J) * np.conj(phase)]), rel=0, abs=1e-13)
        self.assertApproxArray(
            np.array([W[q_idx, idx_bb, idx_aa]]),
            np.array([-J * phase]), rel=0, abs=1e-13)

    def test_b1_reduction_equals_broadcast_onsite_ham_pm(self):
        norb = 2
        nd = norb * norb
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]),
            coeffs={"CoulombInter": np.zeros((1, norb, norb), dtype=complex)})
        onsite_2d = _ted._h2_ham_pm_expected("CoulombIntra", 1.3, norb=norb)
        onsite_4d = onsite_2d.reshape(norb, norb, norb, norb)
        nvol = 6
        onsite_nvol = np.broadcast_to(
            onsite_4d[None, :, :, :, :], (nvol,) + onsite_4d.shape).copy()
        W = W_pm_bond(topo, onsite_nvol, spatial_shape=(3, 2, 1))
        expected = np.broadcast_to(onsite_2d[None, :, :], (nvol, nd, nd))
        self.assertApproxArray(W, expected, rel=0, abs=1e-15)


class TestWPmBondQPhaseCrossPin(ApproxTestCase):
    """Pins ``W_pm_bond``'s q-phase convention against
    ``sc._build_bond_m0_blocks``'s own phase composition (spec: "the SAME
    convention sc._build_bond_m0_blocks' phase uses") on a SHARED fixture:
    a single-direction off-site declaration at R=(1,0,0), (a,b)=(0,1),
    norb=2. Fed to ``_build_bond_m0_blocks`` as a real CoulombInter
    coefficient (whose Hartree ``(00,11)`` matrix element is driven by
    that builder's OWN ``exp(-i q.R)`` phase, extracted directly from its
    real numeric output -- not assumed) and to ``W_pm_bond`` as an
    Exchange coefficient at the SAME R (whose flip-family record carries
    ``exp(-i q.R)`` on the ``(0,0,0),(0,1,1)`` element, spec step 3). Both
    consumers see the IDENTICAL ``kx_array``/``ky_array``/``kz_array``
    mesh (``np.linspace(0, 2*pi, N, endpoint=False)`` -- the same
    construction ``hwave/sc.py``'s own q-grid setup uses ahead of calling
    ``_build_bond_m0_blocks``, e.g. around its ``calc_eliashberg`` driver),
    so this comparison genuinely exercises the SAME q.R value on both
    sides via a REAL call to the shared builder, not two independently
    assumed formulas."""

    def test_qr_matches_build_bond_m0_blocks_phase(self):
        norb = 2
        Nx, Ny, Nz = 5, 1, 1
        kx_array = np.linspace(0.0, 2.0 * np.pi, Nx, endpoint=False)
        ky_array = np.linspace(0.0, 2.0 * np.pi, Ny, endpoint=False)
        kz_array = np.linspace(0.0, 2.0 * np.pi, Nz, endpoint=False)

        C = 0.6 + 0.0j
        bond_set = resolve_interactions(
            {((1, 0, 0), (0, 1)): C}, np.eye(3), norb)
        _S0, C0 = _sc._build_bond_m0_blocks(
            bond_set, {}, {}, norb, kx_array, ky_array, kz_array)
        # Hartree (aa,bb) = (00,11) slot: C0[..., 0, 3] == 2*V_q[0, 1], and
        # V_q[0, 1] is driven ONLY by the declared channel m=R=(1,0,0)'s
        # own [0, 1] matrix entry -- the synthesized reverse partner only
        # populates [1, 0], never [0, 1] -- so C0[q_idx,0,0,0,3] ==
        # 2*C*exp(-i q.R) exactly (C real).
        q_idx = 1
        phase_from_builder = C0[q_idx, 0, 0, 0, 3] / (2.0 * C)

        J = 0.35 - 0.2j
        reverse = np.array([0, 2, 1])
        block = np.zeros((norb, norb), dtype=complex)
        block[0, 1] = J
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]),
            reverse=reverse,
            coeffs={"Exchange": _closed_coeffs(3, norb, reverse, {1: block})})
        onsite = np.zeros((Nx * Ny * Nz, norb, norb, norb, norb),
                           dtype=complex)
        W = W_pm_bond(topo, onsite, spatial_shape=(Nx, Ny, Nz))
        idx_00, idx_11 = 0, 3
        got = W[q_idx, idx_00, idx_11]
        expected = -np.conj(J) * phase_from_builder
        self.assertApprox(got, expected, rel=0, abs=1e-12)


class TestWPmBondChannelZeroEmbeddingPin(ApproxTestCase):
    """Fix round 1 review finding (Task 3): the channel-0 embedding
    (``W[:, 0:nd, 0:nd] = ham_pm_onsite.reshape(nvol, nd, nd)``) had NO
    discriminating numeric pin -- every existing cross-pin fixture uses an
    all-zero on-site block, and ``TestWPmBondStructural``'s B=1 reduction
    test's only nonzero entry sits at the fully-degenerate ``(0,0,0,0)``
    index. At that degenerate index the CORRECT ``(a,b)/(c,d)``
    pair-flattening (row = ``a*norb+b``, col = ``c*norb+d``, exactly what
    a bare ``.reshape(nvol, nd, nd)`` does to a contiguous
    ``(nvol,a,b,c,d)`` tensor -- matching ``_solve_rpa``'s own
    ``ham.reshape(nvol, ndx, ndx)`` consumer convention) is numerically
    IDENTICAL to a WRONG ``(a,c)/(b,d)`` grouping (row = ``a*norb+c``, col
    = ``b*norb+d`` -- what ``onsite.transpose(0,1,3,2,4)`` before reshape
    would produce), so a reviewer-supplied mutant passed every committed
    test.

    This test uses a NON-DIAGONAL-DEGENERATE Hermitian on-site tensor
    (nonzero ONLY at ``(a,b,c,d) = (0,0,1,1)`` and its Hermitian partner
    ``(1,1,0,0)``) where the two groupings land the value at DIFFERENT
    cells (``(0,3)`` vs. ``(1,1)`` -- verified by direct calculation, see
    the module-level mutation-check note below) and compares
    ``W_pm_bond``'s channel-0 block ELEMENTWISE against an expected matrix
    built BY HAND via a plain per-index loop (never via ``.reshape`` --
    validating production against itself would prove nothing)."""

    def test_channel_zero_embedding_pairs_ab_cd_not_ac_bd(self):
        norb = 2
        nd = norb * norb
        Nx, Ny, Nz = 3, 1, 1
        nvol = Nx * Ny * Nz
        V = 0.3 + 0.4j

        onsite = np.zeros((nvol, norb, norb, norb, norb), dtype=complex)
        onsite[:, 0, 0, 1, 1] = V
        onsite[:, 1, 1, 0, 0] = np.conj(V)   # Hermitian partner:
        # W[q,(a,b),(c,d)] must equal conj(W[q,(c,d),(a,b)]), i.e.
        # onsite[a,b,c,d] == conj(onsite[c,d,a,b]).

        # Expected: built by hand, one (a, b, c, d) cell at a time -- the
        # definition of the (a,b)/(c,d) pair-flattening rule, independent
        # of any reshape/transpose call.
        expected = np.zeros((nvol, nd, nd), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                for c in range(norb):
                    for d in range(norb):
                        expected[:, a * norb + b, c * norb + d] = \
                            onsite[:, a, b, c, d]

        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]),
            coeffs={"CoulombInter": np.zeros((1, norb, norb), dtype=complex)})
        W = W_pm_bond(topo, onsite, spatial_shape=(Nx, Ny, Nz))

        assert_approx_array(W, expected, rel=0, abs=1e-15)

        # Self-check that this fixture genuinely discriminates the two
        # groupings (documents the mutant's actual behavior rather than
        # merely asserting it): under (a,c)/(b,d) the SAME (0,0,1,1) entry
        # would land at row=a*norb+c=0*2+1=1, col=b*norb+d=0*2+1=1 -- cell
        # (1, 1) -- not (0, 3), so the correct/wrong mappings disagree at
        # cell (0, 3) (V vs. 0) and at cell (1, 1) (0 vs. V).
        self.assertApprox(expected[0, 0, 3], V, rel=0, abs=0)
        self.assertApprox(expected[0, 1, 1], 0.0, rel=0, abs=0)


class TestWPmBondValidation(ApproxTestCase):
    """Entry-point validation: mesh-injectivity is enforced (via
    ``validate_topology_against_mesh``) and a malformed ``ham_pm_onsite``
    shape is rejected -- both ``ValueError``, never a silent
    misinterpretation."""

    def test_mesh_collision_raises(self):
        # +3 == -3 mod 6 -- the canonical even-mesh self-reversal alias.
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0], [3, 0, 0], [-3, 0, 0]]),
            reverse=np.array([0, 2, 1]),
            coeffs={"CoulombInter": np.zeros((3, 1, 1), dtype=complex)})
        onsite = np.zeros((6 * 4 * 4, 1, 1, 1, 1), dtype=complex)
        with self.assertRaises(ValueError):
            W_pm_bond(topo, onsite, spatial_shape=(6, 4, 4))

    def test_bad_onsite_ndim_rejected(self):
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]),
            coeffs={"CoulombInter": np.zeros((1, 2, 2), dtype=complex)})
        bad = np.zeros((6, 2, 2, 2), dtype=complex)   # ndim=4, not 5
        with self.assertRaises(ValueError):
            W_pm_bond(topo, bad, spatial_shape=(3, 2, 1))

    def test_bad_onsite_nvol_mismatch_rejected(self):
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]),
            coeffs={"CoulombInter": np.zeros((1, 2, 2), dtype=complex)})
        bad = np.zeros((5, 2, 2, 2, 2), dtype=complex)  # nvol=5 != 6
        with self.assertRaises(ValueError):
            W_pm_bond(topo, bad, spatial_shape=(3, 2, 1))

    def test_bad_onsite_non_square_orbital_axes_rejected(self):
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0]]), reverse=np.array([0]),
            coeffs={"CoulombInter": np.zeros((1, 2, 2), dtype=complex)})
        bad = np.zeros((6, 2, 2, 2, 3), dtype=complex)
        with self.assertRaises(ValueError):
            W_pm_bond(topo, bad, spatial_shape=(3, 2, 1))


# =============================================================================
# transverse_bond_bubble_static (Phase W, Task 4)
# =============================================================================
#
# Production tests for hwave.solver.bubble.transverse_bond_bubble_static:
# the static (Omega=0) cross-block bond bubble with the amended
# fwd=block1(down)/rev=block0(up) role contract. The load-bearing pin is
# the V=0 near-machine comparison against tests/test_bond_transverse_ed.py's
# H3 test-local free transverse bond bubble (_transverse_bond_bubble_free,
# imported module-qualified as `_ted`): both sides compose the SAME
# production primitives (bubble._prepare_dense / bubble._bond_pair_full_block)
# with the SAME committed block roles, so agreement is expected at
# near-machine precision (not merely Richardson-extrapolated ED precision,
# which is what H3 itself measures against the ED oracle).


def _transverse_topo_1d(delta_r_1d, reverse, norb):
    """Build a TransverseTopology whose delta_r is the 1-D ring list
    ``delta_r_1d`` (y=z=0 for every channel), for feeding
    ``transverse_bond_bubble_static`` -- only ``delta_r``/``reverse``
    matter to the bubble (``coeffs`` belongs to the VERTEX, W_pm_bond,
    not the bubble); zero-filled coeffs trivially satisfy
    TransverseTopology's Hermitian-closure invariant."""
    B = len(delta_r_1d)
    delta_r = np.zeros((B, 3), dtype=int)
    delta_r[:, 0] = delta_r_1d
    coeffs = {"CoulombInter": np.zeros((B, norb, norb), dtype=complex)}
    return TransverseTopology(
        delta_r=delta_r, reverse=np.array(reverse, dtype=int),
        coeffs=coeffs)


class TestTransverseBondBubbleStaticVZeroPin(ApproxTestCase):
    """The V=0 near-machine pin (spec: "The bond bubble (static)"; plan
    Task 4): ``transverse_bond_bubble_static`` must reproduce
    ``tests/test_bond_transverse_ed.py``'s ``_transverse_bond_bubble_free``
    -- the Phase-H test-local free transverse bond bubble that committed
    the SAME fwd=down/rev=up roles -- at
    ``assert_approx_array(rel=1e-12, abs=1e-13)``, INCLUDING off-diagonal
    ``(m, m')`` cells (both fixtures use the full channel grid, not just
    the diagonal) and the Zeeman-split case (``G_up != G_dn`` -- both
    fixtures below carry a nonzero ``hz``, the same requirement H3's own
    docstring states: a spin-symmetric fixture cannot distinguish the two
    role conventions)."""

    def _check(self, fx, hz, delta_r_1d, reverse, norb, nmat=32):
        green_kw = _ted._h3_free_two_block_green(fx, hz, nmat)
        spatial_shape = (fx.L, 1, 1)
        want = _ted._transverse_bond_bubble_free(
            green_kw, fx.beta, delta_r_1d, norb, spatial_shape)
        topo = _transverse_topo_1d(delta_r_1d, reverse, norb)
        got = bubble_mod.transverse_bond_bubble_static(
            green_kw, None, fx.beta, topo, spatial_shape=spatial_shape)
        assert_approx_array(got, want, rel=1e-12, abs=1e-13)
        return got, want

    def test_norb1_odd_ring_full_channel_grid(self):
        """L=5 ring (odd), norb=1, complex hopping, Zeeman split;
        channels R in {0, +-1, +-2} -- the full 5x5 channel grid,
        exercising every off-diagonal (m, m') cell on an odd ring where R
        and -R are genuinely distinct wrapped channels (mirrors H3's own
        norb=1 fixture)."""
        t = 0.7 * np.exp(0.3j)
        fx = ed_oracle_util.EDFixture(
            L=5, norb=1, t={(0, 0): t}, eps=(0.0,), T=0.5, mu=0.2)
        Rs = [0, 1, -1, 2, -2]
        reverse = [0, 2, 1, 4, 3]
        self._check(fx, 0.15, Rs, reverse, 1)

    def test_norb2_orbital_swap(self):
        """L=3 ring (odd; see H3's scope note on ED tractability), norb=2,
        genuine inter-orbital complex hop, Zeeman split; channels R in
        {0, +-1} -- the multi-orbital off-diagonal (m, m') variant
        (mirrors H3's own norb=2 fixture)."""
        t = {(0, 0): 0.6 + 0.2j, (1, 1): 0.4 - 0.1j,
             (0, 1): 0.15 + 0.05j, (1, 0): 0.15 + 0.05j}
        fx = ed_oracle_util.EDFixture(
            L=3, norb=2, t=t, eps=(0.05, -0.03), T=0.4, mu=0.1)
        Rs = [0, 1, -1]
        reverse = [0, 2, 1]
        self._check(fx, 0.12, Rs, reverse, 2)


class TestTransverseBondBubbleStaticTailEquivalence(ApproxTestCase):
    """``green0_tail=None`` must be numerically equivalent to an explicit
    all-zeros tail array of the same shape/dtype (mirrors
    ``TestBubbleOldVsNewBondStatic``'s/the longitudinal path's own
    None-vs-zeros pin in ``tests/test_bubble_kernel.py``) --
    ``assert_approx_array(rel=0, abs=1e-15)``."""

    def test_none_tail_matches_zero_tail_array(self):
        fx = ed_oracle_util.EDFixture(
            L=4, norb=1, t={(0, 0): 0.5 + 0.1j}, eps=(0.0,), T=0.6, mu=0.1)
        nmat = 16
        green_kw = _ted._h3_free_two_block_green(fx, 0.2, nmat)
        spatial_shape = (fx.L, 1, 1)
        topo = _transverse_topo_1d([0, 1, -1], [0, 2, 1], 1)

        zero_tail = np.zeros_like(green_kw)
        got_none = bubble_mod.transverse_bond_bubble_static(
            green_kw, None, fx.beta, topo, spatial_shape=spatial_shape)
        got_zero = bubble_mod.transverse_bond_bubble_static(
            green_kw, zero_tail, fx.beta, topo, spatial_shape=spatial_shape)
        assert_approx_array(got_zero, got_none, rel=0, abs=1e-15)


class TestTransverseBondBubbleStaticValidation(ApproxTestCase):
    """Entry-point validation: ``nblock != 2``, odd ``nmat``, and a
    mesh-colliding topology are all rejected with ``ValueError`` before
    any FFT work runs."""

    def _small_green(self, nblock, nmat, L=4):
        rng = np.random.default_rng(3)
        shape = (nblock, nmat, L, 1, 1)
        return (rng.normal(size=shape) + 1j * rng.normal(size=shape))

    def test_rejects_nblock_1(self):
        green_kw = self._small_green(1, 8)
        topo = _transverse_topo_1d([0, 1, -1], [0, 2, 1], 1)
        with self.assertRaises(ValueError):
            bubble_mod.transverse_bond_bubble_static(
                green_kw, None, 1.0, topo, spatial_shape=(4, 1, 1))

    def test_rejects_nblock_3(self):
        green_kw = self._small_green(3, 8)
        topo = _transverse_topo_1d([0, 1, -1], [0, 2, 1], 1)
        with self.assertRaises(ValueError):
            bubble_mod.transverse_bond_bubble_static(
                green_kw, None, 1.0, topo, spatial_shape=(4, 1, 1))

    def test_rejects_odd_nmat(self):
        green_kw = self._small_green(2, 7)
        topo = _transverse_topo_1d([0, 1, -1], [0, 2, 1], 1)
        with self.assertRaises(ValueError):
            bubble_mod.transverse_bond_bubble_static(
                green_kw, None, 1.0, topo, spatial_shape=(4, 1, 1))

    def test_mesh_collision_raises(self):
        # +2 == -2 mod 4 -- the canonical even-mesh self-reversal alias.
        green_kw = self._small_green(2, 8, L=4)
        topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0], [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1]),
            coeffs={"CoulombInter": np.zeros((3, 1, 1), dtype=complex)})
        with self.assertRaises(ValueError):
            bubble_mod.transverse_bond_bubble_static(
                green_kw, None, 1.0, topo, spatial_shape=(4, 1, 1))


# =============================================================================
# RPA._run_transverse_bond_pipeline (Phase W, Task 5)
# =============================================================================
#
# The dressing + collapse + the internal pipeline entry
# (RPA._build_ham_pm_onsite / RPA._run_transverse_bond_pipeline in
# rpa.py, appended for this task): on-site vertex assembly ->
# W_pm_bond -> transverse_bond_bubble_static -> the EXISTING _solve_rpa
# dressing -> the (m=0, m'=0) collapse. Anti-vacuity (spec "Phase W --
# implementation + ED campaign"): this entry is called directly by every
# test below, with NO prereq/fallback path.


def _make_ring_ladder_solver(Lx=4, Ly=4, Nmat=8, hz=0.0, U=0.8, T=2.0,
                              norb=2):
    """Minimal ring+ladder RPA solver fixture: ON-SITE CoulombIntra ONLY
    (no off-site declarations anywhere). An Extern field is ALWAYS
    declared (forcing ``solver.spin_mode == "spin-diag"`` and giving
    ``solver.green0`` the required ``nblock=2`` up/down blocks --
    ``transverse_bond_bubble_static`` requires exactly that shape), but
    the default ``hz=0.0`` makes its coefficient zero, so ``G_up ==
    G_down`` EXACTLY (H1 = ham_extern_q * 0 = 0 -> Hnew[0] == Hnew[1]).
    This is deliberate, not incidental: the bond pipeline's ``chi0``
    (``transverse_bond_bubble_static``, fwd=down/rev=up, ``<S+;S+dagger>``)
    and the production plain channel's ``chi0_+-``
    (``_calc_chi0q_transverse``, fwd=up/rev=down, ``<S-;S-dagger>``) are
    PROVABLY the same object only when ``G_up == G_down`` -- off that
    line they are genuinely different physical objects whenever the
    Green's function has off-diagonal orbital structure (a role-swap is
    not a matrix transpose when the two matrices being swapped differ),
    which is exactly what a nonzero ``hz`` combined with the ``norb=2``
    mixing hop below was measured to trigger (diff ~4e-4, pure sign
    flips at the off-diagonal-orbital cells) -- whether production's
    Zeeman-split ``chiq_pm`` genuinely represents ``<S-;S-dagger>`` at
    ``G_up != G_down`` is explicitly the Phase-H/W Task 6 "production
    conjugate-pair granule"'s job, not this gate's. ``hz`` stays a
    parameter (used by the mesh/nblock validation tests below, where
    the value is immaterial) for the same construction to serve both.

    ``norb=2`` (the default) ALSO declares a genuine on-site inter-
    orbital hopping mix, so ``G[a,b]`` is non-diagonal in orbital space
    and ``chiq_pm``'s ``(a,c,b,d)`` tensor is NOT symmetric under an
    arbitrary axis relabelling -- ``norb=1`` cannot discriminate a
    collapse-index mutation (every orbital axis is trivially size-1),
    which is why the collapse mutation check (Task 5) needed this.

    Extends the ``tests/test_bubble_kernel.py`` ``_make_rpa_solver``
    construction pattern (temp Wannier90-like input files, including its
    own ``norb=2`` mixing-hop convention) with ``calc_type=
    'ring+ladder'`` baked into ``info_mode`` at construction time (that
    shared fixture builds an UNSOLVED solver with no ``calc_type``
    parameter, so it cannot drive the production ladder channel), and
    then runs a real ``solver.solve(...)`` -- which populates
    ``green_info["chiq_pm"]`` via the PRODUCTION plain
    ``_build_transverse_channel`` + ``_solve_rpa`` path. Gate W1 compares
    ``_run_transverse_bond_pipeline``'s collapsed output against exactly
    that array.
    """
    import hwave.qlmsio.read_input_k as read_input_k

    t = 0.6 * np.exp(0.25j) if norb >= 2 else 0.6
    mix = 0.35 * np.exp(-0.5j) if norb >= 2 else 0.0
    n_mix = 2 if norb >= 2 else 0
    nr = 2 * norb + n_mix
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
            f.write("{}\n".format(norb))
            for _ in range(norb):
                f.write("0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("gate W1 fixture transfer\n")
            f.write("{}\n{}\n".format(norb, nr))
            f.write(("1 " * nr).strip() + "\n")
            for orb in range(1, norb + 1):
                to = t * orb    # distinct magnitude per orbital
                f.write(" 1 0 0 {o} {o} {re:.12f} {im:.12f}\n".format(
                    o=orb, re=to.real, im=to.imag))
                f.write("-1 0 0 {o} {o} {re:.12f} {im:.12f}\n".format(
                    o=orb, re=np.conj(to).real, im=np.conj(to).imag))
            if norb >= 2:
                f.write(" 0 0 0 1 2 {re:.12f} {im:.12f}\n".format(
                    re=mix.real, im=mix.imag))
                f.write(" 0 0 0 2 1 {re:.12f} {im:.12f}\n".format(
                    re=np.conj(mix).real, im=np.conj(mix).imag))
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("gate W1 fixture CoulombIntra\n")
            f.write("{}\n{}\n".format(norb, norb))
            f.write(("1 " * norb).strip() + "\n")
            for orb in range(1, norb + 1):
                f.write(" 0 0 0 {o} {o} {u:.12f} 0.0\n".format(o=orb, u=U))
        with open(os.path.join(d, "extern.dat"), "w") as f:
            f.write("gate W1 fixture Extern\n")
            f.write("{}\n{}\n".format(norb, norb))
            f.write(("1 " * norb).strip() + "\n")
            for orb in range(1, norb + 1):
                f.write(" 0 0 0 {o} {o} 1.0 0.0\n".format(o=orb))

        interaction = {
            'path_to_input': d,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'Extern': 'extern.dat',
        }
        param = {'T': T, 'mu': 0.0, 'CellShape': [Lx, Ly, 1],
                 'SubShape': [1, 1, 1], 'Nmat': Nmat,
                 'coeff_extern': hz}
        info_input = {'path_to_input': d, 'interaction': interaction}
        info_mode = {'mode': 'RPA', 'param': param,
                     'calc_scheme': 'general', 'calc_type': 'ring+ladder'}

        read_io = read_input_k.QLMSkInput(info_input)
        ham_info = read_io.get_param("ham")
        solver = rpa_mod.RPA(ham_info, {}, info_mode)
        green_info = read_io.get_param("green")

        out_dir = tempfile.mkdtemp(prefix="rpa_gate_w1_")
        solver.solve(green_info, out_dir)
    return solver, green_info


class TestDressingDirectFiniteMatrix(ApproxTestCase):
    """Implementation-independent pin on ``_solve_rpa`` -- exactly the
    call ``_run_transverse_bond_pipeline``'s dressing step makes
    (``_solve_rpa(chi0.reshape(1, nvol, ND, ND), W)``) -- against a
    PLAIN ``numpy.linalg.solve(I + chi0[q] @ W[q], chi0[q])`` reference
    built per q with NO call to ``_solve_rpa`` (spec "The dressing
    (static)"). Random NONSYMMETRIC complex matrices: the ``I + chi0.ham``
    convention has no reason to require symmetry, and a wrong transpose
    orientation would only show up on a nonsymmetric fixture."""

    def test_dressing_matches_plain_solve_per_q(self):
        rng = np.random.default_rng(11)
        nvol, ND = 3, 4
        chi0 = (rng.normal(size=(nvol, ND, ND))
                + 1j * rng.normal(size=(nvol, ND, ND)))
        W = (rng.normal(size=(nvol, ND, ND))
             + 1j * rng.normal(size=(nvol, ND, ND)))

        solver = object.__new__(rpa_mod.RPA)
        solver.lattice = types.SimpleNamespace(nvol=nvol)

        sol = solver._solve_rpa(chi0.reshape(1, nvol, ND, ND), W)
        got = sol[0]

        eye = np.eye(ND, dtype=complex)
        want = np.zeros_like(chi0)
        for q in range(nvol):
            want[q] = np.linalg.solve(eye + chi0[q] @ W[q], chi0[q])

        assert_approx_array(got, want, rel=0, abs=1e-12)


class TestCollapseRule(ApproxTestCase):
    """Elementwise pin on the collapse step (spec "Collapse rule
    (exact)"): ``chiq_pm_static[q, a, c, b, d] = chi_pm_bond[q, (0, a,
    c), (0, b, d)]``. ``W_pm_bond``/``transverse_bond_bubble_static``/
    ``_solve_rpa`` are mocked out so the dressed bond object
    (``chi_pm_bond``) is a FULLY CONTROLLED, non-degenerate random array
    -- the expected collapse is built by DIRECT per-index selection
    (never via a ``.reshape``/``.transpose`` call), so this test does not
    validate production against itself."""

    def test_collapse_picks_m0_subblock_by_direct_indexing(self):
        norb, B, nvol = 2, 2, 3
        nd = norb * norb
        ND = B * nd
        rng = np.random.default_rng(23)
        chi_pm_bond = (rng.normal(size=(nvol, ND, ND))
                       + 1j * rng.normal(size=(nvol, ND, ND)))

        solver = object.__new__(rpa_mod.RPA)
        solver.lattice = types.SimpleNamespace(nvol=nvol, shape=(nvol, 1, 1))
        solver.norb = norb
        solver.fft_workers = 1

        dummy_onsite = np.zeros((nvol, norb, norb, norb, norb),
                                 dtype=complex)
        dummy_topo = object()

        with mock.patch.object(
                rpa_mod.RPA, "_build_ham_pm_onsite",
                return_value=dummy_onsite), \
             mock.patch.object(
                rpa_mod.bond_channels, "W_pm_bond",
                return_value=np.zeros((nvol, ND, ND), dtype=complex)), \
             mock.patch.object(
                rpa_mod.bubble, "transverse_bond_bubble_static",
                return_value=np.zeros((nvol, ND, ND), dtype=complex)), \
             mock.patch.object(
                rpa_mod.RPA, "_solve_rpa",
                return_value=chi_pm_bond.reshape(1, nvol, ND, ND)):
            chi_bond_out, collapsed = solver._run_transverse_bond_pipeline(
                green_kw=None, green0_tail=None, beta=1.0, topo=dummy_topo)

        assert_approx_array(chi_bond_out, chi_pm_bond, rel=0, abs=0)

        expected = np.zeros((nvol, norb, norb, norb, norb), dtype=complex)
        for q in range(nvol):
            for a in range(norb):
                for c in range(norb):
                    for b in range(norb):
                        for d in range(norb):
                            row = a * norb + c    # channel m=0 offset 0
                            col = b * norb + d
                            expected[q, a, c, b, d] = chi_pm_bond[q, row, col]

        assert_approx_array(collapsed, expected, rel=0, abs=0)

        # Self-check that the fixture genuinely discriminates the
        # (a,c)/(b,d) grouping from the (a,b)/(c,d) one this test's
        # mutation check targets: with norb=2 the two groupings read
        # DIFFERENT cells whenever a != c or b != d, e.g. (a,b,c,d) =
        # (0,1,0,0): correct row=a*norb+c=0, col=b*norb+d=2; the (a,b)/
        # (c,d) grouping would read row=a*norb+b=1, col=c*norb+d=0. A
        # REAL runtime discrimination check (the previous version here
        # only compared the two literal integers 0 and 1, which holds
        # unconditionally and exercises neither the fixture nor the
        # array -- vacuous): read BOTH candidate cells from the actual
        # random ``chi_pm_bond`` fixture at q=0 and require them to
        # differ, so this test would genuinely fail to discriminate a
        # collapse-index mutation if the fixture ever degenerated (e.g.
        # an accidentally symmetric array) rather than passing no
        # matter what.
        a_ex, b_ex, c_ex, d_ex = 0, 1, 0, 0
        correct_cell = chi_pm_bond[0, a_ex * norb + c_ex, b_ex * norb + d_ex]
        wrong_cell = chi_pm_bond[0, a_ex * norb + b_ex, c_ex * norb + d_ex]
        self.assertNotEqual(correct_cell, wrong_cell)


class TestGateW1OnsiteReduction(ApproxTestCase):
    """**Gate W1** (spec "Gate W1 (on-site reduction, exact scope)"):
    with ONLY on-site declarations, ``_run_transverse_bond_pipeline``'s
    collapsed static output must equal the Omega=0 slice
    (``l = nmat // 2``) of today's PLAIN ``chiq_pm`` (the production
    ``_build_transverse_channel`` + ``_solve_rpa`` result), at
    ``assert_approx_array(rel=1e-12, abs=1e-13)``. Both sides solve the
    SAME physics through DIFFERENT code paths (the plain path builds its
    on-site vertex/bubble from the FULL solver machinery in one shot;
    the bond pipeline builds ``ham_pm_onsite`` from filtered pre-fold
    declarations and the bubble via the bond-channel machinery with
    B=1), so agreement here is the load-bearing algebraic check that the
    bond lift reduces exactly to the existing plain channel when there
    is nothing off-site to represent.

    Uses the DEFAULT ``_make_ring_ladder_solver`` fixture (``hz=0.0``,
    ``G_up == G_down``): the two sides' bubbles are provably the SAME
    object only on that line (see the fixture's own docstring) -- off it
    they diverge by a measured ~4e-4 that is a genuine convention
    difference between ``<S+;S+dagger>`` and ``<S-;S-dagger>``, not a
    defect gate W1 is scoped to catch (that adjudication is Task 6's
    production conjugate-pair granule)."""

    def test_collapsed_output_matches_plain_chiq_pm_at_omega_zero(self):
        solver, green_info = _make_ring_ladder_solver()
        beta = 1.0 / solver.T

        # B=1 topology: resolve_transverse_topology reads OFF-SITE
        # content only (on-site entries are ignored -- on-site content
        # is represented exclusively through ham_pm_onsite), and this
        # fixture declares none, so every type's coeffs are all-zero at
        # the single mandatory channel 0.
        topo = resolve_transverse_topology(
            solver.ham_info.param_ham, np.eye(3), solver.norb)
        self.assertEqual(len(topo.delta_r), 1)

        chi_pm_bond_static, chiq_pm_static = (
            solver._run_transverse_bond_pipeline(
                solver.green0, solver.green0_tail, beta, topo))

        chiq_pm_plain = green_info["chiq_pm"]
        nmat = chiq_pm_plain.shape[0]
        want = chiq_pm_plain[nmat // 2]

        assert_approx_array(chiq_pm_static, want, rel=1e-12, abs=1e-13)

        # Report the measured margin (spec: "measured margin reported").
        margin = float(np.max(np.abs(chiq_pm_static - want)))
        print("Gate W1 measured max-abs margin: {:.3e}".format(margin))


class TestRunTransverseBondPipelineValidation(ApproxTestCase):
    """Entry-point validation: both failure modes PROPAGATE as
    ``ValueError`` from the underlying calls -- ``_run_transverse_bond_
    pipeline`` adds no guards of its own (no prereq fallback)."""

    def test_mesh_colliding_topology_raises(self):
        solver, _green_info = _make_ring_ladder_solver()
        beta = 1.0 / solver.T
        # +2 == -2 mod 4 -- solver.lattice.shape is (4, 4, 1) by default.
        bad_topo = TransverseTopology(
            delta_r=np.array([[0, 0, 0], [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1]),
            coeffs={"CoulombInter": np.zeros(
                (3, solver.norb, solver.norb), dtype=complex)})
        with self.assertRaises(ValueError):
            solver._run_transverse_bond_pipeline(
                solver.green0, solver.green0_tail, beta, bad_topo)

    def test_nblock_not_two_raises(self):
        solver, _green_info = _make_ring_ladder_solver()
        beta = 1.0 / solver.T
        topo = resolve_transverse_topology(
            solver.ham_info.param_ham, np.eye(3), solver.norb)
        bad_green = solver.green0[:1]          # nblock=1, not 2
        bad_tail = solver.green0_tail[:1]
        with self.assertRaises(ValueError):
            solver._run_transverse_bond_pipeline(
                bad_green, bad_tail, beta, topo)



# =============================================================================
# Phase W, Task 7: the experimental gate (transverse_bond_channels),
# provenance, preflight, and the PairHop off-site warning.
#
# Config key names, rejection wording sources, provenance keys, and the
# config matrix are per docs/superpowers/specs/2026-08-15-bond-transverse-
# design.md "Production surface (experimental gate)"; the wiring pattern
# (config parsing site, prereq validator, resource preflight, npz block,
# stale-option warnings) mirrors sc.py's [eliashberg] bond_channels gate
# (_read_bond_config / _validate_bond_prereqs / _bond_resource_preflight)
# for its transverse sibling. Every test below drives the PUBLIC entry
# (the RPA(...) constructor / solve() / save_results()), never a private
# helper directly, except where a helper IS the object under test
# (_transverse_bond_resource_preflight's own formula).
# =============================================================================

def _make_bond_gate_fixture(*, transverse_bond_channels=None,
                             calc_type="ring+ladder", calc_scheme="general",
                             interactions=None, param_overrides=None,
                             Lx=4, Ly=4, Nmat=32, T=2.0, filling=0.75,
                             input_path='tests/rpa/input'):
    """norb=1 RPA fixture built via the PUBLIC ``RPA(...)`` constructor
    from ``tests/rpa/input``'s existing Wannier90-format files -- the
    SAME fixture family ``tests/test_rpa_ladder.py``'s ``_run_rpa`` (and
    its ``test_offsite_two_body_is_rejected``) uses: on-site CoulombIntra
    (U=4) and, by default, off-site CoulombInter (V=1 at
    R=(+-1,0,0)/(0,+-1,0) -- an ACTIVE transverse-resolved off-site
    declaration the bond gate can represent but the plain ladder's
    representability check (``_check_transverse_representable``)
    rejects) plus Extern (forces ``solver.spin_mode == 'spin-diag'``,
    which the bond pipeline's up/down Green blocks require).

    Does NOT call ``solve()`` -- callers decide (so rejection tests can
    ``assertRaises`` around ``solve()`` itself, the public entry point).

    Returns ``(solver, green_info, out_dir)``.
    """
    import hwave.qlmsio.read_input_k as read_input_k

    if interactions is None:
        interactions = {'CoulombIntra': 'coulombintra.dat',
                         'CoulombInter': 'coulombinter.dat',
                         'Extern': 'extern.dat'}
    interaction_dict = {
        'path_to_input': input_path,
        'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat',
    }
    interaction_dict.update(interactions)

    param = {'T': T, 'filling': filling, 'CellShape': [Lx, Ly, 1],
             'SubShape': [1, 1, 1], 'Nmat': Nmat}
    if transverse_bond_channels is not None:
        param['transverse_bond_channels'] = transverse_bond_channels
    if param_overrides:
        param.update(param_overrides)

    info_mode = {'mode': 'RPA', 'param': param,
                 'calc_scheme': calc_scheme, 'calc_type': calc_type}
    info_input = {'path_to_input': input_path,
                  'interaction': interaction_dict}

    read_io = read_input_k.QLMSkInput(info_input)
    ham_info = read_io.get_param("ham")
    solver = rpa_mod.RPA(ham_info, {}, info_mode)
    green_info = read_io.get_param("green")
    out_dir = tempfile.mkdtemp(prefix="rpa_bond_gate_out_")
    return solver, green_info, out_dir


def _make_spinful_bond_gate_fixture(T=0.5, mu=0.2, Lx=4, Nmat=32, U=0.3,
                                     thop=0.7 + 0.3j, lso=0.35 + 0.15j):
    """norb_phys=1 ``enable_spin_orbital`` ring+ladder fixture with a
    genuine spin-mixing transfer term (off-diagonal generalized-index
    hopping ``lso``), reaching ``solver.spin_mode == 'spinful'`` -- the
    combination the bond gate's prereq validator must reject (spec
    ``spin_mode="spinful"`` -> reject until Phase S). Mirrors
    ``tests/test_rpa_spinful_vertex_exchange.py``'s ``TestVertexExtraction``
    / ``_run_rpa`` construction pattern, the SAME fixture family that
    already established this combination reaches ``spin_mode ==
    'spinful'`` (that file's ``TestSpinConservingLimits`` pins it at
    ``calc_type='ring'``; nothing about the H0 that sets ``spin_mode``
    depends on ``calc_type``).

    Does NOT call ``solve()``. Returns ``(solver, green_info, out_dir)``.
    """
    import hwave.qlmsio.read_input_k as read_input_k

    d = tempfile.mkdtemp(prefix="rpa_bond_gate_spinful_")
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n2\n")
        f.write("0.0 0.0 0.0\n0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n2\n3\n1 1 1\n")
        for i in (1, 2):
            f.write(" 1 0 0 %d %d %.12f %.12f\n"
                     % (i, i, thop.real, thop.imag))
            f.write("-1 0 0 %d %d %.12f %.12f\n"
                     % (i, i, np.conj(thop).real, np.conj(thop).imag))
        f.write(" 0 0 0 1 2 %.12f %.12f\n" % (lso.real, lso.imag))
        f.write(" 0 0 0 2 1 %.12f %.12f\n" % (lso.real, -lso.imag))
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("hdr\n1\n1\n1\n")
        f.write(" 0 0 0 1 1 %.12f 0.0\n" % U)
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}

    param = {"T": T, "mu": mu, "CellShape": [Lx, 1, 1],
             "SubShape": [1, 1, 1], "Nmat": Nmat, "coeff_tail": 1.0,
             "transverse_bond_channels": True}
    info_mode = {"mode": "RPA", "param": param,
                 "enable_spin_orbital": True, "calc_scheme": "general",
                 "calc_type": "ring+ladder"}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
    green_info = io.get_param("green")
    out_dir = tempfile.mkdtemp(prefix="rpa_bond_gate_spinful_out_")
    return solver, green_info, out_dir


def _make_max_shells_fixture(*, transverse_bond_channels=True,
                              max_shells=None, Lx=6, Nmat=8, T=2.0,
                              filling=0.5):
    """norb=1 ring+ladder fixture with off-site CoulombInter declared at
    TWO distinct shells on a Lx=6 chain (``delta_r mod 6`` is injective
    for ``{0, +-1, +-2}``, so no mesh collision): R=+-1 carries a
    NONZERO coefficient (V1=0.4), R=+-2 is DECLARED but exactly zero
    (V2=0.0) -- a truncation that drops only the zero shell is not
    ambiguous (``resolve_transverse_topology`` only refuses a truncation
    that would drop declared NONZERO content), so
    ``transverse_bond_max_shells=1`` has an OBSERVABLE, non-rejected
    effect: it keeps the R=+-1 shell (B=3) and drops the (identically
    zero) R=+-2 shell that the default ``max_shells=None`` keeps (B=5).

    Does NOT call ``solve()``. Returns ``(solver, green_info, out_dir)``.
    """
    import hwave.qlmsio.read_input_k as read_input_k

    d = tempfile.mkdtemp(prefix="rpa_bond_gate_maxshells_")
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n1\n2\n1 1\n"
                 " 1 0 0 1 1 0.5 0.0\n-1 0 0 1 1 0.5 0.0\n")
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("hdr\n1\n1\n1\n 0 0 0 1 1 2.0 0.0\n")
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("hdr\n1\n4\n1 1 1 1\n"
                 " 1 0 0 1 1 0.4 0.0\n-1 0 0 1 1 0.4 0.0\n"
                 " 2 0 0 1 1 0.0 0.0\n-2 0 0 1 1 0.0 0.0\n")
    with open(os.path.join(d, "extern.dat"), "w") as f:
        f.write("hdr\n1\n1\n1\n 0 0 0 1 1 1.0 0.0\n")

    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat",
             "CoulombInter": "coulombinter.dat", "Extern": "extern.dat"}
    param = {"T": T, "filling": filling, "CellShape": [Lx, 1, 1],
             "SubShape": [1, 1, 1], "Nmat": Nmat}
    if transverse_bond_channels is not None:
        param["transverse_bond_channels"] = transverse_bond_channels
    if max_shells is not None:
        param["transverse_bond_max_shells"] = max_shells

    info_mode = {"mode": "RPA", "param": param, "calc_scheme": "general",
                 "calc_type": "ring+ladder"}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
    green_info = io.get_param("green")
    out_dir = tempfile.mkdtemp(prefix="rpa_bond_gate_maxshells_out_")
    return solver, green_info, out_dir


class TestTransverseBondGateConfig(ApproxTestCase):
    """``[mode.param] transverse_bond_channels`` config parsing (spec
    "Production surface"): default false, parsed ONLY under
    ``calc_type='ring+ladder'``, stale-option warnings, and type
    validation -- mirroring sc.py's ``[eliashberg] bond_channels``
    ``_read_bond_config`` pattern for the (eliashberg) longitudinal bond
    gate."""

    ONSITE_ONLY = {'CoulombIntra': 'coulombintra.dat', 'Extern': 'extern.dat'}

    def test_absent_defaults_to_false(self):
        solver, _gi, _out = _make_bond_gate_fixture(
            transverse_bond_channels=None, interactions=self.ONSITE_ONLY)
        self.assertIs(solver.transverse_bond_channels, False)

    def test_gate_absent_equals_gate_false_bitwise(self):
        """Config matrix (spec): "gate absent == gate false (bitwise-
        identical run)"."""
        solver_a, gi_a, out_a = _make_bond_gate_fixture(
            transverse_bond_channels=None, interactions=self.ONSITE_ONLY)
        solver_a.solve(gi_a, out_a)
        solver_b, gi_b, out_b = _make_bond_gate_fixture(
            transverse_bond_channels=False, interactions=self.ONSITE_ONLY)
        solver_b.solve(gi_b, out_b)

        self.assertTrue(np.array_equal(gi_a["chiq"], gi_b["chiq"]))
        self.assertTrue(np.array_equal(gi_a["chiq_pm"], gi_b["chiq_pm"]))
        self.assertNotIn("chiq_pm_bond_static", gi_a)
        self.assertNotIn("chiq_pm_bond_static", gi_b)

    def test_stale_under_calc_type_ring_warns(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            solver, _gi, _out = _make_bond_gate_fixture(
                calc_type="ring", transverse_bond_channels=True,
                interactions=self.ONSITE_ONLY)
        self.assertIs(solver.transverse_bond_channels, False)
        self.assertTrue(any(
            "transverse_bond_channels" in m and "calc_type" in m
            for m in cm.output))

    def test_stale_options_while_gate_false_warn(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _solver, _gi, _out = _make_bond_gate_fixture(
                transverse_bond_channels=False,
                param_overrides={'transverse_bond_max_shells': 1},
                interactions=self.ONSITE_ONLY)
        self.assertTrue(
            any("transverse_bond_max_shells" in m for m in cm.output))

    def test_flag_rejects_non_bool(self):
        with self.assertRaises(ValueError):
            _make_bond_gate_fixture(
                param_overrides={'transverse_bond_channels': 1},
                interactions=self.ONSITE_ONLY)

    def test_max_shells_rejects_negative(self):
        with self.assertRaises(ValueError):
            _make_bond_gate_fixture(
                transverse_bond_channels=True,
                param_overrides={'transverse_bond_max_shells': -1},
                interactions=self.ONSITE_ONLY)

    def test_max_shells_rejects_non_integral(self):
        with self.assertRaises(ValueError):
            _make_bond_gate_fixture(
                transverse_bond_channels=True,
                param_overrides={'transverse_bond_max_shells': 1.5},
                interactions=self.ONSITE_ONLY)

    def test_memory_cap_rejects_non_positive(self):
        with self.assertRaises(ValueError):
            _make_bond_gate_fixture(
                transverse_bond_channels=True,
                param_overrides={'transverse_bond_memory_cap_gb': 0.0},
                interactions=self.ONSITE_ONLY)

    # -----------------------------------------------------------------
    # Case-insensitivity (CLAUDE.md rule; PR #128 had 6 silent-divergence
    # sites in exactly this defect class): every new config key lookup
    # must be case-robust, driven through the PUBLIC entry -- never by
    # inspecting self.param_mod directly, since self.param_mod being a
    # CaseInsensitiveDict is the implementation detail under test, not
    # something the test may lean on as its own proof.
    # -----------------------------------------------------------------

    ACTIVE = {'CoulombIntra': 'coulombintra.dat',
              'CoulombInter': 'coulombinter.dat', 'Extern': 'extern.dat'}

    def test_mixed_case_gate_key_activates_through_public_entry(self):
        solver_lower, gi_lower, out_lower = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver_lower.solve(gi_lower, out_lower)

        solver_mixed, gi_mixed, out_mixed = _make_bond_gate_fixture(
            transverse_bond_channels=None,
            param_overrides={'Transverse_Bond_Channels': True},
            interactions=self.ACTIVE)
        self.assertIs(solver_mixed.transverse_bond_channels, True)
        solver_mixed.solve(gi_mixed, out_mixed)

        self.assertIn("chiq_pm_bond_static", gi_mixed)
        self.assertNotIn("chiq_pm", gi_mixed)
        self.assertTrue(np.array_equal(gi_lower["chiq"], gi_mixed["chiq"]))
        self.assertTrue(np.array_equal(
            gi_lower["chiq_pm_bond_static"], gi_mixed["chiq_pm_bond_static"]))
        self.assertTrue(np.array_equal(
            gi_lower["chiq_pm_static"], gi_mixed["chiq_pm_static"]))

    def test_mixed_case_flag_absent_lowercase_present_does_not_confuse(self):
        # An UPPER-cased false-y flag must behave exactly like the
        # lowercase false-y flag (the documented default) -- guards
        # against a defect where only the TRUE branch was made
        # case-robust.
        solver, _gi, _out = _make_bond_gate_fixture(
            transverse_bond_channels=None,
            param_overrides={'TRANSVERSE_BOND_CHANNELS': False},
            interactions=self.ONSITE_ONLY)
        self.assertIs(solver.transverse_bond_channels, False)

    def test_mixed_case_max_shells_and_memory_cap_keys_thread_through(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=None,
            param_overrides={'TRANSVERSE_bond_channels': True,
                              'Transverse_Bond_Max_Shells': 1,
                              'transverse_bond_MEMORY_cap_gb': 2.5},
            interactions=self.ACTIVE)
        self.assertIs(solver.transverse_bond_channels, True)
        self.assertEqual(solver.transverse_bond_max_shells, 1)
        self.assertEqual(solver.transverse_bond_memory_cap_gb, 2.5)
        solver.solve(gi, out)  # must not raise
        self.assertIn("chiq_pm_bond_static", gi)


class TestTransverseBondGatePrereqs(ApproxTestCase):
    """``_validate_transverse_bond_prereqs``'s REJECT list (spec
    "Production surface"): no active off-site declaration,
    ``spin_mode='spinful'``, and an externally supplied ``chi0q_init``.
    Also proves the representability-check bypass is unreachable when
    the gate is off (spec: "a config test proves the bypass is
    unreachable when the gate is off"), and that it DOES engage on the
    same input with the gate on."""

    ACTIVE = {'CoulombIntra': 'coulombintra.dat',
              'CoulombInter': 'coulombinter.dat', 'Extern': 'extern.dat'}
    ONSITE_ONLY = {'CoulombIntra': 'coulombintra.dat', 'Extern': 'extern.dat'}

    def test_no_active_declaration_rejected(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ONSITE_ONLY)
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("active off-site", str(cm.exception))

    def test_bypass_engages_on_gate_and_solves(self):
        """Gate on: off-site CoulombInter -- unrepresentable by the plain
        vertex -- is REPRESENTED by the bond path instead of being
        rejected."""
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver.solve(gi, out)
        self.assertIn("chiq_pm_bond_static", gi)
        self.assertIn("chiq_pm_static", gi)
        self.assertNotIn("chiq_pm", gi)

    def test_bypass_unreachable_with_gate_off(self):
        """SAME off-site declaration, gate off:
        ``_check_transverse_representable`` still runs and rejects it --
        the bypass branch in ``solve()`` is unreachable when
        ``transverse_bond_channels`` is false."""
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=False, interactions=self.ACTIVE)
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("does not support off-site", str(cm.exception))

    def test_spinful_rejected(self):
        solver, gi, out = _make_spinful_bond_gate_fixture()
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("spin_mode='spinful'", str(cm.exception))

    def test_sublattice_rejected(self):
        """A sublattice run (SubShape < CellShape, has_sublattice=True)
        must be rejected at the gate: ``_validate_transverse_bond_prereqs``
        resolves the topology from the PRE-FOLD ``param_ham_orig``
        declarations but passes the FOLDED ``self.norb`` and the pipeline
        then runs on the supercell mesh -- the original-cell (R, a, b)
        records are never mapped to supercell coordinates, so a
        sublattice run would silently compute the wrong vertices without
        this guard (review fix)."""
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE,
            param_overrides={'SubShape': [2, 1, 1]})
        self.assertTrue(solver.lattice.has_sublattice)
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("sublattice", str(cm.exception))

    def test_external_chi0q_init_rejected(self):
        # Build a spin-free chi0q with a plain ring solve, then feed it
        # back in as chi0q_init to a gate-on ring+ladder solver. The
        # pre-existing spin-diag/external guard (rpa.py, "a same-instance
        # green0 may belong to an OLDER bubble") does not fire for a
        # spin-free chi0q, so this exercises the gate's OWN
        # external-chi0q rejection specifically.
        prep, gi0, out0 = _make_bond_gate_fixture(
            calc_type="ring", transverse_bond_channels=None,
            interactions={'CoulombIntra': 'coulombintra.dat',
                          'CoulombInter': 'coulombinter.dat'})
        prep.solve(gi0, out0)
        self.assertEqual(prep.spin_mode, "spin-free")
        prep.save_results(
            {'path_to_output': out0, 'chi0q': 'chi0q.npz'}, gi0)

        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True,
            interactions={'CoulombIntra': 'coulombintra.dat',
                          'CoulombInter': 'coulombinter.dat'})
        gi.update(solver.read_init(
            {'path_to_input': out0, 'chi0q_init': 'chi0q.npz'}))
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("externally supplied chi0q", str(cm.exception))


class TestTransverseBondMaxShellsThreading(ApproxTestCase):
    """A positive ``transverse_bond_max_shells`` must thread end to end:
    ``_init_transverse_bond_config`` -> ``self.transverse_bond_max_shells``
    -> ``_validate_transverse_bond_prereqs``'s ``resolve_transverse_
    topology(..., max_shells=...)`` call -- with an OBSERVABLE effect on
    the resulting topology/B (not just the rejection paths and the
    default-None case, which is all the earlier config tests covered)."""

    def test_positive_max_shells_truncates_the_resolved_topology(self):
        solver_trunc, gi_trunc, out_trunc = _make_max_shells_fixture(
            max_shells=1)
        solver_trunc.solve(gi_trunc, out_trunc)
        self.assertEqual(solver_trunc.transverse_bond_max_shells, 1)
        topo_trunc = solver_trunc._transverse_bond_topo
        self.assertEqual(len(topo_trunc.delta_r), 3)  # onsite + R=+-1

        # Matches a DIRECT resolve_transverse_topology call with the SAME
        # max_shells on the solver's own (pre-fold) declarations.
        has_sub = getattr(solver_trunc.lattice, "has_sublattice", False)
        interactions = (solver_trunc.ham_info.param_ham_orig if has_sub
                        else solver_trunc.ham_info.param_ham)
        topo_direct = resolve_transverse_topology(
            interactions, np.eye(3), solver_trunc.norb, max_shells=1)
        self.assertTrue(np.array_equal(topo_trunc.delta_r,
                                       topo_direct.delta_r))
        self.assertTrue(np.array_equal(topo_trunc.reverse,
                                       topo_direct.reverse))
        for t in topo_trunc.coeffs:
            self.assertTrue(np.array_equal(topo_trunc.coeffs[t],
                                           topo_direct.coeffs[t]))

        # DIFFERS from the default (max_shells=None, keep every shell):
        # a value silently failing to thread through would produce the
        # SAME (untruncated) B here as above.
        solver_full, gi_full, out_full = _make_max_shells_fixture(
            max_shells=None)
        solver_full.solve(gi_full, out_full)
        self.assertIsNone(solver_full.transverse_bond_max_shells)
        topo_full = solver_full._transverse_bond_topo
        self.assertEqual(len(topo_full.delta_r), 5)  # onsite + R=+-1,+-2
        self.assertGreater(len(topo_full.delta_r), len(topo_trunc.delta_r))

        # And through the npz output: the channel table + provenance key
        # both reflect the truncated (not the default) topology.
        solver_trunc.save_results(
            {'path_to_output': out_trunc, 'chiq': 'chiq.npz'}, gi_trunc)
        data = np.load(os.path.join(out_trunc, 'chiq.npz'),
                       allow_pickle=True)
        self.assertEqual(int(data["transverse_bond_max_shells"]), 1)
        self.assertEqual(tuple(data["transverse_bond_delta_r"].shape),
                          (3, 3))
        self.assertTrue(np.array_equal(
            data["transverse_bond_delta_r"], topo_trunc.delta_r))

    def test_max_shells_dropping_declared_nonzero_content_raises_via_solve(
            self):
        """Public-entry pin (review fix) of the guard
        ``resolve_transverse_topology``'s own truncation logic already
        enforces (directly exercised by
        ``TestResolveTransverseTopology.
        test_max_shells_dropping_a_declared_nonzero_shell_raises``): the
        SAME fixture family as
        ``test_positive_max_shells_truncates_the_resolved_topology``, but
        with ``max_shells=0`` -- which would drop the R=+-1 shell that
        carries a genuinely NONZERO declared CoulombInter coefficient
        (V1=0.4) -- must raise through the full ``solve()`` public entry
        point, not merely when ``resolve_transverse_topology`` is called
        directly."""
        solver, gi, out = _make_max_shells_fixture(max_shells=0)
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("max_shells=0", str(cm.exception))
        self.assertIn("declared nonzero", str(cm.exception))


class TestTransverseBondResourcePreflight(ApproxTestCase):
    """``_transverse_bond_resource_preflight`` (spec "Phase W -- Budget
    (stated, not deferred)"): the ``(3+K_solve)*Nq*ND**2*16`` byte
    estimate against ``transverse_bond_memory_cap_gb`` (REFUSE above the
    cap), and the ``Nq*ND**3 > 1e12`` op-count WARNING (never a
    refusal)."""

    ACTIVE = {'CoulombIntra': 'coulombintra.dat',
              'CoulombInter': 'coulombinter.dat', 'Extern': 'extern.dat'}

    def test_tiny_cap_refuses(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True,
            param_overrides={'transverse_bond_memory_cap_gb': 1e-12},
            interactions=self.ACTIVE)
        with self.assertRaises(ValueError) as cm:
            solver.solve(gi, out)
        self.assertIn("memory_cap_gb", str(cm.exception))

    def test_default_cap_passes_for_tiny_fixture(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver.solve(gi, out)  # must not raise
        self.assertIn("chiq_pm_bond_static", gi)

    def test_formula_matches_documented_estimate(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver.solve(gi, out)
        topo = solver._transverse_bond_topo
        B = len(topo.delta_r)
        ND = B * solver.norb ** 2
        Nq = solver.lattice.nvol
        K_solve = 3
        want_peak = (3 + K_solve) * Nq * ND ** 2 * 16

        solver.transverse_bond_memory_cap_gb = (want_peak * 0.999) / 1.0e9
        with self.assertRaises(ValueError):
            solver._transverse_bond_resource_preflight(topo)

        solver.transverse_bond_memory_cap_gb = (want_peak * 1.001) / 1.0e9
        solver._transverse_bond_resource_preflight(topo)  # must not raise

    def test_op_count_warns_not_refuses(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver.solve(gi, out)
        topo = solver._transverse_bond_topo

        solver.transverse_bond_memory_cap_gb = 1.0e12  # never refuses here
        orig_nvol = solver.lattice.nvol
        solver.lattice.nvol = 10 ** 11  # pushes Nq*ND**3 well past 1e12
        try:
            with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
                solver._transverse_bond_resource_preflight(topo)
        finally:
            solver.lattice.nvol = orig_nvol
        self.assertTrue(any("op-count" in m for m in cm.output))

    def test_preflight_rejection_precedes_any_solve_call(self):
        """Call-order pin (review fix): the resource preflight must run
        BEFORE the (expensive) longitudinal solve, so a cap rejection
        fires before ``_solve_rpa`` is ever invoked -- not merely before
        the bond-resolved solve, but before the plain ring solve too."""
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True,
            param_overrides={'transverse_bond_memory_cap_gb': 1e-12},
            interactions=self.ACTIVE)
        with mock.patch.object(
                rpa_mod.RPA, "_solve_rpa",
                wraps=rpa_mod.RPA._solve_rpa) as spy:
            with self.assertRaises(ValueError) as cm:
                solver.solve(gi, out)
        self.assertIn("memory_cap_gb", str(cm.exception))
        spy.assert_not_called()


class TestTransverseBondGateOutput(ApproxTestCase):
    """End-to-end smoke test on a tiny fixture (anti-vacuity): config
    dict -> ``RPA(...)`` constructor -> ``solve()`` -> ``save_results()``
    -> npz written; keys, shapes, schema metadata, and a value PINNED
    against a direct ``_run_transverse_bond_pipeline`` call on the same
    fixture (spec "Output (static-only)")."""

    ACTIVE = {'CoulombIntra': 'coulombintra.dat',
              'CoulombInter': 'coulombinter.dat', 'Extern': 'extern.dat'}

    def test_npz_keys_shapes_schema_and_cross_pin(self):
        solver, gi, out = _make_bond_gate_fixture(
            transverse_bond_channels=True, interactions=self.ACTIVE)
        solver.solve(gi, out)
        solver.save_results(
            {'path_to_output': out, 'chiq': 'chiq.npz'}, gi)

        data = np.load(os.path.join(out, 'chiq.npz'), allow_pickle=True)

        self.assertNotIn("chiq_pm", data.files)
        for key in ("chiq_pm_bond_static", "chiq_pm_static",
                    "transverse_bond_schema_version",
                    "transverse_output_kind",
                    "transverse_bond_delta_r", "transverse_bond_reverse",
                    "transverse_bond_index_order",
                    "transverse_bond_max_shells",
                    "transverse_spatial_shape", "transverse_q_convention",
                    "transverse_spin_mode", "transverse_normalization"):
            self.assertIn(key, data.files)

        self.assertEqual(int(data["transverse_bond_schema_version"]), 1)
        self.assertEqual(str(data["transverse_output_kind"]), "bond_static")
        self.assertEqual(str(data["transverse_bond_index_order"]),
                          "m*norb**2 + a*norb + b")
        self.assertEqual(str(data["transverse_q_convention"]),
                          "q_d = 2*pi*n_d/N_d, C-order flattening")
        self.assertEqual(str(data["transverse_spin_mode"]), "spin-diag")
        self.assertEqual(str(data["transverse_normalization"]),
                          "per-site, 1/sqrt(Nvol) bilinears")
        self.assertEqual(int(data["transverse_bond_max_shells"]), -1)

        topo = solver._transverse_bond_topo
        B = len(topo.delta_r)
        ND = B * solver.norb ** 2
        Nq = solver.lattice.nvol
        self.assertEqual(tuple(data["chiq_pm_bond_static"].shape),
                          (Nq, ND, ND))
        self.assertEqual(tuple(data["chiq_pm_static"].shape),
                          (Nq, solver.norb, solver.norb,
                           solver.norb, solver.norb))
        self.assertEqual(tuple(data["transverse_bond_delta_r"].shape),
                          (B, 3))
        self.assertEqual(tuple(data["transverse_bond_reverse"].shape),
                          (B,))
        self.assertTrue(np.array_equal(data["transverse_bond_delta_r"],
                                       topo.delta_r))
        self.assertTrue(np.array_equal(data["transverse_bond_reverse"],
                                       topo.reverse))
        self.assertTrue(np.array_equal(
            data["transverse_spatial_shape"], np.array(solver.lattice.shape)))

        # Cross-pin (anti-vacuity): the saved value equals a DIRECT
        # _run_transverse_bond_pipeline call on the same fixture state.
        beta = 1.0 / solver.T
        chi_direct, chiq_direct = solver._run_transverse_bond_pipeline(
            solver.green0, solver.green0_tail, beta, topo)
        self.assertTrue(np.array_equal(
            data["chiq_pm_bond_static"], chi_direct))
        self.assertTrue(np.array_equal(
            data["chiq_pm_static"], chiq_direct))
        self.assertTrue(np.array_equal(
            data["chiq_pm_bond_static"], gi["chiq_pm_bond_static"]))
        self.assertTrue(np.array_equal(
            data["chiq_pm_static"], gi["chiq_pm_static"]))


class TestPairHopOffsiteWarning(ApproxTestCase):
    """``_append_pairhop``'s silent off-site discard now warns, naming
    the ORIGINAL user declarations (spec "Deferred (recorded)": "Off-site
    PairHop physics (warning + tracking issue only).")."""

    def _build(self, d):
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n"
                     "0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n"
                     " 1 0 0 1 1 0.5 0.0\n-1 0 0 1 1 0.5 0.0\n")
        with open(os.path.join(d, "pairhop.dat"), "w") as f:
            # on-site (irvec=(0,0,0)) plus a Hermitian-closed off-site
            # pair (R=+-1 x, both directions declared -- #93 requires
            # Hermitian closure at the reader level).
            f.write("hdr\n1\n3\n1 1 1\n"
                     " 0 0 0 1 1 0.2 0.0\n"
                     " 1 0 0 1 1 0.15 0.0\n"
                     "-1 0 0 1 1 0.15 0.0\n")
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "PairHop": "pairhop.dat"}
        param = {"T": 1.0, "mu": 0.0, "CellShape": [4, 1, 1],
                 "SubShape": [1, 1, 1], "Nmat": 8}
        info_mode = {"mode": "RPA", "param": param,
                     "calc_scheme": "general", "calc_type": "ring"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        return rpa_mod.RPA(io.get_param("ham"), {}, info_mode)

    def test_offsite_pairhop_warns_naming_declarations(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
                self._build(d)
        msgs = "\n".join(cm.output)
        self.assertIn("PairHop", msgs)
        self.assertIn("off-site", msgs)
        # names the two dropped (pre-fold) declarations explicitly
        self.assertIn("(1, 0, 0)", msgs)
        self.assertIn("(-1, 0, 0)", msgs)
        self.assertIn("0.15", msgs)

    def _build_variant(self, d, *, include_offsite):
        """SAME fixture as ``_build`` (transfer/geometry unchanged), but
        with the off-site PairHop pair present or absent, so the two
        variants can be solved and diffed numerically."""
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n"
                     "0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n1\n2\n1 1\n"
                     " 1 0 0 1 1 0.5 0.0\n-1 0 0 1 1 0.5 0.0\n")
        with open(os.path.join(d, "pairhop.dat"), "w") as f:
            if include_offsite:
                f.write("hdr\n1\n3\n1 1 1\n"
                         " 0 0 0 1 1 0.2 0.0\n"
                         " 1 0 0 1 1 0.15 0.0\n"
                         "-1 0 0 1 1 0.15 0.0\n")
            else:
                f.write("hdr\n1\n1\n1\n 0 0 0 1 1 0.2 0.0\n")
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "PairHop": "pairhop.dat"}
        param = {"T": 1.0, "mu": 0.0, "CellShape": [4, 1, 1],
                 "SubShape": [1, 1, 1], "Nmat": 8}
        info_mode = {"mode": "RPA", "param": param,
                     "calc_scheme": "general", "calc_type": "ring"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
        green_info = io.get_param("green")
        return solver, green_info

    def test_default_path_offsite_pairhop_warns_but_output_is_unchanged(
            self):
        """Default-path regression pin (spec adjudication: the base
        solver ALREADY silently discarded off-site PairHop before this
        branch; the new warning is a log-only improvement that fires
        unconditionally, on every path -- including with the
        ``transverse_bond_channels`` gate entirely ABSENT). Proves both
        halves at once: the warning fires, AND ``solve()`` succeeds with
        numeric output BITWISE IDENTICAL to the same run with the
        off-site PairHop declaration removed entirely -- i.e. the
        discard semantics are exactly what they were before this branch,
        only the logging changed."""
        with tempfile.TemporaryDirectory() as d_off, \
                tempfile.TemporaryDirectory() as d_on:
            with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
                solver_off, gi_off = self._build_variant(
                    d_off, include_offsite=True)
            self.assertTrue(any(
                "PairHop" in m and "off-site" in m for m in cm.output))
            # gate absent (never set) == gate false
            self.assertIs(solver_off.transverse_bond_channels, False)

            solver_on, gi_on = self._build_variant(
                d_on, include_offsite=False)

            with tempfile.TemporaryDirectory(
                    prefix="rpa_pairhop_off_") as out_off, \
                    tempfile.TemporaryDirectory(
                        prefix="rpa_pairhop_on_") as out_on:
                solver_off.solve(gi_off, out_off)
                solver_on.solve(gi_on, out_on)

            self.assertTrue(np.array_equal(gi_off["chiq"], gi_on["chiq"]))

    def test_onsite_only_pairhop_does_not_warn(self):
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "geom.dat"), "w") as f:
                f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n"
                         "0.0 0.0 0.0\n")
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write("hdr\n1\n2\n1 1\n"
                         " 1 0 0 1 1 0.5 0.0\n-1 0 0 1 1 0.5 0.0\n")
            with open(os.path.join(d, "pairhop.dat"), "w") as f:
                f.write("hdr\n1\n1\n1\n 0 0 0 1 1 0.2 0.0\n")
            inter = {"path_to_input": d, "Geometry": "geom.dat",
                     "Transfer": "transfer.dat", "PairHop": "pairhop.dat"}
            param = {"T": 1.0, "mu": 0.0, "CellShape": [4, 1, 1],
                     "SubShape": [1, 1, 1], "Nmat": 8}
            info_mode = {"mode": "RPA", "param": param,
                         "calc_scheme": "general", "calc_type": "ring"}
            io = read_input_k.QLMSkInput(
                {"path_to_input": d, "interaction": inter})
            logger = logging.getLogger("hwave.solver.rpa")
            handler = _CollectingHandler()
            logger.addHandler(handler)
            try:
                rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
            finally:
                logger.removeHandler(handler)
            self.assertFalse(any(
                "PairHop declares" in r.getMessage()
                for r in handler.records))


class _CollectingHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        self.records = []

    def emit(self, record):
        self.records.append(record)


if __name__ == "__main__":
    unittest.main()
