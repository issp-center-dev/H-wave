"""The bond topology is channel-agnostic (GitHub issue #181, Tier 3
Phase A): ``BondTopology`` carries the channel set, ``resolve_bond_topology``
resolves it over an explicit ``active_types`` selection, and the Phase W
names (``TransverseTopology``, ``resolve_transverse_topology``,
``transverse_effective_activity``) are thin wrappers/aliases whose
behaviour and diagnostics are unchanged. ``BondSetView`` presents a
topology to the bubble kernel's duck-typed bond-set contract.
"""
import unittest

import numpy as np


def _decl(t, R, a, b, v):
    return {t: {((R, 0, 0), (a, b)): v, ((-R, 0, 0), (b, a)): np.conj(v)}}


class TestResolveBondTopology(unittest.TestCase):

    def test_longitudinal_types_and_order(self):
        from hwave.solver import bond_channels as bc
        inter = _decl("Hund", 1, 0, 0, 0.3)
        topo = bc.resolve_bond_topology(
            inter, np.eye(3), 1, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        self.assertEqual(tuple(topo.coeffs), ("CoulombInter", "Hund", "Ising"))
        self.assertEqual(topo.delta_r.shape, (3, 3))          # 0, +1, -1
        m = [tuple(r) for r in topo.delta_r].index((1, 0, 0))
        self.assertEqual(complex(topo.coeffs["Hund"][m, 0, 0]), 0.3 + 0j)
        self.assertTrue(np.all(topo.coeffs["CoulombInter"] == 0))
        self.assertIsInstance(topo, bc.BondTopology)

    def test_active_types_validated(self):
        from hwave.solver import bond_channels as bc
        for bad in ((), ("CoulombInter", "CoulombInter"), ("Coulomb",),
                    ("PairHop",), ("coulombinter",), "CoulombInter"):
            with self.subTest(active_types=bad):
                with self.assertRaises(ValueError) as cm:
                    bc.resolve_bond_topology({}, np.eye(3), 1, active_types=bad)
                self.assertIn("active_types", str(cm.exception))

    def test_resolvable_type_sets(self):
        from hwave.solver import bond_channels as bc
        self.assertEqual(bc._BOND_RESOLVABLE_TYPES,
                         ("CoulombInter", "Hund", "Ising", "Exchange"))
        self.assertEqual(bc._LONGITUDINAL_ACTIVE_TYPES,
                         ("CoulombInter", "Hund", "Ising"))
        self.assertEqual(bc._TRANSVERSE_ACTIVE_TYPES,
                         ("CoulombInter", "Ising", "Exchange"))

    def test_hund_is_a_shell_source_only_when_active(self):
        """Phase W: Hund contributes no transverse shell. Longitudinal: it
        does. One resolver, two selections."""
        from hwave.solver import bond_channels as bc
        inter = _decl("Hund", 1, 0, 1, 0.2)
        t = bc.resolve_transverse_topology(inter, np.eye(3), 2)
        self.assertEqual(t.delta_r.shape, (1, 3))
        l = bc.resolve_bond_topology(
            inter, np.eye(3), 2, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        self.assertEqual(l.delta_r.shape, (3, 3))

    def test_transverse_wrapper_is_unchanged(self):
        from hwave.solver import bond_channels as bc
        inter = _decl("Exchange", 1, 0, 1, 0.2)
        a = bc.resolve_transverse_topology(inter, np.eye(3), 2)
        b = bc.resolve_bond_topology(
            inter, np.eye(3), 2, active_types=bc._TRANSVERSE_ACTIVE_TYPES)
        self.assertEqual(tuple(a.coeffs), ("CoulombInter", "Ising", "Exchange"))
        for t in a.coeffs:
            self.assertTrue(np.array_equal(a.coeffs[t], b.coeffs[t]))
        self.assertTrue(np.array_equal(a.delta_r, b.delta_r))
        self.assertIs(bc.TransverseTopology, bc.BondTopology)
        with self.assertRaises(ValueError) as cm:
            bc.resolve_transverse_topology(inter, np.eye(3), 0)
        self.assertTrue(str(cm.exception).startswith(
            "resolve_transverse_topology:"), str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            bc.resolve_bond_topology(
                inter, np.eye(3), 0, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        self.assertTrue(str(cm.exception).startswith(
            "resolve_bond_topology:"), str(cm.exception))

    def test_max_shells_guard_uses_the_selected_types(self):
        """A truncation that drops a declared nonzero shell of a SELECTED
        type raises; the same shell of an unselected type is invisible."""
        from hwave.solver import bond_channels as bc
        inter = _decl("Hund", 2, 0, 0, 0.5)
        inter.update(_decl("CoulombInter", 1, 0, 0, 0.1))
        with self.assertRaises(ValueError):
            bc.resolve_bond_topology(
                inter, np.eye(3), 1, max_shells=1,
                active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        t = bc.resolve_transverse_topology(inter, np.eye(3), 1, max_shells=1)
        self.assertEqual(t.delta_r.shape, (3, 3))


class TestBondSetView(unittest.TestCase):

    def test_duck_types_the_bubble_contract(self):
        from hwave.solver import bond_channels as bc
        from hwave.solver import bubble
        topo = bc.resolve_bond_topology(
            _decl("CoulombInter", 1, 0, 0, 0.1), np.eye(3), 1,
            active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        self.assertEqual(view.n_channels, 3)
        self.assertEqual(view.delta_r, ((0, 0, 0), (-1, 0, 0), (1, 0, 0)))
        self.assertEqual(view.reverse, (0, 2, 1))
        self.assertIsInstance(view.n_channels, int)
        self.assertTrue(all(type(c) is int for r in view.delta_r for c in r))
        n, dr = bubble._validate_bond_set(view)       # the kernel's own guard
        self.assertEqual(n, 3)
        self.assertEqual(tuple(tuple(x) for x in dr), view.delta_r)

    def test_rejects_a_broken_topology(self):
        from hwave.solver import bond_channels as bc

        class Broken:
            delta_r = np.array([[0, 0, 0], [1, 0, 0]])
            reverse = np.array([0, 1])          # not geometrically consistent
            coeffs = {"CoulombInter": np.zeros((2, 1, 1), complex)}
        with self.assertRaises(ValueError):
            bc.BondSetView(Broken())


class TestActivityAlias(unittest.TestCase):

    def test_alias_and_zero_topology(self):
        from hwave.solver import bond_channels as bc
        topo = bc.resolve_bond_topology(
            _decl("Ising", 1, 0, 0, 0.0), np.eye(3), 1,
            active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        self.assertFalse(bc.bond_effective_activity(topo))
        self.assertIs(bc.transverse_effective_activity,
                      bc.bond_effective_activity)


if __name__ == "__main__":
    unittest.main()
