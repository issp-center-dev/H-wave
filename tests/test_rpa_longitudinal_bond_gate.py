"""The experimental bond-resolved longitudinal channel of the RPA solver
(``[mode.param] longitudinal_bond_channels = true``; GitHub issue #181,
Tier 3 Phase A): configuration, prerequisites, memory preflight, wiring,
and the ``longitudinal_bond_*`` outputs (spec gates G0 and G4).
"""
import logging
import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as rpa_mod

_IN = "tests/rpa/input_2orb"


def _build(param_extra=None, interactions=None, calc_type="ring",
           mode="RPA", nmat=32, path=_IN, cell=(4, 4, 1)):
    idict = {"path_to_input": path, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    idict.update(interactions or {"CoulombInter": "coulombinter.dat"})
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": list(cell),
           "SubShape": [1, 1, 1], "Nmat": nmat}
    par.update(param_extra or {})
    info = {"mode": mode, "param": par, "enable_spin_orbital": False,
            "calc_scheme": "general", "calc_type": calc_type}
    if mode == "FLEX":
        import hwave.solver.flex as flex_mod
        par.update({"IterationMax": 1, "Mix": 1.0, "EPS": 1})
        del info["calc_type"]
        return flex_mod.FLEX(r.get_param("ham"), {}, info), r
    return rpa_mod.RPA(r.get_param("ham"), {}, info), r


class TestGateConfig(unittest.TestCase):

    def test_defaults_off(self):
        s, _ = _build()
        self.assertFalse(s.longitudinal_bond_channels)
        self.assertIsNone(s.longitudinal_bond_max_shells)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb,
                         rpa_mod.LONGITUDINAL_BOND_MEMORY_CAP_GB_DEFAULT)
        self.assertEqual(rpa_mod.LONGITUDINAL_BOND_MEMORY_CAP_GB_DEFAULT, 8.0)

    def test_on_and_options(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_max_shells": 2,
                       "longitudinal_bond_memory_cap_gb": 1.5})
        self.assertTrue(s.longitudinal_bond_channels)
        self.assertEqual(s.longitudinal_bond_max_shells, 2)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb, 1.5)
        s, _ = _build({"Longitudinal_Bond_Channels": True})   # case-insensitive
        self.assertTrue(s.longitudinal_bond_channels)

    def test_bounds_and_types(self):
        for bad in ({"longitudinal_bond_channels": "yes"},
                    {"longitudinal_bond_channels": 1},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": 0},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": 1.5},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_max_shells": True},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": 0.0},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": float("nan")},
                    {"longitudinal_bond_channels": True,
                     "longitudinal_bond_memory_cap_gb": "8"}):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError) as cm:
                    _build(bad)
                self.assertIn("longitudinal_bond", str(cm.exception))

    def test_dependent_keys_warn_when_off_or_wrong_calc_type(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _build({"longitudinal_bond_max_shells": 1})
        self.assertTrue(any("longitudinal_bond_max_shells" in m
                            for m in cm.output), cm.output)
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            s, _ = _build({"longitudinal_bond_channels": True},
                          calc_type="ring+ladder")
        self.assertFalse(s.longitudinal_bond_channels)
        self.assertTrue(any("calc_type" in m for m in cm.output), cm.output)

    def test_both_gates_true_is_rejected_at_construction(self):
        for ct in ("ring", "ring+ladder"):
            with self.subTest(calc_type=ct):
                with self.assertRaises(ValueError) as cm:
                    _build({"longitudinal_bond_channels": True,
                            "transverse_bond_channels": True}, calc_type=ct)
                self.assertIn("calc_type", str(cm.exception))
        # an explicit false next to a true is fine
        s, _ = _build({"longitudinal_bond_channels": True,
                       "transverse_bond_channels": False})
        self.assertTrue(s.longitudinal_bond_channels)

    def test_info_line_logged(self):
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            _build({"longitudinal_bond_channels": True})
        self.assertTrue(any("longitudinal_bond_channels = True" in m
                            for m in cm.output))

    def test_flex_rejects_only_a_true_key(self):
        with self.assertRaises(ValueError) as cm:
            _build({"longitudinal_bond_channels": True}, mode="FLEX")
        self.assertIn("Phase B", str(cm.exception))
        _build({"longitudinal_bond_channels": False}, mode="FLEX")
        with self.assertRaises(ValueError):
            _build({"longitudinal_bond_channels": "true"}, mode="FLEX")


_EQ2 = "tests/equivalence_input/orb2"


class TestGatePrerequisites(unittest.TestCase):

    def _gate(self, param_extra=None, interactions=None, path=_IN):
        s, _ = _build(dict({"longitudinal_bond_channels": True},
                           **(param_extra or {})), interactions, path=path)
        s.spin_mode = "spin-free"          # what solve() sets before the check
        return s

    def test_refusals(self):
        cases = [
            ({}, {"CoulombInter": "onsite_inter.dat"}, _IN, "declared off-site"),
            ({}, {"PairHop": "offsite_pairhop.dat"}, _EQ2, "PairHop"),
            ({}, {"Exchange": "offsite_exchange.dat"}, _EQ2, "Exchange"),
            ({}, {"CoulombInter": "offsite_coulombinter_interorb.dat",
                  "Exchange": "offsite_exchange.dat"}, _EQ2, "Exchange"),
            ({"SubShape": [2, 1, 1]}, None, _IN, "sublattice"),
        ]
        for extra, inter, path, fragment in cases:
            with self.subTest(fragment=fragment):
                s = self._gate(extra, inter, path)
                with self.assertRaises(ValueError) as cm:
                    s._validate_longitudinal_bond_prereqs()
                self.assertIn("longitudinal_bond_channels", str(cm.exception))
                self.assertIn(fragment, str(cm.exception))

    def test_odd_nmat_is_the_constructor_exit(self):
        """The documented exception to 'every refusal is a ValueError':
        an odd Nmat is refused by the constructor (sys.exit) before the
        gate is ever reached; the gate's own even-Nmat guard is a
        belt-and-braces check for programmatic callers."""
        with self.assertRaises(SystemExit):
            _build({"longitudinal_bond_channels": True, "Nmat": 31})
        s = self._gate()
        s.nmat = 31
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("even", str(cm.exception))

    def test_spin_mode_and_external_chi0q_refused(self):
        s = self._gate()
        s.spin_mode = "spinful"
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("spin", str(cm.exception))
        s = self._gate()
        s._chi0q_external = True
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("chi0q", str(cm.exception))

    def test_accepts_and_returns_topology_and_split(self):
        s = self._gate()
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(tuple(topo.coeffs), ("CoulombInter", "Hund", "Ising"))
        self.assertGreater(topo.delta_r.shape[0], 1)
        self.assertEqual(split.offsite_types, ("CoulombInter",))
        self.assertFalse(split.has_fold)

    def test_declared_but_zero_shell_is_active(self):
        s = self._gate()
        topo, _ = s._validate_longitudinal_bond_prereqs()
        s2 = self._gate()
        for k in list(s2.ham_info.param_ham["CoulombInter"]):
            s2.ham_info.param_ham["CoulombInter"][k] = 0.0
        topo2, _ = s2._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo2.delta_r.shape, topo.delta_r.shape)
        self.assertTrue(all(np.all(a == 0) for a in topo2.coeffs.values()))

    def test_hund_and_ising_shells_count_and_pairlift_is_inert(self):
        s = self._gate({}, {"Hund": "offsite_hund.dat"})
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo.delta_r.shape[0], 3)
        self.assertEqual(split.offsite_types, ("Hund",))
        s = self._gate({}, {"Hund": "offsite_hund.dat"})
        s.ham_info.param_ham["PairLift"] = {((1, 0, 0), (0, 1)): 0.1,
                                            ((-1, 0, 0), (1, 0)): 0.1}
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertTrue(any("PairLift" in m for m in cm.output))
        self.assertEqual(topo.delta_r.shape[0], 3)

    def test_max_shells_truncation_of_a_declared_nonzero_shell_refused(self):
        s = self._gate({"longitudinal_bond_max_shells": 1})
        s.ham_info.param_ham["CoulombInter"][((2, 0, 0), (0, 0))] = 0.05
        s.ham_info.param_ham["CoulombInter"][((-2, 0, 0), (0, 0))] = 0.05
        with self.assertRaises(ValueError):
            s._validate_longitudinal_bond_prereqs()
        s.ham_info.param_ham["CoulombInter"][((2, 0, 0), (0, 0))] = 0.0
        s.ham_info.param_ham["CoulombInter"][((-2, 0, 0), (0, 0))] = 0.0
        topo, _ = s._validate_longitudinal_bond_prereqs()   # zero shell dropped
        self.assertEqual(topo.delta_r.shape[0], 5)          # 4x4: +-x, +-y

    def test_complex_coefficient_refused(self):
        s = self._gate()
        s.ham_info.param_ham["CoulombInter"][((1, 0, 0), (0, 1))] = 0.1 + 0.2j
        s.ham_info.param_ham["CoulombInter"][((-1, 0, 0), (1, 0))] = 0.1 - 0.2j
        with self.assertRaises(ValueError) as cm:
            s._validate_longitudinal_bond_prereqs()
        self.assertIn("complex", str(cm.exception).lower())

    def test_aggregate_coulomb_feeds_the_topology(self):
        s = self._gate({}, {"CoulombInter": "onsite_inter.dat"})
        del s.ham_info.param_ham["CoulombInter"]
        s.ham_info.param_ham["Coulomb"] = {((0, 0, 0), (0, 0)): 1.0,
                                           ((1, 0, 0), (0, 0)): 0.2,
                                           ((-1, 0, 0), (0, 0)): 0.2}
        topo, split = s._validate_longitudinal_bond_prereqs()
        self.assertEqual(topo.delta_r.shape[0], 3)
        self.assertEqual(split.offsite_types, ("CoulombInter",))


class TestGatePreflight(unittest.TestCase):

    def test_estimate_families_and_cap(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e6})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        est = s._longitudinal_bond_resource_preflight(topo)
        nvol, P, nmat = s.lattice.nvol, s.norb ** 2, s.nmat
        ND = topo.delta_r.shape[0] * P
        U = nvol * ND * ND * 16
        self.assertEqual(est["ND"], ND)
        self.assertAlmostEqual(est["peak_solve"], 1.25 * (5 + 2) * U)
        prep = 16 * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = 16 * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        self.assertAlmostEqual(est["bubble_bytes"], 1.25 * max(prep, pair))
        self.assertAlmostEqual(est["peak"], max(est["bubble_bytes"], est["peak_solve"]))
        self.assertEqual(est["cap_bytes"], 1e6 * 1024 ** 3)

    def test_cap_refusal_names_the_shapes(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e-9})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        with self.assertRaises(ValueError) as cm:
            s._longitudinal_bond_resource_preflight(topo)
        msg = str(cm.exception)
        self.assertIn("longitudinal_bond_memory_cap_gb", msg)
        self.assertIn("ND", msg)
        self.assertIn("GiB", msg)

    def test_op_count_warning(self):
        s, _ = _build({"longitudinal_bond_channels": True,
                       "longitudinal_bond_memory_cap_gb": 1e9})
        s.spin_mode = "spin-free"
        topo, _ = s._validate_longitudinal_bond_prereqs()
        s.lattice.nvol = 10 ** 9                    # synthetic: force the warning
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            s._longitudinal_bond_resource_preflight(topo)
        self.assertTrue(any("SVD" in m for m in cm.output), cm.output)


if __name__ == "__main__":
    unittest.main()
