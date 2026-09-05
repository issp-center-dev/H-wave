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


if __name__ == "__main__":
    unittest.main()
