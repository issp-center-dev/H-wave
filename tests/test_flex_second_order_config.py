# tests/test_flex_second_order_config.py
"""The flex_second_order key (spec 2026-09-08 D2, section 3): type/value,
default after auto resolution, applicability refusals, RPA warning, INFO line."""
import logging
import unittest

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"


def _build(param_extra=None, calc_scheme="general", interactions=None, mode="FLEX"):
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
    par.update(param_extra or {})
    info = {"mode": mode, "param": par, "enable_spin_orbital": False, "calc_scheme": calc_scheme}
    if mode == "FLEX":
        import hwave.solver.flex as flex_mod
        return flex_mod.FLEX(r.get_param("ham"), {}, info), r
    import hwave.solver.rpa as rpa_mod
    info["calc_type"] = "ring"
    return rpa_mod.RPA(r.get_param("ham"), {}, info), r


class TestKey(unittest.TestCase):

    def test_default_is_local_and_canonical(self):
        s, _ = _build()
        self.assertEqual(s.flex_second_order, "local")
        s, _ = _build({"flex_second_order": "Takimoto"})
        self.assertEqual(s.flex_second_order, "takimoto")
        s, _ = _build({"FLEX_SECOND_ORDER": "LOCAL"})
        self.assertEqual(s.flex_second_order, "local")

    def test_type_and_value_refused_before_applicability(self):
        for bad in (1, True, "exact", "", None):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError) as cm:
                    _build({"flex_second_order": bad}, calc_scheme="reduced")
                self.assertIn("flex_second_order", str(cm.exception))
                self.assertNotIn("applies to", str(cm.exception))

    def test_reduced_refusals(self):
        # explicit reduced + explicit key
        with self.assertRaises(ValueError) as cm:
            _build({"flex_second_order": "local"}, calc_scheme="reduced")
        self.assertIn('applies to calc_scheme = "general" only', str(cm.exception))
        # auto -> reduced (on-site Hubbard only) + explicit key: names the resolution
        with self.assertRaises(ValueError) as cm:
            _build({"flex_second_order": "local"}, calc_scheme="auto",
                   interactions={"CoulombIntra": "coulombintra.dat"})
        self.assertIn('resolved to "reduced"', str(cm.exception))
        # absent key: silent on both
        s, _ = _build(calc_scheme="reduced")
        self.assertEqual(s.flex_second_order, "local")     # attribute exists, inert
        s, _ = _build(calc_scheme="auto", interactions={"CoulombIntra": "coulombintra.dat"})
        self.assertEqual(s.calc_scheme, "reduced")

    def test_auto_to_general_accepts_the_key(self):
        s, _ = _build({"flex_second_order": "takimoto"}, calc_scheme="auto")
        self.assertEqual(s.calc_scheme, "general")
        self.assertEqual(s.flex_second_order, "takimoto")

    def test_rpa_warns_once_with_the_key_alone(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _build({"flex_second_order": "local"}, mode="RPA")
        hits = [m for m in cm.output if "flex_second_order" in m and "FLEX-only" in m]
        self.assertEqual(len(hits), 1)

    def test_info_line_for_general_runs(self):
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            _build()
        self.assertTrue(any("flex_second_order = local" in m for m in cm.output))
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            _build({"flex_second_order": "takimoto"})
        self.assertTrue(any("flex_second_order = takimoto" in m for m in cm.output))


if __name__ == "__main__":
    unittest.main()
