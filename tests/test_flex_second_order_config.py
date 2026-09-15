# tests/test_flex_second_order_config.py
"""The flex_second_order key (spec 2026-09-08 D2, section 3): type/value,
default after auto resolution, applicability refusals, RPA warning, INFO line."""
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"


def _build(param_extra=None, calc_scheme="general", interactions=None, mode="FLEX", path=_IN2):
    idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
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
        for so in ("local", "takimoto"):
            with self.subTest(second_order=so):
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    _build(None if so == "local" else {"flex_second_order": so})
                self.assertEqual(
                    len([m for m in cm.output if "flex_second_order = {}".format(so) in m]), 1)


class TestSchemeSpelling(unittest.TestCase):
    """calc_scheme is compared case-insensitively everywhere the key acts:
    RPA stores the raw request and FLEX's own dispatch normalises it, so a
    ``calc_scheme = "General"`` run IS a general-scheme run and must build the
    factors, log the INFO line, solve and stamp the provenance exactly like
    the lower-case spelling."""

    def _archives(self, spelling, so):
        with tempfile.TemporaryDirectory() as out:
            with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                s, r = _build({"flex_second_order": so, "Nmat": 8}, calc_scheme=spelling)
            self.assertEqual(
                len([m for m in cm.output if "flex_second_order = {}".format(so) in m]), 1)
            if so == "local":
                self.assertIsNotNone(s._second_order_factors)
            else:
                self.assertIsNone(s._second_order_factors)
            gi = r.get_param("green")
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "chi0q": "chi0q", "chiq": "chiq",
                            "sigma": "sigma", "green": "green"}, gi)
            return {f: dict(np.load(os.path.join(out, f)))
                    for f in ("chi0q.npz", "chiq_s.npz", "sigma.npz", "green.npz")}

    def test_mixed_case_general_builds_solves_and_stamps(self):
        for spelling in ("General", "GENERAL"):
            for so in ("local", "takimoto"):
                with self.subTest(spelling=spelling, second_order=so):
                    for f, z in self._archives(spelling, so).items():
                        self.assertEqual(str(z["flex_second_order"]), so, f)
                        self.assertEqual(int(z["flex_second_order_schema"]), 1, f)

    def test_mixed_case_general_accepts_the_key(self):
        s, _ = _build({"flex_second_order": "takimoto"}, calc_scheme="General")
        self.assertEqual(s.flex_second_order, "takimoto")

    def test_mixed_case_reduced_still_refuses_exchange(self):
        """The inherited RPA validations key off the scheme name too: an
        Exchange (or PairHop) declaration under a mis-cased "Reduced" must hit
        the same refusal as under "reduced" -- the reduced scheme drops those
        vertices entirely, so accepting the request would run a silently
        different Hamiltonian."""
        for itype, fname in (("Exchange", "exchange.dat"), ("PairHop", "pairhop.dat")):
            with tempfile.TemporaryDirectory() as d:
                for f in ("geom.dat", "transfer.dat", "coulombintra.dat"):
                    shutil.copy(os.path.join(_IN2, f), d)
                with open(os.path.join(d, fname), "w") as fw:
                    fw.write("{} in wannier90-like format for uhfk\n2\n1\n 1\n"
                             "   0    0    0    1    2  0.3 0.0\n"
                             "   0    0    0    2    1  0.3 0.0\n".format(itype))
                inter = {itype: fname, "CoulombIntra": "coulombintra.dat"}
                msgs = {}
                for spelling in ("reduced", "Reduced", "REDUCED"):
                    with self.subTest(itype=itype, spelling=spelling):
                        with self.assertRaises(ValueError) as cm:
                            _build(calc_scheme=spelling, interactions=inter, path=d)
                        msgs[spelling] = str(cm.exception)
                        self.assertIn("no density-diagonal content", msgs[spelling])
                        self.assertIn(itype, msgs[spelling])

    def test_mixed_case_auto_resolves_like_auto(self):
        for spelling in ("auto", "Auto", "AUTO"):
            with self.subTest(spelling=spelling):
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    s, _ = _build(calc_scheme=spelling)
                self.assertEqual(s.calc_scheme, "general")
                self.assertEqual(s.calc_scheme_requested, "auto")
                self.assertTrue(s._scheme_resolution.startswith("auto:"), s._scheme_resolution)
                self.assertEqual(
                    len([m for m in cm.output if "flex_second_order = local" in m]), 1)

    def test_mixed_case_reduced_still_refuses_the_key(self):
        for spelling in ("Reduced", "REDUCED"):
            with self.subTest(spelling=spelling):
                with self.assertRaises(ValueError) as cm:
                    _build({"flex_second_order": "local"}, calc_scheme=spelling)
                self.assertIn('applies to calc_scheme = "general" only', str(cm.exception))
        # auto -> reduced keeps its own message (the resolver writes the
        # canonical lower-case name, so only the explicit spelling varies)
        with self.assertRaises(ValueError) as cm:
            _build({"flex_second_order": "local"}, calc_scheme="auto",
                   interactions={"CoulombIntra": "coulombintra.dat"})
        self.assertIn('resolved to "reduced"', str(cm.exception))


if __name__ == "__main__":
    unittest.main()
