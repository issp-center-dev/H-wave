# tests/test_flex_second_order_config.py
"""The flex_second_order key (spec 2026-09-08 D2, section 3): type/value,
default after auto resolution, applicability refusals, RPA warning, INFO line,
and the refusal of the retired "takimoto" value (issue #253)."""
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
        s, _ = _build({"flex_second_order": " Local "})
        self.assertEqual(s.flex_second_order, "local")
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
        s, _ = _build({"flex_second_order": "local"}, calc_scheme="auto")
        self.assertEqual(s.calc_scheme, "general")
        self.assertEqual(s.flex_second_order, "local")

    def test_rpa_warns_once_with_the_key_alone(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _build({"flex_second_order": "local"}, mode="RPA")
        hits = [m for m in cm.output if "flex_second_order" in m and "FLEX-only" in m]
        self.assertEqual(len(hits), 1)

    def test_info_line_for_general_runs(self):
        for extra in (None, {"flex_second_order": "local"}):
            with self.subTest(extra=extra):
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    _build(extra)
                self.assertEqual(
                    len([m for m in cm.output if "flex_second_order = local" in m]), 1)


class TestRetiredTakimoto(unittest.TestCase):
    """``flex_second_order = "takimoto"`` was retired (issue #253): it is
    refused at construction, in any spelling and under any scheme, with a
    message that says why and names ``"local"``."""

    def _refused(self, value, **kw):
        with self.assertRaises(ValueError) as cm:
            _build({"flex_second_order": value}, **kw)
        return str(cm.exception)

    def test_refused_in_any_case_and_names_local(self):
        for value in ("takimoto", "Takimoto", " TAKIMOTO "):
            with self.subTest(value=value):
                msg = self._refused(value)
                self.assertIn('flex_second_order = "takimoto" was retired', msg)
                self.assertIn('flex_second_order = "local"', msg)
                self.assertIn("default", msg)
                self.assertIn("not exact at second order", msg)
                self.assertIn("2.0.0", msg)

    def test_refused_before_the_scheme_applicability(self):
        """The retirement is what the user needs to hear first: it wins over
        the refusal of the key outside the general scheme, and over the
        auto resolution."""
        for scheme, inter in (("reduced", None), ("auto", None),
                              ("auto", {"CoulombIntra": "coulombintra.dat"})):
            with self.subTest(calc_scheme=scheme, interactions=inter):
                msg = self._refused("takimoto", calc_scheme=scheme, interactions=inter)
                self.assertIn("was retired", msg)
                self.assertNotIn("applies to", msg)

    def test_refused_before_the_lattice_is_built(self):
        """The key is parsed before the base constructor, so an otherwise
        broken configuration cannot mask the message."""
        from unittest import mock
        import hwave.solver.rpa as rpa_mod
        with mock.patch.object(rpa_mod, "Lattice",
                               side_effect=AssertionError("lattice built first")):
            msg = self._refused("takimoto")
        self.assertIn("was retired", msg)


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
            self.assertIsNotNone(s._second_order_factors)
            gi = r.get_param("green")
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "chi0q": "chi0q", "chiq": "chiq",
                            "sigma": "sigma", "green": "green"}, gi)
            return {f: dict(np.load(os.path.join(out, f)))
                    for f in ("chi0q.npz", "chiq_s.npz", "sigma.npz", "green.npz")}

    def test_mixed_case_general_builds_solves_and_stamps(self):
        for spelling in ("General", "GENERAL"):
            with self.subTest(spelling=spelling):
                for f, z in self._archives(spelling, "local").items():
                    self.assertEqual(str(z["flex_second_order"]), "local", f)
                    self.assertEqual(int(z["flex_second_order_schema"]), 1, f)

    def test_mixed_case_general_accepts_the_key(self):
        s, _ = _build({"flex_second_order": "Local"}, calc_scheme="General")
        self.assertEqual(s.flex_second_order, "local")

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
