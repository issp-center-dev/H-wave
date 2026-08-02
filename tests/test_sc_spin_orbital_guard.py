"""hwave_sc rejects [mode] enable_spin_orbital (issue #83).

Before the guard, a spin-orbital config RAN to completion and printed
eigenvalues built on four internal inconsistencies (interleaved chi0q
against spin-block consumers -- 27% off; interaction files read at the
spin-orbital dimension; paramagnetic S/C rules on spin-orbital indices;
norb never halved) -- a silent wrong result. Both public entries
(calc_eliashberg and eliashberg_dynamic.solve_dynamic) must refuse.
"""

import os
import subprocess
import sys
import tempfile
import unittest

import hwave.sc as sc


class TestSpinOrbitalGuardUnit(unittest.TestCase):

    def test_truthy_forms_are_rejected(self):
        for label, mode in (("bool", {"enable_spin_orbital": True}),
                            ("string", {"enable_spin_orbital": "true"}),
                            ("cased_key", {"Enable_Spin_Orbital": True})):
            with self.subTest(form=label):
                with self.assertRaises(ValueError) as cm:
                    sc.reject_spin_orbital_mode({"mode": mode})
                self.assertIn("enable_spin_orbital", str(cm.exception))
                self.assertIn("#83", str(cm.exception))

    def test_unrecognized_forms_are_rejected_not_read_as_false(self):
        """A user who misspells 'true' expecting SO must not silently get
        a non-SO run that then misreads their spin-orbital inputs."""
        for label, value in (("misspelled", "ture"), ("garbage", "garbage"),
                             ("int2", 2), ("list", [True]),
                             ("none_explicit", None)):
            with self.subTest(form=label):
                with self.assertRaises(ValueError):
                    sc.reject_spin_orbital_mode(
                        {"mode": {"enable_spin_orbital": value}})

    def test_recognized_false_is_normalized_for_downstream(self):
        """The resolver, not raw truthiness: a string "false" previously
        activated the SO branch of the chi0q convention check and diverged
        inside RPA (truthiness vs == True)."""
        self.assertFalse(sc._resolve_spin_orbital_flag(
            {"mode": {"enable_spin_orbital": "false"}}))
        self.assertTrue(sc._resolve_spin_orbital_flag(
            {"mode": {"Enable_Spin_Orbital": "on"}}))

    def test_loader_forwards_the_resolved_boolean(self):
        """Plumbing pin: _load_chi0q must hand validate_chi0q_index_convention
        the resolved boolean False for a string 'false', not the raw truthy
        string (which activated the SO validation branch)."""
        import numpy as np
        from unittest import mock
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        np.savez(os.path.join(d, "chi0q.npz"),
                 chi0q=np.zeros((2, 1, 1, 1), dtype=complex), momentum_convention="e_plus_ikR")
        inp = {"mode": {"enable_spin_orbital": "false",
                        "param": {"Nmat": 2}},
               "file": {"input": {"path_to_flex_output": d},
                        "output": {"path_to_output": d}}}
        seen = []
        real = sc.validate_chi0q_index_convention

        def spy(data, flag, file_name=""):
            seen.append(flag)
            return real(data, flag, file_name)

        with mock.patch.object(sc, "validate_chi0q_index_convention",
                               side_effect=spy):
            try:
                sc._load_chi0q(inp)
            except ValueError:
                pass  # shape/metadata checks past the spy may still fire
        self.assertEqual(seen, [False])

    def test_falsy_and_absent_forms_pass(self):
        for label, mode in (("absent", {}),
                            ("false", {"enable_spin_orbital": False}),
                            ("string_false",
                             {"enable_spin_orbital": "false"}),
                            ("string_off", {"enable_spin_orbital": "off"})):
            with self.subTest(form=label):
                sc.reject_spin_orbital_mode({"mode": mode})

    def test_guard_survives_python_O(self):
        """Explicit raise, not an assert (bare-assert -O defect class)."""
        code = (
            "import hwave.sc as sc\n"
            "try:\n"
            "    sc.reject_spin_orbital_mode("
            "{'mode': {'enable_spin_orbital': True}})\n"
            "except ValueError:\n"
            "    print('RAISED')\n")
        env = dict(os.environ)
        src = os.path.join(os.path.dirname(os.path.dirname(
            os.path.abspath(__file__))), "src")
        if os.path.isdir(src):
            env["PYTHONPATH"] = src + os.pathsep + env.get("PYTHONPATH", "")
        out = subprocess.run([sys.executable, "-O", "-c", code],
                             capture_output=True, text=True, env=env)
        self.assertIn("RAISED", out.stdout)


class TestSpinOrbitalGuardPublicEntries(unittest.TestCase):

    def test_static_entry_rejects_before_any_work(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n"
                    "1\n 0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n1\n2\n 1 1\n"
                    "   1    0    0    1    1   1.000000   0.0\n"
                    "  -1    0    0    1    1   1.000000   0.0\n")
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra\n1\n1\n 1\n"
                    "   0    0    0    1    1   2.000000   0.0\n")
        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        inp = {
            "mode": {"enable_spin_orbital": True,
                     "param": {"T": 2.0, "filling": 0.5,
                               "CellShape": [4, 4, 1],
                               "SubShape": [1, 1, 1], "Nmat": 8}},
            "file": {"input": {"path_to_input": d,
                               "interaction": {"path_to_input": d,
                                               "Geometry": "geom.dat",
                                               "Transfer": "transfer.dat",
                                               "CoulombIntra":
                                                   "coulombintra.dat"}},
                     "output": {"path_to_output": out}},
            "eliashberg": {"chi0q_mode": "calc", "frequency": "static",
                           "pairing_type": "singlet",
                           "solver_mode": "eigenvalue",
                           "num_eigenvalues": 1},
        }
        with self.assertRaises(ValueError) as cm:
            sc.calc_eliashberg(inp)
        self.assertIn("enable_spin_orbital", str(cm.exception))
        # rejected before any output was produced
        self.assertEqual(os.listdir(out), [])

    def test_static_entry_guard_is_the_first_operation(self):
        """A mode section alone suffices: the guard must fire before any
        required key (mode.param, files) is even read -- a guard moved
        after expensive work would raise KeyError here instead."""
        with self.assertRaises(ValueError) as cm:
            sc.calc_eliashberg({"mode": {"Enable_Spin_Orbital": True}})
        self.assertIn("enable_spin_orbital", str(cm.exception))

    def test_dynamic_entry_rejects_first(self):
        """solve_dynamic is publicly callable; the guard fires before any
        other configuration is read (a minimal dict suffices)."""
        from hwave.solver import eliashberg_dynamic
        with self.assertRaises(ValueError) as cm:
            eliashberg_dynamic.solve_dynamic(
                {"mode": {"enable_spin_orbital": True}})
        self.assertIn("enable_spin_orbital", str(cm.exception))


class TestSpinOrbitalGuardTsweep(unittest.TestCase):
    """hwave_tsweep must reject in preflight, before any FLEX rung."""

    def _base(self, so_mode, run_eli=True):
        base = {"mode": dict(so_mode,
                             param={"CellShape": [2, 2, 1], "Nmat": 8,
                                    "filling": 0.5}),
                "eliashberg": {"solver_mode": "eigenvalue"},
                "file": {"input": {}, "output": {}}}
        cont = {"run_eliashberg": run_eli}
        return base, cont

    def test_preflight_rejects_so_with_eliashberg(self):
        import hwave.tsweep as ts
        for label, so in (("bool", {"enable_spin_orbital": True}),
                          ("string", {"enable_spin_orbital": "true"}),
                          ("cased", {"Enable_Spin_Orbital": True})):
            with self.subTest(form=label):
                base, cont = self._base(so)
                with self.assertRaises(ValueError) as cm:
                    ts.preflight(base, cont)
                self.assertIn("enable_spin_orbital", str(cm.exception))

    def test_preflight_allows_flex_only_and_non_so_sweeps(self):
        import hwave.tsweep as ts
        base, cont = self._base({"enable_spin_orbital": True}, run_eli=False)
        ts.preflight(base, cont)
        base, cont = self._base({"enable_spin_orbital": "false"})
        ts.preflight(base, cont)


if __name__ == "__main__":
    unittest.main()
