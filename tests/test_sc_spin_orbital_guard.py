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

    def test_dynamic_entry_rejects_first(self):
        """solve_dynamic is publicly callable; the guard fires before any
        other configuration is read (a minimal dict suffices)."""
        from hwave.solver import eliashberg_dynamic
        with self.assertRaises(ValueError) as cm:
            eliashberg_dynamic.solve_dynamic(
                {"mode": {"enable_spin_orbital": True}})
        self.assertIn("enable_spin_orbital", str(cm.exception))


if __name__ == "__main__":
    unittest.main()
