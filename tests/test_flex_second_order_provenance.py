# tests/test_flex_second_order_provenance.py
"""Provenance members flex_second_order / flex_second_order_schema in every
general-scheme archive (spec section 3), absent on reduced archives, the
seed reader's WARNING on a kernel mismatch and INFO when absent, and the
off-site warning wording per kernel."""
import os
import tempfile
import unittest

import numpy as np

from tests.test_flex_second_order_config import _build


def _solve_save(s, r, out, files):
    gi = r.get_param("green")
    s.solve(gi, out)
    s.save_results(dict({"path_to_output": out}, **files), gi)
    return gi


class TestProvenance(unittest.TestCase):

    def test_members_in_every_general_archive_and_absent_on_reduced(self):
        with tempfile.TemporaryDirectory() as out:
            for so in ("local", "takimoto"):
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                _solve_save(s, r, out, {"chi0q": "chi0q", "chiq": "chiq", "sigma": "sigma", "green": "green"})
                for f in ("chi0q.npz", "chiq.npz", "chiq_s.npz", "chiq_c.npz", "sigma.npz", "green.npz"):
                    z = np.load(os.path.join(out, f))
                    self.assertEqual(str(z["flex_second_order"]), so, f)
                    self.assertEqual(int(z["flex_second_order_schema"]), 1, f)
                    self.assertEqual(z["flex_second_order"].dtype.kind, "U")
            s, r = _build({"Nmat": 8}, calc_scheme="reduced")
            _solve_save(s, r, out, {"sigma": "sigma", "green": "green"})
            for f in ("sigma.npz", "green.npz", "chiq_s.npz"):
                self.assertNotIn("flex_second_order", np.load(os.path.join(out, f)).files, f)

    def test_seed_warning_and_info(self):
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"flex_second_order": "takimoto", "Nmat": 8})
            _solve_save(s, r, out, {"sigma": "sigma"})
            s2, _ = _build({"flex_second_order": "local", "Nmat": 8})
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})
            self.assertTrue(any("seed computed with flex_second_order = takimoto; this run uses local" in m
                                for m in cm.output))
            z = dict(np.load(os.path.join(out, "sigma.npz")))
            z.pop("flex_second_order"); z.pop("flex_second_order_schema")
            np.savez(os.path.join(out, "legacy.npz"), **z)
            with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "legacy.npz"})
            self.assertTrue(any("not recorded" in m and "reduced scheme or pre-2.1" in m for m in cm.output))
            self.assertFalse(any("WARNING" in m and "seed computed" in m for m in cm.output))

    def test_seed_mismatch_warns_in_both_directions(self):
        with tempfile.TemporaryDirectory() as out:
            for seed_so, run_so in (("takimoto", "local"), ("local", "takimoto")):
                with self.subTest(seed=seed_so, run=run_so):
                    s, r = _build({"flex_second_order": seed_so, "Nmat": 8})
                    _solve_save(s, r, out, {"sigma": "sigma"})
                    os.replace(os.path.join(out, "sigma.npz"),
                               os.path.join(out, "seed_{}.npz".format(seed_so)))
                    s2, _ = _build({"flex_second_order": run_so, "Nmat": 8})
                    with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                        s2.read_init({"path_to_input": out,
                                      "sigma_init": "seed_{}.npz".format(seed_so)})
                    self.assertTrue(any(
                        "seed computed with flex_second_order = {}; this run uses {}".format(
                            seed_so, run_so) in m for m in cm.output), cm.output)

    def test_unreadable_seed_provenance_is_refused(self):
        """A member outside {"local", "takimoto"} or a schema this version does
        not read is a file this run cannot interpret, so it is refused by name
        and value instead of being warned about (or silently accepted)."""
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"flex_second_order": "local", "Nmat": 8})
            _solve_save(s, r, out, {"sigma": "sigma"})
            z = dict(np.load(os.path.join(out, "sigma.npz")))
            np.savez(os.path.join(out, "bad_value.npz"),
                     **dict(z, flex_second_order=np.array("exact", dtype="<U8")))
            np.savez(os.path.join(out, "bad_schema.npz"),
                     **dict(z, flex_second_order_schema=np.int64(2)))
            s2, _ = _build({"flex_second_order": "local", "Nmat": 8})
            with self.assertRaises(ValueError) as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "bad_value.npz"})
            self.assertIn("bad_value.npz", str(cm.exception))
            self.assertIn("unknown flex_second_order", str(cm.exception))
            self.assertIn("exact", str(cm.exception))
            with self.assertRaises(ValueError) as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "bad_schema.npz"})
            self.assertIn("bad_schema.npz", str(cm.exception))
            self.assertIn("unsupported flex_second_order_schema 2", str(cm.exception))
            self.assertIn("this version reads schema 1", str(cm.exception))
            # a reduced-scheme run does not read the members at all
            s3, _ = _build({"Nmat": 8}, calc_scheme="reduced")
            s3.read_init({"path_to_input": out, "sigma_init": "bad_value.npz"})
            s3.read_init({"path_to_input": out, "sigma_init": "bad_schema.npz"})

    def test_malformed_seed_provenance_is_refused(self):
        """The provenance pair must be a single canonical value each: one
        selector string exactly "local" or "takimoto", and one integer schema
        equal to 1. A member that merely converts to those (1.5, True, "1"),
        one that cannot be converted (inf), and one that is not a single value
        (a vector, an empty array) are all files this version cannot
        interpret."""
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"flex_second_order": "local", "Nmat": 8})
            _solve_save(s, r, out, {"sigma": "sigma"})
            z = dict(np.load(os.path.join(out, "sigma.npz")))
            s2, _ = _build({"flex_second_order": "local", "Nmat": 8})
            cases = {
                "schema_float": dict(flex_second_order_schema=np.float64(1.5)),
                "schema_bool": dict(flex_second_order_schema=np.bool_(True)),
                "schema_str": dict(flex_second_order_schema=np.array("1")),
                "schema_inf": dict(flex_second_order_schema=np.float64(np.inf)),
                "schema_vector": dict(flex_second_order_schema=np.array([1, 2])),
                "schema_empty": dict(flex_second_order_schema=np.array([], dtype=np.int64)),
                "selector_vector": dict(
                    flex_second_order=np.array(["local", "unknown"], dtype="<U8")),
                "selector_case": dict(flex_second_order=np.array("Local", dtype="<U8")),
                "selector_empty": dict(flex_second_order=np.array([], dtype="<U8")),
                "selector_number": dict(flex_second_order=np.int64(1)),
            }
            for name, members in cases.items():
                with self.subTest(case=name):
                    f = name + ".npz"
                    np.savez(os.path.join(out, f), **dict(z, **members))
                    with self.assertRaises(ValueError) as cm:
                        s2.read_init({"path_to_input": out, "sigma_init": f})
                    needle = ("unsupported flex_second_order_schema"
                              if name.startswith("schema") else "unknown flex_second_order")
                    self.assertIn(needle, str(cm.exception))
                    self.assertIn(f, str(cm.exception))
            # the canonical pair is still accepted
            s2.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})

    def test_seed_check_runs_under_a_mixed_case_scheme(self):
        """read_init's seed-provenance branch keys off calc_scheme as well, so
        a ``calc_scheme = "General"`` run must still compare the seed's kernel
        with this run's."""
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"flex_second_order": "takimoto", "Nmat": 8})
            _solve_save(s, r, out, {"sigma": "sigma"})
            s2, _ = _build({"flex_second_order": "local", "Nmat": 8}, calc_scheme="General")
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s2.read_init({"path_to_input": out, "sigma_init": "sigma.npz"})
            self.assertTrue(any("seed computed with flex_second_order = takimoto" in m
                                for m in cm.output), cm.output)

    def test_members_in_the_dedicated_bond_archive(self):
        """The bond archive is written by its own branch of save_results, so it
        carries the provenance pair separately from sigma/green."""
        from tests.test_flex_bond_gate import _flex
        with tempfile.TemporaryDirectory() as out:
            s, r = _flex({"longitudinal_bond_output_full": True, "IterationMax": 1,
                          "flex_second_order": "local"}, gate=True)
            self.assertEqual(s.flex_second_order, "local")
            gi = r.get_param("green")
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "sigma": "sigma", "green": "green"}, gi)
            for f in ("longitudinal_bond.npz", "sigma.npz", "green.npz"):
                z = np.load(os.path.join(out, f))
                self.assertEqual(str(z["flex_second_order"]), "local", f)
                self.assertEqual(int(z["flex_second_order_schema"]), 1, f)

    def test_members_in_an_ir_basis_general_run(self):
        """The IR-native general path writes its own sigma/green blocks."""
        try:
            import sparse_ir           # noqa: F401
        except ImportError:
            self.skipTest("sparse-ir not installed")
        from tests.test_flex_ir_general import _make_general_solver
        with tempfile.TemporaryDirectory() as out:
            s, gi = _make_general_solver(64, matsubara_basis="ir", iteration_max=1,
                                         extra_param={"flex_second_order": "local"})
            self.assertTrue(s.use_ir)
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "sigma": "sigma", "green": "green"}, gi)
            for f in ("sigma.npz", "green.npz"):
                z = np.load(os.path.join(out, f))
                self.assertEqual(str(z["flex_second_order"]), "local", f)
                self.assertEqual(int(z["flex_second_order_schema"]), 1, f)

    def test_offsite_warning_wording(self):
        for so, needle in (("local", "exact at second order"), ("takimoto", "is omitted (the same approximation")):
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                s._calc_epsilon_k({})
                G = s._calc_dressed_green(0.5, 0.1, np.zeros((1, 8, 16, 2, 2), complex))
                s._flex_compute_veff_general(s._calc_chi0q(G, np.zeros_like(G), 0.5)[0], s.ham_info.ham_inter_q)
            self.assertTrue(any(needle in m for m in cm.output), (so, cm.output))


if __name__ == "__main__":
    unittest.main()
