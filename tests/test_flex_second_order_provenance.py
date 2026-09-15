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

    #: Every archive a general-scheme FLEX run writes when asked for all of
    #: them. ``chiq_s``/``chiq_c`` are written alongside ``chiq`` rather than
    #: requested by name.
    _ALL = ("chi0q.npz", "chiq.npz", "chiq_s.npz", "chiq_c.npz", "sigma.npz", "green.npz")
    _REQUEST = {"chi0q": "chi0q", "chiq": "chiq", "sigma": "sigma", "green": "green"}

    def test_members_in_every_general_archive(self):
        with tempfile.TemporaryDirectory() as out:
            for so in ("local", "takimoto"):
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                _solve_save(s, r, out, dict(self._REQUEST))
                for f in self._ALL:
                    z = np.load(os.path.join(out, f))
                    self.assertEqual(str(z["flex_second_order"]), so, f)
                    self.assertEqual(int(z["flex_second_order_schema"]), 1, f)
                    self.assertEqual(z["flex_second_order"].dtype.kind, "U")

    def test_absent_on_reduced_and_rpa_archives(self):
        """The members belong to the general FLEX kernel, so no reduced-scheme
        and no RPA archive may carry them.

        Each case runs in its OWN output directory and asks for EVERY archive
        it can write. An absence check that reads a directory a previous run
        has written into is only as good as the assumption that the second run
        rewrites every file the first one left -- which happens to hold here
        (a reduced FLEX solve rewrites ``chiq_s``/``chiq_c`` even though only
        ``sigma``/``green`` were requested) but is not something this test
        should depend on: a stale archive from the general run would otherwise
        be read, and the assertion would fail for a reason that has nothing to
        do with the run under test."""
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"Nmat": 8}, calc_scheme="reduced")
            _solve_save(s, r, out, dict(self._REQUEST))
            written = sorted(os.listdir(out))
            self.assertEqual(written, sorted(self._ALL))      # anti-vacuity: all present
            for f in written:
                z = np.load(os.path.join(out, f))
                self.assertNotIn("flex_second_order", z.files, f)
                self.assertNotIn("flex_second_order_schema", z.files, f)
        with tempfile.TemporaryDirectory() as out:
            s, r = _build({"Nmat": 8}, mode="RPA")
            gi = r.get_param("green")
            s.solve(gi, out)
            s.save_results(dict({"path_to_output": out}, **self._REQUEST), gi)
            written = sorted(os.listdir(out))
            self.assertTrue(written)                          # anti-vacuity
            for f in written:
                z = np.load(os.path.join(out, f))
                self.assertNotIn("flex_second_order", z.files, f)
                self.assertNotIn("flex_second_order_schema", z.files, f)

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
            self.assertTrue(any("not recorded" in m and "H-wave 2.0.0 and earlier" in m
                                for m in cm.output), cm.output)
            self.assertFalse(any("WARNING" in m and "seed computed" in m for m in cm.output))

    def test_warm_start_from_a_provenance_free_seed(self):
        """A seed with NO provenance pair -- what H-wave 2.0.0 and earlier
        wrote -- is usable under either kernel: the solve completes, the
        absence is reported at INFO (not warned about, and not refused), and
        the archives the new run writes carry the pair.

        The end-to-end leg is the point: the seed reader is unit-tested above,
        but nothing else pins that such a seed actually drives a solve to
        completion and that the run stamps its OWN kernel on what it then
        writes."""
        for so in ("local", "takimoto"):
            with self.subTest(kernel=so), tempfile.TemporaryDirectory() as out:
                s, r = _build({"flex_second_order": so, "Nmat": 8})
                _solve_save(s, r, out, {"sigma": "sigma"})
                z = dict(np.load(os.path.join(out, "sigma.npz")))
                self.assertIn("flex_second_order", z)             # anti-vacuity
                z.pop("flex_second_order"); z.pop("flex_second_order_schema")
                legacy = os.path.join(out, "legacy.npz")
                np.savez(legacy, **z)
                s2, r2 = _build({"flex_second_order": so, "Nmat": 8})
                gi = r2.get_param("green")
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    gi.update(s2.read_init({"path_to_input": out,
                                            "sigma_init": "legacy.npz"}))
                self.assertTrue(any("not recorded" in m and "H-wave 2.0.0 and earlier" in m
                                    for m in cm.output), cm.output)
                self.assertFalse(any(m.startswith("WARNING") and "seed computed" in m
                                     for m in cm.output), cm.output)
                np.testing.assert_array_equal(np.asarray(gi["sigma_init"]),
                                              np.asarray(z["sigma"]))
                with tempfile.TemporaryDirectory() as out2:
                    s2.solve(gi, out2)
                    s2.save_results(dict({"path_to_output": out2}, **self._REQUEST), gi)
                    for f in self._ALL:
                        zz = np.load(os.path.join(out2, f))
                        self.assertEqual(str(zz["flex_second_order"]), so, f)
                        self.assertEqual(int(zz["flex_second_order_schema"]), 1, f)

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
