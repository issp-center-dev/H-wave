"""FLEX artifacts written before the bond-placement correction of issue #232
put an off-site two-body row on the bond opposite its Transfer row. Every
artifact a later run consumes carries ``sc_vertex_version`` (3 = the
corrected placement), and a consumer refuses an older artifact exactly when
the declaration is placement-sensitive -- an off-site table that is not
symmetric under ``R -> -R`` on the run's grid. These tests pin the shared
check, the ``sigma_init`` seed gate and the stamps on ``sigma.npz`` /
``green.npz``; the susceptibility and bond-archive gates are pinned in
``tests/test_sc_vertex_adjudicated.py`` and
``tests/test_eliashberg_bond_entry.py``."""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.sc import check_placement_stamp, placement_sensitive_types

_ASYM = {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4}
_SYM = {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
        ((-1, 0, 0), (0, 1)): 0.4, ((1, 0, 0), (1, 0)): 0.4}
_IN2 = "tests/rpa/input_2orb"


class TestPlacementCheck(unittest.TestCase):

    def test_sensitive_types(self):
        self.assertEqual(placement_sensitive_types({"CoulombInter": _ASYM}), ["CoulombInter"])
        self.assertEqual(placement_sensitive_types({"CoulombInter": _SYM}), [])
        # the aggregate Coulomb table carries CoulombInter rows
        self.assertEqual(placement_sensitive_types({"Coulomb": _ASYM}), ["Coulomb"])
        # PairLift has no particle-hole content
        self.assertEqual(placement_sensitive_types({"PairLift": _ASYM}), [])
        # on a 2-site axis +x and -x are the same displacement
        self.assertEqual(placement_sensitive_types({"CoulombInter": _ASYM}, grid=(2, 1, 1)), [])

    def test_stamp_decides_for_a_sensitive_declaration(self):
        decl = {"CoulombInter": _ASYM}
        check_placement_stamp({"sc_vertex_version": np.int64(3)}, "f", decl)
        for data in ({}, {"sc_vertex_version": np.int64(2)}):
            with self.subTest(data=data):
                with self.assertRaisesRegex(ValueError, "#232"):
                    check_placement_stamp(data, "f", decl)

    def test_insensitive_declaration_accepts_a_legacy_file(self):
        check_placement_stamp({}, "f", {"CoulombInter": _SYM})
        check_placement_stamp({"sc_vertex_version": np.int64(2)}, "f", {"CoulombIntra": {}})

    def test_malformed_stamp_is_refused(self):
        for bad in (np.float64(2.5), np.str_("3"), np.bool_(True), np.array([3, 3])):
            with self.subTest(stamp=bad):
                with self.assertRaisesRegex(ValueError, "malformed"):
                    check_placement_stamp({"sc_vertex_version": bad}, "f", {"CoulombInter": _SYM})


def _flex_solver(interaction_file, cell=(4, 4, 1)):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    inter = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombInter": interaction_file}
    io = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": inter})
    mode = {"mode": "FLEX", "calc_scheme": "general",
            "param": {"T": 2.0, "filling": 0.5, "CellShape": list(cell),
                      "SubShape": [1, 1, 1], "Nmat": 32, "IterationMax": 1}}
    return solver_flex.FLEX(io.get_param("ham"), {}, mode), io.get_param("green")


class TestSigmaSeedGate(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _seed(self, stamp):
        members = dict(sigma=np.zeros((1, 32, 16, 4, 4), dtype=complex),
                       cell_shape=np.array([4, 4, 1]), momentum_convention="e_plus_ikR")
        if stamp is not None:
            members["sc_vertex_version"] = np.int64(stamp)
        np.savez(os.path.join(self.tmp, "sigma.npz"), **members)
        return {"path_to_input": self.tmp, "sigma_init": "sigma.npz"}

    def test_legacy_seed_is_refused_for_a_sensitive_declaration(self):
        # tests/rpa/input_2orb/coulombinter.dat declares the inter-orbital
        # rows v_12(-x) / v_21(+x) only: not symmetric under R -> -R
        solver, _ = _flex_solver("coulombinter.dat")
        for stamp in (None, 2):
            with self.subTest(stamp=stamp):
                with self.assertRaisesRegex(ValueError, "#232"):
                    solver.read_init(self._seed(stamp))
        info = solver.read_init(self._seed(3))
        self.assertIn("sigma_init", info)

    def test_legacy_seed_is_accepted_for_an_insensitive_declaration(self):
        solver, _ = _flex_solver("onsite_inter.dat")
        info = solver.read_init(self._seed(None))
        self.assertIn("sigma_init", info)


class TestFlexWritesTheStamp(unittest.TestCase):

    def test_sigma_and_green_carry_version_3(self):
        solver, green = _flex_solver("coulombinter.dat")
        out = tempfile.mkdtemp()
        try:
            solver.solve(green, out)
            solver.save_results({"path_to_output": out, "sigma": "sigma.npz",
                                 "green": "green.npz"}, green)
            for name in ("sigma.npz", "green.npz"):
                with np.load(os.path.join(out, name)) as data:
                    self.assertEqual(int(data["sc_vertex_version"]), 3, name)
        finally:
            shutil.rmtree(out, ignore_errors=True)


class TestGreenLoaderGate(unittest.TestCase):
    """``hwave.sc._load_flex_green`` -- the loader every Eliashberg path
    (static, dynamic, dynamic bond) goes through -- on a real FLEX
    ``green.npz``, stamped and with the stamp removed."""

    @classmethod
    def setUpClass(cls):
        cls.out = tempfile.mkdtemp()
        solver, green = _flex_solver("coulombinter.dat")
        solver.solve(green, cls.out)
        solver.save_results({"path_to_output": cls.out, "green": "green.npz"}, green)
        with np.load(os.path.join(cls.out, "green.npz")) as data:
            members = {k: data[k] for k in data.files}
        members.pop("sc_vertex_version")
        cls.legacy = tempfile.mkdtemp()
        np.savez(os.path.join(cls.legacy, "green.npz"), **members)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.out, ignore_errors=True)
        shutil.rmtree(cls.legacy, ignore_errors=True)

    @staticmethod
    def _load(directory, interactions, cell=(4, 4, 1)):
        import hwave.sc as sc
        inp = {"mode": {"param": {"CellShape": list(cell)}},
               "file": {"output": {"path_to_output": directory}}, "eliashberg": {}}
        return sc._load_flex_green(inp, 2, 4, 4, 1, interactions=interactions)

    def test_stamped_green_is_accepted(self):
        self.assertIsNotNone(self._load(self.out, {"CoulombInter": _ASYM}))

    def test_legacy_green_is_refused_for_a_sensitive_declaration(self):
        with self.assertRaisesRegex(ValueError, "#232"):
            self._load(self.legacy, {"CoulombInter": _ASYM})

    def test_legacy_green_is_accepted_for_an_insensitive_declaration(self):
        self.assertIsNotNone(self._load(self.legacy, {"CoulombInter": _SYM}))
        # +x and -x coincide on a 2-site axis
        self.assertIsNotNone(self._load(self.legacy, {"CoulombInter": _ASYM}, cell=(2, 4, 1)))


if __name__ == "__main__":
    unittest.main()
