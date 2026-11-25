"""Model facts of ``hwave.initial`` (``hwave.initial.facts``)."""
import dataclasses
import hashlib
import os
import shutil
import tempfile
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures
from tests.initial_reference_capture import _read_toml, dump_toml

EXTERN_WARNING = ("Extern is not used by UHFk; give a Zeeman field as on-site "
                  "Transfer in spin-orbital mode")


def _sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_facts_")

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.root, *parts)

    def edit(self, toml, fn):
        doc = _read_toml(toml)
        fn(doc)
        with open(toml, "w") as f:
            f.write(dump_toml(doc))
        return toml


def _two_orbital_chain(dest, *, U0=4.0, subshape=(2, 1, 1), **kw):
    """Two orbitals, U only on orbital 0, hopping on both."""
    transfer = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
                ((1, 0, 0), 1, 1): -0.5, ((-1, 0, 0), 1, 1): -0.5,
                ((0, 0, 0), 0, 1): -0.2, ((0, 0, 0), 1, 0): -0.2}
    param = {"Ncond": 8, "T": 0.1}
    if subshape is not None:
        param["SubShape"] = list(subshape)
    param.update(kw.pop("param", {}))
    return initial_fixtures.write_input(
        dest, cellshape=(4, 1, 1), norb=2, transfer=transfer,
        centers=kw.pop("centers", [(0, 0, 0), (0.5, 0, 0)]),
        interactions={"CoulombIntra": {((0, 0, 0), 0, 0): U0}}, param=param, **kw)


class TestFacts(_Tmp):
    def test_fields_and_frozen(self):
        from hwave.initial import facts
        toml = _two_orbital_chain(self.path("a"))
        f = facts.model_facts(toml)
        self.assertEqual(f.problem.subshape, (2, 1, 1))
        self.assertEqual((f.norb, f.norb_phys, f.spin_orbital, f.iflag_fock), (2, 2, False, True))
        self.assertEqual(f.cellshape, (4, 1, 1))
        self.assertEqual(f.T, 0.1)
        self.assertEqual(f.ncond, 8)
        self.assertEqual(f.ncond_per_cell, Fraction(2))
        self.assertIsNone(f.two_sz)
        self.assertEqual(f.warnings, ())
        with self.assertRaises(dataclasses.FrozenInstanceError):
            f.T = 1.0
        self.assertIsInstance(f.interacting_orbitals, tuple)
        self.assertIsInstance(f.site_groups.groups, tuple)
        # the subshape override reaches the problem
        self.assertEqual(facts.model_facts(toml, subshape=(1, 1, 1)).problem.subshape,
                         (1, 1, 1))

    def test_interacting_orbitals_prefold(self):
        """SubShape 2x1 folds orbital 0 into 0 and 2; the facts read the
        primitive-cell tables."""
        from hwave.initial import facts
        f = facts.model_facts(_two_orbital_chain(self.path("a")))
        self.assertEqual(f.problem.norb, 4)                  # folded
        self.assertEqual(f.interacting_orbitals, (0,))
        self.assertEqual(f.n_c, 1)
        f0 = facts.model_facts(_two_orbital_chain(self.path("b"), U0=0.0))
        self.assertEqual(f0.interacting_orbitals, ())
        self.assertEqual(f0.n_c, 0)

    def test_interacting_orbitals_spin_orbital(self):
        """Two-body indices are physical in SO mode: U on physical orbital 1
        makes spin-orbitals 2 and 3 interacting."""
        from hwave.initial import facts
        transfer = {((1, 0, 0), a, a): -1.0 for a in range(4)}
        transfer.update({((-1, 0, 0), a, a): -1.0 for a in range(4)})
        toml = initial_fixtures.write_input(
            self.path("so"), cellshape=(4, 1, 1), norb=4, transfer=transfer,
            interactions={"CoulombIntra": {((0, 0, 0), 1, 1): 3.0}},
            mode={"enable_spin_orbital": True}, param={"Ncond": 8, "T": 0.1})
        f = facts.model_facts(toml)
        self.assertTrue(f.spin_orbital)
        self.assertEqual((f.norb, f.norb_phys), (4, 2))
        self.assertEqual(f.interacting_orbitals, (2, 3))
        self.assertEqual(f.n_c, 1)

    def test_filling_resolves_like_solver(self):
        from hwave.initial import facts
        from hwave.solver import initial_facade as fac
        variants = {
            "Filling": ({"Filling": 0.25}, 8),
            "filling": ({"filling": 0.25}, 8),
            "round": ({"Filling": 0.3, "Ncond_round_mode": "round"}, 10),
            "round-down": ({"FILLING": 0.3, "ncond_round_mode": "round-down"}, 9),
        }
        for name, (extra, want) in variants.items():
            with self.subTest(variant=name):
                toml = initial_fixtures.square(self.path(name), 4, U=2.0, T=0.1)

                def edit(doc, extra=extra):
                    del doc["mode"]["param"]["Ncond"]
                    doc["mode"]["param"].update(extra)

                self.edit(toml, edit)
                f = facts.model_facts(toml)
                solver, _, _ = fac._build_solver(toml)
                self.assertEqual(f.ncond, solver.param_mod["Ncond"])
                self.assertEqual(f.ncond, want)
                self.assertEqual(f.ncond_per_cell, Fraction(want, 16))
        # 2Sz: (Ncond + 2Sz) / 2 up
        f = facts.model_facts(initial_fixtures.square(self.path("sz"), 4, U=2.0, ncond=15,
                                                      two_sz=1))
        self.assertEqual((f.ncond, f.two_sz), (15, 1))
        self.assertEqual(f.problem.nconds, (8, 7))
        # strict rounding of a non-integer filling: refused like the solver
        toml = initial_fixtures.square(self.path("strict"), 4)
        self.edit(toml, lambda d: (d["mode"]["param"].pop("Ncond"),
                                   d["mode"]["param"].update({"Filling": 0.3})))
        with self.assertRaises(ValueError):
            facts.model_facts(toml)

    def test_extern_warning_text(self):
        from hwave.initial import facts
        toml = initial_fixtures.square(self.path("ext"), 4, U=2.0, subshape=(2, 2, 1))
        with open(self.path("ext", "extern.dat"), "w") as fh:
            fh.write("Extern\n1\n1\n 1\n0 0 0 1 1 0.1 0.0\n")
        self.edit(toml, lambda d: d["file"]["input"]["interaction"].update(
            {"Extern": "extern.dat"}))
        f = facts.model_facts(toml)
        self.assertIn(EXTERN_WARNING, f.warnings)
        self.assertEqual(len(f.warnings), 1)

    def test_no_subshape_warning(self):
        from hwave.initial import facts
        f = facts.model_facts(initial_fixtures.square(self.path("nosub"), 4, U=2.0))
        self.assertEqual(len(f.warnings), 1)
        self.assertIn("SubShape", f.warnings[0])
        self.assertEqual(f.problem.subshape, (4, 4, 1))
        g = facts.model_facts(initial_fixtures.square(self.path("sub"), 4, U=2.0,
                                                      subshape=(2, 2, 1)))
        self.assertEqual(g.warnings, ())
        # a lower-case key is read by the solver (ParamDict): no warning
        toml = initial_fixtures.square(self.path("lower"), 4, U=2.0, subshape=(2, 2, 1))
        self.edit(toml, lambda d: d["mode"]["param"].update(
            {"subshape": d["mode"]["param"].pop("SubShape")}))
        self.assertEqual(facts.model_facts(toml).warnings, ())

    def test_input_sha256(self):
        from hwave.initial import facts
        toml = initial_fixtures.square(self.path("sha"), 4, U=2.0, subshape=(1, 1, 1))
        d = self.path("sha")
        os.makedirs(os.path.join(d, "extra"))
        with open(os.path.join(d, "extra", "geometry_uhf.dat"), "w") as fh:
            fh.write("1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
        with open(os.path.join(d, "extra", "greenone.dat"), "w") as fh:
            fh.write("0 0 0 0 0 0 0 0\n")
        np.savez(os.path.join(d, "extra", "initial.npz"), green=np.zeros(1))
        self.edit(toml, lambda doc: doc["file"]["input"].update(
            {"path_to_input": "extra", "geometry_uhf": "geometry_uhf.dat",
             "onebodyg_uhf": "greenone.dat", "initial": "initial.npz"}))
        f = facts.model_facts(toml)
        want = {"Geometry": _sha(os.path.join(d, "geom.dat")),
                "Transfer": _sha(os.path.join(d, "transfer.dat")),
                "CoulombIntra": _sha(os.path.join(d, "coulombintra.dat")),
                "geometry_uhf": _sha(os.path.join(d, "extra", "geometry_uhf.dat")),
                "onebodyg_uhf": _sha(os.path.join(d, "extra", "greenone.dat"))}
        self.assertEqual(dict(f.input_sha256), want)
        with self.assertRaises(TypeError):
            f.input_sha256["x"] = "y"
        # a missing referenced file is named
        os.remove(os.path.join(d, "extra", "greenone.dat"))
        with self.assertRaisesRegex(ValueError, "onebodyg_uhf"):
            facts.model_facts(toml)


class TestInputPath(_Tmp):
    def test_directory_with_one_uhfk_toml(self):
        from hwave.initial import facts
        toml = initial_fixtures.square(self.path("one"), 4, U=1.0, subshape=(1, 1, 1))
        f = facts.model_facts(self.path("one"))
        self.assertEqual(f.problem.input_files[0].parent.resolve(),
                         facts.model_facts(toml).problem.input_files[0].parent.resolve())

    def test_directory_with_two_tomls_refused(self):
        from hwave.initial import facts
        initial_fixtures.square(self.path("two"), 4, U=1.0, subshape=(1, 1, 1))
        shutil.copy(self.path("two", "input.toml"), self.path("two", "other.toml"))
        with self.assertRaises(ValueError) as cm:
            facts.model_facts(self.path("two"))
        msg = str(cm.exception)
        self.assertIn("input.toml", msg)
        self.assertIn("other.toml", msg)
        # whatever the modes
        self.edit(self.path("two", "other.toml"),
                  lambda d: d["mode"].update({"mode": "UHFr"}))
        with self.assertRaisesRegex(ValueError, "other.toml"):
            facts.model_facts(self.path("two"))

    def test_directory_refusals(self):
        from hwave.initial import facts
        os.makedirs(self.path("empty"))
        with self.assertRaisesRegex(ValueError, "no .*toml"):
            facts.model_facts(self.path("empty"))
        initial_fixtures.square(self.path("uhfr"), 4, U=1.0)
        self.edit(self.path("uhfr", "input.toml"),
                  lambda d: d["mode"].update({"mode": "UHFr"}))
        with self.assertRaisesRegex(ValueError, "UHFr"):
            facts.model_facts(self.path("uhfr"))
        with self.assertRaises(ValueError):
            facts.model_facts(self.path("missing"))


class TestSiteGroups(_Tmp):
    RVEC = np.eye(3)

    def test_periodic_join(self):
        from hwave.initial.facts import cluster_sites
        centers = np.array([[0.0005, 0, 0], [0.9996, 0, 0], [0.5, 0, 0]])
        groups, ambiguous = cluster_sites(centers, self.RVEC, 1e-3)
        self.assertEqual(groups, ((0, 1), (2,)))
        self.assertEqual(ambiguous, ())

    def test_wide_chain_ambiguous(self):
        """0 - 1 - 2 with steps of 0.0008: transitive, 0.0016 wide."""
        from hwave.initial.facts import cluster_sites
        centers = np.array([[0.1, 0, 0], [0.1008, 0, 0], [0.1016, 0, 0], [0.6, 0, 0]])
        groups, ambiguous = cluster_sites(centers, self.RVEC, 1e-3)
        self.assertEqual(groups, ((0, 1, 2), (3,)))
        self.assertEqual(len(ambiguous), 1)
        self.assertIn("wider", ambiguous[0])
        self.assertIn("0, 1, 2", ambiguous[0])

    def test_two_images_of_one_orbital_ambiguous(self):
        """A lattice vector shorter than site_tol joins an orbital with its
        own image; a chain closing around the cell joins two images too."""
        from hwave.initial.facts import cluster_sites
        rvec = np.array([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 5e-4]])
        groups, ambiguous = cluster_sites(np.array([[0.0, 0, 0], [0.5, 0, 0]]), rvec, 1e-3)
        self.assertEqual(groups, ((0,), (1,)))
        self.assertEqual(len(ambiguous), 2)
        self.assertTrue(all("images" in a for a in ambiguous), ambiguous)
        # a ring around a short axis: 0 - 1 - 2 - 0' (period 0.0024 along x)
        rvec = np.array([[0.0024, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
        centers = np.array([[0.0, 0, 0], [1 / 3, 0, 0], [2 / 3, 0, 0]])
        groups, ambiguous = cluster_sites(centers, rvec, 1e-3)
        self.assertEqual(groups, ((0, 1, 2),))
        self.assertTrue(any("images" in a for a in ambiguous), ambiguous)

    def test_model_facts_sources(self):
        from hwave.initial import facts
        toml = _two_orbital_chain(self.path("a"), centers=[(0.0005, 0, 0), (0.9996, 0, 0)])
        f = facts.model_facts(toml)
        self.assertEqual(f.site_groups, facts.SiteGroups("clustered", ((0, 1),), ()))
        # [mode.param] orbital_site wins over the clustering
        self.edit(toml, lambda d: d["mode"]["param"].update({"orbital_site": [0, 1]}))
        f = facts.model_facts(toml)
        self.assertEqual(f.site_groups, facts.SiteGroups("orbital_site", ((0,), (1,)), ()))
        # a user list overrides both
        f = facts.model_facts(toml, orbital_site=[3, 3])
        self.assertEqual(f.site_groups, facts.SiteGroups("user", ((0, 1),), ()))
        with self.assertRaises(ValueError):
            facts.model_facts(toml, orbital_site=[0])


if __name__ == "__main__":
    unittest.main()
