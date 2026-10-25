"""``[mode.param] orbital_site`` (GitHub issue #234).

An ``R = 0`` two-body row between orbitals on different sites of the cell
is an inter-site interaction. With the key it is treated exactly like an
off-site row in the q-only vertices (RPA ring, FLEX general, non-bond
hwave_sc): its density (Hartree) part enters, its exchange crossing does
not. Without the key every orbital is on one site (unchanged results).

The physical check is that the answer no longer depends on the unit cell:
a chain of identical sites with a nearest-neighbour ``V`` is described once
with a one-site cell (``V`` off-site) and once with a two-site cell (one of
the two bonds of each cell is an ``R = 0`` inter-orbital row). With
``orbital_site = [0, 1]`` the two descriptions agree; with the two orbitals
on one site they do not.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver.sites import (cross_site_rows, is_onsite, parse_orbital_site,
                                warn_if_undeclared)


class TestSiteHelpers(unittest.TestCase):

    def test_parse(self):
        self.assertIsNone(parse_orbital_site(None, 2))
        self.assertEqual(parse_orbital_site([0, 1], 2), (0, 1))
        for bad in ([0], [0, 1, 2], "01", [0, -1], [0, 1.0], [True, 0], 3):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    parse_orbital_site(bad, 2)

    def test_is_onsite(self):
        self.assertTrue(is_onsite((0, 0, 0), (0, 1), None))
        self.assertTrue(is_onsite((0, 0, 0), (0, 1), (0, 0)))
        self.assertFalse(is_onsite((0, 0, 0), (0, 1), (0, 1)))
        self.assertTrue(is_onsite((0, 0, 0), (1, 1), (0, 1)))
        self.assertFalse(is_onsite((1, 0, 0), (0, 0), None))

    def test_cross_site_rows_and_warning(self):
        tables = {"CoulombInter": {((0, 0, 0), (0, 1)): 0.3, ((1, 0, 0), (1, 0)): 0.3},
                  "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0}}
        self.assertEqual(cross_site_rows(tables, None), {})
        self.assertEqual(cross_site_rows(tables, (0, 0)), {})
        self.assertEqual(cross_site_rows(tables, (0, 1)),
                         {"CoulombInter": [((0, 0, 0), (0, 1))]})
        with self.assertLogs("qlms", level="WARNING") as logs:
            warn_if_undeclared(tables, None, "RPA")
        self.assertIn("orbital_site", logs.output[0])
        logger = __import__("logging").getLogger("qlms")
        with self.assertLogs(logger, level="WARNING") as logs:
            logger.warning("sentinel")
            warn_if_undeclared(tables, (0, 1), "RPA")
            warn_if_undeclared({"CoulombIntra": tables["CoulombIntra"]}, None, "RPA")
        self.assertEqual(len(logs.output), 1)


# --------------------------------------------------------------------------
# a chain of identical sites, with a one-site and a two-site cell

T_HOP, U, V, NSITE = -1.0, 2.0, 0.6, 8


def _w90(path, title, rows):
    rvecs = sorted({r for r, _, _, _ in rows})
    lines = [title, str(max(max(a, b) for _, a, b, _ in rows) + 1), str(len(rvecs)),
             " ".join(["1"] * len(rvecs))]
    for r, a, b, v in rows:
        lines.append("{} {} {} {} {} {:.12f} 0.0".format(*r, a + 1, b + 1, v))
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _write_chain(d, cell):
    """One-site (``cell=1``) or two-site (``cell=2``) description of the
    same chain; returns the interaction file map."""
    if cell == 1:
        geom = "1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n"
        transfer = [((1, 0, 0), 0, 0, T_HOP), ((-1, 0, 0), 0, 0, T_HOP)]
        intra = [((0, 0, 0), 0, 0, U)]
        inter = [((1, 0, 0), 0, 0, V), ((-1, 0, 0), 0, 0, V)]
    else:
        geom = "2 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0.5 0 0\n"
        bonds = [((0, 0, 0), 0, 1), ((0, 0, 0), 1, 0),       # in-cell bond
                 ((1, 0, 0), 1, 0), ((-1, 0, 0), 0, 1)]      # cross-cell bond
        transfer = [(r, a, b, T_HOP) for r, a, b in bonds]
        intra = [((0, 0, 0), 0, 0, U), ((0, 0, 0), 1, 1, U)]
        inter = [(r, a, b, V) for r, a, b in bonds]
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write(geom)
    _w90(os.path.join(d, "transfer.dat"), "transfer", transfer)
    _w90(os.path.join(d, "coulombintra.dat"), "coulombintra", intra)
    _w90(os.path.join(d, "coulombinter.dat"), "coulombinter", inter)
    return {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat", "CoulombInter": "coulombinter.dat"}


def _run_rpa(d, cell, site, mode="RPA", T=0.5):
    import copy
    import hwave.qlms as qlms
    inter = _write_chain(d, cell)
    param = {"T": T, "CellShape": [NSITE // cell, 1, 1], "SubShape": [1, 1, 1],
             "Nmat": 64, "filling": 0.45, "matsubara_frequency": "center"}
    if mode == "FLEX":
        param.update({"IterationMax": 200, "Mix": 0.5, "EPS": 10})
        param.pop("matsubara_frequency")
    if site is not None:
        param["orbital_site"] = site
    out = os.path.join(d, "output")
    inp = {"log": {"print_level": 0},
           "mode": {"mode": mode, "calc_scheme": "general", "param": param},
           "file": {"input": {"path_to_input": d, "interaction": inter},
                    "output": {"path_to_output": out, "chiq": "chiq", "chi0q": "chi0q"}}}
    qlms.run(input_dict=copy.deepcopy(inp))
    return out


def _uniform_and_staggered(chiq_path, cell):
    """Static charge response of the CHAIN at q = 0 and q = pi, from either
    description (the two-site cell sees both at its Q = 0)."""
    with np.load(chiq_path) as f:
        chiq, fi, wi = f["chiq"], f["freq_index"], f["wavevector_index"]
    l = int(np.where(fi == chiq.shape[0] // 2)[0][0]) if chiq.shape[0] > 1 else 0
    c = chiq[l]
    norb = c.shape[1] // 2
    X = np.einsum("qaabb->qab", c).reshape(-1, 2, norb, 2, norb).sum(axis=(1, 3))
    iq = {tuple(int(x) % (NSITE // cell) for x in w): i for i, w in enumerate(wi)}
    if cell == 1:
        return (X[iq[(0, 0, 0)]][0, 0].real / 1.0,
                X[iq[(NSITE // 2, 0, 0)]][0, 0].real)
    x0 = X[iq[(0, 0, 0)]]
    sign = np.array([1.0, -1.0])
    return (x0.sum().real / 2.0, (sign[:, None] * sign[None, :] * x0).sum().real / 2.0)


class TestCellIndependence(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _dir(self, name):
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        return d

    def test_rpa_ring(self):
        one = _uniform_and_staggered(
            os.path.join(_run_rpa(self._dir("one"), 1, None), "chiq.npz"), 1)
        two_sites = _uniform_and_staggered(
            os.path.join(_run_rpa(self._dir("two"), 2, [0, 1]), "chiq.npz"), 2)
        one_site = _uniform_and_staggered(
            os.path.join(_run_rpa(self._dir("same"), 2, [0, 0]), "chiq.npz"), 2)
        np.testing.assert_allclose(two_sites, one, rtol=1e-10, atol=0)
        # with both orbitals on one site the in-cell bond gets its exchange
        # crossing and the two descriptions disagree (measured 1.1e-4 on the
        # uniform response, against agreement to 1e-15 above)
        gap = max(abs(x - y) for x, y in zip(one_site, one))
        self.assertGreater(gap, 1e-5 * max(abs(x) for x in one))

    def test_ring_plus_ladder_refuses_rows_between_sites(self):
        import hwave.qlms as qlms
        d = self._dir("ladder")
        inter = _write_chain(d, 2)
        inp = {"log": {"print_level": 0},
               "mode": {"mode": "RPA", "calc_scheme": "general", "calc_type": "ring+ladder",
                        "param": {"T": 0.5, "CellShape": [NSITE // 2, 1, 1],
                                  "SubShape": [1, 1, 1], "Nmat": 64, "filling": 0.45,
                                  "orbital_site": [0, 1]}},
               "file": {"input": {"path_to_input": d, "interaction": inter},
                        "output": {"path_to_output": os.path.join(d, "output")}}}
        with self.assertRaisesRegex(ValueError, "orbital_site"):
            qlms.run(input_dict=inp)

    def test_flex_general(self):
        def responses(cell, site):
            out = _run_rpa(self._dir("flex_{}_{}".format(cell, site)), cell, site, mode="FLEX")
            return [_static_density_response(os.path.join(out, name + ".npz"), cell, name)
                    for name in ("chiq_c", "chiq_s")]
        one = responses(1, None)
        two_sites = responses(2, [0, 1])
        one_site = responses(2, [0, 0])
        np.testing.assert_allclose(np.ravel(two_sites), np.ravel(one), rtol=1e-10, atol=0)
        # measured 1.3e-3 on the staggered charge response
        gap = np.max(np.abs(np.ravel(one_site) - np.ravel(one)))
        self.assertGreater(gap, 1e-4 * np.max(np.abs(np.ravel(one))))


def _static_density_response(path, cell, key):
    """Uniform and staggered static density response of the chain from a
    FLEX ``chiq_c`` / ``chiq_s`` file (orbital-pair layout)."""
    with np.load(path) as f:
        c, wi = f[key], f["wavevector_index"]
    X = np.einsum("qaabb->qab", c[c.shape[0] // 2])
    iq = {tuple(int(x) % (NSITE // cell) for x in w): i for i, w in enumerate(wi)}
    if cell == 1:
        return (X[iq[(0, 0, 0)]][0, 0].real, X[iq[(NSITE // 2, 0, 0)]][0, 0].real)
    x0 = X[iq[(0, 0, 0)]]
    sign = np.array([1.0, -1.0])
    return (x0.sum().real / 2.0, (sign[:, None] * sign[None, :] * x0).sum().real / 2.0)


class TestHwaveScVertex(unittest.TestCase):
    """The non-bond hwave_sc vertex: with the two orbitals on different
    sites the in-cell row loses exactly its exchange crossing -- a
    q-independent change confined to the pair-coherence (non-density)
    slots -- and every other element is unchanged."""

    def _sc(self, site):
        from hwave import sc
        d = tempfile.mkdtemp()
        try:
            inter = _write_chain(d, 2)
            interactions = {}
            from hwave.qlmsio import wan90
            for key in ("CoulombIntra", "CoulombInter"):
                interactions[key] = wan90.read_w90(os.path.join(d, inter[key]))
        finally:
            shutil.rmtree(d, ignore_errors=True)
        n = NSITE // 2
        k = np.linspace(0, 2 * np.pi, n, endpoint=False)
        param = {} if site is None else {"orbital_site": site}
        inter_k = sc._build_interaction_k_for_sites(
            k, np.zeros(1), np.zeros(1), interactions, 2, {"mode": {"param": param}})
        return inter_k, sc._build_sc_matrices_all_q(inter_k, 2, n, 1, 1)

    def test_only_the_in_cell_exchange_is_removed(self):
        k_none, (S0, C0) = self._sc(None)
        k_same, (S1, C1) = self._sc([0, 0])
        k_two, (S2, C2) = self._sc([0, 1])
        from hwave.sc import SiteSplitInteractionK
        # every table is split: the chain's inter-cell rows are off-site
        # (issue #239); orbital_site adds the in-cell row between the sites
        for k in (k_none, k_same, k_two):
            self.assertIsInstance(k, SiteSplitInteractionK)
        np.testing.assert_array_equal(S0, S1)
        np.testing.assert_array_equal(C0, C1)
        density = np.zeros(4, dtype=bool)
        density[[0, 3]] = True                    # (00) and (11)
        for M0, M2 in ((S0, S2), (C0, C2)):
            diff = M0 - M2
            self.assertGreater(np.abs(diff).max(), 0.1 * V)
            # density x density block untouched
            np.testing.assert_allclose(diff[..., density, :][..., density], 0, atol=1e-14)
            # the removed part does not depend on q (an R = 0 row)
            np.testing.assert_allclose(diff, np.broadcast_to(diff[:1], diff.shape), atol=1e-14)

    def _tables(self):
        from hwave.qlmsio import wan90
        d = tempfile.mkdtemp()
        try:
            inter = _write_chain(d, 2)
            return {key: wan90.read_w90(os.path.join(d, inter[key]))
                    for key in ("CoulombIntra", "CoulombInter")}
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_single_q_builder_honours_the_split(self):
        from hwave import sc
        n = NSITE // 2
        k = np.linspace(0, 2 * np.pi, n, endpoint=False)
        inter_k = sc._build_interaction_k_for_sites(
            k, np.zeros(1), np.zeros(1), self._tables(), 2,
            {"mode": {"param": {"orbital_site": [0, 1]}}})
        S_all, C_all = sc._build_sc_matrices_all_q(inter_k, 2, n, 1, 1)
        for iq in range(n):
            S1, C1 = sc._build_sc_matrices(inter_k, 2, iq, 0, 0)
            np.testing.assert_allclose(S1, S_all[iq, 0, 0], atol=1e-14)
            np.testing.assert_allclose(C1, C_all[iq, 0, 0], atol=1e-14)

    def test_reduced_flex_warning_names_only_types_that_keep_off_density_content(self):
        # Hund only on the in-cell bond (between the two sites): with the
        # sites declared it has no exchange crossing left; nor has
        # CoulombInter, whose inter-cell rows keep their Hartree part only
        # since issue #239 -- nothing is left for the warning to name
        from hwave import sc
        tables = self._tables()
        tables["Hund"] = {((0, 0, 0), (0, 1)): 0.2, ((0, 0, 0), (1, 0)): 0.2}
        n = NSITE // 2
        k = np.linspace(0, 2 * np.pi, n, endpoint=False)
        inter_k = sc._build_interaction_k_for_sites(
            k, np.zeros(1), np.zeros(1), tables, 2,
            {"mode": {"param": {"orbital_site": [0, 1]}}})
        with self.assertNoLogs("hwave_sc", level="WARNING"):
            sc._warn_reduced_flex_missing_components(inter_k, 2, n, 1, 1)


class TestKeyHandling(unittest.TestCase):

    def test_direct_construction_reads_the_key_case_insensitively(self):
        from hwave.solver.rpa import Interaction, Lattice
        d = tempfile.mkdtemp()
        try:
            inter = _write_chain(d, 2)
            from hwave.qlmsio.read_input_k import QLMSkInput
            ham = QLMSkInput({"path_to_input": d, "interaction": inter}).get_param("ham")
        finally:
            shutil.rmtree(d, ignore_errors=True)
        param = {"T": 0.5, "CellShape": [NSITE // 2, 1, 1], "SubShape": [1, 1, 1],
                 "Nmat": 8, "Orbital_Site": [0, 1]}
        from requests.structures import CaseInsensitiveDict
        lattice = Lattice(CaseInsensitiveDict(param))
        inter_obj = Interaction(lattice, ham, {"mode": "RPA", "param": param})
        self.assertEqual(inter_obj.orbital_site, (0, 1))

    def _ham(self, so=False):
        d = tempfile.mkdtemp()
        try:
            inter = _write_chain(d, 2)
            if so:
                # spin-orbital transfer: interleaved index 2*orb + spin
                rows = [((0, 0, 0), 0, 1), ((0, 0, 0), 1, 0),
                        ((1, 0, 0), 1, 0), ((-1, 0, 0), 0, 1)]
                _w90(os.path.join(d, "transfer.dat"), "transfer",
                     [(r, 2 * a + s, 2 * b + s, T_HOP) for r, a, b in rows for s in range(2)])
                with open(os.path.join(d, "geom.dat"), "w") as f:
                    f.write("2 0 0\n0 1 0\n0 0 1\n4\n0 0 0\n0 0 0\n0.5 0 0\n0.5 0 0\n")
            from hwave.qlmsio.read_input_k import QLMSkInput
            return QLMSkInput({"path_to_input": d, "interaction": inter}).get_param("ham")
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_spin_orbital_labels_run_over_physical_orbitals(self):
        from hwave.solver.rpa import Interaction, Lattice
        from requests.structures import CaseInsensitiveDict
        param = {"T": 0.5, "CellShape": [NSITE // 2, 1, 1], "SubShape": [1, 1, 1], "Nmat": 8}
        mode = {"mode": "RPA", "enable_spin_orbital": True}
        ok = Interaction(Lattice(CaseInsensitiveDict(param)), self._ham(so=True),
                         dict(mode, param=dict(param, orbital_site=[0, 1])))
        self.assertEqual(ok.orbital_site, (0, 1))
        with self.assertRaisesRegex(ValueError, "orbital_site"):
            Interaction(Lattice(CaseInsensitiveDict(param)), self._ham(so=True),
                        dict(mode, param=dict(param, orbital_site=[0, 1, 2, 3])))

    def test_split_is_judged_on_the_pre_fold_table_under_subshape(self):
        from hwave.solver.offsite import split_locality
        from hwave.solver.rpa import Interaction, Lattice
        from requests.structures import CaseInsensitiveDict
        param = {"T": 0.5, "CellShape": [NSITE // 2, 1, 1], "SubShape": [2, 1, 1], "Nmat": 8}
        lattice = Lattice(CaseInsensitiveDict(param))
        inter = Interaction(lattice, self._ham(), {"mode": "RPA", "param": param})
        no_site = split_locality(inter, lattice)
        with_site = split_locality(inter, lattice, site=(0, 1))
        self.assertNotIn(((0, 0, 0), (0, 1)), no_site.offsite_prefold_tbl.get("CoulombInter", {}))
        # without the key the in-cell row stays on-site; the folded on-site
        # table holds its copies in both sub-cells (orbitals 0-1 and 2-3)
        self.assertEqual(sorted(no_site.onsite_tbl["CoulombInter"]),
                         [((0, 0, 0), (0, 1)), ((0, 0, 0), (1, 0)),
                          ((0, 0, 0), (2, 3)), ((0, 0, 0), (3, 2))])
        self.assertFalse(with_site.onsite_tbl.get("CoulombInter"))
        self.assertIn(((0, 0, 0), (0, 1)), with_site.offsite_prefold_tbl["CoulombInter"])
        self.assertIn(((1, 0, 0), (1, 0)), with_site.offsite_prefold_tbl["CoulombInter"])

    def test_hwave_sc_internal_solver_does_not_warn_twice(self):
        from hwave.solver.rpa import Interaction, Lattice
        d = tempfile.mkdtemp()
        try:
            inter = _write_chain(d, 2)
            from hwave.qlmsio.read_input_k import QLMSkInput
            ham = QLMSkInput({"path_to_input": d, "interaction": inter}).get_param("ham")
        finally:
            shutil.rmtree(d, ignore_errors=True)
        param = {"T": 0.5, "CellShape": [NSITE // 2, 1, 1], "SubShape": [1, 1, 1], "Nmat": 8}
        from requests.structures import CaseInsensitiveDict
        logger = __import__("logging").getLogger("qlms")
        with self.assertLogs(logger, level="WARNING") as logs:
            logger.warning("sentinel")
            Interaction(Lattice(CaseInsensitiveDict(param)), ham,
                        {"mode": "RPA", "param": param, "_orbital_site_reported": True})
        self.assertFalse([m for m in logs.output if "orbital_site" in m], logs.output)


if __name__ == "__main__":
    unittest.main()
