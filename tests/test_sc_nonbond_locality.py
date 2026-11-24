"""The non-bond hwave_sc paths keep only the Hartree part of an off-site row
(issue #239, step 3).

The exchange (Fock) crossing of a row with ``R != 0`` pairs a particle-hole
pair on the bond and depends on two momenta; the q-only vertex of the
non-bond paths (``bond_channels = false``) used to collapse it into the
cross slot ``(ab, ab)``, an approximation. It is now left out, as in the RPA
ring and the FLEX general path, and:

* the non-bond spin/charge matrices are exactly the FLEX general path's
  locality-split ones, and exactly the channel-0 block of the bond-resolved
  vertex, for every density type, on- and off-site, on an inter-orbital
  two-orbital chain;
* an off-site Exchange or PairHop, which has no density part and would
  contribute nothing, is refused on both non-bond entries (static and
  dynamic), pointing to ``bond_channels = true``;
* an off-site inter-orbital CoulombInter, Hund or Ising is named in ONE
  warning (its Fock crossing is not included); same-orbital or on-site rows
  are silent;
* with ``orbital_site``, an ``R = 0`` row between orbitals of different
  sites follows the same rules.
"""
import logging
import os
import shutil
import tempfile
import unittest
from types import SimpleNamespace

import numpy as np

import hwave.sc as sc
from hwave.solver import bond_channels as bc
from hwave.solver.offsite import sc_matrices_from_split, split_locality

NORB, NX = 2, 6


def _tables():
    """An inter-orbital two-orbital chain: on-site U and U', and off-site
    CoulombInter, Hund and Ising rows (each with its partner) that are not
    mirror symmetric."""
    t = {"CoulombIntra": {((0, 0, 0), (a, a)): 2.0 for a in range(NORB)},
         "CoulombInter": {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7}}
    for typ, v in (("CoulombInter", 0.4), ("Hund", 0.15), ("Ising", 0.05)):
        rows = t.setdefault(typ, {})
        rows[((1, 0, 0), (0, 1))] = v
        rows[((-1, 0, 0), (1, 0))] = v
        rows[((1, 0, 0), (0, 0))] = 0.5 * v
        rows[((-1, 0, 0), (0, 0))] = 0.5 * v
    return t


def _nonbond(tables):
    k = np.linspace(0, 2 * np.pi, NX, endpoint=False)
    inter_k = sc._build_interaction_k_for_sites(k, np.zeros(1), np.zeros(1), tables, NORB,
                                                {"mode": {"param": {}}})
    return sc._build_sc_matrices_all_q(inter_k, NORB, NX, 1, 1)


class TestNonbondMatrices(unittest.TestCase):

    def test_equal_the_flex_general_and_the_bond_channel0_matrices(self):
        tables = _tables()
        S, C = _nonbond(tables)
        split = split_locality(SimpleNamespace(param_ham=tables),
                               SimpleNamespace(subshape=(1, 1, 1)))
        S_f, C_f = sc_matrices_from_split(split, bc._LONGITUDINAL_ACTIVE_TYPES, NORB, NX, 1, 1)
        topo, view, bsplit = sc._resolve_bond_topology(tables, {"rvec": np.eye(3)}, NORB, None)
        S_b, C_b = sc._build_bond_vertices(topo, bsplit, NORB, (NX, 1, 1))
        nd = NORB * NORB
        S_b = np.asarray(S_b).reshape(NX, 1, 1, *np.shape(S_b)[-2:])[..., :nd, :nd]
        C_b = np.asarray(C_b).reshape(NX, 1, 1, *np.shape(C_b)[-2:])[..., :nd, :nd]
        for name, mine, theirs in (("S flex", S, S_f), ("C flex", C, C_f),
                                   ("S bond", S, S_b), ("C bond", C, C_b)):
            with self.subTest(name):
                # exactly: the same tables, summed in the same order
                np.testing.assert_array_equal(mine, np.asarray(theirs).reshape(mine.shape))

    def test_the_offsite_fock_is_not_collapsed(self):
        """The off-site rows reach the density slots only: the cross slots
        hold the on-site U' alone (q independent)."""
        S, C = _nonbond(_tables())
        for M in (S, C):
            cross = M[..., [1, 2], [1, 2]]
            np.testing.assert_allclose(cross, np.broadcast_to(cross[:1], cross.shape),
                                       atol=1e-14)


# ---------------------------------------------------------------------------
# the entries
# ---------------------------------------------------------------------------

RPA_INPUT = "tests/rpa/input"


def _write_nn(path, value, orbs=(0, 0)):
    rows = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    with open(path, "w") as fw:
        fw.write("generated for test\n1\n{}\n{}\n".format(len(rows), " ".join(["1"] * len(rows))))
        for r in rows:
            fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {:.12f} 0.0\n".format(
                r[0], r[1], r[2], orbs[0] + 1, orbs[1] + 1, value))


class TestEntries(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="hwave_sc_nonbond_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _input(self, extra, frequency="static"):
        inter = {"path_to_input": RPA_INPUT, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
        for t, v in extra.items():
            p = os.path.join(self.tmp, t + ".dat")
            _write_nn(p, v)
            inter[t] = p
        eli = {"chi0q_mode": "calc", "solver_mode": "iteration", "max_iter": 20,
               "convergence_tol": 1e-6}
        if frequency == "dynamic":
            eli = {"chi0q_mode": "flex", "frequency": "dynamic", "solver_mode": "eigenvalue"}
        return {"mode": {"param": {"T": 1.0, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
                                   "Nmat": 32, "filling": 0.5}},
                "file": {"input": {"interaction": inter, "path_to_flex_output": self.tmp},
                         "output": {"path_to_output": self.tmp}},
                "eliashberg": eli}

    def test_offsite_exchange_and_pairhop_are_refused(self):
        for frequency in ("static", "dynamic"):
            for t in ("Exchange", "PairHop"):
                with self.subTest(frequency=frequency, t=t):
                    with self.assertRaisesRegex(ValueError,
                                                r"off-site {}.*bond_channels = true".format(t)):
                        sc.calc_eliashberg(self._input({"CoulombInter": 0.2, t: 0.1}, frequency))

    def test_offsite_interorbital_density_is_warned_once(self):
        tables = _tables()
        with self.assertLogs("hwave_sc", level="WARNING") as logs:
            sc._check_nonbond_offsite(tables, "where")
        hit = [m for m in logs.output if "Fock" in m]
        self.assertEqual(len(hit), 1, logs.output)
        for t in ("CoulombInter", "Hund", "Ising"):
            self.assertIn(t, hit[0])
        self.assertIn("bond_channels = true", hit[0])
        # same-orbital off-site and on-site inter-orbital rows: silent
        quiet = {"CoulombIntra": tables["CoulombIntra"],
                 "CoulombInter": {key: v for key, v in tables["CoulombInter"].items()
                                  if key[1][0] == key[1][1] or key[0] == (0, 0, 0)}}
        with self.assertNoLogs("hwave_sc", level="WARNING"):
            sc._check_nonbond_offsite(quiet, "where")

    def test_cross_site_rows_follow_the_offsite_rules(self):
        """With orbital_site, an R = 0 row between orbitals of different sites
        is off-site: a cross-site Exchange or PairHop is refused, a cross-site
        CoulombInter is named in the Fock warning; the same rows within one
        site are on-site and silent."""
        two_sites, one_site = (0, 1), (0, 0)
        for t in ("Exchange", "PairHop"):
            rows = {t: {((0, 0, 0), (0, 1)): 0.1, ((0, 0, 0), (1, 0)): 0.1}}
            with self.subTest(t=t):
                with self.assertRaisesRegex(ValueError, "off-site {}".format(t)):
                    sc._check_nonbond_offsite(rows, "where", two_sites)
                sc._check_nonbond_offsite(rows, "where", one_site)
        rows = {"CoulombInter": {((0, 0, 0), (0, 1)): 0.5, ((0, 0, 0), (1, 0)): 0.5}}
        with self.assertLogs("hwave_sc", level="WARNING") as logs:
            sc._check_nonbond_offsite(rows, "where", two_sites)
        self.assertEqual(len([m for m in logs.output if "Fock" in m]), 1, logs.output)
        with self.assertNoLogs("hwave_sc", level="WARNING"):
            sc._check_nonbond_offsite(rows, "where", one_site)

    def test_cross_site_exchange_is_refused_by_both_entries(self):
        """End to end: an R = 0 Exchange between the two sites of a dimer
        (orbital_site = [0, 1]) is refused by the static and the dynamic
        non-bond entries."""
        d = os.path.join(self.tmp, "dimer")
        os.makedirs(d)
        shutil.copy(os.path.join("tests/rpa/input_2orb", "geom.dat"), d)
        shutil.copy(os.path.join("tests/rpa/input_2orb", "transfer.dat"), d)
        with open(os.path.join(d, "exchange.dat"), "w") as fw:
            fw.write("in-cell Exchange between the sites\n2\n2\n 1 1\n")
            fw.write("   0    0    0    1    2   0.100000000000   0.000000000000\n")
            fw.write("   0    0    0    2    1   0.100000000000   0.000000000000\n")
        for frequency in ("static", "dynamic"):
            with self.subTest(frequency=frequency):
                inp = self._input({}, frequency)
                inp["file"]["input"]["interaction"] = {
                    "path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "Exchange": "exchange.dat"}
                inp["mode"]["param"]["orbital_site"] = [0, 1]
                with self.assertRaisesRegex(ValueError, "off-site Exchange"):
                    sc.calc_eliashberg(inp)

    def test_bond_path_builds_no_nonbond_table(self):
        """bond_channels = true neither builds the non-bond k-space table nor
        runs its off-site check, on the static and the dynamic route (the
        dynamic one stops at the absent archive, past the point where the
        non-bond entry would have done both)."""
        from unittest import mock
        static = self._input({"CoulombInter": 0.2})
        static["eliashberg"].update({"bond_channels": True, "solver_mode": "eigenvalue",
                                     "num_eigenvalues": 2})
        dynamic = self._input({"CoulombInter": 0.2}, "dynamic")
        dynamic["eliashberg"]["bond_channels"] = True
        with mock.patch.object(sc, "_build_interaction_k_for_sites",
                               side_effect=AssertionError("built on the bond path")), \
                mock.patch.object(sc, "_check_nonbond_offsite",
                                  side_effect=AssertionError("checked on the bond path")):
            sc.calc_eliashberg(static)
            with self.assertRaisesRegex(ValueError, "longitudinal_bond_output_full"):
                sc.calc_eliashberg(dynamic)

    def test_same_orbital_offsite_runs_without_the_warning(self):
        """The single-band fixture's off-site V is same-orbital: no Fock
        warning, and the run completes."""
        with _NoFockWarning(self):
            sc.calc_eliashberg(self._input({"CoulombInter": 0.2}))


class _NoFockWarning(object):
    """Context: no 'Fock' warning on the hwave_sc logger."""

    def __init__(self, case):
        self.case = case

    def __enter__(self):
        self.records = []
        self.handler = logging.Handler()
        self.handler.emit = self.records.append
        logging.getLogger("hwave_sc").addHandler(self.handler)
        return self

    def __exit__(self, *exc):
        logging.getLogger("hwave_sc").removeHandler(self.handler)
        if exc[0] is None:
            msgs = [r.getMessage() for r in self.records if r.levelno >= logging.WARNING]
            self.case.assertFalse([m for m in msgs if "Fock" in m], msgs)
        return False


if __name__ == "__main__":
    unittest.main()
