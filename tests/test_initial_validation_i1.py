"""Validation cases of ``hwave.initial`` for the large tier: V1a, V1b, V2a
and V3 (i), each through ``prepare`` (and, for V1, ``hwave`` and
``collect``).

* **V1** square lattice 16 x 16, hopping 1, U = 4, half filling, SubShape
  2 x 2, Fock on, the default candidate set (every reserved and baseline
  run). Each run is solved in this process from its run directory, as
  ``commands.txt`` would. Frozen references: ``tests/initial_ref/v1a_i1.json``
  and ``tests/initial_ref/v1b.json``, written once by
  ``tests/initial_reference_capture.py --validation`` (model, solver
  settings, H-wave and numpy versions recorded in the files).

  - V1a (T = 0.02): the para row ends para; the Neel row converges, its
    free energy per site below the para and ferro rows by more than 1e-3;
    its moment equals the reference to 1e-6. At T = 1.5 T_N (T_N frozen in
    the reference, from a bisection on max |m_i| > 1e-6, and checked here
    against the independent Stoner root U chi0(Q, T) = 1) every candidate
    converges with max |m_i| < 1e-6.
  - V1b (T = 0, open Fermi shell): the shell |n_x| + |n_y| = 8 holds 30 k
    points (60 states, 30 of them filled; 113 k points below), counted on
    the mesh, and the para seed records it; T_hint = 0.25; the lowest
    converged row's final label is Neel, every row whose final label is
    para or ferro lies above it by more than 1e-3, and statuses, final
    labels and free energies per site of every spin-polarized seed equal
    the reference. A row seeded by a pattern without a spin component (the
    para row and the charge rows) has a deterministic seed: its initial
    state has G_up == G_dn bit for bit and no up-down block, and the two
    C4-partner charge seeds (q = (0, 1/2) and (1/2, 0)), whose target lies
    in the jump of a level crossing at the Fermi level near h = 2^(3/4),
    stop at the same h with the same m to 1e-12. Their runs are not
    compared with the reference: the field cannot break the spin symmetry,
    the Hartree-Fock map preserves the spin-symmetric subspace, and UHFk
    leaves it only through its own round-off (the eigenbasis eigh picks
    inside each spin-degenerate pair and the T = 0 tie order), so whether
    and when such a run leaves it depends on the numpy and BLAS build.
    Measured with the deterministic seeds: on macOS (numpy 2.4 and 1.26)
    the para and the two charge rows reach the Neel state in 51 to 78
    iterations and the q = (1/2, 1/2) charge row does not converge in
    IterationMax; on Linux OpenBLAS (numpy 1.26) all four reach the Neel
    state, after 60 to 267 iterations. Such rows still enter the
    lowest-row assertions when they converge.

  Tolerances: the moment to 1e-6 (the stated V1a tolerance); free energies
  per site to 1e-6, per row against the same row of the reference. EPS =
  10 stops each run at its own distance from the fixed point: rows that
  reach the same state from different seeds differ by up to about 2e-7 per
  site, and one row differs by up to about 3e-8 between numpy 1.26 and 2
  (measured); a different state differs by far more than 1e-3. Green
  functions are never compared: they are not reproducible bit for bit
  across numpy and BLAS builds.

* **V2a** (fast) triangular lattice 12 x 12, half filling, the 120-degree
  literature pattern at K = (2/3, 1/3) on SubShape 2 x 2: listed as not
  representable with the exact reason, the 3 x 3 proposal and the
  rotated-cell note, no instability hints.
* **V3 (i)** a user pattern at q = (7/16, 1/2) on a 64 x 32 lattice
  (fast): on 2 x 2 not representable with the proposal 16 x 2 (which
  divides the lattice); with ``max_denominator = 12`` no cell is proposed.
  Heavy: on ``--subshape 16,2,1`` it is written and checked (the
  enumeration on the 32-site cell takes about 14 s).
"""
import io
import json
import os
import re
import shutil
import tempfile
import unittest

import numpy as np

from tests import initial_fixtures
from tests import initial_reference_capture as capture
from tests.heavy_tests import heavy

#: the free-energy tolerance per site against the frozen references
ENERGY_TOL = 1e-6
#: the moment tolerance of V1a
MOMENT_TOL = 1e-6
#: "ends para" / "disordered": max |m_i| below this
PARA_MOMENT = 1e-6
#: the V1 free-energy margin per site
MARGIN = 1e-3

_Q_IN_LABEL = re.compile(r"q=\(([^)]*)\)")


def _load(path):
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def label_qs(label):
    """The ordering vectors a final label names (``q=(...)`` parts)."""
    return _Q_IN_LABEL.findall(label)


def is_neel(label):
    """A final label of a spin state at q = (1/2, 1/2, 0) only (the
    collinear Neel state in any spin direction)."""
    return label_qs(label) == ["1/2,1/2,0"] and "spin-" in label


def is_ferro(label):
    """A final label of a uniform (q = 0 only) spin state."""
    return label_qs(label) == ["0,0,0"] and "spin-" in label


#: the channels of a pattern without a spin component: their field cannot
#: break the spin symmetry, so a run from such a seed leaves the
#: spin-symmetric subspace only through round-off (see the module docstring)
SPINLESS_CHANNELS = ("para", "charge")


def _tried(cov):
    return [r for r in cov["rows"] if r["tried"]]


def _row(cov, prefix):
    rows = [r for r in _tried(cov) if r["run_id"].startswith(prefix)]
    if len(rows) != 1:
        raise AssertionError("{} tried rows start with {}: {}".format(
            len(rows), prefix, [r["run_id"] for r in rows]))
    return rows[0]


#: the run-id prefixes of the reserved para and ferro rows and of the Neel
#: baseline row (q = (1/2, 1/2), spin z, sign +1) on SubShape 2 x 2
PARA, FERRO, NEEL = "para-", "ferro-", "baseline-q1_2_1_2_0-spin-z"


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_validation_")

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)

    def path(self, *parts):
        return os.path.join(self.root, *parts)


class TestV1(_Tmp):
    """V1a and V1b, prepare -> hwave -> collect on the default large tier."""

    def run_case(self, T, name):
        toml = capture.v1_input(self.path(name + "_src"), T)
        return capture.prepare_run_collect(toml, self.path(name + "_work"))

    def assert_default_set(self, cov):
        families = {r["family"] for r in cov["rows"]}
        self.assertEqual(families, {"para", "ferro", "random", "baseline"})
        self.assertEqual(cov["tier"]["name"], "large")

    @heavy
    def test_v1a_i1(self):
        ref = _load(capture.V1A_REFERENCE)
        self.assertEqual(ref["T"], capture.V1A_T)
        cov = self.run_case(capture.V1A_T, "v1a")
        self.assert_default_set(cov)
        para, ferro, neel = _row(cov, PARA), _row(cov, FERRO), _row(cov, NEEL)
        self.assertEqual(neel["run_id"], ref["neel_run_id"])
        # the para row ends para
        self.assertEqual(para["status"], "converged")
        self.assertLess(capture.max_moment(para), PARA_MOMENT)
        self.assertEqual(para["final_label"], "para")
        # the Neel row converges below para and ferro by > 1e-3 per site
        self.assertEqual(neel["status"], "converged")
        self.assertTrue(is_neel(neel["final_label"]), neel["final_label"])
        self.assertEqual(ferro["status"], "converged")
        for other in (para, ferro):
            self.assertGreater(other["free_energy_per_site"] - neel["free_energy_per_site"],
                               MARGIN, other["run_id"])
        # the moment against the frozen reference
        self.assertAlmostEqual(capture.max_moment(neel), ref["neel_moment"],
                               delta=MOMENT_TOL)

        # T_N: the frozen bisection brackets the independent Stoner root
        self.assertLess(ref["T_N_lower"], ref["T_N"])
        self.assertLessEqual(ref["T_N"] - ref["T_N_lower"], 1e-3)
        t_stoner = capture.stoner_t_n()
        self.assertAlmostEqual(t_stoner, ref["T_N_stoner"], delta=1e-10)
        self.assertLessEqual(ref["T_N_lower"] - 1e-6, t_stoner)
        self.assertLessEqual(t_stoner, ref["T_N"] + 1e-6)
        # at 1.5 T_N every candidate converges to a state without moments
        hot = self.run_case(1.5 * ref["T_N"], "v1a_hot")
        rows = _tried(hot)
        self.assertGreater(len(rows), 3)
        for r in rows:
            self.assertEqual(r["status"], "converged", r["run_id"])
            self.assertLess(capture.max_moment(r), PARA_MOMENT, r["run_id"])

    @heavy
    def test_v1b(self):
        ref = _load(capture.V1B_REFERENCE)
        L = capture.V1_MODEL["L"]
        # the shell, counted on the mesh: e(k) = -2 (cos kx + cos ky), E_F = 0
        n = np.arange(-L // 2 + 1, L // 2 + 1)
        nx, ny = np.meshgrid(n, n, indexing="ij")
        e = -2.0 * (np.cos(2 * np.pi * nx / L) + np.cos(2 * np.pi * ny / L))
        shell = np.abs(nx) + np.abs(ny) == L // 2
        self.assertEqual(int(shell.sum()), 30)
        self.assertTrue(np.all(np.abs(e[shell]) < 1e-12))
        self.assertTrue(np.all(np.abs(e[~shell]) > 0.1))
        below = int((e < -1e-12).sum())
        self.assertEqual(below, 113)
        n_shell_electrons = L * L - 2 * below
        self.assertEqual((2 * int(shell.sum()), n_shell_electrons), (60, 30))

        cov = self.run_case(0.0, "v1b")
        self.assert_default_set(cov)
        self.assertEqual(cov["hint"]["T_hint"], 0.25)
        self.assertEqual(cov["hint"]["t_hint_source"], "rule")
        self.assertEqual(cov["hint"]["W0"], 8.0)
        self.assertEqual(cov["hint"]["L_min"], 16)
        para = _row(cov, PARA)
        self.assertEqual(para["shell"]["h0"], [{"group": "all", "open_shell": True,
                                                "shell_degeneracy": 60,
                                                "shell_fill": 0.5}])
        self.assertEqual(para["h_used"], 0.0)
        self.assertTrue(np.allclose(para["moments_h0"]["spin"], 0.0, atol=1e-12))

        conv = [r for r in _tried(cov) if r["status"] == "converged"]
        lowest = min(conv, key=lambda r: r["free_energy_per_site"])
        self.assertTrue(is_neel(lowest["final_label"]), lowest["final_label"])
        judged = [r for r in conv if r["final_label"] == "para" or is_ferro(r["final_label"])]
        self.assertTrue(judged)       # the ferro row ends ferro
        for r in judged:
            self.assertGreater(r["free_energy_per_site"] - lowest["free_energy_per_site"],
                               MARGIN, r["run_id"])
        # statuses, final labels and energies against the frozen reference
        got = {r["run_id"]: r for r in _tried(cov)}
        self.assertEqual(sorted(got), sorted(ref["rows"]))
        symmetric = {rid for rid, r in got.items() if r["channel"] in SPINLESS_CHANNELS}
        # the para row and the three charge rows
        self.assertEqual(sorted(symmetric), sorted(
            rid for rid in got if rid.startswith(PARA) or "-charge-" in rid))
        self.assertEqual(len(symmetric), 4)
        # their seeds are spin-symmetric bit for bit
        for rid in symmetric:
            with np.load(os.path.join(self.path("v1b_work"), got[rid]["run_dir"],
                                      "initial.npz")) as z:
                for key in ("green", "green_sublattice"):
                    g = z[key]
                    self.assertTrue(np.array_equal(g[:, 0, :, 0, :], g[:, 1, :, 1, :]),
                                    (rid, key))
                    self.assertFalse(np.any(g[:, 0, :, 1, :]) or np.any(g[:, 1, :, 0, :]),
                                     (rid, key))
        # the C4 partners stop at the same crossing h
        partners = [got[rid] for rid in sorted(symmetric)
                    if "-charge-q-0-1-2-0-" in rid or "-charge-q-1-2-0-0-" in rid]
        self.assertEqual(len(partners), 2)
        self.assertEqual(partners[0]["h_used"], partners[1]["h_used"])
        self.assertLess(abs(partners[0]["h_used"] - 2 ** 0.75), 1e-6)
        self.assertAlmostEqual(partners[0]["m_used"], partners[1]["m_used"], delta=1e-12)
        compared = 0
        for rid, want in ref["rows"].items():
            r = got[rid]
            if rid in symmetric:
                # UHFk's round-off decides whether it leaves the
                # spin-symmetric subspace (see the module docstring): not
                # compared
                self.assertIn(r["status"], ("converged", "not_converged"), rid)
                continue
            compared += 1
            self.assertEqual(r["status"], want["status"], rid)
            self.assertEqual(r["final_label"], want["final_label"], rid)
            if want["free_energy_per_site"] is None:
                self.assertIsNone(r["free_energy_per_site"], rid)
            else:
                self.assertAlmostEqual(r["free_energy_per_site"],
                                       want["free_energy_per_site"], delta=ENERGY_TOL,
                                       msg=rid)
        self.assertEqual(compared, 8)


#: the 120-degree spiral at K on its own 3 x 3 cell (a literature row)
TRI_120 = {"patternspec_version": 1, "q": ["2/3", "1/3", "0"], "subshape": [3, 3, 1],
           "phase": 0, "sign": 1, "modulation": "spiral", "axes": ["x", "y"],
           "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "s": [1.0, 0.0]}]}

#: a collinear spin-z pattern at q = (7/16, 1/2) on 16 x 2
V3_PATTERN = {"patternspec_version": 1, "q": ["7/16", "1/2", "0"], "subshape": [16, 2, 1],
              "phase": 0, "sign": 1, "modulation": "cos",
              "entries": [{"orb_a": 0, "orb_b": 0, "r": [0, 0, 0], "cz": [1.0, 0.0]}]}


class TestPrepareValidation(_Tmp):
    """V2a and V3 (i) through ``python -m hwave.initial prepare``."""

    def extra(self, family, pattern, label, citation=None):
        item = {"family": family, "pattern": pattern, "label": label}
        if citation is not None:
            item["citation"] = citation
        path = self.path("extra_{}.json".format(label))
        with open(path, "w") as f:
            json.dump([item], f)
        return path

    def main(self, *argv):
        from hwave.initial import cli
        import contextlib
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = cli.main([str(a) for a in argv])
        return rc, out.getvalue(), err.getvalue()

    def coverage(self, work):
        return _load(os.path.join(work, "coverage.json"))

    def rows_of(self, cov, family):
        return [r for r in cov["rows"] if r["family"] == family]

    def test_v2a_cli(self):
        src = initial_fixtures.triangular(self.path("tri"), 12, U=2.0, T=0.1, ncond=144,
                                          subshape=(2, 2, 1))
        extra = self.extra("literature", TRI_120, "120-degree", "a textbook result")
        work = self.path("WORK")
        rc, _, err = self.main("prepare", os.path.dirname(src), "--out", work,
                               "--extra-candidates", extra, "--max-candidates", "0")
        self.assertEqual(rc, 0, err)
        cov = self.coverage(work)
        rows = self.rows_of(cov, "literature")
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertFalse(row["representable"])
        self.assertFalse(row["tried"])
        self.assertEqual(row["not_tried_reason"], "unrepresentable")
        self.assertEqual(row["reason"], "q_x denominator 3 does not divide B_x = 2")
        self.assertEqual(row["proposed_subshape"], [3, 3, 1])
        self.assertIsNone(row["proposed_cellshape"])
        self.assertIn("√3×√3", row["proposal_note"])
        self.assertEqual(row["citation"], "a textbook result")
        self.assertFalse(os.path.exists(os.path.join(work, row["run_id"])))
        # the large tier: no instability hints
        self.assertEqual(cov["hint"]["status"], "skipped")
        self.assertEqual(cov["hint"]["q_set"], [])
        for r in cov["rows"]:
            self.assertIsNone(r["lambda"])
        # the reserved rows are written and checked
        run_ids = [e["run_id"] for e in _load(os.path.join(work, "manifest.json"))["runs"]]
        self.assertTrue(any(rid.startswith(PARA) for rid in run_ids))
        for rid in run_ids:
            rc, out, err = self.main("check", os.path.join(work, rid, "input.toml"),
                                     "--require-provenance")
            self.assertEqual(rc, 0, err)
            self.assertIn("verified", out)

    def v3_input(self):
        """The doped square lattice of V3: 64 x 32, filling 7/8, U = 6,
        T = 0.02, SubShape 2 x 2; and the extra-candidates file of the
        q = (7/16, 1/2) user pattern."""
        src = initial_fixtures.write_input(
            self.path("sq"), cellshape=(64, 32, 1),
            transfer=initial_fixtures.square_transfer(),
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 6.0}},
            param={"Ncond": 1792, "T": 0.02, "SubShape": [2, 2, 1]})
        return src, self.extra("user", V3_PATTERN, "stripe")

    def test_v3i_cli(self):
        from hwave.initial import cli
        src, extra = self.v3_input()
        # on 2 x 2: listed, not representable, the 16 x 2 proposal
        work = self.path("W22")
        rc, _, err = self.main("prepare", src, "--out", work, "--extra-candidates", extra,
                               "--max-candidates", "0")
        self.assertEqual(rc, 0, err)
        row, = self.rows_of(self.coverage(work), "user")
        self.assertFalse(row["representable"])
        self.assertFalse(row["tried"])
        self.assertEqual(row["not_tried_reason"], "unrepresentable")
        self.assertEqual(row["reason"], "q_x denominator 16 does not divide B_x = 2")
        self.assertEqual(row["proposed_subshape"], [16, 2, 1])
        self.assertIsNone(row["proposed_cellshape"])

        # with max_denominator = 12 (a library keyword of prepare): no cell
        work12 = self.path("W12")
        rc = cli.prepare(src, work12, extra_candidates=extra, max_candidates=0,
                         max_denominator=12, stdout=io.StringIO())
        self.assertEqual(rc, 0)
        row, = self.rows_of(self.coverage(work12), "user")
        self.assertFalse(row["tried"])
        self.assertIsNone(row["proposed_subshape"])
        self.assertIsNone(row["proposed_cellshape"])
        self.assertEqual(row["proposal_note"], "q has a denominator above max_denominator "
                                               "= 12: no cell is proposed")

    @heavy
    def test_v3i_cli_proposed_cell(self):
        # on the proposed --subshape 16,2,1: written and checked
        src, extra = self.v3_input()
        work = self.path("W162")
        rc, _, err = self.main("prepare", src, "--out", work, "--extra-candidates", extra,
                               "--max-candidates", "0", "--subshape", "16,2,1")
        self.assertEqual(rc, 0, err)
        row, = self.rows_of(self.coverage(work), "user")
        self.assertTrue(row["representable"])
        self.assertTrue(row["tried"])
        self.assertEqual(row["subshape"], [16, 2, 1])
        toml = os.path.join(work, row["run_id"], "input.toml")
        rc, out, err = self.main("check", toml, "--require-provenance")
        self.assertEqual(rc, 0, err)
        self.assertIn("verified", out)
        import tomli
        with open(toml, "rb") as f:
            doc = tomli.load(f)
        self.assertEqual(doc["mode"]["param"]["SubShape"], [16, 2, 1])
        self.assertEqual(doc["mode"]["param"]["CellShape"], [64, 32, 1])

if __name__ == "__main__":
    unittest.main()
