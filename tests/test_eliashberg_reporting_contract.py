"""Static-solver reporting contract: the requested pairing channel, the
"both" solver mode, the default solver and the temperature sweep.

* the static non-spin-orbital eigenvalue solver reports the leader of the
  REQUESTED parity sector even when more opposite-parity eigenvalues than
  num_eigenvalues lie above it;
* solver_mode = "both": eigen_converged, the returned value and the exit
  status follow the eigenvalue analysis; the iteration's convergence is
  recorded on its own line;
* an [eliashberg] table without solver_mode runs the eigenvalue solver
  (static and spin-orbital static public entries);
* hwave_tsweep does not accept a NaN leading eigenpair as a result nor chain
  its gap as a seed.
"""
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

import hwave.sc as sc


class TestStaticChannelLeader(unittest.TestCase):
    SHAPE = (1, 1, 8, 1, 1)

    def _sector_basis(self, parity):
        n = int(np.prod(self.SHAPE))
        cols = []
        for j in range(n):
            e = np.zeros(n, complex)
            e[j] = 1.0
            p = sc._project_gap_parity(e.reshape(self.SHAPE),
                                       "singlet" if parity == "even"
                                       else "triplet").ravel()
            for c in cols:
                p = p - np.vdot(c, p) * c
            if np.linalg.norm(p) > 1e-10:
                cols.append(p / np.linalg.norm(p))
        return cols

    def test_requested_sector_leader_below_many_opposite_modes(self):
        """Odd (triplet) eigenvalues 0.9, 0.8, 0.7 lie above the even
        (singlet) leader 0.5; num_eigenvalues = 2 holds odd modes only. The
        singlet request must report 0.5 (dense in-sector oracle)."""
        even, odd = self._sector_basis("even"), self._sector_basis("odd")
        lam_even = [0.5, 0.1, -0.2, 0.05, -0.1][:len(even)]
        lam_odd = [0.9, 0.8, 0.7][:len(odd)]
        n = int(np.prod(self.SHAPE))
        M = np.zeros((n, n), complex)
        for lam, v in zip(lam_even + lam_odd, even + odd):
            M += lam * np.outer(v, v.conj())
        from scipy.sparse.linalg import aslinearoperator
        A = aslinearoperator(M)
        vals, gaps = sc._solve_eigenvalue(None, None, 1, 8, 1, 1,
                                          num_eigenvalues=2,
                                          operator=(A, n))
        vals, gaps, lead, lead_gap, found, _note, _sv = sc._channel_eigenpairs(
            vals, gaps, (A, n), self.SHAPE, "singlet", 2)
        self.assertTrue(found)
        self.assertAlmostEqual(lead.real, max(lam_even), delta=1e-10)
        self.assertTrue(sc._is_gap_parity(lead_gap, "singlet"))


def _legacy(tmp, **eli):
    from tests.test_sc_legacy_golden import _run_legacy
    return _run_legacy(tmp, **eli)


class TestBothModeFollowsTheEigenvalueAnalysis(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _run(self, name, **eli):
        from tests.test_sc_legacy_golden import GOLDEN_DIR, _load_generator
        out = os.path.join(self.tmp, name)
        os.makedirs(out)
        shutil.copy(os.path.join(GOLDEN_DIR, "chi0q.npz"),
                    os.path.join(out, "chi0q.npz"))
        inp = _load_generator().build_input(out)
        inp["eliashberg"]["g2_tail"] = False
        inp["eliashberg"].update(eli)
        lam = sc.calc_eliashberg(inp)
        with open(os.path.join(out, "eigenvalue.dat")) as fh:
            return lam, fh.read()

    def test_failed_iteration_good_eigenanalysis(self):
        # the fixture's 40 shifted steps do not converge
        lam, text = self._run("a", solver_mode="both", num_eigenvalues=6)
        self.assertAlmostEqual(lam, 0.237711293, delta=1e-8)
        self.assertIn("# iteration_converged = false", text)
        self.assertIn("# eigen_converged = true", text)

    def test_good_iteration_failed_eigenanalysis(self):
        real = sc._channel_eigenpairs

        def failed(*a, **k):
            vals, gaps, _lead, _gap, _found, _note, _sv = real(*a, **k)
            return vals, gaps, complex(np.nan, np.nan), None, False, None, None
        with mock.patch.object(sc, "_channel_eigenpairs", failed):
            lam, text = self._run("b", solver_mode="both", num_eigenvalues=6,
                                  max_iter=3000)
        self.assertTrue(np.isnan(lam))
        self.assertIn("# iteration_converged = true", text)
        self.assertIn("# eigen_converged = false", text)


class TestDefaultSolverThroughThePublicEntry(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def test_static(self):
        from tests.test_sc_legacy_golden import GOLDEN_DIR, _load_generator
        out = os.path.join(self.tmp, "s")
        os.makedirs(out)
        shutil.copy(os.path.join(GOLDEN_DIR, "chi0q.npz"),
                    os.path.join(out, "chi0q.npz"))
        inp = _load_generator().build_input(out)
        inp["eliashberg"]["g2_tail"] = False
        inp["eliashberg"].pop("solver_mode")
        with mock.patch.object(sc, "_solve_iteration",
                               side_effect=AssertionError("iteration ran")):
            lam = sc.calc_eliashberg(inp)
        self.assertAlmostEqual(lam, 0.237711293, delta=1e-8)
        with open(os.path.join(out, "eigenvalue.dat")) as fh:
            self.assertIn("# Eigenvalue analysis", fh.read())

    def test_spin_orbital_static(self):
        from tests import test_sc_spin_orbital as S
        files = S._write_model(self.tmp, so=True)
        inp = S._input(self.tmp, files, True, {"chi0q_mode": "calc",
                                               "num_eigenvalues": 4})
        with S._quiet(), mock.patch.object(
                sc, "_so_solve_iteration",
                side_effect=AssertionError("iteration ran")):
            lam = sc.calc_eliashberg(inp)
        self.assertTrue(np.isfinite(lam))
        out = inp["file"]["output"]["path_to_output"]
        with open(os.path.join(out, "eigenvalue.dat")) as fh:
            self.assertIn("# Eigenvalue analysis", fh.read())


class TestTsweepRejectsNaN(unittest.TestCase):

    def test_parse_leading_eig_rejects_nan(self):
        import hwave.tsweep as ts
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        p = os.path.join(tmp, "eigenvalue.dat")
        with open(p, "w") as fw:
            fw.write("# index Re Im |ev| match\n0 nan nan nan 1\n")
        with self.assertRaisesRegex(ValueError, "not finite"):
            ts.parse_leading_eig(p)
        with open(p, "w") as fw:
            fw.write("# eigen_converged=false\n# index Re Im |ev| match\n"
                     "0 0.5 0.0 0.5 1\n")
        with self.assertRaisesRegex(ValueError, "eigen_converged"):
            ts.parse_leading_eig(p)



class TestChannelLeaderFromTheSectorSolve(unittest.TestCase):
    """For a kernel that commutes with parity, the static
    non-spin-orbital solver obtains the reported leader from a solve on the
    requested parity sector's own coordinates, whatever the global set held.
    Through the public static entry with the kernel replaced."""

    N = (1, 1, 4, 4, 1)

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _basis(self, pairing_type):
        n = int(np.prod(self.N))
        cols = []
        for j in range(n):
            e = np.zeros(n, complex)
            e[j] = 1.0
            p = sc._project_gap_parity(e.reshape(self.N), pairing_type).ravel()
            for c in cols:
                p = p - np.vdot(c, p) * c
            if np.linalg.norm(p) > 1e-10:
                cols.append(p / np.linalg.norm(p))
        return cols

    def _kernel(self, even_vals, odd_vals, fill=-0.01, leak=0.0):
        even, odd = self._basis("singlet"), self._basis("triplet")
        n = int(np.prod(self.N))
        M = np.zeros((n, n), complex)
        ev = list(even_vals) + [fill] * (len(even) - len(even_vals))
        ov = list(odd_vals) + [fill] * (len(odd) - len(odd_vals))
        for lam, v in zip(ev + ov, even + odd):
            M += lam * np.outer(v, v.conj())
        if leak:
            M += leak * (np.outer(even[0], odd[0].conj())
                         + np.outer(odd[0], even[0].conj()))
        return M

    def _run(self, M, **eli):
        from scipy.sparse.linalg import aslinearoperator
        from tests.test_sc_legacy_golden import GOLDEN_DIR, _load_generator
        out = os.path.join(self.tmp, "o{}".format(len(os.listdir(self.tmp))))
        os.makedirs(out)
        shutil.copy(os.path.join(GOLDEN_DIR, "chi0q.npz"),
                    os.path.join(out, "chi0q.npz"))
        inp = _load_generator().build_input(out)
        inp["eliashberg"]["g2_tail"] = False
        inp["eliashberg"].update(eli)
        n = M.shape[0]
        with mock.patch.object(sc, "_make_kernel_operator",
                               return_value=(aslinearoperator(M), n)):
            lam = sc.calc_eliashberg(inp)
        with open(os.path.join(out, "eigenvalue.dat")) as fh:
            return lam, fh.read()

    def test_omitted_positive_channel_leader(self):
        """Global which='LM' set [+100 (odd), -10 (even)]; the even (singlet)
        leader +5 is outside it. A leader taken from that set would be -10."""
        M = self._kernel([-10.0, 5.0], [100.0])
        lam, text = self._run(M, solver_mode="eigenvalue", num_eigenvalues=2)
        self.assertAlmostEqual(lam, 5.0, delta=1e-9)
        self.assertIn("# eigen_converged = true", text)
        self.assertIn("parity sector", text)

    def test_negative_channel_spectrum(self):
        """Retained global set [-1, -2] (odd); every even eigenvalue is at or
        below -100. A complement-penalty construction scaled by that set returns NaN."""
        M = self._kernel([-100.0, -120.0, -130.0], [-1.0, -2.0], fill=-150.0)
        lam, text = self._run(M, solver_mode="eigenvalue", num_eigenvalues=2)
        self.assertAlmostEqual(lam, -100.0, delta=1e-8)

    def test_non_commuting_kernel_is_not_solved_in_a_sector(self):
        """A kernel that couples the parities does not commute with parity:
        the sector solve is refused by the parity guard; the reported value is
        the retained eigenpair that is dominantly of the requested parity."""
        from tests.approx_util import _caplog
        M = self._kernel([0.9, 0.2], [0.5], leak=0.05)
        with _caplog("WARNING", logger="hwave_sc") as records:
            lam, text = self._run(M, solver_mode="eigenvalue",
                                  num_eigenvalues=4)
        msg = " ".join(r.message for r in records)
        self.assertIn("does not commute with parity", msg)
        self.assertNotIn("parity sector", text)
        ev, vec = np.linalg.eig(M)
        i = int(np.argmax(ev.real))
        self.assertAlmostEqual(lam, ev[i].real, delta=1e-8)


class TestLegacyModelSmallNumEigenvalues(unittest.TestCase):

    def test_singlet_leader_with_two_eigenvalues(self):
        """Legacy static model, num_eigenvalues = 2: the global set holds the
        repulsive singlet mode -2.752 and a triplet 0.2377; reporting by the
        global set gave -2.752 for the singlet. The parity-sector solve
        reports the singlet leader 0.237711293 (the num_eigenvalues = 10
        value)."""
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        got = _legacy(os.path.join(tmp, "r"), solver_mode="eigenvalue",
                      num_eigenvalues=2)
        lines = [l for l in got["eigenvalue.dat"].decode().splitlines()
                 if l.strip() and not l.startswith("#")]
        self.assertAlmostEqual(float(lines[0]), 0.237711293, delta=1e-8)


def _realistic_triplet(out, mode, **eli):
    """Two-orbital fixture with a small CoulombInter, 32x32, T = 1, Nmat 64,
    triplet: the sector leader 0.0030320961 is threefold degenerate across
    parities in the global solve, whose channel row 0 is -0.00339162472."""
    from tests.test_sc_bond import _two_orbital_input, _write_w90
    os.makedirs(out, exist_ok=True)
    ci = os.path.join(out, "ci2.dat")
    _write_w90(ci, 2, [((1, 0, 0), (0, 1), 0.3), ((-1, 0, 0), (1, 0), 0.3),
                       ((0, 1, 0), (0, 0), 0.2), ((0, -1, 0), (0, 0), 0.2)])
    eli.setdefault("num_eigenvalues", 10)
    inp = _two_orbital_input(out, ci, solver_mode=mode, pairing_type="triplet",
                             max_iter=3000, **eli)
    inp["mode"]["param"]["CellShape"] = [32, 32, 1]
    return sc.calc_eliashberg(inp)


class TestListingLeadsWithTheReportedLeader(unittest.TestCase):
    """Every reader of the static eigenvalue.dat gets the reported leader: it
    is row 0 of the table (match 1) and the "# leading_eigenvalue" line."""

    LEAD = 0.0030320961

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.out = os.path.join(cls.tmp, "eig")
        cls.lam = _realistic_triplet(cls.out, "eigenvalue")
        with open(os.path.join(cls.out, "eigenvalue.dat")) as fh:
            cls.text = fh.read()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_parse_leading_eig_returns_the_reported_value(self):
        import hwave.tsweep as ts
        self.assertAlmostEqual(self.lam, self.LEAD, delta=1e-10)
        re, im, match = ts.parse_leading_eig(
            os.path.join(self.out, "eigenvalue.dat"))
        self.assertAlmostEqual(re, self.lam, delta=1e-10)
        self.assertEqual(match, 1)
        self.assertTrue(ts._eig_parseable(self.out, "eigenvalue.dat"))
        self.assertIn("# leading_eigenvalue = ", self.text)

    def test_row_zero_is_the_leader(self):
        rows = [l.split() for l in self.text.splitlines()
                if l.strip() and not l.startswith("#") and len(l.split()) == 5]
        self.assertEqual(int(rows[0][0]), 0)
        self.assertAlmostEqual(float(rows[0][1]), self.lam, delta=1e-10)
        self.assertEqual(int(rows[0][4]), 1)
        self.assertEqual([int(r[0]) for r in rows], list(range(len(rows))))
        # the old channel row 0 is still listed, below
        self.assertTrue(any(abs(float(r[1]) + 0.00339162472) < 1e-10
                            for r in rows[1:]))

    def test_sweep_result_path(self):
        """hwave_tsweep's result path reads the reported leader from a real
        static file (the solver is replaced by a copy of it)."""
        import copy
        import hwave.tsweep as ts
        text = self.text

        def fake_flex(*, input_dict):
            out = input_dict["file"]["output"]["path_to_output"]
            os.makedirs(out, exist_ok=True)
            np.savez(os.path.join(out, "sigma.npz"), sigma=np.ones((1,)),
                     momentum_convention="e_plus_ikR")
            return {"scf_converged": True, "scf_iterations": 3}

        def fake_eli(input_dict):
            out = input_dict["file"]["output"]["path_to_output"]
            os.makedirs(out, exist_ok=True)
            np.savez(os.path.join(out, "gap_dynamic.npz"), gap=np.ones((1,)))
            with open(os.path.join(out, "eigenvalue.dat"), "w") as fw:
                fw.write(text)
        base = {"mode": {"param": {"CellShape": [2, 2, 1], "Nmat": 16,
                                   "filling": 0.5}},
                "file": {"input": {}, "output": {}},
                "eliashberg": {"frequency": "dynamic",
                               "solver_mode": "eigenvalue"},
                "continuation": {"temperatures": [0.01], "output_dir": "sweep",
                                 "run_eliashberg": True, "warm_start": True,
                                 "seed_gap": True,
                                 "summary_file": "lambda_vs_T.dat"}}
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        with mock.patch("hwave.qlms.run", fake_flex), \
                mock.patch("hwave.sc.calc_eliashberg", fake_eli):
            rows = ts.run(copy.deepcopy(base), base_dir=tmp)
        self.assertEqual(rows[0]["status"], "ok")
        self.assertAlmostEqual(rows[0]["re"], self.lam, delta=1e-10)
        self.assertEqual(rows[0]["match"], 1)

    def test_old_format_is_still_read_from_row_zero(self):
        import hwave.tsweep as ts
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        p = os.path.join(tmp, "eigenvalue.dat")
        with open(p, "w") as fw:
            fw.write("# Iteration eigenvalue\n0.5\n# Eigenvalue analysis\n"
                     "# index  Re  Im  |ev|  match\n   0  0.4 0.0 0.4 1\n"
                     "   1  0.3 0.0 0.3 0\n")
        self.assertEqual(ts.parse_leading_eig(p), (0.4, 0.0, 1))


class TestBothModeThreeValues(unittest.TestCase):

    def test_iteration_leader_and_old_row_zero_differ(self):
        """"both": the iteration's scalar (0.42, mocked), the sector leader
        (+5) and the old global row 0 (-10) are three different values. The
        table and "# leading_eigenvalue" carry +5; the iteration block keeps
        0.42 with its iteration_converged line."""
        import hwave.tsweep as ts
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([-10.0, 5.0], [100.0])
        real = sc._solve_iteration

        def fake_iteration(*a, **k):
            gap, _ev, _conv, n = real(*a, **dict(k, max_iter=2))
            return gap, 0.42, True, 7
        with mock.patch.object(sc, "_solve_iteration", fake_iteration):
            lam, text = t._run(M, solver_mode="both", num_eigenvalues=2)
        self.assertAlmostEqual(lam, 5.0, delta=1e-9)
        lines = text.splitlines()
        i = lines.index("# Iteration eigenvalue")
        self.assertEqual(lines[i + 1], "# iteration_converged = true")
        values = [l for l in lines if l.strip() and not l.startswith("#")]
        self.assertAlmostEqual(float(values[0]), 0.42, delta=1e-12)
        p = os.path.join(t.tmp, sorted(os.listdir(t.tmp))[-1], "eigenvalue.dat")
        re, _im, match = ts.parse_leading_eig(p)
        self.assertAlmostEqual(re, 5.0, delta=1e-9)
        self.assertEqual(match, 1)
        rows = [l.split() for l in values[1:]]
        self.assertAlmostEqual(float(rows[0][1]), 5.0, delta=1e-9)
        self.assertTrue(any(abs(float(r[1]) + 10.0) < 1e-9 for r in rows[1:]))


class TestSectorSolveDetails(unittest.TestCase):

    def test_null_opposite_sector_is_not_refused(self):
        """A kernel whose odd sector is analytically null still commutes
        with parity; the per-sector ratio used to report leakage 1.0."""
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([-10.0, 5.0], [], fill=0.0)
        even = t._basis("singlet")
        M = sum(lam * np.outer(v, v.conj())
                for lam, v in zip([-10.0, 5.0] + [-0.01] * (len(even) - 2),
                                  even))
        lam, text = t._run(M, solver_mode="eigenvalue", num_eigenvalues=1)
        self.assertAlmostEqual(lam, 5.0, delta=1e-9)
        self.assertIn("parity sector", text)

    def test_user_spectral_shift_is_applied_to_the_sector_solve(self):
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([-10.0, 5.0], [100.0])
        lam, text = t._run(M, solver_mode="eigenvalue", num_eigenvalues=2,
                           spectral_shift=20.0)
        self.assertAlmostEqual(lam, 5.0, delta=1e-9)
        self.assertIn("selection LR", text)

    def test_no_false_opposite_channel_warning_on_the_sector_path(self):
        from tests.approx_util import _caplog
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([0.1, 0.05], [100.0, 90.0])
        with _caplog("WARNING", logger="hwave_sc") as records:
            lam, _text = t._run(M, solver_mode="eigenvalue", num_eigenvalues=2)
        self.assertAlmostEqual(lam, 0.1, delta=1e-9)
        self.assertFalse([r for r in records
                          if "belongs to the other channel" in r.message])


class TestBondDiagnosticsLeadingCluster(unittest.TestCase):

    def test_cluster_from_the_sector_spectrum(self):
        """The leading cluster is that of the reported state: with a sector
        spectrum given it is clustered there (doublet 5, 5+1e-9), not in the
        global listing (whose row 0 would describe a different state)."""
        prov = {}
        listing = np.array([5.0, 100.0, -10.0])
        kx = np.linspace(0, 2.0 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2.0 * np.pi, 1, endpoint=False)
        sc._bond_diagnostics_record(prov, None, listing, None, 1, kx, kx, kz,
                                    sector_eigenvalues=np.array(
                                        [5.0, 5.0 + 1e-9, 1.0]))
        self.assertEqual(prov["bond_diagnostics_leading_cluster_size"], 2)
        self.assertEqual(prov["bond_diagnostics_leading_cluster_source"],
                         "channel parity-sector spectrum")

    def test_unavailable_without_a_leader(self):
        prov = {}
        kx = np.linspace(0, 2.0 * np.pi, 4, endpoint=False)
        kz = np.linspace(0, 2.0 * np.pi, 1, endpoint=False)
        sc._bond_diagnostics_record(prov, None, np.array([1.0, 0.5]), None, 1,
                                    kx, kx, kz, leader_available=False)
        self.assertEqual(prov["bond_diagnostics_leading_cluster"], "n/a")


def _assert_rows_are_eigenvectors(tc, M, vals, gaps, rtol=1e-8):
    """Every emitted (lambda, v) pair satisfies A v = lambda v, each judged
    on its own scale: ``||A v - lambda v|| / ||v|| <= rtol |lambda| +
    floor``, with ``floor = 10 n eps ||A||`` the round-off a backward-stable
    eigensolver leaves. A larger unrelated eigenvalue elsewhere in the set
    does not loosen the check of a small one."""
    n = M.shape[0]
    floor = 10.0 * n * np.finfo(float).eps * np.linalg.norm(M, 2)
    for lam, g in zip(vals, gaps):
        v = np.asarray(g).ravel()
        r = np.linalg.norm(M @ v - lam * v) / np.linalg.norm(v)
        tc.assertLessEqual(r, rtol * abs(lam) + floor, (lam, r))


class TestListingWithLeaderKeepsTheGlobalSet(unittest.TestCase):
    """Row 0 of the static table is the reported leader; rows 1.. are the
    global eigen-analysis set exactly as the solver returned it. Nothing is
    merged, so a state may be listed twice (once as row 0) but no vector is
    paired with an eigenvalue it does not have."""

    SHAPE = (1, 1, 4, 1, 1)

    def _check(self, M, vals, gaps, lead, lead_gap):
        out_v, out_g = sc._listing_with_leader(vals, gaps, lead, lead_gap)
        self.assertEqual(len(out_v), len(vals) + 1)
        self.assertEqual(out_v[0], lead)
        np.testing.assert_array_equal(out_g[0],
                                      np.asarray(lead_gap).reshape(self.SHAPE))
        np.testing.assert_array_equal(out_v[1:], vals)
        np.testing.assert_array_equal(out_g[1:], gaps)
        _assert_rows_are_eigenvectors(self, M, out_v, out_g)
        return out_v, out_g

    def test_exactly_degenerate_mixed_basis(self):
        """Eigenvalue 5 shared by one odd and two even states, returned by
        the global solve in a parity-mixed basis: rows 1.. keep all three."""
        rng = np.random.default_rng(3)
        P = TestStaticChannelLeader()
        P.SHAPE = self.SHAPE
        even, odd = P._sector_basis("even"), P._sector_basis("odd")
        vecs = [odd[0], even[0], even[1]]
        others = [(1.0, v) for v in even[2:]] + [(0.5, v) for v in odd[1:]]
        M = 5.0 * sum(np.outer(v, v.conj()) for v in vecs)
        for l, v in others:
            M = M + l * np.outer(v, v.conj())
        mix, _ = np.linalg.qr(rng.standard_normal((3, 3)))
        rows = np.stack(vecs, axis=1) @ mix
        gaps = np.array([rows[:, i].reshape(self.SHAPE) for i in range(3)]
                        + [v.reshape(self.SHAPE) for _l, v in others])
        vals = np.array([5.0, 5.0, 5.0] + [l for l, _v in others], complex)
        out_v, out_g = self._check(M, vals, gaps, 5.0 + 0j,
                                   odd[0].reshape(self.SHAPE))
        V5 = np.stack([g.ravel() for v, g in zip(out_v[1:], out_g[1:])
                       if abs(v - 5.0) < 1e-12], axis=1)
        self.assertEqual(np.linalg.matrix_rank(V5, tol=1e-8), 3)

    def test_non_normal_near_degenerate(self):
        """Eigenvalues 5 and 5 + 1e-9 with eigenvectors 1e-7 apart: the
        distinct state stays listed with its own eigenvalue."""
        n = 4
        S = np.eye(n, dtype=complex)
        S[1, 1] = 1e-7                       # v2 = e0 + 1e-7 e1
        S[0, 1] = 1.0
        D = np.diag([5.0, 5.0 + 1e-9, 0.3, -0.2]).astype(complex)
        M = S @ D @ np.linalg.inv(S)
        v1 = S[:, 0] / np.linalg.norm(S[:, 0])
        v2 = S[:, 1] / np.linalg.norm(S[:, 1])
        gaps = np.array([v2.reshape(self.SHAPE), S[:, 2].reshape(self.SHAPE),
                         S[:, 3].reshape(self.SHAPE)])
        vals = np.array([5.0 + 1e-9, 0.3, -0.2], complex)
        self._check(M, vals, gaps, 5.0 + 0j, v1.reshape(self.SHAPE))

    def test_high_dynamic_range(self):
        """vals = [1, 1e12, 1.5], leader 1: the orthogonal 1.5 state keeps
        eigenvalue 1.5 (an unrelated 1e12 eigenvalue does not widen any
        window), and the leader is listed again in the global rows."""
        Q, _ = np.linalg.qr(np.random.default_rng(5).standard_normal((4, 4)))
        lam = [1.0, 1.0e12, 1.5, 0.25]
        M = (Q * np.array(lam)) @ Q.T
        gaps = np.array([Q[:, i].reshape(self.SHAPE) for i in range(3)],
                        complex)
        vals = np.array(lam[:3], complex)
        out_v, _g = self._check(M, vals, gaps, 1.0 + 0j,
                                Q[:, 0].reshape(self.SHAPE))
        self.assertEqual(int(np.sum(np.abs(out_v - 1.0) < 1e-12)), 2)
        self.assertEqual(int(np.sum(np.abs(out_v - 1.5) < 1e-12)), 1)


class TestStaticTableRowCount(unittest.TestCase):
    """The static table is the reported leader followed by the global
    listing: num_eigenvalues + 1 rows, whether or not the leader's state is
    also in the global set; tsweep still parses the leader."""

    def _rows(self, text):
        return [l.split() for l in text.splitlines()
                if l.strip() and not l.startswith("#") and len(l.split()) == 5]

    def _case(self, num_eigenvalues):
        import hwave.tsweep as ts
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([-10.0, 5.0], [100.0])
        lam, text = t._run(M, solver_mode="eigenvalue",
                           num_eigenvalues=num_eigenvalues)
        self.assertAlmostEqual(lam, 5.0, delta=1e-9)
        rows = self._rows(text)
        self.assertEqual(len(rows), num_eigenvalues + 1)
        self.assertAlmostEqual(float(rows[0][1]), 5.0, delta=1e-9)
        self.assertEqual(int(rows[0][4]), 1)
        p = os.path.join(t.tmp, sorted(os.listdir(t.tmp))[-1],
                         "eigenvalue.dat")
        re, _im, match = ts.parse_leading_eig(p)
        self.assertAlmostEqual(re, 5.0, delta=1e-9)
        self.assertEqual(match, 1)
        return rows

    def test_leader_outside_the_global_set(self):
        rows = self._case(2)                 # global LM set: 100, -10
        self.assertEqual(sum(abs(float(r[1]) - 5.0) < 1e-9 for r in rows), 1)

    def test_leader_also_in_the_global_set(self):
        rows = self._case(3)                 # global LM set: 100, -10, 5
        self.assertEqual(sum(abs(float(r[1]) - 5.0) < 1e-9 for r in rows), 2)


class TestChannelLeaderAlwaysFromTheSectorSolve(unittest.TestCase):
    """Static non-SO, parity-commuting kernel: the reported leader always
    comes from the channel's parity-sector solve, whatever the global
    which='LM' set holds. Counterexamples to taking it from that set
    instead, through the public entry."""

    def _run(self, M, pairing_type, num_eigenvalues):
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        calls = []
        real = sc._parity_sector_basis

        def spy(*a, **k):
            calls.append(a)
            return real(*a, **k)
        with mock.patch.object(sc, "_parity_sector_basis", side_effect=spy):
            lam, text = t._run(M, solver_mode="eigenvalue",
                               num_eigenvalues=num_eigenvalues,
                               pairing_type=pairing_type)
        self.assertGreater(len(calls), 0)
        self.assertIn("own coordinates", text)
        self.assertNotIn("certified", text)
        return lam

    def test_omitted_leader(self):
        """The singlet leader +5 is outside the LM set [+100 (odd), -10]."""
        M = TestChannelLeaderFromTheSectorSolve()._kernel([-10.0, 5.0], [100.0])
        self.assertAlmostEqual(self._run(M, "singlet", 2), 5.0, delta=1e-9)

    def test_negative_listed_channel_pairs_below_the_bound(self):
        """Like the 64x64 triplet: the listed odd pairs are negative (-4.48)
        while the odd sector's largest real part, +4.479, is not listed."""
        M = TestChannelLeaderFromTheSectorSolve()._kernel(
            [5.4, -5.4, 4.5, -4.5], [-4.48, 4.479], fill=-0.01)
        self.assertAlmostEqual(self._run(M, "triplet", 5), 4.479, delta=1e-9)

    def test_large_unrelated_eigenvalue(self):
        """LM set [1e12 (odd), -1 (even)]; the even leader +0.9 is omitted.
        A tolerance scaled by the 1e12 mode would accept -1."""
        M = TestChannelLeaderFromTheSectorSolve()._kernel(
            [-1.0, 0.9], [1.0e12], fill=-0.01)
        # the 1e12 mode leaves ~1e-4 of round-off (1e12 x eps) in the even
        # block, so the value is checked at that level: +0.9, not -1
        self.assertAlmostEqual(self._run(M, "singlet", 2), 0.9, delta=1e-3)

    def test_leader_degenerate_across_parities(self):
        """Even {5, 3}, odd {5}: the global solve may return the degenerate
        5 pair in a parity-mixed basis, which no parity test recognises; the
        singlet leader is 5, not 3."""
        M = TestChannelLeaderFromTheSectorSolve()._kernel(
            [5.0, 3.0], [5.0], fill=-0.01)
        self.assertAlmostEqual(self._run(M, "singlet", 3), 5.0, delta=1e-9)

    def test_realistic_triplet_with_many_eigenvalues(self):
        """Two-orbital 32x32 triplet: its leader 0.0030320961 is threefold
        degenerate across parities; with num_eigenvalues = 32 the listed
        triplet pairs do not contain it."""
        tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, tmp, ignore_errors=True)
        lam = _realistic_triplet(os.path.join(tmp, "o"), "eigenvalue",
                                 num_eigenvalues=32)
        self.assertAlmostEqual(lam, 0.0030320961, delta=1e-10)


class TestStaticTableAtTheDimensionCap(unittest.TestCase):
    """num_eigenvalues >= vec_size - 1: the global solve returns at most
    vec_size - 2 eigenpairs, and the sector solve's leader is prepended to
    however many were returned. The global rows keep their parity order
    (channel parity first) and are not merged."""

    def test_row_count_and_order(self):
        t = TestChannelLeaderFromTheSectorSolve()
        t.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, t.tmp, ignore_errors=True)
        M = t._kernel([-10.0, 5.0], [100.0])
        vec_size = M.shape[0]
        for num in (vec_size - 1, vec_size + 4):
            with self.subTest(num_eigenvalues=num):
                lam, text = t._run(M, solver_mode="eigenvalue",
                                   num_eigenvalues=num)
                self.assertAlmostEqual(lam, 5.0, delta=1e-9)
                rows = [l.split() for l in text.splitlines()
                        if l.strip() and not l.startswith("#")
                        and len(l.split()) == 5]
                self.assertEqual(len(rows), 1 + (vec_size - 2))
                self.assertEqual([int(r[0]) for r in rows],
                                 list(range(len(rows))))
                self.assertAlmostEqual(float(rows[0][1]), 5.0, delta=1e-9)
                self.assertEqual(int(rows[0][4]), 1)
                match = [int(r[4]) for r in rows[1:]]
                # channel-parity (even) rows first, then the opposite parity
                self.assertEqual(match, sorted(match, reverse=True))
                self.assertIn(0, match)
                values = [float(r[1]) for r in rows[1:]]
                for v in (100.0, -10.0, 5.0):
                    self.assertTrue(any(abs(x - v) < 1e-8 for x in values), v)


if __name__ == "__main__":
    unittest.main()
