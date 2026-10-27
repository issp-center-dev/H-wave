"""Diagnostics of a diverging FLEX self-consistency (GitHub issue #282).

A warm start (``sigma_init``) taken across too large a U step can make the
FLEX iteration diverge. The run then stops in the chemical-potential search,
whose message used to blame the filling and the temperature. These tests pin:

* the failed search after a diverged self-energy says "diverged", names the
  likely causes and does NOT give the filling advice;
* a genuinely unreachable target (N(mu) varies, no root) still gives the
  filling/temperature advice;
* the Stoner-factor log line equals an independent computation of
  ``max_q eig[U chi0(q, i nu = 0)]``;
* the diagnostics are log-only: the solver's outputs are bit-identical with
  and without them.
"""

import logging
import re
import unittest

import numpy as np

from hwave.solver import flex_diagnostics as diag

try:
    import sparse_ir  # noqa: F401
    HAVE_SPARSE_IR = True
except ImportError:
    HAVE_SPARSE_IR = False

_LOGGER = "hwave.solver.flex"


def _make_solver(param_extra=None, U=2.0, L=4, T=0.5, Nmat=64,
                 calc_scheme="reduced"):
    """One-orbital square-lattice Hubbard model at half filling."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    param = {'T': T, 'CellShape': [L, L, 1], 'SubShape': [1, 1, 1],
             'Nmat': Nmat, 'IterationMax': 100, 'Mix': 0.3,
             'Ncond': float(L * L)}
    param.update(param_extra or {})
    param = {k: v for k, v in param.items() if v is not None}
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': calc_scheme}
    info_file_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {'path_to_input': 'tests/rpa/input',
                        'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat'},
    }
    read_io = read_input_k.QLMSkInput(info_file_input)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")
    ham_info["CoulombIntra"] = {k: complex(U) for k in ham_info["CoulombIntra"]}
    return solver_flex.FLEX(ham_info, {}, info_mode), green_info


_STONER_RE = re.compile(
    r"Stoner factor max_q eig\[U chi0\(q, i nu = 0\)\]: spin (\S+) at q = "
    r"\(([^)]*)\); charge (\S+) at q")


class TestDivergedWarmStart(unittest.TestCase):

    @unittest.skipUnless(HAVE_SPARSE_IR, "needs the optional sparse-ir package")
    def test_large_u_step_reports_diverged_self_energy(self):
        """U = 2 -> 200 from the converged U = 2 self-energy, with the IR
        bandwidth pinned (ir_wmax = 10): the first map drives Sigma far out
        of the basis, N(mu) no longer depends on mu, and the next mu search
        fails. The message must say "diverged" and must not blame the
        filling or the temperature."""
        ir = {'matsubara_basis': 'ir', 'ir_wmax': 10.0}
        s1, g1 = _make_solver(ir, U=2.0)
        s1.solve(g1, 'tests/flex/output')
        self.assertTrue(s1.scf_converged)

        s2, g2 = _make_solver(dict(ir, Mix=1.0), U=200.0)
        g2["sigma_init"] = g1["sigma"]
        with self.assertRaisesRegex(RuntimeError, "diverged") as cm:
            s2.solve(g2, 'tests/flex/output')
        msg = str(cm.exception)
        self.assertIn("not bracketed", msg)
        self.assertNotIn("Check the target filling", msg)
        self.assertIn("does not depend on mu", msg)
        self.assertIn("smaller U/T step", msg)
        self.assertIn("smaller Mix", msg)
        self.assertIn("cold start", msg)
        self.assertRegex(msg, r"last completed iteration was \d+ with "
                              r"\|dSigma\|/\|Sigma\| = ")
        self.assertRegex(msg, r"Stoner factor spin [0-9.]+")

    def test_non_finite_sigma_reports_diverged(self):
        s, gi = _make_solver(Nmat=32)
        s._calc_epsilon_k(gi)
        nb, nv, nd = s.H0_eigenvalue.shape
        sigma = np.full((nb, s.nmat, nv, nd, nd), np.nan, dtype=np.complex128)
        with self.assertRaisesRegex(RuntimeError, "diverged") as cm, \
                np.errstate(invalid="ignore"):
            s._find_mu_dressed(sigma, 1.0 / s.T, 8.0)
        self.assertIn("self-energy is not finite", str(cm.exception))
        self.assertNotIn("Check the target filling", str(cm.exception))

    def test_unreachable_target_keeps_filling_advice(self):
        """N(mu) runs from 0 to Nstate: a target above Nstate is not
        reachable although N(mu) varies -- the filling message stays."""
        s, gi = _make_solver(Nmat=32)
        s._calc_epsilon_k(gi)
        nb, nv, nd = s.H0_eigenvalue.shape
        sigma = np.zeros((nb, s.nmat, nv, nd, nd), dtype=np.complex128)
        nstate = int(np.prod(s.H0_eigenvalue.shape))
        with self.assertRaisesRegex(RuntimeError,
                                    "Check the target filling/Ncond and temperature") as cm:
            s._find_mu_dressed(sigma, 1.0 / s.T, 10.0 * nstate)
        self.assertNotIn("diverged", str(cm.exception))
        self.assertIn("varies over the bracket", str(cm.exception))


class TestMuFailureMessage(unittest.TestCase):

    def test_flat_is_diverged(self):
        m = diag.mu_failure_message("X", -1e31, 1e31, -8.0, -8.0, 16, True)
        self.assertIn("diverged", m)
        self.assertNotIn("filling/Ncond", m)

    def test_nan_n_is_diverged(self):
        m = diag.mu_failure_message("X", -1.0, 1.0, float("nan"), 3.0, 16, True)
        self.assertIn("N(mu) is not finite", m)

    def test_varying_is_filling(self):
        m = diag.mu_failure_message("X", -1e31, 1e31, -160.0, -144.0, 16, True)
        self.assertIn("filling/Ncond", m)
        self.assertNotIn("diverged", m)

    def test_context_from_monitor(self):
        mon = diag.ScfMonitor(logging.getLogger("test282"))
        mon.record(7, 0.93, {"spin": (1.25, (0.5, 0.5, 0.0))}, last=False)
        m = diag.mu_failure_message("X", -1e31, 1e31, -8.0, -8.0, 16, True, mon)
        self.assertIn("iteration was 7 with |dSigma|/|Sigma| = 9.300e-01", m)
        self.assertIn("spin 1.250000 at q = (0.5, 0.5, 0)", m)


class TestScfMonitor(unittest.TestCase):

    def _records(self, steps, stoner_values):
        lg = logging.getLogger("test282.monitor")
        mon = diag.ScfMonitor(lg)
        with self.assertLogs(lg, level="DEBUG") as cm:
            lg.debug("start")
            for i, (st, lam) in enumerate(zip(steps, stoner_values)):
                mon.record(i + 1, st, {"spin": (lam, (0.5, 0.5, 0.0))},
                           last=i + 1 == len(steps))
        return cm.records, mon

    def test_info_first_and_last_debug_between(self):
        recs, mon = self._records([1.0, 0.3, 0.1, 1e-3], [0.9, 0.9, 0.9, 0.9])
        info = [r for r in recs if r.levelno == logging.INFO]
        debug = [r for r in recs if r.levelno == logging.DEBUG
                 and "Stoner" in r.getMessage()]
        self.assertEqual(len(info), 2)
        self.assertEqual(len(debug), 2)
        self.assertFalse(mon.warned)

    def test_warns_once_on_unstable_large_steps(self):
        recs, mon = self._records([1.0, 0.8, 0.9, 1.2, 2.0, 3.0],
                                  [1.5, 1.6, 1.7, 2.0, 3.0, 4.0])
        warns = [r for r in recs if r.levelno == logging.WARNING]
        self.assertEqual(len(warns), 1)
        self.assertIn("looks divergent", warns[0].getMessage())
        self.assertIn("iteration 3", warns[0].getMessage())
        self.assertTrue(mon.warned)

    def test_no_warning_when_stoner_below_one(self):
        recs, _ = self._records([1.0, 2.0, 3.0, 4.0], [0.99] * 4)
        self.assertFalse([r for r in recs if r.levelno == logging.WARNING])

    def test_no_warning_when_steps_shrink(self):
        recs, _ = self._records([1.0, 0.4, 0.3, 0.2], [1.5] * 4)
        self.assertFalse([r for r in recs if r.levelno == logging.WARNING])

    def test_warns_on_non_finite_step(self):
        recs, _ = self._records([1.0, float("nan")], [0.5, 0.5])
        self.assertEqual(len([r for r in recs if r.levelno == logging.WARNING]), 1)


class TestStonerLog(unittest.TestCase):

    def _first_line(self, solver, gi):
        with self.assertLogs(_LOGGER, level="INFO") as cm:
            solver.solve(gi, 'tests/flex/output')
        lines = [r.getMessage() for r in cm.records if "Stoner factor" in r.getMessage()]
        self.assertTrue(lines)
        mt = _STONER_RE.search(lines[0])
        self.assertIsNotNone(mt, lines[0])
        q = tuple(float(x) for x in mt.group(2).split(","))
        return float(mt.group(1)), q, float(mt.group(3)), lines

    def test_spin_value_matches_independent_bubble(self):
        """One-orbital Hubbard, one cold-start iteration: Sigma = 0, so the
        logged spin value is U max_q chi0(q, 0) of the stored bubble and the
        charge value is -U min_q chi0(q, 0). Both are checked against the
        stored chi0q, and the bubble itself against the Lindhard sum."""
        U, L, T = 2.0, 4, 0.5
        s, gi = _make_solver({'IterationMax': 1, 'mu': 0.0, 'Ncond': None},
                             U=U, L=L, T=T, Nmat=512)
        spin, q, charge, _ = self._first_line(s, gi)
        chi0 = np.asarray(gi["chi0q"])            # (nmat, nvol, nd, nd)
        nmat = chi0.shape[0]
        c00 = chi0[nmat // 2, :, 0, 0]            # spin-free: one orbital
        self.assertAlmostEqual(spin, U * float(np.max(c00.real)), places=5)
        self.assertAlmostEqual(charge, -U * float(np.min(c00.real)), places=5)
        self.assertEqual(q, (0.5, 0.5, 0.0))
        # Lindhard: chi0(q, 0) = -(1/N) sum_k [f(e_k) - f(e_k+q)] / (e_k - e_k+q)
        k = 2 * np.pi * np.arange(L) / L
        kx, ky = np.meshgrid(k, k, indexing="ij")
        e = np.asarray(s.H0_eigenvalue).reshape(L, L)
        f = 1.0 / (np.exp(e / T) + 1.0)
        qx, qy = L // 2, L // 2
        eq = np.roll(np.roll(e, -qx, 0), -qy, 1)
        fq = np.roll(np.roll(f, -qx, 0), -qy, 1)
        de = e - eq
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.where(np.abs(de) > 1e-12, (f - fq) / de, -f * (1 - f) / T)
        lindhard = -r.mean()
        self.assertAlmostEqual(spin, U * lindhard, delta=2e-3 * U * lindhard)

    def test_general_scheme_logs_same_value_for_one_orbital(self):
        sr, gr = _make_solver({'IterationMax': 1}, Nmat=64)
        sg, gg = _make_solver({'IterationMax': 1}, Nmat=64, calc_scheme="general")
        a = self._first_line(sr, gr)
        b = self._first_line(sg, gg)
        self.assertAlmostEqual(a[0], b[0], places=6)
        self.assertAlmostEqual(a[2], b[2], places=6)

    def test_info_lines_first_and_last_only(self):
        s, gi = _make_solver({'IterationMax': 60}, Nmat=64)
        _, _, _, lines = self._first_line(s, gi)
        self.assertTrue(s.scf_converged)
        self.assertGreater(s.scf_iterations, 2)
        self.assertEqual(len(lines), 2)


class TestStonerLogSpinful(unittest.TestCase):

    def test_spinful_equals_spin_free_channels_in_spin_conserving_limit(self):
        """The spin-orbital (spinful) FLEX has one channel with the
        antisymmetrized vertex; for a spin-diagonal model its static
        spectrum is the union of the spin-free spin and charge channels, so
        the logged spinful value is their maximum (first, Sigma = 0 map)."""
        import shutil
        import tempfile
        import tests.test_flex_spin_orbital as so_mod
        vals = {}
        for so in (False, True):
            d = tempfile.mkdtemp(prefix="flex_282_")
            self.addCleanup(shutil.rmtree, d, True)
            s, gi = so_mod._build(d, so_mod._write_model(d, so), so=so, nmat=16,
                                  IterationMax=1)
            with self.assertLogs(_LOGGER, level="INFO") as cm:
                s.solve(gi, d)
            line = [r.getMessage() for r in cm.records
                    if "Stoner factor" in r.getMessage()][0]
            vals[so] = {k: float(v) for k, v in
                        re.findall(r"(spin|charge|spinful) (\S+) at q", line)}
        self.assertEqual(set(vals[True]), {"spinful"})
        self.assertAlmostEqual(vals[True]["spinful"],
                               max(vals[False]["spin"], vals[False]["charge"]),
                               places=5)


class TestBitIdentity(unittest.TestCase):

    def test_outputs_unchanged_by_diagnostics(self):
        """A converging fixture gives bit-identical outputs with the Stoner
        computation and the monitor disabled."""
        import hwave.solver.flex as solver_flex
        s1, g1 = _make_solver({'IterationMax': 60}, Nmat=64)
        s1.solve(g1, 'tests/flex/output')

        orig = solver_flex.FLEX._stoner_factors
        orig_rec = diag.ScfMonitor.record
        try:
            solver_flex.FLEX._stoner_factors = lambda self, chi0q, ch: None
            diag.ScfMonitor.record = lambda self, *a, **k: None
            s2, g2 = _make_solver({'IterationMax': 60}, Nmat=64)
            s2.solve(g2, 'tests/flex/output')
        finally:
            solver_flex.FLEX._stoner_factors = orig
            diag.ScfMonitor.record = orig_rec
        self.assertEqual(s1.scf_iterations, s2.scf_iterations)
        for key in ("sigma", "green", "chi0q", "chiq_s", "chiq_c"):
            self.assertTrue(np.array_equal(g1[key], g2[key]), key)
        self.assertEqual(s1.mu, s2.mu)


if __name__ == "__main__":
    unittest.main()
