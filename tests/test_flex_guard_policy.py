"""The warn-and-continue guard policy of the FLEX Phase B loop (GitHub
issue #199): a transient violation DURING the self-consistency is logged
with its iteration and its location and the iteration continues, the
per-iteration guard diagnostics reach the log, and the final state is
refused all the same."""
import dataclasses
import os
import types
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_OUT = "tests/flex/output"


def _flex(param_extra=None, path=_IN2, inter=None):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if inter is None else inter)
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 1, "Mix": 0.5, "EPS": 8, "mixing_scheme": "linear"}
    par.update(param_extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


def _violating_density(deviation):
    """A stand-in for ``equal_time_density`` whose result violates the
    Hermitian symmetry by ``deviation``: it refuses exactly as the real one
    does when a tolerance is given, and publishes the deviation when
    ``sym_tol`` is ``None``. Building a Green function that violates the
    symmetry THROUGH the solver is not possible from a Hermitian input."""
    from hwave.solver import hartree_fock as hf
    real = hf.equal_time_density

    def fake(green_kw, heff_eig, mu, beta, shape, *, sym_tol=1e-8):
        res = real(green_kw, heff_eig, mu, beta, shape, sym_tol=None)
        if sym_tol is not None and deviation > sym_tol:
            raise ValueError(
                "equal_time_density: the equal-time density violates rho_ab(r) = "
                "conj(rho_ba(-r)) (relative deviation {:.3e} > {:.1e}); the Green "
                "function is not the Green function of a Hermitian problem"
                .format(deviation, sym_tol))
        return dataclasses.replace(res, hermitian_deviation=deviation)

    return real, fake


class _DensityStub:
    """The smallest object ``FLEX._phase_b_density_and_hf`` operates on: a
    free two-orbital band, no chemical-potential target, a zero
    Hartree-Fock map."""

    def __init__(self, policy, tol, nvol=2, norb=2):
        import hwave.solver.flex as flex_mod
        self.s = flex_mod.FLEX.__new__(flex_mod.FLEX)
        self.s.lattice = types.SimpleNamespace(shape=(nvol, 1, 1), nvol=nvol)
        self.s.norb = norb
        self.s.calc_mu = False
        self.s.H0_k = np.zeros((1, nvol, norb, norb), dtype=complex)
        self.s.flex_guard_policy = policy
        self.s.flex_hf_density_tol = tol
        self.s.flex_guard_violations = 0
        self.s._hf_tables = None
        self.static = np.zeros((1, 1, nvol, norb, norb), dtype=complex)
        beta, nmat = 2.0, 8
        iw = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
        eye = np.eye(norb)
        self.green = np.linalg.inv(
            (iw[:, None, None, None]) * eye + 0.0 * self.s.H0_k[0][None])[None]
        self.beta = beta


class TestPhaseBDensityGuardPolicy(unittest.TestCase):

    def _run(self, policy, deviation, tol=1e-8):
        import hwave.solver.flex_hf as flex_hf
        from hwave.solver import hartree_fock as hf
        st = _DensityStub(policy, tol)
        real_density, fake = _violating_density(deviation)
        real_map = flex_hf.hf_map
        flex_hf.hf_map = lambda rho_r, tables, shape, norb: np.zeros(
            (rho_r.shape[0], norb, norb), dtype=complex)
        hf.equal_time_density = fake
        try:
            return st, st.s._phase_b_density_and_hf(st.green, st.static, 0.0, st.beta, 0.0,
                                                    iteration=7)
        finally:
            hf.equal_time_density = real_density
            flex_hf.hf_map = real_map

    def test_warn_logs_the_iteration_and_continues(self):
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            st, (dens, sigma_hf, _heff) = self._run("warn", 1.0e-3)
        msg = "\n".join(cm.output)
        self.assertIn("FLEX iteration 7", msg)
        self.assertIn("rho_ab(r) = conj(rho_ba(-r))", msg)
        self.assertIn("1.000e-03", msg)
        self.assertIn("flex_hf_density_tol", msg)
        self.assertIn('flex_guard_policy = "warn"', msg)
        self.assertEqual(st.s.flex_guard_violations, 1)
        self.assertEqual(dens.hermitian_deviation, 1.0e-3)
        self.assertTrue(np.all(np.isfinite(sigma_hf)))

    def test_warn_below_the_tolerance_is_silent(self):
        import logging
        records = []
        logger = logging.getLogger("hwave.solver.flex")
        h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
        logger.addHandler(h)
        try:
            st, _ = self._run("warn", 1.0e-12)
        finally:
            logger.removeHandler(h)
        self.assertEqual(st.s.flex_guard_violations, 0)
        self.assertFalse([m for m in records if "rho_ba" in m])

    def test_refuse_still_raises(self):
        with self.assertRaises(ValueError) as cm:
            self._run("refuse", 1.0e-3)
        self.assertIn("rho_ab(r)", str(cm.exception))

    def test_the_warning_formats_without_an_iteration(self):
        """``iteration`` is optional, so the message must format when it is
        absent -- a formatting error would be swallowed by logging and
        printed to stderr as "--- Logging error ---"."""
        import contextlib
        import io
        import hwave.solver.flex_hf as flex_hf
        from hwave.solver import hartree_fock as hf
        st = _DensityStub("warn", 1e-8)
        real_density, fake = _violating_density(1.0e-3)
        real_map = flex_hf.hf_map
        flex_hf.hf_map = lambda rho_r, tables, shape, norb: np.zeros(
            (rho_r.shape[0], norb, norb), dtype=complex)
        hf.equal_time_density = fake
        err = io.StringIO()
        try:
            with contextlib.redirect_stderr(err):
                with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                    st.s._phase_b_density_and_hf(st.green, st.static, 0.0, st.beta, 0.0)
        finally:
            hf.equal_time_density = real_density
            flex_hf.hf_map = real_map
        self.assertEqual(err.getvalue(), "")
        self.assertIn("FLEX iteration None", "\n".join(cm.output))
        self.assertEqual(st.s.flex_guard_violations, 1)

    def test_a_raised_tolerance_lets_the_deviation_through_under_refuse(self):
        st, (dens, _sig, _h) = self._run("refuse", 1.0e-3, tol=1.0e-2)
        self.assertEqual(dens.hermitian_deviation, 1.0e-3)
        self.assertEqual(st.s.flex_guard_violations, 0)


def _dress_result(guard_violations, cond_s=1.0e-5, cond_c=2.0e-5, residual=0):
    from hwave.solver.flex_bond import DressResult
    z = np.zeros((1, 1, 1, 1), dtype=complex)
    return DressResult(collapse0=z, collapse_s=z, collapse_c=z,
                       static_s=z[0], static_c=z[0],
                       cond_min_s=cond_s, cond_min_c=cond_c,
                       guard_violations=guard_violations,
                       guard_residual_violations=residual)


class TestFinalBondGuardRefusal(unittest.TestCase):

    def _solver(self, last):
        import hwave.solver.flex as flex_mod
        s = flex_mod.FLEX.__new__(flex_mod.FLEX)
        s.longitudinal_bond_cond_tol = 1.0e-3
        s.flex_guard_policy = "warn"
        s._bond_last = last
        return s

    def test_a_violated_final_map_is_refused(self):
        s = self._solver(_dress_result(3))
        with self.assertRaises(ValueError) as cm:
            s._refuse_final_bond_guard()
        msg = str(cm.exception)
        self.assertIn("3 time(s)", msg)
        self.assertIn("3 conditioning, 0 residual", msg)
        self.assertIn("longitudinal_bond_cond_tol = 1.0e-03", msg)
        self.assertIn("1.000e-05", msg)
        self.assertIn("2.000e-05", msg)
        self.assertIn("instability region", msg)

    def test_a_residual_only_map_is_named_as_such(self):
        """The reduced "static" guard judges the unchecked slices by their
        solve residual; a map that tripped only that must not be reported as
        a conditioning violation."""
        s = self._solver(_dress_result(1, residual=1))
        with self.assertRaises(ValueError) as cm:
            s._refuse_final_bond_guard()
        msg = str(cm.exception)
        self.assertIn("1 time(s)", msg)
        self.assertIn("0 conditioning, 1 residual", msg)

    def test_a_mixed_map_counts_both_kinds(self):
        s = self._solver(_dress_result(5, residual=2))
        with self.assertRaises(ValueError) as cm:
            s._refuse_final_bond_guard()
        self.assertIn("3 conditioning, 2 residual", str(cm.exception))

    def test_a_clean_final_map_passes(self):
        self._solver(_dress_result(0))._refuse_final_bond_guard()
        self._solver(None)._refuse_final_bond_guard()


class TestGuardDiagnosticsInTheLog(unittest.TestCase):

    def test_hartree_fock_run_logs_the_density_guard_line(self):
        s, r = _flex({"flex_hartree_fock": True})
        os.makedirs(_OUT, exist_ok=True)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(r.get_param("green"), _OUT)
        lines = [m for m in cm.output if "guards: density hermitian" in m]
        self.assertEqual(len(lines), 1)
        self.assertIn("(tol 1.0e-08)", lines[0])
        self.assertNotIn("bond cond_min", lines[0])
        # the residuals line is untouched
        self.assertTrue(any("residuals: sigma" in m and "guards" not in m for m in cm.output))

    def test_bond_gate_run_appends_the_conditioning_diagnostics(self):
        s, r = _flex({"flex_hartree_fock": True, "longitudinal_bond_channels": True,
                      "Nmat": 8})
        os.makedirs(_OUT, exist_ok=True)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(r.get_param("green"), _OUT)
        lines = [m for m in cm.output if "guards: density hermitian" in m]
        self.assertEqual(len(lines), 1)
        self.assertIn("bond cond_min spin", lines[0])
        self.assertIn("charge", lines[0])
        self.assertIn("(tol 1.0e-03)", lines[0])

    def test_no_guard_line_without_phase_b(self):
        s, r = _flex({})
        os.makedirs(_OUT, exist_ok=True)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(r.get_param("green"), _OUT)
        self.assertFalse([m for m in cm.output if "guards:" in m])


class TestFinalStateIsNotCoveredByThePolicy(unittest.TestCase):

    def test_the_final_state_density_is_refused_under_warn(self):
        """A solve whose density keeps violating the symmetry ends in a
        refusal even under ``"warn"``: the map is tolerated (and counted
        WHILE the solve runs), the final state is not.

        After the refusal the counter reads 0 again -- a failed solve drops
        every solve-produced member, this provenance one included, so that a
        later reader cannot mistake it for the record of a result."""
        from hwave.solver import hartree_fock as hf
        s, r = _flex({"flex_hartree_fock": True, "flex_guard_policy": "warn",
                      "IterationMax": 1})
        os.makedirs(_OUT, exist_ok=True)
        gi = r.get_param("green")
        seen = []
        s._iteration_hook = lambda d: seen.append(s.flex_guard_violations)
        real_density, fake = _violating_density(1.0e-3)
        hf.equal_time_density = fake
        try:
            with self.assertRaises(ValueError) as cm:
                s.solve(gi, _OUT)
        finally:
            hf.equal_time_density = real_density
        self.assertIn("rho_ab(r)", str(cm.exception))
        self.assertEqual(seen, [1])                           # the map was tolerated, once
        self.assertEqual(s.flex_guard_violations, 0)          # and the failure left nothing
        self.assertNotIn("green", gi)                         # nothing produced survives


class TestIterationMaxZero(unittest.TestCase):
    """No map ran, so there is nothing for the policy to tolerate and
    nothing for the final bond check to judge."""

    def test_hartree_fock_only(self):
        s, r = _flex({"flex_hartree_fock": True, "flex_guard_policy": "warn",
                      "IterationMax": 0})
        os.makedirs(_OUT, exist_ok=True)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(r.get_param("green"), _OUT)
        self.assertFalse([m for m in cm.output if "guards:" in m])
        self.assertEqual(s.flex_guard_violations, 0)

    def test_bond_gate(self):
        s, r = _flex({"flex_hartree_fock": True, "longitudinal_bond_channels": True,
                      "flex_guard_policy": "warn", "Nmat": 8, "IterationMax": 0})
        os.makedirs(_OUT, exist_ok=True)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s.solve(r.get_param("green"), _OUT)
        self.assertFalse([m for m in cm.output if "guards:" in m])
        self.assertEqual(s.flex_guard_violations, 0)
        self.assertIsNone(getattr(s, "_bond_last", None))


class TestGuardProvenance(unittest.TestCase):

    def test_provenance_carries_the_policy_and_the_violation_count(self):
        s, r = _flex({"flex_hartree_fock": True, "flex_guard_policy": "warn"})
        os.makedirs(_OUT, exist_ok=True)
        gi = r.get_param("green")
        s.solve(gi, _OUT)
        prov = s._provenance_block("final_state")
        self.assertEqual(prov["flex_guard_policy"], "warn")
        self.assertEqual(prov["flex_guard_violations"], 0)
        self.assertIsInstance(prov["flex_guard_violations"], int)

    def test_the_violation_count_is_reset_at_solve_entry(self):
        s, r = _flex({"flex_hartree_fock": True, "flex_guard_policy": "warn"})
        os.makedirs(_OUT, exist_ok=True)
        s.flex_guard_violations = 11
        s.solve(r.get_param("green"), _OUT)
        self.assertEqual(s.flex_guard_violations, 0)

    def test_the_bond_archive_records_the_conditioning_tolerance(self):
        s, r = _flex({"flex_hartree_fock": True, "longitudinal_bond_channels": True,
                      "Nmat": 8, "longitudinal_bond_cond_tol": 1.0e-6})
        os.makedirs(_OUT, exist_ok=True)
        gi = r.get_param("green")
        s.solve(gi, _OUT)
        self.assertEqual(float(gi["longitudinal_bond_cond_tol"]), 1.0e-6)


if __name__ == "__main__":
    unittest.main()
