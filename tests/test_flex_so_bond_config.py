"""Configuration and admission of the spin-orbital bond-resolved FLEX
(issue #270, gate C1b-R): ``[mode] enable_spin_orbital = true`` with
``[mode.param] longitudinal_bond_channels = true``.

* :class:`TestFixtures` -- the complex-PairHop fixtures compile and select
  their layouts.
* :class:`TestRefusals` -- every refusal, its exact text, with each key and
  each string value also spelled in mixed case.
* :class:`TestWarningsAndInfo` -- every WARNING / INFO line, logged once.
* :class:`TestCapParsing` -- the memory cap: the non-spin-orbital semantics
  unchanged, the spin-orbital ones (unset = 0.8 x available, 0 disables,
  < 0 refused, no probe).
* :class:`TestAdmission` -- ``FLEX._phase_b_preflight`` builds the
  admission: topology, layout, sector Dyson flag, the over-cap refusal and
  the frequency-batch selection.

The texts are literals here (not read from the solver), so a changed text
fails this module.
"""
import cmath
import contextlib
import logging
import os
import shutil
import tempfile
import tracemalloc
import unittest
from unittest import mock

import numpy as np

from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X

H = ("[mode] enable_spin_orbital = true with [mode.param] longitudinal_bond_channels = true "
     "(FLEX)")
CAP_KEY = "[mode.param] longitudinal_bond_memory_cap_gb"

T_OUTPUT_FULL = (H + " does not write the bond archive in this version (the spin-orbital bond "
                 "archive arrives in a later version); remove longitudinal_bond_output_full "
                 "or set it to false.")
T_OUTPUT_LAYOUT = (H + " does not write the bond archive in this version (the spin-orbital "
                   "bond archive arrives in a later version); remove "
                   "longitudinal_bond_output_layout.")
T_PAIRING = (H + " does not run the in-process pairing in this version (it arrives in a later "
             "version); remove longitudinal_bond_pairing.")
T_STORE = (H + " runs on the host only in this version (the device store arrives in a later "
           "version); set longitudinal_bond_store = \"host\" or \"auto\", or remove it.")
T_GUARD_FREQS = (H + " guards every bosonic frequency (the spin-orbital dressing decomposes "
                 "every block exactly); remove longitudinal_bond_guard_freqs or set it to "
                 "\"all\".")
T_GPU_REQUIRED = (H + " runs on the host only in this version (the device path arrives in a "
                  "later version); remove gpu_required or set it to false.")
T_GPU = H + " runs on the host only in this version; gpu = true is ignored for this run."
T_GUARD_POLICY = (H + ": flex_guard_policy = \"warn\" relaxes the density guard only; the bond "
                  "dressing guard refuses a near-singular block (lower "
                  "longitudinal_bond_cond_tol to admit it).")
T_GUARD_EXACT = (H + ": longitudinal_bond_guard_exact has no effect (the spin-orbital dressing "
                 "decomposes every block exactly).")
T_LAYOUT_IGNORED = ("[mode.param] so_bond_layout is used only with enable_spin_orbital = true "
                    "and longitudinal_bond_channels = true; ignored.")
T_LAYOUT_VALUE = "[mode.param] so_bond_layout must be \"auto\" or \"dense\", got {!r}"
T_LAYOUT_FORCED = ("[mode.param] so_bond_layout = \"sz_packed\" cannot be forced: the S_z-packed "
                   "layout is chosen automatically when the input conserves S_z (no spin-flip "
                   "Transfer entry, no S_z-changing interaction). Use \"auto\" (the default) or "
                   "\"dense\".")
T_GATE_OFF_OFFSITE = (
    "[mode] enable_spin_orbital = true (FLEX) supports on-site and in-cell interactions only "
    "without the bond channels (the antisymmetrized spinful vertex carries the on-site "
    "exchange crossing only); the off-site 'CoulombInter' declaration at irvec=(1, 0, 0), "
    "orbvec=(0, 0) is refused; set [mode.param] longitudinal_bond_channels = true to resolve "
    "their exchange exactly.")
T_B1 = ("spin-orbital bond FLEX: no off-site row; the bond machinery runs with the on-site "
        "channel only (B = 1).")
T_KEY_NAMING = ("longitudinal_bond_channels with enable_spin_orbital: bond channels of the "
                "single spinful channel; 'longitudinal' is historical")
T_PAIRLIFT = ("off-site PairLift is included in the spin-orbital bond vertex (the "
              "non-spin-orbital bond paths omit it: it has no spin/charge content)")
T_NONSO_PAIRLIFT = "are ignored by the bond channels"
T_NB = ("spin-orbital bond FLEX: frequency batch nb = {} (the largest of 16, 8, 4, 2, 1 not "
        "above Nmat that fits the cap)")
T_STORE_INFO = ("bond store: W is written over chibar in place (no later reader of chibar; "
                "chibar_w archive: none)")
T_CAP_RANGE = ("[mode.param] longitudinal_bond_memory_cap_gb must be a finite number >= 0 in "
               "spin-orbital mode (0 disables the check), got -1")
T_CAP_ZERO = ("[mode.param] longitudinal_bond_memory_cap_gb = 0: the memory admission check of "
              "the spin-orbital bond path is disabled")
T_CAP_NO_PROBE = ("no memory limit could be determined; the admission check is disabled; set "
                  "[mode.param] longitudinal_bond_memory_cap_gb to enable it")
T_NONSO_CAP = ("[mode.param] longitudinal_bond_memory_cap_gb must be a finite number > 0 "
               "(binary GiB), got {}")
T_SIDECAR = "needs longitudinal_bond_output_full = true"
T_PACKED_REASON = ("so bond layout: sz_packed (S_z conserved: no spin-flip transfer, no "
                   "S_z-changing interaction)")
T_SOC_REASON = "so bond layout: dense (reason: spin-flip transfer entry R=(-1, 0, 0) (1, 2) = " \
               "(-0.3+0j))"
T_OVERCAP_HEAD = "spin-orbital bond FLEX needs an estimated peak of"
T_REMEDY_AUTO = ("The frequency batch was selected automatically; the smallest batch (nb = 1) "
                 "was tried.")
T_REMEDY_SET = ("longitudinal_bond_freq_batch = 4 is set; a smaller value, or removing it so "
                "that the run chooses, lowers the dressing rows.")
T_REMEDY_SET_1 = "longitudinal_bond_freq_batch = 1 is set (the smallest batch)."

#: a CoulombInter coefficient with an imaginary part (Hermitian-closed)
_COMPLEX_INTER = X.ZEEMAN_SZ_ALL._replace(
    name="ZEEMAN_COMPLEX_INTER", interactions=tuple(
        (t, X._nn_c(0.37 * cmath.exp(0.3j)) if t == "CoulombInter" else rows)
        for t, rows in X.ZEEMAN_SZ_ALL.interactions))


def _mixed(key):
    """``key`` in mixed case (every other letter upper)."""
    return "".join(c.upper() if i % 2 == 0 else c for i, c in enumerate(key))


class _Case(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="flex_so_bond_cfg_")
        self.addCleanup(shutil.rmtree, self.d, True)
        self._n = 0

    def _dir(self):
        self._n += 1
        d = os.path.join(self.d, "c%d" % self._n)
        os.makedirs(d)
        return d

    def build(self, fixture=X.ZEEMAN_SZ_ALL, *, gate=True, so=True, param=None, geom=None):
        d = self._dir()
        input_dict = X.write_flex_case(d, fixture, gate=gate, so=so, param=param)
        if geom is not None:
            with open(os.path.join(d, "geom.dat"), "w") as f:
                f.write(geom)
        return X.build_solver(input_dict)

    @staticmethod
    def preflight(solver, green_info):
        solver._phase_b_reset(green_info)
        solver._phase_b_preflight(green_info)
        return solver._so_admission

    @contextlib.contextmanager
    def capture(self):
        """The messages of every ``hwave`` record (INFO and above) logged inside."""
        msgs = []

        class _H(logging.Handler):
            def emit(self, record):
                msgs.append((record.levelname, record.getMessage()))
        lg = logging.getLogger("hwave")
        h = _H(level=logging.INFO)
        level = lg.level
        lg.addHandler(h)
        lg.setLevel(logging.INFO)
        try:
            yield msgs
        finally:
            lg.removeHandler(h)
            lg.setLevel(level)

    def run_logs(self, fixture=X.ZEEMAN_SZ_ALL, *, preflight=True, **kw):
        with self.capture() as msgs:
            s, gi = self.build(fixture, **kw)
            if preflight:
                self.preflight(s, gi)
        return s, msgs

    def assertOnce(self, msgs, level, text):
        self.assertEqual(msgs.count((level, text)), 1, [m for m in msgs if text[:30] in m[1]])


class TestFixtures(unittest.TestCase):
    def test_complex_pairhop(self):
        from hwave.solver import so_bond
        from hwave.solver.so_records import compile_interaction
        for fx in (X.ZEEMAN_SZ_ALL, X.SOC_ALL_TYPES_CPH):
            with self.subTest(fixture=fx.name):
                tables = X.interaction_tables(fx)
                ph = tables["PairHop"]
                # the Hermitian closure: the -R partner of every row is its conjugate
                for (R, ab), v in ph.items():
                    self.assertEqual(ph[(tuple(-x for x in R), ab[::-1])], complex(v).conjugate())
                self.assertGreater(abs(complex(X.V_CPH).imag), 0.1)
                compiled = compile_interaction(tables, fx.norb_phys)
                rec = [r for r in compiled.records if r["type"] == "PairHop"]
                self.assertTrue(rec)
                self.assertGreater(max(abs(complex(r["v"]).imag) for r in rec), 0.1)
                if fx is X.ZEEMAN_SZ_ALL:
                    self.assertIsNone(so_bond.records_conserve_sz(compiled))
                    self.assertEqual(so_bond.select_layout(compiled, {}, override="auto")[0],
                                     "sz_packed")
        # the replaced rows only
        self.assertEqual(
            [t for t, _ in X.SOC_ALL_TYPES_CPH.interactions],
            [t for t, _ in F.SOC_ALL_TYPES.interactions])
        for (t, a), (_, b) in zip(X.SOC_ALL_TYPES_CPH.interactions,
                                  F.SOC_ALL_TYPES.interactions):
            if t != "PairHop":
                self.assertEqual(a, b)

    def test_onsite_fixture_moved(self):
        import tests.test_sc_spin_orbital_bond as sc_tests
        self.assertIs(sc_tests.ZEEMAN_EXTERN_ONSITE, X.ZEEMAN_EXTERN_ONSITE)
        self.assertTrue(all(r[0] == (0, 0, 0) for _, rows in
                            X.ZEEMAN_EXTERN_ONSITE.interactions for r in rows))

    def test_write_flex_case(self):
        with tempfile.TemporaryDirectory() as d:
            idict = X.write_flex_case(d, X.ZEEMAN_SZ_ALL)
        self.assertNotIn("eliashberg", idict)
        mode = idict["mode"]
        self.assertEqual((mode["mode"], mode["enable_spin_orbital"], mode["calc_scheme"]),
                         ("FLEX", True, "general"))
        par = {k.lower(): v for k, v in mode["param"].items()}
        self.assertNotIn("filling", par)
        self.assertEqual(par["mu"], X.fixture_mu(X.ZEEMAN_SZ_ALL))
        for k, v in (("flex_hartree_fock", True), ("longitudinal_bond_channels", True),
                     ("nmat", 16), ("t", 0.25), ("coeff_tail", 0.0)):
            self.assertEqual(par[k], v, k)


class TestRefusals(_Case):
    """Every refusal of the table, exact text; each key (and string value)
    also in mixed case."""

    def assertRefused(self, text, *, preflight=False, fixture=X.ZEEMAN_SZ_ALL, **kw):
        with self.assertRaises(ValueError) as cm:
            s, gi = self.build(fixture, **kw)
            if preflight:
                self.preflight(s, gi)
        self.assertEqual(str(cm.exception), text)

    def _rows(self):
        return (
            (T_OUTPUT_FULL, {"longitudinal_bond_output_full": True}),
            (T_OUTPUT_LAYOUT, {"longitudinal_bond_output_layout": "sidecar"}),
            (T_OUTPUT_LAYOUT, {"longitudinal_bond_output_layout": "bogus"}),
            (T_PAIRING, {"longitudinal_bond_pairing": "singlet"}),
            (T_PAIRING, {"longitudinal_bond_pairing": "both"}),
            (T_STORE, {"gpu": True, "longitudinal_bond_store": "device"}),
            (T_STORE, {"longitudinal_bond_store": "device"}),
            (T_STORE, {"gpu": False, "longitudinal_bond_store": "device"}),
            (T_GUARD_FREQS, {"longitudinal_bond_guard_freqs": "static"}),
            (T_GPU_REQUIRED, {"gpu_required": True}),
            (T_LAYOUT_VALUE.format("bogus"), {"so_bond_layout": "bogus"}),
            (T_LAYOUT_FORCED, {"so_bond_layout": "sz_packed"}),
        )

    def test_each_row(self):
        for text, param in self._rows():
            mixed = {_mixed(k): (_mixed(v) if isinstance(v, str) else v)
                     for k, v in param.items()}
            for p in (param, mixed):
                exp = text
                if "so_bond_layout" in [k.lower() for k in p] and text.startswith(
                        "[mode.param] so_bond_layout must"):
                    exp = T_LAYOUT_VALUE.format([v for v in p.values()][0])
                with self.subTest(param=p):
                    self.assertRefused(exp, param=p)

    def test_gate_off_offsite(self):
        self.assertRefused(T_GATE_OFF_OFFSITE, preflight=True, fixture=F.SOC_ALL_TYPES,
                           gate=False)

    def test_complex_non_pairhop(self):
        with self.assertRaises(ValueError) as cm:
            self.build(_COMPLEX_INTER)
        self.assertIn("imaginary part", str(cm.exception))
        self.assertIn("CoulombInter", str(cm.exception))

    def test_cap_range(self):
        for key in ("longitudinal_bond_memory_cap_gb", _mixed("longitudinal_bond_memory_cap_gb")):
            with self.subTest(key=key):
                self.assertRefused(T_CAP_RANGE, preflight=True, param={key: -1})

    def test_accepted_values(self):
        """The values the table accepts construct and admit."""
        for p in ({"longitudinal_bond_output_full": False},
                  {"Longitudinal_Bond_Output_Full": False},
                  {"longitudinal_bond_store": "host"}, {"Longitudinal_Bond_Store": "AUTO"},
                  {"longitudinal_bond_guard_freqs": "all"},
                  {"Longitudinal_Bond_Guard_Freqs": " All "},
                  {"longitudinal_bond_pairing": "none"}, {"Longitudinal_Bond_Pairing": "NONE"},
                  {"gpu_required": False}, {"so_bond_layout": "Dense"},
                  {"SO_Bond_Layout": " AUTO "}):
            with self.subTest(param=p):
                s, gi = self.build(param=p)
                self.assertTrue(s._so_bond)
                self.assertIsNotNone(self.preflight(s, gi))

    def test_layout_request_case(self):
        s, gi = self.build(param={"SO_Bond_Layout": "DeNsE"})
        self.assertEqual(s.so_bond_layout, "dense")
        self.assertEqual(self.preflight(s, gi).layout, "dense")
        s, _ = self.build()
        self.assertIsNone(s.so_bond_layout)


class TestWarningsAndInfo(_Case):
    """Every WARNING / INFO line of the table, logged exactly once."""

    def test_each_row(self):
        for key in ("gpu", "GPU"):
            with self.subTest(row="gpu", key=key):
                s, msgs = self.run_logs(param={key: True})
                self.assertOnce(msgs, "WARNING", T_GPU)
                self.assertFalse(s.use_gpu)
        for key, value in (("flex_guard_policy", "warn"), ("Flex_Guard_Policy", "WaRn")):
            with self.subTest(row="flex_guard_policy", key=key):
                _, msgs = self.run_logs(param={key: value})
                self.assertOnce(msgs, "WARNING", T_GUARD_POLICY)
        for key, value in (("longitudinal_bond_guard_exact", "always"),
                           ("Longitudinal_Bond_Guard_Exact", "Decide")):
            with self.subTest(row="guard_exact", key=key):
                _, msgs = self.run_logs(param={key: value})
                self.assertOnce(msgs, "INFO", T_GUARD_EXACT)
        for fx, kw in ((X.ZEEMAN_EXTERN_ONSITE, {"gate": False}),
                       (F.SU2_4X4, {"so": False})):
            for key in ("so_bond_layout", "SO_Bond_Layout"):
                with self.subTest(row="layout ignored", fixture=fx.name, key=key):
                    _, msgs = self.run_logs(fx, param={key: "dense"}, **kw)
                    self.assertOnce(msgs, "WARNING", T_LAYOUT_IGNORED)
        _, msgs = self.run_logs(param={"so_bond_layout": "dense"})
        self.assertNotIn(("WARNING", T_LAYOUT_IGNORED), msgs)
        with self.subTest(row="B = 1"):
            _, msgs = self.run_logs(X.ZEEMAN_EXTERN_ONSITE)
            self.assertOnce(msgs, "INFO", T_B1)
        _, msgs = self.run_logs()
        self.assertNotIn(("INFO", T_B1), msgs)
        for row, text in (("key naming", T_KEY_NAMING), ("layout", T_PACKED_REASON),
                          ("nb", T_NB.format(16)), ("store", T_STORE_INFO)):
            with self.subTest(row=row):
                self.assertOnce(msgs, "INFO", text)
        self.assertNotIn(("INFO", T_PAIRLIFT), msgs)   # ZEEMAN_SZ_ALL has no PairLift
        self.assertFalse(any("Bond-resolved FLEX preflight (ESTIMATE)" in m for _, m in msgs))

    def test_no_sidecar_warning(self):
        """The layout key is refused by its presence, before the non-SO value
        check (whose sidecar WARNING is not logged)."""
        with self.capture() as msgs:
            with self.assertRaises(ValueError):
                self.build(param={"longitudinal_bond_output_layout": "sidecar"})
        self.assertFalse(any(T_SIDECAR in m for _, m in msgs), msgs)

    def test_no_nonso_pairlift_info(self):
        _, msgs = self.run_logs(F.SOC_ALL_TYPES)
        self.assertOnce(msgs, "INFO", T_PAIRLIFT)
        self.assertOnce(msgs, "INFO", T_SOC_REASON)
        self.assertFalse(any(T_NONSO_PAIRLIFT in m for _, m in msgs), msgs)


class TestCapParsing(_Case):
    @staticmethod
    def _probe(nbytes):
        from hwave.solver import host_memory
        return mock.patch.object(host_memory, "available_host_bytes",
                                 return_value=host_memory.HostMemory(
                                     nbytes, "none" if nbytes is None else "test"))

    def test_non_so(self):
        s, _ = self.build(F.SU2_4X4, so=False)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb, 8.0)
        for v in (-1, 0):
            for key in ("longitudinal_bond_memory_cap_gb",
                        _mixed("longitudinal_bond_memory_cap_gb")):
                with self.subTest(value=v, key=key):
                    with self.assertRaises(ValueError) as cm:
                        self.build(F.SU2_4X4, so=False, param={key: v})
                    self.assertEqual(str(cm.exception), T_NONSO_CAP.format(v))

    def test_so(self):
        avail = 10 * 2 ** 30
        s, gi = self.build()
        self.assertIsNone(s.longitudinal_bond_memory_cap_gb)
        with self._probe(avail):
            adm = self.preflight(s, gi)
        self.assertEqual(adm.caps.host_bytes, 0.8 * avail)
        for key in ("longitudinal_bond_memory_cap_gb", _mixed("longitudinal_bond_memory_cap_gb")):
            with self.subTest(value=0, key=key):
                with self.capture() as msgs:
                    s, gi = self.build(param={key: 0})
                    adm = self.preflight(s, gi)
                self.assertOnce(msgs, "WARNING", T_CAP_ZERO)
                self.assertIsNone(adm.caps.host_bytes)
            with self.subTest(value=2, key=key):
                s, gi = self.build(param={key: 2})
                self.assertEqual(self.preflight(s, gi).caps.host_bytes, 2.0 * 2 ** 30)
            with self.subTest(value=-1, key=key):
                s, gi = self.build(param={key: -1})       # kept raw at construction
                self.assertEqual(s.longitudinal_bond_memory_cap_gb, -1)
                with self.assertRaises(ValueError) as cm:
                    self.preflight(s, gi)
                self.assertEqual(str(cm.exception), T_CAP_RANGE)
        with self.subTest(probe=None):
            s, gi = self.build()
            with self.capture() as msgs, self._probe(None):
                adm = self.preflight(s, gi)
            self.assertOnce(msgs, "WARNING", T_CAP_NO_PROBE)
            self.assertIsNone(adm.caps.host_bytes)


class TestAdmission(_Case):
    def test_topology_and_layout(self):
        from hwave.solver import so_bond
        from hwave.solver.so_records import compile_interaction
        for fx, layout, reason, dyson in (
                (F.SOC_ALL_TYPES, "dense", T_SOC_REASON, False),
                (X.ZEEMAN_SZ_ALL, "sz_packed", T_PACKED_REASON, True)):
            with self.subTest(fixture=fx.name):
                s, gi = self.build(fx)
                self.assertTrue(s._so_bond)
                self.assertIsNone(s._so_admission)
                adm = self.preflight(s, gi)
                self.assertIs(s._so_admission, adm)
                ref = so_bond.so_bond_topology(
                    compile_interaction(X.interaction_tables(fx), fx.norb_phys), np.eye(3))
                np.testing.assert_array_equal(adm.topo.delta_r, ref.delta_r)
                np.testing.assert_array_equal(adm.topo.reverse, ref.reverse)
                self.assertEqual(adm.geom.B, 5)
                self.assertEqual((adm.layout, adm.reason), (layout, reason))
                self.assertEqual(adm.sector_map is None, layout == "dense")
                self.assertIs(adm.sector_dyson, dyson)
                self.assertEqual(adm.plan.layout, layout)
                self.assertIn("flex_scf", adm.plan.consumers)
                self.assertEqual(adm.plan.aliases, (("W", "chibar"),))
                self.assertEqual((adm.nb, adm.plan.nb), (16, 16))
                np.testing.assert_array_equal(adm.bond_inputs.delta_r, ref.delta_r)

    def test_dense_override(self):
        s, gi = self.build(param={"so_bond_layout": "dense"})
        adm = self.preflight(s, gi)
        self.assertEqual(adm.layout, "dense")
        self.assertIsNone(adm.sector_map)
        self.assertTrue(adm.sector_dyson)        # the inputs still conserve S_z

    def test_over_cap(self):
        s, gi = self.build(param={"longitudinal_bond_memory_cap_gb": 1e-6})
        with self.assertRaises(ValueError) as cm:
            self.preflight(s, gi)
        msg = str(cm.exception)
        self.assertTrue(msg.startswith(T_OVERCAP_HEAD), msg)
        self.assertIn("(" + CAP_KEY + "). " + T_REMEDY_AUTO + " ", msg)
        for part in ("rows (phase  name  bytes):",
                     "shells (shell  |dr|  max|v|):", "raise " + CAP_KEY,
                     "  resident  collapses  "):
            self.assertIn(part, msg)
        # the admitted-cases table is that of the hwave_sc static path; the
        # FLEX refusal does not append it
        self.assertNotIn("admitted cases", msg)
        self.assertNotIn("hwave_sc static", msg)
        self.assertIsNone(s._so_admission)

    def _spy_plans(self):
        """Patch ``so_bond.plan_so_bond`` with a wrapper that records every
        call (its ``nb``, ``caps``, ``log_summary`` and the plan) and emits
        an unrelated ``so_bond`` record before each call."""
        from hwave.solver import so_bond
        real = so_bond.plan_so_bond
        calls = []

        def spy(**kw):
            so_bond.logger.info("unrelated so_bond record %d", len(calls))
            call = dict(nb=kw.get("nb"), caps=kw.get("caps"),
                        log_summary=kw.get("log_summary", True), plan=None)
            calls.append(call)
            call["plan"] = real(**kw)
            return call["plan"]
        return mock.patch.object(so_bond, "plan_so_bond", side_effect=spy), calls

    def test_over_cap_tries_nb_1(self):
        """The capped plan of an auto-batch over-cap run is the one at the
        smallest batch, nb = 1 (as its remedy sentence says)."""
        patch, calls = self._spy_plans()
        s, gi = self.build(param={"longitudinal_bond_memory_cap_gb": 1e-6})
        with patch as m, self.assertRaises(ValueError) as cm:
            self.preflight(s, gi)
        self.assertTrue(m.called)
        self.assertIn(T_REMEDY_AUTO, str(cm.exception))
        # the trial plans at 16, 8, 4, 2, 1 (unlimited), then the capped one
        self.assertEqual([c["nb"] for c in calls], [16, 8, 4, 2, 1, 1])
        self.assertIsNone(calls[-1]["plan"])            # the capped call refused
        self.assertIsNotNone(calls[-1]["caps"].host_bytes)
        self.assertTrue(all(c["caps"].host_bytes is None for c in calls[:-1]))

    def test_trial_plans_quiet_call_local(self):
        """The trial plans of the batch selection are not logged, the
        admitting plan is (once), and a record another caller emits on the
        ``so_bond`` logger in the meantime is kept: the suppression is a
        per-call option, not a filter on the shared logger."""
        from hwave.solver import so_bond
        s, _ = self.build(param={"longitudinal_bond_memory_cap_gb": 0})
        p8, p4 = self._peaks(s)
        cap = 0.5 * (p8 + p4) / 2 ** 30
        patch, calls = self._spy_plans()
        with patch as m:
            s, msgs = self.run_logs(param={"longitudinal_bond_memory_cap_gb": cap})
        self.assertTrue(m.called)
        self.assertEqual([c["nb"] for c in calls], [16, 8, 4, 4])
        self.assertEqual([c["log_summary"] for c in calls], [False, False, False, True])
        for i in range(len(calls)):
            self.assertOnce(msgs, "INFO", "unrelated so_bond record %d" % i)
        final = calls[-1]["plan"]
        self.assertIs(s._so_admission.plan, final)
        self.assertOnce(msgs, "INFO", final.reason)
        for c in calls[:-1]:
            if c["plan"].reason != final.reason:
                self.assertNotIn(("INFO", c["plan"].reason), msgs)
        self.assertEqual(so_bond.logger.filters, [])

    def test_over_cap_set_batch(self):
        for nb, note in ((4, T_REMEDY_SET), (1, T_REMEDY_SET_1)):
            with self.subTest(nb=nb):
                s, gi = self.build(param={"longitudinal_bond_memory_cap_gb": 1e-6,
                                          "longitudinal_bond_freq_batch": nb})
                s._phase_b_reset(gi)
                tracemalloc.start()
                try:
                    base = tracemalloc.get_traced_memory()[0]
                    tracemalloc.reset_peak()
                    with self.assertRaises(ValueError) as cm:
                        s._phase_b_preflight(gi)
                    peak = tracemalloc.get_traced_memory()[1] - base
                finally:
                    tracemalloc.stop()
                self.assertIn("(" + CAP_KEY + "). " + note + " ", str(cm.exception))
                self.assertLessEqual(peak, 2 ** 20)

    def _peaks(self, s):
        """The plan peaks at nb = 8 and nb = 4 of the run of ``s``, built
        test-side with the admission's rows."""
        from hwave.solver import flex_so_bond, so_bond
        nvol = int(np.prod(s.lattice.shape))
        adm = self.preflight(s, {})
        resident, extra = flex_so_bond.flex_so_bond_rows(
            nmat=s.nmat, nvol=nvol, nd=2, depth=s.anderson_depth, mixing=s.mixing_scheme,
            n_types=sum(1 for t in s._hf_tables.inter_table.values() if t is not None))
        return [so_bond.plan_so_bond(
            nmat=s.nmat, nvol=nvol, B=adm.geom.B, nd=2, consumers=("flex_scf",),
            caps=so_bond.SOBondCaps(None, None, None), nb=nb, layout=adm.layout,
            layout_reason=adm.reason, resident_rows=resident, extra_rows=extra).peak_host
            for nb in (8, 4)]

    def test_batch_selection(self):
        s, _ = self.build(param={"longitudinal_bond_memory_cap_gb": 0})
        p8, p4 = self._peaks(s)
        self.assertGreater(p8, p4)
        cap = 0.5 * (p8 + p4) / 2 ** 30
        _, msgs = self.run_logs(param={"longitudinal_bond_memory_cap_gb": cap})
        self.assertOnce(msgs, "INFO", T_NB.format(4))
        s, gi = self.build(param={"longitudinal_bond_memory_cap_gb": cap})
        adm = self.preflight(s, gi)
        self.assertEqual((adm.nb, adm.plan.nb), (4, 4))
        self.assertLessEqual(adm.plan.peak_host, cap * 2 ** 30)
        # a set batch is used as is, without the selection notice
        s, msgs = self.run_logs(param={"Longitudinal_Bond_Freq_Batch": 8})
        self.assertEqual(s._so_admission.nb, 8)
        self.assertFalse(any("frequency batch nb =" in m for _, m in msgs))

    def test_cond_tol(self):
        s, gi = self.build(param={"Longitudinal_Bond_Cond_Tol": 0.01})
        self.preflight(s, gi)
        self.assertEqual(s.longitudinal_bond_cond_tol, 0.01)

    def test_rectangular_cell(self):
        from hwave.solver import so_bond
        from hwave.solver.so_records import compile_interaction
        # a 4 x 6 mesh, so that the next-nearest neighbours (0, +-2, 0) along
        # y are distinct modulo the mesh
        fx = X.ZEEMAN_SZ_ALL._replace(name="ZEEMAN_SZ_ALL_RECT", shape=(4, 6, 1), interactions=tuple(
            (t, rows + ((((0, 2, 0), (0, 0), 0.05), ((0, -2, 0), (0, 0), 0.05))
                        if t == "CoulombInter" else ()))
            for t, rows in X.ZEEMAN_SZ_ALL.interactions))
        rvec = np.diag([1.0, 2.0, 1.0])
        s, gi = self.build(fx, geom="1 0 0\n0 2 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
        np.testing.assert_array_equal(np.asarray(s.ham_info.param_ham["Geometry"]["rvec"]), rvec)
        adm = self.preflight(s, gi)
        compiled = compile_interaction(X.interaction_tables(fx), 1)
        rect = so_bond.bond_run(compiled, rvec, None)[1]
        square = so_bond.bond_run(compiled, np.eye(3), None)[1]
        for name in ("delta_r", "reverse"):
            np.testing.assert_array_equal(getattr(adm.bond_inputs, name), getattr(rect, name))
        self.assertEqual(adm.bond_inputs.max_shells, rect.max_shells)
        self.assertEqual(np.asarray(rect.delta_r).shape, np.asarray(square.delta_r).shape)
        self.assertFalse(np.array_equal(rect.delta_r, square.delta_r))

    def test_fit_does_not_catch(self):
        from hwave.solver import so_bond
        s, gi = self.build()
        s._phase_b_reset(gi)
        with mock.patch.object(so_bond, "plan_so_bond", side_effect=ValueError("x")) as spy:
            with self.assertRaises(ValueError) as cm:
                s._phase_b_preflight(gi)
        self.assertEqual(spy.call_count, 1)
        self.assertEqual(str(cm.exception), "x")
        self.assertIsNone(s._so_admission)

    def test_non_so_gate_untouched(self):
        """The non-spin-orbital gate: the construction and preflight lines of
        a92b533d (recorded there), one memory estimate, no spin-orbital state."""
        import hwave.solver.flex as flex_mod
        from tests.test_flex_bond_gate import _flex
        recs = []

        class _H(logging.Handler):
            def emit(self, record):
                recs.append((record.levelname, record.getMessage()))
        lg = logging.getLogger("hwave.solver.flex")
        h = _H(level=logging.INFO)
        level = lg.level
        lg.addHandler(h)
        lg.setLevel(logging.INFO)
        orig = flex_mod.FLEX._bond_memory_preflight
        try:
            s, r = _flex({"IterationMax": 1})
            gi = r.get_param("green")
            with mock.patch.object(flex_mod.FLEX, "_bond_memory_preflight", autospec=True,
                                   side_effect=orig) as spy:
                self.preflight(s, gi)
        finally:
            lg.removeHandler(h)
            lg.setLevel(level)
        self.assertEqual(spy.call_count, 1)
        self.assertFalse(s._so_bond)
        self.assertIsNone(s._so_admission)
        self.assertEqual(recs, [("INFO", m) for m in _A92_GATE_LINES])


#: the hwave.solver.flex lines of test_flex_bond_gate._flex({"IterationMax": 1})
#: construction and preflight, recorded at a92b533d
_A92_GATE_LINES = (
    "    flex_second_order = local",
    "    second-order factors: 3 on-site pairs, off-site density, 0.005 MiB",
    "FLEX parameters:",
    "    max_iter        = 1",
    "    mix             = 0.5",
    "    eps             = 1.000000e-12",
    "    sigma Hermitian projection = on",
    "    calc_scheme     = general (explicit)",
    "    flex_hartree_fock = True",
    "    longitudinal_bond_channels (FLEX) = True",
    "    longitudinal_bond_pairing = none",
    "Bond-resolved FLEX preflight (ESTIMATE): B = 5 channels [(0, 0, 0), (-1, 0, 0), "
    "(0, -1, 0), (0, 1, 0), (1, 0, 0)], ND = 20, nvol = 16, Nmat = 8, frequency batch 8; "
    "persistent 0.0014 GiB, peak 0.0075 GiB = 1.25 * (persistent + max phase row) against the "
    "cap 8.0000 GiB\n"
    "  persistent           chibar_W:     0.0008 GiB\n"
    "  persistent           chi_sc_w:     0.0000 GiB\n"
    "  persistent    vertices_static:     0.0005 GiB\n"
    "  persistent              state:     0.0000 GiB\n"
    "  persistent      seed_envelope:     0.0000 GiB\n"
    "  persistent   anderson_history:     0.0000 GiB\n"
    "  persistent      anderson_work:     0.0000 GiB\n"
    "  persistent          collapses:     0.0001 GiB\n"
    "  persistent      eigenpairs_hf:     0.0000 GiB\n"
    "  persistent        flex_arrays:     0.0000 GiB\n"
    "  persistent          hf_tables:     0.0000 GiB\n"
    "  persistent second_order_factors:     0.0000 GiB\n"
    "  phase                green_mu:     0.0000 GiB\n"
    "  phase              density_hf:     0.0000 GiB\n"
    "  phase                  bubble:     0.0003 GiB\n"
    "  phase                dressing:     0.0046 GiB\n"
    "  phase               transport:     0.0002 GiB\n"
    "  phase             convergence:     0.0001 GiB\n"
    "  phase                  mixing:     0.0000 GiB\n"
    "  phase                post_scf:     0.0000 GiB",
)


if __name__ == "__main__":
    unittest.main()
