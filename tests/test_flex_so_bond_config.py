"""Configuration and admission of the spin-orbital bond-resolved FLEX
(issue #270): ``[mode] enable_spin_orbital = true`` with
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

T_OUTPUT_LAYOUT_VALUE = (
    "[mode.param] longitudinal_bond_output_layout must be \"npz\" or \"sidecar\" (\"npz\" "
    "writes the single-file archive, the default; \"sidecar\" writes an index npz plus the "
    "two big channel members as memory-mappable .npy files next to it), got {!r}")
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
# the bond archive (longitudinal_bond_output_full = true)
T_NEED_GREEN = (H + ": longitudinal_bond_output_full = true needs [file.output] green (the bond "
                "archive names the Green function its susceptibility was built from as its "
                "companion green_file; [file.output] keys are case-sensitive); add green to "
                "[file.output].")
T_GREEN_DIR = (H + ": the bond archive '{archive}' and its companion green file '{green}' must "
               "be in the same directory or the green file below it (green_file is stored as a "
               "relative path without '..').")
T_NPZ_TOO_LARGE = (H + ": longitudinal_bond_output_layout = \"npz\" cannot hold the bond archive "
                   "member {member} ({size}, above the small-array bound {bound} = min(1 GiB, "
                   "cap / 8)); set longitudinal_bond_output_layout = \"sidecar\" or remove it "
                   "(the sidecar layout is then chosen).")
T_FORCED = ("spin-orbital bond archive: sidecar layout (the member {member} of {size} is above "
            "the small-array bound {bound}); the susceptibility is written as .npy files next "
            "to the index.")
T_ARCHIVE_KEY_IGNORED = ("[file.output] longitudinal_bond is written only with [mode.param] "
                         "longitudinal_bond_output_full = true; no bond archive is written by "
                         "this run.")
T_NO_ARCHIVE_ITER0 = ("spin-orbital bond archive: IterationMax = 0 runs no map, so no bond "
                      "archive is written.")
T_REMEDY_ARCHIVE = ("longitudinal_bond_output_full = true adds the archive phases "
                    "(archive_bubble, archive_dressing, archive_write); setting it to false "
                    "removes them.")
#: the [mode.param] key of the archive, in mixed case
FULL = {"Longitudinal_Bond_Output_Full": True}

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
            (T_OUTPUT_LAYOUT_VALUE, {"longitudinal_bond_output_layout": "bogus"}),
            (T_PAIRING, {"longitudinal_bond_pairing": "singlet"}),
            (T_PAIRING, {"longitudinal_bond_pairing": "both"}),
            (T_STORE, {"gpu": True, "longitudinal_bond_store": "device"}),
            (T_STORE, {"longitudinal_bond_store": "device"}),
            (T_STORE, {"gpu": False, "longitudinal_bond_store": "device"}),
            (T_GUARD_FREQS, {"longitudinal_bond_guard_freqs": "static"}),
            (T_GPU_REQUIRED, {"gpu_required": True}),
            (T_LAYOUT_VALUE, {"so_bond_layout": "bogus"}),
            (T_LAYOUT_FORCED, {"so_bond_layout": "sz_packed"}),
        )

    def test_each_row(self):
        for text, param in self._rows():
            mixed = {_mixed(k): (_mixed(v) if isinstance(v, str) else v)
                     for k, v in param.items()}
            for p in (param, mixed):
                exp = text
                if text in (T_LAYOUT_VALUE, T_OUTPUT_LAYOUT_VALUE):
                    exp = text.format([v for v in p.values()][0])
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
                  {"longitudinal_bond_output_full": True},
                  {"Longitudinal_Bond_Output_Full": True,
                   "Longitudinal_Bond_Output_Layout": " SideCar "},
                  {"longitudinal_bond_output_layout": "npz"},
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

    def test_sidecar_warning(self):
        """A sidecar layout without longitudinal_bond_output_full writes no
        archive: the WARNING of the non-spin-orbital bond FLEX, once."""
        for p in ({"longitudinal_bond_output_layout": "sidecar"},
                  {"Longitudinal_Bond_Output_Layout": "SiDeCaR"}):
            with self.subTest(param=p):
                with self.capture() as msgs:
                    s, _ = self.build(param=p)
                self.assertEqual(sum(T_SIDECAR in m for lv, m in msgs if lv == "WARNING"), 1,
                                 msgs)
                self.assertEqual(s._so_output_layout, "sidecar")

    def test_no_nonso_pairlift_info(self):
        _, msgs = self.run_logs(F.SOC_ALL_TYPES)
        self.assertOnce(msgs, "INFO", T_PAIRLIFT)
        self.assertOnce(msgs, "INFO", T_SOC_REASON)
        self.assertFalse(any(T_NONSO_PAIRLIFT in m for _, m in msgs), msgs)


class TestArchiveConfig(_Case):
    """The bond archive's configuration (``longitudinal_bond_output_full =
    true``): the layout request, the archive rows of the admission, the
    output-path refusals and the notices."""

    def test_layout_request(self):
        from hwave.solver import so_bond_archive as sa
        rows = (({"longitudinal_bond_output_full": True}, None, "npz"),
                (dict(FULL, LONGITUDINAL_BOND_OUTPUT_LAYOUT=" SIDECAR "), "sidecar", "sidecar"),
                (dict(FULL, longitudinal_bond_output_layout="NpZ"), "npz", "npz"),
                ({"longitudinal_bond_output_full": False}, None, None))
        for fx, member, nbytes in ((X.ZEEMAN_SZ_ALL, "chi_w@0", 409600),
                                   (F.SOC_ALL_TYPES, "chi_w", 1638400)):
            for p, request, layout in rows:
                with self.subTest(fixture=fx.name, param=p):
                    s, gi = self.build(fx, param=dict(p, longitudinal_bond_memory_cap_gb=0))
                    self.assertEqual(s._so_output_layout, request)
                    adm = self.preflight(s, gi)
                    if layout is None:
                        self.assertIsNone(adm.archive)
                        self.assertFalse(adm.plan.keep_chi)
                        self.assertFalse([r for r in adm.plan.rows
                                          if r[0].startswith("archive_")])
                        continue
                    self.assertEqual(adm.archive, sa.ArchiveLayout(layout, False, member, nbytes,
                                                                   float(2 ** 30)))
                    self.assertTrue(adm.plan.keep_chi)
                    self.assertEqual(sorted({r[0] for r in adm.plan.rows
                                             if r[0].startswith("archive_")}),
                                     ["archive_bubble", "archive_dressing", "archive_write"])

    def test_archive_rows_at_the_call_site(self):
        """Every plan of an archive run is built with keep_chi, the archive
        remedy and the archive rows of flex_so_bond_rows(archive=(layout,
        largest inline member))."""
        from hwave.solver import flex_so_bond, so_bond
        from hwave.solver import so_bond_archive as sa
        real_plan, real_rows = so_bond.plan_so_bond, flex_so_bond.flex_so_bond_rows
        for layout in ("npz", "sidecar"):
            with self.subTest(layout=layout):
                plans, rows = [], []

                def plan(**kw):
                    plans.append(kw)
                    return real_plan(**kw)

                def flex_rows(**kw):
                    rows.append(kw)
                    return real_rows(**kw)
                s, gi = self.build(param=dict(FULL, Longitudinal_Bond_Output_Layout=layout,
                                              longitudinal_bond_memory_cap_gb=0))
                with mock.patch.object(so_bond, "plan_so_bond", plan), \
                        mock.patch.object(flex_so_bond, "flex_so_bond_rows", flex_rows):
                    adm = self.preflight(s, gi)
                transfer, extern = s._so_band_tables()[:2]
                compiled = s.ham_info.so_context.compiled
                inline = sa.inline_max_bytes(
                    nmat=16, nvol=16, B=5, nd=2, block_layout="sz_packed", output_layout=layout,
                    n_transfer=len(transfer), n_extern=len(extern),
                    n_rec=len(compiled.records), n_topology=len(compiled.declared_topology))
                self.assertEqual([r["archive"] for r in rows], [(layout, inline)])
                self.assertEqual(len(plans), 2)                  # one trial (nb = 16), final
                for kw in plans:
                    self.assertIs(kw["keep_chi"], True)
                    self.assertEqual(kw["archive_remedy"], T_REMEDY_ARCHIVE)
                    self.assertIn(("archive_write", "archive_write",
                                   min(inline, 16 * 2 ** 20) + 2 ** 20), kw["extra_rows"])
                self.assertEqual(adm.archive.layout, layout)

    def test_need_green(self):
        for info in ({}, {"Green": "green"}):
            with self.subTest(output=info):
                s, _ = self.build(param=FULL)
                with self.assertRaises(ValueError) as cm:
                    s.validate_output_paths(dict(info, path_to_output=self.d), self.d)
                self.assertEqual(str(cm.exception), T_NEED_GREEN)
        # accepted: green given; no archive (output_full false, IterationMax = 0)
        s, _ = self.build(param=FULL)
        s.validate_output_paths({"path_to_output": self.d, "green": "g"}, self.d)
        for p in ({"longitudinal_bond_output_full": False}, dict(FULL, IterationMax=0)):
            with self.subTest(param=p):
                s, _ = self.build(param=p)
                s.validate_output_paths({"path_to_output": self.d}, self.d)

    def test_green_dir(self):
        s, _ = self.build(param=FULL)
        for out in ({"green": "green", "longitudinal_bond": "sub/lb"},
                    {"green": "../g", "longitudinal_bond": "lb"}):
            with self.subTest(output=out):
                with self.assertRaises(ValueError) as cm:
                    s.validate_output_paths(dict(out, path_to_output=self.d), self.d)
                archive = os.path.abspath(os.path.join(self.d, out["longitudinal_bond"] + ".npz"))
                green = os.path.abspath(os.path.join(self.d, out["green"] + ".npz"))
                self.assertEqual(str(cm.exception),
                                 T_GREEN_DIR.format(archive=archive, green=green))
        for out in ({"green": "sub/g", "longitudinal_bond": "lb"},
                    {"green": "g", "longitudinal_bond": "x/../lb"}):
            with self.subTest(output=out):
                s.validate_output_paths(dict(out, path_to_output=self.d), self.d)

    def test_green_dir_shared_rule(self):
        """The output-path check applies the archive's own green_file rule,
        and a green file on another drive (no relative path) is the same
        refusal, not a bare error of os.path.relpath."""
        from hwave.solver import so_bond_archive
        out = {"green": "g", "longitudinal_bond": "lb", "path_to_output": self.d}
        archive = os.path.abspath(os.path.join(self.d, "lb.npz"))
        green = os.path.abspath(os.path.join(self.d, "g.npz"))
        text = T_GREEN_DIR.format(archive=archive, green=green)
        s, _ = self.build(param=FULL)
        s.validate_output_paths(dict(out), self.d)                  # positive control
        calls = []

        def rule(rel, archive_dir):
            calls.append((rel, archive_dir))
            return False
        with mock.patch.object(so_bond_archive, "green_file_ok", rule):
            with self.assertRaises(ValueError) as cm:
                s.validate_output_paths(dict(out), self.d)
        self.assertEqual(calls, [("g.npz", os.path.dirname(archive))])
        self.assertEqual(str(cm.exception), text)

        def other_drive(path, start=os.curdir):
            raise ValueError("path is on mount 'C:', start on mount 'D:'")
        with mock.patch.object(os.path, "relpath", other_drive):
            with self.assertRaises(ValueError) as cm:
                s.validate_output_paths(dict(out), self.d)
        self.assertEqual(str(cm.exception), text)

    def test_no_archive_iter0(self):
        for p in (dict(FULL, IterationMax=0), {"longitudinal_bond_output_full": True,
                                               "iterationmax": 0}):
            with self.subTest(param=p):
                s, msgs = self.run_logs(param=p)
                self.assertOnce(msgs, "INFO", T_NO_ARCHIVE_ITER0)
                self.assertIsNone(s._so_admission.archive)
                self.assertFalse(s._so_admission.plan.keep_chi)
        for p in ({"IterationMax": 0}, FULL):
            with self.subTest(param=p):
                _, msgs = self.run_logs(param=p)
                self.assertNotIn(("INFO", T_NO_ARCHIVE_ITER0), msgs)

    def test_npz_too_large_before_allocation(self):
        """An explicit "npz" with a slot member above the small-array bound
        (cap / 8) is refused by the admission, before any bond array."""
        from hwave.solver import so_bond
        for fx, member, size in ((X.ZEEMAN_SZ_ALL, "chi_w@0", "0.39 MiB"),
                                 (F.SOC_ALL_TYPES, "chi_w", "1.56 MiB")):
            with self.subTest(fixture=fx.name):
                text = T_NPZ_TOO_LARGE.format(member=member, size=size, bound="0.01 MiB")
                param = dict(FULL, Longitudinal_Bond_Output_Layout="NPZ",
                             longitudinal_bond_memory_cap_gb=1e-4)
                s, gi = self.build(fx, param=param)
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
                self.assertEqual(str(cm.exception), text)
                self.assertLessEqual(peak, 2 ** 20)
                s, gi = self.build(fx, param=param)
                with mock.patch.object(so_bond, "open_store", wraps=so_bond.open_store) as spy:
                    with self.assertRaises(ValueError) as cm:
                        s.solve(gi, self.d)
                self.assertEqual(str(cm.exception), text)
                self.assertEqual(spy.call_count, 0)

    def test_forced_sidecar(self):
        """Unset layout and a member above the bound: the sidecar layout and
        one INFO; an explicit choice is never overridden."""
        from hwave.solver import so_bond
        forced = T_FORCED.format(member="chi_w@0", size="0.39 MiB", bound="0.12 MiB")
        for p, layout, info in ((FULL, "sidecar", 1),
                                (dict(FULL, longitudinal_bond_output_layout="Sidecar"),
                                 "sidecar", 0)):
            with self.subTest(param=p):
                with mock.patch.object(so_bond, "small_array_bound",
                                       return_value=float(2 ** 17)):
                    s, msgs = self.run_logs(param=p)
                self.assertEqual(msgs.count(("INFO", forced)), info)
                self.assertEqual(s._so_admission.archive.layout, layout)
                self.assertIs(s._so_admission.archive.forced, bool(info))
        with mock.patch.object(so_bond, "small_array_bound", return_value=float(2 ** 17)):
            with self.assertRaises(ValueError) as cm:
                self.run_logs(param=dict(FULL, longitudinal_bond_output_layout="npz"))
        self.assertEqual(str(cm.exception), T_NPZ_TOO_LARGE.format(
            member="chi_w@0", size="0.39 MiB", bound="0.12 MiB"))
        _, msgs = self.run_logs(param=FULL)
        self.assertFalse([m for _, m in msgs if m.startswith(forced[:40])])

    def test_remedy_archive(self):
        """An over-cap refusal in an archive phase names the key that adds
        the archive phases, after the batch note; without the archive it does
        not."""
        for p, note in ((FULL, T_REMEDY_AUTO),
                        (dict(FULL, longitudinal_bond_freq_batch=4), T_REMEDY_SET)):
            with self.subTest(param=p):
                s, gi = self.build(param=dict(p, longitudinal_bond_memory_cap_gb=1e-6))
                with self.assertRaises(ValueError) as cm:
                    self.preflight(s, gi)
                msg = str(cm.exception)
                self.assertIn(" in phase archive_", msg.split("\n")[0])
                self.assertEqual(msg.count(T_REMEDY_ARCHIVE), 1)
                self.assertIn("(" + CAP_KEY + "). " + note + " " + T_REMEDY_ARCHIVE
                              + " Memory grows", msg)
        s, gi = self.build(param={"longitudinal_bond_memory_cap_gb": 1e-6})
        with self.assertRaises(ValueError) as cm:
            self.preflight(s, gi)
        self.assertNotIn(T_REMEDY_ARCHIVE, str(cm.exception))

    def test_archive_key_ignored(self):
        """[file.output] longitudinal_bond without output_full: one INFO from
        save_results, no archive; not on the non-spin-orbital gate."""
        cases = ((X.ZEEMAN_SZ_ALL, {"Longitudinal_Bond_Output_Full": False}, True, 1),
                 (X.ZEEMAN_SZ_ALL, FULL, True, 0),
                 (F.SU2_4X4, {}, False, 0))
        for fx, p, so, count in cases:
            with self.subTest(fixture=fx.name, param=p):
                d = self._dir()
                inp = X.write_flex_case(d, fx, so=so, param=p,
                                        output={"green": "green", "longitudinal_bond": "lb"})
                info = inp["file"]["output"]
                out = info["path_to_output"]
                os.makedirs(out)
                with self.capture() as msgs:
                    s, gi = X.build_solver(inp)
                    s.validate_output_paths(info, out)
                    s.solve(gi, out)
                    s.save_results(info, gi)
                self.assertEqual(msgs.count(("INFO", T_ARCHIVE_KEY_IGNORED)), count)
                self.assertEqual(os.path.exists(os.path.join(out, "lb.npz")),
                                 bool(so and count == 0))


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
