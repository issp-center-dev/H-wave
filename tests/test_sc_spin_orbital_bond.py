"""hwave_sc in spin-orbital mode with bond-resolved interaction channels
(issue #270): the configuration of the static spin-orbital bond path and
the refusals that come before any file is read."""

import functools
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from requests.structures import CaseInsensitiveDict

import hwave.sc as sc
from hwave.solver import backend, fingerprint, so_bond, so_pairing, spin_orbital
from hwave.sc import (_SO_BOND_CHI0Q_UNUSED, _SO_BOND_DIAGNOSTICS_REFUSAL,
                      _SO_BOND_FLEX_KEYS_UNUSED, _SO_BOND_FLEX_REFUSAL,
                      _SO_BOND_GPU_IGNORED, _SO_BOND_GPU_REQUIRED_REFUSAL,
                      _SO_BOND_GREEN_REFUSAL, _SO_BOND_LAYOUT_IGNORED,
                      _SO_BOND_PRECONDITION_REFUSAL)
from tests import so_bond_entry_fixtures as F
from tests.so_bond_oracles import CELL, closed_row, rel


_STATIC_BOND_TEXT = (
    "[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
    "bond_channels = true is not supported; the spin-orbital pairing vertex is "
    "on-site only. Set bond_channels = false.")


class TestDynamicRefusal(unittest.TestCase):
    """bond_channels = true with frequency = "dynamic" has its own refusal,
    raised by the configuration guard before any file is read (a mode and an
    eliashberg section suffice: a later check would raise KeyError)."""

    def _refusal(self, eli):
        with self.assertRaises(ValueError) as cm:
            sc.calc_eliashberg({"mode": {"Enable_Spin_Orbital": True},
                                "eliashberg": eli})
        return str(cm.exception)

    def test_dynamic_bond_refused(self):
        self.assertEqual(
            so_pairing._SO_DYNAMIC_BOND_REFUSAL,
            "[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
            "bond_channels = true with frequency = \"dynamic\" is not available yet "
            "(issue #270); use frequency = \"static\", or bond_channels = false for "
            "on-site interactions.")
        self.assertEqual(self._refusal({"Frequency": "dynamic", "Bond_Channels": True}),
                         so_pairing._SO_DYNAMIC_BOND_REFUSAL)
        # the static row keeps its own text
        self.assertEqual(self._refusal({"Bond_Channels": True}), _STATIC_BOND_TEXT)


_W = so_pairing._SO_WHAT
_CAP_KEY = "[eliashberg] bond_memory_cap_gb"


class _Warnings(logging.Handler):
    def __init__(self):
        super().__init__(logging.WARNING)
        self.seen = []

    def emit(self, record):
        self.seen.append(record.getMessage())


class TestConfigParsing(unittest.TestCase):
    """The [eliashberg] options of the spin-orbital bond path: every option
    the path does not use is refused with its own text, in a fixed order
    (the first offender wins), keys case-insensitively; options that are
    merely unused get one WARNING each."""

    def _parse(self, eli):
        return sc._read_so_bond_config(CaseInsensitiveDict(eli))

    def _refusal(self, eli):
        with self.assertRaises(ValueError) as cm:
            self._parse(eli)
        return str(cm.exception)

    def test_texts(self):
        self.assertEqual(_SO_BOND_LAYOUT_IGNORED,
                         "[eliashberg] so_bond_layout set but bond_channels=false; it only "
                         "applies to the spin-orbital bond path and is ignored here.")
        self.assertEqual(
            _SO_BOND_GREEN_REFUSAL,
            _W + ": [eliashberg] bond_green is not available with bond_channels in "
            "spin-orbital mode: no run writes a spin-orbital bond-resolved Green function "
            "yet (it arrives with the spin-orbital bond FLEX). Remove bond_green; the Green "
            "function is built from the Transfer (and Extern) input.")
        self.assertEqual(
            _SO_BOND_DIAGNOSTICS_REFUSAL,
            _W + ": [eliashberg] bond_diagnostics is not available in spin-orbital mode "
            "(its character analysis is the singlet/triplet one of the non-spin-orbital "
            "bond path); the gap.npz S_z and spin decompositions describe the leading "
            "state. Remove bond_diagnostics.")
        self.assertEqual(
            _SO_BOND_PRECONDITION_REFUSAL.format(key="bond_precondition_atol"),
            _W + ": [eliashberg] bond_precondition_atol tunes the Hermiticity check of "
            "the non-spin-orbital bond kernel; the spin-orbital kernel is not Hermitian "
            "and has no such check. Remove bond_precondition_atol.")
        self.assertEqual(
            _SO_BOND_FLEX_REFUSAL,
            _W + ": [eliashberg] chi0q_mode = \"flex\" with bond_channels = true is not "
            "available yet: the spin-orbital bond path builds its bond bubble from the "
            "band's Green function and reads no FLEX output. Remove chi0q_mode (or leave "
            "it unset); it is not used on this path.")
        self.assertEqual(
            _SO_BOND_GPU_REQUIRED_REFUSAL,
            _W + ": [eliashberg] gpu_required = true is not available with bond_channels "
            "= true: the spin-orbital bond-resolved kernel runs on the CPU only. Remove "
            "gpu_required (gpu = true alone is ignored with a warning).")
        self.assertEqual(
            _SO_BOND_CHI0Q_UNUSED,
            "[eliashberg] bond_channels=true in spin-orbital mode builds its own "
            "bond-resolved bubble directly from the Green function, so chi0q_mode = '%s' "
            "and chi0q_tensor = '%s' are not used: no chi0q file is read or computed on "
            "this path.")
        self.assertEqual(
            _SO_BOND_GPU_IGNORED,
            "[eliashberg] gpu=true is ignored for bond_channels=true in spin-orbital mode: "
            "the spin-orbital bond-resolved kernel is CPU-only (GPU support is a deferred "
            "follow-up).")
        self.assertEqual(
            _SO_BOND_FLEX_KEYS_UNUSED,
            "[eliashberg] %s: not used for bond_channels=true in spin-orbital mode (the "
            "spin-orbital bond path reads no FLEX output).")

    def test_each_refusal(self):
        on = {"bond_channels": True}
        rows = [
            ("bond_green", {"bond_green": "g.npz"}, _SO_BOND_GREEN_REFUSAL),
            ("bond_diagnostics true", {"bond_diagnostics": True},
             _SO_BOND_DIAGNOSTICS_REFUSAL),
            ("bond_diagnostics false", {"bond_diagnostics": False},
             _SO_BOND_DIAGNOSTICS_REFUSAL),
            ("bond_precondition", {"bond_precondition_dense_limit": False},
             _SO_BOND_PRECONDITION_REFUSAL.format(key="bond_precondition_dense_limit")),
            ("flex", {"chi0q_mode": "flex"}, _SO_BOND_FLEX_REFUSAL),
            ("gpu_required", {"gpu_required": "true", "gpu": True},
             _SO_BOND_GPU_REQUIRED_REFUSAL),
            ("layout", {"so_bond_layout": "sz_packed"},
             '[eliashberg] so_bond_layout = "sz_packed" cannot be forced: the S_z-packed '
             'layout is chosen automatically when the input conserves S_z (no spin-flip '
             'Transfer entry, no S_z-changing interaction). Use "auto" (the default) or '
             '"dense".'),
            ("layout value", {"so_bond_layout": "packed"},
             "[eliashberg] so_bond_layout must be \"auto\" or \"dense\", got 'packed'"),
            ("max_shells", {"bond_max_shells": -1},
             "[eliashberg] bond_max_shells must be >= 0 (shell 0 = the on-site Delta r = 0 "
             "point), got -1"),
        ]
        for raw in (-1.0, -1, float("inf"), float("nan"), True):
            rows.append(("cap %r" % (raw,), {"bond_memory_cap_gb": raw},
                         so_bond.CAP_RANGE_TEXT.format(key=_CAP_KEY, raw=raw)))
        for label, extra, text in rows:
            with self.subTest(row=label):
                self.assertEqual(self._refusal(dict(on, **extra)), text)

    def test_refusal_order(self):
        on = {"bond_channels": True}
        self.assertEqual(self._refusal(dict(on, bond_green="g.npz", bond_diagnostics=True)),
                         _SO_BOND_GREEN_REFUSAL)
        self.assertEqual(
            self._refusal(dict(on, Bond_Precondition_RTOL=1e-3, bond_precondition_atol=0.0)),
            _SO_BOND_PRECONDITION_REFUSAL.format(key="bond_precondition_atol"))
        # every row before the layout one wins over a bad layout
        self.assertEqual(self._refusal(dict(on, chi0q_mode="flex", gpu_required=True,
                                            so_bond_layout="x", bond_max_shells=-1)),
                         _SO_BOND_FLEX_REFUSAL)
        self.assertEqual(self._refusal(dict(on, gpu_required=True, so_bond_layout="x")),
                         _SO_BOND_GPU_REQUIRED_REFUSAL)

    def test_case_insensitive_keys(self):
        self.assertEqual(self._refusal({"Bond_Channels": True, "Bond_Green": "g.npz"}),
                         _SO_BOND_GREEN_REFUSAL)
        self.assertEqual(self._refusal({"Bond_Channels": "on", "Bond_Diagnostics": False}),
                         _SO_BOND_DIAGNOSTICS_REFUSAL)
        self.assertEqual(self._refusal({"BOND_CHANNELS": True, "Chi0q_Mode": "flex"}),
                         _SO_BOND_FLEX_REFUSAL)
        self.assertEqual(self._refusal({"bond_channels": True, "GPU_Required": True}),
                         _SO_BOND_GPU_REQUIRED_REFUSAL)
        cfg = self._parse({"Bond_Channels": True, "SO_Bond_Layout": "Dense"})
        self.assertEqual(cfg.layout_request, "dense")
        cfg = self._parse({"Bond_Channels": True, "Bond_Memory_Cap_GB": 0})
        self.assertEqual(cfg.cap_raw, 0.0)
        cfg = self._parse({"Bond_Channels": True, "Bond_Max_Shells": 2,
                           "Bond_Coeff_Tail": 0.5, "G2_Tail": False,
                           "Bond_Memory_Cap_GB": 3})
        self.assertEqual(cfg, sc._SOBondConfig(True, "auto", 3.0, 2, 0.5, "config", False))

    def test_defaults(self):
        self.assertEqual(self._parse({"bond_channels": True}),
                         sc._SOBondConfig(True, "auto", None, None, 1.0, "default", True))
        self.assertEqual(self._parse({}),
                         sc._SOBondConfig(False, None, None, None, None, None, None))

    def test_bond_off_warnings(self):
        h = _Warnings()
        log = logging.getLogger("hwave_sc")
        log.addHandler(h)
        self.addCleanup(log.removeHandler, h)
        # bond-off: refused keys are never parsed, only named as ignored
        cfg = self._parse({"Bond_Channels": False, "Bond_Green": "g.npz",
                           "bond_diagnostics": "ture", "SO_Bond_Layout": "nonsense"})
        self.assertEqual(cfg, sc._SOBondConfig(False, None, None, None, None, None, None))
        self.assertEqual(h.seen, [
            "[eliashberg] bond_diagnostics, bond_green set but bond_channels=false; the "
            "bond-channel options only apply to bond_channels=true and are ignored here.",
            _SO_BOND_LAYOUT_IGNORED])
        h.seen.clear()
        self._parse({"bond_channels": False})
        self.assertEqual(h.seen, [])

    def test_non_so_cap_unchanged(self):
        with self.assertRaises(ValueError) as cm:
            sc._read_bond_config({"bond_channels": True, "bond_memory_cap_gb": 0})
        self.assertEqual(str(cm.exception),
                         "[eliashberg] bond_memory_cap_gb must be a positive finite "
                         "number, got 0")
        self.assertEqual(sc._read_bond_config({"bond_channels": True})[3], 8.0)

    def test_entry_warnings_once(self):
        log = logging.getLogger("test_sc_spin_orbital_bond.entry")
        log.propagate = False
        h = _Warnings()
        log.addHandler(h)
        self.addCleanup(log.removeHandler, h)
        sc._so_bond_entry_warnings(CaseInsensitiveDict(
            {"Chi0q_Mode": "calc", "GPU": "true", "Flex_Chi_S": "a.npz",
             "flex_chi_c": "b.npz"}), log)
        self.assertEqual(h.seen, [
            _SO_BOND_CHI0Q_UNUSED % ("calc", "unset"),
            _SO_BOND_GPU_IGNORED,
            _SO_BOND_FLEX_KEYS_UNUSED % ("flex_chi_c, Flex_Chi_S",)])
        h.seen.clear()
        sc._so_bond_entry_warnings(CaseInsensitiveDict(
            {"chi0q_tensor": "general", "flex_chi_c": "b.npz"}), log)
        self.assertEqual(h.seen, [_SO_BOND_CHI0Q_UNUSED % ("unset", "general"),
                                  _SO_BOND_FLEX_KEYS_UNUSED % ("flex_chi_c",)])
        h.seen.clear()
        sc._so_bond_entry_warnings(CaseInsensitiveDict({"gpu": False}), log)
        sc._so_bond_entry_warnings(CaseInsensitiveDict(), log)
        self.assertEqual(h.seen, [])

    def test_onsite_run_parses_first(self):
        """A public spin-orbital on-site static run: the spin-orbital bond
        options are parsed before the backend is selected, and a stray
        so_bond_layout is named as ignored once."""
        from tests import test_sc_spin_orbital as S
        d = tempfile.mkdtemp(prefix="sc_so_bond_cfg_")
        self.addCleanup(shutil.rmtree, d, True)
        files = S._write_one_orbital(d)
        inp = S._input(d, files, True, {"chi0q_mode": "calc", "solver_mode": "eigenvalue",
                                         "num_eigenvalues": 2, "SO_Bond_Layout": "dense"},
                       nmat=8, CellShape=[2, 2, 1])
        calls = []
        real_reader, real_backend = sc._read_so_bond_config, backend.get_backend

        def reader(*a, **kw):
            calls.append("read_so_bond_config")
            return real_reader(*a, **kw)

        def get_backend(*a, **kw):
            calls.append("get_backend")
            return real_backend(*a, **kw)

        with S._count_warnings() as warned, \
                mock.patch.object(sc, "_read_so_bond_config", side_effect=reader), \
                mock.patch.object(backend, "get_backend", side_effect=get_backend):
            sc.calc_eliashberg(inp)
        self.assertEqual(calls, ["read_so_bond_config", "get_backend"])
        self.assertEqual(warned.count(_SO_BOND_LAYOUT_IGNORED), 1)
        self.assertIn("eigenvalue.dat", os.listdir(os.path.join(d, "out")))


# ------------------------------------------------------------ run record
_RUN_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond")
#: a 2x2 one-orbital model: nearest-neighbour hopping along x and y (both
#: spins, interleaved spin-orbital indices), on-site U and a
#: nearest-neighbour CoulombInter along x
_RUN_TRANSFER = {}
for _R in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
    for _s in (0, 1):
        _RUN_TRANSFER[(_R, (_s, _s))] = -1.0 + 0j
_RUN_TABLES = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0},
               "CoulombInter": closed_row((1, 0, 0), 0, 0, 0.5)}
_RUN_INPUT = {"mode": {"param": {"CellShape": [2, 2, 1]}}}
_RUN_NMAT = 4


def _context():
    return spin_orbital.build_so_entry_context(None, _RUN_TABLES, 1, _RUN_LOG)


def _run(context=None, bond=None):
    """The run record of the 2x2 model (beta 10, Nmat 4, filling 0.5)."""
    context = _context() if context is None else context
    return sc._so_run_inputs(_RUN_INPUT, _RUN_TRANSFER, None, 1, 10.0, _RUN_NMAT, 0.5,
                             context=context, bond=bond)


class TestRunRecord(unittest.TestCase):
    """The run record of a spin-orbital hwave_sc run: the entry context is
    required, the bond topology is resolved once and carried in the record,
    the gap writer records the bond members, and the archive check takes the
    capability it needs with a fixed precedence for the bond inputs."""

    def test_context_required(self):
        with self.assertRaises(TypeError):
            sc._so_run_inputs(_RUN_INPUT, _RUN_TRANSFER, None, 1, 10.0, _RUN_NMAT, 0.5)
        with self.assertRaises(TypeError):
            sc._SORun(_run().band, _run().state)
        for bad, name in ((None, "NoneType"), ("context", "str")):
            with self.assertRaises(ValueError) as cm:
                sc._so_run_inputs(_RUN_INPUT, _RUN_TRANSFER, None, 1, 10.0, _RUN_NMAT,
                                  0.5, context=bad)
            self.assertEqual(
                str(cm.exception),
                "_so_run_inputs: context must be the run's SOEntryContext "
                "(spin_orbital.build_so_entry_context), got {}".format(name))
        run = _run()
        self.assertIsNone(run.bond)
        self.assertIsInstance(run.context, spin_orbital.SOEntryContext)

    def test_bond_run_built_once(self):
        ctx = _context()
        with mock.patch.object(so_bond, "so_bond_topology",
                               wraps=so_bond.so_bond_topology) as spy:
            b = sc._so_bond_run(ctx, CELL, None)
        self.assertEqual(spy.call_count, 1)
        self.assertIs(spy.call_args.args[0], ctx.compiled)
        self.assertEqual(spy.call_args.kwargs, {"max_shells": None})
        self.assertIs(_run(ctx, bond=b).bond, b)
        self.assertEqual(b.inputs.max_shells, -1)
        self.assertEqual(sc._so_bond_run(ctx, CELL, 2).inputs.max_shells, 2)
        self.assertIs(b.inputs.slot_convention, so_bond.SLOT_CONVENTION)
        self.assertTrue(np.array_equal(b.inputs.delta_r, b.topo.delta_r))
        self.assertTrue(np.array_equal(b.inputs.reverse, b.topo.reverse))
        self.assertEqual([tuple(r) for r in b.inputs.delta_r],
                         [(0, 0, 0), (-1, 0, 0), (1, 0, 0)])
        self.assertEqual(b.fp_resolved_topology, fingerprint.resolved_topology_digest(
            b.topo.delta_r, b.topo.reverse, so_bond.SLOT_CONVENTION))
        # the cutoff refusal of the topology propagates
        with self.assertRaisesRegex(ValueError, r"max_shells=0 would drop"):
            sc._so_bond_run(ctx, CELL, 0)

    def test_writer_bond_members(self):
        ctx = _context()
        b = sc._so_bond_run(ctx, CELL, None)
        run = _run(ctx, bond=b)
        run = run._replace(state=_finite_state(run))
        fp_green = "0" * 64
        m = sc._so_writer_members("sc_so_bond_static", "test", run, fp_green, run.state)
        self.assertEqual(set(m), fingerprint.required_members("gap", "sc_so_bond_static"))
        self.assertEqual(str(m["fp_resolved_topology"]), b.fp_resolved_topology)
        self.assertEqual(int(m["max_shells"]), -1)
        self.assertTrue(np.array_equal(m["delta_r"], b.inputs.delta_r))
        # an explicit bond overrides the run's
        b2 = sc._so_bond_run(ctx, CELL, 1)
        m2 = sc._so_writer_members("sc_so_bond_static", "test", run, fp_green, run.state,
                                   bond=b2.inputs)
        self.assertEqual(int(m2["max_shells"]), 1)
        # the on-site stage takes no bond inputs, and records the sentinel
        with self.assertRaises(ValueError) as cm:
            sc._so_writer_members("sc_so_static", "test", run, fp_green, run.state,
                                  bond=b.inputs)
        self.assertEqual(str(cm.exception), "_so_writer_members: stage sc_so_static is "
                                            "so_onsite and takes no bond inputs")
        m0 = sc._so_writer_members("sc_so_static", "test", run, fp_green, run.state)
        self.assertEqual(set(m0), fingerprint.required_members("gap", "sc_so_static"))
        self.assertEqual(str(m0["fp_resolved_topology"]),
                         fingerprint.ONSITE_RESOLVED_SENTINEL)

    def test_check_archive_precedence(self):
        ctx = _context()
        b = sc._so_bond_run(ctx, CELL, None)
        cases = (
            (_run(ctx, bond=b), dict(need="so_onsite", run_bond=b.inputs),
             "_so_check_archive: run_bond applies to a so_bond consumer only"),
            (_run(ctx), dict(need="so_bond"),
             "_so_check_archive: a so_bond consumer needs the run's bond topology "
             "(so_run.bond)"),
            (_run(ctx, bond=b), dict(need="so_bond",
                                     run_bond=sc._so_bond_run(ctx, CELL, 1).inputs),
             "_so_check_archive: run_bond differs from the run's bond topology "
             "(so_run.bond.inputs)"),
        )
        # the same cutoff, but other channels, another reversal map or another
        # slot convention: each differs from the run's bond topology
        delta_r = b.inputs.delta_r.copy()
        delta_r[-1] = (2, 0, 0)
        reverse = b.inputs.reverse[::-1].copy()
        self.assertFalse(np.array_equal(reverse, b.inputs.reverse))
        for name, other in (
                ("delta_r", fingerprint.BondInputs(delta_r, b.inputs.reverse.copy(),
                                                   so_bond.SLOT_CONVENTION, -1)),
                ("reverse", fingerprint.BondInputs(b.inputs.delta_r.copy(), reverse,
                                                   so_bond.SLOT_CONVENTION, -1)),
                ("slot_convention", fingerprint.BondInputs(
                    b.inputs.delta_r.copy(), b.inputs.reverse.copy(),
                    so_bond.SLOT_CONVENTION.to_string().replace("sr-", "sr+"), -1))):
            cases += ((_run(ctx, bond=b), dict(need="so_bond", run_bond=other),
                       "_so_check_archive: run_bond differs from the run's bond topology "
                       "(so_run.bond.inputs)", name),)
        for run, kw, text, *field in cases:
            with self.subTest(text=text, field=field):
                with mock.patch.object(fingerprint, "check_so_archive") as check:
                    with self.assertRaises(ValueError) as cm:
                        sc._so_check_archive({}, "x.npz", run, "exact", "green",
                                             sc._SO_FLEX_PRODUCER, **kw)
                self.assertEqual(str(cm.exception), text)
                check.assert_not_called()
        # the bond inputs reaching the archive check: the run's by default,
        # an equal explicit one (a copy) is accepted
        same = fingerprint.BondInputs(b.inputs.delta_r.copy(), b.inputs.reverse.copy(),
                                      so_bond.SLOT_CONVENTION, -1)
        for kw in (dict(need="so_bond"), dict(need="so_bond", run_bond=same)):
            with mock.patch.object(fingerprint, "check_so_archive") as check:
                sc._so_check_archive({}, "x.npz", _run(ctx, bond=b), "exact", "green",
                                     sc._SO_FLEX_PRODUCER, **kw)
            self.assertEqual(check.call_args.kwargs["need"], "so_bond")
            self.assertIs(check.call_args.kwargs["run_bond"],
                          kw.get("run_bond", b.inputs))
        with mock.patch.object(fingerprint, "check_so_archive") as check:
            sc._so_check_archive({}, "x.npz", _run(ctx, bond=b), "exact", "green",
                                 sc._SO_FLEX_PRODUCER)
        self.assertEqual(check.call_args.kwargs["need"], "so_onsite")
        self.assertIsNone(check.call_args.kwargs["run_bond"])


def _finite_state(run):
    """The state of ``run`` at a closed filling (mu 0, density = target), as
    a producer writes it."""
    target = run.state.target
    return fingerprint.so_state(beta=run.state.beta, nmat=run.state.nmat, mu=0.0,
                                density=target, calc_mu=True, target=target)


def _onsite_green_archive(test, run):
    """A version-2 green.npz as the spin-orbital FLEX run writes it, built
    from the band and interaction of ``run`` (so the band and interaction
    checks pass and the capability rule is what decides)."""
    nd, nvol = 2, 4
    green = np.zeros((1, _RUN_NMAT, nvol, nd, nd), dtype=np.complex128)
    state = _finite_state(run)
    members = fingerprint.so_provenance_members(
        producer_stage="flex_so_onsite", producer_label="test FLEX producer",
        transfer=run.band.transfer, extern=run.band.extern, meta=run.band.meta,
        fp_green=fingerprint.green_digest(green, fingerprint.GREEN_LAYOUT),
        state=state, interaction=run.context.interaction)
    d = tempfile.mkdtemp(prefix="sc_so_bond_cap_")
    test.addCleanup(shutil.rmtree, d, True)
    path = os.path.join(d, "green.npz")
    np.savez(path, **members, green=green,
             green_layout=np.str_(fingerprint.GREEN_LAYOUT),
             spin_orbital=np.bool_(True), index_convention=np.str_("spin_block"))
    return path


class TestBondCapability(unittest.TestCase):
    """An on-site spin-orbital FLEX Green function is not a bond artifact:
    a bond consumer refuses it on capability, after its band and
    interaction checks pass and before any payload is read; on-site
    consumers accept the same archive."""

    def test_onsite_green_refused_for_bond_consumer(self):
        from tests.test_so_entry_rules import _PayloadReads
        ctx = _context()
        run = _run(ctx, bond=sc._so_bond_run(ctx, CELL, None))
        path = _onsite_green_archive(self, run)
        with _PayloadReads() as spy, np.load(path) as z:
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                sc._so_check_archive(z, path, run, "exact", "green", sc._SO_FLEX_PRODUCER,
                                     need="so_bond")
        self.assertEqual(str(cm.exception),
                         "{} was produced by test FLEX producer (so_onsite), but this run "
                         "needs a so_bond artifact".format(path))
        self.assertEqual(spy.reads, [])

    def test_onsite_consumers_unchanged(self):
        from tests.test_so_entry_rules import _PayloadReads
        ctx = _context()
        for run in (_run(ctx), _run(ctx, bond=sc._so_bond_run(ctx, CELL, None))):
            path = _onsite_green_archive(self, run)
            with _PayloadReads() as spy, np.load(path) as z:
                info = sc._so_check_archive(z, path, run, "exact", "green",
                                            sc._SO_FLEX_PRODUCER)
            self.assertEqual(info.producer_label, "test FLEX producer")
            self.assertEqual(spy.reads, [])


# ------------------------------------------------------------ Green function
_GREEN_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond.green")
_NMAT, _BETA = 32, 2.0


def _read_only(a):
    if isinstance(a, np.ndarray):
        a.flags.writeable = False
    return a


@functools.lru_cache(maxsize=None)
def _green_case(fixture, layout="auto"):
    """The Green function of the spin-orbital bond path of ``fixture``,
    built through the run's own readers from the files of
    :func:`tests.so_bond_entry_fixtures.write_case` (``so_bond_layout =
    layout``): ``(so_run, epsilon_k, config, _SOBondGreen)``. Cached by
    content (the fixture is hashable); every array is read-only."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_green_")
    try:
        inp = F.write_case(d, fixture, eliashberg={"SO_Bond_Layout": layout})
        param = inp["mode"]["param"]
        geom, hr, tables = sc._read_interaction_files(inp)
        norb_phys = spin_orbital.physical_norb(geom["norb"], True)
        ctx = spin_orbital.build_so_entry_context(param, tables, norb_phys, _GREEN_LOG)
        extern, _ = sc._read_so_extern(inp, norb_phys)
        cfg = sc._read_so_bond_config(CaseInsensitiveDict(inp["eliashberg"]))
        beta = sc._coerce_run_beta(param["T"])
        run = sc._so_run_inputs(inp, hr, extern, norb_phys, beta, param["Nmat"],
                                param["filling"], context=ctx,
                                bond=sc._so_bond_run(ctx, geom["rvec"], cfg.max_shells))
        k = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in param["CellShape"]]
        eps = sc._so_epsilon_k(inp, *k, hr, extern, norb_phys)
        out = sc._so_bond_green(eps, run, norb_phys=norb_phys, filling=param["filling"],
                                beta=beta, nmat=param["Nmat"], coeff_tail=cfg.coeff_tail)
    finally:
        shutil.rmtree(d, True)
    for a in (eps, out.eigenvalues, out.carrier.full_sc, out.carrier.deflated_kw,
              out.carrier.green0_tail):
        _read_only(a)
    return run, eps, cfg, out


def _off_sector(G, norb_phys=1):
    """The S_z-off-sector blocks (up-down, down-up) of a spin-block Green
    function ``(nd, nd, ...)``."""
    n = norb_phys
    return np.concatenate([G[:n, n:].ravel(), G[n:, :n].ravel()])


class TestGreen(unittest.TestCase):
    """The Green function of the static spin-orbital bond path on the
    canonical band (gate A2-P: ``test_green_matches_independent_construction``,
    the Green function against an independent construction from the
    fixture's literal tables). A band without spin-flip entries gives an
    exactly block-diagonal G whatever the layout; mu follows the filling per
    spin-orbital, as in the non-spin-orbital run of the same band."""

    def test_green_matches_independent_construction(self):
        """Independent Green function: mu and G(k, i omega_n) of the bond
        path equal a construction from the fixture's literal hopping,
        spin-flip and Zeeman lists."""
        for fx, conserving in ((F.SOC_ALL_TYPES, False), (F.ZEEMAN_EXTERN, True)):
            with self.subTest(fixture=fx.name):
                _, _, cfg, out = _green_case(fx)
                G = out.carrier.full_sc
                G_ind, mu_ind = F.independent_green(fx, nmat=_NMAT, beta=_BETA)
                self.assertEqual(G.shape, (2, 2, 4, 4, 1, _NMAT))
                self.assertEqual(G.shape, G_ind.shape)
                self.assertLessEqual(abs(out.mu - mu_ind), 1e-12)
                self.assertLessEqual(np.abs(G - G_ind).max(), 1e-12 * np.abs(G_ind).max())
                self.assertIs(out.conserving, conserving)
                # the bond path's tail default, never the on-site 0.0
                self.assertEqual(cfg.coeff_tail, 1.0)
                self.assertEqual(out.carrier.coeff_tail_applied, 1.0)
                self.assertEqual(out.carrier.source, "internal")
                self.assertEqual(out.fp_green, sc._so_green_fp(G))
                off = _off_sector(G)
                if conserving:
                    self.assertTrue(np.all(off == 0))
                else:
                    self.assertGreaterEqual(np.linalg.norm(off),
                                            1e-2 * np.linalg.norm(G))

    def test_zeeman_sign_analytic(self):
        """Zeeman sign oracle: +h on spin up, independent of the fixture tables"""
        for fx in (F.ZEEMAN_EXTERN, F.ZEEMAN_TRANSFER):
            with self.subTest(fixture=fx.name):
                out = _green_case(fx)[3]
                kw = dict(fixture=fx, beta=_BETA, nmat=_NMAT)
                self.assertLessEqual(
                    F.zeeman_sign_miss(out.carrier.full_sc, out.mu, h=F.ZEEMAN_H, **kw),
                    1e-12)
                # the flipped sign is visible
                self.assertGreaterEqual(
                    F.zeeman_sign_miss(out.carrier.full_sc, out.mu, h=-F.ZEEMAN_H, **kw),
                    1e-3)

    def test_zeeman_routes_agree(self):
        g_ext = _green_case(F.ZEEMAN_EXTERN)[3]
        g_tr = _green_case(F.ZEEMAN_TRANSFER)[3]
        self.assertLessEqual(rel(g_tr.carrier.full_sc, g_ext.carrier.full_sc), 1e-13)
        self.assertLessEqual(abs(g_tr.mu - g_ext.mu), 1e-13)
        self.assertTrue(g_ext.conserving and g_tr.conserving)

    def test_fp_green_layout_independent(self):
        fx = F.ZEEMAN_EXTERN
        run_d, _, cfg_d, dense = _green_case(fx, "dense")
        run_a, _, cfg_a, auto = _green_case(fx, "auto")
        self.assertEqual((cfg_d.layout_request, cfg_a.layout_request), ("dense", "auto"))
        # the two runs would solve in different layouts
        self.assertEqual(so_bond.select_layout(run_d.context.compiled, run_d.band.transfer,
                                               override="dense")[0], "dense")
        self.assertEqual(so_bond.select_layout(run_a.context.compiled, run_a.band.transfer,
                                               override="auto")[0], "sz_packed")
        # the band's S_z check is evaluated even when the layout is forced dense
        self.assertTrue(dense.conserving and auto.conserving)
        # ... and decides the per-sector construction. Only this spy can see
        # it: LAPACK's Hermitian reduction keeps the exact zeros of a
        # block-diagonal band, so a full diagonalization gives the same
        # exactly block-diagonal G and no value test separates the two.
        eps = _green_case(fx, "dense")[1]
        with mock.patch.object(so_bond, "canonical_band_eigen",
                               wraps=so_bond.canonical_band_eigen) as spy:
            again = sc._so_bond_green(eps, run_d, norb_phys=1, filling=fx.filling,
                                      beta=_BETA, nmat=_NMAT, coeff_tail=cfg_d.coeff_tail)
        self.assertEqual(spy.call_count, 1)
        self.assertIs(spy.call_args.kwargs["conserving"], True)
        # the fingerprints are exact (the archive check compares them); the
        # recomputed values agree to round-off
        self.assertEqual(again.fp_green, dense.fp_green)
        self.assertEqual(dense.fp_green, auto.fp_green)
        self.assertLessEqual(rel(dense.carrier.full_sc, auto.carrier.full_sc), 1e-14)
        self.assertLessEqual(abs(dense.mu - auto.mu), 1e-14)

    def test_carrier_read_only(self):
        """The returned carrier cannot be written in place (its fp_green is
        taken at build time), and coeff_tail must be the resolved value."""
        fx = F.ZEEMAN_EXTERN
        run, eps, cfg, _ = _green_case(fx)
        out = sc._so_bond_green(eps, run, norb_phys=1, filling=fx.filling, beta=_BETA,
                                nmat=_NMAT, coeff_tail=cfg.coeff_tail)
        c = out.carrier
        for name in ("full_sc", "deflated_kw", "green0_tail"):
            with self.subTest(array=name):
                a = getattr(c, name)
                self.assertFalse(a.flags.writeable)
                with self.assertRaises(ValueError):
                    a[(0,) * a.ndim] = 1.0
        self.assertEqual(out.fp_green, sc._so_green_fp(c.full_sc))
        with self.assertRaises(ValueError) as cm:
            sc._so_bond_green(eps, run, norb_phys=1, filling=fx.filling, beta=_BETA,
                              nmat=_NMAT, coeff_tail=None)
        self.assertEqual(str(cm.exception),
                         "_so_bond_green: coeff_tail must be the resolved "
                         "[eliashberg] bond_coeff_tail, got None")

    def test_mu_matches_non_so(self):
        """A band G (x) 1_2: mu of the spin-orbital run (filling per orbital
        per spin, counted over nd = 2 norb_phys spin-orbitals) equals mu of
        the non-spin-orbital run of the same band."""
        fx = F.SU2_4X4
        _, eps, _, out = _green_case(fx)
        hr = {(R, (0, 0)): complex(t) for R, t in fx.hopping}
        k = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in fx.shape]
        w, _ = sc._calc_eigenvalues(sc._build_hamiltonian_k(*k, hr, 1))
        mu_non_so = sc._determine_mu(w, _BETA, fx.filling, 1)
        self.assertLessEqual(abs(out.mu - mu_non_so), 1e-13)
        # counting the spin-orbital eigenvalues over one orbital misses
        self.assertGreater(abs(sc._determine_mu(out.eigenvalues, _BETA, fx.filling, 1)
                               - mu_non_so), 1e-3)

    def test_band_order_from_so_epsilon_k(self):
        """The eigenvalues are those of the spin-block band H(k) built from
        the fixture's literal lists (numpy.linalg.eigvalsh, test side), in
        ascending order at each k, to round-off: for a spin-flip band (full
        diagonalization) and an S_z-conserving one (per sector)."""
        for fx in (F.SOC_ALL_TYPES, F.ZEEMAN_EXTERN):
            with self.subTest(fixture=fx.name):
                out = _green_case(fx)[3]
                self.assertIs(out.conserving, fx is F.ZEEMAN_EXTERN)
                ref = np.linalg.eigvalsh(F._band(fx))           # (Nx, Ny, Nz, nd)
                self.assertEqual(out.eigenvalues.shape, ref.shape)
                self.assertLessEqual(np.abs(out.eigenvalues - ref).max(),
                                     1e-13 * np.abs(ref).max())
                self.assertTrue(np.all(np.diff(out.eigenvalues, axis=-1) >= 0))


class TestWriteCase(unittest.TestCase):
    """The case writer of the public-entry fixtures: overrides replace keys
    case-insensitively section by section, and the non-spin-orbital files
    of a spin-independent band."""

    def _case(self, fixture, **kw):
        d = tempfile.mkdtemp(prefix="sc_so_bond_case_")
        self.addCleanup(shutil.rmtree, d, True)
        return F.write_case(d, fixture, **kw)

    @staticmethod
    def _spellings(section, key):
        return [k for k in section if k.lower() == key.lower()]

    def test_mixed_case_override_replaces_default(self):
        inp = self._case(F.ZEEMAN_EXTERN,
                         param={"nmat": 16, "coeff_extern": 0.1, "FILLING": 0.3},
                         eliashberg={"bond_channels": False, "SOLVER_MODE": "both"},
                         mode={"enable_spin_orbital": True})
        param, eli, mode = inp["mode"]["param"], inp["eliashberg"], inp["mode"]
        for section, key, spelled, value in (
                (param, "nmat", "Nmat", 16), (param, "coeff_extern", "COEFF_EXTERN", 0.1),
                (param, "filling", "filling", 0.3),
                (eli, "bond_channels", "Bond_Channels", False),
                (eli, "solver_mode", "Solver_Mode", "both"),
                (mode, "enable_spin_orbital", "Enable_Spin_Orbital", True)):
            with self.subTest(key=key):
                self.assertEqual(self._spellings(section, key), [spelled])
                self.assertEqual(section[spelled], value)
        self.assertEqual(eli["Frequency"], "static")

    def test_coeff_tail_lands_in_mode_param(self):
        inp = self._case(F.SOC_ALL_TYPES, param={"coeff_tail": 0.5})
        self.assertEqual(inp["mode"]["param"]["coeff_tail"], 0.5)
        self.assertEqual(self._spellings(inp["eliashberg"], "coeff_tail"), [])
        self.assertEqual(self._spellings(inp["mode"], "coeff_tail"), [])
        with self.assertRaises(ValueError) as cm:
            self._case(F.SOC_ALL_TYPES, mode={"Param": {}})
        self.assertEqual(str(cm.exception),
                         "write_case: pass the [mode.param] overrides as param=")

    def test_non_so_files(self):
        fx = F.SU2_4X4
        so = self._case(fx)
        non_so = self._case(fx, mode={"ENABLE_SPIN_ORBITAL": False})
        self.assertEqual(self._spellings(non_so["mode"], "enable_spin_orbital"),
                         ["Enable_Spin_Orbital"])
        self.assertIs(sc._resolve_spin_orbital_flag(non_so), False)
        self.assertIs(sc._resolve_spin_orbital_flag(so), True)
        geom, hr, tables = sc._read_interaction_files(non_so)
        self.assertEqual(geom["norb"], 1)
        self.assertEqual(hr, {(R, (0, 0)): complex(t) for R, t in fx.hopping})
        geom_so, hr_so, tables_so = sc._read_interaction_files(so)
        self.assertEqual(geom_so["norb"], 2)
        self.assertEqual(len(hr_so), 2 * len(hr))
        self.assertEqual(tables, tables_so)
        self.assertEqual(sorted(tables), sorted(t for t, _ in fx.interactions))
        for fx in (F.SOC_ALL_TYPES, F.ZEEMAN_EXTERN, F.ZEEMAN_TRANSFER):
            with self.subTest(fixture=fx.name):
                with self.assertRaises(ValueError) as cm:
                    self._case(fx, mode={"enable_spin_orbital": False})
                self.assertEqual(
                    str(cm.exception),
                    "write_case: fixture {} has spin-dependent Transfer entries or an "
                    "Extern table; it has no non-spin-orbital files".format(fx.name))


if __name__ == "__main__":
    unittest.main()
