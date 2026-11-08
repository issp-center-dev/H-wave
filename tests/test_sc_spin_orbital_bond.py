"""hwave_sc in spin-orbital mode with bond-resolved interaction channels
(issue #270): the configuration of the static spin-orbital bond path and
the refusals that come before any file is read."""

import contextlib
import functools
import gc
import logging
import os
import re
import shutil
import tempfile
import tracemalloc
import unittest
import weakref
from unittest import mock

import numpy as np
from requests.structures import CaseInsensitiveDict

import hwave.sc as sc
# eliashberg_bond is imported here (not first inside a traced phase): the
# phase functions import it lazily, and a first import under tracemalloc
# would count the module's code objects in the phase that triggers it
from hwave.solver import (backend, bond_channels, eliashberg_bond, fingerprint, host_memory,
                          so_bond, so_pairing, spin_orbital)
from hwave.sc import (_SO_BOND_CHI0Q_UNUSED, _SO_BOND_DIAGNOSTICS_REFUSAL,
                      _SO_BOND_FLEX_KEYS_UNUSED, _SO_BOND_FLEX_REFUSAL,
                      _SO_BOND_GPU_IGNORED, _SO_BOND_GPU_REQUIRED_REFUSAL,
                      _SO_BOND_GREEN_REFUSAL, _SO_BOND_LAYOUT_IGNORED,
                      _SO_BOND_PRECONDITION_REFUSAL)
from tests import so_bond_entry_fixtures as F
from tests import so_zero_field as Z
from tests.so_bond_oracles import CELL, closed_row, rel


class TestDynamicRefusal(unittest.TestCase):
    """bond_channels = true with frequency = "dynamic" has its own refusal,
    raised by the configuration guard before any file is read (a mode and an
    eliashberg section suffice: a later check would raise KeyError); with
    frequency = "static" the guard accepts it."""

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
        # the static row is accepted by the configuration guard
        so_pairing.check_so_config({"mode": {"Enable_Spin_Orbital": True},
                                    "eliashberg": {"Bond_Channels": True}})


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
_NMAT, _BETA = 32, 2.0


def _read_only(a):
    if isinstance(a, np.ndarray):
        a.flags.writeable = False
    return a


@functools.lru_cache(maxsize=None)
def _green_case(fixture, layout="auto"):
    """The Green function of the spin-orbital bond path of ``fixture``,
    built through the run's own readers (:func:`tests.so_bond_entry_fixtures.prepare`,
    ``so_bond_layout = layout``) and ``_so_bond_green``: ``(so_run,
    epsilon_k, config, _SOBondGreen)``. Cached by content (the fixture is
    hashable); every array is read-only."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_green_")
    try:
        p = F.prepare(fixture, d, eliashberg={"SO_Bond_Layout": layout})
    finally:
        shutil.rmtree(d, True)
    out = sc._so_bond_green(p.epsilon_k, p.so_run, norb_phys=p.norb_phys, filling=p.filling,
                            beta=p.beta, nmat=p.nmat, coeff_tail=p.cfg.coeff_tail)
    for a in (p.epsilon_k, out.eigenvalues, out.carrier.full_sc, out.carrier.deflated_kw,
              out.carrier.green0_tail):
        _read_only(a)
    return p.so_run, p.epsilon_k, p.cfg, out


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

    def test_non_so_eliashberg_keys_lower_case(self):
        """With spin-orbital mode off every [eliashberg] key is written in
        lower case (the non-spin-orbital reader looks keys up
        case-sensitively); spin-orbital mode keeps the mixed spellings."""
        eli = {"Pairing_Type": "triplet", "G2_Tail": False, "bond_coeff_tail": 0.5}
        non_so = self._case(F.SU2_4X4, mode={"enable_spin_orbital": False}, eliashberg=eli)
        self.assertEqual(non_so["eliashberg"], {
            "frequency": "static", "bond_channels": True, "solver_mode": "eigenvalue",
            "num_eigenvalues": 2, "pairing_type": "triplet", "g2_tail": False,
            "bond_coeff_tail": 0.5})
        so = self._case(F.SU2_4X4, eliashberg=eli)
        self.assertEqual(sorted(so["eliashberg"]), sorted([
            "Frequency", "Bond_Channels", "Solver_Mode", "Num_Eigenvalues", "Pairing_Type",
            "G2_Tail", "bond_coeff_tail"]))

    def test_prepare_is_spin_orbital_only(self):
        d = tempfile.mkdtemp(prefix="sc_so_bond_case_")
        self.addCleanup(shutil.rmtree, d, True)
        with self.assertRaises(ValueError) as cm:
            F.prepare(F.SU2_4X4, d, mode={"Enable_Spin_Orbital": False})
        self.assertEqual(str(cm.exception),
                         "prepare: spin-orbital runs only ([mode] enable_spin_orbital = true)")

    def test_run_case_key(self):
        """The run cache key: overrides frozen case-insensitively, None and
        {} alike, two spellings of one key refused; a layout given twice is
        refused before anything runs."""
        self.assertEqual(F.freeze({"Nmat": 8}), F.freeze({"nmat": 8}))
        self.assertEqual(F.freeze(None), F.freeze({}))
        self.assertEqual(F.freeze({"CellShape": [4, 4, 1], "B": {"X": [1]}}),
                         (("b", (("x", (1,)),)), ("cellshape", (4, 4, 1))))
        self.assertNotEqual(F.freeze({"enable_spin_orbital": False}), F.freeze(None))
        with self.assertRaises(ValueError) as cm:
            F.freeze({"Nmat": 8, "NMAT": 8})
        self.assertEqual(str(cm.exception),
                         "freeze: the keys 'Nmat' and 'NMAT' differ only in case")
        with mock.patch.object(sc, "calc_eliashberg") as run:
            with self.assertRaises(ValueError) as cm:
                F.run_case(F.ZEEMAN_EXTERN, layout="dense",
                           eliashberg={"so_bond_layout": "auto"})
        run.assert_not_called()
        self.assertEqual(str(cm.exception),
                         "run_case: so_bond_layout is set both by layout = 'dense' and in "
                         "eliashberg")

# ------------------------------------------------------------ the static operator
def _bond_off_sector_mask(nd, B):
    """The S_z-off-sector entries of one ``(ND, ND)`` code-slot matrix, from
    the slot labels: row ``(alpha; A, B')``, column ``(beta; C, D)``, off
    where ``S(A) + S(D) != S(B') + S(C)`` (``S(X) = X // (nd / 2)``)."""
    n = nd // 2
    slots = [(m, x, y) for m in range(B) for x in range(nd) for y in range(nd)]
    S = [x // n for x in range(nd)]
    return np.array([[S[a] + S[d] != S[b] + S[c] for (_, c, d) in slots]
                     for (_, a, b) in slots])


_BRIDGE_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond.bridge")
#: the bridge runs pin the configuration resolution: g2_tail and
#: bond_coeff_tail are set away from their defaults
_BRIDGE_ELI = {"G2_Tail": False, "Bond_Coeff_Tail": 0.5}
_BRIDGE = {"dense": F.SOC_ALL_TYPES, "sz_packed": F.ZEEMAN_EXTERN}
_BRIDGE_TOL = 1e-12
_MUTANT_MISS = 1e-6


def _prepare(fixture, eliashberg=None, param=None):
    """:func:`tests.so_bond_entry_fixtures.prepare` of ``fixture`` with the
    bridge configuration (and ``eliashberg`` on top), in a directory removed
    on return."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_bridge_")
    try:
        return F.prepare(fixture, d, eliashberg=dict(_BRIDGE_ELI, **(eliashberg or {})),
                         param=param)
    finally:
        shutil.rmtree(d, True)


@functools.lru_cache(maxsize=None)
def _bridge_prepared(fixture):
    return _prepare(fixture)


def _green_of(p):
    return sc._so_bond_green(p.epsilon_k, p.so_run, norb_phys=p.norb_phys,
                             filling=p.filling, beta=p.beta, nmat=p.nmat,
                             coeff_tail=p.cfg.coeff_tail)


def _admit(p, method="arnoldi"):
    return sc._admit_so_bond_static(
        p.so_run, cfg=p.cfg, cell=p.cell, nmat=p.nmat, nd=2 * p.norb_phys, shape=p.shape,
        solver_mode=p.solver_mode, num_eigenvalues=p.num_eigenvalues, max_iter=p.max_iter,
        spectral_shift=p.spectral_shift, eigenvalue_method=method, log=_BRIDGE_LOG)


def _operator(p, so_green, adm):
    return sc._build_so_bond_operator(p.so_run, adm, so_green, beta=p.beta, nmat=p.nmat,
                                      g2_tail=p.cfg.g2_tail, norb_phys=p.norb_phys)


def _chain(p):
    """The production chain of a prepared run: ``(so_green, admission,
    static)``. A patch that must reach it encloses this whole call."""
    so_green = _green_of(p)
    adm = _admit(p)
    return so_green, adm, _operator(p, so_green, adm)


def _apply(static, v):
    """The production operator on a full flat gap (both layouts)."""
    if static.layout == "sz_packed":
        return static.kernel.matvec(v)
    return static.operator.matvec(v)


def _layout_of(p):
    return so_bond.select_layout(p.so_run.context.compiled, p.so_run.band.transfer,
                                 override=p.cfg.layout_request)[0]


def _gaps(fixture, count=2, seed=11):
    """``count`` seeded random fermion-antisymmetric gaps (full flat)."""
    basis = sc._SOAntisymmetricBasis((fixture.nd, fixture.nd) + tuple(fixture.shape))
    rng = np.random.default_rng(seed)
    return [basis.to_full(rng.standard_normal(basis.dim) + 1j * rng.standard_normal(basis.dim))
            for _ in range(count)]


def _generic_gap(fixture, seed=13):
    n = fixture.nd ** 2 * fixture.nvol
    rng = np.random.default_rng(seed)
    return rng.standard_normal(n) + 1j * rng.standard_normal(n)


def _miss(apply_a, apply_b, gaps):
    """``max_j |A v_j - B v_j| / max_j |B v_j|`` (2-norms)."""
    num = den = 0.0
    for v in gaps:
        b = apply_b(v)
        num = max(num, float(np.linalg.norm(apply_a(v) - b)))
        den = max(den, float(np.linalg.norm(b)))
    return num / den


def _test_side_split(G, coeff_tail, beta):
    """The deflated Green function ``(1, nmat, nvol, nd, nd)`` and its tail
    of a test-side ``G (nd, nd, Nx, Ny, Nz, nmat)``: ``G - c / (i w_n)`` and
    ``c beta / 2`` on the diagonal (the band's eigenvectors are unitary),
    ``None`` for ``c = 0``. This re-implements, on ``G``, the documented
    split of the bubble's inputs (the one ``build_green`` makes on the
    eigenvectors), without importing it: a production deviation from that
    contract is caught, a defect in the contract itself is outside this
    bridge."""
    nd, nmat = G.shape[0], G.shape[-1]
    nvol = int(np.prod(G.shape[2:5]))
    g_kw = np.ascontiguousarray(G.reshape(nd, nd, nvol, nmat).transpose(3, 2, 0, 1))
    if coeff_tail == 0.0:
        return g_kw[None], None
    iw = 1j * F._matsubara(nmat, beta)
    eye = np.eye(nd)
    deflated = g_kw - coeff_tail * eye[None, None] / iw[:, None, None, None]
    tail = np.broadcast_to(0.5 * coeff_tail * beta * eye, (1, nmat, nvol, nd, nd)).astype(complex)
    return deflated[None], tail


def _reference_vertex(p, G, *, coeff_tail):
    """``(V_bare, V_fluct, geom)`` (dense) of the run's interaction on the
    test-side Green function ``G`` through the validated array functions."""
    from tests.so_bond_fixtures import NOCAP
    fx = p.fixture
    nd, shape, nvol = fx.nd, tuple(fx.shape), fx.nvol
    compiled, topo = p.so_run.context.compiled, p.so_run.bond.topo
    geom = so_bond.slot_geometry(topo, shape)
    plan = so_bond.plan_so_bond(nmat=p.nmat, nvol=nvol, B=geom.B, nd=nd,
                                consumers=("pairing_static",), caps=NOCAP, layout="dense")
    gamma = so_bond.so_bond_vertex(compiled, topo, shape, plan=plan)
    deflated, tail = _test_side_split(G, coeff_tail, p.beta)
    ND = geom.B * nd * nd
    X = so_bond.so_bond_bubble(deflated, tail, p.beta, geom, shape, plan=plan, freq="static",
                               sink=np.zeros((nvol, ND, ND), dtype=complex))
    chi = so_bond.so_bond_dress(X[None], gamma, plan=plan, l0=p.nmat // 2,
                                spatial_shape=shape)[0]
    v_bare, v_fluct = so_bond.so_bond_pairing_vertex_parts(gamma, chi, geom, plan=plan)
    return v_bare, v_fluct, geom


def _reference_kernel_of(p, V, G, geom, *, g2_tail):
    from tests.so_bond_fixtures import pair_kernel
    fx = p.fixture
    G2 = bond_channels._g2_from_green(G, p.beta, tail=g2_tail)[..., None]
    return pair_kernel(V[None], G2, geom, fx.nd, tuple(fx.shape), p.beta)


def _bridge_config(p):
    """``(g2_tail, bond_coeff_tail)`` as written in the run's input."""
    eli = CaseInsensitiveDict(p.input_dict["eliashberg"])
    return bool(eli["g2_tail"]), float(eli["bond_coeff_tail"])


_REF_CACHE = {}
_REF_COUNT = {"hit": 0, "miss": 0}


def _ref_key(p):
    """The cache key of :func:`_reference_kernel`. ``K_ref`` is always the
    dense kernel; the run's layout is in the key only so that the run of
    each bridge fixture keeps its own entry."""
    g2_tail, coeff_tail = _bridge_config(p)
    return (p.fixture, _layout_of(p), p.nmat, p.beta, g2_tail, coeff_tail)


def _reference_kernel(p):
    """``K_ref`` of a prepared run: the dense production kernel of the pairing
    vertex built from :func:`tests.so_bond_entry_fixtures.independent_green`
    with the run's written ``g2_tail`` / ``bond_coeff_tail``. Cached per
    :func:`_ref_key` (hits and misses counted)."""
    key = _ref_key(p)
    if key in _REF_CACHE:
        _REF_COUNT["hit"] += 1
        return _REF_CACHE[key]
    _REF_COUNT["miss"] += 1
    g2_tail, coeff_tail = key[4], key[5]
    G, _ = F.independent_green(p.fixture, nmat=p.nmat, beta=p.beta)
    v_bare, v_fluct, geom = _reference_vertex(p, G, coeff_tail=coeff_tail)
    K = _reference_kernel_of(p, v_bare + v_fluct, G, geom, g2_tail=g2_tail)
    _REF_CACHE[key] = K
    return K


def bridge_b_miss(p, static):
    """The input-independent bridge: ``max_j |K_prod v_j - K_ref v_j| / max_j
    |K_ref v_j|`` over two seeded random antisymmetric gaps."""
    K = _reference_kernel(p)
    return _miss(lambda v: _apply(static, v), K.matvec, _gaps(p.fixture))


class _Spy:
    """Wraps ``real``: records the arguments and the return value of every
    call."""

    def __init__(self, real):
        self.real = real
        self.calls = []

    def __call__(self, *args, **kw):
        out = self.real(*args, **kw)
        self.calls.append((args, kw, out))
        return out


def l0_seen(spy):
    """The ``l0`` values the real ``so_bond_dress`` received (``spy``, a
    :class:`_Spy`, wraps it)."""
    return [kw["l0"] for _, kw, _ in spy.calls]


class _PartsSpy:
    """Wraps ``so_bond.so_bond_pairing_vertex_parts``: records its inputs and
    copies of its outputs as returned (the caller combines them in place
    afterwards), and the returned objects themselves."""

    def __init__(self):
        self.real = so_bond.so_bond_pairing_vertex_parts
        self.calls = []

    @staticmethod
    def _copy(a):
        if isinstance(a, so_bond.PackedArrays):
            return so_bond.PackedArrays(a.sector_map, tuple(b.copy() for b in a.blocks))
        return a.copy()

    def __call__(self, gamma_b, chi_b, geom, **kw):
        out = self.real(gamma_b, chi_b, geom, **kw)
        self.calls.append(dict(gamma_b=gamma_b, chi_b=chi_b, out=out,
                               copies=tuple(self._copy(a) for a in out)))
        return out


def _dense(a):
    return a.to_dense() if isinstance(a, so_bond.PackedArrays) else a


def _wiring_kernel(p, parts, so_green):
    """The validated kernel of the production parts: ``pair_kernel((V_bare +
    V_fluct)[None], G2(carrier.full_sc))`` (dense)."""
    vb, vf = (_dense(a) for a in parts.calls[-1]["copies"])
    geom = so_bond.slot_geometry(p.so_run.bond.topo, p.shape)
    return _reference_kernel_of(p, vb + vf, so_green.carrier.full_sc, geom,
                                g2_tail=p.cfg.g2_tail)


class TestKernelBridge(unittest.TestCase):
    """Gate A2-P. (a) wiring: the production operator equals the validated
    kernel built from the same production parts; (b) input-independent: it
    equals the kernel built from a test-side Green function through the
    validated array functions. The mutation table shows which oracle
    catches which wrong input.

    Blind spots and the gates that cover them: both bridges take the
    production Gamma_B, so a vertex filed at -q passes them (the SU(2)
    gate's non-vacuity test); a Zeeman sign flipped in production and in
    the test-side Green function alike passes both (the analytic Zeeman
    oracle of the Green-function tests)."""

    @classmethod
    def setUpClass(cls):
        # the bridge references are built before any production function is
        # patched (row 0 of the mutation table relies on it)
        for fx in _BRIDGE.values():
            _reference_kernel(_bridge_prepared(fx))

    def test_dense_wiring(self):
        p = _bridge_prepared(F.SOC_ALL_TYPES)
        parts, sv, sd = _PartsSpy(), _Spy(so_bond.so_bond_vertex), _Spy(so_bond.so_bond_dress)
        with mock.patch.object(so_bond, "so_bond_vertex", sv), \
                mock.patch.object(so_bond, "so_bond_dress", sd), \
                mock.patch.object(so_bond, "so_bond_pairing_vertex_parts", parts):
            so_green, adm, static = _chain(p)
        self.assertEqual((adm.layout, static.layout), ("dense", "dense"))
        self.assertIsNone(static.kernel)
        self.assertEqual(static.basis.shape, (2, 2) + p.shape)
        self.assertEqual(static.operator.shape, (4 * p.fixture.nvol,) * 2)
        self.assertEqual((len(sv.calls), len(sd.calls), len(parts.calls)), (1, 1, 1))
        # the parts received the vertex and the dressed susceptibility
        self.assertIs(parts.calls[0]["gamma_b"], sv.calls[0][2])
        self.assertIs(sd.calls[0][0][1], sv.calls[0][2])
        self.assertTrue(np.shares_memory(parts.calls[0]["chi_b"], sd.calls[0][2]))
        K = _wiring_kernel(p, parts, so_green)
        gaps = _gaps(p.fixture) + [_generic_gap(p.fixture)]
        self.assertLessEqual(_miss(static.operator.matvec, K.matvec, gaps), _BRIDGE_TOL)
        # the spin-flip band breaks S_z: both residuals are visible
        self.assertGreater(static.r_G, 1e-3)
        self.assertGreater(static.r_V, 1e-6)
        self.assertEqual(static.r_V, so_pairing.bond_sz_conservation_residual(
            parts.calls[0]["out"][0], nd=2, B=adm.geom.B))

    def test_dense_input_independent(self):
        p = _bridge_prepared(F.SOC_ALL_TYPES)
        _, _, static = _chain(p)
        miss = bridge_b_miss(p, static)
        self.assertLessEqual(miss, _BRIDGE_TOL, "bridge (b) dense miss = %.3e" % miss)

    def test_packed_sector_ops_are_the_dense_kernel(self):
        p = _bridge_prepared(F.ZEEMAN_EXTERN)
        parts = _PartsSpy()
        with mock.patch.object(so_bond, "so_bond_pairing_vertex_parts", parts):
            so_green, adm, static = _chain(p)
        self.assertEqual((adm.layout, static.layout), ("sz_packed", "sz_packed"))
        self.assertIsNone(static.operator)
        self.assertIsNone(static.basis)
        self.assertEqual((static.r_V, static.r_G), (0.0, 0.0))
        K = _wiring_kernel(p, parts, so_green)
        n = p.fixture.nd ** 2 * p.fixture.nvol
        rng = np.random.default_rng(17)
        for pp in static.kernel.pp_sectors:
            with self.subTest(pp=pp):
                op = static.kernel.sector_op(pp)
                idx = static.kernel.pp_index(pp)
                x = rng.standard_normal(len(idx)) + 1j * rng.standard_normal(len(idx))

                def embedded(x, idx=idx):
                    phi = np.zeros(n, dtype=complex)
                    phi[idx] = x
                    return K.matvec(phi)[idx]
                self.assertLessEqual(_miss(op.operator.matvec, embedded, [x]), _BRIDGE_TOL)
        miss = bridge_b_miss(p, static)
        self.assertLessEqual(miss, _BRIDGE_TOL)

    def test_non_vacuity(self):
        """Each of the off-site channel blocks of V, V_fluct and the
        spin-off-diagonal Green function moves the reference kernel by at
        least 1e-3 of the operator: the bridge sees them."""
        for layout, fx in _BRIDGE.items():
            p = _bridge_prepared(fx)
            _, adm, static = _chain(p)
            g2_tail, coeff_tail = _bridge_config(p)
            G, _ = F.independent_green(fx, nmat=p.nmat, beta=p.beta)
            v_bare, v_fluct, geom = _reference_vertex(p, G, coeff_tail=coeff_tail)
            npair = fx.nd ** 2
            onsite = [m for m, d in enumerate(geom.delta_r) if not np.any(d)]
            self.assertEqual(len(onsite), 1)
            keep = np.zeros(geom.B * npair, dtype=bool)
            keep[onsite[0] * npair:(onsite[0] + 1) * npair] = True
            V = v_bare + v_fluct
            onsite_only = np.where(keep[:, None] & keep[None, :], V, 0.0)
            variants = {"off-site blocks": (onsite_only, G), "V_fluct": (v_bare, G)}
            if layout == "dense":
                n = fx.norb_phys
                G0 = G.copy()
                G0[:n, n:] = 0.0
                G0[n:, :n] = 0.0
                self.assertGreater(np.abs(G[:n, n:]).max(), 1e-3)
                vb0, vf0, _ = _reference_vertex(p, G0, coeff_tail=coeff_tail)
                variants["spin-off-diagonal G"] = (vb0 + vf0, G0)
            for name, (Vx, Gx) in variants.items():
                with self.subTest(layout=layout, zeroed=name):
                    K = _reference_kernel_of(p, Vx, Gx, geom, g2_tail=g2_tail)
                    miss = _miss(lambda v: _apply(static, v), K.matvec, _gaps(fx))
                    self.assertGreaterEqual(miss, 1e-3, "non-vacuity %s / %s: miss = %.3e"
                                            % (layout, name, miss))

    def test_v_combined_in_place(self):
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                so_green = _green_of(p)
                adm = _admit(p)
                parts = _PartsSpy()
                with mock.patch.object(so_bond, "so_bond_pairing_vertex_parts", parts):
                    gamma_b, V = sc._so_bond_static_vertex(p.so_run, adm, so_green,
                                                           beta=p.beta, nmat=p.nmat)
                v_bare, v_fluct = parts.calls[0]["out"]
                vb0, vf0 = parts.calls[0]["copies"]
                self.assertIs(V, v_bare)
                if layout == "sz_packed":
                    for b, ref, a, f in zip(V.blocks, v_bare.blocks, vb0.blocks, vf0.blocks):
                        self.assertTrue(np.shares_memory(b, ref))
                        np.testing.assert_allclose(b, a + f, rtol=0, atol=1e-15)
                else:
                    self.assertTrue(np.shares_memory(V, v_bare))
                    np.testing.assert_allclose(V, vb0 + vf0, rtol=0, atol=1e-15)

    def test_vertex_from_the_entry_compile(self):
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                sv = _Spy(so_bond.so_bond_vertex)
                with mock.patch.object(so_bond, "so_bond_vertex", sv):
                    _chain(p)
                self.assertEqual(len(sv.calls), 1)
                self.assertIs(sv.calls[0][0][0], p.so_run.context.compiled)
                self.assertIs(sv.calls[0][0][1], p.so_run.bond.topo)

    def test_dress_receives_l0(self):
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                sd = _Spy(so_bond.so_bond_dress)
                with mock.patch.object(so_bond, "so_bond_dress", sd):
                    _chain(p)
                self.assertEqual(l0_seen(sd), [p.nmat // 2])
                X = sd.calls[0][0][0]
                if layout == "sz_packed":
                    self.assertIsInstance(X, so_bond.PackedArrays)
                    self.assertEqual(X.lead_shape, (1, fx.nvol))
                else:
                    self.assertEqual(X.shape[0], 1)
                    self.assertEqual(X.ndim, 4)

    @staticmethod
    def _stiff(real):
        """A wrapper calling ``real`` with the conditioning guard at
        ``cond_tol = inf`` (every block refused)."""
        def wrapper(X, gamma_b, **kw):
            kw["cond_tol"] = np.inf
            return real(X, gamma_b, **kw)
        return wrapper

    def test_conditioning_failure_names_l0(self):
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                with mock.patch.object(so_bond, "so_bond_dress",
                                       self._stiff(so_bond.so_bond_dress)):
                    with self.assertRaises(bond_channels.BondConditioningError) as cm:
                        _chain(p)
                self.assertEqual(cm.exception.l, p.nmat // 2)
                self.assertIn("grid index l={}".format(p.nmat // 2), str(cm.exception))

    def test_method_guard(self):
        text = ("[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
                "eigenvalue_method = 'lanczos' is not available on the spin-orbital bond "
                "path; only \"arnoldi\" is. Set eigenvalue_method = \"arnoldi\" (or remove "
                "the key).")
        self.assertEqual(sc._SO_BOND_METHOD_REFUSAL.format(method=repr("lanczos")), text)
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                _green_of(p)
                sv, sp = _Spy(so_bond.so_bond_vertex), _Spy(so_bond.plan_so_bond)
                with mock.patch.object(so_bond, "so_bond_vertex", sv), \
                        mock.patch.object(so_bond, "plan_so_bond", sp):
                    with self.assertLogs(_BRIDGE_LOG, "INFO") as logs:
                        with self.assertRaises(ValueError) as cm:
                            _admit(p, method="lanczos")
                self.assertEqual(str(cm.exception), text)
                self.assertEqual((len(sv.calls), len(sp.calls)), (0, 0))
                # the layout was selected (its reason logged) before the refusal
                self.assertTrue(any(m.startswith("INFO:%s:so bond layout: %s"
                                                 % (_BRIDGE_LOG.name, layout))
                                    for m in logs.output), logs.output)

    def test_admission_rows_take_the_shift(self):
        """The eigensolver rows follow the run's spectral_shift: the admission
        passes it (and max_iter) to the static rows, and the plan's
        kernel_static eigensolver row is the one of that shift (dense SOC,
        the ARPACK branch of the antisymmetric dimension 28)."""
        rows = {}
        for shift in (None, "auto", 0.5):
            with self.subTest(spectral_shift=shift):
                p = _prepare(F.SOC_ALL_TYPES,
                             eliashberg=None if shift is None else {"Spectral_Shift": shift})
                self.assertEqual(p.spectral_shift, shift)
                spy = _Spy(so_bond.static_so_bond_rows)
                with mock.patch.object(so_bond, "static_so_bond_rows", spy):
                    adm = _admit(p)
                self.assertEqual(len(spy.calls), 1)
                kw = spy.calls[0][1]
                self.assertEqual(kw["spectral_shift"], shift)
                self.assertIs(type(kw["spectral_shift"]), type(shift))
                self.assertEqual(kw["max_iter"], p.max_iter)
                self.assertEqual(adm.layout, "dense")
                want = so_bond.static_so_bond_rows(**kw)
                row = [r for r in adm.plan.rows if r[:2] == ("kernel_static", "eigensolver")]
                self.assertEqual(row, [r for r in want if r[:2] == ("kernel_static",
                                                                     "eigensolver")])
                self.assertEqual(len(row), 1)
                rows[shift] = row[0][2]
        self.assertNotEqual(rows["auto"], rows[0.5])
        self.assertNotEqual(rows[None], rows[0.5])

    def test_carrier_read_only(self):
        """After the operator chain the carrier arrays are still read-only and
        unchanged (no write in place anywhere in the chain)."""
        for layout, fx in _BRIDGE.items():
            with self.subTest(layout=layout):
                p = _bridge_prepared(fx)
                so_green = _green_of(p)
                c = so_green.carrier
                names = [n for n in ("full_sc", "deflated_kw", "green0_tail")
                         if getattr(c, n) is not None]
                self.assertEqual(len(names), 3)
                before = {n: getattr(c, n).copy() for n in names}
                _operator(p, so_green, _admit(p))
                for n in names:
                    a = getattr(c, n)
                    self.assertFalse(a.flags.writeable, n)
                    self.assertEqual(a.tobytes(), before[n].tobytes(), n)
                self.assertEqual(so_green.fp_green, sc._so_green_fp(c.full_sc))

    def test_bond_sz_residual(self):
        nd, B, nvol = 4, 3, 5
        ND = B * nd * nd
        off = _bond_off_sector_mask(nd, B)
        self.assertTrue(off.any() and (~off).any())
        rng = np.random.default_rng(7)
        V = (rng.standard_normal((nvol, ND, ND))
             + 1j * rng.standard_normal((nvol, ND, ND)))
        V[:, off] = 0.0
        self.assertEqual(so_pairing.bond_sz_conservation_residual(V, nd=nd, B=B), 0.0)
        vmax = float(np.abs(V).max())
        i, j = np.argwhere(off)[17]
        V[3, i, j] = 1e-3 * vmax
        r = so_pairing.bond_sz_conservation_residual(V, nd=nd, B=B)
        self.assertLessEqual(abs(r - 1e-3), 1e-15)
        # the residual reads every momentum and every channel block
        V[3, i, j] = 0.0
        last = np.argwhere(off)[-1]
        V[nvol - 1, last[0], last[1]] = -2e-3j * vmax
        self.assertLessEqual(
            abs(so_pairing.bond_sz_conservation_residual(V, nd=nd, B=B) - 2e-3), 1e-15)
        self.assertEqual(so_pairing.bond_sz_conservation_residual(
            np.zeros((nvol, ND, ND), dtype=complex), nd=nd, B=B), 0.0)
        V[1, 0, 0] = np.nan
        self.assertTrue(np.isnan(so_pairing.bond_sz_conservation_residual(V, nd=nd, B=B)))
        with self.assertRaises(ValueError):
            so_pairing.bond_sz_conservation_residual(V[:, :-1], nd=nd, B=B)

    # -------------------------------------------------------- mutation table
    def _row(self, label, fx, patches, *, build=None, reprepare=False):
        """One row: the reference is a cache hit taken before the patch;
        ``patches`` (a list of ``mock.patch.object`` arguments) enclose the
        whole production chain (``prepare`` too when ``reprepare``)."""
        p0 = _bridge_prepared(fx)
        before = dict(_REF_COUNT)
        _reference_kernel(p0)
        self.assertEqual(_REF_COUNT["miss"], before["miss"], label)
        self.assertEqual(_REF_COUNT["hit"], before["hit"] + 1, label)
        with contextlib.ExitStack() as stack:
            for target, name, new in patches:
                stack.enter_context(mock.patch.object(target, name, new))
            p = _prepare(fx) if reprepare else p0
            so_green = _green_of(p)
            if build is not None:
                build(p, so_green)
            adm = _admit(p)
            static = _operator(p, so_green, adm)
        if reprepare:
            self.assertEqual(_ref_key(p), _ref_key(p0))
        return p, so_green, static

    def _stiff_l0(self, p, patches):
        """The conditioning index under ``patches`` (``cond_tol = inf`` at
        the real function)."""
        with contextlib.ExitStack() as stack:
            for target, name, new in patches:
                stack.enter_context(mock.patch.object(target, name, new))
            with self.assertRaises(bond_channels.BondConditioningError) as cm:
                so_green = _green_of(p)
                _operator(p, so_green, _admit(p))
        return cm.exception.l

    def test_mutation_table(self):
        real_bubble = so_bond.so_bond_bubble
        real_g2 = bond_channels._g2_from_green
        real_green = sc._build_bond_green
        real_dress = so_bond.so_bond_dress
        real_extern = sc._read_so_extern
        nmat = _NMAT
        report = []

        def note(label, name, value):
            report.append("%s: %s = %.3e" % (label, name, value))

        # row 0: every oracle passes in both layouts
        for layout, fx in _BRIDGE.items():
            with self.subTest(row="0", layout=layout):
                spy = _Spy(real_dress)
                p, so_green, static = self._row("0/" + layout, fx,
                                                [(so_bond, "so_bond_dress", spy)])
                miss = bridge_b_miss(p, static)
                note("0/" + layout, "bridge (b)", miss)
                self.assertLessEqual(miss, _BRIDGE_TOL)
                self.assertEqual(l0_seen(spy), [nmat // 2])
                spy2 = _Spy(real_dress)
                self.assertEqual(self._stiff_l0(p, [(so_bond, "so_bond_dress",
                                                     self._stiff(spy2))]), nmat // 2)
                if fx is F.ZEEMAN_EXTERN:
                    z = F.zeeman_sign_miss(so_green.carrier.full_sc, so_green.mu, fixture=fx,
                                           h=F.ZEEMAN_H, beta=_BETA, nmat=nmat)
                    note("0/" + layout, "zeeman", z)
                    self.assertLessEqual(z, 1e-12)

        # M1: the bubble takes the Green function at mu + 0.1 (same band)
        with self.subTest(row="M1"):
            fx = F.SOC_ALL_TYPES
            alt, seen = {}, []

            def build(p, so_green):
                ev, evec = so_bond.canonical_band_eigen(p.epsilon_k, p.norb_phys,
                                                        conserving=so_green.conserving)
                alt["c"] = real_green(ev, evec, so_green.mu + 0.1, p.beta, p.nmat,
                                      p.cfg.coeff_tail, None, spin_orbital_mode=True)
                alt["prod"] = so_green.carrier

            def m1(*args, **kw):
                seen.append(args[0] is alt["prod"].deflated_kw
                            and args[1] is alt["prod"].green0_tail)
                return real_bubble(alt["c"].deflated_kw, alt["c"].green0_tail, *args[2:], **kw)
            p, _, static = self._row("M1", fx, [(so_bond, "so_bond_bubble", m1)], build=build)
            self.assertEqual(seen, [True])
            miss = bridge_b_miss(p, static)
            note("M1", "bridge (b)", miss)
            self.assertGreaterEqual(miss, _MUTANT_MISS)

        # M2: g2_tail ignored (forced on; the run sets false)
        with self.subTest(row="M2"):
            seen = []

            def m2(green, beta, tail=False):
                seen.append(tail)
                return real_g2(green, beta, tail=True)
            p, _, static = self._row("M2", F.SOC_ALL_TYPES,
                                     [(bond_channels, "_g2_from_green", m2)])
            self.assertEqual(seen, [False])
            miss = bridge_b_miss(p, static)
            note("M2", "bridge (b)", miss)
            self.assertGreaterEqual(miss, _MUTANT_MISS)

        # M3: bond_coeff_tail ignored (forced to 1.0; the run sets 0.5)
        with self.subTest(row="M3"):
            seen = []

            def m3(*args, **kw):
                seen.append(args[5])
                return real_green(*args[:5], 1.0, *args[6:], **kw)
            p, so_green, static = self._row("M3", F.SOC_ALL_TYPES,
                                            [(sc, "_build_bond_green", m3)])
            self.assertEqual(seen, [0.5])
            self.assertEqual(so_green.carrier.coeff_tail_applied, 1.0)
            miss = bridge_b_miss(p, static)
            note("M3", "bridge (b)", miss)
            self.assertGreaterEqual(miss, _MUTANT_MISS)

        # M4: the bare bubble in place of the dressed susceptibility; the
        # wiring bridge passes (it shares chi_B), bridge (b) does not
        for layout, fx in _BRIDGE.items():
            with self.subTest(row="M4", layout=layout):
                seen = []

                def m4(X, gamma_b, **kw):
                    seen.append(kw["l0"])
                    return X
                parts = _PartsSpy()
                p, so_green, static = self._row(
                    "M4/" + layout, fx, [(so_bond, "so_bond_dress", m4),
                                         (so_bond, "so_bond_pairing_vertex_parts", parts)])
                self.assertEqual(seen, [nmat // 2])
                wiring = _miss(lambda v: _apply(static, v),
                               _wiring_kernel(p, parts, so_green).matvec, _gaps(fx))
                note("M4/" + layout, "wiring bridge", wiring)
                self.assertLessEqual(wiring, _BRIDGE_TOL)
                miss = bridge_b_miss(p, static)
                note("M4/" + layout, "bridge (b)", miss)
                self.assertGreaterEqual(miss, _MUTANT_MISS)

        # M5: a wrong l0 (the spy inside the mutant records what the real
        # function receives)
        for layout, fx in _BRIDGE.items():
            with self.subTest(row="M5", layout=layout):
                seen = []

                def mutant(inner):
                    def m5(X, gamma_b, **kw):
                        seen.append(kw["l0"])
                        kw["l0"] = 0
                        return inner(X, gamma_b, **kw)
                    return m5
                spy = _Spy(real_dress)
                p, _, _ = self._row("M5/" + layout, fx, [(so_bond, "so_bond_dress",
                                                          mutant(spy))])
                self.assertEqual(seen, [nmat // 2])
                self.assertNotEqual(l0_seen(spy), [nmat // 2])
                self.assertEqual(l0_seen(spy), [0])
                l_err = self._stiff_l0(p, [(so_bond, "so_bond_dress",
                                            mutant(self._stiff(real_dress)))])
                note("M5/" + layout, "conditioning index", l_err)
                self.assertNotEqual(l_err, nmat // 2)

        # M6: the Zeeman field reversed at the Extern reader; prepare re-run
        with self.subTest(row="M6"):
            fx = F.ZEEMAN_EXTERN
            seen = []

            def m6(input_dict, norb_phys):
                table, name = real_extern(input_dict, norb_phys)
                seen.append(dict(table))
                return {k: -v for k, v in table.items()}, name
            p, so_green, static = self._row("M6", fx, [(sc, "_read_so_extern", m6)],
                                            reprepare=True)
            self.assertEqual(len(seen), 1)
            self.assertTrue(any(v != 0 for v in seen[0].values()))
            z = F.zeeman_sign_miss(so_green.carrier.full_sc, so_green.mu, fixture=fx,
                                   h=F.ZEEMAN_H, beta=_BETA, nmat=nmat)
            note("M6", "zeeman", z)
            self.assertGreaterEqual(z, 1e-3)
            miss = bridge_b_miss(p, static)
            note("M6", "bridge (b)", miss)
            self.assertGreaterEqual(miss, _MUTANT_MISS)


# ------------------------------------------------------------ the solve
#: the gap archive members of a static spin-orbital bond run, written here
#: literally (the provenance members of the stage come from fingerprint)
_GAP_BASE = frozenset(fingerprint.required_members("gap", "sc_so_bond_static")) | {
    "index_convention", "spin_orbital", "norb_phys", "kx", "ky", "kz", "spin_components",
    "gap", "gap_spin_weights", "eigen_converged", "pp_sector", "gap_status",
    "block_layout", "sz_weight_upup", "sz_weight_dndn", "sz_weight_singlet0",
    "sz_weight_triplet0", "sz_conservation_residual", "sz_conserved"}
_GAP_LISTING = frozenset({"eigenvalues", "eigenvectors", "eigenvector_spin_weights",
                          "eigenvalue_match", "pp_sectors"})
#: the member set of each (solver mode, gap_status) row whose set does not
#: depend on the listing length
_GAP_MEMBERS = {
    ("eigenvalue", "certified"): _GAP_BASE | _GAP_LISTING,
    ("eigenvalue", "uncertified"): _GAP_BASE | _GAP_LISTING,
    ("iteration", "certified"): _GAP_BASE,
    ("iteration", "uncertified"): _GAP_BASE,
    ("iteration", "none"): _GAP_BASE,
    ("both", "certified"): _GAP_BASE | _GAP_LISTING,
}
#: the eigenvalue.dat keys of every static spin-orbital bond run
_BOND_DAT_KEYS = ("bond_channels", "bond_n_channels", "bond_delta_r", "bond_max_shells",
                  "so_bond_layout", "so_bond_layout_reason", "bond_memory_cap", "g2_tail",
                  "bond_coeff_tail", "bond_coeff_tail_source", "approximation",
                  "gap_status", "pp_sector", "selection")
_SOLVE = {"dense": F.SOC_ALL_TYPES, "sz_packed": F.ZEEMAN_EXTERN}


def _dat_keys(text):
    """The ``# key = value`` lines of an eigenvalue file as a dict."""
    out = {}
    for line in text.splitlines():
        if line.startswith("# ") and " = " in line:
            k, v = line[2:].split(" = ", 1)
            out.setdefault(k, v)
    return out


def _solve_prepared(p, out_dir):
    """``sc._solve_so_bond_static`` on the prepared run ``p`` (its own
    [eliashberg] values; the outputs go to ``out_dir``)."""
    eli = p.eli_param
    return sc._solve_so_bond_static(
        p.input_dict, so_run=p.so_run, cfg=p.cfg, epsilon_k=p.epsilon_k, eli_param=eli,
        norb_phys=p.norb_phys, beta=p.beta, nmat=p.nmat, filling=p.filling, shape=p.shape,
        cell=p.cell, kx_array=p.kx, ky_array=p.ky, kz_array=p.kz, output_dir=out_dir,
        gap_file="gap.dat", eigenvalue_file="eigenvalue.dat", solver_mode=p.solver_mode,
        num_eigenvalues=p.num_eigenvalues, eigenvalue_method="arnoldi",
        max_iter=p.max_iter, alpha=eli.get("alpha", 0.5),
        tol=eli.get("convergence_tol", 1.0e-5),
        init_gap_mode=sc._resolve_init_gap(eli.get("init_gap"), "singlet"))


def _solve_case(test, fixture, eliashberg=None):
    """``(value, gap members, eigenvalue.dat text, path of gap.npz)`` of a
    function-level solve of ``fixture`` in a directory removed after the
    test."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_solve_")
    test.addCleanup(shutil.rmtree, d, True)
    p = F.prepare(fixture, d, eliashberg=eliashberg)
    out = os.path.join(d, "solve")
    value = _solve_prepared(p, out)
    path = os.path.join(out, "gap.npz")
    with np.load(path, allow_pickle=False) as z:
        members = {k: z[k] for k in z.files}
    with open(os.path.join(out, "eigenvalue.dat")) as f:
        text = f.read()
    return p, value, members, text, path


class _Calls:
    """A spy that records the arguments of each call and returns what the
    real function returns."""

    def __init__(self, real):
        self.real, self.calls = real, []

    def __call__(self, *args, **kw):
        self.calls.append((args, kw))
        return self.real(*args, **kw)


class TestSolve(unittest.TestCase):
    """The static spin-orbital bond solve and its outputs (gap.npz and
    eigenvalue.dat), one function call on a prepared run per case."""

    def test_solve_function_level(self):
        for layout, fx in _SOLVE.items():
            with self.subTest(layout=layout), _quiet_logs():
                p, value, m, text, path = _solve_case(self, fx)
                self.assertTrue(np.isfinite(value))
                self.assertEqual(str(m["gap_status"]), "certified")
                self.assertEqual(set(m), _GAP_MEMBERS[("eigenvalue", "certified")])
                self.assertEqual(str(m["block_layout"]), layout)
                self.assertEqual(m["pp_sector"].dtype, np.int8)
                self.assertEqual(m["pp_sector"].ndim, 0)
                self.assertEqual(m["pp_sectors"].dtype, np.int8)
                self.assertEqual(m["pp_sectors"].shape, m["eigenvalues"].shape)
                if layout == "dense":
                    self.assertEqual(int(m["pp_sector"]), 2)
                else:
                    self.assertIn(int(m["pp_sector"]), (1, 0, -1))
                self.assertTrue(bool(m["eigen_converged"]))
                self.assertEqual(float(np.real(m["eigenvalues"][0])), value)
                with np.load(path, allow_pickle=False) as z:
                    self.assertEqual(fingerprint.validate_structure(z, path, kind="gap",
                                                                    archive=path),
                                     ("gap", "sc_so_bond_static", 2))
                keys = _dat_keys(text)
                for key in _BOND_DAT_KEYS:
                    self.assertIn(key, keys)
                self.assertNotIn("iteration_converged", keys)
                topo = p.so_run.bond.topo
                self.assertEqual(keys["bond_channels"], "true")
                self.assertEqual(keys["bond_n_channels"], str(len(topo.delta_r)))
                self.assertEqual(keys["bond_delta_r"],
                                 str([[int(x) for x in dr] for dr in topo.delta_r]))
                self.assertEqual(keys["bond_max_shells"], "unset")
                self.assertEqual(keys["so_bond_layout"], layout)
                self.assertTrue(keys["so_bond_layout_reason"].startswith(
                    "so bond layout: %s (" % layout))
                self.assertEqual(keys["g2_tail"], "true")
                self.assertEqual(keys["bond_coeff_tail"], "1.0")
                self.assertEqual(keys["bond_coeff_tail_source"], "default")
                self.assertEqual(keys["approximation"], "static RPA-ladder bond dressing")
                self.assertEqual(keys["gap_status"], "certified")
                self.assertEqual(keys["pp_sector"], str(int(m["pp_sector"])))
                self.assertEqual(keys["producer_stage"], "sc_so_bond_static")
                self.assertEqual(keys["producer_label"], sc._SO_BOND_STATIC_LABEL)
                self.assertNotIn("lambda_pp", keys)
                self.assertEqual(sc._SO_BOND_STATIC_LABEL,
                                 "hwave_sc (spin-orbital, bond, static)")

    def test_zero_gap_weights_nan_silently(self):
        p = _bridge_prepared(F.ZEEMAN_EXTERN)
        so_green = _green_of(p)
        shape = (p.fixture.nd, p.fixture.nd) + tuple(p.shape)
        leg = sc.IterationLeg(float("nan"), None, 2, False, 0, (), "iteration")
        report = sc._so_bond_report("iteration", None, leg, gap_shape=shape)
        self.assertEqual(report.gap_status, "none")
        state = sc._so_own_state(so_green.eigenvalues, p.beta, so_green.mu, p.nmat,
                                 p.so_run.state.target)
        members = sc._so_writer_members("sc_so_bond_static", sc._SO_BOND_STATIC_LABEL,
                                        p.so_run, so_green.fp_green, state,
                                        bond=p.so_run.bond.inputs)
        d = tempfile.mkdtemp(prefix="sc_so_bond_zero_")
        self.addCleanup(shutil.rmtree, d, True)
        import warnings
        with warnings.catch_warnings(), _quiet_logs():
            warnings.simplefilter("error")
            sc._so_output_record(d, "gap.dat", p.norb_phys, None, None, None, None,
                                 p.kx, p.ky, p.kz, sz_residual=0.0,
                                 provenance_members=members, bond_report=report,
                                 block_layout="sz_packed")
        with np.load(os.path.join(d, "gap.npz"), allow_pickle=False) as z:
            m = {k: z[k] for k in z.files}
        self.assertEqual(set(m), _GAP_MEMBERS[("iteration", "none")])
        self.assertEqual(float(np.abs(m["gap"]).max()), 0.0)
        self.assertTrue(np.all(np.isnan(m["gap_spin_weights"])))
        for c in ("upup", "dndn", "singlet0", "triplet0"):
            self.assertEqual(m["sz_weight_" + c].shape, (1,))
            self.assertTrue(np.isnan(m["sz_weight_" + c][0]), c)
        self.assertTrue(np.isfinite(m["sz_conservation_residual"]))
        self.assertEqual(int(m["pp_sector"]), 2)
        self.assertEqual(str(m["gap_status"]), "none")

    def test_solver_mode_legs(self):
        eigen_solver = {"dense": "_so_solve_eigenvalue", "sz_packed": "_so_solve_sectors"}
        iter_solver = {"dense": "_so_solve_iteration", "sz_packed": "_so_iterate_sectors"}
        for layout, fx in _SOLVE.items():
            for mode in ("eigenvalue", "iteration", "both"):
                with self.subTest(layout=layout, mode=mode), _quiet_logs():
                    spies = {name: _Calls(getattr(sc, name)) for name in (
                        "_so_solve_eigenvalue", "_so_solve_iteration", "_so_solve_sectors",
                        "_so_iterate_sectors", "_so_bond_report")}
                    init = _Calls(so_pairing.so_initial_gap)
                    with contextlib.ExitStack() as stack:
                        for name, spy in spies.items():
                            stack.enter_context(mock.patch.object(sc, name, spy))
                        stack.enter_context(mock.patch.object(so_pairing, "so_initial_gap",
                                                              init))
                        # a short iteration: the legs and the member policy are
                        # checked here, not the convergence
                        _, value, m, text, _ = _solve_case(
                            self, fx, eliashberg={"Solver_Mode": mode, "Max_Iter": 20})
                    n = {name: len(spy.calls) for name, spy in spies.items()}
                    run_eigen = mode in ("eigenvalue", "both")
                    run_iter = mode in ("iteration", "both")
                    other = "sz_packed" if layout == "dense" else "dense"
                    self.assertEqual(n[eigen_solver[layout]], int(run_eigen))
                    self.assertEqual(n[iter_solver[layout]], int(run_iter))
                    self.assertEqual(n[eigen_solver[other]], 0)
                    self.assertEqual(n[iter_solver[other]], 0)
                    self.assertEqual(len(init.calls), int(run_iter))
                    self.assertEqual(n["_so_bond_report"], 1)
                    (args, kw), = spies["_so_bond_report"].calls
                    self.assertEqual(args[0], mode)
                    eigen, iteration = args[1], args[2]
                    if run_eigen:
                        self.assertIsInstance(eigen, sc.EigenLeg)
                    else:
                        self.assertIsNone(eigen)
                    if run_iter:
                        self.assertIsInstance(iteration, sc.IterationLeg)
                    else:
                        self.assertIsNone(iteration)
                    status = str(m["gap_status"])
                    # the eigen leg converges; a short iteration may not
                    self.assertEqual(status == "certified", bool(np.isfinite(value)))
                    if run_eigen:
                        self.assertEqual(status, "certified")
                    else:
                        self.assertIn(status, ("certified", "uncertified"))
                    self.assertEqual(set(m), _GAP_MEMBERS[(mode, status)])
                    keys = _dat_keys(text)
                    self.assertEqual("iteration_converged" in keys, mode == "both")
                    # one line, the one after "# Iteration eigenvalue"
                    self.assertEqual(text.count("# iteration_converged = "),
                                     int(mode == "both"))
                    if mode == "both":
                        lines = text.splitlines()
                        at = lines.index("# Iteration eigenvalue")
                        self.assertTrue(lines[at + 1].startswith("# iteration_converged = "))


@contextlib.contextmanager
def _quiet_logs():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


# ------------------------------------------------------------ the public entry
def _nn_rows(v):
    """A nearest-neighbour row ``v`` on orbital 0, on all four bonds."""
    return (((1, 0, 0), (0, 0), v), ((-1, 0, 0), (0, 0), v),
            ((0, 1, 0), (0, 0), v), ((0, -1, 0), (0, 0), v))


#: SOC_ONSITE with every off-site type declared on the nearest-neighbour
#: bonds at exactly zero: B = 5, the same Hamiltonian as SOC_ONSITE
_SOC_DECLARED_ZERO = F.SOC_ONSITE._replace(
    name="SOC_DECLARED_ZERO",
    interactions=F.SOC_ONSITE.interactions + tuple(
        (t, _nn_rows(0.0)) for t in ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop",
                                     "PairLift")))
#: an equal-spin interaction only (Hund on the nearest-neighbour bonds) on
#: the Zeeman band: the pair-S_z 0 sector of the kernel is exactly zero
_ZEEMAN_HUND = F.ZEEMAN_EXTERN._replace(name="ZEEMAN_HUND",
                                        interactions=(("Hund", _nn_rows(0.53)),))
_PACKED_REASON = ("so bond layout: sz_packed (S_z conserved: no spin-flip transfer, no "
                  "S_z-changing interaction)")
_CAP_LINES = ("spin-orbital bond memory cap: ", "no memory limit could be determined")


def _messages(records):
    return [r.getMessage() for r in records]


class TestEntry(unittest.TestCase):
    """The static spin-orbital bond path through the public entry
    (``hwave.sc.calc_eliashberg``): accepted configurations, the log lines,
    degenerate interactions and what the path never touches."""

    def test_b1_runs(self):
        run = F.run_case(F.SOC_ONSITE)
        msgs = _messages(run.records)
        self.assertEqual(msgs.count(
            "bond_channels = true: no off-site row is declared, so B = 1 and the bond "
            "vertex is the on-site spin-orbital vertex"), 1)
        self.assertIn("Bond channels: B = 1, Delta r = [[0, 0, 0]]", msgs)
        self.assertTrue(np.isfinite(run.value))
        self.assertEqual(_dat_keys(run.eigenvalue_dat)["bond_n_channels"], "1")

    def test_all_declared_zero_runs(self):
        # the operator: every off-site block of V_fluct is exactly zero
        d = tempfile.mkdtemp(prefix="sc_so_bond_zero_rows_")
        self.addCleanup(shutil.rmtree, d, True)
        p = F.prepare(_SOC_DECLARED_ZERO, d)
        self.assertEqual(len(p.so_run.bond.topo.delta_r), 5)
        parts = _PartsSpy()
        with mock.patch.object(so_bond, "so_bond_pairing_vertex_parts", parts), _quiet_logs():
            so_green = _green_of(p)
            adm = _admit(p)
            _operator(p, so_green, adm)
        self.assertEqual(adm.layout, "dense")
        v_fluct = _dense(parts.calls[0]["copies"][1])
        npair = p.fixture.nd ** 2
        geom = adm.geom
        onsite = [m for m, dr in enumerate(geom.delta_r) if not np.any(dr)]
        self.assertEqual(len(onsite), 1)
        keep = np.zeros(geom.B * npair, dtype=bool)
        keep[onsite[0] * npair:(onsite[0] + 1) * npair] = True
        off = ~(keep[:, None] & keep[None, :])
        self.assertGreater(np.abs(v_fluct[:, ~off]).max(), 1e-3)
        self.assertEqual(np.abs(v_fluct[:, off]).max(), 0.0)
        # the public entry: the leading eigenvalue of the B = 1 run
        zero, b1 = F.run_case(_SOC_DECLARED_ZERO), F.run_case(F.SOC_ONSITE)
        self.assertEqual(_dat_keys(zero.eigenvalue_dat)["bond_n_channels"], "5")
        self.assertTrue(np.isfinite(b1.value))
        self.assertLessEqual(abs(zero.value - b1.value), 1e-12)

    def test_null_sector_is_not_a_failure(self):
        run = F.run_case(_ZEEMAN_HUND)
        keys = _dat_keys(run.eigenvalue_dat)
        self.assertEqual(keys["so_bond_layout"], "sz_packed")
        self.assertIn("0:null-operator", keys["selection"].split(";"))
        self.assertTrue(bool(run.gap_npz["eigen_converged"]))
        self.assertEqual(str(run.gap_npz["gap_status"]), "certified")
        self.assertTrue(np.isfinite(run.value))

    def test_log_order(self):
        for fx, pairlift in ((F.SOC_ALL_TYPES, True), (F.ZEEMAN_EXTERN, False)):
            with self.subTest(fixture=fx.name):
                run = F.run_case(fx)
                msgs = _messages(run.records)
                topo_line = [i for i, m in enumerate(msgs) if m.startswith("Bond channels: B = ")]
                pl = [i for i, m in enumerate(msgs) if m == sc._SO_BOND_PAIRLIFT_INFO]
                b1 = [i for i, m in enumerate(msgs) if m == sc._SO_BOND_B1_INFO]
                ladder = [i for i, m in enumerate(msgs) if m == sc._SO_BOND_RPA_LADDER]
                cap = [i for i, m in enumerate(msgs) if m.startswith(_CAP_LINES)]
                reason = [i for i, m in enumerate(msgs) if m.startswith("so bond layout: ")]
                for found in (topo_line, ladder, cap, reason):
                    self.assertEqual(len(found), 1, msgs)
                self.assertEqual(b1, [])
                self.assertEqual(len(pl), int(pairlift))
                order = topo_line + pl + ladder + cap + reason
                self.assertEqual(order, sorted(order))
                self.assertEqual(order, list(range(order[0], order[0] + len(order))))
                self.assertEqual(run.records[ladder[0]].levelno, logging.WARNING)
                for i in topo_line + pl + reason:
                    self.assertEqual(run.records[i].levelno, logging.INFO)
        self.assertEqual(
            sc._SO_BOND_RPA_LADDER,
            "[eliashberg] bond_channels=true in spin-orbital mode is a STATIC RPA-ladder bond "
            "dressing on the BARE Green function built from the Transfer (and Extern) input. "
            "It demonstrates the qualitative V-dependence of the pairing eigenvalue; it is NOT "
            "a conserving FLEX result (the bond chi-bar is built here at RPA-ladder level), so "
            "absolute lambda values are not comparable with FLEX/dynamic references.")
        self.assertEqual(
            sc._SO_BOND_PAIRLIFT_INFO,
            "off-site PairLift is included in the spin-orbital bond vertex (the "
            "non-spin-orbital bond paths omit it: it has no spin/charge content)")

    def test_layout_reason_lines(self):
        flip = F.ZEEMAN_EXTERN._replace(
            name="ZEEMAN_FLIP_1E-14", transfer=(((0, 0, 0), (0, 1), 1e-14),
                                                ((0, 0, 0), (1, 0), 1e-14)))
        pairlift = F.ZEEMAN_EXTERN._replace(
            name="ZEEMAN_PAIRLIFT",
            interactions=F.ZEEMAN_EXTERN.interactions + (("PairLift", _nn_rows(0.07)),))
        rows = (
            (F.ZEEMAN_EXTERN, None, "sz_packed", _PACKED_REASON),
            (F.ZEEMAN_EXTERN, {"SO_Bond_Layout": "dense"}, "dense",
             'so bond layout: dense (reason: forced by so_bond_layout = "dense")'),
            (flip, None, "dense", "so bond layout: dense (reason: spin-flip transfer entry "
                                  "R=(0, 0, 0) (1, 2) = (1e-14+0j))"),
            (pairlift, None, "dense", "so bond layout: dense (reason: PairLift declared)"),
        )
        log = logging.getLogger("test_sc_spin_orbital_bond.reason")
        log.propagate = False
        log.setLevel(logging.INFO)
        for fx, eli, layout, reason in rows:
            with self.subTest(fixture=fx.name, eliashberg=eli):
                d = tempfile.mkdtemp(prefix="sc_so_bond_reason_")
                self.addCleanup(shutil.rmtree, d, True)
                p = F.prepare(fx, d, eliashberg=eli)
                h = _Infos()
                log.addHandler(h)
                try:
                    adm = sc._admit_so_bond_static(
                        p.so_run, cfg=p.cfg, cell=p.cell, nmat=p.nmat, nd=fx.nd,
                        shape=p.shape, solver_mode=p.solver_mode,
                        num_eigenvalues=p.num_eigenvalues, max_iter=p.max_iter,
                        spectral_shift=p.spectral_shift, eigenvalue_method="arnoldi", log=log)
                finally:
                    log.removeHandler(h)
                self.assertEqual((adm.layout, adm.layout_reason), (layout, reason))
                self.assertEqual(h.seen.count(reason), 1)

    def test_entry_inventory(self):
        legacy = ("_validate_bond_prereqs", "_bond_memory_estimate", "_bond_resource_preflight",
                  "_resolve_bond_topology", "_build_bond_vertices", "_bond_kernel_parts",
                  "_build_bond_operator", "_validate_bond_green_nmat",
                  "_bond_diagnostics_record", "_check_nonbond_offsite",
                  "_make_kernel_operator", "_so_static_vertex", "_calc_eigenvalues")
        for fx, layout in ((F.SOC_ALL_TYPES, "dense"), (F.ZEEMAN_EXTERN, "sz_packed")):
            with self.subTest(fixture=fx.name):
                with contextlib.ExitStack() as stack:
                    stack.enter_context(mock.patch.object(
                        sc, "_read_bond_config",
                        side_effect=AssertionError("the non-spin-orbital bond reader ran")))
                    spies = {name: stack.enter_context(mock.patch.object(
                        sc, name, wraps=getattr(sc, name))) for name in legacy}
                    get_backend = stack.enter_context(mock.patch.object(
                        backend, "get_backend", wraps=backend.get_backend))
                    topology = stack.enter_context(mock.patch.object(
                        so_bond, "so_bond_topology", wraps=so_bond.so_bond_topology))
                    geometry = stack.enter_context(mock.patch.object(
                        so_bond, "slot_geometry", wraps=so_bond.slot_geometry))
                    vertex = stack.enter_context(mock.patch.object(
                        so_bond, "so_bond_vertex", wraps=so_bond.so_bond_vertex))
                    run_inputs = _Returns(sc._so_run_inputs)
                    stack.enter_context(mock.patch.object(sc, "_so_run_inputs", run_inputs))
                    writer = stack.enter_context(mock.patch.object(
                        sc, "_so_writer_members", wraps=sc._so_writer_members))
                    run = F.run_case(fx, fresh=True, eliashberg={"Gpu": True})
                self.assertTrue(np.isfinite(run.value))
                self.assertEqual(str(run.gap_npz["block_layout"]), layout)
                for name, spy in spies.items():
                    self.assertEqual(spy.call_count, 0, name)
                self.assertEqual(get_backend.call_count, 0)
                self.assertEqual(topology.call_count, 1)
                so_run, = run_inputs.results
                topo, compiled = so_run.bond.topo, so_run.context.compiled
                # the admission and the vertex take the run's topology itself;
                # the bubble re-derives the geometry from the one built from it
                firsts = [call.args[0] for call in geometry.call_args_list]
                self.assertEqual(sum(a is topo for a in firsts), 2)
                for a in firsts:
                    self.assertTrue(a is topo or isinstance(a, so_bond.BondSlotGeometry), a)
                self.assertEqual(vertex.call_count, 1)
                self.assertIs(vertex.call_args.args[0], compiled)
                self.assertIs(vertex.call_args.args[1], topo)
                self.assertEqual(writer.call_count, 1)
                self.assertIs(writer.call_args.kwargs["bond"], so_run.bond.inputs)
                self.assertEqual(_messages(run.records).count(_SO_BOND_GPU_IGNORED), 1)

    def test_gpu_required_refused(self):
        d = tempfile.mkdtemp(prefix="sc_so_bond_gpu_")
        self.addCleanup(shutil.rmtree, d, True)
        inp = F.write_case(d, F.SOC_ALL_TYPES, eliashberg={"GPU_Required": True})
        with mock.patch.object(sc, "_read_interaction_files",
                               side_effect=AssertionError("a file was read")), \
                mock.patch.object(backend, "get_backend",
                                  side_effect=AssertionError("a backend was selected")):
            with self.assertRaises(ValueError) as cm:
                sc.calc_eliashberg(inp)
        self.assertEqual(str(cm.exception), _SO_BOND_GPU_REQUIRED_REFUSAL)
        self.assertFalse(os.path.exists(os.path.join(d, "out")))


class _Returns:
    """A spy that records what each call of the real function returns."""

    def __init__(self, real):
        self.real, self.results = real, []

    def __call__(self, *args, **kw):
        out = self.real(*args, **kw)
        self.results.append(out)
        return out


class _Infos(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.INFO)
        self.seen = []

    def emit(self, record):
        self.seen.append(record.getMessage())


def _config(eli, **mode):
    return {"mode": dict({"Enable_Spin_Orbital": True}, **mode), "eliashberg": eli}


class TestFlip(unittest.TestCase):
    """What the configuration and interaction guards accept with
    ``bond_channels = true`` (function level)."""

    def _refusal(self, fn, *args):
        with self.assertRaises(ValueError) as cm:
            fn(*args)
        return str(cm.exception)

    def test_static_bond_accepted_dynamic_refused(self):
        bond = {"Bond_Channels": True}
        so_pairing.check_so_config(_config(bond))
        so_pairing.check_so_config(_config(dict(bond, Frequency="static")))
        self.assertEqual(self._refusal(so_pairing.check_so_config,
                                       _config(dict(bond, Frequency="dynamic"))),
                         so_pairing._SO_DYNAMIC_BOND_REFUSAL)
        # not used on the static bond path (the entry warns), refused on-site
        flex = ("[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
                "flex_chi_s names a spin or charge channel file, which a spin-orbital FLEX "
                "run does not write (it stores the single spinful channel as chiq). Remove "
                "flex_chi_s and set flex_chiq to the chiq archive of the spin-orbital FLEX "
                "run (default \"chiq.npz\").")
        tensor = ("[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
                  "chi0q_tensor = 'reduced' is not supported; the spin-orbital vertex needs "
                  "the general (4-index) chi0q. Set chi0q_tensor = \"general\" (or \"auto\").")
        for extra, text in (({"Flex_Chi_S": "chis.npz"}, flex),
                            ({"Chi0q_Tensor": "reduced"}, tensor)):
            with self.subTest(extra=extra):
                so_pairing.check_so_config(_config(dict(bond, **extra)))
                self.assertEqual(self._refusal(so_pairing.check_so_config,
                                               _config(dict(extra))), text)
                self.assertEqual(self._refusal(so_pairing.check_so_config,
                                               _config(dict(extra, Bond_Channels=False))),
                                 text)

    def test_chi0q_mode_value_check_kept(self):
        for value in ("FLEX", "Flex"):
            with self.subTest(value=value):
                eli = {"Bond_Channels": True, "Chi0q_Mode": value}
                self.assertEqual(
                    self._refusal(so_pairing.check_so_config, _config(eli)),
                    "[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
                    "chi0q_mode = {!r} is not supported; use \"calc\" (the spinful chi0 is "
                    "computed here), \"load\" (a general chi0q of a spin-orbital RPA run) or "
                    "\"flex\" (the chiq and green of a spin-orbital FLEX run).".format(value))
        eli = {"Bond_Channels": True, "Chi0q_Mode": "flex"}
        so_pairing.check_so_config(_config(eli))
        self.assertEqual(self._refusal(sc._read_so_bond_config, CaseInsensitiveDict(eli)),
                         _SO_BOND_FLEX_REFUSAL)

    def test_offsite_accepted_with_bond_channels(self):
        from tests.test_sc_spin_orbital import OFFSITE_TEXT
        rows = {"CoulombIntra": {((0, 0, 0), (0, 0)): 0.83},
                "CoulombInter": {((1, 0, 0), (0, 0)): 0.37, ((-1, 0, 0), (0, 0)): 0.37}}

        def inp(eli):
            d = _config(eli)
            d["file"] = {"input": {"interaction": {"CoulombIntra": "coulombintra.dat",
                                                   "CoulombInter": "coulombinter.dat"}}}
            return d
        so_pairing.check_so_supported(inp({"Bond_Channels": True}), rows, norb_phys=1)
        so_pairing.check_so_supported(inp({"BOND_CHANNELS": "true", "FREQUENCY": "static"}),
                                      rows, norb_phys=1)
        static_text = OFFSITE_TEXT.format(R=(1, 0, 0), orbs=(1, 1), src="the CoulombInter file")
        for eli in ({}, {"Bond_Channels": False}):
            with self.subTest(eli=eli):
                self.assertEqual(self._refusal(so_pairing.check_so_supported, inp(eli), rows,
                                               None, 1), static_text)
        self.assertEqual(self._refusal(so_pairing.check_so_supported,
                                       inp({"Bond_Channels": True, "Frequency": "dynamic"}),
                                       rows, None, 1),
                         so_pairing._SO_DYNAMIC_BOND_REFUSAL)
        # the orbital-range check still runs on the bond path
        bad = dict(rows, CoulombInter={((1, 0, 0), (0, 1)): 0.37})
        text = self._refusal(so_pairing.check_so_supported, inp({"Bond_Channels": True}), bad,
                             None, 1)
        self.assertIn("outside the physical orbitals 1..1", text)


# ------------------------------------------------------------ memory admission
_ADM_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond.admission")
#: the admission fixtures: 4 x 4, B = 5, nd = 2, Nmat 32, k = 4 (the dense
#: antisymmetric dimension 28 > k + 1: the ARPACK branch)
_ADM_ELI = {"Num_Eigenvalues": 4}
_ADM = {"dense": F.SOC_ALL_TYPES, "sz_packed": F.ZEEMAN_EXTERN}
#: the phases whose rows bound a traced peak (``vertex`` and ``chi_static``
#: have no traced call of their own: ``vertex`` is the resident Gamma_B
#: alone, the ``chi_static`` rows are a subset of the ``pairing_static``
#: rows and run inside the vertex phase)
TRACED_PHASES = ("green_static", "pairing_static", "g2_static", "kernel_static")
_CAP_ZERO_WARNING = ("[eliashberg] bond_memory_cap_gb = 0: the memory admission check of the "
                     "spin-orbital bond path is disabled")
_NO_LIMIT_WARNING = ("no memory limit could be determined; the admission check is disabled; "
                     "set [eliashberg] bond_memory_cap_gb to enable it")


def plan_phases(plan):
    """The phases of ``plan``: every row phase but ``resident``, and the
    closed consumers."""
    return {ph for ph, _, _ in plan.rows} - {"resident"} | set(plan.consumers)


def _solve_leg(p, static, mode):
    """The solve of ``mode`` (``"eigenvalue"`` or ``"iteration"``) on the
    operator ``static``, as the production solve calls it (the raw
    ``spectral_shift`` resolved per leg). Returns the selection of a dense
    eigen solve (else ``None``). This repeats the legs of
    ``sc._solve_so_bond_static`` (its constants and init-gap steps): a
    change there is not seen here, and ``TestAdmission.test_end_to_end_smoke``,
    which traces the production run, is the guard against the drift."""
    nd = 2 * p.norb_phys
    gap_shape = (nd, nd) + tuple(p.shape)
    shift = sc._resolve_spectral_shift(p.spectral_shift, mode)
    if mode == "iteration":
        sigma = so_pairing.so_initial_gap(sc._resolve_init_gap(p.eli_param.get("init_gap"),
                                                               "singlet"),
                                          p.norb_phys, p.kx, p.ky, p.kz)
        sigma = sc._project_gap_parity(sigma, "triplet")
        sigma = sigma / np.linalg.norm(sigma)
    if static.layout == "dense":
        if mode == "iteration":
            sc._so_solve_iteration(static.operator, static.basis, sigma, p.max_iter, 0.5,
                                   1.0e-5, spectral_shift=shift, info_out={})
            return None
        info = {}
        sc._so_solve_eigenvalue(static.operator, static.basis, p.num_eigenvalues, "arnoldi",
                                sigma_shift=None, spectral_shift=shift, info_out=info)
        return info["selection"]
    ops = tuple(static.kernel.sector_op(q) for q in (1, 0, -1))
    if mode == "iteration":
        sc._so_iterate_sectors(ops, sigma, gap_shape=gap_shape, max_iter=p.max_iter,
                               alpha=0.5, tol=1.0e-5, spectral_shift=shift)
    else:
        sc._so_solve_sectors(ops, p.num_eigenvalues, "arnoldi", gap_shape=gap_shape,
                             spectral_shift=shift)
    return None


def traced_phases(p, mode):
    """``(admission, measured, selection)`` of the static chain of the
    prepared run ``p`` traced phase by phase.

    The admission runs untraced; then ONE tracemalloc session covers the
    phase functions in production order (the Green function, the vertex,
    the pair propagator, the kernel and the solve of ``mode``), with
    ``reset_peak`` before each. ``measured[P]`` is the traced peak during
    phase ``P`` of everything allocated since the session started (the
    arrays of earlier phases still alive included); the solve is
    ``measured["kernel_static/solve"]``, bounded by the ``kernel_static``
    rows as well. Gamma_B is released after the vertex phase, as in
    production. Tracing counts only what is allocated after the session
    starts: the arrays of ``p`` (allocated before) are not in ``measured``,
    although the plan's rows include them, so the check errs on the weak
    side."""
    if tracemalloc.is_tracing():
        raise unittest.SkipTest("tracemalloc is already tracing")
    adm = _admit(p)
    measured = {}
    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        so_green = _green_of(p)
        measured["green_static"] = tracemalloc.get_traced_memory()[1]
        tracemalloc.reset_peak()
        gamma_b, V = sc._so_bond_static_vertex(p.so_run, adm, so_green, beta=p.beta,
                                               nmat=p.nmat)
        measured["pairing_static"] = tracemalloc.get_traced_memory()[1]
        del gamma_b
        tracemalloc.reset_peak()
        G2 = sc._so_bond_static_g2(so_green, beta=p.beta, g2_tail=p.cfg.g2_tail)
        measured["g2_static"] = tracemalloc.get_traced_memory()[1]
        tracemalloc.reset_peak()
        static = sc._so_bond_static_kernel(p.so_run, adm, V, G2, so_green, beta=p.beta,
                                           norb_phys=p.norb_phys)
        measured["kernel_static"] = tracemalloc.get_traced_memory()[1]
        tracemalloc.reset_peak()
        selection = _solve_leg(p, static, mode)
        measured["kernel_static/solve"] = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()
    return adm, measured, selection


def check_phase_rows(test, p, mode):
    """Every traced phase of ``p`` within its ``phase_total`` (no margin);
    the ratios and the peak phase are printed on failure. Returns the
    admission."""
    adm, measured, selection = traced_phases(p, mode)
    plan = adm.plan
    peak = max(plan_phases(plan), key=plan.phase_total)
    rows = {key: plan.phase_total(key.split("/")[0]) for key in measured}
    ratios = ", ".join("%s %.3f (%d / %d)" % (key, measured[key] / rows[key], measured[key],
                                              rows[key]) for key in measured)
    msg = "%s %s Nmat %d: peak phase %s; traced / rows: %s" % (
        adm.layout, mode, p.nmat, peak, ratios)
    _ADM_LOG.info(msg)
    test.assertEqual(peak, plan.peak_phase, msg)
    test.assertEqual(set(TRACED_PHASES) | {"vertex", "chi_static"}, plan_phases(plan), msg)
    for key in measured:
        with test.subTest(phase=key):
            test.assertLessEqual(measured[key], rows[key], msg)
    if adm.layout == "dense" and mode == "eigenvalue":
        # the ARPACK branch: the dense fallback is reserved, not run
        test.assertNotEqual(selection, "dense-LR", msg)
    return adm


def _prepare_admission(fixture, eliashberg=None, param=None):
    """:func:`tests.so_bond_entry_fixtures.prepare` of an admission fixture
    (its own [eliashberg] values: no bridge overrides), in a directory
    removed on return."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_admission_")
    try:
        return F.prepare(fixture, d, eliashberg=dict(_ADM_ELI, **(eliashberg or {})),
                         param=param)
    finally:
        shutil.rmtree(d, True)


class TestAdmission(unittest.TestCase):
    """Gate P-mem of the static spin-orbital bond path: the memory
    admission decides before the carrier, and every phase's traced peak
    stays within that phase's rows.

    Traced rows (``tracemalloc`` sees numpy and scipy array buffers, f2py
    work arrays included), which bound the measured peak of their phase:
    ``green_static`` -- ``carrier_build``, ``call_overhead``,
    ``eigendecomposition`` and ``transient``; ``pairing_static`` -- the
    consumer rows (``chibar_static``, ``bubble_static_temp``,
    ``bubble_pair_block``, ``chi_static``, ``dress_batch``, ``v_bare``,
    ``v_fluct``, ``pairing_q_workspace``) and ``carrier``; ``g2_static`` --
    ``carrier``, ``v_static``, ``g2_roll``, ``g2_out``, ``call_overhead``;
    ``kernel_static`` -- ``carrier``, ``v_static``, ``g2``,
    ``kernel_workspace`` and the Python share of ``eigensolver`` (the ARPACK
    arrays, the shift's retained outputs, the gap work).

    Derived rows (not seen by ``tracemalloc``; covered by their arithmetic
    and the forced-failure tests of the bond plan): the native part of the
    bubble row (the FFT scratch), and in ``eigensolver`` the native share of
    the dense fallback at the ``eig`` call (its input copy, the
    eigenvectors, the eigenvalues and the LAPACK ``work`` / ``rwork``
    arrays) and the second Krylov working set (reserved only for an
    antisymmetric dimension above 2000; none here). LAPACK workspace that
    numpy allocates with ``malloc`` inside ``eig`` / ``eigh`` / ``solve``
    and the BLAS buffers are not traced either. A traced peak below its
    rows therefore proves the traced part only; the fallback and the
    second working set are not reached in a passing solve.

    Gamma_B is released after the vertex phase but stays a resident row
    (conservative). The peak phase is the argmax of ``phase_total`` over
    the plan's phases, never assumed. The row set and the row values are
    pinned literally in ``tests/test_so_bond_plan_static.py``; this class
    checks only that the traced allocations stay within them (a ratio
    cannot catch a dropped row)."""

    def test_phase_rows_dense(self):
        check_phase_rows(self, _prepare_admission(_ADM["dense"]), "eigenvalue")

    def test_phase_rows_packed(self):
        adm = check_phase_rows(self, _prepare_admission(_ADM["sz_packed"]), "eigenvalue")
        self.assertEqual(adm.layout, "sz_packed")

    def test_phase_rows_iteration_smoke(self):
        p = _prepare_admission(_ADM["sz_packed"], param={"Nmat": 8},
                               eliashberg={"Solver_Mode": "iteration", "Max_Iter": 1000})
        self.assertIsNone(p.spectral_shift)
        adm = check_phase_rows(self, p, "iteration")
        self.assertEqual((adm.layout, adm.plan.nmat), ("sz_packed", 8))

    def test_peak_phase_consistent(self):
        cases = (("dense", _ADM["dense"], None, None),
                 ("sz_packed", _ADM["sz_packed"], None, None),
                 ("dense", F.SOC_ONSITE, {"Nmat": 8}, "kernel_static"))
        for layout, fx, param, want in cases:
            with self.subTest(fixture=fx.name, param=param):
                plan = _admit(_prepare_admission(fx, param=param)).plan
                self.assertEqual(plan.layout, layout)
                phases = plan_phases(plan)
                self.assertEqual(phases, set(TRACED_PHASES) | {"vertex", "chi_static"})
                peak = max(phases, key=plan.phase_total)
                self.assertEqual(peak, plan.peak_phase)
                self.assertEqual(plan.peak_host, 1.25 * plan.phase_total(peak))
                if want is not None:
                    self.assertEqual(plan.B, 1)
                    self.assertEqual(peak, want)

    def test_end_to_end_smoke(self):
        """One traced session from the admission to the end of
        ``_save_results``, in both layouts: the traced peak of the whole run
        stays within the admitted peak (``peak_host / 1.25``); at the output
        phase the operator (the packed kernel or the dense operator) and the
        Green function carrier are freed (by reference counting alone: the
        garbage collector is off during the run); and the archive
        serializes the report's listing itself, not a second copy."""
        if tracemalloc.is_tracing():
            raise unittest.SkipTest("tracemalloc is already tracing")
        seen = {}
        real_admit, real_green = sc._admit_so_bond_static, sc._so_bond_green
        real_operator = sc._build_so_bond_operator
        real_record, real_save = sc._so_output_record, sc._save_results
        real_savez = np.savez

        def admit(*args, **kw):
            adm = real_admit(*args, **kw)
            seen["plan"] = adm.plan
            tracemalloc.start()
            return adm

        def green(*args, **kw):
            out = real_green(*args, **kw)
            seen["carrier"] = weakref.ref(out.carrier.full_sc)
            return out

        def operator(*args, **kw):
            out = real_operator(*args, **kw)
            seen["operator"] = weakref.ref(out.kernel if out.layout == "sz_packed"
                                           else out.operator)
            return out

        def record(*args, **kw):
            seen["alive"] = {k: seen[k]() is not None for k in ("carrier", "operator")}
            seen["report"] = kw["bond_report"]
            return real_record(*args, **kw)

        def savez(path, **members):
            seen.setdefault("archives", []).append(members)
            return real_savez(path, **members)

        def save(*args, **kw):
            out = real_save(*args, **kw)
            seen["peak"] = tracemalloc.get_traced_memory()[1]
            tracemalloc.stop()
            return out

        patches = ((sc, "_admit_so_bond_static", admit), (sc, "_so_bond_green", green),
                   (sc, "_build_so_bond_operator", operator),
                   (sc, "_so_output_record", record), (np, "savez", savez),
                   (sc, "_save_results", save))
        for layout, fixture in _ADM.items():
            seen.clear()
            with self.subTest(layout=layout):
                enabled = gc.isenabled()
                gc.disable()
                try:
                    with contextlib.ExitStack() as stack:
                        for target, name, fn in patches:
                            stack.enter_context(mock.patch.object(target, name, fn))
                        run = F.run_case(fixture, fresh=True, eliashberg=_ADM_ELI)
                finally:
                    if tracemalloc.is_tracing():
                        tracemalloc.stop()
                    if enabled:
                        gc.enable()
                self.assertTrue(np.isfinite(run.value))
                self.assertEqual(str(run.gap_npz["block_layout"]), layout)
                plan = seen["plan"]
                bound = plan.peak_host / 1.25
                self.assertLessEqual(seen["peak"], bound, "traced %d bytes, bound %d bytes "
                                     "(%.3f)" % (seen["peak"], bound, seen["peak"] / bound))
                self.assertEqual(seen["alive"], {"carrier": False, "operator": False})
                archive, = seen["archives"]
                listing = seen["report"].gaps
                self.assertGreater(len(listing), 0)
                self.assertTrue(np.shares_memory(archive["eigenvectors"], listing),
                                "the archive serializes a second copy of the eigenvector "
                                "listing")

    def test_over_cap_refused_before_carrier(self):
        green = _Calls(sc._build_bond_green)
        vertex = _Calls(so_bond.so_bond_vertex)
        with mock.patch.object(sc, "_build_bond_green", green), \
                mock.patch.object(so_bond, "so_bond_vertex", vertex), _quiet_logs():
            with self.assertRaises(ValueError) as cm:
                F.run_case(_ADM["dense"], fresh=True,
                           eliashberg=dict(_ADM_ELI, Bond_Memory_Cap_GB=1e-6))
        text = str(cm.exception)
        self.assertTrue(text.startswith("spin-orbital bond Eliashberg needs an estimated peak "),
                        text)
        for part in ("Nmat = 32", "Reduce Nmat or the k grid", "raise " + _CAP_KEY,
                     "(" + _CAP_KEY + ")"):
            self.assertIn(part, text)
        self.assertEqual(green.calls, [])
        self.assertEqual(vertex.calls, [])

    def test_kernel_admission_does_not_refuse(self):
        # dense: no memory limit can be determined (the cap is unlimited, one
        # WARNING). The contract is the call count: the kernel is built with
        # its own cap disabled, so its estimate never reads the host memory
        # query (patched to report one byte: a kernel that read it would
        # also refuse the run)
        host_query = mock.Mock(return_value=1)
        with mock.patch.object(eliashberg_bond, "_host_available_bytes", host_query), \
                mock.patch.object(host_memory, "available_host_bytes",
                                  return_value=host_memory.HostMemory(None, "none")):
            run = F.run_case(_ADM["dense"], fresh=True, eliashberg=_ADM_ELI)
        self.assertTrue(np.isfinite(run.value))
        self.assertEqual(host_query.call_count, 0)
        warnings = [r.getMessage() for r in run.records if r.levelno == logging.WARNING]
        self.assertEqual(warnings.count(_NO_LIMIT_WARNING), 1, warnings)
        self.assertEqual(_dat_keys(run.eigenvalue_dat)["bond_memory_cap"],
                         "unlimited (no memory limit could be determined)")
        # packed: a cap one byte above the admitted peak
        plan = _admit(_prepare_admission(_ADM["sz_packed"])).plan
        cap_gb = (plan.peak_host + 1.0) / 2.0 ** 30
        self.assertGreater(cap_gb * 2.0 ** 30, plan.peak_host)
        run = F.run_case(_ADM["sz_packed"], fresh=True,
                         eliashberg=dict(_ADM_ELI, Bond_Memory_Cap_GB=cap_gb))
        self.assertTrue(np.isfinite(run.value))
        self.assertEqual(str(run.gap_npz["block_layout"]), "sz_packed")

    def test_cap_zero_warning(self):
        run = F.run_case(_ADM["dense"], fresh=True, eliashberg=dict(_ADM_ELI,
                                                                     Bond_Memory_Cap_GB=0))
        self.assertTrue(np.isfinite(run.value))
        warnings = [r.getMessage() for r in run.records if r.levelno == logging.WARNING]
        self.assertEqual(warnings.count(_CAP_ZERO_WARNING), 1, warnings)
        self.assertIn("# bond_memory_cap = unlimited ([eliashberg] bond_memory_cap_gb = 0)",
                      run.eigenvalue_dat.splitlines())


# ------------------------------------------------------------ outputs and the reporting contract
#: the gap archive members of every static spin-orbital bond run, written
#: out literally: the provenance members of the stage and the gap-level
#: members
_LIT_BASE = frozenset({
    "delta_r", "fingerprint_version", "fp_band", "fp_capability", "fp_extern", "fp_green",
    "fp_interaction", "fp_meta", "fp_resolved_topology", "fp_topology", "fp_transfer",
    "max_shells", "producer_label", "producer_stage", "producer_version", "prov_cell_shape",
    "prov_coeff_extern", "prov_enable_spin_orbital", "prov_extern_keys", "prov_extern_present",
    "prov_extern_vals", "prov_index_convention", "prov_norb", "prov_norb_phys",
    "prov_rec_cells", "prov_rec_modes", "prov_rec_types", "prov_rec_v", "prov_sub_shape",
    "prov_topology", "prov_transfer_keys", "prov_transfer_vals", "reverse", "slot_convention",
    "state_beta", "state_calc_mu", "state_density", "state_mu", "state_mu_tol", "state_nmat",
    "state_target",
    "index_convention", "spin_orbital", "norb_phys", "kx", "ky", "kz", "spin_components",
    "gap", "gap_spin_weights", "eigen_converged", "pp_sector", "gap_status", "block_layout",
    "sz_weight_upup", "sz_weight_dndn", "sz_weight_singlet0", "sz_weight_triplet0",
    "sz_conservation_residual", "sz_conserved"})
#: ... and with a non-empty eigenvalue listing also these
_LIT_WITH_LISTING = _LIT_BASE | {"eigenvalues", "eigenvectors", "eigenvector_spin_weights",
                                 "eigenvalue_match", "pp_sectors"}
_SZ_COMPONENTS = ("upup", "dndn", "singlet0", "triplet0")
#: the ``# key = value`` keys of eigenvalue.dat of a static spin-orbital bond
#: run in solver_mode = "eigenvalue"
_DAT_KEY_SET = frozenset({
    "approximation", "bond_channels", "bond_coeff_tail", "bond_coeff_tail_source",
    "bond_delta_r", "bond_max_shells", "bond_memory_cap", "bond_n_channels",
    "eigen_converged", "fingerprint_version", "fp_band", "fp_capability", "fp_extern",
    "fp_green", "fp_interaction", "fp_meta", "fp_resolved_topology", "fp_topology",
    "fp_transfer", "g2_tail", "gap_layout", "gap_npz", "gap_spin_weights", "gap_status",
    "index_convention", "leading_eigenvalue", "match_column", "max_shells", "pp_sector",
    "producer_label", "producer_stage", "producer_version", "prov_coeff_extern",
    "prov_enable_spin_orbital", "prov_extern_present", "prov_index_convention", "prov_norb",
    "prov_norb_phys", "selection", "slot_convention", "so_bond_layout",
    "so_bond_layout_reason", "spin_decomposition", "spin_orbital", "state_beta",
    "state_calc_mu", "state_density", "state_mu", "state_mu_tol", "state_nmat",
    "state_target", "weight_columns"})
#: the cap of the runs whose bond_memory_cap line is compared exactly (an
#: unset cap names the momentary available memory)
_FIXED_CAP = {"Bond_Memory_Cap_GB": 2.0}
_UNSET_CAP = re.compile(r"^\d+\.\d{3} GiB \(0\.8 x available \d+\.\d{3} GiB, source: .+\)$"
                        r"|^unlimited \(no memory limit could be determined\)$")


def _dat_count(text, key):
    """How many ``# key = ...`` lines an eigenvalue file has."""
    return sum(1 for line in text.splitlines() if line.startswith("# %s = " % key))


def _exit_status(value):
    """The exit status of the ``hwave_sc`` command whose run returned
    ``value`` (``main`` with the run replaced by that value)."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_exit_")
    try:
        path = os.path.join(d, "in.toml")
        with open(path, "w") as f:
            f.write("[mode]\nmode = \"SC\"\n")
        with mock.patch.object(sc, "calc_eliashberg", return_value=value), \
                mock.patch("sys.argv", ["hwave_sc", "-q", path]), \
                mock.patch("logging.basicConfig"), _quiet_logs():
            try:
                sc.main()
            except SystemExit as exc:
                return exc.code
        return 0
    finally:
        shutil.rmtree(d, True)


def _policy_rows(shape):
    """The ten rows of the gap archive policy: ``(name, solver_mode, eigen
    leg, iteration leg, expected gap_status, expected members, expected gap,
    S_z weight length, finite weights)``; the legs are hand-made."""
    rng = np.random.default_rng(29)

    def gap():
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    g = [gap() for _ in range(6)]
    nan = float("nan")
    zeros = np.zeros(shape, dtype=complex)
    e_ok = sc.EigenLeg(np.array([0.5, 0.2], complex), np.array([g[0], g[1]]),
                       np.array([0, 1], np.int8), True, 0.5 + 0j, "+1:dense-LR;0:dense-LR", ())
    e_unconv = sc.EigenLeg(np.array([0.2], complex), np.array([g[1]]), np.array([1], np.int8),
                           False, complex(nan), "+1:dense-LR;0:not-converged", (0,))
    e_sentinel = sc.EigenLeg.from_dense(np.full(2, nan + 0j), np.array([g[2], g[3]]),
                                        selection="not-converged")
    e_empty = sc.EigenLeg(np.zeros(0, complex), np.zeros((0,) + shape, complex),
                          np.zeros(0, np.int8), False, complex(nan),
                          "+1:empty;0:empty;-1:empty", ())
    i_ok = sc.IterationLeg(0.45, g[4], 0, True, 12, (), "iteration")
    i_unconv = sc.IterationLeg(nan, g[5], -1, False, 50, (-1,), "iteration")
    i_none = sc.IterationLeg(nan, None, 2, False, 0, (), "iteration")
    L, B = _LIT_WITH_LISTING, _LIT_BASE
    return (
        ("eigenvalue converged", "eigenvalue", e_ok, None, "certified", L, g[0], 2, True),
        ("eigenvalue finite unconverged", "eigenvalue", e_unconv, None, "uncertified", L, g[1],
         1, True),
        ("eigenvalue dense NaN sentinel", "eigenvalue", e_sentinel, None, "none", L, zeros, 1,
         False),
        ("eigenvalue packed empty", "eigenvalue", e_empty, None, "none", B, zeros, 1, False),
        ("iteration converged", "iteration", None, i_ok, "certified", B, g[4], 1, True),
        ("iteration unconverged", "iteration", None, i_unconv, "uncertified", B, g[5], 1, True),
        ("iteration no gap", "iteration", None, i_none, "none", B, zeros, 1, False),
        ("both eigen converged", "both", e_ok, i_unconv, "certified", L, g[0], 2, True),
        ("both eigen failed, iteration gap", "both", e_sentinel, i_ok, "uncertified", L, g[4],
         1, True),
        ("both no candidate", "both", e_empty, i_none, "none", B, zeros, 1, False),
    )


#: the SO on-site static run (bond channels off)
_ONSITE_ELI = {"Bond_Channels": False, "Chi0q_Mode": "calc"}


def _null_sector_op(op):
    """``op`` with its kernel replaced by the exact zero."""
    from scipy.sparse.linalg import LinearOperator
    n = op.basis.n
    zero = LinearOperator((n, n), matvec=lambda x: np.zeros(n, complex),
                          rmatvec=lambda y: np.zeros(n, complex), dtype=complex)
    return sc.SectorOp(op.pp_sector, zero, op.basis, True)


def _empty_sector_op(kernel, op):
    """``op`` with an empty fermion-antisymmetric gap space."""
    from scipy.sparse.linalg import LinearOperator
    basis = sc._SOAntisymmetricBasis.from_subset(kernel.gap_shape, np.zeros(0, np.int64),
                                                 reverse=sc._reverse_k_and_orbital)
    zero = LinearOperator((0, 0), matvec=lambda x: np.zeros(0, complex),
                          rmatvec=lambda y: np.zeros(0, complex), dtype=complex)
    return sc.SectorOp(op.pp_sector, zero, basis, True)


@contextlib.contextmanager
def patched_sectors(null=(), empty=()):
    """The S_z-packed kernel's sector operators with the sectors ``null``
    replaced by the exact zero and the sectors ``empty`` by an empty gap
    space."""
    real = so_bond.PackedBondPairKernel.sector_op

    def sector_op(self, p):
        op = real(self, p)
        if p in null:
            return _null_sector_op(op)
        if p in empty:
            return _empty_sector_op(self, op)
        return op
    with mock.patch.object(so_bond.PackedBondPairKernel, "sector_op", sector_op):
        yield


#: the antisymmetric dimension of the Zeeman fixture's pair-S_z 0 sector
#: (the leader's; the +1 and -1 sectors have 6)
ZEEMAN_DIM_PP0 = 16
_REAL_DENSE_EIGEN = sc._so_dense_reduced_eigen


def eigen_fails(dim):
    """The dense sector diagonalization of the ``dim``-dimensional sector
    raises (the sector is failed, "not-converged")."""
    def fake(red, d, k, seed_vec=None):
        if d == dim:
            raise np.linalg.LinAlgError("forced")
        return _REAL_DENSE_EIGEN(red, d, k, seed_vec)
    return mock.patch.object(sc, "_so_dense_reduced_eigen", fake)


def empty_gap_warnings(run):
    """How many times ``run`` logged the empty gap space WARNING."""
    return sum(1 for r in run.records
               if r.levelno == logging.WARNING and r.getMessage() == sc._SO_BOND_EMPTY_GAP)


def reporting_cell(test, fixture, mode, *, layout="auto", patches=()):
    """``(run, keys, exit status)`` of a fresh public-entry run of
    ``fixture`` in ``mode`` with the context managers ``patches`` active."""
    with contextlib.ExitStack() as stack:
        for p in patches:
            stack.enter_context(p)
        run = F.run_case(fixture, fresh=True, layout=layout, eliashberg={"Solver_Mode": mode})
    return run, _dat_keys(run.eigenvalue_dat), _exit_status(run.value)


def check_reporting_cell(test, run, keys, status, want):
    """One cell of the reporting matrix: ``want`` = ``(reported, gap_status,
    pp_sector, eigen_converged, iteration_converged, exit status)``;
    ``reported`` is a number (compared to the eigenvalue tolerance of
    run-to-run compares) or NaN, ``iteration_converged`` None outside
    ``"both"``."""
    from tests import so_run_compare
    reported, gap_status, pp, eigen_ok, iter_ok, exit_status = want
    m = run.gap_npz
    if np.isnan(reported):
        test.assertTrue(np.isnan(run.value), run.value)
    else:
        test.assertTrue(so_run_compare.eigenvalues_close([run.value], [reported]),
                        (run.value, reported))
    test.assertEqual(str(m["gap_status"]), gap_status)
    test.assertEqual(keys["gap_status"], gap_status)
    test.assertEqual(int(m["pp_sector"]), pp)
    test.assertEqual(keys["pp_sector"], str(pp))
    test.assertEqual(bool(m["eigen_converged"]), eigen_ok)
    test.assertEqual(keys["eigen_converged"], "true" if eigen_ok else "false")
    if iter_ok is None:
        test.assertEqual(_dat_count(run.eigenvalue_dat, "iteration_converged"), 0)
    else:
        test.assertEqual(_dat_count(run.eigenvalue_dat, "iteration_converged"), 1)
        test.assertEqual(keys["iteration_converged"], "true" if iter_ok else "false")
    test.assertEqual(status, exit_status)
    if gap_status == "none":
        test.assertEqual(np.count_nonzero(m["gap"]), 0)
    else:
        test.assertGreater(np.count_nonzero(m["gap"]), 0)


class TestOutputs(unittest.TestCase):
    """The outputs of the static spin-orbital bond solve: the gap archive's
    member policy, its structural check, the eigenvalue.dat keys, the empty
    gap space, and the on-site path staying out of the bond code."""

    def test_gap_member_policy(self):
        p = _bridge_prepared(F.ZEEMAN_EXTERN)
        so_green = _green_of(p)
        shape = (p.fixture.nd, p.fixture.nd) + tuple(p.shape)
        state = sc._so_own_state(so_green.eigenvalues, p.beta, so_green.mu, p.nmat,
                                 p.so_run.state.target)
        members = sc._so_writer_members("sc_so_bond_static", sc._SO_BOND_STATIC_LABEL,
                                        p.so_run, so_green.fp_green, state,
                                        bond=p.so_run.bond.inputs)
        # the literal sets are the ones the run's own policy names
        self.assertEqual(_LIT_BASE, _GAP_BASE)
        self.assertEqual(_LIT_WITH_LISTING, _GAP_BASE | _GAP_LISTING)
        rows = _policy_rows(shape)
        self.assertEqual(len(rows), 10)
        for name, mode, eigen, iteration, status, want, want_gap, n_sz, finite in rows:
            with self.subTest(row=name):
                report = sc._so_bond_report(mode, eigen, iteration, gap_shape=shape)
                self.assertEqual(report.gap_status, status)
                d = tempfile.mkdtemp(prefix="sc_so_bond_policy_")
                self.addCleanup(shutil.rmtree, d, True)
                match = (np.ones(len(report.vals), dtype=bool) if len(report.vals) else None)
                with _quiet_logs():
                    sc._so_output_record(d, "gap.dat", p.norb_phys, None, None, None, match,
                                         p.kx, p.ky, p.kz, sz_residual=0.0,
                                         provenance_members=members, bond_report=report,
                                         block_layout="sz_packed")
                with np.load(os.path.join(d, "gap.npz"), allow_pickle=False) as z:
                    m = {k: z[k] for k in z.files}
                self.assertEqual(set(m), want)
                self.assertEqual(str(m["gap_status"]), status)
                self.assertTrue(np.array_equal(m["gap"], want_gap))
                self.assertEqual(m["gap_spin_weights"].shape, (4,))
                for c in _SZ_COMPONENTS:
                    w = m["sz_weight_" + c]
                    self.assertEqual(w.shape, (n_sz,), c)
                    if finite:
                        self.assertTrue(np.all(np.isfinite(w)), c)
                    else:
                        self.assertTrue(np.all(np.isnan(w)), c)
                if finite:
                    self.assertTrue(np.all(np.isfinite(m["gap_spin_weights"])))
                else:
                    self.assertTrue(np.all(np.isnan(m["gap_spin_weights"])))
                if "pp_sectors" in m:
                    self.assertEqual(m["pp_sectors"].dtype, np.int8)
                    self.assertEqual(m["pp_sectors"].shape, (len(report.vals),))

    def test_pp_labels_agree_with_sz_weights(self):
        """Through the entry, every listed ``pp_sectors`` label and the
        written gap's ``pp_sector`` carry the pair S_z of their content: +1
        all up-up, -1 all down-down, 0 all S_z = 0 (singlet plus triplet), to
        1e-12; a dense result is labelled 2. The packed runs are the B = 1
        Zeeman listing (every sector listed) and the default Zeeman run
        (B = 5); the dense one is the B = 1 spin-orbit coupled listing."""
        def content(w, p):
            return {1: w["upup"], -1: w["dndn"], 0: w["singlet0"] + w["triplet0"]}[p]

        b1 = {"Bond_Channels": True, "Bond_Coeff_Tail": 0, "Num_Eigenvalues": _B1_NEV}
        for fx, eli in ((ZEEMAN_EXTERN_ONSITE, b1), (F.ZEEMAN_EXTERN, None)):
            with self.subTest(fixture=fx.name):
                g = F.run_case(fx, eliashberg=eli).gap_npz
                self.assertEqual(str(g["block_layout"]), "sz_packed")
                labels = [int(p) for p in g["pp_sectors"]]
                if fx is ZEEMAN_EXTERN_ONSITE:
                    self.assertEqual(set(labels), {1, 0, -1})
                for i, p in enumerate(labels):
                    w = {c: float(g["sz_weight_" + c][i]) for c in _SZ_COMPONENTS}
                    self.assertGreaterEqual(content(w, p), 1.0 - 1e-12, (i, p, w))
                pp = int(g["pp_sector"])
                w = so_pairing.so_sz_decomposition(g["gap"], fx.norb_phys)
                self.assertGreaterEqual(content(w, pp), 1.0 - 1e-12, (pp, w))
        g = F.run_case(F.SOC_ONSITE, eliashberg=b1).gap_npz
        self.assertEqual(str(g["block_layout"]), "dense")
        self.assertEqual(set(int(p) for p in g["pp_sectors"]), {2})
        self.assertEqual(int(g["pp_sector"]), 2)

    def test_packed_iteration_note(self):
        """The eigenvalue.dat note of an S_z-packed power iteration describes
        the reported sector's iteration: its Rayleigh residual (a number) and
        its own step count."""
        seen = []
        real = sc._so_iterate_sectors

        def spy(*args, **kw):
            seen.append(real(*args, **kw))
            return seen[-1]
        with mock.patch.object(sc, "_so_iterate_sectors", spy):
            run = F.run_case(F.ZEEMAN_EXTERN, fresh=True, param={"Nmat": 8},
                             eliashberg={"Solver_Mode": "iteration"})
        it, = seen
        self.assertTrue(it.converged)
        lead = next(x for x in it.leaders if x.pp_sector == it.pp_sector)
        note, = [ln for ln in run.eigenvalue_dat.splitlines() if "solver_mode=" in ln]
        resid = re.search(r"with residual \|\|K v - lambda v\|\|/\|\|v\|\| = (\S+); ", note)
        self.assertIsNotNone(resid, note)
        self.assertEqual(resid.group(1), "{:.3e}".format(lead.info["rayleigh_residual"]))
        self.assertTrue(note.endswith("the iteration converged after {} steps.".format(
            lead.n_iter)), note)

    def test_archive_passes_its_structural_check(self):
        # the run's own gap.npz (the run cache keeps the members only)
        for layout, fx in _SOLVE.items():
            with self.subTest(layout=layout):
                d = tempfile.mkdtemp(prefix="sc_so_bond_archive_")
                self.addCleanup(shutil.rmtree, d, True)
                input_dict = F.write_case(d, fx)
                with _quiet_logs():
                    sc.calc_eliashberg(input_dict)
                path = os.path.join(input_dict["file"]["output"]["path_to_output"],
                                    "gap.npz")
                with np.load(path, allow_pickle=False) as z:
                    self.assertEqual(str(z["block_layout"]), layout)
                    self.assertEqual(
                        fingerprint.validate_structure(z, path, kind="gap", archive=path),
                        ("gap", "sc_so_bond_static", 2))

    def test_eigenvalue_dat_keys(self):
        delta_r = "[[0, 0, 0], [-1, 0, 0], [0, -1, 0], [0, 1, 0], [1, 0, 0]]"
        reasons = {
            "dense": "so bond layout: dense (reason: spin-flip transfer entry R=(-1, 0, 0) "
                     "(1, 2) = (-0.3+0j))",
            "sz_packed": _PACKED_REASON}
        selections = {"dense": "LR_retry", "sz_packed": "+1:dense-LR;0:dense-LR;-1:dense-LR"}
        for layout, fx in _SOLVE.items():
            with self.subTest(layout=layout):
                run = F.run_case(fx, eliashberg=_FIXED_CAP)
                text = run.eigenvalue_dat
                keys = _dat_keys(text)
                d = tempfile.mkdtemp(prefix="sc_so_bond_keys_")
                self.addCleanup(shutil.rmtree, d, True)
                p = F.prepare(fx, d, eliashberg=_FIXED_CAP)
                with _quiet_logs():
                    adm = _admit(p)
                m = run.gap_npz
                want = {
                    "bond_channels": "true",
                    "bond_n_channels": "5",
                    "bond_delta_r": delta_r,
                    "bond_max_shells": "unset",
                    "so_bond_layout": layout,
                    "so_bond_layout_reason": reasons[layout],
                    "bond_memory_cap": adm.cap_summary,
                    "g2_tail": "true",
                    "bond_coeff_tail": "1.0",
                    "bond_coeff_tail_source": "default",
                    "approximation": "static RPA-ladder bond dressing",
                    "gap_status": "certified",
                    "pp_sector": str(int(m["pp_sector"])),
                    "producer_stage": "sc_so_bond_static",
                    "producer_label": "hwave_sc (spin-orbital, bond, static)",
                    "eigen_converged": "true",
                }
                for key, value in want.items():
                    self.assertEqual(keys.get(key), value, key)
                    self.assertEqual(_dat_count(text, key), 1, key)
                self.assertEqual(adm.cap_summary,
                                 "2.000 GiB ([eliashberg] bond_memory_cap_gb = 2.0)")
                self.assertEqual(int(m["pp_sector"]), 2 if layout == "dense" else 0)
                self.assertEqual(keys["selection"], selections[layout])
                self.assertEqual(set(keys), _DAT_KEY_SET)
                self.assertEqual(_dat_count(text, "iteration_converged"), 0)
                self.assertNotIn("lambda_pp", keys)
                # the unset cap of the default run names the available memory
                default = _dat_keys(F.run_case(fx).eigenvalue_dat)["bond_memory_cap"]
                self.assertRegex(default, _UNSET_CAP)

    def test_onsite_path_never_enters_bond_code(self):
        solve = mock.Mock(side_effect=AssertionError("the bond solve was entered"))
        admit = mock.Mock(side_effect=AssertionError("the bond admission was entered"))
        with mock.patch.object(sc, "_solve_so_bond_static", solve), \
                mock.patch.object(sc, "_admit_so_bond_static", admit):
            run = F.run_case(F.SOC_ONSITE, fresh=True, eliashberg=_ONSITE_ELI)
        self.assertEqual((solve.call_count, admit.call_count), (0, 0))
        self.assertTrue(np.isfinite(run.value))
        self.assertNotIn("block_layout", run.gap_npz)
        self.assertNotIn("bond_channels", _dat_keys(run.eigenvalue_dat))

    def test_empty_gap_space_warns(self):
        with patched_sectors(empty=(1, 0, -1)):
            run = F.run_case(F.ZEEMAN_EXTERN, fresh=True,
                             eliashberg={"Solver_Mode": "iteration"})
        self.assertEqual(
            sc._SO_BOND_EMPTY_GAP,
            "spin-orbital bond pairing: no pair-S_z sector has a fermion-antisymmetric gap "
            "space on this k grid; no gap is reported (gap_status = \"none\")")
        msgs = [r.getMessage() for r in run.records if r.levelno == logging.WARNING]
        self.assertEqual(msgs.count(sc._SO_BOND_EMPTY_GAP), 1, msgs)
        m = run.gap_npz
        self.assertTrue(np.isnan(run.value))
        self.assertEqual(str(m["gap_status"]), "none")
        self.assertEqual(set(m), _LIT_BASE)
        self.assertEqual(np.count_nonzero(m["gap"]), 0)
        self.assertTrue(np.all(np.isnan(m["gap_spin_weights"])))
        for c in _SZ_COMPONENTS:
            self.assertEqual(m["sz_weight_" + c].shape, (1,))
            self.assertTrue(np.isnan(m["sz_weight_" + c][0]), c)
        # no iteration ran: the note says so (no step count)
        note = ("# solver_mode='iteration': no pair-S_z sector has a fermion-antisymmetric "
                "gap space on this k grid, so no power iteration ran; NO eigenvalue is "
                "reported (the value below is NaN) and the gap written is zeros "
                "(gap_status = \"none\").")
        lines = run.eigenvalue_dat.splitlines()
        self.assertIn(note, lines)
        self.assertNotIn("steps", run.eigenvalue_dat)
        # the eigenvalue analysis of an empty space: the same one WARNING
        with patched_sectors(empty=(1, 0, -1)):
            run = F.run_case(F.ZEEMAN_EXTERN, fresh=True,
                             eliashberg={"Solver_Mode": "eigenvalue"})
        self.assertEqual(str(run.gap_npz["gap_status"]), "none")
        self.assertEqual(_dat_keys(run.eigenvalue_dat)["selection"],
                         "+1:empty;0:empty;-1:empty")
        self.assertEqual(empty_gap_warnings(run), 1)

    def test_failed_sector_is_not_an_empty_gap_space(self):
        # pair-S_z 0 failed, +1 and -1 empty: no gap is reported, but the gap
        # space is not empty, so the empty-space WARNING stays silent
        with patched_sectors(empty=(1, -1)), eigen_fails(ZEEMAN_DIM_PP0):
            run = F.run_case(F.ZEEMAN_EXTERN, fresh=True,
                             eliashberg={"Solver_Mode": "eigenvalue"})
        self.assertEqual(str(run.gap_npz["gap_status"]), "none")
        self.assertTrue(np.isnan(run.value))
        self.assertEqual(_dat_keys(run.eigenvalue_dat)["selection"],
                         "+1:empty;0:not-converged;-1:empty")
        self.assertEqual(empty_gap_warnings(run), 0)


class TestReportingThroughEntry(unittest.TestCase):
    """Solver modes and failures through the public entry: a forced state of
    the sector solvers reaches gap.npz, eigenvalue.dat and the exit status
    (the full matrix is the heavy module's)."""

    def test_state_rows_smoke(self):
        # a null pair-S_z +1 sector next to converged ones, in "both": the
        # leader stays the unpatched run's
        lam = F.run_case(F.ZEEMAN_EXTERN).value
        run, keys, status = reporting_cell(self, F.ZEEMAN_EXTERN, "both",
                                           patches=(patched_sectors(null=(1,)),))
        check_reporting_cell(self, run, keys, status, (lam, "certified", 0, True, True, 0))
        self.assertEqual(keys["selection"], "+1:null-operator;0:dense-LR;-1:dense-LR")


# ------------------------------------------------------------ SU(2) reduction
def _su2_inputs(fixture, layout):
    """``(p, geom_info, interactions)``: the prepared spin-orbital run of
    ``fixture`` with the bridge configuration (``g2_tail = false``,
    ``bond_coeff_tail = 0.5``) in ``layout``, and the Geometry and the
    interaction tables of its files as the production reader returns them
    (physical orbitals: the same tables a non-spin-orbital run reads)."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_su2_")
    try:
        eli = dict(_BRIDGE_ELI)
        if layout != "auto":
            eli["SO_Bond_Layout"] = layout
        p = F.prepare(fixture, d, eliashberg=eli)
        geom_info, _, interactions = sc._read_interaction_files(p.input_dict)
    finally:
        shutil.rmtree(d, True)
    return p, geom_info, interactions


def nonso_kernels(test, p, geom_info, interactions):
    """``({pairing_type: K}, topo)``: the non-spin-orbital static bond
    operators ``sc._build_bond_operator`` (one call per pairing type) on the
    spin-up block of the run's band (the band is spin independent: asserted),
    with mu of the non-spin-orbital search, the run's ``g2_tail`` and
    ``bond_coeff_tail``, and the topology the non-spin-orbital resolver
    builds from the same tables."""
    n = p.norb_phys
    eps = p.epsilon_k
    test.assertEqual(float(np.abs(eps[:n, n:]).max()), 0.0)
    test.assertEqual(float(np.abs(eps[n:, n:] - eps[:n, :n]).max()), 0.0)
    ev, evec = sc._calc_eigenvalues(np.ascontiguousarray(eps[:n, :n]))
    mu = sc._determine_mu(ev, p.beta, p.filling, n)
    topo, _, split = sc._resolve_bond_topology(interactions, geom_info, n, None)
    carrier = sc._build_bond_green(ev, evec, mu, p.beta, p.nmat, p.cfg.coeff_tail, None,
                                   spin_orbital_mode=False)
    test.assertEqual(carrier.coeff_tail_applied, 0.5)
    ops = {}
    for eta in ("singlet", "triplet"):
        (A, size), _, _ = sc._build_bond_operator(
            topo, split, carrier, n, p.kx, p.ky, p.kz, p.beta, eta,
            bond_max_shells=None, bond_memory_cap_gb=8.0, g2_tail=p.cfg.g2_tail)
        test.assertEqual(size, n * n * int(np.prod(p.shape)))
        ops[eta] = A
    return ops, topo


def parity_bases(norb, shape):
    """``(E_+, E_-, n_self)`` of the non-spin-orbital gap ``(norb, norb, Nx,
    Ny, Nz)`` under ``P: (a, b, k) -> (b, a, -k)``: ``E_+`` has the columns
    ``(e_u + e_Pu) / sqrt 2`` for ``u < Pu`` and ``e_u`` for a self-paired
    ``u``, ``E_-`` the columns ``(e_u - e_Pu) / sqrt 2`` for ``u < Pu``."""
    gshape = (norb, norb) + tuple(shape)
    n = int(np.prod(gshape))
    partner = sc._reverse_k_and_orbital(np.arange(n).reshape(gshape)).ravel()
    plus, minus = [], []
    for u in range(n):
        pu = int(partner[u])
        if u == pu:
            c = np.zeros(n)
            c[u] = 1.0
            plus.append(c)
        elif u < pu:
            c = np.zeros(n)
            c[u] = c[pu] = 1.0 / np.sqrt(2.0)
            plus.append(c)
            c = np.zeros(n)
            c[u], c[pu] = 1.0 / np.sqrt(2.0), -1.0 / np.sqrt(2.0)
            minus.append(c)
    n_self = int(np.sum(partner == np.arange(n)))
    return np.array(plus).T, np.array(minus).T, n_self


def _projected(A, E):
    """``E^dagger A E``, one matvec per column."""
    return E.conj().T @ np.column_stack([np.asarray(A.matvec(E[:, j].astype(complex)))
                                         for j in range(E.shape[1])])


def so_reduced_matrices(static):
    """The dense reduced spin-orbital kernel: ``[basis.reduce(operator)]``
    (dense layout) or the reduced matrix of every pair-S_z sector
    (``sector_op(p)``, S_z-packed layout)."""
    if static.layout == "sz_packed":
        out = []
        for pp in static.kernel.pp_sectors:
            op = static.kernel.sector_op(pp)
            out.append(op.basis.reduce(op.operator).matmat(
                np.eye(op.basis.dim, dtype=complex)))
        return out
    return [static.basis.reduce(static.operator).matmat(
        np.eye(static.basis.dim, dtype=complex))]


def su2_match(test, fixture, layout):
    """The SU(2) matching of the spin-orbital bond kernel of ``fixture`` in
    ``layout`` against the non-spin-orbital singlet and triplet bond kernels
    (the conditions are asserted; returns the measures): ``mu`` = the
    eigenvalues of ``K_s = E_+^dagger K_singlet E_+`` once and of ``K_t =
    E_-^dagger K_triplet E_-`` three times, matched one to one against the
    eigenvalues ``lambda`` of the reduced spin-orbital kernel (the direct sum
    of the pair-S_z sectors in the packed layout) by ``match_spectrum``."""
    import scipy.linalg
    from tests.so_bond_fixtures import KAPPA_MAX, dense_spectrum, match_spectrum
    with _quiet_logs():
        p, geom_info, interactions = _su2_inputs(fixture, layout)
        _, _, static = _chain(p)
        K, topo = nonso_kernels(test, p, geom_info, interactions)
    test.assertEqual(static.layout, layout if layout != "auto" else "sz_packed")
    test.assertEqual(bond_channels.non_su2_bonds(topo), [])
    n = p.norb_phys
    n_phi = n * n * int(np.prod(p.shape))
    Ep, Em, n_self = parity_bases(n, p.shape)
    test.assertEqual(n_self, n * so_bond.n_trim(p.shape))
    Ks, Kt = _projected(K["singlet"], Ep), _projected(K["triplet"], Em)
    mats = so_reduced_matrices(static)
    full = sc._SOAntisymmetricBasis((2 * n, 2 * n) + tuple(p.shape))
    test.assertEqual(full.dim, 2 * n_phi - n_self)
    test.assertEqual(sum(m.shape[0] for m in mats), full.dim)
    test.assertEqual(Ks.shape[0] + 3 * Kt.shape[0], full.dim)
    K_so = scipy.linalg.block_diag(*mats)
    spec = dense_spectrum(K_so)
    spec_s, spec_t = np.linalg.eigvals(Ks), np.linalg.eigvals(Kt)
    mu = np.concatenate([spec_s] + [spec_t] * 3)
    test.assertTrue(np.all(np.isfinite(spec.vals)) and np.all(np.isfinite(mu)))
    test.assertLessEqual(float(spec.kappa.max()), KAPPA_MAX)
    # every spin-orbital eigenpair: ||K v - lambda v|| <= 1e-10 ||K||_F ||v||
    knorm = np.linalg.norm(K_so)
    for i in range(len(spec.vals)):
        v = spec.right[:, i]
        test.assertLessEqual(np.linalg.norm(K_so @ v - spec.vals[i] * v),
                             1e-10 * knorm * np.linalg.norm(v), i)
    cols, ratio = match_spectrum(spec, mu)
    miss = float(np.abs(mu - spec.vals[cols]).max())
    out = {"ratio": ratio, "miss": miss / spec.scale, "kappa": float(spec.kappa.max()),
           "lead_so": float(spec.vals.real.max()),
           "lead_nonso": float(np.concatenate([spec_s, spec_t]).real.max())}
    if static.layout == "sz_packed":
        # each pair-S_z sector on its own: +-1 (equal spin) carry the triplet
        # spectrum, 0 the singlet and the S_z = 0 triplet
        want = {1: spec_t, 0: np.concatenate([spec_s, spec_t]), -1: spec_t}
        test.assertEqual(tuple(static.kernel.pp_sectors), (1, 0, -1))
        for pp, M in zip(static.kernel.pp_sectors, mats):
            test.assertEqual(M.shape[0], len(want[pp]), pp)
            sector = dense_spectrum(M)
            test.assertLessEqual(float(sector.kappa.max()), KAPPA_MAX, pp)
            out["ratio_pp%+d" % pp] = match_spectrum(sector, want[pp])[1]
    return out


def _vertex_at_minus_q(calls):
    """``so_bond.so_bond_vertex`` read at ``-q`` (both layouts); each call is
    appended to ``calls``."""
    from tests.so_bond_oracles import mesh_negate
    real = so_bond.so_bond_vertex

    def flipped(compiled, topo, shape, **kw):
        g = real(compiled, topo, shape, **kw)
        neg = mesh_negate(shape)
        calls.append(shape)
        if isinstance(g, so_bond.PackedArrays):
            return so_bond.PackedArrays(g.sector_map,
                                        tuple(np.ascontiguousarray(b[neg]) for b in g.blocks))
        return np.ascontiguousarray(g[neg])
    return flipped


def _su2_mixed_fixture():
    """The two-orbital 6x1x1 chain of the non-spin-orbital bond pairing
    oracle (mirror- and transpose-asymmetric hopping, spin independent) with
    the spin-rotation invariant interaction on mixed bonds of the
    spin-orbital pairing-vertex SU(2) reduction (inter-orbital CoulombInter
    different at +R and -R; Hund, Ising and Exchange = Hund - 2 Ising; real
    PairHop; on-site CoulombIntra)."""
    from tests import test_sc_bond_pairing_oracle as nonso_oracle
    from tests import test_so_bond_pairing_oracle as so_oracle
    hopping, transfer = [], []
    for (r, a, b), t in sorted(nonso_oracle.TRANSFER[2].items()):
        R = (r, 0, 0)
        if (a, b) == (0, 0):
            hopping.append((R, t))
        else:
            transfer.extend((R, (2 * a + s, 2 * b + s), t) for s in (0, 1))
    interactions = tuple(
        (typ, tuple(sorted((tuple(R), tuple(ab), v) for (R, ab), v in table.items())))
        for typ, table in so_oracle.SU2_MIXED.items())
    return F.EntryFixture("SU2_MIXED_6", (6, 1, 1), 2, tuple(hopping), tuple(transfer), (),
                          0.0, interactions, 0.4)


#: the S_z-packed and the forced dense layout of an S_z-conserving input
_SU2_LAYOUTS = ("auto", "dense")


class TestSU2(unittest.TestCase):
    """Gate A2-SU2: with spin-independent hopping and a spin-rotation
    invariant interaction (CoulombIntra, CoulombInter, PairHop and ``J_Ex =
    J_H - 2 J_I`` on every bond), the spectrum of the spin-orbital bond
    kernel is that of the non-spin-orbital singlet kernel on the even pairs
    once and of the triplet kernel on the odd pairs three times, in both
    layouts; the public runs report ``max(lambda_s, lambda_t)`` with the
    configured ``g2_tail`` and ``bond_coeff_tail`` (the configuration
    resolution).

    Blind spots and the gates that cover them: the S_z-conserving fixtures
    never fill the spin-flip slots of the vertex (the one-channel cross-path
    gate and the input-independent bridge on the spin-orbit coupled fixture
    cover the dense spin-flip layout); a relabelling of the pair-S_z
    sectors that keeps their spectra passes (the entry-level check of the
    reported sector labels against the S_z weights pins the labels); none
    of the gates runs past the RPA instability (every fixture is below
    it)."""

    def test_fast_both_layouts(self):
        for layout in _SU2_LAYOUTS:
            with self.subTest(layout=layout):
                m = su2_match(self, F.SU2_4X4, layout)
                self.assertLessEqual(max(v for k, v in m.items() if k.startswith("ratio")),
                                     1.0, m)

    def test_end_to_end_lambda(self):
        eli = {"Bond_Channels": True, "G2_Tail": False, "Bond_Coeff_Tail": 0.5}
        nonso = {p: F.run_case(F.SU2_4X4, mode={"enable_spin_orbital": False},
                               eliashberg=dict(eli, Pairing_Type=p))
                 for p in ("singlet", "triplet")}
        for p, run in nonso.items():
            lines = run.eigenvalue_dat.splitlines()
            for line in ("# eigen_converged = true", "# g2_tail = False",
                         "# bond_coeff_tail_applied = 0.5", "# bond_coeff_tail_source = config"):
                self.assertIn(line, lines, p)
            self.assertIn("'{}' parity sector".format(p), run.eigenvalue_dat)
            self.assertTrue(np.isfinite(run.value), p)
        lam = max(nonso["singlet"].value, nonso["triplet"].value)
        # the leader of this fixture is the triplet: the comparison is not
        # decided by one channel alone being right
        self.assertGreater(nonso["triplet"].value, nonso["singlet"].value)
        for layout in ("dense", "auto"):
            with self.subTest(layout=layout):
                run = F.run_case(F.SU2_4X4, layout=layout, eliashberg=eli)
                lines = run.eigenvalue_dat.splitlines()
                for line in ("# g2_tail = false", "# bond_coeff_tail = 0.5",
                             "# bond_coeff_tail_source = config"):
                    self.assertIn(line, lines)
                self.assertEqual(str(run.gap_npz["block_layout"]),
                                 "dense" if layout == "dense" else "sz_packed")
                self.assertTrue(bool(run.gap_npz["eigen_converged"]))
                self.assertTrue(np.isfinite(run.value))
                self.assertLessEqual(abs(run.value - lam), 1e-8 * max(1.0, abs(lam)),
                                     (run.value, lam))

    def test_non_vacuity(self):
        for layout in _SU2_LAYOUTS:
            with self.subTest(layout=layout):
                calls = []
                with mock.patch.object(so_bond, "so_bond_vertex", _vertex_at_minus_q(calls)):
                    m = su2_match(self, F.SU2_4X4, layout)
                self.assertEqual(len(calls), 1)
                self.assertGreaterEqual(m["miss"], 1e-3, m)

    def test_norb2_both_layouts(self):
        # two orbitals on the 6x1x1 mixed-bond chain (inter-orbital rows,
        # CoulombInter different at +R and -R)
        fx = _su2_mixed_fixture()
        for layout in _SU2_LAYOUTS:
            with self.subTest(layout=layout):
                m = su2_match(self, fx, layout)
                self.assertLessEqual(max(v for k, v in m.items() if k.startswith("ratio")),
                                     1.0, m)


# ------------------------------------------------------------ one channel
def _onsite_rows(fixture, name):
    """``fixture`` with its R = 0 interaction rows only (B = 1)."""
    return fixture._replace(name=name, interactions=tuple(
        (t, tuple(r for r in rows if r[0] == (0, 0, 0)))
        for t, rows in fixture.interactions if any(r[0] == (0, 0, 0) for r in rows)))


#: the band and field of ZEEMAN_EXTERN with its on-site interaction rows only
ZEEMAN_EXTERN_ONSITE = _onsite_rows(F.ZEEMAN_EXTERN, "ZEEMAN_EXTERN_ONSITE")
#: num_eigenvalues >= dim - 1 = 27 of the 4x4 antisymmetric space: both
#: paths diagonalize densely (and so does every pair-S_z sector)
_B1_NEV = 27
_ONSITE_DENSE_INFO = ("spin-orbital: num_eigenvalues = 27 >= dim - 1 = 27 of the "
                      "antisymmetric space; dense diagonalization")


class TestCrossPathB1(unittest.TestCase):
    """Gate A2-B1: the spin-orbital bond path with on-site rows only (B = 1)
    and no tail on the bubble's Green function (``bond_coeff_tail = 0`` and
    the RPA ``coeff_tail = 0``; ``g2_tail`` is on, the default, on both
    paths) equals the spin-orbital on-site path (``chi0q_mode = "calc"``) --
    the leading eigenvalue, the leading gap and the whole 27-value listing
    (every solve dense, so the listings are complete up to one value) -- on
    a spin-flip SOC band (dense layout) and on a Zeeman band (S_z-packed
    layout). The two paths share no kernel code (bubble, vertex, kernel
    operator): the RPA chi0q against the bond bubble, the on-site pairing
    vertex against the bond pairing vertex parts, the on-site kernel operator
    against the bond pair kernel. The leader alone does not see the bare
    vertex of an on-site interaction (k independent, it couples only to the
    local pair amplitude); the listing does.

    Blind spots and the gates that cover them: the two paths share the
    band, mu, the antisymmetric basis, the dense eigensolver and the
    reporting (validated by the exact-diagonalization gates of those
    functions); B = 1 has no off-site block (the kernel bridges and the
    SU(2) gate cover them); none of the gates runs past the RPA instability
    (every fixture is below it)."""

    def test_b1_matches_onsite_path(self):
        from tests import so_run_compare
        for fx, layout in ((F.SOC_ONSITE, "dense"), (ZEEMAN_EXTERN_ONSITE, "sz_packed")):
            with self.subTest(fixture=fx.name):
                self.assertTrue(all(r[0] == (0, 0, 0) for _, rows in fx.interactions
                                    for r in rows))
                a = F.run_case(fx, eliashberg={"Bond_Channels": True, "Bond_Coeff_Tail": 0,
                                               "Num_Eigenvalues": _B1_NEV})
                b = F.run_case(fx, eliashberg={"Bond_Channels": False, "Chi0q_Mode": "calc",
                                               "Num_Eigenvalues": _B1_NEV})
                ga, gb = a.gap_npz, b.gap_npz
                self.assertEqual(str(ga["producer_stage"]), "sc_so_bond_static")
                self.assertEqual(str(gb["producer_stage"]), "sc_so_static")
                self.assertEqual(str(ga["block_layout"]), layout)
                self.assertEqual(np.asarray(ga["delta_r"]).tolist(), [[0, 0, 0]])
                for key in ("state_beta", "state_nmat", "state_target"):
                    self.assertEqual(ga[key], gb[key], key)
                self.assertLessEqual(abs(float(ga["state_mu"]) - float(gb["state_mu"])), 1e-12)
                # the solves: dense branches on both paths
                sel = [l.split(" = ", 1)[1] for l in a.eigenvalue_dat.splitlines()
                       if l.startswith("# selection = ")]
                self.assertEqual(len(sel), 1)
                for part in sel[0].split(";"):
                    self.assertEqual(part.split(":")[-1], "dense-LR", sel)
                self.assertTrue(any(r.getMessage().startswith(_ONSITE_DENSE_INFO)
                                    for r in b.records))
                for g in (ga, gb):
                    self.assertTrue(bool(g["eigen_converged"]))
                    self.assertEqual(len(g["eigenvalues"]), _B1_NEV)
                lam = b.value
                self.assertTrue(np.isfinite(lam) and lam > 0.0)
                self.assertLessEqual(abs(a.value - lam), 1e-10 * max(1.0, abs(lam)),
                                     (a.value, lam))
                self.assertTrue(so_run_compare.same_leading_vector(
                    ga["eigenvalues"], ga["gap"], gb["gap"], vecs_a=ga["eigenvectors"],
                    vecs_b=gb["eigenvectors"]))
                # the whole listing, both sorted by descending real part
                la, lb = (np.asarray(g["eigenvalues"]) for g in (ga, gb))
                la, lb = (v[np.argsort(-v.real, kind="stable")] for v in (la, lb))
                self.assertTrue(so_run_compare.eigenvalues_close(la, lb, tol=1e-10),
                                float(np.abs(la - lb).max()))


# ------------------------------------------------------------ zero field
#: the field sequence of the zero-field limit (solved in this order)
ZF_H = (0.0, 0.5e-3, 1e-3, 2e-3, 4e-3)


def _x_rows(rows):
    """The on-site and +-x rows of an interaction table."""
    return tuple(r for r in rows if r[0] in ((0, 0, 0), (1, 0, 0), (-1, 0, 0)))


#: the 4x2 mesh with the hopping of SU2_4X4 and its spin-rotation invariant
#: interaction on site and on the x bond only (B = 3: +-y channels would
#: coincide on L_y = 2), filling 0.35 (no particle-hole symmetry), the field
#: as an Extern file (H1(R = 0) = 1, coeff_extern = h; 0 here)
ZF_4X2 = F.EntryFixture(
    "ZF_4X2", (4, 2, 1), 1, F.SU2_4X4.hopping, (), (((0, 0, 0), (0, 0), 1.0),), 0.0,
    tuple((t, _x_rows(rows)) for t, rows in F.SU2_4X4.interactions), 0.35)


def _at_field(fixture, h):
    """``fixture`` with the field ``h`` (spin up ``+h``): the Extern
    coefficient, or ``+-h`` on site in Transfer for a fixture without an
    Extern table."""
    name = "{}_h{!r}".format(fixture.name, h)
    if fixture.extern:
        return fixture._replace(name=name, coeff_extern=h)
    n = fixture.norb_phys
    rows = tuple(((0, 0, 0), (2 * a + s, 2 * a + s), (h if s == 0 else -h))
                 for a in range(n) for s in (0, 1))
    return fixture._replace(name=name, transfer=fixture.transfer + rows)


def reduced_spectrum_inputs(static):
    """``(K, E)``: the dense reduced spin-orbital kernel (the direct sum of
    the pair-S_z sectors in the S_z-packed layout, :func:`so_reduced_matrices`)
    and the isometry ``E`` from its coordinates to the full flat gap."""
    import scipy.linalg
    mats = so_reduced_matrices(static)
    if static.layout == "sz_packed":
        bases = [static.kernel.sector_op(pp).basis for pp in static.kernel.pp_sectors]
    else:
        bases = [static.basis]
    E = np.hstack([b.scatter(b.to_full(np.eye(b.dim, dtype=complex))) for b in bases])
    return scipy.linalg.block_diag(*mats), E


def zero_field_chain(test, fixture, layout, h, *, oracle=True):
    """``(K, E)`` (:func:`reduced_spectrum_inputs`) of the static operator of
    ``fixture`` at the field ``h`` in ``layout`` (``"auto"`` = S_z-packed):
    its own prepared run, band, Green function, bubble, dressing, vertex,
    G2 and kernel; with ``oracle`` the Green function is checked to carry
    ``+h`` on spin up (the analytic oracle, one physical orbital)."""
    fx = _at_field(fixture, h)
    d = tempfile.mkdtemp(prefix="sc_so_bond_zf_")
    try:
        eli = {} if layout == "auto" else {"SO_Bond_Layout": layout}
        with _quiet_logs():
            p = F.prepare(fx, d, eliashberg=eli)
    finally:
        shutil.rmtree(d, True)
    with _quiet_logs():
        so_green, _, static = _chain(p)
    test.assertEqual(static.layout, "sz_packed" if layout == "auto" else layout)
    if oracle and fx.norb_phys == 1:
        test.assertLessEqual(F.zeeman_sign_miss(so_green.carrier.full_sc, so_green.mu,
                                                fixture=fx, h=h, beta=p.beta, nmat=p.nmat),
                             1e-12, h)
    return reduced_spectrum_inputs(static)


def zero_field_spectra(test, fixture, layout, H=ZF_H, *, oracle=True):
    """``(lams, vecs0, Ks)``: every eigenvalue of the reduced kernel at each
    field of ``H`` (rows, in the order of ``H``), the eigenvectors at
    ``H[0]`` in the full flat gap (columns) and the reduced kernels."""
    lams, vecs0, Ks = [], None, []
    for k, h in enumerate(H):
        K, E = zero_field_chain(test, fixture, layout, h, oracle=oracle)
        vals, vecs = np.linalg.eig(K)
        test.assertTrue(np.all(np.isfinite(vals)), h)
        lams.append(vals)
        Ks.append(K)
        if k == 0:
            vecs0 = E @ vecs
    return np.array(lams), vecs0, Ks


def layouts_match(test, K_dense, K_packed):
    """The largest ratio of the matching rule between the spectra of the
    dense and the S_z-packed reduced kernels (asserted <= 1)."""
    from tests.so_bond_fixtures import KAPPA_MAX, dense_spectrum, match_spectrum
    spec = dense_spectrum(K_dense)
    test.assertLessEqual(float(spec.kappa.max()), KAPPA_MAX)
    ratio = match_spectrum(spec, np.linalg.eigvals(K_packed))[1]
    test.assertLessEqual(ratio, 1.0)
    return ratio


def zero_field_gate(test, fixture):
    """The zero-field limit of ``fixture`` in both layouts and the two
    layouts matched at each field; returns the measures."""
    out, Ks = {}, {}
    shape = (fixture.nd, fixture.nd) + tuple(fixture.shape)
    for layout in ("dense", "auto"):
        with test.subTest(layout=layout):
            lams, vecs0, Ks[layout] = zero_field_spectra(test, fixture, layout)
            out[layout] = Z.zero_field_measures(lams, vecs0, shape, fixture.norb_phys, ZF_H,
                                                kernels=Ks[layout][:2])
            Z.assert_zero_field(test, out[layout])
    out["match"] = max(layouts_match(test, a, b) for a, b in zip(Ks["dense"], Ks["auto"]))
    return out


#: the field of SU2_4X4 as an Extern file (coeff_extern = h; 0 here),
#: filling 0.4 with t2 (no particle-hole symmetry)
ZF_4X4 = F.SU2_4X4._replace(name="ZF_4X4", extern=(((0, 0, 0), (0, 0), 1.0),))
ZF_4X4_B1 = _onsite_rows(ZF_4X4, "ZF_4X4_B1")


def _zf_norb2():
    """The two-orbital square lattice of the weak-field gate of the
    spin-orbital on-site path on 4x4 (``tests/test_sc_spin_orbital.py``):
    nearest-neighbour hopping -1.0 (orbital 0), -0.7 (orbital 1) and 0.15
    (inter-orbital), on-site 0.3 on orbital 1, all spin independent; the
    on-site Kanamori interaction (CoulombIntra 1.6 and 1.9, CoulombInter
    0.5, Hund 0.12, Ising 0.04, Exchange = Hund - 2 Ising, PairHop 0.05);
    on the nearest-neighbour bonds CoulombInter 0.2 on both orbitals and
    the inter-orbital Hund 0.06, Ising 0.02, Exchange = Hund - 2 Ising; the
    field as Extern on both orbitals; filling 0.35."""
    nn = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
    o = (0, 0, 0)
    transfer = []

    def add(R, a, b, t):
        transfer.extend((R, (2 * a + s, 2 * b + s), t) for s in (0, 1))
    add(o, 1, 1, 0.3)
    for R in nn:
        add(R, 1, 1, -0.7)
        add(R, 0, 1, 0.15)
        add(R, 1, 0, 0.15)
    offd = ((0, 1), (1, 0))

    def rows(v, pairs, bonds=(o,)):
        return tuple((R, ab, v) for R in bonds for ab in pairs)
    interactions = (
        ("CoulombIntra", (((o), (0, 0), 1.6), ((o), (1, 1), 1.9))),
        ("CoulombInter", rows(0.5, offd) + rows(0.2, ((0, 0), (1, 1)), nn)),
        ("Hund", rows(0.12, offd) + rows(0.06, offd, nn)),
        ("Ising", rows(0.04, offd) + rows(0.02, offd, nn)),
        ("Exchange", rows(0.12 - 2 * 0.04, offd) + rows(0.06 - 2 * 0.02, offd, nn)),
        ("PairHop", rows(0.05, offd)),
    )
    return F.EntryFixture("ZF_4X4_NORB2", (4, 4, 1), 2, tuple((R, -1.0) for R in nn),
                          tuple(transfer), (((o), (0, 0), 1.0), ((o), (1, 1), 1.0)), 0.0,
                          interactions, 0.35)


ZF_4X4_NORB2 = _zf_norb2()
ZF_4X4_NORB2_B1 = _onsite_rows(ZF_4X4_NORB2, "ZF_4X4_NORB2_B1")


def entry_spectra(test, fixture, layout, H=ZF_H):
    """``(lams, vecs0)`` of the public runs of ``fixture`` at each field of
    ``H`` in ``layout``: the complete listings (``num_eigenvalues`` above
    the antisymmetric dimension, every solve dense) and the listed
    eigenvectors at ``H[0]`` (columns, full flat gap)."""
    dim = sc._SOAntisymmetricBasis((fixture.nd, fixture.nd) + tuple(fixture.shape)).dim
    lams, vecs0 = [], None
    for h in H:
        with _quiet_logs():
            run = F.run_case(_at_field(fixture, h), layout=layout,
                             eliashberg={"Num_Eigenvalues": dim + 2})
        g = run.gap_npz
        test.assertEqual(str(g["block_layout"]), "dense" if layout == "dense" else "sz_packed")
        test.assertTrue(bool(g["eigen_converged"]), h)
        vals = np.asarray(g["eigenvalues"])
        test.assertEqual(vals.shape, (dim,))
        test.assertTrue(np.all(np.isfinite(vals)), h)
        lams.append(vals)
        if vecs0 is None:
            vecs0 = np.asarray(g["eigenvectors"]).reshape(dim, -1).T
    return np.array(lams), vecs0


def entry_zero_field_gate(test, fixture):
    """The zero-field limit of ``fixture`` through the public entry in both
    layouts; at each field the two listings agree one to one (assignment,
    ``so_run_compare.eigenvalues_close``). Returns the measures."""
    from scipy.optimize import linear_sum_assignment
    from tests import so_run_compare
    shape = (fixture.nd, fixture.nd) + tuple(fixture.shape)
    out, lams = {}, {}
    for layout in ("dense", "auto"):
        with test.subTest(fixture=fixture.name, layout=layout):
            lams[layout], vecs0 = entry_spectra(test, fixture, layout)
            # the isolation of the selected clusters needs the kernels: the
            # operator of the same fixture at the first two fields
            kernels = [zero_field_chain(test, fixture, layout, h)[0] for h in ZF_H[:2]]
            out[layout] = Z.zero_field_measures(lams[layout], vecs0, shape, fixture.norb_phys,
                                                ZF_H, kernels=kernels)
            Z.assert_zero_field(test, out[layout])
    out["layouts"] = 0.0
    for a, b in zip(lams["dense"], lams["auto"]):
        r, c = linear_sum_assignment(np.abs(a[:, None] - b[None, :]))
        out["layouts"] = max(out["layouts"], float(np.max(
            np.abs(a[r] - b[c]) / np.maximum(1.0, np.abs(a[r])))))
        test.assertTrue(so_run_compare.eigenvalues_close(a[r], b[c]), out)
    return out


_NN = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))


def _zf_mixed_pair():
    """A two-orbital 4x4 model with on-site rows only (B = 1) whose leading
    clusters include a singlet and an S_z = 0 triplet 8.9e-4 apart: the field
    couples them at first order (``|C| = 0.047``), so ``h_4 |C| / gap`` is
    about 0.2 and the quadratic fit does not describe the pair, although
    their slopes alone call them isolated. Spin-independent hopping (the
    hopping of SU2_4X4 on orbital 0, -0.8 on orbital 1, inter-orbital -0.2
    on site and 0.15 along x, on-site 0.3 on orbital 1), CoulombIntra 0.83
    on both orbitals, the inter-orbital CoulombInter 0.5, the field as
    Extern on both orbitals, filling 0.4."""
    o = (0, 0, 0)
    transfer = []

    def add(R, a, b, t):
        transfer.extend((R, (2 * a + s, 2 * b + s), t) for s in (0, 1))
    add(o, 1, 1, 0.3)
    for R in _NN:
        add(R, 1, 1, -0.8)
    add(o, 0, 1, -0.2)
    add(o, 1, 0, -0.2)
    add((1, 0, 0), 0, 1, 0.15)
    add((-1, 0, 0), 1, 0, 0.15)
    interactions = (("CoulombIntra", ((o, (0, 0), 0.83), (o, (1, 1), 0.83))),
                    ("CoulombInter", ((o, (0, 1), 0.5), (o, (1, 0), 0.5))))
    return F.EntryFixture("ZF_MIXED_PAIR_B1", (4, 4, 1), 2, F.SU2_4X4.hopping,
                          tuple(transfer), ((o, (0, 0), 1.0), (o, (1, 1), 1.0)), 0.0,
                          interactions, 0.4)


ZF_MIXED_PAIR_B1 = _zf_mixed_pair()


class TestZeroField(unittest.TestCase):
    """Gate A2-0S3b: on an SU(2) model at a filling without particle-hole
    symmetry, a weak Zeeman field (Extern, ``coeff_extern = h``, ``+h`` on
    spin up) moves the eigenvalues of the spin-orbital bond kernel linearly:
    an SU(2) triplet (degeneracy 3m) splits into m slopes near zero (pair
    S_z = 0) and m pairs of opposite slopes (up-up / down-down), a singlet
    does not move at first order. ``h = 0, 0.5e-3, 1e-3, 2e-3, 4e-3`` are
    solved in increasing order, each with its own band, Green function,
    bubble, vertex and kernel; the paths are tracked by assignment and
    fitted by ``a h + b h^2`` (``tests/so_zero_field.py``). A selected
    cluster must also be isolated from every other cluster against the
    field-induced mixing (``r_ij <= 0.03``), else the fixture is rejected.
    Both layouts; at each field the two spectra agree by the matching rule.

    Blind spots and the gates that cover them: the weak-field conditions do
    not see the sign of the field (the analytic Green-function oracle pins
    it at every field); a field that reaches only one of the bubble and the
    pair propagator, and Gamma_B read at -q, keep the symmetry these
    conditions test (the input-independent kernel bridge sees the first,
    the SU(2) gate's non-vacuity test the second); no gate here looks past
    the RPA pole (every fixture is below the instability)."""

    def test_operator_both_layouts(self):
        m = zero_field_gate(self, ZF_4X2)
        self.assertLessEqual(m["match"], 1.0, m)

    def test_mixed_pair_is_rejected(self):
        """A fixture whose selected singlet / S_z = 0 triplet pair is mixed by
        the field fails the isolation rule (``r`` about 0.21, above 0.03) and
        so is rejected; the slope rule alone selects the pair."""
        lams, vecs0, Ks = zero_field_spectra(self, ZF_MIXED_PAIR_B1, "dense", oracle=False)
        shape = (ZF_MIXED_PAIR_B1.nd, ZF_MIXED_PAIR_B1.nd) + tuple(ZF_MIXED_PAIR_B1.shape)
        m = Z.zero_field_measures(lams, vecs0, shape, ZF_MIXED_PAIR_B1.norb_phys, ZF_H,
                                  kernels=Ks[:2])
        self.assertGreater(m["isolation"], Z.ISOLATION_MAX, m)
        with self.assertRaisesRegex(AssertionError, "the fixture is rejected: a selected "
                                                    "cluster is mixed by the field"):
            Z.assert_zero_field(self, m)

    def test_non_vacuity(self):
        """The field written as ``coeff_extern H1 (x) 1`` (sigma^z dropped) is
        absorbed by mu: no eigenvalue moves and the slope floor (condition
        2) fails, in both layouts (the Green-function oracle is off here).
        The replacement band is asserted to have received the Extern table
        at every nonzero field."""
        real = sc._so_epsilon_k
        fields = []

        def no_sigma_z(input_dict, kx, ky, kz, hr, extern_table, norb_phys):
            eps = real(input_dict, kx, ky, kz, hr, None, norb_phys)
            coeff = sc._so_coeff_extern(input_dict)
            fields.append((coeff, extern_table is not None))
            if extern_table is not None:
                n = norb_phys
                h1 = coeff * sc._build_hamiltonian_k(kx, ky, kz, extern_table, n)
                eps[:n, :n] += h1
                eps[n:, n:] += h1
            return eps
        shape = (ZF_4X2.nd, ZF_4X2.nd) + tuple(ZF_4X2.shape)
        for layout in ("dense", "auto"):
            with self.subTest(layout=layout):
                fields.clear()
                with mock.patch.object(sc, "_so_epsilon_k", no_sigma_z):
                    lams, vecs0, _ = zero_field_spectra(self, ZF_4X2, layout, oracle=False)
                self.assertEqual(sorted(set(fields)), [(h, True) for h in ZF_H])
                m = Z.zero_field_measures(lams, vecs0, shape, ZF_4X2.norb_phys, ZF_H)
                self.assertLess(m["floor"], 1e-6, m)


# ------------------------------------------------------------ non-SU(2) bonds
#: each S_z-conserving non-SU(2) type alone on the nearest-neighbour bonds
NONSU2_TYPES = (("Ising", 0.29), ("Hund", 0.53), ("Exchange", 0.31))


def nonsu2_fixture(typ, v):
    """4x4, one orbital, the spin-independent hopping of SU2_4X4, on-site
    CoulombIntra 0.83 and ``typ = v`` alone on the nearest-neighbour bonds,
    no field, no PairLift."""
    return F.EntryFixture("NONSU2_" + typ.upper(), F.SU2_4X4.shape, 1, F.SU2_4X4.hopping, (),
                          (), 0.0, (("CoulombIntra", (((0, 0, 0), (0, 0), 0.83),)),
                                    (typ, tuple((R, (0, 0), v) for R in _NN))), 0.4)


def _pair_spin_coords(basis, n):
    """The reduced coordinates of ``basis`` (whole gap ``(nd, nd, ...)``,
    spin-block index) by pair S_z: ``{+1: up-up, -1: down-down, 0: the
    rest}``, from the representative gap index of each coordinate."""
    A, B = np.unravel_index(basis.rep, basis.shape)[:2]
    sa, sb = A >= n, B >= n
    return {1: np.flatnonzero(~sa & ~sb), -1: np.flatnonzero(sa & sb),
            0: np.flatnonzero(sa != sb)}


def pair_sz_blocks(test, static, n):
    """``{pp: K_pp}``: the reduced kernel of each pair-S_z sector. Dense
    layout: ``basis.reduce`` restricted to the test-side coordinates of the
    sector (the blocks between sectors vanish: asserted at 1e-12 of the
    reduced kernel); S_z-packed layout: ``sector_op(pp)`` reduced."""
    if static.layout == "sz_packed":
        mats = so_reduced_matrices(static)
        return dict(zip(static.kernel.pp_sectors, mats))
    M = so_reduced_matrices(static)[0]
    coords = _pair_spin_coords(static.basis, n)
    test.assertEqual(sum(len(c) for c in coords.values()), M.shape[0])
    out = {pp: M[np.ix_(c, c)] for pp, c in coords.items()}
    off = M.copy()
    for c in coords.values():
        off[np.ix_(c, c)] = 0.0
    test.assertLessEqual(float(np.abs(off).max()), 1e-12 * float(np.abs(M).max()))
    return out


def nonsu2_measures(test, blocks, Ks, Kt):
    """The S_z-conserving non-SU(2) oracle on the pair-S_z blocks of one
    layout: the matching ratios of the spectra of ``K_+1`` and ``K_-1``
    against the non-spin-orbital triplet kernel ``K_t`` on the odd pairs
    (asserted <= 1, ``kappa_c <= KAPPA_MAX``), and ``d``, the largest
    assigned distance of the spectrum of ``K_0`` from that of ``K_s (+)
    K_t``, over ``max(1, max|lambda|)``."""
    import scipy.linalg
    from tests.so_bond_fixtures import KAPPA_MAX, dense_spectrum, match_spectrum
    spec_t = dense_spectrum(Kt)
    test.assertLessEqual(float(spec_t.kappa.max()), KAPPA_MAX)
    out = {}
    for pp in (1, -1):
        vals = np.linalg.eigvals(blocks[pp])
        test.assertTrue(np.all(np.isfinite(vals)), pp)
        test.assertEqual(len(vals), len(spec_t.vals), pp)
        out["ratio_pp%+d" % pp] = match_spectrum(spec_t, vals)[1]
    ref = dense_spectrum(scipy.linalg.block_diag(Ks, Kt))
    mu = np.linalg.eigvals(blocks[0])
    test.assertTrue(np.all(np.isfinite(mu)) and np.all(np.isfinite(ref.vals)))
    test.assertEqual(len(mu), len(ref.vals))
    cols, _ = match_spectrum(ref, mu)
    out["d"] = float(np.abs(mu - ref.vals[cols]).max()) / ref.scale
    out["lead"] = float(np.concatenate([np.linalg.eigvals(b) for b in blocks.values()]).real.max())
    return out


def nonsu2_inputs(test, fixture):
    """``(p, K, topo)``: the prepared spin-orbital run of ``fixture`` (bridge
    configuration, automatic layout) and the non-spin-orbital singlet and
    triplet bond kernels on the same files (:func:`nonso_kernels`)."""
    with _quiet_logs():
        p, geom_info, interactions = _su2_inputs(fixture, "auto")
        K, topo = nonso_kernels(test, p, geom_info, interactions)
    test.assertNotEqual(bond_channels.non_su2_bonds(topo), [])
    return p, K, topo


def nonso_reduced(p, K):
    """``(K_s, K_t)``: the non-spin-orbital singlet kernel on the even pairs
    and the triplet kernel on the odd pairs."""
    Ep, Em, _ = parity_bases(p.norb_phys, p.shape)
    return _projected(K["singlet"], Ep), _projected(K["triplet"], Em)


def nonsu2_gate(test, p, statics, Ks, Kt):
    """The three checks of the S_z-conserving non-SU(2) gate on the dense
    and the S_z-packed operator ``statics`` of one run; returns the
    measures."""
    from tests.so_bond_fixtures import KAPPA_MAX, dense_spectrum, match_spectrum
    test.assertEqual(_layout_of(p), "sz_packed")
    test.assertEqual(statics["dense"].layout, "dense")
    test.assertEqual(statics["sz_packed"].layout, "sz_packed")
    test.assertLessEqual(statics["dense"].r_V, 1e-12)
    test.assertEqual(statics["dense"].r_G, 0.0)
    out = {"r_V": statics["dense"].r_V}
    full = {}
    for layout, static in statics.items():
        blocks = pair_sz_blocks(test, static, p.norb_phys)
        m = nonsu2_measures(test, blocks, Ks, Kt)
        out[layout] = m
        test.assertLessEqual(max(m["ratio_pp+1"], m["ratio_pp-1"]), 1.0, (layout, m))
        test.assertGreater(m["d"], 1e-3, (layout, m))
        full[layout] = np.concatenate([np.linalg.eigvals(b) for b in blocks.values()])
    # the whole listing of the two layouts, one to one
    spec = dense_spectrum(so_reduced_matrices(statics["dense"])[0])
    test.assertLessEqual(float(spec.kappa.max()), KAPPA_MAX)
    out["layouts"] = match_spectrum(spec, full["sz_packed"])[1]
    test.assertLessEqual(out["layouts"], 1.0, out)
    return out


def _dense_cfg(p):
    """``p`` with the dense layout requested (same band, files and config)."""
    return p._replace(cfg=p.cfg._replace(layout_request="dense"))


def nonsu2_production(test, p):
    """The dense and the S_z-packed production operators of ``p`` on one
    Green function (each layout its own admission, bubble, vertex, G2 and
    kernel)."""
    with _quiet_logs():
        so_green = _green_of(p)
        statics = {}
        for layout, q in (("dense", _dense_cfg(p)), ("sz_packed", p)):
            statics[layout] = _operator(q, so_green, _admit(q))
    return statics


def independent_chain(p, layout):
    """The static operator of ``p`` in ``layout`` built from scratch, without
    the production operator builder: its own Green function
    (``_so_bond_green``), its own plan (``plan_so_bond`` with the static
    rows of the run's ``solver_mode``, ``num_eigenvalues``, ``max_iter`` and
    ``spectral_shift``; the ``select_layout`` reason for the packed
    layout), vertex, bubble into its own sink, dressing at ``l0 = nmat //
    2``, pairing-vertex parts (``V = V_bare + V_fluct`` in new arrays),
    ``G2`` and kernel. Returns ``(static, arrays)`` with the bubble, chi and
    V in ``arrays``."""
    from scipy.sparse.linalg import LinearOperator
    from tests.so_bond_fixtures import NOCAP
    compiled, topo, shape = p.so_run.context.compiled, p.so_run.bond.topo, p.shape
    nd, nvol = 2 * p.norb_phys, int(np.prod(shape))
    so_green = _green_of(p)
    carrier = so_green.carrier
    geom = so_bond.slot_geometry(topo, shape)
    reason = None
    if layout == "sz_packed":
        chosen, reason = so_bond.select_layout(compiled, p.so_run.band.transfer, override="auto")
        assert chosen == "sz_packed", chosen
    plan = so_bond.plan_so_bond(
        nmat=p.nmat, nvol=nvol, B=geom.B, nd=nd, consumers=("pairing_static",), caps=NOCAP,
        layout=layout, layout_reason=reason,
        extra_rows=so_bond.static_so_bond_rows(
            nmat=p.nmat, nvol=nvol, B=geom.B, nd=nd, layout=layout,
            solver_mode=p.solver_mode, num_eigenvalues=p.num_eigenvalues,
            max_iter=p.max_iter, spatial_shape=shape, spectral_shift=p.spectral_shift))
    gamma = so_bond.so_bond_vertex(compiled, topo, shape, plan=plan)
    l0 = p.nmat // 2
    if layout == "sz_packed":
        sm = so_bond.sz_sector_map(geom, nd)
        X = so_bond.so_bond_bubble(carrier.deflated_kw, carrier.green0_tail, p.beta, geom,
                                   shape, plan=plan, freq="static",
                                   sink=so_bond.PackedArrays.zeros(sm, (nvol,)))
        cb = so_bond.so_bond_dress(so_bond.PackedArrays(sm, tuple(b[None] for b in X.blocks)),
                                   gamma, plan=plan, l0=l0, spatial_shape=shape)
        chi = so_bond.PackedArrays(sm, tuple(b[0] for b in cb.blocks))
        vb, vf = so_bond.so_bond_pairing_vertex_parts(gamma, chi, geom, plan=plan)
        V = so_bond.PackedArrays(sm, tuple(a + b for a, b in zip(vb.blocks, vf.blocks)))
    else:
        ND = geom.B * nd * nd
        X = so_bond.so_bond_bubble(carrier.deflated_kw, carrier.green0_tail, p.beta, geom,
                                   shape, plan=plan, freq="static",
                                   sink=np.zeros((nvol, ND, ND), dtype=complex))
        chi = so_bond.so_bond_dress(X[None], gamma, plan=plan, l0=l0, spatial_shape=shape)[0]
        vb, vf = so_bond.so_bond_pairing_vertex_parts(gamma, chi, geom, plan=plan)
        V = vb + vf
    G2 = bond_channels._g2_from_green(carrier.full_sc, p.beta, tail=p.cfg.g2_tail)[..., None]
    if layout == "sz_packed":
        kernel = so_bond.PackedBondPairKernel(
            so_bond.PackedArrays(sm, tuple(b[None] for b in V.blocks)), G2, geom, sm,
            plan=plan, spatial_shape=shape, beta=p.beta, nfreq=1, host_cap=None)
        static = sc._SOBondStatic(layout, reason, geom, plan, None, None, kernel, 0.0,
                                  float(so_bond.green_leak(carrier.full_sc, p.norb_phys,
                                                           axes=(0, 1))))
    else:
        vertex = eliashberg_bond.PairVertexUniform(
            eliashberg_bond.ArrayBlockSource({"V": V[None]}, nd * nd), "V", geom.B, nd * nd,
            1, nvol)
        kernel = eliashberg_bond.BondPairKernel(
            vertex, G2, geom, xp=np, spatial_shape=shape, norb=nd, beta=p.beta, nfreq=1,
            V_inst=None, residency="host", host_cap=np.inf)
        n = nd * nd * nvol
        static = sc._SOBondStatic(
            layout, reason, geom, plan, LinearOperator((n, n), matvec=kernel.matvec,
                                                       dtype=complex),
            sc._SOAntisymmetricBasis((nd, nd) + shape), None,
            float(so_pairing.bond_sz_conservation_residual(V, nd=nd, B=geom.B)),
            float(so_bond.green_leak(carrier.full_sc, p.norb_phys, axes=(0, 1))))
    return static, {"bubble": X, "chi": chi, "V": V, "green": carrier.full_sc}


def _parts(a):
    return a.blocks if isinstance(a, so_bond.PackedArrays) else (a,)


def shares_any(a, b):
    """Whether any block of ``a`` shares memory with any block of ``b``."""
    return any(np.shares_memory(x, y) for x in _parts(a) for y in _parts(b))


class TestNonSU2Bonds(unittest.TestCase):
    """Gate A2-nonSU2: an S_z-conserving but not spin-rotation invariant bond
    (Ising, Hund or Exchange alone on the nearest-neighbour bonds, with the
    on-site CoulombIntra), in both layouts. Premises, asserted: the hopping
    is spin independent, there is no PairLift and no field (h = 0), the
    non-spin-orbital resolver lists the bonds as not SU(2) symmetric, and
    the dense pairing vertex conserves pair S_z (``r_V <= 1e-12``), so the
    kernel is block diagonal in pair S_z and the automatic layout is the
    S_z-packed one. The equal-spin pairs (pair S_z = +-1) couple only
    through the longitudinal block, where the non-spin-orbital spin and
    charge vertices are exact: each of the +1 and -1 sectors has the
    spectrum of the non-spin-orbital triplet kernel on the odd pairs
    (matching rule). The pair S_z = 0 sector does not have the SU(2)
    pattern: its spectrum differs from that of the singlet and the triplet
    kernel together by more than 1e-3. The two layouts' whole spectra agree
    by the matching rule. For Hund and Exchange the two layouts share the
    band and the Green function (each builds its own bubble, vertex, G2 and
    kernel); for Ising the two chains are built from scratch apart, from
    their own prepared runs, and share no array.

    Blind spots and the gates that cover them: the transverse (spin-flip
    particle-hole) halves of Exchange and CoulombIntra never reach the
    equal-spin pairing elements, so the pair S_z = +-1 match is exactly
    blind to them (the SU(2) gate catches a flipped Exchange crossing);
    Gamma_B read at -q keeps the equal-spin spectra of these
    inversion-symmetric fixtures (the SU(2) gate's non-vacuity test); none
    of the gates runs past the RPA instability (every fixture is below
    it)."""

    def test_each_type_both_layouts(self):
        for typ, v in NONSU2_TYPES[1:]:
            with self.subTest(type=typ):
                p, K, _ = nonsu2_inputs(self, nonsu2_fixture(typ, v))
                Ks, Kt = nonso_reduced(p, K)
                nonsu2_gate(self, p, nonsu2_production(self, p), Ks, Kt)

    def test_layouts_built_independently(self):
        typ, v = NONSU2_TYPES[0]
        fx = nonsu2_fixture(typ, v)
        p, K, _ = nonsu2_inputs(self, fx)
        Ks, Kt = nonso_reduced(p, K)
        with _quiet_logs():
            p_dense = _su2_inputs(fx, "dense")[0]
            p_packed = _su2_inputs(fx, "auto")[0]
        self.assertEqual(p_dense.cfg.layout_request, "dense")
        chains = {"dense": independent_chain(p_dense, "dense"),
                  "sz_packed": independent_chain(p_packed, "sz_packed")}
        arrays = {k: c[1] for k, c in chains.items()}
        for name in ("bubble", "chi", "V", "green"):
            self.assertFalse(shares_any(arrays["dense"][name], arrays["sz_packed"][name]), name)
        statics = {k: c[0] for k, c in chains.items()}
        nonsu2_gate(self, p_packed, statics, Ks, Kt)
        # each chain is the production operator of its own run
        prod = nonsu2_production(self, p)
        for layout, static in statics.items():
            for a, b in zip(so_reduced_matrices(static), so_reduced_matrices(prod[layout])):
                self.assertLessEqual(rel(a, b), 1e-12, layout)

    def test_non_vacuity(self):
        """The bare bubble in place of the dressed susceptibility in V_fluct
        (chi_B = chibar_B) breaks the equal-spin match in both layouts."""
        typ, v = NONSU2_TYPES[0]
        p, K, _ = nonsu2_inputs(self, nonsu2_fixture(typ, v))
        Ks, Kt = nonso_reduced(p, K)
        calls = []

        def undressed(X, gamma_b, **kw):
            calls.append(kw["l0"])
            return X
        with mock.patch.object(so_bond, "so_bond_dress", undressed):
            statics = nonsu2_production(self, p)
        self.assertEqual(calls, [p.nmat // 2] * 2)
        for layout, static in statics.items():
            with self.subTest(layout=layout):
                m = nonsu2_measures(self, pair_sz_blocks(self, static, p.norb_phys), Ks, Kt)
                self.assertGreater(min(m["ratio_pp+1"], m["ratio_pp-1"]), 1.0, m)


# ------------------------------------------------------------ refusals and warnings
#: the entry stages of a static spin-orbital bond run, in call order: the
#: bond options, the input files, the compiled interaction, the bond
#: topology, the memory admission and the Green function (the first bond
#: array)
_ENTRY_STAGES = (
    (sc, "_read_so_bond_config"), (sc, "_read_interaction_files"),
    (spin_orbital, "build_so_entry_context"), (sc, "_so_bond_run"),
    (sc, "_admit_so_bond_static"), (sc, "_so_bond_green"))
_STAGE_NAMES = tuple(name for _, name in _ENTRY_STAGES)
#: the non-spin-orbital bond readers and guards, and the backend selection
#: (which a bond run skips)
_NON_SO_STAGES = ((sc, "_read_bond_config"), (sc, "_validate_bond_prereqs"),
                  (sc, "_check_nonbond_offsite"), (backend, "get_backend"))
#: fragments of the non-spin-orbital bond refusals
_NON_SO_FRAGMENTS = ("bond_channels=true is incompatible", "bond_channels=true requires",
                     "must be a positive finite number")
_SO_CONFIG = "[mode] enable_spin_orbital = true (hwave_sc, issue #83): [eliashberg] "
_SAME_OPERATOR = (
    "{t} coefficient at R=(-1, 0, 0) orbitals (1, 1) has an imaginary part (-0.1). In "
    "H-wave's declaration convention the entry and its partner (−R,b,a) multiply the same "
    "operator, so the imaginary part has no operator content; only PairHop carries a "
    "physical complex phase (complex couplings beyond PairHop: issue #302). Remove the "
    "imaginary part.")
_COMPLEX_EXCHANGE = (
    "Exchange coefficient at R=(-1, 0, 0) orbitals (1, 1) has an imaginary part (-0.1). The "
    "partner entry (−R,b,a) multiplies the Hermitian conjugate, so an imaginary part would "
    "be a Dzyaloshinskii-Moriya-type term; H-wave's declaration convention folds it, and "
    "that extension is tracked in issue #302. Remove the imaginary part.")
#: the first paragraph of the over-cap refusal of the 4 x 4, B = 5, Nmat 32
#: fixtures at a cap of 1e-6 GiB
_OVER_CAP_HEAD = (
    "spin-orbital bond Eliashberg needs an estimated peak of 0.00 GiB in phase "
    "pairing_static (B = 5 channels, ND = B·(2·norb_phys)² = 20, nvol = 16, Nmat = 32; "
    "largest row in that phase: dress_batch = 0.00 GiB) but the cap is 0.00 GiB ([eliashberg] "
    "bond_memory_cap_gb). Memory grows as Nmat · nvol · norb_phys⁴ (the bubble's per-shift "
    "workspace) and Nmat · nvol · norb_phys² (the Green function); independent of Nmat, as "
    "nvol · B² · norb_phys⁴ (the vertex arrays and the kernel workspace); the eigensolver rows "
    "grow with the reduced gap dimension dim (about nvol · (2·norb_phys)² / 2): linearly times "
    "num_eigenvalues for Arnoldi, as dim² for the dense fallback (dim <= 2000), with max_iter "
    "for the power iteration, and they depend on the layout and spectral_shift. Reduce Nmat "
    "or the k grid, or raise [eliashberg] bond_memory_cap_gb (0 disables the check); pruning "
    "bond shells helps only for shells the per-shell table below shows with negligible "
    "coupling.")
#: the table headers of that refusal, in order (the dense layout marked as
#: the run's)
_OVER_CAP_TABLES = (
    "rows (phase  name  bytes):",
    "admitted cases (machine | path | admitted):",
    "  S_z-packed layout (S_z-conserving input: Zeeman, no SOC, no PairLift)",
    "  dense layout (S_z broken by SOC or PairLift)  <- this run",
    "shells (shell  |dr|  max|v|):")


def _complex_row_fixture(typ):
    """SOC_ONSITE with a complex nearest-neighbour ``typ`` row along x (its
    Hermitian partner at -R)."""
    return F.SOC_ONSITE._replace(
        name="SOC_COMPLEX_" + typ.upper(),
        interactions=F.SOC_ONSITE.interactions + (
            (typ, (((1, 0, 0), (0, 0), 0.3 + 0.1j), ((-1, 0, 0), (0, 0), 0.3 - 0.1j))),))


def _refusal_rows():
    """``(label, fixture, eliashberg, text, first stage not entered)`` of
    every refusal of the static spin-orbital bond path through the public
    entry; ``text`` ``None``: the function-level admission's refusal."""
    from tests.test_sc_spin_orbital import OFFSITE_TEXT
    rows = []
    for typ in ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift"):
        text = _COMPLEX_EXCHANGE if typ == "Exchange" else _SAME_OPERATOR.format(t=typ)
        rows.append(("complex " + typ, _complex_row_fixture(typ), None, text, "_so_bond_run"))
    rows += [
        ("layout sz_packed", F.SOC_ALL_TYPES, {"SO_Bond_Layout": "sz_packed"},
         '[eliashberg] so_bond_layout = "sz_packed" cannot be forced: the S_z-packed layout '
         'is chosen automatically when the input conserves S_z (no spin-flip Transfer entry, '
         'no S_z-changing interaction). Use "auto" (the default) or "dense".',
         "_read_interaction_files"),
        ("layout x", F.SOC_ALL_TYPES, {"So_Bond_Layout": "x"},
         "[eliashberg] so_bond_layout must be \"auto\" or \"dense\", got 'x'",
         "_read_interaction_files"),
    ]
    for raw in (-1, float("nan")):
        rows.append(("cap %r" % raw, F.SOC_ALL_TYPES, {"Bond_Memory_Cap_GB": raw},
                     "[eliashberg] bond_memory_cap_gb must be a finite number >= 0 in "
                     "spin-orbital mode (0 disables the check), got %r" % raw,
                     "_read_interaction_files"))
    rows += [
        ("over cap", F.SOC_ALL_TYPES, {"Bond_Memory_Cap_GB": 1e-6}, None, "_so_bond_green"),
        ("over cap, declared-zero outer shell", _SOC_DECLARED_ZERO,
         {"Bond_Memory_Cap_GB": 1e-6}, None, "_so_bond_green"),
        ("max_shells drops a used shell", F.SOC_ALL_TYPES, {"Bond_Max_Shells": 0},
         "so_bond_topology: max_shells=0 would drop bond shell(s) used by nonzero "
         "interaction records: shell 1 with displacements [(-1, 0, 0), (0, -1, 0), "
         "(0, 1, 0), (1, 0, 0)]. Use max_shells >= 1, or remove those declarations.",
         "_admit_so_bond_static"),
        ("dynamic", F.SOC_ALL_TYPES, {"Frequency": "dynamic"},
         so_pairing._SO_DYNAMIC_BOND_REFUSAL, "_read_so_bond_config"),
        ("chi0q_mode flex", F.SOC_ALL_TYPES, {"Chi0q_Mode": "flex"}, _SO_BOND_FLEX_REFUSAL,
         "_read_interaction_files"),
        ("chi0q_mode flex + flex_chi_s", F.SOC_ALL_TYPES,
         {"Chi0q_Mode": "flex", "Flex_Chi_S": "chi_s.npz"}, _SO_BOND_FLEX_REFUSAL,
         "_read_interaction_files"),
    ]
    for value in ("FLEX", "Flex"):
        rows.append(("chi0q_mode %s" % value, F.SOC_ALL_TYPES, {"Chi0q_Mode": value},
                     _SO_CONFIG + "chi0q_mode = {!r} is not supported; use \"calc\" (the "
                     "spinful chi0 is computed here), \"load\" (a general chi0q of a "
                     "spin-orbital RPA run) or \"flex\" (the chiq and green of a "
                     "spin-orbital FLEX run).".format(value), "_read_so_bond_config"))
    for value in (True, False):
        rows.append(("bond_diagnostics %s" % value, F.SOC_ALL_TYPES,
                     {"Bond_Diagnostics": value}, _SO_BOND_DIAGNOSTICS_REFUSAL,
                     "_read_interaction_files"))
    rows += [
        ("bond_precondition_atol", F.SOC_ALL_TYPES, {"Bond_Precondition_ATOL": 1e-12},
         _SO_BOND_PRECONDITION_REFUSAL.format(key="Bond_Precondition_ATOL"),
         "_read_interaction_files"),
        ("bond_green", F.SOC_ALL_TYPES, {"Bond_Green": "green.npz"}, _SO_BOND_GREEN_REFUSAL,
         "_read_interaction_files"),
        ("gpu_required", F.SOC_ALL_TYPES, {"GPU_Required": True, "Gpu": True},
         _SO_BOND_GPU_REQUIRED_REFUSAL, "_read_interaction_files"),
        ("off-site rows, bond_channels = false", F.SOC_ALL_TYPES, {"Bond_Channels": False},
         OFFSITE_TEXT.format(R=(1, 0, 0), orbs=(1, 1), src="the CoulombInter file"),
         "build_so_entry_context"),
        ("eigenvalue_method lanczos", F.SOC_ALL_TYPES, {"Eigenvalue_Method": "lanczos"},
         _SO_CONFIG + "eigenvalue_method = 'lanczos' is not supported; only \"arnoldi\" is. "
         "Set eigenvalue_method = \"arnoldi\" (or remove the key).", "_read_so_bond_config"),
    ]
    return rows


def _admission_refusal(fixture, eliashberg):
    """The refusal text of the function-level admission of ``fixture`` (the
    production readers up to the band, then ``_admit_so_bond_static``)."""
    d = tempfile.mkdtemp(prefix="sc_so_bond_refusal_")
    try:
        with _quiet_logs():
            p = F.prepare(fixture, d, eliashberg=eliashberg)
            sc._admit_so_bond_static(
                p.so_run, cfg=p.cfg, cell=p.cell, nmat=p.nmat, nd=fixture.nd, shape=p.shape,
                solver_mode=p.solver_mode, num_eigenvalues=p.num_eigenvalues,
                max_iter=p.max_iter, spectral_shift=p.spectral_shift,
                eigenvalue_method="arnoldi", log=sc.logger)
    except ValueError as err:
        return str(err)
    finally:
        shutil.rmtree(d, True)
    raise AssertionError("the admission of {} was not refused".format(fixture.name))


class TestRefusalMatrix(unittest.TestCase):
    """Gate A2-R (refusals) of the static spin-orbital bond path through the
    public entry (``hwave.sc.calc_eliashberg``), keys in mixed case: every
    row is refused with its exact text, before any bond array is allocated
    (the entry stages after the refusing one are never entered, the Green
    function never), no output is written, and no row reaches a
    non-spin-orbital reader, guard or the backend selection, nor any of
    their texts."""

    def _spied(self, fixture, eliashberg, *, refused=True):
        """``(refusal text or None, entered stages, non-SO calls, vertex
        calls, output directory)`` of a run of ``fixture`` through the entry
        with every stage spied."""
        d = tempfile.mkdtemp(prefix="sc_so_bond_refusal_")
        self.addCleanup(shutil.rmtree, d, True)
        inp = F.write_case(d, fixture, eliashberg=eliashberg)
        text = None
        with contextlib.ExitStack() as stack:
            stack.enter_context(_quiet_logs())
            stages = {name: stack.enter_context(mock.patch.object(
                mod, name, wraps=getattr(mod, name))) for mod, name in _ENTRY_STAGES}
            non_so = {name: stack.enter_context(mock.patch.object(
                mod, name, wraps=getattr(mod, name))) for mod, name in _NON_SO_STAGES}
            vertex = stack.enter_context(mock.patch.object(
                so_bond, "so_bond_vertex", wraps=so_bond.so_bond_vertex))
            if refused:
                with self.assertRaises(ValueError) as cm:
                    sc.calc_eliashberg(inp)
                text = str(cm.exception)
            else:
                sc.calc_eliashberg(inp)
        entered = {name: spy.call_count for name, spy in stages.items()}
        calls = {name: spy.call_count for name, spy in non_so.items()}
        return text, entered, calls, vertex.call_count, os.path.join(d, "out")

    def _refuse(self, fixture, eliashberg):
        text, entered, non_so, vertex, out = self._spied(fixture, eliashberg)
        self.assertFalse(os.path.exists(out))
        # an on-site (bond-off) run selects its backend; a bond run never does
        bond_off = any(k.lower() == "bond_channels" and v is False
                       for k, v in (eliashberg or {}).items())
        for name, n in non_so.items():
            self.assertEqual(n, int(bond_off and name == "get_backend"), name)
        self.assertEqual(vertex, 0)
        return text, entered

    def test_spies_see_a_valid_run(self):
        """The control of the spies: a valid bond run through the same spy
        set enters every stage (the Green function and the vertex included)
        and no non-spin-orbital reader, guard or backend selection."""
        text, entered, non_so, vertex, out = self._spied(F.ZEEMAN_EXTERN, None,
                                                         refused=False)
        self.assertIsNone(text)
        self.assertTrue(all(n == 1 for n in entered.values()), entered)
        self.assertEqual(vertex, 1)
        self.assertEqual(set(non_so.values()), {0}, non_so)
        self.assertTrue(os.path.exists(os.path.join(out, "gap.npz")))

    def test_flex_refusal_text(self):
        """The chi0q_mode = "flex" refusal, written out (the rows compare
        with the constant)."""
        self.assertEqual(_SO_BOND_FLEX_REFUSAL, (
            _SO_CONFIG + "chi0q_mode = \"flex\" with bond_channels = true is not available "
            "yet: the spin-orbital bond path builds its bond bubble from the band's Green "
            "function and reads no FLEX output. Remove chi0q_mode (or leave it unset); it is "
            "not used on this path."))

    def test_rows(self):
        rows = _refusal_rows()
        self.assertEqual(len(rows), 24)
        for label, fixture, eli, text, first_skipped in rows:
            with self.subTest(row=label):
                got, entered = self._refuse(fixture, eli)
                if text is None:
                    text = _admission_refusal(fixture, eli)
                self.assertEqual(got, text)
                cut = _STAGE_NAMES.index(first_skipped)
                self.assertEqual({n: entered[n] > 0 for n in _STAGE_NAMES},
                                 {n: i < cut for i, n in enumerate(_STAGE_NAMES)}, entered)
                for fragment in _NON_SO_FRAGMENTS:
                    self.assertNotIn(fragment, got)

    def test_over_cap_texts(self):
        """The over-cap refusal (both fixtures of the matrix): its first
        paragraph, and the per-shell table with, for the declared-zero
        nearest-neighbour shell, the clause naming the cutoff that removes
        it."""
        for fixture, tail in (
                (F.SOC_ALL_TYPES, "shells (shell  |dr|  max|v|):\n  1  1  0.31"),
                (_SOC_DECLARED_ZERO,
                 "shells (shell  |dr|  max|v|):\n  1  1  0\n1 outer shell(s) declared only "
                 "with zero coupling; [eliashberg] bond_max_shells = 0 removes them")):
            with self.subTest(fixture=fixture.name):
                text, _ = self._refuse(fixture, {"Bond_Memory_Cap_GB": 1e-6})
                self.assertEqual(text.split("\n", 1)[0], _OVER_CAP_HEAD)
                lines = text.split("\n")
                at = [lines.index(h) for h in _OVER_CAP_TABLES]
                self.assertEqual(at, sorted(at))
                self.assertEqual(sum(ln.endswith("  <- this run") for ln in lines), 1)
                self.assertNotIn("store slot", text)
                for path in ("FLEX", "dynamic", "A100"):
                    self.assertNotIn(path, text)
                self.assertTrue(text.endswith("\n" + tail), text)
                self.assertEqual(text.count("declared only with zero coupling"),
                                 int(fixture is _SOC_DECLARED_ZERO))


#: a complex nearest-neighbour PairHop row along x and y (its Hermitian
#: partner at -R): accepted, S_z-conserving
_COMPLEX_PAIRHOP = (("PairHop", (((1, 0, 0), (0, 0), 0.2 + 0.15j),
                                 ((-1, 0, 0), (0, 0), 0.2 - 0.15j),
                                 ((0, 1, 0), (0, 0), 0.2 + 0.15j),
                                 ((0, -1, 0), (0, 0), 0.2 - 0.15j))),)
#: the Zeeman fixture (S_z-packed) with the complex PairHop row
ZEEMAN_PAIRHOP = F.ZEEMAN_EXTERN._replace(
    name="ZEEMAN_COMPLEX_PAIRHOP", interactions=F.ZEEMAN_EXTERN.interactions + _COMPLEX_PAIRHOP)
#: the warning rows: ``(label, param, eliashberg, WARNING text)``
WARNING_ROWS = (
    ("orbital_site", {"Orbital_Site": [1]}, {}, spin_orbital.ORBITAL_SITE_WARNING),
    ("gpu", {}, {"Gpu": True}, _SO_BOND_GPU_IGNORED),
    ("chi0q_tensor", {}, {"Chi0q_Tensor": "reduced"},
     _SO_BOND_CHI0Q_UNUSED % ("unset", "reduced")),
    ("flex_chi_c", {}, {"Flex_Chi_C": "x.npz"}, _SO_BOND_FLEX_KEYS_UNUSED % ("Flex_Chi_C",)),
    ("sigma_shift (packed)", {}, {"Sigma_Shift": 0.1}, sc._SIGMA_SHIFT_IGNORED),
)
#: the cheap state point of the warning runs (a fixed cap: no cap WARNING)
_WARN_PARAM = {"Nmat": 8}
_WARN_ELI = {"Bond_Memory_Cap_GB": 2.0}


def warning_messages(run):
    return [r.getMessage() for r in run.records if r.levelno == logging.WARNING]


def check_completed(test, run, layout="sz_packed"):
    """The run completed: a finite leading eigenvalue, both outputs, the
    expected layout."""
    test.assertTrue(np.isfinite(run.value))
    test.assertIsNotNone(run.eigenvalue_dat)
    test.assertEqual(str(run.gap_npz["block_layout"]), layout)


class TestWarningMatrix(unittest.TestCase):
    """Gate A2-R (warnings) of the static spin-orbital bond path through the
    public entry: one run carrying every accepted-but-unused option (keys in
    mixed case) and a complex PairHop row completes, and each WARNING
    appears exactly once (the RPA-ladder WARNING of every bond run
    included, nothing else); each row also completes alone with its own
    WARNING once, the complex PairHop row with none; and ``so_bond_layout``
    with ``bond_channels = false`` is named as ignored once on the on-site
    path."""

    def test_one_completing_run(self):
        param, eli = dict(_WARN_PARAM), dict(_WARN_ELI)
        for _, p, e, _ in WARNING_ROWS:
            param.update(p)
            eli.update(e)
        run = F.run_case(ZEEMAN_PAIRHOP, fresh=True, param=param, eliashberg=eli)
        check_completed(self, run)
        want = [text for _, _, _, text in WARNING_ROWS] + [sc._SO_BOND_RPA_LADDER]
        self.assertEqual(sorted(warning_messages(run)), sorted(want))

    def test_each_warning_row(self):
        rows = WARNING_ROWS + (("complex PairHop", {}, {}, None),)
        for label, param, eli, text in rows:
            with self.subTest(row=label):
                fx = ZEEMAN_PAIRHOP if text is None else F.ZEEMAN_EXTERN
                run = F.run_case(fx, fresh=True, param=dict(_WARN_PARAM, **param),
                                 eliashberg=dict(_WARN_ELI, **eli))
                check_completed(self, run)
                want = [sc._SO_BOND_RPA_LADDER] + ([] if text is None else [text])
                self.assertEqual(sorted(warning_messages(run)), sorted(want))

    def test_layout_key_with_bond_off(self):
        run = F.run_case(F.SOC_ONSITE, fresh=True, param=_WARN_PARAM,
                         eliashberg={"Bond_Channels": False, "Chi0q_Mode": "calc",
                                     "SO_Bond_Layout": "dense"})
        self.assertTrue(np.isfinite(run.value))
        self.assertIsNotNone(run.eigenvalue_dat)
        self.assertNotIn("block_layout", run.gap_npz)
        warned = warning_messages(run)
        self.assertEqual(warned.count(sc._SO_BOND_LAYOUT_IGNORED), 1, warned)
        self.assertNotIn(sc._SO_BOND_RPA_LADDER, warned)

    def test_sigma_shift_text_shared(self):
        """One sigma_shift WARNING text for every eigenvalue solve that takes
        sigma_shift: the ARPACK path, the dense spin-orbital branch and the
        S_z-packed sector solve (the row above)."""
        from scipy.sparse.linalg import aslinearoperator
        self.assertEqual(sc._SIGMA_SHIFT_IGNORED, (
            "[eliashberg] sigma_shift is ignored for eigenvalue_method='arnoldi' (it targets "
            "the shift-invert methods); use spectral_shift to bias the arnoldi selection."))
        n = 12
        A = np.diag(np.arange(1.0, n + 1)).astype(complex)
        op = aslinearoperator(A)
        basis = sc._SOAntisymmetricBasis((2, 2, 3, 1, 1))
        calls = (
            ("ARPACK", lambda: sc._solve_leading(lambda: (op, n), n, "arnoldi",
                                                 num_eigenvalues=2, sigma_shift=0.1)),
            ("dense spin-orbital", lambda: sc._so_solve_eigenvalue(
                aslinearoperator(np.eye(basis.n, dtype=complex)), basis, basis.dim, "arnoldi",
                sigma_shift=0.1)))
        for name, call in calls:
            with self.subTest(solve=name), self.assertLogs(sc.logger, "WARNING") as cm:
                call()
            self.assertIn(sc._SIGMA_SHIFT_IGNORED, [r.getMessage() for r in cm.records])


if __name__ == "__main__":
    unittest.main()
