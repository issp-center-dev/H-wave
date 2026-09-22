"""Memory preflight of the bond-resolved FLEX gate (spec 2026-09-06
section 3.6): the named-buffer table, the batch selection, the refusal
naming every row, the operation warnings, the tracer invariant and the
subprocess RSS validation."""
import json
import logging
import os
import subprocess
import sys
import tempfile
import tracemalloc
import unittest

import numpy as np

from tests.heavy_tests import heavy

_IN2 = "tests/rpa/input_2orb"


def _kw(**over):
    kw = dict(nmat=8, nvol=16, norb=2, B=3, depth=0, output_full=False, split_seed=False,
              n_types=1, freq_batch=None, cap_gb=8.0, mixing="linear")
    kw.update(over)
    return kw


class TestEstimate(unittest.TestCase):

    def test_rows_follow_the_table(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        est = estimate_bond_memory(**_kw(depth=4, mixing="anderson", output_full=True, split_seed=True,
                                         n_types=2))
        nmat, nvol, norb, B = 8, 16, 2, 3
        P = norb ** 2; ND = B * P
        U = nmat * nvol * ND * ND * 16; G = nmat * nvol * P * 16; C = nmat * nvol * P * P * 16
        S = nvol * ND * ND * 16; H = nvol * P * 16
        pr = est["persistent_rows"]
        self.assertEqual(pr["chibar_W"], 2 * U)
        self.assertEqual(pr["chi_sc_w"], 2 * U)
        self.assertEqual(pr["vertices_static"], 5 * S)
        self.assertEqual(pr["state"], G + H)
        self.assertEqual(pr["seed_envelope"], 2 * G + H)
        self.assertEqual(pr["anderson_history"], 2 * 4 * 2 * G)
        self.assertEqual(pr["anderson_work"], 2 * 3 * 2 * G)
        self.assertEqual(pr["collapses"], 3 * C)
        self.assertEqual(pr["eigenpairs_hf"], 5 * H)
        self.assertEqual(pr["flex_arrays"], 5 * G)
        self.assertEqual(pr["hf_tables"], 2 * H)
        # the second-order factor pack: this call passes no factor bytes
        # (the default, i.e. flex_second_order = "takimoto", which compiles none)
        self.assertEqual(pr["second_order_factors"], 0)
        self.assertEqual(est["persistent"], sum(pr.values()))
        ph = est["phase_rows"]
        self.assertEqual(ph["green_mu"], 6 * G + H)
        self.assertEqual(ph["density_hf"], 2 * G + 35 * H)
        prep = 16 * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = 16 * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        self.assertEqual(ph["bubble"], max(prep, pair))
        self.assertEqual(ph["dressing"], 6 * est["nb"] * nvol * ND * ND * 16)
        self.assertEqual(ph["transport"], 4 * G + 4 * C)
        self.assertEqual(ph["convergence"], 7 * G + H)
        self.assertEqual(ph["mixing"], 4 * G)
        self.assertEqual(ph["post_scf"], G + U)
        self.assertEqual(est["nb"], nmat)
        self.assertAlmostEqual(est["peak"], 1.25 * (est["persistent"] + max(ph.values())))

    def test_nb_selection_and_refusal(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        big = estimate_bond_memory(**_kw())
        nvol, ND = 16, 12
        per_batch = 6 * nvol * ND * ND * 16
        others = max(v for k, v in big["phase_rows"].items() if k != "dressing")
        # a cap that admits exactly three batches worth of dressing
        cap = 1.25 * (big["persistent"] + max(others, 3 * per_batch)) / 1024 ** 3
        est = estimate_bond_memory(**_kw(cap_gb=cap * (1 + 1e-9)))
        self.assertEqual(est["nb"], 3 if 3 * per_batch >= others else 8)
        est = estimate_bond_memory(**_kw(cap_gb=cap * (1 + 1e-9), freq_batch=2))
        self.assertEqual(est["nb"], 2)
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(**_kw(cap_gb=cap * (1 + 1e-9), freq_batch=8))
        self.assertIn("longitudinal_bond_freq_batch", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(**_kw(cap_gb=1e-6))
        msg = str(cm.exception)
        for row in ("chibar_W", "green_mu", "dressing", "transport", "post_scf", "GiB"):
            self.assertIn(row, msg)

    def test_operation_formulas(self):
        from hwave.solver.flex_bond import dressing_ops, transport_ops
        self.assertEqual(dressing_ops(8, 16, 12), 2 * 8 * 16 * 12 ** 3)
        P = 4
        self.assertAlmostEqual(transport_ops(3, 8, 16, 2),
                               9 * 8 * 16 * (P ** 2 * np.log2(8 * 16) + P ** 3))


def _solver(param_extra=None, inter=None):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if inter is None else inter)
    r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 8, "IterationMax": 2, "Mix": 0.5, "EPS": 1e-12,
           "flex_hartree_fock": True, "longitudinal_bond_channels": True}
    par.update(param_extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


def _preimport():
    """Import everything the solve imports lazily, so that module objects do
    not count as solve memory."""
    import scipy.fft, scipy.linalg  # noqa: F401
    from hwave.solver import (flex_bond, flex_hf, flex_mixing, hartree_fock, bubble,  # noqa: F401
                              bond_channels, offsite, matsubara, backend, vertex_table)
    from hwave.solver import _sc_matrices_myo  # noqa: F401
    import hwave.sc  # noqa: F401


class _Tracer:
    """``solver._phase_tracer``: per-phase tracemalloc peaks measured from
    the traced memory at SOLVE entry (``base``), so every buffer the solve
    allocates -- persistent ones included, also those a phase allocates
    for the first time such as the Anderson stacks -- counts, while the
    process's pre-solve state does not. The invariant compares each phase
    with persistent + its row."""

    def __init__(self, base=0):
        self.peaks = {}
        self.base = base

    def __call__(self, name):
        tracer = self

        class _Ctx:
            def __enter__(self_):
                tracemalloc.reset_peak()
                self_.start = tracemalloc.get_traced_memory()[0]
                return self_

            def __exit__(self_, *exc):
                peak = tracemalloc.get_traced_memory()[1]
                tracer.peaks[name] = max(tracer.peaks.get(name, 0), peak - tracer.base)
                return False
        return _Ctx()


class TestSolverPreflight(unittest.TestCase):

    def test_operation_warnings_are_logged(self):
        import hwave.solver.flex_bond as fb
        s, r = _solver()
        old = fb._DRESSING_OPS_WARN, fb._TRANSPORT_OPS_WARN
        fb._DRESSING_OPS_WARN, fb._TRANSPORT_OPS_WARN = 1.0, 1.0
        try:
            with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
                s._phase_b_reset({})
                s._phase_b_preflight(r.get_param("green"))
        finally:
            fb._DRESSING_OPS_WARN, fb._TRANSPORT_OPS_WARN = old
        text = "\n".join(cm.output)
        self.assertIn("dense", text)
        self.assertIn("transport", text)

    def _run_traced(self, param_extra, green_extra=None):
        s, r = _solver(param_extra)
        gi = r.get_param("green")
        gi.update(green_extra or {})
        _preimport()
        tracer = _Tracer()
        s._phase_tracer = tracer
        tracemalloc.start()
        try:
            with tempfile.TemporaryDirectory() as out:
                tracer.base = tracemalloc.get_traced_memory()[0]
                s.solve(gi, out)
        finally:
            tracemalloc.stop()
        return s, tracer

    def test_preflight_rows_match_tracemalloc_per_phase(self):
        for mixing in ("linear", "anderson"):
            with self.subTest(mixing=mixing):
                s, tracer = self._run_traced({"mixing_scheme": mixing, "anderson_depth": 3,
                                              "longitudinal_bond_output_full": True})
                est = s._bond_est
                self.assertTrue(tracer.peaks, "no phase was traced")
                for name, measured in tracer.peaks.items():
                    self.assertIn(name, est["phase_rows"])
                    bound = est["persistent"] + est["phase_rows"][name]
                    self.assertLessEqual(measured, bound,
                                         "phase {} peak since solve entry {} > persistent {} + row {}".format(
                                             name, measured, est["persistent"], est["phase_rows"][name]))
                self.assertGreaterEqual(len(tracer.peaks), 7)

    def test_tracer_invariant_split_warm_start(self):
        import hwave.solver.flex as flex_mod
        s, r = _solver({"IterationMax": 1})
        gi = r.get_param("green")
        with tempfile.TemporaryDirectory() as out:
            s.solve(gi, out)
            s.save_results({"path_to_output": out, "sigma": "sigma.npz"}, gi)
            s2, r2 = _solver({"IterationMax": 1})
            gi2 = r2.get_param("green")
            gi2.update(s2.read_init({"path_to_input": out, "sigma_init": "sigma.npz"}))
            _preimport()
            tracer = _Tracer()
            s2._phase_tracer = tracer
            tracemalloc.start()
            try:
                tracer.base = tracemalloc.get_traced_memory()[0]
                s2.solve(gi2, out)
            finally:
                tracemalloc.stop()
        est = s2._bond_est
        self.assertEqual(est["persistent_rows"]["seed_envelope"] > 0, True)
        for name, measured in tracer.peaks.items():
            self.assertLessEqual(measured, est["persistent"] + est["phase_rows"][name], name)


_RSS_SCRIPT = r"""
import json, resource, sys, tempfile
import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.flex as flex_mod
solve = sys.argv[1] == "solve"
idict = {"path_to_input": "tests/rpa/input_2orb", "Geometry": "geom.dat", "Transfer": "transfer.dat",
         "CoulombInter": "coulombinter.dat"}
r = read_input_k.QLMSkInput({"path_to_input": "tests/rpa/input_2orb", "interaction": idict})
par = {"T": 2.0, "filling": 0.5, "CellShape": [6, 6, 1], "SubShape": [1, 1, 1], "Nmat": 64,
       "IterationMax": 3, "Mix": 0.5, "EPS": 1e-12, "flex_hartree_fock": True,
       "longitudinal_bond_channels": True, "longitudinal_bond_output_full": True,
       "mixing_scheme": "anderson", "anderson_depth": 8}
info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
s = flex_mod.FLEX(r.get_param("ham"), {}, info)
gi = r.get_param("green")
peak = None
if solve:
    with tempfile.TemporaryDirectory() as out:
        s.solve(gi, out)
        s.save_results({"path_to_output": out, "sigma": "sigma.npz", "green": "green.npz"}, gi)
    peak = s._bond_est["peak"]
rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
if sys.platform != "darwin":
    rss *= 1024
print(json.dumps({"rss": rss, "peak": peak}))
"""


class TestSubprocessRSS(unittest.TestCase):
    #: incremental regression bound (gated peak RSS minus the calibrated
    #: baseline), calibrated 2026-09-06 on the 6x6, Nmat=64, depth-8,
    #: output_full case (measured: baseline 80 MiB, gated 238 MiB, estimate
    #: 186 MiB, increment 158 MiB); 20 % slack is applied in the assertion.
    INCREMENT_BOUND_BYTES = 165 * 1024 ** 2

    @heavy
    def test_subprocess_rss_below_estimate_plus_allowance(self):
        env = dict(os.environ, PYTHONPATH="src:.")

        def _run(kind):
            out = subprocess.run([sys.executable, "-B", "-c", _RSS_SCRIPT, kind], env=env,
                                 capture_output=True, text=True, check=True)
            line = [ln for ln in out.stdout.splitlines() if ln.startswith("{")][-1]
            return json.loads(line)
        base = _run("baseline")
        gated = _run("solve")
        logging.getLogger("qlms").info("RSS baseline %.1f MiB, gated %.1f MiB, estimate %.1f MiB",
                                       base["rss"] / 2 ** 20, gated["rss"] / 2 ** 20,
                                       gated["peak"] / 2 ** 20)
        self.assertLess(gated["rss"], gated["peak"] + base["rss"])
        self.assertLess(gated["rss"] - base["rss"], 1.2 * self.INCREMENT_BOUND_BYTES)


class TestDeviceMemoryTable(unittest.TestCase):
    _KW = dict(nmat=64, nvol=16, norb=2, B=3, depth=4, output_full=False, split_seed=False,
               n_types=1, freq_batch=None, cap_gb=200.0, mixing="anderson")

    def test_no_device_table_without_available_bytes(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        est = estimate_bond_memory(**self._KW)
        self.assertNotIn("device_rows", est)

    #: the device table's symbols at ``_KW`` (spec 4.6 plus the two
    #: persistent rows that are resident for the whole Phase B solve)
    _IT = 16
    _NVOL = 16
    _P = 4
    _ND = 3 * _P
    _NMAT = 64

    @classmethod
    def _sym(cls):
        it, nvol, P, ND, nmat = cls._IT, cls._NVOL, cls._P, cls._ND, cls._NMAT
        return (nvol * ND * ND * it,                 # S
                nmat * nvol * P * it,                # G
                nmat * nvol * P * P * it,            # C
                7 * nvol * ND * ND * it)             # dressing per batch row

    @classmethod
    def _bubble(cls):
        """The bubble's device temporaries -- the same expression as the
        host row, which was derived for that allocation pattern."""
        it, nvol, P, ND, nmat = cls._IT, cls._NVOL, cls._P, cls._ND, cls._NMAT
        prep = it * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        return max(prep, pair)

    @classmethod
    def _guard(cls, nb):
        """The conditioning guard's device row (issue #197), incremental
        to the dressing row at the same batch (under
        ``guard_freqs = "all"``, the default; factor 4, raised from 3 to
        keep headroom over the measured ~3x peak)."""
        return 4 * nb * cls._NVOL * cls._ND * cls._ND * cls._IT

    @classmethod
    def _need(cls, nb, factor_bytes=0):
        S, G, C, per = cls._sym()
        # 5 G the SCF loop's arrays, one more G the bubble's tail (issue #196)
        return 1.25 * (5 * S + 5 * G + G + factor_bytes
                       + max(per * nb + cls._guard(nb), 6 * C, cls._bubble()))

    def test_device_rows_and_selection(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        S, G, C, per = self._sym()
        est = estimate_bond_memory(device_available=int(2.0e8), **self._KW)
        self.assertEqual(est["device_rows"]["vertices_static"], 5 * S)
        # the SCF loop's device Green-function/sigma arrays outlive every bond
        # phase, so they are a persistent device row like the host's
        self.assertEqual(est["device_rows"]["flex_arrays"], 5 * G)
        self.assertEqual(est["device_rows"]["second_order_factors"], 0)
        self.assertEqual(est["device_rows"]["transport"], 6 * C)
        # the bubble runs on the device too (issue #196), so its temporaries
        # are a device phase row with the host row's expression
        self.assertEqual(est["device_rows"]["bubble"], self._bubble())
        # ... and its tail is a solve-scoped device copy, so a persistent row
        self.assertEqual(est["device_rows"]["green0_tail"], G)
        nb = est["device_nb"]
        self.assertEqual(est["device_rows"]["dressing"], per * nb)
        self.assertEqual(est["device_rows"]["guard"], self._guard(nb))
        self.assertLessEqual(self._need(nb), 0.9 * 2.0e8)
        self.assertTrue(nb == self._NMAT or self._need(nb + 1) > 0.9 * 2.0e8)
        nb_host = estimate_bond_memory(device_available=None, **self._KW)["nb"]
        self.assertEqual(est["nb"], min(nb_host, est["device_nb"]))

    def test_guard_freqs_static_caps_the_row_at_one_slice(self):
        """Under ``guard_freqs = "static"`` only the zero-frequency slice
        (``nvol`` blocks) is ever guarded per batch, so the row does not
        scale with the selected ``nb``."""
        from hwave.solver.flex_bond import estimate_bond_memory
        S, G, C, per = self._sym()
        est = estimate_bond_memory(device_available=int(2.0e8), guard_freqs="static", **self._KW)
        self.assertEqual(est["device_rows"]["guard"], 4 * self._NVOL * self._ND * self._ND * self._IT)

    def test_guard_disabled_allocates_nothing(self):
        """``guard_enabled=False`` (the conditioning floor is off,
        ``cond_tol is None``) allocates no device guard row at all."""
        from hwave.solver.flex_bond import estimate_bond_memory
        est = estimate_bond_memory(device_available=int(2.0e8), guard_enabled=False, **self._KW)
        self.assertEqual(est["device_rows"]["guard"], 0)

    def test_second_order_mirror_tightens_the_device_batch(self):
        """``flex_second_order = "local"`` compiles a factor pack that is
        mirrored on the device for the whole solve, so it has to lower the
        admissible batch (and, at a tight enough cap, refuse the run)."""
        from hwave.solver.flex_bond import estimate_bond_memory
        cap = 1.1e7
        avail = int(cap / 0.9)
        bare = estimate_bond_memory(device_available=avail, **self._KW)
        fb = 2_000_000
        packed = estimate_bond_memory(device_available=avail, factor_bytes=fb, **self._KW)
        self.assertEqual(packed["device_rows"]["second_order_factors"], fb)
        self.assertLess(packed["device_nb"], bare["device_nb"])
        self.assertLessEqual(self._need(packed["device_nb"], fb), cap)
        # and a pack large enough to eat the remaining headroom refuses at nb = 1
        tight = int(self._need(1) * 1.01 / 0.9)
        estimate_bond_memory(device_available=tight, **self._KW)      # fits without the mirror
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=tight,
                                 factor_bytes=int(self._need(1)), **self._KW)
        self.assertIn("second_order_factors", str(cm.exception))

    def test_refusal_names_the_phase(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=1000, **self._KW)
        self.assertIn("device", str(cm.exception).lower())
        self.assertTrue("dressing" in str(cm.exception) or "transport" in str(cm.exception))

    def test_the_tail_row_can_decide_admission(self):
        """``BondDeviceContext`` allocates the tail AFTER the reading the
        device table is written against, so a run that fits without it and
        not with it has to be refused rather than fail mid-solve."""
        from hwave.solver.flex_bond import estimate_bond_memory
        S, G, C, per = self._sym()
        phase = max(per * 1 + self._guard(1), 6 * C, self._bubble())
        without = 1.25 * (5 * S + 5 * G + phase)
        with_tail = 1.25 * (5 * S + 5 * G + G + phase)
        self.assertGreater(with_tail, without)
        cap = 0.5 * (without + with_tail)
        avail = int(cap / 0.9)
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=avail, **self._KW)
        self.assertIn("green0_tail", str(cm.exception))

    def test_explicit_batch_validated_against_both_tables(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        # Pick device_available so the nb=64 device need clears the cap (but
        # nb=1 still fits, so this exercises the explicit-batch check rather
        # than the generic phase-at-nb=1 refusal).
        it = 16; nvol = 16; P = 4; ND = 3 * P; nmat = 64
        S = nvol * ND * ND * it; C = nmat * nvol * P * P * it
        dev_need_64 = 1.25 * (5 * S + max(7 * 64 * nvol * ND * ND * it, 6 * C))
        device_available = int(dev_need_64 / 0.9 * 0.8)
        kw = dict(self._KW, freq_batch=64)
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=device_available, **kw)
        msg = str(cm.exception)
        self.assertIn("longitudinal_bond_freq_batch", msg)
        self.assertIn("device", msg.lower())


class TestHostBubbleRowOnTheDevicePath(unittest.TestCase):
    """With the bubble on the device the HOST no longer allocates its
    temporaries: it holds one finished channel-pair block at a time on its
    way into the store, i.e. ``C`` bytes (issue #196). The host-only
    preflight call, which runs before the backend is known, keeps the CPU
    expression."""

    _KW = dict(nmat=64, nvol=16, norb=2, B=3, depth=4, output_full=False, split_seed=False,
               n_types=1, freq_batch=None, cap_gb=200.0, mixing="anderson")

    @classmethod
    def _sym(cls):
        it, nvol, P, ND, nmat = 16, 16, 4, 3 * 4, 64
        prep = it * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        return max(prep, pair), nmat * nvol * P * P * it

    def test_host_only_call_keeps_the_cpu_expression(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        cpu_row, _C = self._sym()
        est = estimate_bond_memory(**self._KW)
        self.assertEqual(est["phase_rows"]["bubble"], cpu_row)

    def test_with_a_device_reading_the_host_row_is_one_block(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        cpu_row, C = self._sym()
        self.assertNotEqual(cpu_row, C)
        est = estimate_bond_memory(device_available=2 ** 62, **self._KW)
        self.assertEqual(est["phase_rows"]["bubble"], C)


class TestDevicePhasePrecedence(unittest.TestCase):
    """Which phase row a refusal names when two of them are exactly equal:
    bubble, then dressing, then transport."""

    def test_helper_breaks_ties_in_order(self):
        from hwave.solver.flex_bond import _largest_dev_phase
        self.assertEqual(_largest_dev_phase({"bubble": 5, "dressing": 5, "transport": 5}),
                         "bubble")
        self.assertEqual(_largest_dev_phase({"bubble": 4, "dressing": 5, "transport": 5}),
                         "dressing")
        self.assertEqual(_largest_dev_phase({"bubble": 4, "dressing": 4, "transport": 5}),
                         "transport")
        self.assertEqual(_largest_dev_phase({"bubble": 6, "dressing": 4, "transport": 5}),
                         "bubble")

    def test_a_real_bubble_transport_tie_names_the_bubble(self):
        """norb = 1, B = 4, Nmat = 40: ``pair`` and ``6 C`` coincide."""
        from hwave.solver.flex_bond import estimate_bond_memory
        kw = dict(nmat=40, nvol=8, norb=1, B=4, depth=0, output_full=False, split_seed=False,
                  n_types=1, freq_batch=None, cap_gb=200.0, mixing="linear")
        it, nvol, P, ND, nmat = 16, 8, 1, 4, 40
        pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        self.assertEqual(pair, 6 * nmat * nvol * P * P * it)
        self.assertGreater(pair, 7 * nvol * ND * ND * it)
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=1000, **kw)
        self.assertIn("'bubble'", str(cm.exception))

    def test_a_real_dressing_transport_tie_names_the_dressing(self):
        """norb = 2, B = 6, Nmat = 42: the dressing batch and ``6 C``
        coincide, and both are above the bubble."""
        from hwave.solver.flex_bond import estimate_bond_memory
        kw = dict(nmat=42, nvol=4, norb=2, B=6, depth=0, output_full=False, split_seed=False,
                  n_types=1, freq_batch=None, cap_gb=200.0, mixing="linear")
        it, nvol, P, ND, nmat = 16, 4, 4, 24, 42
        pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        dressing = 7 * nvol * ND * ND * it
        self.assertEqual(dressing, 6 * nmat * nvol * P * P * it)
        self.assertGreater(dressing, pair)
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=1000, **kw)
        self.assertIn("'dressing'", str(cm.exception))


class TestDeviceBubbleRow(unittest.TestCase):
    """A shape whose BUBBLE is the largest of the three device phase rows
    (issue #196). ``B = 5``, ``norb = 1``, ``Nmat = 44``, ``nvol = 8``:
    the bubble's per-pair buffers grow with ``Nmat`` while the dressing
    batch (plus its incremental guard row, issue #197) does not, and the
    transport's ``6 C`` stays below them at ``norb = 1``. ``Nmat`` is 44,
    not 42: the guard row's factor was raised from 3 to 4 (issue #197
    review), which would otherwise put the guarded dressing row above the
    bubble at ``Nmat = 42``; 44 is the minimal integer that keeps the
    bubble dominant."""

    _KW = dict(nmat=44, nvol=8, norb=1, B=5, depth=0, output_full=False, split_seed=False,
               n_types=1, freq_batch=None, cap_gb=200.0, mixing="linear")

    _IT, _NVOL, _P, _ND, _NMAT, _B = 16, 8, 1, 5, 44, 5

    @classmethod
    def _rows(cls):
        it, nvol, P, ND, nmat = cls._IT, cls._NVOL, cls._P, cls._ND, cls._NMAT
        prep = it * nvol * (ND ** 2 + 4 * P * (nmat + 2))
        pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
        return dict(bubble=max(prep, pair),
                    transport=6 * nmat * nvol * P * P * it,
                    dressing1=7 * 1 * nvol * ND * ND * it,
                    guard1=4 * 1 * nvol * ND * ND * it,
                    # 5 S vertices, 5 G the SCF arrays, one more G the tail
                    persistent=5 * nvol * ND * ND * it + 6 * nmat * nvol * P * it)

    def test_the_fixture_really_is_bubble_dominated(self):
        r = self._rows()
        self.assertGreater(r["bubble"], r["transport"])
        self.assertGreater(r["bubble"], r["dressing1"] + r["guard1"])

    def test_device_need_takes_the_bubble_row(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        r = self._rows()
        need1 = 1.25 * (r["persistent"] + r["bubble"])
        need2 = 1.25 * (r["persistent"] + 2 * (r["dressing1"] + r["guard1"]))
        self.assertGreater(need2, need1)          # nb = 2 is dressing-dominated
        # admit nb = 1 and refuse nb = 2, so the selected batch is 1 and the
        # bubble is the phase row the need is built on
        avail = int((need1 + need2) / 2 / 0.9)
        est = estimate_bond_memory(device_available=avail, **self._KW)
        self.assertEqual(est["device_nb"], 1)
        self.assertEqual(est["device_rows"]["bubble"], r["bubble"])
        self.assertAlmostEqual(est["device_need"], need1, delta=1e-6 * need1)

    def test_refusal_names_the_bubble_when_it_is_the_largest_row(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        with self.assertRaises(ValueError) as cm:
            estimate_bond_memory(device_available=1000, **self._KW)
        msg = str(cm.exception)
        self.assertIn("device", msg.lower())
        self.assertIn("'bubble'", msg)

    def test_the_table_print_carries_the_row(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        r = self._rows()
        need1 = 1.25 * (r["persistent"] + r["bubble"])
        est = estimate_bond_memory(device_available=int(need1 * 1.5 / 0.9), **self._KW)
        self.assertIn("bubble", est["device_table"])


if __name__ == "__main__":
    unittest.main()
