"""The static RPA pole of the bond paths (issue #345, item 1).

The conditioning guard of the bond dressing refuses a denominator
``1 -/+ chibar W`` only within ``cond_tol`` of its pole. Past the pole (the
largest real eigenvalue of ``-sign chibar(q, 0) W`` above 1) the denominator
is well conditioned again, the dressed susceptibility is no longer positive
and the pairing eigenvalue describes no physical state, yet the run used to
finish without a word. Every bond path now computes that indicator, logs it
at INFO and, past the pole, logs ONE WARNING and writes the provenance keys
``rpa_past_pole``, ``rpa_pole_indicator_<channel>`` and ``rpa_pole_q_<channel>``.
The pairing eigenvalue and ``eigen_converged`` are unchanged, and a run below
the pole writes nothing new.

Each path has a pair of cheap fixtures straddling the pole. The checks:

* below: no new WARNING and outputs identical to a run with the check
  stubbed out (no new keys);
* above: exactly one new WARNING naming the channel, the indicator and q,
  the keys present, the pairing eigenvalue identical to the stubbed run;
* an independent oracle: past the pole exactly when the Hermitian part of
  the dressed static susceptibility has a negative eigenvalue;
* the conditioning guard score of the above fixture is above ``cond_tol``,
  so it is exactly the case the guard does not see.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

COND_TOL = 1.0e-3
POLE_TEXT = "past its pole"


class _Records(logging.Handler):
    def __init__(self):
        super().__init__(level=logging.INFO)
        self.records = []

    def emit(self, record):
        self.records.append(record)


@contextlib.contextmanager
def _captured():
    root = logging.getLogger()
    h = _Records()
    level = root.level
    root.addHandler(h)
    root.setLevel(logging.INFO)
    try:
        yield h.records
    finally:
        root.removeHandler(h)
        root.setLevel(level)


def _warnings(records):
    return [r.getMessage() for r in records if r.levelno >= logging.WARNING]


def _pole_warnings(records):
    return [m for m in _warnings(records) if POLE_TEXT in m]


@contextlib.contextmanager
def _stubbed():
    """The pole check switched off: the outputs of the code before it."""
    from hwave.solver import bond_channels as bc

    def stub(indicators, **kw):
        return bc.RpaPoleReport({})
    with mock.patch.object(bc, "check_rpa_pole", stub):
        yield


@contextlib.contextmanager
def _failing_diagnostic():
    """The pole diagnostic's eigenvalue solve raising, as LAPACK can."""
    from hwave.solver import bond_channels as bc

    def boom(m):
        raise np.linalg.LinAlgError("Eigenvalues did not converge")
    with mock.patch.object(bc, "_pole_eigvals", boom):
        yield


def _modes(stub, failing):
    st = contextlib.ExitStack()
    if stub:
        st.enter_context(_stubbed())
    if failing:
        st.enter_context(_failing_diagnostic())
    return st


def _min_herm_eig(chi):
    chi = np.asarray(chi)
    h = 0.5 * (chi + np.swapaxes(chi.conj(), -1, -2))
    return float(np.linalg.eigvalsh(h).min())


def _cond_score(chibar, W, sign):
    from hwave.solver import bond_channels as bc
    mat = np.eye(W.shape[-1]) + sign * (np.asarray(chibar) @ np.asarray(W))
    return bc.bond_conditioning_score(mat.reshape(-1, 1, 1, *mat.shape[-2:]))[0]


# ---------------------------------------------------------------------------
# the helper
# ---------------------------------------------------------------------------

class TestIndicatorHelper(unittest.TestCase):

    def test_one_by_one_block_is_u_chi0(self):
        """``ND = 1``: the spin indicator is ``U max_q chi0(q)``, the charge
        one ``-U min_q chi0(q)``, at the q of that extremum."""
        from hwave.solver import bond_channels as bc
        shape = (2, 2, 1)
        chi0 = np.array([0.1, 0.4, 0.25, 0.2]).reshape(4, 1, 1).astype(complex)
        for U in (0.5, 2.5, 4.0):
            W = np.full((4, 1, 1), U, dtype=complex)
            lam, q, _ = bc.static_pole_indicator(chi0, W, "spin", shape)
            self.assertAlmostEqual(lam, U * 0.4, delta=1e-14)
            self.assertEqual(q, (0.0, 0.5, 0.0))
            lam, q, _ = bc.static_pole_indicator(chi0, W, "charge", shape)
            self.assertAlmostEqual(lam, -U * 0.1, delta=1e-14)
            self.assertEqual(q, (0.0, 0.0, 0.0))

    def test_symmetry_equivalent_maxima_give_the_lowest_q(self):
        """Degenerate maxima (q and its images) differing by round-off give
        the lowest flat q index, so two equivalent runs report one q."""
        from hwave.solver import bond_channels as bc
        for eps in (0.0, 1e-15, -1e-15):
            chi0 = np.array([0.1, 0.4, 0.4 + eps, 0.2]).reshape(4, 1, 1).astype(complex)
            W = np.full((4, 1, 1), 2.0, dtype=complex)
            lam, q, _ = bc.static_pole_indicator(chi0, W, "spin", (2, 2, 1))
            self.assertEqual(q, (0.0, 0.5, 0.0))

    def _problem(self, seed, scale):
        rng = np.random.default_rng(seed)
        nvol, nd = 4, 3
        a = rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd))
        chibar = a @ np.swapaxes(a.conj(), -1, -2) / nd + 0.1 * np.eye(nd)
        b = rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd))
        W = 0.5 * (b + np.swapaxes(b.conj(), -1, -2))
        return chibar, W * scale

    def test_inversion_of_the_dressed_susceptibility(self):
        """``chi = [1 + sign chibar W]^-1 chibar``: the eigenvalues mu of
        ``-sign chi W`` are ``lambda / (1 - lambda)``; inverting them gives
        the indicator of chibar back, on both sides of the pole."""
        from hwave.solver import bond_channels as bc
        shape = (4, 1, 1)
        seen = set()
        for channel in ("spin", "charge", "spinful"):
            chibar, W0 = self._problem(3, 1.0)
            lam0, _, _ = bc.static_pole_indicator(chibar, W0, channel, shape)
            for target in (0.6, 1.4):
                W = W0 * (target / lam0)
                lam, q, _ = bc.static_pole_indicator(chibar, W, channel, shape)
                chi, _ = bc.dress_batch(chibar[None], W, channel, l0=0, nmat=0,
                                        spatial_shape=shape, cond_tol=None)
                lam_d, q_d, _ = bc.pole_indicator_from_dressed(chi[0], W, channel, shape)
                self.assertAlmostEqual(lam, target, delta=1e-10)
                self.assertAlmostEqual(lam_d, lam, delta=1e-10)
                self.assertEqual(q_d, q)
                # the independent oracle: the dressed chi turns indefinite
                self.assertEqual(_min_herm_eig(chi[0]) < 0, lam >= 1.0)
                seen.add((channel, lam >= 1.0))
        self.assertEqual(len(seen), 6)

    def test_a_complex_pair_does_not_count(self):
        """Only a REAL eigenvalue reaching 1 makes the denominator singular."""
        from hwave.solver import bond_channels as bc
        rot = np.array([[1.5, -1.0], [1.0, 1.5]], dtype=complex)    # 1.5 +- 1i
        chibar = np.eye(2, dtype=complex)[None]
        W = -rot[None]                                                 # charge: eig(-chibar W)
        lam, q, max_re = bc.static_pole_indicator(chibar, W, "charge", (1, 1, 1))
        self.assertTrue(np.isnan(lam))
        self.assertAlmostEqual(max_re, 1.5, delta=1e-12)
        rep = bc.check_rpa_pole({"charge": (lam, q, max_re)}, where="test",
                                consequence="none.", log=logging.getLogger("t"))
        self.assertFalse(rep.past)
        self.assertEqual(rep.npz_members(), {})

    def test_report_members_only_past_the_pole(self):
        from hwave.solver import bond_channels as bc
        log = logging.getLogger("hwave.test.pole")
        below = {"spin": (0.99, (0.5, 0.5, 0.0), 0.99), "charge": (0.2, (0.0, 0.0, 0.0), 0.2)}
        with self.assertLogs(log, level="INFO") as cm:
            rep = bc.check_rpa_pole(below, where="test", consequence="none.", log=log)
        self.assertFalse(rep.past)
        self.assertEqual(rep.npz_members(), {})
        self.assertEqual(rep.text_members(), {})
        self.assertEqual([r.levelname for r in cm.records], ["INFO"])
        self.assertIn("spin 0.990000 at q = (0.5, 0.5, 0)", cm.output[0])
        above = dict(below, charge=(1.25, (0.5, 0.5, 0.0), 1.25))
        with self.assertLogs(log, level="INFO") as cm:
            rep = bc.check_rpa_pole(above, where="test", consequence="Some consequence.",
                                    log=log)
        self.assertTrue(rep.past)
        self.assertEqual([r.levelname for r in cm.records], ["INFO", "WARNING"])
        msg = cm.records[1].getMessage()
        for part in ("charge", "1.250000", "(0.5, 0.5, 0)", "Some consequence.", POLE_TEXT,
                     "rpa_past_pole = true"):
            self.assertIn(part, msg)
        self.assertNotIn("spin channel", msg)
        m = rep.npz_members()
        self.assertEqual(set(m), {"rpa_past_pole", "rpa_pole_indicator_spin",
                                  "rpa_pole_q_spin", "rpa_pole_indicator_charge",
                                  "rpa_pole_q_charge"})
        self.assertIs(bool(m["rpa_past_pole"]), True)
        self.assertEqual(float(m["rpa_pole_indicator_charge"]), 1.25)
        np.testing.assert_array_equal(m["rpa_pole_q_charge"], [0.5, 0.5, 0.0])
        self.assertEqual(set(rep.npz_members("longitudinal_bond_")),
                         {"longitudinal_bond_" + k for k in m})
        t = rep.text_members()
        self.assertEqual(t["rpa_past_pole"], "true")
        self.assertEqual(t["rpa_pole_q_charge"], "(0.5, 0.5, 0)")
        self.assertEqual(float(t["rpa_pole_indicator_charge"]), 1.25)
        # the threshold is inclusive
        rep = bc.check_rpa_pole({"spin": (1.0, (0.0, 0.0, 0.0), 1.0)}, where="t",
                                consequence=".", log=log)
        self.assertTrue(rep.past)


class TestFlexMapIndicatorOnTheDevice(unittest.TestCase):
    """The FLEX map computes its indicators from the bare static chibar
    before the W batch overwrites it, on a host or a (fake) device store
    alike, with one static slice copied to the host on the device path."""

    def test_host_and_device_store_agree(self):
        from hwave.solver import bond_channels as bc, flex_bond as fb
        from tests.fake_device import fake_cupy
        from tests.test_bond_half_grid import _hermitian_in_frequency
        from tests.test_flex_bond_dressing import _problem
        from tests.test_flex_bond_gpu import _onsite_factors_on
        import types
        nmat, nvol, nd, B = 8, 4, 4, 3
        ND = B * nd
        _x, S, C = _problem(nmat=nmat, nvol=nvol, nd=nd, B=B, seed=5)
        x = 0.05 * _hermitian_in_frequency(nmat, nvol, ND, seed=5)
        want = {ch: bc.static_pole_indicator(x[nmat // 2], W, ch, (nvol, 1, 1))
                for ch, W in (("spin", S), ("charge", C))}
        view = types.SimpleNamespace(n_channels=B)
        got = {}
        with fake_cupy() as cp:
            for name, sx in (("host", np), ("device", cp)):
                with fb.BondBlockStore(nmat, nvol, ND, nd, ("chibar",), {"W": "chibar"},
                                       xp=sx) as st, \
                        fb.BondDeviceContext.for_view(cp, S, C, view, 2) as dev:
                    st.put_freq_batch("chibar", 0, nmat, x)
                    res = fb.dress_and_build_w(st, dev, nb=3, output_full=False, nmat=nmat,
                                               nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                               factors=_onsite_factors_on(cp),
                                               pole_vertices=(S, C))
                    got[name] = res.pole
                    # without the vertices nothing is computed
                    st.put_freq_batch("chibar", 0, nmat, x)
                    res0 = fb.dress_and_build_w(st, dev, nb=3, output_full=False, nmat=nmat,
                                                nvol=nvol, nd=nd, spatial_shape=(nvol, 1, 1),
                                                factors=_onsite_factors_on(cp))
                    self.assertIsNone(res0.pole)
        for name in ("host", "device"):
            self.assertEqual(set(got[name]), {"spin", "charge"})
            for ch in want:
                self.assertAlmostEqual(got[name][ch][0], want[ch][0], delta=1e-12)
                self.assertEqual(got[name][ch][1], want[ch][1])


# ---------------------------------------------------------------------------
# non-SO static (hwave_sc, bond_channels = true)
# ---------------------------------------------------------------------------

class _DressSpy:
    """Records the static chibar, W and dressed chi of every bond dressing."""

    def __init__(self):
        from hwave.solver import bond_channels as bc
        self.calls = []
        self._real = bc.dress_batch

    def __call__(self, cb, W, channel, **kw):
        chi, cond = self._real(cb, W, channel, **kw)
        self.calls.append((channel, np.array(cb[0]), np.array(W),
                           None if chi is None else np.array(chi[0]), cond))
        return chi, cond


def _run_nonso_static(U, V, stub=False, failing=False):
    """The 4x4 one-band hwave_sc bond run (tests/rpa/input band, T = 1,
    Nmat = 128, on-site U and a nearest-neighbour V); returns (eigenvalue.dat
    text, records, spy calls)."""
    import hwave.sc as sc
    from hwave.solver import bond_channels as bc
    from tests.test_sc_bond import _write_w90, _nn_entries, _base_input
    d = tempfile.mkdtemp(prefix="hwave_pole_")
    try:
        ci = os.path.join(d, "ci.dat")
        _write_w90(ci, 1, _nn_entries(V))
        cu = os.path.join(d, "cu.dat")
        _write_w90(cu, 1, [((0, 0, 0), (0, 0), U)])
        out = os.path.join(d, "out")
        inp = _base_input(out, coulomb_inter=ci, bond_channels=True, solver_mode="eigenvalue")
        inp["file"]["input"]["interaction"]["CoulombIntra"] = cu
        spy = _DressSpy()
        with contextlib.ExitStack() as st:
            st.enter_context(_modes(stub, failing))
            st.enter_context(mock.patch.object(bc, "dress_batch", spy))
            records = st.enter_context(_captured())
            sc.calc_eliashberg(inp)
        with open(os.path.join(out, "eigenvalue.dat")) as f:
            text = f.read()
        return text, records, spy.calls
    finally:
        shutil.rmtree(d, ignore_errors=True)


def _comment_keys(text, sep=" = "):
    keys = {}
    for line in text.splitlines():
        if line.startswith("# ") and sep in line:
            k, v = line[2:].split(sep, 1)
            keys[k.strip()] = v.strip()
    return keys


def _drop_pole_lines(text):
    return "".join(l for l in text.splitlines(True) if "rpa_p" not in l)


def _stable_text(text):
    """``text`` without the lines that record the memory available at run
    time (the spin-orbital memory admission), which differ run to run."""
    return "".join(l for l in text.splitlines(True) if "available" not in l)


class TestNonSOStatic(unittest.TestCase):
    """spin: U = 5.3464 (indicator 0.990) / 5.5085 (1.020) at V = 0.1;
    charge: V = 0.78 (0.970) / 0.82 (1.030) at U = 1; both at q = (pi, pi)."""

    def _pair(self, U, V, channel, below_value, above):
        text, rec, calls = _run_nonso_static(U, V)
        text0, rec0, _ = _run_nonso_static(U, V, stub=True)
        from hwave.solver import bond_channels as bc
        static = {c[0]: c for c in calls}
        ch, cb, W, chi, cond = static[channel]
        sign = bc._DRESS_CHANNELS[channel]
        lam, q, _ = bc.static_pole_indicator(cb, W, channel, (4, 4, 1))
        self.assertAlmostEqual(lam, below_value, delta=5e-4)
        self.assertEqual(q, (0.5, 0.5, 0.0))
        # chibar(q, 0) is Hermitian positive definite, so the spectrum is real
        self.assertLess(np.abs(cb - np.swapaxes(cb.conj(), -1, -2)).max(), 1e-12)
        self.assertGreater(_min_herm_eig(cb), 0.0)
        # the oracle and the guard gap
        self.assertEqual(_min_herm_eig(chi) < 0, above)
        self.assertGreater(_cond_score(cb, W, sign), COND_TOL)
        self.assertGreater(cond, COND_TOL)
        # INFO always
        info = [r.getMessage() for r in rec if "RPA pole indicator" in r.getMessage()]
        self.assertEqual(len(info), 1)
        self.assertIn("{} {:.6f}".format(channel, lam), info[0])
        keys = _comment_keys(text)
        if not above:
            self.assertEqual(text, text0)
            self.assertEqual(_warnings(rec), _warnings(rec0))
            self.assertFalse(any(k.startswith("rpa_p") for k in keys))
            return
        self.assertEqual(_pole_warnings(rec0), [])
        new = [m for m in _warnings(rec) if m not in _warnings(rec0)]
        self.assertEqual(len(new), 1)
        for part in (channel, "{:.6f}".format(lam), "(0.5, 0.5, 0)", POLE_TEXT,
                     "pairing eigenvalue"):
            self.assertIn(part, new[0])
        self.assertEqual(keys["rpa_past_pole"], "true")
        self.assertAlmostEqual(float(keys["rpa_pole_indicator_" + channel]), lam, delta=1e-7)
        self.assertEqual(keys["rpa_pole_q_" + channel], "(0.5, 0.5, 0)")
        # everything else, the eigenvalues and eigen_converged included, is the
        # stubbed run's byte for byte
        self.assertEqual(_drop_pole_lines(text), text0)
        self.assertIn("# eigen_converged = true", text)

    def test_spin_below(self):
        self._pair(5.3464, 0.1, "spin", 0.990, False)

    def test_spin_above(self):
        self._pair(5.5085, 0.1, "spin", 1.020, True)

    def test_charge_below(self):
        self._pair(1.0, 0.78, "charge", 0.970, False)

    def test_charge_above(self):
        self._pair(1.0, 0.82, "charge", 1.030, True)


# ---------------------------------------------------------------------------
# SO static (dense and S_z-packed)
# ---------------------------------------------------------------------------

def _so_fixture(soc, V, N):
    from tests import so_bond_entry_fixtures as F
    return F.SOC_ALL_TYPES._replace(
        name="POLE_{}_{}_{}".format(soc, V, N), shape=(N, N, 1),
        transfer=(F._SOC if soc else ()),
        interactions=(("CoulombIntra", (((0, 0, 0), (0, 0), 0.5),)),
                      ("CoulombInter", F._nn(V))))


def _run_so(soc, V, N, stub=False, failing=False):
    """(RunResult, [(chibar, W, chi)] per dressed block)."""
    from tests import so_bond_entry_fixtures as F
    from hwave.solver import bond_channels as bc, so_bond as sb
    blocks = []
    real_dense, real_packed = bc.dress_batch, sb._dress_packed

    def dense(cb, W, channel, **kw):
        chi, cond = real_dense(cb, W, channel, **kw)
        blocks.append((np.array(cb[0]), np.array(W), np.array(chi[0])))
        return chi, cond

    def packed(X, g, **kw):
        out = real_packed(X, g, **kw)
        for s in range(3):
            if X.blocks[s].shape[-1]:
                blocks.append((np.array(X.blocks[s][0]), np.array(g.blocks[s]),
                               np.array(out.blocks[s][0])))
        return out
    with contextlib.ExitStack() as st:
        st.enter_context(_modes(stub, failing))
        st.enter_context(mock.patch.object(bc, "dress_batch", dense))
        st.enter_context(mock.patch.object(sb, "_dress_packed", packed))
        r = F.run_case(_so_fixture(soc, V, N), fresh=True)
    return r, blocks


class TestSOStatic(unittest.TestCase):
    """dense: the 6x6 Rashba band, U = 0.5, V = 0.52 (0.962) / 0.60 (1.130);
    S_z-packed: the 4x4 band without SOC, V = 0.50 (0.894) / 0.60 (1.0985)."""

    def _pair(self, soc, V, N, layout, value, above):
        from hwave.solver import bond_channels as bc
        r, blocks = _run_so(soc, V, N)
        r0, _ = _run_so(soc, V, N, stub=True)
        self.assertEqual(str(r.gap_npz["block_layout"]), layout)
        lams = [bc.static_pole_indicator(cb, W, "spinful", (N, N, 1)) for cb, W, _ in blocks]
        lam = max(x[0] for x in lams)
        self.assertAlmostEqual(lam, value, delta=5e-4)
        # Hermitian positive definite static bubble: eig(-chibar Gamma) has the
        # spectrum of the Hermitian -L^H Gamma L, chibar = L L^H
        for cb, W, chi in blocks:
            self.assertLess(np.abs(cb - np.swapaxes(cb.conj(), -1, -2)).max(), 1e-12)
            L = np.linalg.cholesky(cb)
            herm = -np.swapaxes(L.conj(), -1, -2) @ W @ L
            ev_h = np.sort(np.linalg.eigvalsh(0.5 * (herm + np.swapaxes(herm.conj(), -1, -2))))
            ev = np.linalg.eigvals(-(cb @ W))
            np.testing.assert_allclose(np.sort(ev.real), ev_h, atol=1e-10)
            self.assertGreater(_cond_score(cb, W, +1.0), COND_TOL)
        self.assertEqual(min(_min_herm_eig(chi) for _, _, chi in blocks) < 0, above)
        info = [x.getMessage() for x in r.records if "RPA pole indicator" in x.getMessage()]
        self.assertEqual(len(info), 1)
        self.assertEqual(r.value, r0.value)
        new = [m for m in _warnings(r.records) if m not in _warnings(r0.records)]
        pole_keys = {k for k in r.gap_npz if k.startswith("rpa_p")}
        if not above:
            self.assertEqual(new, [])
            self.assertEqual(pole_keys, set())
            self.assertEqual(_stable_text(r.eigenvalue_dat), _stable_text(r0.eigenvalue_dat))
            self.assertEqual(set(r.gap_npz), set(r0.gap_npz))
            for k in r.gap_npz:
                np.testing.assert_array_equal(r.gap_npz[k], r0.gap_npz[k], err_msg=k)
            return
        self.assertEqual(len(new), 1)
        for part in ("spinful", "{:.6f}".format(lam), "(0.5, 0.5, 0)", POLE_TEXT):
            self.assertIn(part, new[0])
        self.assertEqual(pole_keys, {"rpa_past_pole", "rpa_pole_indicator_spinful",
                                     "rpa_pole_q_spinful"})
        self.assertIs(bool(r.gap_npz["rpa_past_pole"]), True)
        self.assertAlmostEqual(float(r.gap_npz["rpa_pole_indicator_spinful"]), lam, delta=1e-10)
        np.testing.assert_array_equal(r.gap_npz["rpa_pole_q_spinful"], [0.5, 0.5, 0.0])
        keys = _comment_keys(r.eigenvalue_dat)
        self.assertEqual(keys["rpa_past_pole"], "true")
        self.assertEqual(_stable_text(_drop_pole_lines(r.eigenvalue_dat)),
                         _stable_text(r0.eigenvalue_dat))
        for k in r0.gap_npz:
            np.testing.assert_array_equal(r.gap_npz[k], r0.gap_npz[k], err_msg=k)
        self.assertIs(bool(r.gap_npz["eigen_converged"]), bool(r0.gap_npz["eigen_converged"]))

    def test_dense_below(self):
        self._pair(True, 0.52, 6, "dense", 0.962, False)

    def test_dense_above(self):
        self._pair(True, 0.60, 6, "dense", 1.130, True)

    def test_packed_below(self):
        self._pair(False, 0.50, 4, "sz_packed", 0.894, False)

    def test_packed_above(self):
        self._pair(False, 0.60, 4, "sz_packed", 1.0985, True)


# ---------------------------------------------------------------------------
# FLEX bond final state, its in-process pairing and the archive reader
# ---------------------------------------------------------------------------

_FLEX_BAND = "tests/rpa/input"


def _write_uv(d, U, V):
    from tests.test_sc_bond import _write_w90, _nn_entries
    os.makedirs(d, exist_ok=True)
    for n in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(_FLEX_BAND, n), d)
    _write_w90(os.path.join(d, "coulombintra.dat"), 1, [((0, 0, 0), (0, 0), U)])
    _write_w90(os.path.join(d, "coulombinter.dat"), 1, _nn_entries(V))


def _interaction(d):
    return {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat", "CoulombInter": "coulombinter.dat"}


_NMAT = 16
_T = 0.5


def _flex_input(d, out, itmax, bond=True, extra=None):
    par = {"T": _T, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": _NMAT, "IterationMax": itmax, "Mix": 0.3, "EPS": 8,
           "flex_hartree_fock": True, "mixing_scheme": "linear",
           "longitudinal_bond_channels": bond}
    if bond:
        par.update(longitudinal_bond_output_full=True, longitudinal_bond_pairing="singlet")
    par.update(extra or {})
    return {"log": {"print_level": 0},
            "mode": {"mode": "FLEX", "calc_scheme": "general", "param": par},
            "file": {"input": {"path_to_input": "", "interaction": _interaction(d)},
                     "output": {"path_to_output": out, "chiq_s": "chiq_s", "chiq_c": "chiq_c",
                                "green": "green", "longitudinal_bond": "longitudinal_bond"}}}


def _sc_dynamic_input(d, flex_out, out):
    return {"mode": {"param": {"T": _T, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
                               "Nmat": _NMAT, "filling": 0.5}},
            "file": {"input": {"path_to_flex_output": flex_out, "interaction": _interaction(d)},
                     "output": {"path_to_output": out}},
            "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic", "bond_channels": True,
                           "pairing_type": "singlet", "solver_mode": "eigenvalue",
                           "num_eigenvalues": 3}}


def _npz(path):
    with np.load(path, allow_pickle=True) as z:
        return {k: z[k] for k in z.files}


def _read(path):
    with open(path) as f:
        return f.read()


class _FlexCase:
    """One FLEX bond run (with the in-process singlet pairing and the full
    archive) and one hwave_sc dynamic run on its archive, with and without
    the pole check."""

    def __init__(self, U, V, itmax, stub=False, failing=False, extra=None):
        import hwave.qlms as qlms
        import hwave.sc as sc
        self.tmp = tempfile.mkdtemp(prefix="hwave_pole_flex_")
        d = os.path.join(self.tmp, "in")
        _write_uv(d, U, V)
        self.flex = os.path.join(self.tmp, "flex")
        self.sc = os.path.join(self.tmp, "sc")
        with contextlib.ExitStack() as st:
            st.enter_context(_modes(stub, failing))
            self.flex_records = st.enter_context(_captured())
            qlms.run(input_dict=_flex_input(d, self.flex, itmax, extra=extra))
        with contextlib.ExitStack() as st:
            st.enter_context(_modes(stub, failing))
            self.sc_records = st.enter_context(_captured())
            self.sc_value = sc.calc_eliashberg(_sc_dynamic_input(d, self.flex, self.sc))
        self.arch = _npz(os.path.join(self.flex, "longitudinal_bond.npz"))
        self.chiq_s = _npz(os.path.join(self.flex, "chiq_s.npz"))
        self.pair = _npz(os.path.join(self.flex, "eliashberg_bond_singlet.npz"))
        self.pair_dat = _read(os.path.join(self.flex, "eigenvalue_bond_singlet.dat"))
        self.dyn = _npz(os.path.join(self.sc, "gap_dynamic.npz"))
        self.dyn_dat = _read(os.path.join(self.sc, "eigenvalue.dat"))

    def close(self):
        shutil.rmtree(self.tmp, ignore_errors=True)


_FLEX_CASES = {}


def _flex_case(key, *args, stub=False, failing=False, extra=None):
    k = (key, stub, failing, tuple(sorted((extra or {}).items())))
    if k not in _FLEX_CASES:
        _FLEX_CASES[k] = _FlexCase(*args, stub=stub, failing=failing, extra=extra)
    return _FLEX_CASES[k]


def tearDownModule():
    for c in _FLEX_CASES.values():
        c.close()
    _FLEX_CASES.clear()


_FLEX_FIXTURES = {
    # converged FLEX states, charge channel at q = (pi, pi)
    "charge_below": ((4.0, 1.3, 200), "charge", 0.9202),
    "charge_above": ((4.0, 1.6, 200), "charge", 1.3999),
    # the first map at V = 0 (IterationMax = 1), spin channel
    "spin_above": ((5.0, 0.0, 1), "spin", 1.1316),
}

_PREFIX = "longitudinal_bond_"


class TestFlexBond(unittest.TestCase):

    def _cases(self, name):
        args, channel, value = _FLEX_FIXTURES[name]
        return (_flex_case(name, *args), _flex_case(name, *args, stub=True), channel, value)

    def _indicator(self, case, channel):
        """The producer indicator recomputed from the archived dressed static
        chi (the inversion) -- and the independent oracle on the same chi."""
        from hwave.solver import bond_channels as bc
        chi = case.arch["longitudinal_bond_chi_" + channel[0]]
        W = case.arch["S_bond" if channel == "spin" else "C_bond"]
        return bc.pole_indicator_from_dressed(chi, W, channel, (4, 4, 1)), _min_herm_eig(chi)

    def _check(self, name, above):
        case, ref, channel, value = self._cases(name)
        (lam, q, _), min_eig = self._indicator(case, channel)
        self.assertAlmostEqual(lam, value, delta=5e-4)
        self.assertEqual(min_eig < 0, above)
        self.assertGreater(float(case.arch["longitudinal_bond_cond_min_" + channel[0]]), COND_TOL)
        # the final-iteration Stoner line carries the bond indicator
        stoner = [r.getMessage() for r in case.flex_records
                  if "RPA pole indicator (final map)" in r.getMessage()]
        self.assertFalse(any("Stoner factor" in r.getMessage() for r in case.flex_records))
        self.assertTrue(stoner)
        self.assertIn("{} {:.6f}".format(channel, lam), stoner[-1])
        # the pairing eigenvalues are the stubbed run's
        self.assertEqual(float(case.pair["eigenvalue"]), float(ref.pair["eigenvalue"]))
        self.assertEqual(case.sc_value, ref.sc_value)
        flex_new = [m for m in _warnings(case.flex_records) if m not in _warnings(ref.flex_records)]
        sc_new = [m for m in _warnings(case.sc_records) if m not in _warnings(ref.sc_records)]
        pole_keys = {
            "arch": {k for k in case.arch if "rpa_p" in k},
            "chiq_s": {k for k in case.chiq_s if "rpa_p" in k},
            "pair": {k for k in case.pair if "rpa_p" in k},
            "dyn": {k for k in case.dyn if "rpa_p" in k}}
        if not above:
            self.assertEqual(flex_new, [])
            self.assertEqual(sc_new, [])
            self.assertEqual(pole_keys, dict.fromkeys(pole_keys, set()))
            for got, want in ((case.arch, ref.arch), (case.chiq_s, ref.chiq_s),
                              (case.pair, ref.pair), (case.dyn, ref.dyn)):
                self.assertEqual(set(got), set(want))
                for k in want:
                    if k != "bond_archive":        # the archive's (temporary) path
                        np.testing.assert_array_equal(got[k], want[k], err_msg=k)
            self.assertEqual(case.pair_dat, ref.pair_dat)
            self.assertEqual(case.dyn_dat, ref.dyn_dat)
            return
        self.assertEqual(len(flex_new), 1)
        self.assertEqual(len(sc_new), 1)
        for msg in (flex_new[0], sc_new[0]):
            for part in (channel, "{:.6f}".format(lam), "(0.5, 0.5, 0)", POLE_TEXT):
                self.assertIn(part, msg)
        bare = {"rpa_past_pole", "rpa_pole_indicator_spin", "rpa_pole_q_spin",
                "rpa_pole_indicator_charge", "rpa_pole_q_charge"}
        self.assertEqual(pole_keys["arch"], {_PREFIX + k for k in bare})
        self.assertEqual(pole_keys["chiq_s"], {_PREFIX + k for k in bare})
        self.assertEqual(pole_keys["pair"], bare)
        self.assertEqual(pole_keys["dyn"], bare)
        producer = float(case.arch[_PREFIX + "rpa_pole_indicator_" + channel])
        # the producer's value (from chibar) and the reader's (the inversion of
        # the archived dressed chi) agree
        self.assertAlmostEqual(float(case.dyn["rpa_pole_indicator_" + channel]), producer,
                               delta=1e-10)
        self.assertAlmostEqual(lam, producer, delta=1e-10)
        self.assertEqual(float(case.pair["rpa_pole_indicator_" + channel]), producer)
        np.testing.assert_array_equal(case.arch[_PREFIX + "rpa_pole_q_" + channel],
                                      [0.5, 0.5, 0.0])
        self.assertIn("rpa_past_pole=true", case.pair_dat)
        self.assertIn("rpa_past_pole=true", case.dyn_dat)
        self.assertEqual("".join(l for l in case.pair_dat.splitlines(True) if "rpa_p" not in l),
                         ref.pair_dat)
        self.assertEqual("".join(l for l in case.dyn_dat.splitlines(True) if "rpa_p" not in l),
                         ref.dyn_dat)
        self.assertIs(bool(case.pair["eigen_converged"]), bool(ref.pair["eigen_converged"]))
        # every other member is the stubbed run's
        for got, want in ((case.arch, ref.arch), (case.chiq_s, ref.chiq_s),
                          (case.pair, ref.pair), (case.dyn, ref.dyn)):
            self.assertEqual({k for k in got if "rpa_p" not in k}, set(want))
            for k in want:
                if k != "bond_archive":
                    np.testing.assert_array_equal(got[k], want[k], err_msg=k)

    def test_charge_below(self):
        self._check("charge_below", False)

    def test_charge_above(self):
        self._check("charge_above", True)

    def test_spin_above(self):
        self._check("spin_above", True)

    def test_spin_at_zero_v_is_the_non_bond_stoner_factor(self):
        """At V = 0 the bond spin indicator is the non-bond FLEX Stoner factor
        of the same first map (an independent code path)."""
        import hwave.qlms as qlms
        case, _, channel, _ = self._cases("spin_above")
        (lam, _, _), _ = self._indicator(case, "spin")
        d = os.path.join(case.tmp, "in")
        out = os.path.join(case.tmp, "nonbond")
        with _captured() as rec:
            qlms.run(input_dict=_flex_input(d, out, 1, bond=False))
        line = [r.getMessage() for r in rec if "Stoner factor" in r.getMessage()][-1]
        value = float(line.split("spin ")[1].split()[0])
        self.assertAlmostEqual(lam, value, delta=1e-6)




# ---------------------------------------------------------------------------
# review round 1: adversarial spectra, a failing diagnostic, the half grid
# ---------------------------------------------------------------------------

class TestRealnessAfterTheInversion(unittest.TestCase):
    """mu / (1 + mu) next to mu = -1: realness is judged on the inverted
    value, so a negligible imaginary part of mu that the inversion blows up
    is not a real indicator."""

    def _inv(self, mus):
        from hwave.solver import bond_channels as bc
        ev = np.asarray(mus, dtype=complex).reshape(1, -1)
        return bc._pole_reduce(ev, lambda m: m / (1.0 + m), (1, 1, 1))

    def test_real_mu_on_either_side_of_minus_one(self):
        eps = 1e-6
        lam, q, _ = self._inv([-1.0 - eps, 0.5])
        self.assertAlmostEqual(lam, (1.0 + eps) / eps, delta=1e-3)     # past
        self.assertEqual(q, (0.0, 0.0, 0.0))
        lam, _, _ = self._inv([-1.0 + eps, 0.5])
        self.assertAlmostEqual(lam, 1.0 / 3.0, delta=1e-12)             # below
        self.assertLess(self._inv([-1.0 + eps])[0], -1e5)

    def test_tiny_imaginary_part_next_to_minus_one_is_complex(self):
        from hwave.solver import bond_channels as bc
        lam, _, max_re = self._inv([-1.0 + 1e-12j, 0.5])
        self.assertAlmostEqual(lam, 1.0 / 3.0, delta=1e-12)
        self.assertAlmostEqual(max_re, 1.0, delta=1e-6)
        lam, q, _ = self._inv([-1.0 + 1e-12j])
        self.assertTrue(np.isnan(lam))
        self.assertIsNone(q)
        rep = bc.check_rpa_pole({"charge": self._inv([-1.0 + 1e-12j, 0.5])}, where="t",
                                consequence=".", log=logging.getLogger("t"))
        self.assertFalse(rep.past)
        self.assertEqual(rep.npz_members(), {})

    def test_mu_exactly_minus_one_is_undetermined(self):
        """mu = -1 is lambda = +inf or -inf: on neither side of the pole, so
        a stable remainder does not certify the channel below it, while a
        remainder past the pole still is past it."""
        from hwave.solver import bond_channels as bc
        for mus in ([-1.0, 0.5], [-1.0]):
            lam, q, max_re = self._inv(mus)
            self.assertTrue(np.isnan(lam) and q is None and np.isnan(max_re), mus)
        lam, q, _ = self._inv([-1.0, -2.0])          # mu = -2 is lambda = 2
        self.assertAlmostEqual(lam, 2.0, delta=1e-12)
        self.assertEqual(q, (0.0, 0.0, 0.0))
        with _captured() as rec:
            rep = bc.check_rpa_pole({"spin": self._inv([-1.0, 0.5])}, where="t",
                                    consequence=".", log=logging.getLogger("t"))
        self.assertFalse(rep.past)
        self.assertEqual(rep.npz_members(), {})
        self.assertEqual(len([m for m in _warnings(rec) if "could not be computed" in m]), 1)

    def test_complex_pairs_on_either_side_of_real_tol(self):
        from hwave.solver import bond_channels as bc
        from hwave.solver import flex_diagnostics as diag
        t = 1.2
        for frac, real in ((0.5, True), (2.0, False)):
            d = frac * diag.REAL_TOL * t
            pair = [t + 1j * d, t - 1j * d]
            # direct: the pair is the spectrum of -sign chibar W
            lam, _, _ = bc._pole_reduce(np.array([pair]), lambda e: e, (1, 1, 1))
            self.assertEqual(np.isfinite(lam), real, (frac, lam))
            # inverted: mu chosen so that mu / (1 + mu) is that pair
            mus = [x / (1.0 - x) for x in pair]
            lam, _, _ = self._inv(mus)
            self.assertEqual(np.isfinite(lam), real, (frac, lam))
            if real:
                self.assertAlmostEqual(lam, t, delta=1e-9)


class TestReportedQAttainsTheValue(unittest.TestCase):
    """The q tie never pairs a q below the pole with a maximum past it: the
    reported value is the reported q's own."""

    def test_tie_across_the_threshold(self):
        from hwave.solver import bond_channels as bc
        ev = np.array([[1.0 - 4e-10], [1.0 + 1e-10]], dtype=complex)
        lam, q, _ = bc._pole_reduce(ev, lambda e: e, (2, 1, 1))
        self.assertEqual(lam, 1.0 + 1e-10)
        self.assertEqual(q, (0.5, 0.0, 0.0))

    def test_value_is_the_q_own(self):
        from hwave.solver import bond_channels as bc
        ev = np.array([[0.7 - 4e-10], [0.7]], dtype=complex)
        lam, q, _ = bc._pole_reduce(ev, lambda e: e, (2, 1, 1))
        self.assertEqual(q, (0.0, 0.0, 0.0))
        self.assertEqual(lam, 0.7 - 4e-10)


class TestScfContextLabel(unittest.TestCase):
    """A failed post-SCF step quotes the last iteration's indicator under the
    path's own name: the bond path's is no Stoner factor."""

    def test_label(self):
        from hwave.solver import flex_diagnostics as diag
        m = diag.ScfMonitor(logging.getLogger("t"))
        m.iteration, m.step = 3, 0.1
        m.stoner = {"spin": (0.5, (0.0, 0.0, 0.0), 0.5)}
        self.assertIn("and Stoner factor ", m.context())
        m.context_label = "RPA pole indicator (final map)"
        self.assertIn("and RPA pole indicator (final map) ", m.context())
        self.assertNotIn("Stoner", m.context())


class TestOneWarningWithAnUnavailableChannel(unittest.TestCase):
    """One channel past the pole and the other unavailable: ONE WARNING that
    names both, and keys only for the channel that was computed."""

    def test_mixed(self):
        from hwave.solver import bond_channels as bc
        with _captured() as rec:
            rep = bc.check_rpa_pole({"spin": bc.POLE_UNAVAILABLE,
                                     "charge": (1.2, (0.5, 0.5, 0.0), 1.2)},
                                    where="t", consequence=".", log=logging.getLogger("t"))
        msgs = _warnings(rec)
        self.assertEqual(len(msgs), 1, msgs)
        self.assertIn(POLE_TEXT, msgs[0])
        self.assertIn("could not be computed", msgs[0])
        self.assertIn("spin channel(s)", msgs[0])
        self.assertTrue(rep.past)
        self.assertEqual(set(rep.npz_members()),
                         {"rpa_past_pole", "rpa_pole_indicator_charge", "rpa_pole_q_charge"})


class TestPackedSectors(unittest.TestCase):
    """The S_z-packed check takes every sector (a maximum in a +-1 sector
    is found) and treats one failing sector as undetermined."""

    def _run(self, per_sector):
        import types
        from hwave import sc
        from hwave.solver import bond_channels as bc
        from hwave.solver import so_bond
        sm = so_bond.sz_sector_map(types.SimpleNamespace(reverse=np.array([0])), 2)
        X = so_bond.PackedArrays.zeros(sm, (1,))
        G = so_bond.PackedArrays.zeros(sm, (1,))
        vals = iter(per_sector[name] for name in sm.sector_order)
        with mock.patch.object(bc, "static_pole_indicator",
                               side_effect=lambda *a, **k: next(vals)), \
                _captured() as rec:
            rep = sc._so_bond_static_pole(X, G, (1, 1, 1))
        return rep, rec

    def test_maximum_in_a_transverse_sector(self):
        q = (0.0, 0.0, 0.0)
        for name in ("+1", "-1"):
            per = {"0": (0.3, q, 0.3), "+1": (0.2, q, 0.2), "-1": (0.1, q, 0.1)}
            per[name] = (1.4, q, 1.4)
            rep, rec = self._run(per)
            self.assertTrue(rep.past, name)
            self.assertEqual(rep.indicators["spinful"][0], 1.4)
            msgs = _pole_warnings(rec)
            self.assertEqual(len(msgs), 1)
            self.assertIn("S_z sector " + name, msgs[0])

    def test_one_failing_sector(self):
        from hwave.solver import bond_channels as bc
        q = (0.0, 0.0, 0.0)
        rep, rec = self._run({"0": (0.3, q, 0.3), "+1": bc.POLE_UNAVAILABLE,
                              "-1": (0.1, q, 0.1)})
        self.assertFalse(rep.past)
        self.assertEqual(len([m for m in _warnings(rec) if "could not be computed" in m]), 1)
        rep, rec = self._run({"0": (0.3, q, 0.3), "+1": bc.POLE_UNAVAILABLE,
                              "-1": (1.1, q, 1.1)})
        self.assertTrue(rep.past)
        self.assertEqual(len(_warnings(rec)), 1)

    def test_no_real_eigenvalue_is_not_a_failure(self):
        rep, rec = self._run({"0": (float("nan"), None, 0.2),
                              "+1": (float("nan"), None, 0.4),
                              "-1": (float("nan"), None, 0.1)})
        self.assertFalse(rep.past)
        self.assertEqual(_warnings(rec), [])


class TestDiagnosticFailureDoesNotEndARun(unittest.TestCase):
    """A failing eigenvalue solve of the diagnostic leaves the run complete
    and its numbers unchanged: one WARNING, no pole provenance."""

    UNAVAILABLE = "could not be computed"

    def _one_notice(self, records):
        msgs = [m for m in _warnings(records) if self.UNAVAILABLE in m]
        self.assertEqual(len(msgs), 1, msgs)
        self.assertEqual(_pole_warnings(records), [])

    def test_unit(self):
        from hwave.solver import bond_channels as bc
        chibar = np.eye(2, dtype=complex)[None]
        with _failing_diagnostic():
            got = bc.static_pole_indicator(chibar, chibar, "spin", (1, 1, 1))
            got_d = bc.pole_indicator_from_dressed(chibar, chibar, "spin", (1, 1, 1))
        for g in (got, got_d):
            self.assertTrue(np.isnan(g[0]) and g[1] is None and np.isnan(g[2]))

    def test_nonso_static(self):
        text, rec, _ = _run_nonso_static(5.5085, 0.1, failing=True)
        text0, _, _ = _run_nonso_static(5.5085, 0.1, stub=True)
        self._one_notice(rec)
        self.assertEqual(text, text0)

    def _so(self, soc, V, N):
        r, _ = _run_so(soc, V, N, failing=True)
        r0, _ = _run_so(soc, V, N, stub=True)
        self._one_notice(r.records)
        self.assertEqual(r.value, r0.value)
        self.assertEqual(set(r.gap_npz), set(r0.gap_npz))
        self.assertEqual(_stable_text(r.eigenvalue_dat), _stable_text(r0.eigenvalue_dat))

    def test_so_dense(self):
        self._so(True, 0.60, 6)

    def test_so_packed(self):
        self._so(False, 0.60, 4)

    def test_flex_final_map_and_archive_reader(self):
        args = _FLEX_FIXTURES["charge_above"][0]
        case = _flex_case("charge_above", *args, failing=True)
        ref = _flex_case("charge_above", *args, stub=True)
        self._one_notice(case.flex_records)
        self._one_notice(case.sc_records)
        self.assertEqual(float(case.pair["eigenvalue"]), float(ref.pair["eigenvalue"]))
        self.assertEqual(case.sc_value, ref.sc_value)
        for got, want in ((case.arch, ref.arch), (case.chiq_s, ref.chiq_s),
                          (case.pair, ref.pair), (case.dyn, ref.dyn)):
            self.assertEqual(set(got), set(want))
        self.assertEqual(case.pair_dat, ref.pair_dat)
        self.assertEqual(case.dyn_dat, ref.dyn_dat)


class TestHalfGridStaticIndex(unittest.TestCase):
    """Under longitudinal_bond_half_grid the static slice is stored row 0;
    the final map's indicators equal the full grid's."""

    def test_half_equals_full(self):
        args = _FLEX_FIXTURES["charge_above"][0]
        full = _flex_case("charge_above", *args)
        half = _flex_case("charge_above", *args, extra={"longitudinal_bond_half_grid": True})
        self.assertTrue(bool(half.arch["longitudinal_bond_half_grid"]))
        for ch in ("spin", "charge"):
            a = float(full.arch[_PREFIX + "rpa_pole_indicator_" + ch])
            b = float(half.arch[_PREFIX + "rpa_pole_indicator_" + ch])
            self.assertAlmostEqual(a, b, delta=1e-12 * max(1.0, abs(a)))
            np.testing.assert_array_equal(full.arch[_PREFIX + "rpa_pole_q_" + ch],
                                          half.arch[_PREFIX + "rpa_pole_q_" + ch])
        self.assertTrue(bool(half.arch[_PREFIX + "rpa_past_pole"]))


if __name__ == "__main__":
    unittest.main()
