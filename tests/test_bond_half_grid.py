"""The half bosonic grid of the bond-resolved FLEX store (issue #223).

Every bosonic object of the bond map satisfies, at the SAME q,
``X(q, l) = X(q, nmat - l)^dagger`` on the ``ND x ND`` pair matrix, so with
``[mode.param] longitudinal_bond_half_grid = true`` the block store keeps
only the ``nmat // 2 + 1`` slices ``l(h) = (h + nmat // 2) mod nmat`` (the
static slice ``h = 0`` and the Nyquist slot ``l = 0`` included) and the
transport, the collapses and the published archives rebuild the rest.

The oracle is the full-grid solve: on every fixture below the half-grid
solve must reproduce the self-energy, the Green function, mu, the channel-0
collapses, every ``longitudinal_bond_*`` member and every archive member
to 1e-12 relative (the regrouped frequency batches move the round-off).
The guard statistics count the solved slices only and are compared
against that count instead. The mutation tests break the bookkeeping one
way at a time (the Nyquist slot, the transpose of the dagger, the pair
swap, the partner index) and require the oracle to fail.
"""
import contextlib
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np

from hwave.solver import flex_bond

_IN2 = "tests/rpa/input_2orb"
_TOL = 1.0e-12

# interaction tables (Rx, Ry, Rz, a, b, re, im) of the two-orbital fixture
_U = [(0, 0, 0, 1, 1, 1.5, 0.0), (0, 0, 0, 2, 2, 1.5, 0.0)]
_HUND = [(0, 0, 0, 1, 2, 0.3, 0.0), (0, 0, 0, 2, 1, 0.3, 0.0)]
# inter-orbital, NOT symmetric under R -> -R
_VASYM = [(1, 0, 0, 1, 2, 0.5, 0.0), (-1, 0, 0, 2, 1, 0.5, 0.0),
          (-1, 0, 0, 1, 2, 0.2, 0.0), (1, 0, 0, 2, 1, 0.2, 0.0),
          (0, 1, 0, 1, 1, 0.4, 0.0), (0, -1, 0, 1, 1, 0.4, 0.0)]
_EXCH = [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0)]
_PH = [(0, 1, 0, 1, 1, 0.25, 0.0), (0, -1, 0, 1, 1, 0.25, 0.0)]
# Hermitian but time-reversal breaking hopping (t(k)* != t(-k))
_TCPLX = [(0, 1, 0, 1, 1, 1.0, 0.4), (0, -1, 0, 1, 1, 1.0, -0.4),
          (0, 1, 0, 2, 2, 1.0, -0.3), (0, -1, 0, 2, 2, 1.0, 0.3),
          (0, 0, 0, 1, 2, 0.5, 0.3), (0, 0, 0, 2, 1, 0.5, -0.3),
          (-1, 0, 0, 1, 2, 0.5, -0.2), (1, 0, 0, 2, 1, 0.5, 0.2)]

#: the fixtures: None = the committed input_2orb CoulombInter table
_CASES = {
    "UV": dict(inter=None),
    "U+Hund+V_tail": dict(inter={"CoulombIntra": _U, "Hund": _HUND, "CoulombInter": _VASYM},
                          extra={"coeff_tail": 1.0}),
    "Exchange+PairHop": dict(inter={"CoulombIntra": _U, "CoulombInter": _VASYM,
                                    "Exchange": _EXCH, "PairHop": _PH}),
    "complex_t": dict(inter={"CoulombIntra": _U, "CoulombInter": _VASYM}, transfer=_TCPLX),
    "complex_t+Exchange+PairHop": dict(inter={"CoulombIntra": _U, "CoulombInter": _VASYM,
                                              "Exchange": _EXCH, "PairHop": _PH},
                                       transfer=_TCPLX),
}


def _table(name, rows, norb=2):
    Rs = sorted({r[:3] for r in rows})
    s = "{} half-grid test\n{}\n{}\n".format(name, norb, len(Rs))
    s += " " + " ".join(["1"] * len(Rs)) + "\n"
    for r in rows:
        s += "  {:d} {:d} {:d} {:d} {:d} {:.12f} {:.12f}\n".format(*r)
    return s


@contextlib.contextmanager
def _input(case):
    """A reader of the committed two-orbital geometry with the case's
    tables (written to a temporary directory when not the committed one)."""
    import hwave.qlmsio.read_input_k as read_input_k
    spec = _CASES[case]
    if spec["inter"] is None and spec.get("transfer") is None:
        idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat"}
        yield read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
        return
    d = tempfile.mkdtemp()
    try:
        shutil.copy(os.path.join(_IN2, "geom.dat"), d)
        if spec.get("transfer") is None:
            shutil.copy(os.path.join(_IN2, "transfer.dat"), d)
        else:
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write(_table("Transfer", spec["transfer"]))
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for k, rows in spec["inter"].items():
            fn = k.lower() + ".dat"
            with open(os.path.join(d, fn), "w") as f:
                f.write(_table(k, rows))
            idict[k] = fn
        yield read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    finally:
        shutil.rmtree(d, ignore_errors=True)


_BASE = {"T": 1.0, "filling": 0.45, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
         "Nmat": 8, "IterationMax": 3, "Mix": 0.3, "EPS": 1e-14,
         "flex_hartree_fock": True, "longitudinal_bond_channels": True,
         "mixing_scheme": "linear", "longitudinal_bond_freq_batch": 3}

_OUTFILES = {"sigma": "sigma.npz", "green": "green.npz", "chi0q": "chi0q.npz",
             "longitudinal_bond": "longitudinal_bond.npz"}


def _solver(case, extra, half):
    import hwave.solver.flex as flex_mod
    par = dict(_BASE)
    par.update(_CASES[case].get("extra", {}))
    par.update(extra or {})
    if half is not None:
        par["longitudinal_bond_half_grid"] = half
    par = {k: v for k, v in par.items() if v is not None}
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
            "calc_scheme": "general"}
    return info, flex_mod.FLEX


def _run(case, extra=None, half=True, interval=False, save=True):
    """One gated solve; returns ``(solver, green_info, {file: {member: array}})``."""
    from hwave.solver import bond_channels as bc
    info, FLEX = _solver(case, extra, half)
    with _input(case) as r:
        s = FLEX(r.get_param("ham"), {}, info)
        gi = r.get_param("green")
        files = {}
        with contextlib.ExitStack() as st:
            if interval:
                # the resolution "auto" makes on a device (issue #276 tests)
                st.enter_context(mock.patch.dict(bc._GUARD_AUTO, {"numpy": "interval"}))
            out = st.enter_context(tempfile.TemporaryDirectory())
            s.solve(gi, out)
            if save:
                s.save_results(dict(_OUTFILES, path_to_output=out), gi)
                for fn in sorted(os.listdir(out)):
                    if fn.endswith(".npz"):
                        with np.load(os.path.join(out, fn), allow_pickle=False) as z:
                            files[fn] = {m: np.array(z[m]) for m in z.files}
    return s, gi, files


#: members that legitimately differ: the guard statistics count the solved
#: slices (checked separately), the provenance key exists on the half run only
_COUNTS = ("longitudinal_bond_guard_exact_blocks_total",
           "longitudinal_bond_guard_exact_blocks_max",
           "longitudinal_bond_guard_blocks_per_iteration")


def _is_roundoff_diagnostic(name):
    """Members whose VALUE is a round-off measurement (the Hermiticity error
    of the density, the SCF residuals): compared absolutely."""
    leaf = name.split(":")[-1]
    return leaf.endswith("_error") or leaf.endswith("_residual")


def _rel(a, b, name=""):
    a, b = np.asarray(a), np.asarray(b)
    if _is_roundoff_diagnostic(name) and a.shape == b.shape and a.dtype.kind in "fc":
        return float(np.abs(a - b).max()) if a.size else 0.0
    if a.shape != b.shape:
        return np.inf
    if a.size == 0:
        return 0.0
    if a.dtype.kind not in "fc":
        return 0.0 if np.array_equal(a, b) else np.inf
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(b))):
        return np.inf
    scale = max(float(np.abs(a).max()), 1e-300)
    return float(np.abs(a - b).max()) / scale


def _compare(full, half):
    """``{name: relative deviation}`` over everything the two runs publish."""
    sf, gf, ff = full
    sh, gh, fh = half
    dev = {"mu": abs(float(sf.mu) - float(sh.mu)) / max(1.0, abs(float(sf.mu)))}
    keys = set(gf) | set(gh)
    for k in sorted(keys):
        if k in _COUNTS or k in ("longitudinal_bond_half_grid", "longitudinal_bond_nb"):
            continue
        a, b = gf.get(k), gh.get(k)
        if not isinstance(a, (np.ndarray, np.generic)) and not isinstance(
                b, (np.ndarray, np.generic)):
            continue
        if a is None or b is None:
            dev["green_info:" + k] = np.inf
            continue
        dev["green_info:" + k] = _rel(a, b, k)
    if set(ff) != set(fh):
        dev["files"] = np.inf
    for fn in ff:
        ma, mb = ff[fn], fh.get(fn, {})
        extra = set(mb) - set(ma)
        if extra - {"longitudinal_bond_half_grid"}:
            dev[fn + ":extra"] = np.inf
        for m in ma:
            if m in _COUNTS or m == "longitudinal_bond_nb":
                continue
            if m not in mb:
                dev[fn + ":" + m] = np.inf
                continue
            dev[fn + ":" + m] = _rel(ma[m], mb[m], m)
    return dev


class _OracleCase(unittest.TestCase):

    def _assert_equal_runs(self, case, extra=None, interval=False, save=True):
        full = _run(case, extra, half=False, interval=interval, save=save)
        half = _run(case, extra, half=True, interval=interval, save=save)
        dev = _compare(full, half)
        bad = {k: v for k, v in dev.items() if not v <= _TOL}
        self.assertEqual(bad, {}, "full vs half deviations above {}: {}".format(_TOL, bad))
        sf, gf, ff = full
        sh, gh, fh = half
        # the oracle really compared the dynamic content
        self.assertIn("sigma", gf)
        self.assertIn("chiq_s", gf)
        self.assertIn("longitudinal_bond_chi_s", gf)
        nmat = int(sf.nmat)
        ns = nmat // 2 + 1
        # the published batch: the pinned one, clamped to the stored axis
        self.assertEqual(int(gh["longitudinal_bond_nb"]), min(int(gf["longitudinal_bond_nb"]), ns))
        # the half run stored nmat/2+1 slices, published the provenance key
        # and the full-grid archive
        self.assertTrue(sh.longitudinal_bond_half_grid)
        # the mirrored slices are guarded on their own when time reversal is
        # broken (complex hopping), and always under "warn" or "static"
        self.assertIs(sh._bond_last.mirror_guarded,
                      case.startswith("complex_t") or sf.longitudinal_bond_guard_freqs != "all"
                      or sf.flex_guard_policy != "refuse")
        self.assertNotIn("longitudinal_bond_half_grid", gf)
        self.assertTrue(bool(gh["longitudinal_bond_half_grid"]))
        if save:
            self.assertTrue(bool(fh["chiq_s.npz"]["longitudinal_bond_half_grid"]))
        if "longitudinal_bond.npz" in ff:
            self.assertEqual(fh["longitudinal_bond.npz"]["chi_s_w"].shape[0], nmat)
            self.assertTrue(bool(fh["longitudinal_bond.npz"]["longitudinal_bond_half_grid"]))
        # cond_min where defined
        for ch in ("s", "c"):
            a = float(gf["longitudinal_bond_cond_min_" + ch])
            b = float(gh["longitudinal_bond_cond_min_" + ch])
            if np.isfinite(a):
                self.assertLessEqual(abs(a - b), _TOL * max(1.0, abs(a)), ch)
            else:
                self.assertFalse(np.isfinite(b), ch)
        # the guard statistics count the solved slices: ns / nmat of the
        # full grid's under "all", the same single static slice under "static"
        gb_f = int(gf["longitudinal_bond_guard_blocks_per_iteration"])
        gb_h = int(gh["longitudinal_bond_guard_blocks_per_iteration"])
        if sf.longitudinal_bond_guard_freqs == "all":
            self.assertEqual(gb_h * nmat, gb_f * ns)
        else:
            self.assertEqual(gb_h, gb_f)
        return dev


class TestIndexMaps(unittest.TestCase):

    def test_grid_index_and_slice_count(self):
        from hwave.solver import half_grid as hg
        self.assertEqual(hg.nslices(2), 2)
        self.assertEqual(hg.nslices(8), 5)
        np.testing.assert_array_equal(hg.grid_index(2), [1, 0])
        np.testing.assert_array_equal(hg.grid_index(4), [2, 3, 0])
        np.testing.assert_array_equal(hg.grid_index(8), [4, 5, 6, 7, 0])
        np.testing.assert_array_equal(hg.grid_index(8, 3, 5), [7, 0])
        # nu_h = 2 pi T h: the grid frequency 2 l - nmat is 2 h (mod 2 nmat)
        ls = hg.grid_index(16)
        np.testing.assert_array_equal((2 * ls - 16) % 32, (2 * np.arange(9)) % 32)
        for bad in (0, 3, 7, -2):
            with self.assertRaises(ValueError):
                hg.nslices(bad)


def _hermitian_in_frequency(nmat, nvol, n, seed=0):
    """A random full-grid stack with X(l) = X(nmat - l)^dagger (endpoints
    l = 0 and nmat/2 Hermitian)."""
    rng = np.random.default_rng(seed)
    h0 = nmat // 2
    x = (rng.standard_normal((nmat, nvol, n, n))
         + 1j * rng.standard_normal((nmat, nvol, n, n)))
    for l in (0, h0):
        x[l] = 0.5 * (x[l] + np.conj(x[l]).swapaxes(-1, -2))
    for l in range(1, h0):
        x[l] = np.conj(x[nmat - l]).swapaxes(-1, -2)
    return x


class TestStore(unittest.TestCase):

    def _check(self, nmat, nvol=3, nd=4, B=3):
        from hwave.solver import half_grid as hg
        ND = B * nd
        full = _hermitian_in_frequency(nmat, nvol, ND, seed=nmat)
        np.testing.assert_array_equal(hg.expand(hg.half_from_full(full, nmat), nmat), full)
        with flex_bond.BondBlockStore(nmat, nvol, ND, nd, ("chibar",), {"W": "chibar"},
                                      half_grid=True) as st, \
                flex_bond.BondBlockStore(nmat, nvol, ND, nd, ("chibar",),
                                         {"W": "chibar"}) as ref:
            self.assertEqual(st.nslices, nmat // 2 + 1)
            self.assertEqual(st.get_freq_batch("chibar", 0, st.nslices).shape[0],
                             nmat // 2 + 1)
            np.testing.assert_array_equal(st.grid_index(0, st.nslices), hg.grid_index(nmat))
            np.testing.assert_array_equal(ref.grid_index(0, nmat), np.arange(nmat))
            for a in range(B):
                for b in range(B):
                    blk = np.ascontiguousarray(full[:, :, a * nd:(a + 1) * nd,
                                                    b * nd:(b + 1) * nd])
                    st.put_pair_fullgrid("chibar", a, b, blk)
                    ref.put_pair_fullgrid("chibar", a, b, blk)
            for a in range(B):
                for b in range(B):
                    want = full[:, :, a * nd:(a + 1) * nd, b * nd:(b + 1) * nd]
                    got = st.get_pair_full("W", a, b, xp=np)
                    self.assertTrue(got.flags.c_contiguous)
                    np.testing.assert_array_equal(got, want)
                    got_ref = ref.get_pair_full("W", a, b, xp=np)
                    np.testing.assert_array_equal(got_ref, want)
                    np.testing.assert_array_equal(ref.get_pair_full("W", a, b), want)
            np.testing.assert_array_equal(st.expand_full("chibar"), full)
            out = np.zeros_like(full)
            self.assertIs(st.expand_full("chibar", out=out, chunk=1), out)
            np.testing.assert_array_equal(out, full)
            # on the full grid expand_full is the stored array itself
            self.assertIs(ref.expand_full("chibar"), ref._arrays["chibar"])
            np.testing.assert_array_equal(ref.expand_full("chibar"), full)

    def test_round_trip(self):
        for nmat in (2, 4, 8, 16):
            with self.subTest(nmat=nmat):
                self._check(nmat)

    def test_nmat_two_has_no_mirrored_slice(self):
        from hwave.solver import half_grid as hg
        x = _hermitian_in_frequency(2, 1, 2)
        half = hg.half_from_full(x, 2)
        np.testing.assert_array_equal(half, x[[1, 0]])

    def test_fill_mirror_and_expand_with_transposed_pair(self):
        from hwave.solver import half_grid as hg
        full = _hermitian_in_frequency(8, 2, 6, seed=3)
        part = full.copy()
        part[1:4] = np.nan
        np.testing.assert_array_equal(hg.fill_mirror(part, 8), full)
        # pair (0, 1) of a 2 x 2 block layout reads the mirror of pair (1, 0)
        ab = hg.half_from_full(full[:, :, :3, 3:], 8)
        ba = hg.half_from_full(full[:, :, 3:, :3], 8)
        np.testing.assert_array_equal(hg.expand(ab, 8, ba=ba), full[:, :, :3, 3:])
        self.assertFalse(np.array_equal(hg.expand(ab, 8), full[:, :, :3, 3:]))

    def test_shape_errors(self):
        from hwave.solver import half_grid as hg
        with self.assertRaises(ValueError):
            hg.expand(np.zeros((4, 1, 2, 2)), 8)
        with self.assertRaises(ValueError):
            hg.half_from_full(np.zeros((6, 1, 2, 2)), 8)
        with self.assertRaises(ValueError):
            flex_bond.BondBlockStore(7, 1, 2, 1, ("chibar",), half_grid=True)


class TestDressBatchGridIndices(unittest.TestCase):
    """``dress_batch(ls=...)`` reports the TRUE grid index of a half-grid batch."""

    def _batch(self, nb=3, nvol=2, ND=2):
        cb = np.zeros((nb, nvol, ND, ND), dtype=complex)
        W = np.zeros((nvol, ND, ND), dtype=complex)
        for i in range(ND):
            W[:, i, i] = 1.0
        return cb, W

    def test_refusal_names_the_true_l(self):
        from hwave.solver import bond_channels as bc
        nmat = 8
        cb, W = self._batch()
        # the charge denominator 1 + chibar W is singular in row 1 at q = 1
        cb[1, 1] = -np.eye(2)
        ls = np.array([7, 0, 1])            # wraps through the Nyquist slot
        with self.assertRaises(bc.BondConditioningError) as cm:
            bc.dress_batch(cb, W, "charge", l0=7, nmat=nmat, spatial_shape=(2, 1, 1), ls=ls)
        self.assertEqual(cm.exception.l, 0)
        self.assertIn("grid index l=0", str(cm.exception))
        self.assertIn("bosonic Matsubara index {}".format(-nmat), str(cm.exception))

    def test_static_guard_row_is_found_by_grid_index(self):
        from hwave.solver import bond_channels as bc
        cb, W = self._batch()
        cb[2, 0] = -np.eye(2)       # the row whose l is nmat // 2 = 4
        with self.assertRaises(bc.BondConditioningError) as cm:
            bc.dress_batch(cb, W, "charge", l0=6, nmat=8, spatial_shape=(2, 1, 1),
                           guard_freqs="static", ls=np.array([6, 7, 4]))
        self.assertEqual(cm.exception.l, 4)
        # no row at the static frequency: nothing is conditioning-checked
        chi, cond = bc.dress_batch(cb * 0, W, "charge", l0=5, nmat=8, spatial_shape=(2, 1, 1),
                                   guard_freqs="static", ls=np.array([5, 6, 7]))
        self.assertIsNone(cond)

    def test_length_mismatch_is_refused(self):
        from hwave.solver import bond_channels as bc
        cb, W = self._batch()
        with self.assertRaisesRegex(ValueError, "ls must give"):
            bc.dress_batch(cb, W, "spin", l0=0, nmat=8, spatial_shape=(2, 1, 1),
                           ls=np.array([4, 5]))


class TestConfig(unittest.TestCase):

    def _raw(self, param, gate=True, mode_extra=None):
        from hwave.solver.flex import FLEX
        p = {"flex_hartree_fock": True, "longitudinal_bond_channels": gate, "Nmat": 8}
        p.update(param)
        info = {"param": p, "calc_scheme": "general"}
        info.update(mode_extra or {})
        return FLEX._parse_phase_b_keys(info)

    def test_default_false_and_case_insensitive_key(self):
        self.assertFalse(self._raw({})["longitudinal_bond_half_grid"])
        self.assertTrue(self._raw({"longitudinal_bond_half_grid": True})
                        ["longitudinal_bond_half_grid"])
        self.assertTrue(self._raw({"Longitudinal_Bond_Half_Grid": True})
                        ["longitudinal_bond_half_grid"])
        self.assertTrue(self._raw({"longitudinal_bond_half_grid": np.bool_(True)})
                        ["longitudinal_bond_half_grid"])

    def test_only_booleans(self):
        for v in ("true", 1, 0, None, 1.0):
            with self.subTest(v=v):
                with self.assertRaisesRegex(ValueError, "longitudinal_bond_half_grid must be a "
                                                        "boolean"):
                    self._raw({"longitudinal_bond_half_grid": v})

    def test_stale_without_the_gate(self):
        from hwave.solver.flex import FLEX
        self.assertIn("longitudinal_bond_half_grid", FLEX._PHASE_B_BOND_KEYS)
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            out = self._raw({"longitudinal_bond_half_grid": "garbage"}, gate=False)
        self.assertFalse(out["longitudinal_bond_half_grid"])
        self.assertIn("longitudinal_bond_half_grid", "\n".join(cm.output))
        # neither switch on: the same stale warning
        from hwave.solver.flex import FLEX as F
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            out = F._parse_phase_b_keys({"param": {"longitudinal_bond_half_grid": True}})
        self.assertFalse(out["longitudinal_bond_half_grid"])
        self.assertIn("longitudinal_bond_half_grid", "\n".join(cm.output))

    def test_refused_with_pairing(self):
        for p in ("singlet", "Triplet", "both"):
            with self.subTest(pairing=p):
                with self.assertRaisesRegex(ValueError, "does not support the in-process "
                                                        "pairing yet"):
                    self._raw({"longitudinal_bond_half_grid": True,
                               "longitudinal_bond_pairing": p})
        # pairing "none" and a half grid of false are accepted
        self._raw({"longitudinal_bond_half_grid": True, "longitudinal_bond_pairing": "None"})
        self._raw({"longitudinal_bond_half_grid": False, "longitudinal_bond_pairing": "singlet",
                   "IterationMax": 1})

    def test_refused_with_spin_orbital(self):
        for key in ("enable_spin_orbital", "Enable_Spin_Orbital"):
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, "enable_spin_orbital"):
                    self._raw({"longitudinal_bond_half_grid": True}, mode_extra={key: True})
        self._raw({"longitudinal_bond_half_grid": True},
                  mode_extra={"enable_spin_orbital": False})

    def test_rpa_warns_it_is_flex_only(self):
        from tests.test_flex_bond_gate import _rpa
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            _rpa({"longitudinal_bond_half_grid": True})
        self.assertTrue(any("longitudinal_bond_half_grid" in m and "FLEX-only" in m
                            for m in cm.output), cm.output)


class TestNonHermitianVertexRefused(unittest.TestCase):
    """The input layer cannot produce a non-Hermitian S/C: every interaction
    table must be Hermitian-closed (the reader refuses ``X_ab(R) !=
    conj(X_ba(-R))``), off-site coefficients must be real, and a complex
    Hermitian-closed on-site table gives Hermitian vertices (or is refused by
    the Hartree-Fock map). The refusal is therefore exercised on a vertex
    built non-Hermitian behind the reader's back."""

    def test_direct_check(self):
        from hwave.solver.flex import FLEX
        S = np.zeros((2, 3, 3), dtype=complex)
        S[1, 0, 1] = 1.0
        S[1, 1, 0] = 1.0
        FLEX._check_bond_vertices_hermitian(S, S)
        C = S.copy()
        C[1, 0, 1] = 1.0 + 1e-9j
        with self.assertRaisesRegex(ValueError, r"needs Hermitian bond vertices.*C\(q\)"):
            FLEX._check_bond_vertices_hermitian(S, C)

    def test_solve_refuses_before_the_first_map(self):
        from hwave.solver import bond_channels as bc
        real = bc.build_sc_bond_channel

        def skewed(topo, V0, label, **kw):
            out = real(topo, V0, label, **kw)
            if label == "S":
                out = out.copy()
                out[0, 0, -1] += 1e-6j
            return out

        for half in (True, False):
            with self.subTest(half=half):
                with mock.patch.object(bc, "build_sc_bond_channel", skewed):
                    if half:
                        with self.assertRaisesRegex(ValueError, "needs Hermitian bond "
                                                                "vertices"):
                            _run("UV", {"IterationMax": 1}, half=True, save=False)
                    else:
                        # the full grid does not need the symmetry
                        _run("UV", {"IterationMax": 1}, half=False, save=False)


class TestBitIdentityWhenOff(unittest.TestCase):
    """``longitudinal_bond_half_grid = false`` (and the key absent) never
    reaches the half-grid code: every helper of :mod:`half_grid` raises if
    called, and the two runs agree.

    The two runs execute the same code, so they are compared to round-off,
    not bitwise: threaded FFT/BLAS reductions are not guaranteed to be
    bit-reproducible between two runs in one process (observed at 1e-16
    under pytest with coverage in CI).  Bit identity of the key-false path
    against the previous revision is checked separately on archived
    outputs."""

    def test_off_never_touches_the_half_grid(self):
        from hwave.solver import half_grid as hg

        def boom(*a, **k):
            raise AssertionError("half-grid helper called on a full-grid run")

        with mock.patch.multiple(hg, expand=boom, fill_mirror=boom, half_from_full=boom,
                                 grid_index=boom, describe_slices=boom):
            absent = _run("UV", {"longitudinal_bond_output_full": True}, half=None)
            off = _run("UV", {"longitudinal_bond_output_full": True}, half=False)
        sa, ga, fa = absent
        so, go, fo = off
        self.assertFalse(so.longitudinal_bond_half_grid)
        def same(a, b, msg):
            a, b = np.asarray(a), np.asarray(b)
            if np.issubdtype(a.dtype, np.inexact):
                np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-14, err_msg=msg)
            else:
                np.testing.assert_array_equal(a, b, err_msg=msg)

        self.assertEqual(sorted(ga), sorted(go))
        for k in ga:
            if isinstance(ga[k], np.ndarray):
                same(ga[k], go[k], k)
        self.assertEqual(sorted(fa), sorted(fo))
        for fn in fa:
            self.assertEqual(sorted(fa[fn]), sorted(fo[fn]))
            self.assertNotIn("longitudinal_bond_half_grid", fa[fn])
            for m in fa[fn]:
                same(fa[fn][m], fo[fn][m], fn + ":" + m)
        same(float(sa.mu), float(so.mu), "mu")


class TestMemory(unittest.TestCase):

    _KW = dict(nmat=64, nvol=16, norb=2, B=5, depth=0, split_seed=False, n_types=1,
               freq_batch=None, cap_gb=1e6, mixing="linear", w_in_place=True)

    def test_store_rows_scale(self):
        for of in (False, True):
            with self.subTest(output_full=of):
                full = flex_bond.estimate_bond_memory(output_full=of, **self._KW)
                half = flex_bond.estimate_bond_memory(output_full=of, half_grid=True, **self._KW)
                U = 64 * 16 * (5 * 4) ** 2 * 16
                Uh = 33 * 16 * (5 * 4) ** 2 * 16
                self.assertEqual(full["persistent_rows"]["chibar_W"], U)
                self.assertEqual(half["persistent_rows"]["chibar_W"], Uh)
                self.assertEqual(half["persistent_rows"]["chi_sc_w"], 2 * Uh if of else 0)
                self.assertEqual(half["nslices"], 33)
                self.assertEqual(full["nslices"], 64)
                self.assertLessEqual(half["nb"], 33)
                self.assertEqual(half["dressing_ops"] * 64, full["dressing_ops"] * 33)
                # the transport rebuilds full-grid blocks: unchanged row
                self.assertEqual(half["phase_rows"]["transport"],
                                 full["phase_rows"]["transport"])
                if of:
                    # the publish expansion: two full members, one half alive
                    G = 64 * 16 * 4 * 16
                    self.assertEqual(half["phase_rows"]["post_scf"], G + 2 * U - Uh)

    def test_pinned_batch_is_clamped_to_the_stored_axis(self):
        est = flex_bond.estimate_bond_memory(output_full=False, half_grid=True,
                                             **dict(self._KW, freq_batch=64))
        self.assertEqual(est["nb"], 33)
        est = flex_bond.estimate_bond_memory(output_full=False, half_grid=True,
                                             **dict(self._KW, freq_batch=5))
        self.assertEqual(est["nb"], 5)

    def test_device_rows_select_within_the_stored_axis(self):
        est = flex_bond.estimate_bond_memory(output_full=False, half_grid=True,
                                             device_available=2 ** 50, **self._KW)
        self.assertEqual(est["device_nb"], 33)
        self.assertEqual(est["nb"], 33)

    def test_solver_preflight_and_transfer_volume(self):
        from hwave.solver.flex import FLEX
        info, _ = _solver("UV", {"IterationMax": 0, "longitudinal_bond_freq_batch": None}, True)
        with _input("UV") as r:
            s = FLEX(r.get_param("ham"), {}, info)
            with tempfile.TemporaryDirectory() as out:
                s.solve(r.get_param("green"), out)
        self.assertTrue(s._bond_est_kwargs["half_grid"])
        est = s._bond_est
        self.assertEqual(est["nslices"], 5)
        self.assertLessEqual(s._bond_nb, 5)
        with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
            s._log_bond_transfer_volume(est)
        line = [m for m in cm.output if "transfer volume" in m][0]
        u = 5 * est["nvol"] * est["ND"] ** 2 * 16
        S_b, C_b, G_b = est["S_bytes"], est["C_bytes"], est["G_bytes"]
        h2d = u + est["B"] ** 2 * C_b
        d2h = u + u * 3 + 3 * C_b + 2 * S_b + G_b
        gib = 1024.0 ** 3
        self.assertIn("host to device {:.4f} GiB".format(h2d / gib), line)
        self.assertIn("device to host {:.4f} GiB".format(d2h / gib), line)

    def test_device_admission_with_a_fake_device(self):
        import types
        import hwave.solver.flex as flex_mod
        info, FLEX = _solver("UV", {"IterationMax": 0, "longitudinal_bond_freq_batch": None},
                             True)
        with _input("UV") as r:
            s = FLEX(r.get_param("ham"), {}, info)
            with tempfile.TemporaryDirectory() as out:
                s.solve(r.get_param("green"), out)
        fake_xp = types.SimpleNamespace(__name__="cupy")
        with mock.patch.object(flex_mod._bk, "device_available_bytes", lambda: 2 ** 50):
            with self.assertLogs("hwave.solver.flex", level="INFO"):
                s._bond_device_admission(fake_xp)
        self.assertEqual(s._bond_est["device_nb"], 5)
        self.assertEqual(s._bond_nb, 5)


class TestOracleMixingAndGuard(_OracleCase):
    """The committed U + V fixture across the mixing schemes, guard modes and
    output switches (Nmat 8, batch 3: the second batch holds l = 7 and the
    Nyquist slot l = 0)."""

    def test_mixing_guard_freqs_output_full(self):
        for mixing in ("linear", "anderson"):
            for gf in ("all", "static"):
                for of in (True, False):
                    with self.subTest(mixing=mixing, guard_freqs=gf, output_full=of):
                        self._assert_equal_runs("UV", {
                            "mixing_scheme": mixing, "longitudinal_bond_guard_freqs": gf,
                            "longitudinal_bond_output_full": of})

    def test_guard_off(self):
        """The disabled guard (``cond_tol = None``) is a state of the map
        functions only -- the configuration requires a floor in (0, 1) -- so
        it is compared at the map level: on every map of a real solve the
        captured bubble is dressed and transported once on a full-grid and
        once on a half-grid store, with the guard off."""
        from hwave.solver import half_grid as hg
        real_dress = flex_bond.dress_and_build_w
        real_tr = flex_bond.calc_self_energy_bond
        seen = {"maps": 0}

        def spy_dress(store, dev, **kw):
            nmat, nvol, nd = kw["nmat"], kw["nvol"], kw["nd"]
            ND = store.ND
            cb = np.array(store.get_freq_batch("chibar", 0, nmat))
            kw0 = dict(kw, cond_tol=None, output_full=True)
            names = ("chibar", "chi_s_w", "chi_c_w")
            f = flex_bond.BondBlockStore(nmat, nvol, ND, nd, names, {"W": "chibar"})
            h = flex_bond.BondBlockStore(nmat, nvol, ND, nd, names, {"W": "chibar"},
                                         half_grid=True)
            f.put_freq_batch("chibar", 0, nmat, cb)
            h.put_freq_batch("chibar", 0, h.nslices, hg.half_from_full(cb, nmat))
            rf = real_dress(f, dev, **kw0)
            rh = real_dress(h, dev, **dict(kw0, nb=2))
            self.assertEqual((rf.guard_blocks, rh.guard_blocks), (0, 0))
            for k in ("collapse0", "collapse_s", "collapse_c", "static_s", "static_c"):
                self.assertLessEqual(_rel(getattr(rf, k), getattr(rh, k)), _TOL, k)
            for name in ("W", "chi_s_w", "chi_c_w"):
                self.assertLessEqual(_rel(f.expand_full(name), h.expand_full(name)), _TOL, name)
            # kept alive for the transport of the same map (released there)
            seen["store_pair"] = (f, h)
            seen["maps"] += 1
            return real_dress(store, dev, **kw)

        def spy_tr(store, green_kw, beta, view, shape, norb, workers, xp=np):
            f, h = seen.pop("store_pair")
            with f, h:
                sf = real_tr(f, green_kw, beta, view, shape, norb, workers, xp=xp)
                sh = real_tr(h, green_kw, beta, view, shape, norb, workers, xp=xp)
            self.assertLessEqual(_rel(sf, sh), _TOL)
            seen["sigma_dev"] = max(seen.get("sigma_dev", 0.0), _rel(sf, sh))
            return real_tr(store, green_kw, beta, view, shape, norb, workers, xp=xp)

        with mock.patch.object(flex_bond, "dress_and_build_w", side_effect=spy_dress), \
                mock.patch.object(flex_bond, "calc_self_energy_bond", side_effect=spy_tr):
            _run("UV", {"IterationMax": 2}, half=False, save=False)
        self.assertEqual(seen["maps"], 2)

    def test_guard_exact_decide_and_always(self):
        # the interval guard the device resolves to, on the host (issue #276):
        # "decide" brackets the earlier maps and is forced exact on the last
        for ge in ("decide", "always"):
            for gf in ("all", "static"):
                with self.subTest(guard_exact=ge, guard_freqs=gf):
                    self._assert_equal_runs("UV", {
                        "longitudinal_bond_guard_exact": ge, "longitudinal_bond_guard_freqs": gf,
                        "longitudinal_bond_output_full": True}, interval=True)

    def test_pinned_batches(self):
        # nmat 8 (5 slices): 1, 2 (batch [4, 5) is the Nyquist slot alone),
        # 4, 5 (one batch from nu = 0 through the Nyquist slot) and 8 (> the
        # stored axis, clamped); nmat 16 (9 slices): 5 straddles l = 15 -> 0
        for nmat, nb in ((8, 1), (8, 2), (8, 4), (8, 5), (8, 8), (16, 5), (16, 9)):
            with self.subTest(nmat=nmat, nb=nb):
                self._assert_equal_runs("UV", {"Nmat": nmat, "longitudinal_bond_freq_batch": nb,
                                               "longitudinal_bond_output_full": True})


class TestOracleFixtures(_OracleCase):
    """U + Hund + V with the Green-function tail, off-site Exchange + PairHop
    (the exact bond-enlarged second order) and complex, time-reversal
    breaking hopping, at Nmat 8 and 16."""

    def _sweep(self, case):
        for nmat in (8, 16):
            for mixing in ("linear", "anderson"):
                with self.subTest(case=case, nmat=nmat, mixing=mixing):
                    self._assert_equal_runs(case, {"Nmat": nmat, "mixing_scheme": mixing,
                                                   "longitudinal_bond_output_full": True})

    def test_hund_with_tail(self):
        self._sweep("U+Hund+V_tail")

    def test_exchange_pairhop(self):
        self._sweep("Exchange+PairHop")

    def test_complex_hopping(self):
        self._sweep("complex_t")

    def test_rectangular_cell(self):
        # nx != ny pins the C-order -q map of the time-reversal check against
        # the lattice: a wrong flattening would switch the skip off for the
        # real-hopping fixture (asserted in _assert_equal_runs)
        for case in ("UV", "Exchange+PairHop", "complex_t"):
            with self.subTest(case=case):
                self._assert_equal_runs(case, {"CellShape": [4, 3, 1],
                                               "longitudinal_bond_output_full": True})

    def test_complex_hopping_with_exchange_pairhop_static_guard(self):
        # the remaining cross: complex t with the bond second order of an
        # off-site Exchange
        self._assert_equal_runs("complex_t+Exchange+PairHop", {
            "Nmat": 8, "mixing_scheme": "anderson", "longitudinal_bond_guard_freqs": "static",
            "longitudinal_bond_output_full": False})


class TestArchiveRoundTrip(unittest.TestCase):
    """A half-grid run through ``qlms.run`` (the TOML entry, so the key is
    accepted by the [mode.param] placement check) writes a FULL-grid bond
    archive that ``hwave_sc`` pairs to the same leading eigenvalue as the
    archive of the full-grid run."""

    def _flex(self, tmp, half):
        import hwave.qlms as qlms
        from tests.test_eliashberg_bond_entry import _IN1
        par = {"T": 0.5, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
               "Nmat": 32, "IterationMax": 200, "Mix": 0.3, "EPS": 8,
               "flex_hartree_fock": True, "longitudinal_bond_channels": True,
               "longitudinal_bond_output_full": True, "mixing_scheme": "linear",
               "longitudinal_bond_freq_batch": 4, "longitudinal_bond_half_grid": half}
        inp = {"mode": {"mode": "FLEX", "calc_scheme": "general", "param": par},
               "file": {"input": {"path_to_input": "",
                                  "interaction": {"path_to_input": _IN1, "Geometry": "geom.dat",
                                                  "Transfer": "transfer.dat",
                                                  "CoulombIntra": "coulombintra.dat",
                                                  "CoulombInter": "coulombinter.dat"}},
                        "output": {"path_to_output": tmp, "chiq_s": "chiq_s", "chiq_c": "chiq_c",
                                   "green": "green", "longitudinal_bond": "longitudinal_bond"}}}
        qlms.run(input_dict=inp)

    def test_hwave_sc_on_a_half_grid_archive(self):
        import hwave.sc as sc
        from tests.test_eliashberg_bond_entry import _sc_input
        lam, chi = {}, {}
        for half in (False, True):
            flex_dir = tempfile.mkdtemp()
            out = tempfile.mkdtemp()
            try:
                self._flex(flex_dir, half)
                with np.load(os.path.join(flex_dir, "longitudinal_bond.npz")) as z:
                    self.assertEqual(z["chi_s_w"].shape[0], 32)
                    self.assertEqual("longitudinal_bond_half_grid" in z.files, half)
                    chi[half] = np.array(z["chi_s_w"])
                lam[half] = sc.calc_eliashberg(_sc_input(flex_dir, out, 0.5, 32, (4, 4, 1)))
            finally:
                shutil.rmtree(flex_dir, ignore_errors=True)
                shutil.rmtree(out, ignore_errors=True)
        self.assertLessEqual(_rel(chi[False], chi[True]), _TOL)
        self.assertTrue(np.isfinite(lam[False]))
        self.assertLessEqual(abs(lam[True] - lam[False]), 1e-10 * max(1.0, abs(lam[False])))


# ---------------------------------------------------------------------------
# the mirrored slices' own conditioning (time reversal broken / preserved)
# ---------------------------------------------------------------------------

def _score(m):
    """The guard score of a block stack, from numpy.linalg.svd (independent
    of the production scorer): min(smin/smax, smin/max(1, smax))."""
    sv = np.linalg.svd(m, compute_uv=False)
    smax, smin = sv[..., 0], sv[..., -1]
    return np.minimum(smin / smax, smin / np.maximum(1.0, smax))


def _unitary(rng, n):
    q, _ = np.linalg.qr(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)))
    return q


_NMAT, _NVOL, _ND = 8, 4, 4          # spatial (4, 1, 1); norb = 2, one channel
_LP, _QP = 6, 1                       # the planted STORED block (nu > 0)


def _planted_map(trs, delta=1e-3, seed=330):
    """A synthetic full-grid bubble with X(q, l) = X(q, nmat - l)^dagger and a
    real symmetric spin vertex W (V(-q) = V(q)^T). One stored block (q = 1,
    l = 6) has the spin denominator 1 - X W = E with score ``delta``; its
    mirror (q = 1, l = 2) is 1 - X^dagger W, whose score is ~500x smaller
    (W has condition number ~30). ``trs``: also impose X(-q) = X(q)^T (time
    reversal), which puts the mirror's score at the stored block (q = 3,
    l = 6); otherwise every other block is random and small."""
    from hwave.solver import half_grid as hg
    rng = np.random.default_rng(seed)
    n = _ND
    O, _ = np.linalg.qr(rng.normal(size=(n, n)))
    W = O @ np.diag([1.0, 0.6, 0.3, 0.03]) @ O.T                 # real symmetric
    E = _unitary(rng, n) @ np.diag([1.0, 0.9, 0.8, delta]) @ _unitary(rng, n)
    M = -(E - np.eye(n)) @ np.linalg.inv(W)                      # spin: 1 - M W = E
    x = 0.02 * (rng.normal(size=(_NMAT, _NVOL, n, n)) + 1j * rng.normal(size=(_NMAT, _NVOL, n, n)))
    mq = hg.minus_q_index((_NVOL, 1, 1))
    stored = hg.grid_index(_NMAT)
    for l in stored:
        if trs:
            x[l] = 0.5 * (x[l] + x[l][mq].swapaxes(-1, -2))
        if l in (0, _NMAT // 2):
            x[l] = 0.5 * (x[l] + np.conj(x[l]).swapaxes(-1, -2))
            if trs:
                x[l] = 0.5 * (x[l] + x[l][mq].swapaxes(-1, -2))
    x[_LP, _QP] = M
    if trs:
        x[_LP, mq[_QP]] = M.T
    hg.fill_mirror(x, _NMAT)
    S = np.broadcast_to(W, (_NVOL, n, n)).astype(complex).copy()
    C = np.broadcast_to(0.1 * np.eye(n), (_NVOL, n, n)).astype(complex).copy()
    return x, S, C


def _map_run(x, S, C, half, cond_tol, nb=3, **kw):
    """One dress_and_build_w on a full-grid or half-grid store; returns
    ``("ok", res, warnings)`` or ``("refused", exc)``."""
    import types
    from hwave.solver import bond_channels as bc
    from hwave.solver import half_grid as hg
    from tests.test_flex_bond_dressing import _onsite_factors
    factors = _onsite_factors()
    view = types.SimpleNamespace(n_channels=_ND // 4)
    with flex_bond.BondBlockStore(_NMAT, _NVOL, _ND, 4, ("chibar",), {"W": "chibar"},
                                  half_grid=half) as st, \
            flex_bond.BondDeviceContext.for_view(np, S, C, view, 2) as dev:
        st.put_freq_batch("chibar", 0, st.nslices,
                          hg.half_from_full(x, _NMAT) if half else x)
        try:
            with mock.patch.object(bc.logger, "warning") as warn:
                res = flex_bond.dress_and_build_w(
                    st, dev, nb=nb, output_full=False, nmat=_NMAT, nvol=_NVOL, nd=4,
                    spatial_shape=(_NVOL, 1, 1), factors=factors, cond_tol=cond_tol, **kw)
        except (bc.BondConditioningError, ValueError) as exc:
            return ("refused", exc)
        recs = sorted((c.args[1], c.args[3], c.args[4]) for c in warn.call_args_list
                      if "grid index" in c.args[0])
        return ("ok", res, recs)


class TestMirrorGuardBrokenTimeReversal(unittest.TestCase):
    """The worst block sits at nu < 0 (q = 1, l = 2) and its stored partner
    passes: the half grid must refuse / warn exactly as the full grid."""

    TOL = 3e-5

    def setUp(self):
        self.x, self.S, self.C = _planted_map(trs=False)
        x = self.x
        mat = np.eye(_ND) - x @ self.S[None]                         # spin, every (l, q)
        sc = _score(mat)
        self.scores = sc
        # preconditions of the fixture: exactly one block below TOL, at nu < 0
        below = np.argwhere(sc < self.TOL)
        self.assertEqual(below.tolist(), [[2, _QP]])
        self.assertGreater(sc[_LP, _QP], 10 * self.TOL)
        charge = _score(np.eye(_ND) + x @ self.C[None])
        self.assertGreater(charge.min(), 0.1)

    def test_refuse(self):
        full = _map_run(self.x, self.S, self.C, False, self.TOL)
        half = _map_run(self.x, self.S, self.C, True, self.TOL)
        self.assertEqual(full[0], "refused")
        self.assertEqual(half[0], "refused")
        ef, eh = full[1], half[1]
        self.assertEqual((eh.channel, eh.l, eh.iq), (ef.channel, ef.l, ef.iq))
        self.assertEqual((eh.l, eh.iq), (2, _QP))
        self.assertAlmostEqual(eh.worst, ef.worst, delta=1e-12)
        self.assertAlmostEqual(eh.worst, float(self.scores[2, _QP]), delta=1e-12)
        self.assertIn("bosonic Matsubara index {}".format(2 * 2 - _NMAT), str(eh))

    def test_warn(self):
        full = _map_run(self.x, self.S, self.C, False, self.TOL, guard_policy="warn")
        half = _map_run(self.x, self.S, self.C, True, self.TOL, guard_policy="warn")
        self.assertEqual(full[0], "ok")
        self.assertEqual(half[0], "ok")
        self.assertEqual(full[2], [("spin", 2, 1)])
        self.assertEqual(half[2], full[2])
        self.assertEqual(half[1].guard_violations, full[1].guard_violations)
        self.assertTrue(half[1].mirror_guarded)
        # cond_min covers the whole grid: the mirror's minimum, as on the full grid
        self.assertAlmostEqual(half[1].cond_min_s, full[1].cond_min_s, delta=1e-12)
        self.assertAlmostEqual(full[1].cond_min_s, float(self.scores.min()), delta=1e-12)
        # the guard statistics count the stored slices only
        self.assertEqual(half[1].guard_blocks * _NMAT, full[1].guard_blocks * 5)

    def test_no_violation_matches(self):
        full = _map_run(self.x, self.S, self.C, False, 1e-7)
        half = _map_run(self.x, self.S, self.C, True, 1e-7)
        self.assertEqual((full[0], half[0]), ("ok", "ok"))
        self.assertAlmostEqual(half[1].cond_min_s, full[1].cond_min_s, delta=1e-12)
        self.assertAlmostEqual(half[1].cond_min_c, full[1].cond_min_c, delta=1e-12)
        # the minimum is NOT at the static slice
        l, q = np.unravel_index(np.argmin(self.scores), self.scores.shape)
        self.assertNotEqual(l, _NMAT // 2)

    def test_static_residual(self):
        """Under guard_freqs = "static" the nu < 0 slices are checked by their
        own solve residual: the residual is scaled (only for this test) so
        that only the largest one, the planted mirror block's, exceeds the
        1e-6 threshold; full and half grid must refuse (or warn) there."""
        from hwave.solver import bond_channels as bc
        real = bc.solve_residual
        rec = []

        def record(mat, chi, cb):
            r = real(mat, chi, cb)
            rec.append(np.array(r))
            return r

        with mock.patch.object(bc, "solve_residual", record):
            out = _map_run(self.x, self.S, self.C, False, self.TOL, nb=_NMAT,
                           guard_freqs="static")
        self.assertEqual(out[0], "ok")
        r = rec[0].reshape(_NMAT, _NVOL)                 # spin, the one full-grid batch
        order = np.argsort(r.ravel())[::-1]
        self.assertEqual(np.unravel_index(order[0], r.shape), (2, _QP))
        r1, r2 = r.ravel()[order[0]], r.ravel()[order[1]]
        self.assertGreater(r1, 4 * r2)
        k = 1e-6 / np.sqrt(r1 * r2)

        def scaled(mat, chi, cb):
            return real(mat, chi, cb) * k

        for policy in ("refuse", "warn"):
            with self.subTest(policy=policy):
                with mock.patch.object(bc, "solve_residual", scaled):
                    full = _map_run(self.x, self.S, self.C, False, self.TOL,
                                    guard_freqs="static", guard_policy=policy)
                    half = _map_run(self.x, self.S, self.C, True, self.TOL,
                                    guard_freqs="static", guard_policy=policy)
                if policy == "refuse":
                    self.assertEqual((full[0], half[0]), ("refused", "refused"))
                    for e in (full[1], half[1]):
                        self.assertIn("grid index l=2) and q-point index (1, 0, 0)", str(e))
                else:
                    self.assertEqual((full[0], half[0]), ("ok", "ok"))
                    self.assertEqual(half[1].guard_residual_violations, 1)
                    self.assertEqual(full[1].guard_residual_violations, 1)

    def test_every_violation_once_per_frequency(self):
        """A floor above most scores, one frequency per batch: every slice
        reports its worst block once, the self-mirrored l = 0 and nmat/2
        included -- the same records and count as the full grid (they are
        not guarded a second time as their own mirrors)."""
        full = _map_run(self.x, self.S, self.C, False, 0.99, nb=1, guard_policy="warn")
        half = _map_run(self.x, self.S, self.C, True, 0.99, nb=1, guard_policy="warn")
        self.assertEqual((full[0], half[0]), ("ok", "ok"))
        ls = sorted({r[1] for r in full[2]})
        self.assertEqual(ls, list(range(_NMAT)))
        self.assertEqual(half[2], full[2])
        self.assertEqual(half[1].guard_violations, full[1].guard_violations)

    def test_mutation_without_mirror_guard_passes_silently(self):
        """Mutation: the mirror guard disabled (the relation reported as
        holding). The half grid then passes where the full grid refuses, so
        test_refuse above would fail."""
        from hwave.solver import half_grid as hg
        calls = {"n": 0}

        def holds(*a, **k):
            calls["n"] += 1
            return {"holds": True, "deviation": 0.0, "chibar_max": 1.0, "vertex_max": 1.0}

        with mock.patch.object(hg, "time_reversal_check", holds):
            half = _map_run(self.x, self.S, self.C, True, self.TOL)
        self.assertGreater(calls["n"], 0)
        self.assertEqual(half[0], "ok")


class TestMirrorSkipWithTimeReversal(unittest.TestCase):
    """With time reversal every mirrored score is a stored one, (q, l') <->
    (-q, l): the skip is taken and the refusal / minimum is the same score."""

    def setUp(self):
        self.x, self.S, self.C = _planted_map(trs=True)

    def test_mirror_scores_equal_stored_scores(self):
        from hwave.solver import half_grid as hg
        x, S, C = self.x, self.S, self.C
        mq = hg.minus_q_index((_NVOL, 1, 1))
        holds, dev = hg.mirror_scores_are_stored(hg.half_from_full(x, _NMAT), S, C,
                                                 (_NVOL, 1, 1))
        self.assertTrue(holds)
        self.assertLessEqual(dev, 1e-15)
        for V, sign in ((S, -1.0), (C, 1.0)):
            sc = _score(np.eye(_ND) + sign * (x @ V[None]))
            for l in range(1, _NMAT // 2):
                np.testing.assert_allclose(sc[l], sc[_NMAT - l][mq], rtol=0, atol=1e-12)
        # and the broken fixture really is broken
        xb, Sb, Cb = _planted_map(trs=False)
        self.assertFalse(hg.mirror_scores_are_stored(hg.half_from_full(xb, _NMAT), Sb, Cb,
                                                     (_NVOL, 1, 1))[0])

    def test_a_vertex_breaking_time_reversal_is_detected(self):
        """chibar alone time-reversal symmetric is not enough: a complex
        Hermitian vertex (V(q)^T != V(-q)) must switch the mirror guard on."""
        from hwave.solver import half_grid as hg
        x, S, C = self.x, self.S.copy(), self.C
        S[:, 0, 1] += 0.05j
        S[:, 1, 0] -= 0.05j
        holds, dev = hg.mirror_scores_are_stored(hg.half_from_full(x, _NMAT), S, C,
                                                 (_NVOL, 1, 1))
        self.assertFalse(holds)
        self.assertGreater(dev, 1e-3)
        out = _map_run(x, S, C, True, 1e-7)
        self.assertEqual(out[0], "ok")
        self.assertTrue(out[1].mirror_guarded)

    def test_skip_taken_and_same_refusal_score(self):
        sc = _score(np.eye(_ND) - self.x @ self.S[None])
        tol = 3e-5
        self.assertEqual(sorted(map(tuple, np.argwhere(sc < tol).tolist())),
                         [(2, _QP), (_LP, 3)])
        full = _map_run(self.x, self.S, self.C, False, tol)
        half = _map_run(self.x, self.S, self.C, True, tol)
        self.assertEqual((full[0], half[0]), ("refused", "refused"))
        # the full grid meets the nu < 0 block first; the half grid finds its
        # time-reversal partner (-q, nmat - l) among the stored slices and
        # reports the refusal at the TRUE location, with the same score
        self.assertEqual((full[1].l, full[1].iq), (2, _QP))
        self.assertEqual((half[1].l, half[1].iq), (2, _QP))
        self.assertEqual(half[1].q, full[1].q)
        self.assertAlmostEqual(half[1].worst, full[1].worst, delta=1e-12)
        self.assertIn("grid index l=2) and q-point index (1, 0, 0)", str(half[1]))
        ok_f = _map_run(self.x, self.S, self.C, False, 1e-7)
        ok_h = _map_run(self.x, self.S, self.C, True, 1e-7)
        self.assertFalse(ok_h[1].mirror_guarded)
        self.assertAlmostEqual(ok_h[1].cond_min_s, ok_f[1].cond_min_s, delta=1e-12)
        self.assertAlmostEqual(ok_h[1].cond_min_c, ok_f[1].cond_min_c, delta=1e-12)

    def test_skip_is_logged_per_solver_decision(self):
        # solver level: real hopping skips, complex hopping guards the mirror
        for case, guarded in (("UV", False), ("complex_t", True)):
            with self.subTest(case=case):
                with self.assertLogs("hwave.solver.flex", level="INFO") as cm:
                    s, _gi, _f = _run(case, {"IterationMax": 2}, half=True, save=False)
                self.assertIs(s._bond_last.mirror_guarded, guarded)
                lines = [m for m in cm.output if "longitudinal_bond_half_grid (map" in m]
                self.assertEqual(len(lines), 1, lines)
                self.assertIn("guarded on their own" if guarded else "not guarded separately",
                              lines[0])


def _count_mirror_guards():
    """A patch counting the dress_batch calls of the mirror guard
    (guard_only=True) and the counter."""
    from hwave.solver import bond_channels as bc
    real = bc.dress_batch
    calls = {"n": 0}

    def spy(*a, **k):
        if k.get("guard_only"):
            calls["n"] += 1
        return real(*a, **k)

    return mock.patch.object(bc, "dress_batch", spy), calls


class TestSkipOnlyUnderRefuseAll(unittest.TestCase):
    """With time reversal the mirror is skipped only under guard_freqs =
    "all" with the "refuse" policy; "warn" and "static" always guard it, so
    their records and counts are the full grid's."""

    def setUp(self):
        self.x, self.S, self.C = _planted_map(trs=True)

    def _half(self, tol, **kw):
        p, calls = _count_mirror_guards()
        with p:
            out = _map_run(self.x, self.S, self.C, True, tol, **kw)
        return out, calls["n"]

    def test_refuse_all_skips(self):
        out, n = self._half(1e-7)
        self.assertEqual(out[0], "ok")
        self.assertFalse(out[1].mirror_guarded)
        self.assertEqual(n, 0)

    def test_warn_guards_the_mirror(self):
        for tol in (3e-5, 0.99):
            with self.subTest(tol=tol):
                full = _map_run(self.x, self.S, self.C, False, tol, nb=1, guard_policy="warn")
                half, n = self._half(tol, nb=1, guard_policy="warn")
                self.assertEqual((full[0], half[0]), ("ok", "ok"))
                self.assertTrue(half[1].mirror_guarded)
                self.assertGreater(n, 0)
                self.assertEqual(half[1].guard_violations, full[1].guard_violations)
                self.assertEqual(half[2], full[2])
                self.assertGreater(len(full[2]), 1)

    def test_static_guards_the_mirror(self):
        from hwave.solver import bond_channels as bc
        real = bc.solve_residual

        def huge(mat, chi, cb):
            return real(mat, chi, cb) * 0 + 1.0       # every block "violates"

        with mock.patch.object(bc, "solve_residual", huge):
            full = _map_run(self.x, self.S, self.C, False, 1e-7, nb=1, guard_freqs="static",
                            guard_policy="warn")
            half, n = self._half(1e-7, nb=1, guard_freqs="static", guard_policy="warn")
            ref, n_ref = self._half(1e-7, nb=1, guard_freqs="static")
        self.assertTrue(half[1].mirror_guarded)
        self.assertGreater(n, 0)
        self.assertEqual(half[1].guard_residual_violations, full[1].guard_residual_violations)
        self.assertEqual(full[1].guard_residual_violations, 2 * _NMAT)
        self.assertEqual(half[2], full[2])
        # "static" under "refuse": refused like the full grid (at the first
        # batch, the static slot, which is its own mirror)
        self.assertEqual(ref[0], "refused")
        self.assertEqual(n_ref, 0)

    def test_static_refusal_only_at_negative_frequencies(self):
        """Under "static" with "refuse" the residual exceeds its threshold
        only for the slices 0 < l < nmat/2 (nu < 0): the full grid refuses,
        and so must the half grid, which can see those slices only through
        its own mirror solves (a skip under "static" would pass)."""
        from hwave.solver import bond_channels as bc
        real_db, real_res = bc.dress_batch, bc.solve_residual

        def db(*a, **k):
            ls = k.get("ls")
            l0 = k["l0"]
            if ls is not None:
                ls = np.asarray(ls)
                neg = bool(np.all((ls > 0) & (ls < _NMAT // 2)))
            else:
                neg = 0 < l0 < _NMAT // 2
            if neg:
                with mock.patch.object(bc, "solve_residual",
                                       lambda m, c, b: real_res(m, c, b) * 0 + 1.0):
                    return real_db(*a, **k)
            return real_db(*a, **k)

        with mock.patch.object(bc, "dress_batch", db):
            full = _map_run(self.x, self.S, self.C, False, 1e-7, nb=1, guard_freqs="static")
            half = _map_run(self.x, self.S, self.C, True, 1e-7, nb=1, guard_freqs="static")
        self.assertEqual(full[0], "refused")
        self.assertEqual(half[0], "refused")
        for out in (full, half):
            l = int(str(out[1]).split("grid index l=")[1].split(")")[0])
            self.assertTrue(0 < l < _NMAT // 2, str(out[1]))

    def test_endpoint_worst_with_larger_batches(self):
        """The batch worst at the static slot (nb = 3; review experiment):
        under "warn" the mirror is guarded, so the records are those of a run
        with the mirror guard forced; under "refuse" the decision is the
        full grid's (the first location reported can differ)."""
        from hwave.solver import half_grid as hg
        x = self.x.copy()
        W = self.S[0].real
        x[_NMAT // 2, 0] = (1 - 1e-6) * np.linalg.inv(W)
        forced = {"holds": False, "deviation": 1.0, "chibar_max": 1.0, "vertex_max": 1.0,
                  "hermitian_deviation": 0.0}
        for policy in ("warn", "refuse"):
            with self.subTest(policy=policy):
                full = _map_run(x, self.S, self.C, False, 3e-5, nb=3, guard_policy=policy)
                half = _map_run(x, self.S, self.C, True, 3e-5, nb=3, guard_policy=policy)
                with mock.patch.object(hg, "time_reversal_check", lambda *a, **k: forced):
                    guarded = _map_run(x, self.S, self.C, True, 3e-5, nb=3,
                                       guard_policy=policy)
                self.assertEqual(half[0], full[0])
                self.assertEqual(guarded[0], full[0])
                if policy == "warn":
                    self.assertEqual(half[2], guarded[2])
                    self.assertIn(("spin", 2, _QP), half[2])
                    self.assertEqual(half[1].guard_violations, guarded[1].guard_violations)


class TestViolationRecordSchema(unittest.TestCase):
    """A direct dress_batch caller sees the documented record keys only."""

    def test_keys(self):
        from hwave.solver import bond_channels as bc
        keys = {"channel", "l", "q", "kind", "value", "iteration"}
        x, S, _C = _planted_map(trs=False)
        v = []
        with mock.patch.object(bc.logger, "warning"):
            bc.dress_batch(x[:3], S, "spin", l0=0, nmat=_NMAT, spatial_shape=(_NVOL, 1, 1),
                           cond_tol=0.99, guard_policy="warn", violations=v)
            real = bc.solve_residual
            with mock.patch.object(bc, "solve_residual",
                                   lambda m, c, b: real(m, c, b) * 0 + 1.0):
                bc.dress_batch(x[:3], S, "spin", l0=0, nmat=_NMAT,
                               spatial_shape=(_NVOL, 1, 1), guard_freqs="static",
                               guard_policy="warn", violations=v)
        self.assertEqual({r["kind"] for r in v}, {"conditioning", "residual"})
        for r in v:
            self.assertEqual(set(r), keys)


class TestNearMissTimeReversal(unittest.TestCase):
    """Relation (3) holds only to ~5e-13: the mirrored score then differs from
    its stored partner's by a relative ~1.6e-8 (the block is nearly singular).
    A floor placed between the two must give the full grid's decision."""

    def _maps(self, q_pert=3, sign=1.0):
        from hwave.solver import half_grid as hg
        x, S, C = _planted_map(trs=True)
        rng = np.random.default_rng(7)
        z = rng.normal(size=(_ND, _ND)) + 1j * rng.normal(size=(_ND, _ND))
        x = x.copy()
        x[_LP, q_pert] += sign * 5e-13 * np.abs(x).max() / np.abs(z).max() * z
        hg.fill_mirror(x, _NMAT)
        check = hg.time_reversal_check(hg.half_from_full(x, _NMAT), S, C, (_NVOL, 1, 1))
        self.assertTrue(check["holds"])
        self.assertGreater(check["deviation"], 1e-13)
        sc = _score(np.eye(_ND) - x @ S[None])
        stored, mirror = float(sc[_LP, 3]), float(sc[2, _QP])
        # the two scores differ (by ~1.6e-8 relative here), far more than
        # the 1e-12 of the relation itself
        self.assertGreater(abs(stored - mirror) / mirror, 1e-10)
        return x, S, C, stored, mirror

    def test_same_decision_both_orders(self):
        orders = set()
        for sign in (1.0, -1.0):
            x, S, C, stored, mirror = self._maps(sign=sign)
            orders.add(stored < mirror)
            tol = 0.5 * (stored + mirror)
            with self.subTest(stored_below=stored < mirror):
                full = _map_run(x, S, C, False, tol)
                half = _map_run(x, S, C, True, tol)
                self.assertEqual(full[0], "refused")
                self.assertEqual(half[0], full[0])
                self.assertEqual((half[1].l, half[1].iq), (full[1].l, full[1].iq))
                self.assertAlmostEqual(half[1].worst, full[1].worst, delta=1e-15)
                f = _map_run(x, S, C, False, tol, nb=1, guard_policy="warn")
                h = _map_run(x, S, C, True, tol, nb=1, guard_policy="warn")
                self.assertEqual(h[2], f[2])
                self.assertEqual(h[1].guard_violations, f[1].guard_violations)
        # both orders were exercised: the floor between the two scores once
        # catches only the stored block, once only the mirrored one
        self.assertEqual(orders, {True, False})

    def test_mutation_without_margin(self):
        """The margin removed: the skip then decides from the stored score
        alone and the half grid's decision differs from the full grid's."""
        from hwave.solver import half_grid as hg
        x, S, C, stored, mirror = self._maps()
        tol = 0.5 * (stored + mirror)
        full = _map_run(x, S, C, False, tol)
        with mock.patch.object(hg, "decision_margin", lambda *a, **k: 0.0):
            half = _map_run(x, S, C, True, tol)
        same = half[0] == full[0] and (half[0] != "refused" or (
            (half[1].l, half[1].iq, half[1].worst) == (full[1].l, full[1].iq, full[1].worst)))
        self.assertFalse(same)

    def test_margin_bound(self):
        from hwave.solver import half_grid as hg
        chk = {"deviation": 5e-13, "chibar_max": 20.0, "vertex_max": 1.0,
               "hermitian_deviation": 0.0}
        m = hg.decision_margin(chk, 4, 1e-3)
        self.assertAlmostEqual(m, 4 * (2 * 16 * 5e-13 * 20.0) * (1 + 1e3), delta=1e-18)
        self.assertEqual(hg.decision_margin(dict(chk, deviation=0.0), 4, 1e-3),
                         hg.DECISION_MARGIN_FLOOR)
        # the Hermiticity deviation of the vertices enters the bound too
        self.assertAlmostEqual(hg.decision_margin(dict(chk, deviation=0.0,
                                                       hermitian_deviation=5e-13), 4, 1e-3), m,
                               delta=1e-18)


class TestMirrorGuardDecide(unittest.TestCase):
    """guard_exact = "decide" with the interval guard on the mirrored rows:
    a batch the interval passes returns a bracket from guard_only without a
    solve, and the map's cond_min bracket holds the exact whole-grid minimum."""

    def test_bracket_and_early_return(self):
        from hwave.solver import bond_channels as bc
        x, S, C = _planted_map(trs=False, delta=0.5)
        sc = np.minimum(_score(np.eye(_ND) - x @ S[None]).min(),
                        _score(np.eye(_ND) + x @ C[None]).min())
        solves = {"n": 0}
        real = bc._threaded.solve

        def count(*a, **k):
            solves["n"] += 1
            return real(*a, **k)

        with mock.patch.object(bc._threaded, "solve", count):
            out = _map_run(x, S, C, True, 1e-7, guard_method="interval", guard_exact="decide")
        self.assertEqual(out[0], "ok")
        res = out[1]
        self.assertTrue(res.mirror_guarded)
        # only the stored batches were solved (2 channels x 2 batches)
        self.assertEqual(solves["n"], 4)
        brackets = [b for b in (res.cond_bracket_s, res.cond_bracket_c) if b is not None]
        self.assertTrue(brackets)
        for low, up, _i in brackets:
            self.assertLessEqual(low, up)
        self.assertLessEqual(min(b[0] for b in brackets), sc + 1e-12)


class TestSeedSymmetry(unittest.TestCase):
    """Under the half grid a sigma_init seed must satisfy
    Sigma(k, -i w) = Sigma(k, i w)^dagger."""

    def test_perturbed_seed_is_refused(self):
        from hwave.solver.flex import FLEX
        with tempfile.TemporaryDirectory() as tmp:
            info, _ = _solver("UV", {}, True)
            with _input("UV") as r:
                s = FLEX(r.get_param("ham"), {}, info)
                gi = r.get_param("green")
                s.solve(gi, tmp)
                s.save_results({"path_to_output": tmp, "sigma": "sigma.npz"}, gi)
                path = os.path.join(tmp, "sigma.npz")
                with np.load(path) as z:
                    data = {k: np.array(z[k]) for k in z.files}
                bad = dict(data)
                pert = np.zeros_like(data["sigma_fluct"])
                pert[:, 0] += 1e-6 * np.eye(pert.shape[-1])
                bad["sigma_fluct"] = data["sigma_fluct"] + pert
                bad["sigma"] = data["sigma"] + pert
                bad_path = os.path.join(tmp, "bad.npz")
                np.savez(bad_path, **bad)
                expected = (1, s.nmat, s.lattice.nvol, s.norb, s.norb)
                good_env, bad_env = s._read_sigma(path), s._read_sigma(bad_path)
                s._assemble_static_seed({"sigma_init_envelope": good_env}, expected)
                with self.assertRaisesRegex(ValueError, "sigma_fluct\\(k, -i w\\)"):
                    s._assemble_static_seed({"sigma_init_envelope": bad_env}, expected)
                # the full grid does not need the symmetry
                info_f, _ = _solver("UV", {}, False)
                sf = FLEX(r.get_param("ham"), {}, info_f)
                sf._assemble_static_seed({"sigma_init_envelope": sf._read_sigma(bad_path)},
                                         expected)


# ---------------------------------------------------------------------------
# mutations: each breaks the half-grid bookkeeping one way; the oracle fails
# ---------------------------------------------------------------------------

def _mut_expand(kind):
    """A mutated :func:`half_grid.expand` (and fill_mirror) of one kind."""
    calls = {"n": 0}

    def expand(ab, nmat, ba=None, out=None, chunk=None):
        calls["n"] += 1
        h0 = int(nmat) // 2
        ba = ab if ba is None else ba
        if out is None:
            out = np.empty((nmat,) + tuple(ab.shape[1:]), dtype=ab.dtype)
        out[h0:] = ab[:h0]
        out[0] = 0 if kind == "nyquist" else ab[h0]
        for l in range(1, h0):
            h = h0 - l
            if kind == "partner":
                h = min(h + 1, h0)                  # nmat - l - 1 instead of nmat - l
            x = np.conj(ba[h])
            out[l] = x if kind == "no_transpose" else x.swapaxes(-1, -2)
        return out

    def fill_mirror(full, nmat):
        calls["n"] += 1
        h0 = int(nmat) // 2
        if kind == "nyquist":
            full[0] = 0
        for l in range(1, h0):
            src = nmat - l - (1 if kind == "partner" else 0)
            x = np.conj(full[src])
            full[l] = x if kind == "no_transpose" else x.swapaxes(-1, -2)
        return full

    return expand, fill_mirror, calls


class TestMutations(unittest.TestCase):
    """Each mutation must make the full-vs-half oracle fail (after asserting
    that the mutated code really ran)."""

    def _dev(self, case="UV", extra=None):
        # the density guard of issue #199 refuses most mutated runs outright
        # (a broken mirror makes the equal-time density non-Hermitian); it
        # is relaxed here so that the deviation itself is what fails
        extra = dict({"longitudinal_bond_output_full": True, "flex_hf_density_tol": 1e3},
                     **(extra or {}))
        full = _run(case, extra, half=False)
        return full, extra

    def _mutated_worst(self, patches, calls, case="UV", extra=None):
        full, extra = self._dev(case, extra)
        with contextlib.ExitStack() as st:
            for p in patches:
                st.enter_context(p)
            half = _run(case, extra, half=True)
        self.assertGreater(calls["n"], 0, "the mutation did not run")
        dev = _compare(full, half)
        return max(dev.values()), dev

    def _expand_mutation(self, kind):
        from hwave.solver import half_grid as hg
        expand, fill, calls = _mut_expand(kind)
        worst, dev = self._mutated_worst(
            [mock.patch.object(hg, "expand", expand), mock.patch.object(hg, "fill_mirror", fill)],
            calls)
        self.assertGreater(worst, 1e-6, (kind, dev))
        # the self-energy itself moves, not only the published members
        self.assertGreater(dev["green_info:sigma"], 1e-8, kind)

    def test_control_unmutated_reimplementation_passes(self):
        # the mutation harness itself is faithful: kind None reproduces
        from hwave.solver import half_grid as hg
        expand, fill, calls = _mut_expand(None)
        worst, dev = self._mutated_worst(
            [mock.patch.object(hg, "expand", expand), mock.patch.object(hg, "fill_mirror", fill)],
            calls)
        self.assertLessEqual(worst, _TOL, dev)

    def test_drop_the_nyquist_slot(self):
        self._expand_mutation("nyquist")

    def test_conjugate_without_transpose(self):
        self._expand_mutation("no_transpose")

    def test_partner_off_by_one(self):
        self._expand_mutation("partner")

    def test_miss_the_pair_swap(self):
        # the mirror of pair (alpha, beta) read from (alpha, beta) itself
        from hwave.solver import half_grid as hg
        calls = {"n": 0}

        def get_pair_full(self, name, alpha, beta, xp=None):
            if not self.half_grid:
                return real(self, name, alpha, beta, xp=xp)
            calls["n"] += 1
            ab = self.get_pair(name, alpha, beta)
            full = hg.expand(ab, self.nmat, ba=ab)
            return full if xp is None else flex_bond._bk.to_device(full, xp)

        real = flex_bond.BondBlockStore.get_pair_full
        worst, dev = self._mutated_worst(
            [mock.patch.object(flex_bond.BondBlockStore, "get_pair_full", get_pair_full)], calls)
        self.assertGreater(dev["green_info:sigma"], 1e-6, dev)

    def test_store_without_the_nyquist_slot(self):
        # the store itself drops the Nyquist slot: the half-grid run cannot
        # rebuild l = 0 and must not pass silently
        from hwave.solver import half_grid as hg
        calls = {"n": 0}
        real_hff = hg.half_from_full

        def half_from_full(full, nmat):
            calls["n"] += 1
            out = real_hff(full, nmat)
            out[-1] = 0
            return out

        worst, dev = self._mutated_worst(
            [mock.patch.object(hg, "half_from_full", half_from_full)], calls)
        self.assertGreater(worst, 1e-6, dev)


if __name__ == "__main__":
    unittest.main()
