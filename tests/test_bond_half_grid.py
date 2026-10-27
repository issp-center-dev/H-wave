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
    called, and the two runs are bit-identical."""

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
        self.assertEqual(sorted(ga), sorted(go))
        for k in ga:
            if isinstance(ga[k], np.ndarray):
                np.testing.assert_array_equal(ga[k], go[k], err_msg=k)
        self.assertEqual(sorted(fa), sorted(fo))
        for fn in fa:
            self.assertEqual(sorted(fa[fn]), sorted(fo[fn]))
            self.assertNotIn("longitudinal_bond_half_grid", fa[fn])
            for m in fa[fn]:
                np.testing.assert_array_equal(fa[fn][m], fo[fn][m], err_msg=fn + ":" + m)
        self.assertEqual(float(sa.mu), float(so.mu))


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
