"""The bare-H0 hint temperature of ``hwave.initial`` (``hwave.initial.hints``,
design section 6.4, R1)."""
import dataclasses
import math
import os
import shutil
import tempfile
import time
import unittest

from tests import initial_fixtures

WIDE = ("wide Wannier window: T_hint uses the active-window span W_act = {}, "
        "not the all-band span W_all = {}")
NO_MESH = ("no k-mesh axis longer than 2: T_hint = W0/8 is a fixed cap, not a "
           "mesh-derived scale; give --t-hint")


def _chain(t, a=0):
    return {((1, 0, 0), a, a): -t, ((-1, 0, 0), a, a): -t}


def _spin_transfer_text(norb, entries):
    """A 9-column (spin-resolved) Transfer file of ``{(R, a, b, s, t): v}``
    (0-based), Hermitian-closed."""
    closed = {}
    for (r, a, b, s, t), v in entries.items():
        closed[(r, a, b, s, t)] = closed.get((r, a, b, s, t), 0) + v
        neg = tuple(-x for x in r)
        if (neg, b, a, t, s) != (r, a, b, s, t):
            closed[(neg, b, a, t, s)] = closed.get((neg, b, a, t, s), 0) + complex(v).conjugate()
    lines = ["Transfer in wannier90-like format for uhfk w90spin/1",
             str(norb), "1", " 1"]
    for (r, a, b, s, t), v in closed.items():
        v = complex(v)
        lines.append("{} {} {} {} {} {} {} {:.15g} {:.15g}".format(
            *r, a + 1, b + 1, s, t, v.real, v.imag))
    return "\n".join(lines) + "\n"


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_thint_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def load(self, **kw):
        from hwave.solver import initial_facade as fac
        return fac.load_problem(initial_fixtures.write_input(self.dest(), **kw))

    def hint(self, problem, **kw):
        from hwave.initial import hints
        return hints.t_hint_bare(problem, **kw)

    def assertRel(self, a, b, rtol=1e-12):
        self.assertLessEqual(abs(a - b), rtol * max(1.0, abs(b)), (a, b))


class TestTHint(_Tmp):
    def square(self, L, **kw):
        from hwave.solver import initial_facade as fac
        return fac.load_problem(initial_fixtures.square(self.dest(), L, **kw))

    def plaquette(self, T=0.0):
        """A 2 x 2 x 1 lattice (one bond per axis: on a length-2 axis R and
        -R are the same bond)."""
        return self.load(cellshape=(2, 2, 1),
                         transfer={((1, 0, 0), 0, 0): -1.0, ((0, 1, 0), 0, 0): -1.0},
                         param={"Ncond": 4, "T": T})

    def test_v1b_value(self):
        h = self.hint(self.square(16))
        self.assertEqual(h.W0, 8.0)
        self.assertEqual(h.W_act, 8.0)
        self.assertEqual(h.W_all, 8.0)
        self.assertEqual(h.L_min, 16)
        self.assertEqual(h.T_hint, 0.25)
        self.assertEqual(h.t_hint_source, "rule")
        self.assertEqual(h.t_hint_rule, "mesh")
        self.assertEqual(h.W_source, "active_window")
        self.assertEqual(h.rule, "T_hint = max(W0/(2*L_min), 1e-4*W0)")
        self.assertEqual(len(h.E_F), 1)
        self.assertLessEqual(abs(h.E_F[0]), 1e-12)
        self.assertEqual(h.D, (8.0,))
        self.assertEqual(h.W_act_groups, (8.0,))
        self.assertEqual(h.warnings, ())
        self.assertIn("1/(4*T_hint)", h.effect)
        self.assertIn("0.25", h.effect)
        with self.assertRaises(dataclasses.FrozenInstanceError):
            h.T_hint = 1.0

    def test_t_positive_input(self):
        h = self.hint(self.square(4, T=0.1))
        self.assertEqual(h.T_hint, 0.1)
        self.assertEqual(h.t_hint_source, "input")
        self.assertIsNone(h.t_hint_rule)
        # the inputs of the rule are still computed and recorded
        self.assertEqual(h.L_min, 4)
        self.assertEqual(h.W0, 8.0)
        # a T > 0 input with no mesh axis: no cap, so no cap warning
        h2 = self.hint(self.plaquette(T=0.1))
        self.assertEqual(h2.T_hint, 0.1)
        self.assertIsNone(h2.L_min)
        self.assertIsNone(h2.t_hint_rule)
        self.assertEqual(h2.warnings, ())

    def test_user_t_hint_precedence(self):
        for T in (0.0, 0.1):
            with self.subTest(T=T):
                p = self.square(4, T=T)
                h = self.hint(p, t_user=0.05)
                self.assertEqual(h.T_hint, 0.05)
                self.assertEqual(h.t_hint_source, "user")
                self.assertIsNone(h.t_hint_rule)
                # the problem (and the runs) keep the input T
                self.assertEqual(p.T, T)
        p = self.square(4)
        for bad in (0.0, -0.1, float("nan"), float("inf"), True, "0.1"):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError):
                    self.hint(p, t_user=bad)

    def test_origin_invariance(self):
        C = 7.3
        # a two-orbital chain with hybridisation off half filling, so E_F
        # is not at a symmetric point
        base = dict(_chain(1.0, 0)); base.update(_chain(0.5, 1))
        base[((0, 0, 0), 0, 1)] = -0.3
        base[((0, 0, 0), 1, 0)] = -0.3
        base[((0, 0, 0), 1, 1)] = 0.7
        shifted = dict(base)
        shifted[((0, 0, 0), 0, 0)] = C
        shifted[((0, 0, 0), 1, 1)] = 0.7 + C
        kw = dict(cellshape=(12, 1, 1), norb=2, param={"Ncond": 10, "T": 0.0})
        h0 = self.hint(self.load(transfer=base, **kw))
        h1 = self.hint(self.load(transfer=shifted, **kw))
        for name in ("W0", "W_act", "W_all", "T_hint"):
            self.assertRel(getattr(h1, name), getattr(h0, name))
        for a, b in zip(h1.D, h0.D):
            self.assertRel(a, b)
        for a, b in zip(h1.W_act_groups, h0.W_act_groups):
            self.assertRel(a, b)
        for a, b in zip(h1.E_F, h0.E_F):
            self.assertRel(a, b + C)
        self.assertEqual(h1.L_min, h0.L_min)

    def two_band(self, ligand):
        """A flat ligand level and a half-filled chain of width 4."""
        transfer = dict(_chain(1.0, 0))
        transfer[((0, 0, 0), 1, 1)] = ligand
        return self.load(cellshape=(16, 1, 1), norb=2, transfer=transfer,
                         param={"Ncond": 48, "T": 0.0})

    def test_active_window(self):
        h = self.hint(self.two_band(-22.0))
        self.assertAlmostEqual(h.W_act, 4.0, delta=1e-9)
        self.assertAlmostEqual(h.W_all, 24.0, delta=1e-9)
        self.assertAlmostEqual(h.W0, 4.0, delta=1e-9)
        self.assertAlmostEqual(h.T_hint, 4.0 / 32, delta=1e-12)
        self.assertAlmostEqual(h.D[0], 4.0, delta=1e-9)
        self.assertEqual(h.warnings, (WIDE.format(4, 24),))
        h2 = self.hint(self.two_band(-40.0))
        self.assertEqual(h2.T_hint, h.T_hint)
        self.assertAlmostEqual(h2.W_all, 42.0, delta=1e-9)
        self.assertEqual(h2.warnings, (WIDE.format(4, 42),))

    def test_band_insulator_window(self):
        # lower band [-4, -2] full, upper band [2, 4] empty, a ligand at -30
        transfer = dict(_chain(0.5, 0)); transfer.update(_chain(0.5, 1))
        transfer[((0, 0, 0), 0, 0)] = -3.0
        transfer[((0, 0, 0), 1, 1)] = 3.0
        transfer[((0, 0, 0), 2, 2)] = -30.0
        h = self.hint(self.load(cellshape=(16, 1, 1), norb=3, transfer=transfer,
                                param={"Ncond": 64, "T": 0.0}))
        self.assertLessEqual(abs(h.E_F[0]), 1e-12)
        # the Fermi bands are the bands bordering the gap
        self.assertAlmostEqual(h.D[0], 8.0, delta=1e-9)
        self.assertAlmostEqual(h.W_act, 8.0, delta=1e-9)
        self.assertAlmostEqual(h.W_all, 34.0, delta=1e-9)
        self.assertAlmostEqual(h.T_hint, 8.0 / 32, delta=1e-12)
        self.assertEqual(h.warnings, (WIDE.format(8, 34),))

    def spin_chain(self, *, t_up, t_dn, offset, ncond, two_sz, subshape=None):
        entries = {((1, 0, 0), 0, 0, 0, 0): -t_up, ((1, 0, 0), 0, 0, 1, 1): -t_dn}
        if offset:
            entries[((0, 0, 0), 0, 0, 1, 1)] = offset
        return self.load(cellshape=(16, 1, 1), norb=1, transfer=_chain(1.0),
                         param=dict({"Ncond": ncond, "2Sz": two_sz, "T": 0.0},
                                    **({"SubShape": list(subshape)} if subshape else {})),
                         extra_files={"transfer.dat": _spin_transfer_text(1, entries)})

    def test_2sz_offset_and_empty_sector(self):
        h0 = self.hint(self.spin_chain(t_up=1.0, t_dn=1.0, offset=0.0, ncond=16, two_sz=0))
        h1 = self.hint(self.spin_chain(t_up=1.0, t_dn=1.0, offset=30.0, ncond=16, two_sz=0))
        self.assertEqual(len(h0.E_F), 2)
        self.assertRel(h1.T_hint, h0.T_hint)
        self.assertRel(h1.W_act, 4.0)
        self.assertRel(h1.E_F[1], h0.E_F[1] + 30.0)
        self.assertAlmostEqual(h1.W_all, 34.0, delta=1e-9)
        # an empty down sector (2Sz = Ncond): its wider band never enters
        h2 = self.hint(self.spin_chain(t_up=1.0, t_dn=3.0, offset=0.0, ncond=8, two_sz=8))
        self.assertEqual(len(h2.E_F), 2)
        self.assertIsNone(h2.E_F[1])
        self.assertIsNone(h2.D[1])
        self.assertIsNone(h2.W_act_groups[1])
        self.assertAlmostEqual(h2.W_act_groups[0], 4.0, delta=1e-9)
        self.assertAlmostEqual(h2.W_act, 4.0, delta=1e-9)
        self.assertAlmostEqual(h2.W_all, 12.0, delta=1e-9)

    def test_full_or_empty_everywhere_uses_w_all(self):
        # a full band: no window anywhere, W_act = W_all
        h = self.hint(self.square(4, ncond=32))
        self.assertEqual(h.E_F, (None,))
        self.assertEqual(h.W_act, h.W_all)
        self.assertEqual(h.W_act, 8.0)

    def test_floor(self):
        # a single flat level: W_act = 0, W0 is the floor
        h = self.hint(self.load(cellshape=(4, 1, 1), transfer={((0, 0, 0), 0, 0): 0.5},
                                param={"Ncond": 4, "T": 0.0}))
        self.assertEqual(h.W_act, 0.0)
        self.assertEqual(h.W0, 1e-3)
        self.assertEqual(h.W_source, "floor")
        self.assertEqual(h.T_hint, max(1e-3 / 8, 1e-7))

    def test_structural_axes(self):
        ladder = dict(_chain(1.0))
        # one rung bond on the length-2 axis (R and -R are the same bond there)
        ladder[((0, 1, 0), 0, 0)] = -1.0
        h = self.hint(self.load(cellshape=(64, 2, 1), transfer=ladder,
                                param={"Ncond": 128, "T": 0.0}))
        self.assertEqual(h.L_min, 64)
        self.assertEqual(h.t_hint_rule, "mesh")
        self.assertRel(h.T_hint, h.W0 / 128, 1e-15)
        bilayer = initial_fixtures.square_transfer()
        bilayer[((0, 0, 1), 0, 0)] = -1.0
        h = self.hint(self.load(cellshape=(32, 32, 2), transfer=bilayer,
                                param={"Ncond": 2048, "T": 0.0}))
        self.assertEqual(h.L_min, 32)
        self.assertRel(h.T_hint, h.W0 / 64, 1e-15)
        h = self.hint(self.plaquette())
        self.assertIsNone(h.L_min)
        self.assertEqual(h.T_hint, h.W0 / 8)
        self.assertEqual(h.t_hint_rule, "no_mesh_axis_cap")
        self.assertEqual(h.warnings, (NO_MESH,))
        self.assertEqual(h.rule, "T_hint = W0/8 (no k-mesh axis longer than 2)")

    def test_subshape_same_bands(self):
        # the folded cell has the same spectrum, so the same T_hint
        h1 = self.hint(self.square(8))
        h2 = self.hint(self.square(8, subshape=(2, 2, 1)))
        for name in ("W0", "W_act", "W_all", "T_hint"):
            self.assertRel(getattr(h2, name), getattr(h1, name))
        self.assertEqual(h2.L_min, 8)
        # the primitive-cell bands are the spectrum of the solver's T(k)
        import numpy as np
        from hwave.initial import hints
        for p in (self.square(8, subshape=(2, 2, 1)),
                  self.spin_chain(t_up=1.0, t_dn=3.0, offset=30.0, ncond=16, two_sz=0,
                                  subshape=(4, 1, 1))):
            mine = np.sort(np.concatenate([e.ravel() for e in hints.group_bands(p)]))
            ref = np.sort(np.linalg.eigvalsh(np.asarray(p.ham_trans_k)).ravel())
            np.testing.assert_allclose(mine, ref, rtol=0, atol=1e-12)
        # per-group bands of a fixed 2Sz on a folded cell
        kw = dict(t_up=1.0, t_dn=3.0, offset=30.0, ncond=16, two_sz=0)
        g1 = self.hint(self.spin_chain(**kw))
        g2 = self.hint(self.spin_chain(subshape=(4, 1, 1), **kw))
        self.assertEqual(len(g2.E_F), 2)
        for name in ("E_F", "D", "W_act_groups"):
            for a, b in zip(getattr(g2, name), getattr(g1, name)):
                self.assertRel(a, b)
        self.assertRel(g2.W_act, 12.0)
        self.assertRel(g2.T_hint, g1.T_hint)

    def test_spin_orbital_mode(self):
        # the square lattice in the interleaved layout 2*orb + spin
        transfer = {(r, s, s): v for (r, a, b), v in
                    initial_fixtures.square_transfer().items() for s in (0, 1)}
        p = self.load(cellshape=(8, 8, 1), norb=2, transfer=transfer,
                      param={"Ncond": 64, "T": 0.0},
                      mode={"enable_spin_orbital": True})
        self.assertTrue(p.enable_spin_orbital)
        h = self.hint(p)
        self.assertEqual(h.W0, 8.0)
        self.assertEqual(h.T_hint, 0.5)
        self.assertEqual(h.W_act_groups, (8.0,))
        p2 = self.load(cellshape=(8, 8, 1), norb=2, transfer=transfer,
                       param={"Ncond": 64, "T": 0.0, "SubShape": [2, 2, 1]},
                       mode={"enable_spin_orbital": True})
        self.assertRel(self.hint(p2).T_hint, 0.5)


if __name__ == "__main__":
    unittest.main()
