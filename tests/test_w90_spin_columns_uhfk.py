#!/usr/bin/env python3
"""UHFk consuming the spin-resolved (9-column) Transfer format, issue #323
phase 2: the acceptance cases of the request.

1. spin-diagonal, spin-independent: the 9-column file gives the same H(k)
   and the same converged result as the 7-column file;
2. the same model accepts a fixed 2Sz;
3. a spin-diagonal but spin-asymmetric (Zeeman-like) file accepts 2Sz, a
   spin-flip (SOC) file refuses it;
4. an SOC file matches the equivalent enable_spin_orbital (interleaved)
   input after the index reorder;
5. a sublattice fold keeps the down-spin block of a spin-resolved table.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver import spin_orbital

MARKED = "Transfer in wannier90-like format for uhfk w90spin/1"


def _write_geom(d, norb):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(norb))
        for i in range(norb):
            f.write("{} 0 0\n".format(0.1 * i))


def _closed(entries):
    """Hermitian closure of {(R, a, b, s, t): v} (1-based a, b)."""
    out = {}
    for (r, a, b, s, t), v in entries.items():
        out[(r, a, b, s, t)] = out.get((r, a, b, s, t), 0) + v
        neg = tuple(-x for x in r)
        if (neg, b, a, t, s) == (r, a, b, s, t):
            continue                      # self-conjugate on-site term
        out[(neg, b, a, t, s)] = out.get((neg, b, a, t, s), 0) + np.conj(v)
    return out


def _write_spin_transfer(d, norb, entries):
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(MARKED, norb))
        for (r, a, b, s, t), v in _closed(entries).items():
            f.write("{} {} {} {} {} {} {} {:.12f} {:.12f}\n".format(
                *r, a, b, s, t, v.real, v.imag))


def _write_plain_transfer(d, norb, entries):
    """{(R, a, b): v} written as 7 columns (Hermitian-closed)."""
    closed = {}
    for (r, a, b), v in entries.items():
        closed[(r, a, b)] = closed.get((r, a, b), 0) + v
        neg = tuple(-x for x in r)
        if (neg, b, a) == (r, a, b):
            continue
        closed[(neg, b, a)] = closed.get((neg, b, a), 0) + np.conj(v)
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n{}\n1\n 1\n".format(norb))
        for (r, a, b), v in closed.items():
            f.write("{} {} {} {} {} {:.12f} {:.12f}\n".format(
                *r, a, b, v.real, v.imag))


def _write_coulomb_intra(d, norb, u):
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n{}\n1\n 1\n".format(norb))
        for a in range(1, norb + 1):
            f.write("0 0 0 {} {} {:.6f} 0.0\n".format(a, a, u))


def _solver(d, cellshape, *, subshape=(1, 1, 1), twosz=None, ncond=2,
            enable_spin_orbital=False, coulomb=False, iteration_max=1):
    # SubShape defaults to CellShape (fully folded) in this fork, so the
    # unfolded cases say so explicitly
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.uhfk as uhfk_module
    files = {'path_to_input': d, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat'}
    if coulomb:
        files['CoulombIntra'] = 'coulombintra.dat'
    read_io = read_input_k.QLMSkInput({'path_to_input': d,
                                       'interaction': files})
    param = {'Ncond': ncond, 'T': 0.2, 'CellShape': cellshape,
             'IterationMax': iteration_max, 'EPS': 8, 'Mix': 0.5,
             'RndSeed': 1, 'strict_hermite': True}
    param['SubShape'] = list(subshape)
    if twosz is not None:
        param['2Sz'] = twosz
    info_mode = {'mode': 'UHFk', 'param': param,
                 'enable_spin_orbital': enable_spin_orbital}
    solver = uhfk_module.UHFk(read_io.get_param("ham"),
                              {"print_level": 0, "print_step": 100}, info_mode)
    return solver, read_io.get_param("green")


def _prepare(solver):
    """The solve() prologue up to the block detection."""
    solver._make_ham_trans()
    solver._make_ham_inter()
    solver._detect_blocks()


# two-orbital chain: on-site splitting, intra- and inter-orbital hopping
_ORB = {((0, 0, 0), 1, 1): -0.3, ((0, 0, 0), 2, 2): 0.3,
        ((1, 0, 0), 1, 1): -1.0, ((1, 0, 0), 2, 2): -0.7,
        ((1, 0, 0), 1, 2): -0.2 + 0.05j}


def _spin_independent(orb):
    return {(r, a, b, s, s): v for (r, a, b), v in orb.items() for s in (0, 1)}


class TestSpinIndependentEquivalence(unittest.TestCase):
    def setUp(self):
        self.d7 = tempfile.mkdtemp(prefix="w90spin_7_")
        self.d9 = tempfile.mkdtemp(prefix="w90spin_9_")
        for d in (self.d7, self.d9):
            _write_geom(d, 2)
            _write_coulomb_intra(d, 2, 2.0)
        _write_plain_transfer(self.d7, 2, _ORB)
        _write_spin_transfer(self.d9, 2, _spin_independent(_ORB))

    def tearDown(self):
        for d in (self.d7, self.d9):
            shutil.rmtree(d, ignore_errors=True)

    def test_case1_same_hk_and_same_converged_result(self):
        s7, g7 = _solver(self.d7, [4, 1, 1], coulomb=True, iteration_max=200)
        s9, g9 = _solver(self.d9, [4, 1, 1], coulomb=True, iteration_max=200)
        self.assertFalse(s7.spin_resolved_transfer)
        self.assertTrue(s9.spin_resolved_transfer)
        s7._make_ham_trans()
        s9._make_ham_trans()
        np.testing.assert_allclose(s9.ham_trans, s7.ham_trans, atol=1e-12)
        with tempfile.TemporaryDirectory() as o7, \
                tempfile.TemporaryDirectory() as o9:
            s7.solve(g7, o7)
            s9.solve(g9, o9)
        self.assertAlmostEqual(s9.physics["Ene"]["Total"],
                               s7.physics["Ene"]["Total"], places=10)
        np.testing.assert_allclose(s9.Green, s7.Green, atol=1e-10)

    def test_case2_fixed_2sz_accepted(self):
        s9, _ = _solver(self.d9, [4, 1, 1], coulomb=True, twosz=0)
        _prepare(s9)
        self.assertFalse(s9.sz_free)
        # two pure-spin groups, one electron each
        self.assertEqual(sorted(s9.group_nconds), [1, 1])
        self.assertEqual(len(set(s9.block_to_group)), 2)

    def test_case5_sublattice_fold_matches_plain(self):
        s7, _ = _solver(self.d7, [4, 1, 1], subshape=[2, 1, 1])
        s9, _ = _solver(self.d9, [4, 1, 1], subshape=[2, 1, 1])
        s7._make_ham_trans()
        s9._make_ham_trans()
        self.assertEqual(s9.ham_trans.shape, (2, 8, 8))
        np.testing.assert_allclose(s9.ham_trans, s7.ham_trans, atol=1e-12)


class TestSpinDependent(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="w90spin_dep_")
        _write_geom(self.d, 2)
        _write_coulomb_intra(self.d, 2, 2.0)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_case3_zeeman_accepts_2sz_and_keeps_sectors(self):
        h = 0.4
        entries = _spin_independent(_ORB)
        entries[((0, 0, 0), 1, 1, 0, 0)] += -h    # up: -h on orbital 1
        entries[((0, 0, 0), 1, 1, 1, 1)] += +h    # dn: +h
        _write_spin_transfer(self.d, 2, entries)
        s, _ = _solver(self.d, [4, 1, 1], coulomb=True, twosz=2)
        _prepare(s)
        self.assertEqual(sorted(s.group_nconds), [0, 2])
        # the field sits on the spin-diagonal blocks only
        s._make_ham_trans()
        np.testing.assert_allclose(s.ham_trans[:, :2, 2:], 0.0, atol=1e-12)
        self.assertAlmostEqual(
            (s.ham_trans[0, 0, 0] - s.ham_trans[0, 2, 2]).real, -2 * h, places=12)

    def test_case3_spin_flip_refuses_2sz(self):
        entries = _spin_independent(_ORB)
        entries[((0, 0, 0), 1, 2, 0, 1)] = 0.15 + 0.1j     # on-site SOC
        _write_spin_transfer(self.d, 2, entries)
        s, _ = _solver(self.d, [4, 1, 1], coulomb=True, twosz=0)
        with self.assertRaisesRegex(ValueError, "spin-flip"):
            _prepare(s)
        # and the same file runs Sz-free
        s, _ = _solver(self.d, [4, 1, 1], coulomb=True)
        _prepare(s)
        self.assertTrue(s.sz_free)

    def test_strict_hermite_names_orbital_and_spin(self):
        # a spin-flip entry without its Hermitian partner: the strict check
        # must fail and name the entry as (orbital, spin), not as the
        # spin-block composite index
        entries = _spin_independent(_ORB)
        _write_spin_transfer(self.d, 2, entries)
        with open(os.path.join(self.d, "transfer.dat"), "a") as f:
            f.write("0 0 0 1 2 0 1 0.150000000000 0.100000000000\n")
        with self.assertLogs("qlms", level="ERROR") as cm, \
                self.assertRaises(SystemExit):
            _solver(self.d, [4, 1, 1])
        text = "\n".join(cm.output)
        self.assertIn("Hermiticity check failed", text)
        self.assertIn("(orbital, spin) = (1, 0), (2, 1)", text)

    def test_case5_fold_keeps_down_block(self):
        entries = _spin_independent(_ORB)
        entries[((1, 0, 0), 1, 1, 1, 1)] = -0.5         # dn hops differently
        _write_spin_transfer(self.d, 2, entries)
        s, _ = _solver(self.d, [4, 1, 1], subshape=[2, 1, 1])
        s._make_ham_trans()
        norb = s.norb                                    # 4 after the fold
        up = s.ham_trans[:, :norb, :norb]
        dn = s.ham_trans[:, norb:, norb:]
        self.assertGreater(np.abs(up - dn).max(), 0.1)
        np.testing.assert_allclose(s.ham_trans[:, :norb, norb:], 0.0, atol=1e-12)


class TestSOCEquivalenceWithSpinOrbitalMode(unittest.TestCase):
    """Case 4: a spin-resolved SOC file and the interleaved spin-orbital
    file of the same model give the same H(k) up to the index order."""

    def setUp(self):
        self.d9 = tempfile.mkdtemp(prefix="w90spin_soc9_")
        self.dso = tempfile.mkdtemp(prefix="w90spin_socso_")

    def tearDown(self):
        for d in (self.d9, self.dso):
            shutil.rmtree(d, ignore_errors=True)

    def test_same_hk_after_reorder(self):
        norb = 2
        entries = _spin_independent(_ORB)
        entries[((0, 0, 0), 1, 2, 0, 1)] = 0.15 + 0.1j      # on-site SOC
        entries[((1, 0, 0), 2, 1, 1, 0)] = 0.08              # Rashba-like
        _write_geom(self.d9, norb)
        _write_spin_transfer(self.d9, norb, entries)
        # the same table folded to interleaved indices 2*orb + spin, with
        # the doubled geometry spin-orbital mode expects
        _write_geom(self.dso, 2 * norb)
        with open(os.path.join(self.dso, "transfer.dat"), "w") as f:
            f.write("Transfer SOI format\n{}\n1\n 1\n".format(2 * norb))
            for (r, a, b, s, t), v in _closed(entries).items():
                ia = spin_orbital.interleaved_index(a - 1, s) + 1
                ib = spin_orbital.interleaved_index(b - 1, t) + 1
                f.write("{} {} {} {} {} {:.12f} {:.12f}\n".format(
                    *r, ia, ib, v.real, v.imag))
        s9, _ = _solver(self.d9, [4, 1, 1])
        sso, _ = _solver(self.dso, [4, 1, 1], enable_spin_orbital=True)
        s9._make_ham_trans()
        sso._make_ham_trans()
        np.testing.assert_allclose(
            s9.ham_trans,
            spin_orbital.interleaved_to_spinblock(sso.ham_trans, norb),
            atol=1e-12)


if __name__ == "__main__":
    unittest.main()
