#!/usr/bin/env python3
"""Shared sublattice fold of R-space Green / Hamiltonian arrays.

``fold.fold_layout`` / ``fold.unfold_layout`` are the one definition of the
layout change between the original lattice ``(Lvol, ns, norb_orig, ns,
norb_orig)`` and the sublattice (supercell) basis ``(Nvol, ns, norb, ns,
norb)`` with ``norb = norb_orig * Bvol``.  They must reproduce, bit for bit:

* the d72f35f6 ``UHFk._reshape_green`` / ``UHFk._deflate_green`` loops (frozen
  verbatim in ``tests/uhfk_legacy_d72f35f6.py``), both conventions;
* ``RPA._reshape_green`` (both conventions), the solver's own vectorised copy;
* the folded Transfer table that ``fold.reshape_interaction`` builds
  (``hamiltonian=True`` on a primitive H(R)).
"""
import types
import unittest

import numpy as np

from hwave.solver import fold
from hwave.solver import hamiltonian_k
from hwave.solver.rpa import RPA
from hwave.solver.uhfk import UHFk

from tests import uhfk_legacy_d72f35f6 as legacy

# (CellShape, SubShape, norb_orig, ns); the last one is shaped like the
# enable_spin_orbital layout (ns = 1, spin folded into the orbital index).
CASES = [
    ((4, 4, 1), (2, 2, 1), 1, 2),
    ((4, 3, 1), (2, 3, 1), 2, 2),
    ((6, 4, 2), (3, 2, 2), 1, 2),
    ((2, 2, 1), (2, 2, 1), 4, 1),
]


def _layout(cellshape, subshape, norb_orig, ns):
    shape = tuple(L // B for L, B in zip(cellshape, subshape))
    bvol = int(np.prod(subshape))
    return types.SimpleNamespace(
        cellshape=tuple(cellshape), subshape=tuple(subshape), shape=shape,
        norb_orig=norb_orig, norb=norb_orig * bvol, ns=ns,
        lvol=int(np.prod(cellshape)), nvol=int(np.prod(shape)))


def _random_distinct(rng, shape):
    g = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    # distinct entries, so a wrong index cannot hide behind equal values
    assert np.unique(g).size == g.size
    return g


def _kw(lay):
    return dict(cellshape=lay.cellshape, subshape=lay.subshape,
                norb_orig=lay.norb_orig, ns=lay.ns)


def _uhfk_stub(lay):
    solver = UHFk.__new__(UHFk)
    for key in ("cellshape", "subshape", "shape", "norb_orig", "norb", "ns"):
        setattr(solver, key, getattr(lay, key))
    return solver


def _rpa_stub(lay):
    stub = types.SimpleNamespace(
        lattice=types.SimpleNamespace(cellshape=lay.cellshape,
                                      subshape=lay.subshape,
                                      shape=lay.shape),
        norb_orig=lay.norb_orig, norb=lay.norb, ns=lay.ns)
    return stub


def _assert_bitwise(test, a, b):
    test.assertEqual(a.shape, b.shape)
    test.assertEqual(a.dtype, b.dtype)
    # array_equal is value equality; compare the bytes too, so -0.0 vs +0.0
    # or a dtype change cannot pass
    test.assertTrue(np.array_equal(a, b))
    test.assertEqual(np.ascontiguousarray(a).tobytes(),
                     np.ascontiguousarray(b).tobytes())


class TestFoldLayoutMatchesLegacy(unittest.TestCase):
    def test_fold_matches_reshape_green_bitwise(self):
        rng = np.random.default_rng(364)
        for case in CASES:
            lay = _layout(*case)
            with self.subTest(case=case):
                g = _random_distinct(
                    rng, (lay.lvol, lay.ns, lay.norb_orig, lay.ns,
                          lay.norb_orig))
                got = fold.fold_layout(g, **_kw(lay))
                _assert_bitwise(self, got, legacy.reshape_green(lay, g))

    def test_unfold_matches_deflate_green_bitwise(self):
        rng = np.random.default_rng(365)
        for case in CASES:
            lay = _layout(*case)
            gs = _random_distinct(
                rng, (lay.nvol, lay.ns, lay.norb, lay.ns, lay.norb))
            for ham in (False, True):
                with self.subTest(case=case, hamiltonian=ham):
                    got = fold.unfold_layout(gs, hamiltonian=ham, **_kw(lay))
                    _assert_bitwise(
                        self, got,
                        legacy.deflate_green(lay, gs, hamiltonian=ham))

    def test_fold_matches_rpa_reshape_green_bitwise(self):
        # the Hamiltonian convention has no UHFk loop; RPA's copy carries it
        rng = np.random.default_rng(366)
        for case in CASES:
            lay = _layout(*case)
            g = _random_distinct(
                rng, (lay.lvol, lay.ns, lay.norb_orig, lay.ns, lay.norb_orig))
            for ham in (False, True):
                with self.subTest(case=case, hamiltonian=ham):
                    got = fold.fold_layout(g, hamiltonian=ham, **_kw(lay))
                    ref = RPA._reshape_green(_rpa_stub(lay), g,
                                             hamiltonian=ham)
                    self.assertTrue(np.array_equal(got, ref))

    def test_uhfk_methods_delegate_bitwise(self):
        rng = np.random.default_rng(367)
        for case in CASES:
            lay = _layout(*case)
            solver = _uhfk_stub(lay)
            g = _random_distinct(
                rng, (lay.lvol, lay.ns, lay.norb_orig, lay.ns, lay.norb_orig))
            gs = _random_distinct(
                rng, (lay.nvol, lay.ns, lay.norb, lay.ns, lay.norb))
            with self.subTest(case=case):
                _assert_bitwise(self, solver._reshape_green(g),
                                legacy.reshape_green(lay, g))
                _assert_bitwise(self, solver._deflate_green(gs),
                                legacy.deflate_green(lay, gs))
                _assert_bitwise(
                    self, solver._deflate_green(gs, hamiltonian=True),
                    legacy.deflate_green(lay, gs, hamiltonian=True))

    def test_real_input_matches_legacy(self):
        # the legacy loops write into a complex128 buffer; a real input
        # (e.g. a real transfer table) must come out the same way
        rng = np.random.default_rng(368)
        lay = _layout((4, 3, 1), (2, 3, 1), 2, 2)
        g = rng.standard_normal(
            (lay.lvol, lay.ns, lay.norb_orig, lay.ns, lay.norb_orig))
        gs = rng.standard_normal(
            (lay.nvol, lay.ns, lay.norb, lay.ns, lay.norb))
        _assert_bitwise(self, fold.fold_layout(g, **_kw(lay)),
                        legacy.reshape_green(lay, g))
        for ham in (False, True):
            _assert_bitwise(
                self, fold.unfold_layout(gs, hamiltonian=ham, **_kw(lay)),
                legacy.deflate_green(lay, gs, hamiltonian=ham))


class TestFoldLayoutProperties(unittest.TestCase):
    def test_round_trip(self):
        rng = np.random.default_rng(369)
        for case in CASES:
            lay = _layout(*case)
            g = _random_distinct(
                rng, (lay.lvol, lay.ns, lay.norb_orig, lay.ns, lay.norb_orig))
            for ham in (False, True):
                with self.subTest(case=case, hamiltonian=ham):
                    gs = fold.fold_layout(g, hamiltonian=ham, **_kw(lay))
                    back = fold.unfold_layout(gs, hamiltonian=ham,
                                              **_kw(lay))
                    _assert_bitwise(self, back, g)
                    # and the fold of the unfold of a folded array is itself:
                    # fold(hamiltonian=True) is the inverse of
                    # unfold(hamiltonian=True) on sublattice-periodic input
                    _assert_bitwise(
                        self,
                        fold.fold_layout(back, hamiltonian=ham, **_kw(lay)),
                        gs)

    def test_conventions_differ_off_diagonal(self):
        # guard against both flags silently mapping to one layout
        rng = np.random.default_rng(370)
        lay = _layout((4, 4, 1), (2, 2, 1), 1, 2)
        g = _random_distinct(
            rng, (lay.lvol, lay.ns, lay.norb_orig, lay.ns, lay.norb_orig))
        a = fold.fold_layout(g, hamiltonian=False, **_kw(lay))
        b = fold.fold_layout(g, hamiltonian=True, **_kw(lay))
        self.assertFalse(np.array_equal(a, b))
        # they agree on the in-cell diagonal blocks (ri == rj)
        n = lay.norb_orig
        for ir in range(lay.norb // n):
            sl = slice(ir * n, (ir + 1) * n)
            self.assertTrue(np.array_equal(a[:, :, sl, :, sl],
                                           b[:, :, sl, :, sl]))

    def test_hamiltonian_fold_matches_folded_transfer(self):
        # fold_layout(h_r, hamiltonian=True) of the primitive Transfer table
        # is the supercell Transfer table fold.reshape_interaction builds
        # (the call UHFk._reshape_interaction makes for the Transfer).
        rng = np.random.default_rng(371)
        for cellshape, subshape, norb in [((6, 4, 2), (3, 2, 2), 1),
                                          ((4, 4, 1), (2, 2, 1), 2),
                                          ((6, 4, 2), (3, 2, 2), 2)]:
            lay = _layout(cellshape, subshape, norb, 1)
            # one entry per distinct R mod L (no wrap-around collision, so
            # the accumulation order cannot differ), every orbital pair
            table = {}
            for rx in range(-(cellshape[0] // 2) + 1, cellshape[0] // 2 + 1):
                for ry in range(-(cellshape[1] // 2) + 1,
                                cellshape[1] // 2 + 1):
                    for rz in range(cellshape[2]):
                        for a in range(norb):
                            for b in range(norb):
                                table[((rx, ry, rz), (a, b))] = complex(
                                    rng.standard_normal(),
                                    rng.standard_normal())
            with self.subTest(cellshape=cellshape, subshape=subshape,
                              norb=norb):
                h_r = hamiltonian_k.dense_table(table, cellshape, norb)
                h_r = h_r.reshape(lay.lvol, 1, norb, 1, norb)
                got = fold.fold_layout(h_r, hamiltonian=True, **_kw(lay))

                folded = fold.reshape_interaction(
                    table, lay.subshape, lay.shape, norb_so_orig=norb,
                    norb_phys_orig=norb, enable_spin_orbital=False,
                    drop_spin_block=True)
                ref = hamiltonian_k.dense_table(folded, lay.shape, lay.norb)
                ref = ref.reshape(lay.nvol, 1, lay.norb, 1, lay.norb)
                _assert_bitwise(self, got, ref)
                # and the Green convention is NOT the Hamiltonian table
                self.assertFalse(np.array_equal(
                    fold.fold_layout(h_r, **_kw(lay)), ref))

    def test_wrong_shape_raises(self):
        lay = _layout((4, 4, 1), (2, 2, 1), 2, 2)
        # an already folded array (norb instead of norb_orig) is refused
        folded = np.zeros((lay.lvol, lay.ns, lay.norb, lay.ns, lay.norb),
                          dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "norb_orig"):
            fold.fold_layout(folded, **_kw(lay))
        with self.assertRaisesRegex(ValueError, "norb_orig"):
            _uhfk_stub(lay)._reshape_green(folded)
        # unfold of an unfolded array is refused too
        orig = np.zeros((lay.lvol, lay.ns, lay.norb_orig, lay.ns,
                         lay.norb_orig), dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, "norb_orig"):
            fold.unfold_layout(orig, **_kw(lay))
        with self.assertRaisesRegex(ValueError, "norb_orig"):
            _uhfk_stub(lay)._deflate_green(orig)

    def test_subshape_must_divide_cellshape(self):
        with self.assertRaisesRegex(ValueError, "subshape"):
            fold.fold_layout(np.zeros((12, 1, 1, 1, 1)), cellshape=(4, 3, 1),
                             subshape=(2, 2, 1), norb_orig=1, ns=1)


if __name__ == "__main__":
    unittest.main()
