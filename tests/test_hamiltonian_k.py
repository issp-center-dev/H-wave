#!/usr/bin/env python3
"""hwave.solver.hamiltonian_k (refactoring plan 2.1): unit behaviour, and
element-by-element equivalence with the three solver copies it will
replace -- UHFk._make_ham_trans (normal and spin-orbital mode, with and
without a sublattice fold), rpa.Interaction._make_ham_trans (normal and
spin-orbital) and sc._build_hamiltonian_k.

The solver copies are untouched by the plan-2.1 PR; these tests are the
gate for switching them (2.2-2.4).
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver import hamiltonian_k as hk
from hwave.solver import spin_orbital


def _sb(orb, spin, norb):
    """Spin-block composite index spin*norb + orb (the 9-column fold, #323)."""
    return spin * norb + orb


def _hermitian_table(rng, norb, bonds):
    """Random table with both members of every Hermitian pair
    ((R, a, b), (-R, b, a)) present; ``bonds`` are the R vectors to use."""
    tbl = {}
    for irvec in bonds:
        for a in range(norb):
            for b in range(norb):
                v = complex(rng.normal(), rng.normal())
                neg = tuple(-r for r in irvec)
                if irvec == (0, 0, 0) and a == b:
                    v = complex(v.real, 0.0)
                tbl[(irvec, (a, b))] = tbl.get((irvec, (a, b)), 0) + v
                tbl[(neg, (b, a))] = tbl.get((neg, (b, a)), 0) + np.conj(v)
    return tbl


def _explicit_k(table, shape, dim):
    """sum_R e^{+ikR} T(R) by hand on the FFT grid, C-ordered k."""
    nx, ny, nz = shape
    out = np.zeros((nx, ny, nz, dim, dim), dtype=np.complex128)
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for (irvec, (a, b)), v in table.items():
                    phase = 2 * np.pi * (ix * irvec[0] / nx + iy * irvec[1] / ny
                                         + iz * irvec[2] / nz)
                    out[ix, iy, iz, a, b] += v * np.exp(1j * phase)
    return out.reshape(nx * ny * nz, dim, dim)


class TestDenseTable(unittest.TestCase):
    def test_negative_r_wraps_like_numpy_indexing(self):
        tab = hk.dense_table({((-1, 0, 0), (0, 0)): 2.0}, (4, 1, 1), 1)
        self.assertEqual(tab[3, 0, 0, 0, 0], 2.0)
        self.assertEqual(np.count_nonzero(tab), 1)

    def test_wraparound_images_accumulate(self):
        # R = +1 and R = -3 are the same slot on a 4-site ring
        tab = hk.dense_table({((1, 0, 0), (0, 0)): 1.0,
                              ((-3, 0, 0), (0, 0)): 0.5}, (4, 1, 1), 1)
        self.assertEqual(tab[1, 0, 0, 0, 0], 1.5)

    def test_out_of_range_raises_by_default(self):
        with self.assertRaisesRegex(ValueError, r"outside \[0, 1\)"):
            hk.dense_table({((0, 0, 0), (0, 1)): 1.0}, (2, 1, 1), 1)

    def test_out_of_range_drop_warns_once(self):
        tbl = {((0, 0, 0), (0, 0)): -1.0, ((0, 0, 0), (1, 1)): 0.5,
               ((0, 0, 0), (0, 1)): 0.2}
        with self.assertLogs("qlms", level="WARNING") as cm:
            tab = hk.dense_table(tbl, (2, 1, 1), 1, out_of_range="drop")
        self.assertEqual(len(cm.output), 1)
        self.assertIn("2 entries", cm.output[0])
        self.assertEqual(tab[0, 0, 0, 0, 0], -1.0)
        self.assertEqual(np.count_nonzero(tab), 1)

    def test_bad_policy_rejected(self):
        with self.assertRaises(ValueError):
            hk.dense_table({}, (1, 1, 1), 1, out_of_range="ignore")


class TestToKAndEmbed(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(12)

    def test_to_k_is_plus_ikr_sum_on_fft_grid(self):
        shape, norb = (4, 3, 1), 2
        tbl = _hermitian_table(self.rng, norb,
                               [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, -1, 0)])
        got = hk.to_k(hk.dense_table(tbl, shape, norb))
        np.testing.assert_allclose(got, _explicit_k(tbl, shape, norb),
                                   atol=1e-12)

    def test_embed_spin_is_kron_identity(self):
        m = self.rng.normal(size=(3, 2, 2)) + 1j * self.rng.normal(size=(3, 2, 2))
        got = hk.embed_spin(m)
        for k in range(3):
            np.testing.assert_allclose(got[k], np.kron(np.eye(2), m[k]))

    def test_spin_block_table_equals_embedded_orbital_table(self):
        # a spin-independent table written with explicit spins (the
        # spin-resolved 9-column Transfer, folded to s*norb + a) must build
        # the same H(k) as the 7-column table embedded over spin
        shape, norb = (4, 1, 1), 2
        orb = _hermitian_table(self.rng, norb, [(0, 0, 0), (1, 0, 0)])
        spin = {}
        for (irvec, (a, b)), v in orb.items():
            for s in (0, 1):
                spin[(irvec, (_sb(a, s, norb), _sb(b, s, norb)))] = v
        np.testing.assert_allclose(
            hk.one_body_k(spin, shape, 2 * norb),
            hk.one_body_k(orb, shape, norb, spin_embed=True), atol=1e-12)

    def test_spin_flip_entry_lands_in_off_diagonal_block(self):
        norb = 2
        key = ((0, 0, 0), (_sb(1, 0, norb), _sb(0, 1, norb)))
        h = hk.one_body_k({key: 0.3 + 0.1j}, (1, 1, 1), 2 * norb)
        self.assertEqual(h[0, 1, 2], 0.3 + 0.1j)   # (a=1, up) -> (a=0, dn)
        self.assertEqual(np.count_nonzero(h), 1)


# ---------------------------------------------------------------------------
# equivalence with the solver copies
# ---------------------------------------------------------------------------

def _write_geom(d, norb):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(norb))
        for i in range(norb):
            f.write("{} 0 0\n".format(0.1 * i))


def _write_transfer(d, norb_file, table, extra_rows=()):
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n{}\n1\n 1\n".format(norb_file))
        for (irvec, (a, b)), v in table.items():
            f.write("{} {} {} {} {} {:.12f} {:.12f}\n".format(
                *irvec, a + 1, b + 1, v.real, v.imag))
        for row in extra_rows:
            f.write(row + "\n")


def _input(d):
    import hwave.qlmsio.read_input_k as read_input_k
    return read_input_k.QLMSkInput({'path_to_input': d, 'interaction': {
        'path_to_input': d, 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat'}}).get_param("ham")


class TestEquivalenceUHFk(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(7)
        self.d = tempfile.mkdtemp(prefix="hk_uhfk_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _solver(self, cellshape, subshape, enable_spin_orbital):
        import hwave.solver.uhfk as uhfk_module
        param = {'Ncond': 2, 'T': 0.5, 'CellShape': cellshape,
                 'SubShape': subshape, 'IterationMax': 1, 'EPS': 8,
                 'Mix': 0.5, 'RndSeed': 1}
        info_mode = {'mode': 'UHFk', 'param': param,
                     'enable_spin_orbital': enable_spin_orbital}
        s = uhfk_module.UHFk(_input(self.d), {"print_level": 0,
                                              "print_step": 100}, info_mode)
        s._make_ham_trans()
        return s

    def _check_normal(self, subshape):
        norb = 2
        tbl = _hermitian_table(self.rng, norb, [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        _write_geom(self.d, norb)
        # a spin-block row (index 3 > norb) that normal mode drops
        _write_transfer(self.d, norb, tbl, extra_rows=["0 0 0 3 3 0.7 0.0"])
        s = self._solver([4, 3, 1], subshape, False)
        # UHFk drops the spin-block row once at the boundary
        # (_drop_spin_block_transfer), so the table the builder sees is
        # in range in both geometries
        got = hk.one_body_k(s.param_ham["Transfer"], s.shape, s.norb,
                            spin_embed=True)
        np.testing.assert_allclose(got, s.ham_trans, atol=1e-12)

    def test_normal_mode_no_fold(self):
        self._check_normal([1, 1, 1])

    def test_normal_mode_folded(self):
        self._check_normal([2, 1, 1])

    def _check_spin_orbital(self, subshape):
        nd = 4     # two physical orbitals, interleaved 2*orb + spin
        tbl = _hermitian_table(self.rng, nd, [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        _write_geom(self.d, nd)
        _write_transfer(self.d, nd, tbl)
        s = self._solver([4, 3, 1], subshape, True)
        got = hk.one_body_k(s.param_ham["Transfer"], s.shape, s.nd)
        np.testing.assert_allclose(got, s.ham_trans, atol=1e-12)

    def test_spin_orbital_mode_no_fold(self):
        self._check_spin_orbital([1, 1, 1])

    def test_spin_orbital_mode_folded(self):
        self._check_spin_orbital([2, 1, 1])


class TestEquivalenceRPA(unittest.TestCase):
    def setUp(self):
        self.rng = np.random.default_rng(3)
        self.d = tempfile.mkdtemp(prefix="hk_rpa_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _solver(self, cellshape, subshape, enable_spin_orbital):
        import hwave.solver.rpa as rpa_module
        info_mode = {'mode': 'RPA', 'calc_scheme': 'general',
                     'enable_spin_orbital': enable_spin_orbital,
                     'param': {'T': 0.5, 'Ncond': 2, 'CellShape': cellshape,
                               'SubShape': subshape, 'Nmat': 16}}
        return rpa_module.RPA(_input(self.d), {"print_level": 0,
                                               "print_step": 100}, info_mode)

    def _check_normal(self, subshape):
        norb = 2
        tbl = _hermitian_table(self.rng, norb, [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        _write_geom(self.d, norb)
        _write_transfer(self.d, norb, tbl, extra_rows=["0 0 0 3 3 0.7 0.0"])
        s = self._solver([4, 3, 1], subshape, False)
        ham = s.ham_info
        shape = s.lattice.shape
        if subshape == [1, 1, 1]:
            with self.assertLogs("qlms", level="WARNING"):
                tab_r = hk.dense_table(ham.param_ham["Transfer"], shape,
                                       ham.norb, out_of_range="drop")
        else:
            tab_r = hk.dense_table(ham.param_ham["Transfer"], shape, ham.norb)
        np.testing.assert_allclose(tab_r.reshape(ham.ham_trans_r.shape),
                                   ham.ham_trans_r, atol=1e-12)
        np.testing.assert_allclose(hk.to_k(tab_r), ham.ham_trans_q, atol=1e-12)

    def test_normal_mode_no_fold(self):
        self._check_normal([1, 1, 1])

    def test_normal_mode_folded(self):
        self._check_normal([2, 1, 1])

    def _check_spin_orbital(self, subshape):
        nd = 4
        tbl = _hermitian_table(self.rng, nd, [(0, 0, 0), (1, 0, 0), (0, 1, 0)])
        _write_geom(self.d, nd)
        _write_transfer(self.d, nd, tbl)
        s = self._solver([4, 3, 1], subshape, True)
        ham = s.ham_info
        # RPA works in spin-block order: reorder after the layout-agnostic
        # build, exactly where the solver copy remaps each index
        got = spin_orbital.interleaved_to_spinblock(
            hk.one_body_k(ham.param_ham["Transfer"], s.lattice.shape, 2 * ham.norb),
            ham.norb)
        np.testing.assert_allclose(got, ham.ham_trans_q, atol=1e-12)

    def test_spin_orbital_mode_no_fold(self):
        self._check_spin_orbital([1, 1, 1])

    def test_spin_orbital_mode_folded(self):
        self._check_spin_orbital([2, 1, 1])


class TestEquivalenceSC(unittest.TestCase):
    def test_build_hamiltonian_k_matches(self):
        from hwave import sc
        rng = np.random.default_rng(5)
        shape, norb = (4, 3, 2), 3
        tbl = _hermitian_table(rng, norb, [(0, 0, 0), (1, 0, 0), (0, -1, 1),
                                           (2, 1, 0)])
        karr = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in shape]
        eps = sc._build_hamiltonian_k(*karr, tbl, norb)       # (p, p, nx, ny, nz)
        got = hk.to_k(hk.dense_table(tbl, shape, norb)).reshape(
            *shape, norb, norb).transpose(3, 4, 0, 1, 2)
        np.testing.assert_allclose(got, eps, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
