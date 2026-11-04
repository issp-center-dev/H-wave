#!/usr/bin/env python3
"""RPA one-body construction through hamiltonian_k (refactoring plan 2.4).

Interaction._make_ham_trans and the Extern build now go through
hamiltonian_k; a normal-mode spin-orbital-format file has its spin block
dropped once at the boundary (Interaction.__init__), with one warning, in
every geometry -- it used to be skipped silently (issue #303) -- and the
spin-orbital branch refuses an out-of-range index with the builder's
message.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np


def _count_ignore_warnings(records):
    return sum(1 for r in records if "These terms are ignored" in r.getMessage())


def _write_chain(d, norb_file, rows, geom_norb=1, extern_rows=None):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(geom_norb))
        for i in range(geom_norb):
            f.write("{} 0 0\n".format(0.1 * i))
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("t\n{}\n1\n 1\n".format(norb_file))
        for row in rows:
            f.write(row + "\n")
    if extern_rows is not None:
        with open(os.path.join(d, "extern.dat"), "w") as f:
            f.write("h\n{}\n1\n 1\n".format(norb_file))
            for row in extern_rows:
                f.write(row + "\n")


def _rpa(d, subshape, *, enable_spin_orbital=False, extern=False,
         coeff_extern=0.0):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as rpa_module
    files = {'path_to_input': d, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat'}
    if extern:
        files['Extern'] = 'extern.dat'
    read_io = read_input_k.QLMSkInput({'path_to_input': d,
                                       'interaction': files})
    param = {'T': 0.5, 'Ncond': 1, 'CellShape': [4, 1, 1],
             'SubShape': subshape, 'Nmat': 16}
    if coeff_extern:
        param['coeff_extern'] = coeff_extern
    info_mode = {'mode': 'RPA', 'calc_scheme': 'general',
                 'enable_spin_orbital': enable_spin_orbital, 'param': param}
    return rpa_module.RPA(read_io.get_param("ham"),
                          {"print_level": 0, "print_step": 100}, info_mode)


_CHAIN = [" 1 0 0 1 1 -1.0 0.0", "-1 0 0 1 1 -1.0 0.0"]
_SPIN_BLOCK = [" 1 0 0 2 2 0.7 0.0", "-1 0 0 2 2 0.7 0.0"]


class TestNormalModeSpinBlockDroppedAtBoundary(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="rpa_hk_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _check(self, subshape):
        _write_chain(self.d, 2, _CHAIN + _SPIN_BLOCK)
        with self.assertLogs("qlms", level="WARNING") as cm:
            s = _rpa(self.d, subshape)
        self.assertEqual(_count_ignore_warnings(cm.records), 1)
        ham = s.ham_info
        # the pre-fold table no longer carries the spin block
        src = ham.param_ham_orig if s.lattice.has_sublattice else ham.param_ham
        self.assertTrue(all(o == 0 for (_, ov) in src["Transfer"] for o in ov))
        # H0(k) equals the one built from the same chain written without
        # the spin-block rows (a plain 7-column file, num_wann = 1)
        d_ref = tempfile.mkdtemp(prefix="rpa_hk_ref_")
        try:
            _write_chain(d_ref, 1, _CHAIN)
            ref = _rpa(d_ref, subshape).ham_info
        finally:
            shutil.rmtree(d_ref, ignore_errors=True)
        self.assertEqual(ham.ham_trans_q.shape, ref.ham_trans_q.shape)
        np.testing.assert_allclose(ham.ham_trans_q, ref.ham_trans_q, atol=1e-12)
        np.testing.assert_allclose(ham.ham_trans_r, ref.ham_trans_r, atol=1e-12)
        self.assertGreater(np.abs(ham.ham_trans_q).max(), 1.0)   # the chain is there

    def test_unfolded_warns_once(self):
        self._check([1, 1, 1])

    def test_folded_warns_once(self):
        self._check([2, 1, 1])

    def test_extern_spin_block_dropped_once(self):
        _write_chain(self.d, 2, _CHAIN,
                     extern_rows=[" 0 0 0 1 1 0.25 0.0", " 0 0 0 2 2 0.9 0.0"])
        with self.assertLogs("qlms", level="WARNING") as cm:
            s = _rpa(self.d, [1, 1, 1], extern=True, coeff_extern=1.0)
        self.assertEqual(_count_ignore_warnings(cm.records), 1)
        np.testing.assert_allclose(s.ham_info.ham_extern_q[:, 0, 0], 0.25,
                                   atol=1e-12)


class TestSpinOrbitalBranch(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="rpa_hk_so_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def test_out_of_range_index_raises(self):
        # geometry norb = 2 (one physical orbital); index 3 is outside
        _write_chain(self.d, 2, _CHAIN + [" 0 0 0 3 3 0.5 0.0"], geom_norb=2)
        with self.assertRaisesRegex(ValueError, r"outside \[0, 2\)"):
            _rpa(self.d, [1, 1, 1], enable_spin_orbital=True)

    def test_interleaved_file_lands_in_spin_block_order(self):
        # two physical orbitals: interleaved index 3 (1-based) is
        # (orbital 1, up) and must sit at spin-block position 1
        rows = [" 0 0 0 3 3 0.4 0.0", " 0 0 0 1 1 -0.4 0.0",
                " 0 0 0 2 2 -0.4 0.0", " 0 0 0 4 4 0.4 0.0"]
        _write_chain(self.d, 4, rows, geom_norb=4)
        s = _rpa(self.d, [1, 1, 1], enable_spin_orbital=True)
        h = s.ham_info.ham_trans_q[0]
        np.testing.assert_allclose(np.diag(h).real, [-0.4, 0.4, -0.4, 0.4],
                                   atol=1e-12)


if __name__ == "__main__":
    unittest.main()
