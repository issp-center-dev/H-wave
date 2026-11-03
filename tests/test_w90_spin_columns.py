#!/usr/bin/env python3
"""Spin-resolved (9-column) Wannier90-like Transfer format, issue #323, phase 1.

The reader folds a marked 9-column file to spin-block keys (s*norb + a) and
refuses every layout the format rules exclude; the k-space input reader
records the format beside the table; UHFk and RPA refuse the table until
they consume it, and refuse it outright together with enable_spin_orbital.
"""
import os
import shutil
import tempfile
import unittest

from hwave.qlmsio import wan90


MARKED = "Transfer in wannier90-like format for uhfk w90spin/1"
PLAIN = "Transfer in wannier90-like format for uhfk"


def _write(lines, suffix="_tr.dat"):
    fd, path = tempfile.mkstemp(suffix=suffix)
    os.close(fd)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _read(lines, **kw):
    path = _write(lines)
    try:
        return wan90.read_w90(path, **kw)
    finally:
        os.remove(path)


# two physical orbitals, one R-point, all four spin blocks of the (1,2) pair
_TWO_ORB_SPIN = [
    MARKED, "2", "1", "1",
    " 1 0 0  1 2  0 0   0.50 0.25",
    " 1 0 0  1 2  1 1  -0.50 0.00",
    " 1 0 0  1 2  0 1   0.10 0.00",
    " 1 0 0  1 2  1 0   0.10 0.00",
]


class TestReaderSpinColumns(unittest.TestCase):
    def test_marked_file_folds_to_spin_block_keys(self):
        data = _read(_TWO_ORB_SPIN, allow_spin_columns=True)
        # (a, s) -> s*2 + a with a 0-based: (1,2,s,t) -> (0 + 2s, 1 + 2t)
        self.assertEqual(data[((1, 0, 0), (0, 1))], 0.50 + 0.25j)
        self.assertEqual(data[((1, 0, 0), (2, 3))], -0.50)
        self.assertEqual(data[((1, 0, 0), (0, 3))], 0.10)
        self.assertEqual(data[((1, 0, 0), (2, 1))], 0.10)
        self.assertEqual(len(data), 4)

    def test_header_marker_accessor(self):
        path = _write(_TWO_ORB_SPIN)
        try:
            self.assertEqual(wan90.read_w90_header_marker(path), ("w90spin", 1))
        finally:
            os.remove(path)
        path = _write([PLAIN, "1", "1", "1", "0 0 0 1 1 1.0 0.0"])
        try:
            self.assertIsNone(wan90.read_w90_header_marker(path))
        finally:
            os.remove(path)

    def test_marked_file_refused_by_default(self):
        with self.assertRaisesRegex(ValueError, "accepted only for the Transfer"):
            _read(_TWO_ORB_SPIN)

    def test_unknown_marker_version_refused(self):
        lines = [MARKED.replace("w90spin/1", "w90spin/2")] + _TWO_ORB_SPIN[1:]
        with self.assertRaisesRegex(ValueError, "unknown spin-resolved format version"):
            _read(lines, allow_spin_columns=True)

    def test_unmarked_nine_columns_refused(self):
        lines = [PLAIN] + _TWO_ORB_SPIN[1:]
        with self.assertRaisesRegex(ValueError, "no 'w90spin/1' marker"):
            _read(lines, allow_spin_columns=True)

    def test_marked_seven_columns_refused(self):
        lines = [MARKED, "1", "1", "1", "0 0 0 1 1 1.0 0.0"]
        with self.assertRaisesRegex(ValueError, "expected 9"):
            _read(lines, allow_spin_columns=True)

    def test_mixed_rows_refused(self):
        lines = _TWO_ORB_SPIN + [" 0 0 0  1 1  1.0 0.0"]
        with self.assertRaisesRegex(ValueError, "line 9 has 7 columns"):
            _read(lines, allow_spin_columns=True)

    def test_eight_and_ten_columns_refused(self):
        for row in ("0 0 0 1 1 0 1.0 0.0", "0 0 0 1 1 0 0 1.0 0.0 7"):
            with self.assertRaisesRegex(ValueError, "columns, expected"):
                _read([PLAIN, "1", "1", "1", row], allow_spin_columns=True)
        with self.assertRaisesRegex(ValueError, "columns, expected"):
            _read([MARKED, "1", "1", "1", "0 0 0 1 1 0 1.0 0.0"],
                  allow_spin_columns=True)

    def test_spin_index_outside_two_flavors_refused(self):
        lines = [MARKED, "1", "1", "1", "0 0 0 1 1 0 2 1.0 0.0"]
        with self.assertRaisesRegex(ValueError, "N_f = 2"):
            _read(lines, allow_spin_columns=True)

    def test_orbital_index_beyond_header_norb_refused(self):
        lines = [MARKED, "1", "1", "1", "0 0 0 1 2 0 0 1.0 0.0"]
        with self.assertRaisesRegex(ValueError, "outside 1..1"):
            _read(lines, allow_spin_columns=True)

    def test_nonunit_degeneracy_refused_in_spin_format(self):
        lines = [MARKED, "1", "2", "1 2",
                 "0 0 0 1 1 0 0 1.0 0.0", "1 0 0 1 1 0 0 1.0 0.0"]
        with self.assertRaisesRegex(ValueError, "degeneracy"):
            _read(lines, allow_spin_columns=True)

    def test_duplicate_full_key_refused_but_spin_pairs_distinct(self):
        # same (R, a, b) with different (s, t) is fine (the test above);
        # repeating (R, a, b, s, t) is a duplicate
        lines = _TWO_ORB_SPIN + [" 1 0 0  1 2  0 0   0.50 0.25"]
        with self.assertRaises(SystemExit):
            _read(lines, allow_spin_columns=True)

    def test_plain_file_unchanged(self):
        data = _read([PLAIN, "1", "1", "1", "0 0 0 1 1 -1.0 0.5"],
                     allow_spin_columns=True)
        self.assertEqual(data, {((0, 0, 0), (0, 0)): -1.0 + 0.5j})


def _write_chain(d, transfer_lines, coulomb_lines=None):
    """One physical orbital on a chain; transfer/coulombinter given as lines."""
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("  1.0   0.0   0.0\n  0.0   1.0   0.0\n  0.0   0.0   1.0\n")
        f.write("1\n   0.0   0.0   0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("\n".join(transfer_lines) + "\n")
    if coulomb_lines is not None:
        with open(os.path.join(d, "coulombinter.dat"), "w") as f:
            f.write("\n".join(coulomb_lines) + "\n")


_SPIN_CHAIN = [
    MARKED, "1", "2", "1 1",
    " 1 0 0  1 1  0 0  -1.0 0.0", "-1 0 0  1 1  0 0  -1.0 0.0",
    " 1 0 0  1 1  1 1  -1.0 0.0", "-1 0 0  1 1  1 1  -1.0 0.0",
]
_PLAIN_CHAIN = [PLAIN, "1", "2", "1 1",
                " 1 0 0  1 1  -1.0 0.0", "-1 0 0  1 1  -1.0 0.0"]


class TestInputReaderAndSolverRefusal(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="w90spin_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _input(self, with_coulomb=False):
        import hwave.qlmsio.read_input_k as read_input_k
        files = {'path_to_input': self.d, 'Geometry': 'geom.dat',
                 'Transfer': 'transfer.dat'}
        if with_coulomb:
            files['CoulombInter'] = 'coulombinter.dat'
        return read_input_k.QLMSkInput(
            {'path_to_input': self.d, 'interaction': files})

    def test_input_reader_records_format(self):
        _write_chain(self.d, _SPIN_CHAIN)
        ham = self._input().get_param("ham")
        self.assertEqual(ham.source_formats["transfer"], ("w90spin", 1))
        self.assertEqual(ham.source_dimensions["transfer"], 1)
        # folded keys: spin-block (s*1 + 0) = s
        self.assertIn(((1, 0, 0), (1, 1)), ham["Transfer"])
        _write_chain(self.d, _PLAIN_CHAIN)
        ham = self._input().get_param("ham")
        self.assertIsNone(ham.source_formats["transfer"])

    def test_two_body_file_in_spin_format_refused(self):
        _write_chain(self.d, _PLAIN_CHAIN, coulomb_lines=[
            "CoulombInter w90spin/1", "1", "1", "1",
            " 1 0 0  1 1  0 0  0.5 0.0"])
        with self.assertRaisesRegex(ValueError, "accepted only for the Transfer"):
            self._input(with_coulomb=True)

    def _uhfk(self, enable_spin_orbital=False):
        import hwave.solver.uhfk as uhfk_module
        param = {'Ncond': 1, 'T': 0.5, 'CellShape': [4, 1, 1],
                 'SubShape': [1, 1, 1],
                 'IterationMax': 1, 'EPS': 8, 'Mix': 0.5, 'RndSeed': 1}
        info_mode = {'mode': 'UHFk', 'param': param,
                     'enable_spin_orbital': enable_spin_orbital}
        return uhfk_module.UHFk(self._input().get_param("ham"),
                                {"print_level": 0, "print_step": 100}, info_mode)

    def test_uhfk_consumes_spin_resolved_table(self):
        import numpy as np
        _write_chain(self.d, _SPIN_CHAIN)
        solver = self._uhfk()
        self.assertTrue(solver.spin_resolved_transfer)
        solver._make_ham_trans()
        # both spins carry the same -t chain: T(k) = -2 cos k on each spin
        self.assertEqual(solver.ham_trans.shape, (4, 2, 2))
        k = 2 * np.pi * np.arange(4) / 4
        np.testing.assert_allclose(solver.ham_trans[:, 0, 0], -2 * np.cos(k),
                                   atol=1e-12)
        np.testing.assert_allclose(solver.ham_trans[:, 1, 1], -2 * np.cos(k),
                                   atol=1e-12)
        np.testing.assert_allclose(solver.ham_trans[:, 0, 1], 0.0, atol=1e-12)

    def test_uhfk_refuses_clash_with_enable_spin_orbital(self):
        _write_chain(self.d, _SPIN_CHAIN)
        with self.assertRaisesRegex(ValueError, "enable_spin_orbital = true"):
            self._uhfk(enable_spin_orbital=True)

    def test_uhfk_plain_file_unaffected(self):
        _write_chain(self.d, _PLAIN_CHAIN)
        solver = self._uhfk()
        self.assertEqual(solver.norb_orig, 1)
        self.assertFalse(solver.enable_spin_orbital)

    def test_rpa_consumes_spin_resolved_table(self):
        import hwave.solver.rpa as rpa_module
        _write_chain(self.d, _SPIN_CHAIN)
        info_mode = {'mode': 'RPA', 'calc_scheme': 'reduced', 'param': {
            'T': 0.5, 'Ncond': 1, 'CellShape': [4, 1, 1],
            'SubShape': [1, 1, 1], 'Nmat': 16}}
        solver = rpa_module.RPA(self._input().get_param("ham"),
                                {"print_level": 0, "print_step": 100}, info_mode)
        self.assertTrue(solver.ham_info.spin_resolved_transfer)
        self.assertTrue(solver.ham_info.spinful_one_body)
        self.assertFalse(solver.ham_info.enable_spin_orbital)
        # nd x nd in spin-block order, both spins carrying the -t chain
        self.assertEqual(solver.ham_info.ham_trans_q.shape, (4, 2, 2))

    def test_rpa_refuses_clash_with_enable_spin_orbital(self):
        import hwave.solver.rpa as rpa_module
        _write_chain(self.d, _SPIN_CHAIN)
        info_mode = {'mode': 'RPA', 'calc_scheme': 'general',
                     'enable_spin_orbital': True, 'param': {
                         'T': 0.5, 'Ncond': 1, 'CellShape': [4, 1, 1],
                         'SubShape': [1, 1, 1], 'Nmat': 16}}
        with self.assertRaisesRegex(ValueError, "enable_spin_orbital = true"):
            rpa_module.RPA(self._input().get_param("ham"),
                           {"print_level": 0, "print_step": 100}, info_mode)

    def test_flex_hf_bare_transfer_refuses(self):
        from hwave.solver import flex_hf
        _write_chain(self.d, _SPIN_CHAIN)
        with self.assertRaisesRegex(ValueError, "accepted only for the Transfer"):
            flex_hf._bare_transfer_k(os.path.join(self.d, "transfer.dat"),
                                     [2, 1, 1], 1)


if __name__ == "__main__":
    unittest.main()
