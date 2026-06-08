#!/usr/bin/env python3
"""Tests for the Wannier90-format reader (src/hwave/qlmsio/wan90.py).

Covers:
- degeneracy (ndegen) weighting: H(k) = sum_R e^{ikR} H(R) / ndegen(R)
- declared-count validation
- duplicate (R, orbital) detection
"""

import os
import tempfile
import unittest

from hwave.qlmsio import wan90


def _write_hr(lines):
    fd, path = tempfile.mkstemp(suffix="_hr.dat")
    os.close(fd)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


class TestReadW90Degeneracy(unittest.TestCase):
    def test_degeneracy_divides_matrix_elements(self):
        # num_wann=1, nrpts=2; degeneracies [1, 2].
        # The second R-point value must be divided by its degeneracy 2.
        lines = [
            "comment line",
            "1",
            "2",
            "1 2",
            "0 0 0 1 1 4.0 0.0",
            "1 0 0 1 1 6.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            data = wan90.read_w90(path)
        finally:
            os.remove(path)

        self.assertAlmostEqual(data[((0, 0, 0), (0, 0))].real, 4.0, places=12)
        # 6.0 / ndegen(=2) = 3.0
        self.assertAlmostEqual(data[((1, 0, 0), (0, 0))].real, 3.0, places=12)

    def test_all_ones_degeneracy_unchanged(self):
        # H-wave's own files always write degeneracy 1: values must be untouched.
        lines = [
            "comment line",
            "1",
            "3",
            "1 1 1",
            "-1 0 0 1 1 -1.0 0.0",
            "0 0 0 1 1 2.0 0.0",
            "1 0 0 1 1 -1.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            data = wan90.read_w90(path)
        finally:
            os.remove(path)

        self.assertAlmostEqual(data[((-1, 0, 0), (0, 0))].real, -1.0, places=12)
        self.assertAlmostEqual(data[((0, 0, 0), (0, 0))].real, 2.0, places=12)
        self.assertAlmostEqual(data[((1, 0, 0), (0, 0))].real, -1.0, places=12)


class TestReadW90Validation(unittest.TestCase):
    def test_duplicate_entry_is_rejected(self):
        # Two lines for the same (R, m, n) must be detected and rejected.
        lines = [
            "comment line",
            "1",
            "1",
            "1",
            "0 0 0 1 1 4.0 0.0",
            "0 0 0 1 1 5.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            with self.assertRaises(SystemExit):
                wan90.read_w90(path)
        finally:
            os.remove(path)

    def test_degeneracy_count_mismatch_is_rejected(self):
        # Declares nr=3 R-grid points but the degeneracy line has only 2 entries.
        lines = [
            "comment line",
            "1",
            "3",
            "1 1",
            "0 0 0 1 1 4.0 0.0",
            "1 0 0 1 1 6.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            with self.assertRaises(SystemExit):
                wan90.read_w90(path)
        finally:
            os.remove(path)

    def test_more_distinct_rpoints_than_declared_is_rejected(self):
        # With non-unit degeneracy the positional mapping is strict: declaring
        # nr=2 but listing three distinct R-vectors is malformed.
        lines = [
            "comment line",
            "1",
            "2",
            "1 2",
            "0 0 0 1 1 4.0 0.0",
            "1 0 0 1 1 6.0 0.0",
            "2 0 0 1 1 8.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            with self.assertRaises(SystemExit):
                wan90.read_w90(path)
        finally:
            os.remove(path)

    def test_all_unit_degeneracy_allows_more_rpoints_than_nr(self):
        # H-wave's loose format may declare nr=1 with a single '1' degeneracy
        # yet list several R-vectors; with all degeneracies 1 this is accepted.
        lines = [
            "comment line",
            "1",
            "1",
            "1",
            "0 0 0 1 1 4.0 0.0",
            "1 0 0 1 1 6.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            data = wan90.read_w90(path)
        finally:
            os.remove(path)
        self.assertAlmostEqual(data[((1, 0, 0), (0, 0))].real, 6.0, places=12)

    def test_nonunit_degeneracy_with_sparse_data_is_rejected(self):
        # Non-unit degeneracy requires a dense listing so the positional
        # R -> ndegen mapping is well-defined. nr=3 with a degeneracy of 2 but
        # only 2 distinct R-vectors present must be rejected rather than guessed.
        lines = [
            "comment line",
            "1",
            "3",
            "1 2 1",
            "0 0 0 1 1 4.0 0.0",
            "1 0 0 1 1 6.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            with self.assertRaises(SystemExit):
                wan90.read_w90(path)
        finally:
            os.remove(path)

    def test_sparse_format_with_fewer_rpoints_is_accepted(self):
        # H-wave's sparse format: nr=9 R-grid points (degeneracy all 1) but only
        # the nonzero hoppings are listed. This must NOT be rejected.
        lines = [
            "comment line",
            "1",
            "9",
            "1 1 1 1 1 1 1 1 1",
            "-1 0 0 1 1 -1.0 0.0",
            "0 -1 0 1 1 -1.0 0.0",
            "0 1 0 1 1 -1.0 0.0",
            "1 0 0 1 1 -1.0 0.0",
        ]
        path = _write_hr(lines)
        try:
            data = wan90.read_w90(path)
        finally:
            os.remove(path)
        self.assertEqual(len(data), 4)
        self.assertAlmostEqual(data[((1, 0, 0), (0, 0))].real, -1.0, places=12)


if __name__ == "__main__":
    unittest.main()
