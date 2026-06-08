#!/usr/bin/env python3
"""Robustness tests for the real-space input reader (read_input.py)."""

import os
import tempfile
import unittest

from hwave.qlmsio.read_input import QLMSInput


def _write(lines):
    fd, path = tempfile.mkstemp(suffix=".def")
    os.close(fd)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _stub(file_names):
    inp = QLMSInput.__new__(QLMSInput)
    inp.file_names = file_names
    return inp


class TestReadGreen(unittest.TestCase):
    def test_single_row_onebodyg_is_2d(self):
        # A one-entry OneBodyG file must still yield a 2D (1, 4) array so that
        # downstream iteration over rows works.
        path = _write([
            "comment", "ngreen 1", "h2", "h3", "h4",
            "0 0 1 0",
        ])
        inp = _stub({"onebodyg": path})
        try:
            data = inp._read_green("onebodyg")
        finally:
            os.remove(path)
        self.assertEqual(data.ndim, 2)
        self.assertEqual(data.shape, (1, 4))


class TestLoadTextRobustness(unittest.TestCase):
    def test_trailing_blank_line_is_ignored(self):
        path = _write([
            "comment", "ncoulomb 2", "h2", "h3", "h4",
            "0 0 4.0",
            "1 1 4.0",
            "",          # trailing blank line must not crash
        ])
        inp = _stub({})
        try:
            data = inp._load_text(path, "real")
        finally:
            os.remove(path)
        self.assertEqual(len(data), 2)

    def test_inconsistent_column_count_is_rejected(self):
        path = _write([
            "comment", "ntrans 2", "h2", "h3", "h4",
            "0 0 0 0 1.0 0.0",
            "1 0 0 1.0 0.0",   # one index column short
        ])
        inp = _stub({})
        try:
            with self.assertRaises(SystemExit):
                inp._load_text(path, "complex")
        finally:
            os.remove(path)

    def test_missing_file_is_handled(self):
        inp = _stub({})
        with self.assertRaises(SystemExit):
            inp._load_text("/nonexistent/path/to/file.def", "real")


if __name__ == "__main__":
    unittest.main()
