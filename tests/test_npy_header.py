#!/usr/bin/env python3
"""The shared .npy header reader and the contracts its two callers rely on.

`hwave.solver.npy_header` exists so the Eliashberg frequency probe and the
bond preflight's Green peek answer "what shape is this array?" the same
way, without loading the array and without numpy's private API (numpy 2.4
removed the private helpers both used).

The cases below are the ones that distinguish a correct reader from a
lenient one. Version 3.0 matters most: numpy writes it whenever a dtype
carries field names outside latin-1, publishes no reader for it, and its
2.0 reader ACCEPTS malformed 3.0 headers that `np.load` rejects -- which
would let the probe report a frequency count for a file that cannot be
loaded.
"""

import io
import os
import tempfile
import unittest
import zipfile
import warnings

import numpy as np

from hwave.solver import npy_header


def _npy_bytes(version, header_text, encoding="utf-8"):
    """A synthetic .npy stream carrying `header_text` at `version`."""
    hb = header_text.encode(encoding)
    pad = 64 - ((10 + len(hb) + 1) % 64)
    body = hb + b" " * pad + b"\n"
    return (b"\x93NUMPY" + bytes(version)
            + len(body).to_bytes(4, "little") + body)


def _member(array):
    """`array` written through numpy, opened as its .npy member."""
    d = tempfile.mkdtemp(prefix="npy_header_")
    path = os.path.join(d, "a.npz")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        np.savez(path, k=array)
    return path


class TestReadNpyHeaderShape(unittest.TestCase):
    def test_reads_an_ordinary_version_1_0_header(self):
        path = _member(np.zeros((7, 2), dtype=complex))
        with zipfile.ZipFile(path) as z, z.open("k.npy") as f:
            self.assertEqual(npy_header.read_npy_header_shape(f), (7, 2))

    def test_reads_a_genuine_version_3_0_header(self):
        # Non-latin-1 field names are what makes numpy choose 3.0; the
        # pre-existing 3.0 tests force the version on an ASCII dtype and
        # so cannot tell a utf-8 reader from a latin-1 one.
        dt = np.dtype([("温度", "c16"), ("μ", "f8")])
        path = _member(np.zeros((12, 3), dtype=dt))
        with zipfile.ZipFile(path) as z, z.open("k.npy") as f:
            shape = npy_header.read_npy_header_shape(f)
        self.assertEqual(shape, (12, 3))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(np.load(path)["k"].shape, shape)

    def test_rejects_a_3_0_header_numpy_itself_would_reject(self):
        # numpy's 2.0 reader applies Python-2 compatibility filtering and
        # accepts this, returning (1, 24, 1); np.load raises. Delegating
        # 3.0 to that reader would report a shape for an unloadable file.
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, "
            "'shape': (1L, 24L, 1L)}")
        with self.assertRaises(ValueError):
            npy_header.read_npy_header_shape(io.BytesIO(raw))

    def test_rejects_a_3_0_header_that_is_not_valid_utf8(self):
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, 'shape': (4,)}")
        # Corrupt one header byte into an invalid utf-8 lead byte.
        i = raw.index(b"descr")
        raw = raw[:i] + b"\xff" + raw[i + 1:]
        with self.assertRaises(ValueError):
            npy_header.read_npy_header_shape(io.BytesIO(raw))

    def test_rejects_a_truncated_3_0_header(self):
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, 'shape': (4,)}")
        with self.assertRaises(ValueError):
            npy_header.read_npy_header_shape(io.BytesIO(raw[:20]))

    def test_rejects_a_3_0_shape_that_is_not_integers(self):
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, "
            "'shape': ('a', 2)}")
        with self.assertRaises(ValueError):
            npy_header.read_npy_header_shape(io.BytesIO(raw))

    def test_unknown_future_version_raises_the_documented_exception(self):
        raw = b"\x93NUMPY" + bytes((9, 9)) + b"\x00" * 4
        with self.assertRaises(npy_header.UnsupportedNpyHeaderVersion):
            npy_header.read_npy_header_shape(io.BytesIO(raw))


class TestCallerContracts(unittest.TestCase):
    """Each caller needs a DIFFERENT failure mode; both are pinned here."""

    def test_eliashberg_probe_returns_none_for_an_unreadable_header(self):
        # Its contract is to never raise: the authoritative loader below
        # it must be the one to produce the diagnostic.
        from hwave.solver import eliashberg_dynamic as ed

        d = tempfile.mkdtemp(prefix="npy_header_bad_")
        path = os.path.join(d, "chiq_s.npz")
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, "
            "'shape': (1L, 24L, 1L)}")
        with zipfile.ZipFile(path, "w") as z:
            z.writestr("chiq_s.npy", raw + b"\x00" * 64)
        self.assertIsNone(
            ed._npz_freq_size(path, ("chiq_s", "chiq"), axis=0))

    def test_sc_peek_raises_its_own_exception_for_an_unknown_version(self):
        # Its contract is the opposite: an unknown version must surface,
        # so the caller's except branch can turn it into a clear error
        # rather than silently preflighting on a guess.
        import hwave.sc as sc

        raw = b"\x93NUMPY" + bytes((9, 9)) + b"\x00" * 4
        with self.assertRaises(sc._UnsupportedNpyHeaderVersion):
            sc._read_npy_header_shape(io.BytesIO(raw))


if __name__ == "__main__":
    unittest.main()
