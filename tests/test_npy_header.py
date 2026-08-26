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
        # Assert numpy really emitted 3.0: if its version-selection rule
        # ever changes, this test would otherwise slip onto the delegated
        # 1.0/2.0 path and keep passing while covering nothing.
        with zipfile.ZipFile(path) as z, z.open("k.npy") as f:
            self.assertEqual(tuple(np.lib.format.read_magic(f)), (3, 0))
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

    def test_rejects_a_truncated_length_field(self):
        # A 3.0 preamble is 12 bytes (6 magic + 2 version + 4 length), so
        # raw[:12] keeps the length field intact and would only re-test
        # the empty-body path the previous test covers. Cut inside the
        # length field itself.
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, 'shape': (4,)}")
        with self.assertRaises(ValueError) as cm:
            npy_header.read_npy_header_shape(io.BytesIO(raw[:11]))
        self.assertIn("length field", str(cm.exception))

    def test_refuses_a_huge_declared_length_without_reading_the_body(self):
        """The length is four attacker-supplied bytes, so it must be
        checked BEFORE the body read -- a probe exists precisely to avoid
        reading data, and must not be talked into allocating 4 GiB."""

        class _Recorder(io.BytesIO):
            def __init__(self, data):
                super().__init__(data)
                self.requested = []

            def read(self, n=-1):
                self.requested.append(n)
                return super().read(n)

        raw = (b"\x93NUMPY" + bytes((3, 0))
               + (0xFFFFFFFF).to_bytes(4, "little"))
        fh = _Recorder(raw)
        with self.assertRaises(ValueError):
            npy_header.read_npy_header_shape(fh)
        self.assertNotIn(0xFFFFFFFF, fh.requested)
        self.assertTrue(all(n <= npy_header._MAX_HEADER_LEN
                            for n in fh.requested if isinstance(n, int)
                            and n > 0))

    #: Malformed headers numpy rejects in its HEADER validation, so the
    #: two parsers can be required to agree on them across versions.
    _REJECTED_BY_NUMPYS_HEADER_CHECK = {
        "not a dict": "[1, 2, 3]",
        "missing descr": "{'fortran_order': False, 'shape': (4,)}",
        "missing fortran_order": "{'descr': '<c16', 'shape': (4,)}",
        "extra key": ("{'descr': '<c16', 'fortran_order': False, "
                      "'shape': (4,), 'extra': 1}"),
        "fortran_order not bool": ("{'descr': '<c16', 'fortran_order': 0, "
                                   "'shape': (4,)}"),
        "invalid descr": ("{'descr': 'not-a-dtype', "
                          "'fortran_order': False, 'shape': (4,)}"),
    }

    #: Headers this parser rejects MORE strictly than numpy's header
    #: check does. numpy tests dimensions with isinstance(x, int) only,
    #: which admits a negative value and -- since bool subclasses int --
    #: True; np.load then fails later, while READING THE DATA, in a way
    #: that depends on how much payload follows and on the numpy
    #: version. So no agreement with np.load is asserted for these: an
    #: earlier version of this test did assert it and passed locally
    #: while failing on every CI runner.
    _REJECTED_MORE_STRICTLY_THAN_NUMPY = {
        "negative dimension": ("{'descr': '<c16', 'fortran_order': False, "
                               "'shape': (-4,)}"),
        "bool dimension": ("{'descr': '<c16', 'fortran_order': False, "
                           "'shape': (True,)}"),
    }

    def test_rejects_headers_np_load_rejects(self):
        """Refusing what the loader refuses is the point of this parser:
        reporting a shape for a file that cannot be loaded would turn a
        clear load error into a wrong frequency count or an unrelated
        preflight failure."""
        for name, header_text in self._REJECTED_BY_NUMPYS_HEADER_CHECK.items():
            with self.subTest(case=name):
                raw = _npy_bytes((3, 0), header_text)
                with self.assertRaises(ValueError):
                    npy_header.read_npy_header_shape(io.BytesIO(raw))
                d = tempfile.mkdtemp(prefix="npy_header_bad_")
                path = os.path.join(d, "a.npy")
                with open(path, "wb") as f:
                    f.write(raw + b"\x00" * 128)
                with self.assertRaises((ValueError, TypeError)):
                    np.load(path)

    def test_rejects_dimensions_that_are_not_real_array_sizes(self):
        """Deliberately stricter than numpy's header check.

        Both callers use the returned shape as a SIZE -- a frequency
        count, a memory budget -- so a negative or boolean dimension
        would be arithmetic on a value that is not an array size.
        Refusing it early is safe: the caller falls back to the
        authoritative loader, which raises the real error.
        """
        for name, header_text in self._REJECTED_MORE_STRICTLY_THAN_NUMPY.items():
            with self.subTest(case=name):
                raw = _npy_bytes((3, 0), header_text)
                with self.assertRaises(ValueError):
                    npy_header.read_npy_header_shape(io.BytesIO(raw))

    def test_accepts_a_multibyte_header_numpy_accepts(self):
        """Version 3.0 exists FOR names outside latin-1, so the size
        limit must bound the DECODED length as numpy's does. A header of
        a few thousand CJK characters exceeds 10,000 BYTES while staying
        well under 10,000 characters; np.load reads it, so this must
        too -- a byte-based limit would reject exactly the files this
        version was introduced to carry."""
        dt = np.dtype([("温" * 4000, "c16")])
        path = _member(np.zeros((3,), dtype=dt))
        with zipfile.ZipFile(path) as z, z.open("k.npy") as f:
            f.read(8)                      # past magic + version
            byte_len = int.from_bytes(f.read(4), "little")
            char_len = len(f.read(byte_len).decode("utf-8"))
        self.assertGreater(byte_len, npy_header._MAX_HEADER_CHARS)
        self.assertLessEqual(char_len, npy_header._MAX_HEADER_CHARS)
        with zipfile.ZipFile(path) as z, z.open("k.npy") as f:
            self.assertEqual(npy_header.read_npy_header_shape(f), (3,))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(np.load(path)["k"].shape, (3,))

    def test_accepts_a_four_byte_per_character_header_near_the_byte_bound(self):
        """The pre-read byte bound is 4x the character limit because utf-8
        spends at most four bytes per code point. Astral-plane names hit
        that worst case, so a header near 40,000 bytes but under 10,000
        characters must still be read -- a bound trimmed to, say, 30,000
        would reject a file np.load loads, and the ~12,000-byte CJK case
        above would not notice."""
        name = "\U0001F300" * 9000            # 4 bytes each
        text = ("{'descr': [('" + name + "', '<c16')], "
                "'fortran_order': False, 'shape': (2,)}")
        hb = text.encode("utf-8")
        self.assertLessEqual(len(text), npy_header._MAX_HEADER_CHARS)
        self.assertGreater(len(hb), 3 * npy_header._MAX_HEADER_CHARS)
        raw = (b"\x93NUMPY" + bytes((3, 0))
               + len(hb).to_bytes(4, "little") + hb)
        self.assertEqual(
            npy_header.read_npy_header_shape(io.BytesIO(raw)), (2,))
        d = tempfile.mkdtemp(prefix="npy_header_astral_")
        path = os.path.join(d, "a.npy")
        with open(path, "wb") as f:
            f.write(raw + b"\x00" * (16 * 2))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            self.assertEqual(np.load(path).shape, (2,))

    def test_decoded_length_limit_is_the_boundary(self):
        """At the limit the header is read; one character over, refused
        -- the same boundary numpy draws."""
        for extra, should_raise in ((0, False), (1, True)):
            with self.subTest(over_by=extra):
                base = "{'descr': '<c16', 'fortran_order': False, "
                tail = "'shape': (4,)}"
                pad = npy_header._MAX_HEADER_CHARS - len(base) - len(tail)
                text = base + " " * (pad + extra) + tail
                self.assertEqual(
                    len(text), npy_header._MAX_HEADER_CHARS + extra)
                # Build the stream directly: _npy_bytes pads to a 64-byte
                # boundary, which would push the decoded length past the
                # value under test.
                hb = text.encode("utf-8")
                raw = (b"\x93NUMPY" + bytes((3, 0))
                       + len(hb).to_bytes(4, "little") + hb)
                fh = io.BytesIO(raw)
                if should_raise:
                    with self.assertRaises(ValueError):
                        npy_header.read_npy_header_shape(fh)
                else:
                    self.assertEqual(
                        npy_header.read_npy_header_shape(fh), (4,))
                # Lock the boundary to numpy's, not merely to itself: a
                # limit that drifted off by one would still look
                # self-consistent here without this.
                d = tempfile.mkdtemp(prefix="npy_header_bound_")
                path = os.path.join(d, "a.npy")
                with open(path, "wb") as f:
                    f.write(raw + b"\x00" * (16 * 4))
                if should_raise:
                    with self.assertRaises(ValueError):
                        np.load(path)
                else:
                    self.assertEqual(np.load(path).shape, (4,))

    def test_rejects_mixed_key_types_without_leaking_typeerror(self):
        """sc's probe catches only OSError/ValueError, so a header whose
        keys cannot be ordered against each other must still come back
        as ValueError rather than escaping as TypeError."""
        raw = _npy_bytes(
            (3, 0),
            "{'descr': '<c16', 'fortran_order': False, "
            "'shape': (4,), 0: 'x'}")
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
