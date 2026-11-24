#!/usr/bin/env python3
"""chi0q index-convention validation on the read side.

RPA writes ``index_convention="spin_block"`` into chi0q.npz / chiq.npz so that
spin-orbital (SO) consumers do not silently mix it with UHFk's interleaved
``2*orb+spin`` ordering. The read paths must enforce this: a stale chi0q file
produced before the SO convention fix (commit 9dd9a21) carries no
``index_convention`` key, and feeding it into an SO calculation would silently
use the wrong ordering. These tests pin the fail-fast behaviour for both
consumers:

  * ``hwave.solver.rpa._read_chi0q`` (the ``chi0q_init`` resume path), via the
    shared helper ``validate_chi0q_index_convention``, and
  * ``hwave.sc._load_chi0q`` (the Eliashberg downstream consumer).
"""

import os
import tempfile
import unittest

import numpy as np

from hwave.solver.rpa import validate_chi0q_index_convention
from hwave.sc import _load_chi0q


def _write_npz(path, with_convention, convention="spin_block",
               spin_orbital=None):
    kwargs = dict(chi0q=np.zeros((1, 1, 2, 2), dtype=complex))
    if with_convention:
        kwargs["index_convention"] = convention
    if spin_orbital is not None:
        kwargs["spin_orbital"] = spin_orbital
    np.savez(path, **kwargs)


class TestValidateIndexConvention(unittest.TestCase):
    """Shared helper used by the RPA chi0q_init resume path."""

    def _load(self, with_convention, convention="spin_block"):
        tmp = tempfile.mkdtemp()
        fn = os.path.join(tmp, "chi0q.npz")
        _write_npz(fn, with_convention, convention)
        return np.load(fn), fn

    def test_non_so_skips_validation(self):
        # Without SO, the convention is irrelevant: a file lacking the key loads.
        data, fn = self._load(with_convention=False)
        validate_chi0q_index_convention(data, enable_spin_orbital=False,
                                        file_name=fn)  # must not raise

    def test_so_missing_convention_rejected(self):
        data, fn = self._load(with_convention=False)
        with self.assertRaises(ValueError):
            validate_chi0q_index_convention(data, enable_spin_orbital=True,
                                            file_name=fn)

    def test_so_wrong_convention_rejected(self):
        data, fn = self._load(with_convention=True, convention="interleaved")
        with self.assertRaises(ValueError):
            validate_chi0q_index_convention(data, enable_spin_orbital=True,
                                            file_name=fn)

    def test_so_spin_block_accepted(self):
        data, fn = self._load(with_convention=True, convention="spin_block")
        validate_chi0q_index_convention(data, enable_spin_orbital=True,
                                        file_name=fn)  # must not raise


class TestLoadChi0qSC(unittest.TestCase):
    """hwave.sc._load_chi0q (Eliashberg downstream consumer)."""

    def _input_dict(self, tmp, enable_spin_orbital):
        return {
            "file": {"output": {"path_to_output": tmp, "chi0q": "chi0q"}},
            "mode": {"enable_spin_orbital": enable_spin_orbital},
        }

    def test_so_missing_convention_rejected(self):
        tmp = tempfile.mkdtemp()
        _write_npz(os.path.join(tmp, "chi0q.npz"), with_convention=False)
        with self.assertRaises(ValueError):
            _load_chi0q(self._input_dict(tmp, enable_spin_orbital=True))

    def test_so_spin_block_accepted(self):
        tmp = tempfile.mkdtemp()
        _write_npz(os.path.join(tmp, "chi0q.npz"), with_convention=True,
                   spin_orbital=True)
        chi0q, _ = _load_chi0q(self._input_dict(tmp, enable_spin_orbital=True))
        self.assertEqual(chi0q.shape, (1, 1, 2, 2))

    def test_so_needs_spin_orbital_marker(self):
        # every RPA run writes index_convention = "spin_block", so the tag
        # alone does not identify a spin-orbital chi0q
        for marker, pattern in ((None, "no spin_orbital marker"),
                                (False, "spin_orbital = False"),
                                ("true", "malformed spin_orbital"),
                                (np.array([True]), "malformed spin_orbital")):
            with self.subTest(marker=marker):
                tmp = tempfile.mkdtemp()
                _write_npz(os.path.join(tmp, "chi0q.npz"), with_convention=True,
                           spin_orbital=marker)
                with self.assertRaisesRegex(ValueError, pattern):
                    _load_chi0q(self._input_dict(tmp, enable_spin_orbital=True))

    def test_non_so_ignores_spin_orbital_marker(self):
        tmp = tempfile.mkdtemp()
        _write_npz(os.path.join(tmp, "chi0q.npz"), with_convention=True,
                   spin_orbital=False)
        chi0q, _ = _load_chi0q(self._input_dict(tmp, enable_spin_orbital=False))
        self.assertEqual(chi0q.shape, (1, 1, 2, 2))

    def test_non_so_missing_convention_accepted(self):
        tmp = tempfile.mkdtemp()
        _write_npz(os.path.join(tmp, "chi0q.npz"), with_convention=False)
        chi0q, _ = _load_chi0q(self._input_dict(tmp, enable_spin_orbital=False))
        self.assertEqual(chi0q.shape, (1, 1, 2, 2))


if __name__ == "__main__":
    unittest.main()
