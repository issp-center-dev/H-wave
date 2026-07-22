"""Tests for the ``orbitalidx.def`` parser."""
from __future__ import annotations

import sys, os, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

from tools._uhfk_to_mvmc.orbitalidx_reader import (
    parse_orbitalidx_def,
    OrbitalidxFormatError,
)


def _write_orbitalidx(tmp_path, body, n_orb_idx, complex_type):
    """Helper to write a minimal orbitalidx.def into tmp_path."""
    content = (
        "====================================\n"
        f"NOrbitalIdx {n_orb_idx}\n"
        f"ComplexType {complex_type}\n"
        "====================================\n"
        "====================================\n"
        f"{body}"
    )
    path = os.path.join(tmp_path, "orbitalidx.def")
    with open(path, "w") as fp:
        fp.write(content)
    return path


def test_parse_pbc_no_sign_column():
    """PBC mode: 3-column form (i, j, orbital_idx) + optimize-flag lines."""
    with tempfile.TemporaryDirectory() as tmp:
        body = (
            "0 0 0\n"
            "0 1 1\n"
            "1 0 1\n"
            "1 1 0\n"
            "0 1\n"
            "1 1\n"
        )
        path = _write_orbitalidx(tmp, body, n_orb_idx=2, complex_type=1)
        info = parse_orbitalidx_def(path)

        assert info["n_orbital_idx"] == 2
        assert info["complex_type"] == 1
        assert info["mapping"][(0, 0)] == (0, 1)
        assert info["mapping"][(0, 1)] == (1, 1)
        assert info["mapping"][(1, 0)] == (1, 1)
        assert info["mapping"][(1, 1)] == (0, 1)
        assert info["optimize_flags"] == {0: 1, 1: 1}
        assert info["nsite"] == 2


def test_parse_apbc_with_sign_column():
    """APBC mode: 4-column form (i, j, orbital_idx, sign)."""
    with tempfile.TemporaryDirectory() as tmp:
        body = (
            "0 0 0  1\n"
            "0 1 1  1\n"
            "1 0 1 -1\n"
            "1 1 0  1\n"
            "0 1\n"
            "1 1\n"
        )
        path = _write_orbitalidx(tmp, body, n_orb_idx=2, complex_type=1)
        info = parse_orbitalidx_def(path)

        assert info["has_sign_column"]
        assert info["mapping"][(1, 0)] == (1, -1)
        assert info["mapping"][(0, 0)] == (0, 1)


def test_missing_norbitalidx_raises():
    """A file without NOrbitalIdx header must raise."""
    with tempfile.TemporaryDirectory() as tmp:
        bad = os.path.join(tmp, "orbitalidx.def")
        with open(bad, "w") as fp:
            fp.write("ComplexType 1\n0 0 0\n")
        with pytest.raises(OrbitalidxFormatError, match="NOrbitalIdx"):
            parse_orbitalidx_def(bad)


def test_mixed_3_and_4_column_raises():
    """Mixing 3-column and 4-column body lines is invalid."""
    with tempfile.TemporaryDirectory() as tmp:
        body = (
            "0 0 0\n"
            "0 1 1  1\n"
            "1 0 1\n"
            "1 1 0\n"
            "0 1\n"
            "1 1\n"
        )
        path = _write_orbitalidx(tmp, body, n_orb_idx=2, complex_type=1)
        with pytest.raises(OrbitalidxFormatError, match="mixed"):
            parse_orbitalidx_def(path)
