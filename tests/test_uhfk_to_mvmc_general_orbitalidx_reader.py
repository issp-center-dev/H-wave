"""Tests for orbitalidx_general_reader (InOrbitalGeneral 6-column parser)."""
from __future__ import annotations

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.orbitalidx_general_reader import (
    OrbitalidxFormatError,
    detect_orbitalidx_format,
    parse_orbitalidx_general_def,
)


def _write(text):
    tmp = tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    )
    tmp.write(text)
    tmp.close()
    return tmp.name


def _general_body(nsite):
    """Return a minimal valid orbitalidx_general.def body for the given
    Nsite (Ns=nsite). Every unordered spin-orbital pair (all_i < all_j)
    is its own class with sign=+1."""
    total = 2 * nsite * nsite - nsite
    lines = [
        "======================",
        f"NOrbitalIdx  {total}",
        "ComplexType 1",
        "======================",
        "== i_spn_j_spn_OrbitalIdx ==",
        "======================",
    ]
    idx = 0
    for all_i in range(2 * nsite):
        for all_j in range(all_i + 1, 2 * nsite):
            i, spn_i = all_i % nsite, all_i // nsite
            j, spn_j = all_j % nsite, all_j // nsite
            lines.append(f"{i} {spn_i} {j} {spn_j} {idx} 1")
            idx += 1
    for k in range(total):
        lines.append(f"{k} 1")
    return "\n".join(lines) + "\n"


def test_parse_orbitalidx_general_def_returns_expected_metadata():
    """Minimal 2-site general def parses NOrbitalIdx, ComplexType, mapping."""
    nsite = 2
    path = _write(_general_body(nsite))
    info = parse_orbitalidx_general_def(path)
    total = 2 * nsite * nsite - nsite  # = 6
    assert info["n_orbital_idx"] == total
    assert info["complex_type"] == 1
    assert info["nsite"] == nsite
    assert info["format"] == "general"
    assert info["has_sign_column"] is True
    # Sample mapping entry: (all_i=0, all_j=1) -> (0, +1)
    assert info["mapping"][(0, 1)] == (0, 1)


def test_detect_orbitalidx_format_general_6col():
    """6-column body → 'general'."""
    path = _write(_general_body(nsite=2))
    assert detect_orbitalidx_format(path) == "general"


def test_detect_orbitalidx_format_antiparallel_3col():
    """3-column PBC body maps to 'antiparallel'."""
    body = (
        "======================\n"
        "NOrbitalIdx  4\n"
        "ComplexType 1\n"
        "======================\n"
        "== i_j_OrbitalIdx ==\n"
        "======================\n"
        "0 0 0\n0 1 1\n1 0 2\n1 1 3\n"
        "0 1\n1 1\n2 1\n3 1\n"
    )
    path = _write(body)
    assert detect_orbitalidx_format(path) == "antiparallel"


def test_detect_orbitalidx_format_antiparallel_4col():
    """4-column APBC body maps to 'antiparallel'."""
    body = (
        "======================\n"
        "NOrbitalIdx  4\n"
        "ComplexType 1\n"
        "======================\n"
        "== i_j_OrbitalIdx ==\n"
        "======================\n"
        "0 0 0 1\n0 1 1 1\n1 0 2 1\n1 1 3 1\n"
        "0 1\n1 1\n2 1\n3 1\n"
    )
    path = _write(body)
    assert detect_orbitalidx_format(path) == "antiparallel"


def test_detect_orbitalidx_format_invalid_arity_raises():
    """5-column body (non-existent format) → OrbitalidxFormatError."""
    body = (
        "======================\n"
        "NOrbitalIdx  1\n"
        "ComplexType 1\n"
        "======================\n"
        "0 0 1 0 0\n"
        "0 1\n"
    )
    path = _write(body)
    with pytest.raises(OrbitalidxFormatError, match="arity"):
        detect_orbitalidx_format(path)


def test_parse_orbitalidx_general_rejects_upper_triangle_violation():
    """all_i >= all_j on any row → OrbitalidxFormatError."""
    body = (
        "======================\n"
        "NOrbitalIdx  6\n"
        "ComplexType 1\n"
        "======================\n"
        "== i_spn_j_spn_OrbitalIdx ==\n"
        "======================\n"
        # Second row has all_i (=1+0*2=1) == all_j (=1+0*2=1) → invalid
        "0 0 1 0 0 1\n"
        "1 0 1 0 1 1\n"
        "0 0 0 1 2 1\n"
        "1 0 0 1 3 1\n"
        "0 1 1 1 4 1\n"
        "0 1 0 1 5 1\n"
        "0 1\n1 1\n2 1\n3 1\n4 1\n5 1\n"
    )
    path = _write(body)
    with pytest.raises(OrbitalidxFormatError, match="all_i"):
        parse_orbitalidx_general_def(path)


def test_parse_orbitalidx_general_rejects_invalid_sign():
    """sign not in {-1, +1} → OrbitalidxFormatError."""
    # Corrupt the first mapping row's sign column (` 0 1` → ` 0 2`) so we
    # actually hit the mapping-row sign check rather than the ComplexType
    # header (which also ends in " 1\n").
    body = _general_body(nsite=2).replace(" 0 1\n", " 0 2\n", 1)
    path = _write(body)
    with pytest.raises(OrbitalidxFormatError, match="sign"):
        parse_orbitalidx_general_def(path)


def test_parse_orbitalidx_general_rejects_row_count_mismatch():
    """Row count != 2*Ns^2 - Ns → OrbitalidxFormatError."""
    body = _general_body(nsite=2)
    # Drop the first mapping line (splitlines index 6 for nsite=2: after
    # the 6-line header block). This leaves 5 mapping rows against the
    # expected 2*Ns^2 - Ns = 6.
    lines = body.splitlines(keepends=True)
    body = "".join(lines[:6] + lines[7:])
    path = _write(body)
    with pytest.raises(OrbitalidxFormatError, match="rows"):
        parse_orbitalidx_general_def(path)


def test_parse_orbitalidx_general_rejects_header_only_file():
    """Header-only file with no mapping rows raises OrbitalidxFormatError
    (regression guard; earlier revision fell through to a raw ValueError
    on `int('NOrbitalIdx')` when data_start was never advanced)."""
    body = (
        "======================\n"
        "NOrbitalIdx  0\n"
        "ComplexType 1\n"
        "======================\n"
    )
    path = _write(body)
    with pytest.raises(OrbitalidxFormatError, match="no mapping rows"):
        parse_orbitalidx_general_def(path)
