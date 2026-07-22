"""Tests for orbitalidx_general_emitter (all-unique-classes fallback).

The emitter is triggered under SOC + SubShape > [1, 1, 1] to bypass
StdFace's class over-grouping. The tests verify that:
  1. The emitted file has the expected NOrbitalIdx count.
  2. The emitted file round-trips through
     ``parse_orbitalidx_general_def`` without errors.
  3. Every upper-triangle (all_i, all_j) pair receives exactly one
     unique class idx (i.e., no class is shared across pairs).
  4. All emitted signs are +1.
  5. Invalid arguments raise ``ValueError``.
"""
from __future__ import annotations

import os
import sys
import tempfile

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.orbitalidx_general_emitter import (
    emit_orbitalidx_all_unique,
)
from tools._uhfk_to_mvmc.orbitalidx_general_reader import (
    detect_orbitalidx_format,
    parse_orbitalidx_general_def,
)


@pytest.mark.parametrize("nsite,expected", [
    (1, 1),   # 1 * (2*1 - 1) = 1
    (2, 6),   # 2 * (2*2 - 1) = 6
    (3, 15),  # 3 * (2*3 - 1) = 15
    (4, 28),  # 4 * (2*4 - 1) = 28
])
def test_emit_orbitalidx_all_unique_class_count(nsite, expected):
    """Emitter returns ``nsite * (2*nsite - 1)`` classes = 2*Ns^2 - Ns."""
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        n_classes = emit_orbitalidx_all_unique(nsite, path)
        assert n_classes == expected
    finally:
        os.unlink(path)


def test_emit_orbitalidx_all_unique_reads_back_via_parser():
    """The emitted file parses through the general reader and reports
    the same Nsite, class count, and every mapping entry has sign +1."""
    nsite = 4
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        n_classes = emit_orbitalidx_all_unique(nsite, path)
        # Format detection agrees the emitted file is 6-column general.
        assert detect_orbitalidx_format(path) == "general"
        info = parse_orbitalidx_general_def(path)
        assert info["nsite"] == nsite
        assert info["n_orbital_idx"] == n_classes
        assert info["complex_type"] == 1
        # Every mapping entry has sign +1 (no antisymmetry-encoded rows).
        signs = [sg for _, sg in info["mapping"].values()]
        assert set(signs) == {1}
    finally:
        os.unlink(path)


def test_emit_orbitalidx_all_unique_covers_full_upper_triangle():
    """Every (all_i, all_j) with all_i < all_j (0 <= all < 2*Ns) receives
    exactly one class idx, and class idxs run 0..n_classes-1 with no
    duplicates or gaps."""
    nsite = 5
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        n_classes = emit_orbitalidx_all_unique(nsite, path)
        info = parse_orbitalidx_general_def(path)
        expected_pairs = {
            (a, b)
            for a in range(2 * nsite)
            for b in range(a + 1, 2 * nsite)
        }
        assert set(info["mapping"].keys()) == expected_pairs
        class_idxs = [idx for idx, _ in info["mapping"].values()]
        assert set(class_idxs) == set(range(n_classes))
        # optimize flag one per class.
        assert len(info["optimize_flags"]) == n_classes
    finally:
        os.unlink(path)


def test_emit_orbitalidx_all_unique_complex_type_0():
    """complex_type=0 emits ComplexType 0 in the header (the parser
    accepts either 0 or 1 at the parse layer; the bridge itself
    additionally rejects ComplexType != 1 in the CLI, but the emitter
    should honor whichever value the caller passed)."""
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        emit_orbitalidx_all_unique(2, path, complex_type=0)
        info = parse_orbitalidx_general_def(path)
        assert info["complex_type"] == 0
    finally:
        os.unlink(path)


def test_emit_orbitalidx_all_unique_rejects_bad_nsite():
    """nsite < 1 raises ValueError."""
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        with pytest.raises(ValueError, match="nsite"):
            emit_orbitalidx_all_unique(0, path)
        with pytest.raises(ValueError, match="nsite"):
            emit_orbitalidx_all_unique(-3, path)
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_emit_orbitalidx_all_unique_rejects_bad_complex_type():
    """complex_type not in {0, 1} raises ValueError."""
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        with pytest.raises(ValueError, match="complex_type"):
            emit_orbitalidx_all_unique(2, path, complex_type=2)
        with pytest.raises(ValueError, match="complex_type"):
            emit_orbitalidx_all_unique(2, path, complex_type=-1)
    finally:
        if os.path.exists(path):
            os.unlink(path)


def test_emit_orbitalidx_all_unique_aggregate_class_consistency():
    """End-to-end: emitting the all-unique file and feeding it to
    ``aggregate_general_orbital_params`` with a randomly constructed
    antisymmetric F matrix trivially passes the class-consistency check
    (each class contains exactly one signed F entry)."""
    import numpy as np
    from tools._uhfk_to_mvmc.general_output_writer import (
        aggregate_general_orbital_params,
    )
    nsite = 3
    with tempfile.NamedTemporaryFile(
        "w", suffix=".def", delete=False, encoding="utf-8"
    ) as tmp:
        path = tmp.name
    try:
        n_classes = emit_orbitalidx_all_unique(nsite, path, complex_type=1)
        info = parse_orbitalidx_general_def(path)
        # Random complex antisymmetric F (deterministic seed).
        rng = np.random.default_rng(7919)
        M = rng.standard_normal((2 * nsite, 2 * nsite)) + \
            1j * rng.standard_normal((2 * nsite, 2 * nsite))
        F = M - M.T
        params = aggregate_general_orbital_params(
            F, info["mapping"], n_classes,
            epsilon_noise=0.0, complex_type=1,
            rng=rng, class_consistency_tol=1e-12,
        )
        assert params.shape == (n_classes,)
        # Each class idx corresponds to exactly one (all_i, all_j)
        # entry with sign=+1, so params[idx] should equal F[all_i, all_j].
        for (all_i, all_j), (idx, sign) in info["mapping"].items():
            assert np.isclose(params[idx], sign * F[all_i, all_j])
    finally:
        os.unlink(path)
