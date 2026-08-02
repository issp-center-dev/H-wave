"""Pinning tests for the --debug-writer flag in tools/uhfk_to_mvmc.py.

The G0-writer-check producer contract requires --debug-writer to dump
F_pre_noise.npz + F_post_aggregate.npz. Under --epsilon-noise 0 the two
dumps are bit-identical (identity assertion in G0). Under the shipping
1e-8 rank-lift, they differ by the noise term.
"""
from __future__ import annotations

import argparse
import os

import numpy as np
import pytest

from tools.uhfk_to_mvmc import _debug_dump_writer_frames


def _make_args(tmp_path, debug_writer: bool):
    return argparse.Namespace(
        output=str(tmp_path / "zqp_orbital_uhfk.dat"),
        debug_writer=debug_writer,
    )


def _make_mapping_and_params(n_pairs: int = 3):
    """Construct a small synthetic (mapping, params) pair.

    Returns
    -------
    mapping : dict[(all_i, all_j)] -> (idx, sign)
    params  : (n_pairs,) complex128
    F_pre   : (2*Nsite, 2*Nsite) complex128 antisymmetric matrix that the
        mapping+params reconstruct exactly (before noise).
    """
    Nsite = 3
    dim = 2 * Nsite
    # Three distinct upper-triangle (all_i, all_j) pairs.
    coords = [(0, 1), (0, 4), (2, 3)]
    params = np.array(
        [0.1 + 0.2j, -0.5 + 0.3j, 0.7 - 0.4j], dtype=np.complex128
    )
    mapping = {}
    F_pre = np.zeros((dim, dim), dtype=np.complex128)
    for idx, ((all_i, all_j), sign) in enumerate(zip(coords, [1, -1, 1])):
        mapping[(all_i, all_j)] = (idx, sign)
        F_pre[all_i, all_j] = sign * params[idx]
        F_pre[all_j, all_i] = -F_pre[all_i, all_j]
    return mapping, params, F_pre


def test_debug_dump_writer_frames_noop_without_flag(tmp_path):
    """Without --debug-writer set, the helper writes nothing."""
    mapping, params, F_pre = _make_mapping_and_params()
    args = _make_args(tmp_path, debug_writer=False)
    _debug_dump_writer_frames(args, F_pre, params, mapping)
    assert not (tmp_path / "F_pre_noise.npz").exists()
    assert not (tmp_path / "F_post_aggregate.npz").exists()


def test_debug_dump_writer_frames_writes_both_files(tmp_path):
    """--debug-writer set: both F dumps land next to --output."""
    mapping, params, F_pre = _make_mapping_and_params()
    args = _make_args(tmp_path, debug_writer=True)
    _debug_dump_writer_frames(args, F_pre, params, mapping)
    assert (tmp_path / "F_pre_noise.npz").exists()
    assert (tmp_path / "F_post_aggregate.npz").exists()


def test_debug_dump_frames_bit_identical_under_zero_noise(tmp_path):
    """G0 identity: F_pre_noise == F_post_aggregate exactly
    when the mapping+params round-trip matches F_pre (i.e. the noise
    branch of aggregate_general_orbital_params was disabled)."""
    mapping, params, F_pre = _make_mapping_and_params()
    args = _make_args(tmp_path, debug_writer=True)
    _debug_dump_writer_frames(args, F_pre, params, mapping)
    F_pre_saved = np.load(tmp_path / "F_pre_noise.npz")["F"]
    F_post_saved = np.load(tmp_path / "F_post_aggregate.npz")["F"]
    assert np.array_equal(F_pre_saved, F_pre)
    # Under zero-noise the round-trip is EXACT.
    assert np.array_equal(F_pre_saved, F_post_saved), (
        "F_pre_noise must equal F_post_aggregate under --epsilon-noise 0"
    )


def test_debug_dump_frames_post_aggregate_carries_noise(tmp_path):
    """Under shipping rank-lift noise, the round-trip introduces the
    noise term on the reconstructed F_post_aggregate, while F_pre_noise
    stays clean."""
    mapping, params, F_pre = _make_mapping_and_params()
    # Simulate the noise-added params by perturbing them.
    noise_amp = 1e-6
    rng = np.random.default_rng(42)
    params_noisy = params + (
        rng.uniform(-noise_amp, noise_amp, size=params.shape)
        + 1j * rng.uniform(-noise_amp, noise_amp, size=params.shape)
    )
    args = _make_args(tmp_path, debug_writer=True)
    _debug_dump_writer_frames(args, F_pre, params_noisy, mapping)
    F_pre_saved = np.load(tmp_path / "F_pre_noise.npz")["F"]
    F_post_saved = np.load(tmp_path / "F_post_aggregate.npz")["F"]
    # F_pre_noise unchanged.
    assert np.array_equal(F_pre_saved, F_pre)
    # F_post_aggregate carries the noise; differs from F_pre by the
    # noise amplitude.
    diff = np.max(np.abs(F_post_saved - F_pre_saved))
    assert diff > noise_amp * 0.1, (
        f"expected post-aggregate noise scale ~{noise_amp}; got diff={diff}"
    )
    assert diff < noise_amp * 10.0, (
        f"unexpectedly large drift; got diff={diff}"
    )


def test_debug_dump_frames_post_aggregate_is_antisymmetric(tmp_path):
    """The reconstructed F_post_aggregate must remain antisymmetric —
    the antisymmetric fill in _debug_dump_writer_frames MUST use the
    same all_i / all_j / sign / params round-trip that parse_emitted_F
    applies to the emitted zqp file."""
    mapping, params, F_pre = _make_mapping_and_params()
    args = _make_args(tmp_path, debug_writer=True)
    _debug_dump_writer_frames(args, F_pre, params, mapping)
    F_post_saved = np.load(tmp_path / "F_post_aggregate.npz")["F"]
    assert np.allclose(F_post_saved, -F_post_saved.T, atol=1e-15)


def test_debug_dump_frames_ignores_negative_idx(tmp_path):
    """Mapping entries with idx < 0 (structurally-forced-zero classes)
    MUST be skipped; F_post_aggregate stays zero at those slots."""
    Nsite = 3
    dim = 2 * Nsite
    mapping = {
        (0, 1): (0, 1),
        (0, 4): (1, -1),
        (2, 5): (-1, 1),  # idx = -1 -> skip; must NOT read params[-1]
    }
    params = np.array([0.3 + 0.1j, -0.2 + 0.0j], dtype=np.complex128)
    F_pre = np.zeros((dim, dim), dtype=np.complex128)
    F_pre[0, 1] = 0.3 + 0.1j
    F_pre[1, 0] = -F_pre[0, 1]
    F_pre[0, 4] = -1 * (-0.2 + 0.0j)
    F_pre[4, 0] = -F_pre[0, 4]
    args = _make_args(tmp_path, debug_writer=True)
    _debug_dump_writer_frames(args, F_pre, params, mapping)
    F_post_saved = np.load(tmp_path / "F_post_aggregate.npz")["F"]
    # (2, 5) slot MUST stay zero.
    assert F_post_saved[2, 5] == 0.0
    assert F_post_saved[5, 2] == 0.0
    # (0, 1) and (0, 4) round-trip correctly.
    assert F_post_saved[0, 1] == pytest.approx(0.3 + 0.1j, abs=1e-15)
    assert F_post_saved[0, 4] == pytest.approx(0.2 + 0.0j, abs=1e-15)
