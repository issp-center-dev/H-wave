"""Tests for output_writer (zqp_orbital_uhfk.dat text format)."""
from __future__ import annotations

import sys, os, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

from tools._uhfk_to_mvmc.output_writer import (
    aggregate_orbital_params,
    write_zqp_orbital,
)


def test_aggregate_orbital_params_translation_symmetric():
    """2-site translation-symmetric: f_{ij} = f(j-i) ⇒ 2 unique params."""
    F = np.array([
        [0.5 + 0.0j, 0.3 + 0.2j],
        [0.3 + 0.2j, 0.5 + 0.0j],
    ], dtype=np.complex128)
    mapping = {
        (0, 0): (0, 1),
        (0, 1): (1, 1),
        (1, 0): (1, 1),
        (1, 1): (0, 1),
    }
    # Disable rank-lift noise (default 1e-6) so the asserted exact
    # aggregation values are not perturbed; the noise path itself is
    # exercised by test_uhfk_to_mvmc_pfaffian_density.py's APBC case.
    params = aggregate_orbital_params(
        F, mapping, n_orbital_idx=2, epsilon_noise=0,
    )

    assert params[0] == pytest.approx(0.5 + 0.0j)
    assert params[1] == pytest.approx(0.3 + 0.2j)


def test_write_zqp_orbital_format():
    """Output file must put NOrbitalIdx on line 2 (mVMC fscanf reads
    sscanf("%s %d") from line 2). Body lines = "idx real imag"."""
    with tempfile.TemporaryDirectory() as tmp:
        params = np.array([1.0 + 0j, -2.5 - 0.5j], dtype=np.complex128)
        out_path = os.path.join(tmp, "zqp_orbital_uhfk.dat")
        write_zqp_orbital(out_path, params)

        with open(out_path) as fp:
            lines = [ln.rstrip("\n") for ln in fp]

        # 5 header lines + 2 body lines = 7 total
        assert len(lines) == 7
        assert lines[0].startswith("==")
        # Line 2 (0-indexed 1): mVMC sscanf("%s %d") reads this
        toks_l2 = lines[1].split()
        assert toks_l2[0] == "NOrbitalIdx"
        assert int(toks_l2[1]) == 2
        # Lines 3, 4, 5 (0-indexed 2, 3, 4): banners / content markers
        for header_idx in (2, 3, 4):
            assert lines[header_idx].startswith("==") or "OrbitalIdx" in lines[header_idx], (
                f"line {header_idx} expected banner or content marker, got {lines[header_idx]!r}"
            )

        # Body
        toks0 = lines[5].split()
        assert int(toks0[0]) == 0
        assert float(toks0[1]) == pytest.approx(1.0)
        assert float(toks0[2]) == pytest.approx(0.0)
        toks1 = lines[6].split()
        assert int(toks1[0]) == 1
        assert float(toks1[1]) == pytest.approx(-2.5)
        assert float(toks1[2]) == pytest.approx(-0.5)


def test_aggregate_rank_lift_noise_default_is_complex_bounded():
    """Default ``epsilon_noise=1e-6`` perturbs each param by less than
    that amplitude in BOTH real and imag parts when the orbitalidx
    ``complex_type=1`` is used (the bridge default for APBC)."""
    F = np.array([
        [0.5 + 0.0j, 0.3 + 0.2j],
        [0.3 + 0.2j, 0.5 + 0.0j],
    ], dtype=np.complex128)
    mapping = {
        (0, 0): (0, 1),
        (0, 1): (1, 1),
        (1, 0): (1, 1),
        (1, 1): (0, 1),
    }
    params_clean = aggregate_orbital_params(
        F, mapping, n_orbital_idx=2, epsilon_noise=0,
    )
    params_noisy = aggregate_orbital_params(
        F, mapping, n_orbital_idx=2,
        epsilon_noise=1.0e-6, complex_type=1,
    )
    # Each component perturbed by at most 0.5 * 1e-6 (uniform on [-0.5, 0.5])
    deltas = params_noisy - params_clean
    assert np.all(np.abs(deltas.real) <= 1.0e-6), deltas
    assert np.all(np.abs(deltas.imag) <= 1.0e-6), deltas
    # Non-zero: the noise actually fired
    assert np.any(deltas.real != 0.0)
    assert np.any(deltas.imag != 0.0)


def test_aggregate_rank_lift_noise_real_complex_type_zero():
    """When ``complex_type=0`` (real variational parameters), noise must
    be added to the real part only; imag stays zero so the output
    remains a real Fij."""
    F = np.array([
        [0.5 + 0.0j, 0.3 + 0.0j],
        [0.3 + 0.0j, 0.5 + 0.0j],
    ], dtype=np.complex128)
    mapping = {
        (0, 0): (0, 1),
        (0, 1): (1, 1),
        (1, 0): (1, 1),
        (1, 1): (0, 1),
    }
    params = aggregate_orbital_params(
        F, mapping, n_orbital_idx=2,
        epsilon_noise=1.0e-6, complex_type=0,
    )
    assert np.all(params.imag == 0.0), params.imag
    # Real part still perturbed
    assert np.any(params.real != 0.5)
    assert np.any(params.real != 0.3)


def test_aggregate_rank_lift_noise_reproducible_with_seeded_rng():
    """Two runs with the same RNG seed must produce identical output."""
    F = np.eye(4, dtype=np.complex128) * 0.25
    mapping = {(i, j): (i * 4 + j, 1) for i in range(4) for j in range(4)}
    rng1 = np.random.default_rng(123)
    rng2 = np.random.default_rng(123)
    p1 = aggregate_orbital_params(
        F, mapping, n_orbital_idx=16, epsilon_noise=1.0e-6, rng=rng1,
    )
    p2 = aggregate_orbital_params(
        F, mapping, n_orbital_idx=16, epsilon_noise=1.0e-6, rng=rng2,
    )
    np.testing.assert_array_equal(p1, p2)
