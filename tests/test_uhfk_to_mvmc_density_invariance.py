"""Verify density matrix is invariant under PairProduct gauges
(F scaling, column permutation, complex APBC orientation).

See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np

from tools._uhfk_to_mvmc.density_check import density_from_amplitudes


def test_density_invariant_under_global_A_scaling():
    """Scaling A by an overall unitary phase must not change G^sigma."""
    rng = np.random.default_rng(seed=42)
    Ns = 8
    Nocc = 3
    A = rng.standard_normal((Ns, Nocc)).astype(np.complex128)
    A += 1j * rng.standard_normal((Ns, Nocc))
    Q, _ = np.linalg.qr(A)
    A_orth = Q[:, :Nocc]

    G = density_from_amplitudes(A_orth)
    phase = np.exp(1j * 0.37)
    G_scaled = density_from_amplitudes(A_orth * phase)
    np.testing.assert_allclose(G, G_scaled, atol=1e-12)


def test_density_invariant_under_column_unitary_mix():
    """A → A @ U for unitary U (mixing occupied orbitals) leaves G^sigma
    unchanged (occupied subspace identity)."""
    rng = np.random.default_rng(seed=7)
    Ns = 8
    Nocc = 3
    A = rng.standard_normal((Ns, Nocc)).astype(np.complex128)
    A += 1j * rng.standard_normal((Ns, Nocc))
    Q, _ = np.linalg.qr(A)
    A_orth = Q[:, :Nocc]

    M = rng.standard_normal((Nocc, Nocc)).astype(np.complex128)
    M += 1j * rng.standard_normal((Nocc, Nocc))
    U, _ = np.linalg.qr(M)

    G = density_from_amplitudes(A_orth)
    G_mixed = density_from_amplitudes(A_orth @ U)
    np.testing.assert_allclose(G, G_mixed, atol=1e-12)
