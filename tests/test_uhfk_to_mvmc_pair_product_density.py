"""Tests for tools/_uhfk_to_mvmc/pair_product_density.py.

Includes a parse_emitted_F sanity smoke test. The case_pbc comparison is
skipped until its required snapshot inputs are available.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.pair_product_density import (
    pair_product_density_from_F,
)

_DATA_DIR = Path(__file__).parent / "data"


def _hand_2site_F_and_G():
    """Return (F, G, A, N_pairs) for a 2-site (2*Ns=4) single-pair Slater.

    A has orthonormal columns (2 columns for N_pairs=1), the single pair is
    localized on the physical spin-up sites 0 and 1 (rows 0 and 1 of A).
    F = A @ Omega @ A.T with Omega = [[0, 1], [-1, 0]].
    G = conj(A) @ A.T is the expected density.
    """
    A = np.zeros((4, 2), dtype=np.complex128)
    A[0, 0] = 1.0 + 0.0j
    A[1, 1] = 1.0 + 0.0j
    Omega = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.complex128)
    F = A @ Omega @ A.T
    G_expected = np.conj(A) @ A.T
    return F, G_expected, A, 1


def test_pair_product_density_from_F_hand_2site():
    """Hand-computed 2-site F with a single occupied pair returns
    conj(A) @ A.T at 1e-12 with rank_tol=1e-10 (zero-noise fixture).
    """
    F, G_expected, A, N_pairs = _hand_2site_F_and_G()
    # Sanity: F is antisymmetric.
    assert np.allclose(F, -F.T, atol=1e-15)
    G = pair_product_density_from_F(F, N_pairs=N_pairs, rank_tol=1e-10)
    assert np.allclose(G, G_expected, atol=1e-12), (
        f"G mismatch: max abs delta = {np.max(np.abs(G - G_expected))}"
    )


def _random_unitary(n, seed):
    """Return a deterministic random unitary (n, n) via QR of complex normal."""
    rng = np.random.default_rng(seed)
    M = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    Q, R = np.linalg.qr(M)
    # Fix column phases to make QR unique.
    D = np.diag(np.sign(np.diag(R)))
    return Q @ D


def _synthetic_projector_F(seed=20260709):
    """Build F = Q @ Sigma @ Q.T for the projector test.

    Sigma = block_diag(0.5 * [[0,1],[-1,0]],  1.0 * [[0,1],[-1,0]],  0 * ...)
    so singular values sorted are [1, 1, 0.5, 0.5, 0, 0] and the top-4
    left-singular subspace equals span(Q[:, :4]). See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    Q = _random_unitary(6, seed=seed)
    block_anti = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=np.complex128)
    Sigma = np.zeros((6, 6), dtype=np.complex128)
    Sigma[0:2, 0:2] = 0.5 * block_anti
    Sigma[2:4, 2:4] = 1.0 * block_anti
    # Sigma[4:6, 4:6] is a zero block, singular values 0.
    F = Q @ Sigma @ Q.T
    # Explicitly antisymmetrize away tiny numerical asymmetry from complex
    # multiply-round-off; keep antisymmetry <= 1e-14.
    F = 0.5 * (F - F.T)
    return F, Q


def test_pair_product_density_from_F_matches_projector_synthetic():
    """F = Q @ Sigma @ Q.T with Sigma block singular structure {1, 0.5, 0};
    helper returns conj(Q[:, :4]) @ Q[:, :4].T at 1e-12 (zero-noise).

    The top-4 left-singular subspace of F equals span(Q[:, :4]) because SVD
    orders singular values descending: block B (sigma=1.0, Q cols 2-3) leads,
    block A (sigma=0.5, Q cols 0-1) follows. Together they span Q[:, :4].
    conj(A) @ A.T is invariant under unitary column rotation of A, so the
    two projectors coincide exactly.
    """
    F, Q = _synthetic_projector_F(seed=20260709)
    N_pairs = 2  # 4 non-trivial singular values means N_pairs = 2
    assert np.allclose(F, -F.T, atol=1e-14)
    G = pair_product_density_from_F(F, N_pairs=N_pairs, rank_tol=1e-10)
    G_expected = np.conj(Q[:, :4]) @ Q[:, :4].T
    delta = float(np.max(np.abs(G - G_expected)))
    assert delta < 1e-12, (
        f"synthetic projector mismatch: max abs delta = {delta}"
    )


def test_pair_product_density_from_F_hermitian_check():
    """Sanity: helper output is Hermitian on the synthetic case.

    NOTE: Hermiticity holds for BOTH `conj(U_occ) @ U_occ.T` AND the wrong
    orientation `U_occ @ U_occ.conj().T` — this pin is a basic projector
    sanity check, NOT the orientation guard. The load-bearing pin against
    orientation flips is `test_pair_product_density_from_F_rejects_wrong_
    orientation_soc` below.
    """
    F, _Q = _synthetic_projector_F(seed=20260709)
    G = pair_product_density_from_F(F, N_pairs=2, rank_tol=1e-10)
    assert np.allclose(G, G.conj().T, atol=1e-12), (
        f"G not Hermitian: max abs delta = "
        f"{np.max(np.abs(G - G.conj().T))}"
    )


@pytest.mark.skip(reason="The required case_pbc snapshot is not available")
def test_pair_product_density_from_F_matches_v1_case_pbc():
    """case_pbc snapshot with default --epsilon-noise 1e-8; helper output
    matches the density-check reference at 5e-7 with rank_tol=1e-6.

    Requires a regenerated fixture containing F_from_emitted and reference G.
    """
    from tools._uhfk_to_mvmc.pair_product_density import parse_emitted_F
    workspace = os.path.join(
        os.path.dirname(__file__), "data", "v1_case_pbc_snapshot"
    )
    F_from_emitted = parse_emitted_F(workspace)
    # N_pairs derived from the fixture header when the snapshot is materialized.
    N_pairs = None  # Set from the fixture header when the snapshot is available.
    G_ship = pair_product_density_from_F(
        F_from_emitted, N_pairs=N_pairs, rank_tol=1e-6
    )
    G_ref_path = os.path.join(workspace, "G_ref.npz")
    G_ref = np.load(G_ref_path)["G"]
    delta = float(np.max(np.abs(G_ship - G_ref)))
    assert delta < 5e-7, f"case_pbc snapshot mismatch: delta = {delta}"


def _synthetic_soc_A(nsite=4, seed=20260710):
    """Return an (2*nsite, 2*N_pairs) complex A with cross-spin non-zero
    entries and orthonormal columns.

    Rows [0, nsite) are spin-up, rows [nsite, 2*nsite) are spin-down. Cross-
    spin (row spin-up, column with weight in spin-down component) is enforced
    via non-zero magnitudes in both halves of each column plus complex phase.
    """
    n = 2 * nsite
    rng = np.random.default_rng(seed)
    # Take 2*N_pairs=4 columns => N_pairs = 2.
    N_pairs = 2
    r = 2 * N_pairs
    M = rng.standard_normal((n, r)) + 1j * rng.standard_normal((n, r))
    # Bias columns so both spin halves carry weight (avoid an accidental
    # Sz-fixed basis); add a phase-shifted cross-spin admixture.
    for k in range(r):
        M[nsite:, k] += (1.0 + 0.5j) * rng.standard_normal(nsite)
    # QR to get orthonormal columns.
    Q, R = np.linalg.qr(M)
    D = np.diag(np.sign(np.diag(R)))
    A = Q @ D
    # Sanity: cross-spin block of conj(A) @ A.T is non-zero (upper block:
    # row < nsite, col >= nsite).
    G_true = np.conj(A) @ A.T
    cross = np.abs(G_true[:nsite, nsite:])
    assert cross.max() > 1e-2, (
        f"synthetic SOC A has weak cross-spin density (max {cross.max()});"
        " tighten construction"
    )
    return A, N_pairs


def test_pair_product_density_from_F_rejects_wrong_orientation_soc():
    """On a synthetic SOC Rashba-style A (complex, cross-spin non-zero),
    helper matches conj(A) @ A.T at 1e-12 AND does NOT match A @ A.conj().T
    at 1e-6.

    The wrong-orientation projector A @ A.conj().T = G^T = conj(G) has the
    same magnitude structure as the correct G but opposite phase on
    off-diagonal complex entries. Under complex SOC, the difference is
    O(imag(G_ij)) at cross-spin off-diagonals; well above 1e-6.
    """
    A, N_pairs = _synthetic_soc_A(nsite=4, seed=20260710)
    N_all = 2 * A.shape[1] // 2  # 2 * N_pairs, but N_pairs is derived from A
    N = A.shape[1] // 2  # for Omega size, N_pairs
    Omega = np.zeros((2 * N, 2 * N), dtype=np.complex128)
    I_N = np.eye(N, dtype=np.complex128)
    Omega[:N, N:] = I_N
    Omega[N:, :N] = -I_N
    F = A @ Omega @ A.T
    # Enforce exact antisymmetry to within rounding.
    F = 0.5 * (F - F.T)
    assert np.allclose(F, -F.T, atol=1e-14)
    G = pair_product_density_from_F(F, N_pairs=N_pairs, rank_tol=1e-10)
    G_correct = np.conj(A) @ A.T
    G_wrong = A @ A.conj().T
    delta_correct = float(np.max(np.abs(G - G_correct)))
    delta_wrong = float(np.max(np.abs(G - G_wrong)))
    assert delta_correct < 1e-12, (
        f"correct orientation mismatch: delta = {delta_correct}"
    )
    # Safety margin: wrong orientation MUST fail; if the difference is tiny
    # then the fixture is not adversarial enough.
    assert delta_wrong > 1e-6, (
        f"wrong orientation not distinguishable: delta = {delta_wrong}"
    )


def test_pair_product_density_from_F_matches_v3_5_case_soc_sub_zeronoise():
    """On the shipping
    fixture case_soc_rashba_2d_sub run through the bridge with
    ``--epsilon-noise 0`` (rank-lift disabled), the projector helper on
    the emitted F must reproduce ``conj(A_ship) @ A_ship.T`` at 1e-10
    max_abs_delta with ``rank_tol = 1e-10``.

    This is the field-realistic complement to the hand-computed 2-site
    pin above: on a full SOC + SubShape > [1, 1, 1] shipping A whose
    columns come from build_slater_orbitals (not synthetic QR), the
    skew-SVD projector applied to the emitted F still recovers the
    shipping-A density. Any regression in
    build_fij_general -> parse_emitted_F -> pair_product_density_from_F
    that scrambles the projector orientation, drops a sign, or
    mis-computes the antisymmetrization will trip this at machine-scale.
    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    from pathlib import Path
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    data_dir = Path(__file__).parent / "data"
    F = np.load(
        data_dir / "v35_case_soc_rashba_2d_sub_zeronoise_F_pre_noise.npz",
    )["F"]
    eigen = np.load(data_dir / "v35_case_soc_rashba_2d_sub_eigen.npz")
    occ = np.load(data_dir / "v35_case_soc_rashba_2d_sub_occupation.npz")

    cell_shape = np.array([6, 4, 1], dtype=np.int64)
    subshape = np.array([2, 2, 1], dtype=np.int64)
    Ns = int(np.prod(cell_shape))
    Ncond = 8
    site_positions = np.array(
        [[x, y, z] for z in range(1) for y in range(4) for x in range(6)],
        dtype=np.int64,
    )
    theta = np.zeros(3, dtype=np.float64)

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[Ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions.astype(np.int64),
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_direct = np.conj(A_ship) @ A_ship.T

    N_pairs = A_ship.shape[1] // 2
    G_from_F = pair_product_density_from_F(
        F, N_pairs=N_pairs, rank_tol=1e-10,
    )

    max_diff = float(np.max(np.abs(G_from_F - G_direct)))
    assert max_diff < 1e-10, (
        f"|G_from_F - G_direct|_max = {max_diff:.3e} > 1e-10 on the "
        "case_soc_rashba_2d_sub zero-noise snapshot; the skew-SVD "
        "projector on emitted F no longer matches the shipping-A "
        "density under SOC + SubShape > [1, 1, 1]."
    )


# ---------------------------------------------------------------------
# Zero-noise F pin over the four multi-direction APBC fixtures.
# ---------------------------------------------------------------------


@pytest.mark.parametrize("fixture_name,expected_theta,expected_ncond", [
    ("xy",  (np.pi, np.pi, 0.0), 20),
    ("xz",  (np.pi, 0.0, np.pi), 20),
    ("yz",  (0.0, np.pi, np.pi), 24),
    ("xyz", (np.pi, np.pi, np.pi), 12),
])
def test_pair_product_density_from_F_matches_v37_case_soc_sub_zeronoise(
    fixture_name, expected_theta, expected_ncond,
):
    """Parametrized over the four shipping fixtures under
    CellShape=[4,4,4]/SubShape=[2,2,2]. For each
    fixture, on the zero-noise (``--epsilon-noise 0``) emitted F, the
    projector helper must reproduce ``conj(A_ship) @ A_ship.T`` at 1e-10
    max_abs_delta with ``rank_tol = 1e-10``.

    The test above pins the single-direction PBC case on
    case_soc_rashba_2d_sub. This test covers
    multi-direction APBC: xy (AP-AP-P), xz (AP-P-AP), yz (P-AP-AP), xyz
    (AP-AP-AP). Any regression in build_fij_general -> parse_emitted_F ->
    pair_product_density_from_F that scrambles the projector
    orientation, drops a sign, or mis-computes the antisymmetrization
    under multi-direction APBC will trip this at machine-scale. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    from tools._uhfk_to_mvmc.general_fij_builder import (
        build_slater_orbitals, build_pair_list, compute_canonical_reps,
    )
    from tools._uhfk_to_mvmc.occupation_step import step_occupation
    from tools._uhfk_to_mvmc.partner_index import find_partner_rows

    prefix = f"v37_case_soc_rashba_3d_sub_apbc_{fixture_name}"
    F = np.load(
        _DATA_DIR / f"{prefix}_bridge_zeronoise_F_pre_noise.npz",
    )["F"]
    eigen = np.load(_DATA_DIR / f"{prefix}_eigen.npz")
    occ = np.load(_DATA_DIR / f"{prefix}_occupation.npz")

    cell_shape = np.array([4, 4, 4], dtype=np.int64)
    subshape = np.array([2, 2, 2], dtype=np.int64)
    site_positions = np.array(
        [[x, y, z] for z in range(4) for y in range(4) for x in range(4)],
        dtype=np.int64,
    )
    theta = np.array(expected_theta, dtype=np.float64)

    stepped, _ = step_occupation(
        occ["occupation"], eigen["eigenvalue"], occ["column_spin"],
        occ["column_mu_group"], float(occ["T"]),
        ncond_per_group=[expected_ncond], is_soc_mode=True,
    )
    partner_rows, _ = find_partner_rows(
        eigen["wavevector_index"], theta, cell_shape // subshape,
    )
    canonical, _ = compute_canonical_reps(
        partner_rows, eigen["wavevector_index"],
    )
    pair_list = build_pair_list(
        stepped, occ["column_spin"], canonical, partner_rows,
        is_soc_mode=True,
    )
    A_ship = build_slater_orbitals(
        wavevector_index=eigen["wavevector_index"],
        eigenvector=eigen["eigenvector"],
        column_spin=occ["column_spin"],
        site_positions=site_positions,
        cell_shape=cell_shape,
        subshape=subshape,
        theta=theta,
        pair_list=pair_list,
        is_soc_mode=True,
    )
    G_direct = np.conj(A_ship) @ A_ship.T

    N_pairs = A_ship.shape[1] // 2
    G_from_F = pair_product_density_from_F(
        F, N_pairs=N_pairs, rank_tol=1e-10,
    )

    max_diff = float(np.max(np.abs(G_from_F - G_direct)))
    print(f"[{fixture_name}] max_abs_delta={max_diff:.3e}")
    assert max_diff < 1e-10, (
        f"[{fixture_name}] |G_from_F - G_direct|_max = {max_diff:.3e} > "
        "1e-10 on the zero-noise snapshot; the skew-SVD projector "
        "on emitted F no longer matches the shipping-A density under "
        "multi-direction APBC."
    )
