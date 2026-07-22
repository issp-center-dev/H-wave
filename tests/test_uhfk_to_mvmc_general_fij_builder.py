"""Tests for general_fij_builder InOrbitalGeneral F construction."""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.general_fij_builder import (
    compute_canonical_reps,
    validate_general_prerequisites,
)
from tools._uhfk_to_mvmc.general_fij_builder import build_pair_list
from tools._uhfk_to_mvmc.general_fij_builder import (
    build_fij_general,
    build_slater_orbitals,
)
from tools._uhfk_to_mvmc.fij_builder import build_amplitudes, build_fij_phys
from tools._uhfk_to_mvmc.partner_index import find_partner_rows


def _klist(n):
    return np.roll(np.arange(n) - (n // 2), -(n // 2))


def test_compute_canonical_reps_1d_pbc_l4():
    """PBC L=4: partners are (0, 0), (1, -1)=(1, 3), (2, 2), (3, 1)=(3, 1).
    Self-pair rows: 0, 2. Non-self: {1, 3}. Canonical for {1, 3} is min-lex
    of wavevector_index tuple: (1,0,0) < (-1,0,0)=(3-4=-1) in signed form."""
    wv = np.array([[v, 0, 0] for v in _klist(4)], dtype=np.int64)
    # For PBC, partner_rows for translation-invariance is m ≡ -n mod L per direction
    # (theta=0 gives shift=0). See
    # ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    # _klist(4) = [0, 1, -2, -1]. partner mapping under theta=0:
    # row 0 (n=0) -> partner (m=0) -> row 0 (self)
    # row 1 (n=1) -> partner (m=-1) -> row 3
    # row 2 (n=-2) -> partner (m=2 ≡ -2) -> row 2 (self)
    # row 3 (n=-1) -> partner (m=1) -> row 1
    partner_rows = np.array([0, 3, 2, 1], dtype=np.int64)

    canonical, self_pairs = compute_canonical_reps(partner_rows, wv)
    canonical_set = set(canonical)
    self_set = set(self_pairs)

    assert self_set == {0, 2}
    # Non-self {1, 3}: canonical is the lex-smaller wavevector_index tuple.
    # wv[1] = (1, 0, 0), wv[3] = (-1, 0, 0). Lex order treats -1 < 1, so
    # canonical for {1, 3} is row 3.
    assert canonical_set == {0, 2, 3}
    assert len(canonical) == 3  # 2 self + 1 non-self representative


def test_validate_general_prerequisites_rejects_mixed_block():
    """column_spin containing -1 (mixed block, C case) → ValueError."""
    stepped = np.zeros((2, 4), dtype=np.float64)
    column_spin = np.array([0, 1, -1, 1], dtype=np.int64)
    partner_rows = np.array([0, 1], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    with pytest.raises(ValueError, match="mixed-block"):
        validate_general_prerequisites(
            Ncond=2, stepped_occupation=stepped, column_spin=column_spin,
            partner_rows=partner_rows, wavevector_index=wavevector_index,
        )


def test_validate_rejects_mixed_column_spin_when_not_soc_mode():
    """With is_soc_mode=False, column_spin=-1 triggers ValueError."""
    with pytest.raises(ValueError, match="mixed-block"):
        validate_general_prerequisites(
            Ncond=2,
            stepped_occupation=np.array([[1.0, 0.0], [0.0, 1.0]]),
            column_spin=np.array([-1, -1], dtype=np.int64),
            partner_rows=np.array([1, 0], dtype=np.int64),
            wavevector_index=np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64),
            is_soc_mode=False,
        )


def test_validate_accepts_mixed_column_spin_when_soc_mode():
    """With is_soc_mode=True, column_spin=-1 is allowed."""
    validate_general_prerequisites(
        Ncond=2,
        stepped_occupation=np.array([[1.0, 0.0], [0.0, 1.0]]),
        column_spin=np.array([-1, -1], dtype=np.int64),
        partner_rows=np.array([1, 0], dtype=np.int64),
        wavevector_index=np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64),
        is_soc_mode=True,
    )


def test_validate_rejects_unknown_column_spin_under_soc_mode():
    """SOC path: column_spin values outside {-1, 0, 1} rejected."""
    with pytest.raises(ValueError, match="not recognized"):
        validate_general_prerequisites(
            Ncond=2,
            stepped_occupation=np.array([[1.0, 0.0], [0.0, 1.0]]),
            column_spin=np.array([5, -3], dtype=np.int64),
            partner_rows=np.array([1, 0], dtype=np.int64),
            wavevector_index=np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64),
            is_soc_mode=True,
        )


def test_validate_soc_rejects_non_self_imbalance():
    """SOC rejects non-self (k, p) with unequal occupations.

    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    # Ncond=2 matches sum=2 so we pass the Ncond parity + sum guards and
    # reach the SOC balance check. Row 0 has 2 occupied columns, row 1
    # has 0, and partners are {0<->1} (non-self) so the imbalance must
    # trigger.
    with pytest.raises(ValueError, match="imbalance"):
        validate_general_prerequisites(
            Ncond=2,
            stepped_occupation=np.array([
                [1.0, 1.0, 0.0, 0.0],  # k=0: 2 occupied
                [0.0, 0.0, 0.0, 0.0],  # k=1: 0 occupied
            ]),
            column_spin=np.array([-1, -1, -1, -1], dtype=np.int64),
            partner_rows=np.array([1, 0], dtype=np.int64),
            wavevector_index=np.array(
                [[0, 0, 0], [1, 0, 0]], dtype=np.int64,
            ),
            is_soc_mode=True,
        )


def test_validate_soc_rejects_self_pair_odd():
    """SOC rejects a self-pair canonical k with odd n_occ(k).

    Two self-pair canonical blocks: k=0 has 3 occupied (ODD -> reject);
    k=1 has 1 occupied. Total Ncond=4 is even so the Ncond parity guard
    is satisfied and the SOC self-pair parity check fires. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    with pytest.raises(ValueError, match="odd"):
        validate_general_prerequisites(
            Ncond=4,
            stepped_occupation=np.array([
                [1.0, 1.0, 1.0, 0.0],  # k=0 self-pair: 3 occupied (ODD)
                [1.0, 0.0, 0.0, 0.0],  # k=1 self-pair: 1 occupied
            ]),
            column_spin=np.array([-1, -1, -1, -1], dtype=np.int64),
            partner_rows=np.array([0, 1], dtype=np.int64),
            wavevector_index=np.array(
                [[0, 0, 0], [1, 0, 0]], dtype=np.int64,
            ),
            is_soc_mode=True,
        )


def test_validate_general_prerequisites_rejects_odd_ncond():
    stepped = np.zeros((2, 4), dtype=np.float64)
    stepped[0, 0] = 1.0
    stepped[0, 2] = 1.0
    stepped[1, 0] = 1.0
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([0, 1], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    with pytest.raises(ValueError, match="odd"):
        validate_general_prerequisites(
            Ncond=3, stepped_occupation=stepped, column_spin=column_spin,
            partner_rows=partner_rows, wavevector_index=wavevector_index,
        )


def test_validate_general_prerequisites_rejects_excess_up_imbalance():
    """Non-self canonical block with n_excess_up_k != n_excess_up_p → reject."""
    # nvol_folded = 2, nd = 4 (2 up + 2 down cols). Partner: 0<->1 (non-self).
    # Minimally correct Ncond=2 (even, matches sum) so we actually reach the
    # excess-imbalance check rather than tripping the odd-Ncond guard first.
    stepped = np.zeros((2, 4), dtype=np.float64)
    stepped[0, 0] = 1.0  # up @ row 0, col 0
    stepped[0, 1] = 1.0  # up @ row 0, col 1
    # No down anywhere → n_cross_kd=0, n_cross_dk=0
    # #up@0 = 2, #up@1 = 0 → excess_up_k=2, excess_up_p=0 → violates
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([1, 0], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    with pytest.raises(ValueError, match="excess"):
        validate_general_prerequisites(
            Ncond=2, stepped_occupation=stepped, column_spin=column_spin,
            partner_rows=partner_rows, wavevector_index=wavevector_index,
        )


def test_validate_general_prerequisites_accepts_balanced_sz_free():
    """Sz-fixed 2Sz=0 closed shell: N up at k pairs with N down at partner(k)."""
    stepped = np.zeros((2, 4), dtype=np.float64)
    stepped[0, 0] = 1.0  # up @ row 0
    stepped[1, 2] = 1.0  # down @ row 1 = partner(0)
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([1, 0], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    validate_general_prerequisites(
        Ncond=2, stepped_occupation=stepped, column_spin=column_spin,
        partner_rows=partner_rows, wavevector_index=wavevector_index,
    )  # No exception


def test_build_pair_list_non_self_all_cross_2sz_zero():
    """Sz-fixed 2Sz=0: canonical block (0, partner=1). #up_k=#up_p=1,
    #down_k=#down_p=1, no excess. Expect 2 cross pairs:
      1. (up@k, down@partner)
      2. (up@partner, down@k)
    Ordering follows the canonical pair construction documented in
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    stepped = np.zeros((2, 4), dtype=np.float64)
    stepped[0, 0] = 1.0  # up @ row 0
    stepped[0, 2] = 1.0  # down @ row 0
    stepped[1, 0] = 1.0  # up @ row 1
    stepped[1, 2] = 1.0  # down @ row 1
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([1, 0], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    # Canonical of {0, 1} = row 0 (lex tie-break: (0,0,0) < (1,0,0)).
    canonical, _ = compute_canonical_reps(partner_rows, wavevector_index)
    assert canonical == [0]

    pairs = build_pair_list(stepped, column_spin, canonical, partner_rows)
    assert len(pairs) == 2
    # First: (up@k=0, down@partner=1)
    assert pairs[0] == {"alpha": (0, 0, "up"), "beta": (1, 2, "down")}
    # Second: (up@partner=1, down@k=0)
    assert pairs[1] == {"alpha": (1, 0, "up"), "beta": (0, 2, "down")}


def test_build_pair_list_non_self_same_spin_excess_2sz_positive():
    """Sz-fixed 2Sz=2 in a 2-row toy: canonical (0, 1). #up_k=1, #up_p=1,
    #down anywhere = 0. Expect 1 same-spin up-up pair (no cross)."""
    stepped = np.zeros((2, 4), dtype=np.float64)
    stepped[0, 0] = 1.0
    stepped[1, 0] = 1.0
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([1, 0], dtype=np.int64)
    wavevector_index = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    canonical, _ = compute_canonical_reps(partner_rows, wavevector_index)
    pairs = build_pair_list(stepped, column_spin, canonical, partner_rows)
    assert len(pairs) == 1
    assert pairs[0] == {"alpha": (0, 0, "up"), "beta": (1, 0, "up")}


def test_build_pair_list_self_pair_cross_only():
    """Self-pair k=0 with 1 up + 1 down occupied → single (up, down)
    cross pair, no excess (each spin count is 1, min = 1)."""
    stepped = np.zeros((1, 4), dtype=np.float64)
    stepped[0, 0] = 1.0  # up
    stepped[0, 2] = 1.0  # down
    column_spin = np.array([0, 0, 1, 1], dtype=np.int64)
    partner_rows = np.array([0], dtype=np.int64)  # self
    wavevector_index = np.array([[0, 0, 0]], dtype=np.int64)
    canonical, self_pairs = compute_canonical_reps(partner_rows, wavevector_index)
    assert canonical == [0] and self_pairs == [0]
    pairs = build_pair_list(stepped, column_spin, canonical, partner_rows)
    assert len(pairs) == 1
    assert pairs[0] == {"alpha": (0, 0, "up"), "beta": (0, 2, "down")}


def test_build_fij_general_antisymmetric():
    """F built from any Slater must satisfy F[i, j] = -F[j, i]."""
    # 2-site 1-orbital, synthetic single occupied pair.
    Ns = 2
    Ne = 2
    A = np.zeros((2 * Ns, Ne), dtype=np.complex128)
    # column 0: up orbital with amplitudes (0.7, 0.5) on the up rows
    A[0, 0] = 0.7
    A[1, 0] = 0.5
    # column 1: down orbital with amplitudes (0.3, 0.9) on the down rows
    A[2, 1] = 0.3
    A[3, 1] = 0.9

    F = build_fij_general(A)
    assert F.shape == (2 * Ns, 2 * Ns)
    for i in range(2 * Ns):
        for j in range(2 * Ns):
            assert F[i, j] == pytest.approx(-F[j, i], abs=1e-14)


def test_build_fij_general_matches_v21_antiparallel_closed_shell_pbc():
    """Sz-fixed 2Sz=0 PBC L=4 Ne=2 closed shell: General F[i↑, j↓] block
    must equal A_up @ A_down.T.

    Reuses the free-particle synthetic eigenvector that already covers
    the same case in test_uhfk_to_mvmc_fij_builder.py. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    L = 4
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array([[i, 0, 0] for i in range(L)], dtype=np.float64)
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)
    # nd = 2 (single orbital: 1 up + 1 down); v[k, 0, 0] = 1, v[k, 1, 1] = 1
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)

    stepped_occupation = np.zeros((L, 2), dtype=np.float64)
    # Occupy k=0 (both up and down) — closed-shell singlet at k=0
    n0 = list(wv[:, 0]).index(0)
    stepped_occupation[n0, 0] = 1.0
    stepped_occupation[n0, 1] = 1.0

    # Reference
    A_up_v21, A_down_v21 = build_amplitudes(
        wv, eigenvector, stepped_occupation, column_spin, column_mu_group,
        site_positions, norb_orig=1, theta=theta, L=L_vec,
    )
    F_v21 = build_fij_phys(A_up_v21, A_down_v21)  # (L, L)

    # General
    cell_shape = np.array([L, 1, 1], dtype=np.int64)
    subshape = np.array([1, 1, 1], dtype=np.int64)
    partner_rows, _ = find_partner_rows(wv, theta, L_vec)
    canonical, _ = compute_canonical_reps(partner_rows, wv)
    pairs = build_pair_list(stepped_occupation, column_spin, canonical, partner_rows)
    A = build_slater_orbitals(
        wv, eigenvector, column_spin,
        site_positions.astype(np.int64), cell_shape, subshape, theta, pairs,
    )
    F_gen = build_fij_general(A)
    assert F_gen.shape == (2 * L, 2 * L)
    # F_gen[i↑, j↓] block (top-right L x L) must match F
    F_gen_ud = F_gen[:L, L:]
    np.testing.assert_allclose(F_gen_ud, F_v21, atol=1e-12)
    # up-up and down-down blocks must be exactly 0 (no excess)
    np.testing.assert_allclose(F_gen[:L, :L], np.zeros((L, L)), atol=1e-14)
    np.testing.assert_allclose(F_gen[L:, L:], np.zeros((L, L)), atol=1e-14)


def test_build_fij_general_populates_up_up_block_for_2sz_positive():
    """Synthetic 2Sz=2 case (Ne_up=2, Ne_down=0) with up occupied at
    partner k=0 self-pair (must be even) + non-self (k=+π/4, -π/4)
    should populate F[up, up] block via same-spin excess pair."""
    L = 4
    wv = np.array([[v, 0, 0] for v in _klist(L)], dtype=np.int64)
    site_positions = np.array([[i, 0, 0] for i in range(L)], dtype=np.int64)
    theta = np.array([0.0, 0.0, 0.0], dtype=np.float64)
    L_vec = np.array([L, 1, 1], dtype=np.int64)
    eigenvector = np.zeros((L, 2, 2), dtype=np.complex128)
    for k in range(L):
        eigenvector[k, 0, 0] = 1.0
        eigenvector[k, 1, 1] = 1.0
    column_spin = np.array([0, 1], dtype=np.int64)

    # PBC L=4 partner map: rows are wv=[0, 1, -2, -1]. Partner under theta=0:
    #   row 0 (k=0)  ↔ row 0  (self)
    #   row 1 (k=1)  ↔ row 3  (k=-1)
    #   row 2 (k=-2) ↔ row 2  (self)
    # Occupy up at rows 1 and 3 (non-self canonical block).
    stepped = np.zeros((L, 2), dtype=np.float64)
    stepped[1, 0] = 1.0
    stepped[3, 0] = 1.0

    partner_rows, _ = find_partner_rows(wv, theta, L_vec)
    canonical, _ = compute_canonical_reps(partner_rows, wv)
    # Canonical for {1, 3}: wv[1]=(1,0,0), wv[3]=(-1,0,0). Lex smaller = -1 → row 3.
    # Self-pair rows 0 and 2 also in canonical.
    assert 3 in canonical
    pairs = build_pair_list(stepped, column_spin, canonical, partner_rows)
    # Only one same-spin up-up pair expected (from the non-self block).
    assert len(pairs) == 1
    assert pairs[0]["alpha"][2] == "up" and pairs[0]["beta"][2] == "up"

    cell_shape = np.array([L, 1, 1], dtype=np.int64)
    subshape = np.array([1, 1, 1], dtype=np.int64)
    A = build_slater_orbitals(
        wv, eigenvector, column_spin, site_positions,
        cell_shape, subshape, theta, pairs,
    )
    F = build_fij_general(A)
    # Up-up block (top-left L x L) should be non-zero somewhere.
    up_up = F[:L, :L]
    assert np.max(np.abs(up_up)) > 1e-3
    # Up-down block should be zero (no cross pair).
    assert np.max(np.abs(F[:L, L:])) < 1e-14
    # Antisymmetry: F.T == -F
    np.testing.assert_allclose(F.T, -F, atol=1e-14)


def test_build_slater_orbitals_soc_spin_block_permutation():
    """Under is_soc_mode=True, A rows are indexed r_phys + spin * Ns_phys
    (site-major, spin-minor). Choose Ns_phys=2 so two sites share the same
    intra-cell orbital component and any orbital-based row aliasing
    would collapse them; assert every distinct row is populated.

    The two pair members here pull from different k_rows (k=0 and k=1)
    with the same column index 0 — the new SOC branch reads eigenvector
    at the specified k_row, so populating (k=0, col=0) at column 0 of A
    and (k=1, col=0) at column 1 of A must not silently fall back to
    summing across k.
    """
    import numpy as np
    from tools._uhfk_to_mvmc.general_fij_builder import build_slater_orbitals

    Ns_phys = 2
    # 2 folded k-points, 2 spin-orbital rows per k (a_folded=0, spin=0/1)
    eigenvector = np.zeros((2, 2, 2), dtype=complex)
    # Inject amplitudes so each (r_phys, spin) row gets a unique amplitude
    eigenvector[0, 0, 0] = 1.0 + 0.0j  # k=0, packed idx 0 (a_folded=0, s=0)
    eigenvector[0, 1, 0] = 2.0 + 0.0j  # k=0, packed idx 1 (a_folded=0, s=1)
    eigenvector[1, 0, 0] = 3.0 + 0.0j  # k=1, packed idx 0
    eigenvector[1, 1, 0] = 4.0 + 0.0j  # k=1, packed idx 1
    column_spin = np.array([-1, -1], dtype=np.int64)  # mixed
    site_R_int = np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64)
    A = build_slater_orbitals(
        wavevector_index=np.array([[0, 0, 0], [1, 0, 0]], dtype=np.int64),
        eigenvector=eigenvector,
        column_spin=column_spin,
        site_positions=site_R_int,
        cell_shape=np.array([2, 1, 1], dtype=np.int64),
        subshape=np.array([1, 1, 1], dtype=np.int64),
        theta=np.zeros(3, dtype=np.float64),
        # single pair between (k=0, col=0) and (k=1, col=0)
        pair_list=[((0, 0), (1, 0))],
        is_soc_mode=True,
    )
    # A shape: (2 * Ns_phys, 2 * len(pair_list)) = (4, 2)
    assert A.shape == (4, 2)
    # Every row of A must have a nonzero amplitude - collapse would leave
    # a row all-zero.
    for row in range(4):
        assert np.any(np.abs(A[row]) > 1e-12), (
            f"row {row} = r_phys + spin * Ns_phys collapsed; "
            "check permutation formula"
        )


def test_build_pair_list_v3_spin_aware_regression():
    """A/B fixtures produce byte-for-byte identical pair lists to the
    frozen goldens; guards against the SOC refactor accidentally routing
    the non-SOC path through the SOC rule. The pair-list inputs are loaded
    from tracked snapshots under ``tests/data/``, so a clean-checkout CI run
    without the gitignored ``work/`` directory still exercises the test.
    """
    import json
    from pathlib import Path

    repo_root = Path(__file__).resolve().parent.parent
    data_dir = repo_root / "tests/data"
    goldens = json.loads(
        (data_dir / "v3_pair_list_regression.json").read_text()
    )
    for case in ("case_pbc_sz2", "case_zeeman_sz_free"):
        inputs = np.load(data_dir / f"v3_pair_list_input_{case}.npz")
        partner_rows, _ = find_partner_rows(
            inputs["wavevector_index"], np.zeros(3),
            inputs["cell_shape"],
        )
        canonical, _ = compute_canonical_reps(
            partner_rows, inputs["wavevector_index"]
        )
        pair_list = build_pair_list(
            inputs["occupation"], inputs["column_spin"],
            canonical, partner_rows,
        )

        def to_jsonable(pair):
            if isinstance(pair, dict):
                return {
                    k: (
                        list(v)
                        if isinstance(v, (list, tuple, np.ndarray))
                        else v
                    )
                    for k, v in pair.items()
                }
            if isinstance(pair, (list, tuple)):
                return [int(x) for x in pair]
            raise TypeError(f"pair type {type(pair)!r}")

        observed = [to_jsonable(p) for p in pair_list]
        assert observed == goldens[case], (
            f"{case}: pair list changed!\n"
            f"Observed:\n{observed}\n"
            f"Golden:\n{goldens[case]}"
        )


def test_build_pair_list_soc_non_self_canonical_alpha_pairing():
    """Under is_soc_mode=True, non-self canonical block pairs the
    alpha-th occupied column at k with the alpha-th at partner(k)
    as documented in ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    Each pair member carries its k_row so
    build_slater_orbitals can read the eigenvector at the correct k.
    """
    stepped = np.array(
        [
            [1.0, 1.0, 0.0, 0.0],  # k=0: cols 0, 1 occupied
            [1.0, 1.0, 0.0, 0.0],  # k=1: cols 0, 1 occupied
        ],
        dtype=np.float64,
    )
    column_spin = np.array([-1, -1, -1, -1], dtype=np.int64)
    canonical = [0]
    partner_rows = np.array([1, 0], dtype=np.int64)
    pairs = build_pair_list(
        stepped, column_spin, canonical, partner_rows,
        is_soc_mode=True,
    )
    # 2 occupied columns at k, 2 at partner -> 2 alpha-alpha pairs;
    # each pair carries ((k_row_alpha, col_alpha), (k_row_beta, col_beta)).
    assert len(pairs) == 2
    assert pairs[0] == ((0, 0), (1, 0))
    assert pairs[1] == ((0, 1), (1, 1))


def test_build_pair_list_soc_self_pair_odd_rejected():
    """SOC self canonical block with odd n_occ raises ValueError.

    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    stepped = np.array([[1.0, 1.0, 1.0, 0.0]], dtype=np.float64)
    column_spin = np.array([-1, -1, -1, -1], dtype=np.int64)
    canonical = [0]
    partner_rows = np.array([0], dtype=np.int64)
    with pytest.raises(ValueError, match="odd"):
        build_pair_list(
            stepped, column_spin, canonical, partner_rows,
            is_soc_mode=True,
        )


def test_build_pair_list_soc_non_self_imbalance_rejected():
    """SOC non-self canonical block with n_occ(k) != n_occ(partner)
    raises ValueError. See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    stepped = np.array(
        [
            [1.0, 1.0, 0.0, 0.0],  # k=0: 2 occupied
            [1.0, 0.0, 0.0, 0.0],  # k=1: 1 occupied - imbalance
        ],
        dtype=np.float64,
    )
    column_spin = np.array([-1, -1, -1, -1], dtype=np.int64)
    canonical = [0]
    partner_rows = np.array([1, 0], dtype=np.int64)
    with pytest.raises(ValueError, match="imbalance"):
        build_pair_list(
            stepped, column_spin, canonical, partner_rows,
            is_soc_mode=True,
        )


def test_build_pair_list_soc_self_pair_consecutive_columns():
    """SOC self canonical k: consecutive occupied columns get paired
    (0, 1), (2, 3), ... . Each pair member carries its
    k_row (both == self-pair k here).

    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    stepped = np.array([[1.0, 1.0, 1.0, 1.0]], dtype=np.float64)
    column_spin = np.array([-1, -1, -1, -1], dtype=np.int64)
    canonical = [0]
    partner_rows = np.array([0], dtype=np.int64)
    pairs = build_pair_list(
        stepped, column_spin, canonical, partner_rows,
        is_soc_mode=True,
    )
    assert pairs == [((0, 0), (0, 1)), ((0, 2), (0, 3))]
