"""Tests for occupation_step (T=0 step function from SCF occupation)."""
from __future__ import annotations

import numpy as np
import pytest

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc.occupation_step import (
    step_occupation,
    OccupationGuardError,
)


def test_step_occupation_t0_preserves_zero_one():
    """At T=0 the saved occupation is already 0/1; step_occupation must
    return it unchanged and report Ne consistent with the input."""
    nvol = 4
    nd = 2  # 1 up, 1 down (single orbital, Sz-fixed)
    # Block 0 (column 0, up): occupations [1, 1, 0, 0] over k=0..3 → Ne_up=2
    # Block 1 (column 1, down): same shape → Ne_down=2
    occupation = np.array([
        [1.0, 1.0],
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
    ], dtype=np.float64)
    eigenvalue = np.array([
        [-2.0, -2.0],
        [-1.0, -1.0],
        [+1.0, +1.0],
        [+2.0, +2.0],
    ], dtype=np.float64)
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)
    T = 0.0
    ncond_per_group = [2, 2]  # Ne_up = Ne_down = 2

    occ_step, summary = step_occupation(
        occupation, eigenvalue, column_spin, column_mu_group,
        T, ncond_per_group,
    )

    np.testing.assert_array_equal(occ_step, occupation)
    assert summary["ne_per_group"] == [2, 2]


def test_fractional_residual_fails_fast():
    """Two or more fractional occupations → raise."""
    nvol, nd = 4, 2
    eigenvalue = np.array([
        [-1.0, -1.0],
        [-0.01, -0.01],
        [+0.01, +0.01],
        [+1.0, +1.0],
    ], dtype=np.float64)
    # T=0.1 leaves f near 0.5 for the inner two k → 2 fractional per group
    occupation = np.array([
        [0.99, 0.99],
        [0.55, 0.55],
        [0.45, 0.45],
        [0.01, 0.01],
    ], dtype=np.float64)
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)
    T = 0.1
    ncond_per_group = [2, 2]

    with pytest.raises(OccupationGuardError, match="fractional"):
        step_occupation(
            occupation, eigenvalue, column_spin, column_mu_group,
            T, ncond_per_group,
        )


def test_fermi_level_degeneracy_fails_fast():
    """Equal HOMO and LUMO eigenvalues → raise."""
    nvol, nd = 4, 2
    eigenvalue = np.array([
        [-1.0, -1.0],
        [+0.0, +0.0],
        [+0.0, +0.0],  # degenerate at Fermi level (HOMO = LUMO = 0)
        [+1.0, +1.0],
    ], dtype=np.float64)
    occupation = np.array([
        [1.0, 1.0],
        [1.0, 1.0],
        [0.0, 0.0],
        [0.0, 0.0],
    ], dtype=np.float64)
    column_spin = np.array([0, 1], dtype=np.int64)
    column_mu_group = np.array([0, 1], dtype=np.int64)

    with pytest.raises(OccupationGuardError, match="degeneracy"):
        step_occupation(
            occupation, eigenvalue, column_spin, column_mu_group,
            0.0, [2, 2],
        )


def test_sz_free_column_spin_minus_one_fails_fast():
    """column_spin = -1 (Sz-free mixed block) must raise."""
    occupation = np.zeros((2, 2))
    eigenvalue = np.zeros((2, 2))
    column_spin = np.array([-1, -1], dtype=np.int64)  # mixed
    column_mu_group = np.array([0, 0], dtype=np.int64)

    with pytest.raises(OccupationGuardError, match="Sz-free"):
        step_occupation(
            occupation, eigenvalue, column_spin, column_mu_group,
            0.0, [2],
        )


def test_step_occupation_accepts_all_mixed_columns_under_soc_mode():
    """SOC accepts column_spin = -1 and occupies the lowest Ncond states.

    Note: SCF occupation values are 0/1 (T=0 Slater determinant) so
    the fractional-residual guard passes; the intent here is
    to verify the SOC-mode dispatch bypasses the ``column_spin < 0``
    reject and that eigenvalue ordering selects the lowest ``Ncond``
    states regardless of spin label.
    """
    # SCF occupation "misaligns" with the eig ordering below, so the
    # assertions below verify that step_occupation uses eig ordering
    # (not the input occupation) once the SOC-mode guard is passed.
    occ = np.array([[0.0, 0.0], [1.0, 1.0]])
    eig = np.array([[0.0, 0.1], [0.2, 0.3]])
    col_spin = np.array([-1, -1], dtype=np.int64)
    col_mu = np.array([0, 0], dtype=np.int64)
    stepped, summary = step_occupation(
        occ, eig, col_spin, col_mu, T=0.0,
        ncond_per_group=[2],
        is_soc_mode=True,
    )
    # Lowest 2 eigenvalues (0.0 at [0,0] and 0.1 at [0,1]) get occ=1
    assert stepped[0, 0] == 1.0
    assert stepped[0, 1] == 1.0
    assert stepped[1, 0] == 0.0
    assert stepped[1, 1] == 0.0
    assert summary["ne_per_group"] == [2]


def test_step_occupation_rejects_mixed_columns_without_soc_mode():
    """column_spin = -1 without is_soc_mode still raises."""
    with pytest.raises(OccupationGuardError, match="column_spin = -1"):
        step_occupation(
            np.zeros((2, 2)), np.zeros((2, 2)),
            np.array([-1, -1], dtype=np.int64),
            np.array([0, 0], dtype=np.int64),
            0.0, [2],
            is_soc_mode=False,
        )
