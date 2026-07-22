"""Unit tests for tools/_uhfk_to_mvmc/boundary_input.py.

It covers:
  - absent-key -> all-PBC default (matches UHFk)
  - container-shape / element-type checks (bridge-side)
  - delegation to hwave.solver._apbc_phase.normalize_boundary_condition
  - eigen.npz twist_offset consistency guard
"""
import math

import numpy as np
import pytest

from tools._uhfk_to_mvmc.boundary_input import (
    BoundaryInputError,
    check_eigen_twist_consistency,
    normalize_boundary_condition_list,
)


def test_absent_key_defaults_to_all_pbc():
    theta = normalize_boundary_condition_list(None)
    assert theta == (0.0, 0.0, 0.0)


def test_valid_apbc_list_passes_through():
    theta = normalize_boundary_condition_list(["ap", "periodic", "periodic"])
    assert theta == pytest.approx((math.pi, 0.0, 0.0))


def test_rejects_non_list_container():
    with pytest.raises(BoundaryInputError, match="list or tuple"):
        normalize_boundary_condition_list("periodic")


def test_rejects_non_string_element():
    with pytest.raises(BoundaryInputError, match="index 1"):
        normalize_boundary_condition_list(["periodic", 0, "periodic"])


def test_rejects_wrong_outer_length_from_hwave_helper():
    with pytest.raises(ValueError):
        # Delegated to H-wave helper (raises ValueError, not BoundaryInputError)
        normalize_boundary_condition_list(["periodic", "periodic"])


def test_eigen_twist_consistency_passes_on_match():
    theta = (math.pi, 0.0, 0.0)
    twist_saved = np.array([0.5, 0.0, 0.0], dtype=np.float64)
    # twist_offset(theta_d = pi) = 0.5 (see hwave/_apbc_phase.twist_offset)
    check_eigen_twist_consistency(theta, twist_saved)


def test_eigen_twist_consistency_rejects_mismatch():
    theta = (0.0, 0.0, 0.0)
    twist_saved = np.array([0.5, 0.0, 0.0], dtype=np.float64)
    with pytest.raises(BoundaryInputError, match="twist_offset"):
        check_eigen_twist_consistency(theta, twist_saved)


def test_eigen_twist_consistency_missing_twist_offset_is_ok_when_pbc():
    # Legacy eigen.npz without twist_offset field; only PBC path accepts.
    theta = (0.0, 0.0, 0.0)
    check_eigen_twist_consistency(theta, None)


def test_eigen_twist_consistency_missing_twist_offset_rejects_apbc():
    theta = (math.pi, 0.0, 0.0)
    with pytest.raises(BoundaryInputError, match="does not contain twist_offset"):
        check_eigen_twist_consistency(theta, None)
