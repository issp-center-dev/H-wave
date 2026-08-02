# tests/test_apbc_phase_normalize.py
import math

import pytest

from hwave.solver._apbc_phase import normalize_boundary_condition


def test_accepts_periodic_variants():
    for bc in (
        ["periodic", "periodic", "periodic"],
        ["p", "P", "Periodic"],
        [" periodic ", "p", "P"],
    ):
        theta = normalize_boundary_condition(bc)
        assert theta == (0.0, 0.0, 0.0)


def test_accepts_antiperiodic_variants():
    theta = normalize_boundary_condition(["antiperiodic", "AP", "ap"])
    assert theta == pytest.approx((math.pi, math.pi, math.pi))


def test_accepts_mixed_pbc_apbc():
    theta = normalize_boundary_condition(["ap", "periodic", "p"])
    assert theta == pytest.approx((math.pi, 0.0, 0.0))


def test_rejects_open_and_ambiguous():
    for bc in (
        ["open", "periodic", "periodic"],
        ["aperiodic", "periodic", "periodic"],
        ["", "periodic", "periodic"],
        ["antiperiodicX", "periodic", "periodic"],
        ["peri odic", "periodic", "periodic"],
    ):
        with pytest.raises(ValueError):
            normalize_boundary_condition(bc)


def test_rejects_wrong_length():
    with pytest.raises(ValueError):
        normalize_boundary_condition(["periodic", "periodic"])
    with pytest.raises(ValueError):
        normalize_boundary_condition(["p", "p", "p", "p"])


def test_boundary_theta_is_radians_not_twist_offset():
    """``normalize_boundary_condition`` MUST
    return radians (`0` or `pi`), NOT the dimensionless
    ``twist_offset = theta / (2*pi)`` (which would be `0` or `0.5`). This
    pin blocks a bridge implementation that accidentally passes
    ``eigen.npz["twist_offset"]`` into ``gauge_lift(boundary_theta=...)``
    and produces silently wrong APBC densities.

    See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    from hwave.solver._apbc_phase import twist_offset

    theta = normalize_boundary_condition(["antiperiodic", "periodic", "periodic"])
    # Radians in {0, pi} — the value that gauge_lift and trans_emit want.
    assert theta[0] == pytest.approx(math.pi, abs=1e-12)
    assert theta[1] == 0.0
    assert theta[2] == 0.0
    # In particular NOT 0.5 (which would be the twist_offset value under
    # `twist_offset(theta) = theta / (2*pi) = pi / (2*pi) = 0.5`).
    assert theta[0] != pytest.approx(0.5, abs=1e-6)

    # And confirm the twist_offset -> theta relationship so any future
    # implementer chasing the "dimensionless" branch has an explicit
    # anchor: theta = 2 * pi * twist_offset.
    to = twist_offset(theta)
    assert to[0] == pytest.approx(0.5, abs=1e-12)
    assert to[1] == 0.0
    assert to[2] == 0.0
    # Round-trip: 2 * pi * twist_offset(theta) == theta.
    assert 2.0 * math.pi * to[0] == pytest.approx(theta[0], abs=1e-12)
