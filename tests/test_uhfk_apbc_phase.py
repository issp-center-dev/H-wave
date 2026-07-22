"""Unit tests for the APBC gauge-phase helper."""
import numpy as np
import pytest

from hwave.solver._apbc_phase import (
    normalize_boundary_condition,
    transfer_phase,
    inverse_gauge_phase,
    twist_offset,
    sublattice_offset,
    full_lattice_displacement,
)


def test_normalize_lowercase_full_words():
    assert normalize_boundary_condition(["periodic", "antiperiodic", "periodic"]) == (0.0, np.pi, 0.0)


def test_normalize_short_and_case():
    assert normalize_boundary_condition(["P", "AP", "p"]) == (0.0, np.pi, 0.0)
    assert normalize_boundary_condition(["Periodic", "AntiPeriodic", "ap"]) == (0.0, np.pi, np.pi)


def test_normalize_rejects_wrong_length():
    with pytest.raises(ValueError, match="length"):
        normalize_boundary_condition(["periodic", "periodic"])


def test_normalize_rejects_unknown_value():
    with pytest.raises(ValueError, match="unknown"):
        normalize_boundary_condition(["periodic", "foo", "periodic"])


def test_normalize_rejects_non_string():
    with pytest.raises(TypeError):
        normalize_boundary_condition([0, 1, 0])


def test_transfer_phase_all_periodic_is_one():
    L = np.array([4, 4, 1])
    theta = np.array([0.0, 0.0, 0.0])
    assert transfer_phase(np.array([1, 2, 0]), theta, L) == pytest.approx(1.0)


def test_transfer_phase_antiperiodic_unit_displacement():
    L = np.array([4, 1, 1])
    theta = np.array([np.pi, 0.0, 0.0])
    expected = np.exp(1j * np.pi / 4)
    assert transfer_phase(np.array([1, 0, 0]), theta, L) == pytest.approx(expected)


def test_transfer_phase_signed_displacement():
    L = np.array([4, 1, 1])
    theta = np.array([np.pi, 0.0, 0.0])
    pos = transfer_phase(np.array([1, 0, 0]), theta, L)
    neg = transfer_phase(np.array([-1, 0, 0]), theta, L)
    assert neg == pytest.approx(np.conj(pos))


def test_transfer_phase_multi_direction():
    L = np.array([4, 4, 1])
    theta = np.array([np.pi, np.pi, 0.0])
    expected = np.exp(1j * np.pi * (1.0 / 4.0 + 2.0 / 4.0))
    assert transfer_phase(np.array([1, 2, 0]), theta, L) == pytest.approx(expected)


def test_inverse_gauge_phase_diagonal_is_one():
    L = np.array([4, 4, 1])
    theta = np.array([np.pi, np.pi, 0.0])
    r = np.array([2, 3, 0])
    assert inverse_gauge_phase(r, r, theta, L) == pytest.approx(1.0)


def test_inverse_gauge_phase_value_and_sign():
    L = np.array([4, 4, 1])
    theta = np.array([np.pi, 0.0, 0.0])
    r_i = np.array([3, 1, 0])
    r_j = np.array([1, 1, 0])
    expected = np.exp(1j * np.pi * (3 - 1) / 4)
    assert inverse_gauge_phase(r_i, r_j, theta, L) == pytest.approx(expected)
    assert inverse_gauge_phase(r_j, r_i, theta, L) == pytest.approx(np.conj(expected))


def test_twist_offset_periodic_zero():
    assert tuple(twist_offset([0.0, 0.0, 0.0])) == (0.0, 0.0, 0.0)


def test_twist_offset_antiperiodic_half():
    assert tuple(twist_offset([np.pi, 0.0, np.pi])) == (0.5, 0.0, 0.5)


# ---- sublattice helpers (sublattice + APBC) ---------------------------------


def _hwave_encode_orbit(orig_orb, bx, by, bz, norb_orig, Bx, By, Bz, spin=0):
    """Mirror of UHFk._reshape_orbit_ / _reshape_orbit_spin.

    Used in tests below to generate input orbital indices the same way uhfk.py
    folds them, so the recovery via sublattice_offset is exercised on the
    actual encoding (not a hand-rolled clone of the formula).
    """
    if spin == 0:
        return orig_orb + norb_orig * (bx + Bx * (by + By * bz))
    return orig_orb + norb_orig * (bx + Bx * (by + By * (bz + Bz * spin)))


@pytest.mark.parametrize(
    "Bx,By,Bz,bx,by,bz,orig_orb,norb_orig",
    [
        (2, 1, 1, 0, 0, 0, 0, 1),
        (2, 1, 1, 1, 0, 0, 0, 1),
        (2, 2, 1, 1, 1, 0, 0, 1),
        (3, 2, 1, 2, 1, 0, 1, 2),
        (2, 2, 2, 1, 0, 1, 0, 1),
    ],
)
def test_sublattice_offset_roundtrip_no_spin(
    Bx, By, Bz, bx, by, bz, orig_orb, norb_orig
):
    """Encode (bx, by, bz) with H-wave's formula then recover with the helper."""
    aa = _hwave_encode_orbit(orig_orb, bx, by, bz, norb_orig, Bx, By, Bz)
    assert sublattice_offset(aa, norb_orig, (Bx, By, Bz)) == (bx, by, bz)


@pytest.mark.parametrize(
    "Bx,By,Bz,bx,by,bz,orig_orb,norb_orig,spin",
    [
        (2, 1, 1, 1, 0, 0, 0, 2, 0),
        (2, 1, 1, 1, 0, 0, 1, 2, 1),  # spin=1 must still recover the spatial part
        (2, 2, 1, 1, 1, 0, 0, 2, 1),
    ],
)
def test_sublattice_offset_roundtrip_spin_orbital_strips_spin(
    Bx, By, Bz, bx, by, bz, orig_orb, norb_orig, spin
):
    """Spin lives at a higher position than (bx, by, bz); recovery must ignore it."""
    aa = _hwave_encode_orbit(orig_orb, bx, by, bz, norb_orig, Bx, By, Bz, spin=spin)
    assert sublattice_offset(aa, norb_orig, (Bx, By, Bz)) == (bx, by, bz)


def test_full_lattice_displacement_subshape_111_matches_irvec():
    """With SubShape = [1, 1, 1], Delta = irvec (no sublattice offset)."""
    # Only one sublattice site (0, 0, 0), so alpha = beta = 0 with norb_orig = 1.
    delta = full_lattice_displacement(
        irvec=(1, -1, 0), alpha=0, beta=0, norb_orig=1, subshape=(1, 1, 1)
    )
    assert delta == (1, -1, 0)


def test_full_lattice_displacement_same_sublattice_site():
    """Same sublattice site at both ends -> Delta = irvec * SubShape."""
    norb_orig = 1
    Bx, By, Bz = 2, 2, 1
    # Encode source = target = sublattice site (1, 0, 0), orig_orb = 0
    alpha = _hwave_encode_orbit(0, 1, 0, 0, norb_orig, Bx, By, Bz)
    beta = alpha
    delta = full_lattice_displacement(
        irvec=(1, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta == (1 * Bx + 0, 0, 0)


def test_full_lattice_displacement_cross_sublattice():
    """Source and target on different sublattice sites within the supercell."""
    norb_orig = 1
    Bx, By, Bz = 2, 2, 1
    # source at sublattice site (0, 0, 0), target at sublattice site (1, 1, 0)
    alpha = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz)
    beta = _hwave_encode_orbit(0, 1, 1, 0, norb_orig, Bx, By, Bz)
    # Same supercell -> irvec = (0, 0, 0). Full displacement = (1, 1, 0).
    delta = full_lattice_displacement(
        irvec=(0, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta == (1, 1, 0)


def test_full_lattice_displacement_boundary_crossing_with_sublattice():
    """irvec = +1 supercell + cross-sublattice offset must compose additively."""
    norb_orig = 1
    Bx, By, Bz = 2, 1, 1
    # source at sublattice (1, 0, 0), target at sublattice (0, 0, 0)
    # one supercell forward -> full Delta = 1*2 + (0 - 1) = 1
    alpha = _hwave_encode_orbit(0, 1, 0, 0, norb_orig, Bx, By, Bz)
    beta = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz)
    delta = full_lattice_displacement(
        irvec=(1, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta == (1, 0, 0)


def test_full_lattice_displacement_z_direction_only():
    """SubShape stretches only in z; source/target differ in z within supercell."""
    norb_orig = 1
    Bx, By, Bz = 1, 1, 2
    # source at sublattice (0, 0, 1), target at sublattice (0, 0, 0)
    # in-supercell -> irvec_z = 0, Delta_z = 0*Bz + (0 - 1) = -1
    alpha = _hwave_encode_orbit(0, 0, 0, 1, norb_orig, Bx, By, Bz)
    beta = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz)
    delta = full_lattice_displacement(
        irvec=(0, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta == (0, 0, -1)


def test_full_lattice_displacement_z_boundary_crossing():
    """irvec_z = +1 supercell + cross-sublattice z-offset."""
    norb_orig = 1
    Bx, By, Bz = 1, 1, 2
    # source at sublattice (0, 0, 1), target at sublattice (0, 0, 0), one
    # supercell forward in z -> Delta_z = 1*2 + (0 - 1) = 1
    alpha = _hwave_encode_orbit(0, 0, 0, 1, norb_orig, Bx, By, Bz)
    beta = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz)
    delta = full_lattice_displacement(
        irvec=(0, 0, 1), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta == (0, 0, 1)


def test_full_lattice_displacement_full_3d_sublattice():
    """All three directions have SubShape > 1 and cross-sublattice offsets."""
    norb_orig = 1
    Bx, By, Bz = 2, 2, 2
    # source at sublattice (1, 1, 1), target at sublattice (0, 0, 0)
    # in-supercell -> Delta = (0 - 1, 0 - 1, 0 - 1) = (-1, -1, -1)
    alpha = _hwave_encode_orbit(0, 1, 1, 1, norb_orig, Bx, By, Bz)
    beta = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz)
    delta_in = full_lattice_displacement(
        irvec=(0, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta_in == (-1, -1, -1)

    # Same with irvec = (1, -1, 1):
    # Delta = (1*2 + (0-1), -1*2 + (0-1), 1*2 + (0-1)) = (1, -3, 1)
    delta_out = full_lattice_displacement(
        irvec=(1, -1, 1), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    assert delta_out == (1, -3, 1)


def test_full_lattice_displacement_z_with_spin_orbital():
    """Spin-orbital encoding must NOT shift the z component."""
    norb_orig = 2  # spin-orbital uses orig_orb in [0, norb_orig)
    Bx, By, Bz = 1, 1, 2
    # source at sublattice (0, 0, 1) with spin=1, target at sublattice (0, 0, 0) with spin=0
    alpha = _hwave_encode_orbit(0, 0, 0, 1, norb_orig, Bx, By, Bz, spin=1)
    beta = _hwave_encode_orbit(0, 0, 0, 0, norb_orig, Bx, By, Bz, spin=0)
    delta = full_lattice_displacement(
        irvec=(0, 0, 0), alpha=alpha, beta=beta,
        norb_orig=norb_orig, subshape=(Bx, By, Bz),
    )
    # Spin contributes to alpha/beta but lives ABOVE Bz; sublattice_offset
    # masks it out via (subl_idx // (Bx*By)) % Bz, so z-component stays correct.
    assert delta == (0, 0, -1)
