"""APBC under SubShape > [1, 1, 1]: SubShape-independence + analytic match.

These tests drive the full _init_interaction -> _make_ham_trans pipeline on
object.__new__(UHFk) stubs. _init_interaction injects the gauge phase into
Transfer in its original signed-irvec representation, then sublattice-folds.
We check:

1. SubShape choice is internal: the eigenvalue spectrum of H(k) for the SAME
   physical problem must be identical regardless of SubShape.
2. Sublattice-folded spectra match the analytic free-fermion APBC closed
   form in 1D, 2D, and 3D (including a fully-sublatticed 3D case with
   antiperiodic boundary in every direction).
"""
import numpy as np
import pytest

from hwave.solver.uhfk import UHFk


def make_full_stub(cellshape, subshape, transfer_dict,
                   boundary_theta=(0.0, 0.0, 0.0), norb_orig=1):
    """UHFk stub minimal for _init_interaction + _make_ham_trans.

    Skips the rest of __init__ but populates everything those two methods
    read. Spin-symmetric (ns=2, enable_spin_orbital=False); single original
    orbital by default.
    """
    s = object.__new__(UHFk)
    s.cellshape = tuple(cellshape)
    s.cellvol = cellshape[0] * cellshape[1] * cellshape[2]
    s.subshape = tuple(subshape)
    s.subvol = subshape[0] * subshape[1] * subshape[2]
    s.has_sublattice = s.subvol > 1
    s.shape = (
        cellshape[0] // subshape[0],
        cellshape[1] // subshape[1],
        cellshape[2] // subshape[2],
    )
    s.nvol = s.shape[0] * s.shape[1] * s.shape[2]
    s.norb_orig = norb_orig
    s.norb = norb_orig * s.subvol
    s.ns = 2
    s.nd = s.norb * s.ns
    s.enable_spin_orbital = False
    # _init_orbit is skipped by this stub, so mirror the attribute it sets and
    # that _reshape_interaction reads. Non-SO here, so the physical-orbital
    # count equals norb_orig.
    s.norb_phys_orig = norb_orig
    s.boundary_theta = tuple(boundary_theta)
    s.boundary_periodic = all(t == 0.0 for t in boundary_theta)

    s.param_ham = {
        "Geometry": {
            "norb": norb_orig,
            "rvec": np.eye(3, dtype=np.float64),
            "center": np.zeros((norb_orig, 3), dtype=np.float64),
        },
        "Transfer": dict(transfer_dict),
    }
    return s


def spectrum(s):
    """Run the production pipeline and return sorted eigenvalue spectrum."""
    s._init_interaction()
    s._make_ham_trans()
    eigs = []
    for hk in s.ham_trans:
        eigs.extend(np.linalg.eigvalsh(hk).tolist())
    return np.sort(np.asarray(eigs, dtype=np.float64))


# -------- analytic free-fermion APBC spectra (reference) --------------------


def expected_1d(L, theta_x, t=1.0):
    spec = []
    for n in range(L):
        eps = -2.0 * t * np.cos((2 * np.pi * n + theta_x) / L)
        spec.extend([eps, eps])  # spin-degenerate
    return np.sort(np.asarray(spec))


def expected_2d(Lx, Ly, theta_x, theta_y, t=1.0):
    spec = []
    for nx in range(Lx):
        for ny in range(Ly):
            kx = (2 * np.pi * nx + theta_x) / Lx
            ky = (2 * np.pi * ny + theta_y) / Ly
            eps = -2.0 * t * (np.cos(kx) + np.cos(ky))
            spec.extend([eps, eps])
    return np.sort(np.asarray(spec))


def expected_3d(Lx, Ly, Lz, tx, ty, tz, t=1.0):
    spec = []
    for nx in range(Lx):
        for ny in range(Ly):
            for nz in range(Lz):
                kx = (2 * np.pi * nx + tx) / Lx
                ky = (2 * np.pi * ny + ty) / Ly
                kz = (2 * np.pi * nz + tz) / Lz
                eps = -2.0 * t * (np.cos(kx) + np.cos(ky) + np.cos(kz))
                spec.extend([eps, eps])
    return np.sort(np.asarray(spec))


def nn_1d(t=1.0):
    return {((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t}


def nn_2d(t=1.0):
    return {
        ((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t,
        ((0, 1, 0), (0, 0)): -t, ((0, -1, 0), (0, 0)): -t,
    }


def nn_3d(t=1.0):
    return {
        ((1, 0, 0), (0, 0)): -t, ((-1, 0, 0), (0, 0)): -t,
        ((0, 1, 0), (0, 0)): -t, ((0, -1, 0), (0, 0)): -t,
        ((0, 0, 1), (0, 0)): -t, ((0, 0, -1), (0, 0)): -t,
    }


# -------- 1) SubShape-independence ------------------------------------------


@pytest.mark.parametrize(
    "subshape", [[1, 1, 1], [2, 1, 1], [4, 1, 1]]
)
def test_subshape_independence_1d_apbc(subshape):
    """Same physical 1D L=4 APBC -> same spectrum across SubShape choices."""
    L = 4
    s = make_full_stub(
        [L, 1, 1], subshape, nn_1d(),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    got = spectrum(s)
    want = expected_1d(L, np.pi)
    np.testing.assert_allclose(got, want, atol=1e-12)


@pytest.mark.parametrize(
    "subshape", [[1, 1, 1], [2, 1, 1], [2, 2, 1], [4, 4, 1]]
)
def test_subshape_independence_2d_apbc(subshape):
    """Same physical 2D 4x4 APBC -> same spectrum across SubShape choices."""
    Lx, Ly = 4, 4
    s = make_full_stub(
        [Lx, Ly, 1], subshape, nn_2d(),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    got = spectrum(s)
    want = expected_2d(Lx, Ly, np.pi, 0.0)
    np.testing.assert_allclose(got, want, atol=1e-12)


# -------- 2) Sublattice + APBC analytic match (all directions) --------------


def test_sublattice_apbc_1d_full_supercell():
    """SubShape = CellShape (degenerate single-supercell fold) under APBC."""
    L = 6
    s = make_full_stub(
        [L, 1, 1], [L, 1, 1], nn_1d(),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    got = spectrum(s)
    want = expected_1d(L, np.pi)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_sublattice_apbc_2d_supercell_with_z_subshape():
    """2D APBC in x, with non-trivial SubShape spanning x and y."""
    Lx, Ly = 4, 4
    s = make_full_stub(
        [Lx, Ly, 1], [2, 2, 1], nn_2d(),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    got = spectrum(s)
    want = expected_2d(Lx, Ly, np.pi, 0.0)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_sublattice_apbc_3d_all_directions_with_z_subshape():
    """3D fully-APBC with SubShape > 1 in every direction (including z)."""
    Lx, Ly, Lz = 4, 4, 4
    s = make_full_stub(
        [Lx, Ly, Lz], [2, 2, 2], nn_3d(),
        boundary_theta=(np.pi, np.pi, np.pi),
    )
    got = spectrum(s)
    want = expected_3d(Lx, Ly, Lz, np.pi, np.pi, np.pi)
    np.testing.assert_allclose(got, want, atol=1e-12)


def test_sublattice_apbc_mixed_directions_with_z():
    """Mixed: APBC in y and z, PBC in x; non-trivial SubShape in y and z."""
    Lx, Ly, Lz = 2, 4, 4
    s = make_full_stub(
        [Lx, Ly, Lz], [1, 2, 2], nn_3d(),
        boundary_theta=(0.0, np.pi, np.pi),
    )
    got = spectrum(s)
    want = expected_3d(Lx, Ly, Lz, 0.0, np.pi, np.pi)
    np.testing.assert_allclose(got, want, atol=1e-12)
