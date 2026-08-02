"""APBC input guards: SubShape and non-density-interaction rejection."""
import pytest

from hwave.solver.uhfk import UHFk


def _make_stub(param_mod, param_ham=None):
    s = object.__new__(UHFk)
    s.param_mod = dict(param_mod)
    s.param_ham = dict(param_ham) if param_ham is not None else {}
    return s


def _base_mod(boundary=None, subshape=None, cell=(4, 1, 1)):
    pm = {"CellShape": list(cell)}
    if boundary is not None:
        pm["BoundaryCondition"] = boundary
    if subshape is not None:
        pm["SubShape"] = subshape
    return pm


# ---- _init_lattice + BoundaryCondition path ---------------------------------


def test_default_bc_marks_periodic_and_zero_theta():
    s = _make_stub(_base_mod(boundary=None))
    s._init_lattice()
    assert s.boundary_periodic is True
    assert s.boundary_theta == (0.0, 0.0, 0.0)


def test_apbc_with_explicit_subshape_111_is_accepted():
    s = _make_stub(_base_mod(
        boundary=["antiperiodic", "periodic", "periodic"],
        subshape=[1, 1, 1],
    ))
    s._init_lattice()
    assert s.boundary_periodic is False


def test_apbc_with_omitted_subshape_is_accepted():
    """Default SubShape (= CellShape) is acceptable under APBC.

    The gauge phase is applied to Transfer in its pre-fold signed-irvec
    representation, so any sublattice choice is valid (including the
    degenerate single-supercell SubShape = CellShape case).
    """
    s = _make_stub(_base_mod(boundary=["antiperiodic", "periodic", "periodic"]))
    s._init_lattice()
    assert s.boundary_periodic is False


def test_apbc_with_nontrivial_subshape_is_accepted():
    """SubShape > [1, 1, 1] combined with APBC is supported."""
    s = _make_stub(_base_mod(
        boundary=["antiperiodic", "periodic", "periodic"],
        subshape=[2, 1, 1],
    ))
    s._init_lattice()
    assert s.boundary_periodic is False
    assert s.subshape == (2, 1, 1)


def test_bad_boundary_length_raises():
    s = _make_stub(_base_mod(
        boundary=["antiperiodic", "periodic"], subshape=[1, 1, 1]
    ))
    with pytest.raises(ValueError, match="length"):
        s._init_lattice()


def test_unknown_boundary_value_raises():
    s = _make_stub(_base_mod(
        boundary=["antiperiodic", "twisted", "periodic"], subshape=[1, 1, 1]
    ))
    with pytest.raises(ValueError, match="unknown"):
        s._init_lattice()


# ---- _check_apbc_interactions path -----------------------------------------


def _apbc_stub(non_density_term=None, density_term=None):
    """Stub that has run through _init_lattice (APBC active)."""
    s = object.__new__(UHFk)
    s.boundary_periodic = False
    s.boundary_theta = (3.141592653589793, 0.0, 0.0)
    s.param_ham = {}
    if non_density_term is not None:
        s.param_ham[non_density_term] = {((1, 0, 0), (0, 0, 0, 0)): 0.5}
    if density_term is not None:
        s.param_ham[density_term] = {((1, 0, 0), (0, 0)): 1.0}
    return s


@pytest.mark.parametrize("term", ["PairHop", "PairLift"])
def test_apbc_rejects_non_density_interaction(term):
    s = _apbc_stub(non_density_term=term)
    with pytest.raises(ValueError, match=term):
        s._check_apbc_interactions()


def test_pbc_allows_non_density_interaction():
    s = _apbc_stub(non_density_term="PairHop")
    s.boundary_periodic = True
    s.boundary_theta = (0.0, 0.0, 0.0)
    s._check_apbc_interactions()  # must not raise


@pytest.mark.parametrize(
    "term", ["CoulombIntra", "CoulombInter", "Coulomb", "Hund", "Ising", "Exchange"]
)
def test_apbc_allows_density_interaction(term):
    s = _apbc_stub(density_term=term)
    s._check_apbc_interactions()  # must not raise


def test_apbc_allows_pairhop_with_empty_dict():
    """An empty PairHop entry must not trip the guard (no terms => no risk)."""
    s = object.__new__(UHFk)
    s.boundary_periodic = False
    s.boundary_theta = (3.141592653589793, 0.0, 0.0)
    s.param_ham = {"PairHop": {}}
    s._check_apbc_interactions()  # must not raise
