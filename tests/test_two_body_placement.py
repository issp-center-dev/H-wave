"""A two-body row and a Transfer row with the same ``(R, a, b)`` act on the
same bond (issue #232).

Both file formats follow Wannier90 / RESPACK: a row ``(R, a, b)`` places
orbital ``a`` in the original cell and orbital ``b`` in the cell displaced by
``R`` -- ``<w_{a,0}|H|w_{b,R}>`` for ``hr.dat`` and ``W_ab(R)`` for
``zvo_ur.dat``. Only the RELATIVE placement of the two kinds of rows is
physical (a global mirror ``R -> -R`` changes labels, not results), so every
check below is mirror-invariant:

* UHFk against an independent real-space Hartree-Fock solver that contains
  no H-wave code, compared through the sorted magnitudes of the inter-orbital
  bond amplitudes;
* RPA (ring) against the static Hartree linear response of an independent
  real-space model, after fixing normalization and the momentum frame on the
  free response, where the two frames related by the mirror fit equally.

The models are chosen so that the two placements are different Hamiltonians:
an inter-orbital off-site row whose table is not symmetric under ``R -> -R``,
with a transpose-asymmetric hopping.
"""
import os
import tempfile
import unittest

import numpy as np
from scipy.optimize import brentq

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa
from hwave.solver.uhfk import UHFk


# ---------------------------------------------------------------------------
# UHFk: independent real-space Hartree-Fock
# ---------------------------------------------------------------------------

_L, _NORB, _T, _NCOND = 4, 2, 0.05, 8
_TRANSFER = {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
             (1, 0, 1): -1.0, (-1, 1, 0): -1.0}          # (r, a, b): value
#: ONE inter-orbital bond: orbital 0 at the origin, orbital 1 at +x -- the
#: bond of the strong Transfer row (+1, 0, 1) -- declared with its partner.
_ROWS = {(1, 0, 1): 0.4, (-1, 1, 0): 0.4}


def _idx(j, o, s):
    return ((j % _L) * _NORB + o) * 2 + s


def _quartic_terms(itype):
    """``(coeff, p, q, r, s)`` for ``coeff c+_p c_q c+_r c_s`` -- the operator
    each declared row stands for, with the counting the solvers use: a row
    and its partner name the SAME density operator for the density types
    (each row carries half of it), while for Exchange they are two distinct
    operators and for PairHop each other's Hermitian conjugate."""
    terms = []
    for (r, a, b), v in _ROWS.items():
        for j in range(_L):
            ja, jb = j, j + r                       # a at the origin, b at +r
            if itype == "CoulombInter":
                terms += [(v / 2, _idx(ja, a, s), _idx(ja, a, s), _idx(jb, b, t), _idx(jb, b, t))
                          for s in range(2) for t in range(2)]
            elif itype == "Hund":
                terms += [(-v / 2, _idx(ja, a, s), _idx(ja, a, s), _idx(jb, b, s), _idx(jb, b, s))
                          for s in range(2)]
            elif itype == "Ising":
                terms += [(v / 2 * (1 - 2 * s) * (1 - 2 * t), _idx(ja, a, s), _idx(ja, a, s),
                           _idx(jb, b, t), _idx(jb, b, t)) for s in range(2) for t in range(2)]
            elif itype == "Exchange":
                terms.append((v, _idx(ja, a, 0), _idx(jb, b, 0), _idx(jb, b, 1), _idx(ja, a, 1)))
            elif itype == "PairHop":
                terms.append((v, _idx(ja, a, 0), _idx(jb, b, 0), _idx(ja, a, 1), _idx(jb, b, 1)))
            else:
                raise ValueError(itype)
    return terms


def _reference_bonds(itype):
    """Converged ``<c+_i c_j>`` of the independent solver, reduced to the
    sorted magnitudes of the two inter-orbital nearest-neighbour bonds."""
    m = _L * _NORB * 2
    h0 = np.zeros((m, m))
    for (r, a, b), v in _TRANSFER.items():           # v c+_{j,a} c_{j+r,b}
        for j in range(_L):
            for s in range(2):
                h0[_idx(j, a, s), _idx(j + r, b, s)] += v
    terms = _quartic_terms(itype)

    def mean_field(rho):                             # d<H_int>/d rho, Hartree + Fock
        f = np.zeros((m, m), complex)
        for v, i1, i2, i3, i4 in terms:
            f[i1, i2] += v * rho[i3, i4]
            f[i3, i4] += v * rho[i1, i2]
            f[i1, i4] -= v * rho[i3, i2]
            f[i3, i2] -= v * rho[i1, i4]
            if i2 == i3:
                f[i1, i4] += v
        return 0.5 * (f + f.conj().T)

    rho = np.zeros((m, m), complex)
    for _ in range(5000):
        w, vec = np.linalg.eigh(h0 + mean_field(rho))
        mu = brentq(lambda x: (1 / (np.exp((w - x) / _T) + 1)).sum() - _NCOND,
                    w.min() - 5, w.max() + 5)
        new = (vec.conj() / (np.exp((w - mu) / _T) + 1)) @ vec.T
        if np.abs(new - rho).max() < 1e-13:
            rho = new
            break
        rho = 0.5 * rho + 0.5 * new
    else:
        raise AssertionError("independent Hartree-Fock reference did not converge")
    return sorted([abs(rho[_idx(0, 0, 0), _idx(1, 1, 0)]),
                   abs(rho[_idx(1, 0, 0), _idx(0, 1, 0)])])


def _uhfk_bonds(itype):
    ham = {"Geometry": {"norb": _NORB, "rvec": np.eye(3), "center": np.zeros((_NORB, 3)),
                        "degree": np.ones(_NORB)},
           "Transfer": {((r, 0, 0), (a, b)): v for (r, a, b), v in _TRANSFER.items()},
           itype: {((r, 0, 0), (a, b)): v for (r, a, b), v in _ROWS.items()}}
    param = {"CellShape": [_L, 1, 1], "SubShape": [1, 1, 1], "Ncond": _NCOND, "2Sz": None,
             "IterationMax": 3000, "eps": 1e-12, "Mix": 0.5, "RndSeed": 7, "T": _T}
    solver = UHFk(ham, {"print_level": 0}, {"mode": "UHFk"}, param)
    solver.solve({"initial_mode": "random"}, ".")
    green = solver.Green
    return sorted([abs(green[1, 0, 0, 0, 1]), abs(green[-1, 0, 0, 0, 1])])


class TestUHFkPlacement(unittest.TestCase):
    TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop")

    def test_reference_sees_the_placement(self):
        """Anti-vacuity: the interaction moves the strong bond away from its
        free value by far more than the tolerance below, for every type."""
        free = _reference_bonds_free()
        for itype in self.TYPES:
            with self.subTest(type=itype):
                self.assertGreater(abs(_reference_bonds(itype)[1] - free[1]), 1e-3)

    def test_converged_state_matches_independent_hartree_fock(self):
        for itype in self.TYPES:
            with self.subTest(type=itype):
                np.testing.assert_allclose(_uhfk_bonds(itype), _reference_bonds(itype),
                                           rtol=0, atol=1e-8)


def _reference_bonds_free():
    saved = dict(_ROWS)
    try:
        for k in _ROWS:
            _ROWS[k] = 0.0
        return _reference_bonds("CoulombInter")
    finally:
        _ROWS.update(saved)


# ---------------------------------------------------------------------------
# RPA: independent real-space Hartree linear response
# ---------------------------------------------------------------------------

_RL, _RT, _RMU, _RV = 3, 0.5, 0.2, 0.3
_HOP = {(0, 0): 0.7 + 0.3j, (1, 1): 0.5 - 0.2j, (0, 1): 0.25 + 0.15j, (1, 0): -0.4 + 0.05j}
_EPS = (0.10, -0.05)
_RM = _RL * _NORB * 2


def _ridx(j, o, s):
    return ((j % _RL) * _NORB + o) * 2 + s


def _run_rpa(v, nmat=2048, bond_channels=False):
    """Density-density static chiq[q, a, b] (spins summed) of the RPA solver,
    or with ``bond_channels`` the charge channel of the bond-resolved
    longitudinal susceptibility (ring + the off-site exchange crossing)."""
    nd = 2 * _NORB
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n2\n3\n1 1 1\n")
            for (a, b), t in _HOP.items():
                c = np.conj(_HOP[(b, a)])
                f.write(" 1 0 0 %d %d %.15f %.15f\n" % (a + 1, b + 1, t.real, t.imag))
                f.write("-1 0 0 %d %d %.15f %.15f\n" % (a + 1, b + 1, c.real, c.imag))
            for o, e in enumerate(_EPS):
                f.write(" 0 0 0 %d %d %.15f 0\n" % (o + 1, o + 1, e))
        with open(os.path.join(d, "coulombinter.dat"), "w") as f:
            f.write("hdr\n2\n2\n1 1\n 1 0 0 1 2 %.15f 0\n-1 0 0 2 1 %.15f 0\n" % (v, v))
        inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat"}
        io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
        param = {"T": _RT, "mu": _RMU, "CellShape": [_RL, 1, 1], "SubShape": [1, 1, 1],
                 "Nmat": nmat, "coeff_tail": 1.0}
        if bond_channels:
            param["longitudinal_bond_channels"] = True
        solver = solver_rpa.RPA(io.get_param("ham"), {},
                                {"mode": "RPA", "calc_scheme": "general", "param": param})
        green = io.get_param("green")
        os.makedirs(os.path.join(d, "out"))
        solver.solve(green, os.path.join(d, "out"))
        if bond_channels:
            cc = np.asarray(green["longitudinal_bond_chiq_c_static"])    # (q, a, b, c, d)
            return np.array([[[cc[q, a, a, b, b] for b in range(_NORB)]
                              for a in range(_NORB)] for q in range(_RL)])
        arr = np.asarray(green["chiq"])
    nf = arr.shape[0]
    x = arr.reshape(nf, _RL, nd, nd, nd, nd)[nf // 2]
    sb = lambda o, s: s * _NORB + o
    return np.array([[[sum(x[q, sb(a, s), sb(a, s), sb(b, t), sb(b, t)]
                           for s in range(2) for t in range(2))
                       for b in range(_NORB)] for a in range(_NORB)] for q in range(_RL)])


def _free_hamiltonian(mirror):
    """The Transfer rows as ``t c+_{j,a} c_{j+1,b}``, or the mirror image."""
    h = np.zeros((_RM, _RM), complex)
    for j in range(_RL):
        for (a, b), t in _HOP.items():
            for s in range(2):
                p, q = ((_ridx(j + 1, a, s), _ridx(j, b, s)) if mirror
                        else (_ridx(j, a, s), _ridx(j + 1, b, s)))
                h[p, q] += t
                h[q, p] += np.conj(t)
        for o, e in enumerate(_EPS):
            for s in range(2):
                h[_ridx(j, o, s), _ridx(j, o, s)] += e
    return h


def _density(h):
    w, u = np.linalg.eigh(h)
    return (u.conj() / (np.exp((w - _RMU) / _RT) + 1)) @ u.T


def _occupation(rho, j, o):
    return sum(rho[_ridx(j, o, s), _ridx(j, o, s)] for s in range(2)).real


def _mean_field(rho, v, bond, fock):
    """Mean field of ``v n_{j0, orb 0} n_{j1, orb 1}`` on ``bond(j) = (j0, j1)``:
    Hartree, plus the bond's exchange (Fock) crossing when ``fock``."""
    terms = [(_ridx(j0, 0, s), _ridx(j1, 1, t))
             for j0, j1 in (bond(j) for j in range(_RL)) for s in range(2) for t in range(2)]
    f = np.zeros((_RM, _RM), complex)
    for p, q in terms:                               # v c+_p c_p c+_q c_q
        f[p, p] += v * rho[q, q]
        f[q, q] += v * rho[p, p]
        if fock:
            f[p, q] -= v * rho[q, p]
            f[q, p] -= v * rho[p, q]
    return 0.5 * (f + f.conj().T)


def _linear_response(h, v, bond, eps=1e-6, fock=False):
    """``chi_ab(q) = d<n_a(q)>/dh`` for ``H' = h sum_j e^{-iqj} n_{j,b}``, the
    mean-field shift of the free state subtracted so the reference state is
    the free one (bare-G RPA): Hartree only is the ring, Hartree + Fock adds
    the bond's exchange ladder."""
    rho0 = _density(h)
    counter = _mean_field(rho0, v, bond, fock)

    def scf(pert):
        rho = rho0
        for _ in range(2000):
            new = _density(h + pert + _mean_field(rho, v, bond, fock) - counter)
            if np.abs(new - rho).max() < 1e-14:
                return new
            rho = 0.5 * rho + 0.5 * new
        raise AssertionError("independent linear-response reference did not converge")

    def field(weights, b):
        p = np.zeros((_RM, _RM), complex)
        for j in range(_RL):
            for s in range(2):
                p[_ridx(j, b, s), _ridx(j, b, s)] = weights[j]
        return p

    out = np.zeros((_RL, _NORB, _NORB), complex)
    for q in range(_RL):
        k = 2 * np.pi * q / _RL
        for b in range(_NORB):
            rc = scf(eps * field(np.cos(k * np.arange(_RL)), b))
            rs = scf(-eps * field(np.sin(k * np.arange(_RL)), b))
            for a in range(_NORB):
                na = lambda r: sum(np.exp(1j * k * j) * _occupation(r, j, a) for j in range(_RL))
                out[q, a, b] = ((na(rc) - na(rho0)) + 1j * (na(rs) - na(rho0))) / eps
    return out


class _ResponseFrame:
    """Shared machinery: fix normalization and the momentum frame on the free
    response (both mirror-related frames fit there), then compare the
    interacting solver output with the independent response for the two
    placements of the declared row."""
    BOND_CHANNELS = False
    FOCK = False

    @staticmethod
    def _model(mirror):
        return _free_hamiltonian(mirror)

    @classmethod
    def _fit_frame(cls, free):
        fits = []
        for mirror in (False, True):
            for qsign in (1, -1):
                m = _linear_response(cls._model(mirror), 0.0, lambda j: (j, j + 1))
                m = m[[(qsign * q) % _RL for q in range(_RL)]]
                c = np.vdot(m.ravel(), free.ravel()) / np.vdot(m.ravel(), m.ravel())
                fits.append((np.abs(free - c * m).max() / np.abs(free).max(), mirror, qsign, c))
        return min(fits, key=lambda f: f[0])

    @classmethod
    def setUpClass(cls):
        cls.free = _run_rpa(1e-12 if not cls.BOND_CHANNELS else 1e-9,
                            bond_channels=cls.BOND_CHANNELS)
        cls.full = _run_rpa(_RV, bond_channels=cls.BOND_CHANNELS)
        cls.frame = cls._fit_frame(cls.free)

    def _deviation(self, same_bond):
        err, mirror, qsign, c = self.frame
        if same_bond == mirror:
            bond = lambda j: (j + 1, j)       # orbital 0 at j+1, orbital 1 at j
        else:
            bond = lambda j: (j, j + 1)
        m = c * _linear_response(self._model(mirror), _RV, bond, fock=self.FOCK)
        m = m[[(qsign * q) % _RL for q in range(_RL)]]
        return np.abs(self.full - m).max() / np.abs(self.full - self.free).max()

    def test_frame_is_fixed_by_the_free_response(self):
        self.assertLess(self.frame[0], 1e-6)

    def test_matches_the_same_bond_response(self):
        self.assertLess(self._deviation(same_bond=True), 1e-5)

    def test_opposite_bond_is_rejected(self):
        """Anti-vacuity: the other placement misses by orders of magnitude."""
        self.assertGreater(self._deviation(same_bond=False), 0.1)


class TestRPARingPlacement(_ResponseFrame, unittest.TestCase):
    """The ring (``chiq``) is the Hartree-only response."""


class TestRPABondChannelPlacement(_ResponseFrame, unittest.TestCase):
    """``longitudinal_bond_channels``: the ring plus the bond's exchange
    crossing is the Hartree + Fock response -- this pins the bond-channel
    resolvers to the same placement as the ring's density vertex."""
    BOND_CHANNELS = True
    FOCK = True


_SF = 0.13 + 0.05j     # on-site spin flip on orbital 0, so the SO path is spinful


def _run_rpa_spin_orbital(v, nmat=2048):
    """The same model in ``enable_spin_orbital`` mode (Transfer indices
    ``2 * orb + spin``, the interaction file in physical orbitals) with an
    on-site spin flip; density-density static chiq[q, a, b], spins summed
    (the spin-orbital chiq slots are spin-block, ``s * norb + orb``)."""
    nd = 2 * _NORB
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % nd + "0 0 0\n" * nd)
        rows = []
        for (a, b), t in _HOP.items():
            c = np.conj(_HOP[(b, a)])
            for s in range(2):
                rows += [(1, 2 * a + s, 2 * b + s, t), (-1, 2 * a + s, 2 * b + s, c)]
        for o, e in enumerate(_EPS):
            rows += [(0, 2 * o + s, 2 * o + s, e) for s in range(2)]
        rows += [(0, 0, 1, _SF), (0, 1, 0, np.conj(_SF))]
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n%d\n3\n1 1 1\n" % nd)
            for r, i, j, x in rows:
                f.write("%2d 0 0 %d %d %.15f %.15f\n" % (r, i + 1, j + 1, x.real, x.imag))
        with open(os.path.join(d, "coulombinter.dat"), "w") as f:
            f.write("hdr\n2\n2\n1 1\n 1 0 0 1 2 %.15f 0\n-1 0 0 2 1 %.15f 0\n" % (v, v))
        inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat"}
        io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
        param = {"T": _RT, "mu": _RMU, "CellShape": [_RL, 1, 1], "SubShape": [1, 1, 1],
                 "Nmat": nmat, "coeff_tail": 1.0}
        solver = solver_rpa.RPA(io.get_param("ham"), {},
                                {"mode": "RPA", "calc_scheme": "general",
                                 "enable_spin_orbital": True, "param": param})
        green = io.get_param("green")
        os.makedirs(os.path.join(d, "out"))
        solver.solve(green, os.path.join(d, "out"))
        arr = np.asarray(green["chiq"])
    nf = arr.shape[0]
    x = arr.reshape(nf, _RL, nd, nd, nd, nd)[nf // 2]
    sb = lambda o, s: s * _NORB + o
    return np.array([[[sum(x[q, sb(a, s), sb(a, s), sb(b, t), sb(b, t)]
                           for s in range(2) for t in range(2))
                       for b in range(_NORB)] for a in range(_NORB)] for q in range(_RL)])


class TestRPASpinOrbitalPlacement(_ResponseFrame, unittest.TestCase):
    """``enable_spin_orbital``: the same placement as the spin-free ring."""

    @classmethod
    def setUpClass(cls):
        cls.free = _run_rpa_spin_orbital(1e-12)
        cls.full = _run_rpa_spin_orbital(_RV)
        cls.frame = cls._fit_frame(cls.free)

    @staticmethod
    def _model(mirror):
        h = _free_hamiltonian(mirror)
        for j in range(_RL):
            h[_ridx(j, 0, 0), _ridx(j, 0, 1)] += _SF
            h[_ridx(j, 0, 1), _ridx(j, 0, 0)] += np.conj(_SF)
        return h


class TestDocumentedRPAEquation(unittest.TestCase):
    """The RPA equation of ``docs/*/source/algorithm/rpa.rst`` against the
    same independent real-space response, in the documented conventions
    only -- no H-wave code is involved:

    * Transfer row ``(R, a, b)``: ``t c+_{j,a} c_{j+R,b}``, so ``H(k) =
      sum_R e^{+ikR} t(R)`` and ``c_k = N^{-1/2} sum_R e^{-ikR} c_R``;
    * ``rho(q) = sum_k c+_{k+q} c_k = sum_R e^{+iqR} n_R`` and the
      documented ``X^{aa,bb}(q) ~ <rho_a(-q); rho_b(q)>``;
    * the documented ``W^{bb,aa}_q = sum_r X_ab(r) e^{+iqr}`` for a
      density-density row ``(r, a, b)`` (first pair at ``-q``, second at
      ``+q``).

    The documented vertex matches the field pair of ``X0`` with the
    ``+q`` pair of ``W`` and the response pair with the ``-q`` pair:
    ``W_hat(q)[row b, col a] = W^{aa,bb}_q``, the two pairs exchanged --
    the same numbers as ``W^{bb,aa}_{-q}``. ``X = [I + X0 W_hat]^{-1} X0``
    then holds exactly (with the analytic scale ``1 / N_site`` of the
    unnormalized ``rho``), while the positional assignment of earlier
    editions, ``W_hat(q)[row b, col a] = W^{bb,aa}_q``, misses by order
    one (issue #232)."""

    @classmethod
    def setUpClass(cls):
        h = _free_hamiltonian(False)
        bond = lambda j: (j, j + 1)            # row (+1, 0, 1): n_{j,0} n_{j+1,1}

        def documented_x(v):
            out = _linear_response(h, v, bond)
            # out[q] = d<sum_j e^{+iqj} n_ja>/dh for H' = h sum_j e^{-iqj} n_jb;
            # the documented X(q) = <rho(-q); rho(q)> is -out at -q
            return -np.array([out[(-q) % _RL] for q in range(_RL)])

        cls.x0, cls.x = documented_x(0.0), documented_x(_RV)
        cls.qs = 2.0 * np.pi * np.arange(_RL) / _RL

    def _deviation(self, sign):
        def w_hat(q):                           # W^{bb,aa}_q at [row b, col a]
            w = np.zeros((2, 2), complex)
            w[1, 0] = _RV * np.exp(+1j * q)     # X_01(+1)
            w[0, 1] = _RV * np.exp(-1j * q)     # the partner X_10(-1)
            return w
        # sign = -1: W^{bb,aa}_{-q}, i.e. W^{aa,bb}_q (the documented,
        # pair-exchanged vertex); sign = +1: the positional assignment
        lam = 1.0 / _RL
        dev = max(np.abs(np.linalg.solve(np.eye(2) + lam * self.x0[q] @ w_hat(sign * self.qs[q]),
                                         self.x0[q]) - self.x[q]).max() for q in range(_RL))
        return dev / np.abs(self.x - self.x0).max()

    def test_the_documented_equation_holds(self):
        self.assertLess(self._deviation(-1), 1e-5)

    def test_the_positional_pair_assignment_is_rejected(self):
        self.assertGreater(self._deviation(+1), 0.1)


if __name__ == "__main__":
    unittest.main()
