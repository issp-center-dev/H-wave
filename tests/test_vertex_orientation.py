#!/usr/bin/env python3

"""Pin the orbital orientation of the interaction vertex against exact
calculations that use no convention from this codebase.

Issue #96. ``rpa.py::_make_ham_inter`` and ``sc.py::_build_interaction_k`` both
turn an interaction file into a momentum-space vertex, and they disagreed: for
a bond geometry that is not symmetric under ``a <-> b`` the two produced
matrices that are orbital transposes of each other. For a real, R-symmetric
interaction they coincide element-wise, which is why it survived.

H-wave's own paper (arXiv:2308.00324) settles the convention: Eq.(12) defines
the tensor ``W_ij^{aa'bb'}`` with the first pair at site i and the second at
site j, Eq.(16) puts the first pair at momentum ``+q``, and Eq.(21) defines the
matrix used in the RPA equation as ``[W(q)]^{ab} = W_q^{ba}``. For a
density-density term ``W^{aabb} = V_ab(R)`` that gives
``[W(q)]^{(aa),(bb)} = V_ba(q)``: the matrix is ``V(q)^T``.

The two tests below check that independently of the paper, by exact
calculation on an explicit Hamiltonian:

* ``TestRPAVertexOrientation`` puts the interaction between OPPOSITE SPINS
  only. There is then no exchange diagram -- the two densities are different
  species -- and the cross-spin response vanishes identically at V = 0, so its
  leading term is exactly one bubble-vertex-bubble chain with nothing else
  mixed in.

* ``TestPairingVertexOrientation`` computes the bare pair-scattering amplitude
  in the two-electron, zero-total-momentum sector. That is exact at first
  order and needs no diagrams at all.

Both select ``V(q)^T``. These are slow-ish but small; keep them.
"""

import itertools
import unittest

import numpy as np

np.seterr(all="ignore")     # spurious Accelerate matmul warnings

NSITE = 3          # >= 3 so that q and -q are distinct
NORB = 2
BETA = 3.0

# v_ab(R) with v_ab(R) != v_ba(R); every bond declared from both ends so the
# Hamiltonian stays Hermitian
BONDS = {(+1, 0, 1): 1.0, (-1, 1, 0): 1.0,
         (+1, 1, 0): 0.35, (-1, 0, 1): 0.35,
         (+1, 0, 0): 0.5, (-1, 0, 0): 0.5}

T_INTRA = {0: 0.9 * np.exp(0.35j), 1: 0.6 * np.exp(-0.55j)}
T_INTER = 0.45 * np.exp(0.8j)
ONSITE = {0: 0.15, 1: -0.25}


def v_of_q(q):
    """V_ab(q) = sum_R V_ab(R) e^{-iqR}."""
    M = np.zeros((NORB, NORB), dtype=complex)
    for (R, a, b), v in BONDS.items():
        M[a, b] += v * np.exp(-1j * q * R)
    return M


# --------------------------------------------------------------------------
# 1. the RPA vertex, from the cross-spin response
# --------------------------------------------------------------------------

NSPIN = 2
NUP = NDN = 3


def _mode(i, a, s):
    return (i * NORB + a) * NSPIN + s


def _sector():
    ups = [_mode(i, a, 0) for i in range(NSITE) for a in range(NORB)]
    dns = [_mode(i, a, 1) for i in range(NSITE) for a in range(NORB)]
    out = []
    for su in itertools.combinations(ups, NUP):
        for sd in itertools.combinations(dns, NDN):
            m = 0
            for p in su + sd:
                m |= 1 << p
            out.append(m)
    return sorted(out)


STATES = _sector()
INDEX = {m: i for i, m in enumerate(STATES)}
DIM = len(STATES)


def _sgn(mask, p):
    return -1.0 if bin(mask & ((1 << p) - 1)).count("1") % 2 else 1.0


def _hop(mask, p, q):
    if not (mask >> q) & 1:
        return None
    s = _sgn(mask, q)
    m = mask & ~(1 << q)
    if (m >> p) & 1:
        return None
    return m | (1 << p), s * _sgn(m, p)


def _hamiltonian(V):
    H = np.zeros((DIM, DIM), dtype=complex)
    for idx, mask in enumerate(STATES):
        for s in range(NSPIN):
            for i in range(NSITE):
                for a in range(NORB):
                    p = _mode(i, a, s)
                    if (mask >> p) & 1:
                        H[idx, idx] += ONSITE[a]
                    for d, amp in ((+1, T_INTRA[a]), (-1, np.conj(T_INTRA[a]))):
                        r = _hop(mask, _mode((i + d) % NSITE, a, s), p)
                        if r:
                            H[INDEX[r[0]], idx] += amp * r[1]
                    other = 1 - a
                    amp = T_INTER if a == 0 else np.conj(T_INTER)
                    r = _hop(mask, _mode(i, other, s), p)
                    if r:
                        H[INDEX[r[0]], idx] += amp * r[1]
        # interaction between OPPOSITE SPINS ONLY -> no exchange diagram
        for (R, a, b), v in BONDS.items():
            for i in range(NSITE):
                if ((mask >> _mode(i, a, 0)) & 1 and
                        (mask >> _mode((i + R) % NSITE, b, 1)) & 1):
                    H[idx, idx] += V * v
    return H


def _rho(q, a, spin):
    d = np.array([sum(np.exp(-1j * q * i)
                      for i in range(NSITE) if (m >> _mode(i, a, spin)) & 1)
                  for m in STATES], dtype=complex)
    return np.diag(d)


def _static_chi(E, U, opA, opB):
    """int_0^beta dtau <A(tau) B(0)>, at zero bosonic frequency."""
    w = np.exp(-BETA * (E - E.min()))
    A = U.conj().T @ opA @ U
    B = U.conj().T @ opB @ U
    dE = E[:, None] - E[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        kern = (w[None, :] - w[:, None]) / dE
    deg = np.abs(dE) < 1e-10
    kern[deg] = (BETA * w[:, None] * np.ones_like(dE))[deg]
    return np.sum(A * B.T * kern) / w.sum()


def _measure(V, qs):
    E, U = np.linalg.eigh(_hamiltonian(V))
    uu = np.zeros((len(qs), NORB, NORB), dtype=complex)
    ud = np.zeros((len(qs), NORB, NORB), dtype=complex)
    for iq, q in enumerate(qs):
        rp = {(a, s): _rho(q, a, s) for a in range(NORB) for s in range(NSPIN)}
        rm = {(a, s): _rho(-q, a, s) for a in range(NORB) for s in range(NSPIN)}
        for a in range(NORB):
            for b in range(NORB):
                uu[iq, a, b] = _static_chi(E, U, rp[(a, 0)], rm[(b, 0)])
                ud[iq, a, b] = _static_chi(E, U, rp[(a, 0)], rm[(b, 1)])
    return uu, ud


class TestRPAVertexOrientation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
        cls.chi0, cls.ud0 = _measure(0.0, cls.qs)

    def test_cross_spin_response_vanishes_without_interaction(self):
        """The premise: at V = 0 the two spin species are uncorrelated, so the
        first-order term is the whole story at O(V)."""
        self.assertLess(float(np.max(np.abs(self.ud0[1:]))), 1e-12)

    def test_rpa_vertex_is_the_orbital_transpose(self):
        dV = 1.0e-4
        _, ud = _measure(dV, self.qs)
        lead = (ud - self.ud0) / dV
        for iq in range(1, len(self.qs)):     # q = 0 carries a disconnected part
            v = v_of_q(self.qs[iq])
            c0 = self.chi0[iq]
            res = {}
            for name, M in (("V", v), ("VT", v.T)):
                pred = c0 @ M @ c0
                c = np.vdot(pred, lead[iq]) / np.vdot(pred, pred)
                res[name] = (np.max(np.abs(lead[iq] - c * pred))
                             / np.max(np.abs(lead[iq])))
            self.assertLess(res["VT"], 1e-4, "q index {}".format(iq))
            self.assertGreater(res["V"], 1e-1, "q index {}".format(iq))


# --------------------------------------------------------------------------
# 2. the pairing vertex, from the bare pair-scattering amplitude
# --------------------------------------------------------------------------

class TestPairingVertexOrientation(unittest.TestCase):
    """<k' a up, -k' b down| H_int |k a up, -k b down>, exact at first order.

    The orbitals cannot change (n_a is diagonal in orbital), so the amplitude
    is a matrix in (a, b) and can be compared with V(q) directly.
    """

    @staticmethod
    def _pair(k, a, b):
        up = np.exp(1j * k * np.arange(NSITE)) / np.sqrt(NSITE)
        dn = np.exp(-1j * k * np.arange(NSITE)) / np.sqrt(NSITE)
        return np.outer(up, dn)

    @staticmethod
    def _apply(psi, a, b):
        out = np.zeros_like(psi)
        for i in range(NSITE):
            for j in range(NSITE):
                coeff = 0.0
                for (dR, c, d), v in BONDS.items():
                    if (dR % NSITE) == (j - i) % NSITE and (c, d) == (a, b):
                        coeff += v
                    if (dR % NSITE) == (i - j) % NSITE and (c, d) == (b, a):
                        coeff += v
                out[i, j] = coeff * psi[i, j]
        return out

    def test_pairing_vertex_is_the_orbital_transpose(self):
        ks = 2.0 * np.pi * np.arange(NSITE) / NSITE
        checked = 0
        for ik, k in enumerate(ks):
            for ikp, kp in enumerate(ks):
                q = kp - k
                Vq = v_of_q(q)
                if np.allclose(Vq, Vq.T, atol=1e-12):
                    continue                  # q = 0: the two coincide
                amp = np.zeros((NORB, NORB), dtype=complex)
                for a in range(NORB):
                    for b in range(NORB):
                        amp[a, b] = np.vdot(self._pair(kp, a, b),
                                            self._apply(self._pair(k, a, b),
                                                        a, b))
                res = {}
                for name, M in (("V", Vq), ("VT", Vq.T)):
                    c = np.vdot(M, amp) / np.vdot(M, M)
                    res[name] = (np.max(np.abs(amp - c * M))
                                 / np.max(np.abs(amp)))
                self.assertLess(res["VT"], 1e-10)
                self.assertGreater(res["V"], 1e-1)
                checked += 1
        self.assertGreater(checked, 0, "no discriminating (k, k') pair")


# --------------------------------------------------------------------------
# 3. the two builders must agree
# --------------------------------------------------------------------------

class TestBuildersAgree(unittest.TestCase):
    """`rpa._make_ham_inter` and `sc._build_interaction_k` must produce the
    same matrix, in the orientation the two tests above selected."""

    def test_sc_builder_matches_the_transposed_convention(self):
        import hwave.sc as sc

        qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
        vr = {((R, 0, 0), (a, b)): v for (R, a, b), v in BONDS.items()}
        ik = sc._build_interaction_k(qs, np.array([0.0]), np.array([0.0]),
                                     {"CoulombInter": vr}, NORB)
        built = ik["CoulombInter"].transpose(2, 3, 4, 0, 1).reshape(
            NSITE, NORB, NORB)
        want = np.array([v_of_q(q).T for q in qs])
        self.assertGreater(np.max(np.abs(want - want.transpose(0, 2, 1))), 1e-6)
        np.testing.assert_allclose(built, want, rtol=0.0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
