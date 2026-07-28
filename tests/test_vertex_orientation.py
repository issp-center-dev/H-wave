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

_ERRSTATE = None


def setUpModule():
    """numpy's Accelerate-backed complex matmul emits spurious FP warnings on
    the +-1 operator matrices below.  Scope the suppression to this module so
    it cannot mask warnings in tests that run later in the same process."""
    global _ERRSTATE
    _ERRSTATE = np.seterr(all="ignore")


def tearDownModule():
    if _ERRSTATE is not None:
        np.seterr(**_ERRSTATE)

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
            # ...and the coefficient itself, not just the orientation: with the
            # unnormalised rho_a(q) = sum_i e^{-iq r_i} n_ia used here the chain
            # is exactly -(1/N_site) chi0 V^T chi0.  A uniform double-count of
            # the interaction would survive the free-scale fit above.
            np.testing.assert_allclose(
                lead[iq], -(1.0 / NSITE) * (c0 @ v.T @ c0),
                rtol=0.0, atol=2e-4 * np.max(np.abs(lead[iq])),
                err_msg="q index {}".format(iq))


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
                # the coefficient too: (2/N_site) V(q)^T for this fixture --
                # 2 because both (up,down) orderings of the density-density
                # term contribute, 1/N_site from the plane-wave normalisation.
                np.testing.assert_allclose(
                    amp, (2.0 / NSITE) * Vq.T, rtol=0.0, atol=1e-12)
                checked += 1
        self.assertGreater(checked, 0, "no discriminating (k, k') pair")


# --------------------------------------------------------------------------
# 3. the two builders must agree
# --------------------------------------------------------------------------

class TestBuildersAgreeOnTheDensityBlock(unittest.TestCase):
    """`rpa._make_ham_inter` and `sc._build_interaction_k` must produce the
    same matrix, checked on the CoulombInter density block.

    This compares against the RPA builder itself -- constructing a real
    `Interaction` and reducing its `ham_inter_q` -- rather than against a
    locally assumed orientation, so it tests the agreement the fix is about
    and not the assumption behind it.

    SCOPE. The reduction ``einsum("ksasatbtb->ksatb", ...)`` keeps only the
    elements whose two indices agree within each particle-hole pair, and the
    ``[:, :NORB, :NORB]`` slice takes the up/up spin block. That is the block
    a density-density interaction lives in, and it is what this series
    determined. Exchange, PairLift and PairHop occupy other spin or
    off-density slots and are NOT covered here -- their agreement with
    `rpa.py` is not asserted anywhere in this file.
    """

    def test_sc_builder_matches_the_rpa_builder_for_coulombinter(self):
        import os
        import tempfile

        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.sc as sc
        from hwave.solver.rpa import Interaction, Lattice

        with tempfile.TemporaryDirectory(prefix="hwave_i96_") as d:
            with open(os.path.join(d, "geom.dat"), "w") as fh:
                fh.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n2\n"
                         "   0.0 0.0 0.0\n   0.5 0.0 0.0\n")
            rows = ["   0    0    0    1    1   0.0 0.0",
                    "   1    0    0    1    1  -1.0 0.0",
                    "  -1    0    0    1    1  -1.0 0.0",
                    "   1    0    0    2    2  -1.0 0.0",
                    "  -1    0    0    2    2  -1.0 0.0"]
            with open(os.path.join(d, "transfer.dat"), "w") as fh:
                fh.write("t\n2\n9\n 1 1 1 1 1 1 1 1 1\n"
                         + "\n".join(rows) + "\n")
            ci = ["   {}    0    0    {}    {}   {:.12f}   0.0".format(
                      R, a + 1, b + 1, v)
                  for (R, a, b), v in sorted(BONDS.items())]
            with open(os.path.join(d, "coulombinter.dat"), "w") as fh:
                fh.write("CoulombInter\n2\n9\n 1 1 1 1 1 1 1 1 1\n"
                         + "\n".join(ci) + "\n")

            interaction = {"path_to_input": d, "Geometry": "geom.dat",
                           "Transfer": "transfer.dat",
                           "CoulombInter": "coulombinter.dat"}
            ham_param = read_input_k.QLMSkInput(
                {"path_to_input": d, "interaction": interaction}
            ).get_param("ham")
            info_mode = {"mode": "RPA",
                         "param": {"T": 1.0, "mu": 0.0,
                                   "CellShape": [NSITE, 1, 1],
                                   "SubShape": [1, 1, 1], "Nmat": 4},
                         "calc_scheme": "reduced"}
            lattice = Lattice(info_mode["param"])
            inter = Interaction(lattice, ham_param, info_mode)

            qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
            _, _, interactions = sc._read_interaction_files(
                {"file": {"input": {"path_to_input": d,
                                    "interaction": interaction}}})
            ik = sc._build_interaction_k(qs, np.array([0.0]), np.array([0.0]),
                                         interactions, NORB)

        ns = 2
        nd = NORB * ns
        rpa_V = np.einsum(
            "ksasatbtb->ksatb",
            inter.ham_inter_q.reshape(lattice.nvol, *(ns, NORB) * 4)
        ).reshape(lattice.nvol, nd, nd)[:, :NORB, :NORB]
        sc_V = ik["CoulombInter"].transpose(2, 3, 4, 0, 1).reshape(
            NSITE, NORB, NORB)

        # non-vacuous: the fixture is orbital-asymmetric, so the two
        # orientations really are different matrices
        self.assertGreater(
            np.max(np.abs(sc_V - sc_V.transpose(0, 2, 1))), 1e-6)
        np.testing.assert_allclose(sc_V, rpa_V, rtol=0.0,
                                   atol=1e-12 * np.max(np.abs(rpa_V)))
        # ...and both are the orientation the exact tests above selected
        want = np.array([v_of_q(q).T for q in qs])
        np.testing.assert_allclose(sc_V, want, rtol=0.0, atol=1e-12)


class TestAllInteractionTypesTransposed(unittest.TestCase):
    """`_to_k` is type-agnostic, so the transpose reaches every interaction it
    converts -- CoulombIntra, CoulombInter, Hund, Exchange, Ising, PairLift and
    PairHop -- but only CoulombInter is exercised above. The committed
    Kanamori fixtures are all orbital-symmetric, so a wrong orientation for the
    others would not show up anywhere.

    Asymmetric input is reachable: the readers do not validate orbital symmetry
    (#93), and while the FLEX general path rejects OFF-SITE two-body terms it
    accepts an asymmetric ON-SITE orbital matrix, which then reaches the MYO
    S/C builder.
    """

    # every type that rpa.py places through _append_inter, i.e. with the
    # density-density slots (b, b, a, a) the determination above is for
    TRANSPOSED = ("CoulombIntra", "CoulombInter", "Hund", "Exchange", "Ising",
                  "PairLift")
    # PairHop is placed by _append_pairhop with the slots (b, a, a, b) instead,
    # so the density-density result does not carry over and it is left alone
    NOT_TRANSPOSED = ("PairHop",)

    def test_density_types_are_stored_pair_transposed(self):
        import hwave.sc as sc

        qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
        vr = {((R, 0, 0), (a, b)): v for (R, a, b), v in BONDS.items()}
        want = np.array([v_of_q(q).T for q in qs])
        self.assertGreater(np.max(np.abs(want - want.transpose(0, 2, 1))), 1e-6)

        for itype in self.TRANSPOSED:
            with self.subTest(interaction=itype):
                ik = sc._build_interaction_k(
                    qs, np.array([0.0]), np.array([0.0]), {itype: vr}, NORB)
                got = ik[itype].transpose(2, 3, 4, 0, 1).reshape(
                    NSITE, NORB, NORB)
                np.testing.assert_allclose(got, want, rtol=0.0, atol=1e-12)

    def test_pairhop_is_deliberately_left_alone(self):
        """PairHop keeps the pre-#96 placement, and that is NOT a claim that it
        is correct.

        `rpa.py` places PairHop through `_append_pairhop`, with the slots
        (b, a, a, b) instead of `_append_inter`'s density-density (b, b, a, a),
        so the determination this series rests on does not transfer to it.  It
        is left untouched because changing it would be a guess; whether it
        should be transposed is open, and tracked separately.

        This test therefore pins the STATUS QUO only.  It deliberately does not
        assert agreement with `rpa.py` -- that is exactly what is unresolved.

        The fixture is ON-SITE: `_append_pairhop` discards every ``irvec !=
        (0,0,0)``, so an off-site PairHop fixture would be vacuous on the RPA
        side and must not be used to reason about the two builders.
        """
        import hwave.sc as sc

        P = {(0, 1): 1.0, (1, 0): 0.35}          # on-site, asymmetric
        vr = {((0, 0, 0), (a, b)): v for (a, b), v in P.items()}
        want = np.zeros((NORB, NORB), dtype=complex)
        for (a, b), v in P.items():
            want[a, b] = v                        # NOT transposed
        self.assertGreater(np.max(np.abs(want - want.T)), 1e-6)

        qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
        for itype in self.NOT_TRANSPOSED:
            with self.subTest(interaction=itype):
                ik = sc._build_interaction_k(
                    qs, np.array([0.0]), np.array([0.0]), {itype: vr}, NORB)
                got = ik[itype][:, :, 0, 0, 0]    # on-site: q-independent
                np.testing.assert_allclose(got, want, rtol=0.0, atol=1e-12)

    def test_asymmetric_on_site_input_reaches_the_sc_builders(self):
        """Guard the reachability claim rather than assuming it: an asymmetric
        ON-SITE inter-orbital entry survives into the S/C matrices, and makes
        the charge matrix non-symmetric under the orbital-pair transpose."""
        import hwave.sc as sc
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo

        kx = np.linspace(0, 2 * np.pi, 2, endpoint=False)
        kz = np.array([0.0])
        sym = {((0, 0, 0), (0, 1)): 1.0, ((0, 0, 0), (1, 0)): 1.0}
        asym = {((0, 0, 0), (0, 1)): 1.0, ((0, 0, 0), (1, 0)): 0.4}

        def charge(vr, builder):
            ik = sc._build_interaction_k(kx, kx, kz, {"CoulombInter": vr}, NORB)
            return builder(ik, NORB, 2, 2, 1)[1]

        for builder in (sc._build_sc_matrices_all_q, build_sc_matrices_myo):
            with self.subTest(builder=builder.__name__):
                c_sym = charge(sym, builder)
                c_asym = charge(asym, builder)
                self.assertTrue(
                    np.allclose(c_sym, np.swapaxes(c_sym, -1, -2)),
                    "symmetric input must give a pair-symmetric charge matrix")
                self.assertFalse(
                    np.allclose(c_asym, np.swapaxes(c_asym, -1, -2)),
                    "asymmetric on-site input does reach the S/C builders, so "
                    "the orientation is observable there")


if __name__ == "__main__":
    unittest.main()
