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

    # the types rpa.py places through _append_inter, i.e. with the
    # density-density slots (b, b, a, a) the determination above is for
    TRANSPOSED = ("CoulombIntra", "CoulombInter", "Hund", "Exchange", "Ising",
                  "PairLift")
    # PairHop is placed by _append_pairhop with the slots (b, a, a, b) instead,
    # so the density-density result does not carry over. It has its own
    # determination, TestPairHopEndToEnd below, which selects the UNtransposed
    # placement.
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

    def test_pairhop_is_not_stored_pair_transposed(self):
        """PairHop keeps the plain ``[a, b]`` placement, which
        TestPairHopEndToEnd measures against exact diagonalization.

        The fixture is ON-SITE: ``rpa.py::_append_pairhop`` discards every
        ``irvec != (0,0,0)``, so an off-site PairHop fixture would be vacuous on
        the RPA side and must not be used to reason about the two builders.
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


# --------------------------------------------------------------------------
# 5. PairHop: its own determination (issue #100)
# --------------------------------------------------------------------------

# PairHop is the one type `rpa.py` does not route through `_append_inter`, and
# in the S/C matrices it lands in the pair-ANTIdiagonal Case 4 rather than the
# density-density Case 1/3. The argument above therefore says nothing about it,
# and it needs two separate measurements:
#
#   1. what chi0's pair index means -- because relabelling a pair (x,y)->(y,x)
#      on one side turns a pair-diagonal vertex into a pair-antidiagonal one,
#      so the vertex alone cannot answer the question;
#   2. where the PairHop vertex sits when the pair index labels a bilinear.
#
# Both are done by exact diagonalization of an explicit Hamiltonian. They are
# INPUTS to the answer, not the answer: combining them by hand is exactly where
# this went wrong once (the vertex sits BETWEEN two chi0 factors, so chi0's
# ROW-pair reversal lands on the vertex's COLUMN side, not its row side), and
# two derivations of the same shape reached opposite conclusions. The placement
# is therefore settled by TestPairHopEndToEnd, which measures the whole chain
# and needs neither of them.

PH_NSITE = 3
PH_NORB = 2
PH_BETA = 3.0
PH_MU = 0.2
PH_NTAU = 96

PH_T_INTRA = {0: 0.9 * np.exp(0.35j), 1: 0.6 * np.exp(-0.55j)}
PH_T_INTER = 0.45 * np.exp(0.8j)
PH_ONSITE = {0: 0.15, 1: -0.25}


class TestChi0PairIndexMeaning(unittest.TestCase):
    """`rpa.py::_calc_chi0q` builds, in the general branch,

        chi0[a,c,b,d] = G[a,b] * G_rev[d,c] * sgn

    Which BILINEAR does each flattened pair stand for? Measure it: compute the
    exact correlator of the generalised densities

        A_{ab}(q) = sum_i e^{-i q r_i} c+_{i a} c_{i b}

    by grand-canonical diagonalization, evaluate the kernel above with G from
    the SAME diagonalization, and search all 24 index permutations. Only one
    spin species is needed (the Hamiltonian does not couple spins and A acts
    inside one Fock space), so the trace is 2**6 = 64 dimensional and exact.
    """

    NMODE = PH_NSITE * PH_NORB

    @classmethod
    def _fock(cls):
        dim = 1 << cls.NMODE
        c = []
        for p in range(cls.NMODE):
            op = np.zeros((dim, dim))
            for m in range(dim):
                if (m >> p) & 1:
                    sg = -1.0 if bin(m & ((1 << p) - 1)).count("1") % 2 else 1.0
                    op[m & ~(1 << p), m] = sg
            c.append(op)
        return c

    @classmethod
    def _run(cls, ntau):
        mode = lambda i, a: i * PH_NORB + a
        C = cls._fock()
        CD = [op.T.conj() for op in C]
        dim = 1 << cls.NMODE

        h1 = np.zeros((cls.NMODE, cls.NMODE), dtype=complex)
        for i in range(PH_NSITE):
            for a in range(PH_NORB):
                h1[mode(i, a), mode(i, a)] += PH_ONSITE[a] - PH_MU
                for d, amp in ((+1, PH_T_INTRA[a]),
                               (-1, np.conj(PH_T_INTRA[a]))):
                    h1[mode((i + d) % PH_NSITE, a), mode(i, a)] += amp
            h1[mode(i, 1), mode(i, 0)] += PH_T_INTER
            h1[mode(i, 0), mode(i, 1)] += np.conj(PH_T_INTER)
        assert np.allclose(h1, h1.conj().T)

        H = np.zeros((dim, dim), dtype=complex)
        for p in range(cls.NMODE):
            for q in range(cls.NMODE):
                if h1[p, q] != 0:
                    H += h1[p, q] * (CD[p] @ C[q])
        E, U = np.linalg.eigh(H)
        E -= E.min()
        w = np.exp(-PH_BETA * E)
        Z = w.sum()
        rot = lambda op: U.conj().T @ op @ U

        taus = np.arange(ntau) * PH_BETA / ntau
        expo = np.exp(np.multiply.outer(taus, E[:, None] - E[None, :]))

        def corr(A, B):
            """<A(tau) B(0)> on the tau grid."""
            return np.einsum("m,tmn,mn,nm->t", w, expo, A, B) / Z

        cr = [rot(op) for op in C]
        cdr = [rot(op) for op in CD]

        # G_{ab}(r, tau) = -<T c_{r a}(tau) c+_{0 b}(0)>
        G = np.zeros((ntau, PH_NSITE, PH_NORB, PH_NORB), dtype=complex)
        for r in range(PH_NSITE):
            for a in range(PH_NORB):
                for b in range(PH_NORB):
                    G[:, r, a, b] = -corr(cr[mode(r, a)], cdr[mode(0, b)])

        def bilinear(q, a, b):
            op = np.zeros((dim, dim), dtype=complex)
            for i in range(PH_NSITE):
                op += np.exp(-1j * q * i) * (CD[mode(i, a)] @ C[mode(i, b)])
            return op

        qs = 2.0 * np.pi * np.arange(PH_NSITE) / PH_NSITE
        ed = np.zeros((PH_NSITE,) + (PH_NORB,) * 4, dtype=complex)
        for iq, q in enumerate(qs):
            pl = {(a, b): rot(bilinear(q, a, b))
                  for a in range(PH_NORB) for b in range(PH_NORB)}
            mi = {(a, b): rot(bilinear(-q, a, b))
                  for a in range(PH_NORB) for b in range(PH_NORB)}
            for a in range(PH_NORB):
                for b in range(PH_NORB):
                    for c in range(PH_NORB):
                        for d in range(PH_NORB):
                            ed[iq, a, b, c, d] = corr(
                                pl[(a, b)], mi[(c, d)]).sum() * PH_BETA / ntau

        sgn = np.full(ntau, -1.0)
        sgn[0] = 1.0
        Grev = G[np.ix_((-np.arange(ntau)) % ntau,
                        (-np.arange(PH_NSITE)) % PH_NSITE)]
        rt = (G[:, :, :, None, :, None]
              * sgn[:, None, None, None, None, None]
              * Grev[:, :, None, :, None, :])           # [l,r,a,d,b,c]
        rt = rt.transpose(0, 1, 2, 5, 4, 3)             # -> [l,r,a,c,b,d]
        phase = np.exp(-1j * np.multiply.outer(qs, np.arange(PH_NSITE)))
        code = -np.einsum("qr,lracbd->qacbd", phase, rt) * (PH_BETA / ntau)
        return cls._ranked(ed, code)

    @classmethod
    def setUpClass(cls):
        cls.coarse = cls._run(PH_NTAU)
        cls.fine = cls._run(2 * PH_NTAU)

    @staticmethod
    def _ranked(ed, code):
        out = []
        for perm in itertools.permutations(range(4)):
            # q = 0 carries the disconnected part <A><A>, which the code's
            # bubble kernel does not contain; it cannot be compared.
            cand = np.transpose(code, (0,) + tuple(1 + p for p in perm))[1:]
            ref = ed[1:]
            den = np.vdot(cand, cand)
            if abs(den) < 1e-30:
                continue
            s = np.vdot(cand, ref) / den
            out.append((np.max(np.abs(ref - s * cand)) / np.max(np.abs(ref)),
                        perm))
        out.sort()
        return out

    def test_row_pair_is_reversed_and_column_pair_is_not(self):
        ranked = self.coarse
        best, runner_up = ranked[0], ranked[1]
        # code[a,c,b,d] ~ ED[c,a,b,d] = <A_{ca}(q) A_{bd}(-q)>
        self.assertEqual(best[1], (1, 0, 2, 3),
                         "chi0's pair index does not mean what #100 measured")
        self.assertLess(best[0], 1e-2)
        self.assertGreater(runner_up[0], 20 * best[0],
                           "the fixture does not separate the permutations")

    def test_the_residual_is_discretisation_not_disagreement(self):
        """Halving the tau spacing halves the winner's residual, so it is the
        grid and not a mismatch. Without this the threshold above would be a
        tolerance picked to make the test pass."""
        coarse, fine = self.coarse[0], self.fine[0]
        self.assertEqual(coarse[1], fine[1])
        self.assertLess(fine[0], 0.6 * coarse[0])
        self.assertGreater(fine[0], 0.4 * coarse[0])


class TestPairHopVertexPlacement(unittest.TestCase):
    """H_PH = sum_{aa'} P_{aa'} A^up_{aa'} A^down_{aa'} is a product of an UP
    and a DOWN bilinear, so the opposite-spin construction of
    TestRPAVertexOrientation applies once the densities are generalised to
    orbital-off-diagonal ones: the cross-spin response vanishes at P = 0 and
    its leading term is exactly one bubble-vertex-bubble chain,

        chi^ud(q) = -P X0(q) . Gamma . X0(q) + O(P^2),

    in the orbital-PAIR basis. Inverting it gives Gamma directly.

    P must be Hermitian-closed (P_{a'a} = conj(P_{aa'})) for H to be Hermitian;
    P is then Hermitian as a matrix, so P^T = conj(P) and the transpose is
    still a different matrix whenever P is complex. That is what makes the
    fixture discriminating.
    """

    P = {(0, 1): 1.0 + 0.3j, (1, 0): 1.0 - 0.3j}
    NPAIR = PH_NORB * PH_NORB

    @staticmethod
    def _measure(amp, q):
        """(X0, chi^ud) in the pair basis, at PairHop strength `amp`."""
        mode = lambda i, a, s: (i * PH_NORB + a) * 2 + s
        H = np.zeros((DIM, DIM), dtype=complex)
        for idx, mask in enumerate(STATES):
            for s in range(2):
                for i in range(NSITE):
                    for a in range(NORB):
                        p = mode(i, a, s)
                        if (mask >> p) & 1:
                            H[idx, idx] += ONSITE[a]
                        for d, t in ((+1, T_INTRA[a]),
                                     (-1, np.conj(T_INTRA[a]))):
                            r = _hop(mask, mode((i + d) % NSITE, a, s), p)
                            if r:
                                H[INDEX[r[0]], idx] += t * r[1]
                        t = T_INTER if a == 0 else np.conj(T_INTER)
                        r = _hop(mask, mode(i, 1 - a, s), p)
                        if r:
                            H[INDEX[r[0]], idx] += t * r[1]
        if amp:
            for (a, ap), v in TestPairHopVertexPlacement.P.items():
                for i in range(NSITE):
                    up = np.zeros((DIM, DIM), dtype=complex)
                    dn = np.zeros((DIM, DIM), dtype=complex)
                    for idx, mask in enumerate(STATES):
                        r = _hop(mask, mode(i, a, 0), mode(i, ap, 0))
                        if r:
                            up[INDEX[r[0]], idx] += r[1]
                        r = _hop(mask, mode(i, a, 1), mode(i, ap, 1))
                        if r:
                            dn[INDEX[r[0]], idx] += r[1]
                    # up and dn commute (different species, even parity), so
                    # up @ dn is the operator itself -- symmetrising would put
                    # a spurious factor of 2 into the measured coefficient.
                    H += amp * v * (up @ dn)
        assert np.allclose(H, H.conj().T), "H_PH is Hermitian only if P is"
        E, U = np.linalg.eigh(H)

        def bil(qq, a, b, s):
            op = np.zeros((DIM, DIM), dtype=complex)
            for idx, mask in enumerate(STATES):
                for i in range(NSITE):
                    r = _hop(mask, mode(i, a, s), mode(i, b, s))
                    if r:
                        op[INDEX[r[0]], idx] += np.exp(-1j * qq * i) * r[1]
            return op

        n = TestPairHopVertexPlacement.NPAIR
        uu = np.zeros((n, n), dtype=complex)
        ud = np.zeros((n, n), dtype=complex)
        pl = {(a, b, s): bil(q, a, b, s)
              for a in range(NORB) for b in range(NORB) for s in range(2)}
        mi = {(a, b, s): bil(-q, a, b, s)
              for a in range(NORB) for b in range(NORB) for s in range(2)}
        for a in range(NORB):
            for b in range(NORB):
                for c in range(NORB):
                    for d in range(NORB):
                        i, j = a * NORB + b, c * NORB + d
                        uu[i, j] = _static_chi(E, U, pl[(a, b, 0)],
                                               mi[(c, d, 0)])
                        ud[i, j] = _static_chi(E, U, pl[(a, b, 0)],
                                               mi[(c, d, 1)])
        return uu, ud

    @classmethod
    def setUpClass(cls):
        q = 2.0 * np.pi / NSITE
        cls.X0, cls.ud0 = cls._measure(0.0, q)
        # CENTRAL difference: the O(P^2) term cancels, so the remaining
        # truncation is O(dP^2) rather than the O(dP) of a forward difference
        # (which leaves ~1e-5 here). test_the_truncation_is_quadratic below
        # measures that scaling instead of asserting it.
        inv = np.linalg.inv(cls.X0)
        cls.cond_X0 = float(np.linalg.cond(cls.X0))

        def gamma_at(dP):
            _, up = cls._measure(dP, q)
            _, dn = cls._measure(-dP, q)
            # chi^ud = -X0 Gamma X0  =>  Gamma = -X0^-1 chi^ud X0^-1
            return -inv @ ((up - dn) / (2.0 * dP)) @ inv * NSITE

        cls.gamma = gamma_at(1.0e-4)
        cls.gamma_coarse = gamma_at(2.0e-4)

    def test_cross_spin_response_vanishes_without_pairhop(self):
        """The premise: at P = 0 the leading term is the whole story."""
        self.assertLess(float(np.max(np.abs(self.ud0))), 1e-12)

    def test_vertex_is_diagonal_in_the_bilinear_pair(self):
        """Gamma connects A^up_{ab} to A^down_{ab} -- the SAME pair -- so it is
        diagonal in the bilinear pair index, not antidiagonal."""
        off = self.gamma - np.diag(np.diag(self.gamma))
        # the floor is the O(dP^2) remainder of the central difference and the
        # conditioning of the X0 inversion (cond(X0) is checked below)
        self.assertLess(float(np.max(np.abs(off))),
                        1e-7 * float(np.max(np.abs(self.gamma))))

    def test_the_truncation_is_quadratic_and_x0_is_well_conditioned(self):
        """Halving dP cuts the off-diagonal leakage by about four, so the
        tolerances above are set by a measured truncation and not chosen to
        make the test pass."""
        fine = float(np.max(np.abs(self.gamma - np.diag(np.diag(self.gamma)))))
        coarse = float(np.max(np.abs(
            self.gamma_coarse - np.diag(np.diag(self.gamma_coarse)))))
        self.assertGreater(coarse / fine, 3.0)
        self.assertLess(coarse / fine, 5.0)
        self.assertLess(self.cond_X0, 1.0e3)

    def test_the_diagonal_carries_p_untransposed_with_unit_coefficient(self):
        """Gamma[(a,b),(a,b)] = P[a,b] exactly, sign and normalisation
        included -- not merely proportional. The complex fixture is what
        separates P from P^T = conj(P)."""
        got = np.array([self.gamma[a * NORB + b, a * NORB + b]
                        for (a, b) in self.P])
        want = np.array([self.P[(a, b)] for (a, b) in self.P])
        np.testing.assert_allclose(got, want, rtol=0.0,
                                   atol=1e-7 * float(np.max(np.abs(want))))
        # and P^T is genuinely excluded
        wt = np.array([self.P[(b, a)] for (a, b) in self.P])
        self.assertGreater(np.max(np.abs(got - wt)),
                           1e-2 * float(np.max(np.abs(got))))



class TestPairHopEndToEnd(unittest.TestCase):
    """The determination that actually settles #100.

    No index reasoning anywhere between the input file and the answer. The
    chain under test is the one the solver runs:

        interaction dict -> _build_interaction_k -> _build_sc_matrices_all_q
                                                                     -> S
        chi0 <- rpa.py's own kernel, fed with an EXACT Green function
        prediction:  chi^ud = -chi0 . S . chi0        (linear order)

    compared against the cross-spin response of the SAME Hamiltonian by exact
    diagonalization. It is run for both orientations of the input matrix, and
    whichever the chain needs in order to reproduce the exact answer is what
    ``_to_k`` must store.

    Why linear order is exactly one chain: H_PH couples an UP bilinear to a
    DOWN one, so the cross-spin response vanishes identically at P = 0 and its
    derivative is a single bubble-vertex-bubble term.

    Why chi^ud: Case 4 sets S = C, and with the MYO signs
    chi_s = [I - chi0 Us]^-1 chi0 and chi_c = [I + chi0 Uc]^-1 chi0, so
    chi^ud = (chi_c - chi_s)/2 = -chi0 . S . chi0 + O(P^2).

    Grand canonical throughout, because rpa.py's chi0 is. Every operator here
    except the Green function conserves (N_up, N_down), so the trace splits
    into sectors of dimension at most C(6,3)^2 = 400.
    """

    NTAU = 192
    MU = 0.2
    DP = 1.0e-4
    # Hermitian-closed so that H_PH is Hermitian. P is then Hermitian as a
    # matrix, so P^T = conj(P) is a DIFFERENT matrix and the fixture
    # discriminates -- a real asymmetric P would make H non-Hermitian.
    P_IN = {(0, 1): 1.0 + 0.3j, (1, 0): 1.0 - 0.3j}

    @classmethod
    def _h1(cls):
        n = NSITE * NORB
        m = lambda i, a: i * NORB + a
        h = np.zeros((n, n), dtype=complex)
        for i in range(NSITE):
            for a in range(NORB):
                h[m(i, a), m(i, a)] += ONSITE[a] - cls.MU
                for d, t in ((+1, T_INTRA[a]), (-1, np.conj(T_INTRA[a]))):
                    h[m((i + d) % NSITE, a), m(i, a)] += t
            h[m(i, 1), m(i, 0)] += T_INTER
            h[m(i, 0), m(i, 1)] += np.conj(T_INTER)
        assert np.allclose(h, h.conj().T)
        return h

    @staticmethod
    def _cdc(mask, p, q):
        if not (mask >> q) & 1:
            return None
        sg = -1.0 if bin(mask & ((1 << q) - 1)).count("1") % 2 else 1.0
        mm = mask & ~(1 << q)
        if (mm >> p) & 1:
            return None
        sg *= -1.0 if bin(mm & ((1 << p) - 1)).count("1") % 2 else 1.0
        return mm | (1 << p), sg

    @classmethod
    def _green(cls):
        """G_{ab}(r, tau) = -<T c_{r a}(tau) c+_{0 b}(0)>.

        One spin species suffices: at P = 0 the grand-canonical density matrix
        factorises over spin.
        """
        nmode = NSITE * NORB
        dim = 1 << nmode
        C = []
        for p in range(nmode):
            op = np.zeros((dim, dim))
            for st in range(dim):
                if (st >> p) & 1:
                    sg = (-1.0 if bin(st & ((1 << p) - 1)).count("1") % 2
                          else 1.0)
                    op[st & ~(1 << p), st] = sg
            C.append(op)
        CD = [c.T.conj() for c in C]
        h = cls._h1()
        H = np.zeros((dim, dim), dtype=complex)
        for p in range(nmode):
            for q in range(nmode):
                if h[p, q]:
                    H += h[p, q] * (CD[p] @ C[q])
        E, U = np.linalg.eigh(H)
        E -= E.min()
        w = np.exp(-BETA * E)
        taus = np.arange(cls.NTAU) * BETA / cls.NTAU
        expo = np.exp(np.multiply.outer(taus, E[:, None] - E[None, :]))
        m = lambda i, a: i * NORB + a
        G = np.zeros((cls.NTAU, NSITE, NORB, NORB), dtype=complex)
        for r in range(NSITE):
            for a in range(NORB):
                for b in range(NORB):
                    A = U.conj().T @ C[m(r, a)] @ U
                    B = U.conj().T @ CD[m(0, b)] @ U
                    G[:, r, a, b] = -np.einsum(
                        "m,tmn,mn,nm->t", w, expo, A, B) / w.sum()
        return G

    @classmethod
    def _chi0_code(cls, G, qs):
        """rpa.py::_calc_chi0q's general branch at zero bosonic frequency:
        chi0[a,c,b,d] = G[a,b] * G_rev[d,c] * sgn, G_rev(l,r) = G(-l, -r)."""
        nt = cls.NTAU
        sgn = np.full(nt, -1.0)
        sgn[0] = 1.0
        Grev = G[np.ix_((-np.arange(nt)) % nt, (-np.arange(NSITE)) % NSITE)]
        rt = (G[:, :, :, None, :, None]
              * sgn[:, None, None, None, None, None]
              * Grev[:, :, None, :, None, :])         # [l,r,a,d,b,c]
        rt = rt.transpose(0, 1, 2, 5, 4, 3)           # -> [l,r,a,c,b,d]
        phase = np.exp(-1j * np.multiply.outer(qs, np.arange(NSITE)))
        out = -np.einsum("qr,lracbd->qacbd", phase, rt) * (BETA / nt)
        return out.reshape(NSITE, NORB ** 2, NORB ** 2)

    @classmethod
    def _response(cls, amp, qs, same_spin=False):
        """int dtau <A^u_{ab}(q,tau) A^s_{cd}(-q,0)>, grand canonical."""
        nm = NSITE * NORB
        mode = lambda i, a, s: (i * NORB + a) * 2 + s
        h = cls._h1()
        npair = NORB * NORB
        num = np.zeros((NSITE, npair, npair), dtype=complex)
        Z = 0.0
        ups = [mode(i, a, 0) for i in range(NSITE) for a in range(NORB)]
        dns = [mode(i, a, 1) for i in range(NSITE) for a in range(NORB)]
        for nu in range(nm + 1):
            for nd in range(nm + 1):
                states = []
                for su in itertools.combinations(ups, nu):
                    for sd in itertools.combinations(dns, nd):
                        m = 0
                        for p in su + sd:
                            m |= 1 << p
                        states.append(m)
                states.sort()
                idx = {st: i for i, st in enumerate(states)}
                D = len(states)

                def bil(q, a, b, s):
                    op = np.zeros((D, D), dtype=complex)
                    for j, mask in enumerate(states):
                        for i in range(NSITE):
                            r = cls._cdc(mask, mode(i, a, s), mode(i, b, s))
                            if r:
                                op[idx[r[0]], j] += np.exp(-1j * q * i) * r[1]
                    return op

                H = np.zeros((D, D), dtype=complex)
                for s in range(2):
                    for p in range(nm):
                        for q in range(nm):
                            if not h[p, q]:
                                continue
                            for j, mask in enumerate(states):
                                r = cls._cdc(mask,
                                             mode(p // NORB, p % NORB, s),
                                             mode(q // NORB, q % NORB, s))
                                if r:
                                    H[idx[r[0]], j] += h[p, q] * r[1]
                if amp:
                    # H_PH = sum_{aa'} P_{aa'} A^u_{aa'} A^d_{aa'}; the two
                    # bilinears commute, so no symmetrisation and no stray 2.
                    for (a, ap), v in cls.P_IN.items():
                        for i in range(NSITE):
                            u = np.zeros((D, D), dtype=complex)
                            w_ = np.zeros((D, D), dtype=complex)
                            for j, mask in enumerate(states):
                                r = cls._cdc(mask, mode(i, a, 0),
                                             mode(i, ap, 0))
                                if r:
                                    u[idx[r[0]], j] += r[1]
                                r = cls._cdc(mask, mode(i, a, 1),
                                             mode(i, ap, 1))
                                if r:
                                    w_[idx[r[0]], j] += r[1]
                            H += amp * v * (u @ w_)
                assert np.allclose(H, H.conj().T), "H must be Hermitian"
                E, U = np.linalg.eigh(H)
                # _h1 already carries -MU on the diagonal, so H is H0 - MU N
                # and NO further fugacity factor may be applied here.
                w = np.exp(-BETA * E)
                dE = E[:, None] - E[None, :]
                with np.errstate(divide="ignore", invalid="ignore"):
                    kern = (w[None, :] - w[:, None]) / dE
                deg = np.abs(dE) < 1e-10
                kern[deg] = (BETA * w[:, None] * np.ones_like(dE))[deg]
                sb = 0 if same_spin else 1
                for iq, q in enumerate(qs):
                    A = [U.conj().T @ bil(q, a, b, 0) @ U
                         for a in range(NORB) for b in range(NORB)]
                    B = [U.conj().T @ bil(-q, a, b, sb) @ U
                         for a in range(NORB) for b in range(NORB)]
                    for i, Ar in enumerate(A):
                        for j, Br in enumerate(B):
                            num[iq, i, j] += np.sum(Ar * Br.T * kern)
                Z += w.sum()
        return num / Z

    @staticmethod
    def _reverse_row_pair(arr):
        rev = np.array([(x % NORB) * NORB + (x // NORB)
                        for x in range(NORB ** 2)])
        return arr[:, rev, :]

    @classmethod
    def setUpClass(cls):
        cls.qs = 2.0 * np.pi * np.arange(NSITE) / NSITE
        cls.chi0 = cls._chi0_code(cls._green(), cls.qs)
        cls.bubble = cls._reverse_row_pair(cls._response(0.0, cls.qs,
                                                         same_spin=True))
        cls.ud0 = cls._response(0.0, cls.qs)
        cls.lead = cls._reverse_row_pair(
            (cls._response(+cls.DP, cls.qs) - cls._response(-cls.DP, cls.qs))
            / (2.0 * cls.DP))

    @staticmethod
    def _worst(pred, ref):
        out = 0.0
        for iq in range(1, NSITE):     # q = 0 carries a disconnected part
            s = np.vdot(pred[iq], ref[iq]) / np.vdot(pred[iq], pred[iq])
            out = max(out, np.max(np.abs(ref[iq] - s * pred[iq]))
                      / np.max(np.abs(ref[iq])))
        return out

    def test_cross_spin_response_vanishes_without_pairhop(self):
        """The premise: at P = 0 the leading term is the whole story. q = 0 is
        excluded because it carries the disconnected part <A^u><A^d>."""
        self.assertLess(float(np.max(np.abs(self.ud0[1:]))), 1e-12)

    def test_the_code_kernel_reproduces_the_exact_bubble(self):
        """Locate the residual before using it: chi0 built from rpa.py's own
        kernel matches the exact same-spin bubble to the imaginary-time
        discretisation, so any larger residual downstream is a real
        disagreement and not the grid."""
        self.assertLess(self._worst(self.chi0, self.bubble), 1e-2)

    def _sc_matrix(self, mat):
        import hwave.sc as sc

        k = np.array([0.0])
        vr = {((0, 0, 0), ab): v for ab, v in mat.items()}
        ik = sc._build_interaction_k(k, k, k, {"PairHop": vr}, NORB)
        S, C = sc._build_sc_matrices_all_q(ik, NORB, 1, 1, 1)
        np.testing.assert_allclose(S, C, rtol=0.0, atol=1e-14)
        return S[0, 0, 0]

    def test_the_untransposed_placement_reproduces_exact_diagonalization(self):
        P = self.P_IN
        PT = {(b, a): v for (a, b), v in P.items()}
        keep = self._sc_matrix(P)     # _to_k leaves PairHop alone -> S = P
        flip = self._sc_matrix(PT)

        # say out loud which placement each one is, so the test cannot silently
        # invert if `_to_k` changes
        self.assertAlmostEqual(complex(keep[0 * NORB + 1, 1 * NORB + 0]),
                               P[(0, 1)], places=12)
        self.assertAlmostEqual(complex(flip[0 * NORB + 1, 1 * NORB + 0]),
                               P[(1, 0)], places=12)

        r_keep = self._worst(
            np.array([-self.chi0[i] @ keep @ self.chi0[i]
                      for i in range(NSITE)]), self.lead)
        r_flip = self._worst(
            np.array([-self.chi0[i] @ flip @ self.chi0[i]
                      for i in range(NSITE)]), self.lead)
        self.assertLess(r_keep, 1e-2,
                        "S[(a,b),(b,a)] = P[a,b] must reproduce the exact "
                        "cross-spin response")
        self.assertGreater(r_flip, 1e-1,
                           "the transposed placement must be excluded")


if __name__ == "__main__":
    unittest.main()
