"""The bond-resolved PAIRING kernels put a declared two-body row on the bond
of the Transfer row with the same ``(R, a, b)`` (issue #232).

``tests/test_two_body_placement.py`` pins this for the particle-hole paths.
Here the same question is asked of the particle-particle (Eliashberg) use of
the bond channels: the bare Cooper block and the fluctuation term of the
static kernel (``bond_channels.make_bond_kernel_parts``) and the dynamic
kernel (``eliashberg_bond.BondPairKernel``, with its bare term
``instantaneous_vertex``).

The reference contains no H-wave code. The ``NX``-site chain with ``NORB``
orbitals is treated as ONE cell with ``NX * NORB`` orbitals, so that every
interaction is "local"; the pairing vertex is built in the spin-orbital basis
from the antisymmetrized bare interaction ``A`` (first order) and its two
static particle-hole second-order diagrams (direct and crossed) evaluated with
the exact cluster Green function on the same finite Matsubara grid. The gap
equation is ``Delta_pq = 1/2 sum_rs Gamma_pq,rs kappa_sr`` with the pair
amplitude ``kappa_pq = <c_p c_q>``; the kernels act on the ``(up, down)``
block.

The model is chosen so that the two placements are different Hamiltonians:
an inter-orbital CoulombInter table that is not symmetric under ``R -> -R``
and a transpose-asymmetric hopping. Each check is two-sided: the production
kernel matches the documented placement to round-off and misses the opposite
placement by O(1). The production second-order term is made exactly second
order by dressing with the BARE bond bubble (``chi_s = chi_c = chi_bar``), so
it uses the bubble the production code builds from the same Green function.
"""
import types
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from tests.eliashberg_bond_fixtures import _channel0_from_gate, _tables_from_pairs

NX, NORB, BETA, NMAT, MU = 6, 2, 2.0, 64, 0.15
ND = NORB * NORB
#: Transfer rows (r, a, b): t c+_{j,a} c_{j+r,b}; H(k) = sum_r e^{+ikr} t(r).
#: The inter-orbital bond (0 at the origin, 1 at +x) is strong, its mirror
#: (0 at the origin, 1 at -x) carries no hopping.
TRANSFER = {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
            (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5}
#: CoulombInter rows (r, a, b): each declared with its partner (-r, b, a), so
#: the bond (0 at the origin, 1 at +x) carries 0.4 and its mirror 0.1.
ROWS = {(1, 0, 1): 0.4, (-1, 1, 0): 0.4, (-1, 0, 1): 0.1, (1, 1, 0): 0.1}

_K = 2 * np.pi * np.arange(NX) / NX
_WN = (2 * np.arange(NMAT) + 1 - NMAT) * np.pi / BETA


def _green_sc():
    """Bare G(k, i w_n) in the sc.py layout (norb, norb, Nx, 1, 1, nmat)."""
    h = np.zeros((NX, NORB, NORB), complex)
    for (r, a, b), t in TRANSFER.items():
        h[:, a, b] += t * np.exp(1j * _K * r)
    g = np.empty((NORB, NORB, NX, 1, 1, NMAT), complex)
    for n, w in enumerate(_WN):
        g[:, :, :, 0, 0, n] = np.linalg.inv(
            (1j * w + MU) * np.eye(NORB)[None] - h).transpose(1, 2, 0)
    return g


# ---------------------------------------------------------------------------
# independent reference: the chain as one cell, spin-orbital vertex
# ---------------------------------------------------------------------------

def _cluster(sign):
    """Hopping matrix and density-density matrix ``W`` (``H_int = 1/2
    sum_{I != J} W_IJ n_I n_J``, spin summed) of the chain. ``sign = +1`` puts
    a row ``(r, a, b)`` on the bond of the Transfer row ``(r, a, b)`` (``a``
    at ``j``, ``b`` at ``j + r``); ``sign = -1`` on the opposite bond. As in
    the solvers, a row and its partner name the same density operator and
    each carries half of it."""
    M = NX * NORB

    def idx(j, a):
        return (j % NX) * NORB + a

    h = np.zeros((M, M), complex)
    for (r, a, b), t in TRANSFER.items():
        for j in range(NX):
            h[idx(j, a), idx(j + r, b)] += t
    W = np.zeros((M, M))
    for (r, a, b), v in ROWS.items():
        for j in range(NX):
            p, q = idx(j, a), idx(j + sign * r, b)
            W[p, q] += v / 2
            W[q, p] += v / 2
    return h, W


def _reference_vertices(sign):
    """``(Gamma1, Gamma2)`` as spin-orbital tensors ``[p, q, r, s]`` (out
    ``p, q``; in ``r, s``): the antisymmetrized bare interaction and the
    static second-order particle-hole vertex (direct minus its ``r <-> s``
    crossed partner)."""
    h, W = _cluster(sign)
    M = h.shape[0]
    P2 = 2 * M                                    # spin-orbital p = 2 I + s
    U = np.zeros((P2, P2, P2, P2))
    for I in range(M):
        for J in range(M):
            if W[I, J] != 0:
                for s in range(2):
                    for t in range(2):
                        U[2 * I + s, 2 * J + t, 2 * I + s, 2 * J + t] += W[I, J]
    A = U - U.transpose(0, 1, 3, 2)
    # L[6, 7, 5, 8] = T sum_n G_67(i w_n) G_85(i w_n)
    L = np.zeros((P2,) * 4, complex)
    for w in _WN:
        G = np.kron(np.linalg.inv((1j * w + MU) * np.eye(M) - h), np.eye(2))
        L += np.einsum("ab,cd->abdc", G, G) / BETA
    # direct: Gd[1,2,3,4] = sum A[1,5,3,6] L[6,7,5,8] A[7,2,8,4]
    left = A.transpose(0, 2, 3, 1).reshape(P2 * P2, P2 * P2)
    mid = L.transpose(0, 2, 1, 3).reshape(P2 * P2, P2 * P2)
    right = A.transpose(0, 2, 1, 3).reshape(P2 * P2, P2 * P2)
    Gd = (left @ mid @ right).reshape(P2, P2, P2, P2).transpose(0, 2, 1, 3)
    return A, Gd - Gd.transpose(0, 1, 3, 2)


def _real_space_pairs(X):
    """``kappa_{(j,a),(j',b)} = F_ab(j' - j) = (1/N) sum_k e^{-ik(j'-j)} X_ab(k)``
    for ``X_ab(k) = <c_{k a up} c_{-k b dn}>``."""
    M = NX * NORB
    kap = np.zeros((M, M), complex)
    for j in range(NX):
        for jp in range(NX):
            kap[j * NORB:(j + 1) * NORB, jp * NORB:(jp + 1) * NORB] = (
                X * np.exp(-1j * _K * (jp - j))[None, None]).mean(-1)
    return kap


def _reference_gap_real_space(Gam, X):
    """``Delta_{(0,a),(d,b)}`` (up, down) for the pair field ``X``:
    ``Delta_pq = 1/2 sum_rs Gamma_pq,rs kappa_sr``."""
    M = NX * NORB
    P2 = 2 * M
    kap = _real_space_pairs(X)
    kso = np.zeros((P2, P2), complex)
    kso[0::2, 1::2] = kap               # <c_{I up} c_{J dn}>
    kso[1::2, 0::2] = -kap.T            # <c_{J dn} c_{I up}> = -<c_{I up} c_{J dn}>
    D = 0.5 * np.einsum("pqrs,sr->pq", Gam, kso)[0::2, 1::2]
    return np.stack([D[0:NORB, d * NORB:(d + 1) * NORB] for d in range(NX)], -1)


def _to_k(Dr):
    """``Delta_ab(k) = sum_d e^{+ikd} Delta_ab(d)``, the Fourier convention of
    the pairing term ``sum_k Delta_ab(k) c+_{k a up} c+_{-k b dn}``."""
    return np.einsum("abd,kd->abk", Dr, np.exp(1j * np.outer(_K, np.arange(NX))))


# ---------------------------------------------------------------------------
# production inputs
# ---------------------------------------------------------------------------

def _production_inputs():
    inter = {}
    for (r, a, b), v in ROWS.items():
        inter.setdefault(("CoulombInter", (r, 0, 0)), np.zeros((NORB, NORB)))[a, b] = v
    whole, onsite, offsite = _tables_from_pairs(inter, 0.0, NORB)
    S0, C0 = _channel0_from_gate(whole, onsite, offsite, NORB, (NX, 1, 1))
    bond_set = bc.resolve_interactions(
        {((r, 0, 0), (a, b)): v for (r, a, b), v in ROWS.items()}, np.eye(3), NORB)
    topo = bc.resolve_bond_topology(
        {"CoulombInter": {((r, 0, 0), (a, b)): complex(v) for (r, a, b), v in ROWS.items()}},
        np.eye(3), NORB, active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
    return bond_set, topo, S0, C0


def _pair_field(rng, eta, nfreq=None):
    """A random gap of definite spin parity: ``phi_ab(k) = +-phi_ba(-k)``."""
    shape = (NORB, NORB, NX) + (() if nfreq is None else (nfreq,))
    phi = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    sgn = 1.0 if eta == "singlet" else -1.0
    mk = (-np.arange(NX)) % NX
    if nfreq is None:
        return 0.5 * (phi + sgn * phi[:, :, mk].transpose(1, 0, 2))
    return 0.5 * (phi + sgn * phi[:, :, mk, ::-1].transpose(1, 0, 2, 3))


def _static_x(g, phi):
    """``X_cd(k) = T sum_n G_ce(k, w_n) G_df(-k, -w_n) phi_ef(k)``."""
    gm = np.roll(g[:, :, ::-1, :, :, ::-1], 1, axis=2)[:, :, :, 0, 0]   # G(-k, -w)
    return np.einsum("cekn,dfkn,efk->cdk", g[:, :, :, 0, 0], gm, phi) / BETA


def _dynamic_x(g, phi_w):
    """The frequency sum of the dynamic pair field ``F(k, n) = G2(k, n) phi(k, n)``."""
    gm = np.roll(g[:, :, ::-1, :, :, ::-1], 1, axis=2)[:, :, :, 0, 0]
    return np.einsum("cekn,dfkn,efkn->cdk", g[:, :, :, 0, 0], gm, phi_w) / BETA


class _PlacementCase(object):
    """Shared reference data; subclasses compare one production kernel."""

    @classmethod
    def setUpClass(cls):
        cls.g = _green_sc()
        cls.ref = {s: _reference_vertices(s) for s in (+1, -1)}
        cls.bond_set, cls.topo, cls.S0, cls.C0 = _production_inputs()

    def _assert_two_sided(self, got, X, order, label):
        def deviation(sign):
            want = _to_k(_reference_gap_real_space(self.ref[sign][order - 1], X))
            return np.linalg.norm(got - want) / np.linalg.norm(want)
        same, opposite = deviation(+1), deviation(-1)
        self.assertLess(same, 1e-10, "{}: same-bond deviation {:.3e}".format(label, same))
        self.assertGreater(opposite, 1e-1,
                           "{}: opposite-bond deviation {:.3e}".format(label, opposite))


class TestReferenceSeesThePlacement(_PlacementCase, unittest.TestCase):
    def test_placements_differ_at_both_orders(self):
        X = _static_x(self.g, _pair_field(np.random.default_rng(0), "singlet"))
        for order in (1, 2):
            a = _reference_gap_real_space(self.ref[+1][order - 1], X)
            b = _reference_gap_real_space(self.ref[-1][order - 1], X)
            self.assertGreater(np.linalg.norm(a - b) / np.linalg.norm(a), 0.3)

    def test_first_order_is_the_bare_density_interaction(self):
        """``Delta_{I up, J dn} = -W_IJ kappa_IJ``: the reference's own
        normalization, derived by hand from ``H_int``."""
        X = _static_x(self.g, _pair_field(np.random.default_rng(1), "triplet"))
        _, W = _cluster(+1)
        want = -W * _real_space_pairs(X)
        got = _reference_gap_real_space(self.ref[+1][0], X)
        for d in range(NX):
            np.testing.assert_allclose(got[:, :, d], want[0:NORB, d * NORB:(d + 1) * NORB],
                                       atol=1e-14)


class TestStaticBondKernelPlacement(_PlacementCase, unittest.TestCase):
    """``make_bond_kernel_parts``: the bare Cooper block (``V^pp`` from
    ``bare_bond_vertices``) and the fluctuation term with the bare bubble."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        nq = (NX, 1, 1, ND, ND)
        cls.S_b, cls.C_b, cls.Vs, cls.Vt = bc.bare_bond_vertices(
            cls.bond_set, cls.S0.reshape(nq), cls.C0.reshape(nq), NORB)
        cls.chib = bc.bond_bubble(cls.g, cls.bond_set, BETA)

    def _parts(self, eta):
        _, A_fl, A_pp, _ = bc.make_bond_kernel_parts(
            self.chib, self.chib, self.S_b, self.C_b, self.Vs, self.Vt, self.g,
            self.bond_set, eta, BETA)
        return A_fl, A_pp

    def test_bare_cooper_block(self):
        for i, eta in enumerate(("singlet", "triplet")):
            phi = _pair_field(np.random.default_rng(10 + i), eta)
            _, A_pp = self._parts(eta)
            got = A_pp.matvec(phi.ravel()).reshape(NORB, NORB, NX)
            self._assert_two_sided(got, _static_x(self.g, phi), 1, "V^pp " + eta)

    def test_fluctuation_term(self):
        for i, eta in enumerate(("singlet", "triplet")):
            phi = _pair_field(np.random.default_rng(20 + i), eta)
            A_fl, _ = self._parts(eta)
            got = A_fl.matvec(phi.ravel()).reshape(NORB, NORB, NX)
            self._assert_two_sided(got, _static_x(self.g, phi), 2, "K^fl " + eta)


class TestDynamicBondKernelPlacement(_PlacementCase, unittest.TestCase):
    """``BondPairKernel`` with the FLEX gate's own ``S``/``C`` bond vertices.
    The pair-leg geometry does not depend on frequency, so the vertex is the
    static second-order one, flat in frequency; the kernel then maps a
    frequency-dependent gap to ``-(1/N) sum_k' Gamma(k, k') sum_n F(k', n)``."""

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.S = bc.build_sc_bond_channel(cls.topo, cls.S0, "S")
        cls.C = bc.build_sc_bond_channel(cls.topo, cls.C0, "C")
        cls.chib = bc.bond_bubble(cls.g, cls.bond_set, BETA).reshape(
            NX, cls.S.shape[-1], cls.S.shape[-1])

    def _kernel(self, gam, V_inst):
        from hwave.solver.eliashberg_bond import (ArrayBlockSource, BondPairKernel,
                                                  PairVertexUniform)
        from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
        B = int(self.topo.delta_r.shape[0])
        flat = np.ascontiguousarray(np.broadcast_to(gam[None], (NMAT,) + gam.shape))
        vert = PairVertexUniform(source=ArrayBlockSource({"G": flat}, ND), slot="G",
                                 B=B, nd=ND, nmat=NMAT, nvol=NX)
        return BondPairKernel(vert, calc_g2_dynamic(self.g, BETA), bc.BondSetView(self.topo),
                              xp=np, spatial_shape=(NX, 1, 1), norb=NORB, beta=BETA,
                              nfreq=NMAT, V_inst=V_inst, residency="host",
                              host_cap=10 ** 12, device_cap=10 ** 12)

    def _apply(self, K, phi_w):
        out = K.matvec(phi_w[:, :, :, None, None, :].ravel()).reshape(
            NORB, NORB, NX, NMAT)
        np.testing.assert_allclose(out, np.broadcast_to(out[..., :1], out.shape),
                                   atol=1e-12)            # flat vertex, flat output
        return out[..., 0]

    def test_fluctuation_term(self):
        from hwave.solver.eliashberg_bond import COEFF
        for i, eta in enumerate(("singlet", "triplet")):
            cs, cc = COEFF[eta]
            gam = cs * self.S @ self.chib @ self.S + cc * self.C @ self.chib @ self.C
            phi_w = _pair_field(np.random.default_rng(30 + i), eta, NMAT)
            got = self._apply(self._kernel(gam, None), phi_w)
            self._assert_two_sided(got, _dynamic_x(self.g, phi_w), 2, "Gamma " + eta)

    def test_bare_term(self):
        from hwave.solver.eliashberg_bond import instantaneous_vertex
        zero = np.zeros((NX,) + self.S.shape[1:], complex)
        for i, eta in enumerate(("singlet", "triplet")):
            V_inst = instantaneous_vertex(self.S, self.C, ND, eta, (NX, 1, 1))
            phi_w = _pair_field(np.random.default_rng(40 + i), eta, NMAT)
            got = self._apply(self._kernel(zero, V_inst), phi_w)
            self._assert_two_sided(got, _dynamic_x(self.g, phi_w), 1, "V_inst " + eta)


class TestGapBondProjection(_PlacementCase, unittest.TestCase):
    def test_psi_is_the_pair_at_plus_delta_r(self):
        """``gap_bond_projection``'s ``psi_m[a, b]`` is the real-space gap of
        the pair with ``a`` in the origin cell and ``b`` at ``+delta_r[m]``."""
        from hwave.solver.eliashberg_bond import gap_bond_projection
        X = _static_x(self.g, _pair_field(np.random.default_rng(50), "singlet"))
        Dr = _reference_gap_real_space(self.ref[+1][1], X)          # (a, b, d)
        view = bc.BondSetView(self.topo)
        psi = gap_bond_projection(_to_k(Dr)[:, :, :, None, None, None], view, (NX, 1, 1))
        for m, R in enumerate(view.delta_r):
            np.testing.assert_allclose(psi[m, :, :, 0], Dr[:, :, int(R[0]) % NX], atol=1e-13)


if __name__ == "__main__":
    unittest.main()
