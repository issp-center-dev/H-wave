"""The static ``hwave_sc`` bond pairing kernel at FIRST order in an off-site
coupling, against an independent spin-orbital construction (issue #239).

Since #239 the static bond path resolves every off-site ``CoulombInter``,
``Hund``, ``Ising``, ``Exchange`` and ``PairHop`` row on its bond channel.
Here the ``O(coupling)`` coefficient of the static pairing kernel built by
the production chain (``sc._resolve_bond_topology`` ->
``sc._build_bond_vertices`` -> ``bond_channels.dress_batch`` ->
``sc._bond_kernel_parts``, through ``sc._build_bond_operator`` where the
model allows it) is compared with the first-order pairing vertex of the SAME
Hamiltonian written out in the spin-orbital basis of a finite chain.

The reference contains no H-wave vertex code. The ``NX``-site chain is
treated as one cell with ``NX * norb`` orbitals; each interaction row is
expanded into its quartic operator terms ``v c+_p c_q c+_r c_s`` by the
term builders the exact-diagonalization campaigns use (the operator
conventions the bond-resolved longitudinal vertex was adjudicated against,
``tests/test_bond_longitudinal_ed.py``), normal ordered, antisymmetrized,
and contracted with the pair amplitude: ``Delta_pq = 1/2 sum_rs A_pq,rs
kappa_sr`` -- the construction of ``tests/test_bond_pairing_placement.py``,
generalized from density-density terms to arbitrary quartic ones. The singlet
is the ``(up, down)`` block, the triplet the equal-spin ``(up, up)`` block
(for a spin-rotation-invariant interaction the triplet components agree; for
an Sz-only coupling such as Hund or Ising alone they do not, and the
kernel's triplet is the equal-spin one). The term builders place a row
``(R, a, b)`` with orbital ``a`` on site ``j + R`` and ``b`` on site ``j``
(their fixture frame is the mirror of the file frame), so a file row
``(R, a, b)`` is expanded with the builder at ``-R``.

A declared row and its Hermitian partner ``(-R, b, a)`` name ONE operator
(the file reader requires both): the reference expands the row once.

Each check is run for both signs of the coupling and for both orientations
of the declared row (``+R`` and ``-R``), at ``norb = 1`` (the static path's
supported case) and on the inter-orbital ``norb = 2`` bond of the placement
test. The bare part ``K^pp`` must equal the reference at every coupling, and
the first-order coefficient of the FULL static kernel (bare plus dressed
fluctuation part, by a Richardson-extrapolated central difference) must equal
it too. Where the reference does not vanish identically it is asserted to be
non-zero (the check is not vacuous).
"""
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import bond_channels as bc
from tests.ed_oracle_util import EDFixture
from tests.test_bond_longitudinal_ed import _terms_pairhop_offsite
from tests.test_bond_transverse_ed import (_terms_exchange_offsite, _terms_hund_offsite,
                                           _terms_ising_offsite)

NX, BETA, NMAT, MU = 6, 2.0, 64, 0.15
TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop")
_K = 2 * np.pi * np.arange(NX) / NX
_WN = (2 * np.arange(NMAT) + 1 - NMAT) * np.pi / BETA

#: Transfer rows (r, a, b): t c+_{j,a} c_{j+r,b}. norb = 1: a plain chain;
#: norb = 2: the mirror- and transpose-asymmetric chain of the placement test.
TRANSFER = {
    1: {(1, 0, 0): -1.0, (-1, 0, 0): -1.0, (0, 0, 0): 0.1},
    2: {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
        (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5},
}
#: the declared row (a, b) per orbital count
ROW = {1: (0, 0), 2: (0, 1)}


def _terms(t, fx, a, b, R, v):
    """Quartic terms of the file row ``(R, a, b)`` of type ``t`` at value
    ``v``, in the file frame (``a`` on site ``j``, ``b`` on site ``j + R``)."""
    if t == "CoulombInter":
        terms = []
        for j in range(fx.L):
            for s1 in range(2):
                for s2 in range(2):
                    p, r = fx.mode(j, a, s1), fx.mode((j + R) % fx.L, b, s2)
                    terms.append((p, p, r, r, v))
        return terms
    builder = {"Hund": _terms_hund_offsite, "Ising": _terms_ising_offsite,
               "Exchange": _terms_exchange_offsite, "PairHop": _terms_pairhop_offsite}[t]
    return builder(fx, a, b, -R, v)


def _green_sc(norb):
    """Bare G(k, i w_n) in the sc.py layout (norb, norb, Nx, 1, 1, nmat)."""
    h = np.zeros((NX, norb, norb), complex)
    for (r, a, b), t in TRANSFER[norb].items():
        h[:, a, b] += t * np.exp(1j * _K * r)
    g = np.empty((norb, norb, NX, 1, 1, NMAT), complex)
    for n, w in enumerate(_WN):
        g[:, :, :, 0, 0, n] = np.linalg.inv(
            (1j * w + MU) * np.eye(norb)[None] - h).transpose(1, 2, 0)
    return g


def _antisymmetrized(terms, P2):
    """``A[p, q, r, s]`` (out ``p, q``; in ``r, s``) of the two-body part of
    ``sum v c+_p c_q c+_r c_s`` = ``sum v c+_p c+_r c_s c_q`` (+ one-body),
    written as ``1/2 sum U[a, b, c, d] c+_a c+_b c_d c_c`` with ``U`` symmetric
    under the exchange of the two particles, then antisymmetrized."""
    U = np.zeros((P2,) * 4, complex)
    for (p, q, r, s, v) in terms:
        U[p, r, q, s] += v
        U[r, p, s, q] += v
    return U - U.transpose(0, 1, 3, 2)


def _reference_gap(norb, A, X, eta):
    """``Delta_ab(k)`` of the first-order gap equation for the pair field
    ``X_ab(k)``: the ``(up, down)`` pair ``<c_{k a up} c_{-k b dn}>`` for the
    singlet, the equal-spin pair ``<c_{k a up} c_{-k b up}>`` for the triplet
    (see the module docstring)."""
    M = NX * norb
    kap = np.zeros((M, M), complex)
    for j in range(NX):
        for jp in range(NX):
            kap[j * norb:(j + 1) * norb, jp * norb:(jp + 1) * norb] = (
                X * np.exp(-1j * _K * (jp - j))[None, None]).mean(-1)
    kso = np.zeros((2 * M, 2 * M), complex)
    if eta == "singlet":
        kso[0::2, 1::2] = kap               # <c_{I up} c_{J dn}>
        kso[1::2, 0::2] = -kap.T            # <c_{J dn} c_{I up}>
        D = 0.5 * np.einsum("pqrs,sr->pq", A, kso)[0::2, 1::2]
    else:
        kso[0::2, 0::2] = kap               # <c_{I up} c_{J up}>, antisymmetric
        D = 0.5 * np.einsum("pqrs,sr->pq", A, kso)[0::2, 0::2]
    Dr = np.stack([D[0:norb, d * norb:(d + 1) * norb] for d in range(NX)], -1)
    return np.einsum("abd,kd->abk", Dr, np.exp(1j * np.outer(_K, np.arange(NX))))


def _static_x(g, phi):
    """``X_cd(k) = T sum_n G_ce(k, w_n) G_df(-k, -w_n) phi_ef(k)``."""
    gm = np.roll(g[:, :, ::-1, :, :, ::-1], 1, axis=2)[:, :, :, 0, 0]
    return np.einsum("cekn,dfkn,efk->cdk", g[:, :, :, 0, 0], gm, phi) / BETA


def _pair_field(rng, norb, eta):
    """A random gap of definite spin parity: ``phi_ab(k) = +-phi_ba(-k)``."""
    phi = rng.standard_normal((norb, norb, NX)) + 1j * rng.standard_normal((norb, norb, NX))
    sgn = 1.0 if eta == "singlet" else -1.0
    mk = (-np.arange(NX)) % NX
    return 0.5 * (phi + sgn * phi[:, :, mk].transpose(1, 0, 2))


def _declaration(t, R, a, b, v):
    """The Hermitian-closed declaration the file reader produces for the row
    ``(R, a, b)`` of type ``t`` at value ``v``."""
    return {t: {((R, 0, 0), (a, b)): v, ((-R, 0, 0), (b, a)): np.conj(v)}}


def _production_parts(norb, interactions, eta, g):
    """``(A, A_fl, A_pp)`` of the static bond kernel built by the production
    chain from ``interactions`` and the Green function ``g``."""
    shape = (NX, 1, 1)
    topo, view, split = sc._resolve_bond_topology(interactions, {"rvec": np.eye(3)}, norb, None)
    if norb == 1:
        kx = np.linspace(0, 2 * np.pi, NX, endpoint=False)
        carrier = sc._BondGreen(full_sc=g, deflated_kw=sc._green_sc_to_canonical(g),
                                green0_tail=None, coeff_tail_requested=0.0,
                                coeff_tail_applied=0.0, source="external_npz")
        (A, _), _, attr = sc._build_bond_operator(
            topo, split, carrier, norb, kx, np.zeros(1), np.zeros(1), BETA, eta,
            bond_max_shells=None, bond_memory_cap_gb=8.0, g2_tail=False)
        return A, attr["fl"], attr["pp"]
    # the two-orbital chain is not inversion symmetric, so the Hermitian-path
    # preconditions of the full operator do not apply; the same chain below
    # the operator assembly
    chi_bar = bc.bond_bubble(g, view, BETA).reshape(NX, view.n_channels * norb ** 2, -1)
    S, C = sc._build_bond_vertices(topo, split, norb, shape)
    dressed = ((ch, bc.dress_batch(chi_bar[None], W, ch, l0=0, nmat=0,
                                   spatial_shape=shape)[0][0])
               for ch, W in (("spin", S), ("charge", C)))
    A, A_fl, A_pp, _ = sc._bond_kernel_parts(S, C, dressed, g, view, norb, shape, BETA, eta,
                                             topo=topo)
    return A, A_fl, A_pp


#: (type, pairing channel) pairs whose first-order pairing the static bond
#: kernel reproduces exactly: all of them. The singlet of an off-site Exchange
#: or PairHop comes from ``eliashberg_bond.offsite_singlet_bare_vertex`` (it
#: was missing before #239), every other one from the channel-0 blocks.
EXACT = {(t, eta) for t in TYPES for eta in ("singlet", "triplet")}

#: (type, channel) pairs whose first-order reference vanishes identically:
#: a same-spin coupling has no (up, down) pairing, a spin-flip or pair
#: transfer no equal-spin one
ZERO_REFERENCE = {("Hund", "singlet"), ("Exchange", "triplet"), ("PairHop", "triplet")}


class _OracleCase(object):

    H = 1.0e-3      # central-difference step of the full-kernel coefficient

    def _check(self, norb, t, eta):
        g = _green_sc(norb)
        a, b = ROW[norb]
        fx = EDFixture(NX, norb, {}, [0.0] * norb, 1.0, 0.0)
        P2 = 2 * NX * norb
        rng = np.random.default_rng(7)
        for R in (+1, -1):
            # first-order coefficient of the reference: the vertex is linear
            A_ref = _antisymmetrized(_terms(t, fx, a, b, R, 1.0), P2)
            phi = _pair_field(rng, norb, eta)
            X = _static_x(g, phi)
            want = _reference_gap(norb, A_ref, X, eta)
            if (t, eta) in ZERO_REFERENCE:
                self.assertLess(np.abs(want).max(), 1e-14)
            else:
                # the check is not vacuous
                self.assertGreater(np.linalg.norm(want), 1e-2 * np.linalg.norm(X))
            scale = max(np.linalg.norm(want), np.linalg.norm(X))
            full = {}
            for J in (self.H, -self.H, self.H / 2, -self.H / 2):
                A, A_fl, A_pp = _production_parts(norb, _declaration(t, R, a, b, J), eta, g)
                pp = A_pp.matvec(phi.ravel()).reshape(norb, norb, NX)
                err = np.linalg.norm(pp / J - want) / scale
                self.assertLess(err, 1e-10, (norb, t, R, eta, J, "bare", err))
                full[J] = A.matvec(phi.ravel()).reshape(norb, norb, NX)
            # Richardson-extrapolated central difference: the fluctuation part
            # starts at O(J^2), even in J, so it cancels; the O(J^3) remainder
            # is removed to O(h^4) by the extrapolation
            d1 = (full[self.H] - full[-self.H]) / (2 * self.H)
            d2 = (full[self.H / 2] - full[-self.H / 2]) / self.H
            coeff = (4.0 * d2 - d1) / 3.0
            err = np.linalg.norm(coeff - want) / scale
            self.assertLess(err, 1e-8, (norb, t, R, eta, "full", err))

    def _check_all(self, t, eta):
        for norb in (1, 2):
            self._check(norb, t, eta)


class TestFirstOrderPairingOracle(_OracleCase, unittest.TestCase):

    def test_exact_types_and_channels(self):
        for t, eta in sorted(EXACT):
            for norb in (1, 2):
                with self.subTest(t=t, eta=eta, norb=norb):
                    self._check(norb, t, eta)

    def test_offsite_type_moves_the_single_orbital_result(self):
        """At ``norb = 1`` every newly bond-resolved type changes the static
        kernel at the coupling used in the issue's comparison (0.3 on top of
        a nearest-neighbour V); a type that were silently dropped would leave
        it unchanged."""
        g = _green_sc(1)
        rng = np.random.default_rng(3)
        base = _declaration("CoulombInter", 1, 0, 0, 0.5)
        base["CoulombIntra"] = {((0, 0, 0), (0, 0)): 1.0}
        for eta in ("singlet", "triplet"):
            phi = _pair_field(rng, 1, eta).ravel()
            ref = _production_parts(1, base, eta, g)[0].matvec(phi)
            for t in TYPES[1:]:
                inter = dict(base, **_declaration(t, 1, 0, 0, 0.3))
                got = _production_parts(1, inter, eta, g)[0].matvec(phi)
                self.assertGreater(np.linalg.norm(got - ref), 1e-3 * np.linalg.norm(ref),
                                   (t, eta))


class TestLambdaPpDerivative(unittest.TestCase):
    """The derivative at ``J = 0`` of the bare-part attribution
    ``lambda_pp`` (``bond_channels.attribute_lambda``) of a single-band model
    with ``U`` and a nearest-neighbour ``V``, along the nearest-neighbour
    form factor of the channel (``cos k`` singlet, ``sin k`` triplet, and the
    on-site constant for the singlet of a PairHop, which moves on-site pairs; the
    leading eigenvector of this short chain's triplet kernel is a higher
    harmonic with almost no nearest-neighbour weight), equals the Rayleigh
    quotient of the reference first-order vertex on that vector -- for every
    type and channel the kernel reproduces exactly (``EXACT``). The derivative is a
    Richardson-extrapolated central difference over ``J = +-0.01, +-0.02``
    through the whole production chain."""

    def test_derivative_matches_the_reference(self):
        from hwave.solver.bond_channels import attribute_lambda, pair_weight
        g = _green_sc(1)
        W = pair_weight(g, BETA)
        fx = EDFixture(NX, 1, {}, [0.0], 1.0, 0.0)
        base = _declaration("CoulombInter", 1, 0, 0, 0.5)
        base["CoulombIntra"] = {((0, 0, 0), (0, 0)): 2.0}
        n = NX
        for eta in ("singlet", "triplet"):
            for t in TYPES:
                if (t, eta) not in EXACT:
                    continue
                # an off-site PairHop moves ON-SITE pairs, so its singlet
                # acts on the on-site (constant) form factor
                v0 = (np.ones(NX) if (t, eta) == ("PairHop", "singlet") else
                      np.cos(_K) if eta == "singlet" else np.sin(_K)).astype(complex)
                A_ref = _antisymmetrized(_terms(t, fx, 0, 0, 1, 1.0), 2 * NX)

                def ref_op(v, A_ref=A_ref, eta=eta):
                    phi = sc._project_gap_parity(np.asarray(v).reshape(1, 1, NX, 1, 1),
                                                 eta)[:, :, :, 0, 0]
                    return _reference_gap(1, A_ref, _static_x(g, phi), eta).ravel()

                from scipy.sparse.linalg import LinearOperator
                O = LinearOperator((n, n), matvec=ref_op, dtype=complex)
                r1 = attribute_lambda(v0, W, O, O, op_full=O)["lambda"]
                lam = {}
                for J in (0.01, -0.01, 0.02, -0.02):
                    inter = dict(base)
                    inter.update(_declaration(t, 1, 0, 0, J)
                                 if t != "CoulombInter" else
                                 {"CoulombInter": {k: v + J for k, v in
                                                   base["CoulombInter"].items()}})
                    _, A_fl, A_pp = _production_parts(1, inter, eta, g)
                    lam[J] = attribute_lambda(v0, W, A_pp, A_fl)["lambda_pp"]
                d1 = (lam[0.02] - lam[-0.02]) / 0.04
                d2 = (lam[0.01] - lam[-0.01]) / 0.02
                deriv = (4.0 * d2 - d1) / 3.0
                with self.subTest(t=t, eta=eta):
                    self.assertAlmostEqual(deriv, r1, delta=1e-8 * max(1.0, abs(r1)))
                    if (t, eta) not in ZERO_REFERENCE:
                        self.assertGreater(abs(r1), 1e-3)


if __name__ == "__main__":
    unittest.main()
