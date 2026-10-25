"""The bare (first-order) TRIPLET pairing vertex of an on-site inter-orbital
interaction (issue #262).

An inter-orbital density-type term (CoulombInter ``U'``, Hund, Ising at
``R = 0``, ``a != b``) sits in the spin/charge matrices twice: its Hartree
content in the density slot ``(aa, bb)`` and its Fock content in the cross slot
``(ab, ab)``. In the pairing channel the cross slot is the exchange image of
the density slot, so ``1/2 (C - S)`` over both counted the term twice in the
triplet channel. The checks below compare the first-order kernel of every
route that builds the term -- the chi0q route (``_compute_vertices_general``
at a zero bubble), the FLEX route (``_compute_vertices_flex`` at zero
susceptibilities, both conventions) and the bond-resolved static kernel --
with the independent spin-orbital construction of
``tests/test_sc_bond_pairing_oracle.py``, on its inter-orbital,
mirror-asymmetric two-orbital chain, for every type and both channels, and
for off-site rows on the non-bond routes, which collapse them onto the same
slots.
"""
import unittest

import numpy as np

import hwave.sc as sc
import tests.test_sc_bond_pairing_oracle as O

NORB = 2
ONSITE_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop")
#: the reduced (Kuroki) builder refuses Exchange and PairHop
KUROKI_TYPES = ("CoulombInter", "Hund", "Ising")


def _onsite_declaration(t, v):
    return {t: {((0, 0, 0), (0, 1)): v, ((0, 0, 0), (1, 0)): v}}


def _reference(t, eta, g, rng):
    """``(phi, Delta_ref)``: a random gap of the channel's parity and the
    first-order reference gap of the on-site row ``(0, 0, 1)`` of ``t`` at unit
    coupling."""
    fx = O.EDFixture(O.NX, NORB, {}, [0.0] * NORB, 1.0, 0.0)
    A = O._antisymmetrized(O._terms(t, fx, 0, 1, 0, 1.0), 2 * O.NX * NORB)
    phi = O._pair_field(rng, NORB, eta)
    return phi, O._reference_gap(NORB, A, O._static_x(g, phi), eta)


class TestTripletBareVertex(unittest.TestCase):

    def setUp(self):
        self.g = O._green_sc(NORB)
        self.kx = np.linspace(0, 2 * np.pi, O.NX, endpoint=False)
        self.G2 = sc._calc_g2(self.g, O.BETA, tail=False)

    def _assert_matches(self, got, want, label):
        scale = max(np.linalg.norm(want), 1e-2)
        self.assertLess(np.linalg.norm(got - want) / scale, 1e-10, label)

    def _nonbond(self, vertex):
        A, _ = sc._make_kernel_operator(vertex, self.G2, NORB, O.NX, 1, 1)
        return A

    def test_flex_route_both_conventions(self):
        rng = np.random.default_rng(5)
        z = np.zeros((O.NX, 1, 1, NORB ** 2, NORB ** 2), complex)
        for conv, types in (("myo", ONSITE_TYPES), ("kuroki", KUROKI_TYPES)):
            for t in types:
                inter_k = sc._build_interaction_k(self.kx, np.zeros(1), np.zeros(1),
                                                  _onsite_declaration(t, 1.0), NORB)
                for eta in ("singlet", "triplet"):
                    with self.subTest(conv=conv, t=t, eta=eta):
                        phi, want = _reference(t, eta, self.g, rng)
                        V = sc._compute_vertices_flex(z, z, inter_k, NORB, O.NX, 1, 1,
                                                      pairing_type=eta, convention=conv)
                        got = self._nonbond(V).matvec(phi.ravel()).reshape(want.shape)
                        self._assert_matches(got, want, (conv, t, eta))

    def test_chi0q_route(self):
        rng = np.random.default_rng(6)
        nmat = 2
        chi0q = np.zeros((NORB,) * 4 + (O.NX, 1, 1, nmat), complex)
        for t in ONSITE_TYPES:
            inter_k = sc._build_interaction_k(self.kx, np.zeros(1), np.zeros(1),
                                              _onsite_declaration(t, 1.0), NORB)
            for eta in ("singlet", "triplet"):
                with self.subTest(t=t, eta=eta):
                    phi, want = _reference(t, eta, self.g, rng)
                    V = sc._compute_vertices_general(chi0q, inter_k, NORB, O.NX, 1, 1, nmat,
                                                     pairing_type=eta)
                    got = self._nonbond(V).matvec(phi.ravel()).reshape(want.shape)
                    self._assert_matches(got, want, (t, eta))

    def test_bond_static_kernel(self):
        """The bare part of the bond-resolved static kernel (its channel-0
        block goes through the same term)."""
        rng = np.random.default_rng(7)
        for t in ONSITE_TYPES:
            for eta in ("singlet", "triplet"):
                with self.subTest(t=t, eta=eta):
                    phi, want = _reference(t, eta, self.g, rng)
                    J = 1.0e-3
                    decl = _onsite_declaration(t, J)
                    # the bond path needs a declared off-site row; a zero one
                    decl.setdefault("CoulombInter", {}).update(
                        {((1, 0, 0), (0, 1)): 0.0, ((-1, 0, 0), (1, 0)): 0.0})
                    _, _, A_pp = O._production_parts(NORB, decl, eta, self.g)
                    got = A_pp.matvec(phi.ravel()).reshape(want.shape) / J
                    self._assert_matches(got, want, (t, eta))

    def test_offsite_rows_on_the_nonbond_routes(self):
        """Built from a plain table (no locality split), the vertex collapses
        an off-site row onto the q-only vertex, its Fock content into the
        cross slots -- what the non-bond path did before issue #239: the
        triplet bare term of an off-site density-type row is exact once those
        slots are left out (it was off by up to 3x). The singlet of an
        off-site Exchange or PairHop is not representable there (the non-bond
        entries refuse it since #239) and is not checked."""
        rng = np.random.default_rng(9)
        z = np.zeros((O.NX, 1, 1, NORB ** 2, NORB ** 2), complex)
        chi0q = np.zeros((NORB,) * 4 + (O.NX, 1, 1, 2), complex)
        fx = O.EDFixture(O.NX, NORB, {}, [0.0] * NORB, 1.0, 0.0)
        cases = [(t, eta) for t in ONSITE_TYPES for eta in ("singlet", "triplet")
                 if not (t in ("Exchange", "PairHop") and eta == "singlet")]
        for R in (1, -1):
            for t, eta in cases:
                A = O._antisymmetrized(O._terms(t, fx, 0, 1, R, 1.0), 2 * O.NX * NORB)
                phi = O._pair_field(rng, NORB, eta)
                want = O._reference_gap(NORB, A, O._static_x(self.g, phi), eta)
                inter_k = sc._build_interaction_k(self.kx, np.zeros(1), np.zeros(1),
                                                  O._declaration(t, R, 0, 1, 1.0), NORB)
                for route in ("flex", "general"):
                    with self.subTest(R=R, t=t, eta=eta, route=route):
                        V = (sc._compute_vertices_flex(z, z, inter_k, NORB, O.NX, 1, 1,
                                                       pairing_type=eta, convention="myo")
                             if route == "flex" else
                             sc._compute_vertices_general(chi0q, inter_k, NORB, O.NX, 1, 1, 2,
                                                          pairing_type=eta))
                        got = self._nonbond(V).matvec(phi.ravel()).reshape(want.shape)
                        self._assert_matches(got, want, (R, t, eta, route))

    def test_the_triplet_reference_is_not_vacuous(self):
        """The density-type triplet references are nonzero, so the triplet
        checks above would catch the doubled term."""
        rng = np.random.default_rng(8)
        for t in KUROKI_TYPES:
            phi, want = _reference(t, "triplet", self.g, rng)
            self.assertGreater(np.linalg.norm(want), 1e-2 * np.linalg.norm(phi), t)


if __name__ == "__main__":
    unittest.main()
