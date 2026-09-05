"""The bond-resolved longitudinal S/C vertex builder and the per-channel
dressing helper (GitHub issue #181, Tier 3 Phase A).

``build_sc_bond_channel`` places the Tier-1 locality-split pair-space
matrix verbatim as the channel-0 block and adds, for every off-site
reversal orbit, the bond-DIAGONAL Fock entries ``w_t * Re v`` at
``(m, a, b)`` and ``(reverse[m], b, a)``, ``w_t`` the type's adjudicated
``cross`` coefficient from ``vertex_table`` (S or C column). Gate G1:
for CoulombInter it must equal the Eliashberg module's
``bare_bond_vertices`` exactly.
"""
import unittest

import numpy as np

from hwave.solver import bond_channels as bc
from hwave.solver.vertex_table import sc_coefficients


def _topo(decls, norb, types=bc._LONGITUDINAL_ACTIVE_TYPES):
    return bc.resolve_bond_topology(decls, np.eye(3), norb, active_types=types)


def _two_sided(t, R, a, b, v):
    return {t: {((R, 0, 0), (a, b)): v, ((-R, 0, 0), (b, a)): np.conj(v)}}


def _sc0(decl, norb, L, types=("CoulombInter", "Hund", "Ising")):
    """(S0, C0) at (L, nd, nd) from the locality split of ``decl`` on an
    L-site ring, via the production owner of the channel-0 block."""
    from hwave.solver.offsite import split_locality, sc_matrices_from_split

    class H:
        param_ham = decl
        param_ham_orig = None

        def _reshape_interaction(self, tbl, so):
            return tbl

    class Lat:
        subshape = (1, 1, 1)

    S0, C0 = sc_matrices_from_split(split_locality(H(), Lat()), types, norb, L, 1, 1)
    nd = norb * norb
    return S0.reshape(L, nd, nd), C0.reshape(L, nd, nd)


class TestBuildScBondChannel(unittest.TestCase):

    def test_channel_zero_verbatim_and_bond_diagonal_fock(self):
        nvol, norb = 5, 2
        nd = norb * norb
        rng = np.random.default_rng(0)
        S0 = rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd))
        topo = _topo(_two_sided("Hund", 1, 0, 1, 0.3), norb)
        S = bc.build_sc_bond_channel(topo, S0, "S")
        C = bc.build_sc_bond_channel(topo, S0, "C")
        B = topo.delta_r.shape[0]
        self.assertEqual(S.shape, (nvol, B * nd, B * nd))
        self.assertEqual(S.dtype, np.complex128)
        self.assertTrue(np.array_equal(S[:, :nd, :nd], S0))
        self.assertTrue(np.array_equal(C[:, :nd, :nd], S0))
        s_h, c_h = sc_coefficients("Hund", "cross")     # (-1, +1)
        m_plus = [tuple(r) for r in topo.delta_r].index((1, 0, 0))
        m_minus = int(topo.reverse[m_plus])
        I1 = m_plus * nd + 0 * norb + 1                  # (m, a=0, b=1)
        I2 = m_minus * nd + 1 * norb + 0                 # (reverse m, b=1, a=0)
        for I in (I1, I2):
            np.testing.assert_allclose(S[:, I, I], s_h * 0.3)
            np.testing.assert_allclose(C[:, I, I], c_h * 0.3)
        mask = np.ones((B * nd, B * nd), bool)
        mask[:nd, :nd] = False
        mask[I1, I1] = mask[I2, I2] = False
        self.assertEqual(np.abs(S[:, mask]).max(), 0.0)
        self.assertEqual(np.abs(C[:, mask]).max(), 0.0)

    def test_same_orbital_bond_writes_two_slots(self):
        topo = _topo(_two_sided("Ising", 1, 0, 0, 0.5), 1)
        S = bc.build_sc_bond_channel(topo, np.zeros((3, 1, 1), complex), "S")
        C = bc.build_sc_bond_channel(topo, np.zeros((3, 1, 1), complex), "C")
        s_i, c_i = sc_coefficients("Ising", "cross")       # (+1, -1)
        np.testing.assert_allclose(np.diagonal(S, axis1=1, axis2=2)[0],
                                   [0.0, s_i * 0.5, s_i * 0.5])
        np.testing.assert_allclose(np.diagonal(C, axis1=1, axis2=2)[0],
                                   [0.0, c_i * 0.5, c_i * 0.5])

    def test_types_selects_and_orders_nothing_else(self):
        decl = _two_sided("Hund", 1, 0, 0, 0.3)
        decl.update(_two_sided("Ising", 1, 0, 0, 0.2))
        topo = _topo(decl, 1)
        W0 = np.zeros((3, 1, 1), complex)
        both = bc.build_sc_bond_channel(topo, W0, "S")
        hund = bc.build_sc_bond_channel(topo, W0, "S", types=("Hund",))
        ising = bc.build_sc_bond_channel(topo, W0, "S", types=("Ising",))
        np.testing.assert_allclose(both, hund + ising)
        np.testing.assert_allclose(np.diagonal(hund, axis1=1, axis2=2)[0],
                                   [0.0, -0.3, -0.3])
        np.testing.assert_allclose(np.diagonal(ising, axis1=1, axis2=2)[0],
                                   [0.0, 0.2, 0.2])
        empty = bc.build_sc_bond_channel(topo, W0, "S", types=())
        self.assertEqual(np.abs(empty).max(), 0.0)         # channel 0 only

    def test_matches_bare_bond_vertices_for_coulombinter(self):
        """G1: CoulombInter through the new builder equals the Eliashberg
        objects VERBATIM on a norb=1 two-shell topology, on the case-M
        inter-orbital bond (norb=2) and on a mixed same/inter-orbital
        declaration. Both channel-0 builders (``_build_bond_m0_blocks``
        and the Tier-1 locality split) and both bubbles
        (``bond_bubble``/``bubble.bond_bubble_static``, whose channel-0
        block is the solver's static ``chi0q`` slice bit-for-bit) share
        one pair frame."""
        from hwave.sc import _build_bond_m0_blocks
        L = 5
        for norb, entries in ((1, [(1, 0, 0, 0.7), (2, 0, 0, 0.2)]),
                              (2, [(1, 0, 1, 0.4)]),
                              (2, [(1, 0, 0, 0.3), (1, 0, 1, 0.4), (1, 1, 1, 0.1)])):
            with self.subTest(norb=norb, entries=entries):
                decl = {}
                for (R, a, b, v) in entries:
                    decl[((R, 0, 0), (a, b))] = v
                    decl[((-R, 0, 0), (b, a))] = np.conj(v)
                bond_set = bc.resolve_interactions(decl, np.eye(3), norb=norb)
                kx = 2.0 * np.pi * np.arange(L) / L
                S0q, C0q = _build_bond_m0_blocks(
                    bond_set, {}, {}, norb, kx, np.array([0.0]), np.array([0.0]))
                S_ref, C_ref, _, _ = bc.bare_bond_vertices(bond_set, S0q, C0q, norb)
                S_ref = S_ref[:, 0, 0]
                C_ref = C_ref[:, 0, 0]
                S0, C0 = _sc0({"CoulombInter": decl}, norb, L)
                topo = _topo({"CoulombInter": decl}, norb)
                self.assertEqual([tuple(r) for r in topo.delta_r],
                                 [tuple(r) for r in bond_set.delta_r])
                S_new, C_new = bc.W_sc_bond(topo, S0, C0)
                np.testing.assert_allclose(S_new, S_ref, rtol=0, atol=1e-13)
                np.testing.assert_allclose(C_new, C_ref, rtol=0, atol=1e-13)

    def test_hermitian_per_q_and_transpose_under_q_reversal(self):
        """Structure on a 5-point ring with an a != b record: W(q) is
        Hermitian at every q, and W(-q) == W(q)^T (real coefficients:
        the density slot (aa,bb) carries V_ab(q) = V_ab(-q)^*, the bond-
        diagonal blocks are real and q-independent)."""
        L = 5
        decl = _two_sided("Hund", 1, 0, 1, 0.3)
        decl.update(_two_sided("CoulombInter", 1, 0, 0, 0.2))
        S0, C0 = _sc0(decl, 2, L)
        S, C = bc.W_sc_bond(_topo(decl, 2), S0, C0)
        for q in range(L):
            qm = (-q) % L
            np.testing.assert_allclose(S[q].conj().T, S[q], atol=1e-13)
            np.testing.assert_allclose(C[q].conj().T, C[q], atol=1e-13)
            np.testing.assert_allclose(S[q].T, S[qm], atol=1e-13)
            np.testing.assert_allclose(C[q].T, C[qm], atol=1e-13)
        self.assertGreater(np.abs(C[1] - C[1].T).max(), 0.1)   # not trivially symmetric

    def test_validation(self):
        z = np.zeros((3, 1, 1), complex)
        topo_c = _topo(_two_sided("CoulombInter", 1, 0, 0, 0.1 + 0.2j), 1)
        with self.assertRaises(ValueError) as cm:                # complex coefficient
            bc.build_sc_bond_channel(topo_c, z, "S")
        self.assertIn("complex", str(cm.exception))
        bc.build_sc_bond_channel(topo_c, z, "S", imag_tol=1.0)   # tolerance honoured
        topo = _topo(_two_sided("CoulombInter", 1, 0, 0, 0.1), 1)
        with self.assertRaises(ValueError):
            bc.build_sc_bond_channel(topo, z, "X")
        with self.assertRaises(ValueError):                      # wrong W0 shape
            bc.build_sc_bond_channel(topo, np.zeros((3, 2, 2), complex), "S")
        with self.assertRaises(ValueError):                      # non-finite W0
            bc.build_sc_bond_channel(topo, np.full((3, 1, 1), np.nan), "S")
        with self.assertRaises(ValueError):                      # unknown type
            bc.build_sc_bond_channel(topo, z, "S", types=("Exchange",))
        with self.assertRaises(ValueError):                      # duplicate type
            bc.build_sc_bond_channel(topo, z, "S", types=("CoulombInter", "CoulombInter"))
        with self.assertRaises(ValueError):                      # bad imag_tol
            bc.build_sc_bond_channel(topo, z, "S", imag_tol=-1.0)
        with self.assertRaises(ValueError):                      # mutated topology
            bad = bc.BondTopology(topo.delta_r, topo.reverse, topo.coeffs)
            bad.coeffs["CoulombInter"] = np.zeros((7, 1, 1), complex)
            bc.build_sc_bond_channel(bad, z, "S")


class TestDressChannel(unittest.TestCase):

    def _random_problem(self, nvol=6, ND=3, seed=1):
        rng = np.random.default_rng(seed)
        chi_bar = 0.1 * (rng.normal(size=(nvol, ND, ND))
                         + 1j * rng.normal(size=(nvol, ND, ND)))
        W = 0.5 * rng.normal(size=(nvol, ND, ND))
        return chi_bar, W

    def test_matches_direct_formula_and_reports_cond(self):
        chi_bar, W = self._random_problem()
        chi, cond = bc.dress_channel(chi_bar, W, "spin", spatial_shape=(3, 2, 1))
        ND = W.shape[1]
        mat = np.eye(ND)[None] - chi_bar @ W
        np.testing.assert_array_equal(chi, np.linalg.solve(mat, chi_bar))
        sv = np.linalg.svd(mat, compute_uv=False)
        expected = float(np.min(np.minimum(
            sv[:, -1] / sv[:, 0], sv[:, -1] / np.maximum(1.0, sv[:, 0]))))
        self.assertIsInstance(cond, float)
        self.assertAlmostEqual(cond, expected, places=14)
        chi_c, cond_c = bc.dress_channel(chi_bar, W, "charge", spatial_shape=(3, 2, 1))
        np.testing.assert_array_equal(
            chi_c, np.linalg.solve(np.eye(ND)[None] + chi_bar @ W, chi_bar))
        self.assertGreater(cond_c, 0.0)
        # cond_tol=None: no guard, no score
        _, none = bc.dress_channel(chi_bar, W, "spin", spatial_shape=(3, 2, 1), cond_tol=None)
        self.assertIsNone(none)

    def test_dress_bond_is_two_dress_channel_calls(self):
        chi_bar, W = self._random_problem(nvol=4, ND=2)
        S = W
        C = -0.3 * W
        chi_s, chi_c = bc.dress_bond(chi_bar.reshape(2, 2, 1, 2, 2),
                                     S.reshape(2, 2, 1, 2, 2),
                                     C.reshape(2, 2, 1, 2, 2))
        cs, _ = bc.dress_channel(chi_bar, S, "spin", spatial_shape=(2, 2, 1))
        cc, _ = bc.dress_channel(chi_bar, C, "charge", spatial_shape=(2, 2, 1))
        np.testing.assert_array_equal(chi_s.reshape(4, 2, 2), cs)
        np.testing.assert_array_equal(chi_c.reshape(4, 2, 2), cc)
        # and both equal the former inline formula bit-for-bit
        I = np.eye(2)[None]
        np.testing.assert_array_equal(chi_s.reshape(4, 2, 2),
                                      np.linalg.solve(I - chi_bar @ S, chi_bar))
        np.testing.assert_array_equal(chi_c.reshape(4, 2, 2),
                                      np.linalg.solve(I + chi_bar @ C, chi_bar))

    def test_refusal_names_channel_and_q(self):
        chi_bar = np.zeros((2, 1, 1), complex)
        chi_bar[1, 0, 0] = 1.0
        W = np.ones((2, 1, 1))                        # q=1: 1 - 1 = 0 -> singular
        with self.assertRaises(ValueError) as cm:
            bc.dress_channel(chi_bar, W, "spin", spatial_shape=(2, 1, 1))
        self.assertIn("spin", str(cm.exception))
        self.assertIn("(1, 0, 0)", str(cm.exception))
        # the charge channel of the same data is regular
        bc.dress_channel(chi_bar, W, "charge", spatial_shape=(2, 1, 1))

    def test_validation(self):
        chi_bar, W = self._random_problem(nvol=4, ND=2)
        with self.assertRaises(ValueError):
            bc.dress_channel(chi_bar, W, "up", spatial_shape=(2, 2, 1))
        with self.assertRaises(ValueError):           # nvol != prod(spatial_shape)
            bc.dress_channel(chi_bar, W, "spin", spatial_shape=(3, 1, 1))
        with self.assertRaises(ValueError):           # shape mismatch
            bc.dress_channel(chi_bar, W[:, :1, :1], "spin", spatial_shape=(2, 2, 1))
        with self.assertRaises(ValueError):           # not square
            bc.dress_channel(chi_bar[:, :, :1], W[:, :, :1], "spin", spatial_shape=(2, 2, 1))


if __name__ == "__main__":
    unittest.main()
