"""Spin-orbital Eliashberg (hwave_sc) pairing vertex tests (#83).

The vertex under test is :func:`hwave.solver.so_pairing.so_pairing_vertex`:
from the antisymmetrized spinful vertex ``Gamma`` (``rpa.spinful_bubble_vertex``,
the vertex the spin-orbital RPA and FLEX resum) and the static spinful
susceptibility ``chi = [1 + chi0 Gamma]^-1 chi0`` it builds the pairing vertex
``V`` of the spin-orbital gap equation in the layout ``sc._make_kernel_operator``
consumes (``Delta_ad = -sum_bc V_abcd F_bc``).

* ``TestFirstOrderOracle``: the bare part against an independent first-order
  construction in the spin-orbital basis of a finite chain with complex
  spin-flip hopping (the quartic operator terms of every on-site type,
  antisymmetrized and contracted with a general antisymmetric pair
  amplitude over all spin blocks).
* ``TestSU2Limit``: with spin-independent hopping and a spin-rotation
  invariant interaction, the kernel acting on ``phi (x) i sigma_y`` (and on the
  three triplet components ``phi (x) sigma_j i sigma_y``) equals the mapped
  action of the existing non-spin-orbital singlet (triplet) kernel, bare and
  fluctuation parts together. Hund and Exchange alone conserve S_z but are not
  rotation invariant; for them only the equal-spin triplet (the non-SO
  triplet's own definition) is an exact limit. Their singlet and ``d_z``
  (opposite-spin) components are not compared: the non-spin-orbital singlet
  and triplet formulas assume spin-rotation invariance and are only
  approximate for these terms (issue #88), while the spin-orbital vertex
  treats them exactly. For Hund alone, which acts between equal spins only,
  the spin-orbital kernel returns exactly zero opposite-spin pairing.
"""
import shutil
import tempfile
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import bubble, rpa
import tests.test_flex_spin_orbital as fso
import tests.test_sc_bond_pairing_oracle as O
from tests.test_bond_transverse_ed import _terms_pairlift_offsite
from tests.so_bond_fixtures import onsite_decl as _onsite_decl, rpa_spinful_gamma

NX, BETA, NMAT, MU = O.NX, O.BETA, O.NMAT, O.MU
NP = 2            # physical orbitals
ND = 2 * NP       # spin-orbitals, spin-block order s * NP + a
_K = O._K
_WN = O._WN

#: Pauli matrices and the spin structure of the singlet / triplet gaps
_IY = np.array([[0, 1], [-1, 0]], complex)
_SX = np.array([[0, 1], [1, 0]], complex)
_SY = np.array([[0, -1j], [1j, 0]], complex)
_SZ = np.diag([1.0, -1.0]).astype(complex)
_TRIPLET = {"z": _SZ @ _IY, "x": _SX @ _IY, "y": _SY @ _IY}


def _gamma(decl):
    """``Gamma (NX, ND^2, ND^2)`` of the spin-orbital RPA built from the
    declarations ``{type: {((R, 0, 0), (a, b)): v}}`` (physical orbital
    indices, as the interaction files are read in spin-orbital mode)."""
    return rpa_spinful_gamma(decl, NP, (NX, 1, 1), beta=BETA, value_format="%.12f")


def _onsite_terms(t, fx, v):
    """Quartic terms ``v c+_p c_q c+_r c_s`` of :func:`_onsite_decl` on the
    ``NX``-site chain (independent of every H-wave vertex builder)."""
    if t == "CoulombIntra":
        return [(fx.mode(j, a, 0), fx.mode(j, a, 0), fx.mode(j, a, 1), fx.mode(j, a, 1),
                 v * (1.0 if a == 0 else 0.8))
                for j in range(NX) for a in range(NP)]
    if t == "PairLift":
        return _terms_pairlift_offsite(fx, 0, 1, 0, v)
    return O._terms(t, fx, 0, 1, 0, v)


def _green(h):
    """``G(k, i w_n)`` in the sc.py layout ``(n, n, NX, 1, 1, NMAT)``."""
    n = h.shape[1]
    g = np.empty((n, n, NX, 1, 1, NMAT), complex)
    for i, w in enumerate(_WN):
        g[:, :, :, 0, 0, i] = np.linalg.inv(
            (1j * w + MU) * np.eye(n)[None] - h).transpose(1, 2, 0)
    return g


def _h_orbital():
    """The mirror- and transpose-asymmetric two-orbital chain of the oracle,
    ``(NX, NP, NP)``."""
    h = np.zeros((NX, NP, NP), complex)
    for (r, a, b), t in O.TRANSFER[2].items():
        h[:, a, b] += t * np.exp(1j * _K * r)
    return h


def _h_spin_orbital(soc):
    """The same chain in the spin-orbital basis (spin-block order) plus a
    complex spin-flip (Rashba-like) hopping of strength ``soc``."""
    h = np.zeros((NX, ND, ND), complex)
    h0 = _h_orbital()
    h[:, :NP, :NP] = h0
    h[:, NP:, NP:] = h0
    rng = np.random.default_rng(11)
    for (r, a, b) in ((1, 0, 1), (1, 0, 0), (0, 0, 1), (1, 1, 1)):
        lam = soc * (rng.standard_normal() + 1j * rng.standard_normal())
        hop = np.zeros((ND, ND), complex)
        hop[a, NP + b] = lam                      # c+_{j, a, up} c_{j+r, b, dn}
        hk = hop[None] * np.exp(1j * _K * r)[:, None, None]
        h += hk + hk.conj().transpose(0, 2, 1)
    return h


def _mode_of_spin_block():
    """EDFixture's in-site mode index ``2 a + s`` of spin-block index
    ``s NP + a``."""
    return np.array([2 * (A % NP) + A // NP for A in range(ND)])


def _reference_gap_so(A, X):
    """``Delta_AB(k)`` (spin-block) of the first-order gap equation
    ``D_pq = 1/2 sum_rs A_pq,rs kappa_sr`` for the general (all spin blocks)
    antisymmetric spin-orbital pair amplitude ``X_AB(k)``: the
    construction of ``O._reference_gap`` without the restriction to the
    ``(up, down)`` or ``(up, up)`` block."""
    perm = _mode_of_spin_block()
    M = NX * ND
    kap = np.zeros((M, M), complex)
    for j in range(NX):
        for jp in range(NX):
            blk = (X * np.exp(-1j * _K * (jp - j))[None, None]).mean(-1)
            kap[np.ix_(j * ND + perm, jp * ND + perm)] = blk
    D = 0.5 * np.einsum("pqrs,sr->pq", A, kap)
    Dr = np.stack([D[np.ix_(perm, d * ND + perm)] for d in range(NX)], -1)
    return np.einsum("abd,kd->abk", Dr, np.exp(1j * np.outer(_K, np.arange(NX))))


def _antisymmetric_gap(rng):
    """A random spin-orbital gap with ``Delta_AB(k) = -Delta_BA(-k)``."""
    phi = (rng.standard_normal((ND, ND, NX, 1, 1))
           + 1j * rng.standard_normal((ND, ND, NX, 1, 1)))
    return sc._project_gap_parity(phi, "triplet")[:, :, :, 0, 0]


def _static_bubble(g):
    """Static (zero bosonic frequency) general bubble ``(NX, n, n, n, n)`` in
    the ``[a, c, b, d]`` layout, from the sc-layout Green function."""
    gk = g[:, :, :, 0, 0, :].transpose(3, 2, 0, 1)[None]
    chi0 = bubble.dense_bubble(gk, None, BETA, spatial_shape=(NX, 1, 1), scheme="general")
    return chi0[0, NMAT // 2]


class TestFirstOrderOracle(unittest.TestCase):
    """The bare part, with spin-orbit coupling, against the independent
    first-order spin-orbital construction, for every on-site type."""

    TYPES = (("CoulombIntra", 1.0), ("CoulombInter", 1.0), ("Hund", 1.0),
             ("Ising", 1.0), ("Exchange", 1.0), ("PairHop", 0.7 + 0.4j),
             ("PairLift", 1.0))

    def test_bare_part_matches_spin_orbital_reference(self):
        from hwave.solver.so_pairing import so_pairing_vertex, so_pairing_vertex_parts
        g = _green(_h_spin_orbital(0.3))
        G2 = sc._calc_g2(g, BETA, tail=False)
        fx = O.EDFixture(NX, NP, {}, [0.0] * NP, 1.0, 0.0)
        rng = np.random.default_rng(1)
        zero = np.zeros((NX, ND * ND, ND * ND), complex)
        # the Green function has spin-flip content: the spin-orbit coupling
        # is load bearing
        self.assertGreater(np.linalg.norm(g[:NP, NP:]), 1e-2 * np.linalg.norm(g))
        for t, v in self.TYPES:
            with self.subTest(t=t):
                A = O._antisymmetrized(_onsite_terms(t, fx, v), 2 * NX * NP)
                phi = _antisymmetric_gap(rng)
                X = O._static_x(g, phi)
                want = _reference_gap_so(A, X)
                self.assertGreater(np.linalg.norm(want), 1e-2 * np.linalg.norm(X))
                gam = _gamma(_onsite_decl(t, v))
                V_bare, V_fl = so_pairing_vertex_parts(gam, zero, ND, (NX, 1, 1))
                self.assertEqual(V_bare.shape, (ND,) * 4 + (NX, 1, 1))
                self.assertEqual(np.abs(V_fl).max(), 0.0)
                for V in (V_bare, so_pairing_vertex(gam, zero, ND, (NX, 1, 1))):
                    K, _ = sc._make_kernel_operator(V, G2, ND, NX, 1, 1)
                    got = K.matvec(phi.ravel()).reshape(want.shape)
                    err = np.linalg.norm(got - want) / np.linalg.norm(want)
                    self.assertLess(err, 1e-10, (t, err))


class TestSU2Limit(unittest.TestCase):
    """Spin-independent hopping: the spin-orbital kernel reproduces the
    non-spin-orbital singlet and triplet kernels (``sc._compute_vertices_general``
    with the same chi0), bare plus fluctuation part."""

    @classmethod
    def setUpClass(cls):
        h = _h_orbital()
        cls.g = _green(h)
        hso = np.zeros((NX, ND, ND), complex)
        hso[:, :NP, :NP] = h
        hso[:, NP:, NP:] = h
        cls.gso = _green(hso)
        cls.G2 = sc._calc_g2(cls.g, BETA, tail=False)
        cls.G2so = sc._calc_g2(cls.gso, BETA, tail=False)
        chi0 = _static_bubble(cls.g)
        cls.chi0q = chi0.transpose(1, 2, 3, 4, 0).reshape(NP, NP, NP, NP, NX, 1, 1, 1)
        cls.chi0so = _static_bubble(cls.gso).reshape(NX, ND * ND, ND * ND)

    def _check(self, decl, components):
        """``components``: the (channel, spin matrix) pairs to compare."""
        from hwave.solver.so_pairing import so_pairing_vertex
        kx = np.linspace(0, 2 * np.pi, NX, endpoint=False)
        inter_k = sc._build_interaction_k(kx, np.zeros(1), np.zeros(1), decl, NP)
        gam = _gamma(decl)
        chi = np.linalg.solve(np.eye(ND * ND)[None] + self.chi0so @ gam, self.chi0so)
        # the dressing is not negligible: the check covers the fluctuation part
        self.assertGreater(np.linalg.norm(chi - self.chi0so), 1e-2 * np.linalg.norm(self.chi0so))
        Kso, _ = sc._make_kernel_operator(so_pairing_vertex(gam, chi, ND, (NX, 1, 1)),
                                          self.G2so, ND, NX, 1, 1)
        rng = np.random.default_rng(2)
        for eta, spin in components:
            V = sc._compute_vertices_general(self.chi0q, inter_k, NP, NX, 1, 1, 1,
                                             pairing_type=eta, static_index=0)
            K, _ = sc._make_kernel_operator(V, self.G2, NP, NX, 1, 1)
            phi = sc._project_gap_parity(
                rng.standard_normal((NP, NP, NX, 1, 1))
                + 1j * rng.standard_normal((NP, NP, NX, 1, 1)), eta)
            out = K.matvec(phi.ravel()).reshape(NP, NP, NX)
            want = np.einsum("st,abk->satbk", _SPIN[spin], out).reshape(ND, ND, NX)
            gap = np.einsum("st,abk->satbk", _SPIN[spin], phi[:, :, :, 0, 0])
            got = Kso.matvec(gap.ravel()).reshape(ND, ND, NX)
            with self.subTest(eta=eta, spin=spin):
                self.assertGreater(np.linalg.norm(want), 1e-3 * np.linalg.norm(phi))
                self.assertLess(np.linalg.norm(got - want), 1e-12 * np.linalg.norm(want))

    def _check_su2(self, decl):
        self._check(decl, [("singlet", "s")] + [("triplet", c) for c in "zxy"])

    def test_reproduces_singlet_and_triplet_kernels_coulomb_intra(self):
        self._check_su2(_onsite_decl("CoulombIntra", 1.2))

    def test_reproduces_singlet_and_triplet_kernels_coulomb_inter(self):
        self._check_su2(_onsite_decl("CoulombInter", 0.6))

    def test_reproduces_singlet_and_triplet_kernels_pairhop(self):
        self._check_su2(_onsite_decl("PairHop", 0.5))

    def test_reproduces_singlet_and_triplet_kernels_hund_plus_exchange(self):
        """Hund and Exchange at equal J: the rotation-invariant Kanamori
        Hund coupling."""
        decl = dict(_onsite_decl("Hund", 0.3))
        decl.update(_onsite_decl("Exchange", 0.3))
        self._check_su2(decl)

    def test_reproduces_singlet_and_triplet_kernels_all(self):
        decl = {}
        for t, v in (("CoulombIntra", 1.2), ("CoulombInter", 0.6), ("Hund", 0.25),
                     ("Exchange", 0.25), ("PairHop", 0.2)):
            decl.update(_onsite_decl(t, v))
        self._check_su2(decl)

    def test_equal_spin_triplet_of_hund_and_exchange_alone(self):
        """Hund or Exchange alone conserves S_z only: the equal-spin triplet
        (``d_x``, ``d_y``) -- the non-SO triplet's definition -- is still an
        exact limit."""
        for t, v in (("Hund", 0.3), ("Exchange", 0.25)):
            with self.subTest(t=t):
                self._check(_onsite_decl(t, v), [("triplet", "x"), ("triplet", "y")])

    def test_hund_alone_has_no_opposite_spin_pairing(self):
        """Hund alone couples equal spins only, so the spin-orbital kernel maps
        the singlet and ``d_z`` gaps (opposite-spin blocks) to exactly zero,
        while the equal-spin ``d_x`` gap gets a nonzero response. The
        non-spin-orbital singlet kernel is not compared here (issue #88)."""
        from hwave.solver.so_pairing import so_pairing_vertex
        gam = _gamma(_onsite_decl("Hund", 0.3))
        chi = np.linalg.solve(np.eye(ND * ND)[None] + self.chi0so @ gam, self.chi0so)
        Kso, _ = sc._make_kernel_operator(so_pairing_vertex(gam, chi, ND, (NX, 1, 1)),
                                          self.G2so, ND, NX, 1, 1)
        rng = np.random.default_rng(5)
        for spin in "szx":
            phi = rng.standard_normal((NP, NP, NX)) + 1j * rng.standard_normal((NP, NP, NX))
            gap = np.einsum("st,abk->satbk", _SPIN[spin], phi)
            out = Kso.matvec(gap.ravel()).reshape(2, NP, 2, NP, NX)
            with self.subTest(spin=spin):
                if spin == "x":
                    self.assertGreater(np.abs(out[0, :, 0]).max(), 1e-3)
                else:
                    self.assertEqual(np.abs(out).max(), 0.0)


_SPIN = {"s": _IY, "z": _TRIPLET["z"], "x": _TRIPLET["x"], "y": _TRIPLET["y"]}


class TestInterface(unittest.TestCase):

    def test_parts_sum_and_layout(self):
        from hwave.solver.so_pairing import so_pairing_vertex, so_pairing_vertex_parts
        rng = np.random.default_rng(4)
        gam = _gamma(_onsite_decl("CoulombIntra", 1.0))
        chi = rng.standard_normal((NX, ND * ND, ND * ND)) + 0j
        Vb, Vf = so_pairing_vertex_parts(gam, chi, ND, (NX, 1, 1))
        np.testing.assert_array_equal(so_pairing_vertex(gam, chi, ND, (NX, 1, 1)), Vb + Vf)
        # fluctuation part: -Gamma chi Gamma read in the kernel's (a, b, c, d)
        want = -(gam @ chi @ gam).reshape(NX, ND, ND, ND, ND)
        np.testing.assert_allclose(np.moveaxis(Vf[..., 0, 0], -1, 0), want, atol=1e-14)

    def test_rejects_bad_shapes(self):
        from hwave.solver.so_pairing import so_pairing_vertex
        gam = np.zeros((NX, ND * ND, ND * ND), complex)
        with self.assertRaises(ValueError):
            so_pairing_vertex(gam, gam[:, :1], ND, (NX, 1, 1))
        with self.assertRaises(ValueError):
            so_pairing_vertex(gam, gam, ND, (NX + 1, 1, 1))
        with self.assertRaises(ValueError):
            so_pairing_vertex(gam, gam, ND + 1, (NX, 1, 1))

    def test_rejects_momentum_dependent_vertex(self):
        """An off-site interaction makes Gamma depend on q; its exchange
        crossing then depends on k + k' and is not representable here."""
        from hwave.solver.so_pairing import so_pairing_vertex
        gam = _gamma(_onsite_decl("CoulombIntra", 1.0))
        gam = gam * (1.0 + 0.1 * np.cos(_K))[:, None, None]
        with self.assertRaisesRegex(ValueError, "on-site"):
            so_pairing_vertex(gam, np.zeros_like(gam), ND, (NX, 1, 1))


class TestSharedVertex(unittest.TestCase):
    """Delegation pin, not an oracle: the antisymmetrized spinful vertex has
    ONE builder, shared by the spin-orbital FLEX and the Eliashberg solver,
    so the two cannot drift apart. Both sides call the same code, so this
    test says nothing about whether that vertex is correct; its correctness
    is checked by the exact-diagonalization gate
    (``tests/test_so_pairing_ed.py``) and the SU(2) limits."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="so_vertex_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_matches_flex_spinful_vertex(self):
        s, _ = fso._build(self.d, fso._write_model(self.d, True, soc=0.2))
        nvol, nd = s.lattice.nvol, 2 * s.norb
        ham = np.asarray(s.ham_info.ham_inter_q).reshape(nvol, nd, nd, nd, nd)
        got = rpa.spinful_bubble_vertex(ham, s.ham_info.spinful_exchange())
        want = np.asarray(s._spinful_vertex(np))
        self.assertEqual(got.shape, (nvol, nd * nd, nd * nd))
        np.testing.assert_array_equal(got, want)


if __name__ == "__main__":
    unittest.main()
