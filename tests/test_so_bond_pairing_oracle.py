"""The static pairing vertex of the spinful bond vertex
(``so_bond.so_bond_pairing_vertex_parts``) in the bond pairing kernel.

From the code-slot vertex ``Gamma_B`` and its dressed susceptibility
``chi_B = [1 + chibar_B Gamma_B]^{-1} chibar_B`` the parts are

    V_bare = 0.25 (Gamma_B + Gamma_B^P),     V_fluct(q) = -Gamma_B(q) chi_B(q) Gamma_B(q),

both ``(nvol, ND, ND)`` in the block layout of
``eliashberg_bond.BondPairKernel`` with ``norb = nd = 2 * norb_phys``
spin-orbital modes and no instantaneous vertex.

* ``TestFirstOrderSOBondPairing``: the bare kernel, and the O(coupling)
  coefficient of the full static kernel, against the independent first-order
  construction of ``tests/test_sc_bond_pairing_oracle.py`` written out in the
  spin-orbital basis of the chain (complex spin-flip hopping, every off-site
  type at both orientations, a general antisymmetric gap over all spin
  blocks). The fluctuation part starts at second order and cancels in the
  central difference.
* ``TestSU2ReductionOfVFluct``: with spin-independent hopping and a
  spin-rotation invariant interaction on mixed bonds, the fluctuation kernel
  on ``phi (x) i sigma_y`` and ``phi (x) sigma_z i sigma_y`` equals the
  non-spin-orbital singlet and triplet fluctuation kernels built from the
  same bubble.
"""
import dataclasses
import tracemalloc
import unittest

import numpy as np

import hwave.sc as sc
import tests.test_sc_bond_pairing_oracle as O
from hwave.solver import bond_channels as bc
from hwave.solver.so_bond import (SLOT_CONVENTION, BondSlotGeometry, SOBondCaps, gamma_p,
                                  plan_so_bond, slot_geometry, so_bond_bubble, so_bond_dress,
                                  so_bond_pairing_vertex_parts, so_bond_topology,
                                  so_bond_vertex)
from hwave.solver.so_records import compile_interaction
from tests.heavy_tests import heavy
from tests.so_bond_fixtures import exact, green_kw, pair_kernel
from tests.so_bond_oracles import CELL, build_terms, closed_row, mesh_negate, rel

NX, BETA, NMAT, MU = O.NX, O.BETA, O.NMAT, O.MU
SHAPE = (NX, 1, 1)
NOCAP = SOBondCaps(None, None, None)
_K = O._K
_WN = O._WN

#: the off-site types and their coefficients (PairHop complex)
OFFSITE_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift", "PairHop")
COEFF = dict(zip(OFFSITE_TYPES, (0.37, 0.53, 0.71, 0.29, 0.61, 0.43 + 0.27j)))
U_INTRA = 0.83
#: the declared off-site row (a, b) per orbital count
ROW = {1: (0, 0), 2: (0, 1), 3: (0, 2)}
#: first-order tolerances: the bare kernel at every coupling (matvec round-off
#: through the pair-bubble sum) and the extrapolated central difference (its
#: O(h^4) remainder)
BARE_TOL = 1e-10
RICHARDSON_TOL = 1e-8
SU2_TOL = 1e-12

#: in-cell three-orbital loop 0 -> 1 -> 2 -> 0 with complex amplitudes, plus
#: an intra-orbital chain and on-site energies (rows (r, a, b): t c+_{j,a} c_{j+r,b})
_LOOP_3 = {(0, 0, 1): 0.3 + 0.2j, (0, 1, 2): 0.25 - 0.15j, (0, 2, 0): -0.2 + 0.3j}
TRANSFER_3 = {(0, 0, 0): -0.2, (0, 1, 1): 0.3, (0, 2, 2): 0.1,
              (1, 0, 0): -1.0, (-1, 0, 0): -1.0, (1, 1, 1): -0.7, (-1, 1, 1): -0.7,
              (1, 2, 2): -0.5, (-1, 2, 2): -0.5, (1, 0, 2): -0.4, (-1, 2, 0): -0.4}
for (_r, _a, _b), _t in _LOOP_3.items():
    TRANSFER_3[(_r, _a, _b)] = _t
    TRANSFER_3[(-_r, _b, _a)] = np.conj(_t)
TRANSFER = {1: O.TRANSFER[1], 2: O.TRANSFER[2], 3: TRANSFER_3}

#: Pauli structure of the singlet and the S_z = 0 triplet gap
ISY = np.array([[0, 1], [-1, 0]], dtype=complex)
SZ_ISY = np.diag([1.0, -1.0]).astype(complex) @ ISY


def _h_orbital(norb):
    h = np.zeros((NX, norb, norb), complex)
    for (r, a, b), t in TRANSFER[norb].items():
        h[:, a, b] += t * np.exp(1j * _K * r)
    return h


def _h_spin_orbital(norb, soc):
    """The chain in the spin-orbital basis (spin-block order ``s * norb + a``)
    plus a complex spin-flip hopping of strength ``soc``."""
    nd = 2 * norb
    h0 = _h_orbital(norb)
    h = np.zeros((NX, nd, nd), complex)
    h[:, :norb, :norb] = h0
    h[:, norb:, norb:] = h0
    rng = np.random.default_rng(11)
    for (r, a, b) in sorted({(1, 0, norb - 1), (1, 0, 0), (0, 0, norb - 1), (1, norb - 1, 0)}):
        lam = soc * (rng.standard_normal() + 1j * rng.standard_normal())
        hop = np.zeros((nd, nd), complex)
        hop[a, norb + b] = lam                     # c+_{j, a, up} c_{j+r, b, dn}
        hk = hop[None] * np.exp(1j * _K * r)[:, None, None]
        h += hk + hk.conj().transpose(0, 2, 1)
    return h


def _green(h):
    """``G(k, i w_n)`` in the sc.py layout ``(n, n, NX, 1, 1, NMAT)``."""
    n = h.shape[1]
    g = np.empty((n, n, NX, 1, 1, NMAT), complex)
    for i, w in enumerate(_WN):
        g[:, :, :, 0, 0, i] = np.linalg.inv((1j * w + MU) * np.eye(n)[None] - h).transpose(1, 2, 0)
    return g


_GREEN = {}


def _green_so(norb):
    """The spin-orbital Green function with spin-flip hopping (cached, read-only)."""
    if norb not in _GREEN:
        g = _green(_h_spin_orbital(norb, 0.3))
        g.flags.writeable = False
        _GREEN[norb] = g
    return _GREEN[norb]


def _tables(norb, t, R, scale):
    """The Hermitian-closed declaration of type ``t`` on the row
    ``(R, a, b)`` (``CoulombIntra``: every orbital on-site) at
    ``scale`` times its coefficient."""
    if t == "CoulombIntra":
        return {t: {((0, 0, 0), (a, a)): U_INTRA * scale for a in range(norb)}}
    a, b = ROW[norb]
    return {t: closed_row(R, a, b, COEFF[t] * scale)}


def _reference_gap_so(A, X, norb):
    """``Delta_AB(k)`` (spin-block, ``(nd, nd, NX)``) of the first-order gap
    equation ``D_pq = 1/2 sum_rs A_pq,rs kappa_sr`` for the antisymmetric
    spin-orbital pair amplitude ``X_AB(k)``; ``A`` is over the interleaved
    chain modes ``nd j + 2 a + s``."""
    nd = 2 * norb
    perm = np.array([2 * (A_ % norb) + A_ // norb for A_ in range(nd)])
    M = NX * nd
    kap = np.zeros((M, M), complex)
    for j in range(NX):
        for jp in range(NX):
            blk = (X * np.exp(-1j * _K * (jp - j))[None, None]).mean(-1)
            kap[np.ix_(j * nd + perm, jp * nd + perm)] = blk
    D = 0.5 * np.einsum("pqrs,sr->pq", A, kap)
    Dr = np.stack([D[np.ix_(perm, d * nd + perm)] for d in range(NX)], -1)
    return np.einsum("abd,kd->abk", Dr, np.exp(1j * np.outer(_K, np.arange(NX))))


def _antisymmetric_gap(rng, nd):
    """A random gap with ``Delta_AB(k) = -Delta_BA(-k)``, ``(nd, nd, NX)``."""
    phi = rng.standard_normal((nd, nd, NX)) + 1j * rng.standard_normal((nd, nd, NX))
    mk = (-np.arange(NX)) % NX
    return 0.5 * (phi - phi[:, :, mk].transpose(1, 0, 2))


def _reference(norb, t, R):
    """``(want, scale, phi)``: the first-order coefficient of the gap equation
    for the unit-coefficient declaration, on the seed-7 gap."""
    nd = 2 * norb
    g = _green_so(norb)
    A = O._antisymmetrized(build_terms(_tables(norb, t, R, 1.0), NX, norb), NX * nd)
    phi = _antisymmetric_gap(np.random.default_rng(7), nd)
    X = O._static_x(g, phi)
    want = _reference_gap_so(A, X, norb)
    return want, max(np.linalg.norm(want), np.linalg.norm(X)), phi


def _static_kernel(V, green_sc, view, nd, shape, beta):
    """The static bond pairing kernel of the code-slot vertex ``V``
    ``(nvol, B nd^2, B nd^2)`` on a one-point frequency axis, no
    instantaneous vertex."""
    G2 = bc._g2_from_green(green_sc, beta, tail=False)
    return pair_kernel(V[None], G2[..., None], view, nd, shape, beta)


def _apply(kernel, nd):
    return lambda phi: kernel.matvec(np.asarray(phi).ravel()).reshape(nd, nd, NX)


def _chain(tables, norb, g):
    """``(Gamma_B, chi_B, geom, topo, plan)`` of the production chain:
    ``so_bond_vertex``, the static ``so_bond_bubble`` of ``g`` and its
    dressing ``so_bond_dress``."""
    nd = 2 * norb
    c = compile_interaction(tables, norb)
    topo = so_bond_topology(c, CELL)
    geom = slot_geometry(topo, SHAPE)
    ND = geom.B * nd * nd
    plan = plan_so_bond(nmat=NMAT, nvol=NX, B=geom.B, nd=nd, consumers=("pairing_static",),
                        caps=NOCAP)
    gb = so_bond_vertex(c, topo, SHAPE, plan=plan)
    chibar = so_bond_bubble(green_kw(g), None, BETA, topo, SHAPE, plan=plan, freq="static",
                            sink=np.empty((NX, ND, ND), dtype=np.complex128))
    chi = so_bond_dress(chibar[None], gb, plan=plan, l0=NMAT // 2, spatial_shape=SHAPE)[0]
    return gb, chi, chibar, geom, topo, plan


def _production(norb, t, R, J):
    """``(K_pp, K)``: the bare and the full static kernel of the declaration
    at ``J`` times its coefficient."""
    nd = 2 * norb
    g = _green_so(norb)
    gb, chi, _, geom, topo, plan = _chain(_tables(norb, t, R, J), norb, g)
    vb, vf = so_bond_pairing_vertex_parts(gb, chi, geom, plan=plan)
    view = bc.BondSetView(topo)
    return (_apply(_static_kernel(vb, g, view, nd, SHAPE, BETA), nd),
            _apply(_static_kernel(vb + vf, g, view, nd, SHAPE, BETA), nd))


class TestFirstOrderSOBondPairing(unittest.TestCase):
    """Gate A1-O: the bare spinful bond pairing
    kernel and the O(coupling) coefficient of the full static kernel against
    the independent first-order spin-orbital construction (tolerances 1e-10
    and 1e-8). The fluctuation part is O(J^2) and cancels in the central
    difference; it is gated by ``TestSU2ReductionOfVFluct``."""

    H = 1.0e-3

    def _check(self, norb, t, R):
        want, scale, phi = _reference(norb, t, R)
        self.assertGreater(np.linalg.norm(want), 1e-2 * scale)          # not vacuous
        full = {}
        for J in (self.H, -self.H, self.H / 2, -self.H / 2):
            K_pp, K = _production(norb, t, R, J)
            err = np.linalg.norm(K_pp(phi) / J - want) / scale
            self.assertLess(err, BARE_TOL, (norb, t, R, J, "bare", err))
            full[J] = K(phi)
        d1 = (full[self.H] - full[-self.H]) / (2 * self.H)
        d2 = (full[self.H / 2] - full[-self.H / 2]) / self.H
        err = np.linalg.norm((4 * d2 - d1) / 3 - want) / scale
        self.assertLess(err, RICHARDSON_TOL, (norb, t, R, "full", err))

    def test_norb_2_fast(self):
        for t in ("Exchange", "PairLift", "PairHop"):
            for R in (+1, -1):
                with self.subTest(t=t, R=R):
                    self._check(2, t, R)
        with self.subTest(t="CoulombIntra"):
            self._check(2, "CoulombIntra", 0)

    @heavy
    def test_all_types_norb_1_2_3(self):
        for norb in (1, 2, 3):
            for t in OFFSITE_TYPES:
                for R in (+1, -1):
                    with self.subTest(norb=norb, t=t, R=R):
                        self._check(norb, t, R)
            with self.subTest(norb=norb, t="CoulombIntra"):
                self._check(norb, "CoulombIntra", 0)


# ---------------------------------------------------------------------------
# The parts helper: bounded memory and refusals
# ---------------------------------------------------------------------------

#: a spin-rotation invariant interaction on mixed bonds (norb_phys = 2):
#: inter-orbital CoulombInter different at +R and -R, Hund + Ising + Exchange
#: with J_Exchange = J_Hund - 2 J_Ising, real PairHop, on-site CoulombIntra
SU2_MIXED = {
    "CoulombIntra": {((0, 0, 0), (0, 0)): 0.83, ((0, 0, 0), (1, 1)): 0.83},
    "CoulombInter": {**closed_row(1, 0, 1, 0.37), **closed_row(-1, 0, 1, 0.22)},
    "Hund": closed_row(1, 0, 1, 0.53),
    "Ising": closed_row(1, 0, 1, 0.11),
    "Exchange": closed_row(1, 0, 1, 0.53 - 2 * 0.11),
    "PairHop": closed_row(1, 0, 1, 0.43),
}
_TOPO_6 = so_bond_topology(compile_interaction(SU2_MIXED, 2), CELL)
#: the channels 0, +-1 on the 6-cell chain (B = 3, nd = 4: ND = 48)
GEOM_6 = slot_geometry(_TOPO_6, SHAPE)
#: the channels 0, +-1, +-2 on the same chain (B = 5)
GEOM_6_B5 = BondSlotGeometry.from_arrays(
    np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0], [2, 0, 0], [-2, 0, 0]]),
    np.array([0, 2, 1, 4, 3]), SHAPE)


def _plan():
    return plan_so_bond(nmat=NMAT, nvol=NX, B=3, nd=4, consumers=("pairing_static",),
                        caps=NOCAP)


def _gamma_b_fixture():
    return so_bond_vertex(compile_interaction(SU2_MIXED, 2), _TOPO_6, SHAPE, plan=_plan())


def _chi_fixture():
    rng = np.random.default_rng(9)
    return rng.standard_normal((NX, 48, 48)) + 1j * rng.standard_normal((NX, 48, 48))


class TestPairingParts(unittest.TestCase):
    """The parts helper: the formulas, bounded peak memory, inputs left
    untouched, and its refusals before any allocation."""

    def test_parts_are_the_formulas(self):
        g, chi = _gamma_b_fixture(), _chi_fixture()
        vb, vf = so_bond_pairing_vertex_parts(g, chi, GEOM_6, plan=_plan())
        self.assertEqual((vb.shape, vb.dtype, vf.shape, vf.dtype),
                         (g.shape, np.complex128, g.shape, np.complex128))
        gp = gamma_p(g, GEOM_6, out=np.empty_like(g))
        self.assertLessEqual(rel(vb, 0.25 * (g + gp)), 1e-15)
        self.assertLessEqual(rel(vf, -(g @ chi @ g)), 1e-14)

    def test_parts_peak_is_bounded(self):
        g, chi = _gamma_b_fixture(), _chi_fixture()
        g0, chi0 = g.copy(), chi.copy()
        plan = _plan()
        tracemalloc.start()
        self.addCleanup(tracemalloc.stop)
        vb, vf = so_bond_pairing_vertex_parts(g, chi, GEOM_6, plan=plan)
        peak = tracemalloc.get_traced_memory()[1]
        self.assertLessEqual(peak, vb.nbytes + vf.nbytes + 2 * 16 * 48 * 48 + 16 * 1024)
        self.assertTrue(np.array_equal(g, g0))
        self.assertTrue(np.array_equal(chi, chi0))

    def test_mismatched_geometry_and_aliasing_refused(self):
        g, chi = _gamma_b_fixture(), _chi_fixture()
        plan = _plan()
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the admission plan has B = 3 but the call has "
                "B = 5 ({})".format(plan.reason))):
            so_bond_pairing_vertex_parts(g, chi, GEOM_6_B5, plan=plan)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: chi_b shares memory with gamma_b")):
            so_bond_pairing_vertex_parts(g, g, GEOM_6, plan=plan)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: chi_b shares memory with gamma_b")):
            so_bond_pairing_vertex_parts(g, g[:, :, ::-1], GEOM_6, plan=plan)
        with self.assertRaisesRegex(ValueError, exact("gamma_p: out shares memory with gamma_b")):
            gamma_p(g, GEOM_6, out=g)
        with self.assertRaisesRegex(ValueError, exact("gamma_p: out shares memory with gamma_b")):
            gamma_p(g, GEOM_6, out=g[:, :, ::-1])
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: chi_b must be a complex128 array of shape "
                "(6, 48, 48), got complex128 (6, 48, 47)")):
            so_bond_pairing_vertex_parts(g, chi[:, :, :47], GEOM_6, plan=plan)

    def test_non_default_convention_refused(self):
        conv = dataclasses.replace(SLOT_CONVENTION, phase_row=0)
        geom = BondSlotGeometry.from_arrays(GEOM_6.delta_r, GEOM_6.reverse, SHAPE,
                                            convention=conv)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: slot convention "
                "sr-:sc-:q-:nu+:rowrev1:colrev0:phr0:phc0 is not the pinned "
                "sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0")):
            so_bond_pairing_vertex_parts(_gamma_b_fixture(), _chi_fixture(), geom, plan=_plan())

    def test_gamma_p_refusal(self):
        g = _gamma_b_fixture()
        g[1, 0, 1] += 1e-6
        with self.assertRaisesRegex(ValueError, r"^so_bond_pairing_vertex_parts: Gamma_B violates "
                                                r"the particle-exchange symmetry \(\|Gamma\^P - "
                                                r"Gamma\| = \d\.\d{3}e[-+]\d{2} > 1e-13 "
                                                r"\|Gamma\|\)$"):
            so_bond_pairing_vertex_parts(g, np.zeros_like(g), GEOM_6, plan=_plan())

    def test_non_finite_inputs_refused_before_allocation(self):
        """A NaN or inf in ``gamma_b`` or ``chi_b`` is refused, naming the
        argument, before the outputs are allocated (the traced peak stays
        below one output array)."""
        plan = _plan()
        for arg in ("gamma_b", "chi_b"):
            for bad in (np.nan, np.inf, -np.inf):
                with self.subTest(arg=arg, bad=bad):
                    g, chi = _gamma_b_fixture(), _chi_fixture()
                    (g if arg == "gamma_b" else chi)[1, 0, 1] = bad
                    tracemalloc.start()
                    try:
                        base = tracemalloc.get_traced_memory()[0]
                        with self.assertRaisesRegex(ValueError, exact(
                                "so_bond_pairing_vertex_parts: {} has non-finite entries"
                                .format(arg))):
                            so_bond_pairing_vertex_parts(g, chi, GEOM_6, plan=plan)
                        peak = tracemalloc.get_traced_memory()[1] - base
                    finally:
                        tracemalloc.stop()
                    self.assertLess(peak, g.nbytes)


# ---------------------------------------------------------------------------
# SU(2) reduction of the fluctuation part
# ---------------------------------------------------------------------------

NP = 2


def _pair_field(seed, eta):
    """A random ``(NP, NP, NX)`` gap of the parity of ``eta``
    (``phi_ab(k) = +-phi_ba(-k)``)."""
    return O._pair_field(np.random.default_rng(seed), NP, eta)


def _spin_gap(phi, spin):
    """``Delta_{(s,a),(t,b)}(k) = spin_st phi_ab(k)``, spin-block ``(2 NP, 2 NP, NX)``."""
    return np.einsum("st,abk->satbk", spin, phi).reshape(2 * NP, 2 * NP, NX)


def _fluct_kernels(tables, gamma_q_flip=False):
    """``(so_fl, nonso_fl, pre)``: the spin-orbital fluctuation kernel
    (``V_fluct`` of the parts helper; with ``gamma_q_flip`` the outer
    vertices are read at ``-q``), the non-spin-orbital fluctuation kernels
    ``{eta: A_fl}`` (``sc._bond_kernel_parts``) built from the same bubble,
    and the preconditions of the comparison for the caller to assert:
    ``non_su2_bonds`` (must be empty), ``same_channels`` (the two topologies
    agree) and ``bubble_rel`` (the spin-up block of the spinful bubble
    against the orbital bubble)."""
    g = _green(_h_orbital(NP))
    gso = np.zeros((2 * NP, 2 * NP) + g.shape[2:], complex)
    gso[:NP, :NP] = g
    gso[NP:, NP:] = g
    gb, chi, chibar, geom, topo, plan = _chain(tables, NP, gso)
    if gamma_q_flip:
        gf = gb[mesh_negate(SHAPE)]
        vf = -(gf @ chi @ gf)
    else:
        _, vf = so_bond_pairing_vertex_parts(gb, chi, geom, plan=plan)
    so_fl = _apply(_static_kernel(vf, gso, bc.BondSetView(topo), 2 * NP, SHAPE, BETA), 2 * NP)

    topo_n, view_n, split = sc._resolve_bond_topology(tables, {"rvec": np.eye(3)}, NP, None)
    pre = {"non_su2_bonds": bc.non_su2_bonds(topo_n),
           "same_channels": np.array_equal(np.asarray(topo_n.delta_r), np.asarray(topo.delta_r))}
    S, C = sc._build_bond_vertices(topo_n, split, NP, SHAPE)
    B = geom.B
    # the spin-up block of the spinful code bubble is the orbital bubble
    t = chibar.reshape(NX, B, 2, NP, 2, NP, B, 2, NP, 2, NP)
    cb = np.ascontiguousarray(t[:, :, 0, :, 0, :, :, 0, :, 0, :]).reshape(NX, B * NP * NP, -1)
    pre["bubble_rel"] = rel(cb, bc.bond_bubble(g, view_n, BETA).reshape(cb.shape))
    dressed = [(ch, bc.dress_batch(cb[None], W, ch, l0=0, nmat=0, spatial_shape=SHAPE)[0][0])
               for ch, W in (("spin", S), ("charge", C))]
    nonso = {}
    for eta in ("singlet", "triplet"):
        _, A_fl, _, _ = sc._bond_kernel_parts(S, C, iter(dressed), g, view_n, NP, SHAPE, BETA,
                                              eta, topo=topo_n)
        nonso[eta] = (lambda op: lambda phi: op.matvec(phi.ravel()).reshape(NP, NP, NX))(A_fl)
    return so_fl, nonso, pre


class TestSU2ReductionOfVFluct(unittest.TestCase):
    """The SU(2) reduction of ``V_fluct``: with
    spin-independent hopping and a spin-rotation invariant interaction on
    mixed bonds (``bond_channels.non_su2_bonds == []``), the spinful
    fluctuation kernel on ``phi (x) i sigma_y`` and ``phi (x) sigma_z i
    sigma_y`` equals the non-spin-orbital singlet and triplet fluctuation
    kernels ``c_s S chi_s S + c_c C chi_c C`` dressed from the same bubble;
    reading the outer ``Gamma_B`` at ``-q`` misses."""

    def test_su2_reduction_of_v_fluct(self):
        so_fl, nonso_fl, pre = _fluct_kernels(SU2_MIXED)
        self.assertEqual(pre["non_su2_bonds"], [])
        self.assertTrue(pre["same_channels"])
        self.assertLessEqual(pre["bubble_rel"], 1e-13)
        res = {}
        for eta, spin in (("singlet", ISY), ("triplet", SZ_ISY)):
            phi = _pair_field(7, eta)
            want_orb = nonso_fl[eta](phi)
            self.assertGreater(np.abs(want_orb).max(), 1e-3 * np.abs(phi).max(), eta)
            res[eta] = rel(so_fl(_spin_gap(phi, spin)), _spin_gap(want_orb, spin))
        if max(res.values()) > SU2_TOL:
            self.fail("SU(2) reduction residuals (max-norm relative): {}".format(res))
        phi_s = _pair_field(7, "singlet")
        flipped = _fluct_kernels(SU2_MIXED, gamma_q_flip=True)[0]
        self.assertGreaterEqual(rel(flipped(_spin_gap(phi_s, ISY)), so_fl(_spin_gap(phi_s, ISY))),
                                1e-3)


if __name__ == "__main__":
    unittest.main()
