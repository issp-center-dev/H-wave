"""Gate C1b-SU2 (issue #270): the spin-orbital bond FLEX in the SU(2) limit
against the non-spin-orbital bond FLEX.

With a spin-independent band and a spin-rotation invariant interaction
(``J_Exchange = J_Hund - 2 J_Ising`` on every bond) the spinful map must
reduce to the non-spin-orbital one: ``Sigma_SO = Sigma_nonSO (x) 1_2``
(off-spin blocks exactly zero), and the spinful susceptibility projects onto
the non-spin-orbital charge and spin susceptibilities. Staged (R13):

* :class:`TestProjection` -- :func:`su2_project` on the code slots (R25),
  pinned on the bare bubble and on a dressed susceptibility of the on-site
  one-band Hubbard table (charge ``chibar / (1 + U chibar)``, spin
  ``chibar / (1 - U chibar)``); the label-form transverse pairing is
  identically zero.
* :class:`TestArraySU2` -- array level on ``SU2_MIXED`` (two orbitals,
  mixed bonds, a 6-cell chain, ``B = 3``): (i) the second order (both
  dressings truncated) and (ii) the full map of
  ``flex_so_bond.so_bond_flex_map`` against the non-spin-orbital
  production chain (bond bubble, ``flex_bond.dress_and_build_w`` with its
  ``bond_w2`` second order, ``calc_self_energy_bond``), both layouts; the
  projections of the full ``chi_B`` onto the non-spin-orbital ``chi_s`` /
  ``chi_c`` of every bond block.
* :class:`TestEntrySU2` -- the FLEX entry on ``SU2_4X4`` at a fixed chemical
  potential, one map: Sigma_HF, Sigma_fluct, G, the projections of the
  published chi0q / chiq; the content is off-site and carries the bond
  exchange; a truncated spin-orbital dressing misses (a broken ``P_B``
  misses at array level, two orbitals; at one orbital in the SU(2) limit it
  is invisible).
* :class:`TestConverged` -- (iii) the converged fixed points (heavy; its
  fast smoke is :meth:`TestEntrySU2.test_full_map`).

Both fixtures carry off-site Exchange and PairHop, so the non-spin-orbital
arm runs its ``bond_w2`` second order (asserted). Mutations patch module
attributes looked up at call time; every reference is computed before its
patch is entered.
"""
import math
import tempfile
import unittest
from typing import NamedTuple
from unittest import mock

import numpy as np

from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X
from tests.heavy_tests import heavy
from tests.so_bond_oracles import CELL, rel

#: the gate tolerance of the one-map and array comparisons
TOL = 1e-12
#: the tolerance of the projection pins (one dressing, no transport)
TOL_PROJ = 1e-13
#: the miss floor of a deliberate break
MISS = 1e-3
#: the margin by which a charge/spin swap must miss the projection pin
SWAP_MISS = 1e-2
#: the ring (third order and up) share of Sigma_fluct over its second order
RING_SHARE = 0.1
#: the off-site content floor of the non-spin-orbital Sigma_fluct
CONTENT = 0.1
#: the floor of the bond-exchange effect on Sigma_fluct
EXCHANGE_EFFECT = 1e-2
#: the SCF threshold of the converged comparison: max(10 eps, 1e-8), EPS = 12
CONV_TOL = max(10 * 1e-12, 1e-8)


# --------------------------------------------------------------------------
# The SU(2) projection of a spinful code-slot susceptibility
# --------------------------------------------------------------------------

class SU2Projection(NamedTuple):
    """The SU(2) channels of a spinful code-slot susceptibility, each
    ``lead + (B norb**2, B norb**2)`` on the pair slots ``(m; a, b)``."""
    charge: np.ndarray
    spin_z: np.ndarray
    transverse_ud: np.ndarray
    transverse_du: np.ndarray


def _spin_blocks(chi, norb):
    """``blk(s1, s2, s3, s4)``: the block ``chi[(m; s1 a, s2 b), (m'; s3 c,
    s4 d)]`` of a code-slot array (rank 4 ``(..., ND, ND)`` or the published
    rank-6 channel-0 ``(nmat, nvol, nd, nd, nd, nd)``), reshaped to ``lead +
    (B norb**2, B norb**2)``."""
    chi = np.asarray(chi)
    nd = 2 * int(norb)
    P = nd * nd
    if chi.ndim == 6:
        if chi.shape[2:] != (nd,) * 4:
            raise ValueError("su2_project: a rank-6 chi needs the trailing shape {}, got {}"
                             .format((nd,) * 4, chi.shape))
        # [l, q, a, c, b, d] = chi4[l, q, a nd + c, b nd + d] (R7)
        chi = chi.reshape(chi.shape[:2] + (P, P))
    ND = chi.shape[-1]
    if chi.ndim < 2 or chi.shape[-2] != ND or ND % P:
        raise ValueError("su2_project: chi must be (..., ND, ND) with ND a multiple of "
                         "(2 norb)**2 = {}, got {}".format(P, chi.shape))
    B = ND // P
    lead = chi.shape[:-2]
    t = chi.reshape(lead + (B, 2, norb, 2, norb, B, 2, norb, 2, norb))
    n = B * norb * norb

    def blk(s1, s2, s3, s4):
        return t[..., :, s1, :, s2, :, :, s3, :, s4, :].reshape(lead + (n, n))
    return blk


def su2_project(chi, norb):
    """The SU(2) projection of a spinful code-slot susceptibility (R25: the
    code-slot order as stored), ``sigma_up = +1``, ``sigma_dn = -1``:

    * ``charge[(m;a,b),(m';c,d)] = 1/2 sum_{s,s'} chi[(m; s a, s b), (m'; s' c, s' d)]``;
    * ``spin_z = 1/2 sum_{s,s'} sigma_s sigma_s' chi[...same...]``;
    * ``transverse_ud = chi[(m; up a, dn b), (m'; up c, dn d)]``,
      ``transverse_du = chi[(m; dn a, up b), (m'; dn c, up d)]``.

    ``chi``: rank 4 ``(..., ND, ND)`` code slots or the published rank-6
    channel-0 ``chiq`` ``(nmat, nvol, nd, nd, nd, nd)`` (reshaped to
    ``(nmat, nvol, nd**2, nd**2)`` first). Returns :class:`SU2Projection`."""
    blk = _spin_blocks(chi, norb)
    up, dn = 0, 1
    sig = {up: 1.0, dn: -1.0}
    charge = 0.5 * sum(blk(s, s, s2, s2) for s in (up, dn) for s2 in (up, dn))
    spin_z = 0.5 * sum(sig[s] * sig[s2] * blk(s, s, s2, s2) for s in (up, dn) for s2 in (up, dn))
    return SU2Projection(charge, spin_z, blk(up, dn, up, dn), blk(dn, up, dn, up))


def label_form_blocks(chi, norb):
    """The label-form transverse pairing ``chi[(up a, dn b), (dn c, up d)]``
    and ``chi[(dn a, up b), (up c, dn d)]``: identically zero in the code
    slots (R25)."""
    blk = _spin_blocks(chi, norb)
    return blk(0, 1, 1, 0), blk(1, 0, 0, 1)


# --------------------------------------------------------------------------
# Array-level construction (the chain of tests/test_sc_bond_pairing_oracle.py)
# --------------------------------------------------------------------------

NX, BETA, NMAT, MU = 6, 2.0, 16, 0.15
SHAPE = (NX, 1, 1)
NB = 4
#: Transfer rows (r, a, b) of the chain (tests/test_sc_bond_pairing_oracle.TRANSFER)
CHAIN_TRANSFER = {
    1: {(1, 0, 0): -1.0, (-1, 0, 0): -1.0, (0, 0, 0): 0.1},
    2: {(0, 0, 0): -0.4, (0, 1, 1): 0.6, (0, 0, 1): -0.3, (0, 1, 0): -0.3,
        (1, 0, 1): -1.0, (-1, 1, 0): -1.0, (1, 0, 0): -0.5, (-1, 0, 0): -0.5},
}
#: the on-site one-band coupling of the projection pin
U_ONSITE = 0.7


def chain_green(norb):
    """``(g, g_so)``: the bare chain Green function ``(nmat, nx, norb,
    norb)`` and its spinful spin-block form ``(nmat, nx, 2 norb, 2 norb)``
    (the same ``g`` on both spins, off-spin blocks zero)."""
    k = 2 * np.pi * np.arange(NX) / NX
    wn = (2 * np.arange(NMAT) + 1 - NMAT) * np.pi / BETA
    h = np.zeros((NX, norb, norb), complex)
    for (r, a, b), t in CHAIN_TRANSFER[norb].items():
        h[:, a, b] += t * np.exp(1j * k * r)
    g = np.empty((NMAT, NX, norb, norb), complex)
    for i, w in enumerate(wn):
        g[i] = np.linalg.inv((1j * w + MU) * np.eye(norb)[None] - h)
    g_so = np.zeros((NMAT, NX, 2 * norb, 2 * norb), complex)
    g_so[..., :norb, :norb] = g
    g_so[..., norb:, norb:] = g
    return g, g_so


def _transfer_table(norb):
    """The interleaved spinful Transfer table of the chain (spin-diagonal)."""
    return {((r, 0, 0), (2 * a + s, 2 * b + s)): t
            for (r, a, b), t in CHAIN_TRANSFER[norb].items() for s in (0, 1)}


class SOSide(NamedTuple):
    res: object          # flex_so_bond.SOBondFlexMap
    chibar: object       # the full dense chibar_B (w_in_place = False), else None
    plan: object
    gamma: object


def so_side(tables, norb, layout, *, w_in_place=None):
    """One production map :func:`flex_so_bond.so_bond_flex_map` of the chain
    (``layout``, ``nb = NB``, the dressing's default guard); with
    ``w_in_place = False`` (dense) the full ``chibar_B`` is read back from
    its own slot."""
    from hwave.solver import flex_so_bond, so_bond
    from hwave.solver.so_records import compile_interaction
    nd = 2 * norb
    c = compile_interaction(tables, norb)
    topo = so_bond.so_bond_topology(c, CELL)
    geom = so_bond.slot_geometry(topo, SHAPE)
    reason, sm = None, None
    if layout == "sz_packed":
        selected, reason = so_bond.select_layout(c, _transfer_table(norb), override="auto")
        assert selected == "sz_packed", reason
        sm = so_bond.sz_sector_map(topo, nd)
    plan = so_bond.plan_so_bond(nmat=NMAT, nvol=NX, B=geom.B, nd=nd, consumers=("flex_scf",),
                                caps=so_bond.SOBondCaps(None, None, None), layout=layout,
                                layout_reason=reason, nb=NB, w_in_place=w_in_place)
    gamma = so_bond.so_bond_vertex(c, topo, SHAPE, plan=plan)
    G = chain_green(norb)[1][None].copy()
    chibar = None
    with so_bond.open_store(plan, sm) as st:
        res = flex_so_bond.so_bond_flex_map(st, gamma, G, G, None, BETA, geom, plan=plan,
                                            cond_tol=so_bond.bond_channels._BOND_COND_FLOOR)
        if w_in_place is False:
            assert layout == "dense"
            chibar = np.array(st.get_freq_batch("chibar", 0, NMAT), copy=True)
    return SOSide(res, chibar, plan, gamma)


def so_dress_full(side):
    """The full dressed ``chi_B`` of a dense :func:`so_side` (``w_in_place
    = False``) by the production dressing, batch by batch."""
    from hwave.solver import so_bond
    return np.concatenate([
        so_bond.so_bond_dress(side.chibar[l0:l0 + NB], side.gamma, plan=side.plan, l0=l0,
                              spatial_shape=SHAPE)
        for l0 in range(0, NMAT, NB)])


class NonSOSide(NamedTuple):
    sigma: np.ndarray    # (nmat, nx, norb, norb)
    chibar: np.ndarray   # (nmat, nx, ND, ND)
    chi_s: np.ndarray
    chi_c: np.ndarray
    dress: object        # flex_bond.DressResult
    n_w2: int            # the number of bond_w2 factors


def nonso_side(tables, norb):
    """The non-spin-orbital production chain of the same tables: the bond
    view and vertices (``sc._resolve_bond_topology``,
    ``sc._build_bond_vertices``), the bond bubble, the dressing loop
    ``flex_bond.dress_and_build_w`` (its ring and the ``bond_w2`` second
    order) and ``calc_self_energy_bond``."""
    import hwave.sc as sc
    from hwave.solver import bond_channels as bc
    from hwave.solver import flex_bond
    from hwave.solver.second_order import build_bond_w2_factors
    topo, view, split = sc._resolve_bond_topology(tables, {"rvec": np.eye(3)}, norb, None)
    assert bc.non_su2_bonds(topo) == []
    S, C = sc._build_bond_vertices(topo, split, norb, SHAPE)
    ND = S.shape[-1]
    w2 = build_bond_w2_factors(split, norb, topo.delta_r, SHAPE)
    g = chain_green(norb)[0][None].copy()
    names = ("chibar", "W", "chi_s_w", "chi_c_w")
    with flex_bond.BondBlockStore(NMAT, NX, ND, norb * norb, names) as st, \
            flex_bond.BondDeviceContext.for_view(np, S, C, view, norb) as dev:
        flex_bond.assemble_bubble(st, g, None, BETA, view, SHAPE, 1)
        chibar = np.array(st.get_freq_batch("chibar", 0, NMAT), copy=True)
        res = flex_bond.dress_and_build_w(st, dev, nb=NB, output_full=True, nmat=NMAT,
                                          nvol=NX, nd=norb * norb, spatial_shape=SHAPE,
                                          bond_w2=w2)
        sigma = flex_bond.calc_self_energy_bond(st, g, BETA, view, SHAPE, norb, 1)[0]
        chi_s = np.array(st.get_freq_batch("chi_s_w", 0, NMAT), copy=True)
        chi_c = np.array(st.get_freq_batch("chi_c_w", 0, NMAT), copy=True)
    return NonSOSide(sigma, chibar, chi_s, chi_c, res, len(w2))


def nonso_bubble(tables, norb):
    """The non-spin-orbital bond bubble ``(nmat, nx, ND, ND)`` of the chain."""
    import hwave.sc as sc
    from hwave.solver import flex_bond
    _, view, _ = sc._resolve_bond_topology(tables, {"rvec": np.eye(3)}, norb, None)
    ND = int(view.n_channels) * norb * norb
    g = chain_green(norb)[0][None].copy()
    with flex_bond.BondBlockStore(NMAT, NX, ND, norb * norb, ("chibar",)) as st:
        flex_bond.assemble_bubble(st, g, None, BETA, view, SHAPE, 1)
        return st.detach("chibar")


def _dress_identity(X_, gamma_b, **kwargs):
    """``so_bond.so_bond_dress`` truncated to its input (chi := chibar)."""
    from hwave.solver import so_bond
    if isinstance(X_, so_bond.PackedArrays):
        return so_bond.PackedArrays(X_.sector_map, tuple(b.copy() for b in X_.blocks))
    return np.array(X_, copy=True)


def _dress_batch_identity(cb, V, channel, **kwargs):
    """``bond_channels.dress_batch`` truncated to its input (chi := chibar)."""
    return cb.copy(), 1.0


def _spin_rows(tc, sigma_so, ref, norb, label, tol=TOL):
    """Both spin blocks of ``sigma_so`` equal ``ref`` within ``tol``; the
    off-spin blocks are exactly zero. Returns the two deviations."""
    up = rel(sigma_so[..., :norb, :norb], ref)
    dn = rel(sigma_so[..., norb:, norb:], ref)
    for v in (up, dn):
        tc.assertTrue(np.isfinite(v), label)
    tc.assertGreater(float(np.abs(ref).max()), 0.0, label)
    tc.assertLessEqual(up, tol, label + " up-up")
    tc.assertLessEqual(dn, tol, label + " dn-dn")
    tc.assertFalse(np.any(sigma_so[..., :norb, norb:]), label + " up-dn")
    tc.assertFalse(np.any(sigma_so[..., norb:, :norb]), label + " dn-up")
    return up, dn


# --------------------------------------------------------------------------
# The projection pin
# --------------------------------------------------------------------------

class TestProjection(unittest.TestCase):
    """:func:`su2_project` on the bare and a dressed spinful bubble of the
    on-site one-band Hubbard chain (``B = 1``, production bubble and
    dressing)."""

    @classmethod
    def setUpClass(cls):
        from hwave.solver import flex_bond, so_bond
        from hwave.solver.so_records import compile_interaction
        tables = {"CoulombIntra": {((0, 0, 0), (0, 0)): U_ONSITE}}
        c = compile_interaction(tables, 1)
        topo = so_bond.so_bond_topology(c, CELL)
        B = so_bond.slot_geometry(topo, SHAPE).B
        assert B == 1
        plan = so_bond.plan_so_bond(nmat=NMAT, nvol=NX, B=B, nd=2, consumers=("chi_dynamic",),
                                    caps=so_bond.SOBondCaps(None, None, None), layout="dense")
        gamma = so_bond.so_bond_vertex(c, topo, SHAPE, plan=plan)
        g_so = chain_green(1)[1][None].copy()
        with flex_bond.BondBlockStore(NMAT, NX, plan.ND, 4, ("chibar",)) as st:
            so_bond.so_bond_bubble(g_so, None, BETA, topo, SHAPE, plan=plan, freq="all", sink=st)
            cls.chibar = st.detach("chibar")
        cls.chi = so_bond.so_bond_dress(cls.chibar, gamma, plan=plan, l0=0, spatial_shape=SHAPE)
        cls.x = nonso_bubble(tables, 1)                  # (nmat, nx, 1, 1)

    def test_bare_and_dressed(self):
        x = self.x
        rows = []
        for label, chi in (("bare", self.chibar), ("dressed", self.chi)):
            p = su2_project(chi, 1)
            if label == "bare":
                charge_ref = spin_ref = x
                swap = None
            else:
                charge_ref = x / (1 + U_ONSITE * x)
                spin_ref = x / (1 - U_ONSITE * x)
                swap = (rel(p.charge, spin_ref), rel(p.spin_z, charge_ref))
            dev = {"charge": rel(p.charge, charge_ref), "spin_z": rel(p.spin_z, spin_ref),
                   "transverse_ud": rel(p.transverse_ud, p.spin_z),
                   "transverse_du": rel(p.transverse_du, p.spin_z)}
            rows.append((label, dev, swap))
            with self.subTest(chi=label):
                for k, v in dev.items():
                    self.assertTrue(np.isfinite(v), k)
                    self.assertLessEqual(v, TOL_PROJ, k)
                self.assertGreater(float(np.abs(p.spin_z).max()), 0.0)
                if swap is not None:                     # a charge/spin swap misses
                    for v in swap:
                        self.assertTrue(np.isfinite(v))
                        self.assertGreaterEqual(v, SWAP_MISS)
                # R25: the label-form pairing is identically zero in code slots
                for blk in label_form_blocks(chi, 1):
                    self.assertFalse(np.any(blk))
                # the rank-6 channel-0 form gives the same projection
                nd = 2
                six = np.asarray(chi).reshape(NMAT, NX, nd, nd, nd, nd)
                for a, b in zip(su2_project(six, 1), p):
                    self.assertTrue(np.array_equal(a, b))
        print("\nC1b-SU2 projection pin:\n" + "\n".join(
            "  %-8s " % lab + "  ".join("%s %.2e" % kv for kv in dev.items())
            + ("" if sw is None else "  swap misses %.2e %.2e" % sw)
            for lab, dev, sw in rows))

    def test_refusals(self):
        with self.assertRaises(ValueError):
            su2_project(np.zeros((2, 3, 5, 5)), 1)
        with self.assertRaises(ValueError):
            su2_project(np.zeros((2, 3, 2, 2, 2, 1)), 1)


# --------------------------------------------------------------------------
# Array level: SU2_MIXED on the 6-cell chain
# --------------------------------------------------------------------------

class TestArraySU2(unittest.TestCase):
    """C1b-SU2 (i) and (ii) at array level on ``SU2_MIXED`` (norb 2, B = 3,
    Nmat 16)."""

    NORB = 2
    TABLES = X.SU2_MIXED

    @classmethod
    def setUpClass(cls):
        cls.nonso = nonso_side(cls.TABLES, cls.NORB)          # the references first
        cls.so = {lay: so_side(cls.TABLES, cls.NORB, lay) for lay in ("dense", "sz_packed")}
        cls.so_full = so_side(cls.TABLES, cls.NORB, "dense", w_in_place=False)

    def test_second_order(self):
        from hwave.solver import bond_channels as bc
        from hwave.solver import so_bond
        n = self.NORB
        so_spy = mock.Mock(side_effect=_dress_identity)
        nonso_spy = mock.Mock(side_effect=_dress_batch_identity)
        # the non-spin-orbital second order (its ring removed) first
        with mock.patch.object(bc, "dress_batch", nonso_spy):
            ref = nonso_side(self.TABLES, n)
        self.assertGreater(nonso_spy.call_count, 0)                # applied
        self.assertEqual(so_spy.call_count, 0)
        self.assertGreater(ref.n_w2, 0)                            # the bond_w2 branch
        rows = []
        for layout in ("dense", "sz_packed"):
            before = so_spy.call_count
            with mock.patch.object(so_bond, "so_bond_dress", so_spy):
                side = so_side(self.TABLES, n, layout)
            self.assertEqual(so_spy.call_count - before, math.ceil(NMAT / NB))   # applied
            with self.subTest(layout=layout):
                up, dn = _spin_rows(self, side.res.sigma_fluct[0], ref.sigma, n,
                                    "second order " + layout)
                rows.append((layout, up, dn))
        # non-vacuity: the ring share of the full Sigma_fluct over its second order
        share = float(np.abs(self.nonso.sigma - ref.sigma).max() / np.abs(ref.sigma).max())
        print("\nC1b-SU2 (i) second order (up-up, dn-dn): " + "  ".join(
            "%s %.2e %.2e" % r for r in rows) + "  ring share %.3f" % share)
        self.assertTrue(np.isfinite(share))
        self.assertGreaterEqual(share, RING_SHARE)

    def test_full_map(self):
        n = self.NORB
        ref = self.nonso
        self.assertGreater(ref.n_w2, 0)                            # the bond_w2 branch
        rows = []
        for layout, side in sorted(self.so.items()):
            with self.subTest(layout=layout):
                up, dn = _spin_rows(self, side.res.sigma_fluct[0], ref.sigma, n,
                                    "full map " + layout)
                rows.append((layout, up, dn))
        # the second-slot map (the projections' source) is the same map
        self.assertLessEqual(rel(self.so_full.res.sigma_fluct, self.so["dense"].res.sigma_fluct),
                             TOL)
        print("\nC1b-SU2 (ii) full map (up-up, dn-dn): " + "  ".join(
            "%s %.2e %.2e" % r for r in rows))

    def test_projections(self):
        n = self.NORB
        ref = self.nonso
        side = self.so_full
        chi = so_dress_full(side)
        P = (2 * n) ** 2
        rows = {}
        bare = su2_project(side.chibar, n)
        rows["bare charge"] = rel(bare.charge, ref.chibar)
        rows["bare spin_z"] = rel(bare.spin_z, ref.chibar)
        p = su2_project(chi, n)
        rows["charge"] = rel(p.charge, ref.chi_c)
        rows["spin_z"] = rel(p.spin_z, ref.chi_s)
        rows["transverse_ud"] = rel(p.transverse_ud, ref.chi_s)
        rows["transverse_du"] = rel(p.transverse_du, ref.chi_s)
        # the map's published channel-0 blocks project onto the non-SO collapses
        p0 = su2_project(side.res.chi_c0, n)
        rows["c0 chibar"] = rel(su2_project(side.res.chi0_c0, n).charge, ref.dress.collapse0)
        rows["c0 charge"] = rel(p0.charge, ref.dress.collapse_c)
        rows["c0 spin_z"] = rel(p0.spin_z, ref.dress.collapse_s)
        rows["c0 = chi_B channel 0"] = rel(side.res.chi_c0, chi[..., :P, :P])
        print("\nC1b-SU2 (ii) projections: " + "  ".join("%s %.2e" % kv for kv in rows.items()))
        for k, v in rows.items():
            self.assertTrue(np.isfinite(v), k)
            self.assertLessEqual(v, TOL, k)
        # the charge and spin channels really differ (a swap would miss)
        self.assertGreaterEqual(rel(p.charge, ref.chi_s), SWAP_MISS)
        for arr in (side.chibar, chi):
            for blk in label_form_blocks(arr, n):
                self.assertFalse(np.any(blk))

    def test_mutation(self):
        """``so_bond.p_b`` replaced by the identity (the non-SO ring's
        ``pin`` orientation is then no longer matched) misses."""
        from hwave.solver import so_bond
        from tests.test_flex_so_bond_entry import p_b_identity
        n = self.NORB
        ref = self.nonso.sigma                                  # the reference first
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                spy = mock.Mock(wraps=p_b_identity)
                self.assertEqual(spy.call_count, 0)
                with mock.patch.object(so_bond, "p_b", spy):
                    side = so_side(self.TABLES, n, layout)
                self.assertGreater(spy.call_count, 0)           # the break was applied
                miss = rel(side.res.sigma_fluct[0][..., :n, :n], ref)
                print("C1b-SU2 array mutation P_B identity (%s) miss %.3e" % (layout, miss))
                self.assertTrue(np.isfinite(miss))
                self.assertGreaterEqual(miss, MISS)


# --------------------------------------------------------------------------
# The FLEX entry: SU2_4X4
# --------------------------------------------------------------------------

FX = F.SU2_4X4
DENSE = {"so_bond_layout": "dense"}
#: the Matsubara count of the entry runs
NMAT_ENTRY = X.FLEX_NMAT
#: the interaction types of the bond exchange content (J_H = J_I = J_Ex = 0 without them)
EXCHANGE_TYPES = ("Hund", "Ising", "Exchange")


def _run(fixture, *, so, param=None):
    """The :class:`tests.so_bond_flex_fixtures.RunResult` of the bond-gate
    FLEX run of :func:`tests.so_bond_flex_fixtures.write_flex_case`
    (``so``; ``param`` overrides ``[mode.param]``)."""
    with tempfile.TemporaryDirectory(prefix="so_bond_su2_") as d:
        return X.run_input(X.write_flex_case(d, fixture, so=so, param=param))


def _one_map(fixture, *, so, param=None):
    res = _run(fixture, so=so, param=dict(param or {}, IterationMax=1))
    if len(res.payloads) != 1:
        raise RuntimeError("_one_map: {} payloads".format(len(res.payloads)))
    return res


def _converge(fixture, *, so, param=None):
    par = dict(X.CONVERGE_PARAM)
    par.update(param or {})
    res = _run(fixture, so=so, param=par)
    if not res.solver.scf_converged:
        raise RuntimeError("_converge: {} (so {}) did not converge in {} iterations".format(
            fixture.name, so, res.solver.scf_iterations))
    return res


def _nonso_checks(tc, res):
    """The non-spin-orbital arm: the bond gate with its bond_w2 second
    order, SU(2) bonds."""
    from hwave.solver import bond_channels as bc
    s = res.solver
    tc.assertFalse(getattr(s, "_so_bond", False))
    tc.assertIsNotNone(s._bond_w2)
    tc.assertGreater(len(s._bond_w2), 0)
    tc.assertEqual(bc.non_su2_bonds(s._bond_topo), [])


class TestEntrySU2(unittest.TestCase):
    """C1b-SU2 (i)-(ii) at the FLEX entry: one map of ``SU2_4X4`` at a fixed
    chemical potential, spin-orbital (packed and dense) against non-spin-orbital."""

    @classmethod
    def setUpClass(cls):
        cls.nonso = _one_map(FX, so=False)                     # the reference first
        cls.so = {"sz_packed": _one_map(FX, so=True), "dense": _one_map(FX, so=True, param=DENSE)}

    def test_full_map(self):
        _nonso_checks(self, self.nonso)
        q = self.nonso.payloads[0]
        g_n = self.nonso.green_info["green"]
        norb = FX.norb_phys
        rows = []
        for layout, res in sorted(self.so.items()):
            self.assertEqual(res.solver._so_admission.layout, layout)
            p = res.payloads[0]
            self.assertEqual(p["mu"], q["mu"])
            with self.subTest(layout=layout):
                dev = {}
                dev["Sigma_HF"] = _spin_rows(self, p["static_new"], q["static_new"], norb,
                                             "Sigma_HF")
                dev["Sigma_fluct"] = _spin_rows(self, p["fluct_new"], q["fluct_new"], norb,
                                                "Sigma_fluct")
                dev["G"] = _spin_rows(self, res.green_info["green"], g_n, norb, "G")
                pairs = q["chiq_s"].shape[2] ** 2
                cn = {k: q[k].reshape(NMAT_ENTRY, -1, pairs, pairs)
                      for k in ("chi0q", "chiq_s", "chiq_c")}
                bare = su2_project(p["chi0q"], norb)
                pr = su2_project(p["chiq_s"], norb)
                proj = {"chi0q charge": rel(bare.charge, cn["chi0q"]),
                        "chi0q spin_z": rel(bare.spin_z, cn["chi0q"]),
                        "chiq charge": rel(pr.charge, cn["chiq_c"]),
                        "chiq spin_z": rel(pr.spin_z, cn["chiq_s"]),
                        "chiq transverse_ud": rel(pr.transverse_ud, cn["chiq_s"]),
                        "chiq transverse_du": rel(pr.transverse_du, cn["chiq_s"])}
                for k, v in proj.items():
                    self.assertTrue(np.isfinite(v), k)
                    self.assertLessEqual(v, TOL, k)
                for blk in label_form_blocks(p["chiq_s"], norb):
                    self.assertFalse(np.any(blk))
                # the charge and spin channels really differ
                self.assertGreaterEqual(rel(cn["chiq_c"], cn["chiq_s"]), SWAP_MISS)
                rows.append((layout, dev, proj))
        print("\nC1b-SU2 entry one map:")
        for layout, dev, proj in rows:
            print("  %-9s " % layout + "  ".join("%s %.2e/%.2e" % (k, *v) for k, v in dev.items())
                  + "\n            " + "  ".join("%s %.2e" % kv for kv in proj.items()))

    def test_non_vacuity(self):
        q = self.nonso.payloads[0]
        fluct = q["fluct_new"]
        off_site = float(np.abs(fluct - fluct.mean(axis=-3, keepdims=True)).max()
                         / np.abs(fluct).max())
        no_ex = FX._replace(name="SU2_4X4_NO_EXCHANGE", interactions=tuple(
            (t, rows) for t, rows in FX.interactions if t not in EXCHANGE_TYPES))
        alt = _one_map(no_ex, so=False).payloads[0]
        effect = rel(alt["fluct_new"], fluct)
        effect_hf = rel(alt["static_new"], q["static_new"])
        print("\nC1b-SU2 non-vacuity: off-site share of Sigma_fluct %.3e; bond exchange "
              "effect on Sigma_fluct %.3e (Sigma_HF %.3e)" % (off_site, effect, effect_hf))
        for v in (off_site, effect):
            self.assertTrue(np.isfinite(v))
        self.assertGreaterEqual(off_site, CONTENT)
        self.assertGreaterEqual(effect, EXCHANGE_EFFECT)

    def test_mutation(self):
        """A truncated dressing (chi := chibar, the second order) in the
        spin-orbital run misses the non-spin-orbital map. (``P_B`` replaced
        by the identity is invisible here: at one physical orbital in the
        SU(2) limit it maps every block onto an equal one; its break row is
        :meth:`TestArraySU2.test_mutation`, two orbitals.)"""
        from hwave.solver import so_bond
        ref = self.nonso.payloads[0]["fluct_new"]               # the reference (non-SO) first
        norb = FX.norb_phys
        for layout, param in (("sz_packed", None), ("dense", DENSE)):
            with self.subTest(layout=layout):
                spy = mock.Mock(side_effect=_dress_identity)
                self.assertEqual(spy.call_count, 0)
                with mock.patch.object(so_bond, "so_bond_dress", spy):
                    res = _one_map(FX, so=True, param=param)
                self.assertEqual(res.solver._so_admission.layout, layout)
                self.assertGreater(spy.call_count, 0)           # the break was applied
                miss = rel(res.payloads[0]["fluct_new"][..., :norb, :norb], ref)
                print("C1b-SU2 entry mutation chi := chibar (%s) miss %.3e" % (layout, miss))
                self.assertTrue(np.isfinite(miss))
                self.assertGreaterEqual(miss, MISS)


class TestConverged(unittest.TestCase):
    """C1b-SU2 (iii): the converged fixed points (heavy; fast smoke
    :meth:`TestEntrySU2.test_full_map`)."""

    @heavy
    def test_fixed_points(self):
        n = _converge(FX, so=False)
        _nonso_checks(self, n)
        norb = FX.norb_phys
        b = n.green_info
        for layout, param in (("sz_packed", None), ("dense", DENSE)):
            with self.subTest(layout=layout):
                s = _converge(FX, so=True, param=param)
                self.assertEqual(s.solver._so_admission.layout, layout)
                a = s.green_info
                sig = _spin_rows(self, a["sigma"], b["sigma"], norb, "Sigma", tol=CONV_TOL)
                grn = _spin_rows(self, a["green"], b["green"], norb, "G", tol=CONV_TOL)
                mu = abs(a["physics"]["mu"] - b["physics"]["mu"])
                ncond = abs(a["physics"]["NCond"] - b["physics"]["NCond"])
                print("C1b-SU2 converged %s (%d / %d iterations; residuals sigma %.1e / %.1e) "
                      "Sigma %.2e/%.2e  G %.2e/%.2e  mu %.2e  NCond %.2e" % (
                          layout, s.solver.scf_iterations, n.solver.scf_iterations,
                          s.solver.scf_sigma_residual, n.solver.scf_sigma_residual,
                          *sig, *grn, mu, ncond))
                for v in (mu, ncond):
                    self.assertTrue(np.isfinite(v))
                    self.assertLessEqual(v, CONV_TOL)


if __name__ == "__main__":
    unittest.main()
