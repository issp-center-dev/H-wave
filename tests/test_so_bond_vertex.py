"""The spinful bond vertex ``Gamma_B`` of the spin-orbital interaction
(``so_bond.so_bond_vertex``), its exchange image ``so_bond.gamma_p`` and the
admission plan (``so_bond.plan_so_bond``).

``Gamma_B(q)`` is an ``ND x ND`` matrix per momentum, ``ND = B * nd**2``,
in the code slots of ``so_bond.SLOT_CONVENTION``:
``Gamma_B(q)[(m;x,y),(m';z,w)] = e^{+iq.dr_m'} Gamma(q)[(-dr_m; x,y),
(+dr_m'; w,z)]``, ``Gamma`` the antisymmetrized (direct plus crossed) label-form
interaction, ``H_int = (1/4N) sum_q Gamma_IJ(q) rho_I(q) rho_J(-q)`` up to
one-body terms.
"""
import dataclasses
import tracemalloc
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import flex_bond
from hwave.solver.so_bond import (SLOT_CONVENTION, BondSlotGeometry, PackedArrays, SOBondCaps,
                                  SOBondPlan, _file_records, estimate_packed_kernel_workspace_bytes,
                                  estimate_sector_solver_bytes, gamma_p, plan_so_bond,
                                  slot_geometry, so_bond_bubble, so_bond_dress,
                                  so_bond_pairing_vertex_parts, so_bond_topology, so_bond_vertex,
                                  sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import exact, onsite_decl, production_gamma_b, rpa_spinful_gamma
from tests.so_bond_oracles import (BASE, CELL, COMPILER_TABLES, TYPES, TYPES_OFFSITE, closed_row,
                                   gamma_code_to_label, gamma_label_to_code, mesh_momenta,
                                   mesh_negate, rel)

NOCAP = SOBondCaps(None, None, None)
SIGN = SLOT_CONVENTION.slot_sign
DR3 = np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]], dtype=np.int64)
REV3 = np.array([0, 2, 1], dtype=np.int64)
GEOM = BondSlotGeometry.from_arrays(DR3, REV3, (5, 1, 1))
_rng = np.random.default_rng(5270)
G_RAND = _rng.normal(size=(5, 12, 12)) + 1j * _rng.normal(size=(5, 12, 12))

#: real coefficients of the off-site types (complex PairHop is tested separately)
_REAL = {"CoulombInter": 0.37, "Hund": 0.53, "Ising": 0.71, "Exchange": 0.29,
         "PairLift": 0.61, "PairHop": 0.43}


#: one on-site table per type (PairHop complex)
ONSITE_TABLES = {t: onsite_decl(t, BASE[t]) for t in TYPES}
#: every off-site type on in-cell inter-orbital (cross-site) R = 0 rows, at once
CROSS_SITE_R0 = {t: {((0, 0, 0), (0, 1)): 0.7 * BASE[t],
                     ((0, 0, 0), (1, 0)): 0.7 * np.conj(BASE[t])}
                 for t in TYPES_OFFSITE}


def _bond_row(t, R, norb, v=None):
    """A mirror-asymmetric declaration of type ``t`` on the bond ``+-R``:
    ``(R, 0, norb-1) = v`` and, for ``norb >= 2``, ``(R, norb-1, 0) = 0.6 v``
    (each with its Hermitian partner)."""
    v = _REAL[t] if v is None else v
    tbl = dict(closed_row(R, 0, norb - 1, v))
    if norb >= 2:
        tbl.update(closed_row(R, norb - 1, 0, 0.6 * v))
    return {t: tbl}


def _plan(topo, nd, shape, *, consumers=("vertex",)):
    return plan_so_bond(nmat=2, nvol=int(np.prod(shape)), B=len(topo.delta_r), nd=nd,
                        consumers=consumers, caps=NOCAP)


def _plan_B(B):
    return plan_so_bond(nmat=2, nvol=5, B=B, nd=4, consumers=("vertex",), caps=NOCAP)


def _gamma_b_geom(tables, norb, shape):
    topo = so_bond_topology(compile_interaction(tables, norb), CELL)
    return production_gamma_b(tables, norb, topo, shape), slot_geometry(topo, shape)


def _gamma_b(tables, norb, shape):
    return _gamma_b_geom(tables, norb, shape)[0]


def _project_charge_sz(g, norb):
    """The longitudinal projection of ``Gamma_B`` onto the non-spin-orbital
    bond slots ``(m; a, b)`` (index ``m * norb**2 + a * norb + b``):
    ``C = 1/2 sum_{s,s'} Gamma[(m; as, bs), (m'; cs', ds')]`` and
    ``S = -1/2 sum_{s,s'} sigma_s sigma_s' Gamma[...]``."""
    nvol, ND = g.shape[0], g.shape[1]
    nd = 2 * norb
    B = ND // (nd * nd)
    t = g.reshape(nvol, B, 2, norb, 2, norb, B, 2, norb, 2, norb)
    diag = np.einsum("qmsasbntctd->qmsabntcd", t)                 # (q, m, s, a, b, m', s', c, d)
    sig = np.array([1.0, -1.0])
    C = 0.5 * diag.sum(axis=(2, 6))
    S = -0.5 * np.einsum("qmsabntcd,s,t->qmabncd", diag, sig, sig)
    n = B * norb * norb
    return S.reshape(nvol, n, n), C.reshape(nvol, n, n)


def _non_so_SC(tables, norb, shape):
    """``S``, ``C`` of the non-spin-orbital static bond path and its channels."""
    topo, _view, split = sc._resolve_bond_topology(tables, {"rvec": np.eye(3)}, norb, None)
    S, C = sc._build_bond_vertices(topo, split, norb, shape)
    return S, C, topo


def _slots(g):
    """``(G[(up up),(dn dn)], G[(dn dn),(up up)], G[(dn up),(dn up)],
    G[(up dn),(up dn)])`` of a one-band ``(4, 4)`` block (real parts)."""
    return tuple(float(g[i, j].real) for i, j in ((0, 3), (3, 0), (2, 2), (1, 1)))


def _label_neg(geom):
    return mesh_negate(tuple(geom.spatial_shape))


def _i1(gl, geom):
    """``max|Gamma_IJ(q) - Gamma_JI(-q)| / max|Gamma|`` (label form)."""
    neg = _label_neg(geom)
    return rel(gl, np.transpose(gl[neg], (0, 2, 1)))


def _i2(gl, geom):
    """``Gamma_{(-d';W,Z),(-d;Y,X)}(q) = e^{+iq.(d-d')} Gamma_{(d;X,Y),(d';Z,W)}(q)*``,
    the max-norm relative residual (label form)."""
    nvol, NL, _ = gl.shape
    B = geom.B
    nd = int(round(np.sqrt(NL // B)))
    rev = np.asarray(geom.reverse)
    dr = np.asarray(geom.delta_r, dtype=float)
    idx = np.arange(NL)
    m, x, y = idx // (nd * nd), (idx // nd) % nd, idx % nd
    f = rev[m] * nd * nd + y * nd + x                     # (d; X, Y) -> (-d; Y, X)
    ph = np.exp(1j * mesh_momenta(tuple(geom.spatial_shape)) @ dr.T)[:, m]   # e^{iq.d}
    lhs = gl[:, f[None, :], f[:, None]]                   # [(J~), (I~)] placed at [I, J]
    rhs = ph[:, :, None] * ph.conj()[:, None, :] * gl.conj()
    return rel(lhs, rhs)


def _i3(c, geom):
    """The crossed filing of the records against minus the direct filing of
    their re-paired forms ``c+_1 c_4 c+_3 c_2``.

    A crossed-vs-direct consistency check of the one filer, and as such
    structurally self-referential: both sides go through ``_file_records``,
    so a slot-mapping error common to the direct and crossed filings cancels.
    The slot mapping itself is pinned independently by the B = 1, projection
    and Hubbard tests."""
    nvol, ND = geom.nvol, geom.B * (2 * c.norb_phys) ** 2
    crossed = np.zeros((nvol, ND, ND), complex)
    _file_records(crossed, c.records, geom, parts=("crossed",))
    rep = c.records.copy()
    rep["modes"] = c.records["modes"][:, [0, 3, 2, 1]]
    rep["cells"] = c.records["cells"][:, [0, 3, 2, 1]]
    direct = np.zeros_like(crossed)
    _file_records(direct, rep, geom, parts=("direct",))
    return rel(crossed, -direct)


def _gamma_p_phase_free(g, geom):
    """The phase-free, channel-preserving exchange image
    ``Gamma_B(-q)[(m'; w,z), (m; y,x)]`` (the reading that must miss)."""
    nvol, ND, _ = g.shape
    nd = int(round(np.sqrt(ND // geom.B)))
    idx = np.arange(ND)
    m, x, y = idx // (nd * nd), (idx // nd) % nd, idx % nd
    f = m * nd * nd + y * nd + x
    return g[_label_neg(geom)][:, f[None, :], f[:, None]]


def _static_pairing_chain_observations(p):
    """The boundary observations of the static pairing chain on the plan's
    own dimensions: vertex, static bubble, dressing of one slice, parts."""
    common = dict(nvol=p.nvol, B=p.B, nd=p.nd, ND=p.ND, layout="dense")
    return [("so_bond_vertex", dict(common)),
            ("so_bond_bubble", dict(common, freq="static", nmat=p.nmat)),
            ("so_bond_dress", dict(nvol=p.nvol, ND=p.ND, batch=1)),
            ("so_bond_pairing_vertex_parts", dict(nvol=p.nvol, ND=p.ND, B=p.B))]


class TestReductions(unittest.TestCase):
    """Gate A1-U (a), (b), (c), (h): reductions of Gamma_B to the existing
    spinful vertex, to the non-spin-orbital bond vertices and to the
    one-band Hubbard vertex."""

    def test_b1_matches_spinful_vertex(self):
        for t, tbl in ONSITE_TABLES.items():
            with self.subTest(type=t):
                g = _gamma_b(tbl, 2, (4, 1, 1))
                self.assertEqual(g.shape, (4, 16, 16))
                self.assertLessEqual(rel(g, rpa_spinful_gamma(tbl, 2, (4, 1, 1))), 1e-13)

    def test_b1_self_partner_rows(self):
        for t in ("Hund", "Ising", "Exchange", "PairHop", "PairLift"):
            with self.subTest(type=t):
                tbl = {t: {((0, 0, 0), (0, 0)): 0.4, ((0, 0, 0), (1, 1)): 0.7}}
                g = _gamma_b(tbl, 2, (4, 1, 1))
                ref = rpa_spinful_gamma(tbl, 2, (4, 1, 1))
                if t in ("Hund", "PairLift"):
                    # same-orbital on-site Hund is one-body (n_s n_s = n_s) and
                    # same-orbital PairLift vanishes (c+_up c_dn c+_up c_dn = 0):
                    # direct and crossed filings cancel exactly in both builders
                    self.assertEqual(np.count_nonzero(g), 0)
                    self.assertEqual(np.count_nonzero(ref), 0)
                else:
                    self.assertLessEqual(rel(g, ref), 1e-13)

    def test_crosssite_r0_rows_in_channel_0(self):
        g = _gamma_b(CROSS_SITE_R0, 2, (4, 1, 1))
        self.assertEqual(g.shape, (4, 16, 16))                         # B = 1: channel 0 only
        self.assertLessEqual(rel(g, rpa_spinful_gamma(CROSS_SITE_R0, 2, (4, 1, 1),
                                                        orbital_site=None)), 1e-13)

    def test_longitudinal_projection_equals_non_so(self):
        self._projection(norbs=(1, 2))

    def test_longitudinal_projection_norb3(self):
        self._projection(norbs=(3,))

    def _projection(self, norbs):
        for t in ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop"):
            for norb in norbs:
                for R in (+1, -1):
                    with self.subTest(type=t, norb=norb, R=R):
                        tbl = _bond_row(t, R, norb)
                        g, geom = _gamma_b_geom(tbl, norb, (6, 1, 1))
                        S_so, C_so = _project_charge_sz(g, norb)
                        S, C, topo = _non_so_SC(tbl, norb, (6, 1, 1))
                        np.testing.assert_array_equal(geom.delta_r, topo.delta_r)
                        self.assertLessEqual(rel(S_so, S), 1e-13)
                        self.assertLessEqual(rel(C_so, C), 1e-13)

    def test_one_band_hubbard(self):
        g = _gamma_b({"CoulombIntra": {((0, 0, 0), (0, 0)): 1.3}}, 1, (3, 1, 1))
        for q in range(3):
            self.assertEqual(_slots(g[q]), (1.3, 1.3, -1.3, -1.3))
            self.assertEqual(np.count_nonzero(g[q]), 4)


class TestIdentities(unittest.TestCase):
    """Gate A1-U (d), (e), (f): the label identities I1 (particle exchange),
    I2 (Hermiticity) and I3 (crossing antisymmetry), the exchange image
    Gamma^P = Gamma_B, Hermiticity of Gamma_B(q), complex off-site PairHop and
    off-site PairLift."""

    def test_i1_i2_i3_gamma_p_hermitian(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        topo = so_bond_topology(c, CELL)
        geom = slot_geometry(topo, (5, 1, 1))
        g = so_bond_vertex(c, topo, (5, 1, 1), plan=_plan(topo, 4, (5, 1, 1)))
        self.assertEqual(g.shape, (5, 5 * 16, 5 * 16))
        gl = gamma_code_to_label(g, geom, 4, SIGN)
        self.assertLessEqual(_i1(gl, geom), 1e-13)
        self.assertLessEqual(_i2(gl, geom), 1e-13)
        self.assertLessEqual(_i3(c, geom), 1e-13)
        out = np.empty_like(g)
        self.assertIs(gamma_p(g, geom, out=out), out)
        self.assertLessEqual(rel(out, g), 1e-13)                         # Gamma^P = Gamma_B
        self.assertGreaterEqual(rel(_gamma_p_phase_free(g, geom), g), 1e-3)
        self.assertLessEqual(max(rel(g[q].conj().T, g[q]) for q in range(5)), 1e-13)

    def test_i_checks_detect_violations(self):
        """The identity residuals are not vacuous: a random code array
        violates I1 and I2, and Gamma^P differs from it."""
        gl = gamma_code_to_label(G_RAND, GEOM, 2, SIGN)
        self.assertGreaterEqual(_i1(gl, GEOM), 1e-3)
        self.assertGreaterEqual(_i2(gl, GEOM), 1e-3)
        self.assertGreaterEqual(rel(gamma_p(G_RAND, GEOM, out=np.empty_like(G_RAND)), G_RAND),
                                1e-3)

    def test_gamma_p_explicit_entries(self):
        """``Gamma^P(q)[(m;x,y),(m';z,w)] = e^{+iq.(dr_m' - dr_m)}
        Gamma_B(-q)[(rev m'; w,z), (rev m; y,x)]`` entry by entry."""
        out = gamma_p(G_RAND, GEOM, out=np.empty_like(G_RAND))
        rev, nd = [0, 2, 1], 2
        for (q, m, x, y, mp, z, w) in ((1, 1, 0, 1, 2, 1, 0), (3, 2, 1, 1, 0, 0, 1),
                                       (4, 1, 1, 0, 1, 0, 0)):
            dq = (-q) % 5
            ph = np.exp(2j * np.pi * q * (DR3[mp][0] - DR3[m][0]) / 5)
            want = ph * G_RAND[dq, rev[mp] * 4 + w * nd + z, rev[m] * 4 + y * nd + x]
            got = out[q, m * 4 + x * nd + y, mp * 4 + z * nd + w]
            self.assertAlmostEqual(got, want, places=14)

    def test_gamma_p_refusals(self):
        with self.assertRaisesRegex(ValueError, exact("gamma_p: out shares memory with gamma_b")):
            gamma_p(G_RAND, GEOM, out=G_RAND)
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: out must be a complex128 array of shape (5, 12, 12), got complex128 "
                "(5, 12, 11)")):
            gamma_p(G_RAND, GEOM, out=np.empty((5, 12, 11), complex))
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: gamma_b has shape (5, 10, 10); the geometry needs (5, B*nd**2, "
                "B*nd**2) with B = 3 and nvol = 5")):
            gamma_p(np.zeros((5, 10, 10), complex), GEOM, out=np.empty((5, 10, 10), complex))

    def test_label_round_trip(self):
        back = gamma_label_to_code(gamma_code_to_label(G_RAND, GEOM, 2, SIGN), GEOM, 2, SIGN)
        self.assertLessEqual(rel(back, G_RAND), 1e-15)

    def test_complex_pairhop_offsite(self):
        g, geom = _gamma_b_geom(_bond_row("PairHop", 1, 2, v=0.43 + 0.27j), 2, (5, 1, 1))
        self.assertGreaterEqual(np.abs(g.imag).max(), 1e-3 * np.abs(g).max())
        gl = gamma_code_to_label(g, geom, 4, SIGN)
        self.assertLessEqual(_i2(gl, geom), 1e-13)
        self.assertLessEqual(_i1(gl, geom), 1e-13)
        self.assertLessEqual(max(rel(g[q].conj().T, g[q]) for q in range(5)), 1e-13)

    def test_offsite_pairlift(self):
        g = _gamma_b(_bond_row("PairLift", 1, 2, v=0.61), 2, (5, 1, 1))
        self.assertGreater(np.abs(g).max(), 0.1)
        S, C = _project_charge_sz(g, 2)
        self.assertEqual(np.abs(S).max(), 0.0)
        self.assertEqual(np.abs(C).max(), 0.0)


class TestMeshAndPlan(unittest.TestCase):
    """Gate A1-U (g) and the admission plan: mesh-coincident channels, the
    consumer closure, the rows, and refusals before any large allocation."""

    def test_mesh_coincident_refused(self):
        with self.assertRaisesRegex(ValueError, r"both wrap to the same displacement "
                                                r"\(1, 0, 0\) modulo the mesh"):
            _gamma_b(_bond_row("CoulombInter", 1, 1), 1, (2, 1, 1))

    def test_label_outside_geometry_refused(self):
        c = compile_interaction(_bond_row("CoulombInter", 2, 1), 1)
        with self.assertRaisesRegex(ValueError, exact(
                "_file_records: label displacement (-2, 0, 0) of a CoulombInter record is not "
                "a channel of the bond geometry")):
            _file_records(np.zeros((5, 12, 12), complex), c.records, GEOM)

    def test_consumer_closure(self):
        p = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("pairing_static",), caps=NOCAP)
        self.assertEqual(p.consumers, ("vertex", "chi_static", "pairing_static"))
        for helper, obs in _static_pairing_chain_observations(p):
            p.require(helper, **obs)                                    # no ValueError
        d = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_dynamic",), caps=NOCAP)
        self.assertEqual(d.consumers, ("vertex", "chi_dynamic"))
        self.assertEqual(d.slots, ("chibar",))
        d.require("so_bond_bubble", freq="all", nvol=5, B=3, nd=2, ND=12, layout="dense",
                  nmat=8, placement="host", slots=("chibar",))
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_bubble: the admission plan has no row 'chibar' (closed consumers: "
                "vertex, chi_static, pairing_static; {})".format(p.reason))):
            p.require("so_bond_bubble", freq="all", nvol=5, B=3, nd=2, ND=12, layout="dense",
                      nmat=8, placement="host", slots=("chibar",))

    def test_non_default_convention_refused(self):
        conv = dataclasses.replace(SLOT_CONVENTION, phase_row=0)
        geom = BondSlotGeometry.from_arrays(DR3, REV3, (5, 1, 1), convention=conv)
        with self.assertRaisesRegex(ValueError, exact(
                "gamma_p: slot convention sr-:sc-:q-:nu+:rowrev1:colrev0:phr0:phc0 is not the "
                "pinned sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0")):
            gamma_p(np.zeros((5, 12, 12), complex), geom, out=np.empty((5, 12, 12), complex))
        with self.assertRaisesRegex(ValueError, exact(
                "_file_records: slot convention sr-:sc-:q-:nu+:rowrev1:colrev0:phr0:phc0 is not "
                "the pinned sr-:sc-:q-:nu+:rowrev1:colrev0:phr-:phc0")):
            _file_records(np.zeros((5, 12, 12), complex),
                          compile_interaction(_bond_row("Hund", 1, 1), 1).records, geom)

    def test_rows_unique_after_closure(self):
        for consumers in (("vertex", "chi_static", "pairing_static"), ("pairing_static",),
                          ("chi_static", "chi_dynamic", "pairing_static")):
            p = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=consumers, caps=NOCAP)
            keys = [(ph, n) for ph, n, _ in p.rows]
            self.assertEqual(len(keys), len(set(keys)))
        a = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("pairing_static",), caps=NOCAP)
        b = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("vertex", "pairing_static",
                                                               "chi_static", "vertex"), caps=NOCAP)
        self.assertEqual(a.rows, b.rows)

    def test_dress_capacity(self):
        st = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("pairing_static",), caps=NOCAP)
        dy = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_dynamic",), caps=NOCAP)
        self.assertEqual((st.dress_capacity, dy.dress_capacity), (1, 8))
        self.assertEqual(dict((n, b) for _, n, b in st.rows)["dress_batch"], 7 * 16 * 1 * 5 * 12 * 12)
        self.assertEqual(dict((n, b) for _, n, b in dy.rows)["dress_batch"], 7 * 16 * 8 * 5 * 12 * 12)
        dy.require("so_bond_dress", nvol=5, ND=12, batch=8)
        for bad in (0, -1):
            with self.assertRaisesRegex(ValueError, exact(
                    "so_bond_dress: the batch length must be a positive integer, got {}"
                    .format(bad))):
                dy.require("so_bond_dress", nvol=5, ND=12, batch=bad)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: a batch of 2 slices exceeds the plan's dress_capacity = 1 "
                "({})".format(st.reason))):
            st.require("so_bond_dress", nvol=5, ND=12, batch=2)

    def test_plan_rows_b1(self):
        p = plan_so_bond(nmat=8, nvol=4, B=1, nd=2, consumers=("vertex", "chi_static"), caps=NOCAP)
        rows = dict((n, b) for _, n, b in p.rows)
        self.assertEqual(rows["gamma_b"], 16 * 4 * 16)
        # 16 nvol (3 nmat P^2 + 3 nmat P + 8 P), P = nd^2 = 4
        self.assertEqual(rows["bubble_pair_block"], 16 * 4 * (3 * 8 * 16 + 3 * 8 * 4 + 8 * 4))
        self.assertEqual(p.nb, 8)
        self.assertEqual((p.nslices, p.ND, p.layout, p.placement), (8, 4, "dense", "host"))
        self.assertIsInstance(p.rows, tuple)
        self.assertIsInstance(p.aliases, tuple)
        self.assertIsInstance(p.slots, tuple)
        resident = 16 * 4 * 16
        chi_static = 16 * 4 * 16 * 3 + 32768 + 7 * 16 * 4 * 16
        self.assertEqual(p.peak_host, 1.25 * (resident + chi_static))
        self.assertEqual(plan_so_bond(nmat=40, nvol=4, B=1, nd=2, consumers=("vertex",),
                                      caps=NOCAP).nb, 16)
        self.assertTrue(dataclasses.is_dataclass(p) and isinstance(p, SOBondPlan))
        with self.assertRaises(dataclasses.FrozenInstanceError):
            p.nb = 3

    def test_static_pairing_peak(self):
        """B = 3, nd = 4 (ND = 48), nvol = 8, nmat = 8, pairing only."""
        p = plan_so_bond(nmat=8, nvol=8, B=3, nd=4, consumers=("pairing_static",), caps=NOCAP)
        full = 16 * 8 * 48 * 48                         # 294912: one (nvol, ND, ND) array
        pair_block = 16 * 8 * (3 * 8 * 256 + 3 * 8 * 16 + 8 * 16)    # 851968, P = nd^2 = 16
        dress = 7 * 16 * 1 * 8 * 48 * 48                # 2064384 (dress_capacity = 1)
        workspace = 2 * 16 * 48 * 48                    # 73728
        pairing_rows = {"chibar_static": full, "bubble_static_temp": full,
                        "bubble_pair_block": pair_block, "chi_static": full,
                        "dress_batch": dress, "v_bare": full, "v_fluct": full,
                        "pairing_q_workspace": workspace}
        self.assertEqual({n: b for ph, n, b in p.rows if ph == "pairing_static"}, pairing_rows)
        self.assertEqual(sum(pairing_rows.values()), 4464640)
        # the pairing phase (4464640) exceeds the static-susceptibility one (3801088);
        # the resident row is gamma_b = full
        self.assertEqual(p.peak_host, 1.25 * (full + 4464640))
        self.assertEqual(p.peak_host, 5949440.0)

    def test_static_chain_peak_within_plan(self):
        """The traced peak of the static chain (vertex, static bubble,
        dressing, pairing parts) stays within the plan's unpadded estimate
        ``peak_host / 1.25``. One physical orbital and nmat = 256 make the
        bubble's per-shift pair block (which carries the Green-function
        stacks besides the bubble block itself) the dominant row."""
        shape, nmat, norb = (8, 1, 1), 256, 1
        tables = {"CoulombIntra": {((0, 0, 0), (0, 0)): 1.0},
                  "CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}}
        c = compile_interaction(tables, norb)
        topo = so_bond_topology(c, CELL)
        geom = slot_geometry(topo, shape)
        nd, nvol = 2 * norb, int(np.prod(shape))
        ND = geom.B * nd * nd
        plan = plan_so_bond(nmat=nmat, nvol=nvol, B=geom.B, nd=nd,
                            consumers=("pairing_static",), caps=NOCAP)
        rng = np.random.default_rng(31)
        g = 0.01 * (rng.standard_normal((1, nmat, nvol, nd, nd))
                    + 1j * rng.standard_normal((1, nmat, nvol, nd, nd)))
        tracemalloc.start()
        self.addCleanup(tracemalloc.stop)
        base = tracemalloc.get_traced_memory()[0]
        gb = so_bond_vertex(c, topo, shape, plan=plan)
        sink = np.empty((nvol, ND, ND), dtype=np.complex128)
        so_bond_bubble(g, None, 10.0, topo, shape, plan=plan, freq="static", sink=sink)
        chi = so_bond_dress(sink[None], gb, plan=plan, l0=nmat // 2, spatial_shape=shape,
                            cond_tol=None)[0]
        so_bond_pairing_vertex_parts(gb, chi, geom, plan=plan)
        measured = tracemalloc.get_traced_memory()[1] - base
        self.assertLessEqual(measured, plan.peak_host / 1.25)

    def test_bubble_phase_peak_within_its_rows(self):
        """The traced peak of the bubble alone stays within its own rows:
        the static bubble within ``bubble_pair_block + bubble_static_temp``
        (the sink is the caller's), the every-slice bubble within
        ``bubble_pair_block`` (the store is preallocated). One physical
        orbital and nmat = 256 make the Green-function stacks of a shift a
        visible part of the row."""
        shape, nmat, norb = (8, 1, 1), 256, 1
        tables = {"CoulombIntra": {((0, 0, 0), (0, 0)): 1.0},
                  "CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}}
        topo = so_bond_topology(compile_interaction(tables, norb), CELL)
        geom = slot_geometry(topo, shape)
        nd, nvol = 2 * norb, int(np.prod(shape))
        ND = geom.B * nd * nd
        plan = plan_so_bond(nmat=nmat, nvol=nvol, B=geom.B, nd=nd,
                            consumers=("chi_static", "chi_dynamic"), caps=NOCAP)
        rows = {(ph, n): b for ph, n, b in plan.rows}
        g = 0.01 * (np.random.default_rng(37).standard_normal((1, nmat, nvol, nd, nd)) + 0j)
        sink = np.empty((nvol, ND, ND), dtype=np.complex128)

        def traced(call):
            tracemalloc.start()
            try:
                base = tracemalloc.get_traced_memory()[0]
                call()
                return tracemalloc.get_traced_memory()[1] - base
            finally:
                tracemalloc.stop()

        static = traced(lambda: so_bond_bubble(g, None, 10.0, topo, shape, plan=plan,
                                               freq="static", sink=sink))
        self.assertLessEqual(static, rows[("chi_static", "bubble_pair_block")]
                             + rows[("chi_static", "bubble_static_temp")])
        with flex_bond.BondBlockStore(nmat, nvol, ND, nd * nd, ("chibar",)) as store:
            every = traced(lambda: so_bond_bubble(g, None, 10.0, topo, shape, plan=plan,
                                                  freq="all", sink=store))
        self.assertLessEqual(every, rows[("chi_dynamic", "bubble_pair_block")])

    def test_plan_mismatch_and_cap_refused_before_allocation(self):
        c = compile_interaction(COMPILER_TABLES, 2)
        topo = so_bond_topology(c, CELL)
        tracemalloc.start()
        self.addCleanup(tracemalloc.stop)
        # the peak each refused call adds over what is already traced: the
        # test runner's own bookkeeping (subtest reports, its progress line)
        # is not the code under test

        def start():
            tracemalloc.reset_peak()
            return tracemalloc.get_traced_memory()[0]

        def added(base):
            return tracemalloc.get_traced_memory()[1] - base
        peaks = []
        for field, bad, plan_v, call_v in (("nd", _plan(topo, 2, (5, 1, 1)), 2, 4),
                                           ("nvol", _plan(topo, 4, (4, 1, 1)), 4, 5),
                                           ("B", _plan_B(1), 1, 5)):
            with self.subTest(field=field):
                base = start()
                with self.assertRaisesRegex(ValueError, exact(
                        "so_bond_vertex: the admission plan has {f} = {p} but the call has "
                        "{f} = {c} ({r})".format(f=field, p=plan_v, c=call_v, r=bad.reason))):
                    so_bond_vertex(c, topo, (5, 1, 1), plan=bad)
                peaks.append(added(base))
        head = ("spin-orbital bond FLEX needs an estimated peak of 0.00 GiB in phase vertex "
                "(B = 5 channels, ND = B·(2·norb_phys)² = 80, nvol = 5, Nmat = 2; one store "
                "slot = 16·Nmat·nvol·ND² bytes = 0.00 GiB) but the cap is 0.00 GiB "
                "(caps.host_bytes). Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce Nmat "
                "or the k grid, or prune bond shells.\n")
        base = start()
        with self.assertRaises(ValueError) as cm:
            plan_so_bond(nmat=2, nvol=5, B=5, nd=4, consumers=("vertex",),
                         caps=SOBondCaps(1024.0, 1e-6, None))
        peaks.append(added(base))
        msg = str(cm.exception)
        self.assertTrue(msg.startswith(head), msg)
        self.assertIn("resident  gamma_b  512000", msg)
        self.assertEqual(len(peaks), 4)
        self.assertLess(max(peaks), 100 * 1024)

    def test_plan_argument_refusals(self):
        for kw, msg in (
                (dict(consumers=("vertex", "flex")),
                 "plan_so_bond: unknown consumer 'flex'; the consumers are vertex, chi_static, "
                 "chi_dynamic, pairing_static"),
                (dict(consumers=()), "plan_so_bond: consumers must name at least one of vertex, "
                                     "chi_static, chi_dynamic, pairing_static"),
                (dict(half_grid=True), "plan_so_bond: the half bosonic grid is not available on "
                                       "the spin-orbital bond path"),
                (dict(placement_request="device"),
                 "plan_so_bond: placement_request must be 'host' on the spin-orbital bond path, "
                 "got 'device'"),
                (dict(nd=3), "plan_so_bond: nd must be a positive even integer (2·norb_phys), "
                             "got 3"),
                (dict(nb=0), "plan_so_bond: nb must be a positive integer, got 0"),
                (dict(nmat=3), "plan_so_bond: nmat must be a positive even integer, got 3"),
                (dict(caps=SOBondCaps(-1.0, -1.0, None)),
                 "plan_so_bond: caps.host_bytes must be None (unlimited) or a finite number of "
                 "bytes >= 0, got -1.0"),
                (dict(caps=SOBondCaps(float("nan"), None, None)),
                 "plan_so_bond: caps.host_bytes must be None (unlimited) or a finite number of "
                 "bytes >= 0, got nan")):
            with self.subTest(kw=kw):
                args = dict(nmat=8, nvol=5, B=3, nd=2, consumers=("vertex",), caps=NOCAP)
                args.update(kw)
                with self.assertRaisesRegex(ValueError, exact(msg)):
                    plan_so_bond(**args)

    def test_cap_admits_and_unlimited_logs(self):
        p = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("pairing_static",),
                         caps=SOBondCaps(1e9, 1.0, None))
        self.assertEqual((p.cap_host, p.cap_host_raw, p.cap_device), (1e9, 1.0, None))
        self.assertLessEqual(p.peak_host, 1e9)
        with self.assertLogs("hwave.solver.so_bond", level="INFO") as cm:
            q = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("vertex",), caps=NOCAP)
        self.assertIsNone(q.cap_host)
        self.assertIn(q.reason, "\n".join(cm.output))


def _sigma_sq(B, nd):
    """``sum_sigma n_sigma^2`` of the S_z-packed layout: ``n_0 = 2 B norb^2``,
    ``n_{+-1} = B norb^2``."""
    norb = nd // 2
    n0, n1 = 2 * B * norb * norb, B * norb * norb
    return n0 * n0 + 2 * n1 * n1


_PACKED_WHY = "so bond layout: sz_packed (S_z conserved: no spin-flip transfer, " \
              "no S_z-changing interaction)"


def _packed_plan(**kw):
    args = dict(nmat=8, nvol=8, B=3, nd=4, consumers=("pairing_static",), caps=NOCAP,
                layout="sz_packed", layout_reason=_PACKED_WHY)
    args.update(kw)
    return plan_so_bond(**args)


class TestPackedPlan(unittest.TestCase):
    """Construction check of the admission plan in the S_z-packed layout:
    the rows scale with ``sum_sigma n_sigma^2 = (6/16) ND^2``, the pair block
    and the static bubble temporary do not, the layout is observed at the
    dressing and the pairing parts, and the refusals."""

    def _scaled(self, p, S):
        nvol, nmat = p.nvol, p.nmat
        return {"gamma_b": 16 * nvol * S, "chibar_static": 16 * nvol * S,
                "chi_static": 16 * nvol * S, "v_bare": 16 * nvol * S, "v_fluct": 16 * nvol * S,
                "dress_batch": 7 * 16 * p.dress_capacity * nvol * S,
                "chibar": 16 * nmat * nvol * S, "chi": 16 * nmat * nvol * S,
                "pairing_q_workspace": 2 * 16 * S, "vertex_batch_workspace": 2 * 16 * S}

    def test_sigma_sq_is_six_sixteenths(self):
        for B, nd in ((1, 2), (3, 4), (5, 6)):
            ND = B * nd * nd
            self.assertEqual(16 * _sigma_sq(B, nd), 6 * ND * ND)

    def test_packed_rows(self):
        for B, nd in ((1, 2), (3, 4), (5, 2)):
            for consumers, keep in ((("pairing_static",), False),
                                    (("chi_static", "chi_dynamic"), True)):
                with self.subTest(B=B, nd=nd, consumers=consumers):
                    kw = dict(nmat=8, nvol=8, B=B, nd=nd, consumers=consumers, caps=NOCAP,
                              keep_chi=keep, vertex_batch="chi_dynamic" in consumers)
                    pk = plan_so_bond(layout="sz_packed", layout_reason=_PACKED_WHY, **kw)
                    de = plan_so_bond(**kw)
                    S, ND = _sigma_sq(B, nd), B * nd * nd
                    self.assertEqual((pk.layout, pk.layout_reason), ("sz_packed", _PACKED_WHY))
                    self.assertEqual([(ph, n) for ph, n, _ in pk.rows],
                                     [(ph, n) for ph, n, _ in de.rows])
                    scaled_p, scaled_d = self._scaled(pk, S), self._scaled(de, ND * ND)
                    dense_rows = {(ph, n): b for ph, n, b in de.rows}
                    seen = set()
                    for ph, n, b in pk.rows:
                        if n in scaled_p:
                            self.assertEqual(b, scaled_p[n], (ph, n))
                            self.assertEqual(dense_rows[(ph, n)], scaled_d[n], (ph, n))
                        else:
                            self.assertIn(n, ("bubble_pair_block", "bubble_static_temp"))
                            self.assertEqual(b, dense_rows[(ph, n)], (ph, n))
                        seen.add(n)
                    self.assertIn("bubble_pair_block", seen)
                    self.assertLess(pk.peak_host, de.peak_host)

    def test_packed_rows_cover_every_phase(self):
        p = _packed_plan(consumers=("pairing_static", "chi_dynamic"), keep_chi=True,
                         vertex_batch=True)
        names = {n for _, n, _ in p.rows}
        self.assertEqual(names, {"gamma_b", "chibar_static", "bubble_static_temp",
                                 "bubble_pair_block", "dress_batch", "chi_static", "chibar",
                                 "chi", "vertex_batch_workspace", "v_bare", "v_fluct",
                                 "pairing_q_workspace"})
        rows = {(ph, n): b for ph, n, b in p.rows}
        self.assertEqual(rows[("chi_static", "bubble_static_temp")], 16 * 8 * 48 * 48)
        self.assertEqual(rows[("chi_dynamic", "vertex_batch_workspace")], 2 * 16 * 864)
        self.assertEqual(rows[("resident", "gamma_b")], 16 * 8 * 864)
        # the resident vertex row is the allocation of the packed array
        self.assertEqual(rows[("resident", "gamma_b")],
                         PackedArrays.zeros(sz_sector_map(GEOM, 4), (8,)).nbytes)

    def test_reason_and_layout_reason(self):
        d = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("vertex",), caps=NOCAP)
        self.assertEqual(d.layout_reason,
                         "so bond layout: dense (reason: not selected by the caller)")
        self.assertIn("; dense layout, host placement;", d.reason)
        p = _packed_plan()
        self.assertIn("; sz_packed layout, host placement;", p.reason)
        forced = 'so bond layout: dense (reason: forced by so_bond_layout = "dense")'
        f = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("vertex",), caps=NOCAP,
                         layout="dense", layout_reason=forced)
        self.assertEqual(f.layout_reason, forced)

    def test_layout_argument_refusals(self):
        for kw, msg in (
                (dict(layout="Sz_packed", layout_reason=_PACKED_WHY),
                 "plan_so_bond: layout must be 'dense' or 'sz_packed', got 'Sz_packed'"),
                (dict(layout="sz_packed"),
                 "plan_so_bond: layout = sz_packed needs its layout_reason (the reason line "
                 "of the layout selection)"),
                (dict(layout="sz_packed", layout_reason="so bond layout: dense (reason: x)"),
                 "plan_so_bond: layout_reason must be a string starting with 'so bond layout: "
                 "sz_packed (', got 'so bond layout: dense (reason: x)'"),
                (dict(layout_reason=3),
                 "plan_so_bond: layout_reason must be a string starting with 'so bond layout: "
                 "dense (', got 3"),
                (dict(vertex_batch=True),
                 "plan_so_bond: vertex_batch needs the consumer chi_dynamic (the vertex batch "
                 "reads dressed frequency batches)")):
            with self.subTest(kw=kw):
                args = dict(nmat=8, nvol=5, B=3, nd=2, consumers=("pairing_static",),
                            caps=NOCAP)
                args.update(kw)
                with self.assertRaisesRegex(ValueError, exact(msg)):
                    plan_so_bond(**args)

    def test_layout_observed_at_dress_and_parts(self):
        pk = _packed_plan(nvol=5, consumers=("pairing_static",))
        de = plan_so_bond(nmat=8, nvol=5, B=3, nd=4, consumers=("pairing_static",), caps=NOCAP)
        # the dressing and the parts observe the layout as an optional field (dense when absent)
        de.require("so_bond_dress", nvol=5, ND=48, batch=1)
        de.require("so_bond_dress", nvol=5, ND=48, batch=1, layout="dense")
        pk.require("so_bond_dress", nvol=5, ND=48, batch=1, layout="sz_packed")
        pk.require("so_bond_pairing_vertex_parts", nvol=5, ND=48, B=3, layout="sz_packed")
        for helper, obs in (("so_bond_dress", dict(nvol=5, ND=48, batch=1)),
                            ("so_bond_pairing_vertex_parts", dict(nvol=5, ND=48, B=3))):
            with self.subTest(helper=helper):
                with self.assertRaisesRegex(ValueError, exact(
                        "{}: the admission plan has layout = sz_packed but the call has "
                        "layout = dense ({})".format(helper, pk.reason))):
                    pk.require(helper, **obs)
                with self.assertRaisesRegex(ValueError, exact(
                        "{}: the admission plan has layout = dense but the call has "
                        "layout = sz_packed ({})".format(helper, de.reason))):
                    de.require(helper, layout="sz_packed", **obs)
        with self.assertRaisesRegex(ValueError, exact(
                "SOBondPlan.require: so_bond_dress observes the fields nvol, ND, batch, "
                "layout (optional), got ND, batch, nvol, other")):
            de.require("so_bond_dress", nvol=5, ND=48, batch=1, other=1)

    def test_packed_x_with_dense_plan_refused(self):
        de = plan_so_bond(nmat=8, nvol=5, B=3, nd=4, consumers=("pairing_static",), caps=NOCAP)
        pk = _packed_plan(nvol=5, consumers=("pairing_static",))
        sm = sz_sector_map(GEOM, 4)
        xp = PackedArrays.zeros(sm, (1, 5))
        gb = np.zeros((5, 48, 48), dtype=np.complex128)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: the admission plan has layout = dense but the call has "
                "layout = sz_packed ({})".format(de.reason))):
            so_bond_dress(xp, gb, plan=de, l0=4, spatial_shape=(5, 1, 1))
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_dress: the admission plan has layout = sz_packed but the call has "
                "layout = dense ({})".format(pk.reason))):
            so_bond_dress(np.zeros((1, 5, 48, 48), complex), gb, plan=pk, l0=4,
                          spatial_shape=(5, 1, 1))
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_parts: the admission plan has layout = sz_packed but "
                "the call has layout = dense ({})".format(pk.reason))):
            so_bond_pairing_vertex_parts(gb, np.zeros_like(gb), GEOM, plan=pk)

    def test_vertex_batch_require(self):
        p = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_dynamic",), caps=NOCAP,
                         vertex_batch=True)
        self.assertEqual(dict((n, b) for _, n, b in p.rows)["vertex_batch_workspace"],
                         2 * 16 * 12 * 12)
        p.require("so_bond_pairing_vertex_batch", nvol=5, ND=12, B=3, layout="dense", batch=8)
        q = plan_so_bond(nmat=8, nvol=5, B=3, nd=2, consumers=("chi_dynamic",), caps=NOCAP)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_batch: the admission plan has no row "
                "'vertex_batch_workspace' (closed consumers: vertex, chi_dynamic; {})"
                .format(q.reason))):
            q.require("so_bond_pairing_vertex_batch", nvol=5, ND=12, B=3, layout="dense",
                      batch=8)
        with self.assertRaisesRegex(ValueError, exact(
                "so_bond_pairing_vertex_batch: a batch of 9 slices exceeds the plan's batch "
                "capacity (dress_capacity = 8) ({})".format(p.reason))):
            p.require("so_bond_pairing_vertex_batch", nvol=5, ND=12, B=3, layout="dense",
                      batch=9)

    def test_packed_cap_refusal(self):
        with self.assertRaises(ValueError) as cm:
            _packed_plan(caps=SOBondCaps(1024.0, 1e-6, None))
        msg = str(cm.exception)
        self.assertTrue(msg.startswith("spin-orbital bond Eliashberg needs an estimated peak"),
                        msg)
        self.assertIn("one store slot = 16·Nmat·nvol·Σn_σ² bytes (Σn_σ² = (6/16)·ND²) = ", msg)
        self.assertIn("resident  gamma_b  {}".format(16 * 8 * 864), msg)


def _solver_formula(dims, k, n_full_gap, n_sector_vec):
    """The sector solver's estimate, written out: per sector the larger of
    its branch (dense or ARPACK) and the null-operator row, the maximum over
    sectors, plus the candidate and kept dense gaps, the reduced vectors and
    the fixed call allowance."""
    per = []
    for d, n in zip(dims, n_sector_vec):
        if d <= max(2 * k + 1, 20):
            branch = 16 * 7 * d * d + 16 * 3 * n * d + 16 * 4 * d
        else:
            ncv = min(max(2 * k + 1, 20), d)
            work = 16 * (ncv * d + 3 * ncv * ncv + 14 * ncv + 3 * d)
            # the default (unshifted) solve keeps the spectral-radius set and the
            # first eigenvalues; no second working set at these dimensions
            branch = (work + 24 * n + 6 * 16 * n + 16 * 3 * k * (d + n)
                      + 16 * (min(6, d - 2) * (d + 1) + k))
        null = 16 * (min(k, d) * (d + n) + n_full_gap)
        per.append(max(branch, null))
    return (max(per) + 16 * (k + sum(min(k, d) for d in dims)) * n_full_gap
            + 16 * sum(min(k, d) * d for d in dims) + 16 * 1024)


class TestSectorEstimators(unittest.TestCase):
    """Construction check of the two memory estimators of the S_z-sector
    eigen solve: their values equal the written-out formulas on both
    branches and both null-operator regimes, and the argument refusals."""

    def test_sector_solver_formula(self):
        cases = []
        for k in (1, 10):
            for dims in ((7,), (20,), (21,), (60,), (7, 20, 21, 60), (0, 21, 3)):
                for extra in (0, 5, 4000):
                    nvec = tuple(2 * d + extra for d in dims)
                    # a long gap: the null-operator row (one scattered gap) wins
                    for pad in (3, 10 ** 6):
                        cases.append((dims, k, sum(nvec) + pad, nvec))
        for dims, k, nfull, nvec in cases:
            with self.subTest(dims=dims, k=k, nvec=nvec):
                got = estimate_sector_solver_bytes(list(dims), k, nfull, n_sector_vec=nvec)
                self.assertIsInstance(got, int)
                self.assertEqual(got, _solver_formula(dims, k, nfull, nvec))

    def test_sector_solver_branches(self):
        # k = 1: dim 20 is dense, 21 ARPACK (ncv = 20); k = 10: 21 dense, 60 ARPACK (ncv = 21)
        tail = 16 * 1024
        self.assertEqual(estimate_sector_solver_bytes([20], 1, 40, n_sector_vec=[40]),
                         max(16 * 7 * 400 + 16 * 3 * 40 * 20 + 16 * 80, 16 * (60 + 40))
                         + 16 * 2 * 40 + 16 * 20 + tail)
        # the retained sets: the spectral_shift="auto" spectral-radius solve (6 pairs) and the
        # eigenvalues of the first which='LM' solve
        arp = (16 * (20 * 21 + 3 * 400 + 14 * 20 + 3 * 21) + 24 * 42 + 6 * 16 * 42
               + 16 * 3 * (21 + 42) + 16 * (6 * 22 + 1))
        self.assertEqual(estimate_sector_solver_bytes([21], 1, 42, n_sector_vec=[42]),
                         max(arp, 16 * (63 + 42)) + 16 * 2 * 42 + 16 * 21 + tail)
        self.assertEqual(estimate_sector_solver_bytes([21], 10, 42, n_sector_vec=[42]),
                         16 * 7 * 441 + 16 * 3 * 42 * 21 + 16 * 84 + 16 * 20 * 42
                         + 16 * 10 * 21 + tail)
        arp = (16 * (21 * 60 + 3 * 441 + 14 * 21 + 3 * 60) + 24 * 120 + 6 * 16 * 120
               + 16 * 30 * (60 + 120) + 16 * (6 * 61 + 10))
        self.assertEqual(estimate_sector_solver_bytes([60], 10, 120, n_sector_vec=[120]),
                         arp + 16 * 20 * 120 + 16 * 10 * 60 + tail)

    def test_sector_solver_refusals(self):
        name = "estimate_sector_solver_bytes"
        for args, kw, msg in (
                (([], 1, 4), dict(n_sector_vec=[]),
                 "dims must be a non-empty sequence of integers >= 0, got []"),
                (([2.0], 1, 4), dict(n_sector_vec=[4]),
                 "dims must be a non-empty sequence of integers >= 0, got [2.0]"),
                (([-1], 1, 4), dict(n_sector_vec=[4]),
                 "dims must be a non-empty sequence of integers >= 0, got [-1]"),
                (([2, 1], 1, 8), dict(n_sector_vec=[4]),
                 "n_sector_vec must be a sequence of 2 integers >= 0 (one per dims entry), "
                 "got [4]"),
                (([2], 1, 8), dict(n_sector_vec=[True]),
                 "n_sector_vec must be a sequence of 1 integers >= 0 (one per dims entry), "
                 "got [True]"),
                (([3], 1, 8), dict(n_sector_vec=[5]),
                 "dims[0] = 3 exceeds n_sector_vec[0] / 2 = 5 / 2"),
                (([2], 0, 8), dict(n_sector_vec=[4]), "k must be an integer >= 1, got 0"),
                (([2], True, 8), dict(n_sector_vec=[4]), "k must be an integer >= 1, got True"),
                (([2], 1, 3), dict(n_sector_vec=[4]),
                 "n_full_gap must be an integer >= sum(n_sector_vec) = 4, got 3"),
                (([2], 1, 4.0), dict(n_sector_vec=[4]),
                 "n_full_gap must be an integer >= sum(n_sector_vec) = 4, got 4.0")):
            with self.subTest(msg=msg):
                with self.assertRaisesRegex(ValueError, exact("{}: {}".format(name, msg))):
                    estimate_sector_solver_bytes(*args, **kw)

    def test_kernel_workspace_formula(self):
        for B in (1, 3, 5):
            for counts in ((1, 2, 1), (4, 8, 4), (9, 18, 9)):
                for nvol, nfreq in ((1, 1), (6, 4)):
                    with self.subTest(B=B, counts=counts, nvol=nvol, nfreq=nfreq):
                        sq = sum(n * n for n in counts)
                        self.assertEqual(
                            estimate_packed_kernel_workspace_bytes(B, counts, nvol, nfreq),
                            {"hoisted_blocks": 16 * B * B * sq * nvol * nfreq,
                             "g2_sectors": 16 * sq * nvol * nfreq,
                             "hoist_workspace": 16 * 6 * max(counts) ** 2 * nvol * nfreq,
                             "g2_gather": 16 * max(counts) ** 2 * nvol * nfreq,
                             "matvec_workspace": 16 * (max(3 + B, 6) * max(counts) * nvol
                                                       * nfreq + 3 * 8192),
                             "call_overhead": 16 * 1024})

    def test_kernel_workspace_refusals(self):
        name = "estimate_packed_kernel_workspace_bytes"
        for args, msg in (
                ((0, (1, 2, 1), 1, 1), "B must be an integer >= 1, got 0"),
                ((1.0, (1, 2, 1), 1, 1), "B must be an integer >= 1, got 1.0"),
                ((1, (1, 2), 1, 1), "pp_pair_counts must be three positive integers, got (1, 2)"),
                ((1, (1, 0, 1), 1, 1),
                 "pp_pair_counts must be three positive integers, got (1, 0, 1)"),
                ((1, (1, 2, 1), 0, 1), "nvol must be an integer >= 1, got 0"),
                ((1, (1, 2, 1), 1, False), "nfreq must be an integer >= 1, got False")):
            with self.subTest(msg=msg):
                with self.assertRaisesRegex(ValueError, exact("{}: {}".format(name, msg))):
                    estimate_packed_kernel_workspace_bytes(*args)


if __name__ == "__main__":
    unittest.main()
