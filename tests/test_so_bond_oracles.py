"""The table-driven operator builder and the exact-diagonalization oracle of
the spin-orbital bond-resolved vertex.

``tests/so_bond_oracles.py`` is the independent reference of every later
spin-orbital bond gate, so it is pinned here against existing, already
validated pieces and never against the production compiler:

* it imports neither ``hwave.solver.so_records`` nor ``hwave.solver.so_bond``;
* its operator builder reproduces the term builders of the earlier
  exact-diagonalization tests (in their frame, which is the global mirror of
  the file frame), and misses them in the unmirrored frame;
* its blocked sparse Jordan-Wigner diagonalization reproduces a dense one,
  and its free correlator is the analytic Lindhard function in the label form
  (spin-block mode index; the interleaved reading misses);
* the label filing rule of a quartic monomial, used with the analytic
  Lindhard function as ``-chibar(q) Gamma(-q) chibar(q)``, reproduces the
  first-order vertex part of the exact bond correlators, and four wrong
  readings of the rule miss.
"""

import ast
import os
import unittest

import numpy as np

from tests import ed_oracle_util
from tests.ed_oracle_util import EDFixture
from tests.so_bond_oracles import (
    BASE, COMPILER_TABLES, FIXTURE_I, FIXTURE_II, FIXTURE_X, FIXTURES, JWED, NEEDS_RICHARDSON,
    SOC_II, TRANSFER_II, TRUNCATION_GUARD, TYPES, TYPES_OFFSITE, FixtureSpec, Monomial, apply_terms,
    build_h1, build_terms, closed_row, ed_cached, ed_mode, expand_tables, file_gamma_label,
    first_order_prediction, fock_states, free_density, free_static_corr, freeze_table,
    gamma_code_to_label, gamma_label_to_code, geom_of, label_list, lindhard_static, rel,
    sketch_set, sketch_value, truncation, vertex_part_sketches)
from tests.test_bond_longitudinal_ed import _terms_pairhop_offsite
from tests.test_bond_transverse_ed import (_terms_exchange_offsite, _terms_hund_offsite,
                                           _terms_ising_offsite, _terms_pairlift_offsite)
import tests.test_so_pairing_ed as so_pairing_ed
import tests.test_so_pairing_vertex as so_pairing_vertex

_ORACLE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "so_bond_oracles.py")
_EXISTING = {"Hund": _terms_hund_offsite, "Ising": _terms_ising_offsite,
             "Exchange": _terms_exchange_offsite, "PairLift": _terms_pairlift_offsite,
             "PairHop": _terms_pairhop_offsite}


def _imported_modules(tree):
    """Every module an ``import`` statement names, with ``from a import b``
    also contributing ``a.b`` (so ``from hwave.solver import so_bond`` is seen)."""
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            mod = ("." * node.level) + (node.module or "")
            names.add(mod)
            names.update("{}.{}".format(mod, alias.name) for alias in node.names)
    return names


# ---------------------------------------------------------------------------
# Frame pin helpers
# ---------------------------------------------------------------------------

def _mirror_asymmetric_row(t):
    """``X_01(+1) = v``, ``X_01(-1) = 0.6 v`` with their partners: the row
    differs from its mirror image and from its orbital transpose."""
    v = BASE[t]
    return {**closed_row(1, 0, 1, v), **closed_row(-1, 0, 1, 0.6 * v)}


def _existing_builder(t, tbl, *, mirror, L=3, norb=2):
    """The existing exact-diagonalization term builders on every bond of
    ``tbl`` (one call per partner pair, full coefficient). They put orbital
    ``a`` of row ``(R, a, b)`` on site ``j + R`` and ``b`` on site ``j``;
    ``mirror=True`` calls them at ``-R``."""
    fx = EDFixture(L, norb, {}, [0.0] * norb, 1.0, 0.0)
    done, terms = set(), []
    for (R, (a, b)) in sorted(tbl):
        partner = (tuple(-x for x in R), (b, a))
        if (R, (a, b)) in done:
            continue
        done.update({(R, (a, b)), partner})
        p = np.conj if t == "PairHop" else (lambda z: z)
        v = 0.5 * (tbl[(R, (a, b))] + p(tbl.get(partner, 0.0)))
        r = -R[0] if mirror else R[0]
        if t == "CoulombInter":
            terms += ed_oracle_util.canonical_density_terms(fx, [(a, b, r, v)])
        else:
            terms += _EXISTING[t](fx, a, b, r, v)
    return terms


def _action(terms, nmode):
    """``H |psi>`` on four fixed random vectors in each of the N = 1, 2, 3
    sectors, concatenated."""
    rng = np.random.default_rng(31)
    out = []
    for n in (1, 2, 3):
        states = fock_states(nmode, n)
        vecs = rng.normal(size=(len(states), 4)) + 1j * rng.normal(size=(len(states), 4))
        out.append(apply_terms(terms, states, vecs).ravel())
    return np.concatenate(out)


# ---------------------------------------------------------------------------
# Jordan-Wigner helpers
# ---------------------------------------------------------------------------

def _dense_static_corr(h1, terms, mu, beta, A, B):
    """Dense reference: the full ``2^nmode`` Fock space, Kronecker-product
    Jordan-Wigner operators, one ``eigh``; returns the Kubo correlator and
    the one-body density."""
    nmode = h1.shape[0]
    z, eye = np.diag([1.0, -1.0]), np.eye(2)
    a = np.array([[0.0, 1.0], [0.0, 0.0]])
    C = []
    for m in range(nmode):
        op = np.array([[1.0]])
        # mode m is bit m of the state; the kron puts mode nmode-1 leftmost
        for k in reversed(range(nmode)):
            op = np.kron(op, z if k < m else (a if k == m else eye))
        C.append(op.astype(complex))
    CD = [c.conj().T for c in C]
    P = np.array([[CD[p] @ C[q] for q in range(nmode)] for p in range(nmode)])
    num = sum(P[m, m] for m in range(nmode))
    H = np.einsum("pq,pqab->ab", h1, P)
    for (p, q, r, s, v) in terms:
        H = H + v * P[p, q] @ P[r, s]
    E, V = np.linalg.eigh(H - mu * num)
    E = E - E.min()
    w = np.exp(-beta * E)
    w /= w.sum()
    dE = E[None, :] - E[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(dE) > 1e-10, (w[:, None] - w[None, :]) / dE, beta * w[:, None])

    def op(M):
        return V.conj().T @ np.einsum("pq,pqab->ab", M, P) @ V
    Ae, Be = op(A), op(B)
    corr = np.sum(Ae * Be.T * K) - beta * np.sum(w * np.diag(Ae)) * np.sum(w * np.diag(Be))
    n = np.einsum("pqab,ba->pq", P, (V * w) @ V.conj().T)
    return corr, n


def _blocked_vs_dense():
    """L = 4, norb_phys = 1 (8 modes): the blocked sparse ED against the
    dense one, on a generic interacting Hamiltonian and random one-body
    operators."""
    rng = np.random.default_rng(41)
    tbl = {"CoulombIntra": {((0, 0, 0), (0, 0)): 0.9},
           "CoulombInter": closed_row(1, 0, 0, 0.4),
           "Exchange": closed_row(1, 0, 0, 0.3),
           "PairHop": closed_row(1, 0, 0, 0.2 + 0.35j),
           "PairLift": closed_row(1, 0, 0, 0.25)}
    transfer = {}
    for (R, A, B, t) in ((0, 0, 1, 0.2 - 0.1j), (1, 0, 0, -1.0), (1, 1, 1, -0.8),
                         (1, 0, 1, 0.3 + 0.4j), (0, 0, 0, 0.1), (0, 1, 1, -0.2)):
        transfer[((R, 0, 0), (A, B))] = transfer.get(((R, 0, 0), (A, B)), 0) + t
        key = ((-R, 0, 0), (B, A))
        if key != ((R, 0, 0), (A, B)):
            transfer[key] = transfer.get(key, 0) + np.conj(t)
    h1 = build_h1(transfer, 4, 1)
    terms = build_terms(tbl, 4, 1)
    A = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    B = rng.normal(size=(8, 8)) + 1j * rng.normal(size=(8, 8))
    ed = JWED(h1, terms, mu=0.2, beta=1.5)
    want_corr, want_n = _dense_static_corr(h1, terms, 0.2, 1.5, A, B)
    return max(abs(ed.static_corr(A, B) - want_corr) / abs(want_corr), rel(ed.density(), want_n))


#: L = 2, norb_phys = 2 (8 modes), the generic two-orbital band of fixture (ii):
#: the spin-block and the interleaved in-cell index differ here
_FX_L2_NP2 = FixtureSpec.make(2, 2, 2.0, 0.8, dict(FIXTURE_II.transfer_dict), {}, zeeman=0.1)


def _interleaved_labels(labels, norb):
    """The labels read with ``X`` as the interleaved in-cell index ``2 a + s``
    (re-expressed in the spin-block index the ED helpers take)."""
    def sb(X):
        a, s = divmod(X, 2)
        return s * norb + a
    return [(d, sb(X), sb(Y)) for (d, X, Y) in labels]


def _free_dev(fixture, *, qs=(1, 2), interleaved_labels=False):
    """Free correlator from the exact diagonalization at zero coupling
    against ``u^T chibar(q) v`` with the analytic label-form Lindhard
    function, two random sketches per q."""
    geom = fixture.geom
    labels = label_list(geom, fixture.nd)
    if interleaved_labels:
        labels = _interleaved_labels(labels, fixture.norb_phys)
    chibar = lindhard_static(*fixture.band, geom)
    ed = ed_cached(fixture, 0.0)
    rng = np.random.default_rng(43)
    got, want = [], []
    for q in qs:
        for _ in range(2):
            u = rng.normal(size=len(labels)) + 1j * rng.normal(size=len(labels))
            v = rng.normal(size=len(labels)) + 1j * rng.normal(size=len(labels))
            got.append(sketch_value(ed, u, v, q, labels, fixture.L, fixture.norb_phys))
            want.append(u @ chibar[q] @ v)
    return rel(np.array(got), np.array(want))


def _filing_check(fixture, *, qs, sketches):
    """``(dev, variants, targeted)``: the first-order vertex part of the
    exact bond correlators against ``-chibar(q) Gamma(-q) chibar(q)`` with
    ``Gamma`` filed from the table-driven monomials; the same for four wrong
    readings of the filing; and each targeted sketch's ``|ref|`` relative to
    the largest random one."""
    name = next(k for k, f in FIXTURES.items() if f == fixture)
    sk = sketch_set(fixture, sketches)
    ref = vertex_part_sketches(fixture, sketches=sketches, qs=qs,
                               richardson=NEEDS_RICHARDSON[name])
    geom, nd = fixture.geom, fixture.nd
    chibar = lindhard_static(*fixture.band, geom)
    mono = expand_tables(fixture.tables_dict, fixture.norb_phys)

    def dev(gamma):
        pred = first_order_prediction(chibar, gamma, geom.spatial_shape, sk, qs)
        return float(np.abs(pred - ref).max() / np.abs(ref).max())
    good = file_gamma_label(mono, geom, nd)
    variants = {"q_flipped": dev(file_gamma_label(mono, geom, nd, q_sign=-1)),
                "direct_only": dev(file_gamma_label(mono, geom, nd, crossed=False)),
                "crossed_plus": dev(file_gamma_label(mono, geom, nd, crossed_sign=+1.0)),
                "half": dev(0.5 * good)}
    n_random = 8
    scale = np.abs(ref[:, :n_random]).max()
    targeted = [float(np.abs(ref[:, i]).max() / scale) for i in range(n_random, len(sk))]
    return dev(good), variants, targeted


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestOracleIndependence(unittest.TestCase):
    def test_imports(self):
        with open(_ORACLE_PATH, encoding="utf-8") as fh:
            names = _imported_modules(ast.parse(fh.read()))
        self.assertFalse({"hwave.solver.so_records", "hwave.solver.so_bond"} & names)
        self.assertFalse({n for n in names if n.endswith("so_records") or n.endswith("so_bond")})

    def test_import_scan_sees_every_form(self):
        tree = ast.parse("import hwave.solver.so_bond\n"
                         "from hwave.solver import so_records\n"
                         "from hwave.solver.so_bond import SLOT_CONVENTION\n")
        names = _imported_modules(tree)
        self.assertLessEqual({"hwave.solver.so_bond", "hwave.solver.so_records"}, names)


class TestExpandTables(unittest.TestCase):
    def test_modes_are_spin_block(self):
        # norb_phys = 2, Exchange row (R, 0, 1): a up = 0, b up = 1, b dn = 3, a dn = 2
        m = expand_tables({"Exchange": closed_row(1, 0, 1, 0.29)}, 2)
        self.assertIn(Monomial(0.29, (0, 1, 3, 2), ((0, 0, 0), (1, 0, 0), (1, 0, 0), (0, 0, 0))), m)
        self.assertEqual(ed_mode(1, 3, 2), 4 * 1 + 2 * 1 + 1)

    def test_type_keys_case_insensitive(self):
        a = expand_tables({"exchange": closed_row(1, 0, 1, 0.29)}, 2)
        b = expand_tables({"Exchange": closed_row(1, 0, 1, 0.29)}, 2)
        self.assertEqual(a, b)

    def test_unknown_type_refused(self):
        with self.assertRaises(ValueError) as cm:
            expand_tables({"Hubbard": closed_row(1, 0, 1, 0.29)}, 2)
        self.assertEqual(str(cm.exception), "expand_tables: unknown interaction type 'Hubbard'")


class TestOracleBuilderFrame(unittest.TestCase):
    def test_oracle_builder_frame(self):
        # norb_phys = 2, 3 cells, 12 modes
        for t in TYPES_OFFSITE:
            tbl = _mirror_asymmetric_row(t)
            mine = _action(build_terms({t: tbl}, 3, 2), 12)
            mirrored = _action(_existing_builder(t, tbl, mirror=True), 12)
            unmirrored = _action(_existing_builder(t, tbl, mirror=False), 12)
            self.assertLessEqual(rel(mine, mirrored), 1e-13, t)
            self.assertGreaterEqual(rel(mine, unmirrored), 1e-3, t)

    def test_onsite_rows(self):
        L = so_pairing_vertex.NX
        fx = EDFixture(L, 2, {}, [0.0, 0.0], 1.0, 0.0)
        for t in TYPES:
            mine = _action(build_terms(so_pairing_vertex._onsite_decl(t, BASE[t]), L, 2), 4 * L)
            want = _action(so_pairing_vertex._onsite_terms(t, fx, BASE[t]), 4 * L)
            self.assertLessEqual(rel(mine, want), 1e-13, t)


class TestFixtures(unittest.TestCase):
    def test_fixture_ii_band_is_the_pairing_chain(self):
        self.assertEqual(TRANSFER_II, so_pairing_ed.TRANSFER)
        self.assertEqual(SOC_II, so_pairing_ed.SOC)

    def test_fixtures_are_hashable_content_keys(self):
        self.assertEqual(hash(FIXTURE_I), hash(FixtureSpec(*FIXTURE_I)))
        self.assertNotEqual(FIXTURE_I, FIXTURE_I._replace(zeeman=0.15))
        self.assertEqual(set(FIXTURE_X), set(FIXTURES))
        self.assertEqual(set(NEEDS_RICHARDSON), set(FIXTURES))

    def test_shared_tables_hermitian_closed(self):
        from hwave.solver import declarations
        for t, tbl in COMPILER_TABLES.items():
            declarations.validate_hermitian_closure(t, tbl)
        for f in FIXTURES.values():
            for t, tbl in f.tables_dict.items():
                declarations.validate_hermitian_closure(t, tbl)

    def test_freeze_refuses_duplicate_keys(self):
        with self.assertRaises(ValueError) as cm:
            freeze_table([(((1, 0, 0), (0, 1)), 0.1), ((np.array([1, 0, 0]), (0, 1)), 0.2)])
        self.assertEqual(str(cm.exception),
                         "freeze_table: entry R=(1, 0, 0) orbitals (0, 1) given twice")


class TestJWED(unittest.TestCase):
    def test_blocked_matches_dense(self):
        self.assertLessEqual(_blocked_vs_dense(), 1e-12)

    def test_free_one_body_matches_ed(self):
        rng = np.random.default_rng(47)
        h1 = FIXTURE_I.h1
        A = rng.normal(size=h1.shape) + 1j * rng.normal(size=h1.shape)
        B = rng.normal(size=h1.shape) + 1j * rng.normal(size=h1.shape)
        ed = JWED(h1, [], mu=FIXTURE_I.mu, beta=FIXTURE_I.beta)
        want = ed.static_corr(A, B)
        self.assertLessEqual(abs(free_static_corr(h1, FIXTURE_I.mu, FIXTURE_I.beta, A, B) - want)
                             / abs(want), 1e-12)
        self.assertLessEqual(rel(free_density(h1, FIXTURE_I.mu, FIXTURE_I.beta), ed.density()),
                             1e-12)

    def test_free_correlator_is_the_lindhard(self):
        self.assertLessEqual(_free_dev(FIXTURE_I), 1e-10)

    def test_free_correlator_norb2_index_pin(self):
        # the two-cell mesh has q in {0, 1} only
        self.assertLessEqual(_free_dev(_FX_L2_NP2, qs=(0, 1)), 1e-10)
        self.assertGreaterEqual(_free_dev(_FX_L2_NP2, qs=(0, 1), interleaved_labels=True), 1e-3)


class TestFilingRuleAgainstED(unittest.TestCase):
    """The label filing rule (direct ``v e^{+iq.t}`` and crossed
    ``-v e^{+iq.t}`` with both orders) against the first-order vertex part
    of the exact bond correlators; test-side filer only."""

    def test_label_rule_matches_ed(self):
        dev, variants, targeted = _filing_check(FIXTURE_I, qs=(1, 2), sketches="fast")
        self.assertLessEqual(dev, 1e-6)
        for name in ("q_flipped", "direct_only", "crossed_plus", "half"):
            self.assertGreaterEqual(variants[name], 1e-3, name)
        self.assertEqual(len(targeted), 3)
        self.assertTrue(all(t >= 1e-2 for t in targeted), targeted)

    def test_step_truncation_margin(self):
        tr = truncation(FIXTURE_I)
        self.assertEqual(NEEDS_RICHARDSON["i"], tr > TRUNCATION_GUARD)
        self.assertLessEqual(tr, TRUNCATION_GUARD)

    def test_unpinned_fixture_needs_an_explicit_step(self):
        sub = FIXTURE_I._replace(zeeman=0.05)
        with self.assertRaises(ValueError) as cm:
            vertex_part_sketches(sub, sketches="fast", qs=(1,))
        self.assertEqual(str(cm.exception),
                         "vertex_part_sketches: no pinned step for this fixture; pass x=")


class TestSketchSet(unittest.TestCase):
    def test_sector_pairs(self):
        # PairLift in the run: (sigma, sigma); without it: (sigma, -sigma)
        geom = FIXTURE_I.geom
        # norb_phys = 1: the spin-block index is the spin (0 up, 1 down); the
        # S_z change of rho_XY is s_Y - s_X
        sector = np.array([Y - X for (_, X, Y) in label_list(geom, 2)])
        no_lift = FIXTURE_I._replace(tables=tuple(
            (t, tb) for t, tb in FIXTURE_I.tables if t != "PairLift"))
        for fx, rule in ((FIXTURE_I, +1), (no_lift, -1)):
            sk = sketch_set(fx, "fast")
            self.assertEqual(len(sk), 11)
            for (u, v), s in zip(sk[8:], (0, +1, -1)):
                self.assertEqual(set(sector[np.abs(u) > 0]), {s})
                self.assertEqual(set(sector[np.abs(v) > 0]), {rule * s})
        self.assertEqual(len(sketch_set(FIXTURE_I, "full")), 14)


class TestGammaSlotMap(unittest.TestCase):
    def test_round_trip_and_explicit_entries(self):
        geom, nd = FIXTURE_I.geom, 2
        rng = np.random.default_rng(53)
        shape = (5, geom.B * nd * nd, geom.B * nd * nd)
        g = rng.normal(size=shape) + 1j * rng.normal(size=shape)
        for sign in ((-1, -1), (1, -1), (-1, 1), (1, 1)):
            code = gamma_label_to_code(g, geom, nd, sign)
            self.assertLessEqual(rel(gamma_code_to_label(code, geom, nd, sign), g), 1e-15)
        # (sr, sc) = (-1, -1): Gamma_B(q)[(m;x,y),(m';z,w)]
        #   = e^{+iq.dr_m'} Gamma(q)[(-dr_m; x,y), (+dr_m'; w,z)]
        code = gamma_label_to_code(g, geom, nd, (-1, -1))
        lab = {k: i for i, k in enumerate(label_list(geom, nd))}
        dr = [tuple(d) for d in geom.delta_r]
        for (q, m, x, y, mp, z, w) in ((2, 1, 0, 1, 4, 1, 0), (3, 3, 1, 1, 2, 0, 1)):
            want = np.exp(2j * np.pi * q * dr[mp][0] / 5) * g[
                q, lab[(tuple(-c for c in dr[m]), x, y)], lab[(dr[mp], w, z)]]
            got = code[q, m * nd * nd + x * nd + y, mp * nd * nd + z * nd + w]
            self.assertAlmostEqual(got, want, places=14)

    def test_wrong_shape_refused(self):
        with self.assertRaises(ValueError) as cm:
            gamma_label_to_code(np.zeros((5, 4, 4)), FIXTURE_I.geom, 2, (-1, -1))
        self.assertEqual(str(cm.exception),
                         "gamma_label_to_code: expected shape (5, 20, 20), got (5, 4, 4)")


class TestGeometry(unittest.TestCase):
    def test_labels_of_fixture_i(self):
        g = geom_of(FIXTURE_I)
        self.assertEqual([tuple(d) for d in g.delta_r],
                         [(0, 0, 0), (-1, 0, 0), (1, 0, 0), (-2, 0, 0), (2, 0, 0)])
        self.assertEqual(list(g.reverse), [0, 2, 1, 4, 3])

    def test_missing_displacement_refused(self):
        geom = geom_of(FIXTURE_II)
        mono = expand_tables({"CoulombInter": closed_row(2, 0, 1, 0.3)}, 2)
        with self.assertRaises(ValueError) as cm:
            file_gamma_label(mono, geom, 4)
        self.assertEqual(str(cm.exception),
                         "file_gamma_label: displacement (-2, 0, 0) is not a label of the geometry")


if __name__ == "__main__":
    unittest.main()
