"""Off-site Exchange in the bond-resolved FLEX (GitHub issue #181).

An off-site Exchange ``J c+_{a up} c_{b up} c+_{b dn} c_{a dn}`` (= ``-J S+_a S-_b``)
has no Hartree content; its first order is the bond (Fock) crossing, and its
second order mixes spin-flip and density vertices across bonds. With an
off-site Exchange declared, the second order of EVERY block of the bond map
is built from :func:`hwave.solver.second_order.build_bond_w2_factors`: every
record pair of the antisymmetrised interaction at weight 1/2 in every
placement of the bond-enlarged basis.

What is pinned here:

* the factorised kernel ``1/2 sum A chibar B`` against a test-side
  enumeration of the same kernel in the bond-enlarged basis (written from
  the independent real-space oracle's records,
  ``tests/test_second_order_oracle.py``), element by element;
* the second-order self-energy of ONE production bond map (frozen bare
  ``G``) against that oracle's EXACT second order, through the production
  transport, for Exchange alone and mixed with on-site and off-site partners
  on the same and on other bonds -- including the Heisenberg bond -- and,
  as the falsification, the previous mixed-block rule missing it;
* the SU(2) criterion ``J_Exchange = J_Hund - 2 J_Ising`` of a bond (exact
  diagonalization) and the warning it drives;
* the policy: without the bond channels an off-site Exchange is refused.
"""
import itertools
import logging
import os
import shutil
import tempfile
import types
import unittest

import numpy as np

import tests.test_flex_bond_sopt as sopt
from tests.test_second_order_factors import _split_for
from tests.test_second_order_oracle import (_BETA, _SHAPE, _bare_green, _coords, _g,
                                            _gamma_entries, _site, _site_of, oracle_records,
                                            oracle_sigma2)

NORB = 2
UP, DN = 0, 1
EXCH = [(1, 0, 0, 1, 2, 1.0, 0.0), (-1, 0, 0, 2, 1, 1.0, 0.0)]
_KANAMORI = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.8, 0.0), (0, 0, 0, 2, 2, 0.8, 0.0)],
             "CoulombInter": [(0, 0, 0, 1, 2, 0.5, 0.0), (0, 0, 0, 2, 1, 0.5, 0.0),
                              (0, 1, 0, 1, 1, 0.7, 0.0), (0, -1, 0, 1, 1, 0.7, 0.0)],
             "Hund": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)],
             "Exchange": EXCH + [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)]}
CASES = {
    "X": {"Exchange": EXCH},
    "I+X (Heisenberg)": {"Ising": [r[:5] + (-0.5, 0.0) for r in EXCH], "Exchange": EXCH},
    "V+X same bond": {"CoulombInter": [r[:5] + (0.7, 0.0) for r in EXCH], "Exchange": EXCH},
    "U+X": {"CoulombIntra": [(0, 0, 0, 1, 1, 0.8, 0.0), (0, 0, 0, 2, 2, 0.8, 0.0)],
            "Exchange": EXCH},
    "V other bond+X": {"CoulombInter": [(0, 1, 0, 1, 1, 0.7, 0.0), (0, -1, 0, 1, 1, 0.7, 0.0)],
                       "Exchange": EXCH},
    "Kanamori+V+X": _KANAMORI,
    # off-site PairHop (issue #181 Phase C): alone, with density terms and
    # with Exchange on the same bond, and with an on-site U
    "PH": {"PairHop": [r[:5] + (0.6, 0.0) for r in EXCH]},
    "V+PH same bond": {"CoulombInter": [r[:5] + (0.7, 0.0) for r in EXCH],
                       "PairHop": [r[:5] + (0.6, 0.0) for r in EXCH]},
    "X+PH same bond": {"Exchange": EXCH, "PairHop": [r[:5] + (0.6, 0.0) for r in EXCH]},
    "U+PH other bond": {"CoulombIntra": [(0, 0, 0, 1, 1, 0.8, 0.0), (0, 0, 0, 2, 2, 0.8, 0.0)],
                        "PairHop": [(0, 1, 0, 1, 1, 0.6, 0.0), (0, -1, 0, 1, 1, 0.6, 0.0)]},
    # off-site PairLift on another bond (inert at first order on the
    # paramagnet, but present at second order)
    "PL other bond+X": {"PairLift": [(0, 1, 0, 1, 2, 0.4, 0.0), (0, -1, 0, 2, 1, 0.4, 0.0)],
                        "Exchange": EXCH},
}


def _delta_r(rows):
    shells = {(0, 0, 0)}
    for rr in rows.values():
        for r in rr:
            shells.add(tuple(r[:3]))
            shells.add(tuple(-x for x in r[:3]))
    return sorted(shells, key=lambda t: (t != (0, 0, 0), t))


def _enumerated_w2(G_kw, beta, recs, delta_r, norb=NORB, shape=_SHAPE):
    """The exact second-order kernel in the bond-enlarged basis, enumerated in
    real space from the oracle's records: every record pair at weight 1/2;
    the external pair ``(alpha; c, a)`` is ``p = (0, a, up)``,
    ``s = (-R_alpha, c, sigma_s)``, ``(beta; d, b)`` is ``p' = (j, b, up)``,
    ``s' = (j - R_beta, d, sigma_s)`` at ``R = -j`` -- the placement the bond
    transport reads (``G_cd(R + R_beta - R_alpha)``)."""
    from hwave.solver import backend as _bk, matsubara as _ms
    nx, ny, nz = shape
    nvol = nx * ny * nz
    nmat = G_kw.shape[1]
    nd = norb * norb
    B = len(delta_r)
    g_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw[0].reshape(nmat, nvol * nd), axis=0).reshape(nmat, nx, ny, nz, nd),
        axes=(1, 2, 3), workers=1).reshape(nmat, nvol, norb, norb)
    rev = np.zeros_like(g_rt)
    for t in range(nmat):
        for i in range(nvol):
            xi, yi, zi = _coords(i, shape)
            rev[t, i] = g_rt[(-t) % nmat, _site(-xi, -yi, -zi, shape)]
    sgn = -np.ones(nmat)
    sgn[0] = 1.0

    def G(P, Q, table, flip):
        if (P // norb) % 2 != (Q // norb) % 2:
            return None
        xi = _coords(_site_of(P, norb), shape)
        xj = _coords(_site_of(Q, norb), shape)
        d = [(b - a) if flip else (a - b) for a, b in zip(xi, xj)]
        return table[:, _site(*d, shape=shape), P % norb, Q % norb]

    Gam = _gamma_entries(recs)
    v1, v2 = {}, {}
    for (p, q, r, s), (x, _off) in Gam.items():
        # Gamma_{p q, r s} keyed by (p, s); read as Gamma_{r' s', p' q'} keyed by (p', s')
        v1.setdefault((p, s), []).append((q, r, x))
        v2.setdefault((r, q), []).append((p, s, x))
    W = np.zeros((nmat, nvol, B * nd, B * nd), complex)
    for j2 in range(nvol):
        xj, yj, zj = _coords(j2, shape)
        rsite = _site(-xj, -yj, -zj, shape)
        for al, be in itertools.product(range(B), repeat=2):
            Ra, Rb = delta_r[al], delta_r[be]
            ssite = _site(-Ra[0], -Ra[1], -Ra[2], shape)
            spsite = _site(xj - Rb[0], yj - Rb[1], zj - Rb[2], shape)
            for a, b, c, d in itertools.product(range(norb), repeat=4):
                acc = np.zeros(nmat, complex)
                for ss in (UP, DN):
                    p, s = _g(0, a, UP, norb), _g(ssite, c, ss, norb)
                    pp, sp = _g(j2, b, UP, norb), _g(spsite, d, ss, norb)
                    for (q, r, x1) in v1.get((p, s), []):
                        for (rp, qp, x2) in v2.get((pp, sp), []):
                            gf = G(r, rp, g_rt, False)
                            gr = G(qp, q, rev, True)
                            if gf is not None and gr is not None:
                                acc += 0.5 * x1 * x2 * gf * gr
                W[:, rsite, al * nd + c * norb + a, be * nd + d * norb + b] += -acc * sgn / beta
    Wq = _bk.spatial_fftn(W.reshape(nmat, nx, ny, nz, (B * nd) ** 2), axes=(1, 2, 3), workers=1)
    return _ms.tau_to_boson(Wq.reshape(nmat, nvol * (B * nd) ** 2), axis=0).reshape(
        nmat, nvol, B * nd, B * nd)


class TestFactorisedKernel(unittest.TestCase):

    def test_factors_equal_the_enumerated_bond_kernel(self):
        from hwave.solver import bubble
        from hwave.solver.second_order import build_bond_w2_factors
        s0, _ = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0), (0, 0, 0, 2, 2, 1.0, 0.0)]})
        G = _bare_green(s0, _BETA)
        nd = NORB * NORB
        for name, rows in CASES.items():
            with self.subTest(case=name):
                delta_r = _delta_r(rows)
                want = _enumerated_w2(G, _BETA, oracle_records(rows, NORB), delta_r)
                _, split = _split_for(rows)
                fac = build_bond_w2_factors(split, NORB, delta_r, _SHAPE)
                B = len(delta_r)
                view = types.SimpleNamespace(n_channels=B, delta_r=np.asarray(delta_r))
                chib = np.zeros((G.shape[1], G.shape[2], B * nd, B * nd), complex)
                for (m, mp), blk in bubble._iter_bond_dynamic(G, None, _BETA, view,
                                                              spatial_shape=_SHAPE):
                    chib[:, :, m * nd:(m + 1) * nd, mp * nd:(mp + 1) * nd] = blk
                got = sum(0.5 * A[None] @ chib @ Bf[None] for _, A, Bf in fac)
                ref = np.abs(want).max()
                self.assertGreater(ref, 1e-3)
                self.assertLess(np.abs(got - want).max(), 1e-12 * ref)


def _write_inputs(d, lam, rows_by_type):
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(sopt._IN2, f), d)
    rows_by_type = dict(rows_by_type)
    rows_by_type.setdefault("CoulombIntra", [(0, 0, 0, 1, 1, 0.0, 0.0), (0, 0, 0, 2, 2, 0.0, 0.0)])
    inter = {}
    for t, rows in rows_by_type.items():
        sopt._write_wan(os.path.join(d, t + ".dat"), t, NORB,
                        [tuple(r[:5]) + (r[5] * lam, r[6] * lam) for r in rows])
        inter[t] = t + ".dat"
    return inter


def _one_map(lam, rows, legacy=False):
    """Sigma_fluct of ONE production bond map at coupling scale ``lam`` from the
    frozen bare G (``legacy``: the mixed-block rule instead of the factors)."""
    from hwave.solver import flex_bond
    with tempfile.TemporaryDirectory() as d:
        s, gi = sopt._solver(d, _write_inputs(d, lam, rows), True)
        s._phase_b_reset(gi)
        s._phase_b_preflight(gi)
        s._calc_epsilon_k({})
        beta = 1.0 / s.T
        mu = s.mu_value
        _, tail = s._calc_green(beta, mu)
        nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
        nd = norb * norb
        G = s._calc_dressed_green(beta, mu, np.zeros((1, nmat, nvol, norb, norb), complex))
        B = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, tail, beta, s._bond_view, sopt._SHAPE, 1)
            with flex_bond.BondDeviceContext.for_view(np, s._bond_S, s._bond_C,
                                                      s._bond_view, norb) as dev:
                flex_bond.dress_and_build_w(
                    store, dev, nb=nmat, output_full=False, nmat=nmat, nvol=nvol, nd=nd,
                    spatial_shape=sopt._SHAPE, factors=s._second_order_factors,
                    bond_w2=None if legacy else s._bond_w2_on(np))
            sig = flex_bond.calc_self_energy_bond(store, G, beta, s._bond_view, sopt._SHAPE,
                                                  norb, 1)
    return sig[0], G


def _second_order_coefficient(rows, legacy=False, h=3e-3):
    """O(lambda^2) of Sigma_fluct: symmetric stencil, Richardson over h, h/2, h/4."""
    vals = []
    for hh in (h, h / 2, h / 4):
        fp, G = _one_map(hh, rows, legacy)
        fm, _ = _one_map(-hh, rows, legacy)
        f0, _ = _one_map(0.0, rows, legacy)
        vals.append((fp + fm - 2 * f0) / (2 * hh * hh))
    r1 = (4 * vals[1] - vals[0]) / 3
    r2 = (4 * vals[2] - vals[1]) / 3
    return (16 * r2 - r1) / 15, G


class TestSecondOrderOfTheBondMap(unittest.TestCase):

    def test_exact_second_order(self):
        for name, rows in CASES.items():
            with self.subTest(case=name):
                c, G = _second_order_coefficient(rows)
                orc = oracle_sigma2(G, 1.0 / sopt._T, oracle_records(rows, NORB), NORB, "exact")[0]
                self.assertGreater(np.abs(orc).max(), 1e-3)
                self.assertLess(sopt._rel(c, orc), 1e-10, name)

    def test_the_previous_mixed_block_rule_misses_the_cross_terms(self):
        """Falsification: the mixed-block rule 1/4 (S chibar S + C chibar C)
        with zero bond-bond blocks gets the Heisenberg bond's second order
        wrong (it carries only half of the density x Exchange cross term)."""
        rows = CASES["I+X (Heisenberg)"]
        c, G = _second_order_coefficient(rows, legacy=True)
        orc = oracle_sigma2(G, 1.0 / sopt._T, oracle_records(rows, NORB), NORB, "exact")[0]
        self.assertGreater(sopt._rel(c, orc), 1e-4)


class TestSU2Criterion(unittest.TestCase):

    def test_exact_diagonalization_pins_the_ratio(self):
        """A bond commutes with the total S- exactly when J_E = J_H - 2 J_I."""
        from tests import ed_oracle_util as eu
        from tests.test_bond_transverse_ed import (_terms_exchange_offsite, _terms_hund_offsite,
                                                   _terms_ising_offsite)
        fx = eu.EDFixture(3, 1, {(0, 0): -1.0}, [0.0], 0.5, 0.0)
        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        Sm = sum(CD[fx.mode(j, 0, 1)] @ C[fx.mode(j, 0, 0)] for j in range(fx.L))
        for jI, jH, jE in ((0.2, 0.0, -0.4), (0.0, 0.3, 0.3), (0.1, 0.3, 0.1), (0.2, 0.0, -0.2),
                           (0.0, 0.3, -0.3)):
            terms = (_terms_ising_offsite(fx, 0, 0, 1, jI) + _terms_hund_offsite(fx, 0, 0, 1, jH)
                     + _terms_exchange_offsite(fx, 0, 0, 1, jE))
            H = eu.h_int_from_terms(fx, terms)
            dev = np.abs(H @ Sm - Sm @ H).max()
            if abs(jE - (jH - 2 * jI)) < 1e-12:
                self.assertLess(dev, 1e-14, (jI, jH, jE))
            else:
                self.assertGreater(dev, 0.1, (jI, jH, jE))

    def test_criterion_on_the_split(self):
        import hwave.solver.flex as flex_mod
        bond = [(1, 0, 0, 1, 2), (-1, 0, 0, 2, 1)]
        for rows, n_bad in (
                ({"Ising": [b + (0.2, 0.0) for b in bond],
                  "Exchange": [b + (-0.4, 0.0) for b in bond]}, 0),
                ({"Hund": [b + (0.3, 0.0) for b in bond],
                  "Exchange": [b + (0.3, 0.0) for b in bond]}, 0),
                ({"Exchange": [b + (0.3, 0.0) for b in bond]}, 2),
                ({"Ising": [b + (0.2, 0.0) for b in bond]}, 2),
                ({"CoulombInter": [b + (0.5, 0.0) for b in bond]}, 0)):
            with self.subTest(rows=sorted(rows)):
                from hwave.solver import bond_channels as bc
                _, split = _split_for(rows)
                topo = bc.resolve_bond_topology(split.whole_tbl, np.eye(3), NORB,
                                                active_types=bc._LONGITUDINAL_BOND_TYPES)
                self.assertEqual(len(flex_mod.FLEX._non_su2_bonds(topo)), n_bad)


class TestPolicy(unittest.TestCase):

    def _solver(self, rows, gate, extra=None):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s, gi = sopt._solver(d, _write_inputs(d, 1.0, rows), gate)
        for k, v in (extra or {}).items():
            setattr(s, k, v)
        return s, gi

    def test_without_the_bond_channels_an_offsite_exchange_is_refused(self):
        for t in ("Exchange", "PairHop"):
            with self.subTest(type=t):
                s, gi = self._solver({t: [r[:5] + (0.2, 0.0) for r in EXCH]}, gate=False)
                with self.assertRaisesRegex(ValueError, "off-site {} .*longitudinal_bond_channels "
                                                        "= true".format(t)):
                    s._validate_phase_b_interactions()

    def test_gate_accepts_and_warns_only_off_the_su2_ratio(self):
        log = logging.getLogger("hwave.solver.flex")
        for rows, warns in (({"Ising": [r[:5] + (0.1, 0.0) for r in EXCH],
                              "Exchange": [r[:5] + (-0.2, 0.0) for r in EXCH]}, False),
                            ({"Exchange": [r[:5] + (0.2, 0.0) for r in EXCH]}, True)):
            with self.subTest(rows=sorted(rows)):
                s, gi = self._solver(rows, gate=True)
                records = []
                h = logging.Handler(level=logging.WARNING)
                h.emit = lambda rec: records.append(rec.getMessage())
                log.addHandler(h)
                try:
                    with tempfile.TemporaryDirectory() as out:
                        s.solve(gi, out)
                finally:
                    log.removeHandler(h)
                self.assertTrue(np.all(np.isfinite(np.asarray(gi["green"]))))
                self.assertEqual(any("SU(2)" in m for m in records), warns, records)
                if warns:      # 1-based orbitals, as in the interaction file
                    self.assertTrue(any("orbitals=(1, 2)" in m or "orbitals=(2, 1)" in m
                                        for m in records), records)

    def test_degenerate_onsite_row_is_refused(self):
        from hwave.solver.second_order import DegenerateRowError
        from hwave.solver.second_order import build_bond_w2_factors
        rows = {"Exchange": [r[:5] + (0.2, 0.0) for r in EXCH] + [(0, 0, 0, 1, 1, 0.3, 0.0)]}
        # the solver refuses the row at construction, where the q-only
        # second-order kernel is compiled ...
        with self.assertRaises(DegenerateRowError):
            self._solver(rows, gate=True)
        # ... and the bond-enlarged factors refuse it themselves (the split
        # is taken from a solver whose q-only compile is skipped)
        from unittest import mock
        from hwave.solver import second_order
        stub = types.SimpleNamespace(triples=(), A_off=None, vpair=None, nbytes=0)
        with mock.patch.object(second_order, "build_factors", return_value=stub):
            _, split = _split_for(rows)
        with self.assertRaises(DegenerateRowError):
            build_bond_w2_factors(split, NORB, _delta_r(rows), _SHAPE)

    def test_memory_cap_refuses_before_the_factors_are_allocated(self):
        from unittest import mock
        from hwave.solver import second_order
        s, gi = self._solver({"Exchange": [r[:5] + (0.2, 0.0) for r in EXCH]}, gate=True,
                             extra={"longitudinal_bond_memory_cap_gb": 1e-9})
        s._phase_b_reset(gi)
        with mock.patch.object(second_order, "build_bond_w2_factors",
                               side_effect=AssertionError("allocated before admission")):
            with self.assertRaisesRegex(ValueError, "GiB"):
                s._phase_b_preflight(gi)

    def test_reverse_declared_su2_bond_is_not_flagged(self):
        """Ising declared on one direction, Exchange on the reverse partner:
        the closed topology is SU(2) symmetric."""
        import hwave.solver.flex as flex_mod
        from hwave.solver import bond_channels as bc
        _, split = _split_for({"Ising": [(1, 0, 0, 1, 2, 0.2, 0.0), (-1, 0, 0, 2, 1, 0.2, 0.0)],
                               "Exchange": [(-1, 0, 0, 2, 1, -0.4, 0.0), (1, 0, 0, 1, 2, -0.4, 0.0)]})
        topo = bc.resolve_bond_topology(split.whole_tbl, np.eye(3), NORB,
                                        active_types=bc._LONGITUDINAL_BOND_TYPES)
        self.assertEqual(flex_mod.FLEX._non_su2_bonds(topo), [])

    def test_truncated_pairlift_shell_is_refused(self):
        rows = {"Exchange": [r[:5] + (0.2, 0.0) for r in EXCH],
                "PairLift": [(2, 0, 0, 1, 2, 0.3, 0.0), (-2, 0, 0, 2, 1, 0.3, 0.0)]}
        s, gi = self._solver(rows, gate=True, extra={"longitudinal_bond_max_shells": 1})
        s._phase_b_reset(gi)
        with self.assertRaisesRegex(ValueError, "PairLift row at irvec=.*longitudinal_bond_max_shells"):
            s._phase_b_preflight(gi)

    def test_no_qonly_offsite_factors_with_bond_exchange(self):
        s, gi = self._solver({"Exchange": [r[:5] + (0.2, 0.0) for r in EXCH]}, gate=True)
        f = s._second_order_factors
        self.assertIsNone(f.A_off)
        self.assertIsNone(f.vpair)

    def test_operation_count_includes_the_exact_second_order(self):
        from hwave.solver.flex_bond import dressing_ops
        self.assertEqual(dressing_ops(8, 16, 12), 2.0 * 8 * 16 * 12 ** 3)
        self.assertEqual(dressing_ops(8, 16, 12, w2_triples=3), 2.0 * 4 * 8 * 16 * 12 ** 3)

    def test_density_only_bond_map_keeps_the_previous_second_order(self):
        s, gi = self._solver({"CoulombInter": [r[:5] + (0.3, 0.0) for r in EXCH]}, gate=True)
        s._phase_b_reset(gi)
        s._phase_b_preflight(gi)
        self.assertIsNone(s._bond_w2)


if __name__ == "__main__":
    unittest.main()
