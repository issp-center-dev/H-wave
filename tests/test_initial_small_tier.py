"""The small tier of ``hwave.initial``: the eigen-solve of the hints
(``hints.solve_small``), the hint eigenpairs on the para bands
(``hints.hint_eigenpairs``), the instability rows and the orbital family
of the enumerator, the interim tier rule, and the para reference of the
command line (seeds, labels and collect share it).

Physics anchors:

* V1a (square 16x16, half filling, U = 4, T = 0.02): the leading spin
  eigenvalue is at (1/2, 1/2);
* V2b (triangular 12x12, Ncond = 88, T = 0.1, U = 2): one orbital, so the
  spin eigenvalue is ``U chi0`` with the bare Lindhard ``chi0`` per spin;
  its maximum is at K or K' and K is a strict local maximum on the index
  grid; the negative controls (half-filled triangular, half-filled square)
  behave as recorded.
"""
import contextlib
import io
import json
import os
import shutil
import tempfile
import unittest
from fractions import Fraction
from unittest import mock

import numpy as np

from tests import initial_fixtures
from tests.heavy_tests import heavy


def _facts(toml, subshape=None):
    from hwave.initial.facts import model_facts
    return model_facts(toml, subshape=subshape)


def _quiet():
    return contextlib.redirect_stdout(io.StringIO())


class _Tmp(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="hwave-small-")
        self.n = 0

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def dest(self, name=None):
        self.n += 1
        return os.path.join(self.tmp, name or "m{}".format(self.n))


def _random_psd(rng, n, rank=None):
    rank = n if rank is None else rank
    a = rng.standard_normal((n, rank)) + 1j * rng.standard_normal((n, rank))
    return a @ a.conj().T / n


def _random_herm(rng, n):
    a = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    return (a + a.conj().T) / 2


# ---------------------------------------------------------------------------
# the eigen-solve
# ---------------------------------------------------------------------------

class TestSolveSmall(unittest.TestCase):

    def test_symmetric_form_equals_K(self):
        """Every pair of the symmetric form is an eigenpair of K = -chi0
        Gamma (1e-12), the eigenvalues are the positive eigenvalues of K, and
        ||x|| = 1; also for a rank-deficient chi0."""
        from hwave.initial.hints import solve_small
        rng = np.random.default_rng(5)
        for rank in (8, 5):
            with self.subTest(rank=rank):
                chi = _random_psd(rng, 8, rank)
                gam = _random_herm(rng, 8)
                res = solve_small(chi, gam, top_k=8)
                K = -chi @ gam
                w = np.linalg.eigvals(K)
                pos = np.sort(w.real[(w.real > 1e-9) & (np.abs(w.imag) < 1e-9)])[::-1]
                lams = np.array([lam for lam, _ in res.pairs])
                np.testing.assert_allclose(lams, pos, atol=1e-12, rtol=0)
                for lam, x in res.pairs:
                    self.assertAlmostEqual(np.linalg.norm(x), 1.0, delta=1e-13)
                    self.assertLessEqual(np.linalg.norm(K @ x - lam * x),
                                         1e-12 * max(1.0, abs(lam)))
                self.assertTrue(res.gamma_self_adjoint)
                self.assertEqual(res.lambda_imag, (0.0,) * len(res.pairs))
                self.assertEqual(res.q_status, "ok" if len(pos) >= 8 else
                                 "fewer than top_k positive eigenvalues")

    def test_negative_only_no_positive(self):
        """A positive-definite Gamma gives a negative K: no pair, never a
        null vector."""
        from hwave.initial.hints import solve_small
        rng = np.random.default_rng(1)
        chi = _random_psd(rng, 6)
        g = _random_psd(rng, 6) + np.eye(6)
        res = solve_small(chi, g)
        self.assertEqual(res.q_status, "no positive eigenvalue")
        self.assertEqual(res.pairs, ())
        self.assertEqual(res.n_positive, 0)
        # the null space of chi0 coupled by Gamma (a defective lambda = 0
        # block, chi0 = diag(1, 0), Gamma = sigma^x) gives no pair either
        res = solve_small(np.diag([1.0, 0.0]), np.array([[0, 1], [1, 0]], dtype=complex))
        self.assertEqual(res.q_status, "no positive eigenvalue")
        # chi0 = 0
        res = solve_small(np.zeros((3, 3)), -np.eye(3))
        self.assertEqual(res.q_status, "no positive eigenvalue")

    def test_fewer_than_topk(self):
        from hwave.initial.hints import solve_small
        res = solve_small(np.diag([1.0, 0.5, 0.2]), np.diag([-2.0, 1.0, 1.0]), top_k=3)
        self.assertEqual(res.q_status, "fewer than top_k positive eigenvalues")
        self.assertEqual(len(res.pairs), 1)
        self.assertAlmostEqual(res.pairs[0][0], 2.0, delta=1e-14)
        np.testing.assert_allclose(np.abs(res.pairs[0][1]), [1, 0, 0], atol=1e-14)
        res = solve_small(np.diag([1.0, 0.5, 0.2]), np.diag([-2.0, -1.0, -1.0]), top_k=2)
        self.assertEqual(res.q_status, "ok")
        self.assertEqual([round(lam, 12) for lam, _ in res.pairs], [2.0, 0.5])
        with self.assertRaises(ValueError):
            solve_small(np.eye(2), np.eye(2), top_k=0)

    def test_sector_merge_ties_order(self):
        """Sectors merged before top_k: sorted by lambda, ties by the sector
        order (charge, S^z) and then the index; a degenerate set crossing the
        cutoff is kept whole."""
        from hwave.initial.hints import solve_small
        chi = {"charge": np.eye(1), "S^z": np.eye(2)}
        gam = {"charge": -2.0 * np.eye(1), "S^z": np.diag([-2.0, -1.0])}
        res = solve_small(chi, gam, top_k=1)
        self.assertEqual(res.sectors, ("charge", "S^z"))
        self.assertEqual([lam for lam, _ in res.pairs], [2.0, 2.0])
        self.assertEqual(res.groups, (0, 0))
        self.assertEqual(res.n_positive, 3)
        self.assertEqual(res.q_status, "ok")
        res = solve_small(chi, gam, top_k=3)
        self.assertEqual(res.sectors, ("charge", "S^z", "S^z"))
        self.assertEqual(res.groups, (0, 0, 1))
        # S^z first when its eigenvalue is larger
        gam2 = {"charge": -1.0 * np.eye(1), "S^z": np.diag([-3.0, -1.0])}
        res = solve_small(chi, gam2, top_k=2)
        self.assertEqual(res.sectors, ("S^z", "charge", "S^z"))
        with self.assertRaises(ValueError):
            solve_small({"charge": np.eye(1)}, {"S^z": np.eye(1)})

    def test_degenerate_canonical_basis(self):
        """A degenerate set gets the canonical basis of its span: the same
        vectors whatever basis the eigensolver returns, a single pair a
        fixed phase (real positive pivot)."""
        from hwave.initial.hints import solve_small
        rng = np.random.default_rng(7)
        q, _ = np.linalg.qr(rng.standard_normal((4, 4)) + 1j * rng.standard_normal((4, 4)))
        d = np.diag([2.0, 2.0, 1.0, -1.0])
        th = 0.7
        r = np.array([[np.cos(th), -np.sin(th) * 1j], [-np.sin(th) * 1j, np.cos(th)]])
        q2 = q.copy()
        q2[:, :2] = q[:, :2] @ r
        a = solve_small(np.eye(4), -(q @ d @ q.conj().T), top_k=3)
        b = solve_small(np.eye(4), -(q2 @ d @ q2.conj().T), top_k=3)
        for (la, xa), (lb, xb) in zip(a.pairs, b.pairs):
            self.assertAlmostEqual(la, lb, delta=1e-12)
            np.testing.assert_allclose(xa, xb, atol=1e-10)
        self.assertEqual(a.groups, (0, 0, 1))
        # the phase of a single pair: its largest entry is real positive
        x = a.pairs[2][1]
        k = int(np.argmax(np.abs(x)))
        self.assertAlmostEqual(x[k].imag, 0.0, delta=1e-14)
        self.assertGreater(x[k].real, 0.0)
        c = solve_small(np.eye(4), -(q @ d @ q.conj().T) * np.exp(0j), top_k=3)
        np.testing.assert_allclose(c.pairs[2][1], x, atol=1e-12)

    def test_phi0_branch_at_negative_real(self):
        """A real negative tr AA gives phi0 = -pi/2 whatever the sign of a
        round-off imaginary part (bond patterns at q = -q, where e^{iq.r} =
        -1 carries round-off), so the phase in a name is reproducible."""
        import math
        from hwave.initial.patterns import Entry, phi0
        q = (Fraction(1, 2), 0, 0)
        for eps in (1e-15, -1e-15, 0.0):
            A = (Entry(0, 0, (1, 0, 0), c0=1.0 + eps * 1j),
                 Entry(0, 0, (-1, 0, 0), c0=1.0))
            self.assertEqual(phi0(A, q), -math.pi / 2)

    def test_not_self_adjoint_general_solve(self):
        """A non-Hermitian Gamma: a warning, the general solve of K with the
        imaginary parts reported, gamma_self_adjoint false."""
        from hwave.initial.hints import solve_small
        g = -np.array([[2.0, 0.1], [0.0, 1.0]])
        with self.assertLogs("qlms.initial.hints", level="WARNING"):
            res = solve_small(np.eye(2), g)
        self.assertFalse(res.gamma_self_adjoint)
        self.assertTrue(res.warnings)
        self.assertEqual([round(lam, 12) for lam, _ in res.pairs], [2.0, 1.0])
        self.assertEqual(res.lambda_imag, (0.0, 0.0))
        K = -g
        for lam, x in res.pairs:
            self.assertLessEqual(np.linalg.norm(K @ x - lam * x), 1e-12)


# ---------------------------------------------------------------------------
# the hint eigenpairs on the para bands
# ---------------------------------------------------------------------------

def _para(toml):
    from hwave.initial import hints
    f = _facts(toml)
    with _quiet():
        para = hints.para_hf(f, hints.t_hint_bare(f.problem))
    return f, para


def _spin_lambda(eigen):
    """lambda_spin(q): the largest S^z eigenvalue of each q (0 when none)."""
    out = {}
    for q, res in eigen.items():
        out[q] = max([lam for (lam, _x), s in zip(res.pairs, res.sectors) if s == "S^z"],
                     default=0.0)
    return out


def _P(lam, m1, m2, L):
    """lambda(m1, m2) exceeds its 8 index-grid neighbours and lambda(Gamma)."""
    def at(a, b):
        return lam[(Fraction(a % L, L), Fraction(b % L, L), Fraction(0))]
    v = at(m1, m2)
    nb = [at(m1 + i, m2 + j) for i in (-1, 0, 1) for j in (-1, 0, 1) if (i, j) != (0, 0)]
    return v > max(nb) and v > at(0, 0)


class TestHintEigenpairs(_Tmp):

    def eigen(self, toml):
        from hwave.initial import cells, hints
        f, para = _para(toml)
        self.assertEqual(para.status, "converged", para.skip_reason)
        P = hints.pair_space(f)
        return f, para, P, hints.hint_eigenpairs(f, para, P,
                                                 cells.commensurate_q(f.cellshape))

    def test_v1a_i2_leading_q(self):
        """V1a: the leading spin eigenvalue is at (1/2, 1/2), and it leads the
        hint ranking."""
        from hwave.initial.enumerate import hint_pairs
        toml = initial_fixtures.square(self.dest(), 16, U=4.0, T=0.02, subshape=(2, 2, 1),
                                       param={"IterationMax": 2000})
        _f, _para_, P, eigen = self.eigen(toml)
        self.assertEqual(P.sectors, ("charge", "S^z"))
        lam = _spin_lambda(eigen)
        best = max(lam, key=lam.get)
        self.assertEqual(best, (Fraction(1, 2), Fraction(1, 2), Fraction(0)))
        top = hint_pairs(eigen)[0]
        self.assertEqual((top.q, top.sector, top.rank), (best, "S^z", 0))
        self.assertEqual(len(eigen), 256)

    def test_v2b_property(self):
        """V2b: argmax lambda_spin in {K, K'}; P(K); the frozen values
        (U chi0 per spin at K, Gamma, M) and lambda / chi0 = U on the whole
        mesh; the negative controls."""
        from hwave.initial import hints
        L = 12
        toml = initial_fixtures.triangular(self.dest(), L, U=2.0, T=0.1, ncond=88,
                                           subshape=(1, 1, 1), param={"IterationMax": 2000})
        f, para, P, eigen = self.eigen(toml)
        lam = _spin_lambda(eigen)
        K = (Fraction(2, 3), Fraction(1, 3), Fraction(0))
        Kp = (Fraction(1, 3), Fraction(2, 3), Fraction(0))
        best = max(lam.values())
        self.assertIn(max(lam, key=lam.get), (K, Kp))
        self.assertAlmostEqual(lam[K], best, delta=1e-12)
        self.assertAlmostEqual(lam[Kp], best, delta=1e-12)
        self.assertTrue(_P(lam, 8, 4, L))
        frozen = {K: 2.0 * 0.2269, (Fraction(0),) * 3: 2.0 * 0.2023,
                  (Fraction(1, 2), Fraction(0), Fraction(0)): 2.0 * 0.1699}
        for q, v in frozen.items():
            self.assertAlmostEqual(lam[q], v, delta=2e-4)
        # one orbital: lambda_spin(q) = U chi0(q) per spin, at every q
        up = int(np.nonzero((P.labels == [0, 0, 0, 0, 0, 0, 0]).all(axis=1))[0][0])
        for q in eigen:
            chi = hints.lindhard(para.bands, q, P)
            self.assertAlmostEqual(lam[q], 2.0 * float(chi[up, up].real), delta=1e-10)
        # the leading pair at K is a spin-z pattern on orbital 0
        amp = eigen[K].amplitudes[0]
        self.assertEqual([(e.orb_a, e.orb_b, e.r) for e in amp], [(0, 0, (0, 0, 0))])
        self.assertEqual((amp[0].c0, amp[0].cx, amp[0].cy), (0, 0, 0))
        # negative controls (same mesh, T 0.1)
        half = initial_fixtures.triangular(self.dest(), L, U=2.0, T=0.1, ncond=144,
                                           subshape=(1, 1, 1), param={"IterationMax": 2000})
        self.assertFalse(_P(_spin_lambda(self.eigen(half)[3]), 8, 4, L))
        sq = initial_fixtures.square(self.dest(), L, U=2.0, T=0.1, ncond=144,
                                     subshape=(1, 1, 1), param={"IterationMax": 2000})
        lam_sq = _spin_lambda(self.eigen(sq)[3])
        self.assertFalse(_P(lam_sq, 8, 4, L))
        self.assertTrue(_P(lam_sq, 6, 6, L))

    def test_hint_record_fields(self):
        """The vertex fields of the hint record: every mesh q keyed by its q
        string, the method, the self-adjointness."""
        from hwave.initial import hints
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2, subshape=(2, 2, 1),
                                       param={"IterationMax": 2000})
        f, para, _P_, eigen = self.eigen(toml)
        v = hints.vertex_fields(eigen)
        self.assertEqual(v["method"], "lindhard_static")
        self.assertIs(v["gamma_self_adjoint"], True)
        self.assertEqual(len(v["q_set"]), 16)
        self.assertEqual(v["q_set"][0], "0,0,0")
        self.assertEqual(set(v["q_status"]), set(v["q_set"]))
        self.assertTrue(set(v["q_status"].values())
                        <= {"ok", "fewer than top_k positive eigenvalues",
                            "no positive eigenvalue"})
        rec = hints.hint_record(hints.t_hint_bare(f.problem), para, vertex=v)
        self.assertEqual(rec["status"], "computed")

    def test_amplitude_entries(self):
        """The conjugate field on the labels becomes PatternSpec entries:
        X[s, t] = O_m per bond slot, c^mu = tr(sigma^mu X)."""
        from hwave.initial import hints
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2, subshape=(1, 1, 1))
        P = hints.pair_space(_facts(toml))
        O = np.zeros(P.size, dtype=complex)
        lab = [tuple(x) for x in P.labels.tolist()]
        O[lab.index((0, 0, 0, 0, 0, 0, 0))] = 0.5
        O[lab.index((0, 1, 0, 1, 0, 0, 0))] = -0.5
        O[lab.index((0, 0, 0, 1, 0, 0, 0))] = 0.25j
        e, = hints.amplitude_entries(P, O)
        self.assertEqual((e.orb_a, e.orb_b, e.r), (0, 0, (0, 0, 0)))
        self.assertAlmostEqual(e.c0, 0.0)
        self.assertAlmostEqual(e.cz, 1.0)
        # X[0, 1] = 0.25i: c^x = tr(sigma^x X) = 0.25i, c^y = tr(sigma^y X)
        # = i * 0.25i = -0.25
        self.assertAlmostEqual(e.cx, 0.25j)
        self.assertAlmostEqual(e.cy, -0.25)
        with self.assertRaises(ValueError):
            hints.amplitude_entries(P, np.zeros(P.size))


# ---------------------------------------------------------------------------
# the enumerator: instability rows, orbital family, tiers
# ---------------------------------------------------------------------------

def _bond_model(dest, *, T=0.1, iteration_max=2000):
    """Square 8x8 with U = 2 and a nearest-neighbour CoulombInter V = 1
    (Fock on), Ncond = 48: hints with bond patterns the baseline lacks."""
    V = {((1, 0, 0), 0, 0): 1.0, ((-1, 0, 0), 0, 0): 1.0,
         ((0, 1, 0), 0, 0): 1.0, ((0, -1, 0), 0, 0): 1.0}
    return initial_fixtures.write_input(
        dest, cellshape=(8, 8, 1), transfer=initial_fixtures.square_transfer(),
        interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0}, "CoulombInter": V},
        param={"Ncond": 48, "T": T, "SubShape": [2, 2, 1], "IterationMax": iteration_max})


def _key(r):
    return (r["run_id"], r["tried"], r["not_tried_reason"], r["representable"])


class TestInstabilityRows(_Tmp):

    def cset(self, toml, **kw):
        from hwave.initial.enumerate import candidate_states
        kw.setdefault("n_random", 0)
        kw.setdefault("tier", "small")
        with _quiet():
            return candidate_states(toml, **kw)

    def test_hint_rows_only_order_never_prune(self):
        """The hints change only the instability rows: every other row is the
        same at another T_hint and in the large tier (one orbital: no orbital
        rows). Every hint pair is listed; realised rows follow the rank; the
        listed rows carry their eigenvalue and reason."""
        toml = _bond_model(self.dest())
        a = self.cset(toml, quota_instability=4)
        b = self.cset(toml, quota_instability=4, t_hint=0.3)
        large = self.cset(toml, tier="large")
        self.assertEqual(a.tier, "small")
        non = [[_key(r) for r in c.rows if r["family"] != "instability"]
               for c in (a, b, large)]
        self.assertEqual(non[0], non[1])
        self.assertEqual(non[0], non[2])
        hint = a.model["hint"]
        self.assertEqual(hint["status"], "computed")
        n_pairs = sum(1 for r in a.rows if r["family"] == "instability")
        inst = [r for r in a.rows if r["family"] == "instability"]
        self.assertGreater(n_pairs, 0)
        # every row has its rank and eigenvalue; every rank appears
        ranks = {r["rank"] for r in inst}
        self.assertEqual(ranks, set(range(max(ranks) + 1)))
        for r in inst:
            self.assertGreater(r["lambda"], 0.0)
            self.assertEqual(r["lambda_imag"], 0.0)
        # the realised rows (tried or not) come in rank order
        realised = [r for r in inst if r["representable"]
                    and r["not_tried_reason"] != "quota: instability"]
        self.assertEqual([r["rank"] for r in realised],
                         sorted(r["rank"] for r in realised))
        tried = [r for r in inst if r["tried"]]
        self.assertTrue(tried)
        self.assertLessEqual(len(tried), 4)
        listed = [r for r in inst if r["not_tried_reason"] == "quota: instability"]
        self.assertTrue(listed)
        self.assertGreaterEqual(min(r["rank"] for r in listed),
                                min(r["rank"] for r in tried))
        unrep = [r for r in inst if r["not_tried_reason"] == "unrepresentable"]
        self.assertTrue(unrep)
        for r in unrep:
            self.assertFalse(r["representable"])
            self.assertEqual(list(r["subshape"]), [2, 2, 1])
            q = [Fraction(x) for x in r["q"]]
            self.assertEqual(list(r["proposed_subshape"]), [x.denominator for x in q])
        # T_hint changes the eigenvalues, so the instability rows may differ,
        # but not the rows of the other families (checked above)
        self.assertNotEqual([r["lambda"] for r in inst][:3],
                            [r["lambda"] for r in b.rows if r["family"] == "instability"][:3])

    def test_quota_filled_after_dedup(self):
        """Duplicates of baseline rows consume no instability quota: more
        pairs are realised until the quota of the cell is reached."""
        toml = _bond_model(self.dest())
        for quota in (2, 6):
            with self.subTest(quota=quota):
                c = self.cset(toml, quota_instability=quota)
                tried = [r for r in c.rows if r["family"] == "instability" and r["tried"]]
                self.assertEqual(len(tried), quota)
                dups = [r for r in c.rows if r["family"] == "instability"
                        and (r["not_tried_reason"] or "").startswith("duplicate of baseline")]
                self.assertTrue(dups)

    def test_spiral_and_phase_rows(self):
        """Under a passed SU(2) check an S^z pair at q not equal to -q gives
        two phase rows and one spiral; at q = -q its -O row is listed."""
        toml = initial_fixtures.chain(self.dest(), 12, U=2.0, T=0.2, ncond=8,
                                      subshape=(3, 1, 1), param={"IterationMax": 2000})
        c = self.cset(toml, quota_instability=12, max_candidates=64)
        self.assertEqual(c.model["su2"]["status"], "passed")
        third = [r for r in c.rows if r["family"] == "instability"
                 and r["q_string"] == "1/3,0,0" and r["channel"] in ("spin-z", "spin-spiral")
                 and r["rank"] == min(x["rank"] for x in c.rows
                                      if x["family"] == "instability"
                                      and x["q_string"] == "1/3,0,0"
                                      and x["channel"] == "spin-z")]
        self.assertEqual(sorted(r["channel"] for r in third),
                         ["spin-spiral", "spin-z", "spin-z"])
        zero = [r for r in c.rows if r["family"] == "instability"
                and r["q_string"] == "0,0,0" and r["representable"]
                and r["not_tried_reason"] != "quota: instability"]
        signs = {r["sign"] for r in zero}
        self.assertEqual(signs, {1, -1})

    def test_orbital_family(self):
        """The small tier adds n_a - n_b per site, uniform and with the
        (q, eta) signs (family orbital, baseline quota); the large tier has
        none."""
        from hwave.initial.enumerate import orbital_pairs
        toml = initial_fixtures.two_orbital(self.dest(), 8, delta=0.5, U=(2.0, 2.0),
                                            T=0.2, subshape=(2, 1, 1),
                                            param={"IterationMax": 2000})
        f = _facts(toml)
        self.assertEqual(orbital_pairs(f), (((0, 1),),))
        small = self.cset(toml, max_candidates=64)
        large = self.cset(toml, tier="large", max_candidates=64)
        orb = [r for r in small.rows if r["family"] == "orbital"]
        self.assertTrue(orb)
        self.assertFalse([r for r in large.rows if r["family"] == "orbital"])
        self.assertEqual({r["q_string"] for r in orb}, {"0,0,0", "1/2,0,0"})
        cand = {c.candidate_id: c for c in small.candidates}
        uni = [cand[r["candidate_id"]] for r in orb
               if r["q_string"] == "0,0,0" and r["sign"] == 1]
        self.assertTrue(uni)
        ent = {(e.orb_a, e.c0) for e in uni[0].pattern.entries}
        self.assertEqual(ent, {(0, 1.0), (1, -1.0)})
        for r in orb:
            self.assertEqual(r["channel"], "charge")
        # the orbital rows count toward the baseline quota
        capped = self.cset(toml, max_candidates=0)
        for r in capped.rows:
            if r["family"] == "orbital":
                self.assertFalse(r["tried"])
                self.assertTrue(r["not_tried_reason"].startswith(("quota: orbital",
                                                                  "duplicate of")))

    def test_spin_orbital_and_fixed_sz(self):
        """The small tier in spin-orbital mode (one sector, physical-orbital
        entries) and with a fixed 2Sz (spin-diagonal pair space, S^z
        patterns only): hints computed, the coverage document validates (its
        su2 record has a null method when the check is not applicable)."""
        from hwave.initial.report import coverage_document
        so = initial_fixtures.so_chain(self.dest(), 8, U=2.0, T=0.2, soc=0.3,
                                       param={"IterationMax": 2000, "SubShape": [2, 1, 1]})
        sz = initial_fixtures.square(self.dest(), 6, U=2.0, T=0.2, subshape=(2, 2, 1),
                                     two_sz=0, param={"IterationMax": 2000})
        for toml, status in ((so, "not applicable (spin-orbital)"),
                             (sz, "not applicable (2Sz fixed)")):
            with self.subTest(status=status):
                c = self.cset(toml)
                doc = coverage_document(c, phase_marker="pre_run")
                self.assertEqual(doc["model"]["su2"]["status"], status)
                self.assertIsNone(doc["model"]["su2"]["method"])
                self.assertEqual(doc["hint"]["status"], "computed")
                inst = [r for r in c.rows if r["family"] == "instability"]
                self.assertTrue(inst)
                cand = {x.candidate_id: x for x in c.candidates}
                for r in inst:
                    for e in cand[r["candidate_id"]].pattern.entries:
                        self.assertLess(max(e.orb_a, e.orb_b), 1)
                        if toml == sz:
                            self.assertEqual((e.cx, e.cy), (0, 0))
                self.assertFalse([r for r in c.rows if r["channel"] == "spin-spiral"])

    def test_large_tier_fixed_sz_coverage(self):
        """The large tier with a fixed 2Sz gives a valid coverage document:
        the spin-rotation record has null max_residual and method (the S_z
        check is in sz_residual)."""
        from hwave.initial.report import coverage_document
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2, subshape=(2, 2, 1),
                                       two_sz=0)
        c = self.cset(toml, tier="large")
        doc = coverage_document(c, phase_marker="pre_run")
        su2 = doc["model"]["su2"]
        self.assertEqual(su2["status"], "not applicable (2Sz fixed)")
        self.assertIsNone(su2["method"])
        self.assertIsNone(su2["max_residual"])
        self.assertIsNotNone(su2["sz_residual"])

    def test_rank_and_ids_independent_of_eigensolver(self):
        """Eigenvalues equal by symmetry (q and -q, lattice images, the
        charge and S^z bond pairs at q = 0, two equivalent orbitals) are
        ordered by their degenerate set and the tie keys, not by round-off:
        random phases on every eigenvector of the eigensolver leave the
        ids, ranks, tried set and patterns unchanged."""
        def orb(dest):
            t = {}
            for a in range(2):
                t[((1, 0, 0), a, a)] = -1.0
                t[((-1, 0, 0), a, a)] = -1.0
            t[((0, 0, 0), 0, 1)] = -0.3
            t[((0, 0, 0), 1, 0)] = -0.3
            return initial_fixtures.write_input(
                dest, cellshape=(8, 1, 1), norb=2, transfer=t,
                interactions={"CoulombIntra": {((0, 0, 0), a, a): 0.5 for a in range(2)},
                              "CoulombInter": {((0, 0, 0), 0, 1): 1.5,
                                               ((0, 0, 0), 1, 0): 1.5}},
                param={"Ncond": 8, "T": 0.2, "SubShape": [4, 1, 1], "IterationMax": 3000})

        def summary(c):
            cand = {x.candidate_id: x for x in c.candidates}
            out = []
            for r in c.rows:
                p = cand[r["candidate_id"]].pattern
                ent = [(e.orb_a, e.orb_b, e.r) + tuple(
                    round(v, 8) + 0.0 for z in (e.c0, e.cx, e.cy, e.cz)
                    for v in (complex(z).real, complex(z).imag)) for e in p.entries]
                out.append((r["run_id"], r["tried"], r["not_tried_reason"], r["rank"],
                            None if p.phase is None else round(p.phase, 8) + 0.0, ent))
            return out

        eigh = np.linalg.eigh
        for name, toml in (("bond", _bond_model(self.dest())), ("orbital", orb(self.dest()))):
            with self.subTest(model=name):
                plain = summary(self.cset(toml, max_candidates=64))
                for seed in (3, 11):
                    rng = np.random.default_rng(seed)

                    def phased(a, *args, **kw):
                        w, v = eigh(a, *args, **kw)
                        if np.iscomplexobj(v):
                            v = v * np.exp(1j * rng.uniform(0, 2 * np.pi, v.shape[-1]))
                        else:
                            v = v * rng.choice([-1.0, 1.0], v.shape[-1])
                        return w, v
                    with mock.patch.object(np.linalg, "eigh", phased):
                        got = summary(self.cset(toml, max_candidates=64))
                    self.assertEqual(got, plain)

    def test_growth_memo_several_subshapes(self):
        """The instability budget grows on one cell while another cell keeps
        its realisations: the ids stay unique and the candidate set equals
        the one without the realisation memo (SubShapes 4x4 then 2x2)."""
        from hwave.initial import enumerate as en
        toml = _bond_model(self.dest())
        orig = en.instability_candidates

        def run(memo_on):
            def ic(*a, **k):
                if not memo_on:
                    k["memo"] = None
                return orig(*a, **k)
            with mock.patch.object(en, "instability_candidates", ic):
                c = self.cset(toml, subshapes=[(4, 4, 1), (2, 2, 1)], max_candidates=64)
            return [(r["run_id"], r["tried"], r["not_tried_reason"], r["rank"])
                    for r in c.rows]

        with_memo = run(True)
        ids = [r[0] for r in with_memo]
        self.assertEqual(len(ids), len(set(ids)))
        self.assertEqual(with_memo, run(False))

    def test_degenerate_clusters(self):
        from hwave.initial.hints import degenerate_clusters
        self.assertEqual(degenerate_clusters([1.0, 2.0, 2.0 + 1e-12, 1.0 - 1e-13, 0.5]),
                         [1, 0, 0, 1, 2])
        self.assertEqual(degenerate_clusters([]), [])

    def test_para_not_converged_skips_hints(self):
        """IterationMax = 1: the para step does not converge (the Fock term
        of V renormalises the bonds), so the hints are skipped with its
        reason and no instability row is made; the other rows are those of
        a converged para step."""
        toml = _bond_model(self.dest(), iteration_max=1)
        c = self.cset(toml)
        hint = c.model["hint"]
        self.assertEqual(hint["status"], "skipped")
        self.assertIn("did not converge", hint["skip_reason"])
        self.assertFalse([r for r in c.rows if r["family"] == "instability"])
        ok = self.cset(_bond_model(self.dest()))
        self.assertEqual([_key(r) for r in c.rows],
                         [_key(r) for r in ok.rows if r["family"] != "instability"])

    def test_para_argument_checks(self):
        from hwave.initial import hints
        from hwave.initial.enumerate import candidate_states
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2, subshape=(2, 2, 1),
                                       param={"IterationMax": 2000})
        with self.assertRaises(ValueError):
            candidate_states(toml, tier="large", para_check=True)
        f = _facts(toml)
        with _quiet():
            para = hints.para_hf(f, hints.t_hint_bare(f.problem, t_user=0.5))
        with self.assertRaises(ValueError):
            candidate_states(toml, tier="large", para=para)
        # a para of another T_hint is refused
        with self.assertRaisesRegex(ValueError, "T_hint"):
            candidate_states(toml, tier="small", para=para, n_random=0)
        # a matching one is used as given
        c = candidate_states(toml, tier="small", para=para, t_hint=0.5, n_random=0)
        self.assertEqual(c.model["hint"]["status"], "computed")
        self.assertEqual(c.model["hint"]["T_hint"], 0.5)


class TestTierRule(_Tmp):

    def choose(self, toml, tier="auto", nbytes=1 << 40, **kw):
        from hwave.initial import enumerate as en
        from hwave.solver import initial_facade
        f = _facts(toml, subshape=(1, 1, 1))
        mem = initial_facade.HostMemory(nbytes, "test" if nbytes is not None else "none")
        with mock.patch.object(initial_facade, "available_host_bytes", return_value=mem):
            return en.choose_tier(f, tier, **kw)

    def test_interim_auto_rule(self):
        """auto depends on the model only: small within the fixed size limits
        (|P|, N_q, N_q |P|^3), large beyond one of them, the same tier for
        every host memory, missing included; the host memory only refuses a
        small tier that does not fit."""
        from hwave.initial import enumerate as en
        from hwave.initial.report import INTERIM_TIER_RULE
        toml = initial_fixtures.square(self.dest(), 8, U=2.0, T=0.2, subshape=(2, 2, 1))
        small = self.choose(toml)
        self.assertEqual((small.name, small.chosen_by), ("small", "auto (interim rule)"))
        self.assertEqual(small.warnings, ())
        pre = small.record["preflight"]
        self.assertEqual(pre["small_bytes"],
                         16 * (4 ** 2 + 64 * 3 * 4 + 64 * 2 ** 2 + 64 * 2 * 4))
        self.assertEqual(small.record["estimate"]["bytes"], pre["small_bytes"])
        self.assertEqual(small.record["rule"], INTERIM_TIER_RULE)
        self.assertEqual(small.record["constants"]["small_max_pairs"], 400)
        self.assertEqual(small.record["constants"]["small_max_q"], 4096)
        need = pre["small_bytes"]
        # the same choice for every host memory, missing included
        for nbytes in (1 << 40, 4 * need - 1, need, None):
            with self.subTest(nbytes=nbytes):
                c = self.choose(toml, nbytes=nbytes)
                self.assertEqual((c.name, c.warnings), ("small", ()))
        # beyond a fixed limit: large on every host
        for name, value in (("SMALL_TIER_MAX_PAIRS", 3), ("SMALL_TIER_MAX_Q", 63),
                            ("SMALL_TIER_MAX_OPERATIONS", 64 * 4 ** 3 - 1)):
            for nbytes in (1 << 40, need, None):
                with self.subTest(limit=name, nbytes=nbytes), \
                        mock.patch.object(en, name, value):
                    c = self.choose(toml, nbytes=nbytes)
                    self.assertEqual((c.name, c.warnings),
                                     ("large", (en.SMALL_TIER_LIMIT_WARNING,)))
                    self.assertIsNone(c.record["estimate"]["bytes"])
        # the host memory only refuses a small tier that does not fit
        for tier in ("auto", "small"):
            with self.assertRaisesRegex(ValueError, "exceeds the available host memory"):
                self.choose(toml, tier, nbytes=need - 1)
        # explicit tiers
        with mock.patch.object(en, "SMALL_TIER_MAX_PAIRS", 3):
            over = self.choose(toml, "small")
        self.assertEqual(over.name, "small")
        self.assertTrue(over.warnings)
        self.assertEqual(self.choose(toml, "small", nbytes=None).name, "small")
        large = self.choose(toml, "large")
        self.assertEqual((large.name, large.chosen_by), ("large", "user"))
        with self.assertRaises(NotImplementedError):
            self.choose(toml, "medium")
        with self.assertRaises(ValueError):
            self.choose(toml, "huge")

    def test_tier_record_in_coverage(self):
        """candidate_states records the choice; the coverage document carries
        the estimate and the interim rule."""
        from hwave.initial.enumerate import candidate_states
        from hwave.initial.report import INTERIM_TIER_RULE, coverage_document
        from hwave.solver import initial_facade
        toml = initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2, subshape=(2, 2, 1),
                                       param={"IterationMax": 2000})
        from hwave.initial import enumerate as en
        for nbytes, pairs, name in ((1 << 40, 400, "small"), (None, 400, "small"),
                                    (None, 3, "large")):
            mem = initial_facade.HostMemory(nbytes, "test")
            with mock.patch.object(initial_facade, "available_host_bytes",
                                   return_value=mem), \
                    mock.patch.object(en, "SMALL_TIER_MAX_PAIRS", pairs), _quiet():
                c = candidate_states(toml, n_random=0)
            self.assertEqual(c.tier, name)
            doc = coverage_document(c, phase_marker="pre_run")
            self.assertEqual(doc["tier"]["name"], name)
            self.assertEqual(doc["tier"]["rule"], INTERIM_TIER_RULE)
            self.assertEqual(doc["hint"]["status"],
                             "computed" if name == "small" else "skipped")
            self.assertIsNotNone(doc["tier"]["preflight"])
            if name == "small":
                self.assertEqual(doc["tier"]["estimate"]["bytes"],
                                 doc["tier"]["preflight"]["small_bytes"])
            else:
                self.assertIn("exceeds the preview limits",
                              " ".join(doc["model"]["warnings"]))


# ---------------------------------------------------------------------------
# the command line: one para reference for seeds, labels and collect
# ---------------------------------------------------------------------------

class TestCliPara(_Tmp):

    def main(self, *argv):
        from hwave.initial import cli
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            rc = cli.main([str(a) for a in argv])
        return rc, out.getvalue(), err.getvalue()

    @heavy
    def test_prepare_small_para_reference(self):
        """prepare --tier small --para-check keeps the para reference in
        WORK/para_hf (SHA-256 in the manifest), seeds on it, labels against
        its h = 0 seed (the para row is "para"), and collect reads the same
        reference; a changed para file is refused."""
        from tests.initial_reference_capture import run_solver
        # the Fock term of V renormalises the bonds of the para reference, so
        # its h = 0 seed differs from the bare paramagnet
        src = os.path.dirname(_bond_model(self.dest()))
        work = os.path.join(self.tmp, "WORK")
        rc, out, err = self.main("prepare", src, "--out", work, "--tier", "small",
                                 "--para-check", "--max-candidates", "2")
        self.assertEqual(rc, 0, err)
        with open(os.path.join(work, "manifest.json")) as f:
            man = json.load(f)
        para = man["para"]
        self.assertEqual(para["status"], "converged")
        self.assertEqual(para["green"], "para_hf/green.npz")
        self.assertTrue(os.path.isfile(os.path.join(work, "para_hf", "green.npz")))
        self.assertTrue(os.path.isdir(os.path.join(work, "para_hf", "check")))
        with open(os.path.join(work, "coverage.json")) as f:
            cov = json.load(f)
        self.assertEqual(cov["tier"]["name"], "small")
        self.assertEqual(cov["tier"]["chosen_by"], "user")
        self.assertEqual(cov["hint"]["status"], "computed")
        tried = [r for r in cov["rows"] if r["tried"]]
        self.assertTrue(tried)
        self.assertEqual({r["seed_base"] for r in tried}, {"para_hf"})
        para_row, = [r for r in tried if r["family"] == "para"]
        self.assertEqual(para_row["label"], "para")
        # against the bare paramagnet the same seed carries the bond
        # renormalisation
        from hwave.initial import cli
        from hwave.initial.facts import model_facts
        import numpy as _np
        run_toml = os.path.join(work, para_row["run_id"], "input.toml")
        f = model_facts(run_toml)
        with _np.load(os.path.join(work, para_row["run_id"], "initial.npz")) as d:
            G = _np.array(d["green_sublattice"])
        self.assertIn("bond-charge", cli._state_label(G, f, {}))
        # run the para candidate only, move WORK and delete the sources, then
        # collect: the moved para reference gives the final label para
        run_solver(os.path.join(work, para_row["run_id"], "input.toml"))
        moved = os.path.join(self.tmp, "moved", "WORK")
        os.makedirs(os.path.dirname(moved))
        shutil.move(work, moved)
        shutil.rmtree(src)
        work = moved
        rc, out, err = self.main("collect", work)
        self.assertEqual(rc, 0, err)
        with open(os.path.join(work, "coverage.json")) as f:
            post = json.load(f)
        row, = [r for r in post["rows"] if r["run_id"] == para_row["run_id"]]
        self.assertEqual(row["status"], "converged")
        self.assertEqual(row["final_label"], "para")
        # the para reference is part of the checked state of WORK: a changed
        # green or record, or a manifest without the para entry, is refused
        man_path = os.path.join(work, "manifest.json")
        for name in ("green.npz", "para_hf.json", None):
            with self.subTest(changed=name):
                keep = {}
                targets = [man_path] if name is None else [
                    os.path.join(work, "para_hf", name)]
                for t in targets:
                    with open(t, "rb") as f:
                        keep[t] = f.read()
                if name is None:
                    with open(man_path) as f:
                        m = json.load(f)
                    m["para"] = None
                    with open(man_path, "w") as f:
                        json.dump(m, f)
                else:
                    with open(targets[0], "ab") as f:
                        f.write(b"\0" if name.endswith(".npz") else b" ")
                rc, out, err = self.main("collect", work)
                self.assertEqual(rc, 1)
                self.assertIn("para_hf/" + name if name else "no para reference", err)
                for t, data in keep.items():
                    with open(t, "wb") as f:
                        f.write(data)
        rc, out, err = self.main("collect", work)
        self.assertEqual(rc, 0, err)

    def test_para_check_with_large_refused(self):
        src = os.path.dirname(initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2,
                                                      subshape=(2, 2, 1)))
        work = os.path.join(self.tmp, "W")
        rc, _, err = self.main("prepare", src, "--out", work, "--tier", "large",
                               "--para-check")
        self.assertEqual(rc, 1)
        self.assertIn("--para-check", err)
        self.assertFalse(os.path.exists(work))
        rc, _, err = self.main("prepare", src, "--out", work, "--tier", "medium")
        self.assertEqual(rc, 1)
        self.assertIn("later release", err)

    def test_large_tier_has_no_para_entry(self):
        src = os.path.dirname(initial_fixtures.square(self.dest(), 4, U=2.0, T=0.2,
                                                      subshape=(2, 2, 1)))
        work = os.path.join(self.tmp, "W")
        rc, _, err = self.main("prepare", src, "--out", work, "--tier", "large",
                               "--max-candidates", "0", "--subshape", "2,2,1")
        self.assertEqual(rc, 0, err)
        with open(os.path.join(work, "manifest.json")) as f:
            self.assertIsNone(json.load(f)["para"])
        self.assertFalse(os.path.exists(os.path.join(work, "para_hf")))


if __name__ == "__main__":
    unittest.main()
