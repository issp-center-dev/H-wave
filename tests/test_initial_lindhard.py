"""The static Lindhard response on a pair space: ``static_lindhard`` and
``apply_lindhard`` of the solver facade, and the pair-space wrappers
``hints.lindhard`` / ``hints.lindhard_apply`` (sectors).

Three independent checks pin the result:

* the label-form oracle ``tests/so_bond_oracles.lindhard_static`` (its own
  ``eigh`` and Fermi function, one mu), mapped onto the pair-space labels
  by the relation derived in :func:`_oracle_on_labels`;
* a real-space linear response on the finite lattice: the exact density
  change of a q-modulated one-body field, from the eigenbasis of the whole
  lattice Hamiltonian (no k space). This pins the label coordinates, the
  anchoring of the phase and the sign ``delta rho = chi0 (-delta H)``;
* the V2b triangular values (chi0 per spin, mesh mean) against the recorded
  values and against the analytic dispersion.
"""
import os
import shutil
import tempfile
import tracemalloc
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures
from tests import so_bond_oracles as oracles


def _facts(toml, subshape=None):
    from hwave.initial.facts import model_facts
    return model_facts(toml, subshape=subshape)


def _failed_su2():
    from hwave.initial.symmetry import Su2Result
    return Su2Result(status="failed", max_residual=1.0, method=None, sz_residual=None)


def _hermitian_shift(p, seed, scale=0.05):
    """``H0(k) + Sigma'[G]`` for a random block-projected Hermitian G: a
    generic Hamiltonian in the model's blocks."""
    from hwave.initial import symmetry
    from hwave.solver import initial_facade as fa
    rng = np.random.default_rng(seed)
    gs = (p.nvol, p.ns, p.norb, p.ns, p.norb)
    g = rng.standard_normal(gs) + 1j * rng.standard_normal(gs)
    g = symmetry.hermitian_green(g, p.shape).reshape(p.nvol, p.nd, p.nd)
    blk = np.full(p.nd, -1)
    for k, b in enumerate(p.blocks.block_info):
        blk[np.asarray(b)] = k
    g[:, blk[:, None] != blk[None, :]] = 0
    g = scale * g.reshape(gs)
    return np.asarray(p.ham_trans_k) + fa.hf_potential(
        g, fa.hf_tables(p), include_fock=p.iflag_fock, blocks=p.blocks)


def _bands(p, H, T=None):
    from hwave.solver import initial_facade as fa
    return fa.bands_of(H, p.blocks, p.T if T is None else T, shape=p.shape, ns=p.ns,
                       norb=p.norb, ene_cutoff=p.ene_cutoff)


def _nd(p, a, s):
    from hwave.initial import hints
    return np.asarray(hints._nd_index(p, np.asarray(a), np.asarray(s)), dtype=np.int64)


def _qindex(q, shape):
    return int(oracles.mesh_index(np.array([int(Fraction(x) * n) for x, n in zip(q, shape)]),
                                  shape))


def _eiq(q, r):
    t = sum((Fraction(x) * int(y) for x, y in zip(q, r)), Fraction(0)) % 1
    return np.exp(2j * np.pi * float(t))


class _Geom:
    def __init__(self, delta_r, shape):
        self.delta_r = np.asarray(delta_r, dtype=np.int64)
        self.spatial_shape = tuple(shape)


def _oracle_on_labels(p, labels, H, mu, T, qs):
    """``lindhard_static`` mapped onto the labels.

    The oracle is ``chi[(d; X, Y), (d'; Z, W)] = (1/N) sum_k e^{ik.d}
    e^{i(k+q).d'} sum U_Ya(k) U*_Za(k) U_Wb(k+q) U*_Xb(k+q) F``. The pair-space
    response of the label ``m = (alpha, beta, r)`` to ``n = (alpha', beta',
    r')`` is ``(1/N) sum_k e^{i(k+q).(r - r')} sum u_alpha,l(k+q)
    u*_alpha',l(k+q) u_beta',m(k) u*_beta,m(k) F`` (l a band at k+q, m at k), so
    ``chi0[m, n] = e^{-iq.r'} chi[(-r'; alpha', beta'), (r; beta, alpha)]``.
    """
    r = labels[:, 4:7]
    D = np.unique(np.concatenate([r, -r]), axis=0)
    geom = _Geom(D, p.shape)
    nd = p.nd
    chi = oracles.lindhard_static(H, mu, 1.0 / T, geom)
    pos = {tuple(int(v) for v in d): k for k, d in enumerate(D)}
    al = _nd(p, labels[:, 0], labels[:, 1])
    be = _nd(p, labels[:, 2], labels[:, 3])
    row = np.array([pos[tuple(int(v) for v in -x)] for x in r]) * nd * nd + al * nd + be
    col = np.array([pos[tuple(int(v) for v in x)] for x in r]) * nd * nd + be * nd + al
    out = {}
    for q in qs:
        c = chi[_qindex(q, p.shape)]
        ph = np.array([_eiq(q, x) for x in r])
        out[q] = c[np.ix_(row, col)].T * np.conj(ph)[None, :]
    return out


class _Base(unittest.TestCase):
    def setUp(self):
        import logging
        self.tmp = tempfile.mkdtemp()
        self._n = 0
        logging.getLogger("qlms").setLevel(logging.WARNING)
        self.addCleanup(logging.getLogger("qlms").setLevel, logging.NOTSET)

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def dest(self):
        self._n += 1
        return os.path.join(self.tmp, "case{}".format(self._n))

    def golden(self, **kw):
        return initial_fixtures.golden(self.dest(), **kw)

    def model(self, toml, seed=3):
        from hwave.initial import hints
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f, su2=_failed_su2())
        H = _hermitian_shift(p, seed)
        return f, p, P, H


QS = ((0, 0, 0), (Fraction(1, 3), 0, 0), (Fraction(1, 2), Fraction(1, 3), 0),
      (Fraction(5, 6), Fraction(2, 3), 0))


class TestOracle(_Base):

    def test_vs_oracle(self):
        """Normal mode (one mu), fixed 2Sz (two mu: the down block is shifted
        by mu_up - mu_dn for the one-mu oracle, which leaves the eigenvectors
        and xi unchanged) and spin-orbital mode, 1e-10."""
        from hwave.solver import initial_facade as fa
        cases = {
            "normal": dict(cellshape=(6, 3, 1)),
            "2Sz": dict(cellshape=(6, 3, 1), two_sz=2, types=("CoulombIntra", "CoulombInter",
                                                              "Ising")),
            "spin-orbital": dict(cellshape=(6, 3, 1), spin_orbital=True),
        }
        for name, kw in cases.items():
            with self.subTest(name):
                f, p, P, H = self.model(self.golden(**kw))
                b = _bands(p, H)
                Ho = np.array(H)
                if name == "2Sz":
                    self.assertEqual(len(b.mu), 2)
                    self.assertGreater(abs(b.mu[0] - b.mu[1]), 1e-3)
                    for blk, g in zip(p.blocks.block_info, p.blocks.block_to_group):
                        idx = np.asarray(blk)
                        Ho[:, idx, idx] += b.mu[0] - b.mu[g]
                else:
                    self.assertEqual(len(b.mu), 1)
                want = _oracle_on_labels(p, P.labels, Ho, b.mu[0], b.T, QS)
                for q in QS:
                    got = fa.static_lindhard(b, q, P.labels)
                    scale = max(1.0, float(np.max(np.abs(want[q]))))
                    self.assertLessEqual(float(np.max(np.abs(got - want[q]))), 1e-10 * scale,
                                         msg=str(q))
                    self.assertGreater(float(np.max(np.abs(got))), 1e-3)


class TestRealSpace(_Base):

    def _lattice(self, p, H):
        """The whole-lattice one-body matrix ``h[(R1, i), (R2, j)]`` of ``H =
        sum c^dag_i(R1) h c_j(R2)``, ``h = (1/N) sum_k e^{ik.(R1 - R2)}
        H(k)``, with cells in C order."""
        N = p.nvol
        R = oracles.mesh_grid(p.shape)
        kv = oracles.mesh_momenta(p.shape)
        ph = np.exp(1j * (R @ kv.T))                  # (R, k): e^{ik.R}
        h = np.einsum("ak,kij,bk->aibj", ph, np.asarray(H), np.conj(ph)) / N
        return h.reshape(N * p.nd, N * p.nd), R

    def test_real_space_response(self):
        """``x = -chi0 h`` for a q-modulated field ``delta H = sum_R sum_n h_n
        e^{iq.R} c^dag_{alpha_n}(R + r_n) c_{beta_n}(R)``, with ``x_m = (1/N)
        sum_R e^{-iq.R} delta <c^dag_{beta_m}(R) c_{alpha_m}(R + r_m)>`` from the
        exact real-space response; the lattice matrix is first checked against
        the solver's Green function."""
        from hwave.solver import initial_facade as fa
        for name, kw in {"normal": {}, "spin-orbital": {"spin_orbital": True}}.items():
            with self.subTest(name):
                f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1), **kw), seed=5)
                b = _bands(p, H)
                self.assertEqual(len(b.mu), 1)
                h, R = self._lattice(p, H)
                N, nd = p.nvol, p.nd
                e, U = np.linalg.eigh(h)
                fe = 1.0 / (1.0 + np.exp((e - b.mu[0]) / b.T))
                # equilibrium: G[r]_ab = <c^dag_a(r) c_b(0)> = rho0[(0, b), (r, a)]
                rho0 = (U * fe) @ U.conj().T
                G = fa.green_function(H, p.T, p.blocks, shape=p.shape, ns=p.ns, norb=p.norb,
                                      ene_cutoff=p.ene_cutoff).reshape(N, nd, nd)
                G_rs = rho0.reshape(N, nd, N, nd)[0].transpose(1, 2, 0)   # (r, a, b)
                np.testing.assert_allclose(G_rs, G, atol=1e-12)
                cell = {tuple(int(v) for v in x): k for k, x in enumerate(R)}

                def at(Rv):
                    return cell[tuple(int(Rv[i]) % p.shape[i] for i in range(3))]

                al = _nd(p, P.labels[:, 0], P.labels[:, 1])
                be = _nd(p, P.labels[:, 2], P.labels[:, 3])
                rr = P.labels[:, 4:7]
                F = oracles._pair_kernel(e[:, None] - b.mu[0], e[None, :] - b.mu[0],
                                         1.0 / b.T)
                rng = np.random.default_rng(17)
                for q in ((Fraction(1, 3), Fraction(1, 3), 0), (Fraction(1, 2), Fraction(2, 3), 0),
                          (0, 0, 0)):
                    hv = rng.standard_normal(P.size) + 1j * rng.standard_normal(P.size)
                    V = np.zeros((N * nd, N * nd), dtype=np.complex128)
                    for n in range(P.size):
                        for c, Rv in enumerate(R):
                            V[at(Rv + rr[n]) * nd + al[n], c * nd + be[n]] += hv[n] * _eiq(q, Rv)
                    drho = -U @ ((U.conj().T @ V @ U) * F) @ U.conj().T
                    x = np.zeros(P.size, dtype=np.complex128)
                    for m in range(P.size):
                        for c, Rv in enumerate(R):
                            x[m] += np.conj(_eiq(q, Rv)) * drho[at(Rv + rr[m]) * nd + al[m],
                                                                c * nd + be[m]]
                    x /= N
                    chi = fa.static_lindhard(b, q, P.labels)
                    scale = max(1.0, float(np.max(np.abs(x))))
                    self.assertLessEqual(float(np.max(np.abs(x + chi @ hv))), 1e-12 * scale,
                                         msg=str(q))
                    self.assertLessEqual(
                        float(np.max(np.abs(x + fa.apply_lindhard(b, q, P.labels, hv)))),
                        1e-12 * scale)


class TestProperties(_Base):

    def test_psd_and_hermitian(self):
        from hwave.solver import initial_facade as fa
        for name, kw in {"normal": {}, "spin-orbital": {"spin_orbital": True},
                         "2Sz": {"two_sz": 2, "types": ("CoulombIntra", "Hund")}}.items():
            with self.subTest(name):
                f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1), **kw))
                b = _bands(p, H)
                for q in ((0, 0, 0), (Fraction(1, 2), Fraction(1, 3), 0),
                          (Fraction(1, 6), Fraction(2, 3), 0)):
                    chi = fa.static_lindhard(b, q, P.labels)
                    self.assertEqual(chi.shape, (P.size, P.size))
                    self.assertEqual(chi.dtype, np.complex128)
                    np.testing.assert_array_equal(chi, chi.conj().T)
                    w = np.linalg.eigvalsh(chi)
                    self.assertGreater(w[-1], 0.0)
                    self.assertGreaterEqual(w[0], -1e-13 * w[-1])

    def test_apply_equals_matrix(self):
        from hwave.solver import initial_facade as fa
        for name, kw in {"normal": {}, "spin-orbital": {"spin_orbital": True}}.items():
            with self.subTest(name):
                f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1), **kw))
                b = _bands(p, H)
                rng = np.random.default_rng(23)
                for q in ((0, 0, 0), (Fraction(5, 6), Fraction(1, 3), 0)):
                    chi = fa.static_lindhard(b, q, P.labels)
                    v = rng.standard_normal(P.size) + 1j * rng.standard_normal(P.size)
                    got = fa.apply_lindhard(b, q, P.labels, v)
                    self.assertEqual(got.shape, (P.size,))
                    want = chi @ v
                    self.assertLessEqual(float(np.max(np.abs(got - want))),
                                         1e-12 * float(np.max(np.abs(want))))
                    V = rng.standard_normal((P.size, 3)) + 1j * rng.standard_normal((P.size, 3))
                    got = fa.apply_lindhard(b, q, P.labels, V)
                    self.assertEqual(got.shape, (P.size, 3))
                    want = chi @ V
                    self.assertLessEqual(float(np.max(np.abs(got - want))),
                                         1e-12 * float(np.max(np.abs(want))))

    def test_r_representative_free(self):
        """``e^{i(k+q).r}`` is periodic on the mesh: a label's r shifted by a
        mesh period gives the same response."""
        from hwave.solver import initial_facade as fa
        f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1)))
        b = _bands(p, H)
        lab = np.array(P.labels)
        lab[:, 4] += 6
        lab[:, 5] -= 3
        q = (Fraction(1, 6), Fraction(1, 3), 0)
        np.testing.assert_allclose(fa.static_lindhard(b, q, lab),
                                   fa.static_lindhard(b, q, P.labels), atol=1e-14)

    def test_apply_does_not_form_dense(self):
        """The matrix-free product never holds a ``|P|^2`` array: 2048
        synthetic labels on a 16x16 mesh, peak traced memory below a quarter of
        the dense matrix."""
        from hwave.solver import initial_facade as fa
        p = _facts(initial_fixtures.square(self.dest(), 16, U=1.0, T=0.2,
                                           subshape=(1, 1, 1))).problem
        b = _bands(p, np.asarray(p.ham_trans_k))
        ii = np.arange(2048)
        r = ii // 4
        lab = np.column_stack([np.zeros_like(ii), (ii // 2) % 2, np.zeros_like(ii), ii % 2,
                               r % 16, r // 16, np.zeros_like(ii)]).astype(np.int64)
        self.assertEqual(len({tuple(x) for x in lab}), 2048)
        v = np.random.default_rng(1).standard_normal(2048).astype(np.complex128)
        q = (Fraction(1, 4), Fraction(1, 2), 0)
        tracemalloc.start()
        try:
            fa.apply_lindhard(b, q, lab, v)
            _, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        self.assertLess(peak, 0.25 * 16 * 2048 ** 2)


def _tri_eps(k):
    """The triangular dispersion of the fixture (bonds (1,0), (0,1), (-1,1)
    and negatives, coefficient -1)."""
    k1, k2 = k[:, 0], k[:, 1]
    return -2.0 * (np.cos(k1) + np.cos(k2) + np.cos(k2 - k1))


class TestV2b(_Base):

    def test_v2b_values(self):
        """12x12 triangular, Ncond = 88, T = 0.1, U = 2, one orbital: the
        para-HF bands are the bare bands shifted by U n/2, and chi0 per spin is
        the mesh mean of F(eps(k+q), eps(k)). Probe values 1e-4; the analytic
        dispersion 1e-10."""
        from scipy import optimize
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        toml = initial_fixtures.triangular(self.dest(), 12, U=2.0, T=0.1, ncond=88,
                                           subshape=(1, 1, 1))
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f)
        n_sigma = 88 / 144 / 2
        H = np.asarray(p.ham_trans_k) + 2.0 * n_sigma * np.eye(p.nd)
        b = _bands(p, H)
        up = int(np.nonzero((P.labels == [0, 0, 0, 0, 0, 0, 0]).all(axis=1))[0][0])
        kv = oracles.mesh_momenta(p.shape)
        eps = _tri_eps(kv)
        mu = optimize.brentq(lambda m: np.sum(1.0 / (1.0 + np.exp((eps - m) / 0.1))) - 44.0,
                             -7.0, 4.0, xtol=1e-15)
        add = oracles.mesh_add(p.shape)
        recorded = {"K": ((Fraction(2, 3), Fraction(1, 3), 0), 0.2269),
                 "Gamma": ((0, 0, 0), 0.2023),
                 "M": ((Fraction(1, 2), 0, 0), 0.1699)}
        for name, (q, value) in recorded.items():
            with self.subTest(name):
                chi = fa.static_lindhard(b, q, P.labels)
                got = float(chi[up, up].real)
                self.assertAlmostEqual(got, value, delta=1e-4)
                kq = add[_qindex(q, p.shape)]
                F = oracles._pair_kernel(eps[kq] - mu, eps - mu, 10.0)
                self.assertAlmostEqual(got, float(np.mean(F)), delta=1e-10)
        # K is the maximum of chi0 per spin over the mesh (the V2b claim)
        vals = {}
        for m1 in range(12):
            for m2 in range(12):
                q = (Fraction(m1, 12), Fraction(m2, 12), 0)
                vals[(m1, m2)] = float(fa.apply_lindhard(
                    b, q, P.labels, np.eye(P.size)[:, up])[up].real)
        best = max(vals.values())
        self.assertAlmostEqual(vals[(8, 4)], best, delta=1e-12)


class TestDegenerateAndZeroT(_Base):

    def test_degenerate_limit(self):
        """q = 0 on one orbital: only the degenerate pairs l = m at the same k
        contribute, each with -f'(eps - mu) = f (1 - f) / T."""
        from hwave.solver import initial_facade as fa
        from hwave.initial import hints
        f = _facts(initial_fixtures.chain(self.dest(), 8, U=1.0, T=0.2, ncond=6,
                                          subshape=(1, 1, 1)))
        p = f.problem
        b = _bands(p, np.asarray(p.ham_trans_k))
        P = hints.pair_space(f)
        chi = fa.static_lindhard(b, (0, 0, 0), P.labels)
        eps = -2.0 * np.cos(2.0 * np.pi * np.arange(8) / 8)
        fm = 1.0 / (1.0 + np.exp((eps - b.mu[0]) / 0.2))
        up = int(np.nonzero((P.labels == [0, 0, 0, 0, 0, 0, 0]).all(axis=1))[0][0])
        self.assertAlmostEqual(float(chi[up, up].real), float(np.mean(fm * (1 - fm) / 0.2)),
                               delta=1e-13)

    def test_kernel_continuous_at_tolerance(self):
        """The difference quotient just outside the degeneracy tolerance
        agrees with the derivative used inside it."""
        from hwave.solver import initial_facade as fa
        T, mu = 0.1, 0.0
        for x in (-0.15, 0.0, 0.07):
            f = 1.0 / (1.0 + np.exp((x - mu) / T))
            for d in (3e-11, 2e-10, 1e-8):
                y = x + d
                fy = 1.0 / (1.0 + np.exp((y - mu) / T))
                k = fa._lindhard_kernel(np.array([x]), np.array([f]), np.array([y]),
                                        np.array([fy]), T)
                self.assertAlmostEqual(float(k[0]), f * (1 - f) / T, delta=1e-5)

    def test_zero_t_closed_shell(self):
        """T = 0 with a closed shell: degenerate pairs have equal occupations
        and contribute 0 (the T -> 0 limit of -f'); the result equals that
        of a small T. With a fixed 2Sz the up and down blocks are filled to
        different levels: a filled up level and an empty down level of the
        same energy belong to different blocks and do not pair."""
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        cases = {"Sz free": dict(ncond=6), "2Sz": dict(ncond=4, param={"2Sz": 2})}
        for name, kw in cases.items():
            with self.subTest(name):
                toml = initial_fixtures.chain(self.dest(), 6, U=1.0, T=0.0,
                                              subshape=(1, 1, 1), **kw)
                f = _facts(toml)
                p = f.problem
                P = hints.pair_space(f)
                H = np.asarray(p.ham_trans_k)
                b0 = _bands(p, H, T=0.0)
                b1 = _bands(p, H, T=1e-3)
                self.assertEqual(len(b0.mu), 1 if name == "Sz free" else 2)
                for q in ((0, 0, 0), (Fraction(1, 3), 0, 0), (Fraction(1, 2), 0, 0),
                          (Fraction(2, 3), 0, 0)):
                    c0 = fa.static_lindhard(b0, q, P.labels)
                    c1 = fa.static_lindhard(b1, q, P.labels)
                    np.testing.assert_allclose(c0, c1, atol=1e-12)
                    np.testing.assert_allclose(
                        fa.apply_lindhard(b0, q, P.labels, np.ones(P.size)),
                        c0 @ np.ones(P.size), atol=1e-12)

    def test_zero_t_open_shell_refused(self):
        """T = 0 with an open shell: a degenerate pair in one block with
        different occupations has an unbounded response."""
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        toml = initial_fixtures.chain(self.dest(), 6, U=1.0, T=0.0, ncond=3,
                                      subshape=(1, 1, 1))
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f)
        b = _bands(p, np.asarray(p.ham_trans_k), T=0.0)
        for q in ((Fraction(1, 3), 0, 0), (Fraction(2, 3), 0, 0)):
            with self.assertRaisesRegex(ValueError, "open Fermi shell"):
                fa.static_lindhard(b, q, P.labels)
            with self.assertRaisesRegex(ValueError, "open Fermi shell"):
                fa.apply_lindhard(b, q, P.labels, np.ones(P.size))


    def test_zero_t_open_shell_off_labels(self):
        """T = 0 with an open shell only on an orbital outside the labels
        (orbital 1 does not interact and has its own block): the response on
        the labels is finite and equals that of a small T."""
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        toml = initial_fixtures.two_orbital(self.dest(), 6, delta=10.0, U=(1.0, 0.0),
                                            T=0.0, ncond=15, subshape=(1, 1, 1))
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f)
        self.assertTrue(np.all(P.labels[:, 0] == 0) and np.all(P.labels[:, 2] == 0))
        H = np.asarray(p.ham_trans_k)
        b0 = _bands(p, H, T=0.0)
        b1 = _bands(p, H, T=1e-3)
        lab1 = np.array(P.labels[P.labels[:, 1] == P.labels[:, 3]])
        lab1[:, 0] = lab1[:, 2] = 1
        for q in ((Fraction(1, 3), 0, 0), (Fraction(2, 3), 0, 0)):
            np.testing.assert_allclose(fa.static_lindhard(b0, q, P.labels),
                                       fa.static_lindhard(b1, q, P.labels), atol=1e-12)
            with self.assertRaisesRegex(ValueError, "open Fermi shell"):
                fa.static_lindhard(b0, q, lab1)

    def test_zero_t_nearly_open_shell(self):
        """At T = 0 a pair split by less than LINDHARD_ZERO_T_GAP with
        different occupations is refused (nan), not given 1/|split|; a wider
        split keeps the difference quotient."""
        from hwave.solver import initial_facade as fa
        for d in (2e-10, 1e-8, 0.5 * fa.LINDHARD_ZERO_T_GAP):
            k = fa._lindhard_kernel(np.array([0.0]), np.array([1.0]), np.array([d]),
                                    np.array([0.0]), 0.0)
            self.assertTrue(np.isnan(k[0]))
        k = fa._lindhard_kernel(np.array([0.0]), np.array([1.0]), np.array([1e-3]),
                                np.array([0.0]), 0.0)
        self.assertAlmostEqual(float(k[0]), 1e3, delta=1e-9)
        k = fa._lindhard_kernel(np.array([0.0]), np.array([1.0]), np.array([1e-8]),
                                np.array([1.0]), 0.0)
        self.assertEqual(float(k[0]), 0.0)


class TestArguments(_Base):

    def test_refusals(self):
        from hwave.solver import initial_facade as fa
        f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1)))
        b = _bands(p, H)
        with self.assertRaisesRegex(ValueError, "not on the k-mesh"):
            fa.static_lindhard(b, (Fraction(1, 5), 0, 0), P.labels)
        with self.assertRaisesRegex(ValueError, "not on the k-mesh"):
            fa.static_lindhard(b, (0.3, 0, 0), P.labels)
        with self.assertRaisesRegex(ValueError, r"\(\|P\|, 7\)"):
            fa.static_lindhard(b, (0, 0, 0), P.labels[:, :6])
        bad = np.array(P.labels)
        bad[0, 0] = 2
        with self.assertRaisesRegex(ValueError, "orbital"):
            fa.static_lindhard(b, (0, 0, 0), bad)
        bad = np.array(P.labels)
        bad[0, 1] = 2
        with self.assertRaisesRegex(ValueError, "spin"):
            fa.static_lindhard(b, (0, 0, 0), bad)
        with self.assertRaisesRegex(ValueError, "v must be"):
            fa.apply_lindhard(b, (0, 0, 0), P.labels, np.ones(P.size + 1))
        # floats on the mesh and strings are accepted
        np.testing.assert_allclose(fa.static_lindhard(b, (0.5, 1 / 3, 0), P.labels),
                                   fa.static_lindhard(b, ("1/2", "1/3", 0), P.labels),
                                   atol=1e-15)


    def test_empty_labels(self):
        from hwave.solver import initial_facade as fa
        p = _facts(initial_fixtures.chain(self.dest(), 6, T=0.1, subshape=(1, 1, 1))).problem
        b = _bands(p, np.asarray(p.ham_trans_k))
        lab = np.zeros((0, 7), dtype=np.int64)
        self.assertEqual(fa.static_lindhard(b, (0, 0, 0), lab).shape, (0, 0))
        self.assertEqual(fa.apply_lindhard(b, (0, 0, 0), lab, np.zeros(0)).shape, (0,))


class TestPairSpaceWrappers(_Base):

    def test_full_space(self):
        from hwave.initial import hints
        from hwave.solver import initial_facade as fa
        f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1)))
        b = _bands(p, H)
        q = (Fraction(1, 2), Fraction(1, 3), 0)
        np.testing.assert_array_equal(hints.lindhard(b, q, P),
                                      fa.static_lindhard(b, q, P.labels))
        v = np.arange(P.size) + 0.5j
        np.testing.assert_array_equal(hints.lindhard_apply(b, q, P, v),
                                      fa.apply_lindhard(b, q, P.labels, v))

    def test_sectors(self):
        """Under a passed spin-rotation check (para bands with SU(2)) chi0 does
        not couple the charge and S^z sectors, and the two sector blocks hold
        the spectrum of chi0 on the spin-diagonal labels."""
        from hwave.initial import hints
        toml = initial_fixtures.golden(self.dest(), cellshape=(6, 3, 1),
                                       types=("CoulombIntra", "CoulombInter"))
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f)
        self.assertEqual(P.sectors, (hints.SECTOR_CHARGE, hints.SECTOR_SZ))
        b = _bands(p, np.asarray(p.ham_trans_k))
        q = (Fraction(1, 6), Fraction(1, 3), 0)
        full = hints.lindhard(b, q, P)
        uc, uz = (u.toarray() for u in P.sector_basis)
        np.testing.assert_allclose(uc.conj().T @ full @ uz, 0, atol=1e-13)
        cc = hints.lindhard(b, q, P, sector=hints.SECTOR_CHARGE)
        cz = hints.lindhard(b, q, P, sector=hints.SECTOR_SZ)
        np.testing.assert_allclose(cc, uc.conj().T @ full @ uc, atol=1e-14)
        diag = np.nonzero(P.labels[:, 1] == P.labels[:, 3])[0]
        want = np.linalg.eigvalsh(full[np.ix_(diag, diag)])
        got = np.sort(np.concatenate([np.linalg.eigvalsh(cc), np.linalg.eigvalsh(cz)]))
        np.testing.assert_allclose(got, want, atol=1e-13)
        v = np.random.default_rng(2).standard_normal(cz.shape[0]) + 0j
        np.testing.assert_allclose(hints.lindhard_apply(b, q, P, v, sector=hints.SECTOR_SZ),
                                   cz @ v, atol=1e-13)

    def test_foreign_bands_refused(self):
        from hwave.initial import hints
        f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1)))
        p2 = _facts(self.golden(cellshape=(4, 3, 1), types=("CoulombIntra",))).problem
        b2 = _bands(p2, np.asarray(p2.ham_trans_k))
        with self.assertRaisesRegex(ValueError, "bands"):
            hints.lindhard(b2, (0, 0, 0), P)
        so = _facts(self.golden(cellshape=(6, 3, 1), spin_orbital=True)).problem
        b3 = _bands(so, np.asarray(so.ham_trans_k))
        with self.assertRaisesRegex(ValueError, "bands"):
            hints.lindhard(b3, (0, 0, 0), P)

    def test_other_block_structure_refused(self):
        """Bands of a fixed 2Sz (spin blocks) against a pair space built
        without it (spin-flip labels): refused, not silently zero."""
        from hwave.initial import hints
        types = ("CoulombIntra", "CoulombInter")
        f, p, P, H = self.model(self.golden(cellshape=(6, 3, 1), types=types))
        self.assertTrue(np.any(P.labels[:, 1] != P.labels[:, 3]))
        p2 = _facts(self.golden(cellshape=(6, 3, 1), types=types, two_sz=2)).problem
        b2 = _bands(p2, np.asarray(p2.ham_trans_k))
        for fn in (lambda: hints.lindhard(b2, (0, 0, 0), P),
                   lambda: hints.lindhard_apply(b2, (0, 0, 0), P, np.ones(P.size))):
            with self.assertRaisesRegex(ValueError, "block structure"):
                fn()


if __name__ == "__main__":
    unittest.main()
