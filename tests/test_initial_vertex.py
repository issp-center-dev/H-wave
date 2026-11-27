"""The Hartree-Fock vertex of ``hwave.initial``: ``pair_space``, ``vertex``,
``vertex_sparse`` and ``vertex_self_adjoint``.

The oracle is the solver's own map. At q = 0 the vertex applied to the
rho-components of a random Hermitian Green function must give the
potential ``hf_potential`` returns, slot by slot, and nothing outside the
pair space. At q != 0 the model is loaded on a SubShape that makes q
commensurate, a q-modulated rho is written into the folded Green function,
and the complex-linear extension of ``hf_potential`` (its value on the
Hermitian part plus i times its value on the anti-Hermitian part) must be
``(Gamma(q) x)_m e^{iq.R}`` at every anchor R. This pins the label
coordinates, the anchoring of the phases and the PairHop linear read
together.
"""
import itertools
import shutil
import tempfile
import tracemalloc
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures

TYPES = tuple(initial_fixtures.GOLDEN_TABLES)


def _facts(toml, subshape=None):
    from hwave.initial.facts import model_facts
    return model_facts(toml, subshape=subshape)


def _failed_su2():
    from hwave.initial.symmetry import Su2Result
    return Su2Result(status="failed", max_residual=1.0, method=None, sz_residual=None)


def _passed_su2():
    from hwave.initial.symmetry import Su2Result
    return Su2Result(status="passed", max_residual=0.0, method=None, sz_residual=None)


def _block_projected_random_green(p, rng):
    from hwave.initial import symmetry
    gs = (p.nvol, p.ns, p.norb, p.ns, p.norb)
    g = rng.standard_normal(gs) + 1j * rng.standard_normal(gs)
    g = symmetry.hermitian_green(g, p.shape).reshape(p.nvol, p.nd, p.nd)
    blk = np.full(p.nd, -1)
    for k, b in enumerate(p.blocks.block_info):
        blk[np.asarray(b)] = k
    g[:, blk[:, None] != blk[None, :]] = 0
    return g.reshape(gs)


def _sigma_r(p, g):
    from hwave.solver import initial_facade
    sig = initial_facade.hf_potential(g, p.tables, include_fock=p.iflag_fock, blocks=p.blocks)
    return np.fft.fftn(sig.reshape(*p.shape, p.nd, p.nd), axes=(0, 1, 2),
                       norm="forward").reshape(p.nvol, p.nd, p.nd)


def _slot_indices(p, labels):
    """``(r index of -r, i, j)`` of each label: the Hamiltonian-convention
    potential of label m is ``Sigma_H(-r)[i, j]``; its rho-component is
    ``G(-r)[j, i]``."""
    from hwave.initial import hints
    neg = hints._neg_r(p.shape)
    rid = np.ravel_multi_index(tuple(labels[:, 4 + i] % p.shape[i] for i in range(3)),
                               p.shape)
    i = np.asarray(hints._nd_index(p, labels[:, 0], labels[:, 1]))
    j = np.asarray(hints._nd_index(p, labels[:, 2], labels[:, 3]))
    return neg[rid], i, j


def _dagger(g, shape):
    """``g^dag`` of a stored array ``(nvol, nd, nd)``: ``conj(g(-r))^T``."""
    from hwave.initial import hints
    return np.conj(np.swapaxes(g[hints._neg_r(shape)], -1, -2))


def _qdot(q, r):
    return sum((Fraction(x) * int(y) for x, y in zip(q, r)), Fraction(0)) % 1


def _eiq(q, r):
    return np.exp(2j * np.pi * float(_qdot(q, r)))


class _Base(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self._n = 0

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def golden(self, **kw):
        self._n += 1
        return initial_fixtures.golden("{}/m{}".format(self.tmp, self._n), **kw)

    def assert_identity(self, toml, *, tol=1e-13):
        """q = 0: Gamma x equals hf_potential slot by slot, on random
        Hermitian G and (complex-linear extension) on a random non-Hermitian
        rho."""
        from hwave.initial import hints
        f = _facts(toml)
        p = f.problem
        P = hints.pair_space(f, su2=_failed_su2())
        gamma = hints.vertex(f, (0, 0, 0), P)
        ri, i, j = _slot_indices(p, P.labels)
        rng = np.random.default_rng(7)
        for _ in range(2):
            g = _block_projected_random_green(p, rng)
            gm = g.reshape(p.nvol, p.nd, p.nd)
            sig = _sigma_r(p, g)
            x = gm[ri, j, i]
            want = sig[ri, i, j]
            got = gamma @ x
            scale = max(1.0, float(np.max(np.abs(want))))
            self.assertLessEqual(float(np.max(np.abs(got - want))), tol * scale)
            rest = sig.copy()
            rest[ri, i, j] = 0
            self.assertLessEqual(float(np.max(np.abs(rest))), tol * scale)
        # complex-linear extension: rho = H1 + i H2
        g1 = _block_projected_random_green(p, rng)
        g2 = _block_projected_random_green(p, rng)
        x = g1.reshape(p.nvol, p.nd, p.nd)[ri, j, i] + 1j * g2.reshape(p.nvol, p.nd, p.nd)[ri, j, i]
        want = (_sigma_r(p, g1) + 1j * _sigma_r(p, g2))[ri, i, j]
        self.assertLessEqual(float(np.max(np.abs(gamma @ x - want))),
                             tol * max(1.0, float(np.max(np.abs(want)))))
        return P, gamma

    def assert_supercell(self, toml, q, B, *, tol=1e-12):
        """q != 0 through the solver map of the folded cell ``B``."""
        from hwave.initial import hints
        from hwave.solver import initial_facade
        f = _facts(toml)
        p1 = f.problem
        P = hints.pair_space(f, su2=_failed_su2())
        gamma = hints.vertex(f, q, P)
        pB = initial_facade.load_problem(toml, subshape=B)
        nphys, norb_orig = p1.norb_phys, p1.norb_orig
        Nsup = pB.shape
        so = pB.enable_spin_orbital

        def ir(o):
            return o[0] + B[0] * (o[1] + B[1] * o[2])

        def idx(s, a, o):
            if so:
                return 2 * (a + nphys * ir(o)) + s
            return s * pB.norb + a + norb_orig * ir(o)

        offsets = list(itertools.product(range(B[0]), range(B[1]), range(B[2])))
        rng = np.random.default_rng(11)
        x = rng.standard_normal(P.size) + 1j * rng.standard_normal(P.size)
        dG = np.zeros(tuple(Nsup) + (pB.nd, pB.nd), dtype=np.complex128)
        placed = []
        for n, (a, s, b, t, *r) in enumerate(P.labels):
            for oa in offsets:
                R = [oa[k] - r[k] for k in range(3)]
                D = tuple((R[k] // B[k]) % Nsup[k] for k in range(3))
                oR = tuple(R[k] % B[k] for k in range(3))
                ph = _eiq(q, R)
                i, jj = idx(s, a, oa), idx(t, b, oR)
                self.assertEqual(dG[D + (jj, i)], 0)
                dG[D + (jj, i)] = x[n] * ph
                placed.append((n, D, i, jj, ph))
        dG = dG.reshape(pB.nvol, pB.nd, pB.nd)
        dd = _dagger(dG, pB.shape)
        gs = (pB.nvol, pB.ns, pB.norb, pB.ns, pB.norb)
        h1 = ((dG + dd) / 2).reshape(gs)
        h2 = ((dG - dd) / 2j).reshape(gs)
        sig = (_sigma_r(pB, h1) + 1j * _sigma_r(pB, h2)).reshape(tuple(Nsup) + (pB.nd, pB.nd))
        y = gamma @ x
        scale = max(1.0, float(np.max(np.abs(y))))
        covered = np.zeros(sig.shape, dtype=bool)
        for n, D, i, jj, ph in placed:
            self.assertLessEqual(abs(sig[D + (i, jj)] - y[n] * ph), tol * scale,
                                 msg="label {} anchor phase {}".format(P.labels[n], ph))
            covered[D + (i, jj)] = True
        self.assertLessEqual(float(np.max(np.abs(sig[~covered]), initial=0.0)), tol * scale)
        return P, gamma


class TestPairHopExtension(_Base):
    """The complex-linear extension of the solver's map is exact."""

    def test_pairhop_extension_identity(self):
        for fock in (True, False):
            with self.subTest(fock=fock):
                self.assert_identity(self.golden(types=("PairHop",), fock=fock))
                self.assert_identity(self.golden(fock=fock))

    def test_identity_modes(self):
        for kw in ({"two_sz": 0}, {"spin_orbital": True}, {"spin_orbital": True, "fock": False},
                   {"inter_orbital": False}):
            with self.subTest(**kw):
                self.assert_identity(self.golden(**kw))

    def test_supercell_q_phases(self):
        cases = (
            ({}, (Fraction(1, 3), 0, 0), (3, 1, 1)),
            ({}, (Fraction(1, 2), Fraction(1, 3), 0), (2, 3, 1)),
            ({"fock": False}, (Fraction(2, 3), Fraction(1, 3), 0), (3, 3, 1)),
            ({"types": ("PairHop",)}, (Fraction(1, 3), Fraction(2, 3), 0), (3, 3, 1)),
            ({"types": ("CoulombInter",)}, (Fraction(1, 6), 0, 0), (6, 1, 1)),
            ({"spin_orbital": True}, (Fraction(1, 3), Fraction(1, 3), 0), (3, 3, 1)),
            ({"two_sz": 0}, (Fraction(1, 3), 0, 0), (3, 1, 1)),
        )
        for kw, q, B in cases:
            with self.subTest(kw=kw, q=q):
                self.assert_supercell(self.golden(**kw), q, B)


class TestPairSpace(_Base):

    def _label_set(self, P):
        return {tuple(int(v) for v in row) for row in P.labels}

    def test_fock_off_pair_space(self):
        from hwave.initial import hints
        tables = {"CoulombIntra": initial_fixtures.GOLDEN_TABLES["CoulombIntra"],
                  "PairHop": {((0, 0, 0), 0, 1): 0.2 + 0.1j}}
        onsite = {(a, s, a, s, 0, 0, 0) for a in range(2) for s in range(2)}
        pairhop = {(a, s, 1 - a, s, 0, 0, 0) for a in range(2) for s in range(2)}
        # Transfer connects the orbitals (without Fock only the spins are
        # separate blocks): the PairHop slots stay
        f = _facts(self.golden(tables=tables, fock=False))
        self.assertEqual(len(f.problem.blocks.block_info), 2)
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertEqual(self._label_set(P), onsite | pairhop)
        # no inter-orbital Transfer: the blocks separate the orbitals and the
        # solver's mask removes the inter-orbital PairHop slots
        f = _facts(self.golden(tables=tables, fock=False, inter_orbital=False))
        self.assertEqual(len(f.problem.blocks.block_info), 4)
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertEqual(self._label_set(P), onsite)
        self.assertFalse(P.include_fock)

    def test_fock_on_adds_exchange_slots(self):
        from hwave.initial import hints
        bond = {((1, 0, 0), 0, 0): 0.3}
        f = _facts(self.golden(tables={"CoulombInter": bond}, fock=False))
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertEqual(self._label_set(P), {(0, s, 0, s, 0, 0, 0) for s in range(2)})
        f = _facts(self.golden(tables={"CoulombInter": bond}, fock=True))
        P = hints.pair_space(f, su2=_failed_su2())
        want = {(0, s, 0, s, 0, 0, 0) for s in range(2)}
        want |= {(0, s, 0, t, rx, 0, 0) for s in range(2) for t in range(2) for rx in (1, -1)}
        self.assertEqual(self._label_set(P), want)

    def test_two_sz_masks_transverse(self):
        from hwave.initial import hints
        f = _facts(self.golden(two_sz=0))
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertTrue(np.all(P.labels[:, 1] == P.labels[:, 3]))
        self.assertGreater(P.size, 0)

    def test_labels_layout(self):
        from hwave.initial import hints
        f = _facts(self.golden(cellshape=(6, 4, 1)))
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertEqual(P.labels.shape[1], 7)
        self.assertEqual(P.labels.dtype, np.int64)
        self.assertFalse(P.labels.flags.writeable)
        # centred representative: -N/2 < r <= N/2
        for i, n in enumerate(f.problem.shape):
            self.assertTrue(np.all(P.labels[:, 4 + i] > -(n / 2)))
            self.assertTrue(np.all(P.labels[:, 4 + i] <= n // 2))
        # the PairHop bond (2, 0, 0) (1, 0) and its partner are present
        s = self._label_set(P)
        self.assertIn((0, 0, 1, 0, 2, 0, 0), s)
        self.assertIn((1, 0, 0, 0, -2, 0, 0), s)
        # spin-orbital mode: the same physical labels
        f2 = _facts(self.golden(cellshape=(6, 4, 1), spin_orbital=True))
        P2 = hints.pair_space(f2, su2=_failed_su2())
        self.assertTrue(np.all(P2.labels[:, [0, 2]] < 2))
        self.assertEqual(P2.sectors, (hints.SECTOR_FULL,))

    def test_reverse_closed(self):
        from hwave.initial import hints
        for kw in ({}, {"fock": False}, {"two_sz": 0}, {"spin_orbital": True},
                   {"types": ("PairHop",)}, {"types": ("PairLift",)}):
            with self.subTest(**kw):
                f = _facts(self.golden(**kw))
                p = f.problem
                out, inp, coef, _rid, _kind = hints._kept_terms(p, f)
                self.assertTrue(np.all(coef != 0))
                self.assertEqual(set(out.tolist()), set(inp.tolist()))
                P = hints.pair_space(f, su2=_failed_su2())
                s = {tuple(int(v) for v in row) for row in P.labels}
                for (a, si, b, t, rx, ry, rz) in s:
                    nb = [(-v) for v in (rx, ry, rz)]
                    nb = [v - n if v > n // 2 else (v + n if v <= -(n / 2) else v)
                          for v, n in zip(nb, p.shape)]
                    self.assertIn((b, t, a, si, *nb), s)

    def test_interacting_orbitals_only(self):
        from hwave.initial import hints
        f = _facts(self.golden(tables={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0}}))
        self.assertEqual(f.interacting_orbitals, (0,))
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertTrue(np.all(P.labels[:, [0, 2]] == 0))
        self.assertEqual(P.size, 4)


class TestVertex(_Base):

    def _q_points(self, shape):
        qs = [(0, 0, 0), (Fraction(1, shape[0]), 0, 0)]
        if shape[1] > 1:
            qs.append((Fraction(1, shape[0]), Fraction(1, shape[1]), 0))
            qs.append((Fraction(shape[0] // 2, shape[0]), Fraction(shape[1] - 1, shape[1]), 0))
        return qs

    def test_gamma_hermitian(self):
        from hwave.initial import hints
        cases = [{"types": (name,)} for name in TYPES] + [{}]
        for base in cases:
            for extra in ({}, {"fock": False}, {"two_sz": 0}, {"spin_orbital": True}):
                kw = dict(base, **extra)
                with self.subTest(**kw):
                    f = _facts(self.golden(**kw))
                    P = hints.pair_space(f, su2=_failed_su2())
                    for q in self._q_points(f.problem.shape):
                        g = hints.vertex(f, q, P)
                        ok, res = hints.vertex_self_adjoint(g)
                        self.assertTrue(ok, msg="q={} residual {}".format(q, res))
                        ok, res = hints.vertex_self_adjoint(hints.vertex_sparse(f, q, P))
                        self.assertTrue(ok, msg="sparse q={} residual {}".format(q, res))

    def test_self_adjoint_detects_defect(self):
        from hwave.initial import hints
        import scipy.sparse
        g = np.array([[1.0, 2.0 + 1e-6j], [2.0, 0.5]])
        self.assertFalse(hints.vertex_self_adjoint(g)[0])
        self.assertFalse(hints.vertex_self_adjoint(scipy.sparse.csr_array(g))[0])
        h = np.array([[1.0, 2.0 + 1j], [2.0 - 1j, 0.5]])
        ok, res = hints.vertex_self_adjoint(h)
        self.assertTrue(ok)
        self.assertEqual(res, 0.0)

    def test_su2_sectors_equal_full(self):
        from hwave.initial import hints, symmetry
        L = 4
        tables = {"CoulombIntra": {((0, 0, 0), 0, 0): 2.0, ((0, 0, 0), 1, 1): 1.0},
                  "CoulombInter": {((1, 0, 0), 0, 0): 0.4, ((0, 0, 0), 0, 1): 0.6,
                                   ((0, 1, 0), 1, 0): 0.3},
                  "PairHop": {((0, 0, 0), 0, 1): 0.2}}
        f = _facts(self.golden(cellshape=(L, L, 1), tables=tables))
        su2 = symmetry.su2_check(f.problem)
        self.assertEqual(su2.status, "passed")
        P = hints.pair_space(f, su2=su2)
        self.assertEqual(P.sectors, (hints.SECTOR_CHARGE, hints.SECTOR_SZ))
        self.assertEqual(hints.pair_space(f).sectors, P.sectors)
        # every slot here carries the transverse labels too
        diag = P.labels[:, 1] == P.labels[:, 3]
        self.assertEqual(2 * P.sector_dim(hints.SECTOR_SZ), int(np.sum(~diag)))
        for k, name in enumerate(P.sectors):
            sl = P.sector_labels[k]
            self.assertEqual(sl.shape, (P.sector_dim(name), 6))
            self.assertTrue(np.all(sl[:, 1] == (0 if name == hints.SECTOR_CHARGE else 3)))
        uc, uz = (b.toarray() for b in P.sector_basis)
        for q in self._q_points(f.problem.shape):
            g = hints.vertex(f, q, P)
            gc = hints.vertex(f, q, P, sector=hints.SECTOR_CHARGE)
            gz = hints.vertex(f, q, P, sector=hints.SECTOR_SZ)
            np.testing.assert_allclose(gc, uc.conj().T @ g @ uc, rtol=0, atol=1e-14)
            # the sectors do not mix
            self.assertLessEqual(float(np.max(np.abs(uc.conj().T @ g @ uz))), 1e-12)
            full = np.sort(np.linalg.eigvalsh(g))
            parts = np.sort(np.concatenate([np.linalg.eigvalsh(gc)]
                                           + [np.linalg.eigvalsh(gz)] * 3))
            nz_full = full[np.abs(full) > 1e-9]
            nz_parts = parts[np.abs(parts) > 1e-9]
            self.assertEqual(nz_full.size, nz_parts.size)
            np.testing.assert_allclose(nz_full, nz_parts, rtol=0, atol=1e-12)

    def test_no_sectors_without_su2(self):
        from hwave.initial import hints
        for kw, su2 in (({}, None), ({"two_sz": 0}, None), ({"spin_orbital": True}, None)):
            with self.subTest(**kw):
                f = _facts(self.golden(**kw))
                P = hints.pair_space(f, su2=su2)
                self.assertEqual(P.sectors, (hints.SECTOR_FULL,))
                self.assertIsNone(P.sector_basis[0])
                self.assertEqual(P.sector_dim(hints.SECTOR_FULL), P.size)

    def test_sparse_equals_dense(self):
        from hwave.initial import hints, symmetry
        import scipy.sparse
        cases = [{"types": (name,)} for name in TYPES] + [{}, {"fock": False},
                                                          {"two_sz": 0}, {"spin_orbital": True}]
        for kw in cases:
            with self.subTest(**kw):
                f = _facts(self.golden(**kw))
                P = hints.pair_space(f, su2=_failed_su2())
                for q in self._q_points(f.problem.shape):
                    s = hints.vertex_sparse(f, q, P)
                    self.assertIsInstance(s, scipy.sparse.csr_array)
                    self.assertEqual(s.dtype, np.complex128)
                    self.assertTrue(s.has_sorted_indices)
                    self.assertTrue(s.has_canonical_format)
                    np.testing.assert_allclose(s.toarray(), hints.vertex(f, q, P),
                                               rtol=0, atol=1e-14)
        tables = {"CoulombIntra": {((0, 0, 0), 0, 0): 2.0},
                  "CoulombInter": {((1, 0, 0), 0, 0): 0.4, ((0, 1, 0), 0, 1): 0.3},
                  "PairHop": {((0, 0, 0), 0, 1): 0.2}}
        f = _facts(self.golden(cellshape=(4, 4, 1), tables=tables))
        su2 = symmetry.su2_check(f.problem)
        self.assertEqual(su2.status, "passed")
        P = hints.pair_space(f, su2=su2)
        nnz = None
        for q in self._q_points(f.problem.shape):
            for name in P.sectors:
                s = hints.vertex_sparse(f, q, P, sector=name)
                self.assertTrue(s.has_canonical_format)
                self.assertEqual(s.shape, (P.sector_dim(name),) * 2)
                np.testing.assert_allclose(s.toarray(), hints.vertex(f, q, P, sector=name),
                                           rtol=0, atol=1e-14)
            full = hints.vertex_sparse(f, q, P)
            if nnz is None:
                nnz = full.nnz
            self.assertEqual(full.nnz, nnz)  # structural: the same for every q

    def test_sparse_memory_is_nnz_bound(self):
        """At |P| ~ 2000 the sparse assembly allocates O(nnz), far below the
        dense |P|^2 array."""
        from hwave.initial import hints
        Lx, Ly = 25, 20
        bonds = {}
        for rx in range(-(Lx // 2) + 1, Lx // 2 + 1):
            for ry in range(-(Ly // 2) + 1, Ly // 2 + 1):
                if (rx, ry) != (0, 0) and (-rx, -ry, 0) not in [k[0] for k in bonds]:
                    bonds[((rx, ry, 0), 0, 0)] = 0.01
        toml = initial_fixtures.write_input(
            "{}/big".format(self.tmp), cellshape=(Lx, Ly, 1),
            transfer=initial_fixtures.square_transfer(),
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 1.0},
                          "CoulombInter": initial_fixtures.hermitian_closed(bonds)},
            param={"Ncond": Lx * Ly, "T": 0.1, "SubShape": [1, 1, 1],
                   "trustme_interaction_range": True})
        f = _facts(toml)
        P = hints.pair_space(f, su2=_failed_su2())
        self.assertGreaterEqual(P.size, 1900)
        q = (Fraction(1, Lx), Fraction(2, Ly), 0)
        hints.vertex_sparse(f, q, P)  # warm the imports
        tracemalloc.start()
        try:
            s = hints.vertex_sparse(f, q, P)
            _cur, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        dense = 16 * P.size ** 2
        msg = "peak {} B, |P| {}, nnz {}".format(peak, P.size, s.nnz)
        self.assertLess(peak, dense / 20, msg=msg)
        self.assertLess(peak, 400 * max(s.nnz, P.size) + 2_000_000, msg=msg)


class TestVertexArguments(_Base):

    def test_q_off_mesh_refused(self):
        from hwave.initial import hints
        f = _facts(self.golden())
        P = hints.pair_space(f, su2=_failed_su2())
        with self.assertRaisesRegex(ValueError, "not on the k-mesh"):
            hints.vertex(f, (Fraction(1, 4), 0, 0), P)
        with self.assertRaisesRegex(ValueError, "not on the k-mesh"):
            hints.vertex_sparse(f, (Fraction(1, 4), 0, 0), P)

    def test_foreign_pair_space_refused(self):
        from hwave.initial import hints
        f = _facts(self.golden())
        f2 = _facts(self.golden(cellshape=(6, 4, 1)))
        P2 = hints.pair_space(f2, su2=_failed_su2())
        with self.assertRaisesRegex(ValueError, "another model"):
            hints.vertex(f, (0, 0, 0), P2)
        f3 = _facts(self.golden(fock=False))
        P3 = hints.pair_space(f3, su2=_failed_su2())
        with self.assertRaisesRegex(ValueError, "another model"):
            hints.vertex_sparse(f, (0, 0, 0), P3)

    def test_unknown_sector_refused(self):
        from hwave.initial import hints
        f = _facts(self.golden())
        P = hints.pair_space(f, su2=_failed_su2())
        with self.assertRaisesRegex(ValueError, "unknown sector"):
            hints.vertex(f, (0, 0, 0), P, sector=hints.SECTOR_SZ)

    def test_folded_facts_use_the_primitive_cell(self):
        from hwave.initial import hints
        toml = self.golden(param={"SubShape": [3, 1, 1]})
        fB = _facts(toml)
        f1 = _facts(toml, subshape=(1, 1, 1))
        PB = hints.pair_space(fB, su2=_failed_su2())
        P1 = hints.pair_space(f1, su2=_failed_su2())
        np.testing.assert_array_equal(PB.labels, P1.labels)
        q = (Fraction(1, 3), 0, 0)
        np.testing.assert_array_equal(hints.vertex(fB, q, PB), hints.vertex(f1, q, P1))


if __name__ == "__main__":
    unittest.main()
