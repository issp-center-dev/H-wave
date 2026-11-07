"""The sector building blocks of the spin-orbital eigensolver in ``hwave.sc``.

* ``_SOAntisymmetricBasis.from_subset(shape, indices, reverse=None)``: the
  fermion-antisymmetric basis restricted to an index subset closed under the
  involution ``R``; its ``to_full`` / ``to_reduced`` act on the sector
  vector (length ``len(indices)``) and ``scatter`` places a sector vector in
  the full flat gap. A subset that is not closed is refused.
* ``basis.reduce(A, adjoint=True)`` attaches ``rmatvec = E^dagger K^dagger E``;
  the default call builds the forward operator only, as before.
* ``SectorOp(pp_sector, operator, basis, has_adjoint)``.
* ``_so_dense_reduced_eigen(red, dim, k, seed_vec)``: the dense branch of
  ``_so_solve_eigenvalue`` as a helper that returns the reduced matrix too
  and logs nothing.
* ``_so_solve_sectors`` / ``_so_iterate_sectors``: the eigen solve and the
  power iteration sector by sector, with the failure semantics of the state
  table (empty, null, not converged, residual failure, none certified,
  mixed), and ``_so_bond_report``: the one report of both layouts.
"""
import logging
import re
import unittest
from unittest import mock

import numpy as np
from scipy.sparse.linalg import LinearOperator

from hwave.solver.eliashberg import static as sc_static
from hwave.solver.eliashberg import eigen as sc_eigen
from hwave.solver.eliashberg import so_solve as sc_so
from hwave.solver.eliashberg import dynamic as eliashberg_dynamic
from tests.so_bond_fixtures import exact, pp_flat_index, quiet, subset_bases
from tests.so_bond_fixtures import sector_op_of as _sector

SHAPE = (4, 4, 3, 1, 1)                     # (nd, nd, Nx, Ny, Nz), nd = 4
SHAPE_W = (4, 4, 3, 1, 1, 4)                # with nfreq = 4


def _pp_index(shape):
    """Per pair-S_z sector (+1, 0, -1) the ascending flat gap indices of the
    pairs (b, c) with S_z(b) + S_z(c) = p (spin-block modes)."""
    return pp_flat_index(shape[0], int(np.prod(shape[2:])))


def _random_operator(n, seed):
    rng = np.random.default_rng(seed)
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    return M, LinearOperator((n, n), matvec=lambda x: M @ x, rmatvec=lambda y: M.conj().T @ y,
                             matmat=lambda X: M @ X, dtype=complex)


class TestFromSubset(unittest.TestCase):
    """Construction check: the subset basis on the pair-S_z sectors of a
    static and a frequency-resolved gap."""

    def _cases(self):
        return ((SHAPE, None), (SHAPE_W, eliashberg_dynamic._reverse_kw_and_orbital))

    def test_partition_and_coordinates(self):
        for shape, rev in self._cases():
            full = sc_eigen._SOAntisymmetricBasis(shape, reverse=rev)
            n_full = int(np.prod(shape))
            dims = 0
            for idx in _pp_index(shape):
                b = sc_eigen._SOAntisymmetricBasis.from_subset(shape, idx, reverse=rev)
                np.testing.assert_array_equal(b.indices, idx)
                self.assertEqual(b.n, len(idx))
                self.assertEqual(b.shape, tuple(shape))
                dims += b.dim
                rng = np.random.default_rng(len(idx))
                c = rng.normal(size=b.dim) + 1j * rng.normal(size=b.dim)
                v = b.to_full(c)
                self.assertEqual(v.shape, (len(idx),))
                np.testing.assert_allclose(b.to_reduced(v), c, rtol=0, atol=1e-15 * np.abs(c).max())
                g = b.scatter(v)
                self.assertEqual(g.shape, (n_full,))
                np.testing.assert_array_equal(g[idx], v)
                self.assertEqual(np.count_nonzero(np.delete(g, idx)), 0)
                # the scattered gap is fermion-antisymmetric on the full space
                np.testing.assert_allclose(full.to_full(full.to_reduced(g)), g, rtol=0,
                                           atol=1e-15 * np.abs(g).max())
                # matrix form: E has orthonormal columns
                E = b.to_full(np.eye(b.dim))
                self.assertEqual(E.shape, (len(idx), b.dim))
                np.testing.assert_allclose(E.conj().T @ E, np.eye(b.dim), atol=1e-15)
            self.assertEqual(dims, full.dim)

    def test_not_closed_refused(self):
        idx = _pp_index(SHAPE)[1]
        # drop the partner of the first non-self-paired index
        full = sc_eigen._SOAntisymmetricBasis(SHAPE)
        i = next(int(j) for j in idx if full.partner[j] != j)
        sub = np.setdiff1d(idx, [full.partner[i]])
        with self.assertRaisesRegex(ValueError, exact(
                "the subset of the antisymmetric basis is not closed under the involution "
                "(index {} maps to {})".format(i, int(full.partner[i])))):
            sc_eigen._SOAntisymmetricBasis.from_subset(SHAPE, sub)

    def test_bad_indices_refused(self):
        idx = _pp_index(SHAPE)[0]
        for bad in (idx[::-1], np.concatenate([idx, idx[:1]]), np.array([-1]),
                    np.array([int(np.prod(SHAPE))]), idx.reshape(1, -1)):
            with self.assertRaises(ValueError):
                sc_eigen._SOAntisymmetricBasis.from_subset(SHAPE, bad)

    def test_full_basis_scatter(self):
        full = sc_eigen._SOAntisymmetricBasis(SHAPE)
        self.assertIsNone(full.indices)
        v = np.arange(full.n, dtype=complex)
        np.testing.assert_array_equal(full.scatter(v), v)


class TestReduceAdjoint(unittest.TestCase):
    """Construction check: ``reduce(A, adjoint=True)`` carries the adjoint of
    the reduced operator; the default call carries none."""

    def test_adjoint(self):
        idx = _pp_index(SHAPE)[1]
        b = sc_eigen._SOAntisymmetricBasis.from_subset(SHAPE, idx)
        M, A = _random_operator(len(idx), 3)
        red = b.reduce(A, adjoint=True)
        E = b.to_full(np.eye(b.dim))
        R = E.conj().T @ M @ E
        cols = np.array([red.rmatvec(e) for e in np.eye(b.dim)]).T
        np.testing.assert_allclose(cols, R.conj().T, rtol=0, atol=1e-13 * np.abs(R).max())
        fwd = np.array([red.matvec(e) for e in np.eye(b.dim)]).T
        np.testing.assert_allclose(fwd, R, rtol=0, atol=1e-13 * np.abs(R).max())

    def test_default_has_no_adjoint(self):
        b = sc_eigen._SOAntisymmetricBasis(SHAPE)
        _, A = _random_operator(b.n, 4)
        red = b.reduce(A)
        self.assertEqual(red.shape, (b.dim, b.dim))
        with self.assertRaises(NotImplementedError):
            red.rmatvec(np.ones(b.dim, dtype=complex))

    def test_sector_op(self):
        b = sc_eigen._SOAntisymmetricBasis(SHAPE)
        _, A = _random_operator(b.n, 5)
        op = sc_so.SectorOp(0, A, b, True)
        self.assertEqual(op._fields, ("pp_sector", "operator", "basis", "has_adjoint"))
        self.assertIs(op.basis, b)


class TestDenseReducedEigen(unittest.TestCase):
    """Construction check: the extracted dense branch returns what
    ``_so_solve_eigenvalue`` returns on its dense branch, and logs nothing."""

    def test_matches_solver_and_is_silent(self):
        idx = _pp_index(SHAPE)[0]
        b = sc_eigen._SOAntisymmetricBasis.from_subset(SHAPE, idx)
        _, A = _random_operator(len(idx), 6)
        red = b.reduce(A)
        k = b.dim
        with self.assertLogs(sc_static.logger, level="INFO") as cm:
            vals0, vecs0 = sc_so._so_solve_eigenvalue(A, b, k, "arnoldi", map_back=False)
        self.assertTrue(any("dense diagonalization" in r for r in cm.output))
        logger = logging.getLogger()
        with self.assertNoLogs(logger, level="DEBUG"):
            vals, vecs, M = sc_so._so_dense_reduced_eigen(red, b.dim, k, None)
        np.testing.assert_array_equal(vals, vals0)
        np.testing.assert_array_equal(vecs, vecs0)
        np.testing.assert_array_equal(M, red.matmat(np.eye(b.dim, dtype=complex)))
        seed = np.random.default_rng(1).normal(size=b.dim) + 0j
        vals1, vecs1, _ = sc_so._so_dense_reduced_eigen(red, b.dim, b.dim - 1, seed)
        vals2, vecs2 = sc_so._so_solve_eigenvalue(A, b, b.dim - 1, "arnoldi", seed_vec=seed,
                                               map_back=False)
        np.testing.assert_array_equal(vals1, vals2)
        np.testing.assert_array_equal(vecs1, vecs2)


class TestClosureAndDims(unittest.TestCase):
    """Construction check: the pair-S_z sectors of a gap with one physical
    orbital on a 5 x 3 mesh partition the antisymmetric space."""

    def test_closure_and_dims(self):
        shape = (2, 2, 5, 3, 1)
        full = sc_eigen._SOAntisymmetricBasis(shape)
        subs = [sc_eigen._SOAntisymmetricBasis.from_subset(shape, idx) for idx in _pp_index(shape)]
        self.assertEqual(sum(b.dim for b in subs), full.dim)
        idx = _pp_index(shape)[1]
        i = next(int(j) for j in idx if full.partner[j] != j)
        with self.assertRaisesRegex(ValueError, exact(
                "the subset of the antisymmetric basis is not closed under the involution "
                "(index {} maps to {})".format(i, int(full.partner[i])))):
            sc_eigen._SOAntisymmetricBasis.from_subset(shape, np.setdiff1d(idx, [full.partner[i]]))


# --- the sector solvers --------------------------------------------------------------

GAP = (4, 4, 16, 1, 1)                      # full gap; antisymmetric dim 124
_FULL = sc_eigen._SOAntisymmetricBasis(GAP)
_EMPTY_MSG = ("spin-orbital: the antisymmetric gap space is empty on this grid; "
              "no eigenvalue is reported")


def _bases(dims):
    """Subset bases of the requested dimensions on disjoint orbits of GAP."""
    return subset_bases(GAP, dims)


def _normal(eigs, seed=0):
    """A normal matrix with the eigenvalues ``eigs``."""
    d = len(eigs)
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return (Q * np.asarray(eigs, dtype=complex)) @ Q.conj().T


def _scalar_sector(p, basis, lam):
    """``K = lam * 1`` on the sector vector (its reduced spectrum is ``lam``
    with every rounding identical across sectors)."""
    n = basis.n
    op = LinearOperator((n, n), dtype=complex, matvec=lambda x: lam * np.asarray(x).ravel(),
                        rmatvec=lambda y: np.conj(lam) * np.asarray(y).ravel())
    return sc_so.SectorOp(p, op, basis, True)


def _lift(basis, c):
    return basis.scatter(basis.to_full(np.asarray(c)))


def _unit_gap(basis, i):
    g = _lift(basis, np.eye(basis.dim)[i])
    g /= np.linalg.norm(g)
    return g.reshape(GAP)


def _seed():
    rng = np.random.default_rng(11)
    c = rng.normal(size=_FULL.dim) + 1j * rng.normal(size=_FULL.dim)
    return _FULL.to_full(c).reshape(GAP)


def _solve(ops, k=10, **kw):
    return sc_so._so_solve_sectors(ops, k, "arnoldi", gap_shape=GAP, **kw)


def _iterate(ops, sigma_init=None, **kw):
    kw.setdefault("spectral_shift", "auto")
    return sc_so._so_iterate_sectors(ops, _seed() if sigma_init is None else sigma_init,
                                  gap_shape=GAP, max_iter=400, alpha=0.0, tol=1e-9, **kw)


def _report(mode, eigen=None, iteration=None):
    e = None if eigen is None else sc_so.EigenLeg.from_sectors(eigen)
    i = None if iteration is None else sc_so.IterationLeg.from_sectors(iteration)
    return sc_so._so_bond_report(mode, e, i, gap_shape=GAP)


_REAL_SOLVE_LEADING = sc_eigen._solve_leading


def _failing_arnoldi(dim):
    """``_solve_leading`` with the eigen solve of a ``dim``-dimensional
    operator replaced by the not-converged sentinel."""
    def fake(make_operator, vec_size, solver_mode, *a, **kw):
        if solver_mode == "arnoldi" and vec_size == dim:
            k = kw.get("num_eigenvalues", 10)
            w = np.eye(vec_size, k, dtype=complex)
            v = np.full(k, np.nan + 0j)
            return v[0], w[:, 0], {"eigenvalues": v, "eigenvectors": w,
                                   "selection": "not-converged"}
        return _REAL_SOLVE_LEADING(make_operator, vec_size, solver_mode, *a, **kw)
    return fake


def _last_iterate(dim):
    rng = np.random.default_rng(dim)
    v = rng.normal(size=dim) + 1j * rng.normal(size=dim)
    return v / np.linalg.norm(v)


def _failing_iteration(*dims):
    """``_solve_leading`` with the power iteration of operators of the
    dimensions ``dims`` replaced by a non-converged run ending at
    ``_last_iterate(dim)``."""
    def fake(make_operator, vec_size, solver_mode, *a, **kw):
        if solver_mode == "iteration" and vec_size in dims:
            return float("nan"), _last_iterate(vec_size), {"converged": False, "n_iter": 7}
        return _REAL_SOLVE_LEADING(make_operator, vec_size, solver_mode, *a, **kw)
    return fake


_REAL_DENSE = sc_so._so_dense_reduced_eigen


def _perturbed_dense(dim, which=1):
    """The dense branch with the eigenvalue of pair ``which`` of a
    ``dim``-dimensional sector moved by 1e-6 (its vector no longer matches)."""
    def fake(red, d, k, seed_vec=None):
        vals, vecs, M = _REAL_DENSE(red, d, k, seed_vec)
        if d == dim:
            vals = vals.copy()
            vals[which] += 1e-6
        return vals, vecs, M
    return fake


def _messages(cm):
    return [r.getMessage() for r in cm.records]


class TestStateTable(unittest.TestCase):
    """Gate: the failure semantics of the sector solvers and the report."""

    def test_row1_empty(self):
        b0, b1 = _bases((0, 5))
        ops = [_sector(1, b0, np.zeros((0, 0))), _sector(0, b1, _normal([3, 1, 0.5, -1, -2]))]
        s = _solve(ops)
        self.assertEqual(s.leaders[0].selection, "empty")
        self.assertTrue(s.leaders[0].certified)
        self.assertEqual((s.leaders[0].dim, s.leaders[0].n_returned), (0, 0))
        self.assertEqual(s.failed_sectors, ())
        self.assertTrue(s.converged)
        self.assertAlmostEqual(s.leading, 3.0, delta=1e-12)
        np.testing.assert_array_equal(s.pp_sector, np.zeros(5, np.int8))
        self.assertEqual(s.selection, "+1:empty;0:dense-LR")
        it = _iterate(ops)
        self.assertEqual([lead.status for lead in it.leaders], ["empty", "converged"])
        self.assertEqual(it.failed_sectors, ())
        self.assertTrue(it.converged)
        self.assertAlmostEqual(it.value, 3.0, delta=1e-7)
        self.assertEqual(it.pp_sector, 0)

    def test_row2_null_operator(self):
        for dim in (7, 40):
            (b,) = _bases((dim,))
            op = _sector(1, b, np.zeros((dim, dim)))
            s1, s2 = _solve([op], k=10), _solve([op], k=10)
            self.assertEqual(s1.leaders[0].selection, "null-operator")
            self.assertEqual(s1.failed_sectors, ())
            self.assertTrue(s1.converged)
            self.assertEqual((s1.leaders[0].norm_kind, s1.leaders[0].norm_estimate), ("exact", 0.0))
            n = min(10, dim)
            self.assertTrue(np.array_equal(s1.vals, np.zeros(n)))
            self.assertTrue(np.array_equal(s1.pp_sector, np.full(n, 1, np.int8)))
            for i in range(n):
                want = _unit_gap(b, i)
                self.assertTrue(np.array_equal(s1.gaps[i], want))
                self.assertAlmostEqual(np.linalg.norm(s1.gaps[i]), 1.0, delta=1e-15)
            self.assertTrue(np.array_equal(s1.gaps, s2.gaps))
            with mock.patch.object(sc_so, "_solve_leading", wraps=sc_eigen._solve_leading) as spy:
                it = _iterate([op])
            self.assertEqual(spy.call_count, 0)
            lead = it.leaders[0]
            self.assertEqual((lead.status, lead.value, lead.n_iter), ("null", 0.0, 0))
            self.assertTrue(np.array_equal(lead.gap, _unit_gap(b, 0)))
            self.assertTrue(it.converged)
            self.assertEqual((it.value, it.pp_sector, it.selection), (0.0, 1, "null-operator"))
            self.assertTrue(np.array_equal(it.gap, _unit_gap(b, 0)))

    def test_row3_not_converged(self):
        # row 3: a sector whose eigen solve does not converge is failed (its
        # pairs are not listed, leading = NaN); the other sectors still list
        # their certified pairs, which the report writes as "uncertified"
        b1, b0 = _bases((22, 5))
        ops = [_sector(1, b1, _normal(np.concatenate([[5.0], np.linspace(-1, 1, 21)]), 1)),
               _sector(0, b0, _normal([0.9, 0.5, 0.1, -0.3, -1], 2))]
        with mock.patch.object(sc_so, "_solve_leading", _failing_arnoldi(22)), \
                self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve(ops)
        self.assertEqual(_messages(cm), ["spin-orbital sector +1: the eigen solve did not "
                                         "converge; the sector is failed."])
        self.assertEqual(s.leaders[0].selection, "not-converged")
        self.assertFalse(s.leaders[0].certified)
        self.assertEqual(s.failed_sectors, (1,))
        self.assertFalse(s.converged)
        self.assertTrue(np.isnan(s.leading))
        np.testing.assert_array_equal(s.pp_sector, np.zeros(5, np.int8))
        self.assertAlmostEqual(s.vals[0], 0.9, delta=1e-12)
        r = _report("eigenvalue", s)
        self.assertTrue(np.isnan(r.reported))
        self.assertFalse(r.eigen_converged)
        self.assertEqual((r.gap_status, r.pp_sector), ("uncertified", 0))
        self.assertTrue(np.array_equal(r.gap, s.gaps[0]))
        self.assertEqual(r.failed_sectors, (1,))
        # "both": the eigen leg failed -> NaN, the iteration's gap
        it = _iterate(ops)
        self.assertTrue(it.converged)
        r = _report("both", s, it)
        self.assertTrue(np.isnan(r.reported))
        self.assertFalse(r.eigen_converged)
        self.assertTrue(r.iteration_converged)
        self.assertEqual((r.gap_status, r.pp_sector), ("uncertified", it.pp_sector))
        self.assertTrue(np.array_equal(r.gap, it.gap))
        self.assertEqual((r.failed_sectors, r.iteration_failed_sectors), ((1,), ()))

    def test_row4_residual_any_pair(self):
        b1, b0 = _bases((5, 4))
        M1 = _normal([2.0, 1.0, 0.5, -0.5, -1.0], 3)
        ops = [_sector(1, b1, M1), _sector(0, b0, _normal([0.7, 0.2, -0.1, -0.6], 4))]
        with mock.patch.object(sc_so, "_so_dense_reduced_eigen", _perturbed_dense(5, which=1)), \
                self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve(ops)
        self.assertEqual(len(cm.records), 1)
        self.assertTrue(re.match(r"^spin-orbital sector \+1: a returned eigenpair fails the "
                                 r"residual check \(\d\.\d{3}e[+-]\d+ > 1e-10 relative to the "
                                 r"scale \d\.\d{6}e[+-]\d+\); the sector is failed\.$",
                                 _messages(cm)[0]), _messages(cm)[0])
        lead = s.leaders[0]
        self.assertEqual(lead.selection, "residual-failed")
        self.assertFalse(lead.certified)
        self.assertEqual(lead.n_returned, 5)
        self.assertGreater(lead.residual, 1e-8)
        self.assertEqual(s.failed_sectors, (1,))
        self.assertFalse(s.converged)
        np.testing.assert_array_equal(s.pp_sector, np.zeros(4, np.int8))
        # the leading pair of the failed sector alone passes: the failure is
        # the non-leading pair's
        red = b1.reduce(ops[0].operator)
        vals, vecs, M = _REAL_DENSE(red, 5, 10)
        tol = 1e-10 * max(1.0, abs(vals[0]), np.linalg.norm(M, 2))
        self.assertLessEqual(np.linalg.norm(M @ vecs[:, 0] - vals[0] * vecs[:, 0]), tol)
        r = _report("eigenvalue", s)
        self.assertTrue(np.isnan(r.reported))
        self.assertEqual((r.gap_status, r.pp_sector), ("uncertified", 0))

    def test_row5_all_empty(self):
        ops = [_sector(p, b, np.zeros((0, 0))) for p, b in zip((1, 0, -1), _bases((0, 0, 0)))]
        with self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve(ops)
        self.assertEqual(_messages(cm), [_EMPTY_MSG])
        self.assertEqual(s.vals.shape, (0,))
        self.assertEqual(s.gaps.shape, (0,) + GAP)
        self.assertEqual((s.pp_sector.shape, s.pp_sector.dtype), ((0,), np.int8))
        self.assertFalse(s.converged)
        self.assertTrue(np.isnan(s.leading))
        self.assertEqual(s.failed_sectors, ())
        with self.assertNoLogs(sc_static.logger, level="WARNING"):
            it = _iterate(ops)          # the empty-space WARNING is the eigen leg's
        self.assertTrue(np.isnan(it.value))
        self.assertFalse(it.converged)
        self.assertIsNone(it.gap)
        self.assertEqual((it.pp_sector, it.n_iter, it.n_iter_per_sector), (2, 0, (0, 0, 0)))
        for mode in ("eigenvalue", "iteration", "both"):
            r = _report(mode, s if mode != "iteration" else None, it if mode != "eigenvalue" else None)
            self.assertTrue(np.isnan(r.reported), mode)
            self.assertFalse(r.eigen_converged, mode)
            self.assertEqual((r.gap_status, r.pp_sector), ("none", 2), mode)
            self.assertEqual(r.gap.shape, GAP)
            self.assertEqual(np.count_nonzero(r.gap), 0)
            self.assertEqual(r.iteration_converged, False if mode == "both" else None)

    def test_row6_no_certified_pair(self):
        b1, b0 = _bases((5, 0))
        ops = [_sector(1, b1, _normal([2.0, 1.0, 0.5, -0.5, -1.0], 3)),
               _sector(0, b0, np.zeros((0, 0)))]
        with mock.patch.object(sc_so, "_so_dense_reduced_eigen", _perturbed_dense(5)), \
                self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve(ops)
        self.assertNotIn(_EMPTY_MSG, _messages(cm))
        self.assertEqual((s.vals.shape, s.gaps.shape), ((0,), (0,) + GAP))
        self.assertFalse(s.converged)
        self.assertEqual(s.failed_sectors, (1,))
        with mock.patch.object(sc_so, "_solve_leading", _failing_iteration(5)):
            it = _iterate(ops)
        self.assertEqual(it.failed_sectors, (1,))
        r = _report("eigenvalue", s)
        self.assertTrue(np.isnan(r.reported))
        self.assertEqual((r.gap_status, r.pp_sector), ("none", 2))
        self.assertEqual(np.count_nonzero(r.gap), 0)
        last = _lift(b1, _last_iterate(5)).reshape(GAP)
        for mode in ("iteration", "both"):
            r = _report(mode, s if mode == "both" else None, it)
            self.assertTrue(np.isnan(r.reported), mode)
            self.assertFalse(r.eigen_converged, mode)
            self.assertEqual((r.gap_status, r.pp_sector), ("uncertified", 1), mode)
            self.assertTrue(np.array_equal(r.gap, last), mode)

    def test_row7_mixed_null_and_converged(self):
        bn, b0 = _bases((4, 5))
        for lam, want, want_p in ((-0.4, 0.0, 1), (0.3, 0.3, 0)):
            eigs = [lam, -2.0, -2.5, -3.0, -3.5]
            ops = [_sector(1, bn, np.zeros((4, 4))), _sector(0, b0, _normal(eigs, 5))]
            s = _solve(ops)
            self.assertTrue(s.converged)
            self.assertAlmostEqual(s.leading, want, delta=1e-12)
            self.assertEqual(int(s.pp_sector[0]), want_p)
            self.assertEqual(len(s.vals), 9)
            r = _report("eigenvalue", s)
            self.assertAlmostEqual(r.reported, want, delta=1e-12)
            self.assertTrue(r.eigen_converged)
            self.assertEqual(r.gap_status, "certified")
            it = _iterate(ops)
            self.assertTrue(it.converged)
            self.assertAlmostEqual(it.value, want, delta=1e-7)
            self.assertEqual(it.pp_sector, want_p)
            self.assertEqual(it.selection, "null-operator" if want_p == 1 else "iteration")

    def test_both_mode(self):
        b1, b0 = _bases((5, 6))
        ops = [_sector(1, b1, _normal([0.8, 0.1, -0.2, -0.5, -1], 6)),
               _sector(0, b0, _normal([1.5, 0.4, 0.2, -0.3, -0.7, -1.2], 7))]
        s = _solve(ops)
        with mock.patch.object(sc_so, "_solve_leading", _failing_iteration(5)):
            it = _iterate(ops)
        self.assertFalse(it.converged)
        r = _report("both", s, it)
        self.assertAlmostEqual(r.reported, 1.5, delta=1e-12)
        self.assertTrue(r.eigen_converged)
        self.assertFalse(r.iteration_converged)
        self.assertEqual((r.gap_status, r.pp_sector), ("certified", 0))
        self.assertTrue(np.array_equal(r.gap, s.gaps[0]))
        self.assertEqual((r.failed_sectors, r.iteration_failed_sectors), ((), (1,)))
        self.assertTrue(np.array_equal(r.vals, s.vals))

    def test_gap_status(self):
        g = np.ones(GAP, complex)
        e_empty = sc_so.EigenLeg(np.zeros(0, complex), np.zeros((0,) + GAP, complex),
                              np.zeros(0, np.int8), False, complex(np.nan), "", ())
        i_empty = sc_so.IterationLeg(float("nan"), None, 2, False, 0, (), "iteration")
        i_fail = sc_so.IterationLeg(float("nan"), g, 0, False, 9, (0,), "iteration")
        cases = (("eigenvalue", e_empty, None, "none"), ("iteration", None, i_empty, "none"),
                 ("both", e_empty, i_empty, "none"), ("iteration", None, i_fail, "uncertified"),
                 ("both", e_empty, i_fail, "uncertified"))
        for mode, e, i, want in cases:
            r = sc_so._so_bond_report(mode, e, i, gap_shape=GAP)
            self.assertEqual(r.gap_status, want, (mode, want))
            self.assertEqual(r.gap.shape, GAP)
            self.assertEqual(r.gap_status == "certified", bool(np.isfinite(r.reported)))

    def test_report_reductions(self):
        g = [np.full(GAP, i + 1, complex) for i in range(5)]
        E_OK = sc_so.EigenLeg(np.array([0.5, 0.2], complex), np.array([g[0], g[1]]),
                           np.array([0, 1], np.int8), True, 0.5 + 0j, "+1:dense-LR;0:dense-LR", ())
        E_FAIL = sc_so.EigenLeg(np.array([0.2], complex), np.array([g[1]]), np.array([1], np.int8),
                             False, complex(np.nan), "+1:dense-LR;0:not-converged", (0,))
        I_OK = sc_so.IterationLeg(0.49, g[2], 0, True, 10, (), "iteration")
        I_FAIL = sc_so.IterationLeg(float("nan"), g[3], -1, False, 50, (-1,), "iteration")
        I_NULL = sc_so.IterationLeg(0.0, g[4], 1, True, 0, (), "null-operator")
        nan = float("nan")
        # (mode, eigen, iteration) -> (reported, eigen_converged, iteration_converged,
        #  gap, gap_status, pp_sector, selection, failed_sectors, iteration_failed_sectors)
        table = [
            ("eigenvalue", E_OK, None, (0.5, True, None, 0, "certified", 0, E_OK.selection, (), None)),
            ("eigenvalue", E_FAIL, None, (nan, False, None, 1, "uncertified", 1, E_FAIL.selection,
                                          (0,), None)),
            ("iteration", None, I_OK, (0.49, True, None, 2, "certified", 0, "iteration", (), ())),
            ("iteration", None, I_FAIL, (nan, False, None, 3, "uncertified", -1, "iteration",
                                         (-1,), (-1,))),
            ("iteration", None, I_NULL, (0.0, True, None, 4, "certified", 1, "null-operator", (),
                                         ())),
        ]
        for it, ic, gi in ((I_OK, True, 2), (I_FAIL, False, 3), (I_NULL, True, 4)):
            table.append(("both", E_OK, it, (0.5, True, ic, 0, "certified", 0, E_OK.selection, (),
                                             it.failed_sectors)))
            table.append(("both", E_FAIL, it, (nan, False, ic, gi, "uncertified", it.pp_sector,
                                               E_FAIL.selection, (0,), it.failed_sectors)))
        for mode, e, i, want in table:
            r = sc_so._so_bond_report(mode, e, i, gap_shape=GAP)
            got = (r.reported, r.eigen_converged, r.iteration_converged, r.gap, r.gap_status,
                   r.pp_sector, r.selection, r.failed_sectors, r.iteration_failed_sectors)
            msg = (mode, e and e.selection, i and i.selection)
            if np.isnan(want[0]):
                self.assertTrue(np.isnan(got[0]), msg)
            else:
                self.assertEqual(got[0], want[0], msg)
            self.assertEqual(got[1:3], want[1:3], msg)
            self.assertTrue(np.array_equal(got[3], g[want[3]]), msg)
            self.assertEqual(got[4:], want[4:], msg)
            if mode == "iteration":
                self.assertEqual((r.vals.shape, r.gaps.shape, r.pp_sectors.shape),
                                 ((0,), (0,) + GAP, (0,)))
                self.assertEqual(r.pp_sectors.dtype, np.int8)
            else:
                self.assertTrue(np.array_equal(r.vals, e.vals))
                self.assertTrue(np.array_equal(r.gaps, e.gaps))
                self.assertTrue(np.array_equal(r.pp_sectors, e.pp_sector))
            self.assertIsInstance(r.reported, float)

    def test_sector_ops_validation(self):
        b1, b0 = _bases((3, 3))
        op = _sector(2, b1, np.eye(3))
        for solver, call in (("_so_solve_sectors", lambda ops, **kw: _solve(ops, **kw)),
                             ("_so_iterate_sectors", lambda ops, **kw: _iterate(ops, **kw))):
            with self.assertRaisesRegex(ValueError, exact(
                    "{}: pp_sector must be +1, 0 or -1, got 2".format(solver))):
                call([op])
            with self.assertRaisesRegex(ValueError, exact(
                    "{}: pp_sector 0 appears twice".format(solver))):
                call([_sector(0, b1, np.eye(3)), _sector(0, b0, np.eye(3))])
            e = _bases((0, 0, 0))
            empty = [_sector(p, b, np.zeros((0, 0))) for p, b in zip((1, 0, -1), e)]
            null = [_sector(1, b1, np.zeros((3, 3))), _sector(-1, b0, np.zeros((3, 3)))]
            for ops in (empty, null):
                with self.assertRaisesRegex(ValueError, "spectral_shift must be a positive finite"):
                    call(ops, spectral_shift=-1.0)
        bm, bp, bz = _bases((3, 4, 5))
        ops = [_sector(-1, bm, _normal([1, 2, 3], 1)), _sector(1, bp, _normal([1, 2, 3, 4], 2)),
               _sector(0, bz, _normal([1, 2, 3, 4, 5], 3))]
        self.assertEqual([lead.pp_sector for lead in _solve(ops).leaders], [1, 0, -1])
        self.assertEqual([lead.pp_sector for lead in _iterate(ops).leaders], [1, 0, -1])


class TestSectorIteration(unittest.TestCase):
    """Gate: the per-sector power iteration and its reduction."""

    def test_reduction(self):
        b1, b0, bm = _bases((5, 6, 3))
        ops = [_sector(1, b1, _normal([1.0, 0.3, -0.2, -0.4, -0.5], 8)),
               _sector(0, b0, _normal([2.0, 0.5, 0.3, -0.1, -0.4, -0.6], 9)),
               _sector(-1, bm, np.zeros((3, 3)))]
        it = _iterate(ops)
        self.assertTrue(it.converged)
        self.assertAlmostEqual(it.value, 2.0, delta=1e-7)
        self.assertEqual(it.pp_sector, 0)
        self.assertEqual(it.selection, "iteration")
        self.assertTrue(np.array_equal(it.gap, it.leaders[1].gap))
        self.assertEqual(it.gap.shape, GAP)
        n1, n0 = it.leaders[0].n_iter, it.leaders[1].n_iter
        self.assertGreater(min(n1, n0), 0)
        self.assertEqual(it.n_iter_per_sector, (n1, n0, 0))
        self.assertEqual(it.n_iter, max(n1, n0))
        # the gap is the leading eigenvector of its sector, lifted
        red = b0.reduce(ops[1].operator)
        c = b0.to_reduced(it.gap.ravel()[b0.indices])
        self.assertLess(np.linalg.norm(red.matvec(c) - 2.0 * c), 1e-6 * np.linalg.norm(c))

    def test_iteration_leaders(self):
        b1, b0, bm = _bases((0, 4, 5))
        ops = [_sector(1, b1, np.zeros((0, 0))), _sector(0, b0, np.zeros((4, 4))),
               _sector(-1, bm, _normal([1.0, 0.5, 0.2, -0.3, -0.6], 10))]
        with mock.patch.object(sc_so, "_solve_leading", _failing_iteration(5)):
            it = _iterate(ops)
        e, z, f = it.leaders
        self.assertEqual((e.pp_sector, e.status, e.gap, e.n_iter), (1, "empty", None, 0))
        self.assertTrue(np.isnan(e.value))
        self.assertEqual((z.pp_sector, z.status, z.value, z.n_iter), (0, "null", 0.0, 0))
        self.assertTrue(np.array_equal(z.gap, _unit_gap(b0, 0)))
        self.assertEqual((f.pp_sector, f.status, f.n_iter), (-1, "failed", 7))
        self.assertTrue(np.isnan(f.value))
        self.assertTrue(np.array_equal(f.gap, _lift(bm, _last_iterate(5)).reshape(GAP)))
        self.assertEqual(it.n_iter_per_sector, (0, 0, 7))

    def test_dim_one_exact(self):
        b1, b0 = _bases((1, 1))
        ops = [_sector(1, b1, [[0.25]]), _sector(0, b0, [[-0.5]])]
        with mock.patch.object(sc_so, "_solve_leading", wraps=sc_eigen._solve_leading) as spy:
            it = _iterate(ops)
        self.assertEqual(spy.call_count, 0)
        self.assertEqual([lead.status for lead in it.leaders], ["converged", "converged"])
        self.assertEqual([lead.n_iter for lead in it.leaders], [0, 0])
        self.assertAlmostEqual(it.value, 0.25, delta=1e-15)
        self.assertEqual(it.pp_sector, 1)
        self.assertTrue(np.array_equal(it.gap, _unit_gap(b1, 0)))
        self.assertTrue(it.converged)

    def test_dim_one_non_finite_is_failed(self):
        b1, b0 = _bases((1, 1))
        for bad in (np.nan, np.inf, complex(1.0, np.inf)):
            with self.subTest(element=bad):
                ops = [_scalar_sector(1, b1, bad), _scalar_sector(0, b0, 0.25)]
                with self.assertLogs(sc_static.logger, level="WARNING") as cm, \
                        np.errstate(invalid="ignore"):
                    it = _iterate(ops)
                self.assertEqual(len(cm.records), 1)
                self.assertIn("one-dimensional sector is not finite", _messages(cm)[0])
                lead = it.leaders[0]
                self.assertEqual((lead.status, lead.n_iter), ("failed", 0))
                self.assertTrue(np.isnan(lead.value))
                self.assertTrue(np.array_equal(lead.gap, _unit_gap(b1, 0)))
                self.assertEqual(it.failed_sectors, (1,))
                self.assertFalse(it.converged)
                self.assertTrue(np.isnan(it.value))
                s = _solve([_scalar_sector(0, b0, 0.25)])
                for mode in ("iteration", "both"):
                    r = _report(mode, s if mode == "both" else None, it)
                    self.assertEqual(r.iteration_failed_sectors, (1,), mode)
                    if mode == "iteration":
                        self.assertTrue(np.isnan(r.reported))
                        self.assertFalse(r.eigen_converged)
                        self.assertEqual((r.gap_status, r.pp_sector), ("uncertified", 1))
                    else:
                        self.assertFalse(r.iteration_converged)
                        self.assertTrue(r.eigen_converged)

    def test_non_finite_converged_value_is_failed(self):
        """A sector (dim > 1) whose power iteration reports convergence with
        a non-finite value is failed, also when it holds the true leader."""
        b1, b0 = _bases((5, 6))
        ops = [_sector(1, b1, _normal([0.8, 0.1, -0.2, -0.5, -1], 6)),
               _sector(0, b0, _normal([1.5, 0.4, 0.2, -0.3, -0.7, -1.2], 7))]
        s = _solve(ops)
        self.assertEqual(int(s.pp_sector[0]), 0)           # sector 0 holds the leader

        for bad in (float("nan"), float("inf")):
            def fake(make_operator, vec_size, solver_mode, *a, **kw):
                if solver_mode == "iteration" and vec_size == 6:
                    return bad, _last_iterate(6), {"converged": True, "n_iter": 3}
                return _REAL_SOLVE_LEADING(make_operator, vec_size, solver_mode, *a, **kw)
            with self.subTest(value=bad):
                with mock.patch.object(sc_so, "_solve_leading", fake):
                    it = _iterate(ops)
                lead = it.leaders[1]
                self.assertEqual((lead.pp_sector, lead.status, lead.n_iter), (0, "failed", 3))
                self.assertTrue(np.isnan(lead.value))
                self.assertEqual(it.failed_sectors, (0,))
                self.assertFalse(it.converged)
                self.assertTrue(np.isnan(it.value))
                for mode in ("iteration", "both"):
                    r = _report(mode, s if mode == "both" else None, it)
                    self.assertEqual(r.iteration_failed_sectors, (0,), mode)
                    if mode == "iteration":
                        self.assertFalse(r.eigen_converged)
                        self.assertTrue(np.isnan(r.reported))
                    else:
                        self.assertFalse(r.iteration_converged)

    def test_report_requires_a_finite_iteration_value(self):
        g = np.ones(GAP, complex)
        e_ok = sc_so.EigenLeg(np.array([0.5], complex), np.array([g]), np.array([0], np.int8),
                           True, 0.5 + 0j, "0:dense-LR", ())
        for bad in (float("nan"), float("inf"), complex(0.1, np.inf)):
            leg = sc_so.IterationLeg(bad, g, 0, True, 3, (), "iteration")
            r = sc_so._so_bond_report("iteration", None, leg, gap_shape=GAP)
            self.assertTrue(np.isnan(r.reported), bad)
            self.assertFalse(r.eigen_converged, bad)
            self.assertEqual(r.gap_status, "uncertified", bad)
            r = sc_so._so_bond_report("both", e_ok, leg, gap_shape=GAP)
            self.assertFalse(r.iteration_converged, bad)
            self.assertEqual((r.reported, r.gap_status), (0.5, "certified"), bad)

    def test_equal_leaders_deterministic(self):
        b1, b0, bm = _bases((1, 1, 1))
        ops = [_scalar_sector(-1, bm, 0.7), _scalar_sector(0, b0, 0.5), _scalar_sector(1, b1, 0.7)]
        it = _iterate(ops)
        self.assertEqual(it.pp_sector, 1)
        ops = [_scalar_sector(1, b1, 0.7), _scalar_sector(-1, bm, 0.7 + 0.1j)]
        it = _iterate(ops)
        self.assertEqual(it.pp_sector, -1)
        self.assertEqual(np.real(it.value), np.real(it.leaders[0].value))

    def test_failed_sector(self):
        b1, b0, bm = _bases((5, 6, 4))
        ops = [_sector(1, b1, _normal([3.0, 0.3, -0.2, -0.4, -0.5], 8)),
               _sector(0, b0, _normal([2.0, 0.5, 0.3, -0.1, -0.4, -0.6], 9)),
               _sector(-1, bm, _normal([1.0, 0.5, 0.3, -0.1], 12))]
        with mock.patch.object(sc_so, "_solve_leading", _failing_iteration(6)):
            it = _iterate(ops)
        self.assertEqual(it.failed_sectors, (0,))
        self.assertFalse(it.converged)
        self.assertTrue(np.isnan(it.value))
        self.assertEqual(it.pp_sector, 0)
        self.assertTrue(np.array_equal(it.gap, _lift(b0, _last_iterate(6)).reshape(GAP)))
        r = _report("iteration", None, it)
        self.assertTrue(np.isnan(r.reported))
        self.assertFalse(r.eigen_converged)
        self.assertEqual((r.gap_status, r.pp_sector, r.failed_sectors), ("uncertified", 0, (0,)))
        # two failed sectors: the gap is the first one's in the order (+1, 0, -1)
        with mock.patch.object(sc_so, "_solve_leading", _failing_iteration(4, 6)):
            it = _iterate(ops[::-1])
        self.assertEqual(it.failed_sectors, (0, -1))
        self.assertEqual(it.pp_sector, 0)
        self.assertTrue(np.array_equal(it.gap, _lift(b0, _last_iterate(6)).reshape(GAP)))


class TestSectorSolve(unittest.TestCase):
    """Gate: the per-sector eigen solve, its certification and its union."""

    def test_arnoldi_only(self):
        (b,) = _bases((3,))
        for method in ("shift-invert-gmres", "subspace"):
            with self.assertRaisesRegex(ValueError, exact(
                    "_so_solve_sectors: the per-sector solve supports eigenvalue_method = "
                    "\"arnoldi\" only, got {}".format(method))):
                sc_so._so_solve_sectors([_sector(1, b, np.eye(3))], 2, method, gap_shape=GAP)

    def test_dense_branch_rule(self):
        retry_eigs = np.concatenate([-np.arange(1.0, 11.0), 0.01 * np.arange(1.0, 13.0)])
        lm = np.sort(retry_eigs[np.argsort(-np.abs(retry_eigs))[:10]])
        # the real trigger of the largest-real-part re-solve: the 10
        # largest-magnitude eigenvalues all have real part below their
        # smallest magnitude
        self.assertLess(lm.real.max(), np.abs(lm).min())
        cases = ((20, np.linspace(-1, 2, 20), {}, "dense-LR", 2.0),
                 (21, np.linspace(-1, 2, 21), {}, "dense-LR", 2.0),
                 (22, np.linspace(0.1, 2.2, 22), {}, "LM", 2.2),
                 (22, retry_eigs, {}, "LR_retry", 0.12),
                 (22, np.linspace(-2, 1, 22), {"spectral_shift": 2.5}, "LR", 1.0))
        for dim, eigs, kw, sel, lead in cases:
            (b,) = _bases((dim,))
            M = _normal(eigs, dim)
            if sel == "LR_retry":
                ev = np.linalg.eigvals(M)
                top = ev[np.argsort(-np.abs(ev))[:10]]
                self.assertLess(top.real.max(), np.abs(top).min())
            if sel == "LR_retry":
                with self.assertLogs(sc_static.logger, level="WARNING") as cm:
                    s = _solve([_sector(0, b, M)], k=10, **kw)
                self.assertEqual(len(cm.records), 1)
                self.assertIn("re-solving for the largest real part", _messages(cm)[0])
            else:
                s = _solve([_sector(0, b, M)], k=10, **kw)
            self.assertEqual(s.leaders[0].selection, sel, (dim, sel))
            self.assertTrue(s.leaders[0].certified, (dim, sel))
            self.assertEqual(s.leaders[0].n_returned, min(10, dim))
            self.assertEqual(len(s.vals), min(10, dim))
            self.assertEqual(s.leaders[0].norm_kind, "exact" if sel == "dense-LR" else "adjoint")
            self.assertAlmostEqual(s.leading, lead, delta=1e-9, msg=(dim, sel))
            self.assertEqual(s.selection, "0:" + sel)

    def test_tie_break_and_order(self):
        b1, b0, bm = _bases((3, 3, 3))
        ops = [_scalar_sector(-1, bm, 1.0), _scalar_sector(0, b0, 1.0), _scalar_sector(1, b1, 1.0)]
        s = _solve(ops)
        np.testing.assert_array_equal(s.pp_sector, np.repeat(np.array([1, 0, -1], np.int8), 3))
        ops = [_scalar_sector(1, b1, 1.0), _scalar_sector(-1, bm, 1.0 + 0.5j)]
        s = _solve(ops, k=4)
        np.testing.assert_array_equal(s.pp_sector, np.array([-1, -1, -1, 1], np.int8))
        self.assertTrue(s.converged)
        self.assertAlmostEqual(s.leading, 1.0 + 0.5j, delta=1e-14)
        # the union is ordered by descending real part across sectors
        b1, b0 = _bases((4, 4))
        ops = [_sector(1, b1, _normal([0.9, 0.3, -0.2, -1.0], 1)),
               _sector(0, b0, _normal([1.1, 0.5, -0.1, -2.0], 2))]
        s = _solve(ops, k=5)
        np.testing.assert_allclose(s.vals.real, [1.1, 0.9, 0.5, 0.3, -0.1], atol=1e-12)
        np.testing.assert_array_equal(s.pp_sector, np.array([0, 1, 0, 1, 0], np.int8))
        # the gaps are the sector eigenvectors scattered to the full gap
        for lam, gap, p in zip(s.vals, s.gaps, s.pp_sector):
            b = b1 if p == 1 else b0
            self.assertEqual(np.count_nonzero(np.delete(gap.ravel(), b.indices)), 0)
            c = b.to_reduced(gap.ravel()[b.indices])
            red = b.reduce(ops[0 if p == 1 else 1].operator)
            self.assertLess(np.linalg.norm(red.matvec(c) - lam * c), 1e-12)

    def test_seed_split(self):
        b1, b0, bm = _bases((5, 30, 4))
        ops = [_sector(1, b1, _normal([1, 0.5, 0.2, -0.1, -0.3], 1)),
               _sector(0, b0, _normal(np.linspace(-1, 2, 30), 2)),
               _sector(-1, bm, _normal([0.4, 0.1, -0.2, -0.5], 3))]
        seed = _seed().ravel().copy()
        seed[bm.indices] = 0.0
        seed = seed.reshape(GAP)
        with mock.patch.object(sc_so, "_so_dense_reduced_eigen", wraps=sc_so._so_dense_reduced_eigen) as d, \
                mock.patch.object(sc_so, "_so_solve_eigenvalue", wraps=sc_so._so_solve_eigenvalue) as a, \
                self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve(ops, k=4, seed_vec=seed)
        # a random seed: no selected eigenvector is near it
        self.assertEqual(len(cm.records), 1)
        self.assertIn("the seed's branch is not among the selected eigenpairs", _messages(cm)[0])
        self.assertTrue(s.converged)
        dense_seeds = [c.args[3] if len(c.args) > 3 else c.kwargs.get("seed_vec")
                       for c in d.call_args_list]
        self.assertEqual(len(dense_seeds), 2)
        np.testing.assert_array_equal(dense_seeds[0], b1.to_reduced(seed.ravel()[b1.indices]))
        self.assertIsNone(dense_seeds[1])
        self.assertEqual(a.call_count, 1)
        np.testing.assert_array_equal(a.call_args.kwargs["seed_vec"],
                                      b0.to_reduced(seed.ravel()[b0.indices]))
        self.assertIs(a.call_args.args[0], ops[1].operator)
        self.assertIs(a.call_args.args[1], b0)
        self.assertFalse(a.call_args.kwargs["map_back"])
        # iteration: a zero seed part starts from the fixed random vector
        b1, b0 = _bases((5, 6))
        ops = [_sector(1, b1, _normal([1, 0.5, 0.2, -0.1, -0.3], 1)),
               _sector(0, b0, _normal([2.0, 0.5, 0.3, -0.1, -0.4, -0.6], 2))]
        seed = _seed().ravel().copy()
        seed[b0.indices] = 0.0
        with mock.patch.object(sc_so, "_solve_leading", wraps=sc_eigen._solve_leading) as spy:
            _iterate(ops, sigma_init=seed.reshape(GAP))
        inits = [c.kwargs["init_vec"] for c in spy.call_args_list]
        np.testing.assert_array_equal(inits[0], b1.to_reduced(seed[b1.indices]))
        rng = np.random.default_rng(0)
        np.testing.assert_array_equal(inits[1], rng.standard_normal(6) + 1j * rng.standard_normal(6))

    def test_round_off_seed_part_is_no_seed(self):
        b1, b0, bm = _bases((5, 6, 4))
        ops = [_sector(1, b1, _normal([1, 0.5, 0.2, -0.1, -0.3], 1)),
               _sector(0, b0, _normal([2.0, 0.5, 0.3, -0.1, -0.4, -0.6], 2)),
               _sector(-1, bm, _normal([0.4, 0.1, -0.2, -0.5], 3))]
        base = np.zeros(int(np.prod(GAP)), complex)
        base[b0.indices] = _seed().ravel()[b0.indices]
        noise = np.zeros_like(base)
        noise[b1.indices] = _seed().ravel()[b1.indices]
        noise *= np.linalg.norm(base) / np.linalg.norm(noise)
        real = sc_so._so_dense_reduced_eigen
        for scale, seeded_plus in ((1e-12, False), (1e-9, False), (1e-7, True)):
            with self.subTest(scale=scale):
                seed = base + scale * noise
                with mock.patch.object(sc_so, "_so_dense_reduced_eigen", wraps=real) as d, \
                        self.assertLogs(sc_static.logger, level="INFO") as cm:
                    sc_static.logger.info("start")
                    s = _solve(ops, k=4, seed_vec=seed.reshape(GAP))
                self.assertTrue(s.converged)
                seeds = [c.args[3] for c in d.call_args_list]
                self.assertEqual(len(seeds), 3)
                np.testing.assert_array_equal(seeds[1], b0.to_reduced(seed[b0.indices]))
                self.assertIsNone(seeds[2])                     # exactly zero: no INFO
                info = [m for m in _messages(cm) if "seed part" in m]
                if seeded_plus:
                    np.testing.assert_array_equal(seeds[0], b1.to_reduced(seed[b1.indices]))
                    self.assertEqual(info, [])
                else:
                    self.assertIsNone(seeds[0])
                    self.assertEqual(len(info), 1)
                    self.assertTrue(re.match(
                        r"^spin-orbital sector \+1: the seed part \(norm \d\.\d{3}e[+-]\d+, "
                        r"\d\.\d{3}e[+-]\d+ of the seed\) is round-off; the sector is solved "
                        r"unseeded\.$", info[0]), info[0])

    def test_eigen_leg_dim_one_non_finite(self):
        """A one-dimensional sector with a non-finite element fails its dense
        diagonalization (LinAlgError): the sector is failed and the run is
        not converged in "eigenvalue" and "both" mode."""
        b1, b0 = _bases((1, 4))
        for bad in (np.nan, np.inf, complex(1.0, np.inf)):
            with self.subTest(element=bad):
                ops = [_scalar_sector(1, b1, bad),
                       _sector(0, b0, _normal([0.7, 0.2, -0.1, -0.6], 4))]
                with self.assertLogs(sc_static.logger, level="WARNING") as cm, \
                        np.errstate(invalid="ignore"):
                    s = _solve(ops)
                self.assertTrue(any("dense diagonalization failed" in m for m in _messages(cm)),
                                _messages(cm))
                self.assertEqual(s.leaders[0].selection, "not-converged")
                self.assertEqual(s.failed_sectors, (1,))
                self.assertFalse(s.converged)
                self.assertTrue(np.isnan(s.leading))
                with np.errstate(invalid="ignore"), quiet():
                    it = _iterate(ops)
                for mode in ("eigenvalue", "both"):
                    r = _report(mode, s, it if mode == "both" else None)
                    self.assertTrue(np.isnan(r.reported), mode)
                    self.assertFalse(r.eigen_converged, mode)
                    self.assertEqual(r.failed_sectors, (1,), mode)

    def test_non_finite_seed_refused(self):
        b1, b0 = _bases((5, 6))
        ops = [_sector(1, b1, _normal([1, 0.5, 0.2, -0.1, -0.3], 1)),
               _sector(0, b0, _normal([2.0, 0.5, 0.3, -0.1, -0.4, -0.6], 2))]
        for bad in (np.nan, np.inf):
            seed = _seed().copy()
            seed.ravel()[b0.indices[0]] = bad
            with self.subTest(entry=bad), self.assertRaisesRegex(ValueError, exact(
                    "_so_solve_sectors: seed_vec has non-finite entries")):
                _solve(ops, k=4, seed_vec=seed)

    def test_seeded_union_orders_by_overlap(self):
        b1, b0 = _bases((5, 6))
        M1 = _normal([1.0, 0.5, 0.2, -0.1, -0.3], 1)
        M0 = _normal([2.0, 0.6, 0.3, -0.1, -0.4, -0.6], 2)
        ops = [_sector(1, b1, M1), _sector(0, b0, M0)]
        w1, v1 = np.linalg.eigh(M1)                 # normal and Hermitian: real spectrum
        w0, v0 = np.linalg.eigh(M0)
        a = _lift(b1, v1[:, np.argmin(abs(w1 - 0.5))])
        b = _lift(b0, v0[:, np.argmin(abs(w0 - 0.3))])
        seed = (0.8 * a / np.linalg.norm(a) + 0.6 * b / np.linalg.norm(b)).reshape(GAP)
        unseeded = _solve(ops, k=6)
        with self.assertNoLogs(sc_static.logger, level="WARNING"):
            s = _solve(ops, k=6, seed_vec=seed)
        self.assertTrue(s.converged)
        # the selection is unchanged; only the order follows the seed
        np.testing.assert_allclose(np.sort(s.vals.real), np.sort(unseeded.vals.real), atol=1e-12)
        np.testing.assert_allclose(s.vals[:2].real, [0.5, 0.3], atol=1e-12)
        np.testing.assert_array_equal(s.pp_sector[:2], np.array([1, 0], np.int8))
        self.assertAlmostEqual(s.leading.real, 0.5, delta=1e-12)
        self.assertAlmostEqual(unseeded.leading.real, 2.0, delta=1e-12)
        # a seed outside every selected branch: warned, still ordered by overlap
        with self.assertLogs(sc_static.logger, level="WARNING") as cm:
            _solve(ops, k=6, seed_vec=_seed())
        self.assertEqual(len(cm.records), 1)
        self.assertTrue(re.match(
            r"^spin-orbital: the largest overlap of a selected eigenvector with the seed is "
            r"0\.\d{3} \(< 0\.5\); the seed's branch is not among the selected "
            r"eigenpairs\.$", _messages(cm)[0]), _messages(cm)[0])

    def test_arpack_residual_failure(self):
        (b,) = _bases((30,))
        op = _sector(0, b, _normal(np.linspace(-1, 2, 30), 6))
        real = sc_so._so_solve_eigenvalue

        def perturbed(*a, **kw):
            vals, vecs = real(*a, **kw)
            vals = np.array(vals)
            vals[1] += 1e-6
            return vals, vecs
        with mock.patch.object(sc_so, "_so_solve_eigenvalue", perturbed), \
                self.assertLogs(sc_static.logger, level="WARNING") as cm:
            s = _solve([op], k=4)
        lead = s.leaders[0]
        self.assertEqual((lead.selection, lead.norm_kind, lead.certified),
                         ("residual-failed", "adjoint", False))
        self.assertGreater(lead.residual, 1e-8)
        self.assertEqual((s.failed_sectors, s.converged), ((0,), False))
        self.assertEqual(len(s.vals), 0)
        self.assertEqual(len(cm.records), 1)
        self.assertIn("fails the residual check", _messages(cm)[0])

    def test_adjoint_estimate_not_converged_is_logged(self):
        (b,) = _bases((30,))
        op = _sector(1, b, _normal(np.linspace(-1, 2, 30), 7))
        real = sc_so._so_sector_norm_estimate

        def unconverged(red, dim, adjoint):
            est, steps, change, _ = real(red, dim, adjoint)
            return est, steps, change, False
        with mock.patch.object(sc_so, "_so_sector_norm_estimate", unconverged), \
                self.assertLogs(sc_static.logger, level="INFO") as cm:
            s = _solve([op], k=4)
        self.assertTrue(s.converged)
        self.assertEqual(s.leaders[0].norm_kind, "adjoint")
        pat = (r"^spin-orbital sector \+1: operator-norm estimate adjoint after \d+ steps "
               r"\(relative change \d\.\d\de[+-]\d+\); residual scale \d\.\d{6}e[+-]\d+$")
        self.assertEqual(sum(bool(re.match(pat, m)) for m in _messages(cm)), 1)

    def test_zero_seed_part_on_arpack(self):
        b1, b0 = _bases((5, 30))
        ops = [_sector(1, b1, _normal([1, 0.5, 0.2, -0.1, -0.3], 1)),
               _sector(0, b0, _normal(np.linspace(-1, 2, 30), 2))]
        seed = _seed().ravel().copy()
        seed[b0.indices] = 0.0
        unseeded = _solve(ops, k=4)
        for s_vec in (seed, np.zeros_like(seed)):
            with mock.patch.object(sc_so, "_so_solve_eigenvalue",
                                   wraps=sc_so._so_solve_eigenvalue) as a, quiet():
                s = _solve(ops, k=4, seed_vec=s_vec.reshape(GAP))
            self.assertEqual(a.call_count, 1)
            self.assertIsNone(a.call_args.kwargs["seed_vec"])
            np.testing.assert_allclose(np.sort(s.vals.real), np.sort(unseeded.vals.real),
                                       atol=1e-12)
        # an all-zero seed is no seed: the unseeded order
        np.testing.assert_array_equal(s.pp_sector, unseeded.pp_sector)

    def test_basis_shape_is_the_gap_shape(self):
        (b,) = _bases((3,))
        with self.assertRaisesRegex(ValueError, exact(
                "_so_solve_sectors: the basis of sector +1 lives on a gap of shape "
                "{}, not {}".format(GAP, GAP + (1,)))):
            sc_so._so_solve_sectors([_sector(1, b, np.eye(3))], 2, "arnoldi",
                                 gap_shape=GAP + (1,))

    def test_norm_estimate_and_adjoint_contract(self):
        (b,) = _bases((30,))
        rng = np.random.default_rng(4)
        U, _ = np.linalg.qr(rng.normal(size=(30, 30)) + 1j * rng.normal(size=(30, 30)))
        V, _ = np.linalg.qr(rng.normal(size=(30, 30)) + 1j * rng.normal(size=(30, 30)))
        sv = np.concatenate([[5.0], np.linspace(1.0, 0.1, 29)])
        M = (U * sv) @ V.conj().T
        exact_norm = np.linalg.norm(M, 2)
        s = _solve([_sector(1, b, M)], k=4)
        lead = s.leaders[0]
        self.assertEqual(lead.norm_kind, "adjoint")
        self.assertTrue(lead.certified)
        self.assertLessEqual(abs(lead.norm_estimate - exact_norm), 1e-3 * exact_norm)
        with self.assertLogs(sc_static.logger, level="INFO") as cm:
            s = _solve([_sector(1, b, M, adjoint=False)], k=4)
        lead = s.leaders[0]
        self.assertEqual(lead.norm_kind, "forward")
        self.assertGreater(lead.norm_estimate, 0.0)
        self.assertLessEqual(lead.norm_estimate, exact_norm * (1 + 1e-12))
        pat = (r"^spin-orbital sector \+1: operator-norm estimate forward after \d+ steps "
               r"\(relative change \d\.\d\de[+-]\d+\); residual scale \d\.\d{6}e[+-]\d+$")
        self.assertEqual(sum(bool(re.match(pat, m)) for m in _messages(cm)), 1)
        # the default reduction carries no adjoint
        with self.assertRaises(NotImplementedError):
            b.reduce(_sector(1, b, M).operator).rmatvec(np.ones(30, complex))


class TestReportAdapter(unittest.TestCase):
    """Construction check: the dense legs pass the existing solvers' results
    through the one report."""

    def test_dense_leg_is_pass_through(self):
        basis = sc_eigen._SOAntisymmetricBasis(SHAPE)
        _, A = _random_operator(basis.n, 21)
        with self.assertLogs(sc_static.logger, level="WARNING") as cm:
            vals, gaps = sc_so._so_solve_eigenvalue(A, basis, 6, "arnoldi")
        self.assertEqual(len(cm.records), 1)
        self.assertIn("re-solving for the largest real part", _messages(cm)[0])
        leg = sc_so.EigenLeg.from_dense(vals, gaps, selection="LM")
        self.assertIs(leg.vals, vals)
        self.assertIs(leg.gaps, gaps)
        self.assertEqual(leg.failed_sectors, ())
        r = sc_so._so_bond_report("eigenvalue", leg, None, gap_shape=SHAPE)
        self.assertTrue(np.array_equal(r.vals, vals))
        self.assertTrue(np.array_equal(r.gap, gaps[0]))
        self.assertEqual((r.reported, r.eigen_converged, r.pp_sector),
                         (float(vals[0].real), True, 2))
        np.testing.assert_array_equal(r.pp_sectors, np.full(len(vals), 2, np.int8))
        self.assertEqual((r.gap_status, r.selection), ("certified", "LM"))

    def test_dense_matches_existing_contract(self):
        """Construction check: the documented rule of the dense spin-orbital
        run (reported value, convergence flags and written gap)."""
        g_eig = [np.full(SHAPE, 1.0 + 0j), np.full(SHAPE, 2.0 + 0j)]
        g_it = np.full(SHAPE, 3.0 + 0j)
        for eig_ok in (True, False):
            vals = np.array([0.4, 0.1]) + 0j if eig_ok else np.full(2, np.nan + 0j)
            eleg = sc_so.EigenLeg.from_dense(vals, np.array(g_eig), selection="LM")
            self.assertEqual(eleg.converged, eig_ok)
            for it_ok in (True, False):
                ileg = sc_so.IterationLeg.from_dense(g_it, 0.39 if it_ok else float("nan"), it_ok,
                                                  12, {"converged": it_ok})
                self.assertEqual((ileg.pp_sector, ileg.selection, ileg.failed_sectors),
                                 (2, "iteration", ()))
                for mode in ("eigenvalue", "iteration", "both"):
                    r = sc_so._so_bond_report(mode, eleg if mode != "iteration" else None,
                                           ileg if mode != "eigenvalue" else None,
                                           gap_shape=SHAPE)
                    if mode == "iteration":
                        want = (0.39 if it_ok else np.nan, it_ok, None, g_it)
                    elif mode == "eigenvalue":
                        # a NaN listing (failed-solve sentinel) is no candidate
                        want = (0.4 if eig_ok else np.nan, eig_ok, None,
                                g_eig[0] if eig_ok else np.zeros(SHAPE, complex))
                    else:
                        want = (0.4 if eig_ok else np.nan, eig_ok, it_ok,
                                g_eig[0] if eig_ok else g_it)
                    msg = (mode, eig_ok, it_ok)
                    np.testing.assert_equal(r.reported, want[0], err_msg=str(msg))
                    self.assertEqual((r.eigen_converged, r.iteration_converged), want[1:3], msg)
                    self.assertTrue(np.array_equal(r.gap, want[3]), msg)
                    self.assertEqual(r.pp_sector, 2, msg)
                    status = ("certified" if np.isfinite(want[0])
                              else "none" if mode == "eigenvalue" else "uncertified")
                    self.assertEqual(r.gap_status, status, msg)


if __name__ == "__main__":
    unittest.main()
