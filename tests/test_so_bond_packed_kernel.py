"""The bond pairing kernel on the pair-S_z sectors of the S_z-packed vertex.

``so_bond.PackedBondPairKernel`` applies the bond-resolved pairing kernel

``out_AD(k) = -(1/N) sum_{alpha, beta} e^{-ik.R_beta} sum_{k'}
V_{alpha beta; ABCD}(k - k') e^{-ik'.R_alpha} F_BC(k')``

(``F = G2 phi``) of a packed vertex (``PackedArrays`` or a
``PackedBondBlockStore`` slot) one pair-S_z sector ``p = S_z(B) + S_z(C)``
at a time. ``sector_op(p)`` returns the ``hwave.sc.SectorOp`` of sector
``p``: a forward and adjoint operator on the sector vector and the cached
fermion-antisymmetric basis of the sector's gap components.

Every result is compared with the dense production kernel
(``eliashberg_bond.BondPairKernel`` on ``V.to_dense()``), gathered at the
test-side pair lists (``tests.so_bond_oracles.pp_sector_lists``).

The eigenvalue set of the physical static kernel (the packed production
chain of ``tests.so_bond_fixtures.static_pair_kernels``) is compared with
the dense spectrum by the matching rule of ``tests.so_bond_fixtures``: an
optimal assignment with ``|lambda - mu| <= 1e-8 max(1, max|lambda|) max(1,
kappa_c)``, ``kappa_c`` the condition number of the eigenvalue's cluster
(fixtures with ``kappa_c > 1e4`` are rejected), and eigenvectors compared
as cluster projectors.
"""
import contextlib
import functools
import unittest
from unittest import mock

import numpy as np
from scipy.sparse.linalg import LinearOperator

from hwave import sc
from hwave.solver.eliashberg import so_solve as sc_so
from hwave.solver import eliashberg_dynamic
from hwave.solver.so_bond import (PackedArrays, PackedBondBlockStore, PackedBondPairKernel,
                                  estimate_packed_kernel_workspace_bytes, plan_so_bond)
from tests.so_bond_fixtures import (KAPPA_MAX, MATCH_TOL, NOCAP, assert_recomputed,
                                    cluster_projector, dense_spectrum, exact, match_spectrum,
                                    pp_flat_index, quiet, static_pair_kernels)
from tests.so_bond_fixtures import FIXTURE_NAMES as NAMES
from tests.so_bond_fixtures import cvec as _cvec
from tests.so_bond_fixtures import gap_shape as _gap_shape
from tests.so_bond_fixtures import kernel_setup as _setup
from tests.so_bond_fixtures import random_kernels as _kernels
from tests.so_bond_oracles import (FIXTURE_LZSZ, FIXTURE_LZSZ_GRID, FIXTURE_SZ_GRID,
                                   FIXTURE_SZ_GRID_H0, FIXTURE_SZ_GRID_SYM, FIXTURE_SZ_RING,
                                   PP_ORDER, pp_sector_lists)

TOL = 1e-13
#: (fixture, nfreq) of the forward and adjoint gates
CASES = ((FIXTURE_SZ_GRID, 1), (FIXTURE_LZSZ, 1), (FIXTURE_LZSZ, 4))


def _pp_index(fx, nfreq):
    """Per pp sector ``(+1, 0, -1)`` the ascending flat gap indices (test side)."""
    return pp_flat_index(fx.nd, fx.nvol * nfreq)


def _reverse(nfreq):
    return sc._reverse_k_and_orbital if nfreq == 1 else eliashberg_dynamic._reverse_kw_and_orbital


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / np.max(np.abs(b)))


class TestPackedKernel(unittest.TestCase):
    """Gate P-alg (pairing kernel): the packed kernel equals the dense
    production kernel on ``V.to_dense()`` (<= 1e-13 relative, max-norm) on
    the full gap and per pair-S_z sector; the dense kernel maps a sector
    onto itself (its output outside the sector <= 1e-13 of its maximum)."""

    def test_matvec_equals_dense(self):
        for fx, nfreq in CASES:
            packed, dense, _, _ = _kernels(fx, nfreq)
            shape = _gap_shape(fx, nfreq)
            full = sc._SOAntisymmetricBasis(shape, reverse=_reverse(nfreq))
            rng = np.random.default_rng(31)
            n = int(np.prod(shape))
            gaps = (_cvec(rng, n), full.to_full(_cvec(rng, full.dim)), _cvec(rng, n))
            for t, phi in enumerate(gaps):
                with self.subTest(fixture=NAMES[fx], nfreq=nfreq, gap=t):
                    want = dense.matvec(phi)
                    self.assertLessEqual(_rel(packed.matvec(phi), want), TOL)

    def test_sector_op_equals_dense_gathered(self):
        for fx, nfreq in CASES:
            packed, dense, _, _ = _kernels(fx, nfreq)
            n = int(np.prod(_gap_shape(fx, nfreq)))
            rng = np.random.default_rng(37)
            for p, idx in zip(PP_ORDER, _pp_index(fx, nfreq)):
                with self.subTest(fixture=NAMES[fx], nfreq=nfreq, p=p):
                    op = packed.sector_op(p)
                    x = _cvec(rng, len(idx))
                    phi = np.zeros(n, dtype=complex)
                    phi[idx] = x
                    out = dense.matvec(phi)
                    self.assertLessEqual(_rel(op.operator.matvec(x), out[idx]), TOL)
                    # conservation and non-vacuity of the split
                    scale = float(np.max(np.abs(out)))
                    self.assertGreater(float(np.max(np.abs(out[idx]))), 1e-3)
                    self.assertLessEqual(float(np.max(np.abs(np.delete(out, idx)))),
                                         TOL * scale)

    def test_store_slot_equals_packed_arrays(self):
        fx, nfreq = FIXTURE_LZSZ, 4
        geom, sm, plan = _setup(fx)
        packed, _, V, G2 = _kernels(fx, nfreq)
        with PackedBondBlockStore(nfreq, fx.nvol, sm, ("V",)) as store:
            store.put_freq_batch("V", 0, nfreq, V)
            k = PackedBondPairKernel(store, G2, geom, sm, plan=plan,
                                     spatial_shape=fx.spatial_shape, beta=fx.beta, nfreq=nfreq)
        phi = _cvec(np.random.default_rng(41), int(np.prod(_gap_shape(fx, nfreq))))
        assert_recomputed(self, k.matvec(phi), packed.matvec(phi))


class TestPackedKernelAdjoint(unittest.TestCase):
    """Gate P-alg (adjoint): ``|<y, K x> - <K^H y, x>| <= 1e-13 ||y|| ||K x||``
    per sector and on the full gap; on the ring the sector ``rmatvec`` equals
    the conjugate transpose of the dense kernel's sector matrix, built column
    by column (<= 1e-13 relative, max-norm)."""

    def test_adjoint_identity(self):
        for fx, nfreq in CASES:
            packed, _, _, _ = _kernels(fx, nfreq)
            rng = np.random.default_rng(43)
            for p in PP_ORDER:
                with self.subTest(fixture=NAMES[fx], nfreq=nfreq, p=p):
                    A = packed.sector_op(p).operator
                    x, y = _cvec(rng, A.shape[0]), _cvec(rng, A.shape[0])
                    Kx = A.matvec(x)
                    lhs, rhs = np.vdot(y, Kx), np.vdot(A.rmatvec(y), x)
                    self.assertLessEqual(abs(lhs - rhs),
                                         TOL * np.linalg.norm(y) * np.linalg.norm(Kx))
            n = int(np.prod(_gap_shape(fx, nfreq)))
            x, y = _cvec(rng, n), _cvec(rng, n)
            Kx = packed.matvec(x)
            self.assertLessEqual(abs(np.vdot(y, Kx) - np.vdot(packed.rmatvec(y), x)),
                                 TOL * np.linalg.norm(y) * np.linalg.norm(Kx))

    def test_rmatvec_is_dense_conjugate_transpose(self):
        fx = FIXTURE_SZ_RING
        for nfreq in (1, 4):
            packed, dense, _, _ = _kernels(fx, nfreq)
            n = int(np.prod(_gap_shape(fx, nfreq)))
            for p, idx in zip(PP_ORDER, _pp_index(fx, nfreq)):
                with self.subTest(nfreq=nfreq, p=p):
                    M = np.empty((len(idx), len(idx)), dtype=complex)
                    for c, j in enumerate(idx):
                        phi = np.zeros(n, dtype=complex)
                        phi[j] = 1.0
                        M[:, c] = dense.matvec(phi)[idx]
                    A = packed.sector_op(p).operator
                    R = np.array([A.rmatvec(e) for e in np.eye(len(idx))]).T
                    self.assertLessEqual(_rel(R, M.conj().T), TOL)


class TestPackedKernelSectors(unittest.TestCase):
    """Construction check: the sector API (indices, shapes, the cached
    closed bases) and the memory table."""

    def test_sector_bases(self):
        for fx, nfreq in ((FIXTURE_SZ_GRID, 1), (FIXTURE_LZSZ, 1), (FIXTURE_LZSZ, 4)):
            packed, _, _, _ = _kernels(fx, nfreq)
            shape = _gap_shape(fx, nfreq)
            self.assertEqual(packed.gap_shape, shape)
            self.assertEqual(packed.pp_sectors, (1, 0, -1))
            full = sc._SOAntisymmetricBasis(shape, reverse=_reverse(nfreq))
            dims = 0
            for p, idx, pairs in zip(PP_ORDER, _pp_index(fx, nfreq), pp_sector_lists(fx.nd)):
                with self.subTest(fixture=NAMES[fx], nfreq=nfreq, p=p):
                    np.testing.assert_array_equal(packed.pp_index(p), idx)
                    self.assertEqual(packed.sector_shape(p),
                                     (len(pairs),) + fx.spatial_shape
                                     + ((nfreq,) if nfreq > 1 else ()))
                    op = packed.sector_op(p)
                    self.assertIsInstance(op, sc.SectorOp)
                    self.assertEqual((op.pp_sector, op.has_adjoint), (p, True))
                    b = op.basis
                    np.testing.assert_array_equal(b.indices, idx)
                    self.assertEqual(b.shape, shape)
                    self.assertEqual(op.operator.shape, (len(idx), len(idx)))
                    # closed under the involution: the sector partner map is
                    # an involution and agrees with the full one
                    np.testing.assert_array_equal(b.partner[b.partner], np.arange(len(idx)))
                    np.testing.assert_array_equal(idx[b.partner], full.partner[idx])
                    dims += b.dim
                    again = packed.sector_op(p)
                    self.assertIs(again.basis, b)
            self.assertEqual(dims, full.dim)

    def test_bad_sector(self):
        packed = _kernels(FIXTURE_SZ_RING, 1)[0]
        for bad in (2, True, "0", 0.0):
            with self.assertRaisesRegex(ValueError, exact(
                    "PackedBondPairKernel: pp_sector must be +1, 0 or -1, got {!r}".format(bad))):
                packed.sector_op(bad)

    def test_memory_table(self):
        for fx, nfreq in CASES:
            packed = _kernels(fx, nfreq)[0]
            geom, sm, _ = _setup(fx)
            counts = [len(l) for l in pp_sector_lists(fx.nd)]
            want = estimate_packed_kernel_workspace_bytes(geom.B, counts, fx.nvol, nfreq)
            self.assertEqual(packed.memory_table(), want)
            # each packed vertex entry lies in exactly one pp sector: the
            # hoisted blocks hold the packed vertex once
            self.assertEqual(packed.hoisted_nbytes, want["hoisted_blocks"])
            self.assertEqual(want["hoisted_blocks"],
                             nfreq * fx.nvol * 16 * sum(n * n for n in sm.sizes))


class TestPackedKernelRefusals(unittest.TestCase):
    """Construction check: the refusals of the constructor."""

    def _args(self, fx=FIXTURE_SZ_RING, nfreq=1):
        geom, sm, plan = _setup(fx)
        _, _, V, G2 = _kernels(fx, nfreq)
        return dict(store=V, G2=G2, view=geom, sector_map=sm, plan=plan,
                    spatial_shape=fx.spatial_shape, beta=fx.beta, nfreq=nfreq)

    def _refused(self, msg, **change):
        kw = self._args()
        kw.update(change)
        with self.assertRaisesRegex(ValueError, exact(msg)):
            PackedBondPairKernel(kw.pop("store"), kw.pop("G2"), kw.pop("view"),
                                 kw.pop("sector_map"), **kw)

    def test_g2_pp_leak_refused(self):
        kw = self._args()
        G2 = kw["G2"].copy()
        G2[0, 0, 0, 1, 2, 0, 0, 0] = 0.5 - 0.25j      # pair (0, 0) from pair (0, 1)
        self._refused("PackedBondPairKernel: G2 couples pair S_z sectors (max |entry| = "
                      "{:.3e}); the packed kernel needs a Green function built per spin "
                      "sector".format(abs(0.5 - 0.25j)), G2=G2)

    def test_dense_plan_and_dimensions(self):
        fx = FIXTURE_SZ_RING
        geom, sm, plan = _setup(fx)
        dense_plan = plan_so_bond(nmat=2, nvol=fx.nvol, B=geom.B, nd=fx.nd,
                                  consumers=("vertex",), caps=NOCAP)
        self._refused("PackedBondPairKernel: the admission plan has layout = 'dense'; the "
                      "packed kernel needs 'sz_packed'", plan=dense_plan)
        self._refused("PackedBondPairKernel: the admission plan has nvol = {} but "
                      "spatial_shape (6, 1, 1) has 6 k points".format(fx.nvol),
                      spatial_shape=(6, 1, 1))
        other = plan_so_bond(nmat=2, nvol=fx.nvol, B=geom.B + 1, nd=fx.nd,
                             consumers=("vertex",), caps=NOCAP, layout="sz_packed",
                             layout_reason=plan.layout_reason)
        self._refused("PackedBondPairKernel: the admission plan has ND = {} but the sector map "
                      "has ND = {}".format(other.ND, sm.ND), plan=other)
        view = type("V", (), {"n_channels": geom.B - 1, "delta_r": geom.delta_r[:-1]})()
        self._refused("PackedBondPairKernel: view.n_channels = {} but the sector map has "
                      "B = {}".format(geom.B - 1, geom.B), view=view)

    def test_store_and_backend(self):
        fx = FIXTURE_SZ_RING
        geom, sm, _ = _setup(fx)
        V = self._args()["store"]
        self._refused("PackedBondPairKernel: the vertex must have leading shape (nfreq, nvol) = "
                      "(1, {}), got {}".format(fx.nvol, (2, fx.nvol)),
                      store=PackedArrays(sm, tuple(np.concatenate([b, b]) for b in V.blocks)))
        self._refused("PackedBondPairKernel: the vertex must be a PackedBondBlockStore or "
                      "PackedArrays, got ndarray", store=V.to_dense())
        self._refused("PackedBondPairKernel: G2 must have shape {}, got {}".format(
            (2,) * 4 + fx.spatial_shape + (1,), (2,) * 4 + fx.spatial_shape),
            G2=self._args()["G2"][..., 0])
        self._refused("PackedBondPairKernel: the S_z-packed kernel runs on the host (xp must be "
                      "numpy), got cupy", xp=type("M", (), {"__name__": "cupy"})())
        self._refused("PackedBondPairKernel: residency must be 'host', got 'device'",
                      residency="device")
        table = _kernels(fx, 1)[0].memory_table()
        need = sum(table.values())
        self._refused("PackedBondPairKernel: the kernel needs {} bytes (hoisted_blocks = {}, "
                      "g2_sectors = {}, hoist_workspace = {}, g2_gather = {}, "
                      "matvec_workspace = {}, call_overhead = {}) but host_cap is 1000 bytes"
                      .format(need, table["hoisted_blocks"], table["g2_sectors"],
                              table["hoist_workspace"], table["g2_gather"],
                              table["matvec_workspace"], table["call_overhead"]),
                      host_cap=1000)
        with PackedBondBlockStore(1, fx.nvol, sm, ("W",)) as store:
            with self.assertRaises(KeyError):
                kw = self._args()
                kw.pop("store")
                PackedBondPairKernel(store, kw.pop("G2"), kw.pop("view"), kw.pop("sector_map"),
                                     **kw)

    def test_host_cap_values(self):
        for bad in (float("nan"), float("inf"), -1.0, "lots"):
            with self.subTest(host_cap=bad):
                self._refused("PackedBondPairKernel: host_cap must be None (unlimited) or a "
                              "finite number of bytes >= 0, got {!r}".format(bad), host_cap=bad)
        # None is the only unlimited value; a finite cap above the table admits
        kw = self._args()
        need = sum(_kernels(FIXTURE_SZ_RING, 1)[0].memory_table().values())
        for cap in (None, need, float(need)):
            PackedBondPairKernel(kw["store"], kw["G2"], kw["view"], kw["sector_map"],
                                 plan=kw["plan"], spatial_shape=kw["spatial_shape"],
                                 beta=kw["beta"], nfreq=kw["nfreq"], host_cap=cap)

    def test_full_gap_size(self):
        packed = _kernels(FIXTURE_SZ_RING, 1)[0]
        n = int(np.prod(packed.gap_shape))
        with self.assertRaisesRegex(ValueError, exact(
                "PackedBondPairKernel.matvec: the gap must have {} entries, got {}".format(
                    n, n + 1))):
            packed.matvec(np.zeros(n + 1, dtype=complex))


# ---------------------------------------------------------------------------
# the eigenvalue set of the physical static kernel, packed vs dense
# ---------------------------------------------------------------------------

#: eigenpairs per solve (the sector dense branch holds dim <= max(2k + 1, 20))
K_EIG = 10
#: power iteration (both layouts): steps, mixing, stop on the iterate change
ITER = dict(max_iter=5000, alpha=0.5, tol=1.0e-10)
#: the iteration gap against the dense one (sine of the angle between them):
#: the iterate change stops at ITER["tol"], the distance to the eigenvector is
#: that over the spectral gap of the shifted iteration, < 1e3 tol here
ITER_GAP_TOL = 1.0e-7
#: the tie tolerance of the sector union (relative, floor 1)
TIE_TOL = 1.0e-12
#: the dense branch of a sector solve and the LR_retry tolerance of the
#: unseeded which='LM' set (both as in the solver)
DENSE_MAX = max(2 * K_EIG + 1, 20)
LM_TOL = 1.0e-10
#: the smallest relative gap that makes a selected eigenvalue set independent
#: of round-off (the K-th against the (K+1)-th value of the selection key)
SET_GAP = 1.0e-6


def _isometry(gap_shape, indices=None):
    """Test side: the isometry ``E`` ``(n_full, d)`` onto the fermion-
    antisymmetric static gaps ``Delta_AB(k) = -Delta_BA(-k)`` supported on
    the flat indices ``indices`` (all when None): one column ``(e_i - e_{R
    i}) / sqrt(2)`` per orbit ``{i, R i}``, ``i`` the smaller index, in
    ascending ``i``; ``R: (A, B, k) -> (B, A, -k)``."""
    shape = tuple(gap_shape)
    n = int(np.prod(shape))
    a, b, x, y, z = np.unravel_index(np.arange(n), shape)
    partner = np.ravel_multi_index((b, a, -x % shape[2], -y % shape[3], -z % shape[4]), shape)
    sel = np.arange(n) if indices is None else np.asarray(indices)
    assert np.all(np.isin(partner[sel], sel)), "the index set is not closed"
    reps = sel[sel < partner[sel]]
    E = np.zeros((n, len(reps)))
    E[reps, np.arange(len(reps))] = 2 ** -0.5
    E[partner[reps], np.arange(len(reps))] = -(2 ** -0.5)
    return E


def _columns(apply, E):
    return np.column_stack([apply(E[:, j].astype(complex)) for j in range(E.shape[1])])


class _EigenCase:
    """The physical static kernel of a fixture with its dense spectra
    (test side): ``E`` and the dense reduced matrix ``D = E^T K E`` of the
    production dense kernel, its :class:`DenseSpectrum` ``full``; per pp
    sector ``p`` the flat indices, ``E_p`` (full rows), ``T_p = E^T E_p``
    and the spectrum ``spec[p]`` of ``E_p^T K_p E_p`` with ``K_p`` the
    packed ``sector_op(p)``."""

    def __init__(self, fx):
        self.fx = fx
        self.packed, self.dense = static_pair_kernels(fx)
        self.shape = (fx.nd, fx.nd) + fx.spatial_shape
        self.E = _isometry(self.shape)
        self.full = dense_spectrum(self.E.T @ _columns(self.dense.matvec, self.E))
        self.idx, self.T, self.spec = {}, {}, {}
        for p, idx in zip(PP_ORDER, _pp_index(fx, 1)):
            Ep = _isometry(self.shape, idx)[idx]
            A = self.packed.sector_op(p).operator
            self.idx[p] = idx
            self.T[p] = self.E[idx].T @ Ep
            self.spec[p] = dense_spectrum(Ep.T @ _columns(A.matvec, Ep))

    def ops(self):
        return [self.packed.sector_op(p) for p in self.packed.pp_sectors]

    def dense_operator(self):
        n = int(np.prod(self.shape))
        return LinearOperator((n, n), matvec=lambda v: self.dense.matvec(np.asarray(v).ravel()),
                              dtype=complex)

    def dense_eigen(self):
        """``(vals, gaps, selection)`` of the dense leg: the production
        ``_so_solve_eigenvalue`` on the dense kernel, unseeded (cached by
        the fixture)."""
        return _dense_eigen(self.fx)

    def union(self):
        """``(mu, tags)``: the sector eigenvalues in the order ``(+1, 0,
        -1)`` and ``(p, index in the sector)`` of each."""
        mu = np.concatenate([self.spec[p].vals for p in PP_ORDER])
        tags = [(p, j) for p in PP_ORDER for j in range(len(self.spec[p].vals))]
        return mu, tags

    def leader(self):
        """``(index in full, its sector)`` of the largest real part of the
        dense spectrum; the sector holds the eigenvalue matched to it."""
        i = int(np.argmax(self.full.vals.real))
        mu, tags = self.union()
        cols, _ = match_spectrum(self.full, mu)
        return i, tags[int(np.flatnonzero(cols == i)[0])][0]

    def lift(self, x):
        """Full reduced coordinates -> the flat gap."""
        return self.E @ x


@functools.lru_cache(maxsize=None)
def _eigen_case(fx):
    return _EigenCase(fx)


@functools.lru_cache(maxsize=None)
def _dense_eigen(fx):
    case = _eigen_case(fx)
    info = {}
    with quiet():
        vals, gaps = sc._so_solve_eigenvalue(case.dense_operator(),
                                             sc._SOAntisymmetricBasis(case.shape), K_EIG,
                                             "arnoldi", info_out=info)
    return vals, gaps, info["selection"]


def _sine(g, h):
    """The sine of the angle between two gaps, ``||g' - h' (h'^H g')||_2``
    of the normalized ``g'``, ``h'`` (the component of ``g'`` orthogonal to
    ``h'``). Unlike ``sqrt(1 - cos^2)`` it has no ``sqrt(eps)`` floor: two
    gaps equal to round-off give a round-off sine."""
    g, h = np.ravel(g), np.ravel(h)
    g = g / np.linalg.norm(g)
    h = h / np.linalg.norm(h)
    return float(np.linalg.norm(g - h * np.vdot(h, g)))


def _ties(vals, i, j):
    return abs(vals[i] - vals[j]) <= TIE_TOL * max(1.0, abs(vals[i]), abs(vals[j]))


class TestEigenSet(unittest.TestCase):
    """Gate P-alg (eigenvalue set) on the 5x3 grid (B = 5, norb_phys = 1)
    in the Zeeman field 0.15 and at zero field with spin-symmetric hopping:
    the static kernel of the packed production chain splits into the pair-S_z
    sectors (antisymmetric dimensions 7 / 15 / 7); the union of the dense
    sector spectra equals the dense reduced spectrum of the production
    ``BondPairKernel`` by the matching rule, the eigenvectors agree as
    cluster projectors (<= 1e-8 max(1, kappa_c)); the per-sector solve
    (dense branch in every sector) lists the top-10 eigenvalues by real
    part; a degenerate +1 / -1 pair is listed +1 first; the field splits
    the two sectors' leaders by >= 1e-3; and the reports of the packed and
    the dense legs agree in the three solver modes."""

    FIXTURES = (FIXTURE_SZ_GRID, FIXTURE_SZ_GRID_H0)

    def test_sector_dimensions(self):
        for fx in self.FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                case = _eigen_case(fx)
                full = sc._SOAntisymmetricBasis(case.shape)
                dims = [case.packed.sector_op(p).basis.dim for p in PP_ORDER]
                self.assertEqual(dims, [case.T[p].shape[1] for p in PP_ORDER])
                self.assertEqual(dims, [7, 15, 7])
                self.assertEqual(sum(dims), full.dim)
                self.assertEqual(case.E.shape[1], full.dim)

    def test_union_matches_dense_spectrum(self):
        for fx in self.FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                case = _eigen_case(fx)
                full = case.full
                self.assertLessEqual(full.kappa.max(), KAPPA_MAX)
                for p in PP_ORDER:
                    self.assertLessEqual(case.spec[p].kappa.max(), KAPPA_MAX)
                mu, tags = case.union()
                self.assertEqual(len(mu), len(full.vals))
                cols, ratio = match_spectrum(full, mu)
                self.assertLessEqual(ratio, 1.0)
                for c in full.clusters:
                    Pd = cluster_projector(full.right, full.left, c)
                    Pp = np.zeros_like(Pd)
                    for p in PP_ORDER:
                        J = [j for (q, j), col in zip(tags, cols) if q == p and col in c]
                        if J:
                            sp = case.spec[p]
                            Pp += case.T[p] @ cluster_projector(sp.right, sp.left, J) @ case.T[p].T
                    kappa = full.kappa[c[0]]
                    self.assertLessEqual(np.linalg.norm(Pp - Pd, 2),
                                         MATCH_TOL * max(1.0, kappa), c)

    def _solve(self, case, **kw):
        with quiet():
            return sc._so_solve_sectors(case.ops(), K_EIG, "arnoldi", gap_shape=case.shape, **kw)

    def test_sector_solve_lists_the_top_real_parts(self):
        for fx in self.FIXTURES:
            with self.subTest(fixture=NAMES[fx]):
                case = _eigen_case(fx)
                full = case.full
                s = self._solve(case)
                self.assertEqual(s.selection, "+1:dense-LR;0:dense-LR;-1:dense-LR")
                self.assertEqual([(l.pp_sector, l.dim) for l in s.leaders],
                                 [(1, 7), (0, 15), (-1, 7)])
                self.assertTrue(s.converged)
                self.assertTrue(all(l.certified and l.norm_kind == "exact" for l in s.leaders))
                self.assertEqual(len(s.vals), K_EIG)
                _, ratio = match_spectrum(full, s.vals)
                self.assertLessEqual(ratio, 1.0)
                bound = MATCH_TOL * full.scale * max(1.0, full.kappa.max())
                want = np.sort(full.vals.real)[::-1][:K_EIG]
                np.testing.assert_allclose(s.vals.real, want, rtol=0, atol=bound)
                # the leader: its sector and its gap (a simple eigenvalue)
                i, p_lead = case.leader()
                c = next(c for c in full.clusters if i in c)
                self.assertEqual(len(c), 1)
                self.assertEqual(int(s.pp_sector[0]), p_lead)
                self.assertLessEqual(_sine(s.gaps[0], case.lift(full.right[:, i])),
                                     MATCH_TOL * max(1.0, full.kappa[i]))

    def test_degenerate_pairs_plus_first(self):
        case = _eigen_case(FIXTURE_SZ_GRID_H0)
        s = self._solve(case)
        lead = {l.pp_sector: l.value for l in s.leaders}
        self.assertTrue(_ties([lead[1], lead[-1]], 0, 1), lead)
        pairs = 0
        for i in range(len(s.vals)):
            for j in range(i + 1, len(s.vals)):
                if _ties(s.vals, i, j) and {int(s.pp_sector[i]), int(s.pp_sector[j])} == {1, -1}:
                    pairs += 1
                    self.assertEqual(int(s.pp_sector[i]), 1, (i, j, s.vals[i], s.vals[j]))
        self.assertGreaterEqual(pairs, 2)

    def test_field_splits_the_sectors(self):
        for fx in (FIXTURE_SZ_GRID, FIXTURE_SZ_GRID_SYM):
            with self.subTest(fixture=NAMES[fx]):
                lead = {l.pp_sector: l.value for l in self._solve(_eigen_case(fx)).leaders}
                self.assertGreaterEqual(abs(lead[1] - lead[-1]), 1e-3)

    def test_reports_agree(self):
        case = _eigen_case(FIXTURE_SZ_GRID)
        full = case.full
        i_lead, p_lead = case.leader()
        lam = full.vals[i_lead]
        tol = MATCH_TOL * full.scale * max(1.0, full.kappa[i_lead])
        A = case.dense_operator()
        basis = sc._SOAntisymmetricBasis(case.shape)
        rng = np.random.default_rng(47)
        init = case.lift(rng.normal(size=basis.dim) + 1j * rng.normal(size=basis.dim))
        init = init.reshape(case.shape)
        vals, gaps, selection = case.dense_eigen()
        with quiet():
            gap, value, conv, n_iter = sc._so_solve_iteration(A, basis, init, ITER["max_iter"],
                                                              ITER["alpha"], ITER["tol"],
                                                              spectral_shift="auto")
            it_p = sc._so_iterate_sectors(case.ops(), init, gap_shape=case.shape,
                                          spectral_shift="auto", **ITER)
        legs = {"dense": (sc.EigenLeg.from_dense(vals, gaps, selection=selection),
                          sc.IterationLeg.from_dense(gap, value, conv, n_iter, {})),
                "packed": (sc.EigenLeg.from_sectors(self._solve(case)),
                           sc.IterationLeg.from_sectors(it_p))}
        for mode in ("eigenvalue", "iteration", "both"):
            reports = {}
            for layout, (eig, it) in legs.items():
                reports[layout] = sc._so_bond_report(
                    mode, None if mode == "iteration" else eig,
                    None if mode == "eigenvalue" else it, gap_shape=case.shape)
            rd, rp = reports["dense"], reports["packed"]
            with self.subTest(mode=mode):
                self.assertEqual((rd.gap_status, rp.gap_status), ("certified", "certified"))
                self.assertEqual((rd.pp_sector, rp.pp_sector), (2, p_lead))
                self.assertLessEqual(abs(rp.reported - lam.real), tol)
                self.assertLessEqual(abs(rd.reported - lam.real), tol)
                gap_tol = (ITER_GAP_TOL if mode == "iteration"
                           else MATCH_TOL * max(1.0, full.kappa[i_lead]))
                self.assertLessEqual(_sine(rp.gap, rd.gap), gap_tol)
                self.assertLessEqual(_sine(rp.gap, case.lift(full.right[:, i_lead])), gap_tol)
                if mode != "iteration":
                    for r in (rd, rp):
                        self.assertEqual(len(r.vals), K_EIG)
                        self.assertLessEqual(match_spectrum(full, r.vals)[1], 1.0)
                self.assertEqual((rd.eigen_converged, rp.eigen_converged), (True, True))
                if mode == "both":
                    self.assertEqual((rd.iteration_converged, rp.iteration_converged),
                                     (True, True))
        self.assertTrue(it_p.converged)
        self.assertLessEqual(abs(np.real(it_p.value) - np.real(value)), tol)


@contextlib.contextmanager
def _spy_sector_solves():
    """Record every ``sc._so_solve_eigenvalue`` call: the basis, the seed it
    received (reduced coordinates or None), the returned eigenvalues and the
    selection."""
    calls = []
    real = sc_so._so_solve_eigenvalue

    def spy(A, basis, num, method, **kw):
        vals, vecs = real(A, basis, num, method, **kw)
        info = kw.get("info_out")
        calls.append({"basis": basis, "seed": kw.get("seed_vec"), "vals": np.array(vals),
                      "selection": None if info is None else info.get("selection")})
        return vals, vecs
    with mock.patch.object(sc_so, "_so_solve_eigenvalue", spy):
        yield calls


def _top(values, key, k):
    """Indices of the ``k`` largest ``key(values)`` and the relative gap of
    the k-th to the (k+1)-th key."""
    keys = key(values)
    order = np.argsort(-keys, kind="stable")
    gap = (keys[order[k - 1]] - keys[order[k]]) / max(1.0, float(np.abs(values).max()))
    return set(order[:k].tolist()), float(gap)


class TestEigenSetArpack(unittest.TestCase):
    """Gate P-alg (eigenvalue set, ARPACK) on the 5x3 grid with norb_phys
    = 2 (B = 5; Zeeman 0.15 plus the ``i lambda s_z`` inter-orbital hopping;
    pp dimensions 29 / 60 / 29, all above the dense branch at k = 10) and a
    repulsive-dominant interaction: every sector runs on ARPACK -- unseeded
    ``which='LM'`` where its set holds the largest real part, the
    largest-real re-solve (``"LR_retry"``) in the sector where it cannot --
    with the adjoint-based residual check; the per-sector eigenvalues equal
    the selected sets of the dense sector spectra and the listed leaders the
    dense reduced spectrum of the full kernel (matching rule); a seed equal
    to the dense leading gap is split per sector into exactly what each
    sector solve receives and reproduces the leader."""

    @classmethod
    def setUpClass(cls):
        cls.case = case = _eigen_case(FIXTURE_LZSZ_GRID)
        # the selection each unseeded which='LM' sector solve must make
        cls.expected = {}
        for p in PP_ORDER:
            v = case.spec[p].vals
            lm, _ = _top(v, np.abs, K_EIG)
            top = v[sorted(lm)]
            retry = top.real.max() < np.abs(top).min() - LM_TOL * max(1.0, np.abs(top).max())
            cls.expected[p] = "LR_retry" if retry else "LM"

    def test_fixture_preconditions(self):
        case = self.case
        self.assertEqual([case.T[p].shape[1] for p in PP_ORDER], [29, 60, 29])
        self.assertTrue(all(case.T[p].shape[1] > DENSE_MAX for p in PP_ORDER))
        self.assertLessEqual(case.full.kappa.max(), KAPPA_MAX)
        negative = []
        for p in PP_ORDER:
            sp = case.spec[p]
            self.assertLessEqual(sp.kappa.max(), KAPPA_MAX)
            v = sp.vals
            if v[np.argmax(np.abs(v))].real < 0:
                negative.append(p)
            # the selected sets are independent of round-off
            lm, gap_lm = _top(v, np.abs, K_EIG)
            self.assertGreaterEqual(gap_lm, SET_GAP, p)
            top = v[sorted(lm)]
            margin = (np.abs(top).min() - top.real.max()) / max(1.0, np.abs(top).max())
            self.assertGreaterEqual(abs(margin), SET_GAP, p)
            if self.expected[p] == "LR_retry":
                self.assertGreaterEqual(_top(v, np.real, K_EIG)[1], SET_GAP, p)
        self.assertTrue(negative)                       # repulsive-dominant
        self.assertEqual([self.expected[p] for p in PP_ORDER], ["LM", "LR_retry", "LM"])

    def _sector_of(self, basis):
        return next(p for p in PP_ORDER if self.case.packed.sector_op(p).basis is basis)

    def test_every_sector_on_arpack(self):
        case = self.case
        with self.assertLogs("hwave_sc", level="WARNING") as logs, \
                _spy_sector_solves() as calls:
            s = sc._so_solve_sectors(case.ops(), K_EIG, "arnoldi", gap_shape=case.shape)
        retries = [m for m in logs.output if "re-solving for the largest real part" in m]
        self.assertEqual(len(retries), 1, logs.output)
        self.assertTrue(s.converged)
        self.assertEqual([l.selection for l in s.leaders],
                         [self.expected[p] for p in PP_ORDER])
        for lead in s.leaders:
            self.assertEqual((lead.norm_kind, lead.certified), ("adjoint", True))
            self.assertLessEqual(lead.residual, 1e-10)
            self.assertEqual(lead.n_returned, K_EIG)
        self.assertEqual(len(calls), 3)
        for call in calls:
            p = self._sector_of(call["basis"])
            with self.subTest(p=p):
                self.assertIsNone(call["seed"])
                self.assertEqual(call["selection"], self.expected[p])
                sp = case.spec[p]
                key = np.abs if self.expected[p] == "LM" else np.real
                want, _ = _top(sp.vals, key, K_EIG)
                cols, ratio = match_spectrum(sp, call["vals"])
                self.assertLessEqual(ratio, 1.0)
                self.assertEqual(set(cols.tolist()), want)

    def _expected_listing(self):
        """The top-K by real part of the union of the selected sector sets."""
        vals = []
        for p in PP_ORDER:
            v = self.case.spec[p].vals
            key = np.abs if self.expected[p] == "LM" else np.real
            vals.extend(v[sorted(_top(v, key, K_EIG)[0])])
        vals = np.array(vals)
        return vals[np.argsort(-vals.real, kind="stable")[:K_EIG]]

    def test_leaders_match_the_full_spectrum(self):
        case = self.case
        full = case.full
        with quiet():
            s = sc._so_solve_sectors(case.ops(), K_EIG, "arnoldi", gap_shape=case.shape)
        vals, gaps, selection = case.dense_eigen()
        _, ratio = match_spectrum(full, s.vals)
        self.assertLessEqual(ratio, 1.0)
        bound = MATCH_TOL * full.scale * max(1.0, full.kappa.max())
        np.testing.assert_allclose(s.vals.real, self._expected_listing().real, rtol=0,
                                   atol=bound)
        i_lead, p_lead = case.leader()
        tol = MATCH_TOL * full.scale * max(1.0, full.kappa[i_lead])
        self.assertLessEqual(abs(s.leading - full.vals[i_lead]), tol)
        self.assertEqual(int(s.pp_sector[0]), p_lead)
        self.assertLessEqual(_sine(s.gaps[0], case.lift(full.right[:, i_lead])),
                             MATCH_TOL * max(1.0, full.kappa[i_lead]))
        # the reports of the two layouts
        rd = sc._so_bond_report("eigenvalue", sc.EigenLeg.from_dense(
            vals, gaps, selection=selection), None, gap_shape=case.shape)
        rp = sc._so_bond_report("eigenvalue", sc.EigenLeg.from_sectors(s), None,
                                gap_shape=case.shape)
        self.assertEqual((rd.gap_status, rp.gap_status), ("certified", "certified"))
        self.assertEqual((rd.pp_sector, rp.pp_sector), (2, p_lead))
        self.assertLessEqual(abs(rp.reported - rd.reported), tol)
        self.assertLessEqual(_sine(rp.gap, rd.gap), MATCH_TOL * max(1.0, full.kappa[i_lead]))

    def test_seed_split(self):
        case = self.case
        full = case.full
        i_lead, p_lead = case.leader()
        vals, gaps, _ = case.dense_eigen()
        with quiet():
            # a seeded which='LM' run tracks its branch and is never re-solved
            # for the largest real part, so the seeded solve takes the shift
            with _spy_sector_solves() as calls:
                s = sc._so_solve_sectors(case.ops(), K_EIG, "arnoldi", gap_shape=case.shape,
                                         spectral_shift="auto", seed_vec=gaps[0])
        seed = gaps[0].ravel()
        self.assertEqual(s.selection, "+1:LR;0:LR;-1:LR")
        self.assertEqual(len(calls), 3)
        for call in calls:
            p = self._sector_of(call["basis"])
            with self.subTest(p=p):
                Ep = _isometry(case.shape, case.idx[p])
                part = Ep.T @ seed
                if p == p_lead:
                    self.assertLessEqual(np.abs(call["seed"] - part).max(),
                                         1e-14 * np.abs(seed).max())
                else:
                    # the leading gap lies in its own sector: the other parts
                    # are round-off of the dense solve and are no seed
                    self.assertLessEqual(np.linalg.norm(part), 1e-8 * np.linalg.norm(seed))
                    self.assertIsNone(call["seed"])
        tol = MATCH_TOL * full.scale * max(1.0, full.kappa[i_lead])
        self.assertLessEqual(abs(s.leading - full.vals[i_lead]), tol)
        self.assertLessEqual(abs(s.leading - vals[0]), tol)
        self.assertEqual(int(s.pp_sector[0]), p_lead)
        self.assertLessEqual(_sine(s.gaps[0], gaps[0]), MATCH_TOL * max(1.0, full.kappa[i_lead]))


if __name__ == "__main__":
    unittest.main()
