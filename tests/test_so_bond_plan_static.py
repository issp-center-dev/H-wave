"""The admission rows of the static spin-orbital bond path.

* :data:`tests.so_bond_oracles.FIXTURE_SZ_GRID_4X2`: a 4x2 S_z-conserving
  mesh with interactions along x only (B = 3 channels: on site, +-x) and
  even extents on both axes, so ``N_TRIM = 4`` antisymmetry-fixed points.
* ``so_bond.n_trim`` and ``so_bond.pp_sector_dims``: the fixed points of
  ``k -> -k`` and the pair-S_z sector sizes of a static gap, against the
  antisymmetric basis and the packed kernel's sectors.
* ``so_bond.static_so_bond_rows`` and its estimators
  (``estimate_dense_layout_solver_bytes``, ``estimate_sector_iteration_bytes``,
  ``estimate_sector_solver_bytes(..., dense_fallback=True)``): every row
  equals its written-out formula; the eigensolver row follows the branch
  the solver takes and reserves the failure-only dense fallback of an
  ARPACK branch (``dim <= DENSE_FALLBACK_MAX = 2000``).
* ``plan_so_bond(..., extra_rows=...)``: the peak is the maximum of
  ``1.25 SOBondPlan.phase_total(P)`` over the consumers and every other
  phase of the rows (``SOBondPlan.peak_phase``); the plans without extra
  rows are unchanged. ``SOBondPlan.require("so_bond_static_kernel", ...)``
  checks the gap length ``n_gap = nd^2 nvol``.
* The traced peaks (``tracemalloc``) of the solver branches, the dense
  fallback after an immediate and after a late ARPACK failure included, are
  within their rows. ``tracemalloc`` sees the Python allocations (numpy
  arrays, the ARPACK work arrays scipy allocates); the allocations inside
  ``np.linalg.eig`` (its input copy, LAPACK's eigenvectors and eigenvalues,
  the ``zgeev`` workspace) are native. At the fallback the two moments are
  kept apart: while ``np.linalg.eig`` runs, the traced bytes plus the host
  LAPACK's native bytes (``zgeev`` workspace query) are at most the row;
  after it returns (the native buffers freed, a third ``dim x dim`` array --
  the reordered eigenvectors -- alive), the traced peak is at most the row
  minus its native share plus ``16 dim^2``. The native share covers the host
  LAPACK's sizes arithmetically.
"""
import logging
import tracemalloc
import unittest
from unittest import mock

import numpy as np
import scipy.linalg.lapack as lapack
from scipy.sparse.linalg import ArpackNoConvergence, LinearOperator
from scipy.sparse.linalg import eigs as real_eigs

from hwave import sc
from hwave.solver import so_bond, spin_orbital
from hwave.solver.so_bond import SOBondCaps, plan_so_bond, slot_geometry, so_bond_topology
from hwave.solver.so_records import compile_interaction
from tests import so_bond_entry_fixtures as F
from tests.so_bond_fixtures import (FIXTURE_NAMES, NOCAP, chain_plan, exact, quiet,
                                    random_kernels, sector_op_of, subset_bases)
from tests.so_bond_oracles import (CELL, FIXTURE_LZSZ_GRID, FIXTURE_SZ_GRID, FIXTURE_SZ_GRID_4X2,
                                   PP_ORDER)

logger = logging.getLogger(__name__)


def _n_trim(shape):
    """``prod_i (2 if N_i is even else 1)`` (test side)."""
    return int(np.prod([2 if n % 2 == 0 else 1 for n in shape]))


class TestFixture4x2(unittest.TestCase):
    """Construction check of the 4x2 fixture: three channels, accepted by
    the slot geometry on its own mesh, four fixed points."""

    def test_fixture_4x2_topology(self):
        fx = FIXTURE_SZ_GRID_4X2
        self.assertEqual(fx.spatial_shape, (4, 2, 1))
        self.assertEqual(fx.norb_phys, 1)
        topo = so_bond_topology(compile_interaction(fx.tables_dict, fx.norb_phys), CELL)
        self.assertEqual(len(topo.delta_r), 3)
        self.assertEqual(sorted(tuple(int(x) for x in d) for d in topo.delta_r),
                         [(-1, 0, 0), (0, 0, 0), (1, 0, 0)])
        geom = slot_geometry(topo, fx.spatial_shape)
        self.assertEqual(geom.B, 3)
        self.assertEqual(geom.nvol, 8)
        self.assertEqual(_n_trim(fx.spatial_shape), 4)


# ---------------------------------------------------------------------------
# the literal formulas (test side)
# ---------------------------------------------------------------------------

OV = 16 * 1024


def _F(d):
    """The failure-only dense fallback of an ARPACK solve (dim <= 2000):
    four dim^2 arrays and two vectors, the zgeev workspace (max(130 dim, 4352)
    complex entries) and rwork (2 dim reals)."""
    return 16 * (4 * d * d + 2 * d) + 16 * max(130 * d, 4352) + 8 * 2 * d if d <= 2000 else 0


def _native_share(d):
    """The native part of ``_F(d)`` (inside ``np.linalg.eig``): the input
    copy, LAPACK's eigenvectors and eigenvalues, ``work`` and ``rwork``."""
    return 16 * (2 * d * d + d) + 16 * max(130 * d, 4352) + 8 * 2 * d


def _host_native(d):
    """The native bytes of ``np.linalg.eig`` on a complex ``d x d`` matrix with
    this host's LAPACK (its optimal ``zgeev`` workspace)."""
    lwork = int(lapack.zgeev_lwork(d, compute_vl=0, compute_vr=1)[0].real)
    return 16 * (2 * d * d + d + lwork) + 8 * 2 * d


def _retained(d, k, kind):
    """The outputs kept into the largest-real solve: none for a numeric
    shift, the ``min(6, d - 2)`` eigenpairs of the spectral-radius solve for
    "auto", and those plus the ``k`` first which='LM' eigenvalues unshifted."""
    radius_set = 16 * min(6, d - 2) * (d + 1)
    return {"numeric": 0, "auto": radius_set, "none": radius_set + 16 * k}[kind]


def _arpack(d, k, n, kind="none"):
    """The ARPACK branch for the spectral_shift kind ``kind``, with a second
    Krylov working set above d = 2000 (the placeholder which='LM' solve after
    a late failure, while the failed solve's arrays are alive)."""
    ncv = min(max(2 * k + 1, 20), d)
    work = 16 * (ncv * d + 3 * ncv * ncv + 14 * ncv + 3 * d)
    return ((2 if d > 2000 else 1) * work + 24 * n + 6 * 16 * n
            + 16 * 3 * k * (d + n) + _retained(d, k, kind))


def _direct(d, n):
    return 16 * (7 * d * d + 3 * n * d + 4 * d)


def _dense_solver(d, k, n, mode, max_iter):
    eig = ((_direct(d, n) if k >= d - 1 else _arpack(d, k, n) + _F(d))
           + 16 * k * n + 16 * k * d + 24 * n + OV)
    it = 16 * (6 * d + 4 * n) + 24 * n + 80 * max_iter + OV
    return {"eigenvalue": eig, "iteration": it, "both": eig + it}[mode]


def _sector_solver(dims, k, n_full, nvec, fallback):
    per = []
    for d, n in zip(dims, nvec):
        if d <= max(2 * k + 1, 20):
            branch = _direct(d, n)
        else:
            branch = _arpack(d, k, n) + (_F(d) if fallback else 0)
        per.append(max(branch, 16 * (min(k, d) * (d + n) + n_full)))
    return (max(per) + 16 * (k + sum(min(k, d) for d in dims)) * n_full
            + 16 * sum(min(k, d) * d for d in dims) + OV)


def _sector_iteration(dims, n_full, nvec, max_iter):
    return (max(16 * (6 * d + 3 * n) for d, n in zip(dims, nvec)) + 16 * 4 * n_full
            + 80 * max_iter + OV)


def _sectors(nd, shape):
    """``(pair_counts, n_sector_vec, dims)`` of the pp sectors (+1, 0, -1)."""
    norb, nvol, trim = nd // 2, int(np.prod(shape)), _n_trim(shape)
    counts = (norb ** 2, 2 * norb ** 2, norb ** 2)
    equal = (norb ** 2 * nvol - norb * trim) // 2
    return counts, tuple(c * nvol for c in counts), (equal, norb ** 2 * nvol, equal)


def _green_transient(nmat, nvol, nd):
    """The build temporaries (two k-resolved nd x nd arrays and one
    Green-sized temporary below numpy's 256 KiB in-place threshold) or the
    fp_green digest (three gathered blocks of at most 16 frequencies or 64
    MiB, plus their finiteness mask, beyond the G of headroom over the
    resident carrier), whichever is larger."""
    G, sl, chunk = 16 * nmat * nvol * nd * nd, 16 * nvol * nd * nd, 64 * 2 ** 20
    if sl <= chunk:
        block = sl * min(16, nmat, max(1, chunk // sl))
    else:
        block = 16 * nd * nd * max(1, chunk // (16 * nd * nd))
    return max(2 * sl + min(G, 256 * 1024), 3 * block + block // 8 - G)


def _literal_rows(nmat, shape, B, nd, layout, mode, k, max_iter):
    nvol = int(np.prod(shape))
    norb = nd // 2
    G, T4, n, ND = 16 * nmat * nvol * nd * nd, 16 * nvol * nd ** 4, nd * nd * nvol, B * nd * nd
    S = ND * ND if layout == "dense" else (2 * B * norb ** 2) ** 2 + 2 * (B * norb ** 2) ** 2
    counts, nvec, dims = _sectors(nd, shape)
    if layout == "dense":
        work = 16 * nvol * ND * ND + 3 * T4 + (3 + B) * 16 * n + 16 * nvol * nd * nd
        solver = _dense_solver(sum(dims), k, n, mode, max_iter)
    else:
        sq, m = sum(c * c for c in counts), max(counts)
        work = (16 * B * B * sq * nvol + 16 * sq * nvol + 16 * 6 * m * m * nvol
                + 16 * m * m * nvol + 16 * (max(3 + B, 6) * m * nvol + 3 * 8192) + OV)
        solver = ((_sector_solver(dims, k, n, nvec, True) if mode != "iteration" else 0)
                  + (_sector_iteration(dims, n, nvec, max_iter) if mode != "eigenvalue" else 0))
    return (("green_static", "carrier_build", 4 * G), ("green_static", "call_overhead", OV),
            ("green_static", "eigendecomposition", 16 * nvol * nd * nd + 8 * nvol * nd),
            ("green_static", "transient", _green_transient(nmat, nvol, nd)),
            ("pairing_static", "carrier", 3 * G),
            ("g2_static", "carrier", 3 * G), ("g2_static", "v_static", 16 * nvol * S),
            ("g2_static", "g2_roll", G), ("g2_static", "g2_out", 2 * T4),
            ("g2_static", "call_overhead", OV),
            ("kernel_static", "carrier", 3 * G), ("kernel_static", "v_static", 16 * nvol * S),
            ("kernel_static", "g2", T4), ("kernel_static", "kernel_workspace", work),
            ("kernel_static", "eigensolver", solver))


def _rows(nmat, shape, B, nd, layout, mode="eigenvalue", k=10, max_iter=1000):
    return so_bond.static_so_bond_rows(
        nmat=nmat, nvol=int(np.prod(shape)), B=B, nd=nd, layout=layout, solver_mode=mode,
        num_eigenvalues=k, max_iter=max_iter, spatial_shape=shape)


def _static_plan(nmat, shape, B, nd, layout, caps=NOCAP, extra=None, **kw):
    rows = _rows(nmat, shape, B, nd, layout, **kw) if extra is None else extra
    return plan_so_bond(nmat=nmat, nvol=int(np.prod(shape)), B=B, nd=nd,
                        consumers=("pairing_static",), caps=caps, layout=layout,
                        layout_reason=None if layout == "dense" else PACKED_REASON,
                        extra_rows=rows)


PACKED_REASON = "so bond layout: sz_packed (S_z conserved)"


class TestSectorSizes(unittest.TestCase):
    """Construction check of ``n_trim`` and ``pp_sector_dims`` against the
    antisymmetric basis and the packed kernel's sectors."""

    def test_n_trim_and_dims(self):
        for shape, trim in (((4, 4, 1), 4), ((5, 3, 1), 1), ((6, 1, 1), 2), ((4, 2, 1), 4)):
            with self.subTest(shape=shape):
                self.assertEqual(so_bond.n_trim(shape), trim)
                self.assertEqual(so_bond.n_trim(shape), _n_trim(shape))
                for nd in (2, 4):
                    got = so_bond.pp_sector_dims(nd, shape)
                    self.assertEqual(tuple(got), _sectors(nd, shape))
                    self.assertEqual(sum(got.dims),
                                     sc._SOAntisymmetricBasis((nd, nd) + shape).dim)
        got = so_bond.pp_sector_dims(2, (4, 4, 1))
        self.assertEqual(got.dims, (6, 16, 6))
        self.assertEqual(got.pair_counts, (1, 2, 1))
        self.assertEqual(got.n_sector_vec, (16, 32, 16))

    def test_pp_sector_dims_match_kernel(self):
        for fx in (FIXTURE_SZ_GRID, FIXTURE_SZ_GRID_4X2, FIXTURE_LZSZ_GRID):
            with self.subTest(fx=FIXTURE_NAMES[fx]):
                kernel = random_kernels(fx, 1)[0]
                got = so_bond.pp_sector_dims(fx.nd, fx.spatial_shape)
                for i, p in enumerate(PP_ORDER):
                    self.assertEqual(got.pair_counts[i] * fx.nvol,
                                     int(np.prod(kernel.sector_shape(p))))
                    self.assertEqual(got.n_sector_vec[i], kernel.sector_op(p).basis.n)
                    self.assertEqual(got.dims[i], kernel.sector_op(p).basis.dim)


class TestArguments(unittest.TestCase):
    """Construction check of the argument refusals of the static rows and
    their estimators (one exact text per rule)."""

    def test_argument_validation(self):
        dense = so_bond.estimate_dense_layout_solver_bytes
        it = so_bond.estimate_sector_iteration_bytes
        sv = so_bond.estimate_sector_solver_bytes
        kw = dict(nmat=8, nvol=16, B=5, nd=2, layout="dense", solver_mode="eigenvalue",
                  num_eigenvalues=10, max_iter=1000, spatial_shape=(4, 4, 1))

        def rows(**over):
            return lambda: so_bond.static_so_bond_rows(**dict(kw, **over))
        cases = (
            (lambda: so_bond.n_trim((4, True, 1)),
             "n_trim: spatial_shape must be three integers >= 1, got (4, True, 1)"),
            (lambda: so_bond.n_trim((4, 4)),
             "n_trim: spatial_shape must be three integers >= 1, got (4, 4)"),
            (lambda: so_bond.n_trim((4, 0, 1)),
             "n_trim: spatial_shape must be three integers >= 1, got (4, 0, 1)"),
            (lambda: so_bond.pp_sector_dims(3, (4, 4, 1)),
             "pp_sector_dims: nd must be an even integer >= 2, got 3"),
            (lambda: so_bond.pp_sector_dims(True, (4, 4, 1)),
             "pp_sector_dims: nd must be an integer >= 2, got True"),
            (lambda: so_bond.pp_sector_dims(2, (4.0, 4, 1)),
             "pp_sector_dims: spatial_shape must be three integers >= 1, got (4.0, 4, 1)"),
            (lambda: dense(60, 4, 100, "eigenvalue", max_iter=1),
             "estimate_dense_layout_solver_bytes: n must be an integer >= 120, got 100"),
            (lambda: dense(28.0, 4, 64, "eigenvalue", max_iter=1),
             "estimate_dense_layout_solver_bytes: dim must be an integer >= 1, got 28.0"),
            (lambda: dense(28, True, 64, "eigenvalue", max_iter=1),
             "estimate_dense_layout_solver_bytes: k must be an integer >= 1, got True"),
            (lambda: dense(28, 4, 64, "power", max_iter=1),
             "estimate_dense_layout_solver_bytes: solver_mode must be 'eigenvalue', "
             "'iteration' or 'both', got 'power'"),
            (lambda: dense(28, 4, 64, "iteration", max_iter=0),
             "estimate_dense_layout_solver_bytes: max_iter must be an integer >= 1, got 0"),
            (lambda: it([6], 16, n_sector_vec=[16], max_iter=1000.0),
             "estimate_sector_iteration_bytes: max_iter must be an integer >= 1, got 1000.0"),
            (lambda: it([6], 16, n_sector_vec=[16], max_iter=True),
             "estimate_sector_iteration_bytes: max_iter must be an integer >= 1, got True"),
            (lambda: it([9], 16, n_sector_vec=[16], max_iter=10),
             "estimate_sector_iteration_bytes: dims[0] = 9 exceeds n_sector_vec[0] / 2 = "
             "16 / 2"),
            (lambda: it([6], 15, n_sector_vec=[16], max_iter=10),
             "estimate_sector_iteration_bytes: n_full_gap must be an integer >= "
             "sum(n_sector_vec) = 16, got 15"),
            (lambda: sv([6], 4, 16, n_sector_vec=[16], dense_fallback=1),
             "estimate_sector_solver_bytes: dense_fallback must be a bool, got 1"),
            (rows(nvol=15),
             "static_so_bond_rows: nvol must equal prod(spatial_shape) = 16, got 15"),
            (rows(layout="packed"),
             "static_so_bond_rows: layout must be 'dense' or 'sz_packed', got 'packed'"),
            (rows(nmat=7), "static_so_bond_rows: nmat must be an even integer >= 2, got 7"),
            (rows(nmat=32.0), "static_so_bond_rows: nmat must be an integer >= 2, got 32.0"),
            (rows(B=3.0), "static_so_bond_rows: B must be an integer >= 1, got 3.0"),
            (rows(nd=3), "static_so_bond_rows: nd must be an even integer >= 2, got 3"),
            (rows(num_eigenvalues=True),
             "static_so_bond_rows: num_eigenvalues must be an integer >= 1, got True"),
            (rows(max_iter=0), "static_so_bond_rows: max_iter must be an integer >= 1, got 0"),
            (rows(solver_mode="power"),
             "static_so_bond_rows: solver_mode must be 'eigenvalue', 'iteration' or 'both', "
             "got 'power'"),
            (rows(spatial_shape=(4, 4, True)),
             "static_so_bond_rows: spatial_shape must be three integers >= 1, got (4, 4, True)"),
            (lambda: _static_plan(8, (4, 4, 1), 5, 2, "dense",
                                  extra=(("resident", "carrier", 16),)),
             "plan_so_bond: extra_rows must hold (phase, name, bytes) tuples with non-empty "
             "string phase and name, a phase other than 'resident' and bytes an integer >= 0, "
             "got ('resident', 'carrier', 16)"),
            (lambda: _static_plan(8, (4, 4, 1), 5, 2, "dense",
                                  extra=(["kernel_static", "g2", 16],)),
             "plan_so_bond: extra_rows must hold (phase, name, bytes) tuples with non-empty "
             "string phase and name, a phase other than 'resident' and bytes an integer >= 0, "
             "got ['kernel_static', 'g2', 16]"),
            (lambda: _static_plan(8, (4, 4, 1), 5, 2, "dense",
                                  extra=(("kernel_static", "", 16),)),
             "plan_so_bond: extra_rows must hold (phase, name, bytes) tuples with non-empty "
             "string phase and name, a phase other than 'resident' and bytes an integer >= 0, "
             "got ('kernel_static', '', 16)"),
            (lambda: _static_plan(8, (4, 4, 1), 5, 2, "dense",
                                  extra=(("kernel_static", "g2", 16.0),)),
             "plan_so_bond: extra_rows must hold (phase, name, bytes) tuples with non-empty "
             "string phase and name, a phase other than 'resident' and bytes an integer >= 0, "
             "got ('kernel_static', 'g2', 16.0)"))
        for call, msg in cases:
            with self.subTest(msg=msg):
                with self.assertRaisesRegex(ValueError, exact(msg)):
                    call()


class TestRows(unittest.TestCase):
    """Construction check of the static rows: every row equals its
    written-out formula in both layouts and every solver mode, the
    eigensolver row follows the branch the solver takes (and reserves the
    failure-only dense fallback of an ARPACK branch), and the plan's peak is
    taken over every phase of its rows."""

    def test_rows_dense_and_packed(self):
        for shape, nd in (((4, 4, 1), 2), ((5, 3, 1), 4)):
            dim = sum(_sectors(nd, shape)[2])
            for k in (2, 10, dim - 1):
                for layout in ("dense", "sz_packed"):
                    for mode in ("eigenvalue", "iteration", "both"):
                        with self.subTest(shape=shape, nd=nd, k=k, layout=layout, mode=mode):
                            got = _rows(8, shape, 5, nd, layout, mode, k, 1000)
                            self.assertEqual(got, _literal_rows(8, shape, 5, nd, layout, mode,
                                                                k, 1000))
                            self.assertTrue(all(type(b) is int for _, _, b in got))

    def test_dense_solver_row_branches(self):
        dense = so_bond.estimate_dense_layout_solver_bytes
        tail = lambda d, k, n: 16 * k * n + 16 * k * d + 24 * n + OV  # noqa: E731
        # k >= dim - 1: the direct dense branch only
        self.assertEqual(dense(28, 27, 64, "eigenvalue", max_iter=1000),
                         _direct(28, 64) + tail(28, 27, 64))
        # ARPACK with the fallback reserved
        self.assertEqual(dense(60, 4, 128, "eigenvalue", max_iter=1000),
                         _arpack(60, 4, 128) + 16 * (4 * 3600 + 120) + 16 * 130 * 60
                         + 8 * 120 + tail(60, 4, 128))
        self.assertEqual(dense(2000, 10, 4000, "eigenvalue", max_iter=1000),
                         _arpack(2000, 10, 4000) + 16 * (4 * 2000 ** 2 + 4000)
                         + 16 * 130 * 2000 + 8 * 4000
                         + tail(2000, 10, 4000))
        # above the fallback bound: ARPACK alone
        self.assertEqual(dense(2001, 10, 4002, "eigenvalue", max_iter=1000),
                         _arpack(2001, 10, 4002) + tail(2001, 10, 4002))
        self.assertEqual(dense(60, 4, 128, "iteration", max_iter=1000),
                         16 * (6 * 60 + 4 * 128) + 24 * 128 + 80 * 1000 + OV)
        self.assertEqual(dense(60, 4, 128, "both", max_iter=7),
                         dense(60, 4, 128, "eigenvalue", max_iter=7)
                         + dense(60, 4, 128, "iteration", max_iter=7))
        self.assertEqual(so_bond.DENSE_FALLBACK_MAX, sc._DENSE_FALLBACK_MAX)
        self.assertEqual(so_bond.DENSE_FALLBACK_MAX, 2000)
        # the 32 x 32 grid at nd = 4, B = 5, Nmat = 256: the kernel phase (dim 8184,
        # ARPACK without a fallback) stays below the pairing phase, which a cap of
        # 16.8 GiB admits
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                plan = _static_plan(256, (32, 32, 1), 5, 4, layout)
                self.assertEqual(plan.peak_phase, "pairing_static")
                self.assertLess(plan.phase_total("kernel_static"),
                                plan.phase_total("pairing_static"))
                _static_plan(256, (32, 32, 1), 5, 4, layout,
                             caps=SOBondCaps(16.8 * 2 ** 30, 16.8, None))

    def test_sector_fallback_row(self):
        sv = so_bond.estimate_sector_solver_bytes
        off = sv([60], 4, 120, n_sector_vec=[120])
        self.assertEqual(off, 130304 + 16 * (6 * 61 + 4))
        self.assertEqual(sv([60], 4, 120, n_sector_vec=[120], dense_fallback=True),
                         off + 16 * (4 * 3600 + 120) + 16 * 130 * 60 + 8 * 120)
        self.assertEqual(sv([60], 4, 120, n_sector_vec=[120], dense_fallback=np.bool_(True)),
                         off + 16 * (4 * 3600 + 120) + 16 * 130 * 60 + 8 * 120)
        for k in (1, 10):
            for dims in ((7,), (20,), (21,), (60,), (7, 20, 21, 60), (0, 21, 3), (2001,)):
                for extra in (0, 5, 4000):
                    nvec = tuple(2 * d + extra for d in dims)
                    for pad in (3, 10 ** 6):
                        nfull = sum(nvec) + pad
                        with self.subTest(dims=dims, k=k, nvec=nvec, pad=pad):
                            for fb in (False, True):
                                self.assertEqual(
                                    sv(list(dims), k, nfull, n_sector_vec=nvec,
                                       dense_fallback=fb),
                                    _sector_solver(dims, k, nfull, nvec, fb))

    def test_fallback_bound_in_rows(self):
        """The rows reserve the dense fallback up to reduced dimension 2000 and
        not above, in both layouts (nd = 2, k = 10)."""
        def solver_row(shape, layout):
            rows = _rows(2, shape, 1, 2, layout)
            self.assertEqual(rows, _literal_rows(2, shape, 1, 2, layout, "eigenvalue", 10, 1000))
            return rows[-1][2]
        # dense: dim = 2 nvol - N_TRIM never equals 2000 (an even nvol has
        # N_TRIM >= 2 and no nvol near 1000 has enough even factors), so the
        # sides are 1998 and 2001
        tail = lambda d, n: 16 * 10 * n + 16 * 10 * d + 24 * n + OV  # noqa: E731
        for shape, d, fallback in (((1000, 1, 1), 1998, True), ((7, 143, 1), 2001, False)):
            with self.subTest(layout="dense", dim=d):
                self.assertEqual(sum(_sectors(2, shape)[2]), d)
                n = 4 * int(np.prod(shape))
                got = solver_row(shape, "dense")
                self.assertEqual(got - _arpack(d, 10, n) - tail(d, n),
                                 _F(d) if fallback else 0)
                self.assertEqual(fallback, got - _arpack(d, 10, n) - tail(d, n) > 0)
        # packed: the S_z = 0 sector has dim nvol, the largest sector
        for shape, dims, fallback in (((40, 50, 1), (998, 2000, 998), True),
                                      ((3, 667, 1), (1000, 2001, 1000), False)):
            with self.subTest(layout="sz_packed", dims=dims):
                _, nvec, got_dims = _sectors(2, shape)
                self.assertEqual(got_dims, dims)
                n = 4 * int(np.prod(shape))
                got = solver_row(shape, "sz_packed")
                rest = (16 * (10 + 30) * n + 16 * 10 * sum(dims) + OV)
                self.assertEqual(got - rest, max(_arpack(dims[1], 10, nvec[1])
                                                 + (_F(dims[1]) if fallback else 0),
                                                 _arpack(dims[0], 10, nvec[0]) + _F(dims[0])))

    def test_retry_working_set_row(self):
        """The ARPACK branch by spectral_shift: the outputs kept into the
        largest-real solve (none for a number, the spectral-radius set for
        "auto", that set and the first eigenvalues unshifted), and a second
        Krylov working set above the dense-fallback bound for every shift;
        both estimators and the static rows. The values accepted are those
        the solver accepts."""
        dense = so_bond.estimate_dense_layout_solver_bytes
        sv = so_bond.estimate_sector_solver_bytes
        tail = lambda d, k, n: 16 * k * n + 16 * k * d + 24 * n + OV  # noqa: E731
        kinds = ((None, "none"), ("none", "none"), ("NONE", "none"), ("auto", "auto"),
                 (1.0, "numeric"), (2, "numeric"), (np.float64(0.5), "numeric"))
        for d, k in ((60, 4), (1998, 1), (2000, 2), (2001, 1), (2044, 1), (2044, 10)):
            n = 2 * d + 8
            for shift, kind in kinds:
                with self.subTest(dim=d, k=k, spectral_shift=shift):
                    self.assertEqual(dense(d, k, n, "eigenvalue", max_iter=10,
                                           spectral_shift=shift),
                                     _arpack(d, k, n, kind) + _F(d) + tail(d, k, n))
                    self.assertEqual(sv([d], k, n, n_sector_vec=[n], dense_fallback=True,
                                        spectral_shift=shift),
                                     _arpack(d, k, n, kind) + _F(d)
                                     + 16 * (k + k) * n + 16 * k * d + OV)
        # the three kinds written out at dim 60, k 4 (no second set) and dim 2044,
        # k 1 (the second set, the same for every kind)
        work60 = 16 * (20 * 60 + 3 * 400 + 14 * 20 + 3 * 60)
        base60 = work60 + 24 * 128 + 6 * 16 * 128 + 16 * 3 * 4 * (60 + 128) + _F(60)
        for shift, retained in ((1.0, 0), ("auto", 16 * 6 * 61), (None, 16 * (6 * 61 + 4))):
            with self.subTest(dim=60, spectral_shift=shift):
                self.assertEqual(dense(60, 4, 128, "eigenvalue", max_iter=10,
                                       spectral_shift=shift),
                                 base60 + retained + tail(60, 4, 128))
        work2044 = 16 * (20 * 2044 + 3 * 400 + 14 * 20 + 3 * 2044)
        base2044 = (2 * work2044 + 24 * 4096 + 6 * 16 * 4096 + 16 * 3 * (2044 + 4096))
        for shift, retained in ((1.0, 0), ("auto", 16 * 6 * 2045),
                                (None, 16 * (6 * 2045 + 1))):
            with self.subTest(dim=2044, spectral_shift=shift):
                self.assertEqual(dense(2044, 1, 4096, "eigenvalue", max_iter=10,
                                       spectral_shift=shift),
                                 base2044 + retained + tail(2044, 1, 4096))
        # the default is the unshifted solve; the direct dense branch and the
        # iteration do not depend on the shift
        self.assertEqual(dense(60, 4, 128, "eigenvalue", max_iter=10),
                         dense(60, 4, 128, "eigenvalue", max_iter=10, spectral_shift=None))
        self.assertEqual(dense(28, 27, 64, "both", max_iter=10, spectral_shift=1.0),
                         dense(28, 27, 64, "both", max_iter=10, spectral_shift="auto"))
        # the static rows pass it through: only the eigensolver row moves
        shape = (8, 8, 1)               # an ARPACK sector (dim 64) in both layouts
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                kw = dict(nmat=8, nvol=64, B=5, nd=2, layout=layout, solver_mode="eigenvalue",
                          num_eigenvalues=2, max_iter=10, spatial_shape=shape)
                rows = {shift: so_bond.static_so_bond_rows(spectral_shift=shift, **kw)
                        for shift in (None, "auto", 1.5)}
                self.assertEqual(so_bond.static_so_bond_rows(**kw), rows[None])
                self.assertEqual(rows[None], _literal_rows(8, shape, 5, 2, layout,
                                                           "eigenvalue", 2, 10))
                self.assertEqual(rows[None][:-1], rows["auto"][:-1])
                self.assertEqual(rows[None][:-1], rows[1.5][:-1])
                self.assertEqual(rows[None][-1][2] - rows["auto"][-1][2], 16 * 2)
                if layout == "dense":
                    # one ARPACK space of dim 2 nvol - N_TRIM = 124
                    self.assertEqual(rows["auto"][-1][2] - rows[1.5][-1][2],
                                     16 * 6 * 125)
                else:
                    # the largest ARPACK sector (dim 64) sets the maximum
                    self.assertEqual(rows["auto"][-1][2] - rows[1.5][-1][2],
                                     16 * 6 * 65)
        too_big = np.array(10 ** 1000, dtype=object)       # float() overflows
        bad_values = (True, np.bool_(True), 0, -1.0, float("nan"), float("inf"), 10 ** 1000,
                      too_big, np.longdouble("inf"), "Auto", "x", "1.0", [1.0])
        for bad in bad_values:
            for name, call in (
                    ("estimate_dense_layout_solver_bytes",
                     lambda: dense(60, 4, 128, "eigenvalue", max_iter=10, spectral_shift=bad)),
                    ("estimate_sector_solver_bytes",
                     lambda: sv([60], 4, 120, n_sector_vec=[120], spectral_shift=bad)),
                    ("static_so_bond_rows",
                     lambda: so_bond.static_so_bond_rows(
                         nmat=8, nvol=64, B=5, nd=2, layout="dense",
                         solver_mode="eigenvalue", num_eigenvalues=2, max_iter=10,
                         spatial_shape=shape, spectral_shift=bad))):
                with self.subTest(name=name, spectral_shift=repr(bad)[:40]):
                    with self.assertRaisesRegex(ValueError, exact(
                            "{}: spectral_shift must be None, \"none\", \"auto\" or a "
                            "positive finite number, got {!r}".format(name, bad))):
                        call()
        # the same values as the solver: accepted and refused alike (np.bool_
        # is refused here although float() converts it; bool is refused by both)
        for value in (None, "none", "NONE", "auto", 1.0, 2, np.float64(0.5), np.int32(3),
                      np.longdouble(2.0)) + bad_values:
            if isinstance(value, np.bool_):
                continue
            with self.subTest(mirror=repr(value)[:40]):
                try:
                    sc._validate_spectral_shift(sc._normalize_spectral_shift(value), "arnoldi")
                    solver_ok = True
                except ValueError:
                    solver_ok = False
                try:
                    dense(60, 4, 128, "eigenvalue", max_iter=10, spectral_shift=value)
                    rows_ok = True
                except ValueError:
                    rows_ok = False
                self.assertEqual(rows_ok, solver_ok)

    def test_fallback_covers_second_working_set(self):
        """At every dimension <= 2000 on an ARPACK branch (k <= dim - 2) the
        dense fallback row is at least one ARPACK working set, so the second
        set (which no solve there starts while a failed one is alive) needs no
        row of its own."""
        for k in (1, 2, 10):
            for d in range(k + 2, so_bond.DENSE_FALLBACK_MAX + 1):
                if so_bond._dense_fallback_bytes(d) < so_bond._arpack_working_set_bytes(d, k):
                    self.fail("dim {} k {}: fallback {} < working set {}".format(
                        d, k, so_bond._dense_fallback_bytes(d),
                        so_bond._arpack_working_set_bytes(d, k)))
            self.assertEqual(so_bond._arpack_retry_bytes(so_bond.DENSE_FALLBACK_MAX, k), 0)
            self.assertEqual(so_bond._arpack_retry_bytes(so_bond.DENSE_FALLBACK_MAX + 1, k),
                             so_bond._arpack_working_set_bytes(so_bond.DENSE_FALLBACK_MAX + 1,
                                                               k))

    def test_peak_over_row_phases(self):
        # norb_phys = 1 on 4 x 4, one channel, every eigenpair: the direct dense
        # branch (dim 28) makes the kernel phase the largest
        shape, dim = (4, 4, 1), 28
        rows = _rows(2, shape, 1, 2, "dense", k=dim - 1)
        plan = _static_plan(2, shape, 1, 2, "dense", extra=rows)
        resident = sum(b for ph, _, b in plan.rows if ph == "resident")
        total = resident + sum(b for ph, _, b in rows if ph == "kernel_static")
        self.assertEqual(plan.phase_total("kernel_static"), total)
        self.assertEqual(plan.peak_phase, "kernel_static")
        self.assertEqual(plan.peak_host, 1.25 * total)
        consumers_only = max(plan.phase_total(c) for c in plan.consumers)
        self.assertGreater(plan.peak_host, 1.25 * consumers_only)
        phases = ("vertex", "chi_static", "pairing_static", "green_static", "g2_static",
                  "kernel_static")
        best = max(phases, key=plan.phase_total)
        self.assertEqual(best, plan.peak_phase)
        self.assertEqual(plan.peak_host, 1.25 * max(plan.phase_total(p) for p in phases))
        self.assertEqual(plan.phase_total("vertex"), resident)
        self.assertIn("in phase kernel_static", plan.reason)
        for bad in ("resident", "chi_dynamic", "kernel", 3):
            with self.subTest(phase=bad):
                with self.assertRaisesRegex(ValueError, exact(
                        "SOBondPlan.phase_total: the plan has no phase {!r}".format(bad))):
                    plan.phase_total(bad)
        cap = 1.1 * total
        with self.assertRaises(ValueError) as cm:
            _static_plan(2, shape, 1, 2, "dense", extra=rows,
                         caps=SOBondCaps(cap, cap / 2 ** 30, None))
        msg = str(cm.exception)
        self.assertTrue(msg.startswith("spin-orbital bond Eliashberg needs an estimated peak "
                                       "of {:.2f} GiB in phase kernel_static ".format(
                                           1.25 * total / 2 ** 30)), msg)
        self.assertIn("  kernel_static  eigensolver  {}".format(rows[-1][2]), msg)

    def test_green_static_rows(self):
        """The green_static rows over Nmat, both branches of the transient
        row, and the digest block of a frequency slice beyond CHUNK_BYTES."""
        won = set()
        for shape in ((4, 4, 1), (16, 16, 1)):
            for nd in (2, 4):
                for nmat in (2, 8, 16, 32, 64, 512):
                    with self.subTest(shape=shape, nd=nd, nmat=nmat):
                        nvol = int(np.prod(shape))
                        got = [r for r in _rows(nmat, shape, 3, nd, "sz_packed")
                               if r[0] == "green_static"]
                        self.assertEqual(got, [r for r in _literal_rows(
                            nmat, shape, 3, nd, "sz_packed", "eigenvalue", 10, 1000)
                            if r[0] == "green_static"])
                        sl, G = 16 * nvol * nd * nd, 16 * nmat * nvol * nd * nd
                        won.add("build" if got[3][2] == 2 * sl + min(G, 256 * 1024)
                                else "digest")
        self.assertEqual(won, {"build", "digest"})
        shape = (2048, 2048, 1)                 # one frequency slice: 256 MiB
        got = [r for r in _rows(2, shape, 1, 2, "sz_packed") if r[0] == "green_static"]
        G = 16 * 2 * 2048 * 2048 * 4
        self.assertEqual(got[3], ("green_static", "transient",
                                  max(2 * G // 2 + 256 * 1024,
                                      3 * 64 * 2 ** 20 + 8 * 2 ** 20 - G)))
        self.assertEqual(got[3], ("green_static", "transient",
                                  _green_transient(2, 2048 * 2048, 2)))

    def test_existing_plans_unchanged(self):
        sets = ((("vertex",), {}), (("chi_static",), {}), (("chi_dynamic",), {}),
                (("chi_dynamic",), dict(keep_chi=True, vertex_batch=True)),
                (("pairing_static",), {}), (("chi_dynamic", "pairing_static"), {}))
        for fx in FIXTURE_NAMES:
            for layout in ("dense", "sz_packed"):
                for consumers, kw in sets:
                    with self.subTest(fx=FIXTURE_NAMES[fx], layout=layout, consumers=consumers,
                                      kw=kw):
                        with quiet():
                            plan = chain_plan(fx, layout, consumers, 4, **kw)
                            same = chain_plan(fx, layout, consumers, 4, extra_rows=(), **kw)
                        # new code against new code: extra_rows=() adds nothing; the
                        # behaviour before extra_rows existed is pinned by the
                        # consumer-only recomputation below
                        self.assertEqual(plan, same)
                        self.assertTrue({ph for ph, _, _ in plan.rows}
                                        <= set(plan.consumers) | {"resident"})
                        # the consumer-only rule of the earlier plans
                        resident = sum(b for ph, _, b in plan.rows if ph == "resident")
                        old = [(1.25 * (resident + sum(b for ph, _, b in plan.rows if ph == c)),
                                c) for c in plan.consumers]
                        peak, phase = max(old, key=lambda p: p[0])
                        self.assertEqual(plan.peak_host, peak)
                        self.assertEqual(plan.peak_phase, phase)
                        self.assertIn("estimated peak {:.3f} GiB in phase {};".format(
                            peak / 2 ** 30, phase), plan.reason)

    def test_require_kernel_boundary(self):
        shape, B, nd = (4, 4, 1), 5, 2
        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                plan = _static_plan(8, shape, B, nd, layout)
                obs = dict(nvol=16, B=B, nd=nd, ND=B * nd * nd, layout=layout)
                plan.require("so_bond_static_kernel", n_gap=64, **obs)
                with self.assertRaisesRegex(ValueError, exact(
                        "so_bond_static_kernel: the admission plan has n_gap = nd^2 nvol = 64 "
                        "but the call has n_gap = 63 ({})".format(plan.reason))):
                    plan.require("so_bond_static_kernel", n_gap=63, **obs)
                with self.assertRaisesRegex(ValueError, exact(
                        "SOBondPlan.require: so_bond_static_kernel observes the fields nvol, "
                        "B, nd, ND, layout, n_gap, got B, ND, layout, nd, nvol")):
                    plan.require("so_bond_static_kernel", **obs)
                bare = _static_plan(8, shape, B, nd, layout, extra=())
                with self.assertRaisesRegex(ValueError, exact(
                        "so_bond_static_kernel: the admission plan has no row 'v_static' "
                        "(closed consumers: vertex, chi_static, pairing_static; {})".format(
                            bare.reason))):
                    bare.require("so_bond_static_kernel", n_gap=64, **obs)


# ---------------------------------------------------------------------------
# traced peaks of the solvers against their rows
# ---------------------------------------------------------------------------

def _trace(fn):
    """``(result, peak)``: ``fn()`` under ``tracemalloc``, the peak above the
    allocations alive before the call."""
    if tracemalloc.is_tracing():
        raise unittest.SkipTest("tracemalloc is already tracing")
    tracemalloc.start()
    try:
        out = fn()
        return out, tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()


def _log_ratio(what, traced, row):
    logger.info("%s: traced %d bytes, row %d bytes, ratio %.3f", what, traced, row,
                traced / row)


def _normal(eigs, seed):
    """A normal matrix with the eigenvalues ``eigs``."""
    d = len(eigs)
    rng = np.random.default_rng(seed)
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return (Q * np.asarray(eigs, dtype=complex)) @ Q.conj().T


def _random(d, seed):
    rng = np.random.default_rng(seed)
    return _normal(rng.normal(size=d) + 0.3j * rng.normal(size=d), seed)


def _slow(d, seed):
    """Leading eigenvalues 1 and 0.97, the rest inside |z| < 0.5: the power
    iteration contracts slowly."""
    rng = np.random.default_rng(seed)
    rest = 0.5 * rng.uniform(size=d - 2) * np.exp(2j * np.pi * rng.uniform(size=d - 2))
    return _normal(np.concatenate([[1.0, 0.97], rest]), seed)


def _full_op(shape, M):
    """``(basis, operator)``: the full antisymmetric basis of ``shape`` and
    an operator whose reduced form is ``M`` (forward and adjoint products,
    no ``matmat``, as the kernel's)."""
    basis = sc._SOAntisymmetricBasis(shape)
    return basis, sector_op_of(0, basis, M, matmat=False).operator


def _init(basis, seed):
    rng = np.random.default_rng(seed)
    return basis.to_full(rng.normal(size=basis.dim) + 1j * rng.normal(size=basis.dim))


#: the dense layout's ARPACK case (dim 60) and direct dense case (dim 12)
ARPACK_SHAPE, DIRECT_SHAPE = (2, 2, 8, 4, 1), (2, 2, 4, 2, 1)
#: the packed cases: a gap of 256 entries holding an ARPACK, a dense and a null sector
PACKED_SHAPE, PACKED_DIMS = (4, 4, 16, 1, 1), (60, 7, 20)
MAX_ITER = 1000


def _packed_ops(make):
    bases = subset_bases(PACKED_SHAPE, PACKED_DIMS)
    ops = [sector_op_of(p, b, np.zeros((b.dim, b.dim)) if i == 2 else make(b.dim, 20 + i),
                        matmat=False)
           for i, (p, b) in enumerate(zip(PP_ORDER, bases))]
    return ops, [b.dim for b in bases], [b.n for b in bases]


def _eig_spy(calls):
    """``np.linalg.eig`` recording the input shape and the traced bytes alive
    when it is called."""
    real = np.linalg.eig

    def spy(a):
        calls.append((np.shape(a), tracemalloc.get_traced_memory()[0]))
        return real(a)
    return spy


def _forced(*a, **kw):
    """An ARPACK solve that fails at once."""
    raise ArpackNoConvergence("forced", np.zeros(0), np.zeros((0, 0)))


def _late(*a, **kw):
    """The real ARPACK solve stopped after one restart: every work array is
    allocated before it fails."""
    return real_eigs(*a, **dict(kw, maxiter=1))


def _negative_leading(d, seed):
    """Six eigenvalues near -5 lead in magnitude, the rest lie in [0, 0.5):
    with k = 1 the unshifted which='LM' solve finds a negative eigenvalue and
    re-solves for the largest real part (spectral_shift="auto")."""
    rng = np.random.default_rng(seed)
    return _normal(np.concatenate([-5.0 - rng.uniform(size=6),
                                   0.5 * rng.uniform(size=d - 6)]), seed)


def _diagonal_op(basis, ev):
    """A kernel on the gap of ``basis`` whose reduced form is ``diag(ev)``
    (forward and adjoint products; cheap at any dimension)."""
    def fwd(x):
        return basis.to_full(ev * basis.to_reduced(np.asarray(x).ravel()))

    def adj(x):
        return basis.to_full(np.conj(ev) * basis.to_reduced(np.asarray(x).ravel()))
    return LinearOperator((basis.n, basis.n), dtype=complex, matvec=fwd, rmatvec=adj)


def _which_recorded(lr_solve, whiches):
    """``eigs`` recording the ``which`` of each call; the which='LR' solve
    runs ``lr_solve``, every other one the real solve."""
    def solve(*a, **kw):
        whiches.append(kw.get("which"))
        return (lr_solve if kw.get("which") == "LR" else real_eigs)(*a, **kw)
    return solve


class TestSolverMemory(unittest.TestCase):
    """Gate P-mem on the eigensolver rows: the traced peak of each solver
    branch the static spin-orbital path can take is at most its row -- the
    dense layout's ARPACK branch (dim 60, k = 4), direct dense branch (dim
    12, k = 11; with the native bytes of its ``np.linalg.eig`` added at the
    moment of the call) and power iteration (1000 steps at most, a slowly
    contracting operator); the packed sector solve with an ARPACK, a dense
    and a null sector, and the packed iteration. After an ARPACK failure,
    immediate or late (the real solve stopped after one restart, its work
    arrays kept alive by the exception), the dense fallback of both layouts
    is measured at its two moments: while ``np.linalg.eig`` runs, the traced
    bytes plus the host LAPACK's native bytes are at most the row (the ratio
    to the row without the fallback is logged); after it returns, the
    traced peak (three ``dim x dim`` arrays, no native share) is at most the
    row minus its native share plus ``16 dim^2``. The host LAPACK's native
    bytes are at most that share. The cases include a late failure at dim
    161, k = 1, where the ARPACK term leaves little slack. The operators are
    built before the trace starts."""

    def test_solver_rows_bound_traced(self):
        dense = so_bond.estimate_dense_layout_solver_bytes
        for shape, k, label in ((ARPACK_SHAPE, 4, "arpack"), (DIRECT_SHAPE, 11, "direct")):
            with self.subTest(branch=label):
                basis, A = _full_op(shape, _random(sc._SOAntisymmetricBasis(shape).dim, 3))
                info, calls = {}, []
                with mock.patch.object(np.linalg, "eig", side_effect=_eig_spy(calls)), quiet():
                    _, peak = _trace(lambda: sc._so_solve_eigenvalue(
                        A, basis, k, "arnoldi", info_out=info))
                d = basis.dim
                self.assertEqual(info["selection"] == "dense-LR", label == "direct")
                self.assertEqual(d, 60 if label == "arpack" else 12)
                row = dense(d, k, basis.n, "eigenvalue", max_iter=MAX_ITER)
                _log_ratio("dense layout, {} branch".format(label), peak, row)
                self.assertLessEqual(peak, row)
                if label == "direct":
                    self.assertEqual([c[0] for c in calls], [(d, d)])
                    at_eig = calls[0][1] + 16 * (d * d + d) + _host_native(d)
                    _log_ratio("dense layout, direct branch at eig (traced + native)",
                               at_eig, row)
                    self.assertLessEqual(at_eig, row)
                else:
                    self.assertEqual(calls, [])
        with self.subTest(branch="iteration"):
            basis, A = _full_op(ARPACK_SHAPE, _slow(60, 4))
            init, info = _init(basis, 5), {}
            with quiet():
                _, peak = _trace(lambda: sc._so_solve_iteration(
                    A, basis, init, MAX_ITER, 0.0, 1e-10, info_out=info))
            self.assertGreaterEqual(info["n_iter"], 200)
            row = dense(basis.dim, 4, basis.n, "iteration", max_iter=MAX_ITER)
            _log_ratio("dense layout, iteration ({} steps)".format(info["n_iter"]), peak, row)
            self.assertLessEqual(peak, row)
        with self.subTest(branch="packed eigenvalue"):
            ops, dims, nvec = _packed_ops(_random)
            n_full = int(np.prod(PACKED_SHAPE))
            with quiet():
                s, peak = _trace(lambda: sc._so_solve_sectors(ops, 4, "arnoldi",
                                                              gap_shape=PACKED_SHAPE))
            self.assertTrue(s.converged)
            self.assertEqual([lead.selection for lead in s.leaders][1:],
                             ["dense-LR", "null-operator"])
            self.assertIn(s.leaders[0].selection, ("LM", "LR_retry"))
            row = so_bond.estimate_sector_solver_bytes(dims, 4, n_full, n_sector_vec=nvec,
                                                       dense_fallback=True)
            _log_ratio("packed sector solve", peak, row)
            self.assertLessEqual(peak, row)
        with self.subTest(branch="packed iteration"):
            ops, dims, nvec = _packed_ops(_slow)
            n_full = int(np.prod(PACKED_SHAPE))
            init = _init(sc._SOAntisymmetricBasis(PACKED_SHAPE), 6).reshape(PACKED_SHAPE)
            with quiet():
                it, peak = _trace(lambda: sc._so_iterate_sectors(
                    ops, init, gap_shape=PACKED_SHAPE, max_iter=MAX_ITER, alpha=0.0,
                    tol=1e-10))
            self.assertGreaterEqual(it.n_iter, 200)
            row = so_bond.estimate_sector_iteration_bytes(dims, n_full, n_sector_vec=nvec,
                                                          max_iter=MAX_ITER)
            _log_ratio("packed iteration ({} steps)".format(it.n_iter), peak, row)
            self.assertLessEqual(peak, row)

    def _check_fallback(self, what, run, row, d, selection):
        """Run ``run`` with ``np.linalg.eig`` spied, then the bounds of the
        fallback while ``np.linalg.eig`` runs and after it returns."""
        calls = []
        with mock.patch.object(np.linalg, "eig", side_effect=_eig_spy(calls)), quiet():
            out, peak = _trace(run)
        # the fallback ran: one dense diagonalization of the reduced operator
        self.assertEqual(selection(out), "dense-LR")
        self.assertEqual([c[0] for c in calls], [(d, d)])
        # after eig: the matrix, the eigenvectors and their reordered copy
        after_eig = row - _native_share(d) + 16 * d * d
        self.assertLessEqual(_host_native(d), _native_share(d))
        _log_ratio(what + ", traced against the row's share after eig", peak, after_eig)
        self.assertLessEqual(peak, after_eig)
        at_eig = calls[0][1] + 16 * (d * d + d) + _host_native(d)
        _log_ratio(what + ", traced + native at eig", at_eig, row)
        _log_ratio(what + ", traced + native at eig, without the fallback row", at_eig,
                   row - _F(d))
        self.assertLessEqual(at_eig, row)

    def test_dense_fallback_traced(self):
        for fail, failure in ((_forced, "immediate"), (_late, "late")):
            with self.subTest(layout="dense", failure=failure):
                basis, A = _full_op(ARPACK_SHAPE, _random(60, 7))
                info = {}
                row = so_bond.estimate_dense_layout_solver_bytes(60, 4, basis.n, "eigenvalue",
                                                                 max_iter=MAX_ITER)
                with mock.patch.object(sc, "eigs", side_effect=fail) as patched:
                    self._check_fallback(
                        "dense layout, {} failure".format(failure),
                        lambda: sc._so_solve_eigenvalue(A, basis, 4, "arnoldi",
                                                        info_out=info),
                        row, 60, lambda _: info["selection"])
                self.assertTrue(patched.called)
            with self.subTest(layout="sz_packed", failure=failure):
                bases = subset_bases(PACKED_SHAPE, (60,))
                ops = [sector_op_of(1, bases[0], _random(60, 8), matmat=False)]
                n_full = int(np.prod(PACKED_SHAPE))
                row = so_bond.estimate_sector_solver_bytes([60], 4, n_full,
                                                           n_sector_vec=[bases[0].n],
                                                           dense_fallback=True)
                with mock.patch.object(sc, "eigs", side_effect=fail) as patched:
                    self._check_fallback(
                        "packed sector, {} failure".format(failure),
                        lambda: sc._so_solve_sectors(ops, 4, "arnoldi",
                                                     gap_shape=PACKED_SHAPE),
                        row, 60, lambda s: s.leaders[0].selection)
                self.assertTrue(patched.called)
        # k = 1 at dim 161: the k-dependent ARPACK term leaves little slack, so
        # the bound after eig is the one that holds
        with self.subTest(layout="dense", failure="late", dim=161):
            shape = (2, 2, 9, 9, 1)
            basis, A = _full_op(shape, _random(161, 11))
            self.assertEqual(basis.dim, 161)
            info = {}
            row = so_bond.estimate_dense_layout_solver_bytes(161, 1, basis.n, "eigenvalue",
                                                             max_iter=MAX_ITER)
            with mock.patch.object(sc, "eigs", side_effect=_late) as patched:
                self._check_fallback(
                    "dense layout dim 161, late failure",
                    lambda: sc._so_solve_eigenvalue(A, basis, 1, "arnoldi", info_out=info),
                    row, 161, lambda _: info["selection"])
            self.assertTrue(patched.called)

    def test_shift_estimate_traced(self):
        """The sets kept alive into the largest-real solve: the which='LM'
        spectral-radius solve of spectral_shift="auto" (k = 1, min(6, dim - 2) = 6 pairs) and,
        for the unshifted solve's re-solve, the first which='LM' set as well;
        solved, and after an immediate or a late failure of the which='LR'
        solve (the dense fallback), in both layouts (dim 60). A numeric
        shift (no spectral-radius solve, no retry) is bounded by its row
        without the second ARPACK working set."""
        n_full = int(np.prod(PACKED_SHAPE))
        for shift, make, order, label in (("auto", _random, ["LM", "LR"], "LR"),
                                          (None, _negative_leading, ["LM", "LM", "LR"],
                                           "LR_retry"),
                                          (1.0, _random, ["LR"], "LR")):
            for lr_solve, failure in ((real_eigs, "none"), (_forced, "immediate"),
                                      (_late, "late")):
                for layout in ("dense", "sz_packed"):
                    with self.subTest(spectral_shift=shift, failure=failure, layout=layout):
                        info, whiches = {}, []
                        if layout == "dense":
                            basis, A = _full_op(ARPACK_SHAPE, make(60, 7))
                            row = so_bond.estimate_dense_layout_solver_bytes(
                                60, 1, basis.n, "eigenvalue", max_iter=MAX_ITER,
                                spectral_shift=shift)

                            def run():
                                return sc._so_solve_eigenvalue(
                                    A, basis, 1, "arnoldi", spectral_shift=shift,
                                    info_out=info)

                            def selection(_):
                                return info["selection"]
                        else:
                            bases = subset_bases(PACKED_SHAPE, (60,))
                            ops = [sector_op_of(1, bases[0], make(60, 8), matmat=False)]
                            row = so_bond.estimate_sector_solver_bytes(
                                [60], 1, n_full, n_sector_vec=[bases[0].n],
                                dense_fallback=True, spectral_shift=shift)

                            def run():
                                return sc._so_solve_sectors(ops, 1, "arnoldi",
                                                            gap_shape=PACKED_SHAPE,
                                                            spectral_shift=shift)

                            def selection(s):
                                return s.leaders[0].selection
                        what = "{} layout, spectral_shift={}, {} failure".format(
                            layout, shift, failure)
                        with mock.patch.object(sc, "eigs",
                                               side_effect=_which_recorded(lr_solve, whiches)):
                            if failure == "none":
                                with quiet():
                                    out, peak = _trace(run)
                                self.assertEqual(selection(out), label)
                                _log_ratio(what, peak, row)
                                self.assertLessEqual(peak, row)
                            else:
                                self._check_fallback(what, run, row, 60, selection)
                        self.assertEqual(whiches, order)

    def test_retry_above_fallback_bound_traced(self):
        """Above the dense-fallback bound (reduced dim 2044, k = 1) a late
        failure runs a second which='LM' solve for placeholder vectors while
        the failed solve's arrays are alive, for every shift: "auto" (the
        spectral-radius solve, the failed LR solve, the placeholder solve), a
        numeric shift (the failed LR solve, the placeholder solve) and the
        unshifted solve (its first which='LM' solve failing late, then the
        placeholder solve). The traced peak is within the row of that shift,
        which reserves the second working set, in both layouts."""
        shape = (2, 2, 32, 32, 1)
        basis = sc._SOAntisymmetricBasis(shape)
        d = basis.dim
        self.assertEqual(d, 2044)
        rng = np.random.default_rng(3)
        ev = rng.normal(size=d) + 0.3j * rng.normal(size=d)
        A = _diagonal_op(basis, ev)
        for shift, order, failing in (("auto", ["LM", "LR", "LM"], "LR"),
                                      (1.0, ["LR", "LM"], "LR"),
                                      (None, ["LM", "LM"], "first")):
            for layout in ("dense", "sz_packed"):
                with self.subTest(spectral_shift=shift, layout=layout):
                    info, whiches = {}, []
                    if layout == "dense":
                        row = so_bond.estimate_dense_layout_solver_bytes(
                            d, 1, basis.n, "eigenvalue", max_iter=MAX_ITER,
                            spectral_shift=shift)

                        def run():
                            return sc._so_solve_eigenvalue(A, basis, 1, "arnoldi",
                                                           spectral_shift=shift,
                                                           info_out=info)

                        def selection(_):
                            return info["selection"]
                    else:
                        op = sc.SectorOp(0, A, basis, True)
                        row = so_bond.estimate_sector_solver_bytes(
                            [d], 1, basis.n, n_sector_vec=[basis.n], dense_fallback=True,
                            spectral_shift=shift)

                        def run():
                            return sc._so_solve_sectors([op], 1, "arnoldi", gap_shape=shape,
                                                        spectral_shift=shift)

                        def selection(s):
                            return s.leaders[0].selection
                    if failing == "LR":
                        solve = _which_recorded(_late, whiches)
                    else:
                        # the first solve fails late, the placeholder solve runs
                        def solve(*a, **kw):
                            whiches.append(kw.get("which"))
                            return (_late if len(whiches) == 1 else real_eigs)(*a, **kw)
                    with mock.patch.object(sc, "eigs", side_effect=solve), quiet():
                        out, peak = _trace(run)
                    self.assertEqual(selection(out), "not-converged")
                    self.assertEqual(whiches, order)
                    _log_ratio("{} layout dim {}, spectral_shift={}, late failure".format(
                        layout, d, shift), peak, row)
                    self.assertLessEqual(peak, row)

    def test_fallback_native_share(self):
        # the row's native share of the fallback covers this host's LAPACK at the
        # dimensions the traced tests use, at small dimensions (where the optimal
        # zgeev workspace of some builds has a fixed floor) and up to the bound
        self.assertEqual(so_bond._dense_fallback_bytes(28),
                         16 * (4 * 784 + 56) + 16 * 4352 + 8 * 56)
        for d in (60, 78, 124, 598, 1998, 2000):
            with self.subTest(dim=d):
                self.assertEqual(_F(d) - _native_share(d), 16 * (2 * d * d + d))
                self.assertLessEqual(_host_native(d), _native_share(d))

    @unittest.skipUnless(hasattr(lapack, "zgeev_lwork"),
                         "scipy.linalg.lapack has no zgeev_lwork workspace query")
    def test_small_dim_workspace_floor(self):
        for d in range(16, 41):
            with self.subTest(dim=d):
                lwork = int(lapack.zgeev_lwork(d, compute_vl=0, compute_vr=1)[0].real)
                self.assertLessEqual(lwork, max(130 * d, 4352))
                self.assertLessEqual(_host_native(d), _native_share(d))

def _zeeman_green_inputs(shape, nmat):
    """``(epsilon_k, so_run)`` of the Zeeman fixture (field as Extern) on
    ``shape``, from its literal lists (no files)."""
    fx = F.ZEEMAN_EXTERN
    hr = {}
    for R, t in fx.hopping:
        for s in (0, 1):
            hr[(R, (s, s))] = complex(t)
    extern = {(R, ab): complex(v) for R, ab, v in fx.extern}
    tables = {t: {(R, ab): v for R, ab, v in rows} for t, rows in fx.interactions}
    inp = {"mode": {"param": {"CellShape": list(shape), "coeff_extern": fx.coeff_extern}}}
    ctx = spin_orbital.build_so_entry_context(None, tables, 1, logger)
    run = sc._so_run_inputs(inp, hr, extern, 1, 2.0, nmat, fx.filling, context=ctx)
    k = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in shape]
    with quiet():
        eps = sc._so_epsilon_k(inp, *k, hr, extern, 1)
    return eps, run


def _two_orbital_green_inputs(shape, nmat):
    """``(epsilon_k, so_run)`` of a two-orbital (nd = 4) S_z-conserving band
    on ``shape``: the Zeeman fixture's hopping on orbital 0, 0.8 of it on
    orbital 1, an on-site inter-orbital hopping 0.2 and the field as Extern
    on both orbitals."""
    fx = F.ZEEMAN_EXTERN
    hr = {}
    for a, scale in ((0, 1.0), (1, 0.8)):
        for R, t in fx.hopping:
            for s in (0, 1):
                hr[(R, (2 * a + s, 2 * a + s))] = complex(t) * scale
    for s in (0, 1):
        hr[((0, 0, 0), (s, 2 + s))] = hr[((0, 0, 0), (2 + s, s))] = 0.2 + 0j
    extern = {((0, 0, 0), (a, a)): 1.0 + 0j for a in (0, 1)}
    tables = {"CoulombIntra": {((0, 0, 0), (a, a)): 1.0 for a in (0, 1)}}
    inp = {"mode": {"param": {"CellShape": list(shape), "coeff_extern": fx.coeff_extern}}}
    ctx = spin_orbital.build_so_entry_context(None, tables, 2, logger)
    run = sc._so_run_inputs(inp, hr, extern, 2, 2.0, nmat, fx.filling, context=ctx)
    k = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in shape]
    with quiet():
        eps = sc._so_epsilon_k(inp, *k, hr, extern, 2)
    return eps, run


class TestGreenMemory(unittest.TestCase):
    """Gate P-mem on the green_static rows: the traced peak of the Green
    function build of the static spin-orbital bond path (eigendecomposition,
    mu, carrier and its fp_green digest) is at most the phase's rows, with no
    margin, on small and large grids and on the digest-dominated ones
    (Nmat <= 16, and a 32x32 grid whose digest blocks reach G/2), and on a
    two-orbital band (nd = 4)."""

    def test_green_static_row_traced(self):
        fx = F.ZEEMAN_EXTERN
        for shape, nmat in (((4, 4, 1), 32), ((16, 16, 1), 64), ((8, 8, 1), 16),
                            ((32, 32, 1), 32)):
            with self.subTest(shape=shape, nmat=nmat):
                eps, run = _zeeman_green_inputs(shape, nmat)

                def build():
                    return sc._so_bond_green(eps, run, norb_phys=1, filling=fx.filling,
                                             beta=2.0, nmat=nmat, coeff_tail=1.0)
                build()                         # first-call allocations
                out, peak = _trace(build)
                self.assertEqual(out.carrier.full_sc.shape, (2, 2) + shape + (nmat,))
                row = sum(b for ph, _, b in _rows(nmat, shape, 3, 2, "sz_packed")
                          if ph == "green_static")
                _log_ratio("green_static {} nmat {}".format(shape, nmat), peak, row)
                self.assertLessEqual(peak, row)
                G = 16 * nmat * int(np.prod(shape)) * 4
                if shape == (4, 4, 1):
                    # the carrier build alone does not bound it
                    self.assertGreater(peak, 4 * G + OV)
        with self.subTest(shape=(6, 6, 1), nmat=32, nd=4):
            shape, nmat = (6, 6, 1), 32
            eps, run = _two_orbital_green_inputs(shape, nmat)

            def build4():
                return sc._so_bond_green(eps, run, norb_phys=2, filling=fx.filling, beta=2.0,
                                         nmat=nmat, coeff_tail=1.0)
            build4()
            out, peak = _trace(build4)
            self.assertEqual(out.carrier.full_sc.shape, (4, 4) + shape + (nmat,))
            self.assertTrue(out.conserving)
            row = sum(b for ph, _, b in _rows(nmat, shape, 3, 4, "sz_packed")
                      if ph == "green_static")
            _log_ratio("green_static {} nmat {} nd 4".format(shape, nmat), peak, row)
            self.assertLessEqual(peak, row)


if __name__ == "__main__":
    unittest.main()
