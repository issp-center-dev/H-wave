"""Peak host memory of the S_z-packed spin-orbital bond layout.

Every measurement is an in-process ``tracemalloc`` trace started after the
inputs are built; it is compared with the admission row (or estimate) that
is meant to bound it:

* the static chain (``so_bond_vertex``, the static ``so_bond_bubble``,
  ``so_bond_dress``, ``so_bond_pairing_vertex_parts``) in both layouts against
  the plan of the consumer ``pairing_static`` (``plan.peak_host / 1.25`` is
  the largest phase sum of the rows);
* the packed store of a ``chi_dynamic`` plan against its ``chibar`` row;
* the sector eigen solve (``hwave.sc._so_solve_sectors``) against
  ``so_bond.estimate_sector_solver_bytes`` (plus the kernel's
  ``matvec_workspace`` on the real kernel, whose products run inside the
  solve);
* the construction and the sector products of ``so_bond.PackedBondPairKernel``
  against ``so_bond.estimate_packed_kernel_workspace_bytes``;
* the packed ``so_bond.p_b`` against its one-slice workspace.

The ratio of a trace to its row is logged (a row far above its trace hides
nothing, but says the row is loose); only the upper bounds are asserted.

``tracemalloc`` sees the allocations made through Python's allocators
(numpy arrays included), not the buffers LAPACK, ARPACK or the FFT library
allocate in C; the rows that stand for those (the eigensolver work arrays,
the FFT intermediates) are derived from the library documentation, not
measured. The allowances are the same for every numpy version.
"""
import logging
import tracemalloc
import types
import unittest

import numpy as np

from hwave import sc
from hwave.solver import so_bond
from hwave.solver.so_bond import (PackedArrays, PackedBondPairKernel, SOBondCaps,
                                  estimate_packed_kernel_workspace_bytes,
                                  estimate_sector_solver_bytes, open_store, p_b, plan_so_bond,
                                  select_layout, slot_geometry, so_bond_bubble, so_bond_dress,
                                  so_bond_pairing_vertex_parts, so_bond_topology, so_bond_vertex,
                                  sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_fixtures import (canonical_green, exact, quiet, sector_op_of,
                                    static_pair_kernels, subset_bases)
from tests.so_bond_oracles import CELL, FIXTURE_LZSZ_GRID, PP_ORDER, sz_of_mode

logger = logging.getLogger(__name__)
NOCAP = SOBondCaps(None, None, None)
#: the Matsubara count of the static chain: small, so that the rows of whole
#: (nvol, ND, ND) objects dominate the bubble's per-shift pair block
NMAT = 4


def _trace(fn):
    """``(result, peak)``: ``fn()`` under ``tracemalloc``; ``peak`` is the
    traced peak above the allocations that were alive before the call."""
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


# ---------------------------------------------------------------------------
# the static chain on the 5 x 3 grid with norb_phys = 2 (B = 5, ND = 80)
# ---------------------------------------------------------------------------

class _Chain:
    """The inputs of the static chain of :data:`FIXTURE_LZSZ_GRID`."""

    def __init__(self):
        fx = FIXTURE_LZSZ_GRID
        self.fx = fx
        self.compiled = compile_interaction(fx.tables_dict, fx.norb_phys)
        self.topo = so_bond_topology(self.compiled, CELL)
        self.geom = slot_geometry(self.topo, fx.spatial_shape)
        layout, self.reason = select_layout(self.compiled, fx.transfer_dict, override="auto")
        assert layout == "sz_packed", self.reason
        self.sm = sz_sector_map(self.topo, fx.nd)
        self.green = canonical_green(fx, NMAT)[None].copy()

    def plan(self, layout, consumers):
        return plan_so_bond(nmat=NMAT, nvol=self.fx.nvol, B=self.geom.B, nd=self.fx.nd,
                            consumers=consumers, caps=NOCAP, layout=layout,
                            layout_reason=self.reason if layout == "sz_packed" else None)

    def run(self, layout, plan):
        """Vertex, static bubble, dressing (a batch of one) and the two
        parts, every intermediate kept alive to the end."""
        fx, shape = self.fx, self.fx.spatial_shape
        gamma = so_bond_vertex(self.compiled, self.topo, shape, plan=plan)
        if layout == "dense":
            chibar = np.empty((fx.nvol, plan.ND, plan.ND), dtype=np.complex128)
        else:
            chibar = PackedArrays.zeros(self.sm, (fx.nvol,))
        so_bond_bubble(self.green, None, fx.beta, self.topo, shape, plan=plan, freq="static",
                       sink=chibar)
        if layout == "dense":
            chi = so_bond_dress(chibar[None], gamma, plan=plan, l0=NMAT // 2,
                                spatial_shape=shape)[0]
        else:
            batch = PackedArrays(self.sm, tuple(b[None] for b in chibar.blocks))
            chi = so_bond_dress(batch, gamma, plan=plan, l0=NMAT // 2, spatial_shape=shape)
            chi = PackedArrays(self.sm, tuple(b[0] for b in chi.blocks))
        parts = so_bond_pairing_vertex_parts(gamma, chi, self.geom, plan=plan)
        return gamma, chibar, chi, parts


class TestStaticChainMemory(unittest.TestCase):
    """Gate P-mem on the static chain of the 5x3 grid with norb_phys = 2 (B
    = 5, ND = 80) at Nmat = 4: the traced peak of the packed chain is at most
    the packed plan's largest phase sum (``peak_host / 1.25``) and below the
    traced peak of the dense chain (the dense chain stays within its own
    plan as well). The ratio traced / row is logged (expected >= 0.3)."""

    @classmethod
    def setUpClass(cls):
        cls.chain = _Chain()
        cls.peaks, cls.plans = {}, {}
        for layout in ("dense", "sz_packed"):
            plan = cls.chain.plan(layout, ("pairing_static",))
            with quiet():
                _, peak = _trace(lambda: cls.chain.run(layout, plan))
            cls.peaks[layout], cls.plans[layout] = peak, plan

    def test_packed_within_plan(self):
        row = self.plans["sz_packed"].peak_host / 1.25
        _log_ratio("static chain, sz_packed", self.peaks["sz_packed"], row)
        self.assertLessEqual(self.peaks["sz_packed"], row)

    def test_dense_within_plan(self):
        row = self.plans["dense"].peak_host / 1.25
        _log_ratio("static chain, dense", self.peaks["dense"], row)
        self.assertLessEqual(self.peaks["dense"], row)

    def test_packed_below_dense(self):
        self.assertLess(self.peaks["sz_packed"], self.peaks["dense"])


class TestStoreRow(unittest.TestCase):
    """Construction check: the packed store of a ``chi_dynamic`` plan holds
    exactly its ``chibar`` row, 6/16 of the dense row."""

    def test_store_bytes(self):
        chain = _Chain()
        rows = {}
        for layout in ("dense", "sz_packed"):
            plan = chain.plan(layout, ("chi_dynamic",))
            rows[layout] = dict((name, b) for ph, name, b in plan.rows if ph == "chi_dynamic")
            if layout == "sz_packed":
                self.assertEqual(plan.slots, ("chibar",))
                with open_store(plan, chain.sm) as store:
                    nbytes = store.nbytes
        self.assertEqual(nbytes, rows["sz_packed"]["chibar"])
        self.assertEqual(16 * nbytes, 6 * rows["dense"]["chibar"])


# ---------------------------------------------------------------------------
# the sector eigen solve
# ---------------------------------------------------------------------------

def _gap(B):
    """A static gap ``(nd = 4, 4, 16 B, 1, 1)``: its length grows with the
    channel count ``B`` of the grids the solve is sized for."""
    return (4, 4, 16 * B, 1, 1)


def _matrix(d, seed):
    """A normal ``d x d`` matrix with complex random eigenvalues."""
    rng = np.random.default_rng(seed)
    ev = rng.normal(size=d) + 0.3j * rng.normal(size=d)
    Q, _ = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))
    return (Q * ev) @ Q.conj().T


#: (dims, k, null sectors): the dense branch (dim <= max(2k + 1, 20)), the
#: ARPACK branch (k = 10: dim > 21) and null operators; dims 200 need B >= 3
SOLVER_CASES = (((7,), 1, ()), ((7,), 10, ()), ((20,), 1, ()), ((20,), 10, ()),
                ((7, 20, 7), 10, ()),
                ((60,), 10, ()), ((29, 60, 29), 10, ()), ((200,), 10, ()), ((200,), 1, ()),
                ((20,), 10, (0,)), ((200,), 10, (0,)), ((29, 60), 10, (1,)))


class TestSectorSolverMemory(unittest.TestCase):
    """Gate P-mem on the sector eigen solve: the traced peak of
    ``_so_solve_sectors`` is at most ``estimate_sector_solver_bytes`` on
    synthetic sectors (the dense branch at dims 7 and 20 with k = 1 and 10;
    ARPACK at dims 60, 29/60/29 and 200; null operators, alone and next to an
    ARPACK sector) on gaps of three lengths, and on the real static kernel of
    the 5x3 grid with norb_phys = 2 (pp dimensions 29 / 60 / 29, ARPACK at k =
    10 including the largest-real re-solve) at most the estimate plus the
    kernel's ``matvec_workspace``. The sector operators are built before
    the trace starts. Ratios are logged. (On the real kernel the adjoint
    operator-norm estimate makes 40 of the 504 kernel applications of the
    solve, at most 60 per sector.)"""

    def test_synthetic_sectors(self):
        for B in (1, 3, 5):
            shape = _gap(B)
            full_dim = sc._SOAntisymmetricBasis(shape).dim
            for dims, k, null in SOLVER_CASES:
                if sum(dims) > full_dim:
                    continue
                with self.subTest(B=B, dims=dims, k=k, null=null):
                    bases = subset_bases(shape, dims)
                    # forward and adjoint products only, as the kernel's sector_op
                    ops = [sector_op_of(p, b, np.zeros((b.dim, b.dim)) if i in null
                                        else _matrix(b.dim, 10 * B + i), matmat=False)
                           for i, (p, b) in enumerate(zip(PP_ORDER, bases))]
                    row = estimate_sector_solver_bytes(
                        [b.dim for b in bases], k, int(np.prod(shape)),
                        n_sector_vec=[b.n for b in bases])
                    with quiet():
                        s, peak = _trace(lambda: sc._so_solve_sectors(
                            ops, k, "arnoldi", gap_shape=shape))
                    self.assertTrue(s.converged)
                    want = ["null-operator" if i in null else
                            "dense-LR" if d <= max(2 * k + 1, 20) else "arpack"
                            for i, d in enumerate(dims)]
                    got = [lead.selection if lead.selection in ("null-operator", "dense-LR")
                           else "arpack" for lead in s.leaders]
                    self.assertEqual(got, want)
                    _log_ratio("sector solve B={} dims={} k={} null={}".format(
                        B, dims, k, null), peak, row)
                    self.assertLessEqual(peak, row)

    def test_real_kernel(self):
        packed = static_pair_kernels(FIXTURE_LZSZ_GRID)[0]
        ops = [packed.sector_op(p) for p in PP_ORDER]
        dims = [op.basis.dim for op in ops]
        self.assertEqual(dims, [29, 60, 29])
        mw = packed.memory_table()["matvec_workspace"]
        n_full = int(np.prod(packed.gap_shape))
        k = 10
        row = estimate_sector_solver_bytes(dims, k, n_full,
                                           n_sector_vec=[op.basis.n for op in ops])
        with quiet():
            s, peak = _trace(lambda: sc._so_solve_sectors(
                ops, k, "arnoldi", gap_shape=packed.gap_shape))
        self.assertTrue(s.converged)
        self.assertNotIn("dense-LR", s.selection)
        self.assertIn("LR_retry", s.selection)
        _log_ratio("sector solve, real kernel, k={}".format(k), peak, row + mw)
        self.assertLessEqual(peak, row + mw)


# ---------------------------------------------------------------------------
# the packed pairing kernel
# ---------------------------------------------------------------------------

#: channel displacements and reverse maps of B = 1, 3 and 5 bond channels
CHANNELS = {1: ([(0, 0, 0)], [0]),
            3: ([(0, 0, 0), (1, 0, 0), (-1, 0, 0)], [0, 2, 1]),
            5: ([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)], [0, 2, 1, 4, 3])}


def _kernel_inputs(B, nfreq, norb=2, shape=(5, 3, 1), seed=0):
    """``(V, G2, view, sector_map, plan)``: a random packed vertex and a
    random pair propagator that conserves the pair S_z exactly."""
    nd, nvol = 2 * norb, int(np.prod(shape))
    delta_r, reverse = CHANNELS[B]
    view = types.SimpleNamespace(n_channels=B, delta_r=delta_r, reverse=np.array(reverse))
    sm = sz_sector_map(view, nd)
    plan = plan_so_bond(nmat=2, nvol=nvol, B=B, nd=nd, consumers=("vertex",), caps=NOCAP,
                        layout="sz_packed", layout_reason="so bond layout: sz_packed (test)")
    rng = np.random.default_rng(seed)
    V = PackedArrays(sm, tuple(rng.normal(size=(nfreq, nvol, n, n))
                               + 1j * rng.normal(size=(nfreq, nvol, n, n)) for n in sm.sizes))
    g2 = (rng.normal(size=(nd,) * 4 + shape + (nfreq,))
          + 1j * rng.normal(size=(nd,) * 4 + shape + (nfreq,)))
    pp = np.array([[int(round(sz_of_mode(b, norb) + sz_of_mode(c, norb))) for c in range(nd)]
                   for b in range(nd)])
    g2[pp[:, None, :, None] != pp[None, :, None, :]] = 0
    return V, g2, view, sm, plan


#: rows of the kernel's two phases: its construction and one sector product
BUILD_ROWS = ("hoisted_blocks", "g2_sectors", "hoist_workspace", "g2_gather", "call_overhead")
APPLY_ROWS = ("matvec_workspace", "call_overhead")


class TestKernelWorkspace(unittest.TestCase):
    """Gate P-mem on the packed pairing kernel (norb_phys = 2 on the 5x3
    grid, B = 1, 3 and 5 channels; nfreq = 1, and 4 for B = 1; and B = 3 on
    a 12x6 grid, where ``matvec_workspace`` is several times the call
    overhead, so the bound tests its scaling term): the traced
    peak of the constructor (hoisting) is at most the build rows of
    ``estimate_packed_kernel_workspace_bytes`` (the hoisted blocks, the
    per-sector pair propagator, the hoisting workspace of one channel pair,
    the gathered pair propagator before its copy and the call overhead); the
    traced peak of one forward and one adjoint product of each sector
    operator at most ``matvec_workspace`` plus the call overhead. Ratios are
    logged."""

    def test_build_and_products(self):
        for B, nfreq, shape in ((1, 1, (5, 3, 1)), (3, 1, (5, 3, 1)), (5, 1, (5, 3, 1)),
                                (1, 4, (5, 3, 1)), (3, 1, (12, 6, 1))):
            with self.subTest(B=B, nfreq=nfreq, shape=shape):
                V, g2, view, sm, plan = _kernel_inputs(B, nfreq, shape=shape)
                K, peak = _trace(lambda: PackedBondPairKernel(
                    V, g2, view, sm, plan=plan, spatial_shape=shape, beta=2.0, nfreq=nfreq))
                table = K.memory_table()
                self.assertEqual(table, estimate_packed_kernel_workspace_bytes(
                    B, [K.sector_shape(p)[0] for p in PP_ORDER], K.nvol, nfreq))
                self.assertEqual(set(table), set(BUILD_ROWS + APPLY_ROWS))
                if shape == (12, 6, 1):
                    self.assertGreaterEqual(table["matvec_workspace"],
                                            4 * table["call_overhead"])
                row = sum(table[r] for r in BUILD_ROWS)
                _log_ratio("kernel build B={} nfreq={} shape={}".format(B, nfreq, shape), peak,
                           row)
                self.assertLessEqual(peak, row)
                row = sum(table[r] for r in APPLY_ROWS)
                for p in PP_ORDER:
                    op = K.sector_op(p).operator
                    x = np.random.default_rng(p + 2).normal(size=op.shape[0]) + 0j
                    for how in ("matvec", "rmatvec"):
                        _, peak = _trace(lambda: getattr(op, how)(x))
                        _log_ratio("kernel {} B={} nfreq={} shape={} p={}".format(
                            how, B, nfreq, shape, p), peak, row)
                        self.assertLessEqual(peak, row, (p, how))


class TestPBWorkspace(unittest.TestCase):
    """Gate P-mem on the packed ``p_b`` into a caller's ``out``: its traced
    peak is at most one leading slice of the largest sector block (its
    workspace) plus the call overhead."""

    def test_p_b_packed(self):
        sm = _Chain().sm
        rng = np.random.default_rng(3)
        lead = (NMAT, FIXTURE_LZSZ_GRID.nvol)
        M = PackedArrays(sm, tuple(rng.normal(size=lead + (n, n)) + 0j for n in sm.sizes))
        out = PackedArrays.zeros(sm, lead)
        _, peak = _trace(lambda: p_b(M, sm.nd, out=out))
        row = 16 * lead[1] * max(sm.sizes) ** 2 + so_bond._CALL_OVERHEAD
        _log_ratio("p_b packed", peak, row)
        self.assertLessEqual(peak, row)


class TestEstimatorByHand(unittest.TestCase):
    """Construction check: one multi-sector case of the solver estimate
    evaluated by hand (dims 7 and 21, k = 1, sector vectors 16 and 44, a
    gap of 64 entries)."""

    def test_two_sectors(self):
        # the dense branch (dim 7 <= 20): 16 (7 * 49 + 3 * 16 * 7 + 4 * 7)
        dense = 16 * 707                                          # 11312
        # the ARPACK branch (dim 21 > 20, ncv = 20):
        # 16 (20 * 21 + 3 * 400 + 14 * 20 + 3 * 21) + 24 * 44 + 6 * 16 * 44
        # + 16 * 3 * 1 * (21 + 44), plus the retained sets: the spectral-radius solve of
        # spectral_shift="auto" (6 vectors and values) and the one eigenvalue of
        # the first which='LM' solve, 16 (6 * 22 + 1)
        arpack = 16 * 1963 + 1056 + 4224 + 3120 + 2128            # 41936
        # the null rows: 16 (1 * (7 + 16) + 64) and 16 (1 * (21 + 44) + 64)
        null = max(16 * 87, 16 * 129)
        worst = max(dense, arpack, null)
        # candidate and kept gaps 16 (1 + 2) 64, reduced vectors 16 (7 + 21)
        want = worst + 16 * 3 * 64 + 16 * 28 + 16 * 1024
        self.assertEqual(want, 61840)
        self.assertEqual(estimate_sector_solver_bytes([7, 21], 1, 64, n_sector_vec=[16, 44]),
                         want)


class TestEstimatorRefusals(unittest.TestCase):
    """Construction check: the estimators refuse each invalid argument
    (exact texts; the cases not listed with their formulas)."""

    def test_sector_solver(self):
        name = "estimate_sector_solver_bytes"
        for args, kw, msg in (
                ((np.zeros((1, 2), dtype=int), 1, 8), dict(n_sector_vec=[4, 4]),
                 "dims must be a non-empty sequence of integers >= 0, got "
                 + repr(np.zeros((1, 2), dtype=int))),
                (("2", 1, 8), dict(n_sector_vec=[4]),
                 "dims must be a non-empty sequence of integers >= 0, got '2'"),
                (([2], 1, 8), dict(n_sector_vec=[-4]),
                 "n_sector_vec must be a sequence of 1 integers >= 0 (one per dims entry), "
                 "got [-4]"),
                (([2], 1.5, 8), dict(n_sector_vec=[4]), "k must be an integer >= 1, got 1.5"),
                (([2], 1, True), dict(n_sector_vec=[4]),
                 "n_full_gap must be an integer >= sum(n_sector_vec) = 4, got True")):
            with self.subTest(msg=msg):
                with self.assertRaisesRegex(ValueError, exact("{}: {}".format(name, msg))):
                    estimate_sector_solver_bytes(*args, **kw)

    def test_kernel_workspace(self):
        name = "estimate_packed_kernel_workspace_bytes"
        for args, msg in (
                ((True, (1, 2, 1), 1, 1), "B must be an integer >= 1, got True"),
                ((1, (1, 2.0, 1), 1, 1),
                 "pp_pair_counts must be three positive integers, got (1, 2.0, 1)"),
                ((1, (1, 2, 1, 1), 1, 1),
                 "pp_pair_counts must be three positive integers, got (1, 2, 1, 1)"),
                ((1, (1, 2, 1), 1.0, 1), "nvol must be an integer >= 1, got 1.0"),
                ((1, (1, 2, 1), 1, 0), "nfreq must be an integer >= 1, got 0")):
            with self.subTest(msg=msg):
                with self.assertRaisesRegex(ValueError, exact("{}: {}".format(name, msg))):
                    estimate_packed_kernel_workspace_bytes(*args)


if __name__ == "__main__":
    unittest.main()
