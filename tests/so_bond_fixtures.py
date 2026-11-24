"""Shared test fixtures of the spin-orbital vertex tests that need the
production readers (``tests/so_bond_oracles.py`` stays free of solver
imports).

* :func:`onsite_decl` -- the Hermitian-closed on-site declaration of a type.
* :func:`rpa_spinful_gamma` -- the antisymmetrized spinful vertex of the
  spin-orbital RPA (``rpa.spinful_bubble_vertex`` of ``ham_inter_q`` and
  ``spinful_exchange()``), read from interaction files.
* :func:`production_gamma_b` -- the production ``Gamma_B`` of a table set
  through a vertex-only admission plan.
* :func:`green_kw`, :func:`pair_kernel` -- the Green-function layout of
  the bond bubble and the bond pairing kernel of code-slot vertex blocks.
* :func:`exact` -- an anchored regex of an exact refusal message.
* :func:`assert_recomputed` -- two computations of the same arithmetic
  agree to round-off (:data:`RECOMPUTE_RTOL`).
* :func:`chibar_static`, :func:`production_gamma_label`,
  :func:`production_first_order` -- the production spinful bond vertex in
  the label form of an exact-diagonalization fixture and its first-order
  correlator ``u^T [-chibar(q) Gamma(-q) chibar(q)] v``.
* :func:`canonical_green` -- the free Green function of an S_z-conserving
  fixture built per spin block (exact zeros between the spin blocks).
* :func:`static_pair_inputs`, :func:`static_pair_kernels` -- the static
  pairing vertex and pair propagator of an S_z-conserving fixture through
  the packed production chain, and its kernel, packed and dense.
* :func:`dense_spectrum`, :func:`cluster_projector`, :func:`match_spectrum`
  -- the eigenvalue matching rule: optimal assignment, cluster-level
  condition numbers from the left and right eigenvectors.
* the packed production chain of the conserving fixtures (``chain_*``,
  the bubbles, the dressing, the pairing parts and batch), the packed
  pairing kernel on random inputs (``kernel_setup``, ``random_kernels``),
  the FLEX ring of the ``P_B`` gate (``pb_*``) and sector operators of a
  given reduced matrix (:func:`subset_bases`, :func:`sector_op_of`),
  shared by the packed-layout tests and the index-map mutation tests.
"""
import contextlib
import functools
import logging
import os
import re
import shutil
import tempfile

from typing import NamedTuple

import numpy as np
import scipy.linalg
from scipy.optimize import linear_sum_assignment
from scipy.sparse.linalg import LinearOperator

import hwave.qlmsio.read_input_k as read_input_k
from hwave import sc
from hwave.solver import bond_channels, flex_bond, so_bond
from hwave.solver import eliashberg_bond as eb
from hwave.solver import rpa
from hwave.solver.so_bond import (SLOT_CONVENTION, PackedArrays, PackedBondPairKernel,
                                  SOBondCaps, plan_so_bond, select_layout, slot_geometry,
                                  so_bond_topology, so_bond_vertex, sz_sector_map)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import (BASE, CELL, FIXTURE_LZSZ, FIXTURE_LZSZ_GRID, FIXTURE_SZ_GRID,
                                   FIXTURE_SZ_GRID_4X2, FIXTURE_SZ_GRID_H0, FIXTURE_SZ_GRID_SYM,
                                   FIXTURE_SZ_RING, TYPES_OFFSITE, OracleGeometry, closed_row,
                                   expand_tables, file_gamma_label, first_order_prediction,
                                   gamma_code_to_label, lindhard_static, mesh_index,
                                   ph_sector_lists, pp_sector_lists, random_green, sz_of_mode)


@contextlib.contextmanager
def quiet():
    """Silence log records up to WARNING inside the block."""
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _wfile(path, head, n, rows, value_format):
    line = "  %d %d %d %d %d  " + value_format + "  " + value_format + "\n"
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(head, n))
        for r in rows:
            f.write(line % r)


def onsite_decl(t, v):
    """The Hermitian-closed on-site declaration of type ``t`` (two physical
    orbitals): the diagonal rows of CoulombIntra, the inter-orbital pair
    ``(0, 1)`` / ``(1, 0)`` of every other type."""
    if t == "CoulombIntra":
        return {t: {((0, 0, 0), (0, 0)): v, ((0, 0, 0), (1, 1)): 0.8 * v}}
    return {t: {((0, 0, 0), (0, 1)): v, ((0, 0, 0), (1, 0)): np.conj(v)}}


def rpa_spinful_gamma(decl, norb, shape, *, beta=2.0, value_format="%.17g",
                      orbital_site=None):
    """``Gamma (nvol, nd**2, nd**2)`` (``nd = 2 * norb``) of the spin-orbital
    RPA built from the declarations ``{type: {((Rx, Ry, Rz), (a, b)): v}}``
    (physical orbital indices, as the interaction files are read in
    spin-orbital mode) on the mesh ``shape``. ``value_format`` is the printf
    format of the coefficients in the files; ``orbital_site`` is passed to
    the solver parameters when given."""
    nd = 2 * norb
    d = tempfile.mkdtemp(prefix="so_vertex_reader_")
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % nd + "0 0 0\n" * nd)
        _wfile(os.path.join(d, "transfer.dat"), "Transfer", nd,
               [(0, 0, 0, p + 1, p + 1, 0.1, 0.0) for p in range(nd)], value_format)
        inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for t, tbl in decl.items():
            _wfile(os.path.join(d, t + ".dat"), t, norb,
                   [(R[0], R[1], R[2], a + 1, b + 1, complex(v).real, complex(v).imag)
                    for (R, (a, b)), v in tbl.items()], value_format)
            inter[t] = t + ".dat"
        io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
        par = {"T": 1.0 / beta, "mu": 0.0, "CellShape": list(shape),
               "SubShape": [1, 1, 1], "Nmat": 4}
        if orbital_site is not None:
            par["orbital_site"] = orbital_site
        with quiet():
            s = rpa.RPA(io.get_param("ham"), {}, {
                "mode": "RPA", "param": par, "enable_spin_orbital": True,
                "calc_scheme": "general"})
        nvol = int(np.prod(shape))
        ham = np.asarray(s.ham_info.ham_inter_q).reshape(nvol, nd, nd, nd, nd)
        return rpa.spinful_bubble_vertex(ham, s.ham_info.spinful_exchange())
    finally:
        shutil.rmtree(d, ignore_errors=True)


def exact(msg):
    """The regex that matches exactly ``msg``."""
    return "^" + re.escape(msg) + "$"


#: two computations of the same arithmetic on different buffers agree to this
#: (relative, max-norm): the SIMD loops of numpy < 2 depend on the memory
#: alignment, so a recomputed float array is not reproducible bitwise
RECOMPUTE_RTOL = 1.0e-15


def assert_recomputed(tc, got, want, msg=None, rtol=RECOMPUTE_RTOL):
    """``max|got - want| <= rtol max|want|`` for arrays or ``PackedArrays``
    (block by block against the largest block entry), shapes equal. Use a
    bitwise compare only where an array is moved, not recomputed."""
    if isinstance(want, PackedArrays):
        tc.assertIsInstance(got, PackedArrays, msg)
        pairs = list(zip(got.blocks, want.blocks))
        tc.assertEqual(len(got.blocks), len(want.blocks), msg)
    else:
        pairs = [(np.asarray(got), np.asarray(want))]
    scale = max(float(np.max(np.abs(w))) if w.size else 0.0 for _, w in pairs)
    for g, w in pairs:
        tc.assertEqual(g.shape, w.shape, msg)
        if w.size:
            tc.assertLessEqual(float(np.max(np.abs(g - w))), rtol * scale, msg)


def production_gamma_b(tables, norb, topo, shape):
    """The production ``Gamma_B`` ``(nvol, ND, ND)`` of the interaction
    ``tables`` (``norb`` physical orbitals) on the channels of ``topo`` (any
    object with ``delta_r`` and ``reverse``) and the mesh ``shape``, built
    through a vertex-only admission plan."""
    c = compile_interaction(tables, norb)
    plan = plan_so_bond(nmat=2, nvol=int(np.prod(shape)), B=len(topo.delta_r), nd=2 * norb,
                        consumers=("vertex",), caps=SOBondCaps(None, None, None))
    return so_bond_vertex(c, topo, shape, plan=plan)


def green_kw(g):
    """sc layout ``(nd, nd, Nx, 1, 1, nmat)`` -> the bond-bubble layout
    ``(1, nmat, Nx, nd, nd)`` (a chain along x)."""
    return np.ascontiguousarray(g[:, :, :, 0, 0, :].transpose(3, 2, 0, 1))[None]


def pair_kernel(blocks, G2, view, nd, shape, beta):
    """The bond pairing kernel of the code-slot vertex blocks ``(nfreq, nvol,
    B nd^2, B nd^2)`` on the uniform bosonic grid of the pair bubble ``G2``
    ``(nd, nd, nd, nd, Nx, Ny, Nz, nfreq)``, no instantaneous vertex (a
    frequency-independent part sits in every block)."""
    nvol = int(np.prod(shape))
    B = int(view.n_channels)
    nfreq = int(blocks.shape[0])
    vertex = eb.PairVertexUniform(eb.ArrayBlockSource({"V": blocks}, nd * nd), "V", B,
                                  nd * nd, nfreq, nvol)
    return eb.BondPairKernel(vertex, G2, view, xp=np, spatial_shape=shape, norb=nd,
                             beta=beta, nfreq=nfreq, V_inst=None, residency="host",
                             host_cap=np.inf)


@functools.lru_cache(maxsize=None)
def chibar_static(fixture):
    """The analytic static label bubble of a fixture's band (read-only
    cache by content)."""
    return lindhard_static(*fixture.band, fixture.geom)


@functools.lru_cache(maxsize=None)
def production_gamma_label(fixture, frame):
    """Production ``Gamma_B`` of ``fixture``'s tables, mapped to the label
    form by ``gamma_code_to_label`` on the production channels and embedded
    into the label geometry of ``frame`` (a fixture with the same mesh and
    band; channels matched by displacement modulo the mesh)."""
    nd, shape = fixture.nd, fixture.spatial_shape
    tables = fixture.tables_dict
    topo = so_bond_topology(compile_interaction(tables, fixture.norb_phys), CELL)
    gb = production_gamma_b(tables, fixture.norb_phys, topo, shape)
    geom_p = slot_geometry(topo, shape)
    gl_p = gamma_code_to_label(gb, geom_p, nd, SLOT_CONVENTION.slot_sign)
    geom = frame.geom
    chan = {int(mesh_index(np.asarray(d), shape)): m for m, d in enumerate(geom.delta_r)}
    npair = nd * nd
    idx = np.concatenate([chan[int(mesh_index(np.asarray(d), shape))] * npair + np.arange(npair)
                          for d in geom_p.delta_r])
    NL = len(geom.delta_r) * npair
    out = np.zeros((fixture.nvol, NL, NL), dtype=np.complex128)
    out[:, idx[:, None], idx[None, :]] = gl_p
    out.flags.writeable = False
    return out


def production_first_order(fixture, frame, sketches, qs):
    """``u^T [-chibar(q) Gamma(-q) chibar(q)] v`` ``(len(qs), len(sketches))``
    of the production vertex of ``fixture``'s tables on the labels and band
    of ``frame``."""
    return first_order_prediction(chibar_static(frame), production_gamma_label(fixture, frame),
                                  frame.spatial_shape, sketches, qs)


@functools.lru_cache(maxsize=None)
def canonical_green(fixture, nmat):
    """``G(k, i w_n) = [i w_n + mu - H_s(k)]^{-1}`` ``(nmat, nvol, nd, nd)``
    (spin-block index, ``w_n = (2n + 1 - nmat) pi T``) of an S_z-conserving
    fixture, inverted per spin block ``s``: the entries between the spin
    blocks are exactly zero. Refuses a band with a nonzero spin-off-diagonal
    entry. Read-only (cached by content)."""
    h1_k, mu, beta = fixture.band
    h1_k = np.asarray(h1_k, dtype=np.complex128)
    nvol, nd, _ = h1_k.shape
    norb = nd // 2
    up, dn = slice(0, norb), slice(norb, nd)
    if np.any(h1_k[:, up, dn] != 0) or np.any(h1_k[:, dn, up] != 0):
        raise ValueError("canonical_green: the band has spin-off-diagonal entries")
    w = (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
    eye = np.eye(norb)
    out = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
    for blk in (up, dn):
        A = (1j * w + mu)[:, None, None, None] * eye - h1_k[None, :, blk, blk]
        out[:, :, blk, blk] = np.linalg.inv(A)
    out.flags.writeable = False
    return out


@functools.lru_cache(maxsize=None)
def static_pair_inputs(fixture, nmat=8):
    """``(V, G2, geom, sector_map, plan)``: the static pairing vertex of an
    S_z-conserving fixture through the packed production chain -- the packed
    ``Gamma_B`` (``so_bond_vertex``), the static packed bubble of
    :func:`canonical_green` (``so_bond_bubble``), its dressing
    (``so_bond_dress``, a batch of one at ``l0 = nmat // 2``) and ``V =
    V_bare + V_fluct`` (``so_bond_pairing_vertex_parts``) as
    ``PackedArrays`` of leading shape ``(1, nvol)`` -- with the pair
    propagator ``sc._calc_g2`` of the same Green function, the slot
    geometry, the sector map and a packed vertex plan. Cached by content."""
    nd, shape, nvol = fixture.nd, fixture.spatial_shape, fixture.nvol
    c = compile_interaction(fixture.tables_dict, fixture.norb_phys)
    topo = so_bond_topology(c, CELL)
    geom = slot_geometry(topo, shape)
    layout, reason = select_layout(c, fixture.transfer_dict, override="auto")
    if layout != "sz_packed":
        raise ValueError("static_pair_inputs: the fixture does not conserve S_z ({})"
                         .format(reason))

    def plan(consumer):
        return plan_so_bond(nmat=nmat, nvol=nvol, B=geom.B, nd=nd, consumers=(consumer,),
                            caps=SOBondCaps(None, None, None), layout=layout,
                            layout_reason=reason)
    sm = sz_sector_map(topo, nd)
    gamma = so_bond_vertex(c, topo, shape, plan=plan("vertex"))
    green = canonical_green(fixture, nmat)
    chibar = PackedArrays.zeros(sm, (nvol,))
    so_bond.so_bond_bubble(green[None].copy(), None, fixture.beta, topo, shape,
                           plan=plan("chi_static"), freq="static", sink=chibar)
    chi = so_bond.so_bond_dress(PackedArrays(sm, tuple(b[None] for b in chibar.blocks)), gamma,
                                plan=plan("chi_static"), l0=nmat // 2, spatial_shape=shape)
    chi = PackedArrays(sm, tuple(b[0] for b in chi.blocks))
    v_bare, v_fluct = so_bond.so_bond_pairing_vertex_parts(gamma, chi, geom,
                                                           plan=plan("pairing_static"))
    V = PackedArrays(sm, tuple((a + b)[None] for a, b in zip(v_bare.blocks, v_fluct.blocks)))
    g_kw = green.reshape((nmat,) + shape + (nd, nd)).transpose(4, 5, 1, 2, 3, 0)
    G2 = sc._calc_g2(np.ascontiguousarray(g_kw), fixture.beta)[..., None]
    return V, G2, geom, sm, plan("vertex")


@functools.lru_cache(maxsize=None)
def static_pair_kernels(fixture, nmat=8):
    """``(packed, dense)``: the static pairing kernel of
    :func:`static_pair_inputs`. ``packed`` is the ``PackedBondPairKernel``
    of ``V``, ``dense`` the production ``BondPairKernel`` of
    ``V.to_dense()`` (:func:`pair_kernel`). Cached by content."""
    V, G2, geom, sm, plan = static_pair_inputs(fixture, nmat)
    shape = fixture.spatial_shape
    packed = PackedBondPairKernel(V, G2, geom, sm, plan=plan, spatial_shape=shape,
                                  beta=fixture.beta, nfreq=1)
    dense = pair_kernel(V.to_dense(), G2, geom, fixture.nd, shape, fixture.beta)
    return packed, dense


#: the eigenvalue matching rule: ``|lambda - mu| <= MATCH_TOL max(1, max|lambda|)
#: max(1, kappa_c)``; clusters join at ``MATCH_TOL max(1, max|lambda|)``
MATCH_TOL = 1.0e-8
#: the largest cluster condition number a fixture may have
KAPPA_MAX = 1.0e4


class DenseSpectrum(NamedTuple):
    """The complete spectrum of a dense matrix: eigenvalues ``vals``, right
    and left eigenvectors (columns, ``Y^H M = Lambda Y^H``), the clusters
    (``bond_channels.cluster_eigenvalues`` at ``MATCH_TOL scale``), the
    cluster condition number ``kappa`` of each eigenvalue's cluster and
    ``scale = max(1, max|lambda|)``."""
    vals: np.ndarray
    right: np.ndarray
    left: np.ndarray
    clusters: tuple
    kappa: np.ndarray
    scale: float


def cluster_projector(right, left, idx):
    """The spectral projector ``X_c (Y_c^H X_c)^{-1} Y_c^H`` of the
    eigenvectors ``idx`` (right ``X``, left ``Y``)."""
    X, Y = right[:, list(idx)], left[:, list(idx)]
    return X @ np.linalg.solve(Y.conj().T @ X, Y.conj().T)


def dense_spectrum(M):
    """:class:`DenseSpectrum` of the dense matrix ``M`` (``scipy.linalg.eig``
    with the left eigenvectors); ``kappa`` of a cluster is the 2-norm of its
    spectral projector (``||x|| ||y|| / |y^H x|`` for a simple
    eigenvalue)."""
    vals, left, right = scipy.linalg.eig(np.asarray(M), left=True, right=True)
    scale = max(1.0, float(np.abs(vals).max())) if vals.size else 1.0
    clusters = tuple(tuple(int(i) for i in c)
                     for c in bond_channels.cluster_eigenvalues(vals, deg_tol=MATCH_TOL * scale))
    kappa = np.empty(len(vals))
    for c in clusters:
        kappa[list(c)] = np.linalg.norm(cluster_projector(right, left, c), 2)
    return DenseSpectrum(vals, right, left, clusters, kappa, scale)


def match_spectrum(ref, mu):
    """Assign each value of ``mu`` to a distinct eigenvalue of the
    :class:`DenseSpectrum` ``ref`` (``linear_sum_assignment`` on ``|lambda -
    mu|``, complex). Returns ``(cols, ratio)``: ``cols[j]`` the eigenvalue of
    ``ref`` assigned to ``mu[j]``, and the largest ``|lambda - mu| /
    (MATCH_TOL scale max(1, kappa_c))`` (the rule holds when ``ratio <=
    1``)."""
    mu = np.asarray(mu, dtype=complex).ravel()
    if mu.size > ref.vals.size:
        raise ValueError("match_spectrum: {} values against {} eigenvalues".format(
            mu.size, ref.vals.size))
    cost = np.abs(mu[:, None] - ref.vals[None, :])
    rows, cols = linear_sum_assignment(cost)
    out = np.empty(mu.size, dtype=np.int64)
    out[rows] = cols
    bound = MATCH_TOL * ref.scale * np.maximum(1.0, ref.kappa[out])
    ratio = float(np.max(cost[np.arange(mu.size), out] / bound)) if mu.size else 0.0
    return out, ratio


# ---------------------------------------------------------------------------
# the S_z-packed production chain of the conserving fixtures (shared by the
# packed-layout tests and the index-map mutation tests)
# ---------------------------------------------------------------------------

NOCAP = SOBondCaps(None, None, None)
#: the subtest labels of the conserving fixtures
FIXTURE_NAMES = {FIXTURE_SZ_RING: "ring", FIXTURE_SZ_GRID: "grid", FIXTURE_LZSZ: "lzsz",
                 FIXTURE_SZ_GRID_H0: "grid h=0", FIXTURE_SZ_GRID_SYM: "grid spin-symmetric",
                 FIXTURE_LZSZ_GRID: "lzsz grid", FIXTURE_SZ_GRID_4X2: "grid 4x2"}


@functools.lru_cache(maxsize=None)
def chain_setup(fx):
    """``(compiled, topo, geom)`` of a fixture's interaction."""
    c = compile_interaction(fx.tables_dict, fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    return c, topo, slot_geometry(topo, fx.spatial_shape)


def chain_plan(fx, layout, consumers=("vertex",), nmat=2, **kw):
    """The admission plan of a fixture in ``layout`` (``"sz_packed"`` takes
    the selection's reason line and asserts the selection)."""
    c, topo, _ = chain_setup(fx)
    reason = None
    if layout == "sz_packed":
        selected, reason = select_layout(c, fx.transfer_dict, override="auto")
        assert selected == "sz_packed", (FIXTURE_NAMES.get(fx), reason)
    return plan_so_bond(nmat=nmat, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=consumers, caps=NOCAP, layout=layout,
                        layout_reason=reason, **kw)


@functools.lru_cache(maxsize=None)
def chain_gamma(fx, layout):
    """The production ``Gamma_B`` of a fixture in ``layout``."""
    c, topo, _ = chain_setup(fx)
    return so_bond_vertex(c, topo, fx.spatial_shape, plan=chain_plan(fx, layout))


def chain_lists(fx):
    """The test-side ph sector lists of a fixture's channels."""
    return ph_sector_lists(chain_setup(fx)[2].B, fx.nd)


def chain_sector_map(fx):
    return sz_sector_map(chain_setup(fx)[2], fx.nd)


def random_packed(sm, lead, seed):
    """Seeded complex random ``PackedArrays`` of leading shape ``lead`` (an
    int is one leading axis)."""
    lead = (lead,) if isinstance(lead, (int, np.integer)) else tuple(lead)
    rng = np.random.default_rng(seed)
    return PackedArrays(sm, tuple(rng.normal(size=lead + (n, n))
                                  + 1j * rng.normal(size=lead + (n, n)) for n in sm.sizes))


def copy_packed(p):
    return PackedArrays(p.sector_map, tuple(b.copy() for b in p.blocks))


#: Matsubara counts of the bubble fixtures
BUBBLE_NMAT = {FIXTURE_SZ_GRID: 8, FIXTURE_LZSZ: 8}


def bubble_green(fx):
    """The bubble input ``(1, nmat, nvol, nd, nd)``: the per-spin-block G."""
    return canonical_green(fx, BUBBLE_NMAT[fx])[None].copy()


def static_bubble(fx, layout, g, sink=None):
    """The production static bubble of ``g`` in ``layout`` (into ``sink``)."""
    topo = chain_setup(fx)[1]
    plan = chain_plan(fx, layout, ("chi_static",), BUBBLE_NMAT[fx])
    ND = plan.ND
    if sink is None:
        sink = (np.empty((fx.nvol, ND, ND), dtype=np.complex128) if layout == "dense"
                else PackedArrays.zeros(chain_sector_map(fx), (fx.nvol,)))
    out = so_bond.so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan,
                                 freq="static", sink=sink)
    assert out is sink
    return sink


def every_slice_bubble(fx, layout, g):
    """``(chibar, pairs)``: the whole slot (dense array or a copy of the
    packed blocks) and ``get_pair_full`` of every channel pair."""
    topo = chain_setup(fx)[1]
    plan = chain_plan(fx, layout, ("chi_dynamic",), BUBBLE_NMAT[fx])
    sm = chain_sector_map(fx)
    with so_bond.open_store(plan, sm if layout == "sz_packed" else None) as store:
        self_out = so_bond.so_bond_bubble(g, None, fx.beta, topo, fx.spatial_shape, plan=plan,
                                          freq="all", sink=store)
        assert self_out is None
        if layout == "dense":
            chibar = store.get_freq_batch("chibar", 0, plan.nmat).copy()
        else:
            chibar = PackedArrays(sm, tuple(b.copy() for b in store.packed("chibar").blocks))
        pairs = {(a, b): store.get_pair_full("chibar", a, b)
                 for a in range(plan.B) for b in range(plan.B)}
    return chibar, pairs


@functools.lru_cache(maxsize=None)
def production_bubbles(fx, freq):
    """``(dense, packed)`` production bubbles of the canonical G
    (``freq = "static"`` or ``"all"``)."""
    run = static_bubble if freq == "static" else every_slice_bubble
    return run(fx, "dense", bubble_green(fx)), run(fx, "sz_packed", bubble_green(fx))


#: the batch rows of the frequency-batch comparison (grid index l0 = 2, nb = 4)
BATCH = slice(2, 6)
#: relative size and seed of the sector-preserving perturbation of the bubbles
PERTURB = (0.05, 2711)


def packed_lead(p, lead):
    """The packed blocks of ``p`` indexed by ``lead`` along the leading axes (copies)."""
    return PackedArrays(p.sector_map, tuple(np.array(b[lead]) for b in p.blocks))


def dress_inputs(fx, kind, scale=1.0):
    """``(X_dense, X_packed, gamma_dense, gamma_packed, l0)``: a batch of the
    packed production bubbles of the canonical G (``"static"``: the static
    bubble as a batch of one; ``"batch"``: four slices of the every-slice
    bubble), each sector block plus a seeded random block of
    ``PERTURB[0]`` times its largest entry (this breaks the lattice and
    +1/-1 symmetries that make blocks tie in the guard's argmin, and keeps
    the input sector diagonal), its dense image, and the production vertex
    times ``scale``."""
    if kind == "static":
        p, l0 = packed_lead(production_bubbles(fx, "static")[1], np.s_[None]), BUBBLE_NMAT[fx] // 2
    else:
        p, l0 = packed_lead(production_bubbles(fx, "all")[1][0], BATCH), BATCH.start
    rng = np.random.default_rng(PERTURB[1])
    for b in p.blocks:
        b += PERTURB[0] * np.abs(b).max() * (rng.normal(size=b.shape)
                                             + 1j * rng.normal(size=b.shape))
    gd, gp = chain_gamma(fx, "dense"), chain_gamma(fx, "sz_packed")
    gp = PackedArrays(gp.sector_map, tuple(scale * b for b in gp.blocks))
    return p.to_dense(), p, scale * gd, gp, l0


def dress_plan(fx, layout, kind):
    return chain_plan(fx, layout, ("chi_static",) if kind == "static" else ("chi_dynamic",),
                      BUBBLE_NMAT[fx])


def dress(fx, layout, X, gamma, l0, kind, **kw):
    """``(chi, diagnostics)`` of the production dressing in ``layout``."""
    diag = {}
    chi = so_bond.so_bond_dress(X, gamma, plan=dress_plan(fx, layout, kind), l0=l0,
                                spatial_shape=fx.spatial_shape, diagnostics=diag, **kw)
    return chi, diag


@functools.lru_cache(maxsize=None)
def dressed(fx, kind):
    """``(X_dense, X_packed, gamma_dense, gamma_packed, chi_dense,
    chi_packed)``: the dressing inputs of :func:`dress_inputs` and the
    production dressing in both layouts (static: the batch of one taken
    at its single slice)."""
    xd, xp, gd, gp, l0 = dress_inputs(fx, kind)
    chi_d, _ = dress(fx, "dense", xd, gd, l0, kind)
    chi_p, _ = dress(fx, "sz_packed", xp, gp, l0, kind)
    if kind == "static":
        xd, xp, chi_d, chi_p = xd[0], packed_lead(xp, 0), chi_d[0], packed_lead(chi_p, 0)
    return xd, xp, gd, gp, chi_d, chi_p


def parts_plan(fx, layout):
    return chain_plan(fx, layout, ("pairing_static",), BUBBLE_NMAT[fx])


def batch_plan(fx, layout, **kw):
    return chain_plan(fx, layout, ("chi_dynamic",), BUBBLE_NMAT[fx], vertex_batch=True, **kw)


def pairing_parts(fx, layout, gamma, chi):
    """``(V_bare, V_fluct)`` of the production static pairing parts."""
    return so_bond.so_bond_pairing_vertex_parts(gamma, chi, chain_setup(fx)[2],
                                                plan=parts_plan(fx, layout))


def vertex_batch(fx, layout, gamma, Y, out=None, plan=None):
    """The production pairing vertex batch written into ``out`` (a fresh
    array or ``PackedArrays`` when None)."""
    if out is None:
        out = (np.empty(Y.shape, dtype=np.complex128) if layout == "dense"
               else PackedArrays.zeros(Y.sector_map, Y.lead_shape))
    plan = batch_plan(fx, layout) if plan is None else plan
    got = so_bond.so_bond_pairing_vertex_batch(gamma, Y, chain_setup(fx)[2], plan=plan,
                                               out=out)
    assert got is out
    return out


# -- the packed pairing kernel on random inputs ------------------------------

@functools.lru_cache(maxsize=None)
def kernel_setup(fx):
    """``(geom, sector_map, plan)`` of a fixture (packed vertex plan)."""
    c, topo, geom = chain_setup(fx)
    layout, reason = select_layout(c, fx.transfer_dict, override="auto")
    assert layout == "sz_packed", (FIXTURE_NAMES[fx], reason)
    plan = plan_so_bond(nmat=2, nvol=fx.nvol, B=geom.B, nd=fx.nd, consumers=("vertex",),
                        caps=NOCAP, layout=layout, layout_reason=reason)
    return geom, sz_sector_map(topo, fx.nd), plan


def random_vertex(sm, nfreq, nvol, seed):
    """Seeded random packed vertex of leading shape ``(nfreq, nvol)``."""
    rng = np.random.default_rng(seed)
    return PackedArrays(sm, tuple(rng.normal(size=(nfreq, nvol, n, n))
                                  + 1j * rng.normal(size=(nfreq, nvol, n, n))
                                  for n in sm.sizes))


def pp_of_pair(nd):
    """``p = S_z(b) + S_z(c)`` of every pair ``(b, c)`` (test side)."""
    norb = nd // 2
    return np.array([[int(round(sz_of_mode(b, norb) + sz_of_mode(c, norb))) for c in range(nd)]
                     for b in range(nd)])


def random_g2(fx, nfreq, seed):
    """A random pair propagator ``G2[i, l, j, m]`` (output pair ``(i, j)``,
    input pair ``(l, m)``) that conserves the pair S_z exactly."""
    nd = fx.nd
    rng = np.random.default_rng(seed)
    shape = (nd,) * 4 + fx.spatial_shape + (nfreq,)
    g2 = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    pp = pp_of_pair(nd)
    off = pp[:, None, :, None] != pp[None, :, None, :]
    g2[off] = 0
    return g2


def gap_shape(fx, nfreq):
    return (fx.nd, fx.nd) + fx.spatial_shape + ((nfreq,) if nfreq > 1 else ())


def pp_flat_index(nd, grid):
    """Per pp sector ``(+1, 0, -1)`` the ascending flat gap indices of a gap
    ``(nd, nd) + grid points`` (test side)."""
    return [(l[:, None] * grid + np.arange(grid)[None, :]).ravel() for l in pp_sector_lists(nd)]


@functools.lru_cache(maxsize=None)
def random_kernels(fx, nfreq, seed=5):
    """``(packed, dense, V, G2)`` on a random vertex and pair propagator."""
    geom, sm, plan = kernel_setup(fx)
    V = random_vertex(sm, nfreq, fx.nvol, seed)
    G2 = random_g2(fx, nfreq, seed + 1)
    packed = PackedBondPairKernel(V, G2, geom, sm, plan=plan, spatial_shape=fx.spatial_shape,
                                  beta=fx.beta, nfreq=nfreq)
    dense = pair_kernel(V.to_dense(), G2, geom, fx.nd, fx.spatial_shape, fx.beta)
    return packed, dense, V, G2


def cvec(rng, n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


# -- the FLEX ring of the P_B gate (a chain of L = 5 sites, norb_phys = 1) -----

PB_L, PB_NMAT, PB_BETA, PB_ND = 5, 4, 1.7, 2
PB_SHAPE = (PB_L, 1, 1)


def _pb_tables(nnn, pairlift):
    """CoulombIntra on site and every off-site type on the nearest bond (and
    at 0.47 of it on the next-nearest bond), PairHop complex."""
    out = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"]}}
    for t in TYPES_OFFSITE:
        if t == "PairLift" and not pairlift:
            continue
        v = BASE[t] * (1 + 0.6j if t == "PairHop" else 1.0)
        out[t] = dict(closed_row(1, 0, 0, v))
        if nnn:
            out[t].update(closed_row(2, 0, 0, 0.47 * v))
    return out


#: all seven types; the S_z-conserving subset (PairLift dropped)
PB_TABLES_ALL = {nnn: _pb_tables(nnn, True) for nnn in (False, True)}
PB_TABLES_SZ = {nnn: _pb_tables(nnn, False) for nnn in (False, True)}


def _pb_label_geometry(nnn):
    """The oracle's labels, listed independently of the production channels
    (another order; the functional does not depend on it)."""
    dr = [(1, 0, 0), (0, 0, 0), (-1, 0, 0)] + ([(-2, 0, 0), (2, 0, 0)] if nnn else [])
    rev = [dr.index(tuple(-x for x in d)) for d in dr]
    return OracleGeometry(dr, rev, PB_SHAPE)


PB_GEOM_L = {nnn: _pb_label_geometry(nnn) for nnn in (False, True)}


def pb_gamma_label(tables, geom):
    return file_gamma_label(expand_tables(tables, 1), geom, PB_ND)


def pb_green(seed, sz):
    """``(nmat, nvol, nd, nd)``; with ``sz`` the spin-off-diagonal entries
    are set to exactly zero."""
    g = random_green(seed, PB_NMAT, PB_L, PB_ND)
    if sz:
        g[..., 0, 1] = 0
        g[..., 1, 0] = 0
    return g


@functools.lru_cache(maxsize=None)
def pb_setup(nnn, sz):
    """``(compiled, topo, geom)`` of the P_B chain."""
    tables = (PB_TABLES_SZ if sz else PB_TABLES_ALL)[nnn]
    c = compile_interaction(tables, 1)
    topo = so_bond_topology(c, CELL)
    return c, topo, slot_geometry(topo, PB_SHAPE)


def pb_plan(nnn, sz, layout):
    c, topo, _ = pb_setup(nnn, sz)
    reason = None
    if layout == "sz_packed":
        selected, reason = select_layout(c, {}, override="auto")
        assert selected == "sz_packed", reason
    return plan_so_bond(nmat=PB_NMAT, nvol=PB_L, B=len(topo.delta_r), nd=PB_ND,
                        consumers=("chi_dynamic",), caps=NOCAP, layout=layout,
                        layout_reason=reason)


def pb_sigma(source, green, geom):
    """The production self-energy of the ``W`` slot of ``source``."""
    return flex_bond.calc_self_energy_bond(source, green[None].copy(), PB_BETA, geom, PB_SHAPE,
                                           PB_ND, 1)[0]


def pb_sigma_prod_packed(nnn, seed):
    """The production self-energy through the S_z-packed chain: packed
    Gamma_B, the bubble in a packed store, the ring per sector, ``p_b`` of
    the packed ring into the store's ``W`` slot, read by the self-energy."""
    c, topo, geom = pb_setup(nnn, True)
    plan = pb_plan(nnn, True, "sz_packed")
    sm = sz_sector_map(topo, PB_ND)
    gam = so_bond_vertex(c, topo, PB_SHAPE, plan=plan)
    assert isinstance(gam, PackedArrays)
    g = pb_green(seed, True)
    with so_bond.PackedBondBlockStore(PB_NMAT, PB_L, sm, ("chibar", "W")) as st:
        so_bond.so_bond_bubble(g[None].copy(), None, PB_BETA, topo, PB_SHAPE, plan=plan,
                               freq="all", sink=st)
        X = st.packed("chibar")
        M = PackedArrays(sm, tuple(0.5 * gb[None] @ xb @ gb[None]
                                   for gb, xb in zip(gam.blocks, X.blocks)))
        W = st.packed("W")
        assert so_bond.p_b(M, PB_ND, out=W) is W
        return pb_sigma(st, g, geom)


# -- sector operators of a given reduced matrix ------------------------------

def subset_bases(shape, dims):
    """Subset bases of the requested dimensions on disjoint orbits of the
    antisymmetric basis of ``shape``."""
    full = sc._SOAntisymmetricBasis(shape)
    out, start = [], 0
    for d in dims:
        reps = full.rep[start:start + d]
        start += d
        idx = np.sort(np.concatenate([reps, full.partner[reps]])).astype(np.int64)
        out.append(sc._SOAntisymmetricBasis.from_subset(shape, idx))
    assert start <= full.dim, (shape, dims)
    return out


def sector_op_of(p, basis, M, *, adjoint=True, matmat=True):
    """A ``SectorOp`` whose reduced operator is ``M``: ``K = E M E^dagger``
    on the sector vector, with its adjoint precomputed when ``adjoint``;
    ``matmat = False`` gives forward and adjoint products only (as the
    kernel's ``sector_op``)."""
    E = basis.to_full(np.eye(basis.dim, dtype=complex))
    K = E @ np.asarray(M, dtype=complex) @ E.conj().T
    n = basis.n
    kw = dict(matvec=lambda x: K @ np.asarray(x).ravel())
    if matmat:
        kw["matmat"] = lambda X: K @ X
    if adjoint:
        KH = np.ascontiguousarray(K.conj().T)
        kw["rmatvec"] = lambda y: KH @ np.asarray(y).ravel()
    return sc.SectorOp(p, LinearOperator((n, n), dtype=complex, **kw), basis, adjoint)


# ---------------------------------------------------------------------------
# Second order: grids, the on-site control and the Green functions
# ---------------------------------------------------------------------------

from tests.so_bond_oracles import FIXTURES, fixture_name, free_green, rel  # noqa: E402,F401

#: Matsubara grids of the fast and the full-grid gates
C1A_NMAT_FAST = 8
C1A_NMAT_FULL = 64
#: seed of the random (non-factorized) Green function ``"rand"``
C1A_SEED = 2701
#: the on-site control: FIXTURE_I's band with CoulombIntra only, on the
#: single-channel (B = 1) label geometry
C1A_ONSITE_TABLES = {"CoulombIntra": {((0, 0, 0), (0, 0)): 0.83}}
C1A_ONSITE_GEOM = OracleGeometry([(0, 0, 0)], [0], (5, 1, 1))

_C1A_KINDS = ("band", "rand", "herm")


@functools.lru_cache(maxsize=None)
def _c1a_green_cached(name, kind, nmat):
    fx = FIXTURES[name]
    if kind == "band":
        h1_k, mu, beta = fx.band
        g = free_green(h1_k, mu, beta, nmat)
    else:
        g = random_green(C1A_SEED, nmat, fx.nvol, fx.nd)
        if kind == "herm":
            g = 0.5 * (g + np.conj(np.swapaxes(g[::-1], -1, -2)))
    g = np.ascontiguousarray(g, dtype=np.complex128)
    g.flags.writeable = False
    return g


def c1a_green(fx, kind, nmat):
    """The Green function ``(nmat, nvol, nd, nd)`` of a pinned fixture (no
    tail): ``"band"`` = the free Green function of ``fx.band``, ``"rand"`` =
    ``random_green(C1A_SEED, ...)``, ``"herm"`` = ``1/2 [G_r(k, n) +
    G_r(k, nmat - 1 - n)^dagger]``. Each call returns a fresh copy."""
    name = fixture_name(fx)
    if name is None:
        raise ValueError("c1a_green: the fixture is not a pinned (named) fixture")
    if kind not in _C1A_KINDS:
        raise ValueError("c1a_green: unknown kind {!r} (expected one of {})".format(
            kind, _C1A_KINDS))
    return _c1a_green_cached(name, kind, int(nmat)).copy()


def win(nmat):
    """The inner eighth of a frequency grid of ``nmat`` points."""
    return slice(nmat // 2 - nmat // 16, nmat // 2 + nmat // 16)


import tests.so_bond_oracles as so_bond_oracles  # noqa: E402


def c1a_gamma_label(fx, tables, geom, **filing):
    """The label vertex ``Gamma(q)`` of ``tables`` on the labels of ``geom``
    (``file_gamma_label`` of the documented operators; ``filing`` selects a
    wrong filing). Looked up through the module attribute
    ``so_bond_oracles.file_gamma_label`` at call time."""
    return so_bond_oracles.file_gamma_label(
        so_bond_oracles.expand_tables(tables, fx.norb_phys), geom, fx.nd, **filing)


def c1a_onsite_flex_sigma(green, *, rpa_return):
    """The production second order ``(nmat, nvol, nd, nd)`` of one
    spin-orbital FLEX map (``calc_scheme = "general"``, ``flex_hartree_fock
    = true``) with the on-site control tables on FIXTURE_I's ring, at the
    Green function ``green`` ``(nmat, 5, 2, 2)``.

    The inputs (identity cell, FIXTURE_I's Transfer rows in 1-based
    interleaved indices, ``C1A_ONSITE_TABLES`` with norb 1) are written to a
    temporary directory that is removed before returning. The bubble is
    ``_calc_chi0q(G[None], 0, beta)``; ``_flex_compute_veff_spinful`` runs
    with ``_solve_rpa`` replaced by ``rpa_return(chi0q)`` (with the identity,
    ``chi := chibar`` and ``V = P[1/2 Gamma chibar Gamma]`` exactly, the
    second order); the result is ``_calc_self_energy_general(G[None], V,
    beta)[0]``. ``rpa_return`` is a required callable ``chi0q -> chi``."""
    import hwave.solver.flex as flex_mod
    from unittest import mock
    from tests.test_second_order_factors import _write_wan
    from tests.so_bond_oracles import FIXTURE_I
    fx = FIXTURE_I
    beta = fx.beta
    green = np.asarray(green, dtype=np.complex128)
    nmat = green.shape[0]
    d = tempfile.mkdtemp(prefix="c1a_onsite_")
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n2\n" + "0 0 0\n" * 2)
        rows = [(R[0], R[1], R[2], A + 1, B + 1, complex(v).real, complex(v).imag)
                for (R, (A, B)), v in sorted(fx.transfer_dict.items())]
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", 2, rows)
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for t, tbl in C1A_ONSITE_TABLES.items():
            crows = [(R[0], R[1], R[2], a + 1, b + 1, complex(v).real, complex(v).imag)
                     for (R, (a, b)), v in sorted(tbl.items())]
            _write_wan(os.path.join(d, t + ".dat"), t, 1, crows)
            idict[t] = t + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 1.0 / beta, "mu": fx.mu, "CellShape": list(fx.spatial_shape),
               "SubShape": [1, 1, 1], "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
               "flex_hartree_fock": True}
        with quiet():
            s = flex_mod.FLEX(r.get_param("ham"), {}, {
                "mode": "FLEX", "param": par, "enable_spin_orbital": True,
                "calc_scheme": "general"})
            s._calc_epsilon_k({})
            Gk = green[None]
            chi0 = s._calc_chi0q(Gk, np.zeros_like(Gk), beta)[0]
            with mock.patch.object(s, "_solve_rpa",
                                   side_effect=lambda chi0q, ham: rpa_return(chi0q)):
                _, v, _ = s._flex_compute_veff_spinful(chi0)
            return s._calc_self_energy_general(Gk, v, beta)[0]
    finally:
        shutil.rmtree(d, ignore_errors=True)


# ---------------------------------------------------------------------------
# The production chain of the second-order (and third-order) ring (gate C1a-P)
# ---------------------------------------------------------------------------

def _pair_fields(geom, nd):
    """``(I, m, x, y)`` of the flat slot ``I = m nd**2 + x nd + y``."""
    npair = nd * nd
    I = np.arange(len(geom.delta_r) * npair)
    return I, I // npair, (I % npair) // nd, I % nd


def _read(M, f):
    """``W[..., I, J] = M[..., f(I), f(J)]``."""
    return M[..., f[:, None], f]


def _reading_reversed(M, geom, nd):
    _, m, x, y = _pair_fields(geom, nd)
    return _read(M, np.asarray(geom.reverse)[m] * nd * nd + y * nd + x)


def _reading_reversed_phased(sign):
    def reading(M, geom, nd):
        _, m, _, _ = _pair_fields(geom, nd)
        ph = np.asarray(geom.phase)[:, m]                       # (nvol, ND)
        p = ph if sign > 0 else 1.0 / ph
        return _reading_reversed(M, geom, nd) * (p[:, :, None] * np.conj(p)[:, None, :])
    return reading


def _reading_channel0_identity(M, geom, nd):
    I, m, x, y = _pair_fields(geom, nd)
    return _read(M, np.where(m == 0, I, m * nd * nd + y * nd + x))


#: the pair readings ``(M, geom, nd) -> W`` of the ring ``M``
#: ``(nmat, nvol, ND, ND)`` on the slot geometry; ``"p_b"`` is accepted (the
#: module attribute ``so_bond.p_b`` is looked up at call time), the others
#: are rejected readings ``W[..., I, J] = M[..., f(I), f(J)]`` (the phased
#: ones times ``e^{+-iq.(dr_m(I) - dr_m(J))}``)
READINGS = {
    "p_b": lambda M, geom, nd: so_bond.p_b(M, nd),
    "reversed": _reading_reversed,
    "reversed_plus": _reading_reversed_phased(+1),
    "reversed_minus": _reading_reversed_phased(-1),
    "no_swap": lambda M, geom, nd: M,
    "channel0_identity": _reading_channel0_identity,
}


def c1a_sigma_prod(fx, tables, green, *, reading="p_b", ring_factor=0.5, ring_order=2):
    """The production self-energy ``(nmat, nvol, nd, nd)`` of the ring of
    ``tables`` on ``fx``'s mesh at the Green function ``green`` ``(nmat,
    nvol, nd, nd)`` (dense layout, no caps).

    ``Gamma_B`` (``so_bond_vertex``) and ``chibar_B`` (``so_bond_bubble``
    into a store, detached) give the ring ``M = ring_factor Gamma_B chibar_B
    Gamma_B`` (``ring_order = 2``), ``Gamma_B chibar_B Gamma_B chibar_B
    Gamma_B`` (``3``; coefficient +1, the NEGATIVE of the third-order part
    ``-Gamma chibar Gamma chibar Gamma`` of the solver's ``Gamma (chi -
    chibar/2) Gamma``) or its negative control with the first ``chibar_B``
    transposed (``"3T"``); ``W = READINGS[reading](M, geom,
    nd)`` on the slot geometry ``so_bond.slot_geometry(topo, shape)``
    (looked up at call time), and ``flex_bond.calc_self_energy_bond`` of
    ``W``. A geometry that is not a ``BondSlotGeometry`` with the channels
    of the topology is refused before the transport."""
    if ring_order not in (2, 3, "3T"):
        raise ValueError("c1a_sigma_prod: unknown ring_order {!r}".format(ring_order))
    if ring_order != 2 and ring_factor != 0.5:
        raise ValueError("c1a_sigma_prod: ring_factor applies to ring_order=2 only")
    nmat = green.shape[0]
    c = compile_interaction(tables, fx.norb_phys)
    topo = so_bond_topology(c, CELL)
    geom = so_bond.slot_geometry(topo, fx.spatial_shape)
    if not (isinstance(geom, so_bond.BondSlotGeometry)
            and np.array_equal(geom.delta_r, topo.delta_r)):
        raise ValueError("c1a_sigma_prod: the transport needs the slot geometry of the topology")
    plan = plan_so_bond(nmat=nmat, nvol=fx.nvol, B=len(topo.delta_r), nd=fx.nd,
                        consumers=("chi_dynamic",), caps=NOCAP, layout="dense")
    gamma_b = so_bond_vertex(c, topo, fx.spatial_shape, plan=plan)
    with flex_bond.BondBlockStore(nmat, fx.nvol, plan.ND, fx.nd ** 2, ("chibar",)) as st:
        so_bond.so_bond_bubble(green[None].copy(), None, fx.beta, topo, fx.spatial_shape,
                               plan=plan, freq="all", sink=st)
        chibar = st.detach("chibar")
    assert chibar.shape == (nmat, fx.nvol, plan.ND, plan.ND), chibar.shape
    g = gamma_b[None]
    if ring_order == 2:
        t1 = g @ chibar
        M = t1 @ g
        del t1
        M *= ring_factor
    else:
        t1 = g @ (chibar if ring_order == 3 else np.swapaxes(chibar, -2, -1))
        t2 = t1 @ g
        del t1
        t3 = t2 @ chibar
        del t2
        M = t3 @ g
        del t3
    del chibar
    W = READINGS[reading](M, geom, fx.nd)
    del M
    return flex_bond.calc_self_energy_bond(eb.ArrayBlockSource({"W": W}, fx.nd ** 2),
                                           green[None].copy(), fx.beta, geom, fx.spatial_shape,
                                           fx.nd, 1)[0]


def c1a_hf_param_ham(tables, spatial_shape):
    """A ``param_ham`` of ``tables`` for ``build_flex_hf_tables``:
    ``{type: {((Rx mod Nx, Ry mod Ny, Rz mod Nz), (a, b)): complex}}`` with
    0-based orbitals and the values of keys that coincide on the mesh summed.
    This is the structure ``QLMSkInput`` hands the solver except that the
    reader keeps the raw (possibly negative) displacements ``R`` while this
    helper reduces them modulo the mesh; ``build_flex_hf_tables`` accepts
    either and gives the same Hartree-Fock map."""
    shape = tuple(int(n) for n in spatial_shape)
    out = {}
    for itype, tbl in tables.items():
        rows = {}
        for (R, (a, b)), v in tbl.items():
            R = (int(R), 0, 0) if isinstance(R, (int, np.integer)) else tuple(int(x) for x in R)
            key = (tuple(r % n for r, n in zip(R, shape)), (int(a), int(b)))
            rows[key] = rows.get(key, 0.0) + complex(v)
        out[itype] = rows
    return out


def c1a_rho_r(rho_modes, spatial_shape, norb_phys):
    """The real-space spin-block density ``(nvol, nd, nd)`` of the
    interleaved-mode density ``rho_modes[p, q] = <c+_p c_q>``:
    ``rho_r[r, X, Y] = rho[ed_mode(r, X), ed_mode(0, Y)] = <c+_X(r) c_Y(0)>``,
    the convention of ``hartree_fock.accumulate_hf`` (``r`` the C-order cell
    index)."""
    rho = np.asarray(rho_modes)
    nvol = int(np.prod(spatial_shape))
    nd = 2 * int(norb_phys)
    X = np.arange(nd)
    rows = so_bond_oracles.ed_mode(np.arange(nvol)[:, None], X[None, :], norb_phys)  # [r, X]
    cols = so_bond_oracles.ed_mode(0, X, norb_phys)                                   # [Y]
    if rho.shape != (nvol * nd, nvol * nd):
        raise ValueError("c1a_rho_r: density {0} is not ({1}, {1})".format(rho.shape, nvol * nd))
    return np.ascontiguousarray(rho[rows[:, :, None], cols[None, None, :]], dtype=np.complex128)
