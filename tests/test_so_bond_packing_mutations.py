"""Index-map mutations of the S_z-packed layout are detected.

Every packed consumer reads the sector map through ``so_bond._consumer_map``
(the consumers ``so_bond._MAP_CONSUMERS``). A mutation replaces the map
that ONE consumer reads (``dataclasses.replace`` of one field) while every
other consumer reads the true map; the test asserts that the consumer was
called and that the field it read differs from the true one, then that the
packed result misses the dense production result restricted to the
test-side sector lists (``tests.so_bond_oracles.ph_sector_lists``) by at
least 1e-3, while the unmutated run of the same comparison passes.

A mutation applied to EVERY consumer at once is invisible when it is a
consistent relabeling of the slots inside each sector: the packed chain
then equals the dense chain restricted to the relabeled lists, which is
why each injection above is one-sided.

The refusal of an archive whose stored ``sz_sector_index`` is permuted is
tested with the archive members (``tests/test_so_bond_packed_archive.py``).
"""
import contextlib
import dataclasses
import unittest
from unittest import mock

import numpy as np

from hwave.solver import so_bond
from hwave.solver.so_bond import (PackedArrays, PackedBondPairKernel, PPGather, SectorCross,
                                  gamma_p, open_store, so_bond_bubble, so_bond_vertex)
from tests.so_bond_fixtures import (BATCH, BUBBLE_NMAT, PB_BETA, PB_GEOM_L, PB_ND,
                                    PB_TABLES_SZ, exact, pair_kernel, pb_gamma_label, pb_green,
                                    static_pair_inputs)
from tests.so_bond_fixtures import FIXTURE_NAMES as NAMES
from tests.so_bond_fixtures import bubble_green as _green
from tests.so_bond_fixtures import chain_gamma as _gamma
from tests.so_bond_fixtures import chain_lists as _lists
from tests.so_bond_fixtures import chain_plan as _plan
from tests.so_bond_fixtures import chain_sector_map as _sm
from tests.so_bond_fixtures import chain_setup as _setup
from tests.so_bond_fixtures import cvec as _cvec
from tests.so_bond_fixtures import dress as _dress
from tests.so_bond_fixtures import dress_inputs as _dress_inputs
from tests.so_bond_fixtures import dressed as _dressed
from tests.so_bond_fixtures import gap_shape as _gap_shape
from tests.so_bond_fixtures import kernel_setup as _kernel_setup
from tests.so_bond_fixtures import pairing_parts as _parts
from tests.so_bond_fixtures import pb_setup as _pb_setup
from tests.so_bond_fixtures import pb_sigma_prod_packed as _sigma_prod_packed
from tests.so_bond_fixtures import production_bubbles as _bubbles
from tests.so_bond_fixtures import random_g2 as _random_g2
from tests.so_bond_fixtures import random_kernels as _kernels
from tests.so_bond_fixtures import random_packed as _random_packed
from tests.so_bond_fixtures import vertex_batch as _batch
from tests.so_bond_oracles import (FIXTURE_LZSZ, FIXTURE_SZ_GRID, phi2_gradient, rel,
                                   restrict)

#: every mutated gate misses by at least this; every unmutated gate passes at TOL
MISS = 1e-3
TOL = 1e-13
#: the pre-asserted relative difference of the two blocks a swap exchanges
DIFFER = 1e-2
#: inject at every consumer
ALL = "every consumer"
_TRUE_MAP = so_bond._consumer_map


@contextlib.contextmanager
def _inject(consumer, mutate):
    """Patch ``so_bond._consumer_map``: consumer ``consumer`` (or every
    consumer for :data:`ALL`) reads ``mutate(true map)``, the others the
    true map. Yields the list of ``(name, map read)`` calls."""
    calls = []

    def hooked(name, sector_map):
        sm = _TRUE_MAP(name, sector_map)
        if consumer in (name, ALL):
            sm = mutate(sm)
        calls.append((name, sm))
        return sm
    with mock.patch.object(so_bond, "_consumer_map", hooked):
        yield calls


def _same(a, b):
    """Whether two map fields (arrays, tuples, named tuples) are equal."""
    if isinstance(a, np.ndarray) or isinstance(b, np.ndarray):
        return np.array_equal(np.asarray(a), np.asarray(b))
    if isinstance(a, tuple) and isinstance(b, tuple):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    return a == b


# -- mutations (each changes one field of the true map) ----------------------

def _size(sm, sigma):
    return int(sm.sizes[sm.sigmas.index(sigma)])


def _roll_local(sigma):
    """``local`` cyclically shifted inside sector ``sigma``."""
    def mutate(sm):
        slots = np.flatnonzero(sm.sector_of == sigma)
        local = sm.local.copy()
        local[slots] = (local[slots] + 1) % _size(sm, sigma)
        return dataclasses.replace(sm, local=local)
    return mutate, "local"


def _swap_sector_of(sm):
    """``sector_of`` with the sectors +1 and -1 exchanged."""
    s = sm.sector_of
    return dataclasses.replace(sm, sector_of=np.where(s == 1, -1, np.where(s == -1, 1, s))
                               .astype(np.int8))


def _gather_map(fn):
    def mutate(sm):
        return dataclasses.replace(sm, pp_gather=tuple(fn(sm, p, g) for p, g in
                                                       zip(sm.pp_order, sm.pp_gather)))
    return mutate, "pp_gather"


def _roll_row_local(sm, p, g):
    """The rows read from sector 0 cyclically shifted inside sector 0."""
    n = _size(sm, 0)
    in0 = (g.sector == 0)[None]
    return g._replace(row_local=np.where(in0, (g.row_local + 1) % n, g.row_local))


def _swap_pp_sector(sm, p, g):
    s = g.sector
    return g._replace(sector=np.where(s == 1, -1, np.where(s == -1, 1, s)).astype(np.int8))


#: a sector value no ph sector has: the gather never reads such an entry
_NO_SECTOR = 2


def _drop_pp0_spin_flip(sm, p, g):
    """The pp 0 entries of the ph sectors +1 and -1 are not read."""
    if p != 0:
        return g
    return g._replace(sector=np.where(np.abs(g.sector) == 1, _NO_SECTOR, g.sector)
                      .astype(np.int8))


def _cross_map(field, fn):
    def mutate(sm):
        return dataclasses.replace(sm, **{field: tuple(
            fn(sigma, cross) for sigma, cross in zip(sm.sigmas, getattr(sm, field)))})
    return mutate, field


def _same_sector(sigma, cross):
    """Read target sector ``sigma`` from sector ``sigma`` (true permutation)."""
    return SectorCross(sigma, cross.perm)


def _identity_reading(sigma, cross):
    """Read target sector ``sigma`` from itself, unpermuted."""
    return SectorCross(sigma, np.arange(len(cross.perm)))


def _roll_perm_0(sigma, cross):
    """The permutation onto sector 0 cyclically shifted."""
    return cross._replace(perm=np.roll(cross.perm, 1)) if sigma == 0 else cross


def _drop_minus_one(sm):
    return dataclasses.replace(sm, sector_order=("0", "+1"))


def _table(rows):
    return "\n" + "\n".join("  %-58s %.3e %s %.0e" % r for r in rows)


class _Gate(unittest.TestCase):

    def _verdict(self, rows):
        """``rows``: ``(label, value, op, bound)`` with op ``"<="`` (an
        unmutated gate) or ``">="`` (a mutated gate or a pre-assertion)."""
        bad = [r for r in rows if not (r[1] <= r[3] if r[2] == "<=" else r[1] >= r[3])]
        self.assertEqual(bad, [], _table(rows))

    def _applied(self, calls, consumer, field, true_map):
        """The mutated consumer was called and read a field that differs."""
        read = [sm for name, sm in calls if name == consumer]
        self.assertTrue(read, "consumer {!r} was not called ({})".format(
            consumer, sorted({n for n, _ in calls})))
        for sm in read:
            self.assertFalse(_same(getattr(sm, field), getattr(true_map, field)),
                             (consumer, field))
        for name, sm in calls:
            if name != consumer:
                self.assertIs(sm, true_map, name)


# -- the gates ---------------------------------------------------------------

def _gamma_miss(fx):
    """Gamma_B packed (production filing) vs the dense Gamma_B restricted."""
    c, topo, _ = _setup(fx)
    p = so_bond_vertex(c, topo, fx.spatial_shape, plan=_plan(fx, "sz_packed"))
    return max(rel(b, w) for b, w in zip(p.blocks, restrict(_gamma(fx, "dense"), _lists(fx))))


def _store_chibar(fx):
    """The every-slice packed bubble of the canonical G (copied blocks)."""
    topo = _setup(fx)[1]
    plan = _plan(fx, "sz_packed", ("chi_dynamic",), BUBBLE_NMAT[fx])
    with open_store(plan, _sm(fx)) as st:
        so_bond_bubble(_green(fx), None, fx.beta, topo, fx.spatial_shape, plan=plan,
                       freq="all", sink=st)
        return [b.copy() for b in st.packed("chibar").blocks]


def _kernel_miss(V, G2, fx, nfreq, dense, setup):
    geom, sm, plan = setup
    k = PackedBondPairKernel(V, G2, geom, sm, plan=plan, spatial_shape=fx.spatial_shape,
                             beta=fx.beta, nfreq=nfreq)
    phi = _cvec(np.random.default_rng(31), int(np.prod(_gap_shape(fx, nfreq))))
    return rel(k.matvec(phi), dense.matvec(phi))


def _static_setup(fx):
    V, G2, geom, sm, plan = static_pair_inputs(fx)
    return V, G2, (geom, sm, plan)


class TestMutations(_Gate):
    """Gate P-mut: one-sided injections at one consumer each. Every mutated
    gate misses by >= 1e-3 while the unmutated gate passes (<= 1e-13); a
    swap of two blocks is pre-asserted non-vacuous (the exchanged blocks
    differ by >= 1e-2 relative)."""

    def test_gamma_packer_within_sector(self):
        """The Gamma_B packer reads local positions permuted inside sector 0."""
        rows = []
        for fx in (FIXTURE_SZ_GRID, FIXTURE_LZSZ):
            sm = _sm(fx)
            rows.append(("%s unmutated" % NAMES[fx], _gamma_miss(fx), "<=", TOL))
            mutate, field = _roll_local(0)
            with _inject("gamma_packer", mutate) as calls:
                rows.append(("%s local permuted in sector 0" % NAMES[fx], _gamma_miss(fx),
                             ">=", MISS))
            self._applied(calls, "gamma_packer", field, sm)
        self._verdict(rows)

    def test_kernel_within_sector(self):
        """The kernel gather reads sector-0 rows permuted inside sector 0."""
        rows = []
        for fx, nfreq in ((FIXTURE_SZ_GRID, 1), (FIXTURE_LZSZ, 4)):
            _, dense, V, G2 = _kernels(fx, nfreq)
            setup = _kernel_setup(fx)
            rows.append(("%s nfreq %d unmutated" % (NAMES[fx], nfreq),
                         _kernel_miss(V, G2, fx, nfreq, dense, setup), "<=", TOL))
            mutate, field = _gather_map(_roll_row_local)
            with _inject("kernel_gather", mutate) as calls:
                rows.append(("%s nfreq %d row_local permuted in sector 0" % (NAMES[fx], nfreq),
                             _kernel_miss(V, G2, fx, nfreq, dense, setup), ">=", MISS))
            self._applied(calls, "kernel_gather", field, setup[1])
        self._verdict(rows)

    def test_kernel_sector_swap(self):
        """The kernel gather reads +1 entries from block -1 and vice
        versa, on the physical static vertex of the L_z S_z fixture. The
        swap substitutes V@-1 for V@+1 (both blocks in their local order,
        local i <-> (m; a, b) in both); Gamma@+1 = Gamma@-1, so the
        difference is the fluctuation part driven by chi."""
        fx = FIXTURE_LZSZ
        V, G2, setup = _static_setup(fx)
        dense = pair_kernel(V.to_dense(), G2, setup[0], fx.nd, fx.spatial_shape, fx.beta)
        rows = [("lzsz rel |V@+1 - V@-1| (pre-assertion)", rel(V[+1], V[-1]), ">=", DIFFER),
                ("lzsz static unmutated", _kernel_miss(V, G2, fx, 1, dense, setup), "<=", TOL)]
        mutate, field = _gather_map(_swap_pp_sector)
        with _inject("kernel_gather", mutate) as calls:
            rows.append(("lzsz static +1 <-> -1 in pp_gather.sector",
                         _kernel_miss(V, G2, fx, 1, dense, setup), ">=", MISS))
        self._applied(calls, "kernel_gather", field, setup[1])
        self._verdict(rows)

    def test_kernel_drop_pp0_spin_flip(self):
        """The pp 0 gather skips the entries of the ph sectors +1 and -1.
        The skipped entries of the gather workspace are never written; the
        workspace is zero-initialized during the construction so that they
        read as zeros (a deterministic drop)."""
        rows = []
        for fx, nfreq in ((FIXTURE_SZ_GRID, 1), (FIXTURE_LZSZ, 4)):
            _, dense, V, G2 = _kernels(fx, nfreq)
            setup = _kernel_setup(fx)
            g0 = setup[1].pp_gather[setup[1].pp_order.index(0)]
            self.assertTrue(np.any(np.abs(g0.sector) == 1))           # entries to drop exist
            rows.append(("%s nfreq %d unmutated" % (NAMES[fx], nfreq),
                         _kernel_miss(V, G2, fx, nfreq, dense, setup), "<=", TOL))
            mutate, field = _gather_map(_drop_pp0_spin_flip)
            with _inject("kernel_gather", mutate) as calls, \
                    mock.patch.object(np, "empty", np.zeros):
                miss = _kernel_miss(V, G2, fx, nfreq, dense, setup)
            rows.append(("%s nfreq %d pp 0 entries of sectors +-1 dropped" % (NAMES[fx], nfreq),
                         miss, ">=", MISS))
            self._applied(calls, "kernel_gather", field, setup[1])
        self._verdict(rows)

    def test_w_same_sector(self):
        """``p_b`` reads block sigma of W from block sigma of M (the true
        permutation). Pre-asserted: P(M@+1) and P(M@-1) differ (P the true
        permutation onto +1), which is what the mutation substitutes."""
        fx = FIXTURE_LZSZ
        xd, xp, gd, gp, Yd, Yp = _dressed(fx, "batch")
        Md = gd[None] @ (Yd - 0.5 * xd) @ gd[None]
        Mp = PackedArrays(xp.sector_map, tuple(
            g[None] @ (y - 0.5 * x) @ g[None] for g, y, x in zip(gp.blocks, Yp.blocks,
                                                               xp.blocks)))
        sm = Mp.sector_map
        perm = sm.pb_cross[sm.sigmas.index(1)].perm

        def P(b):
            return np.take(np.take(b, perm, axis=-2), perm, axis=-1)

        want = restrict(so_bond.p_b(Md, fx.nd), _lists(fx))

        def miss():
            W = so_bond.p_b(Mp, fx.nd)
            return max(rel(b, w) for b, w in zip(W.blocks, want))

        rows = [("lzsz rel |P(M@+1) - P(M@-1)| (pre-assertion)", rel(P(Mp[+1]), P(Mp[-1])),
                 ">=", DIFFER),
                ("lzsz unmutated", miss(), "<=", TOL)]
        mutate, field = _cross_map("pb_cross", _same_sector)
        with _inject("w_cross", mutate) as calls:
            rows.append(("lzsz W@sigma read from M@sigma", miss(), ">=", MISS))
        self._applied(calls, "w_cross", field, sm)
        self._verdict(rows)

    def test_w_identity_reading_misses_pb_pin(self):
        """The packed FLEX chain with ``W@sigma = M@sigma`` (identity
        reading at ``p_b``) misses the second-order functional's gradient,
        which the true reading matches to 1e-12 (seed 3, nearest bond)."""
        g = pb_green(3, True)
        S_or = phi2_gradient(g, PB_BETA, pb_gamma_label(PB_TABLES_SZ[False], PB_GEOM_L[False]),
                             PB_GEOM_L[False])
        rows = [("Phi_2 gradient, packed chain, unmutated",
                 rel(_sigma_prod_packed(False, 3), S_or), "<=", 1e-12)]
        mutate, field = _cross_map("pb_cross", _identity_reading)
        with _inject("w_cross", mutate) as calls:
            rows.append(("Phi_2 gradient, packed chain, identity reading",
                         rel(_sigma_prod_packed(False, 3), S_or), ">=", MISS))
        self._applied(calls, "w_cross", field, _pb_map())
        self._verdict(rows)

    def test_store_scatter_within_sector(self):
        """The packed store scatters the pair blocks with local positions
        permuted inside sector +1 (every-slice bubble, L_z S_z fixture)."""
        fx = FIXTURE_LZSZ
        sm = _sm(fx)
        want = restrict(_bubbles(fx, "all")[0][0], _lists(fx))

        def miss():
            return max(rel(b, w) for b, w in zip(_store_chibar(fx), want))

        rows = [("lzsz unmutated", miss(), "<=", TOL)]
        mutate, field = _roll_local(1)
        with _inject("store_scatter", mutate) as calls:
            rows.append(("lzsz local permuted in sector +1", miss(), ">=", MISS))
        self._applied(calls, "store_scatter", field, sm)
        self._verdict(rows)

    def test_store_scatter_sector_swap(self):
        """The packed store scatters the +1 pairs into block -1 and vice
        versa (every-slice bubble, L_z S_z fixture, whose +1 and -1 bubble
        blocks differ)."""
        fx = FIXTURE_LZSZ
        sm = _sm(fx)
        want = restrict(_bubbles(fx, "all")[0][0], _lists(fx))

        def miss():
            return max(rel(b, w) for b, w in zip(_store_chibar(fx), want))

        rows = [("lzsz rel |chibar@+1 - chibar@-1| (pre-assertion)", rel(want[1], want[2]),
                 ">=", DIFFER),
                ("lzsz unmutated", miss(), "<=", TOL)]
        with _inject("store_scatter", _swap_sector_of) as calls:
            rows.append(("lzsz sector_of +1 <-> -1", miss(), ">=", MISS))
        self._applied(calls, "store_scatter", "sector_of", sm)
        self._verdict(rows)

    def _gamma_p_rows(self, mutate, field, label):
        rows = []
        for fx in (FIXTURE_SZ_GRID, FIXTURE_LZSZ):
            geom, sm = _setup(fx)[2], _sm(fx)
            x = _random_packed(sm, fx.nvol, 4)
            dense = x.to_dense()
            want = restrict(gamma_p(dense, geom, out=np.empty_like(dense)), _lists(fx))

            def miss():
                out = gamma_p(x, geom, out=PackedArrays.zeros(sm, (fx.nvol,)))
                return max(rel(b, w) for b, w in zip(out.blocks, want))

            rows += [("%s rel |x@+1 - x@-1| (pre-assertion)" % NAMES[fx], rel(x[+1], x[-1]),
                      ">=", DIFFER),
                     ("%s unmutated" % NAMES[fx], miss(), "<=", TOL)]
            with _inject("gamma_p", mutate) as calls:
                rows.append(("%s %s" % (NAMES[fx], label), miss(), ">=", MISS))
            self._applied(calls, "gamma_p", field, sm)
        return rows

    def test_gamma_p_cross(self):
        """The exchange image reads target sector sigma from sector sigma
        (true permutation), on random block-diagonal input whose +1 and -1
        blocks differ (on Gamma_B itself the reading is invisible:
        Gamma@+1 = Gamma@-1)."""
        mutate, field = _cross_map("gp_cross", _same_sector)
        self._verdict(self._gamma_p_rows(mutate, field, "gp_cross source = sigma"))

    def test_gamma_p_within_sector(self):
        """The exchange image onto sector 0 reads a cyclically shifted
        permutation (random block-diagonal input)."""
        mutate, field = _cross_map("gp_cross", _roll_perm_0)
        self._verdict(self._gamma_p_rows(mutate, field, "gp_cross perm shifted in sector 0"))

    def test_bare_term_exchange_reading(self):
        """The bare term ``1/4 (Gamma + Gamma^P)`` of the pairing parts and
        of the vertex batch with the exchange image read without the slot
        permutation (each sector read from itself, unpermuted, at
        ``gamma_p``), on a random vertex; the mutated image is pre-asserted
        to differ from the true one, and the exchange-symmetry refusal is
        bypassed (bound set to infinity) in both layouts."""
        fx = FIXTURE_LZSZ
        geom, sm, lists = _setup(fx)[2], _sm(fx), _lists(fx)
        gp_ = _random_packed(sm, fx.nvol, 4)
        gd = gp_.to_dense()
        _, _, _, _, Yd, Yp = _dressed(fx, "batch")
        cp = PackedArrays(sm, tuple(b[0].copy() for b in Yp.blocks))
        cd = Yd[0].copy()
        image = restrict(gamma_p(gd, geom, out=np.empty_like(gd)), lists)
        mutate, field = _cross_map("gp_cross", _identity_reading)
        with _inject("gamma_p", mutate) as calls:
            mutated = gamma_p(gp_, geom, out=PackedArrays.zeros(sm, (fx.nvol,)))
        self._applied(calls, "gamma_p", field, sm)
        with mock.patch.object(so_bond, "_EXCHANGE_RTOL", np.inf):
            want_parts = restrict(_parts(fx, "dense", gd, cd)[0], lists)
            want_batch = restrict(_batch(fx, "dense", gd, Yd), lists)

            def miss_parts():
                vb = _parts(fx, "sz_packed", gp_, cp)[0]
                return max(rel(b, w) for b, w in zip(vb.blocks, want_parts))

            def miss_batch():
                V = _batch(fx, "sz_packed", gp_, Yp)
                return max(rel(b, w) for b, w in zip(V.blocks, want_batch))

            label = "exchange image read without the slot permutation"
            rows = [("lzsz rel |mutated image - Gamma^P| (pre-assertion)",
                     max(rel(b, w) for b, w in zip(mutated.blocks, image)), ">=", DIFFER),
                    ("lzsz parts V_bare unmutated", miss_parts(), "<=", TOL),
                    ("lzsz batch V unmutated", miss_batch(), "<=", TOL)]
            with _inject("gamma_p", mutate) as calls:
                rows.append(("lzsz parts V_bare, " + label, miss_parts(), ">=", MISS))
            self._applied(calls, "gamma_p", field, sm)
            with _inject("gamma_p", mutate) as calls:
                rows.append(("lzsz batch V, " + label, miss_batch(), ">=", MISS))
            self._applied(calls, "gamma_p", field, sm)
        self._verdict(rows)

    def test_dress_drops_a_sector(self):
        """The dressing reads a sector order without "-1": that sector is
        not dressed and stays zero in the result (built on the map of X),
        and a packed array on the mutated map itself is refused."""
        fx = FIXTURE_LZSZ
        xd, xp, gd, gp, l0 = _dress_inputs(fx, "static")
        want = restrict(_dress(fx, "dense", xd, gd, l0, "static")[0], _lists(fx))

        def miss():
            chi = _dress(fx, "sz_packed", xp, gp, l0, "static")[0]
            return max(rel(b, w) for b, w in zip(chi.blocks, want))

        rows = [("lzsz unmutated", miss(), "<=", TOL)]
        with _inject("dress_sectors", _drop_minus_one) as calls:
            rows.append(("lzsz sector_order without '-1'", miss(), ">=", MISS))
        self._applied(calls, "dress_sectors", "sector_order", xp.sector_map)
        self._verdict(rows)
        with self.assertRaisesRegex(ValueError, exact(
                "PackedArrays: sector_map.sector_order must be ('0', '+1', '-1'), got "
                "('0', '+1')")):
            PackedArrays.zeros(_drop_minus_one(xp.sector_map), (1,))


def _pb_map():
    """The sector map of the packed FLEX chain of ``_sigma_prod_packed``."""
    return so_bond.sz_sector_map(_pb_setup(False, True)[1], PB_ND)


# -- the shared mutation -------------------------------------------------------

def _relabeled(sm, reverse, order):
    """The sector map after the bijection that relabels the slots inside
    each sector: packed position ``offsets[s] + i`` holds the dense slot
    ``index_s[order[s][i]]`` (``index_s`` the ascending slots of sector
    ``s``). Every field is recomputed from that slot list by its
    definition: ``local`` inverts it, ``pb_cross`` / ``gp_cross`` are the
    slot maps ``(m; x, y) -> (m; y, x)`` / ``(reverse[m]; y, x)`` read
    through the new ``local``, and ``pp_gather`` reads ``V_{ab; ABCD}``
    at the new local positions of ``(a; A, B)`` and ``(b; C, D)``."""
    nd, B = sm.nd, sm.B
    npair = nd * nd
    lists = [sm.index[sm.offsets[s]:sm.offsets[s + 1]][np.asarray(order[s])] for s in range(3)]
    local = np.empty(sm.ND, dtype=np.int64)
    for l in lists:
        local[l] = np.arange(len(l))

    def cross(f):
        out = []
        for s, sigma in enumerate(sm.sigmas):
            tgt = f(lists[s])
            assert set(sm.sector_of[tgt].tolist()) == {-sigma}
            out.append(SectorCross(-sigma, local[tgt]))
        return tuple(out)

    def split(I):
        return I // npair, (I % npair) // nd, I % nd

    def swap(I):
        m, x, y = split(I)
        return m * npair + y * nd + x

    def image(I):
        m, x, y = split(I)
        return np.asarray(reverse)[m] * npair + y * nd + x

    gathers = []
    chan = np.arange(B)[:, None, None] * npair
    for pairs in sm.pp_pairs:
        A, D = np.divmod(pairs, nd)                  # output pair (A, D)
        Bm, C = np.divmod(pairs, nd)                 # input pair (B, C)
        row = chan + (A[:, None] * nd + Bm[None, :])[None]
        col = chan + (C[None, :] * nd + D[:, None])[None]
        gathers.append(PPGather(sm.sector_of[row[0]], local[row], local[col]))
    return dataclasses.replace(sm, index=np.concatenate(lists), local=local,
                               pb_cross=cross(swap), gp_cross=cross(image),
                               pp_gather=tuple(gathers)), lists


class TestSharedMutation(_Gate):
    """Gate P-mut (shared mutation): the same consistent relabeling read by
    every consumer is invisible -- each packed result of the chain equals
    the dense production result restricted to the RELABELED sector lists
    (<= 1e-13), and the kernel equals the dense kernel; this is why each
    mutation of the one-sided gates is injected at one consumer only."""

    def test_identity_relabeling_is_the_map(self):
        """Construction check: the test-side recomputation with the identity
        order reproduces every field of the production map."""
        fx = FIXTURE_LZSZ
        sm = _sm(fx)
        same, _ = _relabeled(sm, _setup(fx)[2].reverse, [np.arange(n) for n in sm.sizes])
        for f in dataclasses.fields(sm):
            self.assertTrue(_same(getattr(same, f.name), getattr(sm, f.name)), f.name)

    def test_shared_map_mutation_is_invisible(self):
        fx = FIXTURE_LZSZ
        c, topo, geom = _setup(fx)
        sm = _sm(fx)
        rng = np.random.default_rng(14)
        order = [rng.permutation(n) for n in sm.sizes]
        relabeled, lists = _relabeled(sm, geom.reverse, order)
        for f in ("index", "local", "pb_cross", "gp_cross", "pp_gather"):
            self.assertFalse(_same(getattr(relabeled, f), getattr(sm, f)), f)
        nfreq = BATCH.stop - BATCH.start
        G2 = _random_g2(fx, nfreq, 6)

        # the dense chain
        gd = _gamma(fx, "dense")
        xd_all = _bubbles(fx, "all")[0][0]
        xd = xd_all[BATCH]
        Yd = _dress(fx, "dense", xd, gd, BATCH.start, "batch")[0]
        Wd = so_bond.p_b(gd[None] @ (Yd - 0.5 * xd) @ gd[None], fx.nd)
        Vd = _batch(fx, "dense", gd, Yd)
        kd = pair_kernel(Vd, G2, geom, fx.nd, fx.spatial_shape, fx.beta)

        def miss(packed, dense):
            return max(rel(b, restrict(dense, [l])[0]) for b, l in zip(packed.blocks, lists))

        rows = []
        with _inject(ALL, lambda m: relabeled) as calls:
            gp_ = so_bond_vertex(c, topo, fx.spatial_shape, plan=_plan(fx, "sz_packed"))
            rows.append(("Gamma_B", miss(gp_, gd), "<=", TOL))
            img = gamma_p(gp_, geom, out=PackedArrays.zeros(sm, (fx.nvol,)))
            rows.append(("Gamma^P", miss(img, gamma_p(gd, geom, out=np.empty_like(gd))),
                         "<=", TOL))
            xp = PackedArrays(sm, tuple(b[BATCH].copy() for b in _store_chibar(fx)))
            rows.append(("chibar (store, every slice)", miss(xp, xd), "<=", TOL))
            Yp = _dress(fx, "sz_packed", xp, gp_, BATCH.start, "batch")[0]
            rows.append(("chi (dressed batch)", miss(Yp, Yd), "<=", TOL))
            Mp = PackedArrays(sm, tuple(g[None] @ (y - 0.5 * x) @ g[None] for g, y, x in
                                        zip(gp_.blocks, Yp.blocks, xp.blocks)))
            rows.append(("W = P_B(M)", miss(so_bond.p_b(Mp, fx.nd), Wd), "<=", TOL))
            Vp = _batch(fx, "sz_packed", gp_, Yp)
            rows.append(("V (vertex batch)", miss(Vp, Vd), "<=", TOL))
            rows.append(("kernel matvec", _kernel_miss(Vp, G2, fx, nfreq, kd, _kernel_setup(fx)),
                         "<=", TOL))
        self.assertEqual({n for n, _ in calls}, set(so_bond._MAP_CONSUMERS))
        self.assertTrue(all(m is relabeled for _, m in calls))
        self._verdict(rows)


if __name__ == "__main__":
    unittest.main()
