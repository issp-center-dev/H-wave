"""The bond map ``P_B`` of the spin-orbital FLEX effective interaction.

``so_bond.p_b(M, nd)`` builds ``W = P_B[M]`` from the screened ring
``M = 1/2 Gamma_B chibar_B Gamma_B`` of the code slots:

``P_B[M](q)[(m; x, y), (m'; z, w)] = M(q)[(m; y, x), (m'; w, z)]``,

the orbital swap of both pair indices within each channel (no channel
reversal, no phase). In the S_z-packed layout the swap sends the ph sector
``sigma`` onto ``-sigma``: block ``sigma`` of ``W`` is read from block
``-sigma`` of ``M``.

``W`` is pinned by the second-order functional of the label vertex,
``Phi_2[G] = -(T/8) sum_{l, q} Tr[(chibar(q, i nu_l) Gamma(-q))**2]``: the
production self-energy ``flex_bond.calc_self_energy_bond`` of ``W`` equals
``beta dPhi_2/dG`` (``tests.so_bond_oracles.phi2_gradient``, an oracle
that never imports the production compiler or the bond module).
"""
import functools
import types
import unittest

import numpy as np

from hwave.solver import eliashberg_bond as eb
from hwave.solver import flex_bond
from hwave.solver.so_bond import PackedArrays, so_bond_bubble, so_bond_vertex, sz_sector_map
from hwave.solver import so_bond
from tests.so_bond_fixtures import PB_BETA as BETA
from tests.so_bond_fixtures import PB_GEOM_L as GEOM_L
from tests.so_bond_fixtures import PB_L as L
from tests.so_bond_fixtures import PB_ND as ND
from tests.so_bond_fixtures import PB_NMAT as NMAT
from tests.so_bond_fixtures import PB_SHAPE as SHAPE
from tests.so_bond_fixtures import PB_TABLES_ALL as TABLES_ALL
from tests.so_bond_fixtures import PB_TABLES_SZ as TABLES_SZ
from tests.so_bond_fixtures import exact
from tests.so_bond_fixtures import pb_gamma_label as gamma_label
from tests.so_bond_fixtures import pb_green as _green
from tests.so_bond_fixtures import pb_plan as _plan
from tests.so_bond_fixtures import pb_setup as _setup
from tests.so_bond_fixtures import pb_sigma as _sigma
from tests.so_bond_fixtures import pb_sigma_prod_packed as _sigma_prod_packed
from tests.so_bond_fixtures import random_packed as _random_packed
from tests.so_bond_oracles import ph_sector_lists, phi2_gradient, rel, restrict

SEEDS = ((3, False), (5, False), (3, True))           # (Green seed, next-nearest shell)


def _relmax(a, b):
    return rel(a, b)


def _off_sector_sigma(S):
    """The spin-off-diagonal entries of a self-energy (norb_phys = 1)."""
    return np.stack([S[..., 0, 1], S[..., 1, 0]])


@functools.lru_cache(maxsize=None)
def _ring_dense(nnn, sz, seed):
    """``M = 1/2 Gamma_B chibar_B Gamma_B`` (nmat, nvol, ND, ND), production."""
    c, topo, _ = _setup(nnn, sz)
    plan = _plan(nnn, sz, "dense")
    gam = so_bond_vertex(c, topo, SHAPE, plan=plan)
    with flex_bond.BondBlockStore(NMAT, L, plan.ND, ND * ND, ("chibar",)) as st:
        so_bond_bubble(_green(seed, sz)[None].copy(), None, BETA, topo, SHAPE, plan=plan,
                       freq="all", sink=st)
        X = np.array(st.get_freq_batch("chibar", 0, NMAT))
    return 0.5 * gam[None] @ X @ gam[None]


def _reversed_reading(M, nd, topo):
    """The wrong reading ``(reverse[m]; y, x)`` of both pair indices."""
    npair = nd * nd
    I = np.arange(M.shape[-1])
    m, x, y = I // npair, (I % npair) // nd, I % nd
    f = np.asarray(topo.reverse)[m] * npair + y * nd + x
    return M[..., f[:, None], f]


def _sigma_prod(nnn, seed, reading):
    """The production self-energy of ``W = reading(M)`` (dense layout, all
    seven types)."""
    _, topo, geom = _setup(nnn, False)
    M = _ring_dense(nnn, False, seed)
    if reading == "p_b":
        W = so_bond.p_b(M, ND)
    elif reading == "reversed":
        W = _reversed_reading(M, ND, topo)
    else:                                                  # no swap
        W = M
    return _sigma(eb.ArrayBlockSource({"W": W}, ND * ND), _green(seed, False), geom)


def _table(rows):
    return "\n".join("  seed %d nnn %-5s  %-9s %.3e" % r for r in rows)


class TestPBPin(unittest.TestCase):
    """Gate P-PB: the bond map P_B is the per-channel orbital swap.

    The production self-energy of ``W = P_B[1/2 Gamma_B chibar_B Gamma_B]``
    (``flex_bond.calc_self_energy_bond``) equals ``beta dPhi_2/dG`` of the
    second-order functional ``Phi_2 = -(T/8) sum Tr[(chibar Gamma)**2]``
    to 1e-12 (max-norm over the whole array); the readings ``(reverse[m];
    y, x)`` and "no swap" miss it by at least 1e-1. The gate pins ``P_B``
    given the production ``Gamma_B``: the oracle's vertex is filed from the
    interaction tables independently, so an error of ``Gamma_B`` would fail
    it too, but ``Gamma_B`` itself is pinned by the exact-diagonalization
    gates of the vertex.

    Dense: all seven interaction types (Ising included) and a random Green
    function with every spin entry. Packed: the S_z-conserving subset
    (PairLift dropped) and a spin-diagonal Green function, where the
    oracle's spin-off-diagonal self-energy vanishes as well.
    """

    def test_dense_every_type(self):
        rows, fails = [], []
        for seed, nnn in SEEDS:
            S_or = phi2_gradient(_green(seed, False), BETA, gamma_label(TABLES_ALL[nnn],
                                                                        GEOM_L[nnn]),
                                 GEOM_L[nnn])
            for reading in ("p_b", "reversed", "no swap"):
                r = _relmax(_sigma_prod(nnn, seed, reading), S_or)
                rows.append((seed, nnn, reading, r))
                if (r > 1e-12) if reading == "p_b" else not r >= 1e-1:
                    fails.append(rows[-1])
        self.assertEqual(fails, [], "\n" + _table(rows))

    def test_packed_conserving_variant(self):
        rows, fails = [], []
        for seed, nnn in SEEDS:
            g = _green(seed, True)
            S_or = phi2_gradient(g, BETA, gamma_label(TABLES_SZ[nnn], GEOM_L[nnn]), GEOM_L[nnn])
            S_pk = _sigma_prod_packed(nnn, seed)
            r = _relmax(S_pk, S_or)
            rows.append((seed, nnn, "packed", r))
            if not (r <= 1e-12
                    and np.abs(_off_sector_sigma(S_pk)).max() == 0.0
                    and np.abs(_off_sector_sigma(S_or)).max() <= 1e-13 * np.abs(S_or).max()):
                fails.append(rows[-1] + (np.abs(_off_sector_sigma(S_pk)).max(),
                                         np.abs(_off_sector_sigma(S_or)).max()))
        self.assertEqual(fails, [], "\n" + _table(rows))


def _swap_formula(M, nd):
    """``P_B`` by explicit loops over the code slots ``(m; x, y)``."""
    ND = M.shape[-1]
    npair = nd * nd
    src = np.empty(ND, dtype=np.int64)
    for m in range(ND // npair):
        for x in range(nd):
            for y in range(nd):
                src[m * npair + x * nd + y] = m * npair + y * nd + x
    out = np.empty_like(M)
    for I in range(ND):
        for J in range(ND):
            out[..., I, J] = M[..., src[I], src[J]]
    return out


def _random(shape, seed):
    rng = np.random.default_rng(seed)
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


REV3 = types.SimpleNamespace(reverse=np.array([0, 2, 1]))
REV2 = types.SimpleNamespace(reverse=np.array([1, 0]))


class TestPBFormula(unittest.TestCase):
    """Construction check: ``p_b`` is the index formula in both layouts, the
    packed result is the restriction of the dense one, and every refusal
    text is exact (raised before anything is written)."""

    def test_dense_formula(self):
        for B, nd, lead in ((3, 2, (2, 3)), (2, 4, (3,)), (1, 2, ())):
            M = _random(lead + (B * nd * nd,) * 2, 11 + B)
            want = _swap_formula(M, nd)
            self.assertTrue(np.array_equal(so_bond.p_b(M, nd), want))
            out = np.full_like(M, np.nan)
            self.assertIs(so_bond.p_b(M, nd, out=out), out)
            self.assertTrue(np.array_equal(out, want))
            # a non-contiguous out and a non-contiguous M
            buf = np.full(lead + (B * nd * nd, 2 * B * nd * nd), np.nan, dtype=np.complex128)
            out2 = buf[..., ::2]
            self.assertIs(so_bond.p_b(np.asfortranarray(M), nd, out=out2), out2)
            self.assertTrue(np.array_equal(out2, want))
            # an involution
            self.assertTrue(np.array_equal(so_bond.p_b(want, nd), M))

    def test_packed_is_restricted_dense(self):
        for topo, nd, lead in ((REV3, 2, (2, 3)), (REV2, 4, (3,))):
            sm = sz_sector_map(topo, nd)
            Mp = _random_packed(sm, lead, 29)
            want = restrict(so_bond.p_b(Mp.to_dense(), nd), ph_sector_lists(sm.B, nd))
            got = so_bond.p_b(Mp, nd)
            self.assertIsInstance(got, PackedArrays)
            self.assertIs(got.sector_map, sm)
            for b, w in zip(got.blocks, want):
                self.assertTrue(np.array_equal(b, w))
            out = PackedArrays(sm, tuple(np.full_like(b, np.nan) for b in Mp.blocks))
            self.assertIs(so_bond.p_b(Mp, nd, out=out), out)
            for b, w in zip(out.blocks, want):
                self.assertTrue(np.array_equal(b, w))
            back = so_bond.p_b(got, nd)
            for b, w in zip(back.blocks, Mp.blocks):
                self.assertTrue(np.array_equal(b, w))

    def test_refusals(self):
        nd_text = "p_b: nd must be a positive even integer (nd = 2*norb_phys), got {!r}"
        M = _random((2, 12, 12), 1)
        for bad in (0, 3, -2, True, 2.0, "2"):
            with self.assertRaisesRegex(ValueError, exact(nd_text.format(bad))):
                so_bond.p_b(M, bad)
        dense_text = ("p_b: M must be a complex128 (..., ND, ND) array with ND a positive "
                      "multiple of nd**2 = 4, got {}")
        for bad, got in ((M.astype(np.complex64), "complex64 (2, 12, 12)"),
                         (M.real.copy(), "float64 (2, 12, 12)"),
                         (_random((2, 10, 10), 2), "complex128 (2, 10, 10)"),
                         (_random((2, 12, 8), 3), "complex128 (2, 12, 8)"),
                         (_random((12,), 4), "complex128 (12,)"),
                         (np.zeros((0, 0), dtype=np.complex128), "complex128 (0, 0)"),
                         ([[1j]], "list")):
            with self.assertRaisesRegex(ValueError, exact(dense_text.format(got))):
                so_bond.p_b(bad, 2)
        out_text = "p_b: out must be a complex128 array of shape (2, 12, 12), got {}"
        for bad, got in ((np.empty((2, 12, 12), dtype=np.complex64), "complex64 (2, 12, 12)"),
                         (np.empty((1, 12, 12), dtype=np.complex128), "complex128 (1, 12, 12)"),
                         ([0], "list")):
            with self.assertRaisesRegex(ValueError, exact(out_text.format(got))):
                so_bond.p_b(M, 2, out=bad)
        with self.assertRaisesRegex(ValueError, exact("p_b: out shares memory with M")):
            so_bond.p_b(M, 2, out=M)
        with self.assertRaisesRegex(ValueError, exact("p_b: out shares memory with M")):
            so_bond.p_b(M, 2, out=M[::-1])

        sm = sz_sector_map(REV3, 2)
        Mp = _random_packed(sm, (2,), 5)
        mixed = "p_b: M and out must both be PackedArrays or both dense arrays"
        with self.assertRaisesRegex(ValueError, exact(mixed)):
            so_bond.p_b(Mp, 2, out=np.empty((2, 12, 12), dtype=np.complex128))
        with self.assertRaisesRegex(ValueError, exact(mixed)):
            so_bond.p_b(M, 2, out=PackedArrays.zeros(sm, (2,)))
        with self.assertRaisesRegex(ValueError, exact(
                "p_b: the sector map of M has nd = 2, the call has nd = 4")):
            so_bond.p_b(Mp, 4)
        with self.assertRaisesRegex(ValueError, exact(nd_text.format(3))):
            so_bond.p_b(Mp, 3)
        other = sz_sector_map(types.SimpleNamespace(reverse=np.array([0, 1, 2])), 2)
        sentinel = PackedArrays(other, tuple(np.full((2, n, n), 7.0 + 0j) for n in other.sizes))
        with self.assertRaisesRegex(ValueError, exact(
                "p_b: the sector map of out is not the sector map of M")):
            so_bond.p_b(Mp, 2, out=sentinel)
        self.assertTrue(all(np.all(b == 7.0) for b in sentinel.blocks))
        with self.assertRaisesRegex(ValueError, exact(
                "p_b: the packed out has leading shape (3,), M has (2,)")):
            so_bond.p_b(Mp, 2, out=PackedArrays.zeros(sm, (3,)))
        with self.assertRaisesRegex(ValueError, exact("p_b: out shares memory with M")):
            so_bond.p_b(Mp, 2, out=Mp)
        alias = PackedArrays(sm, (Mp.blocks[0].copy(), Mp.blocks[1].copy(), Mp.blocks[1]))
        with self.assertRaisesRegex(ValueError, exact("p_b: out shares memory with M")):
            so_bond.p_b(Mp, 2, out=alias)


if __name__ == "__main__":
    unittest.main()
