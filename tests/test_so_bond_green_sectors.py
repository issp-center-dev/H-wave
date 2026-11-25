"""The canonical per-spin-sector band Green function and the borrowed-G
projection of the spin-orbital bond path.

* ``so_bond.canonical_band_eigen`` diagonalizes the spin-block band
  ``epsilon_k`` per spin block when the band conserves S_z (structural
  zeros between the blocks; per k a stable ascending sort, the eigenvalue
  order of ``sc._calc_eigenvalues``) and refuses a band read in the wrong
  index order; on a non-conserving band it is ``sc._calc_eigenvalues``.
* ``so_bond.green_leak`` / ``so_bond.project_green_sectors`` measure and
  remove the S_z-off-sector entries of a Green function.
* ``so_bond.project_borrowed_green`` uses a borrowed Green function as is
  (no leak), projects a round-off leak with one WARNING, and refuses a
  larger leak as an incompatible Hamiltonian.
"""
import unittest
from unittest import mock

import numpy as np

from hwave.solver.eliashberg import so_solve as sc_so
from hwave.solver.eliashberg import bond as sc_bond
from hwave.solver.eliashberg import band as sc_band
from hwave.solver import so_bond, spin_orbital
from hwave.solver.fingerprint import green_digest
from hwave.solver.so_records import compile_interaction
from hwave.solver.so_bond import (PackedArrays, SOBondCaps, band_conserves_sz,
                                  canonical_band_eigen, green_leak, plan_so_bond,
                                  project_borrowed_green, project_green_sectors, select_layout,
                                  so_bond_bubble, so_bond_topology, sz_sector_map)
from tests.so_bond_fixtures import assert_recomputed, canonical_green, exact, quiet
from tests.so_bond_oracles import (BASE, CELL, FIXTURE_LZSZ, FIXTURE_SZ_GRID, closed_row)

NMAT = 8
LAYOUT = "b_f_k_o_o_v1"
FIELDS = (0.0, 0.1)

WARN = ("the Green function '%s' has S_z-off-sector entries at the round-off level (leak = "
        "max off-sector |G| / max |G| = %.3e <= 1e-12: numerical leakage from a "
        "non-canonical producer); they are projected to exact zeros")
REFUSE = ("the Green function '%s' breaks S_z (leak = max off-sector |G| / max |G| = %.3e > "
          "1e-12) although this run's input conserves S_z; it was produced from a different "
          "Hamiltonian. Regenerate it from this run's input.")
INTERLEAVED = ("canonical_band_eigen: epsilon_k has S_z-off-sector entries (max %.3e) in the "
               "spin-block reading although the Transfer table conserves S_z (is it in the "
               "interleaved order?)")


def _extern(norb, h):
    """The Extern table of a Zeeman field ``h`` (``None`` at zero field)."""
    return {((0, 0, 0), (a, a)): 1.0 for a in range(norb)} if h else None


def _halves(norb):
    """The spin-up and spin-down halves of the spin-block index."""
    return slice(0, norb), slice(norb, 2 * norb)


def _epsilon_k(fixture, h):
    """The spin-block band of ``fixture``'s Transfer table with the Zeeman
    field ``h`` as an Extern term, built by ``sc._so_epsilon_k``."""
    norb = fixture.norb_phys
    ks = [np.linspace(0, 2.0 * np.pi, n, endpoint=False) for n in fixture.spatial_shape]
    extern = _extern(norb, h)
    with quiet():
        return sc_so._so_epsilon_k({"mode": {"param": {"coeff_extern": h}}}, *ks,
                                fixture.transfer_dict, extern, norb)


def _green(fixture, eigen):
    """The bond-path band Green function ``(1, nmat, nvol, nd, nd)`` of an
    eigen-decomposition (``sc._build_bond_green``, internal build)."""
    ev, evec = eigen
    return sc_bond._build_bond_green(ev, evec, fixture.mu, fixture.beta, NMAT, 0.0, None,
                                spin_orbital_mode=True).deflated_kw


def _off(G, norb):
    """The S_z-off-sector entries of ``G`` (spin-block, last two axes)."""
    up, dn = _halves(norb)
    return np.concatenate([G[..., up, dn].ravel(), G[..., dn, up].ravel()])


def _band_green(fixture, h, tables, override):
    """The band Green function of the run (Transfer, Extern field ``h``,
    interaction ``tables``, ``so_bond_layout = override``): ``conserving``
    comes from ``band_conserves_sz``; the layout is selected as the run
    would, but the Green function never reads it. Returns ``(G, layout)``."""
    norb = fixture.norb_phys
    extern = _extern(norb, h)
    conserving = band_conserves_sz(fixture.transfer_dict, extern, norb).conserving
    layout, _ = select_layout(compile_interaction(tables, norb), fixture.transfer_dict,
                              override=override)
    eps = _epsilon_k(fixture, h)
    return _green(fixture, canonical_band_eigen(eps, norb, conserving=conserving)), layout


class TestCanonicalBandEigen(unittest.TestCase):
    """Gate P-G (i): the canonical band Green function of a conserving band,
    for h = 0 and h = 0.1, norb_phys 1 (5x3 grid) and 2 (chain)."""

    def test_eigenvalues_match_full_diagonalization(self):
        for fixture in (FIXTURE_SZ_GRID, FIXTURE_LZSZ):
            for h in FIELDS:
                with self.subTest(norb=fixture.norb_phys, h=h):
                    eps = _epsilon_k(fixture, h)
                    ev, evec = canonical_band_eigen(eps, fixture.norb_phys, conserving=True)
                    ref, _ = sc_band._calc_eigenvalues(eps)
                    self.assertEqual(ev.shape, ref.shape)
                    self.assertEqual(evec.shape, ref.shape + (ref.shape[-1],))
                    self.assertLessEqual(float(np.max(np.abs(ev - ref))), 1e-13)
                    self.assertTrue(np.all(np.diff(ev, axis=-1) >= 0))
                    # the eigenvectors diagonalize the band in this order
                    H = eps.transpose(2, 3, 4, 0, 1)
                    resid = H @ evec - evec * ev[..., None, :]
                    self.assertLessEqual(float(np.max(np.abs(resid))), 1e-13)
                    # structural zeros: no eigenvector mixes the spin blocks
                    norb = fixture.norb_phys
                    up, dn = _halves(norb)
                    pure_up = np.all(evec[..., dn, :] == 0, axis=-2)
                    pure_dn = np.all(evec[..., up, :] == 0, axis=-2)
                    self.assertTrue(np.all(pure_up ^ pure_dn))

    def test_stable_sort_of_ties(self):
        """Equal eigenvalues keep the block order (spin up first): a stable
        ascending sort with the matching column permutation."""
        eps = np.zeros((4, 4, 2, 1, 1), dtype=np.complex128)
        eps[[0, 1, 2, 3], [0, 1, 2, 3], 0, 0, 0] = [1.0, 0.0, 1.0, -1.0]
        eps[[0, 1, 2, 3], [0, 1, 2, 3], 1, 0, 0] = [0.5, 0.5, 0.5, 0.5]
        ev, evec = canonical_band_eigen(eps, 2, conserving=True)
        np.testing.assert_array_equal(ev[0, 0, 0], [-1.0, 0.0, 1.0, 1.0])
        np.testing.assert_array_equal(np.argmax(np.abs(evec[0, 0, 0]), axis=0), [3, 1, 0, 2])
        np.testing.assert_array_equal(np.argmax(np.abs(evec[1, 0, 0]), axis=0), [0, 1, 2, 3])

    def test_green_is_canonical(self):
        for fixture in (FIXTURE_SZ_GRID, FIXTURE_LZSZ):
            for h in FIELDS:
                with self.subTest(norb=fixture.norb_phys, h=h):
                    norb = fixture.norb_phys
                    eps = _epsilon_k(fixture, h)
                    G = _green(fixture, canonical_band_eigen(eps, norb, conserving=True))
                    self.assertTrue(np.all(_off(G, norb) == 0))
                    self.assertEqual(green_leak(G, norb), 0.0)
                    ref = canonical_green(fixture._replace(zeeman=h), NMAT)
                    self.assertLessEqual(float(np.max(np.abs(G[0] - ref))), 1e-13)

    def test_full_diagonalization_leaks(self):
        """The interleaved full diagonalization of the same band leaks S_z at
        round-off: the reason for the canonical rule. (norb_phys = 2: with one
        orbital the interleaved and the spin-block orders coincide.)
        Measured at beta = 2: 2.7e-16 (h = 0) and 2.5e-16 (h = 0.1)."""
        fixture, norb = FIXTURE_LZSZ, FIXTURE_LZSZ.norb_phys
        for h in FIELDS:
            with self.subTest(h=h):
                eps = spin_orbital.spinblock_to_interleaved(_epsilon_k(fixture, h), norb,
                                                            axes=(0, 1))
                G_il = _green(fixture, sc_band._calc_eigenvalues(eps))
                G = spin_orbital.interleaved_to_spinblock(G_il, norb, axes=(-2, -1))
                leak = green_leak(G, norb)
                self.assertGreater(leak, 0.0)
                self.assertLessEqual(leak, 1e-12)

    def test_green_does_not_read_the_layout(self):
        """The forced-dense and the automatic layout give the same Green
        function (to round-off of the recomputation); so does a PairLift
        interaction on a spin-diagonal band, whose layout is dense while the
        band conserves S_z."""
        fixture = FIXTURE_SZ_GRID
        tables = fixture.tables_dict
        pairlift = {"CoulombIntra": {((0, 0, 0), (0, 0)): BASE["CoulombIntra"]},
                    "PairLift": closed_row(1, 0, 0, BASE["PairLift"])}
        for h in FIELDS:
            with self.subTest(h=h):
                g_auto, l_auto = _band_green(fixture, h, tables, "auto")
                g_dense, l_dense = _band_green(fixture, h, tables, "dense")
                g_pl, l_pl = _band_green(fixture, h, pairlift, "auto")
                self.assertEqual((l_auto, l_dense, l_pl), ("sz_packed", "dense", "dense"))
                assert_recomputed(self, g_dense, g_auto)
                assert_recomputed(self, g_pl, g_auto)
                self.assertTrue(np.all(_off(g_pl, 1) == 0))

    def test_not_conserving_is_full_diagonalization(self):
        eps = _epsilon_k(FIXTURE_LZSZ, 0.1)
        eps = eps.copy()
        eps[0, 2] += 0.05
        eps[2, 0] += 0.05
        ev, evec = canonical_band_eigen(eps, 2, conserving=False)
        ref_ev, ref_evec = sc_band._calc_eigenvalues(eps)
        assert_recomputed(self, ev, ref_ev)
        assert_recomputed(self, evec, ref_evec)


class TestCanonicalBandEigenRefusals(unittest.TestCase):
    """Construction check: refusals of canonical_band_eigen."""

    def test_interleaved_band_refused(self):
        eps = spin_orbital.spinblock_to_interleaved(_epsilon_k(FIXTURE_LZSZ, 0.1), 2,
                                                    axes=(0, 1))
        up, dn = _halves(2)
        off = max(float(np.max(np.abs(eps[up, dn]))), float(np.max(np.abs(eps[dn, up]))))
        self.assertGreater(off, 0.0)
        with self.assertRaisesRegex(ValueError, exact(INTERLEAVED % off)):
            canonical_band_eigen(eps, 2, conserving=True)

    def test_non_finite_band_refused(self):
        """A non-finite entry is refused as such, also where it would read
        as an S_z-off-sector entry."""
        for where in ((0, 0), (0, 2)):
            with self.subTest(entry=where):
                eps = _epsilon_k(FIXTURE_LZSZ, 0.1).copy()
                eps[where + (0, 0, 0)] = np.nan
                with self.assertRaisesRegex(ValueError, exact(
                        "canonical_band_eigen: epsilon_k has non-finite entries")):
                    canonical_band_eigen(eps, 2, conserving=True)

    def test_shape_and_flag(self):
        eps = _epsilon_k(FIXTURE_SZ_GRID, 0.0)
        with self.assertRaisesRegex(ValueError, "nd, nd, Nx, Ny, Nz"):
            canonical_band_eigen(eps[..., 0], 1, conserving=True)
        with self.assertRaisesRegex(ValueError, "2 \\* norb_phys"):
            canonical_band_eigen(eps, 2, conserving=True)
        check = band_conserves_sz(FIXTURE_SZ_GRID.transfer_dict, None, 1)
        with self.assertRaisesRegex(ValueError, "conserving must be a bool"):
            canonical_band_eigen(eps, 1, conserving=check)


def _packed_bubble(fixture, G):
    """The packed static bond bubble of ``G`` (``so_bond_bubble``)."""
    c = compile_interaction(fixture.tables_dict, fixture.norb_phys)
    topo = so_bond_topology(c, CELL)
    layout, reason = select_layout(c, fixture.transfer_dict, override="auto")
    sm = sz_sector_map(topo, fixture.nd)
    plan = plan_so_bond(nmat=NMAT, nvol=fixture.nvol, B=len(topo.delta_r), nd=fixture.nd,
                        consumers=("chi_static",), caps=SOBondCaps(None, None, None),
                        layout=layout, layout_reason=reason)
    sink = PackedArrays.zeros(sm, (fixture.nvol,))
    so_bond_bubble(np.array(G), None, fixture.beta, topo, fixture.spatial_shape, plan=plan,
                   freq="static", sink=sink)
    return sink


def _dense_bubble(fixture, G):
    """The dense static bond bubble of ``G`` and its S_z-off-sector mask."""
    c = compile_interaction(fixture.tables_dict, fixture.norb_phys)
    topo = so_bond_topology(c, CELL)
    B, nd = len(topo.delta_r), fixture.nd
    ND = B * nd * nd
    plan = plan_so_bond(nmat=NMAT, nvol=fixture.nvol, B=B, nd=nd, consumers=("chi_static",),
                        caps=SOBondCaps(None, None, None))
    sink = np.zeros((fixture.nvol, ND, ND), dtype=np.complex128)
    so_bond_bubble(np.array(G), None, fixture.beta, topo, fixture.spatial_shape, plan=plan,
                   freq="static", sink=sink)
    so = sz_sector_map(topo, nd).sector_of
    return sink, so[:, None] != so[None, :]


class TestBorrowedGreen(unittest.TestCase):
    """Gate P-G (ii)/(iii): a borrowed Green function with a round-off
    S_z leak is projected (one WARNING), the results equal the clean run;
    a 1e-6 leak is refused."""

    @classmethod
    def setUpClass(cls):
        cls.fixture = FIXTURE_SZ_GRID
        eps = _epsilon_k(cls.fixture, 0.1)
        cls.clean = _green(cls.fixture, canonical_band_eigen(eps, 1, conserving=True))
        cls.clean.flags.writeable = False

    def _perturbed(self, rel):
        G = np.array(self.clean)
        rng = np.random.default_rng(0)
        z = rng.normal(size=G[..., 0:1, 1:2].shape) + 1j * rng.normal(size=G[..., 0:1, 1:2].shape)
        scale = rel * float(np.max(np.abs(G)))
        G[..., 0:1, 1:2] += scale * z / np.max(np.abs(z))
        G[..., 1:2, 0:1] += scale * np.conj(z) / np.max(np.abs(z))
        return G

    def test_clean_used_as_is(self):
        with self.assertNoLogs(so_bond.logger, level="WARNING"):
            used, leak = project_borrowed_green(self.clean, norb_phys=1, source="g.npz")
        self.assertIs(used, self.clean)
        self.assertEqual(leak, 0.0)

    def test_round_off_leak_projected(self):
        G = self._perturbed(1e-14)
        leak = green_leak(G, 1)
        self.assertTrue(0.0 < leak <= 1e-12)
        with self.assertLogs(so_bond.logger, level="WARNING") as cm:
            used, got = project_borrowed_green(G, norb_phys=1, source="g.npz")
        self.assertEqual(cm.output, ["WARNING:%s:%s" % (so_bond.logger.name,
                                                         WARN % ("g.npz", leak))])
        self.assertEqual(got, leak)
        self.assertTrue(np.all(_off(used, 1) == 0))
        np.testing.assert_array_equal(used, self.clean)
        self.assertEqual(green_digest(used, LAYOUT), green_digest(self.clean, LAYOUT))
        self.assertTrue(np.all(_off(G, 1) != 0))         # the input is not modified
        a = _packed_bubble(self.fixture, used)
        b = _packed_bubble(self.fixture, self.clean)
        scale = max(float(np.max(np.abs(x))) for x in b.blocks)
        for x, y in zip(a.blocks, b.blocks):
            self.assertLessEqual(float(np.max(np.abs(x - y))), 1e-12 * scale)

    def test_projection_is_not_vacuous(self):
        """With the projection disabled the perturbed Green function reaches
        the bubble: the dense bubble has S_z-off-sector entries and the
        packed bubble refuses it."""
        G = self._perturbed(1e-14)
        with mock.patch.object(so_bond, "project_green_sectors",
                               lambda g, norb_phys, *, axes=(-2, -1): g), quiet():
            used, _ = project_borrowed_green(G, norb_phys=1, source="g.npz")
        self.assertIs(used, G)
        dense, mask = _dense_bubble(self.fixture, used)
        off = max(float(np.max(np.abs(dense[q][mask]))) for q in range(dense.shape[0]))
        self.assertGreater(off, 1e-3 * 1e-14 * float(np.max(np.abs(G))))
        with self.assertRaisesRegex(ValueError, "static bond bubble has S_z-off-sector"):
            _packed_bubble(self.fixture, used)

    def test_incompatible_refused(self):
        G = self._perturbed(1e-6)
        leak = green_leak(G, 1)
        with self.assertRaisesRegex(ValueError, exact(REFUSE % ("g.npz", leak))):
            project_borrowed_green(G, norb_phys=1, source="g.npz")


class TestGreenSectorHelpers(unittest.TestCase):
    """Construction check: green_leak and project_green_sectors on other
    axes, the zero array and non-finite input."""

    def test_axes(self):
        rng = np.random.default_rng(1)
        G = rng.normal(size=(4, 3, 4, 2)) + 0j          # spin-orbital axes 0 and 2
        moved = np.moveaxis(G, (0, 2), (-2, -1))
        ref = float(np.max(np.abs(_off(moved, 2)))) / float(np.max(np.abs(G)))
        self.assertEqual(green_leak(G, 2, axes=(0, 2)), ref)
        P = project_green_sectors(G, 2, axes=(0, 2))
        self.assertIsNot(P, G)
        pm = np.moveaxis(P, (0, 2), (-2, -1))
        self.assertTrue(np.all(_off(pm, 2) == 0))
        up, dn = _halves(2)
        np.testing.assert_array_equal(pm[..., up, up], moved[..., up, up])
        np.testing.assert_array_equal(pm[..., dn, dn], moved[..., dn, dn])
        self.assertEqual(green_leak(P, 2, axes=(0, 2)), 0.0)

    def test_zero_and_bad_input(self):
        self.assertEqual(green_leak(np.zeros((2, 2)), 1), 0.0)
        with self.assertRaisesRegex(ValueError, "2 \\* norb_phys"):
            green_leak(np.zeros((3, 3)), 1)
        G = np.zeros((2, 2), dtype=complex)
        G[0, 0] = np.nan
        self.assertTrue(np.isnan(green_leak(G, 1)))
        with self.assertRaisesRegex(ValueError, "non-finite"):
            project_borrowed_green(G, norb_phys=1, source="g.npz")


if __name__ == "__main__":
    unittest.main()
