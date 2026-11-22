"""Tests of the UHFk readers of ``hwave.results.uhf``: the bands, the Green
function, the order parameters and the spin-orbital density helper."""
import shutil
import unittest
from fractions import Fraction

import numpy as np

import hwave.results as results
from hwave.results import (
    Bands, GreenR, MissingOutputError, UnsupportedLayoutError, open_run,
    read_bands, read_energy, read_onebodyg, read_scf_status, read_uhfk_green,
    read_uhfk_order_parameters, so_density_from_green, sum_rule_tolerance)
from tests.hf_square_oracle import neel_gap_equation
from tests.results_fixtures import (
    SO_XY_FIELDS, TAU_CENTRE, TAU_RHO, assert_deep_frozen, copy_fixture,
    fixture, write_uhfk_geometry)

HALF = Fraction(1, 2)
ZERO = Fraction(0)

SIGMA = (np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]]),
         np.array([[1, 0], [0, -1]]))


def spin_of_slot(rho):
    """``S_a = 1/2 sum_{st} (sigma_a)_{st} rho_{st}`` of a Green slot
    ``rho[s, t] = <c^dag_s c_t>``, written out independently of the
    package."""
    return np.array([0.5 * np.sum(s * rho).real for s in SIGMA])


def uhfk_sum_rule_bound(toml):
    """``sum_rule_tolerance`` of a UHFk run that performed SCF steps: the
    Green function has ``nvol * nd^2`` elements, ``nd = ns * norb *
    subvol``."""
    run = open_run(toml)
    status = read_scf_status(run)
    if not status.iterations:
        raise AssertionError("the run performed no SCF step")
    green = read_uhfk_green(run)
    member = green.sublattice if green.sublattice is not None else green
    nvol, ns, norb = member.green.shape[:3]
    nd = ns * norb
    return sum_rule_tolerance(status.residual, mix=run.param("Mix"),
                              g_size=nvol * nd * nd, nvol=nvol, nd=nd)


def q_index(op, q):
    return op.q_frac.index(tuple(Fraction(c) for c in q))


class TestExports(unittest.TestCase):
    def test_names(self):
        for name in ("Bands", "read_bands", "GreenR", "read_uhfk_green",
                     "read_uhfk_order_parameters", "OneBodyG", "read_onebodyg",
                     "so_density_from_green"):
            self.assertIn(name, results.__all__)
            self.assertTrue(hasattr(results, name))


class TestSoHelper(unittest.TestCase):
    def test_so_helper_index_map(self):
        norb_phys, subvol = 2, 2
        nd = 2 * norb_phys * subvol
        p, q = np.meshgrid(np.arange(nd), np.arange(nd), indexing="ij")
        g0 = 100 * p + q + 1j * (p - q)
        rho = so_density_from_green(g0, norb_phys, subvol)
        self.assertEqual(rho.shape, (norb_phys * subvol, 2, 2))
        for j in range(norb_phys * subvol):
            for s in range(2):
                for t in range(2):
                    self.assertEqual(rho[j, s, t], g0[2 * j + s, 2 * j + t])

    def test_so_helper_shape_checked(self):
        with self.assertRaisesRegex(ValueError, "G0"):
            so_density_from_green(np.zeros((6, 6)), 2, 2)
        with self.assertRaisesRegex(ValueError, "G0"):
            so_density_from_green(np.zeros((8, 7)), 2, 2)
        with self.assertRaisesRegex(ValueError, "norb_phys_orig"):
            so_density_from_green(np.zeros((8, 8)), 0, 2)


class TestGroupedOrbitals(unittest.TestCase):
    """``uhfk_2orb_tau``: two orbitals at the same fractional centre
    (0.25, 0.5, 0), CellShape = SubShape = [2, 1, 1], a synthetic
    green.npz."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfk_2orb_tau")

    def expected_cell_spins(self, ir):
        return sum(spin_of_slot(TAU_RHO[(ir, k)]) for k in range(2))

    def test_grouped_multi_orbital_fractional_tau(self):
        op = read_uhfk_order_parameters(self.toml, site_of_orbital=[0, 0])
        np.testing.assert_allclose(op.site_positions,
                                   [(0.25, 0.5, 0.0), (1.25, 0.5, 0.0)], atol=0)
        spins = np.array([self.expected_cell_spins(ir) for ir in range(2)])
        np.testing.assert_allclose(op.spin, spins, atol=1e-14)
        charges = [sum(np.trace(TAU_RHO[(ir, k)]).real for k in range(2))
                   for ir in range(2)]
        np.testing.assert_allclose(op.charge, charges, atol=1e-14)
        # M(q) = N_site^-1 sum_s S_s e^{2 pi i q.r_s}, q = (1/2, 0, 0)
        positions = np.array([(0.25, 0.5, 0.0), (1.25, 0.5, 0.0)])
        expected = sum(spins[s] * np.exp(2j * np.pi * 0.5 * positions[s, 0])
                       for s in range(2)) / 2
        self.assertGreater(np.max(np.abs(expected.imag)), 1e-3)   # complex
        np.testing.assert_allclose(op.spin_fourier[q_index(op, (HALF, 0, 0))],
                                   expected, atol=1e-14)
        self.assertEqual(op.cell_shape, (2, 1, 1))
        self.assertEqual(op.lattice_shape, (2, 1, 1))
        self.assertEqual(op.n_cell, 1)
        self.assertEqual(op.source.name, "green.npz")
        assert_deep_frozen(op)

    def test_without_grouping_every_orbital_is_a_site(self):
        op = read_uhfk_order_parameters(self.toml)
        self.assertEqual(op.spin.shape, (4, 3))
        np.testing.assert_array_equal(op.site_of_orbital, [0, 1, 2, 3])
        np.testing.assert_allclose(
            op.site_positions,
            [TAU_CENTRE, TAU_CENTRE, (1.25, 0.5, 0.0), (1.25, 0.5, 0.0)], atol=0)
        # entry j = k + 2 ir
        for ir in range(2):
            for k in range(2):
                np.testing.assert_allclose(op.spin[k + 2 * ir],
                                           spin_of_slot(TAU_RHO[(ir, k)]),
                                           atol=1e-14)
        assert_deep_frozen(op)

    def test_site_positions_replace_the_centres(self):
        op = read_uhfk_order_parameters(self.toml, site_of_orbital=[0, 0],
                                        site_positions=[[0.1, 0.2, 0.0]])
        np.testing.assert_allclose(op.site_positions,
                                   [(0.1, 0.2, 0.0), (1.1, 0.2, 0.0)], atol=1e-15)
        with self.assertRaisesRegex(ValueError, "site_positions"):
            read_uhfk_order_parameters(self.toml, site_of_orbital=[0, 0],
                                       site_positions=[[0.1, 0.2, 0.0]] * 2)

    def test_distinct_centres_cannot_be_grouped(self):
        toml = copy_fixture("uhfk_2orb_tau")
        write_uhfk_geometry(toml.parent / "geom.dat",
                            [TAU_CENTRE, (0.75, 0.5, 0.0)])
        with self.assertRaisesRegex(ValueError, r"site 0.*orbitals \[0, 1\]"):
            read_uhfk_order_parameters(toml, site_of_orbital=[0, 0])
        # as separate sites, or with explicit positions, they are fine
        op = read_uhfk_order_parameters(toml)
        np.testing.assert_allclose(op.site_positions[1], (0.75, 0.5, 0.0))
        read_uhfk_order_parameters(toml, site_of_orbital=[0, 0],
                                   site_positions=[TAU_CENTRE])

    def test_bad_site_of_orbital(self):
        for bad in ([0, 2], [0], [0, 1, 1], [-1, 0], [0.5, 0]):
            with self.subTest(bad=bad), \
                    self.assertRaisesRegex(ValueError, "site_of_orbital"):
                read_uhfk_order_parameters(self.toml, site_of_orbital=bad)

    def test_green_member_layouts(self):
        green = read_uhfk_green(self.toml)
        self.assertIsInstance(green, GreenR)
        self.assertEqual(green.green.shape, (2, 2, 2, 2, 2))
        self.assertEqual(green.sublattice.green.shape, (1, 2, 4, 2, 4))
        self.assertEqual(green.cell_shape, (2, 1, 1))
        self.assertEqual(green.sublattice.cell_shape, (1, 1, 1))
        np.testing.assert_array_equal(green.r_index, [(0, 0, 0), (1, 0, 0)])
        self.assertFalse(green.spin_orbital)
        assert_deep_frozen(green)


class TestGeometryRequired(unittest.TestCase):
    def test_missing_geometry_file(self):
        toml = copy_fixture("uhfk_u8")
        (toml.parent / "geom.dat").unlink()
        with self.assertRaises(MissingOutputError) as cm:
            read_uhfk_order_parameters(toml)
        self.assertEqual(cm.exception.path.name, "geom.dat")

    def test_geometry_norb_must_match_green(self):
        toml = copy_fixture("uhfk_u8")
        write_uhfk_geometry(toml.parent / "geom.dat", [(0, 0, 0), (0.5, 0, 0)])
        with self.assertRaisesRegex(UnsupportedLayoutError, "Geometry"):
            read_uhfk_order_parameters(toml)

    def test_so_entries_of_one_orbital_share_a_centre(self):
        toml = copy_fixture("uhfk_so_xy")
        write_uhfk_geometry(toml.parent / "geom.dat",
                            [(0, 0, 0), (0.5, 0, 0), (0, 0, 0), (0, 0, 0)])
        with self.assertRaisesRegex(UnsupportedLayoutError, "entries 0 and 1"):
            read_uhfk_order_parameters(toml)


class TestSpinOrbital(unittest.TestCase):
    """``uhfk_so_xy``: two physical orbitals in on-site fields along n_0
    (in the x-y plane) and n_1 (in the x-z plane)."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfk_so_xy")

    def test_so_field_directions_phase(self):
        op = read_uhfk_order_parameters(self.toml)
        self.assertEqual(op.spin.shape, (2, 3))
        for k in range(2):
            for a in range(3):
                with self.subTest(orbital=k, axis=a):
                    self.assertAlmostEqual(op.spin[k, a], 0.5 * SO_XY_FIELDS[k][a],
                                           delta=1e-8)
        energy = read_energy(self.toml)
        tol = uhfk_sum_rule_bound(self.toml)
        totals = 2 * op.n_cell * op.uniform
        for value, total in zip(totals, (energy.sx, energy.sy, energy.sz)):
            self.assertLessEqual(abs(value - total), tol)
        self.assertEqual(op.meta["solver"], "UHFk")
        self.assertTrue(op.meta["spin_orbital"])
        assert_deep_frozen(op)

    def test_green_and_bands(self):
        green = read_uhfk_green(self.toml)
        self.assertTrue(green.spin_orbital)
        self.assertEqual(green.green.shape, (1, 1, 4, 1, 4))
        self.assertIsNone(green.sublattice)
        bands = read_bands(self.toml)
        self.assertTrue(bands.spin_orbital)
        self.assertEqual(bands.energies.shape, (1, 4))


class TestNeelU8(unittest.TestCase):
    """``uhfk_u8``: tests/uhfk/CoulombIntra (4 x 4, U = 8, SubShape 2 x 2)
    from its own start: a collinear Neel state with a tilted axis."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfk_u8")
        cls.m, cls.e = neel_gap_equation(4, 8)

    def test_oracle_self_check(self):
        self.assertAlmostEqual(self.m, 0.4474190440247483, delta=1e-12)
        self.assertAlmostEqual(self.e, -7.38962582304832, delta=1e-10)
        m, e = neel_gap_equation(16, 2)
        self.assertAlmostEqual(m, 0.1917429137967485, delta=1e-12)
        self.assertAlmostEqual(e, -291.4365477135228, delta=1e-9)

    def test_golden_neel_u8(self):
        op = read_uhfk_order_parameters(self.toml)
        self.assertAlmostEqual(op.staggered_norm, self.m, delta=1e-8)
        self.assertAlmostEqual(read_energy(self.toml).total, self.e, delta=1e-8)
        self.assertLess(op.uniform_norm, 1e-10)
        # the axis is tilted: Sz alone is not the observable
        np.testing.assert_allclose(np.abs(op.staggered.real),
                                   (0.113051, 0.0, 0.432901), atol=1e-5)
        self.assertEqual(op.spin_ordering_vector, (HALF, HALF, ZERO))
        self.assertEqual(op.cell_shape, (2, 2, 1))
        self.assertEqual(op.lattice_shape, (4, 4, 1))
        self.assertEqual(op.n_cell, 4)
        np.testing.assert_allclose(op.site_positions,
                                   [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)])
        self.assertEqual(op.meta["green_member"], "green_sublattice")
        assert_deep_frozen(op)

    def test_cross_solver(self):
        from hwave.results import read_uhfr_order_parameters
        positions = np.array([(i % 4, i // 4, 0) for i in range(16)], dtype=float)
        uhfr = read_uhfr_order_parameters(fixture("uhfr_u8"),
                                          site_positions=positions,
                                          lattice_shape=(4, 4, 1))
        uhfk = read_uhfk_order_parameters(self.toml)
        self.assertAlmostEqual(uhfr.staggered_norm, uhfk.staggered_norm, delta=1e-8)

    def test_green_layout(self):
        green = read_uhfk_green(self.toml)
        self.assertEqual(green.green.shape, (16, 2, 1, 2, 1))
        self.assertEqual(green.cell_shape, (4, 4, 1))
        # r = rz + Nz (ry + Ny rx)
        self.assertEqual(tuple(green.r_index[1]), (0, 1, 0))
        self.assertEqual(tuple(green.r_index[4]), (1, 0, 0))
        sub = green.sublattice
        self.assertEqual(sub.green.shape, (4, 2, 4, 2, 4))
        self.assertEqual(sub.cell_shape, (2, 2, 1))
        self.assertEqual(tuple(sub.r_index[2]), (1, 0, 0))
        # the deflated member at R = (1, 0, 0) is the sublattice entry
        # aa = a + norb * ir with ir = irx + Bx (iry + By irz) = 1
        np.testing.assert_allclose(green.green[4, :, 0, :, 0],
                                   sub.green[0, :, 1, :, 0], atol=0)
        np.testing.assert_allclose(green.green[1, :, 0, :, 0],
                                   sub.green[0, :, 2, :, 0], atol=0)
        self.assertEqual(green.meta["legacy_notes"], [])
        assert_deep_frozen(green)

    def test_bands(self):
        bands = read_bands(self.toml, with_vectors=True)
        self.assertIsInstance(bands, Bands)
        self.assertEqual(bands.energies.shape, (4, 8))
        np.testing.assert_allclose(bands.k_cart, bands.k_index @ bands.k_unit)
        for k in range(4):
            v = bands.eigenvectors[k]
            np.testing.assert_allclose(v.conj().T @ v, np.eye(8), atol=1e-10)
        self.assertTrue(bands.single_reservoir)
        self.assertAlmostEqual(bands.mu, 4.0, delta=1e-8)
        self.assertIsNone(bands.mu_groups)
        self.assertFalse(bands.spin_orbital)
        assert_deep_frozen(bands)
        self.assertIsNone(read_bands(self.toml).eigenvectors)


class TestRotationInvariance(unittest.TestCase):
    def test_rotation_invariance(self):
        """``uhfk_u8_random``: from a random start the Neel state comes out
        along another axis, with the same staggered moment."""
        m, _ = neel_gap_equation(4, 8)
        op = read_uhfk_order_parameters(fixture("uhfk_u8_random"))
        self.assertGreater(read_scf_status(fixture("uhfk_u8_random")).iterations, 1)
        self.assertAlmostEqual(op.staggered_norm, m, delta=1e-8)
        golden = read_uhfk_order_parameters(fixture("uhfk_u8"))
        axis = np.abs(op.staggered.real) / op.staggered_norm
        golden_axis = np.abs(golden.staggered.real) / golden.staggered_norm
        self.assertGreater(np.linalg.norm(axis - golden_axis), 1e-3)


class TestLargeLattice(unittest.TestCase):
    def test_u2_l16(self):
        toml = fixture("uhfk_u2_l16")
        m, e = neel_gap_equation(16, 2)
        op = read_uhfk_order_parameters(toml)
        self.assertAlmostEqual(op.staggered_norm, m, delta=1e-7)
        self.assertAlmostEqual(read_energy(toml).total, e, delta=1e-7)
        self.assertEqual(op.n_cell, 64)


class TestSumRules(unittest.TestCase):
    def check_sum_rules(self, name):
        toml = fixture(name)
        op = read_uhfk_order_parameters(toml)
        energy = read_energy(toml)
        tol = uhfk_sum_rule_bound(toml)
        n_site = op.charge.shape[0]
        self.assertLessEqual(abs(op.n_cell * op.charge.sum() - energy.ncond), tol)
        totals = n_site * op.n_cell * op.uniform
        for value, total in zip(totals, (energy.sx, energy.sy, energy.sz)):
            self.assertLessEqual(abs(value - total), tol)
        return op, energy, tol

    def test_sum_rules_within_residual_bound(self):
        for name in ("uhfk_u8", "uhfk_cdw"):
            with self.subTest(name=name):
                self.check_sum_rules(name)

    def test_ferromagnet_axis(self):
        """``uhfk_fm_axis``: a saturated ferromagnet along the generic axis
        of its start."""
        op, energy, tol = self.check_sum_rules("uhfk_fm_axis")
        self.assertEqual(open_run(fixture("uhfk_fm_axis")).param("Mix"), 0.5)
        u = np.array([0.6, 0.2537, 0.7734])
        u = u / np.linalg.norm(u)
        np.testing.assert_allclose(op.uniform / op.uniform_norm, u, atol=1e-8)
        self.assertEqual(op.n_cell, 16)
        self.assertEqual(op.meta["green_member"], "green")

    def test_u0_finite_t(self):
        """``uhfk_u0_t005``: at U = 0 the Hamiltonian does not depend on
        spin, so no spin order may be reported (a guard against false
        positives of the reader). T > 0 avoids the degenerate Fermi shell of
        T = 0, where an antiferromagnetic start ends at Sz_total = 1."""
        toml = fixture("uhfk_u0_t005")
        op = read_uhfk_order_parameters(toml)
        self.assertLess(op.uniform_norm, 1e-8)
        self.assertLess(op.staggered_norm, 1e-8)
        self.assertEqual(read_scf_status(toml).status, "converged")


class TestChargeOrder(unittest.TestCase):
    """``uhfk_cdw``: docs/tutorial/Hubbard/UHFk (V = 4)."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfk_cdw")

    def test_cdw(self):
        op = read_uhfk_order_parameters(self.toml)
        np.testing.assert_allclose(sorted(op.charge, reverse=True),
                                   [1.990014, 1.990014, 0.009986, 0.009986],
                                   atol=1e-5)
        self.assertAlmostEqual(op.charge_order_amplitude, 0.990014, delta=1e-5)
        self.assertEqual(op.charge_ordering_vector, (HALF, HALF, ZERO))
        bands = read_bands(self.toml)
        self.assertAlmostEqual(bands.homo, 0.1598, delta=1e-3)
        self.assertAlmostEqual(bands.lumo, 31.8402, delta=1e-3)
        self.assertAlmostEqual(bands.gap, 31.6804, delta=1e-3)
        self.assertAlmostEqual(bands.gap, bands.lumo - bands.homo, delta=1e-12)
        self.assertAlmostEqual(bands.mu, 16.0, delta=1e-6)
        self.assertIn("not a band gap", bands.meta["homo_lumo_note"])

    def test_bands_raw_order_and_stable_sort(self):
        bands = read_bands(self.toml)
        with np.load(bands.source) as data:
            np.testing.assert_array_equal(bands.energies, data["eigenvalue"])
        ordered = np.take_along_axis(bands.energies, bands.sort_index, 1)
        self.assertTrue(np.all(np.diff(ordered, axis=1) >= 0))
        # a tie keeps the archive order
        toml = copy_fixture("uhfk_cdw")
        path = toml.parent / "output" / "eigen.npz"
        with np.load(path) as data:
            members = {k: data[k] for k in data.files}
        members["eigenvalue"] = np.tile([1.0, 0.0, 1.0, 0.0, 2.0, 0.0, -1.0, 1.0],
                                        (4, 1))
        np.savez(path, **members)
        tied = read_bands(toml)
        np.testing.assert_array_equal(tied.sort_index[0], [6, 1, 3, 5, 0, 2, 7, 4])

    def test_non_integer_ncond_gives_no_homo(self):
        toml = copy_fixture("uhfk_cdw")
        path = toml.parent / "output" / "energy.dat"
        path.write_text(path.read_text().replace("NCond   = 16.0",
                                                 "NCond   = 15.5"))
        bands = read_bands(toml)
        self.assertIsNone(bands.homo)
        self.assertIsNone(bands.lumo)
        self.assertIsNone(bands.gap)
        self.assertIn("NCond", bands.meta["homo_lumo_unavailable"])


class TestFixedSz(unittest.TestCase):
    def test_2sz_fixed(self):
        bands = read_bands(fixture("uhfk_2sz_fixed"))
        self.assertIs(bands.single_reservoir, False)
        self.assertIsNone(bands.mu)
        self.assertEqual(list(bands.mu_groups), [0, 1])
        self.assertAlmostEqual(bands.mu_groups[0], 2.9930814020896315, delta=1e-12)
        self.assertIsNone(bands.homo)
        self.assertIsNone(bands.lumo)
        self.assertIsNone(bands.gap)
        self.assertIn("homo_lumo_unavailable", bands.meta)
        assert_deep_frozen(bands)


class TestLegacy(unittest.TestCase):
    def test_legacy_tutorial_green(self):
        toml = fixture("uhfk_cdw_legacy")
        with self.assertRaisesRegex(UnsupportedLayoutError, r"legacy=True"):
            read_uhfk_green(toml)
        with self.assertRaisesRegex(UnsupportedLayoutError, r"legacy=True"):
            read_uhfk_order_parameters(toml)
        green = read_uhfk_green(toml, legacy=True)
        self.assertIsNone(green.green)
        self.assertEqual(green.sublattice.green.shape, (4, 2, 4, 2, 4))
        members = [note["member"] for note in green.meta["legacy_notes"]]
        self.assertEqual(members, ["momentum_convention", "green"])
        op = read_uhfk_order_parameters(toml, legacy=True)
        self.assertEqual(op.meta["green_member"], "green_sublattice")
        self.assertEqual(op.meta["legacy_notes"], green.meta["legacy_notes"])
        self.assertAlmostEqual(op.charge_order_amplitude, 0.99, delta=1e-2)
        assert_deep_frozen(green)

    def test_legacy_eigen(self):
        toml = copy_fixture("uhfk_cdw")
        path = toml.parent / "output" / "eigen.npz"
        with np.load(path) as data:
            members = {k: data[k] for k in data.files if k != "momentum_convention"}
        np.savez(path, **members)
        with self.assertRaisesRegex(UnsupportedLayoutError, r"legacy=True"):
            read_bands(toml)
        bands = read_bands(toml, legacy=True)
        self.assertEqual(bands.meta["legacy_notes"][0]["assumed"], "e_plus_ikR")


class TestZeroSteps(unittest.TestCase):
    def test_left_over_outputs_refused(self):
        toml = copy_fixture("uhfk_u8_zero_steps")
        shutil.rmtree(toml.parent / "output")
        shutil.copytree(fixture("uhfk_u8").parent / "output", toml.parent / "output")
        (toml.parent / "output" / "scf_status.dat").write_text(
            (fixture("uhfk_u8_zero_steps").parent / "output" / "scf_status.dat")
            .read_text())
        for reader in (read_bands, read_uhfk_green, read_uhfk_order_parameters):
            with self.subTest(reader=reader.__name__), \
                    self.assertRaisesRegex(MissingOutputError, "not from this run"):
                reader(toml)


class TestOneBodyG(unittest.TestCase):
    def test_onebodyg_known_issue(self):
        toml = copy_fixture("uhfk_u8")
        text = toml.read_text()
        self.assertIn("[file.output]\n", text)
        toml.write_text(text.replace("[file.output]\n",
                                     '[file.output]\n  onebodyg = "greenone.dat"\n'))
        (toml.parent / "output" / "greenone.dat").write_text(
            "  0   0   0   0  5.0e-01 0.0e+00\n")
        with self.assertRaisesRegex(UnsupportedLayoutError, "#357") as cm:
            read_onebodyg(toml)
        self.assertIn("#357", cm.exception.hint)
        self.assertEqual(cm.exception.path.name, "greenone.dat")

    def test_onebodyg_not_configured(self):
        with self.assertRaises(MissingOutputError):
            read_onebodyg(fixture("uhfk_u8"))

    def test_uhfr_run_refused(self):
        for reader in (read_bands, read_uhfk_green, read_uhfk_order_parameters,
                       read_onebodyg):
            with self.subTest(reader=reader.__name__), \
                    self.assertRaisesRegex(UnsupportedLayoutError, "not a UHFk run"):
                reader(fixture("uhfr_tutorial"))


if __name__ == "__main__":
    unittest.main()
