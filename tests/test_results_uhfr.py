"""Tests of the UHFr readers of ``hwave.results.uhf``: the eigenstates, the
Green function and the order parameters."""
import re
import unittest

import numpy as np

import hwave.results as results
from hwave.results import (
    MissingOutputError, UHFrEigen, UHFrGreen, UnsupportedLayoutError, open_run,
    read_energy, read_scf_status, read_uhfr_eigen, read_uhfr_green,
    read_uhfr_order_parameters, sum_rule_tolerance)
from tests.hf_square_oracle import neel_gap_equation
from tests.results_fixtures import (
    REPO, assert_deep_frozen, copy_fixture, fixture, run_input, write_greenone)

#: the 4 x 4 cluster of tests/uhfr/coulombintra: site i at (i mod 4, i // 4)
#: in units of the lattice vectors
U8_POSITIONS = np.array([(i % 4, i // 4, 0) for i in range(16)], dtype=float)
U8_LATTICE = (4, 4, 1)


def sum_rule_bound(run, n_site):
    """``sum_rule_tolerance`` of a UHFr run that performed SCF steps."""
    status = read_scf_status(run)
    if not status.iterations:
        raise AssertionError("the run performed no SCF step; there is no "
                             "residual to bound the sum rules with")
    return sum_rule_tolerance(status.residual, mix=run.param("Mix"),
                              g_size=4 * n_site ** 2, nvol=1, nd=2 * n_site)


def green_lines(path):
    return [line for line in path.read_text().splitlines() if line.strip()]


class TestExports(unittest.TestCase):
    def test_names(self):
        for name in ("UHFrEigen", "read_uhfr_eigen", "UHFrGreen",
                     "read_uhfr_green", "read_uhfr_order_parameters"):
            self.assertIn(name, results.__all__)
            self.assertTrue(hasattr(results, name))


class TestTutorial(unittest.TestCase):
    """docs/tutorial/Hubbard/UHFr: Nsite 8, 2Sz = 0, the 16 diagonal
    entries in greenone.def."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfr_tutorial")
        cls.hwave_run = open_run(cls.toml)

    def test_eigen_blocks(self):
        eigen = read_uhfr_eigen(self.toml)
        self.assertIsInstance(eigen, UHFrEigen)
        self.assertEqual(list(eigen.energies), ["spin-up", "spin-down"])
        for block in ("spin-up", "spin-down"):
            self.assertEqual(eigen.energies[block].shape, (8,))
            self.assertTrue(np.all(np.diff(eigen.energies[block]) >= 0))
        self.assertIsNone(eigen.eigenvectors)
        self.assertFalse(eigen.with_vectors)
        self.assertEqual([p.name for p in eigen.sources],
                         ["spin-up_eigen.npz", "spin-down_eigen.npz"])
        assert_deep_frozen(eigen)

    def test_eigen_with_vectors(self):
        eigen = read_uhfr_eigen(self.hwave_run, with_vectors=True)
        self.assertTrue(eigen.with_vectors)
        for block in ("spin-up", "spin-down"):
            v = eigen.eigenvectors[block]
            self.assertEqual(v.shape, (8, 8))
            np.testing.assert_allclose(v.conj().T @ v, np.eye(8), atol=1e-10)
        assert_deep_frozen(eigen)

    def test_green(self):
        green = read_uhfr_green(self.toml)
        self.assertIsInstance(green, UHFrGreen)
        self.assertEqual(green.green.shape, (16, 16))
        self.assertEqual(green.pairs.shape, (16, 4))
        self.assertEqual(int(green.mask.sum()), 16)
        np.testing.assert_array_equal(
            green.pairs, [(i, s, i, s) for i in range(8) for s in (0, 1)])
        # the values are the file's, at index site + N * spin
        first = green_lines(green.source)[3].split()
        i, s = int(first[0]), int(first[1])
        self.assertEqual(green.green[i + 8 * s, i + 8 * s],
                         complex(float(first[4]), float(first[5])))
        self.assertFalse(green.mask[0, 8])
        self.assertEqual(green.green[0, 8], 0)
        assert_deep_frozen(green)

    def test_tutorial_fixed_2sz(self):
        op = read_uhfr_order_parameters(self.toml)
        energy = read_energy(self.toml)
        tol = sum_rule_bound(self.hwave_run, 8)
        self.assertLessEqual(abs(op.charge.sum() - energy.ncond), tol)
        self.assertLessEqual(abs(op.charge.sum() - self.hwave_run.param("Ncond")), tol)
        self.assertLessEqual(abs(8 * op.uniform[2] - energy.sz), tol)
        self.assertIs(op.transverse_available, False)
        np.testing.assert_array_equal(op.spin[:, :2], 0.0)
        self.assertEqual(op.n_cell, 1)
        self.assertEqual(op.source, read_uhfr_green(self.toml).source)
        assert_deep_frozen(op)

    def test_offdiagonal_not_required_fixed_2sz(self):
        green = read_uhfr_green(self.toml)
        self.assertFalse(any(s != t for _, s, _, t in green.pairs.tolist()))
        read_uhfr_order_parameters(self.toml)       # does not raise

    def test_energy_has_sx_sy(self):
        energy = read_energy(self.toml)
        self.assertEqual(energy.sx, 0.0)
        self.assertEqual(energy.sy, 0.0)
        self.assertEqual(set(energy.mu), {"spin-up", "spin-down"})
        u8 = read_energy(fixture("uhfr_u8"))
        self.assertIsNotNone(u8.sx)
        self.assertIsNotNone(u8.sy)
        self.assertEqual(list(u8.values)[-2:], ["Sx", "Sy"])


class TestTiltedNeel(unittest.TestCase):
    """tests/uhfr/coulombintra: the 4 x 4 Hubbard cluster at U = 8, a
    collinear Neel state with a tilted axis; every (i, s, i, t) entry."""

    @classmethod
    def setUpClass(cls):
        cls.toml = fixture("uhfr_u8")
        cls.hwave_run = open_run(cls.toml)
        cls.m, cls.e = neel_gap_equation(4, 8)

    def test_oracle_values(self):
        self.assertAlmostEqual(self.m, 0.4474190440247483, places=13)
        self.assertAlmostEqual(self.e, -7.38962582304832, places=10)
        self.assertLessEqual(abs(read_energy(self.toml).total - self.e), 1e-8)

    def test_trans_def_is_the_square_lattice(self):
        # the precondition of the positions: every hopping links sites one
        # lattice step apart (periodic modulo 4)
        lines = (self.toml.parent / "trans.def").read_text().splitlines()[5:]
        rows = [line.split() for line in lines if line.strip()]
        self.assertEqual(len(rows), 128)
        for row in rows:
            i, j = int(row[0]), int(row[2])
            d = (U8_POSITIONS[i] - U8_POSITIONS[j]) % 4
            d = np.where(d > 2, d - 4, d)
            self.assertEqual(np.abs(d).sum(), 1, row)

    def test_eigen_sz_free(self):
        eigen = read_uhfr_eigen(self.toml, with_vectors=True)
        self.assertEqual(list(eigen.energies), ["sz-free"])
        self.assertEqual(eigen.energies["sz-free"].shape, (32,))
        self.assertEqual(eigen.eigenvectors["sz-free"].shape, (32, 32))
        self.assertEqual(eigen.sources[0].name, "sz-free_eigen.dat.npz")
        assert_deep_frozen(eigen)

    def test_green(self):
        green = read_uhfr_green(self.toml)
        self.assertEqual(green.green.shape, (32, 32))
        self.assertEqual(int(green.mask.sum()), 64)
        g = np.where(green.mask, green.green, 0)
        np.testing.assert_allclose(g, g.conj().T, atol=1e-12)
        assert_deep_frozen(green)

    def test_uhfr_tilted_neel(self):
        op = read_uhfr_order_parameters(self.toml, site_positions=U8_POSITIONS,
                                        lattice_shape=U8_LATTICE)
        self.assertLessEqual(abs(op.staggered_norm - self.m), 1e-8)
        self.assertLessEqual(abs(op.mean_local_moment - self.m), 1e-8)
        self.assertEqual(op.cell_shape, U8_LATTICE)
        self.assertEqual(op.lattice_shape, U8_LATTICE)
        self.assertEqual(op.n_cell, 1)
        self.assertIs(op.transverse_available, True)
        # the axis is tilted: the transverse entries carry part of the moment
        self.assertGreater(float(np.max(np.abs(op.spin[:, :2]))), 1e-3)
        assert_deep_frozen(op)

    def test_totals(self):
        op = read_uhfr_order_parameters(self.toml)
        energy = read_energy(self.toml)
        tol = sum_rule_bound(self.hwave_run, 16)
        self.assertLessEqual(abs(op.charge.sum() - energy.ncond), tol)
        for axis, total in enumerate((energy.sx, energy.sy, energy.sz)):
            self.assertLessEqual(abs(16 * op.uniform[axis] - total), tol)

    def test_spin_orientation(self):
        # G[a, b] = <c^dag_a c_b>: a spin along +y on site 0 has
        # G[0, N] = i/2 and G[N, 0] = -i/2 (the order.py convention); the
        # file lines are "i s j t Re Im"
        toml = copy_fixture("uhfr_u8")
        path = toml.parent / "output" / "green.dat"
        values = {(0, 0, 0, 0): (0.5, 0.0), (0, 1, 0, 1): (0.5, 0.0),
                  (0, 0, 0, 1): (0.0, 0.5), (0, 1, 0, 0): (0.0, -0.5)}
        lines = []
        for line in green_lines(path):
            key = tuple(int(x) for x in line.split()[:4])
            if key in values:
                line = "{} {} {} {} {} {}".format(*key, *values[key])
            lines.append(line)
        path.write_text("\n".join(lines) + "\n")
        op = read_uhfr_order_parameters(toml)
        np.testing.assert_allclose(op.spin[0], (0.0, 0.5, 0.0), atol=1e-15)
        self.assertEqual(op.charge[0], 1.0)

    def test_without_lattice_shape(self):
        for kwargs in ({}, {"site_positions": U8_POSITIONS},
                       {"lattice_shape": U8_LATTICE}):
            op = read_uhfr_order_parameters(self.toml, **kwargs)
            self.assertIsNone(op.staggered, kwargs)
            self.assertIsNone(op.staggered_norm, kwargs)
            self.assertIsNone(op.spin_fourier, kwargs)
            self.assertIsNone(op.q_frac, kwargs)
            self.assertEqual(op.charge.shape, (16,))
            self.assertLessEqual(abs(op.mean_local_moment - self.m), 1e-8)
            expected = kwargs.get("lattice_shape")
            self.assertEqual(op.cell_shape, expected, kwargs)
            self.assertEqual(op.lattice_shape, expected, kwargs)
            self.assertEqual(op.n_cell, 1)
            self.assertIn("fourier_unavailable", op.meta)


class TestRequiredEntries(unittest.TestCase):
    """The OneBodyG entries the order parameters need, and the layout
    checks of green.dat."""

    def test_missing_diagonal_entry(self):
        # the tutorial's own greenone.def lists sites 0-4 only
        toml = copy_fixture("uhfr_tutorial")
        (toml.parent / "greenone.def").write_bytes(
            (REPO / "docs" / "tutorial" / "Hubbard" / "UHFr" / "greenone.def")
            .read_bytes())
        run_input(toml)
        read_uhfr_green(toml)                       # the file itself reads
        with self.assertRaises(MissingOutputError) as ctx:
            read_uhfr_order_parameters(toml)
        message = str(ctx.exception)
        for site in (5, 6, 7):
            for s in (0, 1):
                self.assertIn("({0} {1} {0} {1})".format(site, s), message)
        self.assertNotIn("(4 0 4 0)", message)
        self.assertIn("add these lines to the OneBodyG file 'greenone.def' and "
                      "rerun", message)

    def test_missing_offdiagonal_entry_sz_free(self):
        toml = copy_fixture("uhfr_u8")
        write_greenone(toml.parent / "greenone.def",
                       [(i, s, i, t) for i in range(16) for s in (0, 1)
                        for t in (0, 1) if (i, s, t) != (3, 1, 0)])
        run_input(toml)
        self.assertEqual(len(green_lines(toml.parent / "output" / "green.dat")), 63)
        with self.assertRaises(MissingOutputError) as ctx:
            read_uhfr_order_parameters(toml)
        message = str(ctx.exception)
        self.assertIn("(3 1 3 0)", message)
        self.assertEqual(len(re.findall(r"\(\d+ \d \d+ \d\)", message)), 1)

    def test_many_missing_listed_up_to_ten(self):
        # only the diagonal entries in a sz-free run: 32 are missing
        toml = copy_fixture("uhfr_u8")
        path = toml.parent / "output" / "green.dat"
        keep = [line for line in green_lines(path)
                if line.split()[1] == line.split()[3]]
        path.write_text("\n".join(keep) + "\n")
        with self.assertRaises(MissingOutputError) as ctx:
            read_uhfr_order_parameters(toml)
        message = str(ctx.exception)
        self.assertEqual(len(re.findall(r"\(\d+ \d \d+ \d\)", message)), 10)
        self.assertIn("and 22 more", message)
        self.assertIn("(0 0 0 1)", message)

    def test_onebodyg_key_in_any_case(self):
        # UHFr reads [file.input] OneBodyG in any case (qlms.py:155-159)
        toml = copy_fixture("uhfr_tutorial")
        text = toml.read_text()
        toml.write_text(text.replace('OneBodyG = "greenone.def"',
                                     'onebodyG = "greenone.def"'))
        path = toml.parent / "output" / "green.dat"
        path.write_text("\n".join(green_lines(path)[:-1]) + "\n")
        with self.assertRaises(MissingOutputError) as ctx:
            read_uhfr_order_parameters(toml)
        self.assertIn("(7 1 7 1)", str(ctx.exception))
        self.assertIn("'greenone.def'", str(ctx.exception))

    def test_duplicate_pair_rejected(self):
        toml = copy_fixture("uhfr_tutorial")
        path = toml.parent / "output" / "green.dat"
        lines = green_lines(path)
        path.write_text("\n".join(lines + [lines[0]]) + "\n")
        with self.assertRaises(UnsupportedLayoutError) as ctx:
            read_uhfr_green(toml)
        message = str(ctx.exception)
        self.assertIn("(0 0 0 0)", message)
        self.assertIn("lines 1 and 17", message)

    def test_bad_lines_rejected(self):
        toml = copy_fixture("uhfr_tutorial")
        path = toml.parent / "output" / "green.dat"
        lines = green_lines(path)
        for bad in ("0 0 0 0 1.0", "0 0 x 0 1.0 0.0", "8 0 8 0 1.0 0.0",
                    "0 2 0 2 1.0 0.0", "-1 0 0 0 1.0 0.0",
                    "0 0 0 0 nan 0.0", "0.5 0 0 0 1.0 0.0"):
            path.write_text("\n".join(lines[1:] + [bad]) + "\n")
            with self.assertRaises(UnsupportedLayoutError, msg=bad) as ctx:
                read_uhfr_green(toml)
            self.assertIn("line 16", str(ctx.exception), bad)

    def test_green_needs_onebodyg(self):
        # without OneBodyG the solver writes no green.dat (uhfr.py:1302-1314)
        toml = copy_fixture("uhfr_tutorial")
        toml.write_text(toml.read_text().replace('OneBodyG = "greenone.def"\n', ""))
        (toml.parent / "output" / "green.dat").unlink()
        with self.assertRaises(MissingOutputError) as ctx:
            read_uhfr_green(toml)
        self.assertIn("OneBodyG", ctx.exception.hint)

    def test_zero_steps(self):
        toml = fixture("uhfr_tutorial_zero_steps")
        self.assertEqual(read_scf_status(toml).iterations, 0)
        for reader in (read_uhfr_green, read_uhfr_order_parameters):
            with self.assertRaises(MissingOutputError) as ctx:
                reader(toml)
            self.assertIn("no SCF step", ctx.exception.hint)
        with self.assertRaises(MissingOutputError):
            read_uhfr_eigen(toml)

    def test_eigen_not_configured(self):
        toml = copy_fixture("uhfr_tutorial")
        toml.write_text(toml.read_text().replace('eigen = "eigen"\n', ""))
        with self.assertRaises(MissingOutputError):
            read_uhfr_eigen(toml)

    def test_eigen_bad_shape(self):
        toml = copy_fixture("uhfr_tutorial")
        path = toml.parent / "output" / "spin-up_eigen.npz"
        np.savez(path, eigenvalue=np.zeros(7), eigenvector=np.eye(7))
        with self.assertRaises(UnsupportedLayoutError):
            read_uhfr_eigen(toml)

    def test_not_uhfr(self):
        with self.assertRaises(UnsupportedLayoutError):
            read_uhfr_green(fixture("uhfk_cdw"))


if __name__ == "__main__":
    unittest.main()
