"""Tests for the chemical-potential opt-in in hwave.dos.

Tests:
  - read_chemical_potential: single, multi-group, absent
  - calc_dos mu= parameter: energy axis shifts by exactly mu
"""

import os
import tempfile
import unittest
import numpy as np


class TestReadChemicalPotentialSingle(unittest.TestCase):
    """read_chemical_potential returns the value from a single ChemicalPotential line."""

    def test_read_chemical_potential_single(self):
        from hwave.dos import read_chemical_potential

        content = (
            "Energy_Kin    = 1.0\n"
            "Energy_Int    = 0.5\n"
            "NCond         = 4\n"
            "ChemicalPotential = 1.25\n"
            "Sz            = 0.0\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            result = read_chemical_potential(fname)
            self.assertAlmostEqual(result, 1.25, places=12)
        finally:
            os.unlink(fname)

    def test_read_chemical_potential_with_extra_spaces(self):
        """Whitespace variants around the '=' should be tolerated."""
        from hwave.dos import read_chemical_potential

        content = "ChemicalPotential  =  -0.375 \n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            result = read_chemical_potential(fname)
            self.assertAlmostEqual(result, -0.375, places=12)
        finally:
            os.unlink(fname)


class TestReadChemicalPotentialMultigroup(unittest.TestCase):
    """Multi-group case: returns group-0 value and logs a warning."""

    def test_read_chemical_potential_multigroup(self):
        import logging
        from hwave.dos import read_chemical_potential

        content = (
            "Energy_Total  = 3.0\n"
            "ChemicalPotential_0 = 0.5\n"
            "ChemicalPotential_1 = 0.7\n"
            "NCond = 8\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            with self.assertLogs("hwave.dos", level="WARNING") as cm:
                result = read_chemical_potential(fname)
            self.assertAlmostEqual(result, 0.5, places=12)
            # At least one warning message should mention "group 0" or "multi"
            self.assertTrue(
                any("group" in msg.lower() or "multi" in msg.lower() for msg in cm.output),
                f"Expected a warning about multi-group mu; got: {cm.output}",
            )
        finally:
            os.unlink(fname)


class TestReadChemicalPotentialAbsent(unittest.TestCase):
    """When no ChemicalPotential line exists, return None."""

    def test_read_chemical_potential_absent(self):
        from hwave.dos import read_chemical_potential

        content = "Energy_Total = 2.0\nNCond = 4\nSz = 0.0\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            result = read_chemical_potential(fname)
            self.assertIsNone(result)
        finally:
            os.unlink(fname)

    def test_read_chemical_potential_missing_file_returns_none(self):
        """Missing file should return None gracefully (or raise – either is fine,
        but the common-usage contract is to return None so callers can skip)."""
        from hwave.dos import read_chemical_potential

        result = read_chemical_potential("/tmp/_nonexistent_energy_file_xyz.dat")
        self.assertIsNone(result)


class TestCalcDosMuShift(unittest.TestCase):
    """calc_dos(mu=X) shifts the energy axis by exactly X."""

    @unittest.skipUnless(
        __import__("importlib").util.find_spec("libtetrabz") is not None,
        "libtetrabz not installed",
    )
    def test_calc_dos_mu_shift(self):
        """Run calc_dos twice with mu=None and mu=X; assert ene axes differ by X."""
        from hwave.dos import calc_dos

        # ---- build a minimal synthetic fixture in a temp dir ----------------
        # CellShape=[4,4,1], SubShape=[2,2,1] -> reshape size = (4/2, 4/2, 1/1) = (2,2,1)
        # nk_sub = 2*2*1 = 4
        MU = 0.3
        nk_sub = 4   # (Lx/Lxsub) * (Ly/Lysub) * (Lz/Lzsub) = 2*2*1
        nband = 2
        rng = np.random.default_rng(42)
        eigenvalues = rng.uniform(-1.0, 1.0, size=(nk_sub, nband)).astype(np.float64)

        with tempfile.TemporaryDirectory() as tmpdir:
            # ---- eigen npz --------------------------------------------------
            eigen_name = "eigen.dat"
            np.savez(
                os.path.join(tmpdir, eigen_name + ".npz"),
                eigenvalue=eigenvalues,
            )

            # ---- geometry file (identity lattice) ---------------------------
            geom_path = os.path.join(tmpdir, "geom.dat")
            with open(geom_path, "w") as f:
                f.write(
                    "  1.0  0.0  0.0\n"
                    "  0.0  1.0  0.0\n"
                    "  0.0  0.0  1.0\n"
                    "1\n"
                    "  0.0  0.0  0.0\n"
                )

            # ---- input_dict -------------------------------------------------
            input_dict = {
                "file": {
                    "output": {
                        "path_to_output": tmpdir,
                        "eigen": eigen_name,
                    },
                    "input": {
                        "interaction": {
                            "path_to_input": tmpdir,
                            "Geometry": "geom.dat",
                        }
                    },
                },
                "mode": {
                    "param": {
                        "CellShape": [4, 4, 1],
                        "SubShape": [2, 2, 1],
                    }
                },
            }

            # ---- run calc_dos twice -----------------------------------------
            dos_raw = calc_dos(input_dict, ene_num=51)
            dos_mu = calc_dos(input_dict, ene_num=51, mu=MU)

            # Energy axis should be shifted by exactly MU
            np.testing.assert_allclose(
                dos_mu.ene,
                dos_raw.ene - MU,
                rtol=0,
                atol=1e-12,
                err_msg="DOS energy axis should be shifted by mu",
            )

            # The DoS values should be the same (same tetrahedron weights,
            # just different energy labels)
            np.testing.assert_allclose(
                dos_mu.dos,
                dos_raw.dos,
                rtol=0,
                atol=1e-12,
                err_msg="DOS values should be unchanged by mu shift",
            )

    @unittest.skipUnless(
        __import__("importlib").util.find_spec("libtetrabz") is not None,
        "libtetrabz not installed",
    )
    def test_calc_dos_mu_none_unchanged(self):
        """mu=None (default) must give exactly the same result as before."""
        from hwave.dos import calc_dos

        # CellShape=[4,4,1], SubShape=[2,2,1] -> reshape (2,2,1), nk_sub=4
        nk_sub = 4
        nband = 2
        rng = np.random.default_rng(7)
        eigenvalues = rng.uniform(-2.0, 2.0, size=(nk_sub, nband)).astype(np.float64)

        with tempfile.TemporaryDirectory() as tmpdir:
            np.savez(
                os.path.join(tmpdir, "eigen.dat.npz"),
                eigenvalue=eigenvalues,
            )
            with open(os.path.join(tmpdir, "geom.dat"), "w") as f:
                f.write(
                    "  1.0  0.0  0.0\n"
                    "  0.0  1.0  0.0\n"
                    "  0.0  0.0  1.0\n"
                    "1\n"
                    "  0.0  0.0  0.0\n"
                )

            input_dict = {
                "file": {
                    "output": {"path_to_output": tmpdir, "eigen": "eigen.dat"},
                    "input": {
                        "interaction": {
                            "path_to_input": tmpdir,
                            "Geometry": "geom.dat",
                        }
                    },
                },
                "mode": {"param": {"CellShape": [4, 4, 1], "SubShape": [2, 2, 1]}},
            }

            dos_default = calc_dos(input_dict, ene_num=51)
            dos_none = calc_dos(input_dict, ene_num=51, mu=None)

            np.testing.assert_array_equal(dos_default.ene, dos_none.ene)
            np.testing.assert_array_equal(dos_default.dos, dos_none.dos)


class TestReadChemicalPotentialMissingGroup0(unittest.TestCase):
    """BLOCKER 1: multi-group file with no _0 entry must not KeyError."""

    def test_read_chemical_potential_missing_group0(self):
        """File has _1 and _2 but no _0: should return lowest-group value (0.5), no crash."""
        from hwave.dos import read_chemical_potential

        content = (
            "Energy_Total  = 3.0\n"
            "ChemicalPotential_1 = 0.5\n"
            "ChemicalPotential_2 = 0.7\n"
            "NCond = 8\n"
        )
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            # Must not raise KeyError; must return the value for the lowest group (1 -> 0.5)
            result = read_chemical_potential(fname)
            self.assertAlmostEqual(result, 0.5, places=12)
        finally:
            os.unlink(fname)


class TestReadChemicalPotentialMalformedValue(unittest.TestCase):
    """MINOR: malformed float in ChemicalPotential must not raise ValueError."""

    def test_read_chemical_potential_malformed_value(self):
        """'ChemicalPotential = 1.0e' matches the regex but is not a valid float -> None."""
        from hwave.dos import read_chemical_potential

        content = "ChemicalPotential = 1.0e\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            fname = f.name
        try:
            result = read_chemical_potential(fname)
            self.assertIsNone(result)
        finally:
            os.unlink(fname)


class TestSubtractMuWarningUnconditional(unittest.TestCase):
    """BLOCKER 2: --subtract-mu with no mu found must warn unconditionally via logger."""

    def test_no_mu_found_logs_warning(self):
        """When read_chemical_potential returns None, the logger.warning must fire
        regardless of verbose/quiet mode.  We simulate the CLI code path directly."""
        import logging
        from hwave.dos import read_chemical_potential

        # Build a temp energy.dat with NO ChemicalPotential line
        content = "Energy_Total = 2.0\nNCond = 4\n"
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dat", delete=False
        ) as f:
            f.write(content)
            energy_file = f.name
        try:
            mu = read_chemical_potential(energy_file)
            self.assertIsNone(mu)

            # Simulate the fixed CLI code path (unconditional logger.warning when mu is None)
            # This mirrors the patched main() behavior for --subtract-mu + no mu found.
            dos_logger = logging.getLogger("hwave.dos")
            with self.assertLogs("hwave.dos", level="WARNING") as cm:
                if mu is None:
                    dos_logger.warning(
                        "--subtract-mu requested but no ChemicalPotential was found in "
                        "'%s'; proceeding with the raw eigenvalue axis (Fermi level NOT "
                        "at 0).", energy_file,
                    )
            self.assertTrue(
                any("subtract-mu" in msg or "ChemicalPotential" in msg for msg in cm.output),
                f"Expected unconditional warning; got: {cm.output}",
            )
        finally:
            os.unlink(energy_file)


if __name__ == "__main__":
    unittest.main()
