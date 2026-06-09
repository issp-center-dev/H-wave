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


if __name__ == "__main__":
    unittest.main()
