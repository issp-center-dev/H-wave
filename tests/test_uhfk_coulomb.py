#!/usr/bin/env python3
"""The combined 'Coulomb' interaction (RESPACK zvo_ur.dat style) for UHFk.

A single 'Coulomb' input is split by the solver into an on-site same-orbital
part (CoulombIntra) and the remainder (CoulombInter). Feeding the existing
CoulombInter data (all off-site, so entirely the inter part) through the
'Coulomb' keyword must reproduce the dedicated CoulombInter result.
"""

import os
import unittest

import numpy as np
import tomli

import hwave.qlms


def _read_energy_total(filename):
    with open(filename) as fh:
        for line in fh.read().splitlines():
            k, v = line.split("=")
            if k.strip().lower() == "energy_total":
                return float(v)
    raise KeyError("Energy_Total not found")


class TestUHFkCombinedCoulomb(unittest.TestCase):
    def test_coulomb_keyword_matches_coulombinter(self):
        cur = os.getcwd()
        os.chdir("tests/uhfk/CoulombInter")
        try:
            with open("input.toml", "rb") as f:
                params = tomli.load(f)

            inter = params["file"]["input"]["interaction"]
            coulomb_file = inter.pop("CoulombInter")
            inter["Coulomb"] = coulomb_file
            params["file"]["output"]["energy"] = "energy_coulomb.dat"

            hwave.qlms.run(input_dict=params)

            result = _read_energy_total("output/energy_coulomb.dat")
            expect = _read_energy_total("output_ref/energy.dat")
        finally:
            os.chdir(cur)

        self.assertTrue(np.isclose(result, expect, rtol=0.0, atol=1.0e-8))

    def _run_with_conflicting_coulomb(self, explicit_key):
        cur = os.getcwd()
        os.chdir("tests/uhfk/CoulombInter")
        try:
            with open("input.toml", "rb") as f:
                params = tomli.load(f)

            inter = params["file"]["input"]["interaction"]
            coulomb_file = inter.pop("CoulombInter")
            inter["Coulomb"] = coulomb_file
            if explicit_key == "CoulombIntra":
                # the CoulombInter file contains off-site entries, which
                # the read-time validation (issue #93) would reject under
                # the CoulombIntra rule BEFORE the conflict check this
                # test pins -- point the explicit key at a valid on-site
                # same-orbital file instead
                import tempfile
                tmp = tempfile.mkdtemp()
                self.addCleanup(__import__("shutil").rmtree, tmp, True)
                cpath = os.path.join(tmp, "coulombintra_conflict.dat")
                with open(cpath, "w") as fobj:
                    fobj.write("CoulombIntra for the conflict test\n"
                               "2\n1\n 1\n"
                               "   0    0    0    1    1   1.0   0.0\n"
                               "   0    0    0    2    2   1.0   0.0\n")
                inter[explicit_key] = cpath
            else:
                inter[explicit_key] = coulomb_file  # ambiguous co-spec
            params["file"]["output"]["energy"] = "energy_conflict.dat"

            with self.assertRaises(SystemExit):
                hwave.qlms.run(input_dict=params)
        finally:
            os.chdir(cur)

    def test_coulomb_combined_with_coulombinter_is_rejected(self):
        # 'Coulomb' together with explicit CoulombInter is ambiguous.
        self._run_with_conflicting_coulomb("CoulombInter")

    def test_coulomb_combined_with_coulombintra_is_rejected(self):
        # 'Coulomb' together with explicit CoulombIntra is ambiguous.
        self._run_with_conflicting_coulomb("CoulombIntra")




class TestUHFkIsingFactor(unittest.TestCase):
    """The historical /4 read an Ising file as J S^z S^z; the documented
    Hamiltonian is J (n_up - n_down)(n_up - n_down) (issue #106, decided
    with the maintainer). Equivalence pin: the corrected reader with the
    COUPLING SCALED BY 1/4 must reproduce the historical reference
    energies exactly -- proving the change is precisely the removal of
    the factor, not a different physics edit."""

    OLD_REFERENCE = {
        "Energy_Total": -10.74418063479572,
        "Energy_Band": -8.794240820599656,
        "Energy_CoulombIntra": -5.358774927711394,
        "Energy_Ising": 3.40883511351533,
    }

    def test_quarter_scaled_file_reproduces_the_old_reference(self):
        import shutil
        import tempfile
        import tomli
        import hwave.qlms

        src = "tests/uhfk/Ising"
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        for f in ("geom.dat", "transfer.dat", "coulombintra.dat",
                  "input.toml"):
            shutil.copy(os.path.join(src, f), d)
        # input.toml references files under output_ref (initial green)
        shutil.copytree(os.path.join(src, "output_ref"),
                        os.path.join(d, "output_ref"))
        with open(os.path.join(src, "ising.dat")) as f:
            lines = f.read().splitlines()
        out = lines[:4]
        for line in lines[4:]:
            p = line.split()
            out.append("  {} {} {} {} {}  {:.12f}  {:.12f}".format(
                *p[:5], float(p[5]) / 4.0, float(p[6]) / 4.0))
        with open(os.path.join(d, "ising.dat"), "w") as f:
            f.write("\n".join(out) + "\n")
        cur = os.getcwd()
        os.chdir(d)
        try:
            with open("input.toml", "rb") as f:
                params = tomli.load(f)
            hwave.qlms.run(input_dict=params)
            with open(os.path.join(
                    params["file"]["output"]["path_to_output"],
                    "energy.dat")) as f:
                text = f.read()
        finally:
            os.chdir(cur)
        got = {}
        for line in text.splitlines():
            if "=" in line:
                k, v = line.split("=")
                got[k.strip()] = float(v)
        for k, want in self.OLD_REFERENCE.items():
            # SCF-convergence tolerance, not bitwise: the seeded initial
            # green in output_ref was regenerated with the corrected
            # factor, so the iteration path differs at the EPS level
            self.assertAlmostEqual(got[k], want, places=7, msg=k)


if __name__ == "__main__":
    unittest.main()
