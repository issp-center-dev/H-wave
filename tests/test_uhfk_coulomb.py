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
            inter[explicit_key] = coulomb_file  # ambiguous co-specification
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


if __name__ == "__main__":
    unittest.main()
