"""Regression test: UHFr block splitting must fill states by energy.

When block detection splits a spin sector into disconnected sub-blocks, the
electrons of that sector must still occupy the lowest eigenvalues of the
whole sector (a single Fermi level per sector), not a per-sub-block quota.

System: 3 sites, no interaction.
  site 0: on-site energy +10 (isolated)
  sites 1,2: on-site energy -10, hopping matrix element -1 (levels -11, -9)
Ncond = 2, 2Sz = 2, T = 0: both electrons are spin-up and must occupy the
levels -11 and -9, giving total energy -20.  A size-proportional split
instead pins one electron on the +10 site (energy -1).
"""
import os
import tempfile
import unittest

import numpy as np

import hwave.qlms


def _write_trans(filename):
    # H = - sum t_{isjt} c^dag_{is} c_{jt}  (Transfer coeff = -1)
    rows = []
    for s in range(2):
        rows.append((0, s, 0, s, -10.0))   # on-site +10
        rows.append((1, s, 1, s, +10.0))   # on-site -10
        rows.append((2, s, 2, s, +10.0))   # on-site -10
        rows.append((1, s, 2, s, +1.0))    # hop -1
        rows.append((2, s, 1, s, +1.0))
    with open(filename, "w") as f:
        f.write("========================\n")
        f.write("NTransfer {}\n".format(len(rows)))
        f.write("========================\n")
        f.write("========i_j_s_tijs======\n")
        f.write("========================\n")
        for i, s, j, t, v in rows:
            f.write("{} {} {} {} {:.12f} 0.0\n".format(i, s, j, t, v))


def _read_energy(filename):
    tbl = {}
    with open(filename, "r") as f:
        for line in f.read().splitlines():
            if "=" in line:
                k, v = line.split("=")
                tbl[k.strip().lower()] = float(v)
    return tbl


class TestUHFr2SzSpinFlipRejected(unittest.TestCase):
    def test_spin_flip_transfer_with_fixed_2sz_raises(self):
        """A one-body spin-flip term breaks Sz conservation, so a fixed 2Sz
        is ill-defined.  UHFk rejects it; UHFr must too (previously the
        cross-spin Hamiltonian entries were silently zeroed every SCF step,
        so the two solvers disagreed on the same input, one of them
        silently)."""
        cur = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                rows = []
                for s in range(2):
                    rows.append((0, s, 1, s, 1.0))
                    rows.append((1, s, 0, s, 1.0))
                # spin-flip on site 0: up <-> down
                rows.append((0, 0, 0, 1, 0.3))
                rows.append((0, 1, 0, 0, 0.3))
                with open("trans.def", "w") as f:
                    f.write("========================\n")
                    f.write("NTransfer {}\n".format(len(rows)))
                    f.write("========================\n")
                    f.write("========i_j_s_tijs======\n")
                    f.write("========================\n")
                    for i, s, j, t, v in rows:
                        f.write("{} {} {} {} {:.12f} 0.0\n".format(
                            i, s, j, t, v))
                params = {
                    "log": {"print_level": 1, "print_step": 1},
                    "mode": {
                        "mode": "UHFr",
                        "param": {
                            "Nsite": 2,
                            "Ncond": 2,
                            "2Sz": 0,
                            "IterationMax": 100,
                            "EPS": 8,
                            "Mix": 0.5,
                            "RndSeed": 123456789,
                            "T": 0.0,
                        },
                    },
                    "file": {
                        "input": {
                            "path_to_input": "",
                            "interaction": {"Trans": "trans.def"},
                        },
                        "output": {
                            "path_to_output": "output",
                            "energy": "energy.dat",
                        },
                    },
                }
                with self.assertRaises(ValueError):
                    hwave.qlms.run(input_dict=params)
            finally:
                os.chdir(cur)


class TestUHFr2SzNoiseTolerated(unittest.TestCase):
    def test_noise_level_spin_flip_does_not_abort(self):
        """Noise-level cross-spin entries (relative ~1e-11 of the hopping
        scale) must be masked, not fatal: the guard compares against the
        transfer scale."""
        cur = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                rows = []
                for s in range(2):
                    rows.append((0, s, 1, s, 1.0))
                    rows.append((1, s, 0, s, 1.0))
                # numerical-noise spin flip, 11 orders below the hopping
                rows.append((0, 0, 0, 1, 1.0e-11))
                rows.append((0, 1, 0, 0, 1.0e-11))
                with open("trans.def", "w") as f:
                    f.write("========================\n")
                    f.write("NTransfer {}\n".format(len(rows)))
                    f.write("========================\n")
                    f.write("========i_j_s_tijs======\n")
                    f.write("========================\n")
                    for i, s, j, t, v in rows:
                        f.write("{} {} {} {} {:.14g} 0.0\n".format(
                            i, s, j, t, v))
                params = {
                    "log": {"print_level": 1, "print_step": 1},
                    "mode": {
                        "mode": "UHFr",
                        "param": {
                            "Nsite": 2,
                            "Ncond": 2,
                            "2Sz": 0,
                            "IterationMax": 100,
                            "EPS": 8,
                            "Mix": 0.5,
                            "RndSeed": 123456789,
                            "T": 0.0,
                        },
                    },
                    "file": {
                        "input": {
                            "path_to_input": "",
                            "interaction": {"Trans": "trans.def"},
                        },
                        "output": {
                            "path_to_output": "output",
                            "energy": "energy.dat",
                        },
                    },
                }
                hwave.qlms.run(input_dict=params)  # must NOT raise
                self.assertTrue(os.path.exists("output/energy.dat"))
            finally:
                os.chdir(cur)


class TestUHFrBlockOccupation(unittest.TestCase):
    def test_sector_electrons_fill_lowest_levels_across_sub_blocks(self):
        cur = os.getcwd()
        with tempfile.TemporaryDirectory() as tmp:
            os.chdir(tmp)
            try:
                _write_trans("trans.def")
                params = {
                    "log": {"print_level": 1, "print_step": 1},
                    "mode": {
                        "mode": "UHFr",
                        "param": {
                            "Nsite": 3,
                            "Ncond": 2,
                            "2Sz": 2,
                            "IterationMax": 100,
                            "EPS": 8,
                            "Mix": 0.5,
                            "RndSeed": 123456789,
                            "T": 0.0,
                        },
                    },
                    "file": {
                        "input": {
                            "path_to_input": "",
                            "interaction": {"Trans": "trans.def"},
                        },
                        "output": {
                            "path_to_output": "output",
                            "energy": "energy.dat",
                        },
                    },
                }
                hwave.qlms.run(input_dict=params)
                result = _read_energy("output/energy.dat")
                self.assertTrue(
                    np.isclose(result["energy_total"], -20.0,
                               rtol=0.0, atol=1.0e-8),
                    "expected ground-state energy -20.0 (both spin-up "
                    "electrons in the low-energy sub-block), got {}".format(
                        result["energy_total"]))
            finally:
                os.chdir(cur)


if __name__ == "__main__":
    unittest.main()
