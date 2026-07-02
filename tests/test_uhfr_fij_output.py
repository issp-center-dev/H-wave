"""Regression test: UHFr fij/eigen output after sub-block detection.

save_results looks up the fixed green_list keys ('sz-free', 'spin-up',
'spin-down').  Block detection must not rename these keys even when a block
splits into disconnected sub-blocks, otherwise a standard mVMC workflow
(Trans-only system, fij output) dies with KeyError after the SCF loop and
the per-key eigen npz files change names.
"""
import os
import tempfile
import unittest

import numpy as np

import hwave.qlms


def _write_trans(filename):
    # Two decoupled sites (diagonal on-site energies): the sz-free block
    # splits into 4 disconnected sub-blocks.
    rows = []
    for s in range(2):
        rows.append((0, s, 0, s, 1.0))    # on-site -1
        rows.append((1, s, 1, s, -1.0))   # on-site +1
    with open(filename, "w") as f:
        f.write("========================\n")
        f.write("NTransfer {}\n".format(len(rows)))
        f.write("========================\n")
        f.write("========i_j_s_tijs======\n")
        f.write("========================\n")
        for i, s, j, t, v in rows:
            f.write("{} {} {} {} {:.12f} 0.0\n".format(i, s, j, t, v))


class TestUHFrFijOutput(unittest.TestCase):
    def test_fij_and_eigen_written_with_standard_keys(self):
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
                            "Nsite": 2,
                            "Ncond": 2,
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
                            "eigen": "eigen.dat",
                            "fij": "fij.dat",
                        },
                    },
                }
                hwave.qlms.run(input_dict=params)
                self.assertTrue(
                    os.path.exists("output/sz-free_fij.dat"),
                    "fij output must be written under the 'sz-free' key")
                self.assertTrue(
                    os.path.exists("output/sz-free_eigen.dat.npz"),
                    "eigen npz must keep the 'sz-free' key name")
                # eigenvalues must cover the whole block, sorted ascending
                data = np.load("output/sz-free_eigen.dat.npz")
                ev = data["eigenvalue"]
                self.assertEqual(ev.shape, (4,))
                np.testing.assert_array_equal(ev, np.sort(ev))
            finally:
                os.chdir(cur)


if __name__ == "__main__":
    unittest.main()
