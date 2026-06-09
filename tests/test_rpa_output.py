#!/usr/bin/env python3
"""RPA output: the transverse susceptibility chiq_pm must be written.

For calc_type == "ring+ladder" the solver computes the transverse channel
chi_+-(q) and stores it in green_info["chiq_pm"]; save_results must persist it
so the (otherwise wasted) computation is accessible.
"""

import os
import tempfile
import unittest

import numpy as np

from hwave.solver.rpa import RPA


class TestRPASaveChiqPm(unittest.TestCase):
    def test_chiq_pm_is_written_to_chiq_npz(self):
        solver = RPA.__new__(RPA)
        solver._init_wavevec = lambda: None  # bypass heavy wavevector setup
        solver.calc_chiq = True
        solver.freq_index = np.array([0])
        solver.kvec = np.zeros((1, 3))
        solver.wavenum_table = np.zeros((1, 3), dtype=int)

        green_info = {
            "chiq": np.ones((1, 1, 1, 1), dtype=complex),
            "chiq_pm": np.full((1, 1, 1, 1), 2.0, dtype=complex),
        }

        tmpdir = tempfile.mkdtemp()
        try:
            info_outputfile = {"path_to_output": tmpdir, "chiq": "chiq"}
            solver.save_results(info_outputfile, green_info)

            data = np.load(os.path.join(tmpdir, "chiq.npz"))
            self.assertIn("chiq_pm", data.files)
            self.assertAlmostEqual(data["chiq_pm"].flat[0].real, 2.0, places=10)
            # spin-orbital axis ordering must be recorded for downstream consumers
            self.assertIn("index_convention", data.files)
            self.assertEqual(str(data["index_convention"]), "spin_block")
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
