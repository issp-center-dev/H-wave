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
        solver.nmat = 1
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


class TestRPAChi0qPassthroughMetadata(unittest.TestCase):
    """A chi0q loaded via chi0q_init and re-saved must keep the frequency
    metadata of the file that PRODUCED it: stamping the current run's
    freq_index/nmat mislabels the axis for the strict hwave_sc loader."""

    def _make_solver(self):
        solver = RPA.__new__(RPA)
        solver._init_wavevec = lambda: None
        solver.calc_chiq = False
        solver.nmat = 16
        solver.freq_index = np.arange(16)
        solver.kvec = np.zeros((1, 3))
        solver.wavenum_table = np.zeros((1, 3), dtype=int)
        return solver

    def _save_and_load(self, solver, tmpdir):
        green_info = {"chi0q": np.ones((5, 1, 1, 1), dtype=complex)}
        solver.save_results({"path_to_output": tmpdir, "chi0q": "chi0q"},
                            green_info)
        return np.load(os.path.join(tmpdir, "chi0q.npz"))

    def test_passthrough_keeps_producing_run_metadata(self):
        solver = self._make_solver()
        solver._chi0q_init_meta = {"freq_index": np.arange(2, 7), "nmat": 8}
        tmpdir = tempfile.mkdtemp()
        try:
            data = self._save_and_load(solver, tmpdir)
            np.testing.assert_array_equal(data["freq_index"], np.arange(2, 7))
            self.assertEqual(int(data["nmat"]), 8)
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_passthrough_of_legacy_file_describes_axis_without_nmat(self):
        # the documented file format guarantees a freq_index key, so a
        # legacy passthrough must still write one -- describing the stored
        # axis (0..n-1) without fabricating an nmat claim (downstream then
        # treats it explicitly as ambiguous unless Nmat matches)
        solver = self._make_solver()
        solver._chi0q_init_meta = {"freq_index": None, "nmat": None}
        tmpdir = tempfile.mkdtemp()
        try:
            data = self._save_and_load(solver, tmpdir)
            np.testing.assert_array_equal(data["freq_index"], np.arange(5))
            self.assertNotIn("nmat", data.files)
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_own_chi0q_uses_current_run_metadata(self):
        solver = self._make_solver()
        solver._chi0q_init_meta = None  # chi0q computed by this run
        tmpdir = tempfile.mkdtemp()
        try:
            data = self._save_and_load(solver, tmpdir)
            np.testing.assert_array_equal(data["freq_index"], np.arange(16))
            self.assertEqual(int(data["nmat"]), 16)
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_chiq_from_chi0q_init_keeps_producing_run_metadata(self):
        # chiq is computed FROM the chi0q input, so its frequency axis is
        # inherited from the chi0q_init file, not the current run
        solver = self._make_solver()
        solver.calc_chiq = True
        solver._chi0q_init_meta = {"freq_index": np.arange(2, 7), "nmat": 8}
        tmpdir = tempfile.mkdtemp()
        try:
            green_info = {"chiq": np.ones((5, 1, 1, 1), dtype=complex)}
            solver.save_results({"path_to_output": tmpdir, "chiq": "chiq"},
                                green_info)
            data = np.load(os.path.join(tmpdir, "chiq.npz"))
            np.testing.assert_array_equal(data["freq_index"], np.arange(2, 7))
            self.assertEqual(int(data["nmat"]), 8)
        finally:
            import shutil
            shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
