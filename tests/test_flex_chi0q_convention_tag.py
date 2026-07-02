"""Regression test: FLEX chi0q.npz must carry the index_convention tag.

FLEX chi0q is produced by the same spin-block RPA internals as RPA chi0q,
and the chi0q consumers (RPA.read_chi0q, hwave_sc._load_chi0q) reject
untagged files in spin-orbital mode as "predating the convention fix".
FLEX.save_results must therefore write index_convention='spin_block' just
like RPA.save_results does.
"""
import os
import tempfile
import unittest

import numpy as np

from hwave.solver.flex import FLEX
from hwave.solver.rpa import validate_chi0q_index_convention


def _make_flex_stub(nmat, freq_index):
    flex = FLEX.__new__(FLEX)
    flex._init_wavevec = lambda: None
    flex.nmat = nmat
    flex.freq_index = freq_index
    flex.kvec = np.eye(3)
    flex.wavenum_table = np.zeros((2, 3))
    flex._flex_general = False
    return flex


class TestFlexChi0qConventionTag(unittest.TestCase):
    def test_chi0q_npz_carries_spin_block_tag(self):
        flex = _make_flex_stub(4, np.arange(4))

        with tempfile.TemporaryDirectory() as tmp:
            flex.save_results(
                {"path_to_output": tmp, "chi0q": "chi0q"},
                {"chi0q": np.zeros((4, 2, 2, 2), dtype=np.complex128)})
            data = np.load(os.path.join(tmp, "chi0q.npz"))
            self.assertIn(
                "index_convention", data,
                "FLEX chi0q.npz must carry the index_convention tag")
            self.assertEqual(str(data["index_convention"]), "spin_block")
            # the SO-mode validator must accept a FLEX-written file
            validate_chi0q_index_convention(
                data, enable_spin_orbital=True, file_name="chi0q.npz")

    def test_chi0q_npz_freq_metadata_matches_full_grid(self):
        # FLEX never applies the matsubara_frequency filter to its outputs
        # (unlike RPA), so the saved freq_index must describe the FULL grid
        # actually stored, not the restricted user option; otherwise the
        # strict hwave_sc loader rejects a perfectly valid file.
        flex = _make_flex_stub(4, np.array([2]))  # e.g. "center"

        with tempfile.TemporaryDirectory() as tmp:
            flex.save_results(
                {"path_to_output": tmp, "chi0q": "chi0q"},
                {"chi0q": np.zeros((4, 2, 2, 2), dtype=np.complex128)})
            data = np.load(os.path.join(tmp, "chi0q.npz"))
            np.testing.assert_array_equal(data["freq_index"], np.arange(4))
            self.assertEqual(int(data["nmat"]), 4)


if __name__ == "__main__":
    unittest.main()
