"""Spin-orbital Eliashberg (hwave_sc) pairing vertex tests (#83)."""
import shutil
import tempfile
import unittest

import numpy as np

import tests.test_flex_spin_orbital as fso


class TestSharedVertex(unittest.TestCase):
    """The antisymmetrized spinful vertex has ONE builder, shared by the
    spin-orbital FLEX and the Eliashberg solver."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="so_vertex_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_matches_flex_spinful_vertex(self):
        from hwave.solver import rpa
        s, _ = fso._build(self.d, fso._write_model(self.d, True, soc=0.2))
        nvol, nd = s.lattice.nvol, 2 * s.norb
        ham = np.asarray(s.ham_info.ham_inter_q).reshape(nvol, nd, nd, nd, nd)
        got = rpa.spinful_bubble_vertex(ham, s.ham_info.spinful_exchange())
        want = np.asarray(s._spinful_vertex(np))
        self.assertEqual(got.shape, (nvol, nd * nd, nd * nd))
        np.testing.assert_array_equal(got, want)


if __name__ == "__main__":
    unittest.main()
