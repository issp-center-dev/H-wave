"""Regression test: wrap-around hopping must accumulate in H(k) scatter.

On a small periodic lattice (e.g. L=2 in x), the hops R=+1 and R=-1 are
physically distinct bonds that map onto the same array slot of the
real-space transfer table (index -1 wraps to L-1).  The scatter
``tab_r[R] = v`` silently overwrites one bond with the other, halving the
band width: for a 1-orbital ring with t=1 the true H(k) = 2 t cos(k) gives
bands {-2, +2} on the L=2 grid, not {-1, +1}.

Both UHFk and RPA build H(k) this way; they must agree with the analytic
result (and hence with each other).
"""
import os
import shutil
import tempfile
import unittest

import numpy as np


def _write_l2_chain(d):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("  1.0   0.0   0.0\n  0.0   1.0   0.0\n  0.0   0.0   1.0\n")
        f.write("1\n   0.0   0.0   0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer L2\n1\n2\n 1 1\n")
        f.write("   1    0    0    1    1  1.000000  0.0\n")
        f.write("  -1    0    0    1    1  1.000000  0.0\n")


class TestWraparoundHoppingUHFk(unittest.TestCase):
    def test_l2_ring_band_is_pm_2t(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.uhfk as uhfk_module

        d = tempfile.mkdtemp(prefix="uhfk_l2_")
        try:
            _write_l2_chain(d)
            info_mode = {
                'mode': 'UHFk',
                'param': {'Ncond': 1, 'T': 0.5, 'CellShape': [2, 1, 1],
                          'IterationMax': 1, 'EPS': 8, 'Mix': 0.5,
                          'RndSeed': 1,
                          # +-1 hops on L=2 exceed the strict range check
                          # (width 3 > 2); this test targets the wrap-around
                          # accumulation of the relaxed path
                          'trustme_interaction_range': True},
            }
            info_input = {
                'path_to_input': d,
                'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                                'Transfer': 'transfer.dat'},
            }
            read_io = read_input_k.QLMSkInput(info_input)
            info_log = {"print_level": 0, "print_step": 100}
            solver = uhfk_module.UHFk(read_io.get_param("ham"), info_log,
                                      info_mode)
            solver._make_ham_trans()
            bands = np.sort(np.unique(np.round(
                np.linalg.eigvalsh(solver.ham_trans).real.ravel(), 8)))
        finally:
            shutil.rmtree(d, ignore_errors=True)

        np.testing.assert_allclose(
            bands, [-2.0, 2.0], atol=1e-8,
            err_msg="L=2 ring with t=1: H(k)=2t cos(k) gives bands +-2")


class TestWraparoundExternRPA(unittest.TestCase):
    def test_extern_wraparound_bonds_are_summed(self):
        from hwave.solver.rpa import Interaction

        class _FakeLattice:
            shape = (2, 1, 1)
            nvol = 2

        obj = object.__new__(Interaction)
        obj.lattice = _FakeLattice()
        obj.norb = 1
        obj.enable_spin_orbital = False
        obj.param_ham = {
            "Transfer": {((1, 0, 0), (0, 0)): 1.0,
                         ((-1, 0, 0), (0, 0)): 1.0},
            "Extern": {((1, 0, 0), (0, 0)): 0.5,
                       ((-1, 0, 0), (0, 0)): 0.5},
        }
        obj._make_ham_trans()

        # R=+1 and R=-1 wrap onto slot 1 on an L=2 ring and must be summed
        np.testing.assert_allclose(obj.ham_trans_r[1, 0, 0], 2.0, atol=1e-12)
        np.testing.assert_allclose(
            obj.ham_extern_r[1, 0, 0], 1.0, atol=1e-12,
            err_msg="Extern wrap-around bonds must accumulate like Transfer")


class TestWraparoundHoppingRPA(unittest.TestCase):
    def test_l2_ring_band_is_pm_2t(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as solver_rpa

        d = tempfile.mkdtemp(prefix="rpa_l2_")
        try:
            _write_l2_chain(d)
            info_mode = {
                'mode': 'RPA',
                'param': {'T': 1.0, 'mu': 0.0, 'CellShape': [2, 1, 1],
                          'Nmat': 16},
                'calc_scheme': 'reduced',
            }
            info_input = {
                'path_to_input': d,
                'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                                'Transfer': 'transfer.dat'},
            }
            read_io = read_input_k.QLMSkInput(info_input)
            solver = solver_rpa.RPA(read_io.get_param("ham"), {}, info_mode)
            hq = solver.ham_info.ham_trans_q
            bands = np.sort(np.unique(np.round(
                np.linalg.eigvalsh(hq).real.ravel(), 8)))
        finally:
            shutil.rmtree(d, ignore_errors=True)

        np.testing.assert_allclose(
            bands, [-2.0, 2.0], atol=1e-8,
            err_msg="L=2 ring with t=1: H(k)=2t cos(k) gives bands +-2")


if __name__ == "__main__":
    unittest.main()
