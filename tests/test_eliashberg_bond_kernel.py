"""Bond-resolved dynamic pairing kernel (spec 10.1)."""
import unittest

import numpy as np

from tests.eliashberg_bond_fixtures import physical_fixture


def _dev(fx):
    from hwave.solver.flex_bond import BondDeviceContext
    return BondDeviceContext(np, fx["S"], fx["C"])


class TestAccumulator(unittest.TestCase):
    def test_accumulator_stages_agree(self):                      # 10.1.7
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, COEFF
        fx = physical_fixture(norb=2, shape=(4, 2, 1), nmat=4)
        for eta in ("singlet", "triplet"):
            with _dev(fx) as dev:
                acc = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                            nvol=fx["nvol"], nd=fx["nd"],
                                            spatial_shape=fx["spatial_shape"])
                acc.add_dressed(fx["store"], cond_tol=1e-6)
                v1 = acc.finish()[eta]
                acc2 = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                             nvol=fx["nvol"], nd=fx["nd"],
                                             spatial_shape=fx["spatial_shape"])
                acc2.add_channel("spin", fx["store"], "chi_s_w")
                acc2.add_channel("charge", fx["store"], "chi_c_w")
                v2 = acc2.finish()[eta]
                g1 = v1.source.get_freq_batch(v1.slot, 0, fx["nmat"])
                g2 = v2.source.get_freq_batch(v2.slot, 0, fx["nmat"])
                np.testing.assert_array_equal(g1, g2)
                cs, cc = COEFF[eta]
                S, C = fx["S"], fx["C"]
                ref = cs * (S[None] @ fx["chi_s"] @ S[None]) + cc * (C[None] @ fx["chi_c"] @ C[None])
                np.testing.assert_allclose(g1, ref, rtol=0, atol=1e-13)
                # a missing stage is refused
                acc3 = PairVertexAccumulator(dev, pairing_types=(eta,), nb=3, nmat=fx["nmat"],
                                             nvol=fx["nvol"], nd=fx["nd"],
                                             spatial_shape=fx["spatial_shape"])
                acc3.add_channel("spin", fx["store"], "chi_s_w")
                with self.assertRaisesRegex(ValueError, "charge stage"):
                    acc3.finish()

    def test_uniform_refuses_two_types_and_non_finite(self):
        from hwave.solver.eliashberg_bond import PairVertexAccumulator, ArrayBlockSource
        fx = physical_fixture(norb=1, shape=(2, 2, 1), nmat=4)
        with _dev(fx) as dev:
            with self.assertRaisesRegex(ValueError, "exactly one pairing type"):
                PairVertexAccumulator(dev, pairing_types=("singlet", "triplet"), nb=2,
                                      nmat=4, nvol=4, nd=1, spatial_shape=(2, 2, 1))
            bad = fx["chi_s"].copy()
            bad[1, 0, 0, 0] = np.nan
            src = ArrayBlockSource({"chi_s_w": bad}, 1)
            acc = PairVertexAccumulator(dev, pairing_types=("singlet",), nb=2, nmat=4, nvol=4,
                                        nd=1, spatial_shape=(2, 2, 1))
            with self.assertRaisesRegex(ValueError, r"non-finite .* batch \[0, 2\)"):
                acc.add_channel("spin", src, "chi_s_w")


class TestAdmission(unittest.TestCase):
    def test_admission_table(self):                              # 10.1.11
        from hwave.solver.eliashberg_bond import estimate_pair_memory
        tab = estimate_pair_memory(nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=3,
                                   num_eigenvalues=2, residency="auto", ir=False, L_B=0,
                                   n_channels=1, in_process=False)
        nd, ND, vec = 4, 12, 4 * 16 * 8 * 16
        blk = 8 * 16 * 16 * 16
        self.assertEqual(tab.rows["G2"], 16 * 16 * 8 * 16)
        self.assertEqual(tab.rows["gap_work"], (3 + 3) * vec)
        self.assertEqual(tab.rows["hoisted_blocks"], 9 * blk)
        self.assertEqual(tab.rows["stream_workspace"], 3 * blk)
        self.assertEqual(tab.rows["eigen_vectors"], (2 * 2 + 1 + 2) * vec)
        self.assertEqual(tab.rows["vertex_slot"], 8 * 16 * ND * ND * 16)   # post-processing uniform
        h, d = tab.need("device")
        self.assertEqual(d, int(1.25 * (tab.rows["G2"] + tab.rows["gap_work"] + 9 * blk + 3 * blk)))
        self.assertEqual(h, int(1.25 * (tab.rows["vertex_slot"] + tab.rows["eigen_vectors"])))
        big = 10 ** 12
        self.assertEqual(tab.choose("auto", host_cap=big, device_cap=big), "device")
        self.assertEqual(tab.choose("auto", host_cap=big, device_cap=tab.need("host")[1]), "host")
        self.assertEqual(tab.choose("auto", host_cap=tab.need("stream")[0],
                                    device_cap=tab.need("stream")[1]), "stream")
        with self.assertRaisesRegex(MemoryError, "stream"):
            tab.choose("auto", host_cap=1, device_cap=1)
        with self.assertRaisesRegex(MemoryError, "device"):
            tab.choose("device", host_cap=big, device_cap=1)
        tab_ip = estimate_pair_memory(nmat=8, ntau=0, nfreq=8, nvol=16, norb=2, B=3,
                                      num_eigenvalues=2, residency="auto", ir=False, L_B=0,
                                      n_channels=1, in_process=True)
        self.assertEqual(tab_ip.rows["vertex_slot"], 0)
        tab_ir = estimate_pair_memory(nmat=8, ntau=5, nfreq=6, nvol=16, norb=2, B=3,
                                      num_eigenvalues=2, residency="auto", ir=True, L_B=4,
                                      n_channels=2, in_process=True)
        self.assertEqual(tab_ir.rows["coefficients_build"], 2 * 2 * 5 * 16 * ND * ND * 16)
        self.assertEqual(tab_ir.rows["coefficients"], 2 * 5 * 16 * ND * ND * 16)
        self.assertEqual(tab_ir.rows["hoisted_blocks"], 9 * 5 * 16 * 16 * 16)


if __name__ == "__main__":
    unittest.main()
