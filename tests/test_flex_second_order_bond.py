"""D5: the bond gate's channel-0 second order is the local kernel under
'local' (the rev-19 expression under 'takimoto'); G0 (declared-zero
off-site = standalone Hartree-Fock FLEX) holds under both values; the
memory table carries the factor row."""
import tempfile
import unittest

import numpy as np

from tests.test_flex_bond_gate import _flex, _collect, _zeroed_offsite_copy


class TestBondGate(unittest.TestCase):

    def test_g0_under_both_kernels(self):
        for so in ("local", "takimoto"):
            with self.subTest(so=so), tempfile.TemporaryDirectory() as inp, tempfile.TemporaryDirectory() as out:
                _zeroed_offsite_copy(inp)
                recs = {}
                for gate in (True, False):
                    s, r = _flex({"flex_second_order": so}, path=inp, gate=gate)
                    recs[gate] = _collect(s, r.get_param("green"), out)
                self.assertGreaterEqual(len(recs[True]), 2)
                self.assertEqual(len(recs[True]), len(recs[False]))
                for a, b in zip(recs[True], recs[False]):
                    for k in ("static_new", "fluct_new", "chi0q", "chiq_s", "chiq_c"):
                        np.testing.assert_allclose(a[k], b[k], rtol=0, atol=1e-12, err_msg=(so, k))

    def test_channel0_block_equals_the_standalone_kernel(self):
        from hwave.solver import flex_bond
        from hwave.solver.second_order import dense_w2
        s, r = _flex({"flex_second_order": "local", "IterationMax": 1})
        gi = r.get_param("green")
        s._phase_b_reset(gi); s._phase_b_preflight(gi); s._calc_epsilon_k({})
        beta = 0.5
        nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
        nd = norb * norb
        G = s._calc_dressed_green(beta, 0.1, np.zeros((1, nmat, nvol, norb, norb), complex))
        B = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, None, beta, s._bond_view, (4, 4, 1), 1)
            cb = np.array(store.get_freq_batch("chibar", 0, nmat))
            flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on, C_on=s._bond_C_on,
                                        nb=nmat, output_full=False, nmat=nmat, nvol=nvol, nd=nd,
                                        spatial_shape=(4, 4, 1), factors=s._second_order_factors,
                                        second_order="local")
            W00 = np.array(store.get_freq_batch("W", 0, nmat))[:, :, :nd, :nd]
        S, C = s._bond_S, s._bond_C
        I = np.eye(B * nd)
        chi_s = np.linalg.solve(I - cb @ S, cb); chi_c = np.linalg.solve(I + cb @ C, cb)
        ring = (1.5 * S @ (chi_s - cb) @ S + 0.5 * C @ (chi_c - cb) @ C)[:, :, :nd, :nd]
        w2 = dense_w2(cb[:, :, :nd, :nd], s._second_order_factors)
        np.testing.assert_allclose(W00, ring + w2, rtol=1e-12, atol=1e-13)

    def test_local_without_factors_is_refused(self):
        from hwave.solver.flex_bond import BondBlockStore, dress_and_build_w
        nmat, nvol, nd, B = 2, 2, 1, 1
        S = np.zeros((nvol, B * nd, B * nd), complex)
        with BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            with self.assertRaises(ValueError) as cm:
                dress_and_build_w(store, S, S, S_on=S, C_on=S, nb=nmat, output_full=False,
                                  nmat=nmat, nvol=nvol, nd=nd, spatial_shape=(2, 1, 1),
                                  factors=None, second_order="local")
        self.assertIn("flex_second_order", str(cm.exception))

    def test_memory_row(self):
        from hwave.solver.flex_bond import estimate_bond_memory
        est = estimate_bond_memory(nmat=8, nvol=16, norb=2, B=3, depth=0, output_full=False, split_seed=False,
                                   n_types=1, freq_batch=None, cap_gb=8.0, mixing="linear", factor_bytes=1234)
        self.assertEqual(est["persistent_rows"]["second_order_factors"], 1234)
        self.assertEqual(est["persistent"], sum(est["persistent_rows"].values()))

    def test_memory_row_is_fed_from_the_solver(self):
        s, r = _flex({"flex_second_order": "local", "IterationMax": 1})
        gi = r.get_param("green")
        s._phase_b_reset(gi)
        s._phase_b_preflight(gi)
        self.assertGreater(s._second_order_factors.nbytes, 0)
        self.assertEqual(s._bond_est["persistent_rows"]["second_order_factors"],
                         s._second_order_factors.nbytes)


if __name__ == "__main__":
    unittest.main()
