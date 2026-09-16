"""D5: the bond gate's channel-0 second order is the local kernel under
'local' (the rev-19 expression under 'takimoto'); G0 (declared-zero
off-site = standalone Hartree-Fock FLEX) holds under both values; the
effective interaction keeps the frequency-reflection Hermiticity D-7 once
the mixed second-order blocks carry the pair permutation of spec
2026-09-16 R3 (issue #192); the memory table carries the factor row."""
import contextlib
import tempfile
import unittest
from unittest import mock

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

    def _bond_w(self, identity_permutation=False):
        """``(W, nd)``: the FULL bond-resolved ``W`` (every block,
        ``(nmat, nvol, ND, ND)``) of the module's fixture and the pair
        dimension of one channel block, built through the same recipe as
        :meth:`test_channel0_block_equals_the_standalone_kernel`.

        The fixture is ``tests/rpa/input_2orb``'s ``coulombinter.dat``, whose
        off-site content includes the INTER-ORBITAL rows ``v_12(-x) = 1`` /
        ``v_21(+x) = 1``: each bond channel's ``(norb, norb)`` coefficient
        block is orbital-off-diagonal and NOT symmetric under the orbital
        swap, which is what makes the mixed blocks' pair permutation visible
        here at all (it is the identity on an orbital-diagonal bond).

        With ``identity_permutation`` the permutation is replaced by the
        identity -- the pre-R3 behaviour -- so the caller can see whether it
        is load-bearing."""
        from hwave.solver import flex_bond
        s, r = _flex({"flex_second_order": "local", "IterationMax": 1})
        gi = r.get_param("green")
        s._phase_b_reset(gi); s._phase_b_preflight(gi); s._calc_epsilon_k({})
        beta = 0.5
        nmat, nvol, norb = s.nmat, s.lattice.nvol, s.norb
        nd = norb * norb
        G = s._calc_dressed_green(beta, 0.1, np.zeros((1, nmat, nvol, norb, norb), complex))
        B = s._bond_view.n_channels
        patch = (mock.patch.object(flex_bond, "mixed_pair_permutation",
                                   lambda nb, ndd, nrb: np.arange(nb * ndd))
                 if identity_permutation else contextlib.nullcontext())
        with flex_bond.BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, None, beta, s._bond_view, (4, 4, 1), 1)
            with patch:
                flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on,
                                            C_on=s._bond_C_on, nb=nmat, output_full=False,
                                            nmat=nmat, nvol=nvol, nd=nd, spatial_shape=(4, 4, 1),
                                            factors=s._second_order_factors,
                                            second_order="local")
            return np.array(store.get_freq_batch("W", 0, nmat)), nd

    def test_w_hermiticity_and_the_correction_is_load_bearing(self):
        """D-7 on the FULL effective interaction: ``W(q, l)^dagger =
        W(q, -l)`` on the store's bosonic frequency index, where ``-l`` is
        ``(nmat - l) % nmat`` (the convention pinned by
        ``tests/test_flex_bond_sigma_transport.py``'s reflection and by the
        local kernel's own ``W2[(nmat - l) % nmat]`` symmetry).

        This is what carries ``Sigma(k, i w)^dagger = Sigma(k, -i w)``
        through the bond transport, so a pair permutation applied to only
        ONE of the two pair axes -- or to a non-Hermitian combination --
        would break it even while the coefficient gates still passed on a
        real-coefficient fixture. The second leg asserts the permutation is
        load-bearing here: with it replaced by the identity (the pre-R3
        behaviour, issue #192) the MIXED blocks move by more than 1e-3 of
        their own size, and nothing outside them moves at all.

        Why the mixed blocks' own size is the yardstick for that second
        leg and ``max|W|`` is not: ``max|W|`` is the channel-0 block, the
        RPA-resummed ring, and the mixed blocks are some three orders
        smaller than it on every fixture (the bond-channel entries of the
        vertex are the bond-DIAGONAL ``v_ab(R)``, while channel 0 carries
        the full q-dependent density slot). Measured here: ``max|W|``
        1.9e+1, mixed blocks 1.6e-2, and the permutation moves the mixed
        blocks by 3.4e-2 of themselves -- 2.8e-5 of ``max|W|``. Normalising
        the second leg by ``max|W|`` would therefore be asking for a ratio
        no fixture can produce; the third assertion (nothing outside the
        mixed strips moves) is what makes this normalisation honest."""
        W, nd = self._bond_w()
        nmat = W.shape[0]
        scale = np.abs(W).max()
        self.assertGreater(scale, 1e-6)                              # anti-vacuity
        dev = max(np.abs(W[l].conj().swapaxes(-1, -2) - W[(-l) % nmat]).max()
                  for l in range(nmat))
        self.assertLess(dev / scale, 1e-12,
                        "the bond-resolved W is not Hermitian under the frequency "
                        "reflection: {:.3e} of its own size".format(dev / scale))
        W_identity, _ = self._bond_w(identity_permutation=True)
        delta = np.abs(W - W_identity)
        mixed = max(np.abs(W[:, :, :nd, nd:]).max(), np.abs(W[:, :, nd:, :nd]).max())
        self.assertGreater(mixed, 1e-6 * scale,                      # anti-vacuity
                           "the mixed blocks of W are empty on this fixture")
        moved = delta.max() / mixed
        self.assertGreater(moved, 1e-3,
                           "the mixed blocks' pair permutation does not change W on this "
                           "fixture ({:.3e} of the mixed blocks' own size), so the "
                           "Hermiticity assertion above says nothing about it".format(moved))
        # and it moves ONLY the mixed blocks: the channel-0 block and the
        # bond-bond blocks are outside the mask the permutation is applied
        # through, so they must be untouched
        outside = delta.copy()
        outside[:, :, :nd, nd:] = 0.0
        outside[:, :, nd:, :nd] = 0.0
        self.assertLess(outside.max(), 1e-14 * scale,
                        "the pair permutation changed W outside the mixed (channel-0 x bond) "
                        "blocks by {:.3e} of max|W|".format(outside.max() / scale))

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
