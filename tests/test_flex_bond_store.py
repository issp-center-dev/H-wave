"""BondBlockStore and the dynamic bubble assembly (spec 2026-09-06
sections 3.1 and 3.5): pair/frequency-batch views, detach/release,
rollback, the static-slice pin against bond_bubble_static, and the
every-frequency channel-0 identity with the general bubble under a
nonzero coeff_tail."""
import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"


def _flex(param_extra):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombInter": "coulombinter.dat"}
    r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
    par.update(param_extra)
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


class TestBondBlockStore(unittest.TestCase):

    def test_views_detach_release(self):
        from hwave.solver.flex_bond import BondBlockStore
        nmat, nvol, nd, B = 4, 3, 4, 2
        with BondBlockStore(nmat, nvol, B * nd, nd, ("chibar", "W")) as store:
            blk = np.arange(nmat * nvol * nd * nd, dtype=float).reshape(nmat, nvol, nd, nd) + 0j
            store.put_pair("chibar", 1, 0, blk)
            np.testing.assert_array_equal(store.get_pair("chibar", 1, 0), blk)
            full = store.get_freq_batch("chibar", 1, 3)
            self.assertEqual(full.shape, (2, nvol, B * nd, B * nd))
            self.assertTrue(full.flags.c_contiguous)
            np.testing.assert_array_equal(full[:, :, nd:, :nd], blk[1:3])
            store.put_freq_batch("W", 0, 4, np.ones((4, nvol, B * nd, B * nd), complex))
            w = store.detach("W")
            self.assertTrue(w.flags.owndata)
            with self.assertRaises(KeyError):
                store.get_pair("W", 0, 0)
        with self.assertRaises(RuntimeError):
            store.get_pair("chibar", 0, 0)                       # released
        self.assertEqual(np.abs(w).max(), 1.0)                  # detached array survives

    def test_release_on_exception_and_rollback(self):
        from hwave.solver.flex_bond import BondBlockStore
        try:
            with BondBlockStore(2, 2, 2, 1, ("chibar",)) as store:
                raise ValueError("boom")
        except ValueError:
            pass
        with self.assertRaises(RuntimeError):
            store.get_freq_batch("chibar", 0, 1)
        with self.assertRaises(ValueError):
            BondBlockStore(2, 2, 2, 1, ("chibar", "chibar"))    # duplicate name refused


class TestBubbleAssembly(unittest.TestCase):

    def _prepared(self, coeff_tail):
        s, r = _flex({"coeff_tail": coeff_tail})
        gi = r.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s.solve(gi, "tests/flex/output")           # gate off: builds green0 / tails
        return s, gi

    def test_static_slice_equals_bond_bubble_static(self):
        from hwave.solver import bond_channels as bc, bubble
        from hwave.solver.flex_bond import BondBlockStore, assemble_bubble
        s, gi = self._prepared(0.0)
        topo = bc.resolve_bond_topology(s.ham_info.param_ham, np.eye(3), s.norb,
                                        active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        nd, nvol, nmat = s.norb ** 2, s.lattice.nvol, s.nmat
        ND = view.n_channels * nd
        with BondBlockStore(nmat, nvol, ND, nd, ("chibar",)) as store:
            assemble_bubble(store, s.green0, s.green0_tail, 0.5, view, (4, 4, 1), 1)
            stat = bubble.bond_bubble_static(s.green0, s.green0_tail, 0.5, view, spatial_shape=(4, 4, 1))
            np.testing.assert_allclose(store.get_freq_batch("chibar", nmat // 2, nmat // 2 + 1)[0],
                                       np.asarray(stat), rtol=0, atol=1e-13)

    def test_channel_zero_equals_general_bubble_at_every_frequency_with_tail(self):
        from hwave.solver import bond_channels as bc
        from hwave.solver.flex_bond import BondBlockStore, assemble_bubble
        s, gi = self._prepared(1.0)
        topo = bc.resolve_bond_topology(s.ham_info.param_ham, np.eye(3), s.norb,
                                        active_types=bc._LONGITUDINAL_ACTIVE_TYPES)
        view = bc.BondSetView(topo)
        nd, nvol, nmat, norb = s.norb ** 2, s.lattice.nvol, s.nmat, s.norb
        ND = view.n_channels * nd
        # the loop's own green_scf / green0_tail: bare G minus the tail (first iteration)
        beta = 0.5
        green_kw = s._calc_dressed_green(beta, s.mu, np.zeros((1, nmat, nvol, norb, norb), complex))
        iomega = (np.arange(nmat) * 2 + 1 - nmat) * np.pi / beta
        ev = s.H0_eigenvector
        VVt = ev @ np.conj(ev).swapaxes(-2, -1)
        tail_w = ((1.0 / (1j * iomega))[None, :, None, None, None] * VVt[:, None])
        green_scf = green_kw - tail_w
        ref = s._calc_chi0q(green_scf, s.green0_tail, beta)[0]       # (nmat, nvol, norb, norb, norb, norb)
        with BondBlockStore(nmat, nvol, ND, nd, ("chibar",)) as store:
            assemble_bubble(store, green_scf, s.green0_tail, beta, view, (4, 4, 1), 1)
            ch0 = store.get_freq_batch("chibar", 0, nmat)[:, :, :nd, :nd]
        np.testing.assert_allclose(ch0.reshape(nmat, nvol, norb, norb, norb, norb), ref, rtol=0, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
