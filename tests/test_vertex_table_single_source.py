#!/usr/bin/env python3
"""Both vertex builders derive from the one adjudicated table.

The pair-space S/C matrices (hwave.sc) and the ring's spin-resolved
cross-slot correction (hwave.solver.rpa) encode the same physics in two
representations, tied by the channel decomposition W_same = (C - S)/2,
W_cross = (C + S)/2. Historically the two drifted (#96, #104, #113);
they now read one table (hwave.solver.vertex_table). These tests measure
the consistency from the ASSEMBLED objects, not the table -- so a future
edit that bypasses the table on either side still fails here.
"""

import unittest

import numpy as np


class TestAssembledConsistency(unittest.TestCase):

    def _sc_cross(self, itype, J=0.7):
        """S and C at the (01, 01) cross slot from the sc builders."""
        import hwave.sc as sc

        k = np.array([0.0])
        entries = {((0, 0, 0), (0, 1)): J, ((0, 0, 0), (1, 0)): J}
        ik = sc._build_interaction_k(k, k, k, {itype: entries}, 2)
        S, C = sc._build_sc_matrices_all_q(ik, 2, 1, 1, 1)
        return S[0, 0, 0, 1, 1], C[0, 0, 0, 1, 1]

    def _rpa_fierz(self, itype, J=0.7):
        """Same-spin and cross-spin cross-slot entries of the ring's
        Fierz tensor, read from the assembled solver."""
        import hwave.solver.rpa as rpa_mod

        param_ham = {'Geometry': {'norb': 2, 'rvec': np.eye(3),
                                  'center': [[0, 0, 0]] * 2},
                     'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                     itype: {((0, 0, 0), (0, 1)): J,
                             ((0, 0, 0), (1, 0)): J}}
        info = {'mode': 'RPA', 'param': {'T': 2.0, 'filling': 0.5,
                'CellShape': [1, 1, 1], 'SubShape': [1, 1, 1], 'Nmat': 4},
                'enable_spin_orbital': False, 'calc_scheme': 'general',
                'calc_type': 'ring'}
        sv = rpa_mod.RPA(param_ham, {}, info)
        f = sv.ham_info.ham_fierz_q
        if f is None:
            return 0.0, 0.0
        f = np.asarray(f).reshape(-1, 4, 4, 4, 4)[0]
        # spin-major nd = 4: up = {0: 0, 1: 1}, down = {0: 2, 1: 3};
        # entries sit at (a, b, a, b) per spin block (see the appender)
        same = f[0, 1, 0, 1]          # (up a, up b, up a, up b)
        cross = f[2, 3, 0, 1]         # (dn a, dn b, up a, up b)
        return same, cross

    def test_ring_fierz_is_the_channel_decomposition_of_sc(self):
        """For every type with cross-slot content: measure S/C from the
        sc builders and the same-spin/cross-spin entries from the ring,
        and require W_same == (C - S)/2, W_cross == (C + S)/2."""
        for itype in ("CoulombInter", "Hund", "Ising", "Exchange"):
            with self.subTest(interaction=itype):
                S, C = self._sc_cross(itype)
                same, cross = self._rpa_fierz(itype)
                self.assertAlmostEqual(
                    complex(same), complex((C - S) / 2.0), places=12,
                    msg="%s same-spin" % itype)
                self.assertAlmostEqual(
                    complex(cross), complex((C + S) / 2.0), places=12,
                    msg="%s cross-spin" % itype)

    def test_types_without_cross_content_have_no_fierz(self):
        for itype in ("PairHop", "PairLift", "CoulombIntra"):
            with self.subTest(interaction=itype):
                entries = ({((0, 0, 0), (0, 0)): 0.7,
                            ((0, 0, 0), (1, 1)): 0.7}
                           if itype == "CoulombIntra" else
                           {((0, 0, 0), (0, 1)): 0.7,
                            ((0, 0, 0), (1, 0)): 0.7})
                import hwave.solver.rpa as rpa_mod

                param_ham = {'Geometry': {'norb': 2, 'rvec': np.eye(3),
                                          'center': [[0, 0, 0]] * 2},
                             'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                             itype: entries}
                info = {'mode': 'RPA',
                        'param': {'T': 2.0, 'filling': 0.5,
                                  'CellShape': [1, 1, 1],
                                  'SubShape': [1, 1, 1], 'Nmat': 4},
                        'enable_spin_orbital': False,
                        'calc_scheme': 'general', 'calc_type': 'ring'}
                sv = rpa_mod.RPA(param_ham, {}, info)
                self.assertIsNone(sv.ham_info.ham_fierz_q)


if __name__ == "__main__":
    unittest.main()
