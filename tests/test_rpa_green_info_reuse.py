#!/usr/bin/env python3
"""Reusing a green_info dict across RPA solves (issue #109).

A solve stores its chi0q back into green_info. Handing that dict to a
second (fresh) solver is ordinary library usage -- e.g. scanning
interaction parameters on a fixed bubble -- and used to crash with
AttributeError: the externally-supplied-chi0q branch reached the
spin_mode dispatch before anything set spin_mode on the new instance
(the file-based chi0q_init route runs the shape detection in
_read_chi0q; the in-memory route ran none).

The reused bubble is the one the first solver computed, so the second
solve must also REPRODUCE the first solver's chiq exactly.
"""

import os
import unittest

import numpy as np


def _make(scheme, interactions, path='tests/rpa/input_2orb'):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as rpa_mod

    idict = {'path_to_input': path, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat'}
    idict.update(interactions)
    r = read_input_k.QLMSkInput({'path_to_input': path,
                                 'interaction': idict})
    par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
           'SubShape': [1, 1, 1], 'Nmat': 32}
    sv = rpa_mod.RPA(r.get_param("ham"), {},
                     {'mode': 'RPA', 'param': par,
                      'enable_spin_orbital': False,
                      'calc_scheme': scheme, 'calc_type': 'ring'})
    return sv, r


class TestGreenInfoReuse(unittest.TestCase):

    def _roundtrip(self, scheme):
        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)

        sv1, r1 = _make(scheme, inter)
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')
        self.assertIn('chi0q', g)
        chiq_first = np.array(g['chiq'], copy=True)

        sv2, _ = _make(scheme, inter)
        sv2.solve(g, 'tests/rpa/output')   # crashed with AttributeError
        self.assertEqual(sv2.spin_mode, sv1.spin_mode)
        np.testing.assert_array_equal(np.asarray(g['chiq']), chiq_first)

    def test_general_scheme_roundtrip(self):
        self._roundtrip('general')

    def test_reduced_scheme_roundtrip(self):
        self._roundtrip('reduced')

    def test_squashed_scheme_roundtrip(self):
        self._roundtrip('squashed')

    def test_spin_diag_blocks_are_detected_in_memory(self):
        """A two-block (spin-diag) bubble handed over in memory must select
        spin-diag mode, exactly as the file-based chi0q_init route does."""
        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)

        sv1, r1 = _make('general', inter)
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')
        chiq_free = np.array(g['chiq'], copy=True)

        sv2, r2 = _make('general', inter)
        g2 = r2.get_param("green")
        g2['chi0q'] = np.stack([np.asarray(g['chi0q'])] * 2, axis=0)
        sv2.solve(g2, 'tests/rpa/output')
        self.assertEqual(sv2.spin_mode, "spin-diag")
        # two identical blocks describe the same spin-symmetric system
        np.testing.assert_array_equal(np.asarray(g2['chiq']), chiq_free)


if __name__ == "__main__":
    unittest.main()
