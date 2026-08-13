#!/usr/bin/env python3
"""calc_scheme='squashed' is rejected at construction (issue #144).

squashed computed the same susceptibility as reduced at several times the
cost; its extra spin-resolved output slots were structurally zero. Removed
in 2.0. The legacy PROVENANCE token stays accepted (see
tests/test_rpa_green_info_reuse.py): rejection of the configuration and
acceptance of old files' metadata are separate concerns.
"""
import unittest


def _ham(norb=1):
    import numpy as np
    return {'Geometry': {'norb': norb, 'rvec': np.eye(3),
                         'center': [[0, 0, 0]] * norb},
            'Transfer': {((0, 0, 0), (0, 0)): 1.0},
            'CoulombIntra': {((0, 0, 0), (0, 0)): 4.0}}


def _mode(mode, extra=None):
    par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 1, 1],
           'SubShape': [1, 1, 1], 'Nmat': 4}
    if extra:
        par.update(extra)
    return {'mode': mode, 'param': par, 'enable_spin_orbital': False,
            'calc_scheme': 'squashed'}


class TestSquashedRejected(unittest.TestCase):

    def test_rpa_rejects_squashed_at_construction(self):
        import hwave.solver.rpa as rpa_mod
        with self.assertRaises(ValueError) as cm:
            rpa_mod.RPA(_ham(), {}, _mode('RPA'))
        self.assertIn("calc_scheme='squashed' was removed in H-wave 2.0",
                      str(cm.exception))
        self.assertIn("reduced", str(cm.exception))
        self.assertIn("(l,q,a,b)", str(cm.exception))

    def test_flex_rejects_squashed_at_construction(self):
        import hwave.solver.flex as flex_mod
        with self.assertRaises(ValueError) as cm:
            flex_mod.FLEX(_ham(), {},
                          _mode('FLEX', {'IterationMax': 1, 'Mix': 1.0}))
        self.assertIn("calc_scheme='squashed' was removed in H-wave 2.0",
                      str(cm.exception))


class TestAutoNeverResolvesToSquashed(unittest.TestCase):
    """auto must keep resolving to reduced/general only (regression)."""

    def _scheme_for(self, ham):
        import hwave.solver.rpa as rpa_mod
        info = _mode('RPA')
        info['calc_scheme'] = 'auto'
        return rpa_mod.RPA(ham, {}, info).calc_scheme

    def test_density_only_resolves_to_reduced(self):
        self.assertEqual(self._scheme_for(_ham()), 'reduced')

    def test_exchange_resolves_to_general(self):
        ham = _ham(norb=2)
        ham['Exchange'] = {((0, 0, 0), (0, 1)): 0.5,
                          ((0, 0, 0), (1, 0)): 0.5}
        self.assertEqual(self._scheme_for(ham), 'general')

    def test_ring_ladder_resolves_to_general(self):
        import hwave.solver.rpa as rpa_mod
        info = _mode('RPA')
        info['calc_scheme'] = 'auto'
        info['calc_type'] = 'ring+ladder'
        self.assertEqual(rpa_mod.RPA(_ham(), {}, info).calc_scheme,
                         'general')


if __name__ == "__main__":
    unittest.main()
