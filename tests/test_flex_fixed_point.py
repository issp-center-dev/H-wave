#!/usr/bin/env python3
"""The converged FLEX solution is a fixed point of the SCF map.

Feed a converged self-energy back as the warm start and run a SINGLE
iteration with full mixing: the output must reproduce the converged
sigma and susceptibilities. This pins the consistency of the whole
pipeline -- if the first-iteration path ever diverges from the in-loop
path (normalization drift, an asymmetric sigma update, a vertex applied
differently on entry), the fixed point moves and these tests fail.

The matrix covers every interaction type FLEX accepts, under every
calc_scheme. The reduced/squashed schemes REJECT Exchange and PairHop
(one policy since #107: their vertex has no density-diagonal content,
so those schemes would drop them entirely); the rejection is pinned
here alongside the fixed-point cells.
"""

import os
import unittest

import numpy as np


def _run(scheme, interactions, path, itmax, eps, mix, sigma_init=None):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod

    idict = {'path_to_input': path, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat'}
    idict.update(interactions)
    r = read_input_k.QLMSkInput({'path_to_input': path,
                                 'interaction': idict})
    filling = 0.75 if 'input_2orb' not in path else 0.5
    pf = {'T': 2.0, 'filling': filling, 'CellShape': [4, 4, 1],
          'SubShape': [1, 1, 1], 'Nmat': 32,
          'IterationMax': itmax, 'Mix': mix, 'EPS': eps}
    fx = flex_mod.FLEX(r.get_param("ham"), {},
                       {'mode': 'FLEX', 'param': pf,
                        'enable_spin_orbital': False,
                        'calc_scheme': scheme})
    g = r.get_param("green")
    if sigma_init is not None:
        g['sigma_init'] = sigma_init
    fx.solve(g, 'tests/rpa/output')
    return g


CELLS = [
    ("U 1orb", 'tests/rpa/input',
     {'CoulombIntra': 'coulombintra.dat'}),
    ("V offsite 1orb", 'tests/rpa/input',
     {'CoulombInter': 'coulombinter.dat'}),
    ("U+V+J 2orb", 'tests/rpa/input_2orb',
     {'CoulombIntra': 'coulombintra.dat',
      'CoulombInter': 'onsite_inter.dat',
      'Hund': 'hund_onsite.dat'}),
    ("Ising 2orb", 'tests/rpa/input_2orb',
     {'Ising': 'hund_onsite.dat'}),
    ("Exchange 2orb", 'tests/rpa/input_2orb',
     {'Exchange': 'hund_onsite.dat'}),
    ("PairHop 2orb", 'tests/rpa/input_2orb',
     {'PairHop': 'hund_onsite.dat'}),
]

ATOL = 1.0e-8


class TestFlexFixedPoint(unittest.TestCase):

    def _check_cell(self, scheme, name, path, interactions):
        os.makedirs('tests/rpa/output', exist_ok=True)
        converged = _run(scheme, interactions, path,
                         itmax=300, eps=10, mix=0.5)
        sigma = np.array(converged['sigma'], copy=True)
        restart = _run(scheme, interactions, path,
                       itmax=1, eps=1, mix=1.0, sigma_init=sigma)
        ds = float(np.max(np.abs(np.asarray(restart['sigma']) - sigma)))
        self.assertLess(
            ds, ATOL,
            "%s / %s: one full-mixing sweep moved the converged sigma by "
            "%.3e -- the SCF map is inconsistent with itself"
            % (scheme, name, ds))
        for key in ('chiq_s', 'chiq_c'):
            d = float(np.max(np.abs(np.asarray(restart[key])
                                    - np.asarray(converged[key]))))
            self.assertLess(
                d, ATOL,
                "%s / %s: %s of the warm-started sweep differs from the "
                "converged one by %.3e" % (scheme, name, key, d))

    def _check_or_reject(self, scheme, name, path, interactions):
        if name in ("Exchange 2orb", "PairHop 2orb"):
            # ONE policy since #107: the density-diagonal schemes reject
            # Exchange and PairHop for both solvers -- their adjudicated
            # vertex has no density-diagonal content, so the schemes
            # would drop them entirely (zero effect, not an
            # approximation). Previously reduced raised SystemExit,
            # squashed silently accepted, FLEX warned.
            with self.assertRaises(ValueError) as cm:
                self._check_cell(scheme, name, path, interactions)
            self.assertIn('general', str(cm.exception))
            return
        self._check_cell(scheme, name, path, interactions)

    def test_reduced(self):
        for name, path, interactions in CELLS:
            with self.subTest(interaction=name):
                self._check_or_reject('reduced', name, path, interactions)

    def test_squashed(self):
        for name, path, interactions in CELLS:
            with self.subTest(interaction=name):
                self._check_or_reject('squashed', name, path, interactions)

    def test_general(self):
        for name, path, interactions in CELLS:
            with self.subTest(interaction=name):
                self._check_cell('general', name, path, interactions)


if __name__ == "__main__":
    unittest.main()
