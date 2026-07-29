#!/usr/bin/env python3
"""Off-site density-density interactions in the general (full-vertex) FLEX.

The general path used to reject every off-site entry, citing a feared mismatch
between the S/C construction's k-grid and chi0's FFT q-grid. The fear was
unfounded -- both are the C-ordered 2*pi*i/n grid -- and the density-density
family (CoulombInter, Hund, Ising) needs nothing beyond the e^{-iqR} phases
that `_build_interaction_k` already carries: the vertex is V(q) at the density
slots, ordinary extended-Hubbard RPA at one iteration.

The clean grid test is ONE orbital: there are no inter-orbital pair slots
there, so FLEX at IterationMax=1 must reproduce the RPA ring exactly -- any
grid misalignment would appear directly. Exchange and PairHop stay rejected
off-site (their particle-hole pair is non-local, so no q-only vertex can
represent it), as does CoulombIntra (UHFk reads only its r = 0 component, so
accepting it here would create a cross-solver semantics divergence).
"""

import os
import unittest

import numpy as np


def _run_pair(path, interactions, cell, filling, flex_iters=1):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as rpa_mod
    import hwave.solver.flex as flex_mod

    par = {'T': 2.0, 'filling': filling, 'CellShape': cell,
           'SubShape': [1, 1, 1], 'Nmat': 32}

    def io():
        idict = {'path_to_input': path, 'Geometry': 'geom.dat',
                 'Transfer': 'transfer.dat'}
        idict.update(interactions)
        return read_input_k.QLMSkInput({'path_to_input': path,
                                        'interaction': idict})

    os.makedirs('tests/rpa/output', exist_ok=True)
    i1 = io()
    rpa = rpa_mod.RPA(i1.get_param("ham"), {},
                      {'mode': 'RPA', 'param': dict(par),
                       'enable_spin_orbital': False,
                       'calc_scheme': 'general', 'calc_type': 'ring'})
    gr = i1.get_param("green")
    rpa.solve(gr, 'tests/rpa/output')

    i2 = io()
    pf = dict(par)
    pf.update({'IterationMax': flex_iters, 'Mix': 1.0, 'EPS': 1})
    fx = flex_mod.FLEX(i2.get_param("ham"), {},
                       {'mode': 'FLEX', 'param': pf,
                        'enable_spin_orbital': False,
                        'calc_scheme': 'general'})
    gf = i2.get_param("green")
    fx.solve(gf, 'tests/rpa/output')
    return gr, gf


class TestOffsiteGeneralFLEX(unittest.TestCase):

    def test_one_orbital_offsite_v_matches_the_rpa_ring_exactly(self):
        """At one orbital there are no inter-orbital slots, so any difference
        from the ring would be a grid or phase error. Non-cubic included so
        both spatial axes are exercised."""
        for cell in ([4, 4, 1], [4, 6, 1]):
            with self.subTest(cell=cell):
                gr, gf = _run_pair('tests/rpa/input',
                                   {'CoulombInter': 'coulombinter.dat'},
                                   cell, filling=0.75)
                np.testing.assert_allclose(
                    np.asarray(gf['chi0q']), np.asarray(gr['chi0q']),
                    rtol=0.0, atol=1e-12)
                chiq = np.asarray(gr['chiq'])
                zz = chiq[:, :, 0, 0, 0, 0] - chiq[:, :, 0, 0, 1, 1]
                cc = chiq[:, :, 0, 0, 0, 0] + chiq[:, :, 0, 0, 1, 1]
                np.testing.assert_allclose(
                    np.asarray(gf['chiq_s'])[:, :, 0, 0, 0, 0], zz,
                    rtol=0.0, atol=1e-12)
                np.testing.assert_allclose(
                    np.asarray(gf['chiq_c'])[:, :, 0, 0, 0, 0], cc,
                    rtol=0.0, atol=1e-12)

    def test_two_orbital_offsite_v_differs_only_at_interorbital_slots(self):
        """chi0 must still be exact; the channel difference against the ring
        must be confined to the inter-orbital pair slots -- the same
        vertex-content difference as on-site U' (#104), where exact
        diagonalization already showed the FLEX side to be the correct one.
        Anything appearing elsewhere would be a genuine off-site defect."""
        gr, gf = _run_pair('tests/rpa/input_2orb',
                           {'CoulombInter': 'coulombinter.dat'},
                           [4, 4, 1], filling=0.5)
        np.testing.assert_allclose(
            np.asarray(gf['chi0q']), np.asarray(gr['chi0q']),
            rtol=0.0, atol=1e-12)
        norb, npair = 2, 4
        chiq = np.asarray(gr['chiq'])
        nw, nq = chiq.shape[0], chiq.shape[1]
        zz = np.zeros((nw, nq, npair, npair), dtype=complex)
        for a in range(norb):
            for c in range(norb):
                for b in range(norb):
                    for d in range(norb):
                        zz[:, :, a*norb+c, b*norb+d] = (
                            chiq[:, :, a, c, b, d]
                            - chiq[:, :, a, c, norb+b, norb+d])
        diff = np.abs(zz - np.asarray(gf['chiq_s']).reshape(
            nw, nq, npair, npair))
        mask = np.zeros((npair, npair), dtype=bool)
        for i in (1, 2):
            for j in (1, 2):
                mask[i, j] = True
        self.assertGreater(float(np.max(diff[:, :, mask])), 1e-3,
                           "the inter-orbital slots must differ while #104 "
                           "is open")
        self.assertLess(float(np.max(diff[:, :, ~mask])), 1e-3,
                        "any difference outside the inter-orbital slots "
                        "would be an off-site defect, not #104")

    def test_offsite_exchange_and_coulombintra_are_still_rejected(self):
        for itype in ("Exchange", "CoulombIntra"):
            with self.subTest(interaction=itype):
                with self.assertRaises(ValueError):
                    _run_pair('tests/rpa/input_2orb',
                              {itype: 'coulombinter.dat'},
                              [4, 4, 1], filling=0.5)

    def test_multi_iteration_offsite_v_runs_and_stays_finite(self):
        """The self-energy path receives a q-dependent V_eff once off-site
        entries are allowed; three iterations must run and produce finite
        sigma and susceptibilities."""
        _, gf = _run_pair('tests/rpa/input',
                          {'CoulombInter': 'coulombinter.dat'},
                          [4, 4, 1], filling=0.75, flex_iters=3)
        for key in ('sigma', 'chiq_s', 'chiq_c'):
            self.assertTrue(np.all(np.isfinite(np.asarray(gf[key]))),
                            "%s must be finite" % key)


if __name__ == "__main__":
    unittest.main()
