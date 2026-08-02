#!/usr/bin/env python3
"""RPA == FLEX at one shot, across calc_schemes (#107 phase 1, #104 fix).

RPA is the one-shot limit of FLEX (one iteration, no self-energy feedback,
bare Green function). These tests pin the measured equivalence matrix:

* reduced / squashed: equal to ~1e-13 in every interaction cell (measured
  before the #104 fix and unchanged by it -- the Fierz slots sit off the
  density diagonal those schemes read).
* general: the multi-orbital on-site cells (inter-orbital CoulombInter,
  Hund, combined U+V+J) disagreed by up to 1.2e-2 while #104 was open --
  the ring's longitudinal W lacked the adjudicated Fierz cross-slot
  content -- and must now be element-complete equal.
* spin modes: the corrected W is shared by every spin mode, so a
  spin-symmetric system solved in spin-diag mode (two identical bubble
  blocks) must reproduce the spin-free result bit for bit.

The spin/charge mapping (chi_s = chi_uu - chi_ud, chi_c = chi_uu + chi_ud on
the spin blocks of RPA's chiq) was verified on the CoulombIntra control,
where the two solvers agreed both before and after the fix.
"""

import os
import unittest

import numpy as np


def _read(path, interactions):
    import hwave.qlmsio.read_input_k as read_input_k

    idict = {'path_to_input': path, 'Geometry': 'geom.dat',
             'Transfer': 'transfer.dat'}
    idict.update(interactions)
    return read_input_k.QLMSkInput({'path_to_input': path,
                                    'interaction': idict})


def _run(scheme, solver, path, interactions, filling):
    import hwave.solver.rpa as rpa_mod
    import hwave.solver.flex as flex_mod

    r = _read(path, interactions)
    ham = r.get_param("ham")
    par = {'T': 2.0, 'filling': filling, 'CellShape': [4, 4, 1],
           'SubShape': [1, 1, 1], 'Nmat': 32}
    if solver == 'rpa':
        sv = rpa_mod.RPA(ham, {}, {'mode': 'RPA', 'param': par,
                                   'enable_spin_orbital': False,
                                   'calc_scheme': scheme,
                                   'calc_type': 'ring'})
    else:
        pf = dict(par)
        pf.update({'IterationMax': 1, 'Mix': 1.0, 'EPS': 1})
        sv = flex_mod.FLEX(ham, {}, {'mode': 'FLEX', 'param': pf,
                                     'enable_spin_orbital': False,
                                     'calc_scheme': scheme})
    g = r.get_param("green")
    os.makedirs('tests/rpa/output', exist_ok=True)
    sv.solve(g, 'tests/rpa/output')
    return g


class TestGeneralSchemeOneShot(unittest.TestCase):
    """The cells that carried the #104 deficiency, element-complete."""

    CASES = [
        ("inter-orbital V", {'CoulombInter': 'onsite_inter.dat'}),
        ("Hund", {'Hund': 'hund_onsite.dat'}),
        ("U+V+J", {'CoulombIntra': 'coulombintra.dat',
                   'CoulombInter': 'onsite_inter.dat',
                   'Hund': 'hund_onsite.dat'}),
    ]

    def test_multiorbital_onsite_cells_match_flex(self):
        norb = 2
        for name, inter in self.CASES:
            with self.subTest(interaction=name):
                gr = _run('general', 'rpa', 'tests/rpa/input_2orb',
                          inter, 0.5)
                gf = _run('general', 'flex', 'tests/rpa/input_2orb',
                          inter, 0.5)
                np.testing.assert_allclose(
                    np.asarray(gf['chi0q']), np.asarray(gr['chi0q']),
                    rtol=0.0, atol=1e-12)
                chiq = np.asarray(gr['chiq'])
                cs = np.asarray(gf['chiq_s'])
                cc = np.asarray(gf['chiq_c'])
                recon = np.zeros_like(chiq)
                same = 0.5 * (cc + cs)
                diff = 0.5 * (cc - cs)
                for s1 in (0, 1):
                    for s2 in (0, 1):
                        blk = same if s1 == s2 else diff
                        recon[:, :,
                              s1*norb:(s1+1)*norb, s1*norb:(s1+1)*norb,
                              s2*norb:(s2+1)*norb, s2*norb:(s2+1)*norb] = blk
                np.testing.assert_allclose(
                    chiq, recon, rtol=0.0, atol=1e-12,
                    err_msg="%s: the ring must reproduce the FLEX channel "
                            "solve element-complete (this cell was 1e-2 off "
                            "while #104 was open)" % name)


class TestReducedSquashedOneShot(unittest.TestCase):
    """reduced/squashed were equal before the #104 fix and must stay equal:
    the Fierz slots sit off the density diagonal these schemes read."""

    def _blocks(self, cq, norb):
        if cq.ndim == 8:
            # squashed layout: chiq[w, k, s1, s1', a, s2, s2', b]
            uu = cq[:, :, 0, 0, :, 0, 0, :]
            ud = cq[:, :, 0, 0, :, 1, 1, :]
        else:
            uu = cq[:, :, :norb, :norb]
            ud = cq[:, :, :norb, norb:]
        return uu, ud

    def test_matrix_cells(self):
        cases = [
            ('tests/rpa/input', {'CoulombInter': 'coulombinter.dat'},
             0.75, 1),
            ('tests/rpa/input_2orb', {'CoulombIntra': 'coulombintra.dat',
                                      'CoulombInter': 'onsite_inter.dat',
                                      'Hund': 'hund_onsite.dat'},
             0.5, 2),
        ]
        for scheme in ('reduced', 'squashed'):
            for path, inter, filling, norb in cases:
                with self.subTest(scheme=scheme, path=path):
                    gr = _run(scheme, 'rpa', path, inter, filling)
                    gf = _run(scheme, 'flex', path, inter, filling)
                    uu, ud = self._blocks(np.asarray(gr['chiq']), norb)
                    cs = np.asarray(gf['chiq_s'])[:, :, :norb, :norb]
                    cc = np.asarray(gf['chiq_c'])[:, :, :norb, :norb]
                    np.testing.assert_allclose(uu - ud, cs,
                                               rtol=0.0, atol=1e-12)
                    np.testing.assert_allclose(uu + ud, cc,
                                               rtol=0.0, atol=1e-12)


class TestFierzDomain(unittest.TestCase):
    """The correction applies exactly on the adjudicated domain: on-site
    inter-orbital entries of the PRE-fold table. Folding maps an off-site
    bond to (0,0,0) between supercell orbitals, so judging locality on the
    folded table smuggled unadjudicated off-site content into the
    correction (a one-orbital +-x bond folded with SubShape=[2,1,1]
    produced a spurious Fierz tensor of magnitude V)."""

    def _ham_info(self, norb, entries, sub, key='CoulombInter'):
        import hwave.solver.rpa as rpa_mod

        param_ham = {'Geometry': {'norb': norb, 'rvec': np.eye(3),
                                  'center': [[0, 0, 0]] * norb},
                     'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                     key: entries}
        info = {'mode': 'RPA', 'param': {'T': 2.0, 'filling': 0.5,
                'CellShape': [4, 1, 1], 'SubShape': list(sub), 'Nmat': 4},
                'enable_spin_orbital': False, 'calc_scheme': 'general',
                'calc_type': 'ring'}
        return rpa_mod.RPA(param_ham, {}, info).ham_info

    def test_folded_offsite_bond_is_a_fierz_noop(self):
        hi = self._ham_info(1, {((1, 0, 0), (0, 0)): 0.7,
                                ((-1, 0, 0), (0, 0)): 0.7}, [2, 1, 1])
        self.assertIsNone(hi.ham_fierz_q)

    def test_folded_onsite_interorbital_keeps_the_correction(self):
        hi = self._ham_info(2, {((0, 0, 0), (0, 1)): 0.7,
                                ((0, 0, 0), (1, 0)): 0.7}, [2, 1, 1])
        self.assertIsNotNone(hi.ham_fierz_q)
        self.assertAlmostEqual(float(np.max(np.abs(hi.ham_fierz_q))),
                               0.7, places=12)

    def test_intra_only_and_one_orbital_have_no_fierz(self):
        hi = self._ham_info(1, {((0, 0, 0), (0, 0)): 4.0}, [1, 1, 1],
                            key='CoulombIntra')
        self.assertIsNone(hi.ham_fierz_q)
        hi = self._ham_info(1, {((1, 0, 0), (0, 0)): 0.7,
                                ((-1, 0, 0), (0, 0)): 0.7}, [1, 1, 1])
        self.assertIsNone(hi.ham_fierz_q)

    def test_folded_offsite_aggregate_is_a_fierz_noop(self):
        # the aggregate path must also judge locality pre-fold: an off-site
        # same-orbital bond arriving through 'Coulomb' folds to (0,0,0)
        # between supercell orbitals just like the explicit table
        hi = self._ham_info(1, {((1, 0, 0), (0, 0)): 0.7,
                                ((-1, 0, 0), (0, 0)): 0.7}, [2, 1, 1],
                            key='Coulomb')
        self.assertIsNone(hi.ham_fierz_q)

    def test_fierz_support_sits_on_the_cross_slots_only(self):
        # exact support per type: on-site inter-orbital entries at 2
        # orbitals, no folding. nd = 4 (spin-major: 0=(up,0), 1=(up,1),
        # 2=(dn,0), 3=(dn,1)); the pair-diagonal entries (a,b,a,b) carry
        # the adjudicated coefficient -- same-spin for CoulombInter (-V),
        # Hund (+J) and Ising (-I), cross-spin for Exchange (+J') -- and
        # nothing else is populated
        J = 0.7
        up = {0: 0, 1: 1}
        dn = {0: 2, 1: 3}
        same = [(up, up), (dn, dn)]
        cross = [(up, dn), (dn, up)]
        table = {
            'CoulombInter': (-J, same),
            'Hund': (+J, same),
            'Ising': (-J, same),
            'Exchange': (+J, cross),
        }
        for key, (coeff, blocks) in table.items():
            with self.subTest(interaction=key):
                hi = self._ham_info(2, {((0, 0, 0), (0, 1)): J,
                                        ((0, 0, 0), (1, 0)): J},
                                    [1, 1, 1], key=key)
                f = np.asarray(hi.ham_fierz_q).reshape(4, 4, 4, 4, 4)[0]
                want = np.zeros((4, 4, 4, 4), dtype=complex)
                for s_bra, s_ket in blocks:
                    for a, b in ((0, 1), (1, 0)):
                        want[s_bra[a], s_bra[b],
                             s_ket[a], s_ket[b]] = coeff
                np.testing.assert_allclose(
                    f, want, atol=1e-13,
                    err_msg="%s Fierz support" % key)

    def test_aggregate_coulomb_fierz_equals_explicit(self):
        # under folding the aggregate path must split the PRE-fold table
        explicit = self._ham_info(2, {((0, 0, 0), (0, 1)): 0.7,
                                      ((0, 0, 0), (1, 0)): 0.7}, [2, 1, 1])
        aggregate = self._ham_info(
            2, {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7,
                ((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 2.0},
            [2, 1, 1], key='Coulomb')
        np.testing.assert_allclose(aggregate.ham_fierz_q,
                                   explicit.ham_fierz_q, atol=1e-13)


class TestSpinModeConsistency(unittest.TestCase):

    def test_spin_diag_matches_spin_free_on_symmetric_input(self):
        """A spin-symmetric system solved in spin-diag mode (two identical
        bubble blocks through the chi0q_init route) must reproduce the
        spin-free chiq bit for bit. The corrected longitudinal W is one
        tensor shared by every spin mode, so the #104 fix applies to
        spin-resolved runs identically -- this is what a channel-level
        patch of the spin-free path alone could not have guaranteed."""
        import hwave.solver.rpa as rpa_mod

        inter = {'CoulombInter': 'onsite_inter.dat'}
        gr = _run('general', 'rpa', 'tests/rpa/input_2orb', inter, 0.5)
        chi0q_free = np.asarray(gr['chi0q'])
        chiq_free = np.asarray(gr['chiq'])

        import tempfile

        out = tempfile.mkdtemp()
        np.savez(os.path.join(out, 'chi0q_diag_su2test.npz'),
                 chi0q=np.stack([chi0q_free, chi0q_free], axis=0), momentum_convention="e_plus_ikR")
        try:
            r = _read('tests/rpa/input_2orb', inter)
            par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                   'SubShape': [1, 1, 1], 'Nmat': 32}
            rpa = rpa_mod.RPA(r.get_param("ham"), {},
                              {'mode': 'RPA', 'param': par,
                               'enable_spin_orbital': False,
                               'calc_scheme': 'general',
                               'calc_type': 'ring'})
            g = r.get_param("green")
            g.update(rpa.read_init({'path_to_input': out,
                                    'chi0q_init': 'chi0q_diag_su2test.npz'}))
            rpa.solve(g, out)
            self.assertEqual(rpa.spin_mode, "spin-diag")
            np.testing.assert_array_equal(np.asarray(g['chiq']), chiq_free)
        finally:
            import shutil

            shutil.rmtree(out, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
