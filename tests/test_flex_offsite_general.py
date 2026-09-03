#!/usr/bin/env python3
"""Off-site interactions in the general (full-vertex) FLEX: the measured policy.

Accepted off-site (#181, Tier 1): CoulombInter, Hund and Ising -- same-orbital
or inter-orbital, with or without sublattice folding. Each enters as its
Hartree vertex V_ab(q) on the density (aa,bb) slots ONLY (the locality-split
S/C builder, hwave.solver._sc_matrices_myo); the exchange crossing of an
off-site term is a non-local particle-hole pair, not representable by a
q-only matrix, and is absent -- exactly the RPA ring's reading, so for every
accepted class FLEX at IterationMax=1 is element-complete equal to the ring:
every spin-orbital component of chiq at every frequency and q, reconstructed
from the FLEX channels, to <= 1e-13. Locality is judged on the PRE-fold
declarations, so the answer does not depend on SubShape. The solver warns
once per solve that only the Hartree vertex of off-site terms enters
(bond-resolved treatment: #181 Tier 3).

Before the split the shared builder wrote off-site content into the cross
(ab,ab) slots and its on-site `l1 != l3` density gate deleted the same-orbital
off-site Hund/Ising: 3.0e-1 / 4.6e-1 against the ring at one orbital,
2.4e-2 .. 3.4e-2 for the inter-orbital classes, 1.3e-1 under folding.

Still rejected: off-site Exchange (its local regrouping -J S+_i S-_j is
transverse; the longitudinal content is unadjudicated -- #181 Tier 2; the ring
returns the bare bubble for it, so "equal to the ring" would be agreement at
zero vertex effect) and off-site PairHop (no local-bilinear regrouping at
all; the ring drops it, #157). An off-site CoulombIntra row is refused by the
reader (#93).
"""

import logging
import os
import unittest

import numpy as np


def _run_pair(path, interactions, cell, filling, flex_iters=1,
              inject=None, sub=(1, 1, 1)):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as rpa_mod
    import hwave.solver.flex as flex_mod

    par = {'T': 2.0, 'filling': filling, 'CellShape': cell,
           'SubShape': list(sub), 'Nmat': 32}

    def io():
        idict = {'path_to_input': path, 'Geometry': 'geom.dat',
                 'Transfer': 'transfer.dat'}
        idict.update(interactions)
        r = read_input_k.QLMSkInput({'path_to_input': path,
                                     'interaction': idict})
        ham = r.get_param("ham")
        if inject:
            for itype, key, val in inject:
                ham[itype][key] = val
        return r, ham

    os.makedirs('tests/rpa/output', exist_ok=True)
    i1, ham1 = io()
    rpa = rpa_mod.RPA(ham1, {},
                      {'mode': 'RPA', 'param': dict(par),
                       'enable_spin_orbital': False,
                       'calc_scheme': 'general', 'calc_type': 'ring'})
    gr = i1.get_param("green")
    rpa.solve(gr, 'tests/rpa/output')

    i2, ham2 = io()
    pf = dict(par)
    pf.update({'IterationMax': flex_iters, 'Mix': 1.0, 'EPS': 1})
    fx = flex_mod.FLEX(ham2, {},
                       {'mode': 'FLEX', 'param': pf,
                        'enable_spin_orbital': False,
                        'calc_scheme': 'general'})
    gf = i2.get_param("green")
    fx.solve(gf, 'tests/rpa/output')
    return gr, gf


def _assert_element_complete_equal(test, gr, gf, norb, atol=1e-12):
    """Reconstruct RPA's full spin-orbital chiq from the FLEX channels and
    compare EVERY element -- including the spin-mixed pair blocks, which must
    vanish for the reconstruction to be complete."""
    np.testing.assert_allclose(np.asarray(gf['chi0q']),
                               np.asarray(gr['chi0q']),
                               rtol=0.0, atol=atol)
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
    np.testing.assert_allclose(chiq, recon, rtol=0.0, atol=atol)


class TestOffsiteGeneralFLEX(unittest.TestCase):

    def test_one_orbital_offsite_v_matches_the_rpa_ring_exactly(self):
        """One orbital, off-site V, element-complete, on 4x4 and non-cubic
        4x6. (The in-loop chemical potential was also checked bit-for-bit
        identical between the two solvers; the `mu` attribute FLEX exposes
        afterwards belongs to its post-loop dressed bookkeeping and has no
        RPA counterpart.)"""
        for cell in ([4, 4, 1], [4, 6, 1]):
            with self.subTest(cell=cell):
                gr, gf = _run_pair('tests/rpa/input',
                                   {'CoulombInter': 'coulombinter.dat'},
                                   cell, filling=0.75)
                _assert_element_complete_equal(self, gr, gf, norb=1)

    def test_one_orbital_offsite_v_3d_with_a_z_bond(self):
        """A bond along z on a 4x4x2 lattice pins the third axis of the grid
        contract, which no committed fixture exercises."""
        gr, gf = _run_pair(
            'tests/rpa/input', {'CoulombInter': 'coulombinter.dat'},
            [4, 4, 2], filling=0.75,
            inject=[('CoulombInter', ((0, 0, 1), (0, 0)), 0.3),
                    ('CoulombInter', ((0, 0, -1), (0, 0)), 0.3)])
        _assert_element_complete_equal(self, gr, gf, norb=1)

    def test_two_orbital_sameorb_offsite_v_matches_exactly(self):
        """TWO orbitals with same-orbital (a == b) off-site bonds: the
        accepted class must agree with the ring in EVERY element -- there is
        no masked region. This is what makes the acceptance class sharp: any
        Fierz-slot contamination from the off-site entries would fail here,
        where the earlier confinement-masked test could not see it."""
        for cell in ([4, 4, 1], [4, 6, 1]):
            with self.subTest(cell=cell):
                gr, gf = _run_pair('tests/rpa/input_2orb',
                                   {'CoulombInter': 'offsite_sameorb.dat'},
                                   cell, filling=0.5)
                _assert_element_complete_equal(self, gr, gf, norb=2)

    def test_one_orbital_offsite_hund_and_ising_match_the_ring(self):
        """Same-orbital off-site Hund / Ising (every one-orbital model): a
        physical density-density term on the (00,00) slot. The shared
        builder's on-site `l1 != l3` gate used to delete it (measured
        3.0e-1 / 4.6e-1 against the ring); the locality split keeps it."""
        for itype in ('Hund', 'Ising'):
            with self.subTest(itype=itype):
                gr, gf = _run_pair('tests/rpa/input',
                                   {itype: 'coulombinter.dat'},
                                   [4, 4, 1], filling=0.75)
                _assert_element_complete_equal(self, gr, gf, norb=1)
                # anti-vacuity: the ring's vertex has a real effect here
                chiq = np.asarray(gr['chiq'])
                self.assertGreater(
                    np.max(np.abs(chiq[:, :, :1, :1, :1, :1]
                                  - np.asarray(gr['chi0q']))), 1e-3)

    def test_two_orbital_offsite_interorbital_classes_match_the_ring(self):
        """Off-site inter-orbital CoulombInter / Hund / Ising: the density
        (Hartree) half is q-only and ring-identical once the off-site
        content stays out of the cross (ab,ab) slots, whose particle-hole
        pair is non-local for R != 0 (measured 2.4e-2..3.4e-2 with the
        on-site slot map, <= 7e-15 with the locality split). The 2-orbital
        coulombinter.dat also carries ON-site a != b entries, so it pins
        that the split keeps the on-site Fierz content intact."""
        cases = [
            ('tests/rpa/input_2orb', {'CoulombInter': 'coulombinter.dat'}),
            ('tests/equivalence_input/orb2', {'Hund': 'offsite_hund.dat'}),
            ('tests/equivalence_input/orb2', {'Ising': 'offsite_ising.dat'}),
            ('tests/equivalence_input/orb2',
             {'CoulombInter': 'offsite_coulombinter_interorb.dat'}),
        ]
        for path, interactions in cases:
            with self.subTest(interaction=list(interactions)[0], path=path):
                gr, gf = _run_pair(path, interactions, [4, 4, 1],
                                   filling=0.5)
                _assert_element_complete_equal(self, gr, gf, norb=2)
                chiq = np.asarray(gr['chiq'])
                self.assertGreater(
                    np.max(np.abs(chiq[:, :, :2, :2, :2, :2]
                                  - np.asarray(gr['chi0q']))), 1e-3)

    def test_offsite_under_sublattice_folding_matches_the_ring(self):
        """Locality is judged on the PRE-fold declarations (the RPA ring's
        reading): a same-orbital +-x bond folded with SubShape=[2,1,1]
        becomes an intra-supercell inter-orbital entry, which the folded
        table cannot tell from on-site input. Reading it as on-site put
        the bond's Fock crossing into the cross slot and the solvers
        differed by 1.3e-1; the answer must not depend on SubShape.
        SubShape=[4,1,1] maps EVERY x displacement to (0,0,0)."""
        cases = [
            ('tests/rpa/input_2orb', {'CoulombInter': 'offsite_sameorb.dat'}),
            ('tests/equivalence_input/orb2', {'Hund': 'offsite_hund.dat'}),
            ('tests/equivalence_input/orb2', {'Ising': 'offsite_ising.dat'}),
            # same-orbital off-site Hund/Ising folded (the one-orbital
            # fixture's +-x/+-y bonds read as Hund / Ising)
            ('tests/rpa/input', {'Hund': 'coulombinter.dat'}),
            ('tests/rpa/input', {'Ising': 'coulombinter.dat'}),
        ]
        for path, interactions in cases:
            norb_phys = 1 if path == 'tests/rpa/input' else 2
            for sub in ((2, 1, 1), (4, 1, 1)):
                with self.subTest(interaction=list(interactions)[0],
                                  path=path, sub=sub):
                    gr, gf = _run_pair(path, interactions, [4, 4, 1],
                                       filling=0.5, sub=sub)
                    norb_folded = int(np.asarray(gr['chi0q']).shape[2])
                    self.assertEqual(norb_folded, norb_phys * sub[0])
                    _assert_element_complete_equal(self, gr, gf,
                                                   norb=norb_folded)
                    chiq = np.asarray(gr['chiq'])
                    n = norb_folded
                    self.assertGreater(
                        np.max(np.abs(chiq[:, :, :n, :n, :n, :n]
                                      - np.asarray(gr['chi0q']))), 1e-3)

    def test_folded_key_collision_between_the_two_parts_is_summed(self):
        """A displacement that is a full lattice period (R = (4,0,0) on a
        4-cell axis; internal table -- the reader-side wrap is the
        (n-1,0,0) form tested elsewhere) is judged off-site by the
        irvec == 0 rule, exactly as the RPA ring judges it, yet it folds
        onto the SAME key as the on-site entry of the same orbital pair.
        The whole-table k-space array must carry the SUM of the two
        (folding the combined table), not whichever part was merged
        last; the ring is the reference, and the on-site entry must
        matter."""
        wrapped = [('CoulombInter', ((4, 0, 0), (0, 1)), 0.2),
                   ('CoulombInter', ((-4, 0, 0), (1, 0)), 0.2)]
        gr, gf = _run_pair('tests/rpa/input_2orb',
                           {'CoulombInter': 'onsite_inter.dat'},
                           [4, 4, 1], filling=0.5, sub=(2, 1, 1),
                           inject=wrapped)
        _assert_element_complete_equal(self, gr, gf, norb=4)
        # anti-vacuity: the on-site declaration changes the answer
        gr_off_only, _ = _run_pair(
            'tests/rpa/input_2orb', {'CoulombInter': 'onsite_inter.dat'},
            [4, 4, 1], filling=0.5, sub=(2, 1, 1),
            inject=wrapped + [('CoulombInter', ((0, 0, 0), (0, 1)), 0.0),
                              ('CoulombInter', ((0, 0, 0), (1, 0)), 0.0)])
        self.assertGreater(np.max(np.abs(np.asarray(gr['chiq'])
                                         - np.asarray(gr_off_only['chiq']))),
                           1e-3)

        # Hund / Ising: a full-period SAME-orbital entry folds onto the
        # (a,a) key, whose density content the shared builder's on-site
        # gate deletes; with the split it comes through the off-site part,
        # as the ring (no gate) carries it
        for itype in ('Hund', 'Ising'):
            with self.subTest(itype=itype):
                wrapped_same = [(itype, ((4, 0, 0), (0, 0)), 0.15),
                                (itype, ((-4, 0, 0), (0, 0)), 0.15)]
                gr, gf = _run_pair('tests/equivalence_input/orb2',
                                   {itype: 'offsite_%s.dat' % itype.lower()},
                                   [4, 4, 1], filling=0.5, sub=(2, 1, 1),
                                   inject=wrapped_same)
                _assert_element_complete_equal(self, gr, gf, norb=4)
                gr0, _ = _run_pair('tests/equivalence_input/orb2',
                                   {itype: 'offsite_%s.dat' % itype.lower()},
                                   [4, 4, 1], filling=0.5, sub=(2, 1, 1),
                                   inject=[(t, k, 0.0)
                                           for t, k, _ in wrapped_same])
                self.assertGreater(
                    np.max(np.abs(np.asarray(gr['chiq'])
                                  - np.asarray(gr0['chiq']))), 1e-3)

    def test_aggregate_coulomb_full_period_entry_under_folding_follows_the_ring(self):
        """Aggregate `Coulomb` with a same-orbital displacement equal to a
        full lattice period (internal table) under folding: the reader's
        folded table maps it to R = 0 same-orbital, which split_coulomb
        then classifies as CoulombIntra (the U slot); the pre-fold split
        classifies the same entry as off-site CoulombInter. The RPA ring
        reads the folded classification, so the whole-table array must
        be built from the folded-then-normalised table (element-complete
        parity), while the locality split still judges the pre-fold
        declarations."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        import hwave.solver.flex as flex_mod

        table = {((0, 0, 0), (0, 0)): 3.0, ((0, 0, 0), (1, 1)): 3.0,
                 ((0, 0, 0), (0, 1)): 0.5, ((0, 0, 0), (1, 0)): 0.5,
                 ((4, 0, 0), (0, 0)): 0.2, ((-4, 0, 0), (0, 0)): 0.2,
                 ((1, 0, 0), (1, 1)): 0.1, ((-1, 0, 0), (1, 1)): 0.1}
        par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
               'SubShape': [2, 1, 1], 'Nmat': 32}

        def io():
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat'}})
            ham = r.get_param("ham")
            ham['Coulomb'] = dict(table)
            return r, ham

        os.makedirs('tests/rpa/output', exist_ok=True)
        r1, ham1 = io()
        rpa = rpa_mod.RPA(ham1, {}, {'mode': 'RPA', 'param': dict(par),
                                     'enable_spin_orbital': False,
                                     'calc_scheme': 'general',
                                     'calc_type': 'ring'})
        gr = r1.get_param("green")
        rpa.solve(gr, 'tests/rpa/output')
        r2, ham2 = io()
        pf = dict(par)
        pf.update({'IterationMax': 1, 'Mix': 1.0, 'EPS': 1})
        fx = flex_mod.FLEX(ham2, {}, {'mode': 'FLEX', 'param': pf,
                                      'enable_spin_orbital': False,
                                      'calc_scheme': 'general'})
        gf = r2.get_param("green")
        fx.solve(gf, 'tests/rpa/output')
        _assert_element_complete_equal(self, gr, gf, norb=4)
        # anti-vacuity: the full-period entry changes the ring's answer
        r3, ham3 = io()
        ham3['Coulomb'][((4, 0, 0), (0, 0))] = 0.0
        ham3['Coulomb'][((-4, 0, 0), (0, 0))] = 0.0
        rpa3 = rpa_mod.RPA(ham3, {}, {'mode': 'RPA', 'param': dict(par),
                                      'enable_spin_orbital': False,
                                      'calc_scheme': 'general',
                                      'calc_type': 'ring'})
        gr3 = r3.get_param("green")
        rpa3.solve(gr3, 'tests/rpa/output')
        self.assertGreater(np.max(np.abs(np.asarray(gr['chiq'])
                                         - np.asarray(gr3['chiq']))), 1e-3)

    def test_aggregate_coulomb_under_folding_equals_the_explicit_split(self):
        """Aggregate `Coulomb` is normalised into CoulombIntra + CoulombInter
        BEFORE the locality split and the per-part folding; the result
        must be bit-identical to the explicit declaration of the same
        tables (nothing lost or double-counted by folding each part
        separately), and the off-site part must have an effect."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        onsite_u = {((0, 0, 0), (0, 0)): 3.0, ((0, 0, 0), (1, 1)): 3.0}
        bond = {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                ((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}

        def _chiq(tables):
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat'}})
            ham = r.get_param("ham")
            ham.update(tables)
            pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                  'SubShape': [2, 1, 1], 'Nmat': 32,
                  'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
            fx = flex_mod.FLEX(ham, {}, {'mode': 'FLEX', 'param': pf,
                                         'enable_spin_orbital': False,
                                         'calc_scheme': 'general'})
            os.makedirs('tests/rpa/output', exist_ok=True)
            gf = r.get_param("green")
            fx.solve(gf, 'tests/rpa/output')
            return np.asarray(gf['chiq_s']), np.asarray(gf['chiq_c'])

        agg = _chiq({'Coulomb': {**onsite_u, **bond}})
        explicit = _chiq({'CoulombIntra': onsite_u, 'CoulombInter': bond})
        onsite_only = _chiq({'CoulombIntra': onsite_u})
        for a, e in zip(agg, explicit):
            np.testing.assert_array_equal(a, e)
        self.assertGreater(np.max(np.abs(explicit[1] - onsite_only[1])), 1e-3)

    def test_offsite_input_warns_that_only_the_hartree_vertex_enters(self):
        """Every off-site two-body term enters the general path as its
        Hartree (density-slot) vertex only; the exchange crossing is not
        representable by a q-only matrix and is absent (the ring makes the
        same approximation; bond-resolved treatment is #181 Tier 3). The
        solver says so once per solve, naming the types; purely on-site
        input stays silent."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        def _solve(path, interactions, sub=(1, 1, 1)):
            idict = {'path_to_input': path, 'Geometry': 'geom.dat',
                     'Transfer': 'transfer.dat'}
            idict.update(interactions)
            r = read_input_k.QLMSkInput({'path_to_input': path,
                                         'interaction': idict})
            pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                  'SubShape': list(sub), 'Nmat': 32,
                  'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
            fx = flex_mod.FLEX(r.get_param("ham"), {},
                               {'mode': 'FLEX', 'param': pf,
                                'enable_spin_orbital': False,
                                'calc_scheme': 'general'})
            gf = r.get_param("green")
            fx.solve(gf, 'tests/rpa/output')
            fx._test_green_info = gf
            return fx

        os.makedirs('tests/rpa/output', exist_ok=True)
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            fx = _solve('tests/rpa/input_2orb',
                        {'CoulombInter': 'coulombinter.dat',
                         'Hund': 'offsite_hund.dat'})
        hits = [m for m in cm.output if 'exchange' in m.lower()
                and 'off-site' in m.lower()]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertIn('CoulombInter', hits[0])
        self.assertIn('Hund', hits[0])

        # once per SOLVE, not once per solver instance: a second solve on
        # the same object (the S/C matrices are cached across the SCF
        # iterations of one solve) must say it again
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            fx.solve(fx._test_green_info, 'tests/rpa/output')
        self.assertEqual(
            len([m for m in cm.output if 'off-site' in m.lower()]), 1,
            cm.output)

        # folded: the warning names the PRE-fold declarations' types
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            _solve('tests/rpa/input_2orb',
                   {'CoulombInter': 'offsite_sameorb.dat'}, sub=(4, 1, 1))
        self.assertEqual(
            len([m for m in cm.output if 'off-site' in m.lower()]), 1,
            cm.output)

        # on-site only: no such warning (assertNoLogs is 3.10+; count)
        logger = logging.getLogger('hwave.solver.flex')
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        logger.addHandler(handler)
        try:
            _solve('tests/rpa/input_2orb', {'CoulombInter': 'onsite_inter.dat'})
        finally:
            logger.removeHandler(handler)
        self.assertFalse([r for r in records
                          if 'off-site' in r.getMessage().lower()])

    def test_wrapped_declaration_reads_like_its_signed_form(self):
        """A -x bond may be declared as (n-1, 0, 0) on an n-cell lattice.
        The symmetrisation reverses displacements modulo the grid (roll+flip,
        as uhfk.py does), so the wrapped and signed forms of the same bond
        must give bit-identical chiq; a sign-flipped key lookup would miss
        the wrapped partner and silently halve the coefficient."""
        signed = []
        wrapped = [('CoulombInter', ((-1, 0, 0), (0, 0)), 0.0),
                   ('CoulombInter', ((3, 0, 0), (0, 0)), 1.0)]
        gr1, _ = _run_pair('tests/rpa/input',
                           {'CoulombInter': 'coulombinter.dat'},
                           [4, 4, 1], filling=0.75, inject=signed)
        gr2, _ = _run_pair('tests/rpa/input',
                           {'CoulombInter': 'coulombinter.dat'},
                           [4, 4, 1], filling=0.75, inject=wrapped)
        np.testing.assert_array_equal(np.asarray(gr1['chiq']),
                                      np.asarray(gr2['chiq']))

    def test_aggregate_coulomb_equals_the_explicit_declaration(self):
        """An aggregate 'Coulomb' table means CoulombIntra (r = 0 diagonal)
        plus CoulombInter (everything else). The general path used to drop
        the key silently -- zero vertex, chiq_s off by 1e-1 -- because
        neither the guard nor the k-space builder knew the aggregate name.
        Now it must be bit-identical to the same input declared explicitly."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        def run(aggregate):
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat',
                                 'CoulombIntra': 'coulombintra.dat',
                                 'CoulombInter': 'onsite_inter.dat'}})
            ham = r.get_param("ham")
            if aggregate:
                agg = {}
                agg.update(ham.pop('CoulombIntra'))
                agg.update(ham.pop('CoulombInter'))
                ham['Coulomb'] = agg
            pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                  'SubShape': [1, 1, 1], 'Nmat': 32,
                  'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
            fx = flex_mod.FLEX(ham, {}, {'mode': 'FLEX', 'param': pf,
                                         'enable_spin_orbital': False,
                                         'calc_scheme': 'general'})
            gf = r.get_param("green")
            fx.solve(gf, 'tests/rpa/output')
            return gf

        os.makedirs('tests/rpa/output', exist_ok=True)
        ge, ga = run(False), run(True)
        for key in ('chiq_s', 'chiq_c'):
            np.testing.assert_array_equal(np.asarray(ge[key]),
                                          np.asarray(ga[key]))

    def test_aggregate_coulomb_offsite_follows_the_same_policy(self):
        """Off-site entries arriving through the aggregate table must take
        the same route as explicit CoulombInter: an a != b off-site bond
        declared as `Coulomb` gives chiq_s bit-identical to the explicit
        CoulombInter declaration of the same bond (and the bond must have
        an effect, so the comparison cannot pass vacuously)."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        def _chiq_s(key, table):
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat'}})
            ham = r.get_param("ham")
            if table is not None:
                ham[key] = table
            pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                  'SubShape': [1, 1, 1], 'Nmat': 32,
                  'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
            fx = flex_mod.FLEX(ham, {}, {'mode': 'FLEX', 'param': pf,
                                         'enable_spin_orbital': False,
                                         'calc_scheme': 'general'})
            os.makedirs('tests/rpa/output', exist_ok=True)
            gf = r.get_param("green")
            fx.solve(gf, 'tests/rpa/output')
            # an off-site CoulombInter has charge-channel content only
            # (density (S, C) = (0, +2)), so the effect check is on chiq_c
            return np.asarray(gf['chiq_s']), np.asarray(gf['chiq_c'])

        bond = {((1, 0, 0), (0, 1)): 0.7, ((-1, 0, 0), (1, 0)): 0.7}
        via_coulomb = _chiq_s('Coulomb', bond)
        via_inter = _chiq_s('CoulombInter', bond)
        none = _chiq_s('CoulombInter', None)
        np.testing.assert_array_equal(via_coulomb[0], via_inter[0])
        np.testing.assert_array_equal(via_coulomb[1], via_inter[1])
        self.assertGreater(np.max(np.abs(via_inter[1] - none[1])), 1e-3)

    def test_missing_prefold_table_fails_closed(self):
        """If folding is active but the pre-fold table is unavailable, the
        guard must refuse rather than fall back to the folded table (the
        folded table is exactly what the bypass exploited)."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        r = read_input_k.QLMSkInput(
            {'path_to_input': 'tests/rpa/input_2orb',
             'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                             'Geometry': 'geom.dat',
                             'Transfer': 'transfer.dat',
                             'CoulombInter': 'offsite_sameorb.dat'}})
        pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
              'SubShape': [2, 1, 1], 'Nmat': 32,
              'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
        fx = flex_mod.FLEX(r.get_param("ham"), {},
                           {'mode': 'FLEX', 'param': pf,
                            'enable_spin_orbital': False,
                            'calc_scheme': 'general'})
        if hasattr(fx.ham_info, 'param_ham_orig'):
            del fx.ham_info.param_ham_orig
        os.makedirs('tests/rpa/output', exist_ok=True)
        with self.assertRaises(ValueError) as cm:
            fx.solve(r.get_param("green"), 'tests/rpa/output')
        self.assertIn("param_ham_orig", str(cm.exception))

    def test_one_sided_onsite_hund_reads_as_the_mean(self):
        """The ring's declaration symmetrisation applies to every type, not
        just CoulombInter: a one-sided on-site Hund entry injected into the
        internal TABLE after the read (file input must be Hermitian-closed
        since #93) must give bit-identical chiq to the halved both-ends
        table."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        def run(entries):
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat'}})
            ham = r.get_param("ham")
            ham['Hund'] = dict(entries)
            par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                   'SubShape': [1, 1, 1], 'Nmat': 32}
            rpa = rpa_mod.RPA(ham, {}, {'mode': 'RPA', 'param': par,
                                        'enable_spin_orbital': False,
                                        'calc_scheme': 'general',
                                        'calc_type': 'ring'})
            gr = r.get_param("green")
            rpa.solve(gr, 'tests/rpa/output')
            return np.asarray(gr['chiq'])

        os.makedirs('tests/rpa/output', exist_ok=True)
        one = run({((0, 0, 0), (0, 1)): 0.6})
        both = run({((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3})
        np.testing.assert_array_equal(one, both)

    def test_rejected_offsite_classes(self):
        """Everything outside the measured-equal class must be rejected."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        cases = [
            # off-site Exchange: no adjudicated longitudinal S/C content
            # (its local regrouping -J S+_i S-_j is transverse; #181 Tier 2)
            ({'Exchange': 'coulombinter.dat'}, 'tests/rpa/input', [1, 1, 1]),
            ({'Exchange': 'offsite_exchange.dat'},
             'tests/equivalence_input/orb2', [1, 1, 1]),
            # off-site PairHop: no local-bilinear particle-hole regrouping
            # exists at all (RPA drops it too, #157)
            ({'PairHop': 'offsite_pairhop.dat'},
             'tests/equivalence_input/orb2', [1, 1, 1]),
            # CoulombIntra is on-site same-orbital by definition; an
            # off-site row under that key is refused by the READER (#93)
            ({'CoulombIntra': 'coulombinter.dat'}, 'tests/rpa/input',
             [1, 1, 1]),
            # folding must not hide a rejected class: with SubShape=[4,1,1]
            # every +-x displacement canonicalizes to (0,0,0) in the folded
            # table, so the guard must scan the PRE-fold table -- before
            # that fix folded off-site input was accepted and diverged
            ({'Exchange': 'offsite_exchange.dat'},
             'tests/equivalence_input/orb2', [4, 1, 1]),
        ]
        os.makedirs('tests/rpa/output', exist_ok=True)
        for interactions, path, sub in cases:
            with self.subTest(interactions=list(interactions)[0], sub=sub):
                idict = {'path_to_input': path, 'Geometry': 'geom.dat',
                         'Transfer': 'transfer.dat'}
                idict.update(interactions)

                def _construct_and_solve():
                    # the rejection may now fire at READ time (issue #93:
                    # e.g. an off-site file configured as CoulombIntra
                    # violates its on-site same-orbital rule) or at the
                    # FLEX off-site guard; both are loud ValueErrors and
                    # both verdicts are what this test pins
                    r = read_input_k.QLMSkInput({'path_to_input': path,
                                                 'interaction': idict})
                    pf = {'T': 2.0, 'filling': 0.5,
                          'CellShape': [4, 4, 1],
                          'SubShape': sub, 'Nmat': 32,
                          'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
                    fx = flex_mod.FLEX(r.get_param("ham"), {},
                                       {'mode': 'FLEX', 'param': pf,
                                        'enable_spin_orbital': False,
                                        'calc_scheme': 'general'})
                    fx.solve(r.get_param("green"), 'tests/rpa/output')

                with self.assertRaises(ValueError):
                    _construct_and_solve()

    def test_one_sided_offsite_declaration_matches_the_ring_exactly(self):
        """A one-sided internal TABLE (only +R; built after the read, since
        file input must be closed per #93) means the reversal-symmetric
        operator v n(i) n(i+R), whose exact vertex is v cos(qR). Both solvers
        now read it that way -- the S/C builders since #114, the ring since
        rpa.py got the same symmetrisation (before that the ring kept a
        one-sided v e^{-iqR}, off by 1.2e-2 in chiq, the #106 discrepancy) --
        so the element-complete equivalence must hold for one-sided input,
        and it must equal the v/2-both-ends reading bit for bit."""
        one_sided = [('CoulombInter', ((-1, 0, 0), (0, 0)), 0.0),
                     ('CoulombInter', ((0, -1, 0), (0, 0)), 0.0)]
        gr, gf = _run_pair('tests/rpa/input',
                           {'CoulombInter': 'coulombinter.dat'},
                           [4, 4, 1], filling=0.75, inject=one_sided)
        _assert_element_complete_equal(self, gr, gf, norb=1)
        halved = [('CoulombInter', ((s * d[0], s * d[1], 0), (0, 0)), 0.5)
                  for d in ((1, 0), (0, 1)) for s in (1, -1)]
        gr2, _ = _run_pair('tests/rpa/input',
                           {'CoulombInter': 'coulombinter.dat'},
                           [4, 4, 1], filling=0.75, inject=halved)
        np.testing.assert_array_equal(np.asarray(gr['chiq']),
                                      np.asarray(gr2['chiq']))

    def test_multi_iteration_offsite_v_runs_and_stays_finite(self):
        """Three iterations exercise the q-dependent V_eff through the
        self-energy path; the outputs must stay finite."""
        _, gf = _run_pair('tests/rpa/input',
                          {'CoulombInter': 'coulombinter.dat'},
                          [4, 4, 1], filling=0.75, flex_iters=3)
        for key in ('sigma', 'chiq_s', 'chiq_c'):
            self.assertTrue(np.all(np.isfinite(np.asarray(gf[key]))),
                            "%s must be finite" % key)




class TestAggregateCoulombPreservesKeyCase(unittest.TestCase):
    """An interaction declared with a non-canonical key must still be
    applied when the aggregate ``Coulomb`` input is used.

    FLEX's general path normalises the aggregate ``Coulomb`` table into
    ``CoulombIntra``/``CoulombInter`` before scanning it. That rebuild used
    a dict comprehension, which returned a plain dict and so discarded the
    reader's ``CaseInsensitiveDict``. The reader stores every table under
    the spelling the USER wrote, so afterwards a ``hund`` key was invisible
    to both the off-site guard and the interaction builder and the term was
    dropped with no error, no warning, and a plausible-looking result.

    The comparison below would pass vacuously if the interaction were
    dropped in BOTH runs, so it also requires the term to CHANGE the
    result: canonical-vs-absent must differ.
    """

    def _run(self, hund_key, with_hund=True, hund_file="hund_onsite.dat",
             file_dir=None):
        import shutil
        import tempfile

        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        src = "tests/rpa/input_2orb"
        d = tempfile.mkdtemp(prefix="flex_keycase_")
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for f in ("geom.dat", "transfer.dat"):
            shutil.copy(os.path.join(src, f), d)
        if file_dir is not None:
            src = file_dir
        # Aggregate Coulomb: r=0 orbital-diagonal entries only, so
        # split_coulomb sends all of it to CoulombIntra.
        with open(os.path.join(d, "coulomb.dat"), "w") as f:
            f.write("aggregate Coulomb\n2\n1\n 1\n"
                    "   0    0    0    1    1   4.000000000000"
                    "   0.000000000000\n"
                    "   0    0    0    2    2   4.000000000000"
                    "   0.000000000000\n")
        shutil.copy(os.path.join(src, hund_file), d)

        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "Coulomb": "coulomb.dat"}
        if with_hund:
            inter[hund_key] = hund_file
        io = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
               "SubShape": [1, 1, 1], "Nmat": 16,
               "IterationMax": 1, "Mix": 1.0, "EPS": 1}
        solver = flex_mod.FLEX(
            io.get_param("ham"), {},
            {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
             "calc_scheme": "general"})
        green_info = io.get_param("green")
        out = tempfile.mkdtemp(prefix="flex_keycase_out_")
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        solver.solve(green_info, out)
        return np.asarray(green_info["chiq_s"])

    def test_non_canonical_key_is_applied_like_the_canonical_one(self):
        absent = self._run("Hund", with_hund=False)
        canonical = self._run("Hund")
        lowercase = self._run("hund")

        # Non-vacuity first: the term must actually do something, or the
        # equality below would hold simply because it is dropped twice.
        self.assertGreater(
            float(np.max(np.abs(canonical - absent))), 1e-6,
            "the fixture's Hund term has no effect, so this test could "
            "not tell an applied interaction from a dropped one")

        self.assertTrue(
            np.array_equal(lowercase, canonical),
            "a lowercase 'hund' key gives a different result from 'Hund' "
            "when the aggregate Coulomb input is used: max|diff| = {}"
            .format(float(np.max(np.abs(lowercase - canonical)))))

    def test_non_canonical_key_still_reaches_the_offsite_guard(self):
        """The same container loss hid entries from the off-site GUARD as
        well as from the interaction builder, and an unsupported off-site
        declaration that the guard cannot see completes silently instead
        of raising. The test above exercises the builder with an on-site
        term; this one exercises the off-site route both ways: a
        supported off-site type (Hund) must be APPLIED under either
        spelling (identically, and with an effect), and a rejected one
        (Exchange) must RAISE under either spelling.
        """
        canonical = self._run("Hund", hund_file="offsite_hund.dat")
        lower = self._run("hund", hund_file="offsite_hund.dat")
        absent = self._run("Hund", with_hund=False)
        np.testing.assert_array_equal(canonical, lower)
        self.assertGreater(np.max(np.abs(canonical - absent)), 1e-3)

        for key in ("Exchange", "exchange"):
            with self.subTest(key=key):
                with self.assertRaises(ValueError) as cm:
                    self._run(key, hund_file="offsite_exchange.dat",
                              file_dir="tests/equivalence_input/orb2")
                self.assertIn("off-site", str(cm.exception))

if __name__ == "__main__":
    unittest.main()
