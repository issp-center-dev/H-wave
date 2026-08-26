#!/usr/bin/env python3
"""Off-site interactions in the general (full-vertex) FLEX: the measured policy.

Accepted off-site: CoulombInter with SAME-orbital pairs (a == b), without
sublattice folding. For exactly that class, FLEX at IterationMax=1 is
element-complete equal to the RPA ring -- every spin-orbital component of chiq
at every frequency and q, reconstructed from the FLEX channels, to ~1e-14 --
because the vertex is V(q) on the density slots alone and the S/C k-grid is
the same C-ordered 2*pi*i/n grid as chi0's FFT axis.

Everything else stays rejected, each for a measured reason:

* CoulombInter with a != b off-site, and Hund / Ising off-site: the MYO S/C
  builder applies the full on-site Kanamori slot mapping, which feeds the
  Fierz (Case 2) inter-orbital slots; off-site, the particle-hole pair behind
  those slots is non-local and not representable by a q-only matrix. Off-site
  Hund / Ising differ from the ring by 3e-2 / 7e-2 even at ONE orbital, where
  no inter-orbital slot exists to blame.
* Off-site combined with sublattice folding: folding turns part of an a == b
  bond into intra-cell inter-orbital coupling and the solvers then differ by
  2e-2 (the folded analogue of the #104 content); deferred to #107.
* Exchange / PairHop off-site (non-local pair), CoulombIntra off-site (#106).
"""

import os
import unittest

import numpy as np


def _run_pair(path, interactions, cell, filling, flex_iters=1,
              inject=None):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as rpa_mod
    import hwave.solver.flex as flex_mod

    par = {'T': 2.0, 'filling': filling, 'CellShape': cell,
           'SubShape': [1, 1, 1], 'Nmat': 32}

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
        """Off-site entries arriving through the aggregate table must hit the
        same guard as explicit CoulombInter: a != b off-site is rejected."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        r = read_input_k.QLMSkInput(
            {'path_to_input': 'tests/rpa/input_2orb',
             'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                             'Geometry': 'geom.dat',
                             'Transfer': 'transfer.dat'}})
        ham = r.get_param("ham")
        ham['Coulomb'] = {((1, 0, 0), (0, 1)): 0.7,
                          ((-1, 0, 0), (1, 0)): 0.7}
        pf = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
              'SubShape': [1, 1, 1], 'Nmat': 32,
              'IterationMax': 1, 'Mix': 1.0, 'EPS': 1}
        fx = flex_mod.FLEX(ham, {}, {'mode': 'FLEX', 'param': pf,
                                     'enable_spin_orbital': False,
                                     'calc_scheme': 'general'})
        os.makedirs('tests/rpa/output', exist_ok=True)
        with self.assertRaises(ValueError):
            fx.solve(r.get_param("green"), 'tests/rpa/output')

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
            # a != b off-site CoulombInter (the 2orb fixture has such bonds)
            ({'CoulombInter': 'coulombinter.dat'}, 'tests/rpa/input_2orb',
             [1, 1, 1]),
            ({'Hund': 'coulombinter.dat'}, 'tests/rpa/input', [1, 1, 1]),
            ({'Ising': 'coulombinter.dat'}, 'tests/rpa/input', [1, 1, 1]),
            ({'Exchange': 'coulombinter.dat'}, 'tests/rpa/input', [1, 1, 1]),
            ({'CoulombIntra': 'coulombinter.dat'}, 'tests/rpa/input',
             [1, 1, 1]),
            # the accepted class, but with sublattice folding
            ({'CoulombInter': 'offsite_sameorb.dat'}, 'tests/rpa/input_2orb',
             [2, 1, 1]),
            # folding that maps the bond direction to a single supercell:
            # every off-site displacement canonicalizes to (0,0,0) in the
            # folded table, so the guard must scan the PRE-fold table --
            # before that fix this input was accepted and diverged
            ({'CoulombInter': 'offsite_sameorb.dat'}, 'tests/rpa/input_2orb',
             [4, 1, 1]),
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

    def _run(self, hund_key, with_hund=True, hund_file="hund_onsite.dat"):
        import shutil
        import tempfile

        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod

        src = "tests/rpa/input_2orb"
        d = tempfile.mkdtemp(prefix="flex_keycase_")
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        for f in ("geom.dat", "transfer.dat"):
            shutil.copy(os.path.join(src, f), d)
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
        of raising. The test above exercises the builder (its Hund term is
        on-site); this one exercises the guard.
        """
        with self.assertRaises(ValueError) as cm:
            self._run("Hund", hund_file="offsite_hund.dat")
        self.assertIn("off-site", str(cm.exception))

        with self.assertRaises(ValueError) as cm:
            self._run("hund", hund_file="offsite_hund.dat")
        self.assertIn("off-site", str(cm.exception))

if __name__ == "__main__":
    unittest.main()
