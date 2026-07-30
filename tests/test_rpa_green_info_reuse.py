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


class TestValidationStrictness(unittest.TestCase):

    def test_junk_shapes_are_rejected_even_under_optimized_python(self):
        """The validation must not be built on `assert`: under python -O a
        malformed array was silently accepted and assigned a spin mode
        (verified during review with shape (3, 99, 7, 8))."""
        import subprocess
        import sys
        import textwrap

        code = textwrap.dedent("""
            import sys, numpy as np
            sys.path.insert(0, %r)
            import hwave.solver.rpa as rpa_mod
            param_ham = {'Geometry': {'norb': 1, 'rvec': np.eye(3),
                                      'center': [[0, 0, 0]]},
                         'Transfer': {((0, 0, 0), (0, 0)): 1.0},
                         'CoulombIntra': {((0, 0, 0), (0, 0)): 4.0}}
            info = {'mode': 'RPA', 'param': {'T': 2.0, 'filling': 0.5,
                    'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                    'Nmat': 16},
                    'enable_spin_orbital': False, 'calc_scheme': 'reduced',
                    'calc_type': 'ring'}
            sv = rpa_mod.RPA(param_ham, {}, info)
            try:
                sv._validate_chi0q_shape(
                    np.zeros((3, 99, 7, 8), dtype=complex), source="probe")
                print("ACCEPTED")
            except ValueError:
                print("REJECTED")
        """) % (os.path.abspath("src"),)
        r = subprocess.run([sys.executable, "-O", "-c", code],
                           capture_output=True, text=True)
        self.assertIn("REJECTED", r.stdout,
                      "junk chi0q accepted under python -O: %s %s"
                      % (r.stdout, r.stderr[-300:]))

    def test_full_range_spin_diag_input_does_not_log_partial_range(self):
        """The partial-frequency check must read the frequency axis: a
        spin-diag array carries the spin-block axis first (shape[0] == 2),
        which used to be misreported as 'partial range ... 2 in 32'."""
        import logging

        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)
        sv1, r1 = _make('general', inter)
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')

        sv2, r2 = _make('general', inter)
        g2 = r2.get_param("green")
        g2['chi0q'] = np.stack([np.asarray(g['chi0q'])] * 2, axis=0)
        with self.assertLogs("hwave.solver.rpa", level=logging.INFO) as cm:
            sv2.solve(g2, 'tests/rpa/output')
        self.assertFalse(
            any("partial range" in m for m in cm.output),
            "full-range spin-diag input misreported as partial: %s"
            % [m for m in cm.output if "partial" in m])

    def test_spinful_shaped_input_selects_spinful_mode(self):
        """An in-memory chi0q with nd == norb * ns must dispatch to the
        spinful branch (the file route already did); a zero bubble solves
        to a zero chiq without touching the Green-function machinery."""
        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)
        sv, r = _make('general', inter)
        g = r.get_param("green")
        nvol = sv.lattice.nvol
        nd = sv.nd
        g['chi0q'] = np.zeros((sv.nmat, nvol, nd, nd, nd, nd),
                              dtype=complex)
        sv.solve(g, 'tests/rpa/output')
        self.assertEqual(sv.spin_mode, "spinful")
        self.assertTrue(np.all(np.asarray(g['chiq']) == 0.0))


class TestFileRouteExceptionType(unittest.TestCase):

    def test_malformed_chi0q_init_file_raises_valueerror(self):
        """The public failure type for a malformed chi0q_init file is
        ValueError (standardized by the #109 validation rework; it was a
        bare AssertionError before, which python -O silently disabled)."""
        import tempfile

        sv, _ = _make('general', {'CoulombInter': 'onsite_inter.dat'})
        d = tempfile.mkdtemp()
        np.savez(os.path.join(d, 'chi0q_bad.npz'),
                 chi0q=np.zeros((3, 99, 7, 8), dtype=complex))
        with self.assertRaises(ValueError) as cm:
            sv.read_init({'path_to_input': d, 'chi0q_init': 'chi0q_bad.npz'})
        msg = str(cm.exception)
        self.assertIn('chi0q_bad.npz', msg)
        self.assertIn('(3, 99, 7, 8)', msg)

    def test_partial_range_spin_diag_reports_the_frequency_axis(self):
        """When a spin-diag bubble genuinely carries fewer frequencies than
        Nmat, the partial-range message must report the frequency-axis
        length, not the spin-block count."""
        import logging

        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)
        sv1, r1 = _make('general', inter)
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')

        sv2, r2 = _make('general', inter)
        g2 = r2.get_param("green")
        half = np.asarray(g['chi0q'])[:16]
        g2['chi0q'] = np.stack([half, half], axis=0)
        with self.assertLogs("hwave.solver.rpa", level=logging.INFO) as cm:
            sv2.solve(g2, 'tests/rpa/output')   # completes on 16 frequencies
        self.assertEqual(np.asarray(g2['chiq']).shape[0], 16)
        partial = [m for m in cm.output if "partial range" in m]
        self.assertTrue(partial, "expected a partial-range message")
        self.assertIn("16 in 32", partial[0])


class TestReusedFrequencyMetadata(unittest.TestCase):
    """save_results after a reused-chi0q solve must label the frequency
    axis of the run that PRODUCED the bubble. Before the metadata
    inheritance, a 9-frequency partial bubble reused under a full-range
    config was saved with freq_index 0..31 and nmat 32 -- a mislabeling
    this PR made reachable by fixing the reuse crash."""

    def _make_ranged(self, freq_range):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        r = read_input_k.QLMSkInput(
            {'path_to_input': 'tests/rpa/input_2orb',
             'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                             'Geometry': 'geom.dat',
                             'Transfer': 'transfer.dat',
                             'CoulombInter': 'onsite_inter.dat'}})
        par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
               'SubShape': [1, 1, 1], 'Nmat': 32}
        if freq_range is not None:
            par['matsubara_frequency'] = freq_range
        sv = rpa_mod.RPA(r.get_param("ham"), {},
                         {'mode': 'RPA', 'param': par,
                          'enable_spin_orbital': False,
                          'calc_scheme': 'general', 'calc_type': 'ring'})
        return sv, r

    def test_partial_bubble_keeps_its_producing_axis(self):
        import shutil
        import tempfile

        os.makedirs('tests/rpa/output', exist_ok=True)
        sv1, r1 = self._make_ranged([12, 20])
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')
        self.assertIn('chi0q_freq_meta', g)

        sv2, _ = self._make_ranged(None)
        sv2.solve(g, 'tests/rpa/output')
        d = tempfile.mkdtemp()
        try:
            sv2.save_results({'path_to_output': d, 'chi0q': 'chi0q.npz',
                              'chiq': 'chiq.npz'}, g)
            for f in ('chi0q.npz', 'chiq.npz'):
                z = np.load(os.path.join(d, f))
                np.testing.assert_array_equal(z['freq_index'],
                                              np.arange(12, 21))
                self.assertEqual(int(z['nmat']), 32)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_untagged_partial_bubble_gets_the_ambiguous_label(self):
        """Without producing-run metadata (a hand-built green_info), a
        partial bubble must be labeled 0..n-1 with NO nmat claim -- the
        same explicit ambiguity as a legacy chi0q_init file -- rather
        than the reusing run's full axis."""
        import shutil
        import tempfile

        os.makedirs('tests/rpa/output', exist_ok=True)
        sv1, r1 = self._make_ranged([12, 20])
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')
        g.pop('chi0q_freq_meta')

        sv2, _ = self._make_ranged(None)
        sv2.solve(g, 'tests/rpa/output')
        d = tempfile.mkdtemp()
        try:
            sv2.save_results({'path_to_output': d, 'chi0q': 'chi0q.npz'}, g)
            z = np.load(os.path.join(d, 'chi0q.npz'))
            np.testing.assert_array_equal(z['freq_index'], np.arange(9))
            self.assertNotIn('nmat', z)
        finally:
            shutil.rmtree(d, ignore_errors=True)


class TestRoundThreeHardening(unittest.TestCase):
    """Adversarial-review items: stale provenance, stale results, the
    spin-mode collision, and the ladder/external-bubble hazard."""

    def _solved(self, scheme='general'):
        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)
        sv, r = _make(scheme, inter)
        g = r.get_param("green")
        sv.solve(g, 'tests/rpa/output')
        return sv, g

    def test_stale_metadata_on_a_replaced_bubble_is_rejected(self):
        """Replacing chi0q while leaving the previous solve's metadata must
        fail loudly, not silently mislabel the saved output."""
        sv1, g = self._solved()
        # replace with a differently-sized bubble, keep the stale meta
        g['chi0q'] = np.asarray(g['chi0q'])[:8]
        sv2, _ = _make('general', {'CoulombInter': 'onsite_inter.dat'})
        with self.assertRaises(ValueError) as cm:
            sv2.solve(g, 'tests/rpa/output')
        self.assertIn('stale', str(cm.exception))

    def test_failed_validation_leaves_no_stale_results(self):
        """A malformed reused chi0q must not leave the previous solve's
        chiq/chiq_pm in green_info (they would be saved as this run's)."""
        sv1, g = self._solved()
        g['chi0q'] = np.zeros((3, 5), dtype=complex)
        g.pop('chi0q_freq_meta')
        sv2, _ = _make('general', {'CoulombInter': 'onsite_inter.dat'})
        with self.assertRaises(ValueError):
            sv2.solve(g, 'tests/rpa/output')
        self.assertNotIn('chiq', g)
        self.assertNotIn('chiq_pm', g)

    def test_producer_norb_mismatch_is_rejected(self):
        """Shape-only detection cannot distinguish a spin-free two-orbital
        bubble from a spinful one-orbital one; the recorded producer
        identity must catch the cross-configuration reuse."""
        sv1, g = self._solved()          # produced with norb = 2
        sv2, r2 = _make('general', {'CoulombIntra': 'coulombintra.dat'},
                        path='tests/rpa/input')   # a norb = 1 solver
        g2 = r2.get_param("green")
        g2['chi0q'] = g['chi0q']
        g2['chi0q_freq_meta'] = g['chi0q_freq_meta']
        with self.assertRaises(ValueError) as cm:
            sv2.solve(g2, 'tests/rpa/output')
        self.assertIn('norb', str(cm.exception))

    def test_spin_diag_ladder_rejects_an_external_bubble(self):
        """The spin-diag transverse channel needs the Green functions that
        produced the bubble; a same-instance green0 may belong to an older
        one, so an externally supplied chi0q must be rejected there."""
        import hwave.solver.rpa as rpa_mod
        import hwave.qlmsio.read_input_k as read_input_k

        inter = {'CoulombInter': 'onsite_inter.dat'}
        os.makedirs('tests/rpa/output', exist_ok=True)
        r = read_input_k.QLMSkInput(
            {'path_to_input': 'tests/rpa/input_2orb',
             'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                             'Geometry': 'geom.dat',
                             'Transfer': 'transfer.dat',
                             'CoulombInter': 'onsite_inter.dat'}})
        par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
               'SubShape': [1, 1, 1], 'Nmat': 32}
        sv = rpa_mod.RPA(r.get_param("ham"), {},
                         {'mode': 'RPA', 'param': par,
                          'enable_spin_orbital': False,
                          'calc_scheme': 'general',
                          'calc_type': 'ring+ladder'})
        g = r.get_param("green")
        sv.solve(g, 'tests/rpa/output')      # internal: fine
        # replace the bubble with a spin-diag-shaped external one; the
        # instance still holds green0 from the FIRST bubble
        g['chi0q'] = np.stack([np.asarray(g['chi0q'])[:, :, :2, :2, :2, :2]
                               if np.asarray(g['chi0q']).ndim == 6
                               else np.asarray(g['chi0q'])] * 2, axis=0)
        g.pop('chi0q_freq_meta', None)
        with self.assertRaises(ValueError) as cm:
            sv.solve(g, 'tests/rpa/output')
        self.assertIn('externally supplied', str(cm.exception))

    def test_chained_reuse_keeps_the_producing_axis(self):
        """A bubble passed through TWO consumers must still carry the
        producing run's axis (the write-back covers chains, including
        bubbles that entered through the chi0q_init file route)."""
        import shutil
        import tempfile
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        def make_ranged(fr):
            r = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': {'path_to_input': 'tests/rpa/input_2orb',
                                 'Geometry': 'geom.dat',
                                 'Transfer': 'transfer.dat',
                                 'CoulombInter': 'onsite_inter.dat'}})
            par = {'T': 2.0, 'filling': 0.5, 'CellShape': [4, 4, 1],
                   'SubShape': [1, 1, 1], 'Nmat': 32}
            if fr:
                par['matsubara_frequency'] = fr
            sv = rpa_mod.RPA(r.get_param("ham"), {},
                             {'mode': 'RPA', 'param': par,
                              'enable_spin_orbital': False,
                              'calc_scheme': 'general',
                              'calc_type': 'ring'})
            return sv, r

        os.makedirs('tests/rpa/output', exist_ok=True)
        sv1, r1 = make_ranged([12, 20])
        g = r1.get_param("green")
        sv1.solve(g, 'tests/rpa/output')
        sv2, _ = make_ranged(None)
        sv2.solve(g, 'tests/rpa/output')     # first consumer
        sv3, _ = make_ranged(None)
        sv3.solve(g, 'tests/rpa/output')     # second consumer, via write-back
        d = tempfile.mkdtemp()
        try:
            sv3.save_results({'path_to_output': d, 'chiq': 'chiq.npz'}, g)
            z = np.load(os.path.join(d, 'chiq.npz'))
            np.testing.assert_array_equal(z['freq_index'], np.arange(12, 21))
            self.assertEqual(int(z['nmat']), 32)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
