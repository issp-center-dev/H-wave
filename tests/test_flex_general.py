#!/usr/bin/env python3

"""Tests for the general-mode (paramagnetic full-vertex) FLEX guards.

These tests cover Task 1: ``calc_scheme="general"`` must be accepted for
paramagnetic (spin-free) FLEX and rejected (fail-fast ValueError) for
``enable_spin_orbital`` (at construction) and for non-spin-free spin modes
(spin-diag / spinful, checked inside ``solve()``).

Tests must be run from the repository root (they use relative paths like
``tests/rpa/input``).
"""

import os
import unittest
import numpy as np


def _make_solver(mode_cls, Lx=8, Ly=8, Nmat=64, T=2.0, mu=0.0,
                 calc_scheme="reduced", extra_mode=None,
                 extra_interactions=None):
    """Build a FLEX/RPA solver from the 1-orbital ``tests/rpa/input`` data.

    Replicates the body of ``tests/test_flex.py``'s ``TestFLEX._make_solver``
    helper so these tests are self-contained.
    """
    info_log = {}
    info_mode = {
        'mode': mode_cls,
        'param': {
            'T': T,
            'mu': mu,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
        },
        'calc_scheme': calc_scheme,
    }
    if extra_mode:
        info_mode.update(extra_mode)

    info_file = {
        'input': {
            'path_to_input': 'tests/rpa/input',
            'interaction': {
                'path_to_input': 'tests/rpa/input',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        },
        'output': {
            'path_to_output': 'tests/flex/output',
        },
    }
    if extra_interactions:
        info_file['input']['interaction'].update(extra_interactions)

    os.makedirs('tests/flex/output', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file['input'])
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    if mode_cls == "FLEX":
        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    else:
        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

    return solver, green_info, info_file


class TestFLEXGeneralGuards(unittest.TestCase):
    """Guards for calc_scheme='general' FLEX (v1, paramagnetic only)."""

    def test_general_accepted_for_spin_free(self):
        """calc_scheme='general' must construct without raising (spin-free)."""
        solver, green_info, _ = _make_solver("FLEX", calc_scheme="general")
        self.assertIsNotNone(solver)
        self.assertTrue(solver._flex_general)

    def test_general_rejected_for_enable_spin_orbital(self):
        """calc_scheme='general' + enable_spin_orbital must raise at construct."""
        with self.assertRaises(ValueError):
            _make_solver("FLEX", calc_scheme="general",
                         extra_mode={'enable_spin_orbital': True})

    def test_spin_mode_guard_in_solve(self):
        """The general + non-spin-free guard in solve() must raise.

        Constructing a genuine spin-diag input would require a magnetic-field /
        spin-dependent transfer setup that the 1-orbital fixture does not
        provide.  The guard runs in solve() right after _calc_epsilon_k (which
        is what determines spin_mode from H0), so we stub _calc_epsilon_k to
        report a spin-diag Hamiltonian and confirm solve() fails fast with a
        ValueError before doing any real FLEX work.
        """
        solver, green_info, info_file = _make_solver(
            "FLEX", calc_scheme="general")
        self.assertTrue(solver._flex_general)

        def fake_epsilon_k(gi):
            solver.spin_mode = "spin-diag"

        solver._calc_epsilon_k = fake_epsilon_k
        with self.assertRaises(ValueError):
            solver.solve(green_info, info_file['output']['path_to_output'])

    def test_reduced_still_works(self):
        """Regression: calc_scheme='reduced' FLEX still constructs."""
        solver, green_info, _ = _make_solver("FLEX", calc_scheme="reduced")
        self.assertIsNotNone(solver)
        self.assertFalse(solver._flex_general)

    def test_general_rejected_for_ring_ladder(self):
        """calc_type='ring+ladder' forces scheme='general' in RPA, but FLEX
        general is ring-only — it must still be rejected (not silently accepted
        as a plain general run)."""
        with self.assertRaises(ValueError):
            _make_solver("FLEX", calc_scheme="general",
                         extra_mode={'calc_type': 'ring+ladder'})

    def test_general_rejected_for_spinful_at_solve(self):
        """The solve()-time guard must also reject spin_mode='spinful'."""
        solver, green_info, info_file = _make_solver(
            "FLEX", calc_scheme="general")

        def fake_epsilon_k(gi):
            solver.spin_mode = "spinful"

        solver._calc_epsilon_k = fake_epsilon_k
        with self.assertRaises(ValueError):
            solver.solve(green_info, info_file['output']['path_to_output'])


class TestFLEXGeneralWarningGating(unittest.TestCase):
    """The density-density reduction warning must fire for reduced/squashed but
    be suppressed for general (where the off-diagonal vertices are kept)."""

    def _construct(self, scheme):
        """Construct a 2-orbital FLEX while pretending an Exchange interaction is
        present. ``self.ham_info`` is a fresh ``Interaction`` built inside
        ``RPA.__init__`` (not the passed object), so patch the class method
        ``Interaction.has_interaction_exchange`` rather than an instance — this
        also avoids wiring an extra interaction file."""
        from unittest.mock import patch
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex
        import hwave.solver.rpa as solver_rpa
        info_input = {
            'path_to_input': 'tests/rpa/input_2orb',
            'interaction': {
                'path_to_input': 'tests/rpa/input_2orb',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            },
        }
        ham = read_input_k.QLMSkInput(info_input).get_param("ham")
        info_mode = {
            'mode': 'FLEX',
            'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                      'SubShape': [1, 1, 1], 'Nmat': 32},
            'calc_scheme': scheme,
        }
        with patch.object(solver_rpa.Interaction, 'has_interaction_exchange',
                          return_value=True):
            return solver_flex.FLEX(ham, {}, info_mode)

    def test_warning_retained_for_squashed(self):
        # 'squashed' (not 'reduced': RPA errors on reduced+exchange, rpa.py:643).
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            self._construct('squashed')
        self.assertTrue(any('density-density' in m for m in cm.output))

    def test_warning_suppressed_for_general(self):
        with self.assertNoLogs('hwave.solver.flex', level='WARNING'):
            self._construct('general')


if __name__ == '__main__':
    unittest.main()
