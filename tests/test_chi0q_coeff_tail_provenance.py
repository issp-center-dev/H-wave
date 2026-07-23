#!/usr/bin/env python3

"""coeff_tail provenance for chi0q files (issue #80).

The Matsubara high-frequency tail correction ``coeff_tail`` changes chi0q at
O(1) (lambda differed by 50% in the issue-#80 reproduction), yet the produced
``chi0q.npz`` did not record which value was used: a ``chi0q_mode="calc"`` run
whose config lacks ``coeff_tail`` silently recomputes DIFFERENT physics than a
``"load"`` run reading a tail-corrected file.

Contract pinned here:
- RPA and FLEX ``save_results`` stamp ``coeff_tail`` into chi0q.npz.
- An RPA run that passes a ``chi0q_init`` file through re-saves the INPUT
  file's coeff_tail (or omits the key if the input did not record one) --
  stamping the current run's value would mislabel the data, mirroring the
  freq_index/nmat provenance handling.
- ``hwave_sc``'s chi0q loader and RPA's chi0q_init reader warn when the
  file's recorded coeff_tail differs from the config's effective value.
- Files without the key (older builds) load silently as before.

Tests must be run from the repository root (they use ``tests/rpa/input``).
"""

import logging
import os
import tempfile

import numpy as np
import pytest


def _make_rpa(coeff_tail=None, Lx=4, Ly=4, Nmat=8):
    """Build a 1-orbital RPA solver from the ``tests/rpa/input`` fixture."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    info_input = {
        'path_to_input': 'tests/rpa/input',
        'interaction': {
            'path_to_input': 'tests/rpa/input',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }
    param = {'T': 2.0, 'mu': 0.0,
             'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1], 'Nmat': Nmat}
    if coeff_tail is not None:
        param['coeff_tail'] = coeff_tail
    info_mode = {'mode': 'RPA', 'param': param, 'calc_scheme': 'general'}
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    solver = sol_rpa.RPA(ham, {}, info_mode)
    solver._init_wavevec()
    return solver


def _save_chi0q(solver, green_info):
    d = tempfile.mkdtemp()
    solver.save_results({'path_to_output': d, 'chi0q': 'chi0q'}, green_info)
    with np.load(os.path.join(d, 'chi0q.npz')) as data:
        return {k: data[k] for k in data.files}


def _stub_chi0q(solver):
    nvol = solver.lattice.nvol
    return np.zeros((solver.nmat, nvol, 1, 1), dtype=complex)


# ---------------------------------------------------------------------------
# writers
# ---------------------------------------------------------------------------

def test_rpa_save_records_coeff_tail():
    solver = _make_rpa(coeff_tail=1.0)
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert 'coeff_tail' in saved
    assert float(saved['coeff_tail']) == 1.0


def test_rpa_save_records_default_coeff_tail():
    solver = _make_rpa()          # default coeff_tail = 0.0
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert float(saved['coeff_tail']) == 0.0


def test_rpa_resave_keeps_input_files_coeff_tail():
    """chi0q_init pass-through: the re-saved file must carry the INPUT file's
    coeff_tail, not the current run's."""
    solver = _make_rpa(coeff_tail=0.0)
    solver._chi0q_init_meta = {
        'freq_index': np.arange(solver.nmat),
        'nmat': solver.nmat,
        'coeff_tail': 1.0,
    }
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert float(saved['coeff_tail']) == 1.0


def test_rpa_resave_omits_key_when_input_lacked_it():
    """A legacy chi0q_init file without the key: re-saving must not fabricate
    a coeff_tail claim."""
    solver = _make_rpa(coeff_tail=1.0)
    solver._chi0q_init_meta = {
        'freq_index': np.arange(solver.nmat),
        'nmat': solver.nmat,
        'coeff_tail': None,
    }
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert 'coeff_tail' not in saved


def test_flex_save_records_coeff_tail():
    from tests.test_flex_general import _make_general_flex
    flex = _make_general_flex(norb=2)
    flex.coeff_tail = 1.0
    d = tempfile.mkdtemp()
    nvol = flex.lattice.nvol
    chi0 = np.zeros((flex.nmat, nvol, 2, 2, 2, 2), dtype=complex)
    flex.save_results({'path_to_output': d, 'chi0q': 'chi0q'},
                      {'chi0q': chi0})
    with np.load(os.path.join(d, 'chi0q.npz')) as data:
        assert float(data['coeff_tail']) == 1.0


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------

def _write_chi0q_file(d, nmat=8, nvol=16, coeff_tail=None):
    kwargs = {}
    if coeff_tail is not None:
        kwargs['coeff_tail'] = coeff_tail
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((nmat, nvol, 1, 1), dtype=complex),
             freq_index=np.arange(nmat), nmat=nmat, **kwargs)


def _sc_input(d, nmat=8, coeff_tail=None):
    param = {'T': 2.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
             'Nmat': nmat, 'filling': 0.5}
    if coeff_tail is not None:
        param['coeff_tail'] = coeff_tail
    return {'file': {'output': {'path_to_output': d}},
            'mode': {'param': param}, 'eliashberg': {}}


def _warnings(caplog):
    return [r.getMessage() for r in caplog.records
            if r.levelno >= logging.WARNING]


def test_sc_load_chi0q_warns_on_coeff_tail_mismatch(caplog):
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d, coeff_tail=1.0)
    with caplog.at_level(logging.WARNING, logger='hwave_sc'):
        sc._load_chi0q(_sc_input(d))          # config default 0.0 != file 1.0
    assert any('coeff_tail' in m for m in _warnings(caplog))


def test_sc_load_chi0q_silent_on_matching_coeff_tail(caplog):
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d, coeff_tail=1.0)
    with caplog.at_level(logging.WARNING, logger='hwave_sc'):
        sc._load_chi0q(_sc_input(d, coeff_tail=1.0))
    assert not any('coeff_tail' in m for m in _warnings(caplog))


def test_sc_load_chi0q_silent_on_legacy_file(caplog):
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d)                       # no coeff_tail key
    with caplog.at_level(logging.WARNING, logger='hwave_sc'):
        sc._load_chi0q(_sc_input(d, coeff_tail=1.0))
    assert not any('coeff_tail' in m for m in _warnings(caplog))


def test_rpa_read_chi0q_warns_on_mismatch_and_keeps_provenance(caplog):
    solver = _make_rpa(coeff_tail=0.0)
    d = tempfile.mkdtemp()
    nvol = solver.lattice.nvol
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex),
             freq_index=np.arange(solver.nmat), nmat=solver.nmat,
             coeff_tail=1.0)
    with caplog.at_level(logging.WARNING, logger='hwave.solver.rpa'):
        solver._read_chi0q(os.path.join(d, 'chi0q.npz'))
    assert any('coeff_tail' in m for m in _warnings(caplog))
    assert solver._chi0q_init_meta['coeff_tail'] == 1.0


def test_rpa_read_chi0q_silent_on_match(caplog):
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    nvol = solver.lattice.nvol
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex),
             freq_index=np.arange(solver.nmat), nmat=solver.nmat,
             coeff_tail=1.0)
    with caplog.at_level(logging.WARNING, logger='hwave.solver.rpa'):
        solver._read_chi0q(os.path.join(d, 'chi0q.npz'))
    assert not any('coeff_tail' in m for m in _warnings(caplog))
    assert solver._chi0q_init_meta['coeff_tail'] == 1.0
