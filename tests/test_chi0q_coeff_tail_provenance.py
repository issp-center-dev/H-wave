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


def _flex_ir_stub(write_densified):
    """General FLEX solver posing as an IR run (no sparse-ir dependency)."""
    from types import SimpleNamespace
    from tests.test_flex_general import _make_general_flex
    flex = _make_general_flex(norb=2)
    flex.coeff_tail = 1.0          # configured but IGNORED on the IR path
    flex.use_ir = True
    flex.write_densified = write_densified
    if not write_densified:
        # sparse-node output path reads the bosonic IR axis metadata
        flex._ir_axB = SimpleNamespace(
            freq_n=np.arange(0, 8, 2), beta=0.5, wmax=10.0, eps=1e-8,
            L=4, statistics="B")
    return flex


@pytest.mark.parametrize('write_densified', [True, False])
def test_flex_ir_save_omits_coeff_tail(write_densified):
    """On the IR path FLEX bypasses the uniform-grid tail machinery
    (flex.py: aa = 0.0 if self.use_ir else self.coeff_tail), so the file must
    NOT claim the configured coeff_tail was applied -- omit the key for both
    IR-native and densified-IR output."""
    flex = _flex_ir_stub(write_densified)
    d = tempfile.mkdtemp()
    nvol = flex.lattice.nvol
    chi0 = np.zeros((flex.nmat, nvol, 2, 2, 2, 2), dtype=complex)
    flex.save_results({'path_to_output': d, 'chi0q': 'chi0q'},
                      {'chi0q': chi0})
    with np.load(os.path.join(d, 'chi0q.npz')) as data:
        assert 'coeff_tail' not in data.files


def test_rpa_read_then_resave_preserves_coeff_tail():
    """End-to-end pass-through: _read_chi0q on a tagged file, then
    save_results -- the INPUT file's coeff_tail must survive the round trip
    even though the current run uses a different value."""
    solver = _make_rpa(coeff_tail=0.0)
    d = tempfile.mkdtemp()
    nvol = solver.lattice.nvol
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex),
             freq_index=np.arange(solver.nmat), nmat=solver.nmat,
             coeff_tail=1.0, tail_endpoint="branch_mean_v1",
             momentum_convention="e_plus_ikR")
    solver._read_chi0q(os.path.join(d, 'chi0q.npz'))
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert float(saved['coeff_tail']) == 1.0
    assert str(saved['tail_endpoint']) == _MARKER


def test_rpa_read_then_resave_omits_for_legacy_file():
    """End-to-end pass-through of a legacy file WITHOUT the key: re-saving
    must not fabricate a coeff_tail claim from the current run's config."""
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    nvol = solver.lattice.nvol
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex),
             freq_index=np.arange(solver.nmat), nmat=solver.nmat, momentum_convention="e_plus_ikR")
    solver._read_chi0q(os.path.join(d, 'chi0q.npz'))
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert 'coeff_tail' not in saved


# ---------------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------------

def _write_chi0q_file(d, nmat=8, nvol=16, coeff_tail=None, marker=True):
    kwargs = {}
    if coeff_tail is not None:
        kwargs['coeff_tail'] = coeff_tail
        if marker:
            # new files carry the endpoint marker; pass marker=False to
            # model a pre-#134 file
            kwargs['tail_endpoint'] = "branch_mean_v1"
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((nmat, nvol, 1, 1), dtype=complex),
             freq_index=np.arange(nmat), nmat=nmat, **kwargs, momentum_convention="e_plus_ikR")


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
             coeff_tail=1.0, tail_endpoint="branch_mean_v1",
             momentum_convention="e_plus_ikR")
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
             coeff_tail=1.0, tail_endpoint="branch_mean_v1",
             momentum_convention="e_plus_ikR")
    with caplog.at_level(logging.WARNING, logger='hwave.solver.rpa'):
        solver._read_chi0q(os.path.join(d, 'chi0q.npz'))
    assert not any('coeff_tail' in m for m in _warnings(caplog))
    assert solver._chi0q_init_meta['coeff_tail'] == 1.0


# ---------------------------------------------------------------------------
# equal-time endpoint convention marker (issue #134)
# ---------------------------------------------------------------------------
#
# The #134 fix changes what a NONZERO-coeff_tail chi0q file contains at
# O(1/Nmat): pre-fix files carry the one-sided-endpoint error and are not
# comparable with a current recomputation at the same Nmat. Files now stamp
# tail_endpoint = "branch_mean_v1"; loaders REJECT a nonzero-tail file
# without the marker (fail-closed, like the momentum-convention gate) and
# any unrecognized marker value. Zero-tail files are exempt: the endpoint
# treatment is vacuous when the tail machinery is off.

_MARKER = "branch_mean_v1"


def _write_rpa_chi0q(solver, d, **extra):
    nvol = solver.lattice.nvol
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex),
             freq_index=np.arange(solver.nmat), nmat=solver.nmat,
             momentum_convention="e_plus_ikR", **extra)
    return os.path.join(d, 'chi0q.npz')


def test_rpa_read_rejects_premarker_nonzero_tail_file():
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=1.0)
    with pytest.raises(ValueError, match="tail_endpoint"):
        solver._read_chi0q(fn)


def test_rpa_read_rejects_unknown_endpoint_marker():
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=1.0,
                          tail_endpoint="branch_mean_v999")
    with pytest.raises(ValueError, match="tail_endpoint"):
        solver._read_chi0q(fn)


def test_rpa_read_accepts_marked_nonzero_tail_file():
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=1.0, tail_endpoint=_MARKER)
    solver._read_chi0q(fn)
    assert solver._chi0q_init_meta['tail_endpoint'] == _MARKER


def test_rpa_read_accepts_zero_tail_file_without_marker():
    solver = _make_rpa(coeff_tail=0.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=0.0)
    solver._read_chi0q(fn)


def test_rpa_read_accepts_legacy_file_without_tail_key():
    """No coeff_tail key at all: unknown tail, unchanged acceptance."""
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d)
    solver._read_chi0q(fn)


def test_rpa_fresh_save_stamps_endpoint_marker():
    solver = _make_rpa(coeff_tail=1.0)
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert str(saved['tail_endpoint']) == _MARKER


def test_rpa_resave_propagates_endpoint_marker():
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=1.0, tail_endpoint=_MARKER)
    solver._read_chi0q(fn)
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert str(saved['tail_endpoint']) == _MARKER


def test_rpa_resave_does_not_fabricate_marker_for_zero_tail_legacy():
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=0.0)
    solver._read_chi0q(fn)
    saved = _save_chi0q(solver, {'chi0q': _stub_chi0q(solver)})
    assert 'tail_endpoint' not in saved


def test_sc_load_rejects_premarker_nonzero_tail_file():
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d, coeff_tail=1.0, marker=False)
    with pytest.raises(ValueError, match="tail_endpoint"):
        sc._load_chi0q(_sc_input(d, coeff_tail=1.0))


def test_sc_load_accepts_marked_nonzero_tail_file():
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d, coeff_tail=1.0)
    sc._load_chi0q(_sc_input(d, coeff_tail=1.0))


def test_sc_load_accepts_zero_tail_file_without_marker():
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    _write_chi0q_file(d, coeff_tail=0.0)
    sc._load_chi0q(_sc_input(d, coeff_tail=0.0))


def test_flex_save_stamps_endpoint_marker():
    """Uniform-grid FLEX chi0q.npz must carry the marker next to
    coeff_tail; the IR path omits both (tail machinery bypassed)."""
    from tests.test_flex_general import _make_general_flex
    flex = _make_general_flex(norb=2)
    flex.coeff_tail = 1.0
    d = tempfile.mkdtemp()
    nvol = flex.lattice.nvol
    chi0 = np.zeros((flex.nmat, nvol, 2, 2, 2, 2), dtype=complex)
    flex.save_results({'path_to_output': d, 'chi0q': 'chi0q'},
                      {'chi0q': chi0})
    with np.load(os.path.join(d, 'chi0q.npz')) as f:
        assert float(f['coeff_tail']) == 1.0
        assert str(f['tail_endpoint']) == _MARKER


@pytest.mark.parametrize('write_densified', [True, False])
def test_flex_ir_save_omits_endpoint_marker(write_densified):
    """IR output must not claim an endpoint treatment that never ran."""
    flex = _flex_ir_stub(write_densified)
    d = tempfile.mkdtemp()
    nvol = flex.lattice.nvol
    chi0 = np.zeros((flex.nmat, nvol, 2, 2, 2, 2), dtype=complex)
    flex.save_results({'path_to_output': d, 'chi0q': 'chi0q'},
                      {'chi0q': chi0})
    with np.load(os.path.join(d, 'chi0q.npz')) as f:
        assert 'tail_endpoint' not in f.files


# ---------------------------------------------------------------------------
# coeff_tail config validation (deep-review should_fix)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("bad", [float('nan'), float('inf'),
                                 -float('inf'), "abc", True, [1.0],
                                 "1.0", np.str_("1.0"), np.bool_(True),
                                 np.array(1.0), np.array([1.0])])
def test_rpa_rejects_pathological_coeff_tail(bad):
    """Type-strict: numeric strings (quoted TOML numbers), booleans and
    arrays -- 0-d included -- are rejected, not coerced."""
    with pytest.raises(ValueError, match="coeff_tail"):
        _make_rpa(coeff_tail=bad)


@pytest.mark.parametrize("ok", [0.0, 1.0, 0.5, -1.0, 2, np.float64(1.0)])
def test_rpa_accepts_finite_real_coeff_tail(ok):
    solver = _make_rpa(coeff_tail=ok)
    assert solver.coeff_tail == float(ok)


# ---------------------------------------------------------------------------
# in-memory reuse route (round-7 review: it bypassed the endpoint gate)
# ---------------------------------------------------------------------------

def _mem_chi0q(solver):
    """Full-frequency general-scheme bubble (the _make_rpa fixture uses
    calc_scheme='general': one orbital, six axes)."""
    nvol = solver.lattice.nvol
    return np.zeros((solver.nmat, nvol, 1, 1, 1, 1), dtype=complex)


def _mem_meta(solver, chi0q, coeff_tail, tail_endpoint):
    from hwave.solver.rpa import _chi0q_fingerprint
    meta = {"freq_index": list(range(solver.nmat)), "nmat": solver.nmat,
            "coeff_tail": coeff_tail, "norb": int(solver.norb),
            "calc_scheme": solver.calc_scheme,
            "fingerprint": _chi0q_fingerprint(chi0q)}
    if tail_endpoint is not None:
        meta["tail_endpoint"] = tail_endpoint
    return meta


def _solve_with_mem_chi0q(solver, chi0q, meta):
    d = tempfile.mkdtemp()
    gi = {"chi0q": chi0q}
    if meta is not None:
        gi["chi0q_freq_meta"] = meta
    solver.solve(gi, d)
    return gi


def test_mem_reuse_rejects_premarker_nonzero_tail():
    solver = _make_rpa(coeff_tail=1.0)
    chi0q = _mem_chi0q(solver)
    meta = _mem_meta(solver, chi0q, 1.0, None)
    with pytest.raises(ValueError, match="tail_endpoint"):
        _solve_with_mem_chi0q(solver, chi0q, meta)


def test_mem_reuse_rejects_unknown_marker():
    solver = _make_rpa(coeff_tail=1.0)
    chi0q = _mem_chi0q(solver)
    meta = _mem_meta(solver, chi0q, 1.0, "branch_mean_v999")
    with pytest.raises(ValueError, match="tail_endpoint"):
        _solve_with_mem_chi0q(solver, chi0q, meta)


def test_mem_reuse_accepts_marked_and_resaves_unchanged():
    solver = _make_rpa(coeff_tail=0.0)      # config differs from producer
    chi0q = _mem_chi0q(solver)
    meta = _mem_meta(solver, chi0q, 1.0, _MARKER)
    gi = _solve_with_mem_chi0q(solver, chi0q, meta)
    saved = _save_chi0q(solver, gi)
    assert float(saved['coeff_tail']) == 1.0
    assert str(saved['tail_endpoint']) == _MARKER


def test_mem_reuse_untagged_full_freq_resaves_without_fabrication():
    """A full-Nmat external chi0q WITHOUT metadata previously fell
    through with _chi0q_init_meta unset, and save_results stamped the
    CURRENT run's coeff_tail and endpoint marker onto data of unknown
    provenance (round-7 review)."""
    solver = _make_rpa(coeff_tail=1.0)
    chi0q = _mem_chi0q(solver)              # full frequency count
    gi = _solve_with_mem_chi0q(solver, chi0q, None)
    saved = _save_chi0q(solver, gi)
    assert 'coeff_tail' not in saved
    assert 'tail_endpoint' not in saved


def test_mem_reuse_end_to_end_chain_preserves_provenance():
    """Legitimate chained workflow (round-8 review): solver 1 COMPUTES a
    coeff_tail=1.0 bubble (no hand-built metadata), solver 2 reuses the
    same green_info, and the save from solver 2 must still carry the
    producing run's coeff_tail and endpoint marker."""
    import hwave.qlmsio.read_input_k as read_input_k
    d = tempfile.mkdtemp()
    solver1 = _make_rpa(coeff_tail=1.0)
    gi = read_input_k.QLMSkInput({
        'path_to_input': 'tests/rpa/input',
        'interaction': {'path_to_input': 'tests/rpa/input',
                        'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat'},
    }).get_param("green")
    solver1.solve(gi, d)
    assert gi["chi0q_freq_meta"]["tail_endpoint"] == _MARKER

    solver2 = _make_rpa(coeff_tail=1.0)
    solver2.solve(gi, d)                    # in-memory reuse route
    saved = _save_chi0q(solver2, gi)
    assert float(saved['coeff_tail']) == 1.0
    assert str(saved['tail_endpoint']) == _MARKER


@pytest.mark.parametrize("bad", [True, "1.0", float('nan'), float('inf')])
def test_mem_reuse_rejects_malformed_tail_provenance(bad):
    """Foreign provenance must not smuggle a malformed coeff_tail past
    the endpoint gate (round-8 review: float() normalized these)."""
    solver = _make_rpa(coeff_tail=1.0)
    chi0q = _mem_chi0q(solver)
    meta = _mem_meta(solver, chi0q, bad, _MARKER)
    with pytest.raises(ValueError, match="coeff_tail"):
        _solve_with_mem_chi0q(solver, chi0q, meta)


@pytest.mark.parametrize("bad", [True, "1.0", float('nan'), float('inf')])
def test_rpa_read_rejects_malformed_tail_provenance(bad):
    solver = _make_rpa(coeff_tail=1.0)
    d = tempfile.mkdtemp()
    fn = _write_rpa_chi0q(solver, d, coeff_tail=bad, tail_endpoint=_MARKER)
    with pytest.raises(ValueError, match="coeff_tail"):
        solver._read_chi0q(fn)


@pytest.mark.parametrize("bad", [True, "1.0", float('nan'), float('inf')])
def test_sc_load_rejects_malformed_tail_provenance(bad):
    import hwave.sc as sc
    d = tempfile.mkdtemp()
    np.savez(os.path.join(d, 'chi0q.npz'),
             chi0q=np.zeros((8, 16, 1, 1), dtype=complex),
             freq_index=np.arange(8), nmat=8, coeff_tail=bad,
             tail_endpoint="branch_mean_v1",
             momentum_convention="e_plus_ikR")
    with pytest.raises(ValueError, match="coeff_tail"):
        sc._load_chi0q(_sc_input(d, coeff_tail=1.0))
