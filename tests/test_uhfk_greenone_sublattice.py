"""Physical one-body output must retain both sites' sublattice offsets.

The independent oracle uses only eigen.npz and physical coordinates. The
legacy-output comparison checks symmetric 4x4 unfolded-path invariance.
An anisotropic 4x3 cell checks displacement packing and initial_uhf restarts;
its Peierls variant distinguishes the displacement sign.

At T=0.01 the unfolded fixtures have partially occupied degenerate levels,
and the stripe ([2,1,1]) fixtures have thermal corrections above 1e-10.
Their oracle therefore uses Fermi occupations from the saved chemical
potentials. The gapped [2,2,1] AF fixtures also check integer band filling.
Diagonal complex hopping alone leaves the [2,2,1] AF density real; an
additional Peierls hopping fixture explicitly tests complex AF matrix entries.
"""

import itertools
import os
from pathlib import Path
import tempfile

import numpy as np
import pytest

from hwave.qlmsio.read_input_k import QLMSkInput
from hwave.solver.uhfk import UHFk


CELL = np.array([4, 4, 1])
POSITIONS = np.array([(x, y, 0) for y in range(4) for x in range(4)])


def _write_w90(path, title, entries):
    lines = [title, '1', str(len(entries)), ' '.join(['1'] * len(entries))]
    for (x, y, z), value in entries:
        lines.append(f'{x} {y} {z} 1 1 {value.real:.16e} {value.imag:.16e}')
    path.write_text('\n'.join(lines) + '\n')


def _run_case(directory, subshape, complex_hopping, peierls_phase=0.0, *,
              cell=(4, 4, 1), ty=-1.0, interaction=4.0, ncond=16,
              temperature=0.01, initial_uhf=None, iteration_max=3000):
    """Run the real input reader, self-consistency loop, and output writers."""
    directory = Path(os.path.relpath(directory))
    cell = np.array(cell)
    positions = np.array([(x, y, z) for z in range(cell[2])
                          for y in range(cell[1]) for x in range(cell[0])])
    volume = len(positions)
    raw = directory / 'raw'
    output = directory / 'output'
    raw.mkdir(parents=True)
    output.mkdir()
    (raw / 'geom.dat').write_text('1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n')
    entries = [((x, y, 0), (-1.0 if x else ty) * np.exp(1j * peierls_phase * x))
               for x, y in [(-1, 0), (1, 0), (0, -1), (0, 1)]]
    if complex_hopping:
        entries += [((1, 1, 0), 0.3 + 0.2j), ((-1, -1, 0), 0.3 - 0.2j)]
    _write_w90(raw / 'transfer.dat', 'Transfer', entries)
    _write_w90(raw / 'coulombintra.dat', 'CoulombIntra', [((0, 0, 0), complex(interaction))])
    geometry = ('1 0 0\n0 1 0\n0 0 1\n0 0 0\n'
                f'{cell[0]} 0 0\n0 {cell[1]} 0\n0 0 {cell[2]}\n')
    geometry += ''.join(f'{x} {y} {z} 0\n' for x, y, z in positions)
    (raw / 'geometry_uhf.dat').write_text(geometry)
    pairs = list(itertools.product(range(volume), range(2), range(volume), range(2)))
    requests = f'=====\nNCisAjs {len(pairs)}\n=====\n=====\n=====\n'
    requests += ''.join(f'{i} {s} {j} {t}\n' for i, s, j, t in pairs)
    (raw / 'OneBodyG.dat').write_text(requests)
    sub = np.array(subshape)
    nf = int(np.prod(sub))
    nv = int(np.prod(cell // sub))
    initial = np.zeros((nv, 2, nf, 2, nf), dtype=complex)
    for y in range(sub[1]):
        for x in range(sub[0]):
            fo = x + sub[0] * y
            moment = 0.4 * (-1) ** (x + y) if nf > 1 else 0.0
            initial[0, 0, fo, 0, fo] = 0.5 + moment
            initial[0, 1, fo, 1, fo] = 0.5 - moment
    np.savez(raw / 'initial.npz', **{
        'green_sublattice' if nf > 1 else 'green': initial})
    initial_input = ({'initial': 'initial.npz'} if initial_uhf is None else
                     {'initial_uhf': os.path.relpath(initial_uhf, raw)})
    reader = QLMSkInput({
        'path_to_input': str(raw), **initial_input,
        'geometry_uhf': 'geometry_uhf.dat', 'onebodyg_uhf': 'OneBodyG.dat',
        'interaction': {'path_to_input': str(raw), 'Geometry': 'geom.dat',
                        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat'},
    })
    solver = UHFk(reader.get_param('ham'), {'print_level': 0}, {
        'mode': 'UHFk', 'param': {
            'CellShape': cell.tolist(), 'SubShape': list(subshape),
            'Ncond': ncond, '2Sz': 0, 'T': temperature, 'EPS': 12,
            'IterationMax': iteration_max, 'Mix': 0.5, 'RndSeed': 1,
        },
    })
    green_info = reader.get_param('green')
    initial_density = solver._initial_green(green_info).copy()
    solver.solve(green_info, str(output))
    if iteration_max > 1:
        assert solver.physics['Rest'] < 1e-12
    solver.save_results({
        'path_to_output': str(output), 'green': 'green.npz',
        'eigen': 'eigen.npz', 'energy': 'energy.dat', 'onebodyg': 'greenone.dat',
    }, green_info)
    data = np.loadtxt(output / 'greenone.dat')
    density = np.empty((volume, 2, volume, 2), dtype=complex)
    indices = data[:, :4].astype(int)
    density[tuple(indices.T)] = data[:, 4] + 1j * data[:, 5]
    energy = dict(line.split('=', 1) for line in (output / 'energy.dat').read_text().splitlines())
    energy = {key.strip(): float(value) for key, value in energy.items()}
    with np.load(output / 'eigen.npz') as archive:
        eigen = {key: archive[key] for key in archive.files}
    with np.load(output / 'green.npz') as archive:
        green = archive['green_sublattice' if nf > 1 else 'green']
    return {'density': density, 'eigen': eigen, 'green': green,
            'energy': energy, 'sub': sub, 'output': output, 'cell': cell,
            'positions': positions, 'temperature': temperature, 'ncond': ncond,
            'initial_density': initial_density, 'peierls_phase': peierls_phase}


@pytest.fixture(scope='module')
def cases():
    with tempfile.TemporaryDirectory(prefix='greenone_', dir='.') as directory:
        yield {(sub, hopping): _run_case(
            Path(directory) / f'sub{sub[0]}{sub[1]}{sub[2]}_complex{int(hopping)}',
            sub, hopping)
            for sub in [(2, 2, 1), (1, 1, 1), (2, 1, 1)]
            for hopping in [False, True]}


def _reconstruct(case, thermal=True, phase_sign=1):
    """Build physical Bloch waves without using any UHFk method."""
    sub = case['sub']
    shape = case['cell'] // sub
    cells = case['positions'] // sub
    within = case['positions'] % sub
    orbitals = within[:, 0] + sub[0] * (within[:, 1] + sub[1] * within[:, 2])
    nf = int(np.prod(sub))
    eigen = case['eigen']
    values, vectors = eigen['eigenvalue'], eigen['eigenvector']
    momenta = 2 * np.pi * eigen['wavevector_index'] / shape
    phases = np.exp(phase_sign * 1j * (cells @ momenta.T)) / np.sqrt(np.prod(shape))
    density = np.zeros((len(cells), 2, len(cells), 2), dtype=complex)
    for spin in range(2):
        energies = values[:, spin * nf:(spin + 1) * nf]
        occupations = np.zeros(energies.size)
        occupations[np.argsort(energies, axis=None, kind='stable')[:case['ncond'] // 2]] = 1
        occupations = occupations.reshape(energies.shape)
        if thermal:
            mu = case['energy'][f'ChemicalPotential_{spin}']
            occupations = 1 / (1 + np.exp(np.clip((energies - mu) / case['temperature'], -700, 700)))
        for k, band in itertools.product(range(len(momenta)), range(nf)):
            psi = vectors[k, spin * nf + orbitals, spin * nf + band] * phases[:, k]
            density[:, spin, :, spin] += occupations[k, band] * np.outer(psi.conj(), psi)
    return density


@pytest.mark.parametrize('sub', [(2, 2, 1), (1, 1, 1), (2, 1, 1)])
@pytest.mark.parametrize('hopping', [False, True], ids=['real', 'complex'])
def test_independent_reconstruction(cases, sub, hopping):
    case = cases[sub, hopping]
    expected = _reconstruct(case)
    error = np.max(np.abs(case['density'] - expected))
    cutoff_error = np.max(np.abs(expected - _reconstruct(case, thermal=False)))
    print(f'{sub=} {hopping=} reconstruction={error:.12e} cutoff_error={cutoff_error:.12e} '
          f'max_imag={np.max(np.abs(expected.imag)):.12e}')
    np.testing.assert_allclose(case['density'], expected, atol=1e-10, rtol=0)


@pytest.mark.parametrize('hopping', [False, True], ids=['real', 'complex'])
def test_gapped_af_integer_band_reconstruction(cases, hopping):
    case = cases[(2, 2, 1), hopping]
    np.testing.assert_allclose(case['density'], _reconstruct(case, thermal=False),
                               atol=1e-10, rtol=0)


def test_folded_displacement_sign(cases):
    # With only two cells per direction, r == -r modulo the shape. The
    # [2,1,1] case has four cells along y and really discriminates the sign.
    case = cases[(2, 1, 1), True]
    expected = _reconstruct(case)
    wrong_sign = _reconstruct(case, phase_sign=-1)
    assert np.max(np.abs(expected.imag)) > 0.008
    wrong_error = np.max(np.abs(case['density'] - wrong_sign))
    print(f'Wrong displacement sign error={wrong_error:.12e}')
    assert wrong_error > 0.15
    np.testing.assert_allclose(case['density'], expected, atol=1e-10, rtol=0)


def test_complex_af_reconstruction_and_hermiticity(tmp_path):
    case = _run_case(tmp_path, (2, 2, 1), True, peierls_phase=0.4)
    expected = _reconstruct(case, thermal=False)
    assert np.max(np.abs(expected.imag)) > 0.01
    error = np.max(np.abs(case['density'] - expected))
    print(f'Peierls AF reconstruction={error:.12e}, '
          f'max_imag={np.max(np.abs(expected.imag)):.12e}')
    np.testing.assert_allclose(case['density'], expected, atol=1e-10, rtol=0)
    density = case['density'].reshape(32, 32)
    np.testing.assert_allclose(density, density.conj().T, atol=1e-12, rtol=0)


def test_af_sum_rules(cases):
    case = cases[(2, 2, 1), False]
    diagonal = case['density'].reshape(32, 32).diagonal().real.reshape(16, 2)
    particles = diagonal.sum()
    spin = (diagonal[:, 0] - diagonal[:, 1]).sum() / 2
    print(f'AF N={particles:.15g}, NCond={case["energy"]["NCond"]:.15g}, '
          f'Sz={spin:.15g}, energy Sz={case["energy"]["Sz"]:.15g}')
    assert abs(particles - case['energy']['NCond']) < 1e-8
    assert abs(spin - case['energy']['Sz']) < 1e-8


def test_af_sublattice_resolution(cases):
    case = cases[(2, 2, 1), False]
    diagonal = case['density'].reshape(32, 32).diagonal().real.reshape(16, 2)
    spin = (diagonal[:, 0] - diagonal[:, 1]) / 2
    within = POSITIONS % case['sub']
    orbitals = within[:, 0] + 2 * within[:, 1]
    green = case['green']
    reference = (green[0, 0, orbitals, 0, orbitals].real
                 - green[0, 1, orbitals, 1, orbitals].real) / 2
    signs = (-1.) ** (POSITIONS[:, 0] + POSITIONS[:, 1])
    print(f'AF site Sz={spin[:4]}, staggered={np.mean(signs * spin):.15g}, '
          f'folded staggered={np.mean(signs * reference):.15g}')
    assert np.min(signs * spin) > 0.3
    np.testing.assert_allclose(spin, reference, atol=1e-10, rtol=0)
    assert abs(np.mean(signs * (spin - reference))) < 1e-10


@pytest.mark.parametrize('sub', [(2, 2, 1), (1, 1, 1), (2, 1, 1)])
def test_complex_hermiticity(cases, sub):
    density = cases[sub, True]['density'].reshape(32, 32)
    np.testing.assert_allclose(density, density.conj().T, atol=1e-12, rtol=0)


@pytest.mark.parametrize('hopping', [False, True], ids=['real', 'complex'])
def test_unfolded_path_invariance(cases, hopping):
    """Unfolded-path invariance: reproduce the pre-change writer to 1e-14."""
    case = cases[(1, 1, 1), hopping]
    before = np.empty((16, 2, 16, 2), dtype=complex)
    for i, s, j, t in itertools.product(range(16), range(2), range(16), range(2)):
        dx, dy, dz = (POSITIONS[j] - POSITIONS[i] + CELL) % CELL
        value = case['green'][dx + 4 * (dy + 4 * dz), s, 0, t, 0]
        before[i, s, j, t] = complex(float(f'{value.real:.12e}'), float(f'{value.imag:.12e}'))
    np.testing.assert_allclose(case['density'], before, atol=1e-14, rtol=0)


@pytest.fixture(scope='module', params=[0.0, 0.4], ids=['real', 'peierls'])
def anisotropic_case(request):
    with tempfile.TemporaryDirectory(prefix='greenone_rect_', dir='.') as directory:
        yield _run_case(Path(directory), (1, 1, 1), False,
                        peierls_phase=request.param, cell=(4, 3, 1), ty=-0.4,
                        interaction=2.0, ncond=6, temperature=0.05)


def test_anisotropic_reconstruction(anisotropic_case):
    case = anisotropic_case
    expected = _reconstruct(case)
    error = np.max(np.abs(case['density'] - expected))
    print(f'Anisotropic max_imag={np.max(np.abs(expected.imag)):.12e}, '
          f'reconstruction={error:.12e}')
    np.testing.assert_allclose(case['density'], expected, atol=1e-10, rtol=0)
    if case['peierls_phase'] == 0.0:
        bond_x = case['density'][0, 0, 1, 0]
        bond_y = case['density'][0, 0, 4, 0]
        print(f'Anisotropic x={bond_x:.12e}, y={bond_y:.12e}, '
              f'ratio={abs(bond_x / bond_y):.12e}')
        np.testing.assert_allclose(bond_x, 0.249944, atol=5e-7, rtol=0)
        np.testing.assert_allclose(bond_y, 0.000084, atol=5e-7, rtol=0)
        assert abs(bond_x) > 100 * abs(bond_y)
    else:
        # Real hopping is inversion symmetric and cannot discriminate d from -d.
        assert np.max(np.abs(expected.imag)) > 0.07
        wrong_error = np.max(np.abs(case['density'] - _reconstruct(case, phase_sign=-1)))
        print(f'Anisotropic opposite displacement error={wrong_error:.12e}')
        assert wrong_error > 0.1


@pytest.mark.parametrize('with_header', [False, True], ids=['headerless', 'legacy_header'])
def test_initial_uhf_round_trip(anisotropic_case, tmp_path, with_header):
    case = anisotropic_case
    initial_uhf = case['output'] / 'greenone.dat'
    if with_header:
        # Preserve the old five-line input format as well as direct output reuse.
        intermediate = tmp_path / 'intermediate'
        intermediate.mkdir()
        initial_uhf = intermediate / 'initial_with_header.dat'
        initial_uhf.write_text('=====\nNInitial 576\n=====\n=====\n=====\n'
                               + (case['output'] / 'greenone.dat').read_text())
    restarted = _run_case(tmp_path, (1, 1, 1), False, peierls_phase=case['peierls_phase'],
                          cell=(4, 3, 1), ty=-0.4, interaction=2.0, ncond=6,
                          temperature=0.05, iteration_max=1,
                          initial_uhf=initial_uhf)
    input_error = np.max(np.abs(restarted['initial_density'] - case['green']))
    error = np.max(np.abs(restarted['green'] - case['green']))
    print(f'Round trip input={input_error:.12e}, after_one_iteration={error:.12e}')
    np.testing.assert_allclose(restarted['initial_density'], case['green'], atol=1e-10, rtol=0)
    np.testing.assert_allclose(restarted['green'], case['green'], atol=1e-8, rtol=0)


def test_initial_uhf_rejects_sublattice(tmp_path, caplog):
    initial_uhf = tmp_path / 'initial_uhf.dat'
    initial_uhf.write_text('0 0 0 0 0.5 0\n0 1 0 1 0.5 0\n')
    with pytest.raises(SystemExit) as error:
        _run_case(tmp_path, (2, 1, 1), False, initial_uhf=initial_uhf)
    assert error.value.code == 1
    assert ('initial_uhf requires SubShape = [1, 1, 1]; give a sublattice-ordered '
            'initial state through initial (green.npz) instead') in caplog.text
