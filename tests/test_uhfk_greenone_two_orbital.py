"""UHFk greenone.dat and initial_uhf on a two-orbital model, checked by value.

The model is the non-interacting two-orbital chain of
tests/uhfk_orientation/norb2 (asymmetric complex hopping, T_ab(R) differs
from T_ab(-R) and has no inversion symmetry). Without interactions UHFk is
exact, so every element of greenone.dat must equal the exact one-body
density matrix <c^dag_{R_i a s} c_{R_j b t}> from a brute-force
diagonalisation of the real-space Hamiltonian, built independently of any
UHFk code from the documented transfer convention
H = sum T_ab(R) c^dag_a(R_i) c_b(R_i + R).

The cells are [4,1,1] (sign check), [3,2,1] (a 3x2 cell distinguishes
x-fast from z-fast packing of the displacement) and [3,4,3] with extra,
different complex hoppings along y and z (every axis is extended, so a
swap of any two packing roles is caught). The runs go through the public
entry hwave.qlms.run and the initial_uhf files through QLMSkInput, so the
geometry_uhf, onebodyg_uhf and initial_uhf readers are exercised as well.

The initial_uhf tests feed three files back into UHFk and compare the
resulting initial Green function with green.npz of the converged run:
UHFk's own greenone.dat, and the UHFr outputs green (headerless) and
initial (five-line header) of the same cluster.
"""

import itertools
import os

import numpy as np
import pytest

import hwave.qlms as qlms
from hwave.qlmsio import wan90
from hwave.qlmsio.read_input_k import QLMSkInput
from hwave.solver.uhfk import UHFk


FIXTURE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "uhfk_orientation", "norb2")
NORB = 2
TEMPERATURE = 0.05


# Hoppings along y and z for the 3D case (R, a, b, T_ab(R)); each comes with
# its Hermitian partner T_ba(-R) = conj(T_ab(R)). They differ from each
# other and from the x hoppings, and are not inversion symmetric.
EXTRA_3D = [((0, 1, 0), 1, 1, -0.30 + 0.10j), ((0, 1, 0), 1, 2, -0.12 + 0.07j),
            ((0, 0, 1), 2, 2, -0.18 - 0.06j), ((0, 0, 1), 2, 1, -0.05 + 0.09j)]


def _sites(lx, ly, lz):
    """Site index i -> (x, y, z, orbital); i = a + NORB * (x + lx * (y + ly * z))."""
    return {a + NORB * (x + lx * (y + ly * z)): (x, y, z, a)
            for z in range(lz) for y in range(ly) for x in range(lx)
            for a in range(NORB)}


def _write_transfer(path, extra):
    with open(os.path.join(FIXTURE, "transfer.dat")) as fr:
        lines = fr.read().splitlines()
    rows = lines[4:]
    for rvec, a, b, value in extra:
        for r, i, j, v in ((rvec, a, b, value),
                           (tuple(-x for x in rvec), b, a, value.conjugate())):
            rows.append("{} {} {} {} {} {:.16e} {:.16e}".format(
                *r, i, j, v.real, v.imag))
    nrvec = len({tuple(row.split()[:3]) for row in rows})
    with open(path, "w") as fw:
        fw.write("\n".join([lines[0], lines[1], str(nrvec), " ".join(["1"] * nrvec)]
                           + rows) + "\n")


def _write_inputs(directory, lx, ly, lz):
    with open(os.path.join(FIXTURE, "geom.dat")) as src, \
            open(os.path.join(directory, "geom.dat"), "w") as dst:
        dst.write(src.read())
    _write_transfer(os.path.join(directory, "transfer.dat"),
                    EXTRA_3D if lz > 1 else [])
    sites = _sites(lx, ly, lz)
    with open(os.path.join(directory, "geometry_uhf.dat"), "w") as fw:
        fw.write("1 0 0\n0 1 0\n0 0 1\n0 0 0\n")
        fw.write(f"{lx} 0 0\n0 {ly} 0\n0 0 {lz}\n")
        for i in range(len(sites)):
            fw.write("{} {} {} {}\n".format(*sites[i]))
    pairs = list(itertools.product(range(len(sites)), range(2),
                                   range(len(sites)), range(2)))
    with open(os.path.join(directory, "greenone.def"), "w") as fw:
        fw.write(f"=====\nNCisAjs {len(pairs)}\n=====\n=====\n=====\n")
        fw.writelines("{} {} {} {}\n".format(i, s, j, t) for i, s, j, t in pairs)
    return sites, pairs


def _exact_density(directory, shape, ncond, sites):
    """<c^dag_I c_J> with I = i + N s, by diagonalising the cluster H."""
    lx, ly, lz = shape
    nsite = len(sites)
    index = {v: k for k, v in sites.items()}
    hamiltonian = np.zeros((2 * nsite, 2 * nsite), dtype=complex)
    for (rvec, (a, b)), value in wan90.read_w90(
            os.path.join(directory, "transfer.dat")).items():
        for x, y, z in itertools.product(range(lx), range(ly), range(lz)):
            i = index[(x, y, z, a)]
            j = index[((x + rvec[0]) % lx, (y + rvec[1]) % ly,
                       (z + rvec[2]) % lz, b)]
            for s in range(2):
                hamiltonian[i + nsite * s, j + nsite * s] += value
    assert np.allclose(hamiltonian, hamiltonian.conj().T, atol=1e-14)
    energies, vectors = np.linalg.eigh(hamiltonian)
    low, high = energies.min() - 10, energies.max() + 10
    for _ in range(200):
        mu = 0.5 * (low + high)
        count = (1 / (np.exp((energies - mu) / TEMPERATURE) + 1)).sum()
        low, high = (mu, high) if count < ncond else (low, mu)
    occupation = 1 / (np.exp((energies - mu) / TEMPERATURE) + 1)
    return hamiltonian, np.einsum("il,l,jl->ij", vectors.conj(), occupation, vectors)


def _load_rows(path):
    data = np.loadtxt(path, ndmin=2)
    return {tuple(int(v) for v in row[:4]): row[4] + 1j * row[5] for row in data}


def _uhfk_input(directory, extra=None):
    return {"path_to_input": directory, "geometry_uhf": "geometry_uhf.dat",
            "onebodyg_uhf": "greenone.def", **(extra or {}),
            "interaction": {"path_to_input": directory, "Geometry": "geom.dat",
                            "Transfer": "transfer.dat"}}


def _uhfk_mode(shape, ncond):
    return {"mode": "UHFk", "flag_fock": True,
            "param": {"T": TEMPERATURE, "Ncond": ncond, "IterationMax": 50,
                      "EPS": 12, "Mix": 1.0, "RndSeed": 1,
                      "CellShape": list(shape), "SubShape": [1, 1, 1]}}


def _write_uhfr_inputs(directory, hamiltonian, pairs, nsite):
    rows = [(I % nsite, I // nsite, J % nsite, J // nsite, -hamiltonian[I, J])
            for I in range(2 * nsite) for J in range(2 * nsite)
            if abs(hamiltonian[I, J]) > 0]       # UHFr Transfer carries coeff -1
    with open(os.path.join(directory, "trans.def"), "w") as fw:
        fw.write(f"====\nNTransfer {len(rows)}\n====\n====\n====\n")
        for i, s, j, t, v in rows:
            fw.write(f"{i} {s} {j} {t} {v.real:.17e} {v.imag:.17e}\n")
    with open(os.path.join(directory, "greenone_uhfr.def"), "w") as fw:
        fw.write(f"====\nNCisAjs {len(pairs)}\n====\n====\n====\n")
        fw.writelines("{} {} {} {}\n".format(*k) for k in pairs)


@pytest.fixture(scope="module", params=[(4, 1, 1), (3, 2, 1), (3, 4, 3)],
                ids=["4x1x1", "3x2x1", "3x4x3"])
def case(request, tmp_path_factory):
    shape = request.param
    name = "x".join(str(n) for n in shape)
    ncond = 6 * int(np.prod(shape)) // 4
    root = str(tmp_path_factory.mktemp(f"greenone2orb_{name}"))
    sites, pairs = _write_inputs(root, *shape)
    nsite = len(sites)
    hamiltonian, exact = _exact_density(root, shape, ncond, sites)

    qlms.run(input_dict={
        "log": {"print_level": 0, "print_step": 1000},
        "mode": _uhfk_mode(shape, ncond),
        "file": {"input": _uhfk_input(root),
                 "output": {"path_to_output": os.path.join(root, "out_k"),
                            "energy": "energy.dat", "green": "green",
                            "onebodyg": "greenone.dat"}}})

    _write_uhfr_inputs(root, hamiltonian, pairs, nsite)
    qlms.run(input_dict={
        "log": {"print_level": 0, "print_step": 1000},
        "mode": {"mode": "UHFr", "param": {
            "Nsite": nsite, "Ncond": ncond, "T": TEMPERATURE,
            "IterationMax": 50, "EPS": 12, "Mix": 1.0, "RndSeed": 1}},
        "file": {"input": {"path_to_input": root, "OneBodyG": "greenone_uhfr.def",
                           "interaction": {"Trans": "trans.def"}},
                 "output": {"path_to_output": os.path.join(root, "out_r"),
                            "energy": "energy.dat", "green": "green.dat",
                            "initial": "initial.dat"}}})

    with np.load(os.path.join(root, "out_k", "green.npz")) as archive:
        green = archive["green"]
    return {"shape": shape, "name": name, "ncond": ncond, "root": root, "sites": sites,
            "pairs": pairs, "nsite": nsite, "exact": exact, "green": green,
            "onebodyg": _load_rows(os.path.join(root, "out_k", "greenone.dat")),
            "uhfr": _load_rows(os.path.join(root, "out_r", "green.dat"))}


def _exact(case, i, s, j, t):
    nsite = case["nsite"]
    return case["exact"][i + nsite * s, j + nsite * t]


def _max_error(case, values, reference):
    return max(abs(values[k] - reference(*k)) for k in case["pairs"])


def test_model_is_asymmetric(case):
    # Guard: the fixture must make the sign and packing observable.
    exact = case["exact"]
    assert np.max(np.abs(exact.imag)) > 0.05
    assert np.max(np.abs(exact - exact.T)) > 0.05


def test_onebodyg_matches_exact_density(case):
    error = _max_error(case, case["onebodyg"], lambda *k: _exact(case, *k))
    print(f"{case['name']}: onebodyg vs exact = {error:.3e}")
    assert error < 1e-10


def test_onebodyg_matches_uhfr_cluster(case):
    # UHFr on the same cluster defines the meaning of the i s j t Re Im row.
    assert _max_error(case, case["uhfr"], lambda *k: _exact(case, *k)) < 1e-10
    assert _max_error(case, case["onebodyg"], lambda *k: case["uhfr"][k]) < 1e-10


def test_previous_layout_is_rejected(case):
    # The previous writer read G at R_j - R_i, packed x-fast. On this model it
    # is far from the exact value, so the comparison above is discriminating.
    (lx, ly, lz), sites, green = case["shape"], case["sites"], case["green"]

    def previous(i, s, j, t):
        (xi, yi, zi, a), (xj, yj, zj, b) = sites[i], sites[j]
        r = (xj - xi) % lx + lx * ((yj - yi) % ly + ly * ((zj - zi) % lz))
        return green[r, s, a, t, b]

    error = max(abs(previous(*k) - _exact(case, *k)) for k in case["pairs"])
    print(f"{case['name']}: previous layout vs exact = {error:.3e}")
    assert error > 0.1


def _initial_from(case, directory, file_name):
    reader = QLMSkInput(_uhfk_input(directory, {"initial_uhf": file_name}))
    solver = UHFk(reader.get_param("ham"), {"print_level": 0},
                  _uhfk_mode(case["shape"], case["ncond"]))
    return solver._initial_green(reader.get_param("green"))


@pytest.mark.parametrize("source", ["uhfk_greenone", "uhfr_green", "uhfr_initial"])
def test_initial_uhf_round_trip(case, source):
    """initial_uhf reproduces green.npz from each real-space file."""
    root = case["root"]
    file_name = {"uhfk_greenone": os.path.join("out_k", "greenone.dat"),
                 "uhfr_green": os.path.join("out_r", "green.dat"),
                 "uhfr_initial": os.path.join("out_r", "initial.dat")}[source]
    initial = _initial_from(case, root, file_name)
    error = np.max(np.abs(initial - case["green"]))
    print(f"{case['name']} {source}: initial vs green.npz = {error:.3e}")
    # greenone.dat carries 13 significant digits; the UHFr files full precision.
    assert error < (1e-11 if source == "uhfk_greenone" else 1e-10)


def test_initial_uhf_is_inverse_of_writer(case):
    """Full-precision exact values placed by initial_uhf give green.npz."""
    reader = QLMSkInput(_uhfk_input(case["root"]))
    solver = UHFk(reader.get_param("ham"), {"print_level": 0},
                  _uhfk_mode(case["shape"], case["ncond"]))
    values = {k: _exact(case, *k) for k in case["pairs"]}
    geometry = {"site2vec": case["sites"]}
    initial = solver._initial_green_uhf(values, geometry)
    assert np.max(np.abs(initial - case["green"])) < 1e-10
    # And the onebodyg writer applied to that array returns the same values.
    solver.Green = initial
    path = os.path.join(case["root"], "rewritten_greenone.dat")
    solver._save_greenone(path, {"onebodyg_uhf": case["pairs"],
                                 "geometry_uhf": geometry})
    rewritten = _load_rows(path)
    assert max(abs(rewritten[k] - values[k]) for k in case["pairs"]) < 1e-11
