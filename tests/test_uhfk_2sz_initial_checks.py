"""UHFk input checks that used to be silent (#371, #372).

#371: with 2Sz fixed, the spin-block structure removes every
spin-off-diagonal element of the mean-field Hamiltonian, so an interaction
term that changes S_z (PairLift) drops out of the calculation. The run goes
on, as before, and logs one WARNING that names the term.

#372: ``initial_uhf`` uses only rows whose site i is in the origin cell; the
other rows are counted in one INFO line. Two rows that give the same
element must agree to 1e-10 (absolute, on the complex value), otherwise
the input is refused.
"""
import logging
import os
import tempfile

import numpy as np
import pytest
import tomli

from hwave.qlmsio import read_input_k
from hwave.qlmsio.read_input_k import QLMSkInput, DUPLICATE_ROW_TOL
from hwave.solver import hartree_fock
from hwave.solver.uhfk import UHFk


SZ_WARNING = "do not conserve S_z and are ignored under a fixed 2Sz"


def _solver(case, two_sz, iteration_max=None, interaction=None):
    """UHFk on tests/uhfk/<case> through the real reader, with 2Sz set
    (``None`` removes it: Sz-free mode); ``interaction`` replaces the
    ``[file.input.interaction]`` file keys other than Geometry/Transfer."""
    cur = os.getcwd()
    os.chdir(os.path.join("tests", "uhfk", case))
    try:
        with open("input.toml", "rb") as f:
            params = tomli.load(f)
        if two_sz is None:
            params["mode"]["param"].pop("2Sz", None)
        else:
            params["mode"]["param"]["2Sz"] = two_sz
        if iteration_max is not None:
            params["mode"]["param"]["IterationMax"] = iteration_max
        params["log"]["print_level"] = 0
        if interaction is not None:
            files = params["file"]["input"]["interaction"]
            for key in list(files):
                if key not in ("path_to_input", "Geometry", "Transfer"):
                    del files[key]
            files.update(interaction)
        reader = read_input_k.QLMSkInput(params["file"]["input"])
        solver = UHFk(reader.get_param("ham"), params["log"], params["mode"])
        return solver, reader.get_param("green")
    finally:
        os.chdir(cur)


def _detect(case, two_sz, interaction=None):
    solver, _ = _solver(case, two_sz, interaction=interaction)
    solver._make_ham_trans()
    solver._make_ham_inter()
    solver._detect_blocks()
    return solver


def _sz_warnings(caplog):
    return [r for r in caplog.records
            if r.levelno == logging.WARNING and SZ_WARNING in r.getMessage()]


# ---------------------------------------------------------------- #371

def test_sz_breaking_types_selection_rule():
    """Only PairLift's spin table changes S_z (p + q != r + w)."""
    tables = hartree_fock.build_interaction_tables(
        {"Transfer": {}}, 1, (1, 1, 1))
    spin = {}
    inter = {}
    for name in hartree_fock.HF_TYPE_ORDER:
        param = {name: {((0, 0, 0), (0, 0)): 1.0}}
        t = hartree_fock.build_interaction_tables(param, 1, (1, 1, 1))
        inter[name] = t.inter_table[name]
        spin[name] = t.spin_table[name]
    assert hartree_fock.sz_breaking_types(inter, spin) == ["PairLift"]
    assert hartree_fock.sz_breaking_types(tables.inter_table,
                                          tables.spin_table) == []
    # a declared but all-zero PairLift table contributes nothing
    inter["PairLift"] = np.zeros_like(inter["PairLift"])
    assert hartree_fock.sz_breaking_types(inter, spin) == []


def test_pairlift_with_fixed_2sz_warns_once_per_run(caplog):
    solver, green_info = _solver("PairLift", 0, iteration_max=1)
    with tempfile.TemporaryDirectory() as tmp, \
            caplog.at_level(logging.WARNING):
        solver.solve(green_info, tmp)
    warnings = _sz_warnings(caplog)
    assert len(warnings) == 1
    assert "PairLift" in warnings[0].getMessage()
    assert "2Sz" in warnings[0].getMessage()


def test_pairlift_contribution_is_removed_by_the_2sz_blocks():
    """The warning is true: with 2Sz fixed, the Hamiltonian is the same with
    and without PairLift even for a Green function with up-down parts,
    while the unmasked PairLift term is nonzero there."""
    solver = _detect("PairLift", 0)
    rng = np.random.default_rng(7)
    shape = (solver.nvol, 2, solver.norb, 2, solver.norb)
    green = 0.1 * rng.standard_normal(shape) + 0j
    solver.Green = green
    solver._make_ham()
    with_pl = solver.ham.copy()
    raw = np.zeros_like(with_pl)
    hartree_fock.accumulate_hf(
        raw, green, {k: (v if k == "PairLift" else None)
                     for k, v in solver.inter_table.items()},
        solver.spin_table, solver.shape, include_fock=solver.iflag_fock)
    assert np.max(np.abs(raw)) > 1e-3
    solver.inter_table["PairLift"] = None
    solver._make_ham()
    np.testing.assert_array_equal(with_pl, solver.ham)


@pytest.mark.parametrize("case, interaction", [
    ("CoulombIntra", None),
    ("Exchange", None),
    # tests/uhfk/Hund has no data files: read the Exchange case's
    # two-site table as a Hund term instead
    ("Exchange", {"Hund": "exchange.dat"}),
], ids=["CoulombIntra", "Exchange", "Hund"])
def test_sz_conserving_terms_with_fixed_2sz_do_not_warn(case, interaction,
                                                         caplog):
    with caplog.at_level(logging.WARNING):
        solver = _detect(case, 0, interaction)
    if interaction is not None:
        assert solver.inter_table["Hund"] is not None
    assert _sz_warnings(caplog) == []


def test_pairlift_with_free_2sz_does_not_warn(caplog):
    with caplog.at_level(logging.WARNING):
        _detect("PairLift", None)
    assert _sz_warnings(caplog) == []


# ---------------------------------------------------------------- #372

def _geom(positions):
    return {"site2vec": {n: (x, y, z, a)
                         for n, (x, y, z, a) in enumerate(positions)}}


@pytest.fixture
def chain(tmp_path):
    """A 4x1x1 one-orbital chain (SubShape [1, 1, 1]) read from files, with
    the given initial_uhf rows."""
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "geom.dat").write_text("1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
    (raw / "transfer.dat").write_text(
        "Transfer\n1\n2\n1 1\n-1 0 0 1 1 -1.0 0.0\n1 0 0 1 1 -1.0 0.0\n")
    (raw / "geometry_uhf.dat").write_text(
        "1 0 0\n0 1 0\n0 0 1\n0 0 0\n4 0 0\n0 1 0\n0 0 1\n"
        + "".join("{} 0 0 0\n".format(x) for x in range(4)))

    def make(rows):
        (raw / "initial_uhf.dat").write_text(rows)
        reader = QLMSkInput({
            "path_to_input": os.path.relpath(raw),
            "initial_uhf": "initial_uhf.dat",
            "geometry_uhf": "geometry_uhf.dat",
            "interaction": {"path_to_input": os.path.relpath(raw),
                            "Geometry": "geom.dat",
                            "Transfer": "transfer.dat"}})
        solver = UHFk(reader.get_param("ham"), {"print_level": 0}, {
            "mode": "UHFk", "param": {
                "CellShape": [4, 1, 1], "SubShape": [1, 1, 1], "Ncond": 4,
                "2Sz": 0, "T": 0.1, "EPS": 12, "IterationMax": 1,
                "Mix": 0.5, "RndSeed": 1}})
        return solver, reader.get_param("green")
    return make


def _all_rows(value_of):
    lines = []
    for i in range(4):
        for s in range(2):
            for j in range(4):
                for t in range(2):
                    v = value_of(i, s, j, t)
                    lines.append("{} {} {} {} {!r} {!r}\n".format(
                        i, s, j, t, v.real, v.imag))
    return "".join(lines)


def _translation_invariant(i, s, j, t):
    d = (i - j) % 4
    return complex(0.1 * (d + 1) + 0.01 * s, 0.02 * t - 0.03 * d)


def test_rows_outside_origin_cell_are_counted(chain, caplog):
    solver, green_info = chain(_all_rows(_translation_invariant))
    with caplog.at_level(logging.INFO):
        green = solver._initial_green(green_info)
    info = [r.getMessage() for r in caplog.records
            if "outside the origin cell were not used" in r.getMessage()]
    # 3 of the 4 sites are outside the origin cell: 3 * 2 * 4 * 2 rows
    assert len(info) == 1
    assert info[0].startswith("initial_uhf: 48 distinct (i, s, j, t) row(s)")
    # values: G[d, s, 0, t, 0] from the origin row (i = 0, j = -d mod 4)
    for d in range(4):
        for s in range(2):
            for t in range(2):
                assert green[d, s, 0, t, 0] == _translation_invariant(
                    0, s, (-d) % 4, t)


def test_origin_rows_only_log_nothing(chain, caplog):
    rows = "".join(line for line in _all_rows(_translation_invariant)
                   .splitlines(True) if line.startswith("0 "))
    solver, green_info = chain(rows)
    with caplog.at_level(logging.INFO):
        solver._initial_green(green_info)
    assert "outside the origin cell" not in caplog.text


def test_equal_duplicate_rows_are_accepted(chain):
    row = "0 0 1 0 0.25 -0.5\n"
    solver, green_info = chain(row + row)
    assert green_info["initial_uhf"] == {(0, 0, 1, 0): 0.25 - 0.5j}
    green = solver._initial_green(green_info)
    assert green[3, 0, 0, 0, 0] == 0.25 - 0.5j


def test_conflicting_repeated_row_is_refused(chain, caplog):
    rows = "0 0 1 0 0.25 0\n0 0 1 0 0.5 0\n"
    with pytest.raises(SystemExit) as error:
        chain(rows)
    assert error.value.code == 1
    assert "row 0 0 1 0 appears twice with different values" in caplog.text
    assert "(0.25+0j)" in caplog.text and "(0.5+0j)" in caplog.text


def test_duplicate_within_tolerance_is_accepted(chain):
    rows = "0 0 1 0 0.25 0\n0 0 1 0 {!r} 0\n".format(
        0.25 + 0.5 * DUPLICATE_ROW_TOL)
    solver, green_info = chain(rows)
    solver._initial_green(green_info)


def test_distinct_rows_mapping_to_one_element(chain, caplog):
    """Site 4 at x = 5 is site 1 (x = 1) after the modulo L_x = 4, so the
    rows (0 0 1 0) and (0 0 4 0) give the same element."""
    solver, _ = chain("0 0 0 0 0.5 0\n")
    positions = [(x, 0, 0, 0) for x in range(4)] + [(5, 0, 0, 0)]
    geom = _geom(positions)
    same = {(0, 0, 1, 0): 0.25 + 0j, (0, 0, 4, 0): 0.25 + 0j}
    green = solver._initial_green_uhf(same, geom)
    assert green[3, 0, 0, 0, 0] == 0.25
    differ = {(0, 0, 1, 0): 0.25 + 0j, (0, 0, 4, 0): 0.75 + 0j}
    with pytest.raises(SystemExit) as error:
        solver._initial_green_uhf(differ, geom)
    assert error.value.code == 1
    msg = caplog.text
    assert "rows (0 0 1 0) and (0 0 4 0) give the same element" in msg
    assert "G[s=0, a=0; t=0, b=0](d=[3, 0, 0])" in msg
    assert "(0.25+0j)" in msg and "(0.75+0j)" in msg


def test_repeated_rows_cannot_drift_past_the_tolerance(chain, caplog):
    """Each repeat is compared with the first value: 0, 0.75e-10 and
    1.5e-10 are refused although neighbouring values agree."""
    step = 0.75 * DUPLICATE_ROW_TOL
    rows = "".join("0 0 1 0 {!r} 0\n".format(k * step) for k in range(3))
    with pytest.raises(SystemExit):
        chain(rows)
    assert "row 0 0 1 0 appears twice with different values" in caplog.text


def test_element_collisions_cannot_drift_past_the_tolerance(chain, caplog):
    solver, _ = chain("0 0 0 0 0.5 0\n")
    geom = _geom([(x, 0, 0, 0) for x in range(4)] + [(5, 0, 0, 0), (9, 0, 0, 0)])
    step = 0.75 * DUPLICATE_ROW_TOL
    rows = {(0, 0, j, 0): complex(k * step) for k, j in enumerate((1, 4, 5))}
    with pytest.raises(SystemExit):
        solver._initial_green_uhf(rows, geom)
    assert "rows (0 0 1 0) and (0 0 5 0) give the same element" in caplog.text


def _legacy_initial_green_uhf(solver, info, geom):
    """The construction before #372, verbatim apart from names."""
    lx, ly, lz = solver.cellshape
    lvol, norb, ns = solver.cellvol, solver.norb, solver.ns

    def _pack_site(ix, iy, iz):
        return iz + lz * (iy + ly * ix)
    tbl = {idx: (_pack_site(ix, iy, iz), a)
           for idx, (ix, iy, iz, a) in geom["site2vec"].items()}
    green = np.zeros((lvol, lvol, ns, norb, ns, norb), dtype=np.complex128)
    for (i, s, j, t), v in info.items():
        isite, a = tbl[i]
        _, b = tbl[j]
        ix, iy, iz, _ = geom["site2vec"][i]
        jx, jy, jz, _ = geom["site2vec"][j]
        green[isite, _pack_site((ix - jx) % lx, (iy - jy) % ly,
                                (iz - jz) % lz), s, a, t, b] = v
    return green[0]


def test_valid_input_is_bitwise_unchanged(chain):
    """A full-lattice file (every i): the same bytes as before."""
    rng = np.random.default_rng(3)
    solver, green_info = chain(_all_rows(
        lambda i, s, j, t: complex(*rng.standard_normal(2))))
    new = solver._initial_green(green_info)
    old = _legacy_initial_green_uhf(solver, green_info["initial_uhf"],
                                    green_info["geometry_uhf"])
    assert new.dtype == old.dtype and new.shape == old.shape
    assert new.tobytes() == old.tobytes()
