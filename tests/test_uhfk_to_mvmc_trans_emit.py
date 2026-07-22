"""Tests for tools/_uhfk_to_mvmc/trans_emit.py."""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

from tools._uhfk_to_mvmc import trans_emit
from tools._uhfk_to_mvmc.trans_emit import (
    TransEmitError,
    _unpack_soi,
    emit_trans_def,
    parse_hwave_transfer,
)

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _write_transfer(path, rows):
    """Write a minimal Wannier90-like Transfer.dat with the given rows.

    ``rows`` is a list of ``(rx, ry, rz, iWan, jWan, re, im)`` tuples.
    """
    n = len(rows)
    with open(path, "w") as fw:
        fw.write("test Transfer.dat header\n")
        fw.write("1\n")  # num_wann placeholder
        fw.write(f"{n}\n")
        ones = ["1"] * n
        # ndegen: 15 per line.
        for i in range(0, n, 15):
            fw.write(" ".join(ones[i:i + 15]) + "\n")
        for (rx, ry, rz, iWan, jWan, re, im) in rows:
            fw.write(
                f"  {rx:4d} {ry:4d} {rz:4d} {iWan:4d} {jWan:4d}  "
                f"{re:.12f}  {im:.12f}\n"
            )


@pytest.mark.parametrize("a, s", [(0, 0), (0, 1)])
def test_unpack_soi_roundtrip(a, s):
    """norb_orig=1 permits only (a, s) = (0, 0), (0, 1).

    See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    """
    iWan = 2 * a + s + 1
    assert _unpack_soi(iWan) == (a, s)


def test_unpack_soi_rejects_nonpositive():
    with pytest.raises(TransEmitError, match="must be >= 1"):
        _unpack_soi(0)


def test_parse_hwave_transfer_reads_all_rows(tmp_path):
    """Verify parse strips header/nr/ndegen and returns all data rows."""
    src = tmp_path / "Transfer.dat"
    rows = [
        (1, 0, 0, 1, 1, -1.0, 0.0),
        (-1, 0, 0, 1, 1, -1.0, 0.0),
        (0, 1, 0, 1, 2, 0.0, 0.5),
    ]
    _write_transfer(str(src), rows)
    entries = parse_hwave_transfer(str(src))
    assert len(entries) == 3
    assert entries[0] == (1, 0, 0, 1, 1, complex(-1.0, 0.0))
    assert entries[2] == (0, 1, 0, 1, 2, complex(0.0, 0.5))


def test_parse_hwave_transfer_rejects_short_row(tmp_path):
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("header\n")
        fw.write("1\n")
        fw.write("1\n")
        fw.write("1\n")  # ndegen line
        fw.write("  0 0 0 1 1  0.5\n")  # only 6 tokens
    with pytest.raises(TransEmitError, match="expected 7 columns"):
        parse_hwave_transfer(str(src))


def test_parse_hwave_transfer_divides_by_ndegen(tmp_path):
    """Wannier90-format Transfer.dat with non-unit ndegen must divide each
    coefficient by the R-point degeneracy, matching H-wave's
    ``hwave.qlmsio.wan90.read_w90``. Fixture stores ``raw = target *
    ndegen[R]`` so the parser must return exactly ``target`` on every row.
    """
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("test Wannier90 header with non-unit ndegen\n")
        fw.write("2\n")               # num_wann (unused)
        fw.write("3\n")               # nr = 3 distinct R-vectors
        fw.write("1 2 1\n")           # ndegen for R0, R1, R2
        # Dense listing: each R appears with raw = 2.0 * ndegen[R] so
        # every returned value must be exactly 2.0.
        fw.write("   0   0   0   1   1   2.0   0.0\n")   # R0, deg=1
        fw.write("   1   0   0   1   1   4.0   0.0\n")   # R1, deg=2
        fw.write("  -1   0   0   1   1   2.0   0.0\n")   # R2, deg=1
    entries = parse_hwave_transfer(str(src))
    assert len(entries) == 3
    for entry in entries:
        _, _, _, _, _, val = entry
        assert val == pytest.approx(complex(2.0, 0.0)), (
            f"got {val}, expected 2.0 (post-ndegen division)"
        )


def test_parse_hwave_transfer_all_unit_ndegen_no_op(tmp_path):
    """The existing case_soc_rashba_2d_nosub fixture uses all-ones ndegen,
    so division must be a no-op and pre-existing test values pass through
    unchanged. Guards the 0.20% E2E delta against silent regression.
    """
    src = tmp_path / "Transfer.dat"
    rows = [
        (1, 0, 0, 1, 1, -1.0, 0.0),
        (0, 1, 0, 1, 2, 0.0, 0.5),
    ]
    _write_transfer(str(src), rows)
    entries = parse_hwave_transfer(str(src))
    assert entries[0] == (1, 0, 0, 1, 1, complex(-1.0, 0.0))
    assert entries[1] == (0, 1, 0, 1, 2, complex(0.0, 0.5))


def test_parse_hwave_transfer_rejects_non_unit_sparse(tmp_path):
    """Non-unit ndegen paired with a sparse listing (fewer distinct Rs
    than declared) must fail, matching H-wave's read_w90."""
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("header\n")
        fw.write("1\n")
        fw.write("2\n")               # nr = 2 distinct R-vectors declared
        fw.write("1 2\n")             # non-unit ndegen
        # Only one distinct R present -> sparse.
        fw.write("   0   0   0   1   1   1.0   0.0\n")
    with pytest.raises(TransEmitError, match="non-unit ndegen"):
        parse_hwave_transfer(str(src))


def test_parse_hwave_transfer_rejects_ndegen_count_mismatch(tmp_path):
    """Declared nr does not match the number of ndegen ints found."""
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("header\n")
        fw.write("1\n")
        fw.write("3\n")               # nr = 3
        fw.write("1 1\n")             # only 2 ints (rows = 1, but count = 2)
        fw.write("   0   0   0   1   1   1.0   0.0\n")
    with pytest.raises(TransEmitError, match="ndegen count mismatch"):
        parse_hwave_transfer(str(src))


def test_parse_hwave_transfer_rejects_duplicate_entries(tmp_path):
    """Duplicate (R, iWan, jWan) entries are rejected, matching H-wave's
    ``hwave.qlmsio.wan90.read_w90`` duplicate-entry contract. Without this
    check ``emit_trans_def`` would double-count the hopping in trans.def
    and silently corrupt the SOC Hamiltonian.
    """
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("test header\n")
        fw.write("2\n")               # num_wann (unused)
        fw.write("1\n")               # nr = 1
        fw.write("1\n")               # ndegen
        fw.write("   0   0   0   1   1   1.0   0.0\n")
        fw.write("   0   0   0   1   1   1.0   0.0\n")  # duplicate key
    with pytest.raises(TransEmitError, match="duplicate"):
        parse_hwave_transfer(str(src))


def test_parse_hwave_transfer_skips_blank_and_comment_rows(tmp_path):
    src = tmp_path / "Transfer.dat"
    with open(str(src), "w") as fw:
        fw.write("header\n")
        fw.write("1\n")
        fw.write("2\n")
        fw.write("1 1\n")
        fw.write("  0 0 0 1 1  1.0  0.0\n")
        fw.write("\n")  # blank
        fw.write("# comment inside data block\n")
        fw.write("  1 0 0 1 1  -1.0  0.0\n")
    entries = parse_hwave_transfer(str(src))
    assert len(entries) == 2


def test_emit_trans_def_expands_all_sites(tmp_path):
    """2x2x1 lattice, 1 Transfer entry at R=(1,0,0), up->dn ->
    4 rows (one per source site), each with i != j and swapped spins."""
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    # iWan = 1 (a=0, s=0=up), jWan = 2 (a=0, s=1=dn), val = +i (Rashba-like).
    _write_transfer(str(src), [(1, 0, 0, 1, 2, 0.0, 1.0)])
    emit_trans_def(str(src), [2, 2, 1], str(out))
    data_rows = _read_trans_def_rows(str(out))
    assert len(data_rows) == 4  # one per source site
    for i, s, j, t, re, im in data_rows:
        assert s == 1 and t == 0  # H-wave spin endpoints are swapped
        assert i != j
        # PBC: for source site (ix, iy, 0), target = ((ix+1) mod 2, iy, 0).
        Lx = 2
        ix = i % Lx
        iy = (i // Lx) % 2
        jx = j % Lx
        jy = (j // Lx) % 2
        assert (jx, jy) == ((ix + 1) % Lx, iy)


def test_emit_trans_def_spin_diagonal_preserves_minus_real_value(tmp_path):
    """The derived rule reduces to the byte-verified ``-v`` for real
    spin-diagonal entries because swapping equal spins is a no-op.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    # 1x1x1 lattice, iWan=jWan=1 (spin-diagonal), val=+1.0.
    _write_transfer(str(src), [(0, 0, 0, 1, 1, 1.0, 0.0)])
    emit_trans_def(str(src), [1, 1, 1], str(out))
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 1
    _, _, _, _, re, im = rows[0]
    assert re == pytest.approx(-1.0)
    assert im == pytest.approx(0.0)


def test_emit_trans_def_spin_diagonal_conjugates_complex_value(tmp_path):
    """Complex diagonal coefficients follow the same general ``-conj(v)`` rule."""
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(str(src), [(0, 0, 0, 1, 1, 0.3, 0.4)])
    emit_trans_def(str(src), [1, 1, 1], str(out))
    rows = _read_trans_def_rows(str(out))
    _, s, _, t, re, im = rows[0]
    assert (s, t) == (0, 0)
    assert complex(re, im) == pytest.approx(-0.3 + 0.4j)


def test_emit_trans_def_spin_symmetric_complex_soc_swaps_and_conjugates(
    tmp_path,
):
    """A spin-symmetric SOC block exposes the old ``+v`` bug.

    Both off-diagonal H-wave spin entries carry ``0.3+0.4j``. The general
    mapping swaps the spin labels and emits ``-conj(v)``; the row order pins
    which source Transfer entry produced each emitted row. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    value = 0.3 + 0.4j
    _write_transfer(str(src), [
        (0, 0, 1, 1, 2, value.real, value.imag),
        (0, 0, 1, 2, 1, value.real, value.imag),
    ])
    emit_trans_def(str(src), [1, 1, 1], str(out))
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 2
    assert [(row[1], row[3]) for row in rows] == [(1, 0), (0, 1)]
    emitted_values = [complex(row[4], row[5]) for row in rows]
    assert emitted_values == pytest.approx([-value.conjugate()] * 2)
    assert all(emitted != pytest.approx(value) for emitted in emitted_values), (
        "regression must differ from the old spin-preserving +v rule"
    )


def test_emit_trans_def_v36_rashba_matches_old_plus_v_matrix(tmp_path):
    """The x/y Rashba spin matrices hide the old rule's defect.

    At fixed R they obey ``v[t,s] = -conj(v[s,t])``: the x pair has equal
    imaginary entries and the y pair has antisymmetric real entries.
    Therefore swapping spin endpoints and emitting ``-conj(v)`` produces
    exactly the same matrix as the old spin-preserving ``+v`` rule.

    This pins MATRIX equivalence only. Comparing dicts deliberately
    discards row order and signed-zero formatting, because the emitted
    text is NOT byte-identical to the old rule's: off-diagonal rows carry
    swapped spin labels and the partner row supplies the transposed
    entry. Byte stability across this change is not a contract; see the
    module docstring in ``trans_emit.py`` and
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    transfer_rows = [
        (1, 0, 0, 1, 2, 0.0, +0.5),  # +x: up -> down
        (1, 0, 0, 2, 1, 0.0, +0.5),  # +x: down -> up
        (0, 1, 0, 1, 2, -0.5, 0.0),  # +y: up -> down
        (0, 1, 0, 2, 1, +0.5, 0.0),  # +y: down -> up
    ]
    _write_transfer(str(src), transfer_rows)
    cell_shape = [2, 2, 1]
    emit_trans_def(str(src), cell_shape, str(out))
    rows = _read_trans_def_rows(str(out))

    new_matrix = {
        (i, s, j, t): complex(re, im) for i, s, j, t, re, im in rows
    }
    old_plus_v_matrix = {}
    Lx, Ly, _ = cell_shape
    for rx, ry, _rz, i_wan, j_wan, re, im in transfer_rows:
        for iy in range(Ly):
            for ix in range(Lx):
                i_site = ix + Lx * iy
                j_site = (ix + rx) % Lx + Lx * ((iy + ry) % Ly)
                old_plus_v_matrix[
                    (i_site, i_wan - 1, j_site, j_wan - 1)
                ] = complex(re, im)

    assert new_matrix == pytest.approx(old_plus_v_matrix)


def test_emit_trans_def_matches_vmcdry_format_prefix(tmp_path):
    """The header must mirror vmcdry.out's five-line preamble so mVMC
    accepts the file interchangeably."""
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(str(src), [(0, 0, 0, 1, 1, 1.0, 0.0)])
    emit_trans_def(str(src), [1, 1, 1], str(out))
    with open(str(out)) as fp:
        lines = fp.readlines()
    assert lines[0].startswith("=")
    assert "NTransfer" in lines[1]
    assert "1" in lines[1]  # count
    assert lines[2].startswith("=")
    assert "i_j_s_tijs" in lines[3]
    assert lines[3].startswith("=")
    assert lines[4].startswith("=")


def test_emit_trans_def_rejects_non_norb1_orbital(tmp_path):
    """With norb_orig=1, iWan=3 (a_phys=1) must fail loudly.

    See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(str(src), [(0, 0, 0, 3, 1, 1.0, 0.0)])
    with pytest.raises(TransEmitError, match="norb_orig == 1"):
        emit_trans_def(str(src), [1, 1, 1], str(out))


def test_emit_trans_def_rejects_bad_cell_shape(tmp_path):
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(str(src), [(0, 0, 0, 1, 1, 1.0, 0.0)])
    with pytest.raises(TransEmitError, match="length 3"):
        emit_trans_def(str(src), [4, 4], str(out))


def test_validate_trans_def_entries_rejects_unsupported_orbital():
    entries = [(0, 0, 0, 3, 1, 1.0 + 0.0j)]

    with pytest.raises(TransEmitError, match="norb_orig == 1"):
        trans_emit.validate_trans_def_entries(entries, [1, 1, 1])


def test_validate_trans_def_entries_accepts_supported_entries():
    entries = [
        (1, 0, 0, 1, 2, -1.0 + 0.5j),
        (-1, 0, 0, 2, 1, -1.0 - 0.5j),
    ]

    trans_emit.validate_trans_def_entries(entries, [4, 1, 1])


def test_emit_trans_def_from_entries_materializes_generator(tmp_path):
    entry = (0, 0, 0, 1, 1, 1.0 + 0.0j)
    entries = (item for item in [entry])
    out = tmp_path / "trans.def"

    trans_emit.emit_trans_def_from_entries(
        entries, [1, 1, 1], str(out),
    )

    lines = out.read_text().splitlines()
    assert lines[1].split() == ["NTransfer", "1"]
    assert len(lines[5:]) == 1


def test_emit_trans_def_pbc_wraps_across_edge(tmp_path):
    """For a Lx=3 lattice with R=(1,0,0), ix=2 -> jx = (2+1) mod 3 = 0."""
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(str(src), [(1, 0, 0, 1, 1, -1.0, 0.0)])
    emit_trans_def(str(src), [3, 1, 1], str(out))
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 3
    # Expect (i=0, j=1), (i=1, j=2), (i=2, j=0).
    pairs = sorted((r[0], r[2]) for r in rows)
    assert pairs == [(0, 1), (1, 2), (2, 0)]


def _make_soc_fixture(tmp, *, transfer_body):
    """Assemble a minimal SOC bridge fixture at ``tmp`` and return the
    argv suffix (paths only). Mirrors the skeleton in
    ``tests/test_uhfk_to_mvmc_dispatch.py::_write_fixture`` but is
    inlined here so the trans_emit test module has no cross-test
    imports. Writes an ``input.toml`` with
    ``enable_spin_orbital = true`` and PBC.
    """
    Nsite = 2
    Ncond = 2
    with open(os.path.join(tmp, "input.toml"), "w") as fp:
        fp.write(
            "[mode.param]\n"
            f"Ncond = {Ncond}\n"
            "2Sz = 0\n"
            "T = 0.0\n"
            f"CellShape = [{Nsite}, 1, 1]\n"
            "SubShape  = [1, 1, 1]\n"
            'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
            "enable_spin_orbital = true\n"
        )
    eigenvalue = np.array([[0.0, 1.0], [0.5, 1.5]], dtype=np.float64)
    np.savez(
        os.path.join(tmp, "eigen.npz"),
        eigenvalue=eigenvalue,
        eigenvector=np.eye(2, dtype=np.complex128).reshape(1, 2, 2)
        .repeat(Nsite, axis=0),
        wavevector_unit=np.eye(3, dtype=np.float64),
        wavevector_index=np.array(
            [[v, 0, 0] for v in [0, 1]], dtype=np.int64,
        ),
        twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
    )
    occ = np.zeros((Nsite, 2), dtype=np.float64)
    occ[0, 0] = 1.0
    occ[0, 1] = 1.0
    np.savez(
        os.path.join(tmp, "occupation.npz"),
        occupation=occ,
        mu=np.array([0.0, 0.0], dtype=np.float64),
        T=np.float64(0.0),
        column_spin=np.array([0, 1], dtype=np.int64),
        column_mu_group=np.array([0, 1], dtype=np.int64),
    )
    with open(os.path.join(tmp, "geometry_uhf.dat"), "w") as fp:
        fp.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
        fp.write("0.0 0.0 0.0\n")
        fp.write(f"{Nsite} 0 0\n0 1 0\n0 0 1\n")
        for i in range(Nsite):
            fp.write(f"{i} 0 0 0\n")
    # Minimal 6-column General orbitalidx.def.
    total = 2 * Nsite * Nsite - Nsite
    lines = [
        "======================",
        f"NOrbitalIdx  {total}",
        "ComplexType 1",
        "======================",
        "== i_spn_j_spn_OrbitalIdx ==",
        "======================",
    ]
    idx = 0
    for all_i in range(2 * Nsite):
        for all_j in range(all_i + 1, 2 * Nsite):
            i, spn_i = all_i % Nsite, all_i // Nsite
            j, spn_j = all_j % Nsite, all_j // Nsite
            lines.append(f"{i} {spn_i} {j} {spn_j} {idx} 1")
            idx += 1
    for k in range(total):
        lines.append(f"{k} 1")
    with open(os.path.join(tmp, "orbitalidx.def"), "w") as fp:
        fp.write("\n".join(lines) + "\n")
    transfer_path = os.path.join(tmp, "Transfer.dat")
    with open(transfer_path, "w") as fp:
        fp.write(transfer_body)
    return transfer_path


def _run_soc_bridge(tmp_path):
    """Run the product CLI against a fixture rooted at ``tmp_path``."""
    return subprocess.run(
        [
            sys.executable, "tools/uhfk_to_mvmc.py",
            "--input", str(tmp_path / "input.toml"),
            "--eigen", str(tmp_path / "eigen.npz"),
            "--occupation", str(tmp_path / "occupation.npz"),
            "--geometry", str(tmp_path / "geometry_uhf.dat"),
            "--orbitalidx", str(tmp_path / "orbitalidx.def"),
            "--output", str(tmp_path / "zqp_orbital_uhfk.dat"),
            "--no-check-density",
            "--epsilon-noise", "0",
            "--transfer", str(tmp_path / "Transfer.dat"),
            "--emit-trans", str(tmp_path / "trans.def.bridge"),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


def _make_soc_sublattice_fixture(tmp_path, *, transfer_body=None):
    """Copy the tracked SOC + SubShape fixture into ``tmp_path``."""
    fixture_dir = os.path.join(
        REPO_ROOT,
        "tests/validation/uhfk_mvmc_pairproduct/"
        "case_soc_rashba_2d_sub_apbc",
    )
    data_dir = os.path.join(REPO_ROOT, "tests/data")
    for name in ("input.toml", "geometry_uhf.dat"):
        shutil.copyfile(
            os.path.join(fixture_dir, name),
            tmp_path / name,
        )
    shutil.copyfile(
        os.path.join(data_dir, "v36_case_soc_rashba_2d_sub_apbc_eigen.npz"),
        tmp_path / "eigen.npz",
    )
    shutil.copyfile(
        os.path.join(
            data_dir,
            "v36_case_soc_rashba_2d_sub_apbc_occupation.npz",
        ),
        tmp_path / "occupation.npz",
    )
    if transfer_body is None:
        shutil.copyfile(
            os.path.join(fixture_dir, "Transfer.dat"),
            tmp_path / "Transfer.dat",
        )
    else:
        (tmp_path / "Transfer.dat").write_text(transfer_body)

    # The source fixture's orbitalidx files are generated E2E artifacts and
    # intentionally not tracked. The CLI only needs a valid General-format
    # input here; SOC + SubShape replaces it with --emit-orbitalidx.
    nsite = 24
    total = nsite * (2 * nsite - 1)
    lines = [
        "======================",
        f"NOrbitalIdx  {total}",
        "ComplexType 1",
        "======================",
        "== i_spn_j_spn_OrbitalIdx ==",
        "======================",
    ]
    idx = 0
    for all_i in range(2 * nsite):
        for all_j in range(all_i + 1, 2 * nsite):
            i, spn_i = all_i % nsite, all_i // nsite
            j, spn_j = all_j % nsite, all_j // nsite
            lines.append(f"{i} {spn_i} {j} {spn_j} {idx} 1")
            idx += 1
    lines.extend(f"{idx} 1" for idx in range(total))
    (tmp_path / "orbitalidx.def").write_text("\n".join(lines) + "\n")


def _run_soc_sublattice_bridge(tmp_path, **path_overrides):
    """Run all SOC + SubShape output producers against ``tmp_path``."""
    paths = {
        "input": tmp_path / "input.toml",
        "eigen": tmp_path / "eigen.npz",
        "occupation": tmp_path / "occupation.npz",
        "geometry": tmp_path / "geometry_uhf.dat",
        "orbitalidx": tmp_path / "orbitalidx.def",
        "output": tmp_path / "zqp_orbital_uhfk.dat",
        "transfer": tmp_path / "Transfer.dat",
        "emit_trans": tmp_path / "trans.def.bridge",
        "emit_orbitalidx": tmp_path / "orbitalidxgen.def.bridge",
    }
    paths.update(path_overrides)
    return subprocess.run(
        [
            sys.executable, "tools/uhfk_to_mvmc.py",
            "--input", str(paths["input"]),
            "--eigen", str(paths["eigen"]),
            "--occupation", str(paths["occupation"]),
            "--geometry", str(paths["geometry"]),
            "--orbitalidx", str(paths["orbitalidx"]),
            "--output", str(paths["output"]),
            "--no-check-density",
            "--epsilon-noise", "0",
            "--debug-writer",
            "--transfer", str(paths["transfer"]),
            "--emit-trans", str(paths["emit_trans"]),
            "--emit-orbitalidx", str(paths["emit_orbitalidx"]),
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )


@pytest.mark.parametrize(
    "input_name",
    [
        "input.toml",
        "eigen.npz",
        "occupation.npz",
        "geometry_uhf.dat",
        "orbitalidx.def",
        "Transfer.dat",
    ],
)
def test_cli_rejects_emit_orbitalidx_matching_any_input(tmp_path, input_name):
    """A requested output must never overwrite a direct CLI input."""
    _make_soc_sublattice_fixture(tmp_path)
    input_path = tmp_path / input_name
    original = input_path.read_bytes()

    result = _run_soc_sublattice_bridge(
        tmp_path, emit_orbitalidx=input_path,
    )

    assert result.returncode == 2, result.stderr
    assert "refusing to overwrite an input" in result.stderr
    assert input_path.read_bytes() == original
    for name in (
        "F_pre_noise.npz",
        "F_post_aggregate.npz",
        "zqp_orbital_uhfk.dat",
        "trans.def.bridge",
    ):
        assert not (tmp_path / name).exists(), f"unexpected output: {name}"


@pytest.mark.parametrize("alias_kind", ["symlink", "hardlink"])
def test_cli_rejects_emit_orbitalidx_aliasing_transfer_input(
    tmp_path, alias_kind,
):
    """Existing symlink and hardlink aliases must not bypass preflight."""
    _make_soc_sublattice_fixture(tmp_path)
    transfer_path = tmp_path / "Transfer.dat"
    alias_path = tmp_path / "transfer-alias.dat"
    if alias_kind == "symlink":
        alias_path.symlink_to(transfer_path)
    else:
        os.link(transfer_path, alias_path)
    original = transfer_path.read_bytes()

    result = _run_soc_sublattice_bridge(
        tmp_path, emit_orbitalidx=alias_path,
    )

    assert result.returncode == 2, result.stderr
    assert "refusing to overwrite an input" in result.stderr
    assert transfer_path.read_bytes() == original


def test_cli_rejects_colliding_outputs_before_writing(tmp_path):
    """All active output destinations must be mutually distinct."""
    _make_soc_sublattice_fixture(tmp_path)
    shared_output = tmp_path / "shared-output.dat"

    result = _run_soc_sublattice_bridge(
        tmp_path,
        output=shared_output,
        emit_orbitalidx=shared_output,
    )

    assert result.returncode == 2, result.stderr
    assert (
        "output destinations resolve to the same path or file"
        in result.stderr
    )
    for name in (
        "shared-output.dat",
        "F_pre_noise.npz",
        "F_post_aggregate.npz",
        "trans.def.bridge",
    ):
        assert not (tmp_path / name).exists(), f"unexpected output: {name}"


def test_cli_rejects_debug_output_matching_input(tmp_path):
    """The fixed --debug-writer paths participate in collision checks."""
    _make_soc_sublattice_fixture(tmp_path)
    eigen_path = tmp_path / "F_pre_noise.npz"
    (tmp_path / "eigen.npz").replace(eigen_path)
    original = eigen_path.read_bytes()

    result = _run_soc_sublattice_bridge(tmp_path, eigen=eigen_path)

    assert result.returncode == 2, result.stderr
    assert "refusing to overwrite an input" in result.stderr
    assert eigen_path.read_bytes() == original
    for name in (
        "F_post_aggregate.npz",
        "orbitalidxgen.def.bridge",
        "zqp_orbital_uhfk.dat",
        "trans.def.bridge",
    ):
        assert not (tmp_path / name).exists(), f"unexpected output: {name}"


def test_cli_validates_transfer_before_soc_sublattice_outputs(tmp_path):
    """A rejected Transfer.dat must not overwrite any earlier producer."""
    _make_soc_sublattice_fixture(
        tmp_path,
        transfer_body=(
            "test header\n"
            "1\n"
            "1\n"
            "1\n"
            "  1 0 0 1 1  -1.0 0.0\n"
        ),
    )
    sentinels = {
        "orbitalidxgen.def.bridge": b"PREEXISTING_ORBITALIDX\n",
        "F_pre_noise.npz": b"PREEXISTING_F_PRE\n",
        "F_post_aggregate.npz": b"PREEXISTING_F_POST\n",
    }
    for name, content in sentinels.items():
        (tmp_path / name).write_bytes(content)

    result = _run_soc_sublattice_bridge(tmp_path)

    assert result.returncode == 2, result.stderr
    assert "ERROR (transfer_hermiticity)" in result.stderr
    assert "missing partner" in result.stderr
    for name, content in sentinels.items():
        assert (tmp_path / name).read_bytes() == content, (
            f"{name} was overwritten before Transfer.dat validation"
        )
    assert not (tmp_path / "zqp_orbital_uhfk.dat").exists()
    assert not (tmp_path / "trans.def.bridge").exists()


def test_cli_soc_sublattice_emits_all_requested_outputs(tmp_path):
    """Hermitian input still produces orbital, debug, zqp, and trans files."""
    _make_soc_sublattice_fixture(tmp_path)

    result = _run_soc_sublattice_bridge(tmp_path)

    assert result.returncode == 0, result.stderr
    for name in (
        "orbitalidxgen.def.bridge",
        "F_pre_noise.npz",
        "F_post_aggregate.npz",
        "zqp_orbital_uhfk.dat",
        "trans.def.bridge",
    ):
        assert (tmp_path / name).is_file(), f"missing output: {name}"


def test_cli_rejects_nonhermitian_transfer_without_final_outputs(tmp_path):
    """The product CLI must fail closed before emitting the bridge pair."""
    _make_soc_fixture(
        str(tmp_path),
        transfer_body=(
            "test header\n"
            "1\n"
            "1\n"
            "1\n"
            "  1 0 0 1 1  -1.0 0.0\n"
        ),
    )

    result = _run_soc_bridge(tmp_path)

    assert result.returncode == 2, result.stderr
    assert "ERROR (transfer_hermiticity)" in result.stderr
    assert "missing partner" in result.stderr
    assert not (tmp_path / "zqp_orbital_uhfk.dat").exists()
    assert not (tmp_path / "trans.def.bridge").exists()


def test_cli_rejects_empty_transfer_without_final_outputs(tmp_path):
    """A bridge output with no transfer terms is not a valid product."""
    _make_soc_fixture(
        str(tmp_path),
        transfer_body="test header\n1\n0\n",
    )

    result = _run_soc_bridge(tmp_path)

    assert result.returncode == 2, result.stderr
    assert "ERROR (transfer_hermiticity)" in result.stderr
    assert "bridge output requires at least one transfer term" in result.stderr
    assert not (tmp_path / "zqp_orbital_uhfk.dat").exists()
    assert not (tmp_path / "trans.def.bridge").exists()


def test_cli_accepts_hermitian_transfer(tmp_path):
    """A valid Hermitian pair continues through both bridge outputs."""
    _make_soc_fixture(
        str(tmp_path),
        transfer_body=(
            "test header\n"
            "1\n"
            "2\n"
            "1 1\n"
            "  1 0 0 1 1  -1.0 0.0\n"
            " -1 0 0 1 1  -1.0 0.0\n"
        ),
    )

    result = _run_soc_bridge(tmp_path)

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "zqp_orbital_uhfk.dat").is_file()
    assert (tmp_path / "trans.def.bridge").is_file()


def test_atomic_output_no_partial_state_on_failure(tmp_path):
    """An unsupported SOC orbital must fail before any requested output.

    Pre-existing zqp, trans, and debug artifacts all survive the failed call.
    """
    tmp = str(tmp_path)
    # This real same-index on-site term passes Hermiticity, but iWan=jWan=3
    # decodes to physical orbital 1, outside the single-orbital scope. See
    # ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
    transfer_body = (
        "test header\n"
        "3\n"
        "1\n"
        "1\n"
        "  0 0 0 3 3  0.5 0.0\n"
    )
    _make_soc_fixture(tmp, transfer_body=transfer_body)

    zqp_path = os.path.join(tmp, "zqp_orbital_uhfk.dat")
    trans_path = os.path.join(tmp, "trans.def.bridge")
    # Pre-existing sentinel files at the two output destinations.
    sentinel_zqp = b"PREEXISTING_ZQP_CONTENT\n"
    sentinel_trans = b"PREEXISTING_TRANS_CONTENT\n"
    sentinel_f_pre = b"PREEXISTING_F_PRE\n"
    sentinel_f_post = b"PREEXISTING_F_POST\n"
    with open(zqp_path, "wb") as fp:
        fp.write(sentinel_zqp)
    with open(trans_path, "wb") as fp:
        fp.write(sentinel_trans)
    (tmp_path / "F_pre_noise.npz").write_bytes(sentinel_f_pre)
    (tmp_path / "F_post_aggregate.npz").write_bytes(sentinel_f_post)

    result = subprocess.run(
        [
            sys.executable, "tools/uhfk_to_mvmc.py",
            "--input", os.path.join(tmp, "input.toml"),
            "--eigen", os.path.join(tmp, "eigen.npz"),
            "--occupation", os.path.join(tmp, "occupation.npz"),
            "--geometry", os.path.join(tmp, "geometry_uhf.dat"),
            "--orbitalidx", os.path.join(tmp, "orbitalidx.def"),
            "--output", zqp_path,
            "--no-check-density",
            "--epsilon-noise", "0",
            "--debug-writer",
            "--transfer", os.path.join(tmp, "Transfer.dat"),
            "--emit-trans", trans_path,
        ],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert result.returncode != 0, result.stderr
    assert "trans_emit" in result.stderr, result.stderr
    # Pre-existing files at the destinations must be untouched.
    with open(zqp_path, "rb") as fp:
        assert fp.read() == sentinel_zqp, (
            "zqp_orbital_uhfk.dat corrupted despite SOC emit failure"
        )
    with open(trans_path, "rb") as fp:
        assert fp.read() == sentinel_trans, (
            "trans.def corrupted despite SOC emit failure"
        )
    assert (tmp_path / "F_pre_noise.npz").read_bytes() == sentinel_f_pre
    assert (tmp_path / "F_post_aggregate.npz").read_bytes() == sentinel_f_post
    # This rejection happens during preflight, so the absence of temps here
    # confirms that preflight did not enter the atomic-write transaction.
    stray = [
        n for n in os.listdir(tmp)
        if n.startswith(".uhfk_to_mvmc.") and n.endswith(".tmp")
    ]
    assert stray == [], f"preflight unexpectedly created temps: {stray}"


def test_atomic_output_cleans_created_temps_when_replace_fails(
    tmp_path, monkeypatch,
):
    """A post-mkstemp commit fault cleans temps and preserves outputs."""
    from tools import uhfk_to_mvmc

    transfer_body = (
        "test header\n"
        "1\n"
        "2\n"
        "1 1\n"
        "  1 0 0 1 1  -1.0 0.0\n"
        " -1 0 0 1 1  -1.0 0.0\n"
    )
    _make_soc_fixture(str(tmp_path), transfer_body=transfer_body)
    zqp_path = tmp_path / "zqp_orbital_uhfk.dat"
    trans_path = tmp_path / "trans.def.bridge"
    sentinel_zqp = b"PREEXISTING_ZQP_CONTENT\n"
    sentinel_trans = b"PREEXISTING_TRANS_CONTENT\n"
    zqp_path.write_bytes(sentinel_zqp)
    trans_path.write_bytes(sentinel_trans)

    created_temps = []
    real_mkstemp = tempfile.mkstemp

    def recording_mkstemp(*args, **kwargs):
        fd, path = real_mkstemp(*args, **kwargs)
        created_temps.append(path)
        return fd, path

    def fail_replace(_source, _destination):
        raise OSError("injected os.replace failure")

    monkeypatch.setattr(
        uhfk_to_mvmc.tempfile, "mkstemp", recording_mkstemp,
    )
    monkeypatch.setattr(uhfk_to_mvmc.os, "replace", fail_replace)

    with pytest.raises(OSError, match="injected os.replace failure"):
        uhfk_to_mvmc.main(
            [
                "--input", str(tmp_path / "input.toml"),
                "--eigen", str(tmp_path / "eigen.npz"),
                "--occupation", str(tmp_path / "occupation.npz"),
                "--geometry", str(tmp_path / "geometry_uhf.dat"),
                "--orbitalidx", str(tmp_path / "orbitalidx.def"),
                "--output", str(zqp_path),
                "--no-check-density",
                "--epsilon-noise", "0",
                "--transfer", str(tmp_path / "Transfer.dat"),
                "--emit-trans", str(trans_path),
            ]
        )

    assert len(created_temps) == 2
    assert all(not os.path.exists(path) for path in created_temps)
    assert zqp_path.read_bytes() == sentinel_zqp
    assert trans_path.read_bytes() == sentinel_trans


def _read_trans_def_rows(path):
    """Parse mVMC trans.def emitted by ``emit_trans_def`` and return the
    data rows as list of ``(i, s, j, t, re, im)`` tuples."""
    with open(path) as fp:
        lines = fp.readlines()
    # First 5 lines are the header (mirrors vmcdry). Data starts at line 5.
    rows = []
    for line in lines[5:]:
        toks = line.split()
        if not toks:
            continue
        i, s, j, t = int(toks[0]), int(toks[1]), int(toks[2]), int(toks[3])
        re, im = float(toks[4]), float(toks[5])
        rows.append((i, s, j, t, re, im))
    return rows


# ---------------------------------------------------------------------------
# SOC + APBC boundary_theta wrap-phase behavior. See
# ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
# ---------------------------------------------------------------------------


def test_emit_trans_def_no_gauge_when_boundary_theta_none(tmp_path):
    """Default (boundary_theta=None) is byte-identical to explicit
    all-PBC (boundary_theta=(0, 0, 0)).
    """
    src = tmp_path / "Transfer.dat"
    _write_transfer(src, [
        (1, 0, 0, 1, 1, -1.0, 0.0),
        (-1, 0, 0, 1, 1, -1.0, 0.0),
        (1, 0, 0, 1, 2, 0.0, 0.5),
        (-1, 0, 0, 1, 2, 0.0, -0.5),
    ])
    out_default = tmp_path / "trans_default.def"
    out_pbc = tmp_path / "trans_pbc.def"
    emit_trans_def(str(src), [4, 1, 1], str(out_default))
    emit_trans_def(
        str(src), [4, 1, 1], str(out_pbc), boundary_theta=(0.0, 0.0, 0.0),
    )
    with open(out_default) as fp:
        default_bytes = fp.read()
    with open(out_pbc) as fp:
        pbc_bytes = fp.read()
    assert default_bytes == pbc_bytes


def test_emit_trans_def_apbc_flips_sign_for_boundary_crossing_spin_diag(
    tmp_path,
):
    """1D L=4 lattice, R=(1,0,0), spin-diagonal T=-1.0 under APBC(x).
    Interior bonds (ix in {0, 1, 2}) emit mVMC trans.def with the
    spin-diagonal reduction ``-conj(T) = +1.0`` alone.
    The boundary bond (ix=3, jx=0) additionally acquires the AP wrap
    sign (-1), so its emitted trans value is -1.0. Matches the sign
    pattern StdFace bakes into vmcdry's trans.def under phase0 = 180
    (see tests/validation/uhfk_mvmc_pairproduct/case_apbc/).
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(src, [(1, 0, 0, 1, 1, -1.0, 0.0)])
    emit_trans_def(
        str(src), [4, 1, 1], str(out),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 4
    # Row tuple layout after mapping: (i_site, s_tgt, j_site, s_src, re, im).
    # Key by source flat index.
    by_source = {r[0]: r for r in rows}
    # ix=0 -> jx=1 (interior): trans = -T_hwave * wrap(=1) = +1
    assert by_source[0][2] == 1
    assert by_source[0][4] == pytest.approx(+1.0)
    # ix=1 -> jx=2 (interior)
    assert by_source[1][2] == 2
    assert by_source[1][4] == pytest.approx(+1.0)
    # ix=2 -> jx=3 (interior)
    assert by_source[2][2] == 3
    assert by_source[2][4] == pytest.approx(+1.0)
    # ix=3 -> jx=0 (boundary, wraps once): trans = -T_hwave * wrap(=-1) = -1
    assert by_source[3][2] == 0
    assert by_source[3][4] == pytest.approx(-1.0)


def test_emit_trans_def_apbc_flips_sign_for_boundary_crossing_negative_R(
    tmp_path,
):
    """R = (-1, 0, 0) hop under APBC(x). The boundary crossing is at
    ix=0 (jx = -1 wraps to Lx-1 = 3). Every other source lands
    interior. Guards the negative-wrap branch of the wrap-count phase.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(src, [(-1, 0, 0, 1, 1, -1.0, 0.0)])
    emit_trans_def(
        str(src), [4, 1, 1], str(out),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 4
    # Row tuple layout after mapping: (i_site, s_tgt, j_site, s_src, re, im).
    by_source = {r[0]: r for r in rows}
    # ix=0 -> jx=3 (boundary, wraps once negatively): trans = -T * (-1) = -1
    assert by_source[0][2] == 3
    assert by_source[0][4] == pytest.approx(-1.0)
    # ix=1 -> jx=0 (interior)
    assert by_source[1][2] == 0
    assert by_source[1][4] == pytest.approx(+1.0)
    # ix=2 -> jx=1 (interior)
    assert by_source[2][2] == 1
    assert by_source[2][4] == pytest.approx(+1.0)
    # ix=3 -> jx=2 (interior)
    assert by_source[3][2] == 2
    assert by_source[3][4] == pytest.approx(+1.0)


def test_emit_trans_def_apbc_flips_sign_on_rashba_off_diagonal(tmp_path):
    """Rashba s != t entries also acquire the AP wrap sign at boundary
    crossings. Under APBC(x), a hop with val = +0.5i (up -> dn) at
    R = (1, 0, 0) becomes -0.5i on the boundary-crossing source
    (ix = Lx-1 -> jx = 0). The general rule swaps it to dn -> up and
    gives ``-conj(+0.5i) = +0.5i`` before the wrap phase is applied.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(src, [(1, 0, 0, 1, 2, 0.0, 0.5)])
    emit_trans_def(
        str(src), [4, 1, 1], str(out),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 4
    # Row tuple layout after mapping: (i_site, s_tgt, j_site, s_src, re, im).
    by_source = {r[0]: r for r in rows}
    # ix=0..2 (interior): (re, im) = (0, +0.5)
    for ix in (0, 1, 2):
        _, s, _, t, re, im = by_source[ix]
        assert s == 1 and t == 0
        assert re == pytest.approx(0.0)
        assert im == pytest.approx(+0.5)
    # ix=3 (boundary): sign-flipped to (0, -0.5)
    _, s, _, t, re, im = by_source[3]
    assert s == 1 and t == 0
    assert re == pytest.approx(0.0)
    assert im == pytest.approx(-0.5)


def test_emit_trans_def_apbc_preserves_hermitian(tmp_path):
    """A Hermitian Transfer.dat pair (T[R] = conj(T[-R]) for the same
    orbital pair) stays Hermitian after the emitter applies wrap
    phases under APBC. Concretely for R = (1, 0, 0), val = +0.5i and
    R = (-1, 0, 0), val = -0.5i (the Rashba x-hop pair), the emitted
    mVMC trans.def rows must satisfy
    ``trans(i, s, j, t, re, im) = trans(j, t, i, s, re, -im)`` on
    every (i, j) pair — the trans.def encoding of Hermiticity.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(src, [
        (1, 0, 0, 1, 2, 0.0, 0.5),   # up -> dn, R = +1, val = +0.5i
        (-1, 0, 0, 2, 1, 0.0, -0.5),  # dn -> up, R = -1, val = -0.5i
    ])
    emit_trans_def(
        str(src), [4, 1, 1], str(out),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    rows = _read_trans_def_rows(str(out))
    # Build a dict keyed on (i, s, j, t).
    by_key = {(r[0], r[1], r[2], r[3]): (r[4], r[5]) for r in rows}
    # Every (i, s, j, t) row must have a conjugate partner at (j, t, i, s).
    for (i, s, j, t), (re, im) in by_key.items():
        assert (j, t, i, s) in by_key, (
            f"missing conjugate partner for ({i},{s},{j},{t})"
        )
        conj_re, conj_im = by_key[(j, t, i, s)]
        assert conj_re == pytest.approx(re)
        assert conj_im == pytest.approx(-im)


def test_emit_trans_def_apbc_matches_case_apbc_pattern(tmp_path):
    """1D L=8 chain, single spin-diagonal NN hop T=-1.0 under APBC(x).
    The 16 emitted rows must match case_apbc's vmcdry-produced
    trans.def pattern: interior bonds carry +1.0, boundary bonds
    (jumping across the L-1 -> 0 seam in either direction) carry
    -1.0. Ground-truth check that the bridge's APBC output is
    interchangeable with StdFace's phase0 = 180 output.
    """
    src = tmp_path / "Transfer.dat"
    out = tmp_path / "trans.def"
    _write_transfer(src, [
        (1, 0, 0, 1, 1, -1.0, 0.0),
        (-1, 0, 0, 1, 1, -1.0, 0.0),
    ])
    emit_trans_def(
        str(src), [8, 1, 1], str(out),
        boundary_theta=(np.pi, 0.0, 0.0),
    )
    rows = _read_trans_def_rows(str(out))
    assert len(rows) == 16
    boundary_pairs = {(0, 7), (7, 0)}
    interior_count = 0
    boundary_count = 0
    for i, s, j, t, re, im in rows:
        assert s == 0 and t == 0
        assert im == pytest.approx(0.0)
        if (i, j) in boundary_pairs:
            assert re == pytest.approx(-1.0)
            boundary_count += 1
        else:
            assert re == pytest.approx(+1.0)
            interior_count += 1
    assert boundary_count == 2
    assert interior_count == 14
