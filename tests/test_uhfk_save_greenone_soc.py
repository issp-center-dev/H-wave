"""SOC-mode _save_greenone: (i, s, j, t) rows for all (s, t) combinations.

Two coverage points:

1. ``test_save_greenone_soc_emits_all_st_combinations`` verifies that under
   ``enable_spin_orbital = true`` the greenone writer emits rows for every
   ``(s, t) in {0, 1}^2`` that ``OneBodyG.dat`` requests.
2. ``test_save_greenone_backward_compat_non_soc`` locks in the non-SOC
   path: only ``s == t`` rows are emitted.

Both tests drive H-wave end-to-end through a temporary input directory so
the SOC packing convention (``2 * a + s``) is exercised as it appears in
production.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent


def _write_geometry_dat(path: Path, norbit: int) -> None:
    """Write ``Geometry.dat``: primitive vectors + Norbit + zero positions.

    ``norbit`` is the count H-wave reads from ``Geometry.dat`` (the total
    spin-orbital count under SOC, otherwise the physical orbital count).
    """
    lines = [
        "  1.000000000000   0.000000000000   0.000000000000\n",
        "  0.000000000000   1.000000000000   0.000000000000\n",
        "  0.000000000000   0.000000000000   1.000000000000\n",
        "{}\n".format(norbit),
    ]
    for _ in range(norbit):
        lines.append(
            "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
        )
    path.write_text("".join(lines))


def _write_geometry_uhf_dat(path: Path, nsites: int) -> None:
    """Write ``geometry_uhf.dat`` for a 1D ``nsites`` chain, one orbital
    per unit cell. Cell shape is ``(nsites, 1, 1)``."""
    lines = [
        "  1.000000000000   0.000000000000   0.000000000000\n",
        "  0.000000000000   1.000000000000   0.000000000000\n",
        "  0.000000000000   0.000000000000   1.000000000000\n",
        "  0.000000000000   0.000000000000   0.000000000000\n",
        "  {} 0 0\n".format(nsites),
        "  0 1 0\n",
        "  0 0 1\n",
    ]
    for i in range(nsites):
        lines.append("  {} 0 0 0\n".format(i))
    path.write_text("".join(lines))


def _write_wannier90_like(path: Path, header: str, entries):
    """Write a Wannier90-like interaction file consumed by ``read_w90``.

    ``entries`` is an iterable of ``(irx, iry, irz, iorb1_1based,
    iorb2_1based, re, im)`` tuples. Each entry has degeneracy 1 (the
    all-ones ndegen block puts ``read_w90`` on the fast path).
    """
    entries = list(entries)
    n = len(entries)
    lines = [
        "{}\n".format(header),
        "1\n",           # num_wann placeholder (unused by read_w90)
        "{}\n".format(n),
        " ".join(["1"] * n) + "\n",
    ]
    for (irx, iry, irz, o1, o2, re, im) in entries:
        lines.append(
            "  {:4d} {:4d} {:4d} {:4d} {:4d}  {:.12f}  {:.12f}\n".format(
                irx, iry, irz, o1, o2, re, im,
            )
        )
    path.write_text("".join(lines))


def _write_transfer_dat(path: Path, entries) -> None:
    _write_wannier90_like(
        path, "Transfer in wannier90-like format for uhfk", entries,
    )


def _write_coulomb_intra_dat(path: Path, norbit_phys: int) -> None:
    """U=0 on every physical orbital of the primitive cell so the SCF
    converges immediately (non-interacting reference)."""
    entries = [(0, 0, 0, o + 1, o + 1, 0.0, 0.0) for o in range(norbit_phys)]
    _write_wannier90_like(
        path, "CoulombIntra in wannier90-like format for uhfk", entries,
    )


def _write_onebodyg_dat(path: Path, entries) -> None:
    """Write ``OneBodyG.dat`` listing ``(i, s, j, t)`` requests."""
    entries = list(entries)
    lines = [
        "===============================\n",
        "NCisAjs         {}\n".format(len(entries)),
        "===============================\n",
        "======== Green functions ======\n",
        "===============================\n",
    ]
    for (i, s, j, t) in entries:
        lines.append("    {}     {}     {}     {}\n".format(i, s, j, t))
    path.write_text("".join(lines))


def _write_input_toml(path: Path, *, enable_spin_orbital: bool) -> None:
    """Minimal 2-site 1D chain input.toml. The SOC flag is toggled via
    ``enable_spin_orbital``; everything else is identical."""
    lines = [
        "[log]\n",
        "  print_level = 1\n",
        "  print_step  = 20\n",
        "\n",
        "[mode]\n",
        '  mode = "UHFk"\n',
        "  flag_fock = false\n",
    ]
    if enable_spin_orbital:
        lines.append("  enable_spin_orbital = true\n")
    lines.extend([
        "\n",
        "[mode.param]\n",
        "  Ncond        = 2\n",
        "  IterationMax = 100\n",
        "  EPS          = 14\n",
        "  Mix          = 0.5\n",
        "  RndSeed      = 12345\n",
        "  T            = 0.0\n",
        "  CellShape    = [3, 1, 1]\n",
        "  SubShape     = [1, 1, 1]\n",
        '  BoundaryCondition = ["periodic", "periodic", "periodic"]\n',
        "\n",
        "[file]\n",
        "[file.input]\n",
        '  path_to_input = ""\n',
        '  geometry_uhf = "geometry_uhf.dat"\n',
        '  onebodyg_uhf = "OneBodyG.dat"\n',
        "[file.input.interaction]\n",
        '  path_to_input = "./"\n',
        '  Geometry     = "Geometry.dat"\n',
        '  Transfer     = "Transfer.dat"\n',
        '  CoulombIntra = "CoulombIntra.dat"\n',
        "\n",
        "[file.output]\n",
        '  path_to_output = "output"\n',
        '  energy   = "energy.dat"\n',
        '  eigen    = "eigen.npz"\n',
        '  green    = "green.npz"\n',
        '  onebodyg = "greenone.dat"\n',
    ])
    path.write_text("".join(lines))


def _write_minimal_soc_input(tmp_path: Path) -> Path:
    """SOC-mode fixture: 3-site 1D chain, 1 physical orbital / cell,
    Norbit=2 (two spin-orbitals per cell), NN hopping ``t = 1`` diagonal
    in spin, U=0. ``CellShape = [3, 1, 1]`` so both ``R = +ex`` and ``R
    = -ex`` fit within the range check (width 3 <= 3). OneBodyG.dat
    requests every ``(i, s, j, t)`` row on the ``i = j = 0`` diagonal
    site (4 rows spanning all ``(s, t) in {0, 1}^2``)."""
    _write_input_toml(tmp_path / "input.toml", enable_spin_orbital=True)
    _write_geometry_dat(tmp_path / "Geometry.dat", norbit=2)
    _write_geometry_uhf_dat(tmp_path / "geometry_uhf.dat", nsites=3)
    # Wannier90 SOI ordering: file index = 2 * orb + spin + 1.
    # For orb=0 the two spin-orbitals are file indices 1 (up) and 2 (down).
    # NN hopping H = -t sum_{i, s} (c_{i+1, s}^dagger c_{i, s} + h.c.),
    # diagonal in spin, encoded as sigma matrix diag(-t, -t).
    _write_transfer_dat(tmp_path / "Transfer.dat", [
        (+1, 0, 0, 1, 1, -1.0, 0.0),   # up-up, R = +ex
        (-1, 0, 0, 1, 1, -1.0, 0.0),   # up-up, R = -ex (hermitian conjugate)
        (+1, 0, 0, 2, 2, -1.0, 0.0),   # down-down, R = +ex
        (-1, 0, 0, 2, 2, -1.0, 0.0),   # down-down, R = -ex
    ])
    _write_coulomb_intra_dat(tmp_path / "CoulombIntra.dat", norbit_phys=1)
    # Cover every (s, t) combination on the on-site block at i=j=0.
    entries = [(0, s, 0, t) for s in range(2) for t in range(2)]
    _write_onebodyg_dat(tmp_path / "OneBodyG.dat", entries)
    return tmp_path


def _write_minimal_non_soc_input(tmp_path: Path) -> Path:
    """Non-SOC fixture: 3-site 1D chain, 1 physical orbital / cell,
    Norbit=1, NN hopping ``t = 1``, U=0. OneBodyG.dat requests only the
    ``s == t`` rows."""
    _write_input_toml(tmp_path / "input.toml", enable_spin_orbital=False)
    _write_geometry_dat(tmp_path / "Geometry.dat", norbit=1)
    _write_geometry_uhf_dat(tmp_path / "geometry_uhf.dat", nsites=3)
    _write_transfer_dat(tmp_path / "Transfer.dat", [
        (+1, 0, 0, 1, 1, -1.0, 0.0),
        (-1, 0, 0, 1, 1, -1.0, 0.0),
    ])
    _write_coulomb_intra_dat(tmp_path / "CoulombIntra.dat", norbit_phys=1)
    # Only diagonal (s == t) rows on the i=j=0 block.
    entries = [(0, s, 0, s) for s in range(2)]
    _write_onebodyg_dat(tmp_path / "OneBodyG.dat", entries)
    return tmp_path


def _run_hwave(input_dir: Path) -> None:
    """Drive ``hwave.qlms.main`` in a subprocess. Raises with combined
    stdout/stderr if the run fails so the test message is actionable."""
    env = os.environ.copy()
    env["PYTHONPATH"] = str(REPO_ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    script = (
        "import numpy as np\n"
        "if not hasattr(np, 'float'):\n"
        "    np.float = float\n"
        "import sys\n"
        "sys.argv = ['hwave', 'input.toml']\n"
        "from hwave.qlms import main\n"
        "main()\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=str(input_dir), env=env, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "H-wave failed (rc={rc}):\nstdout:\n{out}\nstderr:\n{err}".format(
                rc=result.returncode, out=result.stdout, err=result.stderr,
            )
        )


def _parse_greenone(path: Path):
    """Parse ``greenone.dat`` into ``(i, s, j, t, re, im)`` tuples."""
    if not path.exists():
        raise FileNotFoundError(str(path))
    rows = []
    with path.open() as fp:
        for line in fp:
            toks = line.split()
            if len(toks) < 6:
                continue
            try:
                i, s, j, t = (int(toks[0]), int(toks[1]),
                              int(toks[2]), int(toks[3]))
                re, im = float(toks[4]), float(toks[5])
            except ValueError:
                continue
            rows.append((i, s, j, t, re, im))
    return rows


def test_save_greenone_soc_emits_all_st_combinations(tmp_path):
    """SCF with ``enable_spin_orbital = true``: ``_save_greenone`` writes
    rows for every ``(s, t) in {0, 1}^2`` requested in ``OneBodyG.dat``.

    This test exercises the writer through ``UHFk.solve``.
    """
    input_dir = _write_minimal_soc_input(tmp_path)
    _run_hwave(input_dir)
    rows = _parse_greenone(input_dir / "output" / "greenone.dat")
    # OneBodyG.dat lists 4 rows (2^2 = 4 (s, t) pairs at i = j = 0)
    assert len(rows) == 4, "expected 4 rows, got {}".format(len(rows))
    unique_st_pairs = {(int(r[1]), int(r[3])) for r in rows}
    assert unique_st_pairs == {(0, 0), (0, 1), (1, 0), (1, 1)}


def test_save_greenone_backward_compat_non_soc(tmp_path):
    """SCF with ``enable_spin_orbital = false``: existing ``s == t``
    emission path is unchanged (only diagonal spin rows in
    ``greenone.dat``)."""
    input_dir = _write_minimal_non_soc_input(tmp_path)
    _run_hwave(input_dir)
    rows = _parse_greenone(input_dir / "output" / "greenone.dat")
    assert len(rows) == 2, "expected 2 rows, got {}".format(len(rows))
    for r in rows:
        assert int(r[1]) == int(r[3]), (
            "non-SOC greenone.dat has s != t: {}".format(r)
        )
