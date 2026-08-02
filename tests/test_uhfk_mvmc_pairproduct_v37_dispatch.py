"""CLI allowlist dispatch tests.

Pins the CLI's shift from a blanket 'multi-direction APBC rejected'
message to the frozenset-based allowlist that admits the supported lattice
shapes while rejecting novel combinations.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest


_REPO = Path(__file__).resolve().parents[1]
_CLI = str(_REPO / "tools" / "uhfk_to_mvmc.py")


def _prepare_stub_workspace(tmp_path, boundary_condition, cell_shape,
                             sub_shape):
    """Write a minimal input.toml + stub geometry / transfer / coulomb
    files sufficient to reach the allowlist check. The CLI will still
    fail downstream on missing eigen.npz etc., which is fine — we only
    care whether the reject predicate fires."""
    bc_toml = ", ".join(f'"{s}"' for s in boundary_condition)
    (tmp_path / "input.toml").write_text(
        f"""
[log]
print_level = 0

[mode]
mode = "UHFk"
enable_spin_orbital = true

[mode.param]
Ncond = 4
CellShape = {list(cell_shape)}
SubShape = {list(sub_shape)}
BoundaryCondition = [{bc_toml}]
IterationMax = 1
EPS = 8
Mix = 0.5
RndSeed = 12345
T = 0.0

[file]
[file.input]
path_to_input = ""
geometry_uhf = "geometry_uhf.dat"
[file.input.interaction]
path_to_input = "./"
Geometry     = "Geometry.dat"
Transfer     = "Transfer.dat"
CoulombIntra = "CoulombIntra.dat"

[file.output]
path_to_output = "output"
"""
    )
    # Minimal 1-site geometry stubs; enough content that early parsers
    # do not choke before reaching the allowlist check.
    Nsites_x, Nsites_y, Nsites_z = cell_shape
    site_lines = "\n".join(
        f"  {x} {y} {z} 0"
        for z in range(Nsites_z)
        for y in range(Nsites_y)
        for x in range(Nsites_x)
    )
    (tmp_path / "geometry_uhf.dat").write_text(
        "  1.000000000000   0.000000000000   0.000000000000\n"
        "  0.000000000000   1.000000000000   0.000000000000\n"
        "  0.000000000000   0.000000000000   1.000000000000\n"
        "  0.000000000000   0.000000000000   0.000000000000\n"
        f"  {Nsites_x} 0 0\n  0 {Nsites_y} 0\n  0 0 {Nsites_z}\n"
        + site_lines + "\n"
    )
    (tmp_path / "Geometry.dat").write_text(
        "  1.000000000000   0.000000000000   0.000000000000\n"
        "  0.000000000000   1.000000000000   0.000000000000\n"
        "  0.000000000000   0.000000000000   1.000000000000\n"
        "  0.000000000000   0.000000000000   0.000000000000\n"
        f"  {Nsites_x} 0 0\n  0 {Nsites_y} 0\n  0 0 {Nsites_z}\n"
        "1\n"
        "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
    )
    (tmp_path / "Transfer.dat").write_text("Transfer\n1\n0\n1\n1\n")
    (tmp_path / "CoulombIntra.dat").write_text("CoulombIntra\n1\n0\n1\n1\n")


def _run_cli(tmp_path):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(_REPO / "src") + os.pathsep + str(_REPO)
    return subprocess.run(
        [sys.executable, _CLI,
         "--input", str(tmp_path / "input.toml"),
         "--eigen", str(tmp_path / "eigen.npz"),
         "--occupation", str(tmp_path / "occupation.npz"),
         "--geometry", str(tmp_path / "geometry_uhf.dat"),
         "--orbitalidx", str(tmp_path / "orbitalidx.def"),
         "--output", str(tmp_path / "zqp.dat"),
         "--no-check-density"],
        capture_output=True, text=True, env=env,
    )


def test_v37_single_direction_apbc_on_new_lattice_rejected(tmp_path):
    """[4,4,4]/[2,2,2] with single-direction APBC (e.g., x-only) MUST
    be rejected: the allowlist covers only xy/xz/yz/xyz masks on this
    lattice, while (1,0,0) is allowed on [6,4,1]/[2,2,1]."""
    _prepare_stub_workspace(
        tmp_path,
        boundary_condition=["antiperiodic", "periodic", "periodic"],
        cell_shape=(4, 4, 4), sub_shape=(2, 2, 2),
    )
    res = _run_cli(tmp_path)
    assert res.returncode == 2, res.stderr
    assert "unsupported SOC + APBC + SubShape combination" in res.stderr, res.stderr


def test_v37_xy_apbc_on_new_lattice_accepted_by_predicate(tmp_path):
    """[4,4,4]/[2,2,2] with AP-AP-P MUST NOT trigger the allowlist
    reject; the CLI may still fail later on missing eigen.npz."""
    _prepare_stub_workspace(
        tmp_path,
        boundary_condition=["antiperiodic", "antiperiodic", "periodic"],
        cell_shape=(4, 4, 4), sub_shape=(2, 2, 2),
    )
    res = _run_cli(tmp_path)
    assert "unsupported SOC + APBC + SubShape combination" not in res.stderr, res.stderr


def test_v36_shipping_fixture_still_accepted(tmp_path):
    """[6,4,1]/[2,2,1] with AP-P-P (case_soc_rashba_2d_sub_apbc)
    MUST remain accepted."""
    _prepare_stub_workspace(
        tmp_path,
        boundary_condition=["antiperiodic", "periodic", "periodic"],
        cell_shape=(6, 4, 1), sub_shape=(2, 2, 1),
    )
    res = _run_cli(tmp_path)
    assert "unsupported SOC + APBC + SubShape combination" not in res.stderr, res.stderr
