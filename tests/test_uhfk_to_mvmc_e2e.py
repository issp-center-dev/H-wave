"""End-to-end bridge tests: H-wave UHFk SCF → bridge → density check.
Uses the existing case_1d_hubbard fixture from apbc_complexuhf.
"""
from __future__ import annotations

import os, sys, shutil, subprocess, tempfile
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
FIXTURE = os.path.join(
    REPO_ROOT, "tests/validation/apbc_complexuhf/case_1d_hubbard"
)


def _run_hwave(case_dir, boundary, ncond, subshape=(1, 1, 1)):
    """Copy case_dir, patch input.toml, run H-wave UHFk. Return work dir.

    ``subshape`` overrides the SubShape line in input.toml when > [1,1,1].
    """
    tmp = tempfile.mkdtemp(prefix="hwave_bridge_e2e_")
    for fn in os.listdir(case_dir):
        src = os.path.join(case_dir, fn)
        if os.path.isfile(src):
            shutil.copy(src, tmp)
    toml_path = os.path.join(tmp, "input.toml")
    with open(toml_path) as fp:
        lines = fp.readlines()
    new_lines = []
    saw_flag_fock = False
    saw_occupation = False
    for ln in lines:
        s = ln.strip()
        if s.startswith("BoundaryCondition"):
            new_lines.append(
                '  BoundaryCondition = ["'
                + '", "'.join(boundary)
                + '"]\n'
            )
        elif s.startswith("Ncond"):
            new_lines.append(f"  Ncond = {ncond}\n")
        elif s.startswith("EPS"):
            # Tighten SCF convergence so the density check 1e-10 tol holds.
            new_lines.append("  EPS          = 14\n")
        elif s.startswith("SubShape"):
            sub_str = ", ".join(str(int(x)) for x in subshape)
            new_lines.append(f"  SubShape     = [{sub_str}]\n")
        elif s.startswith("flag_fock"):
            new_lines.append("  flag_fock = false\n")
            saw_flag_fock = True
        else:
            if "occupation" in s and "=" in s:
                saw_occupation = True
            new_lines.append(ln)
    # Inject flag_fock = false in [mode] (NOT [mode.param]). H-wave's
    # uhfk.py:82 reads ``flag_fock`` from ``info_mode`` directly, i.e.
    # the outer [mode] section, per docs/.../uhfk/.../index_config.rst.
    # Putting it under [mode.param] is silently ignored.
    if not saw_flag_fock:
        out = []
        injected = False
        for ln in new_lines:
            if ln.strip().startswith("[mode]") and not injected:
                out.append(ln)
                out.append("  flag_fock = false\n")
                injected = True
                continue
            out.append(ln)
        if not injected:
            out.insert(0, "[mode]\n  flag_fock = false\n")
        new_lines = out
    # Add occupation under [file.output] if not present
    if not saw_occupation:
        out2 = []
        in_fileout = False
        injected = False
        for ln in new_lines:
            out2.append(ln)
            if ln.strip() == "[file.output]":
                in_fileout = True
            elif in_fileout and ln.strip().startswith("[") and not injected:
                out2.insert(-1, '  occupation = "occupation.npz"\n')
                injected = True
                in_fileout = False
        if not injected:
            out2.append('  occupation = "occupation.npz"\n')
        new_lines = out2
    with open(toml_path, "w") as fp:
        fp.writelines(new_lines)

    # hwave.qlms has no __main__ entry, and qlmsio.wan90.read_geometry uses
    # the legacy ``np.float`` alias removed in numpy >= 1.20. Install a
    # process-local shim and call main() directly, matching the recipe in
    # tests/validation/apbc_complexuhf/run.sh.
    runner = (
        "import numpy as np\n"
        "if not hasattr(np, 'float'):\n"
        "    np.float = float\n"
        "import sys\n"
        "sys.argv = ['hwave', 'input.toml']\n"
        "from hwave.qlms import main\n"
        "main()\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", runner],
        cwd=tmp, capture_output=True, text=True,
    )
    if result.returncode != 0:
        # leave tmp in place to aid post-mortem
        raise RuntimeError(
            f"hwave.qlms failed in {tmp}: stdout={result.stdout} "
            f"stderr={result.stderr}"
        )
    return tmp


def _write_orbitalidx_pbc(path, L):
    with open(path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx {L}\n")
        fp.write("ComplexType 1\n")
        fp.write("======================\n")
        fp.write("======================\n")
        for i in range(L):
            for j in range(L):
                fp.write(f"{i} {j} {(j - i) % L}\n")
        for k in range(L):
            fp.write(f"{k} 1\n")


def test_e2e_pbc_case_1d_hubbard_density_check_passes():
    """PBC L=8, Ncond=2, flag_fock=False: bridge → density check OK."""
    work = _run_hwave(
        FIXTURE,
        boundary=("periodic", "periodic", "periodic"),
        ncond=2,
    )
    try:
        L = 8
        orbitalidx_path = os.path.join(work, "orbitalidx.def")
        _write_orbitalidx_pbc(orbitalidx_path, L)

        onebodyg_path = os.path.join(work, "output", "greenone.dat")
        assert os.path.isfile(onebodyg_path), (
            f"missing onebodyg output at {onebodyg_path}"
        )
        out = os.path.join(work, "zqp_orbital_uhfk.dat")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(work, "input.toml"),
             "--eigen", os.path.join(work, "output", "eigen.npz"),
             "--occupation", os.path.join(work, "output", "occupation.npz"),
             "--geometry", os.path.join(work, "geometry_uhf.dat"),
             "--orbitalidx", orbitalidx_path,
             "--output", out,
             "--check-density",
             "--onebodyg-uhf", onebodyg_path],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode == 0, (
            f"bridge failed: stdout={result.stdout} stderr={result.stderr}"
        )
        assert os.path.isfile(out)
        assert "density check OK" in result.stdout
    finally:
        shutil.rmtree(work)


def _write_orbitalidx_apbc(path, L):
    """APBC orbitalidx.def: 4-column form with sign = -1 when j < i (the
    pair crosses the boundary cut). Translation idx = (j - i) mod L."""
    with open(path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx {L}\n")
        fp.write("ComplexType 1\n")
        fp.write("======================\n")
        fp.write("======================\n")
        for i in range(L):
            for j in range(L):
                idx = (j - i) % L
                sign = -1 if j < i else 1
                fp.write(f"{i} {j} {idx} {sign:+d}\n")
        for k in range(L):
            fp.write(f"{k} 1\n")


def test_e2e_apbc_case_1d_hubbard_density_check_passes():
    """APBC L=8, Ncond=4, flag_fock=False: bridge → density check OK."""
    work = _run_hwave(
        FIXTURE,
        boundary=("antiperiodic", "periodic", "periodic"),
        ncond=4,
    )
    try:
        L = 8
        orbitalidx_path = os.path.join(work, "orbitalidx.def")
        _write_orbitalidx_apbc(orbitalidx_path, L)

        onebodyg_path = os.path.join(work, "output", "greenone.dat")
        assert os.path.isfile(onebodyg_path)
        out = os.path.join(work, "zqp_orbital_uhfk.dat")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(work, "input.toml"),
             "--eigen", os.path.join(work, "output", "eigen.npz"),
             "--occupation", os.path.join(work, "output", "occupation.npz"),
             "--geometry", os.path.join(work, "geometry_uhf.dat"),
             "--orbitalidx", orbitalidx_path,
             "--output", out,
             "--check-density",
             "--onebodyg-uhf", onebodyg_path],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode == 0, (
            f"bridge failed: stdout={result.stdout} stderr={result.stderr}"
        )
        assert "density check OK" in result.stdout
    finally:
        shutil.rmtree(work)


def test_e2e_apbc_subshape_2_density_check_passes():
    """APBC L=8 SubShape=[2,1,1] Ncond=4: bridge --check-density.

    Exercises the folded-eigen and unfold path against H-wave's
    greenone.dat at 1e-10.
    """
    work = _run_hwave(
        FIXTURE,
        boundary=("antiperiodic", "periodic", "periodic"),
        ncond=4,
        subshape=(2, 1, 1),
    )
    try:
        L = 8
        orbitalidx_path = os.path.join(work, "orbitalidx.def")
        _write_orbitalidx_apbc(orbitalidx_path, L)

        onebodyg_path = os.path.join(work, "output", "greenone.dat")
        out = os.path.join(work, "zqp_orbital_uhfk.dat")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(work, "input.toml"),
             "--eigen", os.path.join(work, "output", "eigen.npz"),
             "--occupation", os.path.join(work, "output", "occupation.npz"),
             "--geometry", os.path.join(work, "geometry_uhf.dat"),
             "--orbitalidx", orbitalidx_path,
             "--output", out,
             "--check-density",
             "--onebodyg-uhf", onebodyg_path,
             "--epsilon-noise", "0"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode == 0, (
            f"stdout={result.stdout} stderr={result.stderr}"
        )
        assert "density check OK" in result.stdout
    finally:
        shutil.rmtree(work)
