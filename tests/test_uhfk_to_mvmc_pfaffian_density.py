"""Deterministic Pfaffian-emulating density check for the bridge's F.

The existing :file:`tests/test_uhfk_to_mvmc_density_check.py` validates
the amplitudes A_up / A_down ahead of orbitalidx aggregation.
mVMC does NOT consume A_up / A_down; it reads ``zqp_orbital_uhfk.dat``,
multiplies each value by the sign column of ``orbitalidx.def`` to
reconstruct F_ij, and evaluates a Pfaffian on the resulting matrix.

This test exercises that exact roundtrip:
1. Build a tiny H-wave / orbitalidx fixture in memory.
2. Run the bridge end-to-end to produce ``zqp_orbital_uhfk.dat``.
3. Reconstruct ``F_mvmc[i, j] = params[idx] * sign`` as mVMC's
   ``slater.c:85`` does.
4. Enumerate all fixed-particle configurations, sum ``|det F_mvmc[ups,
   downs]|^2``, and confirm the per-site density matches the H-wave
   ``greenone.dat`` reference to 1e-10.

The test isolates the bridge's correctness from any mVMC projection
quirks (NSPGaussLeg, NMPTrans, QPTrans), so a future regression in F
construction or orbitalidx sign aggregation fails the bridge CI before
the mVMC E2E harness even runs.
"""
from __future__ import annotations

import itertools
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
FIXTURE = os.path.join(
    REPO_ROOT, "tests/validation/apbc_complexuhf/case_1d_hubbard"
)


def _run_hwave_uhfk(work, ncond=2, boundary="periodic", subshape=(1, 1, 1)):
    """Run H-wave on a L=8 1D Hubbard fixture with the bridge's
    occupation key enabled. ``boundary`` is the x-direction boundary
    condition (``"periodic"`` or ``"antiperiodic"``). ``subshape``
    overrides the SubShape line in ``input.toml``. Returns work."""
    fixture = FIXTURE
    for fn in os.listdir(fixture):
        src = os.path.join(fixture, fn)
        if os.path.isfile(src):
            shutil.copy(src, work)

    toml_path = os.path.join(work, "input.toml")
    with open(toml_path) as fp:
        lines = fp.readlines()
    out_lines = []
    saw_occupation = False
    for ln in lines:
        s = ln.strip()
        if s.startswith("BoundaryCondition"):
            out_lines.append(f'  BoundaryCondition = ["{boundary}", "periodic", "periodic"]\n')
        elif s.startswith("Ncond"):
            out_lines.append(f"  Ncond        = {ncond}\n")
        elif s.startswith("EPS"):
            out_lines.append("  EPS          = 14\n")
        elif s.startswith("SubShape"):
            sub_str = ", ".join(str(int(x)) for x in subshape)
            out_lines.append(f"  SubShape     = [{sub_str}]\n")
        else:
            if "occupation" in s and "=" in s:
                saw_occupation = True
            out_lines.append(ln)
    if not saw_occupation:
        out_lines2 = []
        injected = False
        for ln in out_lines:
            out_lines2.append(ln)
            if ln.strip() == "[file.output]" and not injected:
                out_lines2.append('  occupation = "occupation.npz"\n')
                injected = True
        out_lines = out_lines2
    # flag_fock=false in [mode]
    out_lines2 = []
    injected = False
    for ln in out_lines:
        out_lines2.append(ln)
        if ln.strip() == "[mode]" and not injected:
            out_lines2.append("  flag_fock = false\n")
            injected = True
    with open(toml_path, "w") as fp:
        fp.writelines(out_lines2)

    env = os.environ.copy()
    env["PYTHONPATH"] = os.path.join(REPO_ROOT, "src") + ":" + env.get("PYTHONPATH", "")
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
        cwd=work, env=env, capture_output=True, text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"H-wave failed: {result.stderr}")
    return work


def _write_translation_orbitalidx(path, nsite):
    """Write a per-(i, j) unique orbitalidx.def (NOrbitalIdx = nsite^2)
    with ComplexType=1 and no APBC sign column (3-column form)."""
    with open(path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx {nsite * nsite}\n")
        fp.write("ComplexType 1\n")
        fp.write("======================\n")
        fp.write("======================\n")
        for i in range(nsite):
            for j in range(nsite):
                fp.write(f"{i} {j} {i * nsite + j}\n")
        for k in range(nsite * nsite):
            fp.write(f"{k} 1\n")


def _read_orbitalidx_mapping(path):
    mapping = {}
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("==") or line.startswith("NOrbitalIdx") or line.startswith("ComplexType"):
                continue
            toks = line.split()
            if len(toks) >= 4:
                try:
                    i, j, idx, sgn = int(toks[0]), int(toks[1]), int(toks[2]), int(toks[3])
                    mapping[(i, j)] = (idx, sgn)
                except ValueError:
                    continue
            elif len(toks) == 3:
                try:
                    i, j, idx = int(toks[0]), int(toks[1]), int(toks[2])
                    mapping[(i, j)] = (idx, 1)
                except ValueError:
                    continue
    return mapping


def _read_zqp_orbital(path):
    params = {}
    with open(path) as fp:
        for line in fp:
            line = line.strip()
            if not line or line.startswith("==") or line.startswith("NOrbitalIdx"):
                continue
            toks = line.split()
            if len(toks) >= 3:
                try:
                    idx = int(toks[0])
                    params[idx] = float(toks[1]) + 1j * float(toks[2])
                except ValueError:
                    continue
    return params


def _reconstruct_F_mvmc(orbidx_path, zqp_path, nsite):
    mapping = _read_orbitalidx_mapping(orbidx_path)
    params = _read_zqp_orbital(zqp_path)
    F = np.zeros((nsite, nsite), dtype=complex)
    for (i, j), (idx, sgn) in mapping.items():
        F[i, j] = params[idx] * sgn
    return F


def _enumerate_density(F, ne_up, ne_down, nsite):
    up_cfgs = list(itertools.combinations(range(nsite), ne_up))
    dn_cfgs = list(itertools.combinations(range(nsite), ne_down))
    norm = 0.0
    occ_up = np.zeros(nsite)
    occ_dn = np.zeros(nsite)
    for ups in up_cfgs:
        F_rows = F[list(ups), :]
        for dns in dn_cfgs:
            amp = np.linalg.det(F_rows[:, list(dns)])
            w = abs(amp) ** 2
            if w < 1e-30:
                continue
            norm += w
            for i in ups:
                occ_up[i] += w
            for j in dns:
                occ_dn[j] += w
    if norm > 0:
        occ_up /= norm
        occ_dn /= norm
    return occ_up, occ_dn, norm


def _parse_hwave_diag(path):
    diag = {}
    with open(path) as fp:
        for line in fp:
            toks = line.split()
            if len(toks) < 6:
                continue
            try:
                i, s, j, t = int(toks[0]), int(toks[1]), int(toks[2]), int(toks[3])
                if i == j and s == t:
                    diag[(i, s)] = float(toks[4])
            except ValueError:
                continue
    return diag


def test_bridge_F_reconstructs_correct_density_pbc():
    """The F written by the bridge to zqp_orbital_uhfk.dat must, after
    orbitalidx sign roundtrip, evaluate via Pfaffian (enumeration over
    fixed-N configs) to the same per-site density as H-wave's
    greenone.dat. Catches sign / aggregation regressions independent of
    mVMC's projection options."""
    nsite = 8
    ne_up = 1  # Ne_per_spin = 1 -> k=0 only (non-degenerate, PBC)
    ne_dn = 1
    ne_total = ne_up + ne_dn

    with tempfile.TemporaryDirectory() as tmp:
        # 1. Run H-wave UHFk
        hwave = os.path.join(tmp, "hwave")
        os.makedirs(hwave)
        _run_hwave_uhfk(hwave, ncond=ne_total, boundary="periodic")

        # 2. Write a unique-idx orbitalidx.def (no signs)
        mvmc = os.path.join(tmp, "mvmc")
        os.makedirs(mvmc)
        orbidx_path = os.path.join(mvmc, "orbitalidx.def")
        _write_translation_orbitalidx(orbidx_path, nsite)

        # 3. Run the bridge
        out_path = os.path.join(mvmc, "zqp_orbital_uhfk.dat")
        env = os.environ.copy()
        env["PYTHONPATH"] = REPO_ROOT + ":" + env.get("PYTHONPATH", "")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(hwave, "input.toml"),
             "--eigen", os.path.join(hwave, "output", "eigen.npz"),
             "--occupation", os.path.join(hwave, "output", "occupation.npz"),
             "--geometry", os.path.join(hwave, "geometry_uhf.dat"),
             "--orbitalidx", orbidx_path,
             "--output", out_path,
             "--check-density",
             "--onebodyg-uhf", os.path.join(hwave, "output", "greenone.dat"),
             # Disable rank-lift noise so the Pfaffian density check at
             # tol 1e-10 is not overwhelmed by the default 1e-6 noise.
             "--epsilon-noise", "0"],
            cwd=REPO_ROOT, env=env, capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"bridge failed: stdout={result.stdout}, stderr={result.stderr}"
        )

        # 4. Reconstruct F as mVMC's slater.c does, enumerate density
        F = _reconstruct_F_mvmc(orbidx_path, out_path, nsite)
        occ_up, occ_dn, norm = _enumerate_density(F, ne_up, ne_dn, nsite)

        # 5. Compare against H-wave greenone. The shipped OneBodyG.dat
        # only queries (i=0, j=0..7) for both spins, so we have H-wave
        # ground truth for the site-0 diagonal plus all off-diagonal
        # Green G[0, j].
        diag = _parse_hwave_diag(os.path.join(hwave, "output", "greenone.dat"))
        ref_up0 = diag.get((0, 0))
        ref_dn0 = diag.get((0, 1))
        assert ref_up0 is not None, "H-wave greenone missing (0, 0, 0, 0)"
        assert ref_dn0 is not None, "H-wave greenone missing (0, 1, 0, 1)"

        assert abs(occ_up[0] - ref_up0) < 1e-10, (
            f"site 0 up: Pfaffian density {occ_up[0]:.6e} != "
            f"H-wave {ref_up0:.6e}"
        )
        assert abs(occ_dn[0] - ref_dn0) < 1e-10, (
            f"site 0 down: Pfaffian density {occ_dn[0]:.6e} != "
            f"H-wave {ref_dn0:.6e}"
        )
        # All other sites should be translation-equivalent to site 0
        # (paramagnetic homogeneous Slater on a translation-symmetric
        # 1D chain). Fail if any site density drifts.
        for i in range(nsite):
            assert abs(occ_up[i] - ref_up0) < 1e-10, (
                f"site {i} up: density {occ_up[i]:.6e} not translation-"
                f"equivalent to site 0 ({ref_up0:.6e}); bridge sign "
                f"aggregation is non-uniform"
            )
            assert abs(occ_dn[i] - ref_dn0) < 1e-10, (
                f"site {i} down: density {occ_dn[i]:.6e} not "
                f"translation-equivalent to site 0 ({ref_dn0:.6e})"
            )


def test_bridge_F_reconstructs_correct_density_apbc():
    """Same Pfaffian roundtrip test for APBC L=8, Ne=4 (Ne_per_spin=2,
    the inner (k=+/-pi/8) shell completely filled).

    This is the exact configuration where mVMC's actual <H> evaluation
    misbehaves (see ``tests/validation/uhfk_mvmc_pairproduct/README.md``
    case_apbc note), so locking in that the BRIDGE's F is still
    mathematically correct keeps the bug surface bounded to mVMC's
    projection options."""
    nsite = 8
    ne_up = 2  # closed shell at k=+/-pi/8
    ne_dn = 2
    ne_total = ne_up + ne_dn

    with tempfile.TemporaryDirectory() as tmp:
        hwave = os.path.join(tmp, "hwave")
        os.makedirs(hwave)
        _run_hwave_uhfk(hwave, ncond=ne_total, boundary="antiperiodic")

        mvmc = os.path.join(tmp, "mvmc")
        os.makedirs(mvmc)
        orbidx_path = os.path.join(mvmc, "orbitalidx.def")
        # APBC orbitalidx needs the 4th sign column. Bridge fail-fasts
        # without it; build a translation-symmetric one with all signs +1
        # (the bridge will absorb the boundary phase into F values).
        with open(orbidx_path, "w") as fp:
            fp.write("======================\n")
            fp.write(f"NOrbitalIdx {nsite * nsite}\n")
            fp.write("ComplexType 1\n")
            fp.write("======================\n")
            fp.write("======================\n")
            for i in range(nsite):
                for j in range(nsite):
                    fp.write(f"{i} {j} {i * nsite + j}  1\n")
            for k in range(nsite * nsite):
                fp.write(f"{k} 1\n")

        out_path = os.path.join(mvmc, "zqp_orbital_uhfk.dat")
        env = os.environ.copy()
        env["PYTHONPATH"] = REPO_ROOT + ":" + env.get("PYTHONPATH", "")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(hwave, "input.toml"),
             "--eigen", os.path.join(hwave, "output", "eigen.npz"),
             "--occupation", os.path.join(hwave, "output", "occupation.npz"),
             "--geometry", os.path.join(hwave, "geometry_uhf.dat"),
             "--orbitalidx", orbidx_path,
             "--output", out_path,
             "--check-density",
             "--onebodyg-uhf", os.path.join(hwave, "output", "greenone.dat"),
             # Disable rank-lift noise so the Pfaffian density check at
             # tol 1e-10 is not overwhelmed by the default 1e-6 noise.
             "--epsilon-noise", "0"],
            cwd=REPO_ROOT, env=env, capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"bridge failed: stdout={result.stdout}, stderr={result.stderr}"
        )

        F = _reconstruct_F_mvmc(orbidx_path, out_path, nsite)
        occ_up, occ_dn, norm = _enumerate_density(F, ne_up, ne_dn, nsite)

        diag = _parse_hwave_diag(os.path.join(hwave, "output", "greenone.dat"))
        ref_up0 = diag[(0, 0)]
        ref_dn0 = diag[(0, 1)]
        assert abs(occ_up[0] - ref_up0) < 1e-10
        assert abs(occ_dn[0] - ref_dn0) < 1e-10
        for i in range(nsite):
            assert abs(occ_up[i] - ref_up0) < 1e-10, (
                f"site {i} up: density {occ_up[i]:.6e} not translation-"
                f"equivalent to site 0 ({ref_up0:.6e}) under APBC"
            )


def test_bridge_F_reconstructs_correct_density_apbc_subshape_2():
    """Same Pfaffian roundtrip test for APBC L=8 with SubShape=[2,1,1]
    (folded lattice has 4 cells x 2 sublattice = 8 physical sites).
    Ne=4 gives a closed shell."""
    nsite = 8
    ne_up = 2
    ne_dn = 2
    ne_total = ne_up + ne_dn

    with tempfile.TemporaryDirectory() as tmp:
        hwave = os.path.join(tmp, "hwave")
        os.makedirs(hwave)
        _run_hwave_uhfk(hwave, ncond=ne_total, boundary="antiperiodic",
                         subshape=(2, 1, 1))

        mvmc = os.path.join(tmp, "mvmc")
        os.makedirs(mvmc)
        orbidx_path = os.path.join(mvmc, "orbitalidx.def")
        with open(orbidx_path, "w") as fp:
            fp.write("======================\n")
            fp.write(f"NOrbitalIdx {nsite * nsite}\n")
            fp.write("ComplexType 1\n")
            fp.write("======================\n")
            fp.write("======================\n")
            for i in range(nsite):
                for j in range(nsite):
                    fp.write(f"{i} {j} {i * nsite + j}  1\n")
            for k in range(nsite * nsite):
                fp.write(f"{k} 1\n")

        out_path = os.path.join(mvmc, "zqp_orbital_uhfk.dat")
        env = os.environ.copy()
        env["PYTHONPATH"] = REPO_ROOT + ":" + env.get("PYTHONPATH", "")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", os.path.join(hwave, "input.toml"),
             "--eigen", os.path.join(hwave, "output", "eigen.npz"),
             "--occupation", os.path.join(hwave, "output", "occupation.npz"),
             "--geometry", os.path.join(hwave, "geometry_uhf.dat"),
             "--orbitalidx", orbidx_path,
             "--output", out_path,
             "--check-density",
             "--onebodyg-uhf", os.path.join(hwave, "output", "greenone.dat"),
             "--epsilon-noise", "0"],
            cwd=REPO_ROOT, env=env, capture_output=True, text=True,
        )
        assert result.returncode == 0, (
            f"bridge failed: stdout={result.stdout}, stderr={result.stderr}"
        )

        F = _reconstruct_F_mvmc(orbidx_path, out_path, nsite)
        occ_up, occ_dn, norm = _enumerate_density(F, ne_up, ne_dn, nsite)

        diag = _parse_hwave_diag(os.path.join(hwave, "output", "greenone.dat"))
        ref_up0 = diag[(0, 0)]
        ref_dn0 = diag[(0, 1)]
        assert abs(occ_up[0] - ref_up0) < 1e-10
        assert abs(occ_dn[0] - ref_dn0) < 1e-10
        for i in range(nsite):
            assert abs(occ_up[i] - ref_up0) < 1e-10
