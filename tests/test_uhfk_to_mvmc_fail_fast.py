"""Fail-fast tests for the bridge CLI."""
from __future__ import annotations

import sys, os, subprocess, tempfile
sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir))

import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _write_minimal_inputs(tmp, subshape=(1, 1, 1)):
    """Write a minimal valid input set; returns paths dict."""
    paths = {}
    paths["input"] = os.path.join(tmp, "input.toml")
    with open(paths["input"], "w") as fp:
        sub_str = ", ".join(str(int(x)) for x in subshape)
        fp.write(
            "[mode.param]\n"
            "Ncond = 2\n"
            "2Sz = 0\n"
            "T = 0.0\n"
            "CellShape = [4, 1, 1]\n"
            f"SubShape  = [{sub_str}]\n"
            'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
        )
    paths["eigen"] = os.path.join(tmp, "eigen.npz")
    np.savez(paths["eigen"],
        eigenvalue=np.zeros((4, 2), dtype=np.float64),
        eigenvector=np.zeros((4, 2, 2), dtype=np.complex128),
        wavevector_unit=np.eye(3, dtype=np.float64),
        wavevector_index=np.array(
            [[v, 0, 0] for v in [0, 1, -2, -1]], dtype=np.int64
        ),
        twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
    )
    paths["occupation"] = os.path.join(tmp, "occupation.npz")
    np.savez(paths["occupation"],
        occupation=np.zeros((4, 2), dtype=np.float64),
        mu=np.array([0.0, 0.0], dtype=np.float64),
        T=np.float64(0.0),
        column_spin=np.array([0, 1], dtype=np.int64),
        column_mu_group=np.array([0, 1], dtype=np.int64),
    )
    paths["geometry"] = os.path.join(tmp, "geometry_uhf.dat")
    with open(paths["geometry"], "w") as fp:
        # 3 unit_vec lines + 1 degree line + 3 cell_vec lines + Ns site lines.
        fp.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
        fp.write("0.0 0.0 0.0\n")
        fp.write("4 0 0\n0 1 0\n0 0 1\n")
        for i in range(4):
            fp.write(f"{i} 0 0 0\n")
    paths["orbitalidx"] = os.path.join(tmp, "orbitalidx.def")
    with open(paths["orbitalidx"], "w") as fp:
        fp.write("======================\n")
        fp.write("NOrbitalIdx 4\n")
        fp.write("ComplexType 1\n")
        fp.write("======================\n")
        fp.write("======================\n")
        # 4x4 mapping (3-column = PBC)
        for i in range(4):
            for j in range(4):
                fp.write(f"{i} {j} {(j - i) % 4}\n")
        for k in range(4):
            fp.write(f"{k} 1\n")
    paths["output"] = os.path.join(tmp, "zqp_orbital_uhfk.dat")
    return paths


def test_subshape_not_dividing_cellshape_fails_fast():
    """SubShape[d] must divide CellShape[d]. E.g. SubShape=[3,1,1] on
    CellShape=[4,1,1] must be rejected fail-fast."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(3, 1, 1))
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"],
             "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"],
             "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0
        assert "SubShape" in result.stderr and (
            "divide" in result.stderr or "does not divide" in result.stderr
        )


def test_eigenvector_shape_mismatch_fails_fast():
    """If the eigen.npz shape does not match nd = 2 * norb_orig * subvol,
    the CLI must raise before building amplitudes."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(2, 1, 1))
        # Overwrite eigen.npz with the wrong nd (2 instead of 4)
        np.savez(paths["eigen"],
            eigenvalue=np.zeros((2, 2), dtype=np.float64),
            eigenvector=np.zeros((2, 2, 2), dtype=np.complex128),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array([[v, 0, 0] for v in [0, -1]], dtype=np.int64),
            twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
        )
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"],
             "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"],
             "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0
        assert "shape" in result.stderr.lower() or "nd" in result.stderr.lower()


def test_apbc_without_sign_column_fails_fast():
    """APBC + orbitalidx.def lacking 4-column sign → reject.

    The eigen.npz twist_offset is updated in lockstep with the input.toml
    BoundaryCondition so the eigen twist consistency check
    does not preempt the sign-column reject this test is targeting.
    """
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        with open(paths["input"], "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 2\n2Sz = 0\nT = 0.0\n"
                "CellShape = [4, 1, 1]\nSubShape  = [1, 1, 1]\n"
                'BoundaryCondition = ["antiperiodic", "periodic", "periodic"]\n'
            )
        # Rewrite eigen.npz with twist_offset matching the APBC BC above.
        np.savez(
            paths["eigen"],
            eigenvalue=np.zeros((4, 2), dtype=np.float64),
            eigenvector=np.zeros((4, 2, 2), dtype=np.complex128),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array(
                [[v, 0, 0] for v in [0, 1, -2, -1]], dtype=np.int64
            ),
            twist_offset=np.array([0.5, 0.0, 0.0], dtype=np.float64),
        )
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"], "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"], "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0, f"stderr={result.stderr}"
        assert "sign column" in result.stderr or "sign" in result.stderr


def test_multi_orbital_geometry_fails_fast():
    """Geometry with norb > 1 must be rejected."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        # Overwrite geometry with 2 orbitals per cell (2 lines per R)
        with open(paths["geometry"], "w") as fp:
            fp.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
            fp.write("0.0 0.0 0.0\n")
            fp.write("4 0 0\n0 1 0\n0 0 1\n")
            for i in range(4):
                fp.write(f"{i} 0 0 0\n")
                fp.write(f"{i} 0 0 1\n")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"], "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"], "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0, f"stderr={result.stderr}"
        assert "norb" in result.stderr.lower() or "orbital" in result.stderr.lower()


def test_sz_free_column_spin_fails_fast():
    """column_spin = -1 in occupation.npz → reject."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        np.savez(paths["occupation"],
            occupation=np.zeros((4, 2), dtype=np.float64),
            mu=np.array([0.0], dtype=np.float64),
            T=np.float64(0.0),
            column_spin=np.array([-1, -1], dtype=np.int64),
            column_mu_group=np.array([0, 0], dtype=np.int64),
        )
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"], "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"], "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0, f"stderr={result.stderr}"
        assert "Sz-free" in result.stderr or "column_spin" in result.stderr


def test_cli_failure_leaves_no_output_artifact():
    """A validation failure must not leave a readable output artifact.

    Use the inexpensive Sz-free fail-fast path and verify that
    ``zqp_orbital_uhfk.dat`` remains absent.
    """
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        np.savez(
            paths["occupation"],
            occupation=np.zeros((4, 2), dtype=np.float64),
            mu=np.array([0.0], dtype=np.float64),
            T=np.float64(0.0),
            column_spin=np.array([-1, -1], dtype=np.int64),
            column_mu_group=np.array([0, 0], dtype=np.int64),
        )
        assert not os.path.exists(paths["output"])
        result = subprocess.run(
            [
                sys.executable, "tools/uhfk_to_mvmc.py",
                "--input", paths["input"], "--eigen", paths["eigen"],
                "--occupation", paths["occupation"],
                "--geometry", paths["geometry"],
                "--orbitalidx", paths["orbitalidx"],
                "--output", paths["output"],
                "--no-check-density",
            ],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0
        assert not os.path.exists(paths["output"]), (
            "CLI left a zqp_orbital_uhfk.dat artifact even though "
            "validation rejected the inputs "
            f"(stderr: {result.stderr!r})"
        )


def test_general_path_rejects_ncond_odd():
    """CLI rejects odd Ncond with clear message when routed to General.

    ``derive_ne_per_group`` requires ``(Ncond +/- 2Sz)`` to be even, so
    the plan's "no 2Sz" wording is expressed here as ``2Sz = 1`` (odd)
    which yields integer ``Ne_up = 2``, ``Ne_down = 1`` and still routes
    to the General branch (``is_antiparallel_metadata`` needs ``2Sz==0``).
    """
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        # Overwrite input.toml with odd Ncond; 2Sz=1 keeps Ne_up/Ne_down
        # integer while disabling the antiparallel dispatch flag.
        with open(paths["input"], "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 3\n"
                "2Sz = 1\n"
                "T = 0.0\n"
                "CellShape = [4, 1, 1]\n"
                "SubShape  = [1, 1, 1]\n"
                'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
            )
        # Replace eigen.npz with non-degenerate eigenvalues so
        # step_occupation's Fermi-level degeneracy guard does not fire
        # before validate_general_prerequisites gets to reject odd Ncond.
        np.savez(
            paths["eigen"],
            eigenvalue=np.array(
                [[-2.0, -1.5], [-0.5, 0.0], [0.5, 1.0], [1.5, 2.0]],
                dtype=np.float64,
            ),
            eigenvector=np.zeros((4, 2, 2), dtype=np.complex128),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array(
                [[v, 0, 0] for v in [0, 1, -2, -1]], dtype=np.int64
            ),
            twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
        )
        # Overwrite orbitalidx.def with 6-column General format.
        with open(paths["orbitalidx"], "w") as fp:
            nsite = 4
            total = 2 * nsite * nsite - nsite
            fp.write("======================\n")
            fp.write(f"NOrbitalIdx  {total}\n")
            fp.write("ComplexType 1\n")
            fp.write("======================\n")
            fp.write("== i_spn_j_spn_OrbitalIdx ==\n")
            fp.write("======================\n")
            idx = 0
            for all_i in range(2 * nsite):
                for all_j in range(all_i + 1, 2 * nsite):
                    i, spn_i = all_i % nsite, all_i // nsite
                    j, spn_j = all_j % nsite, all_j // nsite
                    fp.write(f"{i} {spn_i} {j} {spn_j} {idx} 1\n")
                    idx += 1
            for k in range(total):
                fp.write(f"{k} 1\n")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"],
             "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"],
             "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0
        assert "odd" in result.stderr or "even" in result.stderr


def test_general_path_rejects_mixed_block_column_spin():
    """column_spin == -1 (mixed block, C case) is rejected in General."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = _write_minimal_inputs(tmp, subshape=(1, 1, 1))
        np.savez(
            paths["occupation"],
            occupation=np.zeros((4, 2), dtype=np.float64),
            mu=np.array([0.0], dtype=np.float64),
            T=np.float64(0.0),
            column_spin=np.array([-1, -1], dtype=np.int64),
            column_mu_group=np.array([0, 0], dtype=np.int64),
        )
        # Point orbitalidx.def to 6-column so dispatch enters General branch.
        with open(paths["orbitalidx"], "w") as fp:
            nsite = 4
            total = 2 * nsite * nsite - nsite
            fp.write("======================\n")
            fp.write(f"NOrbitalIdx  {total}\n")
            fp.write("ComplexType 1\n")
            fp.write("======================\n")
            fp.write("== i_spn_j_spn_OrbitalIdx ==\n")
            fp.write("======================\n")
            idx = 0
            for all_i in range(2 * nsite):
                for all_j in range(all_i + 1, 2 * nsite):
                    i, spn_i = all_i % nsite, all_i // nsite
                    j, spn_j = all_j % nsite, all_j // nsite
                    fp.write(f"{i} {spn_i} {j} {spn_j} {idx} 1\n")
                    idx += 1
            for k in range(total):
                fp.write(f"{k} 1\n")
        result = subprocess.run(
            [sys.executable, "tools/uhfk_to_mvmc.py",
             "--input", paths["input"], "--eigen", paths["eigen"],
             "--occupation", paths["occupation"],
             "--geometry", paths["geometry"],
             "--orbitalidx", paths["orbitalidx"],
             "--output", paths["output"],
             "--no-check-density"],
            capture_output=True, text=True, cwd=REPO_ROOT,
        )
        assert result.returncode != 0
        assert "mixed block" in result.stderr
