"""Tests for tools/uhfk_to_mvmc.py CLI format-first dispatch.

Verifies the four-case matrix:
  (is_antiparallel_metadata, orbitalidx_format) → routing / rejection

See ``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile

import numpy as np
import pytest

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))


def _minimal_general_orbitalidx(nsite=2):
    total = 2 * nsite * nsite - nsite
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
    for k in range(total):
        lines.append(f"{k} 1")
    return "\n".join(lines) + "\n"


def _minimal_antiparallel_orbitalidx(nsite=2):
    total = nsite * nsite
    lines = [
        "======================",
        f"NOrbitalIdx  {total}",
        "ComplexType 1",
        "======================",
        "== i_j_OrbitalIdx ==",
        "======================",
    ]
    idx = 0
    for i in range(nsite):
        for j in range(nsite):
            lines.append(f"{i} {j} {idx}")
            idx += 1
    for k in range(total):
        lines.append(f"{k} 1")
    return "\n".join(lines) + "\n"


def _write_fixture(tmp, sub_shape=(1, 1, 1), two_sz=0, column_spin=(0, 1)):
    """Emit input.toml + eigen.npz + occupation.npz + geometry_uhf.dat
    minimal enough for CLI to reach the dispatch branch. The mapping
    file path is filled in by the caller.

    Eigenvalues are set to distinct non-degenerate values so
    step_occupation's Fermi-level degeneracy guard does not trip.
    """
    Nsite = 2
    Ncond = 2
    with open(os.path.join(tmp, "input.toml"), "w") as fp:
        sub_str = ", ".join(str(int(x)) for x in sub_shape)
        fp.write(
            "[mode.param]\n"
            f"Ncond = {Ncond}\n"
            f"2Sz = {two_sz}\n"
            "T = 0.0\n"
            f"CellShape = [{Nsite}, 1, 1]\n"
            f"SubShape  = [{sub_str}]\n"
            'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
        )
    # Distinct eigenvalues per (nvol_row, column) so no Fermi-level
    # degeneracy fires. Column 0 (mu-group 0): 0.0 at row 0 < 0.5 at row 1.
    # Column 1 (mu-group 1): 1.0 at row 0 < 1.5 at row 1.
    eigenvalue = np.array(
        [[0.0, 1.0], [0.5, 1.5]], dtype=np.float64
    )
    np.savez(
        os.path.join(tmp, "eigen.npz"),
        eigenvalue=eigenvalue,
        eigenvector=np.eye(2, dtype=np.complex128).reshape(1, 2, 2).repeat(Nsite, axis=0),
        wavevector_unit=np.eye(3, dtype=np.float64),
        wavevector_index=np.array([[v, 0, 0] for v in [0, 1]], dtype=np.int64),
        twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
    )
    # Occupation: single up + single down at row 0 satisfies pair-closure for
    # PBC self-pair at k=0.
    occ = np.zeros((Nsite, 2), dtype=np.float64)
    occ[0, 0] = 1.0
    occ[0, 1] = 1.0
    np.savez(
        os.path.join(tmp, "occupation.npz"),
        occupation=occ,
        mu=np.array([0.0, 0.0], dtype=np.float64),
        T=np.float64(0.0),
        column_spin=np.array(column_spin, dtype=np.int64),
        column_mu_group=np.array([0, 1], dtype=np.int64),
    )
    with open(os.path.join(tmp, "geometry_uhf.dat"), "w") as fp:
        fp.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
        fp.write("0.0 0.0 0.0\n")
        fp.write("2 0 0\n0 1 0\n0 0 1\n")
        for i in range(Nsite):
            fp.write(f"{i} 0 0 0\n")


def _run_cli(tmp, orbitalidx_body, extra_args=()):
    orbidx = os.path.join(tmp, "orbitalidx.def")
    with open(orbidx, "w") as fp:
        fp.write(orbitalidx_body)
    out = os.path.join(tmp, "zqp_orbital_uhfk.dat")
    return subprocess.run(
        [sys.executable, "tools/uhfk_to_mvmc.py",
         "--input", os.path.join(tmp, "input.toml"),
         "--eigen", os.path.join(tmp, "eigen.npz"),
         "--occupation", os.path.join(tmp, "occupation.npz"),
         "--geometry", os.path.join(tmp, "geometry_uhf.dat"),
         "--orbitalidx", orbidx,
         "--output", out,
         "--no-check-density",
         "--epsilon-noise", "0",
         *extra_args],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )


def test_dispatch_antiparallel_metadata_plus_antiparallel_format_uses_v21_path():
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        res = _run_cli(tmp, _minimal_antiparallel_orbitalidx(nsite=2))
        assert res.returncode == 0, res.stderr


def test_dispatch_antiparallel_metadata_plus_general_format_forced_general():
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        # Should succeed via forced-General branch and possibly emit a
        # WARNING (pair closure holds, so no warning is expected here).
        assert res.returncode == 0, res.stderr


def test_dispatch_not_antiparallel_plus_antiparallel_format_rejected():
    """2Sz not provided (Sz-free) + antiparallel 3-column orbitalidx → reject."""
    with tempfile.TemporaryDirectory() as tmp:
        # Write input.toml without 2Sz (defaults treated as null in dispatch).
        _write_fixture(tmp, two_sz=0)
        # Overwrite input.toml to remove 2Sz.
        with open(os.path.join(tmp, "input.toml"), "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 2\n"
                "T = 0.0\n"
                "CellShape = [2, 1, 1]\n"
                "SubShape  = [1, 1, 1]\n"
                'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
            )
        res = _run_cli(tmp, _minimal_antiparallel_orbitalidx(nsite=2))
        assert res.returncode == 2, res.stderr
        assert "General" in res.stderr or "orbitalidx_general" in res.stderr


# ---------------------------------------------------------------------------
# Bridge boundary_input helper wiring
# ---------------------------------------------------------------------------


def _rewrite_input_toml(tmp, bc_entries=None, enable_soc=False, two_sz="0"):
    """Overwrite input.toml. ``bc_entries`` = None -> omit BoundaryCondition
    (absent-key default path). ``bc_entries`` = list -> emit the given list.
    """
    lines = [
        "[mode.param]",
        "Ncond = 2",
        f"2Sz = {two_sz}",
        "T = 0.0",
        "CellShape = [2, 1, 1]",
        "SubShape  = [1, 1, 1]",
    ]
    if bc_entries is not None:
        bc_str = ", ".join(f'"{b}"' for b in bc_entries)
        lines.append(f"BoundaryCondition = [{bc_str}]")
    if enable_soc:
        lines.append("enable_spin_orbital = true")
    with open(os.path.join(tmp, "input.toml"), "w") as fp:
        fp.write("\n".join(lines) + "\n")


def _rewrite_eigen_npz(tmp, twist_offset, *, apbc_shifted=False):
    """Overwrite eigen.npz with configurable twist_offset. When ``apbc_shifted``
    is True, place eigenvalue minima so pair-closure succeeds under APBC(x).
    """
    Nsite = 2
    if apbc_shifted:
        eigenvalue = np.array([[0.0, 1.5], [0.5, 1.0]], dtype=np.float64)
    else:
        eigenvalue = np.array([[0.0, 1.0], [0.5, 1.5]], dtype=np.float64)
    kwargs = dict(
        eigenvalue=eigenvalue,
        eigenvector=np.eye(2, dtype=np.complex128).reshape(1, 2, 2).repeat(Nsite, axis=0),
        wavevector_unit=np.eye(3, dtype=np.float64),
        wavevector_index=np.array([[v, 0, 0] for v in [0, 1]], dtype=np.int64),
    )
    if twist_offset is not None:
        kwargs["twist_offset"] = np.array(twist_offset, dtype=np.float64)
    np.savez(os.path.join(tmp, "eigen.npz"), **kwargs)


@pytest.mark.parametrize(
    "alias", ["AP", "ap", "  antiperiodic  ", "ANTIPERIODIC"]
)
def test_dispatch_soc_apbc_alias_accepted(alias):
    """enable_spin_orbital = true + APBC alias -> reaches General branch.

    The boundary alias helper canonicalizes every AP spelling to theta =
    pi, and dispatch routes to the General-SOC branch. A downstream failure
    (e.g. a tiny fixture without a valid
    Rashba Transfer.dat) is acceptable here because the point is the
    dispatch gate. See ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    End-to-end validation lives in
    ``tests/validation/uhfk_mvmc_pairproduct/case_soc_rashba_2d_nosub_apbc``.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp)
        _rewrite_input_toml(
            tmp, bc_entries=[alias, "periodic", "periodic"], enable_soc=True,
        )
        _rewrite_eigen_npz(
            tmp, twist_offset=[0.5, 0.0, 0.0], apbc_shifted=True,
        )
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        # The reject must NOT be a boundary parsing error.
        assert "ERROR (boundary)" not in res.stderr, res.stderr


@pytest.mark.parametrize(
    "alias", ["ap", "AP", "  antiperiodic  ", "ANTIPERIODIC"]
)
def test_dispatch_nonsoc_apbc_alias_uniform(alias):
    """Non-SOC + APBC alias -> byte-identical output to canonical 'ap'."""
    canonical_bytes = None
    canonical_stderr = None
    canonical_rc = None
    with tempfile.TemporaryDirectory() as tmp0:
        _write_fixture(tmp0)
        _rewrite_input_toml(
            tmp0, bc_entries=["ap", "periodic", "periodic"],
        )
        _rewrite_eigen_npz(
            tmp0, twist_offset=[0.5, 0.0, 0.0], apbc_shifted=True,
        )
        res0 = _run_cli(tmp0, _minimal_general_orbitalidx(nsite=2))
        canonical_rc = res0.returncode
        canonical_stderr = res0.stderr
        out_path = os.path.join(tmp0, "zqp_orbital_uhfk.dat")
        if os.path.exists(out_path):
            with open(out_path, "rb") as fp:
                canonical_bytes = fp.read()

    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp)
        _rewrite_input_toml(
            tmp, bc_entries=[alias, "periodic", "periodic"],
        )
        _rewrite_eigen_npz(
            tmp, twist_offset=[0.5, 0.0, 0.0], apbc_shifted=True,
        )
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        assert res.returncode == canonical_rc, (res.returncode, res.stderr)
        out_path = os.path.join(tmp, "zqp_orbital_uhfk.dat")
        if canonical_bytes is not None:
            with open(out_path, "rb") as fp:
                assert fp.read() == canonical_bytes


@pytest.mark.parametrize(
    "ambiguous",
    ["open", "aperiodic", "", "foo", "peri odic", "antiperiodicX"],
)
def test_dispatch_rejects_ambiguous_boundary_pre_dispatch(ambiguous):
    """Invalid BoundaryCondition string -> CLI fails before dispatch."""
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp)
        _rewrite_input_toml(
            tmp, bc_entries=[ambiguous, "periodic", "periodic"],
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        assert res.returncode != 0, res.stderr
        assert "ERROR (boundary)" in res.stderr, res.stderr
        assert not os.path.exists(
            os.path.join(tmp, "zqp_orbital_uhfk.dat")
        ), "output file leaked despite pre-dispatch reject"


def test_dispatch_absent_boundary_defaults_to_pbc():
    """Omitted BoundaryCondition -> byte-identical to explicit all-PBC."""
    with tempfile.TemporaryDirectory() as tmp_ref:
        _write_fixture(tmp_ref)
        _rewrite_input_toml(
            tmp_ref, bc_entries=["periodic", "periodic", "periodic"],
        )
        res_ref = _run_cli(tmp_ref, _minimal_general_orbitalidx(nsite=2))
        assert res_ref.returncode == 0, res_ref.stderr
        with open(
            os.path.join(tmp_ref, "zqp_orbital_uhfk.dat"), "rb"
        ) as fp:
            ref_bytes = fp.read()

    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp)
        _rewrite_input_toml(tmp, bc_entries=None)  # omit key entirely
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        assert res.returncode == 0, res.stderr
        with open(
            os.path.join(tmp, "zqp_orbital_uhfk.dat"), "rb"
        ) as fp:
            assert fp.read() == ref_bytes


def test_dispatch_rejects_eigen_twist_mismatch():
    """input.toml APBC vs eigen.npz PBC twist_offset -> pre-dispatch reject."""
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp)  # writes twist_offset = [0, 0, 0]
        _rewrite_input_toml(
            tmp, bc_entries=["antiperiodic", "periodic", "periodic"],
        )
        # keep twist_offset from _write_fixture ([0,0,0]) → mismatch
        res = _run_cli(tmp, _minimal_general_orbitalidx(nsite=2))
        assert res.returncode != 0, res.stderr
        assert "ERROR (eigen twist)" in res.stderr, res.stderr
        assert not os.path.exists(
            os.path.join(tmp, "zqp_orbital_uhfk.dat")
        ), "output file leaked despite twist mismatch"


def test_normalize_boundary_condition_list_shape():
    """Direct unit-level shape checks for the helper (belt-and-braces)."""
    from tools._uhfk_to_mvmc.boundary_input import (
        BoundaryInputError, normalize_boundary_condition_list,
    )
    with pytest.raises(BoundaryInputError):
        normalize_boundary_condition_list("periodic")
    with pytest.raises(BoundaryInputError):
        normalize_boundary_condition_list(["periodic", None, "periodic"])
    with pytest.raises(BoundaryInputError):
        normalize_boundary_condition_list([1, 2, 3])


# ---------------------------------------------------------------------------
# is_soc_mode dispatch flag and six-case matrix
# ---------------------------------------------------------------------------


def test_dispatch_soc_mode_forces_general_format():
    """enable_spin_orbital = true + 3-col antiparallel orbitalidx -> reject.

    The 3/4-column antiparallel format does not
    carry spin-off-diagonal classes, so SOC dispatch requires the
    6-column General format. The CLI must reject before entering any
    downstream builder and emit an error mentioning ``6-column``.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        _rewrite_input_toml(
            tmp, bc_entries=["periodic", "periodic", "periodic"],
            enable_soc=True,
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        res = _run_cli(tmp, _minimal_antiparallel_orbitalidx(nsite=2))
        assert res.returncode != 0, res.stderr
        assert "6-column" in res.stderr, res.stderr
        assert not os.path.exists(
            os.path.join(tmp, "zqp_orbital_uhfk.dat")
        ), "output file leaked despite reject"


def test_dispatch_soc_mode_general_ok():
    """enable_spin_orbital = true + 6-col orbitalidx -> reaches General-SOC.

    The ``(*, general, True)`` cell of the dispatch matrix
    routes to the General branch with ``is_soc_mode=True``. Tasks 5-10
    add SOC-specific behavior inside the receiving modules. The dispatch
    reaching the General branch is verified by (a) no SOC+antiparallel
    ``6-column`` reject and (b) a successful run with the ``wrote ...
    General params`` stdout line. The SOC path additionally requires
    ``--transfer`` and ``--emit-trans`` under SOC, so the test provides a
    minimal Transfer.dat and asks the bridge to emit trans.def alongside
    ``zqp_orbital_uhfk.dat``.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        _rewrite_input_toml(
            tmp, bc_entries=["periodic", "periodic", "periodic"],
            enable_soc=True,
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        # The SOC path emits trans.def; provide a minimal
        # Wannier90-like Transfer.dat with two spin-diagonal entries
        # (R = (+1, 0, 0) and R = (-1, 0, 0)) so the emitter has data
        # to walk without inventing physics beyond the dispatch test's
        # scope.
        transfer_path = os.path.join(tmp, "Transfer.dat")
        with open(transfer_path, "w") as fp:
            fp.write("test Transfer.dat\n")
            fp.write("1\n")
            fp.write("2\n")
            fp.write("1 1\n")
            fp.write("   1  0  0  1  1  -1.0  0.0\n")
            fp.write("  -1  0  0  1  1  -1.0  0.0\n")
        emit_trans_path = os.path.join(tmp, "trans.def.bridge")
        res = _run_cli(
            tmp, _minimal_general_orbitalidx(nsite=2),
            extra_args=(
                "--transfer", transfer_path,
                "--emit-trans", emit_trans_path,
            ),
        )
        assert res.returncode == 0, res.stderr
        # No SOC+antiparallel reject fired.
        assert "6-column" not in res.stderr, res.stderr
        # Dispatch actually reached the General branch (writes params).
        assert "General params" in res.stdout, res.stdout
        assert os.path.exists(
            os.path.join(tmp, "zqp_orbital_uhfk.dat")
        )
        # SOC dispatch also emits trans.def alongside.
        assert os.path.exists(emit_trans_path), (
            "SOC bridge should emit trans.def when --emit-trans is set"
        )


def test_dispatch_soc_subshape_no_longer_rejects():
    """enable_spin_orbital = true + SubShape > [1, 1, 1] is supported.

    Downstream errors (the minimal fixture is not a
    converged SCF, so the density check or Slater emit may still error)
    are permissible; only the SubShape dispatch guard is under test. See
    ``docs/en/source/algorithm/uhfk_to_mvmc.rst``.
    """
    Nsite = 4
    subvol = 2
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "input.toml"), "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 2\n"
                "2Sz = 0\n"
                "T = 0.0\n"
                f"CellShape = [{Nsite}, 1, 1]\n"
                "SubShape  = [2, 1, 1]\n"
                'BoundaryCondition = ["periodic", "periodic", "periodic"]\n'
                "enable_spin_orbital = true\n"
            )
        nd = 2 * subvol  # spinful, subvol columns per spin.
        L_folded = Nsite // subvol
        np.savez(
            os.path.join(tmp, "eigen.npz"),
            eigenvalue=np.zeros((L_folded, nd), dtype=np.float64),
            eigenvector=np.eye(nd, dtype=np.complex128).reshape(1, nd, nd)
            .repeat(L_folded, axis=0),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array(
                [[v, 0, 0] for v in range(L_folded)], dtype=np.int64,
            ),
            twist_offset=np.array([0.0, 0.0, 0.0], dtype=np.float64),
        )
        occ = np.zeros((Nsite, 2), dtype=np.float64)
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
        transfer_path = os.path.join(tmp, "Transfer.dat")
        with open(transfer_path, "w") as fp:
            fp.write("test Transfer.dat\n")
            fp.write("1\n")
            fp.write("2\n")
            fp.write("1 1\n")
            fp.write("   1  0  0  1  1  -1.0  0.0\n")
            fp.write("  -1  0  0  1  1  -1.0  0.0\n")
        emit_trans_path = os.path.join(tmp, "trans.def.bridge")
        emit_orbitalidx_path = os.path.join(tmp, "orbitalidxgen.def.bridge")
        res = _run_cli(
            tmp, _minimal_general_orbitalidx(nsite=Nsite),
            extra_args=(
                "--transfer", transfer_path,
                "--emit-trans", emit_trans_path,
                "--emit-orbitalidx", emit_orbitalidx_path,
            ),
        )
        # The dispatch guards must not fire: with --emit-orbitalidx supplied
        # and no antiperiodic direction, SOC + SubShape > [1, 1, 1] is
        # supported, so any failure has to come from a later stage.
        assert "requires --emit-orbitalidx" not in res.stderr, res.stderr
        assert (
            "unsupported SOC + APBC + SubShape combination" not in res.stderr
        ), res.stderr


def test_dispatch_soc_single_apbc_arbitrary_subshape_rejected_by_v37_allowlist():
    """The allowlist contains two exact shipping shapes:
      (a) CellShape=[6,4,1] / SubShape=[2,2,1]
      (b) CellShape=[4,4,4] / SubShape=[2,2,2]

    An arbitrary non-shipping shape like CellShape=[4,1,1] /
    SubShape=[2,1,1] with single-direction APBC MUST be rejected
    pre-dispatch with the allowlist's REJECT_MESSAGE.
    """
    Nsite = 4
    subvol = 2
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "input.toml"), "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 2\n"
                "2Sz = 0\n"
                "T = 0.0\n"
                f"CellShape = [{Nsite}, 1, 1]\n"
                "SubShape  = [2, 1, 1]\n"
                'BoundaryCondition = ["antiperiodic", "periodic", "periodic"]\n'
                "enable_spin_orbital = true\n"
            )
        nd = 2 * subvol  # spinful, subvol columns per spin.
        L_folded = Nsite // subvol
        np.savez(
            os.path.join(tmp, "eigen.npz"),
            eigenvalue=np.zeros((L_folded, nd), dtype=np.float64),
            eigenvector=np.eye(nd, dtype=np.complex128).reshape(1, nd, nd)
            .repeat(L_folded, axis=0),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array(
                [[v, 0, 0] for v in range(L_folded)], dtype=np.int64,
            ),
            # APBC in x -> twist_offset[0] = 0.5.
            twist_offset=np.array([0.5, 0.0, 0.0], dtype=np.float64),
        )
        occ = np.zeros((Nsite, 2), dtype=np.float64)
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
        transfer_path = os.path.join(tmp, "Transfer.dat")
        with open(transfer_path, "w") as fp:
            fp.write("test Transfer.dat\n")
            fp.write("1\n")
            fp.write("2\n")
            fp.write("1 1\n")
            fp.write("   1  0  0  1  1  -1.0  0.0\n")
            fp.write("  -1  0  0  1  1  -1.0  0.0\n")
        emit_trans_path = os.path.join(tmp, "trans.def.bridge")
        emit_orbitalidx_path = os.path.join(tmp, "orbitalidxgen.def.bridge")
        res = _run_cli(
            tmp, _minimal_general_orbitalidx(nsite=Nsite),
            extra_args=(
                "--transfer", transfer_path,
                "--emit-trans", emit_trans_path,
                "--emit-orbitalidx", emit_orbitalidx_path,
            ),
        )
        from tools._uhfk_to_mvmc.allowlist_predicate import REJECT_MESSAGE

        # Single-direction APBC on a non-shipping CellShape/SubShape
        # (here [4,1,1]/[2,1,1]) is not one of the two allowlisted
        # (apbc_mask, sub_shape, cell_shape) triples, so it must be
        # rejected pre-dispatch with REJECT_MESSAGE.
        assert res.returncode == 2, (
            "Supported allowlist: single-direction APBC on a non-shipping "
            "CellShape/SubShape must be rejected pre-dispatch. Got "
            f"returncode={res.returncode}; stderr={res.stderr!r}"
        )
        assert REJECT_MESSAGE in res.stderr, res.stderr


# ---------------------------------------------------------------------------
# SOC output preflight
# ---------------------------------------------------------------------------


def _write_soc_transfer_dat(path):
    """Minimal Wannier90-like Transfer.dat with two spin-diagonal entries.

    Matches the fixture used by ``test_dispatch_soc_mode_general_ok`` so
    the SOC preflight tests don't couple to the trans_emit parser's
    Wannier90 details beyond header + one ndegen line + two data rows.
    """
    with open(path, "w") as fp:
        fp.write("test Transfer.dat\n")
        fp.write("1\n")
        fp.write("2\n")
        fp.write("1 1\n")
        fp.write("   1  0  0  1  1  -1.0  0.0\n")
        fp.write("  -1  0  0  1  1  -1.0  0.0\n")


def _run_cli_soc_with_paths(
    tmp, orbitalidx_body, *, output_path, emit_trans_path,
):
    """Invoke the CLI with explicit --output and --emit-trans (SOC mode).

    Used by the preflight regression tests because ``_run_cli`` hard-codes
    ``--output`` to ``tmp/zqp_orbital_uhfk.dat`` and the same-path /
    directory-target tests need to point both flags at custom locations.
    """
    orbidx = os.path.join(tmp, "orbitalidx.def")
    with open(orbidx, "w") as fp:
        fp.write(orbitalidx_body)
    transfer_path = os.path.join(tmp, "Transfer.dat")
    _write_soc_transfer_dat(transfer_path)
    return subprocess.run(
        [
            sys.executable, "tools/uhfk_to_mvmc.py",
            "--input", os.path.join(tmp, "input.toml"),
            "--eigen", os.path.join(tmp, "eigen.npz"),
            "--occupation", os.path.join(tmp, "occupation.npz"),
            "--geometry", os.path.join(tmp, "geometry_uhf.dat"),
            "--orbitalidx", orbidx,
            "--output", output_path,
            "--no-check-density",
            "--epsilon-noise", "0",
            "--transfer", transfer_path,
            "--emit-trans", emit_trans_path,
        ],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )


def test_dispatch_rejects_soc_output_same_path():
    """SOC with --output == --emit-trans -> pre-write reject.

    Two independent ``os.replace`` calls are not
    truly atomic. If both flags target the same file, the second call
    would silently overwrite the first commit; the preflight must catch
    this before either temp file is opened.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        _rewrite_input_toml(
            tmp, bc_entries=["periodic", "periodic", "periodic"],
            enable_soc=True,
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        shared_path = os.path.join(tmp, "collision.dat")
        res = _run_cli_soc_with_paths(
            tmp, _minimal_general_orbitalidx(nsite=2),
            output_path=shared_path, emit_trans_path=shared_path,
        )
        assert res.returncode != 0, res.stderr
        assert "resolve to the same path" in res.stderr, res.stderr
        # The reject fires before any write; the shared path must not
        # have been created and no stray .tmp file may leak.
        assert not os.path.exists(shared_path), (
            "output file leaked despite same-path reject"
        )
        stray = [
            n for n in os.listdir(tmp)
            if n.startswith(".uhfk_to_mvmc.") and n.endswith(".tmp")
        ]
        assert stray == [], f"stray tmp files not cleaned up: {stray}"


def test_dispatch_rejects_soc_output_directory_target():
    """SOC with --output pointing at an existing directory -> reject.

    ``os.replace`` cannot atomically overwrite a
    directory with a file, so the preflight rejects this up front.
    Confirms the reject fires without a partial write and that the
    directory is untouched.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        _rewrite_input_toml(
            tmp, bc_entries=["periodic", "periodic", "periodic"],
            enable_soc=True,
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        dir_target = os.path.join(tmp, "existing_dir")
        os.makedirs(dir_target)
        emit_trans_path = os.path.join(tmp, "trans.def.bridge")
        res = _run_cli_soc_with_paths(
            tmp, _minimal_general_orbitalidx(nsite=2),
            output_path=dir_target, emit_trans_path=emit_trans_path,
        )
        assert res.returncode != 0, res.stderr
        assert "existing directory" in res.stderr, res.stderr
        # The directory must survive intact.
        assert os.path.isdir(dir_target), "target directory clobbered"
        # No emit-trans output leaked either.
        assert not os.path.exists(emit_trans_path), (
            "trans.def leaked despite directory-target reject"
        )


def test_dispatch_rejects_soc_emit_trans_directory_target():
    """SOC with --emit-trans pointing at an existing directory -> reject.

    Complement of ``test_dispatch_rejects_soc_output_directory_target``:
    the preflight runs the same directory check on both --output and
    --emit-trans, so pointing --emit-trans at a directory must also
    fail-fast.
    """
    with tempfile.TemporaryDirectory() as tmp:
        _write_fixture(tmp, two_sz=0, column_spin=(0, 1))
        _rewrite_input_toml(
            tmp, bc_entries=["periodic", "periodic", "periodic"],
            enable_soc=True,
        )
        _rewrite_eigen_npz(tmp, twist_offset=[0.0, 0.0, 0.0])
        dir_target = os.path.join(tmp, "existing_dir")
        os.makedirs(dir_target)
        output_path = os.path.join(tmp, "zqp_orbital_uhfk.dat")
        res = _run_cli_soc_with_paths(
            tmp, _minimal_general_orbitalidx(nsite=2),
            output_path=output_path, emit_trans_path=dir_target,
        )
        assert res.returncode != 0, res.stderr
        assert "existing directory" in res.stderr, res.stderr
        assert os.path.isdir(dir_target), "target directory clobbered"
        assert not os.path.exists(output_path), (
            "zqp output leaked despite emit-trans directory reject"
        )


def test_dispatch_soc_multi_apbc_subshape_rejects_pre_dispatch():
    """SOC + multi-direction APBC
    (xy mask) + SubShape=[2,2,1] on CellShape=[4,4,1] is not one of the
    two validated (apbc_mask, sub_shape, cell_shape) triples in
    ``tools._uhfk_to_mvmc.allowlist_predicate`` (the shipping shapes are
    CellShape=[6,4,1]/SubShape=[2,2,1] and
    CellShape=[4,4,4]/SubShape=[2,2,2]), so it remains pre-dispatch
    rejected. Complements the shape-narrowing pin
    ``test_dispatch_soc_single_apbc_arbitrary_subshape_rejected_by_v37_allowlist``.
    """
    Nsite = 4
    subvol = 2
    with tempfile.TemporaryDirectory() as tmp:
        with open(os.path.join(tmp, "input.toml"), "w") as fp:
            fp.write(
                "[mode.param]\n"
                "Ncond = 2\n"
                "2Sz = 0\n"
                "T = 0.0\n"
                f"CellShape = [{Nsite}, {Nsite}, 1]\n"
                "SubShape  = [2, 2, 1]\n"
                'BoundaryCondition = ["antiperiodic", "antiperiodic", "periodic"]\n'
                "enable_spin_orbital = true\n"
            )
        nd = 2 * subvol
        L_folded = Nsite // subvol
        np.savez(
            os.path.join(tmp, "eigen.npz"),
            eigenvalue=np.zeros((L_folded, nd), dtype=np.float64),
            eigenvector=np.eye(nd, dtype=np.complex128).reshape(1, nd, nd)
            .repeat(L_folded, axis=0),
            wavevector_unit=np.eye(3, dtype=np.float64),
            wavevector_index=np.array(
                [[v, 0, 0] for v in range(L_folded)], dtype=np.int64,
            ),
            twist_offset=np.array([0.5, 0.5, 0.0], dtype=np.float64),
        )
        occ = np.zeros((Nsite * Nsite, 2), dtype=np.float64)
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
            fp.write(f"{Nsite} 0 0\n0 {Nsite} 0\n0 0 1\n")
            for iy in range(Nsite):
                for ix in range(Nsite):
                    fp.write(f"{ix} {iy} 0 0\n")
        transfer_path = os.path.join(tmp, "Transfer.dat")
        with open(transfer_path, "w") as fp:
            fp.write("test Transfer.dat\n")
            fp.write("1\n")
            fp.write("2\n")
            fp.write("1 1\n")
            fp.write("   1  0  0  1  1  -1.0  0.0\n")
            fp.write("  -1  0  0  1  1  -1.0  0.0\n")
        emit_trans_path = os.path.join(tmp, "trans.def.bridge")
        emit_orbitalidx_path = os.path.join(tmp, "orbitalidxgen.def.bridge")
        res = _run_cli(
            tmp, _minimal_general_orbitalidx(nsite=Nsite * Nsite),
            extra_args=(
                "--transfer", transfer_path,
                "--emit-trans", emit_trans_path,
                "--emit-orbitalidx", emit_orbitalidx_path,
            ),
        )
        from tools._uhfk_to_mvmc.allowlist_predicate import REJECT_MESSAGE

        assert res.returncode == 2, (
            "Supported allowlist: SOC + APBC + SubShape combination not "
            "covered by tools._uhfk_to_mvmc.allowlist_predicate must "
            f"fire pre-dispatch. Got returncode={res.returncode}; "
            f"stderr={res.stderr!r}"
        )
        assert REJECT_MESSAGE in res.stderr, res.stderr
