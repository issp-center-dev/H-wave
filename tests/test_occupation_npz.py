"""Tests for the UHFk._save_occupation() output file.

The new occupation.npz captures the converged SCF Fermi-Dirac weights,
the per-mu-group chemical potential, the SCF temperature, and the per
eigenvector-column metadata (spin character and mu-group index). The
file is consumed downstream by the H-wave -> mVMC PairProduct bridge
to project finite-T occupations onto Slater determinants without
re-running the SCF.

The tests drive UHFk through its normal SCF entry point so the saved state
matches what production runs emit. The bridge usage is documented in
``docs/en/source/uhfk/tools/uhfk_to_mvmc.rst``.
"""
from __future__ import annotations

import logging
import os
import tempfile

import numpy as np
import pytest

from hwave.qlmsio.read_input_k import QLMSkInput
from hwave.solver.uhfk import UHFk


# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

# Wannier90-like input files for a tiny 1D Hubbard chain. The shapes are
# chosen so the SCF converges in a handful of iterations and so each
# eigenvalue spectrum is non-degenerate at the Fermi level (no T=0
# tie-breaker ambiguity).
_GEOMETRY_DAT = (
    "  1.000000000000   0.000000000000   0.000000000000\n"
    "  0.000000000000   1.000000000000   0.000000000000\n"
    "  0.000000000000   0.000000000000   1.000000000000\n"
    "1\n"
    "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
)

# 1D nearest-neighbor hopping, t = 1 (Hermitian: +R and -R both present).
_TRANSFER_DAT = (
    "Transfer in wannier90-like format for uhfk\n"
    "1\n"
    "2\n"
    " 1 1\n"
    "  -1    0    0    1    1  -1.000000000000   0.000000000000\n"
    "   1    0    0    1    1  -1.000000000000   0.000000000000\n"
)

# On-site Hubbard U on the single orbital.
def _coulomb_intra_dat(U: float) -> str:
    return (
        "CoulombIntra in wannier90-like format for uhfk\n"
        "1\n"
        "1\n"
        " 1\n"
        "   0    0    0    1    1   {:.12f}   0.000000000000\n".format(U)
    )


def _write_inputs(tmp_path: str, U: float = 2.0) -> None:
    """Write the minimal Wannier90-like fixture files into tmp_path."""
    with open(os.path.join(tmp_path, "Geometry.dat"), "w") as f:
        f.write(_GEOMETRY_DAT)
    with open(os.path.join(tmp_path, "Transfer.dat"), "w") as f:
        f.write(_TRANSFER_DAT)
    with open(os.path.join(tmp_path, "CoulombIntra.dat"), "w") as f:
        f.write(_coulomb_intra_dat(U))


def _run_minimal_uhfk_for_occupation(
    tmp_path: str,
    *,
    T: float = 0.0,
    two_sz: int = 0,
    cell_l: int = 4,
    U: float = 2.0,
    ncond: int = 4,
) -> UHFk:
    """Run a small UHFk SCF and return the converged solver instance.

    The fixture is a 1D Hubbard chain of length ``cell_l`` with a single
    orbital per cell. Default is half filling (``ncond = cell_l``) with
    on-site repulsion ``U`` -- the CoulombIntra Hartree term couples
    up<->down on the same orbital, so the block detector produces a
    single mixed block (one mu-group). To exercise the Sz-fixed
    two-mu-group branch, call with ``U = 0`` (no Coulomb coupling ->
    two separate blocks) and pick a filling that leaves the Fermi level
    non-degenerate -- e.g. ``cell_l = 4`` with ``ncond = 2`` puts one
    electron per spin at the unique lowest eigenvalue ``E = -2t``.
    """
    _write_inputs(tmp_path, U=U)
    info_inputfile = {
        "path_to_input": "",
        "interaction": {
            "path_to_input": tmp_path,
            "Geometry": "Geometry.dat",
            "Transfer": "Transfer.dat",
            "CoulombIntra": "CoulombIntra.dat",
        },
    }
    read_io = QLMSkInput(info_inputfile)
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    info_log = {"print_level": 0, "print_step": 100}
    info_mode = {
        "mode": "UHFk",
        "param": {
            "Ncond": ncond,
            "2Sz": two_sz,
            "IterationMax": 200,
            "EPS": 8,
            "Mix": 0.5,
            "RndSeed": 12345,
            "T": T,
            "CellShape": [cell_l, 1, 1],
            "SubShape": [1, 1, 1],
        },
    }
    solver = UHFk(ham_info, info_log, info_mode)
    solver.solve(green_info, tmp_path)
    return solver


# ---------------------------------------------------------------------------
# Step 1: failing test (T=0, Sz-fixed two-block)
# ---------------------------------------------------------------------------

def test_save_occupation_writes_required_keys_t0():
    """occupation.npz must carry the 5 keys with the documented shapes/dtypes."""
    with tempfile.TemporaryDirectory() as tmp:
        solver = _run_minimal_uhfk_for_occupation(tmp, T=0.0, two_sz=0)
        occ_path = os.path.join(tmp, "occupation.npz")
        solver._save_occupation(occ_path)

        npz = np.load(occ_path, allow_pickle=False)
        assert set(npz.files) == {
            "occupation", "mu", "T", "column_spin", "column_mu_group"
        }
        # Shapes
        nd = solver.nd
        nvol = solver.nvol
        assert npz["occupation"].shape == (nvol, nd)
        assert npz["column_spin"].shape == (nd,)
        assert npz["column_mu_group"].shape == (nd,)
        assert npz["T"].shape == ()  # scalar
        # Dtypes
        assert npz["occupation"].dtype == np.float64
        assert npz["mu"].dtype == np.float64
        assert npz["T"].dtype == np.float64
        assert npz["column_spin"].dtype.kind == "i"
        assert npz["column_mu_group"].dtype.kind == "i"
        # T=0 occupations are 0 or 1 (no fractional weights)
        occ = npz["occupation"]
        unique_vals = np.unique(occ)
        assert set(unique_vals.tolist()).issubset({0.0, 1.0}), (
            "T=0 occupation should be 0/1, got values {}".format(unique_vals)
        )
        # Total electron count equals Ncond (here = cell_l = 4)
        np.testing.assert_allclose(occ.sum(), 4.0)
        # mu records the SCF chemical potential value
        np.testing.assert_allclose(npz["T"], 0.0)


def test_save_occupation_szfixed_has_two_mu_groups():
    """Sz-fixed UHF (2Sz=0) without spin-coupling interaction yields 2 mu-groups.

    Use ``U = 0`` so the block detector does NOT merge the up and down
    sectors via the CoulombIntra Hartree term (which is what couples
    them in the U>0 case). With two separate blocks, ``2Sz = 0`` gives
    one mu-group per spin sector.
    """
    with tempfile.TemporaryDirectory() as tmp:
        # Quarter filling (ncond=2 of total 4 sites x 2 spins = 8 states):
        # one electron per spin at k=0 (E=-2t). The single occupied state
        # per spin is unique and non-degenerate at the Fermi level, so
        # the T=0 step occupation is deterministic.
        solver = _run_minimal_uhfk_for_occupation(
            tmp, T=0.0, two_sz=0, U=0.0, ncond=2
        )
        occ_path = os.path.join(tmp, "occupation.npz")
        solver._save_occupation(occ_path)

        npz = np.load(occ_path, allow_pickle=False)
        # 2Sz=0 with no spin-coupling interaction keeps up and down in
        # separate blocks -> two mu-groups.
        assert npz["mu"].shape == (2,)
        # column_spin must mark each column as either up (0) or down (1).
        assert set(npz["column_spin"].tolist()) == {0, 1}
        # No "mixed" columns in Sz-fixed mode with decoupled blocks.
        assert -1 not in npz["column_spin"].tolist()
        # column_mu_group must reference both groups.
        assert set(npz["column_mu_group"].tolist()) == {0, 1}
        # Each column's mu-group must agree with its spin character: in
        # Sz-fixed mode block 0 is up (column_spin=0) -> mu-group 0,
        # block 1 is down (column_spin=1) -> mu-group 1.
        for n in range(solver.nd):
            spin = int(npz["column_spin"][n])
            grp = int(npz["column_mu_group"][n])
            assert spin == grp, (
                "spin/mu-group inconsistency at column {}: "
                "spin={}, mu_group={}".format(n, spin, grp)
            )
        # Each spin sector occupies ncond/2 = 1 state out of cell_l = 4.
        occ = npz["occupation"]
        n_up = occ[:, npz["column_spin"] == 0].sum()
        n_down = occ[:, npz["column_spin"] == 1].sum()
        np.testing.assert_allclose(n_up, 1.0)
        np.testing.assert_allclose(n_down, 1.0)


def test_save_occupation_t_positive_records_fractional_weights():
    """T>0 SCF leaves fractional Fermi-Dirac weights in the saved occupation."""
    # Pick T large enough relative to the bandwidth (4t for 1D NN hopping)
    # that some occupation falls in (0.01, 0.99). T=0.5 is comfortably
    # inside the smearing regime for t=1.
    with tempfile.TemporaryDirectory() as tmp:
        solver = _run_minimal_uhfk_for_occupation(
            tmp, T=0.5, two_sz=0
        )
        occ_path = os.path.join(tmp, "occupation.npz")
        solver._save_occupation(occ_path)

        npz = np.load(occ_path, allow_pickle=False)
        np.testing.assert_allclose(npz["T"], 0.5)
        occ = npz["occupation"]
        # At least one entry is fractional (not just 0/1).
        assert np.any((occ > 0.01) & (occ < 0.99)), (
            "T>0 SCF should produce some fractional occupations, got {}"
            .format(np.unique(occ))
        )
        # Total particle count still matches Ncond.
        np.testing.assert_allclose(occ.sum(), 4.0, atol=1.0e-6)


def test_save_occupation_called_from_save_results_when_keyword_present():
    """save_results() must write occupation.npz when the keyword is provided."""
    with tempfile.TemporaryDirectory() as tmp:
        solver = _run_minimal_uhfk_for_occupation(tmp, T=0.0, two_sz=0)
        info_outputfile = {
            "path_to_output": tmp,
            "occupation": "occupation.npz",
            "green": "",       # suppress green save (no green_info needed)
        }
        solver.save_results(info_outputfile, green_info={})
        assert os.path.exists(os.path.join(tmp, "occupation.npz"))


def test_save_occupation_not_called_when_keyword_absent():
    """Without the keyword, save_results() must NOT write occupation.npz.

    Backward compatibility: existing input.toml files that do not
    mention the new keyword keep their current output layout.
    """
    with tempfile.TemporaryDirectory() as tmp:
        solver = _run_minimal_uhfk_for_occupation(tmp, T=0.0, two_sz=0)
        info_outputfile = {
            "path_to_output": tmp,
            "green": "",  # suppress green save
        }
        solver.save_results(info_outputfile, green_info={})
        assert not os.path.exists(os.path.join(tmp, "occupation.npz"))
