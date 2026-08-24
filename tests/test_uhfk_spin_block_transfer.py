"""Regression tests: spin-block Transfer entries must not crash _check_hermite.

In normal mode (enable_spin_orbital = False) a Transfer table in Wannier90
spin-orbital format may carry orbital indices up to 2*norb; _check_orbital_index
deliberately accepts them and _make_ham_trans drops them with a warning.  But
_check_hermite used to scatter every Transfer entry into an array of shape
(nx,ny,nz,norb,norb), so any spin-block index >= norb raised an unhandled
IndexError during __init__ -- the documented "ignored with warning" path was
unreachable.  With a sublattice fold (where self.norb is inflated by subvol)
the entries fit into the array instead, but landed in sublattice-orbital slots
they do not belong to, so a dropped term could spuriously fail the strict
hermiticity check.

The fix excludes the spin-block entries from the hermiticity check, matching
_make_ham_trans dropping them from the Hamiltonian the check guards.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np


def _write_soi_format_chain(d, spin_block_lines):
    """One physical orbital; transfer file in spin-orbital format (num_wann=2)."""
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("  1.0   0.0   0.0\n  0.0   1.0   0.0\n  0.0   0.0   1.0\n")
        f.write("1\n   0.0   0.0   0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer SOI-format chain\n2\n1\n 1\n")
        f.write("   0    0    0    1    1  -1.000000  0.0\n")
        for line in spin_block_lines:
            f.write(line)


def _build_solver(d, cellshape, subshape=None):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.uhfk as uhfk_module

    param = {'Ncond': 1, 'T': 0.5, 'CellShape': cellshape,
             'IterationMax': 1, 'EPS': 8, 'Mix': 0.5, 'RndSeed': 1,
             'strict_hermite': True}
    if subshape is not None:
        param['SubShape'] = subshape
    info_mode = {'mode': 'UHFk', 'param': param}
    info_input = {
        'path_to_input': d,
        'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                        'Transfer': 'transfer.dat'},
    }
    read_io = read_input_k.QLMSkInput(info_input)
    info_log = {"print_level": 0, "print_step": 100}
    return uhfk_module.UHFk(read_io.get_param("ham"), info_log, info_mode)


class TestUHFkSpinBlockTransfer(unittest.TestCase):
    def test_spin_block_entry_does_not_crash_hermite_check(self):
        # index 2 (1-based) is the spin block of the single physical orbital;
        # the entry is deliberately non-hermitian (complex onsite) to pin down
        # that dropped entries are excluded from the check, not checked in a
        # larger array
        d = tempfile.mkdtemp(prefix="uhfk_spinblock_")
        try:
            _write_soi_format_chain(
                d, ["   0    0    0    2    2  0.500000  0.3\n"])
            try:
                solver = _build_solver(d, [2, 1, 1], subshape=[1, 1, 1])
            except (IndexError, SystemExit) as e:
                self.fail(
                    "UHFk init must ignore spin-block Transfer entries in the "
                    "hermiticity check, but raised {!r}".format(e))
            # the constructed Hamiltonian must contain only the orbital block,
            # duplicated over spin: T(k) = -1 * identity(2)
            solver._make_ham_trans()
            np.testing.assert_allclose(
                solver.ham_trans,
                np.broadcast_to(-np.eye(2), solver.ham_trans.shape),
                atol=1e-12,
                err_msg="spin-block Transfer entries must be dropped from H(k)")
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_spin_block_entry_not_misfolded_into_sublattice_orbitals(self):
        # with the default SubShape (= CellShape) self.norb is inflated by the
        # fold, so a spin-block index used to land inside the hermiticity-check
        # array -- in a sublattice-orbital slot it does not belong to.  An
        # unpaired hop at r=+1 misread as "orbital 1" is non-hermitian there
        # and used to abort the strict check spuriously.
        d = tempfile.mkdtemp(prefix="uhfk_spinblock_fold_")
        try:
            _write_soi_format_chain(
                d, ["   1    0    0    2    2  0.700000  0.0\n"])
            try:
                _build_solver(d, [3, 1, 1])
            except (IndexError, SystemExit) as e:
                self.fail(
                    "spin-block Transfer entries must not be folded into "
                    "sublattice orbital slots of the hermiticity check, but "
                    "UHFk init raised {!r}".format(e))
        finally:
            shutil.rmtree(d, ignore_errors=True)

if __name__ == "__main__":
    unittest.main()
