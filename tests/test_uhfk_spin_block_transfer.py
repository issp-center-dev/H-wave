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

fold.reshape_interaction had the worse variant of the same defect: its
normal-mode fold (a + norb_phys_orig * cellidx) relabeled spin-block indices
onto genuine sublattice orbitals, so the entries survived _make_ham_trans's
index guard and corrupted H(k) with spurious inter-orbital hoppings instead of
being dropped.  The fold now skips them too (warning once, since the solvers'
own spin-block warning can no longer fire for folded tables).
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

    def test_spin_block_entry_not_folded_into_hamiltonian(self):
        # spin-block hops at R=+-1 fold (with the physical stride) onto the
        # folded-orbital indices of neighboring cells, so they used to survive
        # _make_ham_trans's index guard and enter H(k) as spurious
        # inter-orbital hoppings; the correct H(k) keeps only the orbital
        # block: -1 * identity at every k
        d = tempfile.mkdtemp(prefix="uhfk_spinblock_fold_ham_")
        try:
            _write_soi_format_chain(
                d, ["   1    0    0    2    2  0.700000  0.0\n",
                    "  -1    0    0    2    2  0.700000  0.0\n"])
            solver = _build_solver(d, [6, 1, 1], subshape=[3, 1, 1])
            for orbvec in (ov for (_, ov) in solver.param_ham["Transfer"]):
                self.assertTrue(
                    all(0 <= o < solver.norb for o in orbvec),
                    "folded Transfer must not carry spin-block indices: "
                    "{}".format(orbvec))
            solver._make_ham_trans()
            np.testing.assert_allclose(
                solver.ham_trans,
                np.broadcast_to(-np.eye(solver.nd), solver.ham_trans.shape),
                atol=1e-12,
                err_msg="spin-block Transfer entries must not fold into H(k)")
        finally:
            shutil.rmtree(d, ignore_errors=True)


def _count_ignore_warnings(records):
    return sum(1 for r in records if "These terms are ignored" in r.getMessage())


class TestSpinBlockIgnoreWarningCount(unittest.TestCase):
    """The ignore-warning is emitted exactly once per run in every geometry.

    Without a fold it comes from _make_ham_trans; with one (including the
    default SubShape = CellShape) it comes from the fold, and _make_ham_trans
    then sees no spin-block entry, so it must not fire a second time.
    """

    def _count(self, cellshape, subshape):
        d = tempfile.mkdtemp(prefix="uhfk_spinblock_warn_")
        try:
            _write_soi_format_chain(
                d, ["   1    0    0    2    2  0.700000  0.0\n",
                    "  -1    0    0    2    2  0.700000  0.0\n"])
            with self.assertLogs("qlms", level="WARNING") as cm:
                solver = _build_solver(d, cellshape, subshape=subshape)
                solver._make_ham_trans()
            return _count_ignore_warnings(cm.records)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_no_fold_warns_once(self):
        # CellShape must exceed twice the hopping range (_check_cellsize)
        self.assertEqual(self._count([4, 1, 1], [1, 1, 1]), 1)

    def test_default_subshape_warns_once(self):
        self.assertEqual(self._count([3, 1, 1], None), 1)

    def test_proper_subshape_warns_once(self):
        self.assertEqual(self._count([6, 1, 1], [3, 1, 1]), 1)


class TestFoldSpinBlockDropIsOneBodyOnly(unittest.TestCase):
    """Only one-body tables (Transfer, Extern) may carry a droppable spin block.

    Two-body tables use physical orbital indices in every mode, so an index
    >= norb_phys_orig there is invalid input, never an ignorable spin term.
    The fold must fail closed on it instead of dropping it with the (wrong)
    advice to enable spin-orbital mode -- before this guard, a folded RPA run
    with such an entry completed with that warning, where the unfolded run
    (and the pre-fix code) crashed.
    """

    def _fold(self, ham, *, norb_so_orig, norb_phys_orig, drop_spin_block,
              enable_spin_orbital=False):
        from hwave.solver import fold
        return fold.reshape_interaction(
            ham, (2, 1, 1), (2, 1, 1),
            norb_so_orig=norb_so_orig, norb_phys_orig=norb_phys_orig,
            enable_spin_orbital=enable_spin_orbital,
            drop_spin_block=drop_spin_block)

    def test_two_body_out_of_range_index_raises_normal_mode(self):
        ham = {((1, 0, 0), (1, 1)): 0.3, ((-1, 0, 0), (1, 1)): 0.3}
        with self.assertRaisesRegex(ValueError, "orbital index"):
            self._fold(ham, norb_so_orig=1, norb_phys_orig=1,
                       drop_spin_block=False)

    def test_two_body_out_of_range_index_raises_spin_orbital_mode(self):
        # SO mode: geometry norb is 2 (spin-orbital), one physical orbital;
        # two-body tables are still physical-indexed, so index 1 is invalid
        ham = {((1, 0, 0), (1, 1)): 0.3, ((-1, 0, 0), (1, 1)): 0.3}
        with self.assertRaisesRegex(ValueError, "orbital index"):
            self._fold(ham, norb_so_orig=2, norb_phys_orig=1,
                       drop_spin_block=False)

    def test_two_body_highest_valid_index_is_retained(self):
        ham = {((1, 0, 0), (1, 1)): 0.3, ((-1, 0, 0), (1, 1)): 0.3}
        for norb_so_orig in (2, 4):   # normal mode / SO mode geometry
            out = self._fold(ham, norb_so_orig=norb_so_orig, norb_phys_orig=2,
                             drop_spin_block=False)
            self.assertTrue(
                np.isclose(sum(out.values()), 2 * sum(ham.values())),
                "a valid highest physical index must survive the fold "
                "(norb_so_orig={})".format(norb_so_orig))

    def test_one_body_spin_block_is_dropped_with_one_warning(self):
        ham = {((0, 0, 0), (0, 0)): -1.0,
               ((1, 0, 0), (1, 1)): 0.7, ((-1, 0, 0), (1, 1)): 0.7}
        with self.assertLogs("qlms", level="WARNING") as cm:
            out = self._fold(ham, norb_so_orig=1, norb_phys_orig=1,
                             drop_spin_block=True)
        self.assertEqual(_count_ignore_warnings(cm.records), 1)
        # the onsite term folds onto the two sublattice orbitals (0 and 1);
        # nothing from the spin block (index 1 in the ORIGINAL basis, which
        # the physical stride would have relabeled onto folded index 1 as
        # well) may leak in: total weight is exactly the folded onsite term
        self.assertTrue(all(o in (0, 1) for (_, ov) in out for o in ov))
        self.assertTrue(np.isclose(sum(out.values()), -2.0),
                        "spin-block hops leaked into the folded table: "
                        "{}".format(out))


def _write_rpa_chain(d, coul_orb):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n1\n0 0 0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("t\n1\n2\n 1 1\n 1 0 0 1 1 -1.0 0.0\n-1 0 0 1 1 -1.0 0.0\n")
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("V\n1\n2\n 1 1\n")
        f.write(" 1 0 0 {o} {o} 0.3 0.0\n-1 0 0 {o} {o} 0.3 0.0\n".format(
            o=coul_orb))


class TestRPAFoldRejectsOutOfRangeTwoBodyIndex(unittest.TestCase):
    """Public RPA path: a folded run must not silently drop an invalid
    CoulombInter orbital index (the unfolded run crashes on it)."""

    def _build(self, coul_orb, subshape):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_module
        d = tempfile.mkdtemp(prefix="rpa_spinblock_")
        try:
            _write_rpa_chain(d, coul_orb)
            info_input = {'path_to_input': d, 'interaction': {
                'path_to_input': d, 'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombInter': 'coulombinter.dat'}}
            info_mode = {'mode': 'RPA', 'calc_scheme': 'reduced', 'param': {
                'T': 0.5, 'Ncond': 1, 'CellShape': [4, 1, 1],
                'SubShape': subshape, 'Nmat': 16}}
            read_io = read_input_k.QLMSkInput(info_input)
            return rpa_module.RPA(read_io.get_param("ham"),
                                  {"print_level": 0, "print_step": 100},
                                  info_mode)
        finally:
            shutil.rmtree(d, ignore_errors=True)

    def test_out_of_range_index_under_fold_raises(self):
        with self.assertRaisesRegex(ValueError, "orbital index"):
            self._build(2, [2, 1, 1])

    def test_valid_index_under_fold_constructs(self):
        solver = self._build(1, [2, 1, 1])
        self.assertIn("CoulombInter", solver.ham_info.param_ham)


class TestRPAFoldDropsExternSpinBlock(unittest.TestCase):
    """Extern is one-body: a spin-block entry in a folded RPA run is dropped
    with the single ignore-warning, and the physical entry survives."""

    def test_extern_spin_block_dropped_physical_kept(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_module
        d = tempfile.mkdtemp(prefix="rpa_extern_spinblock_")
        try:
            _write_rpa_chain(d, 1)
            with open(os.path.join(d, "extern.dat"), "w") as f:
                f.write("h\n2\n1\n 1\n")
                f.write(" 0 0 0 1 1 0.25 0.0\n")     # physical onsite field
                f.write(" 0 0 0 2 2 0.90 0.0\n")     # spin block -> dropped
            info_input = {'path_to_input': d, 'interaction': {
                'path_to_input': d, 'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat', 'Extern': 'extern.dat'}}
            # explicit scheme: with no two-body term the chi0q-only path
            # requires one (calc_scheme is a top-level [mode] key)
            info_mode = {'mode': 'RPA', 'calc_scheme': 'reduced', 'param': {
                'T': 0.5, 'Ncond': 1, 'CellShape': [4, 1, 1],
                'SubShape': [2, 1, 1], 'Nmat': 16, 'coeff_extern': 1.0}}
            read_io = read_input_k.QLMSkInput(info_input)
            with self.assertLogs("qlms", level="WARNING") as cm:
                solver = rpa_module.RPA(
                    read_io.get_param("ham"),
                    {"print_level": 0, "print_step": 100}, info_mode)
            self.assertEqual(_count_ignore_warnings(cm.records), 1)
            ext = solver.ham_info.param_ham["Extern"]
            self.assertTrue(all(o in (0, 1) for (_, ov) in ext for o in ov),
                            "folded Extern must carry only sublattice-orbital "
                            "indices: {}".format(ext))
            # the physical onsite field appears once per sublattice cell
            self.assertTrue(np.isclose(sum(ext.values()), 2 * 0.25),
                            "spin-block Extern entry leaked into the fold: "
                            "{}".format(ext))
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
