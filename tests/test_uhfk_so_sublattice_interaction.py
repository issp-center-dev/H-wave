"""Regression test: UHFk spin-orbital + sublattice + two-body interaction.

In enable_spin_orbital mode the two-body interaction files use PHYSICAL
orbital indices [0, norb/2), and _make_ham_inter sizes its tables with
norb_phys = norb_phys_orig * subvol (layout a + norb_phys_orig * cellidx).
The sublattice fold of the interaction terms must therefore use the
physical-orbital stride norb_phys_orig, not the spin-orbital count
norb_orig: folding with the wrong stride crashes with IndexError (subvol=2)
or silently misplaces interaction terms.

The check runs the full SCF with an on-site U and asserts the converged
total energy is sublattice-invariant (SubShape [1,1,1] vs [2,1,1]), which
holds because the zero initial Green keeps the trajectory symmetric.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np


def _write_inputs(d, norb_phys=1):
    nso = 2 * norb_phys  # spin-orbital count (geometry norb in SO mode)
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("  1.0   0.0   0.0\n  0.0   1.0   0.0\n  0.0   0.0   1.0\n")
        f.write("{}\n".format(nso))
        for _ in range(nso):
            f.write("   0.0   0.0   0.0\n")
    rows = []
    # inter-cell spin-diagonal hops in x for all SO states
    for so in range(1, nso + 1):
        rows.append((1, 0, 0, so, so, 1.0))
        rows.append((-1, 0, 0, so, so, 1.0))
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer SO\n{}\n1\n 1\n".format(nso))
        for (rx, ry, rz, a, b, v) in rows:
            f.write("  {:3d} {:3d} {:3d} {:3d} {:3d}  {:.6f}  0.0\n".format(
                rx, ry, rz, a, b, v))
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n{}\n1\n 1\n".format(norb_phys))
        for p in range(1, norb_phys + 1):
            f.write("   0    0    0  {0:3d}  {0:3d}   2.000000   0.0\n".format(p))


def _run_energy(d, subshape, cellshape, norb_phys):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.uhfk as uhfk_module

    info_mode = {
        'mode': 'UHFk',
        'param': {'Ncond': 2 * norb_phys, 'T': 0.1, 'CellShape': cellshape,
                  'SubShape': subshape, 'IterationMax': 500, 'EPS': 10,
                  'Mix': 0.5, 'RndSeed': 1},
        'enable_spin_orbital': True,
    }
    info_input = {
        'path_to_input': d,
        'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                        'Transfer': 'transfer.dat',
                        'CoulombIntra': 'coulombintra.dat'},
    }
    read_io = read_input_k.QLMSkInput(info_input)
    info_log = {"print_level": 0, "print_step": 100}
    solver = uhfk_module.UHFk(read_io.get_param("ham"), info_log, info_mode)
    out = tempfile.mkdtemp(prefix="uhfk_so_out_")
    try:
        solver.solve(read_io.get_param("green"), out)
    finally:
        shutil.rmtree(out, ignore_errors=True)
    return solver.physics["Ene"]["Total"].real


class TestUHFkSOSublatticeInteraction(unittest.TestCase):
    def _check_invariant(self, norb_phys, cellshape, subshape):
        d = tempfile.mkdtemp(prefix="uhfk_so_int_")
        try:
            _write_inputs(d, norb_phys=norb_phys)
            e1 = _run_energy(d, [1, 1, 1], cellshape, norb_phys)
            e2 = _run_energy(d, subshape, cellshape, norb_phys)
        finally:
            shutil.rmtree(d, ignore_errors=True)
        self.assertTrue(
            np.isclose(e1, e2, rtol=0.0, atol=1e-8),
            "SO+interaction energy must be sublattice-invariant: "
            "E[1,1,1]={} vs E{}={}".format(e1, subshape, e2))

    def test_energy_invariant_1orb_subvol2(self):
        self._check_invariant(norb_phys=1, cellshape=[4, 1, 1],
                              subshape=[2, 1, 1])

    def test_energy_invariant_2orb_subvol2(self):
        self._check_invariant(norb_phys=2, cellshape=[4, 1, 1],
                              subshape=[2, 1, 1])


class TestUHFkSORejects2Sz(unittest.TestCase):
    def test_2sz_fixed_rejected_in_spin_orbital_mode(self):
        """Sz is not a good quantum number when spin is folded into the
        orbital index, so a fixed 2Sz cannot be enforced: block detection
        would silently assign the total Ncond to the spin-mixed blocks and
        ignore the constraint.  Fail fast instead."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.uhfk as uhfk_module

        d = tempfile.mkdtemp(prefix="uhfk_so_2sz_")
        try:
            _write_inputs(d, norb_phys=1)
            info_mode = {
                'mode': 'UHFk',
                'param': {'Ncond': 2, '2Sz': 0, 'T': 0.1,
                          'CellShape': [4, 1, 1], 'IterationMax': 100,
                          'EPS': 8, 'Mix': 0.5, 'RndSeed': 1},
                'enable_spin_orbital': True,
            }
            info_input = {
                'path_to_input': d,
                'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                                'Transfer': 'transfer.dat',
                                'CoulombIntra': 'coulombintra.dat'},
            }
            read_io = read_input_k.QLMSkInput(info_input)
            info_log = {"print_level": 0, "print_step": 100}
            with self.assertRaises(ValueError):
                uhfk_module.UHFk(read_io.get_param("ham"), info_log,
                                 info_mode)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
