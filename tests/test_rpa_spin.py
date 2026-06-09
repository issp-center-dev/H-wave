#!/usr/bin/env python3

import os, sys
import unittest
import copy

import numpy as np

class TestRPAspin(unittest.TestCase):
    def run_test(self, params):
        Lx, Ly, Nmat = 16,16,128

        spin_orbital = params.get("enable_spin_orbital", False)
        calc_scheme  = params.get("calc_scheme", "general")
        params_inter = params.get("inter", {})

        #----------------------------------------------------------------
        # input_dict = toml.load("tests/rpa/input.toml")
        info_log = {}
        info_mode = {
            'mode': 'RPA',
            'param': {
                'T': 2.0,
                #'mu': 0.0,
                'filling': 0.75,
                'CellShape': [Lx,Ly,1],
                'SubShape': [1,1,1],
                'Nmat': Nmat,
            },
            'enable_spin_orbital': spin_orbital,
            'calc_scheme': calc_scheme,
        }
        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input',
                'interaction': {
                    'path_to_input': 'tests/rpa/input',
                    'Geometry': 'geom_so.dat' if spin_orbital else 'geom.dat',
                    'CoulombIntra': 'coulombintra.dat',
                    'CoulombInter': 'coulombinter.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
                'chi0q': 'chi0q.npz',
            },
        }

        os.makedirs(info_file['output']['path_to_output'], exist_ok=True)

        # overwrite by params
        info_file['input']['interaction'].update(params_inter)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")

        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        solver.solve(green_info, info_file['output']['path_to_output'])

        # chi0q = green_info["chi0q"]
        chiq = copy.deepcopy(green_info["chiq"])

        solver.save_results(info_file['output'], green_info)

        #----------------------------------------------------------------
        # load chi0q and retry

        green_info.update( solver.read_init({
            'path_to_input': 'tests/rpa/output',
            'chi0q_init': 'chi0q.npz',
        }) )

        solver.solve(green_info, info_file['output']['path_to_output'])

        chiq_new = copy.deepcopy(green_info["chiq"])

        self.assertTrue(np.allclose(chiq, chiq_new), "check chiq")
        # self.assertTrue(True)
        pass

    #----------------------------------------------------------------
    def test_general_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer.dat',
            }
        })

    def test_general_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_general_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
            }
        })

    def test_general_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
            }
        })

    def test_general_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_general_spinful(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital.dat',
            }
        })

    #----------------------------------------------------------------
    def test_reduced_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer.dat',
            }
        })

    def test_reduced_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_reduced_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
            }
        })

    def test_reduced_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
            }
        })

    def test_reduced_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_reduced_spinful(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital.dat',
            }
        })

    #----------------------------------------------------------------
    def test_squashed_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer.dat',
            }
        })

    def test_squashed_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_squashed_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
            }
        })

    def test_squashed_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
            }
        })

    def test_squashed_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_squashed_spinful(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital.dat',
            }
        })

    # def test_squashed_spinful_with_interaction(self):
    #     self.run_test({
    #         'enable_spin_orbital': True,
    #         'calc_scheme': 'squashed',
    #         'inter': {
    #             'Transfer': 'transfer_spin_orbital.dat',
    #             'CoulombIntra': 'coulombintra.dat',
    #             'CoulombInter': 'coulombinter.dat',
    #         }
    #     })


class TestRPASublatticeInvarianceSO(unittest.TestCase):
    """Spin-orbital mode must be sublattice-invariant.

    Folding a physical lattice into a supercell (SubShape != [1,1,1]) describes
    the same physics, so the band structure and the chi0q spectrum must be
    identical with and without folding. This exercises the SO sublattice
    reshaping of both the transfer (SO indices, spin-decoded) and the
    interactions (physical-orbital indices) together -- a path the other RPA SO
    tests skip by always using SubShape=[1,1,1].
    """

    def _write_inputs(self, d):
        # 1 physical orbital (SO count = 2); SO transfer has on-site spin-flip + inter-cell hop
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("  1.0   0.0   0.0\n  0.0   1.0   0.0\n  0.0   0.0   1.0\n")
            f.write("2\n   0.0   0.0   0.0\n   0.0   0.0   0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer SO\n2\n6\n 1 1 1 1 1 1\n")
            f.write("   0    0    0    1    2  0.300000  0.0\n")  # on-site up<->down
            f.write("   0    0    0    2    1  0.300000  0.0\n")
            f.write("   1    0    0    1    1  1.000000  0.0\n")  # inter-cell, up
            f.write("  -1    0    0    1    1  1.000000  0.0\n")
            f.write("   1    0    0    2    2  1.000000  0.0\n")  # inter-cell, down
            f.write("  -1    0    0    2    2  1.000000  0.0\n")
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            # interaction uses the physical orbital index
            f.write("CoulombIntra\n1\n1\n 1\n   0    0    0    1    1   2.000000   0.0\n")

    def _run(self, d, subshape):
        info_mode = {
            'mode': 'RPA',
            'param': {'T': 1.0, 'mu': 0.0, 'CellShape': [2, 1, 1],
                      'SubShape': subshape, 'Nmat': 16},
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
        }
        info_input = {
            'path_to_input': d,
            'interaction': {'path_to_input': d, 'Geometry': 'geom.dat',
                            'Transfer': 'transfer.dat',
                            'CoulombIntra': 'coulombintra.dat'},
        }
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as solver_rpa
        read_io = read_input_k.QLMSkInput(info_input)
        solver = solver_rpa.RPA(read_io.get_param("ham"), {}, info_mode)
        green_info = read_io.get_param("green")
        solver.solve(green_info, os.path.join(d, "out"))

        bands = np.sort(np.ravel(solver.H0_eigenvalue.real))
        chi0 = green_info["chi0q"]
        nd = chi0.shape[-1]
        chi0_eig = np.sort(
            np.linalg.eigvals(chi0.reshape(-1, nd, nd)).real.ravel())
        return bands, chi0_eig

    def test_band_and_chi0q_sublattice_invariant(self):
        import tempfile
        import shutil
        d = tempfile.mkdtemp(prefix="rpa_so_sub_")
        try:
            self._write_inputs(d)
            b1, c1 = self._run(d, [1, 1, 1])
            b2, c2 = self._run(d, [2, 1, 1])
        finally:
            shutil.rmtree(d, ignore_errors=True)

        # same number of states / chi0q entries, just re-blocked
        self.assertEqual(len(b1), len(b2))
        self.assertEqual(len(c1), len(c2))
        np.testing.assert_allclose(
            b1, b2, atol=1e-8,
            err_msg="SO band structure must be sublattice-invariant")
        np.testing.assert_allclose(
            c1, c2, atol=1e-6,
            err_msg="SO chi0q spectrum must be sublattice-invariant")


if __name__ == '__main__':
    unittest.main()
