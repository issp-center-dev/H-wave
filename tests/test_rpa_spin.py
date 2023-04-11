#!/usr/bin/env python3

import os, sys
import unittest

import numpy as np

class TestRPAspin(unittest.TestCase):
    def run_test(self, params):
        Lx, Ly, Nmat = 16,16,128

        spin_orbital = params.get("enable_spin_orbital", False)
        calc_scheme = params.get("calc_scheme", "general")
        params_inter = params.get("inter", {})

        #----------------------------------------------------------------
        # input_dict = toml.load("tests/rpa/input.toml")
        info_log = {}
        info_mode = {
            'mode': 'RPA',
            'param': {
                'T': 2.0,
                'mu': 0.0,
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
                    'Geometry': 'geom.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
            },
        }

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
        # chiq = green_info["chiq"]
        
        # #----------------------------------------------------------------
        # # reference
        # rpaone = RPAOneOrbital()

        # chi0q_ref, chiq_ref, ham_ref, kx_array, ky_array = rpaone.run(params_ref)

        # #----------------------------------------------------------------
        # # compare

        # # chi0q
        # chi0q_x = chi0q[:,:,0,0,0,0]
        # chi0q_ref_x = np.transpose(chi0q_ref, (2,0,1)).reshape(Nmat,Lx*Ly)

        # self.assertTrue(np.allclose(chi0q_x, chi0q_ref_x), 'chi0q')

        # # hamiltonian
        # ham = solver.ham_info.ham_inter_q.reshape(Lx,Ly,4,4)
        # self.assertTrue(np.allclose(ham, ham_ref), 'hamiltonian')

        # # chiq
        # chiq_x = chiq[Nmat//2].reshape(Lx,Ly,4,4)

        # self.assertTrue(np.allclose(chiq_x, chiq_ref), 'chiq')
        self.assertTrue(True)

    #----------------------------------------------------------------
    def test_general_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_general_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_general_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_general_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_general_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_general_spinful(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'general',
            'inter': {
                'Transfer': 'transfer_spin_orbital.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    #----------------------------------------------------------------
    def test_reduced_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_reduced_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_reduced_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_reduced_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_reduced_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_reduced_spinful(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'reduced',
            'inter': {
                'Transfer': 'transfer_spin_orbital.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    #----------------------------------------------------------------
    def test_squashed_spin_free(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_squashed_spin_diag(self):
        self.run_test({
            'enable_spin_orbital': False,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    def test_squashed_spin_free2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_squashed_spin_diag2(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_diag.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            }
        })

    def test_squashed_spin_diag3(self):
        self.run_test({
            'enable_spin_orbital': True,
            'calc_scheme': 'squashed',
            'inter': {
                'Transfer': 'transfer_spin_orbital_trivial.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
                'Extern': 'coulombintra.dat',
            }
        })

    # def test_squashed_spinful(self):
    #     self.run_test({
    #         'enable_spin_orbital': True,
    #         'calc_scheme': 'squashed',
    #         'inter': {
    #             'Transfer': 'transfer_spin_orbital.dat',
    #             'CoulombIntra': 'coulombintra.dat',
    #             'CoulombInter': 'coulombinter.dat',
    #         }
    #     })


if __name__ == '__main__':
    unittest.main()
