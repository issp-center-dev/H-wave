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
                    'Geometry': 'geom.dat',
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
