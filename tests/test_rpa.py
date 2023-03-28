#!/usr/bin/env python3

import os, sys
import unittest

import numpy as np

class RPAOneOrbital:
    def __init__(self):
        pass

    def _get_epsilon(self, kx_array, ky_array, t, t1):
        epsilon_k = np.zeros((np.shape(kx_array)[0], np.shape(ky_array)[0]))
        for idx, kx in enumerate(kx_array):
            for idy, ky in enumerate(ky_array):
                epsilon_k[idx][idy] = 2.0*t*(np.cos(kx)+np.cos(ky))+2.0*t1*np.cos(kx+ky)
        return epsilon_k

    def _get_green(self, Nx, Ny, Nmat, myu, beta, epsilon_k):
        iomega = np.array([(2.0*i+1.0-Nmat)*np.pi for i in range(Nmat)])/beta
        green_kw = np.zeros((Nx, Ny, Nmat), dtype=np.complex128)
        for idx in range(Nx):
            for idy in range(Ny):
                for idz in range(Nmat):
                    green_kw[idx][idy][idz] = 1.0/(1J*iomega[idz]-(epsilon_k[idx][idy] - myu))
        return green_kw

    def _calc_chi0q(self, Nx, Ny, Nmat, green_a_kw, green_b_kw, beta):
        from numpy.fft import ifft2
        #q2r FFT
        def _green_q2r_fft(green_kw):
            green_rw = np.zeros(green_kw.shape, dtype=np.complex128)
            for idz in range(green_kw.shape[2]):
                g_q_w = green_kw[:,:,idz]
                g_r_w = ifft2(g_q_w)
                green_rw[:, :, idz] = g_r_w
            return green_rw
    
        green_a_rw = _green_q2r_fft(green_a_kw)
        green_b_rw = _green_q2r_fft(green_b_kw)
        
        from numpy.fft import fft
        #mat2tau FFT
        def _green_mat2tau_fft(green_rw):
            green_rt = np.zeros(green_rw.shape, dtype=np.complex128)        
            for idx in range(green_rt.shape[0]):
                for idy in range(green_rt.shape[1]):
                    green_rt[idx][idy] = fft(green_rw[idx][idy])
                    delta_omega = 1.0/Nmat -1.0
                    for idz in range(Nmat):
                        green_rt[idx][idy][idz] *= np.exp(-1J*np.pi*idz*delta_omega)
            return green_rt

        green_a_rt = _green_mat2tau_fft(green_a_rw)
        green_b_rt = _green_mat2tau_fft(green_b_rw)
    
        #Calculate Chi(r, t)
        chi_rt = np.zeros((Nx, Ny, Nmat), dtype=np.complex128)
        for idz in range(Nmat):
            for idx in range(Nx):
                for idy in range(Ny):
                    idx_inv = (Nx-idx)%Nx
                    idy_inv = (Ny-idy)%Ny
                    idz_inv = (Nmat-idz)%Nmat
                    if idz == 0:
                        chi_rt[idx][idy][idz] = green_a_rt[idx][idy][idz] * green_b_rt[idx_inv][idy_inv][idz_inv]
                    else:
                        chi_rt[idx][idy][idz] = -green_a_rt[idx][idy][idz] * green_b_rt[idx_inv][idy_inv][idz_inv]

        #r2q FFT
        chi_qt = np.zeros((Nx, Ny, Nmat), dtype=np.complex128)
        from numpy.fft import fft2
        for idz in range(Nmat):
            chir_t = chi_rt[:,:,idz]
            chiq_t = fft2(chir_t)
            chi_qt[:, :, idz] = chiq_t

        chi_qw = np.zeros((Nx, Ny, Nmat), dtype=np.complex128)
        #tau2mat FFT
        from numpy.fft import ifft
        delta_omega = -1.0
        for idx in range(Nx):
            for idy in range(Ny):
                for idz in range(Nmat):
                    chi_qt[idx][idy][idz] *= np.exp(1J*np.pi*idz*delta_omega)
                chi_qw[idx][idy] = ifft(chi_qt[idx][idy])
        return -chi_qw/beta

    def run(self, params):
        t       = params.get("t",         1.0)
        t1      = params.get("t1",        0.5)
        U       = params.get("U",         0.0)
        V       = params.get("V",         0.0)
        J_hund  = params.get("J_hund",    0.0)
        J_ex    = params.get("J_ex",      0.0)
        J_ising = params.get("J_ising",   0.0)
        J_pair  = params.get("J_pair",    0.0)

        Lx      = params.get("Lx",        16)
        Ly      = params.get("Ly",        16)
        Nmat    = params.get("Nmat",      128)
        myu     = params.get("myu",       0.0)
        beta    = params.get("beta",      0.5)

        kx_array = np.linspace(0, 2.*np.pi, Lx, endpoint=False)
        ky_array = np.linspace(0, 2.*np.pi, Ly, endpoint=False)

        epsilon_k = self._get_epsilon(kx_array, ky_array, t, t1)

        green_kw = self._get_green(Lx, Ly, Nmat, myu, beta, epsilon_k)

        chi0q = self._calc_chi0q(Lx, Ly, Nmat, green_kw, green_kw, beta)

        # os.makedirs(output_dir, exist_ok=True)
        # with open(os.path.join(output_dir, "chiq_fft.dat"), "w") as fw:
        #     pass

        ham = np.zeros((Lx,Ly,4,4), dtype=np.float64)
        chi = np.zeros((Lx,Ly,4,4), dtype=np.float64)
        
        for idqx, kx in enumerate(kx_array):
            for idqy, ky in enumerate(ky_array):
                chi0 = chi0q[idqx][idqy][Nmat//2].real

                Vq = 2.0*V*(np.cos(kx)+np.cos(ky))
                J_hund_q = -2.0*J_hund*(np.cos(kx)+np.cos(ky))
                J_ex_q = -2.0*J_ex*(np.cos(kx)+np.cos(ky))
                J_ising_q = 2.0*J_ising*(np.cos(kx)+np.cos(ky))
                J_pair_q = 2.0*J_pair*(np.cos(kx)+np.cos(ky))
            
                I = np.identity(4)
                X0 = np.matrix([[chi0, 0, 0, 0],[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, chi0]])
                W = np.matrix([[Vq+J_hund_q+J_ising_q, 0, 0, U+Vq-J_ising_q],
#                               [0, J_pair_q, J_ex_q, 0],
#                               [0, J_ex_q, J_pair_q, 0],
                               [0, J_ex_q, J_pair_q, 0],
                               [0, J_pair_q, J_ex_q, 0],
                               [U+Vq-J_ising_q, 0, 0 , Vq+J_hund_q+J_ising_q]])
                X = (I+X0*W).I*X0

                ham[idqx,idqy,:,:] = W
                chi[idqx,idqy,:,:] = X
                
                chiss = X[0,0]
                chisbs = X[0,3]

                # fw.write("{} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(
                #     kx, ky,
                #     chi0q[idqx][idqy][Nmat//2].real,
                #     X[0, 0], X[0, 1], X[0, 2], X[0, 3],
                #     X[1, 1], X[1, 2], X[1, 3],
                #     X[2, 2], X[2, 3],
                #     X[3, 3]))            
            
        return chi0q, chi, ham, kx_array, ky_array

class TestRPA(unittest.TestCase):
    def run_test(self, params, params_ref):
        Lx, Ly, Nmat = 16,16,128

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
        }
        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input',
                'interaction': {
                    'path_to_input': 'tests/rpa/input',
                    'Geometry': 'geom.dat',
                    'Transfer': 'transfer.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
            },
        }

        # overwrite by params
        info_file['input']['interaction'].update(params)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")

        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        solver.solve(green_info, info_file['output']['path_to_output'])

        chi0q = green_info["chi0q"]
        chiq = green_info["chiq"]
        
        #----------------------------------------------------------------
        # reference
        rpaone = RPAOneOrbital()

        chi0q_ref, chiq_ref, ham_ref, kx_array, ky_array = rpaone.run(params_ref)

        #----------------------------------------------------------------
        # compare

        # chi0q
        chi0q_x = chi0q[:,:,0,0,0,0]
        chi0q_ref_x = np.transpose(chi0q_ref, (2,0,1)).reshape(Nmat,Lx*Ly)

        self.assertTrue(np.allclose(chi0q_x, chi0q_ref_x), 'chi0q')

        # hamiltonian
        ham = solver.ham_info.ham_inter_q.reshape(Lx,Ly,4,4)
        self.assertTrue(np.allclose(ham, ham_ref), 'hamiltonian')

        # chiq
        chiq_x = chiq[Nmat//2].reshape(Lx,Ly,4,4)

        self.assertTrue(np.allclose(chiq_x, chiq_ref), 'chiq')

    def test_U_and_V(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat', 'CoulombInter': 'coulombinter.dat' },
            { 't1': 0.5, 'U': 4.0, 'V': 1.0 }
        )
        
    def test_U_and_Hund(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat', 'Hund': 'coulombinter.dat' },
            { 't1': 0.5, 'U': 4.0, 'J_hund': 1.0 }
        )
        
    def test_U_and_Exchange(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat', 'Exchange': 'coulombinter.dat' },
            { 't1': 0.5, 'U': 4.0, 'J_ex': 1.0 }
        )
        
    def test_U_and_Ising(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat', 'Ising': 'coulombinter.dat' },
            { 't1': 0.5, 'U': 4.0, 'J_ising': 1.0 }
        )
        
    def test_U_and_PairLift(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat', 'PairLift': 'coulombinter.dat' },
            { 't1': 0.5, 'U': 4.0, 'J_pair': 1.0 }
        )
        
    def test_U_and_all(self):
        self.run_test(
            { 'CoulombIntra': 'coulombintra.dat',
              'CoulombInter': 'coulombinter.dat',
              'Hund': 'coulombinter.dat',
              'Exchange': 'coulombinter.dat',
              'Ising': 'coulombinter.dat',
              'PairLift': 'coulombinter.dat',
             },
            { 't1': 0.5, 'U': 4.0,
              'V': 1.0,
              'J_hund': 1.0,
              'J_ex': 1.0,
              'J_ising': 1.0,
              'J_pair': 1.0,
             }
        )

if __name__ == '__main__':
    unittest.main()
