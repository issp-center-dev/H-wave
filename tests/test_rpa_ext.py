#!/usr/bin/env python3

import os, sys
import unittest

import numpy as np

class RPAExternal:
    def __init__(self):
        pass

    def _get_epsilon(self, kx_array, ky_array, t, t1, H):
        epsilon_k = np.zeros((2, np.shape(kx_array)[0], np.shape(ky_array)[0]), dtype=np.complex128)
        for idx, kx in enumerate(kx_array):
            for idy, ky in enumerate(ky_array):
                epsilon_k[0][idx][idy] = 2.0*t*(np.cos(kx)+np.cos(ky))+2.0*t1*np.cos(kx+ky) - H
                epsilon_k[1][idx][idy] = 2.0*t*(np.cos(kx)+np.cos(ky))+2.0*t1*np.cos(kx+ky) + H
        return epsilon_k

    def _get_green(self, Nx, Ny, Nmat, myu, beta, epsilon_k):
        iomega = np.array([(2.0*i+1.0-Nmat)*np.pi for i in range(Nmat)])/beta
        green_kw = np.zeros((Nx, Ny, Nmat), dtype=np.complex128)
        for idx in range(Nx):
            for idy in range(Ny):
                for idz in range(Nmat):
                    green_kw[idx][idy][idz] = 1.0/(1J*iomega[idz]-(epsilon_k[idx][idy] - myu))
        return green_kw

    def _calc_chi0q(self, Nx, Ny, Nmat, green_a_kw, green_b_kw, myu, beta):
        from numpy.fft import fft
        #mat2tau FFT
        def _green_mat2tau_fft(green_rw):
            green_rt = np.zeros(green_rw.shape, dtype=np.complex64)        
            for idx in range(green_rt.shape[0]):
                for idy in range(green_rt.shape[1]):
                    green_rt[idx][idy] = fft(green_rw[idx][idy])
                    delta_omega = 1.0/Nmat -1.0
                    for idz in range(Nmat):
                        green_rt[idx][idy][idz] *= np.exp(-1J*np.pi*idz*delta_omega)
            return green_rt
    
        green_a_kt = _green_mat2tau_fft(green_a_kw)
        green_b_kt = _green_mat2tau_fft(green_b_kw)

        from numpy.fft import ifft2
        #q2r FFT
        def _green_q2r_fft(green_kt):
            green_rt = np.zeros(green_kt.shape, dtype=np.complex64)
            for idz in range(green_kt.shape[2]):
                g_q_t = green_kt[:,:,idz]
                g_r_t = ifft2(g_q_t)
                green_rt[:, :, idz] = g_r_t
            return green_rt
    
        green_a_rt = _green_q2r_fft(green_a_kt)
        green_b_rt = _green_q2r_fft(green_b_kt)
    
        #Calculate Chi(r, t)
        chi_rt = np.zeros((Nx, Ny, Nmat), dtype=np.complex64)
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
        chi_qt = np.zeros((Nx, Ny, Nmat), dtype=np.complex64)
        from numpy.fft import fft2
        for idz in range(Nmat):
            chir_t = chi_rt[:,:,idz]
            chiq_t = fft2(chir_t)
            chi_qt[:, :, idz] = chiq_t

        chi_qw = np.zeros((Nx, Ny, Nmat), dtype=np.complex64)
        #tau2mat FFT
        from numpy.fft import ifft
        delta_omega = -1.0
        for idx in range(Nx):
            for idy in range(Ny):
                for idz in range(Nmat):
                    chi_qt[idx][idy][idz] *= np.exp(1J*np.pi*idz*delta_omega)
                chi_qw[idx][idy] = ifft(chi_qt[idx][idy])
            
        #return chi0q
        return -chi_qw/beta

    def _calc_chi0q_exact(self, kx_array, ky_array, epsilon_k, myu, beta, omega=0.0):    
        fermi = 1.0/(1.0+np.exp(beta*(epsilon_k-myu)))
        chi0q_exact = np.zeros((np.shape(kx_array)[0], np.shape(ky_array)[0]))
        Nx = np.shape(kx_array)[0]
        Ny = np.shape(ky_array)[0]
        eta=1e-12
        for idqx, qx in enumerate(kx_array):
            for idqy, qy in enumerate(ky_array):
                for idkx, kx in enumerate(kx_array):
                    for idky, ky in enumerate(ky_array):
                        kx_plus_qx = (idqx+idkx)%Nx
                        ky_plus_qy = (idqy+idky)%Ny
                        chi0q_exact[idqx][idqy] += -(fermi[kx_plus_qx][ky_plus_qy]-fermi[idkx][idky])/(epsilon_k[kx_plus_qx][ky_plus_qy]- epsilon_k[idkx][idky]+1J*eta)
        return chi0q_exact/(Nx*Ny).real

    def calc_chi0q_direct(self, Nx, Ny, Nmat, idqx, idqy, green_a_kw, green_b_kw, myu, beta):
        chi0q_ = 0.0
        idmat=0
        for idx in range(Nx):
            for idy in range(Ny):
                for idz in range(Nmat):
                    chi0q_ += green_a_kw[int(idx+idqx)%Nx][int(idy+idqy)%Ny][int(idz+idmat)%Nmat]*green_b_kw[idx][idy][idz]
        return -chi0q_/(Nx*Ny*beta)

    def run(self, params):
        t     = params.get("t", 1.0)
        t1    = params.get("t1", 0.5)
        U     = params.get("U", 4.0)
        V     = params.get("V", 1.0)
        H     = params.get("H", 0.1)

        Lx    = params.get("Lx", 16)
        Ly    = params.get("Ly", 16)
        Nmat  = params.get("Nmat", 64)
        myu   = params.get("myu", 0.0)
        beta  = params.get("beta", 1.0)

        kx_array = np.linspace(0, 2.*np.pi, Lx, endpoint=False)
        ky_array = np.linspace(0, 2.*np.pi, Ly, endpoint=False)

        epsilon_k = self._get_epsilon(kx_array, ky_array, t, t1, H)

        green_kw = np.zeros((2, Lx, Ly, Nmat), dtype=np.complex64)

        chi0q = np.zeros((2, Lx, Ly, Nmat), dtype=np.complex64)
        for idx in range(2):
            green_kw[idx] = self._get_green(Lx, Ly, Nmat, myu, beta, epsilon_k[idx])
            chi0q[idx] = self._calc_chi0q(Lx, Ly, Nmat, green_kw[idx], green_kw[idx],  myu, beta)
    
        # with open("chiq_fft_beta{}.dat".format(beta), "w") as fw:
        #     pass

        chiq = np.zeros((Lx,Ly,2,2), dtype=np.complex128)

        I = np.eye(2)
        for idqx, kx in enumerate(kx_array):
            for idqy, ky in enumerate(ky_array):
                _chi0 = chi0q[:, idqx, idqy, Nmat//2].real
                #X0_mat = np.array([[_chi0[0], 0 ], [0, _chi0[1]]])
                X0_mat = np.matrix([[_chi0[0], 0 ], [0, _chi0[1]]])
                Vq = 2.0*V*(np.cos(kx)+np.cos(ky))
                Vss = Vq
                Vsbs = U+Vq
                V_mat = np.array([[Vss, Vsbs], [Vsbs, Vss]])
                X_mat = np.dot(np.linalg.inv(I + X0_mat*V_mat), X0_mat)

                chiq[idqx,idqy,:,:] = X_mat
                # fw.write("{} {} {} {} {} {}\n".format(kx, ky, X_mat[0][0], X_mat[0][1], X_mat[1][0], X_mat[1][1]))

        return chi0q, chiq


class TestRPAExternal(unittest.TestCase):
    def run_test(self, params = {}):

        Lx = params.get('Lx', 16)
        Ly = params.get('Ly', 16)
        Nmat = params.get('Nmat', 64)
        T = params.get('T', 0.1)
        H = params.get('H', 2.0)

        #----------------------------------------------------------------
        # input_dict = toml.load("tests/rpa/input.toml")
        info_log = {}
        info_mode = {
            'mode': 'RPA',
            'param': {
                'T': T,
                'mu': 0.0,
                'CellShape': [Lx,Ly,1],
                'SubShape': [1,1,1],
                'Nmat': Nmat,
                'coeff_extern': -H,
            },
            'calc_scheme': 'general',
        }
        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input',
                'interaction': {
                    'path_to_input': 'tests/rpa/input',
                    'Geometry': 'geom.dat',
                    'Transfer': 'transfer.dat',
                    'Extern': 'extern.dat',
                    'CoulombIntra': 'coulombintra.dat',
                    'CoulombInter': 'coulombinter.dat',
                },
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
            },
        }

        # # overwrite by params
        # info_file['input']['interaction'].update(params)

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

        params_ref = {}
        params_ref['Lx'] = Lx
        params_ref['Ly'] = Ly
        params_ref['Nmat'] = Nmat
        params_ref['beta'] = 1.0/T
        params_ref['H'] = H
        
        rpa_ref = RPAExternal()
        chi0q_ref, chiq_ref = rpa_ref.run(params_ref)

        #----------------------------------------------------------------
        # compare

        # chi0q
        self.assertEqual(chi0q.shape, (2,Nmat,Lx*Ly,1,1,1,1))

        x = chi0q.reshape(2,Nmat,Lx,Ly)
        y = np.transpose(chi0q_ref,(0,3,1,2))

        self.assertTrue(np.allclose(x,y), "chi0q")

        # chiq
        self.assertEqual(chiq.shape, (Nmat,Lx*Ly,2,2,2,2))
        
        x = chiq[Nmat//2].reshape(Lx,Ly,4,4)
        y = chiq_ref.reshape(Lx,Ly,2,2)

        self.assertTrue(np.allclose(x[:,:,0,0], y[:,:,0,0]), "chi0_00")
        self.assertTrue(np.allclose(x[:,:,3,3], y[:,:,1,1]), "chi0_11")
        self.assertTrue(np.allclose(x[:,:,0,3], y[:,:,0,1]), "chi0_01")
        self.assertTrue(np.allclose(x[:,:,3,0], y[:,:,1,0]), "chi0_10")

        pass

    def test_external(self):
        self.run_test()

if __name__ == '__main__':
    unittest.main()
