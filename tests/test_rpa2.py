#!/usr/bin/env python3

import os, sys
import unittest

import numpy as np
import toml

class RPATwoOrbital:
    def __init__(self):
        pass

    def _get_epsilon(self, kx_array, ky_array, t, t1):
        epsilon_k = np.zeros((2, 2, np.shape(kx_array)[0], np.shape(ky_array)[0]), dtype=np.complex128)
        for idx, kx in enumerate(kx_array):
            for idy, ky in enumerate(ky_array):
                epsilon_k[0][0][idx][idy] = 2.0*t*np.cos(ky)
                epsilon_k[1][1][idx][idy] = 2.0*t*np.cos(ky)
                epsilon_k[0][1][idx][idy] = t1*(1.0+np.exp(1J*kx))
                epsilon_k[1][0][idx][idy] = t1*(1.0+np.exp(-1J*kx))
        #XXX
        epsilon_k[np.abs(epsilon_k) < 1.0e-8] = 0
                
        return epsilon_k

    def _get_green(self, Nx, Ny, Nmat, myu, beta, epsilon_k):
        iomega = np.array([(2.0*i+1.0-Nmat)*np.pi for i in range(Nmat)])/beta
        green_kw = np.zeros((2, 2, 2, Nx, Ny, Nmat), dtype=np.complex128)
        #diagonalize epsilon_k

        #XXX
        ew = np.zeros((Nx,Ny,2),dtype=np.complex128)
        ev = np.zeros((Nx,Ny,2,2),dtype=np.complex128)

        for idx in range(Nx):
            for idy in range(Ny):
                #eig = np.linalg.eig(epsilon_k[:,:,idx,idy])
                eig = np.linalg.eigh(epsilon_k[:,:,idx,idy])
                ene = eig[0]
                vec = eig[1]

                #XXX
                ew[idx,idy] = eig[0]
                ev[idx,idy] = eig[1]
                
                for iorb in range(2):
                    for jorb in range(2):
                        for igamma in range(2):
                            factor = vec[iorb, igamma]*np.conj(vec[jorb, igamma])
                            ene_gamma = ene[igamma]
                            for idz in range(Nmat):
                                green_kw[igamma][iorb][jorb][idx][idy][idz] = factor/(1J*iomega[idz]-(ene_gamma - myu))

        #XXX
        self.ew = ew
        self.ev = ev
        return np.sum(green_kw, axis=0)

    def _calc_chi0q(self, Nx, Ny, Nmat, green_a_kw, green_b_kw, myu, beta):
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

        #return chi0q
        return -chi_qw/beta

    def run(self, params):
        t     = params.get("t", 1.0)
        t1    = params.get("t1", 0.5)
        U     = params.get("U", 4.0)
        V     = params.get("V", 1.0)

        Lx    = params.get("Lx", 8)
        Ly    = params.get("Ly", 16)
        Nmat  = params.get("Nmat", 64)
        myu   = params.get("myu", 0.0)
        beta  = params.get("beta", 0.5)

        kx_array = np.linspace(0, 2.*np.pi, Lx, endpoint=False)
        ky_array = np.linspace(0, 2.*np.pi, Ly, endpoint=False)

        epsilon_k = self._get_epsilon(kx_array, ky_array, t, t1)

        # eig = np.linalg.eig(epsilon_k[:,:,1,0])

        green_kw = self._get_green(Lx, Ly, Nmat, myu, beta, epsilon_k)

        # Traditional way
        # chi0q = np.zeros((2, 2, Lx, Ly, Nmat), dtype=np.complex128)
        # for iorb in range(2):
        #     for jorb in range(2):
        #         g1 = np.zeros((Lx, Ly, Nmat), dtype = np.complex128)
        #         g2 = np.zeros((Lx, Ly, Nmat), dtype = np.complex128)
        #         for igamma in range(2):
        #             g1 += green_kw[igamma][iorb][jorb]
        #             g2 += green_kw[igamma][jorb][iorb]
        #         chi0q[iorb][jorb] = calc_chi0q(Lx, Ly, Nmat, g1,  g2,  myu, beta)

        chi0q = np.zeros((2, 2, 2, 2, Lx, Ly, Nmat), dtype=np.complex128)
        for iorb in range(2):
            for jorb in range(2):
                for iorb1 in range(2):
                    for jorb1 in range(2):
                        chi0q[iorb][iorb1][jorb][jorb1] += self._calc_chi0q(Lx, Ly, Nmat, green_kw[iorb][jorb], green_kw[jorb1][iorb1], myu, beta)

        # with open("chiq_fft_beta{}.dat".format(beta), "w") as fw:
        #     pass

        ham = np.zeros((Lx,Ly,8,8), dtype=np.complex128)
        chi = np.zeros((Lx,Ly,8,8), dtype=np.complex128)

        for idqx, kx in enumerate(kx_array):
            for idqy, ky in enumerate(ky_array):
                I = np.identity(8, dtype=np.complex128)
                #chi0 = chi0q[:,:,:,:, idqx, idqy,int(Nmat/2)].real
                chi0 = chi0q[:,:,:,:, idqx, idqy,Nmat//2]
                #Vxq = V*(1.0+np.exp(-1J*kx))
                Vxq = V*(1.0+np.exp(+1J*kx))
                Vyq = 2.0*V*np.cos(ky)
                # X0 = np.matrix(([chi0[0][0][0][0], 0, 0, chi0[0][0][1][1], 0, 0, 0, 0],
                #                 [0, 0, 0, 0, 0, 0, 0, 0],
                #                 [0, 0, 0, 0, 0, 0, 0, 0],
                #                 [chi0[1][1][0][0], 0, 0, chi0[1][1][1][1], 0, 0, 0, 0],
                #                 [0, 0, 0, 0, chi0[0][0][0][0], 0, 0, chi0[0][0][1][1]],
                #                 [0, 0, 0, 0, 0, 0, 0, 0],
                #                 [0, 0, 0, 0, 0, 0, 0, 0],
                #                 [0, 0, 0, 0, chi0[1][1][0][0], 0, 0, chi0[1][1][1][1]]), dtype=np.complex128)
                X0 = np.zeros((8,8),dtype=np.complex128)
                for ia in range(4):
                    ia1 = ia//2
                    ia2 = ia%2
                    for ib in range(4):
                        ib1 = ib//2
                        ib2 = ib%2
                        X0[ia,ib] = chi0[ia1,ia2,ib1,ib2]
                        X0[ia+4,ib+4] = chi0[ia1,ia2,ib1,ib2]
                
                W =  np.matrix([[Vyq, 0, 0, np.conj(Vxq), U+Vyq, 0, 0, np.conj(Vxq)],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [Vxq, 0, 0, Vyq, Vxq, 0, 0, U+Vyq],
                                [U+Vyq, 0, 0, np.conj(Vxq), Vyq, 0, 0, np.conj(Vxq)],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [0, 0, 0, 0, 0, 0, 0, 0],
                                [Vxq, 0, 0, U+Vyq, Vxq, 0, 0, Vyq]], dtype=np.complex128)
                X = (I+X0@W).I*X0

                ham[idqx,idqy,:,:] = W
                chi[idqx,idqy,:,:] = X

                # fw.write("{} {} ".format(kx, ky))
                # for i in range(8):
                #     for j in range(8):
                #         fw.write("{} ".format(X[i,j].real))
                # fw.write("\n")

        return chi0q, chi, ham, kx_array, ky_array, green_kw, epsilon_k


class TestRPATwoOrbital(unittest.TestCase):
    def run_test(self, params, params_ref):

        Lx, Ly, Nmat = 8,16,64
        T = 1.0

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
            },
        }
        info_file = {
            'input': {
                'path_to_input': 'tests/rpa/input_2orb',
                'interaction': {
                    'path_to_input': 'tests/rpa/input_2orb',
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

        params_ref['Lx'] = params_ref.get('Lx', Lx)
        params_ref['Ly'] = params_ref.get('Ly', Ly)
        params_ref['Nmat'] = params_ref.get('Nmat', Nmat)
        params_ref['beta'] = params_ref.get('beta', 1.0/T)
        
        rpatwo = RPATwoOrbital()
        chi0q_ref, chiq_ref, ham_ref, kxarr, kyarr, green_ref, eps_ref = rpatwo.run(params_ref)

        #----------------------------------------------------------------
        # compare

        # epsilon_k
        self.assertEqual(solver.ham_info.ham_trans_q.shape, (Lx,Ly,1,2,2), "epsk.shape")
        
        epsk = solver.ham_info.ham_trans_q.reshape(Lx,Ly,2,2)
        epsk_ref = np.transpose(eps_ref, (2,3,0,1))

        self.assertTrue(np.allclose(epsk,epsk_ref),"epsk")

        # eigenvalues and eigenvectors
        ew = solver.H0_eigenvalue.reshape(Lx,Ly,2)
        ev = solver.H0_eigenvector.reshape(Lx,Ly,2,2)

        ew_ref = rpatwo.ew
        ev_ref = rpatwo.ev

        for ix in range(Lx):
            for iy in range(Ly):
                eig = ew[ix,iy]
                evec = ev[ix,iy]
                eig_ref = ew_ref[ix,iy]
                evec_ref = ev_ref[ix,iy]

                if np.isclose(eig_ref[0],eig_ref[1]):
                    ng0,ng1 = 0,0
                    for ig in range(2):
                        e0 = evec_ref[:,0]
                        e1 = evec_ref[:,1]
                        ee = evec[:,ig]
                        if not np.isclose(eig[ig],eig_ref[0]):
                            print("eigenvalue mismatch",ix,iy,eig[ig],eig_ref[0])
                            self.assertTrue(False,"eig mismatch")
                        else:
                            if np.allclose(ee,e0) or np.allclose(ee,-e0):
                                ng0 += 1
                            elif np.allclose(ee,e1) or np.allclose(ee,-e1):
                                ng1 += 1
                            else:
                                print(f"eigenvector {ig} not match",ix,iy,ee,e0,e1)
                                self.assertTrue(False,"evec mismatch")
                    if ng0 == 1 and ng1 == 1:
                        pass
                    else:
                        print("eigenvector set mismatch",ix,iy,ng0,ng1)
                        self.assertTrue(False,"ng0ng1")
                else:
                    ng0,ng1 = 0,0
                    for ig in range(2):
                        e0 = evec_ref[:,0]
                        e1 = evec_ref[:,1]
                        ee = evec[:,ig]
                        if np.isclose(eig[ig],eig_ref[0]):
                            if np.allclose(ee,e0) or np.allclose(ee,-e0):
                                ng0 += 1
                            else:
                                print(f"eigenvector {ig} not match for vec0",ix,iy,ee,e0)
                                self.assertTrue(False,"evec mismatch 0")
                        elif np.isclose(eig[ig],eig_ref[1]):
                            if np.allclose(ee,e1) or np.allclose(ee,-e1):
                                ng1 += 1
                            else:
                                print(f"eigenvector {ig} not match for vec1",ix,iy,ee,e1)
                                self.assertTrue(False,"evec mismatch 1")
                        else:
                            print(f"eigenvalue {ig} not found",ix,iy,eig[ig],eig_ref[0],eig_ref[1])
                            self.assertTrue(False,"eig mismatch")
                    if ng0 == 1 and ng1 == 1:
                        pass
                    else:
                        print("eigenvector set mismatch",ix,iy,ng0,ng1)
                        self.assertTrue(False,"evec ng0ng1")

        # green function
        green = solver.green0

        x=green.reshape(Nmat,Lx,Ly,2,2)
        # green_ref[a][b][kx][ky][wn] -> [wn,kx,ky,a,b]
        y=np.transpose(green_ref,(4,2,3,0,1))

        self.assertTrue(np.allclose(x,y), "green function")

        # chi0q
        x = chi0q.reshape(Nmat,Lx,Ly,2,2,2,2)
        y = np.transpose(chi0q_ref,(6,4,5,0,1,2,3))

        self.assertTrue(np.allclose(x,y), "chi0q")

        # hamiltonian
        x = solver.ham_info.ham_inter_q.reshape(Lx,Ly,16,16)

        # W^{aa's,bb't} \delta{ss'} \delta{tt'} -> W^{(as)(a's'),(bt)(b't')}
        y = np.einsum('ijsabtcd,su,tv->ijsaubvctd',
                      ham_ref.reshape(Lx,Ly,*(2,2,2)*2),
                      np.identity(2),
                      np.identity(2)).reshape(Lx,Ly,16,16)

        self.assertTrue(np.allclose(x,y), "hamiltonian")

        # chiq
        x = chiq[Nmat//2].reshape(Lx,Ly,*(2,2)*4)
        #XXX
        #x[np.abs(x)<1.0e-8]=0
        
        y = np.einsum('ijsabtcd,su,tv->ijsaubvctd',
                      chiq_ref.reshape(Lx,Ly,*(2,2,2)*2),
                      np.identity(2),
                      np.identity(2)).reshape(Lx,Ly,*(2,2)*4)
        #XXX
        #y[np.abs(y)<1.0e-8]=0

        #for ss in [ [0,0,0,0], [1,1,1,1], [0,1,0,1], [0,1,1,0], [0,0,1,1] ]:
        for s in range(16):
            ss = [ s%2, (s//2)%2, (s//4)%2, (s//8)%2 ]
            # print(ss)
            for ix in range(Lx):
                for iy in range(Ly):
                    xx = x[ix,iy,ss[0],:,ss[1],:,ss[2],:,ss[3],:].reshape(4,4)
                    yy = y[ix,iy,ss[0],:,ss[1],:,ss[2],:,ss[3],:].reshape(4,4)
                    # for ia in range(4):
                    #     for ib in range(4):
                    #         print(ix,iy,ia,ib,xx[ia,ib],yy[ia,ib])
                    self.assertTrue(np.allclose(xx,yy), "chi_{}{}{}{} k={},{}".format(*ss,ix,iy))

        pass

    def test_2orbital(self):
        self.run_test({
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
        }, {})
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
