"""UHFk spin-orbital ``_make_ham`` against the interaction loop it replaced
(#83): the SO branch now runs the shared ``hartree_fock.accumulate_hf`` in
the virtual spin-block frame. The pre-refactor loop and its frame
conversions are carried here VERBATIM, and the two Hamiltonians must agree
bit for bit for every interaction type, with and without the Fock term.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

NP = 3
NSO = 2 * NP


def _legacy_so_to_virtual_green(self, G_so):
    nvol = G_so.shape[0]
    norb_phys = self.norb_phys
    G_flat = G_so.reshape(nvol, self.nd, self.nd)
    G_virt = np.zeros((nvol, 2, norb_phys, 2, norb_phys), dtype=np.complex128)
    for s in range(2):
        for t in range(2):
            G_virt[:, s, :, t, :] = G_flat[:, s::2, t::2]
    return G_virt


def _legacy_virtual_ham_to_so(self, H_virt):
    nvol = H_virt.shape[0]
    norb_phys = self.norb_phys
    nd = self.nd
    H_so = np.zeros((nvol, nd, nd), dtype=np.complex128)
    for s in range(2):
        for t in range(2):
            H_so[:, s::2, t::2] = H_virt[:, s, :, t, :]
    return H_so


def _legacy_make_ham(self):
    """The enable_spin_orbital branch of UHFk._make_ham before #83."""
    nx, ny, nz = self.shape
    nvol = self.nvol
    nd = self.nd
    norb_inter = self.norb_phys
    nd_virt = 2 * norb_inter
    gab_r = _legacy_so_to_virtual_green(self, self.Green)
    gbb = np.diagonal(gab_r, axis1=2, axis2=4)[0, :, :, :]
    ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
    ham += self.ham_trans
    if self.enable_spin_orbital:
        for type in ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising', 'PairLift', 'Exchange']:
            if self.inter_table[type] is not None:
                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb_inter,norb_inter)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # non-cross term
                #   sum_r J_{ab}(r) G_{bb,uv}(0) Spin{s,u,v,t}
                hh0 = np.einsum('uvb, suvt -> stb', gbb, spin)
                hh1 = np.einsum('rab, stb -> rsta', jab_r, hh0)
                hh2 = np.einsum('rsta, ab -> rsatb', hh1, np.eye(norb_inter, norb_inter))
                hh3 = np.sum(hh2, axis=0)  # shape: (2, norb_inter, 2, norb_inter)

                if self.enable_spin_orbital:
                    hh3_nd = _legacy_virtual_ham_to_so(self,
                        np.broadcast_to(hh3.reshape(1, 2, norb_inter, 2, norb_inter),
                                       (nvol, 2, norb_inter, 2, norb_inter)).copy()
                    )
                else:
                    hh3_nd = np.broadcast_to(hh3.reshape(nd, nd), (nvol, nd, nd))

                ham += hh3_nd

                # cross term
                #   - sum_r J_{ab}(r) G_{ba,uv}(r) Spin{s,u,t,v} e^{ikr}
                if self.iflag_fock:
                    hh4 = np.einsum('rab, rubva, sutv -> rsatb', jab_r, gab_r, spin, optimize=True)

                    #   fourier transform: sum_r (*) e^{ikr}
                    hh5 = np.fft.ifftn(hh4.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')

                    if self.enable_spin_orbital:
                        ham -= _legacy_virtual_ham_to_so(self,
                            hh5.reshape(nvol, 2, norb_inter, 2, norb_inter)
                        )
                    else:
                        ham -= hh5.reshape(nvol, nd, nd)

        # interaction term: PairHop type
        for type in ['PairHop']:
            if self.inter_table[type] is not None:
                # coefficient of interaction term J_{ab}(r)
                jab_r = self.inter_table[type].reshape(nvol,norb_inter,norb_inter)
                # and its spin combination  Spin(s1,s2,s3,s4)
                spin = self.spin_table[type]

                # non-cross and cross term
                #   + sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,v,t} e^{ikr}
                #   - sum_r J_{ab}(r) G_{ab,uv}(-r) Spin{s,u,t,v} e^{ikr}

                if self.iflag_fock:
                    hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                    hh2 = np.einsum('rvbua, sutv -> rsbta', np.conjugate(gab_r), spin)
                    hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, (hh1 - hh2))
                    hh4 = np.fft.ifftn(hh3.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')
                else:
                    hh1 = np.einsum('rvbua, suvt -> rsbta', np.conjugate(gab_r), spin)
                    hh3 = np.einsum('rab, rsbta -> rsatb', jab_r, hh1)
                    hh4 = np.fft.ifftn(hh3.reshape(nx,ny,nz,nd_virt,nd_virt), axes=(0,1,2), norm='forward')

                if self.enable_spin_orbital:
                    ham += _legacy_virtual_ham_to_so(self,
                        hh4.reshape(nvol, 2, norb_inter, 2, norb_inter)
                    )
                else:
                    ham += hh4.reshape(nvol, nd, nd)


    if len(self.block_info) > 1:
        mask = np.zeros((nd, nd), dtype=bool)
        for blk in self.block_info:
            idx = np.array(blk)
            mask[np.ix_(idx, idx)] = True
        ham[:, ~mask] = 0.0
    return ham


def _write(d):
    rng = np.random.default_rng(11)
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % NSO)
        f.write("0 0 0\n" * NSO)
    H0 = rng.normal(size=(NSO, NSO)) + 1j * rng.normal(size=(NSO, NSO))
    H0 = 0.2 * (H0 + H0.conj().T)
    T = 0.3 * (rng.normal(size=(NSO, NSO)) + 1j * rng.normal(size=(NSO, NSO))) + np.eye(NSO)
    rows = []
    for a in range(NSO):
        for b in range(NSO):
            rows += [(0, a, b, H0[a, b]), (1, a, b, T[a, b]), (-1, a, b, T[b, a].conjugate())]
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n%d\n1\n 1\n" % NSO)
        for R, a, b, v in rows:
            f.write(" %d 0 0 %d %d %.12f %.12f\n" % (R, a + 1, b + 1, v.real, v.imag))
    files = {}
    offd = [(a, b) for a in range(NP) for b in range(NP) if a != b]
    spec = {"CoulombIntra": ([(a, a, 1.0 + 0.1 * a) for a in range(NP)], []),
            "CoulombInter": ([(a, b, 0.5) for a, b in offd],
                             [(a, b, 0.2 + 0.05 * a + 0.01 * b) for a in range(NP) for b in range(NP)])}
    for nm, v in (("Hund", 0.15), ("Exchange", 0.12), ("PairHop", 0.07),
                  ("Ising", 0.05), ("PairLift", 0.03)):
        spec[nm] = ([(a, b, v) for a, b in offd], [(a, b, 0.5 * v) for a, b in offd])
    for nm, (on, off) in spec.items():
        fn = nm.lower() + ".dat"
        with open(os.path.join(d, fn), "w") as f:
            f.write("%s\n%d\n1\n 1\n" % (nm, NP))
            for a, b, v in on:
                f.write(" 0 0 0 %d %d %.12f 0.0\n" % (a + 1, b + 1, v))
            for a, b, v in off:
                f.write(" 1 0 0 %d %d %.12f 0.0\n" % (a + 1, b + 1, v))
                f.write(" -1 0 0 %d %d %.12f 0.0\n" % (b + 1, a + 1, v))
        files[nm] = fn
    return files


class TestSOMakeHamMatchesTheLegacyLoop(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="uhfk_so_legacy_")
        self.files = _write(self.d)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _solver(self, fock, sub):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.uhfk as uhfk
        inter = {"path_to_input": self.d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        inter.update(self.files)
        io = read_input_k.QLMSkInput({"path_to_input": self.d, "interaction": inter})
        mode = {"mode": "UHFk", "enable_spin_orbital": True, "flag_fock": fock,
                "param": {"Ncond": NSO, "T": 0.2, "CellShape": [6, 1, 1],
                          "SubShape": [sub, 1, 1], "RndSeed": 3}}
        s = uhfk.UHFk(io.get_param("ham"), {"print_level": 0}, mode)
        s._make_ham_trans()
        s._make_ham_inter()
        s._detect_blocks()
        s.Green = s._initial_green(io.get_param("green"))
        rng = np.random.default_rng(5)
        g = rng.normal(size=s.Green.shape) + 1j * rng.normal(size=s.Green.shape)
        s.Green = g * 10.0 ** rng.uniform(-6, 2, size=g.shape)   # wide dynamic range
        return s

    def test_bit_identical(self):
        for fock in (True, False):
            for sub in (1, 2):
                with self.subTest(fock=fock, sub=sub):
                    s = self._solver(fock, sub)
                    want = _legacy_make_ham(s)
                    s._make_ham()
                    np.testing.assert_array_equal(s.ham, want)


if __name__ == "__main__":
    unittest.main()
