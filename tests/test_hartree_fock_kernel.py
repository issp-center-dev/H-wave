"""Bit-identity oracle for the Hartree-Fock kernel extracted from UHFk
(#181 Phase B, spec 2026-09-06 section 2.2).

``_LegacyUHFk`` carries the PRE-refactor bodies of ``UHFk._make_ham_inter`` and
``UHFk._make_ham`` VERBATIM (extracted mechanically from develop 703ca32c before
the refactor; do not edit them). The new ``hwave.solver.hartree_fock`` kernel
must reproduce their arithmetic bit for bit on every interaction type, with
and without the Fock term, on complex, non-symmetric paramagnetic and
non-paramagnetic densities.
"""
import logging
import unittest

import numpy as np

from hwave.qlmsio import wan90
from hwave.solver.kgrid import reverse_fft_axes

logger = logging.getLogger(__name__)


def exit(code):  # the legacy body calls exit(1) on the ambiguous Coulomb case
    raise SystemExit(code)


class _LegacyUHFk:
    """Just enough state for the two legacy methods (normal mode only)."""

    def __init__(self, param_ham, norb, shape, iflag_fock):
        self.param_ham = param_ham
        self.shape = tuple(shape)
        self.nvol = int(np.prod(shape))
        self.norb = norb
        self.ns = 2
        self.nd = 2 * norb
        self.enable_spin_orbital = False
        self.iflag_fock = iflag_fock
        self.block_info = [list(range(self.nd))]
        self.ham_trans = np.zeros((self.nvol, self.nd, self.nd), dtype=np.complex128)

    def _make_ham_inter(self):
        logger.debug(">>> _make_ham_inter")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        # In spin-orbital mode, interaction arrays use physical orbital count
        norb     = self.norb_phys if self.enable_spin_orbital else self.norb
        nd       = self.nd

        #----------------
        # interaction table
        #----------------
        self.inter_table = {}
        self.spin_table = {}

        #----------------
        # Coulomb Intra and Coulomb Inter
        #----------------
        if 'Coulomb' in self.param_ham.keys():
            # The aggregate 'Coulomb' input already provides both the intra and
            # inter parts; combining it with explicit CoulombIntra/CoulombInter
            # is ambiguous (the explicit terms would be silently dropped).
            if ('CoulombIntra' in self.param_ham.keys()
                    or 'CoulombInter' in self.param_ham.keys()):
                logger.error(
                    "Coulomb cannot be specified together with "
                    "CoulombIntra or CoulombInter")
                exit(1)

            # assume zvo_ur.dat
            # divide into r=0 (coulomb intra) and r!=0 (coulomb inter)
            # via the decomposition shared with RPA/FLEX
            coulomb_intra, coulomb_inter = wan90.split_coulomb(
                self.param_ham["Coulomb"])

            # coulomb intra: diagonal part of r=0 cell
            uab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            # coulomb inter: off-diagonal part
            vab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in coulomb_intra.items():
                uab_r[(*irvec,*orbvec)] += v
            for (irvec,orbvec), v in coulomb_inter.items():
                vab_r[(*irvec,*orbvec)] += v

            # coulomb intra
            # interaction coeffs
            self.inter_table["CoulombIntra"] = uab_r # r=0 component
            # spin combination
            self.spin_table["CoulombIntra"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["CoulombIntra"][0,1,1,0] = 1
            self.spin_table["CoulombIntra"][1,0,0,1] = 1

            #XXX
            # n_up n_up = n_up -> include in one-body term

            # coulomb inter
            # J~ab(r) = Jab(r) + Jba(-r)
            vba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(vab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["CoulombInter"] = (vab_r + vba)/2
            # spin combination
            self.spin_table["CoulombInter"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["CoulombInter"][0,0,0,0] = 1
            self.spin_table["CoulombInter"][1,1,1,1] = 1
            self.spin_table["CoulombInter"][0,1,1,0] = 1
            self.spin_table["CoulombInter"][1,0,0,1] = 1

        else:
            self.inter_table["CoulombIntra"] = None
            self.inter_table["CoulombInter"] = None

            # coulomb intra
            if 'CoulombIntra' in self.param_ham.keys():
                uab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

                # only r=0 and a=b component
                for (irvec,orbvec), v in self.param_ham["CoulombIntra"].items():
                    alpha, beta = orbvec
                    if irvec == (0,0,0) and alpha == beta:
                        uab_r[(*irvec, *orbvec)] += v

                # interaction coeffs
                self.inter_table["CoulombIntra"] = uab_r
                # spin combination
                self.spin_table["CoulombIntra"] = np.zeros((2,2,2,2), dtype=int)
                self.spin_table["CoulombIntra"][0,1,1,0] = 1
                self.spin_table["CoulombIntra"][1,0,0,1] = 1

            if 'CoulombInter' in self.param_ham.keys():
                vab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

                for (irvec,orbvec), v in self.param_ham["CoulombInter"].items():
#                    if irvec != (0,0,0):
                        vab_r[(*irvec, *orbvec)] += v

                vba = np.conjugate(
                    np.transpose(
                        reverse_fft_axes(vab_r, (0, 1, 2)),
                        (0,1,2,4,3)
                    )
                )

                # interaction coeffs
                self.inter_table["CoulombInter"] = (vab_r + vba)/2
                # spin combination
                self.spin_table["CoulombInter"] = np.zeros((2,2,2,2), dtype=int)
                self.spin_table["CoulombInter"][0,0,0,0] = 1
                self.spin_table["CoulombInter"][1,1,1,1] = 1
                self.spin_table["CoulombInter"][0,1,1,0] = 1
                self.spin_table["CoulombInter"][1,0,0,1] = 1

        #----------------
        # Hund
        #----------------
        if 'Hund' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Hund"].items():
                jab_r[(*irvec, *orbvec)] += v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(jab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs : -J^{Hund} by convention
            self.inter_table["Hund"] = -(jab_r + jba)/2
            # spin combination
            self.spin_table["Hund"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Hund"][0,0,0,0] = 1
            self.spin_table["Hund"][1,1,1,1] = 1
        else:
            self.inter_table["Hund"] = None

        #----------------
        # Ising
        #----------------
        if 'Ising' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Ising"].items():
                jab_r[(*irvec, *orbvec)] += v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(jab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs -- no 1/4: the documented Hamiltonian is
            # J (n_up - n_down)(n_up - n_down), and the RPA/FLEX vertex
            # content was adjudicated by exact diagonalization against
            # exactly that operator (#106); the historical /4 read the
            # file as J S^z S^z instead, so the same Ising file meant
            # couplings differing by 4 between UHFk and RPA/FLEX
            self.inter_table["Ising"] = (jab_r + jba)/2
            # spin combination
            self.spin_table["Ising"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Ising"][0,0,0,0] = 1
            self.spin_table["Ising"][1,1,1,1] = 1
            self.spin_table["Ising"][0,1,1,0] = -1
            self.spin_table["Ising"][1,0,0,1] = -1
        else:
            self.inter_table["Ising"] = None

        #----------------
        # PairLift
        #----------------
        if 'PairLift' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["PairLift"].items():
                jab_r[(*irvec, *orbvec)] += v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(jab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["PairLift"] = (jab_r + jba)/2
            # spin combination
            self.spin_table["PairLift"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["PairLift"][0,0,1,1] = 1
            self.spin_table["PairLift"][1,1,0,0] = 1
        else:
            self.inter_table["PairLift"] = None

        #----------------
        # Exchange
        #----------------
        if 'Exchange' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["Exchange"].items():
                jab_r[(*irvec, *orbvec)] += v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(jab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs : -J^{Ex} by convention
            self.inter_table["Exchange"] = -(jab_r + jba)/2
            # spin combination
            self.spin_table["Exchange"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["Exchange"][0,1,0,1] = 1
            self.spin_table["Exchange"][1,0,1,0] = 1
        else:
            self.inter_table["Exchange"] = None

        #----------------
        # PairHop
        #----------------
        if 'PairHop' in self.param_ham.keys():
            jab_r = np.zeros((nx,ny,nz,norb,norb), dtype=np.complex128)

            for (irvec,orbvec), v in self.param_ham["PairHop"].items():
                jab_r[(*irvec, *orbvec)] += v

            # J~ab(r) = Jab(r) + Jba(-r)
            jba = np.conjugate(
                np.transpose(
                    reverse_fft_axes(jab_r, (0, 1, 2)),
                    (0,1,2,4,3)
                )
            )

            # interaction coeffs
            self.inter_table["PairHop"] = (jab_r + jba)/2
            # spin combination
            self.spin_table["PairHop"] = np.zeros((2,2,2,2), dtype=int)
            self.spin_table["PairHop"][0,1,1,0] = 1
            self.spin_table["PairHop"][1,0,0,1] = 1
        else:
            self.inter_table["PairHop"] = None

    def _make_ham(self):
        logger.debug(">>> _make_ham")

        nx,ny,nz = self.shape
        nvol     = self.nvol
        nd       = self.nd
        ns       = self.ns

        # In spin-orbital mode, use virtual (ns=2, norb_phys) form for interactions
        if self.enable_spin_orbital:
            norb_inter = self.norb_phys
            nd_virt = 2 * norb_inter
            gab_r = self._so_to_virtual_green(self.Green)
        else:
            norb_inter = self.norb
            nd_virt = nd
            gab_r = self.Green

        # green function G_{ab,st}(r) : gab_r(r,s,a,t,b)
        # diagonal part G_{bb,st}(r=0) : gbb(s,t,b)
        gbb = np.diagonal(gab_r, axis1=2, axis2=4)[0,:,:,:]

        # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
        ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

        # transfer term  T_{ab}(k) (note convention)
        logger.debug("Transfer")
        ham += self.ham_trans

        # interaction term: Coulomb type
        for type in ['CoulombIntra', 'CoulombInter', 'Hund', 'Ising', 'PairLift', 'Exchange']:
            if self.inter_table[type] is not None:
                logger.debug(type)

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
                    hh3_nd = self._virtual_ham_to_so(
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
                        ham -= self._virtual_ham_to_so(
                            hh5.reshape(nvol, 2, norb_inter, 2, norb_inter)
                        )
                    else:
                        ham -= hh5.reshape(nvol, nd, nd)

        # interaction term: PairHop type
        for type in ['PairHop']:
            if self.inter_table[type] is not None:
                logger.debug(type)

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
                    ham += self._virtual_ham_to_so(
                        hh4.reshape(nvol, 2, norb_inter, 2, norb_inter)
                    )
                else:
                    ham += hh4.reshape(nvol, nd, nd)

        # Enforce block structure: zero out cross-block entries
        if len(self.block_info) > 1:
            mask = np.zeros((nd, nd), dtype=bool)
            for blk in self.block_info:
                idx = np.array(blk)
                ix = np.ix_(idx, idx)
                mask[ix] = True
            ham[:, ~mask] = 0.0

        # store
        self.ham = ham


def _random_density(rng, nvol, norb, paramagnetic):
    """Hermitian-closed real-space density rho[r, s, a, t, b] with
    rho(r) = conj(rho(-r))^T on the (s,a),(t,b) axes."""
    nd = 2 * norb
    g = rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd))
    if paramagnetic:
        blk = g[:, :norb, :norb]
        g = np.zeros_like(g)
        g[:, :norb, :norb] = blk
        g[:, norb:, norb:] = blk
    return g.reshape(nvol, 2, norb, 2, norb)


def _cases():
    rng = np.random.default_rng(7)
    shape = (2, 2, 1)
    nvol = 4
    one = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0}}
    kanamori = {
        "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
        "CoulombInter": {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7},
        "Hund": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
        "Exchange": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
        "PairHop": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3},
    }
    offsite = {
        "CoulombInter": {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                         ((0, 1, 0), (0, 0)): 0.2, ((0, -1, 0), (0, 0)): 0.2},
        "Hund": {((1, 0, 0), (0, 0)): 0.1, ((-1, 0, 0), (0, 0)): 0.1},
        "Ising": {((1, 0, 0), (0, 1)): 0.15, ((-1, 0, 0), (1, 0)): 0.15},
        "PairLift": {((1, 0, 0), (0, 1)): 0.05, ((-1, 0, 0), (1, 0)): 0.05},
        "Exchange": {((0, 1, 0), (0, 1)): 0.12, ((0, -1, 0), (1, 0)): 0.12},
        "PairHop": {((1, 0, 0), (0, 1)): 0.08, ((-1, 0, 0), (1, 0)): 0.08},
    }
    aggregate = {"Coulomb": {((0, 0, 0), (0, 0)): 2.0, ((1, 0, 0), (0, 1)): 0.3,
                             ((-1, 0, 0), (1, 0)): 0.3}}
    return [
        (one, 1, shape, _random_density(rng, nvol, 1, True)),
        (kanamori, 2, shape, _random_density(rng, nvol, 2, True)),
        (offsite, 2, shape, _random_density(rng, nvol, 2, True)),
        (offsite, 2, shape, _random_density(rng, nvol, 2, False)),
        (aggregate, 2, shape, _random_density(rng, nvol, 2, True)),
    ]


class TestKernelBitIdentity(unittest.TestCase):

    def test_tables_equal_legacy(self):
        from hwave.solver import hartree_fock as hf
        for param_ham, norb, shape, _ in _cases():
            with self.subTest(types=sorted(param_ham)):
                got = hf.build_interaction_tables(param_ham, norb, shape)
                leg = _LegacyUHFk(param_ham, norb, shape, True)
                leg._make_ham_inter()
                self.assertEqual(set(got.inter_table), set(leg.inter_table))
                for k in leg.inter_table:
                    if leg.inter_table[k] is None:
                        self.assertIsNone(got.inter_table[k], k)
                        continue
                    self.assertTrue(np.array_equal(leg.inter_table[k], got.inter_table[k]), k)
                    self.assertTrue(np.array_equal(leg.spin_table[k], got.spin_table[k]), k)

    def test_accumulate_equal_legacy_both_fock_settings(self):
        from hwave.solver import hartree_fock as hf
        rng = np.random.default_rng(3)
        for param_ham, norb, shape, rho in _cases():
            tabs = hf.build_interaction_tables(param_ham, norb, shape)
            nvol = int(np.prod(shape))
            nd = 2 * norb
            ham0 = rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd))
            for fock in (True, False):
                with self.subTest(types=sorted(param_ham), fock=fock):
                    leg = _LegacyUHFk(param_ham, norb, shape, fock)
                    leg._make_ham_inter()
                    leg.ham_trans = ham0.copy()
                    leg.Green = rho
                    leg._make_ham()
                    out = ham0.copy()
                    hf.accumulate_hf(out, rho, tabs.inter_table, tabs.spin_table, shape,
                                     include_fock=fock)
                    self.assertTrue(np.array_equal(out, leg.ham))

    def test_discarded_report_in_fixed_type_order(self):
        from hwave.solver import hartree_fock as hf
        ph = {"CoulombIntra": {((1, 0, 0), (0, 0)): 1.0, ((0, 0, 0), (0, 1)): 0.5,
                               ((0, 0, 0), (0, 0)): 2.0}}
        tabs = hf.build_interaction_tables(ph, 2, (2, 2, 1))
        self.assertEqual(tabs.discarded,
                         (("CoulombIntra", (1, 0, 0), (0, 0), 1.0),
                          ("CoulombIntra", (0, 0, 0), (0, 1), 0.5)))
        self.assertEqual(tabs.inter_table["CoulombIntra"][0, 0, 0, 0, 0], 2.0)
        self.assertEqual(hf.build_interaction_tables({}, 1, (2, 2, 1)).discarded, ())

    def test_ambiguous_coulomb_is_a_valueerror(self):
        from hwave.solver import hartree_fock as hf
        with self.assertRaises(ValueError):
            hf.build_interaction_tables({"Coulomb": {}, "CoulombIntra": {}}, 1, (2, 2, 1))

    def test_helpers(self):
        from hwave.solver import hartree_fock as hf
        a = np.array([[[1.0, 2.0 + 1j], [2.0 - 1j, 3.0]]])
        ok, err = hf.is_hermitian_batch(a)
        self.assertTrue(ok); self.assertEqual(err, 0.0)
        ok, err = hf.is_hermitian_batch(a + np.array([[[0, 1e-3], [0, 0]]]))
        self.assertFalse(ok); self.assertGreater(err, 1e-4)
        self.assertTrue(issubclass(hf.NonFiniteError, FloatingPointError))


class TestUHFkSolverBitIdentity(unittest.TestCase):
    """The production UHFk (normal mode) builds exactly the legacy
    Hamiltonian on real fixtures, iteration by iteration."""

    def _run(self, case, fock):
        import os, tempfile, tomli
        from hwave.qlmsio import read_input_k
        from hwave.solver.uhfk import UHFk
        cur = os.getcwd()
        os.chdir(os.path.join("tests", "uhfk", case))
        try:
            with open("input.toml", "rb") as f:
                params = tomli.load(f)
            params["file"]["input"].pop("initial", None)
            params["mode"]["param"]["flag_fock"] = fock
            params["mode"]["param"]["IterationMax"] = 3
            read_io = read_input_k.QLMSkInput(params["file"]["input"])
            ham_info = read_io.get_param("ham")
            green_info = read_io.get_param("green")
            solver = UHFk(ham_info, params["log"], params["mode"])
            seen = []
            orig = solver._make_ham

            def spy():
                orig()
                leg = _LegacyUHFk(solver.param_ham, solver.norb, solver.shape, solver.iflag_fock)
                leg._make_ham_inter()
                leg.ham_trans = solver.ham_trans
                leg.Green = solver.Green
                leg.block_info = solver.block_info
                leg._make_ham()
                seen.append(np.array_equal(solver.ham, leg.ham))
            solver._make_ham = spy
            with tempfile.TemporaryDirectory() as tmp:
                solver.solve(green_info, tmp)
            return seen
        finally:
            os.chdir(cur)

    def test_every_fixture_every_iteration(self):
        for case in ("CoulombIntra", "CoulombInter", "Ising", "Exchange", "PairHop", "PairLift", "2sz0"):
            for fock in (True, False):
                with self.subTest(case=case, fock=fock):
                    seen = self._run(case, fock)
                    self.assertTrue(seen and all(seen), (case, fock, seen))


if __name__ == "__main__":
    unittest.main()
