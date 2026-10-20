"""The general FLEX ring is contracted in the conserving pair orientation
(issue #251).

A FLEX self-energy is Phi-derivable: its order-n ring term is the functional
derivative of ``Phi_n[G] = sum_q Tr (chi0[G] U)^n``. Using ONLY the solver's
own routines -- the bubble ``_calc_chi0q``, the effective interaction
``_calc_veff_general`` (default ``flex_second_order = "local"``) and the
self-energy kernel
``_calc_self_energy_general`` -- this is checked as the statement that the
ratio ``dPhi_n / Tr[Sigma_n dG]`` is the SAME for every direction ``dG`` (a
wrongly oriented kernel gives an O(1) spread). The Green function and the
directions are random complex arrays with no symmetry, so no orientation can
pass by coincidence. The vertices are the solver's own spin and charge
matrices of a Kanamori model with complex hopping (their Hermitian closure is
what lets the two G legs of the bubble combine into ONE kernel orientation;
a structureless random matrix has no conserving orientation at all).

The order-n term is isolated inside the production assembly by feeding it
truncated susceptibilities: ``chi = chi0`` gives the second-order kernel and
``chi = chi0 +- chi0 U chi0`` adds exactly the third-order ring term.

Before issue #251 the ring (third order and up) was contracted in the
transposed orientation of the exact second order (``second_order.
accumulate_batch``): the ratio spread was 0.5-3.9 at third order against
1e-10 now (measured on this fixture).
"""
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

NORB = 2
NMAT = 16


def _write_inputs(d):
    rng = np.random.default_rng(3)
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % NORB)
        f.write("0 0 0\n" * NORB)
    H0 = rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB))
    H0 = 0.3 * (H0 + H0.conj().T)
    T = 0.4 * (rng.normal(size=(NORB, NORB)) + 1j * rng.normal(size=(NORB, NORB))) + np.eye(NORB)
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n%d\n1\n 1\n" % NORB)
        for a in range(NORB):
            for b in range(NORB):
                for R, v in ((0, H0[a, b]), (1, T[a, b]), (-1, T[b, a].conjugate())):
                    f.write(" %d 0 0 %d %d %.12f %.12f\n" % (R, a + 1, b + 1, v.real, v.imag))
    files = {}
    offd = [(a, b) for a in range(NORB) for b in range(NORB) if a != b]
    spec = {"CoulombIntra": [(a, a, 1.0) for a in range(NORB)],
            "CoulombInter": [(a, b, 0.5) for a, b in offd],
            "Hund": [(a, b, 0.15) for a, b in offd],
            "Exchange": [(a, b, 0.15) for a, b in offd],
            "PairHop": [(a, b, 0.07) for a, b in offd]}
    for nm, rows in spec.items():
        fn = nm.lower() + ".dat"
        with open(os.path.join(d, fn), "w") as f:
            f.write("%s\n%d\n1\n 1\n" % (nm, NORB))
            for a, b, v in rows:
                f.write(" 0 0 0 %d %d %.12f 0.0\n" % (a + 1, b + 1, v))
        files[nm] = fn
    return files


class TestRingOrientation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod
        cls.d = tempfile.mkdtemp(prefix="flex_ring_orient_")
        files = _write_inputs(cls.d)
        inter = {"path_to_input": cls.d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        inter.update(files)
        io = read_input_k.QLMSkInput({"path_to_input": cls.d, "interaction": inter})
        par = {"T": 0.5, "filling": 0.4, "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": NMAT, "IterationMax": 0, "orbital_site": [0] * NORB}
        prev = logging.root.manager.disable
        logging.disable(logging.WARNING)
        try:
            cls.s = flex_mod.FLEX(io.get_param("ham"), {}, {
                "mode": "FLEX", "param": par, "enable_spin_orbital": False,
                "calc_scheme": "general"})
        finally:
            logging.disable(prev)
        cls.beta = 1.0 / cls.s.T
        cls.nvol = cls.s.lattice.nvol
        cls.ndx = NORB * NORB
        rng = np.random.default_rng(1)
        shape = (1, NMAT, cls.nvol, NORB, NORB)
        cls.G = 0.3 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))
        cls.dirs = [rng.normal(size=shape) + 1j * rng.normal(size=shape) for _ in range(4)]
        c0 = cls.s._calc_chi0q(cls.G, np.zeros_like(cls.G), cls.beta)[0]
        _, Us, Uc = cls.s._inflate_chi0q_and_ham_general(c0, cls.s.ham_info.ham_inter_q)
        cls.Us = np.asarray(Us).reshape(cls.nvol, cls.ndx, cls.ndx)
        cls.Uc = np.asarray(Uc).reshape(cls.nvol, cls.ndx, cls.ndx)
        cls.Z = np.zeros_like(cls.Us)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def _chi0(self, G):
        c = self.s._calc_chi0q(G, np.zeros_like(G), self.beta)[0]
        return c                                    # (nmat, nvol, no, no, no, no)

    def _phi(self, G, U, n):
        c = self._chi0(G).reshape(NMAT, self.nvol, self.ndx, self.ndx)
        m = c @ U[None]
        p = m.copy()
        for _ in range(n - 1):
            p = p @ m
        return np.trace(p, axis1=-2, axis2=-1).sum()

    def _ratio_spread(self, sigma, U, n):
        ratios = []
        for dG in self.dirs:
            h = 1e-6
            dphi = (self._phi(self.G + h * dG, U, n) - self._phi(self.G - h * dG, U, n)) / (2 * h)
            ratios.append(dphi / np.einsum("gfkmn,gfknm->", sigma, dG))
        r = np.array(ratios)
        return float(np.max(np.abs(r - r.mean())) / abs(r.mean()))

    def _veff(self, chi_s, chi_c, Us, Uc, kernel):
        factors = self.s._second_order_factors if kernel == "local" else None
        return self.s._calc_veff_general(self._chi0(self.G), chi_s, chi_c, Us, Uc,
                                         factors=factors, second_order=kernel)

    def _ring3(self, channel, kernel, U):
        """Sigma of the third-order ring term of one channel, isolated by
        differencing two production assemblies."""
        c0 = self._chi0(self.G)
        sh = c0.shape
        c2 = c0.reshape(NMAT, self.nvol, self.ndx, self.ndx)
        sign = 1.0 if channel == "spin" else -1.0     # chi_s = [1 - chi0 Us]^-1 chi0, chi_c with +
        third = (c2 + sign * (c2 @ U[None] @ c2)).reshape(sh)
        if channel == "spin":
            a = self._veff(third, c0, U, self.Z, kernel)
            b = self._veff(c0, c0, U, self.Z, kernel)
        else:
            a = self._veff(c0, third, self.Z, U, kernel)
            b = self._veff(c0, c0, self.Z, U, kernel)
        return self.s._calc_self_energy_general(self.G, a - b, self.beta)

    def test_third_order_ring_is_phi_derivable(self):
        # the default "local" kernel; the legacy "takimoto" kernel keeps the
        # earlier orientation on purpose (it reproduces earlier releases)
        for channel, U in (("spin", self.Us), ("charge", self.Uc)):
            with self.subTest(channel=channel):
                sigma = self._ring3(channel, "local", U)
                self.assertGreater(float(np.abs(sigma).max()), 1e-6)
                self.assertLess(self._ratio_spread(sigma, U, 3), 1e-7)

class TestBondRingOrientation(unittest.TestCase):
    """The bond-resolved channels (longitudinal_bond_channels) obey the same
    rule in their enlarged basis: one production map
    (``FLEX._phase_b_bond_map``: bond bubble, dressing, effective
    interaction, bond transport) with the vertices scaled by ``x``; the part
    odd in ``x`` is the ring's odd orders (every second-order term is even
    in ``x`` or does not scale at all), and its leading ``x^3`` coefficient
    must be the functional derivative of
    ``1.5 Tr (chibar S)^3 - 0.5 Tr (chibar C)^3``."""

    @classmethod
    def setUpClass(cls):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod
        cls.d = tempfile.mkdtemp(prefix="flex_bond_ring_orient_")
        files = _write_inputs(cls.d)
        with open(os.path.join(cls.d, "coulombinter.dat"), "a") as f:
            for a in range(NORB):
                for b in range(NORB):
                    v = 0.2 + 0.05 * a + 0.02 * b
                    f.write(" 1 0 0 %d %d %.12f 0.0\n -1 0 0 %d %d %.12f 0.0\n"
                            % (a + 1, b + 1, v, b + 1, a + 1, v))
        inter = {"path_to_input": cls.d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        inter.update(files)
        io = read_input_k.QLMSkInput({"path_to_input": cls.d, "interaction": inter})
        par = {"T": 0.5, "filling": 0.4, "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": NMAT, "IterationMax": 0, "orbital_site": [0] * NORB,
               "flex_hartree_fock": True, "longitudinal_bond_channels": True}
        prev = logging.root.manager.disable
        logging.disable(logging.WARNING)
        try:
            s = flex_mod.FLEX(io.get_param("ham"), {}, {
                "mode": "FLEX", "param": par, "enable_spin_orbital": False,
                "calc_scheme": "general"})
            s.solve(io.get_param("green"), cls.d)
            s._phase_b_prepare_vertices()
        finally:
            logging.disable(prev)
        cls.s = s
        cls.beta = 1.0 / s.T
        rng = np.random.default_rng(1)
        shape = (1, NMAT, s.lattice.nvol, NORB, NORB)
        cls.G = 0.3 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))
        cls.dirs = [rng.normal(size=shape) + 1j * rng.normal(size=shape) for _ in range(3)]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def _store(self):
        from hwave.solver import flex_bond
        s = self.s
        names, aliases = s._bond_store_layout()
        nd = NORB * NORB
        return flex_bond.BondBlockStore(NMAT, s.lattice.nvol, s._bond_view.n_channels * nd,
                                        nd, names, aliases)

    def _sigma(self, x):
        from hwave.solver import flex_bond
        s = self.s
        dev = flex_bond.BondDeviceContext.for_view(
            np, x * s._bond_S, x * s._bond_C, x * s._bond_S_on, x * s._bond_C_on,
            s._bond_view, NORB, green0_tail=np.zeros_like(self.G))
        with self._store() as store, dev:
            return np.asarray(s._phase_b_bond_map(store, dev, self.G, self.G, self.beta, 1)[3])

    def _phi3(self, G):
        from hwave.solver import flex_bond
        s = self.s
        with self._store() as store:
            flex_bond.assemble_bubble(store, G, np.zeros_like(G), self.beta, s._bond_view,
                                      tuple(int(v) for v in s.lattice.shape), 1)
            cb = np.array(store.get_freq_batch("chibar", 0, NMAT))
        tot = 0.0
        for U, w in ((s._bond_S, 1.5), (s._bond_C, -0.5)):
            m = cb @ np.asarray(U)[None]
            tot += w * np.trace(m @ m @ m, axis1=-2, axis2=-1).sum()
        return tot

    def test_third_order_bond_ring_is_phi_derivable(self):
        h = 0.02
        odd = lambda x: 0.5 * (self._sigma(x) - self._sigma(-x))
        sigma3 = (32.0 * odd(h) - odd(2 * h)) / (24.0 * h ** 3)   # x^5 removed
        r = []
        for dG in self.dirs:
            e = 1e-6
            dphi = (self._phi3(self.G + e * dG) - self._phi3(self.G - e * dG)) / (2 * e)
            r.append(dphi / np.einsum("gfkmn,gfknm->", sigma3, dG))
        r = np.array(r)
        self.assertLess(float(np.max(np.abs(r - r.mean())) / abs(r.mean())), 1e-5)


class TestBruteForceReferenceOrientation(unittest.TestCase):
    """The literal-loop MYO reference of tests/flex_bruteforce_ref.py obeys
    the same rule: its ring is conserving only with both pairs of the
    matrix product transposed into Eq. (3)'s ``V_{mu m, nu n}`` slots."""

    def test_reference_pair_is_phi_derivable_only_transposed(self):
        from tests.flex_bruteforce_ref import chi0_bruteforce, sigma_bruteforce
        no, nk, nm, T = 2, 3, 2, 0.5
        rng = np.random.default_rng(0)
        G = rng.normal(size=(no, no, nk, nm)) + 1j * rng.normal(size=(no, no, nk, nm))
        ndx = no * no
        u = rng.normal(size=(ndx, ndx)) + 1j * rng.normal(size=(ndx, ndx))
        u4 = (u + u.conj().T).reshape(no, no, no, no)
        # the Hermitian closure of a physical vertex: U[(a c),(b d)] = conj U[(c a),(d b)]
        U = (0.5 * (u4 + u4.transpose(1, 0, 3, 2).conj())).reshape(ndx, ndx)

        def slices(G):
            c = chi0_bruteforce(G, T, nk)
            return [c[:, :, :, :, q, iv].reshape(ndx, ndx) for q in range(nk) for iv in range(nm)]

        def phi(G, n):
            return sum(np.trace(np.linalg.matrix_power(m @ U, n)) for m in slices(G))

        dirs = [rng.normal(size=G.shape) + 1j * rng.normal(size=G.shape) for _ in range(3)]
        for n in (2, 3):
            spreads = {}
            for transposed in (False, True):
                ms = slices(G)
                V = np.zeros((no, no, no, no, nk, nm), complex)
                for i, (q, iv) in enumerate((q, iv) for q in range(nk) for iv in range(nm)):
                    v4 = (U @ np.linalg.matrix_power(ms[i] @ U, n - 1)).reshape(no, no, no, no)
                    V[:, :, :, :, q, iv] = v4.transpose(1, 0, 3, 2) if transposed else v4
                S = sigma_bruteforce(G, V, T, nk)
                r = []
                for dG in dirs:
                    h = 1e-6
                    dphi = (phi(G + h * dG, n) - phi(G - h * dG, n)) / (2 * h)
                    r.append(dphi / np.einsum("mnkw,nmkw->", S, dG))
                r = np.array(r)
                spreads[transposed] = float(np.max(np.abs(r - r.mean())) / abs(r.mean()))
            with self.subTest(order=n):
                self.assertLess(spreads[True], 1e-7)
                self.assertGreater(spreads[False], 1e-2)   # the untransposed wiring is not


if __name__ == "__main__":
    unittest.main()
