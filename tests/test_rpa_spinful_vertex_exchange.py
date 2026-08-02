#!/usr/bin/env python3

"""Spinful RPA antisymmetrized interaction vertex (issue #137).

The spinful (enable_spin_orbital) general solve resums a single vertex
tensor; ring-only content left every spin-flip pair slot of chiq at the
bare bubble. The vertex is now the antisymmetrized bare particle-hole
vertex Gamma = D + crossed(D)|on-site.

Ground truths used here, in decreasing independence:
  1. exact diagonalization of a genuinely spinful interacting lattice
     model (2 sites x 2 orbitals, dim 256): dchi/dV at V -> 0 must
     equal the RPA first order once the Hartree-Fock self-energy
     insertion (computed from the same free model) is removed;
  2. the vertex EXTRACTED from the solver itself,
     W = chi0^-1 (dchi/dU) chi0^-1, against the analytic D + X slots;
  3. the spin-conserving limit: density slots must match the spin-free
     (non-spin-orbital ring) solve, spin-flip slots the adjudicated
     ring+ladder transverse series, as the spin flip eps -> 0.

Tests must run from the repository root.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


T = 0.5
BETA = 1.0 / T
MU = 0.2
THOP = 0.7 + 0.3j
LSO = 0.35 + 0.15j


def _write_geom(d, nd):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % nd)
        for _ in range(nd):
            f.write("0.0 0.0 0.0\n")


def _run_rpa(d, inter, lx, nmat, so=True, scheme="general",
             calc_type="ring", exchange=None):
    param = {"T": T, "mu": MU, "CellShape": [lx, 1, 1],
             "SubShape": [1, 1, 1], "Nmat": nmat, "coeff_tail": 1.0}
    if exchange is not None:
        param["spinful_vertex_exchange"] = exchange
    info_mode = {"mode": "RPA", "param": param,
                 "enable_spin_orbital": so, "calc_scheme": scheme,
                 "calc_type": calc_type}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
    gi = io.get_param("green")
    out = os.path.join(d, "out")
    os.makedirs(out, exist_ok=True)
    solver.solve(gi, out)
    return solver, gi


class TestVertexExtraction(unittest.TestCase):
    """Extract W = chi0^-1 (dchi/dU) chi0^-1 from the solver on a
    1-orbital spinful chain with CoulombIntra U: it must be the
    antisymmetrized vertex D + X (in units of U: -1 on the density
    slots (uu),(dd)+(dd),(uu); +1 on the spin-flip pair diagonal
    (ud),(ud)+(du),(du)), q-independent. With
    spinful_vertex_exchange = false only D must remain (the pre-#137
    ring-only content)."""

    LX = 4
    NMAT = 256

    def _write(self, d, u):
        _write_geom(d, 2)
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n2\n3\n1 1 1\n")
            for i in (1, 2):
                f.write(" 1 0 0 %d %d %.12f %.12f\n"
                        % (i, i, THOP.real, THOP.imag))
                f.write("-1 0 0 %d %d %.12f %.12f\n"
                        % (i, i, np.conj(THOP).real, np.conj(THOP).imag))
            f.write(" 0 0 0 1 2 %.12f %.12f\n" % (LSO.real, LSO.imag))
            f.write(" 0 0 0 2 1 %.12f %.12f\n" % (LSO.real, -LSO.imag))
            f.write(" 0 0 0 1 1 0.10 0.0\n")
            f.write(" 0 0 0 2 2 -0.10 0.0\n")
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat"}
        if u != 0.0:
            with open(os.path.join(d, "coulombintra.dat"), "w") as f:
                f.write("hdr\n1\n1\n1\n")
                f.write(" 0 0 0 1 1 %.12f 0.0\n" % u)
            inter["CoulombIntra"] = "coulombintra.dat"
        return inter

    def _chiq(self, u, exchange):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = self._write(tmp.name, u)
        _, gi = _run_rpa(tmp.name, inter, self.LX, self.NMAT,
                         exchange=exchange)
        key = "chi0q" if u == 0.0 else "chiq"
        arr = np.asarray(gi[key])
        nf = arr.shape[0]
        return arr.reshape(nf, self.LX, 4, 4)[nf // 2]

    def _extract(self, exchange):
        U1, U2 = 0.02, 0.04
        chi0 = self._chiq(0.0, exchange)
        d1 = (self._chiq(U1, exchange) - chi0) / U1
        d2 = (self._chiq(U2, exchange) - chi0) / U2
        dchi = 2 * d1 - d2                       # Richardson
        W = np.zeros((self.LX, 4, 4), dtype=complex)
        for q in range(self.LX):
            inv = np.linalg.inv(chi0[q])
            W[q] = inv @ dchi[q] @ inv
        return W

    def test_antisymmetrized_vertex(self):
        D = np.zeros((4, 4))
        D[0, 3] = D[3, 0] = -1.0                 # (uu),(dd) direct
        X = np.zeros((4, 4))
        X[1, 1] = X[2, 2] = +1.0                 # (ud),(ud) exchange
        W = self._extract(True)
        for q in range(self.LX):
            self.assertLess(np.abs(W[q] - D - X).max(), 2e-3)

    def test_ring_only_optout(self):
        D = np.zeros((4, 4))
        D[0, 3] = D[3, 0] = -1.0
        W = self._extract(False)
        for q in range(self.LX):
            self.assertLess(np.abs(W[q] - D).max(), 2e-3)

    def test_flag_rejects_non_bool(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = self._write(tmp.name, 0.1)
        with self.assertRaises(ValueError):
            _run_rpa(tmp.name, inter, self.LX, 16, exchange="true")


class TestSpinConservingLimits(unittest.TestCase):
    """As the spin flip eps -> 0, for EVERY interaction type: the
    density slots of the spinful chiq must converge to the spin-free
    (non-spin-orbital ring) result, and the spin-flip slots to the
    adjudicated ring+ladder transverse series. PairLift is the
    documented exception on the transverse side: ring+ladder drops its
    transverse content ('outside both used blocks'), while the
    antisymmetrized spinful vertex keeps it -- exact diagonalization
    sides with the spinful result (TestFirstOrderAgainstED), so only a
    loose bound records the known ring+ladder gap here."""

    LX = 4
    NORB = 2
    NMAT = 256
    EPS = 1e-4
    T00, T11, T01 = 0.7 + 0.3j, 0.5 - 0.2j, 0.25 + 0.15j
    E0, E1 = 0.10, -0.05
    V = 0.3

    TYPES = {
        "CoulombIntra": [(0, 0)],
        "CoulombInter": [(0, 1), (1, 0)],
        "Hund": [(0, 1), (1, 0)],
        "Ising": [(0, 1), (1, 0)],
        "Exchange": [(0, 1), (1, 0)],
        "PairLift": [(0, 1), (1, 0)],
        "PairHop": [(0, 1), (1, 0)],
    }

    def _write(self, d, so, eps, tname):
        nd = 2 * self.NORB if so else self.NORB
        _write_geom(d, nd)
        hop = {(0, 0): self.T00, (1, 1): self.T11,
               (0, 1): self.T01, (1, 0): self.T01}
        def idx(o, s):
            return (2 * o + s if so else o) + 1
        smax = 2 if so else 1
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n%d\n3\n1 1 1\n" % nd)
            for (a, b), t in hop.items():
                for s in range(smax):
                    f.write(" 1 0 0 %d %d %.12f %.12f\n"
                            % (idx(a, s), idx(b, s), t.real, t.imag))
                    f.write("-1 0 0 %d %d %.12f %.12f\n"
                            % (idx(a, s), idx(b, s),
                               np.conj(hop[(b, a)]).real,
                               np.conj(hop[(b, a)]).imag))
            for o, e in ((0, self.E0), (1, self.E1)):
                for s in range(smax):
                    f.write(" 0 0 0 %d %d %.12f 0.0\n"
                            % (idx(o, s), idx(o, s), e))
            if so and eps != 0.0:
                z = eps * (1.0 + 0.5j)
                for o in range(self.NORB):
                    f.write(" 0 0 0 %d %d %.12g %.12g\n"
                            % (idx(o, 0), idx(o, 1), z.real, z.imag))
                    f.write(" 0 0 0 %d %d %.12g %.12g\n"
                            % (idx(o, 1), idx(o, 0), z.real, -z.imag))
        fn = tname.lower() + ".dat"
        ents = self.TYPES[tname]
        with open(os.path.join(d, fn), "w") as f:
            f.write("hdr\n%d\n%d\n%s\n" % (self.NORB, len(ents),
                                           " ".join("1" * len(ents))))
            for a, b in ents:
                f.write(" 0 0 0 %d %d %.12f 0.0\n" % (a + 1, b + 1, self.V))
        return {"path_to_input": d, "Geometry": "geom.dat",
                "Transfer": "transfer.dat", tname: fn}

    def _run(self, so, eps, tname, calc_type="ring"):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = self._write(tmp.name, so, eps, tname)
        return _run_rpa(tmp.name, inter, self.LX, self.NMAT, so=so,
                        calc_type=calc_type)

    def test_all_types(self):
        LX, NORB = self.LX, self.NORB
        ND = 2 * NORB
        for tname in self.TYPES:
            with self.subTest(type=tname):
                _, g0 = self._run(True, 0.0, tname)
                s1, g1 = self._run(True, self.EPS, tname)
                self.assertEqual(s1.spin_mode, "spinful")
                c0 = np.asarray(g0["chiq"])
                c1 = np.asarray(g1["chiq"])
                nf = c0.shape[0]
                r0 = c0.reshape(nf, LX, ND, ND, ND, ND)[nf // 2]
                r1 = c1.reshape(nf, LX, ND, ND, ND, ND)[nf // 2]
                mask = np.abs(r0) > 1e-12
                self.assertLess(np.abs((r1 - r0)[mask]).max(), 1e-7)
                _, gpm = self._run(False, 0.0, tname,
                                   calc_type="ring+ladder")
                pm_ref = np.asarray(gpm["chiq_pm"])
                nfp = pm_ref.shape[0]
                pm_ref = pm_ref.reshape(nfp, LX, NORB, NORB, NORB,
                                        NORB)[nfp // 2]
                up, dn = slice(0, NORB), slice(NORB, ND)
                pm1 = r1[:, up, dn, up, dn]
                dev = np.abs(pm1 - pm_ref).max()
                if tname == "PairLift":
                    self.assertLess(dev, 2e-2)
                else:
                    self.assertLess(dev, 1e-7)


class TestFirstOrderAgainstED(unittest.TestCase):
    """Exact diagonalization oracle on a genuinely spinful interacting
    model: 2 sites x 2 orbitals x 2 spins (dim 256), complex hoppings,
    orbital mixing and a LARGE spin flip. dchi/dV at V -> 0 minus the
    Hartree-Fock self-energy-insertion piece must equal the RPA first
    order. This pins the crossing construction for a same-orbital
    density type (CoulombIntra), a spin-exchange type (Exchange) and
    the type ring+ladder cannot represent (PairLift)."""

    LX = 2
    NORB = 2
    NSO = 4
    NMODE = 8
    DIM = 256
    NMAT = 512
    T00, T11, T01 = 0.7 + 0.3j, 0.5 - 0.2j, 0.25 + 0.15j
    E0, E1 = 0.10, -0.05

    @classmethod
    def _mode(cls, j, o, s):
        return cls.NSO * j + 2 * o + s

    @classmethod
    def setUpClass(cls):
        H1 = np.zeros((cls.NMODE, cls.NMODE), dtype=complex)
        hop = {(0, 0): cls.T00, (1, 1): cls.T11,
               (0, 1): cls.T01, (1, 0): cls.T01}
        for j in range(cls.LX):
            jp = (j + 1) % cls.LX
            for (a, b), t in hop.items():
                for s in range(2):
                    H1[cls._mode(jp, a, s), cls._mode(j, b, s)] += t
                    H1[cls._mode(j, a, s), cls._mode(jp, b, s)] += \
                        np.conj(hop[(b, a)])
            for o, e in ((0, cls.E0), (1, cls.E1)):
                for s in range(2):
                    H1[cls._mode(j, o, s), cls._mode(j, o, s)] += e
            for o in range(cls.NORB):
                H1[cls._mode(j, o, 0), cls._mode(j, o, 1)] += LSO
                H1[cls._mode(j, o, 1), cls._mode(j, o, 0)] += np.conj(LSO)
        cls.H1 = H1
        ops = []
        for p in range(cls.NMODE):
            c = np.zeros((cls.DIM, cls.DIM))
            for st in range(cls.DIM):
                if (st >> p) & 1:
                    sign = (-1) ** bin(st & ((1 << p) - 1)).count("1")
                    c[st ^ (1 << p), st] = sign
            ops.append(c.astype(complex))
        cls.C = ops
        cls.CD = [c.conj().T for c in ops]
        cls.hop = hop

    def _h_int(self, tname, v):
        C, CD, mode = self.C, self.CD, self._mode
        H = np.zeros((self.DIM, self.DIM), dtype=complex)
        for j in range(self.LX):
            if tname == "CoulombIntra":
                H += v * (CD[mode(j, 0, 0)] @ C[mode(j, 0, 0)]
                          @ CD[mode(j, 0, 1)] @ C[mode(j, 0, 1)])
            elif tname == "Exchange":
                for a, b in ((0, 1), (1, 0)):
                    H += v * (CD[mode(j, a, 0)] @ C[mode(j, b, 0)]
                              @ CD[mode(j, b, 1)] @ C[mode(j, a, 1)])
            elif tname == "PairLift":
                A = (CD[mode(j, 0, 0)] @ C[mode(j, 0, 1)]
                     @ CD[mode(j, 1, 0)] @ C[mode(j, 1, 1)])
                H += v * (A + A.conj().T)
            else:
                raise ValueError(tname)
        return H

    def _chi_ed(self, tname, v, h1=None):
        C, CD = self.C, self.CD
        if h1 is None:
            h1 = self.H1
            Hint = self._h_int(tname, v)
        else:
            Hint = 0.0
        with np.errstate(all="ignore"):
            H = np.zeros((self.DIM, self.DIM), dtype=complex)
            for p in range(self.NMODE):
                for q in range(self.NMODE):
                    if h1[p, q] != 0:
                        H += h1[p, q] * (CD[p] @ C[q])
            H = H + Hint
            N = sum(CD[p] @ C[p] for p in range(self.NMODE))
            E, V = np.linalg.eigh(H - MU * N)
        E = E - E.min()
        w = np.exp(-BETA * E)
        w /= w.sum()
        ND = self.NSO
        O = {}
        with np.errstate(all="ignore"):
            for qi in range(self.LX):
                for a in range(ND):
                    for c in range(ND):
                        op = np.zeros((self.DIM, self.DIM), dtype=complex)
                        for j in range(self.LX):
                            ph = np.exp(2j * np.pi * qi * j / self.LX)
                            op += ph * (CD[self.NSO * j + a]
                                        @ C[self.NSO * j + c])
                        O[(qi, a, c)] = V.conj().T @ op @ V
        dE = E[None, :] - E[:, None]
        dw = w[:, None] - w[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            K = np.where(np.abs(dE) > 1e-10, dw / dE, BETA * w[:, None])
        out = np.zeros((self.LX, ND, ND, ND, ND), dtype=complex)
        for qi in range(self.LX):
            qneg = (-qi) % self.LX
            for a in range(ND):
                for c in range(ND):
                    A = O[(qi, a, c)]
                    for b in range(ND):
                        for d in range(ND):
                            B = O[(qneg, d, b)]
                            val = (K * A * B.T).sum()
                            if qi == 0:
                                val -= BETA * (w * np.diag(A)).sum() \
                                    * (w * np.diag(B)).sum()
                            out[qi, a, c, b, d] = val
        return out / self.LX

    def _hf_h1(self, tname, v):
        """H1 + first-order Hartree-Fock self-energy of the interaction
        (Wick with the free density matrix)."""
        E, W = np.linalg.eigh(self.H1 - MU * np.eye(self.NMODE))
        f = 1.0 / (np.exp(np.clip(BETA * E, -500, 500)) + 1.0)
        rho = ((W * f) @ W.conj().T).T
        S = np.zeros((self.NMODE, self.NMODE), dtype=complex)
        mode = self._mode

        def add(p, q, r, s, v_):
            # E = v [rho_pq rho_rs + rho_ps (d_qr - rho_rq)];
            # H_MF = sum (dE/d rho_xy) c^+_x c_y
            S[p, q] += v_ * rho[r, s]
            S[r, s] += v_ * rho[p, q]
            S[p, s] += v_ * ((1.0 if q == r else 0.0) - rho[r, q])
            S[r, q] += -v_ * rho[p, s]

        for j in range(self.LX):
            if tname == "CoulombIntra":
                add(mode(j, 0, 0), mode(j, 0, 0),
                    mode(j, 0, 1), mode(j, 0, 1), v)
            elif tname == "Exchange":
                for a, b in ((0, 1), (1, 0)):
                    add(mode(j, a, 0), mode(j, b, 0),
                        mode(j, b, 1), mode(j, a, 1), v)
            elif tname == "PairLift":
                add(mode(j, 0, 0), mode(j, 0, 1),
                    mode(j, 1, 0), mode(j, 1, 1), v)
                add(mode(j, 1, 1), mode(j, 1, 0),
                    mode(j, 0, 1), mode(j, 0, 0), v)
        S = 0.5 * (S + S.conj().T)
        return self.H1 + S

    def _run_rpa(self, tname, v):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        _write_geom(d, self.NSO)
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n%d\n3\n1 1 1\n" % self.NSO)
            def idx(o, s):
                return 2 * o + s + 1
            for (a, b), t in self.hop.items():
                for s in range(2):
                    f.write(" 1 0 0 %d %d %.12f %.12f\n"
                            % (idx(a, s), idx(b, s), t.real, t.imag))
                    f.write("-1 0 0 %d %d %.12f %.12f\n"
                            % (idx(a, s), idx(b, s),
                               np.conj(self.hop[(b, a)]).real,
                               np.conj(self.hop[(b, a)]).imag))
            for o, e in ((0, self.E0), (1, self.E1)):
                for s in range(2):
                    f.write(" 0 0 0 %d %d %.12f 0.0\n"
                            % (idx(o, s), idx(o, s), e))
            for o in range(self.NORB):
                f.write(" 0 0 0 %d %d %.12f %.12f\n"
                        % (idx(o, 0), idx(o, 1), LSO.real, LSO.imag))
                f.write(" 0 0 0 %d %d %.12f %.12f\n"
                        % (idx(o, 1), idx(o, 0), LSO.real, -LSO.imag))
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat"}
        if v != 0.0:
            fn = tname.lower() + ".dat"
            with open(os.path.join(d, fn), "w") as f:
                if tname == "CoulombIntra":
                    f.write("hdr\n%d\n1\n1\n" % self.NORB)
                    f.write(" 0 0 0 1 1 %.12f 0.0\n" % v)
                else:
                    f.write("hdr\n%d\n2\n1 1\n" % self.NORB)
                    f.write(" 0 0 0 1 2 %.12f 0.0\n" % v)
                    f.write(" 0 0 0 2 1 %.12f 0.0\n" % v)
            inter[tname] = fn
        _, gi = _run_rpa(d, inter, self.LX, self.NMAT)
        key = "chi0q" if v == 0.0 else "chiq"
        arr = np.asarray(gi[key])
        nf = arr.shape[0]
        return arr.reshape(nf, self.LX, self.NSO, self.NSO, self.NSO,
                           self.NSO)[nf // 2]

    def _to_rpa_layout(self, x):
        """ED (interleaved 2*o+s) -> solver layout (spin-block combined
        index, conj, q-reversal); the mapping is pinned by the V = 0
        calibration assert in test_first_order."""
        i2sb = [None] * self.NSO
        for o in range(self.NORB):
            for s in range(2):
                i2sb[s * self.NORB + o] = 2 * o + s
        y = x[:, i2sb][:, :, i2sb][:, :, :, i2sb][:, :, :, :, i2sb]
        return y.conj()[[(-q) % self.LX for q in range(self.LX)]]

    def test_first_order(self):
        V1, V2 = 0.02, 0.04
        chi0_rpa = self._run_rpa("CoulombIntra", 0.0)
        chi0_ed = self._chi_ed("CoulombIntra", 0.0)
        cal = np.abs(self._to_rpa_layout(chi0_ed) - chi0_rpa).max()
        self.assertLess(cal, 1e-5)          # layout + free-model pin
        for tname in ("PairLift", "Exchange", "CoulombIntra"):
            with self.subTest(type=tname):
                d_ed = self._to_rpa_layout(
                    2 * (self._chi_ed(tname, V1) - chi0_ed) / V1
                    - (self._chi_ed(tname, V2) - chi0_ed) / V2)
                d_se = self._to_rpa_layout(
                    2 * (self._chi_ed(tname, 0.0,
                                      h1=self._hf_h1(tname, V1))
                         - chi0_ed) / V1
                    - (self._chi_ed(tname, 0.0,
                                    h1=self._hf_h1(tname, V2))
                       - chi0_ed) / V2)
                d_rpa = (2 * (self._run_rpa(tname, V1) - chi0_rpa) / V1
                         - (self._run_rpa(tname, V2) - chi0_rpa) / V2)
                dev = np.abs((d_ed - d_se) - d_rpa).max()
                scale = np.abs(d_ed - d_se).max()
                self.assertGreater(scale, 0.05)   # oracle is nontrivial
                self.assertLess(dev, 5e-4)


if __name__ == "__main__":
    unittest.main()
