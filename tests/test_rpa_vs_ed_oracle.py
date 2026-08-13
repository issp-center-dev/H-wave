#!/usr/bin/env python3

"""Exact diagonalization as the ground truth for the RPA solver.

Coverage: the spin-free general scheme for CoulombIntra, CoulombInter,
Exchange and PairHop (real and complex), plus the reduced scheme for
CoulombIntra. Hund and Ising, the spin-diagonal and spinful paths and
the GPU backend are adjudicated elsewhere
(tests/test_rpa_spinful_vertex_exchange.py carries the spinful
exact-diagonalization oracle) or, for the interaction-tensor
conventions, by the fixed-point pin at the bottom of this module.

The susceptibility conventions of this code base (which orbital pair
sits in which slot, which index of a pair carries the creation
operator, the momentum sign) form an exact-symmetry family that makes
"scan the layouts until something matches" adjudications ambiguous: at
zero coupling several members agree, and they split only once the
interaction is switched on. This module therefore fixes the frame ONCE,
by DERIVATION from the documented equations, and then compares each
mode's first-order response against exact diagonalization of the
Hamiltonian the input-file specification defines.

The frame
---------
The documented bare susceptibility

    X0^{aa',bb'}(q) = -(T/N) sum_k G^{ab}(k+q) G^{b'a'}(k)

is the correlator of the bilinears

    O1 = c^+_{a'} c_{a}  at -q      and     O2 = c^+_{b} c_{b'}  at +q,

so a Lehmann evaluation storing ed[a, c, b, d] = <(c^+_a c_c);
(c^+_d c_b)> maps onto the solver's slots by swapping the two indices
of each pair. That single derived map (never a scan) is
``_to_solver_slots`` below; it is verified against the solver's own
chi0q at zero coupling, where it must hold to round-off.

The main fixture uses two sites, where q and -q coincide and the
momentum assignment above is therefore untestable; a separate
three-site calibration (``TestMomentumFrame``) pins it where the two
differ.

What is compared
----------------
The FIRST-ORDER response in the coupling. Exact diagonalization
contains both the vertex correction and the Hartree-Fock self-energy
insertion; the ring resummation contains only the former, so the
insertion piece (computed from the same free model) is subtracted from
the exact side. What remains differs from the ring prediction only by
the exchange wiring, which the ring form omits by construction -- so
the comparison is restricted to the slots the ring can populate. The
solver is additionally required not to populate any slot the exact
response leaves empty, and the number of jointly populated slots is
pinned per declaration type; the exchange-only slots themselves are
recorded by those counts rather than by an explicit mask.

Tests must be run from the repository root.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


LX, NORB, NS = 2, 2, 2
ND = NORB * NS                      # generalized orbital (spin-block)
NMODE = LX * NORB * NS
DIM = 1 << NMODE
T, MU = 0.5, 0.2
BETA = 1.0 / T
V1, V2 = 0.02, 0.04                 # Richardson pair for d/dV at V -> 0
HOP = {(0, 0): 0.7 + 0.3j, (1, 1): 0.5 - 0.2j,
       (0, 1): 0.25 + 0.15j, (1, 0): 0.25 + 0.15j}
E_ONSITE = (0.10, -0.05)
PHASE = np.exp(0.4j)                # complex-declaration phase


def _mode(j, o, s):
    """ED single-particle mode (site-major, interleaved 2*orb+spin)."""
    return 4 * j + 2 * o + s


def _sb(o, s):
    """Solver generalized-orbital index (spin-block, spin*norb+orb)."""
    return s * NORB + o


def _build_h1():
    h = np.zeros((NMODE, NMODE), dtype=complex)
    for j in range(LX):
        jp = (j + 1) % LX
        for (a, b), t in HOP.items():
            for s in range(2):
                h[_mode(jp, a, s), _mode(j, b, s)] += t
                h[_mode(j, a, s), _mode(jp, b, s)] += np.conj(HOP[(b, a)])
        for o, e in enumerate(E_ONSITE):
            for s in range(2):
                h[_mode(j, o, s), _mode(j, o, s)] += e
    return h


H1 = _build_h1()


def _annihilators():
    ops = []
    for p in range(NMODE):
        c = np.zeros((DIM, DIM))
        for st in range(DIM):
            if (st >> p) & 1:
                sign = (-1) ** bin(st & ((1 << p) - 1)).count("1")
                c[st ^ (1 << p), st] = sign
        ops.append(c.astype(complex))
    return ops


C = _annihilators()
CD = [c.conj().T for c in C]


# --------------------------------------------------------------------
# interaction definitions: the operator content the input-file
# specification assigns to each declaration, written out explicitly
# --------------------------------------------------------------------

def _n(j, o, s):
    return CD[_mode(j, o, s)] @ C[_mode(j, o, s)]


def _h_int(kind, v, phase=1.0):
    """Hamiltonian of one declaration type at coupling v.

    The operators follow the input-file specification (and the
    real-space UHF InterAll transforms): CoulombIntra U n_up n_dn;
    CoulombInter V n_a n_b; Exchange J S^+_a S^-_b; PairLift
    J c^+_{a up} c_{a dn} c^+_{b up} c_{b dn}; PairHop
    P c^+_{a up} c_{b up} c^+_{a dn} c_{b dn}. Complex declarations are
    Hermitian-closed, the partner carrying the conjugate.
    """
    H = np.zeros((DIM, DIM), dtype=complex)
    coef = v * phase
    for j in range(LX):
        if kind == "CoulombIntra":
            H += v * (_n(j, 0, 0) @ _n(j, 0, 1))
        elif kind == "CoulombInter":
            na = _n(j, 0, 0) + _n(j, 0, 1)
            nb = _n(j, 1, 0) + _n(j, 1, 1)
            H += v * (na @ nb)
        elif kind == "Exchange":
            for a, b in ((0, 1), (1, 0)):
                H += v * (CD[_mode(j, a, 0)] @ C[_mode(j, b, 0)]
                          @ CD[_mode(j, b, 1)] @ C[_mode(j, a, 1)])
        elif kind == "PairLift":
            A = (CD[_mode(j, 0, 0)] @ C[_mode(j, 0, 1)]
                 @ CD[_mode(j, 1, 0)] @ C[_mode(j, 1, 1)])
            H += v * (A + A.conj().T)
        elif kind in ("PairHop", "PairHopComplex"):
            A = (CD[_mode(j, 0, 0)] @ C[_mode(j, 1, 0)]
                 @ CD[_mode(j, 0, 1)] @ C[_mode(j, 1, 1)])
            H += coef * A + np.conj(coef) * A.conj().T
        else:
            raise ValueError(kind)
    return H


def _hf_h1(kind, v, phase=1.0):
    """H1 plus the first-order Hartree-Fock self-energy of the same
    interaction, by Wick contraction with the free density matrix."""
    ev, W = np.linalg.eigh(H1 - MU * np.eye(NMODE))
    f = 1.0 / (np.exp(np.clip(BETA * ev, -500, 500)) + 1.0)
    rho = ((W * f) @ W.conj().T).T          # rho[p, q] = <c^+_p c_q>
    S = np.zeros((NMODE, NMODE), dtype=complex)

    def add(p, q, r, s, c_):
        # E = c [rho_pq rho_rs + rho_ps (delta_qr - rho_rq)];
        # H_MF = sum (dE / d rho_xy) c^+_x c_y
        S[p, q] += c_ * rho[r, s]
        S[r, s] += c_ * rho[p, q]
        S[p, s] += c_ * ((1.0 if q == r else 0.0) - rho[r, q])
        S[r, q] += -c_ * rho[p, s]

    coef = v * phase
    for j in range(LX):
        if kind == "CoulombIntra":
            add(_mode(j, 0, 0), _mode(j, 0, 0),
                _mode(j, 0, 1), _mode(j, 0, 1), v)
        elif kind == "CoulombInter":
            for s1 in range(2):
                for s2 in range(2):
                    add(_mode(j, 0, s1), _mode(j, 0, s1),
                        _mode(j, 1, s2), _mode(j, 1, s2), v)
        elif kind == "Exchange":
            for a, b in ((0, 1), (1, 0)):
                add(_mode(j, a, 0), _mode(j, b, 0),
                    _mode(j, b, 1), _mode(j, a, 1), v)
        elif kind == "PairLift":
            add(_mode(j, 0, 0), _mode(j, 0, 1),
                _mode(j, 1, 0), _mode(j, 1, 1), v)
            add(_mode(j, 1, 1), _mode(j, 1, 0),
                _mode(j, 0, 1), _mode(j, 0, 0), v)
        elif kind in ("PairHop", "PairHopComplex"):
            add(_mode(j, 0, 0), _mode(j, 1, 0),
                _mode(j, 0, 1), _mode(j, 1, 1), coef)
            add(_mode(j, 1, 0), _mode(j, 0, 0),
                _mode(j, 1, 1), _mode(j, 0, 1), np.conj(coef))
        else:
            raise ValueError(kind)
    return H1 + 0.5 * (S + S.conj().T)


def _chi_ed(hint=None, h1=None):
    """Static chi[q, a, c, b, d] = <(c^+_a c_c)(q); (c^+_d c_b)(-q)>,
    connected, by the Lehmann representation."""
    H = np.zeros((DIM, DIM), dtype=complex)
    h1m = H1 if h1 is None else h1
    for p in range(NMODE):
        for q in range(NMODE):
            if h1m[p, q] != 0:
                H += h1m[p, q] * (CD[p] @ C[q])
    if hint is not None:
        H = H + hint
    N = sum(CD[p] @ C[p] for p in range(NMODE))
    with np.errstate(all="ignore"):
        ev, V = np.linalg.eigh(H - MU * N)
    ev = ev - ev.min()
    w = np.exp(-BETA * ev)
    w /= w.sum()
    O = {}
    for qi in range(LX):
        for a in range(ND):
            for c in range(ND):
                oa, sa = a % NORB, a // NORB
                oc, sc = c % NORB, c // NORB
                op = np.zeros((DIM, DIM), dtype=complex)
                for j in range(LX):
                    ph = np.exp(2j * np.pi * qi * j / LX)
                    op += ph * (CD[_mode(j, oa, sa)] @ C[_mode(j, oc, sc)])
                O[(qi, a, c)] = V.conj().T @ op @ V
    dE = ev[None, :] - ev[:, None]
    dw = w[:, None] - w[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(dE) > 1e-10, dw / dE, BETA * w[:, None])
    out = np.zeros((LX, ND, ND, ND, ND), dtype=complex)
    for qi in range(LX):
        qn = (-qi) % LX
        for a in range(ND):
            for c in range(ND):
                A = O[(qi, a, c)]
                for b in range(ND):
                    for d in range(ND):
                        B = O[(qn, d, b)]
                        val = (K * A * B.T).sum()
                        if qi == 0:
                            val -= BETA * (w * np.diag(A)).sum() \
                                * (w * np.diag(B)).sum()
                        out[qi, a, c, b, d] = val
    return out / LX


def _to_solver_slots(x):
    """DERIVED map from the Lehmann layout to the solver's slots.

    The documented X0^{aa',bb'} is the correlator of c^+_{a'} c_{a} with
    c^+_{b} c_{b'}, while ``_chi_ed`` stores <(c^+_a c_c); (c^+_d c_b)>:
    the two differ by swapping the indices inside each pair. No scan.
    """
    return np.transpose(x, (0, 2, 1, 4, 3))


def _ed_first_order(kind, phase=1.0):
    """Vertex part of d(chi)/dV at V -> 0, in the solver's slots: the
    exact response minus the Hartree-Fock insertion piece."""
    base = _chi_ed()
    d_full = (2 * (_chi_ed(_h_int(kind, V1, phase)) - base) / V1
              - (_chi_ed(_h_int(kind, V2, phase)) - base) / V2)
    d_ins = (2 * (_chi_ed(h1=_hf_h1(kind, V1, phase)) - base) / V1
             - (_chi_ed(h1=_hf_h1(kind, V2, phase)) - base) / V2)
    return _to_solver_slots(d_full - d_ins), _to_solver_slots(base)


# --------------------------------------------------------------------
# solver side
# --------------------------------------------------------------------

_FILES = {"CoulombIntra": [(0, 0)], "CoulombInter": [(0, 1), (1, 0)],
          "Exchange": [(0, 1), (1, 0)], "PairLift": [(0, 1), (1, 0)],
          "PairHop": [(0, 1), (1, 0)],
          "PairHopComplex": [(0, 1), (1, 0)]}


def _write_inputs(d, kind, v, phase=1.0):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % NORB)
        for _ in range(NORB):
            f.write("0.0 0.0 0.0\n")
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n3\n1 1 1\n" % NORB)
        for (a, b), t in HOP.items():
            f.write(" 1 0 0 %d %d %.12f %.12f\n"
                    % (a + 1, b + 1, t.real, t.imag))
            f.write("-1 0 0 %d %d %.12f %.12f\n"
                    % (a + 1, b + 1, np.conj(HOP[(b, a)]).real,
                       np.conj(HOP[(b, a)]).imag))
        for o, e in enumerate(E_ONSITE):
            f.write(" 0 0 0 %d %d %.12f 0.0\n" % (o + 1, o + 1, e))
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    ents = _FILES[kind]
    name = "PairHop" if kind == "PairHopComplex" else kind
    fn = name.lower() + ".dat"
    with open(os.path.join(d, fn), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n" % (NORB, len(ents),
                                       " ".join("1" * len(ents))))
        for i, (a, b) in enumerate(ents):
            val = v * (phase if i == 0 else np.conj(phase))
            f.write(" 0 0 0 %d %d %.12f %.12f\n"
                    % (a + 1, b + 1, val.real, val.imag))
    inter[name] = fn
    return inter


def _run_rpa(kind, v, scheme="general", phase=1.0, nmat=512):
    tmp = tempfile.TemporaryDirectory()
    try:
        inter = _write_inputs(tmp.name, kind, v, phase)
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": nmat,
                               "coeff_tail": 1.0},
                     "calc_scheme": scheme}
        io = read_input_k.QLMSkInput(
            {"path_to_input": tmp.name, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        arr = np.asarray(gi["chiq"])
    finally:
        tmp.cleanup()
    nf = arr.shape[0]
    if scheme == "reduced":
        return arr.reshape(nf, LX, ND, ND)[nf // 2]
    return arr.reshape(nf, LX, ND, ND, ND, ND)[nf // 2]


def _solver_first_order(kind, scheme, phase=1.0):
    c0 = _run_rpa(kind, 1e-9, scheme, phase)
    c1 = _run_rpa(kind, V1, scheme, phase)
    c2 = _run_rpa(kind, V2, scheme, phase)
    return 2 * (c1 - c0) / V1 - (c2 - c0) / V2


class TestFrameCalibration(unittest.TestCase):
    """The derived map must reproduce the solver's own bare bubble."""

    def test_chi0_matches_exact_lindhard(self):
        _, base = _ed_first_order("CoulombIntra")
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = _write_inputs(tmp.name, "CoulombIntra", 1e-9)
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 4096,
                               "coeff_tail": 1.0},
                     "calc_scheme": "general"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": tmp.name, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        chi0 = np.asarray(gi["chi0q"])
        nf = chi0.shape[0]
        # spin-free bubble is stored in orbital space; the exact one is
        # spin-diagonal in the pairing implied by the documented formula
        chi0 = chi0.reshape(nf, LX, NORB, NORB, NORB, NORB)[nf // 2]
        ref = base[:, :NORB, :NORB, :NORB, :NORB]
        self.assertLess(np.abs(chi0 - ref).max(), 1e-6)


class TestMomentumFrame(unittest.TestCase):
    """Pin the momentum assignment where q and -q actually differ.

    The main fixture has two sites, so every grid momentum satisfies
    q = -q. Even on a larger lattice a single band cannot settle the
    question: the static Lindhard kernel is symmetric under exchanging
    its two states, which makes chi0(q) even in q for one orbital. The
    q-odd content lives in the orbital off-diagonal slots, so this
    calibration uses a three-site ring with TWO orbitals and complex
    inter-orbital hopping, and compares the solver's bare bubble with
    an analytic Lindhard evaluation written with the documented
    assignment; the reversed assignment must be excluded."""

    L3 = 3
    HOP = {(0, 0): 0.7 + 0.3j, (1, 1): 0.5 - 0.2j,
           (0, 1): 0.25 + 0.15j, (1, 0): 0.25 + 0.15j}
    EPS = (0.10, -0.05)

    def _hk(self, k):
        h = np.zeros((2, 2), dtype=complex)
        ph = np.exp(2j * np.pi * k / self.L3)
        for (a, b), t in self.HOP.items():
            h[a, b] += t * ph + np.conj(self.HOP[(b, a)]) * np.conj(ph)
        for o, e in enumerate(self.EPS):
            h[o, o] += e
        return h

    def _lindhard(self, reversed_frame):
        ev = np.zeros((self.L3, 2))
        U = np.zeros((self.L3, 2, 2), dtype=complex)
        for k in range(self.L3):
            w, v = np.linalg.eigh(self._hk(k))
            ev[k], U[k] = w - MU, v

        def f(e):
            return 1.0 / (np.exp(np.clip(BETA * e, -500, 500)) + 1.0)

        out = np.zeros((self.L3, 2, 2, 2, 2), dtype=complex)
        for q in range(self.L3):
            for k in range(self.L3):
                kq = (k - q) % self.L3 if reversed_frame \
                    else (k + q) % self.L3
                for m in range(2):
                    for n in range(2):
                        de = ev[kq][m] - ev[k][n]
                        if abs(de) < 1e-12:
                            w0 = -BETA * f(ev[k][n]) * (1 - f(ev[k][n]))
                        else:
                            w0 = (f(ev[kq][m]) - f(ev[k][n])) / de
                        out[q] += (-w0 / self.L3) * np.einsum(
                            "a,b,B,A->aAbB", U[kq][:, m], U[kq][:, m].conj(),
                            U[k][:, n], U[k][:, n].conj())
        return out

    def test_documented_assignment(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n2\n"
                    "0.0 0.0 0.0\n0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n2\n3\n1 1 1\n")
            for (a, b), t in self.HOP.items():
                f.write(" 1 0 0 %d %d %.12f %.12f\n"
                        % (a + 1, b + 1, t.real, t.imag))
                f.write("-1 0 0 %d %d %.12f %.12f\n"
                        % (a + 1, b + 1, np.conj(self.HOP[(b, a)]).real,
                           np.conj(self.HOP[(b, a)]).imag))
            for o, e in enumerate(self.EPS):
                f.write(" 0 0 0 %d %d %.12f 0.0\n" % (o + 1, o + 1, e))
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat"}
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU,
                               "CellShape": [self.L3, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 2048,
                               "coeff_tail": 1.0},
                     "calc_scheme": "general"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        chi0 = np.asarray(gi["chi0q"])
        nf = chi0.shape[0]
        got = chi0.reshape(nf, self.L3, 2, 2, 2, 2)[nf // 2]
        ref = self._lindhard(False)
        rev = self._lindhard(True)
        # the two assignments must be distinguishable at all
        self.assertGreater(np.abs(ref - rev).max(), 1e-3)
        self.assertLess(np.abs(got - ref).max(), 1e-4)
        self.assertGreater(np.abs(got - rev).max(), 1e-3)


class TestRPAGeneralVsED(unittest.TestCase):
    """Density and pair-crossing declarations, general scheme, against
    exact diagonalization: the solver must reproduce the exact
    first-order vertex response on the slots the ring form can
    populate, and must not populate any slot the exact response leaves
    empty. The number of jointly populated slots is pinned per type so
    that a regression which silently DROPS ring-representable slots
    cannot pass by shrinking the comparison mask."""

    CASES = [("CoulombIntra", 1.0), ("CoulombInter", 1.0),
             ("Exchange", 1.0), ("PairHop", 1.0),
             ("PairHopComplex", PHASE)]

    # relative bound: the solver evaluates its bubble on a finite
    # Matsubara grid while the exact side sums analytically, which
    # leaves a truncation difference of order 1e-3 relative at the
    # Nmat used here (measured 1.0e-4 on a 0.11 response); the
    # structural errors this module exists to catch are O(1) relative.
    RTOL = 2e-3

    # jointly populated slot counts, measured; a drop means the solver
    # stopped representing part of the channel
    SHARED = {"CoulombIntra": 64, "CoulombInter": 128, "Exchange": 64,
              "PairHop": 64, "PairHopComplex": 64}

    def test_first_order(self):
        for kind, phase in self.CASES:
            with self.subTest(type=kind):
                ed, _ = _ed_first_order(kind, phase)
                sol = _solver_first_order(kind, "general", phase)
                shared = (np.abs(ed) > 1e-6) & (np.abs(sol) > 1e-6)
                self.assertEqual(int(shared.sum()), self.SHARED[kind])
                scale = np.abs(ed).max()
                self.assertGreater(scale, 1e-3)
                self.assertLess(np.abs((ed - sol)[shared]).max(),
                                self.RTOL * scale)
                # the solver must not populate slots the exact response
                # leaves empty
                only_sol = (np.abs(sol) > 1e-6) & (np.abs(ed) <= 1e-6)
                self.assertLess(np.abs(sol[only_sol]).max(initial=0.0),
                                self.RTOL * scale)

    def test_pairlift_content_is_outside_the_longitudinal_ring(self):
        """PairLift couples spin-FLIP bilinears, and the spin-free
        longitudinal path carries only spin-conserving pair slots (the
        bare bubble is built in orbital space and inflated
        spin-diagonally). The transverse (+-) assembly does not pick
        this channel up either -- PairLift falls outside the two spin
        blocks it uses, so its ladder vertex vanishes identically -- so
        NO current spin-free path represents it; the spinful
        antisymmetrized vertex does (issue #137). The solver is
        therefore silent here, and the whole exact response must sit in
        spin-flip slots: a recorded channel limitation, so that a
        future change which starts populating them is noticed."""
        ed, _ = _ed_first_order("PairLift")
        sol = _solver_first_order("PairLift", "general")
        self.assertLess(np.abs(sol).max(), 1e-6)
        self.assertGreater(np.abs(ed).max(), 1e-3)
        flip = np.zeros((ND, ND, ND, ND), dtype=bool)
        for a in range(ND):
            for ap in range(ND):
                for b in range(ND):
                    for bp in range(ND):
                        if (a // NORB != ap // NORB) or \
                                (b // NORB != bp // NORB):
                            flip[a, ap, b, bp] = True
        self.assertLess(np.abs(ed[:, ~flip]).max(), 1e-6)


class TestReducedMatchesGeneralDiagonal(unittest.TestCase):
    """For a purely density-diagonal interaction the reduced scheme is
    the density projection of the general one, so it inherits the same
    exact-diagonalization agreement.

    This does NOT extend to interactions carrying pair-crossing content
    (an on-site CoulombInter, through its Fierz slots): the reduced
    scheme resums inside the density subspace only, dropping the
    cross-pair intermediates the general resummation keeps, and the two
    then differ already at first order (measured 8.9e-3 on a 0.11
    response). That is the documented size reduction of the scheme, not
    a defect, so it is deliberately not asserted here."""

    def test_density_projection(self):
        for kind in ("CoulombIntra",):
            with self.subTest(type=kind):
                gen = _solver_first_order(kind, "general")
                red = _solver_first_order(kind, "reduced")
                proj = np.einsum("qaabb->qab", gen)
                np.testing.assert_allclose(red, proj, atol=1e-10)
                ed, _ = _ed_first_order(kind)
                ed_proj = np.einsum("qaabb->qab", ed)
                self.assertLess(np.abs(red - ed_proj).max(),
                                2e-3 * np.abs(ed_proj).max())


class TestComplexPhaseIsNotConjugated(unittest.TestCase):
    """Direct regression for the bubble-pair convention (issue #139):
    a complex PairHop must produce the susceptibility of the DECLARED
    Hamiltonian, not of its complex conjugate. Reading the interaction
    without the intra-pair transpose returns the conjugate model, which
    this test separates by more than three orders of magnitude."""

    def test_declared_not_conjugate(self):
        sol = _solver_first_order("PairHopComplex", "general", PHASE)
        ed_decl, _ = _ed_first_order("PairHopComplex", PHASE)
        ed_conj, _ = _ed_first_order("PairHopComplex", np.conj(PHASE))
        shared = (np.abs(ed_decl) > 1e-6) & (np.abs(sol) > 1e-6)
        dev_decl = np.abs((ed_decl - sol)[shared]).max()
        shared_c = (np.abs(ed_conj) > 1e-6) & (np.abs(sol) > 1e-6)
        dev_conj = np.abs((ed_conj - sol)[shared_c]).max()
        self.assertLess(dev_decl, 2e-3 * np.abs(ed_decl).max())
        self.assertGreater(dev_conj, 1e-2)
        # the two models must genuinely differ, or the test is vacuous
        self.assertGreater(np.abs(ed_decl - ed_conj).max(), 1e-2)


class TestConversionIsAFixedPointForRealDeclarations(unittest.TestCase):
    """Backward-compatibility pin, asserted directly on the tensors
    rather than inferred from unchanged fixtures: for every declaration
    type with REAL Hermitian-closed coefficients the conversion is the
    identity, so no previously computed result can move. Only a complex
    pair-crossing declaration is affected, and there it must conjugate
    the pair-crossing slots."""

    REAL_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising",
                  "Exchange", "PairLift", "PairHop")

    def _vertices(self, kind, phase=1.0, offsite=False):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % NORB)
            for _ in range(NORB):
                f.write("0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n%d\n3\n1 1 1\n" % NORB)
            for (a, b), t in HOP.items():
                f.write(" 1 0 0 %d %d %.12f %.12f\n"
                        % (a + 1, b + 1, t.real, t.imag))
                f.write("-1 0 0 %d %d %.12f %.12f\n"
                        % (a + 1, b + 1, np.conj(HOP[(b, a)]).real,
                           np.conj(HOP[(b, a)]).imag))
            for o, e in enumerate(E_ONSITE):
                f.write(" 0 0 0 %d %d %.12f 0.0\n" % (o + 1, o + 1, e))
        fn = kind.lower() + ".dat"
        ents = ([(1, 0, 1), (-1, 1, 0)] if offsite
                else ([(0, 0, 0)] if kind == "CoulombIntra"
                      else [(0, 0, 1), (0, 1, 0)]))
        with open(os.path.join(d, fn), "w") as f:
            f.write("hdr\n%d\n%d\n%s\n" % (NORB, len(ents),
                                            " ".join("1" * len(ents))))
            for i, (r, a, b) in enumerate(ents):
                val = 0.3 * (phase if i == 0 else np.conj(phase))
                f.write("%3d 0 0 %d %d %.12f %.12f\n"
                        % (r, a + 1, b + 1, val.real, val.imag))
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", kind: fn}
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 16},
                     "calc_scheme": "general"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        w = np.asarray(solver.ham_info.ham_inter_q)
        fz = getattr(solver.ham_info, "ham_fierz_q", None)
        return w, (None if fz is None else np.asarray(fz))

    # on-site declarations with a != b carry the longitudinal Fierz
    # (cross-slot) correction; the others do not build the tensor at all
    FIERZ_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange")

    def test_real_declarations_are_unchanged(self):
        for kind in self.REAL_TYPES:
            with self.subTest(type=kind):
                w, fz = self._vertices(kind)
                np.testing.assert_array_equal(
                    solver_rpa._to_bubble_pair_convention(w), w)
                self.assertGreater(np.abs(w).max(), 1e-12)
                if kind in self.FIERZ_TYPES:
                    # presence asserted, so that a regression dropping
                    # the tensor cannot pass by making it None
                    self.assertIsNotNone(fz)
                    self.assertGreater(np.abs(fz).max(), 1e-12)
                    np.testing.assert_array_equal(
                        solver_rpa._to_bubble_pair_convention(fz), fz)
                else:
                    self.assertTrue(fz is None or np.abs(fz).max() == 0.0)

    def test_offsite_density_declaration_is_unchanged(self):
        """Off-site density content is a fixed point as well. The
        property is structural -- the occupied slots are pair-diagonal,
        where swapping a pair's two indices does nothing -- so it holds
        for whatever bond phase W(q) carries (on this two-site lattice
        the phases happen to be real)."""
        w, _ = self._vertices("CoulombInter", offsite=True)
        self.assertGreater(np.abs(w).max(), 1e-12)
        np.testing.assert_array_equal(
            solver_rpa._to_bubble_pair_convention(w), w)
        nd = w.shape[-1]
        occupied = np.argwhere(np.abs(w.reshape(-1, *(nd,) * 4)) > 1e-12)
        for _, a, ap, b, bp in occupied:
            self.assertEqual((a, b), (ap, bp))

    def test_complex_pairhop_is_conjugated_on_crossing_slots(self):
        w, _ = self._vertices("PairHop", phase=PHASE)
        conv = solver_rpa._to_bubble_pair_convention(w)
        # mask taken from the OCCUPIED slots, not from the slots that
        # happened to move: a partial conversion must fail here
        occupied = np.abs(w) > 1e-12
        self.assertEqual(int(occupied.sum()), 4 * LX)
        np.testing.assert_allclose(conv[occupied],
                                   np.conj(w[occupied]), atol=1e-12)
        # and the declaration really is complex, or the pin is vacuous
        self.assertGreater(np.abs(w[occupied].imag).max(), 1e-3)


class TestFlexVertexNeedsNoConversion(unittest.TestCase):
    """FLEX's general path builds its spin/charge matrices directly in
    the bubble convention, so it must NOT be put through the helper as
    well: a future "let us make FLEX consistent" change that applied it
    a second time would conjugate a complex declaration back.

    The pin builds the S/C matrices through the very helpers
    ``FLEX._inflate_chi0q_and_ham_general`` calls and compares the slot
    carrying a complex PairHop against the converted RPA tensor.
    """

    def test_flex_sc_matches_converted_rpa_slot(self):
        from hwave.sc import _build_interaction_k
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = _write_inputs(tmp.name, "PairHopComplex", 1.0, PHASE)
        io = read_input_k.QLMSkInput(
            {"path_to_input": tmp.name, "interaction": inter})
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 16},
                     "calc_scheme": "general"}
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        no = solver.norb

        # FLEX general vertex, built exactly as the solver builds it
        kx = np.linspace(0, 2.0 * np.pi, LX, endpoint=False)
        k0 = np.linspace(0, 2.0 * np.pi, 1, endpoint=False)
        inter_k = _build_interaction_k(kx, k0, k0,
                                       solver.ham_info.param_ham, no)
        Us, Uc = build_sc_matrices_myo(inter_k, no, LX, 1, 1)
        Us = np.asarray(Us).reshape(LX, no * no, no * no)
        Uc = np.asarray(Uc).reshape(LX, no * no, no * no)

        # RPA vertex after the conversion
        conv = solver_rpa._to_bubble_pair_convention(
            np.asarray(solver.ham_info.ham_inter_q)).reshape(
            LX, ND, ND, ND, ND)

        # the declaration lives on the up/down pair-crossing slot
        rpa_slot = conv[0, _sb(0, 0), _sb(1, 0), _sb(1, 1), _sb(0, 1)]
        self.assertGreater(abs(rpa_slot.imag), 1e-3)
        # ... and on FLEX's orbital pair (0,1),(1,0) slot
        flex_slot = Us[0, 0 * no + 1, 1 * no + 0]
        np.testing.assert_allclose(flex_slot, rpa_slot, atol=1e-12)
        np.testing.assert_allclose(Uc[0, 0 * no + 1, 1 * no + 0],
                                   rpa_slot, atol=1e-12)
        # double-converting FLEX would land on the conjugate: assert the
        # two are actually distinguishable, or the pin is vacuous
        self.assertGreater(abs(flex_slot - np.conj(flex_slot)), 1e-3)


class TestTransverseComplexPairHop(unittest.TestCase):
    """The transverse (ladder) channel adjudicated DIRECTLY against
    exact diagonalization for a complex PairHop.

    The ring solve consumes the interaction in the bubble-pair
    convention (issue #139), but the transverse assembly re-pairs the
    tensor itself and must receive the HAMILTONIAN convention. The two
    routes are easy to conflate -- relating the spinful spin-flip slice
    to chiq_pm through an assumed pair-Hermiticity relation suggests the
    opposite answer -- so this pins the transverse channel on its own:
    chiq_pm's first-order response is compared against the exact
    transverse response of the declared Hamiltonian.

    The Lehmann tensor built here and the solver's spin-flip block are
    related by a PAIR-SLOT swap (see ``_ed_first_order_pm``), which is a
    different map from the longitudinal one -- reusing the longitudinal
    helper here is exactly the mistake this docstring exists to prevent.
    """

    def _chi_pm_ed(self, hint=None, h1=None):
        """Static transverse chi[q, a, c, b, d] over ORBITAL indices,
        built from the spin-flip bilinears directly."""
        H = np.zeros((DIM, DIM), dtype=complex)
        h1m = H1 if h1 is None else h1
        for p_ in range(NMODE):
            for q_ in range(NMODE):
                if h1m[p_, q_] != 0:
                    H += h1m[p_, q_] * (CD[p_] @ C[q_])
        if hint is not None:
            H = H + hint
        N = sum(CD[p_] @ C[p_] for p_ in range(NMODE))
        with np.errstate(all="ignore"):
            ev, V = np.linalg.eigh(H - MU * N)
        ev = ev - ev.min()
        w = np.exp(-BETA * ev)
        w /= w.sum()
        # O_pm[a, c] = sum_j e^{iqr} c^+_{j a up} c_{j c dn}
        # O_mp[d, b] = sum_j e^{-iqr} c^+_{j d dn} c_{j b up}
        Opm, Omp = {}, {}
        for qi in range(LX):
            for a in range(NORB):
                for c in range(NORB):
                    op = np.zeros((DIM, DIM), dtype=complex)
                    om = np.zeros((DIM, DIM), dtype=complex)
                    for j in range(LX):
                        ph = np.exp(2j * np.pi * qi * j / LX)
                        op += ph * (CD[_mode(j, a, 0)] @ C[_mode(j, c, 1)])
                        om += np.conj(ph) * (CD[_mode(j, a, 1)]
                                             @ C[_mode(j, c, 0)])
                    Opm[(qi, a, c)] = V.conj().T @ op @ V
                    Omp[(qi, a, c)] = V.conj().T @ om @ V
        dE = ev[None, :] - ev[:, None]
        dw = w[:, None] - w[None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            K = np.where(np.abs(dE) > 1e-10, dw / dE, BETA * w[:, None])
        out = np.zeros((LX, NORB, NORB, NORB, NORB), dtype=complex)
        for qi in range(LX):
            for a in range(NORB):
                for c in range(NORB):
                    A = Opm[(qi, a, c)]
                    for b in range(NORB):
                        for d in range(NORB):
                            B = Omp[(qi, d, b)]
                            out[qi, a, c, b, d] = (K * A * B.T).sum()
        return out / LX

    def _ed_first_order_pm(self, phase):
        base = self._chi_pm_ed()
        full = (2 * (self._chi_pm_ed(_h_int("PairHop", V1, phase)) - base)
                / V1
                - (self._chi_pm_ed(_h_int("PairHop", V2, phase)) - base)
                / V2)
        ins = (2 * (self._chi_pm_ed(h1=_hf_h1("PairHop", V1, phase))
                    - base) / V1
               - (self._chi_pm_ed(h1=_hf_h1("PairHop", V2, phase))
                  - base) / V2)
        # The Lehmann tensor above is
        #   R[a,c,b,d] = <(c+_{a up} c_{c dn})(q); (c+_{d dn} c_{b up})(-q)>
        # while the solver's spin-flip block is
        #   L[a,c,b,d] = <(c+_{c dn} c_{a up})(-q); (c+_{b up} c_{d dn})(q)>.
        # Static Kubo symmetry gives R[a,c,b,d] = L[b,d,a,c], so the
        # conversion is the PAIR-SLOT swap -- not the intra-pair swap the
        # longitudinal map uses (the two coincide only under accidental
        # symmetries, which a two-site fixture has).
        return np.transpose(full - ins, (0, 3, 4, 1, 2))

    def _solver_pm(self, v, phase):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        inter = _write_inputs(tmp.name, "PairHopComplex", v, phase)
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "mu": MU, "CellShape": [LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 512,
                               "coeff_tail": 1.0},
                     "calc_scheme": "general", "calc_type": "ring+ladder"}
        io = read_input_k.QLMSkInput(
            {"path_to_input": tmp.name, "interaction": inter})
        solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
        gi = io.get_param("green")
        out = os.path.join(tmp.name, "out")
        os.makedirs(out, exist_ok=True)
        solver.solve(gi, out)
        arr = np.asarray(gi["chiq_pm"])
        nf = arr.shape[0]
        return arr.reshape(nf, LX, NORB, NORB, NORB, NORB)[nf // 2]

    def test_pair_slot_map_is_the_one_kubo_symmetry_selects(self):
        """Pin the CHOICE of map algebraically.

        On the two-site fixture below the pair-slot map and the
        longitudinal intra-pair map happen to agree to round-off
        (every momentum is self-inverse there), so the numerical
        comparison cannot show which one is right. The relation they
        differ on is exact and needs no lattice: for the transverse
        tensor, static Kubo symmetry gives R[a,c,b,d] = L[b,d,a,c],
        i.e. the solver frame is reached by swapping the two PAIR
        slots. A generic tensor carrying that symmetry satisfies the
        relation under the pair-slot map only."""
        rng = np.random.RandomState(139)
        n = NORB
        r = (rng.randn(LX, n, n, n, n) + 1j * rng.randn(LX, n, n, n, n))
        # impose the Kubo/pair-Hermiticity structure R[a,c,b,d] taken to
        # the solver frame L by the pair-slot swap
        L = np.transpose(r, (0, 3, 4, 1, 2))
        np.testing.assert_allclose(
            np.transpose(L, (0, 3, 4, 1, 2)), r, atol=1e-14)
        # the longitudinal intra-pair map does NOT reproduce it
        wrong = np.transpose(r, (0, 2, 1, 4, 3))
        self.assertGreater(np.abs(wrong - L).max(), 1e-3)

    def test_ladder_uses_the_hamiltonian_convention(self):
        c0 = self._solver_pm(1e-9, PHASE)
        sol = (2 * (self._solver_pm(V1, PHASE) - c0) / V1
               - (self._solver_pm(V2, PHASE) - c0) / V2)
        ed_decl = self._ed_first_order_pm(PHASE)
        ed_conj = self._ed_first_order_pm(np.conj(PHASE))
        shared = (np.abs(ed_decl) > 1e-6) & (np.abs(sol) > 1e-6)
        # pinned so that a solver returning only a couple of correct
        # elements cannot satisfy the comparison by shrinking the mask
        self.assertEqual(int(shared.sum()), 32)
        scale = np.abs(ed_decl).max()
        self.assertGreater(scale, 1e-3)
        dev = np.abs((ed_decl - sol)[shared]).max()
        self.assertLess(dev, 5e-3 * scale)
        # and it must not invent slots the exact response leaves empty,
        # nor drop ones it fills
        only_sol = (np.abs(sol) > 1e-6) & (np.abs(ed_decl) <= 1e-6)
        self.assertLess(np.abs(sol[only_sol]).max(initial=0.0),
                        5e-3 * scale)
        missing = (np.abs(ed_decl) > 1e-3 * scale) & (np.abs(sol) <= 1e-6)
        self.assertEqual(int(missing.sum()), 0)
        # the conjugate model must be clearly excluded, or the pin says
        # nothing about the orientation
        shared_c = (np.abs(ed_conj) > 1e-6) & (np.abs(sol) > 1e-6)
        self.assertGreater(np.abs((ed_conj - sol)[shared_c]).max(),
                           1e-1 * scale)


class TestComplexPairHopGpuMatchesCpu(unittest.TestCase):
    """The conversion must behave identically on the device path. The
    committed GPU tests use real declarations, on which it is the
    identity, so they cannot see it."""

    def test_cpu_gpu(self):
        try:
            import cupy
            cupy.zeros(1)
        except Exception:
            self.skipTest("cupy unavailable or no usable CUDA device")
        out = {}
        for gpu in (False, True):
            tmp = tempfile.TemporaryDirectory()
            self.addCleanup(tmp.cleanup)
            inter = _write_inputs(tmp.name, "PairHopComplex", 0.3, PHASE)
            info_mode = {"mode": "RPA",
                         "param": {"T": T, "mu": MU,
                                   "CellShape": [LX, 1, 1],
                                   "SubShape": [1, 1, 1], "Nmat": 64,
                                   "coeff_tail": 1.0, "gpu": gpu},
                         "calc_scheme": "general"}
            io = read_input_k.QLMSkInput(
                {"path_to_input": tmp.name, "interaction": inter})
            solver = solver_rpa.RPA(io.get_param("ham"), {}, info_mode)
            gi = io.get_param("green")
            o = os.path.join(tmp.name, "out")
            os.makedirs(o, exist_ok=True)
            solver.solve(gi, o)
            out[gpu] = np.asarray(gi["chiq"])
        np.testing.assert_allclose(out[True], out[False], atol=1e-10)


if __name__ == "__main__":
    unittest.main()
