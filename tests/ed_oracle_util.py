#!/usr/bin/env python3

"""Lehmann-representation exact-diagonalization machinery, extracted from
tests/test_rpa_vs_ed_oracle.py behind a round-off gate (#151).

This module carries the mechanics that are shared across every ED oracle
fixture: building the free one-body Hamiltonian and annihilation operators
for an L-site, norb-orbital, spin-1/2 chain, diagonalizing a Hamiltonian in
the Lehmann representation, evaluating the static connected density-density
correlator, mapping it into the solver's slot convention, and Richardson-
extrapolating a first-order response in a coupling constant.

``EDFixture`` is fully lazy: constructing one costs nothing beyond storing
its parameters. The exponential-in-Hilbert-space work (building the
annihilation operators, diagonalizing) happens only when a caller asks for
it, and is cached on the instance thereafter.
"""

import numpy as np


class EDFixture:
    """A two-spin, ``norb``-orbital, ``L``-site tight-binding chain, as an
    exact-diagonalization fixture.

    Parameters
    ----------
    L : int
        Number of sites (periodic chain).
    norb : int
        Number of orbitals per site.
    t : dict[(int, int), complex]
        Nearest-neighbour hopping amplitudes, keyed by orbital pair
        ``(a, b)``, both directions present (``t[(b, a)]`` is read for the
        Hermitian-conjugate hop). Spin-independent.
    eps : sequence[float]
        On-site energy per orbital, length ``norb``.
    T : float
        Temperature.
    mu : float
        Chemical potential.

    All heavy state (annihilation operators, the one-body Hamiltonian) is
    built lazily on first use and cached.
    """

    def __init__(self, L, norb, t, eps, T, mu):
        self.L = L
        self.norb = norb
        self.t = t
        self.eps = eps
        self.T = T
        self.mu = mu
        self.beta = 1.0 / T
        self.nd = norb * 2                  # generalized orbital (spin-block)
        self.nmode = L * self.nd
        self.dim = 1 << self.nmode
        self._h1 = None
        self._c = None

    def mode(self, j, o, s):
        """Single-particle mode index (site-major, interleaved 2*orb+spin)."""
        return self.nd * j + 2 * o + s

    def build_h1(self):
        """One-body Hamiltonian in mode space, (nmode, nmode). Cached."""
        if self._h1 is not None:
            return self._h1
        h = np.zeros((self.nmode, self.nmode), dtype=complex)
        for j in range(self.L):
            jp = (j + 1) % self.L
            for (a, b), tv in self.t.items():
                for s in range(2):
                    h[self.mode(jp, a, s), self.mode(j, b, s)] += tv
                    h[self.mode(j, a, s), self.mode(jp, b, s)] += \
                        np.conj(self.t[(b, a)])
            for o, e in enumerate(self.eps):
                for s in range(2):
                    h[self.mode(j, o, s), self.mode(j, o, s)] += e
        self._h1 = h
        return h

    def annihilators(self):
        """Dense annihilation operators, one per mode, as (dim, dim)
        matrices in the occupation-number basis. Cached."""
        if self._c is not None:
            return self._c
        dim = self.dim
        ops = []
        for p in range(self.nmode):
            c = np.zeros((dim, dim))
            for st in range(dim):
                if (st >> p) & 1:
                    sign = (-1) ** bin(st & ((1 << p) - 1)).count("1")
                    c[st ^ (1 << p), st] = sign
            ops.append(c.astype(complex))
        self._c = ops
        return ops


def _diagonalize(fx, hint=None, h1=None):
    """Diagonalize H1(+hint) - mu*N in the many-body Hilbert space.

    Returns (ev, w, V): eigenvalues shifted to a zero ground state, the
    normalized Boltzmann weights, and the eigenvector matrix.
    """
    C = fx.annihilators()
    CD = [c.conj().T for c in C]
    h1m = fx.build_h1() if h1 is None else h1
    H = np.zeros((fx.dim, fx.dim), dtype=complex)
    for p in range(fx.nmode):
        for q in range(fx.nmode):
            if h1m[p, q] != 0:
                H += h1m[p, q] * (CD[p] @ C[q])
    if hint is not None:
        H = H + hint
    N = sum(CD[p] @ C[p] for p in range(fx.nmode))
    with np.errstate(all="ignore"):
        ev, V = np.linalg.eigh(H - fx.mu * N)
    ev = ev - ev.min()
    w = np.exp(-fx.beta * ev)
    w /= w.sum()
    return ev, w, V


def lehmann_kernel(fx, hint=None, h1=None):
    """The static Lehmann kernel K[m, n] = (w_m - w_n)/(e_m - e_n), with
    the thermal-limit value beta*w_m on (near-)degenerate pairs, for
    H1(+hint) - mu*N."""
    ev, w, _V = _diagonalize(fx, hint=hint, h1=h1)
    dE = ev[None, :] - ev[:, None]
    dw = w[:, None] - w[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(dE) > 1e-10, dw / dE, fx.beta * w[:, None])
    return K


def chi_connected(fx, hint=None, h1=None):
    """Static chi[q, a, c, b, d] = <(c^+_a c_c)(q); (c^+_d c_b)(-q)>,
    connected, by the Lehmann representation, for H1(+hint) - mu*N.

    ``a``, ``c``, ``b``, ``d`` are generalized-orbital (spin-block, per
    ``fx.nd``) indices, NOT physical-orbital indices.
    """
    C = fx.annihilators()
    CD = [c.conj().T for c in C]
    ev, w, V = _diagonalize(fx, hint=hint, h1=h1)
    O = {}
    for qi in range(fx.L):
        for a in range(fx.nd):
            for c in range(fx.nd):
                oa, sa = a % fx.norb, a // fx.norb
                oc, sc = c % fx.norb, c // fx.norb
                op = np.zeros((fx.dim, fx.dim), dtype=complex)
                for j in range(fx.L):
                    ph = np.exp(2j * np.pi * qi * j / fx.L)
                    op += ph * (CD[fx.mode(j, oa, sa)] @ C[fx.mode(j, oc, sc)])
                O[(qi, a, c)] = V.conj().T @ op @ V
    dE = ev[None, :] - ev[:, None]
    dw = w[:, None] - w[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        K = np.where(np.abs(dE) > 1e-10, dw / dE, fx.beta * w[:, None])
    out = np.zeros((fx.L, fx.nd, fx.nd, fx.nd, fx.nd), dtype=complex)
    for qi in range(fx.L):
        qn = (-qi) % fx.L
        for a in range(fx.nd):
            for c in range(fx.nd):
                A = O[(qi, a, c)]
                for b in range(fx.nd):
                    for d in range(fx.nd):
                        B = O[(qn, d, b)]
                        val = (K * A * B.T).sum()
                        if qi == 0:
                            val -= fx.beta * (w * np.diag(A)).sum() \
                                * (w * np.diag(B)).sum()
                        out[qi, a, c, b, d] = val
    return out / fx.L


def to_solver_slots(x):
    """DERIVED map from the Lehmann layout to the solver's slots (swap the
    two indices inside each pair; see test_rpa_vs_ed_oracle's module
    docstring for the derivation)."""
    return np.transpose(x, (0, 2, 1, 4, 3))


def richardson(fn_of_v, v1):
    """Two-point Richardson extrapolation of d(fn)/dv at v -> 0:

        2*(fn(v1) - fn(0))/v1 - (fn(2*v1) - fn(0))/(2*v1)

    ED SIDE ONLY -- the solver-side prediction is exactly linear in the
    coupling and has no Richardson stencil (plan review, must-fix 7)."""
    f0 = fn_of_v(0.0)
    f1 = fn_of_v(v1)
    f2 = fn_of_v(2 * v1)
    return 2 * (f1 - f0) / v1 - (f2 - f0) / (2 * v1)


def canonical_density_terms(fx, pairs):
    """CANONICAL list. pairs = [(a, b, R, coeff)] with a, b PHYSICAL
    ORBITAL indices (0..norb-1, NEVER generalized spin-orbital indices --
    round-3 fix: the two domains are now named and distinct throughout;
    generalized indices appear ONLY in correlator layouts). One entry per
    unordered class. Spin expansion, exact:
      - a == b and R == 0 (CoulombIntra U): per site, the single term
        U * n_{j,a,up} n_{j,a,dn}  (matches _h_int).
      - otherwise (CoulombInter V, on- or off-site, any orbitals): per
        site, ALL FOUR spin pairings
        V * n_{j,a,sigma} n_{j+R,b,sigma'} for sigma, sigma' in
        {up, dn} x {up, dn}  (the spin-summed density product).
    Element-level pins for U, single-orbital V, and case-M inter-orbital
    V are part of the Task-3 tests. The one-body t matrix (orbital
    domain) lifts to the spin-orbital H1 as t (x) I_spin via mode(j,o,s).

    Returns a list of (p, q, r, s, coeff) quartic terms, each representing
    coeff * c^+_p c_q c^+_r c_s (here always a density-density product,
    p == q and r == s), suitable for ``h_int_from_terms`` and
    ``hf_h1_from_terms``.
    """
    terms = []
    for (a, b, R, coeff) in pairs:
        if a == b and R == 0:
            for j in range(fx.L):
                p = fx.mode(j, a, 0)
                r = fx.mode(j, a, 1)
                terms.append((p, p, r, r, coeff))
        else:
            for j in range(fx.L):
                jp = (j + R) % fx.L
                for s1 in range(2):
                    for s2 in range(2):
                        p = fx.mode(j, a, s1)
                        r = fx.mode(jp, b, s2)
                        terms.append((p, p, r, r, coeff))
    return terms


def h_int_from_terms(fx, terms):
    """Build sum_terms coeff * c^+_p c_q c^+_r c_s in the many-body
    Hilbert space, for ``terms`` as produced by
    ``canonical_density_terms`` (or any other (p, q, r, s, coeff)
    quartic list over mode indices)."""
    C = fx.annihilators()
    CD = [c.conj().T for c in C]
    H = np.zeros((fx.dim, fx.dim), dtype=complex)
    for (p, q, r, s, coeff) in terms:
        H += coeff * (CD[p] @ C[q] @ CD[r] @ C[s])
    return H


def hf_h1_from_terms(fx, terms):
    """H1 plus the first-order Hartree-Fock self-energy of the interaction
    ``terms`` describe, by Wick contraction with the free density matrix
    (generalized ``add()`` engine, works for any quartic term list --
    on-site or displaced across bonds).
    """
    h1 = fx.build_h1()
    ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
    f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
    rho = ((W * f) @ W.conj().T).T          # rho[p, q] = <c^+_p c_q>
    S = np.zeros((fx.nmode, fx.nmode), dtype=complex)

    def add(p, q, r, s, c_):
        # E = c [rho_pq rho_rs + rho_ps (delta_qr - rho_rq)];
        # H_MF = sum (dE / d rho_xy) c^+_x c_y
        S[p, q] += c_ * rho[r, s]
        S[r, s] += c_ * rho[p, q]
        S[p, s] += c_ * ((1.0 if q == r else 0.0) - rho[r, q])
        S[r, q] += -c_ * rho[p, s]

    for (p, q, r, s, coeff) in terms:
        add(p, q, r, s, coeff)
    return h1 + 0.5 * (S + S.conj().T)
