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


def _static_kernel(ev, w, beta):
    """The static Lehmann kernel K[m, n] = (w_m - w_n)/(e_n - e_m), with
    the thermal-limit value beta*w_m on (near-)degenerate pairs.

    Shared by ``lehmann_kernel`` and ``chi_connected`` -- one copy so the
    two can never drift apart.
    """
    dE = ev[None, :] - ev[:, None]
    dw = w[:, None] - w[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(np.abs(dE) > 1e-10, dw / dE, beta * w[:, None])


def lehmann_kernel(fx, hint=None, h1=None):
    """The static Lehmann kernel K[m, n] = (w_m - w_n)/(e_n - e_m), with
    the thermal-limit value beta*w_m on (near-)degenerate pairs, for
    H1(+hint) - mu*N."""
    ev, w, _V = _diagonalize(fx, hint=hint, h1=h1)
    return _static_kernel(ev, w, fx.beta)


def chi_connected(fx, hint=None, h1=None):
    """Static chi[q, a, c, b, d] = <(c^+_a c_c)(q); (c^+_d c_b)(-q)>,
    connected, by the Lehmann representation, for H1(+hint) - mu*N.

    ``a``, ``c``, ``b``, ``d`` are generalized-orbital (spin-block, per
    ``fx.nd``) indices, NOT physical-orbital indices.

    The disconnected piece <O><O'> is subtracted at EVERY q (not just
    q=0): under translation invariance the q!=0 bilinears have zero
    thermal average and the subtraction is a no-op there, but a
    Hamiltonian passed in via ``hint``/``h1`` need not be translation
    invariant, and an unconditional subtraction is the only form that
    stays correct for that case too.
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
    K = _static_kernel(ev, w, fx.beta)
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


# ---------------------------------------------------------------------------
# Sector-block ED engine (#151, Task 4)
# ---------------------------------------------------------------------------

def _apply_ops(state, ops):
    """Apply a sequence of second-quantized (mode, dagger) operators to the
    Fock ket ``state`` (an int bitmask over mode indices), IN THE GIVEN
    ORDER -- ``ops[0]`` is the operator applied FIRST, i.e. the RIGHTMOST
    factor of the operator string (matching matrix-vector application
    order: ``CD[p] @ C[q] @ |state>`` applies ``C[q]`` first). Returns
    ``(new_state, sign)``, or ``None`` if any step is annihilation on an
    empty slot or creation on an occupied one. Sign convention matches
    ``EDFixture.annihilators()`` exactly (same ``bin(...).count("1")``
    parity rule), so this and the dense path can never silently diverge on
    fermion sign."""
    sign = 1
    for mode, dagger in ops:
        bit = (state >> mode) & 1
        if dagger and bit:
            return None
        if not dagger and not bit:
            return None
        sign *= (-1) ** bin(state & ((1 << mode) - 1)).count("1")
        state ^= (1 << mode)
    return state, sign


def _apply_pq(state, p, q):
    """c^dag_p c_q |state>, matching h_int_from_terms's CD[p] @ C[q]
    matrix-product order (q, the annihilation/rightmost factor, applied
    first)."""
    return _apply_ops(state, ((q, False), (p, True)))


def _apply_pqrs(state, p, q, r, s):
    """c^dag_p c_q c^dag_r c_s |state>, matching h_int_from_terms's
    CD[p] @ C[q] @ CD[r] @ C[s] order."""
    return _apply_ops(state, ((s, False), (r, True), (q, False), (p, True)))


def _apply_pair(state, p, q):
    """c_p c_q |state> (both annihilation; q, the rightmost factor,
    applied first) -- matches Delta_{R,a,b}(q)'s c_{...,down} c_{...,up}
    order when p is the down-spin mode and q the up-spin mode."""
    return _apply_ops(state, ((q, False), (p, False)))


class SectorED:
    """Exact diagonalization by (N_up, N_dn) symmetry block, for a
    Hamiltonian H1(+terms) - mu*(N_up+N_dn) that conserves both particle
    numbers separately -- true of the free chain always, of the CANONICAL
    density ``terms`` always, and of an ``h1`` built by
    ``hf_h1_from_terms`` up to floating-point noise on cross-spin entries
    (the free density matrix it Wick-contracts against is itself exactly
    spin-block-diagonal, since ``build_h1`` never mixes spins; see
    ``_SECTOR_FLOOR`` below, which turns any non-noise-level violation
    into a loud failure rather than a silently wrong block).

    This replaces the dense ``(dim, dim)`` path for systems too large to
    hold as dense operators: case M (L=3, norb=2) has dim=4096, and 12
    dense annihilation matrices alone (``EDFixture.annihilators()``) would
    be multiple GB. The largest symmetry block here is
    (N_up, N_dn) = (3, 3), dim 400 -- dense ``eigh`` per block is trivial.

    Boltzmann weights and cross-sector Lehmann denominators use the
    GRAND-CANONICAL spectrum: every sector's K = H_sector -
    mu*(N_up+N_dn) is diagonalized independently, but the eigenvalues are
    shifted by the GLOBAL minimum (over every sector, not each sector's
    own) before the Boltzmann weights are normalized over the FULL Fock
    space -- exactly what the dense ``_diagonalize`` does with one big
    matrix, just computed block by block.

    Bond and density bilinears whose two legs share creation/annihilation
    spin are block-DIAGONAL in (N_up, N_dn); density bilinears with
    opposite-spin legs (the transverse/spin-flip generalized-orbital
    entries ``chi_connected`` also has to produce) and the pair operator
    (always N_up,N_dn -> N_up-1,N_dn-1) are block-OFF-diagonal. Two
    Lehmann-sum helpers cover both cases with NO special-casing for the
    trivial (0, 0) shift: ``_lehmann_transpose`` (the (K*A*B.T) contraction
    ``chi_connected`` needs) and ``_lehmann_dagger`` (the <A;B^dagger>
    contraction the bond/pair correlators need) both build their kernel by
    concatenating the two involved sectors' (grand-canonical) eigenvalues
    and weights and slicing the cross block of the ordinary
    ``_static_kernel`` -- when the two sectors coincide this is elementwise
    identical to ``_static_kernel`` applied to that one sector directly
    (both blocks are read off the very same underlying (ev, w) arrays by
    index), so the in-sector case falls out of the general formula for
    free.

    Operators are built PER (q, correlator) call and discarded once
    consumed -- never the full cross product of q and channel held at
    once. One density bilinear's cross-sector element count already sums
    to ~8.5e5 across case M's 49 sectors (``sum_k C(6,k)^2 = C(12,6)``,
    squared over the independent N_up/N_dn sums); holding all of them for
    every (q, a, c) combination simultaneously would be several hundred
    MB, so at most a "one q's worth of one operator role" cache is kept
    (see ``chi_connected``'s ``b_side`` and ``bond_correlator``/
    ``pair_correlator``'s ``built`` lists), never the full set.
    """

    _SECTOR_FLOOR = 1e-9   # max tolerated (N_up, N_dn)-non-conserving H
                            # matrix element before it is treated as a real
                            # bug rather than floating-point noise.

    def __init__(self, fx, terms=(), h1=None):
        self.fx = fx
        self.terms = list(terms)
        self.h1 = fx.build_h1() if h1 is None else h1
        self._build_sectors()

    # -- sector bookkeeping ------------------------------------------------

    def _build_sectors(self):
        fx = self.fx
        up_mask = sum(1 << p for p in range(fx.nmode) if p % 2 == 0)
        full_mask = (1 << fx.nmode) - 1
        sector_states = {}
        for state in range(fx.dim):
            nup = bin(state & up_mask).count("1")
            ndn = bin(state & ~up_mask & full_mask).count("1")
            sector_states.setdefault((nup, ndn), []).append(state)
        sector_index = {
            key: {s: i for i, s in enumerate(states)}
            for key, states in sector_states.items()
        }
        self._sector_states = sector_states
        self._sector_index = sector_index

        raw_ev, V = {}, {}
        for key, states in sector_states.items():
            Hloc = self._sector_hamiltonian(states, sector_index[key])
            nup, ndn = key
            Kloc = Hloc - fx.mu * (nup + ndn) * np.eye(len(states))
            with np.errstate(all="ignore"):
                ev_k, V_k = np.linalg.eigh(Kloc)
            raw_ev[key] = ev_k
            V[key] = V_k
        gmin = min(ev_k.min() for ev_k in raw_ev.values())
        ev = {key: ev_k - gmin for key, ev_k in raw_ev.items()}
        w_unnorm = {key: np.exp(-fx.beta * ev_k) for key, ev_k in ev.items()}
        z = sum(w_k.sum() for w_k in w_unnorm.values())
        w = {key: w_k / z for key, w_k in w_unnorm.items()}
        self._ev, self._w, self._V = ev, w, V

    def _sector_hamiltonian(self, states, index_map):
        """Dense (M, M) H = H1 + terms restricted to one (N_up, N_dn)
        sector, via direct bit-level second quantization (never through
        ``fx.annihilators()``'s dense (dim, dim) matrices -- that dense
        path is exactly what the block engine exists to avoid)."""
        fx = self.fx
        M = len(states)
        Hloc = np.zeros((M, M), dtype=complex)
        dropped = 0.0
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                coeff = self.h1[p, q]
                if coeff == 0:
                    continue
                for col, state in enumerate(states):
                    r = _apply_pq(state, p, q)
                    if r is None:
                        continue
                    new_state, sign = r
                    row = index_map.get(new_state)
                    if row is None:
                        dropped = max(dropped, abs(coeff))
                        continue
                    Hloc[row, col] += coeff * sign
        for (p, q, r, s, coeff) in self.terms:
            for col, state in enumerate(states):
                res = _apply_pqrs(state, p, q, r, s)
                if res is None:
                    continue
                new_state, sign = res
                row = index_map.get(new_state)
                if row is None:
                    dropped = max(dropped, abs(coeff))
                    continue
                Hloc[row, col] += coeff * sign
        if dropped > self._SECTOR_FLOOR:
            raise AssertionError(
                "SectorED: h1/terms have a (N_up, N_dn)-non-conserving "
                "matrix element of magnitude {:.3e} -- SectorED requires an "
                "exactly number-conserving Hamiltonian".format(dropped))
        return Hloc

    # -- generic bilinear machinery -----------------------------------------

    def _build_operator(self, mode_terms):
        """mode_terms: [(p, q, weight), ...] for sum_j weight_j *
        c^dag_p c_q. Returns (delta, op): delta = (dNup, dNdn) this
        bilinear imparts (the SAME for every entry -- a site sum keeps
        orbital/spin fixed, only the site index varies, so this is
        computed once and asserted consistent); op maps each ket sector
        key PRESENT in the fixture to its dense block in the SECTOR
        EIGENBASIS (bra sector = ket sector + delta; sectors with no such
        bra partner are simply absent from the dict, not zero-filled)."""
        delta = None
        for (p, q, _weight) in mode_terms:
            d = ((p % 2 == 0) - (q % 2 == 0), (p % 2 == 1) - (q % 2 == 1))
            if delta is None:
                delta = d
            elif d != delta:
                raise ValueError("mode_terms mix incompatible sector shifts")
        op = {}
        for key_ket, states_ket in self._sector_states.items():
            key_bra = (key_ket[0] + delta[0], key_ket[1] + delta[1])
            states_bra = self._sector_states.get(key_bra)
            if states_bra is None:
                continue
            idx_bra = self._sector_index[key_bra]
            block = np.zeros((len(states_bra), len(states_ket)), dtype=complex)
            for col, state in enumerate(states_ket):
                for (p, q, weight) in mode_terms:
                    r = _apply_pq(state, p, q)
                    if r is None:
                        continue
                    new_state, sign = r
                    block[idx_bra[new_state], col] += weight * sign
            op[key_ket] = self._V[key_bra].conj().T @ block @ self._V[key_ket]
        return delta, op

    def _build_pair_operator(self, mode_terms):
        """mode_terms: [(p_dn, q_up, weight), ...] for sum_j weight_j *
        c_{p_dn} c_{q_up} (both annihilation; matches Delta_{R,a,b}(q)'s
        c_down c_up order). Always delta = (-1, -1). Same block/eigenbasis
        convention as ``_build_operator``."""
        delta = (-1, -1)
        op = {}
        for key_ket, states_ket in self._sector_states.items():
            key_bra = (key_ket[0] - 1, key_ket[1] - 1)
            states_bra = self._sector_states.get(key_bra)
            if states_bra is None:
                continue
            idx_bra = self._sector_index[key_bra]
            block = np.zeros((len(states_bra), len(states_ket)), dtype=complex)
            for col, state in enumerate(states_ket):
                for (p, q, weight) in mode_terms:
                    r = _apply_pair(state, p, q)
                    if r is None:
                        continue
                    new_state, sign = r
                    block[idx_bra[new_state], col] += weight * sign
            op[key_ket] = self._V[key_bra].conj().T @ block @ self._V[key_ket]
        return delta, op

    def _thermal_avg(self, op, delta):
        if delta != (0, 0):
            return 0.0
        total = 0.0
        for key, block in op.items():
            total += (self._w[key] * np.diag(block)).sum()
        return total

    def _pair_kernel(self, key_m, key_n):
        """K[m, n] across sectors key_m (rows) and key_n (cols), via the
        SAME ``_static_kernel`` the dense path uses: concatenate the two
        sectors' (grand-canonical) eigenvalues/weights into one array and
        slice the cross block. When key_m == key_n this is elementwise
        identical to ``_static_kernel`` on that one sector alone (both
        halves of the concatenation read off the same (ev, w) by index),
        so the in-sector case is not special-cased."""
        ev_m, w_m = self._ev[key_m], self._w[key_m]
        ev_n, w_n = self._ev[key_n], self._w[key_n]
        ev_cat = np.concatenate([ev_m, ev_n])
        w_cat = np.concatenate([w_m, w_n])
        M = len(ev_m)
        return _static_kernel(ev_cat, w_cat, self.fx.beta)[:M, M:]

    def _lehmann_transpose(self, opA, deltaA, opB, deltaB):
        """sum_{m,n} K[m,n] A[m,n] B[n,m], for A: ket n -> bra n+deltaA,
        B: ket m -> bra m+deltaB. Nonzero only where deltaB == -deltaA
        (the only way A's bra sector can be B's ket sector, and vice
        versa); zero otherwise, with no special-casing needed."""
        if deltaB[0] != -deltaA[0] or deltaB[1] != -deltaA[1]:
            return 0.0
        total = 0.0
        for key_n, blockA in opA.items():
            key_m = (key_n[0] + deltaA[0], key_n[1] + deltaA[1])
            blockB = opB.get(key_m)
            if blockB is None:
                continue
            K = self._pair_kernel(key_m, key_n)
            total += (K * blockA * blockB.T).sum()
        return total

    def _lehmann_dagger(self, opA, deltaA, opB, deltaB):
        """sum_{m,n} K[m,n] A[m,n] conj(B[m,n]), i.e. <A; B^dagger>; both A
        and B map ket -> ket+delta, so a nonzero result needs
        deltaA == deltaB."""
        if deltaA != deltaB:
            return 0.0
        total = 0.0
        for key_n, blockA in opA.items():
            key_m = (key_n[0] + deltaA[0], key_n[1] + deltaA[1])
            blockB = opB.get(key_n)
            if blockB is None:
                continue
            K = self._pair_kernel(key_m, key_n)
            total += (K * blockA * np.conj(blockB)).sum()
        return total

    def _channel_orbitals(self, chan):
        if self.fx.norb == 1:
            (R,) = chan
            return R, 0, 0
        R, a, b = chan
        return R, a, b

    # -- public correlators --------------------------------------------

    def chi_connected(self):
        """Same layout/normalization as the module-level dense
        ``chi_connected``: chi[q, a, c, b, d], generalized-orbital
        indices, 1/L, disconnected subtraction at every q over every
        nonzero-average bilinear (including the mixed-spin/transverse
        generalized-orbital entries, computed exactly via the cross-sector
        Lehmann sum -- unlike bond/pair bilinears these are NOT always
        in-sector, since a, c range over both spins)."""
        fx = self.fx
        out = np.zeros((fx.L, fx.nd, fx.nd, fx.nd, fx.nd), dtype=complex)
        for qi in range(fx.L):
            qn = (-qi) % fx.L
            # One q's worth of the (d, b) role, built ONCE and reused
            # across every (a, c) -- never the full (q, a, c, b, d) product
            # held at once (see class docstring, case-M memory).
            b_side = {}
            for d in range(fx.nd):
                od, sd = d % fx.norb, d // fx.norb
                for b in range(fx.nd):
                    ob, sb = b % fx.norb, b // fx.norb
                    mode_terms = [
                        (fx.mode(j, od, sd), fx.mode(j, ob, sb),
                         np.exp(2j * np.pi * qn * j / fx.L))
                        for j in range(fx.L)]
                    delta, op = self._build_operator(mode_terms)
                    b_side[(d, b)] = (delta, op, self._thermal_avg(op, delta))
            for a in range(fx.nd):
                oa, sa = a % fx.norb, a // fx.norb
                for c in range(fx.nd):
                    oc, sc = c % fx.norb, c // fx.norb
                    mode_terms = [
                        (fx.mode(j, oa, sa), fx.mode(j, oc, sc),
                         np.exp(2j * np.pi * qi * j / fx.L))
                        for j in range(fx.L)]
                    deltaA, opA = self._build_operator(mode_terms)
                    avgA = self._thermal_avg(opA, deltaA)
                    for b in range(fx.nd):
                        for d in range(fx.nd):
                            deltaB, opB, avgB = b_side[(d, b)]
                            val = self._lehmann_transpose(opA, deltaA, opB, deltaB)
                            val -= fx.beta * avgA * avgB
                            out[qi, a, c, b, d] = val
        return out / fx.L

    def bond_correlator(self, channels):
        """Xph[q, sigma, sigma', I, J], connected, of B_{R,a,b,sigma}(q) =
        sum_j e^{iqj} c^dag_{j+R,a,sigma} c_{j,b,sigma} against the
        conjugate leg B^dagger_{R',a',b',sigma'}(q) (SAME q on both legs --
        the physical <B; B^dagger> susceptibility, NOT chi_connected's
        independent-(-q)-leg tensor; see the module docstring). ``I``, ``J``
        index the EXPLICIT ordered ``channels`` list ``[(R, a, b), ...]``
        (``[(R,)]`` at norb=1). Both legs of B share creation/annihilation
        spin, so every channel operator is block-diagonal (delta ==
        (0, 0)) regardless of sigma -- bond bilinears never leave their
        sector, unlike chi_connected's mixed-spin entries. 1/L-normalized,
        matching chi_connected (see the (m=0, a=b) internal-consistency
        pin in tests/test_bond_vs_ed_oracle.py)."""
        fx = self.fx
        n_i = len(channels)
        out = np.zeros((fx.L, 2, 2, n_i, n_i), dtype=complex)
        for qi in range(fx.L):
            built = []   # built[I][sigma] = (delta, op, avg)
            for chan in channels:
                R, a, b = self._channel_orbitals(chan)
                row = []
                for sigma in range(2):
                    mode_terms = [
                        (fx.mode((j + R) % fx.L, a, sigma),
                         fx.mode(j, b, sigma),
                         np.exp(2j * np.pi * qi * j / fx.L))
                        for j in range(fx.L)]
                    delta, op = self._build_operator(mode_terms)
                    row.append((delta, op, self._thermal_avg(op, delta)))
                built.append(row)
            for i in range(n_i):
                for sigma in range(2):
                    deltaA, opA, avgA = built[i][sigma]
                    for j in range(n_i):
                        for sigmap in range(2):
                            deltaB, opB, avgB = built[j][sigmap]
                            val = self._lehmann_dagger(opA, deltaA, opB, deltaB)
                            val -= fx.beta * avgA * np.conj(avgB)
                            out[qi, sigma, sigmap, i, j] = val
        return out / fx.L

    def pair_correlator(self, channels):
        """Xpp[q, i, j], connected, of Delta_{R,a,b}(q) = sum_j e^{iqj}
        c_{j+R,a,down} c_{j,b,up} against Delta^dagger(q); the ANNIHILATOR
        Delta maps (N_up, N_dn) -> (N_up-1, N_dn-1), always, so every
        channel's operator couples exactly two adjacent sectors via the
        same cross-sector machinery ``chi_connected``'s mixed-spin entries
        use. ``i``, ``j`` index the EXPLICIT ordered ``channels`` list
        ``[(R, a, b), ...]`` (``[(R,)]`` at norb=1). NO disconnected piece
        -- <Delta> is exactly 0 by number conservation, so it is not
        computed at all (rather than subtracted as a measured zero). 1/L-
        normalized, matching chi_connected/bond_correlator."""
        fx = self.fx
        n_i = len(channels)
        out = np.zeros((fx.L, n_i, n_i), dtype=complex)
        for qi in range(fx.L):
            built = []
            for chan in channels:
                R, a, b = self._channel_orbitals(chan)
                mode_terms = [
                    (fx.mode((j + R) % fx.L, a, 1), fx.mode(j, b, 0),
                     np.exp(2j * np.pi * qi * j / fx.L))
                    for j in range(fx.L)]
                built.append(self._build_pair_operator(mode_terms))
            for i in range(n_i):
                deltaA, opA = built[i]
                for j in range(n_i):
                    deltaB, opB = built[j]
                    out[qi, i, j] = self._lehmann_dagger(opA, deltaA, opB, deltaB)
        return out / fx.L
