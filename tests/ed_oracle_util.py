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
        self._kernel_cache = {}   # (key_m, key_n) -> _pair_kernel result;
                                   # the kernel depends only on the sector
                                   # pair, not on which (a, c, b, d)/channel
                                   # asked for it, so every caller in a
                                   # chi_connected()/bond_correlator()/
                                   # pair_correlator() run shares one copy
                                   # per sector pair instead of rebuilding
                                   # it per cell (review finding, #151 Task
                                   # 4 fix loop: was previously rebuilt
                                   # O(nd^4) times per q on case M).
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
        so the in-sector case is not special-cased.

        Depends ONLY on the (key_m, key_n) sector pair -- not on which
        bilinear or which (q, a, c, b, d)/channel cell is asking for it --
        so it is memoized on the instance (review finding: chi_connected's
        nd**4 cells per q were each rebuilding this from scratch, most of
        them for a sector pair some earlier cell had already built)."""
        cache_key = (key_m, key_n)
        cached = self._kernel_cache.get(cache_key)
        if cached is not None:
            return cached
        ev_m, w_m = self._ev[key_m], self._w[key_m]
        ev_n, w_n = self._ev[key_n], self._w[key_n]
        ev_cat = np.concatenate([ev_m, ev_n])
        w_cat = np.concatenate([w_m, w_n])
        M = len(ev_m)
        K = _static_kernel(ev_cat, w_cat, self.fx.beta)[:M, M:]
        self._kernel_cache[cache_key] = K
        return K

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
        pin, both same-spin and cross-spin, in
        tests/test_bond_vs_ed_oracle.py).

        The (j+R) % fx.L site-offset direction and which leg carries the
        +q vs -q phase are NOT independently re-derived or cross-checked
        by anything in this module -- they are pinned by construction
        against ``chi_connected`` only at R=0 (where the offset is a
        no-op). The R != 0 orientation is deliberately left to Task 5's
        ``ed_to_solver_bond_map`` pin (full (m, m') matrix, zero coupling,
        orientation-exact against ``bond_bubble``), which is therefore the
        LOAD-BEARING check for this convention, not a redundant one."""
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
        normalized, matching chi_connected/bond_correlator.

        NOT independently numerically pinned by this module: the
        Xpp[q] = Xpp[q]^dagger check in
        tests/test_bond_vs_ed_oracle.py's pair_correlator smoke test is a
        property of ``_lehmann_dagger``'s real-symmetric-kernel structure
        (holds for ANY pair of operators sharing a sector shift, by a
        K[m,n] = K[n,m] argument -- see that test's comment), not evidence
        that the Delta operator itself (site offset, up/down leg
        assignment, annihilation-order sign) is built correctly. The
        cross-sector Lehmann sum here also uses ``mu`` through each
        sector's grand-canonical K, but every OTHER cross-sector path this
        module currently exercises (chi_connected's mixed-spin slots) has
        Delta-N = 0, where mu cancels between the two legs -- pair_correlator
        is the only in-repo path where mu enters a Delta-N != 0 Lehmann
        denominator, and nothing here checks that it is being subtracted
        correctly. Numerical correctness of this (N_up-1, N_dn-1) pair path
        -- including the mu-sensitive cross-sector denominator -- is
        PINNED BY TASK 7's pin 3b (eigenbasis_pair_bubble vs
        pair_correlator at V=0), which is therefore load-bearing and must
        NOT be weakened to a diagnostic-only comparison."""
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


# ---------------------------------------------------------------------------
# Adjudication machinery (#151, Task 6): the shared verdict rule and the
# structural-zero support propagation, used by ph (Task 6/8) and pp
# (Task 7/8) alike.
# ---------------------------------------------------------------------------

def projected_structural_zero_mask(free_support, effective_vertex_support):
    """ONE support-propagation interface for ph AND pp (round-5 fix: no
    boolean projector propagation -- ``I OR P`` cannot represent the
    triplet's FIXED-POINT CANCELLATION, where ``P[i,i] = 1`` makes
    ``Q_t[i,i] = 0`` exactly, e.g. the on-site (R=0, a=b) pair; the caller
    supplies instead the exact analytic support of the EFFECTIVE vertex).

    ALL inputs are BOOLEAN index-level supports -- never thresholded
    numeric matrices: ``free_support`` is the bubble's (ND, ND) reachability
    from the topology (all-True for these dense free-fermion fixtures,
    stated explicitly by every caller); ``effective_vertex_support`` is, for
    ph, the boolean support of the analytic ``dW`` (dS/dC); for pp channel
    eta, the boolean support of the ANALYTIC projected vertex
    ``Q_eta (dVpp) Q_eta``, constructed symbolically per channel with the
    fixed-point cancellation carried exactly.

    Propagated support = ``free @ effective_vertex @ free`` as BOOLEAN
    products (OR-AND semiring, i.e. graph reachability through two hops --
    NOT numpy's arithmetic ``@`` on bool arrays, which sums rather than
    ORs); ``zero_mask = ~propagated``. An all-False propagated support is a
    zero-only granule (every cell structurally zero).

    Both inputs must be the same (ND, ND) boolean shape; this function does
    not know or care whether ``ND`` is q-resolved -- callers that have a
    per-q effective-vertex support reduce it (structural OR over q) before
    calling, since the analytic support is a property of the topology, not
    of any single numeric q sample (a coefficient like ``cos(qR)`` is
    structurally nonzero even though it vanishes at isolated q).
    """
    free_support = np.asarray(free_support)
    effective_vertex_support = np.asarray(effective_vertex_support)
    if free_support.dtype != bool or effective_vertex_support.dtype != bool:
        # Review finding: silently casting a NUMERIC array (e.g. a raw dW
        # matrix passed by mistake instead of its boolean support) via
        # dtype=bool would treat every nonzero float as True and 0.0 as
        # False -- accidentally "working" for an exact-zero matrix while
        # defeating the whole point of the ANALYTIC (never thresholded)
        # support contract this function documents. Refuse anything that
        # is not already boolean rather than coercing it.
        raise ValueError(
            "projected_structural_zero_mask: free_support (dtype={}) and "
            "effective_vertex_support (dtype={}) must both be BOOLEAN "
            "arrays -- this function never thresholds a numeric matrix; "
            "the caller supplies the analytic support directly".format(
                free_support.dtype, effective_vertex_support.dtype))
    if free_support.ndim != 2 or free_support.shape[0] != free_support.shape[1]:
        raise ValueError(
            "projected_structural_zero_mask: free_support must be a square "
            "(ND, ND) matrix, got shape {}".format(free_support.shape))
    if free_support.shape != effective_vertex_support.shape:
        raise ValueError(
            "projected_structural_zero_mask: free_support {} and "
            "effective_vertex_support {} must share the same (ND, ND) "
            "boolean shape".format(free_support.shape,
                                    effective_vertex_support.shape))

    def _bool_matmul(a, b):
        # (a @ b)[i, j] = OR_k a[i, k] AND b[k, j]. A sum of 0/1 ANDs is > 0
        # iff at least one AND is 1, so integer matmul + threshold computes
        # the OR-AND boolean semiring product exactly.
        return (a.astype(np.int64) @ b.astype(np.int64)) > 0

    propagated = _bool_matmul(
        _bool_matmul(free_support, effective_vertex_support), free_support)
    return ~propagated


def adjudicate_granule(D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat,
                        zero_mask, label):
    """Spec tolerance rule, corrected sides and AUTHORITATIVE estimates
    (plan review round 2): the compared values are the FINER estimates,
        D_ed   = D_ed_vhalf      (smaller Richardson step)
        D_pred = D_pred_2nmat    (larger Matsubara window)
    and the coarser companions exist only to measure convergence:
        delta_rich = max|D_ed_v1 - D_ed_vhalf|        (ED side only)
        delta_nmat = max|D_pred_nmat - D_pred_2nmat|  (solver side only)
        tol = 10 * max(delta_rich, delta_nmat, 1e-13)
    zero cells (zero_mask): |D_pred| <= tol AND |D_ed| <= tol
    bearing cells: |D_pred - D_ed| <= tol
    power: max over bearing |D_pred| >= 100*tol else INCONCLUSIVE
    ZERO-ONLY granule (zero_mask.all()): all cells zero-validated, power
    test SKIPPED, status 'PASS-ZERO' (amended spec).

    ``zero_mask`` may be given at the (ND, ND) external-leg shape (the
    common case: a structural support that does not depend on q) and is
    broadcast onto ``D_pred``'s full (q, ND, ND) grid; a q-resolved mask is
    also accepted directly at the full shape.

    Returns dict(label, delta_rich, delta_nmat, tol, max_signal, status,
    failures, pred_full, bearing_mask) -- the CANONICAL FULL-GRID record
    (plan review round 3: a bearing-cells-only flattened column cannot be
    aligned across directions). pred_full is D_pred on the canonical cell
    grid, flattened in the FIXED order (q ascending, external leg I
    row-major, external leg J row-major, then real part followed by
    imaginary part of each cell -- one interleaved (..., 2) axis
    flattened last); bearing_mask is ~zero_mask DUPLICATED onto that
    size-2 component axis before the identical flattening (round-4 fix:
    without the duplication its length is half of pred_full's).
    sensitivity_rank aligns columns by this shared grid with no
    remapping. The record IS the Task-6/7/8 -> Task-9 interface. Prints
    the audit tuple. Pure function; no global state.
    """
    D_ed_v1 = np.asarray(D_ed_v1)
    D_ed_vhalf = np.asarray(D_ed_vhalf)
    D_pred_nmat = np.asarray(D_pred_nmat)
    D_pred_2nmat = np.asarray(D_pred_2nmat)
    zero_mask = np.asarray(zero_mask, dtype=bool)

    if D_ed_v1.shape != D_ed_vhalf.shape:
        raise ValueError(
            "adjudicate_granule[{}]: D_ed_v1 shape {} != D_ed_vhalf shape "
            "{}".format(label, D_ed_v1.shape, D_ed_vhalf.shape))
    if D_pred_nmat.shape != D_pred_2nmat.shape:
        raise ValueError(
            "adjudicate_granule[{}]: D_pred_nmat shape {} != D_pred_2nmat "
            "shape {}".format(label, D_pred_nmat.shape, D_pred_2nmat.shape))
    if D_ed_vhalf.shape != D_pred_2nmat.shape:
        raise ValueError(
            "adjudicate_granule[{}]: ED grid shape {} != solver grid shape "
            "{}".format(label, D_ed_vhalf.shape, D_pred_2nmat.shape))

    # Review finding (must_fix): a NaN/Inf anywhere compares False against
    # every numeric threshold below (|x| > tol, |x| <= tol alike), so a
    # broken (non-finite) upstream computation would silently slip through
    # every zero/bearing check and could even land a zero-only granule on
    # PASS-ZERO. Fail loudly instead, exactly as Task 5's Pin 2b does for
    # the same class of hazard.
    for name, arr in (("D_ed_v1", D_ed_v1), ("D_ed_vhalf", D_ed_vhalf),
                       ("D_pred_nmat", D_pred_nmat),
                       ("D_pred_2nmat", D_pred_2nmat)):
        if not np.all(np.isfinite(arr)):
            raise ValueError(
                "adjudicate_granule[{}]: {} contains a non-finite (NaN/Inf) "
                "value -- refusing to adjudicate a broken input".format(
                    label, name))

    # AUTHORITATIVE sides (plan review round 2).
    D_ed = D_ed_vhalf
    D_pred = D_pred_2nmat

    delta_rich = float(np.max(np.abs(D_ed_v1 - D_ed_vhalf)))
    delta_nmat = float(np.max(np.abs(D_pred_nmat - D_pred_2nmat)))
    tol = 10.0 * max(delta_rich, delta_nmat, 1e-13)

    if zero_mask.shape == D_pred.shape:
        zero_full = zero_mask
    else:
        zero_full = np.broadcast_to(zero_mask, D_pred.shape)
    bearing_full = ~zero_full

    diff = np.abs(D_pred - D_ed)
    zero_bad = zero_full & ((np.abs(D_pred) > tol) | (np.abs(D_ed) > tol))
    bearing_bad = bearing_full & (diff > tol)

    failures = []
    for idx in np.argwhere(zero_bad):
        idx = tuple(int(i) for i in idx)
        failures.append(("zero", idx, complex(D_pred[idx]), complex(D_ed[idx]),
                          float(diff[idx])))
    for idx in np.argwhere(bearing_bad):
        idx = tuple(int(i) for i in idx)
        failures.append(("bearing", idx, complex(D_pred[idx]), complex(D_ed[idx]),
                          float(diff[idx])))

    zero_only = bool(zero_full.all())
    has_bearing_signal = bool(bearing_full.any())
    max_signal = float(np.max(np.abs(D_pred[bearing_full]))) \
        if has_bearing_signal else 0.0

    if failures:
        status = "FAIL"
    elif zero_only:
        status = "PASS-ZERO"
        max_signal = 0.0   # power test skipped for a zero-only granule
    elif max_signal >= 100.0 * tol:
        status = "PASS"
    else:
        status = "INCONCLUSIVE"

    pred_full = np.stack([D_pred.real, D_pred.imag], axis=-1).ravel().copy()
    bearing_mask = np.stack([bearing_full, bearing_full], axis=-1).ravel().copy()

    print("adjudicate[{}]: delta_rich={:.6e} delta_nmat={:.6e} tol={:.6e} "
          "max_signal={:.6e} status={} n_failures={}".format(
              label, delta_rich, delta_nmat, tol, max_signal, status,
              len(failures)))

    return dict(label=label, delta_rich=delta_rich, delta_nmat=delta_nmat,
                tol=tol, max_signal=max_signal, status=status,
                failures=failures, pred_full=pred_full,
                bearing_mask=bearing_mask)


# ---------------------------------------------------------------------------
# pp reference path (#151, Task 7): the exchange projector (fermionic
# exchange symmetry, independent of bond_channels) and the independent
# finite-Nmat pp bubble reference (direct sum over the one-body eigenbasis).
# ---------------------------------------------------------------------------

def exchange_projector(channels):
    """``(P, Q_s, Q_t)`` over the ordered-pair basis ``channels`` (a list of
    ``(R, a, b)`` tuples -- ALWAYS the full triple, never the ``(R,)``
    norb=1 shorthand ``SectorED.bond_correlator``/``pair_correlator`` accept,
    since orbital indices are meaningful for a pair operator even at
    norb=1), built from FERMIONIC EXCHANGE SYMMETRY ALONE and independent of
    ``hwave.solver.bond_channels`` (design doc, "Pair operators are
    projected before comparison": the singlet/triplet resolution of the
    ordered-spin pair ``Delta_{R,a,b}(q) = sum_j e^{iqj} c_{j+R,a,down}
    c_{j,b,up}`` is fixed by requiring the OVERALL two-electron amplitude be
    symmetric (singlet) / antisymmetric (triplet) under simultaneous
    ``R -> -R``, orbital swap ``a<->b`` -- the two ways of describing the
    same bond from either endpoint).

    ``P[i, j] = 1`` iff ``channels[j] == (-R, b, a)`` for
    ``channels[i] == (R, a, b)``; ``Q_s = (I + P)/2``, ``Q_t = (I - P)/2``.
    ``channels`` must be UNIQUE and CLOSED under the involution
    ``(R, a, b) -> (-R, b, a)`` (round-4 fix: the invariant was implicit) --
    raises ``ValueError`` naming the offending/missing channel otherwise.
    ``P`` is idempotent-involutive (``P @ P = I``), so ``Q_s``/``Q_t`` are
    complementary orthogonal projectors (``Q_s + Q_t = I``, ``Q_s^2 = Q_s``,
    ``Q_t^2 = Q_t``, ``Q_s Q_t = 0``) for any closed, unique channel list.
    """
    channels = list(channels)
    index_of = {}
    for i, ch in enumerate(channels):
        if ch in index_of:
            raise ValueError(
                "exchange_projector: channels must be unique -- {} appears "
                "at both index {} and index {}".format(
                    ch, index_of[ch], i))
        index_of[ch] = i
    n = len(channels)
    P = np.zeros((n, n), dtype=complex)
    for i, ch in enumerate(channels):
        R, a, b = ch
        partner = (-R, b, a)
        j = index_of.get(partner)
        if j is None:
            raise ValueError(
                "exchange_projector: channels is not closed under the "
                "involution (R,a,b)->(-R,b,a) -- channel {} is missing its "
                "partner {}".format(ch, partner))
        P[i, j] = 1.0
    Id = np.eye(n, dtype=complex)
    Q_s = 0.5 * (Id + P)
    Q_t = 0.5 * (Id - P)
    return P, Q_s, Q_t


def eigenbasis_pair_bubble(fx, channels, nmat):
    """The INDEPENDENT finite-Nmat pp bubble reference: a direct Matsubara
    sum over the one-body EIGENBASIS (never through momentum space or
    ``hwave.solver.bond_channels._g2_from_green`` -- a genuinely different
    computational route from ``production_pair_bubble``, compared at pin 3a).

    Derivation (Wick's theorem on the free reference state, verified
    numerically against ``SectorED.pair_correlator`` at pin 3b): for
    ``Delta_{R,a,b}(q=0) = sum_j c_{j+R,a,down} c_{j,b,up}``, the only
    nonzero contraction of ``<Delta_{R,a,b}; Delta_{R',c,d}^dagger>`` pairs
    the two DOWN legs and the two UP legs (down/up never contract in a
    spin-conserving free theory); converting the resulting imaginary-time
    integral to a Matsubara sum gives, per site pair ``(j, j')``::

        X0_pp[(R,a,b),(R',c,d)]
            = (1/(L*beta)) sum_n sum_{j,j'}
                  G[(j+R,a),(j'+R',c)](iwn) * G[(j,b),(j',d)](-iwn)

    with ``G`` the free single-spin-sector Green function (spin-independent
    ``H1``, so up and down share one ``G``), built here directly from the
    one-body eigenbasis: ``G(iwn) = U @ diag(1/(i*wn + mu - eps)) @
    U^dagger`` for ``eps, U = eigh(H1_one_spin_sector)``. ``G(-iwn[n])`` is
    read off the SAME array at the mirrored index ``nmat-1-n`` (the centered
    fermionic Matsubara grid satisfies ``iwn[nmat-1-n] = -iwn[n]`` exactly,
    an algebraic identity, not an approximation).

    ``channels``: a list of ``(R, a, b)`` triples (the full ordered-pair
    basis, matching ``exchange_projector``'s convention -- NOT
    ``SectorED``'s ``(R,)`` norb=1 shorthand; the caller maps between the two
    when calling ``SectorED.pair_correlator``).
    """
    # Review finding (must_fix): the centered fermionic grid
    # iwn = (2n+1-nmat)*pi/beta is only antisymmetric (iwn[nmat-1-n] =
    # -iwn[n], the identity this function's G(-iwn) shortcut relies on) for
    # an EVEN nmat; an odd/zero/negative nmat silently builds a grid that is
    # not a fermionic Matsubara grid at all (e.g. containing iwn=0), and
    # production_pair_bubble (same formula, via free_green) would silently
    # agree with it at pin 3a while both compute the wrong quantity.
    # Round-2 review finding (should_fix): int(nmat) on a NaN/inf/None/
    # non-numeric value raises ValueError/TypeError/OverflowError instead of
    # the intended, consistent ValueError -- normalize through float() first
    # and reject non-finite values explicitly before the integer/parity check.
    try:
        nmat_f = float(nmat)
    except (TypeError, ValueError):
        raise ValueError(
            "eigenbasis_pair_bubble: nmat must be a positive even integer "
            "(the centered fermionic Matsubara grid) -- got {!r}".format(nmat))
    if (not np.isfinite(nmat_f) or nmat_f != int(nmat_f) or nmat_f <= 0
            or int(nmat_f) % 2 != 0):
        raise ValueError(
            "eigenbasis_pair_bubble: nmat must be a positive even integer "
            "(the centered fermionic Matsubara grid) -- got {!r}".format(nmat))
    nmat = int(nmat_f)
    L, norb = fx.L, fx.norb
    idx_up = [fx.mode(j, o, 0) for j in range(L) for o in range(norb)]
    h1_spin = fx.build_h1()[np.ix_(idx_up, idx_up)]
    eps, U = np.linalg.eigh(h1_spin)
    beta = fx.beta
    iwn = (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi / beta
    denom = 1j * iwn[:, None] + fx.mu - eps[None, :]
    ginv = 1.0 / denom
    green = np.einsum('pa,na,qa->npq', U, ginv, np.conj(U))   # (nmat, M, M)
    green_minus = green[::-1]                                  # G(-iwn[n])

    def _mode2(j, o):
        return j * norb + o

    channels = list(channels)
    n_ch = len(channels)
    X = np.zeros((n_ch, n_ch), dtype=complex)
    for i, (R, a, b) in enumerate(channels):
        for jx, (Rp, c, d) in enumerate(channels):
            acc = 0.0 + 0.0j
            for j in range(L):
                p1 = _mode2((j + R) % L, a)
                p2 = _mode2(j, b)
                for jp in range(L):
                    q1 = _mode2((jp + Rp) % L, c)
                    q2 = _mode2(jp, d)
                    acc += np.sum(green[:, p1, q1] * green_minus[:, p2, q2])
            X[i, jx] = acc / (L * beta)
    return X
