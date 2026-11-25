"""Pattern realisation, phases and signs (``hwave.initial.patterns``;
design sections 5.1 "Rows per eigenvector" and 6.3; validation R2).

``realise`` turns a PatternSpec into the Hamiltonian-convention ``O(k)``
of a folded cell. A bond slot ``(a, b, r)`` multiplies ``c^dag_a(R + r)
c_b(R)`` (R = the cell of the annihilation operator); the field is the
anchored modulation ``O(R) = (e^{i(q.R+phi)} A + e^{-i(q.R+phi)} A^dag)/2``
with ``(A^dag)_ab(r) = e^{-iq.r} conj(A_ba(-r))``. The oracle below builds
the same operator on the whole lattice with explicit loops, and reads the
Hamiltonian-convention blocks off it, so it shares no code with
``realise``.
"""
import cmath
import math
import os
import shutil
import tempfile
import time
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures

_S0 = np.eye(2, dtype=complex)
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
_SZ = np.array([[1, 0], [0, -1]], dtype=complex)
_SIG = {"c0": _S0, "cx": _SX, "cy": _SY, "cz": _SZ}

_CHAIN = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0}


def _qdot(q, R):
    """``q.R`` modulo 1, exactly."""
    return float(sum(Fraction(x) * int(r) for x, r in zip(q, R)) % 1)


def _block(entry):
    """The 2x2 spin block of a 'cos' entry, ``1/2 sum_mu c^mu sigma^mu``."""
    return 0.5 * sum(getattr(entry, k) * _SIG[k] for k in _SIG)


def _oracle(spec, facts, subshape, *, normalise=True):
    """The Hamiltonian-convention ``O(k)`` of a 'cos' or 'spiral' spec on
    the folded cell ``subshape``, built on the whole lattice.

    ``M[(cell_i, a, s), (cell_j, b, t)]`` is the coefficient of
    ``c^dag_{a s}(cell_i) c_{b t}(cell_j)``; the folded block is
    ``O_H(D)[aa, bb] = M[(ra, a, s), (D*B + rb, b, t)]`` (``aa`` in
    supercell 0)."""
    p = facts.problem
    L = tuple(int(x) for x in p.cellshape)
    B = tuple(int(x) for x in subshape)
    nphys = facts.norb_phys
    so = facts.spin_orbital
    cells = [(x, y, z) for x in range(L[0]) for y in range(L[1]) for z in range(L[2])]
    cidx = {c: i for i, c in enumerate(cells)}
    n = len(cells) * nphys * 2

    def g(cell, a, s):
        cell = tuple(c % l for c, l in zip(cell, L))
        return (cidx[cell] * nphys + a) * 2 + s

    M = np.zeros((n, n), dtype=complex)
    q = spec.q
    terms = []        # (blocks per entry, phase)
    if spec.modulation == "cos":
        terms.append(({e.key: _block(e) for e in spec.entries}, spec.phase))
    else:
        sig = {"x": _SX, "y": _SY, "z": _SZ}
        for ax, shift in zip(spec.axes, (0.0, -math.pi / 2)):
            terms.append(({e.key: e.s * sig[ax] for e in spec.entries},
                          spec.phase + shift))
    for amp, phi in terms:
        for R in cells:
            ph = cmath.exp(1j * (2 * math.pi * _qdot(q, R) + phi))
            for (a, b, r), X in amp.items():
                Rc = tuple(R[i] + r[i] for i in range(3))
                # A: c^dag_a(R + r) c_b(R)
                for s in range(2):
                    for t in range(2):
                        M[g(Rc, a, s), g(R, b, t)] += 0.5 * ph * X[s, t]
                # A^dag: (A^dag)_ba(-r) = e^{iq.r} conj(A_ab(r)), slot
                # c^dag_b(R - r) c_a(R)
                Rm = tuple(R[i] - r[i] for i in range(3))
                adj = cmath.exp(2j * math.pi * _qdot(q, r)) * np.conj(X).T
                for t in range(2):
                    for s in range(2):
                        M[g(Rm, b, t), g(R, a, s)] += 0.5 * np.conj(ph) * adj[t, s]
    M *= spec.sign
    assert np.allclose(M, M.conj().T, atol=1e-13)
    # fold: supercell D, in-cell offsets
    N = tuple(l // b for l, b in zip(L, B))
    bvol = B[0] * B[1] * B[2]
    norb_b = nphys * bvol * (2 if so else 1)
    nd = 2 * nphys * bvol
    offs = [(x, y, z) for z in range(B[2]) for y in range(B[1]) for x in range(B[0])]

    def fidx(ir, a, s):
        if so:
            return (2 * a + s) + 2 * nphys * ir
        return s * (nphys * bvol) + a + nphys * ir

    OH = np.zeros(N + (nd, nd), dtype=complex)
    for D in [(x, y, z) for x in range(N[0]) for y in range(N[1]) for z in range(N[2])]:
        for ia, ra in enumerate(offs):
            for ib, rb in enumerate(offs):
                cj = tuple(D[i] * B[i] + rb[i] for i in range(3))
                for a in range(nphys):
                    for b in range(nphys):
                        for s in range(2):
                            for t in range(2):
                                OH[D][fidx(ia, a, s), fidx(ib, b, t)] = \
                                    M[g(ra, a, s), g(cj, b, t)]
    del norb_b
    ok = np.fft.ifftn(OH, axes=(0, 1, 2), norm="forward").reshape(-1, nd, nd)
    if normalise:
        ok = ok / math.sqrt(np.sum(np.abs(OH) ** 2) / bvol)
    return ok


def _to_r(o_k, shape):
    nd = o_k.shape[-1]
    return np.fft.fftn(o_k.reshape(tuple(shape) + (nd, nd)), axes=(0, 1, 2),
                       norm="forward")


def _site_values(o_k, facts, subshape):
    """The on-site, spin-up diagonal of ``O(R = 0)`` per in-cell site (one
    orbital, normal mode)."""
    shape = tuple(c // b for c, b in zip(facts.cellshape, subshape))
    o_r = _to_r(o_k, shape)[0, 0, 0]
    bvol = int(np.prod(subshape))
    return np.array([o_r[i, i].real for i in range(bvol)])


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_phases_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def chain(self, L, subshape, *, inter=None, ising=None, U=1.0, norb=1,
              transfer=None, centers=None, so=False, ncond=None):
        from hwave.initial import facts
        nphys = norb // 2 if so else norb
        interactions = {"CoulombIntra": {((0, 0, 0), a, a): U for a in range(nphys)}}
        if inter:
            interactions["CoulombInter"] = {((1, 0, 0), 0, 0): inter,
                                            ((-1, 0, 0), 0, 0): inter}
        if ising:
            interactions["Ising"] = {((1, 0, 0), 0, 0): ising, ((-1, 0, 0), 0, 0): ising}
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(L, 1, 1), norb=norb,
            transfer=_CHAIN if transfer is None else transfer, centers=centers,
            interactions=interactions,
            param={"Ncond": L * norb if ncond is None else ncond, "T": 0.1,
                   "SubShape": list(subshape)},
            mode={"enable_spin_orbital": True} if so else None)
        return facts.model_facts(toml)


def _onsite(c0=0, cx=0, cy=0, cz=0, a=0):
    from hwave.initial.patterns import Entry
    return Entry(a, a, (0, 0, 0), c0=c0, cx=cx, cy=cy, cz=cz)


def _spec(entries, q, subshape, *, phase=0.0, sign=1, modulation="cos", axes=None):
    from hwave.initial.patterns import PatternSpec
    return PatternSpec(q=tuple(str(Fraction(x)) for x in q), subshape=subshape, phase=phase,
                       sign=sign, modulation=modulation, axes=axes, entries=tuple(entries))


# --------------------------------------------------------------------------
# realise against the whole-lattice oracle
# --------------------------------------------------------------------------

class TestRealise(_Tmp):
    def test_matches_oracle_onsite_and_bond(self):
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        f = self.chain(4, (1, 1, 1))
        spec = _spec([Entry(0, 0, (1, 0, 0), c0=0.3 + 0.7j, cy=0.2 - 0.1j),
                      _onsite(c0=0.4, cz=-0.9)], ("1/4", "0", "0"), (4, 1, 1),
                     phase=0.37, sign=-1)
        got = patterns.realise(spec, f, (4, 1, 1))
        np.testing.assert_allclose(got, _oracle(spec, f, (4, 1, 1)), atol=1e-13)

    def test_normalised_per_primitive_cell(self):
        from hwave.initial import patterns
        f = self.chain(4, (1, 1, 1))
        for B in ((2, 1, 1), (4, 1, 1)):
            o = patterns.realise(_spec([_onsite(cz=1.0)], ("1/2", "0", "0"), B), f, B)
            bvol = int(np.prod(B))
            self.assertAlmostEqual(float(np.sum(np.abs(o) ** 2)) / o.shape[0] / bvol,
                                   1.0, places=12)

    def test_transfer_as_pattern(self):
        """The Hamiltonian convention is the solver's: a Hermitian Transfer
        written as a q = 0 pattern (a row ``(R, a, b): v`` is the slot
        ``(a, b, -R)``, since a slot follows the Green layout ``G[r] =
        <c^dag(r) c(0)>``) realises to ``ham_trans_k`` itself, on the
        primitive and on a folded cell."""
        from hwave.initial import facts, patterns
        from hwave.initial.patterns import Entry
        t = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
             ((1, 0, 0), 0, 1): -0.3 + 0.2j, ((-1, 0, 0), 1, 0): -0.3 - 0.2j,
             ((0, 0, 0), 0, 1): 0.1j, ((0, 0, 0), 1, 0): -0.1j, ((0, 0, 0), 1, 1): 0.4}
        for B in ((1, 1, 1), (2, 1, 1)):
            toml = initial_fixtures.write_input(
                self.dest(), cellshape=(4, 1, 1), norb=2, transfer=t,
                interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 1.0}},
                param={"Ncond": 4, "T": 0.1, "SubShape": list(B)})
            f = facts.model_facts(toml)
            A = [Entry(a, b, tuple(-x for x in R), c0=2 * v) for (R, a, b), v in t.items()]
            got = patterns.realise(_spec(A, (0, 0, 0), (1, 1, 1)), f, B, normalise=False)
            np.testing.assert_allclose(got, f.problem.ham_trans_k, atol=1e-14)

    def test_q_zero_equals_fold_layout(self):
        """At q = 0 the realised pattern is the primitive one folded by the
        solver's shared Hamiltonian-convention fold."""
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        from hwave.solver import fold
        f = self.chain(4, (1, 1, 1))
        spec = _spec([Entry(0, 0, (1, 0, 0), c0=0.5 + 0.2j, cx=0.3),
                      _onsite(cz=0.7)], ("0", "0", "0"), (1, 1, 1), phase=0.2)
        prim = patterns.realise(spec, f, (1, 1, 1))
        o_r = _to_r(prim, (4, 1, 1)).reshape(4, 2, 1, 2, 1)
        folded = fold.fold_layout(o_r, cellshape=(4, 1, 1), subshape=(2, 1, 1),
                                  norb_orig=1, ns=2, hamiltonian=True)
        want = np.fft.ifftn(folded.reshape(2, 1, 1, 4, 4), axes=(0, 1, 2),
                            norm="forward").reshape(2, 4, 4)
        got = patterns.realise(spec, f, (2, 1, 1))
        np.testing.assert_allclose(got, want, atol=1e-13)

    def test_hermitian(self):
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        f = self.chain(6, (1, 1, 1))
        spec = _spec([Entry(0, 0, (1, 0, 0), c0=0.3 + 0.7j, cx=1j)], ("1/3", "0", "0"),
                     (3, 1, 1), phase=0.1)
        o = patterns.realise(spec, f, (3, 1, 1))
        np.testing.assert_allclose(o, np.conj(np.swapaxes(o, -1, -2)), atol=1e-14)

    def test_anchored_adjoint_q_bond(self):
        """The R0 bond-sign pin at the pattern level: a complex bond
        amplitude at q = (1/2, 0, 0). The A^dag slot carries
        ``e^{-iq.r} conj(A_ba(-r))``: at q = 1/2 and r = 1 this is a sign,
        so the wrong rule (no phase, or e^{+iq.r} conjugated) gives a
        different operator."""
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        f = self.chain(4, (1, 1, 1))
        spec = _spec([Entry(0, 0, (1, 0, 0), c0=0.6 + 0.8j)], ("1/2", "0", "0"), (2, 1, 1),
                     phase=0.3)
        got = patterns.realise(spec, f, (2, 1, 1), normalise=False)
        np.testing.assert_allclose(got, _oracle(spec, f, (2, 1, 1), normalise=False),
                                   atol=1e-14)
        # the bond c^dag(R+1) c(R) for R = 0 (in-cell site 1 <- site 0):
        # (e^{i phi} A)/2 from A, plus the A^dag slot of R' = 1 on the same
        # bond read backwards (c^dag(0) c(1)) gives e^{-i(pi+phi)} e^{i pi}...
        o_r = _to_r(got, (2, 1, 1))
        a = 0.6 + 0.8j
        # O_H(D=0)[site1 up, site0 up] = coefficient of c^dag_1 c_0, i.e.
        # slot (R=0, r=+1): (e^{i phi} A)/2 with the spin block A = a/2
        self.assertAlmostEqual(complex(o_r[0, 0, 0][1, 0]), cmath.exp(0.3j) * a / 4, places=13)
        # O_H(D=0)[site0 up, site1 up] = coefficient of c^dag_0 c_1: the
        # A^dag slot with annihilation cell R = 1, r = -1:
        # e^{-i(pi + phi)} e^{-iq.(-1)} conj(a/2) / 2 = e^{-i phi} conj(a) / 4
        self.assertAlmostEqual(complex(o_r[0, 0, 0][0, 1]),
                               cmath.exp(-0.3j) * np.conj(a) / 4, places=13)
        # the wrong rule (no e^{-iq.r}) would give the opposite sign here
        self.assertGreater(abs(complex(o_r[0, 0, 0][0, 1])
                               + cmath.exp(-0.3j) * np.conj(a) / 4), 0.1)

    def test_explicit_minus_r_no_doubling(self):
        """A user pattern given with both +r and -r entries is two entries
        of A: the realised O is the sum of the two, nothing is added or
        merged."""
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        f = self.chain(4, (1, 1, 1))
        e1 = Entry(0, 0, (1, 0, 0), c0=0.5 + 0.1j)
        e2 = Entry(0, 0, (-1, 0, 0), c0=0.2 - 0.4j)
        q = ("1/4", "0", "0")
        pair = patterns.realise(_spec([e1, e2], q, (4, 1, 1)), f, (4, 1, 1), normalise=False)
        one = patterns.realise(_spec([e1], q, (4, 1, 1)), f, (4, 1, 1), normalise=False)
        two = patterns.realise(_spec([e2], q, (4, 1, 1)), f, (4, 1, 1), normalise=False)
        np.testing.assert_allclose(pair, one + two, atol=1e-14)
        np.testing.assert_allclose(pair, _oracle(_spec([e1, e2], q, (4, 1, 1)), f, (4, 1, 1),
                                                 normalise=False), atol=1e-14)

    def test_spin_orbital_physical_index(self):
        """In spin-orbital mode an entry on physical orbital 1 lands on the
        interleaved slots 2 and 3 (spin_orbital.interleaved_index)."""
        from hwave.initial import patterns
        so_t = {(r, 2 * a + s, 2 * a + s): v for (r, a, _b), v in
                {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
                 ((1, 0, 0), 1, 1): -1.0, ((-1, 0, 0), 1, 1): -1.0}.items()
                for s in range(2)}
        f = self.chain(4, (1, 1, 1), norb=4, transfer=so_t, so=True, ncond=8)
        o = patterns.realise(_spec([_onsite(cz=1.0, a=1)], ("0", "0", "0"), (1, 1, 1)), f,
                             (1, 1, 1))
        o_r = _to_r(o, (4, 1, 1))[0, 0, 0]
        d = np.real(np.diag(o_r))
        np.testing.assert_allclose(d, [0, 0, 1 / math.sqrt(2), -1 / math.sqrt(2)], atol=1e-14)
        np.testing.assert_allclose(o_r - np.diag(np.diag(o_r)), 0, atol=1e-14)

    def test_random_same_seed_same_pattern(self):
        from hwave.initial import patterns
        from hwave.initial.patterns import PatternSpec
        f = self.chain(4, (2, 1, 1))
        spec = PatternSpec(q=None, subshape=(2, 1, 1), phase=None, sign=None,
                           modulation="random", seed=7)
        o1 = patterns.realise(spec, f, (2, 1, 1))
        o2 = patterns.realise(spec, f, (2, 1, 1))
        np.testing.assert_array_equal(o1, o2)
        o3 = patterns.realise(PatternSpec(q=None, subshape=(2, 1, 1), phase=None, sign=None,
                                          modulation="random", seed=8), f, (2, 1, 1))
        self.assertGreater(float(np.max(np.abs(o1 - o3))), 1e-3)
        self.assertAlmostEqual(float(np.sum(np.abs(o1) ** 2)) / o1.shape[0] / 2, 1.0, places=12)
        # cell-resolved: the two cells of B differ
        o_r = _to_r(o1, (2, 1, 1))[0, 0, 0]
        self.assertGreater(abs(o_r[0, 0] - o_r[1, 1]), 1e-3)
        np.testing.assert_allclose(o1, np.conj(np.swapaxes(o1, -1, -2)), atol=1e-14)

    def test_random_interacting_orbitals_only(self):
        """A random pattern lives on the interacting orbitals I_c (a ligand
        with no two-body term stays zero) and within the solver blocks."""
        from hwave.initial import patterns
        from hwave.initial.patterns import PatternSpec
        f = self.chain_with_ligand()
        self.assertEqual(f.interacting_orbitals, (0,))
        o = patterns.realise(PatternSpec(q=None, subshape=(1, 1, 1), phase=None, sign=None,
                                         modulation="random", seed=3), f, (1, 1, 1))
        o_r = _to_r(o, f.problem.shape)
        norb = f.problem.norb
        lig = [1, norb + 1]
        self.assertEqual(float(np.max(np.abs(o_r[..., lig, :]))), 0.0)
        self.assertEqual(float(np.max(np.abs(o_r[..., :, lig]))), 0.0)

    def chain_with_ligand(self):
        from hwave.initial import facts
        t = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
             ((0, 0, 0), 0, 1): -0.5, ((0, 0, 0), 1, 0): -0.5}
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=2, transfer=t,
            centers=[(0, 0, 0), (0.5, 0, 0)],
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0}},
            param={"Ncond": 4, "T": 0.1, "SubShape": [1, 1, 1]})
        return facts.model_facts(toml)

    def test_random_within_blocks(self):
        """Two interacting orbitals that nothing connects form separate
        blocks: the random pattern has no entry between them."""
        from hwave.initial import facts, patterns
        from hwave.initial.patterns import PatternSpec
        t = {((1, 0, 0), a, a): -1.0 for a in range(2)}
        t.update({((-1, 0, 0), a, a): -1.0 for a in range(2)})
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=2, transfer=t,
            interactions={"CoulombIntra": {((0, 0, 0), a, a): 2.0 for a in range(2)}},
            param={"Ncond": 4, "T": 0.1, "SubShape": [1, 1, 1]})
        f = facts.model_facts(toml)
        self.assertGreater(len(f.problem.blocks.block_info), 1)
        o = patterns.realise(PatternSpec(q=None, subshape=(1, 1, 1), phase=None, sign=None,
                                         modulation="random", seed=2), f, (1, 1, 1))
        block_of = np.empty(f.problem.nd, dtype=int)
        for k, blk in enumerate(f.problem.blocks.block_info):
            block_of[np.asarray(blk)] = k
        cross = block_of[:, None] != block_of[None, :]
        self.assertEqual(float(np.max(np.abs(o[:, cross]))), 0.0)
        self.assertGreater(float(np.max(np.abs(o))), 0.1)

    def test_random_2sz_spin_diagonal(self):
        from hwave.initial import facts, patterns
        from hwave.initial.patterns import PatternSpec
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), transfer=_CHAIN,
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0}},
            param={"Ncond": 4, "T": 0.1, "SubShape": [2, 1, 1], "2Sz": 0})
        f = facts.model_facts(toml)
        o = patterns.realise(PatternSpec(q=None, subshape=(2, 1, 1), phase=None, sign=None,
                                         modulation="random", seed=1), f, (2, 1, 1))
        n = f.problem.norb
        self.assertEqual(float(np.max(np.abs(o[:, :n, n:]))), 0.0)
        self.assertEqual(float(np.max(np.abs(o[:, n:, :n]))), 0.0)

    def test_not_representable_refused(self):
        from hwave.initial import patterns
        f = self.chain(4, (1, 1, 1))
        with self.assertRaisesRegex(ValueError, "not representable"):
            patterns.realise(_spec([_onsite(c0=1)], ("1/3", "0", "0"), (3, 1, 1)), f,
                             (2, 1, 1))


# --------------------------------------------------------------------------
# phi0, phase rows, signs
# --------------------------------------------------------------------------

class TestPhases(_Tmp):
    def test_phi0_rule(self):
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        # on-site: tr AA = e^{2 i theta} tr(X X)
        for theta in (0.0, 0.3, -1.2, 2.0):
            A = (_onsite(c0=cmath.exp(1j * theta)),)
            phi = patterns.phi0(A, ("0", "0", "0"))
            self.assertAlmostEqual(cmath.exp(2j * (phi + theta)).real, 1.0, places=12)
        # a lone bond slot without its -r partner: tr AA = 0 -> phi0 = 0
        A = (Entry(0, 0, (1, 0, 0), c0=0.3 + 0.4j),)
        self.assertEqual(abs(patterns.trAA(A, ("0", "0", "0"))), 0.0)
        self.assertEqual(patterns.phi0(A, ("1/4", "0", "0")), 0.0)
        # sigma^x versus sigma^y weight: tr AA of (sx + i sy)/2 vanishes
        A = (_onsite(cx=1.0, cy=1j),)
        self.assertLess(abs(patterns.trAA(A, ("0", "0", "0"))), 1e-15)
        self.assertEqual(patterns.phi0(A, ("0", "0", "0")), 0.0)
        # the bilinear is the explicit bond-slot sum, with e^{iq.r}
        A = (Entry(0, 0, (1, 0, 0), c0=0.6), Entry(0, 0, (-1, 0, 0), c0=0.2j))
        q = ("1/3", "0", "0")
        x, y = 0.3, 0.1j        # the spin blocks are c0/2 * 1
        want = 2 * (x * y * cmath.exp(2j * math.pi / 3) + y * x * cmath.exp(-2j * math.pi / 3))
        self.assertGreater(abs(want), 0.01)
        self.assertAlmostEqual(patterns.trAA(A, q), want, places=14)

    def test_phi0_maximises_the_hermitian_norm(self):
        """At q = -q the per-cell norm of O(phi) is that of the R = 0 cell,
        and phi0 maximises it (the anchored adjoint is the one realise
        uses)."""
        from hwave.initial import patterns
        from hwave.initial.patterns import Entry
        f = self.chain(4, (1, 1, 1))
        A = (Entry(0, 0, (1, 0, 0), c0=0.6 + 0.3j), Entry(0, 0, (-1, 0, 0), c0=0.1 - 0.5j),
             _onsite(c0=0.2 + 0.4j))
        q = ("1/2", "0", "0")
        phi0 = patterns.phi0(A, q)
        grid = np.linspace(-math.pi / 2, math.pi / 2, 721)
        norms = [float(np.sum(np.abs(patterns.realise(_spec(A, q, (2, 1, 1), phase=p), f,
                                                      (2, 1, 1), normalise=False)) ** 2))
                 for p in grid]
        best = grid[int(np.argmax(norms))]
        self.assertLess(abs(cmath.exp(2j * best) - cmath.exp(2j * phi0)), 0.01)

    def test_two_phases_q_quarter(self):
        from hwave.initial import patterns
        f = self.chain(4, (4, 1, 1))
        A = (_onsite(c0=1.0),)
        q = ("1/4", "0", "0")
        rows = patterns.phase_rows(A, q, sign_map=False)
        self.assertEqual(len(rows), 2)
        self.assertEqual([r.neg for r in rows], [False, False])
        vals = [_site_values(patterns.realise(_spec(A, q, (4, 1, 1), phase=r.phase), f,
                                              (4, 1, 1)), f, (4, 1, 1)) for r in rows]
        scale = float(np.max(np.abs(vals[0])))
        np.testing.assert_allclose(vals[0] / scale, [1, 0, -1, 0], atol=0.01)
        np.testing.assert_allclose(vals[1] / scale, [.71, -.71, -.71, .71], atol=0.01)
        self.assertIsNone(rows[1].equivalent_to)
        # the same rows with a sign map (even n: -O is the n/2 translate)
        self.assertEqual(patterns.phase_rows(A, q, sign_map=True), rows)

    def test_q_third_tr_even_no_map(self):
        """A one-orbital chain with CoulombInter V, a charge pattern at
        q = 1/3 (n = 3, TR-even, no sign map)."""
        from hwave.initial import patterns, symmetry
        f = self.chain(6, (3, 1, 1), inter=0.5)
        A = (_onsite(c0=1.0),)
        q = ("1/3", "0", "0")
        O = patterns.realise(_spec(A, q, (3, 1, 1)), f, (3, 1, 1))
        su2 = symmetry.su2_check(f.problem)
        tr = symmetry.tr_check(f.problem)
        sign_map = symmetry.sign_map_verified(O, problem=f.problem, su2=su2, tr=tr)
        self.assertFalse(sign_map)
        A = patterns.fix_sign(A, q)
        rows = patterns.phase_rows(A, q, sign_map=sign_map)
        self.assertEqual(len(rows), 2)
        phi0 = patterns.phi0(A, q)
        self.assertAlmostEqual(rows[0].phase, phi0, places=15)
        self.assertAlmostEqual(rows[1].phase, phi0 + math.pi / 3, places=15)
        self.assertIsNone(rows[0].equivalent_to)
        self.assertEqual(rows[1].equivalent_to, "-O(phi0), translated")
        self.assertEqual([r.neg for r in rows], [False, False])      # no row re-signed
        vals = [_site_values(patterns.realise(_spec(A, q, (3, 1, 1), phase=r.phase), f,
                                              (3, 1, 1)), f, (3, 1, 1)) for r in rows]
        scale = float(np.max(np.abs(vals[0])))
        m1, m2 = vals[0] / scale, vals[1] / scale
        np.testing.assert_allclose(m1, [1, -.5, -.5], atol=1e-12)
        np.testing.assert_allclose(m2, [.5, -1, .5], atol=1e-12)
        self.assertAlmostEqual(float(np.sum(m1 ** 3)), 0.75, delta=1e-12)
        self.assertAlmostEqual(float(np.sum(m2 ** 3)), -0.75, delta=1e-12)
        # the node-centred phase is not emitted
        node = _site_values(patterns.realise(_spec(A, q, (3, 1, 1), phase=phi0 + math.pi / 6),
                                             f, (3, 1, 1)), f, (3, 1, 1)) / scale
        np.testing.assert_allclose(node, [.87, -.87, 0], atol=0.01)
        self.assertAlmostEqual(float(np.sum(node ** 3)), 0.0, delta=1e-12)
        for r in rows:
            self.assertGreater(abs(math.sin(r.phase - (phi0 + math.pi / 6))), 0.1)
        # the second row is the translate of -O(phi0)
        O2 = patterns.realise(_spec(A, q, (3, 1, 1), phase=rows[1].phase), f, (3, 1, 1))
        self.assertEqual(symmetry.find_map(-O, O2, problem=f.problem, su2=su2, tr=tr),
                         "translation")

    def test_q_third_tr_odd_with_map(self):
        """A TR-odd spin pattern at q = 1/3 with 2Sz free, a TR-invariant
        H_base and SU(2) failing (Ising): time reversal maps O to -O, so the
        rows are phi0 and phi0 + pi/6."""
        from hwave.initial import patterns, symmetry
        f = self.chain(6, (3, 1, 1), ising=0.3)
        A = (_onsite(cz=1.0),)
        q = ("1/3", "0", "0")
        O = patterns.realise(_spec(A, q, (3, 1, 1)), f, (3, 1, 1))
        su2 = symmetry.su2_check(f.problem)
        tr = symmetry.tr_check(f.problem)
        self.assertEqual(su2.status, "failed")
        self.assertTrue(tr.invariant)
        sign_map = symmetry.sign_map_verified(O, problem=f.problem, su2=su2, tr=tr)
        self.assertTrue(sign_map)
        self.assertEqual(symmetry.find_map(O, -O, problem=f.problem, su2=su2, tr=tr),
                         "time_reversal")
        rows = patterns.phase_rows(A, q, sign_map=sign_map)
        phi0 = patterns.phi0(A, q)
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(rows[0].phase, phi0, places=15)
        self.assertAlmostEqual(rows[1].phase, phi0 + math.pi / 6, places=15)
        self.assertEqual([r.equivalent_to for r in rows], [None, None])
        self.assertEqual([r.neg for r in rows], [False, False])

    def test_q_half_third_order_six(self):
        from hwave.initial import cells, patterns
        q = ("1/2", "1/3", "0")
        self.assertEqual(cells.order(q), 6)
        A = (_onsite(c0=1.0),)
        for sign_map in (False, True):
            rows = patterns.phase_rows(A, q, sign_map=sign_map)
            self.assertEqual([r.phase for r in rows], [0.0, math.pi / 6])
            self.assertEqual([r.neg for r in rows], [False, False])
            self.assertEqual([r.equivalent_to for r in rows], [None, None])

    def test_self_conjugate_rows(self):
        """q = -q (q = 0 included): one phase, or phi0 and phi0 + pi/2 when
        tr AA vanishes; each with its -O row unless a map merges them."""
        from hwave.initial import patterns
        A = (_onsite(c0=1.0),)
        for q in (("0", "0", "0"), ("1/2", "1/2", "0")):
            rows = patterns.phase_rows(A, q, sign_map=True)
            self.assertEqual([(r.phase, r.neg) for r in rows], [(0.0, False)])
            rows = patterns.phase_rows(A, q, sign_map=False)
            self.assertEqual([(r.phase, r.neg) for r in rows], [(0.0, False), (0.0, True)])
        A = (_onsite(cx=1.0, cy=1j),)          # |tr AA| = 0
        rows = patterns.phase_rows(A, ("0", "0", "0"), sign_map=False)
        self.assertEqual([(r.phase, r.neg) for r in rows],
                         [(0.0, False), (0.0, True), (math.pi / 2, False),
                          (math.pi / 2, True)])
        rows = patterns.phase_rows(A, ("0", "0", "0"), sign_map=True)
        self.assertEqual([(r.phase, r.neg) for r in rows],
                         [(0.0, False), (math.pi / 2, False)])

    def test_fix_sign(self):
        from hwave.initial import patterns
        q = ("0", "0", "0")
        A = (_onsite(cz=-2.0), _onsite(c0=1.0, a=1))
        got = patterns.fix_sign(A, q)
        self.assertEqual([e.cz for e in got], [2.0, 0.0])
        self.assertEqual([e.c0 for e in got], [0.0, -1.0])
        self.assertIs(type(got), tuple)
        # already positive: unchanged
        self.assertEqual(patterns.fix_sign(got, q), got)
        # a raw eigensolver phase: the phi0 row is made positive
        A = (_onsite(cz=-2.0 * cmath.exp(0.7j)),)
        got = patterns.fix_sign(A, q)
        phi = patterns.phi0(got, q)
        o0 = (cmath.exp(1j * phi) * got[0].cz + cmath.exp(-1j * phi) * np.conj(got[0].cz)) / 2
        self.assertGreater(o0.real, 0)
        # ties: the lowest index (orbital 0 before orbital 1) after rounding
        A = (_onsite(cz=-1.0), _onsite(cz=1.0 + 1e-13, a=1))
        got = patterns.fix_sign(A, q)
        self.assertEqual([e.cz.real for e in got], [1.0, -(1.0 + 1e-13)])


# --------------------------------------------------------------------------
# degenerate subspaces
# --------------------------------------------------------------------------

class TestCanonicalBasis(unittest.TestCase):
    def test_canonical_basis_unitary_invariance(self):
        from hwave.initial import patterns
        rng = np.random.default_rng(11)
        n, g = 12, 3
        X = rng.normal(size=(n, g)) + 1j * rng.normal(size=(n, g))
        V, _ = np.linalg.qr(X)
        W = patterns.canonical_basis(V)
        for seed in range(5):
            r = np.random.default_rng(100 + seed)
            U, _ = np.linalg.qr(r.normal(size=(g, g)) + 1j * r.normal(size=(g, g)))
            np.testing.assert_allclose(patterns.canonical_basis(V @ U), W, atol=1e-10)
        # orthonormal, spanning the same subspace
        np.testing.assert_allclose(W.conj().T @ W, np.eye(g), atol=1e-12)
        np.testing.assert_allclose(V @ (V.conj().T @ W), W, atol=1e-12)

    def test_pivot_rule(self):
        """Pivots are the largest-norm coordinates of the deflated V^H, ties
        by the lowest index; the first column is the normalised projection of
        the first pivot's unit vector, its pivot entry real positive."""
        from hwave.initial import patterns
        V = np.zeros((4, 2), dtype=complex)
        V[1, 0] = 1j
        V[3, 1] = -1
        W = patterns.canonical_basis(V)
        np.testing.assert_allclose(W[:, 0], [0, 1, 0, 0], atol=1e-15)
        np.testing.assert_allclose(W[:, 1], [0, 0, 0, 1], atol=1e-15)

    def test_deflation(self):
        """Two parallel largest rows: after the first pivot the second is
        deflated to zero, so the next pivot is the next independent row."""
        from hwave.initial import patterns
        V = np.zeros((6, 2), dtype=complex)
        V[:2, 0] = 1 / math.sqrt(2)
        V[2:, 1] = 0.5
        U, _ = np.linalg.qr(np.array([[1, 2j], [0.5, -1]]))
        W = patterns.canonical_basis(V @ U)
        np.testing.assert_allclose(W[:, 0], V[:, 0], atol=1e-12)
        np.testing.assert_allclose(W[:, 1], V[:, 1], atol=1e-12)

    def test_one_vector(self):
        from hwave.initial import patterns
        v = np.array([[0.1], [-0.6j], [0.3]], dtype=complex)
        v /= np.linalg.norm(v)
        W = patterns.canonical_basis(v)
        self.assertAlmostEqual(complex(W[1, 0]), abs(v[1, 0]), places=15)


# --------------------------------------------------------------------------
# spirals and the rows of one eigenvector
# --------------------------------------------------------------------------

class TestSpiral(_Tmp):
    def test_spiral_random_global_phase(self):
        """An SU(2)-passing spin eigenvector at q = 1/3 with a random global
        phase: exactly 3 rows, 2 collinear + 1 spiral, the spiral equal to
        ``cos(q.R) s (x) sigma^x + sin(q.R) s (x) sigma^y`` for real s."""
        from hwave.initial import patterns, symmetry
        from hwave.initial.patterns import Entry
        f = self.chain(6, (3, 1, 1))
        su2 = symmetry.su2_check(f.problem)
        self.assertEqual(su2.status, "passed")
        q = ("1/3", "0", "0")
        s = {(0, 0, (0, 0, 0)): 0.9, (0, 0, (1, 0, 0)): 0.3, (0, 0, (-1, 0, 0)): 0.3}
        theta = float(np.random.default_rng(5).uniform(0, 2 * math.pi))
        z = cmath.exp(1j * theta)
        # A = e^{i theta} s (x) sigma^z: c^z = 2 s
        A = tuple(Entry(a, b, r, cz=2 * v * z) for (a, b, r), v in s.items())
        O = patterns.realise(_spec(A, q, (3, 1, 1)), f, (3, 1, 1))
        self.assertTrue(symmetry.sign_map_verified(O, problem=f.problem, su2=su2,
                                                   tr=symmetry.tr_check(f.problem)))
        rows = patterns.pattern_rows(A, q, sign_map=True, spiral=True, subshape=(3, 1, 1))
        self.assertEqual(len(rows), 3)
        self.assertEqual([r.spiral for r in rows], [False, False, True])
        self.assertEqual([r.spec.modulation for r in rows], ["cos", "cos", "spiral"])
        spiral = rows[2].spec
        self.assertEqual(spiral.axes, ("x", "y"))
        self.assertAlmostEqual(spiral.phase, rows[0].spec.phase, places=15)
        got = patterns.realise(spiral, f, (3, 1, 1))
        real_s = tuple(Entry(a, b, r, s=v) for (a, b, r), v in s.items())
        want = patterns.realise(_spec(real_s, q, (3, 1, 1), modulation="spiral",
                                      axes=("x", "y")), f, (3, 1, 1))
        np.testing.assert_allclose(got, want, atol=1e-12)
        # the hand form cos(q.R) s sigma^x + sin(q.R) s sigma^y on the lattice
        hand = _oracle(_spec(real_s, q, (3, 1, 1), modulation="spiral", axes=("x", "y")),
                       f, (3, 1, 1))
        np.testing.assert_allclose(got, hand, atol=1e-12)
        o_r = _to_r(got, (2, 1, 1))[0, 0, 0]
        n = 3
        kappa = o_r[0, n].real / 0.9          # the normalisation of the on-site s
        self.assertGreater(kappa, 0)
        for R in range(3):
            c, sn = math.cos(2 * math.pi * R / 3), math.sin(2 * math.pi * R / 3)
            blk = np.array([[o_r[R, R], o_r[R, n + R]], [o_r[n + R, R], o_r[n + R, n + R]]])
            np.testing.assert_allclose(blk, kappa * 0.9 * (c * _SX + sn * _SY), atol=1e-12)
        # independent of the global phase
        A2 = tuple(Entry(a, b, r, cz=2 * v * cmath.exp(2.1j)) for (a, b, r), v in s.items())
        rows2 = patterns.pattern_rows(A2, q, sign_map=True, spiral=True, subshape=(3, 1, 1))
        np.testing.assert_allclose(patterns.realise(rows2[2].spec, f, (3, 1, 1)), got,
                                   atol=1e-12)

    def test_spiral_needs_spin_z_shape(self):
        from hwave.initial import patterns
        with self.assertRaisesRegex(ValueError, "S\\^z"):
            patterns.spiral_spec((_onsite(c0=1.0),), ("1/3", "0", "0"), 0.0)

    def test_spiral_spec_default_cell(self):
        from hwave.initial import patterns
        spec = patterns.spiral_spec((_onsite(cz=2.0),), ("2/3", "1/3", "0"), 0.25)
        self.assertEqual(spec.subshape, (3, 3, 1))
        self.assertEqual(spec.modulation, "spiral")
        self.assertEqual(spec.entries[0].s, 1.0)
        self.assertEqual(spec.phase, 0.25)

    def test_pattern_rows_collinear(self):
        from hwave.initial import patterns
        A = (_onsite(cz=-1.0),)
        rows = patterns.pattern_rows(A, ("1/3", "0", "0"), sign_map=False, subshape=(3, 1, 1))
        self.assertEqual(len(rows), 2)
        self.assertEqual([r.spec.entries[0].cz for r in rows], [1.0, 1.0])   # sign fixed once
        self.assertEqual([r.spec.sign for r in rows], [1, 1])
        self.assertEqual(rows[1].phase_row.equivalent_to, "-O(phi0), translated")
        rows = patterns.pattern_rows(A, ("1/2", "0", "0"), sign_map=False, subshape=(2, 1, 1))
        self.assertEqual([r.spec.sign for r in rows], [1, -1])
        with self.assertRaisesRegex(ValueError, "spiral"):
            patterns.pattern_rows(A, ("1/2", "0", "0"), sign_map=True, spiral=True)


if __name__ == "__main__":
    unittest.main()
