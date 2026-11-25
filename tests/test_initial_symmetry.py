"""Spin-rotation and time-reversal checks and the deduplication maps of
``hwave.initial`` (``hwave.initial.symmetry``; design section 5 step 6, R1).

Every check runs on the real Hartree-Fock map of the solver
(``initial_facade.hf_potential`` on a loaded problem), here on a
two-orbital chain of three cells (``nvol = 3``: a bond interaction needs
both directions ``R = +-1`` declared, and the solver refuses a table wider
than the cell) unless a test needs more cells.

The residuals are those of the check as specified (design section 5 step
6): the largest Frobenius residual over the complete table-support probe
basis, divided by ``max(1, max_e ||Sigma'[e]||)``. On these models they
are exact numbers: 1 for the x and y axes of the Ising model (0 for z),
``sqrt(2)`` for the x and z axes of the revision-3 form on CoulombIntra (0
for its y axis, where ``S_y^T = -S_y`` makes the two forms agree).
"""
import os
import shutil
import tempfile
import time
import unittest

import numpy as np

from tests import initial_fixtures

#: the warning of an S_z check that fails in 2Sz mode (design section 5 step 6)
SZ_WARNING = "the interaction tables break spin conservation about z"

#: two orbitals, spin-independent hopping on both and an on-site mixing
_TRANSFER2 = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0,
              ((1, 0, 0), 1, 1): -0.5, ((-1, 0, 0), 1, 1): -0.5,
              ((0, 0, 0), 0, 1): -0.2, ((0, 0, 0), 1, 0): -0.2}


def _pair(r, a, b, v):
    """``{(r, a, b): v, (-r, b, a): conj(v)}`` (one Hermitian pair)."""
    neg = tuple(-x for x in r)
    return {(tuple(r), a, b): v, (neg, b, a): np.conj(v)}


def _terms(*names, hund=0.8, exchange=0.8):
    """Interaction tables of the two-orbital model by short name."""
    out = {}
    for name in names:
        if name == "intra":
            out["CoulombIntra"] = {((0, 0, 0), 0, 0): 4.0, ((0, 0, 0), 1, 1): 3.0}
        elif name == "inter":            # a bond: orbital 0 to orbital 1 one cell away
            out["CoulombInter"] = {**_pair((1, 0, 0), 0, 1, 1.5),
                                   **_pair((0, 0, 0), 0, 1, 2.0)}
        elif name == "hund":
            out["Hund"] = _pair((0, 0, 0), 0, 1, hund)
        elif name == "exchange":
            out["Exchange"] = _pair((0, 0, 0), 0, 1, exchange)
        elif name == "pairhop":
            out["PairHop"] = _pair((0, 0, 0), 0, 1, 0.6)
        elif name == "pairlift":
            out["PairLift"] = _pair((0, 0, 0), 0, 1, 0.6)
        elif name == "ising":
            out["Ising"] = _pair((0, 0, 0), 0, 1, 1.0)
        else:
            raise KeyError(name)
    return out


def _spin_transfer_text(norb, entries):
    """A 9-column (spin-resolved) Transfer file of ``{(R, a, b, s, t): v}``
    (0-based), Hermitian-closed."""
    closed = {}
    for (r, a, b, s, t), v in entries.items():
        closed[(r, a, b, s, t)] = closed.get((r, a, b, s, t), 0) + v
        neg = tuple(-x for x in r)
        if (neg, b, a, t, s) != (r, a, b, s, t):
            closed[(neg, b, a, t, s)] = closed.get((neg, b, a, t, s), 0) + np.conj(v)
    lines = ["Transfer in wannier90-like format for uhfk w90spin/1",
             str(norb), "1", " 1"]
    for (r, a, b, s, t), v in closed.items():
        v = complex(v)
        lines.append("{} {} {} {} {} {} {} {:.15g} {:.15g}".format(
            *r, a + 1, b + 1, s, t, v.real, v.imag))
    return "\n".join(lines) + "\n"


def _spin_independent(transfer):
    """The 9-column entries of a spin-independent 7-column table
    (one direction of each pair; the writer closes them)."""
    out = {}
    for (r, a, b), v in transfer.items():
        neg = tuple(-x for x in r)
        if (neg, b, a, 0, 0) in out:
            continue
        for s in range(2):
            out[(r, a, b, s, s)] = v
    return out


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_symmetry_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def model(self, *names, fock=True, two_sz=None, transfer=None, cellshape=(3, 1, 1),
              spin_transfer=None, hund=0.8, exchange=0.8, mode=None, norb=2, ncond=6):
        """The loaded two-orbital problem with the named interactions."""
        from hwave.solver import initial_facade as fac
        param = {"Ncond": ncond, "T": 0.1, "SubShape": [1, 1, 1]}
        if two_sz is not None:
            param["2Sz"] = two_sz
        m = dict(mode or {})
        if not fock:
            m["flag_fock"] = False
        extra = None
        if spin_transfer is not None:
            extra = {"transfer.dat": _spin_transfer_text(norb, spin_transfer)}
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=cellshape, norb=norb,
            transfer=_TRANSFER2 if transfer is None else transfer,
            interactions=_terms(*names, hund=hund, exchange=exchange),
            param=param, mode=m or None, extra_files=extra)
        return fac.load_problem(toml)


def _random_green(problem, rng):
    """A generic stored Green function obeying G(r)_ab = conj(G(-r)_ba)."""
    from hwave.initial import symmetry
    nvol, ns, norb = problem.nvol, problem.ns, problem.norb
    g = (rng.normal(size=(nvol, ns, norb, ns, norb))
         + 1j * rng.normal(size=(nvol, ns, norb, ns, norb)))
    return symmetry.hermitian_green(g, problem.shape)


def _spin_ops(norb):
    """``S_mu = sigma^mu / 2 (x) 1_norb`` in the spin-block order."""
    sig = (np.array([[0, 1], [1, 0]], dtype=complex),
           np.array([[0, -1j], [1j, 0]], dtype=complex),
           np.array([[1, 0], [0, -1]], dtype=complex))
    return [np.kron(s / 2, np.eye(norb)) for s in sig]


def _rev3_form(problem):
    """The generator check of design revision 3, kept here only to guard
    against its return: the rho-form identity ``Sigma'[i[S, e]] = i[S,
    Sigma'[e]]`` applied to the STORED probe ``e_G`` without the transpose.
    Returns the relative residual per axis (x, y, z)."""
    from hwave.initial import symmetry
    from hwave.solver import initial_facade as fac
    tables = fac.hf_tables(problem)

    def sig(g):
        return fac.hf_potential(g, tables, include_fock=problem.iflag_fock,
                                blocks=problem.blocks)

    nd = problem.nd
    shp = (problem.nvol, problem.ns, problem.norb, problem.ns, problem.norb)
    res = np.zeros(3)
    norm = 0.0
    for e in symmetry.table_probes(problem):
        s0 = sig(e)
        norm = max(norm, float(np.linalg.norm(s0)))
        e2 = e.reshape(problem.nvol, nd, nd)
        for mu, S in enumerate(_spin_ops(problem.norb)):
            lhs = sig((1j * (S @ e2 - e2 @ S)).reshape(shp))
            rhs = 1j * (S @ s0 - s0 @ S)
            res[mu] = max(res[mu], float(np.linalg.norm(lhs - rhs)))
    return res / max(1.0, norm)


class TestSu2SelfTest(_Tmp):
    """R1: the check itself, on the real ``accumulate_hf``."""

    def test_selftest_intra_fock_passes(self):
        from hwave.initial import symmetry
        r = symmetry.su2_check(self.model("intra"))
        self.assertEqual(r.status, "passed")
        self.assertEqual(r.method, "generator")
        self.assertLessEqual(r.max_residual, 2e-16)
        self.assertIsNone(r.sz_residual)

    def test_selftest_ising_fails(self):
        from hwave.initial import symmetry
        r = symmetry.su2_check(self.model("ising"))
        self.assertEqual(r.status, "failed")
        x, y, z = r.residuals
        np.testing.assert_allclose([x, y], [1.0, 1.0], rtol=1e-12)
        self.assertLessEqual(z, 1e-12)          # Ising conserves S^z

    def test_rev3_form_guard(self):
        p = self.model("intra")
        x, y, z = _rev3_form(p)
        np.testing.assert_allclose([x, z], [np.sqrt(2), np.sqrt(2)], rtol=1e-12)
        self.assertLessEqual(y, 1e-15)
        # while the stored-form check passes on the same model
        from hwave.initial import symmetry
        self.assertEqual(symmetry.su2_check(p).status, "passed")


class TestSu2CaseTable(_Tmp):
    def test_case_table(self):
        from hwave.initial import symmetry
        passing = [("intra",), ("intra", "inter", "hund", "exchange"), ("pairhop",),
                   ("intra", "inter", "hund", "exchange", "pairhop")]
        for names in passing:
            with self.subTest(names=names, fock=True):
                r = symmetry.su2_check(self.model(*names))
                self.assertEqual(r.status, "passed")
                self.assertLessEqual(r.max_residual, 2e-16)
            with self.subTest(names=names, fock=False):
                r = symmetry.su2_check(self.model(*names, fock=False))
                self.assertEqual(r.status, "failed")
                x, y, _ = r.residuals
                self.assertGreater(min(x, y), 0.1)      # O(1)
                self.assertLess(max(x, y), 10.0)
        for names, kw in ((("pairlift",), {}),
                          (("intra", "hund", "exchange"), {"exchange": 0.5})):
            with self.subTest(names=names, **kw):
                r = symmetry.su2_check(self.model(*names, **kw))
                self.assertEqual(r.status, "failed")
                self.assertGreater(r.max_residual, 1e-3)

    def test_spin_dependent_transfer_fails_on_h0(self):
        """R5: a 9-column Transfer with a Zeeman-like or SOC piece breaks
        SU(2) through H0 alone, the interaction tables being SU(2)."""
        from hwave.initial import symmetry
        base = _spin_independent(_TRANSFER2)
        # the 9-column spin-independent file is the 7-column model
        p = self.model("intra", spin_transfer=base)
        self.assertTrue(p.spin_resolved_transfer)
        r = symmetry.su2_check(p)
        self.assertEqual(r.status, "passed")
        self.assertEqual(r.h0_residual, 0.0)
        for extra in ({((0, 0, 0), 0, 0, 0, 0): 0.3, ((0, 0, 0), 0, 0, 1, 1): -0.3},
                      {((0, 0, 0), 0, 1, 0, 0): 0.2j, ((0, 0, 0), 0, 1, 1, 1): -0.2j}):
            with self.subTest(extra=extra):
                entries = dict(base)
                entries.update(extra)
                r = symmetry.su2_check(self.model("intra", spin_transfer=entries))
                self.assertEqual(r.status, "failed")
                self.assertGreater(r.h0_residual, 1e-3)
                # the interaction part alone is equivariant
                self.assertLessEqual(max(r.residuals), 2e-16)

    def test_random_method_above_the_probe_budget(self):
        """``4 N_probe > 2e4`` switches to the random regression; it
        agrees with the generator verdict on both sides."""
        from hwave.initial import symmetry
        for names, want in ((("intra", "inter"), "passed"), (("ising",), "failed")):
            with self.subTest(names=names):
                p = self.model(*names)
                n = len(symmetry.table_probes(p))
                r = symmetry.su2_check(p, max_hf_applications=4 * n - 1)
                self.assertEqual(r.method, "random")
                self.assertEqual(r.status, want)
                self.assertIsNone(r.residuals)
                self.assertEqual(symmetry.su2_check(p).method, "generator")
                # exactly at the budget the generator check still runs
                self.assertEqual(symmetry.su2_check(p, max_hf_applications=4 * n).method,
                                 "generator")
        self.assertEqual(symmetry.MAX_HF_APPLICATIONS, 20000)

    def test_probes_are_a_hermitian_table_support_basis(self):
        from hwave.initial import symmetry
        p = self.model("intra", "inter")
        probes = symmetry.table_probes(p)
        # every probe is a valid stored Green function
        for e in probes:
            np.testing.assert_array_equal(symmetry.hermitian_green(e, p.shape), e)
        # they are linearly independent over the reals
        m = np.array([np.concatenate([e.ravel().real, e.ravel().imag]) for e in probes])
        self.assertEqual(np.linalg.matrix_rank(m), len(probes))
        # and span the slots the tables read: a random Hermitian G with the
        # off-support entries zeroed gives the same potential as its projection
        from hwave.solver import initial_facade as fac
        g = _random_green(p, np.random.default_rng(1))
        coef, *_ = np.linalg.lstsq(m.T, np.concatenate([g.ravel().real, g.ravel().imag]),
                                   rcond=None)
        proj = np.tensordot(coef, np.array(probes), axes=1)
        kw = dict(include_fock=True, blocks=p.blocks)
        np.testing.assert_allclose(fac.hf_potential(proj, p.tables, **kw),
                                   fac.hf_potential(g, p.tables, **kw), atol=1e-12)


class TestFiniteRotation(_Tmp):
    def test_finite_rotation_both_forms(self):
        """theta = 0.37: rho form ``Sigma'_rho[U rho U^dag] = U Sigma'_rho[rho]
        U^dag`` with ``Sigma'_rho[rho] = hf_potential(rho^T)``, and the stored
        form ``Sigma'[U* G U^T] = U Sigma'[G] U^dag``, agree on a generic G."""
        from hwave.initial import symmetry
        from hwave.solver import initial_facade as fac
        theta = 0.37
        for names, passes in ((("intra", "inter", "hund", "exchange", "pairhop"), True),
                              (("ising",), False)):
            p = self.model(*names)
            nd = p.nd
            shp = (p.nvol, p.ns, p.norb, p.ns, p.norb)

            def sig(g):
                return fac.hf_potential(g, p.tables, include_fock=True, blocks=p.blocks)

            rng = np.random.default_rng(7)
            G = _random_green(p, rng)
            n = rng.normal(size=3)
            n /= np.linalg.norm(n)
            U = symmetry.spin_rotation(n, theta, p.norb)
            # stored form, through the module's own rotation of G
            s_res = np.linalg.norm(sig(symmetry.rotate_green(G, U))
                                   - U @ sig(G) @ U.conj().T)
            # rho form, written out here: rho(r) = G(r)^T
            rho = np.swapaxes(G.reshape(p.nvol, nd, nd), -1, -2)
            rho_rot = U @ rho @ U.conj().T

            def sig_rho(r):
                return sig(np.swapaxes(r, -1, -2).reshape(shp))

            r_res = np.linalg.norm(sig_rho(rho_rot) - U @ sig_rho(rho) @ U.conj().T)
            scale = max(1.0, float(np.linalg.norm(sig(G))))
            with self.subTest(names=names):
                np.testing.assert_allclose(s_res, r_res, rtol=1e-10, atol=1e-13)
                if passes:
                    self.assertLessEqual(s_res / scale, 1e-13)
                else:
                    self.assertGreater(s_res / scale, 1e-3)


class TestNotApplicable(_Tmp):
    def test_not_applicable_labels(self):
        from hwave.initial import symmetry
        # 2Sz fixed: only S_z is checked
        r = symmetry.su2_check(self.model("intra", two_sz=0))
        self.assertEqual(r.status, "not applicable (2Sz fixed)")
        self.assertIsNotNone(r.sz_residual)
        self.assertLessEqual(r.sz_residual, 2e-16)
        self.assertEqual(r.warnings, ())
        # spin-orbital mode never qualifies
        so_transfer = {}
        for (rv, a, b), v in _TRANSFER2.items():
            for s in range(2):
                so_transfer[(rv, 2 * a + s, 2 * b + s)] = v
        p = self.model("intra", transfer=so_transfer, norb=4,
                       mode={"enable_spin_orbital": True})
        r = symmetry.su2_check(p)
        self.assertEqual(r.status, "not applicable (spin-orbital)")
        self.assertIsNone(r.max_residual)
        self.assertIsNone(r.method)
        self.assertIsNone(r.sz_residual)

    def test_2sz_sz_failure_warns(self):
        """In 2Sz mode the blocks mask every transverse output, so the
        solver's tables conserve S_z there. Without that mask, PairLift
        (``G^{st}_bb -> Sigma^{st}_aa``, the wrong phase under a rotation
        about z) breaks it: the residual goes to ``sz_residual`` with the
        exact warning and the status does not change."""
        import dataclasses
        from hwave.initial import symmetry
        from hwave.solver import initial_facade as fac
        p = self.model("pairlift", two_sz=0)
        r = symmetry.su2_check(p)
        self.assertEqual(r.status, "not applicable (2Sz fixed)")
        self.assertLessEqual(r.sz_residual, 1e-12)
        self.assertEqual(r.warnings, ())
        one = fac.SolverBlocks((np.arange(p.nd),), (0,), (p.ncond,))
        unmasked = dataclasses.replace(p, blocks=one)
        with self.assertLogs("hwave.initial", level="WARNING") as cm:
            r = symmetry.su2_check(unmasked)
        self.assertEqual(r.status, "not applicable (2Sz fixed)")
        self.assertGreater(r.sz_residual, 1e-3)
        self.assertEqual(r.residuals[:2], (None, None))
        self.assertEqual(r.warnings, (SZ_WARNING,))
        self.assertTrue(any(SZ_WARNING in line for line in cm.output))


class TestTimeReversal(_Tmp):
    def test_tr_invariance(self):
        from hwave.initial import symmetry
        chain = dict(cellshape=(4, 1, 1), ncond=8)
        # real hopping
        r = symmetry.tr_check(self.model("intra", "inter", "hund", **chain))
        self.assertTrue(r.invariant)
        self.assertLessEqual(max(r.h0_residual, r.hf_residual), 1e-12)
        # a complex hopping e^{0.3i} breaks it through H0
        cplx = dict(_TRANSFER2)
        cplx[((1, 0, 0), 0, 0)] = -np.exp(0.3j)
        cplx[((-1, 0, 0), 0, 0)] = -np.exp(-0.3j)
        r = symmetry.tr_check(self.model("intra", transfer=cplx, **chain))
        self.assertFalse(r.invariant)
        self.assertGreater(r.h0_residual, 0.1)
        self.assertLessEqual(r.hf_residual, 1e-12)
        # SOC i lambda sigma^z between the two orbitals is TR-even:
        # 9-column normal mode ...
        lam = 0.25
        entries = _spin_independent(_TRANSFER2)
        entries[((0, 0, 0), 0, 1, 0, 0)] = 1j * lam
        entries[((0, 0, 0), 0, 1, 1, 1)] = -1j * lam
        r = symmetry.tr_check(self.model("intra", spin_transfer=entries, **chain))
        self.assertTrue(r.invariant)
        # ... and spin-orbital mode (interleaved 2*orb + spin)
        so = {}
        for (rv, a, b), v in _TRANSFER2.items():
            for s in range(2):
                so[(rv, 2 * a + s, 2 * b + s)] = v
        for s, sign in ((0, 1), (1, -1)):
            so[((0, 0, 0), 0 + s, 2 + s)] = 1j * lam * sign
            so[((0, 0, 0), 2 + s, 0 + s)] = -1j * lam * sign
        r = symmetry.tr_check(self.model("intra", "inter", transfer=so, norb=4,
                                         mode={"enable_spin_orbital": True}, **chain))
        self.assertTrue(r.invariant)
        self.assertLessEqual(max(r.h0_residual, r.hf_residual), 1e-12)
        # a Zeeman-like 9-column piece is TR-odd
        zee = _spin_independent(_TRANSFER2)
        zee[((0, 0, 0), 0, 0, 0, 0)] = 0.3
        zee[((0, 0, 0), 0, 0, 1, 1)] = -0.3
        r = symmetry.tr_check(self.model("intra", spin_transfer=zee, **chain))
        self.assertFalse(r.invariant)
        self.assertGreater(r.h0_residual, 0.1)
        # 2Sz = 2: the sector is not mapped onto itself
        r = symmetry.tr_check(self.model("intra", two_sz=2, **chain))
        self.assertFalse(r.invariant)
        # 2Sz = 0 is
        r = symmetry.tr_check(self.model("intra", two_sz=0, **chain))
        self.assertTrue(r.invariant)

    def test_tr_equivariance_catches_a_tr_odd_map(self):
        """The HF part is a numerical check of the map, not a reading of
        the tables: the solver's Hartree-Fock map on a complex PairHop
        coefficient is tested as it is."""
        from hwave.initial import symmetry
        r = symmetry.tr_check(self.model("intra", "pairhop", cellshape=(4, 1, 1), ncond=8))
        self.assertTrue(r.invariant)
        self.assertLessEqual(r.hf_residual, 1e-12)


def _spin_pattern(problem, c, *, spin=None, orb=None):
    """``O(R)`` (Hamiltonian convention, folded cell) with ``O = 1/2 sum
    c^mu sigma^mu`` per folded orbital: ``c`` ``(norb, 4)`` gives the
    on-site weights (charge, x, y, z) of each folded orbital; returned in
    k-space ``(nvol, nd, nd)`` as ``realise`` gives it."""
    from hwave.solver import hamiltonian_k
    sig = (np.eye(2, dtype=complex), np.array([[0, 1], [1, 0]], dtype=complex),
           np.array([[0, -1j], [1j, 0]], dtype=complex),
           np.array([[1, 0], [0, -1]], dtype=complex))
    norb = problem.norb
    o_r = np.zeros(tuple(problem.shape) + (2, norb, 2, norb), dtype=complex)
    for a in range(norb):
        m = sum(c[a][mu] * sig[mu] for mu in range(4)) / 2
        o_r[0, 0, 0, :, a, :, a] = m
    nd = 2 * norb
    return hamiltonian_k.to_k(o_r.reshape(tuple(problem.shape) + (nd, nd)))


class TestMaps(_Tmp):
    def square(self, *, U=4.0, fock=True, subshape=(2, 2, 1), L=4):
        from hwave.solver import initial_facade as fac
        toml = initial_fixtures.square(self.dest(), L, U=U, T=0.1, subshape=subshape,
                                       fock=fock)
        return fac.load_problem(toml)

    def test_maps(self):
        from hwave.initial import symmetry
        p = self.square()                       # 2x2 folded cell, one orbital
        su2 = symmetry.su2_check(p)
        tr = symmetry.tr_check(p)
        self.assertEqual(su2.status, "passed")
        self.assertTrue(tr.invariant)
        # Neel at (1/2, 1/2): sub-cell ir = ix + 2 iy, sign (-1)^(ix+iy)
        eta = [1, -1, -1, 1]
        neel = _spin_pattern(p, [[0, 0, 0, e] for e in eta])
        self.assertEqual(symmetry.find_map(neel, -neel, problem=p, su2=su2, tr=tr),
                         "translation")
        self.assertTrue(symmetry.sign_map_verified(neel, problem=p, su2=su2, tr=tr))
        # a TR-odd ferro: translations leave it alone, time reversal flips it
        ferro = _spin_pattern(p, [[0, 0, 0, 1]] * 4)
        self.assertEqual(symmetry.find_map(ferro, ferro, problem=p, su2=su2, tr=tr),
                         "translation")
        self.assertEqual(symmetry.find_map(ferro, -ferro, problem=p, su2=su2, tr=tr),
                         "time_reversal")
        # without time reversal, a global spin rotation still maps it
        no_tr = symmetry.TrResult(invariant=False, h0_residual=1.0, hf_residual=0.0)
        self.assertEqual(symmetry.find_map(ferro, -ferro, problem=p, su2=su2, tr=no_tr),
                         "spin_rotation")
        ferro_x = _spin_pattern(p, [[0, 1, 0, 0]] * 4)
        self.assertEqual(symmetry.find_map(ferro, ferro_x, problem=p, su2=su2, tr=no_tr),
                         "spin_rotation")
        # never assume SU(2): a failed (or inapplicable) check maps nothing
        for status in ("failed", "not applicable (2Sz fixed)"):
            off = symmetry.Su2Result(status=status, max_residual=1.0, method="generator",
                                     sz_residual=None)
            self.assertIsNone(symmetry.find_map(ferro, ferro_x, problem=p, su2=off,
                                                tr=no_tr))
            self.assertFalse(symmetry.sign_map_verified(ferro, problem=p, su2=off,
                                                        tr=no_tr))
        # a translated Neel is found as a translation, not as -O
        self.assertEqual(symmetry.find_map(neel, neel, problem=p, su2=su2, tr=tr),
                         "translation")

    def test_crystal_field_orbital_pattern_has_no_sign_map(self):
        """A crystal-field q = 0 orbital pattern (TR-even, translation
        invariant, spin-independent) is never mapped to -O."""
        from hwave.initial import symmetry
        p = self.model("intra", "inter", "hund", "exchange")
        su2 = symmetry.su2_check(p)
        tr = symmetry.tr_check(p)
        self.assertEqual(su2.status, "passed")
        self.assertTrue(tr.invariant)
        orb = _spin_pattern(p, [[1, 0, 0, 0], [-1, 0, 0, 0]])
        self.assertIsNone(symmetry.find_map(orb, -orb, problem=p, su2=su2, tr=tr))
        self.assertFalse(symmetry.sign_map_verified(orb, problem=p, su2=su2, tr=tr))
        self.assertEqual(symmetry.find_map(orb, orb, problem=p, su2=su2, tr=tr),
                         "translation")

    def test_non_coplanar_minus_needs_time_reversal(self):
        """Three non-coplanar moments: -O is not a proper rotation, so
        without time reversal it is not mapped."""
        from hwave.initial import symmetry
        p = self.square(subshape=(2, 2, 1))
        su2 = symmetry.su2_check(p)
        tr = symmetry.tr_check(p)
        no_tr = symmetry.TrResult(invariant=False, h0_residual=1.0, hf_residual=0.0)
        c = [[0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1], [0, 0.3, 0.2, 0.1]]
        o = _spin_pattern(p, c)
        self.assertIsNone(symmetry.find_map(o, -o, problem=p, su2=su2, tr=no_tr))
        self.assertEqual(symmetry.find_map(o, -o, problem=p, su2=su2, tr=tr),
                         "time_reversal")

    def test_translation_moves_bonds_across_the_cell(self):
        """A dimerised bond pattern on a chain folded by B = 2: the
        translation by one site takes the even bonds (inside the cell) to
        the odd bonds (across the cell boundary, R = +1 for ``1 -> 0``), and
        not to the reversed even bond a wrong carry sign would give."""
        from hwave.initial import symmetry
        from hwave.solver import initial_facade as fac
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(8, 1, 1), transfer={((1, 0, 0), 0, 0): -1.0,
                                                       ((-1, 0, 0), 0, 0): -1.0},
            param={"Ncond": 8, "T": 0.1, "SubShape": [2, 1, 1]})
        p = fac.load_problem(toml)
        self.assertEqual((p.shape, p.norb), ((4, 1, 1), 2))
        nd = p.nd
        off = symmetry.Su2Result(status="failed", max_residual=1.0, method="generator",
                                 sz_residual=None)
        no_tr = symmetry.TrResult(invariant=False, h0_residual=1.0, hf_residual=0.0)

        def bonds(entries):
            o_r = np.zeros((4, 1, 1, nd, nd), dtype=complex)
            for R, i, j, v in entries:
                for s in range(2):
                    o_r[R % 4, 0, 0, s * 2 + i, s * 2 + j] += v
            return np.fft.ifftn(o_r, axes=(0, 1, 2), norm="forward").reshape(4, nd, nd)

        even = bonds([(0, 0, 1, 1.0), (0, 1, 0, 1.0)])
        odd = bonds([(1, 1, 0, 1.0), (-1, 0, 1, 1.0)])
        wrong = bonds([(-1, 1, 0, 1.0), (1, 0, 1, 1.0)])
        self.assertEqual(symmetry.find_map(even, odd, problem=p, su2=off, tr=no_tr),
                         "translation")
        self.assertIsNone(symmetry.find_map(even, wrong, problem=p, su2=off, tr=no_tr))
        # a directed (current-like) pattern keeps its direction
        cur = bonds([(0, 0, 1, 1j), (0, 1, 0, -1j)])
        cur_odd = bonds([(1, 1, 0, 1j), (-1, 0, 1, -1j)])
        self.assertEqual(symmetry.find_map(cur, cur_odd, problem=p, su2=off, tr=no_tr),
                         "translation")
        self.assertIsNone(symmetry.find_map(cur, -cur_odd, problem=p, su2=off, tr=no_tr))
        # time reversal flips a current: allowed only when verified
        tr = symmetry.tr_check(p)
        self.assertTrue(tr.invariant)
        self.assertEqual(symmetry.find_map(cur, -cur_odd, problem=p, su2=off, tr=tr),
                         "time_reversal")
        self.assertTrue(symmetry.sign_map_verified(cur, problem=p, su2=off, tr=tr))
        self.assertFalse(symmetry.sign_map_verified(even, problem=p, su2=off, tr=tr))

    def test_tolerance_is_entrywise_relative(self):
        from hwave.initial import symmetry
        p = self.square()
        su2 = symmetry.su2_check(p)
        tr = symmetry.tr_check(p)
        crys = _spin_pattern(p, [[1, 0, 0, 0], [0.5, 0, 0, 0], [0.5, 0, 0, 0], [1, 0, 0, 0]])
        scale = np.linalg.norm(symmetry._to_r(crys, p))
        bump = np.zeros_like(crys)
        bump[:, 0, 0] = 1e-9 * scale        # a uniform shift of one entry in R = 0
        self.assertIsNone(symmetry.find_map(crys, crys + bump, problem=p, su2=su2, tr=tr))
        bump *= 1e-3
        self.assertEqual(symmetry.find_map(crys, crys + bump, problem=p, su2=su2, tr=tr),
                         "translation")

    def test_shape_mismatch_raises(self):
        from hwave.initial import symmetry
        p = self.square()
        su2 = symmetry.su2_check(p)
        tr = symmetry.tr_check(p)
        o = _spin_pattern(p, [[0, 0, 0, 1]] * 4)
        with self.assertRaises(ValueError):
            symmetry.find_map(o, o[:1], problem=p, su2=su2, tr=tr)


if __name__ == "__main__":
    unittest.main()
