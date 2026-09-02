#!/usr/bin/env python3

"""Regression tests for issue #174: the ``general`` scheme must build the
spinful exchange crossing for EVERY run that can reach the spinful solve
branch, not only for ``enable_spin_orbital = true`` runs.

``RPA._calc_epsilon_k`` sets the spin-orbital branch UNCONDITIONALLY for
the ``trans_mod`` / ``green_init`` npz routes, and ``spin_mode`` is then
read off the supplied H0(k)'s block structure. So a spin-mixing H0 handed
in through an npz reaches ``spin_mode == 'spinful'`` with
``enable_spin_orbital = False``. Before the fix, the precursor of the #137
antisymmetrized crossing (``onsite_r`` in ``Interaction._make_ham_inter``)
was gated on the ``enable_spin_orbital`` flag, so on that route
``ham_spinful_exchange`` was ``None`` and the spinful solve silently fell
back to the pre-#137 ring-only vertex -- measured 1.705e-02 (12.2%
relative) in ``chiq``, confined to the spin-flip pair slots, with no
warning. ``calc_scheme = 'auto'`` promotes INTO that path since #167
(``auto:mixed:trans_mod`` / ``auto:mixed:green_init``).

Each test compares the npz route (route B: ``enable_spin_orbital=False``,
physical Geometry/Transfer, spinful npz) against the ADJUDICATED
spin-orbital run (route A: ``enable_spin_orbital=True``, interleaved SO
Transfer) on IDENTICAL physics. The two routes assemble bit-identical
H0(k), ``chi0q`` and ``ham_inter_q``; the only thing that ever differed
was the vertex crossing, so the comparison is BITWISE.

The reduced scheme is unaffected (the crossing has no density-diagonal
content); that half is recorded in ``tests/test_rpa_spinful_npz_reduced.py``
(issue #161).

The crossing is materialized LAZILY by ``Interaction.spinful_exchange()``
(round-2 review): the dense tensor is (ns*norb)^4 complex128 -- about
1.6 GB at norb = 50 -- and ``ham_spinful_exchange`` keeps a transposed
view of it, so building it in the constructor would pin that allocation
for every solver instance, including the spin-free and spin-diag solves
that never read it. The contract these tests pin is therefore: absent
after CONSTRUCTION in every mode; present after a SPINFUL solve; still
absent after a non-spinful solve.

Tests must run from the repository root.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as solver_rpa


NORB = 2          # physical orbitals
NS = 2
ND = NS * NORB
LX = 8            # 1D chain cells
T = 2.0
FILLING = 0.5
NMAT = 32


# --------------------------------------------------------------------------
# fixture: ONE physical spin-mixing H0, emitted in both index conventions
# --------------------------------------------------------------------------

def _h_blocks():
    """H(R) in SPIN-BLOCK order (index = s*NORB + a) for R = 0, +1, -1.

    The R = 0 block carries a genuinely spin-MIXING (spin off-diagonal,
    non-symmetric, complex) sector -- that is what drives ``spin_mode``
    to ``'spinful'``.
    """
    H = {}
    h0 = np.zeros((ND, ND), dtype=complex)
    h0[0:2, 0:2] = [[0.30, 0.10], [0.10, -0.20]]      # up-up
    h0[2:4, 2:4] = [[-0.10, 0.05], [0.05, 0.25]]      # down-down
    m = np.array([[0.15 + 0.05j, 0.02], [0.03, -0.07 + 0.04j]], dtype=complex)
    h0[0:2, 2:4] = m
    h0[2:4, 0:2] = m.conj().T
    H[0] = h0

    hp = np.zeros((ND, ND), dtype=complex)
    hp[0:2, 0:2] = [[-0.50, 0.07], [0.04, -0.45]]
    hp[2:4, 2:4] = [[-0.40, 0.02], [0.01, -0.35]]
    hp[0:2, 2:4] = [[0.06j, 0.01], [0.02, 0.03]]
    hp[2:4, 0:2] = [[0.0, 0.015], [0.025, 0.0]]
    H[1] = hp
    H[-1] = hp.conj().T
    return H


def _spinblock_to_interleaved(mat):
    """spin-block (s*NORB + a) -> interleaved (2*a + s), on both axes."""
    perm = [(j % NS) * NORB + (j // NS) for j in range(ND)]
    return mat[np.ix_(perm, perm)]


def _w90_header(fh, norb, nr, comment):
    fh.write(comment + "\n")
    fh.write("{}\n".format(norb))
    fh.write("{}\n".format(nr))
    fh.write(" ".join(["1"] * nr) + "\n")


def write_fixture(d):
    """Write geometry / transfer / interaction / npz files into ``d``."""
    H = _h_blocks()

    for name, n in (("geom_so.dat", ND), ("geom_phys.dat", NORB)):
        with open(os.path.join(d, name), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n")
            f.write("{}\n".format(n))
            for _ in range(n):
                f.write("  0.0 0.0 0.0\n")

    # SO Transfer: interleaved indices, the spin-mixing H0 itself
    with open(os.path.join(d, "transfer_so_mix.dat"), "w") as f:
        _w90_header(f, ND, 3, "SO spin-mixing transfer (interleaved)")
        for R in (-1, 0, 1):
            mil = _spinblock_to_interleaved(H[R])
            for i in range(ND):
                for j in range(ND):
                    v = mil[i, j]
                    if v != 0:
                        f.write("{} 0 0 {} {} {:.12f} {:.12f}\n".format(
                            R, i + 1, j + 1, v.real, v.imag))

    # SO Transfer, spin-INDEPENDENT: used by the green_init route, where
    # the spin mixing arrives through the Green function instead, so the
    # SO and non-SO H0 assemblies must coincide before green_init is read
    with open(os.path.join(d, "transfer_so_indep.dat"), "w") as f:
        _w90_header(f, ND, 3, "SO spin-independent transfer (interleaved)")
        for R in (-1, 1):
            for a in range(NORB):
                for s in range(NS):
                    i = 2 * a + s + 1
                    f.write("{} 0 0 {} {} -0.500000 0.0\n".format(R, i, i))

    # physical (non-SO) Transfer: a spin-independent placeholder. On the
    # trans_mod route H0 is REPLACED by the npz, so only its structure
    # matters; on the green_init route it equals transfer_so_indep.
    with open(os.path.join(d, "transfer_phys.dat"), "w") as f:
        _w90_header(f, NORB, 3, "physical 2-orbital diagonal transfer")
        for R in (-1, 1):
            for a in range(NORB):
                f.write("{} 0 0 {} {} -0.500000 0.0\n".format(R, a + 1, a + 1))

    # interactions: PHYSICAL orbital indices in both modes
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        _w90_header(f, NORB, 1, "CoulombIntra")
        for a in range(NORB):
            f.write("0 0 0 {} {} 1.000000 0.0\n".format(a + 1, a + 1))

    with open(os.path.join(d, "coulombinter_onsite.dat"), "w") as f:
        _w90_header(f, NORB, 1, "CoulombInter, on-site cross-orbital")
        f.write("0 0 0 1 2 0.350000 0.0\n")
        f.write("0 0 0 2 1 0.350000 0.0\n")

    with open(os.path.join(d, "pairlift.dat"), "w") as f:
        _w90_header(f, NORB, 1, "PairLift on-site cross-orbital")
        f.write("0 0 0 1 2 0.300000 0.0\n")
        f.write("0 0 0 2 1 0.300000 0.0\n")

    # trans_mod, in both conventions
    tab_sb = np.zeros((LX, ND, ND), dtype=complex)
    tab_il = np.zeros((LX, ND, ND), dtype=complex)
    for R, m in H.items():
        tab_sb[R % LX] += m
        tab_il[R % LX] += _spinblock_to_interleaved(m)
    np.savez(os.path.join(d, "trans_mod_spinblock.npz"), trans_mod=tab_sb)
    np.savez(os.path.join(d, "trans_mod_interleaved.npz"), trans_mod=tab_il)

    # green_init, in both conventions: spin-mixing on-site block
    g = np.zeros((LX, ND, ND), dtype=complex)
    g[0][0:2, 0:2] = [[0.55, 0.05], [0.05, 0.45]]
    g[0][2:4, 2:4] = [[0.40, 0.03], [0.03, 0.50]]
    gm = np.array([[0.12 + 0.04j, 0.02], [0.01, 0.09 - 0.03j]], dtype=complex)
    g[0][0:2, 2:4] = gm
    g[0][2:4, 0:2] = gm.conj().T
    g[1][0:2, 0:2] = 0.05 * np.eye(2)
    g[1][2:4, 2:4] = 0.05 * np.eye(2)
    g[LX - 1] = g[1].conj().T
    np.savez(os.path.join(d, "green_spinblock.npz"), green=g)
    gi = np.zeros_like(g)
    for r in range(LX):
        gi[r] = _spinblock_to_interleaved(g[r])
    np.savez(os.path.join(d, "green_interleaved.npz"), green=gi)


def run_route(d, *, so, interactions, transfer, trans_mod=None,
              green_init=None, calc_scheme="general", extra_param=None):
    """Run one RPA solve through the public input path.

    ``so`` selects the spin-orbital route (SO Geometry/Transfer plus the
    interleaved npz) or the issue-#174 route (physical Geometry/Transfer
    plus the spin-block npz).
    """
    inter = {"path_to_input": d,
             "Geometry": "geom_so.dat" if so else "geom_phys.dat",
             "Transfer": transfer[so]}
    inter.update(interactions)
    input_block = {"path_to_input": d, "interaction": inter}
    if trans_mod is not None:
        input_block["trans_mod"] = trans_mod[so]
    if green_init is not None:
        input_block["green_init"] = green_init[so]

    param = {"T": T, "filling": FILLING, "CellShape": [LX, 1, 1],
             "SubShape": [1, 1, 1], "Nmat": NMAT}
    param.update(extra_param or {})
    info_mode = {"mode": "RPA", "param": param,
                 "enable_spin_orbital": so, "calc_scheme": calc_scheme}

    out = os.path.join(d, "out")
    os.makedirs(out, exist_ok=True)
    rio = read_input_k.QLMSkInput(input_block)
    solver = solver_rpa.RPA(rio.get_param("ham"), {}, info_mode)
    green = rio.get_param("green")
    green.update(solver.read_init(input_block))
    solver.solve(green, out)
    return solver, green


TRANSFER_MIX = {True: "transfer_so_mix.dat", False: "transfer_phys.dat"}
TRANSFER_INDEP = {True: "transfer_so_indep.dat", False: "transfer_phys.dat"}
TRANS_MOD = {True: "trans_mod_interleaved.npz",
             False: "trans_mod_spinblock.npz"}
GREEN_INIT = {True: "green_interleaved.npz", False: "green_spinblock.npz"}


class TestSpinfulNpzGeneralVertex(unittest.TestCase):
    """Issue #174: under ``calc_scheme='general'`` the npz routes must
    resum the SAME antisymmetrized vertex as the adjudicated
    spin-orbital run.

    RED evidence on the unfixed code (develop @ bbb8637f): ``chiq``
    differed by max 1.705e-02 (12.2% relative) on the ``trans_mod``
    route and 1.707e-02 (12.3%) on the ``green_init`` route, while
    ``H0_k`` / ``chi0q`` / ``ham_inter_q`` were already bit-identical --
    isolating the difference to the missing crossing.
    """

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name
        write_fixture(self.d)

    def _assert_routes_agree(self, sA, gA, sB, gB):
        # both routes must really be on the spinful branch, else the
        # comparison would be vacuous
        self.assertEqual(sA.spin_mode, "spinful")
        self.assertEqual(sB.spin_mode, "spinful")
        self.assertTrue(sA.ham_info.enable_spin_orbital)
        self.assertFalse(sB.ham_info.enable_spin_orbital)
        # the inputs to the solve are identical...
        np.testing.assert_array_equal(np.asarray(sA.H0_k),
                                      np.asarray(sB.H0_k))
        np.testing.assert_array_equal(np.asarray(gA["chi0q"]),
                                      np.asarray(gB["chi0q"]))
        np.testing.assert_array_equal(
            np.asarray(sA.ham_info.ham_inter_q),
            np.asarray(sB.ham_info.ham_inter_q))
        # ... so the OUTPUT must be too, BITWISE. Asserted before the
        # structural checks below so that a regression reports the
        # physical divergence (1.7e-02) rather than a None tensor.
        np.testing.assert_array_equal(np.asarray(gA["chiq"]),
                                      np.asarray(gB["chiq"]))
        # the spinful solve MATERIALIZED the crossing on both routes --
        # reading the cache attribute (not calling the builder) is what
        # makes this a pin on consumption, and it must be identical
        self.assertIsNotNone(sA.ham_info.ham_spinful_exchange)
        self.assertIsNotNone(sB.ham_info.ham_spinful_exchange)
        np.testing.assert_array_equal(
            np.asarray(sA.ham_info.ham_spinful_exchange),
            np.asarray(sB.ham_info.ham_spinful_exchange))

    def test_trans_mod_route_matches_spin_orbital_bitwise(self):
        """CoulombIntra, spin-mixing H0 via ``trans_mod``: the
        ``enable_spin_orbital=false`` run must reproduce the
        ``enable_spin_orbital=true`` run bitwise under ``general``."""
        inter = {"CoulombIntra": "coulombintra.dat"}
        sA, gA = run_route(self.d, so=True, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        sB, gB = run_route(self.d, so=False, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        self._assert_routes_agree(sA, gA, sB, gB)

    def test_trans_mod_route_pairlift_matches_spin_orbital_bitwise(self):
        """PairLift + CoulombIntra: the sweep's largest divergence
        (12.3% relative before the fix). PairLift's ring spin table has
        spin-off-diagonal slots, so it exercises a different sourcing
        branch of the crossing than CoulombIntra alone."""
        inter = {"CoulombIntra": "coulombintra.dat",
                 "PairLift": "pairlift.dat"}
        sA, gA = run_route(self.d, so=True, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        sB, gB = run_route(self.d, so=False, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        self._assert_routes_agree(sA, gA, sB, gB)

    def test_green_init_route_matches_spin_orbital_bitwise(self):
        """The ``green_init`` half of the route: ``_calc_epsilon_k`` sets
        the spin-orbital branch for it exactly as for ``trans_mod``, and
        the Fock term of a spin-mixing Green function makes H0 spinful.
        The transfer here is spin-INDEPENDENT and identical in the two
        index conventions, so the assemblies coincide before the npz."""
        inter = {"CoulombIntra": "coulombintra.dat",
                 "PairLift": "pairlift.dat"}
        sA, gA = run_route(self.d, so=True, interactions=inter,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT)
        sB, gB = run_route(self.d, so=False, interactions=inter,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT)
        self._assert_routes_agree(sA, gA, sB, gB)

    def test_auto_scheme_on_npz_route_matches_spin_orbital(self):
        """``calc_scheme='auto'`` promotes the ``trans_mod`` route to
        ``general`` (``auto:mixed:trans_mod``, #167), i.e. the DEFAULT
        configuration chose the defective computation. Pin that the
        promoted run now agrees with the adjudicated one."""
        inter = {"CoulombIntra": "coulombintra.dat",
                 "CoulombInter": "coulombinter_onsite.dat"}
        sA, gA = run_route(self.d, so=True, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        sB, gB = run_route(self.d, so=False, interactions=inter,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD,
                           calc_scheme="auto")
        self.assertEqual(sB.calc_scheme, "general")
        self.assertEqual(sB._scheme_resolution, "auto:mixed:trans_mod")
        self._assert_routes_agree(sA, gA, sB, gB)
    def test_auto_scheme_on_green_init_route_matches_spin_orbital(self):
        """The ``green_init`` analogue of the previous test: ``auto``
        resolves to ``general`` via ``auto:mixed:green_init``.

        PairLift is what makes this fixture spinful at all -- its ring
        spin table has the spin-off-diagonal slots the green_init Fock
        term needs to produce a spin-MIXING H0 -- and it is also a
        conditional type, so it drives the promotion too."""
        inter = {"CoulombIntra": "coulombintra.dat",
                 "PairLift": "pairlift.dat"}
        sA, gA = run_route(self.d, so=True, interactions=inter,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT)
        sB, gB = run_route(self.d, so=False, interactions=inter,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT,
                           calc_scheme="auto")
        self.assertEqual(sB.calc_scheme, "general")
        self.assertEqual(sB._scheme_resolution, "auto:mixed:green_init")
        self._assert_routes_agree(sA, gA, sB, gB)


class TestNonSpinfulRunsUnchanged(unittest.TestCase):
    """No-change pin for issue #174's one-line fix.

    The crossing is now built for every run, but only the
    ``spin_mode == 'spinful'`` branch of ``RPA._solve_impl`` reads
    ``ham_spinful_exchange``; the ``spin-free`` and ``spin-diag``
    branches take ``_fierz_long()`` unconditionally. So building it
    cannot move a non-spinful result.
    """

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name
        write_fixture(self.d)

    def test_spin_free_general_run_is_bitwise_independent_of_exchange(self):
        """A plain spin-free (no npz, spin-independent transfer)
        ``general`` run: ``chiq`` is BITWISE identical with the crossing
        enabled and disabled, and ``spinful_vertex_exchange=false`` is
        exactly the pre-fix code path (the crossing is not consumed), so
        this equality demonstrates the fix changed nothing here."""
        inter = {"CoulombIntra": "coulombintra.dat"}
        s_on, g_on = run_route(self.d, so=False, interactions=inter,
                               transfer=TRANSFER_INDEP)
        s_off, g_off = run_route(
            self.d, so=False, interactions=inter, transfer=TRANSFER_INDEP,
            extra_param={"spinful_vertex_exchange": False})
        self.assertEqual(s_on.spin_mode, "spin-free")
        np.testing.assert_array_equal(np.asarray(g_on["chiq"]),
                                      np.asarray(g_off["chiq"]))
        # the solve never materialized the crossing
        self.assertIsNone(s_on.ham_info.ham_spinful_exchange)
        self.assertFalse(s_on.ham_info._spinful_exchange_built)
        # ... though this input does source one, on demand
        self.assertIsNotNone(s_on.ham_info.spinful_exchange())

    def test_spin_diag_general_run_is_bitwise_independent_of_exchange(self):
        """Same for a spin-DIAGONAL H0 supplied through ``trans_mod``
        (spin-dependent but not spin-mixing): ``spin_mode`` is
        ``'spin-diag'``, which also never reads the crossing."""
        diag = np.zeros((LX, ND, ND), dtype=complex)
        diag[0] = np.diag([0.30, -0.20, -0.10, 0.25])
        diag[1] = np.diag([-0.50, -0.45, -0.40, -0.35])
        diag[LX - 1] = diag[1].conj().T
        np.savez(os.path.join(self.d, "trans_mod_diag.npz"), trans_mod=diag)
        tm = {False: "trans_mod_diag.npz"}
        inter = {"CoulombIntra": "coulombintra.dat"}
        s_on, g_on = run_route(self.d, so=False, interactions=inter,
                               transfer=TRANSFER_INDEP, trans_mod=tm)
        s_off, g_off = run_route(
            self.d, so=False, interactions=inter, transfer=TRANSFER_INDEP,
            trans_mod=tm, extra_param={"spinful_vertex_exchange": False})
        self.assertEqual(s_on.spin_mode, "spin-diag")
        np.testing.assert_array_equal(np.asarray(g_on["chiq"]),
                                      np.asarray(g_off["chiq"]))
        self.assertIsNone(s_on.ham_info.ham_spinful_exchange)
        self.assertFalse(s_on.ham_info._spinful_exchange_built)
        self.assertIsNotNone(s_on.ham_info.spinful_exchange())


class TestExchangeIsNotMaterializedEagerly(unittest.TestCase):
    """Memory regression pin (issue #174 round-2 review).

    The crossing accumulator is ``(ns*norb)^4`` complex128, and
    ``ham_spinful_exchange`` holds a transposed VIEW of it, so the base
    allocation stays alive as long as the attribute does. Building it in
    the constructor would add about 1.6 GB at ``norb = 50`` -- roughly
    doubling the interaction-tensor footprint of a one-cell calculation
    and able to turn a previously working non-spinful input into an OOM.

    ``norb = 12`` here: ``(2*12)^4 * 16 B`` is about 5.3 MB, large enough
    that eager construction is a real cost and cheap enough for the fast
    gate. The pin asserts the named cache, which is non-fragile: it is
    the only reference the solver keeps to that allocation.
    """

    NORB_BIG = 12
    LX = 4

    def _build(self, d):
        nd = 2 * self.NORB_BIG
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n"
                    % self.NORB_BIG)
            for _ in range(self.NORB_BIG):
                f.write("0.0 0.0 0.0\n")
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("hdr\n%d\n2\n1 1\n" % self.NORB_BIG)
            for a in range(1, self.NORB_BIG + 1):
                f.write(" 1 0 0 %d %d -0.500000 0.0\n" % (a, a))
                f.write("-1 0 0 %d %d -0.500000 0.0\n" % (a, a))
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("hdr\n%d\n1\n1\n" % self.NORB_BIG)
            for a in range(1, self.NORB_BIG + 1):
                f.write(" 0 0 0 %d %d 1.000000 0.0\n" % (a, a))
        inter = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat",
                 "CoulombIntra": "coulombintra.dat"}
        info_mode = {"mode": "RPA",
                     "param": {"T": T, "filling": FILLING,
                               "CellShape": [self.LX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 8},
                     "enable_spin_orbital": False,
                     "calc_scheme": "general"}
        rio = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": inter})
        solver = solver_rpa.RPA(rio.get_param("ham"), {}, info_mode)
        return solver, rio, nd

    def test_cache_is_empty_after_construction_and_after_spin_free_solve(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        d = tmp.name
        solver, rio, nd = self._build(d)

        # after CONSTRUCTION: nothing dense is held
        self.assertIsNone(solver.ham_info.ham_spinful_exchange)
        self.assertFalse(solver.ham_info._spinful_exchange_built)

        out = os.path.join(d, "out")
        os.makedirs(out, exist_ok=True)
        green = rio.get_param("green")
        solver.solve(green, out)

        # after a SPIN-FREE solve: still nothing
        self.assertEqual(solver.spin_mode, "spin-free")
        self.assertIsNone(solver.ham_info.ham_spinful_exchange)
        self.assertFalse(solver.ham_info._spinful_exchange_built)

        # and the tensor this input would have allocated is the big one,
        # so the pin above is measuring something real
        X = solver.ham_info.spinful_exchange()
        self.assertIsNotNone(X)
        self.assertEqual(np.asarray(X).shape, (nd, nd, nd, nd))
        self.assertGreater(np.asarray(X).nbytes, 5_000_000)


if __name__ == "__main__":
    unittest.main()
