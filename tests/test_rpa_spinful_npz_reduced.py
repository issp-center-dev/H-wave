#!/usr/bin/env python3

"""Issue #161's ADJUDICATION RECORD: ``calc_scheme = 'reduced'`` with a
spinful one-body Hamiltonian supplied through ``trans_mod`` / ``green_init``.

THE QUESTION (#161). ``rpa.py`` decides ``spin_mode`` from the block
structure of H0(k), independently of ``enable_spin_orbital``: the npz
routes set the spin-orbital branch in ``_calc_epsilon_k``
unconditionally, so a user-supplied spin-mixing ``trans_mod.npz`` /
``green_init.npz`` reaches ``spin_mode == 'spinful'`` with
``enable_spin_orbital = False``, and can be combined with the reduced
scheme -- a state reachable from a public input path. FLEX now rejects
that combination (``FLEX._check_reduced_rejects_spinful``); whether
RPA's own ``chi0q`` / ``chiq`` there are meaningful had not been
adjudicated, and the RPA/FLEX equivalence table has no cell for the
route (its fixture model is TOML plus text inputs).

THE ANSWER RECORDED HERE. The npz route is not a separate computation:
on identical physics it produces BITWISE the same ``chi0q`` and
``chiq`` as the ``enable_spin_orbital = true`` run, which IS
adjudicated -- ``chi0`` against the exact Lindhard function (#133) and
the vertex against exact diagonalization (#137, in the pair convention
fixed by #139). The arrays are identical, so the adjudication transfers
verbatim: nothing about the reduced scheme's meaning changes with how
the spinful H0 arrived. No guard is warranted; the approximation the
reduced scheme makes on a flavour-mixing H0 is the ORDINARY reduced
approximation, and #167's exactness diagnostic already announces it on
this very route (pinned below, with cause token ``trans_mod``).

Route A : ``enable_spin_orbital=True``, SO Geometry/Transfer, interleaved
          npz -- the adjudicated spin-orbital computation.
Route B : ``enable_spin_orbital=False``, physical Geometry/Transfer,
          spin-block npz -- the issue-#161 state.

The ``general`` half of the same comparison is issue #174 and lives in
``tests/test_rpa_spinful_npz_general.py``; there the two routes did NOT
agree until the exchange crossing was built unconditionally. Reduced was
never affected: the crossing has no density-diagonal content, so it
cannot survive ``project_density_pairs``.

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
    """H(R) in SPIN-BLOCK order (index = s*NORB + a) for R = 0, +1, -1,
    with a genuinely spin-MIXING R = 0 sector."""
    H = {}
    h0 = np.zeros((ND, ND), dtype=complex)
    h0[0:2, 0:2] = [[0.30, 0.10], [0.10, -0.20]]
    h0[2:4, 2:4] = [[-0.10, 0.05], [0.05, 0.25]]
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
    H = _h_blocks()

    for name, n in (("geom_so.dat", ND), ("geom_phys.dat", NORB)):
        with open(os.path.join(d, name), "w") as f:
            f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n")
            f.write("{}\n".format(n))
            for _ in range(n):
                f.write("  0.0 0.0 0.0\n")

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

    with open(os.path.join(d, "transfer_so_indep.dat"), "w") as f:
        _w90_header(f, ND, 3, "SO spin-independent transfer (interleaved)")
        for R in (-1, 1):
            for a in range(NORB):
                for s in range(NS):
                    i = 2 * a + s + 1
                    f.write("{} 0 0 {} {} -0.500000 0.0\n".format(R, i, i))

    with open(os.path.join(d, "transfer_phys.dat"), "w") as f:
        _w90_header(f, NORB, 3, "physical 2-orbital diagonal transfer")
        for R in (-1, 1):
            for a in range(NORB):
                f.write("{} 0 0 {} {} -0.500000 0.0\n".format(R, a + 1, a + 1))

    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        _w90_header(f, NORB, 1, "CoulombIntra")
        for a in range(NORB):
            f.write("0 0 0 {} {} 1.000000 0.0\n".format(a + 1, a + 1))

    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        _w90_header(f, NORB, 3, "CoulombInter (off-site, physical indices)")
        for R in (-1, 1):
            for a in range(NORB):
                for b in range(NORB):
                    f.write("{} 0 0 {} {} 0.400000 0.0\n".format(
                        R, a + 1, b + 1))

    with open(os.path.join(d, "pairlift.dat"), "w") as f:
        _w90_header(f, NORB, 1, "PairLift on-site cross-orbital")
        f.write("0 0 0 1 2 0.300000 0.0\n")
        f.write("0 0 0 2 1 0.300000 0.0\n")

    tab_sb = np.zeros((LX, ND, ND), dtype=complex)
    tab_il = np.zeros((LX, ND, ND), dtype=complex)
    for R, m in H.items():
        tab_sb[R % LX] += m
        tab_il[R % LX] += _spinblock_to_interleaved(m)
    np.savez(os.path.join(d, "trans_mod_spinblock.npz"), trans_mod=tab_sb)
    np.savez(os.path.join(d, "trans_mod_interleaved.npz"), trans_mod=tab_il)

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


TRANSFER_MIX = {True: "transfer_so_mix.dat", False: "transfer_phys.dat"}
TRANSFER_INDEP = {True: "transfer_so_indep.dat", False: "transfer_phys.dat"}
TRANS_MOD = {True: "trans_mod_interleaved.npz",
             False: "trans_mod_spinblock.npz"}
GREEN_INIT = {True: "green_interleaved.npz", False: "green_spinblock.npz"}


def run_route(d, *, so, interactions, transfer, trans_mod=None,
              green_init=None, calc_scheme="reduced"):
    """Run one RPA solve through the public input path."""
    inter = {"path_to_input": d,
             "Geometry": "geom_so.dat" if so else "geom_phys.dat",
             "Transfer": transfer[so]}
    inter.update(interactions)
    input_block = {"path_to_input": d, "interaction": inter}
    if trans_mod is not None:
        input_block["trans_mod"] = trans_mod[so]
    if green_init is not None:
        input_block["green_init"] = green_init[so]

    info_mode = {"mode": "RPA",
                 "param": {"T": T, "filling": FILLING,
                           "CellShape": [LX, 1, 1], "SubShape": [1, 1, 1],
                           "Nmat": NMAT},
                 "enable_spin_orbital": so, "calc_scheme": calc_scheme}
    out = os.path.join(d, "out")
    os.makedirs(out, exist_ok=True)
    rio = read_input_k.QLMSkInput(input_block)
    solver = solver_rpa.RPA(rio.get_param("ham"), {}, info_mode)
    green = rio.get_param("green")
    green.update(solver.read_init(input_block))
    solver.solve(green, out)
    return solver, green


class TestReducedSpinfulNpzMatchesSpinOrbital(unittest.TestCase):
    """Issue #161's adjudication record, part 1: the arrays are the same.

    Explicit ``calc_scheme='reduced'`` with a spinful H0 arriving through
    ``trans_mod`` / ``green_init`` and ``enable_spin_orbital = false``
    computes BITWISE the same ``chi0q`` and ``chiq`` as the
    ``enable_spin_orbital = true`` run on identical physics. Because the
    arrays are identical, the npz route INHERITS the spin-orbital mode's
    adjudication -- ``chi0`` = exact Lindhard (#133) and the #137/#139
    vertex work -- and needs no adjudication of its own, and no guard.
    (The investigation measured this equality across nine interaction
    sets; the three below plus the ``green_init`` case are the record
    kept in the suite.)
    """

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name
        write_fixture(self.d)

    def _assert_bitwise(self, sA, gA, sB, gB):
        # both routes must really be spinful, else the record is vacuous
        self.assertEqual(sA.spin_mode, "spinful")
        self.assertEqual(sB.spin_mode, "spinful")
        self.assertTrue(sA.ham_info.enable_spin_orbital)
        self.assertFalse(sB.ham_info.enable_spin_orbital)
        self.assertEqual(sA.calc_scheme, "reduced")
        self.assertEqual(sB.calc_scheme, "reduced")
        np.testing.assert_array_equal(np.asarray(sA.H0_k),
                                      np.asarray(sB.H0_k))
        np.testing.assert_array_equal(np.asarray(gA["chi0q"]),
                                      np.asarray(gB["chi0q"]))
        np.testing.assert_array_equal(np.asarray(gA["chiq"]),
                                      np.asarray(gB["chiq"]))

    def _trans_mod_case(self, interactions):
        sA, gA = run_route(self.d, so=True, interactions=interactions,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        sB, gB = run_route(self.d, so=False, interactions=interactions,
                           transfer=TRANSFER_MIX, trans_mod=TRANS_MOD)
        self._assert_bitwise(sA, gA, sB, gB)

    def test_trans_mod_coulombintra(self):
        """CoulombIntra: the unconditional type, exact under reduced."""
        self._trans_mod_case({"CoulombIntra": "coulombintra.dat"})

    def test_trans_mod_coulombintra_plus_coulombinter(self):
        """CoulombIntra + off-site CoulombInter: CoulombInter is a
        CONDITIONAL type, so this is the case where reduced genuinely
        discards cross-family vertex sectors -- and still both routes
        discard exactly the same ones."""
        self._trans_mod_case({"CoulombIntra": "coulombintra.dat",
                              "CoulombInter": "coulombinter.dat"})

    def test_trans_mod_pairlift_plus_coulombintra(self):
        """PairLift + CoulombIntra: PairLift's ring spin table has
        spin-off-diagonal slots, so it is the interaction most likely to
        expose a spin-index convention mismatch between the interleaved
        and spin-block readings. It does not."""
        self._trans_mod_case({"CoulombIntra": "coulombintra.dat",
                              "PairLift": "pairlift.dat"})

    def test_green_init_route(self):
        """The ``green_init`` half of the route. Here the transfer is
        spin-INDEPENDENT and identical in both index conventions, and
        the spin mixing enters through the initial Green function's Fock
        term -- a different assembly path to the same conclusion."""
        interactions = {"CoulombIntra": "coulombintra.dat",
                        "PairLift": "pairlift.dat"}
        sA, gA = run_route(self.d, so=True, interactions=interactions,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT)
        sB, gB = run_route(self.d, so=False, interactions=interactions,
                           transfer=TRANSFER_INDEP, green_init=GREEN_INIT)
        self._assert_bitwise(sA, gA, sB, gB)


class TestReducedSpinfulNpzExactnessWarning(unittest.TestCase):
    """Issue #161's adjudication record, part 2: the user is told.

    Part of the answer to "is this state acceptable?" is that it is not
    silent. #167's reduced-exactness diagnostic fires on the npz route,
    naming ``trans_mod`` / ``green_init`` as the cause of the flavour
    mixing, so a user who explicitly asks for ``reduced`` with a
    conditional type and a spinful npz is told the scheme is an
    approximation for that input and that ``general`` / ``auto`` is
    exact. No additional guard is warranted on top of it.
    """

    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.d = tmp.name
        write_fixture(self.d)

    def test_warning_fires_with_trans_mod_cause(self):
        interactions = {"CoulombIntra": "coulombintra.dat",
                        "CoulombInter": "coulombinter.dat"}
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as ctx:
            solver, _ = run_route(self.d, so=False,
                                  interactions=interactions,
                                  transfer=TRANSFER_MIX,
                                  trans_mod=TRANS_MOD)
        self.assertEqual(solver.spin_mode, "spinful")
        hits = [m for m in ctx.output
                if "mixes orbital flavour (trans_mod)" in m]
        self.assertTrue(hits, ctx.output)
        self.assertIn("calc_scheme='reduced' with CoulombInter", hits[0])
        self.assertIn("is an APPROXIMATION for this input", hits[0])

    def test_warning_fires_with_green_init_cause(self):
        interactions = {"CoulombIntra": "coulombintra.dat",
                        "PairLift": "pairlift.dat"}
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as ctx:
            solver, _ = run_route(self.d, so=False,
                                  interactions=interactions,
                                  transfer=TRANSFER_INDEP,
                                  green_init=GREEN_INIT)
        self.assertEqual(solver.spin_mode, "spinful")
        hits = [m for m in ctx.output
                if "mixes orbital flavour (green_init)" in m]
        self.assertTrue(hits, ctx.output)


if __name__ == "__main__":
    unittest.main()
