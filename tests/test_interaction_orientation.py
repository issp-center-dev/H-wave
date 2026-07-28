#!/usr/bin/env python3

"""The two paths that turn an interaction file into V(q) must agree on the
orbital orientation.

`rpa.py::Interaction._make_ham_inter` builds the vertex that the susceptibility
is solved with; `sc.py::_build_interaction_k` builds the one the pairing vertex
is assembled from. They read the same file, and they must produce the same
matrix -- not its orbital transpose.

The convention both must follow is fixed by Mizuno, Kobayashi and Suzumura,
"Role of charge fluctuation in Q1D organic superconductor (TMTSF)2ClO4",
arXiv:1010.1084 -- the multi-site extended-Hubbard RPA this code generalizes:

  Eq. (4)   V_ab(q) a+_{k+q s a} a+_{k'-q s' b} a_{k' s' b} a_{k s a}
            -> orbital a is the one carrying momentum transfer +q
  Eq. (13)  chi0_ab(q) = -(T/N) sum_k G_ab(k+q) G_ba(k)
            -> the chi0 convention this code already uses
  Eq. (12)  chi^C = (I + 2 chi0 V + chi0 U)^-1 chi0, V untransposed

The on-site-only sources cited elsewhere (MYO cond-mat/0407094, Kuroki-Aoki
0902.3691) cannot pin this: their vertex matrices are symmetric in the orbital
pair, so the orientation is invisible there.

The fixture below has V_ab(R) != V_ba(R) -- an ordinary zigzag bond
arrangement, with every bond declared from both ends so the Hamiltonian stays
Hermitian. Without that asymmetry V(q) is symmetric and the test is vacuous,
which is asserted.
"""

import os
import tempfile
import unittest

import numpy as np

NX = 4          # >= 3, so that some q differs from -q
NORB = 2

# the same physical bond declared from both ends, twice over, with different
# strengths for the two orbital orderings
BONDS = {(+1, 0, 1): 1.0, (-1, 1, 0): 1.0,
         (+1, 1, 0): 0.3, (-1, 0, 1): 0.3}


def _write_fixture(d):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("  1.0 0.0 0.0\n  0.0 1.0 0.0\n  0.0 0.0 1.0\n2\n"
                "   0.0 0.0 0.0\n   0.5 0.0 0.0\n")
    rows = ["   0    0    0    1    1   0.0 0.0",
            "   1    0    0    1    1  -1.0 0.0",
            "  -1    0    0    1    1  -1.0 0.0",
            "   1    0    0    2    2  -1.0 0.0",
            "  -1    0    0    2    2  -1.0 0.0"]
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("t\n2\n9\n 1 1 1 1 1 1 1 1 1\n" + "\n".join(rows) + "\n")
    ci = ["   {}    0    0    {}    {}   {:.12f}   0.0".format(
              R, a + 1, b + 1, v)
          for (R, a, b), v in sorted(BONDS.items())]
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("CoulombInter\n2\n9\n 1 1 1 1 1 1 1 1 1\n"
                + "\n".join(ci) + "\n")


def _physical_v_of_q():
    """V_ab(q) = sum_R V_ab(R) exp(-i q R), straight from the declared bonds."""
    qs = 2.0 * np.pi * np.arange(NX) / NX
    V = np.zeros((NX, NORB, NORB), dtype=complex)
    for (R, a, b), v in BONDS.items():
        V[:, a, b] += v * np.exp(-1j * qs * R)
    return V


class TestInteractionOrbitalOrientation(unittest.TestCase):

    def _build_both(self, dirpath):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.sc as sc
        from hwave.solver.rpa import Interaction, Lattice

        interaction = {"path_to_input": dirpath,
                       "Geometry": "geom.dat", "Transfer": "transfer.dat",
                       "CoulombInter": "coulombinter.dat"}
        ham_param = read_input_k.QLMSkInput(
            {"path_to_input": dirpath, "interaction": interaction}
        ).get_param("ham")

        info_mode = {"mode": "RPA",
                     "param": {"T": 1.0, "mu": 0.0, "CellShape": [NX, 1, 1],
                               "SubShape": [1, 1, 1], "Nmat": 4},
                     "calc_scheme": "reduced"}
        lattice = Lattice(info_mode["param"])
        inter = Interaction(lattice, ham_param, info_mode)

        ns = 2
        nd = NORB * ns
        ham = np.einsum(
            "ksasatbtb->ksatb",
            inter.ham_inter_q.reshape(lattice.nvol, *(ns, NORB) * 4)
        ).reshape(lattice.nvol, nd, nd)
        rpa_V = ham[:, :NORB, :NORB]          # spin-up/up density block

        qs = 2.0 * np.pi * np.arange(NX) / NX
        _, _, interactions = sc._read_interaction_files(
            {"file": {"input": {"path_to_input": dirpath,
                                "interaction": interaction}}})
        ik = sc._build_interaction_k(qs, np.array([0.0]), np.array([0.0]),
                                     interactions, NORB)
        sc_V = ik["CoulombInter"].transpose(2, 3, 4, 0, 1).reshape(
            NX, NORB, NORB)
        return rpa_V, sc_V

    def test_the_fixture_can_tell_the_orientations_apart(self):
        V = _physical_v_of_q()
        self.assertFalse(np.allclose(V, V.transpose(0, 2, 1)),
                         "V(q) is orbital-symmetric; the test would be vacuous")
        self.assertTrue(np.allclose(V, V.conj().transpose(0, 2, 1)),
                        "V(q) must still be Hermitian")

    def test_both_paths_build_the_same_v_of_q(self):
        phys = _physical_v_of_q()
        with tempfile.TemporaryDirectory(prefix="hwave_i96_") as d:
            _write_fixture(d)
            rpa_V, sc_V = self._build_both(d)

        scale = np.max(np.abs(phys))
        np.testing.assert_allclose(
            sc_V, phys, rtol=0.0, atol=1e-12 * scale,
            err_msg="sc._build_interaction_k must build V_ab(q) with a at +q")
        np.testing.assert_allclose(
            rpa_V, phys, rtol=0.0, atol=1e-12 * scale,
            err_msg="rpa._make_ham_inter must build V_ab(q) with a at +q, "
                    "not its orbital transpose")

    def test_flex_reduced_path_sees_the_same_orientation(self):
        """The FLEX reduced/squashed path consumes the same ``ham_inter_q``.

        ``_inflate_chi0q_and_ham`` contracts it down to the density block, so
        the orientation propagates there too.  The existing exact FLEX
        index-order tests all use on-site interactions, where the transpose is
        invisible.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex

        phys = _physical_v_of_q()
        with tempfile.TemporaryDirectory(prefix="hwave_i96_flex_") as d:
            _write_fixture(d)
            interaction = {"path_to_input": d,
                           "Geometry": "geom.dat", "Transfer": "transfer.dat",
                           "CoulombInter": "coulombinter.dat"}
            ham_param = read_input_k.QLMSkInput(
                {"path_to_input": d, "interaction": interaction}
            ).get_param("ham")
            flex = solver_flex.FLEX(ham_param, {}, {
                "mode": "FLEX",
                "param": {"T": 1.0, "mu": 0.0, "CellShape": [NX, 1, 1],
                          "SubShape": [1, 1, 1], "Nmat": 4},
                "calc_scheme": "reduced"})
            flex.spin_mode = "spin-free"
            nvol = flex.lattice.nvol
            chi0_dummy = np.zeros((flex.nmat, nvol, NORB, NORB),
                                  dtype=complex)
            _, ham = flex._inflate_chi0q_and_ham(
                chi0_dummy, flex.ham_info.ham_inter_q)

        block = np.asarray(ham)[:, :NORB, :NORB]
        scale = np.max(np.abs(phys))
        np.testing.assert_allclose(
            block, phys, rtol=0.0, atol=1e-12 * scale,
            err_msg="the FLEX reduced vertex must carry V_ab(q), not V_ba(q)")

    def test_trans_mod_hartree_term_is_unchanged_by_the_orientation(self):
        """``ham_inter_r`` also feeds ``_calc_trans_mod`` (the ``green_init``
        path).  That contraction sums over ALL r, i.e. it only sees q = 0, and
        for a bond set declared from both ends ``sum_r V_ab(r) == sum_r
        V_ba(r)``.  So this fix must leave ``trans_mod`` alone -- asserted here
        rather than assumed, by running the same contraction on the interaction
        and on its orbital-pair transpose.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        from hwave.solver.rpa import Interaction, Lattice

        with tempfile.TemporaryDirectory(prefix="hwave_i96_tm_") as d:
            _write_fixture(d)
            interaction = {"path_to_input": d,
                           "Geometry": "geom.dat", "Transfer": "transfer.dat",
                           "CoulombInter": "coulombinter.dat"}
            ham_param = read_input_k.QLMSkInput(
                {"path_to_input": d, "interaction": interaction}
            ).get_param("ham")
            info_mode = {"mode": "RPA",
                         "param": {"T": 1.0, "mu": 0.0,
                                   "CellShape": [NX, 1, 1],
                                   "SubShape": [1, 1, 1], "Nmat": 4},
                         "calc_scheme": "reduced"}
            lattice = Lattice(info_mode["param"])
            inter = Interaction(lattice, ham_param, info_mode)

        ns, nvol = 2, lattice.nvol
        nd = NORB * ns
        ww = inter.ham_inter_r.reshape(nvol, nd, nd, nd, nd)
        # the orientation this fix chose vs. the one it replaced
        ww_swapped = ww.transpose(0, 3, 4, 1, 2)

        rng = np.random.default_rng(96)
        gg = (rng.standard_normal((nd, nd))
              + 1j * rng.standard_normal((nd, nd)))
        gg = gg + gg.conj().T          # a Hermitian Green function block

        def hartree(w):
            h1 = np.einsum('rbacd,cd->rab', w, gg)
            h2 = np.einsum('rcdab,dc->rab', w, gg)
            return np.sum(h1 + h2, axis=0) / 2

        np.testing.assert_allclose(
            hartree(ww), hartree(ww_swapped), rtol=0.0, atol=1e-12,
            err_msg="trans_mod's Hartree term must not depend on the "
                    "orbital orientation of the interaction")
        # ...and not vacuously: the two tensors really do differ
        self.assertGreater(np.max(np.abs(ww - ww_swapped)), 1e-6)

    def test_every_density_type_uses_the_same_orientation(self):
        """``_append_inter`` places CoulombIntra, CoulombInter, Hund, Ising,
        PairLift and Exchange through one shared statement, so they must all
        come out in the same orientation.  Checked per type on the same
        asymmetric bond set: the density block must be proportional to V(q) and
        not to its transpose.  The proportionality constant is left free
        because each type enters the spin/charge decomposition with its own
        weight -- only the orientation is under test here.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        from hwave.solver.rpa import Interaction, Lattice

        phys = _physical_v_of_q()
        for itype, fname in (("CoulombInter", "coulombinter.dat"),
                             ("Hund", "hund.dat"),
                             ("Ising", "ising.dat")):
            with self.subTest(interaction=itype):
                with tempfile.TemporaryDirectory(prefix="hwave_i96_t_") as d:
                    _write_fixture(d)
                    if fname != "coulombinter.dat":
                        with open(os.path.join(d, fname), "w") as fh:
                            fh.write(open(os.path.join(
                                d, "coulombinter.dat")).read().replace(
                                    "CoulombInter", itype, 1))
                    interaction = {"path_to_input": d,
                                   "Geometry": "geom.dat",
                                   "Transfer": "transfer.dat",
                                   itype: fname}
                    ham_param = read_input_k.QLMSkInput(
                        {"path_to_input": d, "interaction": interaction}
                    ).get_param("ham")
                    info_mode = {"mode": "RPA",
                                 "param": {"T": 1.0, "mu": 0.0,
                                           "CellShape": [NX, 1, 1],
                                           "SubShape": [1, 1, 1], "Nmat": 4},
                                 "calc_scheme": "reduced"}
                    lattice = Lattice(info_mode["param"])
                    inter = Interaction(lattice, ham_param, info_mode)

                ns = 2
                nd = NORB * ns
                ham = np.einsum(
                    "ksasatbtb->ksatb",
                    inter.ham_inter_q.reshape(lattice.nvol, *(ns, NORB) * 4)
                ).reshape(lattice.nvol, nd, nd)[:, :NORB, :NORB]

                self.assertGreater(np.max(np.abs(ham)), 1e-9,
                                   "{}: vertex is zero".format(itype))
                c = (np.vdot(phys, ham) / np.vdot(phys, phys))
                d_ok = np.max(np.abs(ham - c * phys))
                d_bad = np.max(np.abs(ham - c * phys.transpose(0, 2, 1)))
                self.assertLess(d_ok, 1e-10 * np.max(np.abs(ham)),
                                "{}: not proportional to V(q)".format(itype))
                self.assertGreater(d_bad, 1e-6,
                                   "{}: V(q) and V(q)^T coincide; the check "
                                   "is vacuous".format(itype))


if __name__ == "__main__":
    unittest.main()
