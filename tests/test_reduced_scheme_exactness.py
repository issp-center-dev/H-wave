#!/usr/bin/env python3
"""When is ``calc_scheme='reduced'`` EXACT, and when is it an approximation?

``reduced`` keeps only the density-pair slots of the vertex and of chi
(``density_projection.project_density_pairs``, ``'kaabb->kab'``), while
the adjudicated vertex table (``hwave.solver.vertex_table``) gives
CoulombInter, Hund and Ising a nonzero ``cross`` family in ADDITION to
their ``density`` family. So ``reduced`` discards content for those three
types, and reproduces the density-pair projection of ``general`` only
when the discarded sector cannot re-enter through the RPA series.

That condition holds when the effective one-body Hamiltonian conserves
the ORIGINAL orbital flavour: the discarded Fierz/cross entries live in
unequal-flavour sectors, which a flavour-diagonal H0 leaves unreachable
from a density-pair observable. It fails once the flavours hybridise.

These tests record that boundary as it stands today. They compare the two
schemes DIRECTLY (both requested explicitly), so they say nothing about
which scheme ``auto`` picks -- that is pinned separately below, because
``auto`` currently selects ``reduced`` for exactly the hybridised
multi-orbital case where it is NOT exact. Changing that is a behavioural
and output-shape change deferred to its own design round; these tests are
the honest record of the present state, and a deliberate prompt to update
them when it happens.
"""

import os
import shutil
import tempfile
import unittest

import numpy as np

_SRC = "tests/rpa/input_2orb"

#: A transfer that is diagonal in the orbital index: +-y hops on each
#: orbital, nothing connecting orbital 1 to orbital 2.
_DIAGONAL_TRANSFER = """Transfer, orbital-diagonal only (no hybridisation)
2
4
 1 1 1 1
   0    1    0    1    1  1.0 0.0
   0   -1    0    1    1  1.0 0.0
   0    1    0    2    2  1.0 0.0
   0   -1    0    2    2  1.0 0.0
"""


def _project_density_pairs(chiq_general):
    """The density-pair slots of a general chiq, in the reduced layout.

    ``bubble.contract_general`` stores the trailing axes as ``(a, c, b, d)``
    holding ``G_fwd[a, b] * G_rev[d, c]``; selecting ``[a, a, b, b]`` gives
    ``G_fwd[a, b] * G_rev[b, a]``, which is exactly what
    ``bubble.contract_reduced`` computes. The spin inflation that follows
    uses the same spin-major generalized index as the reduced output, so
    the two are directly comparable.
    """
    return np.einsum("kqaabb->kqab", chiq_general)


class _SchemeComparison(unittest.TestCase):
    def _solve(self, scheme, interactions, hybridised, subshape=(1, 1, 1)):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        d = tempfile.mkdtemp(prefix="reduced_exact_")
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        shutil.copy(os.path.join(_SRC, "geom.dat"), d)
        if hybridised:
            shutil.copy(os.path.join(_SRC, "transfer.dat"), d)
        else:
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write(_DIAGONAL_TRANSFER)

        idict = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat"}
        for key, fname in interactions.items():
            shutil.copy(os.path.join(_SRC, fname), d)
            idict[key] = fname

        reader = read_input_k.QLMSkInput(
            {"path_to_input": d, "interaction": idict})
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
               "SubShape": list(subshape), "Nmat": 16}
        solver = rpa_mod.RPA(
            reader.get_param("ham"), {},
            {"mode": "RPA", "param": par, "enable_spin_orbital": False,
             "calc_scheme": scheme, "calc_type": "ring"})
        green_info = reader.get_param("green")
        out = tempfile.mkdtemp(prefix="reduced_exact_out_")
        self.addCleanup(shutil.rmtree, out, ignore_errors=True)
        solver.solve(green_info, out)
        return solver, np.asarray(green_info["chiq"])


class TestReducedIsExactWhenFlavourIsConserved(_SchemeComparison):
    """A flavour-diagonal H0 leaves the discarded sector unreachable."""

    CASES = {
        "CoulombInter": {"CoulombInter": "onsite_inter.dat"},
        "Hund": {"Hund": "hund_onsite.dat"},
    }

    def test_orbital_diagonal_transfer_gives_exact_agreement(self):
        for name, interactions in self.CASES.items():
            with self.subTest(interaction=name):
                _, general = self._solve("general", interactions, False)
                _, reduced = self._solve("reduced", interactions, False)
                projected = _project_density_pairs(general)
                self.assertEqual(projected.shape, reduced.shape)
                self.assertTrue(
                    np.array_equal(projected, reduced),
                    "reduced should reproduce the density-pair projection "
                    "of general exactly for a flavour-diagonal H0; "
                    "max|diff| = {}".format(
                        float(np.max(np.abs(projected - reduced)))))

    def test_folding_alone_does_not_break_the_agreement(self):
        """Folding mixes SUBLATTICE labels but preserves the original
        orbital flavour, so the discarded unequal-flavour sector stays
        unreachable. (The folded chi0 itself is NOT density-pair block
        diagonal in the supercell index -- a folded one-orbital chain has
        G[0,1] != 0 and hence a nonzero P*chi0*Q -- so the exactness comes
        from flavour conservation, not from any closure of the supercell
        density subspace.)
        """
        interactions = {"CoulombInter": "onsite_inter.dat"}
        _, general = self._solve("general", interactions, False,
                                 subshape=(2, 1, 1))
        _, reduced = self._solve("reduced", interactions, False,
                                 subshape=(2, 1, 1))
        projected = _project_density_pairs(general)
        self.assertTrue(
            np.array_equal(projected, reduced),
            "folding an orbital-diagonal transfer should not break "
            "reduced's exactness; max|diff| = {}".format(
                float(np.max(np.abs(projected - reduced)))))


class TestReducedIsApproximateUnderHybridisation(_SchemeComparison):
    """Inter-orbital hopping makes the discarded sector reachable.

    The magnitudes recorded here are for THIS fixture only. Near an RPA
    instability the resolvent is ill-conditioned and the departure can be
    amplified without a useful upper bound, so these numbers are evidence
    that the effect is structural -- not a bound on its size.
    """

    CASES = {
        "CoulombInter": ({"CoulombInter": "onsite_inter.dat"}, 2.6e-05),
        "Hund": ({"Hund": "hund_onsite.dat"}, 3.8e-05),
    }

    def test_inter_orbital_hopping_makes_the_schemes_differ(self):
        for name, (interactions, recorded) in self.CASES.items():
            with self.subTest(interaction=name):
                _, general = self._solve("general", interactions, True)
                _, reduced = self._solve("reduced", interactions, True)
                diff = float(np.max(np.abs(
                    _project_density_pairs(general) - reduced)))
                self.assertGreater(
                    diff, 1e-8,
                    "the schemes are expected to DIFFER here: reduced "
                    "discards the cross-family vertex content that "
                    "hybridisation makes reachable")
                # Loose bracket around the recorded value: this pins the
                # effect's scale without turning a tolerance drift into a
                # false alarm.
                self.assertLess(diff, 10.0 * recorded)
                self.assertGreater(diff, 0.1 * recorded)


class TestAutoCurrentlySelectsTheApproximation(_SchemeComparison):
    """The present behaviour of ``calc_scheme='auto'``, pinned so that
    changing it is a deliberate act.

    ``auto`` forces ``general`` only for Exchange and PairHop, which carry
    no density-family content at all. For CoulombInter/Hund/Ising it picks
    ``reduced`` regardless of hybridisation -- i.e. exactly the case the
    class above shows is approximate. Correcting this changes chiq from
    four axes to six for affected inputs, so it is deferred to its own
    design round; update this test when that lands.
    """

    def test_auto_picks_reduced_for_a_hybridised_coulombinter_model(self):
        solver, chiq = self._solve(
            "auto", {"CoulombInter": "onsite_inter.dat"}, True)
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertEqual(chiq.ndim, 4)


if __name__ == "__main__":
    unittest.main()
