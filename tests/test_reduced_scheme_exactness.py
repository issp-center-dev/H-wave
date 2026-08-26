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
        "Ising": {"Ising": "hund_onsite.dat"},
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
        unreachable and reduced stays exact.

        The fold DIRECTION decides whether this test means anything.
        ``_DIAGONAL_TRANSFER`` hops along +-y only, so folding along x
        leaves the two sublattice copies dynamically disconnected and
        every cross-sublattice Green element zero (measured: 0.0 for an x
        fold against 4.6e-02 for a y fold). The agreement would then
        follow from the folded density subspace closing after all --
        exactly the explanation this test exists to rule out. Folding
        along y makes those elements nonzero, so the supercell density
        subspace is genuinely NOT closed and the exactness can only come
        from flavour conservation. That precondition is asserted, not
        assumed.
        """
        interactions = {"CoulombInter": "onsite_inter.dat"}
        solver, general = self._solve("general", interactions, False,
                                      subshape=(1, 2, 1))
        _, reduced = self._solve("reduced", interactions, False,
                                 subshape=(1, 2, 1))

        g0 = np.asarray(solver.green0)
        norb = solver.norb
        offdiag = max(
            float(np.max(np.abs(g0[..., a, b])))
            for a in range(norb) for b in range(norb) if a != b)
        self.assertGreater(
            offdiag, 1e-6,
            "this fold left the sublattices disconnected, so the test "
            "cannot distinguish flavour conservation from closure of the "
            "supercell density subspace")

        projected = _project_density_pairs(general)
        diff = float(np.max(np.abs(projected - reduced)))
        # Round-off, not bit-identity: folding changes the array shapes
        # and so the summation order. Measured 1.7e-18 -- thirteen orders
        # below the 2.3e-04 STRUCTURAL departure hybridisation produces,
        # which is the distinction this bound must make.
        self.assertLess(
            diff, 1e-12,
            "folding an orbital-diagonal transfer should not break "
            "reduced's exactness; max|diff| = {}".format(diff))


class TestReducedIsApproximateUnderHybridisation(_SchemeComparison):
    """Inter-orbital hopping makes the discarded sector reachable.

    The magnitudes recorded here are for THIS fixture only. Near an RPA
    instability the resolvent is ill-conditioned and the departure can be
    amplified without a useful upper bound, so these numbers are evidence
    that the effect is structural -- not a bound on its size.
    """

    #: Measured RELATIVE departure: max|proj(general) - reduced| over
    #: max|reduced|. Relative because that is how the finding is stated
    #: and it does not move under an overall rescaling of the fixture.
    CASES = {
        "CoulombInter": ({"CoulombInter": "onsite_inter.dat"}, 2.2572e-04),
        "Hund": ({"Hund": "hund_onsite.dat"}, 3.2763e-04),
        "Ising": ({"Ising": "hund_onsite.dat"}, 3.1681e-04),
    }

    def test_inter_orbital_hopping_makes_the_schemes_differ(self):
        for name, (interactions, recorded) in self.CASES.items():
            with self.subTest(interaction=name):
                _, general = self._solve("general", interactions, True)
                _, reduced = self._solve("reduced", interactions, True)
                scale = float(np.max(np.abs(reduced)))
                self.assertGreater(scale, 0.0)
                relative = float(np.max(np.abs(
                    _project_density_pairs(general) - reduced))) / scale
                self.assertGreater(
                    relative, 1e-8,
                    "the schemes are expected to DIFFER here: reduced "
                    "discards the cross-family vertex content that "
                    "hybridisation makes reachable")
                # Pin the SIZE, not merely the sign: a changed vertex
                # coefficient or normalisation moves this by a factor,
                # which a decade-wide bracket would wave through. 5% sits
                # far above the run-to-run floor (~1e-14, issue #85) and
                # far below any coefficient-scale regression.
                self.assertAlmostEqual(
                    relative / recorded, 1.0, delta=0.05,
                    msg="{}: relative departure {:.4e} is not the "
                        "recorded {:.4e}".format(name, relative, recorded))


class TestAutoCurrentlySelectsTheApproximation(_SchemeComparison):
    """The present behaviour of ``calc_scheme='auto'``, pinned so that
    changing it is a deliberate act.

    ``auto`` forces ``general`` only for Exchange and PairHop, which carry
    no density-family content at all. For CoulombInter/Hund/Ising it picks
    ``reduced`` regardless of hybridisation -- i.e. exactly the case the
    class above shows is approximate. That is issue #167: correcting it
    changes chiq from four axes to six for affected inputs, so it is
    deferred to its own design round; update this test when that lands.
    """

    def test_auto_picks_reduced_for_hybridised_models(self):
        for name, interactions in (
                ("CoulombInter", {"CoulombInter": "onsite_inter.dat"}),
                ("Hund", {"Hund": "hund_onsite.dat"}),
                ("Ising", {"Ising": "hund_onsite.dat"})):
            with self.subTest(interaction=name):
                solver, chiq = self._solve("auto", interactions, True)
                self.assertEqual(solver.calc_scheme, "reduced")
                self.assertEqual(chiq.ndim, 4)


if __name__ == "__main__":
    unittest.main()
