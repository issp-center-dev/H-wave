# tests/test_auto_scheme_exactness.py
"""calc_scheme='auto' is EXACT for the declared input (#167).

Public-entry tests: QLMSkInput -> RPA/FLEX -> solve. The fixture is the
2-orbital 4x4x1 model of tests/rpa/input_2orb (hybridised transfer.dat) or
its orbital-diagonal counterpart written inline.
"""
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

_SRC = "tests/rpa/input_2orb"

_DIAGONAL_TRANSFER = """Transfer, orbital-diagonal only (no hybridisation)
2
4
 1 1 1 1
   0    1    0    1    1  1.0 0.0
   0   -1    0    1    1  1.0 0.0
   0    1    0    2    2  1.0 0.0
   0   -1    0    2    2  1.0 0.0
"""

#: an on-site inter-orbital field (physical indices 1<->2), Hermitian-closed
_OFFDIAG_EXTERN = """Extern, inter-orbital
2
1
 1
   0    0    0    1    2  0.30 0.0
   0    0    0    2    1  0.30 0.0
"""

_DIAG_EXTERN = """Extern, orbital-diagonal
2
1
 1
   0    0    0    1    1  0.30 0.0
   0    0    0    2    2 -0.30 0.0
"""


#: A 1-orbital chain, for the sublattice-folding case: SubShape (2,1,1)
#: makes the POST-fold orbital count 2 while the PRE-FOLD physical count
#: stays 1, so the two candidate index limits disagree.
_CHAIN_GEOM = """  1.000000000000   0.000000000000   0.000000000000
  0.000000000000   1.000000000000   0.000000000000
  0.000000000000   0.000000000000   1.000000000000
1
    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00
"""

_CHAIN_TRANSFER = """Transfer, 1-orbital chain
1
2
 1 1
   1    0    0    1    1  -1.0 0.0
  -1    0    0    1    1  -1.0 0.0
"""

#: same-orbital nearest-neighbour V: a CONDITIONAL type, so the decision
#: turns on the flavour predicate rather than on general_only.
_CHAIN_COULOMBINTER = """CoulombInter, off-site same-orbital
1
2
 1 1
   1    0    0    1    1  0.3 0.0
  -1    0    0    1    1  0.3 0.0
"""

#: Extern declared in the spin-block-extended layout (num_wann = 2 =
#: 2 * norb_phys): index 2 is the SPIN block, which _make_ham_trans
#: discards ("skip spin dependence"). Judged with the pre-fold physical
#: count (1) these rows are out of range; judged with the post-fold count
#: (2) they look like an inter-orbital field and promote spuriously.
_CHAIN_SPIN_BLOCK_EXTERN = """Extern, spin-block extended
2
1
 1
   0    0    0    1    1  0.25 0.0
   0    0    0    1    2  0.20 0.0
   0    0    0    2    1  0.20 0.0
"""


def _project_density_pairs(chiq_general):
    return np.einsum("kqaabb->kqab", chiq_general)


class _Case(unittest.TestCase):
    """Builds a solver from files; ``interactions`` maps reader keys (ANY
    case) to fixture file names in tests/rpa/input_2orb or inline bodies."""

    INLINE = {"diag_transfer": _DIAGONAL_TRANSFER,
              "offdiag_extern": _OFFDIAG_EXTERN,
              "diag_extern": _DIAG_EXTERN,
              # tests/rpa/input_2orb has no Exchange/PairLift file: on-site,
              # Hermitian-closed inter-orbital bodies (physical indices 1<->2)
              "exchange": "Exchange\n2\n1\n 1\n   0 0 0 1 2 0.2 0.0\n   0 0 0 2 1 0.2 0.0\n",
              "pairlift": "PairLift\n2\n1\n 1\n   0 0 0 1 2 0.2 0.0\n   0 0 0 2 1 0.2 0.0\n"}

    def _dir(self):
        d = tempfile.mkdtemp(prefix="auto167_")
        self.addCleanup(shutil.rmtree, d, ignore_errors=True)
        return d

    def _build(self, scheme, interactions, hybridised, *, mode="RPA",
               subshape=(1, 1, 1), calc_type="ring", extern=None,
               coeff_extern=0.0, enable_spin_orbital=False, extra_param=None):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        import hwave.solver.flex as flex_mod
        d = self._dir()
        shutil.copy(os.path.join(_SRC, "geom.dat"), d)
        if hybridised:
            shutil.copy(os.path.join(_SRC, "transfer.dat"), d)
        else:
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write(_DIAGONAL_TRANSFER)
        idict = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat"}
        for key, fname in interactions.items():
            if fname in self.INLINE:
                with open(os.path.join(d, key + ".dat"), "w") as f:
                    f.write(self.INLINE[fname])
                idict[key] = key + ".dat"
            else:
                shutil.copy(os.path.join(_SRC, fname), d)
                idict[key] = fname
        if extern is not None:
            with open(os.path.join(d, "extern.dat"), "w") as f:
                f.write(self.INLINE[extern])
            idict["Extern"] = "extern.dat"
        reader = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
               "SubShape": list(subshape), "Nmat": 16,
               "coeff_extern": coeff_extern}
        par.update(extra_param or {})
        info_mode = {"mode": mode, "param": par,
                     "enable_spin_orbital": enable_spin_orbital,
                     "calc_scheme": scheme, "calc_type": calc_type}
        if mode == "FLEX":
            par.update({"IterationMax": 1, "Mix": 1.0})
            solver = flex_mod.FLEX(reader.get_param("ham"), {}, info_mode)
        else:
            solver = rpa_mod.RPA(reader.get_param("ham"), {}, info_mode)
        return solver, reader.get_param("green")

    def _solve(self, scheme, interactions, hybridised, **kw):
        solver, green_info = self._build(scheme, interactions, hybridised, **kw)
        out = self._dir()
        solver.solve(green_info, out)
        return solver, green_info, out


class TestRPAAutoResolution(_Case):
    """Promotion / non-promotion matrix, resolved through preview_scheme
    (pure) AND confirmed after solve (committed)."""

    def test_construction_leaves_auto_unresolved(self):
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        self.assertEqual(solver.calc_scheme, "auto")
        self.assertEqual(solver.calc_scheme_requested, "auto")
        self.assertIsNone(solver.enable_reduced)
        self.assertIsNone(solver._scheme_resolution)

    def test_diagonal_transfer_keeps_reduced(self):
        for key, fname in (("CoulombInter", "onsite_inter.dat"),
                           ("Hund", "hund_onsite.dat"),
                           ("Ising", "hund_onsite.dat"),
                           ("coulombinter", "onsite_inter.dat")):   # lowercase key
            with self.subTest(interaction=key):
                solver, green_info, _ = self._solve("auto", {key: fname}, False)
                self.assertEqual(solver.calc_scheme, "reduced")
                self.assertEqual(solver._scheme_resolution, "auto:exact:diagonal_transfer")
                self.assertEqual(np.asarray(green_info["chiq"]).ndim, 4)
                self.assertIs(solver.enable_reduced, True)

    def test_hybridised_transfer_promotes(self):
        for key, fname in (("CoulombInter", "onsite_inter.dat"),
                           ("Hund", "hund_onsite.dat"),
                           ("Ising", "hund_onsite.dat")):
            with self.subTest(interaction=key):
                solver, green_info, _ = self._solve("auto", {key: fname}, True)
                self.assertEqual(solver.calc_scheme, "general")
                self.assertEqual(solver._scheme_resolution, "auto:mixed:transfer")
                self.assertEqual(np.asarray(green_info["chiq"]).ndim, 6)

    def test_folded_diagonal_chain_keeps_reduced(self):
        solver, green_info, _ = self._solve("auto", {"CoulombInter": "onsite_inter.dat"},
                                            False, subshape=(1, 2, 1))
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertEqual(solver._scheme_resolution, "auto:exact:folded_diagonal")

    def test_coulombintra_only_never_promotes(self):
        solver, _ = self._build("auto", {"CoulombIntra": "coulombintra.dat"}, True)
        self.assertEqual(solver.preview_scheme(), ("reduced", "auto:no_discarded_content"))

    def test_extern_promotes_only_when_active(self):
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, False,
                                extern="offdiag_extern", coeff_extern=0.0)
        self.assertEqual(solver.preview_scheme(), ("reduced", "auto:exact:diagonal_transfer"))
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, False,
                                extern="offdiag_extern", coeff_extern=0.5)
        self.assertEqual(solver.preview_scheme(), ("general", "auto:mixed:extern"))
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, False,
                                extern="diag_extern", coeff_extern=0.5)
        self.assertEqual(solver.preview_scheme(), ("reduced", "auto:exact:diagonal_transfer"))

    def test_ring_ladder_and_general_only_precedence(self):
        solver, _ = self._build("auto", {"CoulombIntra": "coulombintra.dat"}, False,
                                calc_type="ring+ladder")
        self.assertEqual(solver.preview_scheme(), ("general", "auto:ring_ladder"))
        solver, _ = self._build("auto", {"Exchange": "exchange",
                                         "CoulombInter": "onsite_inter.dat"}, False)
        self.assertEqual(solver.preview_scheme(), ("general", "auto:general_only"))

    def test_preview_is_pure_and_agrees_with_the_commit(self):
        solver, green_info = self._build("auto", {"Hund": "hund_onsite.dat"}, True)
        before = (solver.calc_scheme, solver.enable_reduced, solver._scheme_resolution)
        preview = solver.preview_scheme()
        self.assertEqual((solver.calc_scheme, solver.enable_reduced,
                          solver._scheme_resolution), before)
        solver.solve(green_info, self._dir())
        self.assertEqual(preview, (solver.calc_scheme, solver._scheme_resolution))

    def test_preview_green_info_none_means_no_late_inputs(self):
        solver, _ = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        self.assertEqual(solver.preview_scheme(None)[1], "auto:exact:diagonal_transfer")
        self.assertEqual(solver.preview_scheme({"trans_mod": object()})[1],
                         "auto:mixed:trans_mod")
        self.assertEqual(solver.preview_scheme({"green_init": object()})[1],
                         "auto:mixed:green_init")

    def test_folded_extern_spin_block_rows_do_not_promote(self):
        """The Extern off-diagonal scan reads the PRE-FOLD table, so its
        index limit must be the PRE-FOLD physical orbital count.

        A 1-orbital chain folded onto a 2-cell supercell has post-fold
        norb == 2, so a post-fold limit would admit the spin-block rows
        (file index 2) that ``_make_ham_trans`` discards, and a
        flavour-conserving model would be promoted to 'general'.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        d = self._dir()
        for name, body in (("geom.dat", _CHAIN_GEOM),
                           ("transfer.dat", _CHAIN_TRANSFER),
                           ("coulombinter.dat", _CHAIN_COULOMBINTER),
                           ("extern.dat", _CHAIN_SPIN_BLOCK_EXTERN)):
            with open(os.path.join(d, name), "w") as f:
                f.write(body)
        idict = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat",
                 "Extern": "extern.dat"}
        reader = read_input_k.QLMSkInput({"path_to_input": d,
                                          "interaction": idict})
        info_mode = {"mode": "RPA",
                     "param": {"T": 2.0, "filling": 0.5,
                               "CellShape": [4, 1, 1], "SubShape": [2, 1, 1],
                               "Nmat": 16, "coeff_extern": 0.5},
                     "enable_spin_orbital": False,
                     "calc_scheme": "auto", "calc_type": "ring"}
        solver = rpa_mod.RPA(reader.get_param("ham"), {}, info_mode)
        # the fixture is only meaningful while the two candidate limits
        # differ: post-fold norb == 2, pre-fold physical norb == 1
        self.assertTrue(solver.lattice.has_sublattice)
        self.assertEqual((solver.ham_info.norb, solver.ham_info.norb_orig),
                         (2, 1))
        self.assertEqual(solver.preview_scheme(),
                         ("reduced", "auto:exact:folded_diagonal"))

    def test_show_params_prints_deferred_marker(self):
        solver, _ = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            solver._show_params()
        self.assertTrue(any("auto (resolved at read_init/solve)" in m
                            for m in cm.output), cm.output)

    def test_explicit_schemes_are_resolved_at_construction(self):
        for scheme in ("reduced", "general"):
            solver, _ = self._build(scheme, {"CoulombInter": "onsite_inter.dat"}, True)
            self.assertEqual(solver.calc_scheme, scheme)
            self.assertEqual(solver._scheme_resolution, "explicit")
            self.assertEqual(solver.preview_scheme(), (scheme, "explicit"))


class TestRPALegacyWarning(_Case):
    def test_warning_exactly_when_decision_differs_from_1_0(self):
        # promoted: WARNING naming the reason and the cost consequence
        solver, green_info = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            solver.solve(green_info, self._dir())
        warn = [r for r in cm.records if r.levelno == logging.WARNING
                and "1.0.x" in r.getMessage()]
        self.assertEqual(len(warn), 1, cm.output)
        self.assertIn("rank-6", warn[0].getMessage())
        self.assertIn("auto:mixed:transfer", warn[0].getMessage())
        # not promoted: INFO only
        solver, green_info = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, False)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            solver.solve(green_info, self._dir())
        self.assertFalse([r for r in cm.records if r.levelno == logging.WARNING
                          and "1.0.x" in r.getMessage()], cm.output)
        self.assertTrue(any("auto:exact:diagonal_transfer" in m for m in cm.output))


if __name__ == "__main__":
    unittest.main()
