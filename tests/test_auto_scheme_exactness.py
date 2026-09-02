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

#: Transfer declared in the spin-block-extended layout (num_wann = 2 =
#: 2 * norb_phys) while the model runs in NORMAL (non spin-orbital) mode.
#: The normal branch of ``_make_ham_trans`` keeps only indices < norb (= 1
#: here) and skips the rest, so neither the spin-diagonal row (2,2) nor the
#: spin-block off-diagonal row (1,2) ever enters H0(k).
_CHAIN_SPIN_BLOCK_TRANSFER = """Transfer, spin-block extended (normal mode)
2
2
 1 1
   1    0    0    1    1  -1.0 0.0
  -1    0    0    1    1  -1.0 0.0
   0    0    0    2    2  -0.4 0.0
   1    0    0    1    2   0.7 0.0
  -1    0    0    2    1   0.7 0.0
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
               coeff_extern=0.0, enable_spin_orbital=False, extra_param=None,
               inject_tables=None):
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
        ham = reader.get_param("ham")
        # reader-bypass injection BEFORE construction: FLEX's explicit-reduced
        # diagnostic runs in the constructor, so a post-construction injection
        # could never reach it.
        for tname, tbl in (inject_tables or {}).items():
            ham[tname] = tbl
        if mode == "FLEX":
            par.update({"IterationMax": 1, "Mix": 1.0})
            solver = flex_mod.FLEX(ham, {}, info_mode)
        else:
            solver = rpa_mod.RPA(ham, {}, info_mode)
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

    def test_transfer_spin_block_rows_do_not_promote_in_normal_mode(self):
        """The Transfer off-diagonal scan must mirror ``_make_ham_trans``'s
        index bound.

        A 1-orbital model whose transfer.dat is written in the
        spin-extended layout (num_wann = 2) carries rows on indices >= norb
        that the normal branch of ``_make_ham_trans`` silently skips. They
        are not part of H0(k), so they must not promote the scheme.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        d = self._dir()
        for name, body in (("geom.dat", _CHAIN_GEOM),
                           ("transfer.dat", _CHAIN_SPIN_BLOCK_TRANSFER),
                           ("coulombinter.dat", _CHAIN_COULOMBINTER)):
            with open(os.path.join(d, name), "w") as f:
                f.write(body)
        idict = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat"}
        reader = read_input_k.QLMSkInput({"path_to_input": d,
                                          "interaction": idict})
        info_mode = {"mode": "RPA",
                     "param": {"T": 2.0, "filling": 0.5,
                               "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
                               "Nmat": 16, "coeff_extern": 0.0},
                     "enable_spin_orbital": False,
                     "calc_scheme": "auto", "calc_type": "ring"}
        solver = rpa_mod.RPA(reader.get_param("ham"), {}, info_mode)
        # the fixture is only meaningful while the file really declares rows
        # outside the consumed range
        self.assertEqual(solver.ham_info.norb_orig, 1)
        self.assertTrue(any(max(orbvec) >= 1 for _, orbvec
                            in solver.ham_info.param_ham["Transfer"].keys()))
        self.assertEqual(solver.preview_scheme(),
                         ("reduced", "auto:exact:diagonal_transfer"))

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


class TestAutoIsExact(_Case):
    """The contract, numerically: what auto returns equals the explicit
    exact scheme on the same input (bitwise where the regression suite
    already measures bit-identity; tolerance fallback per
    tests/equivalence_cells._candidate_atol is a recorded calibration event)."""

    CASES = (("CoulombInter", "onsite_inter.dat"),
             ("Hund", "hund_onsite.dat"),
             ("Ising", "hund_onsite.dat"))

    def test_retained_reduced_equals_the_general_density_projection(self):
        for key, fname in self.CASES:
            with self.subTest(interaction=key):
                _, gi_auto, _ = self._solve("auto", {key: fname}, False)
                _, gi_gen, _ = self._solve("general", {key: fname}, False)
                auto = np.asarray(gi_auto["chiq"])
                self.assertEqual(auto.ndim, 4)
                projected = _project_density_pairs(np.asarray(gi_gen["chiq"]))
                self.assertTrue(np.array_equal(projected, auto),
                                "max|diff| = {}".format(
                                    float(np.max(np.abs(projected - auto)))))

    def test_promoted_auto_equals_explicit_general_bitwise(self):
        for key, fname in self.CASES:
            with self.subTest(interaction=key):
                _, gi_auto, _ = self._solve("auto", {key: fname}, True)
                _, gi_gen, _ = self._solve("general", {key: fname}, True)
                self.assertTrue(np.array_equal(np.asarray(gi_auto["chiq"]),
                                               np.asarray(gi_gen["chiq"])))

    def test_folded_diagonal_reduced_equals_general_projection(self):
        _, gi_auto, _ = self._solve("auto", {"CoulombInter": "onsite_inter.dat"},
                                    False, subshape=(1, 2, 1))
        _, gi_gen, _ = self._solve("general", {"CoulombInter": "onsite_inter.dat"},
                                   False, subshape=(1, 2, 1))
        auto = np.asarray(gi_auto["chiq"])
        self.assertEqual(auto.ndim, 4)
        diff = float(np.max(np.abs(_project_density_pairs(np.asarray(gi_gen["chiq"])) - auto)))
        # folding changes summation order: round-off, not bit-identity
        # (measured 1.7e-18 in tests/test_reduced_scheme_exactness.py)
        self.assertLess(diff, 1e-12, diff)

    def test_pairlift_both_sides_of_the_conditional(self):
        # hybridised -> general (conditional + mixed); diagonal -> reduced AND
        # the reduced result equals the general density projection (records
        # the previously unmeasured type)
        # Measured outcome: the diagonal-branch agreement is TRIVIALLY exact
        # (diff 0.0) because PairLift's particle-hole vertex carries no
        # content that reaches the density-pair sector in the spin-free ring
        # (see PAIRLIFT_INERT_WARNING in rpa.py); the row records that fact
        # and the hybridised branch's routing.
        solver, gi, _ = self._solve("auto", {"PairLift": "pairlift"}, True)
        self.assertEqual((solver.calc_scheme, solver._scheme_resolution),
                         ("general", "auto:mixed:transfer"))
        _, gi_auto, _ = self._solve("auto", {"PairLift": "pairlift"}, False)
        _, gi_gen, _ = self._solve("general", {"PairLift": "pairlift"}, False)
        auto = np.asarray(gi_auto["chiq"])
        self.assertEqual(auto.ndim, 4)
        projected = _project_density_pairs(np.asarray(gi_gen["chiq"]))
        diff = float(np.max(np.abs(projected - auto)))
        self.assertTrue(np.array_equal(projected, auto),
                        "PairLift diagonal: reduced vs general projection max|diff| = {}"
                        .format(diff))


class TestReusePathAndFingerprint(_Case):
    """chi0q-reuse (file and in-memory), resolution/read ordering, and the
    fingerprint guard that catches a stale auto-resolved solver being fed a
    different input (#167)."""

    def _write_chi0q(self, scheme, hybridised):
        # RPA.save_results(info_outputfile, green_info) reads
        # info_outputfile["path_to_output"] (rpa.py:2387); np.savez appends .npz
        solver, green_info, out = self._solve(scheme, {"CoulombInter": "onsite_inter.dat"},
                                              hybridised)
        solver.save_results({"path_to_output": out, "chi0q": "chi0q"}, green_info)
        path = os.path.join(out, "chi0q.npz")
        self.assertTrue(os.path.isfile(path), os.listdir(out))
        return path

    def test_four_axis_chi0q_with_promoting_input_is_rejected_with_remediation(self):
        path = self._write_chi0q("reduced", True)
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        with self.assertRaisesRegex(ValueError, "auto:mixed:transfer.*rank-6.*calc_scheme = 'reduced'"):
            solver.read_init({"path_to_input": os.path.dirname(path),
                              "chi0q_init": os.path.basename(path)})

    def test_six_axis_chi0q_with_promoting_input_is_accepted(self):
        path = self._write_chi0q("general", True)
        solver, green_info = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        info = solver.read_init({"path_to_input": os.path.dirname(path),
                                 "chi0q_init": os.path.basename(path)})
        self.assertEqual(solver._scheme_resolution, "auto:mixed:transfer")
        green_info.update(info)
        solver.solve(green_info, self._dir())
        self.assertEqual(np.asarray(green_info["chiq"]).ndim, 6)

    def test_resolution_happens_before_the_chi0q_load(self):
        path = self._write_chi0q("general", True)
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        order = []
        real_resolve, real_read = solver._resolve_auto_scheme, solver._read_chi0q

        def spy_resolve(**kw):
            order.append("resolve"); return real_resolve(**kw)

        def spy_read(fn):
            order.append("read"); return real_read(fn)
        solver._resolve_auto_scheme, solver._read_chi0q = spy_resolve, spy_read
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            solver.read_init({"path_to_input": os.path.dirname(path),
                              "chi0q_init": os.path.basename(path)})
        self.assertEqual(order[:2], ["resolve", "read"])
        msgs = cm.output
        i_est = next(i for i, m in enumerate(msgs) if "chi0q_init load" in m)
        i_read = next(i for i, m in enumerate(msgs) if "read_chi0q" in m)
        self.assertLess(i_est, i_read, "estimate must be logged before the load")

    def test_both_entry_orders_agree(self):
        # read_init first, then solve
        solver_a, gi_a = self._build("auto", {"Hund": "hund_onsite.dat"}, True)
        gi_a.update(solver_a.read_init({}))
        solver_a.solve(gi_a, self._dir())
        # solve directly
        solver_b, gi_b = self._build("auto", {"Hund": "hund_onsite.dat"}, True)
        solver_b.solve(gi_b, self._dir())
        self.assertEqual((solver_a.calc_scheme, solver_a._scheme_resolution),
                         (solver_b.calc_scheme, solver_b._scheme_resolution))
        self.assertTrue(np.array_equal(np.asarray(gi_a["chiq"]), np.asarray(gi_b["chiq"])))

    def test_stale_fingerprint_raises(self):
        solver, green_info = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        solver.read_init({})                       # resolves: no late inputs
        self.assertEqual(solver._scheme_resolution, "auto:exact:diagonal_transfer")
        H = np.asarray(solver.ham_info.ham_trans_q)          # (nvol, norb, norb)
        nvol, norb = H.shape[0], H.shape[1]
        trans_mod = np.einsum("kab,st->ksatb", H, np.eye(2)).reshape(nvol, 2 * norb, 2 * norb)
        green_info["trans_mod"] = trans_mod
        with self.assertRaisesRegex(ValueError, "resolved to 'reduced'"):
            solver.solve(green_info, self._dir())

    def test_mutated_decision_input_is_detected_on_re_resolution(self):
        """The fingerprint records declared types and late-input presence,
        but the decision also reads construction-immutable inputs
        (coeff_extern, the transfer/extern structure). Mutating one after
        resolution violates the #167 immutability invariant; the re-resolve
        must not silently keep the stale scheme."""
        solver, green_info = self._build(
            "auto", {"CoulombInter": "onsite_inter.dat"}, False,
            extern="offdiag_extern", coeff_extern=0.0)
        solver.read_init({})
        self.assertEqual((solver.calc_scheme, solver._scheme_resolution),
                         ("reduced", "auto:exact:diagonal_transfer"))
        solver.ext = 0.5          # the invariant violation
        with self.assertRaises(ValueError) as cm:
            solver.solve(green_info, self._dir())
        msg = str(cm.exception)
        self.assertIn("reduced", msg)
        self.assertIn("auto:exact:diagonal_transfer", msg)
        self.assertIn("general", msg)
        self.assertIn("auto:mixed:extern", msg)
        self.assertIn("coeff_extern", msg)

    def test_trans_mod_promotes_via_public_solve(self):
        solver, green_info = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        H = np.asarray(solver.ham_info.ham_trans_q)
        nvol, norb = H.shape[0], H.shape[1]
        green_info["trans_mod"] = np.einsum("kab,st->ksatb", H, np.eye(2)).reshape(
            nvol, 2 * norb, 2 * norb)
        solver.solve(green_info, self._dir())
        self.assertEqual((solver.calc_scheme, solver._scheme_resolution),
                         ("general", "auto:mixed:trans_mod"))

    def test_spin_orbital_flip_promotes_and_diagonal_keeps_reduced(self):
        # SO geometry: norb (spin-orbital) = 4 for the 2-orbital model; the
        # combined-index-diagonal transfer is the physical one duplicated
        # per spin; a spin flip is an off-diagonal (1,2) entry
        d = self._dir()
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n4\n0 0 0\n0 0 0\n0 0 0\n0 0 0\n")
        diag = ("Transfer SO diag\n4\n2\n 1 1\n" +
                "".join("   0 {r} 0 {i} {i} 1.0 0.0\n".format(r=r, i=i)
                        for r in (1, -1) for i in (1, 2, 3, 4)))
        flip = diag + "   0 0 0 1 2 0.1 0.0\n   0 0 0 2 1 0.1 0.0\n"
        inter = "CoulombInter\n2\n1\n 1\n   0 0 0 1 2 0.5 0.0\n   0 0 0 2 1 0.5 0.0\n"
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod
        for body, expect in ((diag, "auto:exact:diagonal_transfer"),
                             (flip, "auto:mixed:transfer")):
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write(body)
            with open(os.path.join(d, "inter.dat"), "w") as f:
                f.write(inter)
            reader = read_input_k.QLMSkInput({"path_to_input": d, "interaction": {
                "path_to_input": d, "Geometry": "geom.dat",
                "Transfer": "transfer.dat", "CoulombInter": "inter.dat"}})
            solver = rpa_mod.RPA(reader.get_param("ham"), {}, {
                "mode": "RPA", "enable_spin_orbital": True, "calc_scheme": "auto",
                "calc_type": "ring",
                "param": {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
                          # explicit: rpa.py's own default is SubShape ==
                          # CellShape (full folding to a single supercell,
                          # not "no folding"); every other fixture in this
                          # file goes through _Case._build, which sets
                          # SubShape=[1,1,1] by default (test-premise fix,
                          # #167 task 6).
                          "SubShape": [1, 1, 1], "Nmat": 16}})
            self.assertEqual(solver.preview_scheme()[1], expect, body)


class TestAtomicity(_Case):
    """A raise in any phase of _resolve_auto_scheme leaves the solver in the
    AUTO-UNRESOLVED state (never a half-committed one), and a subsequent
    retry with valid inputs still succeeds (#167)."""

    def _assert_unresolved(self, solver):
        self.assertEqual(solver.calc_scheme, "auto")
        self.assertIsNone(solver.enable_reduced)
        self.assertIsNone(solver._scheme_resolution)
        self.assertIsNone(solver._scheme_fingerprint)

    def test_each_failure_phase_leaves_auto_unresolved_and_retry_succeeds(self):
        from hwave.solver import scheme as sch
        from unittest import mock
        solver, green_info = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        phases = {
            "discovery": ("declared_types", ValueError("injected discovery")),
            "predicate": ("flavour_conserved", ValueError("injected predicate")),
            "validation": ("resolve_rpa", lambda *a, **k: ("bogus", "auto:bogus")),
        }
        for name, (target, effect) in phases.items():
            with self.subTest(phase=name):
                kw = {"side_effect": effect} if isinstance(effect, Exception) else {"side_effect": effect}
                with mock.patch.object(sch, target, **kw):
                    with self.assertRaises((ValueError, AssertionError)):
                        solver._resolve_auto_scheme(trans_mod_present=False,
                                                    green_init_present=False)
                self._assert_unresolved(solver)
        # fingerprint phase: _decide_auto_scheme returns a well-formed
        # (chosen, token) pair -- so the deferred-validation assertions on
        # chosen/token would pass -- but a non-iterable `types`. The
        # fingerprint line (fingerprint = (frozenset(types), ...)) sits
        # BEFORE those validation assertions in _resolve_auto_scheme, so
        # frozenset(42) raises TypeError exactly there, never reaching
        # validation.
        with mock.patch.object(solver, "_decide_auto_scheme",
                               return_value=("reduced",
                                             "auto:exact:diagonal_transfer",
                                             42, True)):
            with self.assertRaises(TypeError):
                solver._resolve_auto_scheme(trans_mod_present=False, green_init_present=False)
        self._assert_unresolved(solver)
        # retry
        solver._resolve_auto_scheme(trans_mod_present=False, green_init_present=False)
        self.assertEqual(solver._scheme_resolution, "auto:exact:diagonal_transfer")

    def test_emission_failure_never_propagates(self):
        from unittest import mock
        solver, _ = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        with mock.patch.object(solver, "_emit_auto_resolution_messages",
                               side_effect=RuntimeError("boom")):
            solver._resolve_auto_scheme(trans_mod_present=False, green_init_present=False)
        self.assertEqual(solver._scheme_resolution, "auto:exact:diagonal_transfer")

    def test_unknown_declared_type_fails_closed_via_public_path(self):
        solver, green_info = self._build("auto", {"Hund": "hund_onsite.dat"}, False)
        # inject a foreign table into the reader container (reader-bypass)
        solver.ham_info.param_ham["Kondo"] = {((0, 0, 0), (0, 0)): 1.0}
        with self.assertRaisesRegex(ValueError, "no capability entry for interaction 'Kondo'"):
            solver.solve(green_info, self._dir())
        self._assert_unresolved(solver)

    def test_all_zero_table_is_undeclared_and_non_finite_raises(self):
        solver, _ = self._build("auto", {"Hund": "hund_onsite.dat"}, True)
        solver.ham_info.param_ham["Exchange"] = {((0, 0, 0), (0, 1)): 0.0}
        self.assertEqual(solver.preview_scheme()[1], "auto:mixed:transfer")  # not general_only
        for bad in (float("nan"), float("inf")):
            solver.ham_info.param_ham["Exchange"] = {((0, 0, 0), (0, 1)): bad}
            with self.assertRaisesRegex(ValueError, "non-finite"):
                solver.preview_scheme()


class TestExplicitReducedDiagnostic(_Case):
    def test_one_warning_when_reduced_is_an_approximation(self):
        solver, green_info = self._build("reduced", {"CoulombInter": "onsite_inter.dat"}, True)
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            green_info.update(solver.read_init({}))
            solver.solve(green_info, self._dir())
        hits = [m for m in cm.output if "INPUT-DEPENDENT" in m]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertIn("CoulombInter", hits[0])
        self.assertIn("2.3e-4 / 3.3e-4 / 3.2e-4", hits[0])

    def test_no_warning_when_reduced_is_exact(self):
        solver, green_info = self._build("reduced", {"CoulombInter": "onsite_inter.dat"}, False)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            solver.solve(green_info, self._dir())
        self.assertFalse([m for m in cm.output if "INPUT-DEPENDENT" in m], cm.output)

    def test_dedup_flag_stays_disarmed_across_a_conserved_no_op_call(self):
        # read_init on a diagonal (flavour-conserving) transfer: conditional
        # type (Hund) declared, but conserved -- no warning, and the dedup
        # flag must NOT be armed by this no-op (#167 review finding).
        solver, green_info = self._build("reduced", {"Hund": "hund_onsite.dat"}, False)
        with self.assertLogs("hwave.solver.rpa", level="INFO") as cm:
            green_info.update(solver.read_init({}))
        self.assertFalse([m for m in cm.output if "INPUT-DEPENDENT" in m], cm.output)
        # now diverge the presence flags at the OTHER call point: a trans_mod
        # that mixes flavour (built exactly as
        # TestReusePathAndFingerprint.test_trans_mod_promotes_via_public_solve)
        H = np.asarray(solver.ham_info.ham_trans_q)
        nvol, norb = H.shape[0], H.shape[1]
        green_info["trans_mod"] = np.einsum("kab,st->ksatb", H, np.eye(2)).reshape(
            nvol, 2 * norb, 2 * norb)
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            solver.solve(green_info, self._dir())
        hits = [m for m in cm.output if "INPUT-DEPENDENT" in m]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertIn("trans_mod", hits[0])
        # dedup: once armed (by the actual emission above), a further
        # warning-triggering evaluation emits nothing at all.
        with self.assertRaises(AssertionError):
            with self.assertLogs("hwave.solver.rpa", level="WARNING"):
                solver._emit_reduced_exactness_diagnostic(
                    trans_mod_present=True, green_init_present=False)

    # ---------------------------------------------------------------- #167
    # The explicit-scheme diagnostics are ADVISORY: an explicit reduced run
    # must compute exactly what 1.0.x computed, so a diagnostic that cannot
    # judge the input is skipped (debug), never raised. Fail-closed discovery
    # is kept for the AUTO path, which is where the decision is made.
    def test_unknown_table_never_breaks_an_explicit_reduced_solve(self):
        solver, green_info = self._build("reduced", {"Hund": "hund_onsite.dat"}, True)
        solver.ham_info.param_ham["Kondo"] = {((0, 0, 0), (0, 0)): 1.0}
        with self.assertLogs("hwave.solver.rpa", level="DEBUG") as cm:
            solver.solve(green_info, self._dir())
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertEqual(np.asarray(green_info["chiq"]).ndim, 4)
        self.assertTrue(any("reduced-exactness diagnostic skipped" in m
                            for m in cm.output), cm.output)

    def test_non_finite_table_never_breaks_an_explicit_reduced_solve(self):
        # The reader rejects non-finite input files since #130, so injecting
        # into the already-built table is the honest reproduction of a
        # non-finite entry reaching the diagnostic. The interaction arrays were
        # built at construction, so the SOLVE numerics are unaffected here; the
        # point of the test is only that no ValueError escapes.
        solver, green_info = self._build("reduced", {"Hund": "hund_onsite.dat"}, True)
        key = sorted(solver.ham_info.param_ham["Hund"].keys())[0]
        solver.ham_info.param_ham["Hund"][key] = float("nan")
        with self.assertLogs("hwave.solver.rpa", level="DEBUG") as cm:
            solver.solve(green_info, self._dir())
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertIn("chiq", green_info)
        self.assertTrue(any("reduced-exactness diagnostic skipped" in m
                            for m in cm.output), cm.output)


class TestFLEXAutoResolution(_Case):
    """FLEX's auto rule is H0-INDEPENDENT and decided in the constructor
    (#167): general iff any declared type carries FLEX-forcing vertex
    content, because the dressed Green function can hybridise during the
    SCF iteration no matter how diagonal H0(k) starts out."""

    def test_onsite_coulombinter_diagonal_h0_promotes_h0_independent(self):
        # the FLEX rule is H0-independent: diagonal transfer, still general
        solver, _ = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, False, mode="FLEX")
        self.assertEqual((solver.calc_scheme, solver._scheme_resolution),
                         ("general", "auto:flex_forcing"))
        self.assertIs(solver.enable_reduced, False)
        self.assertTrue(solver._flex_general)

    def test_coulombintra_only_keeps_reduced(self):
        solver, _ = self._build("auto", {"CoulombIntra": "coulombintra.dat"}, True, mode="FLEX")
        self.assertEqual((solver.calc_scheme, solver._scheme_resolution),
                         ("reduced", "auto:no_discarded_content"))

    def test_pairlift_only_keeps_reduced_with_inertness_warning(self):
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            solver, _ = self._build("auto", {"PairLift": "pairlift"}, True, mode="FLEX")
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertTrue(any("PairLift" in m and "exactly zero" in m for m in cm.output),
                        cm.output)

    def test_ring_ladder_rejected_before_resolution(self):
        with self.assertRaisesRegex(ValueError, r"ring\+ladder"):
            self._build("auto", {"CoulombIntra": "coulombintra.dat"}, False,
                        mode="FLEX", calc_type="ring+ladder")

    def test_explicit_reduced_with_forcing_type_warns_at_construction(self):
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            self._build("reduced", {"CoulombInter": "onsite_inter.dat"}, False, mode="FLEX")
        hits = [m for m in cm.output if "regardless of the transfer structure" in m]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertIn("CoulombInter", hits[0])

    def test_explicit_reduced_never_emits_the_rpa_input_dependent_diagnostic(self):
        """FLEX warns ONCE, at construction, unconditionally of H0; RPA's
        H0-conditional diagnostic (whose text cites 'RPA instability' and
        RPA fixture figures) must not also fire on a FLEX solver. Both
        module loggers are captured via their common parent."""
        with self.assertLogs("hwave.solver", level="WARNING") as cm:
            solver, green_info = self._build(
                "reduced", {"CoulombInter": "onsite_inter.dat"}, True, mode="FLEX")
            green_info.update(solver.read_init({}))
            solver.solve(green_info, self._dir())
        hits = [m for m in cm.output if "regardless of the transfer structure" in m]
        self.assertEqual(len(hits), 1, cm.output)
        self.assertFalse([m for m in cm.output if "INPUT-DEPENDENT" in m], cm.output)

    def test_promoted_auto_equals_explicit_general_one_shot(self):
        _, gi_auto, _ = self._solve("auto", {"CoulombInter": "onsite_inter.dat"}, True, mode="FLEX")
        _, gi_gen, _ = self._solve("general", {"CoulombInter": "onsite_inter.dat"}, True, mode="FLEX")
        for key in ("chiq_s", "chiq_c"):
            self.assertTrue(np.array_equal(np.asarray(gi_auto[key]), np.asarray(gi_gen[key])), key)

    def test_unsupported_scheme_name_still_raises_the_actionable_valueerror(self):
        """The step-0 restructuring must not turn an unsupported scheme
        name into an AssertionError. 'AUTO' is the sharp case: the
        inherited _set_scheme compares case-SENSITIVELY, so it never routes
        a mis-cased request through the auto path, and FLEX must report it
        the way it always did."""
        for name in ("bogus", "AUTO"):
            with self.subTest(calc_scheme=name):
                with self.assertRaises(ValueError) as cm:
                    self._build(name, {"CoulombIntra": "coulombintra.dat"},
                                False, mode="FLEX")
                self.assertIn("FLEX requires calc_scheme='reduced' or "
                              "'general'", str(cm.exception))

    def test_read_init_and_solve_are_no_ops_for_the_resolved_state(self):
        solver, green_info = self._build("auto", {"CoulombIntra": "coulombintra.dat"}, True, mode="FLEX")
        green_info.update(solver.read_init({}))
        solver.solve(green_info, self._dir())
        self.assertEqual(solver._scheme_resolution, "auto:no_discarded_content")

    def test_unknown_table_never_breaks_an_explicit_reduced_construction(self):
        """FLEX's explicit-reduced diagnostic runs at construction; an input
        it cannot judge must be skipped, not raised (the computation itself
        stays the 1.0.x one)."""
        with self.assertLogs("hwave.solver.flex", level="DEBUG") as cm:
            solver, _ = self._build(
                "reduced", {"CoulombInter": "onsite_inter.dat"}, False,
                mode="FLEX",
                inject_tables={"Kondo": {((0, 0, 0), (0, 0)): 1.0}})
        self.assertEqual(solver.calc_scheme, "reduced")
        self.assertEqual(solver._scheme_resolution, "explicit")
        self.assertTrue(any("reduced-exactness diagnostic skipped" in m
                            for m in cm.output), cm.output)


class TestSchemeStamp(_Case):
    STAMP = ("calc_scheme", "calc_scheme_requested", "scheme_resolution")

    def _saved(self, mode, scheme, hybridised, outputs):
        solver, green_info, out = self._solve(scheme, {"CoulombInter": "onsite_inter.dat"},
                                              hybridised, mode=mode)
        info_out = {"path_to_output": out}
        info_out.update(outputs)
        solver.save_results(info_out, green_info)
        return out

    def test_rpa_writers_stamp_plain_unicode(self):
        from hwave.solver import scheme as sch
        out = self._saved("RPA", "auto", True, {"chiq": "chiq", "chi0q": "chi0q"})
        for name, expect in (("chiq.npz", ("general", "auto", "auto:mixed:transfer")),
                             ("chi0q.npz", ("general", "auto", "auto:mixed:transfer"))):
            with np.load(os.path.join(out, name), allow_pickle=False) as z:
                for key, val in zip(self.STAMP, expect):
                    self.assertIn(key, z.files, name)
                    self.assertEqual(z[key].dtype.kind, "U", (name, key))
                    self.assertEqual(z[key].item(), val, (name, key))
                self.assertIn(z["scheme_resolution"].item(), sch.RESOLUTION_TOKENS)

    def test_explicit_runs_stamp_explicit(self):
        out = self._saved("RPA", "reduced", False, {"chiq": "chiq"})
        with np.load(os.path.join(out, "chiq.npz"), allow_pickle=False) as z:
            self.assertEqual((z["calc_scheme"].item(), z["calc_scheme_requested"].item(),
                              z["scheme_resolution"].item()), ("reduced", "reduced", "explicit"))

    def test_flex_writers_stamp(self):
        out = self._saved("FLEX", "auto", False, {"chi0q": "chi0q", "chiq_s": "chiq_s", "chiq_c": "chiq_c"})
        for name in ("chi0q.npz", "chiq_s.npz", "chiq_c.npz"):
            with np.load(os.path.join(out, name), allow_pickle=False) as z:
                self.assertEqual(z["scheme_resolution"].item(), "auto:flex_forcing", name)
                self.assertEqual(z["calc_scheme"].item(), "general", name)

    def test_stamped_chi0q_round_trips_through_the_rpa_reader(self):
        out = self._saved("RPA", "general", True, {"chi0q": "chi0q"})
        solver, green_info = self._build("auto", {"CoulombInter": "onsite_inter.dat"}, True)
        green_info.update(solver.read_init({"path_to_input": out, "chi0q_init": "chi0q.npz"}))
        solver.solve(green_info, self._dir())
        self.assertEqual(np.asarray(green_info["chiq"]).ndim, 6)

    def test_stamp_less_pre_2_0_file_is_still_accepted(self):
        out = self._saved("RPA", "general", True, {"chi0q": "chi0q"})
        path = os.path.join(out, "chi0q.npz")
        with np.load(path, allow_pickle=False) as z:
            legacy = {k: z[k] for k in z.files if k not in self.STAMP}
        np.savez(path, **legacy)
        solver, green_info = self._build("general", {"CoulombInter": "onsite_inter.dat"}, True)
        green_info.update(solver.read_init({"path_to_input": out, "chi0q_init": "chi0q.npz"}))
        solver.solve(green_info, self._dir())


if __name__ == "__main__":
    unittest.main()
