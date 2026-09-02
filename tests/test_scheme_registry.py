"""Unit tests for hwave.solver.scheme -- the calc_scheme capability
registry and the auto-resolution rules (#167)."""
import math
import unittest

import numpy as np
from requests.structures import CaseInsensitiveDict


class TestCapabilityRegistryShape(unittest.TestCase):
    def test_every_record_is_named_and_closed(self):
        from hwave.solver import scheme
        for name, cap in scheme.CAPABILITIES.items():
            with self.subTest(name=name):
                self.assertIsInstance(cap, scheme.Capability)
                self.assertEqual(cap._fields,
                                 ("families", "rpa_mode", "flex_forcing",
                                  "sc_legacy_forcing"))
                self.assertIn(cap.rpa_mode, scheme.RPA_MODES)
                self.assertIsInstance(cap.flex_forcing, bool)
                self.assertIsInstance(cap.sc_legacy_forcing, bool)
                self.assertIsInstance(cap.families, frozenset)

    def test_registry_is_read_only_and_names_unique_case_insensitively(self):
        from hwave.solver import scheme
        with self.assertRaises(TypeError):
            scheme.CAPABILITIES["Hund"] = None  # MappingProxyType
        lowered = [k.lower() for k in scheme.CAPABILITIES]
        self.assertEqual(len(lowered), len(set(lowered)))

    def test_frozen_modes_match_the_spec_table(self):
        from hwave.solver import scheme
        expect = {
            "CoulombIntra": ("diag_ok", False, False),
            "CoulombInter": ("conditional", True, True),
            "Hund": ("conditional", True, True),
            "Ising": ("conditional", True, True),
            "Exchange": ("general_only", True, True),
            "PairHop": ("general_only", True, True),
            "PairLift": ("conditional", False, False),
            "Coulomb": ("conditional", True, True),
        }
        self.assertEqual(set(scheme.CAPABILITIES), set(expect))
        for name, (mode, flex, sc) in expect.items():
            cap = scheme.CAPABILITIES[name]
            self.assertEqual((cap.rpa_mode, cap.flex_forcing,
                              cap.sc_legacy_forcing), (mode, flex, sc), name)

    def test_reader_namelist_is_covered_exactly(self):
        """Adding a reader interaction without a capability record (or a
        record the reader cannot produce) fails here."""
        from hwave.solver import scheme
        from hwave.qlmsio.read_input_k import QLMSkInput
        reader_types = {s.lower() for s in QLMSkInput.valid_namelist}
        reader_types -= scheme.NON_INTERACTION_SECTIONS
        self.assertEqual(reader_types,
                         {k.lower() for k in scheme.CAPABILITIES})


class TestStructurallyNonzero(unittest.TestCase):
    def test_empty_and_all_zero_are_undeclared(self):
        from hwave.solver import scheme
        self.assertFalse(scheme.is_structurally_nonzero({}))
        self.assertFalse(scheme.is_structurally_nonzero(
            {((0, 0, 0), (0, 0)): 0.0, ((1, 0, 0), (0, 1)): 0j}))

    def test_numpy_scalars_count(self):
        from hwave.solver import scheme
        self.assertTrue(scheme.is_structurally_nonzero(
            {((0, 0, 0), (0, 0)): np.complex128(1e-300)}))
        self.assertFalse(scheme.is_structurally_nonzero(
            {((0, 0, 0), (0, 0)): np.complex128(0.0)}))

    def test_non_finite_raises_even_after_a_nonzero_entry(self):
        from hwave.solver import scheme
        for bad in (float("nan"), float("inf"), -float("inf"),
                    complex(0.0, float("nan"))):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "non-finite"):
                    scheme.is_structurally_nonzero(
                        {((0, 0, 0), (0, 0)): 1.0, ((1, 0, 0), (0, 0)): bad},
                        source="Hund")


class TestDeclaredTypes(unittest.TestCase):
    def _tables(self, **kw):
        t = CaseInsensitiveDict()
        t["Geometry"] = {"norb": 2}
        t["Transfer"] = {((0, 0, 0), (0, 0)): -1.0}
        for k, v in kw.items():
            t[k] = v
        return t

    def test_lowercase_keys_canonicalize(self):
        from hwave.solver import scheme
        t = self._tables(coulombinter={((0, 0, 0), (0, 1)): 0.5},
                         hund={((0, 0, 0), (0, 1)): 0.1})
        self.assertEqual(scheme.declared_types(t),
                         frozenset({"CoulombInter", "Hund"}))

    def test_sections_and_zero_tables_are_excluded(self):
        from hwave.solver import scheme
        t = self._tables(Extern={((0, 0, 0), (0, 1)): 0.3},
                         Exchange={((0, 0, 0), (0, 1)): 0.0},
                         CoulombIntra={})
        self.assertEqual(scheme.declared_types(t), frozenset())

    def test_unknown_type_fails_closed(self):
        from hwave.solver import scheme
        t = self._tables(Kondo={((0, 0, 0), (0, 0)): 1.0})
        with self.assertRaisesRegex(ValueError,
                                    "no capability entry for interaction 'Kondo'"):
            scheme.declared_types(t)

    def test_plain_dict_input_is_accepted(self):
        from hwave.solver import scheme
        self.assertEqual(
            scheme.declared_types({"PairLift": {((0, 0, 0), (0, 1)): 0.2}}),
            frozenset({"PairLift"}))


class TestLegacyRuleAndTokens(unittest.TestCase):
    def test_legacy_1_0_rule(self):
        from hwave.solver import scheme
        self.assertEqual(scheme.legacy_1_0_resolution(frozenset(), "ring+ladder"), "general")
        self.assertEqual(scheme.legacy_1_0_resolution({"Exchange"}, "ring"), "general")
        self.assertEqual(scheme.legacy_1_0_resolution({"PairHop"}, "ring"), "general")
        for t in ("CoulombIntra", "CoulombInter", "Hund", "Ising", "PairLift", "Coulomb"):
            self.assertEqual(scheme.legacy_1_0_resolution({t}, "ring"), "reduced", t)

    def test_token_vocabulary_is_closed(self):
        from hwave.solver import scheme
        self.assertEqual(scheme.RESOLUTION_TOKENS, frozenset({
            "explicit", "auto:ring_ladder", "auto:general_only",
            "auto:no_discarded_content", "auto:exact:diagonal_transfer",
            "auto:exact:folded_diagonal", "auto:mixed:transfer",
            "auto:mixed:extern", "auto:mixed:trans_mod",
            "auto:mixed:green_init", "auto:flex_forcing"}))


class TestFlavourConserved(unittest.TestCase):
    def _t(self, transfer, extern=None):
        t = CaseInsensitiveDict()
        t["Transfer"] = transfer
        if extern is not None:
            t["extern"] = extern
        return t

    def _call(self, tables, **kw):
        from hwave.solver import scheme
        args = dict(norb_phys=2, coeff_extern=0.0,
                    trans_mod_present=False, green_init_present=False)
        args.update(kw)
        return scheme.flavour_conserved(tables, **args)

    def test_diagonal_transfer_is_conserved(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0, ((1, 0, 0), (1, 1)): 0.5,
                     ((0, 0, 0), (0, 1)): 0.0})   # explicit zero does not promote
        self.assertEqual(self._call(t), (True, "diagonal_transfer"))

    def test_any_nonzero_offdiagonal_promotes_no_tolerance(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0, ((2, 0, 0), (0, 1)): 1e-300})
        self.assertEqual(self._call(t), (False, "transfer"))

    def test_non_finite_transfer_raises(self):
        t = self._t({((0, 0, 0), (0, 0)): float("nan")})
        with self.assertRaisesRegex(ValueError, "non-finite"):
            self._call(t)

    def test_extern_promotes_only_when_active(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0},
                    extern={((0, 0, 0), (0, 1)): 0.2})
        self.assertEqual(self._call(t, coeff_extern=0.0), (True, "diagonal_transfer"))
        self.assertEqual(self._call(t, coeff_extern=0.3), (False, "extern"))

    def test_extern_spin_block_entries_are_ignored(self):
        # indices >= norb_phys never enter H0 (rpa._make_ham_trans skips them)
        t = self._t({((0, 0, 0), (0, 0)): -1.0},
                    extern={((0, 0, 0), (2, 3)): 0.2, ((0, 0, 0), (0, 0)): 0.1})
        self.assertEqual(self._call(t, coeff_extern=1.0), (True, "diagonal_transfer"))

    def test_non_finite_coeff_extern_raises(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0}, extern={((0, 0, 0), (0, 0)): 0.1})
        with self.assertRaisesRegex(ValueError, "coeff_extern"):
            self._call(t, coeff_extern=float("inf"))

    def test_precedence_is_fixed(self):
        t = self._t({((0, 0, 0), (0, 1)): 0.3},        # transfer mixes
                    extern={((0, 0, 0), (0, 1)): 0.2})  # extern mixes
        self.assertEqual(self._call(t, coeff_extern=1.0), (False, "transfer"))
        self.assertEqual(self._call(t, coeff_extern=1.0, green_init_present=True),
                         (False, "green_init"))
        self.assertEqual(self._call(t, coeff_extern=1.0, green_init_present=True,
                                    trans_mod_present=True), (False, "trans_mod"))

    def test_spin_orbital_combined_index(self):
        # combined index 2*orb+spin: (0,1) is a spin flip on orbital 0.
        # The SO Transfer limit is 2*norb_phys (= geometry norb), so index 1
        # is IN range for norb_phys=1 and the flip promotes.
        t = self._t({((0, 0, 0), (0, 0)): -1.0, ((0, 0, 0), (0, 1)): 0.1})
        self.assertEqual(self._call(t, norb_phys=1, enable_spin_orbital=True),
                         (False, "transfer"))
        t = self._t({((0, 0, 0), (0, 0)): -1.0, ((0, 0, 0), (1, 1)): -1.0})
        self.assertEqual(self._call(t, norb_phys=1, enable_spin_orbital=True),
                         (True, "diagonal_transfer"))

    def test_transfer_spin_block_rows_are_ignored_in_normal_mode(self):
        # In normal mode _make_ham_trans keeps only indices < norb, silently
        # skipping the spin-extended rows a num_wann = 2*norb file carries.
        # The predicate must be blind to them too, or a flavour-conserving
        # model would be promoted on a row H0(k) never sees.
        t = self._t({((0, 0, 0), (0, 0)): -1.0, ((0, 0, 0), (2, 3)): 0.5})
        self.assertEqual(self._call(t, norb_phys=2), (True, "diagonal_transfer"))
        # the very same row DOES promote once it is in range
        self.assertEqual(self._call(t, norb_phys=4), (False, "transfer"))

    def test_out_of_range_transfer_non_finite_does_not_raise(self):
        # bounds filter first: an entry the consumer never reads cannot make
        # the predicate refuse to judge the input
        t = self._t({((0, 0, 0), (0, 0)): -1.0,
                     ((0, 0, 0), (2, 3)): float("nan")})
        self.assertEqual(self._call(t, norb_phys=2), (True, "diagonal_transfer"))

    def test_out_of_range_extern_non_finite_does_not_raise(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0},
                    extern={((0, 0, 0), (2, 3)): float("inf"),
                            ((0, 0, 0), (0, 1)): 0.0})
        self.assertEqual(self._call(t, coeff_extern=1.0), (True, "diagonal_transfer"))

    def test_in_range_transfer_non_finite_still_raises(self):
        t = self._t({((0, 0, 0), (0, 0)): -1.0,
                     ((0, 0, 0), (0, 1)): float("nan")})
        with self.assertRaisesRegex(ValueError, "non-finite"):
            self._call(t, norb_phys=2)


class TestResolvers(unittest.TestCase):
    def _rpa(self, types, calc_type="ring", conserved=True, cause="diagonal_transfer",
             has_sublattice=False):
        from hwave.solver import scheme
        return scheme.resolve_rpa(frozenset(types), calc_type, conserved=conserved,
                                  cause=cause, has_sublattice=has_sublattice)

    def test_rpa_decision_order(self):
        self.assertEqual(self._rpa({"CoulombIntra"}, calc_type="ring+ladder"),
                         ("general", "auto:ring_ladder"))
        self.assertEqual(self._rpa({"Exchange", "CoulombInter"}),
                         ("general", "auto:general_only"))
        self.assertEqual(self._rpa({"CoulombInter"}),
                         ("reduced", "auto:exact:diagonal_transfer"))
        self.assertEqual(self._rpa({"Hund"}, has_sublattice=True),
                         ("reduced", "auto:exact:folded_diagonal"))
        self.assertEqual(self._rpa({"Ising"}, conserved=False, cause="transfer"),
                         ("general", "auto:mixed:transfer"))
        self.assertEqual(self._rpa({"PairLift"}, conserved=False, cause="extern"),
                         ("general", "auto:mixed:extern"))
        self.assertEqual(self._rpa({"CoulombIntra"}, conserved=False, cause="transfer"),
                         ("reduced", "auto:no_discarded_content"))
        self.assertEqual(self._rpa(set()), ("reduced", "auto:no_discarded_content"))

    def test_every_rpa_token_is_in_the_vocabulary(self):
        from hwave.solver import scheme
        seen = set()
        for types in ({"Exchange"}, {"CoulombInter"}, {"CoulombIntra"}, set()):
            for ct in ("ring", "ring+ladder"):
                for cons, cause in ((True, "diagonal_transfer"), (False, "transfer"),
                                    (False, "extern"), (False, "trans_mod"),
                                    (False, "green_init")):
                    for sub in (False, True):
                        seen.add(self._rpa(types, ct, cons, cause, sub)[1])
        self.assertTrue(seen <= scheme.RESOLUTION_TOKENS, seen - scheme.RESOLUTION_TOKENS)

    def test_flex_rule(self):
        from hwave.solver import scheme
        self.assertEqual(scheme.resolve_flex({"CoulombInter"}), ("general", "auto:flex_forcing"))
        self.assertEqual(scheme.resolve_flex({"Coulomb"}), ("general", "auto:flex_forcing"))
        self.assertEqual(scheme.resolve_flex({"CoulombIntra"}), ("reduced", "auto:no_discarded_content"))
        self.assertEqual(scheme.resolve_flex({"PairLift"}), ("reduced", "auto:no_discarded_content"))
        self.assertEqual(scheme.resolve_flex(set()), ("reduced", "auto:no_discarded_content"))

    def test_estimate_chi_bytes(self):
        from hwave.solver import scheme
        self.assertEqual(scheme.estimate_chi_bytes("reduced", 16, 4, 3), 16 * 4 * 9 * 16)
        self.assertEqual(scheme.estimate_chi_bytes("general", 16, 4, 3), 16 * 4 * 81 * 16)
        # case-tolerant: an odd-cased explicit scheme must size, not crash
        # (this is a GPU-path pre-flight estimate; raising there would be a
        # new crash class for a request the solver itself accepts)
        self.assertEqual(scheme.estimate_chi_bytes("Reduced", 16, 4, 3),
                         16 * 4 * 9 * 16)
        self.assertEqual(scheme.estimate_chi_bytes("GENERAL", 16, 4, 3),
                         16 * 4 * 81 * 16)
        # ... but anything that is not a resolved scheme still raises
        for bad in ("auto", "AUTO", "bogus", None):
            with self.subTest(scheme=bad):
                with self.assertRaises(ValueError):
                    scheme.estimate_chi_bytes(bad, 16, 4, 3)


if __name__ == "__main__":
    unittest.main()
