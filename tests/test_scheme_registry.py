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


if __name__ == "__main__":
    unittest.main()
