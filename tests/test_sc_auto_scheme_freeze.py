# tests/test_sc_auto_scheme_freeze.py
"""Behaviour-freeze pin for hwave_sc's chi0q_tensor='auto' rule (#167,
contract 4: sc.py's mode is frozen while RPA/FLEX auto change)."""
import unittest


class TestScAutoSchemeFreeze(unittest.TestCase):
    #: (declared interaction keys) -> frozen 1.0.x decision
    CASES = {
        (): "reduced",
        ("path_to_input", "Geometry", "Transfer"): "reduced",
        ("CoulombIntra",): "reduced",
        ("PairLift",): "reduced",                      # frozen: NOT forcing
        ("CoulombInter",): "general",
        ("Hund",): "general",
        ("Ising",): "general",
        ("Exchange",): "general",
        ("PairHop",): "general",
        ("Coulomb",): "general",
        ("coulombinter",): "general",                  # case-insensitive
        ("hund", "CoulombIntra"): "general",
        ("Kondo",): "reduced",                          # unknown key ignored
        ("CoulombIntra", "PairLift"): "reduced",
    }

    def test_frozen_decisions(self):
        from hwave.sc import _auto_chi0q_tensor_scheme
        for keys, expect in self.CASES.items():
            with self.subTest(keys=keys):
                files = {k: "{}.dat".format(k) for k in keys}
                self.assertEqual(_auto_chi0q_tensor_scheme(files), expect)

    def test_every_registry_type_is_pinned(self):
        """A new reader type must get a row here (its sc decision is a
        deliberate choice, not a default)."""
        from hwave.solver import scheme
        pinned = {k for keys in self.CASES for k in keys}
        pinned = {k.lower() for k in pinned}
        for name in scheme.CAPABILITIES:
            self.assertIn(name.lower(), pinned, name)


if __name__ == "__main__":
    unittest.main()
