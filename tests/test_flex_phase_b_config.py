"""Phase B configuration surface of the FLEX solver (spec 2026-09-06
sections 4.1 and D6/D9): the keys, their defaults, the refusal precedence
steps 1-3 (raw key types, domain keys, gate-HF coupling), the stale-key
warning, and the mode/key matrix (FLEX-only keys are warned-and-ignored
outside FLEX)."""
import logging
import os
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN = "tests/rpa/input_2orb"


def _build(param_extra=None, mode="FLEX", calc_scheme="general", interactions=None,
           info_extra=None):
    idict = {"path_to_input": _IN, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": _IN, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
    par.update(param_extra or {})
    info = {"mode": mode, "param": par, "enable_spin_orbital": False, "calc_scheme": calc_scheme}
    info.update(info_extra or {})
    if mode == "FLEX":
        import hwave.solver.flex as flex_mod
        return flex_mod.FLEX(r.get_param("ham"), {}, info), r
    import hwave.solver.rpa as rpa_mod
    info["calc_type"] = "ring"
    return rpa_mod.RPA(r.get_param("ham"), {}, info), r


class TestPhaseBConfig(unittest.TestCase):

    def test_defaults(self):
        s, _ = _build()
        self.assertFalse(s.flex_hartree_fock)
        self.assertFalse(s.longitudinal_bond_channels)
        self.assertFalse(s.longitudinal_bond_output_full)
        self.assertIsNone(s.longitudinal_bond_freq_batch)
        self.assertIsNone(s.longitudinal_bond_max_shells)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb, 8.0)
        self.assertFalse(s._phase_b_active)

    def test_keys_case_insensitive_and_active(self):
        s, _ = _build({"Flex_Hartree_Fock": True, "LONGITUDINAL_bond_channels": True,
                       "longitudinal_bond_output_full": True, "longitudinal_bond_freq_batch": 4,
                       "longitudinal_bond_max_shells": 2, "longitudinal_bond_memory_cap_gb": 3.5})
        self.assertTrue(s.flex_hartree_fock and s.longitudinal_bond_channels)
        self.assertTrue(s.longitudinal_bond_output_full)
        self.assertEqual(s.longitudinal_bond_freq_batch, 4)
        self.assertEqual(s.longitudinal_bond_max_shells, 2)
        self.assertEqual(s.longitudinal_bond_memory_cap_gb, 3.5)
        self.assertTrue(s._phase_b_active)
        s, _ = _build({"flex_hartree_fock": True})
        self.assertTrue(s._phase_b_active)

    def test_refusals(self):
        hf = {"flex_hartree_fock": True}
        gate = {"flex_hartree_fock": True, "longitudinal_bond_channels": True}
        cases = [
            ({"flex_hartree_fock": "yes"}, {}, "boolean"),
            ({"longitudinal_bond_channels": 1}, {}, "boolean"),
            (dict(hf, Nmat=31), {}, "even"),
            (dict(hf, matsubara_basis="ir"), {}, "matsubara_basis"),
            (dict(hf, gpu=True), {}, "gpu"),
            (dict(hf, SubShape=[2, 1, 1]), {}, "sublattice"),
            ({"longitudinal_bond_channels": True}, {}, "flex_hartree_fock"),
            (dict(gate, longitudinal_bond_freq_batch=0), {}, "longitudinal_bond_freq_batch"),
            (dict(gate, longitudinal_bond_freq_batch=33), {}, "longitudinal_bond_freq_batch"),
            (dict(gate, longitudinal_bond_max_shells=0), {}, "longitudinal_bond_max_shells"),
            (dict(gate, longitudinal_bond_memory_cap_gb=0.0), {}, "longitudinal_bond_memory_cap_gb"),
        ]
        for extra, info_extra, fragment in cases:
            with self.subTest(extra=extra):
                with self.assertRaises(ValueError) as cm:
                    _build(extra, info_extra=info_extra)
                self.assertIn(fragment, str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            _build(hf, calc_scheme="reduced")
        self.assertIn("general", str(cm.exception))

    def test_refusals_precede_the_base_constructor(self):
        """Odd Nmat is a ValueError from the FLEX parser (before the base
        constructor's legacy exit) whenever HF or the gate is on; a legacy
        input keeps the legacy exit."""
        with self.assertRaises(ValueError):
            _build({"flex_hartree_fock": True, "Nmat": 31})
        with self.assertRaises(SystemExit):
            _build({"Nmat": 31})

    def test_stale_bond_keys_warn_once_and_are_not_parsed(self):
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            s, _ = _build({"longitudinal_bond_output_full": "junk", "longitudinal_bond_freq_batch": -7,
                           "longitudinal_bond_max_shells": "x", "longitudinal_bond_memory_cap_gb": "y"})
        msgs = [m for m in cm.output if "longitudinal_bond_output_full" in m]
        self.assertEqual(len(msgs), 1)
        for k in ("longitudinal_bond_freq_batch", "longitudinal_bond_max_shells", "longitudinal_bond_memory_cap_gb"):
            self.assertIn(k, msgs[0])
        self.assertFalse(s.longitudinal_bond_output_full)
        self.assertIsNone(s.longitudinal_bond_freq_batch)
        with self.assertLogs("hwave.solver.flex", level="WARNING") as cm:
            _build({"flex_hartree_fock": True, "longitudinal_bond_freq_batch": 3})
        self.assertTrue(any("longitudinal_bond_freq_batch" in m for m in cm.output))

    def test_flex_only_keys_are_ignored_with_a_warning_in_rpa(self):
        with self.assertLogs("hwave.solver.rpa", level="WARNING") as cm:
            s, _ = _build({"flex_hartree_fock": True, "longitudinal_bond_output_full": True,
                           "longitudinal_bond_freq_batch": 2}, mode="RPA")
        msgs = [m for m in cm.output if "FLEX-only" in m]
        self.assertEqual(len(msgs), 1)
        for k in ("flex_hartree_fock", "longitudinal_bond_output_full", "longitudinal_bond_freq_batch"):
            self.assertIn(k, msgs[0])
        self.assertFalse(hasattr(s, "flex_hartree_fock"))

    def test_no_new_warning_for_legacy_inputs(self):
        logger = logging.getLogger("hwave.solver.flex")
        records = []
        h = logging.Handler(); h.emit = lambda rec: records.append(rec.getMessage())
        logger.addHandler(h)
        try:
            _build({})
        finally:
            logger.removeHandler(h)
        self.assertFalse([m for m in records if "longitudinal_bond" in m or "flex_hartree_fock" in m])


if __name__ == "__main__":
    unittest.main()
