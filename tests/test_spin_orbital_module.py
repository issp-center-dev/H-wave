"""The shared spin-orbital conventions (hwave.solver.spin_orbital, #83):
index reorders, the physical-orbital count, the strict switch, and the
public entries that now read them (RPA, UHFk)."""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver import spin_orbital as so


class TestIndexOrders(unittest.TestCase):
    def test_spinblock_order_matches_the_interleaved_slices(self):
        for n in (1, 2, 3, 5):
            rng = np.random.default_rng(n)
            m = rng.normal(size=(4, 2 * n, 2 * n))
            sb = so.interleaved_to_spinblock(m, n).reshape(4, 2, n, 2, n)
            for s in range(2):
                for t in range(2):
                    np.testing.assert_array_equal(sb[:, s, :, t, :], m[:, s::2, t::2])

    def test_round_trip_and_identity_for_one_orbital(self):
        m = np.arange(36.0).reshape(6, 6)
        np.testing.assert_array_equal(
            so.spinblock_to_interleaved(so.interleaved_to_spinblock(m, 3), 3), m)
        m2 = np.arange(4.0).reshape(2, 2)
        np.testing.assert_array_equal(so.interleaved_to_spinblock(m2, 1), m2)

    def test_scatter_and_gather_forms_agree(self):
        for n in (1, 2, 4):
            g = so.spinblock_order(n)
            for j, i in enumerate(g):
                self.assertEqual(so.to_spinblock_index(i, n), j)
                self.assertEqual(so.spin_of(i), j // n)
                self.assertEqual(so.interleaved_index(j % n, j // n), i)

    def test_only_the_requested_axes_move(self):
        m = np.arange(2 * 4 * 4).reshape(2, 4, 4)
        out = so.interleaved_to_spinblock(m, 2, axes=(-1,))
        np.testing.assert_array_equal(out, m[:, :, so.spinblock_order(2)])


class TestPhysicalNorb(unittest.TestCase):
    def test_halves_in_spin_orbital_mode_only(self):
        self.assertEqual(so.physical_norb(6, True), 3)
        self.assertEqual(so.physical_norb(6, False), 6)

    def test_odd_count_raises_naming_the_source(self):
        with self.assertRaisesRegex(ValueError, "even Geometry norb.*3 in geom.dat"):
            so.physical_norb(3, True)

    def test_parity_is_checked_on_the_original_count(self):
        self.assertEqual(so.physical_norb(8, True, check_norb=4), 4)
        with self.assertRaises(ValueError):
            so.physical_norb(6, True, check_norb=3)


class TestResolveFlag(unittest.TestCase):
    def test_values(self):
        for v, want in ((True, True), (False, False), (1, True), (0, False),
                        ("true", True), (" FALSE ", False), ("on", True), ("no", False)):
            with self.subTest(v=v):
                self.assertIs(so.resolve_flag({"enable_spin_orbital": v}), want)

    def test_missing_and_none(self):
        self.assertIs(so.resolve_flag({}), False)
        self.assertIs(so.resolve_flag(None), False)

    def test_key_is_case_insensitive(self):
        self.assertIs(so.resolve_flag({"Enable_Spin_Orbital": True}), True)

    def test_exact_float_zero_and_one_are_accepted(self):
        # UHFk accepted them before the switch was shared
        self.assertIs(so.resolve_flag({"enable_spin_orbital": 1.0}), True)
        self.assertIs(so.resolve_flag({"enable_spin_orbital": 0.0}), False)

    def test_a_string_value_warns(self):
        with self.assertLogs("hwave.solver.spin_orbital", "WARNING") as cm:
            self.assertIs(so.resolve_flag({"enable_spin_orbital": "true"}), True)
        self.assertIn("without quotes", cm.output[0])

    def test_unrecognized_values_raise(self):
        for v in ("ture", 2, 1.5, None):
            with self.subTest(v=v), self.assertRaisesRegex(
                    ValueError, "enable_spin_orbital"):
                so.resolve_flag({"enable_spin_orbital": v})


class TestHistoricalNames(unittest.TestCase):
    def test_rpa_reexports(self):
        import hwave.solver.rpa as rpa
        self.assertIs(rpa.validate_chi0q_index_convention,
                      so.validate_chi0q_index_convention)
        self.assertIs(rpa._so_physical_norb, so.physical_norb)


def _rpa(flag):
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as solver_rpa
    d = os.path.join(os.path.dirname(os.path.abspath(__file__)), "rpa", "input")
    info_mode = {"mode": "RPA",
                 "param": {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": 8},
                 "enable_spin_orbital": flag, "calc_scheme": "general"}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": {
        "path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}})
    return solver_rpa.RPA(io.get_param("ham"), {}, info_mode)


class TestPublicEntries(unittest.TestCase):
    def test_rpa_reads_a_false_string_as_false(self):
        # the raw value used to be kept, so "false" was truthy and RPA
        # halved the single-orbital non-SO geometry (odd count -> refused)
        solver = _rpa("false")
        self.assertIs(solver.ham_info.enable_spin_orbital, False)
        self.assertEqual(solver.ham_info.norb, 1)

    def test_rpa_refuses_a_misspelled_flag(self):
        with self.assertRaisesRegex(ValueError, "enable_spin_orbital"):
            _rpa("ture")

    def test_uhfk_odd_spin_orbital_count_raises(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.uhfk as uhfk
        d = tempfile.mkdtemp(prefix="so_odd_")
        try:
            with open(os.path.join(d, "geom.dat"), "w") as f:
                f.write("1 0 0\n0 1 0\n0 0 1\n3\n0 0 0\n0 0 0\n0 0 0\n")
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write("Transfer\n3\n1\n 1\n")
                for a in range(1, 4):
                    f.write(" 1 0 0 %d %d 1.0 0.0\n -1 0 0 %d %d 1.0 0.0\n" % (a, a, a, a))
            io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": {
                "path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}})
            mode = {"mode": "UHFk", "enable_spin_orbital": True,
                    "param": {"Ncond": 2, "T": 0.1, "CellShape": [4, 1, 1],
                              "SubShape": [1, 1, 1]}}
            with self.assertRaisesRegex(ValueError, "even Geometry norb"):
                uhfk.UHFk(io.get_param("ham"), {"print_level": 0}, mode)
        finally:
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
