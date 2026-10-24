"""A nonzero on-site CoulombInter row of an orbital with itself is refused
at read time, in both readers (the k-space reader used by UHFk/RPA/FLEX and
hwave_sc's own reader): it is a Hubbard-type term that belongs in
CoulombIntra, and every solver used to treat it incompletely."""
import os
import shutil
import tempfile
import unittest

from hwave.solver.declarations import validate_no_onsite_same_orbital_inter


def _write(d, rows, norb=2):
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % norb + "0 0 0\n" * norb)
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("Transfer\n%d\n1\n 1\n" % norb)
        for a in range(1, norb + 1):
            f.write(" 1 0 0 %d %d -1.0 0.0\n -1 0 0 %d %d -1.0 0.0\n" % (a, a, a, a))
    with open(os.path.join(d, "coulombinter.dat"), "w") as f:
        f.write("CoulombInter\n%d\n1\n 1\n" % norb)
        for r in rows:
            f.write(" %d %d %d %d %d %.6f 0.0\n" % r)


_SAME = [(0, 0, 0, 1, 1, 0.5)]
_FINE = [(0, 0, 0, 1, 2, 0.5), (0, 0, 0, 2, 1, 0.5),          # on-site, different orbitals
         (1, 0, 0, 1, 1, 0.3), (-1, 0, 0, 1, 1, 0.3),         # off-site, same orbital
         (0, 0, 0, 2, 2, 0.0)]                                  # declared at zero


class TestValidator(unittest.TestCase):
    def test_refuses_the_nonzero_row_and_names_it(self):
        # the message gives the conversion: U = V in CoulombIntra plus V/2
        with self.assertRaisesRegex(ValueError, r"orbital 1 with itself.*U = V.*V/2.*CoulombIntra"):
            validate_no_onsite_same_orbital_inter(
                "CoulombInter", {((0, 0, 0), (0, 0)): 0.5}, source="x.dat")

    def test_accepts_other_rows_and_other_types(self):
        validate_no_onsite_same_orbital_inter(
            "CoulombInter", {((0, 0, 0), (0, 1)): 0.5, ((1, 0, 0), (0, 0)): 0.3,
                             ((0, 0, 0), (1, 1)): 0.0})
        validate_no_onsite_same_orbital_inter("Hund", {((0, 0, 0), (0, 0)): 0.5})


class TestReaders(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="samev_")

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _k_reader(self, key="CoulombInter"):
        import hwave.qlmsio.read_input_k as read_input_k
        return read_input_k.QLMSkInput({"path_to_input": self.d, "interaction": {
            "path_to_input": self.d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
            key: "coulombinter.dat"}})

    def _sc_reader(self, key="CoulombInter"):
        from hwave import sc
        return sc._read_interaction_files({"file": {"input": {"interaction": {
            "path_to_input": self.d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
            key: "coulombinter.dat"}}}})

    def test_the_key_is_case_insensitive(self):
        _write(self.d, _SAME + [(0, 0, 0, 1, 2, 0.2), (0, 0, 0, 2, 1, 0.2)])
        for reader in (self._k_reader, self._sc_reader):
            with self.subTest(reader=reader.__name__):
                with self.assertRaisesRegex(ValueError, "CoulombIntra"):
                    reader(key="coulombinter")

    def test_the_coulomb_aggregate_routes_the_row_to_coulombintra(self):
        """In the combined ``Coulomb`` file the same row IS the on-site U:
        it is split off into CoulombIntra and accepted, in both readers."""
        _write(self.d, _SAME + [(0, 0, 0, 1, 2, 0.2), (0, 0, 0, 2, 1, 0.2)])
        r = self._k_reader(key="Coulomb")
        # the k-space consumers split the aggregate the same way (wan90.split_coulomb)
        from hwave.qlmsio import wan90
        intra, inter_k = wan90.split_coulomb(r.get_param("ham")["Coulomb"])
        self.assertIn(((0, 0, 0), (0, 0)), intra)
        self.assertNotIn(((0, 0, 0), (0, 0)), inter_k)
        _, _, inter = self._sc_reader(key="Coulomb")
        self.assertIn(((0, 0, 0), (0, 0)), inter["CoulombIntra"])

    def test_both_readers_refuse(self):
        _write(self.d, _SAME + [(0, 0, 0, 1, 2, 0.2), (0, 0, 0, 2, 1, 0.2)])
        for reader in (self._k_reader, self._sc_reader):
            with self.subTest(reader=reader.__name__):
                with self.assertRaisesRegex(ValueError, "CoulombIntra"):
                    reader()

    def test_both_readers_accept_the_other_rows(self):
        _write(self.d, _FINE)
        self._k_reader()
        self._sc_reader()


if __name__ == "__main__":
    unittest.main()
