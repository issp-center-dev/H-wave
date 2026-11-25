"""A coefficient that is not a finite number is refused when it is read.

``nan``, ``inf`` and a literal too large for a float are refused by the
reader of the Wannier90-like files (every file of that format: Transfer,
Extern and the two-body definitions, for UHFk, RPA, FLEX and hwave_sc) and
by the reader of the UHFr definition files, with the file and the line. A
one-body table built in memory is refused by the shared H(k) builder.

Before, UHFk read such a Transfer file, dropped the bond in its block
detection and reported a converged run of another Hamiltonian.
"""
import logging
import os
import tempfile
import unittest

import numpy as np
from requests.structures import CaseInsensitiveDict

from hwave.qlmsio import wan90
from hwave.qlmsio.read_input import QLMSInput
from hwave.qlmsio.read_input_k import QLMSkInput
from hwave.solver import hamiltonian_k
from hwave.solver.uhfk import UHFk


# one orbital, a chain along x: three R-points, degeneracy 1
_HEAD = ["Transfer in wannier90-like format", "1", "3", " 1 1 1"]
_ROWS = ["  -1    0    0    1    1  -1.000000000000  -0.000000000000",
         "   0    0    0    1    1   0.500000000000   0.000000000000",
         "   1    0    0    1    1  -1.000000000000   0.000000000000"]

_SPIN_HEAD = ["Transfer in wannier90-like format w90spin/1", "1", "1", " 1"]
_SPIN_ROWS = ["   0    0    0    1    1    0    0   0.500000000000   0.000000000000",
              "   0    0    0    1    1    1    1  -0.500000000000   0.000000000000"]


def _write(d, name, lines):
    path = os.path.join(d, name)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _with(rows, index, column, token):
    """``rows`` with one token of row ``index`` replaced (``column`` counts
    from the end: -2 the real part, -1 the imaginary part)."""
    out = list(rows)
    fields = out[index].split()
    fields[column] = token
    out[index] = "  ".join(fields)
    return out


class TestWannier90LikeReader(unittest.TestCase):

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name

    def test_a_finite_file_is_read_as_before(self):
        path = _write(self.d, "transfer.dat", _HEAD + _ROWS)
        data = wan90.read_w90(path)
        self.assertEqual(data, {((-1, 0, 0), (0, 0)): complex(-1.0, -0.0),
                                ((0, 0, 0), (0, 0)): complex(0.5, 0.0),
                                ((1, 0, 0), (0, 0)): complex(-1.0, 0.0)})

    def test_every_non_finite_token_is_refused_with_file_and_line(self):
        # the body starts on line 5; row index 2 is line 7
        for token in ("nan", "NaN", "inf", "-inf", "Infinity", "1e999"):
            for column, part in ((-2, "real"), (-1, "imaginary")):
                with self.subTest(token=token, part=part):
                    path = _write(self.d, "transfer.dat",
                                  _HEAD + _with(_ROWS, 2, column, token))
                    with self.assertRaises(ValueError) as cm:
                        wan90.read_w90(path)
                    message = str(cm.exception)
                    self.assertIn("transfer.dat: line 7: non-finite coefficient: ",
                                  message)
                    self.assertIn(token, message)

    def test_the_first_offending_line_is_named(self):
        rows = _with(_with(_ROWS, 0, -2, "nan"), 2, -2, "inf")
        path = _write(self.d, "transfer.dat", _HEAD + rows)
        with self.assertRaisesRegex(
                ValueError, r"transfer\.dat: line 5: non-finite coefficient: nan "):
            wan90.read_w90(path)

    def test_the_spin_resolved_format(self):
        ok = _write(self.d, "ok.dat", _SPIN_HEAD + _SPIN_ROWS)
        self.assertEqual(len(wan90.read_w90(ok, allow_spin_columns=True)), 2)
        bad = _write(self.d, "bad.dat", _SPIN_HEAD + _with(_SPIN_ROWS, 1, -2, "nan"))
        with self.assertRaisesRegex(
                ValueError, r"bad\.dat: line 6: non-finite coefficient: nan "):
            wan90.read_w90(bad, allow_spin_columns=True)


class TestKSpaceInput(unittest.TestCase):
    """Through ``QLMSkInput``, the reader of UHFk, RPA and FLEX: Transfer,
    Extern and a two-body file alike."""

    GEOM = ["1.0 0.0 0.0", "0.0 1.0 0.0", "0.0 0.0 1.0", "1", "0.0 0.0 0.0"]
    ONSITE = ["CoulombIntra in wannier90-like format", "1", "1", " 1",
              "   0    0    0    1    1   4.000000000000   0.000000000000"]

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name
        _write(self.d, "geom.dat", self.GEOM)

    def _read(self, **files):
        info = {"interaction": dict({"path_to_input": self.d, "Geometry": "geom.dat"},
                                    **files)}
        return QLMSkInput(info)

    def test_finite_files_are_read(self):
        _write(self.d, "transfer.dat", _HEAD + _ROWS)
        _write(self.d, "coulombintra.dat", self.ONSITE)
        read = self._read(Transfer="transfer.dat", CoulombIntra="coulombintra.dat")
        self.assertEqual(len(read.get_param("ham")["Transfer"]), 3)

    def test_each_file_with_a_nan_is_refused(self):
        _write(self.d, "transfer.dat", _HEAD + _ROWS)
        _write(self.d, "coulombintra.dat", self.ONSITE)
        _write(self.d, "bad_transfer.dat", _HEAD + _with(_ROWS, 0, -2, "nan"))
        _write(self.d, "bad_extern.dat", _HEAD + _with(_ROWS, 1, -2, "nan"))
        _write(self.d, "bad_coulombintra.dat", _with(self.ONSITE, 4, -2, "nan"))
        for key, name, line in (("Transfer", "bad_transfer.dat", 5),
                                ("Extern", "bad_extern.dat", 6),
                                ("CoulombIntra", "bad_coulombintra.dat", 5)):
            with self.subTest(file=key):
                files = {"Transfer": "transfer.dat", key: name}
                with self.assertRaisesRegex(
                        ValueError, r"{}: line {}: non-finite coefficient: nan ".format(
                            name.replace(".", r"\."), line)):
                    self._read(**files)


class TestRealSpaceDefinitionFiles(unittest.TestCase):
    """The reader of the UHFr ``.def`` files stops with the file and the
    line, in its own way of reporting (an error message and ``exit(1)``)."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.d = self._tmp.name

    def _trans(self, value):
        return ["========================", "NTransfer      2",
                "========================", "========i_j_s_tijs======",
                "========================",
                "    1     0     0     0         1.000000000000000         0.000000000000000",
                "    0     0     1     0         {}        -0.000000000000000".format(value)]

    def _intra(self, value):
        return ["=============================================", "NCoulombIntra      2",
                "=============================================",
                "================== CoulombIntra ================",
                "=============================================",
                "    0         8.000000000000000",
                "    1         {}".format(value)]

    def test_finite_files_are_read(self):
        files = CaseInsensitiveDict({
            "Trans": _write(self.d, "trans.def", self._trans("1.000000000000000")),
            "CoulombIntra": _write(self.d, "coulombintra.def",
                                   self._intra("8.000000000000000"))})
        ham = QLMSInput(files).get_param("ham")
        self.assertEqual(ham["Transfer"][(0, 0, 1, 0)], complex(1.0, -0.0))
        self.assertEqual(ham["CoulombIntra"][(1,)], 8.0)

    def test_a_non_finite_value_stops_the_read(self):
        for key, name, lines in (
                ("Trans", "trans.def", self._trans("nan")),
                ("Trans", "trans.def", self._trans("inf")),
                ("CoulombIntra", "coulombintra.def", self._intra("nan"))):
            with self.subTest(file=name, lines=lines[-1].split()[-2:]):
                path = _write(self.d, name, lines)
                with self.assertLogs("qlms", level=logging.ERROR) as logs:
                    with self.assertRaises(SystemExit) as cm:
                        QLMSInput(CaseInsensitiveDict({key: path}))
                self.assertEqual(cm.exception.code, 1)
                message = "\n".join(logs.output)
                self.assertIn("{}: line 7: non-finite value: ".format(path), message)


class TestTablesBuiltInMemory(unittest.TestCase):
    """A table that did not come from a file is refused by the shared H(k)
    builder, which UHFk, RPA and FLEX go through."""

    TABLE = {((1, 0, 0), (0, 0)): -1.0, ((-1, 0, 0), (0, 0)): -1.0,
             ((0, 1, 0), (0, 0)): -1.0, ((0, -1, 0), (0, 0)): -1.0}

    def test_dense_table_refuses_a_non_finite_coefficient(self):
        hamiltonian_k.dense_table(self.TABLE, (4, 4, 1), 1)
        for value in (float("nan"), float("inf"), complex(0.0, float("nan"))):
            with self.subTest(value=value):
                table = dict(self.TABLE)
                table[((1, 0, 0), (0, 0))] = value
                with self.assertRaisesRegex(
                        ValueError,
                        r"Transfer: non-finite coefficient .* of orbitals "
                        r"\(1, 1\) at R = \(1, 0, 0\)"):
                    hamiltonian_k.dense_table(table, (4, 4, 1), 1)

    def test_an_entry_that_is_dropped_is_not_judged(self):
        # the spin block of a spin-orbital-format table in normal mode
        table = dict(self.TABLE)
        table[((0, 0, 0), (1, 1))] = float("nan")
        logging.disable(logging.WARNING)
        self.addCleanup(logging.disable, logging.NOTSET)
        tab = hamiltonian_k.dense_table(table, (4, 4, 1), 1, out_of_range="drop")
        self.assertTrue(np.all(np.isfinite(tab)))

    def test_uhfk_refuses_instead_of_converging_without_the_bond(self):
        """One orbital on 4 x 4 with ``SubShape = [2, 2, 1]`` and a NaN
        hopping along x. Before, the block detection read the NaN bond as
        no bond, the Hamiltonian was cleared of it, and the run converged
        (total energy -7.05 where the finite model gives -12.57)."""
        nan = float("nan")
        ham = {"Geometry": {"norb": 1, "rvec": np.eye(3),
                            "center": np.zeros((1, 3)), "degree": np.ones(1)},
               "Transfer": {((1, 0, 0), (0, 0)): nan, ((-1, 0, 0), (0, 0)): nan,
                            ((0, 1, 0), (0, 0)): -1.0, ((0, -1, 0), (0, 0)): -1.0},
               "CoulombIntra": {((0, 0, 0), (0, 0)): 4.0}}
        params = {"CellShape": [4, 4, 1], "SubShape": [2, 2, 1], "Ncond": 16,
                  "2Sz": None, "IterationMax": 300, "eps": 1e-8, "Mix": 0.5,
                  "RndSeed": 123456789, "T": 0.0}
        with self.assertRaisesRegex(ValueError, "Transfer: non-finite coefficient"):
            solver = UHFk(ham, {"print_level": 0}, {"mode": "UHFk"}, params)
            solver.solve({"initial_mode": "random"}, ".")


if __name__ == "__main__":
    unittest.main()
