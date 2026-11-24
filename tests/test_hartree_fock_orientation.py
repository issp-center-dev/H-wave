"""The mean-field tables keep a declared off-site row where the reader puts
it: ``inter_table[t][r, a, b]`` is the closed declaration ``(X_ab(r) +
conj(X_ba(-r))) / 2`` itself, with no further reversal (issue #232).

A row ``(r, a, b)`` places orbital ``a`` in the original cell and ``b`` in
the cell displaced by ``r`` -- the bond of a Transfer row ``(r, a, b)``
(Wannier90 / RESPACK). That this table layout realizes that bond in the
mean field is checked physically, against an independent real-space
Hartree-Fock solver, in ``tests/test_two_body_placement.py``; this module
only pins the table layout, so that a reversal step cannot creep back in.
(An earlier revision conjugate-transposed every off-site table here, which
moved each two-body row to the bond opposite its Transfer row.)

TestOneDirectionBondIsRefused at the end covers the READER's side of the
same subject: a bond declared in only one direction never reaches any of
the kernels above."""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver import hartree_fock as hf


def _closed_only(param_ham, norb, shape):
    """The closed tables rebuilt here from the same per-type recipe, so the
    builder is compared against an independent construction."""
    nx, ny, nz = shape
    out = {}
    sign = {"Hund": -1.0, "Exchange": -1.0}
    for t in ("CoulombInter", "Hund", "Ising", "PairLift", "Exchange", "PairHop"):
        if t not in param_ham:
            continue
        tab = np.zeros((nx, ny, nz, norb, norb), dtype=np.complex128)
        for (irvec, orbvec), v in param_ham[t].items():
            tab[(*irvec, *orbvec)] += v
        out[t] = sign.get(t, 1.0) * hf._reverse_closed(tab)
    return out


class TestTablesAreTheClosedDeclaration(unittest.TestCase):
    shape = (3, 3, 1)          # every axis odd: no displacement is its own reverse
    norb = 2

    def _ham(self):
        return {
            "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
            "CoulombInter": {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                             ((0, 1, 0), (0, 0)): 0.2, ((0, -1, 0), (0, 0)): 0.2,
                             ((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7},
            "Hund": {((1, 0, 0), (0, 1)): 0.1, ((-1, 0, 0), (1, 0)): 0.1},
            "Ising": {((1, 0, 0), (0, 1)): 0.15, ((-1, 0, 0), (1, 0)): 0.15},
            "PairLift": {((1, 0, 0), (0, 1)): 0.05, ((-1, 0, 0), (1, 0)): 0.05},
            "Exchange": {((0, 1, 0), (0, 1)): 0.12, ((0, -1, 0), (1, 0)): 0.12},
            "PairHop": {((1, 0, 0), (0, 1)): 0.08 + 0.03j, ((-1, 0, 0), (1, 0)): 0.08 - 0.03j,
                        ((0, 0, 0), (0, 1)): 0.2 + 0.1j, ((0, 0, 0), (1, 0)): 0.2 - 0.1j},
        }

    def test_every_table_is_the_closed_declaration(self):
        ham = self._ham()
        tabs = hf.build_interaction_tables(ham, self.norb, self.shape)
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            with self.subTest(type=t):
                np.testing.assert_array_equal(tabs.inter_table[t], ref)

    def test_the_fixture_distinguishes_a_reversal(self):
        """Non-vacuity: every off-site table here differs from its own
        reversal r -> -r, so a reintroduced reversal step would be caught."""
        closed = _closed_only(self._ham(), self.norb, self.shape)
        for t, ref in closed.items():
            rev = np.roll(np.flip(ref, axis=(0, 1, 2)), shift=1, axis=(0, 1, 2))
            with self.subTest(type=t):
                self.assertFalse(np.allclose(rev, ref), t)

    def test_aggregate_coulomb_equals_the_explicit_split(self):
        """The aggregate ``Coulomb`` table lands in the same tables as the
        explicit CoulombIntra + CoulombInter declaration of the same rows."""
        rows = {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5,
                ((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                ((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7}
        agg = hf.build_interaction_tables({"Coulomb": dict(rows)}, self.norb, self.shape)
        intra = {k: v for k, v in rows.items() if k[0] == (0, 0, 0) and k[1][0] == k[1][1]}
        inter = {k: v for k, v in rows.items() if k not in intra}
        exp = hf.build_interaction_tables({"CoulombIntra": intra, "CoulombInter": inter},
                                          self.norb, self.shape)
        for t in ("CoulombIntra", "CoulombInter"):
            with self.subTest(type=t):
                np.testing.assert_array_equal(agg.inter_table[t], exp.inter_table[t])


class TestOneDirectionBondIsRefused(unittest.TestCase):
    """A two-body bond declared in ONE direction only is rejected at READ
    time, through the public entry point, before any kernel sees it.

    WHY THIS BELONGS BESIDE THE ORIENTATION TESTS. A one-direction table
    looks like the perfect orientation experiment: declare ``v_12(+x)`` and
    nothing else, and see which cell orbital 1 lands in. It is not an
    admissible counterexample, and this test is what says so. The file
    format's two entries ``X_ab(R)`` and ``X_ba(-R)`` are not two couplings
    but ONE written twice (``X_ab(R) = conj(X_ba(-R))``), so a table with
    only one of them does not denote half a bond -- it denotes nothing the
    format can express, and ``declarations.validate_hermitian_closure``
    (issue #93) refuses it rather than letting each solver complete it in
    its own way. Every orientation statement in this module and in the
    modules it names is therefore about which reading a HERMITIAN-CLOSED
    declaration gets; a test that tried to settle it with a lone row would
    fail at the reader and prove nothing about the kernels.

    Driven through ``read_input_k.QLMSkInput`` on a temporary directory --
    the reader's own entry point, not the validator called directly --
    because "before the kernels run" is a property of where the check sits
    in the pipeline, which only the public entry can show."""

    _GEOM = ("  1.000000000000   0.000000000000   0.000000000000\n"
             "  0.000000000000   1.000000000000   0.000000000000\n"
             "  0.000000000000   0.000000000000   1.000000000000\n"
             "2\n"
             "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
             "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n")

    _TRANSFER = ("Transfer in wannier90-like format for uhfk\n2\n2\n1 1\n"
                 "   1    0    0    1    1 -1.000000000000000e+00  0.000000000000000e+00\n"
                 "  -1    0    0    1    1 -1.000000000000000e+00  0.000000000000000e+00\n"
                 "   1    0    0    2    2 -1.000000000000000e+00  0.000000000000000e+00\n"
                 "  -1    0    0    2    2 -1.000000000000000e+00  0.000000000000000e+00\n")

    #: the SAME inter-orbital bond the orientation fixtures carry, with the
    #: -x partner row simply left out
    _ONE_DIRECTION = ("CoulombInter in wannier90-like format for uhfk\n2\n1\n1\n"
                      "   1    0    0    1    2  4.000000000000000e-01"
                      "  0.000000000000000e+00\n")

    #: the closed version of it, which must be ACCEPTED -- the control that
    #: says the refusal is about the missing partner and not about the file
    _BOTH_DIRECTIONS = _ONE_DIRECTION.replace("\n2\n1\n1\n", "\n2\n2\n1 1\n") + (
        "  -1    0    0    2    1  4.000000000000000e-01  0.000000000000000e+00\n")

    def _read(self, inter_text):
        import hwave.qlmsio.read_input_k as read_input_k
        path = tempfile.mkdtemp(prefix="hwave_one_direction_")
        try:
            with open(os.path.join(path, "geom.dat"), "w") as fw:
                fw.write(self._GEOM)
            with open(os.path.join(path, "transfer.dat"), "w") as fw:
                fw.write(self._TRANSFER)
            with open(os.path.join(path, "coulombinter.dat"), "w") as fw:
                fw.write(inter_text)
            idict = {"path_to_input": path, "Geometry": "geom.dat",
                     "Transfer": "transfer.dat", "CoulombInter": "coulombinter.dat"}
            return read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
        finally:
            shutil.rmtree(path, ignore_errors=True)

    def test_a_bond_declared_in_one_direction_is_refused_by_the_reader(self):
        with self.assertRaises(ValueError) as cm:
            self._read(self._ONE_DIRECTION)
        message = str(cm.exception)
        self.assertIn("CoulombInter", message)
        self.assertIn("coulombinter.dat", message)            # names the file
        # names the MISSING PARTNER: the reversed displacement and the
        # swapped orbital pair, in the file's own 1-based numbering
        self.assertIn("NO partner entry", message)
        self.assertIn("R=(-1, 0, 0)", message)
        self.assertIn("orbitals (2, 1)", message)
        # and names the row that was declared, so the operator can find it
        self.assertIn("R=(1, 0, 0) orbitals (1, 2)", message)

    def test_the_closed_declaration_is_accepted(self):
        """The control: the same bond WITH its partner row reads fine, and
        carries both entries. Without this the refusal above could be any
        defect in the temporary fixture."""
        reader = self._read(self._BOTH_DIRECTIONS)
        table = reader.get_param("ham")["CoulombInter"]
        self.assertIn(((1, 0, 0), (0, 1)), table)
        self.assertIn(((-1, 0, 0), (1, 0)), table)
        self.assertAlmostEqual(abs(complex(table[((1, 0, 0), (0, 1))])), 0.4, places=12)


if __name__ == "__main__":
    unittest.main()
