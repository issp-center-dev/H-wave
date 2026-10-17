"""The mean-field tables read an off-site row (r, a, b, v) in the documented
orientation (spec 2026-09-16 D-1): after the reader's closure every off-site
displacement table is conjugate-transposed, which equals the closed table at
-r. On-site tables (r = 0) are untouched; the INFO line of D-8 is emitted
exactly when a table changed.

TestOneDirectionBondIsRefused at the end covers the READER's side of the
same subject: a bond declared in only one direction never reaches any of
the kernels above."""
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.solver import hartree_fock as hf


def _closed_only(param_ham, norb, shape):
    """The tables BEFORE the orientation step: rebuilt here from the same
    per-type recipe (closure only) so the identity new[r] == closed[-r] is
    asserted against an independent construction, not against the function
    under test."""
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


def _reversed_r(tab):
    """closed[-r]: the grid reversal of every displacement axis."""
    return np.roll(np.flip(tab, axis=(0, 1, 2)), shift=1, axis=(0, 1, 2))


#: A norb = 2 declaration whose off-site rows are REAL and ORBITAL-DIAGONAL
#: (a == b): the class the docstring of
#: :func:`hwave.solver.hartree_fock._orient_documented` declares the step
#: never moves. Two orbitals, so the transpose is a real operation and not
#: the identity for structural reasons -- a norb = 1 table would prove
#: nothing about it.
_DIAGONAL_HAM = {
    "CoulombIntra": {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5},
    "CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3,
                     ((1, 0, 0), (1, 1)): 0.25, ((-1, 0, 0), (1, 1)): 0.25},
}


def _random_density(nvol, norb, seed=7):
    """A fixed ``(nvol, 2, norb, 2, norb)`` density for ``accumulate_hf``.

    Seeded, and Hermitian in the combined ``(spin, orbital)`` index at each
    displacement separately. That is NOT the Hermiticity a physical
    real-space density obeys (which relates ``rho(r)`` to ``rho(-r)^H``),
    and it does not have to be: the comparison below is between two TABLE
    sets on the SAME density, so the density itself only has to be fixed
    and nontrivial."""
    rng = np.random.default_rng(seed)
    nd = 2 * norb
    a = (rng.normal(size=(nvol, nd, nd)) + 1j * rng.normal(size=(nvol, nd, nd)))
    a = 0.5 * (a + np.conjugate(np.swapaxes(a, -1, -2)))
    return a.reshape(nvol, 2, norb, 2, norb)


def _mean_field(inter_table, spin_table, rho, shape, norb, include_fock=True):
    nvol = int(np.prod(shape))
    out = np.zeros((nvol, 2 * norb, 2 * norb), dtype=np.complex128)
    return hf.accumulate_hf(out, rho, inter_table, spin_table, shape,
                            include_fock=include_fock)


def _unoriented(tables):
    """``tables`` with the orientation step UNDONE -- the step is an
    involution, so applying it a second time to each oriented type restores
    the closed tables the builder had before it (the same undo
    ``tests/test_flex_second_order_compat._legacy_tables`` performs)."""
    inter = dict(tables.inter_table)
    for t in hf._ORIENTED_TYPES:
        if inter.get(t) is not None:
            inter[t] = hf._orient_documented(inter[t])
    return inter


class TestOrientation(unittest.TestCase):
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

    def test_offsite_tables_equal_the_closed_table_at_minus_r(self):
        ham = self._ham()
        tabs = hf.build_interaction_tables(ham, self.norb, self.shape)
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            with self.subTest(type=t):
                got = tabs.inter_table[t]
                want = _reversed_r(ref)
                want[0, 0, 0] = ref[0, 0, 0]            # r = 0 untouched
                np.testing.assert_array_equal(got, want)
                # non-vacuity: the orientation step changed something off-site
                self.assertFalse(np.array_equal(got, ref), t)

    def test_onsite_tables_are_untouched_including_complex_pairhop(self):
        ham = self._ham()
        tabs = hf.build_interaction_tables(ham, self.norb, self.shape)
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            np.testing.assert_array_equal(tabs.inter_table[t][0, 0, 0], ref[0, 0, 0], t)
        np.testing.assert_array_equal(tabs.inter_table["CoulombIntra"][0, 0, 0],
                                      np.diag([2.0, 1.5]).astype(complex))

    def test_orient_documented_is_an_involution_on_closed_tables(self):
        ham = self._ham()
        closed = _closed_only(ham, self.norb, self.shape)
        for t, ref in closed.items():
            twice = hf._orient_documented(hf._orient_documented(ref))
            np.testing.assert_array_equal(twice, ref, t)

    def test_real_single_orbital_and_real_diagonal_rows_are_invariant(self):
        """Both halves of the invariance claim in ``_orient_documented``'s
        docstring: a real SINGLE-ORBITAL table, and a real ORBITAL-DIAGONAL
        table at norb = 2 -- where the orbital transpose is a real operation
        that is not the identity for structural reasons, so the invariance
        is a statement about the CONTENT of the table and not about its
        shape. The second half also pins the silence of D-8's INFO line:
        nothing moved, so the builder must say nothing."""
        ham = {"CoulombInter": {((1, 0, 0), (0, 0)): 0.3, ((-1, 0, 0), (0, 0)): 0.3}}
        tabs = hf.build_interaction_tables(ham, 1, self.shape)
        closed = _closed_only(ham, 1, self.shape)
        np.testing.assert_array_equal(tabs.inter_table["CoulombInter"], closed["CoulombInter"])

        # norb = 2, orbital-diagonal off-site rows only
        closed2 = _closed_only(_DIAGONAL_HAM, self.norb, self.shape)["CoulombInter"]
        self.assertGreater(np.abs(closed2).max(), 1e-12)            # anti-vacuity
        np.testing.assert_array_equal(hf._orient_documented(closed2), closed2)
        logger = logging.getLogger("qlms.solver.hartree_fock")
        with self.assertNoLogs(logger, level=logging.INFO):
            tabs2 = hf.build_interaction_tables(_DIAGONAL_HAM, self.norb, self.shape)
        np.testing.assert_array_equal(tabs2.inter_table["CoulombInter"], closed2)

    def test_complex_single_orbital_pairhop_is_not_invariant(self):
        ham = {"PairHop": {((1, 0, 0), (0, 0)): 0.1 + 0.05j, ((-1, 0, 0), (0, 0)): 0.1 - 0.05j}}
        tabs = hf.build_interaction_tables(ham, 1, self.shape)
        closed = _closed_only(ham, 1, self.shape)
        self.assertFalse(np.array_equal(tabs.inter_table["PairHop"], closed["PairHop"]))
        np.testing.assert_array_equal(tabs.inter_table["PairHop"], _reversed_r(closed["PairHop"]))

    def test_real_orbital_diagonal_mean_field_does_not_move(self):
        """The declared "never moves" class pinned where a user sees it: the
        MEAN FIELD, not just the table.

        ``accumulate_hf`` on a fixed density gives the same matrix from the
        oriented tables and from the un-oriented ones, bit for bit, on a
        real orbital-diagonal declaration. Taken alone that leg is close to
        a tautology -- on this class the orientation IS the identity, so the
        two table sets are the same arrays and ``accumulate_hf`` cannot tell
        them apart. What the method actually pins is the pair: that the mean
        field of the declared-invariant class does not move WHILE the same
        comparison on the INTER-ORBITAL declaration of ``_ham`` does, which
        is the control on the next lines."""
        nvol = int(np.prod(self.shape))
        rho = _random_density(nvol, self.norb)

        tabs = hf.build_interaction_tables(_DIAGONAL_HAM, self.norb, self.shape)
        oriented = _mean_field(tabs.inter_table, tabs.spin_table, rho,
                               self.shape, self.norb)
        legacy = _mean_field(_unoriented(tabs), tabs.spin_table, rho,
                             self.shape, self.norb)
        self.assertGreater(np.abs(oriented).max(), 1e-12)            # anti-vacuity
        np.testing.assert_array_equal(oriented, legacy)

        # control: the inter-orbital declaration DOES move under the same
        # comparison
        tabs_i = hf.build_interaction_tables(self._ham(), self.norb, self.shape)
        oriented_i = _mean_field(tabs_i.inter_table, tabs_i.spin_table, rho,
                                 self.shape, self.norb)
        legacy_i = _mean_field(_unoriented(tabs_i), tabs_i.spin_table, rho,
                               self.shape, self.norb)
        moved = np.abs(oriented_i - legacy_i).max() / np.abs(oriented_i).max()
        self.assertGreater(
            moved, 1e-3,
            "the orientation step does not move the mean field of the "
            "inter-orbital declaration either ({:.3e} relative), so the "
            "equality above says nothing".format(moved))

    def test_aggregate_coulomb_declaration_is_oriented_too(self):
        """The aggregate ``Coulomb`` type reaches the orientation step.

        ``build_interaction_tables`` is where ``Coulomb`` is split
        (``wan90.split_coulomb``: the ``r = 0`` orbital-diagonal entries are
        ``CoulombIntra``, everything else ``CoulombInter``), i.e. the split
        happens INSIDE the function whose last step orients. The orientation
        loop itself is shared, not per-route, so this is not a route-specific
        step being checked: what the method pins is the COMPOSITION of the
        split with it in the aggregate branch -- a branch that builds its own
        ``uab_r``/``vab_r`` and fills ``inter_table`` separately from the
        explicit one, and could therefore drift from it (a different
        closure, a table written after the loop) without any other test
        noticing. A user who declares one aggregate file must get exactly
        what the explicitly split declaration gives, orientation included,
        which is what the second leg checks."""
        intra = {((0, 0, 0), (0, 0)): 2.0, ((0, 0, 0), (1, 1)): 1.5}
        inter = {((1, 0, 0), (0, 1)): 0.4, ((-1, 0, 0), (1, 0)): 0.4,
                 ((0, 1, 0), (0, 1)): 0.2, ((0, -1, 0), (1, 0)): 0.2}
        aggregate = dict(intra)
        aggregate.update(inter)
        got = hf.build_interaction_tables({"Coulomb": aggregate},
                                          self.norb, self.shape)
        split = hf.build_interaction_tables(
            {"CoulombIntra": intra, "CoulombInter": inter}, self.norb, self.shape)
        for t in ("CoulombIntra", "CoulombInter"):
            np.testing.assert_array_equal(got.inter_table[t],
                                          split.inter_table[t], t)
        # and this declaration really is one the orientation step moves, so
        # the equality above is not two un-oriented tables agreeing
        closed = _closed_only({"CoulombInter": inter}, self.norb, self.shape)
        self.assertFalse(
            np.array_equal(got.inter_table["CoulombInter"], closed["CoulombInter"]),
            "the aggregate Coulomb table was not oriented")
        np.testing.assert_array_equal(
            got.inter_table["CoulombInter"],
            hf._orient_documented(closed["CoulombInter"]))

    def test_info_line_only_when_a_table_changed(self):
        onsite = {"CoulombIntra": {((0, 0, 0), (0, 0)): 2.0},
                  "CoulombInter": {((0, 0, 0), (0, 1)): 0.7, ((0, 0, 0), (1, 0)): 0.7}}
        logger = logging.getLogger("qlms.solver.hartree_fock")
        records = []
        handler = logging.Handler()
        handler.setLevel(logging.INFO)
        handler.emit = records.append
        logger.addHandler(handler)
        old_level = logger.level
        logger.setLevel(logging.INFO)
        try:
            hf.build_interaction_tables(onsite, 2, self.shape)
            self.assertFalse([r for r in records if "documented orientation" in r.getMessage()])
            records.clear()
            hf.build_interaction_tables(self._ham(), 2, self.shape)
            hits = [r for r in records if "documented orientation" in r.getMessage()]
            self.assertEqual(len(hits), 1)
            # (type, r) entries changed: CoulombInter +-x (2), Hund 2, Ising 2, PairLift 2,
            # Exchange +-y (2), PairHop +-x (2); the +-y CoulombInter rows are real and
            # orbital-diagonal -> unchanged
            self.assertIn("12 off-site displacement table(s)", hits[0].getMessage())
            # the user-facing tail, verbatim: it is what points a 2.0.0 user at
            # the migration text, so it is part of the contract of this INFO
            self.assertIn("differ from H-wave 2.0.0 -- see the release note",
                          hits[0].getMessage())
            records.clear()
            diag = {"PairHop": {((1, 0, 0), (0, 0)): 0.1 + 0.05j, ((-1, 0, 0), (0, 0)): 0.1 - 0.05j}}
            hf.build_interaction_tables(diag, 1, self.shape)
            hits = [r for r in records if "documented orientation" in r.getMessage()]
            self.assertEqual(len(hits), 1)
            self.assertIn("2 off-site displacement table(s)", hits[0].getMessage())
        finally:
            logger.removeHandler(handler)
            logger.setLevel(old_level)


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
