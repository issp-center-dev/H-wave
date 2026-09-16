"""Spec 2026-09-16 section 2.4, gate G-UHFk: UHFk END TO END on the
reversed declaration.

An off-site two-body row ``(r, a, b, v)`` is read in the DOCUMENTED
orientation since issue #192 -- ``v n_{j,a} n_{j+r,b}``, orbital ``a`` in
the original cell and ``b`` in the cell displaced by ``r``. The previous
release (H-wave 2.0.0, the reference revision pinned by
``tests/test_flex_second_order_compat.DEVELOP_COMMIT``) read the two cells
the other way round. The compatibility promise that goes with the change
is therefore an EXACT one:

    this branch's UHFk on a declaration ``F``
        == the reference revision's UHFk on ``F^rev``

where ``F^rev`` negates the displacement of every off-site two-body row
and leaves the orbital indices and the on-site rows alone
(``tests/test_flex_hf_scf._reversed_interaction_dir`` writes it). The
in-process, table-level version of the same statement is
``tests/test_flex_second_order_compat`` (``new(F) == legacy(F^rev)``, at
``np.array_equal``); this module is the end-to-end one -- two real
``hwave.qlms.run`` processes, the two source trees, the shipped
``energy.dat`` and ``green.dat.npz``.

WHY THE COMPARISON IS ``np.array_equal`` AND NOT A TOLERANCE
------------------------------------------------------------
Because the two runs do the SAME arithmetic. The orientation step is a
conjugate transpose at fixed ``r`` of the assembled displacement table,
and the reversal of the declaration produces that same table out of the
reference revision's builder -- bit for bit, not to round-off (Task 2's
kernel test pins exactly that with ``np.array_equal``). The tables being
identical, the whole SCF trajectory is identical: same initial Green
(``RndSeed`` is pinned), same mixing, same eigen-decompositions, same
iteration count. A tolerance here would hide a real difference in the
mean field. The only thing that could separate the two processes is a
different numpy/BLAS on the two ``PYTHONPATH``s, which is not a thing on
one machine; if a member ever does differ at round-off, the comparison
below falls back to 1e-12 relative and SAYS which member it was, rather
than quietly widening for everything.

WHICH CASES THE ORIENTATION ACTUALLY MOVES (and why the rest cannot)
--------------------------------------------------------------------
The identity above is asserted for every case. Its anti-vacuity control
-- the reference revision on the UNREVERSED ``F``, which must land
somewhere else -- is only meaningful where the orientation is VISIBLE at
all, and two documented classes of case it is not. Both are structural,
derived rather than observed, and asserted here in the coinciding
direction so that a change making them visible cannot pass unnoticed:

* **Hartree-only** (``flag_fock = false``) for ``CoulombInter`` / ``Hund``
  / ``Ising``: the Hartree term contracts the table only through
  ``sum_r J_ab(r)`` (``hartree_fock.accumulate_hf``: ``hh1 =
  einsum('rab, stb -> rsta', jab_r, hh0)`` then a sum over ``r``). The
  reversal/Hermitian closure ``_reverse_closed`` makes that sum
  HERMITIAN, and the orientation maps it to its own conjugate transpose,
  i.e. to itself. Nothing downstream can see the change. ``PairHop`` is
  NOT in this class: its term (the ``hh6``/``hh7`` block) is not gated by
  ``flag_fock`` and reads the displacement resolved, so the off-site
  complex PairHop moves with the Fock term on or off.

* **A single orbital**: at ``norb = 1`` the transpose is the identity, so
  the oriented table is the COMPLEX CONJUGATE of the original one. With a
  real transfer, ``F^rev``'s whole Hamiltonian is then the complex
  conjugate of ``F``'s, and every real observable ``energy.dat`` carries
  is invariant under complex conjugation -- a single-orbital chain cannot
  make the orientation visible in the energy no matter what phase the
  off-site amplitude carries (the complex PairHop amplitude here does
  move the TABLE, which is what Task 2's kernel gate sees). The
  fixture is kept because the end-to-end identity is worth pinning on the
  single-orbital path too, and because the statement "this one cannot
  move" is itself falsifiable.

The fixtures under ``tests/uhfk_orientation/`` are a ``CellShape =
[4, 1, 1]`` chain with COMPLEX inter-orbital hopping and ``t_01(+x) !=
conj(t_10(+x))`` (the amplitudes of
``tests/test_flex_second_order_ed_chain._fx_complex``), so the band is
neither inversion- nor transpose-symmetric and the two readings of a bond
really are different Hamiltonians. ``coulombinter.dat`` carries an
ASYMMETRIC inter-orbital bond (``v_12(+x) = 0.4 != v_21(+x) = 0.25``),
``hund.dat`` an inter-orbital ``J``, ``pairhop.dat`` a COMPLEX off-site
amplitude; every case also runs an on-site ``coulombintra.dat`` so the
mean field is nontrivial. The ``*_so`` directories hold the same band
written in spin-orbital indices (``so = 2a + s``, spin-diagonal) for the
``enable_spin_orbital = true`` twin of each case: UHFk supports every
two-body type there (``uhfk._check_spin_orbital_compatibility``), the
interaction files keep PHYSICAL orbital indices, and ``Ncond`` counts
total electrons.
"""

import collections
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from tests.heavy_tests import heavy

_FIXTURES = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "uhfk_orientation")

#: What the subprocess runs: one real UHFk calculation, driven exactly as a
#: user drives it (``hwave.qlms.run`` on a parsed input), in the current
#: working directory -- which the runner has made a COPY of the fixture, so
#: no ``output/`` is ever written next to the committed inputs.
_RUN_SCRIPT = ("import json, sys; import hwave.qlms; "
               "hwave.qlms.run(input_dict=json.loads(sys.argv[1]))")

Case = collections.namedtuple(
    "Case", "label fixture norb files twobody spin_orbital fock moves why")


_BOND = {"CoulombIntra": "coulombintra.dat",
         "CoulombInter": "coulombinter.dat", "Hund": "hund.dat"}
_BOND_FILES = ("coulombinter.dat", "hund.dat")
_PAIRHOP = {"CoulombIntra": "coulombintra.dat", "PairHop": "pairhop.dat"}
_PAIRHOP_FILES = ("pairhop.dat",)

_HARTREE_ONLY = ("Hartree-only: the mean field reads the table as "
                 "sum_r J_ab(r), which the Hermitian closure makes "
                 "Hermitian and the orientation therefore leaves alone")
_ONE_ORBITAL = ("norb = 1 with a real transfer: F^rev's Hamiltonian is the "
                "complex conjugate of F's, and energy.dat's observables are "
                "conjugation-invariant")

#: The ten cases of the gate, keyed by label. ``moves`` says whether the
#: orientation is visible in ``energy_total`` at all -- see the module
#: docstring; ``why`` is the derivation for the ones where it is not.
CASES = {c.label: c for c in (
    Case("bond/norb2/fock=on", "norb2", 2, _BOND, _BOND_FILES,
         False, True, True, None),
    Case("bond/norb2/fock=off", "norb2", 2, _BOND, _BOND_FILES,
         False, False, False, _HARTREE_ONLY),
    Case("bond/norb2/fock=on/so", "norb2_so", 2, _BOND, _BOND_FILES,
         True, True, True, None),
    Case("bond/norb2/fock=off/so", "norb2_so", 2, _BOND, _BOND_FILES,
         True, False, False, _HARTREE_ONLY),
    Case("pairhop/norb2/fock=on", "norb2", 2, _PAIRHOP, _PAIRHOP_FILES,
         False, True, True, None),
    Case("pairhop/norb2/fock=off", "norb2", 2, _PAIRHOP, _PAIRHOP_FILES,
         False, False, True, None),
    Case("pairhop/norb2/fock=on/so", "norb2_so", 2, _PAIRHOP, _PAIRHOP_FILES,
         True, True, True, None),
    Case("pairhop/norb2/fock=off/so", "norb2_so", 2, _PAIRHOP, _PAIRHOP_FILES,
         True, False, True, None),
    Case("pairhop/norb1/fock=on", "norb1", 1, _PAIRHOP, _PAIRHOP_FILES,
         False, True, False, _ONE_ORBITAL),
    Case("pairhop/norb1/fock=on/so", "norb1_so", 1, _PAIRHOP, _PAIRHOP_FILES,
         True, True, False, _ONE_ORBITAL),
)}


def _params(case):
    """The input a user would write for ``case``, as the parsed dict
    ``hwave.qlms.run`` takes.

    ``flag_fock`` and ``enable_spin_orbital`` sit in the ``mode`` section,
    NOT in ``mode.param``: ``uhfk.UHFk.__init__`` passes ``info_mode`` (not
    ``info_mode["param"]``) to ``_init_mode``, so a ``flag_fock`` written
    one level down is silently ignored and every case would run with the
    default Fock term on.
    """
    interaction = {"path_to_input": ".", "Geometry": "geom.dat",
                   "Transfer": "transfer.dat"}
    interaction.update(case.files)
    return {
        "log": {"print_level": 0, "print_step": 1000},
        "mode": {"mode": "UHFk",
                 "enable_spin_orbital": case.spin_orbital,
                 "flag_fock": case.fock,
                 "param": {"T": 0.05,
                           "Ncond": 4 if case.norb == 2 else 2,
                           "IterationMax": 1000, "EPS": 12, "Mix": 0.5,
                           "RndSeed": 1,
                           "CellShape": [4, 1, 1], "SubShape": [1, 1, 1]}},
        "file": {"input": {"path_to_input": "", "interaction": interaction},
                 "output": {"path_to_output": "output",
                            "energy": "energy.dat", "green": "green.dat"}},
    }


def _read_energy(path):
    """``energy.dat`` as ``{key: float}`` -- the parse ``tests/test_uhf.py``
    already uses for the reference-data suite, minus its case-insensitive
    container (both sides here are written by the same writer, so the
    spellings match exactly and a difference in the KEY SET is a
    difference worth failing on)."""
    table = {}
    with open(path) as handle:
        lines = handle.read().splitlines()
    for line in lines:
        if not line.strip():
            continue
        key, _, value = line.partition("=")
        table[key.strip()] = float(value)
    return table


def _run(checkout, indir, case):
    """One UHFk run of ``case`` on the declaration in ``indir``, with
    ``checkout``'s source tree, in a COPY of ``indir`` under a temporary
    directory. Returns ``(energy, green)``.

    A subprocess, not an in-process call, for two reasons: the reference
    revision's ``hwave`` cannot be imported into the same interpreter as
    this branch's, and ``hwave.qlms.run`` configures process-wide logging
    and writes relative to the working directory.
    """
    work = tempfile.mkdtemp(prefix="hwave_uhfk_orient_")
    try:
        shutil.copytree(indir, work, dirs_exist_ok=True)
        env = dict(os.environ)
        env["PYTHONPATH"] = os.pathsep.join(
            (os.path.join(checkout, "src"), checkout))
        proc = subprocess.run(
            [sys.executable, "-B", "-c", _RUN_SCRIPT,
             json.dumps(_params(case))],
            cwd=work, env=env, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "UHFk failed for case {} on {} with {}:\n{}\n{}".format(
                    case.label, indir, checkout,
                    proc.stdout[-2000:], proc.stderr[-2000:]))
        energy = _read_energy(os.path.join(work, "output", "energy.dat"))
        with np.load(os.path.join(work, "output", "green.dat.npz"),
                     allow_pickle=True) as archive:
            green = {k: np.array(archive[k]) for k in archive.files}
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return energy, green


def _reversed_dir(case):
    """``F^rev`` of ``case``'s fixture in a fresh temporary directory (the
    caller removes it).

    The writer is ``tests/test_flex_hf_scf._reversed_interaction_dir``, the
    one this branch already uses for the FLEX half of the same contract --
    one implementation of "negate the displacement of every off-site row",
    not two. Its ``copy`` argument carries ``coulombintra.dat`` across
    verbatim: an on-site file has no row whose displacement could be
    negated, so passing it as a file to REVERSE would (rightly) be refused
    as vacuous.
    """
    from tests.test_flex_hf_scf import _reversed_interaction_dir
    return _reversed_interaction_dir(
        os.path.join(_FIXTURES, case.fixture), case.twobody,
        copy=("geom.dat", "transfer.dat", "coulombintra.dat"))


class _OrientationMixin:
    """Shared assertions, mixed into the two gates below (a plain mixin, so
    that no empty ``TestCase`` is discovered for it)."""

    def _assert_same_run(self, a, b, what):
        """Every member of two runs is equal: the key sets, every
        ``energy.dat`` number and every ``green.dat.npz`` array, at
        ``np.array_equal``.

        A member that is merely equal to 1e-12 relative is accepted and
        REPORTED (see the module docstring on why exactness is what is
        expected); anything beyond that fails and names the member.
        """
        (ea, ga), (eb, gb) = a, b
        self.assertEqual(sorted(ea), sorted(eb),
                         "{}: energy.dat key sets differ".format(what))
        self.assertEqual(sorted(ga), sorted(gb),
                         "{}: green.dat.npz member names differ".format(what))
        self.assertGreater(np.abs(ga["green"]).max(), 1e-6,
                           "{}: the Green function is zero, so comparing it "
                           "proves nothing".format(what))
        rounded = []
        for key in sorted(ea):
            if ea[key] == eb[key]:
                continue
            scale = max(abs(ea[key]), 1e-30)
            self.assertLessEqual(
                abs(ea[key] - eb[key]) / scale, 1e-12,
                "{}: energy.dat member {} differs by {:.3e} relative"
                .format(what, key, abs(ea[key] - eb[key]) / scale))
            rounded.append("{} ({:.3e} rel)".format(
                key, abs(ea[key] - eb[key]) / scale))
        for key in sorted(ga):
            if np.array_equal(ga[key], gb[key]):
                continue
            x, y = np.asarray(ga[key]), np.asarray(gb[key])
            self.assertIn(x.dtype.kind, "fc",
                          "{}: non-numeric member {} differs".format(what, key))
            scale = max(float(np.abs(x).max()), 1e-30)
            residual = float(np.abs(x - y).max()) / scale
            self.assertLessEqual(
                residual, 1e-12,
                "{}: green.dat.npz member {} differs by {:.3e} of its own "
                "size".format(what, key, residual))
            rounded.append("{} ({:.3e} rel)".format(key, residual))
        if rounded:
            print("\nNOTE {}: bit-identical except at round-off: {}"
                  .format(what, ", ".join(rounded)))

    def _relative_energy_gap(self, a, b):
        ea, eb = a[0], b[0]
        return abs(ea["Energy_Total"] - eb["Energy_Total"]) / max(
            abs(ea["Energy_Total"]), 1e-30)


class TestUHFkOrientation(_OrientationMixin, unittest.TestCase):
    """G-UHFk: this tree's UHFk on ``F`` equals the reference revision's
    UHFk on ``F^rev``, member for member -- and, wherever the orientation
    is visible at all, differs from the reference revision's UHFk on ``F``.

    Skipped as a whole when the reference checkout is not usable: the
    shared rule of every develop-comparison harness on this branch
    (``tests/test_flex_second_order_compat.develop_checkout``) is that the
    reference tree must be at the NAMED revision and clean, or the
    comparison is against a tree nobody can name.
    """

    def _check(self, *labels):
        from tests.test_flex_second_order_compat import develop_checkout
        here = os.getcwd()
        reference, why = develop_checkout()
        if reference is None:
            self.skipTest(why)
        for label in labels:
            case = CASES[label]
            source = os.path.join(_FIXTURES, case.fixture)
            reversed_dir = _reversed_dir(case)
            try:
                mine = _run(here, source, case)
                theirs = _run(reference, reversed_dir, case)
                control = _run(reference, source, case)
            finally:
                shutil.rmtree(reversed_dir, ignore_errors=True)
            self._assert_same_run(mine, theirs, "case {}".format(label))
            gap = self._relative_energy_gap(mine, control)
            if case.moves:
                self.assertGreater(
                    gap, 1e-6,
                    "case {}: the reference revision's UHFk on the "
                    "UNREVERSED declaration is within {:.3e} relative of "
                    "this tree's, so the identity above says nothing"
                    .format(label, gap))
            else:
                # the coinciding direction, asserted rather than assumed:
                # the reason is structural (see the module docstring), so a
                # change that makes the orientation visible here has to come
                # back through this test.
                self.assertLessEqual(
                    gap, 1e-12,
                    "case {}: the orientation moved a case that cannot see "
                    "it ({:.3e} relative). The derivation is: {}"
                    .format(label, gap, case.why))

    def test_bond_fock_on(self):
        """The inter-orbital ``CoulombInter`` + ``Hund`` bond with the Fock
        term on -- the one case of this module the FAST gate runs."""
        self._check("bond/norb2/fock=on")

    @heavy
    def test_bond_variants(self):
        """The same bond declaration with the Fock term off and in
        spin-orbital mode."""
        self._check("bond/norb2/fock=off", "bond/norb2/fock=on/so",
                    "bond/norb2/fock=off/so")

    @heavy
    def test_pairhop_complex(self):
        """The complex off-site ``PairHop`` amplitude, two orbitals and one,
        Fock term on and off, normal and spin-orbital mode."""
        self._check("pairhop/norb2/fock=on", "pairhop/norb2/fock=off",
                    "pairhop/norb2/fock=on/so", "pairhop/norb2/fock=off/so",
                    "pairhop/norb1/fock=on", "pairhop/norb1/fock=on/so")


class TestFixturesAreOrientationSensitive(_OrientationMixin,
                                          unittest.TestCase):
    """The same statement WITHOUT the reference tree: this tree's UHFk on
    ``F`` against this tree's UHFk on ``F^rev``.

    Its job is to keep the gate above honest when it skips. The identity
    it pins cannot be checked here -- that needs the reference revision --
    but the PREMISE can: that these fixtures are orientation-sensitive at
    all, so that the identity is a statement about the reading of the
    declaration and not about two runs that would agree anyway. One method
    per fixture directory, so a fixture that stops moving is named by the
    failure.
    """

    def _gap(self, label):
        case = CASES[label]
        here = os.getcwd()
        source = os.path.join(_FIXTURES, case.fixture)
        reversed_dir = _reversed_dir(case)
        try:
            return self._relative_energy_gap(_run(here, source, case),
                                             _run(here, reversed_dir, case))
        finally:
            shutil.rmtree(reversed_dir, ignore_errors=True)

    def _assert_moves(self, label):
        gap = self._gap(label)
        self.assertGreater(
            gap, 1e-6,
            "case {}: reversing every off-site displacement moves this "
            "tree's total energy by only {:.3e} relative, so the fixture "
            "cannot see the orientation".format(label, gap))

    def _assert_blind(self, label):
        gap = self._gap(label)
        self.assertLessEqual(
            gap, 1e-12,
            "case {}: reversing the declaration moved a case that cannot "
            "see the orientation ({:.3e} relative). The derivation is: {}"
            .format(label, gap, CASES[label].why))

    def test_norb2_bond(self):
        self._assert_moves("bond/norb2/fock=on")

    def test_norb2_pairhop(self):
        self._assert_moves("pairhop/norb2/fock=on")

    def test_norb2_so_bond(self):
        self._assert_moves("bond/norb2/fock=on/so")

    def test_hartree_only_is_blind(self):
        """``flag_fock = false``: the Hartree term reads the table only
        through ``sum_r J_ab(r)``, which the orientation cannot move."""
        self._assert_blind("bond/norb2/fock=off")

    def test_one_orbital_is_blind(self):
        """``norb = 1`` with a real transfer: ``F^rev``'s Hamiltonian is the
        complex conjugate of ``F``'s, and ``energy.dat`` is conjugation-
        invariant -- in normal and in spin-orbital mode alike."""
        self._assert_blind("pairhop/norb1/fock=on")
        self._assert_blind("pairhop/norb1/fock=on/so")


if __name__ == "__main__":
    unittest.main()
