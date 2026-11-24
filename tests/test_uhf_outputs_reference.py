"""The UHFr/UHFk output files of seven cases against the reference revision.

Every output file of the cases listed in :data:`CASES` (the ``print_check``
trace included) is byte-identical between the reference revision
(``tests.test_flex_second_order_compat.DEVELOP_COMMIT``) and this tree.
This module keeps it that way, with two declared exceptions that later
changes extend:

* :data:`ALLOWED_NEW_FILES`: files this tree may write and the reference
  revision does not;
* :data:`APPENDED_LINES`: files to which this tree may append complete
  lines, one per key, after the reference revision's last line.

Each case runs on both trees in subprocesses (the reference revision's
``hwave`` cannot be imported next to this one), each in a temporary copy of
its input directory, with ``[log] print_check = "check.dat"``.

The cases are the ones checked by hand to be byte-identical between the
two revisions. Not covered: spin-orbital UHFk (its ``Sy`` sign changed
with #350; :class:`TestSpinOrbitalEnergyKeys` pins its ``energy.dat``
key order instead), inputs that use ``initial_uhf``/``geometry_uhf``
(#357), and every case not listed.

CI provisions the reference revision and sets
``HWAVE_REQUIRE_DEVELOP_COMPARISON=1`` (``.github/workflows/
ci-python39.yml``), so a missing reference checkout fails there. Locally
the comparison skips with a reason when ``HWAVE_DEVELOP_CHECKOUT`` does not
name a clean checkout at the reference revision.
"""
import copy
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: Files this tree may write that the reference revision does not:
#: ``scf_status.dat``, the convergence status, written by default.
ALLOWED_NEW_FILES = {"scf_status.dat"}

#: ``{file name: keys}``: lines this tree appends, one per key, in order,
#: after the reference revision's last line of that file: the transverse
#: spin, written by both UHF solvers.
APPENDED_LINES = {"energy.dat": ("Sx", "Sy")}

#: ``{test name suffix: input directory relative to the repository root}``
CASES = {
    "uhfk_coulombintra": "tests/uhfk/CoulombIntra",
    "uhfk_2sz0": "tests/uhfk/2sz0",
    "uhfk_exchange": "tests/uhfk/Exchange",
    "uhfk_coulombinter": "tests/uhfk/CoulombInter",
    "tutorial_hubbard_uhfk": "docs/tutorial/Hubbard/UHFk",
    "tutorial_hubbard_uhfr": "docs/tutorial/Hubbard/UHFr",
    "uhfr_coulombintra": "tests/uhfr/coulombintra",
}

#: The subprocess: read ``input.toml`` in the working directory, add the
#: ``print_check`` trace and run it as ``hwave`` would.
_RUN_SCRIPT = (
    "import tomli, hwave.qlms\n"
    "with open('input.toml', 'rb') as f:\n"
    "    info = tomli.load(f)\n"
    "info.setdefault('log', {})['print_check'] = 'check.dat'\n"
    "hwave.qlms.run(input_dict=info)\n")


def _relative_files(root):
    found = set()
    for here, _dirs, files in os.walk(root):
        for name in files:
            found.add(os.path.relpath(os.path.join(here, name), root))
    return found


def compare_outputs(ref_dir, branch_dir, *, allowed_new_files, appended_lines):
    """The differences between two output directories, as a list of
    messages; empty when the outputs match.

    The file sets must be equal, except that ``branch_dir`` may also hold
    the files in ``allowed_new_files``. Common files are compared as lists
    of complete lines. For a file named in ``appended_lines`` (``{name:
    keys}``), the branch must be the reference followed by ``len(keys)``
    complete lines, the ``j``-th starting with ``keys[j] + " "``. Every
    other common file must be byte-equal.
    """
    problems = []
    ref_files = _relative_files(ref_dir)
    branch_files = _relative_files(branch_dir)
    allowed = set(allowed_new_files)
    for name in sorted(ref_files - branch_files):
        problems.append("{}: written by the reference only".format(name))
    for name in sorted(branch_files - ref_files - allowed):
        problems.append("{}: written by this tree only".format(name))
    for name in sorted(ref_files & branch_files):
        with open(os.path.join(ref_dir, name), "rb") as f:
            ref = f.read().splitlines(keepends=True)
        with open(os.path.join(branch_dir, name), "rb") as f:
            branch = f.read().splitlines(keepends=True)
        keys = appended_lines.get(name)
        if not keys:
            if ref != branch:
                problems.append("{}: differs from the reference".format(name))
            continue
        k = len(keys)
        if len(branch) < k or branch[:-k] != ref:
            problems.append("{}: is not the reference followed by {} lines"
                            .format(name, k))
            continue
        for j, key in enumerate(keys):
            line = branch[len(branch) - k + j]
            if not line.startswith(key.encode() + b" "):
                problems.append("{}: appended line {} does not start with {!r}: {!r}"
                                .format(name, j + 1, key, line))
            if not line.endswith(b"\n"):
                problems.append("{}: appended line {} is not a complete line: {!r}"
                                .format(name, j + 1, line))
    return problems


def _write_lines(d, name, lines):
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, name), "wb") as f:
        f.write(b"".join(lines))


class TestCompareOutputsHelper(unittest.TestCase):
    """:func:`compare_outputs` on hand-written directories."""

    REF = [b"Energy_Total = -1.0\n", b"NCond = 16.0\n",
           b"ChemicalPotential = 0.5\n"]
    KEYS = {"energy.dat": ("Sx", "Sy")}

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.ref = os.path.join(self.d, "ref")
        self.branch = os.path.join(self.d, "branch")
        _write_lines(self.ref, "energy.dat", self.REF)

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _compare(self, branch_lines, allowed=(), appended=None):
        _write_lines(self.branch, "energy.dat", branch_lines)
        return compare_outputs(self.ref, self.branch, allowed_new_files=set(allowed),
                               appended_lines=self.KEYS if appended is None else appended)

    def test_identical_accepted(self):
        self.assertEqual(self._compare(self.REF, appended={}), [])

    def test_trailing_keys_accepted(self):
        self.assertEqual(self._compare(self.REF + [b"Sx = 0.0\n", b"Sy = 0.0\n"]), [])
        self.assertEqual(
            self._compare(self.REF + [b"Sx      = 0.0\n", b"Sy      = 0.0\n"]), [])

    def test_keys_in_middle_rejected(self):
        branch = self.REF[:2] + [b"Sx = 0.0\n", self.REF[2], b"Sy = 0.0\n"]
        problems = self._compare(branch)
        self.assertEqual(len(problems), 1, problems)
        self.assertIn("energy.dat", problems[0])

    def test_keys_replacing_line_rejected(self):
        problems = self._compare(self.REF[:2] + [b"Sx = 0.0\n", b"Sy = 0.0\n"])
        self.assertEqual(len(problems), 1, problems)
        problems = self._compare(self.REF[:2] + [b"Sx = 0.0\n"], appended={"energy.dat": ("Sx",)})
        self.assertEqual(len(problems), 1, problems)

    def test_wrong_key_rejected(self):
        self.assertEqual(len(self._compare(self.REF + [b"Sy = 0.0\n", b"Sx = 0.0\n"])), 2)
        # the key must be followed by a space: "Sxx" is not "Sx"
        self.assertEqual(len(self._compare(self.REF + [b"Sxx = 0.0\n", b"Sy = 0.0\n"])), 1)

    def test_partial_line_rejected(self):
        problems = self._compare(self.REF + [b"Sx = 0.0\n", b"Sy = 0.0"])
        self.assertEqual(len(problems), 1, problems)
        self.assertIn("complete", problems[0])

    def test_byte_change_rejected(self):
        problems = self._compare([self.REF[0], b"NCond = 16.00\n", self.REF[2]], appended={})
        self.assertEqual(len(problems), 1, problems)

    def test_file_sets(self):
        _write_lines(self.ref, "only_ref.dat", [b"x\n"])
        _write_lines(self.branch, "scf_status.dat", [b"x\n"])
        problems = self._compare(self.REF, appended={})
        self.assertEqual(len(problems), 2, problems)
        problems = self._compare(self.REF, allowed={"scf_status.dat"}, appended={})
        self.assertEqual(len(problems), 1, problems)
        self.assertIn("only_ref.dat", problems[0])


def _run_case(tree, source, work):
    """Run the case in ``source`` with ``tree``'s source in a copy of it
    under ``work``; returns the output directory."""
    from tests.test_flex_second_order_compat import (own_tree_subprocess_env,
                                                     reference_subprocess_env)
    # a committed "output" (the tutorials ship one) must not reach the copy
    shutil.copytree(source, work, ignore=shutil.ignore_patterns("output"))
    env = (own_tree_subprocess_env(tree) if tree == _REPO
           else reference_subprocess_env(tree))
    proc = subprocess.run([sys.executable, "-B", "-c", _RUN_SCRIPT],
                          cwd=work, env=env, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError("{} failed with the tree {}:\n{}\n{}".format(
            source, tree, proc.stdout[-2000:], proc.stderr[-2000:]))
    return os.path.join(work, "output")


class TestOutputsMatchReference(unittest.TestCase):
    """This tree's output files of each case against the reference
    revision's, by :func:`compare_outputs`."""

    def _check(self, name):
        from tests.test_flex_second_order_compat import develop_checkout
        reference, why = develop_checkout()
        if reference is None:
            self.skipTest(why)
        source = os.path.join(_REPO, CASES[name])
        tmp = tempfile.mkdtemp(prefix="hwave_uhf_ref_")
        try:
            ref_out = _run_case(reference, source, os.path.join(tmp, "reference"))
            branch_out = _run_case(_REPO, source, os.path.join(tmp, "branch"))
            problems = compare_outputs(ref_out, branch_out,
                                       allowed_new_files=ALLOWED_NEW_FILES,
                                       appended_lines=APPENDED_LINES)
            # the comparison is not vacuous: both traces were written
            self.assertLessEqual({"check.dat", "energy.dat"}, _relative_files(ref_out))
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
        self.assertEqual(problems, [], "{}: {}".format(CASES[name], problems))

    def test_uhfk_coulombintra(self):
        self._check("uhfk_coulombintra")

    def test_uhfk_2sz0(self):
        self._check("uhfk_2sz0")

    def test_uhfk_exchange(self):
        self._check("uhfk_exchange")

    def test_uhfk_coulombinter(self):
        self._check("uhfk_coulombinter")

    def test_tutorial_hubbard_uhfk(self):
        self._check("tutorial_hubbard_uhfk")

    def test_tutorial_hubbard_uhfr(self):
        self._check("tutorial_hubbard_uhfr")

    def test_uhfr_coulombintra(self):
        self._check("uhfr_coulombintra")


class TestSpinOrbitalEnergyKeys(unittest.TestCase):
    """Spin-orbital UHFk is not in :data:`CASES` (its ``Sy`` sign changed
    with #350), so its ``energy.dat`` is pinned here: one site in the field
    ``H = -sigma_x`` (Transfer entry ``(a, b)`` is the coefficient of
    ``c^dag_a c_b``, index ``2*orb + spin + 1``)."""

    KEYS = ["Energy_Total", "Energy_Band", "FreeEnergy_Total", "FreeEnergy_Band",
            "NCond", "Sz", "ChemicalPotential", "Sx", "Sy"]

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self._level = logging.getLogger().level
        logging.getLogger().setLevel(logging.ERROR)

    def tearDown(self):
        logging.getLogger().setLevel(self._level)
        shutil.rmtree(self.d, ignore_errors=True)

    def test_key_order_and_spin(self):
        from hwave import qlms
        with open(os.path.join(self.d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n2\n0 0 0\n0 0 0\n")
        with open(os.path.join(self.d, "transfer.dat"), "w") as f:
            f.write("field\n2\n1\n 1\n0 0 0 1 2 -1.0 0.0\n0 0 0 2 1 -1.0 0.0\n")
        inp = {
            "log": {"print_level": 0},
            "mode": {"mode": "UHFk", "enable_spin_orbital": True, "param": {
                "CellShape": [1, 1, 1], "SubShape": [1, 1, 1], "T": 0.0,
                "Ncond": 1, "IterationMax": 50}},
            "file": {"input": {"path_to_input": self.d, "interaction": {
                "path_to_input": self.d, "Geometry": "geom.dat",
                "Transfer": "transfer.dat"}},
                "output": {"path_to_output": os.path.join(self.d, "out"),
                           "energy": "energy.dat", "green": "green"}},
        }
        qlms.run(input_dict=copy.deepcopy(inp))
        values = []
        with open(os.path.join(self.d, "out", "energy.dat")) as f:
            for line in f:
                key, _, value = line.partition("=")
                values.append((key.strip(), float(value)))
        self.assertEqual([k for k, _ in values], self.KEYS)
        table = dict(values)
        self.assertAlmostEqual(table["Sx"], 0.5, places=10)
        self.assertAlmostEqual(table["Sy"], 0.0, places=10)
        self.assertAlmostEqual(table["Sz"], 0.0, places=10)


if __name__ == "__main__":
    unittest.main()
