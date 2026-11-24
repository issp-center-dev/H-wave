"""The standalone Hartree-Fock-renormalised FLEX loop (spec 2026-09-06
sections 1, 2, 4.2-4.3, 5.6; gates G0-off, G5, D7): gate-off inputs are
identical to develop (exactly with numpy >= 2, to round-off with numpy < 2),
the single-band Hubbard HF term is a chemical-potential shift, the
mean-field seed follows the initial-value trajectory, the W-disabled fixed
point equals paramagnetic UHFk, the reset set, IterationMax = 0 and the provenance block.

TestG0Off's two comparisons need a SECOND source tree at
tests/test_flex_second_order_compat.DEVELOP_COMMIT. CI provisions that
reference revision (a detached worktree under RUNNER_TEMP; see
.github/workflows/ci-python39.yml) and sets
HWAVE_REQUIRE_DEVELOP_COMPARISON=1, so they RUN there and a missing or
unusable reference fails rather than skips. Locally they skip with a
reason when no such checkout is at hand."""
import decimal
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_IN1 = "tests/rpa/input"


def _build(param_extra=None, path=_IN2, interactions=None, mixing="linear"):
    import hwave.solver.flex as flex_mod
    idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    idict.update({"CoulombInter": "coulombinter.dat"} if interactions is None else interactions)
    r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
           "Nmat": 32, "IterationMax": 4, "Mix": 0.5, "EPS": 8, "mixing_scheme": mixing}
    par.update(param_extra or {})
    info = {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"}
    return flex_mod.FLEX(r.get_param("ham"), {}, info), r


_ISSUE282_DIAGNOSTICS = re.compile(
    r"^hwave\.solver\.flex (?:INFO|DEBUG) {3}Stoner factor max_q eig\[U chi0\(q, i nu = 0\)\]: .*\n"
    r"|^hwave\.solver\.flex DEBUG Stoner factor not computed: .*\n"
    r"|^hwave\.solver\.flex WARNING FLEX iteration \d+ looks divergent: .*\n",
    re.MULTILINE)


def _strip_issue282_diagnostics(log_text):
    """Remove the additive diagnostic lines of issue #282 (see the
    develop-comparison test below) from a log in this module's
    "%(name)s %(levelname)s %(message)s" format."""
    return _ISSUE282_DIAGNOSTICS.sub("", log_text)


#: a floating-point literal in a log line: digits with a decimal point
#: and/or an exponent (a plain integer is not matched), optionally signed.
#: It must stand on its own: not preceded by a word character, a dot, a
#: sign, a hyphen or a path separator, and not followed by a word
#: character, a path separator, or a dot, sign or hyphen that continues
#: into a word or number -- so decimal pieces of file names, path
#: components and version strings (file-1.25.dat, run/1.25/output, 1.25-x,
#: v1.26.4, 1.26.4) are left alone, while a sentence-ending period or
#: comma after a value is allowed.
_FLOAT_LITERAL = re.compile(
    r"(?<![\w.+\-/\\])[-+]?(?:\d+\.\d+|\.\d+|\d+(?=[eE][-+]?\d))(?:[eE][-+]?\d+)?"
    r"(?![\w/\\]|[.+\-]\w)")

#: relative tolerance of one printed float of the develop log comparison
#: with numpy < 2 (the same allowance as the archives' ROUNDOFF_RTOL, here
#: relative to the larger magnitude of the two printed values)
LOG_FLOAT_RTOL = 1e-12


#: _FLOAT_LITERAL as the single capture group re.split needs to keep the
#: matched literals in its result
_FLOAT_TOKEN = re.compile("(" + _FLOAT_LITERAL.pattern + ")")


def _log_line_tokens(line):
    """Split a log line into alternating text and float-literal tokens:
    ``[text, float, text, ..., text]`` (re.split with one capture group)."""
    return _FLOAT_TOKEN.split(line)


def _log_lines_differ(la, lb, rtol=LOG_FLOAT_RTOL):
    """Why line ``la`` does not match line ``lb`` token by token, or None.

    Non-numeric text must be identical (so file names, path components and
    version strings, which _FLOAT_LITERAL never matches, compare exactly),
    the float literals must sit at the same places, and each pair must agree
    to ``rtol`` relative to the larger magnitude of the two, compared as
    exact decimal values (so out-of-range literals are compared too)."""
    ta, tb = _log_line_tokens(la), _log_line_tokens(lb)
    if len(ta) != len(tb):
        return "different number of float literals"
    for i, (x, y) in enumerate(zip(ta, tb)):
        if i % 2 == 0:
            if x != y:
                return "text {!r} != {!r}".format(x, y)
        else:
            # compared as exact decimals, not floats: a literal outside the
            # double range would overflow to inf or underflow to 0 and make
            # different values (1e309 vs 1e308, 1e-400 vs 2e-400) compare
            # equal
            dx, dy = decimal.Decimal(x), decimal.Decimal(y)
            if abs(dx - dy) > decimal.Decimal(repr(rtol)) * max(abs(dx), abs(dy)):
                return "float {} != {} beyond rtol {}".format(x, y, rtol)
    return None


def _assert_logs_match(test, la, lb, below_2):
    """Assert that log text ``la`` matches the reference's ``lb``.

    With numpy >= 2 (``below_2`` false) the two must be identical strings.
    With numpy < 2 the last-bit drift of the C-ordered Matsubara transforms
    (#309) can change a printed value (32.0 vs 31.999999999999996), and any
    fixed-digit rounding of both sides has boundaries such a drift can
    straddle; so the logs are compared line by line and token by token
    (:func:`_log_lines_differ`). A failure names the first differing line."""
    if not below_2:
        test.assertEqual(la, lb)
        return
    lines_a, lines_b = la.split("\n"), lb.split("\n")
    for n, (x, y) in enumerate(zip(lines_a, lines_b), 1):
        why = _log_lines_differ(x, y)
        if why is not None:
            test.fail("log line {} differs ({}):\n  here:      {}\n  reference: {}"
                      .format(n, why, x, y))
    if len(lines_a) != len(lines_b):
        test.fail("logs differ in length: {} vs {} lines; first extra line: {!r}".format(
            len(lines_a), len(lines_b),
            (lines_a if len(lines_a) > len(lines_b) else lines_b)[min(len(lines_a),
                                                                      len(lines_b))]))


_G0_SCRIPT = r'''
import os, sys, json, logging, numpy as np
sys.path.insert(0, sys.argv[1] + "/src"); sys.path.insert(0, sys.argv[1])
import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.flex as flex_mod
out = sys.argv[2]
logging.basicConfig(level=logging.DEBUG, filename=os.path.join(out, "log.txt"), filemode="w",
                    format="%(name)s %(levelname)s %(message)s")
p = "tests/rpa/input_2orb"
idict = {"path_to_input": p, "Geometry": "geom.dat", "Transfer": "transfer.dat",
         "CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"}
r = read_input_k.QLMSkInput({"path_to_input": p, "interaction": idict})
par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 32,
       "IterationMax": 3, "Mix": 0.5, "EPS": 8, "mixing_scheme": "anderson",
       "flex_second_order": "local",
       # both orbitals on one site: the behaviour of the reference revision,
       # which ignores the key, stated so that no #234 warning enters the log
       "orbital_site": [0, 0],
       # the reference revision predates the Hermitian projection of the
       # self-energy (issue #281), which changes a symmetric run at round-off
       # only (pinned in tests/test_flex_sigma_hermitian.py); switched off so
       # that this comparison with the reference stays about the code it names
       "flex_sigma_hermitian_projection": False}
s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par, "enable_spin_orbital": False, "calc_scheme": "general"})
gi = r.get_param("green")
s.solve(gi, out)
s.save_results({"path_to_output": out, "sigma": "sigma.npz", "green": "green.npz", "chi0q": "chi0q.npz", "chiq": "chiq.npz"}, gi)
'''


_HF_FIRST_MAP_SCRIPT = r'''
import os, sys, numpy as np
sys.path.insert(0, sys.argv[1] + "/src"); sys.path.insert(0, sys.argv[1])
import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.flex as flex_mod
_root, inp, out = sys.argv[1:4]
idict = {"path_to_input": inp, "Geometry": "geom.dat", "Transfer": "transfer.dat",
         "CoulombInter": "coulombinter.dat"}
r = read_input_k.QLMSkInput({"path_to_input": inp, "interaction": idict})
# IterationMax = 1 with Mix = 1.0 and no seed: the state after the run IS the
# first Hartree-Fock map, taken on the bare (Sigma = 0) density, unmixed
par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 32,
       "IterationMax": 1, "Mix": 1.0, "EPS": 8, "mixing_scheme": "linear",
       "flex_hartree_fock": True, "flex_second_order": "local"}
s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                           "enable_spin_orbital": False, "calc_scheme": "general"})
gi = r.get_param("green")
s.solve(gi, out)
np.savez(os.path.join(out, "static.npz"), sigma_static=np.asarray(gi["sigma_static"]))
'''


def _reversed_interaction_dir(src_dir, files, copy=("geom.dat", "transfer.dat")):
    """A fresh input directory holding ``F^rev``: every listed two-body file
    rewritten with the displacement of each OFF-SITE row negated, and every
    file in ``copy`` (the geometry and the transfer by default) copied
    verbatim. Returns the directory (the caller removes it).

    ``copy`` exists because a declaration is not always geometry + transfer
    + the reversed files: tests/test_uhfk_orientation.py runs a
    ``CoulombIntra`` alongside the off-site types, and an ON-SITE file has
    no row whose displacement could be negated -- passing it in ``files``
    would (rightly) be refused as vacuous below, so it travels in ``copy``.

    ``F^rev`` is what a user of the reference revision's VERTEX paths
    (RPA, the FLEX vertex) has to write to get the Hamiltonian this branch
    reads from ``F``: a row ``(r, a, b, v)`` means ``v n_{j,a} n_{j+r,b}``,
    the bond of a Transfer row ``(r, a, b)`` (issue #232), while those paths
    used to put it on ``v n_{j,b} n_{j+r,a}``; the two readings differ
    exactly by ``r -> -r``. Orbital indices are NOT touched; ``r = 0`` rows
    are copied as they are.

    The wannier90-like header is four lines (title, norb, n_rvec,
    degeneracies) and the rows follow; only the OFF-SITE rows are rewritten,
    in place and column by column, so that a hand diff of F against F^rev
    shows exactly the rows whose displacement was negated. Every
    declared degeneracy must be 1, which is what H-wave's own files write
    and what makes the row order irrelevant -- a file with a real
    Wannier90 degeneracy block would need its header permuted too, so this
    helper refuses one instead of writing a wrong file."""
    dst = tempfile.mkdtemp(prefix="hwave_rev_decl_")
    for f in copy:
        shutil.copy(os.path.join(src_dir, f), dst)
    for f in files:
        with open(os.path.join(src_dir, f)) as fr:
            lines = fr.read().splitlines()
        header, body = lines[:4], lines[4:]
        if any(int(x) != 1 for x in header[3].split()):
            raise ValueError("_reversed_interaction_dir: {} declares a non-unit "
                             "degeneracy block; reversing it would need the header "
                             "permuted too".format(f))
        out, reversed_rows = list(header), 0
        for ln in body:
            parts = ln.split()
            if len(parts) < 7 or not any(int(x) for x in parts[:3]):
                out.append(ln)             # unchanged rows go through verbatim
                continue
            # rewrite the three displacement columns IN PLACE, each into the
            # same FIELD (its leading whitespace included, so a longer value
            # such as -1 for 1 eats a space instead of shifting the line), so
            # that a hand diff of F against F^rev shows the negated rows and
            # nothing else
            rewritten, cursor = "", 0
            for x in parts[:3]:
                end = ln.index(x, cursor) + len(x)
                rewritten += str(-int(x)).rjust(end - cursor)
                cursor = end
            out.append(rewritten + ln[cursor:])
            reversed_rows += 1
        if reversed_rows == 0:
            raise ValueError("_reversed_interaction_dir: {} has no off-site row, so "
                             "F^rev == F and the comparison would be vacuous".format(f))
        with open(os.path.join(dst, f), "w") as fw:
            fw.write("\n".join(out) + "\n")
    return dst


class TestAssertLogsMatch(unittest.TestCase):
    """The log comparator of TestG0Off, on both sides of the numpy gate."""

    def _ok(self, la, lb):
        _assert_logs_match(self, la, lb, below_2=True)

    def _bad(self, la, lb):
        with self.assertRaises(AssertionError):
            _assert_logs_match(self, la, lb, below_2=True)

    def test_last_bit_drift_passes(self):
        self._ok("a\nNcond = 31.999999999999996\nb", "a\nNcond = 32.0\nb")
        self._ok("res 1.0000000000000002e-05.", "res 1e-05.")
        self._ok("x = -0.30000000000000004, y [0.1, 2.0]", "x = -0.3, y [0.1, 2.0]")

    def test_exact_on_numpy_2(self):
        with self.assertRaises(AssertionError):
            _assert_logs_match(self, "Ncond = 31.999999999999996", "Ncond = 32.0",
                               below_2=False)
        _assert_logs_match(self, "Ncond = 32.0", "Ncond = 32.0", below_2=False)

    def test_real_number_change_fails(self):
        self._bad("mu = 1.000000000004", "mu = 1.0")
        self._bad("x -1.5", "x 1.5")
        self._bad("e 1.5e-05", "e 1.5e-06")
        self._bad("mu = 0.5", "mu = 5")

    def test_out_of_double_range_literals(self):
        self._bad("x 1e309", "x 1e308")
        self._bad("x 1e309", "x -1e309")
        self._bad("x 1e-400", "x 2e-400")
        self._ok("x 1e309", "x 1.0e309")
        self._ok("x 1.0e-3", "x 0.001")

    def test_text_difference_fails(self):
        self._bad("hwave INFO converged 1.0", "hwave INFO diverged 1.0")
        self._bad("iteration 3: 1.0", "iteration 4: 1.0")         # integers are text

    def test_path_and_version_difference_fails(self):
        self._bad("read file-1.23456789012345.dat", "read file-1.23456789012346.dat")
        self._bad("read run/1.000000000000001/output", "read run/1.0/output")
        self._bad("numpy v1.26.4", "numpy v1.26.5")
        self._bad("numpy 1.26.4", "numpy 1.26.40")
        self._bad("saved chi_1.5.npz", "saved chi_1.50.npz")

    def test_extra_number_or_line_fails(self):
        self._bad("values 1.0 2.0", "values 1.0")
        self._bad("values 1.0", "values 1.0 x")
        self._bad("a\nb", "a\nb\nc")

    def test_failure_names_the_first_differing_line(self):
        with self.assertRaises(AssertionError) as cm:
            self._ok("a 1.0\nb 2.0\nc 3.0", "a 1.0\nb 2.5\nc 3.5")
        self.assertIn("log line 2", str(cm.exception))

    def test_float_literal_boundaries(self):
        """What _FLOAT_LITERAL takes as a number (the rest is text)."""
        floats = lambda t: _log_line_tokens(t)[1::2]
        self.assertEqual(floats("a 1e-05 b 2.50E+03, c -0.25."), ["1e-05", "2.50E+03", "-0.25"])
        self.assertEqual(floats("[0.1, 2.0]"), ["0.1", "2.0"])
        for text in ("iteration 3 of 10: nvol = 16, norb 2", "shape (4, 4, 1) step 7.",
                     "index [3][-1] k_0 1e x2.5", "numpy v1.26.4", "numpy 1.26.4 and 2.4.6",
                     "read file-1.23456789012345.dat", "read /tmp/run/1.2345/output",
                     "run\\1.5\\x", "saved chi_1.5.npz", "range 1.0000000000000002-x"):
            self.assertEqual(floats(text), [], text)


class TestG0Off(unittest.TestCase):
    """Gate-off inputs: every archive member and the complete log stream
    equal those of the reference revision (the main checkout, pinned by
    tests/test_flex_second_order_compat.DEVELOP_COMMIT): exactly with
    numpy >= 2; with numpy < 2 floating-point values to round-off (archive
    members via _assert_same_to_roundoff, printed floats in the log via
    _assert_logs_match).

    The run below pins flex_second_order = "local" explicitly, on BOTH
    sides (the default of both revisions; spec 2026-09-08 D2). It used to
    pin the "takimoto" kernel, which this tree refuses since issue #253;
    on this on-site declaration the "local" kernel of the reference
    revision, given issue #251's ring orientation by
    tests/test_flex_second_order_compat.reference_script, is the same
    computation as this tree's.

    The declaration it runs is ON-SITE only (onsite_inter.dat +
    coulombintra.dat): an on-site declaration has no row whose placement
    could be read at all, so this guard stays a statement about the Phase B
    gates and not about the placement of off-site rows, which issue #232
    corrected in the vertex builders (the reference revision still puts an
    off-site row on the opposite bond there).

    The gate-ON half of the off-site story is the IDENTITY below: the
    Hartree-Fock map of the same declaration F on both revisions.

    The npz comparison allows this side's archives to carry additional
    members beyond the reference revision's -- every member the reference
    DOES have must still match: exactly with numpy >= 2; with numpy < 2
    integers exactly and floating-point values to round-off (the memory
    layout of the Matsubara transforms changed in #309, which moves the
    last bits of numpy < 2's FFT and BLAS results).
    """

    def test_matches_develop(self):
        from tests.test_flex_second_order_compat import (develop_checkout,
                                                           own_tree_subprocess_env,
                                                           reference_script,
                                                           reference_subprocess_env, moved_stamp,
                                                           _assert_same_to_roundoff,
                                                           NUMPY_BELOW_2)
        here = os.getcwd()
        # one shared rule with the other develop-comparison harnesses: the
        # reference tree must be at the NAMED revision and clean, or the
        # comparison is against something nobody can name
        develop, why = develop_checkout()
        if develop is None:
            self.skipTest(why)
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            for root, out in ((here, a), (develop, b)):
                # _G0_SCRIPT puts root/src on sys.path itself, so PYTHONPATH
                # from the helper is harmless; only the reference invocation
                # needs pytest-cov's subprocess hooks stripped -- this
                # tree's run stays measured
                env = (own_tree_subprocess_env(root) if root == here
                       else reference_subprocess_env(root))
                # the reference run gets issue #251's orientation change only
                script = _G0_SCRIPT if root == here else reference_script(_G0_SCRIPT)
                subprocess.run([sys.executable, "-B", "-c", script, root, out],
                               check=True, cwd=here, env=env)
            for name in ("sigma.npz", "green.npz", "chi0q.npz", "chiq.npz"):
                da, db = np.load(os.path.join(a, name), allow_pickle=True), np.load(os.path.join(b, name), allow_pickle=True)
                # the same treatment as test_flex_second_order_compat.py's
                # _members() comparison: every member the reference revision
                # writes must be present here and equal, while this side may
                # carry members the reference does not have (a later
                # provenance stamp is not a compatibility break). The
                # reference revision does write flex_second_order /
                # flex_second_order_schema -- it is the merge that added
                # them -- so those are compared, not tolerated.
                self.assertEqual(set(db.files) - set(da.files), set(), name)
                self.assertIn("flex_second_order", db.files, name)
                for k in db.files:
                    if moved_stamp(k, db[k], da[k]):        # b = reference, a = here
                        continue
                    _assert_same_to_roundoff(self, da[k], db[k], (name, k))
            with open(os.path.join(a, "log.txt")) as fh:
                la = fh.read().replace(here, "<root>").replace(a, "<out>")
            # Issue #282 added documented, purely additive diagnostic lines
            # that the reference revision cannot write: the per-iteration
            # Stoner factor line (INFO on the first and last iteration,
            # DEBUG otherwise; DEBUG "Stoner factor not computed" if its
            # evaluation fails) and the one-off "looks divergent" WARNING of
            # the iteration monitor. Only those lines are removed, from this
            # side only; every other line must still match byte for byte
            # (with numpy < 2 up to round-off in the printed floats, see below).
            la = _strip_issue282_diagnostics(la)
            with open(os.path.join(b, "log.txt")) as fh:
                lb = fh.read().replace(develop, "<root>").replace(b, "<out>")
            # the whole stream, nothing filtered: the reference revision is
            # the merge that introduced flex_second_order, so it logs the
            # same "flex_second_order = local" INFO line this side does
            # (an earlier version of this guard stripped that line, which
            # would now DELETE a line from one side only)
            self.assertIn("flex_second_order = local", la)
            # exact with numpy >= 2; with numpy < 2 printed floats to
            # round-off, token by token (the last-bit drift of the archives
            # can flip a printed repr, e.g. 32.0 vs 31.999999999999996)
            _assert_logs_match(self, la, lb, NUMPY_BELOW_2)

    def test_first_map_static_equals_develop_on_the_same_declaration(self):
        """Spec 2.4 end to end, as corrected by issue #232: with an OFF-SITE
        declaration F, this branch's first Hartree-Fock map equals the
        reference revision's first map on the SAME declaration F, and differs
        from the reference revision's map on the reversed declaration F^rev.

        An off-site row (r, a, b, v) is read as v n_{j,a} n_{j+r,b}, the bond
        of a Transfer row (r, a, b) (Wannier90 / RESPACK); the reference
        revision predates the reversal step issue #193 added and issue #232
        removed, so it reads the row the same way. Iteration 1 with Mix = 1.0
        and no seed is compared because both sides then map the SAME (bare,
        Sigma = 0) density, so any difference is the interaction reading and
        nothing else; the in-process component version, on the map alone, is
        tests/test_flex_second_order_compat.py::TestCompatibility::
        test_hf_map_reads_the_declaration_as_develop_did.

        np.array_equal, not a tolerance: the two tables are built from the
        same values by the same code. A tolerance here would hide a real
        difference in the map.
        """
        from tests.test_flex_second_order_compat import (develop_checkout,
                                                           own_tree_subprocess_env,
                                                           reference_subprocess_env)
        here = os.getcwd()
        develop, why = develop_checkout()
        if develop is None:
            self.skipTest(why)
        rev = _reversed_interaction_dir(_IN2, ("coulombinter.dat",))
        try:
            # anti-vacuity: the reversed declaration really is a different file
            with open(os.path.join(rev, "coulombinter.dat")) as fh:
                reversed_text = fh.read()
            with open(os.path.join(_IN2, "coulombinter.dat")) as fh:
                original_text = fh.read()
            self.assertNotEqual(reversed_text, original_text)
            with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b, \
                    tempfile.TemporaryDirectory() as c:
                for root, inp, out in ((here, os.path.abspath(_IN2), a),
                                       (develop, os.path.abspath(_IN2), b),
                                       (develop, rev, c)):
                    # same rule as above: only the reference invocations lose
                    # pytest-cov's subprocess hooks
                    env = (own_tree_subprocess_env(root) if root == here
                       else reference_subprocess_env(root))
                    subprocess.run([sys.executable, "-B", "-c", _HF_FIRST_MAP_SCRIPT,
                                    root, inp, out], check=True, cwd=here, env=env)
                mine = np.load(os.path.join(a, "static.npz"))["sigma_static"]
                theirs = np.load(os.path.join(b, "static.npz"))["sigma_static"]
                on_rev = np.load(os.path.join(c, "static.npz"))["sigma_static"]
            scale = np.abs(mine).max()
            self.assertGreater(scale, 1e-6)                          # anti-vacuity
            self.assertEqual(mine.shape, theirs.shape)
            self.assertTrue(np.array_equal(mine, theirs),
                            "first Hartree-Fock map differs from the reference "
                            "revision's on F by {:.3e} of its own size"
                            .format(np.abs(mine - theirs).max() / scale))
            # anti-vacuity: the fixture sees the placement -- the reference
            # revision on the REVERSED declaration lands somewhere else
            gap = np.abs(mine - on_rev).max() / scale
            self.assertGreater(gap, 1e-3,
                               "reversing the declaration does not move this fixture's "
                               "first map ({:.3e} relative); the identity above "
                               "would be vacuous".format(gap))
        finally:
            shutil.rmtree(rev, ignore_errors=True)


class TestStandaloneHF(unittest.TestCase):

    def test_single_band_hubbard_hf_is_a_mu_shift(self):
        """norb=1, on-site U only: Sigma_HF is a constant absorbed by mu, so
        the dressed Green function equals the HF-off one."""
        s_on, r1 = _build({"flex_hartree_fock": True}, path=_IN1, interactions={"CoulombIntra": "coulombintra.dat"})
        s_off, r2 = _build({}, path=_IN1, interactions={"CoulombIntra": "coulombintra.dat"})
        g_on, g_off = r1.get_param("green"), r2.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s_on.solve(g_on, "tests/flex/output"); s_off.solve(g_off, "tests/flex/output")
        st = np.asarray(g_on["sigma_static"])
        self.assertLess(np.abs(st - st[0, 0, 0, 0, 0]).max(), 1e-12)     # k-independent constant
        np.testing.assert_allclose(np.asarray(g_on["green"]), np.asarray(g_off["green"]), rtol=0, atol=1e-10)
        self.assertAlmostEqual(g_on["physics"]["mu"] - g_off["physics"]["mu"], st[0, 0, 0, 0, 0].real, places=9)
        self.assertAlmostEqual(g_on["physics"]["NCond"], g_off["physics"]["NCond"], places=8)

    def test_d7_trajectory_delta_h_is_never_readded(self):
        """A trans_mod seed Delta H with the HF map forced to zero and
        linear mixing m: the static state is (1-m)^k Delta H after k maps."""
        import hwave.solver.flex_hf as flex_hf
        s, r = _build({"flex_hartree_fock": True, "IterationMax": 2, "Mix": 0.25})
        nvol, norb = s.lattice.nvol, s.norb
        shift = np.diag([0.3, -0.2])
        H = np.zeros((nvol, 2, norb, 2, norb), complex)
        base = np.asarray(s.ham_info.ham_trans_q).reshape(nvol, norb, norb)
        H[:, 0, :, 0, :] = base + shift; H[:, 1, :, 1, :] = base + shift
        gi = r.get_param("green"); gi["trans_mod"] = H.reshape(nvol, 2 * norb, 2 * norb)
        real = flex_hf.hf_map
        flex_hf.hf_map = lambda rho_r, tables, shape, norb_: np.zeros((rho_r.shape[0], norb_, norb_), complex)
        seen = []
        s._iteration_hook = lambda d: seen.append(d)
        try:
            s.solve(gi, "tests/flex/output")
        finally:
            flex_hf.hf_map = real
        self.assertEqual(len(seen), 2)
        self.assertEqual(np.abs(seen[0]["static_new"]).max(), 0.0)              # the map is zero
        np.testing.assert_allclose(np.asarray(gi["sigma_static"])[0, 0],
                                   0.75 ** 2 * np.broadcast_to(shift, (nvol, norb, norb)), rtol=0, atol=1e-12)
        self.assertIn("trans_mod", gi)                                          # input preserved

    def test_reset_set_iteration_max_zero_and_provenance(self):
        s, r = _build({"flex_hartree_fock": True})
        gi = r.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s.solve(gi, "tests/flex/output")
        self.assertEqual(s.scf_iterations, 4)
        for k in ("sigma", "sigma_static", "sigma_fluct", "green", "physics", "chi0q", "chiq_s", "chiq_c"):
            self.assertIn(k, gi)
        np.testing.assert_allclose(gi["sigma"], gi["sigma_static"] + gi["sigma_fluct"], atol=1e-14)
        # reuse: a second solver on the same container with IterationMax = 0
        s0, _ = _build({"flex_hartree_fock": True, "IterationMax": 0})
        gi["chiq_s"] = np.zeros(3)                                             # stale, must vanish
        s0.solve(gi, "tests/flex/output")
        for k in ("chi0q", "chiq_s", "chiq_c"):
            self.assertNotIn(k, gi)
        for k in ("sigma", "green", "physics"):
            self.assertIn(k, gi)
        self.assertEqual((s0.scf_iterations, s0.map_iteration, s0.state_iteration), (0, 0, 0))
        self.assertTrue(np.isnan(s0.scf_green_residual))
        with tempfile.TemporaryDirectory() as tmp:
            s0.save_results({"path_to_output": tmp, "sigma": "s.npz", "green": "g.npz", "chiq": "c.npz"}, gi)
            d = np.load(os.path.join(tmp, "s.npz"))
            for k in ("sigma_convention", "sigma_static", "sigma_fluct", "scf_converged", "scf_iterations",
                      "map_iteration", "state_iteration", "scf_sigma_residual", "scf_green_residual",
                      "scf_component_residual", "payload_kind", "hf_density_error", "hf_density_source",
                      "density_target_enforced"):
                self.assertIn(k, d.files, k)
            self.assertEqual(str(d["payload_kind"]), "final_state")
            self.assertEqual(str(d["sigma_convention"]), "split")
            self.assertFalse(os.path.exists(os.path.join(tmp, "c.npz")) and "chiq_s" in np.load(os.path.join(tmp, "c.npz")).files)
        # failure atomicity: a refused re-solve leaves nothing produced behind
        s_bad, _ = _build({"flex_hartree_fock": True})
        gi["trans_mod"] = np.zeros((s.lattice.nvol, 2 * s.norb, 2 * s.norb), complex)
        gi["sigma_init_envelope"] = s._read_sigma.__self__ and None
        gi.pop("sigma_init_envelope")
        import hwave.solver.flex_hf as flex_hf
        env = flex_hf.SigmaSeedEnvelope(sigma=np.zeros(1), sigma_static=None, sigma_fluct=None,
                                        marker="split", ir_meta=None, file_name="x")
        gi["sigma_init_envelope"] = env
        with self.assertRaises(ValueError):
            s_bad.solve(gi, "tests/flex/output")
        for k in ("sigma", "green", "physics", "chi0q"):
            self.assertNotIn(k, gi)
        self.assertFalse(hasattr(s_bad, "physics"))

    def test_hf_off_does_not_write_new_members(self):
        s, r = _build({})
        gi = r.get_param("green")
        os.makedirs("tests/flex/output", exist_ok=True)
        s.solve(gi, "tests/flex/output")
        with tempfile.TemporaryDirectory() as tmp:
            s.save_results({"path_to_output": tmp, "sigma": "s.npz", "chi0q": "c.npz"}, gi)
            self.assertNotIn("sigma_convention", np.load(os.path.join(tmp, "s.npz")).files)
            self.assertNotIn("payload_kind", np.load(os.path.join(tmp, "c.npz")).files)


class TestG5UHFkFixedPoint(unittest.TestCase):

    def test_w_disabled_fixed_point_equals_paramagnetic_uhfk(self):
        """Fluctuations disabled: the converged FLEX Hartree-Fock state equals
        a paramagnetic UHFk run (flag_fock=true, 2Sz=0, same T, per-cell
        normalisation) at the fixed point."""
        import hwave.solver.flex as flex_mod
        from hwave.solver.uhfk import UHFk
        # FLEX: seed sigma_static = HF map of the bare paramagnetic density (through trans_mod)
        s, r = _build({"flex_hartree_fock": True, "IterationMax": 400, "Mix": 0.5, "EPS": 10},
                      path=_IN2, interactions={"CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"})
        real = flex_mod.FLEX._calc_self_energy_general
        flex_mod.FLEX._calc_self_energy_general = lambda self_, g, v, b: np.zeros_like(g)
        gi = r.get_param("green")
        try:
            s.solve(gi, "tests/flex/output")
        finally:
            flex_mod.FLEX._calc_self_energy_general = real
        self.assertTrue(s.scf_converged)
        # UHFk on the same input, paramagnetic (2Sz = 0), T = 2.0, Ncond = filling * 2 * norb * nvol
        nvol, norb = s.lattice.nvol, s.norb
        idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"}
        r2 = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
        par = {"T": 2.0, "Ncond": int(round(0.5 * 2 * norb * nvol)), "2Sz": 0, "flag_fock": True,
               "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "IterationMax": 2000, "EPS": 12, "Mix": 0.5,
               "RndSeed": 1}
        u = UHFk(r2.get_param("ham"), {"print_level": 1, "print_step": 1}, {"mode": "UHFk", "param": par})
        with tempfile.TemporaryDirectory() as tmp:
            u.solve(r2.get_param("green"), tmp)
        G = u.Green                                                   # (nvol, 2, norb, 2, norb)
        # the FLEX final-state density
        from hwave.solver import hartree_fock as hf
        heff = hf.heff_eigenpairs(np.asarray(s.H0_k)[0], gi["sigma_static"][0, 0])
        dens = hf.equal_time_density(gi["green"], heff, s.mu, 1.0 / 2.0, (4, 4, 1))
        np.testing.assert_allclose(G[:, 0, :, 0, :], dens.rho_r, rtol=0, atol=1e-8)
        np.testing.assert_allclose(G[:, 1, :, 1, :], dens.rho_r, rtol=0, atol=1e-8)
        self.assertLess(np.abs(G[:, 0, :, 1, :]).max(), 1e-10)
        h_uhf = np.asarray(u.ham).reshape(nvol, 2, norb, 2, norb)[:, 0, :, 0, :]
        h_flex = np.asarray(s.H0_k)[0] + gi["sigma_static"][0, 0]
        np.testing.assert_allclose(h_flex, h_uhf, rtol=0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
