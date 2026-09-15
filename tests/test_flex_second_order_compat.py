"""Compatibility contract (spec 2026-09-08 section 3): 'takimoto' keeps every
numerical archive member of the general path np.array_equal to develop's;
'local' on CoulombIntra-only input agrees with 'takimoto' to 1e-14 with the
archive-scale floor; the factors are built at construction under 'local'
only; the D7 refusal fires at construction."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

_IN2 = "tests/rpa/input_2orb"
_IN1 = "tests/rpa/input"

_RUN = r'''
import json, os, sys, numpy as np
import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.flex as flex_mod
path, inter_json, out, so, extra_json = sys.argv[1:6]
inter = json.loads(inter_json)
idict = {"path_to_input": path, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
idict.update(inter)
r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 16,
       "IterationMax": 3, "Mix": 0.5, "EPS": 1e-12}
par.update(json.loads(extra_json))
if so != "absent":
    par["flex_second_order"] = so
s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                            "enable_spin_orbital": False, "calc_scheme": "general"})
gi = r.get_param("green")
s.solve(gi, out)
s.save_results({"path_to_output": out, "chi0q": "chi0q", "chiq": "chiq", "sigma": "sigma",
                "green": "green", "energy": "energy.dat"}, gi)
'''

_FILES = ("chi0q.npz", "chiq.npz", "chiq_s.npz", "chiq_c.npz", "sigma.npz", "green.npz")
#: the bond gate adds its dedicated archive to the comparison
_BOND_FILES = _FILES + ("longitudinal_bond.npz",)
#: gate-on parameters: the Hartree-Fock term, the bond channels and the
#: dynamic archive (all three exist on develop as well, #181 Phase B)
_GATE = {"flex_hartree_fock": True, "longitudinal_bond_channels": True,
         "longitudinal_bond_output_full": True}


def _run(checkout, path, inter, out, so, extra=None):
    env = dict(os.environ, PYTHONPATH=os.path.join(checkout, "src") + ":" + checkout)
    subprocess.run([sys.executable, "-B", "-c", _RUN, path, json.dumps(inter), out, so,
                    json.dumps(extra or {})],
                   env=env, check=True, capture_output=True, cwd=checkout)


#: The commit the compatibility harnesses compare against: the merge that
#: precedes this branch (``develop`` at the time the branch was cut). The
#: contract these harnesses assert is BYTE identity with that revision, so
#: the revision has to be named -- comparing against "whatever the reference
#: directory happens to contain" turns a failed comparison into a puzzle and
#: a passing one into nothing at all.
DEVELOP_COMMIT = "59ac623810f2cff841d96bc721c9312b7ad912f7"


def develop_checkout(run=None):
    """``(path, None)`` when the reference checkout is usable, ``(None,
    reason)`` otherwise.

    Usable means: it exists, it is at :data:`DEVELOP_COMMIT`, and its tree
    is clean. A checkout at another revision -- or with local edits -- is
    NOT a reference: the comparison would either fail for reasons that have
    nothing to do with this branch, or pass against a tree nobody can name.
    Either way the answer is to skip and say what was expected.

    ``run`` is the subprocess runner, injectable so the rejections can be
    unit-tested without a second checkout."""
    run = run or subprocess.run
    dev = os.environ.get("HWAVE_DEVELOP_CHECKOUT",
                         os.path.abspath(os.path.join(os.getcwd(), "..", "..", "..")))
    if not os.path.exists(os.path.join(dev, "src", "hwave", "solver", "flex.py")):
        return None, ("reference checkout not found at {} (set HWAVE_DEVELOP_CHECKOUT)"
                      .format(dev))
    try:
        head = run(["git", "-C", dev, "rev-parse", "HEAD"],
                   check=True, capture_output=True, text=True).stdout.strip()
        dirty = run(["git", "-C", dev, "status", "--porcelain"],
                    check=True, capture_output=True, text=True).stdout.strip()
    except (OSError, subprocess.SubprocessError) as exc:
        return None, "cannot read the revision of the reference checkout {}: {}".format(dev, exc)
    if head != DEVELOP_COMMIT:
        return None, ("the reference checkout {} is at {}, but this comparison is against {}"
                      .format(dev, head or "<unknown>", DEVELOP_COMMIT))
    if dirty:
        return None, ("the reference checkout {} has local modifications; this comparison "
                      "needs a clean tree at {}".format(dev, DEVELOP_COMMIT))
    return dev, None


def _members(d, files=_FILES):
    out = {}
    for f in files:
        z = np.load(os.path.join(d, f))
        out[f] = {k: z[k] for k in z.files}
    return out


class TestCompatibility(unittest.TestCase):

    def _develop(self):
        dev, why = develop_checkout()
        if dev is None:
            self.skipTest(why)
        return dev

    def test_takimoto_numerical_members_equal_develop(self):
        dev = self._develop()
        here = os.getcwd()
        compared = 0
        for path, inter, extra, files in (
                (_IN2, {"CoulombInter": "coulombinter.dat"}, None, _FILES),
                (_IN2, {"CoulombInter": "onsite_inter.dat",
                        "CoulombIntra": "coulombintra.dat"}, None, _FILES),
                # the bond gate (flex_hartree_fock + longitudinal_bond_channels):
                # its channel-0 second order is selected by flex_second_order
                # too, so "takimoto" must keep the whole gate-on archive set --
                # the dedicated bond archive included -- equal to develop's
                (_IN2, {"CoulombInter": "coulombinter.dat"}, _GATE, _BOND_FILES)):
            with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
                _run(dev, os.path.abspath(path), inter, a, "absent", extra)
                _run(here, os.path.abspath(path), inter, b, "takimoto", extra)
                ma, mb = _members(a, files), _members(b, files)
                for f in ma:
                    for k, v in ma[f].items():
                        self.assertIn(k, mb[f], (f, k))
                        # "b" too: boolean members such as density_target_enforced
                        # are part of the archive and must not drift either
                        if np.asarray(v).dtype.kind in "fciub":
                            np.testing.assert_array_equal(np.asarray(mb[f][k]), np.asarray(v),
                                                          err_msg=str((f, k)))
                            compared += 1
        # anti-vacuity: the member list must not be empty
        self.assertGreater(compared, 10)

    def test_local_equals_takimoto_on_hubbard_only(self):
        here = os.getcwd()
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            _run(here, os.path.abspath(_IN1), {"CoulombIntra": "coulombintra.dat"}, a, "takimoto")
            _run(here, os.path.abspath(_IN1), {"CoulombIntra": "coulombintra.dat"}, b, "local")
            ma, mb = _members(a), _members(b)
            compared = 0
            for f in ma:
                scale = max([1.0] + [np.abs(np.asarray(v)).max() for v in ma[f].values()
                                     if np.asarray(v).dtype.kind in "fc"])
                for k, v in ma[f].items():
                    if np.asarray(v).dtype.kind in "fc":
                        np.testing.assert_allclose(np.asarray(mb[f][k]), np.asarray(v), rtol=1e-14,
                                                   atol=1e-14 * scale, err_msg=str((f, k)))
                        compared += 1
            self.assertGreater(compared, 5)
            # anti-vacuity: the kernel actually ran on a non-trivial sigma
            self.assertGreater(np.abs(np.asarray(ma["sigma.npz"]["sigma"])).max(), 1e-6)

    def test_calc_veff_general_local_equals_takimoto_on_hubbard(self):
        """The Hubbard-only identity of spec section 3 at the METHOD level:
        with CoulombIntra as the only interaction the local kernel's
        3/2 Us (chi_s - chibar) Us + 1/2 Uc (chi_c - chibar) Uc + W2 equals the
        legacy assembly element-wise. In-process, no subprocess, and it is the
        only direct exercise of _calc_veff_general(second_order="local") --
        every other caller in the suite omits the keyword and takes
        "takimoto"."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod
        beta = 0.5
        for path, norb in ((_IN1, 1), (_IN2, 2)):
            with self.subTest(path=path, norb=norb):
                idict = {"path_to_input": path, "Geometry": "geom.dat",
                         "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
                r = read_input_k.QLMSkInput({"path_to_input": path, "interaction": idict})
                par = {"T": 1.0 / beta, "mu": 0.1, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
                       "Nmat": 8, "IterationMax": 1, "Mix": 1.0, "EPS": 1}
                s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                            "enable_spin_orbital": False,
                                                            "calc_scheme": "general"})
                self.assertEqual(s.flex_second_order, "local")
                self.assertIsNotNone(s._second_order_factors)
                self.assertEqual(s.norb, norb)
                s._calc_epsilon_k({})
                G = s._calc_dressed_green(
                    beta, 0.1, np.zeros((1, s.nmat, s.lattice.nvol, norb, norb), complex))
                chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), beta)[0]
                chi0q, Us, Uc = s._inflate_chi0q_and_ham_general(
                    chi0q_raw, s.ham_info.ham_inter_q)
                chi_s, chi_c = s._solve_channels_general(chi0q, Us, Uc)
                v_thu = s._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc,
                                             second_order="takimoto")
                v_loc = s._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc,
                                             factors=s._second_order_factors,
                                             second_order="local")
                scale = np.abs(v_thu).max()
                self.assertGreater(scale, 1e-6)        # anti-vacuity
                self.assertEqual(v_loc.shape, v_thu.shape)
                np.testing.assert_allclose(v_loc, v_thu, rtol=1e-14, atol=1e-14 * scale)
                # the two branches really are different code paths
                with self.assertRaises(ValueError):
                    s._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc, second_order="local")
                with self.assertRaises(ValueError):
                    s._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc, second_order="nonsense")

    def test_calc_veff_general_batch_path_equals_the_dense_assembly(self):
        """The local branch's frequency batching (spec 2.3): whatever nb the
        batch loop picks, the result must equal the dense assembly

            3/2 Us (chi_s - chibar) Us + 1/2 Uc (chi_c - chibar) Uc + W2(chibar)

        formed on the FULL frequency axis. Two batch shapes are exercised --
        nmat = 28 gives nb = 3 and a final batch of length 1, nmat = 6 gives
        nb = 1 -- on a two-orbital input WITH off-site V (the kernel's
        off-site branch runs) and on one without (on-site factors only)."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_mod
        from hwave.solver.second_order import dense_w2
        beta = 0.5
        for inter, offsite in (({"CoulombInter": "coulombinter.dat"}, True),
                               ({"CoulombIntra": "coulombintra.dat",
                                 "CoulombInter": "onsite_inter.dat"}, False)):
            for nmat in (28, 6):
                with self.subTest(offsite=offsite, nmat=nmat):
                    idict = {"path_to_input": _IN2, "Geometry": "geom.dat",
                             "Transfer": "transfer.dat"}
                    idict.update(inter)
                    r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
                    par = {"T": 1.0 / beta, "mu": 0.1, "CellShape": [4, 4, 1],
                           "SubShape": [1, 1, 1], "Nmat": nmat, "IterationMax": 1,
                           "Mix": 1.0, "EPS": 1}
                    s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                                "enable_spin_orbital": False,
                                                                "calc_scheme": "general"})
                    f = s._second_order_factors
                    self.assertEqual(f.vpair is not None, offsite)
                    self.assertEqual(max(1, s.nmat // 8), 3 if nmat == 28 else 1)
                    s._calc_epsilon_k({})
                    G = s._calc_dressed_green(
                        beta, 0.1, np.zeros((1, s.nmat, s.lattice.nvol, s.norb, s.norb), complex))
                    chi0q, Us, Uc = s._inflate_chi0q_and_ham_general(
                        s._calc_chi0q(G, np.zeros_like(G), beta)[0], s.ham_info.ham_inter_q)
                    chi_s, chi_c = s._solve_channels_general(chi0q, Us, Uc)
                    v = s._calc_veff_general(chi0q, chi_s, chi_c, Us, Uc, factors=f,
                                             second_order="local")
                    ndx = s.norb ** 2
                    shape = (s.nmat, s.lattice.nvol, ndx, ndx)
                    c0 = chi0q.reshape(shape)
                    UsB, UcB = Us[np.newaxis], Uc[np.newaxis]
                    ref = (1.5 * (UsB @ (chi_s.reshape(shape) - c0) @ UsB)
                           + 0.5 * (UcB @ (chi_c.reshape(shape) - c0) @ UcB)
                           + dense_w2(c0, f))
                    scale = np.abs(ref).max()
                    self.assertGreater(scale, 1e-6)                 # anti-vacuity
                    np.testing.assert_allclose(v, ref, rtol=0, atol=1e-13 * scale)

    def test_factors_lifecycle_and_d7_at_construction(self):
        from tests.test_second_order_factors import _split_for
        import hwave.solver.flex as flex_mod
        import hwave.qlmsio.read_input_k as read_input_k
        from hwave.solver.second_order import DegenerateRowError
        s, _ = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0)]})   # takimoto in the helper
        self.assertIsNone(s._second_order_factors)
        self.assertIsNone(s._second_order_device)
        idict = {"path_to_input": _IN2, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "CoulombInter": "coulombinter.dat"}
        r = read_input_k.QLMSkInput({"path_to_input": _IN2, "interaction": idict})
        par = {"T": 2.0, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 8,
               "IterationMax": 1, "Mix": 1.0, "EPS": 1}
        s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                    "enable_spin_orbital": False,
                                                    "calc_scheme": "general"})
        self.assertIsNotNone(s._second_order_factors)
        self.assertFalse(s._second_order_factors.A_on[0].flags.writeable)
        # D7 at construction, 'local' only
        with tempfile.TemporaryDirectory() as d:
            for f in ("geom.dat", "transfer.dat", "coulombintra.dat"):
                shutil.copy(os.path.join(_IN2, f), d)
            with open(os.path.join(d, "hund.dat"), "w") as fw:
                fw.write("Hund in wannier90-like format for uhfk\n2\n1\n 1\n"
                         "   0    0    0    1    1  0.3 0.0\n")
            idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                     "Hund": "hund.dat", "CoulombIntra": "coulombintra.dat"}
            r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
            with self.assertRaises(DegenerateRowError):
                flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                        "enable_spin_orbital": False,
                                                        "calc_scheme": "general"})
            par2 = dict(par, flex_second_order="takimoto")
            s2 = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par2,
                                                         "enable_spin_orbital": False,
                                                         "calc_scheme": "general"})
            self.assertIsNone(s2._second_order_factors)



class TestDevelopCheckoutGuard(unittest.TestCase):
    """:func:`develop_checkout` refuses a reference tree it cannot name.

    The refusals are unit-tested with a stubbed runner rather than a second
    checkout: what matters is that a wrong revision and a dirty tree both
    STOP the comparison and say what was expected, not that git works."""

    class _Result(object):
        def __init__(self, stdout):
            self.stdout = stdout

    def _runner(self, head, dirty):
        def run(cmd, **kw):
            return self._Result(head if "rev-parse" in cmd else dirty)
        return run

    def test_accepts_the_named_revision_with_a_clean_tree(self):
        dev, why = develop_checkout(run=self._runner(DEVELOP_COMMIT + "\n", ""))
        if dev is None:
            # the reference directory itself is missing; that is the one
            # rejection this stub cannot bypass
            self.assertIn("not found", why)
            self.skipTest(why)
        self.assertIsNone(why)

    def test_rejects_a_wrong_revision(self):
        wrong = "0" * 40
        dev, why = develop_checkout(run=self._runner(wrong, ""))
        if dev is None and "not found" in (why or ""):
            self.skipTest(why)
        self.assertIsNone(dev)
        self.assertIn(wrong, why)
        self.assertIn(DEVELOP_COMMIT, why)

    def test_rejects_a_dirty_tree(self):
        dev, why = develop_checkout(
            run=self._runner(DEVELOP_COMMIT, " M src/hwave/solver/flex.py"))
        if dev is None and "not found" in (why or ""):
            self.skipTest(why)
        self.assertIsNone(dev)
        self.assertIn("local modifications", why)
        self.assertIn(DEVELOP_COMMIT, why)

    def test_rejects_a_checkout_git_cannot_read(self):
        def run(cmd, **kw):
            raise OSError("git not found")
        dev, why = develop_checkout(run=run)
        if dev is None and "not found at" in (why or ""):
            self.skipTest(why)
        self.assertIsNone(dev)
        self.assertIn("cannot read the revision", why)


if __name__ == "__main__":
    unittest.main()
