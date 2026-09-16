"""Compatibility contract (spec 2026-09-08 section 3, extended by spec
2026-09-16 section 2.4): 'takimoto' without the Hartree-Fock term keeps
every numerical archive member of the general path np.array_equal to the
reference revision's, and so does 'local' on an input with no off-site
two-body row; the Hartree-Fock map on a declaration F reproduces the
reference revision's map on the REVERSED declaration F^rev exactly; 'local'
on CoulombIntra-only input agrees with 'takimoto' to 1e-14 with the
archive-scale floor; the factors are built at construction under 'local'
only; the D7 refusal fires at construction.

What is NOT here: the paths the interaction-row orientation change moved
(the Hartree-Fock term, the local kernel's off-site vertex, the bond gate).
They cannot be compared to the reference revision on the same declaration
by construction; their numerical state is pinned as committed fixed state
by tests/test_flex_orientation_baseline.py, and their CORRECTNESS by the
ED/oracle modules that file names."""
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


def _run(checkout, path, inter, out, so, extra=None):
    env = dict(os.environ, PYTHONPATH=os.path.join(checkout, "src") + ":" + checkout)
    subprocess.run([sys.executable, "-B", "-c", _RUN, path, json.dumps(inter), out, so,
                    json.dumps(extra or {})],
                   env=env, check=True, capture_output=True, cwd=checkout)


#: The commit the compatibility harnesses compare against: the merge that
#: precedes this branch (``develop`` at the time the branch was cut -- the
#: merge of the local second-order kernel, #191). The contract these
#: harnesses assert is BYTE identity with that revision, so the revision has
#: to be named -- comparing against "whatever the reference directory
#: happens to contain" turns a failed comparison into a puzzle and a passing
#: one into nothing at all.
#:
#: The reference revision already HAS ``flex_second_order`` (it is the merge
#: that added it) and already defaults it to ``"local"``, so every harness
#: below names the kernel on BOTH sides. Leaving the reference side to its
#: default would compare ``"local"`` against ``"takimoto"`` and report the
#: difference between two kernels as a compatibility break.
DEVELOP_COMMIT = "add6dc44d930e1c11cdfb3c70df7ae28bd3b95be"


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


#: the two-body types whose displacement tables ``build_interaction_tables``
#: reads in the documented orientation (spec 2026-09-16 D-1); ``CoulombIntra``
#: is on-site by construction and is not one of them
_ORIENTED_TYPES = ("CoulombInter", "Hund", "Ising", "PairLift", "Exchange", "PairHop")


def _hermitian_density(shape, norb, seed):
    """A deterministic per-spin real-space density ``rho_ab(r)`` obeying the
    convention the Hartree-Fock map requires, ``rho_ab(r) = conj(rho_ba(-r))``
    (it refuses anything else): the transform of a random HERMITIAN k-space
    density."""
    nvol = int(np.prod(shape))
    rng = np.random.default_rng(seed)
    rho_k = (rng.normal(size=(nvol, norb, norb))
             + 1j * rng.normal(size=(nvol, norb, norb)))
    rho_k = 0.5 * (rho_k + np.conjugate(np.swapaxes(rho_k, -1, -2)))
    return np.fft.ifftn(rho_k.reshape(*shape, norb, norb),
                        axes=(0, 1, 2)).reshape(nvol, norb, norb)


def _legacy_tables(tables):
    """``tables`` with the documented-orientation step UNDONE, i.e. the
    tables the reference revision builds from the same declaration.

    :func:`hwave.solver.hartree_fock.orient_documented` is an involution
    (the conjugate transpose at fixed ``r``, with ``r = 0`` untouched), so
    applying it a second time to each oriented type restores the reference
    revision's table exactly -- no second source tree and no
    re-implementation of the old builder."""
    from hwave.solver import hartree_fock as hf
    inter = dict(tables.inter_table)
    for t in _ORIENTED_TYPES:
        if inter.get(t) is not None:
            inter[t] = hf.orient_documented(inter[t])
    return hf.InteractionTables(inter, tables.spin_table, tables.discarded)


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
        """The paths this branch did NOT move: every numerical archive
        member of a ``"takimoto"`` run without the Hartree-Fock term is
        ``np.array_equal`` to the reference revision's, on a fixture with
        INTER-ORBITAL off-site rows and on an on-site one.

        Scope, after the interaction-row orientation change (issue #192,
        spec 2026-09-16 section 2.4): reading an off-site row ``(r, a, b, v)``
        in the documented orientation changes three paths -- the FLEX
        Hartree-Fock term (``flex_hartree_fock``), the local second-order
        kernel's off-site vertex (``flex_second_order = "local"``) and the
        bond gate, which reads the rows through the Hartree-Fock term. The
        legacy kernel ``"takimoto"`` without the Hartree-Fock term is not
        one of them: its vertex is the q-space density block, which the
        orientation of the real-space row does not enter. That is the
        identity asserted here, and its counterpart -- what the moved paths
        DO produce -- is recorded as fixed state in
        ``tests/test_flex_orientation_baseline.py`` (the gate-on case that
        used to sit in this list is its ``smoke`` vector, under both
        kernels) and pinned against develop at the Hartree-Fock map itself
        by :meth:`test_hf_map_equals_develop_on_the_reversed_declaration`
        below and by
        ``tests/test_flex_hf_scf.py::TestG0Off::
        test_first_map_static_equals_develop_on_the_reversed_declaration``.
        """
        dev = self._develop()
        here = os.getcwd()
        compared = 0
        for path, inter, extra, files in (
                (_IN2, {"CoulombInter": "coulombinter.dat"}, None, _FILES),
                (_IN2, {"CoulombInter": "onsite_inter.dat",
                        "CoulombIntra": "coulombintra.dat"}, None, _FILES)):
            with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
                # "takimoto" on BOTH sides: the reference revision defaults
                # the key to "local" (see DEVELOP_COMMIT)
                _run(dev, os.path.abspath(path), inter, a, "takimoto", extra)
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

    def test_local_on_onsite_only_input_equals_develop(self):
        """The other unmoved path: ``flex_second_order = "local"`` on an
        input with NO off-site two-body row is byte-identical to the
        reference revision's.

        The orientation change touches the kernel's OFF-SITE vertex only
        (``second_order.accumulate_batch``'s ``vpair`` terms), so on an
        on-site fixture the exact local kernel must be untouched. Without
        this case the identity above would leave the impression that only
        the legacy kernel is unchanged."""
        dev = self._develop()
        here = os.getcwd()
        inter = {"CoulombInter": "onsite_inter.dat", "CoulombIntra": "coulombintra.dat"}
        with tempfile.TemporaryDirectory() as a, tempfile.TemporaryDirectory() as b:
            _run(dev, os.path.abspath(_IN2), inter, a, "local")
            _run(here, os.path.abspath(_IN2), inter, b, "local")
            ma, mb = _members(a), _members(b)
            compared = 0
            for f in ma:
                for k, v in ma[f].items():
                    self.assertIn(k, mb[f], (f, k))
                    if np.asarray(v).dtype.kind in "fciub":
                        np.testing.assert_array_equal(np.asarray(mb[f][k]), np.asarray(v),
                                                      err_msg=str((f, k)))
                        compared += 1
            self.assertGreater(compared, 10)
            # anti-vacuity: the kernel really ran on a non-trivial sigma
            self.assertGreater(np.abs(np.asarray(ma["sigma.npz"]["sigma"])).max(), 1e-6)

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

    def test_hf_map_equals_develop_on_the_reversed_declaration(self):
        """Spec 2.4: ``hf_map_new(F, rho) == hf_map_develop(F^rev, rho)`` at a
        fixed density (``np.array_equal``).

        This is the COMPATIBILITY statement of the orientation change for the
        mean-field term: nothing is lost, the same physics is now written the
        documented way round. An off-site row ``(r, a, b, v)`` that used to
        mean ``v n_{j,b} n_{j+r,a}`` means ``v n_{j,a} n_{j+r,b}`` here, so a
        user who wants the OLD result reverses the displacement of every
        off-site row (``F -> F^rev``) and gets it back exactly.

        The reference map is reproduced in process rather than through a
        second source tree: it is today's map on the LEGACY tables of
        ``F^rev``, i.e. ``build_flex_hf_tables(F^rev)`` with the orientation
        step undone (:func:`_legacy_tables`; ``orient_documented`` is an
        involution). The end-to-end version of the same statement, against
        the real reference tree and through ``solve``, is
        ``tests/test_flex_hf_scf.py::TestG0Off::
        test_first_map_static_equals_develop_on_the_reversed_declaration``.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        from hwave.solver import flex_hf
        from hwave.solver.kgrid import reverse_fft_axes
        shape, norb = (4, 4, 1), 2
        # the fixture's declaration, read once; its off-site content carries
        # the INTER-ORBITAL rows v_12(-x) = 1 / v_21(+x) = 1, which is what
        # makes F^rev a different declaration at all (reversing the
        # displacement of an orbital-DIAGONAL row is the identity on the
        # Hermitian-closed table)
        idict = {"path_to_input": _IN2, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "CoulombInter": "coulombinter.dat"}
        rows = read_input_k.QLMSkInput(
            {"path_to_input": _IN2, "interaction": idict}).get_param("ham")["CoulombInter"]
        F = {"CoulombInter": dict(rows)}
        F_rev = {"CoulombInter": {((-ir[0], -ir[1], -ir[2]), ov): v
                                  for (ir, ov), v in rows.items()}}
        self.assertNotEqual(F["CoulombInter"], F_rev["CoulombInter"])   # anti-vacuity
        rho = _hermitian_density(shape, norb, seed=20260916)
        # the density really is in the convention hf_map requires
        rev = np.conjugate(np.swapaxes(
            reverse_fft_axes(rho.reshape(*shape, norb, norb), (0, 1, 2)), -1, -2))
        self.assertLess(np.abs(rho - rev.reshape(rho.shape)).max(), 1e-14)
        t_new = flex_hf.build_flex_hf_tables(F, norb, shape)
        t_rev = flex_hf.build_flex_hf_tables(F_rev, norb, shape)
        new = flex_hf.hf_map(rho, t_new, shape, norb)
        reference = flex_hf.hf_map(rho, _legacy_tables(t_rev), shape, norb)
        self.assertGreater(np.abs(new).max(), 1e-6)                    # anti-vacuity
        np.testing.assert_array_equal(new, reference)
        # anti-vacuity of the REVERSAL: today's map on F differs from the
        # reference map on F (the legacy tables of the SAME declaration), so
        # the equality above is a property of F^rev, not of the map
        legacy_on_F = flex_hf.hf_map(rho, _legacy_tables(t_new), shape, norb)
        gap = np.abs(new - legacy_on_F).max() / np.abs(new).max()
        self.assertGreater(gap, 1e-3,
                           "the orientation step does not change this fixture's "
                           "Hartree-Fock map ({:.3e} relative); the equivalence "
                           "above would be vacuous".format(gap))

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
