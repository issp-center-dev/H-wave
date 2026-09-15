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


def _members(d, files=_FILES):
    out = {}
    for f in files:
        z = np.load(os.path.join(d, f))
        out[f] = {k: z[k] for k in z.files}
    return out


class TestCompatibility(unittest.TestCase):

    def _develop(self):
        dev = os.environ.get("HWAVE_DEVELOP_CHECKOUT",
                             os.path.abspath(os.path.join(os.getcwd(), "..", "..", "..")))
        if not os.path.exists(os.path.join(dev, "src", "hwave", "solver", "flex.py")):
            self.skipTest("develop checkout not found (set HWAVE_DEVELOP_CHECKOUT)")
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
                        if np.asarray(v).dtype.kind in "fciu":
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


if __name__ == "__main__":
    unittest.main()
