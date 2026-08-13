#!/usr/bin/env python3

"""Round-off gate for the Lehmann-machinery extraction (#151).

tests/capture_ed_reference.py captured ten arrays from
tests/test_rpa_vs_ed_oracle.py BEFORE the extraction refactor, into
tests/ed_oracle_reference.npz. This module pins every one of those ten
arrays against the refactored ed_oracle_util machinery (and, for the
declaration types the canonical density-term builder does not cover, the
rewired oracle module's own legacy builders), so a coordinated
utility/map change cannot slip through compensating errors.

Also guards the specific regression class the rewiring is prone to: the
legacy names (FX2/C/CD/H1) inside tests/test_rpa_vs_ed_oracle.py are lazy,
built on first use via a module-level ``_fx2_state()``. PEP 562's module
``__getattr__`` only intercepts EXTERNAL attribute access (``oracle.H1``
from outside the module, or after a failed normal lookup); a function
DEFINED INSIDE that module which reads a bare global ``H1`` resolves it
via ``LOAD_GLOBAL`` against the module's own namespace directly, which
``__getattr__`` never sees. If any internal builder in the oracle module
were to regress to a bare global read, it would raise ``NameError`` -- but
only once nothing else had already populated that global, so an in-process
test run (where some earlier test already triggered ``_fx2_state()`` and
usually nothing writes the bare names back) would not catch it. Only a
FRESH interpreter, calling the builders immediately, does. Hence the
subprocess test below.

Tests must be run from the repository root.
"""

import os
import subprocess
import sys
import unittest

import numpy as np

from tests.approx_util import ApproxTestCase, assert_approx_array
import tests.ed_oracle_util as ed_oracle_util
import tests.test_rpa_vs_ed_oracle as oracle


_REF_PATH = os.path.join(os.path.dirname(__file__), "ed_oracle_reference.npz")
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class TestExtractionGate(ApproxTestCase):
    """One test per captured artifact -- ALL ten -- through approx_util
    assertions only (np.testing.assert_allclose is forbidden in new code:
    it uses a looser sum-rule tolerance and treats NaN as equal, see
    tests/approx_util.py)."""

    @classmethod
    def setUpClass(cls):
        cls.U = ed_oracle_util
        cls.fx = ed_oracle_util.EDFixture(
            L=oracle.LX, norb=oracle.NORB, t=oracle.HOP,
            eps=oracle.E_ONSITE, T=oracle.T, mu=oracle.MU)
        with np.load(_REF_PATH) as ref:
            cls.ref = {k: ref[k] for k in ref.files}

    def test_all_ten_keys_present(self):
        expected = {"h1", "hf_CoulombIntra", "hf_CoulombInter", "hf_Exchange",
                    "hf_PairHop", "base_chi", "first_CoulombIntra",
                    "first_CoulombInter", "kernel_sample",
                    "base_chi_slotmapped"}
        self.assertEqual(set(self.ref), expected)

    def test_h1_identical(self):
        got = self.fx.build_h1()
        assert_approx_array(got, self.ref["h1"], rel=0, abs=0.0)

    def test_kernel_identical(self):
        got = self.U.lehmann_kernel(self.fx)      # exposed for the gate
        assert_approx_array(got, self.ref["kernel_sample"], rel=0, abs=1e-13)

    def test_base_chi_identical(self):
        got = self.U.chi_connected(self.fx)
        assert_approx_array(got, self.ref["base_chi"], rel=0, abs=1e-12)

    def test_slot_map_identical(self):
        got = self.U.to_solver_slots(self.ref["base_chi"])
        assert_approx_array(got, self.ref["base_chi_slotmapped"], rel=0, abs=0.0)

    def test_hf_all_kinds_identical(self):
        # CoulombIntra via canonical terms; Exchange/PairHop via the legacy
        # builders re-exported by the rewired oracle module
        terms = self.U.canonical_density_terms(self.fx, [(0, 0, 0, oracle.V1)])
        got = self.U.hf_h1_from_terms(self.fx, terms)
        assert_approx_array(got, self.ref["hf_CoulombIntra"], rel=0, abs=1e-12)

        terms = self.U.canonical_density_terms(self.fx, [(0, 1, 0, oracle.V1)])
        got = self.U.hf_h1_from_terms(self.fx, terms)
        assert_approx_array(got, self.ref["hf_CoulombInter"], rel=0, abs=1e-12)

        got = oracle._hf_h1("Exchange", oracle.V1)
        assert_approx_array(got, self.ref["hf_Exchange"], rel=0, abs=1e-12)

        got = oracle._hf_h1("PairHop", oracle.V1)
        assert_approx_array(got, self.ref["hf_PairHop"], rel=0, abs=1e-12)

    def test_h_int_matches_legacy_na_nb(self):
        # direct pin on H_int itself (not just its first-order response):
        # the (0, 1, R=0, V) canonical class must reproduce the legacy
        # _h_int's na@nb exactly.
        terms = self.U.canonical_density_terms(self.fx, [(0, 1, 0, oracle.V1)])
        assert_approx_array(self.U.h_int_from_terms(self.fx, terms),
                            oracle._h_int("CoulombInter", oracle.V1),
                            rel=0, abs=1e-14)

    def test_h_int_matches_legacy_coulombintra(self):
        terms = self.U.canonical_density_terms(self.fx, [(0, 0, 0, oracle.V1)])
        assert_approx_array(self.U.h_int_from_terms(self.fx, terms),
                            oracle._h_int("CoulombIntra", oracle.V1),
                            rel=0, abs=1e-14)

    def test_first_order_both_kinds_identical(self):
        # CoulombIntra AND CoulombInter, via U.richardson over the
        # canonical-terms vertex function
        cases = {
            "first_CoulombIntra": (0, 0, 0),
            "first_CoulombInter": (0, 1, 0),
        }
        for key, (a, b, r) in cases.items():
            with self.subTest(key=key):
                def vertex(v, a=a, b=b, r=r):
                    terms = self.U.canonical_density_terms(
                        self.fx, [(a, b, r, v)])
                    hint = self.U.h_int_from_terms(self.fx, terms)
                    h1 = self.U.hf_h1_from_terms(self.fx, terms)
                    full = self.U.chi_connected(self.fx, hint=hint)
                    ins = self.U.chi_connected(self.fx, h1=h1)
                    return self.U.to_solver_slots(full - ins)

                got = self.U.richardson(vertex, oracle.V1)
                assert_approx_array(got, self.ref[key], rel=0, abs=1e-10)


class TestFreshImportHasNoNameError(unittest.TestCase):
    """The legacy names (FX2/C/CD/H1) and legacy builders (_h_int,
    _hf_h1, _chi_ed, _ed_first_order) must all work on a FRESH interpreter
    that has done nothing but ``import tests.test_rpa_vs_ed_oracle`` --
    see the module docstring for why this cannot be checked in-process."""

    def test_legacy_builders_work_on_fresh_import(self):
        code = "\n".join([
            "import tests.test_rpa_vs_ed_oracle as m",
            "m._h_int('Exchange', 0.1)",
            "m._h_int('CoulombIntra', 0.1)",
            "m._h_int('CoulombInter', 0.1)",
            "m._h_int('PairLift', 0.1)",
            "m._h_int('PairHop', 0.1)",
            "m._hf_h1('Exchange', 0.1)",
            "m._hf_h1('CoulombIntra', 0.1)",
            "m._hf_h1('CoulombInter', 0.1)",
            "m._hf_h1('PairLift', 0.1)",
            "m._hf_h1('PairHop', 0.1)",
            "m._chi_ed()",
            "m._ed_first_order('CoulombIntra')",
            "m.FX2",
            "m.C",
            "m.CD",
            "m.H1",
            "print('OK')",
        ])
        result = subprocess.run(
            [sys.executable, "-c", code], cwd=_REPO_ROOT,
            capture_output=True, text=True)
        self.assertEqual(
            result.returncode, 0,
            "subprocess import/build failed:\nstdout={}\nstderr={}".format(
                result.stdout, result.stderr))
        self.assertIn("OK", result.stdout)


if __name__ == "__main__":
    unittest.main()
