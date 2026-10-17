"""Unit tests for ``tests.approx_util``: pin the exact ``pytest.approx``
semantics the helper is meant to reproduce (see the module docstring there
for the full table), at the boundary values where a sum-rule-vs-max-rule or
NaN-handling slip would actually show up.

Boundary ("diff == tol") cases are built from binary-exact values (powers of
two) so ``actual - expected`` recovers the tolerance bit-for-bit under IEEE
754 -- not from adding/subtracting the (decimal, not binary-exact) default
1e-6/1e-12 constants, which would make an exact-equality assertion
non-deterministic at the last bit.
"""
import unittest

import numpy as np

from tests.approx_util import (
    ApproxTestCase,
    assert_approx,
    assert_approx_array,
)


class TestAssertApproxDefaults(unittest.TestCase):
    """rel/abs both omitted -> tol = max(1e-6 * |expected|, 1e-12)."""

    def test_absolute_floor_diff_equal_to_tolerance_passes(self):
        # rel * |0| == 0, so the 1e-12 absolute floor is the exact bound here,
        # and 0.0 - x is always exact -> a genuine bit-exact boundary case.
        assert_approx(1e-12, 0.0)

    def test_absolute_floor_diff_just_above_tolerance_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(1.000001e-12, 0.0)

    def test_relative_term_well_within_tolerance_passes(self):
        assert_approx(1.0 + 1e-9, 1.0)

    def test_relative_term_well_outside_tolerance_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(1.0 + 1e-3, 1.0)


class TestAssertApproxRelOnly(unittest.TestCase):
    """rel given, abs omitted -> tol = max(rel * |expected|, 1e-12)."""

    # rel and expected are powers of two so rel * expected (and the
    # subsequent +/- against expected) is exact under IEEE 754.
    _REL = 2.0 ** -4     # 0.0625
    _EXPECTED = 2.0 ** 10  # 1024.0
    _TOL = _REL * _EXPECTED  # 64.0, exact

    def test_diff_equal_to_tolerance_passes(self):
        actual = self._EXPECTED + self._TOL  # 1088.0, exact
        assert_approx(actual, self._EXPECTED, rel=self._REL)

    def test_diff_just_above_tolerance_fails(self):
        actual = self._EXPECTED + self._TOL + 1.0  # unambiguously past the tol
        with self.assertRaises(AssertionError):
            assert_approx(actual, self._EXPECTED, rel=self._REL)

    def test_absolute_floor_still_applies_for_tiny_expected(self):
        # rel * |expected| is far below 1e-12, so 1e-12 is the real bound
        assert_approx(1e-12, 0.0, rel=1e-6)
        with self.assertRaises(AssertionError):
            assert_approx(1.000001e-12, 0.0, rel=1e-6)


class TestAssertApproxAbsOnly(unittest.TestCase):
    """abs given, rel omitted -> tol = EXACTLY abs, no relative term mixed in."""

    _ABS = 2.0 ** -10  # exact binary value
    _EXPECTED = 5.0

    def test_diff_equal_to_tolerance_passes(self):
        actual = self._EXPECTED + self._ABS  # exact (Sterbenz-safe subtraction)
        assert_approx(actual, self._EXPECTED, abs=self._ABS)

    def test_diff_just_above_tolerance_fails(self):
        actual = self._EXPECTED + self._ABS + self._ABS  # 2x the tolerance
        with self.assertRaises(AssertionError):
            assert_approx(actual, self._EXPECTED, abs=self._ABS)

    def test_no_relative_term_leaks_in_for_a_large_expected(self):
        # if a relative term were mixed in, this would pass (1e-6 * 1e9 >> abs);
        # it must NOT, because abs was given explicitly.
        with self.assertRaises(AssertionError):
            assert_approx(1.0e9 + 1.0, 1.0e9, abs=1e-6)


class TestAssertApproxNanAndInf(unittest.TestCase):

    def test_nan_actual_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(float("nan"), 1.0)

    def test_nan_expected_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(1.0, float("nan"))

    def test_nan_actual_and_expected_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(float("nan"), float("nan"))

    def test_same_sign_infinity_passes(self):
        assert_approx(float("inf"), float("inf"))
        assert_approx(float("-inf"), float("-inf"))

    def test_opposite_sign_infinity_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(float("inf"), float("-inf"))

    def test_infinity_vs_finite_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(float("inf"), 1.0)
        with self.assertRaises(AssertionError):
            assert_approx(1.0, float("inf"))

    def test_exact_equality_short_circuits_before_any_tolerance_math(self):
        # rel=0, abs=0 would make the tolerance itself 0; the equality
        # short-circuit must still pass without even reaching that branch.
        assert_approx(3.0, 3.0, rel=0.0, abs=0.0)


class TestAssertApproxComplex(unittest.TestCase):

    def test_complex_within_tolerance_passes(self):
        assert_approx(1.0 + 2.0j, 1.0 + 2.0j + 1e-13)

    def test_complex_outside_tolerance_fails(self):
        with self.assertRaises(AssertionError):
            assert_approx(1.0 + 2.0j + 1.0, 1.0 + 2.0j)


class TestAssertApproxArray(unittest.TestCase):

    def test_all_elements_within_tolerance_passes(self):
        expected = np.array([1.0, 2.0, -3.0])
        assert_approx_array(expected + 1e-13, expected)

    def test_one_element_outside_tolerance_fails(self):
        expected = np.array([1.0, 2.0, -3.0])
        actual = expected.copy()
        actual[1] += 1.0
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected)

    def test_diff_equal_to_tolerance_passes_elementwise(self):
        expected = np.array([2.0 ** 10, 2.0 ** 10])
        tol = (2.0 ** -4) * expected  # 64.0, exact
        actual = expected + tol
        assert_approx_array(actual, expected, rel=2.0 ** -4)

    def test_diff_just_above_tolerance_fails_elementwise(self):
        expected = np.array([2.0 ** 10, 2.0 ** 10])
        tol = (2.0 ** -4) * expected
        actual = expected + tol
        actual[0] += 1.0  # push exactly one element past the tolerance
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected, rel=2.0 ** -4)

    def test_abs_only_no_relative_term_leaks_in(self):
        expected = np.array([1.0e9, 2.0e9])
        actual = expected + 1.0
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected, abs=1e-6)

    def test_nan_element_fails_even_if_others_match(self):
        expected = np.array([1.0, 2.0, 3.0])
        actual = np.array([1.0, float("nan"), 3.0])
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected)

    def test_matching_nan_position_still_fails(self):
        # NaN is never "equal" under pytest.approx semantics, even NaN vs NaN
        expected = np.array([1.0, float("nan")])
        actual = np.array([1.0, float("nan")])
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected)

    def test_same_sign_infinity_elementwise_passes(self):
        expected = np.array([1.0, float("inf")])
        actual = np.array([1.0, float("inf")])
        assert_approx_array(actual, expected)

    def test_opposite_sign_infinity_elementwise_fails(self):
        expected = np.array([1.0, float("inf")])
        actual = np.array([1.0, float("-inf")])
        with self.assertRaises(AssertionError):
            assert_approx_array(actual, expected)

    def test_broadcasts_scalar_expected_against_array_actual(self):
        assert_approx_array(np.array([2.0, 2.0, 2.0]), 2.0)
        with self.assertRaises(AssertionError):
            assert_approx_array(np.array([2.0, 2.5, 2.0]), 2.0)


class TestApproxTestCaseMixin(ApproxTestCase):
    """The unittest.TestCase-bound wrappers must fail via self.fail (i.e.
    self.failureException, normally AssertionError) rather than leaking a
    bare AssertionError from a different code path."""

    def test_assert_approx_passes_silently_within_tolerance(self):
        self.assertApprox(1.0 + 1e-13, 1.0)

    def test_assert_approx_fails_via_the_testcase_failure_path(self):
        with self.assertRaises(self.failureException):
            self.assertApprox(2.0, 1.0)

    def test_assert_approx_nan_fails_via_the_testcase_failure_path(self):
        with self.assertRaises(self.failureException):
            self.assertApprox(float("nan"), 1.0)

    def test_assert_approx_array_passes_silently_within_tolerance(self):
        self.assertApproxArray([1.0, 2.0], [1.0, 2.0 + 1e-13])

    def test_assert_approx_array_fails_via_the_testcase_failure_path(self):
        with self.assertRaises(self.failureException):
            self.assertApproxArray([1.0, 2.0], [1.0, 3.0])

    def test_custom_message_is_used_on_failure(self):
        try:
            self.assertApprox(2.0, 1.0, msg="custom failure text")
        except self.failureException as exc:
            self.assertIn("custom failure text", str(exc))
        else:
            self.fail("expected assertApprox to raise")


if __name__ == "__main__":
    unittest.main()
