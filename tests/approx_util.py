"""Shared pytest.approx-equivalent assertions for the unittest-based H-wave
test suite, plus a caplog-equivalent log capture context manager.

These exist because several test modules were converted from pytest
(``assert x == pytest.approx(y)``, ``caplog.at_level(...)``) to plain
``unittest``. Reproducing ``pytest.approx``'s numeric semantics EXACTLY
matters -- silently loosening (or tightening) a tolerance changes what a
regression test actually catches.

``pytest.approx`` semantics (what this module reproduces):

    tolerance = max(rel * |expected|, abs)
        rel defaults to 1e-6, abs defaults to 1e-12
    rel given, abs omitted   -> tolerance = max(rel * |expected|, 1e-12)
    abs given, rel omitted   -> tolerance = EXACTLY abs (no relative term
                                 mixed in -- this is the case most easily
                                 gotten wrong)
    both omitted              -> tolerance = max(1e-6 * |expected|, 1e-12)
    actual == expected        -> PASS unconditionally (the equality
                                 short-circuit; covers exact equality and
                                 same-sign infinities, e.g. inf == inf)
    NaN anywhere               -> FAIL (``pytest.approx`` has no
                                 ``nan_ok=True`` here, so a NaN actual or
                                 expected value never compares equal)

Why ``numpy.testing.assert_allclose`` is NOT an equivalent substitute:

  * ``assert_allclose`` uses the SUM rule: ``|actual - expected| <= atol +
    rtol * |expected|``. ``pytest.approx`` uses the MAX rule: ``|actual -
    expected| <= max(rel * |expected|, abs)``. For the same ``(rtol, atol)``
    / ``(rel, abs)`` pair the sum rule can be up to 2x LOOSER, so a bare
    ``assert_allclose(rtol=1e-6, atol=1e-12)`` swapped in for
    ``pytest.approx(x)`` silently accepts deviations the original would
    have rejected.
  * ``assert_allclose`` defaults to ``equal_nan=True`` (a NaN actual
    compares equal to a NaN expected, so the assertion PASSES);
    ``pytest.approx`` always FAILS when either side is NaN. A NaN actual
    value -- almost always a sign of a real bug upstream -- can slip past
    ``assert_allclose`` silently where ``pytest.approx`` would have caught
    it.

Use ``assertApprox``/``assertApproxArray`` (via the ``ApproxTestCase``
mixin) or the standalone ``assert_approx``/``assert_approx_array``
functions (for call sites that are not ``unittest.TestCase`` methods)
wherever a test used to rely on ``pytest.approx``. Sites whose ORIGINAL
pytest form already used ``np.testing.assert_allclose`` directly (not
``pytest.approx``) are out of scope -- leave those as they are.
"""
import logging
import math
import operator
import unittest
from contextlib import contextmanager

import numpy as np

_DEFAULT_REL = 1e-6
_DEFAULT_ABS = 1e-12


def _scalar_tolerance(expected, rel, abs_):
    atol = _DEFAULT_ABS if abs_ is None else abs_
    if rel is None and abs_ is not None:
        # abs-only: the tolerance is EXACTLY abs, no relative term mixed in.
        return atol
    relative = (_DEFAULT_REL if rel is None else rel) * operator.abs(expected)
    return max(relative, atol)


def assert_approx(actual, expected, rel=None, abs=None, msg=None):
    """Standalone ``pytest.approx``-equivalent scalar assertion.

    Raises ``AssertionError`` on failure; use this at call sites that are
    not ``unittest.TestCase`` methods (e.g. plain helper functions). See
    the module docstring for the exact tolerance semantics.
    """
    if actual == expected:
        # exact equality short-circuit; also covers inf == inf (same sign)
        return
    diff = operator.abs(actual - expected)
    if math.isinf(operator.abs(expected)):
        # An infinite `expected` that didn't already match via the equality
        # short-circuit above can never be "close" by a relative/absolute
        # tolerance -- the relative term would itself be rel * inf == inf,
        # which would make ANY diff (even inf) compare <= tol. Reject
        # unconditionally instead, matching pytest.approx's own explicit
        # infinite-expected special case.
        raise AssertionError(
            msg or "{!r} != {!r}: expected is infinite and not exactly equal "
                   "to actual".format(actual, expected))
    tol = _scalar_tolerance(expected, rel, abs)
    if math.isnan(diff) or diff > tol:
        # a NaN diff means a NaN actual/expected, or unequal infinities
        # (inf - inf = nan) -- both must FAIL, never silently pass.
        raise AssertionError(
            msg or "{!r} != {!r} within tolerance (diff={!r}, tol={!r})".format(
                actual, expected, diff, tol))


def assert_approx_array(actual, expected, rel=None, abs=None, msg=None):
    """Standalone ``pytest.approx``-equivalent elementwise array assertion.

    Applies the scalar max-rule tolerance elementwise (broadcasting
    ``actual``/``expected`` against each other first) with the same
    equality short-circuit and NaN-fails behavior as ``assert_approx``.
    """
    actual_b, expected_b = np.broadcast_arrays(np.asarray(actual), np.asarray(expected))
    with np.errstate(invalid="ignore", over="ignore"):
        diff = np.abs(actual_b - expected_b)
        atol = _DEFAULT_ABS if abs is None else abs
        if rel is None and abs is not None:
            tol = np.full(diff.shape, atol, dtype=float)
        else:
            relative = (_DEFAULT_REL if rel is None else rel) * np.abs(expected_b)
            tol = np.broadcast_to(np.maximum(relative, atol), diff.shape)
        equal_mask = (actual_b == expected_b)
        # an infinite `expected` that isn't already an exact match (equal_mask)
        # can never be "close" by a tolerance -- see assert_approx's comment.
        inf_expected_mask = np.isinf(np.abs(expected_b))
        nan_mask = np.isnan(diff)
        ok = equal_mask | (~inf_expected_mask & ~nan_mask & (diff <= tol))
    if not np.all(ok):
        bad_idx = tuple(np.argwhere(~ok)[0])
        raise AssertionError(
            msg or "arrays not approx-equal ({} of {} element(s) differ); first "
                   "mismatch at index {}: actual={!r} expected={!r} diff={!r} "
                   "tol={!r}".format(int(np.count_nonzero(~ok)), ok.size, bad_idx,
                                      actual_b[bad_idx], expected_b[bad_idx],
                                      diff[bad_idx], tol[bad_idx]))


class ApproxTestCase(unittest.TestCase):
    """Common ``unittest.TestCase`` mixin adding ``assertApprox`` (scalar)
    and ``assertApproxArray`` (elementwise), both replicating
    ``pytest.approx`` semantics exactly. See the module docstring for the
    tolerance table and why ``np.testing.assert_allclose`` is not an
    equivalent substitute.
    """

    def assertApprox(self, actual, expected, rel=None, abs=None, msg=None):
        try:
            assert_approx(actual, expected, rel=rel, abs=abs, msg=msg)
        except AssertionError as exc:
            self.fail(str(exc))

    def assertApproxArray(self, actual, expected, rel=None, abs=None, msg=None):
        try:
            assert_approx_array(actual, expected, rel=rel, abs=abs, msg=msg)
        except AssertionError as exc:
            self.fail(str(exc))


@contextmanager
def _caplog(level, logger=None):
    """Minimal log-capture context manager, ``caplog.at_level``-equivalent.

    Attaches a handler directly to the named logger (root, matching the
    conventional default, when ``logger`` is None) for the duration of the
    block and yields the list of captured ``LogRecord``\\ s -- each carrying
    a populated ``.message`` attribute, same as the records list a
    log-capture fixture would hand back.

    Both the logger's level AND the handler's own level are set: without
    ``handler.setLevel``, a child logger with its own lower (more verbose)
    level can still propagate records past this handler that
    ``caplog.at_level`` would have filtered out -- propagated records are
    filtered by each handler's OWN level along the way, not by the
    ancestor logger's ``.level`` attribute.
    """
    target = logging.getLogger(logger)
    records = []

    class _Handler(logging.Handler):
        def emit(self, record):
            record.message = record.getMessage()
            records.append(record)

    handler = _Handler()
    handler.setLevel(level)
    old_level = target.level
    target.addHandler(handler)
    target.setLevel(level)
    try:
        yield records
    finally:
        target.removeHandler(handler)
        target.setLevel(old_level)
