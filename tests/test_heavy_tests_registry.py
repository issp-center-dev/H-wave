"""Consistency gate for ``tests/heavy_tests.py``.

The registry is only worth having if it cannot silently drift from the
code, so this module FAILS the FAST gate whenever the two disagree:

* an entry naming a module/class/method that does not exist
  (:meth:`TestHeavyRegistryConsistency.test_every_registry_entry_resolves`)
* a test carrying ``@heavy`` that is absent from the registry
  (:meth:`...test_no_decorated_test_is_missing_from_the_registry`)
* an entry whose test is NOT actually skipped when ``HWAVE_FULL_TESTS``
  is unset (:meth:`...test_registered_tests_are_skipped_in_the_fast_gate`)
* one of the documented fast-gate exceptions having been decorated after
  all (:meth:`...test_documented_exceptions_stay_in_the_fast_gate`)

Everything here works by INSPECTING the loaded test classes -- importing
the named modules and reading ``__unittest_skip__`` /
``__hwave_heavy__`` off the resolved attributes, plus one
``unittest.defaultTestLoader.discover`` walk for the "decorated but
unregistered" direction. Nothing parses source text, so a decorator
applied through an alias, a helper, or a class body still counts.

Discovery is collection-only (``discover`` instantiates ``TestCase``
objects without running them), the same cheap pattern
``test_bubble_kernel.TestGateCollection`` already uses; it must be run
from the repository root, which is this project's standing convention for
the whole suite.

This module is NEVER decorated ``@heavy``: it is precisely the check the
fast gate exists to run.
"""

import importlib
import unittest

from tests import heavy_tests


def _resolve(mod_name, cls_name, meth_name):
    """Return the class-level attribute for ``mod.Class.method``.

    Returns ``(obj, None)`` on success and ``(None, why)`` with a
    human-readable reason when any leg of the path is missing.
    """
    try:
        module = importlib.import_module("tests." + mod_name)
    except ImportError as exc:
        return None, "module tests.{} does not import: {}".format(
            mod_name, exc)
    cls = getattr(module, cls_name, None)
    if cls is None:
        return None, "class {} not found in tests.{}".format(
            cls_name, mod_name)
    if not (isinstance(cls, type) and issubclass(cls, unittest.TestCase)):
        return None, "tests.{}.{} is not a unittest.TestCase".format(
            mod_name, cls_name)
    meth = getattr(cls, meth_name, None)
    if meth is None:
        return None, "method {} not found on tests.{}.{}".format(
            meth_name, mod_name, cls_name)
    if meth_name not in unittest.defaultTestLoader.getTestCaseNames(cls):
        return None, (
            "tests.{}.{}.{} exists but is not collected as a test by "
            "unittest's loader".format(mod_name, cls_name, meth_name))
    return meth, None


def _discovered_heavy_ids():
    """``"module.Class.method"`` for every collected test carrying the
    ``@heavy`` marker, found by walking the discovered suite."""
    suite = unittest.defaultTestLoader.discover("tests", pattern="test_*.py")
    found = set()
    stack = [suite]
    while stack:
        item = stack.pop()
        if isinstance(item, unittest.TestSuite):
            stack.extend(item)
            continue
        if not isinstance(item, unittest.TestCase):
            continue
        cls = type(item)
        name = getattr(item, "_testMethodName", None)
        if name is None:
            continue
        meth = getattr(cls, name, None)
        if meth is None or not getattr(
                meth, heavy_tests.HEAVY_MARKER_ATTR, False):
            continue
        # Strip the "tests." prefix so the id matches the registry's own
        # bare-module-name form under both invocation styles (package
        # relative discovery gives "tests.test_x", a tests/-on-sys.path
        # run gives "test_x").
        mod = cls.__module__
        if mod.startswith("tests."):
            mod = mod[len("tests."):]
        found.add("{}.{}.{}".format(mod, cls.__qualname__, name))
    return found


class TestHeavyRegistryConsistency(unittest.TestCase):
    """The registry and the decorated code must agree, both directions."""

    def test_registry_is_well_formed(self):
        """Rows are 4-tuples of non-empty strings, unique, and sorted."""
        ids = []
        for row in heavy_tests.HEAVY_TESTS:
            self.assertEqual(
                len(row), 4,
                "HEAVY_TESTS rows are (module, class, method, reason): "
                "{!r}".format(row))
            for field in row:
                self.assertIsInstance(field, str, row)
                self.assertTrue(field.strip(), "empty field in {!r}".format(row))
            ids.append(row[:3])
        self.assertEqual(
            len(ids), len(set(ids)),
            "duplicate entries in HEAVY_TESTS: {}".format(
                sorted(i for i in set(ids) if ids.count(i) > 1)))
        self.assertEqual(
            ids, sorted(ids),
            "HEAVY_TESTS must stay sorted by (module, class, method) so "
            "additions produce readable diffs")

    def test_every_registry_entry_resolves(self):
        """No entry may name a test that does not exist."""
        problems = []
        for mod, cls, meth, _reason in heavy_tests.HEAVY_TESTS:
            _obj, why = _resolve(mod, cls, meth)
            if why is not None:
                problems.append(why)
        self.assertEqual(
            problems, [],
            "tests/heavy_tests.py names test(s) that do not exist -- "
            "delete or rename the entry:\n  " + "\n  ".join(problems))

    def test_every_registry_entry_carries_the_decorator(self):
        """An entry whose method lacks ``@heavy`` is a registry lie."""
        undecorated = []
        for mod, cls, meth, _reason in heavy_tests.HEAVY_TESTS:
            obj, why = _resolve(mod, cls, meth)
            if why is not None:
                continue  # reported by test_every_registry_entry_resolves
            if not getattr(obj, heavy_tests.HEAVY_MARKER_ATTR, False):
                undecorated.append("{}.{}.{}".format(mod, cls, meth))
        self.assertEqual(
            undecorated, [],
            "listed in tests/heavy_tests.py but NOT decorated @heavy (so "
            "they still run in the fast gate): {}".format(undecorated))

    def test_no_decorated_test_is_missing_from_the_registry(self):
        """A ``@heavy`` test absent from the registry is undocumented."""
        declared = heavy_tests.heavy_test_ids()
        discovered = _discovered_heavy_ids()
        unregistered = sorted(discovered - declared)
        self.assertEqual(
            unregistered, [],
            "decorated @heavy but missing from HEAVY_TESTS in "
            "tests/heavy_tests.py -- add a row (with its reason): "
            "{}".format(unregistered))
        # The other direction is a stronger statement than
        # test_every_registry_entry_carries_the_decorator: it also catches
        # an entry the loader never collects at all.
        undiscovered = sorted(declared - discovered)
        self.assertEqual(
            undiscovered, [],
            "listed in HEAVY_TESTS but not discovered as a @heavy test by "
            "`unittest discover -s tests`: {}".format(undiscovered))

    def test_registered_tests_are_skipped_in_the_fast_gate(self):
        """With ``HWAVE_FULL_TESTS`` unset, every entry must be skipped;
        with it set, none of them may be."""
        full = heavy_tests.full_suite_enabled()
        wrong = []
        for mod, cls, meth, _reason in heavy_tests.HEAVY_TESTS:
            obj, why = _resolve(mod, cls, meth)
            if why is not None:
                continue  # reported by test_every_registry_entry_resolves
            skipped = bool(getattr(obj, "__unittest_skip__", False))
            if skipped == full:
                wrong.append("{}.{}.{} (skipped={})".format(
                    mod, cls, meth, skipped))
        self.assertEqual(
            wrong, [],
            "with {}={!r} (full_suite_enabled={}) these registered heavy "
            "tests have the WRONG skip state: {}".format(
                heavy_tests.FULL_SUITE_ENV,
                __import__("os").environ.get(heavy_tests.FULL_SUITE_ENV),
                full, wrong))
        if not full:
            # And the skip must be the one this registry installed, not an
            # unrelated skipUnless that happens to sit on the same method.
            for mod, cls, meth, _reason in heavy_tests.HEAVY_TESTS:
                obj, why = _resolve(mod, cls, meth)
                if why is not None:
                    continue
                self.assertEqual(
                    getattr(obj, "__unittest_skip_why__", None),
                    heavy_tests.SKIP_REASON,
                    "{}.{}.{} is skipped for a reason other than the heavy "
                    "gate".format(mod, cls, meth))

    def test_documented_exceptions_stay_in_the_fast_gate(self):
        """The four >=5 s tests kept on purpose must never be decorated."""
        problems = []
        for mod, cls, meth in heavy_tests.FAST_GATE_EXCEPTIONS:
            obj, why = _resolve(mod, cls, meth)
            if why is not None:
                problems.append(why)
                continue
            if getattr(obj, heavy_tests.HEAVY_MARKER_ATTR, False):
                problems.append(
                    "{}.{}.{} is a DOCUMENTED fast-gate exception (see the "
                    "heavy_tests module docstring) but carries @heavy".format(
                        mod, cls, meth))
        self.assertEqual(
            problems, [],
            "fast-gate exception list is out of sync:\n  "
            + "\n  ".join(problems))

    def test_exception_list_and_registry_are_disjoint(self):
        """A test cannot be both opt-in and a documented exception."""
        exceptions = frozenset(
            "{}.{}.{}".format(*row)
            for row in heavy_tests.FAST_GATE_EXCEPTIONS)
        overlap = sorted(exceptions & heavy_tests.heavy_test_ids())
        self.assertEqual(
            overlap, [],
            "listed both in HEAVY_TESTS and in FAST_GATE_EXCEPTIONS: "
            "{}".format(overlap))


class TestFullSuiteEnvParsing(unittest.TestCase):
    """``full_suite_enabled`` reads the environment the documented way."""

    def _with_env(self, value):
        import os
        old = os.environ.get(heavy_tests.FULL_SUITE_ENV)
        if value is None:
            os.environ.pop(heavy_tests.FULL_SUITE_ENV, None)
        else:
            os.environ[heavy_tests.FULL_SUITE_ENV] = value
        try:
            return heavy_tests.full_suite_enabled()
        finally:
            if old is None:
                os.environ.pop(heavy_tests.FULL_SUITE_ENV, None)
            else:
                os.environ[heavy_tests.FULL_SUITE_ENV] = old

    def test_unset_is_off(self):
        self.assertFalse(self._with_env(None))

    def test_falsey_spellings_are_off(self):
        for value in ("", "0", "false", "FALSE", "No", "off", "  0  "):
            self.assertFalse(self._with_env(value), value)

    def test_truthy_spellings_are_on(self):
        for value in ("1", "true", "TRUE", "yes", "on", "anything"):
            self.assertTrue(self._with_env(value), value)
