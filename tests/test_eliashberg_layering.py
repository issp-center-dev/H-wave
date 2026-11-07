"""Layering guard for the sc.py split (dev/sc-split-plan.md, step 0).

``hwave.sc`` is being taken apart into ``hwave/solver/eliashberg/``, one
layer per pull request, each a pure move. Two things must hold throughout:

1. **The shell stays complete.** Every name that the tests or the solver
   modules reach through ``hwave.sc`` (``sc._load_chi0q``,
   ``sc._solve_leading``, ...) is still an attribute of ``hwave.sc`` after a
   move, because the shell re-exports it. The set of names is read from the
   source tree at test time, so a new reference is covered the moment it is
   written and a dropped re-export fails here instead of in some test that
   runs an hour later.

2. **The layer inversion only shrinks.** Modules under ``hwave/solver/`` are
   the layer ``hwave.sc`` is built on; the ones that import ``hwave.sc`` back
   invert that. ``SOLVER_MODULES_IMPORTING_SC`` is the exact set allowed
   today. A module that stops importing ``hwave.sc`` must be removed from the
   set (the test fails otherwise), so the set can only go down; a module that
   starts importing it fails the test.

Both checks parse the sources with ``ast`` rather than importing the
modules, so they run in milliseconds and in the fast gate.
"""
import ast
import os
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SRC = os.path.join(REPO, "src", "hwave")
SOLVER = os.path.join(SRC, "solver")

#: Solver-side modules that import ``hwave.sc`` today (the layer inversion
#: the split removes). Edit this set only to REMOVE an entry.
SOLVER_MODULES_IMPORTING_SC = frozenset({
    # step 1 (leaf layers) released _sc_matrices_myo, offsite and
    # second_order; step 4 (bond, spin-orbital solve) released so_bond;
    # step 6 released eliashberg_bond_io, eliashberg_dynamic,
    # eliashberg_dynamic_so, flex and so_pairing (they read the entry's
    # constants and functions from hwave.solver.eliashberg.static now).
    # The layer inversion is gone: nothing under hwave/solver/ imports
    # hwave.sc.
})

#: The public surface of ``hwave.sc`` that must outlive the split: the
#: console script (``hwave_sc = "hwave.sc:main"``), tsweep and external
#: scripts call these.
PUBLIC_API = ("calc_eliashberg", "main", "reject_spin_orbital_mode")

def _python_files(*roots):
    for root in roots:
        for dirpath, _dirs, files in os.walk(root):
            for f in files:
                if f.endswith(".py"):
                    yield os.path.join(dirpath, f)


def _module_imports_hwave_sc(path):
    """True if the module imports ``hwave.sc`` (any form, anywhere in the
    file, including lazy imports inside functions)."""
    tree = ast.parse(open(path, encoding="utf-8").read(), path)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(a.name == "hwave.sc" for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if node.module == "hwave.sc":
                return True
            if node.module == "hwave" and any(a.name == "sc" for a in node.names):
                return True
            # ``from .. import sc`` inside hwave/solver
            if node.level == 2 and node.module is None and any(a.name == "sc" for a in node.names):
                return True
        elif isinstance(node, ast.Call):
            # importlib.import_module("hwave.sc") / __import__("hwave.sc")
            if node.args and isinstance(node.args[0], ast.Constant) and node.args[0].value == "hwave.sc":
                return True
    return False


def _names_reached_through_sc(paths):
    """Every attribute name read off ``hwave.sc`` in ``paths``.

    Resolved from the AST, not from text: an alias is only counted once
    the file binds it to ``hwave.sc`` (``import hwave.sc as sc``,
    ``from hwave import sc``, ``import hwave.sc as sc_mod``, ...), so a
    local array called ``sc`` in a test, and ``sc._x`` inside a comment or
    docstring, do not register. ``from hwave.sc import a, b`` counts ``a``
    and ``b`` directly.
    """
    names = set()
    for path in paths:
        tree = ast.parse(open(path, encoding="utf-8").read(), path)
        aliases = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for a in node.names:
                    if a.name == "hwave.sc" and a.asname:
                        aliases.add(a.asname)
            elif isinstance(node, ast.ImportFrom):
                if node.module == "hwave.sc":
                    names.update(a.name for a in node.names)
                elif node.module == "hwave":
                    for a in node.names:
                        if a.name == "sc":
                            aliases.add(a.asname or "sc")
        # string patch targets need no alias: mock.patch("hwave.sc.x")
        for node in ast.walk(tree):
            if isinstance(node, ast.Call) and node.args and isinstance(node.args[0], ast.Constant):
                v = node.args[0].value
                if isinstance(v, str) and v.startswith("hwave.sc."):
                    names.add(v[len("hwave.sc."):].split(".")[0])
        if not aliases:
            continue
        for node in ast.walk(tree):
            if (isinstance(node, ast.Attribute)
                    and isinstance(node.value, ast.Name)
                    and node.value.id in aliases):
                names.add(node.attr)
            # getattr(sc, "x") / setattr(sc, "x", v) / mock.patch.object(sc, "x")
            # / monkeypatch.setattr(sc, "x", v): the alias is the first argument
            # and the name a string constant in the second
            elif (isinstance(node, ast.Call) and len(node.args) >= 2
                    and isinstance(node.args[0], ast.Name)
                    and node.args[0].id in aliases
                    and isinstance(node.args[1], ast.Constant)
                    and isinstance(node.args[1].value, str)):
                names.add(node.args[1].value)
    return names


class TestShellStaysComplete(unittest.TestCase):

    def test_every_referenced_name_is_an_attribute_of_hwave_sc(self):
        import hwave.sc as sc

        paths = [p for p in _python_files(HERE, SRC)
                 if os.path.abspath(p) != os.path.abspath(sc.__file__)
                 and os.path.basename(p) != os.path.basename(__file__)]
        referenced = _names_reached_through_sc(paths)
        self.assertGreater(len(referenced), 100,
                           "the scan found suspiciously few references; "
                           "has the import form of the tests changed?")
        missing = sorted(n for n in referenced if not hasattr(sc, n))
        self.assertEqual(missing, [],
                         "names reached through hwave.sc that hwave.sc no longer "
                         "exposes (re-export them from the module they moved to): "
                         + ", ".join(missing))

    def test_public_api_present(self):
        import hwave.sc as sc
        for name in PUBLIC_API:
            self.assertTrue(callable(getattr(sc, name, None)), name)


class TestLayerInversionOnlyShrinks(unittest.TestCase):

    def test_solver_modules_importing_sc_are_exactly_the_allowed_set(self):
        importing = set()
        for path in _python_files(SOLVER):
            mod = os.path.splitext(os.path.basename(path))[0]
            if mod == "__init__":
                continue
            # the split's own package is allowed to be imported BY sc, never
            # to import it: it is checked like any other solver module
            if _module_imports_hwave_sc(path):
                importing.add(mod)
        new = sorted(importing - SOLVER_MODULES_IMPORTING_SC)
        self.assertEqual(new, [],
                         "solver modules that newly import hwave.sc (the solver "
                         "layer must not depend on the entry module): " + ", ".join(new))
        gone = sorted(SOLVER_MODULES_IMPORTING_SC - importing)
        self.assertEqual(gone, [],
                         "these modules no longer import hwave.sc; remove them from "
                         "SOLVER_MODULES_IMPORTING_SC so the set keeps shrinking: "
                         + ", ".join(gone))


if __name__ == "__main__":
    unittest.main()
