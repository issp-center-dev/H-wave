"""Layering guard for the sc.py split (dev/sc-split-plan.md, steps 0 and 7).

``hwave.sc`` was taken apart into ``hwave/solver/eliashberg/``, one layer per
pull request, each a pure move; it is now the entry shell. Two things must
hold from here on:

1. **The shell is the public API and nothing else.** The only names reached
   through ``hwave.sc`` anywhere in the tree are :data:`PUBLIC_API` and
   :data:`PUBLIC_CONSTANTS` (the console script, the temperature sweep and
   external scripts use them); a private name is read from the module that
   owns it. The set of names is read from the source tree at test time, so a
   new ``sc._x`` reference fails here instead of waiting for the shell to
   lose it.

2. **The layer inversion stays gone.** Modules under ``hwave/solver/`` are
   the layer ``hwave.sc`` is built on; none of them imports ``hwave.sc``.
   ``SOLVER_MODULES_IMPORTING_SC`` is the exact set allowed today: empty. A
   module that starts importing ``hwave.sc`` fails the test.

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

#: The public surface of ``hwave.sc`` that outlives the split: the console
#: script (``hwave_sc = "hwave.sc:main"``), tsweep and external scripts call
#: these.
PUBLIC_API = ("calc_eliashberg", "main", "reject_spin_orbital_mode")
PUBLIC_CONSTANTS = ("DEFAULT_SOLVER_MODE", "IN_PROCESS_DEFAULT_SOLVER_MODE",
                    "EXIT_EIGEN_NOT_CONVERGED")

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
            # ``from .. import sc`` / ``from ..sc import x`` inside hwave/solver,
            # and the three-dot forms from inside hwave/solver/eliashberg
            if node.level >= 2 and (
                    (node.module is None and any(a.name == "sc" for a in node.names))
                    or node.module == "sc"):
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


class TestShellIsThePublicApi(unittest.TestCase):

    def test_only_the_public_api_is_reached_through_hwave_sc(self):
        import hwave.sc as sc

        paths = [p for p in _python_files(HERE, SRC)
                 if os.path.abspath(p) != os.path.abspath(sc.__file__)
                 and os.path.basename(p) != os.path.basename(__file__)]
        referenced = _names_reached_through_sc(paths)
        self.assertIn("calc_eliashberg", referenced,
                      "the scan found no reference to hwave.sc.calc_eliashberg; "
                      "has the import form of the tests changed?")
        allowed = set(PUBLIC_API) | set(PUBLIC_CONSTANTS)
        private = sorted(n for n in referenced
                         if n not in allowed and not n.startswith("__"))
        self.assertEqual(private, [],
                         "private names reached through hwave.sc (read them from "
                         "the hwave.solver.eliashberg module that owns them): "
                         + ", ".join(private))

    def test_public_api_present(self):
        import hwave.sc as sc
        for name in PUBLIC_API:
            self.assertTrue(callable(getattr(sc, name, None)), name)
        for name in PUBLIC_CONSTANTS:
            self.assertTrue(hasattr(sc, name), name)
        self.assertEqual(sorted(sc.__all__), sorted(PUBLIC_API + PUBLIC_CONSTANTS))


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
