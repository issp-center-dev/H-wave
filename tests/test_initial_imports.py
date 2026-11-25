"""The import boundary of ``hwave.initial``.

``hwave.initial`` reaches the solver only through
``hwave.solver.initial_facade``: no module of the package imports another
``hwave.solver`` module or ``hwave.qlms``/``hwave.qlmsio``. Later tasks add
modules to the package; this test scans every one of them.
"""
import ast
import os
import unittest

_PKG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                    "src", "hwave", "initial")

#: the only solver-side module ``hwave.initial`` may import
ALLOWED = "hwave.solver.initial_facade"


def _absolute(module, level, package):
    """The absolute name of an import of ``module`` at relative ``level``
    from ``package`` (``hwave.initial`` or a subpackage)."""
    if level == 0:
        return module or ""
    parts = package.split(".")
    base = parts[:len(parts) - (level - 1)]
    return ".".join(base + ([module] if module else []))


def imported_names(path, package="hwave.initial"):
    """Every module name ``path`` imports (``from X import y`` yields
    ``X.y``: ``y`` may be a submodule, and an attribute of a forbidden
    module is forbidden as well)."""
    with open(path, encoding="utf-8") as f:
        tree = ast.parse(f.read(), filename=path)
    names = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.extend(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            mod = _absolute(node.module, node.level, package)
            names.extend("{}.{}".format(mod, alias.name) for alias in node.names)
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
              and node.func.attr == "import_module" and node.args
              and isinstance(node.args[0], ast.Constant)):
            names.append(str(node.args[0].value))
        elif (isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
              and node.func.id == "__import__" and node.args
              and isinstance(node.args[0], ast.Constant)):
            names.append(str(node.args[0].value))
    return names


def _forbidden(name):
    if name == ALLOWED or name.startswith(ALLOWED + "."):
        return False
    if name in ("hwave.solver", "hwave.qlms", "hwave.qlmsio"):
        return True
    return name.startswith(("hwave.solver.", "hwave.qlms.", "hwave.qlmsio.")) \
        or name.startswith("hwave.qlms")


def _modules():
    for dirpath, _, files in os.walk(_PKG):
        for fn in sorted(files):
            if fn.endswith(".py"):
                yield os.path.join(dirpath, fn)


class TestImportBoundary(unittest.TestCase):
    def test_import_boundary(self):
        modules = list(_modules())
        self.assertTrue(modules)
        bad = []
        for path in modules:
            for name in imported_names(path):
                if _forbidden(name):
                    bad.append("{}: {}".format(os.path.relpath(path, _PKG), name))
        self.assertEqual(bad, [], "hwave.initial must import the solver only "
                         "through {}".format(ALLOWED))

    def test_scanner_sees_forbidden_imports(self):
        """The scanner flags each import form (absolute, relative, dotted
        ``from``, ``importlib``) and lets the facade through."""
        import tempfile
        cases = {
            "import hwave.solver.uhfk\n": True,
            "from hwave.solver import uhfk\n": True,
            "from ..solver import uhfk\n": True,
            "from .. import qlms\n": True,
            "from hwave import solver\n": True,
            "import hwave.qlmsio.read_input_k as r\n": True,
            "import importlib\nimportlib.import_module('hwave.solver.rpa')\n": True,
            "from hwave.solver import initial_facade\n": False,
            "from ..solver.initial_facade import load_problem\n": False,
            "import hwave.results\n": False,
            "from . import patterns\n": False,
        }
        d = tempfile.mkdtemp()
        try:
            for src, want in cases.items():
                with self.subTest(src=src):
                    p = os.path.join(d, "m.py")
                    with open(p, "w") as f:
                        f.write(src)
                    got = any(_forbidden(n) for n in imported_names(p))
                    self.assertEqual(got, want)
        finally:
            import shutil
            shutil.rmtree(d, ignore_errors=True)


if __name__ == "__main__":
    unittest.main()
