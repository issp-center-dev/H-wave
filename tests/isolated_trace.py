"""Run a tracemalloc measurement in a fresh Python process (issue #266).

A traced peak taken inside the test process counts every allocation the
measured code makes, including the ones it makes on behalf of state left
behind by EARLIER tests. The case that broke the bond-FLEX preflight
checks: a test that lowers a logger's level and does not restore it makes
the solver emit INFO/DEBUG records that it would otherwise drop, and
pytest's log capture handlers keep every such record (and its formatted
text) until the test ends. Those bytes land inside the measured phase, so
the measurement depended on which modules ran before it, and on the runner
(unittest keeps no records). A running coverage tracer adds its own
bookkeeping the same way.

:func:`run_isolated` calls ``module.function(*args)`` in a new interpreter
started with ``-B`` (no stale bytecode), from the repository root, with
only ``src`` and the repository root on ``PYTHONPATH`` and with the
coverage hand-over variables removed, so neither the earlier tests' state
nor a coverage tracer can reach the measurement. The function must return
a JSON-serializable value; the assertions stay in the calling test.
"""
import json
import os
import subprocess
import sys

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_MARK = "@@ISOLATED_TRACE@@ "

_SCRIPT = r"""
import importlib, json, sys
mod, fn, args = sys.argv[1], sys.argv[2], json.loads(sys.argv[3])
res = getattr(importlib.import_module(mod), fn)(*args)
sys.stdout.write("\n" + {mark!r} + json.dumps(res) + "\n")
""".format(mark=_MARK)


def _clean_env():
    env = {k: v for k, v in os.environ.items()
           if not (k.startswith("COV_CORE_") or k.startswith("COVERAGE_"))}
    env["PYTHONPATH"] = os.pathsep.join([os.path.join(_ROOT, "src"), _ROOT])
    env.setdefault("OPENBLAS_NUM_THREADS", "1")
    return env


def run_isolated(module, function, *args, timeout=600):
    """Return ``module.function(*args)`` evaluated in a fresh process."""
    proc = subprocess.run([sys.executable, "-B", "-c", _SCRIPT, module, function, json.dumps(args)],
                          cwd=_ROOT, env=_clean_env(), capture_output=True, text=True,
                          timeout=timeout)
    lines = [ln for ln in proc.stdout.splitlines() if ln.startswith(_MARK)]
    if proc.returncode != 0 or not lines:
        raise RuntimeError("isolated measurement {}.{} failed (exit {}):\n{}".format(
            module, function, proc.returncode, proc.stderr[-4000:]))
    return json.loads(lines[-1][len(_MARK):])
