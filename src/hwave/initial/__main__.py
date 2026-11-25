"""``python -m hwave.initial``: the commands prepare, check and collect.

The process ends with ``os._exit`` after flushing its streams: the solver's
profiling table, printed when the interpreter shuts down, is not part of
the output of these commands (the solver is left unchanged).
"""
import os
import sys

from hwave.initial.cli import main

if __name__ == "__main__":
    code = main()
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(code)
