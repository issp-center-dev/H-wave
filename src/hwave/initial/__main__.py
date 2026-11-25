"""``python -m hwave.initial``: the commands prepare, check and collect.

The solver's profiling table, printed when the interpreter shuts down, is
not part of the output of these commands: its statistics are dropped
before the normal exit (the solver is left unchanged).
"""
import sys

from hwave.initial.cli import main
from hwave.solver.initial_facade import discard_profile_statistics

if __name__ == "__main__":
    code = main()
    discard_profile_statistics()
    sys.exit(code)
