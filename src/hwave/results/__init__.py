"""Read the outputs of H-wave runs.

Open a run from its input TOML file with :func:`open_run`; ``Run.path(key)``
gives the file the solver wrote for an output key. :func:`read_energy` reads
``energy.dat``; :func:`read_scf_status` gives the convergence status.
"""
from hwave.results.errors import (
    AmbiguousRunError, MissingConfigError, MissingOutputError, ResultsError,
    RunRecordError, UnsupportedLayoutError)
from hwave.results.run import (
    Run, as_run, open_run, read_run_record, run_record_path, write_run_record)
from hwave.results.status import ScfStatus, read_scf_status
from hwave.results.uhf import Energy, read_energy

#: the version of the readers' interpretation of the output files
READER_VERSION = "1"

__all__ = [
    "READER_VERSION",
    "Run", "open_run", "as_run", "read_run_record", "write_run_record",
    "run_record_path",
    "ResultsError", "UnsupportedLayoutError", "MissingOutputError",
    "MissingConfigError", "AmbiguousRunError", "RunRecordError",
    "Energy", "read_energy", "ScfStatus", "read_scf_status",
]
