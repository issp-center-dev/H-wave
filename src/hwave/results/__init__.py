"""Read the outputs of H-wave runs.

Open a run from its input TOML file with :func:`open_run`; ``Run.path(key)``
gives the file the solver wrote for an output key. :func:`read_energy` reads
``energy.dat``; :func:`read_scf_status` gives the convergence status.
:func:`read_uhfr_eigen`, :func:`read_uhfr_green` and
:func:`read_uhfr_order_parameters` read the UHFr eigenstates, Green function
and order parameters. :func:`read_bands`, :func:`read_uhfk_green` and
:func:`read_uhfk_order_parameters` read the UHFk bands, Green function and
order parameters.
:func:`order_parameters_from_density` gives the spin and charge order
parameters of a one-cell density.
"""
from hwave.results.errors import (
    AmbiguousRunError, MissingConfigError, MissingOutputError, ResultsError,
    RunRecordError, UnsupportedLayoutError)
from hwave.results.order import (
    OrderParameters, cosine_amplitude, order_parameters_from_density,
    sum_rule_tolerance)
from hwave.results.run import (
    Run, as_run, open_run, read_run_record, run_record_path, write_run_record)
from hwave.results.status import ScfStatus, read_scf_status
from hwave.results.uhf import (
    Bands, Energy, GreenR, OneBodyG, UHFrEigen, UHFrGreen, read_bands,
    read_energy, read_onebodyg, read_uhfk_green, read_uhfk_order_parameters,
    read_uhfr_eigen, read_uhfr_green, read_uhfr_order_parameters,
    so_density_from_green)

#: the version of the readers' interpretation of the output files
READER_VERSION = "1"

__all__ = [
    "READER_VERSION",
    "Run", "open_run", "as_run", "read_run_record", "write_run_record",
    "run_record_path",
    "ResultsError", "UnsupportedLayoutError", "MissingOutputError",
    "MissingConfigError", "AmbiguousRunError", "RunRecordError",
    "Energy", "read_energy", "ScfStatus", "read_scf_status",
    "UHFrEigen", "read_uhfr_eigen", "UHFrGreen", "read_uhfr_green",
    "read_uhfr_order_parameters",
    "Bands", "read_bands", "GreenR", "read_uhfk_green",
    "read_uhfk_order_parameters", "OneBodyG", "read_onebodyg",
    "so_density_from_green",
    "OrderParameters", "order_parameters_from_density", "sum_rule_tolerance",
    "cosine_amplitude",
]
