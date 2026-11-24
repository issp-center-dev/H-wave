"""Reference snapshot of the NON-spin-orbital RPA and FLEX archives.

The spin-orbital producers (RPA, Hartree-Fock FLEX) write fingerprint and
state members into their archives; a non-spin-orbital run must keep exactly
the member set and values it had before. This module holds the two small
non-spin-orbital runs that ``tests/test_so_provenance.py`` repeats, and,
run as a script, writes their archives into this directory::

    PYTHONPATH=<checkout>/src:<this repo> python3 -B tests/so_provenance_ref/make.py

It was run ONCE, on the tree before the spin-orbital provenance members
were added, to produce the committed ``<run>/<archive>.npz`` files and
``platform.txt`` (the platform the snapshot was made on).
"""
import contextlib
import logging
import os
import platform
import shutil
import sys
import tempfile

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPIN_DIR = os.path.join(os.path.dirname(_HERE), "equivalence_input", "spin")

#: (archive key, file name) of each run
RPA_OUTPUTS = {"chi0q": "chi0q.npz", "chiq": "chiq.npz"}
FLEX_OUTPUTS = {"chi0q": "chi0q.npz", "chiq": "chiq.npz", "green": "green.npz",
                "sigma": "sigma.npz"}


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _run(mode, param, files, outputs, extra_mode, out_dir):
    import hwave.qlms as qlms
    d = tempfile.mkdtemp(prefix="so_prov_ref_")
    try:
        for f in files.values():
            shutil.copy(os.path.join(_SPIN_DIR, f), os.path.join(d, f))
        mode_tbl = {"mode": mode, "enable_spin_orbital": False, "param": dict(param)}
        mode_tbl.update(extra_mode)
        inter = dict(files, path_to_input=d)
        inp = {"log": {"print_level": 0}, "mode": mode_tbl,
               "file": {"input": {"path_to_input": d, "interaction": inter},
                        "output": dict(outputs, path_to_output=out_dir)}}
        with _quiet():
            qlms.run(input_dict=inp)
    finally:
        shutil.rmtree(d, ignore_errors=True)
    return {k: os.path.join(out_dir, v) for k, v in outputs.items()}


def run_rpa_non_so(out_dir):
    """Non-spin-orbital RPA (general, ring) with a Zeeman Extern (spin-diag):
    two orbitals on a 2x2 lattice."""
    param = {"T": 1.0, "filling": 0.4, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1],
             "Nmat": 8, "coeff_extern": 1.0}
    files = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
             "Extern": "extern_zeeman.dat", "CoulombIntra": "coulombintra.dat"}
    return _run("RPA", param, files, RPA_OUTPUTS,
                {"calc_scheme": "general", "calc_type": "ring"}, out_dir)


def run_flex_non_so(out_dir):
    """Non-spin-orbital Hartree-Fock FLEX (general), three maps."""
    param = {"T": 1.0, "filling": 0.4, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1],
             "Nmat": 8, "IterationMax": 3, "Mix": 0.5, "mixing_scheme": "linear",
             "flex_hartree_fock": True}
    files = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombIntra": "coulombintra.dat"}
    return _run("FLEX", param, files, FLEX_OUTPUTS, {"calc_scheme": "general"}, out_dir)


RUNS = {"_run_rpa_non_so": run_rpa_non_so, "_run_flex_non_so": run_flex_non_so}


def platform_tag():
    return "{} {} numpy {}".format(platform.system(), platform.machine(), np.__version__)


def main(dest=_HERE):
    for name, run in RUNS.items():
        out = os.path.join(dest, name.lstrip("_"))
        os.makedirs(out, exist_ok=True)
        run(out)
    with open(os.path.join(dest, "platform.txt"), "w") as f:
        f.write(platform_tag() + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
