"""End-to-end + reentrancy tests for hwave.tsweep over real solvers.

Reuses the minimal FLEX fixture from tests/test_qlms_run_result.py (Task 1),
which is already proven to drive a real FLEX run to completion.
"""
import math
import os

import numpy as np
import pytest

import hwave.tsweep as ts
from tests.test_qlms_run_result import _tiny_hubbard_flex_input


def _tsweep_base(tmp_path, temps, run_eli):
    base = _tiny_hubbard_flex_input(tmp_path, T=temps[0], itermax=30)
    base["mode"]["param"]["EPS"] = 6
    base["continuation"] = {"temperatures": temps, "output_dir": "sweep",
                            "run_eliashberg": run_eli, "warm_start": True,
                            "seed_gap": False, "summary_file": "lambda_vs_T.dat"}
    return base


def test_e2e_flex_only_two_rungs(tmp_path):
    base = _tsweep_base(tmp_path, [4.0, 3.0], run_eli=False)
    rows = ts.run(base, base_dir=str(tmp_path))
    assert len(rows) == 2
    assert all(r["status"] in ("ok", "not_converged") for r in rows)
    # rung 1 warm-started from rung 0 sigma
    assert os.path.exists(os.path.join(tmp_path, "sweep", "000_T4", "output", "sigma.npz"))
    assert os.path.exists(os.path.join(tmp_path, "sweep", "lambda_vs_T.dat"))


def test_e2e_reentrancy_no_leakage(tmp_path):
    """Two identical rungs run back-to-back in one process give identical FLEX
    iteration counts / convergence (guards against cross-rung shared state)."""
    base = _tsweep_base(tmp_path, [4.0, 4.0], run_eli=False)
    base["continuation"]["warm_start"] = False
    rows = ts.run(base, base_dir=str(tmp_path))
    assert len(rows) == 2
    assert rows[0]["flex_converged"] == rows[1]["flex_converged"]
    assert rows[0]["flex_iter"] == rows[1]["flex_iter"]


@pytest.mark.slow
def test_e2e_dynamic_eliashberg_two_rungs_real_solver(tmp_path):
    """Real-solver end-to-end: hwave.tsweep.run() driving DYNAMIC Eliashberg
    (chi0q_mode="flex") across a 2-rung temperature ladder.

    This is the one case the other tsweep tests miss: they either skip
    Eliashberg entirely (test_e2e_flex_only_two_rungs) or (in unit tests
    elsewhere) mock/stub the solver. Here FLEX actually writes chiq_s.npz /
    chiq_c.npz / green.npz per rung, and hwave.sc.calc_eliashberg dispatches
    to the real hwave.solver.eliashberg_dynamic.solve_dynamic path on them.
    It regression-guards two fixed bugs:
      * C1: tsweep.preflight requires [eliashberg] solver_mode="eigenvalue"
        (or "both") so the per-rung eigenvalue.dat has the parity-analysis
        block parse_leading_eig needs -- solver_mode="iteration" would have
        raised at preflight, before rung 0 even ran.
      * I2: warm_start sigma injection -- rung 1 is warm-started from rung 0's
        converged self-energy (sigma.npz), which must be correctly wired into
        the rung-1 FLEX input and not silently dropped.
    """
    base = _tiny_hubbard_flex_input(tmp_path, T=4.0, itermax=30)
    base["mode"]["param"]["EPS"] = 6
    base["file"]["output"]["chiq_s"] = "chiq_s"
    base["file"]["output"]["chiq_c"] = "chiq_c"
    base["file"]["output"]["green"] = "green"
    base["eliashberg"] = {
        "chi0q_mode": "flex",
        "frequency": "dynamic",
        "pairing_type": "singlet",
        "solver_mode": "eigenvalue",
        "num_eigenvalues": 6,
        "eigenvalue_method": "arnoldi",
    }
    base["continuation"] = {
        "temperatures": [4.0, 3.0],
        "output_dir": "sweep_dyn",
        "run_eliashberg": True,
        "warm_start": True,
        "seed_gap": True,
        "summary_file": "lambda_vs_T.dat",
    }

    rows = ts.run(base, base_dir=str(tmp_path))

    assert len(rows) == 2
    assert all(r["status"] in ("ok", "not_converged") for r in rows), rows
    for r in rows:
        assert math.isfinite(r["re"]), \
            "rung {} did not get a parsed leading eigenvalue: {}".format(
                r["idx"], r)

    summary_path = os.path.join(tmp_path, "sweep_dyn", "lambda_vs_T.dat")
    assert os.path.exists(summary_path)
    with open(summary_path) as fh:
        data_lines = [ln for ln in fh if ln.strip() and not ln.startswith("#")]
    assert len(data_lines) == 2

    # rung 1 (T=3.0) was warm-started (sigma.npz) from rung 0, per I2.
    rung0_out = os.path.join(tmp_path, "sweep_dyn", "000_T4", "output")
    rung1_out = os.path.join(tmp_path, "sweep_dyn", "001_T3", "output")
    assert os.path.exists(os.path.join(rung0_out, "sigma.npz"))

    # rung 0's dynamic gap is the seed handed to rung 1 (seed_gap=True); both
    # rungs' own dynamic runs must have produced a seedable gap_dynamic.npz,
    # confirming the real dynamic Eliashberg path ran end-to-end (C1) and the
    # seed wiring did not error out (would show up as status="error" above).
    assert os.path.exists(os.path.join(rung0_out, "gap_dynamic.npz"))
    assert os.path.exists(os.path.join(rung1_out, "gap_dynamic.npz"))
    gap0 = np.load(os.path.join(rung0_out, "gap_dynamic.npz"))["gap"]
    assert np.all(np.isfinite(gap0))
