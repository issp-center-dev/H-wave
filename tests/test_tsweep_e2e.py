"""End-to-end + reentrancy tests for hwave.tsweep over real solvers.

Reuses the minimal FLEX fixture from tests/test_qlms_run_result.py (Task 1),
which is already proven to drive a real FLEX run to completion.
"""
import os

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
