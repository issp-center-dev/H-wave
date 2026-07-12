import os
import numpy as np
import hwave.qlms as qlms


def _tiny_hubbard_flex_input(tmp_path, T=4.0, itermax=2):
    """Minimal single-orbital 2x2 k-space Hubbard FLEX input on disk."""
    d = tmp_path
    (d / "geom.dat").write_text("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n1\n0.0 0.0 0.0\n")
    (d / "transfer.dat").write_text(
        "transfer\n1\n4\n 1 1 1 1\n"
        "  1  0  0  1  1  -1.0  0.0\n -1  0  0  1  1  -1.0  0.0\n"
        "  0  1  0  1  1  -1.0  0.0\n  0 -1  0  1  1  -1.0  0.0\n")
    (d / "coulombintra.dat").write_text("coulombintra\n1\n1\n 1\n0 0 0 1 1 4.0 0.0\n")
    toml = {
        "log": {"print_level": 0},
        "mode": {"mode": "FLEX", "calc_scheme": "reduced",
                 "param": {"T": T, "CellShape": [2, 2, 1], "SubShape": [1, 1, 1],
                           "Nmat": 16, "filling": 0.5, "IterationMax": itermax,
                           "Mix": 0.5, "EPS": 6}},
        "file": {"input": {"path_to_input": "",
                           "interaction": {"path_to_input": str(d),
                                           "Geometry": "geom.dat",
                                           "Transfer": "transfer.dat",
                                           "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": str(d / "output"),
                            "sigma": "sigma"}},
    }
    return toml


def test_run_flex_returns_convergence_dict(tmp_path):
    inp = _tiny_hubbard_flex_input(tmp_path, T=4.0, itermax=2)
    result = qlms.run(input_dict=inp)
    assert isinstance(result, dict)
    assert "scf_converged" in result and isinstance(result["scf_converged"], bool)
    assert "scf_iterations" in result and isinstance(result["scf_iterations"], int)
    assert result["scf_iterations"] >= 1
