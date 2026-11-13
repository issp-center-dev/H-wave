"""Single-record input files must retain their row dimension."""

import os

import numpy as np

from hwave.qlmsio.read_input_k import QLMSkInput


def test_single_row_headerless_initial_uhf(tmp_path):
    (tmp_path / "initial_uhf.dat").write_text("0 1 2 0 0.25 -0.5\n")
    reader = QLMSkInput({
        "path_to_input": os.path.relpath(tmp_path),
        "initial_uhf": "initial_uhf.dat",
    })
    assert reader.get_param("green")["initial_uhf"] == {(0, 1, 2, 0): 0.25 - 0.5j}


def test_single_row_onebodyg_uhf(tmp_path):
    (tmp_path / "greenone.def").write_text(
        "=====\nNCisAjs 1\n=====\n=====\n=====\n0 1 2 0\n")
    reader = QLMSkInput({
        "path_to_input": os.path.relpath(tmp_path),
        "onebodyg_uhf": "greenone.def",
    })
    np.testing.assert_array_equal(reader.get_param("green")["onebodyg_uhf"],
                                  [[0, 1, 2, 0]])
