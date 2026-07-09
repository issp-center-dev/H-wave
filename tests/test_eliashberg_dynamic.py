import numpy as np
import pytest


def _min_input(**eli):
    base = {
        "mode": {"param": {"Nmat": 8, "T": 0.1, "CellShape": [4, 4, 1]}},
        "file": {"output": {"path_to_output": "."}},
        "eliashberg": {"chi0q_mode": "flex"},
    }
    base["eliashberg"].update(eli)
    return base


def test_frequency_default_is_static(monkeypatch):
    import hwave.sc as sc
    assert sc._eliashberg_frequency(_min_input()) == "static"


def test_frequency_invalid_value_raises():
    import hwave.sc as sc
    with pytest.raises(ValueError, match="frequency"):
        sc._eliashberg_frequency(_min_input(frequency="dynamical"))


def test_dynamic_requires_flex():
    import hwave.sc as sc
    inp = _min_input(frequency="dynamic")
    inp["eliashberg"]["chi0q_mode"] = "load"
    with pytest.raises(ValueError, match="flex"):
        sc._validate_dynamic_prereqs(inp)


def test_dynamic_requires_even_nmat():
    import hwave.sc as sc
    inp = _min_input(frequency="dynamic")
    inp["mode"]["param"]["Nmat"] = 7
    with pytest.raises(ValueError, match="even"):
        sc._validate_dynamic_prereqs(inp)


def _write_flex_fixture(tmp_path, nmat=8, norb=1, Nx=2, Ny=2, Nz=1):
    nvol = Nx*Ny*Nz; nd = norb*norb
    rng = np.random.default_rng(3)
    def rc(shape): return rng.standard_normal(shape)+1j*rng.standard_normal(shape)
    np.savez(tmp_path/"chiq_s.npz", chiq_s=rc((nmat, nvol, nd, nd)))
    np.savez(tmp_path/"chiq_c.npz", chiq_c=rc((nmat, nvol, nd, nd)))
    np.savez(tmp_path/"green.npz",  green=rc((1, nmat, nvol, norb, norb)))
    return dict(nmat=nmat, norb=norb, Nx=Nx, Ny=Ny, Nz=Nz)


def test_load_flex_chi_dynamic_shapes(tmp_path):
    from hwave.solver import eliashberg_dynamic as ed
    m = _write_flex_fixture(tmp_path)
    inp = {"mode": {"param": {"Nmat": m["nmat"]}},
           "file": {"output": {"path_to_output": str(tmp_path)}},
           "eliashberg": {"chi0q_mode": "flex"}}
    chis, chic, green, conv = ed.load_flex_chi_dynamic(inp, m["norb"], m["Nx"], m["Ny"], m["Nz"])
    assert chis.shape == (m["Nx"], m["Ny"], m["Nz"], 1, 1, m["nmat"])
    assert green.shape == (m["norb"], m["norb"], m["Nx"], m["Ny"], m["Nz"], m["nmat"])
    assert conv in ("kuroki", "myo")


def test_load_flex_chi_dynamic_grid_mismatch(tmp_path):
    from hwave.solver import eliashberg_dynamic as ed
    m = _write_flex_fixture(tmp_path, nmat=8)
    # overwrite green with a different nmat
    np.savez(tmp_path/"green.npz", green=(np.zeros((1, 6, 4, 1, 1), complex)))
    inp = {"mode": {"param": {"Nmat": 8}},
           "file": {"output": {"path_to_output": str(tmp_path)}},
           "eliashberg": {"chi0q_mode": "flex"}}
    with pytest.raises(ValueError, match="nmat"):
        ed.load_flex_chi_dynamic(inp, 1, 2, 2, 1)


def test_check_memory_aborts_over_limit():
    from hwave.solver import eliashberg_dynamic as ed
    with pytest.raises(MemoryError, match="mem_limit_gb"):
        ed.check_memory(norb=2, Nk=1024, nmat=512, mem_limit_gb=0.01)
    ed.check_memory(norb=1, Nk=4, nmat=8, mem_limit_gb=0)  # disabled: no raise
