"""Channel-decomposition diagnostics that do not require sparse-ir."""

import os

import numpy as np
import pytest

from tests.test_eliashberg_ir import (
    LX,
    LY,
    _eliashberg_input,
    _offsite_input,
    flex_outdir,
    flex_outdir_offsite,
)


@pytest.mark.parametrize("pairing_type", ["singlet", "triplet"])
def test_channel_decomposition_vertex_linear_offsite(flex_outdir_offsite, pairing_type):
    """Spin, charge, and retained bare vertex terms are linearly separable."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed

    inp = _offsite_input(flex_outdir_offsite)
    norb, nx, ny, nz = 1, LX, LY, 1
    chis_w, chic_w, _, convention = ed.load_flex_chi_dynamic(inp, norb, nx, ny, nz)
    kx = np.linspace(0, 2 * np.pi, nx, endpoint=False)
    ky = np.linspace(0, 2 * np.pi, ny, endpoint=False)
    kz = np.linspace(0, 2 * np.pi, nz, endpoint=False)
    _, _, interaction = sc._read_interaction_files(inp)
    inter_k = sc._build_interaction_k(kx, ky, kz, interaction, norb)
    zeros_s, zeros_c = np.zeros_like(chis_w), np.zeros_like(chic_w)

    def vertex(chis, chic):
        return ed.compute_vertices_flex_dynamic(
            chis,
            chic,
            inter_k,
            norb,
            nx,
            ny,
            nz,
            pairing_type=pairing_type,
            convention=convention,
        )

    full = vertex(chis_w, chic_w)
    spin = vertex(chis_w, zeros_c)
    charge = vertex(zeros_s, chic_w)
    bare = vertex(zeros_s, zeros_c)
    np.testing.assert_allclose(full + bare, spin + charge, rtol=1e-10, atol=1e-12)
    if pairing_type == "singlet":
        assert np.max(np.abs(bare)) > 1e-8


@pytest.mark.parametrize("pairing_type", ["singlet", "triplet"])
@pytest.mark.parametrize(
    "extra, expected",
    [
        ({"zero_chi_c": True}, (False, True)),
        ({"zero_chi_s": True}, (True, False)),
        ({"zero_chi_c": True, "zero_chi_s": True}, (True, True)),
        ({"zero_chi_c": "false", "zero_chi_s": "false"}, (False, False)),
    ],
)
def test_channel_decomposition_flags_route_through_solve_dynamic(
    flex_outdir, monkeypatch, pairing_type, extra, expected
):
    """The public uniform solver routes both pairing types without sparse-ir."""
    from hwave.solver import eliashberg_dynamic as ed

    captured = {}
    real_compute = ed.compute_vertices_flex_dynamic

    def capture(chis_w, chic_w, *args, **kwargs):
        captured["zero"] = (
            bool(np.all(chis_w == 0)),
            bool(np.all(chic_w == 0)),
        )
        return real_compute(chis_w, chic_w, *args, **kwargs)

    monkeypatch.setattr(ed, "compute_vertices_flex_dynamic", capture)
    options = dict(extra, pairing_type=pairing_type)
    ed.solve_dynamic(_eliashberg_input(flex_outdir, extra=options))
    assert captured["zero"] == expected
    if any(expected):
        with np.load(os.path.join(flex_outdir, "gap_dynamic.npz")) as data:
            assert bool(data["zero_chi_s"]) is expected[0]
            assert bool(data["zero_chi_c"]) is expected[1]
        with open(os.path.join(flex_outdir, "gap.dat")) as stream:
            header = stream.readline()
        assert "zero_chi_s={}".format(str(expected[0]).lower()) in header
        assert "zero_chi_c={}".format(str(expected[1]).lower()) in header
        with open(os.path.join(flex_outdir, "eigenvalue.dat")) as stream:
            eigenvalue_header = stream.read()
        assert "zero_chi_s={}".format(str(expected[0]).lower()) in eigenvalue_header
        assert "zero_chi_c={}".format(str(expected[1]).lower()) in eigenvalue_header
    else:
        with np.load(os.path.join(flex_outdir, "gap_dynamic.npz")) as data:
            assert "zero_chi_s" not in data and "zero_chi_c" not in data


@pytest.mark.parametrize(
    "flags, ignored",
    [
        ({"zero_chi_s": True}, ["zero_chi_s"]),
        ({"zero_chi_c": True}, ["zero_chi_c"]),
        ({"zero_chi_s": True, "zero_chi_c": True}, ["zero_chi_s", "zero_chi_c"]),
        ({"zero_chi_s": "true"}, ["zero_chi_s"]),
        ({}, []),
        ({"zero_chi_s": False, "zero_chi_c": "false"}, []),
    ],
)
def test_static_channel_flags_warn_they_are_ignored(caplog, flags, ignored):
    """The channel-decomposition flags only affect the dynamic vertex; on the
    static path they are silently inert, so calc_eliashberg must warn."""
    import logging

    import hwave.sc as sc

    with caplog.at_level(logging.WARNING, logger="hwave_sc"):
        sc._warn_if_static_ignores_channel_flags(flags)
    warnings = [
        rec.getMessage()
        for rec in caplog.records
        if rec.levelno >= logging.WARNING
    ]
    if ignored:
        assert any("frequency='static'" in msg for msg in warnings)
        for name in ignored:
            assert any(name in msg for msg in warnings)
    else:
        assert not any("frequency='static'" in msg for msg in warnings)


@pytest.mark.parametrize("freq", [None, "static"])
def test_static_channel_flags_warn_through_calc_eliashberg(caplog, monkeypatch, freq):
    """Integration: the warning is wired into calc_eliashberg's static path
    (not just the helper), for both explicit and default (omitted) static
    frequency. A ``_read_interaction_files`` stub short-circuits right after the
    warning so no fixture files are needed; the warning must already be logged
    by then."""
    import logging

    import hwave.sc as sc

    class _Stop(Exception):
        pass

    def stop(*args, **kwargs):
        raise _Stop()

    monkeypatch.setattr(sc, "_read_interaction_files", stop)

    eli = {"zero_chi_s": True}
    if freq is not None:
        eli["frequency"] = freq
    inp = {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": 8, "filling": 0.5}},
        "file": {"output": {"path_to_output": "."}},
        "eliashberg": eli,
    }
    with caplog.at_level(logging.WARNING, logger="hwave_sc"):
        with pytest.raises(_Stop):
            sc.calc_eliashberg(inp)
    warnings = [
        rec.getMessage()
        for rec in caplog.records
        if rec.levelno >= logging.WARNING
    ]
    assert any(
        "frequency='static'" in msg and "zero_chi_s" in msg for msg in warnings
    )
