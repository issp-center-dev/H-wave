"""Channel-decomposition diagnostics that do not require sparse-ir."""

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
