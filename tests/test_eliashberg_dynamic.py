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


def test_gpu_requested_parses_bool_and_string():
    from hwave.solver import eliashberg_dynamic as ed
    assert ed._gpu_requested({"gpu": True}) is True
    assert ed._gpu_requested({"gpu": False}) is False
    assert ed._gpu_requested({}) is False
    # programmatic string configs: only truthy strings enable GPU
    assert ed._gpu_requested({"gpu": "true"}) is True
    assert ed._gpu_requested({"gpu": "false"}) is False   # not bool("false")==True
    assert ed._gpu_requested({"gpu": "0"}) is False


def test_gpu_true_with_static_frequency_raises():
    """gpu=true is only wired into the dynamic solver; on the static (CPU-only)
    path it must fail loudly rather than silently ignore the flag."""
    import hwave.sc as sc
    inp = _min_input(gpu=True)          # frequency omitted -> static
    with pytest.raises(ValueError, match="gpu"):
        sc.calc_eliashberg(inp)


def _write_flex_fixture(tmp_path, nmat=8, norb=1, Nx=2, Ny=2, Nz=1):
    nvol = Nx*Ny*Nz; nd = norb*norb
    rng = np.random.default_rng(3)
    def rc(shape): return rng.standard_normal(shape)+1j*rng.standard_normal(shape)
    # chi is in orbital-pair space (nd = norb^2), i.e. the general/full-vertex
    # "myo" convention; tag it so the loader treats it as orbital-pair (no
    # spin-orbital block extraction).
    np.savez(tmp_path/"chiq_s.npz", chiq_s=rc((nmat, nvol, nd, nd)),
             chi_convention="myo")
    np.savez(tmp_path/"chiq_c.npz", chiq_c=rc((nmat, nvol, nd, nd)),
             chi_convention="myo")
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


def test_npz_freq_size_reads_header_without_loading(tmp_path):
    """_npz_freq_size returns the stored frequency-axis length from the NPZ
    header (chi: axis 0; green: axis 1) without materializing the array."""
    from hwave.solver import eliashberg_dynamic as ed
    _write_flex_fixture(tmp_path, nmat=12, norb=1, Nx=2, Ny=2, Nz=1)
    assert ed._npz_freq_size(str(tmp_path / "chiq_s.npz"),
                             ("chiq_s", "chiq"), axis=0) == 12
    assert ed._npz_freq_size(str(tmp_path / "green.npz"),
                             ("green",), axis=1) == 12


def test_npz_freq_size_returns_none_on_bad_or_missing(tmp_path):
    """Header probe is best-effort: an absent file, a non-NPZ file, or a missing
    key all return None (so the caller falls back to the config grid and the
    loader raises the real error) rather than a new header-parser failure."""
    from hwave.solver import eliashberg_dynamic as ed
    assert ed._npz_freq_size(str(tmp_path / "nope.npz"),
                             ("chiq_s",), axis=0) is None
    (tmp_path / "garbage.npz").write_bytes(b"not a zip file")
    assert ed._npz_freq_size(str(tmp_path / "garbage.npz"),
                             ("chiq_s",), axis=0) is None
    _write_flex_fixture(tmp_path, nmat=8)
    assert ed._npz_freq_size(str(tmp_path / "chiq_s.npz"),
                             ("absent_key",), axis=0) is None


def test_npz_freq_size_returns_none_if_numpy_header_api_changes(tmp_path, monkeypatch):
    """The probe uses numpy.lib.format internals; if a future numpy renames them
    (AttributeError), the best-effort probe must still return None (fall back to
    the loader) rather than crash."""
    from hwave.solver import eliashberg_dynamic as ed
    from numpy.lib import format as npformat
    _write_flex_fixture(tmp_path, nmat=8)

    def _boom(*a, **k):
        raise AttributeError("simulated numpy internal API change")

    monkeypatch.setattr(npformat, "read_magic", _boom)
    assert ed._npz_freq_size(str(tmp_path / "chiq_s.npz"),
                             ("chiq_s",), axis=0) is None


def test_grid_mismatch_detected_before_loading(monkeypatch, tmp_path):
    """Issue #41: a FLEX file whose stored nmat exceeds the config Nmat must be
    rejected from the NPZ headers BEFORE the full arrays are allocated -- the
    memory guard sizes on the real file, not the smaller config."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    _write_flex_fixture(tmp_path, nmat=16, norb=1, Nx=2, Ny=2, Nz=1)

    called = {"loaded": False}

    def _must_not_run(*a, **k):
        called["loaded"] = True
        raise AssertionError("loader ran before the header-based grid check")

    monkeypatch.setattr(sc, "_load_flex_susceptibilities_full", _must_not_run)

    inp = {"mode": {"param": {"Nmat": 8}},   # config < file (16)
           "file": {"output": {"path_to_output": str(tmp_path)}},
           "eliashberg": {"chi0q_mode": "flex"}}
    with pytest.raises(ValueError, match="nmat"):
        ed.load_flex_chi_dynamic(inp, 1, 2, 2, 1)
    assert called["loaded"] is False


def _write_flex_so_fixture(tmp_path, nmat=8, norb=1, Nx=2, Ny=2, Nz=1):
    """FLEX chi in spin-orbital space (nd_chi = norb*2) to exercise the
    spin-orbital block-expansion path of the static loader."""
    nvol = Nx * Ny * Nz
    nd_so = norb * 2
    rng = np.random.default_rng(5)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    np.savez(tmp_path / "chiq_s.npz", chiq_s=rc((nmat, nvol, nd_so, nd_so)))
    np.savez(tmp_path / "chiq_c.npz", chiq_c=rc((nmat, nvol, nd_so, nd_so)))
    np.savez(tmp_path / "green.npz", green=rc((1, nmat, nvol, norb, norb)))
    return dict(nmat=nmat, norb=norb, Nx=Nx, Ny=Ny, Nz=Nz)


def test_expand_flex_chi_kuroki_norb2_extracts_spin_orbital_block():
    """A norb=2 reduced (kuroki) FLEX chi has nd_chi = norb*ns = 4 = norb^2, so
    the spin-orbital layout is INDISTINGUISHABLE from an orbital-pair layout by
    shape alone. _expand_flex_chi must use the 'kuroki' convention to extract
    the spin-up orbital block and diagonally expand it -- not treat the raw 4x4
    spin-orbital matrix as an orbital-pair susceptibility."""
    import hwave.sc as sc
    norb, ns, Nx, Ny, Nz, nfreq = 2, 2, 2, 2, 1, 3
    nvol, nd_so, nd = Nx * Ny * Nz, norb * 2, norb * norb
    rng = np.random.default_rng(11)
    chi_raw = (rng.standard_normal((nfreq, nvol, nd_so, nd_so))
               + 1j * rng.standard_normal((nfreq, nvol, nd_so, nd_so)))

    out = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention="kuroki")

    # expected: extract up-spin orbital block [:norb, :norb], scatter diagonally
    chi6 = chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_so, nd_so)
    orb = chi6[:, :, :, :, :norb, :norb]
    expected = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
    for l2 in range(norb):
        expected[:, :, :, :, l2::norb, l2::norb] = orb
    np.testing.assert_allclose(out, expected)
    # and it must NOT equal the raw 4x4 (the buggy orbital-pair interpretation)
    assert not np.allclose(out, chi6)


def test_expand_flex_chi_myo_norb2_is_orbital_pair():
    """A norb=2 general (myo) FLEX chi is already orbital-pair (nd_chi=norb^2);
    _expand_flex_chi must pass it through unchanged."""
    import hwave.sc as sc
    norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 3
    nvol, nd = Nx * Ny * Nz, norb * norb
    rng = np.random.default_rng(12)
    chi_raw = (rng.standard_normal((nfreq, nvol, nd, nd))
               + 1j * rng.standard_normal((nfreq, nvol, nd, nd)))

    out = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention="myo")
    np.testing.assert_allclose(out, chi_raw.reshape(nfreq, Nx, Ny, Nz, nd, nd))


def test_expand_flex_chi_norb2_unknown_convention_raises():
    """For the shape-ambiguous norb=2 case, an unknown/typo convention must
    raise -- not silently fall through to the spin-orbital extraction (which
    would corrupt the pairing vertex)."""
    import hwave.sc as sc
    norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 2
    nvol, nd = Nx * Ny * Nz, norb * norb
    chi_raw = np.zeros((nfreq, nvol, nd, nd), dtype=complex)
    with pytest.raises(ValueError, match="chi_convention"):
        sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention="bogus")


def test_expand_flex_chi_dimension_matches_neither_raises():
    """A chi dimension matching neither norb^2 nor norb*ns must raise clearly
    (e.g. norb=3 with nd_chi=5)."""
    import hwave.sc as sc
    norb, Nx, Ny, Nz, nfreq = 3, 2, 1, 1, 2
    nvol, nd_bad = Nx * Ny * Nz, 5
    chi_raw = np.zeros((nfreq, nvol, nd_bad, nd_bad), dtype=complex)
    with pytest.raises(ValueError, match="matches neither"):
        sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, convention="kuroki")


def test_static_loader_untagged_norb2_defaults_kuroki_and_extracts(tmp_path):
    """A legacy untagged norb=2 chi (nd_chi=4) defaults to 'kuroki' and MUST be
    spin-orbital-extracted at the loader level -- the regression the shape-only
    heuristic missed."""
    import hwave.sc as sc
    norb, ns, Nx, Ny, Nz, nmat = 2, 2, 2, 1, 1, 4
    nvol, nd_so, nd = Nx * Ny * Nz, norb * 2, norb * norb
    rng = np.random.default_rng(21)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    # no chi_convention tag -> _read_flex_chi_raw defaults to "kuroki"
    np.savez(tmp_path / "chiq_s.npz", chiq_s=rc((nmat, nvol, nd_so, nd_so)))
    np.savez(tmp_path / "chiq_c.npz", chiq_c=rc((nmat, nvol, nd_so, nd_so)))
    inp = {"mode": {"param": {"Nmat": nmat}},
           "file": {"output": {"path_to_output": str(tmp_path)}},
           "eliashberg": {"chi0q_mode": "flex"}}

    chis, chic, _, conv = sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)
    assert conv == "kuroki"
    assert chis.shape == (Nx, Ny, Nz, nd, nd)
    # verify it is the spin-orbital extraction, not the raw 4x4 passthrough
    raw_s = np.load(tmp_path / "chiq_s.npz")["chiq_s"]
    chi6 = raw_s.reshape(nmat, Nx, Ny, Nz, nd_so, nd_so)[nmat // 2]
    orb = chi6[:, :, :, :norb, :norb]
    expected = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
    for l2 in range(norb):
        expected[:, :, :, l2::norb, l2::norb] = orb
    np.testing.assert_allclose(chis, expected)


def _flex_input(tmp_path, nmat):
    return {"mode": {"param": {"Nmat": nmat}},
            "file": {"output": {"path_to_output": str(tmp_path)}},
            "eliashberg": {"chi0q_mode": "flex"}}


def test_static_flex_loader_matches_full_then_slice(tmp_path):
    """The static FLEX loader must return exactly the static-frequency slice of
    the full-frequency loader (correctness preserved by the slice-first path)."""
    import hwave.sc as sc
    m = _write_flex_so_fixture(tmp_path)
    inp = _flex_input(tmp_path, m["nmat"])
    chis, chic, _, conv = sc._load_flex_susceptibilities(
        inp, m["norb"], m["Nx"], m["Ny"], m["Nz"])
    chis_w, chic_w, _, conv_w = sc._load_flex_susceptibilities_full(
        inp, m["norb"], m["Nx"], m["Ny"], m["Nz"])
    center = m["nmat"] // 2
    np.testing.assert_allclose(chis, chis_w[..., center])
    np.testing.assert_allclose(chic, chic_w[..., center])
    assert conv == conv_w


def test_static_flex_loader_expands_single_frequency(tmp_path, monkeypatch):
    """The static loader must slice the static frequency BEFORE the spin-orbital
    expansion -- expanding one frequency, not the full Nmat axis (the memory
    regression). The full loader still expands the whole axis."""
    import hwave.sc as sc
    m = _write_flex_so_fixture(tmp_path, nmat=8)
    inp = _flex_input(tmp_path, m["nmat"])

    seen = []
    real_expand = sc._expand_flex_chi

    def _spy(chi_raw, *a, **k):
        seen.append(chi_raw.shape[0])    # leading axis = frequency count
        return real_expand(chi_raw, *a, **k)

    monkeypatch.setattr(sc, "_expand_flex_chi", _spy)
    sc._load_flex_susceptibilities(inp, m["norb"], m["Nx"], m["Ny"], m["Nz"])

    assert seen and max(seen) == 1, \
        "static loader expanded {} frequencies (want 1)".format(max(seen))


def test_check_memory_aborts_over_limit():
    from hwave.solver import eliashberg_dynamic as ed
    with pytest.raises(MemoryError, match="mem_limit_gb"):
        ed.check_memory(norb=2, Nk=1024, nmat=512, mem_limit_gb=0.01)
    ed.check_memory(norb=1, Nk=4, nmat=8, mem_limit_gb=0)  # disabled: no raise


def test_memory_guard_runs_before_loading(monkeypatch):
    """The memory guard must abort BEFORE _load_flex_susceptibilities_full
    allocates the full-frequency chi/green arrays (fail-before-allocating):
    a too-tight mem_limit must raise MemoryError without the loader running."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed

    called = {"loaded": False}

    def _must_not_run(*args, **kwargs):
        called["loaded"] = True
        raise AssertionError("loader ran before the memory guard")

    monkeypatch.setattr(sc, "_load_flex_susceptibilities_full", _must_not_run)

    inp = {"mode": {"param": {"Nmat": 512}},
           "file": {"output": {"path_to_output": "."}},
           "eliashberg": {"chi0q_mode": "flex", "mem_limit_gb": 0.01}}
    with pytest.raises(MemoryError, match="mem_limit_gb"):
        ed.load_flex_chi_dynamic(inp, 2, 8, 8, 1)
    assert called["loaded"] is False


def test_dynamic_vertex_matches_static_per_frequency():
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    norb, Nx, Ny, Nz, nmat = 1, 2, 2, 1, 4
    nd = norb*norb
    rng = np.random.default_rng(5)
    chis_w = rng.standard_normal((Nx, Ny, Nz, nd, nd, nmat)) + 0j
    chic_w = rng.standard_normal((Nx, Ny, Nz, nd, nd, nmat)) + 0j
    U = rng.standard_normal((norb, norb, Nx, Ny, Nz)) + 0j
    inter_k = {"CoulombIntra": U}
    Vw = ed.compute_vertices_flex_dynamic(chis_w, chic_w, inter_k, norb,
                                          Nx, Ny, Nz, pairing_type="singlet",
                                          convention="kuroki")
    # compare frequency l against the static routine fed that single slice
    for l in range(nmat):
        Vstat = sc._compute_vertices_flex(chis_w[..., l], chic_w[..., l], inter_k,
                                          norb, Nx, Ny, Nz, pairing_type="singlet",
                                          convention="kuroki")
        assert np.allclose(Vw[..., l], Vstat, atol=1e-12)


def test_kernel_single_k_matches_dense_freq_operator():
    from hwave.solver import eliashberg_dynamic as ed
    from hwave.solver import matsubara as ms
    norb, Nx, Ny, Nz, nmat = 1, 1, 1, 1, 4
    beta = 8.0
    rng = np.random.default_rng(7)
    V = (rng.standard_normal((1, 1, 1, 1, Nx, Ny, Nz, nmat))
         + 1j*rng.standard_normal((1, 1, 1, 1, Nx, Ny, Nz, nmat)))
    G2 = (rng.standard_normal((1, 1, 1, 1, Nx, Ny, Nz, nmat))
          + 1j*rng.standard_normal((1, 1, 1, 1, Nx, Ny, Nz, nmat)))
    phi = (rng.standard_normal((1, 1, Nx, Ny, Nz, nmat))
           + 1j*rng.standard_normal((1, 1, Nx, Ny, Nz, nmat)))
    out = ed.eliashberg_kernel_dynamic(V, G2, phi, norb, beta)
    # dense freq operator at single k (q=0): -(1/1) Bfi diag(Vtau) Bf diag(G2) phi
    I = np.eye(nmat, dtype=complex)
    Bf = np.stack([ms.fermion_to_tau(I[i], axis=0) for i in range(nmat)], 1)
    Bfi = np.stack([ms.tau_to_fermion(I[i], axis=0) for i in range(nmat)], 1)
    Bb = np.stack([ms.boson_to_tau(I[i], axis=0) for i in range(nmat)], 1)
    Vtau = Bb @ V.reshape(nmat)
    op = -(Bfi @ (Vtau[:, None] * Bf) * G2.reshape(nmat)[None, :])
    expected = (op @ phi.reshape(nmat)).reshape(1, 1, Nx, Ny, Nz, nmat)
    assert np.allclose(out, expected, atol=1e-10)


def test_kernel_precomputed_vertex_rt_matches():
    """The pairing vertex's (q, i nu) -> (r, tau) transform is phi-independent,
    so it can be hoisted out of the power-iteration/Arnoldi matvec: applying
    the kernel with the precomputed real-space/tau vertex (``Vs_rt``) must
    match the plain call that transforms ``Vs_q_w`` internally. (Comparison is
    at fp64 round-off, not bit-identity: numpy's FFT/einsum results vary in
    the last bits with buffer alignment, so even two identical kernel calls
    are not bit-reproducible.)"""
    from hwave.solver import eliashberg_dynamic as ed
    norb, Nx, Ny, Nz, nmat = 2, 2, 2, 1, 8
    beta = 7.0
    rng = np.random.default_rng(23)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    V_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    G2_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    phi = rc((norb, norb, Nx, Ny, Nz, nmat))

    ref = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi, norb, beta)
    V_rt = ed.vertex_qw_to_rt(V_w)
    out = ed.eliashberg_kernel_dynamic(None, G2_w, phi, norb, beta,
                                       Vs_rt=V_rt)
    np.testing.assert_allclose(out, ref, rtol=1e-12, atol=1e-12)


def test_kernel_fft_workers_matches_serial():
    """The scipy-parallel spatial FFT path (workers=-1) matches the serial
    numpy path to machine precision. (scipy's FFT backend may differ from
    numpy's in the last bits, hence allclose rather than bit equality.)"""
    from hwave.solver import eliashberg_dynamic as ed
    if ed._SFFT is None:
        pytest.skip("scipy.fft unavailable")
    norb, Nx, Ny, Nz, nmat = 2, 1, 8, 8, 8
    rng = np.random.default_rng(31)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    V = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    G2 = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    phi = rc((norb, norb, Nx, Ny, Nz, nmat))
    serial = ed.eliashberg_kernel_dynamic(
        None, G2, phi, norb, 5.0,
        Vs_rt=ed.vertex_qw_to_rt(V, workers=1), workers=1)
    par = ed.eliashberg_kernel_dynamic(
        None, G2, phi, norb, 5.0,
        Vs_rt=ed.vertex_qw_to_rt(V, workers=-1), workers=-1)
    assert np.max(np.abs(serial - par)) < 1e-11


def test_frequency_inner_is_full_vdot():
    from hwave.solver import eliashberg_dynamic as ed
    rng = np.random.default_rng(11)
    a = rng.standard_normal((2, 2, 2, 2, 1, 3)) + 1j*rng.standard_normal((2, 2, 2, 2, 1, 3))
    b = rng.standard_normal(a.shape) + 1j*rng.standard_normal(a.shape)
    assert np.isclose(ed.frequency_inner(a, b), np.vdot(a, b))


@pytest.mark.parametrize("norb", [1, 2])
def test_g2_dynamic_sums_to_static(norb):
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    Nx, Ny, Nz, nmat = 2, 2, 1, 8
    rng = np.random.default_rng(6)
    # uniquely-populated (non-symmetric in orbital) green so an l-axis
    # transposition cannot pass the norb=2 case
    green = (rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
             + 1j*rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat)))
    beta = 10.0
    g2w = ed.calc_g2_dynamic(green, beta)
    assert g2w.shape == (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    g2_static = sc._calc_g2(green, beta)
    assert np.allclose(g2w.sum(axis=-1), g2_static, atol=1e-12)
    # component check (independent of the sum): reproduce a single element
    green_inv = np.roll(green[:, :, ::-1, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    l2, l5, l3, l6, k, n = 0, min(1, norb-1), min(1, norb-1), 0, 3, 5
    kx, ky, kz = 1, 0, 0
    assert np.isclose(
        g2w[l2, l5, l3, l6, kx, ky, kz, n],
        green[l2, l5, kx, ky, kz, n] * green_inv[l3, l6, kx, ky, kz, n] / beta,
        atol=1e-12)


def _write_geom_transfer_coulomb(input_dir, norb=1):
    import os
    os.makedirs(input_dir, exist_ok=True)
    with open(os.path.join(input_dir, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
        f.write("{}\n".format(norb))
        for _ in range(norb):
            f.write("0.0 0.0 0.0\n")
    with open(os.path.join(input_dir, "transfer.dat"), "w") as f:
        f.write("Transfer\n{}\n5\n 1 1 1 1 1\n".format(norb))
        for orb in range(norb):
            f.write("  1  0  0  {0}  {0}  1.0 0.0\n".format(orb + 1))
            f.write(" -1  0  0  {0}  {0}  1.0 0.0\n".format(orb + 1))
            f.write("  0  1  0  {0}  {0}  1.0 0.0\n".format(orb + 1))
            f.write("  0 -1  0  {0}  {0}  1.0 0.0\n".format(orb + 1))
    with open(os.path.join(input_dir, "coulombintra.dat"), "w") as f:
        f.write("CoulombIntra\n{}\n1\n 1\n".format(norb))
        for orb in range(norb):
            f.write("  0  0  0  {0}  {0}  2.0 0.0\n".format(orb + 1))


def test_solve_dynamic_smoke_returns_finite_lambda(tmp_path):
    """End-to-end smoke: solve_dynamic runs on a tiny generated FLEX fixture
    and returns a finite leading eigenvalue. (Full static-limit / dense-oracle
    acceptance is Task 7.)"""
    import os
    from hwave.solver import eliashberg_dynamic as ed
    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1, Nx=2, Ny=2, Nz=1)

    inp = {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": m["nmat"],
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic",
                       "solver_mode": "iteration", "max_iter": 50},
    }
    lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    assert os.path.exists(os.path.join(output_dir, "eigenvalue.dat"))
    assert os.path.exists(os.path.join(output_dir, "gap_dynamic.npz"))


def test_kernel_cupy_matches_numpy():
    """The dynamic kernel applied to cupy arrays (GPU) must reproduce the
    numpy result to fp64 round-off, and stay on the device."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    from hwave.solver import eliashberg_dynamic as ed
    norb, Nx, Ny, Nz, nmat = 2, 4, 2, 1, 8
    beta = 7.0
    rng = np.random.default_rng(29)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    V_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    G2_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    phi = rc((norb, norb, Nx, Ny, Nz, nmat))

    ref = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi, norb, beta)

    V_rt_g = ed.vertex_qw_to_rt(cupy.asarray(V_w))
    assert isinstance(V_rt_g, cupy.ndarray)
    out_g = ed.eliashberg_kernel_dynamic(None, cupy.asarray(G2_w),
                                         cupy.asarray(phi), norb, beta,
                                         Vs_rt=V_rt_g)
    assert isinstance(out_g, cupy.ndarray)
    np.testing.assert_allclose(cupy.asnumpy(out_g), ref, atol=1e-10)
    # a host-side (numpy) trial gap must also be accepted (the iterative
    # solvers hand the matvec numpy vectors)
    out_h = ed.eliashberg_kernel_dynamic(None, cupy.asarray(G2_w),
                                         phi, norb, beta, Vs_rt=V_rt_g)
    np.testing.assert_allclose(cupy.asnumpy(out_h), ref, atol=1e-10)


def _dynamic_input_dict(input_dir, output_dir, nmat, extra_eliashberg=None):
    eli = {"chi0q_mode": "flex", "frequency": "dynamic",
           "solver_mode": "iteration", "max_iter": 50}
    if extra_eliashberg:
        eli.update(extra_eliashberg)
    return {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": nmat,
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": eli,
    }


def test_solve_dynamic_gpu_falls_back_without_cupy(tmp_path, monkeypatch,
                                                   caplog):
    """[eliashberg] gpu=true on a machine without CuPy must warn and produce
    the identical result via the numpy path."""
    import os
    import logging
    from hwave.solver import backend
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1,
                            Nx=2, Ny=2, Nz=1)

    lam_cpu = ed.solve_dynamic(
        _dynamic_input_dict(input_dir, output_dir, m["nmat"]))

    def _no_cupy():
        raise ImportError("No module named 'cupy'")

    monkeypatch.setattr(backend, "_import_cupy", _no_cupy)
    with caplog.at_level(logging.WARNING, logger="qlms"):
        lam_gpu_flag = ed.solve_dynamic(
            _dynamic_input_dict(input_dir, output_dir, m["nmat"],
                                extra_eliashberg={"gpu": True}))
    assert any("cupy" in rec.message.lower() for rec in caplog.records)
    assert np.isclose(lam_gpu_flag, lam_cpu, atol=1e-12)


def test_solve_dynamic_gpu_end_to_end_matches_cpu(tmp_path):
    """gpu=true with a real CUDA device must reproduce the CPU lambda."""
    cupy = pytest.importorskip("cupy")
    try:
        cupy.zeros(1)
    except Exception:
        pytest.skip("cupy installed but no usable CUDA device")
    import os
    from hwave.solver import eliashberg_dynamic as ed

    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1,
                            Nx=2, Ny=2, Nz=1)

    lam_cpu = ed.solve_dynamic(
        _dynamic_input_dict(input_dir, output_dir, m["nmat"]))
    lam_gpu = ed.solve_dynamic(
        _dynamic_input_dict(input_dir, output_dir, m["nmat"],
                            extra_eliashberg={"gpu": True}))
    assert np.isclose(lam_gpu, lam_cpu, atol=1e-10)


def test_solve_dynamic_requires_green(tmp_path):
    """green_w=None (missing green.npz) must fail-fast with a clear message,
    not an AttributeError (Task 3 review gap / spec section 9)."""
    import os
    from hwave.solver import eliashberg_dynamic as ed
    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1, Nx=2, Ny=2, Nz=1)
    os.remove(os.path.join(output_dir, "green.npz"))  # drop the dressed Green

    inp = {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": m["nmat"],
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic"},
    }
    with pytest.raises(ValueError, match="green.npz"):
        ed.solve_dynamic(inp)


# ---------------------------------------------------------------------------
# Task 7: acceptance tests (convention pins), gauge, and outputs
# ---------------------------------------------------------------------------

def _dense_operator(matvec, vec_size, dtype=complex):
    """Reconstruct the dense matrix of a linear ``matvec`` (flat->flat) by
    applying it to each unit basis vector (one column per call)."""
    K = np.empty((vec_size, vec_size), dtype=dtype)
    e = np.zeros(vec_size, dtype=dtype)
    for j in range(vec_size):
        e[:] = 0.0
        e[j] = 1.0
        K[:, j] = matvec(e)
    return K


@pytest.mark.parametrize("norb", [1, 2])
def test_static_limit_matvec_multi_k(norb):
    """PRIMARY convention pin (independent of the dynamic kernel's own
    self-consistency): on a NON-symmetric mesh (Nx=4, Ny=2) the frequency-flat
    limit of ``eliashberg_kernel_dynamic`` must reproduce, component-by-
    component, the STATIC kernel ``sc._eliashberg_kernel_fft`` evaluated with
    ``G2_static = sum_n G2_w``. This pins the spatial fold and the -(1/N)
    normalization against the static reference, not against a self-extraction.
    Nk=1 would hide both the fold and the 1/N, so the mesh must be multi-k and
    non-symmetric (q = k - k' is then not self-symmetric).

    The ``norb=2`` case additionally pins the multi-orbital orbital einsums
    (``'iljmxyzn,lmxyzn->ijxyzn'`` and ``'abcdxyzt,bcxyzt->adxyzt'``) against the
    INDEPENDENT static 4-index kernel. Because ``V0``/``G2`` are drawn from
    distinct random complex values per orbital-index combination (never
    symmetrized), they are NOT invariant under the orbital-index swaps those
    einsums perform, so an accidental orbital transposition would change the
    result and fail this assertion."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    Nx, Ny, Nz, nmat = 4, 2, 1, 4
    beta = 7.0  # unused by the kernel; G2 is passed raw to both paths
    rng = np.random.default_rng(20 + norb)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    # frequency-flat vertex: build one bosonic slice, broadcast across nmat.
    # ORBITAL-NON-SYMMETRIC: each (l1,l2,l3,l4) combo gets an independent random
    # value, so an accidental orbital transposition WOULD change the result.
    V0 = rc((norb, norb, norb, norb, Nx, Ny, Nz))            # (l1,l2,l3,l4,q)
    V_w = np.broadcast_to(V0[..., np.newaxis],
                          V0.shape + (nmat,)).copy()
    # frequency-resolved pair bubble (varies with n on purpose), also
    # orbital-non-symmetric
    G2_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    # frequency-flat trial gap
    phi0 = rc((norb, norb, Nx, Ny, Nz))                      # (l1,l2,k)
    phi_w = np.broadcast_to(phi0[..., np.newaxis],
                            phi0.shape + (nmat,)).copy()

    out = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi_w, norb, beta)

    # (i) the output must be frequency-flat
    assert np.max(np.abs(out - out[..., :1])) < 1e-10, \
        "frequency-flat input did not give a frequency-flat output"

    # (ii) each frequency slice equals the static matvec with G2_static
    G2_static = G2_w.sum(axis=-1)                            # sum over n
    static_out = sc._eliashberg_kernel_fft(V0, G2_static, phi0, norb)
    for n in range(nmat):
        assert np.max(np.abs(out[..., n] - static_out)) < 1e-10, \
            "dynamic kernel does not match the static kernel in the flat limit"


@pytest.mark.parametrize("norb,Nx,Ny,Nz,nmat", [
    (1, 2, 2, 1, 4),   # dense oracle, norb=1
    (2, 2, 1, 1, 2),   # dense oracle, norb=2 (catches orbital transposition)
])
def test_dense_oracle(norb, Nx, Ny, Nz, nmat):
    """Build the dense operator column-by-column by applying
    ``eliashberg_kernel_dynamic`` to unit basis vectors, take its full spectrum
    with ``scipy.linalg.eig``, and assert the leading eigenvalue (largest real
    part, per ``sc._order_eigenpairs``) equals the shared driver's leading
    eigenvalue. Uses a NONZERO-nu synthetic vertex so the tau-product frequency
    structure is exercised. norb=2 catches an orbital-index transposition."""
    import scipy.linalg
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    from scipy.sparse.linalg import LinearOperator
    rng = np.random.default_rng(31 + norb)
    beta = 5.0

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    V_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))   # nonzero-nu vertex
    G2_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    Nk = Nx * Ny * Nz
    gap_shape = (norb, norb, Nx, Ny, Nz, nmat)
    vec_size = norb * norb * Nk * nmat

    def _matvec(x):
        return ed.eliashberg_kernel_dynamic(
            V_w, G2_w, x.reshape(gap_shape), norb, beta).ravel()

    # dense operator + full spectrum
    K = _dense_operator(_matvec, vec_size)
    dense_vals = scipy.linalg.eig(K, right=False)
    dense_leading = dense_vals[np.argmax(dense_vals.real)]

    # shared driver (same path solve_dynamic uses for the eigenvalue family)
    op = LinearOperator((vec_size, vec_size), matvec=_matvec, dtype=complex)
    lam, _, info = sc._solve_leading(
        lambda: (op, vec_size), vec_size, "arnoldi",
        num_eigenvalues=vec_size - 2)

    assert np.isclose(lam, dense_leading, atol=1e-10), \
        "driver leading eigenvalue {} != dense leading {}".format(
            lam, dense_leading)


@pytest.mark.parametrize("solver_mode", ["arnoldi", "shift-invert-gmres"])
def test_solve_leading_tiny_operator_uses_dense_fallback(solver_mode):
    """ARPACK cannot handle LinearOperator eigenproblems with k >= N - 1.
    The smallest dynamic grid has vec_size=2, so the shared eigenvalue driver
    must use a dense fallback instead of calling eigs."""
    import hwave.sc as sc
    from scipy.sparse.linalg import LinearOperator

    dense = np.array([[2.0, 0.0],
                      [0.0, -3.0]], dtype=complex)
    op = LinearOperator((2, 2), matvec=lambda x: dense @ x, dtype=complex)

    lam, vec, info = sc._solve_leading(
        lambda: (op, 2), 2, solver_mode, num_eigenvalues=10)

    assert np.isclose(lam, 2.0)
    assert info["eigenvalues"].shape == (2,)
    assert np.allclose(dense @ vec, lam * vec)


def test_analytic_flat_diagonal():
    """Closed-form pin: a frequency-flat, orbital-diagonal, k-uniform vertex V0
    with a k-uniform (but frequency-resolved) diagonal pair bubble G2(n) makes
    the k-space kernel rank-one; its single nonzero eigenvalue is the uniform
    gap with lambda = -V0 * sum_n G2(n). (The 1/N spatial fold is cancelled by
    the sum over the Nk-fold uniform mesh.) Assert the driver reproduces it."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    from scipy.sparse.linalg import LinearOperator
    norb, Nx, Ny, Nz, nmat = 1, 2, 2, 1, 4
    beta = 3.0
    Nk = Nx * Ny * Nz
    rng = np.random.default_rng(42)

    v0 = -1.3 + 0.0j                                    # real vertex value
    g = rng.standard_normal(nmat) + 1j * rng.standard_normal(nmat)  # G2(n)
    # V frequency-flat, k-uniform, orbital-diagonal (norb=1 => scalar)
    V_w = np.full((norb, norb, norb, norb, Nx, Ny, Nz, nmat), v0, dtype=complex)
    # G2 k-uniform, frequency-resolved
    G2_w = np.empty((norb, norb, norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
    for n in range(nmat):
        G2_w[..., n] = g[n]

    gap_shape = (norb, norb, Nx, Ny, Nz, nmat)
    vec_size = norb * norb * Nk * nmat

    def _matvec(x):
        return ed.eliashberg_kernel_dynamic(
            V_w, G2_w, x.reshape(gap_shape), norb, beta).ravel()

    op = LinearOperator((vec_size, vec_size), matvec=_matvec, dtype=complex)
    lam, _, _ = sc._solve_leading(
        lambda: (op, vec_size), vec_size, "arnoldi",
        num_eigenvalues=vec_size - 2)

    lam_expected = -v0 * g.sum()
    assert np.isclose(lam, lam_expected, atol=1e-9), \
        "analytic lambda mismatch: {} vs {}".format(lam, lam_expected)


def test_gauge_deterministic():
    """_fix_gauge is invariant to a random global phase: L2-normalize over all
    components, then rotate so the largest-|magnitude| component is real-
    positive (lexicographic tie-break). Applying an arbitrary phase before
    _fix_gauge must reproduce the same array."""
    from hwave.solver import eliashberg_dynamic as ed
    rng = np.random.default_rng(77)
    phi = (rng.standard_normal((2, 2, 3, 2, 1, 4))
           + 1j * rng.standard_normal((2, 2, 3, 2, 1, 4)))

    a = ed._fix_gauge(phi)
    theta = 1.2345
    b = ed._fix_gauge(np.exp(1j * theta) * phi)
    assert np.allclose(a, b, atol=1e-12)

    # gauge property: unit norm, and the largest-|magnitude| entry is real>0
    assert np.isclose(np.linalg.norm(a), 1.0, atol=1e-12)
    flat = a.ravel()
    piv = flat[np.argmax(np.abs(flat))]
    assert piv.real > 0 and abs(piv.imag) < 1e-12

    # explicit equal-magnitude tie-break: two components share the exact max
    # magnitude (|1.0| == |1j| == 1); the argmax-first rule must pick the FIRST
    # in C-order ravel and make IT real-positive (deterministically).
    tie = np.zeros((2, 2, 1, 1, 1, 1), dtype=complex)
    tie[0, 0, 0, 0, 0, 0] = 1.0    # C-order index 0 (first tied component)
    tie[0, 1, 0, 0, 0, 0] = 1j     # C-order index 1 (second tied component)
    g = ed._fix_gauge(tie)
    # deterministic under a global phase rotation
    g2 = ed._fix_gauge(np.exp(1j * 0.7) * tie)
    assert np.allclose(g, g2, atol=1e-12)
    gflat = g.ravel()
    # the FIRST tied component was made real-positive
    assert gflat[0].real > 0 and abs(gflat[0].imag) < 1e-12
    # and it is (one of) the max-|magnitude| entries
    assert np.isclose(abs(gflat[0]), np.max(np.abs(gflat)), atol=1e-12)


def test_outputs_written(tmp_path):
    """solve_dynamic writes gap_dynamic.npz (with the required metadata keys)
    and a gap.dat whose first line starts with '#' and marks frequency=dynamic."""
    import os
    from hwave.solver import eliashberg_dynamic as ed
    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1, Nx=2, Ny=2, Nz=1)

    inp = {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": m["nmat"],
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic",
                       "pairing_type": "singlet",
                       "solver_mode": "iteration", "max_iter": 50},
    }
    lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)

    npz = np.load(os.path.join(output_dir, "gap_dynamic.npz"))
    assert set(npz.files) == {
        "gap", "iomega", "T", "pairing_type", "frequency",
        "eigenvalue", "axis_order", "normalization"}, \
        "unexpected gap_dynamic.npz key set: {}".format(sorted(npz.files))
    assert str(npz["frequency"]) == "dynamic"
    assert npz["iomega"].shape == (m["nmat"],)
    assert npz["gap"].shape == (1, 1, 2, 2, 1, m["nmat"])

    with open(os.path.join(output_dir, "gap.dat")) as fh:
        first = fh.readline()
    assert first.startswith("#")
    assert "frequency=dynamic" in first


# ---------------------------------------------------------------------------
# Task 8: slow end-to-end through the REAL FLEX pipeline + gap symmetry
# ---------------------------------------------------------------------------

def _run_real_flex(input_dir, output_dir, T, Nmat, cell_shape):
    """Drive a REAL FLEX solve through hwave.qlms.run so that genuine
    full-frequency chiq_s.npz / chiq_c.npz and a dressed green.npz are written
    to ``output_dir``. This is deliberately NOT the synthetic
    ``_write_flex_fixture`` (hand-written random npz) used by the fast smoke
    tests -- the e2e must exercise the actual FLEX self-consistency so the
    dynamic Eliashberg reads physically consistent susceptibilities and Green.
    ``green`` is named in [file.output] because FLEX only writes green.npz when
    an output key for it is present.
    """
    import hwave.qlms as qlms
    inp = {
        "mode": {"mode": "FLEX", "calc_scheme": "reduced",
                 "param": {"T": T, "mu": 0.0, "CellShape": cell_shape,
                           "SubShape": [1, 1, 1], "Nmat": Nmat,
                           "IterationMax": 300, "Mix": 0.2, "EPS": 6}},
        "file": {"input": {"path_to_input": "",
                           "interaction": {"path_to_input": input_dir,
                                           "Geometry": "geom.dat",
                                           "Transfer": "transfer.dat",
                                           "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir,
                            "chiq_s": "chiq_s", "chiq_c": "chiq_c",
                            "green": "green"}},
    }
    qlms.run(input_dict=inp)


def _eliashberg_input(input_dir, output_dir, T, Nmat, cell_shape, **eli):
    inp = {
        "mode": {"param": {"T": T, "CellShape": cell_shape,
                           "SubShape": [1, 1, 1], "Nmat": Nmat,
                           "filling": 0.5}},
        "file": {"input": {"interaction": {"path_to_input": input_dir,
                                           "Geometry": "geom.dat",
                                           "Transfer": "transfer.dat",
                                           "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"chi0q_mode": "flex", "pairing_type": "singlet",
                       "solver_mode": "iteration", "max_iter": 300},
    }
    inp["eliashberg"].update(eli)
    return inp


def _gauge_fix_k(vec):
    """Gauge-fix a flat complex gap-over-k vector the same way
    ``ed._fix_gauge`` does: L2-normalize, then rotate the largest-|amplitude|
    component to be real-positive (argmax gives the deterministic C-order
    tie-break). Both eigenvectors are only defined up to a complex scale, so
    this is required before any component-wise comparison."""
    v = np.asarray(vec).astype(complex, copy=True)
    nrm = np.linalg.norm(v)
    if nrm > 0:
        v /= nrm
    piv = v.ravel()[int(np.argmax(np.abs(v)))]
    if piv != 0:
        v /= (piv / abs(piv))
    return v


def _parse_static_gap(path, norb=1):
    """Parse sc._save_results' gap.dat (norb=1: kx ky kz Re Im per row) into a
    flat k-ordered complex vector, in the SAME (ix,iy,iz) C-order that both the
    static writer and the dynamic gap array use."""
    assert norb == 1
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.split()
            rows.append(complex(float(p[3]), float(p[4])))
    return np.array(rows)


@pytest.mark.slow
def test_end_to_end_dynamic_matches_static_symmetry(tmp_path):
    """Slow END-TO-END guarantee through the REAL FLEX pipeline.

    Unlike the fast smoke tests above (which feed a hand-written synthetic
    ``_write_flex_fixture``), this drives the ACTUAL FLEX solver:

      1. generate a tiny 1-orbital square-lattice Hubbard model (U=2.0),
      2. run a real FLEX solve via ``hwave.qlms.run`` -> genuine full-frequency
         ``chiq_s.npz`` / ``chiq_c.npz`` + dressed ``green.npz`` (nmat=8),
      3. run ``sc.calc_eliashberg`` STATIC (frequency omitted/"static") on that
         output -> lambda_static + static gap,
      4. run ``sc.calc_eliashberg`` DYNAMIC (frequency="dynamic") on the SAME
         output -> lambda_dynamic + gap_dynamic.npz,
      5. assert finite lambda, all output files exist, and the dynamic gap's
         low-|w_n| slice reproduces the static gap's k-space sign/nodal pattern.

    Grid/temperature: CellShape=[4,4,1] (Nk=16), Nmat=8, T=0.5, U=2.0. FLEX
    converges in well under a second at these sizes; this is the smallest grid
    that keeps a genuine k-structured gap (an anisotropic s-wave here, values
    ~0.97..1.0 over k) rather than a perfectly flat one.

    Gap-symmetry invariant (robustness vs. rigor):

      * Both the static gap phi_stat(k) and each central-frequency slice of the
        dynamic gap phi_dyn(k, iw_n) are eigenvectors defined only up to a
        complex scale, so each is independently gauge-fixed (``_gauge_fix_k``,
        same convention as ``ed._fix_gauge``) before comparison.
      * PRIMARY (the brief's invariant): the SIGN of every k-component whose
        |amplitude| exceeds 20% of the per-slice max must agree between static
        and dynamic. Thresholding on |amp| is what makes it robust -- the sign
        of a near-zero component is numerical noise, so a strict all-k
        elementwise-sign match would be brittle; the threshold keeps only the
        physically meaningful components. This still catches any k where the
        dynamic low-frequency limit would flip the pairing sign vs. the static
        solver (the classic s- vs d-wave nodal distinction).
      * SECONDARY (extra teeth, so the sign test is not vacuous on a nodeless
        gap): the gauge-fixed dynamic central slice must also lie within 5%
        relative-L2 of the static gap. The genuine frequency renormalization of
        the lowest-|w_n| slice is <1% at these parameters, so 5% is a wide
        margin against linear-algebra backend jitter yet fails on any gross
        k-structure corruption or a single-k sign flip (which would push the
        relative-L2 to O(1)).
      * The two central slices (n = Nmat//2-1 and Nmat//2) must be complex
        conjugates: phi(k,-w_n) = conj phi(k,w_n) is an exact symmetry of the
        linearized gap and holds here to ~1e-16.
    """
    import os
    import hwave.sc as sc

    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)

    T = 0.5
    Nmat = 8
    cell_shape = [4, 4, 1]

    # --- Step 2: REAL FLEX solve -> full-frequency chiq_s/chiq_c + green ---
    _run_real_flex(input_dir, output_dir, T, Nmat, cell_shape)
    for fn in ("chiq_s.npz", "chiq_c.npz", "green.npz"):
        assert os.path.exists(os.path.join(output_dir, fn)), \
            "FLEX did not write {}".format(fn)
    # confirm the FLEX outputs carry the FULL frequency axis (nmat=8)
    assert np.load(os.path.join(output_dir, "chiq_s.npz"))["chiq_s"].shape[0] == Nmat
    assert np.load(os.path.join(output_dir, "green.npz"))["green"].shape[1] == Nmat

    # --- Step 3: STATIC Eliashberg on the real FLEX output ---
    # separate output filenames so the static run's own outputs survive the
    # dynamic run (which writes gap.dat / eigenvalue.dat).
    sc.calc_eliashberg(_eliashberg_input(
        input_dir, output_dir, T, Nmat, cell_shape,
        output_gap="gap_static.dat", output_eigenvalue="eigenvalue_static.dat"))
    assert os.path.exists(os.path.join(output_dir, "gap_static.dat"))
    assert os.path.exists(os.path.join(output_dir, "eigenvalue_static.dat"))
    gap_static = _gauge_fix_k(_parse_static_gap(
        os.path.join(output_dir, "gap_static.dat")))

    # --- Step 4: DYNAMIC Eliashberg on the SAME real FLEX output ---
    lam_dynamic = sc.calc_eliashberg(_eliashberg_input(
        input_dir, output_dir, T, Nmat, cell_shape, frequency="dynamic"))
    assert np.isfinite(lam_dynamic)
    assert os.path.exists(os.path.join(output_dir, "gap.dat"))
    assert os.path.exists(os.path.join(output_dir, "gap_dynamic.npz"))

    npz = np.load(os.path.join(output_dir, "gap_dynamic.npz"))
    gap_dyn = npz["gap"]                        # (1, 1, Nx, Ny, Nz, nmat)
    assert gap_dyn.shape == (1, 1, 4, 4, 1, Nmat)
    assert str(npz["frequency"]) == "dynamic"

    # --- Step 5: gap-symmetry invariant on the two central Matsubara slices ---
    central = (Nmat // 2 - 1, Nmat // 2)
    slices = {}
    for n0 in central:
        gd = _gauge_fix_k(gap_dyn[0, 0, :, :, 0, n0].ravel())
        slices[n0] = gd

        # gauge-fixed slice is (numerically) real
        assert np.max(np.abs(gd.imag)) < 1e-8, \
            "dynamic central slice n0={} not real after gauge fix".format(n0)

        # PRIMARY: sign pattern agrees on the physically meaningful components
        thr = 0.2 * np.max(np.abs(gd))
        mask = np.abs(gd) > thr
        assert np.array_equal(np.sign(gd.real[mask]),
                              np.sign(gap_static.real[mask])), \
            ("dynamic (n0={}) and static gap sign patterns disagree on "
             "|amp|>20%max components".format(n0))

        # SECONDARY: normalized shapes agree within a wide margin
        rel = np.linalg.norm(gd - gap_static) / np.linalg.norm(gap_static)
        assert rel < 0.05, \
            ("dynamic (n0={}) central slice deviates {:.3%} from the static "
             "gap (>5%); k-structure not reproduced".format(n0, rel))

    # exact linearized-gap conjugation symmetry between the two central slices
    a = gap_dyn[0, 0, :, :, 0, central[0]].ravel()
    b = gap_dyn[0, 0, :, :, 0, central[1]].ravel()
    assert (np.linalg.norm(a - np.conj(b)) / np.linalg.norm(a)) < 1e-8, \
        "phi(k,-w_n) != conj phi(k,w_n) across the two central slices"


# ---------------------------------------------------------------------------
# Channel-parity selection for the dynamic (frequency-resolved) gap.
# Regression for: solve_dynamic's eigenvalue path returned ARPACK's raw
# largest-|lambda| WITHOUT the static path's channel-parity filter, so on real
# FLEX data (bp-ICl2 8.3 GPa UV) it reported an even-frequency, odd-(k,orbital)
# mode -- combined parity -1, i.e. the TRIPLET sector -- for a singlet request.
# ---------------------------------------------------------------------------

def _rand_gap(norb=2, Nx=1, Ny=4, Nz=1, nmat=4, seed=0):
    rng = np.random.default_rng(seed)
    shape = (norb, norb, Nx, Ny, Nz, nmat)
    return (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))


def test_reverse_kw_is_an_involution():
    from hwave.solver import eliashberg_dynamic as ed
    g = _rand_gap(seed=1)
    g2 = ed._reverse_kw_and_orbital(ed._reverse_kw_and_orbital(g))
    assert np.allclose(g2, g, atol=1e-12)


def test_parity_projectors_split_even_odd():
    from hwave.solver import eliashberg_dynamic as ed
    g = _rand_gap(seed=2)
    P = ed._reverse_kw_and_orbital
    ge = 0.5 * (g + P(g))   # combined-even  (P ge = ge)   -> singlet sector
    go = 0.5 * (g - P(g))   # combined-odd   (P go = -go)  -> triplet sector
    # even gap is fixed by the singlet projector, annihilated by triplet
    assert np.allclose(ed._project_parity_dynamic(ge, "singlet"), ge, atol=1e-12)
    assert np.linalg.norm(ed._project_parity_dynamic(ge, "triplet")) < 1e-10
    # odd gap is fixed by the triplet projector, annihilated by singlet
    assert np.allclose(ed._project_parity_dynamic(go, "triplet"), go, atol=1e-12)
    assert np.linalg.norm(ed._project_parity_dynamic(go, "singlet")) < 1e-10
    # parity labels
    assert ed._is_parity_dynamic(ge, "singlet")
    assert not ed._is_parity_dynamic(ge, "triplet")
    assert ed._is_parity_dynamic(go, "triplet")
    assert not ed._is_parity_dynamic(go, "singlet")


def test_reorder_promotes_channel_parity_eigenpair():
    """The odd mode has the larger eigenvalue but the singlet request must
    surface the even (combined-parity +1) eigenpair as leading."""
    from hwave.solver import eliashberg_dynamic as ed
    gap_shape = (2, 2, 1, 4, 1, 4)
    P = ed._reverse_kw_and_orbital
    g = _rand_gap(seed=3)
    ge = 0.5 * (g + P(g))              # even/singlet
    go = 0.5 * (g - P(g))              # odd/triplet
    vecs = np.column_stack([go.ravel(), ge.ravel()])   # odd first (col 0)
    vals = np.array([2.0 + 0j, 1.0 + 0j])              # odd has larger lambda
    rvals, rvecs, match = ed._reorder_eigenpairs_by_parity_dynamic(
        vals, vecs, gap_shape, "singlet")
    # even (lambda=1) promoted to leading; match flags reordered with it
    assert np.isclose(rvals[0].real, 1.0)
    assert bool(match[0]) is True
    assert bool(match[1]) is False
    assert np.allclose(rvecs[:, 0], ge.ravel(), atol=1e-12)


def test_solve_dynamic_eigenvalue_writes_match_column(tmp_path):
    """solve_dynamic (eigenvalue path) tags eigenvalue.dat with the
    channel-parity match column and returns an in-sector leading gap."""
    import os
    from hwave.solver import eliashberg_dynamic as ed
    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1, Nx=2, Ny=2, Nz=1)
    inp = {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": m["nmat"],
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic",
                       "pairing_type": "singlet",
                       "solver_mode": "eigenvalue", "num_eigenvalues": 6},
    }
    lam = ed.solve_dynamic(inp)
    assert np.isfinite(lam)
    with open(os.path.join(output_dir, "eigenvalue.dat")) as fh:
        lines = fh.read().splitlines()
    # wiring: the channel-parity match column is written
    assert any("match(1=channel-parity)" in ln for ln in lines)
    # parse the per-eigenvalue table (index Re Im |ev| match)
    rows = [ln.split() for ln in lines
            if ln and not ln.startswith("#") and len(ln.split()) == 5]
    assert rows, "expected a 5-column eigenvalue table with the match flag"
    matches = [int(r[4]) for r in rows]
    res = [float(r[1]) for r in rows]
    # selection: the reported leading lambda is the FIRST channel-parity row
    # if any exists (reordering surfaced it), else the raw leading row.
    if any(matches):
        first_match = matches.index(1)
        assert np.isclose(lam, res[first_match]), \
            "leading lambda is not the first channel-parity eigenvalue"
        # and the reordering put a matching pair at index 0
        assert matches[0] == 1
    else:
        # no in-sector mode on this (non-centrosymmetric) fixture: the driver
        # warns and falls back to the raw leading pair.
        assert np.isclose(lam, res[0])


def test_project_seed_dynamic_guard():
    """Seed projection raises when the seed has no component in the requested
    sector (mirrors sc._solve_iteration's guard), else returns a normalized
    in-sector seed."""
    from hwave.solver import eliashberg_dynamic as ed
    P = ed._reverse_kw_and_orbital
    g = _rand_gap(seed=7)
    ge = 0.5 * (g + P(g))   # even/singlet
    go = 0.5 * (g - P(g))   # odd/triplet
    # even seed survives the singlet projection, in-sector and non-zero, and is
    # returned WITHOUT renormalization (mirrors static _solve_iteration).
    s = ed._project_seed_dynamic(ge, "singlet")
    assert np.linalg.norm(s) > 0
    assert np.allclose(s, ge, atol=1e-12)   # ge is already even -> unchanged
    assert ed._is_parity_dynamic(s, "singlet")
    # a pure-odd seed has no singlet component -> raise
    import pytest
    with pytest.raises(ValueError, match="parity sector"):
        ed._project_seed_dynamic(go, "singlet")


def test_parity_leakage_zero_for_commuting_operator():
    """_parity_leakage ~ 0 for an operator that commutes with parity (identity)
    and O(1) for one that maps into the opposite sector."""
    from hwave.solver import eliashberg_dynamic as ed
    gap_shape = (2, 2, 1, 4, 1, 4)

    class _Op:
        def __init__(self, fn): self.fn = fn
        def matvec(self, x): return self.fn(x)

    ident = _Op(lambda x: x)
    assert ed._parity_leakage(ident, gap_shape, "singlet") < 1e-12

    # elementwise multiply by a P-ODD field d (P d = -d): it maps an even probe
    # to an odd image (and vice versa), i.e. the kernel does NOT commute with P,
    # so every probe leaks fully into the opposite sector.
    P = ed._reverse_kw_and_orbital
    r = _rand_gap(seed=13)
    d = 0.5 * (r - P(r))                 # odd projection -> P d = -d
    def mult(x):
        return (d * x.reshape(gap_shape)).ravel()
    assert ed._parity_leakage(_Op(mult), gap_shape, "singlet") > 0.9


# ---------------------------------------------------------------------------
# Eigenvector-continuation seed (v0) — track one physical branch across an
# exceptional point of the non-Hermitian kernel, where the algebraically
# largest eigenvalue can jump between a real and a complex branch.
# ---------------------------------------------------------------------------

def test_order_by_seed_overlap_picks_overlapping_branch():
    import hwave.sc as sc
    vecs = np.eye(4, dtype=complex)[:, :3]        # e0, e1, e2
    vals = np.array([2.0 + 0j, 0.5 + 0j, 1.0 + 0j])   # largest-real is e0 (2.0)
    seed = vecs[:, 1].copy()                       # overlaps e1 (eigenvalue 0.5)
    v, w = sc._order_by_seed_overlap(vals, vecs, seed)
    assert np.isclose(v[0], 0.5)                    # picks the seeded branch, not 2.0
    assert np.allclose(w[:, 0], vecs[:, 1])
    # a zero seed falls back to real-part ordering
    v2, _ = sc._order_by_seed_overlap(vals, vecs, np.zeros(4, dtype=complex))
    assert np.isclose(v2[0], 2.0)


def _dyn_input(tmp_path, extra_eli):
    import os
    input_dir = str(tmp_path / "input")
    output_dir = str(tmp_path / "output")
    os.makedirs(output_dir, exist_ok=True)
    _write_geom_transfer_coulomb(input_dir, norb=1)
    m = _write_flex_fixture(tmp_path / "output", nmat=8, norb=1, Nx=2, Ny=2, Nz=1)
    eli = {"chi0q_mode": "flex", "frequency": "dynamic",
           "pairing_type": "singlet", "solver_mode": "eigenvalue",
           "num_eigenvalues": 6}
    eli.update(extra_eli)
    return {
        "mode": {"param": {"T": 0.5, "CellShape": [2, 2, 1],
                           "SubShape": [1, 1, 1], "Nmat": m["nmat"],
                           "filling": 0.5}},
        "file": {"input": {"interaction": {
                    "path_to_input": input_dir,
                    "Geometry": "geom.dat", "Transfer": "transfer.dat",
                    "CoulombIntra": "coulombintra.dat"}},
                 "output": {"path_to_output": output_dir}},
        "eliashberg": eli}, output_dir


def test_solve_dynamic_seed_selects_overlapping_eigenpair(tmp_path):
    """Seeding with a saved gap makes the run report the eigenpair that
    overlaps the seed. Seeding from its own leading gap is self-consistent
    (returns the same lambda); seeding from a sub-leading eigenvector selects
    that branch instead of the algebraically-largest one."""
    import os
    from hwave.solver import eliashberg_dynamic as ed
    inp, out = _dyn_input(tmp_path, {})
    lam0 = ed.solve_dynamic(inp)
    gap_path = os.path.join(out, "gap_dynamic.npz")
    assert os.path.exists(gap_path)

    # self-seed: same leading branch -> same lambda
    inp2, out2 = _dyn_input(tmp_path / "b", {"seed_eigenvector": gap_path})
    lam_self = ed.solve_dynamic(inp2)
    assert np.isclose(lam_self, lam0, rtol=1e-6)

    # the run consumed the seed and its written gap overlaps it strongly
    g = np.load(os.path.join(out2, "gap_dynamic.npz"))["gap"].ravel()
    s = np.load(gap_path)["gap"].ravel()
    ov = abs(np.vdot(s / np.linalg.norm(s), g / np.linalg.norm(g)))
    assert ov > 0.99


def test_seed_gap_nmat_mismatch_raises(tmp_path):
    import os
    from hwave.solver import eliashberg_dynamic as ed
    # a seed gap with a different Nmat than the run
    bad = str(tmp_path / "bad_seed.npz")
    np.savez(bad, gap=np.zeros((1, 1, 2, 2, 1, 6), dtype=complex))  # Nmat=6
    inp, _ = _dyn_input(tmp_path, {"seed_eigenvector": bad})        # run Nmat=8
    with pytest.raises(ValueError, match="Nmat"):
        ed.solve_dynamic(inp)
