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


def test_static_limit_matvec_multi_k():
    """PRIMARY convention pin (independent of the dynamic kernel's own
    self-consistency): on a NON-symmetric mesh (Nx=4, Ny=2) the frequency-flat
    limit of ``eliashberg_kernel_dynamic`` must reproduce, component-by-
    component, the STATIC kernel ``sc._eliashberg_kernel_fft`` evaluated with
    ``G2_static = sum_n G2_w``. This pins the spatial fold and the -(1/N)
    normalization against the static reference, not against a self-extraction.
    Nk=1 would hide both the fold and the 1/N, so the mesh must be multi-k and
    non-symmetric (q = k - k' is then not self-symmetric)."""
    import hwave.sc as sc
    from hwave.solver import eliashberg_dynamic as ed
    norb, Nx, Ny, Nz, nmat = 1, 4, 2, 1, 4
    beta = 7.0  # unused by the kernel; G2 is passed raw to both paths
    rng = np.random.default_rng(20)

    def rc(shape):
        return rng.standard_normal(shape) + 1j * rng.standard_normal(shape)

    # frequency-flat vertex: build one bosonic slice, broadcast across nmat
    V0 = rc((norb, norb, norb, norb, Nx, Ny, Nz))            # (l1,l2,l3,l4,q)
    V_w = np.broadcast_to(V0[..., np.newaxis],
                          V0.shape + (nmat,)).copy()
    # frequency-resolved pair bubble (varies with n on purpose)
    G2_w = rc((norb, norb, norb, norb, Nx, Ny, Nz, nmat))
    # frequency-flat trial gap
    phi0 = rc((norb, norb, Nx, Ny, Nz))                      # (l1,l2,k)
    phi_w = np.broadcast_to(phi0[..., np.newaxis],
                            phi0.shape + (nmat,)).copy()

    out = ed.eliashberg_kernel_dynamic(V_w, G2_w, phi_w, norb, beta)

    # (i) the output must be frequency-flat
    assert np.allclose(out, out[..., :1], atol=1e-10), \
        "frequency-flat input did not give a frequency-flat output"

    # (ii) each frequency slice equals the static matvec with G2_static
    G2_static = G2_w.sum(axis=-1)                            # sum over n
    static_out = sc._eliashberg_kernel_fft(V0, G2_static, phi0, norb)
    for n in range(nmat):
        assert np.allclose(out[..., n], static_out, atol=1e-10), \
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
    for key in ("iomega", "T", "pairing_type", "frequency",
                "eigenvalue", "axis_order", "normalization"):
        assert key in npz.files, "missing gap_dynamic.npz key: {}".format(key)
    assert str(npz["frequency"]) == "dynamic"
    assert npz["iomega"].shape == (m["nmat"],)
    assert npz["gap"].shape == (1, 1, 2, 2, 1, m["nmat"])

    with open(os.path.join(output_dir, "gap.dat")) as fh:
        first = fh.readline()
    assert first.startswith("#")
    assert "frequency=dynamic" in first
