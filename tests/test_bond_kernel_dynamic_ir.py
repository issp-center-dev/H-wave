"""Task 7: the IR-axis dynamic bond kernel, the tail estimator, and the
IR round-trip / Hermiticity contracts.

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-
design.md, section 2 (IR-coordinate Hermiticity threshold), section 3.3
(instantaneous per-basis definitions + the tail-estimator pipeline) and
section 3.5 (the IR round-trip contract, steps 1-5, with its two frozen
inequalities). Every number asserted below is defined there; none is tuned
here.

TWO of the spec's frozen contracts are currently NOT met and are recorded as
strict xfails rather than silently weakened (spec section 2 / 3.5 amendment
gates -- see .superpowers/sdd/2026-07-28-dynamic-bond-channels-phaseA/
task-7-report.md for the full measurement tables):

* the ABSOLUTE round-trip budget ``E(N) <= 3*(tail_est_rel + eps_IR)``: the
  dominant IR-vs-uniform difference is the UNIFORM path's own O(1/Nmat)
  frequency-convolution discretization of the FLUCTUATION term, which that
  budget (instantaneous tail + IR representation error) does not model. The
  REFINEMENT half of the same contract -- ``E(2N) <= 0.6*E(N) + eps_IR^max``
  for two successive doublings -- passes, in both the action and the
  eigenvalue metric and in both parity sectors, which is what establishes
  that the two kernels are discretizations of the SAME operator.
* the IR-coordinate Hermiticity threshold ``max(10*eps_IR_roundtrip, 1e-10)``.

What DOES pin the implementation's correctness here, unconditionally:
``test_ir_kernel_equals_scalar_ir_kernel_at_B1`` (exact, 4e-16) plus the
refinement tests.

The IR path is only exercised when ``sparse-ir`` is installed (the same
import-safe skip guard as tests/test_eliashberg_ir.py: a module-level
``pytest.importorskip`` would turn a missing optional dependency into an
ERROR under ``unittest discover``, which CI runs).
"""
import numpy as np
import pytest

from hwave.solver import bond_channels as bc

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

BETA = 5.0


# ---------------------------------------------------------------------------
# Step 1/2: the tail estimator (spec S3.3), on an analytic profile
# ---------------------------------------------------------------------------
def test_tail_estimate_on_exact_inverse_square():
    """The estimator's pipeline is fixed by the spec; on an exact ``1/w^2``
    profile the fit is exact, so the returned ``tail_est`` must reproduce the
    literal window-complement sum it models."""
    beta, nmat = 5.0, 64
    n_t = np.arange(nmat) - nmat // 2
    w = (2 * n_t + 1) * np.pi / beta
    X = (1.0 / w ** 2)[None, None, None, None, :]   # (1,1,1,1,nmat) toy
    est, rel, bad = bc.tail_estimate(X, beta, nmat)
    exact = (1.0 / beta) * sum(
        1.0 / (((2 * m + 1) * np.pi / beta) ** 2)
        for m in range(nmat // 2, 64 * nmat))       # positive tail
    exact *= 2                                      # both signs
    assert not bad
    assert abs(est - exact) / exact < 0.05
    # tail_est_rel is the same quantity over the stored-window sum of the
    # SAME profile (spec S3.3 "Relative form"), so it is checkable in closed
    # form too -- and it must be small here (the window holds most of the
    # 1/w^2 weight).
    denom = (1.0 / beta) * float(np.sum(1.0 / w ** 2))
    assert abs(rel - est / denom) <= 1e-12 * abs(rel)
    assert 0.0 < rel < 0.05


def test_tail_estimate_flags_an_unfittable_profile():
    """The reliability flag (relative outer-shell fit residual > 0.2) is the
    spec's guard against reporting a tail number for data the ``|a|/w^2 +
    |b|/w^3`` model does not describe. A profile that oscillates between
    neighbouring frequencies (an X nowhere near its asymptotic regime) is the
    canonical case; a merely noisy but smooth-on-average profile is NOT
    flagged, which is the estimator behaving as specified (the flag is a
    model-adequacy test at the 0.2 relative-residual level, not a noise
    detector)."""
    beta, nmat = 5.0, 64
    n_t = np.arange(nmat) - nmat // 2
    X = (1.0 + 0.9 * (-1.0) ** n_t)[None, :]
    _, _, bad = bc.tail_estimate(X, beta, nmat)
    assert bad
    rng = np.random.default_rng(11)
    _, _, mild = bc.tail_estimate(rng.random((2, 3, nmat)) + 0.5, beta, nmat)
    assert not mild


def test_tail_estimate_reduces_over_k_and_orbital_pairs():
    """The fitted object is the single scalar profile ``x(n) = max_{k,pairs}
    |X|`` (spec S3.3 "Data fitted"): one fit, never per-k. Broadcasting the
    same 1/w^2 profile over extra leading axes with a k-dependent (<=1)
    modulation must therefore give exactly the max-reduced answer."""
    beta, nmat = 5.0, 32
    n_t = np.arange(nmat) - nmat // 2
    w = (2 * n_t + 1) * np.pi / beta
    prof = 1.0 / w ** 2
    scale = np.array([0.25, 1.0, 0.5])[:, None]
    est_multi, rel_multi, _ = bc.tail_estimate(
        (scale * prof[None, :])[None, :, :], beta, nmat)
    est_one, rel_one, _ = bc.tail_estimate(prof[None, :], beta, nmat)
    assert abs(est_multi - est_one) <= 1e-12 * est_one
    # the denominator is the max-reduced profile's window sum as well
    assert abs(rel_multi - rel_one) <= 1e-12 * rel_one


def test_tail_estimate_rejects_a_mismatched_window():
    with pytest.raises(ValueError, match="last axis"):
        bc.tail_estimate(np.ones((3, 16)), 5.0, 32)
    with pytest.raises(ValueError, match="even number"):
        bc.tail_estimate(np.ones((3, 33)), 5.0, 33)


# ---------------------------------------------------------------------------
# Steps 3-5: the IR kernel
# ---------------------------------------------------------------------------
ir_only = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                             reason="sparse-ir not installed")

WMAX = 10.0          # > the fixture's bandwidth (|e| <~ 2.5) with room to spare
IR_TOL = 1.0e-8
NX = NY = 4
NZ = 1
NMATS = (64, 128, 256)


def _dispersion(Nx=NX, Ny=NY, Nz=NZ, seed=20260728):
    """Real, inversion-symmetric e(k) -- the same construction as
    tests/test_oracle_bond_dynamic._symmetric_green, but with the dispersion
    pinned to a seed and hoisted OUT of the frequency grid: the refinement
    study compares Nmat, 2*Nmat and 4*Nmat, which is only a resolution study
    if the physical system is held fixed (``_symmetric_green`` draws from a
    shared module-level RNG, so calling it three times would give three
    different systems)."""
    rng = np.random.default_rng(seed)
    e = rng.normal(size=(Nx, Ny, Nz))
    return 0.5 * (e + np.roll(e[::-1, ::-1, ::-1], (1, 1, 1), (0, 1, 2)))


def _green_uniform(e, nmat, beta=BETA):
    """G = 1/(i w_n - e(k)): the v1 symmetry class at every Nmat, and a single
    pole per k, so the IR basis represents it to ~1e-9 (measured) and eps_IR
    never masks the quantities under test."""
    n_t = np.arange(nmat) - nmat // 2
    iw = 1j * (2 * n_t + 1) * np.pi / beta
    G = 1.0 / (iw[None, None, None, :] - e[..., None])
    return G[None, None, ...].astype(complex)


def _axes(beta=BETA, wmax=WMAX, tol=IR_TOL):
    """Brief: reuse eliashberg_dynamic's axis-construction helper rather than
    building IRAxis by hand. ``ir_wmax`` is given explicitly, so the (hr,
    inter_k) auto-estimate arguments are unused."""
    from hwave.solver import eliashberg_dynamic as ed
    return ed._ir_axes_for_run({"ir_tol": tol, "ir_wmax": wmax},
                               beta, None, None, 1)


def _compress3(arr, ax, nmat, label):
    """IR-compress an array whose frequency axis is at position 3 (the
    ``(Nx,Ny,Nz,nfreq,ND,ND)`` bond layout); every IRAxis transform contracts
    the LAST axis."""
    from hwave.solver import eliashberg_dynamic as ed
    return np.moveaxis(ed._ir_compress(np.moveaxis(arr, 3, -1), ax, nmat,
                                       label), -1, 3)


def _densify3(arr, ax, nmat):
    return np.moveaxis(ax.eval_to_uniform(
        ax.fit_from_freq(np.moveaxis(arr, 3, -1)), nmat), -1, 3)


_FIXTURES = {}


def _fixture(nmat, pairing, beta=BETA):
    """Uniform and IR builds of the SAME U=4 / NN-V=1 single-band model,
    split into the two operator parts.

    Both kernels are fed the IR-densified chi and green, so the comparison
    isolates the KERNEL ALGEBRA (frequency convolution + the two per-basis
    instantaneous definitions) from the input data's out-of-basis content --
    the same discipline as the scalar gate
    ``tests/test_eliashberg_ir.py::test_ir_matvec_matches_uniform_kernel``:
    the raw uniform-FFT susceptibility additionally carries finite-Nmat
    artifacts (the delta(tau) constant, aliasing images) that no basis of
    bandwidth wmax can represent and that spec S3.5's error budget does not
    model.

    The two parts are built by zeroing the other part's vertex input (Vpp = 0
    -> fluctuation only; chi = 0 -> instantaneous only) so both are produced
    by the SAME production code path under test.
    """
    key = (nmat, pairing, beta)
    if key in _FIXTURES:
        return _FIXTURES[key]
    from hwave.solver import eliashberg_dynamic as ed
    from tests.test_bond_channels_dynamic import _nn_bond_set
    from tests.test_bond_dynamic_hermiticity import _build_vertices

    axF, axB = _axes(beta)
    green_u = _green_uniform(_dispersion(), nmat, beta)
    # eps_IR (spec S3.3, single definition): the green's IR round-trip
    # residual on the run's fermionic axis, full tensor, measured on the
    # run's physical (uniform) green.
    eps_ir = float(np.linalg.norm(
        axF.eval_to_uniform(axF.fit_from_uniform(green_u, nmat), nmat)
        - green_u) / np.linalg.norm(green_u))
    green_n = ed._ir_compress(green_u, axF, nmat, "green")
    green_u = axF.eval_to_uniform(axF.fit_from_freq(green_n), nmat)

    bset = _nn_bond_set()
    S_bond, C_bond, Vpp_s, Vpp_t = _build_vertices(bset, U=4.0,
                                                   shape=(NX, NY, NZ))
    chi_bar = bc.bond_bubble_dynamic(green_u, bset, beta)
    chi_s_u, chi_c_u, _, _ = bc.dress_bond_dynamic(chi_bar, S_bond, C_bond)
    chi_s_n = _compress3(chi_s_u, axB, nmat, "chi_s")
    chi_c_n = _compress3(chi_c_u, axB, nmat, "chi_c")
    chi_s_u = _densify3(chi_s_n, axB, nmat)
    chi_c_u = _densify3(chi_c_n, axB, nmat)

    zv = np.zeros_like(Vpp_s)
    zc_u = np.zeros_like(chi_s_u)
    zc_n = np.zeros_like(chi_s_n)

    def uni(cs, cc, vs, vt):
        return bc.make_bond_kernel_dynamic(cs, cc, S_bond, C_bond, vs, vt,
                                           green_u, bset, pairing, beta)[0]

    def ir(cs, cc, vs, vt):
        return bc.make_bond_kernel_dynamic_ir(cs, cc, S_bond, C_bond, vs, vt,
                                              green_n, bset, pairing, beta,
                                              axF, axB)[0]

    fx = dict(
        axF=axF, axB=axB, beta=beta, nmat=nmat, eps_ir=eps_ir,
        green_u=green_u, green_n=green_n,
        u_full=uni(chi_s_u, chi_c_u, Vpp_s, Vpp_t),
        u_fl=uni(chi_s_u, chi_c_u, zv, zv),
        u_pp=uni(zc_u, zc_u, Vpp_s, Vpp_t),
        i_full=ir(chi_s_n, chi_c_n, Vpp_s, Vpp_t),
        i_fl=ir(chi_s_n, chi_c_n, zv, zv),
        i_pp=ir(zc_n, zc_n, Vpp_s, Vpp_t))
    _FIXTURES[key] = fx
    return fx


def _gap_coeffs(axF, seed=5, npole=6, wpole=2.0, beta=BETA):
    """Spec S3.5 step 1: a probe gap that is IR-representable BY CONSTRUCTION
    (IR coefficients, densified to whatever grid is asked for), drawn with the
    PHYSICAL coefficient weight -- a random Lehmann sum with poles inside
    ``|w| <= wpole``, projected onto the basis.

    Why not white (unit-variance) IR coefficients: measured on this fixture,
    a white draw saturates the basis bandwidth, and the pointwise product
    ``X = G2 * phi`` (a tau-convolution, so the spectral supports ADD) then
    falls outside it. The IR equal-time value ``sum_l X_l u_l(0+)`` is then
    9.7% off the converged window sum and does not improve with Nmat
    (measured 9.70e-2, 9.70e-2, 9.70e-2, 9.76e-2 at Nmat = 64, 128, 256,
    4096), i.e. the comparison would measure the basis bandwidth rather than
    the kernel. With the physical draw the same quantity agrees to 2.6e-6.
    """
    rng = np.random.default_rng(seed)
    w = rng.uniform(-wpole, wpole, size=npole)
    amp = (rng.normal(size=(1, 1, NX, NY, NZ, npole))
           + 1j * rng.normal(size=(1, 1, NX, NY, NZ, npole)))
    iw = 1j * axF.freq_n * np.pi / beta
    return axF.fit_from_freq((amp[..., None] / (iw[None, :] - w[:, None]))
                             .sum(axis=-2))


def _apply(op, x):
    return (op @ x.ravel()).reshape(x.shape)


def _round_trip_action(fx, coeffs):
    """Spec S3.5 steps 2-4, with the PART-AWARE densification of step 3.

    ``densify`` on the IR path means "evaluate the IR-represented object on
    the uniform grid". The operator's output is not a pure IR object: it is
    (IR-representable fluctuation part) + (frequency-FLAT instantaneous part),
    and a frequency-independent constant is exactly what the IR basis cannot
    represent (u_l(i w) ~ 1/w). Fitting the sum destroys the flat part -- so
    the flat part is densified by broadcast and only the fluctuation part is
    fitted. ``test_ir_naive_densification_is_not_a_measurement`` pins that the
    naive alternative is a measurement artifact, not a stricter test.
    """
    axF, nmat = fx["axF"], fx["nmat"]
    phi_d = axF.eval_to_uniform(coeffs, nmat)
    phi_n = axF.eval_to_freq(coeffs)
    Y_uni = _apply(fx["u_full"], phi_d)
    Y_fl = axF.eval_to_uniform(
        axF.fit_from_freq(_apply(fx["i_fl"], phi_n)), nmat)
    Y_pp = _apply(fx["i_pp"], phi_n)
    Y_ir = Y_fl + np.broadcast_to(Y_pp[..., :1], Y_uni.shape)
    E = float(np.linalg.norm(Y_uni - Y_ir) / np.linalg.norm(Y_uni))
    X = _pair_amplitude(fx, phi_d)
    est, rel, bad = bc.tail_estimate(X, fx["beta"], nmat)
    return E, rel, bad


def _pair_amplitude(fx, phi_dense):
    """X_{l3l4}(k', i w_n') = sum_{ef} G2 phi -- the object whose window
    truncation ``tail_estimate`` models (norb = 1 here)."""
    from hwave.solver import eliashberg_dynamic as ed
    G2 = ed.calc_g2_dynamic(fx["green_u"], fx["beta"])
    return G2[0, 0, 0, 0] * phi_dense[0, 0]


def _sector(A, shape, pairing):
    from scipy.sparse.linalg import LinearOperator
    from hwave.solver import eliashberg_dynamic as ed
    n = int(np.prod(shape))

    def mv(v):
        x = ed._project_parity_dynamic(np.asarray(v).reshape(shape), pairing)
        return ed._project_parity_dynamic(
            _apply(A, x), pairing).ravel()
    return LinearOperator((n, n), matvec=mv, dtype=complex)


def _sector_lambda(fx, pairing, which):
    """Leading eigenvalue INSIDE the parity sector (spec S3.5: branches are
    matched by sector, so no eigenpair-sorting rule is needed)."""
    from scipy.sparse.linalg import eigs
    from hwave.solver import eliashberg_dynamic as ed
    axF, nmat = fx["axF"], fx["nmat"]
    c = _gap_coeffs(axF)
    if which == "uniform":
        shape = (1, 1, NX, NY, NZ, nmat)
        v0 = axF.eval_to_uniform(c, nmat)
        op = fx["u_full"]
    else:
        shape = (1, 1, NX, NY, NZ, axF.n_freq)
        v0 = axF.eval_to_freq(c)
        op = fx["i_full"]
    v0 = ed._project_parity_dynamic(v0, pairing).ravel()
    return complex(eigs(_sector(op, shape, pairing), k=1, which="LM", v0=v0,
                        tol=1e-10, return_eigenvectors=False)[0])


def _dense(op, shape):
    size = int(np.prod(shape))
    M = np.zeros((size, size), dtype=complex)
    for i in range(size):
        v = np.zeros(size, dtype=complex)
        v[i] = 1.0
        M[:, i] = op @ v
    return M


def _whitened_residual(K, w):
    sq = np.sqrt(w).ravel()
    Kt = (sq[:, None] * K) / sq[None, :]
    return float(np.linalg.norm(Kt - Kt.conj().T) / np.linalg.norm(Kt))


def _pair_weight(green):
    """w(k, i w_n) = G(k,iw) G(-k,-iw), real and positive in the v1 symmetry
    class (spec section 2). The (-k,-iw) map is the same reverse+roll the
    uniform P0-1 gate uses; it is valid on IR nodes because the node set is
    exactly symmetric under integer negation."""
    wc = (green * np.roll(green[:, :, ::-1, ::-1, ::-1, ::-1],
                          (1, 1, 1), (2, 3, 4)))[0, 0]
    assert np.abs(wc.imag).max() <= 1e-10 * np.abs(wc).max()
    w = wc.real
    assert w.min() >= 1e-12 * w.max()
    return w


# --- correctness: the unconditional pin -------------------------------------
@ir_only
def test_ir_kernel_equals_scalar_ir_kernel_at_B1():
    """B = 1 reduction (spec S3.5): with a single (on-site) bond channel the
    bond kernel must BE the existing scalar IR kernel -- same transport, same
    normalization, same instantaneous convention. Exact equality, so it pins
    the implementation independently of every discretization budget below."""
    from hwave.solver import eliashberg_dynamic as ed
    nmat = 64
    axF, axB = _axes()
    green_u = _green_uniform(_dispersion(), nmat)
    green_n = ed._ir_compress(green_u, axF, nmat, "green")
    bset = bc.resolve_interactions({((0, 0, 0), (0, 0)): 0.0}, np.eye(3), 1)
    assert bset.n_channels == 1
    rng = np.random.default_rng(3)

    # An IR-representable bosonic vertex, fed through the bond algebra as
    # F = 1.5 * S chi_s S with S = 1 (and C = 0), so the bond kernel's
    # fluctuation vertex equals the scalar kernel's pairing vertex.
    cb = (rng.normal(size=(NX, NY, NZ, axB.L))
          + 1j * rng.normal(size=(NX, NY, NZ, axB.L)))
    Vs_n = axB.eval_to_freq(cb)
    S_bond = np.ones((NX, NY, NZ, 1, 1), dtype=complex)
    C_bond = np.zeros((NX, NY, NZ, 1, 1), dtype=complex)
    chi_s = Vs_n[..., None, None] / 1.5
    zero = np.zeros((1, 1), dtype=complex)
    A, vec_size = bc.make_bond_kernel_dynamic_ir(
        chi_s, np.zeros_like(chi_s), S_bond, C_bond, zero, zero, green_n,
        bset, "singlet", BETA, axF, axB)
    assert vec_size == NX * NY * NZ * axF.n_freq

    cphi = (rng.normal(size=(1, 1, NX, NY, NZ, axF.L))
            + 1j * rng.normal(size=(1, 1, NX, NY, NZ, axF.L)))
    phi_n = axF.eval_to_freq(cphi)
    Y_bond = _apply(A, phi_n)
    Y_scalar = ed.eliashberg_kernel_ir(
        ed._ir_vertex_to_rtau(Vs_n[None, None, None, None], axB, axF),
        ed.calc_g2_dynamic(green_n, BETA), phi_n.copy(), axF, BETA)
    np.testing.assert_allclose(Y_bond, Y_scalar, rtol=1e-11, atol=1e-13)


@ir_only
def test_ir_kernel_rejects_uniform_grid_inputs():
    """The input contract (docstring): arrays must already be on the IR
    sampling nodes. A uniform-grid array of the same rank must fail loudly
    rather than be silently misread as node values."""
    from hwave.solver import eliashberg_dynamic as ed
    from tests.test_bond_channels_dynamic import _nn_bond_set
    from tests.test_bond_dynamic_hermiticity import _build_vertices
    nmat = 64
    axF, axB = _axes()
    bset = _nn_bond_set()
    S, C, Vs, Vt = _build_vertices(bset, U=4.0, shape=(NX, NY, NZ))
    green_u = _green_uniform(_dispersion(), nmat)
    green_n = ed._ir_compress(green_u, axF, nmat, "green")
    ND = bset.n_channels
    chi_u = np.zeros((NX, NY, NZ, nmat, ND, ND), dtype=complex)
    chi_n = np.zeros((NX, NY, NZ, axB.n_freq, ND, ND), dtype=complex)
    with pytest.raises(ValueError, match="BOSONIC IR sampling"):
        bc.make_bond_kernel_dynamic_ir(chi_u, chi_u, S, C, Vs, Vt, green_n,
                                       bset, "singlet", BETA, axF, axB)
    with pytest.raises(ValueError, match="ALREADY on the IR nodes"):
        bc.make_bond_kernel_dynamic_ir(chi_n, chi_n, S, C, Vs, Vt, green_u,
                                       bset, "singlet", BETA, axF, axB)
    with pytest.raises(ValueError, match="beta"):
        bc.make_bond_kernel_dynamic_ir(chi_n, chi_n, S, C, Vs, Vt, green_n,
                                       bset, "singlet", 2.0 * BETA, axF, axB)


# --- spec S3.5: the round-trip contract -------------------------------------
@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_kernel_round_trip_refinement(pairing, capsys):
    """Spec S3.5 step 5, REFINEMENT half (frozen):
    ``E(2N) <= 0.6*E(N) + eps_IR^max`` for two successive doublings.

    This is the inequality that establishes the two kernels are
    discretizations of the same operator; the absolute half is the xfail
    below."""
    axF = _axes()[0]
    c = _gap_coeffs(axF)
    E, tails, eps = {}, {}, {}
    for nmat in NMATS:
        fx = _fixture(nmat, pairing)
        E[nmat], tails[nmat], bad = _round_trip_action(fx, c)
        eps[nmat] = fx["eps_ir"]
        assert not bad, (
            "tail estimate unreliable at Nmat={} -- spec S3.3 forbids "
            "recording that as an inconclusive pass".format(nmat))
    with capsys.disabled():
        print("\n  [{}] E(Nmat) = {}\n        tail_est_rel = {}\n"
              "        eps_IR = {}".format(
                  pairing, {k: "%.4e" % v for k, v in E.items()},
                  {k: "%.4e" % v for k, v in tails.items()},
                  {k: "%.3e" % v for k, v in eps.items()}))
    for lo, hi in ((64, 128), (128, 256)):
        eps_max = max(eps[lo], eps[hi])
        assert E[hi] <= 0.6 * E[lo] + eps_max, (
            "no refinement {} -> {}: E = {:.4e} -> {:.4e} (budget {:.4e})"
            .format(lo, hi, E[lo], E[hi], 0.6 * E[lo] + eps_max))


@ir_only
@pytest.mark.xfail(strict=True, reason=(
    "spec S3.5 step 5 ABSOLUTE budget amendment gate: measured E is 28-150x "
    "3*(tail_est_rel + eps_IR) at every resolution, because the dominant "
    "IR-vs-uniform difference is the UNIFORM path's own O(1/Nmat) "
    "frequency-convolution discretization of the FLUCTUATION term -- a term "
    "the budget (instantaneous tail + IR representation error) does not "
    "model. The refinement half of the same contract passes. Do not weaken "
    "the inequality here; amend the spec (task-7-report.md)."))
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_kernel_round_trip_absolute_budget(pairing):
    """Spec S3.5 step 5, ABSOLUTE half (frozen):
    ``E(N) <= 3*(tail_est_rel(N) + eps_IR)`` at each resolution."""
    axF = _axes()[0]
    c = _gap_coeffs(axF)
    for nmat in NMATS:
        fx = _fixture(nmat, pairing)
        E, rel, _ = _round_trip_action(fx, c)
        assert E <= 3.0 * (rel + fx["eps_ir"]), (
            "E({}) = {:.4e} exceeds 3*(tail_est_rel {:.4e} + eps_IR {:.3e})"
            .format(nmat, E, rel, fx["eps_ir"]))


@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_kernel_eigenvalue_refinement(pairing, capsys):
    """The end-to-end EIGENVALUE form of the contract (spec S3.5): an action
    norm does not define an eigenvalue metric. Per sector,
    ``e(N) = |lam_uni(N) - lam_IR| / |lam_uni(N)|`` with the same frozen
    refinement inequality. No densification enters here, so this metric is
    also free of the frequency-flat-component question."""
    axF = _axes()[0]
    c = _gap_coeffs(axF)
    e_rel, tails, eps = {}, {}, {}
    for nmat in NMATS:
        fx = _fixture(nmat, pairing)
        lam_u = _sector_lambda(fx, pairing, "uniform")
        lam_i = _sector_lambda(fx, pairing, "ir")
        e_rel[nmat] = float(abs(lam_u - lam_i) / abs(lam_u))
        tails[nmat] = _round_trip_action(fx, c)[1]
        eps[nmat] = fx["eps_ir"]
    with capsys.disabled():
        print("\n  [{}] e(Nmat) = {}".format(
            pairing, {k: "%.4e" % v for k, v in e_rel.items()}))
    for lo, hi in ((64, 128), (128, 256)):
        eps_max = max(eps[lo], eps[hi])
        assert e_rel[hi] <= 0.6 * e_rel[lo] + eps_max, (
            "no eigenvalue refinement {} -> {}: e = {:.4e} -> {:.4e}"
            .format(lo, hi, e_rel[lo], e_rel[hi]))


@ir_only
@pytest.mark.xfail(strict=True, reason=(
    "same spec S3.5 absolute-budget amendment gate as the action metric: "
    "measured e(N) is 30-90x 3*(tail_est_rel + eps_IR^max) while the "
    "refinement half passes (task-7-report.md)."))
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_kernel_eigenvalue_absolute_budget(pairing):
    """Spec S3.5, eigenvalue form, ABSOLUTE half:
    ``e(N) <= 3*(tail_est_rel(N) + eps_IR^max)``."""
    axF = _axes()[0]
    c = _gap_coeffs(axF)
    eps_max = max(_fixture(n, pairing)["eps_ir"] for n in NMATS)
    for nmat in NMATS:
        fx = _fixture(nmat, pairing)
        lam_u = _sector_lambda(fx, pairing, "uniform")
        lam_i = _sector_lambda(fx, pairing, "ir")
        e_rel = float(abs(lam_u - lam_i) / abs(lam_u))
        rel = _round_trip_action(fx, c)[1]
        assert e_rel <= 3.0 * (rel + eps_max), (
            "e({}) = {:.4e} exceeds 3*(tail_est_rel {:.4e} + eps_IR^max "
            "{:.3e})".format(nmat, e_rel, rel, eps_max))


@ir_only
def test_ir_naive_densification_is_not_a_measurement():
    """Evidence for the part-aware densification of ``_round_trip_action``:
    fitting the FULL IR output (fluctuation + frequency-flat instantaneous
    term) to the IR basis and evaluating on the uniform grid makes the error
    GROW with Nmat, because a constant is not in the span of the IR basis and
    the uniform grid keeps extending to frequencies where the fit has decayed
    to zero while the true answer is still the constant. A monotonically
    growing "error" measures the densification, not the kernel."""
    axF = _axes()[0]
    c = _gap_coeffs(axF)
    naive = []
    for nmat in NMATS:
        fx = _fixture(nmat, "singlet")
        phi_d = axF.eval_to_uniform(c, nmat)
        phi_n = axF.eval_to_freq(c)
        Y_uni = _apply(fx["u_full"], phi_d)
        Y_naive = axF.eval_to_uniform(
            axF.fit_from_freq(_apply(fx["i_full"], phi_n)), nmat)
        naive.append(float(np.linalg.norm(Y_uni - Y_naive)
                           / np.linalg.norm(Y_uni)))
        # ... while the instantaneous part IS exactly flat on both axes, so
        # the broadcast densification is exact by construction.
        for key in ("u_pp", "i_pp"):
            Y = _apply(fx[key], phi_d if key.startswith("u") else phi_n)
            assert (np.abs(Y - Y[..., :1]).max()
                    <= 1e-12 * np.abs(Y).max()), key
    assert naive[0] < naive[1] < naive[2], naive
    # and the part-aware metric goes the other way on the same data
    part_aware = [_round_trip_action(_fixture(n, "singlet"), c)[0]
                  for n in NMATS]
    assert part_aware[0] > part_aware[1] > part_aware[2], part_aware


# --- spec section 2: the IR-coordinate Hermiticity gate ---------------------
@ir_only
@pytest.mark.xfail(strict=True, reason=(
    "spec section 2 amendment gate: the IR-coordinate residual is ~1.0, far "
    "above the frozen max(10*eps_IR_roundtrip, 1e-10) (eps_IR_roundtrip on a "
    "single-pole green is ~1e-16 on the nodes, so the threshold degenerates "
    "to 1e-10). Decomposition (Nmat=64, singlet): instantaneous part 1.31 "
    "(the IR equal-time functional and the constant ket are adjoint in the "
    "physical L^2(0,beta) metric, NOT in the Euclidean sampling metric), "
    "fluctuation part 0.236 shrinking to 0.045 at Nmat=128 (input-resolution "
    "limited; the uniform kernel through the same IR projection measures "
    "0.019). Do not tune the threshold -- amend the spec "
    "(task-7-report.md)."))
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_coordinate_hermiticity_threshold(pairing, capsys):
    """Spec section 2, measurement-coordinates paragraph: the IR operator acts
    on IR SAMPLING VALUES, and IR fit/evaluation is not unitary in the
    Euclidean sampling metric, so the sampling-coordinate residual is reported
    separately from the (machine-precision) uniform one, against the FROZEN
    threshold ``max(10 * eps_IR_roundtrip, 1e-10)``."""
    nmat = 64
    fx = _fixture(nmat, pairing)
    axF, green_n = fx["axF"], fx["green_n"]

    # eps_IR_roundtrip on the axis the operator actually lives on.
    g_rt = axF.eval_to_freq(axF.fit_from_freq(green_n))
    eps_rt = float(np.linalg.norm(g_rt - green_n) / np.linalg.norm(green_n))
    threshold = max(10.0 * eps_rt, 1e-10)

    shape = (1, 1, NX, NY, NZ, axF.n_freq)
    w_n = _pair_weight(green_n)
    res = _whitened_residual(_dense(fx["i_full"], shape), w_n)
    parts = {k: _whitened_residual(_dense(fx[k], shape), w_n)
             for k in ("i_fl", "i_pp")}
    w_u = _pair_weight(fx["green_u"])
    shape_u = (1, 1, NX, NY, NZ, nmat)
    uni = _whitened_residual(_dense(fx["u_full"], shape_u), w_u)
    with capsys.disabled():
        print("\n  [{}] IR-coordinate Hermiticity residual = {:.4e} "
              "(fluctuation {:.4e}, instantaneous {:.4e}); uniform "
              "coordinates {:.4e}; eps_IR_roundtrip = {:.3e}; threshold = "
              "{:.3e}".format(pairing, res, parts["i_fl"], parts["i_pp"],
                              uni, eps_rt, threshold))
    assert res <= threshold, (
        "IR-COORDINATE HERMITICITY GATE FAILED (spec section 2 amendment "
        "gate): residual {:.4e} > max(10*eps_IR_roundtrip {:.3e}, 1e-10) = "
        "{:.3e}".format(res, eps_rt, threshold))
