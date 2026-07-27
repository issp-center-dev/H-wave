"""Task 7: the IR-axis dynamic bond kernel, the tail estimator, and the
IR validation contracts.

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-
design.md, section 2 ("Measurement coordinates", amended 2026-07-28), section
3.3 (instantaneous per-basis definitions + the tail-estimator pipeline) and
section 3.5 ("IR validation", amended 2026-07-28). Every number asserted below
is defined there; none is tuned here.

The amendment was written on the measurements this file makes, so the tests
map one-to-one onto its clauses:

* step 1 -- the normative ANALYTIC two-pole gap generator (seed 20260728), so
  every resolution samples the same underlying function
  (``test_normative_gap_generator_is_analytic_and_ir_representable``), with
  the rejected white-coefficient probe kept as evidence
  (``test_white_ir_coefficients_are_not_a_physical_probe``);
* step 3 -- PART-AWARE densification, with the naive alternative pinned as a
  non-measurement (``test_ir_naive_densification_is_not_a_measurement``);
* step 5 -- three action clauses (instantaneous absolute budget with its
  vacuity guard; refinement; uniform anchor with ``D_uni``) and their two
  eigenvalue twins (with ``d_uni``);
* section 2 -- the Euclidean sampling-coordinate residual is a DIAGNOSTIC,
  never gated; the physical Hermiticity evidence for the IR branch is the
  B = 1 reduction against the oracle-verified scalar ``eliashberg_kernel_ir``
  (measured 4.4e-16).

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
# The IR kernel and the AMENDED (2026-07-28) IR validation contracts
# ---------------------------------------------------------------------------
ir_only = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                             reason="sparse-ir not installed")

WMAX = 10.0          # > the fixture's bandwidth (|e| <~ 2.5) with room to spare
# The E_pp clause compares the two instantaneous conventions against a tail
# that falls like 1/Nmat^2, so the basis must represent X's equal-time value
# BETTER than that tail at the largest tested resolution -- otherwise the
# clause measures basis truncation instead. Measured E_pp/budget ratios:
# ir_tol = 1e-8 (L=21) leaves a ~1e-6 floor and the ratio degrades with
# resolution (0.03, 0.09, 0.92 singlet; 0.39, 0.48, 1.95 triplet -- the last
# one fails); at 1e-10 (L=24) the floor is gone and the ratio is
# resolution-INDEPENDENT (0.04 / 0.38 at all three), and 1e-12 (L=27)
# reproduces those same ratios, i.e. the fixture is converged in the basis.
IR_TOL = 1.0e-10
NX = NY = 4
NZ = 1
NMATS = (64, 128, 256)          # the tested resolution set of spec S3.5
# D_uni(N)/d_uni(N) compare the UNIFORM kernel at N against itself at 2N on the
# shared window, so the largest tested resolution needs its own doubling.
ALL_NMATS = NMATS + (512,)


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


# --- the normative test gap (spec S3.5 step 1, amended) ---------------------
def _gap_poles(Nx=NX, Ny=NY, Nz=NZ, seed=20260728):
    """The two real, inversion-symmetrized pole arrays of the normative
    generator, drawn ONCE from ``default_rng(20260728)``."""
    rng = np.random.default_rng(seed)

    def draw():
        e = rng.normal(size=(Nx, Ny, Nz))
        return 0.5 * (e + np.roll(e[::-1, ::-1, ::-1], (1, 1, 1), (0, 1, 2)))
    return draw(), draw()


def _pwave_factor(Nx=NX, Ny=NY, Nz=NZ):
    """``sin(k_x)`` on the FFT grid: real, odd under the grid inversion
    ``i -> (N-i) % N``, and exactly the harmonic the triplet Cooper vertex
    couples to (``Vpp_t``'s only non-zero block is the antisymmetric
    ``(-x, +x)`` channel pair, measured)."""
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    return np.broadcast_to(np.sin(kx)[:, None, None], (Nx, Ny, Nz)).copy()


def _gap_analytic(nmat, pairing="singlet", beta=BETA):
    """``phi(k, i w) = 1/(i w - e1(k)) - 1/(i w - e2(k))``: an ANALYTIC
    function of omega, so the Nmat / 2*Nmat / 4*Nmat grids sample the same
    underlying object and the refinement inequalities compare resolutions of
    one function rather than independent random vectors (spec S3.5 step 1,
    amended 2026-07-28 -- white IR coefficients are not usable, see
    ``test_white_ir_coefficients_are_not_a_physical_probe``).

    TRIPLET sector: the generator's two pole arrays are inversion-symmetrized,
    so the raw gap is EVEN in k; with an inversion-even green the triplet
    (odd) parity sector is then realized as an ODD-FREQUENCY, even-k gap,
    whose equal-time amplitude ``sum_n X(k, i w_n)`` vanishes identically --
    measured ``||Y^pp_uni||/||Y_uni|| = 9.6e-17``, i.e. exactly the "gap that
    nulls the Cooper term" that spec S3.5's vacuity guard declares an invalid
    fixture for the E_pp clause. The lattice's own p-wave factor ``sin(k_x)``
    is therefore attached in that sector (only there): it is frequency-
    independent, so analyticity in omega, the resolution-independence of the
    generator and its IR representability are all untouched, and it puts the
    probe in the harmonic the triplet Cooper vertex actually couples to.
    """
    e1, e2 = _gap_poles()
    n_t = np.arange(nmat) - nmat // 2
    iw = 1j * (2 * n_t + 1) * np.pi / beta
    phi = (1.0 / (iw[None, None, None, :] - e1[..., None])
           - 1.0 / (iw[None, None, None, :] - e2[..., None]))
    if pairing == "triplet":
        phi = phi * _pwave_factor()[..., None]
    return phi[None, None, ...].astype(complex)


def _gap_projected(axF, nmat, pairing, beta=BETA):
    """Step 1's "then fit -> densify per resolution": the probe actually fed
    to the two kernels is the basis projection of the analytic gap, so the
    comparison carries no out-of-subspace ambiguity.

    The analytic gap is first projected onto the CHANNEL'S combined-parity
    sector (``Delta_{ab}(k,iw) -> Delta_{ba}(-k,-iw)``, even for singlet, odd
    for triplet) -- the subspace the solver actually works in, and the one the
    eigenvalue metric is defined on. This is required, not cosmetic: both
    pole arrays are inversion-symmetrized, so the raw generator is EVEN in k
    and its overlap with the odd-harmonic triplet Cooper vertex vanishes
    identically (measured ||Y^pp_uni||/||Y_uni|| = 9.6e-17), which spec S3.5's
    own vacuity guard declares an invalid fixture for the E_pp clause. The
    projection preserves step 1's resolution-independence: the partner
    frequency of ``n~`` is ``-n~-1``, which lies in the stored window of every
    grid, so all resolutions still sample one analytic function.

    Returns ``(phi_dense, phi_nodes, fit_residual_rel)``.
    """
    from hwave.solver import eliashberg_dynamic as ed
    phi = ed._project_parity_dynamic(
        _gap_analytic(nmat, pairing, beta), pairing)
    coeffs = axF.fit_from_uniform(phi, nmat)
    dense = axF.eval_to_uniform(coeffs, nmat)
    resid = float(np.linalg.norm(dense - phi) / np.linalg.norm(phi))
    return dense, axF.eval_to_freq(coeffs), resid


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
    bandwidth wmax can represent.

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


def _apply(op, x):
    return (op @ x.ravel()).reshape(x.shape)


def _pair_amplitude(fx, phi_dense):
    """X_{l3l4}(k', i w_n') = sum_{ef} G2 phi -- the object whose window
    truncation ``tail_estimate`` models (norb = 1 here)."""
    from hwave.solver import eliashberg_dynamic as ed
    G2 = ed.calc_g2_dynamic(fx["green_u"], fx["beta"])
    return G2[0, 0, 0, 0] * phi_dense[0, 0]


_MEASURED = {}


def _measure(nmat, pairing):
    """All spec S3.5 quantities at one resolution (cached).

    Part-aware densification (step 3, amended): the fluctuation part is
    densified through the basis, the instantaneous part is carried exactly
    (it is flat in the external frequency by construction, S3.3), because a
    frequency-independent constant is precisely what the IR basis cannot
    represent -- see ``test_ir_naive_densification_is_not_a_measurement``.
    """
    key = (nmat, pairing)
    if key in _MEASURED:
        return _MEASURED[key]
    fx = _fixture(nmat, pairing)
    axF = fx["axF"]
    phi_d, phi_n, fit_resid = _gap_projected(axF, nmat, pairing)

    Y_uni = _apply(fx["u_full"], phi_d)
    Y_uni_pp = _apply(fx["u_pp"], phi_d)
    Y_ir_fl = axF.eval_to_uniform(
        axF.fit_from_freq(_apply(fx["i_fl"], phi_n)), nmat)
    Y_ir_pp_nodes = _apply(fx["i_pp"], phi_n)
    Y_ir_pp = np.broadcast_to(Y_ir_pp_nodes[..., :1], Y_uni.shape)
    Y_ir = Y_ir_fl + Y_ir_pp

    est, tail_rel, bad = bc.tail_estimate(_pair_amplitude(fx, phi_d),
                                          fx["beta"], nmat)
    out = dict(
        nmat=nmat, eps_ir=fx["eps_ir"], fit_resid=fit_resid,
        tail_rel=tail_rel, tail_bad=bad,
        Y_uni=Y_uni, Y_uni_pp=Y_uni_pp, Y_ir=Y_ir, Y_ir_pp=Y_ir_pp,
        E=float(np.linalg.norm(Y_uni - Y_ir) / np.linalg.norm(Y_uni)),
        E_pp=float(np.linalg.norm(Y_uni_pp - Y_ir_pp)
                   / np.linalg.norm(Y_uni_pp)),
        pp_weight=float(np.linalg.norm(Y_uni_pp) / np.linalg.norm(Y_uni)),
        lam_uni=_sector_lambda(fx, pairing, "uniform"),
        lam_ir=_sector_lambda(fx, pairing, "ir"))
    out["e_lam"] = float(abs(out["lam_uni"] - out["lam_ir"])
                         / abs(out["lam_uni"]))
    _MEASURED[key] = out
    return out


def _shared_window(Y, nmat, n_small):
    """``Y|_S(n_small)``: the stored window of the coarser grid, a plain slice
    of the finer one (S(N) is a subset of S(2N) under the centered map)."""
    n_t = np.arange(nmat) - nmat // 2
    mask = (n_t >= -(n_small // 2)) & (n_t <= n_small // 2 - 1)
    return Y[..., mask]


def _d_uni(nmat, pairing):
    """``D_uni(Nmat)`` and ``d_uni(Nmat)`` (spec S3.5 step 5, amended): the
    UNIFORM kernel's own self-convergence increment between Nmat and 2*Nmat on
    the shared window -- the correct model of the discretization the IR path
    is compared against."""
    lo = _measure(nmat, pairing)
    hi = _measure(2 * nmat, pairing)
    ref = _shared_window(hi["Y_uni"], 2 * nmat, nmat)
    D = float(np.linalg.norm(lo["Y_uni"] - ref) / np.linalg.norm(ref))
    d = float(abs(lo["lam_uni"] - hi["lam_uni"]) / abs(hi["lam_uni"]))
    return D, d


def _eps_ir_max(pairing):
    """``eps_IR^max``: one number for the whole study, the maximum over the
    tested resolution set {64, 128, 256}."""
    return max(_measure(n, pairing)["eps_ir"] for n in NMATS)


def _sector(A, shape, pairing):
    from scipy.sparse.linalg import LinearOperator
    from hwave.solver import eliashberg_dynamic as ed
    n = int(np.prod(shape))

    def mv(v):
        x = ed._project_parity_dynamic(np.asarray(v).reshape(shape), pairing)
        return ed._project_parity_dynamic(_apply(A, x), pairing).ravel()
    return LinearOperator((n, n), matvec=mv, dtype=complex)


def _sector_lambda(fx, pairing, which):
    """Leading eigenvalue INSIDE the parity sector (spec S3.5: branches are
    matched by sector, so no eigenpair-sorting rule is needed). The IR branch
    runs the complex-ARPACK semantics the amended section 2 makes explicit;
    on this fixture the mode is real to ~1e-12, and the metric uses the
    modulus of the difference either way."""
    from scipy.sparse.linalg import eigs
    from hwave.solver import eliashberg_dynamic as ed
    axF, nmat = fx["axF"], fx["nmat"]
    if which == "uniform":
        shape = (1, 1, NX, NY, NZ, nmat)
        v0 = _gap_projected(axF, nmat, pairing)[0]
        op = fx["u_full"]
    else:
        shape = (1, 1, NX, NY, NZ, axF.n_freq)
        v0 = _gap_projected(axF, nmat, pairing)[1]
        op = fx["i_full"]
    v0 = ed._project_parity_dynamic(v0, pairing).ravel()
    return complex(eigs(_sector(op, shape, pairing), k=1, which="LM", v0=v0,
                        tol=1e-10, return_eigenvectors=False)[0])


def _b1_kernels():
    """The B = 1 reduction data: the bond IR kernel with a single (on-site)
    channel, and the oracle-verified scalar ``eliashberg_kernel_ir`` fed the
    same vertex/green/gap. Shared by the reduction test and by the
    Hermiticity diagnostic (amended section 2 makes this reduction the
    PHYSICAL Hermiticity evidence for the IR branch)."""
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

    cphi = (rng.normal(size=(1, 1, NX, NY, NZ, axF.L))
            + 1j * rng.normal(size=(1, 1, NX, NY, NZ, axF.L)))
    phi_n = axF.eval_to_freq(cphi)
    Y_bond = _apply(A, phi_n)
    Y_scalar = ed.eliashberg_kernel_ir(
        ed._ir_vertex_to_rtau(Vs_n[None, None, None, None], axB, axF),
        ed.calc_g2_dynamic(green_n, BETA), phi_n.copy(), axF, BETA)
    return Y_bond, Y_scalar, vec_size, axF


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
    the implementation independently of every discretization budget below,
    and (amended section 2) it is the PHYSICAL Hermiticity evidence for the
    IR branch."""
    Y_bond, Y_scalar, vec_size, axF = _b1_kernels()
    assert vec_size == NX * NY * NZ * axF.n_freq
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


# --- spec S3.5 step 1 (amended): the normative gap generator ----------------
@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_normative_gap_generator_is_analytic_and_ir_representable(pairing,
                                                                  capsys):
    """Step 1's two requirements: (a) every resolution samples the SAME
    underlying analytic function -- so the shared Matsubara integers carry
    bit-identical values across grids, before and after the parity projection
    -- and (b) its IR fit residual sits at the basis tolerance, so the probe
    carries no out-of-subspace ambiguity into the comparison."""
    from hwave.solver import eliashberg_dynamic as ed
    axF = _axes()[0]
    resids = {}
    for nmat in NMATS:
        _, _, resid = _gap_projected(axF, nmat, pairing)
        resids[nmat] = resid
        assert resid <= 1.0e3 * IR_TOL, (nmat, resid)
    # (a) same analytic function at every resolution: the coarse grid's
    # frequencies are a subset of the fine grid's, and the parity partner of
    # every stored frequency is stored too (n~ <-> -n~-1), so the projection
    # does not break that.
    for lo, hi in ((64, 128), (128, 256)):
        np.testing.assert_array_equal(
            _gap_analytic(lo, pairing),
            _shared_window(_gap_analytic(hi, pairing), hi, lo))
        np.testing.assert_array_equal(
            ed._project_parity_dynamic(_gap_analytic(lo, pairing), pairing),
            _shared_window(ed._project_parity_dynamic(
                _gap_analytic(hi, pairing), pairing), hi, lo))
    with capsys.disabled():
        print("\n  [{}] gap IR fit residual = {}".format(
            pairing, {k: "%.2e" % v for k, v in resids.items()}))


@ir_only
def test_white_ir_coefficients_are_not_a_physical_probe(capsys):
    """Evidence for the amended step 1. White (unit-variance) IR coefficients
    saturate the basis bandwidth; a pointwise Matsubara product is a tau
    convolution, so the spectral supports ADD and ``X = G2 * phi`` falls
    outside the basis. Its IR equal-time value is then ~10% off the uniform
    window sum AND does not improve with resolution -- the comparison would
    measure the basis bandwidth, not the kernel. The normative analytic gap
    keeps the same quantity at the 1e-5 level."""
    from hwave.solver import eliashberg_dynamic as ed
    axF = _axes()[0]
    rng = np.random.default_rng(5)
    white = (rng.normal(size=(1, 1, NX, NY, NZ, axF.L))
             + 1j * rng.normal(size=(1, 1, NX, NY, NZ, axF.L)))

    def equal_time_mismatch(nmat, phi_d, phi_n):
        fx = _fixture(nmat, "singlet")
        G2u = ed.calc_g2_dynamic(fx["green_u"], BETA)[0, 0, 0, 0]
        window_sum = (G2u * phi_d[0, 0]).sum(axis=-1)
        G2n = ed.calc_g2_dynamic(fx["green_n"], BETA)[0, 0, 0, 0]
        X_l = axF.fit_from_freq(G2n * phi_n[0, 0])
        ir_equal_time = BETA * (X_l @ axF.u_zero_plus)
        return float(np.abs(window_sum - ir_equal_time).max()
                     / np.abs(window_sum).max())

    white_mis = {n: equal_time_mismatch(
        n, axF.eval_to_uniform(white, n), axF.eval_to_freq(white))
        for n in NMATS}
    phys_mis = {n: equal_time_mismatch(
        n, *_gap_projected(axF, n, "singlet")[:2])
                for n in NMATS}
    with capsys.disabled():
        print("\n  equal-time mismatch: white {} vs normative {}".format(
            {k: "%.2e" % v for k, v in white_mis.items()},
            {k: "%.2e" % v for k, v in phys_mis.items()}))
    # white: percent-level and NOT converging (ratio ~ 1 across a 4x change)
    assert min(white_mis.values()) > 1e-2
    assert white_mis[256] > 0.5 * white_mis[64]
    # normative: orders of magnitude smaller at every resolution
    assert max(phys_mis.values()) < 1e-4
    assert max(phys_mis.values()) < 0.01 * min(white_mis.values())


# --- spec S3.5 step 5 (amended): the three action clauses -------------------
@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_instantaneous_absolute_budget(pairing, capsys):
    """Amended step 5, INSTANTANEOUS part (absolute):
    ``E_pp(N) <= 3*(tail_est_rel(N) + eps_IR^max)`` -- the clause the tail
    estimator actually models, with the spec's vacuity guard
    ``||Y^pp_uni||_F >= 1e-3 * ||Y_uni||_F`` asserted so a gap that nulls the
    Cooper term cannot make it pass trivially."""
    eps_max = _eps_ir_max(pairing)
    rows = {}
    for nmat in NMATS:
        m = _measure(nmat, pairing)
        assert not m["tail_bad"], (
            "tail estimate unreliable at Nmat={} -- spec S3.3 forbids "
            "recording that as an inconclusive pass".format(nmat))
        assert m["pp_weight"] >= 1e-3, (
            "vacuity guard: ||Y^pp_uni||_F is only {:.3e} of ||Y_uni||_F; "
            "this fixture nulls the Cooper term".format(m["pp_weight"]))
        budget = 3.0 * (m["tail_rel"] + eps_max)
        rows[nmat] = (m["E_pp"], budget, m["pp_weight"])
        assert m["E_pp"] <= budget, (
            "E_pp({}) = {:.4e} exceeds 3*(tail_est_rel {:.4e} + eps_IR^max "
            "{:.3e}) = {:.4e}".format(nmat, m["E_pp"], m["tail_rel"],
                                      eps_max, budget))
    with capsys.disabled():
        print("\n  [{}] E_pp / budget (pp weight): {}".format(
            pairing, {k: "%.3e / %.3e (%.2f)" % v for k, v in rows.items()}))


@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_action_refinement(pairing, capsys):
    """Amended step 5, REFINEMENT: ``E(2N) <= 0.6*E(N) + eps_IR^max`` for two
    successive doublings. This is what establishes that the two kernels are
    discretizations of the same operator."""
    eps_max = _eps_ir_max(pairing)
    E = {n: _measure(n, pairing)["E"] for n in NMATS}
    with capsys.disabled():
        print("\n  [{}] E(Nmat) = {}".format(
            pairing, {k: "%.4e" % v for k, v in E.items()}))
    for lo, hi in ((64, 128), (128, 256)):
        assert E[hi] <= 0.6 * E[lo] + eps_max, (
            "no refinement {} -> {}: E = {:.4e} -> {:.4e} (budget {:.4e})"
            .format(lo, hi, E[lo], E[hi], 0.6 * E[lo] + eps_max))


@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_action_uniform_anchor(pairing, capsys):
    """Amended step 5, UNIFORM ANCHOR:
    ``E(N) <= 3*(D_uni(N) + tail_est_rel(N) + eps_IR^max)`` with
    ``D_uni(N) = ||Y_uni(N) - Y_uni(2N)|_S(N)||_F / ||Y_uni(2N)|_S(N)||_F``.

    The anchor replaces the original (mis-modeled) absolute budget: the
    dominant IR-vs-uniform difference is the uniform path's OWN O(1/Nmat)
    fluctuation truncation, and D_uni measures exactly that."""
    eps_max = _eps_ir_max(pairing)
    rows = {}
    for nmat in NMATS:
        m = _measure(nmat, pairing)
        D, _ = _d_uni(nmat, pairing)
        budget = 3.0 * (D + m["tail_rel"] + eps_max)
        rows[nmat] = (m["E"], D, budget, m["E"] / budget)
        assert m["E"] <= budget, (
            "E({}) = {:.4e} exceeds 3*(D_uni {:.4e} + tail_est_rel {:.4e} + "
            "eps_IR^max {:.3e}) = {:.4e}".format(
                nmat, m["E"], D, m["tail_rel"], eps_max, budget))
    with capsys.disabled():
        print("\n  [{}] E / D_uni / budget (ratio): {}".format(
            pairing,
            {k: "%.3e / %.3e / %.3e (%.2f)" % v for k, v in rows.items()}))


# --- spec S3.5 step 5 (amended): the eigenvalue twins -----------------------
@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_eigenvalue_refinement(pairing, capsys):
    """Eigenvalue twin of the refinement clause, per parity sector:
    ``e(2N) <= 0.6*e(N) + eps_IR^max``. No densification enters this metric,
    so it is also free of the frequency-flat-component question."""
    eps_max = _eps_ir_max(pairing)
    e_rel = {n: _measure(n, pairing)["e_lam"] for n in NMATS}
    with capsys.disabled():
        print("\n  [{}] e(Nmat) = {}  lam_uni = {}".format(
            pairing, {k: "%.4e" % v for k, v in e_rel.items()},
            {n: "%.6f" % _measure(n, pairing)["lam_uni"].real
             for n in NMATS}))
    for lo, hi in ((64, 128), (128, 256)):
        assert e_rel[hi] <= 0.6 * e_rel[lo] + eps_max, (
            "no eigenvalue refinement {} -> {}: e = {:.4e} -> {:.4e}"
            .format(lo, hi, e_rel[lo], e_rel[hi]))


@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_eigenvalue_uniform_anchor(pairing, capsys):
    """Eigenvalue twin of the uniform anchor:
    ``e(N) <= 3*(d_uni(N) + tail_est_rel(N) + eps_IR^max)`` with
    ``d_uni(N) = |lam_uni(N) - lam_uni(2N)| / |lam_uni(2N)|``."""
    eps_max = _eps_ir_max(pairing)
    rows = {}
    for nmat in NMATS:
        m = _measure(nmat, pairing)
        _, d = _d_uni(nmat, pairing)
        budget = 3.0 * (d + m["tail_rel"] + eps_max)
        rows[nmat] = (m["e_lam"], d, budget, m["e_lam"] / budget)
        assert m["e_lam"] <= budget, (
            "e({}) = {:.4e} exceeds 3*(d_uni {:.4e} + tail_est_rel {:.4e} + "
            "eps_IR^max {:.3e}) = {:.4e}".format(
                nmat, m["e_lam"], d, m["tail_rel"], eps_max, budget))
    with capsys.disabled():
        print("\n  [{}] e / d_uni / budget (ratio): {}".format(
            pairing,
            {k: "%.3e / %.3e / %.3e (%.2f)" % v for k, v in rows.items()}))


# --- evidence for the amended step 3 ----------------------------------------
@ir_only
def test_ir_naive_densification_is_not_a_measurement(capsys):
    """Evidence for the part-aware densification of step 3 (amended):
    fitting the FULL IR output (fluctuation + frequency-flat instantaneous
    term) to the IR basis and evaluating on the uniform grid does not converge
    at all -- a constant is not in the span of the IR basis, and the uniform
    grid keeps extending to frequencies where the fit has decayed to zero
    while the true answer is still the constant. The naive number therefore
    stays O(1) and REFUSES the refinement inequality that the part-aware
    metric passes, so it measures the densification, not the kernel. (With a
    white-coefficient probe it grew monotonically, 0.10 -> 0.27 -> 0.57, which
    is the number quoted in the amended spec; with the normative analytic gap
    it is non-monotonic but equally stuck, 0.56 -> 0.70 -> 0.57.)"""
    axF = _axes()[0]
    naive, part_aware = [], []
    for nmat in NMATS:
        fx = _fixture(nmat, "singlet")
        phi_d, phi_n, _ = _gap_projected(axF, nmat, "singlet")
        m = _measure(nmat, "singlet")
        Y_naive = axF.eval_to_uniform(
            axF.fit_from_freq(_apply(fx["i_full"], phi_n)), nmat)
        naive.append(float(np.linalg.norm(m["Y_uni"] - Y_naive)
                           / np.linalg.norm(m["Y_uni"])))
        part_aware.append(m["E"])
        # ... while the instantaneous part IS exactly flat on both axes, so
        # the broadcast densification is exact by construction.
        for key, probe in (("u_pp", phi_d), ("i_pp", phi_n)):
            Y = _apply(fx[key], probe)
            assert (np.abs(Y - Y[..., :1]).max()
                    <= 1e-12 * np.abs(Y).max()), key
    with capsys.disabled():
        print("\n  naive densification E = {}\n  part-aware E        = {}"
              .format(["%.3e" % v for v in naive],
                      ["%.3e" % v for v in part_aware]))
    # stuck at O(1) and failing the refinement inequality at every doubling...
    assert min(naive) > 0.1
    for lo, hi in ((0, 1), (1, 2)):
        assert naive[hi] > 0.6 * naive[lo]
    # ... while the part-aware metric converges monotonically on the same data
    assert part_aware[0] > part_aware[1] > part_aware[2]
    assert part_aware[-1] < 0.01 * naive[-1]


# --- amended section 2: the sampling residual is a DIAGNOSTIC ---------------
@ir_only
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_ir_sampling_hermiticity_is_a_diagnostic_not_a_gate(pairing, capsys):
    """Amended section 2 (2026-07-28): there is NO gated threshold on the
    Euclidean sampling-coordinate residual any more -- it is structural, not
    physical (the equal-time functional is adjoint to a constant ket in the
    physical L^2 metric, not in the Euclidean sampling metric), and the same
    property holds for the pre-existing scalar dynamic IR path, which has
    always run the complex-ARPACK ordering.

    (a) The PHYSICAL Hermiticity evidence for the IR branch is the B = 1
    reduction against the oracle-verified scalar ``eliashberg_kernel_ir``,
    asserted here as well as in its own test.
    (b) The sampling residual is reported only -- asserted finite, never
    gated -- together with its decomposition and the uniform-coordinate
    residual, which is the branch that IS gated (machine precision, P0-1).
    """
    Y_bond, Y_scalar, _, _ = _b1_kernels()
    rel = float(np.abs(Y_bond - Y_scalar).max() / np.abs(Y_scalar).max())
    assert rel <= 1e-11, "physical-metric evidence broken: {:.3e}".format(rel)

    nmat = 64
    fx = _fixture(nmat, pairing)
    axF = fx["axF"]
    shape_n = (1, 1, NX, NY, NZ, axF.n_freq)
    w_n = _pair_weight(fx["green_n"])
    res = _whitened_residual(_dense(fx["i_full"], shape_n), w_n)
    parts = {k: _whitened_residual(_dense(fx[k], shape_n), w_n)
             for k in ("i_fl", "i_pp")}
    w_u = _pair_weight(fx["green_u"])
    uni = _whitened_residual(_dense(fx["u_full"], (1, 1, NX, NY, NZ, nmat)),
                             w_u)
    with capsys.disabled():
        print("\n  [{}] IR sampling-coordinate residual (DIAGNOSTIC, not "
              "gated) = {:.4e} (fluctuation {:.4e}, instantaneous {:.4e}); "
              "uniform coordinates = {:.4e}; B=1 physical evidence = {:.2e}"
              .format(pairing, res, parts["i_fl"], parts["i_pp"], uni, rel))
    assert np.isfinite(res) and np.isfinite(uni)
    assert all(np.isfinite(v) for v in parts.values())
