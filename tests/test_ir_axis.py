"""Unit tests for solver/ir_axis.py (sparse-ir wrapper).

Conventions under test (design doc docs/design/ir-matsubara.md, Sec. 3.1/3.5):
integer Matsubara node indices, node sets exactly symmetric so reflections are
plain axis reversals, pure basis-change transforms contracting the last axis,
eps-accurate uniform-grid round trips for physically smooth objects, and the
beta^- route to equal-time observables.
"""
import numpy as np
import pytest

# Import-safe under BOTH pytest and unittest discovery (the CI runs
# `unittest discover`, where pytest.importorskip at module level would turn
# a missing optional dependency into an import ERROR instead of a skip).
try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

pytestmark = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                                reason="sparse-ir not installed")


def _axis(statistics, beta=50.0, wmax=8.0, eps=1e-8):
    from hwave.solver.ir_axis import IRAxis
    return IRAxis(beta=beta, wmax=wmax, eps=eps, statistics=statistics)


def test_node_sets_are_symmetric_and_sorted():
    for stat, parity in (("F", 1), ("B", 0)):
        ax = _axis(stat)
        n = ax.freq_n
        assert np.all(n % 2 == parity), "F nodes odd / B nodes even integers"
        assert np.array_equal(n, np.sort(n))
        # exact integer symmetry: negation == reversal
        assert np.array_equal(-n[::-1], n)
        t = ax.tau
        assert np.array_equal(np.sort(t), t)
        # reflection is implemented as index reversal: mirror pairs sum to beta
        np.testing.assert_allclose(t + t[::-1], ax.beta,
                                   rtol=0.0, atol=1e-13 * ax.beta)


def test_uniform_roundtrip_greens_function():
    """A physically smooth G (pole structure) on the uniform centered H-wave
    grid must survive uniform -> IR -> uniform within ~eps."""
    beta, nmat = 50.0, 512
    ax = _axis("F", beta=beta)
    rng = np.random.default_rng(1)
    poles = rng.uniform(-6.0, 6.0, size=5)
    weights = rng.uniform(0.1, 1.0, size=5)
    weights /= weights.sum()
    # centered H-wave grid: iw_j = (2j+1-nmat) * pi / beta
    wn = 2 * np.arange(nmat) + 1 - nmat
    iw = 1j * wn * np.pi / beta
    G = (weights[None, :] / (iw[:, None] - poles[None, :])).sum(axis=1)

    coeffs = ax.fit_from_uniform(G, nmat)
    G_back = ax.eval_to_uniform(coeffs, nmat)
    assert np.abs(G_back - G).max() < 1e-7

    # and through the sparse nodes: uniform -> coeffs -> nodes -> coeffs
    g_nodes = ax.eval_to_freq(coeffs)
    coeffs2 = ax.fit_from_freq(g_nodes)
    assert np.abs(coeffs2 - coeffs).max() < 1e-9


def test_freq_tau_composite_transforms_consistent():
    """freq_to_tau followed by tau_to_freq must be the identity on
    IR-representable data (contracting the LAST axis, leading batch dims)."""
    ax = _axis("F")
    rng = np.random.default_rng(2)
    coeffs = rng.standard_normal((3, 4, ax.L)) + 1j * rng.standard_normal((3, 4, ax.L))
    g_nodes = ax.eval_to_freq(coeffs)
    g_tau = ax.freq_to_tau(g_nodes)
    g_nodes2 = ax.tau_to_freq(g_tau)
    np.testing.assert_allclose(g_nodes2, g_nodes, atol=1e-9)
    assert g_tau.shape == (3, 4, ax.n_tau)


def test_beta_minus_gives_fermi_occupation():
    """For G0(iw) = 1/(iw - e), the density from the beta^- evaluation must be
    the Fermi function: -G(beta^-) = f(e) with G(tau) = -<T c(tau) c+>."""
    beta = 50.0
    ax = _axis("F", beta=beta)
    for e in (-2.0, -0.3, 0.4, 3.0):
        g_nodes = 1.0 / (1j * ax.freq_n * np.pi / beta - e)
        coeffs = ax.fit_from_freq(g_nodes)
        n = -float(np.real(coeffs @ ax.u_beta_minus))
        f = 1.0 / (1.0 + np.exp(np.clip(beta * e, -700, 700)))
        assert abs(n - f) < 1e-7, "e={}: n={} f={}".format(e, n, f)


def test_bosonic_eval_on_foreign_tau_nodes():
    """Bosonic coefficients must evaluate exactly on the FERMIONIC tau nodes
    (the kernel's product grid): pin against a direct evaluation of a smooth
    bosonic function chi(iv) = 2g/(v^2+g^2) <-> exp decay in tau."""
    beta = 50.0
    axf = _axis("F", beta=beta)
    axb = _axis("B", beta=beta)
    g = 1.5
    # chi(tau) = exp(-g*tau)+exp(-g*(beta-tau)) normalized-ish (periodic, smooth)
    nu = axb.freq_n * np.pi / beta   # bosonic real frequencies nu_m = 2m pi/beta
    chi_nodes = (2.0 * g / (nu ** 2 + g ** 2)
                 * (1.0 - np.exp(-beta * g)) / 1.0)
    # analytic transform of e^{-g tau} + e^{-g (beta-tau)} on 0..beta:
    # integral_0^beta e^{i nu tau} (e^{-g tau} + e^{-g(beta-tau)}) dtau
    chi_exact_tau = (np.exp(-g * axf.tau) + np.exp(-g * (beta - axf.tau)))
    coeffs = axb.fit_from_freq(chi_nodes)
    chi_on_ftau = axb.eval_to_tau_points(coeffs, axf.tau)
    np.testing.assert_allclose(chi_on_ftau, chi_exact_tau, atol=1e-6)


def test_missing_sparse_ir_raises_actionable(monkeypatch):
    import hwave.solver.ir_axis as ir_axis
    monkeypatch.setattr(ir_axis, "_import_sparse_ir",
                        lambda: (_ for _ in ()).throw(ImportError("nope")))
    with pytest.raises(ImportError, match="sparse-ir"):
        ir_axis.IRAxis(beta=10.0, wmax=5.0, eps=1e-8, statistics="F")


def test_compress_drop_constant_isolates_known_offset(caplog):
    """A representable bosonic function plus a KNOWN frequency-independent
    constant: drop_constant must recover the clean function on the nodes
    (the constant fully isolated), and an unusually large constant must be
    surfaced with a warning (review follow-up: the residual alone cannot
    see it, since the augmented fit absorbs the constant exactly)."""
    import logging
    from hwave.solver import eliashberg_dynamic as ed
    beta, nmat = 50.0, 512
    ax = _axis("B", beta=beta)
    g = 1.7
    nu_u = (2 * np.arange(nmat) - nmat) * np.pi / beta
    chi_clean = 2.0 * g / (nu_u ** 2 + g ** 2) * (1.0 - np.exp(-beta * g))
    const = 0.5 * np.abs(chi_clean).max()          # large, must warn
    chi = (chi_clean + const)[None, :]

    with caplog.at_level(logging.WARNING, logger="qlms"):
        nodes = ed._ir_compress(chi.astype(complex), ax, nmat, "test",
                                drop_constant=True)
    nu_n = ax.freq_n * np.pi / beta
    chi_nodes_clean = 2.0 * g / (nu_n ** 2 + g ** 2) * (1.0 - np.exp(-beta * g))
    np.testing.assert_allclose(nodes[0], chi_nodes_clean, atol=1e-6)
    assert any("unusually large" in rec.message for rec in caplog.records)


# ---------------------------------------------------------------------------
# Stage 3: fit_from_freq_points (values at ARBITRARY integer Matsubara nodes)
# ---------------------------------------------------------------------------

def _gf_at(freq_n, beta, e=0.7):
    """Reference representable function: G(iw_n) = 1/(iw_n - e)."""
    iw = 1j * np.asarray(freq_n) * np.pi / beta
    return 1.0 / (iw - e)


def test_fit_from_freq_points_identity_on_own_nodes():
    """At the axis' own node set the fit must reproduce fit_from_freq."""
    ax = _axis("F")
    vals = _gf_at(ax.freq_n, ax.beta)[None, :]
    c_pts = ax.fit_from_freq_points(vals, ax.freq_n)
    c_own = ax.fit_from_freq(vals)
    np.testing.assert_allclose(c_pts, c_own, atol=1e-12)


def test_fit_from_freq_points_foreign_node_set():
    """Values on a DIFFERENT basis' nodes (a tighter-eps writer -> more
    nodes than this axis' L: overdetermined, the supported direction) must
    recover coefficients that evaluate correctly on the run axis' own
    nodes. The reverse (fewer file nodes than the run L) raises -- covered
    by the underdetermined case in the validation test."""
    ax_run = _axis("F", eps=1e-6)
    ax_file = _axis("F", eps=1e-8)          # tighter basis -> more nodes
    assert not np.array_equal(ax_file.freq_n, ax_run.freq_n)
    assert ax_file.n_freq > ax_run.L
    vals_file = _gf_at(ax_file.freq_n, ax_run.beta)[None, :]
    coeffs = ax_run.fit_from_freq_points(vals_file, ax_file.freq_n)
    vals_run = ax_run.eval_to_freq(coeffs)
    ref_run = _gf_at(ax_run.freq_n, ax_run.beta)[None, :]
    np.testing.assert_allclose(vals_run, ref_run, atol=1e-5)


def test_fit_from_freq_points_validation_errors():
    ax = _axis("F")
    n_ok = ax.freq_n
    vals = _gf_at(n_ok, ax.beta)[None, :]
    # wrong parity (even indices into a fermionic axis)
    with pytest.raises(ValueError, match="parity"):
        ax.fit_from_freq_points(vals, n_ok + 1)
    # not strictly increasing (duplicate)
    bad = n_ok.copy()
    bad[1] = bad[0]
    with pytest.raises(ValueError, match="strictly increasing"):
        ax.fit_from_freq_points(vals, bad)
    # not 1-D
    with pytest.raises(ValueError, match="1-D"):
        ax.fit_from_freq_points(vals, n_ok.reshape(-1, 1))
    # underdetermined (fewer nodes than L)
    few = n_ok[:ax.L - 1]
    with pytest.raises(ValueError, match="at least"):
        ax.fit_from_freq_points(vals[..., :ax.L - 1], few)


def test_fit_from_freq_points_matrix_cached():
    """Second call with the same node set must reuse the cached matrix."""
    ax = _axis("B")
    n = ax.freq_n
    vals = (1.0 / (1.0 + (np.asarray(n) * np.pi / ax.beta) ** 2))[None, :]
    ax.fit_from_freq_points(vals.astype(complex), n)
    n_keys = len(ax._device_m)
    ax.fit_from_freq_points(vals.astype(complex), n)
    assert len(ax._device_m) == n_keys
