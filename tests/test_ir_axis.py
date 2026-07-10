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
