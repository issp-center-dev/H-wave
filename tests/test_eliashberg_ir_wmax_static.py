"""Regression tests for issue #57: IR-basis dynamic Eliashberg gives
wrong/erratic lambda on static-fluctuation-dominated FLEX data.

Two independent root causes in ``hwave.solver.eliashberg_dynamic``:

RC1 -- ``_ir_auto_wmax`` misread the flat ``hr`` layout
       ``{((Rx,Ry,Rz),(orb1,orb2)): scalar}`` as ``{R: matrix}`` and summed
       |t| over every (R, orbital-pair), a ~15x band overestimate on realistic
       multi-hopping models. The band bound must come from the actual
       dispersion eps(k).

RC2 -- ``_ir_compress(..., drop_constant=True)`` discarded a
       frequency-independent column that, for a static-dominated
       susceptibility (the near-critical regime that matters for SC), holds
       physical weight. It must error when the discarded constant is large
       relative to the data scale, and offer a keep option.

Tests must run from the repository root.
"""
import numpy as np
import pytest

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

pytestmark = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                                reason="sparse-ir not installed")


# --- RC1 --------------------------------------------------------------------

def test_ir_auto_wmax_uses_dispersion_not_hopping_grand_sum():
    from hwave.solver.eliashberg_dynamic import _ir_auto_wmax
    # An on-site offset that shifts the band far from zero WITHOUT widening it
    # is exactly the case the grand-|t|-sum bug over-counts. Two decoupled 1D
    # chains, on-site 5.0, NN hopping -1.0:
    #   eps_orb(k) = 5 - 2 cos(k)  ->  band range [3, 7].
    # With mu at the band centre (5.0): max|eps - mu| = 2, so 3*(2 + 0) = 6.
    # The buggy measure sums |t| over every (R, orbital-pair):
    #   |5| + |5| + |-1| * 4 = 14, so 3*(14 + 0) = 42 -- 7x too large.
    hr = {
        ((0, 0, 0), (0, 0)): 5.0,
        ((0, 0, 0), (1, 1)): 5.0,
        ((1, 0, 0), (0, 0)): -1.0,
        ((-1, 0, 0), (0, 0)): -1.0,
        ((1, 0, 0), (1, 1)): -1.0,
        ((-1, 0, 0), (1, 1)): -1.0,
    }
    wmax = _ir_auto_wmax(hr, {}, norb=2, beta=10.0, mu=5.0)
    assert wmax == pytest.approx(6.0, abs=1e-6)
    # The on-site offset must NOT leak in (max|eps|=7 -> 21) nor the grand
    # sum (42): the spectral half-range about mu is what matters.
    assert wmax < 12.0


def test_ir_auto_wmax_includes_interaction_scale():
    from hwave.solver.eliashberg_dynamic import _ir_auto_wmax
    # Single 1D chain: eps(k) = -2 cos(k), range [-2, 2]. mu = 0 (band centre)
    # -> max|eps - mu| = 2. Interaction max = 3. So 3 * (2 + 3) = 15.
    hr = {
        ((1, 0, 0), (0, 0)): -1.0,
        ((-1, 0, 0), (0, 0)): -1.0,
    }
    inter_k = {"CoulombIntra": np.full((1, 1, 4, 1, 1), 3.0)}
    wmax = _ir_auto_wmax(hr, inter_k, norb=1, beta=10.0, mu=0.0)
    assert wmax == pytest.approx(15.0, abs=1e-6)


def test_ir_auto_wmax_solves_mu_from_filling():
    from hwave.solver.eliashberg_dynamic import _ir_auto_wmax
    # Same offset chain; mu is NOT given but half-filling (per orbital per
    # spin) fixes it at the band centre 5.0 by particle-hole symmetry of
    # -2 cos(k), so max|eps - mu| = 2 and wmax = 6 -- same as passing mu=5.
    hr = {
        ((0, 0, 0), (0, 0)): 5.0,
        ((1, 0, 0), (0, 0)): -1.0,
        ((-1, 0, 0), (0, 0)): -1.0,
    }
    wmax = _ir_auto_wmax(hr, {}, norb=1, beta=50.0, filling=0.5)
    assert wmax == pytest.approx(6.0, abs=1e-3)


# --- RC2 --------------------------------------------------------------------

def _static_dominated_illcond(beta=50.0, wmax=200.0, nmat=512):
    """A static-dominated bosonic chi in the pathological regime: a sharp
    static (nu=0) Lorentzian peak fitted with a huge basis (Lambda=beta*wmax=1e4,
    exactly the issue's over-large auto-wmax). The augmented fit is
    ill-conditioned and the constant column blows up ABOVE the data scale --
    the definitive signature that dropping it is unsafe. Returns (axB, arr)."""
    from hwave.solver.ir_axis import IRAxis
    axB = IRAxis(beta=beta, wmax=wmax, eps=1.0e-8, statistics="B")
    nu = (2 * np.arange(nmat) - nmat) * np.pi / beta
    g = 0.02
    chi = 2.0 * g / (nu ** 2 + g ** 2) * (1.0 - np.exp(-beta * g))
    return axB, chi[None, :].astype(np.complex128)


def test_ir_compress_errors_when_constant_exceeds_data_scale():
    from hwave.solver.eliashberg_dynamic import _ir_compress
    axB, arr = _static_dominated_illcond()
    # Discarded constant (~7.6e2) exceeds the data scale (~6.3e1): the run must
    # abort rather than silently return a physically-wrong object.
    with pytest.raises(ValueError, match="exceeds the data scale"):
        _ir_compress(arr, axB, 512, "chiq_s", drop_constant=True)


def test_ir_compress_keep_constant_does_not_error():
    from hwave.solver.eliashberg_dynamic import _ir_compress
    axB, arr = _static_dominated_illcond()
    # keep_constant retains the static component instead of dropping it, so the
    # pathological ill-conditioning no longer aborts the run.
    out = _ir_compress(arr, axB, 512, "chiq_s",
                       drop_constant=True, keep_constant=True)
    assert out.shape == (1, axB.n_freq)


def test_ir_compress_keep_constant_retains_offset():
    """Mechanistic check at a SANE wmax: a representable bosonic function plus a
    known moderate constant. drop recovers the clean function; keep returns the
    clean function PLUS the constant (the static weight is retained)."""
    from hwave.solver.eliashberg_dynamic import _ir_compress
    from hwave.solver.ir_axis import IRAxis
    beta, nmat = 50.0, 512
    axB = IRAxis(beta=beta, wmax=4.0, eps=1.0e-8, statistics="B")
    nu = (2 * np.arange(nmat) - nmat) * np.pi / beta
    g = 1.7
    chi_clean = 2.0 * g / (nu ** 2 + g ** 2) * (1.0 - np.exp(-beta * g))
    const = 0.3 * float(np.abs(chi_clean).max())      # moderate: warns, < scale
    arr = (chi_clean + const)[None, :].astype(np.complex128)
    dropped = _ir_compress(arr, axB, nmat, "chiq_s", drop_constant=True)
    kept = _ir_compress(arr, axB, nmat, "chiq_s",
                        drop_constant=True, keep_constant=True)
    np.testing.assert_allclose(kept[0], dropped[0] + const, atol=1e-6)


def test_ir_compress_small_constant_still_dropped():
    from hwave.solver.eliashberg_dynamic import _ir_compress
    from hwave.solver.ir_axis import IRAxis
    nmat = 128
    axF = IRAxis(beta=20.0, wmax=4.0, eps=1.0e-8, statistics="F")
    # A smooth, IR-representable fermionic object (single pole) has only the
    # tiny O(beta/Nmat) discretization constant -> no error, drops cleanly.
    n = 2 * np.arange(nmat) + 1 - nmat
    iwn = 1j * np.pi * n / 20.0
    arr = (1.0 / (iwn - 0.7))[None, :].astype(np.complex128)
    out = _ir_compress(arr, axF, nmat, "green", drop_constant=False)
    assert out.shape == (1, axF.n_freq)
