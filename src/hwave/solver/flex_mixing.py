"""The provenance-split self-energy state of the Phase B FLEX loop and its
mixing (GitHub issue #181; spec 2026-09-06 sections 2.4 and 5.6).

``sigma = sigma_static + sigma_fluct``: the Hartree-Fock map and the
static mean-field seed are ``static`` (frequency independent, a singleton
frequency axis), the W map is ``fluct``. Both components are mixed with
the SAME affine coefficients: linear mixing component-wise, and a stacked
Anderson mixer whose least-squares system uses REAL coefficients (the
Gram matrix and right-hand side are the real parts of the complex inner
products), so a real affine combination of Hermitian static matrices
stays Hermitian. The legacy ``flex._AndersonMixer`` (total sigma, HF off)
is untouched.
"""
from dataclasses import dataclass

import numpy as np

_GEMV_CHUNK = 1 << 20


def _affine_combination(dX, dR, mix, gamma):
    """``(dX + mix * dR).T @ gamma`` in column chunks: numpy's bundled OpenBLAS
    segfaults in the complex-by-real gemv once a (m, N) operand has more than
    a few million columns per thread (reproduced without H-wave at
    OPENBLAS_NUM_THREADS=1, N=4e6, m=2); the chunks also avoid the full-size
    ``dX + mix * dR`` temporary."""
    n = dX.shape[1]
    out = np.empty(n, dtype=np.complex128)
    for s in range(0, n, _GEMV_CHUNK):
        sl = slice(s, s + _GEMV_CHUNK)
        out[sl] = (dX[:, sl] + mix * dR[:, sl]).T @ gamma
    return out


@dataclass(frozen=True)
class SplitState:
    """``static`` `(nb, 1, nvol, n, n)`, ``fluct`` `(nb, nmat, nvol, n, n)`."""
    static: np.ndarray
    fluct: np.ndarray

    def total(self):
        return self.static + self.fluct           # broadcast over frequency

    @property
    def nmat(self):
        return int(self.fluct.shape[1])


def linear_mix_pair(state, new, mix):
    """``(1 - mix) state + mix new`` on both components."""
    return SplitState(static=(1.0 - mix) * state.static + mix * new.static,
                      fluct=(1.0 - mix) * state.fluct + mix * new.fluct)


class StackedAndersonMixer:
    """Anderson acceleration on the stacked vector
    ``[broadcast(static); fluct]`` with real coefficients (spec 2.4).

    Same update, Tikhonov guard, fallback-to-linear and history restart
    as ``flex._AndersonMixer``; the ``dR``/``dX`` stacks are two
    preallocated ``(depth - 1, N)`` arrays filled in place (spec 3.6).
    """

    def __init__(self, mix, depth):
        self.mix = float(mix)
        self.depth = max(1, int(depth))
        self._xs = []
        self._rs = []
        self._dR = None
        self._dX = None
        self.last_gamma = None

    def stack(self, state):
        nmat = state.nmat
        return np.concatenate([np.repeat(state.static, nmat, axis=1).reshape(-1),
                               state.fluct.reshape(-1)])

    def _unstack(self, vec, like):
        n_static_b = like.static.size * like.nmat
        st = vec[:n_static_b].reshape(like.fluct.shape)[:, :1]
        fl = vec[n_static_b:].reshape(like.fluct.shape)
        return SplitState(static=np.ascontiguousarray(st), fluct=np.ascontiguousarray(fl))

    def _reset(self):
        self._xs = []
        self._rs = []

    def step(self, state, new):
        x = self.stack(state)
        r = self.stack(new) - x
        self._xs.append(x)
        self._rs.append(r)
        if len(self._xs) > self.depth:
            self._xs.pop(0)
            self._rs.pop(0)
        m = len(self._xs)
        self.last_gamma = None
        if m == 1:
            return self._unstack(x + self.mix * r, state)
        N = x.size
        if self._dR is None or self._dR.shape != (self.depth - 1, N):
            self._dR = np.empty((self.depth - 1, N), dtype=np.complex128)
            self._dX = np.empty((self.depth - 1, N), dtype=np.complex128)
        dR = self._dR[:m - 1]
        dX = self._dX[:m - 1]
        for i in range(m - 1):
            np.subtract(self._rs[i + 1], self._rs[i], out=dR[i])
            np.subtract(self._xs[i + 1], self._xs[i], out=dX[i])
        G = (dR.conj() @ dR.T).real
        b = (dR.conj() @ r).real
        lam = 1.0e-10 * max(float(np.trace(G)) / (m - 1), 1.0e-300)
        try:
            gamma = np.linalg.solve(G + lam * np.eye(m - 1), b)
        except np.linalg.LinAlgError:
            gamma = None
        if gamma is not None:
            x_next = x + self.mix * r - _affine_combination(dX, dR, self.mix, gamma)
            if np.all(np.isfinite(x_next)):
                self.last_gamma = gamma
                return self._unstack(x_next, state)
        # singular / non-finite history: plain linear step and restart
        self._reset()
        self._xs.append(x)
        self._rs.append(r)
        return self._unstack(x + self.mix * r, state)


# ---------------------------------------------------------------------------
# Frequency-reversal Hermitian symmetry of the self-energy (issue #281)
# ---------------------------------------------------------------------------
#
# For a Hamiltonian H that is Hermitian -- a Hermitian one-body part H0(k)
# (plus a Hermitian static mean field) and Hermitian two-body terms -- with
# a real chemical potential, the imaginary-time Green function obeys
#
#     G_ab(k, tau)^* = G_ba(k, tau),   0 < tau < beta,
#
# (take the complex conjugate of -Tr[e^{-(beta-tau)K} c_ka e^{-tau K} c_kb^dag]/Z
# with K = H - mu N Hermitian), i.e. G(k, tau) is a Hermitian matrix at
# every k. On the Matsubara axis this reads
#
#     G_ab(k, i w_n)^* = G_ba(k, -i w_n),   i.e.  G(k, i w) = G(k, -i w)^dag,
#
# and the same holds for the self-energy, since
# Sigma(k, i w) = (i w + mu) - H0(k) - G(k, i w)^{-1}. The relation is at
# FIXED k: no time reversal or inversion is involved.
#
# The FLEX map preserves it exactly in exact arithmetic. In (k, tau) space
# every step is a product of matrices that are Hermitian in their index
# pairs: the bubble chi0_{ab,cd}(q, tau) ~ G_ca(k+q, tau) G_bd(k, -tau)
# satisfies chi0_{ab,cd}^* = chi0_{ba,dc}; the RPA/FLEX ladders sum it with
# the bare vertex, which has the same pair symmetry exactly when the
# interaction is Hermitian (the k-space reader rejects a declaration table
# that is not Hermitian-closed, issue #93); and the self-energy contraction
# Sigma_ab(k, tau) ~ sum_q V_{..}(q, tau) G_{..}(k - q, tau) maps such
# products back to a Hermitian G-like object. The Hartree-Fock term
# (Phase B) is a frequency-independent Hermitian matrix. This covers the
# reduced, general and spin-orbital vertices, the bond-gate blocks, and a
# trans_mod / green_init mean field (which only changes the Hermitian H0
# that the solver diagonalises with ``eigh`` -- the Green function is
# rebuilt from those eigenpairs, so it is Hermitian by construction).
#
# The discretised maps keep the relation as well: the uniform grid
# i w_j = (2 j + 1 - Nmat) pi / beta is symmetric (index j <-> Nmat - 1 - j),
# the bosonic grid likewise, and the cyclic FFT convolutions commute with
# the reversal. On the IR basis the fermionic node set is integer-symmetric
# (asserted in ir_axis.IRAxis), and because the basis functions u_l(tau) are
# real, uhat_l(-i w) = uhat_l(i w)^*: projecting the NODE values is the same
# as taking the Hermitian part of every expansion coefficient rho_l (the fit
# matrix P = pinv(E) satisfies P R = conj(P) for the reversal R), so the
# projection is exact there too and the densified output inherits it.
#
# What the map does NOT do is damp an asymmetric part that is already
# present: past a bare Stoner instability the early iterations amplify the
# round-off asymmetric component together with the transient, and linear
# mixing then removes it only at the convergence rate (issue #281). The
# projection below removes it after every update. It changes a symmetric
# iterate only at round-off, so the fixed point is unchanged.


def frequency_reversed_dagger(sigma):
    """``Sigma(k, -i w)^dag`` of a ``(nblock, nfreq, nvol, n, n)`` array
    whose frequency axis (axis 1) is symmetric under reversal (the centered
    uniform grid, the IR node set, or a singleton static axis). Works on
    numpy and cupy arrays."""
    return sigma[:, ::-1].swapaxes(-1, -2).conj()


def _deviation(sigma, rev):
    if not sigma.size:
        return 0.0, 0.0
    a = float(abs(sigma - rev).max())
    scale = float(abs(sigma).max())
    return a, (a / scale if scale > 0.0 else a)


def frequency_hermitian_deviation(sigma):
    """``(abs, rel)``: ``max |Sigma(k, i w) - Sigma(k, -i w)^dag|`` and its
    ratio to ``max |Sigma|`` (``abs`` itself when ``Sigma`` vanishes)."""
    return _deviation(sigma, frequency_reversed_dagger(sigma))


def project_frequency_hermitian(sigma, with_deviation=False):
    """``(Sigma(k, i w) + Sigma(k, -i w)^dag) / 2`` (a new array); with
    ``with_deviation`` also the ``(abs, rel)`` deviation of the INPUT."""
    rev = frequency_reversed_dagger(sigma)
    out = 0.5 * (sigma + rev)
    if with_deviation:
        return out, _deviation(sigma, rev)
    return out


def project_split_state(state):
    """Both components of a :class:`SplitState` projected; the singleton
    static axis makes the static projection ``(S + S^dag) / 2``."""
    return SplitState(static=project_frequency_hermitian(state.static),
                      fluct=project_frequency_hermitian(state.fluct))


def residuals(state, new, G_state, G_new, eps_den=1e-300):
    """``(res_sigma, res_G, res_component)`` of spec 5.6, all before mixing."""
    tot, tot_new = state.total(), new.total()
    # res_sigma: exactly FLEX._calc_convergence (diff / ||sigma_new||, the
    # bare difference when ||sigma_new|| < 1e-30)
    diff = float(np.linalg.norm((tot_new - tot).ravel()))
    norm_new = float(np.linalg.norm(tot_new.ravel()))
    res_sigma = diff if norm_new < 1.0e-30 else diff / norm_new
    res_G = float(np.linalg.norm((np.asarray(G_new) - np.asarray(G_state)).ravel())
                  / max(float(np.linalg.norm(np.asarray(G_state).ravel())), eps_den))
    nmat = state.nmat
    d_static = np.repeat(new.static - state.static, nmat, axis=1).ravel()
    d_fluct = (new.fluct - state.fluct).ravel()
    n_static = np.repeat(new.static, nmat, axis=1).ravel()
    num = np.sqrt(float(np.vdot(d_static, d_static).real + np.vdot(d_fluct, d_fluct).real))
    den = np.sqrt(float(np.vdot(n_static, n_static).real + np.vdot(new.fluct.ravel(), new.fluct.ravel()).real))
    res_component = num / max(den, eps_den)
    return res_sigma, res_G, res_component


class PassCounter:
    """Three consecutive iterations with every residual below ``eps``."""

    def __init__(self, eps, needed=3):
        self.eps = float(eps)
        self.needed = int(needed)
        self.count = 0

    def update(self, res_sigma, res_G, res_component):
        if res_sigma < self.eps and res_G < self.eps and res_component < self.eps:
            self.count += 1
        else:
            self.count = 0
        return (self.count >= self.needed), self.count
