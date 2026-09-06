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
            x_next = x + self.mix * r - (dX + self.mix * dR).T @ gamma
            if np.all(np.isfinite(x_next)):
                self.last_gamma = gamma
                return self._unstack(x_next, state)
        # singular / non-finite history: plain linear step and restart
        self._reset()
        self._xs.append(x)
        self._rs.append(r)
        return self._unstack(x + self.mix * r, state)


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
