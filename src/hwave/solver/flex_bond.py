"""The dynamic bond pipeline of the Phase B FLEX solver (GitHub issue
#181; spec docs/superpowers/specs/2026-09-06-flex-bond-sigma-phase-b-181-design.md,
section 3): the dense bond-block store, the bubble assembly, the batched
dressing, the bond-aware self-energy transport, the memory preflight and
the output helpers.
"""
import numpy as np

from . import backend as _bk
from . import matsubara as _ms
from . import bubble as _bubble


class BondBlockStore:
    """Owner of the dense ``(nmat, nvol, ND, ND)`` bond arrays (spec 3.5).

    A context manager: construction allocates every named array (rolling
    back on a partial failure), ``release()`` runs on every exit of the
    ``with`` block. Pair access (``put_pair`` / ``get_pair``: a view
    ``(nmat, nvol, nd, nd)`` onto channel pair ``(alpha, beta)``) and
    frequency-batch access (``put_freq_batch`` / ``get_freq_batch``: a
    contiguous ``(nb, nvol, ND, ND)`` view) are the only ways in and out;
    ``detach(name)`` transfers ownership of one array to the caller.
    A future streaming store implements the same five methods.
    """

    def __init__(self, nmat, nvol, ND, nd, names):
        names = tuple(names)
        if len(set(names)) != len(names):
            raise ValueError("BondBlockStore: duplicate array names {}".format(names))
        if ND % nd != 0:
            raise ValueError("BondBlockStore: ND={} is not a multiple of nd={}".format(ND, nd))
        self.nmat, self.nvol, self.ND, self.nd = int(nmat), int(nvol), int(ND), int(nd)
        self.B = self.ND // self.nd
        self._arrays = {}
        try:
            for name in names:
                self._arrays[name] = np.zeros((self.nmat, self.nvol, self.ND, self.ND),
                                              dtype=np.complex128)
        except MemoryError:
            self._arrays.clear()
            raise
        self._released = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def _arr(self, name):
        if self._released:
            raise RuntimeError("BondBlockStore: released")
        return self._arrays[name]

    def put_pair(self, name, alpha, beta, block):
        a = self._arr(name)
        nd = self.nd
        a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd] = block

    def get_pair(self, name, alpha, beta):
        a = self._arr(name)
        nd = self.nd
        return a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd]

    def put_freq_batch(self, name, l0, l1, batch):
        self._arr(name)[l0:l1] = batch

    def get_freq_batch(self, name, l0, l1):
        return self._arr(name)[l0:l1]

    def detach(self, name):
        a = self._arr(name)
        del self._arrays[name]
        return a

    def release(self):
        self._arrays.clear()
        self._released = True


def assemble_bubble(store, green_scf, green0_tail, beta, view, spatial_shape, workers):
    """Fill ``store['chibar']`` pair by pair from ``bubble._iter_bond_dynamic``
    (spec 3.1). ``green_scf`` is the TAIL-SUBTRACTED single-block Green
    function the general bubble consumes; ``green0_tail`` its tail."""
    for (alpha, beta_), block in _bubble._iter_bond_dynamic(
            green_scf, green0_tail, beta, view, spatial_shape=tuple(spatial_shape),
            workers=workers):
        store.put_pair("chibar", alpha, beta_, _bk.to_host(block))
        del block


# =============================================================================
# Batched dressing, effective interaction and collapses (spec 3.2-3.3)
# =============================================================================

from dataclasses import dataclass as _dataclass

from . import bond_channels as _bc


@_dataclass(frozen=True)
class DressResult:
    collapse0: np.ndarray      # (nmat, nvol, nd, nd)  channel-0 block of chibar
    collapse_s: np.ndarray     # (nmat, nvol, nd, nd)  channel-0 block of chi_s
    collapse_c: np.ndarray
    static_s: np.ndarray       # (nvol, ND, ND) at Omega = 0
    static_c: np.ndarray
    cond_min_s: float
    cond_min_c: float


def dress_and_build_w(store, S, C, *, nb, output_full, nmat, nvol, nd, spatial_shape,
                      cond_tol=_bc._BOND_COND_FLOOR):
    """The spec 3.2 loop: per frequency batch, dress spin then charge (one
    channel batch alive at a time), consume each into W, the channel-0
    collapses, the static slices (slice assignment into preallocated
    buffers) and, with ``output_full``, the store's ``chi_s_w``/``chi_c_w``."""
    ND = S.shape[-1]
    SpC = S + C
    collapse0 = np.empty((nmat, nvol, nd, nd), dtype=np.complex128)
    collapse_s = np.empty_like(collapse0)
    collapse_c = np.empty_like(collapse0)
    static_s = np.zeros((nvol, ND, ND), dtype=np.complex128)
    static_c = np.zeros((nvol, ND, ND), dtype=np.complex128)
    l_static = nmat // 2
    cond_s = cond_c = np.inf
    for l0 in range(0, nmat, nb):
        l1 = min(nmat, l0 + nb)
        cb = store.get_freq_batch("chibar", l0, l1)
        collapse0[l0:l1] = cb[:, :, :nd, :nd]
        chi_s_b, cs = _bc.dress_batch(cb, S, "spin", l0=l0, nmat=nmat, spatial_shape=spatial_shape,
                                      cond_tol=cond_tol)
        cond_s = min(cond_s, cs if cs is not None else np.inf)
        W_b = 1.5 * (S[None] @ chi_s_b @ S[None])
        collapse_s[l0:l1] = chi_s_b[:, :, :nd, :nd]
        if l0 <= l_static < l1:
            static_s[...] = chi_s_b[l_static - l0]
        if output_full:
            store.put_freq_batch("chi_s_w", l0, l1, chi_s_b)
        del chi_s_b
        chi_c_b, cc = _bc.dress_batch(cb, C, "charge", l0=l0, nmat=nmat, spatial_shape=spatial_shape,
                                      cond_tol=cond_tol)
        cond_c = min(cond_c, cc if cc is not None else np.inf)
        W_b += 0.5 * (C[None] @ chi_c_b @ C[None])
        collapse_c[l0:l1] = chi_c_b[:, :, :nd, :nd]
        if l0 <= l_static < l1:
            static_c[...] = chi_c_b[l_static - l0]
        if output_full:
            store.put_freq_batch("chi_c_w", l0, l1, chi_c_b)
        del chi_c_b
        W_b -= 0.25 * (SpC[None] @ cb @ SpC[None])
        store.put_freq_batch("W", l0, l1, W_b)
        del W_b
    return DressResult(collapse0=collapse0, collapse_s=collapse_s, collapse_c=collapse_c,
                       static_s=static_s, static_c=static_c,
                       cond_min_s=float(cond_s), cond_min_c=float(cond_c))


# =============================================================================
# Bond-aware self-energy transport (spec 3.4)
# =============================================================================

_TRANSPORT_COST_WARN = 1.0e11


def calc_self_energy_bond(store, green_kw, beta, view, shape, norb, workers):
    """Sigma_fluct(k, iw) from the bond-resolved effective interaction ``W``
    in ``store`` (spec 3.4, normative equation with the momentum phases
    e^{+ik.R_alpha} e^{-i(k-q).R_beta} realised as real-space rolls):

        Sigma_ab(k) = T/N sum_{q,nu} sum_{alpha,beta,c,d}
            e^{+ik.R_alpha} W_{(alpha,c,a),(beta,d,b)}(q) e^{-i(k-q).R_beta} G_cd(k-q)

    ``green_kw`` is rank five `(1, nmat, nvol, norb, norb)`; the result has
    the same shape.  With a single on-site channel this reproduces
    `FLEX._calc_self_energy_general` byte for byte."""
    nx, ny, nz = (int(x) for x in shape)
    nvol = nx * ny * nz
    nmat = green_kw.shape[1]
    P = norb
    nd = norb * norb
    B = view.n_channels
    cost = float(B * B) * nmat * nvol * (P ** 2 * np.log2(max(nmat * nvol, 2)) + P ** 3)
    if cost > _TRANSPORT_COST_WARN:
        logger.warning("bond self-energy transport cost estimate %.2e exceeds %.0e "
                       "(B=%d, Nmat=%d, Nvol=%d, norb=%d); expect long iterations",
                       cost, _TRANSPORT_COST_WARN, B, nmat, nvol, norb)
    G_kw = green_kw[0]
    G_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw.reshape(nmat, nvol * P * P), axis=0).reshape(nmat, nx, ny, nz, P * P),
        axes=(1, 2, 3), workers=workers).reshape(nmat, nvol, P, P)
    Sigma_rt = np.zeros((nmat, nvol, P, P), dtype=np.complex128)
    axes = (1, 2, 3)
    for alpha in range(B):
        Ra = np.asarray(view.delta_r[alpha], dtype=int)
        for bt in range(B):
            Rb = np.asarray(view.delta_r[bt], dtype=int)
            blk = np.ascontiguousarray(store.get_pair("W", alpha, bt))               # (nmat, nvol, nd, nd)
            blk_qt = _ms.boson_to_tau(blk.reshape(nmat, nvol * nd * nd), axis=0)
            del blk
            Wab_rt = _bk.spatial_ifftn(blk_qt.reshape(nmat, nx, ny, nz, nd * nd),
                                       axes=axes, workers=workers)
            del blk_qt
            Wab_rt = np.roll(Wab_rt, tuple(-Rb), axis=axes).reshape(nmat, nvol, P, P, P, P)
            A = np.einsum('frcadb,frcd->frab', Wab_rt, G_rt)
            del Wab_rt
            A_rolled = np.roll(A.reshape(nmat, nx, ny, nz, P, P), tuple(Rb - Ra), axis=axes
                               ).reshape(nmat, nvol, P, P)
            del A
            Sigma_rt += A_rolled
            del A_rolled
    del G_rt
    tmp = _bk.spatial_fftn(Sigma_rt.reshape(nmat, nx, ny, nz, P * P), axes=axes, workers=workers)
    del Sigma_rt
    sigma = _ms.tau_to_fermion(tmp.reshape(nmat, nvol * P * P), axis=0)
    del tmp
    sigma = sigma.reshape(1, nmat, nvol, P, P)
    sigma *= 1.0 / beta
    return sigma
