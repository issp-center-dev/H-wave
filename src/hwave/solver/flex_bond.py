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

    @property
    def released(self):
        return self._released


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


from .hartree_fock import NonFiniteError as _NonFiniteError


def _at(iteration):
    return "" if iteration is None else " (SCF iteration {})".format(iteration)


def _dress(cb, V, channel, l0, nmat, spatial_shape, cond_tol, iteration):
    try:
        chi_b, cond = _bc.dress_batch(cb, V, channel, l0=l0, nmat=nmat, spatial_shape=spatial_shape,
                                      cond_tol=cond_tol)
    except ValueError as exc:
        if iteration is None:
            raise
        raise ValueError("{}{}".format(exc, _at(iteration))) from exc
    if not np.all(np.isfinite(chi_b)):
        raise _NonFiniteError("non-finite dressed {} channel in the frequency batch starting "
                              "at l={}{}".format(channel, l0, _at(iteration)))
    return chi_b, cond


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
                      cond_tol=_bc._BOND_COND_FLOOR, iteration=None):
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
        chi_s_b, cs = _dress(cb, S, "spin", l0, nmat, spatial_shape, cond_tol, iteration)
        cond_s = min(cond_s, cs if cs is not None else np.inf)
        W_b = 1.5 * (S[None] @ chi_s_b @ S[None])
        collapse_s[l0:l1] = chi_s_b[:, :, :nd, :nd]
        if l0 <= l_static < l1:
            static_s[...] = chi_s_b[l_static - l0]
        if output_full:
            store.put_freq_batch("chi_s_w", l0, l1, chi_s_b)
        del chi_s_b
        chi_c_b, cc = _dress(cb, C, "charge", l0, nmat, spatial_shape, cond_tol, iteration)
        cond_c = min(cond_c, cc if cc is not None else np.inf)
        W_b += 0.5 * (C[None] @ chi_c_b @ C[None])
        collapse_c[l0:l1] = chi_c_b[:, :, :nd, :nd]
        if l0 <= l_static < l1:
            static_c[...] = chi_c_b[l_static - l0]
        if output_full:
            store.put_freq_batch("chi_c_w", l0, l1, chi_c_b)
        del chi_c_b
        W_b -= 0.25 * (SpC[None] @ cb @ SpC[None])
        if not np.all(np.isfinite(W_b)):
            raise _NonFiniteError("non-finite effective interaction W in the frequency batch "
                                  "[{}, {}){}".format(l0, l1, _at(iteration)))
        store.put_freq_batch("W", l0, l1, W_b)
        del W_b
    return DressResult(collapse0=collapse0, collapse_s=collapse_s, collapse_c=collapse_c,
                       static_s=static_s, static_c=static_c,
                       cond_min_s=float(cond_s), cond_min_c=float(cond_c))


# =============================================================================
# Bond-aware self-energy transport (spec 3.4)
# =============================================================================



def calc_self_energy_bond(store, green_kw, beta, view, shape, norb, workers):
    """Sigma_fluct(k, iw) from the bond-resolved effective interaction ``W``
    in ``store`` (spec 3.4, normative equation):

        Sigma_ab(k) = T/N sum_{q,nu} sum_{alpha,beta,c,d}
            e^{+i(k-q).(R_alpha - R_beta)} W_{(alpha,c,a),(beta,d,b)}(q) G_cd(k-q)

    The bond form factor sits on the internal leg k - q on BOTH ends of W
    (the bubble's pair operator carries its phase on the reversed
    propagator, so chibar_{alpha beta}(q) = -(T/N) sum_k e^{i(k-q).(R_alpha
    - R_beta)} G(k) G(k-q)); with this assignment chibar(q, i nu)^dagger =
    chibar(q, -i nu) and W inherit the plain matrix conjugation symmetry
    and Sigma(k, i w)^dagger = Sigma(k, -i w) holds to round-off.  In real
    space the phase is G(r + R_alpha - R_beta), a roll of G by R_beta - R_alpha.

    ``green_kw`` is rank five `(1, nmat, nvol, norb, norb)`; the result has
    the same shape.  With a single on-site channel this reproduces
    `FLEX._calc_self_energy_general` byte for byte."""
    nx, ny, nz = (int(x) for x in shape)
    nvol = nx * ny * nz
    nmat = green_kw.shape[1]
    P = norb
    nd = norb * norb
    B = view.n_channels
    G_kw = green_kw[0]
    G_rt = _bk.spatial_ifftn(
        _ms.fermion_to_tau(G_kw.reshape(nmat, nvol * P * P), axis=0).reshape(nmat, nx, ny, nz, P * P),
        axes=(1, 2, 3), workers=workers).reshape(nmat, nx, ny, nz, P, P)
    Sigma_rt = np.zeros((nmat, nvol, P, P), dtype=np.complex128)
    axes = (1, 2, 3)
    for alpha in range(B):
        Ra = np.asarray(view.delta_r[alpha], dtype=int)
        for bt in range(B):
            Rb = np.asarray(view.delta_r[bt], dtype=int)
            shift = tuple(int(x) for x in (Rb - Ra))   # e^{+ik'.(R_a - R_b)} G(k') = G(r + R_a - R_b)
            if shift == (0, 0, 0):
                G_sh = G_rt.reshape(nmat, nvol, P, P)
            else:
                G_sh = np.roll(G_rt, shift, axis=axes).reshape(nmat, nvol, P, P)
            blk = np.ascontiguousarray(store.get_pair("W", alpha, bt))               # (nmat, nvol, nd, nd)
            blk_qt = _ms.boson_to_tau(blk.reshape(nmat, nvol * nd * nd), axis=0)
            del blk
            Wab_rt = _bk.spatial_ifftn(blk_qt.reshape(nmat, nx, ny, nz, nd * nd),
                                       axes=axes, workers=workers).reshape(nmat, nvol, P, P, P, P)
            del blk_qt
            A = np.einsum('frcadb,frcd->frab', Wab_rt, G_sh)
            del Wab_rt, G_sh
            Sigma_rt += A
            del A
    del G_rt
    tmp = _bk.spatial_fftn(Sigma_rt.reshape(nmat, nx, ny, nz, P * P), axes=axes, workers=workers)
    del Sigma_rt
    sigma = _ms.tau_to_fermion(tmp.reshape(nmat, nvol * P * P), axis=0)
    del tmp
    sigma = sigma.reshape(1, nmat, nvol, P, P)
    sigma *= 1.0 / beta
    return sigma


# =============================================================================
# Memory preflight (spec 3.6) and output-path resolution (spec 4.2)
# =============================================================================

import os as _os

_DRESSING_OPS_WARN = 1.0e12
_TRANSPORT_OPS_WARN = 1.0e11
_GIB = float(1024 ** 3)


def dressing_ops(nmat, nvol, ND):
    """Operation count of the two dense ND x ND solves over all (l, q)."""
    return 2.0 * nmat * nvol * float(ND) ** 3


def transport_ops(B, nmat, nvol, norb):
    """Operation count of the bond self-energy transport (transforms plus
    the orbital contraction, P^3 = norb^6)."""
    P = norb ** 2
    return float(B * B) * nmat * nvol * (P ** 2 * np.log2(max(nmat * nvol, 2)) + P ** 3)


def estimate_bond_memory(*, nmat, nvol, norb, B, depth, output_full, split_seed, n_types,
                         freq_batch, cap_gb, mixing):
    """The named-buffer lifetime table of spec 3.6 (every row raw), the
    batch selection and the admission decision against ``cap_gb`` (binary
    GiB). Returns a dict with ``persistent_rows``, ``phase_rows`` (at the
    selected ``nb``), ``persistent``, ``nb``, ``peak`` and the symbols;
    raises ``ValueError`` naming every row when even ``nb = 1`` exceeds
    the cap, or when ``freq_batch`` does."""
    nmat, nvol, norb, B = int(nmat), int(nvol), int(norb), int(B)
    depth = int(depth) if mixing == "anderson" else 0
    it = 16
    P = norb * norb
    ND = B * P
    U = nmat * nvol * ND * ND * it
    G = nmat * nvol * P * it
    C = nmat * nvol * P * P * it
    S = nvol * ND * ND * it
    H = nvol * P * it
    persistent_rows = {
        "chibar_W": 2 * U,
        "chi_sc_w": 2 * U if output_full else 0,
        "vertices_static": 5 * S,
        "state": G + H,
        "seed_envelope": (2 * G + H) if split_seed else 0,
        "anderson_history": 2 * depth * 2 * G,
        "collapses": 3 * C,
        "eigenpairs_hf": 5 * H,
        "flex_arrays": 5 * G,
        "hf_tables": int(n_types) * H,
    }
    persistent = sum(persistent_rows.values())
    prep = it * nvol * (ND ** 2 + 4 * P * (nmat + 2))
    pair = it * nvol * (2 * ND ** 2 + 3 * nmat * P ** 2 + 2 * nmat * P + 8 * P)
    per_batch = 6 * nvol * ND * ND * it

    def _phase_rows(nb):
        return {
            "green_mu": 6 * G + H,
            "density_hf": 2 * H + 8 * H + 20 * H + H,
            "bubble": max(prep, pair),
            "dressing": per_batch * nb,
            "transport": 4 * G + 4 * C,
            "convergence": 6 * G + H,
            "mixing": (2 * depth * 2 * G + 4 * G) if mixing == "anderson" else 2 * G,
            "post_scf": G + (U if output_full else C),
        }

    def _peak(nb):
        return 1.25 * (persistent + max(_phase_rows(nb).values()))

    cap_bytes = float(cap_gb) * _GIB

    def _table(nb):
        rows = ["  persistent {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                for k, v in persistent_rows.items()]
        rows += ["  phase      {:>18s}: {:10.4f} GiB".format(k, v / _GIB)
                 for k, v in _phase_rows(nb).items()]
        return "\n".join(rows)

    if freq_batch is not None:
        nb = int(freq_batch)
        if not 1 <= nb <= nmat:
            raise ValueError("longitudinal_bond_freq_batch must be in [1, Nmat={}], got {}"
                             .format(nmat, nb))
        if _peak(nb) > cap_bytes:
            raise ValueError(
                "[mode.param] longitudinal_bond_freq_batch = {} gives an estimated peak of "
                "{:.4f} GiB = 1.25 * (persistent + max phase row) above "
                "longitudinal_bond_memory_cap_gb = {:.4f} GiB (B = {}, ND = {}, nvol = {}, "
                "Nmat = {}); the rows are\n{}\nReduce the batch, the k mesh or Nmat, drop "
                "declared-zero outer shells with longitudinal_bond_max_shells, or raise the cap."
                .format(nb, _peak(nb) / _GIB, cap_bytes / _GIB, B, ND, nvol, nmat, _table(nb)))
    else:
        if _peak(1) > cap_bytes:
            raise ValueError(
                "[mode.param] longitudinal_bond_channels=true: the estimated peak host memory "
                "{:.4f} GiB (frequency batch 1) = 1.25 * (persistent + max phase row) exceeds "
                "longitudinal_bond_memory_cap_gb = {:.4f} GiB (B = {}, ND = {}, nvol = {}, "
                "Nmat = {}); the rows are\n{}\nReduce the k mesh or Nmat, drop declared-zero "
                "outer shells with longitudinal_bond_max_shells, disable "
                "longitudinal_bond_output_full, lower the Anderson depth, or raise the cap."
                .format(_peak(1) / _GIB, cap_bytes / _GIB, B, ND, nvol, nmat, _table(1)))
        nb = 1
        for cand in range(nmat, 0, -1):
            if _peak(cand) <= cap_bytes:
                nb = cand
                break
    phase_rows = _phase_rows(nb)
    return dict(persistent_rows=persistent_rows, phase_rows=phase_rows, persistent=persistent,
                nb=nb, peak=_peak(nb), cap_bytes=cap_bytes, U=U, G_bytes=G, C_bytes=C,
                S_bytes=S, H_bytes=H, B=B, ND=ND, nvol=nvol, nmat=nmat, table=_table(nb),
                dressing_ops=dressing_ops(nmat, nvol, ND),
                transport_ops=transport_ops(B, nmat, nvol, norb))


_NPZ_ARTIFACTS = ("chi0q", "chiq_s", "chiq_c", "chiq", "sigma", "green", "longitudinal_bond")
_DEFAULT_FILES = {"chiq_s": "chiq_s", "chiq_c": "chiq_c",
                  "longitudinal_bond": "longitudinal_bond.npz"}


def resolve_output_paths(info_outputfile, path_to_output, active_outputs):
    """``{name: absolute path}`` for the artifacts in ``active_outputs``
    that will actually be written (spec 4.2): joins with ``path_to_output``,
    applies the ``.npz`` suffix rule of ``numpy.savez`` to the archive
    artifacts only, normalises to absolute paths (symlink equivalence is
    not checked) and refuses two artifacts resolving to one file."""
    out = {}
    for name in active_outputs:
        fn = info_outputfile.get(name, _DEFAULT_FILES.get(name))
        if fn is None:
            continue
        p = _os.path.join(str(path_to_output), str(fn))
        if name in _NPZ_ARTIFACTS and not p.endswith(".npz"):
            p += ".npz"
        out[name] = _os.path.abspath(_os.path.normpath(p))
    seen = {}
    for name, p in out.items():
        if p in seen:
            raise ValueError(
                "[file.output] the artifacts '{}' and '{}' resolve to the same file {}; "
                "give them distinct file names (numpy appends '.npz' to archive names "
                "without that suffix).".format(seen[p], name, p))
        seen[p] = name
    return out
