"""Bond-resolved, frequency-resolved Eliashberg pairing kernel built from the
FLEX bond gate's own vertices and dressed susceptibilities (design spec
2026-09-17-eliashberg-bond-dynamic-design.md). Array-only: no TOML, no files."""
import logging
from dataclasses import dataclass

import numpy as np

from . import backend as _bk
from . import bond_channels as _bc
from . import flex_bond as _fb
from . import matsubara as _ms

logger = logging.getLogger("qlms").getChild("eliashberg_bond")

#: The pairing-type coefficients of the fluctuation vertex (spec 4.2):
#: ``Gamma_eta = c_s S chi_s S + c_c C chi_c C``.
COEFF = {"singlet": (1.5, -0.5), "triplet": (-0.5, -0.5)}

_GIB = 1024.0 ** 3


class ArrayBlockSource:
    """Host arrays ``{name: (nmat, nvol, ND, ND)}`` behind the block-source
    contract of :class:`~hwave.solver.flex_bond.BondBlockStore`."""

    def __init__(self, arrays, nd):
        self.nd = int(nd)
        self._a = {}
        shape = None
        for k, v in arrays.items():
            v = np.asarray(v)
            if v.ndim != 4 or v.shape[2] != v.shape[3] or v.shape[2] % self.nd != 0:
                raise ValueError("ArrayBlockSource: {!r} must be (nmat, nvol, ND, ND) with ND a "
                                 "multiple of nd = {}, got {}".format(k, self.nd, v.shape))
            if shape is not None and v.shape != shape:
                raise ValueError("ArrayBlockSource: shapes differ: {} vs {}".format(shape, v.shape))
            shape = v.shape
            self._a[k] = v
        if shape is None:
            raise ValueError("ArrayBlockSource: no arrays")
        self.nmat, self.nvol, self.ND = int(shape[0]), int(shape[1]), int(shape[2])
        self.B = self.ND // self.nd

    @property
    def names(self):
        return tuple(self._a)

    def get_freq_batch(self, name, l0, l1):
        return self._a[name][l0:l1]

    def put_freq_batch(self, name, l0, l1, block):
        self._a[name][l0:l1] = block

    def get_pair(self, name, alpha, beta):
        nd = self.nd
        return self._a[name][:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd]


@dataclass(frozen=True)
class PairVertexUniform:
    """``Gamma`` on the uniform bosonic grid, living in ``source[slot]``."""
    source: object
    slot: str
    B: int
    nd: int
    nmat: int
    nvol: int


@dataclass(frozen=True)
class PairVertexIR:
    """``Gamma`` compressed onto the bosonic IR basis (spec 4.5)."""
    coeffs: np.ndarray              # (L_B, nvol, ND, ND) host
    const: object                   # (nvol, ND, ND) host, only with ir_keep_static
    axB: object
    fit_residual_rel: np.ndarray    # (B, B) componentwise relative residual


class PairVertexAccumulator:
    """Staged, ADDITIVE construction of ``Gamma = c_s S chi_s S + c_c C chi_c C``
    for one or two pairing types (spec 4.2, 4.5, 5.1), so that the spin and the
    charge susceptibility never need to be resident together.

    ``ir is None``: exactly ONE pairing type; the target is
    ``out_store[out_slot]`` (host), zeroed on construction. ``ir = (axF, axB)``:
    one OR two pairing types; the targets are the per-type, per-stage
    coefficient arrays of spec 4.5 and nothing is written to any store.
    """

    _STAGES = ("spin", "charge")

    def __init__(self, dev, *, pairing_types, nb, nmat, nvol, nd, spatial_shape, ir=None,
                 out_store=None, out_slot="W", ir_keep_static=False, workers=1):
        pairing_types = tuple(pairing_types)
        if not pairing_types:
            raise ValueError("PairVertexAccumulator: no pairing type requested")
        for eta in pairing_types:
            if eta not in COEFF:
                raise ValueError("unknown pairing type {!r}".format(eta))
        if len(set(pairing_types)) != len(pairing_types):
            raise ValueError("PairVertexAccumulator: duplicate pairing types {}"
                             .format(pairing_types))
        self.dev, self.xp = dev, dev.xp
        self.pairing_types = pairing_types
        self.nb, self.nmat, self.nvol, self.nd = int(nb), int(nmat), int(nvol), int(nd)
        if self.nb <= 0:
            raise ValueError("PairVertexAccumulator: nb must be positive, got {}".format(self.nb))
        self.shape = tuple(int(x) for x in spatial_shape)
        self.ND = int(dev.S.shape[-1])
        if self.ND % self.nd != 0:
            raise ValueError("PairVertexAccumulator: the vertex pair dimension ND = {} is not a "
                             "multiple of the channel block size nd = {}".format(self.ND, self.nd))
        self.B = self.ND // self.nd
        self.workers = workers
        self.ir = ir
        self.ir_keep_static = bool(ir_keep_static)
        self._done = set()
        self._resid_mode = False
        self._resid_type = None
        self._resid_seen = set()
        if ir is None:
            if len(self.pairing_types) != 1:
                raise ValueError("PairVertexAccumulator: the uniform grid takes exactly one "
                                 "pairing type, got {}".format(self.pairing_types))
            if out_store is None:
                out_store = ArrayBlockSource(
                    {out_slot: np.zeros((self.nmat, self.nvol, self.ND, self.ND), np.complex128)},
                    self.nd)
            else:
                for l0 in range(0, self.nmat, self.nb):
                    l1 = min(self.nmat, l0 + self.nb)
                    out_store.get_freq_batch(out_slot, l0, l1)[...] = 0.0
            self._target, self._slot = out_store, out_slot
        else:
            axF, axB = ir
            self.axF, self.axB = axF, axB
            fit, ev = axB.uniform_matrices(self.nmat, with_constant=True)
            self._fit, self._ev = np.asarray(fit), np.asarray(ev)      # (nmat, L+1), (L, nmat)
            self.L = int(axB.L)
            self._sol = {eta: {st: np.zeros((self.L + 1, self.nvol, self.ND, self.ND),
                                            np.complex128)
                               for st in self._STAGES} for eta in self.pairing_types}
            self._r = {st: np.zeros((self.B, self.B)) for st in self._STAGES}
            self._a = {st: np.zeros((self.B, self.B)) for st in self._STAGES}

    # -- stages -----------------------------------------------------------
    def add_dressed(self, source, *, cond_tol=_bc._BOND_COND_FLOOR, guard_freqs="all"):
        """In-process stage: per frequency batch of ``source['chibar']`` re-dress
        the spin AND the charge channel with :func:`flex_bond._dress` and add
        both contributions for every requested pairing type."""
        xp, S, C = self.xp, self.dev.S, self.dev.C
        for l0 in range(0, self.nmat, self.nb):
            l1 = min(self.nmat, l0 + self.nb)
            cb = _bk.to_device(source.get_freq_batch("chibar", l0, l1), xp)
            chi_s_b, _ = _fb._dress(cb, S, "spin", l0, self.nmat, self.shape, cond_tol,
                                    None, guard_freqs)
            self._absorb("spin", S[None] @ chi_s_b @ S[None], l0, l1)
            del chi_s_b
            chi_c_b, _ = _fb._dress(cb, C, "charge", l0, self.nmat, self.shape, cond_tol,
                                    None, guard_freqs)
            self._absorb("charge", C[None] @ chi_c_b @ C[None], l0, l1)
            del chi_c_b, cb
        self._done.update(self._STAGES)

    def add_channel(self, channel, source, name):
        """Post-processing stage: add ``c_s S chi S`` (``channel = "spin"``) or
        ``c_c C chi C`` (``"charge"``) from the already dressed ``source[name]``.
        Called once per channel, with one archive member resident at a time."""
        if channel not in self._STAGES:
            raise ValueError("channel must be 'spin' or 'charge', got {!r}".format(channel))
        V = self.dev.S if channel == "spin" else self.dev.C
        for l0 in range(0, self.nmat, self.nb):
            l1 = min(self.nmat, l0 + self.nb)
            chi_b = _bk.to_device(source.get_freq_batch(name, l0, l1), self.xp)
            self._absorb(channel, V[None] @ chi_b @ V[None], l0, l1)
            del chi_b
        self._done.add(channel)

    def _absorb(self, channel, P, l0, l1):
        xp = self.xp
        if not bool(xp.all(xp.isfinite(P))):
            raise ValueError("PairVertexAccumulator: non-finite {} contribution in the frequency "
                             "batch [{}, {})".format(channel, l0, l1))
        idx = 0 if channel == "spin" else 1
        Ph = _bk.to_host(P)
        for eta in self.pairing_types:
            if self._resid_mode and eta != self._resid_type:
                continue
            c = COEFF[eta][idx]
            if self.ir is None:
                self._target.get_freq_batch(self._slot, l0, l1)[...] += c * Ph
            elif not self._resid_mode:
                self._sol[eta][channel] += np.einsum("lvij,lc->cvij", c * Ph, self._fit[l0:l1])
            else:
                self._residual(eta, channel, c * Ph, l0, l1)

    def _residual(self, eta, channel, G, l0, l1):
        """Componentwise per-block misfit of this stage's own fit against this
        stage's own contribution (spec 4.5): ``eval(sol)`` reconstructs exactly
        what the kernel will use -- the smooth coefficients alone when the
        constant is dropped, plus the constant when it is retained."""
        self._resid_seen.add(channel)
        sol = self._sol[eta][channel]
        rec = np.einsum("cvij,cl->lvij", sol[:self.L], self._ev[:, l0:l1])
        if self.ir_keep_static:
            rec = rec + sol[self.L][None]
        nd, B = self.nd, self.B
        for al in range(B):
            for be in range(B):
                sl = (slice(None), slice(None),
                      slice(al * nd, (al + 1) * nd), slice(be * nd, (be + 1) * nd))
                self._r[channel][al, be] = max(self._r[channel][al, be],
                                               float(np.abs(rec[sl] - G[sl]).max()))
                self._a[channel][al, be] = max(self._a[channel][al, be],
                                               float(np.abs(G[sl]).max()))

    # -- finish -------------------------------------------------------------
    def finish(self, *, ir_fit_tol=0.5, stage_callable=None):
        """Close the accumulation. Uniform: a :class:`PairVertexUniform` over the
        target slot. IR: the residual pass of spec 4.5 (replayed through
        ``stage_callable``), the ``ir_fit_tol`` refusal / warning band and the
        frozen :class:`PairVertexIR` per pairing type."""
        for st in self._STAGES:
            if st not in self._done:
                raise ValueError("PairVertexAccumulator: the {} stage was never added".format(st))
        if self.ir is None:
            return {self.pairing_types[0]: PairVertexUniform(
                source=self._target, slot=self._slot, B=self.B, nd=self.nd,
                nmat=self.nmat, nvol=self.nvol)}
        if ir_fit_tol < 0:
            raise ValueError("ir_fit_tol must be >= 0")
        rel = np.full((self.B, self.B), np.nan)
        if ir_fit_tol > 0:
            if stage_callable is None:
                raise ValueError("PairVertexAccumulator.finish: the IR residual pass needs "
                                 "stage_callable (or ir_fit_tol = 0 to skip it)")
            # the residual pass is evaluated on the FIRST pairing type only: the
            # types differ by the scalar coefficients, which cancel in r / a
            self._resid_mode = True
            self._resid_type = self.pairing_types[0]
            self._resid_seen = set()
            try:
                stage_callable(self)
            finally:
                self._resid_mode = False
            # a callable that replays nothing would leave r = a = 0, which reads
            # as a perfect fit AND disables the constant-vs-scale refusal below
            missing = [st for st in self._STAGES if st not in self._resid_seen]
            if len(missing) == len(self._STAGES):
                raise ValueError("PairVertexAccumulator.finish: stage_callable replayed no "
                                 "stage; the IR residual cannot be evaluated")
            if missing:
                raise ValueError("PairVertexAccumulator.finish: stage_callable did not replay "
                                 "the {} stage; the IR residual cannot be evaluated"
                                 .format(missing[0]))
            rel = np.zeros((self.B, self.B))
            for st in self._STAGES:
                # a stage with a = 0 (identically zero contribution) counts as 0
                q = np.where(self._a[st] > 0, self._r[st] / np.where(self._a[st] > 0,
                                                                    self._a[st], 1.0), 0.0)
                rel = np.maximum(rel, q)
            worst = float(rel.max())
            if worst > ir_fit_tol:
                al, be = np.unravel_index(int(rel.argmax()), rel.shape)
                raise ValueError(
                    "IR fit of the bond pairing vertex: componentwise relative residual {:.3e} in "
                    "block ({}, {}) exceeds ir_fit_tol = {:.1e}; raise ir_wmax, increase Nmat of "
                    "the FLEX run, or set ir_fit_tol = 0 to skip the check"
                    .format(worst, al, be, ir_fit_tol))
            if worst >= 0.1 * ir_fit_tol:
                logger.warning("IR fit of the bond pairing vertex: componentwise relative residual "
                               "%.3e is within a decade of ir_fit_tol = %.1e", worst, ir_fit_tol)
            logger.info("IR fit of the bond pairing vertex: max componentwise relative "
                        "residual %.3e", worst)
        else:
            logger.warning("IR fit of the bond pairing vertex: residual check SKIPPED "
                           "(ir_fit_tol = 0)")
        out = {}
        for eta in self.pairing_types:
            sol = self._sol[eta]["spin"]
            sol += self._sol[eta]["charge"]
            self._sol[eta]["charge"] = None
            const = np.ascontiguousarray(sol[self.L]) if self.ir_keep_static else None
            cmax = float(np.abs(sol[self.L]).max())
            if ir_fit_tol > 0:
                scale = float(max(self._a["spin"].max(), self._a["charge"].max()))
                if not self.ir_keep_static and cmax > scale > 0:
                    raise ValueError(
                        "IR compress of the bond pairing vertex ({}): the discarded frequency-"
                        "independent component ({:.3e}) exceeds the data scale ({:.3e}); lower "
                        "ir_wmax, increase Nmat, or set [eliashberg] ir_keep_static_chi = true"
                        .format(eta, cmax, scale))
            logger.info("IR compress of the bond pairing vertex (%s): frequency-independent "
                        "component %.3e %s", eta, cmax,
                        "retained" if self.ir_keep_static else "discarded")
            out[eta] = PairVertexIR(coeffs=np.ascontiguousarray(sol[:self.L]), const=const,
                                    axB=self.axB, fit_residual_rel=rel)
        return out


# =============================================================================
# Memory model and admission (spec 8)
# =============================================================================

@dataclass
class AdmissionTable:
    """The incremental byte rows of the pairing step and the per-residency
    ``(host, device)`` needs derived from them (spec 8)."""
    rows: dict
    residency_need: dict          # residency -> (host_bytes, device_bytes)
    requested: str

    def need(self, residency):
        return self.residency_need[residency]

    def choose(self, requested, host_cap, device_cap):
        order = ("device", "host", "stream") if requested == "auto" else (requested,)
        for r in order:
            h, d = self.residency_need[r]
            if h <= host_cap and d <= device_cap:
                return r
        raise MemoryError(
            "bond pairing kernel: no residency fits (requested {}; host cap {:.3f} GiB, device cap "
            "{:.3f} GiB); needs host/device GiB: {}; rows (GiB): {}".format(
                requested, host_cap / _GIB, device_cap / _GIB,
                {r: (round(h / _GIB, 3), round(d / _GIB, 3))
                 for r, (h, d) in self.residency_need.items()},
                {k: round(v / _GIB, 3) for k, v in self.rows.items()}))


def estimate_pair_memory(*, nmat, ntau, nfreq, nvol, norb, B, num_eigenvalues, residency, ir,
                         L_B, n_channels, in_process, nb=0):
    """The spec-8 incremental memory table of the pairing step. ``L_B`` is the
    bosonic coefficient count (``0`` on the uniform grid), ``n_channels``
    (1 or 2) the number of pairing types whose IR coefficient rows are alive,
    and ``nb`` the dressing batch size of the vertex build (``0`` for a
    kernel-only admission, where no dressing batch is allocated)."""
    nd = norb * norb
    ND = B * nd
    S = ntau if ir else nmat
    vec = nd * nvol * nfreq * 16
    blk = S * nvol * nd * nd * 16
    ncv = 2 * int(num_eigenvalues) + 1
    rows = {
        "vertex_slot": 0 if (ir or in_process) else nmat * nvol * ND * ND * 16,
        "coefficients_build": (2 * (L_B + 1) * nvol * ND * ND * 16 * n_channels) if ir else 0,
        "coefficients": ((L_B + 1) * nvol * ND * ND * 16 * n_channels) if ir else 0,
        # one dressing batch on the device during the vertex build (the FLEX
        # dressing row); the post-processing entry holds one archive member
        # on the host while it feeds that build
        "dressing_workspace": 7 * int(nb) * nvol * ND * ND * 16,
        "source_member": 0 if in_process else nmat * nvol * ND * ND * 16,
        "G2": norb ** 4 * nvol * nfreq * 16,
        "gap_work": (3 + B) * vec,
        "hoisted_blocks": B * B * blk,
        "stream_workspace": 3 * blk,
        "eigen_vectors": (ncv + 2) * vec,
    }
    build_host = max(rows["vertex_slot"], rows["coefficients_build"]) + rows["source_member"]
    # the dressing batch is a BUILD-phase device peak, so it is part of every
    # residency's device need, not only the streaming one
    dev_common = (rows["G2"] + rows["gap_work"] + rows["stream_workspace"]
                  + rows["dressing_workspace"])
    need = {
        "device": (int(1.25 * (build_host + rows["eigen_vectors"])),
                   int(1.25 * (dev_common + rows["hoisted_blocks"]))),
        "host": (int(1.25 * (build_host + rows["hoisted_blocks"] + rows["eigen_vectors"])),
                 int(1.25 * dev_common)),
        "stream": (int(1.25 * (build_host + rows["eigen_vectors"])), int(1.25 * dev_common)),
    }
    return AdmissionTable(rows=rows, residency_need=need, requested=residency)


# =============================================================================
# The pairing kernel on the uniform grid (spec 3, 4.1-4.4)
# =============================================================================

def instantaneous_vertex(S, C, nd, pairing_type, spatial_shape):
    """The frequency-INDEPENDENT pairing vertex of spec 4.2.

    ``0.5 (S_0 + C_0)`` (singlet) / ``0.5 (C_0 - S_0)`` (triplet) built from the
    CHANNEL-0 blocks of the enlarged bond vertices and crossed into the on-site
    ``(l1, l2, l3, l4)`` layout by the shared
    :func:`~hwave.solver.eliashberg_dynamic._instantaneous_vertex` (which
    evaluates ``sc._compute_vertices_flex`` at ``chi_s = chi_c = 0``), so the
    bond path and the on-site path carry bit-identical bare terms.

    Parameters
    ----------
    S, C : ndarray, shape (nvol, ND, ND)
        The enlarged bond spin/charge vertices; only the channel-0 block
        ``[:, :nd, :nd]`` is read.
    nd : int
        ``norb**2``.
    pairing_type : {"singlet", "triplet"}
    spatial_shape : tuple
        ``(Nx, Ny, Nz)``.

    Returns
    -------
    ndarray, shape ``(norb, norb, norb, norb, Nx, Ny, Nz)``
    """
    from . import eliashberg_dynamic as _ed
    nx, ny, nz = (int(x) for x in spatial_shape)
    norb = int(round(nd ** 0.5))
    S0 = np.asarray(S)[:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
    C0 = np.asarray(C)[:, :nd, :nd].reshape(nx, ny, nz, nd, nd)
    return _ed._instantaneous_vertex({}, norb, nx, ny, nz, pairing_type=pairing_type,
                                     convention="myo", sc_matrices=(S0, C0))


def _vinst_as_block(V_inst, nd, nvol):
    """``(norb,)*4 x (Nx, Ny, Nz)`` -> ``(nvol, nd, nd)`` pair-matrix form
    (row pair ``(a, b)``, column pair ``(c, d)``) -- the layout the vertex
    blocks live in, so the bare term can be added to the channel-0 block."""
    return np.ascontiguousarray(
        np.asarray(V_inst).transpose(4, 5, 6, 0, 1, 2, 3).reshape(nvol, nd, nd))


def block_to_rtau_uniform(blk, spatial_shape, norb, workers, xp):
    """One ``(alpha, beta)`` vertex block ``(nmat, nvol, nd, nd)`` on the uniform
    bosonic grid -> ``(a, b, c, d, x, y, z, tau)`` on the module ``xp``.

    The same two transforms the on-site vertex leg
    (:func:`~hwave.solver.eliashberg_dynamic.vertex_qw_to_rt`) applies -- boson
    -> tau on the frequency axis, ``ifftn`` (which carries the single spatial
    fold's ``1/N``) on the spatial axes -- written for the block's own axis
    order and transposed into the kernel's ``(a, b, c, d, x, y, z, t)`` layout.
    """
    nx, ny, nz = spatial_shape
    nmat = blk.shape[0]
    nvol = nx * ny * nz
    nd = norb * norb
    b = _bk.to_device(np.ascontiguousarray(blk) if isinstance(blk, np.ndarray) else blk, xp)
    bt = _ms.boson_to_tau(b.reshape(nmat, nvol * nd * nd), axis=0)
    br = _bk.spatial_ifftn(bt.reshape(nmat, nx, ny, nz, nd * nd), axes=(1, 2, 3),
                           workers=workers)
    return xp.ascontiguousarray(
        br.reshape(nmat, nx, ny, nz, norb, norb, norb, norb).transpose(4, 5, 6, 7, 1, 2, 3, 0))


def _host_available_bytes():
    """Host memory the admission model may spend (the shared probe)."""
    from . import eliashberg_dynamic as _ed
    return _ed._available_ram_bytes()


class BondPairKernel:
    r"""The bond-resolved, frequency-resolved linearized Eliashberg operator
    (spec 3, 4.1-4.4) as a ``matvec`` on a flat gap.

    The gap is ``(norb, norb, Nx, Ny, Nz, nfreq)``; the bond channel index is
    INTERNAL (summed over, never part of the gap). One matvec is

    .. code-block:: text

        F(k, n)        = sum G2(k, n) phi(k, n)                 (orbital pair)
        F_beta(r, tau) = roll(ifft F(r, tau), -R_beta)
        P_alpha        = - sum_beta Gamma_{alpha beta}(r, tau) F_beta(r, tau)
        out(k, n)      = tau->freq fft sum_alpha roll(P_alpha, -R_alpha)

    i.e. ``2B`` rolls and ``B**2`` real-space multiply-accumulates per matvec.
    The roll convention is spec 4.3: ``e^{+ik.R} f(k)`` is ``f(r + R)``, which
    on the FFT grid is ``xp.roll(f, -R, axis=(2, 3, 4))`` -- the same rule the
    static bond kernel (:func:`~hwave.solver.bond_channels.make_bond_kernel`)
    realizes with explicit ``e^{+ik.dr_m}`` phases and the FLEX bond transport
    (``flex_bond.calc_self_energy_bond``) realizes with its ``R_alpha -
    R_beta`` shift.

    The frequency-independent bare vertex ``V_inst`` (spec 4.2) is added to the
    CHANNEL-0 block at every bosonic frequency: on the dense uniform tau grid
    that is exactly the ``delta(tau)`` its Fourier transform represents.

    Normalization follows the on-site uniform kernel
    (:func:`~hwave.solver.eliashberg_dynamic.eliashberg_kernel_dynamic`): the
    ``1/beta`` lives inside ``G2``, the ``-(1/N)`` in the single ``ifftn``, and
    no explicit ``beta`` is applied. The IR branch follows
    :func:`~hwave.solver.eliashberg_dynamic.eliashberg_kernel_ir` instead (the
    IR transforms are physical, so the operator carries one explicit ``beta``).

    Parameters
    ----------
    vertex : PairVertexUniform or PairVertexIR
        The pairing vertex ``Gamma_eta`` (:class:`PairVertexAccumulator`).
    G2 : ndarray
        Pair bubble ``(norb, norb, norb, norb, Nx, Ny, Nz, nfreq)`` on the
        fermionic axis (uniform grid: ``eliashberg_dynamic.calc_g2_dynamic``).
    view : BondSetView
        The bond topology (``n_channels``, ``delta_r``).
    xp : module
        Array backend the matvec runs on (``numpy`` or ``cupy``).
    residency : {"auto", "device", "host", "stream"}
        Where the ``B**2`` transformed vertex blocks live. ``"device"`` and
        ``"host"`` hoist them once; ``"stream"`` rebuilds one block at a time
        inside every matvec. On the numpy backend there is no separate device,
        so ``"auto"`` and ``"device"`` both resolve to ``"host"`` and the
        device cap is the host cap.
    """

    def __init__(self, vertex, G2, view, *, xp, spatial_shape, norb, beta, nfreq, V_inst=None,
                 axF=None, residency="auto", admission=None, host_cap=None, device_cap=None,
                 workers=1):
        self.xp = xp
        self.shape = tuple(int(x) for x in spatial_shape)
        nx, ny, nz = self.shape
        self.nvol = nx * ny * nz
        self.norb, self.nd = int(norb), int(norb) ** 2
        self.beta, self.nfreq, self.workers = float(beta), int(nfreq), workers
        self.view = view
        self.B = int(view.n_channels)
        self.delta_r = [tuple(int(x) for x in r) for r in view.delta_r]
        self.vertex = vertex
        self.is_ir = isinstance(vertex, PairVertexIR)
        self.axF = axF
        if self.is_ir and axF is None:
            raise ValueError("BondPairKernel: the IR vertex needs the fermionic axis axF")
        self.G2 = xp.asarray(G2)
        self.gap_shape = (self.norb, self.norb, nx, ny, nz, self.nfreq)
        self.V_inst = None if V_inst is None else np.asarray(V_inst)
        self._Vinst_rt = None
        self._const_rt = None

        # -- admission (spec 8) --------------------------------------------
        if admission is None:
            admission = estimate_pair_memory(
                nmat=0 if self.is_ir else self.nfreq,
                ntau=(len(axF.tau) if self.is_ir else 0), nfreq=self.nfreq, nvol=self.nvol,
                norb=self.norb, B=self.B, num_eigenvalues=10, residency=residency,
                ir=self.is_ir, L_B=(int(vertex.coeffs.shape[0]) if self.is_ir else 0),
                n_channels=1, in_process=False)
        if host_cap is None:
            host_cap = 0.8 * _host_available_bytes()
        req = residency
        if xp is np:
            # no separate device: "auto"/"device" mean the hoisted host blocks,
            # and the device rows are budgeted against the host cap
            if req in ("auto", "device"):
                req = "host"
            device_cap = host_cap
        elif device_cap is None:
            device_cap = 0.9 * _bk.device_available_bytes()
        self.residency = admission.choose(req, host_cap, device_cap)
        self.admission = admission

        # -- hoisting -------------------------------------------------------
        self._blocks = None
        if self.residency in ("device", "host"):
            mod = xp if self.residency == "device" else np
            self._blocks = {(a, b): self._block_rtau(a, b, mod)
                            for a in range(self.B) for b in range(self.B)}
        if self.is_ir:
            if self.V_inst is not None:
                self._Vinst_rt = xp.asarray(_bk.spatial_ifftn(
                    self.V_inst.astype(complex), axes=(4, 5, 6), workers=workers))
            if vertex.const is not None:
                self._const_rt = {(a, b): self._const_block_r(a, b, xp)
                                  for a in range(self.B) for b in range(self.B)}

    # -- block builders ------------------------------------------------------
    def _raw_block_uniform(self, a, b):
        """The ``(a, b)`` vertex block on the bosonic grid, with the bare term
        added to the channel-0 block at every frequency (spec 4.2)."""
        blk = np.array(self.vertex.source.get_pair(self.vertex.slot, a, b), copy=True)
        if a == 0 and b == 0 and self.V_inst is not None:
            blk += _vinst_as_block(self.V_inst, self.nd, self.nvol)[None]
        return blk

    def _block_rtau(self, a, b, mod):
        """The ``(a, b)`` block in ``(r, tau)`` on module ``mod``."""
        if not self.is_ir:
            return block_to_rtau_uniform(self._raw_block_uniform(a, b), self.shape,
                                         self.norb, self.workers, mod)
        nd = self.nd
        nx, ny, nz = self.shape
        co = self.vertex.coeffs[:, :, a * nd:(a + 1) * nd, b * nd:(b + 1) * nd]  # (L,nvol,nd,nd)
        co = np.moveaxis(co, 0, -1)                                    # (nvol, nd, nd, L)
        vt = self.vertex.axB.eval_to_tau_points(co, self.axF.tau)      # (nvol, nd, nd, ntau)
        vr = _bk.spatial_ifftn(vt.reshape(nx, ny, nz, nd, nd, -1), axes=(0, 1, 2),
                               workers=self.workers)
        arr = vr.reshape(nx, ny, nz, self.norb, self.norb, self.norb,
                         self.norb, -1).transpose(3, 4, 5, 6, 0, 1, 2, 7)
        return mod.ascontiguousarray(mod.asarray(arr))

    def _const_block_r(self, a, b, mod):
        """The retained frequency-independent IR component of the ``(a, b)``
        block, spatially transformed to ``r`` (spec 4.5)."""
        nd = self.nd
        nx, ny, nz = self.shape
        cb = self.vertex.const[:, a * nd:(a + 1) * nd,
                               b * nd:(b + 1) * nd].reshape(nx, ny, nz, nd, nd)
        cr = _bk.spatial_ifftn(cb, axes=(0, 1, 2), workers=self.workers)
        return mod.asarray(cr.reshape(nx, ny, nz, self.norb, self.norb, self.norb,
                                      self.norb).transpose(3, 4, 5, 6, 0, 1, 2))

    def _get_block(self, a, b):
        if self._blocks is not None:
            return self.xp.asarray(self._blocks[(a, b)])
        return self._block_rtau(a, b, self.xp)

    # -- matvec ---------------------------------------------------------------
    def matvec(self, phi_flat):
        """Apply the kernel to a flat HOST gap; returns a flat host array."""
        xp = self.xp
        ax = (2, 3, 4)
        phi = xp.asarray(np.asarray(phi_flat).reshape(self.gap_shape))
        F = xp.einsum("iljmxyzn,lmxyzn->ijxyzn", self.G2, phi)
        F0_r = None
        if self.is_ir:
            F_coeff = self.axF.fit_from_freq(F)
            F_rt = _bk.spatial_ifftn(self.axF.eval_to_tau(F_coeff), axes=ax,
                                     workers=self.workers)
            if self._Vinst_rt is not None or self._const_rt is not None:
                # the delta(tau) terms' tau integral (the bare vertex and the
                # retained IR constant): the UNREGULARIZED Matsubara sum
                # (1/beta) sum_n F(i w_n), i.e. the MIDPOINT 0.5 (F(0^+) +
                # F(0^-)) of the tau = 0 jump, evaluated exactly through the
                # fermionic basis. NOT u_zero_plus -- see
                # eliashberg_dynamic.eliashberg_kernel_ir.
                u0 = xp.asarray(self.axF.u_matsubara_sum)
                F0_r = _bk.spatial_ifftn(F_coeff @ u0, axes=ax, workers=self.workers)
        else:
            F_rt = _bk.spatial_ifftn(_ms.fermion_to_tau(F, axis=-1), axes=ax,
                                     workers=self.workers)
        acc = [None] * self.B
        acc0 = [None] * self.B
        for b in range(self.B):
            Rb = self.delta_r[b]
            F_b = F_rt if Rb == (0, 0, 0) else xp.roll(F_rt, tuple(-x for x in Rb), axis=ax)
            F0_b = None
            if self._const_rt is not None:
                F0_b = F0_r if Rb == (0, 0, 0) else xp.roll(F0_r, tuple(-x for x in Rb),
                                                            axis=ax)
            for a in range(self.B):
                G_ab = self._get_block(a, b)
                P = -xp.einsum("abcdxyzt,bcxyzt->adxyzt", G_ab, F_b)
                acc[a] = P if acc[a] is None else acc[a] + P
                del G_ab, P
                if F0_b is not None:
                    P0 = -xp.einsum("abcdxyz,bcxyz->adxyz", self._const_rt[(a, b)], F0_b)
                    acc0[a] = P0 if acc0[a] is None else acc0[a] + P0
        out_rt = None
        out0_r = None
        for a in range(self.B):
            Ra = self.delta_r[a]
            t = acc[a] if Ra == (0, 0, 0) else xp.roll(acc[a], tuple(-x for x in Ra), axis=ax)
            out_rt = t if out_rt is None else out_rt + t
            if acc0[a] is not None:
                t0 = acc0[a] if Ra == (0, 0, 0) else xp.roll(acc0[a], tuple(-x for x in Ra),
                                                             axis=ax)
                out0_r = t0 if out0_r is None else out0_r + t0
        if self.is_ir:
            out = self.axF.tau_to_freq(_bk.spatial_fftn(out_rt, axes=ax, workers=self.workers))
            if self._Vinst_rt is not None:
                inst_r = -xp.einsum("abcdxyz,bcxyz->adxyz", self._Vinst_rt, F0_r)
                out = out + _bk.spatial_fftn(inst_r, axes=ax, workers=self.workers)[..., None]
            if out0_r is not None:
                out = out + _bk.spatial_fftn(out0_r, axes=ax, workers=self.workers)[..., None]
            out = self.beta * out
        else:
            out = _ms.tau_to_fermion(_bk.spatial_fftn(out_rt, axes=ax, workers=self.workers),
                                     axis=-1)
        return _bk.to_host(out).ravel()

    def memory_table(self):
        """The admission rows (bytes) plus the residency actually chosen."""
        return dict(self.admission.rows, residency=self.residency)


def gap_bond_projection(gap_w, view, spatial_shape):
    """Project a gap onto the bond form factors: ``psi_m = (1/N) sum_k
    e^{-i k . R_m} Delta(k)`` for every channel ``m`` (spec 4.4, 5.1).

    Parameters
    ----------
    gap_w : ndarray, shape (norb, norb, Nx, Ny, Nz, nfreq)
    view : BondSetView
    spatial_shape : tuple

    Returns
    -------
    ndarray, shape ``(B, norb, norb, nfreq)``
    """
    nx, ny, nz = spatial_shape
    gap = np.asarray(gap_w)
    kx = 2 * np.pi * np.arange(nx) / nx
    ky = 2 * np.pi * np.arange(ny) / ny
    kz = 2 * np.pi * np.arange(nz) / nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    out = np.empty((int(view.n_channels),) + gap.shape[:2] + gap.shape[5:], complex)
    for m, R in enumerate(view.delta_r):
        ph = np.exp(-1j * (KX * R[0] + KY * R[1] + KZ * R[2]))
        out[m] = np.einsum("xyz,abxyzn->abn", ph, gap) / (nx * ny * nz)
    return out
