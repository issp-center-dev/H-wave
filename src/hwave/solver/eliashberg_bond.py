"""Bond-resolved, frequency-resolved Eliashberg pairing kernel built from the
FLEX bond gate's own vertices and dressed susceptibilities (design spec
2026-09-17-eliashberg-bond-dynamic-design.md). Array-only: no TOML, no files."""
import logging
from dataclasses import dataclass

import numpy as np

from . import backend as _bk
from . import bond_channels as _bc
from . import flex_bond as _fb

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
    def finish(self, *, ir_fit_tol=1.0e-4, stage_callable=None):
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
            try:
                stage_callable(self)
            finally:
                self._resid_mode = False
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
                         L_B, n_channels, in_process):
    """The spec-8 incremental memory table of the pairing step. ``L_B`` is the
    bosonic coefficient count (``0`` on the uniform grid) and ``n_channels``
    (1 or 2) the number of pairing types whose IR coefficient rows are alive."""
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
        "G2": norb ** 4 * nvol * nfreq * 16,
        "gap_work": (3 + B) * vec,
        "hoisted_blocks": B * B * blk,
        "stream_workspace": 3 * blk,
        "eigen_vectors": (ncv + 2) * vec,
    }
    build_host = max(rows["vertex_slot"], rows["coefficients_build"])
    dev_common = rows["G2"] + rows["gap_work"] + rows["stream_workspace"]
    need = {
        "device": (int(1.25 * (build_host + rows["eigen_vectors"])),
                   int(1.25 * (dev_common + rows["hoisted_blocks"]))),
        "host": (int(1.25 * (build_host + rows["hoisted_blocks"] + rows["eigen_vectors"])),
                 int(1.25 * dev_common)),
        "stream": (int(1.25 * (build_host + rows["eigen_vectors"])), int(1.25 * dev_common)),
    }
    return AdmissionTable(rows=rows, residency_need=need, requested=residency)
