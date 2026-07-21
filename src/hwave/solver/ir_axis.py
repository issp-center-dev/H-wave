"""Sparse-IR Matsubara axis (design: docs/design/ir-matsubara.md, Sec. 3.1).

Wraps the optional ``sparse-ir`` package behind explicit H-wave conventions:

- Frequencies are identified by INTEGER Matsubara indices n with
  ``i w = n * pi / beta`` (n odd for fermions, even for bosons), matching the
  centered H-wave uniform grid ``n_j = 2 j + 1 - Nmat`` (fermionic) /
  ``n_l = 2 l - Nmat`` (bosonic).
- Node sets are EXACTLY symmetric: integer negation of ``freq_n`` and the
  reflection ``tau -> beta - tau`` of ``tau`` are both plain array reversals
  (asserted at construction; tau nodes are snapped to symmetric pair means).
  Frequency-axis reversals in existing solver code therefore keep their
  meaning on IR node arrays.
- Transforms contract the LAST axis and are PHYSICAL: ``tau_to_freq`` is the
  integral over tau (u_l(tau) fits carry the measure) and ``freq_to_tau``
  contains the 1/beta Matsubara sum through the basis. Consumers therefore
  apply NO extra 1/beta at product exits on the IR path -- unlike the
  uniform-FFT path of ``solver/matsubara.py``, whose bare FFTs need explicit
  1/beta factors there.
- Transform matrices are built once on the host; applications are plain
  matmuls, so they work unchanged on numpy or cupy arrays (the caller moves
  matrices to the device via :func:`hwave.solver.backend.array_module_of`
  dispatch -- see :meth:`IRAxis.mats`).
"""
import logging

import numpy as np

from hwave.solver import backend as _bk

_logger = logging.getLogger(__name__)


def _import_sparse_ir():
    """Import sparse_ir (separated out so tests can monkeypatch a missing
    installation)."""
    import sparse_ir
    return sparse_ir


def is_ir_native(data):
    """Discriminator D-1 (design ir-matsubara-stage3.md Sec. 3): a loaded
    npz holds sparse-IR NODE data iff ``matsubara_basis == "ir"``,
    ``frequency_grid == "sparse_ir_nodes"``, AND ``ir_freq_n`` is present.
    Provenance alone is not enough -- densified outputs (e.g. the Stage-1
    ``gap_dynamic.npz``) carry ``matsubara_basis="ir"`` on a uniform-grid
    array.

    Requires no sparse-ir installation (uniform-only readers use it to
    fail fast even where the IR stack is absent)."""
    basis_ir = ("matsubara_basis" in data
                and str(data["matsubara_basis"]) == "ir")
    grid_sparse = ("frequency_grid" in data
                   and str(data["frequency_grid"]) == "sparse_ir_nodes")
    has_freq_n = "ir_freq_n" in data
    native = basis_ir and grid_sparse and has_freq_n
    # ``matsubara_basis='ir'`` alone is valid provenance on a densified
    # uniform array. Any sparse-node marker, however, commits the file to the
    # complete native discriminator; otherwise a malformed native file could
    # fall through to a permissive legacy uniform-grid reader.
    if (grid_sparse or has_freq_n) and not native:
        raise ValueError(
            "incomplete IR-native discriminator: sparse-node data requires "
            "matsubara_basis='ir', frequency_grid='sparse_ir_nodes', and "
            "ir_freq_n together.")
    return native


def ir_native_meta(data):
    """Extract the IR-node metadata of an IR-native npz (Sec. 3 schema)."""
    freq_n = np.asarray(data["ir_freq_n"])
    if freq_n.ndim != 1:
        raise ValueError(
            "IR-native ir_freq_n must be 1-D, got shape {}."
            .format(freq_n.shape))
    if not np.issubdtype(freq_n.dtype, np.integer):
        raise ValueError(
            "IR-native ir_freq_n must contain integer Matsubara indices, "
            "got dtype {}.".format(freq_n.dtype))
    if np.any(np.diff(freq_n) <= 0):
        raise ValueError(
            "IR-native ir_freq_n must be strictly increasing (sorted, no "
            "duplicates).")
    statistics = str(data["ir_statistics"])
    if statistics not in ("F", "B"):
        raise ValueError(
            "IR-native ir_statistics must be 'F' or 'B', got {!r}."
            .format(statistics))
    want_odd = 1 if statistics == "F" else 0
    if np.any(np.abs(freq_n) % 2 != want_odd):
        raise ValueError(
            "IR-native ir_freq_n parity does not match ir_statistics={!r}."
            .format(statistics))
    return {"freq_n": freq_n.astype(np.int64, copy=False),
            "beta": float(data["ir_beta"]),
            "wmax": float(data["ir_wmax"]),
            "tol": float(data["ir_tol"]),
            "L": int(data["ir_L"]),
            "statistics": statistics}


class IRAxis:
    """One statistics sector (fermionic 'F' or bosonic 'B') of the IR axis.

    Attributes
    ----------
    beta, wmax, eps : float
        Basis parameters (Lambda = beta * wmax).
    L : int
        Number of IR coefficients.
    freq_n : ndarray of int
        Sorted, integer-symmetric Matsubara node indices (i w = n pi / beta).
    n_freq, n_tau : int
        Node counts (>= L).
    tau : ndarray of float
        Sorted tau nodes in (0, beta), exactly reversal-symmetric.
    u_beta_minus : ndarray (L,)
        Basis values u_l(beta^-), for equal-time observables
        (n = -Re[coeffs @ u_beta_minus] for a fermionic Green's function).
    """

    def __init__(self, beta, wmax, eps, statistics):
        try:
            sparse_ir = _import_sparse_ir()
        except ImportError:
            raise ImportError(
                "matsubara_basis='ir' requires the optional 'sparse-ir' "
                "package; install it with 'pip install sparse-ir' "
                "(https://sparse-ir.readthedocs.io).")
        if statistics not in ("F", "B"):
            raise ValueError("statistics must be 'F' or 'B', got {!r}"
                             .format(statistics))
        self.beta = float(beta)
        self.wmax = float(wmax)
        self.eps = float(eps)
        self.statistics = statistics

        basis = sparse_ir.FiniteTempBasis(statistics, self.beta, self.wmax,
                                          eps=self.eps)
        self._basis = basis
        self.L = basis.size

        # --- frequency nodes: integer, symmetrized (union with mirror) ---
        wn = np.asarray(
            sparse_ir.MatsubaraSampling(basis).sampling_points, dtype=np.int64)
        wn = np.unique(np.concatenate([wn, -wn]))
        assert np.array_equal(-wn[::-1], wn), \
            "frequency node set not integer-symmetric after symmetrization"
        self.freq_n = wn
        self.n_freq = wn.size

        # --- tau nodes: symmetrized to exact reversal symmetry ---
        tau = np.sort(np.asarray(
            sparse_ir.TauSampling(basis).sampling_points, dtype=np.float64))
        if not np.allclose(tau, self.beta - tau[::-1],
                           atol=1e-14 * self.beta, rtol=0.0):
            # union with the mirror set, then pair below
            tau = np.sort(np.concatenate([tau, self.beta - tau]))
        # Snap each mirror pair to its symmetric mean. The reflection
        # tau -> beta - tau is IMPLEMENTED as a pure index reversal everywhere
        # (never by recomputing beta - tau), so what must hold is the pairing
        # tau_i + tau_{n-1-i} = beta within tight tolerance -- bitwise value
        # symmetry is not required (and not achievable in floating point).
        tau = 0.5 * (tau + (self.beta - tau[::-1]))
        assert np.allclose(tau + tau[::-1], self.beta,
                           rtol=0.0, atol=1e-13 * self.beta), \
            "tau node pairing tau_i + tau_{n-1-i} = beta violated"
        self.tau = tau
        self.n_tau = tau.size

        # --- transform matrices (host, applied to the LAST axis: arr @ M) ---
        # eval_freq (n_freq, L): node values = eval_freq @ coeffs, so the
        # last-axis application uses its transpose; fits use the pinv
        # transposed likewise. pinv(eval)(L, n) -> fit matrix (n, L).
        eval_freq = basis.uhat(self.freq_n).T          # (n_freq, L)
        eval_tau = basis.u(self.tau).T                 # (n_tau, L)
        self._m = {}
        self._m["fit_freq"] = np.ascontiguousarray(np.linalg.pinv(eval_freq).T)  # (n_freq, L)
        self._m["eval_freq"] = np.ascontiguousarray(eval_freq.T)                 # (L, n_freq)
        self._m["fit_tau"] = np.ascontiguousarray(np.linalg.pinv(eval_tau).T)    # (n_tau, L)
        self._m["eval_tau"] = np.ascontiguousarray(eval_tau.T)                   # (L, n_tau)
        # composites: nodes -> nodes through coefficient space
        self._m["freq_to_tau"] = np.ascontiguousarray(
            self._m["fit_freq"] @ self._m["eval_tau"])         # (n_freq, n_tau)
        self._m["tau_to_freq"] = np.ascontiguousarray(
            self._m["fit_tau"] @ self._m["eval_freq"])         # (n_tau, n_freq)

        self.u_beta_minus = basis.u(self.beta * (1.0 - 1e-15))  # (L,)
        # u_l(0^+): equal-time evaluation at the other endpoint, e.g. the
        # anomalous bubble F(tau=0) = (1/beta) sum_nu F(i nu) needed by the
        # instantaneous-vertex term of the dynamic Eliashberg kernel.
        self.u_zero_plus = basis.u(self.beta * 1e-15)            # (L,)

        # device mirrors, filled lazily per array module
        self._device_m = {}

    # -- matrix access with backend dispatch ---------------------------------

    def mats(self, name, xp):
        """Return transform matrix ``name`` on the array module ``xp``."""
        if xp is np:
            return self._m[name]
        key = name
        if key not in self._device_m:
            self._device_m[key] = xp.asarray(self._m[name])
        return self._device_m[key]

    def _apply(self, arr, name):
        xp = _bk.array_module_of(arr)
        return arr @ self.mats(name, xp)

    # -- transforms (all contract the LAST axis) -----------------------------

    def fit_from_freq(self, arr):
        """(..., n_freq) node values -> (..., L) coefficients."""
        return self._apply(arr, "fit_freq")

    def eval_to_freq(self, coeffs):
        """(..., L) coefficients -> (..., n_freq) node values."""
        return self._apply(coeffs, "eval_freq")

    def fit_from_tau(self, arr):
        """(..., n_tau) node values -> (..., L) coefficients."""
        return self._apply(arr, "fit_tau")

    def eval_to_tau(self, coeffs):
        """(..., L) coefficients -> (..., n_tau) node values."""
        return self._apply(coeffs, "eval_tau")

    def freq_to_tau(self, arr):
        """(..., n_freq) -> (..., n_tau), fused fit+evaluate."""
        return self._apply(arr, "freq_to_tau")

    def tau_to_freq(self, arr):
        """(..., n_tau) -> (..., n_freq), fused fit+evaluate."""
        return self._apply(arr, "tau_to_freq")

    def eval_to_tau_points(self, coeffs, tau_points):
        """Evaluate coefficients on ARBITRARY tau points (e.g. this axis is
        bosonic and ``tau_points`` are the fermionic product nodes). The
        evaluation matrix is cached per tau-point set (and per backend), so
        repeated SCF-loop calls with the same node array are matmuls only."""
        tau_points = np.asarray(tau_points, dtype=np.float64)
        key = ("taupts", hash(tau_points.tobytes()))
        if key not in self._device_m:
            self._device_m[key] = np.ascontiguousarray(
                self._basis.u(tau_points))              # (L, npts)
        m = self._device_m[key]
        xp = _bk.array_module_of(coeffs)
        if xp is not np:
            dkey = key + ("dev",)
            if dkey not in self._device_m:
                self._device_m[dkey] = xp.asarray(m)
            m = self._device_m[dkey]
        return coeffs @ m

    def freq_to_tau_points(self, arr, tau_points):
        """(..., n_freq) node values -> values on ARBITRARY tau points
        (fused fit + evaluate, cached like :meth:`eval_to_tau_points`)."""
        return self.eval_to_tau_points(self.fit_from_freq(arr), tau_points)

    def fit_from_freq_points(self, arr, freq_n):
        """(..., npts) values at ARBITRARY integer Matsubara indices
        ``freq_n`` -> (..., L) coefficients (design: Stage-3 addendum
        Sec. 5; the re-ingestion primitive for IR-native files).

        ``freq_n`` must be 1-D, integer, strictly increasing, of the axis'
        parity, with at least L entries whose design matrix has full rank
        -- everything else raises, because a pinv of a deficient system
        would silently amplify noise into physics."""
        freq_n = np.asarray(freq_n)
        if freq_n.ndim != 1:
            raise ValueError(
                "fit_from_freq_points: freq_n must be 1-D, got shape {}"
                .format(freq_n.shape))
        if not np.issubdtype(freq_n.dtype, np.integer):
            raise ValueError(
                "fit_from_freq_points: freq_n must be integer Matsubara "
                "indices, got dtype {}".format(freq_n.dtype))
        want_odd = 1 if self.statistics == "F" else 0
        if np.any(np.abs(freq_n) % 2 != want_odd):
            raise ValueError(
                "fit_from_freq_points: freq_n parity does not match the "
                "{} axis ({} indices required).".format(
                    "fermionic" if want_odd else "bosonic",
                    "odd" if want_odd else "even"))
        if np.any(np.diff(freq_n) <= 0):
            raise ValueError(
                "fit_from_freq_points: freq_n must be strictly increasing "
                "(sorted, no duplicates).")
        if freq_n.size < self.L:
            raise ValueError(
                "fit_from_freq_points: need at least L={} nodes to "
                "determine the coefficients, got {}."
                .format(self.L, freq_n.size))
        freq_n = freq_n.astype(np.int64)
        self._freq_points_fit_matrix(freq_n)   # build + validate + cache
        return self._apply_cached(arr, ("freqpts", hash(freq_n.tobytes())))

    def eval_to_freq_points(self, coeffs, freq_n):
        """(..., L) coefficients -> values at ARBITRARY integer Matsubara
        indices ``freq_n`` (the evaluation half of
        :meth:`fit_from_freq_points`; used for refit residuals at a file's
        own nodes). The node set must have passed
        :meth:`fit_from_freq_points` validation first (it shares the
        cache)."""
        freq_n = np.asarray(freq_n, dtype=np.int64)
        key = ("freqpts_eval", hash(freq_n.tobytes()))
        if key not in self._device_m:
            self._device_m[key] = np.ascontiguousarray(
                self._basis.uhat(freq_n))               # (L, npts)
        return self._apply_cached(coeffs, key)

    def _apply_cached(self, arr, key):
        m = self._device_m[key]
        xp = _bk.array_module_of(arr)
        if xp is not np:
            dkey = key + ("dev",)
            if dkey not in self._device_m:
                self._device_m[dkey] = xp.asarray(m)
            m = self._device_m[dkey]
        return arr @ m

    def _freq_points_fit_matrix(self, freq_n):
        key = ("freqpts", hash(freq_n.tobytes()))
        if key not in self._device_m:
            design = np.ascontiguousarray(
                self._basis.uhat(freq_n).T)             # (npts, L)
            rank = np.linalg.matrix_rank(design)
            if rank < self.L:
                raise ValueError(
                    "fit_from_freq_points: the design matrix for this node "
                    "set is rank-deficient ({} < L={}); the node set cannot "
                    "determine the coefficients.".format(rank, self.L))
            sv = np.linalg.svd(design, compute_uv=False)
            cond = float(sv[0] / sv[-1])
            if cond > 1e10:
                _logger.warning(
                    "fit_from_freq_points: ill-conditioned node set "
                    "(cond=%.2e); the fit may amplify noise.", cond)
            else:
                _logger.debug("fit_from_freq_points: cond=%.2e", cond)
            self._device_m[key] = np.ascontiguousarray(
                np.linalg.pinv(design).T)               # (npts, L)

    # -- uniform-grid (centered H-wave) interface ----------------------------

    def _uniform_n(self, nmat):
        if self.statistics == "F":
            return 2 * np.arange(nmat, dtype=np.int64) + 1 - nmat
        return 2 * np.arange(nmat, dtype=np.int64) - nmat

    def uniform_matrices(self, nmat, with_constant=False):
        """(fit, eval) matrices for the centered uniform H-wave grid:
        fit: (nmat, L) applied as arr(..., nmat) @ fit -> coeffs;
        eval: (L, nmat) applied as coeffs(..., L) @ eval -> uniform values.

        With ``with_constant=True`` the fit design matrix gains a constant
        column (fit shape (nmat, L+1); the last output component is the
        fitted constant). A frequency-independent constant is the known
        O(beta/Nmat) discretization artifact of the uniform-FFT
        susceptibilities (a delta(tau) component, which the IR basis
        correctly cannot represent); the augmented fit isolates it so it can
        be discarded."""
        key = ("uniform", int(nmat), bool(with_constant))
        if key not in self._device_m:
            eval_u = self._basis.uhat(self._uniform_n(nmat)).T   # (nmat, L)
            ev = np.ascontiguousarray(eval_u.T)                   # (L, nmat)
            design = eval_u
            if with_constant:
                design = np.hstack([eval_u,
                                    np.ones((nmat, 1), dtype=eval_u.dtype)])
            fit = np.ascontiguousarray(np.linalg.pinv(design).T)  # (nmat, L[+1])
            self._device_m[key] = (fit, ev)
        return self._device_m[key]

    def fit_from_uniform(self, arr, nmat):
        """(..., nmat) centered-uniform values -> (..., L) coefficients
        (least squares over ALL uniform frequencies)."""
        fit, _ = self.uniform_matrices(nmat)
        xp = _bk.array_module_of(arr)
        if xp is not np:
            fit = xp.asarray(fit)
        return arr @ fit

    def eval_to_uniform(self, coeffs, nmat):
        """(..., L) coefficients -> (..., nmat) centered-uniform values."""
        _, ev = self.uniform_matrices(nmat)
        xp = _bk.array_module_of(coeffs)
        if xp is not np:
            ev = xp.asarray(ev)
        return coeffs @ ev
