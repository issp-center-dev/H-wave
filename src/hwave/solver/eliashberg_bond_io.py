"""Entries, configuration and files of the bond-resolved dynamic pairing kernel."""
import logging
import os
import zipfile
from dataclasses import dataclass

import numpy as np
from numpy.lib import format as _npfmt

logger = logging.getLogger("qlms").getChild("solver").getChild("eliashberg_bond")

_INDEX_ORDER = "I = m*norb**2 + l1*norb + l2"
_FREQ_AXIS = "bosonic l -> 2l - nmat"
_SOLVER_MODES = ("iteration", "eigenvalue", "both")


@dataclass(frozen=True)
class PairingControls:
    pairing_types: tuple
    solver_mode: str = "iteration"
    eigenvalue_method: str = "arnoldi"
    num_eigenvalues: int = 10
    max_iter: int = 1000
    alpha: float = 0.5
    convergence_tol: float = 1.0e-5
    spectral_shift: object = None
    sigma_shift: object = None
    init_gap: object = None
    seed_eigenvector: object = None
    matsubara_basis: str = "uniform"
    ir_wmax: object = None
    ir_tol: float = 1.0e-8
    ir_keep_static_chi: bool = False
    ir_fit_tol: float = 0.1
    fft_workers: int = 1
    bond_memory_cap_gb: object = None

    @classmethod
    def from_param(cls, eli_param, *, pairing_types):
        from . import backend as _bk
        p = dict(eli_param or {})
        for eta in pairing_types:
            if eta not in ("singlet", "triplet"):
                raise ValueError("pairing type must be 'singlet' or 'triplet', got {!r}".format(eta))

        def _pos_int(key, default):
            v = p.get(key, default)
            if isinstance(v, bool) or int(v) != v or int(v) <= 0:
                raise ValueError("[eliashberg] {} must be a positive integer, got {!r}".format(key, v))
            return int(v)

        def _pos_float(key, default, allow_zero=False, allow_none=False):
            v = p.get(key, default)
            if v is None and allow_none:
                return None
            v = float(v)
            if not np.isfinite(v) or v < 0 or (v == 0 and not allow_zero):
                raise ValueError("[eliashberg] {} must be a positive finite number, got {!r}".format(key, p.get(key)))
            return v

        sm = str(p.get("solver_mode", "iteration")).lower()
        if sm not in _SOLVER_MODES:
            raise ValueError("[eliashberg] solver_mode must be one of {}, got {!r}".format(_SOLVER_MODES, sm))
        mb = str(p.get("matsubara_basis", "uniform")).lower()
        if mb not in ("uniform", "ir"):
            raise ValueError("[eliashberg] matsubara_basis must be 'uniform' or 'ir', got {!r}".format(mb))
        return cls(
            pairing_types=tuple(pairing_types), solver_mode=sm,
            eigenvalue_method=str(p.get("eigenvalue_method", "arnoldi")),
            num_eigenvalues=_pos_int("num_eigenvalues", 10), max_iter=_pos_int("max_iter", 1000),
            alpha=_pos_float("alpha", 0.5), convergence_tol=_pos_float("convergence_tol", 1.0e-5),
            spectral_shift=p.get("spectral_shift"), sigma_shift=p.get("sigma_shift"),
            init_gap=p.get("init_gap"), seed_eigenvector=p.get("seed_eigenvector"),
            matsubara_basis=mb, ir_wmax=p.get("ir_wmax"), ir_tol=_pos_float("ir_tol", 1.0e-8),
            ir_keep_static_chi=_bk.as_bool(p.get("ir_keep_static_chi", False)),
            ir_fit_tol=_pos_float("ir_fit_tol", 0.1, allow_zero=True),
            fft_workers=int(p.get("fft_workers", 1)),
            bond_memory_cap_gb=_pos_float("bond_memory_cap_gb", None, allow_none=True))

    def as_eli_param(self):
        """The dict run_leading_eigenproblem / build_seed read."""
        return {"solver_mode": self.solver_mode, "eigenvalue_method": self.eigenvalue_method,
                "num_eigenvalues": self.num_eigenvalues, "max_iter": self.max_iter, "alpha": self.alpha,
                "convergence_tol": self.convergence_tol, "spectral_shift": self.spectral_shift,
                "sigma_shift": self.sigma_shift, "init_gap": self.init_gap,
                "seed_eigenvector": self.seed_eigenvector}


def _npz_member_header(path, name):
    """(shape, dtype) of one npz member from its npy header, without loading it."""
    with zipfile.ZipFile(path) as zf:
        with zf.open(name + ".npy") as fh:
            version = _npfmt.read_magic(fh)
            if version == (1, 0):
                shape, fortran, dtype = _npfmt.read_array_header_1_0(fh)
            elif version == (2, 0):
                shape, fortran, dtype = _npfmt.read_array_header_2_0(fh)
            else:
                raise ValueError(
                    "npz member {!r} in {}: unsupported .npy format version {}".format(name, path, version))
    return tuple(int(x) for x in shape), dtype


class BondArchive:
    def __init__(self, path, small, chi_shape):
        self.path = path
        for k, v in small.items():
            setattr(self, k, v)
        self.chi_shape = chi_shape

    def member(self, name):
        with np.load(self.path) as d:
            arr = np.asarray(d[name])
        if not np.all(np.isfinite(arr)):
            raise ValueError("bond archive {}: member {!r} contains non-finite values".format(self.path, name))
        return arr


def load_bond_archive(path, *, norb, nmat_expected, cell_shape_expected, beta_expected):
    from hwave.solver.rpa import check_momentum_marker
    if not os.path.exists(path):
        raise ValueError("bond archive not found: {} (post-processing needs a FLEX run with "
                         "longitudinal_bond_channels = true and longitudinal_bond_output_full = true)".format(path))
    with np.load(path) as d:
        files = set(d.files)
        schema = int(d["bond_archive_schema"]) if "bond_archive_schema" in files else 0
        if schema == 1 or "S_bond" not in files:
            raise ValueError(
                "bond archive {}: schema {} was written by a version without the bond vertices S and C; "
                "re-run FLEX under this version with longitudinal_bond_channels = true and "
                "longitudinal_bond_output_full = true".format(path, schema))
        if schema != 2:
            raise ValueError("bond archive {}: unsupported archive schema {}".format(path, schema))
        if str(d["index_order"]) != _INDEX_ORDER:
            raise ValueError("bond archive {}: index_order {!r} != {!r}".format(path, str(d["index_order"]), _INDEX_ORDER))
        if str(d["freq_axis"]) != _FREQ_AXIS:
            raise ValueError("bond archive {}: freq_axis {!r} != {!r}".format(path, str(d["freq_axis"]), _FREQ_AXIS))
        if int(d["norb"]) != int(norb):
            raise ValueError("bond archive {}: norb {} != the geometry's {}".format(path, int(d["norb"]), norb))
        nmat = int(d["nmat"])
        if nmat != int(nmat_expected):
            raise ValueError("bond archive {}: nmat {} != this run's Nmat {}".format(path, nmat, nmat_expected))
        cell = tuple(int(x) for x in d["cell_shape"])
        if cell != tuple(int(x) for x in cell_shape_expected):
            raise ValueError("bond archive {}: cell_shape {} != this run's {}".format(path, cell, tuple(cell_shape_expected)))
        beta = float(d["beta"])
        if abs(beta - beta_expected) > 1e-10 * max(1.0, abs(beta_expected)):
            raise ValueError("bond archive {}: beta {} != this run's {}".format(path, beta, beta_expected))
        check_momentum_marker(d, path)
        delta_r = np.asarray(d["delta_r"], dtype=np.int64)
        reverse = np.asarray(d["reverse"], dtype=np.int64)
        types = tuple(str(t) for t in np.asarray(d["types"]).ravel())
        S_bond = np.asarray(d["S_bond"]); C_bond = np.asarray(d["C_bond"])
    B = int(delta_r.shape[0])
    nd = norb * norb
    ND = B * nd
    nvol = int(np.prod(cell))
    if delta_r.ndim != 2 or delta_r.shape[1] != 3 or tuple(delta_r[0]) != (0, 0, 0):
        raise ValueError("bond archive {}: delta_r must be (B, 3) with delta_r[0] = (0, 0, 0)".format(path))
    # Each row must already be its own minimum-image representative modulo
    # cell_shape (component r reduced into the canonical [-m/2, m/2) range
    # per axis of modulus m) -- a raw value differing from that reduction
    # (e.g. R and R + k*L on some axis) is the SAME physical bond stored
    # non-canonically, i.e. it aliases another representative modulo the
    # cell; and two rows may not share one canonical representative either.
    half = np.asarray(cell, dtype=np.int64) // 2
    reduced = np.mod(delta_r, cell)
    reduced = np.where(reduced > half, reduced - np.asarray(cell, dtype=np.int64), reduced)
    if not np.array_equal(reduced, delta_r) \
            or len({tuple(int(x) for x in r) for r in reduced}) != B:
        raise ValueError(
            "bond archive {}: delta_r rows must be pairwise distinct modulo cell_shape {} and "
            "already stored in canonical (minimum-image) form".format(path, cell))
    if reverse.shape != (B,) or reverse[0] != 0 or np.any(reverse < 0) or np.any(reverse >= B) \
            or np.any(reverse[reverse] != np.arange(B)) or np.any(delta_r[reverse] != -delta_r):
        raise ValueError("bond archive {}: reverse is not the reversal involution of delta_r".format(path))
    for name, arr in (("S_bond", S_bond), ("C_bond", C_bond)):
        if arr.shape != (nvol, ND, ND) or not np.all(np.isfinite(arr)):
            raise ValueError("bond archive {}: {} must be finite with shape {}, got {}".format(path, name, (nvol, ND, ND), arr.shape))
    for name in ("chi_s_w", "chi_c_w"):
        shape, dtype = _npz_member_header(path, name)
        if shape != (nmat, nvol, ND, ND) or np.dtype(dtype) != np.dtype(np.complex128):
            raise ValueError("bond archive {}: {} must be complex128 of shape {}, got {} {}".format(path, name, (nmat, nvol, ND, ND), shape, dtype))
    small = dict(S_bond=S_bond, C_bond=C_bond, delta_r=delta_r, reverse=reverse, types=types, B=B, ND=ND,
                 nmat=nmat, nvol=nvol, cell_shape=cell, beta=beta, norb=int(norb))
    return BondArchive(path, small, (nmat, nvol, ND, ND))
