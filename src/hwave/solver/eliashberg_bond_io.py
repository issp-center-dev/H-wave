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
    mod = {tuple(int(x) for x in np.mod(r, cell)) for r in delta_r}
    if len(mod) != B:
        raise ValueError("bond archive {}: two channels of delta_r coincide modulo cell_shape {}".format(path, cell))
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


class _ArchiveView:
    """The ``BondSetView`` contract (``delta_r``, ``reverse``, ``n_channels``)
    rebuilt from the archive arrays: the pairing kernel needs only the bond
    topology, and the archive carries it, so post-processing never re-derives
    the shell enumeration from the interaction files (where a different
    ``bond_max_shells`` would silently give a different channel set)."""

    def __init__(self, delta_r, reverse):
        self.delta_r = tuple(tuple(int(x) for x in r) for r in delta_r)
        self.reverse = tuple(int(x) for x in reverse)
        self.n_channels = len(self.delta_r)


def _frequency_batch(nmat, nvol, ND, host_cap):
    """Frequency batch size of the vertex build: one batch (the seven
    ``nb * nvol * ND**2`` complex workspaces of the spec-8 dressing row) stays
    under a quarter of the host cap. At least 1, at most ``nmat``."""
    per_l = 7 * nvol * ND * ND * 16
    return max(1, min(int(nmat), int(0.25 * host_cap // max(1, per_l))))


def solve_dynamic_bond(input_dict):
    """Bond-resolved dynamic Eliashberg solve from a FLEX bond archive (spec 6).

    The post-processing entry of ``hwave_sc``: reads the enlarged bond
    vertices and the bond-resolved susceptibilities written by a FLEX run with
    ``[mode.param] longitudinal_bond_channels = true`` and
    ``longitudinal_bond_output_full = true``, builds the pairing vertex from
    them (the two archive members consumed ONE AT A TIME), and solves the
    linearized Eliashberg eigenproblem for ONE pairing channel.

    Returns
    -------
    float
        The leading eigenvalue lambda.
    """
    import hwave.sc as sc
    from . import backend as _bk
    from . import eliashberg_bond as _eb
    from . import eliashberg_dynamic as _ed
    from . import flex_bond as _fb

    # publicly callable, so guard the unsupported mode here too (issue #83)
    sc.reject_spin_orbital_mode(input_dict)
    mode_param = input_dict["mode"]["param"]
    T = mode_param["T"]
    beta = sc._coerce_run_beta(T)
    cell = list(mode_param["CellShape"])
    cell = cell + [1] * (3 - len(cell))
    sub = list(mode_param.get("SubShape", cell))
    sub = sub + [1] * (3 - len(sub))
    if sub != [1, 1, 1]:
        raise ValueError(
            "SubShape (sublattice folding) is not supported by the Eliashberg "
            "module: fold the model into the unit cell yourself, or set "
            "SubShape = [1, 1, 1] explicitly.")
    Nx, Ny, Nz = (int(x) for x in cell)
    nvol = Nx * Ny * Nz
    nmat = int(mode_param.get("Nmat", 1024))

    eli_param = input_dict.get("eliashberg", {})
    pairing_type = str(eli_param.get("pairing_type", "singlet")).lower()
    # ONE channel per hwave_sc run: the outputs keep the on-site file names,
    # so a second channel is a second run with its own path_to_output.
    ctl = PairingControls.from_param(eli_param, pairing_types=(pairing_type,))
    use_ir = ctl.matsubara_basis == "ir"
    xp, _gpu_active = _bk.get_backend(_ed._gpu_requested(eli_param), logger=logger,
                                      required=_ed._gpu_required_requested(eli_param))

    geom_info, hr, interactions = sc._read_interaction_files(input_dict)
    norb = int(geom_info["norb"])
    nd = norb * norb

    # --- the archive (validated before any large member is touched) --------
    flex_dir = sc._resolve_flex_dir(input_dict)
    arch_name = str(eli_param.get("flex_bond_archive", "longitudinal_bond.npz"))
    arch_path = arch_name if os.path.isabs(arch_name) else os.path.join(flex_dir, arch_name)
    arch = load_bond_archive(arch_path, norb=norb, nmat_expected=nmat,
                             cell_shape_expected=(Nx, Ny, Nz), beta_expected=beta)
    logger.info("bond pairing: archive %s (B = %d, ND = %d, nmat = %d)",
                arch_path, arch.B, arch.ND, arch.nmat)

    green = sc._load_flex_green(input_dict, norb, Nx, Ny, Nz)
    if green is None:
        raise ValueError(
            "dynamic bond Eliashberg requires the dressed green.npz of the same "
            "FLEX run (the pair bubble G2 is built from it); none was found in "
            "{}. Check [file.input] path_to_flex_output / [eliashberg] flex_green."
            .format(flex_dir))
    if green.shape[-1] != nmat:
        raise ValueError("green.npz has nmat = {} but this run uses Nmat = {}"
                         .format(green.shape[-1], nmat))

    view = _ArchiveView(arch.delta_r, arch.reverse)
    kx, ky, kz = (np.linspace(0.0, 2.0 * np.pi, n, endpoint=False) for n in (Nx, Ny, Nz))

    # --- IR axes (optional) -------------------------------------------------
    axF = axB = None
    if use_ir:
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)
        axF, axB = _ed._ir_axes_for_run({"ir_tol": ctl.ir_tol, "ir_wmax": ctl.ir_wmax},
                                        beta, hr, inter_k, norb, mu=mode_param.get("mu"),
                                        filling=mode_param.get("filling"))
    nfreq = axF.n_freq if use_ir else nmat
    L_B = axB.L if use_ir else 0

    # --- admission (spec 8) -------------------------------------------------
    host_cap = (ctl.bond_memory_cap_gb * _eb._GIB) if ctl.bond_memory_cap_gb \
        else 0.8 * _eb._host_available_bytes()
    device_cap = host_cap if xp is np else 0.9 * _bk.device_available_bytes()
    nb = _frequency_batch(nmat, nvol, arch.ND, host_cap)
    table = _eb.estimate_pair_memory(
        nmat=nmat, ntau=(axF.n_tau if use_ir else 0), nfreq=nfreq, nvol=nvol, norb=norb,
        B=arch.B, num_eigenvalues=ctl.num_eigenvalues, residency="auto", ir=use_ir,
        L_B=L_B, n_channels=1, in_process=False, nb=nb)
    logger.info("bond pairing admission (post-processing, GiB): %s",
                {k: round(v / _eb._GIB, 3) for k, v in table.rows.items()})
    # raises MemoryError with the full table when nothing fits
    residency = table.choose("host" if xp is np else "auto", host_cap, device_cap)

    with _fb.BondDeviceContext(xp, arch.S_bond, arch.C_bond) as dev:
        acc = _eb.PairVertexAccumulator(
            dev, pairing_types=(pairing_type,), nb=nb, nmat=nmat, nvol=nvol, nd=nd,
            spatial_shape=(Nx, Ny, Nz), ir=((axF, axB) if use_ir else None),
            ir_keep_static=ctl.ir_keep_static_chi, workers=ctl.fft_workers)

        def stages(a):
            # one archive member resident at a time (spec 6): the IR residual
            # pass replays this and re-reads each member once more
            for channel, name in (("spin", "chi_s_w"), ("charge", "chi_c_w")):
                arr = arch.member(name)
                a.add_channel(channel, _eb.ArrayBlockSource({name: arr}, nd), name)
                del arr

        stages(acc)
        vertex = acc.finish(ir_fit_tol=ctl.ir_fit_tol, stage_callable=stages)[pairing_type]
        del acc
        V_inst = _eb.instantaneous_vertex(arch.S_bond, arch.C_bond, nd, pairing_type,
                                          (Nx, Ny, Nz))
        g = _ed._ir_compress(green, axF, nmat, "green") if use_ir else green
        G2 = xp.asarray(_ed.calc_g2_dynamic(g, beta))
        del g
        kernel = _eb.BondPairKernel(
            vertex, G2, view, xp=xp, spatial_shape=(Nx, Ny, Nz), norb=norb, beta=beta,
            nfreq=nfreq, V_inst=V_inst, axF=axF, residency=residency, admission=table,
            host_cap=host_cap, device_cap=device_cap, workers=ctl.fft_workers)
        logger.info("bond pairing kernel residency: %s", kernel.residency)
        phi0, seed_vec = _ed.build_seed(ctl.as_eli_param(), pairing_type, norb, kx, ky, kz,
                                        kernel.gap_shape, use_ir, axF, nmat)
        lam, gap_w, eigenvalues_all, eigenvalue_match, note = _ed.run_leading_eigenproblem(
            kernel.matvec, kernel.gap_shape, ctl.as_eli_param(), pairing_type, phi0=phi0,
            seed_vec=seed_vec, use_ir=use_ir, axF=axF, nmat=nmat,
            logger_label="bond pairing kernel", parity_leakage_policy="refuse")
        residency_used = kernel.residency

    # --- outputs (the on-site dynamic file set plus the bond provenance) ----
    out_dir = input_dict["file"]["output"]["path_to_output"]
    os.makedirs(out_dir, exist_ok=True)
    _ed.write_eigenvalue_file(
        os.path.join(out_dir, eli_param.get("output_eigenvalue", "eigenvalue.dat")),
        lam, eigenvalues_all, eigenvalue_match, note,
        header_lines=["bond_channels=true", "residency=" + residency_used])
    extra = {"bond_channels": True, "bond_delta_r": arch.delta_r, "bond_reverse": arch.reverse,
             "bond_archive": arch_path, "bond_residency": residency_used,
             "gap_bond_projection": _eb.gap_bond_projection(gap_w, view, (Nx, Ny, Nz))}
    if use_ir:
        extra.update({"matsubara_basis": "ir", "ir_tol": axF.eps, "ir_wmax": axF.wmax,
                      "ir_L": axF.L, "bond_ir_fit_residual_rel": vertex.fit_residual_rel})
    _ed.write_dynamic_outputs(out_dir, gap_w, lam, T, pairing_type, kx, ky, kz, beta,
                              gap_file=eli_param.get("output_gap", "gap.dat"),
                              extra_meta=extra)
    return lam
