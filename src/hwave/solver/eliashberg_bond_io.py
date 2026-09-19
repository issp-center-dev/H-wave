"""Entries, configuration and files of the bond-resolved dynamic pairing kernel."""
import logging
import os
import secrets
import zipfile
from dataclasses import dataclass

import numpy as np
from numpy.lib import format as _npfmt

logger = logging.getLogger("qlms").getChild("solver").getChild("eliashberg_bond")

_INDEX_ORDER = "I = m*norb**2 + l1*norb + l2"
_FREQ_AXIS = "bosonic l -> 2l - nmat"
_SOLVER_MODES = ("iteration", "eigenvalue", "both")
#: the eigen family sc._solve_leading accepts (see sc.calc_eliashberg)
_EIGENVALUE_METHODS = ("arnoldi", "shift-invert-bicgstab", "shift-invert-gmres",
                       "shift-invert-lgmres")


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
    ir_fit_tol: float = 0.5
    parity_leakage_tol: object = None
    fft_workers: int = 1
    bond_memory_cap_gb: object = None

    @classmethod
    def from_param(cls, eli_param, *, pairing_types):
        p = dict(eli_param or {})
        for eta in pairing_types:
            if eta not in ("singlet", "triplet"):
                raise ValueError("pairing type must be 'singlet' or 'triplet', got {!r}".format(eta))

        def _fft_workers(key, default):
            """``scipy.fft`` worker counts: a positive count, or a negative one
            that ``scipy.fft`` reads as ``os.cpu_count() + 1 + n`` (``-1`` =
            every core, the value the dynamic-solver documentation names) and
            accepts only down to ``-os.cpu_count()``; ``0`` is refused by
            ``scipy.fft`` itself. Checked here so the run fails at input
            parsing, not at the first spatial FFT after the FLEX solve."""
            v = p.get(key, default)
            floor = -(os.cpu_count() or 1)
            if isinstance(v, bool) or int(v) != v or int(v) == 0 or int(v) < floor:
                raise ValueError("[eliashberg] {} must be a positive integer or a negative one "
                                 "down to {} (-1 = use all cores), got {!r}"
                                 .format(key, floor, v))
            return int(v)

        def _pos_int(key, default):
            v = p.get(key, default)
            if isinstance(v, bool) or int(v) != v or int(v) <= 0:
                raise ValueError("[eliashberg] {} must be a positive integer, got {!r}".format(key, v))
            return int(v)

        def _pos_float(key, default, allow_zero=False, allow_none=False):
            v = p.get(key, default)
            if v is None and allow_none:
                return None
            # a bool is an int in Python: "alpha = true" would become 1.0
            if isinstance(v, bool):
                raise ValueError("[eliashberg] {} must be a positive finite number, got {!r}".format(key, p.get(key)))
            v = float(v)
            if not np.isfinite(v) or v < 0 or (v == 0 and not allow_zero):
                raise ValueError("[eliashberg] {} must be a positive finite number, got {!r}".format(key, p.get(key)))
            return v

        def _bool(key, default):
            """A boolean option, refusing unrecognised spellings: read as
            plain truthiness a typo like "ture" would become TRUE and the
            silently wrong branch would run."""
            v = p.get(key, default)
            if isinstance(v, bool):
                return v
            if isinstance(v, str):
                s = v.strip().lower()
                if s in ("true", "yes", "on", "1"):
                    return True
                if s in ("false", "no", "off", "0"):
                    return False
            raise ValueError("[eliashberg] {} must be a boolean, got {!r}".format(key, v))

        sm = str(p.get("solver_mode", "iteration")).lower()
        if sm not in _SOLVER_MODES:
            raise ValueError("[eliashberg] solver_mode must be one of {}, got {!r}".format(_SOLVER_MODES, sm))
        mb = str(p.get("matsubara_basis", "uniform")).lower()
        if mb not in ("uniform", "ir"):
            raise ValueError("[eliashberg] matsubara_basis must be 'uniform' or 'ir', got {!r}".format(mb))
        em = str(p.get("eigenvalue_method", "arnoldi")).lower()
        if em not in _EIGENVALUE_METHODS:
            raise ValueError("[eliashberg] eigenvalue_method must be one of {}, got {!r}"
                             .format(_EIGENVALUE_METHODS, p.get("eigenvalue_method")))
        return cls(
            pairing_types=tuple(pairing_types), solver_mode=sm,
            eigenvalue_method=em,
            num_eigenvalues=_pos_int("num_eigenvalues", 10), max_iter=_pos_int("max_iter", 1000),
            alpha=_pos_float("alpha", 0.5), convergence_tol=_pos_float("convergence_tol", 1.0e-5),
            spectral_shift=p.get("spectral_shift"), sigma_shift=p.get("sigma_shift"),
            init_gap=p.get("init_gap"), seed_eigenvector=p.get("seed_eigenvector"),
            matsubara_basis=mb, ir_wmax=p.get("ir_wmax"), ir_tol=_pos_float("ir_tol", 1.0e-8),
            ir_keep_static_chi=_bool("ir_keep_static_chi", False),
            ir_fit_tol=_pos_float("ir_fit_tol", 0.5, allow_zero=True),
            parity_leakage_tol=_pos_float("parity_leakage_tol", None, allow_zero=True,
                                          allow_none=True),
            fft_workers=_fft_workers("fft_workers", 1),
            bond_memory_cap_gb=_pos_float("bond_memory_cap_gb", None, allow_none=True))

    @property
    def resolved_parity_leakage_tol(self):
        """The parity-probe refusal threshold this run actually uses (4.4).

        ``parity_leakage_tol = None`` (the default) resolves per basis: 1e-8
        on the uniform grid, which real archives reach at machine precision,
        and 2e-2 with ``matsubara_basis = "ir"``, where the IR representation
        of a uniform-FFT archive carries a parity asymmetry of its own that
        decays as ``Nmat^-2`` (measured 7.6e-3 / 1.3e-3 at Nmat 64 / 128 on
        the single-band U = 4, V = 1 run) and is not kernel algebra.
        """
        if self.parity_leakage_tol is not None:
            return float(self.parity_leakage_tol)
        return 2.0e-2 if self.matsubara_basis == "ir" else 1.0e-8

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


def _as_index_array(member, name, path):
    """One integer-valued archive member as ``int64``.

    ``astype(np.int64)`` TRUNCATES, so a float member carrying 0.5 would
    become a DIFFERENT bond topology without a word; the values must be
    exactly integral (and finite) before the cast."""
    arr = np.asarray(member)
    if not np.issubdtype(arr.dtype, np.integer):
        if not np.issubdtype(arr.dtype, np.floating) or not np.all(np.isfinite(arr)) \
                or not np.all(np.mod(arr, 1) == 0):
            raise ValueError("bond archive {}: {} must hold integers, got dtype {} with "
                             "non-integral values".format(path, name, arr.dtype))
    return arr.astype(np.int64)


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
        # every comparison with nan is false, so an archive whose beta is nan
        # would pass the relative check below and the run would proceed on an
        # unknown temperature
        if not np.isfinite(beta) or beta <= 0:
            raise ValueError("bond archive {}: beta must be a finite positive number, got {}"
                             .format(path, beta))
        if not np.isfinite(beta_expected) or beta_expected <= 0:
            raise ValueError("bond archive {}: this run's beta must be a finite positive "
                             "number, got {}".format(path, beta_expected))
        if abs(beta - beta_expected) > 1e-10 * max(1.0, abs(beta_expected)):
            raise ValueError("bond archive {}: beta {} != this run's {}".format(path, beta, beta_expected))
        check_momentum_marker(d, path)
        delta_r = _as_index_array(d["delta_r"], "delta_r", path)
        reverse = _as_index_array(d["reverse"], "reverse", path)
        types = tuple(str(t) for t in np.asarray(d["types"]).ravel())
        S_bond = np.asarray(d["S_bond"]); C_bond = np.asarray(d["C_bond"])
    # the rank check comes BEFORE shape[0] is read as the channel count
    if delta_r.ndim != 2 or delta_r.shape[1] != 3 or delta_r.shape[0] < 1:
        raise ValueError("bond archive {}: delta_r must be (B, 3) with delta_r[0] = (0, 0, 0)".format(path))
    B = int(delta_r.shape[0])
    nd = norb * norb
    ND = B * nd
    nvol = int(np.prod(cell))
    if tuple(delta_r[0]) != (0, 0, 0):
        raise ValueError("bond archive {}: delta_r must be (B, 3) with delta_r[0] = (0, 0, 0)".format(path))
    # NOTE: ``types`` holds the INTERACTION TYPE NAMES the bond vertices were
    # built from (CoulombInter / Hund / Ising), not one entry per channel, so
    # its length is unrelated to B and nothing is checked against B here.
    if not types:
        raise ValueError("bond archive {}: types is empty".format(path))
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


def _device_cap(host_cap, where):
    """The device byte cap of an entry-level admission, with the fallback of
    :func:`~hwave.solver.eliashberg_bond.device_cap_or_host` and this module's
    logger: the device probe is allowed to have no answer (it returns ``None``
    when the query fails, not only when there is no device), and an admission
    that multiplied that would lose the run to a ``TypeError``."""
    from . import eliashberg_bond as _eb
    return _eb.device_cap_or_host(host_cap, where, log=logger)


def _frequency_batch(nmat, nvol, ND, host_cap, device_cap=None):
    """Frequency batch size of the vertex build: at least 1, at most ``nmat``.

    One batch is the seven ``nb * nvol * ND**2`` complex workspaces of the
    spec-8 dressing row, and that batch lives on BOTH sides -- the archive
    member is sliced on the host and the slice is transferred to the array
    module (``PairVertexAccumulator.add_channel`` -> ``to_device``), while
    ``estimate_pair_memory`` charges ``dressing_workspace`` to every
    residency's DEVICE need. Sizing ``nb`` from the host cap alone therefore
    oversizes it on a device backend with plenty of host memory, and admission
    then refuses a configuration a smaller ``nb`` would have passed. So the
    batch is the smaller of the two quarter-cap limits, as the FLEX gate does
    (:mod:`~hwave.solver.flex_bond`). ``device_cap=None`` means one address
    space (numpy), where the host limit already covers both.
    """
    per_l = 7 * nvol * ND * ND * 16
    nb = int(0.25 * host_cap // max(1, per_l))
    if device_cap is not None:
        nb = min(nb, int(0.25 * device_cap // max(1, per_l)))
    return max(1, min(int(nmat), nb))


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
    device_cap = host_cap if xp is np else _device_cap(host_cap, "bond pairing admission")
    nb = _frequency_batch(nmat, nvol, arch.ND, host_cap,
                          device_cap=(None if xp is np else device_cap))
    table = _eb.estimate_pair_memory(
        nmat=nmat, ntau=(axF.n_tau if use_ir else 0), nfreq=nfreq, nvol=nvol, norb=norb,
        B=arch.B, num_eigenvalues=ctl.num_eigenvalues, residency="auto", ir=use_ir,
        L_B=L_B, n_channels=1, in_process=False, nb=nb,
        keep_const=ctl.ir_keep_static_chi and use_ir)
    logger.info("bond pairing admission (post-processing, GiB): %s",
                {k: round(v / _eb._GIB, 3) for k, v in table.rows.items()})
    # raises MemoryError with the full table when nothing fits
    residency = table.choose("host" if xp is np else "auto", host_cap, device_cap,
                             shared=(xp is np))

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
        lam, gap_w, eigenvalues_all, eigenvalue_match, note, leakage, sector_weights = \
            _ed.run_leading_eigenproblem(
                kernel.matvec, kernel.gap_shape, ctl.as_eli_param(), pairing_type, phi0=phi0,
                seed_vec=seed_vec, use_ir=use_ir, axF=axF, nmat=nmat,
                logger_label="bond pairing kernel", parity_leakage_policy="refuse",
                parity_leakage_tol=ctl.resolved_parity_leakage_tol)
        residency_used = kernel.residency

    # --- outputs (the on-site dynamic file set plus the bond provenance) ----
    out_dir = input_dict["file"]["output"]["path_to_output"]
    os.makedirs(out_dir, exist_ok=True)
    # the measured parity leakage is recorded in BOTH outputs (spec 4.4), so a
    # run accepted under a raised parity_leakage_tol carries its own evidence
    head = ["bond_channels=true", "residency=" + residency_used]
    if leakage is not None:
        head.append("parity_leakage={:.6e}".format(leakage))
    _ed.write_eigenvalue_file(
        os.path.join(out_dir, eli_param.get("output_eigenvalue", "eigenvalue.dat")),
        lam, eigenvalues_all, eigenvalue_match, note, header_lines=head,
        sector_weights=sector_weights)
    extra = {"bond_channels": True, "bond_delta_r": arch.delta_r, "bond_reverse": arch.reverse,
             "bond_archive": arch_path, "bond_residency": residency_used,
             "gap_bond_projection": _eb.gap_bond_projection(gap_w, view, (Nx, Ny, Nz))}
    if leakage is not None:
        extra["bond_parity_leakage"] = float(leakage)
    if use_ir:
        extra.update({"matsubara_basis": "ir", "ir_tol": axF.eps, "ir_wmax": axF.wmax,
                      "ir_L": axF.L, "bond_ir_fit_residual_rel": vertex.fit_residual_rel})
    _ed.write_dynamic_outputs(out_dir, gap_w, lam, T, pairing_type, kx, ky, kz, beta,
                              gap_file=eli_param.get("output_gap", "gap.dat"),
                              extra_meta=extra, sector_weights=sector_weights)
    return lam


# =============================================================================
# Entry 2: in-process at the end of the FLEX solve (spec 7)
# =============================================================================

def _inprocess_axes(solver):
    """The fermionic and bosonic IR axes of an in-process pairing run.

    ``ir_wmax`` defaults to three times the sum of the full bandwidth and the
    largest interaction matrix element, both read from the solve's OWN objects
    (the diagonalised H0 and the momentum-space interaction), so the auto value
    does not depend on files this entry never reads.
    """
    from . import backend as _bk
    from .ir_axis import IRAxis
    ctl = solver._pairing_controls
    beta = 1.0 / solver.T
    wmax = ctl.ir_wmax
    if wmax is None:
        ew = _bk.to_host(solver.H0_eigenvalue)
        band = 2.0 * float(np.abs(ew).max())
        inter = solver.ham_info.ham_inter_q
        u = 0.0 if inter is None else float(np.abs(_bk.to_host(inter)).max())
        wmax = 3.0 * (band + u)
        logger.info("bond pairing IR: auto ir_wmax = %.6g", wmax)
    axF = IRAxis(beta=beta, wmax=float(wmax), eps=ctl.ir_tol, statistics="F")
    axB = IRAxis(beta=beta, wmax=float(wmax), eps=ctl.ir_tol, statistics="B")
    return axF, axB


def pairing_preflight(solver):
    """Solve-entry admission (spec 7): refuse BEFORE the first FLEX map when
    even the streaming residency cannot fit next to the FLEX solve's own
    predicted peak.

    The caps are the section-8 caps measured now MINUS what the bond gate has
    already told us it will hold when the pairing step starts
    (``solver._bond_est``). The table is charged the solve's OWN dressing
    batch (``solver._bond_nb``), because the in-process vertex build re-dresses
    ``chibar`` in batches of exactly that size; the post-SCF admission is
    kernel-only (the vertex is built by then) and charges none. Raises
    ``ValueError``; the post-SCF admission re-measures and is a per-channel
    error there, not an exception.
    """
    from . import backend as _bk
    from . import eliashberg_bond as _eb
    ctl = solver._pairing_controls
    nmat, nvol, norb = solver.nmat, solver.lattice.nvol, solver.norb
    B = solver._bond_view.n_channels
    use_ir = ctl.matsubara_basis == "ir"
    axF, axB = _inprocess_axes(solver) if use_ir else (None, None)
    table = _eb.estimate_pair_memory(
        nmat=nmat, ntau=(axF.n_tau if use_ir else 0), nfreq=(axF.n_freq if use_ir else nmat),
        nvol=nvol, norb=norb, B=B, num_eigenvalues=ctl.num_eigenvalues, residency="auto",
        ir=use_ir, L_B=(axB.L if use_ir else 0), n_channels=len(ctl.pairing_types),
        in_process=True, nb=solver._bond_nb,
        keep_const=ctl.ir_keep_static_chi and use_ir)
    host_cap_now = (ctl.bond_memory_cap_gb * _eb._GIB) if ctl.bond_memory_cap_gb \
        else 0.8 * _eb._host_available_bytes()
    # estimate_bond_memory's "peak": 1.25 x (persistent + the largest phase)
    flex_host_peak = float(solver._bond_est["peak"])
    host_cap = host_cap_now - flex_host_peak
    device_cap = host_cap
    if getattr(solver, "use_gpu", False):
        available = _bk.device_available_bytes()
        if available is not None:
            device_need = solver._bond_est.get("device_need")
            if device_need is None:
                # the FLEX device admission runs inside the first map, i.e.
                # after this point; say so rather than implying the 0 below is
                # a measurement
                logger.warning(
                    "longitudinal_bond_pairing preflight: the FLEX device peak is unknown at "
                    "solve entry, charged 0; the post-SCF admission re-measures it and reports "
                    "a shortfall as a per-channel error instead of a refusal here")
                device_need = 0.0
            device_cap = 0.9 * float(available) - float(device_need)
    try:
        residency = table.choose("auto", host_cap, device_cap,
                                 shared=not getattr(solver, "use_gpu", False))
    except MemoryError as exc:
        raise ValueError(
            "longitudinal_bond_pairing: the pairing step would not fit after the FLEX solve "
            "({}); set [eliashberg] matsubara_basis = 'ir' or raise bond_memory_cap_gb"
            .format(exc))
    logger.info("bond pairing preflight at solve entry: residency %s would fit; rows (GiB) %s",
                residency, {k: round(v / _eb._GIB, 3) for k, v in table.rows.items()})


def _record_backstop(green_info, etas, phase, exc):
    """Record a failure the per-phase handlers did not see, without erasing a
    more specific message one of them already stored."""
    msg = "{}: {}: {}".format(phase, type(exc).__name__, exc)
    logger.error("longitudinal_bond_pairing failed outside the per-phase handlers (%s): %s",
                 phase, exc)
    for eta in etas:
        key = "pairing_{}_error".format(eta)
        if key not in green_info:
            green_info[key] = msg


def _requested_types(solver):
    ctl = getattr(solver, "_pairing_controls", None)
    return () if ctl is None else tuple(ctl.pairing_types)


def run_inprocess_pairing(solver, store, dev, green_kw, beta, green_info):
    """The in-process pairing step at the end of the FLEX solve (spec 7).

    NON-THROWING: every ``Exception`` and every device out-of-memory type is
    caught per phase and recorded as ``pairing_<type>_error`` in ``green_info``
    (``KeyboardInterrupt`` / ``SystemExit`` propagate), so a failure here never
    reaches ``_solve_phase_b``'s handler and never costs the FLEX results.
    A failure in a stage SHARED by several channels fails all of them with the
    same root cause; the later per-channel stages fail that channel only.
    """
    # the handler's OWN setup runs inside the guard: an import or a device
    # probe that fails here would otherwise propagate and cost the FLEX
    # results. Python evaluates the ``except`` expression at match time, so
    # rebinding ``caught`` inside the try is what makes the device
    # out-of-memory types part of the match once they are known.
    caught = (Exception,)
    try:
        from . import backend as _bk
        caught = (Exception,) + tuple(_bk._oom_error_types())
        _run_inprocess_pairing(solver, store, dev, green_kw, beta, green_info)
    except caught as exc:
        # the phase handlers below cover every failure the design foresees;
        # this is the backstop that keeps the promise for the one it does not
        _record_backstop(green_info, _requested_types(solver), "pairing", exc)


def _run_inprocess_pairing(solver, store, dev, green_kw, beta, green_info):
    from . import backend as _bk
    from . import bond_channels as _bc
    from . import eliashberg_bond as _eb
    from . import eliashberg_dynamic as _ed
    ctl = solver._pairing_controls
    xp = dev.xp
    nmat, nvol, norb = solver.nmat, solver.lattice.nvol, solver.norb
    nd = norb * norb
    shape = tuple(int(x) for x in solver.lattice.shape)
    Nx, Ny, Nz = shape
    view = solver._bond_view
    use_ir = ctl.matsubara_basis == "ir"
    types = ctl.pairing_types
    label_state = ("last_map_chi / final_green" if solver.scf_converged
                   else "mixed: last-map chi, final green")
    if not solver.scf_converged:
        logger.warning("longitudinal_bond_pairing: the FLEX solve did not converge; the pairing "
                       "eigenvalues are a diagnostic of a mixed state (%s)", label_state)
    oom = tuple(_bk._oom_error_types())
    caught = (Exception,) + oom

    def fail(etas, phase, exc):
        msg = "{}: {}: {}".format(phase, type(exc).__name__, exc)
        logger.error("longitudinal_bond_pairing (%s) failed in phase %s: %s",
                     ",".join(etas), phase, exc)
        for eta in etas:
            green_info["pairing_{}_error".format(eta)] = msg

    try:
        axF, axB = _inprocess_axes(solver) if use_ir else (None, None)
        # the sc.py Green layout, from the SAME array green.npz carries, so the
        # two entries build G2 from bit-identical data (spec 7, test 10.2.3)
        g_host = _bk.to_host(green_kw)[0].reshape(
            nmat, Nx, Ny, Nz, norb, norb).transpose(4, 5, 1, 2, 3, 0).copy()
        if use_ir:
            g_host = _ed._ir_compress(g_host, axF, nmat, "green")
        G2 = xp.asarray(_ed.calc_g2_dynamic(g_host, beta))
        del g_host
        nfreq = axF.n_freq if use_ir else nmat
    except caught as exc:
        fail(types, "setup", exc)
        return

    # IR carries both channels through ONE accumulator (4.5); the uniform grid
    # builds Gamma in the store's W slot, so one channel at a time
    groups = [types] if use_ir else [(eta,) for eta in types]
    for group in groups:
        vertices = None
        try:
            if use_ir:
                # the IR coefficients replace the W slot rather than joining it
                store.release_slot("W")
            acc = _eb.PairVertexAccumulator(
                dev, pairing_types=group, nb=solver._bond_nb, nmat=nmat, nvol=nvol, nd=nd,
                spatial_shape=shape, ir=((axF, axB) if use_ir else None),
                out_store=(None if use_ir else store), out_slot="W",
                ir_keep_static=ctl.ir_keep_static_chi,
                workers=getattr(solver, "fft_workers", 1))

            def stages(a):
                a.add_dressed(store, cond_tol=_bc._BOND_COND_FLOOR,
                              guard_freqs=solver.longitudinal_bond_guard_freqs)

            stages(acc)
            vertices = acc.finish(ir_fit_tol=ctl.ir_fit_tol, stage_callable=stages)
            del acc
        except caught as exc:
            fail(group, "vertex", exc)
            continue
        for eta in group:
            K = None
            try:
                # the previous channel's device arrays are gone by now; give
                # them back to the driver so this admission measures the
                # truth. On the numpy backend there is no device pool to give
                # back to, and a cupy installed next to a CPU run belongs to
                # somebody else -- do not touch it.
                if xp is not np:
                    _bk.free_device_pool()
                host_cap = (ctl.bond_memory_cap_gb * _eb._GIB) if ctl.bond_memory_cap_gb \
                    else 0.8 * _eb._host_available_bytes()
                device_cap = host_cap if xp is np else _device_cap(
                    host_cap, "longitudinal_bond_pairing ({}) admission".format(eta))
                table = _eb.estimate_pair_memory(
                    nmat=nmat, ntau=(axF.n_tau if use_ir else 0), nfreq=nfreq, nvol=nvol,
                    norb=norb, B=view.n_channels, num_eigenvalues=ctl.num_eigenvalues,
                    residency="auto", ir=use_ir, L_B=(axB.L if use_ir else 0),
                    n_channels=1, in_process=True,
                    keep_const=ctl.ir_keep_static_chi and use_ir)
                V_inst = _eb.instantaneous_vertex(solver._bond_S, solver._bond_C, nd, eta, shape)
                K = _eb.BondPairKernel(
                    vertices[eta], G2, view, xp=xp, spatial_shape=shape, norb=norb, beta=beta,
                    nfreq=nfreq, V_inst=V_inst, axF=axF, residency="auto", admission=table,
                    host_cap=host_cap, device_cap=device_cap,
                    workers=getattr(solver, "fft_workers", 1))
                logger.info("longitudinal_bond_pairing (%s): kernel residency %s", eta, K.residency)
                kx, ky, kz = (np.linspace(0.0, 2.0 * np.pi, n, endpoint=False) for n in shape)
                phi0, seed = _ed.build_seed(ctl.as_eli_param(), eta, norb, kx, ky, kz,
                                            K.gap_shape, use_ir, axF, nmat)
                lam, gap_w, evs, match, note, leakage, sector_weights = \
                    _ed.run_leading_eigenproblem(
                        K.matvec, K.gap_shape, ctl.as_eli_param(), eta, phi0=phi0,
                        seed_vec=seed, use_ir=use_ir, axF=axF, nmat=nmat,
                        logger_label="bond pairing kernel ({})".format(eta),
                        parity_leakage_policy="refuse",
                        parity_leakage_tol=ctl.resolved_parity_leakage_tol)
                meta = {"bond_channels": True,
                        "bond_delta_r": np.asarray(view.delta_r, dtype=np.int64),
                        "bond_reverse": np.asarray(view.reverse, dtype=np.int64),
                        "bond_residency": K.residency,
                        "scf_converged": bool(solver.scf_converged),
                        "scf_iterations": int(solver.scf_iterations),
                        "state": label_state,
                        "matsubara_basis": ctl.matsubara_basis,
                        "gap_bond_projection": _eb.gap_bond_projection(gap_w, view, shape),
                        # the SAME two npz keys the file entry writes, so a
                        # reader sees one key set from both entries
                        "gap_sector_weights": np.array(
                            [float(sector_weights[label])
                             for label in _ed._SECTOR_LABELS]),
                        "gap_sector_labels": np.array(_ed._SECTOR_LABELS)}
                if leakage is not None:
                    meta["bond_parity_leakage"] = float(leakage)
                if use_ir:
                    meta.update({"ir_tol": axF.eps, "ir_wmax": axF.wmax, "ir_L": axF.L,
                                 "bond_ir_fit_residual_rel": vertices[eta].fit_residual_rel})
                green_info["pairing_{}_eigenvalue".format(eta)] = float(lam)
                green_info["pairing_{}_eigenvalues".format(eta)] = \
                    None if evs is None else np.asarray(evs)
                green_info["pairing_{}_eigenvalue_match".format(eta)] = \
                    None if match is None else np.asarray(match)
                green_info["pairing_{}_gap".format(eta)] = gap_w
                green_info["pairing_{}_meta".format(eta)] = meta
                green_info["pairing_{}_note".format(eta)] = note
                logger.info("longitudinal_bond_pairing (%s): leading eigenvalue %.8e (%s)",
                            eta, float(lam), label_state)
            except caught as exc:
                fail((eta,), "kernel/solver", exc)
            finally:
                # a FAILED channel must not keep its kernel (and the device
                # blocks it hoisted) alive while the next channel calls
                # free_device_pool and measures what is available
                K = None
        del vertices


def _tmp_name(path):
    """``dir/x.npz`` -> ``dir/.x.<pid>-<token>.tmp.npz``.

    Next to its target (so the publishing rename stays inside one filesystem),
    hidden by the leading dot, and UNIQUE per process and per call: a fixed
    ``x.tmp.npz`` collides whenever two runs share an output directory (a
    parameter sweep, a restarted job), and one run's rename can then publish
    the other's half-written file. The ``.npz`` tail is preserved, so numpy's
    savez suffix rule is still satisfied."""
    head, base = os.path.split(path)
    root, ext = os.path.splitext(base)
    token = secrets.token_hex(4)
    return os.path.join(head, ".{}.{}-{}.tmp{}".format(root, os.getpid(), token, ext))


def write_pairing_outputs(solver, info_outputfile, green_info, path_to_output):
    """The three files of every successful pairing channel (spec 7).

    Called from ``save_results`` AFTER every FLEX artifact. Each channel is
    written under temporary names and then renamed in a fixed order (npz, gap,
    eigenvalue); a write failure removes the temporaries and publishes nothing,
    a rename failure leaves what was already renamed and reports the channel as
    PARTIALLY published. Either way the message lands in
    ``pairing_<type>_error`` and the next channel is still attempted -- this
    function never raises.
    """
    try:
        _write_pairing_outputs(solver, info_outputfile, green_info, path_to_output)
    except Exception as exc:
        _record_backstop(green_info, _requested_types(solver), "outputs", exc)


def _write_pairing_outputs(solver, info_outputfile, green_info, path_to_output):
    from . import eliashberg_dynamic as _ed
    from . import flex_bond as _fb
    ctl = solver._pairing_controls
    shape = tuple(int(x) for x in solver.lattice.shape)
    kx, ky, kz = (np.linspace(0.0, 2.0 * np.pi, n, endpoint=False) for n in shape)
    for eta in ctl.pairing_types:
        if "pairing_{}_error".format(eta) in green_info \
                or "pairing_{}_gap".format(eta) not in green_info:
            continue
        # resolved INSIDE the guard: a name or a result this channel is missing
        # is this channel's failure, not the next channel's
        names, tmps = {}, {}
        try:
            for kind in ("eliashberg_bond", "gap_bond", "eigenvalue_bond"):
                key = "{}_{}".format(kind, eta)
                fn = str(info_outputfile.get(key, _fb._DEFAULT_FILES[key]))
                if kind == "eliashberg_bond" and not fn.endswith(".npz"):
                    fn += ".npz"
                names[kind] = os.path.join(str(path_to_output), fn)
                tmps[kind] = _tmp_name(os.path.abspath(names[kind]))
            for kind in names:
                # a configured name may carry a subdirectory; the temporary
                # lives next to its target, so both need that directory
                os.makedirs(os.path.dirname(tmps[kind]), exist_ok=True)
            meta = green_info["pairing_{}_meta".format(eta)]
            extra = dict(meta)
            evs = green_info.get("pairing_{}_eigenvalues".format(eta))
            if evs is not None:
                extra["eigenvalues_all"] = evs
            _ed.write_dynamic_outputs(
                str(path_to_output), green_info["pairing_{}_gap".format(eta)],
                green_info["pairing_{}_eigenvalue".format(eta)], solver.T, eta,
                kx, ky, kz, 1.0 / solver.T,
                # the ABSOLUTE temporaries: os.path.join(output_dir, abs_path)
                # is abs_path, so a configured name with a subdirectory (or an
                # absolute one) lands where it was asked for and not next to
                # path_to_output under its basename
                gap_file=tmps["gap_bond"], npz_file=tmps["eliashberg_bond"],
                extra_meta=extra)
            _ed.write_eigenvalue_file(
                tmps["eigenvalue_bond"], green_info["pairing_{}_eigenvalue".format(eta)],
                evs, green_info.get("pairing_{}_eigenvalue_match".format(eta)),
                green_info.get("pairing_{}_note".format(eta)),
                header_lines=(
                    ["bond_channels=true",
                     "scf_converged={}".format(str(meta["scf_converged"]).lower()),
                     "state={}".format(meta["state"]),
                     "residency={}".format(meta["bond_residency"])]
                    + (["parity_leakage={:.6e}".format(meta["bond_parity_leakage"])]
                       if "bond_parity_leakage" in meta else [])),
                # the weights travelled here inside meta (they are npz keys);
                # rebuild the dict the header writer takes
                sector_weights=(
                    {str(k): float(v) for k, v in zip(meta["gap_sector_labels"],
                                                      meta["gap_sector_weights"])}
                    if "gap_sector_weights" in meta else None))
        except Exception as exc:
            for t in tmps.values():
                try:
                    os.remove(t)
                except OSError:
                    pass
            green_info["pairing_{}_error".format(eta)] = \
                "outputs: {}: {}".format(type(exc).__name__, exc)
            logger.error("longitudinal_bond_pairing (%s): writing the outputs failed: %s",
                         eta, exc)
            continue
        done = []
        try:
            for kind in ("eliashberg_bond", "gap_bond", "eigenvalue_bond"):
                os.replace(tmps[kind], names[kind])
                done.append(kind)
        except OSError as exc:
            missing = [k for k in tmps if k not in done]
            for k in missing:
                try:
                    os.remove(tmps[k])
                except OSError:
                    pass
            green_info["pairing_{}_error".format(eta)] = \
                "outputs: partially published, missing {}: {}".format(missing, exc)
            logger.error("longitudinal_bond_pairing (%s): PARTIALLY published (missing %s): %s",
                         eta, missing, exc)
            continue
        logger.info("save_results: pairing (%s) outputs %s", eta, ", ".join(names.values()))
