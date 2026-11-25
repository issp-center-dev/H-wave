"""The initial-state writer and its check (design section 9).

:func:`initial_state` seeds one candidate on one ``SubShape`` (design
section 7) and writes the Green function UHFk reads as ``[file.input]
initial``; with a ``site_map`` it writes the ``Initial`` file of a UHFr
input instead. :func:`check_initial_state` checks such a file against the
input that will read it, before ``hwave`` runs.

**Keys** (design section 9.1). ``prod(SubShape) > 1``: ``green_sublattice``
(the member UHFk reads, ``UHFk._read_green_from_data``) and ``green``, its
unfolded view (:func:`initial_facade.unfold_layout`, the solver's
``_deflate_green``; lossy for a seed that breaks the translations of the
primitive cell, and never the member read for it, since the check refuses
a ``SubShape`` mismatch); ``prod(SubShape) = 1``: ``green``. Always the
convention markers ``green_convention = "green_slot_first"`` and
``momentum_convention = "e_plus_ikR"`` in the solver's string form, and
the provenance members :data:`PROVENANCE_KEYS`.

**Pickle-free provenance.** Every provenance member is a 0-d numpy unicode
array holding the canonical JSON of its value (sorted keys, no
whitespace): ``null``, numbers, strings, lists and the ``input_sha256``
object. Files are read with ``allow_pickle=False``; an object member is
refused.

**Why the check matters.** UHFk starts from a zero Green function when
the member it looks for is missing, and reads ``green_sublattice`` of
another ``SubShape`` with the same array shape without a word (2 x 2 ->
4 x 1 on a 4 x 4 lattice). The check refuses both, and every other
mismatch of the provenance with the input, naming the field.

**Bytes.** Members are written in sorted key order into a zip archive
with fixed time stamps, so two writes of one seed give identical bytes
within one numpy version. Across numpy and BLAS builds only the keys and
member contents are comparable (eigh and FFT are not bitwise).

**Where the solver looks** (read_input_k.py: ``os.path.join([file.input]
path_to_input, [file.input] initial)``, relative to the working directory
of ``hwave``, ``path_to_input`` defaulting to ``""``):
:func:`solver_initial_path` and its inverse :func:`initial_key_value`
assume ``hwave`` runs in the TOML's directory. No other code of
``hwave.initial`` joins that path.

**UHFr.** The model input is a UHFk TOML of the same lattice; the
``site_map`` gives, for UHFr site ``i``, its UHFk cell ``(ix, iy, iz)``
in ``CellShape`` and physical orbital ``a``. With the evaluated
``SubShape`` B, cell R lies in supercell ``R_B = R // B`` at in-cell
offset ``ir``, folded orbital ``a' = a + norb_orig * ir``, and
``G_UHFr[(i, s), (j, t)] = <c^dag_{is} c_{jt}> = G[r(R_B,i - R_B,j), s,
a'_i, t, a'_j]`` (the stored layout ``G[r, s, a, t, b] = <c^dag_{sa}(r)
c_{tb}(0)>``).
"""
import hashlib
import json
import math
import numbers
import os
import pathlib
import tempfile
import zipfile
from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np

from hwave.solver import initial_facade

__all__ = [
    "InitialStateError", "ProvenanceMissingError", "CheckResult",
    "INITIAL_STATE_SCHEMA", "PROVENANCE_KEYS", "CONVENTION_MEMBERS",
    "NO_PROVENANCE_WARNING", "SO_UHFR_REFUSAL", "HERMITE_TOL", "EIG_TOL",
    "TRACE_RTOL", "initial_state", "check_initial_state", "check_green",
    "solver_initial_path", "initial_key_value", "uhfr_green", "write_uhfr_initial",
]

#: the provenance schema this module writes and reads
INITIAL_STATE_SCHEMA = 1

#: the provenance members (design section 9.1), in sorted order
PROVENANCE_KEYS = tuple(sorted((
    "initial_state_schema", "run_id", "candidate_id", "label", "q", "phase", "sign",
    "subshape", "cellshape", "norb_orig", "spin_orbital", "index_convention", "two_sz",
    "ncond", "seed", "h_used", "seed_base", "hwave_version", "input_sha256")))

#: the convention markers, in the solver's string form (UHFk._save_green)
CONVENTION_MEMBERS = {"green_convention": "green_slot_first",
                      "momentum_convention": "e_plus_ikR"}

#: the stored index order: spin-block ``G[r, s, a, t, b]`` (normal mode)
#: or interleaved ``2 * orb + spin`` with ``ns = 1`` (spin-orbital mode)
_INDEX_CONVENTION = {False: "spin_block", True: "interleaved"}

#: design section 9.3, verbatim
NO_PROVENANCE_WARNING = ("no provenance: SubShape, SO convention and electron sectors "
                         "cannot be confirmed beyond shape")

#: design section 9.2, verbatim
SO_UHFR_REFUSAL = "UHFr has no spin-orbital mode; use UHFk with enable_spin_orbital"

#: ``max |G(r)_ab - conj(G(-r)_ba)|`` above this is not Hermitian
HERMITE_TOL = 1e-10
#: every eigenvalue of ``G(k)`` must lie in ``[-EIG_TOL, 1 + EIG_TOL]``
EIG_TOL = 1e-10
#: ``|N - Ncond|`` (and each 2Sz sector) within this times ``max(1, Ncond)``
TRACE_RTOL = 1e-8

#: fixed zip time stamps, so equal content gives equal bytes
_ZIP_TIME = (1980, 1, 1, 0, 0, 0)


class InitialStateError(ValueError):
    """An initial-state file, or a seed, is refused; the message names the
    field or the check."""


class ProvenanceMissingError(InitialStateError):
    """``require_provenance=True`` and the file has no provenance members
    (written by an older tool, by hand, or a plain ``hwave`` output)."""


@dataclass(frozen=True)
class CheckResult:
    """The outcome of :func:`check_initial_state` (refusals raise).

    ``verified``: the provenance matches the input. ``key``: the member the
    solver will read. ``warnings``: e.g. the no-provenance warning."""
    verified: bool
    key: str
    warnings: Tuple[str, ...] = ()


# --------------------------------------------------------------------------
# where the solver looks
# --------------------------------------------------------------------------

def _abs(path):
    return pathlib.Path(os.path.abspath(os.fspath(path)))


def _read_doc(config_path):
    import tomli
    with open(config_path, "rb") as f:
        return tomli.load(f)


def _finput(doc):
    """``[file.input]`` by exact key, as ``hwave`` reads it."""
    f = doc.get("file", {})
    f = f.get("input", {}) if isinstance(f, dict) else {}
    return f if isinstance(f, dict) else {}


def _path_to_input(finput):
    v = finput.get("path_to_input", "")
    if not isinstance(v, str):
        raise InitialStateError("[file.input] path_to_input must be a string, got {!r}"
                                .format(v))
    return v


def solver_initial_path(config_path) -> Optional[pathlib.Path]:
    """The file ``hwave config_path`` reads as ``[file.input] initial``,
    run in the TOML's directory: ``os.path.join(<TOML dir>, path_to_input,
    initial)`` (an absolute ``initial`` as it is). ``None`` without an
    ``initial`` key."""
    config_path = _abs(config_path)
    finput = _finput(_read_doc(config_path))
    if "initial" not in finput:
        return None
    initial = finput["initial"]
    if not isinstance(initial, str):
        raise InitialStateError("[file.input] initial must be a string, got {!r}"
                                .format(initial))
    return pathlib.Path(os.path.join(str(config_path.parent),
                                     os.path.join(_path_to_input(finput), initial)))


def initial_key_value(config_path, npz_path) -> str:
    """The ``[file.input] initial`` value under which ``hwave config_path``
    (run in the TOML's directory) reads ``npz_path``: its path relative to
    ``<TOML dir>/path_to_input``, in POSIX form."""
    config_path = _abs(config_path)
    finput = _finput(_read_doc(config_path))
    base = os.path.join(str(config_path.parent), _path_to_input(finput))
    rel = os.path.relpath(os.path.abspath(os.fspath(npz_path)), os.path.abspath(base))
    return pathlib.PurePath(rel).as_posix()


# --------------------------------------------------------------------------
# the checks on a stored Green function (design section 9.2)
# --------------------------------------------------------------------------

def _neg_r(shape):
    """The index of ``-r`` for each C-order index ``r`` (z fastest)."""
    nx, ny, nz = shape
    ix, iy, iz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    return (((-iz) % nz) + nz * (((-iy) % ny) + ny * ((-ix) % nx))).ravel()


def check_green(G, problem, *, what="green") -> None:
    """The design section 9.2 checks of the stored Green function ``G``
    ``(nvol, ns, norb, ns, norb)`` against ``problem`` (an
    :class:`initial_facade.UHFkProblem`): the shape, Hermiticity ``G(r)_ab =
    conj(G(-r)_ba)``, ``0 <= eig G(k) <= 1`` at every k, the trace
    ``Ncond``, the 2Sz sectors ``(Ncond +- 2Sz)/2`` and the solver blocks.
    Raises :class:`InitialStateError` naming the check and ``what``."""
    p = problem
    G = np.asarray(G)
    expected = (p.nvol, p.ns, p.norb, p.ns, p.norb)
    if G.shape != expected:
        mode = "spin-orbital" if p.enable_spin_orbital else "normal"
        raise InitialStateError(
            "{}: shape {} does not fit this input: expected {} (k-grid {} = CellShape {} "
            "/ SubShape {}, ns = {} in {} mode, {} orbitals of the folded cell)".format(
                what, G.shape, expected, list(p.shape), list(p.cellshape),
                list(p.subshape), p.ns, mode, p.norb))
    if not np.all(np.isfinite(G)):
        raise InitialStateError("{}: not finite".format(what))
    nd = p.nd
    g = G.reshape(p.nvol, nd, nd)
    herm = float(np.max(np.abs(g - np.conj(np.swapaxes(g[_neg_r(p.shape)], -1, -2)))))
    if herm > HERMITE_TOL:
        raise InitialStateError(
            "{}: not Hermitian: max |G(r)_ab - conj(G(-r)_ba)| = {:.3g} > {:g}".format(
                what, herm, HERMITE_TOL))
    gk = initial_facade.to_k(G.reshape(tuple(p.shape) + (nd, nd)))
    ev = np.linalg.eigvalsh(gk)
    lo, hi = float(ev.min()), float(ev.max())
    if lo < -EIG_TOL or hi > 1.0 + EIG_TOL:
        raise InitialStateError(
            "{}: eigenvalue bounds: the eigenvalues of G(k) lie in [{:.12g}, {:.12g}], "
            "not within [0, 1] (tolerance {:g})".format(what, lo, hi, EIG_TOL))
    tol = TRACE_RTOL * max(1.0, abs(float(p.ncond)))
    diag = np.real(np.diagonal(g[0])) * p.nvol
    n = float(np.sum(diag))
    if abs(n - float(p.ncond)) > tol:
        raise InitialStateError(
            "{}: trace: the electron number is {:.12g}, the input's Ncond is {}".format(
                what, n, p.ncond))
    if p.two_sz is not None and not p.enable_spin_orbital:
        up = float(np.sum(diag[:p.norb]))
        dn = float(np.sum(diag[p.norb:]))
        want_up = (float(p.ncond) + p.two_sz) / 2
        want_dn = (float(p.ncond) - p.two_sz) / 2
        if abs(up - want_up) > tol or abs(dn - want_dn) > tol:
            raise InitialStateError(
                "{}: 2Sz sectors: spin up {:.12g}, spin down {:.12g}; 2Sz = {} needs {} "
                "and {}".format(what, up, dn, p.two_sz, want_up, want_dn))
    # a component outside the solver blocks above seed.BLOCK_RTOL of the norm
    from hwave.initial.seed import block_check
    total = float(np.linalg.norm(gk))
    msgs = block_check(gk, p.blocks, ns=p.ns, norb=p.norb) if total > 0 else []
    if msgs:
        raise InitialStateError("{}: block check: {}".format(what, "; ".join(msgs)))


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

def _json_value(v):
    if v is None or isinstance(v, (bool, str)):
        return v
    if isinstance(v, numbers.Integral):
        return int(v)
    if isinstance(v, numbers.Real):
        f = float(v)
        if not math.isfinite(f):
            raise ValueError("provenance values must be finite, got {!r}".format(v))
        return f
    if isinstance(v, (list, tuple)):
        return [_json_value(x) for x in v]
    if isinstance(v, dict) or hasattr(v, "items"):
        return {str(k): _json_value(x) for k, x in v.items()}
    raise TypeError("cannot encode {!r} as provenance JSON".format(v))


def _encode(v):
    from hwave.initial.patterns import canonical_json
    return np.array(canonical_json(_json_value(v)))


def _write_npz(path, members):
    """``members`` as an npz archive at ``path``: sorted names, fixed zip
    time stamps, no pickle; written to a temporary file and renamed."""
    path = _abs(path)
    fd, tmp = tempfile.mkstemp(prefix=".initial-", suffix=".npz", dir=str(path.parent))
    os.close(fd)
    try:
        with zipfile.ZipFile(tmp, mode="w", compression=zipfile.ZIP_STORED,
                             allowZip64=True) as zf:
            for name in sorted(members):
                info = zipfile.ZipInfo(name + ".npy", date_time=_ZIP_TIME)
                info.compress_type = zipfile.ZIP_STORED
                info.external_attr = 0o600 << 16
                with zf.open(info, mode="w", force_zip64=True) as fid:
                    np.lib.format.write_array(fid, np.asanyarray(members[name]),
                                              allow_pickle=False)
        # mkstemp creates the file 0600; give it the mode open() would
        umask = os.umask(0)
        os.umask(umask)
        os.chmod(tmp, 0o666 & ~umask)
        os.replace(tmp, str(path))
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise
    return path


# --------------------------------------------------------------------------
# the writer
# --------------------------------------------------------------------------

def _toml_of(input_dir_or_model):
    """The TOML path of the model input (refusing a ``ModelFacts`` and a
    UHFr TOML, naming the route for UHFr)."""
    from hwave.initial.facts import ModelFacts, _resolve_toml
    if isinstance(input_dir_or_model, ModelFacts):
        raise TypeError("initial_state needs the input path (a TOML or a directory): it "
                        "loads the model on the evaluated SubShape")
    if not isinstance(input_dir_or_model, (str, os.PathLike)):
        raise TypeError("input_dir_or_model: a TOML file or an input directory, got {!r}"
                        .format(input_dir_or_model))
    p = pathlib.Path(os.fspath(input_dir_or_model))
    found = sorted(p.glob("*.toml")) if p.is_dir() else [p] if p.is_file() else []
    if len(found) == 1:
        mode = _read_doc(found[0]).get("mode", {})
        if isinstance(mode, dict) and mode.get("mode") == "UHFr":
            raise ValueError(
                "{}: a UHFr input has no lattice geometry; give the UHFk input of the "
                "same lattice as the model, and the UHFr input as uhfr_toml= with a "
                "site_map= (UHFr site -> UHFk cell and orbital)".format(found[0]))
    return _resolve_toml(p)


def _run_of(candidate_or_run, subshape):
    from hwave.initial.candidates import Candidate, CandidateRun
    if isinstance(candidate_or_run, CandidateRun):
        run = candidate_or_run
        if subshape is not None and tuple(int(x) for x in subshape) != tuple(run.subshape):
            raise ValueError(
                "subshape {} conflicts with the run's subshape {} ({}); pass None or the "
                "run's own".format(list(subshape), list(run.subshape), run.run_id))
        return run
    if isinstance(candidate_or_run, Candidate):
        B = candidate_or_run.subshape if subshape is None else subshape
        return CandidateRun(candidate_or_run, tuple(B))
    raise TypeError("candidate_or_run: a Candidate or a CandidateRun, got {!r}"
                    .format(candidate_or_run))


def _seed(facts, run, para):
    """The seed of ``run`` on the cell of ``facts`` (design section 7)."""
    from hwave.initial.cells import representable, to_qfrac
    from hwave.initial.hints import t_hint_bare
    from hwave.initial.patterns import realise
    from hwave.initial.seed import seed_pattern
    if para is not None:
        raise NotImplementedError(
            "para=: seeding on the para Hartree-Fock state arrives in a later release; "
            "this release seeds on the bare H0 (para=None)")
    cand = run.candidate
    spec = cand.pattern
    B = tuple(run.subshape)
    if spec.q is not None:
        ok, reason = representable(to_qfrac(spec.q), B)
        if not ok:
            raise ValueError("{}: not representable on SubShape {}: {}".format(
                run.run_id, list(B), reason))
    p = facts.problem
    O = None if spec.modulation == "none" else realise(spec, facts, B)
    res = seed_pattern(facts, p.ham_trans_k, O, T_run=p.T, W0=t_hint_bare(p).W0,
                       amplitude=cand.amplitude_requested, seed_base="bare")
    if not res.tried:
        raise InitialStateError("{}: refused: {}".format(run.run_id, res.not_tried_reason))
    return res


def _provenance(facts, run, res):
    import hwave
    cand = run.candidate
    spec = cand.pattern
    p = facts.problem
    return {
        "initial_state_schema": INITIAL_STATE_SCHEMA,
        "run_id": run.run_id, "candidate_id": cand.candidate_id, "label": cand.label,
        "q": None if spec.q is None else list(spec.q), "phase": spec.phase,
        "sign": spec.sign, "subshape": list(run.subshape),
        "cellshape": list(facts.cellshape), "norb_orig": int(p.norb_orig),
        "spin_orbital": bool(p.enable_spin_orbital),
        "index_convention": _INDEX_CONVENTION[bool(p.enable_spin_orbital)],
        "two_sz": p.two_sz, "ncond": p.ncond, "seed": spec.seed,
        "h_used": res.h_used, "seed_base": res.seed_base,
        "hwave_version": hwave.__version__,
        "input_sha256": dict(facts.input_sha256),
    }


def initial_state(input_dir_or_model, candidate_or_run, *, out, subshape=None,
                  site_map=None, uhfr_toml=None, para=None) -> pathlib.Path:
    """Seed one candidate and write its initial-state file (design section 9).

    Parameters
    ----------
    input_dir_or_model : path-like
        The UHFk input (a TOML, or a directory holding exactly one); for a
        UHFr target the UHFk input of the same lattice.
    candidate_or_run : Candidate or CandidateRun
        A :class:`Candidate` is evaluated on ``subshape`` (default its
        PatternSpec ``subshape``); a :class:`CandidateRun` on its own
        ``subshape`` (``subshape`` must then be ``None`` or equal, else
        ``ValueError``).
    out : path-like
        The file to write (its directory must exist). UHFk: the npz of
        design section 9.1, read as ``[file.input] initial``; UHFr: the
        ``Initial`` text file.
    site_map : sequence of ``(ix, iy, iz, a)``, optional
        UHFr only, with ``uhfr_toml``: for UHFr site ``i = 0 ... Nsite - 1``
        its UHFk cell in ``CellShape`` and physical orbital; a bijection
        onto the cells times the physical orbitals.
    uhfr_toml : path-like, optional
        UHFr only: the target UHFr input, whose ``Nsite``, ``Ncond`` and
        ``2Sz`` must agree with the map and the model.
    para : None
        Reserved for seeding on the para Hartree-Fock state (a later
        release); ``None`` seeds on the bare ``H0`` (``seed_base =
        "bare"``).

    Returns
    -------
    pathlib.Path
        ``out``.

    Raises
    ------
    InitialStateError
        The seed is refused (block check, no response, the checks of
        design section 9.2) or a spin-orbital model is written for UHFr.
    ValueError
        A conflicting ``subshape``, a pattern q the cell cannot hold, a
        bad ``site_map`` or a UHFr input that does not fit it.
    """
    from hwave.initial.facts import model_facts
    if (site_map is None) != (uhfr_toml is None):
        raise ValueError("site_map and uhfr_toml go together (a UHFr target needs both)")
    toml = _toml_of(input_dir_or_model)
    run = _run_of(candidate_or_run, subshape)
    facts = model_facts(toml, subshape=tuple(run.subshape))
    p = facts.problem
    if site_map is not None:
        if p.enable_spin_orbital:
            raise InitialStateError(SO_UHFR_REFUSAL)
        cells = _check_site_map(site_map, facts)
        _check_uhfr_input(uhfr_toml, facts, len(cells))
    res = _seed(facts, run, para)
    G = np.array(res.G, dtype=np.complex128)
    check_green(G, p, what="the seed of {}".format(run.run_id))
    out = _abs(out)
    if site_map is not None:
        write_uhfr_initial(uhfr_green(G, p, cells), out)
        return out
    members = {k: _encode(v) for k, v in _provenance(facts, run, res).items()}
    for k, v in CONVENTION_MEMBERS.items():
        members[k] = np.array(v)
    if int(np.prod(p.subshape)) > 1:
        members["green_sublattice"] = G
        members["green"] = initial_facade.unfold_layout(
            G, cellshape=p.cellshape, subshape=p.subshape, norb_orig=p.norb_orig, ns=p.ns)
    else:
        members["green"] = G
    return _write_npz(out, members)


# --------------------------------------------------------------------------
# UHFr
# --------------------------------------------------------------------------

def _check_site_map(site_map, facts):
    L = tuple(int(x) for x in facts.cellshape)
    nphys = int(facts.norb_phys)
    if isinstance(site_map, (str, bytes)) or not hasattr(site_map, "__len__"):
        raise ValueError("site_map: a sequence of (ix, iy, iz, a), one per UHFr site, got "
                         "{!r}".format(site_map))
    cells = []
    for i, entry in enumerate(site_map):
        if isinstance(entry, (str, bytes)) or not hasattr(entry, "__len__") \
                or len(entry) != 4:
            raise ValueError("site_map[{}]: (ix, iy, iz, a) expected, got {!r}"
                             .format(i, entry))
        vals = []
        for v in entry:
            if isinstance(v, bool) or not isinstance(v, numbers.Integral):
                raise ValueError("site_map[{}]: integers expected, got {!r}".format(i, entry))
            vals.append(int(v))
        ix, iy, iz, a = vals
        if not (0 <= ix < L[0] and 0 <= iy < L[1] and 0 <= iz < L[2]):
            raise ValueError("site_map[{}]: cell {} is outside CellShape {}".format(
                i, [ix, iy, iz], list(L)))
        if not 0 <= a < nphys:
            raise ValueError("site_map[{}]: orbital {} is not a physical orbital of a cell "
                             "with {} physical orbitals".format(i, a, nphys))
        cells.append((ix, iy, iz, a))
    n = L[0] * L[1] * L[2] * nphys
    if len(set(cells)) != len(cells) or len(cells) != n:
        dup = sorted({c for c in cells if cells.count(c) > 1})
        raise ValueError(
            "site_map must be a bijection onto the {} cells x {} physical orbitals ({} "
            "entries); got {} entries{}".format(
                L[0] * L[1] * L[2], nphys, n, len(cells),
                ", repeated: {}".format(dup[:4]) if dup else ""))
    return tuple(cells)


def _check_uhfr_input(uhfr_toml, facts, nsite):
    doc = _read_doc(uhfr_toml)
    mode = doc.get("mode", {})
    if not isinstance(mode, dict) or mode.get("mode") != "UHFr":
        raise ValueError("{}: uhfr_toml must be a UHFr input ([mode] mode = \"UHFr\")"
                         .format(uhfr_toml))
    # [mode.param] keys are case-insensitive (ParamDict)
    param = {str(k).lower(): v for k, v in (mode.get("param", {}) or {}).items()}
    if param.get("nsite") != nsite:
        raise ValueError("{}: Nsite = {!r} but the site_map has {} sites".format(
            uhfr_toml, param.get("nsite"), nsite))
    if param.get("2sz") != facts.two_sz:
        raise ValueError("{}: 2Sz = {!r} but the model has 2Sz = {!r}".format(
            uhfr_toml, param.get("2sz"), facts.two_sz))
    if "ncond" in param and param["ncond"] != facts.ncond:
        raise ValueError("{}: Ncond = {!r} but the model has Ncond = {!r}".format(
            uhfr_toml, param["ncond"], facts.ncond))


def uhfr_green(G, problem, cells) -> np.ndarray:
    """The UHFr Green function ``(2 Nsite, 2 Nsite)``, index ``i + s *
    Nsite``, of the stored UHFk ``G`` on the cell of ``problem`` for the
    UHFr sites ``cells[i] = (ix, iy, iz, a)``:
    ``G_UHFr[(i, s), (j, t)] = G[r(R_B,i - R_B,j), s, a'_i, t, a'_j]``."""
    p = problem
    B = tuple(p.subshape)
    N = tuple(p.shape)
    cells = np.asarray(cells, dtype=int).reshape(-1, 4)
    R = cells[:, :3]
    RB = R // np.array(B)
    off = R % np.array(B)
    ir = off[:, 0] + B[0] * (off[:, 1] + B[1] * off[:, 2])
    aa = cells[:, 3] + int(p.norb_orig) * ir
    d = (RB[:, None, :] - RB[None, :, :]) % np.array(N)
    r = d[..., 2] + N[2] * (d[..., 1] + N[1] * d[..., 0])
    G = np.asarray(G)
    n = len(cells)
    out = np.zeros((2 * n, 2 * n), dtype=np.complex128)
    for s in range(2):
        for t in range(2):
            out[s * n:(s + 1) * n, t * n:(t + 1) * n] = G[r, s, aa[:, None], t, aa[None, :]]
    return out


def write_uhfr_initial(G_2N, path) -> pathlib.Path:
    """Write ``G_2N`` ``(2N, 2N)`` (index ``i + s * N``) in the UHFr
    ``Initial`` format (the five-line header, then ``i s j t re im`` for
    every nonzero entry, row-major; the layout ``UHFr.save_results``
    writes)."""
    G = np.asarray(G_2N)
    if G.ndim != 2 or G.shape[0] != G.shape[1] or G.shape[0] % 2:
        raise ValueError("G_2N must be (2N, 2N), got shape {}".format(G.shape))
    n = G.shape[0] // 2
    rows = []
    for i1 in range(2 * n):
        for i2 in range(2 * n):
            v = complex(G[i1, i2])
            if abs(v) > 0:
                rows.append("{} {} {} {} {} {}\n".format(
                    i1 % n, i1 // n, i2 % n, i2 // n, v.real, v.imag))
    path = _abs(path)
    with open(path, "w") as f:
        f.write("===============================\n")
        f.write("NCisAjs {}\n".format(len(rows)))
        f.write("===============================\n" * 3)
        f.writelines(rows)
    return path


# --------------------------------------------------------------------------
# the check
# --------------------------------------------------------------------------

def _load_members(path):
    """Every member of the npz at ``path`` (``allow_pickle=False``); an
    object member is refused, naming it."""
    out = {}
    try:
        data = np.load(str(path), allow_pickle=False)
    except (OSError, ValueError, zipfile.BadZipFile) as exc:
        raise InitialStateError("{}: not a readable npz file: {}".format(path, exc)) from None
    with data:
        if not hasattr(data, "files"):
            raise InitialStateError("{}: not an npz archive".format(path))
        for name in data.files:
            try:
                v = data[name]
            except ValueError as exc:
                # an object array needs pickle, which is never used
                raise InitialStateError(
                    "{}: member {!r} is an object array (pickled), which is never read: "
                    "{}".format(path, name, exc)) from None
            out[name] = v
    return out


def _decode(members, key, path):
    if key not in members:
        raise InitialStateError("{}: provenance member {!r} is missing".format(path, key))
    v = members[key]
    if v.shape != () or v.dtype.kind != "U":
        raise InitialStateError(
            "{}: provenance member {!r} must be a 0-d unicode string, got dtype {} shape "
            "{}".format(path, key, v.dtype, v.shape))
    try:
        return json.loads(str(v[()]))
    except ValueError as exc:
        raise InitialStateError("{}: provenance member {!r} is not JSON: {}".format(
            path, key, exc)) from None


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_provenance(members, problem, config_path, path):
    p = problem
    prov = {k: _decode(members, k, path) for k in PROVENANCE_KEYS}
    if prov["initial_state_schema"] != INITIAL_STATE_SCHEMA:
        raise InitialStateError("{}: initial_state_schema {!r} is not {} (this version "
                                "reads schema {})".format(path, prov["initial_state_schema"],
                                                          INITIAL_STATE_SCHEMA,
                                                          INITIAL_STATE_SCHEMA))
    want = {
        "subshape": list(p.subshape), "cellshape": list(p.cellshape),
        "norb_orig": int(p.norb_orig), "spin_orbital": bool(p.enable_spin_orbital),
        "index_convention": _INDEX_CONVENTION[bool(p.enable_spin_orbital)],
        "two_sz": p.two_sz,
    }
    for key in ("spin_orbital", "index_convention", "subshape", "cellshape", "norb_orig",
                "two_sz"):
        if prov[key] != want[key] or type(prov[key]) is not type(want[key]):
            raise InitialStateError(
                "{}: {}: the file was written for {} {}, the input {} has {}".format(
                    path, key, key, json.dumps(prov[key]), config_path,
                    json.dumps(want[key])))
    n = prov["ncond"]
    if isinstance(n, bool) or not isinstance(n, numbers.Real) \
            or abs(float(n) - float(p.ncond)) > 1e-12 * max(1.0, abs(float(p.ncond))):
        raise InitialStateError("{}: ncond: the file was written for Ncond {}, the input {} "
                                "has {}".format(path, json.dumps(n), config_path, p.ncond))
    for key, value in CONVENTION_MEMBERS.items():
        got = members.get(key)
        if got is None or got.shape != () or str(got[()]) != value:
            raise InitialStateError("{}: {}: {!r} expected, got {}".format(
                path, key, value, None if got is None else repr(str(got[()]))))
    warnings = []
    sha = prov["input_sha256"]
    if isinstance(sha, dict):
        now = {}
        for key, f in initial_facade.referenced_input_files(config_path).items():
            now[key] = _sha256(f) if f.is_file() else None
        changed = sorted(k for k in set(now) | set(sha) if now.get(k) != sha.get(k))
        if changed:
            warnings.append("input_sha256: the input files {} differ from those the seed "
                            "was written for".format(", ".join(changed)))
    return warnings


def check_initial_state(config_path, npz_path=None, *, require_provenance=False
                        ) -> CheckResult:
    """Check an initial-state npz against the UHFk input that will read it
    (design section 9.3).

    ``npz_path=None`` checks the file the solver will read
    (:func:`solver_initial_path`). With provenance the members
    ``subshape``, ``cellshape``, ``norb_orig``, ``spin_orbital``,
    ``index_convention``, ``two_sz``, ``ncond`` and the convention markers
    must match the input, and the member the solver reads must be present
    (``green_sublattice`` when ``prod(SubShape) > 1``, else ``green``).
    Every file (with or without provenance) is checked for that member, its
    shape against ``SubShape``, ``norb`` and the spin-orbital layout,
    Hermiticity, the eigenvalue bounds, the trace, the 2Sz sectors and the
    solver blocks (:func:`check_green`). A file without provenance gives
    ``verified = False`` and :data:`NO_PROVENANCE_WARNING`, or with
    ``require_provenance`` a :class:`ProvenanceMissingError`.

    Raises :class:`InitialStateError` on any refusal, naming the field or
    the check (also for an unreadable input, ``ValueError`` from
    :func:`initial_facade.load_problem`).
    """
    config_path = _abs(config_path)
    try:
        problem = initial_facade.load_problem(config_path)
    except ValueError as exc:
        if isinstance(exc, InitialStateError):
            raise
        raise InitialStateError(str(exc)) from None
    finput = _finput(_read_doc(config_path))
    if "initial" in finput and "initial_uhf" in finput:
        raise InitialStateError("{}: [file.input] sets both initial and initial_uhf, which "
                                "the solver refuses; remove initial_uhf".format(config_path))
    if npz_path is None:
        npz_path = solver_initial_path(config_path)
        if npz_path is None:
            raise InitialStateError("{}: [file.input] initial is not set, so the solver "
                                    "reads no initial-state file".format(config_path))
        if not os.path.isfile(npz_path):
            pti = _path_to_input(finput)
            raise InitialStateError(
                "{}: the solver reads [file.input] initial from {} (path_to_input = {} "
                "joined with initial = {}, from the TOML's directory), which resolves to "
                "{}: no such file".format(
                    config_path, os.path.join(pti, finput["initial"]), json.dumps(pti),
                    json.dumps(finput["initial"]), os.path.abspath(npz_path)))
    path = _abs(npz_path)
    if not path.is_file():
        raise InitialStateError("{}: no such file".format(path))
    members = _load_members(path)
    sub = int(np.prod(problem.subshape)) > 1
    has_prov = "initial_state_schema" in members
    warnings = []
    if has_prov:
        warnings += _check_provenance(members, problem, config_path, path)
    if sub and "green_sublattice" in members:
        key, G = "green_sublattice", members["green_sublattice"]
    elif sub and "green" in members and not has_prov:
        # the solver folds 'green' (lossy for a translation-breaking state)
        key = "green"
        try:
            G = initial_facade.fold_layout(
                members["green"], cellshape=problem.cellshape, subshape=problem.subshape,
                norb_orig=problem.norb_orig, ns=problem.ns)
        except ValueError as exc:
            raise InitialStateError("{}: green: shape: {}".format(path, exc)) from None
        if "green_convention" not in members:
            warnings.append("green: no green_convention marker; the solver folds it with "
                            "the offset on the first slot")
    elif not sub and "green" in members:
        key, G = "green", members["green"]
    elif sub and "green" in members:
        raise InitialStateError(
            "{}: green_sublattice is missing: a file with provenance always has it for "
            "SubShape {}; the solver would fold the lossy 'green' view instead".format(
                path, list(problem.subshape)))
    else:
        want = "green_sublattice" if sub else "green"
        raise InitialStateError(
            "{}: {} is missing: the solver reads {} for SubShape {} and, without it, "
            "would start from a zero Green function".format(
                path, want, " or ".join(("green_sublattice", "green")) if sub else want,
                list(problem.subshape)))
    check_green(G, problem, what="{}: {}".format(path, key))
    if not has_prov:
        if require_provenance:
            raise ProvenanceMissingError("{}: {}".format(path, NO_PROVENANCE_WARNING))
        warnings.insert(0, NO_PROVENANCE_WARNING)
    return CheckResult(verified=has_prov, key=key, warnings=tuple(warnings))
