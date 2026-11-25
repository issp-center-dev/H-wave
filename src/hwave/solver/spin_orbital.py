"""Spin-orbital (``[mode] enable_spin_orbital``) conventions shared by the
solvers (GitHub issue #83, step 2 of the spin-orbital roadmap).

In spin-orbital mode the spin is folded into the orbital index:

* the Geometry ``norb`` and the Transfer file count SPIN-ORBITALS, in the
  INTERLEAVED order ``i = 2*orb + spin`` (Wannier90 SOI, UHFk's native
  order);
* the two-body interaction files are written in PHYSICAL orbital indices
  ``0 <= orb < norb // 2``;
* RPA and FLEX work internally, and write their susceptibilities, in the
  SPIN-BLOCK order ``j = spin*norb_phys + orb`` (tagged
  ``index_convention = "spin_block"``); UHFk evaluates its interaction
  terms in the same spin-block ("virtual spin") frame.

Everything that converts between these orders, counts physical orbitals or
reads the switch lives here, so every solver agrees on it.
"""
import logging
from dataclasses import dataclass

import numpy as np
from requests.structures import CaseInsensitiveDict

from . import fingerprint, so_records
from .backend import strict_bool

#: The ``index_convention`` tag RPA/FLEX write into chi0q/chiq files.
INDEX_CONVENTION = "spin_block"

_FLAG = "enable_spin_orbital"

logger = logging.getLogger(__name__)


def resolve_flag(info_mode, name="[mode] enable_spin_orbital", *, allow_float=False):
    """``enable_spin_orbital`` of a ``[mode]`` dict as a real boolean.

    The key lookup is case-insensitive and the value is parsed strictly
    (:func:`hwave.solver.backend.strict_bool`): a misspelled value raises
    instead of silently reading as true or false. A recognized string is
    accepted with a warning. ``allow_float`` additionally accepts an exact
    float 0.0/1.0 -- for UHFk only, whose own parser accepted them before
    the switch was shared; every other caller rejects floats."""
    if info_mode is None:
        return False
    value = CaseInsensitiveDict(info_mode).get(_FLAG, False)
    if allow_float and isinstance(value, float) and value in (0.0, 1.0):
        return bool(value)        # accepted by UHFk before the switch was shared
    if isinstance(value, str):
        logger.warning(
            "%s should be a boolean, not the string %r; write "
            "enable_spin_orbital = true/false without quotes.", name, value)
    return strict_bool(value, name)


def physical_norb(geom_norb, enable_spin_orbital, *, check_norb=None,
                  source="geom.dat"):
    """Physical orbital count from a geometry ``norb``.

    In spin-orbital mode ``geom.dat``'s ``norb`` is the spin-orbital count
    (= 2 * physical orbitals = Wannier90 num_wann), so halve it. Evenness is
    validated on ``check_norb`` (the *original*, pre-sublattice-fold value)
    so the error names the user's actual ``geom.dat`` entry, while the
    returned count is derived from ``geom_norb`` (which may be the post-fold
    value ``orig * subvol``)."""
    if not enable_spin_orbital:
        return geom_norb
    cn = check_norb if check_norb is not None else geom_norb
    if cn % 2 != 0:
        raise ValueError(
            "spin-orbital mode requires an even Geometry norb (the spin-orbital "
            "count = 2 * physical orbitals); got {} in {}".format(cn, source))
    return geom_norb // 2


def _extern_range_sentence(norb_phys):
    return ("In spin-orbital mode Extern uses PHYSICAL orbital indices (norb_phys = "
            "{n}), while the Transfer file uses the 2\u00b7norb_phys = {n2} "
            "spin-orbital indices; write the field per physical orbital (the sign "
            "per spin comes from \u03c3^z).".format(n=norb_phys, n2=2 * norb_phys))


def _extern_index_error(file_name, index_1based, norb_phys):
    return ValueError(
        "Extern file '{f}': orbital index {i} is outside 1..{n}. {rest}".format(
            f=file_name, i=index_1based, n=norb_phys,
            rest=_extern_range_sentence(norb_phys)))


def _extern_header_error(file_name, declared_norb, norb_phys):
    return ValueError(
        "Extern file '{f}': its header declares {d} orbitals, but norb_phys = "
        "{n}. {rest}".format(f=file_name, d=declared_norb, n=norb_phys,
                             rest=_extern_range_sentence(norb_phys)))


def check_extern_indices(table, norb_phys, file_name, declared_norb=None):
    """Refuse a spin-orbital ``Extern`` table outside the physical orbitals.

    In spin-orbital mode the Extern operator is
    ``H = sum H_ab(R) sigma^z_{ss'} c+_{a s} c_{b s'}``: the file is written
    in PHYSICAL orbital indices ``1..norb_phys`` (the spin sign comes from
    ``sigma^z``), not in the ``1..2*norb_phys`` spin-orbital indices of the
    Transfer file. An index outside the physical range used to be dropped
    without a word; it is refused here.

    Parameters
    ----------
    table : dict
        ``{(irvec, (a, b)): value}`` as read by
        :func:`hwave.qlmsio.wan90.read_w90` (0-based ``a``, ``b``), before
        any sublattice folding.
    norb_phys : int
        Physical orbital count (Geometry ``norb // 2``), per original cell.
    file_name : str
        Source path, named in the error.
    declared_norb : int, optional
        The orbital count the file header declares (line 2). Checked first:
        a header other than ``norb_phys`` is refused even when every row is
        in range. ``None`` (a caller without a file) skips it.

    Raises
    ------
    ValueError
        Naming the declared header dimension, or else the first offending
        index (1-based), and both index ranges.
    """
    if declared_norb is not None and int(declared_norb) != norb_phys:
        raise _extern_header_error(file_name, int(declared_norb), norb_phys)
    for (_irvec, orbvec), _v in table.items():
        for i in orbvec:
            if not 0 <= i < norb_phys:
                raise _extern_index_error(file_name, i + 1, norb_phys)


def validate_coeff_extern(raw):
    """``[mode.param] coeff_extern`` of a spin-orbital run as a float: a
    finite real number (a bool is refused: ``float(True) == 1.0`` would
    silently switch the field on). Shared by RPA, FLEX and ``hwave_sc``, so
    one value is accepted or refused alike."""
    import numbers
    if (isinstance(raw, (bool, np.bool_)) or not isinstance(raw, numbers.Real)
            or not np.isfinite(float(raw))):
        raise ValueError(
            "[mode.param] coeff_extern must be a finite real number, got "
            "{!r}".format(raw))
    return float(raw)


def interleaved_index(orb, spin):
    """Interleaved spin-orbital index of ``(orb, spin)``."""
    return 2 * orb + spin


def spin_of(index):
    """Spin (0 = up, 1 = down) of an interleaved index (array-valued too)."""
    return index % 2


def to_spinblock_index(index, norb_phys):
    """Spin-block position of interleaved index ``index``."""
    return (index % 2) * norb_phys + index // 2


def spinblock_index(orb, spin, norb_phys):
    """Spin-block composite index ``spin * norb_phys + orb`` of physical
    orbital ``orb`` with ``spin``: the flattened ``(s, a)`` index of the
    normal-mode ``(nvol, ns, norb, ns, norb)`` arrays, and the key the
    spin-resolved 9-column Transfer reader folds to (issue #323)."""
    return spin * norb_phys + orb


def spinblock_split(index, norb_phys):
    """``(orb, spin)`` of a spin-block composite index (array-valued too)."""
    return index % norb_phys, index // norb_phys


def spin_resolved_transfer_format(param_ham):
    """``(format, version)`` when the Transfer table came from a spin-resolved
    9-column file (``ham_param.source_formats``, recorded by the k-space
    input reader), else ``None``."""
    formats = getattr(param_ham, "source_formats", None) or {}
    return formats.get("transfer")


def _spin_resolved_name_tag(param_ham, fmt):
    files = getattr(param_ham, "source_files", None) or {}
    return files.get("transfer", "Transfer"), "{}/{}".format(*fmt)


def check_spin_resolved_clash(param_ham, enable_spin_orbital):
    """Refuse a spin-resolved (9-column) Transfer together with spin-orbital
    mode: the format carries its own spin columns and physical orbital
    indices (Geometry ``norb`` = physical count), while spin-orbital mode
    folds the spin into interleaved indices and reads the Geometry ``norb``
    as the spin-orbital count. A solver that consumes the spin-resolved
    table (UHFk, issue #323) calls this alone; one that does not yet calls
    :func:`reject_spin_resolved_transfer`, which includes it."""
    fmt = spin_resolved_transfer_format(param_ham)
    if fmt is None or not enable_spin_orbital:
        return
    name, tag = _spin_resolved_name_tag(param_ham, fmt)
    raise ValueError(
        "{}: the spin-resolved Transfer format ({}) carries its own spin "
        "columns and physical orbital indices, which contradicts "
        "[mode] enable_spin_orbital = true (interleaved spin-orbital "
        "indices, Geometry norb = 2 * physical). Use one or the other."
        .format(name, tag))


def check_spin_resolved_norb(param_ham, norb_geom, solver):
    """Refuse a spin-resolved (9-column) Transfer whose line-2 orbital count
    differs from the Geometry ``norb``.

    The reader folds every row to ``s * N + a`` with ``N`` from line 2 of
    the Transfer file, and the solver decodes that key with the Geometry
    ``norb``. When the two differ the spin and orbital indices are mixed up
    without any error (line 2 = 1, Geometry ``norb`` = 2 puts the down spin
    of orbital 1 on the up spin of orbital 2), so the two must agree. The
    7-column format is not affected: its line 2 is not used for indexing.

    Parameters
    ----------
    param_ham : Mapping
        The k-space reader's table container (``source_formats`` /
        ``source_dimensions`` / ``source_files``).
    norb_geom : int
        The Geometry ``norb`` (physical orbitals per unit cell, before any
        sublattice fold).
    solver : str
        Named in the error.

    Raises
    ------
    ValueError
    """
    fmt = spin_resolved_transfer_format(param_ham)
    if fmt is None:
        return
    dims = getattr(param_ham, "source_dimensions", None) or {}
    n_file = dims.get("transfer")
    if n_file is not None and int(n_file) == int(norb_geom):
        return
    name, tag = _spin_resolved_name_tag(param_ham, fmt)
    raise ValueError(
        "{}: line 2 of this spin-resolved Transfer file ({}) gives norb = "
        "{}, but the Geometry norb is {} ({}). In this format line 2 "
        "is the number of PHYSICAL orbitals and must equal the Geometry "
        "norb (the spin is a separate column); the rows are indexed with "
        "it, so a different value would mix up spin and orbital."
        .format(name, tag, n_file, norb_geom, solver))


def spin_free_transfer_table(table, norb_phys):
    """The 7-column table of a spin-resolved Transfer that does not depend
    on spin, else ``None``.

    ``table`` is keyed on the spin-block composite index ``s*norb_phys + a``
    (the reader's fold of the 9-column rows). It is spin-independent when
    no entry couples the two spins (a zero-valued one is ignored) and the
    up-up and down-down blocks are identical, value for value; the result
    is then the up-up block on physical indices -- exactly the table the
    equivalent 7-column file reads to. Any spin dependence (a Zeeman-like
    asymmetry or a spin flip), or an index outside ``[0, 2*norb_phys)``,
    gives ``None``.
    """
    n = int(norb_phys)
    blocks = ({}, {})
    for (irvec, (i, j)), v in table.items():
        if not (0 <= i < 2 * n and 0 <= j < 2 * n):
            return None
        si, sj = i // n, j // n
        if si != sj:
            if v != 0:
                return None
            continue
        blocks[si][(irvec, (i - si * n, j - sj * n))] = v
    up, dn = blocks
    for key in set(up) | set(dn):
        if up.get(key, 0) != dn.get(key, 0):
            return None
    # an entry present in the down block only is zero (checked above)
    return dict(up)


def validate_spin_resolved_closure(table, norb_phys, source="Transfer"):
    """Refuse a spin-resolved Transfer that is not Hermitian-closed.

    In the 9-column format the partner of ``T_{as,bt}(R)`` is
    ``T_{bt,as}(-R)`` and the two must be complex conjugates; a declared
    entry without its partner (a spin-flip row written in one direction
    only, say) is refused, naming both entries as (orbital, spin). The
    same rule and tolerance as the two-body declarations
    (:func:`hwave.solver.declarations.validate_hermitian_closure`).
    ``table`` is keyed on the spin-block index ``s*norb_phys + a``.
    """
    from .declarations import REVERSE_ATOL, REVERSE_RTOL
    n = int(norb_phys)

    def fmt(irvec, i, j):
        (a, s), (b, t) = spinblock_split(i, n), spinblock_split(j, n)
        return "R={} (orbital, spin) = ({}, {}), ({}, {})".format(
            tuple(irvec), a + 1, s, b + 1, t)

    for (irvec, (i, j)), v in table.items():
        if not np.isfinite(complex(v)):
            raise ValueError(
                "{}: the spin-resolved Transfer declares a non-finite "
                "coefficient at {}: {}".format(source, fmt(irvec, i, j), v))
        neg = tuple(-x for x in irvec)
        partner = table.get((neg, (j, i)))
        if partner is None:
            raise ValueError(
                "{}: the spin-resolved Transfer declares {} = {} with NO "
                "partner entry at {}: declare both directions "
                "(T_as,bt(R) = conj T_bt,as(-R))".format(
                    source, fmt(irvec, i, j), v, fmt(neg, j, i)))
        if abs(v - np.conj(partner)) > REVERSE_ATOL + REVERSE_RTOL * abs(v):
            raise ValueError(
                "{}: the spin-resolved Transfer is not Hermitian-closed: "
                "{} = {} but {} = {} (expected its conjugate, {})".format(
                    source, fmt(irvec, i, j), v, fmt(neg, j, i), partner,
                    np.conj(v)))


def private_table_container(param_ham):
    """A copy of the k-space reader's table container that a solver can
    rewrite without touching the caller's: the same tables (shared, not
    copied) and the same record attributes (``source_files``,
    ``source_dimensions``, ``source_formats``; their dicts copied)."""
    new = type(param_ham)()
    for key, value in param_ham.items():
        new[key] = value
    for attr, value in vars(param_ham).items():
        if attr == "_store":
            continue
        setattr(new, attr, dict(value) if isinstance(value, dict) else value)
    return new


def spinful_one_body(ham_info):
    """Whether the one-body term of ``ham_info`` carries a spin index (is
    held as nd x nd): spin-orbital mode or a spin-resolved Transfer.

    The one accessor every consumer uses, so an Interaction-like object
    that predates ``spinful_one_body`` (a test stub, an external caller)
    falls back to its ``enable_spin_orbital`` instead of silently losing
    the spin-orbital gates or raising AttributeError."""
    value = getattr(ham_info, "spinful_one_body", None)
    if value is not None:
        return bool(value)
    return bool(getattr(ham_info, "enable_spin_orbital", False))


def spin_resolved_one_body(ham_info):
    """Whether ``ham_info``'s Transfer came from a spin-resolved (9-column)
    file; False for an Interaction-like object that predates the
    attribute."""
    return bool(getattr(ham_info, "spin_resolved_transfer", False))


def spinful_label(ham_info):
    """How the spinful one-body term was declared, for messages: the
    spin-orbital switch or a spin-resolved (9-column) Transfer (#323)."""
    if getattr(ham_info, "enable_spin_orbital", False):
        return "[mode] enable_spin_orbital = true"
    return "a spin-resolved (9-column) Transfer"


def spinful_term(ham_info):
    """Short form of :func:`spinful_label` for "... does not support X"
    messages; the spin-orbital wording is the historical one the
    equivalence-table reject fragments pin."""
    if getattr(ham_info, "enable_spin_orbital", False):
        return "enable_spin_orbital"
    return "a spin-resolved (9-column) Transfer"


def reject_spin_resolved_transfer(param_ham, solver, enable_spin_orbital):
    """Refuse a spin-resolved (9-column) Transfer table a solver does not
    consume yet, and -- always -- its combination with spin-orbital mode.

    The spin-resolved format carries explicit spin columns and PHYSICAL
    orbital indices (Geometry ``norb`` = physical count); spin-orbital mode
    folds the spin into interleaved indices and reads the Geometry ``norb``
    as the spin-orbital count. The two disagree on what ``norb`` means, so
    they are never combined. Until a solver consumes the spin-resolved
    table (issue #323: UHFk with the H(k) builder switch, RPA / FLEX and
    hwave_sc afterwards) it refuses the table here, right after reading,
    rather than dropping half of it in ``_make_ham_trans``.

    Parameters
    ----------
    param_ham : Mapping
        The k-space reader's table container (``source_formats`` /
        ``source_files`` attributes); a plain dict is treated as 7-column.
    solver : str
        Named in the error.
    enable_spin_orbital : bool
        The resolved ``[mode] enable_spin_orbital``.

    Raises
    ------
    ValueError
    """
    fmt = spin_resolved_transfer_format(param_ham)
    if fmt is None:
        return
    check_spin_resolved_clash(param_ham, enable_spin_orbital)
    name, tag = _spin_resolved_name_tag(param_ham, fmt)
    raise ValueError(
        "{}: the spin-resolved Transfer format ({}) is not yet consumed by "
        "the {} solver. Write the Transfer in the 7-column "
        "form, or fold the spin into the orbital index and set "
        "[mode] enable_spin_orbital = true.".format(name, tag, solver))


def spinblock_order(norb_phys):
    """Gather order ``g`` with ``g[j] =`` the interleaved index of spin-block
    position ``j``: ``M_spinblock = M_interleaved[g][:, g]``."""
    j = np.arange(2 * norb_phys)
    return interleaved_index(j % norb_phys, j // norb_phys)


def interleaved_to_spinblock(mat, norb_phys, axes=(-2, -1)):
    """Reorder the given axes of ``mat`` from interleaved to spin-block
    order (a copy; the identity for ``norb_phys == 1``)."""
    g = spinblock_order(norb_phys)
    out = np.asarray(mat)
    for ax in axes:
        out = np.take(out, g, axis=ax)
    return out


def spinblock_to_interleaved(mat, norb_phys, axes=(-2, -1)):
    """Inverse of :func:`interleaved_to_spinblock`."""
    inv = np.argsort(spinblock_order(norb_phys))
    out = np.asarray(mat)
    for ax in axes:
        out = np.take(out, inv, axis=ax)
    return out


def validate_chi0q_index_convention(data, enable_spin_orbital, file_name=""):
    """Reject a stored chi0q whose spin-orbital index convention is unknown.

    RPA writes ``index_convention="spin_block"`` (spin*norb+orb) into its
    chi0q.npz / chiq.npz. UHFk instead uses the interleaved (2*orb+spin)
    ordering, and chi0q files produced before the SO convention fix
    (commit 9dd9a21) carry no ``index_convention`` key at all. In spin-orbital
    mode the two orderings differ, so silently accepting such a file would mix
    conventions and corrupt the result. Outside spin-orbital mode the
    convention is irrelevant and this is a no-op.

    Parameters
    ----------
    data : Mapping
        Loaded npz (e.g. ``numpy.load(...)``); must support ``in`` and ``[]``.
    enable_spin_orbital : bool
        Whether the consuming calculation runs in spin-orbital mode.
    file_name : str, optional
        Source path, used only for the error message.

    Raises
    ------
    ValueError
        If spin-orbital mode is on and the file lacks a ``spin_block``
        ``index_convention`` marker.
    """
    if not enable_spin_orbital:
        return
    if "index_convention" not in data:
        raise ValueError(
            "chi0q file '{}' lacks an 'index_convention' marker. It predates "
            "the spin-orbital index convention fix and may use the interleaved "
            "(2*orb+spin) ordering; regenerate it with the current RPA solver "
            "(which writes index_convention='spin_block').".format(file_name)
        )
    conv = str(data["index_convention"])
    if conv != INDEX_CONVENTION:
        raise ValueError(
            "chi0q file '{}' has index_convention='{}', but spin-orbital mode "
            "requires 'spin_block' (spin*norb+orb). Regenerate it with the "
            "current RPA solver.".format(file_name, conv)
        )


#: The WARNING a spin-orbital entry logs, once per public run, when
#: ``[mode.param] orbital_site`` is given.
ORBITAL_SITE_WARNING = ("orbital_site is ignored in spin-orbital mode; all R=0 rows are "
                        "treated exactly (see migration note)")

#: The ``[mode]`` dict key under which a public entry hands its
#: :class:`SOEntryContext` to a nested solver.
SO_ENTRY_CONTEXT_KEY = "_so_entry_context"

#: lower-cased table names handed to the interaction compiler (the
#: interaction types and the aggregate ``Coulomb``); other tables of the
#: ``param_ham`` dict (Geometry, Transfer, Extern, Initial, ...) are not
#: interactions
_INTERACTION_KEYS = frozenset(t.lower() for t in so_records.TYPES + ("Coulomb",))


@dataclass
class SOEntryContext:
    """What a public spin-orbital entry establishes once for the whole run.

    Attributes
    ----------
    compiled : so_records.CompiledInteraction
        The interaction compiled from the pre-fold tables.
    orbital_site_warned : bool
        True when ``[mode.param] orbital_site`` was given (and the WARNING
        :data:`ORBITAL_SITE_WARNING` has been logged).
    interaction : fingerprint.InteractionInputs
        The fingerprints of ``compiled`` (``fp_interaction`` over its
        records, ``fp_topology`` over its declared topology), which the
        archive writers record and the archive consumers compare.
    tables_digest : str
        :func:`interaction_tables_digest` of the raw tables ``compiled``
        was compiled from. A nested solver that builds its vertex from its
        own reader tables compares it, so the fingerprints it records
        describe the tables its vertex is built from.
    """

    compiled: so_records.CompiledInteraction
    orbital_site_warned: bool
    interaction: fingerprint.InteractionInputs
    tables_digest: str


#: the refusal of a nested solver whose interaction tables are not the ones
#: its supplied entry context was compiled from
TABLES_MISMATCH_TEXT = (
    "enable_spin_orbital: the interaction tables of this solver differ from "
    "the ones its entry context was compiled from, so the recorded "
    "interaction fingerprint would not describe the vertex. Build the "
    "context and the solver from the same interaction files.")


def interaction_tables_digest(tables):
    """Digest of the raw interaction tables of a spin-orbital run.

    Only the interaction types and the aggregate ``Coulomb`` table are
    covered (keys matched case-insensitively). ``Coulomb`` is split as
    every reader splits it (``wan90.split_coulomb``: the ``R = 0`` diagonal
    goes to CoulombIntra, the rest to CoulombInter) and empty tables are
    left out, so the RPA reader (which keeps the aggregate) and the
    hwave_sc reader (which splits it) of the same files give one digest.
    """
    from hwave.qlmsio import wan90
    canonical = {t.lower(): t for t in so_records.TYPES}
    merged = {}

    def add(name, tbl):
        merged.setdefault(name, {}).update(tbl)

    for key, tbl in tables.items():
        low = str(key).lower()
        if low == "coulomb":
            intra, inter = wan90.split_coulomb(tbl)
            add("CoulombIntra", intra)
            add("CoulombInter", inter)
        elif low in canonical:
            add(canonical[low], tbl)

    def chunks():
        for name in sorted(n for n, t in merged.items() if t):
            yield fingerprint.encode_string(name)
            yield fingerprint.encode_string(fingerprint.table_digest(name, merged[name]))

    return fingerprint.digest("interaction_tables", chunks())


def build_so_entry_context(param, tables, norb_phys, log):
    """The :class:`SOEntryContext` of one public spin-orbital entry.

    Parameters
    ----------
    param : mapping or None
        The ``[mode.param]`` dict. Its ``orbital_site`` key is looked up
        case-insensitively; the value is never parsed or validated, since
        spin-orbital mode treats every ``R = 0`` row exactly. When the key
        is present, :data:`ORBITAL_SITE_WARNING` is logged once on ``log``.
    tables : mapping
        The pre-fold interaction tables ``{type: {((Rx, Ry, Rz), (a, b)):
        v}}`` in physical orbital indices (0-based). Only the interaction
        types and ``Coulomb`` (case-insensitive) are compiled; any other
        table is left out.
    norb_phys : int
        Physical orbitals per cell.
    log : logging.Logger
        The entry's logger.

    Returns
    -------
    SOEntryContext

    Raises
    ------
    ValueError
        Every refusal of :func:`so_records.compile_interaction` (an
        imaginary part on a CoulombInter, Hund, Ising, Exchange or PairLift
        entry, a table that is not Hermitian-closed, ...).
    """
    # A CaseInsensitiveDict (ParamDict included) is asked directly: copying
    # a ParamDict reads every key, so the hwave entry would no longer
    # report the [mode.param] keys the run never read.
    if param is not None and not isinstance(param, CaseInsensitiveDict):
        param = CaseInsensitiveDict(param)
    warned = param is not None and "orbital_site" in param
    if warned:
        log.warning(ORBITAL_SITE_WARNING)
    # The readers own the interaction-key validation (the RPA reader refuses
    # an unknown key; hwave_sc reads a fixed list of types), so the
    # compiler's unknown-type refusal is not reached from this entry.
    interaction = {key: tbl for key, tbl in tables.items()
                   if str(key).lower() in _INTERACTION_KEYS}
    compiled = so_records.compile_interaction(interaction, norb_phys)
    return SOEntryContext(compiled=compiled, orbital_site_warned=warned,
                          interaction=fingerprint.interaction_inputs(compiled),
                          tables_digest=interaction_tables_digest(interaction))


def context_for_nested(info_mode, context):
    """Store ``context`` in ``info_mode[SO_ENTRY_CONTEXT_KEY]`` (a top-level
    key, never under ``param``), so a nested solver built from this
    ``[mode]`` dict reuses it: the interaction is not compiled again and
    the ``orbital_site`` WARNING is not logged again."""
    info_mode[SO_ENTRY_CONTEXT_KEY] = context


def offsite_types(compiled):
    """The sorted names of the types with a nonzero off-site record (its
    four cells not all equal) whose off-site rows enter the spin-orbital
    RPA vertex through their direct part.

    PairHop is never named: the RPA vertex drops an off-site PairHop row
    (with its own WARNING). An ``R = 0`` row between different orbitals is
    not off-site."""
    names = set()
    for r in compiled.records:
        cells = r["cells"]
        if r["type"] != "PairHop" and not (cells == cells[0]).all():
            names.add(str(r["type"]))
    return tuple(sorted(names))
