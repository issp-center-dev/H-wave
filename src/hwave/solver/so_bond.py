"""Spin-orbital bond-resolved two-body objects: the slot convention and the
slot geometry, the bond topology of a compiled interaction, the exact S_z
layout selection (:func:`select_layout`), the canonical per-spin-sector band
eigen-decomposition and the borrowed-Green-function projection
(:func:`canonical_band_eigen`, :func:`project_borrowed_green`), the S_z
sector map and the packed sector blocks (:func:`sz_sector_map`,
:class:`PackedArrays`) with their archive members and checks
(:func:`sector_map_members`, :func:`check_sector_map_members`,
:func:`read_slot_member`), the admission plan (:func:`plan_so_bond`), the
spinful bond vertex (:func:`so_bond_vertex`)
and its exchange image (:func:`gamma_p`), the bond bubble
(:func:`so_bond_bubble`), its RPA dressing (:func:`so_bond_dress`), the
static pairing vertex (:func:`so_bond_pairing_vertex_parts`) and the
pairing kernel of the packed vertex on the pair-S_z sectors
(:class:`PackedBondPairKernel`).

These are array-level building blocks for bond-resolved spin-orbital
calculations (issue #270); no solver entry point calls them yet, and their
signatures may still change when the ``hwave_sc`` and FLEX paths are wired.

In spin-orbital mode (``[mode] enable_spin_orbital``) a bond-resolved
two-body object (the interaction ``Gamma_B``, the bubble ``chibar_B``) is an
``ND x ND`` matrix per momentum, ``ND = B * nd**2``, whose rows and columns
are CODE SLOTS ``(m; x, y)`` (bond channel ``m``, spin-orbital modes ``x``,
``y``; flat position ``m * nd**2 + x * nd + y``). This module fixes what a
code slot means and carries the mesh every phase-aware mapping needs.

Code slots (operator form)
--------------------------
With the bond bilinear ``rho^{(d)}_{XY}(q) = sum_j e^{+iq.j} c+_{X,j}
c_{Y,j+d}``, row slot ``(m; x, y)`` is ``O_{m,xy}(q) = rho^{(sr dr_m)}_{xy}(q)``
and column slot ``(m'; z, w)`` is ``rho^{(sc dr_m')}_{zw}(q)``, with the slot
signs ``(sr, sc) = (-1, -1)``: slot ``(a, b)`` of channel ``m`` has ``a`` in
the origin cell and ``b`` at ``-dr_m``, the geometry of
:class:`hwave.solver.eliashberg_bond.BondPairKernel`. Then

* ``H_int = (1/4N) sum_q sum_rc Gamma_B(q)_rc O_r(q) O_c(q)^dagger`` up to
  one-body terms, for the antisymmetrized ``Gamma_B`` (direct plus crossed
  filing, :func:`so_bond_vertex`); either part filed alone gives ``H_int``
  with ``1/2N``;
* ``chibar_B(q, i nu)_rc = (1/N) int dtau e^{i nu tau}
  <T O_r(q, tau)^dagger O_c(q)>_0`` (the ROW operator is daggered).

Read against the label-form bubble
``chibar_{(d;X,Y),(d';Z,W)}(q, i nu)`` (the correlator of ``rho_{(d;X,Y)}(q)``
and ``rho_{(d';Z,W)}(-q)``), the bubble rule is

``chibar_B(q, i nu)[(m;x,y),(m';z,w)]
= e^{-iq.dr_m} chibar(-q, +i nu)[(+dr_m; y,x), (-dr_m'; z,w)]``,

i.e. the q flip, the row reversal and the row phase all come from
``rho^{(d)}_{XY}(q)^dagger = e^{+iq.d} rho^{(-d)}_{YX}(-q)``; the bosonic
frequency is not flipped. The production bond bubble
(:func:`hwave.solver.bubble.bond_bubble_static`,
``bubble._iter_bond_dynamic``) is pinned to this rule by
``tests/test_so_bond_conventions.py``, which scans every competing reading.

:data:`SLOT_CONVENTION` is the one runtime source of these choices; its
:meth:`SlotConvention.to_string` form is the archive string.
"""

import functools
import logging
import re
import types
from dataclasses import dataclass
from typing import NamedTuple, Optional

import numpy as np

from . import backend as _bk, bond_channels, bubble, flex_bond, matsubara as _ms
from .hartree_fock import NonFiniteError
from .so_records import CompiledInteraction, _displacements

__all__ = ["SlotConvention", "SLOT_CONVENTION", "BondSlotGeometry", "slot_geometry",
           "so_bond_topology", "parse_layout_key", "BandSzCheck", "band_conserves_sz",
           "records_conserve_sz", "select_layout", "canonical_band_eigen", "green_leak",
           "project_green_sectors", "project_borrowed_green", "SECTOR_ORDER", "SectorCross",
           "PPGather",
           "SectorMap", "sz_sector_map", "PackedArrays", "SOBondCaps", "SOBondPlan",
           "plan_so_bond", "estimate_sector_solver_bytes",
           "estimate_packed_kernel_workspace_bytes", "PackedBondBlockStore", "open_store",
           "so_bond_vertex", "gamma_p", "so_bond_bubble", "so_bond_dress",
           "so_bond_pairing_vertex_parts", "so_bond_pairing_vertex_batch", "p_b",
           "PackedBondPairKernel", "sector_map_members", "packed_members",
           "check_sector_map_members", "read_slot_member"]

logger = logging.getLogger(__name__)

_GRAMMAR = ("sr<+|->:sc<+|->:q<+|->:nu<+|->:rowrev<0|1>:colrev<0|1>:"
            "phr<-|0|+>:phc<-|0|+>")
_PATTERN = re.compile(
    r"^sr([+-]):sc([+-]):q([+-]):nu([+-]):rowrev([01]):colrev([01]):phr([-0+]):phc([-0+])$",
    re.ASCII)
_SIGN = {1: "+", -1: "-"}
_PHASE = {1: "+", 0: "0", -1: "-"}


def _is_int(value):
    return isinstance(value, (int, np.integer)) and not isinstance(value, (bool, np.bool_))


@dataclass(frozen=True)
class SlotConvention:
    """The code-slot convention of the bond-resolved bubble and vertex.

    Fields
    ------
    q_sign, nu_sign : +1 or -1
        The code slice ``(q, nu)`` reads the label bubble at
        ``(q_sign q, nu_sign nu)``.
    row_rev, col_rev : 0 or 1
        Whether the row / column mode pair is reversed in the label.
    phase_row, phase_col : -1, 0 or +1
        The code entry carries ``e^{iq.(phase_row dr_m + phase_col dr_m')}``.
    slot_sign : (sr, sc), each +1 or -1
        The operator-form slot signs: row slot ``(m; x, y)`` is
        ``rho^{(sr dr_m)}_{xy}``, column slot ``rho^{(sc dr_m')}``.
    """

    q_sign: int
    nu_sign: int
    row_rev: int
    col_rev: int
    phase_row: int
    phase_col: int
    slot_sign: tuple

    def __post_init__(self):
        for name in ("q_sign", "nu_sign"):
            v = getattr(self, name)
            if not _is_int(v) or int(v) not in (1, -1):
                raise ValueError("SlotConvention: {} must be +1 or -1, got {!r}".format(name, v))
        for name in ("row_rev", "col_rev"):
            v = getattr(self, name)
            if not _is_int(v) or int(v) not in (0, 1):
                raise ValueError("SlotConvention: {} must be 0 or 1, got {!r}".format(name, v))
        for name in ("phase_row", "phase_col"):
            v = getattr(self, name)
            if not _is_int(v) or int(v) not in (-1, 0, 1):
                raise ValueError("SlotConvention: {} must be -1, 0 or +1, got {!r}"
                                 .format(name, v))
        ss = self.slot_sign
        try:
            ok = (len(ss) == 2 and all(_is_int(v) and int(v) in (1, -1) for v in ss))
        except TypeError:
            ok = False
        if not ok:
            raise ValueError("SlotConvention: slot_sign must be a pair of +1/-1, got {!r}"
                             .format(ss))
        for name in ("q_sign", "nu_sign", "row_rev", "col_rev", "phase_row", "phase_col"):
            object.__setattr__(self, name, int(getattr(self, name)))
        object.__setattr__(self, "slot_sign", (int(ss[0]), int(ss[1])))

    def to_string(self):
        """The ASCII archive form
        ``sr<+|->:sc<+|->:q<+|->:nu<+|->:rowrev<0|1>:colrev<0|1>:phr<-|0|+>:phc<-|0|+>``."""
        return "sr{}:sc{}:q{}:nu{}:rowrev{}:colrev{}:phr{}:phc{}".format(
            _SIGN[self.slot_sign[0]], _SIGN[self.slot_sign[1]], _SIGN[self.q_sign],
            _SIGN[self.nu_sign], self.row_rev, self.col_rev, _PHASE[self.phase_row],
            _PHASE[self.phase_col])

    @classmethod
    def from_string(cls, s):
        """Parse :meth:`to_string` output; anything else is a ``ValueError``."""
        match = _PATTERN.fullmatch(s) if isinstance(s, str) else None
        if match is None:
            raise ValueError("slot convention '{}' does not match {}".format(s, _GRAMMAR))
        sr, sc, q, nu, rr, cr, pr, pc = match.groups()
        sign = {"+": 1, "-": -1}
        phase = {"+": 1, "0": 0, "-": -1}
        return cls(q_sign=sign[q], nu_sign=sign[nu], row_rev=int(rr), col_rev=int(cr),
                   phase_row=phase[pr], phase_col=phase[pc], slot_sign=(sign[sr], sign[sc]))

    def label_view(self):
        """The label-index view ``(lambda_r, lambda_c, s_q, s_nu, rowrev,
        colrev, t, p_r, p_c)``: code entry ``[(m;x,y),(m';z,w)](q, nu)`` is
        ``e^{iq.(p_r dr_m + p_c dr_m')} chibar(s_q q, s_nu nu)[(lambda_r dr_m;
        x,y or y,x), (lambda_c dr_m'; z,w or w,z)]`` (``t = 0``: no
        row/column transposition)."""
        return (-self.slot_sign[0], self.slot_sign[1], self.q_sign, self.nu_sign,
                self.row_rev, self.col_rev, 0, self.phase_row, self.phase_col)

    def is_operator_reading(self):
        """Whether the convention is the operator form with a daggered row:
        ``q_sign = -1``, ``nu_sign = +1``, the row (only) reversed, the row
        phase equal to the row slot sign and no column phase."""
        return (self.q_sign == -1 and self.nu_sign == 1 and self.row_rev == 1
                and self.col_rev == 0 and self.phase_row == self.slot_sign[0]
                and self.phase_col == 0)


#: The convention of the production bond bubble and kernel.
SLOT_CONVENTION = SlotConvention(q_sign=-1, nu_sign=+1, row_rev=1, col_rev=0,
                                 phase_row=-1, phase_col=0, slot_sign=(-1, -1))


def _readonly(arr):
    arr.flags.writeable = False
    return arr


@dataclass(frozen=True, eq=False)
class BondSlotGeometry:
    """The bond channels and the mesh of a code-slot array.

    Fields (private, read-only arrays)
    ----------------------------------
    delta_r : (B, 3) int64
        Channel displacements.
    reverse : (B,) int64
        ``delta_r[reverse[m]] == -delta_r[m]``, an involution.
    spatial_shape : (Nx, Ny, Nz)
    q_frac : (nvol, 3) float64
        Momenta in C order as fractions of the reciprocal vectors (grid
        index / shape); ``q = 2 pi q_frac``.
    phase : (nvol, B) complex128
        ``e^{iq.dr_m}``.
    convention : SlotConvention

    Build it with :meth:`from_arrays` (or :func:`slot_geometry`); equality
    is identity.
    """

    delta_r: np.ndarray
    reverse: np.ndarray
    spatial_shape: tuple
    q_frac: np.ndarray
    phase: np.ndarray
    convention: SlotConvention

    @classmethod
    def from_arrays(cls, delta_r, reverse, spatial_shape, convention=SLOT_CONVENTION):
        """Validate and copy the channel arrays and build the q phases.

        Refuses (``ValueError``) a ``delta_r`` that is not an integer
        ``(B, 3)`` array, a ``reverse`` that is not an integer ``(B,)``
        array (ragged or non-numeric input included) or not an involution with
        ``delta_r[reverse[m]] == -delta_r[m]``, a ``spatial_shape`` that is
        not three positive integers, and channels that coincide modulo the
        mesh (the message names both channels and the mesh)."""
        if not isinstance(convention, SlotConvention):
            raise ValueError("BondSlotGeometry: convention must be a SlotConvention, got {!r}"
                             .format(convention))
        try:
            shape = tuple(spatial_shape)
            ok = len(shape) == 3 and all(_is_int(n) and int(n) > 0 for n in shape)
        except TypeError:
            ok = False
        if not ok:
            raise ValueError("BondSlotGeometry: spatial_shape must be three positive integers, "
                             "got {!r}".format(spatial_shape))
        shape = tuple(int(n) for n in shape)

        dr = np.array(delta_r, copy=True)
        if dr.ndim != 2 or dr.shape[1] != 3 or dr.shape[0] < 1 \
                or not np.issubdtype(dr.dtype, np.integer):
            raise ValueError("BondSlotGeometry: delta_r must be an integer (B, 3) array, "
                             "got shape {}".format(dr.shape))
        dr = dr.astype(np.int64)
        B = dr.shape[0]

        bad_reverse = ("BondSlotGeometry: reverse must have shape ({},) with entries in "
                       "[0, {}), got {}")
        try:
            rev = np.array(reverse, copy=True)
        except (TypeError, ValueError):                 # ragged or not array-like
            raise ValueError(bad_reverse.format(B, B, reverse)) from None
        if rev.shape != (B,) or not np.issubdtype(rev.dtype, np.integer) \
                or np.any(rev < 0) or np.any(rev >= B):
            raise ValueError(bad_reverse.format(B, B, rev.tolist()))
        rev = rev.astype(np.int64)
        for m in range(B):
            if rev[rev[m]] != m or np.any(dr[rev[m]] != -dr[m]):
                raise ValueError("BondSlotGeometry: reverse must be an involution with "
                                 "delta_r[reverse[m]] == -delta_r[m]; fails at channel {} "
                                 "(reverse[{}] = {})".format(m, m, int(rev[m])))

        wrapped = np.mod(dr, np.array(shape, dtype=np.int64))
        seen = {}
        for m in range(B):
            key = tuple(int(x) for x in wrapped[m])
            if key in seen:
                o = seen[key]
                raise ValueError("BondSlotGeometry: channels {} (delta_r={}) and {} "
                                 "(delta_r={}) coincide modulo the mesh {}".format(
                                     o, tuple(int(x) for x in dr[o]), m,
                                     tuple(int(x) for x in dr[m]), shape))
            seen[key] = m

        grid = np.indices(shape).reshape(3, -1).T          # C order
        q_frac = grid / np.array(shape, dtype=np.float64)
        phase = np.exp(1j * (2.0 * np.pi * q_frac) @ dr.T.astype(np.float64))
        return cls(delta_r=_readonly(dr), reverse=_readonly(rev), spatial_shape=shape,
                   q_frac=_readonly(np.ascontiguousarray(q_frac, dtype=np.float64)),
                   phase=_readonly(np.ascontiguousarray(phase, dtype=np.complex128)),
                   convention=convention)

    @property
    def nvol(self):
        return int(self.q_frac.shape[0])

    @property
    def B(self):
        return int(self.delta_r.shape[0])

    @property
    def n_channels(self):
        """The channel count ``B``: with :attr:`delta_r` the geometry is a
        bond view of ``eliashberg_bond.BondPairKernel`` and
        ``flex_bond.calc_self_energy_bond``."""
        return self.B


def slot_geometry(topo, spatial_shape, convention=SLOT_CONVENTION):
    """The :class:`BondSlotGeometry` of a bond set on a mesh.

    ``topo`` is any object with ``delta_r`` and ``reverse`` (a
    ``bond_channels.BondTopology``, a ``BondSetView``, or a plain
    namespace). Its channels are rebuilt as a coefficient-free
    ``BondTopology`` and checked by
    ``bond_channels.validate_topology_against_mesh`` (channel 0 on-site,
    the reversal, injectivity modulo the mesh) before the geometry is
    built."""
    rebuilt = bond_channels.BondTopology(topo.delta_r, topo.reverse, {})
    bond_channels.validate_topology_against_mesh(rebuilt, spatial_shape)
    return BondSlotGeometry.from_arrays(rebuilt.delta_r, rebuilt.reverse, spatial_shape,
                                        convention)


def _used_displacements(records):
    """Every particle-hole displacement (and its negative) of the records."""
    used = set()
    for r in records:
        for d in _displacements(tuple(tuple(int(x) for x in c) for c in r["cells"])):
            used.add(d)
            used.add(tuple(-x for x in d))
    return used


def _shell_numbers(resolve, wanted):
    """``{displacement: shell}`` for the off-site displacements ``wanted``.

    ``resolve(k)`` is the topology resolved with ``max_shells = k``; shell
    ``k`` is the set of displacements first kept at ``k`` (the shells
    ``max_shells`` counts), found for ``k = 1, 2, ...`` until every wanted
    displacement has its shell or no shell is left."""
    out = {}
    seen = {(0, 0, 0)}
    total = len(resolve(None).delta_r)
    k = 0
    while set(wanted) - set(out) and len(seen) < total:
        k += 1
        kept = {tuple(int(x) for x in row) for row in resolve(k).delta_r}
        for d in kept - seen:
            out[d] = k
        seen = kept
    return out


def so_bond_topology(compiled, cell, *, max_shells=None):
    """The bond channels of a compiled spin-orbital interaction.

    The channels are the displacements of ``compiled.declared_topology``
    (every particle-hole displacement of every declared monomial,
    zero-valued rows included), resolved by
    ``bond_channels.resolve_bond_topology(..., displacements=)``: channel 0
    on-site, then the off-site shells by length (``cell`` rows are the
    lattice vectors) and lexicographically, closed under reversal.

    ``max_shells`` keeps the first ``max_shells`` off-site shells. A cutoff
    that drops a displacement used by a record of ``compiled.records`` (a
    nonzero coefficient) is refused with a ``ValueError`` naming the shells,
    their displacements and the smallest admissible ``max_shells``; shells
    declared only with zero coupling may be dropped.

    Returns a ``bond_channels.BondTopology`` with ``coeffs == {}``. It
    carries no mesh: :func:`slot_geometry` checks it on one."""
    if not isinstance(compiled, CompiledInteraction):
        raise ValueError("so_bond_topology: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got {}".format(
                             type(compiled).__name__))
    disp = compiled.declared_topology

    def resolve(cutoff):
        return bond_channels.resolve_bond_topology(
            {}, cell, compiled.norb_phys, max_shells=cutoff, active_types=(),
            displacements=disp)

    topo = resolve(max_shells)
    if max_shells is None:
        return topo
    kept = {tuple(int(x) for x in row) for row in topo.delta_r}
    dropped = _used_displacements(compiled.records) - kept
    if dropped:
        shell_of = _shell_numbers(resolve, dropped)
        by_shell = {}
        for d in dropped:
            by_shell.setdefault(shell_of[d], []).append(d)
        parts = ["shell {} with displacements {}".format(k, sorted(by_shell[k]))
                 for k in sorted(by_shell)]
        raise ValueError(
            "so_bond_topology: max_shells={} would drop bond shell(s) used by nonzero "
            "interaction records: {}. Use max_shells >= {}, or remove those "
            "declarations.".format(max_shells, "; ".join(parts), max(by_shell)))
    return topo


# ---------------------------------------------------------------------------
# Layout selection
# ---------------------------------------------------------------------------

#: the values of the ``so_bond_layout`` key
LAYOUT_KEY_VALUES = ("auto", "dense")
_PACKED_REASON = ("so bond layout: sz_packed (S_z conserved: no spin-flip transfer, "
                  "no S_z-changing interaction)")


def parse_layout_key(value, *, key):
    """The ``so_bond_layout`` value, lower case: ``"auto"`` (the default:
    the S_z-packed layout when the input conserves S_z, the dense one
    otherwise) or ``"dense"``. Values are read case-insensitively;
    ``"sz_packed"`` cannot be forced, and anything else (a non-string
    included) is refused (``ValueError``, ``key`` names the setting)."""
    low = value.lower() if isinstance(value, str) else None
    if low in LAYOUT_KEY_VALUES:
        return low
    if low == "sz_packed":
        raise ValueError(
            '{} = "sz_packed" cannot be forced: the S_z-packed layout is chosen automatically '
            'when the input conserves S_z (no spin-flip Transfer entry, no S_z-changing '
            'interaction). Use "auto" (the default) or "dense".'.format(key))
    raise ValueError('{} must be "auto" or "dense", got {!r}'.format(key, value))


class BandSzCheck(NamedTuple):
    """Whether a band conserves S_z; ``entry = (R, (i, j), v)`` is the first
    spin-flip Transfer entry (``None`` when it conserves)."""
    conserving: bool
    entry: Optional[tuple]


def _transfer_key(key, norb_phys):
    """``((Rx, Ry, Rz), (i, j))`` of a Transfer key, indices checked."""
    try:
        R, (i, j) = key
        R = tuple(int(x) for x in R)
        i, j = int(i), int(j)
    except (TypeError, ValueError):
        raise ValueError("band_conserves_sz: a Transfer key must be ((Rx, Ry, Rz), (i, j)), "
                         "got {!r}".format(key)) from None
    nd = 2 * norb_phys
    if len(R) != 3 or not (0 <= i < nd and 0 <= j < nd):
        raise ValueError("band_conserves_sz: Transfer key {!r} is not ((Rx, Ry, Rz), (i, j)) "
                         "with 0 <= i, j < 2·norb_phys = {}".format(key, nd))
    return R, i, j


def band_conserves_sz(transfer, extern, norb_phys):
    """Whether the band of a spin-orbital run conserves S_z.

    ``transfer`` is the Transfer table ``{((Rx, Ry, Rz), (i, j)): v}`` in
    the interleaved spin-orbital index ``i = 2 orb + spin`` (0-based); the
    band conserves S_z iff every spin-off-diagonal entry (``i % 2 != j % 2``)
    is exactly ``0``. ``extern`` (the Extern table in physical orbital
    indices, or ``None``) is accepted and never breaks S_z: its operator
    ``coeff_extern H1 (x) sigma_z`` is spin-diagonal.

    Returns ``BandSzCheck(conserving, entry)``; ``entry = (R, (i, j), v)``
    is the first offending entry in ascending signed key ``(Rx, Ry, Rz, i,
    j)`` (``v`` as a complex number).

    This predicate alone decides whether the Green function is built in the
    canonical per-spin-sector form; it does not read the interaction or the
    ``so_bond_layout`` key. :func:`select_layout` adds the interaction
    records and the override to choose the layout of the bond objects."""
    norb_phys = int(norb_phys)
    bad = []
    for key, v in transfer.items():
        R, i, j = _transfer_key(key, norb_phys)
        if i % 2 != j % 2 and complex(v) != 0:
            bad.append((R + (i, j), complex(v)))
    if not bad:
        return BandSzCheck(True, None)
    k, v = min(bad, key=lambda e: e[0])
    return BandSzCheck(False, (k[:3], (k[3], k[4]), v))


def records_conserve_sz(compiled):
    """The first record of ``compiled.records`` with ``S_z(1) - S_z(2) +
    S_z(3) - S_z(4) != 0`` (spin-block modes: ``S_z = +1/2`` for ``X <
    norb_phys``, ``-1/2`` otherwise), or ``None`` when every record
    conserves S_z."""
    if not isinstance(compiled, CompiledInteraction):
        raise ValueError("records_conserve_sz: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got {}".format(
                             type(compiled).__name__))
    norb = int(compiled.norb_phys)
    for r in compiled.records:
        s = [1 - 2 * (int(x) // norb) for x in r["modes"]]      # 2 S_z
        if s[0] - s[1] + s[2] - s[3] != 0:
            return r
    return None


def select_layout(compiled, transfer_table, *, override):
    """The layout of the spin-orbital bond objects and its reason line.

    ``override`` is the ``so_bond_layout`` value (:func:`parse_layout_key`).
    Precedence: ``"dense"`` forces the dense layout; otherwise the first
    spin-flip Transfer entry (:func:`band_conserves_sz`), then the first
    S_z-changing interaction record (:func:`records_conserve_sz`), selects
    it; with neither the S_z-packed layout ``"sz_packed"`` is selected. The
    selection is exact: an entry of ``1e-14`` selects the dense layout.

    Returns ``(layout, reason)``; the caller logs ``reason`` at INFO."""
    forced = parse_layout_key(override, key="so_bond_layout") == "dense"
    if not isinstance(compiled, CompiledInteraction):
        raise ValueError("select_layout: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got {}".format(
                             type(compiled).__name__))
    if forced:
        return "dense", 'so bond layout: dense (reason: forced by so_bond_layout = "dense")'
    band = band_conserves_sz(transfer_table, None, compiled.norb_phys)
    if not band.conserving:
        R, (i, j), v = band.entry
        return "dense", ("so bond layout: dense (reason: spin-flip transfer entry R=({}, {}, {}) "
                         "({}, {}) = {!r})".format(R[0], R[1], R[2], i + 1, j + 1, complex(v)))
    rec = records_conserve_sz(compiled)
    if rec is not None:
        typ = str(rec["type"])
        if typ == "PairLift":
            return "dense", "so bond layout: dense (reason: PairLift declared)"
        return "dense", ("so bond layout: dense (reason: S_z-changing interaction record of "
                         "type {})".format(typ))
    return "sz_packed", _PACKED_REASON


# ---------------------------------------------------------------------------
# The canonical band Green function and borrowed Green functions
# ---------------------------------------------------------------------------

def _spin_halves(norb_phys):
    """The spin-up and spin-down slices of the spin-block index."""
    return slice(0, norb_phys), slice(norb_phys, 2 * norb_phys)


def _check_norb_phys(helper, norb_phys):
    if not _is_int(norb_phys) or int(norb_phys) < 1:
        raise ValueError("{}: norb_phys must be a positive integer, got {!r}".format(
            helper, norb_phys))
    return int(norb_phys)


def canonical_band_eigen(epsilon_k, norb_phys, *, conserving):
    """Eigen-decomposition of the spin-orbital band, per spin sector when the
    band conserves S_z.

    ``epsilon_k`` is the band ``(nd, nd, Nx, Ny, Nz)`` in the SPIN-BLOCK
    index (``X = s norb_phys + a``, ``nd = 2 norb_phys``), as
    ``sc._so_epsilon_k`` returns it. ``conserving`` must be the bool
    ``band_conserves_sz(...).conserving`` of the run's Transfer table; it
    alone decides the construction (the ``so_bond_layout`` is never read).

    * ``conserving = True``: the spin-up and spin-down blocks are
      diagonalized separately and assembled with structural zeros between
      them, so a Green function built from the result has exactly zero
      S_z-off-sector entries. Per k the ``nd`` eigenvalues are put in
      ascending order by a stable sort (ties keep spin up before spin down)
      with the matching column permutation of the eigenvectors, the order
      of ``sc._calc_eigenvalues``. A band with a nonzero S_z-off-sector
      entry is refused: the table conserves S_z, so the array is not in the
      spin-block order (an interleaved ``2 orb + spin`` array, for example);
      a non-finite band is refused as such first.
    * ``conserving = False``: the full diagonalization of
      ``sc._calc_eigenvalues`` (the same arithmetic).

    Returns ``(eigenvalues (Nx, Ny, Nz, nd), eigenvectors (Nx, Ny, Nz, nd,
    nd))``; ``eigenvectors[..., :, i]`` belongs to ``eigenvalues[..., i]``."""
    norb_phys = _check_norb_phys("canonical_band_eigen", norb_phys)
    if not isinstance(conserving, (bool, np.bool_)):
        raise ValueError("canonical_band_eigen: conserving must be a bool (band_conserves_sz("
                         "...).conserving), got {}".format(type(conserving).__name__))
    eps = np.asarray(epsilon_k)
    if eps.ndim != 5 or eps.shape[0] != eps.shape[1]:
        raise ValueError("canonical_band_eigen: epsilon_k must have shape (nd, nd, Nx, Ny, Nz), "
                         "got {}".format(eps.shape))
    nd = 2 * norb_phys
    if eps.shape[0] != nd:
        raise ValueError("canonical_band_eigen: epsilon_k has nd = {}, expected 2 * norb_phys "
                         "= {}".format(eps.shape[0], nd))
    batch = eps.transpose(2, 3, 4, 0, 1)
    if not conserving:
        return np.linalg.eigh(batch)
    if not np.all(np.isfinite(batch)):
        raise ValueError("canonical_band_eigen: epsilon_k has non-finite entries")
    up, dn = _spin_halves(norb_phys)
    off = max(float(np.max(np.abs(batch[..., up, dn]))),
              float(np.max(np.abs(batch[..., dn, up]))))
    if off != 0.0:
        raise ValueError(
            "canonical_band_eigen: epsilon_k has S_z-off-sector entries (max %.3e) in the "
            "spin-block reading although the Transfer table conserves S_z (is it in the "
            "interleaved order?)" % off)
    w_up, v_up = np.linalg.eigh(batch[..., up, up])
    w_dn, v_dn = np.linalg.eigh(batch[..., dn, dn])
    values = np.concatenate([w_up, w_dn], axis=-1)
    vectors = np.zeros(batch.shape, dtype=np.result_type(v_up, v_dn))
    vectors[..., up, up] = v_up
    vectors[..., dn, dn] = v_dn
    order = np.argsort(values, axis=-1, kind="stable")
    values = np.take_along_axis(values, order, axis=-1)
    vectors = np.take_along_axis(vectors, order[..., None, :], axis=-1)
    return values, vectors


def _green_view(helper, G, norb_phys, axes):
    """``G`` with its two spin-orbital axes moved last (a view)."""
    norb_phys = _check_norb_phys(helper, norb_phys)
    G = np.asarray(G)
    try:
        a, b = axes
        if not (_is_int(a) and _is_int(b) and -G.ndim <= a < G.ndim and -G.ndim <= b < G.ndim):
            raise ValueError
    except (TypeError, ValueError):
        raise ValueError("{}: axes must be two axes of the {}-dimensional Green function, got "
                         "{!r}".format(helper, G.ndim, axes)) from None
    a, b = int(a) % G.ndim, int(b) % G.ndim
    if a == b:
        raise ValueError("{}: axes must be two different axes, got {!r}".format(helper, axes))
    view = np.moveaxis(G, (a, b), (-2, -1))
    nd = 2 * norb_phys
    if view.shape[-2:] != (nd, nd):
        raise ValueError("{}: the spin-orbital axes {!r} have lengths {}, expected 2 * norb_phys "
                         "= {}".format(helper, tuple(axes), view.shape[-2:], nd))
    return view, norb_phys


def green_leak(G, norb_phys, *, axes=(-2, -1)):
    """``max |G_off-sector| / max |G|`` of a spin-block Green function (``0``
    when ``G`` is zero; ``nan`` when ``G`` has a non-finite entry). ``axes``
    are its two spin-orbital axes; the off-sector entries are those between
    the spin-up and spin-down halves."""
    view, norb_phys = _green_view("green_leak", G, norb_phys, axes)
    up, dn = _spin_halves(norb_phys)
    scale = float(np.max(np.abs(view))) if view.size else 0.0
    if not np.isfinite(scale):
        return float("nan")
    if scale == 0.0:
        return 0.0
    off = max(float(np.max(np.abs(view[..., up, dn]))), float(np.max(np.abs(view[..., dn, up]))))
    return off / scale


def project_green_sectors(G, norb_phys, *, axes=(-2, -1)):
    """A copy of the spin-block Green function ``G`` with its S_z-off-sector
    entries set to exact zeros (``axes``: its two spin-orbital axes)."""
    out = np.array(G, copy=True)
    view, norb_phys = _green_view("project_green_sectors", out, norb_phys, axes)
    up, dn = _spin_halves(norb_phys)
    view[..., up, dn] = 0
    view[..., dn, up] = 0
    return out


def project_borrowed_green(G, *, norb_phys, source, axes=(-2, -1)):
    """The Green function a run that conserves S_z uses for a borrowed
    spin-block Green function ``G`` (read from ``source``), and its leak.

    With ``leak = green_leak(G)``:

    * ``leak = 0``: ``G`` itself is used;
    * ``0 < leak <= 1e-12``: numerical leakage from a non-canonical
      producer; ``project_green_sectors(G)`` is used and one WARNING gives
      the measured leak;
    * ``leak > 1e-12`` (or a non-finite ``G``): refused, the Green function
      was produced from a different Hamiltonian.

    Returns ``(green_used, leak)``. The helper checks and projects only:
    comparing the file's recorded digest and recording the digest of
    ``green_used`` stay with the caller."""
    leak = green_leak(G, norb_phys, axes=axes)
    if not np.isfinite(leak):
        raise ValueError("the Green function '%s' has non-finite entries" % (source,))
    if leak == 0.0:
        return G, leak
    if leak <= _LEAK_NUMERICAL:
        logger.warning(
            "the Green function '%s' has S_z-off-sector entries at the round-off level (leak = "
            "max off-sector |G| / max |G| = %.3e <= 1e-12: numerical leakage from a "
            "non-canonical producer); they are projected to exact zeros", source, leak)
        return project_green_sectors(G, norb_phys, axes=axes), leak
    raise ValueError(
        "the Green function '%s' breaks S_z (leak = max off-sector |G| / max |G| = %.3e > "
        "1e-12) although this run's input conserves S_z; it was produced from a different "
        "Hamiltonian. Regenerate it from this run's input." % (source, leak))


# ---------------------------------------------------------------------------
# S_z sectors: the sector map and packed arrays
# ---------------------------------------------------------------------------

#: the ph sectors in storage order, their sigma values and the pp sectors
SECTOR_ORDER = ("0", "+1", "-1")
_SIGMAS = (0, 1, -1)
_PP_ORDER = (1, 0, -1)
#: consumers of the sector map; each reads it only through :func:`_consumer_map`
_MAP_CONSUMERS = ("gamma_packer", "gamma_p", "store_scatter", "dress_sectors", "w_cross",
                  "kernel_gather")
_LEAK_NUMERICAL = 1e-12


class SectorCross(NamedTuple):
    """A cross-sector slot permutation onto one target sector: target local
    position ``i`` reads local position ``perm[i]`` of sector ``source``
    (a sigma value)."""
    source: int
    perm: np.ndarray


class PPGather(NamedTuple):
    """The gather table of one pp sector (``n_p`` gap pairs): for output
    pair ``i = (A, D)`` and input pair ``j = (B, C)``, ``sector[i, j]`` is
    the ph sector (sigma) of ``S_z(A) - S_z(B)`` (equal to that of
    ``S_z(C) - S_z(D)``), ``row_local[a, i, j]`` the local position of slot
    ``(a; A, B)`` and ``col_local[b, i, j]`` that of ``(b; C, D)`` in it."""
    sector: np.ndarray
    row_local: np.ndarray
    col_local: np.ndarray


@dataclass(frozen=True, eq=False)
class SectorMap:
    """Every index map of the S_z-sector packed layout of code-slot arrays.

    Built (and validated) by :func:`sz_sector_map`; the record itself does
    not validate. All arrays are read-only.

    The ph sector of a code slot ``(m; x, y)`` (spin-block modes, ``S_z =
    +1/2`` for ``x < norb_phys``) is ``sigma = S_z(x) - S_z(y)``; an
    S_z-conserving two-body object is block diagonal in it.

    Fields
    ------
    B, nd, norb_phys, ND : int
        Channels, modes ``nd = 2 norb_phys``, ``ND = B nd**2``.
    sector_order : ``("0", "+1", "-1")``; sigmas : ``(0, 1, -1)``
        The storage order of the sectors. A sigma value used as a Python
        index of a 3-sequence in this order addresses its sector.
    sizes : (3,) tuple of int
        ``(2 B norb_phys**2, B norb_phys**2, B norb_phys**2)``.
    offsets : int64 (4,)
        Sector ``s`` (a storage POSITION) occupies
        ``index[offsets[s]:offsets[s + 1]]``. Never index ``offsets`` with a
        sigma value; use ``s = sigmas.index(sigma)``.
    index : int64 (ND,)
        The dense slot ``I = m nd**2 + x nd + y`` of every packed position,
        sectors concatenated, ascending within a sector.
    sector_of : int8 (ND,); local : int64 (ND,)
        Per dense slot: its sigma and its position within its sector.
    pb_cross : 3 :class:`SectorCross`, per target sector
        ``P_B`` (``(m; x, y) -> (m; y, x)``, sigma -> -sigma):
        ``P_B(M)@t = M@source[perm][:, perm]``.
    gp_cross : 3 :class:`SectorCross`, per target sector
        The exchange-image slot map ``(m; x, y) -> (reverse[m]; y, x)``.
    pp_order : ``(1, 0, -1)``
        The pp sectors ``p = S_z(b) + S_z(c)`` of a gap pair ``(b, c)``.
    pp_pairs : 3 int64 arrays, per pp sector
        Ascending flat pairs ``b nd + c``; sizes ``norb_phys**2``,
        ``2 norb_phys**2``, ``norb_phys**2``.
    pp_gather : 3 :class:`PPGather`, per pp sector
        ``V^p_{ab}[i, j] = V@sector[i, j][row_local[a, i, j],
        col_local[b, i, j]]``.
    """

    B: int
    nd: int
    norb_phys: int
    ND: int
    sector_order: tuple
    sigmas: tuple
    sizes: tuple
    offsets: np.ndarray
    index: np.ndarray
    sector_of: np.ndarray
    local: np.ndarray
    pb_cross: tuple
    gp_cross: tuple
    pp_order: tuple
    pp_pairs: tuple
    pp_gather: tuple


def _consumer_map(name, sector_map):
    """The sector map as consumer ``name`` reads it (the identity). Every
    packed consumer obtains its map here, so one consumer's reading can be
    replaced in isolation."""
    if name not in _MAP_CONSUMERS:
        raise ValueError("_consumer_map: unknown consumer {!r}, expected one of {}".format(
            name, _MAP_CONSUMERS))
    return sector_map


def _sector_reverse(topo, nd):
    if not (_is_int(nd) and not isinstance(nd, bool) and int(nd) > 0 and int(nd) % 2 == 0):
        raise ValueError("sz_sector_map: nd must be a positive even integer "
                         "(nd = 2*norb_phys), got {!r}".format(nd))
    rev = np.asarray(topo.reverse)
    B = rev.shape[0] if rev.ndim == 1 else -1
    ok = (rev.ndim == 1 and B >= 1 and np.issubdtype(rev.dtype, np.integer)
          and bool(np.all((rev >= 0) & (rev < B))))
    if ok:
        rev = rev.astype(np.int64)
        ok = bool(np.all(rev[rev] == np.arange(B)))
    if not ok:
        raise ValueError("sz_sector_map: topo.reverse must be an integer (B,) involution with "
                         "entries in [0, B), got {}".format(rev.tolist()))
    return int(nd), rev


def sz_sector_map(topo, nd):
    """The :class:`SectorMap` of the code slots of a bond set.

    ``topo`` carries ``reverse`` (``(B,)``, an involution; a
    ``BondTopology``, a geometry or a namespace); ``nd = 2 norb_phys``.
    Cached by content ``(B, nd, reverse)``: equal inputs return the same
    object."""
    nd, rev = _sector_reverse(topo, nd)
    return _build_sector_map(int(rev.shape[0]), nd, rev.tobytes())


@functools.lru_cache(maxsize=64)
def _build_sector_map(B, nd, rev_bytes):
    rev = np.frombuffer(rev_bytes, dtype=np.int64)
    norb = nd // 2
    npair = nd * nd
    ND = B * npair
    two_sz = np.where(np.arange(nd) < norb, 1, -1)            # 2 S_z per mode
    pair_sigma = ((two_sz[:, None] - two_sz[None, :]) // 2).reshape(npair)
    slot_sigma = np.tile(pair_sigma, B)                         # per dense slot
    lists = [np.flatnonzero(slot_sigma == s).astype(np.int64) for s in _SIGMAS]
    sizes = tuple(int(len(l)) for l in lists)
    offsets = np.concatenate(([0], np.cumsum(sizes))).astype(np.int64)
    index = np.concatenate(lists).astype(np.int64)
    sector_of = slot_sigma.astype(np.int8)
    local = np.empty(ND, dtype=np.int64)
    for l in lists:
        local[l] = np.arange(len(l), dtype=np.int64)

    m = index // npair
    x, y = np.divmod(index % npair, nd)
    swapped = m * npair + y * nd + x                            # (m; y, x)
    crossed = rev[m] * npair + y * nd + x                       # (rev m; y, x)

    def cross(target_slots):
        out = []
        for s, sigma in enumerate(_SIGMAS):
            seg = target_slots[offsets[s]:offsets[s + 1]]
            out.append(SectorCross(-sigma, _readonly(local[seg].copy())))
        return tuple(out)

    pp_pairs, pp_gather = [], []
    pair_p = ((two_sz[:, None] + two_sz[None, :]) // 2).reshape(npair)
    chan = np.arange(B, dtype=np.int64)[:, None, None] * npair
    for p in _PP_ORDER:
        pairs = np.flatnonzero(pair_p == p).astype(np.int64)
        A, D = np.divmod(pairs, nd)                             # output pair i = (A, D)
        Bm, C = np.divmod(pairs, nd)                            # input pair j = (B, C)
        row = A[:, None] * nd + Bm[None, :]                     # pair (A, B) per (i, j)
        col = C[None, :] * nd + D[:, None]                      # pair (C, D) per (i, j)
        pp_pairs.append(_readonly(pairs))
        pp_gather.append(PPGather(_readonly(pair_sigma[row].astype(np.int8)),
                                  _readonly(local[chan + row[None]]),
                                  _readonly(local[chan + col[None]])))

    return SectorMap(B=B, nd=nd, norb_phys=norb, ND=ND, sector_order=SECTOR_ORDER,
                     sigmas=_SIGMAS, sizes=sizes, offsets=_readonly(offsets),
                     index=_readonly(index), sector_of=_readonly(sector_of),
                     local=_readonly(local), pb_cross=cross(swapped), gp_cross=cross(crossed),
                     pp_order=_PP_ORDER, pp_pairs=tuple(pp_pairs), pp_gather=tuple(pp_gather))


def _sector_blocks(dense, sector_map):
    """The diagonal sector blocks of a dense ``(..., ND, ND)`` code-slot
    array (fresh, C-contiguous complex128) and ``max|off-sector|`` (``0.0``
    when every off-sector entry is exactly 0; ``nan`` when one is not
    finite)."""
    sector_of = sector_map.sector_of
    blocks, off_max = [], 0.0
    for s, sigma in enumerate(sector_map.sigmas):
        idx = sector_map.index[sector_map.offsets[s]:sector_map.offsets[s + 1]]
        rows = np.take(dense, idx, axis=-2)
        if off_max == off_max:                          # not yet nan
            off = np.take(rows, np.flatnonzero(sector_of != sigma), axis=-1)
            if off.size and np.any(off != 0):
                v = float(np.max(np.abs(off)))
                off_max = max(off_max, v) if np.isfinite(v) else float("nan")
        blocks.append(np.ascontiguousarray(np.take(rows, idx, axis=-1), dtype=np.complex128))
    return blocks, off_max


def _leak_refusal(source, off_max, scale):
    """The refusal of a dense array with nonzero S_z-off-sector entries:
    ``leak = off_max / scale`` (``max|off-sector| / max|X|``; ``nan`` when
    either is not finite), classified "numerical leakage from a
    non-canonical producer" for ``leak <= 1e-12`` and "incompatible
    Hamiltonian" above (and for ``nan``)."""
    leak = off_max / scale if np.isfinite(scale) else float("nan")
    cls_text = ("numerical leakage from a non-canonical producer"
                if leak <= _LEAK_NUMERICAL else
                "incompatible Hamiltonian (produced from a different, S_z-breaking input)")
    return ValueError(
        "'%s' is dense with nonzero S_z-off-sector entries (measured leak = %.3e; "
        "packing requires exactly 0) although this run conserves S_z: %s. "
        "Regenerate it with this version." % (source, leak, cls_text))


def _block_name(s):
    return SECTOR_ORDER[s] if 0 <= s < len(SECTOR_ORDER) else str(s)


@dataclass(frozen=True, eq=False)
class PackedArrays:
    """A code-slot array in the S_z-sector packed layout.

    ``blocks`` holds the three diagonal sector blocks in the order
    ``sector_map.sector_order``; block ``s`` has shape ``lead_shape +
    (sizes[s], sizes[s])`` and is complex128. The constructor checks these
    invariants (``ValueError`` naming the first violation) and keeps the
    given arrays as they are.

    Ownership: the blocks are fresh allocations only when the object comes
    from :meth:`zeros` or :meth:`from_dense`. Packed arrays handed out by a
    block store hold VIEWS of the store's slots: writes through them reach
    the store, and they must not be used once the slot is released (not
    enforced: the views keep the memory alive but no longer belong to the
    store).
    :meth:`to_dense` always returns a new array (unless ``out`` is given).
    """

    sector_map: SectorMap
    blocks: tuple

    def __post_init__(self):
        sm = self.sector_map
        if not isinstance(sm, SectorMap):
            raise ValueError("PackedArrays: sector_map must be a SectorMap, got {}".format(
                type(sm).__name__))
        if tuple(sm.sector_order) != SECTOR_ORDER:
            raise ValueError("PackedArrays: sector_map.sector_order must be {!r}, got {!r}"
                             .format(SECTOR_ORDER, sm.sector_order))
        if tuple(sm.sigmas) != _SIGMAS:
            raise ValueError("PackedArrays: sector_map.sigmas must be {!r}, got {!r}"
                             .format(_SIGMAS, sm.sigmas))
        if len(sm.sizes) != 3:
            raise ValueError("PackedArrays: sector_map.sizes must have 3 entries, got {!r}"
                             .format(sm.sizes))
        blocks = self.blocks
        if not isinstance(blocks, tuple) or len(blocks) != 3:
            got = ("{} blocks".format(len(blocks)) if isinstance(blocks, tuple)
                   else type(blocks).__name__)
            raise ValueError("PackedArrays: blocks must be a tuple of exactly 3 arrays "
                             "(sectors '0', '+1', '-1'), got {}".format(got))
        lead0 = None
        for s, b in enumerate(blocks):
            name = _block_name(s)
            if not isinstance(b, np.ndarray) or b.dtype != np.complex128:
                got = str(b.dtype) if isinstance(b, np.ndarray) else type(b).__name__
                raise ValueError("PackedArrays: block {} (sector '{}') must be a complex128 "
                                 "ndarray, got {}".format(s, name, got))
            n = int(sm.sizes[s])
            if b.ndim < 2 or b.shape[-2:] != (n, n):
                raise ValueError("PackedArrays: block {} (sector '{}') must end in shape "
                                 "({}, {}), got shape {}".format(s, name, n, n, b.shape))
            if lead0 is None:
                lead0 = b.shape[:-2]
            elif b.shape[:-2] != lead0:
                raise ValueError("PackedArrays: block {} (sector '{}') has leading shape {}, "
                                 "block 0 has {}".format(s, name, b.shape[:-2], lead0))

    @classmethod
    def zeros(cls, sector_map, lead_shape):
        """Fresh zero blocks of leading shape ``lead_shape``."""
        lead = tuple(int(n) for n in lead_shape)
        return cls(sector_map, tuple(np.zeros(lead + (n, n), dtype=np.complex128)
                                     for n in sector_map.sizes))

    @classmethod
    def from_dense(cls, dense, sector_map, *, source):
        """Pack a dense ``(..., ND, ND)`` code-slot array (fresh blocks).

        The off-sector entries must be exactly zero; otherwise the array is
        refused with the measured leak ``max|off-sector| / max|dense|`` and
        its classification (``<= 1e-12``: numerical leakage from a
        non-canonical producer; above: an incompatible Hamiltonian).
        ``source`` names the array in the message."""
        if not isinstance(sector_map, SectorMap):
            raise ValueError("PackedArrays.from_dense: sector_map must be a SectorMap, got "
                             "{}".format(type(sector_map).__name__))
        dense = np.asarray(dense)
        ND = sector_map.ND
        if dense.ndim < 2 or dense.shape[-2:] != (ND, ND):
            raise ValueError("PackedArrays.from_dense: '{}' must end in shape (ND, ND) = "
                             "({}, {}), got shape {}".format(source, ND, ND, dense.shape))
        blocks, off_max = _sector_blocks(dense, sector_map)
        if off_max != 0.0:
            raise _leak_refusal(source, off_max, float(np.max(np.abs(dense))))
        return cls(sector_map, tuple(blocks))

    def __getitem__(self, sigma):
        """The block of sector ``sigma``: ``0``, ``+1``, ``-1`` or ``"0"``,
        ``"+1"``, ``"-1"``."""
        sm = self.sector_map
        if isinstance(sigma, str) and sigma in sm.sector_order:
            return self.blocks[sm.sector_order.index(sigma)]
        if type(sigma) is int and sigma in sm.sigmas:
            return self.blocks[sm.sigmas.index(sigma)]
        raise KeyError("PackedArrays: sector must be 0, +1, -1 or '0', '+1', '-1', got "
                       "{!r}".format(sigma))

    @property
    def lead_shape(self):
        return tuple(self.blocks[0].shape[:-2])

    @property
    def nbytes(self):
        return int(sum(b.nbytes for b in self.blocks))

    def to_dense(self, out=None):
        """The dense ``lead_shape + (ND, ND)`` array (off-sector zeros); a new
        array, or ``out`` (complex128, that shape) overwritten and returned."""
        sm = self.sector_map
        shape = self.lead_shape + (sm.ND, sm.ND)
        if out is None:
            out = np.zeros(shape, dtype=np.complex128)
        else:
            if not isinstance(out, np.ndarray) or out.dtype != np.complex128 \
                    or out.shape != shape:
                raise ValueError("PackedArrays.to_dense: out must be a complex128 array of "
                                 "shape {}, got {} {}".format(
                                     shape, getattr(out, "dtype", type(out).__name__),
                                     getattr(out, "shape", None)))
            out[...] = 0
        for s, b in enumerate(self.blocks):
            idx = sm.index[sm.offsets[s]:sm.offsets[s + 1]]
            out[..., idx[:, None], idx] = b
        return out


# ---------------------------------------------------------------------------
# Archive members of the packed layout
# ---------------------------------------------------------------------------

_LAYOUT_SET = "'dense' or 'sz_packed'"
_MAP_MEMBERS = ("sz_sector_order", "sz_sector_offsets", "sz_sector_index")


def sector_map_members(sector_map):
    """The archive members that declare the slot layout.

    ``None`` (dense layout): ``{"block_layout": "dense"}``. A
    :class:`SectorMap` (packed layout): ``block_layout = "sz_packed"``
    (``<U16``), ``sz_sector_order`` (``<U2 (3,)``, ``["0", "+1", "-1"]``),
    ``sz_sector_offsets`` (int64 ``(4,)``, cumulative sector sizes) and
    ``sz_sector_index`` (int64 ``(ND,)``, the dense slot ``I = m nd**2 + A nd
    + C`` of every packed position). The arrays are fresh copies."""
    if sector_map is None:
        return {"block_layout": np.array("dense", dtype="<U16")}
    if not isinstance(sector_map, SectorMap):
        raise ValueError("sector_map_members: sector_map must be a SectorMap or None, got "
                         "{}".format(type(sector_map).__name__))
    return {"block_layout": np.array("sz_packed", dtype="<U16"),
            "sz_sector_order": np.array(sector_map.sector_order, dtype="<U2"),
            "sz_sector_offsets": np.array(sector_map.offsets, dtype=np.int64),
            "sz_sector_index": np.array(sector_map.index, dtype=np.int64)}


def _packed_key(name, sector):
    return "{}@{}".format(name, sector)


def packed_members(name, packed):
    """The members ``<name>@0``, ``<name>@+1``, ``<name>@-1`` of a
    :class:`PackedArrays` (its blocks, in sector order, not copied)."""
    if not isinstance(name, str) or not name or "@" in name:
        raise ValueError("packed_members: name must be a non-empty string without '@', got "
                         "{!r}".format(name))
    if not isinstance(packed, PackedArrays):
        raise ValueError("packed_members: packed must be a PackedArrays, got {}".format(
            type(packed).__name__))
    return {_packed_key(name, s): b for s, b in zip(packed.sector_map.sector_order,
                                                    packed.blocks)}


def _render(value):
    """A stored or expected entry as printed in a refusal: strings quoted,
    numbers as they are."""
    if isinstance(value, (str, np.str_)):
        return repr(str(value))
    if isinstance(value, np.generic):
        return str(value.item())
    return str(value)


def _member_values(arr, expected_dtype):
    """``arr`` comparable by value with a map member of ``expected_dtype``,
    or ``None`` when its dtype cannot hold such values: an integer member
    accepts any signed or unsigned integer dtype (cast to int64; an
    unsigned value above the int64 range is refused), a unicode member any
    width and byte order (compared as stored, never truncated). Float,
    complex, bool, object and bytes members are refused."""
    want = np.dtype(expected_dtype).kind
    kind = arr.dtype.kind
    if want == "i" and kind in "iu":
        if kind == "u" and arr.size and int(arr.max()) > np.iinfo(np.int64).max:
            return None
        return arr.astype(np.int64, copy=False)
    if want == "U" and kind == "U":
        return arr
    return arr if arr.dtype == expected_dtype else None


def _first_difference(stored, expected):
    """``(position, stored, expected)`` (rendered) of the first difference
    between a stored member (``None`` when absent) and the expected 1-D
    array, or ``None`` when they are equal (by value: see
    :func:`_member_values`). Sentinels: an absent member is ``missing`` at
    position 0; a dtype that cannot hold the values or a wrong shape is
    ``dtype <dt>`` or ``shape <shape>`` on both sides at position 0; a
    stored member shorter
    than the map (with an equal common part) is ``end of member`` at the
    first missing position, a longer one is ``end of map`` on the expected
    side at the first extra position."""
    if stored is None:
        return 0, "missing", _render(expected[0])
    arr = _member_values(np.asarray(stored), expected.dtype)
    if arr is None:
        arr = np.asarray(stored)
        return 0, "dtype {}".format(arr.dtype), "dtype {}".format(expected.dtype)
    if arr.ndim != 1:
        return 0, "shape {}".format(arr.shape), "shape {}".format(expected.shape)
    n = min(arr.shape[0], expected.shape[0])
    differ = np.flatnonzero(arr[:n] != expected[:n])
    if differ.size:
        i = int(differ[0])
        return i, _render(arr[i]), _render(expected[i])
    if arr.shape[0] < expected.shape[0]:
        return n, "end of member", _render(expected[n])
    if arr.shape[0] > expected.shape[0]:
        return n, _render(arr[n]), "end of map"
    return None


def _stored_layout(data):
    """``block_layout`` of an archive as a ``str``, or ``None`` when it is
    absent or not a unicode scalar."""
    if "block_layout" not in data:
        return None
    v = np.asarray(data["block_layout"])
    return str(v[()]) if v.ndim == 0 and v.dtype.kind == "U" else None


def _layout_text(data):
    """The stored ``block_layout`` as printed in a refusal: ``missing``, the
    value (strings quoted) or ``shape <shape>`` for a non-scalar."""
    if "block_layout" not in data:
        return "missing"
    v = np.asarray(data["block_layout"])
    return _render(v[()]) if v.ndim == 0 else "shape {}".format(v.shape)


def check_sector_map_members(data, B, nd, *, file, producer_label):
    """Validate the layout members of an archive ``data`` (a mapping, e.g.
    an ``np.load`` result) against ``B`` channels (``len(delta_r)``) and
    ``nd = 2 norb_phys`` modes; return the layout, ``"dense"`` or
    ``"sz_packed"``.

    For ``"sz_packed"`` the members ``sz_sector_order``,
    ``sz_sector_offsets`` and ``sz_sector_index`` must equal, in shape and
    every entry, the map recomputed from ``(B, nd)`` (it does not depend on
    ``reverse``); the comparison is by value -- any integer dtype for the
    integer members, any unicode width or byte order for
    ``sz_sector_order``, while a float, complex, bool, object or bytes
    member is refused as a dtype. The first difference is refused as

        '<file>' is malformed or truncated: <member> differs from the map
        of (B = <B>, nd = <nd>) at position <i> (<stored> vs <expected>).
        Regenerate it with <producer_label>.

    with these positions and values: an entry that differs is reported at
    its position, strings quoted; an absent member is ``missing`` at
    position 0; a wrong dtype or shape is ``dtype <dt>`` / ``shape
    <shape>`` on both sides at position 0; a truncated member is ``end of
    member`` at its first missing position and a longer one ``end of map``
    (expected side) at its first extra position. ``block_layout`` absent or
    outside the set is reported at position 0 against ``'dense' or
    'sz_packed'`` (``missing``, or the stored value)."""
    if not (_is_int(B) and not isinstance(B, bool) and int(B) > 0):
        raise ValueError("check_sector_map_members: B must be a positive integer, got "
                         "{!r}".format(B))
    B = int(B)
    sm = sz_sector_map(types.SimpleNamespace(reverse=np.arange(B, dtype=np.int64)), nd)

    def malformed(member, i, stored, expected):
        return ValueError(
            "'{}' is malformed or truncated: {} differs from the map of (B = {}, nd = {}) at "
            "position {} ({} vs {}). Regenerate it with {}.".format(
                file, member, B, sm.nd, i, stored, expected, producer_label))

    layout = _stored_layout(data)
    if layout not in LAYOUTS:
        raise malformed("block_layout", 0, _layout_text(data), _LAYOUT_SET)
    if layout == "dense":
        return layout
    expected = sector_map_members(sm)
    for member in _MAP_MEMBERS:
        diff = _first_difference(data[member] if member in data else None, expected[member])
        if diff is not None:
            raise malformed(member, *diff)
    return layout


def _member(data, key, file):
    if key not in data:
        raise ValueError("read_slot_member: '{}' has no member '{}'".format(file, key))
    return np.asarray(data[key])


def read_slot_member(data, name, *, want, sector_map, file):
    """The code-slot member ``name`` of an archive ``data`` in the layout
    ``want`` (``"dense"``: an ndarray ``(..., ND, ND)``; ``"sz_packed"``: a
    :class:`PackedArrays` on ``sector_map``).

    The archive's layout members must have passed
    :func:`check_sector_map_members`; a packed file must in addition carry
    exactly the index map of ``sector_map``.

    * Same layout: the stored arrays as they are.
    * Packed file, dense consumer: expanded slice by slice along the first
      axis; off-sector entries are exactly 0.
    * Dense file, packed consumer: packed slice by slice; every off-sector
      entry must be exactly 0. Otherwise the file is refused with the leak
      of the whole member, ``max|off-sector| / max|member|``, and its
      classification (``<= 1e-12``: numerical leakage from a non-canonical
      producer; above, or not finite: an incompatible Hamiltonian)."""
    if want not in LAYOUTS:
        raise ValueError("read_slot_member: want must be {}, got {!r}".format(_LAYOUT_SET, want))
    layout = _stored_layout(data)
    if layout not in LAYOUTS:
        raise ValueError("read_slot_member: '{}' has block_layout {}, expected {}; check it "
                         "with check_sector_map_members first".format(
                             file, _layout_text(data), _LAYOUT_SET))
    if layout == want == "dense":
        return _member(data, name, file)
    if not isinstance(sector_map, SectorMap):
        raise ValueError("read_slot_member: converting '{}' of '{}' needs the run's sector map, "
                         "got {}".format(name, file, type(sector_map).__name__
                                         if sector_map is not None else None))
    if layout == "sz_packed":
        stored = np.asarray(data["sz_sector_index"]) if "sz_sector_index" in data else None
        if stored is None or not np.array_equal(stored, sector_map.index):
            raise ValueError("read_slot_member: the sector map stored in '{}' differs from this "
                             "run's map (B = {}, nd = {}); check it with "
                             "check_sector_map_members first".format(file, sector_map.B,
                                                                     sector_map.nd))
        packed = PackedArrays(sector_map, tuple(
            _member(data, _packed_key(name, s), file) for s in sector_map.sector_order))
        if want == "sz_packed":
            return packed
        out = np.zeros(packed.lead_shape + (sector_map.ND, sector_map.ND), dtype=np.complex128)
        if not packed.lead_shape:
            return packed.to_dense(out=out)
        for l in range(packed.lead_shape[0]):
            PackedArrays(sector_map, tuple(b[l] for b in packed.blocks)).to_dense(out=out[l])
        return out
    dense = _member(data, name, file)
    ND = sector_map.ND
    if dense.ndim < 2 or dense.shape[-2:] != (ND, ND):
        raise ValueError("read_slot_member: '{}' of '{}' must end in shape (ND, ND) = ({}, {}), "
                         "got shape {}".format(name, file, ND, ND, dense.shape))
    if dense.ndim == 2:
        blocks, off_max = _sector_blocks(dense, sector_map)
        if off_max != 0.0:
            raise _leak_refusal(file, off_max, float(np.max(np.abs(dense))))
        return PackedArrays(sector_map, tuple(blocks))
    out = PackedArrays.zeros(sector_map, dense.shape[:-2])
    off_max, scale = 0.0, 0.0
    for l in range(dense.shape[0]):
        piece = dense[l]
        blocks, off = _sector_blocks(piece, sector_map)
        if off != 0.0:
            off_max = float("nan") if (off != off or off_max != off_max) else max(off_max, off)
        if off_max == 0.0:
            for b, src in zip(out.blocks, blocks):
                b[l] = src
        scale = float(np.max([scale, np.max(np.abs(piece))]))     # nan and inf propagate
    if off_max != 0.0:
        raise _leak_refusal(file, off_max, scale)
    return out


# ---------------------------------------------------------------------------
# The admission plan
# ---------------------------------------------------------------------------

_GIB = float(2 ** 30)
_ENTRY = 16                     # bytes per complex128 entry
#: the consumers of a plan, in canonical order
CONSUMERS = ("vertex", "chi_static", "chi_dynamic", "pairing_static")
#: what each consumer needs besides itself (closed transitively)
_IMPLIES = {"vertex": (), "chi_static": ("vertex",), "chi_dynamic": ("vertex",),
            "pairing_static": ("chi_static",)}
#: the layouts of a plan and the reason line of a dense plan whose caller made no selection
LAYOUTS = ("dense", "sz_packed")
_DEFAULT_LAYOUT_REASON = "so bond layout: dense (reason: not selected by the caller)"


class SOBondCaps(NamedTuple):
    """The memory caps of a plan: ``host_bytes`` (``None`` = unlimited), the
    raw user value ``host_raw`` (diagnostics only) and ``device_bytes``."""
    host_bytes: Optional[float]
    host_raw: Optional[float]
    device_bytes: Optional[float]


#: per helper (and bubble frequency mode): the observed fields compared with
#: the plan, in check order, and the rows the plan must list
_REQUIRE = {
    "so_bond_vertex": (("nvol", "B", "nd", "ND", "layout"), ("gamma_b",)),
    ("so_bond_bubble", "static"): (("nvol", "B", "nd", "ND", "layout", "nmat"),
                                   ("chibar_static", "bubble_static_temp", "bubble_pair_block")),
    ("so_bond_bubble", "all"): (("nvol", "B", "nd", "ND", "layout", "nmat", "placement", "slots"),
                                ("chibar", "bubble_pair_block")),
    "so_bond_dress": (("nvol", "ND", "batch"), ("dress_batch",)),
    "so_bond_pairing_vertex_parts": (("nvol", "ND", "B"),
                                     ("v_bare", "v_fluct", "pairing_q_workspace")),
    "so_bond_pairing_vertex_batch": (("nvol", "ND", "B", "layout", "batch"),
                                     ("vertex_batch_workspace",)),
}
#: optional observed fields (checked after the required ones) and their
#: value when the call does not pass them
_OPTIONAL = {"so_bond_dress": (("layout", "dense"),),
             "so_bond_pairing_vertex_parts": (("layout", "dense"),)}


@dataclass(frozen=True)
class SOBondPlan:
    """The admission plan of the spin-orbital bond helpers.

    Built once by :func:`plan_so_bond` before any large allocation; every
    helper that allocates a large array checks its own boundary against it
    (:meth:`require`) first.

    Fields
    ------
    nmat, nslices, nvol, B, nd, ND : int
        The Green function's Matsubara count (kept also by static plans),
        the bosonic slice count (``nmat``: no half grid), the k points, the
        bond channels, ``nd = 2 * norb_phys`` and ``ND = B * nd**2``.
    nb : int
        The frequency batch size.
    dress_capacity : int
        The largest batch the dressing may take: 1 when the closed consumer
        set has no ``chi_dynamic``, ``nb`` otherwise.
    consumers : tuple of str
        The transitively closed consumer set, canonical order.
    placement, layout, layout_reason, kernel_residency : str
        ``"host"``, ``"dense"`` or ``"sz_packed"`` (with the reason line of
        the selection) and ``"host"``.
    slots : tuple of str
        The store slots (``("chibar",)`` for ``chi_dynamic``, plus ``"chi"``
        with ``keep_chi``).
    aliases : tuple of (str, str)
        In-place slot overwrites (none).
    keep_chi, write_chibar : bool
    small_array_allowed, small_array_bytes : bool, float
        The small-array bound ``min(2**30, cap_host / 8)`` (``2**30`` when
        unlimited); no consumer of this plan uses the small-array helper,
        so it is not allowed.
    rows : tuple of (phase, name, bytes)
        One entry per logical allocation; phase ``"resident"`` lives
        through every phase.
    peak_host, peak_device : float, None
        ``max over phases of 1.25 * (resident + the phase's rows)``; no
        device rows.
    cap_host, cap_host_raw, cap_device : float or None
    reason : str
        The admission summary (logged).
    """

    nmat: int
    nslices: int
    nvol: int
    B: int
    nd: int
    ND: int
    nb: int
    dress_capacity: int
    consumers: tuple
    placement: str
    layout: str
    layout_reason: str
    kernel_residency: str
    slots: tuple
    aliases: tuple
    keep_chi: bool
    write_chibar: bool
    small_array_allowed: bool
    small_array_bytes: float
    rows: tuple
    peak_host: float
    peak_device: Optional[float]
    cap_host: Optional[float]
    cap_host_raw: Optional[float]
    cap_device: Optional[float]
    reason: str

    def require(self, helper, **observed):
        """Check a helper's boundary against the plan before it allocates.

        ``helper`` is ``"so_bond_vertex"``, ``"so_bond_bubble"`` (with the
        observed ``freq``, ``"static"`` or ``"all"``), ``"so_bond_dress"``,
        ``"so_bond_pairing_vertex_parts"`` or
        ``"so_bond_pairing_vertex_batch"``. The observed fields are the
        dimensions visible at that boundary:

        * vertex: ``nvol``, ``B``, ``nd``, ``ND``, ``layout``;
        * static bubble: those and ``nmat``; every-slice bubble: those,
          ``nmat``, ``placement`` and the store ``slots`` (which must hold
          ``"chibar"``);
        * dressing: ``nvol``, ``ND`` and the ``batch`` length
          (``1 <= batch <= dress_capacity``);
        * pairing parts: ``nvol``, ``ND`` and the geometry's ``B`` (with
          ``ND == plan.ND = plan.B * plan.nd**2`` this is
          ``ND == B * plan.nd**2``);
        * the dressing and the pairing parts also observe ``layout``, which
          a call may omit: an omitted ``layout`` is ``"dense"``;
        * pairing vertex batch: ``nvol``, ``ND``, ``B``, ``layout`` and the
          ``batch`` length (the plan needs ``vertex_batch=True``).

        Raises ``ValueError`` naming the first mismatching field (or the
        first missing row) and :attr:`reason`."""
        key = helper
        if helper == "so_bond_bubble":
            freq = observed.pop("freq", None)
            key = (helper, freq)
            if key not in _REQUIRE:
                raise ValueError("SOBondPlan.require: so_bond_bubble needs freq='static' or "
                                 "'all', got {!r}".format(freq))
        if key not in _REQUIRE:
            raise ValueError("SOBondPlan.require: unknown helper {!r}".format(helper))
        fields, needed = _REQUIRE[key]
        given = sorted(observed)
        optional = _OPTIONAL.get(key, ())
        for f, default in optional:
            observed.setdefault(f, default)
        shown = fields + tuple("{} (optional)".format(f) for f, _ in optional)
        fields = fields + tuple(f for f, _ in optional)
        if set(observed) != set(fields):
            raise ValueError("SOBondPlan.require: {} observes the fields {}, got {}".format(
                helper, ", ".join(shown), ", ".join(given)))
        for f in fields:
            got = observed[f]
            if f == "batch":
                if not _is_int(got) or int(got) <= 0:
                    raise ValueError("{}: the batch length must be a positive integer, got {!r}"
                                     .format(helper, got))
                if int(got) > self.dress_capacity:
                    if helper == "so_bond_pairing_vertex_batch":
                        raise ValueError("{}: a batch of {} slices exceeds the plan's batch "
                                         "capacity (dress_capacity = {}) ({})".format(
                                             helper, int(got), self.dress_capacity,
                                             self.reason))
                    raise ValueError("{}: a batch of {} slices exceeds the plan's "
                                     "dress_capacity = {} ({})".format(
                                         helper, int(got), self.dress_capacity, self.reason))
            elif f == "slots":
                if "chibar" not in tuple(got):
                    raise ValueError("{}: the store slots {} do not include 'chibar' ({})"
                                     .format(helper, tuple(got), self.reason))
            elif got != getattr(self, f):
                raise ValueError("{}: the admission plan has {} = {} but the call has {} = {} "
                                 "({})".format(helper, f, getattr(self, f), f, got, self.reason))
        names = {name for _, name, _ in self.rows}
        for name in needed:
            if name not in names:
                raise ValueError("{}: the admission plan has no row '{}' (closed consumers: {}; "
                                 "{})".format(helper, name, ", ".join(self.consumers),
                                              self.reason))


def _positive_int(name, value, *, even=False, what=""):
    if not _is_int(value) or int(value) <= 0 or (even and int(value) % 2):
        raise ValueError("plan_so_bond: {} must be a positive {}integer{}, got {!r}".format(
            name, "even " if even else "", what, value))
    return int(value)


def _close_consumers(consumers):
    if isinstance(consumers, str):
        consumers = (consumers,)
    given = tuple(consumers)
    if not given:
        raise ValueError("plan_so_bond: consumers must name at least one of {}"
                         .format(", ".join(CONSUMERS)))
    closed = set()
    stack = []
    for c in given:
        if c not in _IMPLIES:
            raise ValueError("plan_so_bond: unknown consumer {!r}; the consumers are {}"
                             .format(c, ", ".join(CONSUMERS)))
        stack.append(c)
    while stack:
        c = stack.pop()
        if c not in closed:
            closed.add(c)
            stack.extend(_IMPLIES[c])
    return tuple(c for c in CONSUMERS if c in closed)


def _layout_entries(layout, B, nd):
    """The complex entries of one ``ND x ND`` code-slot matrix in ``layout``:
    ``ND**2`` dense, ``sum_sigma n_sigma**2 = (6/16) ND**2`` S_z-packed
    (``n_0 = 2 B norb**2``, ``n_{+1} = n_{-1} = B norb**2``)."""
    if layout == "dense":
        return (B * nd * nd) ** 2
    norb = nd // 2
    n0, n1 = 2 * B * norb * norb, B * norb * norb
    return n0 * n0 + 2 * n1 * n1


def _plan_rows(consumers, *, nmat, nvol, nd, ND, dress_capacity, keep_chi, layout="dense",
               vertex_batch=False):
    """The row table: ``(phase, name, bytes)``, one entry per logical
    allocation, ``"resident"`` first, then each closed consumer's phase.

    Every row of whole code-slot matrices holds ``S`` entries per matrix
    (:func:`_layout_entries`); the bubble's pair block and the static
    bubble temporary (always a dense ``(nvol, ND, ND)`` array) do not
    depend on the layout."""
    S = _layout_entries(layout, ND // (nd * nd), nd)
    full = _ENTRY * nvol * S                                # one (nvol, ND, ND) object
    dense_full = _ENTRY * nvol * ND * ND                    # one dense (nvol, ND, ND) array
    P = nd * nd
    # the bond bubble's per-shift workspace per k point: its P^2-sized and
    # P-sized frequency-grid buffers (bubble.BOND_BUBBLE_N4_BUFFERS and
    # BOND_BUBBLE_N2_BUFFERS, the counts the bubble's memory contract and
    # sc's estimator use) plus 8 P per-k temporaries
    pair_block = _ENTRY * nvol * (bubble.BOND_BUBBLE_N4_BUFFERS * nmat * P * P
                                  + bubble.BOND_BUBBLE_N2_BUFFERS * nmat * P + 8 * P)
    dress = 7 * _ENTRY * dress_capacity * nvol * S
    rows = [("resident", "gamma_b", full)]
    if "chi_static" in consumers:
        rows += [("chi_static", "chibar_static", full),
                 ("chi_static", "bubble_static_temp", dense_full),
                 ("chi_static", "bubble_pair_block", pair_block),
                 ("chi_static", "dress_batch", dress), ("chi_static", "chi_static", full)]
    if "chi_dynamic" in consumers:
        rows += [("chi_dynamic", "chibar", _ENTRY * nmat * nvol * S)]
        if keep_chi:
            rows += [("chi_dynamic", "chi", _ENTRY * nmat * nvol * S)]
        rows += [("chi_dynamic", "bubble_pair_block", pair_block),
                 ("chi_dynamic", "dress_batch", dress)]
        if vertex_batch:
            rows += [("chi_dynamic", "vertex_batch_workspace", 2 * _ENTRY * S)]
    if "pairing_static" in consumers:
        rows += [("pairing_static", "chibar_static", full),
                 ("pairing_static", "bubble_static_temp", dense_full),
                 ("pairing_static", "bubble_pair_block", pair_block),
                 ("pairing_static", "chi_static", full),
                 ("pairing_static", "dress_batch", dress), ("pairing_static", "v_bare", full),
                 ("pairing_static", "v_fluct", full),
                 ("pairing_static", "pairing_q_workspace", 2 * _ENTRY * S)]
    return tuple((ph, name, int(b)) for ph, name, b in rows)


def _row_table(rows):
    return "rows (phase  name  bytes):\n" + "\n".join(
        "  {}  {}  {}".format(ph, name, b) for ph, name, b in rows)


def plan_so_bond(*, nmat, nvol, B, nd, half_grid=False, consumers, caps,
                 placement_request="host", nb=None, keep_chi=False, write_chibar=False,
                 layout="dense", layout_reason=None, vertex_batch=False):
    """The admission plan of the spin-orbital bond helpers (host placement).

    Parameters
    ----------
    nmat : int
        The Green function's Matsubara count, a positive even integer
        (static plans keep it: the bubble's per-shift pair block scales with
        it).
    nvol, B, nd : int
        k points, bond channels and ``nd = 2 * norb_phys``.
    consumers : iterable of str
        A subset of ``("vertex", "chi_static", "chi_dynamic",
        "pairing_static")``, closed transitively (``pairing_static =>
        chi_static => vertex``, ``chi_dynamic => vertex``).
    caps : SOBondCaps
        ``host_bytes = None`` means unlimited: the peak is only logged; a
        negative or non-finite ``host_bytes`` is refused.
    nb : int, optional
        The frequency batch, default ``min(nmat, 16)``.
    layout : str
        ``"dense"`` (default) or ``"sz_packed"`` (the S_z-sector packed
        layout of the code-slot objects).
    layout_reason : str, optional
        The reason line of the layout selection (``select_layout``); it
        must start with ``"so bond layout: <layout> ("``. ``None`` is
        accepted for the dense layout only and stands for ``"so bond
        layout: dense (reason: not selected by the caller)"``.
    vertex_batch : bool
        Add the row of the per-frequency pairing vertex batch
        (``vertex_batch_workspace = 2 * 16 S`` in the ``chi_dynamic``
        phase); needs the consumer ``chi_dynamic``.

    Rows (16 bytes per complex entry): ``gamma_b = 16 nvol ND^2``
    (resident); the per-shift ``bubble_pair_block = 16 nvol (3 nmat P^2 +
    3 nmat P + 8 P)`` (``P = nd^2``; the 3s are
    ``bubble.BOND_BUBBLE_N4_BUFFERS`` and ``BOND_BUBBLE_N2_BUFFERS``) in
    every bubble or pairing phase; ``chi_static``: the static bubble
    output and its temporary, ``chi_static``; ``chi_dynamic``: the store
    slot ``chibar = 16 nmat nvol ND^2`` (and ``chi`` with ``keep_chi``);
    ``dress_batch = 7 * 16 dress_capacity nvol ND^2`` in every phase that
    dresses; ``pairing_static``: the static bubble output and its temporary,
    the pair block, ``chi_static``, the dressing, ``v_bare``, ``v_fluct``
    (``16 nvol ND^2`` each) and ``pairing_q_workspace = 2 * 16 ND^2``. The
    peak is the maximum over the phases of ``1.25 (resident + the phase's
    rows)``; a peak above ``caps.host_bytes`` is refused with the row table.
    The refusal names the solver as ``Eliashberg`` when the consumers
    include ``pairing_static`` and ``FLEX`` otherwise.

    The ``ND^2`` above is the entry count ``S`` of one ``ND x ND`` code-slot
    matrix: ``S = ND^2`` dense and ``S = sum_sigma n_sigma^2 = (6/16) ND^2``
    S_z-packed, in every row except ``bubble_pair_block`` and
    ``bubble_static_temp`` (dense in both layouts).
    """
    if half_grid:
        raise ValueError("plan_so_bond: the half bosonic grid is not available on the "
                         "spin-orbital bond path")
    if placement_request != "host":
        raise ValueError("plan_so_bond: placement_request must be 'host' on the spin-orbital "
                         "bond path, got {!r}".format(placement_request))
    if not isinstance(caps, SOBondCaps):
        raise ValueError("plan_so_bond: caps must be a SOBondCaps, got {!r}".format(caps))
    if not (isinstance(layout, str) and layout in LAYOUTS):
        raise ValueError("plan_so_bond: layout must be 'dense' or 'sz_packed', got {!r}"
                         .format(layout))
    if layout_reason is None:
        if layout != "dense":
            raise ValueError("plan_so_bond: layout = {} needs its layout_reason (the reason line "
                             "of the layout selection)".format(layout))
        layout_reason = _DEFAULT_LAYOUT_REASON
    head = "so bond layout: {} (".format(layout)
    if not (isinstance(layout_reason, str) and layout_reason.startswith(head)):
        raise ValueError("plan_so_bond: layout_reason must be a string starting with {!r}, "
                         "got {!r}".format(head, layout_reason))
    if caps.host_bytes is not None:
        try:
            cap_ok = float(caps.host_bytes) >= 0.0 and np.isfinite(float(caps.host_bytes))
        except (TypeError, ValueError):
            cap_ok = False
        if not cap_ok:
            raise ValueError("plan_so_bond: caps.host_bytes must be None (unlimited) or a "
                             "finite number of bytes >= 0, got {!r}".format(caps.host_bytes))
    nmat = _positive_int("nmat", nmat, even=True)
    nvol = _positive_int("nvol", nvol)
    B = _positive_int("B", B)
    nd = _positive_int("nd", nd, even=True, what=" (2·norb_phys)")
    nb = min(nmat, 16) if nb is None else _positive_int("nb", nb)
    closed = _close_consumers(consumers)
    vertex_batch = bool(vertex_batch)
    if vertex_batch and "chi_dynamic" not in closed:
        raise ValueError("plan_so_bond: vertex_batch needs the consumer chi_dynamic (the vertex "
                         "batch reads dressed frequency batches)")
    ND = B * nd * nd
    dress_capacity = nb if "chi_dynamic" in closed else 1
    keep_chi, write_chibar = bool(keep_chi), bool(write_chibar)
    rows = _plan_rows(closed, nmat=nmat, nvol=nvol, nd=nd, ND=ND,
                      dress_capacity=dress_capacity, keep_chi=keep_chi, layout=layout,
                      vertex_batch=vertex_batch)
    resident = sum(b for ph, _, b in rows if ph == "resident")
    peaks = [(1.25 * (resident + sum(b for ph, _, b in rows if ph == c)), c) for c in closed]
    peak, phase = max(peaks, key=lambda p: p[0])
    slots = (("chibar",) + (("chi",) if keep_chi else ())) if "chi_dynamic" in closed else ()
    cap = caps.host_bytes
    small = _GIB if cap is None else min(_GIB, float(cap) / 8.0)
    reason = ("spin-orbital bond plan: consumers {}; {} layout, host placement; nmat = {}, "
              "nvol = {}, B = {}, nd = {}, ND = {}, nb = {}, dress_capacity = {}; estimated "
              "peak {:.3f} GiB in phase {}; host cap {}".format(
                  ", ".join(closed), layout, nmat, nvol, B, nd, ND, nb, dress_capacity,
                  peak / _GIB,
                  phase, "unlimited" if cap is None else "{:.3f} GiB".format(cap / _GIB)))
    if cap is not None and peak > cap:
        solver = "Eliashberg" if "pairing_static" in closed else "FLEX"
        slot = ("16·Nmat·nvol·ND² bytes" if layout == "dense" else
                "16·Nmat·nvol·Σn_σ² bytes (Σn_σ² = (6/16)·ND²)")
        raise ValueError(
            "spin-orbital bond {} needs an estimated peak of {:.2f} GiB in phase {} (B = {} "
            "channels, ND = B·(2·norb_phys)² = {}, nvol = {}, Nmat = {}; one store slot = "
            "{} = {:.2f} GiB) but the cap is {:.2f} GiB (caps.host_bytes). "
            "Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce Nmat or the k grid, or prune "
            "bond shells.\n{}".format(
                solver, peak / _GIB, phase, B, ND, nvol, nmat, slot,
                _ENTRY * nmat * nvol * _layout_entries(layout, B, nd) / _GIB, cap / _GIB,
                _row_table(rows)))
    logger.info(reason)
    return SOBondPlan(
        nmat=nmat, nslices=nmat, nvol=nvol, B=B, nd=nd, ND=ND, nb=nb,
        dress_capacity=dress_capacity, consumers=closed, placement="host", layout=layout,
        layout_reason=layout_reason, kernel_residency="host", slots=slots, aliases=(),
        keep_chi=keep_chi, write_chibar=write_chibar, small_array_allowed=False,
        small_array_bytes=small, rows=rows, peak_host=peak, peak_device=None,
        cap_host=None if cap is None else float(cap), cap_host_raw=caps.host_raw,
        cap_device=caps.device_bytes, reason=reason)


# ---------------------------------------------------------------------------
# Memory estimators of the S_z-sector eigen solve
# ---------------------------------------------------------------------------

def _int_at_least(helper, name, value, low):
    if not _is_int(value) or int(value) < low:
        raise ValueError("{}: {} must be an integer >= {}, got {!r}".format(
            helper, name, low, value))
    return int(value)


def _int_list(value):
    """``value`` as a list of Python integers, or ``None`` when it is not a
    one-dimensional sequence of integers."""
    if isinstance(value, np.ndarray):
        if value.ndim != 1:
            return None
        value = value.tolist() if np.issubdtype(value.dtype, np.integer) else list(value)
    if not isinstance(value, (list, tuple)) or not all(_is_int(v) for v in value):
        return None
    return [int(v) for v in value]


#: the Python-level bookkeeping of one call (array headers, index tuples,
#: result records): a fixed allowance that does not scale with the data
_CALL_OVERHEAD = 16 * 1024
#: the entries of one iterator buffer numpy < 2 allocates per operand of a
#: broadcast ``einsum`` contraction (``numpy.BUFSIZE``)
_EINSUM_BUFFER = 8192


def estimate_sector_solver_bytes(dims, k, n_full_gap, *, n_sector_vec):
    """The host bytes of the per-sector eigen solve of the S_z-packed pairing
    kernel, for ``k`` kept eigenpairs.

    ``dims[p]`` is the reduced dimension of sector ``p`` and
    ``n_sector_vec[p]`` the length of its sector vector (``0 <= dims[p] <=
    n_sector_vec[p] / 2``); ``n_full_gap >= sum(n_sector_vec)`` is the
    length of one dense gap. Per sector (16 bytes per complex entry,
    ``n_p = n_sector_vec[p]``):

    * dense branch, ``dim <= max(2k + 1, 20)``: ``16 (7 dim^2 + 3 n_p dim +
      4 dim)`` (the identity input, the copy ``np.linalg.eig`` makes, its
      eigenvectors, the reduced matrix and the 2-norm's copy of it; the
      lifted input, the operator's column-by-column image and its stacked
      form, ``n_p x dim`` each);
    * ARPACK branch: ``16 (ncv dim + 3 ncv^2 + 14 ncv + 3 dim) + 24 n_p +
      6 * 16 n_p + 16 * 3 k (dim + n_p)`` with ``ncv = min(max(2k + 1, 20),
      dim)`` (the Krylov basis, ``workl = 3 ncv^2 + 5 ncv``, ``workev``,
      ``rwork``, the Ritz values, ``workd``, ``resid``; the index arrays and
      the sector-vector temporaries; the returned eigenvectors, those the
      largest-real re-solve keeps from the first solve and the residual
      images of the ``k`` pairs: three reduced and three sector-vector
      columns per pair);
    * null operator: ``16 (min(k, dim) (dim + n_p) + n_full_gap)`` (unit
      vectors, their lifts and one scattered gap).

    A sector counts the larger of its branch and the null-operator row
    (whether an operator vanishes is known only at run time). The result is
    the maximum over sectors plus ``16 (k + sum_p min(k, dim_p))
    n_full_gap`` (the candidate gaps of every sector, alive until the merge,
    and the stacked kept gaps), ``16 sum_p min(k, dim_p) dim_p`` (the reduced
    eigenvectors) and :data:`_CALL_OVERHEAD`. Arguments outside these ranges
    are refused (``ValueError`` naming the argument)."""
    name = "estimate_sector_solver_bytes"
    d = _int_list(dims)
    if not d or any(v < 0 for v in d):
        raise ValueError("{}: dims must be a non-empty sequence of integers >= 0, got {!r}"
                         .format(name, dims))
    n = _int_list(n_sector_vec)
    if n is None or len(n) != len(d) or any(v < 0 for v in n):
        raise ValueError("{}: n_sector_vec must be a sequence of {} integers >= 0 (one per "
                         "dims entry), got {!r}".format(name, len(d), n_sector_vec))
    for p, (dim, n_p) in enumerate(zip(d, n)):
        if 2 * dim > n_p:
            raise ValueError("{}: dims[{}] = {} exceeds n_sector_vec[{}] / 2 = {} / 2".format(
                name, p, dim, p, n_p))
    k = _int_at_least(name, "k", k, 1)
    total_n = sum(n)
    if not _is_int(n_full_gap) or int(n_full_gap) < total_n:
        raise ValueError("{}: n_full_gap must be an integer >= sum(n_sector_vec) = {}, got {!r}"
                         .format(name, total_n, n_full_gap))
    n_full_gap = int(n_full_gap)
    cut = max(2 * k + 1, 20)
    worst = 0
    for dim, n_p in zip(d, n):
        if dim <= cut:
            branch = _ENTRY * (7 * dim * dim + 3 * n_p * dim + 4 * dim)
        else:
            ncv = min(cut, dim)
            branch = (_ENTRY * (ncv * dim + 3 * ncv * ncv + 14 * ncv + 3 * dim)
                      + 24 * n_p + 6 * _ENTRY * n_p + _ENTRY * 3 * k * (dim + n_p))
        null = _ENTRY * (min(k, dim) * (dim + n_p) + n_full_gap)
        worst = max(worst, branch, null)
    kept = sum(min(k, dim) for dim in d)
    return int(worst + _ENTRY * (k + kept) * n_full_gap
               + _ENTRY * sum(min(k, dim) * dim for dim in d) + _CALL_OVERHEAD)


def estimate_packed_kernel_workspace_bytes(B, pp_pair_counts, nvol, nfreq):
    """The host bytes of the S_z-packed pairing kernel's workspace.

    ``pp_pair_counts`` are the gap pairs ``n_p`` of the three pp sectors,
    ``B`` the bond channels, ``nvol`` the k points and ``nfreq`` the
    fermionic frequencies (all positive integers; ``ValueError`` naming the
    argument otherwise). With ``grid = nvol nfreq`` and ``n = max_p n_p``,
    returns a dict (16 bytes per complex entry):

    * ``hoisted_blocks = 16 B^2 sum_p n_p^2 grid`` (the vertex blocks of
      every pp sector, hoisted once);
    * ``g2_sectors = 16 sum_p n_p^2 grid`` (the pair propagator per pp
      sector);
    * ``hoist_workspace = 16 * 6 n^2 grid`` (the hoisting of one channel
      pair: the raw block, its boson-to-tau image, the inverse spatial FFT's
      output and per-axis intermediate, the FFT's internal buffers and, on
      numpy < 2, the contiguous copy its FFT makes of a swapped axis);
    * ``g2_gather = 16 n^2 grid`` (the pair propagator of one sector as
      gathered, before its contiguous copy);
    * ``matvec_workspace = 16 (max(3 + B, 6) n grid + 3 * 8192)`` (one
      matrix-vector product: the ``B`` output accumulators, the transformed
      input, its rolled copy and one leg product; or the transform stage --
      the input, the fermion-to-tau image and its phased copy, the spatial
      FFT's output and intermediate and an axis-swap copy; plus the iterator
      buffers of an ``einsum`` contraction on numpy < 2, one per operand);
    * ``call_overhead`` (:data:`_CALL_OVERHEAD`).

    The constructor needs ``hoisted_blocks + g2_sectors + hoist_workspace +
    g2_gather + call_overhead``; a product on the built kernel needs
    ``matvec_workspace + call_overhead`` on top of the kept blocks.
    """
    name = "estimate_packed_kernel_workspace_bytes"
    B = _int_at_least(name, "B", B, 1)
    counts = _int_list(pp_pair_counts)
    if counts is None or len(counts) != 3 or any(c < 1 for c in counts):
        raise ValueError("{}: pp_pair_counts must be three positive integers, got {!r}".format(
            name, pp_pair_counts))
    nvol = _int_at_least(name, "nvol", nvol, 1)
    nfreq = _int_at_least(name, "nfreq", nfreq, 1)
    sq = sum(c * c for c in counts)
    grid = nvol * nfreq
    n = max(counts)
    return {"hoisted_blocks": _ENTRY * B * B * sq * grid,
            "g2_sectors": _ENTRY * sq * grid,
            "hoist_workspace": _ENTRY * 6 * n * n * grid,
            "g2_gather": _ENTRY * n * n * grid,
            "matvec_workspace": _ENTRY * (max(3 + B, 6) * n * grid + 3 * _EINSUM_BUFFER),
            "call_overhead": _CALL_OVERHEAD}


# ---------------------------------------------------------------------------
# Gamma_B assembly and its exchange image
# ---------------------------------------------------------------------------

def _check_convention(helper, geom):
    conv = getattr(geom, "convention", None)
    if conv != SLOT_CONVENTION:
        shown = conv.to_string() if isinstance(conv, SlotConvention) else repr(conv)
        raise ValueError("{}: slot convention {} is not the pinned {}".format(
            helper, shown, SLOT_CONVENTION.to_string()))


def _nd_of(shape, geom, where, name):
    """``nd`` of a ``(nvol, B nd^2, B nd^2)`` array on ``geom``."""
    nvol, B = geom.nvol, geom.B
    nd = 0
    if len(shape) == 3 and shape[0] == nvol and shape[1] == shape[2] and shape[1] % B == 0:
        nd = int(round(np.sqrt(shape[1] // B)))
        if nd * nd * B != shape[1] or nd % 2 or nd == 0:
            nd = 0
    if nd == 0:
        raise ValueError("{}: {} has shape {}; the geometry needs ({}, B*nd**2, B*nd**2) with "
                         "B = {} and nvol = {}".format(where, name, tuple(shape), nvol, B, nvol))
    return nd


def _sector_name(sigma):
    return SECTOR_ORDER[_SIGMAS.index(int(sigma))]


def _coupling_message(where, typ, s_row, s_col):
    return ("{}: a {} record couples the S_z sectors {} and {} in the code slots; the packed "
            "layout needs an S_z-conserving interaction (select_layout)".format(
                where, typ, _sector_name(s_row), _sector_name(s_col)))


def _first_sector_coupling(records, norb_phys):
    """``(type, sigma_row, sigma_col)`` of the first record whose direct
    filing ``[(.; 1, 2), (.; 4, 3)]`` joins two ph sectors (``S_z(1) -
    S_z(2) != S_z(4) - S_z(3)``; its crossed filing then joins two as
    well), or ``None``."""
    norb = int(norb_phys)
    for r in records:
        s1, s2, s3, s4 = (1 - 2 * (int(x) // norb) for x in r["modes"])   # 2 S_z
        if s1 - s2 != s4 - s3:
            return str(r["type"]), (s1 - s2) // 2, (s4 - s3) // 2
    return None


def _file_records(out, records, geom, *, parts=("direct", "crossed")):
    """Add the records into the code-slot array ``out`` ``(nvol, ND, ND)``.

    A record ``v c+_1(0) c_2(d1) c+_3(t) c_4(t+d2)`` is filed in label form
    (each part alone gives ``H_int = (1/2N) sum_q Gamma_IJ(q) rho_I(q)
    rho_J(-q)``, so direct plus crossed gives it with ``1/4N`` up to one-body
    terms; labels
    ``(d; X, Y)`` of ``rho^{(d)}_{XY}(q) = sum_j e^{+iq.j} c+_{X,j}
    c_{Y,j+d}``)::

        direct:   Gamma[(d1;1,2),(d2;3,4)](q) += v e^{+iq.t},
                  Gamma[(d2;3,4),(d1;1,2)](q) += v e^{-iq.t}
        crossed:  Gamma[(t+d2;1,4),(d1-t;3,2)](q) += -v e^{+iq.t},
                  Gamma[(d1-t;3,2),(t+d2;1,4)](q) += -v e^{-iq.t}

    (the crossed part from ``c+_1 c_2 c+_3 c_4 = -c+_1 c_4 c+_3 c_2 +
    one-body``; the one-body remainders are not part of the vertex). Each
    label entry ``Gamma[(dA; a,b),(dB; c,e)]`` goes to the code slot
    ``[(m; a,b),(m'; e,c)]`` with ``dr_m = -dA``, ``dr_m' = dB`` and the
    phase ``e^{+iq.dr_m'}`` (``geom.phase[:, m']``), i.e.
    ``Gamma_B(q)[(m;x,y),(m';z,w)] = e^{+iq.dr_m'} Gamma(q)[(-dr_m; x,y),
    (+dr_m'; w,z)]``. ``parts`` selects the direct and/or crossed filing. A
    label displacement that is not a channel of ``geom`` (no wrapping) is
    refused.

    ``out`` may also be :class:`PackedArrays` (leading shape ``(nvol,)``):
    each entry is then added at its local positions in the block of its ph
    sector, the map read as the consumer ``"gamma_packer"``; an entry whose
    row and column lie in different sectors is refused (the record does not
    conserve S_z)."""
    _check_convention("_file_records", geom)
    parts = tuple(parts)
    if not parts or set(parts) - {"direct", "crossed"}:
        raise ValueError("_file_records: parts must be a non-empty subset of ('direct', "
                         "'crossed'), got {!r}".format(parts))
    packed = isinstance(out, PackedArrays)
    if packed:
        sm = _consumer_map("gamma_packer", out.sector_map)
        if out.lead_shape != (geom.nvol,) or sm.B != geom.B:
            raise ValueError("_file_records: the packed out has leading shape {} and B = {}; "
                             "the geometry needs ({},) and B = {}".format(
                                 out.lead_shape, sm.B, geom.nvol, geom.B))
        nd = int(sm.nd)
        sector_of, local = sm.sector_of, sm.local
        block_of = {int(sigma): out.blocks[s] for s, sigma in enumerate(sm.sigmas)}
    else:
        if not isinstance(out, np.ndarray) or out.dtype != np.complex128:
            raise ValueError("_file_records: out must be a complex128 ndarray")
        nd = _nd_of(out.shape, geom, "_file_records", "out")
    npair = nd * nd
    chan = {tuple(int(x) for x in d): m for m, d in enumerate(geom.delta_r)}
    rev = geom.reverse
    two_pi_q = 2.0 * np.pi * geom.q_frac
    direct, crossed = "direct" in parts, "crossed" in parts

    def channel(d, typ):
        m = chan.get(tuple(int(x) for x in d))
        if m is None:
            raise ValueError("_file_records: label displacement {} of a {} record is not a "
                             "channel of the bond geometry".format(
                                 tuple(int(x) for x in d), typ))
        return m

    for r in records:
        typ = str(r["type"])
        v = complex(r["v"])
        x1, x2, x3, x4 = (int(x) for x in r["modes"])
        if min(x1, x2, x3, x4) < 0 or max(x1, x2, x3, x4) >= nd:
            raise ValueError("_file_records: a {} record has modes {} outside [0, {})"
                             .format(typ, (x1, x2, x3, x4), nd))
        cells = np.asarray(r["cells"], dtype=np.int64)
        if cells[0].any():
            raise ValueError("_file_records: operator 1 of a {} record must sit in cell 0, "
                             "got {}".format(typ, tuple(int(x) for x in cells[0])))
        d1, t, d2 = cells[1], cells[2], cells[3] - cells[2]
        et = np.exp(1j * (two_pi_q @ t.astype(np.float64)))          # e^{+iq.t}
        entries = []
        if direct:
            entries += [(d1, x1, x2, d2, x3, x4, v * et), (d2, x3, x4, d1, x1, x2, v * et.conj())]
        if crossed:
            entries += [(t + d2, x1, x4, d1 - t, x3, x2, -v * et),
                        (d1 - t, x3, x2, t + d2, x1, x4, -v * et.conj())]
        for dA, a, b, dB, c, e, val in entries:
            m_row = int(rev[channel(dA, typ)])                         # dr_m = -dA
            m_col = channel(dB, typ)                                   # dr_m' = +dB
            if packed:
                row, col = m_row * npair + a * nd + b, m_col * npair + e * nd + c
                s_row, s_col = int(sector_of[row]), int(sector_of[col])
                if s_row != s_col:
                    raise ValueError(_coupling_message("_file_records", typ, s_row, s_col))
                block_of[s_row][:, local[row], local[col]] += val * geom.phase[:, m_col]
                continue
            out[:, m_row * npair + a * nd + b, m_col * npair + e * nd + c] += \
                val * geom.phase[:, m_col]


def so_bond_vertex(compiled, topo, spatial_shape, *, plan):
    """The spinful bond vertex ``Gamma_B`` ``(nvol, ND, ND)`` complex128.

    ``ND = B * nd**2`` (``nd = 2 * norb_phys``), slot ``(m; x, y)`` at
    ``m * nd**2 + x * nd + y`` in the spin-block mode index. Every record of
    ``compiled.records`` is filed by :func:`_file_records`: direct plus
    crossed (the antisymmetrized particle-hole vertex), in the code slots of
    :data:`SLOT_CONVENTION`, with ``H_int = (1/4N) sum_q sum_rc
    Gamma_B(q)_rc O_r(q) O_c(q)^dagger`` up to one-body terms (the
    antisymmetrized vertex counts each pair twice).

    The geometry is built by :func:`slot_geometry` (which refuses channels
    that coincide modulo the mesh), and the boundary is checked against
    ``plan`` (:meth:`SOBondPlan.require`) before ``Gamma_B`` is allocated.

    With ``plan.layout == "sz_packed"`` the result is :class:`PackedArrays`
    ``(nvol, n_sigma, n_sigma)`` on :func:`sz_sector_map`: the diagonal ph
    sector blocks of the same ``Gamma_B`` (an S_z-conserving interaction
    has no entry between two sectors). A record that couples two sectors
    is refused before the blocks are allocated.
    """
    # checked before the geometry is built; the adapter checks again for its
    # own callers (two isinstance tests)
    _check_vertex_args(compiled, plan)
    return _vertex_on_geom(compiled, slot_geometry(topo, spatial_shape), plan=plan)


def _check_vertex_args(compiled, plan):
    if not isinstance(compiled, CompiledInteraction):
        raise ValueError("so_bond_vertex: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got {}".format(
                             type(compiled).__name__))
    if not isinstance(plan, SOBondPlan):
        raise ValueError("so_bond_vertex: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))


def _vertex_on_geom(compiled, geom, *, plan):
    """:func:`so_bond_vertex` on a built geometry (same checks and result)."""
    _check_vertex_args(compiled, plan)
    _check_convention("so_bond_vertex", geom)
    nd = 2 * int(compiled.norb_phys)
    ND = geom.B * nd * nd
    # the vertex layout is the plan's own (the require call checks the
    # other observed fields); the guard against an interaction that couples
    # the S_z sectors is _first_sector_coupling below
    layout = plan.layout
    plan.require("so_bond_vertex", nvol=geom.nvol, B=geom.B, nd=nd, ND=ND, layout=layout)
    if layout == "sz_packed":
        bad = _first_sector_coupling(compiled.records, compiled.norb_phys)
        if bad is not None:
            raise ValueError(_coupling_message("so_bond_vertex", *bad))
        out = PackedArrays.zeros(sz_sector_map(geom, nd), (geom.nvol,))
    else:
        out = np.zeros((geom.nvol, ND, ND), dtype=np.complex128)
    _file_records(out, compiled.records, geom)
    return out


def _negate_index(shape):
    """``neg[q]``: the C-order index of ``-q`` on the mesh."""
    grid = np.indices(shape).reshape(3, -1)
    return np.ravel_multi_index(tuple(np.mod(-grid, np.array(shape)[:, None])), shape)


def gamma_p(gamma_b, geom, *, out):
    """The exchange image of a code-slot vertex, written into ``out``:

    ``Gamma^P(q)[(m;x,y),(m';z,w)] = e^{+iq.(dr_m' - dr_m)}
    Gamma_B(-q)[(rev m'; w,z), (rev m; y,x)]``

    (particle exchange ``Gamma_IJ(q) = Gamma_JI(-q)`` of the label form,
    written in the code slots). A vertex filed by :func:`so_bond_vertex`
    equals its exchange image. ``out`` is a caller-owned complex128 array of
    the shape of ``gamma_b`` that must not share memory with it; the work
    is done per q in two ``ND x ND`` matrices. Returns ``out``.

    ``gamma_b`` and ``out`` may both be :class:`PackedArrays` on the sector
    map of ``geom`` (leading shape ``(nvol,)``): the slot map sends sector
    ``sigma`` onto ``-sigma``, so block ``sigma`` of the image is read from
    block ``-sigma`` of ``gamma_b`` through the exchange-image permutation
    of the map (consumer ``"gamma_p"``), transposed and phased; the work is
    done per q in two ``n_sigma x n_sigma`` matrices."""
    if isinstance(gamma_b, PackedArrays) or isinstance(out, PackedArrays):
        return _gamma_p_packed(gamma_b, geom, out=out)
    if np.shares_memory(out, gamma_b):
        raise ValueError("gamma_p: out shares memory with gamma_b")
    _check_convention("gamma_p", geom)
    if not isinstance(gamma_b, np.ndarray) or gamma_b.dtype != np.complex128:
        raise ValueError("gamma_p: gamma_b must be a complex128 ndarray")
    nd = _nd_of(gamma_b.shape, geom, "gamma_p", "gamma_b")
    if not isinstance(out, np.ndarray) or out.dtype != np.complex128 \
            or out.shape != gamma_b.shape:
        raise ValueError("gamma_p: out must be a complex128 array of shape {}, got {} {}".format(
            gamma_b.shape, getattr(out, "dtype", type(out).__name__),
            tuple(getattr(out, "shape", ()))))
    nvol, ND = gamma_b.shape[0], gamma_b.shape[1]
    npair = nd * nd
    idx = np.arange(ND)
    m = idx // npair
    f = np.asarray(geom.reverse)[m] * npair + (idx % nd) * nd + (idx // nd) % nd
    neg = _negate_index(geom.spatial_shape)
    w1 = np.empty((ND, ND), dtype=np.complex128)
    w2 = np.empty((ND, ND), dtype=np.complex128)
    for q in range(nvol):
        np.take(gamma_b[neg[q]], f, axis=0, out=w1, mode="clip")
        np.take(w1, f, axis=1, out=w2, mode="clip")        # w2[i, j] = Gamma_B(-q)[f(i), f(j)]
        oq = out[q]
        oq[...] = w2.T
        ph = geom.phase[q, m]                                # e^{+iq.dr} per slot
        np.multiply(oq, ph[None, :], out=oq)
        np.multiply(oq, ph.conj()[:, None], out=oq)
    return out


def _same_sector_map(a, b):
    """Whether two sector maps have the same content: ``B``, ``nd``, the
    packed index (it fixes the sizes, offsets and slot lists) and the
    exchange-image permutations (sources and permutations, which carry
    ``reverse``). The other tables are derived from these."""
    if a is b:
        return True
    return (a.B == b.B and a.nd == b.nd and np.array_equal(a.index, b.index)
            and all(x.source == y.source and np.array_equal(x.perm, y.perm)
                    for x, y in zip(a.gp_cross, b.gp_cross)))


def _gamma_p_packed(gamma_b, geom, *, out):
    """:func:`gamma_p` of packed arrays (see there)."""
    if not (isinstance(gamma_b, PackedArrays) and isinstance(out, PackedArrays)):
        raise ValueError("gamma_p: gamma_b and out must both be PackedArrays or both dense "
                         "arrays")
    _check_convention("gamma_p", geom)
    sm = gamma_b.sector_map
    if not _same_sector_map(sm, sz_sector_map(geom, sm.nd)):
        raise ValueError("gamma_p: the sector map of gamma_b (B = {}, nd = {}) is not the map "
                         "of the geometry (B = {})".format(sm.B, sm.nd, geom.B))
    if not _same_sector_map(out.sector_map, sm):
        raise ValueError("gamma_p: the sector map of out is not the sector map of gamma_b")
    for name, arr in (("gamma_b", gamma_b), ("out", out)):
        if arr.lead_shape != (geom.nvol,):
            raise ValueError("gamma_p: the packed {} has leading shape {}; the geometry needs "
                             "({},)".format(name, arr.lead_shape, geom.nvol))
    if any(np.shares_memory(o, g) for o in out.blocks for g in gamma_b.blocks):
        raise ValueError("gamma_p: out shares memory with gamma_b")
    cm = _consumer_map("gamma_p", sm)
    npair = cm.nd * cm.nd
    neg = _negate_index(geom.spatial_shape)
    for s, cross in enumerate(cm.gp_cross):
        src = gamma_b.blocks[cm.sigmas.index(cross.source)]
        dst = out.blocks[s]
        perm = cross.perm
        n = len(perm)
        m = cm.index[cm.offsets[s]:cm.offsets[s + 1]] // npair       # channel per target slot
        w1 = np.empty((n, n), dtype=np.complex128)
        w2 = np.empty((n, n), dtype=np.complex128)
        # mode="clip" lets np.take write into out directly (the default
        # "raise" buffers a copy when out is given); perm is a permutation
        # of the map, always in range, so nothing is clipped
        for q in range(geom.nvol):
            np.take(src[neg[q]], perm, axis=0, out=w1, mode="clip")
            np.take(w1, perm, axis=1, out=w2, mode="clip")  # w2[i, j] = src(-q)[perm i, perm j]
            oq = dst[q]
            oq[...] = w2.T
            ph = geom.phase[q, m]                            # e^{+iq.dr} per target slot
            np.multiply(oq, ph[None, :], out=oq)
            np.multiply(oq, ph.conj()[:, None], out=oq)
    return out


# ---------------------------------------------------------------------------
# The bond map of the effective interaction
# ---------------------------------------------------------------------------

def _describe(arr):
    if isinstance(arr, np.ndarray):
        return "{} {}".format(arr.dtype, arr.shape)
    return type(arr).__name__


def p_b(M, nd, *, out=None):
    """The bond map of the FLEX effective interaction:

    ``P_B[M](q)[(m; x, y), (m'; z, w)] = M(q)[(m; y, x), (m'; w, z)]``,

    the orbital swap of both pair indices within each channel (no channel
    reversal, no phase), applied to every leading index (``q``, bosonic
    frequency). ``W = P_B[1/2 Gamma_B chibar_B Gamma_B]`` is the
    second-order effective interaction read by
    ``flex_bond.calc_self_energy_bond``: its self-energy is ``beta
    dPhi_2/dG`` of ``Phi_2 = -(T/8) sum_{l,q} Tr[(chibar Gamma)**2]``.

    ``M`` is either a dense complex128 ``(..., ND, ND)`` array of code slots
    (``ND`` a positive multiple of ``nd**2``, ``nd = 2 norb_phys`` modes) or
    :class:`PackedArrays` on a sector map with ``nd`` modes. The swap sends
    the ph sector ``sigma`` onto ``-sigma``, so in the packed layout block
    ``sigma`` of the result is read from block ``-sigma`` of ``M`` through
    the swap permutation of the map (consumer ``"w_cross"``).

    ``out`` (optional) receives the result: an array of the layout, shape
    and dtype of ``M`` (packed: the same sector map and leading shape) that
    shares no memory with ``M``; ``None`` allocates it. Every argument is
    checked (``ValueError``) before anything is written. Returns ``out``.

    Workspace beyond ``out``: dense, none for a C-contiguous ``out`` (a
    strided copy) and one copy of ``M`` otherwise; packed, one
    leading slice of one sector block at a time, at most ``16 prod(lead[1:])
    max_sigma n_sigma**2`` bytes (``lead`` the leading shape, ``n_0 =
    ND / 2`` the largest block)."""
    if not (_is_int(nd) and int(nd) > 0 and int(nd) % 2 == 0):
        raise ValueError("p_b: nd must be a positive even integer (nd = 2*norb_phys), got "
                         "{!r}".format(nd))
    nd = int(nd)
    if isinstance(M, PackedArrays):
        return _p_b_packed(M, nd, out)
    if isinstance(out, PackedArrays):
        raise ValueError("p_b: M and out must both be PackedArrays or both dense arrays")
    npair = nd * nd
    if not (isinstance(M, np.ndarray) and M.dtype == np.complex128 and M.ndim >= 2
            and M.shape[-1] == M.shape[-2] and M.shape[-1] > 0 and M.shape[-1] % npair == 0):
        raise ValueError("p_b: M must be a complex128 (..., ND, ND) array with ND a positive "
                         "multiple of nd**2 = {}, got {}".format(npair, _describe(M)))
    if out is None:
        out = np.empty(M.shape, dtype=np.complex128)
    else:
        if not (isinstance(out, np.ndarray) and out.dtype == np.complex128
                and out.shape == M.shape):
            raise ValueError("p_b: out must be a complex128 array of shape {}, got {}".format(
                M.shape, _describe(out)))
        if np.shares_memory(out, M):
            raise ValueError("p_b: out shares memory with M")
    lead = M.shape[:-2]
    split = lead + (M.shape[-1] // npair, nd, nd) * 2
    src = M.reshape(split).swapaxes(-5, -4).swapaxes(-2, -1)    # [(m; y, x), (m'; w, z)]
    if out.flags.c_contiguous:
        np.copyto(out.reshape(split), src)                      # a view: no temporary
    else:
        out[...] = src.reshape(M.shape)
    return out


def _p_b_packed(M, nd, out):
    """:func:`p_b` of packed arrays (see there)."""
    sm = M.sector_map
    if sm.nd != nd:
        raise ValueError("p_b: the sector map of M has nd = {}, the call has nd = {}".format(
            sm.nd, nd))
    lead = M.lead_shape
    if out is None:
        out = PackedArrays(sm, tuple(np.empty(lead + (n, n), dtype=np.complex128)
                                     for n in sm.sizes))
    else:
        if not isinstance(out, PackedArrays):
            raise ValueError("p_b: M and out must both be PackedArrays or both dense arrays")
        if not _same_sector_map(out.sector_map, sm):
            raise ValueError("p_b: the sector map of out is not the sector map of M")
        if out.lead_shape != lead:
            raise ValueError("p_b: the packed out has leading shape {}, M has {}".format(
                out.lead_shape, lead))
        if any(np.shares_memory(o, b) for o in out.blocks for b in M.blocks):
            raise ValueError("p_b: out shares memory with M")
    cm = _consumer_map("w_cross", sm)
    for s, cross in enumerate(cm.pb_cross):
        src = M.blocks[cm.sigmas.index(cross.source)]
        dst = out.blocks[s]
        perm = cross.perm
        # one leading slice at a time through a workspace of one slice;
        # mode="clip" writes into out directly (the in-range permutation of
        # the map is never clipped), also into a non-contiguous dst[i] (a
        # view of the caller's block)
        slices = range(src.shape[0]) if src.ndim > 2 else (Ellipsis,)
        work = np.empty(src.shape[1:] if src.ndim > 2 else src.shape, dtype=np.complex128)
        for i in slices:
            np.take(src[i], perm, axis=-2, out=work, mode="clip")
            np.take(work, perm, axis=-1, out=dst[i], mode="clip")
        del work                    # one sector's workspace alive at a time
    return out


# ---------------------------------------------------------------------------
# The S_z-packed block store
# ---------------------------------------------------------------------------

_HOST_ONLY = "the S_z-packed store is host only"


def _module_name(xp):
    return getattr(xp, "__name__", type(xp).__name__)


class PackedBondBlockStore:
    """Owner of the S_z-packed ``(nmat, nvol, ND, ND)`` bond arrays.

    Every named slot holds the three diagonal sector blocks ``(nslices,
    nvol, n_sigma, n_sigma)`` of :class:`SectorMap` ``sector_map`` (host
    memory, full bosonic grid). The store implements an explicit subset of
    the ``flex_bond.BondBlockStore`` protocol -- the part read by
    ``flex_bond.assemble_bubble``, ``flex_bond.freq_batch_on``,
    ``flex_bond.calc_self_energy_bond``, ``eliashberg_bond.BondPairKernel``
    and this module:

    * attributes ``nmat``, ``nslices`` (``= nmat``), ``nvol``, ``ND``,
      ``nd``, ``B``, ``half_grid`` (``False``), ``resident`` (``"host"``),
      ``layout`` (``"sz_packed"``), ``sector_map``, ``slot_names``,
      ``nbytes``, ``released``. As in ``BondBlockStore``, ``nd`` is the
      PAIR size ``(2 norb_phys)**2`` of one channel block, while
      ``sector_map.nd`` is the mode count ``2 norb_phys``;
    * pair access: :meth:`put_pair` scatters a dense pair block into the
      sector blocks and refuses one with a nonzero entry between two
      sectors; :meth:`get_pair` and :meth:`get_pair_full` return a NEW
      dense pair block (off-sector zeros; writes do not reach the store);
    * frequency-batch access in :class:`PackedArrays`
      (:meth:`put_freq_batch`, :meth:`get_freq_batch`,
      :meth:`get_freq_batch_on`) and the whole slot (:meth:`packed`); the
      returned objects hold VIEWS of the slot (writes reach the store; they
      must not be used once the slot is released -- not enforced, as in
      ``BondBlockStore``);
    * ``grid_index``, ``detach`` (the slot as an owned
      :class:`PackedArrays`; the store then treats the name as released),
      ``release_slot``, ``release``, ``to_host`` (a no-op) and the context
      manager (``release`` on exit).

    Deliberately absent, so that a half-grid or device consumer fails with
    ``AttributeError`` instead of mis-reading packed slots:
    ``put_pair_fullgrid``, ``expand_full``, ``stager`` and ``copy_workers``.
    ``flex_bond.dress_and_build_w``, ``flex_bond.HalfGridDressing`` and
    ``flex_bond.BondDeviceContext`` are not supported on this store.

    ``aliases`` maps extra names onto allocated slots (as in
    ``BondBlockStore``); an alias cannot be detached or released on its
    own. ``half_grid=True`` and an array module other than numpy are
    refused.
    """

    def __init__(self, nmat, nvol, sector_map, names, aliases=None, *, half_grid=False, xp=np):
        if half_grid:
            raise ValueError("PackedBondBlockStore: the half bosonic grid is not available on "
                             "the S_z-packed store")
        if xp is not np:
            raise ValueError("PackedBondBlockStore: {} (xp must be numpy), got {}".format(
                _HOST_ONLY, _module_name(xp)))
        if not isinstance(sector_map, SectorMap):
            raise ValueError("PackedBondBlockStore: sector_map must be a SectorMap, got {}"
                             .format(type(sector_map).__name__))
        names = tuple(names)
        if len(set(names)) != len(names):
            raise ValueError("PackedBondBlockStore: duplicate array names {}".format(names))
        aliases = dict(aliases or {})
        for alias, target in aliases.items():
            if alias in names or target not in names:
                raise ValueError("PackedBondBlockStore: alias {!r} -> {!r} must name a new name "
                                 "for an allocated array {}".format(alias, target, names))
        self._aliases = aliases
        self.sector_map = sector_map
        self.nmat = self.nslices = int(nmat)
        self.nvol = int(nvol)
        self.ND = int(sector_map.ND)
        self.nd = int(sector_map.nd) ** 2
        self.B = int(sector_map.B)
        self.half_grid = False
        self.resident = "host"
        self.layout = "sz_packed"
        self._arrays = {}
        self._released_slots = set()
        try:
            for name in names:
                self._arrays[name] = tuple(
                    np.zeros((self.nslices, self.nvol, n, n), dtype=np.complex128)
                    for n in sector_map.sizes)
        except MemoryError:
            self._arrays.clear()
            raise
        self._released = False

    # -- bookkeeping -------------------------------------------------------
    @property
    def slot_names(self):
        """The allocated slots (released or detached ones excluded) in
        allocation order, then the aliases whose target is still allocated."""
        return tuple(self._arrays) + tuple(a for a, t in self._aliases.items()
                                           if t in self._arrays)

    @property
    def nbytes(self):
        """Bytes of the slots currently allocated (aliases counted once)."""
        return int(sum(b.nbytes for blocks in self._arrays.values() for b in blocks))

    @property
    def released(self):
        return self._released

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def _arr(self, name):
        if self._released:
            raise RuntimeError("PackedBondBlockStore: released")
        name = self._aliases.get(name, name)
        if name in self._released_slots:
            raise KeyError("PackedBondBlockStore: slot {!r} was released".format(name))
        if name not in self._arrays:
            raise KeyError("PackedBondBlockStore: no slot {!r}".format(name))
        return self._arrays[name]

    def _refuse_alias(self, name, what):
        if name in self._aliases:
            raise ValueError("PackedBondBlockStore: cannot {} {!r}, an alias of {!r}".format(
                what, name, self._aliases[name]))

    def _range(self, l0, l1):
        if not (_is_int(l0) and _is_int(l1) and 0 <= int(l0) < int(l1) <= self.nslices):
            raise ValueError("PackedBondBlockStore: the slice range [{}, {}) must lie in "
                             "[0, {}) and be non-empty".format(l0, l1, self.nslices))
        return int(l0), int(l1)

    def _host_only(self, method, xp):
        if not (xp is None or xp is np):
            raise ValueError("PackedBondBlockStore.{}: {} (xp must be None or numpy), got {}"
                             .format(method, _HOST_ONLY, _module_name(xp)))

    # -- pair access -------------------------------------------------------
    def _pair_maps(self, alpha, beta):
        """Per sector: the pair positions ``rp`` (ascending) of that sector
        and the local rows / columns of channel ``alpha`` / ``beta``."""
        B = self.B
        if not (_is_int(alpha) and _is_int(beta) and 0 <= int(alpha) < B
                and 0 <= int(beta) < B):
            raise ValueError("PackedBondBlockStore: channel indices must be integers in [0, {}), "
                             "got ({}, {})".format(B, alpha, beta))
        sm = _consumer_map("store_scatter", self.sector_map)
        P = self.nd
        psig = sm.sector_of[:P]                       # sigma of each pair (channel independent)
        out = []
        for sigma in sm.sigmas:
            rp = np.flatnonzero(psig == sigma)
            out.append((rp, sm.local[int(alpha) * P + rp], sm.local[int(beta) * P + rp]))
        return psig, out

    def put_pair(self, name, alpha, beta, block):
        """Scatter the dense pair block ``(alpha, beta)`` ``(nslices, nvol,
        nd, nd)`` into the sector blocks of slot ``name``. A block with a
        nonzero entry between two sectors is refused before anything is
        written."""
        arrs = self._arr(name)
        psig, maps = self._pair_maps(alpha, beta)
        P = self.nd
        block = np.asarray(block)
        want = (self.nslices, self.nvol, P, P)
        if block.shape != want:
            raise ValueError("PackedBondBlockStore.put_pair: the block of '{}' must have shape "
                             "{}, got {}".format(name, want, block.shape))
        off_mask = psig[:, None] != psig[None, :]
        off_max = 0.0
        for l in range(self.nslices):                 # one slice at a time: small temporaries
            off = block[l][..., off_mask]
            if off.size and np.any(off != 0):
                v = float(np.max(np.abs(off)))
                if not np.isfinite(v):                # NaN or inf: refused as such
                    off_max = v
                    break
                off_max = max(off_max, v)
        if off_max != 0.0:
            raise ValueError(
                "PackedBondBlockStore: the pair block (%d, %d) of '%s' has S_z-off-sector "
                "entries (max |entry| = %.3e): the Green function it was built from is not "
                "block-diagonal in spin; build it per spin sector (canonical_band_eigen) or "
                "project a borrowed one (project_borrowed_green)"
                % (int(alpha), int(beta), name, off_max))
        for arr, (rp, rows, cols) in zip(arrs, maps):
            arr[:, :, rows[:, None], cols] = block[:, :, rp[:, None], rp]

    def get_pair(self, name, alpha, beta):
        """The dense pair block ``(alpha, beta)`` ``(nslices, nvol, nd, nd)``
        of slot ``name``: a NEW array with exact zeros between sectors."""
        arrs = self._arr(name)
        _, maps = self._pair_maps(alpha, beta)
        P = self.nd
        out = np.zeros((self.nslices, self.nvol, P, P), dtype=np.complex128)
        for arr, (rp, rows, cols) in zip(arrs, maps):
            out[:, :, rp[:, None], rp] = arr[:, :, rows[:, None], cols]
        return out

    def get_pair_full(self, name, alpha, beta, xp=None):
        """:meth:`get_pair` (the full grid is stored); ``xp`` must be None or
        numpy."""
        self._host_only("get_pair_full", xp)
        return self.get_pair(name, alpha, beta)

    # -- frequency batches -------------------------------------------------
    def packed(self, name):
        """The whole slot ``name`` as :class:`PackedArrays` of VIEWS,
        leading shape ``(nslices, nvol)``."""
        return PackedArrays(self.sector_map, tuple(self._arr(name)))

    def get_freq_batch(self, name, l0, l1):
        """The slices ``l0 ... l1 - 1`` of slot ``name`` as
        :class:`PackedArrays` of VIEWS."""
        l0, l1 = self._range(l0, l1)
        return PackedArrays(self.sector_map, tuple(b[l0:l1] for b in self._arr(name)))

    def get_freq_batch_on(self, name, l0, l1, xp):
        """:meth:`get_freq_batch`; ``xp`` must be None or numpy."""
        self._host_only("get_freq_batch_on", xp)
        return self.get_freq_batch(name, l0, l1)

    def put_freq_batch(self, name, l0, l1, batch):
        """Copy :class:`PackedArrays` ``batch`` (the store's sector map,
        leading shape ``(l1 - l0, nvol)``) into the slices ``l0 ... l1 - 1``
        of slot ``name``."""
        arrs = self._arr(name)
        l0, l1 = self._range(l0, l1)
        if not (isinstance(batch, PackedArrays)
                and _same_sector_map(batch.sector_map, self.sector_map)
                and batch.lead_shape == (l1 - l0, self.nvol)):
            raise ValueError("PackedBondBlockStore.put_freq_batch: the batch of '{}' must be "
                             "PackedArrays of the store's sector map with leading shape {}"
                             .format(name, (l1 - l0, self.nvol)))
        for dst, src in zip(arrs, batch.blocks):
            dst[l0:l1] = src

    def grid_index(self, s0, s1):
        """The grid index ``l`` of the stored slices ``s0 ... s1 - 1`` (the
        full grid: the slice index itself)."""
        return np.arange(int(s0), int(s1), dtype=np.int64)

    def to_host(self):
        """The identity (the store is host memory)."""
        return None

    # -- ownership ---------------------------------------------------------
    def detach(self, name):
        """Remove slot ``name`` and return it as an owned
        :class:`PackedArrays`."""
        self._refuse_alias(name, "detach")
        arrs = self._arr(name)
        del self._arrays[name]
        self._released_slots.add(name)              # later access: "was released"
        return PackedArrays(self.sector_map, arrs)

    def release_slot(self, name):
        """Free one slot; later access raises ``KeyError``; releasing it
        again is tolerated."""
        if self._released:
            raise RuntimeError("PackedBondBlockStore: released")
        self._refuse_alias(name, "release")
        if name not in self._arrays:
            if name in self._released_slots:
                return
            raise KeyError("PackedBondBlockStore: no slot {!r}".format(name))
        del self._arrays[name]
        self._released_slots.add(name)

    def release(self):
        self._arrays.clear()
        self._released = True


def open_store(plan, sector_map=None):
    """The block store of an admission plan: ``flex_bond.BondBlockStore(
    plan.nslices, plan.nvol, plan.ND, plan.nd**2, plan.slots,
    dict(plan.aliases))`` for a dense plan (no sector map), a
    :class:`PackedBondBlockStore` on ``sector_map`` (required; its ``B`` and
    ``nd`` must be the plan's) for an ``sz_packed`` plan. The plan carries no
    channel ``reverse``; the bubble checks the store's map against the bond
    channels it is called with."""
    if not isinstance(plan, SOBondPlan):
        raise ValueError("open_store: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    if plan.layout == "dense":
        if sector_map is not None:
            raise ValueError("open_store: a dense plan takes no sector map")
        return flex_bond.BondBlockStore(plan.nslices, plan.nvol, plan.ND, plan.nd ** 2,
                                        plan.slots, dict(plan.aliases))
    if sector_map is None:
        raise ValueError("open_store: an sz_packed plan needs the sector map (sz_sector_map)")
    if not isinstance(sector_map, SectorMap):
        raise ValueError("open_store: sector_map must be a SectorMap, got {}".format(
            type(sector_map).__name__))
    if (sector_map.B, sector_map.nd) != (plan.B, plan.nd):
        raise ValueError("open_store: the sector map (B = {}, nd = {}) is not the plan's (B = {}, "
                         "nd = {})".format(sector_map.B, sector_map.nd, plan.B, plan.nd))
    return PackedBondBlockStore(plan.nslices, plan.nvol, sector_map, plan.slots,
                                dict(plan.aliases))


# ---------------------------------------------------------------------------
# Bubble and dressing
# ---------------------------------------------------------------------------

def so_bond_bubble(green_kw, green0_tail, beta, view, spatial_shape, *, plan, freq, sink):
    """The spinful bond bubble ``chibar_B`` in the code slots.

    ``chibar_B(q, i nu)[(m;x,y),(m';z,w)] = e^{-iq.dr_m}
    chibar(-q, +i nu)[(+dr_m; y,x), (-dr_m'; z,w)]`` in terms of the label
    bubble ``chibar_{(d;X,Y),(d';Z,W)}(q, i nu) = -(T/N) sum_{k,n}
    e^{+ik.d} e^{+i(k+q).d'} G_YZ(k, i w_n) G_WX(k+q, i w_n - i nu)``,
    computed by the production bond bubble (``bubble.bond_bubble_static``,
    ``flex_bond.assemble_bubble``) with the general orbital product on the
    spin-orbital modes (``nd = 2 * norb_phys``, spin-block index).

    Parameters
    ----------
    green_kw : ndarray (1, nmat, nvol, nd, nd)
        The Green function with its single block axis (tail-subtracted when
        ``green0_tail`` is given).
    green0_tail : ndarray or None
        Its paired tail (``None`` disables the tail correction).
    beta : float
    view : object with ``delta_r`` and ``reverse``
        The bond channels (a ``BondTopology``, ``BondSetView`` or plain
        namespace).
    spatial_shape : (Nx, Ny, Nz)
    plan : SOBondPlan
    freq : ``"static"`` or ``"all"``
        ``"static"``: the ``nu = 0`` bubble is written into ``sink``, a
        caller-owned complex128 ``(nvol, ND, ND)`` array, which is returned.
        ``"all"``: every bosonic slice is written into the store slot
        ``"chibar"`` of ``sink`` (a ``flex_bond.BondBlockStore`` of the
        full grid, pair size ``nd**2``); returns ``None``.
    sink (S_z-packed layout)
        ``"static"``: a :class:`PackedArrays` of leading shape ``(nvol,)``;
        the bubble is computed dense and packed only if its entries between
        two S_z sectors are exactly zero (otherwise refused, with the leak
        ``max|off-sector| / max|bubble|``, before the sink is written).
        ``"all"``: a :class:`PackedBondBlockStore` (:func:`open_store`); its
        ``put_pair`` scatters each dense pair block and refuses one with an
        off-sector entry. The observed ``layout`` is the sink's
        (``PackedArrays`` / the store's ``layout`` attribute, ``"dense"``
        otherwise) and must be the plan's; a packed sink's sector map must
        be the map of the bond channels (``sz_sector_map(view, nd)`` by
        content, so a map of another ``reverse`` is refused). A Green function that is not
        block-diagonal in spin fails both paths.

    The mesh is validated first (:func:`slot_geometry`: channels that
    coincide modulo the mesh are refused), then the boundary is checked
    against ``plan`` (:meth:`SOBondPlan.require`); both happen before any
    bubble work. This function never allocates the sink. A non-finite
    bubble raises ``hartree_fock.NonFiniteError`` (the static bubble is
    checked before it is written into ``sink``, which is then untouched; the
    every-slice path refuses a non-finite block before storing it). If it
    raises after the bubble has started, the content of a store is
    undefined.
    """
    return _bubble_on_geom(green_kw, green0_tail, beta, slot_geometry(view, spatial_shape),
                           plan=plan, freq=freq, sink=sink)


def _off_sector_max(dense, sector_map):
    """``max |entry|`` of a dense ``(nvol, ND, ND)`` code-slot array outside
    its diagonal sector blocks (one momentum at a time)."""
    so = sector_map.sector_of
    mask = so[:, None] != so[None, :]
    return max((float(np.max(np.abs(dense[q][mask]))) for q in range(dense.shape[0])),
               default=0.0)


def _bubble_on_geom(green_kw, green0_tail, beta, geom, *, plan, freq, sink):
    """:func:`so_bond_bubble` on a built geometry (same checks and result)."""
    _check_convention("so_bond_bubble", geom)
    if not isinstance(plan, SOBondPlan):
        raise ValueError("so_bond_bubble: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    if freq not in ("static", "all"):
        raise ValueError("so_bond_bubble: freq must be 'static' or 'all', got {!r}".format(freq))
    shape = tuple(np.shape(green_kw))
    if len(shape) != 5 or shape[0] != 1 or shape[2] != geom.nvol or shape[3] != shape[4]:
        raise ValueError("so_bond_bubble: green_kw must have shape (1, nmat, {}, nd, nd), got {}"
                         .format(geom.nvol, shape))
    nmat, nvol, nd = int(shape[1]), int(shape[2]), int(shape[3])
    ND = geom.B * nd * nd
    if freq == "static":
        layout = "sz_packed" if isinstance(sink, PackedArrays) else "dense"
    else:
        layout = getattr(sink, "layout", "dense")
    observed = dict(nvol=nvol, B=geom.B, nd=nd, ND=ND, layout=layout, nmat=nmat)
    if freq == "all":
        observed.update(placement=getattr(sink, "resident", None),
                        slots=tuple(getattr(sink, "slot_names", ())))
    plan.require("so_bond_bubble", freq=freq, **observed)
    if layout == "sz_packed":
        sm = sink.sector_map
        what = "static sink" if freq == "static" else "store"
        if freq == "static" and sink.lead_shape != (nvol,):
            raise ValueError("so_bond_bubble: the packed static sink must have leading shape "
                             "({},), got {}".format(nvol, sink.lead_shape))
        if not (isinstance(sm, SectorMap) and sm.nd == nd
                and _same_sector_map(sm, sz_sector_map(geom, nd))):
            raise ValueError("so_bond_bubble: the sector map of the packed {} is not the map of "
                             "the bond channels (B = {}, nd = {}, reverse = {})".format(
                                 what, geom.B, nd, geom.reverse.tolist()))
    bond_set = bond_channels.BondSetView(
        bond_channels.BondTopology(geom.delta_r, geom.reverse, {}))
    if freq == "static" and layout == "sz_packed":
        chibar = bubble.bond_bubble_static(green_kw, green0_tail, beta, bond_set,
                                           spatial_shape=geom.spatial_shape)
        if not np.isfinite(chibar).all():
            raise NonFiniteError("so_bond_bubble: non-finite static bond bubble")
        off_max = _off_sector_max(chibar, sm)
        if off_max != 0.0:
            raise ValueError(
                "so_bond_bubble: the static bond bubble has S_z-off-sector entries (leak = "
                "%.3e): the Green function is not block-diagonal in spin; build it per spin "
                "sector (canonical_band_eigen) or project a borrowed one "
                "(project_borrowed_green)" % (off_max / float(np.max(np.abs(chibar)))))
        # gathered straight into the sink, one momentum at a time (no
        # temporary larger than one sector block of one momentum)
        for s, dst in enumerate(sink.blocks):
            idx = sm.index[sm.offsets[s]:sm.offsets[s + 1]]
            for q in range(nvol):
                dst[q] = chibar[q][idx[:, None], idx]
        return sink
    if freq == "static":
        if not isinstance(sink, np.ndarray) or sink.dtype != np.complex128 \
                or sink.shape != (nvol, ND, ND):
            raise ValueError("so_bond_bubble: the static sink must be a complex128 array of "
                             "shape {}, got {} {}".format(
                                 (nvol, ND, ND), getattr(sink, "dtype", type(sink).__name__),
                                 tuple(getattr(sink, "shape", ()))))
        chibar = bubble.bond_bubble_static(green_kw, green0_tail, beta, bond_set,
                                           spatial_shape=geom.spatial_shape)
        if not np.isfinite(chibar).all():
            raise NonFiniteError("so_bond_bubble: non-finite static bond bubble")
        sink[...] = chibar
        del chibar
        return sink
    want = (nmat, nvol, ND, nd * nd)
    got = tuple(getattr(sink, name, None) for name in ("nmat", "nvol", "ND", "nd"))
    if got != want or getattr(sink, "half_grid", None) is not False:
        raise ValueError("so_bond_bubble: the store must hold the full grid with (nmat, nvol, "
                         "ND, pair size) = {}, got {} (half_grid = {})".format(
                             want, got, getattr(sink, "half_grid", None)))
    flex_bond.assemble_bubble(sink, green_kw, green0_tail, beta, bond_set, geom.spatial_shape,
                              None)
    return None


def so_bond_dress(X, gamma_b, *, plan, l0, spatial_shape, cond_tol=bond_channels._BOND_COND_FLOOR,
                  diagnostics=None):
    """The RPA dressing of a batch of spinful code bubbles:

    ``chi_B(q, i nu) = [1 + chibar_B(q, i nu) Gamma_B(q)]^{-1} chibar_B(q, i nu)``

    per slice and momentum (``bond_channels.dress_batch`` with the channel
    ``"spinful"``: the full antisymmetrized vertex, one plus sign). In the
    code slots this is the bubble-rule image of the label ring
    ``[1 + chibar(q) Gamma(-q)]^{-1} chibar(q)``; no further momentum flip
    enters.

    ``X`` is a batch ``(nb, nvol, ND, ND)`` and so is the result; a static
    caller passes ``X[None]`` with ``l0 = nmat // 2`` and takes ``[0]``.
    ``gamma_b`` is ``(nvol, ND, ND)``. ``l0`` (the grid index of the
    first slice) and ``plan.nmat`` only name the frequency in a refusal. The
    batch length is checked against ``plan.dress_capacity``
    (:meth:`SOBondPlan.require`) before any work, then ``X`` and
    ``gamma_b`` must be finite (``hartree_fock.NonFiniteError`` naming the
    argument otherwise). ``cond_tol`` is the
    conditioning guard of ``dress_batch`` (``None`` disables it).

    S_z-packed layout (a plan of layout ``"sz_packed"``): ``X`` is a
    :class:`PackedArrays` of leading shape ``(nb, nvol)`` and ``gamma_b`` a
    :class:`PackedArrays` of leading shape ``(nvol,)`` with the same sector
    map; the result is a new :class:`PackedArrays` ``(nb, nvol)`` on that
    map. Each ph sector is dressed on its own, ``chi_s = [1 + X_s
    Gamma_s]^{-1} X_s``. The conditioning guard is the dense one: per
    ``(l, q)`` the singular values of the three sector blocks together are
    those of the dense block, so ``sigma_min`` (``sigma_max``) is the
    minimum (maximum) over the sectors and the two criteria of
    ``bond_channels.bond_conditioning_score`` are formed from them. A
    refusal is the ``bond_channels.BondConditioningError`` of
    ``dress_batch`` for the same block, its message followed by ``" (S_z
    sector of sigma_min: <s>; of sigma_max: <s'>)"`` (the second clause
    omitted when the two sectors are the same). Like the dense call, the
    packed guard always decomposes every block exactly on the host: there
    is no ``guard_method``, ``guard_freqs`` or warn policy.

    ``diagnostics`` (a dict, optional) receives ``cond_min``, the guard
    score (``None`` without a guard); packed, also ``sigma_min_sector`` and
    ``sigma_max_sector``, the sectors (``"0"``, ``"+1"``, ``"-1"``) of the
    two singular values at the worst ``(l, q)`` (``None`` without a guard).
    """
    return _dress_impl(X, gamma_b, plan=plan, l0=l0, spatial_shape=spatial_shape,
                       cond_tol=cond_tol, diagnostics=diagnostics)


def _dress_impl(X, gamma_b, *, plan, l0, spatial_shape, cond_tol, diagnostics):
    """:func:`so_bond_dress` (same checks and result, ``diagnostics`` as
    there). The dressing needs no geometry: only ``ND`` (dense) or the
    sector map of ``X`` (packed) and the mesh shape."""
    if not isinstance(plan, SOBondPlan):
        raise ValueError("so_bond_dress: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    if isinstance(X, PackedArrays):
        layout = "sz_packed"
        shape = X.lead_shape + (X.sector_map.ND, X.sector_map.ND)
    else:
        layout = "dense"
        shape = tuple(np.shape(X))
    if len(shape) != 4 or shape[2] != shape[3]:
        raise ValueError("so_bond_dress: X must be a batch (nb, nvol, ND, ND), got shape {}"
                         .format(shape))
    plan.require("so_bond_dress", nvol=shape[1], ND=shape[2], batch=shape[0], layout=layout)
    if isinstance(X, PackedArrays) != isinstance(gamma_b, PackedArrays):
        raise ValueError("so_bond_dress: X and gamma_b must both be PackedArrays or both dense "
                         "arrays")
    if layout == "sz_packed":
        return _dress_packed(X, gamma_b, nmat=plan.nmat, l0=l0, spatial_shape=spatial_shape,
                             cond_tol=cond_tol, diagnostics=diagnostics)
    for arg, arr in (("X", X), ("gamma_b", gamma_b)):
        if not np.isfinite(arr).all():
            raise NonFiniteError("so_bond_dress: {} has non-finite entries".format(arg))
    chi, cond_min = bond_channels.dress_batch(X, gamma_b, "spinful", l0=l0, nmat=plan.nmat,
                                              spatial_shape=spatial_shape, cond_tol=cond_tol)
    if diagnostics is not None:
        diagnostics["cond_min"] = cond_min
    return chi


def _dress_packed(X, gamma_b, *, nmat, l0, spatial_shape, cond_tol, diagnostics):
    """The packed branch of :func:`_dress_impl` (see :func:`so_bond_dress`)."""
    sm = X.sector_map
    nb, nvol = X.lead_shape
    if not _same_sector_map(gamma_b.sector_map, sm):
        raise ValueError("so_bond_dress: the sector map of gamma_b is not the sector map of X")
    if gamma_b.lead_shape != (nvol,):
        raise ValueError("so_bond_dress: the packed gamma_b has leading shape {}; X needs ({},)"
                         .format(gamma_b.lead_shape, nvol))
    Nx, Ny, Nz = (int(x) for x in spatial_shape)
    if Nx * Ny * Nz != nvol:
        raise ValueError("so_bond_dress: prod(spatial_shape) != nvol")
    for arg, arr in (("X", X), ("gamma_b", gamma_b)):
        if not all(np.isfinite(b).all() for b in arr.blocks):
            raise NonFiniteError("so_bond_dress: {} has non-finite entries".format(arg))
    # the sectors to dress, by storage position; the result keeps the map of X
    cm = _consumer_map("dress_sectors", sm)
    sectors = [SECTOR_ORDER.index(name) for name in cm.sector_order]
    mats = {}
    for s in sectors:
        mat = np.matmul(X.blocks[s], gamma_b.blocks[s][None])
        idx = np.arange(mat.shape[-1])
        mat[..., idx, idx] += 1.0
        mats[s] = mat
    cond_min = s_min = s_max = None
    if cond_tol is not None:
        # per (l, q): the extreme singular values over the sectors are those
        # of the dense block (its singular values are the union)
        smin = np.full((len(SECTOR_ORDER), nb * nvol), np.inf)
        smax = np.full((len(SECTOR_ORDER), nb * nvol), -np.inf)
        for s, mat in mats.items():
            n = mat.shape[-1]
            sv = np.linalg.svd(mat.reshape(nb * nvol, n, n), compute_uv=False)
            smin[s], smax[s] = sv[:, -1], sv[:, 0]
        gmin, gmax = smin.min(axis=0), smax.max(axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = gmin / gmax
        ratio = np.where(np.isfinite(ratio), ratio, 0.0)
        pole = gmin / np.maximum(1.0, gmax)
        pole = np.where(np.isfinite(pole), pole, 0.0)
        score = np.minimum(ratio, pole)
        i = int(np.argmin(score))
        cond_min = float(score[i])
        s_min = SECTOR_ORDER[int(np.argmin(smin[:, i]))]
        s_max = SECTOR_ORDER[int(np.argmax(smax[:, i]))]
        if cond_min <= cond_tol:
            l, q = l0 + i // nvol, i % nvol
            qx, rem = divmod(q, Ny * Nz)
            qy, qz = divmod(rem, Nz)
            worst = (cond_min, i, float(ratio[i]), float(pole[i]), float(gmin[i]),
                     float(gmax[i]))
            cause = bond_channels._conditioning_refusal("spinful", worst, cond_tol, 1, 1)
            where = ("sigma_min: {}".format(s_min) if s_min == s_max else
                     "sigma_min: {}; of sigma_max: {}".format(s_min, s_max))
            raise bond_channels.BondConditioningError(
                bond_channels._dress_refusal_text("spinful", l, nmat, (qx, qy, qz))
                + " (S_z sector of {})".format(where),
                channel="spinful", iq=q, q=(qx, qy, qz), l=l, worst=cond_min,
                ratio=worst[2], pole=worst[3], smin=worst[4], smax=worst[5],
                cond_tol=cond_tol) from cause
    out = PackedArrays.zeros(sm, (nb, nvol))
    for s, mat in mats.items():
        out.blocks[s][...] = np.linalg.solve(mat, X.blocks[s])
    if diagnostics is not None:
        diagnostics["cond_min"] = cond_min
        diagnostics["sigma_min_sector"] = s_min
        diagnostics["sigma_max_sector"] = s_max
    return out


# ---------------------------------------------------------------------------
# Static pairing vertex
# ---------------------------------------------------------------------------

#: relative bound of ``|Gamma^P - Gamma|`` (Frobenius norms over all q)
_EXCHANGE_RTOL = 1e-13


def so_bond_pairing_vertex_parts(gamma_b, chi_b, geom, *, plan):
    """The static pairing vertex of the spinful bond vertex, in two parts:

    ``V_bare = 0.25 (Gamma_B + Gamma_B^P)``,
    ``V_fluct(q) = -Gamma_B(q) chi_B(q) Gamma_B(q)``,

    each ``(nvol, ND, ND)`` complex128 in the code slots of
    :data:`SLOT_CONVENTION`, which are the block layout of
    :class:`hwave.solver.eliashberg_bond.BondPairKernel` with ``norb = nd``
    spin-orbital modes (no instantaneous vertex: the bare part lives on all
    channel blocks). ``Gamma_B^P`` is the exchange image (:func:`gamma_p`).

    ``chi_b`` is the DRESSED static susceptibility ``chi_B = [1 + chibar_B
    Gamma_B]^{-1} chibar_B`` (the output of :func:`so_bond_dress`), never the
    bare bubble ``chibar_B``. ``V_fluct`` is the direct particle-hole image
    of the screened vertex with weight one; the bare part already holds both
    exchange images, hence its 1/2 (``0.25 (Gamma + Gamma^P) = 1/2 Gamma``).

    The boundary is checked against ``plan`` (:meth:`SOBondPlan.require`:
    ``nvol``, ``ND`` and the geometry's ``B``) before anything is allocated;
    a ``chi_b`` that shares memory with ``gamma_b`` is refused. The helper
    allocates the two outputs and one per-q workspace of two ``ND x ND``
    matrices, and never writes to ``gamma_b`` or ``chi_b``. A non-finite
    entry in ``gamma_b`` or ``chi_b`` is refused (``ValueError`` naming the
    argument) before the outputs are allocated. ``Gamma_B`` must equal its
    exchange image: ``|Gamma^P - Gamma| > 1e-13 |Gamma|`` (Frobenius norms
    over all q), or a non-finite residual, is a ``ValueError``; the
    symmetrization is checked, not silently applied.

    S_z-packed layout: ``gamma_b`` and ``chi_b`` may both be
    :class:`PackedArrays` of leading shape ``(nvol,)`` on the sector map of
    ``geom`` (a plan of layout ``"sz_packed"``). Both parts are then
    :class:`PackedArrays` on that map, computed sector by sector,
    ``V_bare@s = 0.25 (Gamma@s + Gamma^P@s)`` and ``V_fluct@s(q) =
    -Gamma@s(q) chi@s(q) Gamma@s(q)``, with a per-q workspace of two
    ``n_s x n_s`` matrices; the exchange-symmetry check runs over all
    sectors with the same rule (the Frobenius norms summed over the
    sectors and q).

    Returns ``(V_bare, V_fluct)``.
    """
    name = "so_bond_pairing_vertex_parts"
    if not isinstance(plan, SOBondPlan):
        raise ValueError("{}: plan must be a SOBondPlan (plan_so_bond), got {}".format(
            name, type(plan).__name__))
    if isinstance(gamma_b, PackedArrays) or isinstance(chi_b, PackedArrays):
        return _parts_packed(gamma_b, chi_b, geom, plan=plan)
    if not isinstance(gamma_b, np.ndarray) or gamma_b.dtype != np.complex128 \
            or gamma_b.ndim != 3:
        raise ValueError("{}: gamma_b must be a complex128 (nvol, ND, ND) array".format(name))
    nvol, ND = int(gamma_b.shape[0]), int(gamma_b.shape[1])
    _check_convention(name, geom)
    plan.require(name, nvol=nvol, ND=ND, B=getattr(geom, "B", None))
    _nd_of(gamma_b.shape, geom, name, "gamma_b")
    if not isinstance(chi_b, np.ndarray) or chi_b.dtype != np.complex128 \
            or chi_b.shape != gamma_b.shape:
        raise ValueError("{}: chi_b must be a complex128 array of shape {}, got {} {}".format(
            name, gamma_b.shape, getattr(chi_b, "dtype", type(chi_b).__name__),
            tuple(getattr(chi_b, "shape", ()))))
    if np.shares_memory(chi_b, gamma_b):
        raise ValueError("{}: chi_b shares memory with gamma_b".format(name))
    for arg, arr in (("gamma_b", gamma_b), ("chi_b", chi_b)):
        if not np.isfinite(arr).all():
            raise ValueError("{}: {} has non-finite entries".format(name, arg))

    v_bare = np.empty((nvol, ND, ND), dtype=np.complex128)
    gamma_p(gamma_b, geom, out=v_bare)                       # Gamma^P, in place
    ws1 = np.empty((ND, ND), dtype=np.complex128)
    ws2 = np.empty((ND, ND), dtype=np.complex128)
    asym = 0.0
    norm = 0.0
    for q in range(nvol):
        np.subtract(v_bare[q], gamma_b[q], out=ws1)
        asym += float(np.vdot(ws1, ws1).real)
        norm += float(np.vdot(gamma_b[q], gamma_b[q]).real)
    asym, norm = np.sqrt(asym), np.sqrt(norm)
    # written so that a NaN residual (a non-finite Gamma_B) is refused too
    if not asym <= _EXCHANGE_RTOL * norm:
        raise ValueError("so_bond_pairing_vertex_parts: Gamma_B violates the particle-exchange "
                         "symmetry (|Gamma^P - Gamma| = {:.3e} > {:g} |Gamma|)".format(
                             asym, _EXCHANGE_RTOL))
    v_fluct = np.empty((nvol, ND, ND), dtype=np.complex128)
    for q in range(nvol):
        vb = v_bare[q]
        np.add(vb, gamma_b[q], out=vb)
        vb *= 0.25
        np.matmul(gamma_b[q], chi_b[q], out=ws2)
        vf = v_fluct[q]
        np.matmul(ws2, gamma_b[q], out=vf)
        vf *= -1.0
    return v_bare, v_fluct


def _exchange_residual(pairs):
    """``(|Gamma^P - Gamma|, |Gamma|)``, Frobenius norms over every
    ``(image, vertex)`` pair of ``(.., n, n)`` arrays in ``pairs``, formed
    one ``n x n`` matrix at a time."""
    asym = norm = 0.0
    for img, g in pairs:
        ws = np.empty(g.shape[-2:], dtype=np.complex128)
        for q in range(g.shape[0]):
            np.subtract(img[q], g[q], out=ws)
            asym += float(np.vdot(ws, ws).real)
            norm += float(np.vdot(g[q], g[q]).real)
    return np.sqrt(asym), np.sqrt(norm)


def _refuse_exchange(name, asym, norm):
    # written so that a NaN residual (a non-finite Gamma_B) is refused too
    if not asym <= _EXCHANGE_RTOL * norm:
        raise ValueError("{}: Gamma_B violates the particle-exchange symmetry (|Gamma^P - "
                         "Gamma| = {:.3e} > {:g} |Gamma|)".format(name, asym, _EXCHANGE_RTOL))


def _packed_nvol(name, gamma_b):
    """``nvol`` of a packed ``gamma_b``, whose leading shape must be ``(nvol,)``."""
    lead = gamma_b.lead_shape
    if len(lead) != 1:
        raise ValueError("{}: the packed gamma_b must have leading shape (nvol,), got {}"
                         .format(name, lead))
    return int(lead[0])


def _check_packed_map(name, gamma_b, geom):
    sm = gamma_b.sector_map
    if not _same_sector_map(sm, sz_sector_map(geom, sm.nd)):
        raise ValueError("{}: the sector map of gamma_b (B = {}, nd = {}) is not the map of the "
                         "geometry (B = {})".format(name, sm.B, sm.nd, geom.B))


def _blocks(a):
    """The sector blocks of a packed array; a dense array as one block."""
    return a.blocks if isinstance(a, PackedArrays) else (a,)


def _shares(a, b):
    """Whether any block of ``a`` shares memory with any block of ``b``."""
    return any(np.shares_memory(x, y) for x in _blocks(a) for y in _blocks(b))


def _parts_packed(gamma_b, chi_b, geom, *, plan):
    """The packed branch of :func:`so_bond_pairing_vertex_parts`."""
    name = "so_bond_pairing_vertex_parts"
    if not (isinstance(gamma_b, PackedArrays) and isinstance(chi_b, PackedArrays)):
        raise ValueError("{}: gamma_b and chi_b must both be PackedArrays or both dense arrays"
                         .format(name))
    _check_convention(name, geom)
    nvol = _packed_nvol(name, gamma_b)
    sm = gamma_b.sector_map
    plan.require(name, nvol=nvol, ND=sm.ND, B=getattr(geom, "B", None), layout="sz_packed")
    _check_packed_map(name, gamma_b, geom)
    if not _same_sector_map(chi_b.sector_map, sm):
        raise ValueError("{}: the sector map of chi_b is not the sector map of gamma_b"
                         .format(name))
    if chi_b.lead_shape != gamma_b.lead_shape:
        raise ValueError("{}: the packed chi_b has leading shape {}; gamma_b has {}".format(
            name, chi_b.lead_shape, gamma_b.lead_shape))
    if _shares(chi_b, gamma_b):
        raise ValueError("{}: chi_b shares memory with gamma_b".format(name))
    for arg, arr in (("gamma_b", gamma_b), ("chi_b", chi_b)):
        if not all(np.isfinite(b).all() for b in arr.blocks):
            raise ValueError("{}: {} has non-finite entries".format(name, arg))

    v_bare = PackedArrays.zeros(sm, (nvol,))
    gamma_p(gamma_b, geom, out=v_bare)                       # Gamma^P, in place
    _refuse_exchange(name, *_exchange_residual(zip(v_bare.blocks, gamma_b.blocks)))
    v_fluct = PackedArrays.zeros(sm, (nvol,))
    for vbs, vfs, g, c in zip(v_bare.blocks, v_fluct.blocks, gamma_b.blocks, chi_b.blocks):
        ws = np.empty(g.shape[-2:], dtype=np.complex128)
        for q in range(nvol):
            vb = vbs[q]
            np.add(vb, g[q], out=vb)
            vb *= 0.25
            np.matmul(g[q], c[q], out=ws)
            np.matmul(ws, g[q], out=vfs[q])
            vfs[q] *= -1.0
    return v_bare, v_fluct


def so_bond_pairing_vertex_batch(gamma_b, Y, geom, *, plan, out):
    """The pairing vertex of one dressed frequency batch, written into
    ``out``:

    ``V[l](q) = 1/2 Gamma_B^sym(q) - Gamma_B(q) Y[l](q) Gamma_B(q)``,
    ``1/2 Gamma_B^sym = 0.25 (Gamma_B + Gamma_B^P)``,

    the per-slice sum of the two parts of
    :func:`so_bond_pairing_vertex_parts` with ``chi_b = Y[l]``. ``Y`` is a
    DRESSED batch (the output of :func:`so_bond_dress`), never the bare
    bubble.

    Dense layout: ``gamma_b`` is ``(nvol, ND, ND)``, ``Y`` and ``out`` are
    ``(nb, nvol, ND, ND)`` complex128. S_z-packed layout (a plan of layout
    ``"sz_packed"``): all three are :class:`PackedArrays` on the sector map
    of ``geom``, ``gamma_b`` of leading shape ``(nvol,)``, ``Y`` and
    ``out`` of leading shape ``(nb, nvol)``; the vertex is formed sector by
    sector.

    The boundary is checked against ``plan`` (:meth:`SOBondPlan.require`:
    ``nvol``, ``ND``, the geometry's ``B``, the layout and the batch length
    ``nb <= plan.dress_capacity``; the plan needs ``vertex_batch=True``)
    before anything is written. ``out`` is caller owned (for example a
    store slot's batch view) and is overwritten completely; it must share
    no memory with ``gamma_b`` or ``Y``, which are never written. The
    exchange image ``Gamma^P`` is assembled in ``out[0]`` and the work is
    done per ``(l, q)`` in a workspace of two ``ND x ND`` (packed: ``n_s x
    n_s``) matrices. A non-finite entry of ``gamma_b`` or ``Y`` is refused
    (``ValueError`` naming the argument), and ``Gamma_B`` must equal its
    exchange image as in :func:`so_bond_pairing_vertex_parts`; that check
    runs on the assembled image, so on its refusal ``out[0]`` holds
    ``Gamma^P`` and the other slices of ``out`` are untouched.

    Returns ``out``.
    """
    name = "so_bond_pairing_vertex_batch"
    if not isinstance(plan, SOBondPlan):
        raise ValueError("{}: plan must be a SOBondPlan (plan_so_bond), got {}".format(
            name, type(plan).__name__))
    packed = [isinstance(a, PackedArrays) for a in (gamma_b, Y, out)]
    if any(packed) and not all(packed):
        raise ValueError("{}: gamma_b, Y and out must all be PackedArrays or all dense arrays"
                         .format(name))
    _check_convention(name, geom)
    B = getattr(geom, "B", None)
    if all(packed):
        nvol = _packed_nvol(name, gamma_b)
        sm = gamma_b.sector_map
        nb = Y.lead_shape[0] if len(Y.lead_shape) == 2 else None
        if nb is None or Y.lead_shape[1] != nvol:
            raise ValueError("{}: the packed Y has leading shape {}; gamma_b needs (nb, {})"
                             .format(name, Y.lead_shape, nvol))
        plan.require(name, nvol=nvol, ND=sm.ND, B=B, layout="sz_packed", batch=nb)
        _check_packed_map(name, gamma_b, geom)
        for arg, arr in (("Y", Y), ("out", out)):
            if not _same_sector_map(arr.sector_map, sm):
                raise ValueError("{}: the sector map of {} is not the sector map of gamma_b"
                                 .format(name, arg))
        if out.lead_shape != Y.lead_shape:
            raise ValueError("{}: the packed out has leading shape {}; Y has {}".format(
                name, out.lead_shape, Y.lead_shape))
    else:
        if not isinstance(gamma_b, np.ndarray) or gamma_b.dtype != np.complex128 \
                or gamma_b.ndim != 3:
            raise ValueError("{}: gamma_b must be a complex128 (nvol, ND, ND) array".format(
                name))
        if not (isinstance(Y, np.ndarray) and Y.dtype == np.complex128 and Y.ndim == 4
                and Y.shape[1:] == gamma_b.shape):
            raise ValueError("{}: Y must be a complex128 (nb, nvol, ND, ND) array matching "
                             "gamma_b {}, got {}".format(name, gamma_b.shape, _describe(Y)))
        nvol, ND, nb = int(gamma_b.shape[0]), int(gamma_b.shape[1]), int(Y.shape[0])
        plan.require(name, nvol=nvol, ND=ND, B=B, layout="dense", batch=nb)
        _nd_of(gamma_b.shape, geom, name, "gamma_b")
        if not (isinstance(out, np.ndarray) and out.dtype == np.complex128
                and out.shape == Y.shape):
            raise ValueError("{}: out must be a complex128 array of shape {}, got {}".format(
                name, Y.shape, _describe(out)))
    for arg, arr in (("Y", Y), ("gamma_b", gamma_b)):
        if _shares(out, arr):
            raise ValueError("{}: out shares memory with {}".format(name, arg))
    for arg, arr in (("gamma_b", gamma_b), ("Y", Y)):
        if not all(np.isfinite(b).all() for b in _blocks(arr)):
            raise ValueError("{}: {} has non-finite entries".format(name, arg))

    # Gamma^P into out[0] (a view of the caller's out), checked, then
    # 1/4 (Gamma + Gamma^P) there; the slices l >= 1 read it before slice 0
    # is completed in place
    if all(packed):
        bare = PackedArrays(sm, tuple(b[0] for b in out.blocks))
    else:
        bare = out[0]
    gamma_p(gamma_b, geom, out=bare)
    _refuse_exchange(name, *_exchange_residual(zip(_blocks(bare), _blocks(gamma_b))))
    for bs, g in zip(_blocks(bare), _blocks(gamma_b)):
        for q in range(nvol):
            np.add(bs[q], g[q], out=bs[q])
            bs[q] *= 0.25
    for os_, ys, g in zip(_blocks(out), _blocks(Y), _blocks(gamma_b)):
        ws1 = np.empty(g.shape[-2:], dtype=np.complex128)
        ws2 = np.empty(g.shape[-2:], dtype=np.complex128)
        for l in range(nb - 1, -1, -1):              # slice 0 holds the bare part: last
            for q in range(nvol):
                np.matmul(g[q], ys[l, q], out=ws1)
                np.matmul(ws1, g[q], out=ws2)
                np.subtract(os_[0, q], ws2, out=os_[l, q])
    return out


# ---------------------------------------------------------------------------
# The pairing kernel on the pair-S_z sectors
# ---------------------------------------------------------------------------

_KERNEL = "PackedBondPairKernel"
_ROLL_AXES = (0, 1, 2)


def _rolled(arr, shift):
    """``arr`` rolled by ``shift`` on the three spatial axes (the array
    itself for a zero shift)."""
    return arr if shift == (0, 0, 0) else np.roll(arr, shift, axis=_ROLL_AXES)


class PackedBondPairKernel:
    r"""The bond-resolved pairing kernel of an S_z-packed vertex, applied one
    pair-S_z sector at a time.

    The operator is that of :class:`hwave.solver.eliashberg_bond.BondPairKernel`
    on the uniform grid (no instantaneous vertex):

    .. code-block:: text

        out_AD(k) = -(1/N) sum_{alpha, beta} e^{-ik.R_beta}
            sum_{k'} V_{alpha beta; ABCD}(k - k') e^{-ik'.R_alpha} F_BC(k'),
        F_BC(k) = sum_{EF} G2[B, E, C, F](k) phi_EF(k),

    with the same transforms (boson -> tau and ``ifftn`` for the vertex,
    fermion -> tau and ``ifftn`` for ``F``) and roll conventions. With S_z
    conserved, the pair sector ``p = S_z(B) + S_z(C)`` (pp order ``(+1, 0,
    -1)``) of the gap is mapped onto itself: ``V_{alpha beta; ABCD}`` is
    nonzero only for ``S_z(A) - S_z(B) = S_z(C) - S_z(D)``, i.e. ``S_z(A) +
    S_z(D) = S_z(B) + S_z(C)``, and ``G2`` must not couple two pair sectors
    (refused otherwise, exact zero test).

    Parameters
    ----------
    store : PackedBondBlockStore or PackedArrays
        The vertex: slot ``slot`` of the store, or packed arrays, of leading
        shape ``(nfreq, nvol)`` on ``sector_map`` (the uniform bosonic
        grid; ``nfreq = 1`` for the static kernel).
    G2 : ndarray
        The pair propagator ``(nd, nd, nd, nd, Nx, Ny, Nz, nfreq)``, indexed
        ``G2[B, E, C, F]`` as above (``1/beta`` included, as for
        ``BondPairKernel``).
    view : bond view
        ``n_channels`` (must be the sector map's ``B``) and ``delta_r``.
    sector_map : SectorMap
    plan : SOBondPlan
        Layout ``"sz_packed"``, the mesh's ``nvol`` and the map's ``ND``.
    spatial_shape, beta, nfreq, workers
        As for ``BondPairKernel`` (``beta`` is not applied: ``1/beta`` lives
        in ``G2``).
    xp, residency
        ``numpy`` and ``"host"`` only.
    host_cap : float, optional
        Refuse a kernel whose :meth:`memory_table` total (the build rows
        and the apply rows together) exceeds it; ``None`` (the default) is
        the only unlimited value -- NaN, an infinity or a negative cap is
        refused.

    The vertex is hoisted once: per pp sector ``p`` (``n_p`` gap pairs) the
    entries ``V^p[alpha, beta][i, j] = V_{alpha beta; A B C D}`` (output pair
    ``i = (A, D)``, input pair ``j = (B, C)``) are gathered from the packed
    ph-sector blocks through the map's pp gather tables and transformed to
    ``(r, tau)``, one ``(alpha, beta)`` pair at a time. Every packed entry
    lies in exactly one pp sector, so the hoisted blocks hold the packed
    vertex once (``hoisted_blocks`` of :meth:`memory_table`).

    The gap is flat in the order of ``gap_shape`` (``(nd, nd, Nx, Ny, Nz)``
    for ``nfreq = 1``, else ``(nd, nd, Nx, Ny, Nz, nfreq)``; the same flat
    order as ``BondPairKernel``'s gap). :meth:`sector_op` is the sector
    API.
    """

    def __init__(self, store, G2, view, sector_map, *, plan, spatial_shape, beta, nfreq=1,
                 slot="V", xp=np, residency="host", host_cap=None, workers=1):
        if xp is not np:
            raise ValueError("{}: the S_z-packed kernel runs on the host (xp must be numpy), "
                             "got {}".format(_KERNEL, _module_name(xp)))
        if residency != "host":
            raise ValueError("{}: residency must be 'host', got {!r}".format(_KERNEL, residency))
        if host_cap is not None:
            try:
                cap_ok = float(host_cap) >= 0.0 and np.isfinite(float(host_cap))
            except (TypeError, ValueError):
                cap_ok = False
            if not cap_ok:
                raise ValueError("{}: host_cap must be None (unlimited) or a finite number of "
                                 "bytes >= 0, got {!r}".format(_KERNEL, host_cap))
        if not isinstance(sector_map, SectorMap):
            raise ValueError("{}: sector_map must be a SectorMap, got {}".format(
                _KERNEL, type(sector_map).__name__))
        if not isinstance(plan, SOBondPlan):
            raise ValueError("{}: plan must be a SOBondPlan (plan_so_bond), got {}".format(
                _KERNEL, type(plan).__name__))
        if plan.layout != "sz_packed":
            raise ValueError("{}: the admission plan has layout = {!r}; the packed kernel needs "
                             "'sz_packed'".format(_KERNEL, plan.layout))
        shape = tuple(int(x) for x in spatial_shape)
        if len(shape) != 3 or min(shape) < 1:
            raise ValueError("{}: spatial_shape must be three positive integers, got {!r}"
                             .format(_KERNEL, spatial_shape))
        nvol = int(np.prod(shape))
        nfreq = _int_at_least(_KERNEL, "nfreq", nfreq, 1)
        if plan.nvol != nvol:
            raise ValueError("{}: the admission plan has nvol = {} but spatial_shape {} has {} "
                             "k points".format(_KERNEL, plan.nvol, shape, nvol))
        if plan.ND != sector_map.ND:
            raise ValueError("{}: the admission plan has ND = {} but the sector map has ND = {}"
                             .format(_KERNEL, plan.ND, sector_map.ND))
        B = int(view.n_channels)
        if B != sector_map.B:
            raise ValueError("{}: view.n_channels = {} but the sector map has B = {}".format(
                _KERNEL, B, sector_map.B))
        delta_r = [tuple(int(x) for x in r) for r in view.delta_r]
        if len(delta_r) != B:
            raise ValueError("{}: view.delta_r has {} entries but view.n_channels is {}".format(
                _KERNEL, len(delta_r), B))
        packed = self._vertex(store, slot, sector_map, nfreq, nvol)
        nd = sector_map.nd
        G2 = np.asarray(G2)
        want = (nd,) * 4 + shape + (nfreq,)
        if G2.shape != want:
            raise ValueError("{}: G2 must have shape {}, got {}".format(_KERNEL, want, G2.shape))

        gm = _consumer_map("kernel_gather", sector_map)
        pp_of = np.empty(nd * nd, dtype=np.int64)
        for p, pairs in zip(gm.pp_order, gm.pp_pairs):
            pp_of[pairs] = p
        pp_of = pp_of.reshape(nd, nd)
        off = pp_of[:, None, :, None] != pp_of[None, :, None, :]       # G2[i, l, j, m]
        worst = 0.0
        for i in range(nd):                     # one output mode at a time: small temporaries
            vals = G2[i][off[i]]
            if vals.size and np.any(vals != 0):
                v = float(np.max(np.abs(vals)))
                if not np.isfinite(v):
                    worst = v
                    break
                worst = max(worst, v)
        if worst != 0.0:
            raise ValueError(
                "PackedBondPairKernel: G2 couples pair S_z sectors (max |entry| = %.3e); the "
                "packed kernel needs a Green function built per spin sector" % worst)

        self.sector_map = sector_map
        self.shape = shape
        self.nvol, self.nfreq, self.B, self.nd = nvol, nfreq, B, nd
        self.beta, self.workers = float(beta), workers
        self.delta_r = delta_r
        self.gap_shape = (nd, nd) + shape + ((nfreq,) if nfreq > 1 else ())
        self.pp_sectors = tuple(int(p) for p in gm.pp_order)
        self._pairs = tuple(np.asarray(pairs) for pairs in gm.pp_pairs)
        self._table = estimate_packed_kernel_workspace_bytes(
            B, [len(pairs) for pairs in self._pairs], nvol, nfreq)
        if host_cap is not None:
            need = int(sum(self._table.values()))
            if need > float(host_cap):
                raise ValueError("{}: the kernel needs {} bytes ({}) but host_cap is {} bytes"
                                 .format(_KERNEL, need, ", ".join(
                                     "{} = {}".format(k, v) for k, v in self._table.items()),
                                     host_cap))

        # -- hoisting: per pp sector, the pair-pair blocks in (r, tau) --------
        self._G2, self._V = [], []
        for s, pairs in enumerate(self._pairs):
            o1, o2 = np.divmod(pairs, nd)                               # pair (o1, o2)
            g2 = G2[o1[:, None], o1[None, :], o2[:, None], o2[None, :]]
            self._G2.append(np.ascontiguousarray(np.moveaxis(g2, (0, 1), (-2, -1))))
            self._V.append(self._hoist(packed, gm.pp_gather[s], len(pairs)))
        self._ops = {}

    @staticmethod
    def _vertex(store, slot, sector_map, nfreq, nvol):
        if isinstance(store, PackedBondBlockStore):
            if not _same_sector_map(store.sector_map, sector_map):
                raise ValueError("{}: the store's sector map is not the kernel's".format(_KERNEL))
            packed = store.packed(slot)
        elif isinstance(store, PackedArrays):
            if not _same_sector_map(store.sector_map, sector_map):
                raise ValueError("{}: the vertex's sector map is not the kernel's".format(
                    _KERNEL))
            packed = store
        else:
            raise ValueError("{}: the vertex must be a PackedBondBlockStore or PackedArrays, "
                             "got {}".format(_KERNEL, type(store).__name__))
        if packed.lead_shape != (nfreq, nvol):
            raise ValueError("{}: the vertex must have leading shape (nfreq, nvol) = {}, got {}"
                             .format(_KERNEL, (nfreq, nvol), packed.lead_shape))
        return packed

    def _hoist(self, packed, gather, n_p):
        """``V^p`` ``(B, B, Nx, Ny, Nz, nfreq, n_p, n_p)`` in ``(r, tau)``."""
        nx, ny, nz = self.shape
        B, nfreq = self.B, self.nfreq
        sel = []
        for sigma in self.sector_map.sigmas:
            ii, jj = np.nonzero(gather.sector == sigma)
            if ii.size:
                sel.append((packed[int(sigma)], ii, jj))
        out = np.empty((B, B, nx, ny, nz, nfreq, n_p, n_p), dtype=np.complex128)
        raw = np.empty((nfreq, self.nvol, n_p, n_p), dtype=np.complex128)
        for a in range(B):
            for b in range(B):
                for blk, ii, jj in sel:
                    raw[:, :, ii, jj] = blk[:, :, gather.row_local[a, ii, jj],
                                            gather.col_local[b, ii, jj]]
                bt = _ms.boson_to_tau(raw, axis=0).reshape(nfreq, nx, ny, nz, n_p, n_p)
                br = _bk.spatial_ifftn(bt, axes=(1, 2, 3), workers=self.workers)
                del bt
                out[a, b] = np.moveaxis(br, 0, 3)
                del br
        return out

    # -- sectors -----------------------------------------------------------
    def _pos(self, p):
        if not _is_int(p) or int(p) not in self.pp_sectors:
            raise ValueError("{}: pp_sector must be +1, 0 or -1, got {!r}".format(_KERNEL, p))
        return self.pp_sectors.index(int(p))

    def pp_index(self, p):
        """The ascending flat indices into :attr:`gap_shape` of the
        components ``Delta_BC`` with ``S_z(B) + S_z(C) = p`` (pair-major,
        then the momentum and frequency grid)."""
        pairs = self._pairs[self._pos(p)]
        grid = self.nvol * self.nfreq
        return (pairs[:, None] * grid + np.arange(grid)[None, :]).ravel()

    def sector_shape(self, p):
        """``(n_p, Nx, Ny, Nz)`` (static) or ``(n_p, Nx, Ny, Nz, nfreq)``."""
        n_p = len(self._pairs[self._pos(p)])
        return (n_p,) + self.shape + ((self.nfreq,) if self.nfreq > 1 else ())

    def sector_op(self, p):
        """The :class:`hwave.sc.SectorOp` of pp sector ``p``: the kernel on
        the sector vector (``g[pp_index(p)]``) with ``matvec`` and
        ``rmatvec``, and the fermion-antisymmetric basis of the sector
        (``_SOAntisymmetricBasis.from_subset`` on :attr:`gap_shape`, with
        ``(A, B, k) -> (B, A, -k)`` for ``nfreq = 1`` and ``(A, B, k, n) ->
        (B, A, -k, nfreq - 1 - n)`` otherwise). Built on the first call and
        cached: later calls return the same object."""
        s = self._pos(p)
        op = self._ops.get(s)
        if op is not None:
            return op
        from hwave import sc
        if self.nfreq == 1:
            reverse = sc._reverse_k_and_orbital
        else:
            from . import eliashberg_dynamic
            reverse = eliashberg_dynamic._reverse_kw_and_orbital
        basis = sc._SOAntisymmetricBasis.from_subset(self.gap_shape, self.pp_index(p),
                                                      reverse=reverse)
        from scipy.sparse.linalg import LinearOperator
        n = basis.n
        operator = LinearOperator((n, n), matvec=lambda x: self._apply(s, x),
                                  rmatvec=lambda y: self._apply(s, y, adjoint=True),
                                  dtype=complex)
        op = sc.SectorOp(self.pp_sectors[s], operator, basis, True)
        self._ops[s] = op
        return op

    def _apply(self, s, x, adjoint=False):
        """The kernel (or its adjoint) of sector position ``s`` on a sector
        vector ``(n_p * nvol * nfreq,)``.

        The adjoint reverses every step: the conjugate-transposed pair
        blocks, the opposite rolls, and the transforms' adjoints (``ifftn``
        and ``fftn``, fermion -> tau and tau -> fermion trade places; their
        scale factors cancel in the product)."""
        n_p = len(self._pairs[s])
        grid = self.shape + (self.nfreq,)
        x = np.asarray(x)
        if x.size != n_p * self.nvol * self.nfreq:
            raise ValueError("{}: the sector vector of pp sector {} must have {} entries, got {}"
                             .format(_KERNEL, self.pp_sectors[s], n_p * self.nvol * self.nfreq,
                                     x.size))
        xg = np.moveaxis(x.reshape((n_p,) + grid), 0, -1)      # (Nx, Ny, Nz, nfreq, n_p)
        G2, V = self._G2[s], self._V[s]
        if adjoint:
            F = xg
        else:
            F = np.einsum("...ij,...j->...i", G2, xg)
        F_rt = _bk.spatial_ifftn(_ms.fermion_to_tau(F, axis=3), axes=_ROLL_AXES,
                                 workers=self.workers)
        del F
        # the adjoint accumulates the complex conjugate of each leg,
        # conj(-V^H f) = -V^T conj(f), and conjugates the sums once
        acc = [None] * self.B
        for u in range(self.B):                 # the input leg: alpha, or beta (adjoint)
            R = self.delta_r[u]
            if adjoint:
                F_u = np.conj(_rolled(F_rt, tuple(-r for r in R)))
            else:
                F_u = _rolled(F_rt, R)
            for w in range(self.B):
                if adjoint:
                    P = np.einsum("...ji,...j->...i", V[w, u], F_u)
                else:
                    P = np.einsum("...ij,...j->...i", V[u, w], F_u)
                if acc[w] is None:
                    acc[w] = np.negative(P, out=P)
                else:
                    acc[w] -= P
                del P
            del F_u
        del F_rt
        out_rt = None
        for w in range(self.B):                 # the output leg
            R = self.delta_r[w]
            if adjoint:
                np.conj(acc[w], out=acc[w])
            t = _rolled(acc[w], tuple(-r for r in R) if adjoint else R)
            acc[w] = None
            if out_rt is None:
                out_rt = t
            else:
                out_rt += t
            del t
        out = _ms.tau_to_fermion(_bk.spatial_fftn(out_rt, axes=_ROLL_AXES,
                                                  workers=self.workers), axis=3)
        del out_rt
        if adjoint:
            out = np.conj(np.einsum("...ji,...j->...i", G2, np.conj(out)))
        return np.moveaxis(out, -1, 0).ravel()

    # -- the full gap -----------------------------------------------------
    def _full(self, phi, adjoint, where):
        n = int(np.prod(self.gap_shape))
        phi = np.asarray(phi).ravel()
        if phi.size != n:
            raise ValueError("{}.{}: the gap must have {} entries, got {}".format(
                _KERNEL, where, n, phi.size))
        grid = self.nvol * self.nfreq
        src = phi.reshape(self.nd * self.nd, grid)
        out = np.zeros((self.nd * self.nd, grid), dtype=np.complex128)
        for s, pairs in enumerate(self._pairs):
            out[pairs] = self._apply(s, src[pairs].ravel(), adjoint).reshape(len(pairs), grid)
        return out.ravel()

    def matvec(self, phi_full):
        """The kernel on a flat full gap (the flat order of ``gap_shape``,
        as ``BondPairKernel.matvec``)."""
        return self._full(phi_full, False, "matvec")

    def rmatvec(self, y_full):
        """The adjoint kernel on a flat full gap."""
        return self._full(y_full, True, "rmatvec")

    @property
    def hoisted_nbytes(self):
        """Bytes of the hoisted vertex blocks of all pp sectors."""
        return int(sum(v.nbytes for v in self._V))

    def memory_table(self):
        """The workspace rows of the kernel
        (:func:`estimate_packed_kernel_workspace_bytes`)."""
        return dict(self._table)
