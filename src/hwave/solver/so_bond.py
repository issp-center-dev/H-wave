"""Spin-orbital bond-resolved two-body objects: the slot convention and the
slot geometry, the bond topology of a compiled interaction, the admission
plan (:func:`plan_so_bond`), the spinful bond vertex (:func:`so_bond_vertex`)
and its exchange image (:func:`gamma_p`), the bond bubble
(:func:`so_bond_bubble`), its RPA dressing (:func:`so_bond_dress`) and the
static pairing vertex (:func:`so_bond_pairing_vertex_parts`).

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

import logging
import re
from dataclasses import dataclass
from typing import NamedTuple, Optional

import numpy as np

from . import bond_channels, bubble, flex_bond
from .hartree_fock import NonFiniteError
from .so_records import CompiledInteraction, _displacements

__all__ = ["SlotConvention", "SLOT_CONVENTION", "BondSlotGeometry", "slot_geometry",
           "so_bond_topology", "SOBondCaps", "SOBondPlan", "plan_so_bond", "so_bond_vertex",
           "gamma_p", "so_bond_bubble", "so_bond_dress", "so_bond_pairing_vertex_parts"]

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
# The admission plan
# ---------------------------------------------------------------------------

_GIB = float(2 ** 30)
_ENTRY = 16                     # bytes per complex128 entry
#: the consumers of a plan, in canonical order
CONSUMERS = ("vertex", "chi_static", "chi_dynamic", "pairing_static")
#: what each consumer needs besides itself (closed transitively)
_IMPLIES = {"vertex": (), "chi_static": ("vertex",), "chi_dynamic": ("vertex",),
            "pairing_static": ("chi_static",)}
_LAYOUT_REASON = "dense layout (the S_z-sector packed layout is not available on this path)"


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
}


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
        ``"host"``, ``"dense"`` (with its reason) and ``"host"``.
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
        observed ``freq``, ``"static"`` or ``"all"``), ``"so_bond_dress"`` or
        ``"so_bond_pairing_vertex_parts"``. The observed fields are the
        dimensions visible at that boundary:

        * vertex: ``nvol``, ``B``, ``nd``, ``ND``, ``layout``;
        * static bubble: those and ``nmat``; every-slice bubble: those,
          ``nmat``, ``placement`` and the store ``slots`` (which must hold
          ``"chibar"``);
        * dressing: ``nvol``, ``ND`` and the ``batch`` length
          (``1 <= batch <= dress_capacity``);
        * pairing parts: ``nvol``, ``ND`` and the geometry's ``B`` (with
          ``ND == plan.ND = plan.B * plan.nd**2`` this is
          ``ND == B * plan.nd**2``).

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
        if set(observed) != set(fields):
            raise ValueError("SOBondPlan.require: {} observes the fields {}, got {}".format(
                helper, ", ".join(fields), ", ".join(sorted(observed))))
        for f in fields:
            got = observed[f]
            if f == "batch":
                if not _is_int(got) or int(got) <= 0:
                    raise ValueError("{}: the batch length must be a positive integer, got {!r}"
                                     .format(helper, got))
                if int(got) > self.dress_capacity:
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


def _plan_rows(consumers, *, nmat, nvol, nd, ND, dress_capacity, keep_chi):
    """The row table: ``(phase, name, bytes)``, one entry per logical
    allocation, ``"resident"`` first, then each closed consumer's phase."""
    full = _ENTRY * nvol * ND * ND                          # one (nvol, ND, ND) array
    P = nd * nd
    # the bond bubble's per-shift workspace per k point: its P^2-sized and
    # P-sized frequency-grid buffers (bubble.BOND_BUBBLE_N4_BUFFERS and
    # BOND_BUBBLE_N2_BUFFERS, the counts the bubble's memory contract and
    # sc's estimator use) plus 8 P per-k temporaries
    pair_block = _ENTRY * nvol * (bubble.BOND_BUBBLE_N4_BUFFERS * nmat * P * P
                                  + bubble.BOND_BUBBLE_N2_BUFFERS * nmat * P + 8 * P)
    dress = 7 * _ENTRY * dress_capacity * nvol * ND * ND
    rows = [("resident", "gamma_b", full)]
    if "chi_static" in consumers:
        rows += [("chi_static", "chibar_static", full), ("chi_static", "bubble_static_temp", full),
                 ("chi_static", "bubble_pair_block", pair_block),
                 ("chi_static", "dress_batch", dress), ("chi_static", "chi_static", full)]
    if "chi_dynamic" in consumers:
        rows += [("chi_dynamic", "chibar", _ENTRY * nmat * nvol * ND * ND)]
        if keep_chi:
            rows += [("chi_dynamic", "chi", _ENTRY * nmat * nvol * ND * ND)]
        rows += [("chi_dynamic", "bubble_pair_block", pair_block),
                 ("chi_dynamic", "dress_batch", dress)]
    if "pairing_static" in consumers:
        rows += [("pairing_static", "chibar_static", full),
                 ("pairing_static", "bubble_static_temp", full),
                 ("pairing_static", "bubble_pair_block", pair_block),
                 ("pairing_static", "chi_static", full),
                 ("pairing_static", "dress_batch", dress), ("pairing_static", "v_bare", full),
                 ("pairing_static", "v_fluct", full),
                 ("pairing_static", "pairing_q_workspace", 2 * _ENTRY * ND * ND)]
    return tuple((ph, name, int(b)) for ph, name, b in rows)


def _row_table(rows):
    return "rows (phase  name  bytes):\n" + "\n".join(
        "  {}  {}  {}".format(ph, name, b) for ph, name, b in rows)


def plan_so_bond(*, nmat, nvol, B, nd, half_grid=False, consumers, caps,
                 placement_request="host", nb=None, keep_chi=False, write_chibar=False):
    """The admission plan of the spin-orbital bond helpers (dense, host).

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
    """
    if half_grid:
        raise ValueError("plan_so_bond: the half bosonic grid is not available on the "
                         "spin-orbital bond path")
    if placement_request != "host":
        raise ValueError("plan_so_bond: placement_request must be 'host' on the spin-orbital "
                         "bond path, got {!r}".format(placement_request))
    if not isinstance(caps, SOBondCaps):
        raise ValueError("plan_so_bond: caps must be a SOBondCaps, got {!r}".format(caps))
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
    ND = B * nd * nd
    dress_capacity = nb if "chi_dynamic" in closed else 1
    keep_chi, write_chibar = bool(keep_chi), bool(write_chibar)
    rows = _plan_rows(closed, nmat=nmat, nvol=nvol, nd=nd, ND=ND,
                      dress_capacity=dress_capacity, keep_chi=keep_chi)
    resident = sum(b for ph, _, b in rows if ph == "resident")
    peaks = [(1.25 * (resident + sum(b for ph, _, b in rows if ph == c)), c) for c in closed]
    peak, phase = max(peaks, key=lambda p: p[0])
    slots = (("chibar",) + (("chi",) if keep_chi else ())) if "chi_dynamic" in closed else ()
    cap = caps.host_bytes
    small = _GIB if cap is None else min(_GIB, float(cap) / 8.0)
    reason = ("spin-orbital bond plan: consumers {}; dense layout, host placement; nmat = {}, "
              "nvol = {}, B = {}, nd = {}, ND = {}, nb = {}, dress_capacity = {}; estimated "
              "peak {:.3f} GiB in phase {}; host cap {}".format(
                  ", ".join(closed), nmat, nvol, B, nd, ND, nb, dress_capacity, peak / _GIB,
                  phase, "unlimited" if cap is None else "{:.3f} GiB".format(cap / _GIB)))
    if cap is not None and peak > cap:
        solver = "Eliashberg" if "pairing_static" in closed else "FLEX"
        raise ValueError(
            "spin-orbital bond {} needs an estimated peak of {:.2f} GiB in phase {} (B = {} "
            "channels, ND = B·(2·norb_phys)² = {}, nvol = {}, Nmat = {}; one store slot = "
            "16·Nmat·nvol·ND² bytes = {:.2f} GiB) but the cap is {:.2f} GiB (caps.host_bytes). "
            "Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce Nmat or the k grid, or prune "
            "bond shells.\n{}".format(
                solver, peak / _GIB, phase, B, ND, nvol, nmat,
                _ENTRY * nmat * nvol * ND * ND / _GIB, cap / _GIB, _row_table(rows)))
    logger.info(reason)
    return SOBondPlan(
        nmat=nmat, nslices=nmat, nvol=nvol, B=B, nd=nd, ND=ND, nb=nb,
        dress_capacity=dress_capacity, consumers=closed, placement="host", layout="dense",
        layout_reason=_LAYOUT_REASON, kernel_residency="host", slots=slots, aliases=(),
        keep_chi=keep_chi, write_chibar=write_chibar, small_array_allowed=False,
        small_array_bytes=small, rows=rows, peak_host=peak, peak_device=None,
        cap_host=None if cap is None else float(cap), cap_host_raw=caps.host_raw,
        cap_device=caps.device_bytes, reason=reason)


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
    refused."""
    _check_convention("_file_records", geom)
    parts = tuple(parts)
    if not parts or set(parts) - {"direct", "crossed"}:
        raise ValueError("_file_records: parts must be a non-empty subset of ('direct', "
                         "'crossed'), got {!r}".format(parts))
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
    """
    if not isinstance(compiled, CompiledInteraction):
        raise ValueError("so_bond_vertex: compiled must be a CompiledInteraction "
                         "(so_records.compile_interaction), got {}".format(
                             type(compiled).__name__))
    if not isinstance(plan, SOBondPlan):
        raise ValueError("so_bond_vertex: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    geom = slot_geometry(topo, spatial_shape)
    _check_convention("so_bond_vertex", geom)
    nd = 2 * int(compiled.norb_phys)
    ND = geom.B * nd * nd
    plan.require("so_bond_vertex", nvol=geom.nvol, B=geom.B, nd=nd, ND=ND, layout="dense")
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
    is done per q in two ``ND x ND`` matrices. Returns ``out``."""
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
    geom = slot_geometry(view, spatial_shape)
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
    observed = dict(nvol=nvol, B=geom.B, nd=nd, ND=ND, layout="dense", nmat=nmat)
    if freq == "all":
        observed.update(placement=getattr(sink, "resident", None),
                        slots=tuple(getattr(sink, "slot_names", ())))
    plan.require("so_bond_bubble", freq=freq, **observed)
    bond_set = bond_channels.BondSetView(
        bond_channels.BondTopology(geom.delta_r, geom.reverse, {}))
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


def so_bond_dress(X, gamma_b, *, plan, l0, spatial_shape, cond_tol=bond_channels._BOND_COND_FLOOR):
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
    """
    if not isinstance(plan, SOBondPlan):
        raise ValueError("so_bond_dress: plan must be a SOBondPlan (plan_so_bond), got {}"
                         .format(type(plan).__name__))
    shape = tuple(np.shape(X))
    if len(shape) != 4 or shape[2] != shape[3]:
        raise ValueError("so_bond_dress: X must be a batch (nb, nvol, ND, ND), got shape {}"
                         .format(shape))
    plan.require("so_bond_dress", nvol=shape[1], ND=shape[2], batch=shape[0])
    for arg, arr in (("X", X), ("gamma_b", gamma_b)):
        if not np.isfinite(arr).all():
            raise NonFiniteError("so_bond_dress: {} has non-finite entries".format(arg))
    chi, _ = bond_channels.dress_batch(X, gamma_b, "spinful", l0=l0, nmat=plan.nmat,
                                       spatial_shape=spatial_shape, cond_tol=cond_tol)
    return chi


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

    Returns ``(V_bare, V_fluct)``.
    """
    name = "so_bond_pairing_vertex_parts"
    if not isinstance(plan, SOBondPlan):
        raise ValueError("{}: plan must be a SOBondPlan (plan_so_bond), got {}".format(
            name, type(plan).__name__))
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
