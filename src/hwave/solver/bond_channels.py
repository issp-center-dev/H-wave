"""Bond-channel topology resolution for the bond-resolved Eliashberg extension.

This module builds the single source-of-truth ``ResolvedInteractionSet`` from
the raw R-resolved ``CoulombInter`` dict (as returned by
``hwave.qlmsio.wan90.read_w90``). See
``docs/superpowers/specs/2026-07-25-bond-channels-eliashberg-design.md``
sections 3.0-3.2 for the authoritative algorithm this implements.

Only channel-topology bookkeeping lives here: reversal/Hermiticity closure,
declared-topology inclusion (a channel that appears in the input dict is kept
even at V=0), shell sorting/cutoff, and the immutable ``v_bond`` carrier.
Vertex/bubble construction and any complex-input rejection live elsewhere
(later tasks / ``sc.py`` top level) -- this module stays general.
"""

from dataclasses import dataclass, field
import logging

import numbers
import numpy as np
from scipy.sparse.linalg import LinearOperator

from . import backend as _bk
from . import bubble as _bubble

logger = logging.getLogger("qlms.solver.bond_channels")


@dataclass(frozen=True)
class ResolvedInteractionSet:
    """Reversal-closed, shell-sorted, declared-topology bond channel set.

    Attributes
    ----------
    delta_r : tuple[tuple[int, int, int], ...]
        Bravais-cell displacement per channel. ``delta_r[0] == (0, 0, 0)``.
    v_bond : tuple[np.ndarray, ...]
        ``v_bond[m]`` is a ``(norb, norb)`` complex array holding
        ``V_ab(delta_r[m])``. ``v_bond[0]`` is all-zeros (the on-site Fock
        value is not carried as a bond channel; see spec S3.2 "Delta r=0
        handling"). Every array is marked read-only.
    v_onsite : np.ndarray, shape (norb, norb), complex
        The on-site (``Delta r = 0``) ``V_ab(0)`` block, Hermiticity-closed
        by the SAME reverse/consistency synthesis ``resolve_interactions``
        applies to the ``m != 0`` bond channels: ``V_ab(0) = conj(V_ba(0))``
        when only one orientation is declared; a ``ValueError`` when both are
        declared and disagree beyond ``reverse_atol``/``reverse_rtol``. This
        is the single source of truth for the on-site Fock block consumed by
        ``sc.py._build_bond_m0_blocks`` (review fix I-1: previously ``sc.py``
        read the raw ``CoulombInter`` dict directly for ``R=0``, bypassing
        this closure). Marked read-only. Separate array from ``v_bond[0]``,
        which stays all-zeros per its existing contract.
    reverse : tuple[int, ...]
        ``reverse[m]`` is the channel index of ``-delta_r[m]``.
    n_channels : int
        ``len(delta_r)`` (``B = Z + 1`` in the spec's notation).
    dropped : tuple
        Provenance records (see module docstring) of shells/entries removed
        by ``bond_max_shells`` or synthesized-but-negligible reverse partners.
    """

    delta_r: tuple
    v_bond: tuple
    reverse: tuple
    n_channels: int
    v_onsite: np.ndarray = None
    dropped: tuple = field(default_factory=tuple)


def _shell_length(irvec, lattice_vectors):
    disp = np.asarray(irvec, dtype=float) @ np.asarray(lattice_vectors, dtype=float)
    return float(np.linalg.norm(disp))


def _complete_onsite_hermitian(onsite, norb, reverse_atol, reverse_rtol):
    """Hermiticity-close the on-site (``Delta r = 0``) block.

    Applies the SAME synthesis/consistency rule ``resolve_interactions`` uses
    for the ``m != 0`` bond channels (``V_ab(R) = conj(V_ba(-R))``), specialized
    to ``R = -R = 0``: ``V_ab(0) = conj(V_ba(0))``. If only one orientation is
    declared, the other is synthesized; if both are declared and disagree
    beyond the scale-aware tolerance, raise (spec S3.2 "Reverse closure /
    Hermiticity", review fix I-1).

    Parameters
    ----------
    onsite : dict
        ``{(a, b): value}`` -- the raw ``Delta r = 0`` ``CoulombInter``
        entries (already summed if duplicate keys collapsed upstream).
    norb : int
    reverse_atol, reverse_rtol : float

    Returns
    -------
    np.ndarray, shape (norb, norb), complex
    """
    mat = np.zeros((norb, norb), dtype=complex)
    keys = set(onsite.keys()) | {(b, a) for (a, b) in onsite.keys()}
    for (a, b) in keys:
        if a == b:
            v = onsite.get((a, b))
            if v is not None:
                # R = -R = 0 makes the Hermiticity rule V_ab(0) = conj(V_ba(0))
                # read V_aa(0) = conj(V_aa(0)) on the diagonal, i.e. REAL. An
                # imaginary diagonal cannot be closed -- storing it verbatim
                # would return a non-Hermitian block from a function that
                # promises a Hermiticity-closed one -- so anything beyond the
                # same scale-aware reversal tolerance used for the off-diagonal
                # pairs is refused, and round-off below it is projected away so
                # the result is EXACTLY real.
                tol = reverse_atol + reverse_rtol * abs(v)
                if abs(np.imag(v)) > tol:
                    raise ValueError(
                        "resolve_interactions: on-site V_{}{}(0) = {} has a "
                        "nonzero imaginary part ({}); Hermiticity requires "
                        "V_aa(0) = conj(V_aa(0)), i.e. a REAL diagonal "
                        "(tolerance atol={}, rtol={})".format(
                            a, a, v, np.imag(v), reverse_atol, reverse_rtol,
                        )
                    )
                mat[a, b] = complex(np.real(v), 0.0)
            continue
        v_ab = onsite.get((a, b))
        v_ba = onsite.get((b, a))
        if v_ab is not None and v_ba is not None:
            v_synth = np.conj(v_ba)
            tol = reverse_atol + reverse_rtol * abs(v_ab)
            if abs(v_ab - v_synth) > tol:
                raise ValueError(
                    "resolve_interactions: on-site V_{}{}(0) = {} disagrees "
                    "with conj(V_{}{}(0)) = {} beyond tolerance "
                    "(atol={}, rtol={})".format(
                        a, b, v_ab, b, a, v_synth, reverse_atol, reverse_rtol,
                    )
                )
            # Accepted within tolerance: PROJECT onto the exactly conjugate
            # pair instead of storing both declared values verbatim, which
            # would leave an O(tol) asymmetry in a block this function promises
            # is Hermiticity-closed (review fix R2-2). The symmetric average is
            # taken independently for (a,b) and (b,a); IEEE addition is
            # commutative and negation exact, so the two results are EXACT
            # conjugates of one another.
            mat[a, b] = 0.5 * (v_ab + v_synth)
        elif v_ab is not None:
            mat[a, b] = v_ab
        else:
            mat[a, b] = np.conj(v_ba)
    return mat


def resolve_interactions(
    coulomb_inter,
    lattice_vectors,
    norb,
    *,
    bond_max_shells=None,
    tol_incl=1e-12,
    reverse_atol=1e-10,
    reverse_rtol=1e-8,
):
    """Resolve the raw R-resolved ``CoulombInter`` dict into a
    ``ResolvedInteractionSet``.

    Parameters
    ----------
    coulomb_inter : dict
        ``{(irvec, orbvec): value}`` as returned by
        ``hwave.qlmsio.wan90.read_w90`` (``irvec`` a length-3 int tuple,
        ``orbvec = (a, b)`` 0-based orbital indices, ``value`` complex).
    lattice_vectors : array_like, shape (3, 3)
        Real-space lattice vectors (rows). Used only for Euclidean bond
        length (shell sorting / ``bond_max_shells`` cutoff).
    norb : int
        Number of orbitals; each ``v_bond[m]`` is ``(norb, norb)``.
    bond_max_shells : int or None, optional
        Keep shells ``0..bond_max_shells`` (shell 0 = Delta r = (0,0,0)).
        ``None`` keeps all declared shells. ``bond_max_shells=0`` while any
        *declared* nonzero inter-site ``V`` is present raises ``ValueError``
        (see spec S3.2 "Cutoff") -- "nonzero" here means literally ``!= 0``,
        independent of ``tol_incl``. A negative, non-integral or non-finite
        value also raises ``ValueError`` rather than being clamped/truncated.
    tol_incl : float, optional
        Magnitude tolerance below which a *synthesized* reverse partner is
        treated as negligible for provenance purposes. Declared entries are
        always kept regardless of magnitude (declared topology, not
        magnitude-derived; see spec S3.0/S3.2).
    reverse_atol, reverse_rtol : float, optional
        Scale-aware tolerance for consistency of an explicitly-declared pair
        ``V_ab(R)`` and ``V_ba(-R)`` (Hermiticity check); disagreement beyond
        ``atol + rtol * |value|`` raises ``ValueError``. A pair that passes is
        then PROJECTED onto its symmetric average, so the returned
        ``v_bond``/``v_onsite`` satisfy the reversal/Hermiticity closure
        ``v_bond[m] == v_bond[reverse[m]].conj().T`` and
        ``v_onsite == v_onsite.conj().T`` EXACTLY (to the last bit), not merely
        within the tolerance (review fix R2-2).

    Returns
    -------
    ResolvedInteractionSet
    """
    norb = int(norb)
    dropped = []

    # --- Step 1: group entries by irvec (declared topology, not magnitude) ---
    # A duplicate (irvec, orbvec) key cannot appear here since a Python dict
    # has already deduplicated -- that check lives at the wan90.read_w90
    # parser (spec S3.2 "Duplicates").
    by_irvec = {}
    onsite = {}
    for (irvec, orbvec), value in coulomb_inter.items():
        irvec = tuple(int(x) for x in irvec)
        a, b = int(orbvec[0]), int(orbvec[1])
        if irvec == (0, 0, 0):
            # Delta r=0 entries route to the Hartree/on-site block only, never
            # a Fock bond channel (spec S3.2 "Delta r=0 handling"). Collected
            # here and Hermiticity-closed below (review fix I-1) rather than
            # discarded -- the closure is this module's single source of
            # truth, consumed by sc.py._build_bond_m0_blocks.
            onsite[(a, b)] = complex(value)
            continue
        by_irvec.setdefault(irvec, {})[(a, b)] = complex(value)

    # --- Step 2: reverse closure / Hermiticity synthesis ---
    # For every declared R, ensure -R is present. If both are present,
    # validate consistency V_ab(R) == conj(V_ba(-R)); if only one is present,
    # synthesize the other.
    declared_irvecs = set(by_irvec.keys())
    all_irvecs = set(declared_irvecs)
    for irvec in declared_irvecs:
        neg = tuple(-x for x in irvec)
        all_irvecs.add(neg)

    completed = {}
    visited = set()
    for irvec in all_irvecs:
        if irvec in visited:
            continue
        neg = tuple(-x for x in irvec)
        visited.add(irvec)
        visited.add(neg)

        fwd = by_irvec.get(irvec, {})
        bwd = by_irvec.get(neg, {})

        fwd_mat = np.zeros((norb, norb), dtype=complex)
        bwd_mat = np.zeros((norb, norb), dtype=complex)

        # Build a synthesized version of each direction from the other:
        # V_ab(R) = conj(V_ba(-R))
        synth_fwd_from_bwd = {(a, b): np.conj(v) for (b, a), v in bwd.items()}
        synth_bwd_from_fwd = {(b, a): np.conj(v) for (a, b), v in fwd.items()}

        keys_fwd = set(fwd.keys()) | set(synth_fwd_from_bwd.keys())
        for (a, b) in keys_fwd:
            v_declared = fwd.get((a, b))
            v_synth = synth_fwd_from_bwd.get((a, b))
            if v_declared is not None and v_synth is not None:
                tol = reverse_atol + reverse_rtol * abs(v_declared)
                if abs(v_declared - v_synth) > tol:
                    raise ValueError(
                        "resolve_interactions: V_{}{}({}) = {} disagrees with "
                        "conj(V_{}{}({})) = {} beyond tolerance "
                        "(atol={}, rtol={})".format(
                            a, b, irvec, v_declared,
                            b, a, neg, v_synth,
                            reverse_atol, reverse_rtol,
                        )
                    )
                # Accepted within tolerance -> project onto the exactly
                # conjugate pair (review fix R2-2); see the same step in
                # _complete_onsite_hermitian. The (a,b)-of-R and (b,a)-of-(-R)
                # averages below are exact conjugates of one another, so the
                # returned set satisfies v_bond[m] == v_bond[reverse[m]].conj().T
                # to the last bit rather than only within tolerance.
                fwd_mat[a, b] = 0.5 * (v_declared + v_synth)
            elif v_declared is not None:
                fwd_mat[a, b] = v_declared
            else:
                fwd_mat[a, b] = v_synth
                if abs(v_synth) <= tol_incl:
                    dropped.append(
                        ("synthesized_negligible", irvec, (a, b), v_synth)
                    )

        keys_bwd = set(bwd.keys()) | set(synth_bwd_from_fwd.keys())
        for (a, b) in keys_bwd:
            v_declared = bwd.get((a, b))
            v_synth = synth_bwd_from_fwd.get((a, b))
            if v_declared is not None and v_synth is not None:
                tol = reverse_atol + reverse_rtol * abs(v_declared)
                if abs(v_declared - v_synth) > tol:
                    raise ValueError(
                        "resolve_interactions: V_{}{}({}) = {} disagrees with "
                        "conj(V_{}{}({})) = {} beyond tolerance "
                        "(atol={}, rtol={})".format(
                            a, b, neg, v_declared,
                            b, a, irvec, v_synth,
                            reverse_atol, reverse_rtol,
                        )
                    )
                bwd_mat[a, b] = 0.5 * (v_declared + v_synth)
            elif v_declared is not None:
                bwd_mat[a, b] = v_declared
            else:
                bwd_mat[a, b] = v_synth
                if abs(v_synth) <= tol_incl:
                    dropped.append(
                        ("synthesized_negligible", neg, (a, b), v_synth)
                    )

        completed[irvec] = fwd_mat
        if neg != irvec:
            completed[neg] = bwd_mat
        # if irvec == neg (only possible for (0,0,0), already excluded above)
        # no special handling needed.

    # --- Step 3: sort non-zero-shell channels by (length, lexicographic irvec) ---
    def sort_key(irvec):
        return (_shell_length(irvec, lattice_vectors), irvec)

    nonzero_irvecs = sorted(completed.keys(), key=sort_key)

    # Group into shells: irvecs whose bond length agrees within relative 1e-6.
    shells = []  # list of (length, [irvecs...])
    for irvec in nonzero_irvecs:
        length = _shell_length(irvec, lattice_vectors)
        if shells and abs(length - shells[-1][0]) <= 1e-6 * max(length, shells[-1][0], 1.0):
            shells[-1][1].append(irvec)
        else:
            shells.append((length, [irvec]))

    # --- Step 4: apply bond_max_shells cutoff ---
    # Shell 0 = Delta r = (0,0,0); shell n (n>=1) = nth nearest distinct
    # nonzero-shell length.
    # The ambiguity guard keys on values that are ACTUALLY DECLARED nonzero,
    # not on the tol_incl inclusion threshold (review fix R2-3): tol_incl only
    # governs provenance bookkeeping for SYNTHESIZED partners, so using it here
    # silently discarded a declared-but-tiny inter-site V that the documented
    # rule ("any declared nonzero inter-site V") says is ambiguous.
    has_nonzero_inter_site_v = any(
        value != 0 for entries in by_irvec.values() for value in entries.values()
    )
    if bond_max_shells is not None:
        # Validate BEFORE use (review fix R2-3): a negative value used to be
        # clamped to zero -- silently dropping every inter-site bond AND
        # bypassing the ambiguity guard below -- and a non-integral one used to
        # be truncated by int() without a word.
        try:
            shells_f = float(bond_max_shells)
        except (TypeError, ValueError):
            raise ValueError(
                "resolve_interactions: bond_max_shells must be a non-negative "
                "integer or None, got {!r}".format(bond_max_shells))
        if (not np.isfinite(shells_f) or shells_f < 0
                or shells_f != np.floor(shells_f)):
            raise ValueError(
                "resolve_interactions: bond_max_shells must be a non-negative "
                "integral value (shell 0 = the on-site Delta r = 0 point) or "
                "None, got {!r}".format(bond_max_shells))
        if shells_f == 0 and has_nonzero_inter_site_v:
            raise ValueError(
                "resolve_interactions: bond_max_shells=0 requested but nonzero "
                "inter-site V is declared; this is ambiguous (asks for B=1 "
                "while inter-site V is present). Use bond_channels=false for "
                "the genuinely local model instead."
            )
        # shells 1..bond_max_shells are kept (shell 0 is Delta r=0, always kept
        # separately below); bond_max_shells counts inter-site shells beyond
        # the on-site point, matching "shell 1 = 1st NN".
        n_keep = int(shells_f)
        kept_shells = shells[:n_keep]
        dropped_shells = shells[n_keep:]
        for length, irvecs in dropped_shells:
            for irvec in irvecs:
                dropped.append(("shell_cutoff", irvec, length))
        shells = kept_shells

    ordered_irvecs = [irvec for _, irvecs in shells for irvec in irvecs]

    # --- Step 5: build delta_r / v_bond / reverse ---
    delta_r = [(0, 0, 0)] + ordered_irvecs
    index_of = {dr: i for i, dr in enumerate(delta_r)}

    v_bond_list = [np.zeros((norb, norb), dtype=complex)]
    for irvec in ordered_irvecs:
        arr = np.array(completed[irvec], dtype=complex, copy=True)
        v_bond_list.append(arr)

    reverse = []
    for dr in delta_r:
        neg = tuple(-x for x in dr)
        if neg not in index_of:
            raise ValueError(
                "resolve_interactions: internal error -- channel {} has no "
                "reverse partner {} in the resolved set (reversal closure "
                "invariant violated)".format(dr, neg)
            )
        reverse.append(index_of[neg])

    for arr in v_bond_list:
        arr.flags.writeable = False

    v_onsite = _complete_onsite_hermitian(onsite, norb, reverse_atol, reverse_rtol)
    v_onsite.flags.writeable = False

    return ResolvedInteractionSet(
        delta_r=tuple(delta_r),
        v_bond=tuple(v_bond_list),
        reverse=tuple(reverse),
        n_channels=len(delta_r),
        v_onsite=v_onsite,
        dropped=tuple(dropped),
    )


# ---------------------------------------------------------------------------
# bare_bond_vertices's memory contract  --  read together with
# sc._bond_memory_estimate
# ---------------------------------------------------------------------------
#
# ``bare_bond_vertices`` allocates ``ND x ND`` (NOT q-resolved -- ``S_bond``/
# ``C_bond`` are the q-resolved arrays and are already counted among
# ``sc._BOND_N_Q_ARRAYS``) working buffers while assembling the Cooper
# vertices ``Vpp_s``/``Vpp_t`` (spec S4.5). None of them are ``del``-ed or
# reused today, so all ``BARE_VERTEX_ND2_BUFFERS`` of them are alive
# simultaneously at the function's high-water mark (its very end, just before
# it returns):
#
#   1. ``P``     -- the reversal+orbital-swap permutation matrix
#   2. ``D``     -- diagonal bond-channel interaction (m != 0 only)
#   3. ``Dh``    -- ``D + D^dag``
#   4. ``Id``    -- the ``ND x ND`` identity
#   5. ``Q_s``   -- ``(Id + P)/2`` singlet projector
#   6. ``Q_t``   -- ``(Id - P)/2`` triplet projector
#   7. ``B_s``   -- ``Q_s Dh Q_s``, the bond Cooper block (singlet)
#   8. ``B_t``   -- ``Q_t Dh Q_t``, the bond Cooper block (triplet)
#   9. ``Vpp_s`` -- the returned singlet Cooper vertex (local L_s + B_s)
#   10. ``Vpp_t``-- the returned triplet Cooper vertex (local L_t + B_t)
#
# (The local block temporaries -- ``S0_loc``, ``C0_loc``, ``L_s``, ``L_t``
# and ``_crossing``'s reshape views -- are ``nd x nd = norb**4``-sized, i.e.
# a factor ``B**2`` smaller than ``ND x ND`` for ``B > 1``, so they are not
# separately budgeted.) Imported by ``sc._bond_memory_estimate`` so the two
# sides cannot silently desync; a new persistent ``ND x ND`` buffer here must
# bump this count (and ``tests/test_sc_bond.py`` measures the real peak
# against the budget).
BARE_VERTEX_ND2_BUFFERS = 10


def bare_bond_vertices(bond_set, S0_q, C0_q, norb):
    """Build the bare enlarged vertices ``S_bond``, ``C_bond`` and the bare
    particle-particle (Cooper) vertices ``Vpp_s``, ``Vpp_t``.

    Implements spec S4.3 (fully-indexed S_bond/C_bond element equations) and
    S4.5 (the complete Cooper vertex ``V^pp_eta = Q_eta L_eta Q_eta +
    Q_eta(D+D^dag)Q_eta``). See
    ``docs/superpowers/specs/2026-07-25-bond-channels-eliashberg-design.md``.

    Enlarged bond-major index ``I = m*nd + idx`` with ``nd = norb**2`` and
    ``idx = l1*norb + l2``; ``ND = nd*B`` where ``B = bond_set.n_channels``.

    Parameters
    ----------
    bond_set : ResolvedInteractionSet
        Reversal-closed bond-channel topology (``delta_r``, ``v_bond``,
        ``reverse``, ``n_channels``) from :func:`resolve_interactions`.
    S0_q, C0_q : ndarray, shape (Nx, Ny, Nz, nd, nd)
        The ``m=0`` ph blocks -- the existing Kuroki spin/charge matrices
        (intra ``U``, inter-orbital ``U'``, Hund, Exchange, PairHop) PLUS the
        inter-site Hartree ``2 V_ab(q)`` in the charge ``(aa,bb)`` element.
        Assumed already **Case-2-corrected**: only the on-site ``V_ab(R=0)``
        component stays in the ``m=0`` Fock ``(ab,ab)`` element; the ``R!=0``
        Fock lives in ``bond_set.v_bond`` (spec S4.3 star).
    norb : int
        Number of orbitals.

    Returns
    -------
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        The ``m=0`` sub-block equals ``S0_q``/``C0_q``; the bond-diagonal
        ``m!=0`` Fock element ``(I,I)`` with ``I=m*nd+(l1*norb+l2)`` equals
        ``+Re(V_{l1 l2}(Delta r_m))`` (S) / ``-Re(V_{l1 l2}(Delta r_m))`` (C)
        (spec S4.3, Hermitized per Finding 1 of the #151 ED-adjudication
        campaign -- see the inline comment at the assignment site below for
        the derivation; bit-identical to the un-Hermitized value whenever
        ``V_{l1 l2}(Delta r_m)`` is already real). q-dependent (the Hartree
        lives in the ``m=0`` block).
    Vpp_s, Vpp_t : ndarray, shape (ND, ND), complex
        q-**independent** bare Cooper vertices. The local ``m=0`` block is the
        pair matrix ``L^{s/t}`` obtained from ``(S0 +- C0_local)`` under the
        ``(ac,db)->(ab,cd)`` pair<->ph crossing, with the **inter-site**
        (``R!=0``) Hartree excluded so ``V^pp`` carries only local
        interactions. The ``m!=0`` block is ``Q_eta (D_off + D_off^dag) Q_eta``
        (``= D_off (1 +- P)`` for real inversion-symmetric bonds), with
        ``D_off`` diagonal ``V_ab(Delta r_m)`` (``m!=0`` only) and
        ``P:(m,ab)->(m_bar,ba)`` the reversal+orbital-swap involution. Since
        ``D_off`` is diagonal, ``D_off + D_off^dag`` already reduces
        elementwise to ``2*Re(D_off)`` before the ``Q_eta`` projection below
        -- the SAME Hermitization the ph sector's ``m!=0`` diagonal applies
        explicitly above -- so this pp block is Hermitization-invariant BY
        CONSTRUCTION and is intentionally left untouched by the Finding 1 fix
        (do not also apply ``.real`` to ``D_off`` -- it is already redundant
        and would not change ``Dh``'s value).
    """
    norb = int(norb)
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    delta_r = bond_set.delta_r
    v_bond = bond_set.v_bond
    reverse = bond_set.reverse

    S0_q = np.asarray(S0_q)
    C0_q = np.asarray(C0_q)
    if S0_q.ndim != 5 or S0_q.shape[3:] != (nd, nd):
        raise ValueError(
            "bare_bond_vertices: S0_q must have shape (Nx, Ny, Nz, nd, nd) "
            "with nd=norb**2={}; got {}".format(nd, S0_q.shape))
    if C0_q.shape != S0_q.shape:
        raise ValueError(
            "bare_bond_vertices: C0_q shape {} must match S0_q shape {}".format(
                C0_q.shape, S0_q.shape))
    Nx, Ny, Nz = S0_q.shape[:3]

    # --- internal-interface guard: bond_set MUST be reversal-closed ---
    # (a hand-built non-closed set is a programming error; fail fast rather
    # than silently produce a non-Hermitian Cooper vertex; spec test 10a).
    for m in range(B):
        mr = reverse[m]
        if not (0 <= mr < B) or reverse[mr] != m:
            raise ValueError(
                "bare_bond_vertices: bond_set.reverse is not an involution at "
                "channel {} (reverse={}); the set is not reversal-closed".format(
                    m, reverse))
        neg = tuple(-x for x in delta_r[m])
        if tuple(delta_r[mr]) != neg:
            raise ValueError(
                "bare_bond_vertices: bond_set is not reversal-closed -- "
                "delta_r[reverse[{0}]]={1} != -delta_r[{0}]={2}".format(
                    m, tuple(delta_r[mr]), neg))

    # === S_bond / C_bond (spec S4.3) =======================================
    S_bond = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
    C_bond = np.zeros((Nx, Ny, Nz, ND, ND), dtype=complex)
    # m=0 sub-block = existing Kuroki matrices (Hartree included)
    S_bond[:, :, :, 0:nd, 0:nd] = S0_q
    C_bond[:, :, :, 0:nd, 0:nd] = C0_q
    # m!=0 bond-diagonal Fock: +Re(V_ab(R_m)) (S) / -Re(V_ab(R_m)) (C)
    #
    # Hermitization derivation (Finding 1, #151 ED-adjudication campaign;
    # tests/test_bond_vs_ed_oracle.py::TestNullDirectionSolverSide). A declared
    # channel V_ab(R_m) = v_bond[m][a,b] and its resolve_interactions reversal
    # partner V_ba(-R_m) = v_bond[reverse[m]][b,a] multiply the SAME physical
    # operator sum_j n_{j,a} n_{j+R_m,b}: substituting j' = j+R_m in the
    # partner's sum sum_j n_{j,b} n_{j-R_m,a} reproduces the ORIGINAL sum
    # verbatim (density operators commute), so the two declarations are not
    # independent physical terms but two Fourier-index views of one term whose
    # Hamiltonian coefficient is their SUM, V_ab(R_m) + V_ba(-R_m).
    # resolve_interactions guarantees, EXACTLY (to the last bit, by
    # construction -- see its docstring's "returned set satisfies
    # v_bond[m] == v_bond[reverse[m]].conj().T"), that
    # v_bond[reverse[m]][b,a] == conj(v_bond[m][a,b]); substituting,
    # V_ab(R_m) + V_ba(-R_m) == v_bond[m][a,b] + conj(v_bond[m][a,b])
    # == 2*Re(v_bond[m][a,b]). This is a Hermiticity consequence of the
    # operator algebra, not a special case for a "null" pair: it holds
    # identically whether the closed pair is a synthetic +/-i*eps probe or a
    # genuinely complex, phase-carrying single-sided declaration (whose
    # reverse partner resolve_interactions synthesizes as the exact
    # conjugate) -- a complex CoulombInter value ALWAYS collapses to twice its
    # real part in this density-density Hamiltonian term.
    #
    # This enlarged-index slot I=(m,l1,l2) is never summed with its mirror
    # slot (reverse[m],l2,l1) anywhere downstream (S_bond/C_bond only ever set
    # the (I,I) diagonal, independently per slot) -- unlike the pp sector's
    # D+D^dagger construction below, which projects Q_eta (I +/- P) Q_eta
    # across BOTH mirror slots together and is untouched here. So each ph
    # slot must carry HALF of the physical total, Re(v_bond[m][a,b]), not the
    # full 2*Re(v_bond[m][a,b]) -- using the full sum here would double every
    # already-validated real-coupling S_bond/C_bond entry pinned by this
    # campaign's granules and the #82 golden tests.
    #
    # Before this fix, S_bond/C_bond took vb[l1, l2] RAW (no Hermitization):
    # V_01(+1)=V+i*eps and V_10(-1)=V-i*eps live at DIFFERENT slots and were
    # never combined the way the ED oracle's canonical_density_terms combines
    # them, so S_bond leaked the null direction at magnitude eps. C_bond
    # picked up this SAME Fock effect (HALF its total response) PLUS a
    # separate eps response of equal magnitude from the sc.py Hartree BUG,
    # fixed below (sc._build_bond_m0_blocks summed the raw, un-Hermitized
    # per-shell value into V_ab(q) before this fix).
    #
    # Re(vb[l1, l2]) is bit-identical to vb[l1, l2] whenever the stored value
    # is already real (every direction adjudicated by this campaign, and
    # every #82 golden-test coupling), so this is a provable no-op there.
    for m in range(1, B):
        vb = v_bond[m]
        for l1 in range(norb):
            for l2 in range(norb):
                I = m * nd + l1 * norb + l2
                S_bond[:, :, :, I, I] = vb[l1, l2].real
                C_bond[:, :, :, I, I] = -vb[l1, l2].real

    # === Local Cooper block L^{s/t} (spec S4.5, crossing) ==================
    # The local interactions are q-independent; only the inter-site Hartree is
    # q-dependent. Evaluate at q0=(0,0,0), where every bond phase is 1, and
    # remove the R!=0 Hartree (2 * sum_{m>=1} V_ab(R_m)) from the (aa,bb)
    # charge elements. The on-site R=0 Hartree (= 2 U'_ab) stays: it is a
    # local interaction and yields L^s(ab,ab)=2U' (spec S4.5 table).
    S0_loc = S0_q[0, 0, 0].copy()
    C0_loc = C0_q[0, 0, 0].copy()
    for a in range(norb):
        for b in range(norb):
            hartree_rneq0 = 0.0 + 0.0j
            for m in range(1, B):
                # Same Hermitization as the m!=0 Fock diagonal above (see its
                # comment for the derivation), and for the same reason: this
                # subtraction must exactly cancel the R!=0 Hartree
                # sc._build_bond_m0_blocks adds into C0_q's (aa,bb) element at
                # q=0, which Hermitizes each shell's contribution the
                # identical way (per-shell .real BEFORE the q-phase multiply).
                # Leaving this raw while that Hartree builder is fixed would
                # reintroduce an eps-dependent residual into
                # C0_loc -> L_s/L_t -> Vpp_s/Vpp_t, breaking the pp sector's
                # currently-exact null invariance (Finding 1's report: "Vpp_s/
                # Vpp_t ... holds exactly").
                hartree_rneq0 += v_bond[m][a, b].real
            C0_loc[a * norb + a, b * norb + b] -= 2.0 * hartree_rneq0

    def _crossing(M):
        # L_{(ab),(cd)} = M_{(ac),(db)}: reshape to (l1,l2,l3,l4) and permute
        # so L4[a,b,c,d] = M4[a,c,d,b].
        M4 = M.reshape(norb, norb, norb, norb)
        L4 = M4.transpose(0, 3, 1, 2)
        return L4.reshape(nd, nd).copy()

    L_s = _crossing(S0_loc + C0_loc)   # singlet: (S+C)
    L_t = _crossing(C0_loc - S0_loc)   # triplet: (C-S)

    # === Bond Cooper block Q_eta (D_off + D_off^dag) Q_eta (spec S4.5) =====
    # P:(m,l1,l2) -> (reverse[m], l2, l1). As a permutation matrix,
    # P[I,J]=1 iff state I = P(state J).
    P = np.zeros((ND, ND), dtype=complex)
    for m in range(B):
        mr = reverse[m]
        for l1 in range(norb):
            for l2 in range(norb):
                J = m * nd + l1 * norb + l2
                I = mr * nd + l2 * norb + l1
                P[I, J] = 1.0
    # D_off: diagonal V_ab(R_m), m != 0 only.
    D = np.zeros((ND, ND), dtype=complex)
    for m in range(1, B):
        vb = v_bond[m]
        for l1 in range(norb):
            for l2 in range(norb):
                I = m * nd + l1 * norb + l2
                D[I, I] = vb[l1, l2]
    Dh = D + D.conj().T
    Id = np.eye(ND, dtype=complex)
    Q_s = 0.5 * (Id + P)
    Q_t = 0.5 * (Id - P)
    B_s = Q_s @ Dh @ Q_s
    B_t = Q_t @ Dh @ Q_t

    Vpp_s = np.zeros((ND, ND), dtype=complex)
    Vpp_t = np.zeros((ND, ND), dtype=complex)
    Vpp_s[0:nd, 0:nd] = L_s
    Vpp_t[0:nd, 0:nd] = L_t
    # The bond block is supported purely on the m!=0 sector (D vanishes on
    # m=0 and P does not couple m=0 to m!=0), so this add does not disturb the
    # local L block.
    Vpp_s += B_s
    Vpp_t += B_t

    return S_bond, C_bond, Vpp_s, Vpp_t


# Conditioning floor of the enlarged RPA denominators, applied by dress_bond to
# BOTH the relative ratio sigma_min/sigma_max and the absolute pole distance
# sigma_min/max(1, sigma_max) of the ACTUAL solve matrices, minimized over q
# (see _check_bond_conditioning). Same value -- and, for the relative half, the
# same criterion -- as the milestone off-instability guard,
# tests/test_bond_onari_milestone.py::SIGMA_FLOOR, which pins BOTH the ratio and
# sigma_min above this floor.
_BOND_COND_FLOOR = 1.0e-3


def _check_bond_conditioning(name, mat, cond_tol):
    """Refuse a singular / nearly singular enlarged RPA denominator.

    Returns the guard score (the minimum over q of the smaller of the two
    criteria below) when the block passes, ``None`` when ``cond_tol`` is
    ``None`` (guard disabled); raises ``ValueError`` otherwise (#181 Tier 3
    Phase A: the score is published as ``longitudinal_bond_cond_min_*``).

    ``np.linalg.solve`` raises a bare ``LinAlgError("Singular matrix")`` on an
    exactly singular block -- with no indication of WHICH channel or WHICH
    q-point diverged -- and returns enormous, unreliable numbers just short of
    that, silently. Both are the RPA instability of the bond path, so they get
    one actionable error naming the channel, the q-point and both conditioning
    numbers.

    TWO criteria are applied to each q-block of the ACTUAL solve matrix
    ``I -/+ chi_bar V``; the block is refused when EITHER falls to ``cond_tol``
    or below (the guard score is their minimum, minimized over q):

    1. the RELATIVE ``sigma_min/sigma_max`` ("ratio"), the scale-free
       conditioning number, identical to the milestone off-instability guard
       (``tests/test_bond_onari_milestone.py::SIGMA_FLOOR``);
    2. the ABSOLUTE distance to the pole measured on the natural ``O(1)`` scale
       of ``I -/+ chi_bar V``, ``sigma_min / max(1, sigma_max)`` ("pole
       distance"). Criterion 1 alone is blind whenever the singular values are
       UNIFORMLY small: it is identically 1 for any nonzero ``1x1`` block (the
       reachable ``B = 1, norb = 1`` case of an empty/local interaction set) and
       for any ``eps * I``, so a denominator of ``1e-12`` would pass while the
       dressed vertices came back of order ``1e12``. Normalizing by
       ``max(1, sigma_max)`` rather than by ``sigma_max`` keeps this a genuine
       absolute floor near the identity (where ``sigma_max ~ 1``) without
       double-penalizing a block whose largest singular value is huge --
       criterion 1 already covers that regime.
    """
    if cond_tol is None:
        return None
    Nx, Ny, Nz, ND, _ = mat.shape
    sv = np.linalg.svd(mat.reshape(-1, ND, ND), compute_uv=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sv[:, -1] / sv[:, 0]
    # A zero (or non-finite) largest singular value means the block is the zero
    # matrix -- as singular as it gets.
    ratio = np.where(np.isfinite(ratio), ratio, 0.0)
    # Absolute pole distance on the natural scale of I -/+ chi_bar V.
    pole = sv[:, -1] / np.maximum(1.0, sv[:, 0])
    pole = np.where(np.isfinite(pole), pole, 0.0)
    score = np.minimum(ratio, pole)
    iq = int(np.argmin(score))
    worst = float(score[iq])
    if worst > cond_tol:
        return worst
    qx, rem = divmod(iq, Ny * Nz)
    qy, qz = divmod(rem, Nz)
    raise ValueError(
        "dress_bond: the {} RPA denominator is singular or nearly singular at "
        "the q-point index ({}, {}, {}): sigma_min/sigma_max = {:.3e}, "
        "sigma_min/max(1, sigma_max) = {:.3e}; the smaller of the two is "
        "<= cond_tol = {:.3e} (sigma_min = {:.3e}, sigma_max = {:.3e}). The "
        "bond path has entered the {} instability region, where the dressed "
        "vertices are enormous and numerically meaningless. Reduce the "
        "interaction strength, raise the temperature, refine/reduce the "
        "q-grid, or -- if you deliberately want to study the stiff regime -- "
        "lower cond_tol.".format(
            name, qx, qy, qz, float(ratio[iq]), float(pole[iq]), cond_tol,
            sv[iq, -1], sv[iq, 0], name))


def dress_bond(chi_bar, S_bond, C_bond, cond_tol=_BOND_COND_FLOOR):
    """Dress the enlarged bond bubble into the RPA spin/charge susceptibilities
    at the enlarged bond-major index (spec S4.4).

    Batched linear solves (not explicit inverses), one system per q-point,
    mirroring the exact orientation of the existing ``sc.py`` dressing
    (``_compute_vertices_general``, ``sc.py:1098-1104``) generalized from the
    plain orbital-pair size ``nd`` to the enlarged bond-major size ``ND``::

        chi_s = solve(I_ND - chi_bar @ S_bond, chi_bar)
        chi_c = solve(I_ND + chi_bar @ C_bond, chi_bar)

    ``S_bond``/``C_bond`` are block-diagonal in the bond index ``m`` (see
    :func:`bare_bond_vertices`), but ``chi_bar`` (from :func:`bond_bubble`) is
    not, so the solve mixes bond sectors -- ``chi_s``/``chi_c`` are full
    ``ND x ND`` matrices, not block-diagonal.

    Parameters
    ----------
    chi_bar : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        The enlarged bond bubble (:func:`bond_bubble`).
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        The enlarged bare vertices (:func:`bare_bond_vertices`).
    cond_tol : float or None, optional
        Conditioning floor of the two solve matrices, applied per q-point to
        BOTH the relative ratio ``sigma_min/sigma_max`` and the absolute pole
        distance ``sigma_min/max(1, sigma_max)`` (default ``1e-3``, the
        milestone off-instability floor). The relative ratio alone cannot see
        a uniformly small denominator -- it is identically 1 for a ``1x1``
        block (``B = 1, norb = 1``) and for any ``eps * I`` -- so the absolute
        criterion, measured on the natural ``O(1)`` scale of
        ``I -/+ chi_bar V``, is applied alongside it; see
        :func:`_check_bond_conditioning`. A denominator at or below the floor
        on EITHER criterion raises with the channel, the q-point and both
        measured quantities named -- an exactly singular block would otherwise
        raise a bare ``LinAlgError`` and a nearly singular one would silently
        return enormous, unreliable vertices. ``None`` disables the check.
        This is a BOND-PATH guard only; the legacy ``sc.py`` dressing is
        untouched.

    Returns
    -------
    chi_s, chi_c : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        The dressed spin (minus) and charge (plus) susceptibilities at the
        enlarged bond-major index.

    Raises
    ------
    ValueError
        If either RPA denominator fails the ``cond_tol`` conditioning check.
    """
    chi_bar = np.asarray(chi_bar)
    if chi_bar.ndim != 5 or chi_bar.shape[3] != chi_bar.shape[4]:
        raise ValueError(
            "dress_bond: chi_bar must have shape (Nx, Ny, Nz, ND, ND); got {}"
            .format(chi_bar.shape)
        )
    Nx, Ny, Nz, ND, _ = chi_bar.shape

    S_bond = np.asarray(S_bond)
    C_bond = np.asarray(C_bond)
    if S_bond.shape != chi_bar.shape:
        raise ValueError(
            "dress_bond: S_bond shape {} must match chi_bar shape {}".format(
                S_bond.shape, chi_bar.shape))
    if C_bond.shape != chi_bar.shape:
        raise ValueError(
            "dress_bond: C_bond shape {} must match chi_bar shape {}".format(
                C_bond.shape, chi_bar.shape))

    # Batched RPA solve for all q-points simultaneously (same algebra/
    # orientation as sc.py._compute_vertices_general, generalized nd -> ND).
    # Kept as the legacy schedule -- BOTH denominators are formed and
    # conditioning-checked before either solve -- so the Eliashberg path
    # is unchanged in numbers, error text and order of diagnostics. The
    # production longitudinal bond gate uses the one-channel-at-a-time
    # ``dress_channel`` below instead (#181 Tier 3 Phase A); the two agree
    # bit-for-bit (tests/test_bond_longitudinal_vertex.py).
    I_mat = np.broadcast_to(np.eye(ND, dtype=complex), (Nx, Ny, Nz, ND, ND)).copy()

    mat_s = I_mat - chi_bar @ S_bond
    mat_c = I_mat + chi_bar @ C_bond

    # Off-instability guard BEFORE the solve, so an unstable point is named
    # rather than crashing (exactly singular) or silently producing garbage
    # (nearly singular).
    _check_bond_conditioning("spin", mat_s, cond_tol)
    _check_bond_conditioning("charge", mat_c, cond_tol)

    chi_s = np.linalg.solve(mat_s, chi_bar)
    chi_c = np.linalg.solve(mat_c, chi_bar)

    return chi_s, chi_c


_DRESS_CHANNELS = {"spin": -1.0, "charge": +1.0}


def dress_channel(chi_bar, W, channel, *, spatial_shape, cond_tol=_BOND_COND_FLOOR):
    """Dress ONE RPA channel on the (bond-enlarged) pair basis::

        chi = [I + sign * chi_bar @ W]^{-1} chi_bar,   sign = -1 ("spin"),
                                                              +1 ("charge")

    ``chi_bar``/``W`` are ``(nvol, ND, ND)``; ``spatial_shape`` =
    ``(Nx, Ny, Nz)`` with ``prod == nvol`` (the q-point named in a
    conditioning refusal is decoded from it). ONE denominator is formed,
    conditioning-checked (:func:`_check_bond_conditioning`, criteria
    ``sigma_min/sigma_max`` and ``sigma_min/max(1, sigma_max)``; refused
    at ``<= cond_tol``, disabled by ``cond_tol=None``) and solved with a
    batched ``numpy.linalg.solve`` (no explicit inverse). Returns
    ``(chi, cond_min)``, ``cond_min`` the guard score (``None`` when the
    guard is disabled). Peak transient: the denominator plus the solve's
    LU workspace; the denominator is released before returning.
    """
    if channel not in _DRESS_CHANNELS:
        raise ValueError(
            "dress_channel: channel must be 'spin' or 'charge', got {!r}"
            .format(channel))
    sign = _DRESS_CHANNELS[channel]
    chi_bar = np.asarray(chi_bar)
    W = np.asarray(W)
    if (chi_bar.ndim != 3 or chi_bar.shape[1] != chi_bar.shape[2]
            or W.shape != chi_bar.shape):
        raise ValueError(
            "dress_channel: chi_bar and W must both have shape (nvol, ND, "
            "ND); got {} and {}".format(chi_bar.shape, W.shape))
    try:
        shape = tuple(spatial_shape)
        ok = (len(shape) == 3 and all(
            not isinstance(x, bool) and isinstance(x, (int, np.integer))
            and int(x) > 0 for x in shape))
    except TypeError:
        ok = False
    if not ok:
        raise ValueError(
            "dress_channel: spatial_shape must be three positive integers "
            "(Nx, Ny, Nz), got {!r}".format(spatial_shape))
    Nx, Ny, Nz = (int(x) for x in shape)
    nvol, ND = chi_bar.shape[0], chi_bar.shape[1]
    if Nx * Ny * Nz != nvol:
        raise ValueError(
            "dress_channel: prod(spatial_shape) = {} != nvol = {}".format(
                Nx * Ny * Nz, nvol))
    mat = chi_bar @ W
    if sign < 0:
        np.negative(mat, out=mat)
    idx = np.arange(ND)
    mat[:, idx, idx] += 1.0
    cond_min = _check_bond_conditioning(
        channel, mat.reshape(Nx, Ny, Nz, ND, ND), cond_tol)
    chi = np.linalg.solve(mat, chi_bar)
    del mat
    return chi, cond_min


def _validate_green_beta(green, beta, label):
    """Validate the documented ``green``/``beta`` inputs of the kernel path.

    ``green`` must have the 6-D ``sc.py`` layout ``(norb, norb, Nx, Ny, Nz,
    nmat)`` with equal orbital axes, and ``beta`` must be positive and finite.
    :func:`bond_bubble` already checked this; :func:`make_bond_kernel`,
    :func:`make_bond_kernel_parts` and :func:`pair_weight` used to index
    ``green.shape`` blind, turning a mis-shaped input into an ``IndexError``,
    an opaque einsum/reshape failure, or (for ``beta = 0``) a silent array of
    ``inf``/``nan`` (review fix R2-4).

    Returns the array form of ``green``.
    """
    green = np.asarray(green)
    if green.ndim != 6:
        raise ValueError(
            "{}: green must have shape (norb, norb, Nx, Ny, Nz, nmat); got "
            "ndim={} shape={}".format(label, green.ndim, green.shape))
    if green.shape[0] != green.shape[1]:
        raise ValueError(
            "{}: green's two orbital axes must both equal norb; got shape "
            "{}".format(label, green.shape))
    try:
        beta_f = float(beta)
    except (TypeError, ValueError):
        raise ValueError(
            "{}: beta must be a positive finite number, got {!r}".format(
                label, beta))
    if not np.isfinite(beta_f) or beta_f <= 0.0:
        raise ValueError(
            "{}: beta must be a positive finite number, got {!r}".format(
                label, beta))
    return green


def _g2_from_green(green, beta, tail=False):
    """Static two-particle bubble ``G2[i,j,l,m](k) = (T) sum_n
    G_{ij}(k, iw_n) G_{lm}(-k, -iw_n)`` (``T = 1/beta``).

    This is the bond-kernel implementation (roll+flip for ``G(-k, -w_n)``,
    per-site batched GEMM over the Matsubara axis). With ``tail=True`` it adds
    the same analytic ``1/wn**2`` window-complement correction as
    :func:`hwave.sc._calc_g2`; the default remains the historical raw-window
    behavior for direct callers, while the top-level solver forwards its
    ``[eliashberg] g2_tail`` setting explicitly.
    """
    norb = green.shape[0]
    Nx, Ny, Nz, nmat = green.shape[2], green.shape[3], green.shape[4], green.shape[5]
    nvol = Nx * Ny * Nz
    green_inv = np.roll(
        green[:, :, ::-1, ::-1, ::-1, ::-1], (1, 1, 1), (2, 3, 4))
    A = green.reshape(norb * norb, nvol, nmat)
    Bm = green_inv.reshape(norb * norb, nvol, nmat)
    As = np.moveaxis(A, 1, 0)
    Bs = np.moveaxis(Bm, 1, 0)
    G2 = np.moveaxis(As @ Bs.transpose(0, 2, 1), 0, 2)  # (norb^2, norb^2, nvol)
    G2 = G2.reshape(norb, norb, norb, norb, Nx, Ny, Nz)
    G2 = G2 / beta
    if tail:
        if nmat <= 0 or nmat % 2 != 0:
            raise ValueError(
                "the Matsubara tail correction (g2_tail) requires an even, "
                "positive number of frequencies on the centered fermionic "
                "grid; the Green function has nmat = {}".format(nmat))
        r = 2.0 * np.arange(nmat) + 1.0 - nmat
        coeff = beta * (0.25 - np.sum(1.0 / r**2) / np.pi**2)
        di = np.arange(norb)
        G2[di[:, None], di[:, None], di[None, :], di[None, :]] += coeff
    return G2


def make_bond_kernel(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                     green, bond_set, pairing_type, beta, part="full",
                     g2_tail=False):
    """Build the two-part bond-resolved linearized Eliashberg operator (spec
    S4.5) as a scipy :class:`~scipy.sparse.linalg.LinearOperator`.

    The operator acts on a flattened pairing gap ``phi`` of shape
    ``(norb, norb, Nx, Ny, Nz)`` (the bond channel index is **internal** -- it
    is summed over and never appears in ``phi``).  It has two pieces with
    **different** phase structure (spec S4.5, "Eliashberg operator -- two
    separable-phase parts"):

    Fluctuation part ``K^fl`` (particle-hole ladder, both legs carry ``+``
    phases ``e^{i(k.dr_m + k'.dr_m')}``, convolution in ``q = k - k'``)::

        Y^fl_{l1l2}(k) = -(1/N) sum_{m,m'} e^{i k.dr_m}
            [ F_eta[m,m'] *_q ( e^{i k'.dr_m'} X_{l3l4}(k') ) ](k)

    with the fluctuation vertex (matrix products in the enlarged ``ND`` space at
    each ``q``)::

        F_s(q) = +1.5 S_bond chi_s S_bond - 0.5 C_bond chi_c C_bond
        F_t(q) = -0.5 S_bond chi_s S_bond - 0.5 C_bond chi_c C_bond

    Realized by the **O(B) FFT factorization** (spec S4.5): form ``B`` input
    fields ``e^{i k'.dr_m'} X`` -> IFFT -> real-space MAC with ``F_hat`` -> FFT
    -> ``x e^{i k.dr_m}``; ``2B`` spatial transforms, ``B^2`` real-space
    multiply-accumulates, preserving ``O(N_k log N_k)``.

    Instantaneous / bare particle-particle part ``K^pp`` (Cooper, low-rank, NO
    convolution; cross phase ``e^{i(k.dr_m' - k'.dr_m)}``)::

        A_{m; l3l4}    = (1/N) sum_{k'} e^{-i k'.dr_m} X_{l3l4}(k')
        Y^pp_{l1l2}(k) = -1/2 sum_{m,m'} [V^pp_eta]_{(m,l1l2),(m',l3l4)}
                              e^{i k.dr_m'} A_{m; l3l4}

    In both parts the ``GG phi`` pair bubble is
    ``X_{l3l4}(k') = sum_{l5,l6} G2_{l3,l5,l4,l6}(k') phi_{l5l6}(k')`` with
    ``G2 = (T) sum_n G(k',w_n) G(-k',-w_n)`` (:func:`_g2_from_green`), carrying
    the single factor of ``T = 1/beta``; the single factor of ``1/N`` comes from
    the spatial fold (the FFT convolution for ``K^fl``, the explicit reduction
    for ``K^pp``).  This mirrors the ``-(T/N)`` normalization of the existing
    static solver ``sc._make_kernel_operator`` (which applies the full ``G2``,
    i.e. the ``sqrt(GG)``-symmetrizable pair weight, without an extra ``T`` or
    ``1/N``).

    Parameters
    ----------
    chi_s, chi_c : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        Dressed spin/charge susceptibilities at the enlarged bond-major index
        (:func:`dress_bond`).  ``ND = norb**2 * B``.
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        Bare enlarged spin/charge vertices (:func:`bare_bond_vertices`).
    Vpp_s, Vpp_t : ndarray, shape (ND, ND), complex
        q-independent bare Cooper vertices (:func:`bare_bond_vertices`).
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat), complex
        k-space Matsubara Green function in the ``sc.py`` layout.
    bond_set : ResolvedInteractionSet
        Bond-channel topology (``delta_r``, ``n_channels``).
    pairing_type : {"singlet", "triplet"}
        Selects ``F_s``/``Vpp_s`` (singlet) or ``F_t``/``Vpp_t`` (triplet).
    beta : float
        Inverse temperature (the ``T = 1/beta`` factor in ``G2``).
    g2_tail : bool, optional
        Apply the analytic high-frequency Matsubara correction to ``G2``.
        Direct library calls default to the historical raw-window sum; the
        top-level solver passes its configured value explicitly.

    part : {"full", "fluctuation", "instantaneous"}, optional
        Which piece of ``K = K^fl + K^pp`` the returned operator applies.
        ``"full"`` (the default) is the physical kernel; the two pieces exist
        for the ``lambda = lambda^pp + lambda^fl`` attribution of spec S4.5
        (see :func:`make_bond_kernel_parts` and :func:`attribute_lambda`).

    Returns
    -------
    A : scipy.sparse.linalg.LinearOperator
        The Eliashberg kernel ``K``; ``A.matvec(phi.ravel())`` returns
        ``(K phi).ravel()`` with ``phi`` shaped ``(norb, norb, Nx, Ny, Nz)``.
    vec_size : int
        ``norb**2 * Nx * Ny * Nz``.
    """
    ops, vec_size = _bond_kernel_operators(
        chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set,
        pairing_type, beta, g2_tail=g2_tail)
    if part not in ops:
        raise ValueError(
            "make_bond_kernel: unknown part '{}'. Use one of {}.".format(
                part, sorted(ops)))
    return ops[part], vec_size


def make_bond_kernel_parts(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                           green, bond_set, pairing_type, beta, *,
                           g2_tail=False):
    """Build the bond Eliashberg operator together with its two parts.

    Same arguments as :func:`make_bond_kernel`. The three operators share one
    set of precomputed invariants (the pair bubble ``G2``, the real-space
    fluctuation vertex, the bond phases), so this is cheaper than three
    separate :func:`make_bond_kernel` calls and guarantees the identity
    ``K = K^fl + K^pp`` holds operator-wise (spec S4.5).

    Returns
    -------
    A_full, A_fl, A_pp : scipy.sparse.linalg.LinearOperator
        The full kernel, its fluctuation (particle-hole ladder ``F_eta``) part
        and its instantaneous (bare particle-particle ``1/2 V^pp``) part.
    vec_size : int
    """
    ops, vec_size = _bond_kernel_operators(
        chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t, green, bond_set,
        pairing_type, beta, g2_tail=g2_tail)
    return ops["full"], ops["fluctuation"], ops["instantaneous"], vec_size


def _bond_kernel_operators(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                           green, bond_set, pairing_type, beta, *,
                           g2_tail=False):
    """Shared builder behind :func:`make_bond_kernel` /
    :func:`make_bond_kernel_parts`; returns ``({part: LinearOperator},
    vec_size)`` with parts ``"full"``, ``"fluctuation"``, ``"instantaneous"``.
    """
    if pairing_type == "singlet":
        F_q = 1.5 * (S_bond @ chi_s @ S_bond) - 0.5 * (C_bond @ chi_c @ C_bond)
        Vpp = np.asarray(Vpp_s)
    elif pairing_type == "triplet":
        F_q = -0.5 * (S_bond @ chi_s @ S_bond) - 0.5 * (C_bond @ chi_c @ C_bond)
        Vpp = np.asarray(Vpp_t)
    else:
        raise ValueError(
            "make_bond_kernel: unknown pairing_type '{}'. Use 'singlet' or "
            "'triplet'.".format(pairing_type))

    green = _validate_green_beta(green, beta, "make_bond_kernel")
    norb = green.shape[0]
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nx, Ny, Nz, nmat = green.shape[2], green.shape[3], green.shape[4], green.shape[5]
    N = Nx * Ny * Nz
    delta_r = bond_set.delta_r

    F_q = np.asarray(F_q)
    if F_q.shape != (Nx, Ny, Nz, ND, ND):
        raise ValueError(
            "make_bond_kernel: fluctuation vertex has shape {} but expected "
            "(Nx, Ny, Nz, ND, ND) = {} (ND = norb**2 * n_channels).".format(
                F_q.shape, (Nx, Ny, Nz, ND, ND)))
    if Vpp.shape != (ND, ND):
        raise ValueError(
            "make_bond_kernel: Vpp has shape {} but expected (ND, ND) = "
            "{}.".format(Vpp.shape, (ND, ND)))

    vec_size = nd * N

    # --- invariants (hoisted out of the per-matvec closure) ---------------
    # X_{cd}(k) = sum_{ef} G2p[(cd),(ef),k] phi_{ef}(k), cd = c*norb+d.
    G2 = _g2_from_green(green, beta, tail=g2_tail)       # (i,j,l,m,Nx,Ny,Nz)
    G2p = G2.transpose(0, 2, 1, 3, 4, 5, 6).reshape(nd, nd, N)  # [(c,d),(e,f),k]

    # Real-space fluctuation vertex F_hat[r, m, ab, m', cd] = IFFT_q F_q.
    F_r6 = _bk.spatial_ifftn(F_q, axes=(0, 1, 2)).reshape(
        Nx, Ny, Nz, B, nd, B, nd)

    # Bond phases PH[m](k) = e^{i k.dr_m} on the k-grid.
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    kz = 2.0 * np.pi * np.arange(Nz) / Nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    PH = np.empty((B, Nx, Ny, Nz), dtype=complex)
    for m in range(B):
        dm = delta_r[m]
        PH[m] = np.exp(1j * (KX * dm[0] + KY * dm[1] + KZ * dm[2]))
    PH_k = np.moveaxis(PH, 0, 3)                          # (Nx,Ny,Nz,B)
    PHc = np.conj(PH)                                     # e^{-i k.dr_m}

    Vpp6 = Vpp.reshape(B, nd, B, nd)                      # [m,ab,m',cd]

    def _fluctuation(X_grid):
        """Y^fl: FFT convolution, O(B) factorized (spec S4.5)."""
        Z = PH[:, None, :, :, :] * X_grid[None, :, :, :, :]   # (B,nd,Nx,Ny,Nz)
        Z_r = _bk.spatial_ifftn(Z, axes=(2, 3, 4))
        Z_r_m = np.moveaxis(Z_r, (0, 1), (3, 4))              # (Nx,Ny,Nz,B,nd)
        # W_r[x,y,z,m,ab] = sum_{m',cd} F_hat[...,m,ab,m',cd] Z_r[...,m',cd]
        W_r = np.einsum('xyzaAbB,xyzbB->xyzaA', F_r6, Z_r_m)  # (Nx,Ny,Nz,B,nd)
        W_k = _bk.spatial_fftn(W_r, axes=(0, 1, 2))
        return -np.einsum('xyzm,xyzma->xyza', PH_k, W_k)      # (Nx,Ny,Nz,nd)

    def _instantaneous(X_grid):
        """Y^pp: low-rank bare particle-particle part, no convolution."""
        A_coef = (1.0 / N) * np.einsum('mxyz,cxyz->mc', PHc, X_grid)  # (B,nd)
        # Bcoef[m',ab] = sum_{m,cd} Vpp[(m,ab),(m',cd)] A_coef[m,cd]
        Bcoef = np.einsum('aAbB,aB->bA', Vpp6, A_coef)                # (B,nd)
        return -0.5 * np.einsum('mxyz,ma->xyza', PH, Bcoef)    # (Nx,Ny,Nz,nd)

    def _make_matvec(part):
        def matvec(v):
            phi = np.asarray(v).reshape(norb, norb, Nx, Ny, Nz)
            phi_flat = phi.reshape(nd, N)                # [(e,f),k]
            X = np.einsum('adk,dk->ak', G2p, phi_flat)   # (nd, N)  a=(c,d)
            X_grid = X.reshape(nd, Nx, Ny, Nz)

            if part == "fluctuation":
                Y = _fluctuation(X_grid)
            elif part == "instantaneous":
                Y = _instantaneous(X_grid)
            else:
                Y = _fluctuation(X_grid) + _instantaneous(X_grid)

            Y = Y.reshape(Nx, Ny, Nz, norb, norb)
            Y = np.moveaxis(Y, (3, 4), (0, 1))           # (norb,norb,Nx,Ny,Nz)
            return Y.ravel()
        return matvec

    ops = {part: LinearOperator((vec_size, vec_size),
                                matvec=_make_matvec(part), dtype=complex)
           for part in ("full", "fluctuation", "instantaneous")}
    return ops, vec_size


# ---------------------------------------------------------------------------
# bond bubble memory contract  --  read together with sc._bond_memory_estimate
# ---------------------------------------------------------------------------
#
# The resource preflight (``sc._bond_resource_preflight``) refuses a run whose
# estimated peak exceeds ``bond_memory_cap_gb``. That guard is only worth
# anything if the estimate is not an UNDERCOUNT, so the shared kernel's
# ``bubble.bond_bubble_static`` (the PRODUCTION streaming assembly, the sole
# implementation since the unified-bubble-kernel series' Task 11 legacy
# deletion) is written to a fixed, documented working set: every temporary is
# released (``del`` / in-place rescaling) as soon as it is consumed, and the
# two counts below are imported FROM ``hwave.solver.bubble`` -- which owns
# the canonical definition -- so the value actually used by
# ``sc._bond_memory_estimate`` always matches what the PRODUCTION path
# allocates. ``tests/test_sc_bond.py`` measures the real
# ``bubble.bond_bubble_static`` peak against the budget.
#
# The per-pair high-water mark is inside ``tau_to_boson(chi0_qt)`` in the
# channel-pair loop:
#
#   ``norb**4``-sized, on the ``(N_q, nmat)`` grid (BOND_BUBBLE_N4_BUFFERS;
#   the shared kernel's per-pair generator, ``bubble._iter_bond_channel_
#   pairs``, matches this same count -- its own risk of an extra buffer
#   (the previous pair's block staying alive across ``yield`` in the
#   generator's suspended frame) is closed by an explicit ``del`` right
#   after the ``yield``, documented beside the constant in
#   ``hwave.solver.bubble``):
#     1. ``chi0_qt``            -- the q,tau bubble being transformed
#     2. ``arr * omg``          -- tau_to_boson's internal phase-multiplied copy
#     3. the ``ifft`` output    -- becomes ``chi0_qw`` (rescaled IN PLACE)
#   (``chi0_rt`` is ``del``-ed right after ``spatial_fftn`` consumes it, and
#   the previous iteration's ``chi0_qw`` right after its block is stored, so
#   neither is alive here. ``spatial_fftn``'s own axis-by-axis pipeline holds
#   at most the same three.)
#
#   ``norb**2``-sized, on the ``(N_q, nmat)`` grid (BOND_BUBBLE_N2_BUFFERS,
#   the input ``green`` NOT included -- the preflight budgets that separately
#   as ``carrier_bytes``):
#     1. ``green_rev``          -- G(-r,-tau) (an np.flip VIEW of one np.roll)
#     2. ``G_fwd_sgn``          -- G(r,tau)*sgn(tau)
#     3. one transient          -- covers ``green_rev_shifted`` in the loop and,
#                                  in the setup phase (where no ``norb**4``
#                                  buffer exists yet), the fft temporaries of
#                                  ``fermion_to_tau``/``spatial_ifftn``
#   (``green_kt`` and ``green_rt`` are ``del``-ed as soon as they are consumed.)
#
# The output ``chi_bar`` is ``(N_q, ND, ND)`` and is already one of the
# ``sc._BOND_N_Q_ARRAYS`` q-resolved arrays, so it is not counted twice here.
BOND_BUBBLE_N4_BUFFERS = _bubble.BOND_BUBBLE_N4_BUFFERS
BOND_BUBBLE_N2_BUFFERS = _bubble.BOND_BUBBLE_N2_BUFFERS


def _sc_to_canonical(green):
    """``(p, p, Nx, Ny, Nz, nmat) -> (1, nmat, Nx*Ny*Nz, p, p)``.

    A PRIVATE, local copy of the same transpose/reshape/block-axis-insert
    ``sc._green_sc_to_canonical`` uses -- duplicated rather than imported
    because ``sc.py`` already imports this module (importing back would be
    circular); :func:`bond_bubble` (below) is the only caller.
    """
    p, p2, Nx, Ny, Nz, nmat = green.shape
    if p2 != p:
        raise ValueError(
            "bond_bubble: green's two orbital axes must both equal norb; "
            "got shape {}".format(green.shape))
    transposed = green.transpose(5, 2, 3, 4, 0, 1)      # (nmat,Nx,Ny,Nz,p,p)
    canonical = transposed.reshape(nmat, Nx * Ny * Nz, p, p)[np.newaxis, ...]
    return np.ascontiguousarray(canonical)


def bond_bubble(green, bond_set, beta):
    """Compatibility wrapper for the enlarged bond-channel bubble
    ``chi_bar(q; Delta r, Delta r')``.

    SAME signature and return as the original implementation (numerically
    pinned against it while both existed -- see
    ``tests.test_bubble_kernel.TestBubbleOldVsNewBondStatic``, now a
    self-consistency gate against the shared kernel's other independent
    oracles): converts ``green`` to the shared kernel's canonical layout and
    delegates to :func:`hwave.solver.bubble.bond_bubble_static`. TAIL-OFF
    ALWAYS -- this signature carries no tail argument, so it can only
    reproduce the original raw finite-``nmat`` sum; the production
    tail-corrected bond path (``sc._build_bond_operator``) calls
    ``bubble.bond_bubble_static`` directly with the Green carrier's
    ``deflated_kw``/``green0_tail`` instead of going through this function.
    The ``np.asarray`` coercion below is a CPU-compatibility shim kept for
    this established signature; the production path (sc.py's Green
    carrier into ``bubble.bond_bubble_static``) is device-clean and is
    what carries the GPU claims.

    Parameters
    ----------
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat)
        k-space, Matsubara-frequency Green's function in the ``sc.py``
        layout/normalization (see ``hwave.sc._build_bond_green`` /
        ``hwave.sc._load_flex_green``): ``iomega_n = (2n+1-nmat)*pi/beta``.
    bond_set : ResolvedInteractionSet
        Bond-channel topology (``delta_r``, ``n_channels``) as returned by
        :func:`resolve_interactions`.
    beta : float
        Inverse temperature.

    Returns
    -------
    chi_bar : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        Static (zero bosonic frequency) enlarged bubble; ``ND = norb**2 * B``.
    """
    green = np.asarray(green)
    if green.ndim != 6:
        raise ValueError(
            "bond_bubble: green must have shape (norb, norb, Nx, Ny, Nz, "
            "nmat); got ndim={} shape={}".format(green.ndim, green.shape)
        )
    norb, norb2, Nx, Ny, Nz, nmat = green.shape
    if norb2 != norb:
        raise ValueError(
            "bond_bubble: green's two orbital axes must both equal norb; "
            "got shape {}".format(green.shape)
        )

    canonical = _sc_to_canonical(green)
    chi_bar_flat = _bubble.bond_bubble_static(
        canonical, None, beta, bond_set, spatial_shape=(Nx, Ny, Nz))
    nD = chi_bar_flat.shape[-1]
    return chi_bar_flat.reshape(Nx, Ny, Nz, nD, nD)


# ===========================================================================
# lambda attribution on the symmetrized (Hermitian) kernel  --  spec S4.5
# ===========================================================================
#
# The bond kernel factorizes as ``K = -Gamma W`` where ``W`` is the static pair
# weight ``w = GG`` (multiplication by ``G2`` at each k, applied first inside
# every matvec) and ``Gamma`` is the pairing vertex (the (1/N)-normalized
# k,k'-convolution).  The spec's symmetrized kernel is the ``sqrt(w)``
# similarity
#
#     K~ = -sqrt(W) Gamma sqrt(W) = sqrt(W) K sqrt(W)^{-1},
#
# which is Hermitian (hence has a real spectrum) exactly when ``W`` is
# Hermitian positive semi-definite and ``Gamma`` is Hermitian.  Two identities
# make this usable without ever forming ``sqrt(W)``:
#
#   * ``K~`` is Hermitian  <=>  ``M = W K`` is Hermitian
#     (``sqrt(W) K~ sqrt(W) = W K`` and ``sqrt(W) K~^dag sqrt(W) = K^dag W``),
#     so the runtime residual is measured on ``M``.  The equivalence is EXACT
#     ONLY AT ZERO RESIDUAL: writing ``D = K~ - K~^dag``, one has
#     ``M - M^dag = sqrt(W) D sqrt(W)``, hence
#     ``w_min ||D|| <= ||M - M^dag|| <= w_max ||D||``.
#     At a FINITE tolerance the two criteria therefore differ by ``W``'s
#     conditioning, and a plain absolute tolerance on ``||M - M^dag||`` is
#     effectively ``tol / w_min`` on the quantity that matters (a violation
#     living on a small-``w`` index would pass).  :func:`
#     check_hermitian_preconditions` therefore applies the EXACT two-sided
#     rescaling ``K~ = W^-1/2 M W^-1/2`` (:func:`_metric_inverse_sqrt` /
#     :func:`_refer_to_ktilde`) and compares ``||K~ - K~^dag||`` against
#     ``max(atol, rtol ||K~||)`` -- BOTH sides referred to ``K~``.  Referring
#     only the residual (by the loose ``/w_min`` upper bound) while keeping the
#     ``||M||`` scale would make the relative criterion ``w_min`` times
#     stricter than requested, i.e. tighten it with beta and refuse legitimate
#     kernels; both numbers are reported.
#   * the Rayleigh quotient of ``K~`` at ``v~ = sqrt(W) v`` equals the
#     ``W``-weighted Rayleigh quotient of ``K`` at ``v``:
#     ``<v~|K~|v~>/<v~|v~> = <v|W K|v>/<v|W|v>``,
#     which is what :func:`attribute_lambda` evaluates for each part.

def pair_weight(green, beta, *, g2_tail=False):
    """Static pair weight ``w = GG`` as the Hermitian metric of the kernel.

    This is the SAME object the Eliashberg operator applies to the gap before
    the vertex convolution (``X = G2 phi`` in :func:`make_bond_kernel` and in
    ``sc._make_kernel_operator``), written as a matrix in the orbital-pair
    index so it can be used as an inner-product metric:

        ``W[k][(c,d),(e,f)] = G2_{c,e,d,f}(k)``,  ``G2 = T sum_n G(k,w_n)
        G(-k,-w_n)``.

    Parameters
    ----------
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat), complex
    beta : float

    Returns
    -------
    W : ndarray, shape (Nx, Ny, Nz, nd, nd), complex
        ``nd = norb**2``. For a centrosymmetric spin-degenerate input this is
        Hermitian and positive semi-definite (single band: real and ``>= 0``)
        -- the v1 precondition checked by
        :func:`check_hermitian_preconditions`.
    """
    green = _validate_green_beta(green, beta, "pair_weight")
    norb = green.shape[0]
    nd = norb * norb
    Nx, Ny, Nz = green.shape[2], green.shape[3], green.shape[4]
    G2 = _g2_from_green(green, beta, tail=g2_tail)  # (i,j,l,m,Nx,Ny,Nz)
    W = G2.transpose(0, 2, 1, 3, 4, 5, 6).reshape(nd, nd, Nx, Ny, Nz)
    return np.moveaxis(W, (0, 1), (3, 4)).copy()


def _apply_metric(weight, V):
    """Apply the metric ``W`` to a stack of row vectors ``V`` (``(d, n)``).

    ``weight`` may be

    * ``None`` -- identity metric;
    * a callable ``W(V) -> (d, n)``;
    * a 1-D array of length ``n`` -- a diagonal metric;
    * a 5-D array ``(Nx, Ny, Nz, nd, nd)`` -- the per-k orbital-pair blocks
      returned by :func:`pair_weight`, acting on vectors laid out as
      ``(norb, norb, Nx, Ny, Nz)`` raveled (the gap layout).
    """
    V = np.atleast_2d(np.asarray(V))
    if weight is None:
        return V
    if callable(weight):
        return np.atleast_2d(weight(V))
    W = np.asarray(weight)
    if W.ndim == 1:
        return V * W[np.newaxis, :]
    if W.ndim == 5:
        Nx, Ny, Nz, nd, nd2 = W.shape
        if nd != nd2:
            raise ValueError(
                "_apply_metric: the metric's last two axes must be square; "
                "got shape {}".format(W.shape))
        d, n = V.shape
        if n != nd * Nx * Ny * Nz:
            raise ValueError(
                "_apply_metric: vector length {} does not match the metric "
                "shape {} (expected nd*Nx*Ny*Nz = {})".format(
                    n, W.shape, nd * Nx * Ny * Nz))
        Vg = V.reshape(d, nd, Nx, Ny, Nz)
        out = np.einsum('xyzab,dbxyz->daxyz', W, Vg)
        return out.reshape(d, n)
    raise ValueError(
        "_apply_metric: unsupported metric of shape {}; expected None, a "
        "callable, a 1-D diagonal, or a 5-D (Nx, Ny, Nz, nd, nd) array"
        .format(W.shape))


def _metric_inner(weight, A, B):
    """Gram matrix ``G[i,j] = <a_i, b_j>_W = a_i^dag W b_j`` for row stacks."""
    A = np.atleast_2d(np.asarray(A))
    WB = _apply_metric(weight, B)
    return np.conj(A) @ WB.T


def check_hermitian_preconditions(operator, weight, *, atol=1.0e-8,
                                  rtol=1.0e-8, dense_limit=1024, n_probe=16,
                                  seed=20260726, label="bond", parts=None):
    """Check the v1 runtime preconditions of the Hermitian static path (S4.5).

    The static solver reports ``lambda`` as a real eigenvalue/Rayleigh quotient
    of the symmetrized kernel ``K~ = -sqrt(w) Gamma sqrt(w)``. That is only
    legitimate when the supplied Green function has the inversion/time-reversal
    symmetry which makes the pair weight ``w = GG`` **Hermitian and positive
    semi-definite** (single band: real and ``>= 0``) and the vertex ``Gamma``
    Hermitian. Both are verified here; a violation **raises** -- v1 does not
    fall back to a non-Hermitian biorthogonal solve (documented follow-up).

    Parameters
    ----------
    operator : LinearOperator
        The FULL Eliashberg kernel ``K`` (not one of the parts).
    weight : ndarray or None
        The pair weight metric, see :func:`pair_weight` / :func:`_apply_metric`.
    atol, rtol : float, optional
        A residual ``r`` passes when ``r <= max(atol, rtol * scale)`` with
        ``scale`` the corresponding Frobenius/element scale.

        **Tolerance convention for the KERNEL check (both sides referred to
        ``K~``).** The runtime residual is measured on ``M = W K``, but the
        criterion is about ``K~``; since ``M = sqrt(W) K~ sqrt(W)``, the check
        applies the EXACT two-sided rescaling ``K~ = W^-1/2 M W^-1/2`` and
        compares ``||K~ - K~^dag||_F`` against ``max(atol, rtol ||K~||_F)``.
        Both quantities live on the same object, so ``rtol`` means exactly
        what it says -- a relative asymmetry of the symmetrized kernel -- and
        the criterion does NOT tighten as ``min w`` shrinks with beta (review
        fix I1). The ``M``-referred numbers are reported alongside. When the
        metric carries no spectral information (``None``/callable) the check
        stays on ``M``, which is then ``K~`` itself.
    dense_limit : int, optional
        Build ``K`` densely (``vec_size`` matvecs, ``vec_size**2`` complex
        entries) when ``vec_size <= dense_limit``, giving an EXACT Hermiticity
        residual; the default keeps that below ~16 MB / ~1k matvecs.

        **Cost of the dense branch (review fix I-4).** The ``vec_size``
        matvecs above build the FULL kernel only; when ``parts`` is given
        (the normal production call from ``sc._build_bond_operator``, which
        always passes ``{"pp": ..., "fl": ...}``) the part-wise check below
        densely builds EACH part too, for a TOTAL of ``3 * vec_size`` bond
        kernel matvecs before the eigensolve even starts. With the default
        ``dense_limit = 1024``, a 32x32 single-band production grid
        (``vec_size = norb**2 * Nx * Ny * Nz = 1 * 32 * 32 * 1 = 1024``) sits
        EXACTLY at the limit, so this precondition check alone costs ~3072
        bond-kernel matvecs -- each itself an ``O(B)`` FFT-factorized
        operation (:func:`_bond_kernel_operators`), not a cheap dense
        matrix-vector product. This is a one-time up-front cost per
        Eliashberg solve (not per SCF/eigensolver iteration), but it is not
        negligible at this grid size; ``bond_precondition_dense_limit`` is
        the escape hatch (forces the cheaper randomized-probe estimator
        instead, see below).

        Above ``dense_limit``, a randomized (Hutchinson-style) estimate of the same Frobenius norm is
        used: ``E|a^dag M b|^2 = ||M||_F^2`` for unit-variance complex Gaussian
        probes, from ``n_probe`` probe pairs seeded by ``seed``
        (deterministic run to run). The estimator cannot produce a FALSE
        ALARM -- an exactly Hermitian ``M`` gives ``a^dag(M-M^dag)b = 0`` for
        every probe -- it can only be less sensitive to a small violation than
        the dense residual. Each probe pair costs FOUR matvecs: two for the
        ``M``-referred diagnostic and two for the ``K~``-referred criterion
        (the same probes rescaled by ``W^-1/2``).
    n_probe, seed : int, optional
    label : str, optional
        Prefix used in the error messages.
    parts : dict[str, LinearOperator], optional
        Additive pieces of ``operator`` (e.g. ``{"pp": K_pp, "fl": K_fl}``).
        The S4.5 attribution reports a REAL Rayleigh quotient for each part
        separately, which needs each part -- not only their sum -- Hermitian
        in the ``W`` metric; a violation that cancels between the parts is
        invisible in the full-kernel residual. Checked in the DENSE branch
        only (the probe estimator would cost another ``n_probe`` matvecs per
        part for a diagnostic that is not a hard precondition) and reported as
        ``part_hermiticity_residual_<name>`` /
        ``part_hermiticity_residual_ktilde_<name>``. A violation WARNS rather
        than raising: the reported ``lambda`` (the full kernel) is still
        legitimate, only its pp/fl split becomes questionable -- and
        :func:`attribute_lambda` independently reports the imaginary parts.

    Returns
    -------
    diagnostics : dict
        ``weight_hermiticity_residual``, ``weight_min_eigenvalue``,
        ``weight_scale``, ``kernel_hermiticity_residual`` (absolute Frobenius
        norm of ``W K - (W K)^dag``), ``kernel_hermiticity_residual_ktilde``
        (the EXACT ``||K~ - K~^dag||_F``, the quantity compared against the
        tolerance), ``kernel_hermiticity_scale``/``..._scale_ktilde``,
        ``kernel_hermiticity_relative``/``..._relative_ktilde``,
        ``kernel_hermiticity_referral``, ``method`` ("dense" or "probe"),
        and, when ``parts`` is given and the dense branch ran, one
        ``part_hermiticity_residual*_<name>`` entry per part.

    Raises
    ------
    ValueError
        If the pair weight is not Hermitian / not positive semi-definite, or
        the symmetrized kernel is not Hermitian, within tolerance.
    """
    diagnostics = {}

    # --- 1/2: the pair weight w = GG must be Hermitian and >= 0 ------------
    if weight is None or callable(weight):
        diagnostics["weight_hermiticity_residual"] = 0.0
        diagnostics["weight_min_eigenvalue"] = float("nan")
        diagnostics["weight_scale"] = float("nan")
    else:
        W = np.asarray(weight)
        if W.ndim == 1:
            # diagonal metric: "Hermitian" = real, "PSD" = every entry >= 0
            herm_res = float(np.abs(np.imag(W)).max())
            scale = float(np.abs(W).max())
            min_eig = float(np.real(W).min())
        elif W.ndim == 5:
            herm_res = float(np.abs(W - np.conj(np.swapaxes(W, -1, -2))).max())
            scale = float(np.abs(W).max())
            Wh = 0.5 * (W + np.conj(np.swapaxes(W, -1, -2)))
            min_eig = float(np.linalg.eigvalsh(Wh).min())
        else:
            raise ValueError(
                "check_hermitian_preconditions: unsupported metric shape {}"
                .format(W.shape))
        tol = max(atol, rtol * max(scale, 1.0))
        diagnostics["weight_hermiticity_residual"] = herm_res
        diagnostics["weight_min_eigenvalue"] = min_eig
        diagnostics["weight_scale"] = scale
        if herm_res > tol:
            raise ValueError(
                "[{}] the static pair weight w = GG is NOT Hermitian "
                "(max |w - w^dag| = {:.3e} > tol {:.3e}). The v1 Hermitian "
                "path requires a Green function with the inversion/"
                "time-reversal symmetry that makes w real (single band) / "
                "Hermitian (multi-orbital) and >= 0, so that "
                "K~ = -sqrt(w) Gamma sqrt(w) has a real spectrum. v1 does NOT "
                "fall back to a non-Hermitian biorthogonal solve (deferred "
                "follow-up): supply a symmetric green.npz or disable "
                "bond_channels.".format(label, herm_res, tol))
        if min_eig < -tol:
            raise ValueError(
                "[{}] the static pair weight w = GG is NOT positive "
                "semi-definite (min eigenvalue = {:.3e} < -tol {:.3e}). The "
                "sqrt(w) similarity that symmetrizes the Eliashberg kernel "
                "does not exist, so the reported lambda would not be a real "
                "Rayleigh quotient. v1 errors instead of running a "
                "non-Hermitian biorthogonal solve (deferred follow-up)."
                .format(label, min_eig, tol))

    # --- 3: Hermiticity of the symmetrized kernel, via M = W K ------------
    #
    # Tolerance convention (review fix I1): BOTH the residual and the scale it
    # is compared against are referred to K~, by the EXACT two-sided rescaling
    # K~ = W^-1/2 M W^-1/2 (see the module comment above) -- never by the
    # ||M - M^dag||/w_min upper bound. Referring only the residual (as the
    # first version of this check did) makes the relative criterion w_min
    # times stricter than requested, i.e. it tightens with beta and refuses
    # legitimate kernels at large beta; referring neither lets a violation
    # localized on a small-w index pass. When W carries no spectral
    # information (None / callable metric) the check stays on M, which for
    # W = I IS K~, so the two sides are consistent there too.
    Winv_sqrt = _metric_inverse_sqrt(weight)
    n = operator.shape[0]
    if n <= dense_limit:
        K = np.zeros((n, n), dtype=complex)
        for i in range(n):
            e = np.zeros(n, dtype=complex)
            e[i] = 1.0
            K[:, i] = operator.matvec(e)
        # columns of M = W K are W (K e_j)
        M = _apply_metric(weight, K.T).T
        resid = float(np.linalg.norm(M - np.conj(M.T)))
        scale = float(np.linalg.norm(M))
        N = _refer_to_ktilde(M, Winv_sqrt)
        resid_kt = float(np.linalg.norm(N - np.conj(N.T)))
        scale_kt = float(np.linalg.norm(N))
        method = "dense"
    else:
        rng = np.random.default_rng(seed)
        res_acc = 0.0
        scale_acc = 0.0
        res_kt_acc = 0.0
        scale_kt_acc = 0.0
        for _ in range(int(n_probe)):
            a = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2.0)
            b = (rng.normal(size=n) + 1j * rng.normal(size=n)) / np.sqrt(2.0)
            # <a|M|b> and <b|M|a> with M = W K, then the same probes on
            # K~ = W^-1/2 M W^-1/2 (i.e. the probes rescaled by W^-1/2)
            for tag in ("M", "kt"):
                if tag == "M":
                    u, v = a, b
                else:
                    if Winv_sqrt is None:
                        continue
                    u = _apply_metric(Winv_sqrt, a[np.newaxis, :])[0]
                    v = _apply_metric(Winv_sqrt, b[np.newaxis, :])[0]
                Ku = operator.matvec(u)
                Kv = operator.matvec(v)
                Wu = _apply_metric(weight, u[np.newaxis, :])[0]
                Wv = _apply_metric(weight, v[np.newaxis, :])[0]
                aMb = np.vdot(Wu, Kv)
                bMa = np.vdot(Wv, Ku)
                if tag == "M":
                    res_acc += abs(aMb - np.conj(bMa)) ** 2
                    scale_acc += abs(aMb) ** 2
                else:
                    res_kt_acc += abs(aMb - np.conj(bMa)) ** 2
                    scale_kt_acc += abs(aMb) ** 2
        resid = float(np.sqrt(res_acc / n_probe))
        scale = float(np.sqrt(scale_acc / n_probe))
        if Winv_sqrt is None:
            resid_kt, scale_kt = resid, scale
        else:
            resid_kt = float(np.sqrt(res_kt_acc / n_probe))
            scale_kt = float(np.sqrt(scale_kt_acc / n_probe))
        method = "probe"

    rel = resid / scale if scale > 0.0 else 0.0
    rel_kt = resid_kt / scale_kt if scale_kt > 0.0 else 0.0
    diagnostics["kernel_hermiticity_residual"] = resid
    diagnostics["kernel_hermiticity_scale"] = scale
    diagnostics["kernel_hermiticity_relative"] = rel
    diagnostics["kernel_hermiticity_residual_ktilde"] = resid_kt
    diagnostics["kernel_hermiticity_scale_ktilde"] = scale_kt
    diagnostics["kernel_hermiticity_relative_ktilde"] = rel_kt
    diagnostics["kernel_hermiticity_referral"] = (
        "ktilde" if Winv_sqrt is not None else "M (no spectral metric)")
    diagnostics["method"] = method

    tol_k = max(atol, rtol * scale_kt)
    if resid_kt > tol_k:
        raise ValueError(
            "[{}] the symmetrized Eliashberg kernel K~ = -sqrt(w) Gamma "
            "sqrt(w) is NOT Hermitian: ||K~ - K~^dag||_F = {:.3e} "
            "(K~-scale {:.3e}, relative {:.3e}) exceeds tol {:.3e} "
            "[{} check; the same residual on M = W K is {:.3e} with scale "
            "{:.3e}]. Its spectrum is then not guaranteed real and the "
            "lambda = lambda^pp + lambda^fl attribution as REAL Rayleigh "
            "quotients would be meaningless. v1 errors instead of running a "
            "non-Hermitian biorthogonal solve (documented follow-up)."
            .format(label, resid_kt, scale_kt, rel_kt, tol_k, method,
                    resid, scale))

    # --- 4: part-wise Hermiticity (dense branch only; diagnostic) ---------
    if parts and method == "dense":
        for name in sorted(parts):
            op_part = parts[name]
            Kp = np.zeros((n, n), dtype=complex)
            for i in range(n):
                e = np.zeros(n, dtype=complex)
                e[i] = 1.0
                Kp[:, i] = op_part.matvec(e)
            Mp = _apply_metric(weight, Kp.T).T
            r_p = float(np.linalg.norm(Mp - np.conj(Mp.T)))
            s_p = float(np.linalg.norm(Mp))
            Np = _refer_to_ktilde(Mp, Winv_sqrt)
            r_p_kt = float(np.linalg.norm(Np - np.conj(Np.T)))
            s_p_kt = float(np.linalg.norm(Np))
            diagnostics["part_hermiticity_residual_" + name] = r_p
            diagnostics["part_hermiticity_residual_ktilde_" + name] = r_p_kt
            diagnostics["part_hermiticity_scale_" + name] = s_p
            tol_p = max(atol, rtol * max(s_p_kt, scale_kt))
            if r_p_kt > tol_p:
                logger.warning(
                    "[%s] the '%s' part of the Eliashberg kernel is NOT "
                    "Hermitian in the pair-weight metric "
                    "(||K~_%s - K~_%s^dag||_F = %.3e > tol %.3e) even though "
                    "the FULL kernel passes -- the violation cancels between "
                    "the parts. lambda itself is still a real Rayleigh "
                    "quotient, but its lambda^pp / lambda^fl split (spec "
                    "S4.5) is not trustworthy; check the imaginary parts "
                    "reported by attribute_lambda.",
                    label, name, name, name, r_p_kt, tol_p)

    return diagnostics


def _metric_inverse_sqrt(weight, rcond=None):
    """``W^-1/2`` for the metrics :func:`_apply_metric` understands, or ``None``.

    Returned in the SAME kind/shape as ``weight`` (1-D diagonal or 5-D per-k
    blocks) so it can be fed straight back to :func:`_apply_metric`. It is what
    refers a ``M = W K`` residual EXACTLY to ``K~ = W^-1/2 M W^-1/2`` (review
    fix I1), instead of the conditioning-dependent ``/w_min`` bound.

    Directions whose weight is below ``rcond * w_max`` (default
    ``sqrt(eps) w_max``) carry no pair weight -- the kernel's matching column
    is proportional to it, since ``K = -Gamma W`` -- and are mapped to zero (a
    pseudo-inverse), so a numerically null weight cannot amplify roundoff into
    a false alarm.

    ``None`` is returned for a ``None``/callable metric, whose spectrum is not
    available; the caller then keeps the check on ``M`` (which for ``W = I``
    IS ``K~``).
    """
    if weight is None or callable(weight):
        return None
    W = np.asarray(weight)
    if rcond is None:
        rcond = float(np.sqrt(np.finfo(float).eps))
    if W.ndim == 1:
        w = np.real(np.asarray(W, dtype=complex))
        if w.size == 0:
            return w
        cut = rcond * float(np.abs(w).max())
        out = np.zeros_like(w)
        good = w > cut
        out[good] = 1.0 / np.sqrt(w[good])
        return out
    if W.ndim == 5:
        Wh = 0.5 * (W + np.conj(np.swapaxes(W, -1, -2)))
        evals, evecs = np.linalg.eigh(Wh)
        if evals.size == 0:
            return None
        cut = rcond * float(np.abs(evals).max())
        inv = np.where(evals > cut,
                       1.0 / np.sqrt(np.where(evals > cut, evals, 1.0)), 0.0)
        return (evecs * inv[..., np.newaxis, :]) @ np.conj(
            np.swapaxes(evecs, -1, -2))
    return None


def _refer_to_ktilde(M, Winv_sqrt):
    """``K~ = W^-1/2 M W^-1/2`` for the dense ``M = W K`` (identity if None).

    ``_apply_metric(S, X)`` evaluates ``X S^T``, so ``S M`` is
    ``_apply_metric(S, M.T).T`` and the right multiplication uses
    ``A S = conj(conj(A) S^T)`` (``S`` is Hermitian).
    """
    if Winv_sqrt is None:
        return M
    left = _apply_metric(Winv_sqrt, M.T).T               # S M
    return np.conj(_apply_metric(Winv_sqrt, np.conj(left)))


def attribute_lambda(vec, weight, op_pp, op_fl, op_full=None,
                     atol_imag=1.0e-6):
    """Split ``lambda`` into its bare-pp and fluctuation parts (spec S4.5).

    Evaluates the REAL Rayleigh quotients of the symmetrized kernel,

        ``lambda^X = <v~|K~^X|v~> / <v~|v~> = <v|W K^X|v> / <v|W|v>``
        with ``v~ = sqrt(W) v``,

    for ``X in {pp, fl}``; because ``K = K^pp + K^fl`` the two add up to the
    full quotient exactly. Call
    :func:`check_hermitian_preconditions` first -- the quotients are only real
    (and the split only meaningful) on the Hermitian path.

    Parameters
    ----------
    vec : ndarray
        The converged eigenvector/gap (any shape; raveled internally).
    weight : ndarray or None
        The pair weight metric (:func:`pair_weight`).
    op_pp, op_fl : LinearOperator
        The instantaneous and fluctuation parts
        (:func:`make_bond_kernel_parts`).
    op_full : LinearOperator, optional
        The full kernel. When given, ``lambda`` is its own Rayleigh quotient
        and ``sum_residual`` reports ``|lambda - (lambda^pp + lambda^fl)|``
        (an operator-identity check); otherwise ``lambda`` is the sum.
    atol_imag : float, optional
        Tolerance on the imaginary part of each quotient (a diagnostic; it is
        reported, not enforced -- the enforcement is the precondition check).

    Returns
    -------
    dict with keys ``lambda``, ``lambda_pp``, ``lambda_fl``, ``imag``,
    ``imag_pp``, ``imag_fl``, ``sum_residual``, ``metric_norm``,
    ``imag_within_tol``.
    """
    v = np.asarray(vec).ravel()
    Wv = _apply_metric(weight, v[np.newaxis, :])[0]
    den = np.vdot(v, Wv)
    den_real = float(den.real)
    if den_real <= 0.0:
        raise ValueError(
            "attribute_lambda: the W-norm of the eigenvector is {:.3e} <= 0; "
            "the pair weight is not positive definite on this vector, so the "
            "Rayleigh quotient is not defined. Run "
            "check_hermitian_preconditions first.".format(den_real))

    def _quotient(op):
        q = np.vdot(Wv, op.matvec(v)) / den_real
        return float(q.real), float(q.imag)

    lam_pp, im_pp = _quotient(op_pp)
    lam_fl, im_fl = _quotient(op_fl)
    if op_full is not None:
        lam, im = _quotient(op_full)
    else:
        lam, im = lam_pp + lam_fl, im_pp + im_fl
    return {
        "lambda": lam,
        "lambda_pp": lam_pp,
        "lambda_fl": lam_fl,
        "imag": im,
        "imag_pp": im_pp,
        "imag_fl": im_fl,
        "sum_residual": abs(lam - (lam_pp + lam_fl)),
        "metric_norm": den_real,
        "imag_within_tol": bool(max(abs(im), abs(im_pp), abs(im_fl))
                                <= atol_imag * max(1.0, abs(lam))),
    }


# ===========================================================================
# invariant-subspace tracking across a V sweep  --  spec S7.7
# ===========================================================================

def cluster_eigenvalues(evals, deg_tol=1.0e-3):
    """Group eigenvalues into near-degenerate clusters (spec S7.7 (iii)).

    Sorted by descending real part; consecutive eigenvalues whose COMPLEX
    distance ``|lambda_i - lambda_j|`` is at most ``deg_tol`` join the same
    cluster (single-linkage), so a cluster is the numerically-degenerate
    invariant subspace the tracker follows.

    The distance is the full complex modulus, not ``|Re lambda_i -
    Re lambda_j|``: on the Hermitian static path the spectrum is real and the
    two agree, but a kernel that has drifted off that path (or an ARPACK
    result with a residual imaginary part) can produce eigenvalues sharing a
    real part while being physically distinct, and merging those into one
    "degenerate" subspace would make the tracker follow a span that is not
    invariant.

    Returns
    -------
    list[list[int]]
        Index lists into ``evals``, clusters ordered by descending leading
        real part.
    """
    evals = np.asarray(evals)
    if evals.size == 0:
        return []
    order = list(np.argsort(-evals.real))
    clusters = [[order[0]]]
    for idx in order[1:]:
        if abs(evals[idx] - evals[clusters[-1][-1]]) <= deg_tol:
            clusters[-1].append(idx)
        else:
            clusters.append([idx])
    return [[int(i) for i in c] for c in clusters]


def orthonormalize(vectors, weight=None, rank_tol=1.0e-12):
    """``W``-orthonormal basis of the span of ``vectors`` (rows).

    Uses the eigendecomposition of the Gram matrix ``G_ij = <v_i, v_j>_W``
    (Hermitian positive semi-definite), dropping directions whose Gram
    eigenvalue falls below ``rank_tol`` times the largest -- so a numerically
    rank-deficient set yields a smaller but still orthonormal basis.
    """
    V = np.atleast_2d(np.asarray(vectors, dtype=complex))
    G = _metric_inner(weight, V, V)
    G = 0.5 * (G + np.conj(G.T))
    w, U = np.linalg.eigh(G)
    if w.size == 0 or w.max() <= 0.0:
        return np.zeros((0, V.shape[1]), dtype=complex)
    keep = w > rank_tol * w.max()
    C = U[:, keep] / np.sqrt(w[keep])[np.newaxis, :]
    return C.T @ V


def subspace_similarity(A, B, weight=None):
    """Principal-angle overlap of two subspaces in the ``W`` metric (S7.7 iv).

    Both row stacks are ``W``-orthonormalized; the singular values of the
    overlap matrix ``O_ij = <a_i, b_j>_W`` are the cosines of the principal
    angles. The returned score is the MEAN SQUARED cosine over
    ``min(dim A, dim B)`` -- a basis-invariant subspace similarity in
    ``[0, 1]`` (1 = one subspace contains the other). Comparing whole
    subspaces (rather than one vector) is what makes the tracking survive an
    arbitrary rotation inside a degenerate multiplet.

    Note that the mean (rather than the minimum or the product) makes the
    score comparable between clusters of DIFFERENT dimension, which is what
    :func:`track_subspace` needs -- but it also means EXACT TIES are possible
    (e.g. a seed split evenly over two one-dimensional clusters gives 0.5 for
    both). :func:`track_subspace` resolves such a tie in favour of the FIRST
    cluster in lambda-descending order (it keeps the incumbent best on a
    strict ``>`` comparison, and iterates clusters in the order
    :func:`cluster_eigenvalues` returns them).

    Returns
    -------
    (score, cosines) : (float, ndarray)
    """
    QA = orthonormalize(A, weight)
    QB = orthonormalize(B, weight)
    if QA.shape[0] == 0 or QB.shape[0] == 0:
        return 0.0, np.zeros(0)
    O = _metric_inner(weight, QA, QB)
    s = np.linalg.svd(O, compute_uv=False)
    s = np.clip(s, 0.0, 1.0)
    return float(np.mean(s ** 2)), s


def harmonic_decomposition(vectors, basis, weight=None):
    """Decompose a subspace onto a set of named form factors (spec S7.7).

    For each basis function ``b`` (``W``-normalized) the reported weight is the
    fraction of ``b`` captured by the subspace,
    ``||P_S b||_W^2 = sum_j |<q_j, b>_W|^2`` with ``{q_j}`` a ``W``-orthonormal
    basis of the subspace -- a number in ``[0, 1]``, invariant under any
    rotation inside the subspace.

    Parameters
    ----------
    vectors : ndarray, shape (d, n)
    basis : dict[str, ndarray]
        Named form factors (any shape; raveled). The odd basis of S7.7 is
        ``{sin kx, sin ky, sin kx (cos kx - cos ky),
        sin ky (cos ky - cos kx)}``.
    weight : metric, optional

    Returns
    -------
    dict[str, float]
    """
    Q = orthonormalize(vectors, weight)
    out = {}
    for name, b in basis.items():
        bv = np.asarray(b, dtype=complex).ravel()[np.newaxis, :]
        nb = _metric_inner(weight, bv, bv)[0, 0].real
        if nb <= 0.0:
            out[name] = 0.0
            continue
        bv = bv / np.sqrt(nb)
        if Q.shape[0] == 0:
            out[name] = 0.0
            continue
        ov = _metric_inner(weight, Q, bv)[:, 0]
        out[name] = float(np.sum(np.abs(ov) ** 2))
    return out


def track_subspace(points, *, seed=None, weight=None, deg_tol=1.0e-3,
                   basis=None, capture_tol=0.5, n_ev=None, max_n_ev=None,
                   rank_tol=1.0e-12):
    """Track one invariant subspace across a sequence of solves (spec S7.7).

    **This is the public API for SWEEP DRIVERS, not for a single solve.**
    Tracking is inherently MULTI-POINT: it follows one invariant subspace from
    one ``V`` (or ``T``, or doping) point to the next, so it needs a sequence
    of solves and cannot be evaluated inside a single ``sc.calc_eliashberg``
    run, which has exactly one point. A driver that sweeps ``V`` -- e.g. the
    single-band Onari milestone in ``tests/test_bond_onari_milestone.py``, or
    any user script that loops ``calc_eliashberg`` -- calls this directly with
    the per-point eigenpairs. The SINGLE-POINT pieces of the same analysis (the
    odd-harmonic decomposition of the leading state and the degeneracy
    clustering of its spectrum) are wired into the solver itself behind
    ``[eliashberg] bond_diagnostics = true`` and land in the eigenvalue file's
    provenance block; see ``sc._bond_diagnostics_record``.

    At every point the eigenpairs are clustered into near-degenerate
    invariant subspaces (``deg_tol``); the tracked subspace is the cluster of
    maximum principal-angle overlap (in the ``sqrt(GG)`` metric) with the
    previous point's tracked subspace -- **not** the maximum overlap of a
    single eigenvector, which is ill-defined inside a degenerate multiplet and
    can hop onto a nearby non-degenerate branch. The first point is anchored by
    maximal overlap with ``seed`` (the odd f-seed at ``V = 0``).

    **Tiebreak.** The score is the mean squared principal cosine over
    ``min(dim A, dim B)`` (:func:`subspace_similarity`); clusters are visited
    in the order :func:`cluster_eigenvalues` returns them (descending leading
    real part) and the incumbent best is only replaced on a STRICT
    improvement, so an exact tie is resolved to the first cluster in
    lambda-descending order.

    Parameters
    ----------
    points : sequence
        One entry per sweep point, each either

        * ``(evals, evecs)`` -- ``evecs`` a row stack (or sequence) of
          eigenvectors, ONE PER ROW (``scipy.sparse.linalg.eigs`` returns them
          as COLUMNS: pass ``vecs.T``), each of any shape (raveled here), or
        * a callable ``f(n_ev) -> (evals, evecs)`` -- then ``n_ev`` is
          **adaptive**: if the tracked subspace is not captured
          (best overlap ``< capture_tol``) the point is re-solved with twice as
          many eigenpairs, up to ``max_n_ev``.
    seed : ndarray, optional
        Anchor vector for the first point (e.g. the ``f_x`` gap seed). When
        ``None`` the leading cluster is taken.
    weight : metric, optional
        The ``sqrt(GG)`` metric (:func:`pair_weight`); ``None`` = identity.
    deg_tol : float, optional
        Degeneracy clustering tolerance (spec default ``1e-3``).
    basis : dict[str, ndarray], optional
        Odd harmonics for the per-point decomposition
        (:func:`harmonic_decomposition`).
    capture_tol : float, optional
        Minimum subspace overlap that counts as "captured".
    n_ev : int, optional
        Initial number of eigenpairs requested from callable points
        (default 6, the spec's ``num_eigenvalues`` default).
    max_n_ev : int, optional
        Ceiling for the adaptive escalation (default ``4 * n_ev``).

    Returns
    -------
    list[dict]
        One record per point: ``lambda`` (leading eigenvalue of the tracked
        subspace), ``eigenvalues``, ``dim``, ``cluster`` (indices), ``overlap``
        (subspace similarity with the previous point; with ``seed`` for the
        first), ``principal_cosines``, ``captured``, ``truncated``, ``n_ev``,
        ``harmonics``, and ``vectors`` (the ``W``-orthonormal basis of the
        tracked subspace).

    Notes
    -----
    **Window truncation.** A cluster that reaches the LOW end of the supplied
    eigenvalue window may be INCOMPLETE: the eigensolver was asked for
    ``n_ev`` eigenpairs and a degenerate multiplet can straddle that cut, so
    only part of the invariant subspace is returned. What comes back is then
    an ARBITRARY slice of the true eigenspace -- and with an ARPACK window
    started from an unpinned ``v0`` (``scipy.sparse.linalg.eigs`` reuses a
    process-global saved seed that every earlier call advances) it is a slice
    that depends on how many eigensolves the process happened to run before,
    which makes ``overlap``, ``harmonics`` and ``dim`` non-reproducible. Such
    a cluster is reported with ``truncated = True`` and, for CALLABLE points,
    triggers the same ``n_ev`` escalation as a capture failure -- even when
    the overlap happened to clear ``capture_tol``. Callers that pass fixed
    eigenpairs (a sweep driver holding precomputed spectra) get the flag and
    a warning, and should widen their own window.
    """
    if n_ev is None:
        n_ev = 6
    if max_n_ev is None:
        max_n_ev = 4 * n_ev

    records = []
    prev_basis = None
    if seed is not None:
        prev_basis = np.atleast_2d(np.asarray(seed, dtype=complex).ravel())

    for t, point in enumerate(points):
        n_ev_cur = int(n_ev)
        while True:
            if callable(point):
                evals, evecs = point(n_ev_cur)
            else:
                evals, evecs = point
            evals = np.asarray(evals)
            V = np.atleast_2d(np.asarray(
                [np.asarray(x, dtype=complex).ravel() for x in evecs]))
            n_ev_used = len(evals)

            clusters = cluster_eigenvalues(evals, deg_tol=deg_tol)
            best_score, best_cos, best_cluster = -1.0, np.zeros(0), None
            for cl in clusters:
                if prev_basis is None:
                    # no anchor: take the leading cluster
                    best_score, best_cos, best_cluster = 1.0, np.ones(1), cl
                    break
                score, cos = subspace_similarity(V[cl], prev_basis, weight)
                # strict ">" => an exact tie keeps the incumbent, i.e. the
                # first cluster in lambda-descending order (documented above)
                if score > best_score:
                    best_score, best_cos, best_cluster = score, cos, cl

            captured = bool(best_score >= capture_tol)
            # A cluster that touches the LOW end of the returned window may
            # have been cut in half by the eigensolver's n_ev truncation; see
            # the "Window truncation" note in the docstring. Escalate on that
            # too, not only on a capture failure -- a truncated multiplet can
            # clear capture_tol by luck while its basis (hence overlap,
            # harmonics and dim) is solver-state dependent.
            truncated = bool(
                best_cluster is not None and evals.size > 0
                and int(np.argmin(evals.real)) in best_cluster)
            if ((not captured or truncated) and callable(point)
                    and n_ev_cur < max_n_ev):
                n_ev_cur = min(2 * n_ev_cur, int(max_n_ev))
                continue
            break

        if best_cluster is None:
            raise ValueError(
                "track_subspace: point {} supplied no eigenpairs to track"
                .format(t))

        sub = V[best_cluster]
        Q = orthonormalize(sub, weight, rank_tol=rank_tol)
        lam_cluster = evals[best_cluster]
        record = {
            "lambda": float(np.max(lam_cluster.real)),
            "eigenvalues": lam_cluster,
            "dim": len(best_cluster),
            "cluster": list(best_cluster),
            "overlap": float(best_score),
            "principal_cosines": best_cos,
            "captured": captured,
            "truncated": truncated,
            "n_ev": int(n_ev_used),
            "vectors": Q,
        }
        if basis:
            record["harmonics"] = harmonic_decomposition(sub, basis, weight)
        if not captured:
            logger.warning(
                "track_subspace: point %d -- best subspace overlap %.3f is "
                "below capture_tol %.3f even at n_ev = %d (max_n_ev = %d). "
                "The tracked state may have left the computed eigenvalue "
                "window; increase num_eigenvalues.",
                t, best_score, capture_tol, n_ev_used, max_n_ev)
        if truncated:
            logger.warning(
                "track_subspace: point %d -- the tracked cluster (lambda = "
                "%.6g, dim %d) reaches the low end of the %d returned "
                "eigenpairs, so the degenerate multiplet may be TRUNCATED by "
                "the eigenvalue window. Its basis -- and therefore overlap, "
                "dim and harmonics -- is then an arbitrary slice of the true "
                "invariant subspace and is not reproducible; increase "
                "num_eigenvalues.",
                t, record["lambda"], record["dim"], n_ev_used)
        records.append(record)
        prev_basis = Q if Q.shape[0] > 0 else sub

    return records


# =============================================================================
# Phase W, Task 2 -- the master transverse bond topology
# =============================================================================
#
# APPEND-ONLY REGION. Nothing above this banner is touched by this task (or
# any later Phase W task): `resolve_interactions` and every pre-existing
# function/class in this module stay byte-for-byte identical (the
# append-only contract of
# docs/superpowers/plans/2026-08-16-bond-transverse-phase-w.md, "Global
# Constraints" -- `bare_bond_vertices`/`resolve_interactions` are templates,
# never edited).
#
# Implements the "master transverse topology" of
# docs/superpowers/specs/2026-08-15-bond-transverse-design.md ("The
# enlarged transverse objects (contracts)" -> "The master transverse
# topology", "Canonical orbit rule", "Construction domain (binding for
# steps 2-3)"): `TransverseTopology`, `resolve_transverse_topology`,
# `iter_reversal_orbits`, `validate_topology_against_mesh`,
# `transverse_effective_activity`. This section stays general (topology and
# coefficient bookkeeping only) -- vertex construction (`W_pm_bond`, Phase W
# Task 3) and the bubble (`transverse_bond_bubble_static`, Task 4) are NOT
# built here; they consume this section's objects.

import collections as _collections  # local: keeps the append-only region
                                     # self-contained without touching the
                                     # module's top-of-file import block.

_TRANSVERSE_ACTIVE_TYPES = ("CoulombInter", "Ising", "Exchange")

# Every interaction type the shared resolver can read a bond shell from
# (#181 Tier 3 Phase A). A caller selects a subset as ``active_types``:
# the transverse (spin-flip) channel reads _TRANSVERSE_ACTIVE_TYPES, the
# longitudinal (spin/charge) channel _LONGITUDINAL_ACTIVE_TYPES -- Hund is
# a longitudinal bond source with no transverse content, Exchange the
# reverse (its longitudinal content is adjudicated absent, Tier 2; its
# promotion through the bond basis is a recorded follow-up).
_BOND_RESOLVABLE_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange")
_LONGITUDINAL_ACTIVE_TYPES = ("CoulombInter", "Hund", "Ising")


def _as_readonly_int_array(value, name, *, ndim, last_dim=None):
    """Validate/convert ``value`` into a read-only ``int64`` array for
    :class:`TransverseTopology`'s ``delta_r``/``reverse`` fields.

    Rejects (``ValueError``, never silently clamps/truncates -- the same
    philosophy ``resolve_interactions``'s ``bond_max_shells`` validation
    uses): wrong ``ndim``, a wrong trailing axis length (when ``last_dim``
    is given), a complex dtype, a non-finite value, or a genuinely
    fractional value. A float array whose entries are all exact integers
    (e.g. ``np.array([[0, 0, 0]], dtype=float)``) is accepted and cast.
    """
    arr = np.asarray(value)
    if arr.ndim != ndim:
        raise ValueError(
            "BondTopology: {} must have ndim={}; got shape {}".format(
                name, ndim, arr.shape))
    if last_dim is not None and arr.shape[-1] != last_dim:
        raise ValueError(
            "BondTopology: {} must have last axis of length {}; got "
            "shape {}".format(name, last_dim, arr.shape))
    if np.iscomplexobj(arr):
        raise ValueError(
            "BondTopology: {} must be real/integer-valued; got a "
            "complex dtype ({})".format(name, arr.dtype))
    if arr.size and not np.all(np.isfinite(arr)):
        raise ValueError(
            "BondTopology: {} must be finite; got {}".format(name, arr))
    if not np.issubdtype(arr.dtype, np.integer):
        if arr.size and not np.all(arr == np.round(arr)):
            raise ValueError(
                "BondTopology: {} must be integer-valued (no "
                "fractional part); got {}".format(name, arr))
    out = np.array(arr, dtype=np.int64, copy=True)
    out.flags.writeable = False
    return out


_TransverseTopologyBase = _collections.namedtuple(
    "_TransverseTopologyBase", ["delta_r", "reverse", "coeffs"])


class BondTopology(_TransverseTopologyBase):
    """Alias-safe, closure-validated bond-channel topology for the
    transverse (spin-flip) channel (spec, "The master transverse
    topology").

    Fields
    ------
    delta_r : np.ndarray, shape (B, 3), int64, read-only
        Bravais-cell displacement per channel. ``delta_r[0] == (0, 0, 0)``
        ALWAYS -- the mandatory local channel every topology contains,
        never removed by shell truncation.
    reverse : np.ndarray, shape (B,), int64, read-only
        ``reverse[m]`` is the channel index of ``-delta_r[m]``;
        ``reverse[reverse[m]] == m`` for every ``m`` (an involution) and
        ``reverse[0] == 0``.
    coeffs : dict[str, np.ndarray]
        Per-interaction-type ``(B, norb, norb)`` complex arrays (a plain,
        NEW dict -- mutating the dict passed into the constructor does not
        affect the stored one), all sharing the same ``B``/``norb``. Every
        array's own arrays are read-only, ``coeffs[t][0]`` is EXACTLY zero
        (on-site content is represented exclusively by the caller's
        ``ham_pm_onsite``, never here), and the closure invariant
        ``coeffs[t][reverse[m]] == conj(coeffs[t][m]).T`` holds EXACTLY
        (to the last bit) for every ``m``.

    ``B = len(delta_r)`` (channel count, including channel 0); the enlarged
    dimension a consumer builds from this topology is ``ND = B * norb**2``.

    "Alias-safe", not "immutable": array VALUES cannot be mutated in place
    (``flags.writeable = False``), but ``coeffs`` is an ordinary dict a
    caller could re-key or replace wholesale on an aliased reference --
    every consumer boundary (:func:`validate_topology_against_mesh`)
    therefore re-validates rather than trusting a previously-constructed
    instance is still well-formed.

    Raises
    ------
    ValueError
        On any invariant violation, at construction time (never
        assert -- this is a public constructor, not an internal
        precondition).
    """

    __slots__ = ()

    def __new__(cls, delta_r, reverse, coeffs):
        delta_r = _as_readonly_int_array(delta_r, "delta_r", ndim=2, last_dim=3)
        B = delta_r.shape[0]
        if B < 1:
            raise ValueError(
                "BondTopology: delta_r must carry at least channel 0 "
                "(the mandatory local R=(0,0,0) channel); got shape {}"
                .format(delta_r.shape))
        if not np.array_equal(delta_r[0], np.zeros(3, dtype=np.int64)):
            raise ValueError(
                "BondTopology: channel 0 must be R=(0,0,0) (the "
                "mandatory local channel every topology contains); got "
                "delta_r[0]={}".format(tuple(int(x) for x in delta_r[0])))

        reverse = _as_readonly_int_array(reverse, "reverse", ndim=1)
        if reverse.shape != (B,):
            raise ValueError(
                "BondTopology: reverse must have shape ({},), "
                "matching delta_r's channel count B={}; got {}".format(
                    B, B, reverse.shape))
        if reverse.size and (np.any(reverse < 0) or np.any(reverse >= B)):
            raise ValueError(
                "BondTopology: reverse entries must all be valid "
                "channel indices in [0, {}); got {}".format(
                    B, reverse.tolist()))
        if not np.array_equal(reverse[reverse], np.arange(B, dtype=np.int64)):
            raise ValueError(
                "BondTopology: reverse must be an involution "
                "(reverse[reverse[m]] == m for every channel m); got "
                "reverse={}".format(reverse.tolist()))
        if int(reverse[0]) != 0:
            raise ValueError(
                "BondTopology: reverse[0] must be 0 -- channel 0 "
                "(R=(0,0,0)) is its own reversal partner; got reverse[0]={}"
                .format(int(reverse[0])))
        # Geometric consistency: reverse[m] must genuinely be the channel of
        # -delta_r[m]. This is not separately listed among the spec's
        # constructor-time checks, but reverse's whole MEANING depends on
        # it, and every downstream consumer (the closure check right below,
        # iter_reversal_orbits, W_pm_bond's orbit iteration) silently
        # produces wrong results if it is violated -- the same
        # internal-interface guard bare_bond_vertices applies to
        # ResolvedInteractionSet.reverse (see its "bond_set is not
        # reversal-closed" ValueError above).
        neg = -delta_r
        if not np.array_equal(delta_r[reverse], neg):
            mism = np.any(delta_r[reverse] != neg, axis=1)
            bad = int(np.argmax(mism))
            raise ValueError(
                "BondTopology: reverse is not geometrically "
                "consistent with delta_r -- delta_r[reverse[{0}]]={1} != "
                "-delta_r[{0}]={2}".format(
                    bad, tuple(int(x) for x in delta_r[reverse[bad]]),
                    tuple(int(x) for x in neg[bad])))

        new_coeffs = {}
        common_norb = None
        for type_name, raw_arr in dict(coeffs).items():
            arr = np.array(raw_arr, dtype=complex, copy=True)
            if arr.ndim != 3 or arr.shape[0] != B or arr.shape[1] != arr.shape[2]:
                raise ValueError(
                    "BondTopology: coeffs[{!r}] must have shape "
                    "(B, norb, norb) with B={} (matching delta_r); got {}"
                    .format(type_name, B, arr.shape))
            if common_norb is None:
                common_norb = arr.shape[1]
            elif arr.shape[1] != common_norb:
                raise ValueError(
                    "BondTopology: coeffs arrays disagree on norb -- "
                    "coeffs[{!r}] has norb={} but an earlier type has "
                    "norb={}".format(type_name, arr.shape[1], common_norb))
            if np.any(arr[0] != 0):
                raise ValueError(
                    "BondTopology: coeffs[{!r}][0] (channel 0, the "
                    "on-site R=(0,0,0) channel) must be EXACTLY zero -- "
                    "on-site content is represented exclusively by the "
                    "caller's ham_pm_onsite, never by a topology coeffs "
                    "array; got coeffs[{!r}][0]={}".format(
                        type_name, type_name, arr[0]))
            mism = arr[reverse] - np.conj(np.swapaxes(arr, 1, 2))
            if np.any(mism != 0):
                m_bad = int(np.argmax(np.any(mism != 0, axis=(1, 2))))
                raise ValueError(
                    "BondTopology: coeffs[{!r}] fails the Hermitian "
                    "closure invariant coeffs[reverse[m]] == "
                    "conj(coeffs[m]).T at channel m={} (reverse[m]={}): "
                    "coeffs[{!r}][{}]={} but conj(coeffs[{!r}][{}]).T={}"
                    .format(type_name, m_bad, int(reverse[m_bad]),
                            type_name, m_bad, arr[reverse[m_bad]],
                            type_name, m_bad, np.conj(arr[m_bad]).T))
            arr.flags.writeable = False
            new_coeffs[type_name] = arr

        return super(BondTopology, cls).__new__(
            cls, delta_r, reverse, new_coeffs)


# Phase W name, kept as an alias: the topology carrier is channel-agnostic
# (the CHANNEL is decided by which types the resolver is asked to read and
# by which vertex builder consumes the result; #181 Tier 3 Phase A).
TransverseTopology = BondTopology


def iter_reversal_orbits(topo):
    """Yield ONE representative ``(m, a, b)`` per reversal orbit
    ``{(m, a, b), (reverse[m], b, a)}`` of ``topo`` (spec, "Canonical orbit
    rule"): the representative is the LEXICOGRAPHIC TUPLE-MINIMUM of the
    orbit's (at most two) elements.

    This is the ONE shared helper both the topology layer (this module)
    and the vertex construction (``W_pm_bond``, Phase W Task 3) use to
    enumerate reversal orbits -- unit-tested standalone here, independent
    of either caller, per the spec's "Canonical orbit rule".

    Ranges over EVERY channel ``topo`` carries (``m`` in ``[0, B)``,
    INCLUDING the mandatory on-site channel 0) and every orbital pair
    ``(a, b)`` in ``[0, norb) x [0, norb)``; ``norb`` is read from
    ``topo.coeffs`` (every array there shares one ``(B, norb, norb)``
    shape by :class:`TransverseTopology`'s own constructor validation).
    Callers that only want the OFF-SITE domain (spec, "Construction domain
    (binding for steps 2-3)": steps 2-3 of ``W_pm_bond`` iterate ``R != 0``
    content exclusively) filter this function's output to ``m != 0``
    themselves -- this helper stays domain-agnostic, matching the general
    "Canonical orbit rule" statement.

    At ``m == reverse[m]`` (only possible at ``m == 0``, since a genuine
    ``m != 0`` channel wrapping onto its own reversal on a finite mesh is
    rejected by :func:`validate_topology_against_mesh`, never silently
    produced here) and ``a == b``, the orbit is a SINGLETON (its own
    partner); every other orbit has exactly two elements.

    Parameters
    ----------
    topo : TransverseTopology

    Returns
    -------
    list of (int, int, int)
        One ``(m, a, b)`` representative per orbit, in the deterministic
        order the ``(m, a, b)`` triple-loop first encounters each orbit
        (which is already ascending in the representative itself, since
        the representative is always whichever element of the pair is
        visited first).

    Raises
    ------
    ValueError
        If ``topo.coeffs`` is empty (``norb`` cannot be inferred -- every
        topology :func:`resolve_transverse_topology` returns carries the
        ``CoulombInter``/``Ising``/``Exchange`` keys always, even when
        zero-filled, so this is only reachable from a hand-built
        ``TransverseTopology`` with an empty ``coeffs`` dict).
    """
    if not topo.coeffs:
        raise ValueError(
            "iter_reversal_orbits: topo.coeffs is empty, so norb cannot be "
            "inferred (every TransverseTopology resolve_transverse_topology "
            "returns carries at least the CoulombInter/Ising/Exchange "
            "zero-filled arrays)")
    norb = None
    for type_name, arr in topo.coeffs.items():
        if norb is None:
            norb = arr.shape[1]
        elif arr.shape[1] != norb:
            raise ValueError(
                "iter_reversal_orbits: topo.coeffs arrays disagree on norb "
                "({} vs coeffs[{!r}]'s {})".format(norb, type_name, arr.shape[1]))

    reverse = topo.reverse
    B = len(reverse)
    seen = set()
    reps = []
    for m in range(B):
        mr = int(reverse[m])
        for a in range(norb):
            for b in range(norb):
                key = (m, a, b)
                if key in seen:
                    continue
                partner = (mr, b, a)
                seen.add(key)
                seen.add(partner)
                reps.append(min(key, partner))
    return reps


def validate_topology_against_mesh(topo, spatial_shape, *, arrays=None):
    """The ONE mesh-dependent validation lifecycle point for a
    :class:`TransverseTopology` (spec, "The master transverse topology"):
    invoked by BOTH ``transverse_bond_bubble_static`` (Phase W Task 4) and
    ``W_pm_bond`` (Task 3) at entry, so the injectivity rule below is
    enforced exactly once per call site.

    Validates, in order:

    1. ``spatial_shape``: a length-3 sequence of strictly positive
       integers (``Nx, Ny, Nz``).
    2. (optional) every array in ``arrays`` (a ``{label: ndarray}``
       mapping) has a leading axis equal to ``nvol = prod(spatial_shape)``
       -- the "any independently shaped numerical input's nvol must equal
       prod(spatial_shape)" check the spec's "master transverse topology"
       section requires this function to apply (e.g. a caller's Green
       function or ``ham_pm_onsite`` tensor built on the wrong mesh).
    3. ``topo``'s own invariants, by round-tripping it back through
       :class:`TransverseTopology`'s validating constructor (alias-safe,
       not immutable: ``coeffs`` -- though its ARRAYS are read-only -- is
       itself a plain, mutable dict a caller could have re-keyed since
       ``topo`` was built).
    4. INJECTIVITY of ``delta_r mod spatial_shape`` over ALL channels: any
       two DISTINCT channels whose displacement coincides after periodic
       wrapping (``R`` vs ``R + k*L`` on any axis, e.g. ``R=1`` vs ``R=6``
       on ``L=5``, INCLUDING the even-mesh self-reversal alias
       ``+L/2 == -L/2``) are REJECTED, naming both channels (index and raw
       ``delta_r``) and the mesh. A topology that passes this check can
       never produce duplicate ED operators under
       ``SectorED.bond_correlator_transverse``'s own post-wrap duplicate
       detection (spec).

    Parameters
    ----------
    topo : TransverseTopology
    spatial_shape : sequence of 3 ints
        ``(Nx, Ny, Nz)``.
    arrays : dict[str, array_like], optional
        Additional numerical inputs whose leading axis must equal
        ``nvol = Nx*Ny*Nz``.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        On any of the above.
    """
    try:
        shape_seq = tuple(spatial_shape)
    except TypeError:
        raise ValueError(
            "validate_topology_against_mesh: spatial_shape must be a "
            "length-3 sequence of strictly positive integers; got {!r}"
            .format(spatial_shape))
    if len(shape_seq) != 3:
        raise ValueError(
            "validate_topology_against_mesh: spatial_shape must have "
            "length 3 (Nx, Ny, Nz); got {!r}".format(spatial_shape))
    shape_i = []
    for x in shape_seq:
        try:
            xf = float(x)
        except (TypeError, ValueError):
            raise ValueError(
                "validate_topology_against_mesh: spatial_shape entries "
                "must be integers; got {!r}".format(spatial_shape))
        if not np.isfinite(xf) or xf != np.floor(xf) or xf <= 0:
            raise ValueError(
                "validate_topology_against_mesh: spatial_shape entries "
                "must be strictly positive integers; got {!r}".format(
                    spatial_shape))
        shape_i.append(int(xf))
    Nx, Ny, Nz = shape_i
    nvol = Nx * Ny * Nz

    if arrays:
        for label, value in dict(arrays).items():
            value = np.asarray(value)
            got = value.shape[0] if value.ndim else None
            if value.ndim < 1 or got != nvol:
                raise ValueError(
                    "validate_topology_against_mesh: {} has leading axis "
                    "{} but the mesh {} implies nvol={}".format(
                        label, got, shape_i, nvol))

    # Re-validate topo (alias-safe, not immutable -- see the class
    # docstring); a hand-built or corrupted topology is caught here rather
    # than silently producing wrong results downstream.
    topo = TransverseTopology(topo.delta_r, topo.reverse, topo.coeffs)

    shape_arr = np.array(shape_i, dtype=np.int64)
    delta_r = np.asarray(topo.delta_r)
    B = delta_r.shape[0]
    wrapped = np.mod(delta_r, shape_arr)
    seen = {}
    for m in range(B):
        key = tuple(int(x) for x in wrapped[m])
        if key in seen:
            other = seen[key]
            raise ValueError(
                "validate_topology_against_mesh: channels {} (delta_r={}) "
                "and {} (delta_r={}) collide on the mesh {} -- both wrap to "
                "the same displacement {} modulo the mesh (an alias such as "
                "R vs R+k*L on some axis, or the even-mesh self-reversal "
                "+L/2 == -L/2); a topology whose channels are not "
                "INJECTIVE modulo the mesh cannot be used at this mesh "
                "resolution.".format(
                    other, tuple(int(v) for v in delta_r[other]),
                    m, tuple(int(v) for v in delta_r[m]), shape_i, key))
        seen[key] = m
    return None


def bond_effective_activity(topo):
    """The shared "is this topology doing anything transverse off-site"
    predicate: EXACT nonzero after summation, Hermitian projection and
    shell truncation -- i.e. after everything
    :func:`resolve_transverse_topology` already did to ``topo.coeffs``.

    A topology built entirely from PASS-ZERO declarations (only
    Hund/PairLift, or no off-site transverse-active interaction at all)
    has every ``coeffs`` array identically zero (channel 0 always is, by
    :class:`TransverseTopology`'s own constructor invariant; channel 0 is
    the ONLY channel when nothing off-site was ever declared) and this
    returns ``False``; any genuinely nonzero off-site CoulombInter/
    Ising/Exchange coefficient makes it ``True``.

    This is the ONE shared place that answers the question -- callers
    (``W_pm_bond``, ``transverse_bond_bubble_static``, and any future
    gate/preflight code, Phase W Tasks 3+) use this rather than
    re-deriving their own "is coeffs all zero" check.

    Parameters
    ----------
    topo : TransverseTopology

    Returns
    -------
    bool
    """
    for arr in topo.coeffs.values():
        if np.any(arr != 0):
            return True
    return False


# Phase W name (see BondTopology): one predicate for both channels.
transverse_effective_activity = bond_effective_activity


class BondSetView:
    """A :class:`BondTopology` presented under the bubble kernel's
    duck-typed bond-set contract (``bubble._validate_bond_set``:
    ``n_channels`` a positive int, ``delta_r`` a sequence of integer
    3-tuples, plus ``reverse``): the longitudinal bond bubble
    (``bubble.bond_bubble_static``) and the test-side ED map consume a
    topology through this view instead of a ``ResolvedInteractionSet``
    (#181 Tier 3 Phase A). The topology is re-validated on construction
    (alias-safe, like ``validate_topology_against_mesh``)."""

    def __init__(self, topo):
        topo = BondTopology(delta_r=topo.delta_r, reverse=topo.reverse,
                            coeffs=topo.coeffs)
        self.delta_r = tuple(tuple(int(x) for x in row)
                             for row in np.asarray(topo.delta_r))
        self.reverse = tuple(int(x) for x in np.asarray(topo.reverse))
        self.n_channels = len(self.delta_r)



def _close_offsite_hermitian(by_irvec, norb, reverse_atol, reverse_rtol,
                              type_name, diag="resolve_transverse_topology"):
    """Off-site (``R != 0``) reverse-closure / Hermitian-projection for
    ONE interaction type, per the SAME algorithm ``resolve_interactions``
    (its "Step 2: reverse closure / Hermiticity synthesis", above in this
    module) applies to its ``CoulombInter``-shaped input.

    Deliberately RE-IMPLEMENTED here rather than shared with
    ``resolve_interactions`` -- this task's append-only contract forbids
    editing that function (or extracting a common helper out of it) -- but
    the algorithm itself is the SAME normalization
    :func:`resolve_transverse_topology` must reproduce (carried Task-1
    review finding): for a density-density (commuting-operator) type,
    declaring BOTH ``C`` at ``R`` and ``conj(C)`` at ``-R`` describes the
    SAME physical operator sum TWICE (``bare_bond_vertices``'s Finding-1
    derivation, above, reproduces the argument in full) -- so a
    consistent mirrored pair is Hermitian-PROJECTED onto its exact
    conjugate average (``0.5 * (declared + conj(synthesized))``), never
    stored as either raw declared value, and an INCONSISTENT mirrored pair
    (beyond the scale-aware tolerance) raises rather than silently picking
    one side. A single-direction declaration is stored as declared, with
    its reverse partner synthesized as the exact conjugate. Duplicate
    declarations at the identical ``(type, R, a, b)`` key cannot occur
    here (a plain Python dict already deduplicates that -- same note
    ``resolve_interactions``'s own "Step 1" docstring makes); a caller
    that pre-sums two logical contributions into one dict entry gets that
    pre-summed value carried through this same closure unchanged.

    Parameters
    ----------
    by_irvec : dict[tuple, dict[(int, int), complex]]
        ``{irvec: {(a, b): value}}``, OFF-SITE (``irvec != (0, 0, 0)``)
        entries only (the caller filters out the on-site point).
    norb : int
    reverse_atol, reverse_rtol : float
        Same scale-aware tolerance semantics as ``resolve_interactions``.
    type_name : str
        Only used to make a ``ValueError`` message actionable.

    Returns
    -------
    dict[tuple, np.ndarray]
        ``{irvec: (norb, norb) complex}``, closed under
        ``irvec -> -irvec`` (every key's negation is also a key).
    """
    declared_irvecs = set(by_irvec.keys())
    all_irvecs = set(declared_irvecs)
    for irvec in declared_irvecs:
        all_irvecs.add(tuple(-x for x in irvec))

    completed = {}
    visited = set()
    for irvec in all_irvecs:
        if irvec in visited:
            continue
        neg = tuple(-x for x in irvec)
        visited.add(irvec)
        visited.add(neg)

        fwd = by_irvec.get(irvec, {})
        bwd = by_irvec.get(neg, {})

        fwd_mat = np.zeros((norb, norb), dtype=complex)
        bwd_mat = np.zeros((norb, norb), dtype=complex)

        synth_fwd_from_bwd = {(a, b): np.conj(v) for (b, a), v in bwd.items()}
        synth_bwd_from_fwd = {(b, a): np.conj(v) for (a, b), v in fwd.items()}

        keys_fwd = set(fwd.keys()) | set(synth_fwd_from_bwd.keys())
        for (a, b) in keys_fwd:
            v_declared = fwd.get((a, b))
            v_synth = synth_fwd_from_bwd.get((a, b))
            if v_declared is not None and v_synth is not None:
                tol = reverse_atol + reverse_rtol * abs(v_declared)
                if abs(v_declared - v_synth) > tol:
                    raise ValueError(
                        diag + ": {}_{}{}({}) = {} "
                        "disagrees with conj({}_{}{}({})) = {} beyond "
                        "tolerance (atol={}, rtol={})".format(
                            type_name, a, b, irvec, v_declared,
                            type_name, b, a, neg, v_synth,
                            reverse_atol, reverse_rtol))
                fwd_mat[a, b] = 0.5 * (v_declared + v_synth)
            elif v_declared is not None:
                fwd_mat[a, b] = v_declared
            else:
                fwd_mat[a, b] = v_synth

        keys_bwd = set(bwd.keys()) | set(synth_bwd_from_fwd.keys())
        for (a, b) in keys_bwd:
            v_declared = bwd.get((a, b))
            v_synth = synth_bwd_from_fwd.get((a, b))
            if v_declared is not None and v_synth is not None:
                tol = reverse_atol + reverse_rtol * abs(v_declared)
                if abs(v_declared - v_synth) > tol:
                    raise ValueError(
                        diag + ": {}_{}{}({}) = {} "
                        "disagrees with conj({}_{}{}({})) = {} beyond "
                        "tolerance (atol={}, rtol={})".format(
                            type_name, a, b, neg, v_declared,
                            type_name, b, a, irvec, v_synth,
                            reverse_atol, reverse_rtol))
                bwd_mat[a, b] = 0.5 * (v_declared + v_synth)
            elif v_declared is not None:
                bwd_mat[a, b] = v_declared
            else:
                bwd_mat[a, b] = v_synth

        completed[irvec] = fwd_mat
        if neg != irvec:
            completed[neg] = bwd_mat
        # neg == irvec is impossible here: irvec != (0, 0, 0) always (the
        # caller filters out the on-site point before calling), and
        # -irvec == irvec only at irvec == (0, 0, 0).

    return completed


def _require_transverse_integral(value, label, diag="resolve_transverse_topology"):
    """Coerce ``value`` to ``int`` iff it is an exact integral scalar;
    otherwise raise a ``resolve_transverse_topology``-prefixed
    ``ValueError`` naming ``label``. ``bool`` is rejected even though
    ``isinstance(True, int)`` in Python, since a stray boolean (e.g. an
    ``irvec``/``orbvec`` component that is accidentally a truthiness
    flag) is never a legitimate integral coordinate here.
    """
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(
            diag + ": {} must be an integer, not a "
            "bool, got {!r}".format(label, value))
    # Exact-integer types pass through unchanged (arbitrary precision:
    # no float round-trip, which would silently round values above
    # 2**53 and overflow on huge ints).
    if isinstance(value, (int, np.integer)):
        return int(value)
    # Float types are accepted iff they carry an exact integral value
    # small enough that the value is unambiguous (|v| <= 2**53, the
    # float64 contiguous-integer bound). Everything else -- strings,
    # arrays, None, complex, non-integral or non-finite floats -- is
    # rejected; numeric strings deliberately do NOT coerce.
    if isinstance(value, (float, np.floating)):
        # Validate at the value's NATIVE precision (no float() narrowing:
        # a fractional np.longdouble could round to an integral float64,
        # and a longdouble just above 2**53 could narrow onto the
        # accepted boundary).
        if (not np.isfinite(value) or value != np.floor(value)
                or abs(value) > 2.0**53):
            raise ValueError(
                diag + ": {} must be an integral "
                "value, got {!r}".format(label, value))
        return int(value)
    raise ValueError(
        diag + ": {} must be an integer, got "
        "{!r}".format(label, value))


def _resolve_bond_topology_impl(interactions, cell, norb, *, max_shells,
                                active_types, diag):
    """Shared body of :func:`resolve_bond_topology` and
    :func:`resolve_transverse_topology`; ``diag`` prefixes every
    diagnostic so each public entry reports under its own name."""
    if (isinstance(active_types, str)
            or not isinstance(active_types, (tuple, list))
            or len(active_types) == 0
            or any(not isinstance(t, str) for t in active_types)
            or len(set(active_types)) != len(active_types)
            or any(t not in _BOND_RESOLVABLE_TYPES for t in active_types)):
        raise ValueError(
            "{}: active_types must be a non-empty, duplicate-free sequence "
            "drawn from {}; got {!r} (an aggregate 'Coulomb' table must be "
            "split by the caller first)".format(
                diag, _BOND_RESOLVABLE_TYPES, active_types))
    active_types = tuple(active_types)

    norb = _require_transverse_integral(norb, "norb", diag)
    if norb < 1:
        raise ValueError(
            diag + ": norb must be a positive integer, "
            "got {!r}".format(norb))
    if not hasattr(interactions, "get"):
        raise ValueError(
            diag + ": interactions must be a "
            "dict-like container (mapping interaction type name -> "
            "{(irvec, orbvec): value} sub-dict), got {!r}".format(
                type(interactions)))

    reverse_atol, reverse_rtol = 1e-10, 1e-8

    per_type_by_irvec = {}
    per_type_completed = {}
    union_irvecs = set()
    for type_name in active_types:
        tbl = interactions.get(type_name, {}) or {}
        by_irvec = {}
        for (irvec, orbvec), value in tbl.items():
            try:
                irvec_len = len(irvec)
            except TypeError:
                raise ValueError(
                    diag + ": {} declares an irvec "
                    "key that is not a 3-component vector, got "
                    "{!r}".format(type_name, irvec))
            if irvec_len != 3:
                raise ValueError(
                    diag + ": {} declares an irvec "
                    "key with {} component(s), expected exactly 3, got "
                    "{!r}".format(type_name, irvec_len, irvec))
            irvec = tuple(
                _require_transverse_integral(x, "{} irvec component".format(type_name), diag)
                for x in irvec)

            try:
                orbvec_len = len(orbvec)
            except TypeError:
                raise ValueError(
                    diag + ": {} declares an orbvec "
                    "key that is not a 2-component (a, b) pair, got "
                    "{!r}".format(type_name, orbvec))
            if orbvec_len != 2:
                raise ValueError(
                    diag + ": {} declares an orbvec "
                    "key with {} component(s), expected exactly 2 (a, "
                    "b), got {!r}".format(type_name, orbvec_len, orbvec))
            a = _require_transverse_integral(orbvec[0], "{} orbvec[0]".format(type_name), diag)
            b = _require_transverse_integral(orbvec[1], "{} orbvec[1]".format(type_name), diag)
            if not (0 <= a < norb and 0 <= b < norb):
                raise ValueError(
                    diag + ": {} declares orbital "
                    "index (a={}, b={}) at irvec={} out of range for "
                    "norb={}; orbital indices must satisfy 0 <= a,b < "
                    "norb".format(type_name, a, b, irvec, norb))

            if irvec == (0, 0, 0):
                continue
            if (isinstance(value, (bool, np.bool_))
                    or not isinstance(value, numbers.Number)
                    or not np.isfinite(value)):
                raise ValueError(
                    diag + ": {} declares a coefficient that is not a "
                    "finite number at irvec {}, orbvec {}: got {!r}".format(
                        type_name, irvec, orbvec, value))
            by_irvec.setdefault(irvec, {})[(a, b)] = complex(value)
        per_type_by_irvec[type_name] = by_irvec
        completed = _close_offsite_hermitian(
            by_irvec, norb, reverse_atol, reverse_rtol, type_name, diag)
        per_type_completed[type_name] = completed
        union_irvecs.update(completed.keys())

    # --- shell sort: SAME ordering resolve_interactions uses (m=0 first,
    # then shells by (Euclidean length, lexicographic irvec)) ---
    def sort_key(irvec):
        return (_shell_length(irvec, cell), irvec)

    ordered_irvecs_all = sorted(union_irvecs, key=sort_key)
    shells = []  # list of (length, [irvecs...])
    for irvec in ordered_irvecs_all:
        length = _shell_length(irvec, cell)
        if shells and abs(length - shells[-1][0]) <= 1e-6 * max(
                length, shells[-1][0], 1.0):
            shells[-1][1].append(irvec)
        else:
            shells.append((length, [irvec]))

    if max_shells is not None:
        try:
            shells_f = float(max_shells)
        except (TypeError, ValueError):
            raise ValueError(
                diag + ": max_shells must be a "
                "non-negative integer or None, got {!r}".format(max_shells))
        if (not np.isfinite(shells_f) or shells_f < 0
                or shells_f != np.floor(shells_f)):
            raise ValueError(
                diag + ": max_shells must be a "
                "non-negative integral value (shell 0 is the mandatory "
                "on-site R=(0,0,0) channel, always kept regardless; "
                "max_shells counts OFF-SITE shells beyond it) or None, "
                "got {!r}".format(max_shells))
        n_keep = int(shells_f)
        dropped_shells = shells[n_keep:]

        # Ambiguity guard (review fix, generalizing resolve_interactions's
        # "Step 4" bond_max_shells=0 guard, bond_channels.py above, to
        # every truncation depth): dropping a shell that carries any
        # exactly-nonzero DECLARED coefficient (from CoulombInter, Ising
        # or Exchange) is ALWAYS an error, never silent -- the docstring's
        # "max_shells truncates the UNION of off-site shells" promise
        # would otherwise silently discard content the caller explicitly
        # asked for. Checked against per_type_by_irvec (the RAW declared
        # entries, pre-closure), matching resolve_interactions's own
        # has_nonzero_inter_site_v -- a synthesized-but-negligible
        # reverse partner never triggers this on its own (it is not
        # itself a declared entry), but its declared mirror at the SAME
        # shell (same |R|, since |R| == |-R| always) does.
        dropped_irvecs = {irvec for _, irvecs in dropped_shells
                           for irvec in irvecs}
        if dropped_irvecs:
            offending = []
            for type_name in active_types:
                by_irvec = per_type_by_irvec[type_name]
                for irvec in dropped_irvecs:
                    entries = by_irvec.get(irvec)
                    if entries and any(v != 0 for v in entries.values()):
                        offending.append((type_name, irvec))
            if offending:
                offending.sort(key=lambda x: (x[0], x[1]))
                raise ValueError(
                    diag + ": max_shells={} requested "
                    "but this would drop shell(s) carrying declared "
                    "nonzero off-site content: {}; this is ambiguous "
                    "(asks to truncate away declared coefficients rather "
                    "than merely unresolved/zero-valued shells). Use a "
                    "larger max_shells, or omit the offending "
                    "declarations, for the genuinely truncated model "
                    "instead.".format(max_shells, offending))
        shells = shells[:n_keep]

    ordered_irvecs = [irvec for _, irvecs in shells for irvec in irvecs]

    delta_r = [(0, 0, 0)] + ordered_irvecs
    index_of = {dr: i for i, dr in enumerate(delta_r)}
    B = len(delta_r)

    reverse = [0]
    for irvec in ordered_irvecs:
        neg = tuple(-x for x in irvec)
        if neg not in index_of:
            raise ValueError(
                diag + ": internal error -- channel "
                "{} has no reverse partner {} in the resolved union "
                "(reversal closure invariant violated)".format(irvec, neg))
        reverse.append(index_of[neg])

    coeffs = {}
    for type_name in active_types:
        completed = per_type_completed[type_name]
        arr = np.zeros((B, norb, norb), dtype=complex)
        for m, irvec in enumerate(delta_r):
            if m == 0:
                continue
            block = completed.get(irvec)
            if block is not None:
                arr[m] = block
        coeffs[type_name] = arr

    return BondTopology(delta_r=delta_r, reverse=reverse, coeffs=coeffs)


def resolve_bond_topology(interactions, cell, norb, *, max_shells=None,
                          active_types):
    """Resolve the master bond topology over an explicit ``active_types``
    selection (#181 Tier 3 Phase A): the union of every declared off-site
    shell of those types, channel 0 always ``R=(0,0,0)``, shells ordered
    by length then lexicographically, reversal-closed, ``coeffs`` holding
    EXACTLY ``active_types`` (in that order, zero-filled for a type that
    declares nothing off-site), each ``(B, norb, norb)`` complex and
    Hermitian-closed. Same pre-fold-declaration rule, mirrored-pair
    consistency check and ``max_shells`` truncation guard as
    :func:`resolve_transverse_topology` (whose docstring has the details);
    that function is this one at ``active_types=_TRANSVERSE_ACTIVE_TYPES``.
    Diagnostics are prefixed ``resolve_bond_topology:``.
    """
    return _resolve_bond_topology_impl(
        interactions, cell, norb, max_shells=max_shells,
        active_types=active_types, diag="resolve_bond_topology")


def resolve_transverse_topology(interactions, cell, norb, *, max_shells=None):
    """Resolve the master transverse bond topology (spec, "The master
    transverse topology") from the SAME raw, PRE-FOLD interaction
    declarations ``resolve_interactions`` consumes -- the pre-fold
    locality rule (judge locality on the ORIGINAL declarations, never a
    folded table; ``rpa.py``'s own ``param_ham_orig``-reading builders,
    e.g. ``_append_inter_cross``/``_append_onsite_direct`` above this
    module in the pipeline, apply the identical rule for the same reason).
    Callers pass ``interactions = self.param_ham_orig`` (or the equivalent
    un-folded dict-of-dicts) -- NEVER the post-sublattice-folding
    ``self.param_ham``.

    Collects the UNION of declared off-site displacement vectors from
    EVERY transverse-active type (``CoulombInter``, ``Ising``,
    ``Exchange``; ``Hund``/``PairLift`` and anything else in
    ``interactions`` contribute NO shells and are simply never read --
    "tolerated" per the spec), closes each type's own declarations under
    ``R -> -R`` with the SAME sum-then-Hermitian-project normalization
    ``resolve_interactions`` uses (see :func:`_close_offsite_hermitian`),
    orders channels the SAME way ``resolve_interactions`` does (``m=0``
    first, then shells by ``(Euclidean length, lexicographic irvec)``),
    and returns a :class:`TransverseTopology` whose ``coeffs`` are
    zero-filled, per type, at every channel that type does not declare.

    Parameters
    ----------
    interactions : dict
        A dict-of-dicts container in the SAME shape as
        ``self.param_ham_orig``/``self.param_ham``: keys are interaction
        type names; each sub-dict has the SAME ``{(irvec, orbvec): value}``
        shape ``resolve_interactions``'s own ``coulomb_inter`` parameter
        documents (``irvec`` a length-3 int tuple, ``orbvec = (a, b)``
        0-based orbital indices, ``value`` complex). A missing type key is
        treated exactly like an empty sub-dict; on-site (``irvec ==
        (0, 0, 0)``) entries of every type are ignored here -- on-site
        content is represented exclusively through the caller's
        ``ham_pm_onsite`` (Phase W Task 3's ``W_pm_bond``, spec step 1),
        never through a topology ``coeffs`` channel.
    cell : array_like, shape (3, 3)
        Real-space lattice vectors (rows) -- the SAME quantity
        ``resolve_interactions`` calls ``lattice_vectors``; used
        identically, only for Euclidean shell length / ``max_shells``
        ordering.
    norb : int
        Number of orbitals; every ``coeffs`` array is
        ``(B, norb, norb)``.
    max_shells : int or None, optional
        Truncates the UNION of off-site shells (the ``bond_max_shells``
        semantics ``resolve_interactions`` documents) by WHOLE ``|R|``-
        shells; ``None`` (the default) keeps every declared off-site
        shell. A negative, non-integral or non-finite value raises
        ``ValueError`` (never clamped/truncated). Dropping declared
        content is ALWAYS an error, never silent: if the truncation would
        remove a shell that carries any exactly-nonzero DECLARED
        coefficient (from ``CoulombInter``, ``Ising`` or ``Exchange``),
        this raises ``ValueError`` -- the same ambiguity guard
        ``resolve_interactions``'s ``bond_max_shells=0`` case applies
        (bond_channels.py's "Step 4" above), generalized here to any
        truncation depth: ``max_shells=0`` with declared off-site content
        is simply the ``n_keep=0`` instance of the same rule. A
        truncation that only drops shells with no declared nonzero
        coefficient (e.g. purely synthesized-negligible partners, or
        shells nothing ever declared) remains fine.

    Returns
    -------
    TransverseTopology
        Channel 0 is always ``R=(0,0,0)``; ``coeffs`` holds
        ``"CoulombInter"``, ``"Ising"`` and ``"Exchange"`` keys ALWAYS
        (zero-filled at every channel for a type that declares nothing at
        all off-site), each ``(B, norb, norb)`` complex, Hermitian-closed
        exactly like ``resolve_interactions``'s ``v_bond``.

    Raises
    ------
    ValueError
        If a mirrored declared pair disagrees beyond tolerance, if
        ``max_shells`` is invalid, if ``max_shells`` would drop a shell
        carrying declared nonzero content, or (propagated from
        :class:`TransverseTopology`'s constructor) if the resolved
        topology somehow fails its own invariants.
    """
    return _resolve_bond_topology_impl(
        interactions, cell, norb, max_shells=max_shells,
        active_types=_TRANSVERSE_ACTIVE_TYPES,
        diag="resolve_transverse_topology")


# =============================================================================
# Phase W, Task 3 -- the production bond-resolved transverse vertex W_pm_bond
# =============================================================================
#
# APPEND-ONLY REGION (continued): nothing above this banner is touched by
# this task -- `resolve_interactions`/`bare_bond_vertices` (templates) and
# every Task-2 name (`TransverseTopology`, `resolve_transverse_topology`,
# `iter_reversal_orbits`, `validate_topology_against_mesh`,
# `transverse_effective_activity`) stay byte-for-byte identical.
#
# Implements docs/superpowers/specs/2026-08-15-bond-transverse-design.md,
# "The vertex -- element equations" (steps 1-3, the AMENDED 2026-08-16 flip-
# family assignment) -- the ORDERED-RECORD equations Gate W0
# (tests/test_bond_transverse_ed.py, Task 1) numerically validated against
# exact diagonalization BEFORE this function existed. This function is a
# faithful production transcription of that gate's test-local reference,
# `w_expected_from_records`: same representative-orbit iteration (via the
# shared `iter_reversal_orbits`, filtered to the off-site R != 0 domain),
# same coefficient table, same q-phase convention -- cross-pinned against
# that reference at rel=0/abs=1e-13 on all three W0 ED granule fixtures in
# tests/test_bond_transverse_w.py.

def W_pm_bond(topo, ham_pm_onsite, *, spatial_shape):
    """The bond-resolved transverse (spin-flip) vertex ``W_{+-}(q)``, built
    from ``topo`` (a :class:`TransverseTopology`, spec "The master
    transverse topology") and the CALLER-ASSEMBLED on-site vertex
    ``ham_pm_onsite`` (spec "The vertex -- element equations", steps 1-3).

    Construction (spec, verbatim; also see the module-level Task 1 gate
    ``tests/test_bond_transverse_ed.py``'s ``w_expected_from_records``,
    which encodes the IDENTICAL equations independently and is this
    function's ED-validated oracle):

    1. **Channel-0 (on-site) block** = ``ham_pm_onsite`` verbatim, placed
       (never re-assembled -- ``W_pm_bond`` performs NO assembly of its
       own) into ``W[:, 0:nd, 0:nd]``. ``ham_pm_onsite`` already carries
       the mesh's ``nvol`` leading axis (the shape
       ``_assemble_transverse_vertex`` returns): for genuinely on-site-only
       declarations that tensor is q-INDEPENDENT (constant along the
       leading axis) by construction elsewhere
       (``_check_transverse_representable``), so this step is a placement,
       not a broadcast, despite the shared leading axis -- the B=1
       reduction test pins this: with ``B=1`` (on-site-only topology),
       ``W_pm_bond`` IS ``ham_pm_onsite`` reshaped, bit for bit.
    2. **Cross family** (``CoulombInter`` -> ``-Re(.)``, ``Ising`` ->
       ``+Re(.)``), OFF-SITE (``R != 0``) reversal-orbit representatives
       only (``iter_reversal_orbits(topo)``, filtered to ``m != 0`` --
       that helper's own docstring: "the mandatory local channel 0" is
       included in its general enumeration; the off-site filter belongs
       to every CALLER per the spec's "Construction domain"). Each
       representative ``(m, a, b)`` with coefficient
       ``C = topo.coeffs[type][m][a, b]`` emits the SAME real value into
       BOTH mirrored diagonal target cells, DIRECT placement (AMENDED
       2026-08-16, Task-6 granule adjudication -- the channel carrying
       the declared coefficient is the target, matching the longitudinal
       ``bare_bond_vertices`` Fock-diagonal rule exactly; the PREVIOUS
       ``m <-> reverse[m]`` swap FAILED the multi-orbital off-site
       CoulombInter granule at O(1), invisible at ``a == b`` which is why
       Gate W0's norb=1 fixtures never caught it -- see the spec's "Cross
       family" section for the full derivation):
       ``W[q, (m, a, b), (m, a, b)] += s_type * Re(C)``
       and ``W[q, (reverse[m], b, a), (reverse[m], b, a)] += s_type * Re(C)``
       -- bond-diagonal, q-INDEPENDENT. The per-slot ``Re(.)`` placement
       (never the raw complex value) is deliberate: a Hermitian-closed
       ``+-i*eps`` null-direction pair must leave the diagonal exactly
       real (the structural Hermiticity test with a complex
       ``CoulombInter`` coefficient pins this).
    3. **Flip family** (``Exchange``, ``f_J = -1``), the SAME off-site
       representatives. Each representative ``(m, a, b)`` with
       ``J = topo.coeffs["Exchange"][m][a, b]`` at ``R = topo.delta_r[m]``
       emits the AMENDED (2026-08-16, Gate W0 adjudication) two ordered
       records, q-dependent ONLY here:
       ``W[q, (0,a,a), (0,b,b)] += f_J * conj(J) * exp(-i q.R)`` and
       ``W[q, (0,b,b), (0,a,a)] += f_J * J * exp(+i q.R)``. The PRE-
       amendment draft form (``f_J*J*exp(+iqR)`` at ``(aa,bb)``) failed
       Gate W0's multi-orbital off-site Exchange granule systematically
       (residuals 0.18-0.33); this is the swapped, ED-adjudicated
       assignment.

    The q mesh is ``q = 2*pi*(n_x/N_x, n_y/N_y, n_z/N_z)`` (``spatial_shape
    = (N_x, N_y, N_z)``), C-order flattened -- the SAME convention
    ``sc._build_bond_m0_blocks``'s own phase composition uses (cross-pinned
    by a dedicated test on a shared fixture); ``q.R = 2*pi*sum_d
    n_d*R_d/N_d``. The ``(Nx, Ny, Nz, ...)`` frame used to build this mesh
    is private to this function; the returned array is the canonical
    flattened ``(nvol, ND, ND)`` shape every numerical object in the spec
    uses.

    Parameters
    ----------
    topo : TransverseTopology
        The master transverse bond topology (``resolve_transverse_topology``
        or an equivalent hand-built, invariant-satisfying instance).
    ham_pm_onsite : array_like, shape (nvol, norb, norb, norb, norb),
        complex
        The ALREADY-ASSEMBLED on-site transverse vertex (e.g.
        ``RPA._assemble_transverse_vertex`` applied to the on-site-only
        tensor filtered from the ORIGINAL PRE-FOLD declarations at
        ``irvec == (0, 0, 0)``) -- never re-derived here.
    spatial_shape : sequence of 3 ints
        ``(Nx, Ny, Nz)``; ``nvol = Nx*Ny*Nz`` must match both
        ``ham_pm_onsite``'s leading axis and ``topo``'s own mesh-
        injectivity requirement (``validate_topology_against_mesh``,
        invoked here at entry).

    Returns
    -------
    ndarray, complex128, shape (nvol, ND, ND)
        ``ND = B * norb**2``, ``B = len(topo.delta_r)``.

    Raises
    ------
    ValueError
        Propagated from :func:`validate_topology_against_mesh` (mesh-
        injectivity violation, malformed ``spatial_shape``, a
        ``ham_pm_onsite`` whose leading axis disagrees with ``nvol``, or a
        stale/corrupted ``topo``), or raised directly here if
        ``ham_pm_onsite`` has the wrong ndim/shape or its orbital count
        disagrees with ``topo.coeffs``'s.
    """
    try:
        Nx, Ny, Nz = (int(x) for x in spatial_shape)
    except (TypeError, ValueError):
        raise ValueError(
            "W_pm_bond: spatial_shape must be a length-3 sequence of "
            "integers; got {!r}".format(spatial_shape))
    nvol = Nx * Ny * Nz

    ham_pm_onsite = np.array(ham_pm_onsite, dtype=complex, copy=False)
    if ham_pm_onsite.ndim != 5:
        raise ValueError(
            "W_pm_bond: ham_pm_onsite must have ndim=5 (nvol, norb, norb, "
            "norb, norb); got shape {}".format(ham_pm_onsite.shape))
    norb = ham_pm_onsite.shape[1]
    expected_onsite_shape = (nvol, norb, norb, norb, norb)
    if ham_pm_onsite.shape != expected_onsite_shape:
        raise ValueError(
            "W_pm_bond: ham_pm_onsite must have shape (nvol, norb, norb, "
            "norb, norb) = {} (nvol from spatial_shape={}, norb inferred "
            "from ham_pm_onsite.shape[1]); got {}".format(
                expected_onsite_shape, spatial_shape, ham_pm_onsite.shape))

    # The ONE mesh-dependent validation lifecycle point (spec, "The master
    # transverse topology"): mesh-injectivity of topo.delta_r AND that
    # ham_pm_onsite's leading axis matches nvol.
    validate_topology_against_mesh(
        topo, spatial_shape, arrays={"ham_pm_onsite": ham_pm_onsite})
    # Re-fetch a freshly-validated, alias-safe topology for local use (the
    # same defensive re-construction validate_topology_against_mesh applies
    # internally -- topo.coeffs is a plain dict a caller could have re-keyed
    # between that call returning and this line).
    topo = TransverseTopology(topo.delta_r, topo.reverse, topo.coeffs)

    delta_r = np.asarray(topo.delta_r)
    reverse = np.asarray(topo.reverse)
    B = delta_r.shape[0]
    coeffs = topo.coeffs
    for type_name, arr in coeffs.items():
        if arr.shape[1] != norb:
            raise ValueError(
                "W_pm_bond: topo.coeffs[{!r}] has norb={} but "
                "ham_pm_onsite implies norb={}".format(
                    type_name, arr.shape[1], norb))

    nd = norb * norb
    ND = B * nd
    W = np.zeros((nvol, ND, ND), dtype=complex)

    # Step 1: channel-0 block = the received on-site vertex, placed
    # verbatim (no assembly, no broadcast beyond what ham_pm_onsite's own
    # leading axis already carries).
    W[:, 0:nd, 0:nd] = ham_pm_onsite.reshape(nvol, nd, nd)

    # q mesh (private (Nx, Ny, Nz) frame, C-order flattened), spec's fixed
    # convention -- the same one sc._build_bond_m0_blocks' phase uses.
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    kz = 2.0 * np.pi * np.arange(Nz) / Nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    q_flat = np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=-1)

    # Off-site (R != 0) reversal-orbit representatives -- the shared helper
    # (spec, "Canonical orbit rule") ranges over EVERY channel including 0;
    # steps 2-3 are binding to the off-site construction domain only (spec,
    # "Construction domain (binding for steps 2-3)"), so channel 0 is
    # filtered out here, not inside the shared helper.
    reps = [rep for rep in iter_reversal_orbits(topo) if rep[0] != 0]

    # Step 2: cross family (CoulombInter, Ising) -- bond-diagonal,
    # q-independent, both mirrored diagonal target cells per representative.
    # AMENDED (2026-08-16, Task-6 granule adjudication): DIRECT placement
    # -- content lands at I=(m,a,b), the channel m carrying the declared
    # coefficient (delta_r[m] = +R_declared), matching the longitudinal
    # bare_bond_vertices Fock-diagonal rule exactly (#151-validated
    # precedent). The mirrored partner emits the same Re value at
    # (reverse[m], b, a). The PREVIOUS m<->reverse[m] swap (target at
    # reverse[m] with (a,b) unchanged, target at m with (b,a) swapped)
    # FAILED the multi-orbital off-site CoulombInter granule (L=3, norb=2,
    # a != b) at O(1) (0.12 vs tol 5.6e-4) -- invisible at a=b, which is
    # why Gate W0's norb=1 cross-family fixtures could not catch it; see
    # docs/superpowers/specs/2026-08-15-bond-transverse-design.md, "Cross
    # family", the 2026-08-16 amendment.
    for type_name, s_type in (("CoulombInter", -1.0), ("Ising", 1.0)):
        block_arr = coeffs.get(type_name)
        if block_arr is None:
            continue
        for (m, a, b) in reps:
            val = s_type * float(np.real(block_arr[m, a, b]))
            if val == 0.0:
                continue
            idx1 = m * nd + a * norb + b
            W[:, idx1, idx1] += val
            target2 = int(reverse[m])
            idx2 = target2 * nd + b * norb + a
            W[:, idx2, idx2] += val

    # Step 3: flip family (Exchange) -- the two AMENDED ordered records,
    # q-dependent ONLY through this local-channel phase.
    exch_arr = coeffs.get("Exchange")
    if exch_arr is not None:
        for (m, a, b) in reps:
            J = complex(exch_arr[m, a, b])
            if J == 0.0:
                continue
            R = delta_r[m].astype(float)
            phase = np.exp(1j * (q_flat @ R))
            idx_aa = a * norb + a
            idx_bb = b * norb + b
            W[:, idx_aa, idx_bb] += -1.0 * np.conj(J) * np.conj(phase)
            W[:, idx_bb, idx_aa] += -1.0 * J * phase

    return W


# =============================================================================
# #181 Tier 3 Phase A -- the bond-resolved LONGITUDINAL (spin/charge) vertex
# =============================================================================
#
# APPEND-ONLY REGION (continued). Spec:
# docs/superpowers/specs/2026-09-05-flex-offsite-bond-longitudinal-181-design.md,
# "The longitudinal bond vertices". The channel-0 block is the Tier-1
# locality-split pair-space matrix (hwave.solver.offsite
# .sc_matrices_from_split: on-site Kanamori slots + the off-site Hartree
# vertex V(q) in the density slots) placed VERBATIM; the off-site Fock
# (exchange) crossing that no q-only matrix can carry lives on the
# bond-DIAGONAL of the enlarged index, exactly where the Eliashberg module's
# bare_bond_vertices puts CoulombInter's (+Re v, -Re v) -- generalized here
# to every type through vertex_table's adjudicated ``cross`` row (Gate G1:
# equality with bare_bond_vertices for CoulombInter; Gate G2: exact
# diagonalization for Hund/Ising, tests/test_bond_longitudinal_ed.py).

def build_sc_bond_channel(topo, W0, channel, *, imag_tol=1e-12, types=None):
    """One channel (``"S"`` or ``"C"``) of the bond-resolved longitudinal
    vertex, shape ``(nvol, ND, ND)`` complex128, ``ND = B * norb**2``.

    ``W0`` (``(nvol, nd, nd)``, ``nd = norb**2``) is placed verbatim as
    the channel-0 block. For every type in ``types`` (default: every key
    of ``topo.coeffs``, in order) and every off-site reversal orbit
    ``{(m, a, b), (reverse[m], b, a)}`` (``m != 0``, one representative
    from :func:`iter_reversal_orbits`), the bond-diagonal entries at
    ``I = m*nd + a*norb + b`` and ``I' = reverse[m]*nd + b*norb + a``
    receive ``w_t * Re v_t[m, a, b]`` with ``w_t`` the type's ``cross``
    coefficient (S or C column of ``vertex_table.ADJUDICATED_SC``); both
    members of an orbit carry the same real part by Hermitian closure.
    A coefficient with ``|Im v| > imag_tol`` is refused (real
    coefficients only in this version; the imaginary direction of an
    off-site coefficient is not represented by this bond-diagonal rule).
    Every bond-diagonal entry is q-INDEPENDENT; the q-dependence of the
    vertex sits entirely in ``W0``'s density slots.

    FRAME. The result lives in the general path's pair frame: the one
    ``bubble.bond_bubble_static`` produces (its channel-0 block equals
    the solver's static ``chi0q`` slice to round-off -- measured
    difference 0.0 on the 2-orbital fixture, pinned at 1e-13), the one
    the Tier-1 ``W0`` is built in, and the one the Eliashberg module's
    ``bare_bond_vertices`` uses -- for CoulombInter the two vertex
    builders agree to round-off (Gate G1, pinned at 1e-13 in
    tests/test_bond_longitudinal_vertex.py). Structure: Hermitian at
    every q, and ``W(-q) == W(q)^T`` for real coefficients (NOT
    ``W(q)^dagger == W(-q)``: the density slot (aa,bb) carries
    ``V_ab(q) = V_ab(-q)^*``).
    """
    from hwave.solver.vertex_table import sc_coefficients
    if channel not in ("S", "C"):
        raise ValueError(
            "build_sc_bond_channel: channel must be 'S' or 'C', got {!r}"
            .format(channel))
    # re-validate (alias-safe): coeffs is a plain dict a caller may have
    # re-keyed since the topology was resolved
    topo = BondTopology(delta_r=topo.delta_r, reverse=topo.reverse,
                        coeffs=topo.coeffs)
    keys = list(topo.coeffs)
    if not keys:
        raise ValueError(
            "build_sc_bond_channel: the topology carries no coefficient "
            "arrays, so norb cannot be inferred")
    if types is None:
        types = keys
    else:
        if isinstance(types, str) or not isinstance(types, (list, tuple)):
            raise ValueError(
                "build_sc_bond_channel: types must be a list/tuple of type "
                "names, got {!r}".format(types))
        types = list(types)
    if (any(not isinstance(t, str) for t in types)
            or len(set(types)) != len(types)
            or any(t not in keys for t in types)):
        raise ValueError(
            "build_sc_bond_channel: types must be an ordered, duplicate-"
            "free subset of the topology's coefficient keys {}; got {!r}"
            .format(keys, types))
    try:
        imag_ok = np.isfinite(imag_tol) and imag_tol >= 0
    except TypeError:
        imag_ok = False
    if not imag_ok:
        raise ValueError(
            "build_sc_bond_channel: imag_tol must be a finite number >= 0, "
            "got {!r}".format(imag_tol))
    norb = int(topo.coeffs[keys[0]].shape[1])
    nd = norb * norb
    W0 = np.asarray(W0)
    if W0.dtype.kind not in ("f", "c"):
        raise ValueError(
            "build_sc_bond_channel: W0 must be a real or complex floating "
            "array (it is promoted to complex128); got dtype {}"
            .format(W0.dtype))
    if W0.ndim != 3 or W0.shape[1:] != (nd, nd):
        raise ValueError(
            "build_sc_bond_channel: W0 must have shape (nvol, {0}, {0}) "
            "(nd = norb**2 with norb={1} from the topology); got {2}"
            .format(nd, norb, W0.shape))
    if not np.all(np.isfinite(W0)):
        raise ValueError(
            "build_sc_bond_channel: W0 has non-finite entries")
    nvol = W0.shape[0]
    B = int(np.asarray(topo.delta_r).shape[0])
    ND = B * nd
    W = np.zeros((nvol, ND, ND), dtype=np.complex128)
    W[:, :nd, :nd] = W0
    col = 0 if channel == "S" else 1
    reverse = np.asarray(topo.reverse)
    delta_r = np.asarray(topo.delta_r)
    orbits = [orb for orb in iter_reversal_orbits(topo) if orb[0] != 0]
    for t in types:
        w_t = sc_coefficients(t, "cross")[col]
        arr = topo.coeffs[t]
        for (m, a, b) in orbits:
            v = complex(arr[m, a, b])
            if abs(v.imag) > imag_tol:
                raise ValueError(
                    "build_sc_bond_channel: the off-site {} coefficient at "
                    "channel {} (delta_r={}), orbitals ({}, {}) is complex "
                    "({}); only real off-site coefficients are supported "
                    "by the bond-resolved longitudinal channel in this "
                    "version".format(
                        t, m, tuple(int(x) for x in delta_r[m]), a, b, v))
            if w_t == 0.0 or v.real == 0.0:
                continue
            for I in (m * nd + a * norb + b,
                      int(reverse[m]) * nd + b * norb + a):
                W[:, I, I] += w_t * v.real
    return W


def W_sc_bond(topo, S0, C0, *, imag_tol=1e-12, types=None):
    """``(S_bond, C_bond)``: :func:`build_sc_bond_channel` for both
    channels (tests and Gate G1; production builds one channel at a time
    so that only one ``(nvol, ND, ND)`` vertex is alive)."""
    return (build_sc_bond_channel(topo, S0, "S", imag_tol=imag_tol, types=types),
            build_sc_bond_channel(topo, C0, "C", imag_tol=imag_tol, types=types))
