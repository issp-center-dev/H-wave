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

import numpy as np
from scipy.optimize import nnls
from scipy.sparse.linalg import LinearOperator

from . import backend as _bk
from . import matsubara as _ms

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
        ``+V_{l1 l2}(Delta r_m)`` (S) / ``-V_{l1 l2}(Delta r_m)`` (C)
        (spec S4.3). q-dependent (the Hartree lives in the ``m=0`` block).
    Vpp_s, Vpp_t : ndarray, shape (ND, ND), complex
        q-**independent** bare Cooper vertices. The local ``m=0`` block is the
        pair matrix ``L^{s/t}`` obtained from ``(S0 +- C0_local)`` under the
        ``(ac,db)->(ab,cd)`` pair<->ph crossing, with the **inter-site**
        (``R!=0``) Hartree excluded so ``V^pp`` carries only local
        interactions. The ``m!=0`` block is ``Q_eta (D_off + D_off^dag) Q_eta``
        (``= D_off (1 +- P)`` for real inversion-symmetric bonds), with
        ``D_off`` diagonal ``V_ab(Delta r_m)`` (``m!=0`` only) and
        ``P:(m,ab)->(m_bar,ba)`` the reversal+orbital-swap involution.
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
    # m!=0 bond-diagonal Fock: +V_ab(R_m) (S) / -V_ab(R_m) (C)
    for m in range(1, B):
        vb = v_bond[m]
        for l1 in range(norb):
            for l2 in range(norb):
                I = m * nd + l1 * norb + l2
                S_bond[:, :, :, I, I] = vb[l1, l2]
                C_bond[:, :, :, I, I] = -vb[l1, l2]

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
                hartree_rneq0 += v_bond[m][a, b]
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
        return
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
        return
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


def _bond_cond_scores(mat):
    """Per-point conditioning score for a batch of ``(ND, ND)`` solve
    matrices, ``mat`` shape ``(N, ND, ND)``.

    Mirrors, point for point, the two criteria :func:`_check_bond_conditioning`
    combines into its worst-point guard score (relative ``sigma_min/sigma_max``
    and absolute pole distance ``sigma_min/max(1, sigma_max)``, minimized) --
    see that function's docstring for the rationale. Unlike
    ``_check_bond_conditioning``, which reports only the single worst point and
    raises, this returns the score at every point so callers that need a full
    conditioning MAP (:func:`dress_bond_dynamic`) are not stuck redoing the SVD
    with a different formula.

    PINNING CONTRACT: this is a second, independently-maintained copy of
    ``_check_bond_conditioning``'s ratio/pole formula (unavoidable --
    ``_check_bond_conditioning`` has no per-point return path and is
    byte-invariant). ``tests/test_bond_channels_dynamic.py::
    test_bond_cond_scores_pins_check_bond_conditioning_decision_boundary``
    pins the two formulas together by asserting they agree on every
    raise/no-raise decision; keep that test green if this formula ever
    changes.
    """
    sv = np.linalg.svd(mat, compute_uv=False)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = sv[..., -1] / sv[..., 0]
    ratio = np.where(np.isfinite(ratio), ratio, 0.0)
    pole = sv[..., -1] / np.maximum(1.0, sv[..., 0])
    pole = np.where(np.isfinite(pole), pole, 0.0)
    return np.minimum(ratio, pole)


def dress_bond_dynamic(chi_bar_w, S_bond, C_bond, cond_tol=_BOND_COND_FLOOR):
    """Frequency-batched :func:`dress_bond`: dress the enlarged bond bubble at
    every bosonic Matsubara frequency in one batched solve.

    Identical algebra to :func:`dress_bond` -- ``chi_s = solve(I - chi_bar @
    S_bond, chi_bar)``, ``chi_c = solve(I + chi_bar @ C_bond, chi_bar)`` --
    applied independently at each ``(q, i nu)`` point: ``S_bond``/``C_bond``
    are q-dependent but frequency-independent (spec S4.3), so they are
    broadcast across the bosonic-frequency axis of ``chi_bar_w`` before the
    solve. Every ``(q, i nu)`` slice of the result is exactly what
    ``dress_bond(chi_bar_w[:, :, :, j], S_bond, C_bond)`` would produce for
    that ``j`` (see the slice-equality test) -- this function differs from a
    Python loop over ``dress_bond`` only in doing every ``j`` in one batched
    ``np.linalg.solve`` call.

    Parameters
    ----------
    chi_bar_w : ndarray, shape (Nx, Ny, Nz, nmat, ND, ND), complex
        The frequency-resolved enlarged bond bubble (:func:`bond_bubble_dynamic`).
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        The enlarged bare vertices (:func:`bare_bond_vertices`) -- q-dependent,
        frequency-independent.
    cond_tol : float or None, optional
        Same conditioning floor and same two criteria as :func:`dress_bond`
        (default :data:`_BOND_COND_FLOOR`), applied independently at every
        ``(q, i nu)`` point rather than only at ``q``: the RPA denominator can
        become ill-conditioned at a nonzero bosonic frequency even where its
        static (``i nu_0``) slice is fine. ``None`` disables the check.

    Returns
    -------
    chi_s_w, chi_c_w : ndarray, shape (Nx, Ny, Nz, nmat, ND, ND), complex
        The dressed spin (minus) and charge (plus) susceptibilities at the
        enlarged bond-major index, resolved over the bosonic-frequency axis.
    cond_min_spin, cond_min_charge : ndarray, shape (Nx, Ny, Nz, nmat), float64
        The per-``(q, i nu)`` conditioning score (the same
        ``min(sigma_min/sigma_max, sigma_min/max(1, sigma_max))`` combination
        :func:`_check_bond_conditioning` uses to decide pass/fail) of the spin
        and charge RPA denominators, respectively. Returned regardless of
        ``cond_tol`` (even ``None``) so callers can inspect the full map.

    Raises
    ------
    ValueError
        If either RPA denominator fails the ``cond_tol`` conditioning check at
        any ``(q, i nu)`` point.
    """
    chi_bar_w = np.asarray(chi_bar_w)
    if chi_bar_w.ndim != 6 or chi_bar_w.shape[4] != chi_bar_w.shape[5]:
        raise ValueError(
            "dress_bond_dynamic: chi_bar_w must have shape (Nx, Ny, Nz, nmat, "
            "ND, ND); got {}".format(chi_bar_w.shape)
        )
    Nx, Ny, Nz, nmat, ND, _ = chi_bar_w.shape

    S_bond = np.asarray(S_bond)
    C_bond = np.asarray(C_bond)
    q_shape = (Nx, Ny, Nz, ND, ND)
    if S_bond.shape != q_shape:
        raise ValueError(
            "dress_bond_dynamic: S_bond shape {} must be (Nx, Ny, Nz, ND, ND)"
            " = {} (matching chi_bar_w's q-grid and matrix size)".format(
                S_bond.shape, q_shape))
    if C_bond.shape != q_shape:
        raise ValueError(
            "dress_bond_dynamic: C_bond shape {} must be (Nx, Ny, Nz, ND, ND)"
            " = {} (matching chi_bar_w's q-grid and matrix size)".format(
                C_bond.shape, q_shape))

    # The q-dependent, frequency-independent bare vertices are broadcast
    # across the bosonic axis BY THE MATMUL, not materialised.
    #
    # MEMORY (review fix): the previous form built
    # ``np.broadcast_to(S_bond[..., None, :, :], chi_bar_w.shape)`` and then
    # ``.reshape(N, ND, ND)`` it. That reshape CANNOT be a view -- the
    # broadcast frequency axis has stride 0 -- so numpy silently copied, i.e.
    # two extra FULL ``(N_q * nmat, ND, ND)`` arrays (two frequency-resolved
    # units) that ``sc._bond_memory_estimate``'s dressing budget did not
    # count. numpy's matmul broadcasts the leading batch axes itself, so
    # keeping the 6-D form and reshaping only the RESULT (which is
    # C-contiguous, hence a free view) gives a BIT-IDENTICAL answer at two
    # units less peak. Verified array-wise (max |diff| = 0.0) against the
    # previous form; the identity is likewise materialised only as the
    # ``(ND, ND)`` matrix that broadcasts into the subtraction.
    N = Nx * Ny * Nz * nmat
    chi_flat = chi_bar_w.reshape(N, ND, ND)
    eye = np.eye(ND, dtype=complex)
    mat_s = (eye - chi_bar_w @ S_bond[:, :, :, None, :, :]).reshape(N, ND, ND)
    mat_c = (eye + chi_bar_w @ C_bond[:, :, :, None, :, :]).reshape(N, ND, ND)

    cond_min_spin = _bond_cond_scores(mat_s).reshape(Nx, Ny, Nz, nmat)
    cond_min_charge = _bond_cond_scores(mat_c).reshape(Nx, Ny, Nz, nmat)

    if cond_tol is not None:
        for name, cond_min in (("spin", cond_min_spin),
                                ("charge", cond_min_charge)):
            worst = float(cond_min.min())
            if worst <= cond_tol:
                iq = int(np.argmin(cond_min))
                qx, rem = divmod(iq, Ny * Nz * nmat)
                qy, rem = divmod(rem, Nz * nmat)
                qz, j = divmod(rem, nmat)
                raise ValueError(
                    "dress_bond_dynamic: the {} RPA denominator is singular "
                    "or nearly singular at q-point index ({}, {}, {}), "
                    "bosonic Matsubara index {}: conditioning score = "
                    "{:.3e} <= cond_tol = {:.3e}. The dynamic bond path has "
                    "entered the {} instability region, where the dressed "
                    "vertices are enormous and numerically meaningless. "
                    "Reduce the interaction strength, raise the temperature, "
                    "refine/reduce the q-grid, or -- if you deliberately "
                    "want to study the stiff regime -- lower "
                    "cond_tol.".format(
                        name, qx, qy, qz, j, worst, cond_tol, name))

    chi_s_flat = np.linalg.solve(mat_s, chi_flat)
    chi_c_flat = np.linalg.solve(mat_c, chi_flat)

    chi_s_w = chi_s_flat.reshape(Nx, Ny, Nz, nmat, ND, ND)
    chi_c_w = chi_c_flat.reshape(Nx, Ny, Nz, nmat, ND, ND)

    return chi_s_w, chi_c_w, cond_min_spin, cond_min_charge


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


def _g2_from_green(green, beta):
    """Static two-particle bubble ``G2[i,j,l,m](k) = (T) sum_n
    G_{ij}(k, iw_n) G_{lm}(-k, -iw_n)`` (``T = 1/beta``).

    The SINGLE shared implementation (roll+flip for ``G(-k, -w_n)``, per-site
    batched GEMM over the Matsubara axis) used by both this module
    (:func:`make_bond_kernel` / :func:`pair_weight`) and ``hwave.sc``, which
    imports it as ``bond_channels._g2_from_green`` and exposes it under the
    pre-existing name ``sc._calc_g2`` (review fix: this used to be a verbatim
    copy kept "to avoid a circular import" -- that justification did not
    hold, since ``bond_channels`` never imports ``sc``; only ``sc`` imports
    ``bond_channels``, so hoisting the implementation here and having ``sc``
    delegate to it introduces no cycle).
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
    return G2 / beta


def make_bond_kernel(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                     green, bond_set, pairing_type, beta, part="full"):
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
        pairing_type, beta)
    if part not in ops:
        raise ValueError(
            "make_bond_kernel: unknown part '{}'. Use one of {}.".format(
                part, sorted(ops)))
    return ops[part], vec_size


def make_bond_kernel_parts(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                           green, bond_set, pairing_type, beta):
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
        pairing_type, beta)
    return ops["full"], ops["fluctuation"], ops["instantaneous"], vec_size


def _bond_kernel_operators(chi_s, chi_c, S_bond, C_bond, Vpp_s, Vpp_t,
                           green, bond_set, pairing_type, beta):
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
    G2 = _g2_from_green(green, beta)                      # (i,j,l,m,Nx,Ny,Nz)
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
# make_bond_kernel_dynamic's memory contract  --  read together with
# sc._bond_memory_estimate (mirrors the BOND_BUBBLE_* contract below)
# ---------------------------------------------------------------------------
#
# The dynamic bond kernel is the memory-dominant object of the dynamic path:
# every array below carries the FULL fermionic/bosonic window ``nmat`` on top
# of the ``(N_q, ND, ND)`` footprint the static kernel already had. The
# resource preflight can only refuse an over-budget run if the counts here are
# not an UNDERCOUNT, so the builder is written to a fixed, documented working
# set (every temporary is released -- ``del`` -- as soon as it is consumed)
# and the counts are exported for ``sc._bond_memory_estimate`` to import.
#
# HOW THESE NUMBERS WERE OBTAINED. tracemalloc (which does trace numpy
# buffers -- verified against a bare ``np.empty``), started BEFORE any
# allocation so the live set is a true baseline, then ``reset_peak()`` and one
# build / one matvec. The reported figure is the CUMULATIVE peak over the whole
# phase, not a per-stage transient: a per-stage measurement resets between
# stages and therefore misses arrays that stay live ACROSS a stage boundary,
# which is exactly what dominates here (``Z`` is still alive while
# ``fermion_to_tau`` allocates its own three). Per-stage measurement reports
# ~3.1 units for the matvec; the true cumulative peak is ~4.3. Use the
# cumulative method if these are ever re-derived.
#
# Measured (numpy 1.26), in units of the phase's own buffer size:
#
#   configuration                    build   matvec
#   norb=1 B=5 ND=5   8x8x1 nmat=16   3.58     4.24
#   norb=1 B=5 ND=5   6x6x1 nmat=32   3.29     4.24
#   norb=2 B=3 ND=12  4x4x1 nmat=16   3.23     4.40
#
# The declared constants are the CEIL of the measured maxima, so the preflight
# rounds against itself and can never undercount.
#
# --- build phase, arrays of size (N_q * nmat * ND**2) complex ---------------
# BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS -- high-water while the vertex is hoisted.
# Contributions: the ``S chi S`` / ``C chi C`` sandwich holds {spin term,
# charge term, their difference}, and ``_ms.boson_to_tau(F_q_w, axis=3)`` holds
# {``F_q_w``, its fft output, the phase-multiplied result}. ``spatial_ifftn``
# afterwards holds only input + one output. Measured max 3.58 -> declared 4.
# The caller-owned ``chi_s_w`` / ``chi_c_w`` are the SAME size and are NOT
# counted here; the preflight budgets those separately, exactly as ``green_kw``
# is excluded from BOND_BUBBLE_N2_BUFFERS.
#
# BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT = arrays of that size that stay alive
# for the whole lifetime of the returned LinearOperator:
#   1. ``F_rt``  -- the (r, tau) fluctuation vertex, hoisted out of the
#                   per-matvec closure (it does not depend on phi)
#
# --- per-matvec, arrays of size (N_q * nmat * ND) complex -------------------
# BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS -- high-water during one ``matvec``,
# reached inside ``_ms.fermion_to_tau(Z, axis=3)`` while ``Z`` is still live:
# {``Z``, the fft's internal axis-reordering copy, the fft output, the
# phase-multiplied result}. Measured max 4.40 -> declared 5.
#
# Two structural choices in ``_fluctuation`` keep it at that level rather than
# higher (both measured, see its docstring):
#   * the transform calls are SPLIT, not nested -- nesting keeps the caller's
#     argument alive across the inner call's temporaries and adds a unit;
#   * the real-space multiply-accumulate is a batched MATMUL in the (ND, ND)
#     matrix form, NOT an ``np.einsum`` over the split [m, ab, m', cd] index.
#     einsum buffers the strided operands (4.24 units for that stage alone
#     against the matmul's 1.01) and pushed the whole matvec to 5.45 units even
#     with every transform already split. The MAC, not the transforms, is the
#     stage that decides this constant.
# Reverting either choice silently invalidates the constant and the preflight
# built on it.
#
# NOT counted here, and NOT safely foldable into the count above:
#   * ``phi`` / ``X_w`` (size ``N_q * nmat * nd``) -- smaller by a factor ``B``.
#   * the persistent ``G2p_w`` (size ``N_q * nmat * nd**2``). Its ratio to one
#     matvec unit is ``nd / B``, NOT ``1 / B**2``: it is sub-dominant only when
#     ``B > nd``, and at ``B = 1`` (the single on-site channel) it is ``nd``
#     TIMES a matvec unit. Task 8's preflight must budget ``G2p_w`` as its own
#     explicit term rather than absorbing it into
#     BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS.
BOND_KERNEL_DYNAMIC_VERTEX_BUFFERS = 4
BOND_KERNEL_DYNAMIC_VERTEX_PERSISTENT = 1
BOND_KERNEL_DYNAMIC_MATVEC_BUFFERS = 5


def make_bond_kernel_dynamic(chi_s_w, chi_c_w, S_bond, C_bond, Vpp_s, Vpp_t,
                             green, bond_set, pairing_type, beta):
    r"""Frequency-resolved (dynamic) bond Eliashberg kernel on the UNIFORM
    Matsubara axis (spec 2026-07-27-dynamic-bond-channels-design.md S3.3).

    The frequency extension of :func:`make_bond_kernel`: the gap keeps its
    fermionic Matsubara axis, ``phi_{l1l2}(k, i w_n)``, and the enlarged bond
    vertex acquires a bosonic transfer axis, so the fluctuation part becomes a
    convolution in ``(q, i nu) = (k - k', i w_n - i w_n')`` instead of a purely
    spatial one::

        Y^fl_{ab}(k, n) = -(T/N) sum_{k',n'} sum_{m,m'}
              e^{i k.dr_m} e^{i k'.dr_m'}
              F_eta[(m,ab),(m',cd)](k - k', n - n') X_{cd}(k', n')

        X_{cd}(k', n') = sum_{ef} G2_{c,e,d,f}(k', n') phi_{ef}(k', n')

    with the frequency-resolved pair bubble ``G2(k, n) = (T) G(k, i w_n)
    G(-k, -i w_n)`` (:func:`eliashberg_dynamic.calc_g2_dynamic`, which already
    carries the single ``T = 1/beta``) and the dressed fluctuation vertex,
    now evaluated at every bosonic transfer::

        F_s(q, i nu) = +1.5 S chi_s(q, i nu) S - 0.5 C chi_c(q, i nu) C
        F_t(q, i nu) = -0.5 S chi_s(q, i nu) S - 0.5 C chi_c(q, i nu) C

    The bosonic zero transfer sits at array index ``nmat // 2`` (the
    :func:`bond_bubble_dynamic` map), and the fermionic difference ``n - n'``
    is taken on the stored-window torus -- both are realized implicitly by the
    imaginary-time product below, exactly as in the scalar dynamic kernel
    ``eliashberg_dynamic.eliashberg_kernel_dynamic``.

    The instantaneous (bare Cooper) part is FREQUENCY-INDEPENDENT: the bare
    ``V^pp`` has no bosonic transfer, so the internal ``k', n'`` reduction
    collapses to the window Matsubara sum of ``X`` and the result is flat in
    the external frequency (spec S3.3, uniform normative definition)::

        A_{m; cd}      = (T/N) sum_{k', n'} e^{-i k'.dr_m} X_{cd}(k', n')
        Y^pp_{ab}(k)   = -1/2 sum_{m,m'} [V^pp_eta]_{(m,ab),(m',cd)}
                              e^{i k.dr_m'} A_{m; cd}

    Because ``sum_n calc_g2_dynamic(...) == _g2_from_green(...)`` exactly
    (pinned by ``tests/test_eliashberg_dynamic.py::
    test_g2_dynamic_sums_to_static``), this is literally the static
    :func:`_bond_kernel_operators` instantaneous body applied to the
    frequency-summed ``X`` and broadcast flat over the external frequency
    axis -- no second definition of the Cooper term exists.

    Normalization bookkeeping (same as the scalar dynamic kernel): the single
    ``T = 1/beta`` lives inside ``G2``; the single spatial ``1/N`` comes from
    the ``ifftn`` half of the FFT convolution (numpy's inverse transform
    divides by ``N``, the forward one does not re-multiply); the frequency
    transforms ``boson_to_tau`` / ``fermion_to_tau`` / ``tau_to_fermion`` are
    normalization-free in the same pairing. No extra factors are applied.

    Parameters
    ----------
    chi_s_w, chi_c_w : ndarray, shape (Nx, Ny, Nz, nmat, ND, ND), complex
        Dressed spin/charge susceptibilities on the bosonic Matsubara axis
        (:func:`dress_bond_dynamic`). ``ND = norb**2 * B``.
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND), complex
        Bare enlarged spin/charge vertices (:func:`bare_bond_vertices`); they
        carry no frequency dependence and are broadcast over ``nmat``.
    Vpp_s, Vpp_t : ndarray, shape (ND, ND), complex
        q- and frequency-independent bare Cooper vertices.
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat), complex
        k-space Matsubara Green function in the ``sc.py`` layout.
    bond_set : ResolvedInteractionSet
        Bond-channel topology (``delta_r``, ``n_channels``).
    pairing_type : {"singlet", "triplet"}
    beta : float
        Inverse temperature.

    Returns
    -------
    A : scipy.sparse.linalg.LinearOperator
        ``A.matvec(phi.ravel())`` returns ``(K phi).ravel()`` with ``phi``
        shaped ``(norb, norb, Nx, Ny, Nz, nmat)``.
    vec_size : int
        ``norb**2 * Nx * Ny * Nz * nmat``.
    """
    from . import eliashberg_dynamic as _ed

    green = _validate_green_beta(green, beta, "make_bond_kernel_dynamic")
    norb = green.shape[0]
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nx, Ny, Nz, nmat = (green.shape[2], green.shape[3], green.shape[4],
                        green.shape[5])
    N = Nx * Ny * Nz
    delta_r = bond_set.delta_r

    chi_s_w = np.asarray(chi_s_w)
    chi_c_w = np.asarray(chi_c_w)
    expected = (Nx, Ny, Nz, nmat, ND, ND)
    for name, arr in (("chi_s_w", chi_s_w), ("chi_c_w", chi_c_w)):
        if arr.shape != expected:
            raise ValueError(
                "make_bond_kernel_dynamic: {} has shape {} but expected "
                "(Nx, Ny, Nz, nmat, ND, ND) = {} (ND = norb**2 * "
                "n_channels). The dynamic path needs the frequency-resolved "
                "susceptibility from dress_bond_dynamic, not the static "
                "dress_bond slice.".format(name, arr.shape, expected))
    S_bond = np.asarray(S_bond)
    C_bond = np.asarray(C_bond)
    for name, arr in (("S_bond", S_bond), ("C_bond", C_bond)):
        if arr.shape != (Nx, Ny, Nz, ND, ND):
            raise ValueError(
                "make_bond_kernel_dynamic: {} has shape {} but expected "
                "(Nx, Ny, Nz, ND, ND) = {}.".format(
                    name, arr.shape, (Nx, Ny, Nz, ND, ND)))

    if pairing_type == "singlet":
        Vpp = np.asarray(Vpp_s)
    elif pairing_type == "triplet":
        Vpp = np.asarray(Vpp_t)
    else:
        raise ValueError(
            "make_bond_kernel_dynamic: unknown pairing_type '{}'. Use "
            "'singlet' or 'triplet'.".format(pairing_type))
    if Vpp.shape != (ND, ND):
        raise ValueError(
            "make_bond_kernel_dynamic: Vpp has shape {} but expected "
            "(ND, ND) = {}.".format(Vpp.shape, (ND, ND)))

    vec_size = nd * N * nmat

    # --- hoisted invariants (see the memory contract above) ---------------
    # F_eta(q, i nu): the SAME algebra as the static builder, broadcast over
    # the bosonic axis (S/C are frequency-independent).
    S_b = S_bond[:, :, :, None, :, :]
    C_b = C_bond[:, :, :, None, :, :]
    if pairing_type == "singlet":
        F_q_w = 1.5 * (S_b @ chi_s_w @ S_b) - 0.5 * (C_b @ chi_c_w @ C_b)
    else:
        F_q_w = -0.5 * (S_b @ chi_s_w @ S_b) - 0.5 * (C_b @ chi_c_w @ C_b)

    # (q, i nu) -> (r, tau): the vertex leg of the convolution. boson_to_tau
    # puts the bosonic array on the tau nodes, where the product with the
    # fermionic-tau gap leg IS the circular frequency convolution with the
    # spec's transfer map -- the same pattern as
    # eliashberg_dynamic.vertex_qw_to_rt. spatial_ifftn carries the 1/N.
    F_tau = _ms.boson_to_tau(F_q_w, axis=3)
    del F_q_w
    # Kept in the (ND, ND) matrix form, NOT split into [m, ab, m', cd]: the
    # real-space multiply-accumulate below is exactly the ND-space mat-vec
    # W_I(r, tau) = sum_J F_IJ(r, tau) Z_J(r, tau), which numpy can do as one
    # batched matmul (see _fluctuation for why that matters for memory).
    F_rt = _bk.spatial_ifftn(F_tau, axes=(0, 1, 2))   # (Nx,Ny,Nz,nmat,ND,ND)
    del F_tau

    # Frequency-resolved pair bubble; carries the single T = 1/beta. REUSED
    # from the scalar dynamic path -- no fourth roll+flip copy is introduced.
    G2_w = _ed.calc_g2_dynamic(green, beta)      # (i,j,l,m,Nx,Ny,Nz,nmat)
    G2p_w = G2_w.transpose(0, 2, 1, 3, 4, 5, 6, 7).reshape(nd, nd, N, nmat)
    del G2_w

    # Bond phases PH[m](k) = e^{i k.dr_m} -- purely spatial, broadcast over
    # the frequency axis. Copied verbatim from _bond_kernel_operators.
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

    def _fluctuation(X_grid_w):
        """Y^fl(k, n): spatial FFT convolution x imaginary-time product.

        Two structural choices here are load-bearing for the memory contract
        declared above; both were measured, not guessed (numpy 1.26, in units
        of one ``N_q * nmat * ND`` buffer):

        * The transforms are NOT nested. Each stage is bound to a name and the
          previous one ``del``-ed before the next allocates; nesting them
          (``spatial_ifftn(fermion_to_tau(Z))``) keeps the caller's argument
          alive across the inner call's own two temporaries and adds a unit.
        * The real-space multiply-accumulate is a batched MATMUL in the
          (ND, ND) matrix form, not an ``einsum`` over the split
          ``[m, ab, m', cd]`` index. The two are mathematically identical
          (agreement 3.6e-15) but ``np.einsum`` buffers the strided operands
          and peaks at 4.24 units for that stage alone, where the matmul peaks
          at 1.01 -- this single stage, not the transforms, is what would
          otherwise dominate the whole matvec (measured 5.45 units with the
          einsum, 4.24 with the matmul; every transform split in both).

        Everything is carried in the ``(Nx, Ny, Nz, nmat, B, nd)`` layout so
        the reshape to the ND matrix index is free (no ``moveaxis`` copy) and
        the FFT axes stay leading.
        """
        # X_grid_w is (nd, Nx, Ny, Nz, nmat); the kernel wants (..., nmat, nd).
        X_t = np.moveaxis(X_grid_w, 0, -1)
        # Z_{m,ab}(k, n) = e^{i k.dr_m} X_{ab}(k, n)
        Z = PH_k[:, :, :, None, :, None] * X_t[:, :, :, :, None, :]
        Z_t = _ms.fermion_to_tau(Z, axis=3)
        del Z
        Z_rt = _bk.spatial_ifftn(Z_t, axes=(0, 1, 2))
        del Z_t
        # W_I(r, tau) = sum_J F_IJ(r, tau) Z_J(r, tau), I = (m, ab)
        W_rt = (F_rt @ Z_rt.reshape(Nx, Ny, Nz, nmat, ND, 1))[..., 0]
        del Z_rt
        W_kr = _bk.spatial_fftn(W_rt, axes=(0, 1, 2))
        del W_rt
        W_kw = _ms.tau_to_fermion(W_kr, axis=3)
        del W_kr
        # (Nx,Ny,Nz,nmat,nd)
        return -np.einsum('xyzm,xyztma->xyzta', PH_k,
                          W_kw.reshape(Nx, Ny, Nz, nmat, B, nd))

    def _instantaneous(X_sum):
        """Y^pp(k): VERBATIM the static _bond_kernel_operators body, applied
        to the frequency-summed X (spec S3.3). Returns (Nx,Ny,Nz,nd)."""
        A_coef = (1.0 / N) * np.einsum('mxyz,cxyz->mc', PHc, X_sum)  # (B,nd)
        Bcoef = np.einsum('aAbB,aB->bA', Vpp6, A_coef)               # (B,nd)
        return -0.5 * np.einsum('mxyz,ma->xyza', PH, Bcoef)

    def matvec(v):
        phi = np.asarray(v).reshape(norb, norb, Nx, Ny, Nz, nmat)
        phi_flat = phi.reshape(nd, N, nmat)              # [(e,f), k, n]
        X = np.einsum('adkn,dkn->akn', G2p_w, phi_flat)  # (nd, N, nmat)
        X_grid_w = X.reshape(nd, Nx, Ny, Nz, nmat)
        del X

        Y = _fluctuation(X_grid_w)                       # (Nx,Ny,Nz,nmat,nd)
        # The bare Cooper term has no bosonic transfer: it sees only the
        # window sum of X and is flat in the external frequency.
        Y += _instantaneous(X_grid_w.sum(axis=-1))[:, :, :, None, :]

        Y = Y.reshape(Nx, Ny, Nz, nmat, norb, norb)
        Y = np.moveaxis(Y, (4, 5), (0, 1))   # (norb,norb,Nx,Ny,Nz,nmat)
        return Y.ravel()

    A = LinearOperator((vec_size, vec_size), matvec=matvec, dtype=complex)
    return A, vec_size


def make_bond_kernel_dynamic_ir(chi_s_w, chi_c_w, S_bond, C_bond,
                                Vpp_s, Vpp_t, green_ir, bond_set,
                                pairing_type, beta, axF, axB):
    r"""The IR-axis twin of :func:`make_bond_kernel_dynamic` (spec
    ``2026-07-27-dynamic-bond-channels-design.md`` S3.3, "IR realization").

    Same operator, same normative equations; only the frequency axis changes.
    The gap lives on the fermionic IR SAMPLING frequencies
    (``axF.freq_n``, ``axF.n_freq`` of them) instead of the centered uniform
    window, and every Matsubara transform is an IRAxis matmul instead of a
    phase-twisted FFT. The structure follows
    ``eliashberg_dynamic.eliashberg_kernel_ir`` line by line -- it is the
    scalar (B = 1) special case of this function.

    INPUT CONTRACT (explicit, because both conventions are defensible and the
    wiring must not have to guess): this function takes arrays ALREADY ON THE
    IR SAMPLING NODES and never fits uniform data itself.

    * ``chi_s_w``, ``chi_c_w``: ``(Nx, Ny, Nz, axB.n_freq, ND, ND)`` --
      dressed susceptibilities on the BOSONIC IR sampling frequencies.
    * ``green_ir``: ``(norb, norb, Nx, Ny, Nz, axF.n_freq)`` -- the Green
      function on the FERMIONIC IR sampling frequencies.

    Both are what ``eliashberg_dynamic._ir_compress`` produces from the
    uniform-grid FLEX output, which is exactly how ``solve_dynamic`` already
    prepares ``chis_w``/``green_w`` on the scalar IR path; the caller keeps
    ownership of the fit (including the ``drop_constant`` / ``ir_keep_static_
    chi`` policy for the susceptibilities, which is a data-quality decision,
    not a kernel decision). The static vertices ``S_bond``/``C_bond``/``Vpp``
    are frequency-independent and identical on both axes.

    Normalization. The IR transforms are PHYSICAL (``freq -> tau`` carries the
    1/beta Matsubara sum, ``tau -> freq`` is the integral over tau), while the
    uniform-FFT chain of :func:`make_bond_kernel_dynamic` carries one net
    factor beta; the explicit ``beta`` factor on the return restores the
    identical operator normalization. It multiplies BOTH parts, exactly as in
    ``eliashberg_kernel_ir``: the instantaneous term's equal-time value
    ``X(tau=0+) = (1/beta) sum_n X(i w_n)`` becomes the uniform path's window
    sum ``sum_n X`` after that factor.

    Instantaneous term (spec S3.3, "IR (normative for the IR path)"): the
    infinite-cutoff equal-time value ``sum_l X_l u_l(0+)`` via
    ``axF.u_zero_plus``, NOT a truncated window sum. This is the existing
    scalar IR convention (``eliashberg_kernel_ir`` :1067-1077); X is an
    anomalous (pair) amplitude, hence continuous at tau = 0, so the equal-time
    value needs no 0^+ regularization -- but it does differ from the uniform
    path's window sum by the ordinary ~1/w^2 tail, whose size is what
    :func:`tail_estimate` reports and what spec S3.5's inequalities budget.

    Parameters
    ----------
    chi_s_w, chi_c_w : ndarray, shape (Nx, Ny, Nz, axB.n_freq, ND, ND)
    S_bond, C_bond : ndarray, shape (Nx, Ny, Nz, ND, ND)
    Vpp_s, Vpp_t : ndarray, shape (ND, ND)
    green_ir : ndarray, shape (norb, norb, Nx, Ny, Nz, axF.n_freq)
    bond_set : ResolvedInteractionSet
    pairing_type : {"singlet", "triplet"}
    beta : float
    axF, axB : hwave.solver.ir_axis.IRAxis
        Fermionic / bosonic axes of the run (same beta, wmax, tol); build them
        with ``eliashberg_dynamic._ir_axes_for_run``.

    Returns
    -------
    A : scipy.sparse.linalg.LinearOperator
        ``A.matvec(phi.ravel())`` with ``phi`` shaped
        ``(norb, norb, Nx, Ny, Nz, axF.n_freq)``.
    vec_size : int
        ``norb**2 * Nx * Ny * Nz * axF.n_freq``.
    """
    from . import eliashberg_dynamic as _ed

    green_ir = _validate_green_beta(green_ir, beta,
                                    "make_bond_kernel_dynamic_ir")
    norb = green_ir.shape[0]
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nx, Ny, Nz, nfreq = (green_ir.shape[2], green_ir.shape[3],
                         green_ir.shape[4], green_ir.shape[5])
    N = Nx * Ny * Nz
    delta_r = bond_set.delta_r

    if float(axF.beta) != float(beta) or float(axB.beta) != float(beta):
        raise ValueError(
            "make_bond_kernel_dynamic_ir: the IR axes carry beta = ({}, {}) "
            "but the kernel was given beta = {}. Build both axes for this "
            "run's temperature (eliashberg_dynamic._ir_axes_for_run)."
            .format(axF.beta, axB.beta, beta))
    if nfreq != axF.n_freq:
        raise ValueError(
            "make_bond_kernel_dynamic_ir: green_ir's frequency axis has "
            "length {} but the fermionic IR axis has {} sampling "
            "frequencies. This function takes data ALREADY on the IR nodes "
            "(eliashberg_dynamic._ir_compress), not a uniform grid."
            .format(nfreq, axF.n_freq))

    chi_s_w = np.asarray(chi_s_w)
    chi_c_w = np.asarray(chi_c_w)
    expected = (Nx, Ny, Nz, axB.n_freq, ND, ND)
    for name, arr in (("chi_s_w", chi_s_w), ("chi_c_w", chi_c_w)):
        if arr.shape != expected:
            raise ValueError(
                "make_bond_kernel_dynamic_ir: {} has shape {} but expected "
                "(Nx, Ny, Nz, axB.n_freq, ND, ND) = {} (ND = norb**2 * "
                "n_channels). The IR path needs the susceptibility on the "
                "BOSONIC IR sampling frequencies, not on the uniform grid "
                "and not on the fermionic axis.".format(
                    name, arr.shape, expected))
    S_bond = np.asarray(S_bond)
    C_bond = np.asarray(C_bond)
    for name, arr in (("S_bond", S_bond), ("C_bond", C_bond)):
        if arr.shape != (Nx, Ny, Nz, ND, ND):
            raise ValueError(
                "make_bond_kernel_dynamic_ir: {} has shape {} but expected "
                "(Nx, Ny, Nz, ND, ND) = {}.".format(
                    name, arr.shape, (Nx, Ny, Nz, ND, ND)))

    if pairing_type == "singlet":
        Vpp = np.asarray(Vpp_s)
    elif pairing_type == "triplet":
        Vpp = np.asarray(Vpp_t)
    else:
        raise ValueError(
            "make_bond_kernel_dynamic_ir: unknown pairing_type '{}'. Use "
            "'singlet' or 'triplet'.".format(pairing_type))
    if Vpp.shape != (ND, ND):
        raise ValueError(
            "make_bond_kernel_dynamic_ir: Vpp has shape {} but expected "
            "(ND, ND) = {}.".format(Vpp.shape, (ND, ND)))

    vec_size = nd * N * nfreq

    # --- hoisted invariants -------------------------------------------------
    # F_eta(q, i nu) on the bosonic nodes: the SAME algebra as the uniform
    # builder (S/C are frequency-independent and broadcast over the axis).
    S_b = S_bond[:, :, :, None, :, :]
    C_b = C_bond[:, :, :, None, :, :]
    if pairing_type == "singlet":
        F_q_w = 1.5 * (S_b @ chi_s_w @ S_b) - 0.5 * (C_b @ chi_c_w @ C_b)
    else:
        F_q_w = -0.5 * (S_b @ chi_s_w @ S_b) - 0.5 * (C_b @ chi_c_w @ C_b)

    # (q, bosonic nodes) -> (r, FERMIONIC tau nodes): the vertex leg. The
    # product F*X is anti-periodic, so the tau product lives on the fermionic
    # nodes and the bosonic coefficients are evaluated there exactly -- the
    # _ir_vertex_to_rtau domain rule, which IS the normative treatment of the
    # transfer frequency on the IR path (spec S3.1, "Frequency grids (IR)").
    # IRAxis transforms contract the LAST axis, so the frequency axis is moved
    # out of position 3 and back; the transforms are split (bound, then the
    # input released) for the same memory reason as the uniform kernel.
    F_l = axB.fit_from_freq(np.moveaxis(F_q_w, 3, -1))   # (...,ND,ND,L_B)
    del F_q_w
    F_tau = axB.eval_to_tau_points(F_l, axF.tau)         # (...,ND,ND,n_tau)
    del F_l
    F_tau = np.moveaxis(F_tau, -1, 3)                    # (...,n_tau,ND,ND)
    F_rt = _bk.spatial_ifftn(F_tau, axes=(0, 1, 2))
    del F_tau

    # Frequency-resolved pair bubble ON THE IR NODES; carries the single
    # T = 1/beta. calc_g2_dynamic's G(-k,-i w) is a reverse+roll in k and a
    # plain reversal on the frequency axis -- valid here because the IR node
    # set is EXACTLY symmetric under integer negation (ir_axis.IRAxis
    # construction), which is the same property solve_dynamic's IR path
    # relies on.
    G2_w = _ed.calc_g2_dynamic(green_ir, beta)
    G2p_w = G2_w.transpose(0, 2, 1, 3, 4, 5, 6, 7).reshape(nd, nd, N, nfreq)
    del G2_w

    # Bond phases PH[m](k) = e^{i k.dr_m}: purely spatial, identical to the
    # uniform kernel (copied verbatim from _bond_kernel_operators).
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
    u0 = np.asarray(axF.u_zero_plus)                      # (L_F,)

    def _fluctuation_from_coeffs(X_l):
        """Y^fl on the fermionic IR nodes, from X's IR COEFFICIENTS.

        The bond phase e^{i k.dr_m'} is frequency-independent, so the
        coefficients of Z_{m'} = e^{i k.dr_m'} X are just the phase times X's
        coefficients -- one fit serves both this leg and the instantaneous
        term below (the uniform kernel needs no such sharing: there the
        window sum is a plain reduction of X itself).
        """
        # X_l: (Nx,Ny,Nz,nd,L) -> Z_l: (Nx,Ny,Nz,B,nd,L)
        Z_l = PH_k[:, :, :, :, None, None] * X_l[:, :, :, None, :, :]
        Z_t = np.moveaxis(axF.eval_to_tau(Z_l), -1, 3)    # (...,n_tau,B,nd)
        del Z_l
        Z_rt = _bk.spatial_ifftn(Z_t, axes=(0, 1, 2))
        del Z_t
        n_tau = Z_rt.shape[3]
        # W_I(r, tau) = sum_J F_IJ(r, tau) Z_J(r, tau), I = (m, ab): the same
        # batched ND-space matmul as the uniform kernel (an einsum over the
        # split [m,ab,m',cd] index costs ~4x the peak memory there).
        W_rt = (F_rt @ Z_rt.reshape(Nx, Ny, Nz, n_tau, ND, 1))[..., 0]
        del Z_rt
        W_kr = _bk.spatial_fftn(W_rt, axes=(0, 1, 2))
        del W_rt
        W_kw = np.moveaxis(axF.tau_to_freq(np.moveaxis(W_kr, 3, -1)), -1, 3)
        del W_kr
        return -np.einsum('xyzm,xyztma->xyzta', PH_k,
                          W_kw.reshape(Nx, Ny, Nz, nfreq, B, nd))

    def _instantaneous(X0):
        """Y^pp(k) from the equal-time X(tau=0+) (spec S3.3, IR definition).

        Structurally the uniform kernel's ``_instantaneous`` body with the
        window sum replaced by the infinite-cutoff equal-time value; the
        overall ``beta`` factor on the return turns X(tau=0+) into the
        uniform path's ``sum_n X``.
        """
        A_coef = (1.0 / N) * np.einsum('mxyz,xyzc->mc', PHc, X0)   # (B,nd)
        Bcoef = np.einsum('aAbB,aB->bA', Vpp6, A_coef)             # (B,nd)
        return -0.5 * np.einsum('mxyz,ma->xyza', PH, Bcoef)

    def matvec(v):
        phi = np.asarray(v).reshape(norb, norb, Nx, Ny, Nz, nfreq)
        phi_flat = phi.reshape(nd, N, nfreq)             # [(e,f), k, n]
        X = np.einsum('adkn,dkn->akn', G2p_w, phi_flat)  # (nd, N, n_freq)
        # (Nx,Ny,Nz,nd,n_freq) -> IR coefficients on the last axis
        X_l = axF.fit_from_freq(
            np.moveaxis(X.reshape(nd, Nx, Ny, Nz, nfreq), 0, 3))
        del X

        Y = _fluctuation_from_coeffs(X_l)                # (Nx,Ny,Nz,nfreq,nd)
        # The bare Cooper term has no bosonic transfer: it sees only the
        # equal-time X and is flat in the external frequency.
        Y += _instantaneous(X_l @ u0)[:, :, :, None, :]
        del X_l

        Y = Y.reshape(Nx, Ny, Nz, nfreq, norb, norb)
        Y = np.moveaxis(Y, (4, 5), (0, 1))   # (norb,norb,Nx,Ny,Nz,n_freq)
        return (beta * Y).ravel()

    A = LinearOperator((vec_size, vec_size), matvec=matvec, dtype=complex)
    return A, vec_size


def tail_estimate(X_w, beta, nmat):
    r"""Estimate the Matsubara-window TRUNCATION of the instantaneous term.

    Spec ``2026-07-27-dynamic-bond-channels-design.md`` S3.3, "Cross-basis
    validation", implemented verbatim. The uniform and IR paths use different
    (each self-consistent) normative definitions of the instantaneous
    (delta(tau)) Cooper term -- the uniform one sums the STORED window,
    ``T sum_{n in S} X``, while the IR one evaluates the infinite-cutoff
    equal-time value ``sum_l X_l u_l(0+)``. The difference is the ordinary
    ``~1/w^2`` tail of the pair amplitude X (which is continuous at tau = 0,
    so there is no delta-function/jump ambiguity, only this tail). This
    function estimates its size so the two paths can be compared against a
    quantitative budget instead of a hand-picked tolerance.

    It is an ESTIMATOR, not a bound (spec's words): the returned numbers are
    only as good as the ``|a|/w^2 + |b|/|w|^3`` asymptotic MAGNITUDE model on
    the outer shell, which is exactly what ``unreliable`` reports. The
    second term uses ``|w|``, not ``w``: ``x`` is a nonnegative magnitude
    profile (step 1), so its model must be even in w by construction --
    fitting a SIGNED ``1/w^3`` column against an even target is inconsistent
    on the negative-frequency half of the outer shell (review fix). With
    both columns even (and both strictly positive, so the design is far more
    collinear than the old odd/even pair), the coefficients are recovered
    with NONNEGATIVE least squares (``scipy.optimize.nnls``), not
    unconstrained ``lstsq`` followed by ``abs()``: on real physical data the
    unconstrained fit can return a genuinely negative coefficient (the two
    columns nearly trade off against each other), and forcing ``abs()``
    afterwards assembles a model that was never actually fit -- NNLS solves
    the ``a, b >= 0``-constrained problem directly, so it is never worse
    than any nonnegative guess.

    Pipeline (each step is normative):

    1. *Reduction first* -- one scalar profile ``x(n~) = max_{k, orbital
       pairs} |X(k, i w_n~)|`` over every non-frequency axis. A single fit is
       done on this profile; there are no per-k fits anywhere.
    2. *Outer shell* -- the stored frequencies with ``|n~| > 3*nmat/8`` (the
       outermost quarter of the window, both signs).
    3. *Asymptotic fit* -- nonnegative least squares of ``x(n~) ~ a/w^2 +
       b/|w|^3`` (``a, b >= 0``) over the outer shell.
    4. *Window sets* -- stored ``S = {n~: -nmat/2 <= n~ <= nmat/2 - 1}`` (the
       S3.1 integer map) and tail ``T_bar = {n~: n~ not in S, |n~| <=
       64*nmat}``; note the unstored ``n~ = +nmat/2`` belongs to the tail set.
    5. *Estimator* -- ``tail_est = T * sum_{n~ in T_bar} (|a|/w^2 + |b|/|w|^3)``
       by explicit summation of the fitted model over that closed window.
    6. *Relative form* -- ``tail_est_rel = tail_est / (T * sum_{n~ in S}
       x(n~))``: the denominator is the same-reduction sum over the stored
       set, i.e. the very quantity the tail extends.
    7. *Reliability* -- relative (Euclidean) fit residual over the outer shell
       ``> 0.2`` sets ``unreliable``. The flag is evaluated for the SAME
       ``|a|, |b|`` model that produced ``tail_est``, so it cannot report a
       residual for a model the estimator does not use. In addition, an outer
       shell with fewer than 4 points always sets ``unreliable`` -- with only
       2 or 3 points the 2-parameter fit is under- or exactly-determined
       (residual collapses to ~0 regardless of how well the model actually
       describes the tail), so the residual criterion alone cannot be
       trusted there.

    Parameters
    ----------
    X_w : ndarray
        The pair amplitude ``X_{l3l4}(k', i w_n')`` (or any array whose LAST
        axis is the fermionic Matsubara axis of the centered uniform grid --
        every other axis is max-reduced away in step 1).
    beta : float
        Inverse temperature; ``T = 1/beta`` multiplies both sums.
    nmat : int
        Stored window size; must equal ``X_w.shape[-1]`` (the centered map
        ``n~ = n - nmat//2`` is what defines S, so a mismatch would silently
        shift every frequency).

    Returns
    -------
    tail_est : float
        The estimated ``T sum_{tail} X`` (same units as ``T sum_S X``).
    tail_est_rel : float
        ``tail_est`` relative to the stored-window sum of the profile.
    unreliable : bool
        True when the asymptotic model does not describe the outer shell.
    """
    X_w = np.asarray(X_w)
    if X_w.ndim < 1 or X_w.shape[-1] != int(nmat):
        raise ValueError(
            "tail_estimate: X_w's last axis must be the fermionic Matsubara "
            "axis of length nmat = {}; got shape {}.".format(nmat, X_w.shape))
    nmat = int(nmat)
    if nmat < 8 or nmat % 2 != 0:
        raise ValueError(
            "tail_estimate: nmat must be an even number >= 8 (the centered "
            "grid map and the 3*nmat/8 outer shell both assume it); got {}."
            .format(nmat))
    beta = float(beta)
    if not np.isfinite(beta) or beta <= 0.0:
        raise ValueError(
            "tail_estimate: beta must be a positive finite number, got {!r}"
            .format(beta))
    T = 1.0 / beta

    # 1. reduction to the single scalar profile
    x = np.abs(X_w).reshape(-1, nmat).max(axis=0)
    n_t = np.arange(nmat) - nmat // 2
    w = (2 * n_t + 1) * np.pi / beta

    if not np.any(x > 0.0):
        # No signal: the tail of zero is zero, but there is nothing to fit, so
        # the number carries no information and is flagged as such.
        return 0.0, 0.0, True

    # 2./3. outer shell + asymptotic least squares.
    # The SECOND basis column is 1/|w|**3, not the signed 1/w**3: x is a
    # NONNEGATIVE magnitude profile (max_{k,pairs} |X|), even in w by
    # construction (the centered grid's negative and positive frequencies
    # both contribute magnitude, not signed values), while 1/w**3 is an ODD
    # function of w -- an inconsistent basis for an even target (review fix).
    # 1/|w|**3 matches x's even symmetry exactly.
    #
    # With BOTH columns now even (and strictly positive, and close in shape
    # over one outer shell -- 1/w^2 and 1/|w|^3 differ only by a slowly
    # varying extra 1/|w| factor there), the design is far more collinear
    # than the old (odd, orthogonal-ish) one: an UNCONSTRAINED lstsq can and
    # does return a NEGATIVE coefficient on real physical data even though
    # the true minimum-norm-consistent model has both a, b >= 0 (measured:
    # on the U=4/V=1 milestone fixture at Nmat=64, unconstrained lstsq gives
    # a<0 with |a| comparable to b, and forcing abs(a) afterwards -- the
    # previous code's move -- assembles a model that was never actually fit
    # to the data, one whose residual is ~40x WORSE than the true
    # nonnegative optimum). Since the model is only physically meaningful
    # for a, b >= 0 (a magnitude decays, its coefficient does not flip
    # sign), solve the constrained problem directly with nonnegative least
    # squares instead of unconstrained lstsq + abs(): it minimises the same
    # residual SUBJECT TO a, b >= 0, so it is never worse than any
    # nonnegative choice (including the abs()-of-lstsq guess) and is exact
    # for the normal case (an exact A/w**2 + B/|w|**3 profile with A, B > 0
    # has a, b >= 0 as its unconstrained optimum too, so NNLS reproduces it
    # exactly there as well).
    outer = np.abs(n_t) > 3.0 * nmat / 8.0
    design = np.stack([1.0 / w[outer] ** 2, 1.0 / np.abs(w[outer]) ** 3],
                      axis=1)
    coef, _ = nnls(design, x[outer])
    a = float(coef[0])
    b = float(coef[1])

    # 7. reliability of the model actually used below
    model = a / w[outer] ** 2 + b / np.abs(w[outer]) ** 3
    scale = float(np.linalg.norm(x[outer]))
    resid = float(np.linalg.norm(x[outer] - model))
    n_outer = int(np.count_nonzero(outer))
    unreliable = bool(scale == 0.0 or resid / scale > 0.2 or n_outer < 4)

    # 4./5. explicit summation of the fitted model over the closed tail window
    n_all = np.arange(-64 * nmat, 64 * nmat + 1)
    in_S = (n_all >= -(nmat // 2)) & (n_all <= nmat // 2 - 1)
    n_tail = n_all[~in_S]
    w_tail = (2 * n_tail + 1) * np.pi / beta
    tail_est = float(T * np.sum(a / w_tail ** 2 + b / np.abs(w_tail) ** 3))

    # 6. relative form
    denom = float(T * np.sum(x))
    tail_est_rel = tail_est / denom if denom > 0.0 else float("inf")
    return tail_est, tail_est_rel, unreliable


# ---------------------------------------------------------------------------
# bond_bubble's memory contract  --  read together with sc._bond_memory_estimate
# ---------------------------------------------------------------------------
#
# The resource preflight (``sc._bond_resource_preflight``) refuses a run whose
# estimated peak exceeds ``bond_memory_cap_gb``. That guard is only worth
# anything if the estimate is not an UNDERCOUNT, so :func:`bond_bubble` below
# is written to a fixed, documented working set: every temporary is released
# (``del`` / in-place rescaling) as soon as it is consumed, and the two counts
# here are the number of buffers that may be alive at the high-water mark.
# They are imported by ``sc._bond_memory_estimate`` so the two cannot silently
# desync; a new buffer in ``bond_bubble`` must bump the matching count (and
# ``tests/test_sc_bond.py`` measures the real peak against the budget).
#
# High-water mark: inside ``tau_to_boson(chi0_qt)`` in the channel-pair loop.
#
#   ``norb**4``-sized, on the ``(N_q, nmat)`` grid (BOND_BUBBLE_N4_BUFFERS):
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
#   as ``green_kw``):
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
BOND_BUBBLE_N4_BUFFERS = 3
BOND_BUBBLE_N2_BUFFERS = 3


def bond_bubble(green, bond_set, beta):
    """Compute the enlarged bond-channel bubble ``chi_bar(q; Delta r, Delta r')``.

    Implements spec S4.2 (Onari Eq. 14): the orbital contraction is
    *identical* to the existing ``chi0q`` bubble (mirrored verbatim from
    ``hwave.solver.rpa.RPA._calc_chi0q``: fermionic Matsubara frequency ->
    imaginary time via ``matsubara.fermion_to_tau``, k-space -> real space via
    ``backend.spatial_ifftn``, the ``G(r,tau) * G(-r,-tau) * sgn(tau)``
    product, then back to q-space/bosonic frequency via
    ``backend.spatial_fftn`` + ``matsubara.tau_to_boson``); this function only
    *adds* the bond-channel phase ``e^{i k . (Delta r_m - Delta r_m')}``
    inside the internal k-sum.

    That phase multiplies only the "bare-k" Green's function factor
    (``G_{l4 l2}(k)`` in the ``chi0q`` contraction below), which -- by the
    Fourier shift theorem -- is realized as a *rigid shift* of the same
    ``G(-r,-tau)`` real-space array by ``Delta r_m - Delta r_m'`` (an
    ``np.roll``) before the final spatial FFT. No new normalization or
    Matsubara convention is introduced; both are copied unchanged from
    ``_calc_chi0q``. All ``B**2`` channel-pair transforms are computed
    directly (the ``-q`` symmetry optimization is deferred; see spec S4.2).

    Unlike ``_calc_chi0q``, this function does **not** apply the
    ``green0_tail`` high-frequency Matsubara-tail correction -- that
    correction is deferred to a later task. The static output returned here
    is therefore the raw finite-``nmat`` sum only, not the tail-corrected
    value ``_calc_chi0q`` would produce for the ``m=m'=0`` block.

    Parameters
    ----------
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat)
        k-space, Matsubara-frequency Green's function in the ``sc.py``
        layout/normalization (see ``hwave.sc._calc_green`` /
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
        Bond-major index ``I = m*nd + idx`` with ``nd = norb**2`` and
        ``idx = l1*norb + l2`` (the existing ``sc.py`` orbital-pair
        convention). The ``m=m'=0`` block equals the ordinary ``chi0q``
        (same orbital contraction, phase = 1) within numerical tolerance.
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

    npair = norb * norb
    n_channels = bond_set.n_channels
    nd_enlarged = npair * n_channels

    # --- mirror rpa.py._calc_chi0q exactly (single source of truth for the
    # normalization / Matsubara convention: matsubara.fermion_to_tau and
    # backend.spatial_ifftn/spatial_fftn are the same functions rpa.py calls) ---
    green_kt = _ms.fermion_to_tau(green, axis=-1)
    green_rt = _bk.spatial_ifftn(green_kt, axes=(2, 3, 4))
    # Released as soon as it is consumed: only the buffers listed in
    # ``BOND_BUBBLE_N2_BUFFERS`` / ``BOND_BUBBLE_N4_BUFFERS`` (see the
    # module constants) may be alive at the high-water mark, because
    # ``sc._bond_memory_estimate``'s cap is computed from exactly that list.
    del green_kt

    # G(-r,-tau): identical roll(-1)+flip trick to rpa.py, over both the
    # spatial axes (2,3,4) and the imaginary-time axis (5).
    # (np.flip is a view of the np.roll result, so this is ONE buffer.)
    green_rev = np.flip(
        np.roll(green_rt, -1, axis=(2, 3, 4, 5)), axis=(2, 3, 4, 5)
    )

    sgn = np.full(nmat, -1.0, dtype=complex)
    sgn[0] = 1.0
    sgn = sgn.reshape(1, 1, 1, 1, 1, nmat)
    # G_fwd_sgn[l1,l3](r,tau) = G(r,tau)_{l1,l3} * sgn(tau)
    G_fwd_sgn = green_rt * sgn
    del green_rt          # green_rev does not alias it (np.roll copied)

    chi_bar = np.zeros((Nx, Ny, Nz, nd_enlarged, nd_enlarged), dtype=complex)

    delta_r = bond_set.delta_r
    static_index = nmat // 2
    for m in range(n_channels):
        dm = delta_r[m]
        for mp in range(n_channels):
            dmp = delta_r[mp]
            shift = (dm[0] - dmp[0], dm[1] - dmp[1], dm[2] - dmp[2])

            # Bond phase e^{i k.(Delta r_m - Delta r_m')} multiplying the
            # bare-k Green's function is, by the Fourier shift theorem, a
            # rigid shift of its real-space representation G(-r,-tau) by
            # (Delta r_m - Delta r_m') (an np.roll of the array holding
            # G(-r,-tau) -- see docstring).
            green_rev_shifted = np.roll(green_rev, shift=shift, axis=(2, 3, 4))

            # chi0_rt[l1,l2,l3,l4](r,tau)
            #     = G_fwd[l1,l3](r,tau) * sgn(tau) * G_rev_shifted[l4,l2](r,tau)
            # -- identical orbital contraction to rpa.py's general-mode chi0q
            # (a=l1,b=l2,c=l3,d=l4; no axis is contracted, purely an outer
            # product over orbital indices with an elementwise product over
            # (r,tau)).
            chi0_rt = np.einsum(
                'acxyzt,dbxyzt->abcdxyzt', G_fwd_sgn, green_rev_shifted
            )
            del green_rev_shifted

            chi0_qt = _bk.spatial_fftn(chi0_rt, axes=(4, 5, 6))
            del chi0_rt
            # tau_to_boson allocates one norb**4 temporary internally
            # (``arr * omg``) plus its ifft output, so at this point exactly
            # BOND_BUBBLE_N4_BUFFERS norb**4 arrays are alive; the final
            # rescaling is done IN PLACE so it does not add a fourth.
            chi0_qw = _ms.tau_to_boson(chi0_qt, axis=7)
            del chi0_qt
            chi0_qw *= (-1.0 / beta)

            # Static (zero bosonic frequency) slice only -- this function
            # does not carry a bosonic-frequency axis in its output.
            chi0_static = chi0_qw[..., static_index]  # (l1,l2,l3,l4,x,y,z)
            block = np.moveaxis(chi0_static, (4, 5, 6), (0, 1, 2)).reshape(
                Nx, Ny, Nz, npair, npair
            )

            chi_bar[:, :, :,
                    m * npair:(m + 1) * npair,
                    mp * npair:(mp + 1) * npair] = block
            # Released before the NEXT iteration allocates its own norb**4
            # buffers -- otherwise the previous chi0_qw would still be bound
            # and the real peak would be one norb**4 array above the budget.
            del chi0_static, block, chi0_qw

    return chi_bar


# ---------------------------------------------------------------------------
# bond_bubble_dynamic's memory contract  --  read together with
# sc._bond_memory_estimate(..., dynamic_nmat=...)
# ---------------------------------------------------------------------------
#
# ``bond_bubble_dynamic`` is a code-motion copy of ``bond_bubble`` (see that
# function's own buffer-accounting comment above) that keeps the full
# bosonic-frequency axis in ``chi0_qw`` instead of slicing out the static
# (zero-transfer) index. The per-iteration working set at the high-water
# mark -- the norb**4 buffers alive inside ``tau_to_boson(chi0_qt)`` and the
# norb**2 buffers alive across the channel-pair loop -- is therefore
# UNCHANGED from ``bond_bubble``: each is still one ``(N_q, nmat)``-grid
# array, because the frequency axis is already carried through the whole
# pipeline (``fermion_to_tau`` / ``spatial_ifftn`` / ``tau_to_boson``) even
# in the static function; only the STORED OUTPUT ``chi_bar`` grows, by a
# factor of ``nmat`` (``(N_q, nmat, ND, ND)`` instead of ``(N_q, ND, ND)``),
# because every bosonic-frequency slice is written into the output array
# instead of only the ``static_index`` one. These two constants are kept
# separate from ``BOND_BUBBLE_N4_BUFFERS`` / ``BOND_BUBBLE_N2_BUFFERS`` (even
# though their values are identical) so ``sc._bond_memory_estimate``'s
# dynamic branch cannot silently desync from this function's actual working
# set the way a single shared constant could.
BOND_BUBBLE_DYN_N4_BUFFERS = 3
BOND_BUBBLE_DYN_N2_BUFFERS = 3


def bond_bubble_dynamic(green, bond_set, beta):
    """Compute the frequency-resolved enlarged bond-channel bubble
    ``chi_bar(q, i nu_{j~}; Delta r, Delta r')``.

    Identical to :func:`bond_bubble` (spec S4.2 / Onari Eq. 14; see that
    function's docstring for the full derivation and normalization -- the
    orbital contraction, ``G(r,tau)*G(-r,-tau)*sgn(tau)`` construction, bond
    phase ``e^{i k.(Delta r_m - Delta r_m')}`` realized as a real-space
    ``np.roll``, and the fermion->tau->boson pipeline are all mirrored
    verbatim) except that this function does **not** collapse the bosonic
    Matsubara axis to its zero-transfer (static) slice: the full ``nmat``
    -length bosonic-frequency axis produced by ``matsubara.tau_to_boson`` is
    kept in the output. Zero bosonic transfer (``i nu_{0}``) sits at array
    index ``nmat // 2`` on that axis -- the same convention
    :func:`bond_bubble` uses internally for its ``static_index`` slice, so
    ``bond_bubble_dynamic(...)[..., nmat // 2, :, :]`` equals
    ``bond_bubble(...)`` exactly.

    Like :func:`bond_bubble`, this function does **not** apply the
    ``green0_tail`` high-frequency Matsubara-tail correction -- that
    correction is deferred to a later task (see spec S4.2 / Task 7's IR
    tail estimator). The values returned here are therefore the raw
    finite-``nmat`` sum only.

    Parameters
    ----------
    green : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat)
        k-space, Matsubara-frequency Green's function in the ``sc.py``
        layout/normalization (see ``hwave.sc._calc_green`` /
        ``hwave.sc._load_flex_green``): ``iomega_n = (2n+1-nmat)*pi/beta``.
    bond_set : ResolvedInteractionSet
        Bond-channel topology (``delta_r``, ``n_channels``) as returned by
        :func:`resolve_interactions`.
    beta : float
        Inverse temperature.

    Returns
    -------
    chi_bar_w : ndarray, shape (Nx, Ny, Nz, nmat, ND, ND), complex
        Frequency-resolved enlarged bubble; ``ND = norb**2 * B``. Bond-major
        index ``I = m*nd + idx`` with ``nd = norb**2`` and
        ``idx = l1*norb + l2`` (the existing ``sc.py`` orbital-pair
        convention), same as :func:`bond_bubble`. The bosonic-frequency axis
        (axis 3) has zero transfer at array index ``nmat // 2``; slicing
        that index reproduces :func:`bond_bubble`'s static output.
    """
    green = np.asarray(green)
    if green.ndim != 6:
        raise ValueError(
            "bond_bubble_dynamic: green must have shape (norb, norb, Nx, "
            "Ny, Nz, nmat); got ndim={} shape={}".format(
                green.ndim, green.shape)
        )
    norb, norb2, Nx, Ny, Nz, nmat = green.shape
    if norb2 != norb:
        raise ValueError(
            "bond_bubble_dynamic: green's two orbital axes must both equal "
            "norb; got shape {}".format(green.shape)
        )

    npair = norb * norb
    n_channels = bond_set.n_channels
    nd_enlarged = npair * n_channels

    # --- mirror rpa.py._calc_chi0q exactly (single source of truth for the
    # normalization / Matsubara convention: matsubara.fermion_to_tau and
    # backend.spatial_ifftn/spatial_fftn are the same functions rpa.py calls) ---
    green_kt = _ms.fermion_to_tau(green, axis=-1)
    green_rt = _bk.spatial_ifftn(green_kt, axes=(2, 3, 4))
    # Released as soon as it is consumed: only the buffers listed in
    # ``BOND_BUBBLE_DYN_N2_BUFFERS`` / ``BOND_BUBBLE_DYN_N4_BUFFERS`` (see
    # the module constants) may be alive at the high-water mark, because
    # ``sc._bond_memory_estimate``'s dynamic branch is computed from exactly
    # that list.
    del green_kt

    # G(-r,-tau): identical roll(-1)+flip trick to rpa.py, over both the
    # spatial axes (2,3,4) and the imaginary-time axis (5).
    # (np.flip is a view of the np.roll result, so this is ONE buffer.)
    green_rev = np.flip(
        np.roll(green_rt, -1, axis=(2, 3, 4, 5)), axis=(2, 3, 4, 5)
    )

    sgn = np.full(nmat, -1.0, dtype=complex)
    sgn[0] = 1.0
    sgn = sgn.reshape(1, 1, 1, 1, 1, nmat)
    # G_fwd_sgn[l1,l3](r,tau) = G(r,tau)_{l1,l3} * sgn(tau)
    G_fwd_sgn = green_rt * sgn
    del green_rt          # green_rev does not alias it (np.roll copied)

    chi_bar = np.zeros((Nx, Ny, Nz, nmat, nd_enlarged, nd_enlarged),
                       dtype=complex)

    delta_r = bond_set.delta_r
    for m in range(n_channels):
        dm = delta_r[m]
        for mp in range(n_channels):
            dmp = delta_r[mp]
            shift = (dm[0] - dmp[0], dm[1] - dmp[1], dm[2] - dmp[2])

            # Bond phase e^{i k.(Delta r_m - Delta r_m')} multiplying the
            # bare-k Green's function is, by the Fourier shift theorem, a
            # rigid shift of its real-space representation G(-r,-tau) by
            # (Delta r_m - Delta r_m') (an np.roll of the array holding
            # G(-r,-tau) -- see docstring).
            green_rev_shifted = np.roll(green_rev, shift=shift, axis=(2, 3, 4))

            # chi0_rt[l1,l2,l3,l4](r,tau)
            #     = G_fwd[l1,l3](r,tau) * sgn(tau) * G_rev_shifted[l4,l2](r,tau)
            # -- identical orbital contraction to rpa.py's general-mode chi0q
            # (a=l1,b=l2,c=l3,d=l4; no axis is contracted, purely an outer
            # product over orbital indices with an elementwise product over
            # (r,tau)).
            chi0_rt = np.einsum(
                'acxyzt,dbxyzt->abcdxyzt', G_fwd_sgn, green_rev_shifted
            )
            del green_rev_shifted

            chi0_qt = _bk.spatial_fftn(chi0_rt, axes=(4, 5, 6))
            del chi0_rt
            # tau_to_boson allocates one norb**4 temporary internally
            # (``arr * omg``) plus its ifft output, so at this point exactly
            # BOND_BUBBLE_DYN_N4_BUFFERS norb**4 arrays are alive; the final
            # rescaling is done IN PLACE so it does not add a fourth.
            chi0_qw = _ms.tau_to_boson(chi0_qt, axis=7)
            del chi0_qt
            chi0_qw *= (-1.0 / beta)

            # Full bosonic-frequency axis -- unlike bond_bubble, no
            # static_index slice: every transfer is kept in the output.
            block = np.moveaxis(chi0_qw, (4, 5, 6, 7), (0, 1, 2, 3)).reshape(
                Nx, Ny, Nz, nmat, npair, npair
            )

            chi_bar[:, :, :, :,
                    m * npair:(m + 1) * npair,
                    mp * npair:(mp + 1) * npair] = block
            # Released before the NEXT iteration allocates its own norb**4
            # buffers -- otherwise the previous chi0_qw would still be bound
            # and the real peak would be one norb**4 array above the budget.
            del block, chi0_qw

    return chi_bar


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

def pair_weight(green, beta):
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
    G2 = _g2_from_green(green, beta)          # (i,j,l,m,Nx,Ny,Nz)
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
