"""Declaration symmetrisation, single-sourced (#108 increment 2).

These are TABLE semantics: since issue #93 the file readers reject
non-Hermitian-closed input, so unclosed tables reach these reductions
only through reader-bypassing internal paths (and tests). A declaration
table may carry the same physical bond from both ends:
the two declarations of one coupling are ``(R, a, b)`` and
``(-R, b, a)``. For CoulombInter, Hund, Ising and PairLift the reversed
partner multiplies the SAME operator; for Exchange and PairHop it
multiplies the HERMITIAN CONJUGATE, which for a real Exchange coefficient
has the same effect. For every type except PairHop the table entering any
vertex is the mean

    T~[r, a, b] = (T[r, a, b] + T[-r, b, a]) / 2.

PairHop's partner entry is the HERMITIAN one, so its mean conjugates and
a Hermitian-closed complex declaration keeps its physical phase
(adjudicated in the #100/#105/#106 series; an asymmetric same-operator
declaration means the mean coefficient, and an antisymmetric one is an
identically zero Hamiltonian).

This module owns that content in BOTH representations it is consumed in:

* :func:`symmetrise_dense` -- real space, on a dense
  ``(nx, ny, nz, norb, norb)`` array. Used by the RPA ring's interaction
  assembly (``rpa._make_ham_inter``), which works on declaration tables
  before any Fourier transform. Dense on purpose: table keys may sit in
  a wrapped canonical form (``(n-1, 0, 0)`` for a ``-x`` bond; folded
  tables store canonicalized displacements), where a ``(-R)`` key lookup
  would silently miss the partner and halve the coefficient.
* :func:`symmetrise_k` -- momentum space, on the FFT'd
  ``(norb, norb, nx, ny, nz)`` tables. Used by the pair-space S/C
  builders in :mod:`hwave.sc` (RPA/FLEX/Eliashberg vertex side), whose
  pipeline transforms first and symmetrises second.

The two forms are the SAME reduction: the Fourier transform is linear,
so ``fft(symmetrise_dense(T)) == symmetrise_k(fft(T))`` exactly in
exact arithmetic (the r-space partner at ``(-R, b, a)`` maps to the
orbital transpose AT ``-q``; PairHop's conjugation at fixed ``q`` is the
Fourier image of {conjugate, ``R -> -R``}). That identity is pinned
executably in tests/test_declarations.py; the numerical results differ
only by float summation order, which is why BOTH forms exist here
instead of one delegating to the other through transforms -- each
consumer keeps its historical bit-exact pipeline while the content
lives in one module.

Relation to UHFk: ``uhfk.py`` stores the HERMITIAN mean for every type
(its table is contracted over both ordered orbital pairs, so a complex
same-operator coefficient p appears as p + conj(p) there). The two
conventions give the same physical Hamiltonian; the full
adjudication note lives in :func:`symmetrise_k`'s docstring.
"""

import numpy as np

from hwave.solver.kgrid import reverse_fft_axes

#: the seven documented two-index declaration types. The reduction rule
#: is a PHYSICS dispatch (same-operator vs Hermitian-partner, derived
#: from the operator definitions -- PairHop is the only Hermitian-partner
#: type; PairLift's endpoint reversal is the same operator because its
#: raising bilinears commute), so an unknown key must fail loudly rather
#: than silently receive the same-operator rule.
KNOWN_TYPES = frozenset(("CoulombIntra", "CoulombInter", "Hund",
                         "Exchange", "Ising", "PairLift", "PairHop"))


def symmetrise_dense(arr, hermitian=False):
    """Mean of a dense real-space table with its reversed-bond partner.

    Parameters
    ----------
    arr : ndarray
        Dense table, shape ``(nx, ny, nz, norb, norb)``, index 0 of each
        lattice axis = the on-site/origin cell (FFT order).
    hermitian : bool
        False (default): the same-operator mean
        ``0.5 * (T[r, a, b] + T[-r, b, a])`` -- every type but PairHop.
        True: the Hermitian-partner mean, which conjugates the reversed
        term -- PairHop.

    Returns
    -------
    ndarray
        The symmetrised table, same shape.
    """
    rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                       (0, 1, 2, 4, 3))
    if hermitian:
        rev = np.conjugate(rev)
    return 0.5 * (arr + rev)


def symmetrise_k(inter_k):
    """Momentum-space form of the same reduction, per interaction type.

    The same-operator partner of ``M[b, a](q)`` is the orbital-transposed
    entry AT ``-q`` (the interaction arrays carry ``e^{-iqR}`` phases since
    issue #232; the reduction is sign-agnostic) -- averaging with
    the same-q transpose instead corrupted every off-site interaction,
    collapsing a one-direction bond's ``V(q) = v e^{-iqR}`` to
    ``v cos(qR)`` (measured: the S/C entry vanished outright at
    ``q = pi/2``). The mean used here:

    * every type except PairHop: ``0.5 * (M(q) + M(-q)^T)``. For a
      both-ends declaration this mean is an identity; for an on-site
      complex Hermitian-closed declaration it drops the inert imaginary
      part; for an antisymmetric declaration -- an identically zero
      Hamiltonian -- it gives zero.
    * PairHop: ``0.5 * (M(q) + conj(M(q)^T))`` at the SAME q.
      Conjugation at fixed q is the Fourier image of
      {conjugate, ``R -> -R``} (the #105 derivation), so this is an
      identity for Hermitian-closed input on-site and off-site, and
      preserves the physical complex phase.

    Idempotent, so it is safe that both the all-q and per-q S/C builders
    apply it.

    Relation to UHFk: ``uhfk.py`` stores the HERMITIAN mean -- vba there
    is the CONJUGATED r-reversed transpose (``_make_ham_inter``) -- which
    keeps a complex Hermitian-closed coefficient p intact, while this
    reduction folds it to Re(p). The two tables serve different
    contractions and give the SAME physical Hamiltonian: UHFk sums its
    table over both ordered orbital pairs, so p + conj(p) appears there;
    the S/C builders use each slot exactly once, so the effective real
    coefficient must already be folded in. That the imaginary part of a
    Hermitian-closed same-operator declaration is physically inert was
    adjudicated against exact diagonalization in #105/#106; the slot
    content is #113. For real declarations the two conventions coincide
    identically; the file format itself admits complex entries, and for
    those the equivalence rests on the Hermitian-closure/redundancy
    assumptions stated above.

    Parameters
    ----------
    inter_k : dict
        ``{type: ndarray}`` with arrays of shape
        ``(norb, norb, Nx, Ny, Nz)``.

    Returns
    -------
    dict
        Symmetrised tables, same shapes.
    """
    out = {}
    for itype, M in inter_k.items():
        if itype not in KNOWN_TYPES:
            raise ValueError(
                "unknown interaction type '{}': the symmetrisation rule "
                "is a physics dispatch (same-operator vs Hermitian-"
                "partner) and cannot be guessed for an unrecognised "
                "type".format(itype))
        if itype == "PairHop":
            out[itype] = 0.5 * (M + np.conj(M.transpose(1, 0, 2, 3, 4)))
            continue
        Mrev = reverse_fft_axes(M, (2, 3, 4))
        out[itype] = 0.5 * (M + Mrev.transpose(1, 0, 2, 3, 4))
    return out


#: scale-aware tolerance of the Hermitian-closure validation, matching
#: the bond-channel resolver's reversal tolerance so the two checks
#: agree by construction (issue #93).
REVERSE_ATOL = 1e-10
REVERSE_RTOL = 1e-8


def validate_no_onsite_same_orbital_inter(itype, table, source=""):
    """Reject a nonzero on-site (``R = 0``) ``CoulombInter`` row of an orbital
    with itself.

    ``CoulombInter`` is the density interaction between two different
    orbitals or sites. Each declared row carries half of its operator
    ``V n n`` (the row and its partner together make one bond), so a row
    ``V_aa(0)`` of an orbital with itself is
    ``(V/2) n_a n_a = V n_{a up} n_{a dn} + (V/2) n_a``: a Hubbard term
    ``U = V`` plus a one-body level shift ``V/2``. The solvers do not read
    it that way -- the RPA/FLEX spin/charge vertex gives a ``CoulombInter``
    density row charge content only, and the Hartree-Fock decoupling of
    the same-spin part drops the level shift -- so it is refused; declare
    ``U = V`` in ``CoulombIntra`` and add ``V/2`` to the orbital's on-site
    energy in the Transfer file (or absorb it in the chemical potential).
    A row declared at zero is accepted (it carries no interaction).
    """
    if str(itype).lower() != "coulombinter":
        return
    for (irvec, orbvec), v in table.items():
        if (tuple(int(x) for x in irvec) == (0, 0, 0)
                and int(orbvec[0]) == int(orbvec[1]) and complex(v) != 0):
            raise ValueError(
                "CoulombInter{}: the on-site row of orbital {} with itself (R = (0, 0, 0), "
                "value V = {}) is the Hubbard term U = V plus a one-body level shift V/2 "
                "(each row carries half of V n n), which the solvers do not read as a "
                "CoulombInter density interaction. Declare U = V in CoulombIntra and add V/2 "
                "to the orbital's on-site energy in the Transfer file (or absorb it in the "
                "chemical potential); an on-site interaction between DIFFERENT orbitals "
                "stays in CoulombInter.".format(
                    " in {}".format(source) if source else "", int(orbvec[0]) + 1, v))


def validate_hermitian_closure(itype, table, source=""):
    """Reject a declaration table that is not Hermitian-closed (#93).

    The physical redundancy of the file format is
    ``X_ab(R) = conj(X_ba(-R))``: the two entries denote the same
    coupling (the Hermitian partner for PairHop; for the same-operator
    types a real closed pair, with an inert imaginary part folded by the
    symmetrised reading). A pair that DISAGREES beyond the scale-aware
    tolerance -- including a declared entry whose partner is absent or
    zero -- is not a different physical model; it is a typo or a
    convention mixup, and every solver used to do something different
    with it (silently project, read verbatim, or complete). Fail fast at
    read time instead, naming the file, the entry and both values.

    Keys are the LITERAL file displacements (validation runs before any
    folding, so the wrapped-canonical-key trap does not apply here).
    ``CoulombIntra`` is validated separately: the documented operator is
    on-site and same-orbital, so any entry with ``R != 0`` or ``a != b``
    is rejected outright (UHFk used to drop them silently).
    """
    def _fmt(irvec, a, b):
        return "R={} orbitals ({}, {})".format(tuple(irvec), a + 1, b + 1)

    where = " in {}".format(source) if source else ""
    def _reject_nonfinite(irvec, a, b, v):
        if not np.isfinite(complex(v)):
            raise ValueError(
                "{}{} declares a non-finite coefficient at {}: "
                "{}".format(itype, where, _fmt(irvec, a, b), v))

    if itype == "CoulombIntra":
        for (irvec, (a, b)), v in table.items():
            if tuple(irvec) != (0, 0, 0) or a != b:
                raise ValueError(
                    "CoulombIntra{} declares {} = {}: the documented "
                    "operator is on-site and same-orbital (R = 0, a == "
                    "b); off-site or inter-orbital entries belong in "
                    "CoulombInter".format(where, _fmt(irvec, a, b), v))
            _reject_nonfinite(irvec, a, b, v)
            # self-partner (R = -R = 0, a == b): Hermiticity means REAL
            if abs(np.imag(complex(v))) > (REVERSE_ATOL
                                           + REVERSE_RTOL * abs(v)):
                raise ValueError(
                    "CoulombIntra{} declares a complex diagonal at {}: "
                    "{} -- its Hermitian self-partner requires a real "
                    "value".format(where, _fmt(irvec, a, b), v))
        return
    for (irvec, (a, b)), v in table.items():
        _reject_nonfinite(irvec, a, b, v)
        rkey = ((-irvec[0], -irvec[1], -irvec[2]), (b, a))
        if rkey not in table:
            # membership FIRST (round-1 review): substituting 0 and
            # comparing let a declared |v| below the absolute tolerance
            # pass with no partner at all
            raise ValueError(
                "{}{} declares {} = {} with NO partner entry at R={} "
                "orbitals ({}, {}): declare both directions "
                "(X_ab(R) = conj(X_ba(-R)))".format(
                    itype, where, _fmt(irvec, a, b), v,
                    tuple(-x for x in irvec), b + 1, a + 1))
        partner = table[rkey]
        _reject_nonfinite(rkey[0], b, a, partner)
        tol = REVERSE_ATOL + REVERSE_RTOL * abs(v)
        if abs(v - np.conj(partner)) > tol:
            raise ValueError(
                "{}{} is not Hermitian-closed: {} = {} but the partner "
                "entry R={} orbitals ({}, {}) = {} (expected conj = "
                "{}). The two entries denote the SAME coupling "
                "(X_ab(R) = conj(X_ba(-R))); declare both directions "
                "with matching values -- a disagreement is read as a "
                "typo, not as a different model.".format(
                    itype, where, _fmt(irvec, a, b), v,
                    tuple(-x for x in irvec), b + 1, a + 1, partner,
                    np.conj(partner)))
