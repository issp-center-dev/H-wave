"""The pre-fold locality split of a two-body interaction table
(GitHub issue #181).

One owner for what Tier 1 (PR #182) did inside ``flex.py``: normalise
the aggregate ``Coulomb`` table, judge locality on the PRE-fold
declarations (``param_ham_orig`` under sublattice folding, failing
closed when it is missing), route every declared entry to an on-site or
an off-site part, fold each part separately, and keep the reader's own
folded, normalised table as the WHOLE table. That last point is what
keeps the pre-split summation order bit-identical for the previously
accepted class, and it is why the whole table is never a dict union of
the two folded parts (a full-period displacement folds onto the on-site
key of the same orbital pair) nor the fold of the pre-fold-normalised
table (the reader's fold-then-split classifies such an entry as
CoulombIntra, a split-then-fold as CoulombInter).

The split is DECLARATION-based: every declared off-site entry of every
two-body type is retained -- zero-valued entries, Exchange, PairHop and
PairLift included. Policy (reject / skip / warn) belongs to the callers:
the FLEX general path (``flex.py``) and the RPA bond-resolved
longitudinal gate (``rpa.py``) apply different policies to the same
split.
"""
from dataclasses import dataclass

import numpy as np

_TWO_BODY_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Exchange",
                   "PairHop", "Ising", "PairLift")

# The types whose OFF-site content has a q-only representation in the
# pair-space S/C matrices: the density (aa, bb) slots, i.e. the Hartree
# vertex V(q). Off-site Exchange/PairLift carry no longitudinal content
# (Tier 2, tests/test_offsite_exchange_ed_longitudinal.py); off-site
# PairHop has no local-pair form.
_OFFSITE_DENSITY_TYPES = ("CoulombInter", "Hund", "Ising")


@dataclass(frozen=True)
class LocalitySplit:
    """Result of :func:`split_locality`.

    ``whole_tbl`` is the reader's (folded, under SubShape) table,
    normalised; ``onsite_tbl``/``offsite_tbl`` are the two parts,
    folded separately when folding is active; ``offsite_prefold_tbl`` is
    the off-site part BEFORE folding (for diagnostics that must name the
    declared displacement); ``offsite_types`` lists the declared off-site
    types in declaration order.
    """
    whole_tbl: object
    onsite_tbl: dict
    offsite_tbl: dict
    offsite_prefold_tbl: dict
    offsite_types: tuple
    has_fold: bool


def normalize_coulomb(tbl_dict):
    """Split the aggregate ``Coulomb`` table into ``CoulombIntra`` (the
    r = 0 orbital-diagonal entries) and ``CoulombInter`` (everything
    else) via ``wan90.split_coulomb``, the shared decomposition.

    Without this an aggregate declaration silently produced a ZERO
    vertex in the FLEX general path (measured: chiq_s off by 1e-1
    against the identical explicit declaration). The container CLASS is
    preserved through ``.copy()``: the reader stores each table under
    the spelling the user wrote (``CaseInsensitiveDict``), and a dict
    comprehension here dropped that -- a non-canonically-spelled type
    then vanished from every lookup downstream (measured: 'hund' with an
    aggregate Coulomb declaration gave chiq_s identical to omitting Hund
    entirely, while 'Hund' differed from it by 4.7e-2).
    """
    if "Coulomb" not in tbl_dict:
        return tbl_dict
    if "CoulombIntra" in tbl_dict or "CoulombInter" in tbl_dict:
        raise ValueError(
            "Coulomb cannot be specified together with "
            "CoulombIntra or CoulombInter")
    from hwave.qlmsio import wan90
    intra, inter = wan90.split_coulomb(tbl_dict["Coulomb"])
    out = tbl_dict.copy()
    del out["Coulomb"]
    out["CoulombIntra"] = intra
    out["CoulombInter"] = inter
    return out


def split_locality(ham_info, lattice):
    """Split ``ham_info``'s interaction table by locality, judged on the
    PRE-fold declarations (see the module docstring)."""
    has_fold = tuple(getattr(lattice, "subshape", (1, 1, 1))) != (1, 1, 1)
    # Under folding the split must scan the PRE-fold table:
    # _init_interaction canonicalizes displacements modulo the folded
    # grid, and a folded dimension of size one maps every off-site
    # displacement to (0,0,0) -- e.g. CellShape=[4,4,1] with
    # SubShape=[4,1,1] turns +-x bonds into zero-displacement
    # inter-sublattice entries, which the folded table cannot
    # distinguish from genuinely on-site input.
    scan_ham = ham_info.param_ham
    if has_fold:
        # fail CLOSED: judging locality on the folded table is exactly
        # the bypass this split exists to prevent
        scan_ham = getattr(ham_info, "param_ham_orig", None)
        if scan_ham is None:
            raise ValueError(
                "sublattice folding is active but the pre-fold "
                "interaction table (param_ham_orig) is missing, so "
                "off-site input cannot be validated (the folded "
                "table canonicalizes displacements and can hide "
                "off-site entries).")
    scan_ham = normalize_coulomb(scan_ham)
    onsite_tbl, offsite_tbl = {}, {}
    for itype in _TWO_BODY_TYPES:
        if itype not in scan_ham:
            continue
        for (irvec, orbvec), v in scan_ham[itype].items():
            target = onsite_tbl if tuple(irvec) == (0, 0, 0) else offsite_tbl
            target.setdefault(itype, {})[(irvec, orbvec)] = v
    offsite_prefold_tbl = offsite_tbl
    if has_fold:
        # Fold each part separately: the folded table cannot tell a
        # folded bond from on-site input, so the split had to be made
        # before folding (the pre-fold-locality rule the RPA solver's
        # _append_inter_cross follows for the same reason).
        onsite_tbl = {t: ham_info._reshape_interaction(tbl, False)
                      for t, tbl in onsite_tbl.items()}
        offsite_tbl = {t: ham_info._reshape_interaction(tbl, False)
                       for t, tbl in offsite_tbl.items()}
    whole_tbl = normalize_coulomb(ham_info.param_ham)
    return LocalitySplit(whole_tbl, onsite_tbl, offsite_tbl,
                         offsite_prefold_tbl, tuple(offsite_tbl), has_fold)


def sc_matrices_from_split(split, offsite_types, norb, nx, ny, nz):
    """The Tier-1 locality-split pair-space S/C matrices, each
    ``(nx, ny, nz, norb**2, norb**2)``.

    The whole table keeps every type (its on-site cross/antidiag content
    is read from there); the off-site part handed to the builder is
    FILTERED to ``offsite_types`` (the density types it accepts). The k
    grid is ``linspace(0, 2pi, n, endpoint=False)`` per axis, C-order
    flattened -- the same points, order and flattening as chi0's spatial
    FFT axis, verified by the element-complete equivalence with the RPA
    ring for off-site (q-dependent) entries.
    """
    from hwave.sc import _build_interaction_k
    from hwave.solver._sc_matrices_myo import build_sc_matrices_locality_split
    kx = np.linspace(0, 2.0 * np.pi, nx, endpoint=False)
    ky = np.linspace(0, 2.0 * np.pi, ny, endpoint=False)
    kz = np.linspace(0, 2.0 * np.pi, nz, endpoint=False)
    offsite_density = {t: tbl for t, tbl in split.offsite_tbl.items()
                       if t in offsite_types}
    inter_k = _build_interaction_k(kx, ky, kz, split.whole_tbl, norb)
    inter_k_onsite = _build_interaction_k(kx, ky, kz, split.onsite_tbl, norb)
    inter_k_offsite = _build_interaction_k(kx, ky, kz, offsite_density, norb)
    return build_sc_matrices_locality_split(
        inter_k, inter_k_onsite, inter_k_offsite, norb, nx, ny, nz)
