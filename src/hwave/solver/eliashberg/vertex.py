"""S/C matrices, pairing vertices and the k-space interaction tables of hwave_sc.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 1); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import numpy as np

from hwave.solver.vertex_table import sc_coefficients
from hwave.solver.kgrid import reverse_fft_axes
from hwave.solver.declarations import symmetrise_k
from hwave.solver import bond_channels

logger = logging.getLogger("hwave_sc")


#: ``sc_vertex_version`` from which a FLEX artifact puts every off-site row
#: on the bond of its Transfer row (issue #232).
PLACEMENT_VERSION = 3


#: The two-body types whose off-site placement reaches a FLEX artifact
#: (PairLift has no particle-hole content, see _read_flex_chi_raw).
_PLACEMENT_TYPES = ("Coulomb", "CoulombInter", "Hund", "Ising", "Exchange",
                    "PairHop")


#: Interaction terms whose Kuroki S/C matrices have entries OUTSIDE the
#: density-pair block, keyed by the _build_sc_matrices_all_q case that puts
#: them there.  Case 2 is S/C[(a,b),(a,b)] and case 4 is S/C[(a,b),(b,a)], both
#: with a != b; a reduced run never computes chi on those pair indices,
#: so the corresponding fluctuation dressing is missing from the vertex.
#: Types whose vertex is PARTIALLY represented by a reduced chi: their
#: density-slot content is dressed, their cross-slot content is not.
_REDUCED_FLEX_PARTIAL = ("CoulombInter", "Hund", "Ising")


#: Types with NO density-diagonal vertex content at all
#: (hwave.solver.vertex_table): with a reduced chi nothing of them is
#: dressed, and since the #120 policy no reduced FLEX or RPA run
#: can even be produced with them -- such input is stale or mismatched,
#: and is REJECTED rather than warned about.
_REDUCED_FLEX_REJECTED = ("Exchange", "PairHop")


#: the types whose cross / antidiag (exchange) content the q-only vertex of the
#: non-bond paths cannot hold for a row with R != 0 (issue #239)
_NONBOND_FOCK_TYPES = ("CoulombInter", "Hund", "Ising")


def _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb):
    """Build interaction matrices in k-space.

    Supports CoulombIntra (U), CoulombInter (U'), Hund (J), and Exchange (J').

    Parameters
    ----------
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    interactions : dict
        Interaction parameters keyed by type.
    norb : int
        Number of orbitals.

    Returns
    -------
    inter_k : dict
        Dictionary of interactions in k-space.
        Keys: "CoulombIntra", "CoulombInter", "Hund", "Exchange".
        Values: ndarray of shape (norb, norb, Nx, Ny, Nz).
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    kx_mesh, ky_mesh, kz_mesh = np.meshgrid(
        kx_array, ky_array, kz_array, indexing='ij'
    )

    def _to_k(value_r, transpose=True):
        # The ORBITAL PAIR is stored transposed -- an entry (R, (a, b))
        # lands at [b, a] -- with the phase e^{-iqR} (issue #232; the
        # one-body Hamiltonian keeps e^{+ikR}, see the end of this comment).
        #
        # The interaction is a four-index object, and its MATRIX form carries a
        # pair-index transpose that the one-body Hamiltonian does not. H-wave's
        # own paper (arXiv:2308.00324) makes this explicit: Eq.(12) defines the
        # tensor W_ij^{aa'bb'} with the first pair at site i and the second at
        # site j, Eq.(16) puts the first pair at momentum +q, and Eq.(21) then
        # defines the matrix used in the RPA equation as
        #     [W(q)]^{ab} = W_q^{ba}
        # so for a density-density term W^{aabb} = V_ab(R) the matrix element is
        # [W(q)]^{(aa),(bb)} = V_ba(q), i.e. the matrix is V(q)^T.
        #
        # `rpa.py::_make_ham_inter` has always built that (it stores the entry's
        # first orbital in the SECOND pair slot); this builder did not, so the
        # susceptibility was solved with one orientation and the pairing vertex
        # assembled from the other (issue #96).
        #
        # Both consumers need the transpose, confirmed by exact
        # diagonalization on a bond set with V_ab(R) != V_ba(R):
        #   * the RPA ladder [I + chi0 W]^-1 chi0 -- solved here too, in
        #     _compute_vertices_simple -- selects V^T with a residual that
        #     scales linearly in V (pure O(V^2) truncation), against 98% for V;
        #   * the bare pair-scattering amplitude
        #     <k' a up, -k' b down| H_int |k a up, -k b down>, which is exact at
        #     first order, matches V^T to 6e-16, against 99% for V.
        #
        # NOTE: `_build_hamiltonian_k` is a one-body object and correctly keeps
        # the plain [a, b] placement. Do not "align" this builder to it.
        #
        # PairHop is EXCLUDED (transpose=False below), and unlike the rest that
        # is now a MEASURED result rather than a gap. It is the one type
        # `rpa.py` does not place through `_append_inter`: `_append_pairhop`
        # uses the slots (b, a, a, b) rather than the density-density
        # (b, b, a, a), and in the S/C matrices it lands in the pair-ANTIdiagonal
        # Case 4 rather than Case 1/3, so the density-density determination
        # above says nothing about it.
        #
        # Issue #100 settled it end to end (TestPairHopEndToEnd in
        # tests/test_vertex_orientation.py), with no index reasoning anywhere
        # between the input file and the answer: feed this builder an on-site
        # PairHop matrix, build S from it, form the solver's own linear-order
        # response -chi0 . S . chi0 with chi0 from `rpa.py`'s kernel and an
        # EXACT Green function, and compare against the cross-spin response of
        # the same Hamiltonian by exact diagonalization. The untransposed
        # placement -- i.e. S[(a,b),(b,a)] = P[a,b] -- reproduces it, and the
        # test locates the residual by checking chi0 against the exact bubble
        # in the same run, so what is left is the imaginary-time grid. Refining
        # that grid off-line, 192 -> 384 -> 768 points, takes the residual
        # 2.3e-3 -> 1.1e-3 -> 5.5e-4 (the test itself runs only the first);
        # transposing PairHop sits at 0.55 and does not move.
        #
        # Beware the plausible-looking derivation that says otherwise. The
        # vertex sits BETWEEN two chi0 factors, so its row index lives in
        # chi0's COLUMN convention and its column index in chi0's ROW
        # convention. chi0's row pair is reversed relative to its bilinear and
        # its column pair is not, so the reversal lands on the vertex's COLUMN
        # side. Applying it to the row side instead gives the transpose, and is
        # wrong.
        val_k = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        for (irvec, orbvec), value in value_r.items():
            orb1, orb2 = orbvec
            Rx, Ry, Rz = irvec
            i1, i2 = (orb2, orb1) if transpose else (orb1, orb2)
            # A row (R, a, b) is the bond of a Transfer row (R, a, b): a in
            # the original cell, b in the cell displaced by R (Wannier90 /
            # RESPACK, issue #232). With the slot convention above that
            # operator carries e^{-ik.R}; e^{+ik.R} put it on the opposite
            # bond.
            val_k[i1, i2, :, :, :] += value * np.exp(
                -1j * (kx_mesh * Rx + ky_mesh * Ry + kz_mesh * Rz)
            )
        return val_k

    inter_k = {}
    for itype in ["CoulombIntra", "CoulombInter", "Hund", "Exchange",
                  "Ising", "PairLift", "PairHop"]:
        if itype in interactions:
            inter_k[itype] = _to_k(interactions[itype],
                                   transpose=(itype != "PairHop"))

    return inter_k


def _symmetrise_interactions_k(inter_k):
    """Reduce each interaction to its physical symmetric coefficient.

    Thin delegation: the reduction and its full derivation -- the
    momentum-space form of the reversed-bond partner, the PairHop
    Hermitian rule, and the UHFk-relation adjudication note -- are
    single-sourced in :func:`hwave.solver.declarations.symmetrise_k`
    (#108). Idempotent, so it is safe that both the all-q and per-q
    builders apply it.
    """
    return symmetrise_k(inter_k)


def _accumulate_coeff(dst, coeff, value):
    """dst += coeff * value with IEEE-preserving forms.

    NumPy evaluates ``1.0 * (Inf+0j)`` as ``Inf+NaNj`` (a full complex
    multiply), while a direct add preserves the zero imaginary part, so
    the +-1 coefficients of the vertex table must use direct add/subtract
    to keep the pre-table behavior for non-finite couplings. Zero
    coefficients suppress the contribution entirely (0.0 * Inf is NaN).
    ``dst`` must be an ndarray or writeable view; the update is in place.
    """
    if coeff == 0.0:
        return
    if coeff == 1.0:
        dst += value
    elif coeff == -1.0:
        dst -= value
    else:
        dst += coeff * value


class SiteSplitInteractionK(dict):
    """A k-space interaction table (as returned by
    :func:`_build_interaction_k`) that also carries the locality split
    ``(inter_k_onsite, inter_k_offsite)``: the ``R != 0`` rows (issue #239)
    and, with ``[mode.param] orbital_site`` (#234), the ``R = 0`` rows between
    orbitals of different sites form the off-site part. The S/C builders
    pick the split up from the table itself, so every hwave_sc vertex path
    that receives the table honours the declaration."""

    def __init__(self, data, locality_split):
        super().__init__(data)
        self.locality_split = locality_split


def _offsite_rows(interactions, types, site=None):
    """``{type: [(irvec, orbvec, value), ...]}`` of the nonzero off-site rows:
    ``R != 0``, and with ``orbital_site`` labels ``site`` also the ``R = 0``
    rows between orbitals of different sites (``sites.is_onsite``)."""
    from hwave.solver.sites import is_onsite
    out = {}
    for t in types:
        tbl = interactions.get(t) if hasattr(interactions, "get") else None
        rows = [(irvec, orbvec, v) for (irvec, orbvec), v in (tbl or {}).items()
                if v != 0 and not is_onsite(irvec, orbvec, site)]
        if rows:
            out[t] = rows
    return out


def _nonbond_site(input_dict, norb):
    """The ``[mode.param] orbital_site`` labels (or ``None``)."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver.sites import parse_orbital_site
    param = CaseInsensitiveDict(input_dict.get("mode", {}).get("param", {}) or {})
    return parse_orbital_site(param.get("orbital_site"), norb)


def _build_interaction_k_for_sites(kx_array, ky_array, kz_array, interactions, norb,
                                   input_dict):
    """:func:`_build_interaction_k` for the hwave_sc entries, with the
    locality split of the FLEX general path.

    An ``R != 0`` row keeps its Hartree (density) part only: its exchange
    crossing is left out, as in the RPA ring and the FLEX general path (issue
    #239; it used to be collapsed into the Fock slot of the q-only vertex).
    ``[mode.param] orbital_site`` (#234) makes an ``R = 0`` row between
    orbitals of different sites off-site in the same sense. Without such
    rows the plain table is returned (unchanged results)."""
    from requests.structures import CaseInsensitiveDict
    from hwave.solver.sites import (parse_orbital_site, warn_if_undeclared,
                                    cross_site_rows)
    param = CaseInsensitiveDict(input_dict.get("mode", {}).get("param", {}) or {})
    site = parse_orbital_site(param.get("orbital_site"), norb)
    warn_if_undeclared(interactions, site, "hwave_sc")
    inter_k = _build_interaction_k(kx_array, ky_array, kz_array, interactions, norb)
    rows = cross_site_rows(interactions, site)
    # (the non-bond entries refuse such Exchange / PairHop rows beforehand:
    # _check_nonbond_offsite)
    for t in ("Exchange", "PairHop", "PairLift"):
        if t in rows:
            logger.warning(
                "hwave_sc: %d R = 0 %s row(s) join orbitals of different sites "
                "([mode.param] orbital_site); like off-site %s they have no "
                "density part and are left out of the vertex.", len(rows[t]), t, t)
    # every R != 0 row is off-site: the density part of CoulombInter / Hund /
    # Ising enters through the off-site table, nothing of it through the
    # cross / antidiag slots (which read the on-site table)
    offsite_keys = {t: [key for key in tbl if tuple(int(x) for x in key[0]) != (0, 0, 0)]
                    for t, tbl in interactions.items()}
    offsite_keys = {t: keys for t, keys in offsite_keys.items() if keys}
    if not rows and not offsite_keys:
        return inter_k
    drop = {t: set(rows.get(t, ())) | set(offsite_keys.get(t, ()))
            for t in set(rows) | set(offsite_keys)}
    offsite = {}
    for t, keys in drop.items():
        if t in _NONBOND_FOCK_TYPES:
            # in the source order, as offsite.split_locality builds it: the
            # Fourier sums then accumulate in the same order
            offsite[t] = {key: v for key, v in interactions[t].items() if key in keys}
    onsite = {}
    for t, tbl in interactions.items():
        kept = {key: v for key, v in tbl.items() if key not in drop.get(t, ())}
        if kept:
            onsite[t] = kept
    split = (_build_interaction_k(kx_array, ky_array, kz_array, onsite, norb),
             _build_interaction_k(kx_array, ky_array, kz_array, offsite, norb))
    if rows:
        logger.info("hwave_sc: orbital_site makes %s R = 0 row(s) inter-site (Hartree part "
                    "only).", ", ".join("{} {}".format(len(v), t) for t, v in rows.items()))
    return SiteSplitInteractionK(inter_k, split)


def _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz,
                             _presymmetrised=False, locality_split=None):
    """Build spin (S) and charge (C) interaction matrices for all q-points at once.

    Follows Kuroki et al., Eq.(5) in arXiv:0902.3691:
        S_{l1l2,l3l4}, C_{l1l2,l3l4} for multi-orbital systems.

    Parameters
    ----------
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.
    locality_split : None or (dict, dict), keyword
        Default None: the on-site Kanamori slot map is applied to
        ``inter_k`` as a whole (every existing caller). When given, the
        pair ``(inter_k_onsite, inter_k_offsite)`` -- the k-space forms of
        the ``R == 0`` and ``R != 0`` PRE-fold declarations, separately --
        makes the slot map locality-aware (#181 Tier 1, the general FLEX
        path). ``inter_k`` stays the authoritative whole table (the
        reader's own, folded under SubShape); the two parts are
        locality-ROUTING inputs and need not form a per-type
        decomposition of it (under folding a full-period displacement
        folds onto an on-site key, and an aggregate ``Coulomb`` entry can
        then be CoulombIntra in the whole table while the pre-fold split
        holds it as CoulombInter): the cross (Case 2) and antidiag (Case 4)
        families, whose particle-hole pair is NON-local for ``R != 0``,
        are built from the on-site part only, and the same-orbital
        density elements of Hund / Ising (Case 3 at ``l1 == l3``, where
        an orbital has no ON-site Hund / Ising coupling with itself) take
        the off-site part, which is a physical density-density term.
        Everything else keeps reading ``inter_k`` as a whole, so the
        summation order of every element that existed without the split
        is unchanged (bit-identical output for on-site input and for
        same-orbital off-site CoulombInter, the class accepted before).

    Returns
    -------
    S_all : ndarray
        Spin interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    C_all : ndarray
        Charge interaction matrices, shape (Nx, Ny, Nz, norb^2, norb^2).
    """
    if locality_split is None:
        # an hwave_sc table built with [mode.param] orbital_site (#234)
        # carries its own split (see _build_interaction_k_for_sites)
        locality_split = getattr(inter_k, "locality_split", None)
    nd = norb * norb
    S_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
    C_all = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

    # _presymmetrised is set by the per-q wrapper, which has already
    # symmetrised on the FULL grid: the -q partner of an off-site entry is
    # unreachable from a single-point slice, so re-symmetrising here would
    # average with the wrong (same-q) partner and corrupt off-site input.
    if not _presymmetrised:
        inter_k = _symmetrise_interactions_k(inter_k)
    if locality_split is None:
        inter_k_onsite = inter_k
        inter_k_offsite = None
    else:
        inter_k_onsite, inter_k_offsite = locality_split
        if not _presymmetrised:
            inter_k_onsite = _symmetrise_interactions_k(inter_k_onsite)
            inter_k_offsite = _symmetrise_interactions_k(inter_k_offsite)

    def _get(itype, source=inter_k):
        if source is not None and itype in source:
            return source[itype]  # (norb, norb, Nx, Ny, Nz)
        return None

    U_mat = _get("CoulombIntra")
    Up_mat = _get("CoulombInter")
    J_mat = _get("Hund")
    I_mat = _get("Ising")
    # the families whose pair is non-local off-site read the on-site part
    # (Exchange and PairHop have cross / antidiag content only)
    Up_on = _get("CoulombInter", inter_k_onsite)
    J_on = _get("Hund", inter_k_onsite)
    Jp_on = _get("Exchange", inter_k_onsite)
    I_on = _get("Ising", inter_k_onsite)
    PH_on = _get("PairHop", inter_k_onsite)
    # the same-orbital density elements of Hund / Ising read the off-site
    # part (None without a split: the on-site gate below then applies)
    J_off = _get("Hund", inter_k_offsite)
    I_off = _get("Ising", inter_k_offsite)

    # Build using precomputed index arrays to avoid Python loops
    # for small norb (1-3), the loop overhead is negligible;
    # for larger norb the vectorized approach helps
    l1_arr, l2_arr, l3_arr, l4_arr = np.meshgrid(
        np.arange(norb), np.arange(norb), np.arange(norb), np.arange(norb),
        indexing='ij')
    l1f, l2f, l3f, l4f = l1_arr.ravel(), l2_arr.ravel(), l3_arr.ravel(), l4_arr.ravel()
    idx12 = l1f * norb + l2f
    idx34 = l3f * norb + l4f

    # Case 1: l1 == l2 == l3 == l4
    # Accumulate rather than assign: this element is also reached by Case 3
    # below (which now includes l1 == l3), where an inter-site same-orbital
    # CoulombInter contributes 2 V_aa(q) to the charge channel.
    mask1 = (l1f == l2f) & (l2f == l3f) & (l3f == l4f)
    if U_mat is not None and np.any(mask1):
        sU, cU = sc_coefficients("CoulombIntra", "diag")
        for i in np.where(mask1)[0]:
            _l = l1f[i]
            _accumulate_coeff(S_all[:, :, :, idx12[i], idx34[i]],
                              sU, U_mat[_l, _l])
            _accumulate_coeff(C_all[:, :, :, idx12[i], idx34[i]],
                              cU, U_mat[_l, _l])

    # Case 2: l1==l3, l2==l4, l1!=l2
    # Coefficients come from the single adjudicated table
    # (hwave.solver.vertex_table). Signs, slots and per-type factors were
    # established against exact diagonalization in #113 (Ising sign, the
    # previously missing Hund S term, Exchange moved here from the
    # pair-antidiagonal Case 4); the SU(2) Kanamori consistency check --
    # Hund + Exchange at equal J giving S(ab,ab) without J and
    # C(ab,ab) = -U' + 2J -- is pinned in the tests.
    mask2 = (l1f == l3f) & (l2f == l4f) & (l1f != l2f)
    cross_terms = [(Up_on, "CoulombInter"), (I_on, "Ising"),
                   (J_on, "Hund"), (Jp_on, "Exchange")]
    for i in np.where(mask2)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        for mat, itype in cross_terms:
            if mat is not None:
                sco, cco = sc_coefficients(itype, "cross")
                _accumulate_coeff(s_q, sco, mat[_l1, _l2])
                _accumulate_coeff(c_q, cco, mat[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    # Case 3: l1==l2, l3==l4 -- INCLUDING l1 == l3, but for CoulombInter ONLY.
    # The l1 != l3 exclusion dropped the inter-site same-orbital CoulombInter
    # from the charge diagonal C[(a,a),(a,a)], which must be U_a + 2 V_aa(q)
    # (issue #95); the simple two-index formulation used by chi0q_mode="load"
    # builds exactly that (`Wc = U_k + 2 V_k`, _compute_vertices_simple).
    # Case 1 above writes U_a into the same element, so both accumulate.
    # Hund and Ising stay restricted to l1 != l3 for ON-site input: an
    # orbital has no on-site Hund or Ising coupling with itself, and letting
    # a stray diagonal entry through here would silently move S as well.
    # With a locality split the l1 == l3 element takes the OFF-site part
    # only (an inter-site same-orbital Hund / Ising is a physical
    # density-density term, #181 Tier 1).
    mask3 = (l1f == l2f) & (l3f == l4f)
    for i in np.where(mask3)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l3 = l1f[i], l3f[i]
        if _l1 != _l3:
            density_terms = ((J_mat, "Hund"), (I_mat, "Ising"))
        else:
            density_terms = ((J_off, "Hund"), (I_off, "Ising"))
        for mat, itype in density_terms:
            if mat is not None:
                sco, cco = sc_coefficients(itype, "density")
                _accumulate_coeff(s_q, sco, mat[_l1, _l3])
                _accumulate_coeff(c_q, cco, mat[_l1, _l3])
        if Up_mat is not None:
            sco, cco = sc_coefficients("CoulombInter", "density")
            _accumulate_coeff(s_q, sco, Up_mat[_l1, _l3])
            _accumulate_coeff(c_q, cco, Up_mat[_l1, _l3])
        S_all[:, :, :, idx12[i], idx34[i]] += s_q
        C_all[:, :, :, idx12[i], idx34[i]] += c_q

    # Case 4: l1==l4, l2==l3, l1!=l2
    mask4 = (l1f == l4f) & (l2f == l3f) & (l1f != l2f)
    for i in np.where(mask4)[0]:
        s_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        c_q = np.zeros((Nx, Ny, Nz), dtype=complex)
        _l1, _l2 = l1f[i], l2f[i]
        # Exchange used to sit here; exact diagonalization (issue #113) puts
        # its vertex on the pair-DIAGONAL slot family (Case 2), and end to end
        # the antidiagonal placement produced the right magnitude at the wrong
        # slots in both channels. Only PairHop belongs here (#100/#102).
        if PH_on is not None:
            sco, cco = sc_coefficients("PairHop", "antidiag")
            _accumulate_coeff(s_q, sco, PH_on[_l1, _l2])
            _accumulate_coeff(c_q, cco, PH_on[_l1, _l2])
        S_all[:, :, :, idx12[i], idx34[i]] = s_q
        C_all[:, :, :, idx12[i], idx34[i]] = c_q

    return S_all, C_all


def _build_sc_matrices(inter_k, norb, ix, iy, iz):
    """Spin (S) and charge (C) matrices at a single q-point.

    Delegates to :func:`_build_sc_matrices_all_q` and slices, so there is ONE
    implementation of the S/C content. The previous hand-maintained copy had
    already drifted from the all-q builder once; after the issue #113 vertex
    corrections a second parallel copy would be a liability.
    """
    # Symmetrise on the FULL grid first -- the same-operator partner of an
    # off-site entry lives at -q, so slicing before symmetrising would discard
    # it -- then slice the (small, norb^2 per point) interaction arrays down to
    # the requested q and delegate on a 1x1x1 grid, so the (nd^2 per point) S/C
    # matrices are never built for the whole grid. Indexing goes through
    # np.arange so numpy negative indices keep their usual meaning and
    # out-of-range indices raise instead of returning empty slices.
    def _slice(table):
        out = {}
        for t, M in _symmetrise_interactions_k(table).items():
            jx = np.arange(M.shape[2])[ix]
            jy = np.arange(M.shape[3])[iy]
            jz = np.arange(M.shape[4])[iz]
            out[t] = np.ascontiguousarray(M[:, :, jx:jx+1, jy:jy+1, jz:jz+1])
        return out

    inter_1 = _slice(inter_k)
    split = getattr(inter_k, "locality_split", None)
    split_1 = None if split is None else (_slice(split[0]), _slice(split[1]))
    S_all, C_all = _build_sc_matrices_all_q(inter_1, norb, 1, 1, 1,
                                            _presymmetrised=True,
                                            locality_split=split_1)
    return S_all[0, 0, 0], C_all[0, 0, 0]


def _declarations_partner_closed(interactions, Nx, Ny, Nz, norb):
    """True iff every CoulombIntra/CoulombInter declaration table is
    EXACTLY closed under the reversed-bond partner (R,a,b) <-> (-R,b,a).

    Decided on the RAW tables, before any exponential transform, with
    bitwise comparison: a both-ends declaration carries identical
    literals, and 0.5*(x + x) == x exactly in IEEE, so closure detection
    is exact -- while the k-space arrays of a closed table are only
    ALGEBRAICALLY symmetric (the +R and -R exponentials are evaluated
    independently and differ by roundoff), which is why the decision
    cannot be made after the transform (PR #129 round 5: re-averaging a
    closed configuration changed saved artifacts at the 1e-18 level).
    Non-finite input reads as not closed (explicitly, before any
    comparison), which fails toward symmetrising. The closure test
    compares the table DIRECTLY with its reversed/transposed partner --
    no mean is formed, so a closed coefficient above float64.max / 2
    cannot overflow into a false 'open' (PR #129 round 6 reproduced
    exactly that with the earlier mean-based comparison)."""
    for itype in ("CoulombIntra", "CoulombInter"):
        tbl = interactions.get(itype) if interactions else None
        if not tbl:
            continue
        arr = np.zeros((Nx, Ny, Nz, norb, norb), dtype=complex)
        for (irvec, orbvec), v in tbl.items():
            arr[(irvec[0] % Nx, irvec[1] % Ny, irvec[2] % Nz,
                 *orbvec)] += v
        if not np.all(np.isfinite(arr)):
            return False
        rev = np.transpose(reverse_fft_axes(arr, (0, 1, 2)),
                           (0, 1, 2, 4, 3))
        if not (np.array_equal(arr.real, rev.real)
                and np.array_equal(arr.imag, rev.imag)):
            return False
    return True


def _compute_vertices(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                      pairing_type="singlet", static_index=None,
                      declarations_closed=False):
    """Compute effective pairing interaction V(q).

    Supports two modes:
    - Simple mode (only CoulombIntra/Inter, 2-index chi0q):
      Uses Wc=U+2V, Ws=-U formulation.
    - General mode (with Hund/Exchange, or 4-index chi0q):
      Uses S,C matrices from Kuroki et al.

    When 4-index chi0q is provided (8D array), general mode is always used
    because the full orbital tensor structure must be preserved.

    Pairing types:
    - singlet: V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)
    - triplet: V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S)

    Parameters
    ----------
    chi0q : ndarray
        Bare susceptibility.
        - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
        - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb, Nx, Ny, Nz, nmat : int
        System parameters.
    pairing_type : str
        "singlet" or "triplet". Default "singlet".

    Returns
    -------
    Result depends on mode:
        Simple mode: tuple (Pc_q, Ps_q), each shape (norb, norb, Nx, Ny, Nz).
        General mode: Vs_q, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    # NOTE (#167): this list is the Eliashberg-side vertex-layout choice and
    # is deliberately NOT derived from scheme.CAPABILITIES -- it differs
    # from the chi0q_tensor='auto' rule (no CoulombInter/Coulomb) and is
    # behaviour-frozen. Unifying it is a recorded follow-up decision.
    has_interorbital_vertex = any(k in inter_k for k in
                                  ["Hund", "Exchange", "Ising", "PairHop"])
    chi0q_is_4index = (chi0q.ndim == 8)

    # PairLift does not enter the spin/charge (particle-hole) pairing vertex in
    # either mode (its S=C=0 contribution, verified against the full 4-index
    # RPA), so a configured PairLift term is silently inert. Warn the user.
    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    if chi0q_is_4index or has_interorbital_vertex:
        # General mode: 4-index S,C matrices
        # Required when 4-index chi0q is available or Hund/Exchange present
        return _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                         pairing_type=pairing_type,
                                         static_index=static_index)
    else:
        # Simple mode: backward compatible with original implementation
        # Only used for 2-index chi0q without Hund/Exchange
        if "CoulombInter" in inter_k and norb > 1:
            logger.warning(
                "Multi-orbital CoulombInter in the simple (2-index) vertex mode "
                "uses the reduced density-density approximation: inter-orbital "
                "cross-channel contributions are dropped. Provide a 4-index "
                "chi0q (general mode) for the full Kuroki S/C treatment.")
        return _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                                        pairing_type=pairing_type,
                                        static_index=static_index,
                                        declarations_closed=declarations_closed)


def _compute_vertices_simple(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                             pairing_type="singlet", static_index=None,
                             declarations_closed=False):
    """Compute vertices using simple Wc=U+2V, Ws=-U formulation.

    For singlet:
        Pc = (Wc+Ws)/2 - (1/2) Wc chi_c Wc
        Ps = -Ws + (3/2) Ws chi_s Ws
    For triplet:
        Pc = (Wc-Ws)/2 - (1/2) Wc chi_c Wc   (charge part same sign change)
        Ps = Ws - (1/2) Ws chi_s Ws            (spin part sign flip)

    Returns Pc_q and Ps_q separately for backward compatibility.

    Returns
    -------
    Pc_q : ndarray
        Charge vertex, shape (norb, norb, Nx, Ny, Nz).
    Ps_q : ndarray
        Spin vertex, shape (norb, norb, Nx, Ny, Nz).
    """
    # Symmetrise FIRST (PR #129 round 3): this path read the raw tables,
    # so a one-sided off-site declaration entered with its bare phase
    # (v e^{-iqR} since issue #232) while the
    # ring and the general S/C route read the same Hamiltonian as
    # v cos(qR) -- measured drift 0.7 at q = pi/2 for V(R=+x) = 0.7.
    # SKIPPED when the caller proved the raw declarations partner-closed
    # (round 5): the reduction is only ALGEBRAICALLY idempotent after the
    # exponential transform -- re-averaging a closed configuration mixed
    # independently rounded +R/-R exponentials and changed saved
    # artifacts at the 1e-18 level, violating byte parity for
    # previously-working symmetric runs.
    if not declarations_closed:
        inter_k = _symmetrise_interactions_k(
            {k: inter_k[k] for k in ("CoulombIntra", "CoulombInter")
             if k in inter_k})
    U_k = inter_k.get("CoulombIntra", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))
    V_k = inter_k.get("CoulombInter", np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex))

    # Transpose to batch dimension first: (Nx, Ny, Nz, norb, norb)
    Wc = (U_k + 2.0 * V_k).transpose(2, 3, 4, 0, 1).copy()
    Ws = (-U_k).transpose(2, 3, 4, 0, 1).copy()

    # chi0 at static limit: (Nx, Ny, Nz, norb, norb)
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    chi0_static = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()

    I_mat = np.broadcast_to(np.eye(norb, dtype=complex), (Nx, Ny, Nz, norb, norb)).copy()

    # Batched solve
    mat_s = I_mat + chi0_static @ Ws
    mat_c = I_mat + chi0_static @ Wc

    chis = np.linalg.solve(mat_s, chi0_static)
    chic = np.linalg.solve(mat_c, chi0_static)

    WsChisWs = Ws @ chis @ Ws
    WcChicWc = Wc @ chic @ Wc

    if pairing_type == "singlet":
        Pc_all = (Wc + Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = -Ws + 1.5 * WsChisWs
    elif pairing_type == "triplet":
        Pc_all = (Wc - Ws) / 2.0 - 0.5 * WcChicWc
        Ps_all = Ws - 0.5 * WsChisWs
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Transpose back: (Nx, Ny, Nz, norb, norb) -> (norb, norb, Nx, Ny, Nz)
    Pc_q = Pc_all.transpose(3, 4, 0, 1, 2)
    Ps_q = Ps_all.transpose(3, 4, 0, 1, 2)

    return Pc_q, Ps_q


def _compute_vertices_general(chi0q, inter_k, norb, Nx, Ny, Nz, nmat,
                              pairing_type="singlet", static_index=None):
    """Compute effective pairing interaction using general S,C matrices.

    For singlet (Kuroki et al., arXiv:0902.3691, Eq.(6)):
        V^s = (3/2) S chi_s S - (1/2) C chi_c C + (1/2)(S + C)

    For triplet (Takimoto et al., PRB 69, 104504):
        V^t = -(1/2) S chi_s S - (1/2) C chi_c C + (1/2)(C - S),
    the bare term without its cross elements (ab, ab), a != b
    (:func:`_triplet_bare_vertex`, issue #262).

    Supports both 2-index (reduced) and 4-index (general) chi0q:
    - 2-index: shape (norb, norb, Nx, Ny, Nz, nmat)
      Expanded to diagonal of (norb^2, norb^2) matrix.
    - 4-index: shape (norb, norb, norb, norb, Nx, Ny, Nz, nmat)
      Used directly as (norb^2, norb^2) matrix. This is the correct
      treatment for multi-orbital systems.

    Returns
    -------
    V_q : ndarray
        Effective pairing interaction, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb
    chi0q_is_4index = (chi0q.ndim == 8)

    # Build S, C matrices for all q-points at once: (Nx, Ny, Nz, nd, nd)
    S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)

    # Extract chi0 at static limit for all q-points
    # (zero bosonic frequency; nmat//2 only for a full frequency grid)
    si = nmat // 2 if static_index is None else static_index
    if chi0q_is_4index:
        # (norb, norb, norb, norb, Nx, Ny, Nz, nmat) -> (Nx, Ny, Nz, nd, nd)
        chi0_static = chi0q[:, :, :, :, :, :, :, si].reshape(
            nd, nd, Nx, Ny, Nz).transpose(2, 3, 4, 0, 1).copy()
    else:
        # (norb, norb, Nx, Ny, Nz, nmat) -> expand to (Nx, Ny, Nz, nd, nd)
        chi0_2d = chi0q[:, :, :, :, :, si].transpose(2, 3, 4, 0, 1).copy()
        # chi0_2d shape: (Nx, Ny, Nz, norb, norb)
        # The guard runs for EVERY 2-index input, including norb == 1: the
        # rejection of Exchange/PairHop applies there too (the one-orbital
        # builder gives them zero weight, so accepting them silently omits
        # the interaction -- the norb shortcut below bypassed the helper
        # entirely, round 8). Only the partial-dressing warning is
        # norb-gated, inside the helper.
        _warn_reduced_flex_missing_components(
            inter_k, norb, Nx, Ny, Nz, source="a reduced 2-index chi0q",
            sc_matrices=(S_all, C_all))
        if norb == 1:
            chi0_static = chi0_2d.reshape(Nx, Ny, Nz, 1, 1)
        else:
            # A 2-index (reduced) chi0q is the density-density diagonal
            # of the susceptibility: chi0_2d[a, b] IS chi0_{(a,a),(b,b)} (the
            # matching interaction reduction is einsum('kaabb->kab', ...) in
            # RPA._inflate_chi0q_and_ham).  With the orbital-pair flat index
            # (l1,l2) -> l1*norb + l2 used by _build_sc_matrices_all_q, it
            # therefore belongs at the density-pair positions:
            #
            #     chi0_{(a,a),(b,b)} = chi0_2d[a, b],  everything else zero.
            #
            # The historical placement chi0_{(l1,l2),(l3,l2)} = chi0_2d[l1,l3]
            # (a delta_{l2,l4} scatter, i.e. kron(chi0_2d, I_norb)) read the
            # density-pair index as an ordinary orbital index: it dropped the
            # inter-orbital density coupling chi0_{(0,0),(1,1)} and scattered
            # chi0_2d onto pair indices the reduced scheme never computed.  This
            # is the same defect that _expand_flex_chi carried on the
            # chi0q_mode="flex" route.
            #
            # Off-density rows/columns stay exactly zero.  For the interaction
            # terms that put S/C weight there (CoulombInter, Hund, Ising,
            # Exchange, PairHop) those channels are then undressed -- an honest
            # reflection of what a reduced chi0q contains, rather than a
            # fabricated dressing.
            chi0_static = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)
            dens = np.arange(norb) * norb + np.arange(norb)
            chi0_static[..., dens[:, None], dens[None, :]] = chi0_2d

    # Batched RPA solve for all q-points simultaneously
    # chi_s = [I - chi0 @ S]^{-1} @ chi0
    # chi_c = [I + chi0 @ C]^{-1} @ chi0
    # The solve itself lives in bond_channels.dress_bond, the SINGLE shared
    # dressing helper: the bond path runs the identical algebra at the enlarged
    # bond-major size ND = nd*B, so factoring it out lets the bond reduction
    # test exercise the real production code path instead of a mirror copy.
    # dress_bond performs exactly the same operations in the same order, so
    # this call is bit-identical to the historical inline body.
    #
    # cond_tol=None DISABLES dress_bond's near-singularity guard here: that
    # guard is a BOND-PATH-ONLY policy. This is the legacy default path, which
    # has always simply called np.linalg.solve; an invertible but poorly
    # conditioned RPA denominator must keep returning its (large but finite)
    # result rather than becoming a ValueError for existing users.
    chis, chic = bond_channels.dress_bond(chi0_static, S_all, C_all,
                                          cond_tol=None)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + _triplet_bare_vertex(S_all, C_all, norb)
    else:
        raise ValueError("Unknown pairing_type: '{}'. Use 'singlet' or 'triplet'.".format(
            pairing_type))

    # Reshape (Nx, Ny, Nz, nd, nd) -> (norb, norb, norb, norb, Nx, Ny, Nz)
    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(3, 4, 5, 6, 0, 1, 2)

    return Vs_q


def _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz, sc_matrices=None):
    """Largest |S| or |C| on the blocks a reduced chi cannot dress.

    The reduced susceptibility is zero on every pair index (a,b) with a != b, so
    the vertex is undressed exactly where S or C is nonzero there --
    S/C[(a,b),(a,b)] (case 2) and S/C[(a,b),(b,a)] (case 4) of
    :func:`_build_sc_matrices_all_q`.

    This inspects the COMBINED matrices rather than testing each configured term
    on its own. Those blocks are sums: under the adjudicated slot table
    (#113) case 2 mixes CoulombInter, Ising, Hund AND Exchange, and case 4
    carries PairHop alone. Terms can cancel there (equal CoulombInter and
    Hund cancel case 2 exactly in both channels), and a per-term test would
    then announce missing dressing that does not exist.
    """
    # Reuse the caller's matrices when it already has them: at the full grid
    # these are O(Nq * norb^4) each, and building a second pair while the first
    # is live doubles the allocation for what is only a diagnostic.
    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    nd = norb * norb
    density = np.zeros(nd, dtype=bool)
    density[np.arange(norb) * norb + np.arange(norb)] = True
    off = ~density
    if not off.any():
        return 0.0
    weight = 0.0
    for M in (S_all, C_all):
        # anything outside the density x density sub-block, over all q
        weight = max(weight,
                     float(np.max(np.abs(M[..., off, :]))),
                     float(np.max(np.abs(M[..., :, off]))))
    return weight


def _build_vertex_sc_matrices(convention, inter_k, norb, Nx, Ny, Nz):
    """Build the (S, C) interaction matrices for the given orbital convention.

    Single dispatcher shared by :func:`_compute_vertices_flex` and its
    callers, so a dynamic run can build the frequency-INDEPENDENT matrices
    once and pass them to every per-frequency contraction instead of
    rebuilding ~nmat identical full-grid pairs.
    """
    conv = str(convention).lower()
    if conv == "myo":
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        return build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz)
    if conv == "kuroki":
        return _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
    raise ValueError(
        "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(convention))


def _reject_reduced_flex_unsupported(inter_k, convention="kuroki",
                                     source=None):
    """Reject Exchange/PairHop wherever a Kuroki-convention vertex is built.

    This is the cheap half of :func:`_warn_reduced_flex_missing_components`:
    a dict scan with no S/C construction, so it can (and must) run at EVERY
    Kuroki vertex boundary -- :func:`_compute_vertices_flex` itself, hence
    also ``eliashberg_dynamic._instantaneous_vertex`` and every
    per-frequency dynamic contraction -- making the rejection an enforced
    invariant rather than an ordering artifact of which builder the
    orchestrator happens to call first (round-9 review: the IR
    instantaneous route silently accepted both terms and returned a zero
    vertex when called directly).

    Exchange and PairHop carry NO density-diagonal vertex content
    (hwave.solver.vertex_table), so a Kuroki (reduced-origin) chi dresses
    none of it and the result would silently omit the interaction; the
    general (myo) convention represents them fully and is not restricted.
    """
    if str(convention).lower() != "kuroki":
        return
    rejected = [k for k in _REDUCED_FLEX_REJECTED
                if k in inter_k and np.any(np.asarray(inter_k[k]) != 0)]
    if rejected:
        raise ValueError(
            "the Eliashberg vertex cannot be built from {} together with "
            "{}: those interactions have no density-diagonal vertex "
            "content at all (hwave.solver.vertex_table), so reduced "
            "(density-only) data dresses none of it and the result would "
            "silently omit the interaction. (H-wave's own reduced "
            "runs cannot even be produced with these terms since the "
            "unified scheme policy.) Provide a general (four-index) "
            "susceptibility instead -- for a FLEX source, re-run with "
            "calc_scheme='general'.".format(
                source or ("a REDUCED (calc_scheme='reduced') FLEX "
                           "susceptibility"),
                ", ".join(rejected)))


def _warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz,
                                          convention="kuroki",
                                          source=None, sc_matrices=None):
    """Warn when a reduced (kuroki) FLEX chi cannot support the interaction.

    Call this ONCE per run, from the place that is about to build the pairing
    vertex -- NOT from inside ``_compute_vertices_flex``. That function is
    invoked once per bosonic Matsubara frequency by the dynamic kernel (so the
    warning would repeat ``Nmat`` times, ~1000 in production runs), and again by
    ``eliashberg_dynamic._instantaneous_vertex`` with chi = 0, where the
    message would be doubly misleading. The Exchange/PairHop REJECTION, by
    contrast, is enforced inside ``_compute_vertices_flex`` on every call
    (:func:`_reject_reduced_flex_unsupported`) -- it is a cheap scan and
    must hold on every route, not just the orchestrated one.

    A ``calc_scheme="reduced"`` FLEX run stores only the
    density-density diagonal chi_{(a,a),(b,b)} of the susceptibility.  The
    Kuroki S/C matrices built from CoulombIntra alone live entirely on that
    density-pair block, so the reduced treatment is exact.  CoulombInter,
    Hund and Ising also populate the off-density block S/C[(a,b),(a,b)]
    with a != b -- and there chi is identically zero simply because the
    reduced run never computed it.  Those channels then keep only the bare
    0.5*(S+C) term with no fluctuation dressing: a silent approximation
    rather than a solver error, hence the WARNING.  Exchange and PairHop
    carry NO density-diagonal vertex content at all (vertex_table), so a
    reduced chi dresses nothing of them; combined with the #120 scheme
    policy (no reduced run can be produced with them) that
    combination is REJECTED rather than warned about.

    This is a genuine limitation of the stored data, not of the loader: it
    cannot be repaired on the Eliashberg side.  The universal remedy is a
    general (four-index) susceptibility; for a FLEX source that means
    re-running with ``calc_scheme="general"`` (which stores the full
    orbital-pair chi).
    """
    if str(convention).lower() != "kuroki":
        # Only the reduced route stores a density-only chi; the
        # general (myo) path carries the full orbital-pair susceptibility.
        return
    # The rejection must run BEFORE the norb == 1 shortcut: a one-orbital
    # Exchange/PairHop has zero S/C weight in this builder too, so the
    # interaction would be silently omitted rather than represented --
    # exactly what the rejection exists to prevent (round-7 review; the
    # bypass returned zero vertices for both terms on both routes).
    _reject_reduced_flex_unsupported(inter_k, convention, source)
    if norb <= 1:
        # norb == 1 has no off-density pair index, so no PARTIAL dressing
        # can be missing (the rejection above still applies).
        return
    configured = [k for k in _REDUCED_FLEX_PARTIAL if k in inter_k]
    if not configured:
        return
    # Ask the assembled S/C matrices, not the individual terms: the off-density
    # blocks are sums and the configured terms can cancel there (e.g. equal
    # CoulombInter and Hund cancel case 2 in both channels under the
    # adjudicated table). Only a nonzero combined block means dressing is
    # actually missing.
    if _off_density_sc_weight(inter_k, norb, Nx, Ny, Nz,
                              sc_matrices) == 0.0:
        return
    # Decision above is on the assembled matrices; attribution here is per term,
    # so an inert term that happens to be configured is not named as a cause.
    # Attribute through the same assembled-and-symmetrised path as the
    # decision: a raw declaration that symmetrises to zero (#112) must not
    # be named as a cause.
    def _one_type(k):
        split = getattr(inter_k, "locality_split", None)
        if split is None:
            return {k: inter_k[k]}
        on, off = split
        return SiteSplitInteractionK(
            {k: inter_k[k]},
            ({k: on[k]} if k in on else {}, {k: off[k]} if k in off else {}))

    missing = [k for k in configured
               if _off_density_sc_weight(_one_type(k), norb,
                                         Nx, Ny, Nz) != 0.0] or configured
    logger.warning(
        "The Eliashberg vertex is being built from %s, which carries only the "
        "density-density components chi_{(a,a),(b,b)}, together with "
        "inter-orbital interaction(s) %s. Together those terms give the S/C "
        "matrices nonzero "
        "off-density components S/C[(a,b),(a,b)] and S/C[(a,b),(b,a)] (a != b) "
        "for which the reduced run computed no susceptibility, so those "
        "channels enter the pairing vertex undressed (bare 0.5*(S+C) only) and "
        "the resulting lambda is an approximation. calc_scheme='general' stores "
        "the full orbital-pair chi and removes this gap; note its off-site "
        "support is limited to same-orbital CoulombInter without sublattice "
        "folding -- for a model with off-site INTER-orbital interactions "
        "there is currently no FLEX-dressed vertex without this "
        "approximation.",
        source or ("a REDUCED (calc_scheme='reduced') FLEX "
                   "susceptibility"),
        ", ".join(missing))


def _triplet_bare_vertex(S_all, C_all, norb):
    """The bare (first-order) TRIPLET term of the pair-space vertex,
    ``1/2 (C - S)`` without its cross slots ``(ab, ab)``, ``a != b``.

    An inter-orbital density-type term sits in the spin/charge matrices twice:
    its Hartree content in the density slot ``(aa, bb)`` and its Fock content
    in the cross slot ``(ab, ab)`` (the particle-hole fluctuations need both).
    In the pairing channel the cross slot is the exchange image of the density
    slot, and on a gap of triplet parity the two give the SAME contribution:
    ``1/2 (C - S)`` over both slots counts an on-site ``U'``, Hund or Ising term
    twice (issue #262); an off-site row collapsed onto the q-only vertex (the
    non-bond path before issue #239, and a plain table without the locality
    split) puts its Fock content there too. The density slot alone is exact. In the singlet channel
    the cross slot of these types carries ``S + C = 0`` and the term is right
    as it stands; the types whose pairing lives in the cross slot (Exchange)
    have ``C - S = 0`` there, so nothing is lost by dropping it here. Validated
    against an independent first-order spin-orbital construction for every
    type, on-site and off-site (tests/test_triplet_bare_vertex.py)."""
    bare = 0.5 * (C_all - S_all)
    cross = [a * norb + b for a in range(norb) for b in range(norb) if a != b]
    if cross:
        bare[..., cross, cross] = 0.0
    return bare


def _compute_vertices_flex(chis, chic, inter_k, norb, Nx, Ny, Nz,
                           pairing_type="singlet", convention="kuroki",
                           sc_matrices=None):
    """Compute pairing vertex from pre-computed FLEX susceptibilities.

    Uses the same formula as _compute_vertices_general but takes
    chi_s and chi_c directly instead of computing them from chi0q via RPA.
    The triplet bare term is ``1/2 (C - S)`` without its cross elements
    ``(ab, ab)``, ``a != b`` (:func:`_triplet_bare_vertex`, issue #262).
    This allows using dressed susceptibilities from the FLEX solver.

    ``convention`` records which FLEX path produced the susceptibilities:
    "kuroki" (default; reduced path) or "myo" (general full-vertex path). The
    two S/C builders historically differed in the charge ``C(ab,ab)`` element
    (``-U'+2J`` vs ``-U'+J``); the exact-diagonalization adjudication of the
    per-type vertex content (issue #113: Hund ``+J`` and Exchange ``+J'``
    there) made them IDENTICAL, so the flag no longer selects different
    matrices -- it remains a provenance / layout discriminator, and legacy
    files are guarded by ``sc_vertex_version`` instead. ``chis``/``chic`` are
    expected in the native [a,c,b,d] orbital-pair order regardless
    (issue #78).

    Parameters
    ----------
    chis : ndarray
        Spin susceptibility at static limit, shape (Nx, Ny, Nz, nd, nd)
        where nd = norb^2.
    chic : ndarray
        Charge susceptibility at static limit, same shape.
    inter_k : dict
        Interactions in k-space from _build_interaction_k.
    norb : int
        Number of orbitals.
    Nx, Ny, Nz : int
        Grid dimensions.
    pairing_type : str
        "singlet" or "triplet".
    sc_matrices : tuple of ndarray, optional
        Precomputed ``(S_all, C_all)`` for this convention
        (:func:`_build_vertex_sc_matrices`). The matrices are
        frequency-independent, so a dynamic run passes one pair to every
        per-frequency call instead of rebuilding ~nmat identical full-grid
        pairs (round-9 review).

    Returns
    -------
    Vs_q : ndarray
        Pairing vertex, shape (norb, norb, norb, norb, Nx, Ny, Nz).
    """
    nd = norb * norb

    if "PairLift" in inter_k:
        logger.warning(
            "PairLift is configured but does not contribute to the S/C pairing "
            "vertex (S=C=0); it is ignored in the Eliashberg calculation.")

    # Validate the convention BEFORE the cache branch: with precomputed
    # sc_matrices the dispatcher below is skipped, and an unknown tag would
    # otherwise be silently treated as "not kuroki" = unrestricted
    # (round-10 review reproduced convention="invalid" returning a vertex).
    conv = str(convention).lower()
    if conv not in ("kuroki", "myo"):
        raise ValueError(
            "Unknown convention: '{}'. Use 'kuroki' or 'myo'.".format(
                convention))

    # Enforced at THIS boundary (not only in the orchestrator): every Kuroki
    # vertex construction -- including chi = 0 via _instantaneous_vertex --
    # must reject Exchange/PairHop, or the interaction is silently omitted.
    _reject_reduced_flex_unsupported(inter_k, conv)

    if sc_matrices is not None:
        S_all, C_all = sc_matrices
    else:
        S_all, C_all = _build_vertex_sc_matrices(conv, inter_k,
                                                 norb, Nx, Ny, Nz)

    SChisS = S_all @ chis @ S_all
    CChicC = C_all @ chic @ C_all

    if pairing_type == "singlet":
        Vs_all = 1.5 * SChisS - 0.5 * CChicC + 0.5 * (S_all + C_all)
    elif pairing_type == "triplet":
        Vs_all = -0.5 * SChisS - 0.5 * CChicC + _triplet_bare_vertex(S_all, C_all, norb)
    else:
        raise ValueError("Unknown pairing_type: '{}'".format(pairing_type))

    Vs_q = Vs_all.reshape(Nx, Ny, Nz, norb, norb, norb, norb).transpose(
        3, 4, 5, 6, 0, 1, 2)

    return Vs_q


def _offsite_mirror_asymmetry(table, grid=None):
    """Max |X[R,(a,b)] - X[-R,(a,b)]| over the off-site entries of a raw
    declaration table: zero exactly when moving every off-site row to the
    opposite bond leaves the table unchanged, i.e. when the bond placement
    of issue #232 cannot change any result computed from it.

    With ``grid`` (the ``CellShape`` the run used) the comparison is made
    on that grid -- keys are accumulated modulo the grid, so wrapped or
    aliased displacements, and a displacement that is its own reverse on
    the grid (``+x`` and ``-x`` on a 2-site axis), compare as the vertex
    sees them. Without it the keys are compared literally, a missing mirror
    reading as zero, which fails CLOSED (reported as asymmetric)."""
    if grid is not None:
        grid = tuple(int(n) for n in grid)
        agg = {}
        for (irvec, orbvec), v in table.items():
            r = tuple(int(x) % n for x, n in zip(irvec, grid))
            key = (r, tuple(orbvec))
            agg[key] = agg.get(key, 0.0) + complex(v)
        table = agg
        mirror_of = lambda r: tuple((-x) % n for x, n in zip(r, grid))
        onsite = tuple(0 for _ in grid)
    else:
        mirror_of = lambda r: tuple(-x for x in r)
        onsite = (0, 0, 0)
    worst = 0.0
    for (irvec, orbvec), v in table.items():
        irvec = tuple(int(x) for x in irvec)
        if irvec == onsite:
            continue
        mirror = table.get((mirror_of(irvec), tuple(orbvec)), 0.0)
        diff = abs(complex(v) - complex(mirror))
        if not np.isfinite(diff):
            return float("inf")
        worst = max(worst, diff)
    return worst


def _decode_version_stamp(data, key, path):
    """``int`` value of the scalar stamp ``key`` in an npz, ``None`` when
    absent; anything but one finite, exactly integral number is refused.
    The field is read once, so a malformed stamp reports as such."""
    if key not in getattr(data, "files", data):
        return None
    shown = "<unreadable>"
    try:
        raw = data[key]
        shown = repr(raw)
        arr = np.asarray(raw)
        if arr.dtype.kind not in "iufc" or arr.size != 1:
            raise ValueError("not a numeric scalar")
        val = complex(arr.reshape(())[()])
        if val.imag != 0.0 or not (np.isfinite(val.real)
                                   and float(val.real).is_integer()):
            raise ValueError("not an integer")
    except (TypeError, ValueError, OverflowError):
        raise ValueError("'{}' carries a malformed {} field ({}).".format(path, key, shown))
    return int(val.real)


def placement_sensitive_types(interactions, grid=None):
    """The declared two-body types whose off-site table is not symmetric
    under R -> -R (on ``grid`` when given) -- the declarations for which an
    artifact written before issue #232 describes a different Hamiltonian."""
    if not interactions:
        return []
    return [t for t in _PLACEMENT_TYPES
            if interactions.get(t)
            and _offsite_mirror_asymmetry(interactions[t], grid=grid) > 0.0]


def check_placement_stamp(data, path, interactions, grid=None, what="FLEX file"):
    """Refuse a FLEX artifact written before the bond-placement correction
    of issue #232 (``sc_vertex_version`` absent or below
    :data:`PLACEMENT_VERSION`) when the declaration is placement-sensitive;
    accept it otherwise."""
    ver = _decode_version_stamp(data, "sc_vertex_version", path)
    if ver is not None and ver >= PLACEMENT_VERSION:
        return
    sensitive = placement_sensitive_types(interactions, grid)
    if sensitive:
        raise ValueError(
            "{} '{}' was written before the bond-placement correction of issue #232 "
            "(sc_vertex_version {}), which put an off-site row on the bond opposite its "
            "Transfer row, and {} declare(s) off-site rows that are not symmetric under "
            "R -> -R, so the file describes a different Hamiltonian; regenerate it with "
            "the current code.".format(what, path, "absent" if ver is None else ver,
                                       ", ".join(sensitive)))
