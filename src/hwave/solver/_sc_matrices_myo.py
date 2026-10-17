"""Spin/charge interaction matrices for the full-vertex FLEX path.

Historically this module carried a hand-maintained sibling of
``hwave.sc._build_sc_matrices_all_q`` differing in exactly one element: the
charge (ab, ab) slot, where it used the MYO (cond-mat/0407094) value ``-U'+2J``
against the Kuroki (arXiv:0902.3691) value ``-U'+J`` used by ``sc.py``.

Exact diagonalization of H-wave's documented file operators adjudicated that
difference (issue #113): per interaction type the exact charge vertex at that
slot is ``+J`` from the Hund file term and ``+J`` from the Exchange file term.
For the SU(2) Kanamori combination (Hund + Exchange at equal J) their sum
reproduces the standard literature value ``-U'+2J`` -- so the MYO matrix was
right for the COMBINATION but wrong as a per-type attribution, assigning the
whole ``2J`` to the Hund file entry and nothing to Exchange. With the
per-type values fixed in ``sc.py``, the two builders are identical, and
:func:`build_sc_matrices_myo` simply re-exports the single implementation.

:func:`build_sc_matrices_locality_split` is the FLEX general path's
locality-aware entry (#181, Tier 1). The shared builder applies the on-site
Kanamori slot map to whatever it is handed, with no locality bookkeeping,
which is wrong for off-site input in two ways: an off-site entry lands in
the cross ``(ab,ab)`` / antidiag ``(ab,ba)`` families, whose particle-hole
pair is NON-local for ``R != 0`` (legs on orbital a at 0 and b at R) and
therefore not a function of q alone; and the density-family gate
``l1 != l3`` -- right on-site, where an orbital has no Hund/Ising coupling
with itself -- deletes the physical same-orbital off-site Hund/Ising. The
split entry hands the shared builder the WHOLE table (authoritative for
every density element, in the reader's own entry order) plus two
locality-ROUTING inputs: the on-site part feeds the cross/antidiag
families, the off-site part supplies only the same-orbital ``(aa,aa)``
density elements of Hund/Ising that the on-site gate would delete. Off-site
content therefore reaches the density family only -- the Hartree vertex
``V_ab(q)`` on ``(aa,bb)`` including ``a == b`` -- which is q-representable
exactly and is what the RPA ring carries for off-site input (measured
element-complete equal, tests/test_flex_offsite_general.py). The output is
NOT the sum of two independently built parts. The exchange (Fock) crossing
of an off-site term is deliberately ABSENT: it needs the bond-resolved
vertex (#181, Tier 3). ``hwave.sc``'s own callers (Eliashberg) keep the
shared builder unchanged.
"""

from hwave.sc import _build_sc_matrices_all_q

# The interaction types whose off-site content has a q-only representation
# in the pair-space S/C matrices: density-density on each site, so the
# whole displacement dependence sits in the coefficient V_ab(q) of the
# (aa,bb) density slots. Everything else in the off-site part is a caller
# error (see build_sc_matrices_locality_split). Owned by
# hwave.solver.offsite (one definition, #181 Tier 3); re-exported here for
# the builder's own validation and for existing importers.
from hwave.solver.offsite import _OFFSITE_DENSITY_TYPES  # noqa: E402,F401


def build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz):
    """One implementation: see :func:`hwave.sc._build_sc_matrices_all_q`."""
    return _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)


def build_sc_matrices_locality_split(inter_k, inter_k_onsite, inter_k_offsite,
                                     norb, Nx, Ny, Nz):
    """S/C matrices from a PRE-fold-locality split of the interaction.

    Parameters
    ----------
    inter_k : dict
        k-space interactions (``hwave.sc._build_interaction_k``) built from
        the WHOLE table, in the caller's own entry order. Every slot
        element that exists without a split is read from this, so their
        floating-point summation order is unchanged.
    inter_k_onsite : dict
        Built from the ``R == 0`` declarations only. Feeds the cross
        ``(ab,ab)`` and antidiag ``(ab,ba)`` families, whose particle-hole
        pair is non-local for ``R != 0``.
    inter_k_offsite : dict
        Built from the ``R != 0`` declarations only. Feeds the same-orbital
        density elements ``(aa,aa)`` of Hund / Ising, which the on-site
        gate of the shared builder otherwise deletes. Must contain only
        CoulombInter / Hund / Ising: ``ValueError`` otherwise (fail closed
        -- no off-site type is silently dropped here).
    norb, Nx, Ny, Nz : int

    Returns
    -------
    S_all, C_all : ndarray, shape (Nx, Ny, Nz, norb**2, norb**2)

    Bit-identical to the shared builder for on-site input and for the
    class the general path accepted before the split (same-orbital
    off-site CoulombInter), see tests/test_sc_matrices_locality_split.py.
    """
    bad = [t for t in (inter_k_offsite or {}) if t not in _OFFSITE_DENSITY_TYPES]
    if bad:
        raise ValueError(
            "off-site interaction part contains {}: only {} have a "
            "q-representable (density-slot) off-site vertex in the "
            "pair-space S/C matrices; the caller must reject or route "
            "other off-site types before building.".format(
                ", ".join(sorted(bad)), "/".join(_OFFSITE_DENSITY_TYPES)))
    return _build_sc_matrices_all_q(
        inter_k, norb, Nx, Ny, Nz,
        locality_split=(inter_k_onsite, inter_k_offsite or {}))
