#!/usr/bin/env python3

"""Contract tests and the H1 calibration gate for
``ed_oracle_util.SectorED.bond_correlator_transverse`` (Phase H of the
bond-transverse campaign, issue #110's validation groundwork).

See ``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``
("Phase H -- the transverse ED harness") for the binding contract this
module pins, and ``docs/superpowers/plans/2026-08-15-bond-transverse-
phase-h.md`` for the task-by-task build order. This module (Task 1)
covers:

- the operator-construction contract (shape/dtype, R-wrapping duplicate
  rejection, the mixed-spin sector delta) via
  ``TestBondCorrelatorTransverseContract``;
- H1: the R == 0 channels of ``bond_correlator_transverse`` reproduce the
  bespoke ``_chi_pm_ed`` dense-ED reference from
  ``tests/test_rpa_vs_ed_oracle.py::TestTransverseComplexPairHop`` on the
  SAME 2-site spinful fixture, at round-off (both sides are exact ED on
  the same Hilbert space -- ``assert_approx_array(rel=0, abs=1e-12)``).

- H2: the ED-pinned on-site 7-type transverse table (spec: "The verified
  on-site algebra"), RE-adjudicated through ``bond_correlator_transverse``
  itself rather than through the bespoke ``_chi_pm_ed`` reference H1
  checks against, via ``TestH2OnsiteTableReadjudication``.

H3 (the frame-map derivation and the V=0 bond pin) is a later task in
this same campaign and is not covered here.

Tests must be run from the repository root.
"""

import functools
import unittest

import numpy as np

import tests.test_rpa_vs_ed_oracle as oracle
from tests import ed_oracle_util
from tests.approx_util import assert_approx_array

# Richardson finite-difference step, the #151/campaign-wide standard
# (tests/test_bond_vs_ed_oracle.py's CAMPAIGN_V1) -- reused verbatim so
# H2's granules sit on the same finite-difference footing as every other
# ED-oracle granule in this campaign.
CAMPAIGN_V1 = 0.005

# NOTE: imported as a module (``oracle``), never as bare names -- binding
# ``oracle.TestTransverseComplexPairHop`` at module scope here would make
# unittest's discovery pick it up as a SECOND, duplicate TestCase under
# this module too.


class TestBondCorrelatorTransverseContract(unittest.TestCase):
    """Operator-construction contract, independent of any physics
    comparison: shape/dtype, R-wrapping duplicate rejection, and the
    mixed-spin sector displacement the new method relies on."""

    @staticmethod
    def _small_fx():
        # 2-site, two-orbital, spinful -- fast, disposable fixture for the
        # contract tests (distinct from the H1 fixture below, which must
        # be the EXACT 2-site/2-orbital TestTransverseComplexPairHop
        # fixture). norb=2 so channels use the (R, a, b) 3-tuple form
        # ``_channel_orbitals`` expects (the norb == 1 bare-``(R,)``
        # shorthand is not what the brief's example channels use).
        t = {(0, 0): 0.3 + 0.1j, (1, 1): 0.2 - 0.05j,
             (0, 1): 0.1j, (1, 0): 0.1j}
        return ed_oracle_util.EDFixture(
            L=2, norb=2, t=t, eps=(0.05, -0.02), T=0.5, mu=0.1)

    def test_output_shape_and_dtype(self):
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        channels = [(0, 0, 0), (1, 0, 0)]
        out = se.bond_correlator_transverse(channels)
        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        self.assertEqual(out.dtype, np.complex128)

    def test_duplicate_channels_after_wrapping_rejected(self):
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        # (1, 0, 0) and (1 - L, 0, 0) both wrap (mod L=2) to R=1 -- same
        # channel after wrapping, must be rejected rather than silently
        # double-counted.
        channels = [(1, 0, 0), (1 - fx.L, 0, 0)]
        with self.assertRaises(ValueError):
            se.bond_correlator_transverse(channels)

    def test_mixed_spin_delta_is_plus_minus_one(self):
        # Probe the sector displacement via the same operator-construction
        # API bond_correlator_transverse uses internally: a bilinear whose
        # creation leg is spin UP (mode parity 0) and annihilation leg is
        # spin DOWN (mode parity 1) must impart delta = (+1, -1) -- the
        # up-sector gains one particle, the down-sector loses one.
        fx = self._small_fx()
        se = ed_oracle_util.SectorED(fx)
        mode_terms = [(fx.mode(j, 0, 0), fx.mode(j, 0, 1), 1.0)
                      for j in range(fx.L)]
        delta, _op = se._build_operator(mode_terms)
        self.assertEqual(delta, (1, -1))


class TestH1AgainstBespokeChiPmEd(unittest.TestCase):
    """H1: on R == 0 channels, ``bond_correlator_transverse`` must
    reproduce ``TestTransverseComplexPairHop._chi_pm_ed`` exactly (both
    are exact diagonalization of the SAME Hamiltonian on the SAME 2-site,
    2-orbital spinful Hilbert space -- round-off agreement, not a
    physics-tolerance one)."""

    def test_r0_matches_bespoke_chi_pm_ed(self):
        fx, _C, _CD, h1 = oracle._fx2_state()
        se = ed_oracle_util.SectorED(fx, h1=h1)

        # channels ordered (a outer, b inner) so the flat (I, J) axes
        # reshape directly onto the bespoke reference's (a, c, b, d) axes
        # below: channel (0, a, b) <-> creation orbital a (up),
        # annihilation orbital b (down).
        channels = [(0, a, b) for a in range(oracle.NORB)
                    for b in range(oracle.NORB)]
        got = se.bond_correlator_transverse(channels)
        got = got.reshape(oracle.LX, oracle.NORB, oracle.NORB,
                           oracle.NORB, oracle.NORB)

        bespoke = oracle.TestTransverseComplexPairHop()
        expected = bespoke._chi_pm_ed()  # base Hamiltonian, no interaction

        assert_approx_array(got, expected, rel=0, abs=1e-12)


# ---------------------------------------------------------------------------
# H2: the on-site 7-type re-adjudication, THROUGH the harness
# (bond_correlator_transverse) rather than through the bespoke H1
# reference. Reuses the exact H1 fixture (oracle._fx2_state(), L=2,
# norb=2) and, for 5 of the 7 types, oracle._terms_for's already-ED-
# validated canonical quartic terms (module docstring: "Coverage: the
# spin-free general scheme for CoulombIntra, CoulombInter, Exchange and
# PairHop"; PairLift is exercised as a structural zero in
# TestFierzOffSitePin's spirit and needs no separate validation to serve
# as a PASS-ZERO control here). Ising and Hund are NOT covered by
# oracle._terms_for (its module docstring: "Hund and Ising ... are
# adjudicated elsewhere") and get local term builders below, built
# directly from the file-format-documented Hamiltonian
# (src/hwave/solver/uhfk.py's Hund/Ising comments) rather than reused
# from anywhere else, per the append-only discipline (ed_oracle_util.py
# is not touched this task).
# ---------------------------------------------------------------------------

def _terms_ising_hund(fx, kind, v):
    """Ising/Hund canonical on-site two-orbital (0, 1) quartic terms, at
    coupling ``v``, in the (p, q, r, s, coeff) form ``SectorED``/
    ``hf_h1_from_terms`` consume -- the two types ``oracle._terms_for``
    does not build (see module docstring above). Built directly from the
    file-format-documented Hamiltonian (src/hwave/solver/uhfk.py):

    - Ising:  H = v * (n_{0,up} - n_{0,dn}) * (n_{1,up} - n_{1,dn})
      ("the documented Hamiltonian is J (n_up - n_down)(n_up - n_down)",
      uhfk.py's Ising append comment) -- pure density-density, all four
      spin combinations with sign (+,-,-,+) for (up,up)/(up,dn)/(dn,up)/
      (dn,dn).
    - Hund:   H = -v * (n_{0,up} n_{1,up} + n_{0,dn} n_{1,dn})
      ("interaction coeffs : -J^Hund by convention", uhfk.py's Hund
      append comment) -- same-spin-only density-density, which is
      exactly why Hund's transverse vertex is structurally zero (spec:
      "a same-spin interaction cannot connect the up and down
      propagators of the transverse loop"): a single Wick contraction of
      a same-spin quartic term can never produce a spin-flip (up-
      creation/down-annihilation) bilinear.
    """
    terms = []
    for j in range(fx.L):
        if kind == "Ising":
            for su, su_sign in ((0, 1.0), (1, -1.0)):
                for sv, sv_sign in ((0, 1.0), (1, -1.0)):
                    p = fx.mode(j, 0, su)
                    r = fx.mode(j, 1, sv)
                    terms.append((p, p, r, r, v * su_sign * sv_sign))
        elif kind == "Hund":
            for s in range(2):
                p = fx.mode(j, 0, s)
                r = fx.mode(j, 1, s)
                terms.append((p, p, r, r, -v))
        else:
            raise ValueError(kind)
    return terms


_ORACLE_TERM_TYPES = ("CoulombIntra", "CoulombInter", "Exchange",
                       "PairLift", "PairHop")


def _terms_at(fx, kind, v, phase=1.0):
    """Dispatch to ``oracle._terms_for`` (5 types, already ED-validated
    longitudinally) or the local Ising/Hund builder above."""
    if kind in _ORACLE_TERM_TYPES:
        return oracle._terms_for(fx, kind, v, phase)
    return _terms_ising_hund(fx, kind, v)


@functools.lru_cache(maxsize=None)
def _h2_channels():
    norb = oracle.NORB
    return tuple((0, a, b) for a in range(norb) for b in range(norb))


def _frame_map(x, norb):
    """The Kubo pair-slot-swap transpose, applied exactly as
    ``tests.test_rpa_vs_ed_oracle.TestTransverseComplexPairHop
    ._ed_first_order_pm`` applies it: reshape the flat (q, ND, ND)
    channel-space tensor onto its (q, a, c, b, d) orbital axes (ND =
    norb**2, channel I <-> (a, c), channel J <-> (b, d), matching how
    ``_h2_channels`` orders ``(0, a, b)``), transpose (0, 3, 4, 1, 2),
    and flatten back. This converts a RAW Lehmann-frame tensor (as
    ``bond_correlator_transverse`` returns it) into the frame
    ``_assemble_transverse_vertex``'s own ``ham_pm`` lives in."""
    shape = x.shape
    x5 = x.reshape(shape[0], norb, norb, norb, norb)
    return np.transpose(x5, (0, 3, 4, 1, 2)).reshape(shape)


@functools.lru_cache(maxsize=None)
def _h2_chi0_solver_frame():
    """The V=0 transverse correlator (exact ED, R=0 channels), frame-
    mapped into solver space -- the harness's own stand-in for the
    solver's bare transverse bubble at this on-site fixture."""
    fx, _C, _CD, _H1 = oracle._fx2_state()
    channels = _h2_channels()
    se = ed_oracle_util.SectorED(fx)
    raw = se.bond_correlator_transverse(list(channels))
    return _frame_map(raw, oracle.NORB)


def _h2_ed_derivative_solver_frame(kind, phase=1.0):
    """Richardson-extrapolated d(bond_correlator_transverse)/dv at v -> 0
    (the full-minus-HF-only pattern tests/test_bond_vs_ed_oracle.py's
    ``_unit_direction_raw`` uses for the longitudinal channel, mirrored
    here for the transverse one), frame-mapped into solver space.
    Returns (D_v1, D_vhalf), the two Richardson estimates
    ``adjudicate_granule`` needs to measure ED-side convergence.

    PairLift is special-cased onto the DENSE reference
    (``oracle.TestTransverseComplexPairHop._chi_pm_ed``, the exact
    machinery H1 pins ``bond_correlator_transverse`` against at v=0)
    rather than ``SectorED``: PairLift's declared operator
    ``c^+_{a,up} c_{a,dn} c^+_{b,up} c_{b,dn} + h.c.`` changes
    ``(N_up, N_dn)`` by ``(+2, -2)`` per term (a genuine two-magnon
    creation/annihilation operator -- Sz is not separately conserved,
    only total N), which ``SectorED`` structurally cannot represent (its
    per-(N_up, N_dn)-block Hamiltonian assumption; see the spec's Phase S
    section: "SectorED conserves (N_up, N_dn) and cannot represent
    Sz-breaking fixtures"). This is not a defect in
    ``bond_correlator_transverse`` -- ``SectorED.__init__`` correctly
    raises rather than silently truncating the off-block content -- so
    PairLift alone borrows the FULL (unrestricted) Hilbert-space Lehmann
    machinery, which has no such assumption. ``_chi_pm_ed``'s output
    layout ((q, a, c, b, d), axes matching H1's own reshape) is byte-for-
    byte the same convention ``bond_correlator_transverse``'s
    (q, I, J) output reshapes onto (H1 pins this at v=0 to round-off),
    so the two are drop-in interchangeable here."""
    fx, _C, _CD, _H1 = oracle._fx2_state()
    channels = list(_h2_channels())
    norb = oracle.NORB

    if kind == "PairLift":
        bespoke = oracle.TestTransverseComplexPairHop()

        @functools.lru_cache(maxsize=None)
        def X_of_v(v):
            terms = _terms_at(fx, kind, v, phase)
            hint = ed_oracle_util.h_int_from_terms(fx, terms)
            full = bespoke._chi_pm_ed(hint=hint)
            hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
            hf_only = bespoke._chi_pm_ed(h1=hf_h1)
            diff = (full - hf_only).reshape(fx.L, norb * norb, norb * norb)
            return _frame_map(diff, norb)
    else:
        @functools.lru_cache(maxsize=None)
        def X_of_v(v):
            terms = _terms_at(fx, kind, v, phase)
            full = ed_oracle_util.SectorED(
                fx, terms=terms).bond_correlator_transverse(channels)
            hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
            hf_only = ed_oracle_util.SectorED(
                fx, terms=(), h1=hf_h1).bond_correlator_transverse(channels)
            return _frame_map(full - hf_only, norb)

    D_v1 = ed_oracle_util.richardson(X_of_v, CAMPAIGN_V1)
    D_vhalf = ed_oracle_util.richardson(X_of_v, CAMPAIGN_V1 / 2)
    return D_v1, D_vhalf


def _h2_ham_pm_expected(kind, coeff, norb=2):
    """The ED-pinned on-site transverse table (spec: "The verified
    on-site algebra"; independently corroborated here by
    tests/test_rpa_ladder.py::...
    test_transverse_vertex_matches_exact_diagonalization's captured
    ``ham_pm``, which pins the SAME ``(a, c, b, d)`` axis convention this
    function reproduces -- row channel I <-> (a, c), column channel J
    <-> (b, d), matching ``_frame_map``'s reshape), at UNIT declared
    coupling ``coeff`` on the fixture's orbital pair (0, 1) (CoulombIntra
    uses orbital 0 only, matching ``oracle._terms_for``'s declaration).
    ``coeff`` is the SLOPE-frame value (e.g. 1.0, or a unit-magnitude
    complex phase for PairHop) that ``_h2_ed_derivative_solver_frame``'s
    Richardson derivative is directly comparable to -- not a finite
    coupling."""
    built = np.zeros((norb, norb, norb, norb), dtype=complex)
    if kind == "CoulombIntra":
        built[0, 0, 0, 0] = -coeff
    elif kind == "CoulombInter":
        built[0, 1, 0, 1] = -coeff
        built[1, 0, 1, 0] = -coeff
    elif kind == "Ising":
        built[0, 1, 0, 1] = coeff
        built[1, 0, 1, 0] = coeff
    elif kind == "Exchange":
        built[0, 0, 1, 1] = -coeff
        built[1, 1, 0, 0] = -coeff
    elif kind == "Hund":
        pass                                       # identically zero
    elif kind == "PairLift":
        pass                                       # identically zero
    elif kind == "PairHop":
        built[0, 1, 1, 0] = -coeff
        built[1, 0, 0, 1] = -np.conj(coeff)
    else:
        raise ValueError(kind)
    return built.reshape(norb * norb, norb * norb)


class TestH2OnsiteTableReadjudication(unittest.TestCase):
    """H2: the ED-pinned on-site 7-type transverse table, re-adjudicated
    THROUGH ``bond_correlator_transverse`` (Task 1) rather than through
    the bespoke ``_chi_pm_ed`` reference H1 already checked against --
    the harness's own calibration pin, per
    docs/superpowers/specs/2026-08-15-bond-transverse-design.md ("Phase
    H -- the transverse ED harness").

    For each of the 7 declaration types on the H1 fixture (2-site,
    2-orbital, spinful, complex hoppings -- ``oracle._fx2_state()``), R=0
    channels only: the Richardson-extrapolated first-order response of
    ``bond_correlator_transverse`` (full interaction minus the
    Hartree-Fock self-energy insertion -- the genuine vertex correction,
    the #151 pattern) is frame-mapped by the Kubo pair-slot-swap
    transpose (established on R=0 by
    ``TestTransverseComplexPairHop.test_pair_slot_map_is_the_one_kubo_
    symmetry_selects``) and compared, via ``adjudicate_granule``, against
    ``-chi0 @ ham_pm_expected @ chi0`` -- the SAME sandwich
    ``_assemble_transverse_vertex``'s ``ham_pm`` enters the solver's RPA
    equation through, with ``chi0`` the (also frame-mapped) exact V=0
    transverse correlator and ``ham_pm_expected`` the spec's on-site
    table. Hund and PairLift are PASS-ZERO controls (their transverse
    vertex is identically zero); the other five are expected PASS.

    INVESTIGATION PROTOCOL (binding on any FAIL here): a FAIL does not,
    by itself, indict either the table (ED-proven independently --
    src/hwave/solver/rpa.py's ``_assemble_transverse_vertex`` docstring:
    "verified against exact diagonalization end to end with one common
    scale across all seven types"; re-confirmed by
    tests/test_rpa_ladder.py's
    ``test_transverse_vertex_matches_exact_diagonalization``) or this
    harness (H1-calibrated against the bespoke reference to round-off).
    Before assigning fault, check independently:

    1. Operator phases -- is each channel's creation/annihilation leg
       correctly signed, and does a complex-Hermitian-closed declaration
       (PairHop) carry its conjugate partner on the right slot.
    2. The frame map -- is the (0, 3, 4, 1, 2) pair-slot-swap applied
       consistently to BOTH the ED derivative and chi0 (an inconsistency
       there is invisible on R=0 for the diagonal-only types and shows up
       only where the table is genuinely off-diagonal, i.e. PairHop).
    3. Finite-difference conditioning -- is ``CAMPAIGN_V1`` too large
       (nonlinear contamination pollutes the Richardson estimate) or too
       small (floating-point cancellation dominates); does
       ``delta_rich`` (not ``delta_nmat``, which is structurally zero
       here -- the solver-side sandwich is exact in the coupling, no
       Matsubara truncation) dominate ``tol``.

    Only once all three are independently ruled out does a FAIL indict a
    production defect in ``bond_correlator_transverse`` itself -- the
    expected outcome, going in, is PASS (or PASS-ZERO) on all 7.
    """

    def _check(self, kind, phase=1.0, expected_status="PASS"):
        norb = oracle.NORB
        D_ed_v1, D_ed_vhalf = _h2_ed_derivative_solver_frame(kind, phase)
        chi0 = _h2_chi0_solver_frame()
        ham_pm = _h2_ham_pm_expected(kind, phase, norb=norb)
        # D_pred is EXACT in the coupling (chi0 is exact ED, ham_pm_expected
        # is the exact table value, no Matsubara truncation anywhere) --
        # richardson()'s own docstring: "the solver-side prediction is
        # exactly linear in the coupling and has no Richardson stencil".
        # Passing the SAME exact value at both "nmat" slots makes
        # adjudicate_granule's delta_nmat structurally 0, which is the
        # honest statement here (no solver-side approximation to converge).
        D_pred = -np.einsum('qij,jk,qkl->qil', chi0, ham_pm, chi0)
        # zero_mask is NOT "cells where ham_pm is zero": D_pred is a
        # MATRIX PRODUCT (chi0 @ ham_pm @ chi0), so even a single nonzero
        # ham_pm cell spreads dense signal across most (I, J) cells of
        # D_pred via chi0's own off-diagonal content -- a per-cell
        # ham_pm-sparsity mask flagged that spread as spurious "zero"
        # failures although D_pred and D_ed agreed on it to ~1e-9 (first
        # investigation-protocol pass, before this fix). The only case
        # where D_pred is STRUCTURALLY zero on every cell is ham_pm being
        # the all-zero matrix (Hund/PairLift): 0 sandwiched between any
        # chi0 is 0 everywhere, independent of chi0's structure. For the
        # five active types the whole (q, ND, ND) grid is "bearing" --
        # the full predicted matrix is compared, not just the table's
        # named cell, which is a STRONGER check than the named cell
        # alone.
        all_zero = bool(np.all(ham_pm == 0))
        zero_mask = np.full(ham_pm.shape, all_zero, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred, D_pred, zero_mask,
            "H2/{}".format(kind))
        self.assertEqual(
            rec["status"], expected_status,
            "granule {} status={} (expected {}; delta_rich={:.3e} "
            "delta_nmat={:.3e} tol={:.3e} max_signal={:.3e} "
            "first_failures={}) -- see class docstring's INVESTIGATION "
            "PROTOCOL before assigning fault".format(
                rec["label"], rec["status"], expected_status,
                rec["delta_rich"], rec["delta_nmat"], rec["tol"],
                rec["max_signal"], rec["failures"][:5]))

    def test_coulomb_intra(self):
        self._check("CoulombIntra")

    def test_coulomb_inter(self):
        self._check("CoulombInter")

    def test_ising(self):
        self._check("Ising")

    def test_exchange(self):
        self._check("Exchange")

    def test_hund(self):
        self._check("Hund", expected_status="PASS-ZERO")

    def test_pair_lift(self):
        self._check("PairLift", expected_status="PASS-ZERO")

    def test_pair_hop_complex(self):
        # A genuinely complex P: PairHop is the one on-site type whose
        # ham_pm is NOT symmetric under the pair-slot swap ((0,1,1,0) and
        # (1,0,0,1) carry conjugate, not equal, values), so this is the
        # case that would actually catch a frame-map mistake on-site --
        # mirrors H1/TestTransverseComplexPairHop's choice of PairHop for
        # the same reason.
        self._check("PairHop", phase=oracle.PHASE)


if __name__ == "__main__":
    unittest.main()
