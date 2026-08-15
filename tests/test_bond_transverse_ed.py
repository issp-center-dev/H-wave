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

- H3: the transverse ED-to-solver bond FRAME MAP, derived (not assumed)
  from the ED operator's own Wick contraction, and pinned at V=0 against
  a free transverse bond bubble composed IN THIS TEST MODULE from the
  production kernel primitives (``bubble._prepare_dense``,
  ``bubble._bond_pair_full_block``, PRs #152/#155) -- the test-local
  precursor of Phase W's ``transverse_bond_bubble_static`` -- via
  ``TestH3FrameMapVZero`` and the module-level
  ``transverse_ed_to_solver_map``.

Tests must be run from the repository root.
"""

import functools
import unittest

import numpy as np

import tests.test_rpa_vs_ed_oracle as oracle
from hwave.solver import bubble as _hbubble
from hwave.solver import green as _hgreen
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

    def test_ising_hund_term_builders_match_dense_hamiltonian(self):
        """Standalone Hamiltonian-level pin for ``_terms_ising_hund``
        (fix round 1, review finding: the other 5 types are backed this
        way by ``oracle._h_int``, but Ising/Hund -- the two types this
        module had to build locally -- had no such pin).

        Builds each dense many-body Hamiltonian matrix DIRECTLY from the
        physical operator definitions -- density operators
        (``CD[mode] @ C[mode]``) combined by plain matrix
        multiplication/addition, the SAME construction style
        ``oracle._h_int`` uses for its own types (e.g. its CoulombIntra:
        ``v * (_n(j, 0, 0) @ _n(j, 0, 1))``) -- rather than through the
        generic ``(p, q, r, s, coeff)`` quartic engine
        (``ed_oracle_util.h_int_from_terms``) that
        ``_terms_ising_hund``'s OWN output is fed through elsewhere in
        this module. This is therefore an INDEPENDENT check of the
        formula: it never calls ``_terms_ising_hund`` to build the
        Hamiltonian, only to build the side being checked.

        Convention source (src/hwave/solver/uhfk.py):

        - Ising (``_make_ham_inter``'s Ising block comment): "the
          documented Hamiltonian is J (n_up - n_down)(n_up - n_down)" --
          density operators are diagonal (they commute), so this dense
          form and ``_terms_ising_hund``'s four explicit
          (+,-,-,+)-signed same-site density-density terms are the same
          operator by construction, checked here to round-off.
        - Hund (``_make_ham_inter``'s Hund block comment): "interaction
          coeffs : -J^Hund by convention", with the append machinery
          populating only the SAME-SPIN slots (spin_table entries
          ``(0,0,0,0)``/``(1,1,1,1)`` -- see also rpa.py's
          ``_append_inter``, which drives the actual solver's assembly
          off the identical same-spin-only rule): H_Hund =
          -J * (n_{0,up} n_{1,up} + n_{0,dn} n_{1,dn}).
        """
        fx, C, CD, _H1 = oracle._fx2_state()
        v = 0.37 + 0.0j   # arbitrary nonzero coupling; both Hamiltonians
                          # are real-valued by construction (pure density
                          # operator products)

        def _n(j, o, s):
            m = fx.mode(j, o, s)
            return CD[m] @ C[m]

        def _dense_ising(v):
            H = np.zeros((fx.dim, fx.dim), dtype=complex)
            for j in range(fx.L):
                sz0 = _n(j, 0, 0) - _n(j, 0, 1)
                sz1 = _n(j, 1, 0) - _n(j, 1, 1)
                H = H + v * (sz0 @ sz1)
            return H

        def _dense_hund(v):
            H = np.zeros((fx.dim, fx.dim), dtype=complex)
            for j in range(fx.L):
                H = H - v * (_n(j, 0, 0) @ _n(j, 1, 0)
                             + _n(j, 0, 1) @ _n(j, 1, 1))
            return H

        for kind, dense_fn in (("Ising", _dense_ising),
                                ("Hund", _dense_hund)):
            with self.subTest(kind=kind):
                terms = _terms_ising_hund(fx, kind, v)
                got = ed_oracle_util.h_int_from_terms(fx, terms)
                want = dense_fn(v)
                assert_approx_array(got, want, rel=0, abs=1e-13)

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


# ---------------------------------------------------------------------------
# H3: the transverse ED-to-solver bond frame map, derived from the ED
# operator's own Wick contraction, and pinned at V=0 against a free
# transverse bond bubble composed IN TEST CODE from the production kernel
# primitives (bubble._prepare_dense / bubble._bond_pair_full_block,
# PRs #152/#155) -- the test-local precursor of Phase W's
# transverse_bond_bubble_static.
# ---------------------------------------------------------------------------


def transverse_ed_to_solver_map(channels, delta_r, norb):
    """DERIVED index map from ``SectorED.bond_correlator_transverse``'s
    ``(q, I, J)`` frame to the solver's enlarged transverse-bond-bubble
    slot convention ``m*norb**2 + l1*norb + l2``, where ``(l1, l2)`` is
    the ROW/COLUMN orbital sub-index that
    ``bubble._bond_pair_full_block``'s own ``contract_general`` assembly
    produces INSIDE one ``(m, m')`` channel-pair block -- but only once
    the FORWARD Green factor is fixed to the DOWN block and the REVERSED
    factor to the UP block (see ``_transverse_bond_bubble_free``'s
    docstring for the companion half of this derivation; neither half is
    meaningful without the other, exactly as #151's
    ``ed_to_solver_bond_map`` and ``bond_bubble``'s own ``m``
    convention are a matched pair).

    Derivation (operator-level, mirroring #151's
    ``ed_to_solver_bond_map`` derivation comment in
    tests/test_bond_vs_ed_oracle.py -- read there first).
    ``SectorED.bond_correlator_transverse``'s channel operator is

        S+_{R,a,b}(q) = sum_j e^{+iqj} c^+_{j+R,a,up} c_{j,b,dn}

    which is LITERALLY ``bond_correlator``'s
    ``B_{R,a,b}(q) = sum_j e^{iqj} c^+_{j+R,a} c_{j,b}`` with the
    creation leg fixed to spin up and the annihilation leg fixed to spin
    down -- the site-offset / momentum bookkeeping (the ``(j+R) mod L``
    placement, the ``e^{+iqj}`` phase) is therefore IDENTICAL IN FORM to
    the longitudinal case; only the SPIN FLAVOR each leg carries differs.
    Wick-contracting
    ``<T_tau S+_{R,a,b}(q,tau) S+_{R',c,d}(q)^dagger(0)>`` (only
    up-up and down-down pairings survive, since the free Hamiltonian
    never mixes spin):

        <T_tau c^+_{j+R,a,up}(tau) c_{l+R',c,up}(0)>   -- UP: both legs
            are CREATION operators of the two S+ insertions, so this is
            the REVERSED time order relative to the standard
            G(tau) = -<T c(tau) c^+(0)> convention (creation at tau,
            annihilation at 0) -- the role
            ``kgrid.reverse_fft_axes``'s time/space-reversal trick
            (``prepped.green_rev``) exists to realize;
        <T_tau c_{j,b,dn}(tau) c^+_{l,d,dn}(0)>          -- DOWN: the
            STANDARD forward tau order (annihilation at tau, creation
            at 0) -- the ``prepped.green_rt`` (``green_fwd``) role.

    So the DOWN propagator plays ``_bond_pair_full_block``'s FORWARD
    role and the UP propagator plays its REVERSED role -- the OPPOSITE
    of the naive per-block reading of ``transverse_bubble``'s own
    docstring ("block 0 is G_up ... forward"): that docstring fixes its
    OWN component (a different, R=0-only cross-block object; nothing in
    Phase H or the on-site table pins its fwd/rev choice against THIS
    one), so it does not by itself settle the question here -- this
    module derives the assignment fresh from the ED operator, per the
    task's requirement, and the ``TestH3FrameMapVZero`` pin below
    confirms it numerically (an unswapped ``fwd=up, rev=down`` choice
    fails the pin by O(1), not by a finite-Nmat amount -- see that
    class's docstring).

    Within one ``(m, m')`` channel-pair block,
    ``bubble.contract_general``'s row/column split is
    ``chi0[a',c',b',d'] = g_fwd[a',b'] * g_rev[d',c']``, row = ``(a',c')``,
    column = ``(b',d')``. Matching the Wick contraction above with
    ``g_fwd = G_dn``, ``g_rev = G_up``: the DOWN factor's indices are
    (row channel's ED annihilation orbital ``b_ED``, column channel's ED
    annihilation orbital ``d_ED``) -> ``a' = b_ED``, ``b' = d_ED``; the
    UP factor's indices are (column channel's ED creation orbital
    ``c_ED``, row channel's ED creation orbital ``a_ED``) -> the
    REVERSED factor's first index ``d' = c_ED``, second index
    ``c' = a_ED``. So the ROW pair is ``(a', c') = (b_ED, a_ED)`` --
    ``row_natural = b_ED*norb + a_ED``, an ORBITAL SWAP relative to the
    ED channel's own ``(a, b) = (creation, annihilation)`` ordering --
    structurally the SAME KIND of swap ``ed_to_solver_bond_map``
    established longitudinally (a pair-internal swap, not a channel
    swap), but note the bond-CHANNEL index ``m`` itself is NOT reversed:
    ``delta_r[m]`` enters ``_transverse_bond_bubble_free``'s roll as
    ``+delta_r[m] - delta_r[mp]``, the SAME sign
    ``_bond_pair_full_block``'s established (Pin-2b-verified)
    longitudinal convention uses, with no additional sign flip -- the
    transverse channel/R bookkeeping is UNCHANGED from the longitudinal
    answer; what is genuinely NEW here is the SPIN-BLOCK role swap
    (fwd=down, rev=up) derived above, which has no longitudinal analogue
    (the longitudinal bubble has only one block).

    R=0 reduction (required by the task): at a single channel
    ``delta_r = [0]``, ``slot_to_channel[l1*norb+l2]`` is the index of
    ED channel ``(0, l2, l1)`` -- i.e. ``raw[:, smap][:, :, smap]``
    relabels BOTH the row and column orbital axes by the SAME swap
    permutation ``(l1, l2) -> (l2, l1)``. Because ``raw`` is Hermitian at
    fixed q for ANY channel list (``_lehmann_dagger``'s kernel is real
    and symmetric -- the same structural fact ``pair_correlator``'s
    docstring notes for its own ``Xpp[q] = Xpp[q]^dagger`` check), its
    plain transpose equals its elementwise conjugate:
    ``raw[q,J,I] = conj(raw[q,I,J])``. H2's already-established
    ``_frame_map`` (the Kubo pair-slot-swap, verified there against the
    bespoke ``_chi_pm_ed`` reference) is EXACTLY this transpose
    (``_frame_map(x)[q,J,I] = x[q,I,J]``, confirmed by direct index
    algebra: reshaping ``(q,ND,ND)`` to ``(q,l1,l2,l1',l2')`` and
    transposing ``(0,3,4,1,2)`` swaps the flat row/column labels without
    touching their internal order). ``TestH3FrameMapVZero`` checks
    NUMERICALLY (not just by this argument) that
    ``transverse_ed_to_solver_map``'s single-channel ``smap`` reproduces
    H2's ``_frame_map`` to round-off on the H1/H2 fixture -- confirmed
    (max diff ~2e-16) -- so this map's R=0 sub-map is EXACTLY the
    established Kubo pair-slot-swap, as required.

    This is a derivation, cross-checked empirically by
    ``TestH3FrameMapVZero`` below (not searched for) -- see that class's
    docstring for the measured raw/extrapolated distances and the
    mutation-check finding.

    Parameters
    ----------
    channels : list of (R, a, b) tuples (``(R,)`` at ``norb == 1``) --
        the EXACT ordered list to pass to ``bond_correlator_transverse``.
        Must contain, for every ``delta_r[m]`` and every physical orbital
        pair ``(l1, l2)``, the SWAPPED entry ``(delta_r[m], l2, l1)``.
    delta_r : sequence of int, length B
        The declared 1-D ring channel displacements (``R = delta_r[m]``;
        y/z components are out of scope, mirroring
        ``ed_to_solver_bond_map``'s own scope).
    norb : int

    Returns
    -------
    slot_to_channel : ndarray of int, shape (B * norb**2,)
        ``slot_to_channel[m*norb**2 + l1*norb + l2]`` is the index into
        ``channels`` of ED channel ``(delta_r[m], l2, l1)``.
    """
    nd = norb * norb
    B = len(delta_r)
    index_of = {ch: i for i, ch in enumerate(channels)}
    slot_to_channel = np.full(B * nd, -1, dtype=int)
    for m in range(B):
        R = int(delta_r[m])
        for l1 in range(norb):
            for l2 in range(norb):
                idx = l1 * norb + l2
                key = (R,) if norb == 1 else (R, l2, l1)
                if key not in index_of:
                    raise ValueError(
                        "transverse_ed_to_solver_map: channels is missing "
                        "{!r} (needed for solver slot m={}, (l1,l2)=({},{}))"
                        .format(key, m, l1, l2))
                slot_to_channel[m * nd + idx] = index_of[key]
    return slot_to_channel


def _h3_zeeman_h1(fx, hz):
    """``fx.build_h1()`` plus a per-site, per-orbital ZEEMAN split
    (``+hz`` on spin-up modes -- mode parity 0, ``-hz`` on spin-down --
    parity 1) added on the diagonal, so ``G_up != G_dn``. Required: a
    spin-SYMMETRIC fixture cannot distinguish ``fwd=up,rev=down`` from
    ``fwd=down,rev=up`` (the two Green functions would be numerically
    identical), which would make the frame-map pin below vacuous on the
    one convention it most needs to catch."""
    h1 = fx.build_h1().copy()
    for j in range(fx.L):
        for o in range(fx.norb):
            h1[fx.mode(j, o, 0), fx.mode(j, o, 0)] += hz
            h1[fx.mode(j, o, 1), fx.mode(j, o, 1)] -= hz
    return h1


def _h3_free_two_block_green(fx, hz, nmat):
    """The canonical two-block ``(2, nmat, fx.L, norb, norb)`` free
    Matsubara Green function for ``fx``'s ring, block 0 = spin up
    (``H0(k) + hz``), block 1 = spin down (``H0(k) - hz``) -- the SAME
    Zeeman split ``_h3_zeeman_h1`` adds to the ED Hamiltonian, so the
    solver-side test bubble and the ED reference diagonalize the
    IDENTICAL one-body Hamiltonian. Built through
    ``hwave.solver.green.build_green`` per spin block (brief
    requirement), diagonalizing ``H0(k) +- hz`` via ``np.linalg.eigh``
    per k-point first. ``coeff_tail=0.0`` (no tail correction): the RAW
    bond bubble, O(1/Nmat) convergence -- matching
    tests/test_bond_vs_ed_oracle.py's Pin 2b convention
    (``_richardson_nmat(..., order=1)``) that this pin's own
    Richardson-in-Nmat recipe reuses."""
    L, norb = fx.L, fx.norb
    tmat = np.zeros((norb, norb), dtype=complex)
    for (a, b), v in fx.t.items():
        tmat[a, b] = v
    kx = 2.0 * np.pi * np.arange(L) / L
    phase = np.exp(1j * kx)
    eps_k = (tmat[None, :, :] * phase[:, None, None]
             + tmat.conj().T[None, :, :] * np.conj(phase)[:, None, None])
    eps_k = eps_k + np.diag(np.asarray(fx.eps, dtype=complex))[None, :, :]
    ident = np.eye(norb, dtype=complex)
    eigvals = np.zeros((2, L, norb))
    eigvecs = np.zeros((2, L, norb, norb), dtype=complex)
    for blk, sign in ((0, +1.0), (1, -1.0)):
        hk = eps_k + (sign * hz) * ident[None, :, :]
        ev, V = np.linalg.eigh(hk)
        eigvals[blk] = ev
        eigvecs[blk] = V
    full_kw, _deflated_kw, _tail = _hgreen.build_green(
        eigvals, eigvecs, fx.mu, fx.beta, nmat, 0.0, want_full=True)
    return full_kw


def _transverse_bond_bubble_free(green_kw, beta, delta_r, norb,
                                 spatial_shape, workers=None,
                                 mutate_pair=None):
    """``(nvol, ND, ND)`` static (Omega=0) free transverse bond bubble,
    ``ND = B*norb**2`` -- composed IN TEST CODE from the production
    primitives ``bubble._prepare_dense`` / ``bubble._bond_pair_full_block``
    (PRs #152/#155): one call to ``_prepare_dense``, then per ``(m, m')``
    channel pair, roll the REVERSED block's real-space array by
    ``delta_r[m] - delta_r[mp]`` (the ``_bond_pair_full_block`` /
    ``bond_bubble_static`` roll pattern) and contract. This is the
    test-local precursor of Phase W's ``transverse_bond_bubble_static``.

    FORWARD block = 1 (spin DOWN), REVERSED block = 0 (spin UP) -- the
    DERIVED (not assumed) block assignment; see
    ``transverse_ed_to_solver_map``'s docstring for the Wick-contraction
    argument this composition realizes.

    ``green_kw``: canonical two-block ``(2, nmat, nvol, norb, norb)``,
    block 0 = G_up, block 1 = G_dn (``transverse_bubble``'s own
    labeling). No ``green0_tail`` argument -- RAW bond bubble, O(1/Nmat)
    convergence.

    ``mutate_pair``: ``None`` (normal operation) or an ``(m, mp)`` pair
    whose roll shift is NEGATED -- the Step-4 mutation-check hook,
    documented and exercised manually (not as a committed automated
    toggle; see ``TestH3FrameMapVZero``'s docstring for the measured
    failure magnitude and the discipline this follows,
    tests/test_bubble_kernel.py's established mutation-check pattern:
    edit, measure, revert, ``git diff`` clean)."""
    if green_kw.shape[0] != 2:
        raise ValueError(
            "_transverse_bond_bubble_free: green_kw must carry exactly 2 "
            "blocks (spin up, spin down), got nblock={}".format(
                green_kw.shape[0]))
    prepped = _hbubble._prepare_dense(green_kw, None, beta, spatial_shape,
                                      workers)
    _nblock, nmat, nvol, _nd, _ = green_kw.shape
    npair = norb * norb
    B = len(delta_r)
    ND = B * npair
    sgn = prepped.sgn.reshape((nmat, 1, 1, 1))
    green_fwd_sgn = prepped.green_rt[1] * sgn      # forward = block 1 (down)
    green_rev = prepped.green_rev[0]               # reversed = block 0 (up)
    out = np.zeros((nvol, ND, ND), dtype=complex)
    static_index = nmat // 2
    for m in range(B):
        dr_m = int(delta_r[m])
        for mp in range(B):
            dr_mp = int(delta_r[mp])
            s = dr_m - dr_mp
            if mutate_pair is not None and (m, mp) == mutate_pair:
                s = -s
            block = _hbubble._bond_pair_full_block(
                green_fwd_sgn, green_rev, None, beta, spatial_shape,
                workers, (s, 0, 0))
            out[:, m * npair:(m + 1) * npair,
                   mp * npair:(mp + 1) * npair] = block[static_index]
    return out


def _h3_richardson_nmat(fn, n1):
    """Two-point Richardson extrapolation of ``fn(nmat)`` as
    ``nmat -> infinity``, assuming the leading finite-Nmat error is
    ``O(1/nmat)`` (order=1 -- the RAW, no-tail-correction Matsubara-sum
    convergence rate; same formula and same rate as
    tests/test_bond_vs_ed_oracle.py's ``_richardson_nmat(fn, n1, order=1)``,
    reproduced here rather than imported to keep this module's coupling
    to that one at the level Task 1/2 already established (importing it
    as ``oracle``-style would pull in that module's own heavy fixtures
    for no benefit here)."""
    f1 = fn(n1)
    f2 = fn(2 * n1)
    return 2 * f2 - f1


class TestH3FrameMapVZero(unittest.TestCase):
    """H3 (spec: "Calibration gates" -- H3): the V=0 pin that catches a
    wrong bond reversal or phase convention, which H1/H2 (on-site only)
    cannot. ``bond_correlator_transverse`` at V=0 (a free-fermion
    fixture, no interaction terms) is compared, through
    ``transverse_ed_to_solver_map``, against
    ``_transverse_bond_bubble_free`` -- the test-local composed free
    transverse bond bubble -- on an ODD ring (L=5), COMPLEX hopping,
    NONZERO-R channels including OFF-DIAGONAL ``(m, m')`` cells, plus a
    multi-orbital (norb=2) orbital-swap variant.

    Richardson-in-Nmat (base pair 512/1024, order=1 -- see
    ``_h3_richardson_nmat``): the EXTRAPOLATED distance is asserted
    against ``assert_approx_array(rel=0, abs=2e-9)``. The spec calls for
    "~1e-9"; the measured extrapolated distances below are 8.5e-10
    (norb=1 fixture) and 1.23e-9 (norb=2 fixture) -- the norb=2 case sits
    marginally ABOVE a literal 1e-9, so the ceiling used here is 2e-9 (a
    real, still-tight margin over both measurements, not folded down to
    make a stricter number look achieved). RAW (single Nmat=1024, no
    extrapolation) distances, recorded separately per the #151 Pin-2
    pattern and NEVER folded into the ceiling above:

    - norb=1 fixture (L=5 ring, t=0.7*e^{0.3j}, beta=2, mu=0.2,
      hz=0.15, channels R in {0, +-1, +-2}): raw distance at Nmat=1024
      is 1.979e-4, scaling as C/Nmat with C ~ 0.2026 (measured at
      Nmat=512/1024/2048: 3.958e-4/1.979e-4/9.895e-5, each half the
      previous -- clean O(1/Nmat)); Richardson(512, 1024) extrapolated
      distance 8.49e-10.
    - norb=2 fixture (L=3 ring -- see scope note below -- t including a
      genuine inter-orbital hop, beta=2.5, mu=0.1, hz=0.12, channels R
      in {0, +-1} x orbital pairs {(0,0),(1,1),(0,1),(1,0)}): raw
      distance at Nmat=1024 is 2.474e-4, C ~ 0.2533; Richardson(512,
      1024) extrapolated distance 1.225e-9.

    Scope note on the norb=2 fixture's ring size: the spec's fx5T-free
    description names L=5 for both variants, but ``SectorED`` at L=5,
    norb=2 (nmode=20, largest (Nup,Ndn) sector dim 252, ~20 channels x 5
    q built per test) did not complete in a reasonable time (observed:
    did not return within 2 minutes on a build that returns in seconds
    at L=3) -- the same L=3/norb=2 tractability boundary the spec's own
    Phase W ED-granule section uses for ITS multi-orbital case ("one
    multi-orbital case (L=3, norb=2, ...) for the orbital indexing").
    This module follows that established precedent rather than the
    fx5T-free description's literal L for the orbital-swap variant; the
    ring is still odd, complex-hopping, spinful, and Zeeman-split, and
    the orbital-swap content (the load-bearing NEW check beyond the
    norb=1 fixture) is unaffected by ring size.

    Mutation check (Step 4, manual -- edited, run, reverted; not a
    committed automated toggle, per tests/test_bubble_kernel.py's
    established pattern): negating one bond roll
    (``_transverse_bond_bubble_free``'s ``s = dr_m - dr_mp`` ->
    ``s = -(dr_m - dr_mp)``) at a single off-diagonal ``(m, mp)`` pair
    with nonzero relative shift, at Nmat=1024:

    - norb=1 fixture, pair (m=1, mp=3) (R=+1 vs R=+2, shift -1): distance
      jumps from the baseline 1.979e-4 to 6.038e-2 -- a ~305x increase.
    - norb=2 fixture, pair (m=0, mp=1) (R=0 vs R=+1, shift -1): distance
      jumps from the baseline 2.474e-4 to 2.033e-1 -- a ~822x increase.

    Both are confirmed FAILING at the ceiling used here (``abs=2e-9``)
    by roughly seven orders of magnitude, and at the spec's own "catches
    a wrong bond reversal or phase convention" bar by several orders of
    magnitude on top of that. The negation was applied only inside a
    throwaway probe script during development (never committed to this
    file); ``_transverse_bond_bubble_free`` as committed always uses
    ``s = dr_m - dr_mp`` unconditionally unless ``mutate_pair`` is
    explicitly passed (documented above, not exercised by any test in
    this module), so ``git diff`` on this file carries no trace of the
    mutation.
    """

    def test_r0_submap_reduces_to_kubo_pair_slot_swap(self):
        """The task's explicit reduction requirement: at a single
        R=0 channel, ``transverse_ed_to_solver_map``'s ``smap`` must
        reproduce the ALREADY-ESTABLISHED H1/H2 ``_frame_map`` (the
        Kubo pair-slot-swap, verified in H1/H2 against the bespoke
        ``_chi_pm_ed`` reference and the on-site table) to round-off, on
        the SAME H1/H2 fixture."""
        fx, _C, _CD, h1 = oracle._fx2_state()
        se = ed_oracle_util.SectorED(fx, h1=h1)
        norb = oracle.NORB
        channels = list(_h2_channels())     # [(0, a, b) for a, b in NORB]
        raw = se.bond_correlator_transverse(channels)

        smap = transverse_ed_to_solver_map(channels, [0], norb)
        mapped_h3 = raw[:, smap][:, :, smap]
        mapped_h2 = _frame_map(raw, norb)
        assert_approx_array(mapped_h3, mapped_h2, rel=0, abs=1e-12)

    def _check(self, fx, hz, channels, delta_r, norb, n1=512,
               extrapolated_abs=2e-9):
        h1 = _h3_zeeman_h1(fx, hz)
        se = ed_oracle_util.SectorED(fx, terms=(), h1=h1)
        raw = se.bond_correlator_transverse(channels)
        smap = transverse_ed_to_solver_map(channels, delta_r, norb)
        mapped = raw[:, smap][:, :, smap]

        def solver_of(nmat):
            green_kw = _h3_free_two_block_green(fx, hz, nmat)
            return _transverse_bond_bubble_free(
                green_kw, fx.beta, delta_r, norb, (fx.L, 1, 1))

        rich = _h3_richardson_nmat(solver_of, n1)
        assert_approx_array(rich, mapped, rel=0, abs=extrapolated_abs)

    def test_norb1_odd_ring_full_channel_grid(self):
        """L=5 ring (odd), norb=1, complex hopping, spinful with a
        Zeeman split; channels R in {0, +-1, +-2} -- the FULL 5x5
        channel grid, exercising every off-diagonal (m, m') cell on an
        odd ring where R and -R are genuinely distinct wrapped channels
        (R=+-1 -> {1,4}, R=+-2 -> {2,3})."""
        t = 0.7 * np.exp(0.3j)
        fx = ed_oracle_util.EDFixture(
            L=5, norb=1, t={(0, 0): t}, eps=(0.0,), T=0.5, mu=0.2)
        Rs = [0, 1, -1, 2, -2]
        channels = [(R,) for R in Rs]
        self._check(fx, 0.15, channels, Rs, 1)

    def test_norb2_orbital_swap(self):
        """L=3 ring (odd; see class docstring's scope note), norb=2,
        with a genuine inter-orbital complex hop; channels R in
        {0, +-1} x {(0,0), (1,1), (0,1), (1,0)} -- the (0,1)/(1,0)
        orbital pairs are the load-bearing content: they are the ones an
        un-swapped orbital map (row_natural = a_ED*norb+b_ED instead of
        b_ED*norb+a_ED) gets wrong even when the R/m bookkeeping is
        right."""
        t = {(0, 0): 0.6 + 0.2j, (1, 1): 0.4 - 0.1j,
             (0, 1): 0.15 + 0.05j, (1, 0): 0.15 + 0.05j}
        fx = ed_oracle_util.EDFixture(
            L=3, norb=2, t=t, eps=(0.05, -0.03), T=0.4, mu=0.1)
        Rs = [0, 1, -1]
        orb_pairs = [(0, 0), (1, 1), (0, 1), (1, 0)]
        channels = [(R, a, b) for R in Rs for (a, b) in orb_pairs]
        self._check(fx, 0.12, channels, Rs, 2)


if __name__ == "__main__":
    unittest.main()
