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

import collections
import functools
import types
import unittest

import numpy as np

import tests.test_rpa_vs_ed_oracle as oracle
from hwave.solver import bond_channels
from hwave.solver import bubble as _hbubble
from hwave.solver import green as _hgreen
from hwave.solver import rpa as rpa_mod
from hwave.solver.vertex_table import ring_spin_table
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
            # Phase-H deferred minor (folded into Phase W's Task 1): this
            # dense (dim, dim) matmul is warning-noisy under some BLAS
            # backends (spurious "invalid value encountered" on an
            # all-zero/near-degenerate product) -- the same class of
            # benign warning ed_oracle_util._diagonalize and
            # SectorED._build_sectors already guard with
            # np.errstate(all="ignore") around their own dense linear
            # algebra; the result here is exact integer/half-integer
            # arithmetic on a projector product, never actually NaN/Inf.
            with np.errstate(all="ignore"):
                for j in range(fx.L):
                    sz0 = _n(j, 0, 0) - _n(j, 0, 1)
                    sz1 = _n(j, 1, 0) - _n(j, 1, 1)
                    H = H + v * (sz0 @ sz1)
            return H

        def _dense_hund(v):
            H = np.zeros((fx.dim, fx.dim), dtype=complex)
            with np.errstate(all="ignore"):
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
        the SAME H1/H2 fixture.

        Phase-H deferred minor (folded into Phase W's Task 1, per the
        binding plan's Global Constraints): this ``smap == _frame_map``
        identity is SPIN-SYMMETRY-CONDITIONAL, not a general fact about
        ``transverse_ed_to_solver_map`` -- it holds here because the H1/H2
        fixture (``oracle._fx2_state()``) carries no Zeeman splitting
        (``G_up == G_dn``), the same condition ``_h3_zeeman_h1``'s own
        docstring names as the reason H3 needs an explicit Zeeman term at
        all ("a spin-SYMMETRIC fixture cannot distinguish fwd=up,rev=down
        from fwd=down,rev=up"). Using the spin-symmetric ``_frame_map``
        shortcut is correct HERE, on fx2, and nowhere else in this module
        (H3's own granules below use ``transverse_ed_to_solver_map``
        directly, never ``_frame_map``, precisely because they are
        Zeeman-split). Phase W's production/ED-oracle code imports the
        Zeeman-ROBUST ``transverse_ed_to_solver_map``, never the
        R=0-only, spin-symmetric ``_frame_map`` shortcut this test
        exists to calibrate."""
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


# ---------------------------------------------------------------------------
# Phase W, Task 1: Gate W0 -- the pre-implementation numeric oracle for the
# binding spec's derived W_{+-} element equations
# (docs/superpowers/specs/2026-08-15-bond-transverse-design.md, "The vertex
# -- element equations" + "Gate W0"). This section is TEST-LOCAL and
# PRE-IMPLEMENTATION: none of TransverseTopology/W_pm_bond/
# resolve_transverse_topology exist yet (Tasks 2-3 of the Phase-W plan build
# those AFTER this gate passes) -- `w_expected_from_records` is a standalone
# encoding of the spec's ordered-record equations, built from a lightweight
# test-local topology stand-in (`_TopoLike`) rather than the production
# `TransverseTopology` NamedTuple. A FAIL on any granule below reopens the
# spec section (the discrepancy protocol, #151): never adjust
# `w_expected_from_records` to force a PASS.
# ---------------------------------------------------------------------------

_TopoLike = collections.namedtuple("_TopoLike", ["delta_r", "reverse"])


def _reversal_orbit_representatives(delta_r, reverse, norb):
    """Test-local replica of the spec's shared `iter_reversal_orbits`
    helper (spec, "Canonical orbit rule"): enumerates ONE representative
    ``(m, a, b)`` per reversal orbit ``{(m, a, b), (reverse[m], b, a)}``,
    restricted to the OFF-SITE construction domain (``m != 0`` --
    "steps 2-3 iterate OFF-SITE (R != 0) content EXCLUSIVELY", spec
    "Construction domain"), tuple-minimum representative (matching the
    shared helper's stated rule). Channel 0 (on-site) never
    participates here -- its content lives exclusively in the
    caller-supplied ``onsite_ham_pm`` block of `w_expected_from_records`.
    """
    B = len(delta_r)
    seen = set()
    reps = []
    for m in range(1, B):
        for a in range(norb):
            for b in range(norb):
                key = (m, a, b)
                if key in seen:
                    continue
                partner = (int(reverse[m]), b, a)
                seen.add(key)
                seen.add(partner)
                reps.append(min(key, partner))
    return reps


def _w0_q_mesh_flat(q_mesh):
    """``(nvol, 3)`` flattened q array, ``q_d = 2*pi*n_d/N_d`` with
    C-order flattening -- the SAME convention the spec's "Gate W0"
    section fixes ("q runs on the reciprocal mesh... the SAME
    convention sc._build_bond_m0_blocks' phase uses") and
    ``bond_channels.make_bond_kernel_parts``'s own bond-phase
    construction (``PH[m] = exp(1j*(KX*dm[0]+KY*dm[1]+KZ*dm[2]))`` on
    ``meshgrid(kx, ky, kz, indexing='ij')``, reshaped C-order)
    already uses."""
    Nx, Ny, Nz = q_mesh
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    kz = 2.0 * np.pi * np.arange(Nz) / Nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    return np.stack([KX.ravel(), KY.ravel(), KZ.ravel()], axis=-1)


def w_expected_from_records(topo_like, declarations, q_mesh):
    """Test-local oracle: the spec's ordered-record ``W_{+-}`` equations
    (spec, "The vertex -- element equations"), encoded directly rather
    than via the not-yet-implemented production ``W_pm_bond``. This IS
    gate W0's reference against which this module's ED granules
    adjudicate the spec's derivation, BEFORE any of that production
    code exists.

    Parameters
    ----------
    topo_like : _TopoLike(delta_r, reverse)
        ``delta_r``: int array (B, 3), channel 0 ALWAYS (0, 0, 0).
        ``reverse``: int array (B,), ``reverse[reverse[m]] == m``,
        ``reverse[0] == 0``.
    declarations : dict
        ``declarations["onsite_ham_pm"]``: the ALREADY-ASSEMBLED
        ``(norb, norb, norb, norb)`` on-site vertex (spec step 1: "the
        RECEIVED ham_pm_onsite ... `W_pm_bond` performs NO assembly of
        its own here"), broadcast verbatim into the channel-0 block.
        ``declarations["CoulombInter"]``, ``declarations["Ising"]``,
        ``declarations["Exchange"]``: OFF-SITE-ONLY
        ``dict[m] -> (norb, norb) complex`` coefficient blocks,
        Hermitian-closed exactly like ``resolve_transverse_topology``'s
        ``coeffs`` contract (``coeffs[type][reverse[m]] ==
        conj(coeffs[type][m]).T``) -- entries need only be supplied at
        whichever channel ``_reversal_orbit_representatives`` picks as
        the representative's first element (closure lets the formula
        below recover the mirrored value without a second lookup; see
        the derivation in the docstring body below). Missing types /
        missing channels are treated as zero (PASS-ZERO controls, e.g.
        Hund/PairLift, need not appear at all).
    q_mesh : (Nx, Ny, Nz)
        Spatial shape supplying the q mesh (spec: "q runs on the
        reciprocal mesh q = 2*pi*(n_x/N_x, ...)").

    Returns
    -------
    W : ndarray, complex128, shape (nvol, ND, ND)
        ``ND = B * norb**2``, ``nvol = Nx*Ny*Nz``, C-order flattened --
        the canonical flattened convention every numerical object in
        the spec uses.

    Construction (spec steps 1-3, verbatim):

    1. Channel-0 block = the broadcast on-site ``ham_pm``
       (q-INDEPENDENT).
    2. Cross family (CoulombInter -Re, Ising +Re): for each off-site
       reversal-orbit representative ``(m, a, b)`` with declared
       coefficient ``C = declarations[type][m][a, b]``, BOTH mirrored
       diagonal TARGET cells get the SAME value ``s_type * Re(C)``,
       DIRECT placement (AMENDED 2026-08-16, Task-6 granule
       adjudication -- the second measurement-driven correction of this
       step; see the derivation note below)::

           W[q, (m,a,b), (m,a,b)]                   += s_type * Re(C)
           W[q, (reverse[m],b,a), (reverse[m],b,a)]  += s_type * Re(C)

       (the two target cells are the reversal orbit of the
       REPRESENTATIVE'S OWN channel/orbital pair -- i.e. the target IS
       the representative and its orbit partner, no additional
       ``m <-> reverse[m]`` swap; closure guarantees they carry the
       identical ``Re(.)`` value, so ONE lookup at the representative
       channel ``m`` suffices for both -- derivation below.)

       AMENDED (2026-08-16, Task-6 granule adjudication): the PREVIOUS
       form of this step placed the first entry at ``reverse[m]`` (same
       ``(a,b)``) and the second at ``m`` (swapped ``(b,a)``) -- an
       ``m <-> reverse[m]`` channel swap relative to the direct
       placement above. The multi-orbital off-site CoulombInter granule
       (L=3, norb=2, ``a != b`` -- the ONLY fixture class that can see
       this swap: at ``a == b`` the two placements coincide, which is
       exactly why Gate W0's norb=1 cross-family fixtures never caught
       it) FAILED the swapped placement at O(1) (measured 0.12 against
       ``tol=5.6e-4``) and PASSES, at 5.6e-5, the DIRECT placement above
       -- the SAME rule the ALREADY ED-validated longitudinal
       ``bare_bond_vertices`` Fock-diagonal insertion uses
       (``bond_channels.py``, ``I = m*nd + l1*norb + l2`` gets
       ``Re(v_bond[m][l1,l2])`` directly, no channel swap; #151 Finding
       1). See ``docs/superpowers/specs/2026-08-15-bond-transverse-
       design.md``, "Cross family", the 2026-08-16 amendment, for the
       full derivation and the coordinator ruling.
    3. Flip family (Exchange, ``f_J = -1``): for each off-site
       reversal-orbit representative ``(m, a, b)`` with
       ``J = declarations["Exchange"][m][a, b]`` at ``R = delta_r[m]``::

           W[q, (0,a,a), (0,b,b)] += -conj(J)  * exp(-i q.R)
           W[q, (0,b,b), (0,a,a)] += -J        * exp(+i q.R)

       AMENDED (2026-08-16, gate W0 adjudication): the design doc's
       original draft assigned ``-J*exp(+iqR)`` to the ``(aa,bb))``
       element; this module's own W0 granule (a) (multi-orbital
       off-site Exchange, complex J, non-self-inverse q) FAILED that
       assignment systematically (residuals 0.18-0.33 across every
       tested orientation and R sign) and PASSED, at the Richardson
       noise floor, the SWAPPED assignment above. The spec
       (docs/superpowers/specs/2026-08-15-bond-transverse-design.md,
       "Flip family", the "AMENDED (2026-08-16, GATE W0 ADJUDICATION)"
       paragraph) was updated accordingly: the un-conjugated ``J`` is
       the ``(bb)->(aa)`` scattering amplitude (row=out/col=in, the
       col leg daggered in ``<A;B^dagger>``), landing on
       ``W[(0,b,b),(0,a,a)]`` with phase ``exp(+iq.R)``; the Hermitian
       partner supplies the transposed element. The R=0 (on-site) limit
       is UNCHANGED (both forms coincide there after Hermitian closure
       -- precisely why H2's on-site-only readjudication could not by
       itself catch this, and the off-site granule was required).

    Derivation of step 2's "one lookup, both cells" shortcut (DIRECT
    placement, 2026-08-16 amendment): writing the formula generally as
    "target cell (T, x, y) gets ``s*Re(coeffs[T][x, y])``", the
    representative ``(m, a, b)`` feeds TWO target cells: ``(m, a, b)``
    [value ``s*Re(coeffs[m][a, b])``] and ``(reverse[m], b, a)`` [value
    ``s*Re(coeffs[reverse[m]][b, a])``]. Closure
    (``coeffs[reverse[m]][x, y] = conj(coeffs[m][y, x])``) gives
    ``coeffs[reverse[m]][b, a] = conj(coeffs[m][a, b])``, whose real
    part equals ``Re(coeffs[m][a, b])`` -- so both target cells carry
    the SAME value, computed from the ONE representative-channel
    lookup (the lookup itself is unchanged by the amendment; only WHICH
    cell each of the two values is written to changed).
    """
    delta_r = np.asarray(topo_like.delta_r, dtype=int)
    reverse = np.asarray(topo_like.reverse, dtype=int)
    B = delta_r.shape[0]
    if delta_r.shape != (B, 3):
        raise ValueError(
            "w_expected_from_records: delta_r must have shape (B, 3), got "
            "{}".format(delta_r.shape))
    if reverse.shape != (B,):
        raise ValueError(
            "w_expected_from_records: reverse must have shape (B,), got "
            "{}".format(reverse.shape))
    if not np.array_equal(delta_r[0], np.zeros(3, dtype=int)):
        raise ValueError(
            "w_expected_from_records: channel 0 must be R=(0,0,0), got "
            "{}".format(delta_r[0]))
    if not np.array_equal(reverse[reverse], np.arange(B)):
        raise ValueError(
            "w_expected_from_records: reverse must be an involution "
            "(reverse[reverse[m]] == m for every m)")
    if int(reverse[0]) != 0:
        raise ValueError(
            "w_expected_from_records: reverse[0] must be 0 (channel 0 is "
            "its own reversal partner)")

    onsite_ham_pm = np.asarray(declarations["onsite_ham_pm"], dtype=complex)
    norb = onsite_ham_pm.shape[0]
    if onsite_ham_pm.shape != (norb, norb, norb, norb):
        raise ValueError(
            "w_expected_from_records: onsite_ham_pm must have shape "
            "(norb, norb, norb, norb), got {}".format(onsite_ham_pm.shape))
    nd = norb * norb
    ND = B * nd

    q_flat = _w0_q_mesh_flat(q_mesh)         # (nvol, 3)
    nvol = q_flat.shape[0]

    W = np.zeros((nvol, ND, ND), dtype=complex)

    # Step 1: channel-0 block = the broadcast on-site ham_pm, q-independent.
    W[:, 0:nd, 0:nd] = onsite_ham_pm.reshape(nd, nd)[None, :, :]

    reps = _reversal_orbit_representatives(delta_r, reverse, norb)

    # Step 2: cross family (CoulombInter, Ising), bond-diagonal,
    # q-independent. DIRECT placement (AMENDED 2026-08-16, Task-6
    # granule adjudication -- see the docstring above): the target IS
    # the representative's own (m,a,b), and its orbit partner
    # (reverse[m],b,a) -- no additional m<->reverse[m] swap.
    for kind, s_type in (("CoulombInter", -1.0), ("Ising", 1.0)):
        coeffs_map = declarations.get(kind, {})
        for (m, a, b) in reps:
            block = coeffs_map.get(m)
            if block is None:
                continue
            val = s_type * float(np.real(np.asarray(block)[a, b]))
            if val == 0.0:
                continue
            idx1 = m * nd + a * norb + b
            W[:, idx1, idx1] += val
            target2 = int(reverse[m])
            idx2 = target2 * nd + b * norb + a
            W[:, idx2, idx2] += val

    # Step 3: flip family (Exchange), the two ordered records,
    # q-dependent ONLY through this local-channel phase (spec:
    # "W_pm_bond is q-dependent ONLY through the flip-family
    # local-channel entries"). AMENDED (2026-08-16, gate W0
    # adjudication -- see the docstring above): the un-conjugated J is
    # the (bb)->(aa) scattering amplitude with phase exp(+iq.R); the
    # (aa,bb) element carries conj(J)*exp(-iq.R). The draft's swapped
    # assignment FAILED W0's off-site multi-orbital Exchange granule
    # (residuals 0.18-0.33); this assignment PASSES at the Richardson
    # noise floor.
    exch_map = declarations.get("Exchange", {})
    for (m, a, b) in reps:
        block = exch_map.get(m)
        if block is None:
            continue
        J = complex(np.asarray(block)[a, b])
        if J == 0.0:
            continue
        R = delta_r[m].astype(float)
        phase = np.exp(1j * (q_flat @ R))            # exp(+i q.R), (nvol,)
        idx_aa = a * norb + a
        idx_bb = b * norb + b
        W[:, idx_aa, idx_bb] += -1.0 * np.conj(J) * np.conj(phase)
        W[:, idx_bb, idx_aa] += -1.0 * J * phase

    return W


class TestW0OrderedRecordOracle(unittest.TestCase):
    """Structural self-tests for ``w_expected_from_records`` -- pure
    algebra, no ED involved (the ED-comparison granules are
    ``TestW0Granules`` below). Per Task 1 Step 1's brief: Hermiticity at
    every q, the two Exchange records pinned individually with a
    genuinely complex ``J`` at a non-self-inverse q, and the B=1 on-site
    reduction."""

    _TOPO_OFFSITE = _TopoLike(
        delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]),
        reverse=np.array([0, 2, 1]))

    def test_hermitian_at_every_q(self):
        norb = 2
        J = 0.6 - 0.35j
        declarations = {
            "onsite_ham_pm": np.zeros((norb, norb, norb, norb),
                                       dtype=complex),
            "Exchange": {1: np.array([[0, J], [0, 0]], dtype=complex)},
            "CoulombInter": {1: np.array([[0.3 + 0j, 0], [0, 0]],
                                          dtype=complex)},
            "Ising": {1: np.array([[0.2 + 0j, 0], [0, 0]], dtype=complex)},
        }
        W = w_expected_from_records(self._TOPO_OFFSITE, declarations,
                                     (5, 1, 1))
        got_dagger = np.conj(np.transpose(W, (0, 2, 1)))
        assert_approx_array(got_dagger, W, rel=0, abs=1e-12)

    def test_exchange_entries_pinned_individually(self):
        """A genuinely complex J at a non-self-inverse q (spec, "Gate
        W0": "the two Exchange entries SEPARATELY equal f_J*conj(J)*
        e^{-iqR} and f_J*J*e^{+iqR}" -- the AMENDED (2026-08-16, gate W0
        adjudication) assignment: the un-conjugated J is the (bb)->(aa)
        scattering amplitude, landing on W[(0,b,b),(0,a,a)] with phase
        e^{+iqR}; see `w_expected_from_records`'s own docstring for the
        full amendment note). L=5 (odd) -> every nonzero q is
        non-self-inverse; q_idx=1 (q=2*pi/5) is used here."""
        norb = 2
        J = 0.6 - 0.35j
        declarations = {
            "onsite_ham_pm": np.zeros((norb, norb, norb, norb),
                                       dtype=complex),
            "Exchange": {1: np.array([[0, J], [0, 0]], dtype=complex)},
        }
        Nx = 5
        W = w_expected_from_records(self._TOPO_OFFSITE, declarations,
                                     (Nx, 1, 1))
        q_flat = _w0_q_mesh_flat((Nx, 1, 1))
        q_idx = 1
        q = q_flat[q_idx, 0]
        phase = np.exp(1j * q * 1.0)          # R = delta_r[1] = (1, 0, 0)
        idx_aa, idx_bb = 0, 3                 # channel 0, (0,0) and (1,1)
        assert_approx_array(
            np.array([W[q_idx, idx_aa, idx_bb]]),
            np.array([-np.conj(J) * np.conj(phase)]), rel=0, abs=1e-13)
        assert_approx_array(
            np.array([W[q_idx, idx_bb, idx_aa]]), np.array([-J * phase]),
            rel=0, abs=1e-13)

    def test_b1_reduction_equals_broadcast_onsite_ham_pm(self):
        """Spec step 1: "B=1 reduction: on-site-only declarations give
        B=1 and W_pm_bond IS the broadcast ham_pm -- gate W1's algebraic
        basis." ``onsite_ham_pm`` here is H2's own already-ED-
        readjudicated table (`_h2_ham_pm_expected`, "reuse H2's spied
        builder" per the task brief) at a real CoulombIntra coupling."""
        norb = 2
        topo = _TopoLike(delta_r=np.array([[0, 0, 0]]),
                          reverse=np.array([0]))
        onsite = _h2_ham_pm_expected("CoulombIntra", 1.3, norb=norb)
        declarations = {
            "onsite_ham_pm": onsite.reshape(norb, norb, norb, norb)}
        W = w_expected_from_records(topo, declarations, (3, 2, 1))
        nvol = 6
        broadcast = np.broadcast_to(
            onsite[None, :, :], (nvol,) + onsite.shape)
        assert_approx_array(W, broadcast, rel=0, abs=1e-15)


# ---------------------------------------------------------------------------
# Task 1, Step 2: the three W0 ED granules. New off-site term builders
# (Exchange, Ising) not covered by any existing module -- canonical_
# density_terms (ed_oracle_util.py) already covers off-site CoulombInter
# (plain spin-summed density-density), and is reused directly below.
# ---------------------------------------------------------------------------

def _terms_exchange_offsite(fx, a, b, R, J):
    """Off-site Exchange term, generalizing
    ``tests.test_rpa_vs_ed_oracle._terms_for``'s ON-SITE Exchange
    operator convention (``c^+_{a,up} c_{b,up} c^+_{b,dn} c_{a,dn}``,
    already ED-validated on-site to produce ``ham_pm``'s ``-J`` entry --
    H2's ``test_exchange`` above re-confirms it through this module's
    own harness) to a nonzero bond displacement ``R``: the ``b``-orbital
    legs move to site ``j+R``, the ``a``-orbital legs stay at site
    ``j`` (a bond from site ``j``, orbital ``a``, to site ``j+R``,
    orbital ``b``). ``coeffs["Exchange"][m][a, b] = J`` fed to
    ``w_expected_from_records`` uses EXACTLY this same operator
    convention, so the ``f_J = -1`` coefficient-table sign is consistent
    between the ED side and the oracle side by construction (never
    independently re-derived here -- consistency, not an a-priori
    physical claim, is what gate W0 needs and what
    ``TestW0TermBuilderHamiltonianPins`` below pins at the Hamiltonian
    level). The Hermitian partner is the EXACT operator dagger
    (``conj(J) * (c^+_p c_q c^+_r c_s)^dagger = conj(J) * c^+_s c_r
    c^+_q c_p``), guaranteeing H is Hermitian by construction rather
    than by re-deriving a second, independent physical form for the
    reversed bond."""
    terms = []
    for j in range(fx.L):
        jr = (j + R) % fx.L
        p = fx.mode(j, a, 0)
        q = fx.mode(jr, b, 0)
        r = fx.mode(jr, b, 1)
        s = fx.mode(j, a, 1)
        terms.append((p, q, r, s, J))
        terms.append((s, r, q, p, np.conj(J)))
    return terms


def _terms_ising_offsite(fx, a, b, R, v):
    """Off-site Ising term, generalizing this module's own
    ``_terms_ising_hund``'s ON-SITE (R=0) Ising formula (uhfk.py's
    documented Hamiltonian ``J (n_up - n_down)(n_up - n_down)``) to a
    nonzero bond displacement ``R``: ``H = v * sum_j (n_{j,a,up} -
    n_{j,a,dn}) * (n_{j+R,b,up} - n_{j+R,b,dn})``. Manifestly Hermitian
    for real ``v`` without any mirrored (-R) declaration -- density
    operators at DIFFERENT sites commute, so
    ``(n_{j,a,up}-n_{j,a,dn}) @ (n_{j+R,b,up}-n_{j+R,b,dn})`` is
    Hermitian term by term, and the ``j``-sum over one period already
    visits every bond exactly once (the same reasoning
    ``canonical_density_terms``'s off-site branch relies on)."""
    terms = []
    for j in range(fx.L):
        jr = (j + R) % fx.L
        for su, su_sign in ((0, 1.0), (1, -1.0)):
            for sv, sv_sign in ((0, 1.0), (1, -1.0)):
                p = fx.mode(j, a, su)
                r = fx.mode(jr, b, sv)
                terms.append((p, p, r, r, v * su_sign * sv_sign))
    return terms


class TestW0TermBuilderHamiltonianPins(unittest.TestCase):
    """Hamiltonian-level dense pins (the H2 fix-round pattern -- see
    ``test_ising_hund_term_builders_match_dense_hamiltonian`` above) for
    the two NEW off-site term builders Task 1 needs: neither
    ``canonical_density_terms`` (CoulombInter-shaped, unsigned density-
    density only) nor this module's own on-site ``_terms_ising_hund``
    (R=0 only) covers off-site Exchange or off-site Ising. Each pin uses
    a SMALL tractable fixture (dense ``fx.annihilators()`` -- NOT
    SectorED, since this is a direct-operator-construction check)
    distinct from (and smaller than) the granule fixtures the term
    builders are actually exercised on in ``TestW0Granules`` below (case
    M's dense route is intractable -- see ``ed_oracle_util.SectorED``'s
    own docstring -- which is exactly why the granules use SectorED and
    this pin uses a separate, smaller fixture instead)."""

    def test_offsite_exchange_term_builder_matches_dense_hamiltonian(self):
        fx = ed_oracle_util.EDFixture(
            L=2, norb=2,
            t={(0, 0): 0.3, (1, 1): 0.2, (0, 1): 0.1, (1, 0): 0.1},
            eps=(0.05, -0.02), T=0.5, mu=0.1)
        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        J = 0.4 + 0.3j
        R, a, b = 1, 0, 1

        def _dense():
            H = np.zeros((fx.dim, fx.dim), dtype=complex)
            # np.errstate wrap: this dense (dim, dim) matmul is
            # warning-noisy under some BLAS backends (spurious
            # divide-by-zero/overflow/invalid-value warnings on a
            # near-nilpotent operator product), the same benign-warning
            # class the earlier Ising/Hund dense pin above already
            # guards.
            with np.errstate(all="ignore"):
                for j in range(fx.L):
                    jr = (j + R) % fx.L
                    op = (CD[fx.mode(j, a, 0)] @ C[fx.mode(jr, b, 0)]
                          @ CD[fx.mode(jr, b, 1)] @ C[fx.mode(j, a, 1)])
                    H = H + J * op + np.conj(J) * op.conj().T
            return H

        terms = _terms_exchange_offsite(fx, a, b, R, J)
        got = ed_oracle_util.h_int_from_terms(fx, terms)
        assert_approx_array(got, _dense(), rel=0, abs=1e-13)

    def test_offsite_ising_term_builder_matches_dense_hamiltonian(self):
        fx = ed_oracle_util.EDFixture(
            L=3, norb=1, t={(0, 0): 0.6 + 0.1j}, eps=(0.0,), T=0.4, mu=0.15)
        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        v = 0.37
        R, a, b = 1, 0, 0

        def _n(j, o, s):
            m = fx.mode(j, o, s)
            return CD[m] @ C[m]

        def _dense():
            H = np.zeros((fx.dim, fx.dim), dtype=complex)
            with np.errstate(all="ignore"):
                for j in range(fx.L):
                    jr = (j + R) % fx.L
                    sz_j = _n(j, a, 0) - _n(j, a, 1)
                    sz_jr = _n(jr, b, 0) - _n(jr, b, 1)
                    H = H + v * (sz_j @ sz_jr)
            return H

        terms = _terms_ising_offsite(fx, a, b, R, v)
        got = ed_oracle_util.h_int_from_terms(fx, terms)
        assert_approx_array(got, _dense(), rel=0, abs=1e-13)


def _w0_delta_r_1d(delta_r):
    """This module's granule fixtures are all 1-D rings (y, z always 0
    in the canonical ``(B, 3)`` ``delta_r``); ``bond_correlator_
    transverse``/``transverse_ed_to_solver_map`` (Phase H) both take the
    scalar ring offset ``R`` (their ``(R,)``/``(R, a, b)`` channel
    convention has no y/z component -- see
    ``transverse_ed_to_solver_map``'s own docstring: "y/z components are
    out of scope"), so this extracts just the x-component as a plain int
    list."""
    return [int(r) for r in np.asarray(delta_r)[:, 0]]


def _w0_channels_for(delta_r_1d, norb):
    """The full ``(R, a, b)`` channel list ``transverse_ed_to_solver_map``
    needs (its docstring: "Must contain, for every delta_r[m] and every
    physical orbital pair (l1, l2), the SWAPPED entry (delta_r[m], l2,
    l1)" -- ranging a, b over every combination already contains every
    swap)."""
    if norb == 1:
        return [(R,) for R in delta_r_1d]
    return [(R, a, b) for R in delta_r_1d for a in range(norb)
            for b in range(norb)]


def _w0_ed_derivative_solver_frame(fx, terms_of_v, delta_r, norb):
    """Richardson-extrapolated d(bond_correlator_transverse)/dv at
    v -> 0, full-minus-HF-only (the #151/H2 pattern), frame-mapped into
    the enlarged (``m*norb**2 + a*norb + b``) slot convention via
    ``transverse_ed_to_solver_map`` -- the Zeeman-robust map (this
    module's Phase-H deferred-minor clarification above:
    ``_frame_map`` is a spin-symmetric, R=0-only shortcut, not usable
    here since these granule fixtures have off-site channels, B>1).
    Returns ``(D_v1, D_vhalf)``, the two Richardson estimates
    ``adjudicate_granule`` needs."""
    delta_r_1d = _w0_delta_r_1d(delta_r)
    channels = _w0_channels_for(delta_r_1d, norb)
    smap = transverse_ed_to_solver_map(channels, delta_r_1d, norb)

    @functools.lru_cache(maxsize=None)
    def X_of_v(v):
        terms = terms_of_v(v)
        full = ed_oracle_util.SectorED(
            fx, terms=terms).bond_correlator_transverse(channels)
        hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
        hf_only = ed_oracle_util.SectorED(
            fx, terms=(), h1=hf_h1).bond_correlator_transverse(channels)
        diff = full - hf_only
        return diff[:, smap][:, :, smap]

    D_v1 = ed_oracle_util.richardson(X_of_v, CAMPAIGN_V1)
    D_vhalf = ed_oracle_util.richardson(X_of_v, CAMPAIGN_V1 / 2)
    return D_v1, D_vhalf


def _w0_chi0_solver_frame(fx, delta_r, norb):
    """The V=0 transverse correlator (exact ED, every declared channel),
    frame-mapped into solver space via ``transverse_ed_to_solver_map``
    -- gate W0's ``chi0`` ("chi0 evaluated on the FIXED noninteracting
    Hamiltonian", spec "Gate W0")."""
    delta_r_1d = _w0_delta_r_1d(delta_r)
    channels = _w0_channels_for(delta_r_1d, norb)
    smap = transverse_ed_to_solver_map(channels, delta_r_1d, norb)
    se = ed_oracle_util.SectorED(fx)
    raw = se.bond_correlator_transverse(channels)
    return raw[:, smap][:, :, smap]


def _w0_adjudicate_direction(fx, terms_of_v, topo, declarations, q_mesh,
                              label):
    """One granule direction, end to end: measure the ED-side
    derivative, build gate W0's oracle prediction
    ``-chi0 . w_expected_from_records(...) . chi0`` (chi0 EXACT ED, no
    Matsubara truncation anywhere on the predicted side -- so, exactly
    as H2 documents, ``D_pred`` is passed at BOTH the ``nmat`` and
    ``2*nmat`` slots of ``adjudicate_granule``, making ``delta_nmat``
    structurally 0, the honest statement when there is no solver-side
    approximation to converge), and adjudicate."""
    delta_r = np.asarray(topo.delta_r)
    reverse = np.asarray(topo.reverse)
    norb = np.asarray(declarations["onsite_ham_pm"]).shape[0]
    D_ed_v1, D_ed_vhalf = _w0_ed_derivative_solver_frame(
        fx, terms_of_v, delta_r, norb)
    chi0 = _w0_chi0_solver_frame(fx, delta_r, norb)
    W_expected = w_expected_from_records(topo, declarations, q_mesh)
    D_pred = -np.einsum('qij,qjk,qkl->qil', chi0, W_expected, chi0)
    # zero_mask: the same H2 argument -- D_pred is a MATRIX PRODUCT, so a
    # single nonzero W_expected cell spreads dense signal across most
    # (I, J) cells via chi0's own off-diagonal content. The only
    # structurally-zero case is W_expected being the all-zero matrix.
    all_zero = bool(np.all(W_expected == 0))
    zero_mask = np.full(W_expected.shape, all_zero, dtype=bool)
    return ed_oracle_util.adjudicate_granule(
        D_ed_v1, D_ed_vhalf, D_pred, D_pred, zero_mask, label)


def _w0_direction_set_sv_gate(records, label):
    """Test-local replica of ``ed_oracle_util.sensitivity_rank``'s SVD
    gate (``sv_ratio >= SENS_SV_FLOOR``, ``sigma_min >=
    100*delta_rich``) -- the spec's "Gate W0" campaign-level condition
    on top of each direction's own ``adjudicate_granule`` verdict
    ("the direction-set SVD gate ... whose sensitivity matrix passes
    the #151 SVD gate", spec "Gate W0"). A LOCAL replica, not
    ``ed_oracle_util.sensitivity_rank`` itself: that function's a-priori
    ``SENSITIVITY_EXPECTED_ACTIVE``/``SENSITIVITY_EXPECTED_NULL``
    registry is keyed to the #151 pp/ph campaign's own (fixture,
    channel) labels, and this task's Files list is
    ``tests/test_bond_transverse_ed.py`` alone -- ``ed_oracle_util.py``
    is not touched (nor does it need to be: the math below is the exact
    same formula, just without that unrelated campaign's label
    registry)."""
    names = sorted(records)
    for name in names:
        status = records[name]["status"]
        if status not in ("PASS", "PASS-ZERO"):
            raise ValueError(
                "_w0_direction_set_sv_gate[{}]: direction {!r} has status "
                "{!r} -- only PASS/PASS-ZERO granules feed the SVD gate"
                .format(label, name, status))
    active = [n for n in names if records[n]["status"] != "PASS-ZERO"]
    if not active:
        raise ValueError(
            "_w0_direction_set_sv_gate[{}]: every direction is PASS-ZERO "
            "-- no active column to rank".format(label))

    bearing_arrays = [np.asarray(records[n]["bearing_mask"], dtype=bool)
                       for n in active]
    union_mask = bearing_arrays[0].copy()
    for barr in bearing_arrays[1:]:
        union_mask |= barr
    n_rows = int(np.count_nonzero(union_mask))

    columns = []
    for n in active:
        pred = np.asarray(records[n]["pred_full"], dtype=float)
        bearing = np.asarray(records[n]["bearing_mask"], dtype=bool)
        columns.append(np.where(bearing, pred, 0.0)[union_mask])
    matrix = np.stack(columns, axis=1)
    n_cols = matrix.shape[1]
    delta_rich_max = max(records[n]["delta_rich"] for n in active)

    if n_rows < n_cols:
        sv = np.zeros(0)
        sigma_min, sigma_max, sv_ratio = 0.0, float("nan"), 0.0
    else:
        sv = np.linalg.svd(matrix, compute_uv=False)
        sigma_max, sigma_min = float(sv[0]), float(sv[-1])
        sv_ratio = sigma_min / sigma_max if sigma_max > 0.0 else 0.0

    print("W0 direction-set SVD[{}]: active={} rows={} cols={} "
          "singular_values={} sv_ratio={:.6e} sigma_min={:.6e} "
          "100*delta_rich_max={:.6e}".format(
              label, active, n_rows, n_cols,
              np.array2string(sv, precision=6), sv_ratio, sigma_min,
              100.0 * delta_rich_max))
    return dict(active=active, n_rows=n_rows, n_cols=n_cols,
                singular_values=sv, sv_ratio=sv_ratio, sigma_min=sigma_min,
                sigma_max=sigma_max, delta_rich_max=delta_rich_max)


class TestW0Granules(unittest.TestCase):
    """Gate W0's three ED granules (spec: "Gate W0" + the Phase-W plan's
    Task 1 Step 2): forward comparison of the measured
    d(bond_correlator_transverse)/dv (full-minus-HF-only, frame-mapped)
    against ``-chi0 . w_expected_from_records(...) . chi0``, chi0 = the
    exact V=0 ED correlator -- NO reconstruction, per the spec's binding
    procedure ("FORWARD comparison, no W reconstruction", spec "Gate
    W0"). Each granule runs 2 directions sharing one topology/fixture
    (CoulombInter/Ising: the g1 (R=+-1) and g2 (R=+-2) shells; Exchange:
    two independent complex phases of the SAME off-site bond), so the
    #151 direction-set SVD gate has a genuine (non-trivial-rank)
    sensitivity matrix to rank.

    CRITICAL: a FAIL on any of these is a genuine mismatch between the
    spec's derived W_{+-} equations and exact diagonalization -- STOP,
    do not commit, report verbatim (the discrepancy protocol, #151:
    re-check the encoding against the spec text, the frame map, the
    term builders via their Hamiltonian-level pins, Richardson
    conditioning -- never adjust ``w_expected_from_records`` to force a
    PASS).
    """

    def test_granule_a_multiorbital_offsite_exchange(self):
        """(a) L=3, norb=2, off-site Exchange at R=1 with a != b (a=0,
        b=1) and a genuinely complex J -- the norb=1 Exchange granule
        would collapse the two ordered elements W[(0,aa),(0,bb)] and
        W[(0,bb),(0,aa)] onto the same cell and could not adjudicate
        the record placement (spec, "ED granules"); this one checks
        BOTH independently via the full-grid adjudicate_granule
        comparison. Two directions ("real"/"imag": J = v*1 and
        J = v*1j) on the SAME topology give the SVD gate a rank-2
        matrix to check.

        ADJUDICATION OUTCOME (gate W0 working exactly as designed): the
        design doc's ORIGINAL DRAFT flip-family assignment
        (`W[(0,a,a),(0,b,b)] += -J*e^{+iqR}`,
        `W[(0,b,b),(0,a,a)] += -conj(J)*e^{-iqR}`) FAILED this granule
        systematically -- both directions, `status=FAIL`, residuals
        0.18-0.33 at every nonzero q (q=0 matched exactly, isolating the
        defect to the phase/orientation of the two records, not their
        magnitude or the term builder -- see the investigation recorded
        in `w_expected_from_records`'s docstring and
        docs/superpowers/specs/2026-08-15-bond-transverse-design.md's
        "AMENDED (2026-08-16, GATE W0 ADJUDICATION)" paragraph). The
        finding was adjudicated (coordinator-level) and the spec AMENDED
        to the swapped assignment now encoded in `w_expected_from_records`
        (`W[(0,a,a),(0,b,b)] += -conj(J)*e^{-iqR}`,
        `W[(0,b,b),(0,a,a)] += -J*e^{+iqR}`); with the amendment, this
        granule PASSES at the Richardson noise floor (both directions;
        see the module-level test-run output for the exact numbers)."""
        fx = ed_oracle_util.EDFixture(
            L=3, norb=2,
            t={(0, 0): 0.5 + 0.2j, (1, 1): 0.35 - 0.15j,
               (0, 1): 0.1 + 0.05j, (1, 0): 0.1 + 0.05j},
            eps=(0.05, -0.03), T=0.45, mu=0.1)
        topo = _TopoLike(delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0]]),
                          reverse=np.array([0, 2, 1]))
        q_mesh = (3, 1, 1)
        onsite_ham_pm = np.zeros((2, 2, 2, 2), dtype=complex)
        a, b, R = 0, 1, 1

        records = {}
        for name, phase in (("real", 1.0 + 0.0j), ("imag", 0.0 + 1.0j)):
            def terms_of_v(v, phase=phase):
                return _terms_exchange_offsite(fx, a, b, R, v * phase)

            block = np.zeros((2, 2), dtype=complex)
            block[a, b] = phase
            declarations = {"onsite_ham_pm": onsite_ham_pm,
                             "Exchange": {1: block}}
            records[name] = _w0_adjudicate_direction(
                fx, terms_of_v, topo, declarations, q_mesh,
                "W0/exchange/{}".format(name))

        for name, rec in records.items():
            self.assertEqual(
                rec["status"], "PASS",
                "W0 granule (a) direction {} status={} (delta_rich={:.3e} "
                "tol={:.3e} max_signal={:.3e} first_failures={}) -- STOP, "
                "do not commit; see TestW0Granules docstring's discrepancy "
                "protocol".format(
                    name, rec["status"], rec["delta_rich"], rec["tol"],
                    rec["max_signal"], rec["failures"][:5]))

        gate = _w0_direction_set_sv_gate(records, "W0/exchange")
        self.assertGreaterEqual(
            gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR,
            "W0/exchange direction-set SVD gate: sv_ratio={:.3e} below "
            "floor {:.3e}".format(gate["sv_ratio"],
                                   ed_oracle_util.SENS_SV_FLOOR))
        self.assertGreaterEqual(
            gate["sigma_min"], 100.0 * gate["delta_rich_max"],
            "W0/exchange direction-set SVD gate: sigma_min={:.3e} below "
            "100*delta_rich_max={:.3e}".format(
                gate["sigma_min"], 100.0 * gate["delta_rich_max"]))

    def test_granule_b_complex_offsite_coulomb_inter(self):
        """(b) L=5, norb=1, off-site CoulombInter, real coupling,
        SINGLE-DIRECTION declaration (R=+1/+2 only -- see the
        investigation note below on why NOT also declaring the -R
        mirror as a separate many-body term is the correct construction
        here, not a simplification). Two directions, g1 (R=+1) and g2
        (R=+2), on the SAME topology for the SVD gate. Exact-magnitude
        adjudication (adjudicate_granule's tol sits at the Richardson
        noise floor, ~1e-4 here) already distinguishes the coefficient
        table's -Re(C) from -2*Re(C) or -Re(C)/2 -- any such factor
        error would show up as a clean, well-outside-tolerance ratio.

        Investigation note (recorded per the discrepancy protocol, since
        an earlier draft of this granule used a Hermitian-mirrored (R
        AND -R) declaration and FAILED by a clean, uniform factor of
        ~2.0 at every grid cell): a density-density operator commutes
        with itself at every site pair, so X_{+R} := sum_j n_j n_{j+R}
        and X_{-R} := sum_j n_j n_{j-R} are the SAME operator (reindex
        j -> j-R). Explicitly declaring BOTH "C at R" and "conj(C) at
        -R" as SEPARATE many-body terms therefore builds
        H = C*X_{+R} + conj(C)*X_{+R} = 2*Re(C)*X_{+R} -- twice the
        physical bond strength a single-channel declaration represents
        -- which is exactly the measured ~2x discrepancy. The topology's
        Hermitian CLOSURE (coeffs[reverse[m]] = conj(coeffs[m])) is
        bookkeeping for `w_expected_from_records`'s representative-based
        LOOKUP (Task 1 Step 1's derivation: one lookup already feeds
        BOTH mirrored target cells), not an instruction to additionally
        re-declare the mirror as a second physical Hamiltonian term on
        the ED side -- exactly mirroring how off-site Ising already
        works below (single R declaration, no explicit mirror) and
        `resolve_interactions` (bond_channels.py) itself: when only ONE
        direction is declared in a real file, the STORED value at that
        channel is the RAW declared value (unaveraged), never doubled.
        A genuinely complex C has ZERO net physical effect here in
        EITHER construction (Im(C) always cancels once Hermiticity is
        enforced for a commuting density-density pair), so this granule
        uses a real coupling and lets adjudicate_granule's exact-value
        comparison catch a magnitude or sign error directly."""
        fx = ed_oracle_util.EDFixture(
            L=5, norb=1, t={(0, 0): 0.7 * np.exp(0.3j)}, eps=(0.0,), T=0.5,
            mu=0.2)
        topo = _TopoLike(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                               [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1, 4, 3]))
        q_mesh = (5, 1, 1)
        onsite_ham_pm = np.zeros((1, 1, 1, 1), dtype=complex)
        C = 0.6 + 0.0j

        records = {}
        for name, R, m_pos in (("g1", 1, 1), ("g2", 2, 3)):
            def terms_of_v(v, R=R):
                return ed_oracle_util.canonical_density_terms(
                    fx, [(0, 0, R, v * C)])

            declarations = {
                "onsite_ham_pm": onsite_ham_pm,
                "CoulombInter": {m_pos: np.array([[C]], dtype=complex)},
            }
            records[name] = _w0_adjudicate_direction(
                fx, terms_of_v, topo, declarations, q_mesh,
                "W0/coulombinter/{}".format(name))

        for name, rec in records.items():
            self.assertEqual(
                rec["status"], "PASS",
                "W0 granule (b) direction {} status={} (delta_rich={:.3e} "
                "tol={:.3e} max_signal={:.3e} first_failures={}) -- STOP, "
                "do not commit; see TestW0Granules docstring's discrepancy "
                "protocol".format(
                    name, rec["status"], rec["delta_rich"], rec["tol"],
                    rec["max_signal"], rec["failures"][:5]))

        gate = _w0_direction_set_sv_gate(records, "W0/coulombinter")
        self.assertGreaterEqual(gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)
        self.assertGreaterEqual(gate["sigma_min"],
                                 100.0 * gate["delta_rich_max"])

    def test_granule_c_offsite_ising(self):
        """(c) L=5, norb=1, off-site Ising, real coupling. Two
        directions, g1 (R=+1) and g2 (R=+2), on the SAME topology for
        the SVD gate (Ising's own Hamiltonian is already Hermitian
        without a mirrored declaration, unlike CoulombInter's genuinely
        complex granule above -- see ``_terms_ising_offsite``'s
        docstring)."""
        fx = ed_oracle_util.EDFixture(
            L=5, norb=1, t={(0, 0): 0.55 * np.exp(-0.25j)}, eps=(0.0,),
            T=0.6, mu=0.15)
        topo = _TopoLike(
            delta_r=np.array([[0, 0, 0], [1, 0, 0], [-1, 0, 0],
                               [2, 0, 0], [-2, 0, 0]]),
            reverse=np.array([0, 2, 1, 4, 3]))
        q_mesh = (5, 1, 1)
        onsite_ham_pm = np.zeros((1, 1, 1, 1), dtype=complex)

        records = {}
        for name, R, m_pos in (("g1", 1, 1), ("g2", 2, 3)):
            def terms_of_v(v, R=R):
                return _terms_ising_offsite(fx, 0, 0, R, v)

            declarations = {
                "onsite_ham_pm": onsite_ham_pm,
                "Ising": {m_pos: np.array([[1.0 + 0j]], dtype=complex)},
            }
            records[name] = _w0_adjudicate_direction(
                fx, terms_of_v, topo, declarations, q_mesh,
                "W0/ising/{}".format(name))

        for name, rec in records.items():
            self.assertEqual(
                rec["status"], "PASS",
                "W0 granule (c) direction {} status={} (delta_rich={:.3e} "
                "tol={:.3e} max_signal={:.3e} first_failures={}) -- STOP, "
                "do not commit; see TestW0Granules docstring's discrepancy "
                "protocol".format(
                    name, rec["status"], rec["delta_rich"], rec["tol"],
                    rec["max_signal"], rec["failures"][:5]))

        gate = _w0_direction_set_sv_gate(records, "W0/ising")
        self.assertGreaterEqual(gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)
        self.assertGreaterEqual(gate["sigma_min"],
                                 100.0 * gate["delta_rich_max"])


# =============================================================================
# Task 6: the ED granule campaign -- the PRODUCTION bond-transverse pipeline
# (RPA._run_transverse_bond_pipeline, Task 5) and the PRODUCTION plain-ladder
# channel (RPA._build_transverse_channel, the pre-existing "ring+ladder"
# path) adjudicated against exact diagonalization to FIRST ORDER in the
# declared coupling. Every granule below drives one of these two PRODUCTION
# entry points directly (anti-vacuity, spec "Phase W -- implementation + ED
# campaign": "a test-invocable internal entry ... executes the bond path
# UNCONDITIONALLY ... every granule and gate W1 call THIS entry" -- Task 6
# extends the same discipline to its own campaign).
#
# COORDINATOR RULING (2026-08-16): this module's own first pass through this
# campaign FAILED the multi-orbital off-site CoulombInter granule (below,
# "the regression pin") at O(1) (0.12 vs tol 5.6e-4) against the DRAFT
# W_pm_bond cross-family formula (an m<->reverse[m] channel swap). The
# investigation (recorded in this module's git history and
# .superpowers/sdd/2026-08-16-bond-transverse-phase-w/task-6-report.md) was
# adjudicated CONFIRMED: the spec's "Cross family" step 2 is AMENDED (see
# docs/superpowers/specs/2026-08-15-bond-transverse-design.md) to DIRECT
# placement (no channel swap, matching the longitudinal bare_bond_vertices
# rule exactly), and BOTH `bond_channels.W_pm_bond` and this module's own
# oracle `w_expected_from_records` (Task 1) were corrected identically. The
# granule below is the PERMANENT regression pin for that fix.
#
# Solver-side construction pattern (uniform across every granule): a bare
# `RPA` stub (`object.__new__` + the minimal attributes
# `_run_transverse_bond_pipeline`/`_build_ham_pm_onsite`/
# `_build_transverse_channel` actually read -- the SAME pattern
# tests/test_bond_transverse_w.py's `TestCollapseRule`/
# `TestDressingDirectFiniteMatrix` already use for this module's sibling
# production tests) driven with `green_kw` built DIRECTLY from the ED
# fixture's own eigenbasis via `_h3_free_two_block_green` (this module's own
# Phase-H helper, reused verbatim) -- never through Wannier90 file I/O,
# which would introduce an independent, unverified Transfer-file sign-
# convention risk this campaign has no need to take on. Since RPA's
# non-interacting Green function depends ONLY on the one-body Hamiltonian
# (hopping + any Zeeman field), never on the declared two-body interaction
# values (`_init_interaction`'s `ham_trans_r`/`ham_extern_r` construction,
# rpa.py -- CoulombIntra/CoulombInter/etc. feed ONLY the vertex, never
# green0), `green_kw` is held EXACTLY FIXED while the declared coupling `v`
# varies -- the same "chi0 fixed, only the vertex varies" construction Gate
# W0 uses (tests/test_bond_transverse_ed.py, Task 1), just now measured
# through the REAL production pipeline (real `resolve_transverse_topology`,
# real `W_pm_bond`, real `transverse_bond_bubble_static`, real `_solve_rpa`)
# instead of the test-local oracle formula. Richardson-differencing the
# EXACT (not first-order-truncated) exact solve's output over `v` at
# `v -> 0` reduces algebraically to the SAME first-order quantity
# `-chi0.dW/dv.chi0` the oracle predicts (chi0 is v-independent by
# construction here), so the two constructions are directly comparable.
#
# Two Matsubara resolutions (nmat=1024, 2048; RAW/no-tail bubble,
# `green0_tail=None`, matching the fwd=down/rev=up bond bubble's own
# O(1/nmat) raw convergence) feed `adjudicate_granule`'s own `delta_nmat`
# convergence estimate -- these values were tuned (empirically, during this
# task's own construction) to bring `delta_nmat` comfortably below
# `max_signal/100` for every granule below while staying fast (a few
# milliseconds of `_solve_rpa` per direction; the ED side, at these small
# fixtures -- L<=5, norb<=2 -- is the dominant, still sub-second, cost).
CAMPAIGN_V1_TASK6 = CAMPAIGN_V1


def _task6_bare_bond_solver(fx, norb):
    """Minimal ``RPA`` stub providing exactly what
    ``_run_transverse_bond_pipeline``/``_build_ham_pm_onsite`` read:
    ``norb``, ``ns``, ``fft_workers``, ``lattice.nvol``/``lattice.shape``,
    and ``ham_info.param_ham``/``param_ham_orig`` (both empty -- every
    granule below declares its off-site content directly in a
    ``TransverseTopology`` via ``bond_channels.resolve_transverse_topology``,
    never through the solver's own interaction files, so
    ``_build_ham_pm_onsite()`` correctly returns an all-zero on-site vertex
    for every granule; the on-site reduction itself is gate W1's job, Task
    5). Mirrors tests/test_bond_transverse_w.py's own bare-solver pattern
    (``TestCollapseRule``/``TestDressingDirectFiniteMatrix``)."""
    solver = object.__new__(rpa_mod.RPA)
    solver.norb = norb
    solver.ns = 2
    solver.fft_workers = 1
    solver.ham_info = types.SimpleNamespace(param_ham={}, param_ham_orig={})
    solver.lattice = types.SimpleNamespace(nvol=fx.L, shape=(fx.L, 1, 1))
    return solver


def _task6_bond_pipeline_derivative(fx, hz, norb, nmat, topo_of_v,
                                     v1=CAMPAIGN_V1_TASK6):
    """Richardson-extrapolated d(chi_pm_bond_static)/dv at v -> 0, from TWO
    calls of the PRODUCTION ``RPA._run_transverse_bond_pipeline`` (at v1 and
    2*v1, plus the v=0 baseline) -- ``green_kw`` fixed for the whole sweep
    (see the module docstring above), only ``topo_of_v(v)`` varies."""
    solver = _task6_bare_bond_solver(fx, norb)
    green_kw = _h3_free_two_block_green(fx, hz, nmat)
    beta = fx.beta

    @functools.lru_cache(maxsize=None)
    def Y(v):
        chi_bond, _collapsed = solver._run_transverse_bond_pipeline(
            green_kw, None, beta, topo_of_v(v))
        return chi_bond

    return ed_oracle_util.richardson(Y, v1)


def _task6_adjudicate_bond_direction(fx, norb, terms_of_v, topo_of_v, label,
                                      nmat1=1024, nmat2=2048,
                                      v1=CAMPAIGN_V1_TASK6, zero_mask=None):
    """One granule direction, end to end, through the REAL production bond
    pipeline: Richardson-in-v of ``_run_transverse_bond_pipeline`` at TWO
    Matsubara resolutions (solver-side convergence, ``delta_nmat``) versus
    Gate W0's own ED-side helper ``_w0_ed_derivative_solver_frame``
    (full-minus-HF-only, frame-mapped via ``transverse_ed_to_solver_map`` --
    reused verbatim, unmodified). ``zero_mask=None`` defaults to all-False
    (a dense/bearing comparison -- the default for every ACTIVE granule
    below; PASS-ZERO controls pass an explicit all-True mask)."""
    D_pred_nmat = _task6_bond_pipeline_derivative(
        fx, 0.0, norb, nmat1, topo_of_v, v1)
    D_pred_2nmat = _task6_bond_pipeline_derivative(
        fx, 0.0, norb, nmat2, topo_of_v, v1)
    topo_probe = topo_of_v(1.0)
    delta_r = np.asarray(topo_probe.delta_r)
    D_ed_v1, D_ed_vhalf = _w0_ed_derivative_solver_frame(
        fx, terms_of_v, delta_r, norb)
    if zero_mask is None:
        zero_mask = np.zeros(D_pred_2nmat.shape, dtype=bool)
    return ed_oracle_util.adjudicate_granule(
        D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask, label)


def _terms_hund_offsite(fx, a, b, R, v):
    """Off-site Hund term, generalizing this module's own on-site Hund
    formula (``_terms_ising_hund``, uhfk.py's documented Hamiltonian
    ``H = -J^Hund (n_up n_up + n_dn n_dn)``, same-spin-only density-density)
    to a nonzero bond displacement ``R``:
    ``H = -v * sum_j (n_{j,a,up} n_{j+R,b,up} + n_{j,a,dn} n_{j+R,b,dn})``.
    Same-spin density operators at DIFFERENT sites commute, so this is
    manifestly Hermitian for real ``v`` without a mirrored declaration
    (mirroring ``_terms_ising_offsite``'s own reasoning)."""
    terms = []
    for j in range(fx.L):
        jr = (j + R) % fx.L
        for s in range(2):
            p = fx.mode(j, a, s)
            r = fx.mode(jr, b, s)
            terms.append((p, p, r, r, -v))
    return terms


def _terms_pairlift_offsite(fx, a, b, R, v):
    """Off-site PairLift term (Sz-BREAKING: ``A = c^+_{a,up} c_{a,dn}
    c^+_{b,up} c_{b,dn}`` changes total Sz by +2 -- SectorED's (Nup, Ndown)
    sector blocking cannot represent it; the DENSE route below is required,
    per the Global Constraints' "Fixture bounds" -- "PairLift's quartic
    operator is Sz-BREAKING -- any PairLift ED reference uses the dense
    route, never SectorED"), generalizing the on-site convention
    (``oracle._terms_for``'s PairLift branch, ``H = v*A + conj(v)*A^dagger``)
    to a nonzero bond displacement ``R`` (orbital ``a`` at site ``j``,
    orbital ``b`` at site ``j+R``). The dagger term mirrors
    ``_terms_exchange_offsite``'s own ``(s, r, q, p, conj(coeff))``
    pattern exactly."""
    terms = []
    for j in range(fx.L):
        jr = (j + R) % fx.L
        p = fx.mode(j, a, 0)
        q = fx.mode(j, a, 1)
        r = fx.mode(jr, b, 0)
        s = fx.mode(jr, b, 1)
        terms.append((p, q, r, s, v))
        terms.append((s, r, q, p, np.conj(v)))
    return terms


def _dense_bond_correlator_transverse(fx, channels, hint=None, h1=None):
    """Dense (non-``SectorED``, no (Nup, Ndown) sector blocking) analogue of
    ``SectorED.bond_correlator_transverse``, for Sz-BREAKING interactions
    (PairLift) SectorED cannot represent. Built from the SAME primitives
    ``ed_oracle_util.chi_connected``/``SectorED.bond_correlator_transverse``
    both use (``ed_oracle_util._diagonalize``/``_static_kernel`` -- full
    diagonalization of the ``(fx.dim, fx.dim)`` many-body Hamiltonian, no
    sector shortcut, tractable only at the small fixtures this granule
    uses), with the SAME ``<A;B^dagger>`` pointwise-kernel contraction
    ``SectorED._lehmann_dagger`` applies (``K[m,n]*A[m,n]*conj(B[m,n])``,
    summed -- NOT a matrix product, unlike ``chi_connected``'s own
    ``A @ B.T`` formula, which pairs a DIFFERENT (q, -q) leg convention this
    method does not use: both legs of the transverse bond operator sit at
    the SAME q here, exactly as ``bond_correlator_transverse`` documents).

    Independently self-validated (recorded, not committed as an automated
    test -- consistent with this campaign's established discipline for a
    new correctness-critical primitive): matches
    ``SectorED.bond_correlator_transverse`` to round-off (~4e-16) both at
    V=0 and for a Sz-CONSERVING interacting reference (where SectorED is
    also valid and the two routes can be cross-checked directly) --
    confirming this dense route's contraction convention BEFORE it is used
    for the Sz-breaking PairLift case SectorED cannot check itself against.

    Parameters
    ----------
    fx : EDFixture
    channels : list of (R, a, b)
        Same contract as ``SectorED.bond_correlator_transverse`` (creation
        leg spin UP at site ``j+R`` orbital ``a``, annihilation leg spin
        DOWN at site ``j`` orbital ``b``) -- ALWAYS the full triple (this
        function has no norb=1 shorthand).
    hint, h1 : see ``ed_oracle_util._diagonalize``.

    Returns
    -------
    ndarray, complex128, shape (fx.L, len(channels), len(channels))
        Connected, SAME 1/fx.L normalization as
        ``SectorED.bond_correlator_transverse``.
    """
    C = fx.annihilators()
    CD = [c.conj().T for c in C]
    ev, w, V = ed_oracle_util._diagonalize(fx, hint=hint, h1=h1)
    K = ed_oracle_util._static_kernel(ev, w, fx.beta)
    n_i = len(channels)
    ops_by_q = []
    for qi in range(fx.L):
        built = []
        for (R, a, b) in channels:
            op = np.zeros((fx.dim, fx.dim), dtype=complex)
            for j in range(fx.L):
                ph = np.exp(2j * np.pi * qi * j / fx.L)
                op += ph * (CD[fx.mode((j + R) % fx.L, a, 0)]
                            @ C[fx.mode(j, b, 1)])
            built.append(V.conj().T @ op @ V)
        ops_by_q.append(built)
    out = np.zeros((fx.L, n_i, n_i), dtype=complex)
    for qi in range(fx.L):
        for i in range(n_i):
            A = ops_by_q[qi][i]
            avgA = (w * np.diag(A)).sum()
            for j in range(n_i):
                B = ops_by_q[qi][j]
                avgB = (w * np.diag(B)).sum()
                val = (K * A * np.conj(B)).sum()
                val -= fx.beta * avgA * np.conj(avgB)
                out[qi, i, j] = val
    return out / fx.L


class TestTask6ProductionPipelineGranules(unittest.TestCase):
    """The ED granule campaign for the PRODUCTION bond-resolved transverse
    pipeline (``RPA._run_transverse_bond_pipeline``, spec "Phase W -- ED
    granules"). Every granule below drives that entry directly (anti-
    vacuity). fx5T (spec, verbatim): L=5 spinful ring, norb=1,
    t=0.7*e^{0.3i}, beta=2, mu=0.2.

    CRITICAL: a FAIL on any of these is a genuine mismatch between the
    (now spec-amended) production pipeline and exact diagonalization --
    STOP, do not commit, report verbatim (the discrepancy protocol, #151:
    re-derive before patching; never bend production or the oracle to force
    a PASS). The multi-orbital off-site CoulombInter granule below is a
    PERMANENT regression pin for the 2026-08-16 cross-family placement fix
    (see the module docstring above) -- it must stay green untouched."""

    @staticmethod
    def _fx5t():
        return ed_oracle_util.EDFixture(
            L=5, norb=1, t={(0, 0): 0.7 * np.exp(0.3j)}, eps=(0.0,), T=0.5,
            mu=0.2)

    def test_granule_coulombinter_g1(self):
        """CoulombInter g1 (R=+-1), real coupling C=0.6 (single-direction
        declaration, avoiding the "declare both R and -R separately"
        doubling trap Gate W0 granule (b) already documents). Measured:
        status=PASS, delta_rich=1.5e-05, delta_nmat=4.4e-05, tol=4.4e-04,
        max_signal=8.1e-02."""
        fx = self._fx5t()
        Cc = 0.6 + 0.0j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 0)): v * Cc}},
                np.eye(3), 1)

        def terms_of_v(v):
            return ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 1, v * Cc)])

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/coulombinter_g1")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_coulombinter_g2(self):
        """CoulombInter g2 (R=+-2), real coupling C=0.6. Measured:
        status=PASS, delta_rich=1.2e-05, delta_nmat=4.4e-05, tol=4.4e-04,
        max_signal=8.1e-02."""
        fx = self._fx5t()
        Cc = 0.6 + 0.0j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((2, 0, 0), (0, 0)): v * Cc}},
                np.eye(3), 1)

        def terms_of_v(v):
            return ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 2, v * Cc)])

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/coulombinter_g2")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_coulombinter_g1_g2_sv_gate(self):
        """The direction-set SVD gate over {g1, g2} (mirroring Gate W0's
        own CoulombInter granule (b) pattern exactly: two shells of the
        SAME type, sharing ONE topology so the sensitivity matrix has a
        genuine common grid to rank -- see the module docstring's note on
        why a naive 4-type combined gate is DEGENERATE below). Measured:
        sv_ratio=0.887, sigma_min=0.204, 100*delta_rich_max=1.5e-03."""
        fx = self._fx5t()
        Cc = 0.6 + 0.0j

        def shared_topo(v_g1=0.0, v_g2=0.0):
            # BOTH shells declared in EVERY call (even when one is exactly
            # 0) so every richardson sample point (v=0, v1, 2*v1) shares
            # the IDENTICAL B=5 channel layout -- richardson differences a
            # B=3 baseline against a B=5 perturbed point otherwise
            # (broadcasting silently produces nonsense, measured: a
            # spurious max_signal~117 during this task's own construction).
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 0)): v_g1,
                                   ((2, 0, 0), (0, 0)): v_g2}},
                np.eye(3), 1)

        def topo_g1(v):
            return shared_topo(v_g1=v * Cc)

        def topo_g2(v):
            return shared_topo(v_g2=v * Cc)

        def terms_g1(v):
            return ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 1, v * Cc)])

        def terms_g2(v):
            return ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 2, v * Cc)])

        records = {
            "coulombinter_g1": _task6_adjudicate_bond_direction(
                fx, 1, terms_g1, topo_g1, "task6/sv/coulombinter_g1"),
            "coulombinter_g2": _task6_adjudicate_bond_direction(
                fx, 1, terms_g2, topo_g2, "task6/sv/coulombinter_g2"),
        }
        for name, rec in records.items():
            self.assertEqual(rec["status"], "PASS", (name, rec))

        gate = _w0_direction_set_sv_gate(records, "task6/coulombinter_g1_g2")
        self.assertGreaterEqual(gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)
        self.assertGreaterEqual(gate["sigma_min"],
                                 100.0 * gate["delta_rich_max"])

    def test_granule_ising_g1(self):
        """Ising g1 (R=+-1), real coupling V=1.0. Measured: status=PASS,
        delta_rich=6.1e-07, delta_nmat=7.3e-05, tol=7.3e-04,
        max_signal=1.4e-01."""
        fx = self._fx5t()
        Vi = 1.0 + 0.0j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"Ising": {((1, 0, 0), (0, 0)): v * Vi}}, np.eye(3), 1)

        def terms_of_v(v):
            return _terms_ising_offsite(fx, 0, 0, 1, v * Vi.real)

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/ising_g1")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_exchange_g1(self):
        """Exchange g1 (R=+-1), a=b=0 (single-orbital fx5T -- the spec's
        own note applies: "off-site Exchange with a == b is VALID input, an
        ordinary non-singleton orbit at R != 0"), complex J=0.5-0.3i.
        Measured: status=PASS, delta_rich=4.1e-07, delta_nmat=8.4e-05,
        tol=8.4e-04, max_signal=1.6e-01."""
        fx = self._fx5t()
        J0 = 0.5 - 0.3j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"Exchange": {((1, 0, 0), (0, 0)): v * J0}}, np.eye(3), 1)

        def terms_of_v(v):
            return _terms_exchange_offsite(fx, 0, 0, 1, v * J0)

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/exchange_g1")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_joint_ray_superposition(self):
        """Joint-ray superposition (spec: "two directions declared
        together, first-order additivity"): CoulombInter g1 AND Ising g1
        declared SIMULTANEOUSLY at R=+-1, BOTH scaled by the SAME v --
        exercises W_pm_bond's per-type accumulation at an OVERLAPPING
        channel/orbital cell (both types insert at the identical (m,0,0)
        diagonal target), not just isolated single-type topologies.
        Measured: status=PASS, delta_rich=2.1e-05, delta_nmat=2.9e-05,
        tol=2.9e-04, max_signal=5.4e-02. Additivity independently confirmed
        (not asserted as its own numeric pin here, since adjudicate_granule
        already IS the additivity check by construction -- the ED side's
        own terms are the literal sum of the two individual granules'
        terms): D_pred(joint) vs D_pred(coulombinter_g1)+D_pred(ising_g1)
        agree to 6.6e-07 against a 5.4e-02 scale (measured separately
        during this task's construction)."""
        fx = self._fx5t()
        Cc = 0.6 + 0.0j
        Vi = 1.0 + 0.0j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 0)): v * Cc},
                 "Ising": {((1, 0, 0), (0, 0)): v * Vi}}, np.eye(3), 1)

        def terms_of_v(v):
            a = ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 1, v * Cc)])
            b = _terms_ising_offsite(fx, 0, 0, 1, v * Vi.real)
            return a + b

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/joint_ray")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_hund_offsite_passzero_control(self):
        """PASS-ZERO off-site Hund control, co-declared with the active
        CoulombInter g1 term (spec, Global Constraints: "PASS-ZERO controls
        ... declared ALONGSIDE one active calibrated CoulombInter term, and
        the granule isolates the control's derivative -- the #151 pattern:
        the zero must be measured as a derivative, not as a fallback").

        INVESTIGATION NOTE (recorded per the discrepancy protocol): a first
        construction co-declared Hund at a FIXED coupling scaled
        INDEPENDENTLY of the varying active CoulombInter term (i.e. the
        active term present at its full value even as Hund's coupling -> 0)
        and FAILED -- ED showed a genuine, non-negligible (~0.09 at
        C0=0.6, scaling roughly linearly with C0, present even at C0=0.01)
        first-order-in-Hund-coupling signal. This is NOT a production bug:
        it is a real beyond-RPA (beyond-ladder) many-body correlation
        between the two DIFFERENT interaction types, visible to EXACT
        diagonalization but structurally outside what ANY ladder/RPA
        resummation (production's OWN approximation class) could ever
        capture -- Hund's own single-vertex insertion into the bubble is
        exactly zero (a Sz-conserving, same-spin-only interaction cannot
        connect the up/down propagators of the transverse loop, the
        on-site table's own argument, unaffected by which site pair is
        involved), and confirmed independently at v->0 with NO active term
        present at all (measured ~3e-6, at the Richardson noise floor).
        The FIX: declare BOTH the active CoulombInter g1 term AND the Hund
        control at the SAME v (this method's construction) -- their
        cross-correlation is then second order in v (proportional to
        v*v), invisible to a first-derivative comparison, while the
        topology/Hamiltonian genuinely carries BOTH declarations
        simultaneously (satisfying anti-vacuity: the ED side's Hamiltonian
        is never Hund-alone-at-zero-coupling). Verified: the joint-v ED
        response matches the CoulombInter-g1-alone response to 4.9e-06
        (against a 0.045-0.08 scale) -- Hund's marginal joint-v
        contribution is zero, confirming the production topology's
        structural omission of Hund (it is not one of
        resolve_transverse_topology's transverse-active types at all) is
        physically correct. Measured (this test): status=PASS,
        delta_rich=1.4e-07, delta_nmat=4.4e-05, tol=4.4e-04,
        max_signal=8.1e-02 (== the CoulombInter-g1-alone signal, since
        production's topology represents CoulombInter only -- Hund is
        never read by resolve_transverse_topology's
        _TRANSVERSE_ACTIVE_TYPES union at all)."""
        fx = self._fx5t()
        Cc = 0.6 + 0.0j

        def topo_of_v(v):
            # Production topology: CoulombInter g1 only -- Hund is
            # structurally absent regardless of v (bit-identical topology
            # at every richardson sample point).
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 0)): v * Cc}},
                np.eye(3), 1)

        def terms_of_v(v):
            active = ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 1, v * Cc)])
            control = _terms_hund_offsite(fx, 0, 0, 1, v)
            return active + control

        rec = _task6_adjudicate_bond_direction(
            fx, 1, terms_of_v, topo_of_v, "task6/hund_control")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_granule_pairlift_offsite_passzero_control_dense(self):
        """PASS-ZERO off-site PairLift control, co-declared with an active
        CoulombInter term at a FIXED (v-independent) coupling -- UNLIKE the
        Hund control above, a FIXED active coupling here does NOT
        contaminate the measurement: PairLift's operator changes total Sz
        by +-2 (Sz-BREAKING, spec/Global Constraints: "PairLift's quartic
        operator is Sz-BREAKING -- any PairLift ED reference uses the dense
        route, never SectorED"), so a SINGLE PairLift insertion can never
        contribute to the Sz-changing-by-1 transverse <S+;S+dagger>
        response AT ANY PERTURBATIVE ORDER (an exact spin-selection-rule
        argument, not merely an RPA-order one -- unlike Hund, which
        conserves Sz and has no such protection). Measured directly:
        status=PASS-ZERO, ED derivative magnitude ~1.3e-09 (Richardson
        noise floor), delta_rich=9.2e-09, tol=9.2e-08, production side
        EXACTLY zero (bit-identical across both Matsubara resolutions,
        since PairLift is not a transverse-active type at all).

        MESH CONSTRAINT (recorded, per the spec's Fixture bounds): the
        bond-transverse topology's mesh-injectivity requirement
        (validate_topology_against_mesh, Task 2) REJECTS any off-site
        channel R with 2*R == 0 (mod L) -- at L=2 this rejects EVERY
        off-site shell (R=1 self-collides, R=-1==R=1 mod 2), so a literal
        2-site ring cannot carry ANY off-site declaration through the
        production pipeline at all (a genuine, load-bearing consequence of
        Task 2's own contract, discovered while constructing this granule
        -- not a defect). This granule therefore uses the smallest ODD ring
        (L=3, avoiding the even-mesh self-reversal) that both satisfies
        mesh-injectivity for R=1 AND keeps the DENSE (full, non-sector)
        diagonalization tractable (dim = 2**(3*1*2) = 64) -- "dense", not
        "2-site", is the load-bearing part of the spec's constraint (Sz-
        breaking implies SectorED is unusable; the dense route's cost is
        what bounds the fixture size, not a literal 2-site requirement)."""
        fx = ed_oracle_util.EDFixture(
            L=3, norb=1, t={(0, 0): 0.55 * np.exp(-0.15j)}, eps=(0.0,),
            T=0.55, mu=0.1)
        C0 = 0.5 + 0.0j
        channels_full = [(0, 0, 0), (1, 0, 0), (2, 0, 0)]

        def topo_of_v(v):
            # v is unused: PairLift is not read by resolve_transverse_
            # topology at all, so the topology (carrying only the FIXED
            # active CoulombInter) is bit-identical at every v.
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 0)): C0}}, np.eye(3), 1)

        def terms_of_v(v):
            active = ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 1, C0)])
            control = _terms_pairlift_offsite(fx, 0, 0, 1, v)
            return active + control

        def dense_derivative(v1):
            @functools.lru_cache(maxsize=None)
            def X_of_v(v):
                terms = terms_of_v(v)
                full = _dense_bond_correlator_transverse(
                    fx, channels_full,
                    hint=ed_oracle_util.h_int_from_terms(fx, terms))
                hf_h1 = ed_oracle_util.hf_h1_from_terms(fx, terms)
                hf_only = _dense_bond_correlator_transverse(
                    fx, channels_full, h1=hf_h1)
                return full - hf_only
            return ed_oracle_util.richardson(X_of_v, v1)

        D_ed_v1 = dense_derivative(CAMPAIGN_V1_TASK6)
        D_ed_vhalf = dense_derivative(CAMPAIGN_V1_TASK6 / 2)

        delta_r = np.asarray(topo_of_v(1.0).delta_r)
        delta_r_1d = _w0_delta_r_1d(delta_r)
        # transverse_ed_to_solver_map wants the norb=1 (R,) shorthand
        # channel list; _w0_channels_for's output is built from the SAME
        # delta_r_1d in the SAME order as channels_full (full (R,a,b)
        # triples) above, so the returned smap indices apply unchanged.
        smap = transverse_ed_to_solver_map(
            _w0_channels_for(delta_r_1d, 1), delta_r_1d, 1)
        D_ed_v1 = D_ed_v1[:, smap][:, :, smap]
        D_ed_vhalf = D_ed_vhalf[:, smap][:, :, smap]

        D_pred_1 = _task6_bond_pipeline_derivative(fx, 0.0, 1, 1024, topo_of_v)
        D_pred_2 = _task6_bond_pipeline_derivative(fx, 0.0, 1, 2048, topo_of_v)
        zero_mask = np.ones(D_pred_2.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_1, D_pred_2, zero_mask,
            "task6/pairlift_control_dense")
        self.assertEqual(rec["status"], "PASS-ZERO", rec)

    def test_granule_multiorbital_offsite_coulombinter_regression(self):
        """THE REGRESSION PIN (spec, "Cross family", 2026-08-16 amendment;
        this module's own coordinator-adjudicated finding): the orbital
        indexing guard -- L=3, norb=2, off-site CoulombInter at R=+1 with
        a=0 != b=1, real coupling C=0.55 (single-direction declaration).
        The ONLY fixture class that can see the m<->reverse[m] cross-family
        placement swap (at a == b the two placements coincide -- this is
        exactly why Gate W0's norb=1 CoulombInter/Ising granules could
        never have caught it).

        HISTORY: this granule's FIRST run (against the DRAFT placement,
        target1=reverse[m] with (a,b) unchanged, target2=m with (b,a)
        swapped) FAILED at O(1): status=FAIL, 220 failing cells,
        max_signal=0.119 vs tol=5.6e-04 (a ~0.12 absolute residual at the
        worst cell). Orientation sweep (recorded in
        .superpowers/sdd/2026-08-16-bond-transverse-phase-w/
        task-6-report.md): the SAME topology declaration matches ED under
        the ED-side relabeling (a=1,b=0,R=+1) or (a=0,b=1,R=-1) -- the
        SAME physical bond, confirming the vertex's target cells were
        placed at the WRONG channel/orbital pair, not that the ED-side
        term builder was wrong. Cross-checked against the ALREADY
        ED-validated longitudinal ``bare_bond_vertices`` Fock-diagonal rule
        (``I=(m,l1,l2)`` gets ``Re(v_bond[m][l1,l2])`` DIRECTLY, no
        ``reverse[]`` swap) -- the "direct placement" hypothesis reproduced
        ED to 5.6e-05 (inside tol), confirming the diagnosis. The
        coordinator CONFIRMED this finding and AMENDED the spec (see the
        module docstring above); ``W_pm_bond`` and ``w_expected_from_
        records`` were corrected identically (DIRECT placement: content at
        the channel carrying the declared coefficient, matching the
        longitudinal rule).

        POST-FIX measurement (this test, run against the corrected
        production code): status=PASS, delta_rich=3.4e-06,
        delta_nmat=5.6e-05, tol=5.6e-04, max_signal=1.19e-01 -- this
        granule is now a PERMANENT regression pin (must stay green
        untouched; a future re-introduction of the channel swap would fail
        it immediately and visibly, at O(1))."""
        fx = ed_oracle_util.EDFixture(
            L=3, norb=2,
            t={(0, 0): 0.5 + 0.2j, (1, 1): 0.35 - 0.15j,
               (0, 1): 0.1 + 0.05j, (1, 0): 0.1 + 0.05j},
            eps=(0.05, -0.03), T=0.45, mu=0.1)
        Cc = 0.55 + 0.0j

        def topo_of_v(v):
            return bond_channels.resolve_transverse_topology(
                {"CoulombInter": {((1, 0, 0), (0, 1)): v * Cc}},
                np.eye(3), 2)

        def terms_of_v(v):
            return ed_oracle_util.canonical_density_terms(
                fx, [(0, 1, 1, v * Cc)])

        rec = _task6_adjudicate_bond_direction(
            fx, 2, terms_of_v, topo_of_v,
            "task6/multiorb_coulombinter_REGRESSION")
        self.assertEqual(
            rec["status"], "PASS",
            "REGRESSION: the multi-orbital off-site CoulombInter orbital-"
            "indexing guard FAILED -- status={} (delta_rich={:.3e} "
            "tol={:.3e} max_signal={:.3e} first_failures={}). This granule "
            "adjudicated the 2026-08-16 cross-family placement fix; a FAIL "
            "here means that fix has regressed -- STOP, do not commit, "
            "report verbatim.".format(
                rec["status"], rec["delta_rich"], rec["tol"],
                rec["max_signal"], rec["failures"][:5]))

    def test_granule_multiorbital_offsite_exchange_both_orientations(self):
        """The multi-orbital off-site Exchange granule (spec, "ED
        granules": "a MULTI-ORBITAL OFF-SITE EXCHANGE granule (norb=2,
        a != b, odd ring L=5 [here L=3, matching the CoulombInter
        regression pin's own fixture], evaluated at a non-self-inverse q)
        that checks BOTH ordered local-channel elements ... plus a
        complex/asymmetric-orientation null-direction variant"). L=3,
        norb=2 (the SAME fixture as the regression pin above), off-site
        Exchange at R=+1 with a=0 != b=1. Two directions on the SAME
        topology: "real" (J=v*1, pins the REAL-coefficient orientation) and
        "imag" (J=v*1i -- the complex/null-direction variant: a purely
        imaginary coupling pins the two-record symmetrization
        independently of the real-coefficient direction, since the two
        ordered records W[(0,a,a),(0,b,b)]/W[(0,b,b),(0,a,a)] carry
        conj(J)/J respectively -- a real-only check cannot distinguish
        conj(J) from J). Since the comparison is over the FULL (nvol, ND,
        ND) grid (not just the on-site sub-block), BOTH ordered elements
        are independently exercised by construction -- the norb=1 fx5T
        Exchange granule above collapses them onto one cell and could not.
        Measured: status=PASS for both directions (real:
        delta_rich=8.9e-07, delta_nmat=9.6e-05, tol=9.6e-04,
        max_signal=1.90e-01; imag: delta_rich=3.1e-06, SAME delta_nmat/
        tol/max_signal). Direction-set SVD gate: sv_ratio=1.0 (to
        displayed precision), sigma_min=0.437, 100*delta_rich_max=3.1e-04
        -- a clean, well-conditioned rank-2 matrix (unlike the
        cross-family {CoulombInter, Ising} pair at a SHARED channel, which
        is exactly degenerate -- see the SV-gate note on
        test_granule_coulombinter_g1_g2_sv_gate)."""
        fx = ed_oracle_util.EDFixture(
            L=3, norb=2,
            t={(0, 0): 0.5 + 0.2j, (1, 1): 0.35 - 0.15j,
               (0, 1): 0.1 + 0.05j, (1, 0): 0.1 + 0.05j},
            eps=(0.05, -0.03), T=0.45, mu=0.1)

        records = {}
        for name, phase in (("real", 1.0 + 0.0j), ("imag", 0.0 + 1.0j)):
            def topo_of_v(v, phase=phase):
                return bond_channels.resolve_transverse_topology(
                    {"Exchange": {((1, 0, 0), (0, 1)): v * phase}},
                    np.eye(3), 2)

            def terms_of_v(v, phase=phase):
                return _terms_exchange_offsite(fx, 0, 1, 1, v * phase)

            records[name] = _task6_adjudicate_bond_direction(
                fx, 2, terms_of_v, topo_of_v,
                "task6/multiorb_exchange/{}".format(name))

        for name, rec in records.items():
            self.assertEqual(rec["status"], "PASS", (name, rec))

        gate = _w0_direction_set_sv_gate(records, "task6/multiorb_exchange")
        self.assertGreaterEqual(gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR)
        self.assertGreaterEqual(gate["sigma_min"],
                                 100.0 * gate["delta_rich_max"])


def _sector_s_minus_correlator(se, channels):
    """The conjugate object to ``SectorED.bond_correlator_transverse``'s
    ``<S+;S+dagger>``: ``<S-;S-dagger>``, built the SAME way but with the
    creation/annihilation SPIN ROLES SWAPPED on the bilinear --
    ``S-_{R,a,b}(q) = sum_j e^{+iqj} c^+_{j+R,a,dn} c_{j,b,up}`` (creation
    DOWN, annihilation UP, the opposite of S+'s creation UP/annihilation
    DOWN). Built directly from ``SectorED``'s own
    ``_build_operator``/``_lehmann_dagger``/``_thermal_avg`` primitives (the
    SAME primitives ``bond_correlator_transverse`` uses internally, already
    called directly from this module's own contract tests, e.g.
    ``TestBondCorrelatorTransverseContract.
    test_mixed_spin_delta_is_plus_minus_one``) -- this is a spin-flipped
    SIBLING of that method's construction, not a re-derivation of its
    machinery, chosen (per the production conjugate-pair granule's brief)
    over the "conjugate/transpose the S+ correlator" route: building S-
    DIRECTLY needs no Kubo-relation derivation and is unambiguous by
    construction (the SAME contraction ``bond_correlator_transverse``
    itself already uses, applied to a differently-labeled bilinear)."""
    fx = se.fx
    seen = set()
    wrapped = []
    for chan in channels:
        R, a, b = se._channel_orbitals(chan)
        Rw = R % fx.L
        key = (Rw, a, b)
        if key in seen:
            raise ValueError(
                "_sector_s_minus_correlator: duplicate channel {!r} after "
                "wrapping R mod L (wrapped key {!r})".format(chan, key))
        seen.add(key)
        wrapped.append((Rw, a, b))
    n_i = len(wrapped)
    out = np.zeros((fx.L, n_i, n_i), dtype=complex)
    for qi in range(fx.L):
        built = []
        for (R, a, b) in wrapped:
            mode_terms = [
                (fx.mode((j + R) % fx.L, a, 1), fx.mode(j, b, 0),
                 np.exp(2j * np.pi * qi * j / fx.L))
                for j in range(fx.L)]
            delta, op = se._build_operator(mode_terms)
            built.append((delta, op, se._thermal_avg(op, delta)))
        for i in range(n_i):
            deltaA, opA, avgA = built[i]
            for j in range(n_i):
                deltaB, opB, avgB = built[j]
                val = se._lehmann_dagger(opA, deltaA, opB, deltaB)
                val -= fx.beta * avgA * np.conj(avgB)
                out[qi, i, j] = val
    return out / fx.L


def _hf_h1_from_terms_at(fx, terms, h1_base):
    """Zeeman-aware generalization of ``ed_oracle_util.hf_h1_from_terms``
    (which hardcodes the density-matrix reference to ``fx.build_h1()``,
    i.e. no Zeeman split): the SAME first-order Hartree-Fock self-energy
    Wick-contraction formula (the ``add()`` closure, reproduced verbatim),
    but built from an explicit ``h1_base`` (the Zeeman-split one-body
    Hamiltonian) rather than the bare ``fx.build_h1()`` -- needed because
    the production conjugate-pair granule's ED reference must subtract the
    HF-only response evaluated on the SAME (Zeeman-split) free reference
    the "full" side uses; the module-level helper cannot express this
    (round-trip verified: at hz=0 this reduces to the SAME numbers
    ``hf_h1_from_terms`` gives, since ``h1_base == fx.build_h1()`` there)."""
    ev, W = np.linalg.eigh(h1_base - fx.mu * np.eye(fx.nmode))
    f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
    rho = ((W * f) @ W.conj().T).T
    S = np.zeros((fx.nmode, fx.nmode), dtype=complex)

    def add(p, q, r, s, c_):
        S[p, q] += c_ * rho[r, s]
        S[r, s] += c_ * rho[p, q]
        S[p, s] += c_ * ((1.0 if q == r else 0.0) - rho[r, q])
        S[r, q] += -c_ * rho[p, s]

    for (p, q, r, s, coeff) in terms:
        add(p, q, r, s, coeff)
    return h1_base + 0.5 * (S + S.conj().T)


class TestTask6ProductionConjugatePairGranule(unittest.TestCase):
    """**The production conjugate-pair granule** (the Phase-H ledger's
    MEDIUM, carried into Task 6, spec "Gate W1 (on-site reduction, exact
    scope)": "the G_up != G_down conjugate-pair identity is adjudicated by
    the production conjugate-pair granule (Task 6)"). On a Zeeman-split
    ON-SITE fixture (G_up != G_down), adjudicates the PRODUCTION plain-
    ladder ``chiq_pm`` -- via ``RPA._build_transverse_channel`` +
    ``RPA._solve_rpa``, the pre-existing "ring+ladder" production path, NOT
    the bond pipeline -- against exact diagonalization, using the
    ``<S-;S-dagger>`` identification (``_sector_s_minus_correlator`` above):
    ``_calc_chi0q_transverse``'s own docstring fixes this unambiguously
    (``chi0_+-[a,c,b,d](r,tau) = -G_up[a,b](r,tau)*G_down[d,c](-r,-tau)``,
    i.e. fwd=UP/rev=DOWN -- the OPPOSITE role assignment from the bond
    pipeline's fwd=DOWN/rev=UP, hence the CONJUGATE object per the spec's
    own "Gate W1" note: "the bond pipeline's chi0 ... computes
    <S+;S+dagger> ... the plain channel's chi0_+- ... computes
    <S-;S-dagger>"). A FAIL here is a PRODUCTION finding (per the task
    brief: "STOP, do not commit, report verbatim").

    Solver-side construction: a bare ``RPA`` stub (the SAME
    ``object.__new__`` pattern as the bond-pipeline granules above) with
    ``green0``/``green0_tail`` built directly from the fixture's own
    Zeeman-split eigenbasis (``_h3_free_two_block_green``, reused verbatim
    -- ALREADY the correct block convention here too: block 0 = +hz = up,
    block 1 = -hz = down, matching ``_calc_chi0q_transverse``'s own "block
    0=up, block1=down" contract), and the on-site vertex ``ham_orig`` built
    directly via ``ring_spin_table("CoulombIntra")`` (the SAME building
    block ``RPA._build_ham_pm_onsite`` uses internally) -- independently
    cross-checked against this module's own H2 table
    (``_h2_ham_pm_expected("CoulombIntra", U, norb=1)``, already
    ED-validated) via ``RPA._assemble_transverse_vertex`` before use
    (exact match, confirming the hand-built ``ham_orig`` tensor is correct
    before it feeds the granule). ``green0`` is v-independent (RPA's
    non-interacting Green function never depends on the declared two-body
    coupling -- see the module docstring above), so the SAME "chi0 fixed,
    only the vertex varies" first-order construction applies here too.

    INVESTIGATION NOTE (recorded per the discrepancy protocol -- two bugs
    were found and fixed during this granule's own construction, neither a
    production defect): (1) an early ED-side term builder reused
    ``tests.test_rpa_vs_ed_oracle._terms_for``, which hardcodes that
    module's OWN ``LX=2`` site count internally -- silently building an
    interaction Hamiltonian covering only 2 of this granule's fixture's 3
    sites. Fixed by switching to ``ed_oracle_util.canonical_density_terms``
    (parameterized on ``fx.L`` correctly); the bare bubble (V=0) comparison
    against production ALREADY matched to ~3e-4 (finite-nmat-only)
    throughout, which is what made the CoulombIntra-only mismatch (a
    genuine ~33% gap at hz=0, where the plain path should reduce to a
    trivially self-consistent value) diagnosable as an ED-side term-count
    bug rather than a solver-side one. (2) the HF-only reference for the
    "full-minus-HF-only" subtraction must be built from the ZEEMAN-SPLIT
    one-body Hamiltonian, not the bare one -- ``ed_oracle_util.
    hf_h1_from_terms`` hardcodes the bare reference, hence
    ``_hf_h1_from_terms_at`` above.

    Measured (this test): status=PASS, delta_rich=3.6e-07,
    delta_nmat=7.2e-05, tol=7.2e-04, max_signal=1.34e-01 -- the plain
    bubble+vertex pair IS self-consistent as the ``<S-;S-dagger>`` conjugate
    object at G_up != G_down, for the first time adjudicated against exact
    diagonalization. This PASS resolves the Phase-H MEDIUM finding."""

    def test_zeeman_split_onsite_conjugate_pair(self):
        fx = ed_oracle_util.EDFixture(
            L=3, norb=1, t={(0, 0): 0.6 * np.exp(0.2j)}, eps=(0.0,), T=0.5,
            mu=0.15)
        hz = 0.3
        h1_zeeman = _h3_zeeman_h1(fx, hz)
        norb = 1
        nd = norb * 2
        channels = [(0,)]

        def solver_for_nmat(nmat):
            solver = object.__new__(rpa_mod.RPA)
            solver.norb = norb
            solver.ns = 2
            solver.fft_workers = 1
            solver.lattice = types.SimpleNamespace(
                nvol=fx.L, shape=(fx.L, 1, 1))
            solver.spin_mode = "spin-diag"
            solver.enable_reduced = False
            solver.T = fx.T
            solver.nmat = nmat
            solver.freq_index = np.arange(nmat)
            solver._chi0q_external = False
            green_kw = _h3_free_two_block_green(fx, hz, nmat)
            solver.green0 = green_kw
            solver.green0_tail = np.zeros_like(green_kw)
            return solver

        def ham_orig_of(v):
            ham_r8 = np.zeros((2, norb) * 4, dtype=complex)
            for spinvec, w in ring_spin_table("CoulombIntra").items():
                s1, s2, s3, s4 = spinvec
                orb = (s4, 0, s3, 0, s1, 0, s2, 0)
                ham_r8[orb] += v * w
            return ham_r8.reshape(nd, nd, nd, nd)

        # Cross-check: the hand-built ham_orig assembles to the SAME
        # ED-validated H2 table entry for CoulombIntra (-U) before use.
        probe_solver = object.__new__(rpa_mod.RPA)
        probe_solver.norb = norb
        probe_solver.ns = 2
        probe_solver.lattice = types.SimpleNamespace(
            nvol=1, shape=(1, 1, 1))
        ham_pm_probe = probe_solver._assemble_transverse_vertex(
            ham_orig_of(1.3)[None])
        expected_probe = _h2_ham_pm_expected(
            "CoulombIntra", 1.3, norb=norb)
        assert_approx_array(ham_pm_probe, expected_probe, rel=0, abs=1e-13)

        def pred_at_nmat(nmat, v1=CAMPAIGN_V1_TASK6):
            solver = solver_for_nmat(nmat)

            @functools.lru_cache(maxsize=None)
            def Y(v):
                nvol = fx.L
                ham_orig = np.broadcast_to(
                    ham_orig_of(v)[None], (nvol, nd, nd, nd, nd)).copy()
                chi0q_pm, ham_pm = solver._build_transverse_channel(
                    None, ham_orig)
                chiq_pm = solver._solve_rpa(chi0q_pm, ham_pm)
                return chiq_pm[nmat // 2]

            return ed_oracle_util.richardson(Y, v1)

        D_pred_1 = pred_at_nmat(1024)
        D_pred_2 = pred_at_nmat(2048)

        def terms_of_v(v):
            return ed_oracle_util.canonical_density_terms(fx, [(0, 0, 0, v)])

        def X_of_v(v):
            terms = terms_of_v(v)
            full = _sector_s_minus_correlator(
                ed_oracle_util.SectorED(fx, terms=terms, h1=h1_zeeman),
                channels)
            hf_h1 = _hf_h1_from_terms_at(fx, terms, h1_zeeman)
            hf_only = _sector_s_minus_correlator(
                ed_oracle_util.SectorED(fx, terms=(), h1=hf_h1), channels)
            return full - hf_only

        v1 = CAMPAIGN_V1_TASK6
        D_ed_v1 = ed_oracle_util.richardson(X_of_v, v1).reshape(
            fx.L, norb, norb, norb, norb)
        D_ed_vhalf = ed_oracle_util.richardson(X_of_v, v1 / 2).reshape(
            fx.L, norb, norb, norb, norb)

        zero_mask = np.zeros(D_pred_2.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_1, D_pred_2, zero_mask,
            "task6/production_conjugate_pair")
        self.assertEqual(
            rec["status"], "PASS",
            "PRODUCTION FINDING: the plain-ladder chiq_pm (via "
            "_build_transverse_channel) is NOT self-consistent as the "
            "<S-;S-dagger> conjugate object at G_up != G_down -- status={} "
            "(delta_rich={:.3e} tol={:.3e} max_signal={:.3e} "
            "first_failures={}). STOP, do not commit, report verbatim.".
            format(rec["status"], rec["delta_rich"], rec["tol"],
                   rec["max_signal"], rec["failures"][:5]))


if __name__ == "__main__":
    unittest.main()
