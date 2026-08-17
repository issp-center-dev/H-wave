#!/usr/bin/env python3

"""``TotalNED`` contract and calibration tests (Phase S, Task 1 of the
bond-transverse campaign, issue #110's validation groundwork).

See ``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``
("Phase S -- the full spinful transverse") and
``docs/superpowers/plans/2026-08-17-bond-transverse-phase-s.md`` for the
binding spec/plan this module works under. ``TotalNED``
(``tests/ed_oracle_util.py``) is the total-N-sector ED oracle Task 1
built to represent Sz-BREAKING one-body Hamiltonians (e.g. a spin-orbit
hop ``t_so * c^dagger_{...,up} c_{...,dn}``), which ``SectorED``'s
(N_up, N_dn) blocking cannot host at all.

This module covers:

- construction: on a Sz-breaking ``h1``, ``SectorED`` must reject it
  (``TestTotalNEDConstructionContract``);
- the ``bond_correlator_transverse`` operator-construction contract
  (shape/dtype, R-wrapping duplicate rejection, the mixed-spin sector
  delta) -- mirroring Phase H's
  ``TestBondCorrelatorTransverseContract`` list, but for the total-N
  delta (always 0, never SectorED's ``(+1, -1)``) --
  (``TestBondCorrelatorTransverseContractTotalN``);
- the calibration gate: on a Sz-CONSERVING fixture (with interactions on
  and a R != 0 channel), ``TotalNED.bond_correlator_transverse`` must
  equal ``SectorED.bond_correlator_transverse`` at round-off -- both are
  exact ED on the identical Hilbert space, so any discrepancy is a real
  bug in one of the two engines (``TestTotalNEDCalibration``);
- the Sz-breaking discrimination check: on a genuinely Sz-breaking
  fixture, the mixed-spin (transverse) correlator blocks TotalNED
  computes are independently, measurably NONZERO -- confirming the
  fixture class Gate S0/S2 will build on actually exercises spin mixing,
  not an accidental near-zero response (``TestTotalNEDSzBreakingDiscrimination``).

Tests must be run from the repository root.
"""

import functools
import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as rpa_mod
from tests import ed_oracle_util
from tests.approx_util import assert_approx_array
from tests.test_bond_transverse_ed import (CAMPAIGN_V1, _hf_h1_from_terms_at,
                                            _dense_bond_correlator_transverse,
                                            _terms_exchange_offsite,
                                            _terms_ising_hund,
                                            _terms_pairlift_offsite,
                                            _w0_direction_set_sv_gate)


def _h1_with_spin_orbit(fx, t_so):
    """``fx.build_h1()`` plus an on-site spin-flip hop
    ``t_so * c^dagger_{j,o,up} c_{j,o,dn} + h.c.``, per site/orbital -- a
    Sz-BREAKING one-body term (spin-orbit-like) that ``SectorED`` cannot
    host (its (N_up, N_dn) sector blocking assumes ``h1`` conserves each
    separately) but ``TotalNED`` can (its blocking is by TOTAL N only;
    see ``TotalNED``'s class docstring for the sector-algebra
    derivation: a bilinear of shape ``c^dagger_p c_q`` always conserves
    total N, whatever the spin parities of ``p`` and ``q`` are)."""
    h1 = fx.build_h1().copy()
    for j in range(fx.L):
        for o in range(fx.norb):
            up = fx.mode(j, o, 0)
            dn = fx.mode(j, o, 1)
            h1[up, dn] += t_so
            h1[dn, up] += np.conj(t_so)
    return h1


class TestTotalNEDConstructionContract(unittest.TestCase):
    """Construction on a small spinful fixture with a Sz-BREAKING ``h1``:
    ``SectorED`` must reject it (``_sector_hamiltonian``'s dropped-element
    floor check fires, since a spin-off-diagonal ``h1`` entry moves
    probability between (N_up, N_dn) sectors); ``TotalNED`` must accept it
    (total N is still conserved for the very same bilinear -- see
    ``TotalNED``'s class docstring)."""

    @staticmethod
    def _fx():
        return ed_oracle_util.EDFixture(
            L=2, norb=1, t={(0, 0): 0.4 - 0.1j}, eps=(0.05,), T=0.5, mu=0.1)

    def test_sectored_rejects_sz_breaking_h1(self):
        fx = self._fx()
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        with self.assertRaises(AssertionError):
            ed_oracle_util.SectorED(fx, h1=h1)

    def test_totalned_accepts_sz_breaking_h1(self):
        fx = self._fx()
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        te = ed_oracle_util.TotalNED(fx, h1=h1)   # must not raise
        # Sanity: the sector bookkeeping actually built something usable
        # (every total-N key from 0 to fx.nmode present, weights sum to 1).
        self.assertEqual(sorted(te._sector_states.keys()),
                          list(range(fx.nmode + 1)))
        total_w = sum(w.sum() for w in te._w.values())
        assert_approx_array(total_w, 1.0, rel=0, abs=1e-12)


class TestBondCorrelatorTransverseContractTotalN(unittest.TestCase):
    """Operator-construction contract for
    ``TotalNED.bond_correlator_transverse``, mirroring the Phase-H list
    (``test_bond_transverse_ed.py::TestBondCorrelatorTransverseContract``)
    against ``SectorED``'s IDENTICAL method: shape/dtype, R-wrapping
    duplicate rejection, and the mixed-spin sector displacement -- which
    is ``0`` here (not SectorED's ``(+1, -1)``; see ``TotalNED``'s class
    docstring's sector algebra)."""

    @staticmethod
    def _small_fx():
        # Same fixture shape as Phase H's contract-test fixture (2-site,
        # two-orbital, spinful) -- fast, disposable, norb=2 so channels
        # use the (R, a, b) 3-tuple form ``_channel_orbitals`` expects.
        t = {(0, 0): 0.3 + 0.1j, (1, 1): 0.2 - 0.05j,
             (0, 1): 0.1j, (1, 0): 0.1j}
        return ed_oracle_util.EDFixture(
            L=2, norb=2, t=t, eps=(0.05, -0.02), T=0.5, mu=0.1)

    def test_output_shape_and_dtype(self):
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        channels = [(0, 0, 0), (1, 0, 0)]
        out = te.bond_correlator_transverse(channels)
        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        self.assertEqual(out.dtype, np.complex128)

    def test_duplicate_channels_after_wrapping_rejected(self):
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        # (1, 0, 0) and (1 - L, 0, 0) both wrap (mod L=2) to R=1 -- same
        # channel after wrapping, must be rejected rather than silently
        # double-counted.
        channels = [(1, 0, 0), (1 - fx.L, 0, 0)]
        with self.assertRaises(ValueError):
            te.bond_correlator_transverse(channels)

    def test_mixed_spin_delta_is_zero(self):
        # Probe the sector displacement via the same operator-
        # construction API bond_correlator_transverse uses internally: a
        # bilinear whose creation leg is spin UP (mode parity 0) and
        # annihilation leg is spin DOWN (mode parity 1) must impart
        # delta = 0 under TOTAL-N blocking -- unlike SectorED, where the
        # SAME bilinear gives delta = (+1, -1) (see the class docstring's
        # sector-algebra derivation: total N = N_up + N_dn is invariant
        # under any single creation-plus-annihilation bilinear).
        fx = self._small_fx()
        te = ed_oracle_util.TotalNED(fx)
        mode_terms = [(fx.mode(j, 0, 0), fx.mode(j, 0, 1), 1.0)
                      for j in range(fx.L)]
        delta, _op = te._build_operator(mode_terms)
        self.assertEqual(delta, 0)


class TestTotalNEDCalibration(unittest.TestCase):
    """The calibration gate (Task 1, Step 4): on a Sz-CONSERVING fixture
    WITH interactions on (a real off-site CoulombInter term, not just the
    free chain) and including a R != 0 channel, ``TotalNED.
    bond_correlator_transverse`` must equal ``SectorED.
    bond_correlator_transverse`` at ``assert_approx_array(rel=0,
    abs=1e-12)`` -- both engines diagonalize the IDENTICAL Hamiltonian on
    the IDENTICAL full Hilbert space (TotalNED's total-N blocks are a
    coarser partition of the SAME Fock space SectorED's (N_up, N_dn)
    blocks refine), so any discrepancy beyond round-off is a real bug
    in one of the two sector-block engines, not a physics difference.

    Fixture: L=3, norb=1 ring (dim=64; TotalNED's largest total-N block
    is N=3, dim=C(6,3)=20 -- small and fast), the SAME (t, eps, T, mu) as
    ``test_bond_transverse_ed.py``'s
    ``test_granule_pairlift_offsite_passzero_control_dense`` fixture (a
    known-good parameter set already exercised in this campaign).
    Interaction: real off-site CoulombInter V=0.4 at R=1 (via
    ``canonical_density_terms`` -- V kept real, since a complex
    coefficient on a Hermitian density-density product would break
    Hermiticity of H). Channels: R=0, 1, 2 (norb=1's ``(R,)`` shorthand),
    so R != 0 is exercised on both legs.
    """

    def test_matches_sectored_with_interactions_and_offsite_channel(self):
        fx = ed_oracle_util.EDFixture(
            L=3, norb=1, t={(0, 0): 0.55 * np.exp(-0.15j)}, eps=(0.0,),
            T=0.55, mu=0.1)
        V = 0.4
        terms = ed_oracle_util.canonical_density_terms(fx, [(0, 0, 1, V)])
        channels = [(0,), (1,), (2,)]

        se = ed_oracle_util.SectorED(fx, terms=terms)
        te = ed_oracle_util.TotalNED(fx, terms=terms)

        want = se.bond_correlator_transverse(channels)
        got = te.bond_correlator_transverse(channels)

        self.assertEqual(got.shape, want.shape)
        assert_approx_array(got, want, rel=0, abs=1e-12)


class TestTotalNEDSzBreakingDiscrimination(unittest.TestCase):
    """The Sz-breaking discrimination check (Task 1, Step 5): on a
    genuinely Sz-breaking fixture (on-site spin-orbit hop ``t_so``, see
    ``_h1_with_spin_orbit``), the mixed-spin transverse correlator blocks
    ``TotalNED.bond_correlator_transverse`` computes must be
    independently, measurably NONZERO -- confirming this fixture class
    genuinely exercises spin mixing (rather than landing, by accident, on
    a near-zero response that would make a later Gate S0/S2 pin vacuous).

    Fixture: L=2, norb=1 (dim=16, tiny/fast), t=0.4-0.1j, eps=0.05,
    T=0.5, mu=0.1, on-site spin-orbit t_so=0.3+0.05j (same fixture and
    t_so as ``TestTotalNEDConstructionContract``, so the two tests probe
    the identical Hamiltonian from two angles: constructibility, then
    physical content).

    Measured lower bound (recorded per the brief's "assert a lower bound
    measured and recorded in the docstring"): the R=0 diagonal entry
    ``out[0, 0, 0]`` (q=0, self channel) has ``abs(...) ~ 0.30`` at these
    parameters -- an O(1) fraction of the O(1) coupling scale, not a
    round-off-level residual. The assertion below uses a conservative
    floor of ``1e-3`` (orders of magnitude above any Richardson/round-off
    noise floor used elsewhere in this campaign, e.g. adjudicate_
    granule's ``1e-13`` epsilon), so it is robust to the exact measured
    value while still being a genuine, non-vacuous discrimination check.
    """

    def test_mixed_spin_blocks_are_nonzero(self):
        fx = ed_oracle_util.EDFixture(
            L=2, norb=1, t={(0, 0): 0.4 - 0.1j}, eps=(0.05,), T=0.5, mu=0.1)
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        te = ed_oracle_util.TotalNED(fx, h1=h1)
        channels = [(0,), (1,)]
        out = te.bond_correlator_transverse(channels)

        self.assertEqual(out.shape, (fx.L, len(channels), len(channels)))
        magnitudes = np.abs(out)
        max_mag = float(magnitudes.max())
        # Not vacuous: at least one cell clears a floor far above any
        # round-off/noise scale.
        self.assertGreaterEqual(max_mag, 1e-3, out)
        # Every q slice individually carries a nonzero mixed-spin signal
        # (not just one lucky cell) -- the "blocks are independently
        # nonzero" requirement.
        for qi in range(fx.L):
            self.assertGreaterEqual(float(magnitudes[qi].max()), 1e-3,
                                     (qi, out[qi]))


# ---------------------------------------------------------------------------
# Mandatory addition carried from the Task-1 review (Phase S, Task 2 brief):
# Task 1's checked-in suite never pins the Sz-breaking CONNECTED-SUBTRACTION
# path of ``TotalNED.bond_correlator_transverse`` -- every Task-1 fixture's
# transverse channel happens to have avgA == avgB == 0 (the calibration
# fixture is Sz-CONSERVING, where SectorED's own delta=(+1,-1) structurally
# forces the average to 0; the Sz-breaking discrimination fixture never has
# an on-site term declared, so Task 1's own construction never exercised a
# case where TotalNED's delta==0 average could be nonzero and enter the
# ``val -= fx.beta * avgA * conj(avgB)`` line non-trivially). The Task-1
# reviewer verified OUT-OF-BAND (not committed) that a dense, non-sector
# Lehmann evaluation matches TotalNED to ~1e-16 on a genuinely Sz-breaking
# fixture where the average is measurably nonzero (L=2, t_so=0.3+0.05j;
# <S+>(q=0, R=0) ~ -0.169+0.028j). This is that check, committed, BEFORE
# Gate S0 (below) relies on TotalNED for its own adjudication.
# ---------------------------------------------------------------------------


class TestTotalNEDDenseFullFockCrossPin(unittest.TestCase):
    """Cross-pins ``TotalNED.bond_correlator_transverse`` against
    ``tests.test_bond_transverse_ed._dense_bond_correlator_transverse`` --
    an INDEPENDENT code path (no (N_up, N_dn) OR total-N sector blocking at
    all: one dense diagonalization of the full ``(fx.dim, fx.dim) = (16,
    16)`` Hamiltonian via ``ed_oracle_util._diagonalize``/
    ``_static_kernel``, the SAME dense machinery
    ``TestTransverseComplexPairHop._chi_pm_ed`` and the module-level
    ``ed_oracle_util.chi_connected`` use elsewhere in this campaign) -- on
    the Sz-breaking fixture (L=2, norb=1, t_so=0.3+0.05j, the SAME fixture
    ``TestTotalNEDConstructionContract``/
    ``TestTotalNEDSzBreakingDiscrimination`` use), covering BOTH channels
    ``[(0,), (1,)]`` (R=0 AND R != 0) so every cell of the (fx.L, 2, 2)
    output -- including the R != 0 diagonal and the R=0/R=1 CROSS terms
    (I != J) -- is compared, at V=0 (h1 only -- isolates the ONE-BODY
    ``<S+>`` nonzero-average path the review flagged, with no interaction
    vertex in play at all) AND with an on-site CoulombIntra interaction
    present (V != 0 -- the same "interactions on" combination Gate S0's
    own ``TotalNED`` calls below actually use).

    ``_dense_bond_correlator_transverse``'s own docstring already records
    an informal self-validation against ``SectorED`` on Sz-CONSERVING
    fixtures; this test is the analogous CHECKED-IN pin for the Sz-
    BREAKING case only ``TotalNED`` (not ``SectorED``) can represent.
    """

    @staticmethod
    def _fx():
        return ed_oracle_util.EDFixture(
            L=2, norb=1, t={(0, 0): 0.4 - 0.1j}, eps=(0.05,), T=0.5, mu=0.1)

    def _check(self, terms):
        fx = self._fx()
        h1 = _h1_with_spin_orbit(fx, 0.3 + 0.05j)
        channels = [(0,), (1,)]

        # bond_correlator_transverse/TotalNED accept the norb=1 (R,)
        # shorthand; _dense_bond_correlator_transverse (built for
        # test_bond_transverse_ed's norb>=1 general channels) always wants
        # the full (R, a, b) triple -- expand the shorthand explicitly so
        # both sides receive the SAME channel list under each function's
        # own contract.
        full_channels = [(R, 0, 0) for (R,) in channels]
        hint = ed_oracle_util.h_int_from_terms(fx, terms) if terms else None
        want = _dense_bond_correlator_transverse(
            fx, full_channels, hint=hint, h1=h1)
        got = ed_oracle_util.TotalNED(
            fx, terms=terms, h1=h1).bond_correlator_transverse(channels)

        self.assertEqual(got.shape, want.shape)
        assert_approx_array(got, want, rel=0, abs=1e-12)
        # non-vacuous: confirm the nonzero-average connected-subtraction
        # path this test exists to pin is actually being exercised (not
        # accidentally landing on a near-zero <S+> again)
        self.assertGreater(float(np.abs(want).max()), 1e-3, want)

    def test_free_h1_only(self):
        self._check(terms=())

    def test_with_onsite_coulomb_intra(self):
        fx = self._fx()
        terms = ed_oracle_util.canonical_density_terms(
            fx, [(0, 0, 0, 0.2)])
        self._check(terms=terms)


# ---------------------------------------------------------------------------
# Gate S0 (Phase S, Task 2): pre-implementation adjudication of the spec's
# dressed-tensor transverse extraction equation, BEFORE any production line
# in rpa.py changes (the ordering that caught both Phase-W spec errors,
# applied here to Phase S's central risk).
#
# Binding spec equation (docs/superpowers/specs/2026-08-15-bond-transverse-
# design.md, "Phase S", restated verbatim in the plan's Global Constraints):
# in the spin-block generalized-index convention g = s*norb + o (spin-
# major -- confirmed empirically below, see ``TestGateS0ExtractionAdjudication``'s
# docstring "Operator-phase / index-convention check"), the physical
# transverse response is the AFTER-DRESSING slice
#     chiq_pm[q, a, c, b, d] = chi_dressed[q, a, norb+c, b, norb+d]
# -- the SAME index selection the pre-S spinful branch of
# ``RPA._build_transverse_channel`` applies to chi0 (the #110 site), applied
# to the DRESSED tensor (the RPA-solved chiq, i.e. ``RPA._solve_rpa``'s
# output on the antisymmetrized general vertex) instead.
# ---------------------------------------------------------------------------


def extract_pm_from_dressed(chi_dressed, norb):
    """The spec equation, verbatim: slice the DRESSED spinful general-solve
    tensor ``chi_dressed[freq, q, g1, g2, g3, g4]`` (``g = s*norb + o``,
    spin-major -- ``g < norb`` is spin-up, ``g >= norb`` is spin-down; the
    SAME generalized-index convention ``chi0q_orig`` and ``chiq`` already
    use throughout ``rpa.py`` for ``spin_mode == "spinful"``) at
    ``[..., a, norb+c, b, norb+d]`` -- i.e. the creation-leg orbital ``a``
    stays on the spin-up block, the annihilation-leg orbital ``c`` is read
    off the spin-down block, and likewise ``(b, d)`` for the second
    (conjugate) leg pair.

    This is byte-for-byte the same slice pattern the PRE-S (buggy, #110)
    ``RPA._build_transverse_channel`` spinful branch applies to ``chi0q_orig``
    (the BARE bubble) -- Phase S's production fix (Task 3) applies it to
    ``chi_dressed`` (the RPA-SOLVED tensor) instead; this is the "slice-
    after-solve vs solve-after-slice" distinction the spec calls out as
    the actual #110 fix. This function does not know or care whether its
    input is bare or dressed -- it is pure index selection -- so it is
    reused unchanged as the Task-3 production reference (cross-pinned bit
    for bit against ``RPA._extract_transverse_from_dressed`` there).

    Parameters
    ----------
    chi_dressed : ndarray, shape (nfreq, nvol, nd, nd, nd, nd), nd = 2*norb
        The dressed (RPA-solved) spinful general tensor, spin-major
        generalized-orbital axes.
    norb : int
        The PHYSICAL orbital count (nd = 2*norb).

    Returns
    -------
    ndarray, shape (nfreq, nvol, norb, norb, norb, norb)
    """
    chi_dressed = np.asarray(chi_dressed)
    if chi_dressed.ndim != 6:
        raise ValueError(
            "extract_pm_from_dressed: chi_dressed must be the 6-dim "
            "(nfreq, nvol, nd, nd, nd, nd) general spinful tensor, got "
            "ndim={}".format(chi_dressed.ndim))
    nd = 2 * norb
    if chi_dressed.shape[2:] != (nd, nd, nd, nd):
        raise ValueError(
            "extract_pm_from_dressed: chi_dressed's trailing orbital axes "
            "{} do not match nd=2*norb={} in every slot".format(
                chi_dressed.shape[2:], nd))
    return chi_dressed[..., 0:norb, norb:2 * norb, 0:norb, norb:2 * norb]


# --- the S0 fixture: L=3 (odd -- q and -q are distinct everywhere except
# q=0, required to pin any q-axis convention; L=2 cannot, every momentum
# there satisfies q=-q, see tests/test_rpa_vs_ed_oracle.py::
# TestMomentumFrame's docstring for the same reasoning applied to the
# longitudinal channel), norb_phys=2 (required to pin orbital LEG ORDER:
# at norb_phys=1 every orbital index is 0, so the spec slice and the
# deliberate g2/g4-swap mutation of Task-2 Step 4 read the SAME array
# element and are numerically indistinguishable -- confirmed by direct
# construction during this gate's own development, discarded), a genuine
# Sz-breaking on-site spin-orbit hop t_so (``_h1_with_spin_orbit``,
# already defined above for Task 1), and BOTH CoulombIntra (the spec's
# named "at least one on-site interaction with transverse content") and a
# complex PairHop (needed for the leg-order pin: CoulombIntra's on-site
# vertex table is symmetric under the pair-slot swap -- see
# tests/test_bond_transverse_ed.py::TestH2OnsiteTableReadjudication's
# choice of PairHop for "the one on-site type whose ham_pm is NOT
# symmetric under the pair-slot swap" -- so a CoulombIntra-only probe
# cannot discriminate a leg-order mistake; confirmed here too: with
# CoulombIntra alone the correct and pair-slot-swapped ED references
# agreed to ~1e-14, discarded as vacuous for this purpose).
_S0_L = 3
_S0_NORB = 2
_S0_T, _S0_MU = 0.5, 0.2
_S0_HOP = {(0, 0): 0.7 + 0.3j, (1, 1): 0.5 - 0.2j,
           (0, 1): 0.25 + 0.15j, (1, 0): 0.25 + 0.15j}
_S0_EPS = (0.10, -0.05)
_S0_LSO = 0.3 + 0.05j
_S0_U_UNIT = 0.13
_S0_PAIRHOP_PHASE = np.exp(0.4j)
_S0_NMAT1, _S0_NMAT2 = 1024, 2048


def _s0_sidx(o, s):
    """1-based, interleaved (2*orb+spin) spin-orbital file index -- the
    Transfer file's OWN convention under enable_spin_orbital (confirmed
    directly from ``Interaction._make_ham_trans``'s comment: "The SO
    transfer file uses the interleaved convention (index = 2*orb + spin,
    matching UHFk and the docs), while RPA works internally in spin-block
    order (index = spin*norb + orb)"), NOT the internal spin-major g the
    solver remaps to and ``extract_pm_from_dressed`` slices."""
    return 2 * o + s + 1


def _s0_make_solver(v, nmat):
    """A real spinful RPA solver fixture (production config-dict pattern,
    adapted from tests/test_bond_transverse_w.py::
    _make_spinful_bond_gate_fixture: enable_spin_orbital, calc_scheme=
    "general", calc_type="ring+ladder" -- but here the transverse BOND
    gate stays OFF/absent, so this exercises the PLAIN ladder path Task 3
    will route through the new extraction). CoulombIntra and PairHop are
    declared with PHYSICAL orbital indices (``Interaction._reshape_
    interaction``: "two-body interactions use physical orbital indices"),
    unlike Transfer's spin-orbital file indices above.

    A SINGLE coupling ``v`` scales BOTH CoulombIntra and PairHop together
    from a genuine v=0 (free) baseline -- required for the "full minus
    HF-only" ED pattern below to isolate a well-defined first-order
    vertex response (an earlier construction that varied PairHop alone
    against a FIXED nonzero CoulombIntra background left a spurious
    ~1% mismatch against production that did not shrink with Nmat or the
    Richardson step -- diagnosed as a finite-difference/pipeline-
    construction defect in the gate's own harness, not a production
    finding, and fixed by scaling both couplings by the same v; see the
    class docstring's "Finite-difference conditioning check").

    Does not call ``solve()``; the caller does, exactly once per (v, nmat).
    """
    d = tempfile.mkdtemp(prefix="s0_gate_")
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n"
                 % (2 * _S0_NORB))
        for _ in range(2 * _S0_NORB):
            f.write("0.0 0.0 0.0\n")

    lines = []
    for (a, b), t in _S0_HOP.items():
        for s in range(2):
            lines.append((1, 0, 0, _s0_sidx(a, s), _s0_sidx(b, s), t))
            lines.append((-1, 0, 0, _s0_sidx(a, s), _s0_sidx(b, s),
                          np.conj(_S0_HOP[(b, a)])))
    for o, e in enumerate(_S0_EPS):
        for s in range(2):
            lines.append((0, 0, 0, _s0_sidx(o, s), _s0_sidx(o, s),
                          complex(e)))
    for o in range(_S0_NORB):
        lines.append((0, 0, 0, _s0_sidx(o, 0), _s0_sidx(o, 1), _S0_LSO))
        lines.append((0, 0, 0, _s0_sidx(o, 1), _s0_sidx(o, 0),
                      np.conj(_S0_LSO)))

    nr = len(set((rx, ry, rz) for (rx, ry, rz, i, j, val) in lines))
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n"
                 % (2 * _S0_NORB, nr, " ".join("1" * nr)))
        for (rx, ry, rz, i, j, val) in lines:
            f.write(" %d %d %d %d %d %.12f %.12f\n"
                     % (rx, ry, rz, i, j, val.real, val.imag))

    inter = {"path_to_input": d, "Geometry": "geom.dat",
              "Transfer": "transfer.dat"}
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("hdr\n%d\n1\n1\n" % _S0_NORB)
        f.write(" 0 0 0 1 1 %.12f 0.0\n" % (v * _S0_U_UNIT))
    inter["CoulombIntra"] = "coulombintra.dat"

    ents = [(0, 1), (1, 0)]
    with open(os.path.join(d, "pairhop.dat"), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n"
                 % (_S0_NORB, len(ents), " ".join("1" * len(ents))))
        for i, (a, b) in enumerate(ents):
            val = v * (_S0_PAIRHOP_PHASE if i == 0
                       else np.conj(_S0_PAIRHOP_PHASE))
            f.write(" 0 0 0 %d %d %.12f %.12f\n"
                     % (a + 1, b + 1, val.real, val.imag))
    inter["PairHop"] = "pairhop.dat"

    param = {"T": _S0_T, "mu": _S0_MU, "CellShape": [_S0_L, 1, 1],
              "SubShape": [1, 1, 1], "Nmat": nmat, "coeff_tail": 1.0}
    info_mode = {"mode": "RPA", "param": param, "enable_spin_orbital": True,
                  "calc_scheme": "general", "calc_type": "ring+ladder"}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
    green_info = io.get_param("green")
    return solver, green_info, d


@functools.lru_cache(maxsize=None)
def _s0_chi_dressed(v, nmat):
    """Step 1: run ``solve()`` once (the plain path -- it materializes
    today's WRONG spinful slice into ``green_info["chiq_pm"]``, which this
    function does not read and which is NOT used as a reference anywhere
    in this gate), then separately retrieve, READ-ONLY, off the solved
    instance: ``chi0q`` (``green_info["chi0q"]``, already the full 6-dim
    general tensor for ``spin_mode == "spinful"`` -- unlike the spin-diag/
    spin-free branches, ``_solve_impl`` never reassigns the local ``chi0q``
    for the spinful branch, so the STORED array is exactly the tensor the
    production longitudinal solve receives) and the antisymmetrized
    vertex the solver assembled (``solver.ham_info.ham_inter_q`` +
    ``solver.ham_info.ham_spinful_exchange``, converted via the PRODUCTION
    ``rpa_mod._to_bubble_pair_convention`` helper -- calling the actual
    production function rather than reimplementing its transpose is the
    "poke internals read-only, never reimplement the assembly" rule the
    task brief states). Then dress via ``RPA._solve_rpa`` IN TEST CODE,
    the same two inputs (``chi0q``, this ``ham``) ``_solve_impl``'s own
    ``sol = self._solve_rpa(chi0q, ham)`` line receives for the spinful
    ring+ladder path (rpa.py, spinful branch of the ``if self.calc_chiq``
    block) -- this reproduces ``green_info["chiq"]`` bit for bit (verified
    informally during construction), but is built independently here per
    the brief's instruction rather than read back off ``green_info``.

    ``spinful_vertex_exchange`` defaults True, so ``ham_spinful_exchange``
    is genuinely present here (CoulombIntra is declared and nonzero for
    any ``v != 0``); at ``v == 0`` no interaction is declared at all and
    ``ham_spinful_exchange`` is ``None`` -- treated as an all-zero tensor
    of ``ham_inter_q``'s trailing (nd, nd, nd, nd) shape, matching "no
    interaction, chi_dressed == chi0q" exactly (the RPA equation's
    ``[1 + chi0*0]^-1 * chi0 = chi0`` identity).

    Cached on ``(v, nmat)``: both the correct and the (Step 4) mutated
    extraction read the SAME dressed tensor, so the solve is not repeated
    per extractor.
    """
    solver, green_info, _d = _s0_make_solver(v, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s0_gate_out_"))
    chi0q = np.asarray(green_info["chi0q"])
    ham_inter_q = solver.ham_info.ham_inter_q
    exch = solver.ham_info.ham_spinful_exchange
    if exch is None:
        exch = np.zeros(ham_inter_q.shape[-4:], dtype=complex)
    ham_long = rpa_mod._to_bubble_pair_convention(ham_inter_q + exch[None, ...])
    return solver._solve_rpa(chi0q, ham_long)


def _s0_pred_at_nmat(nmat, v1, extractor):
    """Richardson-extrapolated d(chiq_pm)/dv at v -> 0, via ``extractor``
    applied to the dressed tensor -- the "test-side dressed-then-extracted
    pipeline" the plan's Task-2 Step 3 names. ``chi0q`` never depends on
    the two-body coupling ``v`` (RPA's bare bubble uses only the one-body
    Hamiltonian), but ``chi_dressed(v) = [1 + chi0*ham(v)]^-1 * chi0`` is
    NOT linear in ``v`` even though ``ham(v)`` is (the RPA series resums
    ``ham(v)``, ``ham(v)^2``, ...) -- so, unlike
    tests/test_bond_transverse_ed.py::TestH2OnsiteTableReadjudication's
    analytic ``-chi0 @ ham_pm @ chi0`` sandwich (exactly linear, no
    Richardson stencil needed), the SOLVER side needs the SAME two-point
    Richardson extrapolation the ED side does."""
    @functools.lru_cache(maxsize=None)
    def Y(v):
        chi_dressed = _s0_chi_dressed(v, nmat)
        chiq_pm = extractor(chi_dressed, _S0_NORB)
        nfreq = chiq_pm.shape[0]
        return chiq_pm[nfreq // 2]        # static (Omega=0) component
    return ed_oracle_util.richardson(Y, v1)


def _s0_extract_pm_mutated(chi_dressed, norb):
    """Task-2 Step 4's deliberate mutation: swap the g2/g4 ROLES in the
    spec slice, reading ``chiq_pm[a, c, b, d]`` from
    ``chi_dressed[..., a, norb+d, b, norb+c]`` instead of the correct
    ``chi_dressed[..., a, norb+c, b, norb+d]`` (i.e. which of the two
    spin-down-block orbitals -- ``c`` or ``d`` -- lands on which of the
    two legs is swapped). Implemented as a transpose of the correctly-
    sliced array (equivalent to re-slicing with c/d exchanged, verified
    by direct index algebra) so this and ``extract_pm_from_dressed`` can
    never silently diverge on anything OTHER than the deliberate swap."""
    chi_dressed = np.asarray(chi_dressed)
    y = chi_dressed[..., 0:norb, norb:2 * norb, 0:norb, norb:2 * norb]
    return np.transpose(y, (0, 1, 2, 5, 4, 3))


@functools.lru_cache(maxsize=1)
def _s0_ed_fixture():
    fx = ed_oracle_util.EDFixture(L=_S0_L, norb=_S0_NORB, t=_S0_HOP,
                                    eps=_S0_EPS, T=_S0_T, mu=_S0_MU)
    h1_so = _h1_with_spin_orbit(fx, _S0_LSO)
    return fx, h1_so


def _s0_pairhop_terms(fx, v, phase):
    """(p, q, r, s, coeff) quartic terms for PairHop at coupling ``v``,
    phase ``phase``, on the declared orbital pair (0, 1) -- the SAME
    operator convention tests/test_rpa_vs_ed_oracle.py::_terms_for's
    "PairHop"/"PairHopComplex" branch uses (Hermitian-closed, the second
    term carrying the conjugate coefficient)."""
    coef = v * phase
    terms = []
    for j in range(fx.L):
        terms.append((fx.mode(j, 0, 0), fx.mode(j, 1, 0),
                      fx.mode(j, 0, 1), fx.mode(j, 1, 1), coef))
        terms.append((fx.mode(j, 1, 0), fx.mode(j, 0, 0),
                      fx.mode(j, 1, 1), fx.mode(j, 0, 1), np.conj(coef)))
    return terms


def _s0_channels():
    return [(0, a, b) for a in range(_S0_NORB) for b in range(_S0_NORB)]


def _s0_X_of_v(v):
    """The ED-side "full minus HF-only" difference (the #151 pattern,
    reused throughout this campaign -- see
    tests/test_bond_transverse_ed.py::_h2_ed_derivative_solver_frame and
    ``TestTask6ProductionConjugatePairGranule``): ``TotalNED`` (not
    ``SectorED`` -- the fixture's ``h1`` carries the Sz-breaking spin-
    orbit term) on the FULL interacting Hamiltonian at coupling ``v``,
    minus ``TotalNED`` on the Hartree-Fock-dressed FREE reference at the
    SAME ``v`` (``_hf_h1_from_terms_at``, the Zeeman/spin-orbit-aware
    generalization of ``ed_oracle_util.hf_h1_from_terms`` that Wick-
    contracts against the DECLARED ``h1_so``, not the bare chain --
    required exactly as tests/test_bond_transverse_ed.py::
    TestTask6ProductionConjugatePairGranule's own INVESTIGATION NOTE (2)
    found for its Zeeman fixture). Returns the RAW (Lehmann-frame)
    ``bond_correlator_transverse`` tensor, channel-reshaped onto
    ``(L, norb, norb, norb, norb)`` with axes ``(q, a, c, b, d)`` where
    channel I <-> (a, c) [the FIRST S+ operator's (creation, annihilation)
    orbitals], channel J <-> (b, d) [the SECOND] -- matching
    ``_s0_channels()``'s ``a``-outer/``b``-inner enumeration order exactly
    (the same reshape convention as tests/test_bond_transverse_ed.py::
    _h2_ed_derivative_solver_frame / ``_h2_channels``).
    """
    fx, h1_so = _s0_ed_fixture()
    channels = _s0_channels()
    terms = (ed_oracle_util.canonical_density_terms(
                fx, [(0, 0, 0, v * _S0_U_UNIT)])
              + _s0_pairhop_terms(fx, v, _S0_PAIRHOP_PHASE))
    full = ed_oracle_util.TotalNED(
        fx, terms=terms, h1=h1_so).bond_correlator_transverse(channels)
    hf_h1 = _hf_h1_from_terms_at(fx, terms, h1_so)
    hf_only = ed_oracle_util.TotalNED(
        fx, terms=(), h1=hf_h1).bond_correlator_transverse(channels)
    return (full - hf_only).reshape(fx.L, _S0_NORB, _S0_NORB, _S0_NORB,
                                     _S0_NORB)


def _s0_frame_map(X):
    """The ED-to-solver frame map for this gate, DERIVED and cross-checked
    empirically during this gate's own construction (recorded here since
    neither half is meaningful without the other, exactly as
    tests/test_bond_transverse_ed.py::transverse_ed_to_solver_map's own
    derivation comment structures its writeup):

    1. Kubo pair-slot swap: ``chiq_pm[q, a, c, b, d] = X[q, b, d, a, c]``
       -- the SAME transpose ``(0, 3, 4, 1, 2)``
       tests/test_rpa_vs_ed_oracle.py::TestTransverseComplexPairHop's
       ``_chi_pm_ed``-derived relation and
       tests/test_bond_transverse_ed.py::transverse_ed_to_solver_map's
       R=0 reduction both establish (static Kubo symmetry relates the
       Lehmann-frame ``<A(q); B(q)>``-type tensor to the solver's own
       slot convention by swapping the two OPERATOR-PAIR slots, not the
       legs within a pair) -- NOT the longitudinal intra-pair map
       ``(0, 2, 1, 4, 3)``, which
       ``TestTransverseComplexPairHop.test_pair_slot_map_is_the_one_kubo_
       symmetry_selects`` pins algebraically as the wrong choice for a
       transverse tensor.

    2. A q-axis reversal, ``q -> (-q) % L``, specific to the
       ``enable_spin_orbital``/``spin_mode == "spinful"`` chi0/chiq
       assembly (NOT needed for the plain "general" scheme --
       tests/test_rpa_vs_ed_oracle.py::TestMomentumFrame's L=3
       calibration pins the NON-reversed ``k+q`` assignment for that
       path, and this gate's own Transfer-file wiring was checked against
       that exact test's file-writing convention before assigning fault
       here). Consistent with a previously-recorded finding for this
       exact code path (an SO/non-SO q-axis Fourier-sign asymmetry noted
       during the RPA-SO ED validation campaign that preceded Phase S);
       this gate is the first to PIN it numerically against a genuinely
       Sz-breaking, L=3 fixture (see the class docstring's "Operator-
       phase / index-convention check" for the ablation that isolated
       it: dropping either half of this map below reproduces a large,
       Nmat/Richardson-CONVERGED, non-noise mismatch -- see there for the
       measured numbers).

    Order does not matter (the two axes -- q, and the (a,c)<->(b,d) pair
    swap -- are disjoint), but is fixed here for a single canonical
    reading.
    """
    L = X.shape[0]
    x_qrev = X[[(-q) % L for q in range(L)], ...]
    return np.transpose(x_qrev, (0, 3, 4, 1, 2))


def _s0_ed_derivative():
    """(D_ed_v1, D_ed_vhalf): the Richardson pair (V1, V1/2) of the frame-
    mapped ED reference, for ``adjudicate_granule``'s ED-side convergence
    measurement (``delta_rich``)."""
    D_v1 = _s0_frame_map(ed_oracle_util.richardson(_s0_X_of_v, CAMPAIGN_V1))
    D_vhalf = _s0_frame_map(
        ed_oracle_util.richardson(_s0_X_of_v, CAMPAIGN_V1 / 2))
    return D_v1, D_vhalf


class TestGateS0ExtractionAdjudication(unittest.TestCase):
    """**Gate S0** (Phase S, Task 2): the pre-implementation adjudication
    of the spec's dressed-tensor transverse extraction equation
    (``extract_pm_from_dressed``) against exact diagonalization
    (``TotalNED``), on a genuinely Sz-BREAKING fixture -- BEFORE any
    change to the production ``rpa.py`` line the spec's Task 3 will add
    (``RPA._extract_transverse_from_dressed``). This is the campaign's
    central-risk gate: the same oracle-first discipline that caught two
    spec errors in Phase W, applied here to the equation's leg-order and
    conjugation.

    Pipeline construction (Step 1): a real spinful ``RPA`` solver (gate
    OFF, plain ``calc_type="ring+ladder"``, ``enable_spin_orbital=True``,
    ``calc_scheme="general"``; L=3, norb_phys=2, complex hopping, a
    genuine on-site Sz-breaking spin-orbit hop ``t_so`` on both orbitals,
    CoulombIntra U and a complex PairHop declared together, a SINGLE
    coupling ``v`` scaling both from a genuine v=0 baseline -- see
    ``_s0_make_solver``'s docstring). ``solve()`` runs once per (v, nmat)
    (its ``chiq_pm`` output uses today's WRONG spinful slice and is never
    read here); ``chi0q`` and the antisymmetrized vertex are retrieved
    READ-ONLY off the solved instance and dressed via ``RPA._solve_rpa``
    IN TEST CODE (``_s0_chi_dressed``) -- the SAME two inputs
    ``_solve_impl``'s own spinful-branch ``_solve_rpa`` call receives.

    The reference extraction (Step 2): ``extract_pm_from_dressed`` (module
    level, above) -- the spec equation verbatim, documented with its
    derivation there.

    Operator-phase / index-convention check (part of the discrepancy
    protocol, performed BEFORE trusting any adjudication result): the
    Transfer-file wiring (interleaved spin-orbital file index -> internal
    spin-major ``g``) was independently confirmed by probing
    ``Interaction.ham_trans_q``'s diagonal against 4 distinct declared
    onsite energies on a throwaway norb_phys=2 fixture (g=0,1 <-> spin-up
    orbitals 0,1; g=2,3 <-> spin-down orbitals 0,1 -- exactly
    ``g = s*norb + o``, matching the spec and
    ``Interaction._make_ham_trans``'s own comment). With that pinned, an
    ablation over the 2x2 (pair-slot-swap x q-reversal) choices on THIS
    gate's own fixture found: identity mapping and q-reversal-only both
    leave a max discrepancy of ~0.086 (scale: max|signal| ~0.111) that
    does NOT shrink between Nmat=1024/2048 nor between the V1/V1-halved
    Richardson estimates (i.e. NOT a convergence artifact); applying BOTH
    the pair-slot swap AND the q-reversal (``_s0_frame_map``) collapses
    the discrepancy to the Nmat/Richardson noise floor (~3e-7, see the
    PASS measurement below) -- a >10^5 tolerance-relative gap between the
    wrong and right frame, i.e. unambiguous. (An earlier, CoulombIntra-
    only construction could not discriminate the pair-slot swap at all --
    CoulombIntra's on-site vertex is itself pair-slot-swap symmetric --
    which is why PairHop was added; recorded above in the fixture's own
    docstring.)

    Finite-difference conditioning check: an earlier construction that
    held CoulombIntra FIXED (nonzero) while varying only PairHop's
    coupling left a ~1% (0.00101 / 0.111) mismatch against production
    that did NOT shrink with Nmat -- traced to a genuine methodological
    defect in the gate's own harness (the "full minus HF-only" pattern
    requires BOTH sides to start from a common v=0 FREE baseline; a fixed
    nonzero background makes ED's HF-only reference miss the background
    coupling's OWN vertex correction, which the solver's chi_dressed(v=0)
    already resums via RPA) -- fixed by scaling both couplings by the
    same v (``_s0_make_solver``), not a production finding.

    **Step 3 measurement (the adjudication)**: status=PASS,
    delta_rich=3.30e-07, delta_nmat=1.13e-07, tol=3.30e-06,
    max_signal=1.11e-01 (>> 100*tol -- a genuine, non-vacuous PASS, not
    INCONCLUSIVE). The spec's extraction equation, applied unmodified to
    the dressed tensor and compared through the frame map above, matches
    ``TotalNED`` to Nmat/Richardson-converged round-off on a fixture
    whose mixed-spin (transverse) content Task 1's own discrimination
    check already proved is genuinely, measurably nonzero (not an
    accidental near-zero response).

    **Step 4 measurement (the deliberate mutation)**: swapping the g2/g4
    roles in the extraction (``_s0_extract_pm_mutated`` -- reading
    ``chi_dressed[..., a, norb+d, b, norb+c]`` instead of the correct
    ``chi_dressed[..., a, norb+c, b, norb+d]``) against the SAME
    (unmutated) ED reference: status=FAIL, magnitude
    max|correct - mutated| = 1.09e-01 (essentially the FULL signal scale,
    98% of max_signal) -- the gate genuinely discriminates the mutation,
    confirming Step 3's PASS is not vacuous.
    """

    def test_step3_dressed_extraction_matches_totalned(self):
        D_ed_v1, D_ed_vhalf = _s0_ed_derivative()
        D_pred_nmat = _s0_pred_at_nmat(
            _S0_NMAT1, CAMPAIGN_V1, extract_pm_from_dressed)
        D_pred_2nmat = _s0_pred_at_nmat(
            _S0_NMAT2, CAMPAIGN_V1, extract_pm_from_dressed)

        zero_mask = np.zeros(D_pred_2nmat.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask,
            "S0/dressed_extraction")
        self.assertEqual(
            rec["status"], "PASS",
            "GATE S0 FAILURE: the spec's dressed-tensor transverse "
            "extraction equation does NOT match exact diagonalization on "
            "the Sz-breaking fixture -- status={} (delta_rich={:.3e} "
            "delta_nmat={:.3e} tol={:.3e} max_signal={:.3e} "
            "first_failures={}). This is the campaign's central risk "
            "(#110 spec adjudication): STOP, do not implement the "
            "production change, report verbatim with the discrepancy-"
            "protocol investigation (operator phases, frame map, "
            "finite-difference conditioning -- each independently).".
            format(rec["status"], rec["delta_rich"], rec["delta_nmat"],
                   rec["tol"], rec["max_signal"], rec["failures"][:5]))

    def test_step4_g2_g4_swap_mutation_must_fail(self):
        D_ed_v1, D_ed_vhalf = _s0_ed_derivative()
        D_pred_nmat = _s0_pred_at_nmat(
            _S0_NMAT1, CAMPAIGN_V1, _s0_extract_pm_mutated)
        D_pred_2nmat = _s0_pred_at_nmat(
            _S0_NMAT2, CAMPAIGN_V1, _s0_extract_pm_mutated)

        # magnitude, recorded per the brief's Step 4 (against the
        # correctly-extracted prediction, same Nmat)
        D_pred_correct_2nmat = _s0_pred_at_nmat(
            _S0_NMAT2, CAMPAIGN_V1, extract_pm_from_dressed)
        magnitude = float(
            np.abs(D_pred_correct_2nmat - D_pred_2nmat).max())
        self.assertGreater(
            magnitude, 1e-2,
            "the deliberate g2/g4-swap mutation must produce a "
            "SIZABLE (non-noise) difference from the correct extraction "
            "for this check to be meaningful; measured magnitude={:.3e}".
            format(magnitude))

        zero_mask = np.zeros(D_pred_2nmat.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_pred_nmat, D_pred_2nmat, zero_mask,
            "S0/mutated_g2_g4_swap")
        self.assertNotEqual(
            rec["status"], "PASS",
            "GATE S0 MUTATION-CHECK FAILURE: the deliberate g2/g4-swap "
            "mutation of the extraction equation PASSED the same "
            "granule the correct extraction passes -- the gate is not "
            "discriminating leg order at all (mutation magnitude was "
            "{:.3e}). Per the task brief, a mutation check failing to "
            "fail is itself a STOP.".format(magnitude))


# ---------------------------------------------------------------------------
# Gate S1 (Phase S, Task 3, Step 3): the Sz-conserving reduction. A spinful
# solver (spin_mode == "spinful") on a system whose one-body Hamiltonian
# carries a genuine but TINY on-site spin-orbit hop `t_so = eps` must
# reproduce, as eps -> 0, the IDENTICAL physical system run through the
# established plain path (`_build_transverse_channel`'s spin-free branch,
# which the SAME fixture reaches automatically at eps == 0.0 EXACTLY --
# every declared hopping amplitude and on-site energy below is spin-
# INDEPENDENT, so with no spin-mixing term at all `RPA._init_wavevec`'s own
# auto-detector (rpa.py, "H is spin-free" / "H is spin-diagnoal" /
# "H is spinful" branch) lands on "spin-free", never "spin-diag", for
# t_so == 0 here specifically).
# ---------------------------------------------------------------------------

_S1_L = 3
_S1_NORB = 2
_S1_T, _S1_MU = 0.5, 0.2
_S1_HOP = {(0, 0): 0.6 + 0.2j, (1, 1): 0.4 - 0.1j,
           (0, 1): 0.2 + 0.05j, (1, 0): 0.2 + 0.05j}
_S1_EPS_ONSITE = (0.10, -0.05)
_S1_U = 0.25
_S1_NMAT = 512


def _s1_sidx(o, s):
    """Same interleaved (2*orb+spin) file-index convention as
    ``_s0_sidx`` -- see its docstring."""
    return 2 * o + s + 1


def _s1_make_solver(t_so, nmat):
    """The Gate S1 fixture: IDENTICAL physical system for every value of
    ``t_so`` except the presence/value of a single declared transfer
    block (the on-site spin-flip hop). At ``t_so == 0.0`` the transfer
    file omits that block entirely -- not merely declares it as zero --
    so the two runs ("reference" and "spinful, eps small") differ by a
    literal file-content difference of one declaration, nothing else
    (same hopping, same on-site energies, same CoulombIntra U, same T,
    mu, Nmat). Every declared hopping amplitude/on-site energy is
    spin-INDEPENDENT (the same value written for both ``s=0`` and
    ``s=1`` file rows), so at ``t_so == 0.0`` the auto-detector lands on
    ``spin_mode == "spin-free"`` (not "spin-diag" -- there is no
    up/down asymmetry at all here to trigger that branch instead); at
    any ``t_so != 0.0`` (however tiny, as long as it clears the
    ``np.allclose`` detection floor of ``atol=1e-8``) it lands on
    ``spin_mode == "spinful"``, routing through the NEW
    ``_extract_transverse_from_dressed`` this task adds.

    Does not call ``solve()``; the caller does.
    """
    d = tempfile.mkdtemp(prefix="s1_gate_")
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n"
                 % (2 * _S1_NORB))
        for _ in range(2 * _S1_NORB):
            f.write("0.0 0.0 0.0\n")

    lines = []
    for (a, b), t in _S1_HOP.items():
        for s in range(2):
            lines.append((1, 0, 0, _s1_sidx(a, s), _s1_sidx(b, s), t))
            lines.append((-1, 0, 0, _s1_sidx(a, s), _s1_sidx(b, s),
                          np.conj(_S1_HOP[(b, a)])))
    for o, e in enumerate(_S1_EPS_ONSITE):
        for s in range(2):
            lines.append((0, 0, 0, _s1_sidx(o, s), _s1_sidx(o, s),
                          complex(e)))
    if t_so != 0.0:
        for o in range(_S1_NORB):
            lines.append((0, 0, 0, _s1_sidx(o, 0), _s1_sidx(o, 1), t_so))
            lines.append((0, 0, 0, _s1_sidx(o, 1), _s1_sidx(o, 0),
                          np.conj(t_so)))

    nr = len(set((rx, ry, rz) for (rx, ry, rz, i, j, val) in lines))
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n"
                 % (2 * _S1_NORB, nr, " ".join("1" * nr)))
        for (rx, ry, rz, i, j, val) in lines:
            f.write(" %d %d %d %d %d %.12f %.12f\n"
                     % (rx, ry, rz, i, j, val.real, val.imag))

    inter = {"path_to_input": d, "Geometry": "geom.dat",
              "Transfer": "transfer.dat"}
    with open(os.path.join(d, "coulombintra.dat"), "w") as f:
        f.write("hdr\n%d\n1\n1\n" % _S1_NORB)
        f.write(" 0 0 0 1 1 %.12f 0.0\n" % _S1_U)
    inter["CoulombIntra"] = "coulombintra.dat"

    param = {"T": _S1_T, "mu": _S1_MU, "CellShape": [_S1_L, 1, 1],
              "SubShape": [1, 1, 1], "Nmat": nmat, "coeff_tail": 1.0}
    info_mode = {"mode": "RPA", "param": param, "enable_spin_orbital": True,
                  "calc_scheme": "general", "calc_type": "ring+ladder"}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
    green_info = io.get_param("green")
    return solver, green_info


@functools.lru_cache(maxsize=None)
def _s1_static_pm(t_so, nmat):
    """Runs the fixture at ``t_so``, returns ``(solver.spin_mode,
    chiq_pm[static])`` -- the static (Omega=0) slice of whichever
    ``chiq_pm`` `solve()` populates, through WHICHEVER path production
    routes to for the resulting ``spin_mode`` (spin-free's
    `_build_transverse_channel` at ``t_so == 0.0``,
    `_extract_transverse_from_dressed` at ``t_so != 0.0``)."""
    solver, green_info = _s1_make_solver(t_so, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s1_gate_out_"))
    chiq_pm = np.asarray(green_info["chiq_pm"])
    nfreq = chiq_pm.shape[0]
    return solver.spin_mode, chiq_pm[nfreq // 2]


class TestGateS1SzConservingReduction(unittest.TestCase):
    """**Gate S1** (Phase S, Task 3, Step 3): as the on-site spin-orbit
    hop ``t_so -> 0``, the spinful plain-transverse ``chiq_pm`` (now
    produced by ``_extract_transverse_from_dressed``) must reduce to the
    IDENTICAL physical system's ``chiq_pm`` computed through the
    established plain path (here: `_build_transverse_channel`'s
    spin-free branch, which the fixture reaches automatically at
    ``t_so == 0.0`` -- see ``_s1_make_solver``'s docstring).

    **Measurement** (recorded here, not assumed -- this is what was
    actually run to pick the tolerance below): at ``t_so = eps``, the
    max-elementwise deviation from the ``t_so = 0`` reference was
    measured over ``eps in {1e-3, 1e-4, 1e-5, 1e-6}`` and found to scale
    as eps**2 (NOT linearly) all the way down to eps=1e-6, with no sign
    of hitting a numerical noise floor before that point:

        eps=1e-3  dev=1.321e-07   (dev/eps**2 = 0.132)
        eps=1e-4  dev=1.321e-09   (dev/eps**2 = 0.132)
        eps=1e-5  dev=1.321e-11   (dev/eps**2 = 0.132)
        eps=1e-6  dev=1.325e-13   (dev/eps**2 = 0.132, matches to 4 s.f.)

    (The quadratic-only scaling -- no measurable linear term -- is
    itself consistent with a symmetry: the LEADING correction from a
    single insertion of chi0's O(t_so) cross-block content into the RPA
    geometric series would enter the transverse-projected output
    linearly in general, so its measured absence here is a genuine
    structural finding about this Hermitian, on-site, orbital-diagonal
    ``t_so`` -- not asserted or relied upon beyond what is measured.)

    At ``eps = 1e-6`` (100x above the ``np.allclose`` atol=1e-8 floor
    the spin_mode auto-detector uses, so still unambiguously
    ``spin_mode == "spinful"``), the RAW (non-extrapolated) deviation
    already sits at 1.325e-13 -- inside the brief's starting tolerance
    (rel=1e-12, abs=1e-13: tolerance = max(1e-12*|ref|, 1e-13) with
    ``max|ref| ~ 0.40``, i.e. ~4.0e-13) with a measured margin of ~3x.
    A 2-point Richardson extrapolation in eps (``2*Y(eps/2) - Y(eps)``,
    removing any residual O(eps) term were one to appear) at
    eps=1e-6/eps=5e-7 tightens this further to 6.589e-14 -- margin
    ~6x -- confirming the RAW measurement was not a coincidental
    near-cancellation. Both the raw and the Richardson-extrapolated
    comparisons are asserted below at the brief's starting tolerance
    (rel=1e-12, abs=1e-13); NEITHER needed loosening (the "STOP and
    investigate before loosening" branch was not triggered)."""

    EPS1 = 1e-6

    def test_spinful_reduces_to_established_plain_path(self):
        mode_ref, ref = _s1_static_pm(0.0, _S1_NMAT)
        self.assertEqual(mode_ref, "spin-free",
                          "the t_so=0 reference fixture must resolve to "
                          "spin_mode=='spin-free' (see _s1_make_solver's "
                          "docstring) -- got {}".format(mode_ref))
        scale = float(np.abs(ref).max())
        self.assertGreater(scale, 1e-2, "non-vacuous reference scale")

        mode1, Y1 = _s1_static_pm(self.EPS1, _S1_NMAT)
        mode2, Y2 = _s1_static_pm(self.EPS1 / 2, _S1_NMAT)
        self.assertEqual(mode1, "spinful")
        self.assertEqual(mode2, "spinful")

        dev_raw = float(np.abs(Y1 - ref).max())
        extrap = 2 * Y2 - Y1
        dev_extrap = float(np.abs(extrap - ref).max())

        msg = ("GATE S1 FAILURE: the spinful chiq_pm does not reduce to "
               "the established spin-free plain path as t_so -> 0 -- "
               "dev_raw={:.3e} dev_extrap={:.3e} at t_so={:.1e} "
               "(scale={:.3e}). This is the #137-adjacent claim (\"the "
               "spin-conserving limit reproduces ring+ladder\") applied "
               "to the NEW extraction path; per the task brief, a "
               "measured gap larger than the starting tolerance is a "
               "STOP, not a cue to loosen it.").format(
                   dev_raw, dev_extrap, self.EPS1, scale)
        assert_approx_array(Y1, ref, rel=1e-12, abs=1e-13, msg=msg)
        assert_approx_array(extrap, ref, rel=1e-12, abs=1e-13, msg=msg)


# ---------------------------------------------------------------------------
# Gate S2 (Phase S, Task 3, Step 4): the negative control -- #110 made
# visible. Reproduces the OLD (pre-Phase-S) slice-THEN-solve computation IN
# TEST CODE, on the SAME Sz-breaking S0 fixture Gate S0 adjudicated, and
# shows it FAILS the identical granule the (post-change) production path
# PASSES. Also cross-pins production `chiq_pm` against
# `extract_pm_from_dressed` applied to the test-side dressed tensor, bit for
# bit, on the shared fixture.
# ---------------------------------------------------------------------------


@functools.lru_cache(maxsize=None)
def _s2_production_chiq_pm(v, nmat):
    """Runs the S0 fixture through the PRODUCTION ``solve()`` entry
    point (post-Phase-S) and returns ``green_info["chiq_pm"])`` -- the
    array `_extract_transverse_from_dressed` now populates for this
    (spinful) fixture."""
    solver, green_info, _d = _s0_make_solver(v, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s2_gate_out_"))
    return np.asarray(green_info["chiq_pm"])


def _s2_new_Y(v, nmat):
    arr = _s2_production_chiq_pm(v, nmat)
    return arr[arr.shape[0] // 2]


def _s2_new_pred(nmat, v1):
    return ed_oracle_util.richardson(lambda v: _s2_new_Y(v, nmat), v1)


@functools.lru_cache(maxsize=None)
def _s2_old_extract(v, nmat):
    """Reproduces the OLD (pre-Phase-S, issue #110) slice-THEN-solve
    computation, IN TEST CODE, on the S0 fixture: slice ``chi0`` (the
    BARE bubble, ``green_info["chi0q"]``) at the Sz-conserving block
    FIRST -- the exact index selection the deleted spinful branch of
    ``_build_transverse_channel`` applied to ``chi0q_orig`` -- and THEN
    solve RPA (``RPA._solve_rpa``) on that REDUCED slice, dressed by the
    transverse vertex ``_assemble_transverse_vertex`` extracts from
    ``ham_inter_q`` ALONE (NOT ``ham_inter_q + ham_spinful_exchange`` --
    the deleted branch's own ``ham_orig`` argument, at its
    ``_solve_impl`` call site, was always the bare ``ham_inter_q``; the
    antisymmetrized ``ham_spinful_exchange`` combination is exclusive to
    the LONGITUDINAL solve's own ``ham_long``, never passed to
    ``_build_transverse_channel``)."""
    solver, green_info, _d = _s0_make_solver(v, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s2_gate_out_"))
    chi0q_orig = np.asarray(green_info["chi0q"])
    norb = _S0_NORB
    chi0q_pm_sliced = chi0q_orig[:, :,
                                 0:norb, norb:2 * norb,
                                 0:norb, norb:2 * norb].copy()
    ham_pm = solver._assemble_transverse_vertex(solver.ham_info.ham_inter_q)
    sol_pm_old = solver._solve_rpa(chi0q_pm_sliced, ham_pm)
    return np.asarray(rpa_mod._bk.to_host(sol_pm_old))


def _s2_old_Y(v, nmat):
    arr = _s2_old_extract(v, nmat)
    return arr[arr.shape[0] // 2]


def _s2_old_pred(nmat, v1):
    return ed_oracle_util.richardson(lambda v: _s2_old_Y(v, nmat), v1)


class TestGateS2ProductionRouting(unittest.TestCase):
    """**Gate S2** (Phase S, Task 3, Step 4): the negative control that
    makes issue #110 directly visible in the same granule apparatus
    Gate S0 used. On the IDENTICAL Sz-breaking S0 fixture
    (L=3, norb=2, genuine on-site spin-orbit ``t_so``, CoulombIntra +
    complex PairHop):

    - the OLD (pre-Phase-S) slice-THEN-solve computation
      (``_s2_old_pred``, reproduced here in test code since the
      production branch that used to compute it was deleted in this
      same task) must FAIL the Task-2 granule against ``TotalNED``
      (``_s0_ed_derivative``) -- this is issue #110 itself, made
      directly measurable instead of merely asserted;
    - the NEW production path (``_s2_new_pred``, i.e. the actual
      ``solve()`` output post this task's routing change) must PASS the
      IDENTICAL granule;
    - production ``chiq_pm`` must equal ``extract_pm_from_dressed``
      applied to the test-side dressed tensor (``_s0_chi_dressed``, the
      SAME pipeline Gate S0 built), BIT FOR BIT, on the shared fixture
      (confirming the production routing change did not introduce any
      numerical divergence from the pre-adjudicated reference path, not
      even at the last bit).

    **Measured tuples** (recorded here; both runs share the identical
    ED reference ``D_ed_v1``/``D_ed_vhalf``, so the two ``status``
    fields are directly comparable):

        OLD (slice-then-solve):  status=FAIL, delta_rich=3.300e-07,
            delta_nmat=1.132e-07, tol=3.300e-06, max_signal=1.110e-01,
            n_failures=44
        NEW (production, dressed extraction): status=PASS,
            delta_rich=3.300e-07, delta_nmat=1.134e-07, tol=3.300e-06,
            max_signal=1.111e-01, n_failures=0

    The NEW tuple matches Gate S0's own Step-3 measurement to the
    reported precision (Task 2's ``TestGateS0ExtractionAdjudication``:
    delta_rich=3.30e-07, delta_nmat=1.13e-07, tol=3.30e-06,
    max_signal=1.11e-01) -- expected, since the cross-pin below proves
    production ``chiq_pm`` and Gate S0's own
    ``extract_pm_from_dressed(chi_dressed)`` are the SAME array, bit for
    bit, on this fixture: the production routing change and Gate S0's
    pre-adjudicated pipeline compute the identical quantity.

    Cross-pin (production `chiq_pm` vs `extract_pm_from_dressed`
    applied to `chi_dressed`): ``np.array_equal`` True, max|diff| = 0.0
    -- exact bit-for-bit agreement, not merely close, at
    ``v=1.0, nmat=1024``."""

    def test_old_slice_then_solve_fails_the_s0_granule(self):
        D_ed_v1, D_ed_vhalf = _s0_ed_derivative()
        D_old_nmat = _s2_old_pred(_S0_NMAT1, CAMPAIGN_V1)
        D_old_2nmat = _s2_old_pred(_S0_NMAT2, CAMPAIGN_V1)
        zero_mask = np.zeros(D_old_2nmat.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_old_nmat, D_old_2nmat, zero_mask,
            "S2/old_slice_then_solve")
        self.assertNotEqual(
            rec["status"], "PASS",
            "GATE S2 (negative control) FAILURE: the OLD (pre-Phase-S, "
            "issue #110) slice-then-solve computation PASSED the S0 "
            "granule -- the negative control is not discriminating "
            "issue #110 at all. status={} delta_rich={:.3e} "
            "delta_nmat={:.3e} tol={:.3e} max_signal={:.3e}".format(
                rec["status"], rec["delta_rich"], rec["delta_nmat"],
                rec["tol"], rec["max_signal"]))

    def test_production_passes_the_s0_granule(self):
        D_ed_v1, D_ed_vhalf = _s0_ed_derivative()
        D_new_nmat = _s2_new_pred(_S0_NMAT1, CAMPAIGN_V1)
        D_new_2nmat = _s2_new_pred(_S0_NMAT2, CAMPAIGN_V1)
        zero_mask = np.zeros(D_new_2nmat.shape, dtype=bool)
        rec = ed_oracle_util.adjudicate_granule(
            D_ed_v1, D_ed_vhalf, D_new_nmat, D_new_2nmat, zero_mask,
            "S2/production_dressed_extraction")
        self.assertEqual(
            rec["status"], "PASS",
            "GATE S2 FAILURE: the production (post-Phase-S) chiq_pm does "
            "not pass the SAME granule the OLD slice-then-solve "
            "computation fails -- the #110 fix did not land. status={} "
            "delta_rich={:.3e} delta_nmat={:.3e} tol={:.3e} "
            "max_signal={:.3e} first_failures={}".format(
                rec["status"], rec["delta_rich"], rec["delta_nmat"],
                rec["tol"], rec["max_signal"], rec["failures"][:5]))

    def test_production_matches_reference_extraction_bit_for_bit(self):
        v, nmat = 1.0, _S0_NMAT1
        chiq_pm_prod = _s2_production_chiq_pm(v, nmat)
        chi_dressed = _s0_chi_dressed(v, nmat)
        chiq_pm_ref = extract_pm_from_dressed(chi_dressed, _S0_NORB)
        self.assertEqual(chiq_pm_prod.shape, chiq_pm_ref.shape)
        self.assertTrue(
            np.array_equal(chiq_pm_prod, chiq_pm_ref),
            "production chiq_pm and extract_pm_from_dressed(chi_dressed) "
            "must be BIT FOR BIT identical on the shared S0 fixture -- "
            "max|diff|={:.3e}".format(
                float(np.abs(chiq_pm_prod - chiq_pm_ref).max())))


# ---------------------------------------------------------------------------
# Task 4 (Phase S): the granule campaign -- first-order adjudication of
# PRODUCTION solve()'s chiq_pm (the anti-vacuity rule -- every granule below
# drives the public solve() entry point exactly as Gate S2's production leg
# does, never a hand-assembled chi_dressed) against TotalNED, on the t_so
# fixture FAMILY (a norb=1 variant, reusing Task 1's own Sz-breaking
# construction contract fixture verbatim, and a norb=2 variant, reusing
# Gate S0's own (L, HOP, EPS, LSO, T, MU) verbatim), for EACH on-site
# interaction type with transverse content (spec: "The verified on-site
# algebra", H2's own 7-type table), plus PASS-ZERO controls where a type's
# transverse content genuinely vanishes on this fixture class -- DERIVED
# empirically here, never assumed (PairLift in particular: Sz-breaking, so
# H2's Sz-conserving PASS-ZERO finding is NOT automatically valid once t_so
# lifts the selection rule that produced it -- see
# TestTask4GranuleCampaign's class docstring for the measured verdict).
# ---------------------------------------------------------------------------


def _s4_write_onsite_interaction(d, fname, norb, entries):
    """Write a generic on-site (R=0) two-body interaction declaration file
    in the k-space Wannier90-like format ``QLMSkInput`` reads -- the SAME
    layout ``_s0_make_solver``/``_s1_make_solver`` already write for
    CoulombIntra/PairHop (header line, norb, "1", one nr-flag, then
    "Rx Ry Rz o1 o2 val_re val_im" rows), generalized to an arbitrary list
    of ``(o1_1based, o2_1based, complex_value)`` entries so ONE writer
    serves every on-site interaction TYPE -- confirmed this is the shared
    format across types by direct inspection of
    ``tests/rpa/input_2orb/{onsite_inter.dat,onsite_cplx_hermitian.dat,
    coulombintra.dat}`` (CoulombInter/Exchange/CoulombIntra all use it) and
    ``tests/test_rpa_ladder.py``'s own reuse of one generic
    ``coulombinter.dat`` file under both the ``'Hund'`` and ``'PairLift'``
    interaction-dict keys."""
    with open(os.path.join(d, fname), "w") as f:
        f.write("hdr\n%d\n1\n1\n" % norb)
        for (o1, o2, val) in entries:
            f.write(" 0 0 0 %d %d %.12f %.12f\n"
                     % (o1, o2, val.real, val.imag))


def _s4_make_solver(L, norb, hop, eps, lso, T, mu, kind, v, nmat, phase=1.0,
                     calc_type="ring+ladder", spinful_vertex_exchange=True):
    """Build a spinful RPA solver fixture (``enable_spin_orbital=True``,
    ``calc_scheme="general"``, ``calc_type`` defaulting to
    ``"ring+ladder"`` (override via the ``calc_type`` keyword), transverse
    bond gate absent/OFF -- the plain ladder path) declaring EXACTLY ONE
    on-site interaction type (``kind``) at coupling ``v`` (times ``phase``
    for PairHop, the one complex-capable type here), on a FIXED on-site
    spin-orbit hop ``lso`` present on every orbital (the t_so fixture
    family's own defining feature) plus the given complex hopping/on-site
    energies -- generalizing ``_s0_make_solver``/``_s1_make_solver``'s
    construction pattern to an arbitrary ``(L, norb)`` and a SINGLE
    declared type (rather than S0's fixed CoulombIntra+PairHop pair).

    The interaction file is ALWAYS written (even at ``v == 0.0``, with
    every entry valued exactly ``0.0``) -- matching ``_s0_make_solver``'s
    own pattern, and REQUIRED: ``RPA.calc_chiq`` is
    ``self.ham_info.has_interaction()``, a presence check on the declared
    interaction KEYS, not a check on their VALUES, so omitting the file
    entirely at v=0 leaves ``calc_chiq`` False and ``solve()`` never
    populates ``chiq``/``chiq_pm`` at all (measured directly while
    constructing this fixture: ``green_info`` held only
    ``{"chi0q", "chi0q_freq_meta"}`` with the file omitted). Since the
    ASSEMBLED vertex is all-zero regardless of which interaction TYPE key
    happens to hold the zero-valued declaration, the v=0 solve is still
    fully INDEPENDENT of ``kind`` -- ``_s4_free_chiq_pm`` below exploits
    this to compute it once and share it across every type (the
    runtime-budget optimization this module's docstring's "keep the
    module under ~10 min total" constraint needs).

    Does not call ``solve()``; the caller does."""
    d = tempfile.mkdtemp(prefix="s4_gate_")
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % (2 * norb))
        for _ in range(2 * norb):
            f.write("0.0 0.0 0.0\n")

    def sidx(o, s):
        return 2 * o + s + 1

    lines = []
    for (a, b), t in hop.items():
        for s in range(2):
            lines.append((1, 0, 0, sidx(a, s), sidx(b, s), t))
            lines.append((-1, 0, 0, sidx(a, s), sidx(b, s),
                          np.conj(hop[(b, a)])))
    for o, e in enumerate(eps):
        for s in range(2):
            lines.append((0, 0, 0, sidx(o, s), sidx(o, s), complex(e)))
    for o in range(norb):
        lines.append((0, 0, 0, sidx(o, 0), sidx(o, 1), lso))
        lines.append((0, 0, 0, sidx(o, 1), sidx(o, 0), np.conj(lso)))

    nr = len(set((rx, ry, rz) for (rx, ry, rz, i, j, val) in lines))
    with open(os.path.join(d, "transfer.dat"), "w") as f:
        f.write("hdr\n%d\n%d\n%s\n" % (2 * norb, nr, " ".join("1" * nr)))
        for (rx, ry, rz, i, j, val) in lines:
            f.write(" %d %d %d %d %d %.12f %.12f\n"
                     % (rx, ry, rz, i, j, val.real, val.imag))

    inter = {"path_to_input": d, "Geometry": "geom.dat",
              "Transfer": "transfer.dat"}

    if kind == "CoulombIntra":
        entries = [(1, 1, complex(v))]
    elif kind in ("CoulombInter", "Exchange", "Hund", "PairLift"):
        entries = [(1, 2, complex(v)), (2, 1, complex(v))]
    elif kind == "PairHop":
        cv = v * phase
        entries = [(1, 2, cv), (2, 1, np.conj(cv))]
    else:
        raise ValueError(kind)
    _s4_write_onsite_interaction(d, "onsite.dat", norb, entries)
    inter[kind] = "onsite.dat"

    param = {"T": T, "mu": mu, "CellShape": [L, 1, 1], "SubShape": [1, 1, 1],
              "Nmat": nmat, "coeff_tail": 1.0,
              "spinful_vertex_exchange": spinful_vertex_exchange}
    info_mode = {"mode": "RPA", "param": param, "enable_spin_orbital": True,
                  "calc_scheme": "general", "calc_type": calc_type}
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    solver = rpa_mod.RPA(io.get_param("ham"), {}, info_mode)
    green_info = io.get_param("green")
    return solver, green_info


def _s4_terms_of_v(fx, kind, v, phase=1.0):
    """(p, q, r, s, coeff) quartic terms for ONE on-site type at coupling
    ``v``, matching ``_s4_make_solver``'s file-declared convention exactly:

    - CoulombIntra/CoulombInter: ``canonical_density_terms`` (the SAME
      builder Task 1's own calibration gate uses), on orbital 0 alone
      (CoulombIntra) or the (0, 1) pair (CoulombInter, all four spin
      pairings, matching ``oracle._terms_for``'s on-site CoulombInter
      branch and the shared ``onsite_inter.dat``-style file).
    - Exchange: ``_terms_exchange_offsite(fx, 0, 1, R=0, v)`` (real ``v``)
      -- cross-checked algebraically against ``oracle._terms_for``'s
      on-site Exchange branch: for REAL ``v`` the two are the IDENTICAL
      operator (the two orderings of the (a,b) declaration commute
      because they act on disjoint mode sets, so ``conj(v) == v``
      collapses the only difference between the two builders).
    - Hund: ``_terms_ising_hund`` (H2's own on-site Hund builder, already
      ED-validated at the Hamiltonian level by
      ``test_ising_hund_term_builders_match_dense_hamiltonian``).
    - PairLift: ``_terms_pairlift_offsite(fx, 0, 1, R=0, v)`` (real ``v``,
      same disjoint-mode-set argument as Exchange above collapses its
      ``conj(v)`` partner term to ``oracle._terms_for``'s on-site PairLift
      convention).
    - PairHop: ``_s0_pairhop_terms`` (already defined above for Gate S0,
      genuinely complex-capable via ``phase``).
    """
    if kind == "CoulombIntra":
        return ed_oracle_util.canonical_density_terms(fx, [(0, 0, 0, v)])
    if kind == "CoulombInter":
        return ed_oracle_util.canonical_density_terms(fx, [(0, 1, 0, v)])
    if kind == "Exchange":
        return _terms_exchange_offsite(fx, 0, 1, 0, v)
    if kind == "Hund":
        return _terms_ising_hund(fx, "Hund", v)
    if kind == "PairLift":
        return _terms_pairlift_offsite(fx, 0, 1, 0, v)
    if kind == "PairHop":
        return _s0_pairhop_terms(fx, v, phase)
    raise ValueError(kind)


# --- norb=1 variant: reuses Task 1's OWN Sz-breaking construction-contract
# fixture verbatim (L=2, norb=1, t=0.4-0.1j, eps=0.05, T=0.5, mu=0.1,
# t_so=0.3+0.05j -- already proven to carry a measurably nonzero mixed-spin
# signal by TestTotalNEDSzBreakingDiscrimination). The only on-site type
# that needs just ONE orbital is CoulombIntra (the rest all name an orbital
# PAIR (0, 1) by construction).
_S4_1_L, _S4_1_NORB = 2, 1
_S4_1_T, _S4_1_MU = 0.5, 0.1
_S4_1_HOP = {(0, 0): 0.4 - 0.1j}
_S4_1_EPS = (0.05,)
_S4_1_LSO = 0.3 + 0.05j

# --- norb=2 variant: reuses Gate S0's OWN (L, HOP, EPS, LSO, T, MU)
# verbatim (already proven, by Gate S0's own PASS, to carry a genuinely
# nonzero, correctly-extracted mixed-spin transverse signal under
# production's post-Phase-S routing) for CoulombInter (inter-orbital),
# Exchange, PairHop (complex, reusing Gate S0's own PairHop phase), and the
# Hund/PairLift controls.
_S4_2_L, _S4_2_NORB = _S0_L, _S0_NORB
_S4_2_T, _S4_2_MU = _S0_T, _S0_MU
_S4_2_HOP, _S4_2_EPS, _S4_2_LSO = _S0_HOP, _S0_EPS, _S0_LSO
_S4_2_PAIRHOP_PHASE = _S0_PAIRHOP_PHASE

# A deliberate runtime-budget reduction from Gate S0's Nmat=(1024, 2048):
# the granule campaign runs SIX independent type-isolated fixtures (Gate S0
# ran ONE), so each granule's own Nmat pair is lowered to (512, 1024).
# Measured margins (recorded per-granule in TestTask4GranuleCampaign's
# docstring) confirm this still clears the 100x-signal-over-tol power floor
# by a wide margin for every PASS granule -- the reduction costs precision
# headroom, not adjudication validity.
_S4_NMAT1, _S4_NMAT2 = 512, 1024

# The coupling SCALE multiplying the Richardson finite-difference variable
# ``v`` (which itself ranges over {0, CAMPAIGN_V1/2, CAMPAIGN_V1,
# 2*CAMPAIGN_V1}, the campaign-wide small-step convention) -- an O(1)
# prefactor so the declared file/ED coupling sits at a physically
# reasonable scale (matching Gate S0's own ``_S0_U_UNIT``-style
# construction), kept UNIFORM across every type in this campaign for a
# single, auditable scale.
_S4_UNIT = 0.3


def _s4_scaled(v):
    return v * _S4_UNIT


@functools.lru_cache(maxsize=None)
def _s4_free_chiq_pm(variant, nmat):
    """The v=0 (no interaction declared AT ALL) production ``chiq_pm``,
    SHARED across every on-site TYPE granule on the given fixture variant
    (``"norb1"``/``"norb2"``) -- the free background depends only on
    ``(L, hop, eps, t_so)``, never on which interaction type is being
    probed (``_s4_make_solver`` writes no interaction file at v=0
    regardless of ``kind``), so this is computed ONCE per (variant, nmat)
    rather than once per type, the runtime-budget optimization this
    module's docstring names."""
    if variant == "norb1":
        solver, green_info = _s4_make_solver(
            _S4_1_L, _S4_1_NORB, _S4_1_HOP, _S4_1_EPS, _S4_1_LSO,
            _S4_1_T, _S4_1_MU, "CoulombIntra", 0.0, nmat)
    else:
        solver, green_info = _s4_make_solver(
            _S4_2_L, _S4_2_NORB, _S4_2_HOP, _S4_2_EPS, _S4_2_LSO,
            _S4_2_T, _S4_2_MU, "CoulombInter", 0.0, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s4_gate_out_"))
    chiq_pm = np.asarray(green_info["chiq_pm"])
    return chiq_pm[chiq_pm.shape[0] // 2]


@functools.lru_cache(maxsize=None)
def _s4_type_chiq_pm(variant, kind, v, nmat, phase):
    """Runs the PRODUCTION ``solve()`` entry point (post-Phase-S) for one
    (variant, kind, coupling, Nmat) combination and returns the static
    (Omega=0) slice of ``green_info["chiq_pm"]`` -- the anti-vacuity rule:
    every granule below drives this, never a hand-assembled chi_dressed."""
    if variant == "norb1":
        solver, green_info = _s4_make_solver(
            _S4_1_L, _S4_1_NORB, _S4_1_HOP, _S4_1_EPS, _S4_1_LSO,
            _S4_1_T, _S4_1_MU, kind, v, nmat, phase)
    else:
        solver, green_info = _s4_make_solver(
            _S4_2_L, _S4_2_NORB, _S4_2_HOP, _S4_2_EPS, _S4_2_LSO,
            _S4_2_T, _S4_2_MU, kind, v, nmat, phase)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s4_gate_out_"))
    chiq_pm = np.asarray(green_info["chiq_pm"])
    return chiq_pm[chiq_pm.shape[0] // 2]


def _s4_pred_at_nmat(variant, kind, v1, nmat, phase=1.0):
    """Richardson-extrapolated d(chiq_pm)/dv at v -> 0 through PRODUCTION
    solve() (the solver-side prediction is NOT linear in the declared
    coupling -- the RPA equation resums it -- exactly Gate S0's own
    ``_s0_pred_at_nmat`` reasoning), sharing the v=0 point across every
    type via ``_s4_free_chiq_pm`` (the runtime-budget optimization)."""
    f0 = _s4_free_chiq_pm(variant, nmat)
    f1 = _s4_type_chiq_pm(variant, kind, _s4_scaled(v1), nmat, phase)
    f2 = _s4_type_chiq_pm(variant, kind, _s4_scaled(2 * v1), nmat, phase)
    return 2 * (f1 - f0) / v1 - (f2 - f0) / (2 * v1)


@functools.lru_cache(maxsize=1)
def _s4_norb1_fixture():
    fx = ed_oracle_util.EDFixture(L=_S4_1_L, norb=_S4_1_NORB, t=_S4_1_HOP,
                                    eps=_S4_1_EPS, T=_S4_1_T, mu=_S4_1_MU)
    return fx, _h1_with_spin_orbit(fx, _S4_1_LSO)


@functools.lru_cache(maxsize=1)
def _s4_norb2_fixture():
    fx = ed_oracle_util.EDFixture(L=_S4_2_L, norb=_S4_2_NORB, t=_S4_2_HOP,
                                    eps=_S4_2_EPS, T=_S4_2_T, mu=_S4_2_MU)
    return fx, _h1_with_spin_orbit(fx, _S4_2_LSO)


def _s4_ed_derivative(variant, kind, v1, phase=1.0):
    """(D_v1, D_vhalf): the Richardson pair of TotalNED's full-minus-HF-
    only ``bond_correlator_transverse`` derivative for ONE on-site type,
    on the given fixture variant, frame-mapped via ``_s0_frame_map`` --
    REUSED VERBATIM (interfaces: "REUSE it verbatim"; the map is pure
    index/axis algebra -- q-reversal and the Kubo pair-slot swap -- with
    no dependence on the particular (L, norb) it is applied to)."""
    fx, h1_so = (_s4_norb1_fixture() if variant == "norb1"
                 else _s4_norb2_fixture())
    norb = _S4_1_NORB if variant == "norb1" else _S4_2_NORB
    # bond_correlator_transverse's norb==1 shorthand takes bare (R,) --
    # TotalNED._channel_orbitals rejects the (R, a, b) triple there (it
    # unpacks assuming exactly one element).
    channels = ([(0,)] if norb == 1
                else [(0, a, b) for a in range(norb) for b in range(norb)])

    @functools.lru_cache(maxsize=None)
    def X_of_v(v):
        terms = _s4_terms_of_v(fx, kind, _s4_scaled(v), phase)
        full = ed_oracle_util.TotalNED(
            fx, terms=terms, h1=h1_so).bond_correlator_transverse(channels)
        hf_h1 = _hf_h1_from_terms_at(fx, terms, h1_so)
        hf_only = ed_oracle_util.TotalNED(
            fx, terms=(), h1=hf_h1).bond_correlator_transverse(channels)
        return (full - hf_only).reshape(fx.L, norb, norb, norb, norb)

    D_v1 = _s0_frame_map(ed_oracle_util.richardson(X_of_v, v1))
    D_vhalf = _s0_frame_map(ed_oracle_util.richardson(X_of_v, v1 / 2))
    return D_v1, D_vhalf


def _s4_adjudicate(variant, kind, label, phase=1.0, v1=CAMPAIGN_V1,
                    zero_mask_all=False):
    """One granule: production ``chiq_pm`` (Richardson, two Nmat) vs
    ``TotalNED`` (Richardson, two FD steps), via ``adjudicate_granule``.
    ``zero_mask_all=True`` marks the ENTIRE grid as a zero cell (the H2/
    Task6 PASS-ZERO-control pattern -- used only once a type's content has
    been independently MEASURED as zero, never assumed a priori)."""
    D_ed_v1, D_ed_vhalf = _s4_ed_derivative(variant, kind, v1, phase)
    D_pred_1 = _s4_pred_at_nmat(variant, kind, v1, _S4_NMAT1, phase)
    D_pred_2 = _s4_pred_at_nmat(variant, kind, v1, _S4_NMAT2, phase)
    zero_mask = np.full(D_pred_2.shape, zero_mask_all, dtype=bool)
    return ed_oracle_util.adjudicate_granule(
        D_ed_v1, D_ed_vhalf, D_pred_1, D_pred_2, zero_mask, label)


class TestTask4GranuleCampaign(unittest.TestCase):
    """Task 4, Step 1: the spin-orbit granule campaign. On the t_so
    fixture FAMILY (a norb=1 variant, Task 1's own Sz-breaking
    construction-contract fixture verbatim; a norb=2 variant, Gate S0's
    own (L, HOP, EPS, LSO, T, MU) verbatim), first-order adjudication of
    PRODUCTION ``solve()``'s ``chiq_pm`` against ``TotalNED``, per on-site
    interaction type, via ``adjudicate_granule`` (``_s4_adjudicate``,
    which drives the anti-vacuity rule: every granule below runs the
    public ``solve()`` entry point, never a hand-assembled dressed
    tensor -- mirroring Gate S2's own production leg).

    **MEASURED VERDICT, all SIX types PASS** (delta_rich/delta_nmat/tol/
    max_signal, this module's own measurements):

        norb1/CoulombIntra  delta_rich=1.81e-08 delta_nmat=1.28e-07
            tol=1.28e-06 max_signal=4.87e-02  (margin ~3.8e4 x tol)
        norb2/CoulombInter  delta_rich=1.16e-07 delta_nmat=1.47e-07
            tol=1.47e-06 max_signal=3.95e-02  (margin ~2.7e4 x tol)
        norb2/Exchange      delta_rich=8.80e-09 delta_nmat=1.38e-07
            tol=1.38e-06 max_signal=3.35e-02  (margin ~2.4e4 x tol)
        norb2/PairHop       delta_rich=9.10e-09 delta_nmat=1.36e-07
            tol=1.36e-06 max_signal=3.33e-02  (margin ~2.4e4 x tol)
        norb2/Hund          delta_rich=1.62e-08 delta_nmat=5.56e-10
            tol=1.62e-07 max_signal=1.44e-04  (margin ~8.9e2 x tol)
        norb2/PairLift      delta_rich=1.08e-09 delta_nmat=2.82e-09
            tol=2.82e-08 max_signal=1.32e-03  (margin ~4.7e4 x tol)

    **FINDING (spec amendment request -- the brief's central "derive
    first" question, PLUS one the brief did not anticipate)**: the brief
    named PairLift as needing derivation ("an Sz-breaking quartic that
    TotalNED CAN represent, unlike SectorED... derive first whether its
    transverse first-order content is genuinely zero... or nonzero") and
    named Hund as an ordinary PASS-ZERO control alongside it. MEASURED:
    BOTH types are genuinely, measurably NONZERO on this Sz-breaking
    fixture -- production correctly reproduces both (clean PASS, not
    PASS-ZERO, for either). Hund is NOT a PASS-ZERO control here.

    This is not a production bug -- see the on-site table's own home,
    ``RPA._assemble_transverse_vertex`` (src/hwave/solver/rpa.py), whose
    docstring proves ``ham_pm`` (the SEPARATE, narrower vertex the
    Sz-conserving/reduced ladder path -- ``_build_transverse_channel``'s
    spin-free/spin-diag branches -- actually uses) has Hund's and
    PairLift's entries STRUCTURALLY zero regardless of any one-body
    field: Hund's declared coefficient lives entirely in the "uuuu,
    dddd" blocks of the raw interaction tensor and PairLift's in "uddu,
    duud", NEITHER of which ``_assemble_transverse_vertex`` ever reads
    (it draws only on "uudd, dduu" and "udud, dudu" -- see that method's
    own "Measured block occupancy" table). That argument is about
    ``ham_pm`` alone and says nothing about
    ``_extract_transverse_from_dressed``'s DIFFERENT slice of the FULL
    RPA-DRESSED general tensor: once ``chi0`` itself carries genuine
    cross-spin content (t_so's own doing, already the basis of Gate
    S1/S2), the RPA matrix inversion mixes ALL vertex blocks into the
    sliced output, including Hund's/PairLift's -- there is no remaining
    structural reason for either to stay zero once Sz is broken, and the
    measurement above confirms it is not.

    **THE ISOLATION CONTROL** (``TestTask4TsoZeroControl`` below): re-running
    BOTH Hund and PairLift on the IDENTICAL norb=2 fixture with ``t_so``
    set to EXACTLY ``0.0`` reproduces H2's/production's own established
    on-site table value of PRECISELY zero (production ``chiq_pm``
    BIT-IDENTICAL 0.0 across both Nmat resolutions, ED-side Richardson
    derivative at the noise floor, ~5-6e-09 -- matching H2's own PairLift
    PASS-ZERO magnitude, ~1.3e-09, to within an order of magnitude) --
    isolating ``t_so`` as the SOLE responsible mechanism, not a
    construction defect in this campaign's own fixture/term builders.

    **SPEC AMENDMENT REQUESTED** (recorded here, adjudicated by a human
    coordinator, not self-applied): the spec's on-site table (H2,
    ``docs/superpowers/specs/2026-08-15-bond-transverse-design.md``,
    "The verified on-site algebra") and ``_assemble_transverse_vertex``'s
    own docstring table should note explicitly that "Hund/PairLift: 0" is
    an Sz-CONSERVING-only statement about ``ham_pm`` (the narrow ladder
    vertex), not a statement about the spinful ``_extract_transverse_
    from_dressed`` path once a genuine Sz-breaking one-body term (t_so)
    is present -- where BOTH types carry real, production-correctly-
    reproduced transverse content instead.

    Fixture provenance: norb=1 -- ``TestTotalNEDSzBreakingDiscrimination``'s
    own fixture (L=2, t=0.4-0.1j, eps=0.05, T=0.5, mu=0.1,
    t_so=0.3+0.05j). norb=2 -- Gate S0's own fixture (L=3, complex HOP,
    EPS=(0.10,-0.05), T=0.5, mu=0.2, t_so=0.3+0.05j). Coupling: the
    Richardson finite-difference variable v in {0, CAMPAIGN_V1/2,
    CAMPAIGN_V1, 2*CAMPAIGN_V1} scaled by ``_S4_UNIT=0.3`` before being
    declared (``_s4_scaled``). Nmat=(512, 1024) -- a deliberate reduction
    from Gate S0's (1024, 2048) for the runtime budget (six independent
    fixtures here vs Gate S0's one); every PASS margin above (all >= 887x
    tol) confirms this reduction did not cost adjudication validity.
    PairHop uses Gate S0's own genuinely complex phase
    (``_S4_2_PAIRHOP_PHASE = exp(0.4j)``); every other type uses a real
    coupling (Exchange/PairLift real-``v`` collapse to
    ``oracle._terms_for``'s on-site convention exactly, per
    ``_s4_terms_of_v``'s own docstring).
    """

    def test_norb1_coulombintra(self):
        rec = _s4_adjudicate("norb1", "CoulombIntra",
                              "task4/norb1/CoulombIntra")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_norb2_coulombinter(self):
        rec = _s4_adjudicate("norb2", "CoulombInter",
                              "task4/norb2/CoulombInter")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_norb2_exchange(self):
        rec = _s4_adjudicate("norb2", "Exchange", "task4/norb2/Exchange")
        self.assertEqual(rec["status"], "PASS", rec)

    def test_norb2_pairhop_complex(self):
        rec = _s4_adjudicate("norb2", "PairHop", "task4/norb2/PairHop",
                              phase=_S4_2_PAIRHOP_PHASE)
        self.assertEqual(rec["status"], "PASS", rec)

    def test_norb2_hund(self):
        rec = _s4_adjudicate("norb2", "Hund", "task4/norb2/Hund")
        self.assertEqual(
            rec["status"], "PASS", "measured status={} -- see class "
            "docstring's FINDING: Hund is NOT a PASS-ZERO control on the "
            "Sz-breaking t_so fixture (record={})".format(
                rec["status"], rec))

    def test_norb2_pairlift(self):
        rec = _s4_adjudicate("norb2", "PairLift", "task4/norb2/PairLift")
        self.assertEqual(
            rec["status"], "PASS", "measured status={} -- see class "
            "docstring's FINDING (the brief's own 'derive first' "
            "question): PairLift is NOT a PASS-ZERO control on the "
            "Sz-breaking t_so fixture (record={})".format(
                rec["status"], rec))


# ---------------------------------------------------------------------------
# The isolation control: re-run Hund/PairLift on the IDENTICAL norb=2
# fixture with t_so forced to EXACTLY 0.0, confirming their PASS-ZERO
# behaviour reproduces the H2/production-established on-site table value
# there -- pinning t_so as the mechanism responsible for
# TestTask4GranuleCampaign's nonzero measurement, not a defect in this
# campaign's own term/file-writing construction. A standalone (non-cached-
# shared) construction, deliberately NOT threading through
# _s4_norb2_fixture/_s4_free_chiq_pm/_s4_type_chiq_pm (which hardcode the
# module's own t_so): _s4_make_solver already accepts lso as an explicit
# argument, so this reuses it directly with lso=0.0 rather than mutating
# any shared module state.
# ---------------------------------------------------------------------------


def _s4_tso0_ed_derivative(kind, v1):
    fx = ed_oracle_util.EDFixture(L=_S4_2_L, norb=_S4_2_NORB, t=_S4_2_HOP,
                                    eps=_S4_2_EPS, T=_S4_2_T, mu=_S4_2_MU)
    channels = [(0, a, b) for a in range(_S4_2_NORB)
                for b in range(_S4_2_NORB)]

    @functools.lru_cache(maxsize=None)
    def X_of_v(v):
        terms = _s4_terms_of_v(fx, kind, _s4_scaled(v))
        full = ed_oracle_util.TotalNED(
            fx, terms=terms).bond_correlator_transverse(channels)
        hf_h1 = _hf_h1_from_terms_at(fx, terms, fx.build_h1())
        hf_only = ed_oracle_util.TotalNED(
            fx, terms=(), h1=hf_h1).bond_correlator_transverse(channels)
        return (full - hf_only).reshape(fx.L, _S4_2_NORB, _S4_2_NORB,
                                         _S4_2_NORB, _S4_2_NORB)

    D_v1 = _s0_frame_map(ed_oracle_util.richardson(X_of_v, v1))
    D_vhalf = _s0_frame_map(ed_oracle_util.richardson(X_of_v, v1 / 2))
    return D_v1, D_vhalf


@functools.lru_cache(maxsize=None)
def _s4_tso0_chiq_pm(kind, v, nmat):
    solver, green_info = _s4_make_solver(
        _S4_2_L, _S4_2_NORB, _S4_2_HOP, _S4_2_EPS, 0.0 + 0.0j,
        _S4_2_T, _S4_2_MU, kind, v, nmat)
    solver.solve(green_info, tempfile.mkdtemp(prefix="s4_gate_out_"))
    chiq_pm = np.asarray(green_info["chiq_pm"])
    return chiq_pm[chiq_pm.shape[0] // 2]


def _s4_tso0_pred_at_nmat(kind, v1, nmat):
    f0 = _s4_tso0_chiq_pm(kind, 0.0, nmat)
    f1 = _s4_tso0_chiq_pm(kind, _s4_scaled(v1), nmat)
    f2 = _s4_tso0_chiq_pm(kind, _s4_scaled(2 * v1), nmat)
    return 2 * (f1 - f0) / v1 - (f2 - f0) / (2 * v1)


def _s4_tso0_adjudicate(kind, label, v1=CAMPAIGN_V1):
    D_ed_v1, D_ed_vhalf = _s4_tso0_ed_derivative(kind, v1)
    D_pred_1 = _s4_tso0_pred_at_nmat(kind, v1, _S4_NMAT1)
    D_pred_2 = _s4_tso0_pred_at_nmat(kind, v1, _S4_NMAT2)
    zero_mask = np.ones(D_pred_2.shape, dtype=bool)
    return ed_oracle_util.adjudicate_granule(
        D_ed_v1, D_ed_vhalf, D_pred_1, D_pred_2, zero_mask, label)


class TestTask4TsoZeroControl(unittest.TestCase):
    """The isolation control for ``TestTask4GranuleCampaign``'s Hund/
    PairLift finding: on the IDENTICAL norb=2 fixture with ``t_so`` set to
    EXACTLY ``0.0`` (every other parameter unchanged), BOTH types return
    to PASS-ZERO -- the production-side finite-difference derivative of
    ``chiq_pm`` with respect to the declared coupling (the quantity
    ``_s4_tso0_pred_at_nmat`` builds and ``delta_nmat`` compares across
    resolutions; NOT raw ``chiq_pm`` itself, which is the nonzero v=0
    bare-bubble baseline) is measured BIT-IDENTICAL ``0.0`` at both Nmat
    resolutions (asserted directly, not merely inferred from
    ``delta_nmat=0.0``), ED-side Richardson derivative at the noise
    floor:

        Hund      delta_rich=1.64e-08  ED max|D_vhalf|=5.48e-09
        PairLift  delta_rich=6.44e-09  ED max|D_vhalf|=5.67e-09

    matching H2's own PairLift PASS-ZERO measurement (``ED derivative
    magnitude ~1.3e-09``, ``test_bond_transverse_ed.py``) to within an
    order of magnitude, and reproducing ``_assemble_transverse_vertex``'s
    own documented table entries ("Hund J: 0", "PairLift J: 0") exactly.
    This isolates ``t_so`` as the SOLE mechanism responsible for
    ``TestTask4GranuleCampaign``'s nonzero measurement -- not a defect in
    this campaign's own fixture/term-builder construction."""

    def _assert_derivative_bit_identical_zero(self, kind):
        """Pin the docstring's bit-identity claim directly against the
        production-side array the adjudication is built from, rather than
        trusting ``status == "PASS-ZERO"`` alone to imply it.

        Note what this is checking: production ``chiq_pm`` itself (the
        v=0 baseline, ``f0``) is the nonzero bare-bubble value -- it is
        the coupling-derivative ``_s4_tso0_pred_at_nmat`` builds from it
        (the same finite-difference quantity that feeds ``delta_nmat``)
        that the docstring's "measured BIT-IDENTICAL 0.0" refers to.
        Checking it directly at each Nmat, rather than relying on
        ``delta_nmat == 0.0``, matters because ``delta_nmat`` is only a
        difference between the two resolutions: it would read 0.0 just
        as well if both resolutions agreed on some shared NONZERO value,
        which would not correspond to "Hund/PairLift: 0" at all.
        """
        for nmat in (_S4_NMAT1, _S4_NMAT2):
            D_pred = _s4_tso0_pred_at_nmat(kind, CAMPAIGN_V1, nmat)
            self.assertTrue(
                np.array_equal(D_pred, np.zeros_like(D_pred)),
                "kind={} nmat={}: production derivative is not "
                "bit-identical to 0.0".format(kind, nmat))

    def test_hund_returns_to_pass_zero_at_tso_zero(self):
        rec = _s4_tso0_adjudicate("Hund", "task4/control/tso0/Hund")
        self.assertEqual(rec["status"], "PASS-ZERO", rec)
        self.assertEqual(rec["delta_nmat"], 0.0, rec)
        self._assert_derivative_bit_identical_zero("Hund")

    def test_pairlift_returns_to_pass_zero_at_tso_zero(self):
        rec = _s4_tso0_adjudicate("PairLift", "task4/control/tso0/PairLift")
        self.assertEqual(rec["status"], "PASS-ZERO", rec)
        self.assertEqual(rec["delta_nmat"], 0.0, rec)
        self._assert_derivative_bit_identical_zero("PairLift")


class TestTask4SVDSensitivityGate(unittest.TestCase):
    """Task 4, Step 2: the SVD sensitivity gate across each fixture's
    direction set (the #151 replica pattern, ``_w0_direction_set_sv_gate``
    -- imported from ``test_bond_transverse_ed.py`` rather than
    duplicated, since it is a plain, campaign-label-free utility).

    **norb=2 fixture -- SPLIT into two direction-set gates, a recorded
    rank argument (the "exclude/group" decision Step 2 asks for)**: a
    first construction combined all FIVE norb=2 types (CoulombInter,
    Exchange, PairHop, Hund, PairLift) into ONE 5-column gate and measured
    ``sv_ratio=1.96e-03`` -- technically clearing ``SENS_SV_FLOOR=1e-03``,
    but only by a factor of ~2, an uncomfortably thin margin. The cause is
    NOT a near-degenerate (linearly dependent) direction pair -- the five
    types occupy almost disjoint (a,c,b,d) table cells by construction, so
    they are NOT collinear -- but a raw SCALE mismatch of ~2-3 orders of
    magnitude between the three "directly active" on-site types
    (CoulombInter/Exchange/PairHop, max_signal ~0.03-0.04) and the two
    "t_so-MEDIATED" types this task's own finding surfaced (Hund/
    PairLift, max_signal ~1e-4 to 1e-3): ``sensitivity_rank``/
    ``_w0_direction_set_sv_gate`` never normalize columns, so a scale
    imbalance this large squeezes the smallest singular value even when
    every column is genuinely independent (the SAME failure mode the
    module docstring's "why a naive 4-type combined gate is DEGENERATE"
    note warns about for mixing unrelated topologies in one gate).
    SPLITTING by physical origin -- {CoulombInter, Exchange, PairHop}
    (the established, directly-active on-site types) vs {Hund, PairLift}
    (the newly-surfaced t_so-mediated types) -- gives two well-conditioned
    gates instead of one thin one:

        active3       sv_ratio=8.69e-01  sigma_min=7.22e-02
                       100*delta_rich_max=1.16e-05  (margin ~6200x)
        tso_mediated2 sv_ratio=9.75e-02  sigma_min=2.57e-04
                       100*delta_rich_max=1.62e-06  (margin ~159x)

    Both clear BOTH gates (``sv_ratio >= SENS_SV_FLOOR`` and
    ``sigma_min >= 100*delta_rich_max``) with comfortable margins. No
    direction pair was excluded as exactly degenerate -- none was found;
    the split is a grouping-by-scale choice, not an exclusion.

    **norb=1 fixture**: only ONE on-site type applies (CoulombIntra --
    every other type names an orbital PAIR by construction). A
    single-column "direction set" has no rank to discriminate (SVD of a
    single nonzero column trivially returns ``sv_ratio=1.0``) -- run here
    for completeness/uniformity with the norb=2 gates, but recorded
    explicitly as a DEGENERATE-BY-CONSTRUCTION case, not a genuine
    multi-direction rank test.
    """

    def test_norb2_active_types_direction_set(self):
        records = {
            "CoulombInter": _s4_adjudicate(
                "norb2", "CoulombInter", "task4/sv/norb2/CoulombInter"),
            "Exchange": _s4_adjudicate(
                "norb2", "Exchange", "task4/sv/norb2/Exchange"),
            "PairHop": _s4_adjudicate(
                "norb2", "PairHop", "task4/sv/norb2/PairHop",
                phase=_S4_2_PAIRHOP_PHASE),
        }
        for name, rec in records.items():
            self.assertEqual(rec["status"], "PASS", (name, rec))
        gate = _w0_direction_set_sv_gate(records, "task4/norb2_active3")
        self.assertGreaterEqual(
            gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR,
            "task4/norb2_active3 direction-set SVD gate: sv_ratio={:.3e} "
            "below floor {:.3e}".format(gate["sv_ratio"],
                                         ed_oracle_util.SENS_SV_FLOOR))
        self.assertGreaterEqual(
            gate["sigma_min"], 100.0 * gate["delta_rich_max"],
            "task4/norb2_active3 direction-set SVD gate: sigma_min={:.3e} "
            "below 100*delta_rich_max={:.3e}".format(
                gate["sigma_min"], 100.0 * gate["delta_rich_max"]))

    def test_norb2_tso_mediated_direction_set(self):
        records = {
            "Hund": _s4_adjudicate("norb2", "Hund", "task4/sv/norb2/Hund"),
            "PairLift": _s4_adjudicate(
                "norb2", "PairLift", "task4/sv/norb2/PairLift"),
        }
        for name, rec in records.items():
            self.assertEqual(rec["status"], "PASS", (name, rec))
        gate = _w0_direction_set_sv_gate(records, "task4/norb2_tso_mediated2")
        self.assertGreaterEqual(
            gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR,
            "task4/norb2_tso_mediated2 direction-set SVD gate: "
            "sv_ratio={:.3e} below floor {:.3e}".format(
                gate["sv_ratio"], ed_oracle_util.SENS_SV_FLOOR))
        self.assertGreaterEqual(
            gate["sigma_min"], 100.0 * gate["delta_rich_max"],
            "task4/norb2_tso_mediated2 direction-set SVD gate: "
            "sigma_min={:.3e} below 100*delta_rich_max={:.3e}".format(
                gate["sigma_min"], 100.0 * gate["delta_rich_max"]))

    def test_norb1_single_direction_degenerate_by_construction(self):
        records = {
            "CoulombIntra": _s4_adjudicate(
                "norb1", "CoulombIntra", "task4/sv/norb1/CoulombIntra"),
        }
        self.assertEqual(records["CoulombIntra"]["status"], "PASS")
        gate = _w0_direction_set_sv_gate(records, "task4/norb1_single")
        # A single active column: sv_ratio is trivially 1.0 by
        # construction (one singular value, ratio of itself to itself) --
        # asserted here as a documented, EXPECTED degenerate-rank result,
        # not a genuine rank/conditioning finding.
        self.assertEqual(gate["sv_ratio"], 1.0, gate)
        self.assertEqual(len(gate["singular_values"]), 1, gate)


# ---------------------------------------------------------------------------
# Pre-push hardening round: `spinful_vertex_exchange = false` is a
# ring-only compatibility switch (reproduces the pre-#137 ring-only
# numbers -- the longitudinal solve is NOT antisymmetrized under it).
# Under `calc_type='ring+ladder'` the transverse channel is extracted by
# SLICING that same longitudinal solve (`_extract_transverse_from_dressed`,
# reused verbatim, Gate S2 above), so the opt-out combination would
# silently produce a `chiq_pm` matching neither the legacy ring-only
# transverse assembly (`_assemble_transverse_vertex`) nor the advertised
# full-space dressing. `RPA` now rejects the combination at solve time;
# `calc_type='ring'` never reaches the extraction, so the opt-out keeps
# doing its documented job there and remains accepted.
# ---------------------------------------------------------------------------


class TestSpinfulVertexExchangeOptOutRingLadderRejection(unittest.TestCase):
    _NMAT = 16

    def _solver(self, calc_type, spinful_vertex_exchange):
        return _s4_make_solver(
            _S4_1_L, _S4_1_NORB, _S4_1_HOP, _S4_1_EPS, _S4_1_LSO,
            _S4_1_T, _S4_1_MU, "CoulombIntra", 0.3, self._NMAT,
            calc_type=calc_type,
            spinful_vertex_exchange=spinful_vertex_exchange)

    def test_ring_ladder_rejects_the_opt_out_combination(self):
        solver, green_info = self._solver("ring+ladder", False)
        with self.assertRaises(ValueError) as cm:
            solver.solve(green_info,
                         tempfile.mkdtemp(prefix="s5_reject_out_"))
        msg = str(cm.exception)
        self.assertIn("spinful_vertex_exchange", msg)
        self.assertIn("ring+ladder", msg)
        # a construction/solve-time rejection must not leave a partially
        # populated green_info behind (chiq/chiq_pm are popped at entry
        # and never repopulated before the raise).
        self.assertNotIn("chiq", green_info)
        self.assertNotIn("chiq_pm", green_info)

    def test_ring_mode_still_accepts_the_opt_out_and_differs_from_default(
            self):
        solver_false, gi_false = self._solver("ring", False)
        solver_false.solve(
            gi_false, tempfile.mkdtemp(prefix="s5_ring_false_out_"))
        chiq_false = np.asarray(gi_false["chiq"])
        self.assertTrue(np.all(np.isfinite(chiq_false)))

        solver_true, gi_true = self._solver("ring", True)
        solver_true.solve(
            gi_true, tempfile.mkdtemp(prefix="s5_ring_true_out_"))
        chiq_true = np.asarray(gi_true["chiq"])
        self.assertTrue(np.all(np.isfinite(chiq_true)))

        self.assertEqual(chiq_false.shape, chiq_true.shape)
        # the flag must still do its documented ring-only job: false and
        # true must NOT produce the same longitudinal chiq
        self.assertFalse(
            np.array_equal(chiq_false, chiq_true),
            "spinful_vertex_exchange must still change the ring-only "
            "chiq under calc_type='ring': false and true produced "
            "identical output")


if __name__ == "__main__":
    unittest.main()
