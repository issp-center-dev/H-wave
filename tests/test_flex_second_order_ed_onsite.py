#!/usr/bin/env python3

"""G3: exact diagonalisation of the ON-SITE second order, every interaction type.

On the single-site three-orbital model of
``tests/test_flex_sopt_index_order.py`` (three orbitals with a complex
hopping loop, so no orbital index order or phase can hide), the O(x y)
coefficient of the EXACT fluctuation remainder

    Sigma_ED - Sigma_HF[rho_ED]

must equal the production ``Sigma_fluct`` coefficient of the general FLEX
path with ``flex_second_order = "local"`` on ``CellShape = [1, 1, 1]``.

Why the remainder is the right thing to compare
-----------------------------------------------
The exact self-energy is ``Sigma_HF[G] + Sigma_2[G] + O(v^3)`` with the HF
functional evaluated on the EXACT density, so subtracting
``Sigma_HF[rho_ED]`` removes the first order AND the second-order piece
that the first-order functional picks up from the density response. What is
left at O(v^2) is the second-order skeleton on the bare propagator --
exactly what the production ``Sigma_fluct`` carries at that order (its RPA
ladder starts contributing only at third order).

Independence
------------
The Hamiltonian, the seven operator definitions, the Lehmann evaluation AND
the first-order functional that is subtracted are written HERE, from the
operator table of the design (spec 2.2); nothing on the ED side reads
``hwave.solver.second_order``. The subtracted first order is
:func:`wick_hf_sigma`, the Wick derivative of the ED monomials themselves
(:func:`_h_int_terms`, pinned against the dense operators of
:func:`_h_int` at zero tolerance), so the remainder whose second order the
gate compares is not built with any part of the object under test.

The production first-order functional
(``hf_first_order(compile_onsite(...))``) appears only as a CONVENTION PIN,
in two places, neither of which feeds the second-order comparison:
:meth:`TestWickFunctional.test_independent_functional_equals_the_production_first_order`
(the two functionals agree on several random Hermitian densities, block
diagonal in spin and generic) and
:meth:`TestG3.test_first_order_remainder_vanishes` (the compiled ``Gamma``
IS the first order of this Hamiltonian -- which is what makes the
mirrored-row weight a measurement; see below).

The mirrored-row weight, determined empirically
-----------------------------------------------
The reader delivers both ordered rows ``(a, b)`` and ``(b, a)`` of a
CoulombInter/Hund/Ising/Exchange/PairLift bond, and both rows name the SAME
operator; the solver's reversal closure therefore gives each ordered row
HALF the bond (``_MIRRORED_TYPES`` of ``hwave.solver.second_order``). The
ED Hamiltonian here follows that convention, and
:meth:`TestG3.test_first_order_remainder_vanishes` -- which subtracts the
PRODUCTION functional, precisely so that it can be two-sided (the test-side
functional tracks whichever Hamiltonian it is handed, so it cannot
discriminate a weight) -- PINS it from both sides:
with the half weight the O(v) remainder is <= 2.4e-4 of the self-energy's
own linear size for every type, and with the full weight it is ~0.5 -- half
the first order left unsubtracted, because that Hamiltonian carries twice
the interaction (which would be a factor FOUR at second order) -- for the
four mirrored types that have a first order at all.
``CoulombIntra`` has no transposed partner and ``PairHop``'s transposed row
is its Hermitian partner (a different operator), so neither is halved and
neither moves between the two conventions.

Stencils and tolerances (measured, not guessed)
-----------------------------------------------
The coefficient of a type is extracted with the two-point stencil
``(f(2x) - 2 f(x)) / (2 x^2)`` and of a pair with
``(f(x,y) - f(x,0) - f(0,y)) / (x y)``; both have a leading O(x) error (the
third-order terms), so a two-level refinement ``2 c(x/2) - c(x)`` is
applied, NOT the ``(4 c(x/2) - c(x)) / 3`` of an O(x^2) error -- measured
at ``x = 0.05, nmat = 2048``, worst type over the seven: 1.5e-2 with the
O(x) refinement against 9.2e-2 with the O(x^2) one, and only the former
shrinks as ``x`` is reduced.

The tolerance is the 2e-3 of the brief and is NOT loosened; the stencil is
refined instead. Two error sources set the working point (both measured
with the correct refinement, worst type over the seven):

* the finite ``x`` truncation, which shrinks as ``x^2`` after the
  refinement: 1.5e-2 at ``x = 0.05``, 4.4e-3 at 0.025 and 1.6e-3 at 0.0125
  (``nmat = 2048``), and at 6.25e-3 it is already below the next source
  (5.6e-4, which halving ``x`` again does not improve);
* the production kernel's FINITE Matsubara window, which is a fixed
  fraction of the coefficient (it does not shrink with ``x``) and grows
  towards the edges of the grid: at ``nmat = 4096`` it is 7.4e-4 over the
  innermost eighth of the grid and 1.5e-3 over the brief's inner HALF.

Hence ``nmat = 4096``, ``x = 6.25e-3`` and a comparison window of the
innermost eighth: every entry of the table then agrees to <= 7.4e-4, a
factor 2.7 inside the tolerance. (The inner half would still pass, at 1.5e-3,
but with a margin too thin to be a gate.)
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

from tests.heavy_tests import heavy
from tests import test_flex_sopt_index_order as ref
from tests.test_second_order_factors import _write_wan

C, CD, NORB, BETA, MU, H0 = ref.C, ref.CD, ref.NORB, ref.BETA, ref.MU, ref.H0

UP, DN = 0, 1

#: Types whose transposed row carries the SAME operator, so each declared
#: ordered row is half of the bond (see the module docstring; the weight is
#: pinned by ``test_first_order_remainder_vanishes``, not assumed).
_MIRRORED_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift")

_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange", "PairHop", "PairLift")
_PAIRS = [("CoulombIntra", "CoulombInter"), ("CoulombIntra", "Hund"), ("CoulombIntra", "Ising"),
          ("CoulombIntra", "Exchange"), ("CoulombIntra", "PairHop"), ("CoulombIntra", "PairLift"),
          ("CoulombInter", "Hund")]

#: Cross terms that vanish identically by SPIN algebra (the argument, and
#: the measurement, of ``tests/test_flex_second_order_sopt.py::_KNOWN_ZERO``):
#: ``CoulombIntra``'s only monomial is ``n_{a up} n_{a dn}``, so at second
#: order the U vertex puts one up and one down leg on EACH side of the
#: skeleton and the partner vertex must do the same. On-site ``Hund`` has
#: only same-spin monomials and on-site ``PairLift`` only spin-paired ones
#: (both legs up on one side, both down on the other), so neither can, and
#: both cross terms are zero for every input -- not merely small here.
_KNOWN_ZERO = frozenset({("CoulombIntra", "Hund"), ("CoulombIntra", "PairLift")})

#: Types with NO first order at all (see
#: :meth:`TestG3.test_pairlift_first_order_is_zero`): the first-order gate
#: can confirm that both sides vanish for them, but it cannot discriminate
#: the mirrored-row weight -- zero times two is still zero. PairLift's
#: weight is pinned by the second order alone (a factor two in H is a
#: factor four there).
_NO_FIRST_ORDER = ("PairLift",)

#: Working point of the heavy comparison (see the module docstring).
_NMAT = 4096
_X = 6.25e-3
_WINDOW = _NMAT // 16          # half-width of the comparison window, in frequencies
_TOL = 2.0e-3

#: Anti-vacuity floor on the ED coefficient of a compared entry: a relative
#: comparison against a numerically zero reference would pass for any
#: production value at all. Every non-``_KNOWN_ZERO`` entry measures
#: 6.8e-2 ... 1.4e-1, four orders above it.
_FLOOR = 1.0e-6

#: Half-width of the CONVERGENCE GUARD's window, in frequencies: the inner
#: eighth of the SMALLER of the two grids it compares (``nmat = 2048``), so
#: the two cutoffs are compared on the SAME absolute Matsubara frequencies.
_GUARD_WINDOW = _NMAT // 32

#: Bounds on the drifts the convergence guard measures, as fractions of the
#: coefficient (see :meth:`TestG3.test_extraction_is_converged`). Measured
#: on ``CoulombIntra``: the Matsubara cutoff moves the refined coefficient
#: by 7.2e-4 between ``nmat = 2048`` and 4096 -- the same size as the
#: window error the module docstring records, and a factor 2.8 inside
#: ``_TOL`` -- while halving the stencil step moves it by 5.4e-5 and then
#: 1.4e-5 (ratio 4.01, the O(x^2) the two-level refinement claims).
_GUARD_CUTOFF_BOUND = _TOL / 2.0
_GUARD_STEP_BOUND = _TOL / 10.0

#: Two-sided ceiling for a ``_KNOWN_ZERO`` cross term. Measured at
#: ``x = 6.25e-3``: ED 2.0e-6, production 1.2e-6 (U x J) and ED 1.8e-6,
#: production 3.3e-17 (U x PL), both shrinking as ``x^2`` -- pure stencil
#: truncation against pure-type coefficients of 6.8e-2 ... 1.4e-1.
_ZERO_CEIL = 1.0e-4


@contextlib.contextmanager
def _quiet():
    """Silence the solver's per-map INFO/WARNING narration (one block per
    production map, some seventy maps per heavy method)."""
    lg = logging.getLogger("hwave")
    old = lg.level
    lg.setLevel(logging.ERROR)
    try:
        yield
    finally:
        lg.setLevel(old)


# --------------------------------------------------------------------------
# The model: operators, Hamiltonian, Lehmann evaluation (production-free)
# --------------------------------------------------------------------------

def _n(m, s):
    return CD[ref._so(m, s)] @ C[ref._so(m, s)]


def _h_int(itype, a, b, v):
    """The operator of ONE declared row (spec 2.2 monomials), as a dense matrix."""
    cd = lambda m, s: CD[ref._so(m, s)]
    c = lambda m, s: C[ref._so(m, s)]
    with ref._quiet_matmul():
        if itype == "CoulombIntra":
            return v * (_n(a, UP) @ _n(a, DN))
        if itype == "CoulombInter":
            return v * ((_n(a, UP) + _n(a, DN)) @ (_n(b, UP) + _n(b, DN)))
        if itype == "Hund":
            return -v * (_n(a, UP) @ _n(b, UP) + _n(a, DN) @ _n(b, DN))
        if itype == "Ising":
            return v * ((_n(a, UP) - _n(a, DN)) @ (_n(b, UP) - _n(b, DN)))
        if itype == "Exchange":
            return v * (cd(a, UP) @ c(b, UP) @ cd(b, DN) @ c(a, DN)
                        + cd(a, DN) @ c(b, DN) @ cd(b, UP) @ c(a, UP))
        if itype == "PairHop":
            return v * (cd(a, UP) @ c(b, UP) @ cd(a, DN) @ c(b, DN))
        if itype == "PairLift":
            P = cd(a, UP) @ c(a, DN) @ cd(b, UP) @ c(b, DN)
            return v * (P + P.conj().T)
    raise ValueError(itype)


def _h_int_terms(itype, a, b, v):
    """The SAME monomials as :func:`_h_int`, as a list of
    ``(p, q, r, s, coeff)`` entries meaning ``coeff c^dag_p c_q c^dag_r c_s``
    on the model's mode index ``ref._so(orbital, spin)``.

    :meth:`TestWickFunctional.test_term_list_reproduces_the_ed_operators`
    multiplies these out and compares with :func:`_h_int` at 0 tolerance, so
    the two are the same operator by measurement, not by inspection."""
    so = ref._so
    if itype == "CoulombIntra":
        return [(so(a, UP), so(a, UP), so(a, DN), so(a, DN), v)]
    if itype == "CoulombInter":
        return [(so(a, s1), so(a, s1), so(b, s2), so(b, s2), v)
                for s1 in (UP, DN) for s2 in (UP, DN)]
    if itype == "Hund":
        return [(so(a, s1), so(a, s1), so(b, s1), so(b, s1), -v) for s1 in (UP, DN)]
    if itype == "Ising":
        return [(so(a, s1), so(a, s1), so(b, s2), so(b, s2), v if s1 == s2 else -v)
                for s1 in (UP, DN) for s2 in (UP, DN)]
    if itype == "Exchange":
        return [(so(a, s1), so(b, s1), so(b, 1 - s1), so(a, 1 - s1), v) for s1 in (UP, DN)]
    if itype == "PairHop":
        return [(so(a, UP), so(b, UP), so(a, DN), so(b, DN), v)]
    if itype == "PairLift":
        # ``v (P + P^dagger)`` exactly as :func:`_h_int` writes it -- the
        # SAME ``v`` on both placements, not its conjugate (every PairLift
        # coupling declared in this module is real, and the mirrored-row
        # closure folds a complex one to its real part in any case).
        return [(so(a, UP), so(a, DN), so(b, UP), so(b, DN), v),
                (so(b, DN), so(b, UP), so(a, DN), so(a, UP), v)]
    raise ValueError(itype)


def wick_hf_sigma(nmode, terms, rho):
    """The first-order (Hartree-Fock) self-energy of a quartic term list at
    an arbitrary density ``rho[p, q] = <c^dag_p c_q>``, from the Wick
    derivative alone.

    ``terms`` are ``(p, q, r, s, coeff)`` entries for
    ``coeff c^dag_p c_q c^dag_r c_s``. The functional is
    ``E = coeff [rho_pq rho_rs + rho_ps (delta_qr - rho_rq)]`` and
    ``Sigma_HF = dE / d rho`` symmetrised -- the engine of
    ``ed_oracle_util.hf_h1_from_terms``, evaluated away from the free
    density, which is what the exact-diagonalization gates need.

    This is the TEST SIDE's own first order: it is derived from the ED
    Hamiltonian's monomials and never reads the production compiler, so a
    remainder built with it is independent of ``hwave.solver.second_order``.
    The production functional is compared against it separately
    (:class:`TestWickFunctional`), as a convention pin rather than as part
    of the subtraction."""
    S = np.zeros((nmode, nmode), dtype=complex)
    for (p, q, r, s, c_) in terms:
        S[p, q] += c_ * rho[r, s]
        S[r, s] += c_ * rho[p, q]
        S[p, s] += c_ * ((1.0 if q == r else 0.0) - rho[r, q])
        S[r, q] += -c_ * rho[p, s]
    return 0.5 * (S + S.conj().T)


def _gen_to_mode(rho_gen):
    """A density on the generalised index ``s * NORB + a`` (the spin-block
    order of ``compile_onsite``) re-indexed onto the model's mode index
    ``ref._so(a, s) = a * 2 + s``."""
    M = 2 * NORB
    out = np.zeros((M, M), dtype=complex)
    for sp in range(2):
        for ap in range(NORB):
            for sq in range(2):
                for aq in range(NORB):
                    out[ref._so(ap, sp), ref._so(aq, sq)] = rho_gen[sp * NORB + ap, sq * NORB + aq]
    return out


def _mode_to_gen(m):
    """Inverse re-indexing of :func:`_gen_to_mode`."""
    M = 2 * NORB
    out = np.zeros((M, M), dtype=complex)
    for sp in range(2):
        for ap in range(NORB):
            for sq in range(2):
                for aq in range(NORB):
                    out[sp * NORB + ap, sq * NORB + aq] = m[ref._so(ap, sp), ref._so(aq, sq)]
    return out


def _rows(itype, v):
    """The DECLARED rows of one type, as the reader delivers them: the
    on-site ``CoulombIntra`` on every orbital, and both ordered rows of the
    orbital-1/orbital-2 bond for the two-orbital types (for ``PairHop`` the
    second row is the reverse hop, i.e. the Hermitian partner)."""
    if itype == "CoulombIntra":
        return [(0, 0, 0, m + 1, m + 1, v, 0.0) for m in range(NORB)]
    return [(0, 0, 0, 1, 2, v, 0.0), (0, 0, 0, 2, 1, v, 0.0)]


def _mode_terms(terms, mirrored_weight=0.5, rows_of=None):
    """The quartic term list of the WHOLE declaration (``terms``: a list of
    ``(itype, v)``), with the mirrored-row weight applied exactly as
    :func:`_hamiltonian` applies it to the dense operators."""
    rows_of = rows_of or _rows
    out = []
    for itype, v in terms:
        w = mirrored_weight if itype in _MIRRORED_TYPES else 1.0
        for (rx, ry, rz, a1, b1, vr, vi) in rows_of(itype, v):
            out.extend(_h_int_terms(itype, a1 - 1, b1 - 1, w * complex(vr, vi)))
    return out


def _hamiltonian(terms, mirrored_weight=0.5, rows_of=None):
    """``terms``: list of ``(itype, v)``. Each declared row contributes its
    operator; the rows of the mirrored types (which name the same operator
    twice) carry ``mirrored_weight`` each -- 1/2 for the solver's
    convention, 1 for the full-weight reading the first-order gate refutes.

    ``rows_of`` overrides :func:`_rows` for fixtures with a different
    declaration (the complex ``PairHop`` fixture below)."""
    rows_of = rows_of or _rows
    with ref._quiet_matmul():
        H = np.zeros((ref.DIM, ref.DIM), complex)
        for s in range(2):
            for m in range(NORB):
                for n in range(NORB):
                    if H0[m, n] != 0:
                        H = H + H0[m, n] * (CD[ref._so(m, s)] @ C[ref._so(n, s)])
        for itype, v in terms:
            w = mirrored_weight if itype in _MIRRORED_TYPES else 1.0
            for (rx, ry, rz, a1, b1, vr, vi) in rows_of(itype, v):
                H = H + w * _h_int(itype, a1 - 1, b1 - 1, complex(vr, vi))
        N = sum(CD[p] @ C[p] for p in range(ref.NSO))
    return H - MU * N


def _green_ed(H, iws):
    """Lehmann ``G_{mn}(i omega)`` of the up spin for the exact state."""
    E, V = np.linalg.eigh(H)
    E = E - E.min()
    w = np.exp(-BETA * E)
    Z = w.sum()
    with ref._quiet_matmul():
        cm = [V.conj().T @ C[ref._so(m, UP)] @ V for m in range(NORB)]
    boltz = w[:, None] + w[None, :]
    dE = E[None, :] - E[:, None]
    G = np.zeros((len(iws), NORB, NORB), complex)
    for m in range(NORB):
        for n in range(NORB):
            num = cm[m] * np.conj(cm[n]) * boltz
            for i, iw in enumerate(iws):
                G[i, m, n] = (num / (iw - dE)).sum()
    return G / Z


def _rho_ed(H):
    """``<c+_p c_q>`` of the exact state on the generalised index
    ``p = s * norb + a`` (the spin-block order of ``compile_onsite``)."""
    E, V = np.linalg.eigh(H)
    E = E - E.min()
    w = np.exp(-BETA * E)
    Z = w.sum()
    M = 2 * NORB
    rho = np.zeros((M, M), complex)
    with ref._quiet_matmul():
        for p in range(M):
            for q in range(M):
                sp, mp = divmod(p, NORB)
                sq, mq = divmod(q, NORB)
                op = V.conj().T @ (CD[ref._so(mp, sp)] @ C[ref._so(mq, sq)]) @ V
                rho[p, q] = (np.diag(op) * w).sum() / Z
    return rho


def _onsite_table(terms, rows_of=None):
    """The reader's on-site table ``{type: {((0,0,0),(a,b)): v}}`` of the
    same declared rows the ED Hamiltonian was built from."""
    rows_of = rows_of or _rows
    tbl = {}
    for itype, v in terms:
        tbl.setdefault(itype, {})
        for (rx, ry, rz, a1, b1, vr, vi) in rows_of(itype, v):
            tbl[itype][((0, 0, 0), (a1 - 1, b1 - 1))] = complex(vr, vi)
    return tbl


def _sigma_ed(terms, iws, mirrored_weight=0.5, hf_source="wick", rows_of=None):
    """``(Sigma_ED, Sigma_HF)`` on ``iws``: the exact self-energy from the
    Dyson inversion, and a first-order functional evaluated with the EXACT
    density (the up-spin block).

    ``hf_source`` selects WHOSE first order is subtracted:

    ``"wick"`` (the default, and what the second-order gate uses)
        the TEST SIDE's own :func:`wick_hf_sigma` on the ED monomials. The
        remainder whose O(v^2) coefficient the gate compares is then
        production-free: it is not built with any part of the object under
        test.
    ``"production"``
        ``hf_first_order(compile_onsite(...))``. Used only by
        :meth:`TestG3.test_first_order_remainder_vanishes`, where the point
        IS to confront the production compiler with this Hamiltonian -- a
        convention pin, and the reason that gate can be two-sided in the
        mirrored-row weight at all (the ``"wick"`` functional tracks
        whatever weight the Hamiltonian was built with, so it cannot
        discriminate one).
    """
    H = _hamiltonian(terms, mirrored_weight, rows_of)
    G = _green_ed(H, iws)
    G0 = ref._green0_w(iws)
    sigma = np.array([np.linalg.inv(G0[i]) - np.linalg.inv(G[i]) for i in range(len(iws))])
    if not terms:
        return sigma, np.zeros((NORB, NORB), complex)
    if hf_source == "production":
        from hwave.solver.second_order import compile_onsite, hf_first_order
        Gamma = compile_onsite(_onsite_table(terms, rows_of), NORB)
        return sigma, hf_first_order(Gamma, _rho_ed(H))[:NORB, :NORB]
    if hf_source != "wick":
        raise ValueError("hf_source must be 'wick' or 'production'")
    hf = wick_hf_sigma(2 * NORB, _mode_terms(terms, mirrored_weight, rows_of),
                       _gen_to_mode(_rho_ed(H)))
    return sigma, _mode_to_gen(hf)[:NORB, :NORB]


def _sigma_fluct_ed(terms, iws, mirrored_weight=0.5, hf_source="wick", rows_of=None):
    """``Sigma_ED - Sigma_HF[rho_ED]``; see :func:`_sigma_ed` for
    ``hf_source``."""
    sigma, hf = _sigma_ed(terms, iws, mirrored_weight, hf_source, rows_of)
    return sigma - hf[None]


# --------------------------------------------------------------------------
# The production side
# --------------------------------------------------------------------------

def _production_fluct(terms, nmat):
    """Production ``Sigma_fluct`` of ONE map: the general FLEX path with
    ``flex_second_order = "local"`` on the single site, at the FIXED
    ``mu = MU`` and on the BARE Green function (no SCF, no density update),
    so the propagator is the same ``G0`` the ED expansion is taken around."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp()
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(NORB)
                    + "".join("0 0 0\n" for _ in range(NORB)))
        rows = [(0, 0, 0, m + 1, n + 1, H0[m, n].real, H0[m, n].imag)
                for m in range(NORB) for n in range(NORB) if H0[m, n] != 0]
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", NORB, rows)
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for itype, v in terms:
            _write_wan(os.path.join(d, itype + ".dat"), itype, NORB, _rows(itype, v))
            idict[itype] = itype + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 1.0 / BETA, "mu": MU, "CellShape": [1, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
               "flex_second_order": "local"}
        s = flex_mod.FLEX(r.get_param("ham"), {},
                          {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
                           "calc_scheme": "general"})
        s._calc_epsilon_k({})
        G = s._calc_dressed_green(BETA, MU, np.zeros((1, nmat, 1, NORB, NORB), complex))
        chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), BETA)[0]
        _, v_eff, _, _ = s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q)
        return s._calc_self_energy_general(G, v_eff, BETA)[0, :, 0]
    finally:
        shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------------------------------
# Stencils
# --------------------------------------------------------------------------

def _coeff2(fn, x):
    """Second-order coefficient of one coupling: ``(f(2x) - 2 f(x)) / (2 x^2)``
    (exact for ``f = c1 v + c2 v^2``, with an O(x) error from ``c3``)."""
    return (fn(2 * x) - 2 * fn(x)) / (2 * x * x)


def _coeff11(fn, x, y):
    """Cross coefficient of two couplings. ``f(0, 0) = 0`` on both sides (no
    interaction, no fluctuation self-energy), so it is not evaluated."""
    return (fn(x, y) - fn(x, 0.0) - fn(0.0, y)) / (x * y)


def _refine(coeff, x):
    """Two-level refinement of a stencil whose leading error is O(x)."""
    return 2.0 * coeff(x / 2) - coeff(x)


def _guard_window(nmat):
    """A FIXED absolute Matsubara window, as a slice into an ``nmat`` grid.

    The comparison window of the gate itself is the inner EIGHTH of its own
    grid, so it covers twice the frequency range at ``nmat = 4096`` as at
    2048; a cutoff comparison has to hold the frequencies fixed instead."""
    return slice(nmat // 2 - _GUARD_WINDOW, nmat // 2 + _GUARD_WINDOW)


class _Maps(object):
    """Memoised ED / production maps -- the pure-type points of the pair
    stencils are the same maps the pure-type stencils already asked for."""

    def __init__(self, nmat):
        self.nmat = nmat
        self.iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA
        self._ed, self._pr = {}, {}
        self.runs = 0

    @staticmethod
    def _key(terms):
        return tuple(sorted((t, float(v)) for t, v in terms if v != 0.0))

    def ed(self, terms):
        k = self._key(terms)
        if not k:
            return np.zeros((self.nmat, NORB, NORB), complex)
        if k not in self._ed:
            self._ed[k] = _sigma_fluct_ed(list(k), self.iws)
        return self._ed[k]

    def prod(self, terms):
        k = self._key(terms)
        if not k:
            return np.zeros((self.nmat, NORB, NORB), complex)
        if k not in self._pr:
            self._pr[k] = _production_fluct(list(k), self.nmat)
            self.runs += 1
        return self._pr[k]


#: Phase of the complex ``PairHop`` fixture: a row value ``h * _PH_PHASE``
#: with its Hermitian partner ``conj(h * _PH_PHASE)`` on the transposed row.
_PH_PHASE = 1.0 + 0.5j


def _rows_complex_pairhop(itype, v):
    """The complex ``PairHop`` declaration: ``(1, 2, z)`` and its Hermitian
    partner ``(2, 1, conj z)``, ``z = v * _PH_PHASE``.

    ``PairHop``'s transposed row is the REVERSE hop, i.e. the Hermitian
    conjugate of the row's operator, so the closed table keeps the phase
    (``v P + v^* P^dagger``) -- unlike the same-operator types, whose
    ordered-row sum folds a complex coefficient to its real part. This
    fixture is the one that can see the difference."""
    if itype != "PairHop":
        raise ValueError("the complex fixture declares PairHop only")
    z = complex(v) * _PH_PHASE
    return [(0, 0, 0, 1, 2, z.real, z.imag), (0, 0, 0, 2, 1, z.real, -z.imag)]


def _random_densities(n=4, spin_block_diagonal=True, seed=20260915):
    """``n`` deterministic random Hermitian densities on the generalised
    index ``s * NORB + a``, centred on half filling.

    With ``spin_block_diagonal`` the spin off-diagonal blocks are zeroed
    (the domain spec 2.2's contract (i) states); without it the densities
    are generic Hermitian ones, which is the only way the on-site
    ``PairLift`` first order is nonzero at all."""
    rng = np.random.default_rng(seed)
    M = 2 * NORB
    out = []
    for _ in range(n):
        A = rng.normal(size=(M, M)) + 1j * rng.normal(size=(M, M))
        rho = 0.05 * (A + A.conj().T) + 0.5 * np.eye(M)
        if spin_block_diagonal:
            rho[:NORB, NORB:] = 0.0
            rho[NORB:, :NORB] = 0.0
        out.append(rho)
    return out


class TestWickFunctional(unittest.TestCase):
    """The test side's own operators and first-order functional: what they
    are, and that the production compiler agrees with them.

    The second-order gate subtracts :func:`wick_hf_sigma`, not the
    production functional, so the remainder it compares is independent of
    the object under test. The agreement between the two is asserted HERE,
    as a convention pin."""

    def test_term_list_reproduces_the_ed_operators(self):
        """``_h_int_terms`` multiplied out IS ``_h_int`` -- for every type,
        at a complex coupling, exactly."""
        with ref._quiet_matmul():
            for itype in _TYPES:
                a, b = (0, 0) if itype == "CoulombIntra" else (0, 1)
                v = 0.37 - 0.21j
                dense = np.zeros((ref.DIM, ref.DIM), complex)
                for (p, q, r, t, c_) in _h_int_terms(itype, a, b, v):
                    dense = dense + c_ * (CD[p] @ C[q] @ CD[r] @ C[t])
                with self.subTest(itype=itype):
                    self.assertGreater(np.abs(dense).max(), 1e-3)     # anti-vacuity
                    np.testing.assert_array_equal(dense, _h_int(itype, a, b, v))

    def test_independent_functional_equals_the_production_first_order(self):
        """``wick_hf_sigma`` on the ED monomials equals
        ``hf_first_order(compile_onsite(...))`` on SEVERAL deterministic
        random Hermitian densities -- spin-block-diagonal ones (the domain
        of spec 2.2's contract (i)) and generic ones, which is the only
        family on which the on-site ``PairLift`` first order is nonzero."""
        from hwave.solver.second_order import compile_onsite, hf_first_order
        for block in (True, False):
            for k, rho_gen in enumerate(_random_densities(spin_block_diagonal=block)):
                for itype in _TYPES:
                    terms = [(itype, 0.41)]
                    mine = _mode_to_gen(wick_hf_sigma(2 * NORB, _mode_terms(terms),
                                                      _gen_to_mode(rho_gen)))
                    ref_hf = hf_first_order(compile_onsite(_onsite_table(terms), NORB), rho_gen)
                    with self.subTest(itype=itype, block_diagonal=block, density=k):
                        if not (block and itype in _NO_FIRST_ORDER):
                            self.assertGreater(np.abs(ref_hf).max(), 1e-4)    # anti-vacuity
                        self.assertLess(np.abs(mine - ref_hf).max(),
                                        1e-12 * max(np.abs(ref_hf).max(), 1e-3))

    def test_wick_remainder_has_no_first_order_at_either_weight(self):
        """The remainder the second-order gate subtracts really does start
        at O(v^2) -- and, unlike the production subtraction, it does so for
        EITHER mirrored-row weight, because this functional is the first
        order of whichever Hamiltonian it was handed. That is why the
        weight is pinned by
        :meth:`TestG3.test_first_order_remainder_vanishes` (which uses the
        production functional) and not here."""
        nmat, x = 32, 0.05
        iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA
        for itype in _TYPES:
            for weight in (0.5, 1.0):
                lin = _refine(lambda h: _sigma_fluct_ed([(itype, h)], iws, weight) / h, x)
                sig, _ = _sigma_ed([(itype, x)], iws, weight)
                scale = np.abs(sig).max() / x
                # the stencil's leftover truncation is the O(v^2) term, which
                # the full-weight Hamiltonian carries four times over
                bound = _TOL * (weight / 0.5) ** 2
                with self.subTest(itype=itype, mirrored_weight=weight):
                    self.assertGreater(scale, 1e-3)                   # anti-vacuity
                    self.assertLess(np.abs(lin).max() / scale, bound)

    def test_complex_pairhop_conjugation_is_pinned(self):
        """A COMPLEX ``PairHop`` declaration on the single site: the
        production compiler keeps the phase of the Hermitian partner, and
        the reading that does not is refuted.

        Two-sided. (a) With the Hermitian-partner declaration
        (:func:`_rows_complex_pairhop`) the ED Hamiltonian is Hermitian and
        the production first-order functional leaves no O(v) remainder.
        (b) The non-conjugated declaration -- the same row value on the
        transposed row -- compiles to a DIFFERENT first order at the same
        density (the reversal closure folds it to ``Re z``), so the
        conjugation is measured, not assumed."""
        from hwave.solver.second_order import compile_onsite, hf_first_order
        nmat, x = 32, 0.05
        iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA
        H = _hamiltonian([("PairHop", x)], rows_of=_rows_complex_pairhop)
        self.assertLess(np.abs(H - H.conj().T).max(), 1e-14)          # Hermitian fixture
        self.assertGreater(np.abs(_PH_PHASE.imag), 0.1)               # the phase is real content

        lin = _refine(lambda h: _sigma_fluct_ed(
            [("PairHop", h)], iws, 0.5, "production", _rows_complex_pairhop) / h, x)
        sig, _ = _sigma_ed([("PairHop", x)], iws, 0.5, "production", _rows_complex_pairhop)
        scale = np.abs(sig).max() / x
        self.assertGreater(scale, 1e-3)                               # anti-vacuity
        self.assertLess(np.abs(lin).max() / scale, _TOL,
                        "the complex PairHop first order is not subtracted by the compiled "
                        "Gamma: the Hermitian partner's conjugation is wrong")

        rho = _random_densities(n=1, spin_block_diagonal=False)[0]
        good = _onsite_table([("PairHop", x)], _rows_complex_pairhop)
        z = good["PairHop"][((0, 0, 0), (0, 1))]
        # the REAL-PART reading: what the compiler would produce if it
        # treated PairHop's transposed row as the same operator (the
        # mirrored-row fold) instead of the Hermitian partner. It is a legal
        # Hermitian-closed table, so the guards below do not intercept it.
        folded = {"PairHop": {((0, 0, 0), (0, 1)): complex(z.real),
                              ((0, 0, 0), (1, 0)): complex(z.real)}}
        s_good = hf_first_order(compile_onsite(good, NORB), rho)
        s_folded = hf_first_order(compile_onsite(folded, NORB), rho)
        self.assertGreater(np.abs(s_good).max(), 1e-4)                # anti-vacuity
        self.assertGreater(np.abs(s_good - s_folded).max(), 0.1 * np.abs(s_good).max(),
                           "the real-part reading compiles to the same first order: the "
                           "phase of a complex PairHop carries no content here, so the "
                           "gate above cannot pin the conjugation")
        # and the table that is NOT Hermitian-closed -- the same value on the
        # transposed row -- is refused by name
        unclosed = {"PairHop": {((0, 0, 0), (0, 1)): z, ((0, 0, 0), (1, 0)): z}}
        with self.assertRaises(ValueError) as cm:
            compile_onsite(unclosed, NORB)
        self.assertIn("PairHop", str(cm.exception))
        self.assertIn("Hermitian", str(cm.exception))


class TestG3(unittest.TestCase):
    """Exact diagonalisation of the on-site second order (design 2026-09-08,
    gate G3)."""

    def test_first_order_remainder_vanishes(self):
        """The HF subtraction is exact at first order -- and the mirrored-row
        HALF weight is what makes it so.

        This gate is what licenses the second-order comparison to use the
        production compiler's ``Gamma`` for the subtraction: if ``Gamma``
        were not the first-order functional of THIS Hamiltonian, the O(v)
        remainder would not vanish. It is two-sided: the full-weight
        Hamiltonian (each declared row of a mirrored type at full ``v``)
        must FAIL it, which is what turns the weight from an assumption
        into a measurement.

        Measured (O(v) remainder over the self-energy's own linear size,
        half weight / full weight): CoulombIntra 8.8e-5 (no mirrored row),
        CoulombInter 1.3e-4 / 5.1e-1, Hund 6.1e-5 / 4.9e-1, Ising 2.4e-4 /
        4.9e-1, Exchange 1.8e-4 / 5.1e-1, PairHop 1.9e-4 (Hermitian
        partner, not mirrored), PairLift 7.0e-4 (no first order at all)."""
        nmat = 32
        iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / BETA
        x = 0.05
        for itype in _TYPES:
            for weight, ok in ((0.5, True), (1.0, False)):
                discriminates = itype in _MIRRORED_TYPES and itype not in _NO_FIRST_ORDER
                if weight == 1.0 and not discriminates:
                    # either the two conventions build the same Hamiltonian
                    # (CoulombIntra, PairHop), or the first order vanishes in
                    # both of them (PairLift): nothing to discriminate.
                    continue
                lin = _refine(lambda h: _sigma_fluct_ed([(itype, h)], iws, weight,
                                                       "production") / h, x)
                sig, _ = _sigma_ed([(itype, x)], iws, weight, "production")
                scale = np.abs(sig).max() / x          # the self-energy's own linear size
                with self.subTest(itype=itype, mirrored_weight=weight):
                    self.assertGreater(scale, 1e-3)    # anti-vacuity
                    rel = np.abs(lin).max() / scale
                    if ok:
                        self.assertLess(rel, _TOL,
                                        "{}: the first-order remainder does not vanish with the "
                                        "half weight ({:.2e} of the self-energy)".format(itype, rel))
                    else:
                        self.assertGreater(rel, 0.1,
                                           "{}: the FULL-weight Hamiltonian also passes the "
                                           "first-order gate -- the weight is not pinned"
                                           .format(itype))

    def test_pairlift_first_order_is_zero(self):
        """On-site PairLift has no first order at all on a paramagnetic,
        spin-block-diagonal density: its monomials move a pair of spins,
        so every Hartree and Fock contraction needs a spin-off-diagonal
        ``<c+_{up} c_{dn}>``. Recorded because it is why the first-order
        gate above cannot pin PairLift's weight (both conventions give
        zero) -- only the second order does."""
        from hwave.solver.second_order import compile_onsite, hf_first_order
        tbl = {"PairLift": {((0, 0, 0), (0, 1)): 0.3, ((0, 0, 0), (1, 0)): 0.3}}
        Gamma = compile_onsite(tbl, NORB)
        self.assertGreater(np.abs(Gamma).max(), 1e-3)            # anti-vacuity
        rho = _rho_ed(_hamiltonian([]))
        self.assertLess(np.abs(hf_first_order(Gamma, rho)).max(), 1e-14)

    def test_extraction_is_converged(self):
        """The working point of :meth:`test_every_type_and_pair` is
        converged in BOTH of its discretisations -- the Matsubara cutoff and
        the stencil step -- so the fixed inner window really can stand in
        for an extrapolation.

        The reference test this module's stencils come from
        (``tests/test_flex_second_order_sopt.py``) extrapolates the STEP
        away with a three-level Richardson ladder and leaves the cutoff
        alone; here the cutoff is the larger of the two error sources (the
        production kernel's finite frequency grid is a fixed fraction of the
        coefficient, and it does not shrink with the step), so a ladder in
        the step alone would extrapolate towards the wrong limit. What
        replaces it is this executable guard: the coefficient is extracted
        at two cutoffs and three steps and the drifts are bounded directly.

        Three statements, all on the same absolute frequency window
        (:func:`_guard_window`) so the two cutoffs are comparable:

        * CUTOFF: the refined coefficient at ``nmat = 2048`` and at 4096
          differ by less than :data:`_GUARD_CUTOFF_BOUND` of the
          coefficient;
        * STEP: the refined coefficient at ``_X``, ``_X/2`` and ``_X/4``
          drifts by less than :data:`_GUARD_STEP_BOUND`;
        * ORDER: after the two-level refinement the step error is O(x^2),
          i.e. halving the step shrinks the drift by roughly four.

        ``CoulombIntra`` alone is used -- it is the largest coefficient of
        the table and the one every mixed pair is built on, and the
        discretisation errors are properties of the extraction, not of the
        interaction type. Sixteen maps in all, 3.0 s: below the ~5 s the
        heavy registry is for, so this guard runs in the FAST gate even
        though the comparison it underwrites does not.
        """
        itype = "CoulombIntra"
        with _quiet():
            coeffs = {}
            for nmat in (2048, _NMAT):
                maps = _Maps(nmat)
                win = _guard_window(nmat)
                for h in (_X, _X / 2, _X / 4):
                    ed = _refine(lambda g: _coeff2(lambda v: maps.ed([(itype, v)]), g), h)
                    pr = _refine(lambda g: _coeff2(lambda v: maps.prod([(itype, v)]), g), h)
                    coeffs[(nmat, h)] = (ed[win], pr[win])

        scale = np.abs(coeffs[(_NMAT, _X)][0]).max()
        self.assertGreater(scale, _FLOOR)                            # anti-vacuity
        for k, (ed, pr) in coeffs.items():
            self.assertTrue(np.all(np.isfinite(ed)) and np.all(np.isfinite(pr)),
                            "non-finite coefficient at {}".format(k))

        cutoff = np.abs(coeffs[(_NMAT, _X)][1] - coeffs[(2048, _X)][1]).max() / scale
        self.assertLess(cutoff, _GUARD_CUTOFF_BOUND,
                        "the production coefficient still moves with the Matsubara cutoff "
                        "({:.2e} of the coefficient between nmat 2048 and {})"
                        .format(cutoff, _NMAT))

        d1 = np.abs(coeffs[(_NMAT, _X)][1] - coeffs[(_NMAT, _X / 2)][1]).max() / scale
        d2 = np.abs(coeffs[(_NMAT, _X / 2)][1] - coeffs[(_NMAT, _X / 4)][1]).max() / scale
        self.assertLess(d1, _GUARD_STEP_BOUND, "the step drift is not inside the bound")
        self.assertLess(d2, _GUARD_STEP_BOUND, "the step drift is not inside the bound")
        self.assertGreater(d1, 0.0)                                  # anti-vacuity
        self.assertGreater(d1 / d2, 2.0,
                           "the post-refinement step error does not fall as O(x^2) "
                           "(drifts {:.2e} then {:.2e})".format(d1, d2))
        self.assertLess(d1 / d2, 8.0,
                        "the post-refinement step error falls faster than O(x^2): the "
                        "refinement order recorded in the module docstring is wrong "
                        "(drifts {:.2e} then {:.2e})".format(d1, d2))

    @heavy
    def test_every_type_and_pair(self):
        """G3: the O(x^2) coefficient of every type and the O(x y) coefficient
        of every pair of the covering set agree between exact diagonalisation
        and the production ``Sigma_fluct``."""
        maps = _Maps(_NMAT)
        inner = slice(_NMAT // 2 - _WINDOW, _NMAT // 2 + _WINDOW)
        with _quiet():
            for itype in _TYPES:
                with self.subTest(itype=itype):
                    ed = _refine(lambda h: _coeff2(lambda v: maps.ed([(itype, v)]), h), _X)
                    pr = _refine(lambda h: _coeff2(lambda v: maps.prod([(itype, v)]), h), _X)
                    scale = np.abs(ed[inner]).max()
                    self.assertGreater(scale, _FLOOR)            # anti-vacuity
                    # max-norm RELATIVE error against the coefficient's own
                    # size (an rtol+atol pair would let a large entry hide a
                    # proportionally large error in a small one)
                    self.assertLess(np.abs(pr[inner] - ed[inner]).max() / scale, _TOL,
                                    "second order of {}".format(itype))
            for tx, ty in _PAIRS:
                with self.subTest(pair=(tx, ty)):
                    ed = _refine(lambda h: _coeff11(
                        lambda a, b: maps.ed([(tx, a), (ty, b)]), h, h), _X)
                    pr = _refine(lambda h: _coeff11(
                        lambda a, b: maps.prod([(tx, a), (ty, b)]), h, h), _X)
                    scale = np.abs(ed[inner]).max()
                    if (tx, ty) in _KNOWN_ZERO:
                        # a relative comparison against a zero reference is
                        # vacuous: check the TWO-SIDED zero instead.
                        self.assertLess(scale, _ZERO_CEIL,
                                        "{} x {}: listed as a spin-algebra zero".format(tx, ty))
                        self.assertLess(np.abs(pr[inner]).max(), _ZERO_CEIL,
                                        "{} x {}: the ED cross term is zero but production's "
                                        "is not".format(tx, ty))
                        continue
                    self.assertGreater(scale, _FLOOR)            # anti-vacuity
                    self.assertLess(np.abs(pr[inner] - ed[inner]).max() / scale, _TOL,
                                    "cross term of {} x {}".format(tx, ty))


if __name__ == "__main__":
    unittest.main()
