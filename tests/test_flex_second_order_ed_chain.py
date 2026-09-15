#!/usr/bin/env python3

"""G4: exact diagonalisation of the OFF-SITE second order on a chain.

Task 10 (``tests/test_flex_second_order_ed_onsite.py``, gate G3) settled the
ON-SITE second order on a single site, where every interaction entry is
local and the locality split of ``flex_second_order = "local"`` has nothing
to drop. This module is the off-site half: on periodic chains
(``CellShape = [L, 1, 1]``) with a NEIGHBOUR bond declared, the O(v^2)
coefficient of the exact fluctuation remainder

    Sigma_ED - Sigma_HF[rho_ED]

must equal the production ``Sigma_fluct`` of the general path PLUS the
off/off uncrossed exchange class that the local kernel deliberately drops
(``tests/test_second_order_oracle.py``'s ``"dropped"`` weighting); and with
the bond gate on -- which resums that class -- production alone must equal
the ED remainder.

Two fixtures, and what each adjudicates
--------------------------------------
* ``_fx_chain()`` -- ``L = 4``, one orbital, ``t = -1``, an on-site ``U``
  and a nearest-neighbour ``V``. Asserted: the O(V^2), O(U V) and O(U^2)
  coefficients of ``local + dropped`` and, for ``V``, the gate-on
  coefficient. Off-site Hund and Ising on the same chain are RECORDED
  (spec G2 (c) / G4), not asserted.
* ``_fx_orbital()`` -- ``L = 3``, TWO orbitals, and an ORBITAL-ASYMMETRIC
  inter-orbital bond (``v_01(+x) = v``, ``v_10(+x) = 0.6 v``). This is the
  fixture that adjudicates the two open ledger items of the campaign: the
  orientation of the off-site density vertex (Task 5) and the bond gate's
  deviation on inter-orbital bonds (Task 9). Its verdicts are in
  :class:`TestG4`'s two methods and in the module's measured table below.

Which Hamiltonian the chain ED must be
--------------------------------------
The declared row ``(R, a, b, v)`` of an off-site density type places
orbital ``a`` on the DISPLACED site and orbital ``b`` on the reference one,

    v n_{j+R, a} n_{j, b}

-- the same ``(j + R, a) <- (j, b)`` reading a ``transfer.dat`` row gets
(``EDFixture.build_h1``). This is NOT assumed here: :meth:`TestChain
Hamiltonian.test_ed_hamiltonian_is_the_production_mean_field` compares the
first-order functional of the ED term list against UHFk's OWN mean-field
kernel (``hwave.solver.hartree_fock.accumulate_hf``, Hartree AND Fock, at
the same density) and fails with the opposite (orbital-swapped) placement
by 5.4e-2 / 1.4e-2 of the mean field on the two inter-orbital fixtures,
while the placement above matches at 1e-16. On the single-orbital chain the
two placements are the same Hamiltonian, which is why only the two-orbital
fixture can decide it.

SCOPE OF THAT CLAIM, and of every verdict below it. What is adjudicated
here is the second-order vertex against UHFk's MEAN-FIELD READING of a
declared row: ``accumulate_hf`` is what defines the orbital placement the
chain ED implements, and the gate then asks whether the second order is
consistent with THAT Hamiltonian. It is not an adjudication against a
file-format definition -- the documentation never fixes the sign of
``r_ij`` in an interaction row, so "which site the row's first orbital sits
on" has no independent written answer to appeal to. Both readings are
internally consistent conventions; this module pins that the solver's
first-order and second-order paths agree on ONE of them, and says which.
The k-space bridge the comparison rides on (exponent sign, orbital index
order) is pinned separately, on a band that can see it, by
:meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`
-- without that pin a flipped bridge would silently select the opposite
placement and invert every verdict below.

The mirrored-row HALF weight (each of the two declared orientations of a
bond carries half of it, ``_MIRRORED_TYPES``) is the convention Task 10
established on the on-site types; here it is re-pinned two-sidedly by
:meth:`TestFirstOrder.test_first_order_remainder_vanishes` -- vanishing at
the half weight for every fixture, and 0.48 ... 0.53 of the self-energy's
own linear size at the full weight.

Stencils, working point, tolerance
----------------------------------
The coefficient extraction is Task 10's, imported rather than recopied:
the two-point stencils ``_coeff2`` / ``_coeff11`` with the O(x) refinement
``2 c(x/2) - c(x)``. Working point: ``nmat = 4096``, ``x = 3.125e-3``, and
the innermost eighth of the Matsubara grid -- Task 10's window, every one
of its 512 frequencies on both sides. The two error sources are the same
ones Task 10 measured: the finite ``x`` truncation (which the refinement
leaves at O(x^2)) and the production kernel's finite Matsubara window
(which is a fixed fraction of the coefficient and shrinks as ``1/nmat``:
9.2e-4 at nmat = 1024, 4.3e-4 at 2048, 2.3e-4 at 4096 on the V^2 entry).

Measured at that working point (relative to the ED coefficient's own
maximum over the window; tolerance 2e-3):

    fixture     entry              local + dropped     gate on
    L=4 norb=1  V^2                2.3e-4              2.4e-4
    L=4 norb=1  U V                4.8e-4              --
    L=4 norb=1  U^2                3.2e-4              --
    L=4 norb=1  Hund   (recorded)  2.5e-4              2.5e-4
    L=4 norb=1  Ising  (recorded)  2.2e-4              2.2e-4
    L=3 norb=2  v^2 (asymmetric)   3.5e-4              4.0e-2  <-- recorded

The last row is the campaign's open item: the STANDALONE local path plus
the dropped class reproduces the exact inter-orbital off-site second order
(so the Task 5 orientation of the off-site density vertex is confirmed by
exact diagonalisation), while the BOND GATE misses it by 4.0e-2 -- twenty
times the tolerance, and nearly twice the whole off/off uncrossed class it
is supposed to be resumming (2.2e-2 of the coefficient). That deviation is
RECORDED by a print and routed to a Phase B follow-up; it is not asserted
here, and it is not a second-order kernel finding.
"""
import os
import tempfile
import unittest

import numpy as np

from tests import ed_oracle_util as edu
from tests.heavy_tests import heavy
from tests.test_flex_second_order_ed_onsite import _coeff2, _coeff11, _quiet, _refine
from tests.test_second_order_factors import _write_wan
from tests.test_second_order_oracle import (_MIRRORED_ROW_TYPES as _MIRRORED_TYPES,
                                             oracle_records, oracle_sigma2)

UP, DN = 0, 1

# ``_MIRRORED_TYPES`` above is the oracle's own list, imported rather than
# recopied: the types whose two declared orientations name the SAME bond, so
# each declared row is half of it. Pinned two-sidedly here by
# :meth:`TestFirstOrder.test_first_order_remainder_vanishes`.

#: Working point of the heavy comparisons (see the module docstring).
_NMAT = 4096
_X = 3.125e-3
_WINDOW = _NMAT // 16          # half-width of the comparison window
_TOL = 2.0e-3

#: Anti-vacuity floor on the ED coefficient of a compared entry: a relative
#: comparison against a numerically zero reference would pass for any
#: production value at all. Every entry compared here measures
#: 5.8e-2 ... 3.6e-1.
_FLOOR = 1.0e-6

#: How large the dropped (off/off uncrossed) class must be, relative to the
#: coefficient, for the ``local + dropped`` assertion to be distinguishing
#: at all. Measured: 1.2e-1 for V^2 and 2.3e-1 for Hund on the
#: single-orbital chain, 2.2e-2 on the inter-orbital one.
_DROPPED_FLOOR = 1.0e-3


def _inner():
    """Indices of the comparison window on the production Matsubara grid --
    the innermost eighth, every frequency of it (Task 10's window)."""
    return np.arange(_NMAT // 2 - _WINDOW, _NMAT // 2 + _WINDOW)


def _fx_chain():
    """``L = 4``, one orbital, ``t = -1``: the brief's single-orbital chain."""
    return edu.EDFixture(L=4, norb=1, t={(0, 0): -1.0}, eps=(0.0,), T=0.5, mu=0.2)


def _fx_orbital():
    """``L = 3``, two orbitals, asymmetric intra/inter-orbital hopping and
    distinct on-site energies: the fixture that can see an orbital index
    order (``tests/test_bond_vs_ed_oracle.py``'s case-M shape, at the same
    ``T`` and ``mu``). ``nmode = 12``, Fock dimension 4096, largest
    ``(N_up, N_dn)`` block 400 -- comfortably inside ``SectorED``."""
    return edu.EDFixture(L=3, norb=2,
                         t={(0, 0): -1.0, (1, 1): -0.7, (0, 1): -0.2, (1, 0): -0.2},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


def _fx_complex():
    """``L = 3``, two orbitals, COMPLEX hopping in every orbital channel and
    ``t_01 != conj(t_10)``: a band that is neither inversion-symmetric nor
    transpose-symmetric.

    The two fixtures the gate itself runs on both have real hopping, so
    ``eps(k) = eps(-k)`` and ``eps(k) = eps(k)^T`` hold to round-off there
    and NEITHER the Fourier sign nor the orbital index order of
    :func:`_to_k` can be seen. This fixture exists only to see them
    (:meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`);
    it carries no interaction and is never diagonalised beyond the free
    Green function.

    ``EDFixture.build_h1`` places ``t[(a, b)]`` on ``(j+1, a) <- (j, b)``
    and ``conj(t[(b, a)])`` on the return hop, so ``h1`` is Hermitian for
    ANY pair ``t_01``, ``t_10`` -- the asymmetry below is legal, not a
    broken Hamiltonian, and ``_write_chain_inputs`` writes exactly those
    two amplitudes into ``transfer.dat``."""
    return edu.EDFixture(L=3, norb=2,
                         t={(0, 0): -1.0 + 0.3j, (1, 1): -0.7 - 0.2j,
                            (0, 1): -0.2 + 0.15j, (1, 0): -0.25 - 0.1j},
                         eps=(0.0, 0.3), T=0.5, mu=0.2)


# --------------------------------------------------------------------------
# Declared rows -> the exact Hamiltonian
# --------------------------------------------------------------------------

def _rows_u_v(u, v):
    """On-site ``U`` and nearest-neighbour ``V`` of the single-orbital chain,
    as the reader delivers them (both orientations of the bond)."""
    rows = {}
    if u != 0.0:
        rows["CoulombIntra"] = [(0, 0, 0, 1, 1, u, 0.0)]
    if v != 0.0:
        rows["CoulombInter"] = [(1, 0, 0, 1, 1, v, 0.0), (-1, 0, 0, 1, 1, v, 0.0)]
    return rows


def _rows_offsite(itype, v):
    """One off-site single-orbital bond of ``itype``, both orientations."""
    return {itype: [(1, 0, 0, 1, 1, v, 0.0), (-1, 0, 0, 1, 1, v, 0.0)]}


def _rows_interorbital(v, ratio=0.6):
    """The ORBITAL-ASYMMETRIC inter-orbital bond of the two-orbital chain:
    ``v_01(+x) = v`` and ``v_10(+x) = ratio * v`` with ``ratio != 1``, each
    with its reversed partner. Both orbital assignments of the SAME
    displacement are declared with DIFFERENT amplitudes, so nothing here is
    invariant under the orbital swap the ED adjudicates."""
    return {"CoulombInter": [(1, 0, 0, 1, 2, v, 0.0), (-1, 0, 0, 2, 1, v, 0.0),
                             (1, 0, 0, 2, 1, ratio * v, 0.0), (-1, 0, 0, 1, 2, ratio * v, 0.0)]}


def _rows_interorbital_one(v):
    """ONE inter-orbital bond of the two-orbital chain, both orientations:
    ``v_01(+x) = v`` alone. The orbital swap maps this declaration onto a
    DIFFERENT Hamiltonian (unlike :func:`_rows_interorbital` at ratio 1,
    where the two declared bonds exchange and the swap is a symmetry), so
    this is the shape the orientation check needs."""
    return {"CoulombInter": [(1, 0, 0, 1, 2, v, 0.0), (-1, 0, 0, 2, 1, v, 0.0)]}


def _ed_terms(fx, rows_by_type, mirrored_weight=0.5, swap_orbitals=False):
    """The quartic term list (``(p, q, r, s, coeff)`` for
    ``coeff c^dag_p c_q c^dag_r c_s``, the form ``SectorED`` consumes) of
    the Hamiltonian the DECLARED rows describe.

    A row ``(R, a, b, v)`` of a density type contributes
    ``w v sum_j sum_{s1 s2} sigma(s1, s2) n_{j+R, a, s1} n_{j, b, s2}`` --
    orbital ``a`` on the DISPLACED site (see the module docstring; pinned
    against UHFk's mean-field kernel, not assumed) -- with ``w`` the
    mirrored-row weight and ``sigma`` the type's spin structure: all four
    pairings for ``CoulombInter``, only the same-spin ones with a minus for
    ``Hund``, and the density-DIFFERENCE signs for ``Ising``.
    ``CoulombIntra`` is the on-site ``v n_{j,a,up} n_{j,a,dn}``.

    ``swap_orbitals`` builds the OPPOSITE placement (``a`` on the reference
    site) and exists only so the orientation can be checked two-sidedly.
    """
    terms = []
    for itype, rows in rows_by_type.items():
        w = mirrored_weight if itype in _MIRRORED_TYPES else 1.0
        for (rx, ry, rz, a1, b1, vr, vi) in rows:
            v = w * complex(vr, vi)
            if v == 0.0:
                continue
            a, b, R = a1 - 1, b1 - 1, rx
            if swap_orbitals:
                a, b = b, a
            if itype == "CoulombIntra":
                if (rx, ry, rz) != (0, 0, 0) or a != b:
                    raise ValueError("CoulombIntra must be an on-site same-orbital row")
                for j in range(fx.L):
                    terms.append((fx.mode(j, a, UP), fx.mode(j, a, UP),
                                  fx.mode(j, a, DN), fx.mode(j, a, DN), v))
                continue
            if itype == "CoulombInter":
                def sigma(s1, s2, v=v):
                    return v
            elif itype == "Hund":
                def sigma(s1, s2, v=v):
                    return -v if s1 == s2 else 0.0
            elif itype == "Ising":
                def sigma(s1, s2, v=v):
                    return v if s1 == s2 else -v
            else:
                raise ValueError("chain ED: unsupported interaction type {!r}".format(itype))
            for j in range(fx.L):
                jp = (j + R) % fx.L
                for s1 in range(2):
                    for s2 in range(2):
                        c = sigma(s1, s2)
                        if c == 0.0:
                            continue
                        terms.append((fx.mode(jp, a, s1), fx.mode(jp, a, s1),
                                      fx.mode(j, b, s2), fx.mode(j, b, s2), c))
    return terms


def _hf_sigma(fx, terms, rho):
    """First-order (Hartree-Fock) self-energy of ``terms`` at the density
    ``rho[p, q] = <c^dag_p c_q>``, in mode space.

    The Wick engine is ``ed_oracle_util.hf_h1_from_terms``'s ``add()``
    verbatim (``E = c [rho_pq rho_rs + rho_ps (delta_qr - rho_rq)]`` and
    ``H_MF = sum dE/d rho_xy c^dag_x c_y``), evaluated at an ARBITRARY
    density rather than the free one -- which is the only difference, and
    the reason it is written here: the gates need the functional at the
    EXACT density. :meth:`TestChainHamiltonian
    .test_ed_hamiltonian_is_the_production_mean_field` checks the result
    against UHFk's own kernel, so this is not a second private derivation
    of the mean field but a bridge to the production one."""
    S = np.zeros((fx.nmode, fx.nmode), dtype=complex)

    def add(p, q, r, s, c_):
        S[p, q] += c_ * rho[r, s]
        S[r, s] += c_ * rho[p, q]
        S[p, s] += c_ * ((1.0 if q == r else 0.0) - rho[r, q])
        S[r, q] += -c_ * rho[p, s]

    for (p, q, r, s, coeff) in terms:
        add(p, q, r, s, coeff)
    return 0.5 * (S + S.conj().T)


def _production_hf_sigma_k(fx, rows_by_type, rho):
    """UHFk's OWN mean-field Hamiltonian for the declared rows at the
    density ``rho``, as ``Sigma_HF(k)[(s, a), (t, b)]``.

    ``hwave.solver.hartree_fock.accumulate_hf`` is ``UHFk._make_ham``'s
    interaction loop as a pure function -- the definition of the
    Hamiltonian the solver implements, Hartree AND Fock, for on-site and
    off-site rows alike. It is what makes the first-order gate below
    two-sided in both the mirrored-row weight and the orbital placement,
    exactly as Task 10's gate is two-sided through ``compile_onsite``."""
    from hwave.solver.hartree_fock import accumulate_hf, build_interaction_tables
    norb, L = fx.norb, fx.L
    shape = (L, 1, 1)
    param_ham = {}
    for itype, rows in rows_by_type.items():
        tbl = {}
        for (rx, ry, rz, a1, b1, vr, vi) in rows:
            key = ((rx % L, ry, rz), (a1 - 1, b1 - 1))
            tbl[key] = tbl.get(key, 0.0) + complex(vr, vi)
        param_ham[itype] = tbl
    tabs = build_interaction_tables(param_ham, norb, shape)
    rho_so_r = np.zeros((L, 2, norb, 2, norb), dtype=complex)
    for r in range(L):
        for s in range(2):
            for a in range(norb):
                for t in range(2):
                    for b in range(norb):
                        rho_so_r[r, s, a, t, b] = rho[fx.mode(0, a, s), fx.mode(r, b, t)]
    out = np.zeros((L, 2 * norb, 2 * norb), dtype=complex)
    accumulate_hf(out, rho_so_r, tabs.inter_table, tabs.spin_table, shape,
                  include_fock=True)
    return out


# --------------------------------------------------------------------------
# Mode space -> the solver's (k, orbital) frame
# --------------------------------------------------------------------------

def _k_phases(fx, sign=1.0):
    """``phase[kx, R] = e^{sign i k R}`` on ``k = 2 pi kx / L``. The ONE
    place the Fourier sign of this module's real-space-to-k bridge is
    written; ``sign = -1`` builds the reversed-k candidate that
    :meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`
    has to reject."""
    R = np.arange(fx.L)
    return np.exp(sign * 2j * np.pi * np.outer(R, R) / fx.L)


def _to_k(fx, M, sign=1.0):
    """The UP-spin block of a translation-invariant mode-space matrix, as
    ``A(k)_{ab} = sum_R A[(R, a), (0, b)] e^{i k R}``, ``k = 2 pi kx / L``.

    This is the convention ``EDFixture.build_h1`` and the solver's
    ``epsilon_k[a,b] += t_ab e^{+ikR}`` share (a ``transfer.dat`` row at
    ``R`` carries ``(j + R, a) <- (j, b)``), and
    :meth:`TestChainHamiltonian.test_bare_green_matches_the_production_solver`
    pins the whole chain of it end to end -- the SIGN of the exponent and
    the orbital index order specifically by
    :meth:`TestChainHamiltonian.test_fourier_sign_and_orbital_order_are_pinned`,
    which needs a band that is neither inversion- nor transpose-symmetric to
    see either."""
    norb, L = fx.norb, fx.L
    phase = _k_phases(fx, sign)
    out = np.zeros((L, norb, norb), dtype=complex)
    for kx in range(L):
        for R in range(L):
            for a in range(norb):
                for b in range(norb):
                    out[kx, a, b] += phase[kx, R] * M[fx.mode(R, a, UP), fx.mode(0, b, UP)]
    return out


def _g0_k(fx, iws, sign=1.0):
    """``G0(k, i w) = [(i w + mu) - eps(k)]^{-1}`` of the free chain."""
    ek = _to_k(fx, fx.build_h1(), sign)
    eye = np.eye(fx.norb)
    return np.array([[np.linalg.inv((iw + fx.mu) * eye - ek[kx]) for kx in range(fx.L)]
                     for iw in iws])


def _green_ed_k(ed, fx, iws, sign=1.0):
    """The interacting ``G(k, i w)`` of the up spin, from ``SectorED.green``.

    Only the ``norb`` columns at site 0 are asked for: the fixture is
    translation invariant, so ``G[(R, a), (0, b)]`` already carries every
    displacement, and the restriction turns an ``nmode^2`` Lehmann
    contraction into an ``nmode * norb`` one."""
    L, norb = fx.L, fx.norb
    rows = [fx.mode(j, a, UP) for j in range(L) for a in range(norb)]
    cols = [fx.mode(0, b, UP) for b in range(norb)]
    gm = ed.green(iws, rows=rows, cols=cols)
    phase = _k_phases(fx, sign)
    gk = np.zeros((len(iws), L, norb, norb), dtype=complex)
    for kx in range(L):
        for R in range(L):
            gk[:, kx] += phase[kx, R] * gm[:, R * norb:(R + 1) * norb, :]
    return gk


def _sigma_ed(fx, rows_by_type, iws, mirrored_weight=0.5):
    """``(Sigma_ED, Sigma_HF)`` on ``iws``: the exact self-energy of the
    chain from the Dyson inversion ``G0^{-1} - G^{-1}``, and UHFk's
    first-order functional evaluated with the EXACT density."""
    ed = edu.SectorED(fx, terms=_ed_terms(fx, rows_by_type, mirrored_weight))
    sigma = np.linalg.inv(_g0_k(fx, iws)) - np.linalg.inv(_green_ed_k(ed, fx, iws))
    hf = _production_hf_sigma_k(fx, rows_by_type, ed.density_matrix())
    return sigma, hf[:, :fx.norb, :fx.norb]


def _sigma_fluct_ed(fx, rows_by_type, iws, mirrored_weight=0.5):
    """``Sigma_ED - Sigma_HF[rho_ED]`` -- the exact fluctuation remainder,
    whose O(v^2) coefficient is the bare second-order skeleton."""
    sigma, hf = _sigma_ed(fx, rows_by_type, iws, mirrored_weight)
    return sigma - hf[None]


# --------------------------------------------------------------------------
# The production side
# --------------------------------------------------------------------------

def _write_chain_inputs(d, fx, rows_by_type):
    """``geom.dat`` / ``transfer.dat`` / one file per interaction type
    describing ``fx``'s chain on ``CellShape = [L, 1, 1]``: a one-cell
    geometry with ``norb`` orbitals, the hopping as ``(+1, a, b, t_ab)``
    and its Hermitian partner ``(-1, a, b, conj(t_ba))``, and the on-site
    energies at ``R = 0``."""
    norb = fx.norb
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n{}\n".format(norb))
        for _ in range(norb):
            f.write("0.0 0.0 0.0\n")
    rows = []
    for (a, b), tv in fx.t.items():
        back = np.conj(fx.t[(b, a)])
        rows.append((1, 0, 0, a + 1, b + 1, tv.real, tv.imag))
        rows.append((-1, 0, 0, a + 1, b + 1, back.real, back.imag))
    for o, e in enumerate(fx.eps):
        rows.append((0, 0, 0, o + 1, o + 1, float(np.real(e)), 0.0))
    _write_wan(os.path.join(d, "transfer.dat"), "Transfer", norb, rows)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for itype, irows in rows_by_type.items():
        _write_wan(os.path.join(d, itype.lower() + ".dat"), itype, norb, irows)
        idict[itype] = itype.lower() + ".dat"
    return idict


def _chain_solver(d, fx, rows_by_type, gate, nmat):
    """A FLEX general solver on ``fx``'s chain with ``flex_second_order =
    "local"``, at the FIXED ``mu = fx.mu`` (no SCF, no density update), and
    with the bond gate switched on through its own two parameters when
    ``gate`` -- so the Phase B preflight and vertex preparation run exactly
    as in a real solve."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    idict = _write_chain_inputs(d, fx, rows_by_type)
    reader = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": fx.T, "mu": fx.mu, "CellShape": [fx.L, 1, 1], "SubShape": [1, 1, 1],
           "Nmat": nmat, "IterationMax": 1, "Mix": 1.0, "EPS": 1,
           "flex_second_order": "local"}
    if gate:
        par["flex_hartree_fock"] = True
        par["longitudinal_bond_channels"] = True
    return flex_mod.FLEX(reader.get_param("ham"), {},
                         {"mode": "FLEX", "param": par, "enable_spin_orbital": False,
                          "calc_scheme": "general"})


def _production_sigma(fx, rows_by_type, gate, nmat):
    """Production ``Sigma_fluct`` of ONE map on the BARE Green function at
    the fixed ``mu``, and that Green function.

    ``gate=False`` is the standalone general path
    (``_flex_compute_veff_general`` -> ``_calc_self_energy_general``);
    ``gate=True`` is Phase B's bond-resolved path
    (``dress_and_build_w`` -> ``calc_self_energy_bond``), driven exactly as
    ``tests/test_flex_second_order_sopt.py::_one_map_sigma`` drives it."""
    from hwave.solver import flex_bond
    beta, L, norb = fx.beta, fx.L, fx.norb
    shape = (L, 1, 1)
    with tempfile.TemporaryDirectory() as d, _quiet():
        s = _chain_solver(d, fx, rows_by_type, gate, nmat)
        if gate:
            green_info = {}
            s._phase_b_reset(green_info)
            s._phase_b_preflight(green_info)
        s._calc_epsilon_k({})
        G = s._calc_dressed_green(beta, fx.mu,
                                  np.zeros((1, nmat, L, norb, norb), complex))
        if not gate:
            chi0q_raw = s._calc_chi0q(G, np.zeros_like(G), beta)[0]
            _, v_eff, _, _ = s._flex_compute_veff_general(chi0q_raw, s.ham_info.ham_inter_q)
            return s._calc_self_energy_general(G, v_eff, beta)[0], G
        nd = norb * norb
        nb = s._bond_view.n_channels
        with flex_bond.BondBlockStore(nmat, L, nb * nd, nd, ("chibar", "W")) as store:
            s._phase_b_prepare_vertices()
            flex_bond.assemble_bubble(store, G, np.zeros_like(G), beta, s._bond_view, shape, 1)
            flex_bond.dress_and_build_w(store, s._bond_S, s._bond_C, S_on=s._bond_S_on,
                                        C_on=s._bond_C_on, nb=nmat, output_full=False,
                                        nmat=nmat, nvol=L, nd=nd, spatial_shape=shape,
                                        factors=s._second_order_factors, second_order="local")
            return flex_bond.calc_self_energy_bond(store, G, beta, s._bond_view,
                                                   shape, norb, 1)[0], G


class _Maps(object):
    """Memoised ED / production maps of one fixture: the pure-type points of
    a cross stencil are the same maps the pure-type stencils already asked
    for, and the gate-on and gate-off runs of a type share nothing but the
    ED side."""

    def __init__(self, fx, rows_of, nmat=_NMAT):
        self.fx, self.rows_of, self.nmat = fx, rows_of, nmat
        self.idx = _inner()
        self.iws = 1j * (2 * self.idx + 1 - nmat) * np.pi / fx.beta
        self._ed, self._pr = {}, {}
        self._bare = None

    def _zero(self):
        return np.zeros((len(self.idx), self.fx.L, self.fx.norb, self.fx.norb), complex)

    def ed(self, *args):
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        key = args
        if key not in self._ed:
            self._ed[key] = _sigma_fluct_ed(self.fx, rows, self.iws)
        return self._ed[key]

    def prod(self, *args, **kw):
        gate = kw.pop("gate", False)
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        key = (gate,) + args
        if key not in self._pr:
            self._pr[key] = _production_sigma(self.fx, rows, gate, self.nmat)[0][self.idx]
        return self._pr[key]

    def bare_green(self, rows):
        """The bare Green function the oracle's dropped class is evaluated
        on -- the SAME propagator both sides expand around. It is
        coupling-independent (no SCF, fixed mu), so it is built once from
        whichever map asks for it first and reused."""
        if self._bare is None:
            self._bare = _production_sigma(self.fx, rows, False, self.nmat)[1]
        return self._bare

    def dropped(self, *args):
        """The oracle's off/off UNCROSSED class at the given couplings. The
        oracle's skeleton is EXACTLY quadratic in the couplings, so a
        coefficient needs no stencil: the value at unit coupling IS it."""
        rows = self.rows_of(*args)
        if not rows:
            return self._zero()
        recs = oracle_records(rows, self.fx.norb, (self.fx.L, 1, 1))
        return oracle_sigma2(self.bare_green(rows), self.fx.beta, recs, self.fx.norb,
                             "dropped", (self.fx.L, 1, 1))[0][self.idx]


# --------------------------------------------------------------------------
# Tests
# --------------------------------------------------------------------------

class TestGreenFunction(unittest.TestCase):
    """``SectorED.green`` / ``density_matrix``: the single-particle
    observables the chain gate is built on, pinned where they are known in
    closed form or by an independent dense evaluation."""

    def test_free_chain_green_equals_the_one_body_resolvent(self):
        """With no interaction the cross-sector Lehmann sum must collapse to
        ``[i w - (h1 - mu)]^{-1}`` on every mode pair -- which pins the
        Boltzmann weighting, the grand-canonical denominators (``mu``
        enters them, the sectors differ by one particle) and the global
        partition function all at once."""
        fx = _fx_chain()
        iws = 1j * (2 * np.arange(8) + 1 - 8) * np.pi / fx.beta
        G = edu.green_function(fx, iws)
        h1 = fx.build_h1()
        eye = np.eye(fx.nmode)
        for i, iw in enumerate(iws):
            np.testing.assert_allclose(
                G[i], np.linalg.inv(iw * eye - (h1 - fx.mu * eye)), atol=1e-12)

    def test_single_site_interacting_green_equals_the_dense_lehmann_sum(self):
        """One site, one orbital, ``U n_up n_dn``: the sector engine against
        the DENSE Lehmann sum over the full Fock space
        (``EDFixture.annihilators``), which knows nothing about sectors.
        The free pin above cannot see a sector-pairing error that the
        interaction makes matter; this one can."""
        fx = edu.EDFixture(L=1, norb=1, t={(0, 0): 0.0}, eps=(0.3,), T=0.5, mu=0.2)
        terms = edu.canonical_density_terms(fx, [(0, 0, 0, 0.7)])
        iws = 1j * (2 * np.arange(6) + 1 - 6) * np.pi / fx.beta
        G = edu.green_function(fx, iws, terms)

        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        h1 = fx.build_h1()
        H = edu.h_int_from_terms(fx, terms)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                if h1[p, q] != 0:
                    H = H + h1[p, q] * (CD[p] @ C[q])
        N = sum(CD[p] @ C[p] for p in range(fx.nmode))
        E, V = np.linalg.eigh(H - fx.mu * N)
        E = E - E.min()
        w = np.exp(-fx.beta * E)
        cm = [V.conj().T @ C[p] @ V for p in range(fx.nmode)]
        boltz = w[:, None] + w[None, :]
        dE = E[None, :] - E[:, None]
        ref = np.zeros_like(G)
        for p in range(fx.nmode):
            for q in range(fx.nmode):
                num = cm[p] * np.conj(cm[q]) * boltz
                for i, iw in enumerate(iws):
                    ref[i, p, q] = (num / (iw - dE)).sum()
        ref /= w.sum()
        self.assertGreater(np.abs(ref).max(), 1e-3)            # anti-vacuity
        np.testing.assert_allclose(G, ref, atol=1e-12)

    def test_free_density_matrix_equals_the_fermi_occupation(self):
        """``SectorED.density_matrix`` at zero coupling is the free
        ``<c^dag_p c_q> = [W f(ev) W^dagger]^T`` of ``h1 - mu``, the same
        expression ``hf_h1_from_terms`` Wick-contracts against."""
        fx = _fx_orbital()
        h1 = fx.build_h1()
        ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
        f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
        np.testing.assert_allclose(edu.SectorED(fx).density_matrix(),
                                   ((W * f) @ W.conj().T).T, atol=1e-12)


class TestChainHamiltonian(unittest.TestCase):
    """What the chain ED must be for the comparison to mean anything: the
    solver's own band structure, and the solver's own Hamiltonian."""

    def test_bare_green_matches_the_production_solver(self):
        """The bare ``G(k, i w)`` the solver builds from the chain's
        ``transfer.dat`` equals ``_g0_k``'s ``[(i w + mu) - eps(k)]^{-1}``,
        on BOTH fixtures.

        On these two fixtures the band is inversion- and transpose-symmetric
        (real hopping), so what this pin actually fixes is the
        ``(j + R, a) <- (j, b)`` reading of a transfer row, the ``kx``
        ordering and the site-major ``j * norb + a`` row layout -- checked
        against production rather than re-derived. The exponent SIGN and
        the orbital index order cannot be seen here and are pinned
        separately, on a complex-hopping band, by
        :meth:`test_fourier_sign_and_orbital_order_are_pinned`."""
        nmat = 32
        for fx in (_fx_chain(), _fx_orbital()):
            rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 1e-9, 0.0)]}
            _sig, G = _production_sigma(fx, rows, False, nmat)
            iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
            g0 = _g0_k(fx, iws)
            with self.subTest(L=fx.L, norb=fx.norb):
                self.assertGreater(np.abs(g0).max(), 1e-3)      # anti-vacuity
                np.testing.assert_allclose(G[0], g0, atol=1e-12)
                # and the ED engine's own free Green function is that same
                # object, so the Dyson inversion below subtracts like with like
                np.testing.assert_allclose(_green_ed_k(edu.SectorED(fx), fx, iws),
                                           g0, atol=1e-12)

    def test_fourier_sign_and_orbital_order_are_pinned(self):
        """On a band that is neither inversion- nor transpose-symmetric, the
        solver's own bare ``G(k, i w)`` selects ``_to_k``'s exponent SIGN
        and its orbital index order -- both of which the gate's two
        fixtures are blind to.

        Why this is load-bearing rather than decorative: ``_fx_chain`` and
        ``_fx_orbital`` have real hopping, so ``eps(k) = eps(-k)`` to
        round-off; with the exponent sign flipped every comparison in this
        module would still pass, and
        :meth:`test_ed_hamiltonian_is_the_production_mean_field` would then
        select the OPPOSITE orbital placement for a declared off-site row,
        inverting every downstream verdict with no pin going red. The same
        goes for the orbital index order on a real-symmetric ``eps(k)``.

        Three-sided: the ``e^{+ikR}`` convention matches production at
        1e-12; the reversed-k candidate (``sign = -1``, i.e. the array read
        at ``-kx``) and the orbital-transposed one both deviate by more
        than 1e-2. The two asymmetries are asserted to be O(1) first, so
        neither rejection can be vacuous."""
        nmat = 32
        fx = _fx_complex()
        iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
        eps = _to_k(fx, fx.build_h1())
        rev_idx = (-np.arange(fx.L)) % fx.L
        self.assertGreater(np.abs(eps - eps[rev_idx]).max(), 0.1,
                           "the fixture's band is inversion-symmetric: the Fourier "
                           "sign cannot be seen on it")
        self.assertGreater(np.abs(eps - np.swapaxes(eps, -1, -2)).max(), 0.1,
                           "the fixture's band is transpose-symmetric: the orbital "
                           "index order cannot be seen on it")

        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 1e-9, 0.0)]}
        _sig, G = _production_sigma(fx, rows, False, nmat)
        good = _g0_k(fx, iws)
        self.assertGreater(np.abs(good).max(), 1e-3)             # anti-vacuity
        np.testing.assert_allclose(G[0], good, atol=1e-12)
        # and the ED engine's own free Green function, through the same
        # transform, is that same object
        np.testing.assert_allclose(_green_ed_k(edu.SectorED(fx), fx, iws),
                                   good, atol=1e-12)
        reversed_k = _g0_k(fx, iws, sign=-1.0)
        self.assertGreater(np.abs(reversed_k - G[0]).max(), 1e-2,
                           "the reversed-k candidate also matches production: the "
                           "Fourier sign of _to_k is not pinned")
        transposed = np.swapaxes(good, -1, -2)
        self.assertGreater(np.abs(transposed - G[0]).max(), 1e-2,
                           "the orbital-transposed candidate also matches production: "
                           "the orbital index order of _to_k is not pinned")

    def test_ed_hamiltonian_is_the_production_mean_field(self):
        """The first-order functional of :func:`_ed_terms` equals UHFk's own
        mean-field kernel (``hartree_fock.accumulate_hf``, Hartree AND
        Fock) at the same density, for every fixture's declared rows -- and
        the two variants that would be wrong do NOT.

        This is what fixes the chain Hamiltonian, and it is two-sided in
        both conventions at once:

        * the ORBITAL PLACEMENT of an off-site row. Measured deviation from
          the production mean field with the opposite placement: 0 on the
          single-orbital chain (where the two are the same Hamiltonian),
          5.4e-2 of the mean field on the symmetric inter-orbital bond and
          1.4e-2 on the asymmetric one. This is the convention the
          campaign's second-order comparison rests on, and it is measured
          here rather than assumed.
        * the MIRRORED-ROW weight: at the full weight the off-site part of
          the mean field is doubled, 0.5 ... 1.0 of it.
        """
        cases = [("chain U", _fx_chain(), _rows_u_v(0.3, 0.0), False),
                 ("chain V", _fx_chain(), _rows_u_v(0.0, 0.3), False),
                 ("chain J", _fx_chain(), _rows_offsite("Hund", 0.3), False),
                 ("chain I", _fx_chain(), _rows_offsite("Ising", 0.3), False),
                 ("orbital single bond", _fx_orbital(), _rows_interorbital_one(0.3), True),
                 ("orbital asymmetric", _fx_orbital(), _rows_interorbital(0.3), True)]
        for (name, fx, rows, orbital_sensitive) in cases:
            h1 = fx.build_h1()
            ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
            f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
            rho = ((W * f) @ W.conj().T).T
            ref = _production_hf_sigma_k(fx, rows, rho)
            scale = np.abs(ref).max()
            with self.subTest(case=name):
                self.assertGreater(scale, 1e-3)                  # anti-vacuity
                # the production mean field is spin-diagonal and paramagnetic
                # on this density, so its up block is the whole story
                norb = fx.norb
                np.testing.assert_allclose(ref[:, :norb, norb:], 0.0, atol=1e-12)
                np.testing.assert_allclose(ref[:, norb:, :norb], 0.0, atol=1e-12)
                np.testing.assert_allclose(ref[:, norb:, norb:], ref[:, :norb, :norb],
                                           atol=1e-12)
                up = ref[:, :norb, :norb]
                mine = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows), rho))
                np.testing.assert_allclose(mine, up, atol=1e-11 * scale)
                swapped = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows, swap_orbitals=True), rho))
                if orbital_sensitive:
                    self.assertGreater(np.abs(swapped - up).max(), 1e-3 * scale,
                                       "{}: the orbital placement of the off-site row is not "
                                       "pinned -- both readings match the production mean "
                                       "field".format(name))
                else:
                    np.testing.assert_allclose(swapped, up, atol=1e-11 * scale)
                full = _to_k(fx, _hf_sigma(fx, _ed_terms(fx, rows, mirrored_weight=1.0), rho))
                if any(t in _MIRRORED_TYPES for t in rows):
                    self.assertGreater(np.abs(full - up).max(), 1e-3 * scale,
                                       "{}: the FULL-weight Hamiltonian also matches the "
                                       "production mean field -- the mirrored-row weight is "
                                       "not pinned".format(name))

    def test_onsite_block_matches_the_second_order_compiler(self):
        """For the on-site rows alone, UHFk's mean field is
        ``hf_first_order(compile_onsite(...))`` -- the object Task 10's gate
        is built on. Recorded here so the two ED gates of the campaign are
        known to subtract the SAME first order, not two independent ones."""
        from hwave.solver.second_order import compile_onsite, hf_first_order
        fx = _fx_orbital()
        rows = {"CoulombIntra": [(0, 0, 0, 1, 1, 0.3, 0.0), (0, 0, 0, 2, 2, 0.2, 0.0)],
                "CoulombInter": [(0, 0, 0, 1, 2, 0.15, 0.0), (0, 0, 0, 2, 1, 0.15, 0.0)]}
        h1 = fx.build_h1()
        ev, W = np.linalg.eigh(h1 - fx.mu * np.eye(fx.nmode))
        f = 1.0 / (np.exp(np.clip(fx.beta * ev, -500, 500)) + 1.0)
        rho = ((W * f) @ W.conj().T).T
        ref = _production_hf_sigma_k(fx, rows, rho)
        # k-independent, as an on-site mean field must be
        np.testing.assert_allclose(ref, np.broadcast_to(ref[0], ref.shape), atol=1e-12)
        norb = fx.norb
        tbl = {t: {((0, 0, 0), (r[3] - 1, r[4] - 1)): complex(r[5], r[6]) for r in rr}
               for t, rr in rows.items()}
        rho_gen = np.zeros((2 * norb, 2 * norb), dtype=complex)
        for p in range(2 * norb):
            for q in range(2 * norb):
                sp, ap = divmod(p, norb)
                sq, aq = divmod(q, norb)
                rho_gen[p, q] = rho[fx.mode(0, ap, sp), fx.mode(0, aq, sq)]
        sig1 = hf_first_order(compile_onsite(tbl, norb), rho_gen)
        self.assertGreater(np.abs(sig1).max(), 1e-3)             # anti-vacuity
        np.testing.assert_allclose(ref[0], sig1, atol=1e-11 * np.abs(sig1).max())


class TestFirstOrder(unittest.TestCase):
    """The gate that licenses the second-order comparison: what is left
    after the Hartree-Fock subtraction really does start at O(v^2)."""

    def test_first_order_remainder_vanishes(self):
        """``Sigma_ED - Sigma_HF[rho_ED]`` has NO first order -- with the
        mirrored-row half weight, and only with it.

        Measured (O(v) remainder over the self-energy's own linear size,
        half weight / full weight): chain U 4.0e-5 (no mirrored row),
        chain V 3.3e-4 / 5.4e-1, chain U+V 2.9e-4 / 5.0e-1, chain Hund
        1.3e-4 / 4.8e-1, chain Ising 6.7e-4 / 4.9e-1, orbital asymmetric
        1.9e-4 / 5.3e-1."""
        cases = [("chain U", _fx_chain(), lambda v: _rows_u_v(v, 0.0)),
                 ("chain V", _fx_chain(), lambda v: _rows_u_v(0.0, v)),
                 ("chain U+V", _fx_chain(), lambda v: _rows_u_v(v, v)),
                 ("chain Hund", _fx_chain(), lambda v: _rows_offsite("Hund", v)),
                 ("chain Ising", _fx_chain(), lambda v: _rows_offsite("Ising", v)),
                 ("orbital asymmetric", _fx_orbital(), _rows_interorbital)]
        nmat, x = 32, 0.05
        for (name, fx, rows_of) in cases:
            iws = 1j * (2 * np.arange(nmat) + 1 - nmat) * np.pi / fx.beta
            mirrored = any(t in _MIRRORED_TYPES for t in rows_of(1.0))
            for weight, must_vanish in ((0.5, True), (1.0, False)):
                if weight == 1.0 and not mirrored:
                    continue      # the two conventions build the same Hamiltonian
                lin = _refine(lambda h: _sigma_fluct_ed(fx, rows_of(h), iws, weight) / h, x)
                sig, _hf = _sigma_ed(fx, rows_of(x), iws, weight)
                scale = np.abs(sig).max() / x      # the self-energy's own linear size
                with self.subTest(case=name, mirrored_weight=weight):
                    self.assertGreater(scale, 1e-3)              # anti-vacuity
                    rel = np.abs(lin).max() / scale
                    if must_vanish:
                        self.assertLess(rel, _TOL,
                                        "{}: the first-order remainder does not vanish with the "
                                        "half weight ({:.2e} of the self-energy)".format(name, rel))
                    else:
                        self.assertGreater(rel, 0.1,
                                           "{}: the FULL-weight Hamiltonian also passes the "
                                           "first-order gate -- the weight is not pinned"
                                           .format(name))


class TestG4(unittest.TestCase):
    """Exact diagonalisation of the off-site second order on a chain
    (design 2026-09-15, gate G4)."""

    def _compare(self, label, ed, prod, dropped, gate=None):
        scale = np.abs(ed).max()
        self.assertGreater(scale, _FLOOR, "{}: the ED coefficient is at the floor".format(label))
        np.testing.assert_allclose(prod + dropped, ed, rtol=_TOL, atol=_TOL * scale,
                                   err_msg="{}: local + dropped vs exact diagonalisation"
                                   .format(label))
        if gate is not None:
            np.testing.assert_allclose(gate, ed, rtol=_TOL, atol=_TOL * scale,
                                       err_msg="{}: bond gate vs exact diagonalisation"
                                       .format(label))
        return scale

    @heavy
    def test_coulombinter_chain(self):
        """``L = 4`` chain, on-site ``U`` and nearest-neighbour ``V``: the
        O(V^2), O(U V) and O(U^2) coefficients of the exact remainder equal
        production's ``local`` plus the oracle's dropped class, and the
        gate-on production equals the exact remainder for ``V``.

        The dropped class is load-bearing here (1.2e-1 of the O(V^2)
        coefficient), so the first assertion really does distinguish the
        local weighting from the exact one."""
        maps = _Maps(_fx_chain(), _rows_u_v)
        with _quiet():
            ed_v = _refine(lambda h: _coeff2(lambda v: maps.ed(0.0, v), h), _X)
            pr_v = _refine(lambda h: _coeff2(lambda v: maps.prod(0.0, v), h), _X)
            dr_v = maps.dropped(0.0, 1.0)
            gate_v = _refine(lambda h: _coeff2(lambda v: maps.prod(0.0, v, gate=True), h), _X)
            with self.subTest(entry="V^2"):
                scale = self._compare("V^2", ed_v, pr_v, dr_v, gate=gate_v)
                self.assertGreater(np.abs(dr_v).max(), _DROPPED_FLOOR * scale,
                                   "the dropped class is negligible on this fixture: the "
                                   "assertion would pass against the exact oracle too")

            ed_u = _refine(lambda h: _coeff2(lambda u: maps.ed(u, 0.0), h), _X)
            pr_u = _refine(lambda h: _coeff2(lambda u: maps.prod(u, 0.0), h), _X)
            with self.subTest(entry="U^2"):
                # on-site only: nothing is off-site, so nothing is dropped
                self.assertEqual(np.abs(maps.dropped(1.0, 0.0)).max(), 0.0)
                self._compare("U^2", ed_u, pr_u, 0.0)

            ed_x = _refine(lambda h: _coeff11(lambda u, v: maps.ed(u, v), h, h), _X)
            pr_x = _refine(lambda h: _coeff11(lambda u, v: maps.prod(u, v), h, h), _X)
            dr_x = maps.dropped(1.0, 1.0) - dr_v - maps.dropped(1.0, 0.0)
            with self.subTest(entry="U V"):
                # mixed on/off pairs are never dropped (one entry is on-site,
                # so the pair is always "crossed"); measured 1.7e-18
                self.assertLess(np.abs(dr_x).max(), _FLOOR * np.abs(ed_x).max())
                self._compare("U V", ed_x, pr_x, dr_x)

    @heavy
    def test_interorbital_bond_chain(self):
        """``L = 3``, two orbitals, ORBITAL-ASYMMETRIC inter-orbital bond:
        the adjudication of the campaign's two open items.

        ASSERTED -- the standalone local path plus the dropped class
        reproduces the exact O(v^2) remainder (3.5e-4 of it). Since the ED
        Hamiltonian is pinned to UHFk's own mean field
        (:meth:`TestChainHamiltonian.test_ed_hamiltonian_is_the_production_
        mean_field`), this is an independent confirmation of the ORIENTATION
        of the off-site density vertex: the OPPOSITE orbital placement of
        the same declaration, run through the identical recipe, misses
        ``local + dropped`` by 5.5e-2 -- twenty-seven times the tolerance,
        and the assertion below would fail loudly on it.

        RECORDED -- the bond gate's deviation. Phase B's bond-resolved path
        is supposed to resum exactly the class the local kernel drops
        (2.2e-2 of the coefficient here), and on the single-orbital chain it
        does (2.4e-4). On this inter-orbital bond it misses by 4.0e-2 --
        twenty times the tolerance, and nearly twice the class it is
        recovering, so the deviation is not that class being mishandled but
        a bond-gate defect on inter-orbital off-site bonds. It is printed
        and routed to a Phase B follow-up of issue #181; it is NOT a
        second-order kernel finding, and it is not asserted here."""
        maps = _Maps(_fx_orbital(), _rows_interorbital)
        with _quiet():
            ed = _refine(lambda h: _coeff2(lambda v: maps.ed(v), h), _X)
            pr = _refine(lambda h: _coeff2(lambda v: maps.prod(v), h), _X)
            dr = maps.dropped(1.0)
            gate = _refine(lambda h: _coeff2(lambda v: maps.prod(v, gate=True), h), _X)
            scale = self._compare("inter-orbital v^2", ed, pr, dr)
            self.assertGreater(np.abs(dr).max(), _DROPPED_FLOOR * scale,
                               "the dropped class is negligible on this fixture")
            print("RECORDED bond gate on the inter-orbital off-site bond: relative "
                  "deviation from exact diagonalisation {:.3e} (dropped class {:.3e} of "
                  "the coefficient) -- Phase B follow-up, issue #181"
                  .format(np.abs(gate - ed).max() / scale, np.abs(dr).max() / scale))

    @heavy
    def test_hund_ising_chain_recorded(self):
        """Off-site Hund and Ising on the ``L = 4`` chain, RECORDED (spec
        G2 (c) / G4): the same recipe as the CoulombInter chain, printed
        rather than asserted, because the off-site Hund/Ising second order
        is the follow-up class of issue #181's Tier 3 and not part of this
        gate's asserted surface. Measured: 2.5e-4 (Hund) and 2.2e-4
        (Ising) for both the standalone and the gate-on path."""
        for itype in ("Hund", "Ising"):
            maps = _Maps(_fx_chain(), lambda v, itype=itype: _rows_offsite(itype, v))
            with _quiet():
                ed = _refine(lambda h: _coeff2(lambda v: maps.ed(v), h), _X)
                pr = _refine(lambda h: _coeff2(lambda v: maps.prod(v), h), _X)
                gate = _refine(lambda h: _coeff2(lambda v: maps.prod(v, gate=True), h), _X)
                dr = maps.dropped(1.0)
            scale = np.abs(ed).max()
            self.assertGreater(scale, _FLOOR)                    # anti-vacuity
            print("RECORDED off-site {} on the L = 4 chain: local + dropped {:.3e}, "
                  "gate on {:.3e} (dropped class {:.3e} of the coefficient)"
                  .format(itype, np.abs(pr + dr - ed).max() / scale,
                          np.abs(gate - ed).max() / scale, np.abs(dr).max() / scale))


if __name__ == "__main__":
    unittest.main()
