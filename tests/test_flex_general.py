#!/usr/bin/env python3

"""Tests for the general-mode (paramagnetic full-vertex) FLEX guards.

These tests cover Task 1: ``calc_scheme="general"`` must be accepted for
paramagnetic (spin-free) FLEX and rejected (fail-fast ValueError) for
``enable_spin_orbital`` (at construction) and for non-spin-free spin modes
(spin-diag / spinful, checked inside ``solve()``).

Tests must be run from the repository root (they use relative paths like
``tests/rpa/input``).
"""

import os
import unittest
import numpy as np


def _make_solver(mode_cls, Lx=8, Ly=8, Nmat=64, T=2.0, mu=0.0,
                 calc_scheme="reduced", extra_mode=None,
                 extra_interactions=None):
    """Build a FLEX/RPA solver from the 1-orbital ``tests/rpa/input`` data.

    Replicates the body of ``tests/test_flex.py``'s ``TestFLEX._make_solver``
    helper so these tests are self-contained.
    """
    info_log = {}
    info_mode = {
        'mode': mode_cls,
        'param': {
            'T': T,
            'mu': mu,
            'CellShape': [Lx, Ly, 1],
            'SubShape': [1, 1, 1],
            'Nmat': Nmat,
        },
        'calc_scheme': calc_scheme,
    }
    if extra_mode:
        info_mode.update(extra_mode)

    info_file = {
        'input': {
            'path_to_input': 'tests/rpa/input',
            'interaction': {
                'path_to_input': 'tests/rpa/input',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        },
        'output': {
            'path_to_output': 'tests/flex/output',
        },
    }
    if extra_interactions:
        info_file['input']['interaction'].update(extra_interactions)

    os.makedirs('tests/flex/output', exist_ok=True)

    import hwave.qlmsio.read_input_k as read_input_k
    read_io = read_input_k.QLMSkInput(info_file['input'])
    ham_info = read_io.get_param("ham")
    green_info = read_io.get_param("green")

    if mode_cls == "FLEX":
        import hwave.solver.flex as solver_flex
        solver = solver_flex.FLEX(ham_info, info_log, info_mode)
    else:
        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

    return solver, green_info, info_file


def _make_general_flex():
    """Build a spin-free 2-orbital general FLEX solver.

    Uses the ``tests/rpa/input_2orb`` fixture (norb=2, CoulombIntra +
    CoulombInter) so the MYO S/C matrices are nontrivial 4x4 blocks.  The
    constructor sets ``self.norb``, ``self.lattice`` and ``self.ham_info``
    (with ``self.ham_info.param_ham``), which is all the inflation method
    needs.
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    info_input = {
        'path_to_input': 'tests/rpa/input_2orb',
        'interaction': {
            'path_to_input': 'tests/rpa/input_2orb',
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
        },
    }
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    info_mode = {
        'mode': 'FLEX',
        'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                  'SubShape': [1, 1, 1], 'Nmat': 8},
        'calc_scheme': 'general',
    }
    solver = solver_flex.FLEX(ham, {}, info_mode)
    solver.spin_mode = "spin-free"
    return solver


def _fake_general_chi0q(flex):
    """Build a random complex rank-4 chi0q matching the general-path layout.

    Shape ``(nmat, nvol, norb, norb, norb, norb)`` -- this is what
    ``RPA._calc_chi0q`` returns for the spin-free general scheme after the
    block dimension has been stripped in ``FLEX.solve``.
    """
    rng = np.random.default_rng(0)
    no = flex.norb
    shape = (flex.nmat, flex.lattice.nvol, no, no, no, no)
    return (rng.standard_normal(shape)
            + 1j * rng.standard_normal(shape)).astype(np.complex128)


class TestInflateGeneral(unittest.TestCase):
    """Task 4: rank-4 chi0q pass-through + full MYO S/C inflation."""

    def test_inflate_returns_rank4_chi0q_and_sc_matrices(self):
        flex = _make_general_flex()
        self.assertTrue(flex._flex_general)
        no = flex.norb
        self.assertEqual(no, 2)

        chi0q_raw = _fake_general_chi0q(flex)
        chi0q, Us, Uc = flex._inflate_chi0q_and_ham_general(chi0q_raw, None)

        # chi0q passes through unchanged: rank-4 in orbital space.
        nvol = flex.lattice.nvol
        self.assertEqual(chi0q.ndim, 6)
        self.assertEqual(chi0q.shape,
                         (flex.nmat, nvol, no, no, no, no))
        np.testing.assert_array_equal(chi0q, chi0q_raw)

        # S/C matrices reshaped to (nvol, norb^2, norb^2).
        self.assertEqual(Us.shape, (nvol, no * no, no * no))
        self.assertEqual(Uc.shape, (nvol, no * no, no * no))
        self.assertEqual(Us.shape[-2:], (4, 4))
        self.assertEqual(Uc.shape[-2:], (4, 4))

        # CoulombIntra + CoulombInter present -> nonzero S/C.
        self.assertTrue(np.any(Us != 0.0))
        self.assertTrue(np.any(Uc != 0.0))


class TestFLEXGeneralGuards(unittest.TestCase):
    """Guards for calc_scheme='general' FLEX (v1, paramagnetic only)."""

    def test_general_accepted_for_spin_free(self):
        """calc_scheme='general' must construct without raising (spin-free)."""
        solver, green_info, _ = _make_solver("FLEX", calc_scheme="general")
        self.assertIsNotNone(solver)
        self.assertTrue(solver._flex_general)

    def test_general_rejected_for_enable_spin_orbital(self):
        """calc_scheme='general' + enable_spin_orbital must raise at construct."""
        with self.assertRaises(ValueError):
            _make_solver("FLEX", calc_scheme="general",
                         extra_mode={'enable_spin_orbital': True})

    def test_spin_mode_guard_in_solve(self):
        """The general + non-spin-free guard in solve() must raise.

        Constructing a genuine spin-diag input would require a magnetic-field /
        spin-dependent transfer setup that the 1-orbital fixture does not
        provide.  The guard runs in solve() right after _calc_epsilon_k (which
        is what determines spin_mode from H0), so we stub _calc_epsilon_k to
        report a spin-diag Hamiltonian and confirm solve() fails fast with a
        ValueError before doing any real FLEX work.
        """
        solver, green_info, info_file = _make_solver(
            "FLEX", calc_scheme="general")
        self.assertTrue(solver._flex_general)

        def fake_epsilon_k(gi):
            solver.spin_mode = "spin-diag"

        solver._calc_epsilon_k = fake_epsilon_k
        with self.assertRaises(ValueError):
            solver.solve(green_info, info_file['output']['path_to_output'])

    def test_reduced_still_works(self):
        """Regression: calc_scheme='reduced' FLEX still constructs."""
        solver, green_info, _ = _make_solver("FLEX", calc_scheme="reduced")
        self.assertIsNotNone(solver)
        self.assertFalse(solver._flex_general)

    def test_general_rejected_for_ring_ladder(self):
        """calc_type='ring+ladder' forces scheme='general' in RPA, but FLEX
        general is ring-only — it must still be rejected (not silently accepted
        as a plain general run)."""
        with self.assertRaises(ValueError):
            _make_solver("FLEX", calc_scheme="general",
                         extra_mode={'calc_type': 'ring+ladder'})

    def test_general_rejected_for_spinful_at_solve(self):
        """The solve()-time guard must also reject spin_mode='spinful'."""
        solver, green_info, info_file = _make_solver(
            "FLEX", calc_scheme="general")

        def fake_epsilon_k(gi):
            solver.spin_mode = "spinful"

        solver._calc_epsilon_k = fake_epsilon_k
        with self.assertRaises(ValueError):
            solver.solve(green_info, info_file['output']['path_to_output'])


class TestFLEXGeneralWarningGating(unittest.TestCase):
    """The density-density reduction warning must fire for reduced/squashed but
    be suppressed for general (where the off-diagonal vertices are kept)."""

    def _construct(self, scheme):
        """Construct a 2-orbital FLEX while pretending an Exchange interaction is
        present. ``self.ham_info`` is a fresh ``Interaction`` built inside
        ``RPA.__init__`` (not the passed object), so patch the class method
        ``Interaction.has_interaction_exchange`` rather than an instance — this
        also avoids wiring an extra interaction file."""
        from unittest.mock import patch
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex
        import hwave.solver.rpa as solver_rpa
        info_input = {
            'path_to_input': 'tests/rpa/input_2orb',
            'interaction': {
                'path_to_input': 'tests/rpa/input_2orb',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            },
        }
        ham = read_input_k.QLMSkInput(info_input).get_param("ham")
        info_mode = {
            'mode': 'FLEX',
            'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                      'SubShape': [1, 1, 1], 'Nmat': 32},
            'calc_scheme': scheme,
        }
        with patch.object(solver_rpa.Interaction, 'has_interaction_exchange',
                          return_value=True):
            return solver_flex.FLEX(ham, {}, info_mode)

    def test_warning_retained_for_squashed(self):
        # 'squashed' (not 'reduced': RPA errors on reduced+exchange, rpa.py:643).
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            self._construct('squashed')
        self.assertTrue(any('density-density' in m for m in cm.output))

    def test_warning_suppressed_for_general(self):
        with self.assertNoLogs('hwave.solver.flex', level='WARNING'):
            self._construct('general')


def _kanamori_inter_k(norb=2, U=4.0, Up=2.0, J=0.5, Jp=0.5,
                      Nx=2, Ny=2, Nz=1, with_pairhop=False):
    """Build a constant Kanamori ``inter_k`` dict for the S/C matrix builders.

    Each entry is an ndarray of shape ``(norb, norb, Nx, Ny, Nz)`` matching the
    layout consumed by ``sc._build_sc_matrices_all_q`` / ``_get(itype)``.
    """
    shape = (norb, norb, Nx, Ny, Nz)
    CoulombIntra = np.zeros(shape, dtype=complex)
    CoulombInter = np.zeros(shape, dtype=complex)
    Hund = np.zeros(shape, dtype=complex)
    Exchange = np.zeros(shape, dtype=complex)
    for a in range(norb):
        CoulombIntra[a, a, :] = U
        for b in range(norb):
            if a != b:
                CoulombInter[a, b, :] = Up
                Hund[a, b, :] = J
                Exchange[a, b, :] = Jp
    inter_k = {
        "CoulombIntra": CoulombIntra,
        "CoulombInter": CoulombInter,
        "Hund": Hund,
        "Exchange": Exchange,
    }
    if with_pairhop:
        PairHop = np.zeros(shape, dtype=complex)
        for a in range(norb):
            for b in range(norb):
                if a != b:
                    PairHop[a, b, :] = Jp
        inter_k["PairHop"] = PairHop
    return inter_k


class TestMYOSCMatrices(unittest.TestCase):
    """MYO-convention S/C interaction matrix builder (cond-mat/0407094 Eq.(6))."""

    def test_myo_elements(self):
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        U, Up, J, Jp = 4.0, 2.0, 0.5, 0.5
        ik = _kanamori_inter_k(norb=2, U=U, Up=Up, J=J, Jp=Jp,
                               Nx=2, Ny=2, Nz=1)
        S, C = build_sc_matrices_myo(ik, 2, 2, 2, 1)

        def el(M, l1, l2, l3, l4):
            return M[0, 0, 0, l1 * 2 + l2, l3 * 2 + l4]

        # Case 2 (ab,ab): the MYO-specific charge element
        self.assertAlmostEqual(el(S, 0, 1, 0, 1), Up)
        self.assertAlmostEqual(el(C, 0, 1, 0, 1), -Up + 2 * J)
        # Case 3 (aa,bb)
        self.assertAlmostEqual(el(S, 0, 0, 1, 1), J)
        self.assertAlmostEqual(el(C, 0, 0, 1, 1), 2 * Up - J)
        # Case 4 (ab,ba)
        self.assertAlmostEqual(el(S, 0, 1, 1, 0), Jp)
        self.assertAlmostEqual(el(C, 0, 1, 1, 0), Jp)
        # Case 1 (aaaa)
        self.assertAlmostEqual(el(S, 0, 0, 0, 0), U)
        self.assertAlmostEqual(el(C, 0, 0, 0, 0), U)

    def test_diverges_from_kuroki_only_at_charge_abab(self):
        from hwave.sc import _build_sc_matrices_all_q
        from hwave.solver._sc_matrices_myo import build_sc_matrices_myo
        U, Up, J, Jp = 4.0, 2.0, 0.5, 0.5
        ik = _kanamori_inter_k(norb=2, U=U, Up=Up, J=J, Jp=Jp,
                               Nx=2, Ny=2, Nz=1)
        Sm, Cm = build_sc_matrices_myo(ik, 2, 2, 2, 1)
        Sk, Ck = _build_sc_matrices_all_q(ik, 2, 2, 2, 1)

        # Spin matrices identical.
        np.testing.assert_allclose(Sm, Sk)

        # Charge differs ONLY at the (ab,ab) entries, by exactly +J there.
        diff = Cm - Ck
        # (ab,ab): idx12 == idx34 with l1!=l2; for norb=2 these are
        # (0,1)&(0,1) -> flat (1,1) and (1,0)&(1,0) -> flat (2,2).
        nonzero = [(1, 1), (2, 2)]
        mask = np.zeros((4, 4), dtype=bool)
        for (i, j) in nonzero:
            mask[i, j] = True
        nz = diff[..., mask]
        np.testing.assert_allclose(nz, J)
        zero = diff[..., ~mask]
        np.testing.assert_allclose(zero, 0.0, atol=1.0e-12)


class TestBruteForceRef(unittest.TestCase):
    """Structural sanity checks for the physical-index brute-force reference
    (``tests/flex_bruteforce_ref.py``), which is the independent ground truth
    for the optimized general-mode FLEX path."""

    def test_chi0_shape(self):
        from tests.flex_bruteforce_ref import chi0_bruteforce
        rng = np.random.default_rng(0)
        norb, Nk, nmat = 2, 4, 4
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))
        chi0 = chi0_bruteforce(G, T=1.5, Nk=Nk)
        self.assertEqual(chi0.shape, (2, 2, 2, 2, 4, 4))

    def test_sigma_shape(self):
        from tests.flex_bruteforce_ref import sigma_bruteforce
        rng = np.random.default_rng(1)
        norb, Nk, nmat = 2, 4, 4
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))
        V = rng.standard_normal((norb, norb, norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, norb, norb, Nk, nmat))
        Sig = sigma_bruteforce(G, V, T=1.5, Nk=Nk)
        self.assertEqual(Sig.shape, (2, 2, 4, 4))

    def test_chi0_linearity(self):
        """chi0 is bilinear in G: chi0(2 G) == 4 chi0(G)."""
        from tests.flex_bruteforce_ref import chi0_bruteforce
        rng = np.random.default_rng(2)
        norb, Nk, nmat = 2, 3, 4
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))
        c1 = chi0_bruteforce(G, T=2.0, Nk=Nk)
        c2 = chi0_bruteforce(2 * G, T=2.0, Nk=Nk)
        np.testing.assert_allclose(c2, 4 * c1, atol=1.0e-12)

    def test_sigma_linear_in_V(self):
        """Sigma is linear in V: Sigma(G, 2 V) == 2 Sigma(G, V)."""
        from tests.flex_bruteforce_ref import sigma_bruteforce
        rng = np.random.default_rng(3)
        norb, Nk, nmat = 2, 3, 4
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))
        V = rng.standard_normal((norb, norb, norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, norb, norb, Nk, nmat))
        s1 = sigma_bruteforce(G, V, T=2.0, Nk=Nk)
        s2 = sigma_bruteforce(G, 2 * V, T=2.0, Nk=Nk)
        np.testing.assert_allclose(s2, 2 * s1, atol=1.0e-12)

    def test_sigma_matches_naive_einsum(self):
        """Cross-check sigma_bruteforce against an independent einsum/roll
        implementation of the SAME MYO Eq.(3) -- guards against a loop typo."""
        from tests.flex_bruteforce_ref import sigma_bruteforce
        rng = np.random.default_rng(4)
        norb, Nk, nmat = 2, 4, 4
        T = 1.3
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))
        V = rng.standard_normal((norb, norb, norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, norb, norb, Nk, nmat))

        Sig_ref = sigma_bruteforce(G, V, T=T, Nk=Nk)

        # Independent einsum/roll path. For each (q, iv) shift, build
        # Gs[mu,nu,k,iw] = G[mu,nu,(k-q)%Nk,(iw-iv)%nmat] via np.roll
        # (roll by +q on the k axis and +iv on the iw axis), then contract
        # over (mu, nu) with V[:,m,:,n,q,iv] and accumulate.
        Sig_ein = np.zeros((norb, norb, Nk, nmat), dtype=complex)
        for q in range(Nk):
            for iv in range(nmat):
                Gs = np.roll(np.roll(G, q, axis=2), iv, axis=3)
                # V[mu,m,nu,n] for this (q,iv) ; contract mu,nu with Gs[mu,nu,k,iw]
                Sig_ein += np.einsum('ambn,abki->mnki',
                                     V[:, :, :, :, q, iv], Gs)
        Sig_ein *= (T / Nk)

        np.testing.assert_allclose(Sig_ref, Sig_ein, atol=1.0e-12)

    def test_chi0_matches_naive_einsum(self):
        """Cross-check chi0_bruteforce against an independent einsum/roll
        implementation of the SAME MYO Eq.(5) -- guards against a loop typo
        (closes the coverage gap flagged in the Codex review; sigma already
        had such a cross-check, chi0 did not)."""
        from tests.flex_bruteforce_ref import chi0_bruteforce
        rng = np.random.default_rng(7)
        norb, Nk, nmat = 2, 4, 4
        T = 1.7
        G = rng.standard_normal((norb, norb, Nk, nmat)) \
            + 1j * rng.standard_normal((norb, norb, Nk, nmat))

        chi0_ref = chi0_bruteforce(G, T=T, Nk=Nk)

        # Independent einsum/roll path. For each (q, iv) shift, build
        # Gs[mu,m,k,iw] = G[mu,m,(k+q)%Nk,(iw+iv)%nmat] via np.roll
        # (roll by -q on the k axis and -iv on the iw axis), then contract
        # over (k, iw) with G[n,nu,k,iw]:  chi0[m,n,mu,nu] = -(T/Nk) sum Gs*G.
        chi0_ein = np.zeros((norb, norb, norb, norb, Nk, nmat), dtype=complex)
        for q in range(Nk):
            for iv in range(nmat):
                Gs = np.roll(np.roll(G, -q, axis=2), -iv, axis=3)
                chi0_ein[:, :, :, :, q, iv] = np.einsum('amki,nbki->mnab', Gs, G)
        chi0_ein *= -(T / Nk)

        np.testing.assert_allclose(chi0_ref, chi0_ein, atol=1.0e-12)


if __name__ == '__main__':
    unittest.main()
