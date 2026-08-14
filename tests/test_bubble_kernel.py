"""Gates for the unified particle-hole bubble kernel (Step 2 of the H-wave
bubble-unification series).

Binding spec: ``docs/superpowers/specs/2026-08-14-unified-bubble-kernel-design.md``.
This module grows through the series' tasks; it starts with the Green
builder (``hwave.solver.green.build_green``, extracted from
``RPA._calc_green`` in Task 1) and the dense bubble transport
(``hwave.solver.bubble.dense_bubble``, extracted from ``RPA._calc_chi0q``
in Task 2).
"""

import os
import tempfile
import unittest

import numpy as np

from tests.approx_util import assert_approx_array
from hwave.solver import green as green_mod
from hwave.solver import bubble as bubble_mod
from hwave.solver.bond_channels import resolve_interactions


def _tiny_eig(nvol=4, p=2, seed=7):
    rng = np.random.default_rng(seed)
    ev = rng.uniform(-2.0, 2.0, size=(nvol, p))
    # random unitary per k (QR of a random complex matrix)
    m = rng.normal(size=(nvol, p, p)) + 1j * rng.normal(size=(nvol, p, p))
    q, _ = np.linalg.qr(m)
    return ev, q


def _make_rpa_solver(coeff_tail=0.0, Lx=4, Ly=4, Nmat=8, calc_scheme='general',
                      norb=1, complex_hop=False, hz=0.0):
    """Build a minimal RPA solver with a real eigen-decomposition.

    The default (``norb=1``, ``complex_hop=False``, ``hz=0.0``) reads the
    ``tests/rpa/input`` fixture files (same construction pattern as
    ``tests/test_chi0q_coeff_tail_provenance.py``'s ``_make_rpa`` --
    ``tests/test_rpa_chi0q_guards.py`` only builds a bare
    ``object.__new__`` stub for its shape-guard tests, so this module
    copies the fuller construction pattern used elsewhere in the suite for
    tests that need a real eigen-decomposition).

    ``norb=2``/``complex_hop=True``/``hz!=0`` extend the same pattern with
    hand-written temporary Wannier90-like input files (the
    ``tests/rpa/input`` fixture only has REAL hopping): ``complex_hop``
    gives each orbital a complex NN hopping ``t = 0.7 * exp(0.3j)`` (a
    different magnitude per orbital, so the two orbitals are not
    degenerate); ``hz`` writes an ``Extern`` onsite field and sets
    ``coeff_extern``, which -- for a non-spin-orbital Hamiltonian --
    RPA._calc_epsilon_k splits into two spin-diag blocks (``nblock=2``,
    ``H0 +/- H1``); rpa.py:2778-2797.
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.rpa as sol_rpa

    param = {'T': 2.0, 'mu': 0.0,
             'CellShape': [Lx, Ly, 1], 'SubShape': [1, 1, 1], 'Nmat': Nmat,
             'coeff_tail': coeff_tail}

    if norb == 1 and not complex_hop and hz == 0.0:
        info_input = {
            'path_to_input': 'tests/rpa/input',
            'interaction': {
                'path_to_input': 'tests/rpa/input',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        }
        info_mode = {'mode': 'RPA', 'param': param, 'calc_scheme': calc_scheme}
        ham = read_input_k.QLMSkInput(info_input).get_param("ham")
        solver = sol_rpa.RPA(ham, {}, info_mode)
        solver._init_wavevec()
        return solver

    # Extended fixture: hand-written temp input files, read immediately
    # (QLMSkInput reads eagerly at construction), so the TemporaryDirectory
    # need not outlive this function.
    t = 0.7 * np.exp(0.3j)
    with tempfile.TemporaryDirectory() as d:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n")
            f.write("{}\n".format(norb))
            for _ in range(norb):
                f.write("0.0 0.0 0.0\n")
        # Onsite inter-orbital mixing (Hermitian, R=(0,0,0)): without it
        # H0(k) stays orbital-diagonal for a diagonal-hopping fixture, so
        # its eigenvectors are trivial and G[..., a, b] is zero off the
        # orbital diagonal at every k -- which cannot distinguish
        # contract_general's correct (..., a, c, b, d) axis order from
        # certain wrong permutations (round-1 review: a diagonal G is
        # symmetric under exactly the axis swap a bug would introduce).
        # The mixing makes G genuinely non-diagonal in orbital space.
        mix = 0.4 * np.exp(-0.6j) if norb >= 2 else 0.0
        n_mix = 2 if norb >= 2 else 0
        nr = 2 * norb + n_mix
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("bubble kernel gate fixture\n")
            f.write("{}\n{}\n".format(norb, nr))
            f.write(("1 " * nr).strip() + "\n")
            for orb in range(1, norb + 1):
                to = t * orb  # distinct magnitude per orbital
                f.write(" 1 0 0 {o} {o} {re:.12f} {im:.12f}\n".format(
                    o=orb, re=to.real, im=to.imag))
                f.write("-1 0 0 {o} {o} {re:.12f} {im:.12f}\n".format(
                    o=orb, re=np.conj(to).real, im=np.conj(to).imag))
            if norb >= 2:
                f.write(" 0 0 0 1 2 {re:.12f} {im:.12f}\n".format(
                    re=mix.real, im=mix.imag))
                f.write(" 0 0 0 2 1 {re:.12f} {im:.12f}\n".format(
                    re=np.conj(mix).real, im=np.conj(mix).imag))
        with open(os.path.join(d, "coulombintra.dat"), "w") as f:
            f.write("CoulombIntra gate fixture\n")
            f.write("{}\n{}\n".format(norb, norb))
            f.write(("1 " * norb).strip() + "\n")
            for orb in range(1, norb + 1):
                f.write(" 0 0 0 {o} {o} 1.0 0.0\n".format(o=orb))

        interaction = {
            'path_to_input': d,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        }
        if hz != 0.0:
            with open(os.path.join(d, "extern.dat"), "w") as f:
                f.write("bubble kernel gate fixture extern\n")
                f.write("{}\n{}\n".format(norb, norb))
                f.write(("1 " * norb).strip() + "\n")
                for orb in range(1, norb + 1):
                    f.write(" 0 0 0 {o} {o} 1.0 0.0\n".format(o=orb))
            interaction['Extern'] = 'extern.dat'
            param = dict(param)
            param['coeff_extern'] = hz

        info_input = {'path_to_input': d, 'interaction': interaction}
        info_mode = {'mode': 'RPA', 'param': param, 'calc_scheme': calc_scheme}
        ham = read_input_k.QLMSkInput(info_input).get_param("ham")
        solver = sol_rpa.RPA(ham, {}, info_mode)
        solver._init_wavevec()
        return solver


class TestBuildGreen(unittest.TestCase):

    def test_formula_pin_against_direct_sum(self):
        ev, V = _tiny_eig()
        beta, nmat, mu, ct = 2.0, 8, 0.1, 1.0
        full, defl, tail = green_mod.build_green(ev, V, mu, beta, nmat, ct)
        self.assertEqual(full.shape, (1, nmat, ev.shape[0], 2, 2))
        n, k = 3, 1
        iw = 1j * (2 * n + 1 - nmat) * np.pi / beta
        g_direct = (V[k] * (1.0 / (iw - (ev[k] - mu)))) @ V[k].conj().T
        assert_approx_array(full[0, n, k], g_direct, rel=0, abs=1e-13)
        assert_approx_array(full[0, n, k] - defl[0, n, k],
                             ct * (V[k] @ V[k].conj().T) / iw, rel=0, abs=1e-13)
        assert_approx_array(tail[0, n, k],
                             V[k] @ V[k].conj().T * ct * beta / 2, rel=0, abs=1e-13)

    def test_coeff_tail_zero_aliases(self):
        ev, V = _tiny_eig()
        full, defl, tail = green_mod.build_green(ev, V, 0.0, 2.0, 8, 0.0)
        self.assertIs(full, defl)
        self.assertIsNone(tail)

    def test_matches_rpa_solver_calc_green(self):
        """Oracle: RPA._calc_green on the same eigen-decomposition must
        equal build_green's (deflated, tail) pair exactly (same math,
        transcribed -- not just close)."""
        solver = _make_rpa_solver(coeff_tail=0.5)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green_old, tail_old = solver._calc_green(beta, mu)

        full_new, defl_new, tail_new = green_mod.build_green(
            solver.H0_eigenvalue, solver.H0_eigenvector, mu, beta,
            solver.nmat, solver.coeff_tail)

        assert_approx_array(defl_new, green_old, rel=0, abs=1e-13)
        assert_approx_array(tail_new, tail_old, rel=0, abs=1e-13)

    def test_rejects_unflattened_shapes(self):
        """sc.py's genuinely UNFLATTENED (Nx, Ny, Nz, p) /
        (Nx, Ny, Nz, p, p) layout (ndim 4 / 5) must be rejected -- the
        caller is required to flatten to (nvol, p) / (nvol, p, p) first.

        (A reshape that stays within the accepted ndim pairs -- e.g. to
        (nblock, nvol, p) / (nblock, nvol, p, p) -- is not a case this
        test can use: for nvol=4, p=2 such a reshape is numerically
        self-consistent as a genuine 2-block decomposition of the same
        data (verified: each reshaped (block, k) slice reproduces one
        original (k) slice exactly), so it is legitimately ACCEPTED, not
        a shape to reject. The real invalid case is the unflattened
        spatial layout build_green documents as out of contract.)"""
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev.reshape(2, 2, 1, 2), V.reshape(2, 2, 1, 2, 2),
                                   0.0, 2.0, 8, 0.0)

    def test_rejects_mismatched_p(self):
        ev, V = _tiny_eig(nvol=4, p=2)
        bad_V = np.zeros((4, 3, 3), dtype=complex)
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, bad_V, 0.0, 2.0, 8, 0.0)

    def test_rejects_nonpositive_beta(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 0.0, 8, 0.0)
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, -1.0, 8, 0.0)

    def test_rejects_odd_nmat(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 7, 0.0)

    def test_rejects_nonfinite_coeff_tail(self):
        ev, V = _tiny_eig()
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 8, float('nan'))
        with self.assertRaises(ValueError):
            green_mod.build_green(ev, V, 0.0, 2.0, 8, float('inf'))


class TestBubbleOldVsNewDense(unittest.TestCase):
    """Old-vs-new gate: ``RPA._calc_chi0q`` (OLD) vs
    ``bubble.dense_bubble`` (NEW) on identical inputs, across
    {reduced, general} x {nblock 1, 2} x {tail on (coeff_tail=0.5), off}.
    No dispatch change yet (rpa.py is untouched until Task 6) -- both
    callables exist today, so this is a pure old-vs-new numerics gate
    (Global Constraints: ``assert_approx_array(new, old, rel=1e-6,
    abs=1e-8)``, no bytewise gates, issue #85).

    The fixture is norb=2 with complex NN hopping (``t = 0.7 * exp(0.3j)``,
    a different magnitude per orbital) so the general (non-diagonal
    orbital-pair) branch and the tail's equal-time endpoint correction are
    both exercised non-trivially; ``hz != 0`` (an ``Extern`` onsite field)
    selects the nblock=2 spin-diag path.
    """

    def _oracle(self, coeff_tail, hz):
        """Build one nblock-{1,2} x tail-{on,off} solver instance and
        return (solver, green, tail, beta, spatial_shape).
        ``solver.enable_reduced`` is a plain instance attribute read
        afresh on every ``_calc_chi0q`` call (not fixed by a closure at
        construction time), so the caller toggles it per scheme instead of
        rebuilding the solver."""
        solver = _make_rpa_solver(coeff_tail=coeff_tail, hz=hz, norb=2,
                                   complex_hop=True, Lx=4, Ly=4, Nmat=8)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green, tail = solver._calc_green(beta, mu)
        spatial_shape = tuple(solver.lattice.shape)
        return solver, green, tail, beta, spatial_shape

    def test_old_vs_new_all_cells(self):
        for hz, want_nblock in ((0.0, 1), (0.3, 2)):
            for coeff_tail, tail_label in ((0.0, "off"), (0.5, "on")):
                solver, green, tail, beta, spatial_shape = self._oracle(
                    coeff_tail, hz)
                self.assertEqual(green.shape[0], want_nblock)
                for scheme in ("reduced", "general"):
                    with self.subTest(nblock=want_nblock, tail=tail_label,
                                       scheme=scheme):
                        solver.enable_reduced = (scheme == "reduced")
                        old = solver._calc_chi0q(green, tail, beta)
                        new = bubble_mod.dense_bubble(
                            green, tail, beta, spatial_shape=spatial_shape,
                            scheme=scheme)
                        self.assertEqual(new.shape, old.shape)
                        assert_approx_array(new, old, rel=1e-6, abs=1e-8)

    def test_none_tail_matches_zero_tail_array(self):
        """``green0_tail=None`` (full Green, no tail machinery) and an
        explicit all-zero tail array of the same shape (the legacy
        oracle's only vocabulary -- ``RPA._calc_chi0q`` always requires a
        real array) are two different code paths through
        ``_prepare_dense`` (skip the subtraction entirely vs. subtract an
        all-zero array) that must produce identical "tail off" numerics
        (spec: Green/tail contracts)."""
        for hz, want_nblock in ((0.0, 1), (0.3, 2)):
            solver, green, zero_tail, beta, spatial_shape = self._oracle(
                0.0, hz)
            self.assertEqual(green.shape[0], want_nblock)
            self.assertTrue(bool(np.all(zero_tail == 0)))
            for scheme in ("reduced", "general"):
                with self.subTest(nblock=want_nblock, scheme=scheme):
                    with_none = bubble_mod.dense_bubble(
                        green, None, beta, spatial_shape=spatial_shape,
                        scheme=scheme)
                    with_zero = bubble_mod.dense_bubble(
                        green, zero_tail, beta, spatial_shape=spatial_shape,
                        scheme=scheme)
                    assert_approx_array(with_none, with_zero, rel=0, abs=0)

    def test_rejects_odd_nmat(self):
        solver, green, tail, beta, spatial_shape = self._oracle(0.0, 0.0)
        bad_green = green[:, :-1]
        with self.assertRaises(ValueError):
            bubble_mod.dense_bubble(bad_green, tail[:, :-1], beta,
                                     spatial_shape=spatial_shape,
                                     scheme="general")

    def test_rejects_spatial_shape_mismatch(self):
        solver, green, tail, beta, spatial_shape = self._oracle(0.0, 0.0)
        bad_shape = (spatial_shape[0] + 1,) + spatial_shape[1:]
        with self.assertRaises(ValueError):
            bubble_mod.dense_bubble(green, tail, beta,
                                     spatial_shape=bad_shape,
                                     scheme="general")

    def test_rejects_real_dtype(self):
        solver, green, tail, beta, spatial_shape = self._oracle(0.0, 0.0)
        with self.assertRaises(ValueError):
            bubble_mod.dense_bubble(green.real, tail, beta,
                                     spatial_shape=spatial_shape,
                                     scheme="general")

    def test_rejects_bad_scheme(self):
        solver, green, tail, beta, spatial_shape = self._oracle(0.0, 0.0)
        with self.assertRaises(ValueError):
            bubble_mod.dense_bubble(green, tail, beta,
                                     spatial_shape=spatial_shape,
                                     scheme="bogus")


class TestBubbleOldVsNewIr(unittest.TestCase):
    """Old-vs-new gate: ``FLEX._calc_chi0q_ir`` / ``FLEX._calc_chi0q_general_ir``
    (OLD) vs ``bubble.ir_bubble`` (NEW) on identical inputs, reduced +
    general -- plus the IR axis-compatibility ``ValueError`` cases. No
    dispatch change yet (flex.py is untouched until Task 7).

    Fixtures reuse ``tests/test_flex_ir.py``'s / ``tests/test_flex_ir_general.py``'s
    smallest-``Nmat`` solver-construction helpers (``Nmat=64``) rather than
    re-deriving them here.

    Skips cleanly when ``sparse_ir`` is not importable -- the same guard
    ``tests/test_flex_ir.py`` uses, so ``tests.test_flex_ir`` skipping is
    the canary that this class must also skip.
    """

    @classmethod
    def setUpClass(cls):
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            raise unittest.SkipTest("sparse-ir not installed")

    def test_reduced_old_vs_new(self):
        from tests.test_flex_ir import _make_solver
        T = 0.5
        beta = 1.0 / T
        solver, gi = _make_solver(64, {'matsubara_basis': 'ir'}, T=T)
        solver._calc_epsilon_k(gi)
        solver._ir_setup(beta)
        axF, axB = solver._ir_axF, solver._ir_axB
        self.assertTrue(solver.enable_reduced)

        green, _ = solver._calc_green_ir(beta, 0.0)
        old = solver._calc_chi0q_ir(green, beta)
        spatial_shape = tuple(solver.lattice.shape)

        new = bubble_mod.ir_bubble(green, axF, axB,
                                   spatial_shape=spatial_shape,
                                   scheme="reduced")
        self.assertEqual(new.shape, old.shape)
        assert_approx_array(new, old, rel=1e-6, abs=1e-8)

    def test_general_old_vs_new(self):
        from tests.test_flex_ir_general import _make_general_solver
        T = 2.0
        beta = 1.0 / T
        solver, gi = _make_general_solver(64, "ir", T=T)
        solver._calc_epsilon_k(gi)
        solver._ir_setup(beta)
        axF, axB = solver._ir_axF, solver._ir_axB

        green, _ = solver._calc_green_ir(beta, 0.0)
        old = solver._calc_chi0q_general_ir(green, beta)
        spatial_shape = tuple(solver.lattice.shape)

        new = bubble_mod.ir_bubble(green, axF, axB,
                                   spatial_shape=spatial_shape,
                                   scheme="general")
        self.assertEqual(new.shape, old.shape)
        assert_approx_array(new, old, rel=1e-6, abs=1e-8)

    def test_axis_mismatch_raises_value_error(self):
        from hwave.solver.ir_axis import IRAxis
        beta = 2.0
        axF = IRAxis(beta=beta, wmax=5.0, eps=1e-8, statistics="F")
        axB_bad_wmax = IRAxis(beta=beta, wmax=7.0, eps=1e-8, statistics="B")
        axB_ok = IRAxis(beta=beta, wmax=5.0, eps=1e-8, statistics="B")
        green = np.zeros((1, axF.n_freq, 4, 1, 1), dtype=complex)

        with self.subTest(field="wmax"):
            with self.assertRaises(ValueError):
                bubble_mod.ir_bubble(green, axF, axB_bad_wmax,
                                     spatial_shape=(2, 2, 1), scheme="reduced")

        with self.subTest(field="statistics_swapped"):
            with self.assertRaises(ValueError):
                # axF/axB swapped: axF.statistics == "F" required first arg
                bubble_mod.ir_bubble(green, axB_ok, axF,
                                     spatial_shape=(2, 2, 1), scheme="reduced")

        # sanity: matched axes do NOT raise (isolates the mismatch itself
        # as the cause above, not some other shape/scheme issue)
        bubble_mod.ir_bubble(green, axF, axB_ok, spatial_shape=(2, 2, 1),
                             scheme="reduced")


class TestBubbleOldVsNewBondStatic(unittest.TestCase):
    """Old-vs-new gate: ``bond_channels.bond_bubble`` (OLD, sc-layout
    Green) vs ``bubble.bond_bubble_static`` (NEW, canonical Green) on
    identical inputs -- a 2-orbital NN CoulombInter fixture built exactly
    the way ``tests/test_bond_channels.py``'s bubble tests build it
    (``_nn_square_2orb`` / ``_toy_green_4x4_2orb``, imported directly
    rather than re-derived). Layout conversion between the sc.py Green
    layout (``p, p, Nx, Ny, Nz, nmat``) and the kernel's canonical layout
    (``1, nmat, nvol, p, p``) is an INLINE test helper for now -- Task 8
    adds the named ``sc._green_sc_to_canonical`` adapter.

    PLUS the shared-primitive pin (mirroring
    ``tests/test_bond_vs_ed_oracle.py``'s ``TestPin2aRawChi0`` basis): the
    new kernel's ``m=m'=0`` block against ``RPA._calc_chi0q``'s general
    output at ``coeff_tail=0``, reusing this module's own
    ``_make_rpa_solver`` fixture so the pin stays fast (both sides run
    through the identical ``contract_general``/FFT primitives at
    ``Delta r=0``, so round-off, not the ``rel=1e-6``/``abs=1e-8``
    old-vs-new tolerance, is the right bar).
    """

    @staticmethod
    def _sc_to_canonical(green_sc):
        """``(p, p, Nx, Ny, Nz, nmat) -> (1, nmat, Nx*Ny*Nz, p, p)`` --
        the same transpose/reshape/block-axis-insert Task 8's named
        ``sc._green_sc_to_canonical`` adapter will use, kept as an inline
        helper here per this task's brief."""
        p1, p2, nx, ny, nz, nmat = green_sc.shape
        assert p1 == p2
        g = np.transpose(green_sc, (5, 2, 3, 4, 0, 1))  # (nmat,Nx,Ny,Nz,p,p)
        g = np.ascontiguousarray(g).reshape(nmat, nx * ny * nz, p1, p2)
        return g[np.newaxis, ...]  # (1, nmat, nvol, p, p)

    def test_old_vs_new_2orb_nn(self):
        from hwave.solver.bond_channels import bond_bubble, resolve_interactions
        from tests.test_bond_channels import _nn_square_2orb, _toy_green_4x4_2orb

        green_sc, T, N = _toy_green_4x4_2orb(Nx=4, Ny=4, Nz=1, Nmat=8)
        beta = 1.0 / T
        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
        self.assertEqual(bond_set.n_channels, 5)

        old = bond_bubble(green_sc, bond_set, beta=beta)  # (Nx,Ny,Nz,ND,ND)
        nx, ny, nz, nD, nD2 = old.shape
        old_reshaped = old.reshape(nx * ny * nz, nD, nD2)

        green_canonical = self._sc_to_canonical(green_sc)
        new = bubble_mod.bond_bubble_static(
            green_canonical, None, beta, bond_set,
            spatial_shape=(nx, ny, nz))

        self.assertEqual(new.shape, old_reshaped.shape)
        assert_approx_array(new, old_reshaped, rel=1e-6, abs=1e-8)

    def test_shared_primitive_pin_m0_block(self):
        from hwave.solver.bond_channels import resolve_interactions
        from tests.test_bond_channels import _nn_square_2orb

        solver = _make_rpa_solver(coeff_tail=0.0, hz=0.0, norb=2,
                                   complex_hop=True, Lx=4, Ly=4, Nmat=8)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green, tail = solver._calc_green(beta, mu)
        spatial_shape = tuple(solver.lattice.shape)
        nvol = spatial_shape[0] * spatial_shape[1] * spatial_shape[2]

        solver.enable_reduced = False
        old_general = solver._calc_chi0q(green, tail, beta)
        npair = 2 * 2
        old_m0 = old_general[0, solver.nmat // 2].reshape(nvol, npair, npair)

        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
        new_static = bubble_mod.bond_bubble_static(
            green, None, beta, bond_set, spatial_shape=spatial_shape)
        new_m0 = new_static[:, :npair, :npair]

        assert_approx_array(new_m0, old_m0, rel=0, abs=1e-12)

    def test_tail_on_smoke_runs_and_differs_from_tail_off(self):
        """The tail-on bond bubble has no legacy counterpart (legacy
        ``bond_bubble`` never applied the tail correction); Task 5's
        direct-sum oracle carries the correctness burden for tail-on.
        Here only a smoke check: the tail-on path runs and its output
        differs from the tail-off path fed the SAME (deflated) Green
        function -- isolating whether the endpoint-correction machinery
        actually executes."""
        from hwave.solver.bond_channels import resolve_interactions
        from tests.test_bond_channels import _nn_square_2orb

        solver = _make_rpa_solver(coeff_tail=0.5, hz=0.0, norb=2,
                                   complex_hop=True, Lx=4, Ly=4, Nmat=8)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green_deflated, tail = solver._calc_green(beta, mu)
        spatial_shape = tuple(solver.lattice.shape)
        self.assertTrue(bool(np.any(tail != 0)))

        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
        tail_on = bubble_mod.bond_bubble_static(
            green_deflated, tail, beta, bond_set, spatial_shape=spatial_shape)
        tail_off = bubble_mod.bond_bubble_static(
            green_deflated, None, beta, bond_set, spatial_shape=spatial_shape)

        self.assertEqual(tail_on.shape, tail_off.shape)
        max_abs_diff = float(np.max(np.abs(tail_on - tail_off)))
        self.assertGreater(max_abs_diff, 0.0)

    def test_nblock_2_raises_value_error(self):
        from hwave.solver.bond_channels import resolve_interactions
        from tests.test_bond_channels import _nn_square_2orb

        solver = _make_rpa_solver(coeff_tail=0.0, hz=0.3, norb=2,
                                   complex_hop=True, Lx=4, Ly=4, Nmat=8)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green, tail = solver._calc_green(beta, mu)
        self.assertEqual(green.shape[0], 2)
        spatial_shape = tuple(solver.lattice.shape)

        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
        with self.assertRaises(ValueError):
            bubble_mod.bond_bubble_static(
                green, tail, beta, bond_set, spatial_shape=spatial_shape)

    def test_rejects_spatial_shape_mismatch(self):
        from hwave.solver.bond_channels import resolve_interactions
        from tests.test_bond_channels import _nn_square_2orb

        solver = _make_rpa_solver(coeff_tail=0.0, hz=0.0, norb=2,
                                   complex_hop=True, Lx=4, Ly=4, Nmat=8)
        solver._calc_epsilon_k({})
        beta = 1.0 / solver.T
        mu = solver.mu_value
        green, tail = solver._calc_green(beta, mu)
        spatial_shape = tuple(solver.lattice.shape)
        bad_shape = (spatial_shape[0] + 1,) + spatial_shape[1:]

        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3), norb=2)
        with self.assertRaises(ValueError):
            bubble_mod.bond_bubble_static(
                green, tail, beta, bond_set, spatial_shape=bad_shape)


class TestBondDynamicOracle(unittest.TestCase):
    """Independent oracle for ``bubble._iter_bond_dynamic``: a plain-Python
    Matsubara double sum extending ``tests/test_bond_channels.py``'s
    ``_direct_bond_bubble`` (the static, Omega=0-only recipe at line ~122)
    to EVERY bosonic frequency index ``l``.

    Derivation of the frequency extension (recorded here since it is the
    load-bearing new fact this task adds): the fermionic grid is
    ``iw_n = (2n+1-nmat)*pi/beta`` for ``n = 0 .. nmat-1``
    (``_toy_green_4x4``'s convention); the bosonic Matsubara index at
    generator-slot ``l`` is ``2l - nmat`` (``ir_axis.py``'s ``n_l = 2l -
    Nmat`` convention, pinned by the module docstring / spec's "Output
    contracts"). Because the kernel builds G(tau) from only ``nmat``
    sampled frequencies via an FFT, multiplying two such truncated series
    in tau and transforming back realizes a CIRCULAR (mod-``nmat``, no
    wraparound sign flip) convolution over the SAME ``n=0..nmat-1`` window
    -- not literal analytic continuation to frequencies outside the
    window. Concretely, the forward propagator's frequency index shifts by
    ``l - nmat//2`` (wrapped mod ``nmat``) while the reverse propagator's
    index stays unshifted:

        chi[(m,idx),(mp,idxp)](q,l) =
            -(T/N) sum_{k,n} e^{i k.(dr_m - dr_mp)}
                              G_{l1,l3}(k+q, (n + l - nmat//2) mod nmat)
                              G_{l4,l2}(k, n)

    with idx=l1*norb+l2, idx'=l3*norb+l4 (same bond-major/orbital-minor
    convention as ``_direct_bond_bubble``). This formula was mechanically
    derived and cross-checked against ``bubble.dense_bubble``'s ALREADY
    old-vs-new-gated general-scheme output (``TestBubbleOldVsNewDense``)
    at zero displacement, on both a real single-orbital fixture and a
    fully random-Hamiltonian multi-orbital fixture with no built-in
    k -> -k symmetry (which would otherwise make the shift's sign
    unobservable) -- exact match (<1e-15) at every ``(q, l, orbital)``
    combination confirmed the "+", "mod-nmat, no sign flip" convention
    used above (a "-" shift or a wraparound sign flip both fail on the
    asymmetric fixture). The static (``l = nmat//2``, shift 0) case
    reduces to the existing test's unmodified formula.

    Tiny grid (spec's ambiguity resolution): Nx=2, Ny=Nz=1, p=1, nmat=8 --
    reuses ``_toy_green_4x4`` (it accepts arbitrary Nx/Ny/Nz/Nmat) rather
    than re-deriving the single-orbital tight-binding fixture.
    """

    Nx, Ny, Nz, Nmat = 2, 1, 1, 8

    @staticmethod
    def _materialize_dynamic(green_kw, green0_tail, beta, bond_set,
                             spatial_shape):
        """Test-only helper (spec: "Private interfaces" -- kept in the
        test module, not ``bubble.py``): collects ``_iter_bond_dynamic``'s
        streamed ``((m, mp), block)`` pairs into a dense
        ``(nmat, nvol, ND, ND)`` array, for tiny grids only."""
        n_channels = bond_set.n_channels
        npair = None
        out = None
        for (m, mp), block in bubble_mod._iter_bond_dynamic(
                green_kw, green0_tail, beta, bond_set,
                spatial_shape=spatial_shape):
            if out is None:
                nmat, nvol, npair, _ = block.shape
                out = np.zeros((nmat, nvol, npair * n_channels,
                                npair * n_channels), dtype=np.complex128)
            out[:, :, m * npair:(m + 1) * npair,
                mp * npair:(mp + 1) * npair] = block
        return out

    def _direct_reference(self, green, bond_set, T, N):
        """Single-orbital (``l1=l2=l3=l4=0``) direct sum per (m, mp, l, q)
        -- the formula derived in the class docstring, specialized to
        ``norb=1`` (``idx=idx'=0``, so ``chi`` is just
        ``(n_channels, n_channels, Nmat, Nx, Ny, Nz)``, no orbital axes)."""
        Nx, Ny, Nz, Nmat = self.Nx, self.Ny, self.Nz, self.Nmat
        g = green[0, 0]  # (Nx, Ny, Nz, Nmat)
        B = bond_set.n_channels
        kx_list = 2.0 * np.pi * np.arange(Nx) / Nx
        ky_list = 2.0 * np.pi * np.arange(Ny) / Ny
        kz_list = 2.0 * np.pi * np.arange(Nz) / Nz

        ref = np.zeros((B, B, Nmat, Nx, Ny, Nz), dtype=complex)
        for m in range(B):
            drm = bond_set.delta_r[m]
            for mp in range(B):
                drmp = bond_set.delta_r[mp]
                ddx, ddy, ddz = (drm[0] - drmp[0], drm[1] - drmp[1],
                                 drm[2] - drmp[2])
                for l in range(Nmat):
                    shift = l - Nmat // 2
                    for qx in range(Nx):
                        for qy in range(Ny):
                            for qz in range(Nz):
                                s = 0.0 + 0.0j
                                for kx in range(Nx):
                                    for ky in range(Ny):
                                        for kz in range(Nz):
                                            phase = np.exp(1j * (
                                                kx_list[kx] * ddx
                                                + ky_list[ky] * ddy
                                                + kz_list[kz] * ddz))
                                            kqx = (kx + qx) % Nx
                                            kqy = (ky + qy) % Ny
                                            kqz = (kz + qz) % Nz
                                            n_fwd = (np.arange(Nmat)
                                                     + shift) % Nmat
                                            s += phase * np.sum(
                                                g[kqx, kqy, kqz, n_fwd]
                                                * g[kx, ky, kz, :])
                                ref[m, mp, l, qx, qy, qz] = -(T / N) * s
        return ref

    def test_per_omega_matches_direct_sum(self):
        from tests.test_bond_channels import _toy_green_4x4, _nn_square

        green, T, N = _toy_green_4x4(Nx=self.Nx, Ny=self.Ny, Nz=self.Nz,
                                     Nmat=self.Nmat)
        beta = 1.0 / T
        bond_set = resolve_interactions(_nn_square(0.25), np.eye(3), norb=1)

        g_canon = np.ascontiguousarray(
            np.transpose(green[0, 0], (3, 0, 1, 2))
        ).reshape(1, self.Nmat, self.Nx * self.Ny * self.Nz, 1, 1)

        materialized = self._materialize_dynamic(
            g_canon, None, beta, bond_set,
            spatial_shape=(self.Nx, self.Ny, self.Nz))
        ref = self._direct_reference(green, bond_set, T, N)

        for m in range(bond_set.n_channels):
            for mp in range(bond_set.n_channels):
                for l in range(self.Nmat):
                    with self.subTest(m=m, mp=mp, l=l):
                        code_block = materialized[
                            l, :, m:m + 1, mp:mp + 1].reshape(
                                self.Nx, self.Ny, self.Nz)
                        ref_block = ref[m, mp, l]
                        assert_approx_array(code_block, ref_block,
                                            rel=0, abs=1e-10)

    def test_l_half_nmat_pins_static_output(self):
        """The generator's ``l = nmat // 2`` slice must equal
        ``bond_bubble_static``'s output at ``abs=1e-13`` -- same per-pair
        pipeline (:func:`bubble._bond_pair_full_block`), so this is a
        near-exact numerical pin, not an independent-oracle tolerance."""
        from tests.test_bond_channels import _toy_green_4x4, _nn_square

        green, T, N = _toy_green_4x4(Nx=self.Nx, Ny=self.Ny, Nz=self.Nz,
                                     Nmat=self.Nmat)
        beta = 1.0 / T
        bond_set = resolve_interactions(_nn_square(0.25), np.eye(3), norb=1)
        spatial_shape = (self.Nx, self.Ny, self.Nz)

        g_canon = np.ascontiguousarray(
            np.transpose(green[0, 0], (3, 0, 1, 2))
        ).reshape(1, self.Nmat, self.Nx * self.Ny * self.Nz, 1, 1)

        materialized = self._materialize_dynamic(
            g_canon, None, beta, bond_set, spatial_shape=spatial_shape)
        static = bubble_mod.bond_bubble_static(
            g_canon, None, beta, bond_set, spatial_shape=spatial_shape)

        assert_approx_array(materialized[self.Nmat // 2], static,
                            rel=0, abs=1e-13)

    def test_rejects_nblock_2(self):
        from tests.test_bond_channels import _nn_square
        bond_set = resolve_interactions(_nn_square(0.25), np.eye(3), norb=1)
        green = np.zeros((2, self.Nmat, self.Nx * self.Ny * self.Nz, 1, 1),
                         dtype=complex)
        with self.assertRaises(ValueError):
            list(bubble_mod._iter_bond_dynamic(
                green, None, 2.0, bond_set,
                spatial_shape=(self.Nx, self.Ny, self.Nz)))


class TestBondEndpointShift(unittest.TestCase):
    """Nonzero-displacement, asymmetric-Green, tail-ON gate for the
    equal-time endpoint-mean correction's spatial roll (spec: "Endpoint
    correction under a shift" -- "the SAME spatial roll to the reversed
    branch's jump term as to green_rev"; brief: "an m=m'=0 test cannot see
    this error").

    Fixture: a 2-orbital Hamiltonian with different diagonal dispersions
    per orbital (``eps1`` disperses along x and z, ``eps2`` along y and z,
    at different amplitudes) and a constant complex inter-orbital hopping
    (extends the Task 4 ``_toy_green_4x4_2orb`` idea, but built from an
    explicit Hermitian ``H(k)`` diagonalized via ``green.build_green`` so
    both the FULL and DEFLATED Green functions -- and the endpoint jump
    term -- are available with the exact tail semantics the kernel
    consumes). ``coeff_tail=0.5`` (a fractional value, per the spec's note
    that fractional tail coefficients remain in scope for this path).

    The channel pair compared is ``(m, mp) = (1, 4)``, i.e.
    ``Delta r_m=(-1,0,0)``, ``Delta r_mp=(1,0,0)`` (from
    ``_nn_square_2orb``'s ``resolve_interactions`` topology) -- shift
    ``(-2, 0, 0)``, the LARGEST-magnitude nonzero shift available on this
    topology, chosen because a dropped roll on the endpoint jump term is
    most visible where the phase/shift is least trivial.

    Tolerance derivation (measured, not guessed, AT THIS (m, mp) = (1, 4)
    pair -- the truncation error is shift-dependent, so it was remeasured
    for the actual pair/q-points used here rather than reused from a
    different pair): at this fixture's Nmat=128, the tail-ON
    ``bond_bubble_static`` output differs from the direct FULL-Green
    Matsubara sum (restricted to the same finite n=0..Nmat-1 window -- the
    physically exact answer only in the Nmat -> infinity limit, since the
    coeff_tail correction accelerates but does not eliminate the
    finite-window truncation error) by <= 2e-6 in max-abs at both tested
    q-points (measured 1.7e-6 at q=(0,0,0), 8.6e-7 at q=(1,1,0)), scaling
    ~1/Nmat**2 (measured at Nmat=16/32/64/128/256:
    1.0e-4/2.7e-5/6.8e-6/1.7e-6/4.3e-7, each ~4x the next -- consistent
    with the tail correction's improved but still finite-order
    convergence). The gate uses abs=1e-5, a >5x margin over the measured
    error, which is still several orders of magnitude tighter than the
    scale at which a dropped roll on the jump term shows up (an
    undropped-vs-dropped roll differs by O(1e-2) at this shift -- see the
    mutation-check note in the task report).
    """

    Nx, Ny, Nz, Nmat = 4, 2, 1, 128
    BETA = 1.3
    MU = 0.15
    COEFF_TAIL = 0.5

    @classmethod
    def _build_fixture(cls):
        t = 1.0
        mu1, mu2 = 0.2, -0.35
        t12 = 0.18 + 0.07j
        kx = 2.0 * np.pi * np.arange(cls.Nx) / cls.Nx
        ky = 2.0 * np.pi * np.arange(cls.Ny) / cls.Ny
        kz = 2.0 * np.pi * np.arange(cls.Nz) / cls.Nz
        KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
        eps1 = -2.0 * t * np.cos(KX) - 1.0 * t * np.cos(KZ) - mu1
        eps2 = -2.0 * t * np.cos(KY) - 0.5 * t * np.cos(KZ) - mu2
        nvol = cls.Nx * cls.Ny * cls.Nz
        H = np.zeros((nvol, 2, 2), dtype=complex)
        H[:, 0, 0] = eps1.reshape(nvol)
        H[:, 1, 1] = eps2.reshape(nvol)
        H[:, 0, 1] = t12
        H[:, 1, 0] = np.conj(t12)
        ev, V = np.linalg.eigh(H)
        full, defl, tail = green_mod.build_green(
            ev, V, cls.MU, cls.BETA, cls.Nmat, cls.COEFF_TAIL)
        return full, defl, tail

    @staticmethod
    def _direct_block_full_green(full, bond_set, T, N, spatial_shape, Nmat,
                                 m, mp, qx, qy, qz):
        """Static (Omega=0) direct Matsubara sum over the FULL Green
        function for ONE (m, mp, q) cell, all 2x2 orbital indices --
        restricted (not the whole (m, mp, q) grid) to keep this fast."""
        Nx, Ny, Nz = spatial_shape
        drm = bond_set.delta_r[m]
        drmp = bond_set.delta_r[mp]
        ddx, ddy, ddz = (drm[0] - drmp[0], drm[1] - drmp[1],
                         drm[2] - drmp[2])
        kx_list = 2.0 * np.pi * np.arange(Nx) / Nx
        ky_list = 2.0 * np.pi * np.arange(Ny) / Ny
        kz_list = 2.0 * np.pi * np.arange(Nz) / Nz

        def g_at(l1, l3, kx, ky, kz, n):
            kidx = (kx * Ny + ky) * Nz + kz
            return full[0, n, kidx, l1, l3]

        npair = 4
        block = np.zeros((npair, npair), dtype=complex)
        for l1 in range(2):
            for l2 in range(2):
                for l3 in range(2):
                    for l4 in range(2):
                        s = 0.0 + 0.0j
                        for kx in range(Nx):
                            for ky in range(Ny):
                                for kz in range(Nz):
                                    phase = np.exp(1j * (
                                        kx_list[kx] * ddx
                                        + ky_list[ky] * ddy
                                        + kz_list[kz] * ddz))
                                    kqx = (kx + qx) % Nx
                                    kqy = (ky + qy) % Ny
                                    kqz = (kz + qz) % Nz
                                    for n in range(Nmat):
                                        s += phase * g_at(
                                            l1, l3, kqx, kqy, kqz, n
                                        ) * g_at(l4, l2, kx, ky, kz, n)
                        block[l1 * 2 + l2, l3 * 2 + l4] = -(T / N) * s
        return block

    def test_tail_on_nonzero_shift_matches_direct_sum(self):
        from tests.test_bond_channels import _nn_square_2orb

        full, defl, tail = self._build_fixture()
        spatial_shape = (self.Nx, self.Ny, self.Nz)
        bond_set = resolve_interactions(_nn_square_2orb(), np.eye(3),
                                         norb=2)
        self.assertEqual(bond_set.delta_r[1], (-1, 0, 0))
        self.assertEqual(bond_set.delta_r[4], (1, 0, 0))
        m, mp = 1, 4

        static_tail_on = bubble_mod.bond_bubble_static(
            defl, tail, self.BETA, bond_set, spatial_shape=spatial_shape)

        T = 1.0 / self.BETA
        N = self.Nx * self.Ny * self.Nz
        npair = 4

        for (qx, qy, qz) in [(0, 0, 0), (1, 1, 0)]:
            with self.subTest(q=(qx, qy, qz)):
                kidx = (qx * self.Ny + qy) * self.Nz + qz
                code_block = static_tail_on[
                    kidx, m * npair:(m + 1) * npair,
                    mp * npair:(mp + 1) * npair]
                ref_block = self._direct_block_full_green(
                    full, bond_set, T, N, spatial_shape, self.Nmat,
                    m, mp, qx, qy, qz)
                assert_approx_array(code_block, ref_block, rel=0, abs=1e-5)


if __name__ == "__main__":
    unittest.main()
