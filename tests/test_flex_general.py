#!/usr/bin/env python3

"""Tests for the general-mode (paramagnetic full-vertex) FLEX guards.

These tests cover Task 1: ``calc_scheme="general"`` must be accepted for
paramagnetic (spin-free) FLEX and rejected (fail-fast ValueError) for
``enable_spin_orbital`` (at construction) and for non-spin-free spin modes
(spin-diag / spinful, checked inside ``solve()``).

Tests must be run from the repository root (they use relative paths like
``tests/rpa/input``).
"""

import contextlib
import logging
import os
import unittest
import numpy as np


@contextlib.contextmanager
def _assert_no_warning(testcase, logger_name):
    """Python 3.9-compatible replacement for assertNoLogs (added in 3.10):
    capture records on ``logger_name`` and assert none are WARNING or above."""
    records = []
    handler = logging.Handler()
    handler.emit = records.append
    lg = logging.getLogger(logger_name)
    lg.addHandler(handler)
    try:
        yield
    finally:
        lg.removeHandler(handler)
    warns = [r.getMessage() for r in records if r.levelno >= logging.WARNING]
    testcase.assertEqual(warns, [], "unexpected warnings: {}".format(warns))


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


def _write_2d_2orb_onsite_fixture(dirpath):
    """Materialize a [4,4,1]-compatible 2-orbital fixture with ON-SITE U'.

    Mirrors ``tests/rpa/input_2orb`` (same 2D geom + transfer) but writes an
    ON-SITE (irvec ``0 0 0``) inter-orbital CoulombInter, unlike the committed
    fixture whose CoulombInter is off-site.  The general (full-vertex) path now
    fail-fasts on off-site two-body terms (see ``_inflate_chi0q_and_ham_general``),
    so general-path tests must use on-site interactions.  On-site inter-orbital
    U' still gives nontrivial 4x4 MYO S/C blocks.
    """
    os.makedirs(dirpath, exist_ok=True)
    geom = ("  1.000000000000   0.000000000000   0.000000000000\n"
            "  0.000000000000   1.000000000000   0.000000000000\n"
            "  0.000000000000   0.000000000000   1.000000000000\n"
            "2\n"
            "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
            "    0.500000000000000e+00     0.500000000000000e+00     0.000000000000000e+00\n")
    transfer = ("Transfer in wannier90-like format for uhfk\n"
                "2\n"
                "9\n"
                " 1 1 1 1 1 1 1 1 1\n"
                "   0    1    0    1    1  1.0 0.0\n"
                "   0   -1    0    1    1  1.0 0.0\n"
                "   0    1    0    2    2  1.0 0.0\n"
                "   0   -1    0    2    2  1.0 0.0\n"
                "   0    0    0    1    2  0.5 0.0\n"
                "  -1    0    0    1    2  0.5 0.0\n"
                "   0    0    0    2    1  0.5 0.0\n"
                "   1    0    0    2    1  0.5 0.0\n")
    coulombintra = ("CoulombIntra in wannier90-like format for uhfk\n"
                    "2\n"
                    "1\n"
                    " 1\n"
                    "   0    0    0    1    1   4.000000000000   0.000000000000\n"
                    "   0    0    0    2    2   4.000000000000   0.000000000000\n")
    coulombinter = ("CoulombInter in wannier90-like format for uhfk\n"
                    "2\n"
                    "1\n"
                    " 1\n"
                    "   0    0    0    1    2   1.000000000000   0.000000000000\n"
                    "   0    0    0    2    1   1.000000000000   0.000000000000\n")
    with open(os.path.join(dirpath, "geom.dat"), "w") as f:
        f.write(geom)
    with open(os.path.join(dirpath, "transfer.dat"), "w") as f:
        f.write(transfer)
    with open(os.path.join(dirpath, "coulombintra.dat"), "w") as f:
        f.write(coulombintra)
    with open(os.path.join(dirpath, "coulombinter.dat"), "w") as f:
        f.write(coulombinter)


def _make_general_flex(norb=2):
    """Build a spin-free general FLEX solver.

    For ``norb=2`` (default) uses a self-contained 2-orbital fixture with
    CoulombIntra + ON-SITE CoulombInter (see ``_write_2d_2orb_onsite_fixture``)
    so the MYO S/C matrices are nontrivial 4x4 blocks.  (The committed
    ``tests/rpa/input_2orb`` fixture has OFF-SITE CoulombInter, which the general
    full-vertex path now rejects.)  For ``norb=1`` uses the 1-orbital
    ``tests/rpa/input`` fixture (CoulombIntra only).  The constructor sets
    ``self.norb``, ``self.lattice`` and ``self.ham_info`` (with
    ``self.ham_info.param_ham``), which is all the downstream general-path
    methods need.
    """
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex

    if norb == 1:
        info_input = {
            'path_to_input': 'tests/rpa/input',
            'interaction': {
                'path_to_input': 'tests/rpa/input',
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        }
    else:
        dirpath = os.path.join(tempfile.gettempdir(),
                               "hwave_flex_2d_2orb_onsite")
        _write_2d_2orb_onsite_fixture(dirpath)
        info_input = {
            'path_to_input': dirpath,
            'interaction': {
                'path_to_input': dirpath,
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

        # chi0q is converted RPA->MYO via an orbital-pair transpose; the shape
        # is unchanged (all four orbital legs are norb).
        nvol = flex.lattice.nvol
        self.assertEqual(chi0q.ndim, 6)
        self.assertEqual(chi0q.shape,
                         (flex.nmat, nvol, no, no, no, no))
        np.testing.assert_array_equal(
            chi0q, chi0q_raw.transpose(0, 1, 4, 5, 2, 3))

        # S/C matrices reshaped to (nvol, norb^2, norb^2).
        self.assertEqual(Us.shape, (nvol, no * no, no * no))
        self.assertEqual(Uc.shape, (nvol, no * no, no * no))
        self.assertEqual(Us.shape[-2:], (4, 4))
        self.assertEqual(Uc.shape[-2:], (4, 4))

        # CoulombIntra + CoulombInter present -> nonzero S/C.
        self.assertTrue(np.any(Us != 0.0))
        self.assertTrue(np.any(Uc != 0.0))


def _rpa_general_chi(chi0q, ham, sign):
    """Reference matrix-RPA susceptibility, computed by direct reshape+solve.

    Computes ``chi = [I + sign * chi0q.ham]^{-1} . chi0q`` independently of
    ``_solve_channels_general`` / ``_solve_rpa`` (no block detection, no
    threading) so it is a genuine cross-check.
    """
    nmat, nvol = chi0q.shape[0], chi0q.shape[1]
    no = chi0q.shape[2]
    ndx = no * no
    chi0_2d = chi0q.reshape(nmat, nvol, ndx, ndx)
    ham_2d = ham.reshape(nvol, ndx, ndx)
    eye = np.eye(ndx, dtype=complex)
    mat = eye + sign * (chi0_2d @ ham_2d[np.newaxis])
    sol = np.linalg.solve(mat, chi0_2d)
    return sol.reshape(chi0q.shape)


class TestChiGeneralConsistency(unittest.TestCase):
    """Task 5: general-path chi_s/chi_c equal direct matrix-RPA results."""

    def test_matches_rpa_general(self):
        flex = _make_general_flex()
        chi0, Us, Uc = flex._inflate_chi0q_and_ham_general(
            _fake_general_chi0q(flex), None)
        chi_s, chi_c = flex._solve_channels_general(chi0, Us, Uc)
        chi_s_ref = _rpa_general_chi(chi0, Us, sign=-1)   # [I - chi0 Us]^-1 chi0
        chi_c_ref = _rpa_general_chi(chi0, Uc, sign=+1)   # [I + chi0 Uc]^-1 chi0
        np.testing.assert_allclose(chi_s, chi_s_ref, atol=1e-10)
        np.testing.assert_allclose(chi_c, chi_c_ref, atol=1e-10)


class TestGeneralOutputConvention(unittest.TestCase):
    """The general path computes internally in MYO convention but must expose
    the public chi0q output in the same RPA [a,c,b,d] orbital convention as the
    reduced path, so the saved chi0q key has one consistent meaning."""

    def test_output_chi0q_is_rpa_convention(self):
        flex = _make_general_flex(norb=2)
        flex._myo_sc_cache = None
        chi0_raw = _fake_general_chi0q(flex)
        chi0q_out, v_eff, chi_s, chi_c = \
            flex._flex_compute_veff_general(chi0_raw, None)
        # output chi0q must equal the [a,c,b,d] input (not the MYO transpose)
        np.testing.assert_allclose(
            chi0q_out, chi0_raw, atol=1e-12,
            err_msg="general output chi0q must stay in RPA [a,c,b,d] convention")
        # internal MYO chi0q (the transpose) must differ for a non-symmetric
        # tensor -- confirms the output is genuinely back-converted, not MYO.
        myo = chi0_raw.transpose(0, 1, 4, 5, 2, 3)
        self.assertGreater(np.linalg.norm(chi0q_out - myo), 1e-6)


class TestGeneralChiConventionRoundTrip(unittest.TestCase):
    """General FLEX tags saved chiq_s/chiq_c as MYO; the Eliashberg loader reads
    that tag back so the pairing vertex uses MYO (not Kuroki) S/C matrices."""

    def _save_and_reload_convention(self, flex):
        import tempfile
        import hwave.sc as sc
        d = tempfile.mkdtemp()
        nd = flex.norb * flex.norb
        shape = (flex.nmat, flex.lattice.nvol, nd, nd)
        green_info = {
            "chiq_s": np.zeros(shape, dtype=complex),
            "chiq_c": np.zeros(shape, dtype=complex),
        }
        info_out = {"path_to_output": d, "chiq_s": "chiq_s", "chiq_c": "chiq_c"}
        flex.save_results(info_out, green_info)
        input_dict = {"file": {"output": {"path_to_output": d}},
                      "eliashberg": {}}
        nx, ny, nz = flex.lattice.shape
        _, _, _, conv = sc._load_flex_susceptibilities(
            input_dict, flex.norb, nx, ny, nz)
        return conv

    def test_general_tags_and_reloads_myo(self):
        flex = _make_general_flex(norb=2)
        self.assertEqual(self._save_and_reload_convention(flex), "myo")

    def test_legacy_file_without_tag_defaults_kuroki(self):
        import tempfile
        import hwave.sc as sc
        d = tempfile.mkdtemp()
        # write a legacy chiq file with NO chi_convention field
        nd = 4
        arr = np.zeros((4, 4, nd, nd), dtype=complex)
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=arr)
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=arr)
        input_dict = {"file": {"output": {"path_to_output": d}},
                      "eliashberg": {}}
        _, _, _, conv = sc._load_flex_susceptibilities(input_dict, 2, 2, 2, 1)
        self.assertEqual(conv, "kuroki")

    def test_mismatched_s_c_conventions_raise(self):
        import tempfile
        import hwave.sc as sc
        d = tempfile.mkdtemp()
        nd = 4
        arr = np.zeros((4, 4, nd, nd), dtype=complex)
        # spin says myo, charge says kuroki -> must not silently combine
        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=arr, chi_convention="myo")
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=arr,
                 chi_convention="kuroki")
        input_dict = {"file": {"output": {"path_to_output": d}},
                      "eliashberg": {}}
        with self.assertRaises(ValueError):
            sc._load_flex_susceptibilities(input_dict, 2, 2, 2, 1)

    def test_real_rank6_chi_roundtrips_to_myo_flattening(self):
        """Save REAL general-path (rank-6 MYO) chi_s/chi_c, reload through the
        Eliashberg loader, and confirm the static slice flattens to (nd,nd) in
        the MYO row=(m,n)/col=(mu,nu) order that build_sc_matrices_myo expects."""
        import tempfile
        import hwave.sc as sc
        flex = _make_general_flex(norb=2)
        flex._myo_sc_cache = None
        chi0 = _fake_general_chi0q(flex)
        _, _, chi_s, chi_c = flex._flex_compute_veff_general(chi0, None)
        # real general chi is rank-6 (nmat, nvol, m, n, mu, nu)
        self.assertEqual(chi_s.ndim, 6)
        nx, ny, nz = flex.lattice.shape
        nd = flex.norb * flex.norb

        d = tempfile.mkdtemp()
        flex.save_results(
            {"path_to_output": d, "chiq_s": "chiq_s", "chiq_c": "chiq_c"},
            {"chiq_s": chi_s, "chiq_c": chi_c})
        input_dict = {"file": {"output": {"path_to_output": d}},
                      "eliashberg": {}}
        chis, chic, _, conv = sc._load_flex_susceptibilities(
            input_dict, flex.norb, nx, ny, nz)
        self.assertEqual(conv, "myo")
        self.assertEqual(chis.shape, (nx, ny, nz, nd, nd))
        # expected static slice flattened with row=(m,n), col=(mu,nu)
        center = flex.nmat // 2
        expected = chi_s[center].reshape(nx, ny, nz, nd, nd)
        np.testing.assert_allclose(chis, expected, atol=1e-12)


class TestVeffGeneral(unittest.TestCase):
    """Task 6: MYO fluctuation effective interaction V_eff."""

    def test_single_orbital_matches_reduced_kernel(self):
        U = 3.0
        flex = _make_general_flex(norb=1)
        chi0 = _fake_general_chi0q(flex)        # (nmat, nvol, 1,1,1,1)
        nvol = flex.lattice.nvol
        Us = U * np.ones((nvol, 1, 1), dtype=complex)
        Uc = U * np.ones((nvol, 1, 1), dtype=complex)
        chi_s, chi_c = flex._solve_channels_general(chi0, Us, Uc)
        V = flex._calc_veff_general(chi0, chi_s, chi_c, Us, Uc)
        # reduced fluctuation kernel sandwiched by U: U*(1.5 chi_s + 0.5 chi_c - chi0)*U
        kernel = 1.5 * chi_s + 0.5 * chi_c - chi0          # (nmat,nvol,1,1,1,1)
        Vref = (U * U) * kernel.reshape(V.shape)
        np.testing.assert_allclose(V, Vref, atol=1e-10)

    def test_shape_two_orbital(self):
        flex = _make_general_flex(norb=2)
        chi0, Us, Uc = flex._inflate_chi0q_and_ham_general(
            _fake_general_chi0q(flex), None)
        chi_s, chi_c = flex._solve_channels_general(chi0, Us, Uc)
        V = flex._calc_veff_general(chi0, chi_s, chi_c, Us, Uc)
        self.assertEqual(V.shape, (flex.nmat, flex.lattice.nvol, 4, 4))

    def test_values_two_orbital_vs_elementwise(self):
        rng = np.random.default_rng(99)
        flex = _make_general_flex(norb=2)
        chi0 = _fake_general_chi0q(flex)
        nvol, no = flex.lattice.nvol, flex.norb
        ndx = no * no
        Us = (rng.standard_normal((nvol, ndx, ndx))
              + 1j * rng.standard_normal((nvol, ndx, ndx)))
        Uc = (rng.standard_normal((nvol, ndx, ndx))
              + 1j * rng.standard_normal((nvol, ndx, ndx)))
        chi_s, chi_c = flex._solve_channels_general(chi0, Us, Uc)
        V = flex._calc_veff_general(chi0, chi_s, chi_c, Us, Uc)
        chi0_2d = chi0.reshape(flex.nmat, nvol, ndx, ndx)
        chis_2d = chi_s.reshape(flex.nmat, nvol, ndx, ndx)
        chic_2d = chi_c.reshape(flex.nmat, nvol, ndx, ndx)
        Vref = np.zeros_like(V)
        for f in range(flex.nmat):
            for r in range(nvol):
                Vref[f, r] = (1.5 * Us[r] @ chis_2d[f, r] @ Us[r]
                              + 0.5 * Uc[r] @ chic_2d[f, r] @ Uc[r]
                              - 0.25 * (Us[r] + Uc[r]) @ chi0_2d[f, r] @ (Us[r] + Uc[r]))
        np.testing.assert_allclose(V, Vref, atol=1e-10)


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

    def test_general_rejects_offsite_interaction(self):
        """The general (v1) path supports on-site interactions only; an off-site
        interaction entry (irvec != (0,0,0)) must fail-fast with a ValueError
        when the MYO S/C matrices are built in _inflate_chi0q_and_ham_general."""
        flex = _make_general_flex(norb=2)
        pham = flex.ham_info.param_ham
        # param_ham[itype] is a dict {(irvec, orbvec): value}; inject an off-site
        # CoulombInter entry mirroring that real key structure.  irvec=(1,0,0) is
        # off-site; orbvec=(0,1) is a valid orbital pair.
        key = ((1, 0, 0), (0, 1))
        pham.setdefault("CoulombInter", {})[key] = 1.0
        # ensure a clean cache so the guard (cache-MISS branch) actually runs:
        flex._myo_sc_cache = None
        chi0_raw = _fake_general_chi0q(flex)
        with self.assertRaises(ValueError) as cm:
            flex._inflate_chi0q_and_ham_general(chi0_raw, None)
        # confirm we hit the off-site guard specifically (not some other error)
        self.assertIn("off-site", str(cm.exception).lower())

    def test_general_pairlift_is_inert_and_warns(self):
        """PairLift contributes S=C=0 to the particle-hole spin/charge vertex
        (cf. hwave.sc), so the general path correctly omits it. It must WARN
        (not silently drop, and not error) and the result must be numerically
        identical to the no-PairLift case."""
        # baseline (no PairLift)
        flex0 = _make_general_flex(norb=2)
        flex0._myo_sc_cache = None
        chi0_raw = _fake_general_chi0q(flex0)
        _, base_v, base_s, base_c = \
            flex0._flex_compute_veff_general(chi0_raw, None)

        # with on-site PairLift injected
        flex1 = _make_general_flex(norb=2)
        pham = flex1.ham_info.param_ham
        for (a, b) in ((0, 1), (1, 0)):
            pham.setdefault("PairLift", {})[((0, 0, 0), (a, b))] = 0.3
        flex1._myo_sc_cache = None
        with self.assertLogs('hwave.solver.flex', level='WARNING') as cm:
            _, p_v, p_s, p_c = flex1._flex_compute_veff_general(chi0_raw, None)
        self.assertTrue(any("pairlift" in m.lower() for m in cm.output))
        # inert: identical S/C/V_eff with and without PairLift
        np.testing.assert_allclose(p_v, base_v, atol=1e-12)
        np.testing.assert_allclose(p_s, base_s, atol=1e-12)
        np.testing.assert_allclose(p_c, base_c, atol=1e-12)

    def test_general_iterationmax_zero_no_nameerror(self):
        """IterationMax=0 must not raise NameError in solve(): the SCF loop body
        (which assigns ``diff``) never runs, and the 'did not converge' warning
        branch references ``diff`` — it must be initialized before the loop."""
        flex = _make_general_flex(norb=1)
        flex.max_iter = 0
        green_info = {}
        os.makedirs('tests/flex/output', exist_ok=True)
        # Should warn-and-return (no convergence) without NameError on ``diff``.
        flex.solve(green_info, 'tests/flex/output')


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
        with _assert_no_warning(self, 'hwave.solver.flex'):
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


class TestSelfEnergyGeneral(unittest.TestCase):
    """Task 7: rank-4 orbital self-energy contraction.

    Validates the ONLY new, bug-prone code (the rank-4 orbital contraction
    Sigma_{mn} = sum_{mu,nu} V_{mu m, nu n} G_{mu nu}) in isolation against an
    independent physical-index ground truth.  The frequency/momentum FFT
    transport is reused unchanged from the already-tested reduced path, so it
    is NOT compared here (the brute-force Matsubara wrap is a toy that does not
    match the real fermionic/bosonic FFT phases).
    """

    def test_orbital_contraction_matches_physical(self):
        rng = np.random.default_rng(11)
        norb = 2
        flex = _make_general_flex(norb=2)   # only needed so the method is bound
        # Physical vertex and Green function at a single (r,tau) slice:
        Vphys = (rng.standard_normal((norb, norb, norb, norb))
                 + 1j * rng.standard_normal((norb, norb, norb, norb)))  # (mu,m,nu,n)
        G = rng.standard_normal((norb, norb)) + 1j * rng.standard_normal((norb, norb))  # (mu,nu)
        # Production input: flatten physical V to (1,1,norb^2,norb^2).
        v_rt = Vphys.reshape(1, 1, norb * norb, norb * norb)
        green_rt = G.reshape(1, 1, 1, norb, norb)
        sig = flex._sigma_orbital_contract(v_rt, green_rt)   # (1,1,1,norb,norb)
        # Independent reference: explicit physical loop, indexing Vphys directly.
        sig_ref = np.zeros((norb, norb), dtype=complex)
        for m in range(norb):
            for n in range(norb):
                s = 0j
                for mu in range(norb):
                    for nu in range(norb):
                        s += Vphys[mu, m, nu, n] * G[mu, nu]
                sig_ref[m, n] = s
        np.testing.assert_allclose(sig[0, 0, 0], sig_ref, atol=1e-12)

    def test_orbital_contraction_matches_bruteforce_single_point(self):
        # Cross-check against the committed physical ground truth at Nk=1,nmat=1
        # (pure orbital contraction; T=1 so the (T/Nk) prefactor is 1).
        from tests.flex_bruteforce_ref import sigma_bruteforce
        rng = np.random.default_rng(12)
        norb = 2
        flex = _make_general_flex(norb=2)
        Vphys = (rng.standard_normal((norb, norb, norb, norb))
                 + 1j * rng.standard_normal((norb, norb, norb, norb)))
        G = rng.standard_normal((norb, norb)) + 1j * rng.standard_normal((norb, norb))
        v_rt = Vphys.reshape(1, 1, norb * norb, norb * norb)
        green_rt = G.reshape(1, 1, 1, norb, norb)
        sig = flex._sigma_orbital_contract(v_rt, green_rt)[0, 0, 0]
        # brute-force at one k, one freq:
        Vbf = Vphys.reshape(norb, norb, norb, norb, 1, 1)   # (mu,m,nu,n,q=1,iv=1)
        Gbf = G.reshape(norb, norb, 1, 1)                   # (mu,nu,Nk=1,nmat=1)
        sig_bf = sigma_bruteforce(Gbf, Vbf, T=1.0, Nk=1)[:, :, 0, 0]   # (m,n)
        np.testing.assert_allclose(sig, sig_bf, atol=1e-12)

    def test_full_method_shape(self):
        flex = _make_general_flex(norb=2)
        nblock, nmat, nvol, no = 1, flex.nmat, flex.lattice.nvol, flex.norb
        rng = np.random.default_rng(13)
        green_kw = (rng.standard_normal((nblock, nmat, nvol, no, no))
                    + 1j * rng.standard_normal((nblock, nmat, nvol, no, no)))
        v_eff = (rng.standard_normal((nmat, nvol, no * no, no * no))
                 + 1j * rng.standard_normal((nmat, nvol, no * no, no * no)))
        sig = flex._calc_self_energy_general(green_kw, v_eff, beta=1.0 / flex.T)
        self.assertEqual(sig.shape, (nblock, nmat, nvol, no, no))


class TestGeneralSolveEndToEnd(unittest.TestCase):
    """Task 8: the general FLEX path runs end-to-end through solve().

    Builds a tiny 2-orbital general FLEX (small grid, small Nmat, small
    IterationMax) and runs solve().  Asserts the output arrays are written and
    that sigma is finite (no NaN/Inf) -- a basic sanity that the general path
    actually ran.  Without the solve() dispatch the reduced helper would receive
    the rank-6 general chi0q and either crash on the reduced reshape or produce a
    wrong-shaped sigma; the sigma-shape assertion makes the test depend on the
    dispatch even if no exception is raised.
    """

    def test_runs_and_saves(self):
        import tempfile
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex
        dirpath = os.path.join(tempfile.gettempdir(),
                               "hwave_flex_2d_2orb_onsite")
        _write_2d_2orb_onsite_fixture(dirpath)
        info_input = {
            'path_to_input': dirpath,
            'interaction': {
                'path_to_input': dirpath,
                'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            },
        }
        ham = read_input_k.QLMSkInput(info_input).get_param('ham')
        green = read_input_k.QLMSkInput(info_input).get_param('green')
        info_mode = {
            'mode': 'FLEX',
            'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                      'SubShape': [1, 1, 1], 'Nmat': 8,
                      'IterationMax': 3, 'Mix': 0.5},
            'calc_scheme': 'general',
        }
        solver = solver_flex.FLEX(ham, {}, info_mode)
        self.assertTrue(solver._flex_general)

        green_info = green
        os.makedirs('tests/flex/output', exist_ok=True)
        solver.solve(green_info, 'tests/flex/output')

        for key in ('sigma', 'green', 'chiq_s', 'chiq_c', 'chi0q'):
            self.assertIn(key, green_info)

        # sigma must be finite (no NaN/Inf) -- a basic sanity the path ran.
        self.assertTrue(np.all(np.isfinite(green_info['sigma'])))

        # sigma shape must be the spin-free general shape
        # (nblock=1, Nmat, nvol, norb, norb): this depends on the general
        # dispatch (the reduced path would not produce this for general chi0q).
        nvol = solver.lattice.nvol
        norb = solver.norb
        self.assertEqual(green_info['sigma'].shape,
                         (1, solver.nmat, nvol, norb, norb))

    def test_runs_with_hund_and_exchange(self):
        """End-to-end general solve() with on-site Hund and Exchange present:
        the full-vertex path must build the MYO S/C matrices from these terms
        and converge to a finite sigma (the off-diagonal vertices the reduced
        path would drop are exactly what calc_scheme='general' keeps)."""
        import tempfile
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex
        dirpath = os.path.join(tempfile.gettempdir(),
                               "hwave_flex_2d_2orb_onsite")
        _write_2d_2orb_onsite_fixture(dirpath)
        info_input = {
            'path_to_input': dirpath,
            'interaction': {
                'path_to_input': dirpath,
                'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
                'CoulombInter': 'coulombinter.dat',
            },
        }
        info_mode = {
            'mode': 'FLEX',
            'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                      'SubShape': [1, 1, 1], 'Nmat': 8,
                      'IterationMax': 3, 'Mix': 0.5},
            'calc_scheme': 'general',
        }

        def _solve(with_J):
            # read a FRESH ham/green per solve so injecting J into param_ham does
            # not leak into the other (baseline) run via a shared ham object.
            ham = read_input_k.QLMSkInput(info_input).get_param('ham')
            green = read_input_k.QLMSkInput(info_input).get_param('green')
            s = solver_flex.FLEX(ham, {}, info_mode)
            if with_J:
                # inject on-site Hund and Exchange (orbital off-diagonal) so the
                # general S/C builder uses J / J' -- the terms the reduced path
                # would drop.
                pham = s.ham_info.param_ham
                for (a, b) in ((0, 1), (1, 0)):
                    pham.setdefault("Hund", {})[((0, 0, 0), (a, b))] = 0.2
                    pham.setdefault("Exchange", {})[((0, 0, 0), (a, b))] = 0.1
            s._myo_sc_cache = None
            os.makedirs('tests/flex/output', exist_ok=True)
            s.solve(green, 'tests/flex/output')
            return s, green['sigma']

        solver, sigma_J = _solve(with_J=True)
        _, sigma_noJ = _solve(with_J=False)

        self.assertTrue(np.all(np.isfinite(sigma_J)))
        self.assertEqual(sigma_J.shape,
                         (1, solver.nmat, solver.lattice.nvol,
                          solver.norb, solver.norb))
        # Hund/Exchange must actually enter the MYO S/C vertex: a bug that drops
        # them would give sigma identical to the no-J baseline.
        self.assertGreater(np.linalg.norm(sigma_J - sigma_noJ), 1e-8,
                           "Hund/Exchange must change the general-path self-energy")


def _write_1d_2orb_fixture(dirpath):
    """Materialize a 1D-compatible (x-only hoppings) 2-orbital input fixture.

    The committed ``tests/rpa/input_2orb`` fixture hops in the y direction
    (irvec ``0 +-1 0``), so it is incompatible with CellShape ``[4,1,1]``.  Here
    we write a self-contained fixture with only x-direction hoppings so a 1D
    ``[4,1,1]`` grid (Nk=4) is valid -- giving the small, brute-force-friendly
    problem the convention/pipeline cross-checks need.
    """
    os.makedirs(dirpath, exist_ok=True)
    geom = ("  1.000000000000   0.000000000000   0.000000000000\n"
            "  0.000000000000   1.000000000000   0.000000000000\n"
            "  0.000000000000   0.000000000000   1.000000000000\n"
            "2\n"
            "    0.000000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n"
            "    0.500000000000000e+00     0.000000000000000e+00     0.000000000000000e+00\n")
    # x-only transfer: intra-orbital +-x hops + inter-orbital +-x hops.
    transfer = ("Transfer (x-only) in wannier90-like format for uhfk\n"
                "2\n"
                "8\n"
                " 1 1 1 1 1 1 1 1\n"
                "   1    0    0    1    1  1.0 0.0\n"
                "  -1    0    0    1    1  1.0 0.0\n"
                "   1    0    0    2    2  1.0 0.0\n"
                "  -1    0    0    2    2  1.0 0.0\n"
                "   0    0    0    1    2  0.5 0.0\n"
                "  -1    0    0    1    2  0.5 0.0\n"
                "   0    0    0    2    1  0.5 0.0\n"
                "   1    0    0    2    1  0.5 0.0\n")
    coulombintra = ("CoulombIntra in wannier90-like format for uhfk\n"
                    "2\n"
                    "1\n"
                    " 1\n"
                    "   0    0    0    1    1   4.000000000000   0.000000000000\n"
                    "   0    0    0    2    2   4.000000000000   0.000000000000\n")
    # on-site (irvec 0,0,0) inter-orbital density-density coupling.
    coulombinter = ("CoulombInter in wannier90-like format for uhfk\n"
                    "2\n"
                    "1\n"
                    " 1\n"
                    "   0    0    0    1    2   2.000000000000   0.000000000000\n"
                    "   0    0    0    2    1   2.000000000000   0.000000000000\n")
    with open(os.path.join(dirpath, "geom.dat"), "w") as f:
        f.write(geom)
    with open(os.path.join(dirpath, "transfer.dat"), "w") as f:
        f.write(transfer)
    with open(os.path.join(dirpath, "coulombintra.dat"), "w") as f:
        f.write(coulombintra)
    with open(os.path.join(dirpath, "coulombinter.dat"), "w") as f:
        f.write(coulombinter)


def _make_chi0_general_flex():
    """Build a spin-free 2-orbital general FLEX with CellShape [4,1,1], Nmat=1.

    Used by the chi0-convention and full-pipeline physical cross-checks below.
    Uses a self-contained x-only 2-orbital fixture (see
    ``_write_1d_2orb_fixture``) so the 1D ``[4,1,1]`` (Nk=4) grid is valid.
    """
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    dirpath = os.path.join(tempfile.gettempdir(), "hwave_flex_1d_2orb")
    _write_1d_2orb_fixture(dirpath)
    info_input = {
        'path_to_input': dirpath,
        'interaction': {
            'path_to_input': dirpath,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
        },
    }
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    # Nmat must be even at construction (validator rejects odd / zero); we drop
    # to a single Matsubara point AFTER construction for the nmat=1 cross-check
    # (the brute-force reference works at one frequency point).
    info_mode = {
        'mode': 'FLEX',
        'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 1, 1],
                  'SubShape': [1, 1, 1], 'Nmat': 2},
        'calc_scheme': 'general',
    }
    solver = solver_flex.FLEX(ham, {}, info_mode)
    solver.spin_mode = "spin-free"
    solver.nmat = 1
    return solver


class TestChi0ConventionMatchesMYO(unittest.TestCase):
    """Lock the RPA-vs-MYO orbital-pair transpose as an intentional fact.

    ``RPA._calc_chi0q`` returns the bare bubble in the RPA convention
    ``chi0_opt[..., a, c, b, d]``.  The MYO brute-force reference returns
    ``chi0_bf[m, n, mu, nu, q, iv]``.  The relation found (and locked by the
    fix) is ``chi0_opt[q, a, c, b, d] == chi0_bf[b, d, a, c, q, 0]`` (MYO
    ``(m,n,mu,nu) = (b,d,a,c)``), i.e. the RPA row pair ``(a,c)`` is MYO's
    column ``(mu,nu)`` and the RPA column pair ``(b,d)`` is MYO's row ``(m,n)``.
    """

    def test_calc_chi0q_matches_bruteforce_transposed(self):
        from tests.flex_bruteforce_ref import chi0_bruteforce
        solver = _make_chi0_general_flex()
        nk = solver.lattice.nvol
        self.assertEqual(nk, 4)
        no = solver.norb
        self.assertEqual(no, 2)

        rng = np.random.default_rng(2024)
        shape = (1, 1, nk, no, no)   # (nblock, nmat, k, a, b)
        green_kw = (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)).astype(np.complex128)

        beta = 1.0 / solver.T
        chi0_opt = solver._calc_chi0q(green_kw, np.zeros_like(green_kw), beta)
        # strip the (size-1) spin-block dimension -> (1, nk, no, no, no, no)
        chi0_opt = chi0_opt[0]
        self.assertEqual(chi0_opt.shape, (1, nk, no, no, no, no))

        # Brute-force reference. G_bf[a,b,k,0] = green_kw[0,0,k,a,b].
        G_bf = np.zeros((no, no, nk, 1), dtype=np.complex128)
        for k in range(nk):
            G_bf[:, :, k, 0] = green_kw[0, 0, k]
        chi0_bf = chi0_bruteforce(G_bf, T=solver.T, Nk=nk)

        # chi0_opt[0, q, a, c, b, d] == chi0_bf[b, d, a, c, q, 0]
        for q in range(nk):
            for a in range(no):
                for c in range(no):
                    for b in range(no):
                        for d in range(no):
                            np.testing.assert_allclose(
                                chi0_opt[0, q, a, c, b, d],
                                chi0_bf[b, d, a, c, q, 0],
                                atol=1e-10)


class TestGeneralPipelinePhysical(unittest.TestCase):
    """The FULL optimized general FLEX pipeline must match the physical
    brute-force pipeline at nmat=1 (this is the end-to-end lock for the fix)."""

    def test_pipeline_matches_physical(self):
        from tests.flex_bruteforce_ref import chi0_bruteforce, sigma_bruteforce
        solver = _make_chi0_general_flex()
        nk = solver.lattice.nvol
        no = solver.norb
        self.assertEqual((nk, no), (4, 2))
        beta = 1.0 / solver.T

        rng = np.random.default_rng(7777)
        shape = (1, 1, nk, no, no)
        green_kw = (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)).astype(np.complex128)

        # --- OPTIMIZED general pipeline ---
        chi0_raw = solver._calc_chi0q(
            green_kw, np.zeros_like(green_kw), beta)[0]   # (1,nk,no,no,no,no)
        chi0, Us, Uc = solver._inflate_chi0q_and_ham_general(chi0_raw, None)
        chi_s, chi_c = solver._solve_channels_general(chi0, Us, Uc)
        v_eff = solver._calc_veff_general(chi0, chi_s, chi_c, Us, Uc)
        sigma_opt = solver._calc_self_energy_general(green_kw, v_eff, beta)
        self.assertEqual(sigma_opt.shape, (1, 1, nk, no, no))

        # --- PHYSICAL (paper convention) brute-force pipeline ---
        G_bf = np.zeros((no, no, nk, 1), dtype=np.complex128)
        for k in range(nk):
            G_bf[:, :, k, 0] = green_kw[0, 0, k]
        # chi0_phys axes (m, n, mu, nu, q, iv)
        chi0_phys = chi0_bruteforce(G_bf, solver.T, nk)

        ndx = no * no
        V_phys = np.zeros((no, no, no, no, nk, 1), dtype=np.complex128)
        eye = np.eye(ndx, dtype=np.complex128)
        # The physical reference reuses the production Us/Uc; this is intentional.
        # This test locks the RPA->MYO chi0 transpose + the channel/veff/contraction
        # wiring, NOT the S/C matrix *values* (those are independently pinned by
        # TestMYOSCMatrices.test_myo_elements). The transpose fix is what makes
        # sigma_opt match sigma_phys here regardless of the (shared) S/C values.
        for q in range(nk):
            M0 = chi0_phys[:, :, :, :, q, 0].reshape(ndx, ndx)  # row=(m,n),col=(mu,nu)
            Us_q = Us[q]
            Uc_q = Uc[q]
            chi_s_q = np.linalg.inv(eye - M0 @ Us_q) @ M0
            chi_c_q = np.linalg.inv(eye + M0 @ Uc_q) @ M0
            V_q = (1.5 * Us_q @ chi_s_q @ Us_q
                   + 0.5 * Uc_q @ chi_c_q @ Uc_q
                   - 0.25 * (Us_q + Uc_q) @ M0 @ (Us_q + Uc_q))
            # V_q is the production V_eff at this q: its MYO flatten has
            # row=(mu,m), col=(nu,n) (the same flatten that
            # _sigma_orbital_contract reshapes back to axes (mu,m,nu,n)). That
            # is exactly the V_phys[mu,m,nu,n] layout the brute-force sigma
            # expects, so reshape directly with NO axis permutation.
            V_phys[:, :, :, :, q, 0] = V_q.reshape(no, no, no, no)  # (mu,m,nu,n)

        sigma_phys = sigma_bruteforce(G_bf, V_phys, solver.T, nk)  # (m,n,k,iw)

        max_diff = 0.0
        for k in range(nk):
            for m in range(no):
                for n in range(no):
                    d = abs(sigma_opt[0, 0, k, m, n] - sigma_phys[m, n, k, 0])
                    max_diff = max(max_diff, d)
        print("\nTestGeneralPipelinePhysical max abs diff = {:.3e}".format(
            max_diff))
        for k in range(nk):
            for m in range(no):
                for n in range(no):
                    np.testing.assert_allclose(
                        sigma_opt[0, 0, k, m, n],
                        sigma_phys[m, n, k, 0],
                        atol=1e-8)


def _run_flex_sigma(scheme, *, norb1=True, U=None, extra_interactions=None,
                    Nmat=16, Lx=8, IterationMax=100, Mix=0.3, T=2.0, EPS=8):
    """Build a FLEX solver, run ``solve()`` to convergence, return ``sigma``.

    Deterministic small-but-converged FLEX run used by the limit/regression
    tests.  ``norb1=True`` uses the 1-orbital ``tests/rpa/input`` fixture (via
    ``_make_solver``); ``norb1=False`` uses the 2-orbital
    ``tests/rpa/input_2orb`` fixture (CoulombIntra + CoulombInter, J=0).

    The interaction strength U is fixed by the fixture data files
    (CoulombIntra = 4.0); the ``U`` keyword is accepted for API symmetry but the
    fixtures do not parameterize it, so it is intentionally not applied here.
    Returns the converged ``green_info["sigma"]`` of shape
    ``(1, Nmat, Lx*Lx, norb, norb)``.
    """
    param = {
        'T': T, 'mu': 0.0,
        'CellShape': [Lx, Lx, 1], 'SubShape': [1, 1, 1],
        'Nmat': Nmat, 'IterationMax': IterationMax, 'Mix': Mix, 'EPS': EPS,
    }
    if norb1:
        solver, green_info, info_file = _make_solver(
            "FLEX", Lx=Lx, Ly=Lx, Nmat=Nmat, T=T, calc_scheme=scheme,
            extra_mode={'param': param}, extra_interactions=extra_interactions)
        out = info_file['output']['path_to_output']
        solver.solve(green_info, out)
        return green_info["sigma"]

    # 2-orbital: build directly from a self-contained ON-SITE 2-orbital fixture
    # (mirrors the _construct / TestGeneralSolveEndToEnd path, but WITHOUT
    # patching exchange).  The general path rejects off-site two-body terms, so
    # both the general and squashed runs here use on-site U' for a fair compare.
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    dirpath = os.path.join(tempfile.gettempdir(), "hwave_flex_2d_2orb_onsite")
    _write_2d_2orb_onsite_fixture(dirpath)
    info_input = {
        'path_to_input': dirpath,
        'interaction': {
            'path_to_input': dirpath,
            'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
            'CoulombInter': 'coulombinter.dat',
        },
    }
    if extra_interactions:
        info_input['interaction'].update(extra_interactions)
    ham = read_input_k.QLMSkInput(info_input).get_param('ham')
    green_info = read_input_k.QLMSkInput(info_input).get_param('green')
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': scheme}
    solver = solver_flex.FLEX(ham, {}, info_mode)
    os.makedirs('tests/flex/output_2orb', exist_ok=True)
    solver.solve(green_info, 'tests/flex/output_2orb')
    return green_info["sigma"]


class TestGeneralLimits(unittest.TestCase):
    """Task 9: limit & regression tests for the general FLEX path.

    Step 1 (single-orbital general == reduced) is the key physics lock: for one
    orbital the MYO matrices are 1x1 with Us=Uc=U, and the general fluctuation
    kernel V = 1.5 U chi_s U + 0.5 U chi_c U - 0.25 (2U) chi0 (2U)
    = U^2 (1.5 chi_s + 0.5 chi_c - chi0) is exactly the reduced density-density
    kernel.  The converged self-energies must therefore coincide.  (Verified:
    max abs diff == 0.0, ratio == 1 exactly.)
    """

    def test_single_orbital_equals_reduced(self):
        sig_gen = _run_flex_sigma("general", norb1=True)
        sig_red = _run_flex_sigma("reduced", norb1=True)
        # Exact in the paramagnetic single-orbital Hubbard limit.
        np.testing.assert_allclose(sig_gen, sig_red, atol=1e-8, rtol=1e-6)

    def test_multiorbital_general_differs_from_reduced(self):
        # 2-orbital, CoulombIntra + CoulombInter (U, U'); J=0.  The general path
        # keeps the inter-orbital U' vertices that the reduced/squashed
        # density-density reduction drops, so the converged sigmas must differ.
        sig_gen = _run_flex_sigma("general", norb1=False)
        sig_red = _run_flex_sigma("squashed", norb1=False)
        self.assertGreater(np.linalg.norm(sig_gen - sig_red), 1e-8)
        self.assertTrue(np.all(np.isfinite(sig_gen))
                        and np.all(np.isfinite(sig_red)))
        # Bounded / sane (not a runaway divergence).
        self.assertLess(np.linalg.norm(sig_gen - sig_red), 1e3)

    def test_general_no_density_density_warning(self):
        """General-path construction with an exchange-type interaction present
        must NOT emit the density-density reduction warning.

        This is the focused Hund-suppression counterpart to
        TestFLEXGeneralWarningGating.test_warning_suppressed_for_general: it
        reuses the exact same construction (which patches
        ``Interaction.has_interaction_exchange`` to True so an exchange/Hund-type
        interaction is seen as present), and asserts the warning stays silent on
        the new-physics general path.
        """
        gating = TestFLEXGeneralWarningGating()
        with _assert_no_warning(self, 'hwave.solver.flex'):
            gating._construct('general')


if __name__ == '__main__':
    unittest.main()
