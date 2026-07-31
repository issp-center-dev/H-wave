#!/usr/bin/env python3
"""Tests for the RPA transverse (ladder) channel, calc_type="ring+ladder".

What is covered, and why it looks the way it does:

1. The transverse vertex is pinned per interaction type against values
   established by exact diagonalization (with the O(lambda) self-energy
   removed via the Hartree-Fock response), including the orbital-pair
   symmetrisation convention shared with `uhfk.py` and its per-slot-family
   conjugation rule (plain mean on palindromic pairs, whose declarations
   multiply the same operator; Hermitian mean on heterogeneous pairs, whose
   declarations are Hermitian-partner coefficients -- PairHop).
2. SU(2): `chi_zz = chi_pm` is asserted only where the Hamiltonian is
   genuinely SU(2) symmetric, over the FULL orbital-pair matrix. For on-site
   CoulombInter the identity is exact but currently VIOLATED, because the
   longitudinal ring drops the inter-orbital spin vertex -- that test pins
   the specific #104 signature and must not be "fixed" from the transverse
   side.
3. Off-site input: interactions whose assembled transverse vertex would be
   q-dependent are rejected before the longitudinal solve; off-site Hund and
   PairLift are accepted because their transverse vertex vanishes.
4. Ring-only behaviour and the bare-bubble (Lindhard) checks are unchanged.
"""

import os
import unittest

import numpy as np


class TestRPALadder(unittest.TestCase):
    """Test ladder diagram implementation via SU(2) symmetry."""

    def _run_rpa(self, calc_type="ring", calc_scheme="auto",
                 spin_orbital=False, inter=None, Lx=8, Ly=8, Nmat=64,
                 interactions=None, input_path='tests/rpa/input',
                 T=2.0, filling=0.75):
        """Run RPA and return (solver, green_info).

        Parameters
        ----------
        interactions : dict, optional
            Full interaction dict to use instead of default CoulombIntra.
            If None, uses {'CoulombIntra': 'coulombintra.dat'} plus `inter`.
        inter : dict, optional
            Additional interactions to add to the default.
        input_path : str, optional
            Path to input files directory (default: 'tests/rpa/input').
        T : float, optional
            Temperature (default: 2.0).
        filling : float, optional
            Electron filling (default: 0.75).
        """
        if inter is None:
            inter = {}

        info_log = {}
        info_mode = {
            'mode': 'RPA',
            'param': {
                'T': T,
                'filling': filling,
                'CellShape': [Lx, Ly, 1],
                'SubShape': [1, 1, 1],
                'Nmat': Nmat,
            },
            'enable_spin_orbital': spin_orbital,
            'calc_scheme': calc_scheme,
            'calc_type': calc_type,
        }

        if interactions is None:
            interaction_dict = {
                'path_to_input': input_path,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            }
            interaction_dict.update(inter)
        else:
            interaction_dict = {
                'path_to_input': input_path,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
            }
            interaction_dict.update(interactions)

        info_file = {
            'input': {
                'path_to_input': input_path,
                'interaction': interaction_dict,
            },
            'output': {
                'path_to_output': 'tests/rpa/output',
                'chi0q': 'chi0q.npz',
            },
        }
        os.makedirs(info_file['output']['path_to_output'], exist_ok=True)

        import hwave.qlmsio.read_input_k as read_input_k
        read_io = read_input_k.QLMSkInput(info_file['input'])
        ham_info = read_io.get_param("ham")

        import hwave.solver.rpa as solver_rpa
        solver = solver_rpa.RPA(ham_info, info_log, info_mode)

        green_info = read_io.get_param("green")
        solver.solve(green_info, info_file['output']['path_to_output'])

        return solver, green_info

    def test_transverse_spin_diag_external_chi0q_raises(self):
        """spin-diag transverse channel needs G_up, G_down to build
        chi0_+- = G_up*G_down. With an externally-supplied chi0q (no Green's
        functions), chi0_+- cannot be reconstructed from the longitudinal
        chi0q, so it must raise instead of silently using the up-up block."""
        import types
        import hwave.solver.rpa as rpa_module

        stub = object.__new__(rpa_module.RPA)
        stub.norb, stub.ns, stub.nd = 1, 2, 2
        stub.lattice = types.SimpleNamespace(nvol=4)
        stub.spin_mode = "spin-diag"
        stub.T, stub.nmat = 1.0, 4
        stub.freq_index = np.arange(4)
        stub.enable_reduced = True
        # no green0 attribute -> the external-chi0q fallback path

        nfreq, nvol, norb, nd = 4, 4, 1, 2
        chi0q_orig = np.zeros((2, nfreq, nvol, norb, norb), dtype=complex)
        ham_orig = np.zeros((nvol, nd, nd, nd, nd), dtype=complex)
        with self.assertRaises(ValueError):
            stub._build_transverse_channel(chi0q_orig, ham_orig)

    def test_transverse_spin_free_rejects_reduced_chi0q(self):
        """The spin-free branch once expanded a reduced chi0q with the
        defective delta_{l2,l4} placement; the expansion was dead code
        (ring+ladder forces calc_scheme='general') and is now a loud
        guard. Pin the guard directly, since no normally constructed
        solver can reach it (issue #111).

        All the malformed-shape guards in this method raise ValueError
        consistently: an earlier revision used AssertionError here with
        a claimed invariant-vs-runtime distinction, which a review
        showed was false -- both branches are equally unreachable
        through normal construction (the upstream gate is pinned by
        test_upstream_gate_pins_the_dead_branch_premise below)."""
        import types
        import hwave.solver.rpa as rpa_module

        stub = object.__new__(rpa_module.RPA)
        stub.norb, stub.ns = 2, 2
        stub.lattice = types.SimpleNamespace(nvol=4)
        stub.spin_mode = "spin-free"

        # a reduced (2-index) chi0q: ndim 4, not the rank-4 orbital tensor
        chi0q_reduced = np.zeros((4, 4, 2, 2), dtype=complex)
        ham_orig = np.zeros((4, 4, 4, 4, 4), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._build_transverse_channel(chi0q_reduced, ham_orig)
        self.assertIn("density-pair", str(cm.exception).lower())
        self.assertIn("general", str(cm.exception))

    def test_transverse_spinful_rejects_reduced_chi0q(self):
        """Same guard on the spinful branch: a non-6-dim chi0q cannot
        supply the spin-flip block and must be rejected (issue #111)."""
        import types
        import hwave.solver.rpa as rpa_module

        stub = object.__new__(rpa_module.RPA)
        # ns is ALWAYS 2 in a constructed RPA; a reduced spinful bubble
        # for norb=2 ends in (nd, nd) = (4, 4), so this is the state that
        # could exist just before the upstream scheme validation.
        stub.norb, stub.ns = 2, 2
        stub.lattice = types.SimpleNamespace(nvol=4)
        stub.spin_mode = "spinful"

        chi0q_reduced = np.zeros((4, 4, 4, 4), dtype=complex)
        ham_orig = np.zeros((4, 4, 4, 4, 4), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            stub._build_transverse_channel(chi0q_reduced, ham_orig)
        # pin the INTENDED failure, not just any ValueError: the message
        # must identify the channel, the requirement, and the shape
        msg = str(cm.exception)
        self.assertIn("spinful transverse", msg)
        self.assertIn("requires the full", msg)
        self.assertIn("general-scheme", msg)
        self.assertIn(str(chi0q_reduced.shape), msg)

    def test_upstream_gate_pins_the_dead_branch_premise(self):
        """The two guard tests above exercise states normal construction
        cannot produce. Pin the PREMISE that makes them unreachable --
        the invariant under which the old expansion branches were
        removed: ring+ladder forces the general scheme at construction,
        and every normally reached transverse call receives the full
        rank-4 (6-dim) bubble."""
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        input_path = 'tests/rpa/input'
        info_file_input = {
            'path_to_input': input_path,
            'interaction': {
                'path_to_input': input_path,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        }

        def _construct(calc_scheme):
            info_mode = {
                'mode': 'RPA',
                'param': {'T': 2.0, 'filling': 0.75,
                          'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                          'Nmat': 16},
                'enable_spin_orbital': False,
                'calc_scheme': calc_scheme,
                'calc_type': 'ring+ladder',
            }
            read_io = read_input_k.QLMSkInput(info_file_input)
            return rpa_mod.RPA(read_io.get_param("ham"), {}, info_mode)

        # reduced/squashed exit before any solve
        for scheme in ("reduced", "squashed"):
            with self.subTest(scheme=scheme):
                with self.assertRaises(SystemExit):
                    _construct(scheme)
        # auto resolves to general; general constructs normally
        self.assertEqual(_construct("auto").calc_scheme, "general")
        self.assertEqual(_construct("general").calc_scheme, "general")

        # every normally reached transverse call gets a 6-dim bubble
        seen = []
        original = rpa_mod.RPA._build_transverse_channel

        def spy(self_, chi0q_orig, ham_orig):
            seen.append(chi0q_orig.ndim)
            return original(self_, chi0q_orig, ham_orig)

        rpa_mod.RPA._build_transverse_channel = spy
        try:
            self._run_rpa(calc_type="ring+ladder", Lx=4, Ly=4, Nmat=16)
        finally:
            rpa_mod.RPA._build_transverse_channel = original
        self.assertTrue(seen, "the ladder solve must reach the channel")
        self.assertEqual(seen, [6] * len(seen))

    def test_su2_symmetry_1orb_coulombintra(self):
        """For 1-orbital CoulombIntra, ring+ladder should give chi_zz = chi_+-.

        chi_zz = chi_uu - chi_ud  (from longitudinal RPA)
        chi_+- = transverse RPA susceptibility (from ladder channel)

        SU(2) symmetry requires chi_zz = chi_+- for paramagnetic systems.

        Note: in the density correlation convention used by rpa.py,
        chi_zz = chi_uu - chi_ud (not 2*(chi_uu - chi_ud)),
        and chi_+- = chi_zz (not chi_zz/2).
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2

        # chi_zz from longitudinal channel
        chi_uu = chiq[iw0, :, 0, 0, 0, 0]
        chi_ud = chiq[iw0, :, 0, 0, 1, 1]
        chi_zz = chi_uu - chi_ud

        # chi_+- from transverse channel
        chi_pm = chiq_pm[iw0, :, 0, 0, 0, 0]

        np.testing.assert_allclose(
            chi_zz, chi_pm,
            atol=1e-10,
            err_msg="SU(2) symmetry violated: chi_zz != chi_+- for ring+ladder"
        )

    def test_ring_only_no_transverse(self):
        """For ring-only mode, chiq_pm should not exist."""
        solver, green_info = self._run_rpa(
            calc_type="ring",
            calc_scheme="general",
        )
        self.assertNotIn("chiq_pm", green_info,
                          "Ring-only mode should not compute transverse chi")

    def test_ring_backward_compatibility(self):
        """Ring-only with general scheme should produce non-trivial results."""
        solver, green_info = self._run_rpa(
            calc_type="ring",
            calc_scheme="general",
        )
        chiq = green_info["chiq"]
        nfreq = chiq.shape[0]
        iw0 = nfreq // 2
        chi_uu = chiq[iw0, :, 0, 0, 0, 0]
        self.assertTrue(np.any(np.abs(chi_uu) > 1e-10),
                        "chi_uu should be non-zero")

    def test_transverse_vertex_coulombintra(self):
        """Verify the transverse vertex W_+- = -U for CoulombIntra.

        For H = U n_up n_dn, the Fock exchange gives W_+-(diag) = -U.
        The transverse RPA is: chi_+- = chi0 / (1 + U*chi0)
        which gives the SAME result as chi_s = chi0 / (1 - U*chi0)
        when the sign convention accounts for W_+- = -U.

        SU(2) requires chi_s = chi_+- in density convention.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=32,
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2

        chi_s = chiq[iw0, :, 0, 0, 0, 0] - chiq[iw0, :, 0, 0, 1, 1]
        chi_pm = chiq_pm[iw0, :, 0, 0, 0, 0]

        np.testing.assert_allclose(
            chi_s, chi_pm,
            atol=1e-10,
            err_msg="Transverse vertex: chi_s != chi_+- for CoulombIntra"
        )

    def test_su2_all_matsubara_frequencies(self):
        """SU(2) symmetry should hold at ALL Matsubara frequencies, not just iw0."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder",
            calc_scheme="general",
            Lx=4, Ly=4, Nmat=16,
        )

        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]

        for iw in range(nfreq):
            chi_zz = chiq[iw, :, 0, 0, 0, 0] - chiq[iw, :, 0, 0, 1, 1]
            chi_pm = chiq_pm[iw, :, 0, 0, 0, 0]
            np.testing.assert_allclose(
                chi_zz, chi_pm,
                atol=1e-10,
                err_msg=f"SU(2) violated at Matsubara frequency index {iw}"
            )


    def _assert_su2(self, green_info, norb=1, msg=""):
        """Assert SU(2): chi_zz = chi_+- at every q, over the FULL
        norb^2 x norb^2 orbital-pair matrix.

        This used to sum only the pair-DIAGONAL elements `chiq[..., a,a,a,a]`
        at norb > 1. Those are invariant under the orbital-pair transpose and
        blind to anything living at an off-diagonal pair index, which is
        exactly where the inter-orbital vertex sits -- so the reduced form
        reported agreement while chi_zz and chi_pm differed by 9 percent.
        """
        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]

        nfreq = chiq.shape[0]
        iw0 = nfreq // 2
        npair = norb * norb

        chi_zz = np.zeros((chiq.shape[1], npair, npair), dtype=complex)
        chi_pm = np.zeros_like(chi_zz)
        for a in range(norb):
            for c in range(norb):
                for b in range(norb):
                    for d in range(norb):
                        i, j = a * norb + c, b * norb + d
                        chi_zz[:, i, j] = (
                            chiq[iw0, :, a, c, b, d]
                            - chiq[iw0, :, a, c, norb + b, norb + d])
                        chi_pm[:, i, j] = chiq_pm[iw0, :, a, c, b, d]

        np.testing.assert_allclose(
            chi_zz, chi_pm, atol=1e-10,
            err_msg="SU(2) symmetry violated: {}".format(msg))

    def test_transverse_vertex_matches_exact_diagonalization(self):
        """Pin `ham_pm` for every interaction type against exact
        diagonalization.

        Each value below was obtained by diagonalizing an explicit 3-site,
        2-orbital, spinful model and removing the O(lambda) self-energy exactly
        -- it is the linear response to the Hartree-Fock one-body potential
        built from the non-interacting density matrix, so subtracting it leaves
        the vertex alone. Every case is q-independent to better than 1e-6,
        which is the diagnostic that the extraction is clean: an on-site
        interaction cannot give a q-dependent vertex. The ED-to-code scale is
        one constant for all seven types (3.048 to 3.055, i.e. N_site), which is
        what makes the magnitudes and signs -- not merely the structure -- part
        of the check.

        The table that used to sit in `_build_transverse_channel` had four of
        these wrong and one missing:

            CoulombIntra  -U             table said -U    correct
            CoulombInter  -U'            table said 0     wrong
            Hund           0             table said -J    wrong
            Ising         +J             table said 2J    wrong
            Exchange      -(J+J^T)/2     table said 0     wrong
            PairLift       0             table said 0     correct
            PairHop       -J             no entry

        The orbital pair is symmetrised with the MEAN, which is what an
        interaction file already means here: `uhfk.py` builds every two-body
        table as `(jab_r + jba)/2`, and the two solvers hold identical values
        for symmetric input.
        """
        import hwave.solver.rpa as rpa_mod

        U0, U1, V = 4.0, 4.0, 0.7      # coulombintra.dat, onsite_inter.dat
        expected = {
            # (a, c, b, d) -> value
            "CoulombIntra": {(0, 0, 0, 0): -U0, (1, 1, 1, 1): -U1},
            "CoulombInter": {(0, 1, 0, 1): -V, (1, 0, 1, 0): -V},
            "Hund": {},
            "Exchange": {(0, 0, 1, 1): -V, (1, 1, 0, 0): -V},
            "Ising": {(0, 1, 0, 1): +V, (1, 0, 1, 0): +V},
            "PairLift": {},
            "PairHop": {(0, 1, 1, 0): -V, (1, 0, 0, 1): -V},
        }
        files = {"CoulombIntra": "coulombintra.dat"}

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        for itype, want in expected.items():
            with self.subTest(interaction=itype):
                rpa_mod.RPA._build_transverse_channel = spy
                try:
                    self._run_rpa(
                        calc_type="ring+ladder", calc_scheme="general",
                        input_path='tests/rpa/input_2orb',
                        Lx=4, Ly=4, Nmat=32, T=2.0, filling=0.5,
                        interactions={itype: files.get(itype,
                                                      "onsite_inter.dat")})
                finally:
                    rpa_mod.RPA._build_transverse_channel = original

                w = captured["ham_pm"]
                self.assertTrue(
                    np.allclose(w, w[0:1], atol=1e-12),
                    "an on-site interaction must give a q-independent vertex")
                got = w[0]
                built = np.zeros_like(got)
                for idx, val in want.items():
                    built[idx] = val
                np.testing.assert_allclose(
                    got, built, atol=1e-10,
                    err_msg="{} vertex does not match exact "
                            "diagonalization".format(itype))

    def test_complex_pairhop_keeps_its_full_complex_value(self):
        """A complex Hermitian-closed PairHop must not lose its imaginary part.

        Hermiticity of H requires P_10 = conj(P_01); the exact vertex (#100)
        carries the FULL complex value. The symmetrisation must therefore
        average with the CONJUGATE of the pair swap -- as `uhfk.py` does
        (`jba = conjugate(transpose(...))`). A plain transpose average
        collapses both slots to Re(P), which was measured on an earlier
        revision of this branch: ham_pm = -0.700 + 0.000i for P = 0.7 + 0.4i.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        rpa_mod.RPA._build_transverse_channel = spy
        try:
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"PairHop": "onsite_cplx_hermitian.dat"})
        finally:
            rpa_mod.RPA._build_transverse_channel = original

        got = captured["ham_pm"][0]
        want = np.zeros_like(got)
        want[0, 1, 1, 0] = -(0.7 + 0.4j)
        want[1, 0, 0, 1] = -(0.7 - 0.4j)
        np.testing.assert_allclose(got, want, atol=1e-10)
        # the load-bearing part: the imaginary component survives
        self.assertGreater(float(np.max(np.abs(got.imag))), 0.39)

    def test_complex_density_density_coupling_is_real(self):
        """Complex Hermitian-closed CoulombInter and Ising must yield REAL
        vertices, exactly like Exchange.

        n_a n_b = n_b n_a, so the two declarations multiply the same operator,
        the physical coefficient is the (real) sum, and the declared imaginary
        parts multiply the null operator -- in the many-body Hamiltonian they
        cannot even enter, since a density-density term is diagonal in the
        Fock basis and a complex diagonal would be non-Hermitian.

        This is the case that showed the per-BLOCK conjugation rule was the
        wrong abstraction: density-density types live on the same cross-spin
        block as PairHop but in palindromic slots, and the conjugated mean
        there kept a spurious -0.7 -/+ 0.4i. The rule is per slot family.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        for itype, val in (("CoulombInter", -0.7), ("Ising", +0.7)):
            with self.subTest(interaction=itype):
                rpa_mod.RPA._build_transverse_channel = spy
                try:
                    self._run_rpa(
                        calc_type="ring+ladder", calc_scheme="general",
                        input_path='tests/rpa/input_2orb',
                        Lx=4, Ly=4, Nmat=32, T=2.0, filling=0.5,
                        interactions={itype: "onsite_cplx_hermitian.dat"})
                finally:
                    rpa_mod.RPA._build_transverse_channel = original

                got = captured["ham_pm"][0]
                want = np.zeros_like(got)
                want[0, 1, 0, 1] = want[1, 0, 1, 0] = val
                np.testing.assert_allclose(got, want, atol=1e-10)
                self.assertLess(
                    float(np.max(np.abs(got.imag))), 1e-12,
                    "the density-density coupling (V_01 + V_10)/2 is real; an "
                    "imaginary remnant means the conjugated mean was applied "
                    "to a palindromic slot family")

    def test_complex_exchange_coupling_is_real(self):
        """A complex Hermitian-closed Exchange must yield a REAL vertex.

        The two declarations multiply the SAME operator: X_ab^dagger = X_ba
        and the different-spin bilinears commute, so X_ab = X_ba and
        H = (J_01 + J_10) X. Hermiticity requires that sum to be real, so for
        J = (j, conj(j)) the physical coupling is Re(j) -- here 0.7 -- and any
        surviving imaginary part is an artifact. A conjugated mean on the
        spin-flip block produced exactly that artifact (-0.7 +/- 0.4i); the
        plain mean is correct on this block, in contrast to the cross-spin
        block where PairHop requires the conjugated one.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        rpa_mod.RPA._build_transverse_channel = spy
        try:
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"Exchange": "onsite_cplx_hermitian.dat"})
        finally:
            rpa_mod.RPA._build_transverse_channel = original

        got = captured["ham_pm"][0]
        want = np.zeros_like(got)
        want[0, 0, 1, 1] = want[1, 1, 0, 0] = -0.7
        np.testing.assert_allclose(got, want, atol=1e-10)
        self.assertLess(float(np.max(np.abs(got.imag))), 1e-12,
                        "the Exchange coupling (J_01 + J_10)/2 is real; an "
                        "imaginary remnant means the wrong Hermitian partner "
                        "was used on the spin-flip block")

    def test_a_redundant_declaration_of_zero_gives_a_zero_vertex(self):
        """`V_01 = +1`, `V_10 = -1` denotes an identically zero Hamiltonian,
        because `n_a n_b = n_b n_a`. The vertex must vanish.

        Read unsymmetrised, this produced +1.00 and -1.00 at the two slots --
        i.e. the response depended on how the coefficient was split between two
        entries that denote the same operator. Interaction files are not checked
        for orbital symmetry (#93), so such input is reachable.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        for itype in ("CoulombInter", "Ising", "Exchange"):
            # Exchange too: X_ab^dagger = X_ba are the same Hermitian operator
            # pair, so an antisymmetric J also denotes a zero Hamiltonian.
            with self.subTest(interaction=itype):
                captured.clear()
                rpa_mod.RPA._build_transverse_channel = spy
                try:
                    self._run_rpa(
                        calc_type="ring+ladder", calc_scheme="general",
                        input_path='tests/rpa/input_2orb',
                        Lx=4, Ly=4, Nmat=32, T=2.0, filling=0.5,
                        interactions={itype: "onsite_inter_antisym.dat"})
                finally:
                    rpa_mod.RPA._build_transverse_channel = original
                self.assertLess(
                    float(np.max(np.abs(captured["ham_pm"]))), 1e-12,
                    "a declaration denoting a zero Hamiltonian must give a "
                    "zero vertex")

    def test_exchange_vertex_uses_only_the_symmetric_part(self):
        """Exchange depends only on J_ab + J_ba, and the vertex must too.

        `X_ab^dagger = X_ba`, so `H = sum_ab J_ab (X_ab + X_ba)` and the
        antisymmetric part of J cancels identically. Measured on an asymmetric
        real J = (1.0, 0.35) -- physical, since Exchange is Hermitian for any
        real J -- exact diagonalization requires -(1.0 + 0.35) at BOTH slots,
        while the raw interaction block holds -0.35 and -1.00 separately.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        rpa_mod.RPA._build_transverse_channel = spy
        try:
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"Exchange": "onsite_inter_asym.dat"})
        finally:
            rpa_mod.RPA._build_transverse_channel = original

        got = captured["ham_pm"][0]
        want = np.zeros_like(got)
        want[0, 0, 1, 1] = want[1, 1, 0, 0] = -(1.0 + 0.35) / 2
        np.testing.assert_allclose(got, want, atol=1e-10)

    def test_offsite_two_body_is_rejected(self):
        """The one-orbital `coulombinter.dat` carries OFF-SITE bonds, and the
        transverse channel cannot represent those.

        The transverse pair is ``c+_(i a up) c_(j b down)``; when the
        cross-spin part of the interaction is off-site its two ends sit on
        different sites, the pair is non-local, and `ham_pm` -- which carries
        only on-site orbital pairs -- has nowhere to put it. Measured: the
        extracted vertex is not proportional to V(q) at all (residual 1.0),
        while the longitudinal one is proportional to it to 1e-8.

        Seven SU(2) cases used to run here and passed, but only because the old
        vertex construction happened to cancel to zero for them. They were never
        valid tests of the transverse channel.
        """
        for name, inter in (
            ("CoulombInter", {'CoulombIntra': 'coulombintra.dat',
                              'CoulombInter': 'coulombinter.dat'}),
            ("CoulombInter only", {'CoulombInter': 'coulombinter.dat'}),
            ("Ising", {'CoulombIntra': 'coulombintra.dat',
                       'Ising': 'coulombinter.dat'}),
        ):
            with self.subTest(interaction=name):
                with self.assertRaises(ValueError):
                    self._run_rpa(calc_type="ring+ladder",
                                  calc_scheme="general",
                                  Lx=4, Ly=4, Nmat=32, interactions=inter)

    def test_offsite_pair_swap_respects_the_displacement(self):
        """Swapping the two bilinears of an off-site term also reverses the
        displacement, so the same-operator partner of ``V_ab(+r)`` is
        ``V_ba(-r)``, not ``V_ba(+r)``.

        Two sharp cases pin this. ``V_ab(+x) = +1, V_ba(+x) = -1`` is a
        genuinely NONZERO Hamiltonian -- at fixed q the swap made it cancel to
        an identically zero vertex, slip through the guard, and lose its
        physics silently. ``V_ab(+x) = +1, V_ba(-x) = -1`` IS identically zero
        -- and was rejected. Both are now handled correctly.
        """
        import hwave.solver.rpa as rpa_mod

        # the nonzero Hamiltonian must be rejected (off-site cross-spin)
        with self.assertRaises(ValueError):
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"CoulombInter": "offsite_asym_nonzero.dat"})

        # the identically zero Hamiltonian must be accepted, with a zero vertex
        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        captured.clear()
        rpa_mod.RPA._build_transverse_channel = spy
        try:
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"CoulombInter": "offsite_redundant_zero.dat"})
        finally:
            rpa_mod.RPA._build_transverse_channel = original
        self.assertIn("ham_pm", captured, "the spy must actually have run")
        self.assertLess(float(np.max(np.abs(captured["ham_pm"]))), 1e-12)

    def test_spin_diag_transverse_bubble_matches_the_longitudinal_kernel(self):
        """With G_up == G_dn the transverse bubble must equal the longitudinal
        general chi0 -- both are -G[a,b](r,t) G[d,c](-r,-t), and the
        longitudinal kernel is the one verified against exact diagonalization
        in this series.

        The displacement reversal in `_calc_chi0q_transverse` used to be
        applied AFTER flattening the lattice to nvol, and over the orbital
        axes as well. Both errors are invisible for nd <= 2 (roll(-1)+flip is
        the identity on an axis of length <= 2) and on effectively-1D
        lattices, which is exactly what every existing fixture used. At
        nd = 3 on a 3 x 2 lattice the disagreement was of order unity (2.8).
        """
        import types
        import hwave.solver.rpa as rpa_mod

        rng = np.random.RandomState(7)
        nx, ny, nz, nd, nmat = 3, 2, 1, 3, 8
        nvol = nx * ny * nz
        stub = object.__new__(rpa_mod.RPA)
        stub.lattice = types.SimpleNamespace(shape=(nx, ny, nz), nvol=nvol)
        stub.enable_reduced = False
        stub.fft_workers = 1
        stub.nmat = nmat
        g = (rng.randn(nmat, nvol, nd, nd)
             + 1j * rng.randn(nmat, nvol, nd, nd))
        gkw = np.stack([g, g])                    # G_up == G_dn
        tail = np.zeros_like(gkw)
        pm = stub._calc_chi0q_transverse(gkw, tail, beta=3.0)
        lg = stub._calc_chi0q(gkw, tail, beta=3.0)
        np.testing.assert_allclose(pm, lg[0], rtol=0.0, atol=1e-12)

    def test_a_reused_green_info_does_not_leak_stale_results(self):
        """chiq / chiq_pm must belong to the CURRENT solve.

        save_results writes whatever `chiq_pm` key is present, so a ring-only
        solve on a green_info that carries one from an earlier ring+ladder run
        would silently label the stale transverse data as part of the new
        result; and a solve that fails validation must not leave an earlier
        chiq behind either. The stale keys are injected directly -- feeding
        the previous run's full green_info back in also hands the solver its
        stored chi0q, which exercises the (separately broken) external-chi0q
        path and would test the wrong thing.
        """
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as rpa_mod

        def make(calc_type, interactions):
            interaction_dict = {'path_to_input': 'tests/rpa/input_2orb',
                                'Geometry': 'geom.dat',
                                'Transfer': 'transfer.dat'}
            interaction_dict.update(interactions)
            read_io = read_input_k.QLMSkInput(
                {'path_to_input': 'tests/rpa/input_2orb',
                 'interaction': interaction_dict})
            solver = rpa_mod.RPA(
                read_io.get_param("ham"), {},
                {'mode': 'RPA',
                 'param': {'T': 2.0, 'filling': 0.5,
                           'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
                           'Nmat': 32},
                 'enable_spin_orbital': False,
                 'calc_scheme': "general", 'calc_type': calc_type})
            return solver, read_io

        os.makedirs('tests/rpa/output', exist_ok=True)
        stale = np.zeros((1,), dtype=complex)

        # ring-only: an injected stale chiq_pm must be cleared
        solver, read_io = make("ring", {'CoulombInter': 'onsite_inter.dat'})
        green_info = read_io.get_param("green")
        green_info["chiq_pm"] = stale
        solver.solve(green_info, 'tests/rpa/output')
        self.assertNotIn("chiq_pm", green_info)
        self.assertIn("chiq", green_info)

        # failing validation: neither stale key may survive
        solver2, read_io2 = make("ring+ladder",
                                 {'CoulombInter': 'coulombinter.dat'})
        green_info2 = read_io2.get_param("green")
        green_info2["chiq"] = stale
        green_info2["chiq_pm"] = stale
        with self.assertRaises(ValueError):
            solver2.solve(green_info2, 'tests/rpa/output')
        self.assertNotIn("chiq", green_info2)
        self.assertNotIn("chiq_pm", green_info2)

    def test_tiny_offsite_coupling_is_still_rejected(self):
        """The guard tolerance is relative with NO absolute floor.

        A floor of ``max(scale, 1)`` accepted a purely q-dependent off-site
        vertex whenever its amplitude sat below 1e-10 -- 100 percent q
        variation, silently used by the unsupported calculation. Fourier
        roundoff is ~1e-16 relative, so a relative 1e-10 separates it from
        genuine q-dependence at any amplitude.
        """
        with self.assertRaises(ValueError):
            self._run_rpa(
                calc_type="ring+ladder", calc_scheme="general",
                input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                T=2.0, filling=0.5,
                interactions={"CoulombInter": "offsite_tiny.dat"})

    def test_sublattice_folding_keeps_the_verdicts(self):
        """SubShape != [1,1,1]: the vertex, the guard and the q-reversal all
        act on the FOLDED lattice, so they must keep working when the shape,
        the orbital count and the qrev map all change.

        The pathological off-site declaration stays rejected -- folding turns
        part of the bond intra-cell, but an inter-cell remainder survives at
        any SubShape < CellShape -- and the identically-zero declaration stays
        accepted, because the operator it denotes is zero on any folding.
        """
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        def run(interactions, cell, sub):
            info_mode = {
                'mode': 'RPA',
                'param': {'T': 2.0, 'filling': 0.5, 'CellShape': cell,
                          'SubShape': sub, 'Nmat': 32},
                'enable_spin_orbital': False,
                'calc_scheme': "general", 'calc_type': "ring+ladder",
            }
            interaction_dict = {'path_to_input': 'tests/rpa/input_2orb',
                                'Geometry': 'geom.dat',
                                'Transfer': 'transfer.dat'}
            interaction_dict.update(interactions)
            info_file = {'input': {'path_to_input': 'tests/rpa/input_2orb',
                                   'interaction': interaction_dict},
                         'output': {'path_to_output': 'tests/rpa/output',
                                    'chi0q': 'chi0q.npz'}}
            os.makedirs('tests/rpa/output', exist_ok=True)
            import hwave.qlmsio.read_input_k as read_input_k
            read_io = read_input_k.QLMSkInput(info_file['input'])
            solver = rpa_mod.RPA(read_io.get_param("ham"), {}, info_mode)
            solver.solve(read_io.get_param("green"),
                         info_file['output']['path_to_output'])

        for cell in ([4, 4, 1], [4, 6, 1]):
            with self.subTest(cell=cell):
                # on-site: accepted, q-independent, correct magnitude
                captured.clear()
                rpa_mod.RPA._build_transverse_channel = spy
                try:
                    run({'CoulombInter': 'onsite_inter.dat'}, cell, [2, 2, 1])
                finally:
                    rpa_mod.RPA._build_transverse_channel = original
                w = captured["ham_pm"]
                self.assertLess(float(np.max(np.abs(w - w[0:1]))), 1e-12)
                self.assertAlmostEqual(float(np.max(np.abs(w))), 0.7,
                                       places=10)

                # nonzero off-site: rejected
                with self.assertRaises(ValueError):
                    run({'CoulombInter': 'offsite_asym_nonzero.dat'},
                        cell, [2, 2, 1])

                # identically zero: accepted with a zero vertex
                captured.clear()
                rpa_mod.RPA._build_transverse_channel = spy
                try:
                    run({'CoulombInter': 'offsite_redundant_zero.dat'},
                        cell, [2, 2, 1])
                finally:
                    rpa_mod.RPA._build_transverse_channel = original
                self.assertIn("ham_pm", captured,
                              "the spy must actually have run")
                self.assertLess(
                    float(np.max(np.abs(captured["ham_pm"]))), 1e-12)

    def test_su2_onsite_coulombinter_2orb_equality(self):
        """SU(2) holds exactly for on-site CoulombInter, so chi_zz must
        equal chi_pm. This failed while #104 was open (the longitudinal ring
        dropped the inter-orbital spin vertex, measured at 1e-2 in chiq); the
        adjudicated Fierz cross-slot content restored the identity -- the
        residual dropped from 3e-3 to 4e-34 with the transverse side
        untouched.
        """
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder", calc_scheme="general",
            input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
            T=2.0, filling=0.5,
            interactions={'CoulombInter': 'onsite_inter.dat'},
        )
        self._assert_su2(green_info, norb=2,
                         msg="on-site CoulombInter (2-orbital)")

    def test_combined_types_superpose(self):
        """The vertex of a combined declaration is the sum of the individual
        vertices -- block extraction and symmetrisation are linear, so mixed
        interaction files cannot produce cross-talk between types."""
        import hwave.solver.rpa as rpa_mod

        captured = {}
        original = rpa_mod.RPA._build_transverse_channel

        def spy(inner_self, chi0q_orig, ham_orig):
            out = original(inner_self, chi0q_orig, ham_orig)
            captured["ham_pm"] = np.asarray(out[1])
            return out

        def vertex(interactions):
            rpa_mod.RPA._build_transverse_channel = spy
            try:
                self._run_rpa(
                    calc_type="ring+ladder", calc_scheme="general",
                    input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
                    T=2.0, filling=0.5, interactions=interactions)
            finally:
                rpa_mod.RPA._build_transverse_channel = original
            return captured["ham_pm"].copy()

        parts = [{"CoulombIntra": "coulombintra.dat"},
                 {"CoulombInter": "onsite_inter.dat"},
                 {"Exchange": "onsite_inter_asym.dat"},
                 {"PairHop": "onsite_cplx_hermitian.dat"}]
        combined = {}
        for d in parts:
            combined.update(d)
        total = vertex(combined)
        summed = sum(vertex(d) for d in parts)
        np.testing.assert_allclose(total, summed, atol=1e-12)

    def test_offsite_is_accepted_when_it_cannot_reach_the_vertex(self):
        """Off-site Hund and PairLift are fine, and rejecting them would be
        over-strict.

        The guard is on the cross-spin and spin-flip blocks, not on the whole
        tensor, because those are the only blocks the vertex draws on. Hund is
        same-spin only, so it has no cross-spin block; PairLift's transverse
        vertex vanishes identically. Measured off-site: 1.4e-5 and exactly 0
        respectively, against 0.63 for CoulombInter.
        """
        for name, inter in (
            ("Hund", {'CoulombIntra': 'coulombintra.dat',
                      'Hund': 'coulombinter.dat'}),
            ("PairLift", {'CoulombIntra': 'coulombintra.dat',
                          'PairLift': 'coulombinter.dat'}),
        ):
            with self.subTest(interaction=name):
                solver, green_info = self._run_rpa(
                    calc_type="ring+ladder", calc_scheme="general",
                    Lx=4, Ly=4, Nmat=32, interactions=inter)
                self.assertIn("chiq_pm", green_info)

    def test_su2_onsite_coulombinter_2orb_no_residual_anywhere(self):
        """Retired #104 signature test, inverted: while the issue was open,
        the chi_zz / chi_pm mismatch sat exactly on the U' pair slots (0,1)
        and (1,0) at 3e-3 and nowhere else. With the adjudicated cross-slot
        content in the longitudinal W, the mismatch must be gone EVERYWHERE
        -- including those slots -- not merely reduced."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder", calc_scheme="general",
            input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
            T=2.0, filling=0.5,
            interactions={'CoulombInter': 'onsite_inter.dat'},
        )
        norb = 2
        npair = norb * norb
        chiq = green_info["chiq"]
        chiq_pm = green_info["chiq_pm"]
        iw0 = chiq.shape[0] // 2
        zz = np.zeros((chiq.shape[1], npair, npair), dtype=complex)
        pm = np.zeros_like(zz)
        for a in range(norb):
            for c in range(norb):
                for b in range(norb):
                    for d in range(norb):
                        i, j = a * norb + c, b * norb + d
                        zz[:, i, j] = (chiq[iw0, :, a, c, b, d]
                                       - chiq[iw0, :, a, c, norb + b, norb + d])
                        pm[:, i, j] = chiq_pm[iw0, :, a, c, b, d]
        self.assertLess(float(np.max(np.abs(zz - pm))), 1e-10,
                        "chi_zz and chi_pm must agree in every pair slot")

    def test_su2_onsite_coulombintra_2orb(self):
        """CoulombIntra at two orbitals: SU(2) holds and the ring is correct
        for it, so this one passes -- the control for the case above."""
        solver, green_info = self._run_rpa(
            calc_type="ring+ladder", calc_scheme="general",
            input_path='tests/rpa/input_2orb', Lx=4, Ly=4, Nmat=32,
            T=2.0, filling=0.5,
            interactions={'CoulombIntra': 'coulombintra.dat'},
        )
        self._assert_su2(green_info, norb=2, msg="CoulombIntra (2-orbital)")

    def _compute_lindhard(self, ek, mu, T, qx, qy, Lx, Ly):
        """Compute the Lindhard function chi0_orb(q) analytically.

        Returns the static orbital susceptibility (per spin).
        """
        nvol = Lx * Ly
        beta = 1.0 / T
        KX_flat = np.array([i // Ly for i in range(nvol)])
        KY_flat = np.array([i % Ly for i in range(nvol)])
        xk = ek - mu
        fk = 1.0 / (np.exp(beta * xk) + 1.0)

        chi0 = 0.0
        for ik in range(nvol):
            kpqx = (KX_flat[ik] + qx) % Lx
            kpqy = (KY_flat[ik] + qy) % Ly
            ikpq = kpqx * Ly + kpqy
            dx = xk[ikpq] - xk[ik]
            if abs(dx) > 1e-12:
                chi0 += (fk[ikpq] - fk[ik]) / dx
            else:
                chi0 += (-beta * fk[ik] * (1.0 - fk[ik]))
        chi0 *= -1.0 / nvol
        return chi0

    def _compute_fd_response_real_space(self, ek_real_space_builder, mu, T,
                                         qx, qy, Lx, Ly, dh=1e-6,
                                         field_type="z"):
        """Compute static response via finite-difference in real space.

        Builds the full 2N x 2N real-space Hamiltonian with applied field,
        diagonalizes, and computes the order parameter change.

        Returns chi_zz(q) or chi_xx(q).
        """
        N = Lx * Ly
        beta = 1.0 / T
        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])

        # Phase factor for the applied field: e^{i q.r} in real space
        phase_r = np.zeros(N, dtype=float)
        for ir in range(N):
            phase_r[ir] = np.cos(2 * np.pi * (qx * KX_flat[ir] / Lx
                                                + qy * KY_flat[ir] / Ly))

        def build_H(h):
            H = np.zeros((2 * N, 2 * N))
            # Hopping (from the real-space builder)
            H[:N, :N], H[N:, N:] = ek_real_space_builder()
            # Field
            for i in range(N):
                if field_type == "z":
                    H[i, i] += h * phase_r[i]
                    H[N + i, N + i] -= h * phase_r[i]
                elif field_type == "x":
                    H[i, N + i] += h * phase_r[i]
                    H[N + i, i] += h * phase_r[i]
            return H

        def compute_order_param(H):
            w, v = np.linalg.eigh(H)
            f_occ = 1.0 / (np.exp(beta * (w - mu)) + 1.0)
            rho = v @ np.diag(f_occ) @ v.conj().T
            m = 0.0
            for ir in range(N):
                if field_type == "z":
                    m += phase_r[ir] * (rho[ir, ir] - rho[N + ir, N + ir]).real
                elif field_type == "x":
                    m += phase_r[ir] * (rho[ir, N + ir]
                                        + rho[N + ir, ir]).real
            return m / N

        Hp = build_H(dh)
        Hm = build_H(-dh)
        mp = compute_order_param(Hp)
        mm = compute_order_param(Hm)
        return (mp - mm) / (2 * dh)

    def test_chi0_uniform_matches_lindhard(self):
        """Verify chi0_zz(q=0) from Lindhard matches finite-difference.

        The Lindhard chi0 > 0 (physics convention), while the response
        dM/dh < 0 (field opposes magnetization). So chi0_zz = -dM/dh.
        """
        Lx, Ly = 8, 8
        T = 2.0
        N = Lx * Ly

        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])
        ek = -2.0 * (np.cos(2 * np.pi * KX_flat / Lx)
                      + np.cos(2 * np.pi * KY_flat / Ly))

        from scipy.optimize import brentq
        beta = 1.0 / T

        def filling(mu):
            return 2 * np.mean(1.0 / (np.exp(beta * (ek - mu)) + 1.0)) - 0.75
        mu = brentq(filling, -10, 10)

        # Lindhard: chi0 > 0
        chi0_orb = self._compute_lindhard(ek, mu, T, 0, 0, Lx, Ly)
        chi0_zz = 2 * chi0_orb

        # Analytic response: dM/dh = -chi0_zz
        xk = ek - mu
        fk = 1.0 / (np.exp(beta * xk) + 1.0)
        response = -2 * beta * np.mean(fk * (1 - fk))  # negative

        np.testing.assert_allclose(
            chi0_zz, -response,
            rtol=1e-8,
            err_msg="chi0_zz(q=0): Lindhard vs analytic mismatch"
        )

    def test_chi0_staggered_matches_fd(self):
        """Verify chi0_zz(q=pi,pi) from Lindhard matches real-space FD."""
        Lx, Ly = 6, 6
        T = 2.0
        N = Lx * Ly

        KX_flat = np.array([i // Ly for i in range(N)])
        KY_flat = np.array([i % Ly for i in range(N)])
        ek = -2.0 * (np.cos(2 * np.pi * KX_flat / Lx)
                      + np.cos(2 * np.pi * KY_flat / Ly))

        from scipy.optimize import brentq
        beta = 1.0 / T

        def filling(mu):
            return 2 * np.mean(1.0 / (np.exp(beta * (ek - mu)) + 1.0)) - 0.75
        mu = brentq(filling, -10, 10)

        # Lindhard
        Qx, Qy = Lx // 2, Ly // 2
        chi0_orb = self._compute_lindhard(ek, mu, T, Qx, Qy, Lx, Ly)
        chi0_zz = 2 * chi0_orb

        # Real-space finite-difference
        def build_hopping():
            H = np.zeros((N, N))
            for ir in range(N):
                irx = KX_flat[ir]
                iry = KY_flat[ir]
                for dx, dy in [(1, 0), (-1, 0), (0, 1), (0, -1)]:
                    jx = (irx + dx) % Lx
                    jy = (iry + dy) % Ly
                    jr = jx * Ly + jy
                    H[ir, jr] -= 1.0
            return H, H.copy()

        chi_fd = self._compute_fd_response_real_space(
            build_hopping, mu, T, Qx, Qy, Lx, Ly, dh=1e-6, field_type="z")

        # Lindhard chi0 > 0, FD response dM/dh < 0: chi0_zz = -chi_fd
        np.testing.assert_allclose(
            chi0_zz, -chi_fd,
            rtol=1e-4,
            err_msg="chi0_zz(q=pi,pi): Lindhard vs real-space FD mismatch"
        )


class TestRPAParamValidation(unittest.TestCase):
    """RPA finite-temperature Matsubara setup requires T > 0 and an even Nmat
    (odd Nmat would place omega=0 in the fermionic grid). The chemical-potential
    search also needs 0 < Ncond < Nstate."""

    def _build(self, T=2.0, Nmat=64, Ncond=2):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.rpa as solver_rpa
        info_mode = {
            'mode': 'RPA',
            'param': {'T': T, 'Ncond': Ncond, 'CellShape': [8, 8, 1],
                      'SubShape': [1, 1, 1], 'Nmat': Nmat},
            'calc_scheme': 'general',
        }
        info_input = {
            'path_to_input': 'tests/rpa/input',
            'interaction': {'path_to_input': 'tests/rpa/input',
                            'Geometry': 'geom.dat', 'Transfer': 'transfer.dat',
                            'CoulombIntra': 'coulombintra.dat'},
        }
        read_io = read_input_k.QLMSkInput(info_input)
        return solver_rpa.RPA(read_io.get_param("ham"), {}, info_mode)

    def test_T_zero_rejected(self):
        with self.assertRaises(SystemExit):
            self._build(T=0.0)

    def test_odd_nmat_rejected(self):
        with self.assertRaises(SystemExit):
            self._build(Nmat=5)

    def test_ncond_full_filling_rejected(self):
        # Nstate = nvol * nd = 64 * 2 = 128; Ncond=128 is full filling
        with self.assertRaises(SystemExit):
            self._build(Ncond=128)

    def test_ncond_zero_rejected(self):
        with self.assertRaises(SystemExit):
            self._build(Ncond=0)

    def test_nmat_one_rejected(self):
        with self.assertRaises(SystemExit):
            self._build(Nmat=1)

    def test_min_even_nmat_ok(self):
        self._build(Nmat=2)  # smallest valid (even) grid -> no exit

    def test_valid_params_ok(self):
        self._build()  # T>0, even Nmat, 0<Ncond<Nstate -> no exit


if __name__ == '__main__':
    unittest.main()
