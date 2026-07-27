#!/usr/bin/env python3

"""Reduced (kuroki) FLEX susceptibilities consumed by the Eliashberg solver.

A ``calc_scheme="reduced"`` FLEX run stores ``chiq_s``/``chiq_c`` in the reduced
spin-orbital space, where the matrix index is a DENSITY PAIR: the stored
``X[i, j]`` is ``chi_{(i,i),(j,j)}`` with ``i = s*norb + a`` (see
``FLEX._inflate_chi0q_and_ham``, whose companion interaction reduction is
``einsum('ksasatbtb->ksatb', ...)`` -- the density-density diagonal of the
vertex).

Embedding that matrix into the ``norb^2 x norb^2`` orbital-pair space used by
the pairing vertex therefore has exactly one correct placement:

    out[(a,a), (b,b)] = X[a, b],   everything else zero

i.e. flat indices ``out[a*norb + a, b*norb + b]``.  The historical placement
``out[l1*norb + l2, l3*norb + l2] = X[l1, l3]`` (a delta_{l2,l4} "spectator"
scatter, equivalent to ``kron(X, I_norb)``) reads the density-pair index as an
ordinary orbital index.  It drops the genuine inter-orbital density coupling
``chi_{(0,0),(1,1)}`` and injects X into off-density positions, which changes
the singlet vertex and the static lambda.

Tests must be run from the repository root (they use relative paths like
``tests/rpa/input``).
"""

import os
import tempfile
import unittest

import numpy as np


def _write_2orb_intra_only_fixture(dirpath):
    """A 2-orbital fixture with CoulombIntra ONLY.

    With no inter-orbital two-body term the Kuroki S/C matrices are nonzero
    only on the density-pair block (``_build_sc_matrices_all_q`` case 1), so the
    reduced density-density treatment is EXACT and the reduced pairing vertex
    must agree element-for-element with the 2-index "simple" vertex.  That makes
    the equivalence assertion below a genuine physical check rather than an
    approximation comparison.
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
    with open(os.path.join(dirpath, "geom.dat"), "w") as f:
        f.write(geom)
    with open(os.path.join(dirpath, "transfer.dat"), "w") as f:
        f.write(transfer)
    with open(os.path.join(dirpath, "coulombintra.dat"), "w") as f:
        f.write(coulombintra)


def _make_reduced_flex(dirpath):
    """Build a spin-free ``calc_scheme="reduced"`` FLEX solver on the fixture.

    ``dirpath`` must be a directory owned by the caller (a per-test temporary
    one), so parallel test processes never write and read the same fixture
    files and nothing is left behind afterwards.
    """
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex

    _write_2orb_intra_only_fixture(dirpath)
    info_input = {
        'path_to_input': dirpath,
        'interaction': {
            'path_to_input': dirpath,
            'Geometry': 'geom.dat',
            'Transfer': 'transfer.dat',
            'CoulombIntra': 'coulombintra.dat',
        },
    }
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    info_mode = {
        'mode': 'FLEX',
        'param': {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                  'SubShape': [1, 1, 1], 'Nmat': 8},
        'calc_scheme': 'reduced',
    }
    solver = solver_flex.FLEX(ham, {}, info_mode)
    solver.spin_mode = "spin-free"
    return solver


class TestExpandFlexChiDensityPlacement(unittest.TestCase):
    """``_expand_flex_chi`` must embed the reduced matrix as a density pair."""

    def test_kuroki_block_lands_on_density_pair_positions(self):
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 3
        nvol, nd_so, nd = Nx * Ny * Nz, norb * 2, norb * norb
        rng = np.random.default_rng(11)
        chi_raw = (rng.standard_normal((nfreq, nvol, nd_so, nd_so))
                   + 1j * rng.standard_normal((nfreq, nvol, nd_so, nd_so)))

        out = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz,
                                  convention="kuroki")

        chi6 = chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_so, nd_so)
        X = chi6[:, :, :, :, :norb, :norb]      # spin-up density-pair block
        expected = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                expected[:, :, :, :, a * norb + a, b * norb + b] = X[..., a, b]

        np.testing.assert_allclose(out, expected, atol=1e-14)

    def test_kuroki_keeps_interorbital_density_coupling(self):
        """The specific component the delta_{l2,l4} scatter silently dropped:
        chi_{(0,0),(1,1)} = X[0,1] (flat position [0, nd-1])."""
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nfreq = 2, 1, 1, 1, 1
        nd_so, nd = norb * 2, norb * norb
        chi_raw = np.zeros((nfreq, 1, nd_so, nd_so), dtype=complex)
        chi_raw[0, 0, 0, 1] = 7.0 + 3.0j        # X[0,1]

        out = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz,
                                  convention="kuroki")

        self.assertAlmostEqual(out[0, 0, 0, 0, 0, nd - 1], 7.0 + 3.0j)
        # ...and nowhere else.
        rest = out.copy()
        rest[0, 0, 0, 0, 0, nd - 1] = 0.0
        self.assertAlmostEqual(np.max(np.abs(rest)), 0.0)

    def test_off_density_rows_and_columns_are_zero(self):
        """A reduced chi carries no information about pair indices (a,b) with
        a != b, so those rows/columns must stay exactly zero rather than being
        filled with copies of the density block."""
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nfreq = 3, 2, 1, 1, 2
        nvol, nd_so, nd = Nx * Ny * Nz, norb * 2, norb * norb
        rng = np.random.default_rng(23)
        chi_raw = (rng.standard_normal((nfreq, nvol, nd_so, nd_so))
                   + 1j * rng.standard_normal((nfreq, nvol, nd_so, nd_so)))

        out = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz,
                                  convention="kuroki")

        density = [a * norb + a for a in range(norb)]
        off = [i for i in range(nd) if i not in density]
        self.assertAlmostEqual(np.max(np.abs(out[..., off, :])), 0.0)
        self.assertAlmostEqual(np.max(np.abs(out[..., :, off])), 0.0)


class TestSingleOrbitalIsBitIdentical(unittest.TestCase):
    """The release note claims norb=1 results are BIT-IDENTICAL across this
    change. Substantiate it against the old placement itself, with exact
    equality rather than a tolerance."""

    @staticmethod
    def _old_scatter(chi_full, norb, nfreq, Nx, Ny, Nz):
        """The pre-fix embedding, reproduced verbatim as the reference."""
        nd = norb * norb
        chi_orb = chi_full[:, :, :, :, :norb, :norb]
        out = np.zeros((nfreq, Nx, Ny, Nz, nd, nd), dtype=complex)
        for l2 in range(norb):
            out[:, :, :, :, l2::norb, l2::norb] = chi_orb
        return out

    def test_expand_flex_chi_norb1_matches_old_placement_exactly(self):
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nfreq = 1, 3, 2, 1, 4
        nvol, nd_so = Nx * Ny * Nz, norb * 2
        rng = np.random.default_rng(77)
        chi_raw = (rng.standard_normal((nfreq, nvol, nd_so, nd_so))
                   + 1j * rng.standard_normal((nfreq, nvol, nd_so, nd_so)))

        new = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz,
                                  convention="kuroki")
        old = self._old_scatter(
            chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_so, nd_so),
            norb, nfreq, Nx, Ny, Nz)

        self.assertTrue(np.array_equal(new, old),
                        "norb=1 must be bit-identical to the old placement")

    def test_multi_orbital_genuinely_differs_from_old_placement(self):
        """Guard against the bit-identity test passing vacuously: for norb >= 2
        the two placements must NOT coincide."""
        import hwave.sc as sc

        for norb in (2, 3):
            with self.subTest(norb=norb):
                Nx, Ny, Nz, nfreq = 2, 1, 1, 2
                nvol, nd_so = Nx * Ny * Nz, norb * 2
                rng = np.random.default_rng(80 + norb)
                chi_raw = (rng.standard_normal((nfreq, nvol, nd_so, nd_so))
                           + 1j * rng.standard_normal(
                               (nfreq, nvol, nd_so, nd_so)))
                new = sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz,
                                          convention="kuroki")
                old = self._old_scatter(
                    chi_raw.reshape(nfreq, Nx, Ny, Nz, nd_so, nd_so),
                    norb, nfreq, Nx, Ny, Nz)
                self.assertFalse(np.allclose(new, old))


class TestReducedFlexVertexMatchesLoadVertex(unittest.TestCase):
    """Physical regression: for Sigma=0 the reduced-FLEX pairing vertex must
    equal the vertex the ``load`` path builds from the same chi0q.

    ``chi0q_mode="flex"`` on a reduced run feeds ``chi_s``/``chi_c`` through
    ``_expand_flex_chi`` into ``_compute_vertices_flex`` (norb^2 orbital-pair
    algebra), while ``chi0q_mode="load"`` on the same reduced run feeds chi0q
    into ``_compute_vertices_simple`` (norb x norb algebra).  With CoulombIntra
    only, the S/C matrices are confined to the density-pair block and the two
    formulations are algebraically identical, so the vertices must agree
    exactly:  Vs_flex[a, a, b, b] == (Pc + Ps)[a, b], with Vs_flex zero on every
    off-density component.
    """

    def setUp(self):
        # Per-test fixture dir: parallel test processes must not share it.
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)

    def _build(self):
        import hwave.sc as sc

        dirpath = self._tmp.name
        flex = _make_reduced_flex(dirpath)
        norb = flex.norb
        Nx, Ny, Nz = flex.lattice.shape
        nvol = flex.lattice.nvol
        nd = norb * norb
        nmat = flex.nmat
        si = nmat // 2

        # One random spin-free reduced chi0q feeds BOTH paths.
        rng = np.random.default_rng(7)
        shape = (nmat, nvol, norb, norb)
        chi0_raw = (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)) * 0.05

        # FLEX path: solve the reduced RPA channels exactly as a Sigma=0 run
        # would, giving chi_s/chi_c in the stored (nmat, nvol, nd_so, nd_so)
        # spin-orbital layout.
        _, _, chi_s, chi_c = flex._flex_compute_veff(
            chi0_raw, flex.ham_info.ham_inter_q)
        chi_s = np.asarray(chi_s)
        chi_c = np.asarray(chi_c)

        input_dict = {"file": {"input": {"interaction": {
            "path_to_input": dirpath,
            "Geometry": "geom.dat", "Transfer": "transfer.dat",
            "CoulombIntra": "coulombintra.dat"}}}}
        _, _, interactions = sc._read_interaction_files(input_dict)
        kx = np.linspace(0, 2 * np.pi, Nx, endpoint=False)
        ky = np.linspace(0, 2 * np.pi, Ny, endpoint=False)
        kz = np.linspace(0, 2 * np.pi, Nz, endpoint=False)
        inter_k = sc._build_interaction_k(kx, ky, kz, interactions, norb)

        chis_stat = sc._expand_flex_chi(chi_s[si:si + 1], norb, Nx, Ny, Nz,
                                        "kuroki")[0]
        chic_stat = sc._expand_flex_chi(chi_c[si:si + 1], norb, Nx, Ny, Nz,
                                        "kuroki")[0]

        # load path consumes chi0q as (norb, norb, Nx, Ny, Nz, nmat)
        chi0q_load = chi0_raw.reshape(nmat, Nx, Ny, Nz, norb, norb).transpose(
            4, 5, 1, 2, 3, 0)

        return (sc, inter_k, chis_stat, chic_stat, chi0q_load,
                norb, Nx, Ny, Nz, nmat, si, nd)

    def test_singlet_and_triplet_vertices_match_simple_path(self):
        (sc, inter_k, chis_stat, chic_stat, chi0q_load,
         norb, Nx, Ny, Nz, nmat, si, nd) = self._build()

        for pairing in ("singlet", "triplet"):
            Vs_flex = sc._compute_vertices_flex(
                chis_stat, chic_stat, inter_k, norb, Nx, Ny, Nz,
                pairing_type=pairing, convention="kuroki")
            Pc_q, Ps_q = sc._compute_vertices_simple(
                chi0q_load, inter_k, norb, Nx, Ny, Nz, nmat,
                pairing_type=pairing, static_index=si)
            Vs_load = Pc_q + Ps_q          # (norb, norb, Nx, Ny, Nz)

            for a in range(norb):
                for b in range(norb):
                    np.testing.assert_allclose(
                        Vs_flex[a, a, b, b], Vs_load[a, b],
                        rtol=1e-9, atol=1e-11,
                        err_msg="reduced flex vertex ({}) must match the "
                                "simple/load vertex at density pair "
                                "({},{})".format(pairing, a, b))

    def test_vertex_has_no_off_density_components(self):
        """With CoulombIntra only, S and C live entirely on the density-pair
        block, so the reduced vertex must vanish elsewhere.  A nonzero
        off-density component means chi leaked into pair positions the reduced
        scheme never computed."""
        (sc, inter_k, chis_stat, chic_stat, chi0q_load,
         norb, Nx, Ny, Nz, nmat, si, nd) = self._build()

        Vs_flex = sc._compute_vertices_flex(
            chis_stat, chic_stat, inter_k, norb, Nx, Ny, Nz,
            pairing_type="singlet", convention="kuroki")

        mask = np.ones((norb, norb, norb, norb), dtype=bool)
        for a in range(norb):
            for b in range(norb):
                mask[a, a, b, b] = False
        self.assertAlmostEqual(np.max(np.abs(Vs_flex[mask])), 0.0)


class TestGeneralVertexTwoIndexEmbedding(unittest.TestCase):
    """The same density-pair embedding must hold on the ``load`` route.

    ``_compute_vertices_general`` also accepts a 2-index (reduced) chi0q, which
    it lifts into the norb^2 orbital-pair space itself.  That branch is reached
    by the ordinary two-step workflow: RPA with ``calc_scheme="auto"`` resolves
    to ``"reduced"`` for U + Hund/Ising/PairHop and writes a 2-index chi0q, and
    ``hwave_sc`` then routes to the general S/C formulation because an
    inter-orbital vertex term is present.

    The reduced chi0q carries the same density-pair index as the reduced FLEX
    chi, so it needs the same placement.  Pinning it with a chi0q that is
    EXACTLY density-diagonal makes the 2-index reduction lossless, so the
    4-index and 2-index vertices must agree exactly -- any difference is the
    embedding's fault alone.
    """

    def _run(self, keys):
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nmat = 2, 3, 2, 1, 4
        si = nmat // 2
        rng = np.random.default_rng(3)

        X = (rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))
             + 1j * rng.standard_normal((norb, norb, Nx, Ny, Nz, nmat))) * 0.1
        X4 = np.zeros((norb, norb, norb, norb, Nx, Ny, Nz, nmat), dtype=complex)
        for a in range(norb):
            for b in range(norb):
                X4[a, a, b, b] = X[a, b]

        inter_k = {}
        for k in keys:
            m = rng.standard_normal((norb, norb, Nx, Ny, Nz))
            inter_k[k] = (m + m.transpose(1, 0, 2, 3, 4)).astype(complex)

        out = {}
        for pairing in ("singlet", "triplet"):
            out[pairing] = (
                sc._compute_vertices_general(X4, inter_k, norb, Nx, Ny, Nz,
                                             nmat, pairing_type=pairing,
                                             static_index=si),
                sc._compute_vertices_general(X, inter_k, norb, Nx, Ny, Nz,
                                             nmat, pairing_type=pairing,
                                             static_index=si))
        return out

    def test_density_diagonal_chi0q_gives_identical_vertices(self):
        for keys in (["CoulombIntra"],
                     ["CoulombIntra", "Hund"],
                     ["CoulombIntra", "CoulombInter", "Hund", "Ising"],
                     ["CoulombIntra", "Exchange", "PairHop"]):
            res = self._run(keys)
            for pairing, (V4, V2) in res.items():
                with self.subTest(interactions="+".join(keys), pairing=pairing):
                    np.testing.assert_allclose(
                        V2, V4, rtol=1e-9, atol=1e-11,
                        err_msg="a density-diagonal chi0q loses nothing in the "
                                "2-index reduction, so the vertices must match")


class TestReducedFlexMissingComponentWarning(unittest.TestCase):
    """A reduced FLEX chi cannot dress the off-density S/C blocks that
    inter-orbital interactions create, so consuming the two together must say
    so instead of silently returning an approximate lambda."""

    def _inter_k(self, keys, norb=2, Nx=2, Ny=2, Nz=1):
        return {k: np.ones((norb, norb, Nx, Ny, Nz), dtype=complex)
                for k in keys}

    def test_warns_for_each_unsupported_interorbital_term(self):
        import hwave.sc as sc

        for key in ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop"):
            with self.subTest(interaction=key):
                with self.assertLogs("hwave_sc", level="WARNING") as cm:
                    sc._warn_reduced_flex_missing_components(
                        self._inter_k(["CoulombIntra", key]), 2, 2, 2, 1)
                joined = "\n".join(cm.output)
                self.assertIn(key, joined)
                self.assertIn("calc_scheme='general'", joined)

    def test_silent_for_coulombintra_only(self):
        import hwave.sc as sc
        import logging

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components(
                self._inter_k(["CoulombIntra"]), 2, 2, 2, 1)
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_silent_when_interorbital_term_is_identically_zero(self):
        """A configured term whose inter-orbital entries all vanish -- the V = 0
        or J = 0 endpoint of a coupling/pressure scan -- contributes nothing to
        the off-density S/C blocks, so warning about it would be a false
        positive."""
        import hwave.sc as sc
        import logging

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        inter_k = {"CoulombIntra": np.ones((norb, norb, Nx, Ny, Nz),
                                           dtype=complex)}
        for k in ("CoulombInter", "Hund", "Ising", "Exchange", "PairHop"):
            inter_k[k] = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz)
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_no_warning_when_the_terms_cancel_in_the_combined_block(self):
        """The decision is made on the ASSEMBLED S/C matrices, not term by term.

        Case 4 of _build_sc_matrices_all_q builds S = C = Exchange + PairHop on
        the off-density block, so Exchange = -PairHop cancels it exactly. A
        per-term test would announce missing dressing that does not exist."""
        import hwave.sc as sc
        import logging

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        jp = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        jp[0, 1] = 0.6
        jp[1, 0] = 0.6
        inter_k = {"CoulombIntra": np.ones((norb, norb, Nx, Ny, Nz),
                                           dtype=complex),
                   "Exchange": jp, "PairHop": -jp}
        self.assertEqual(
            sc._off_density_sc_weight(inter_k, norb, Nx, Ny, Nz), 0.0,
            "Exchange = -PairHop must cancel the case-4 off-density block")

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz)
        finally:
            lg.removeHandler(handler)
        self.assertEqual([r.getMessage() for r in records
                          if r.levelno >= logging.WARNING], [])

    def test_warns_when_the_terms_do_not_cancel(self):
        """Same fixture without the cancellation must still warn, so the test
        above cannot pass for the wrong reason."""
        import hwave.sc as sc

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        jp = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        jp[0, 1] = 0.6
        jp[1, 0] = 0.6
        inter_k = {"CoulombIntra": np.ones((norb, norb, Nx, Ny, Nz),
                                           dtype=complex),
                   "Exchange": jp, "PairHop": jp}
        self.assertGreater(
            sc._off_density_sc_weight(inter_k, norb, Nx, Ny, Nz), 0.0)
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            sc._warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz)
        self.assertIn("Exchange", "\n".join(cm.output))

    def test_diagonal_only_interorbital_term_does_not_warn(self):
        """Only the a != b entries build the off-density blocks, so a term whose
        weight sits purely on the orbital diagonal must not warn."""
        import hwave.sc as sc
        import logging

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        diag = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        diag[0, 0] = 1.0
        diag[1, 1] = 1.0

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components({"Hund": diag}, norb,
                                                     Nx, Ny, Nz)
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_silent_for_single_orbital(self):
        """norb=1 has no off-density pair index, so nothing can be missing even
        with an inter-orbital term configured."""
        import hwave.sc as sc
        import logging

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components(
                self._inter_k(["CoulombIntra", "Hund"], 1, 2, 2, 1), 1, 2, 2, 1)
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_myo_convention_does_not_warn(self):
        """The general (myo) FLEX path stores the full orbital-pair chi, so the
        reduced-data warning must not fire there."""
        import hwave.sc as sc
        import logging

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._warn_reduced_flex_missing_components(
                self._inter_k(["CoulombIntra", "Hund"]), 2, 2, 2, 1,
                convention="myo")
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_compute_vertices_flex_itself_stays_silent(self):
        """The warning must NOT live inside _compute_vertices_flex.

        The dynamic kernel calls that function once per bosonic Matsubara
        frequency (nmat ~ 1000 in production), and
        eliashberg_dynamic._zero_chi_vertex calls it again with chi = 0 to
        extract the bare 0.5*(S+C) term. Warning from inside would flood the log
        and fire on a call where the message makes no sense, so the check is
        hoisted to the two callers that run once per calculation."""
        import hwave.sc as sc
        import logging

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        nd = norb * norb
        inter_k = self._inter_k(["CoulombIntra", "Hund"], norb, Nx, Ny, Nz)
        chi = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc._compute_vertices_flex(chi, chi, inter_k, norb, Nx, Ny, Nz,
                                      pairing_type="singlet",
                                      convention="kuroki")
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING
                 and "REDUCED" in r.getMessage()]
        self.assertEqual(warns, [])

    def test_dynamic_vertex_builder_warns_exactly_once(self):
        """The dynamic path must emit the reduced-data warning once for the
        whole run, not once per frequency."""
        from hwave.solver import eliashberg_dynamic as ed
        import logging

        norb, Nx, Ny, Nz, nmat = 2, 2, 2, 1, 6
        nd = norb * norb
        inter_k = self._inter_k(["CoulombIntra", "Hund"], norb, Nx, Ny, Nz)
        chi_w = np.zeros((Nx, Ny, Nz, nd, nd, nmat), dtype=complex)

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            ed.compute_vertices_flex_dynamic(
                chi_w, chi_w, inter_k, norb, Nx, Ny, Nz,
                pairing_type="singlet", convention="kuroki")
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING
                 and "REDUCED" in r.getMessage()]
        self.assertEqual(
            len(warns), 1,
            "expected exactly 1 reduced-data warning for nmat={}, got {}".format(
                nmat, len(warns)))
        self.assertIn("Hund", warns[0])


class TestDiscardedSpinContentIsRefused(unittest.TestCase):
    """The embedding keeps only the spin-UP block. Dropping the rest is lossless
    exactly when it is redundant (down == up AND cross == 0), which holds
    bit-for-bit for a paramagnetic run. Real discarded content is REFUSED, not
    warned about: the eigenvalue would not approximate anything.

    The guard tests the DISCARDED DATA, never the run's spin mode -- unequal
    diagonal blocks do not imply polarization, and equal ones do not imply
    redundancy when the cross blocks are nonzero."""

    @staticmethod
    def _capture(fn):
        import logging
        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            fn()
        finally:
            lg.removeHandler(handler)
        return [r.getMessage() for r in records if r.levelno >= logging.WARNING]

    def _chi(self, norb, Nx, Ny, Nz, nfreq, down_factor=1.0, cross=0.0):
        nvol, nd_so = Nx * Ny * Nz, norb * 2
        rng = np.random.default_rng(9)
        block = (rng.standard_normal((nfreq, nvol, norb, norb))
                 + 1j * rng.standard_normal((nfreq, nvol, norb, norb)))
        chi = np.zeros((nfreq, nvol, nd_so, nd_so), dtype=complex)
        chi[..., :norb, :norb] = block
        chi[..., norb:, norb:] = block * down_factor
        if cross:
            chi[..., :norb, norb:] = block * cross
            chi[..., norb:, :norb] = block * cross
        return chi

    def _expand(self, **kw):
        """Drive the guard where it lives: the loader boundary.

        _expand_flex_chi itself stays a pure layout transform -- policy about
        which data the solver accepts belongs where the data enters it, not in
        the reshape, which is also exercised directly with arbitrary arrays.
        """
        import hwave.sc as sc
        norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 2
        chi = self._chi(norb, Nx, Ny, Nz, nfreq, **kw)
        return lambda: sc._check_spin_block_discarded(chi, norb, "kuroki")

    def test_paramagnetic_chi_passes_silently(self):
        self.assertEqual(self._capture(self._expand()), [])

    def test_unequal_down_block_is_refused(self):
        with self.assertRaises(ValueError) as cm:
            self._expand(down_factor=0.5)()
        self.assertIn("down-spin block differs", str(cm.exception))

    def test_nonzero_cross_spin_blocks_are_refused(self):
        """Equal diagonal blocks, nonzero cross blocks: discarded just the same,
        so this must not slip through."""
        with self.assertRaises(ValueError) as cm:
            self._expand(cross=0.3)()
        self.assertIn("cross-spin blocks are nonzero", str(cm.exception))

    def test_message_does_not_assert_a_spin_mode(self):
        """The npz records no spin_mode, so the message must describe the data
        rather than classify the run."""
        with self.assertRaises(ValueError) as cm:
            self._expand(down_factor=0.5)()
        low = str(cm.exception).lower()
        for claim in ("spin-diag", "spinful"):
            self.assertNotIn(claim, low)

    def test_message_says_what_the_user_can_do(self):
        with self.assertRaises(ValueError) as cm:
            self._expand(down_factor=0.5)()
        msg = str(cm.exception)
        self.assertIn("coeff_extern", msg)
        self.assertIn("chi^(+-)", msg)

    def test_below_the_roundoff_threshold_only_warns(self):
        """A difference within double-precision noise must not abort a
        legitimate paramagnetic run on a backend whose solve is not
        bit-symmetric."""
        import hwave.sc as sc
        norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 2
        chi = self._chi(norb, Nx, Ny, Nz, nfreq)
        scale = float(np.max(np.abs(chi[..., :norb, :norb])))
        chi[..., norb, norb] += 0.1 * sc._SPIN_DISCARD_ROUNDOFF_RATIO * scale
        warns = self._capture(
            lambda: sc._check_spin_block_discarded(chi, norb, "kuroki"))
        self.assertEqual(len(warns), 1)
        self.assertIn("round-off", warns[0])

    def test_above_the_roundoff_threshold_is_refused(self):
        """...but only just below it. An order of magnitude above must raise, so
        the allowance cannot be mistaken for a general tolerance."""
        import hwave.sc as sc
        norb, Nx, Ny, Nz, nfreq = 2, 2, 2, 1, 2
        chi = self._chi(norb, Nx, Ny, Nz, nfreq)
        scale = float(np.max(np.abs(chi[..., :norb, :norb])))
        chi[..., norb, norb] += 10.0 * sc._SPIN_DISCARD_ROUNDOFF_RATIO * scale
        with self.assertRaises(ValueError):
            sc._check_spin_block_discarded(chi, norb, "kuroki")

    def test_threshold_is_anchored_to_machine_precision(self):
        """Pin the provenance of the constant: it is a few hundred ulp, not a
        guess about how weak a physical field can be."""
        import hwave.sc as sc
        self.assertLess(sc._SPIN_DISCARD_ROUNDOFF_RATIO, 1e-12)
        self.assertGreater(sc._SPIN_DISCARD_ROUNDOFF_RATIO,
                           np.finfo(float).eps)

    def test_myo_chi_is_never_checked(self):
        import hwave.sc as sc
        norb, Nx, Ny, Nz, nfreq = 3, 2, 1, 1, 2
        nd = norb * norb
        rng = np.random.default_rng(13)
        chi = (rng.standard_normal((nfreq, Nx * Ny * Nz, nd, nd))
               + 1j * rng.standard_normal((nfreq, Nx * Ny * Nz, nd, nd)))
        self.assertEqual(
            self._capture(
                lambda: sc._check_spin_block_discarded(chi, norb, "myo")), [])


class TestStaticEntryPointWarnsOnce(unittest.TestCase):
    """End-to-end guard on the PUBLIC static route.

    The helper-level and dynamic-builder tests above do not prove that
    ``calc_eliashberg`` itself reaches the check: the call sits in the
    ``chi0q_mode="flex"`` branch, and a future refactor could drop it without
    failing any of them. Run the real entry point on a reduced (kuroki) chi with
    an inter-orbital term and assert exactly one warning."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = self._tmp.name
        _write_2orb_intra_only_fixture(self.dir)
        # On-site inter-orbital V: reaches the off-density S/C blocks.
        with open(os.path.join(self.dir, "coulombinter.dat"), "w") as f:
            f.write("CoulombInter in wannier90-like format for uhfk\n"
                    "2\n1\n 1\n"
                    "   0    0    0    1    2   1.000000000000   0.000000000000\n"
                    "   0    0    0    2    1   1.000000000000   0.000000000000\n")

        self.norb, self.Nx, self.Ny, self.Nz, self.nmat = 2, 4, 4, 1, 8
        nvol, nd_so = self.Nx * self.Ny * self.Nz, self.norb * 2
        rng = np.random.default_rng(41)

        def rc(shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)) * 0.01

        def paramagnetic_chi():
            """Spin-redundant, as a real paramagnetic reduced run stores it:
            identical spin blocks and zero cross blocks. A plain random array
            would (correctly) be refused as spin-resolved."""
            block = rc((self.nmat, nvol, self.norb, self.norb))
            chi = np.zeros((self.nmat, nvol, nd_so, nd_so), dtype=complex)
            chi[..., :self.norb, :self.norb] = block
            chi[..., self.norb:, self.norb:] = block
            return chi

        np.savez(os.path.join(self.dir, "chiq_s.npz"),
                 chiq_s=paramagnetic_chi(), chi_convention="kuroki")
        np.savez(os.path.join(self.dir, "chiq_c.npz"),
                 chiq_c=paramagnetic_chi(), chi_convention="kuroki")
        np.savez(os.path.join(self.dir, "green.npz"),
                 green=rc((1, self.nmat, nvol, self.norb, self.norb)))

    def _input_dict(self):
        return {
            "mode": {"mode": "RPA", "calc_scheme": "reduced",
                     "param": {"T": 2.0, "mu": 0.0,
                               "CellShape": [self.Nx, self.Ny, self.Nz],
                               "SubShape": [1, 1, 1], "Nmat": self.nmat,
                               "filling": 0.5}},
            "file": {
                "input": {"path_to_input": self.dir,
                          "path_to_flex_output": self.dir,
                          "interaction": {
                              "path_to_input": self.dir,
                              "Geometry": "geom.dat",
                              "Transfer": "transfer.dat",
                              "CoulombIntra": "coulombintra.dat",
                              "CoulombInter": "coulombinter.dat"}},
                "output": {"path_to_output": os.path.join(self.dir, "out")}},
            "eliashberg": {"chi0q_mode": "flex", "frequency": "static",
                           "pairing_type": "singlet", "init_gap": "cos",
                           "solver_mode": "eigenvalue",
                           "eigenvalue_method": "arnoldi",
                           "num_eigenvalues": 2,
                           "output_eigenvalue": "eig.dat",
                           "output_gap": "gap.dat"},
        }

    def test_calc_eliashberg_emits_exactly_one_reduced_warning(self):
        import hwave.sc as sc
        import logging

        records = []
        handler = logging.Handler()
        handler.emit = records.append
        lg = logging.getLogger("hwave_sc")
        lg.addHandler(handler)
        try:
            sc.calc_eliashberg(self._input_dict())
        finally:
            lg.removeHandler(handler)

        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING
                 and "REDUCED" in r.getMessage()]
        self.assertEqual(
            len(warns), 1,
            "static calc_eliashberg must emit the reduced-data warning exactly "
            "once, got {}".format(len(warns)))
        self.assertIn("CoulombInter", warns[0])


class TestExactRedundancyHoldsEndToEnd(unittest.TestCase):
    """The exact (zero-tolerance) spin comparison rests on a claim about real
    producers: a paramagnetic reduced FLEX run writes bit-identical spin blocks
    and exactly-zero cross blocks.

    The synthetic tests above construct redundant arrays by hand and so cannot
    substantiate that. Run the actual solver -- with the high-frequency tail
    correction active, and on both the uniform and IR frequency axes -- and
    assert exact equality on the produced output before asserting silence.
    """

    def _run_flex(self, dirpath, extra_param=None):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as solver_flex

        _write_2orb_intra_only_fixture(dirpath)
        info_input = {
            'path_to_input': dirpath,
            'interaction': {
                'path_to_input': dirpath,
                'Geometry': 'geom.dat',
                'Transfer': 'transfer.dat',
                'CoulombIntra': 'coulombintra.dat',
            },
        }
        ham = read_input_k.QLMSkInput(info_input).get_param("ham")
        param = {'T': 2.0, 'mu': 0.0, 'CellShape': [4, 4, 1],
                 'SubShape': [1, 1, 1], 'Nmat': 8,
                 'coeff_tail': 1.0,          # tail correction ACTIVE
                 'IterationMax': 2, 'Mix': 0.5, 'EPS': 6}
        param.update(extra_param or {})
        solver = solver_flex.FLEX(ham, {}, {'mode': 'FLEX', 'param': param,
                                            'calc_scheme': 'reduced'})
        green_info = read_input_k.QLMSkInput(info_input).get_param("green")
        solver.solve(green_info, dirpath)
        return green_info

    def _assert_exactly_redundant(self, chi, norb, what):
        up = chi[..., :norb, :norb]
        down = chi[..., norb:, norb:]
        self.assertTrue(
            np.array_equal(up, down),
            "{}: spin blocks must be bit-identical for a paramagnetic run "
            "(max|up-down| = {:.3e})".format(
                what, float(np.max(np.abs(up - down)))))
        for lo, hi, name in ((slice(None, norb), slice(norb, None), "up-down"),
                             (slice(norb, None), slice(None, norb), "down-up")):
            cross = chi[..., lo, hi]
            self.assertTrue(
                np.array_equal(cross, np.zeros_like(cross)),
                "{}: {} cross block must be exactly zero (max {:.3e})".format(
                    what, name, float(np.max(np.abs(cross)))))

    def _check(self, extra_param, tag):
        import hwave.sc as sc
        import logging

        with tempfile.TemporaryDirectory() as d:
            green_info = self._run_flex(d, extra_param)
            norb = 2
            records = []
            handler = logging.Handler()
            handler.emit = records.append
            lg = logging.getLogger("hwave_sc")
            lg.addHandler(handler)
            try:
                for key in ("chiq_s", "chiq_c"):
                    chi = np.asarray(green_info[key])
                    self._assert_exactly_redundant(chi, norb,
                                                   "{} {}".format(tag, key))
                    sc._check_spin_block_discarded(chi, norb, "kuroki", key)
            finally:
                lg.removeHandler(handler)
            warns = [r.getMessage() for r in records
                     if r.levelno >= logging.WARNING
                     and "spin structure" in r.getMessage()]
            self.assertEqual(warns, [], "{}: unexpected warning".format(tag))

    def test_uniform_axis_with_tail_correction(self):
        self._check(None, "uniform")

    def test_ir_axis_densified(self):
        try:
            import sparse_ir  # noqa: F401
        except ImportError:
            self.skipTest("sparse_ir not installed")
        self._check({'matsubara_basis': 'ir', 'ir_wmax': 20.0}, "ir-densified")


class TestSpinOrbitalShapeErrorNamesTheCause(unittest.TestCase):
    """A spin-orbital FLEX run reaches hwave_sc as a bare shape mismatch.

    FLEX writes chi in its reduced spin-orbital space (nd = norb_phys * ns)
    while hwave_sc reads norb from geom.dat, which in spin-orbital mode holds
    the spin-orbital count -- so nd_chi lands on exactly norb. The dimensions
    alone look like a corrupt file, and the user has no way to tell that the
    real answer is "the Eliashberg vertex is paramagnetic and cannot take this
    model at all". The message must say so.
    """

    def test_message_identifies_spin_orbital_mode(self):
        import hwave.sc as sc

        # geom norb = 4 spin-orbitals (2 physical); FLEX stores nd = 2*2 = 4.
        norb, Nx, Ny, Nz, nfreq = 4, 2, 1, 1, 1
        nd_chi = 4
        chi_raw = np.zeros((nfreq, Nx * Ny * Nz, nd_chi, nd_chi), dtype=complex)

        with self.assertRaises(ValueError) as cm:
            sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, "kuroki")
        msg = str(cm.exception)
        self.assertIn("enable_spin_orbital", msg)
        self.assertIn("does not support", msg)

    def test_unrelated_shape_mismatch_keeps_the_plain_message(self):
        """The hint must not be attached to every shape error."""
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nfreq = 3, 2, 1, 1, 1
        chi_raw = np.zeros((nfreq, Nx * Ny * Nz, 5, 5), dtype=complex)
        with self.assertRaises(ValueError) as cm:
            sc._expand_flex_chi(chi_raw, norb, Nx, Ny, Nz, "kuroki")
        msg = str(cm.exception)
        self.assertIn("matches neither", msg)
        self.assertNotIn("enable_spin_orbital", msg)


class TestLoadersRefuseSpinResolvedInput(unittest.TestCase):
    """The guard now lives at the loader boundary, so pin it THERE.

    The helper-level tests above cannot catch the check being dropped from
    _load_flex_susceptibilities / _load_flex_susceptibilities_full, which is
    exactly the placement that matters. Cover the static and dynamic routes, and
    the dressed Green function, which discards spin blocks independently of chi.
    """

    def _write(self, d, norb=2, nmat=4, Nx=2, Ny=2, Nz=1,
               down_factor=1.0, cross=0.0, green_blocks=1,
               green_down_factor=1.0):
        nvol, nd_so = Nx * Ny * Nz, norb * 2
        rng = np.random.default_rng(3)

        def rc(shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)) * 0.01

        def chi():
            block = rc((nmat, nvol, norb, norb))
            a = np.zeros((nmat, nvol, nd_so, nd_so), dtype=complex)
            a[..., :norb, :norb] = block
            a[..., norb:, norb:] = block * down_factor
            if cross:
                a[..., :norb, norb:] = block * cross
                a[..., norb:, :norb] = block * cross
            return a

        np.savez(os.path.join(d, "chiq_s.npz"), chiq_s=chi(),
                 chi_convention="kuroki")
        np.savez(os.path.join(d, "chiq_c.npz"), chiq_c=chi(),
                 chi_convention="kuroki")
        g0 = rc((nmat, nvol, norb, norb))
        g = np.stack([g0] + [g0 * green_down_factor
                             for _ in range(green_blocks - 1)])
        np.savez(os.path.join(d, "green.npz"), green=g)
        return {"mode": {"param": {"Nmat": nmat}},
                "file": {"output": {"path_to_output": d}},
                "eliashberg": {"chi0q_mode": "flex"}}, norb, Nx, Ny, Nz

    def _run(self, loader, **kw):
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, **kw)
            return getattr(sc, loader)(inp, norb, Nx, Ny, Nz)

    def test_static_loader_accepts_paramagnetic(self):
        self._run("_load_flex_susceptibilities")

    def test_dynamic_loader_accepts_paramagnetic(self):
        self._run("_load_flex_susceptibilities_full")

    def test_static_loader_refuses_unequal_spin_blocks(self):
        with self.assertRaises(ValueError) as cm:
            self._run("_load_flex_susceptibilities", down_factor=0.5)
        self.assertIn("down-spin block differs", str(cm.exception))

    def test_dynamic_loader_refuses_unequal_spin_blocks(self):
        with self.assertRaises(ValueError) as cm:
            self._run("_load_flex_susceptibilities_full", down_factor=0.5)
        self.assertIn("down-spin block differs", str(cm.exception))

    def test_static_loader_refuses_nonzero_cross_blocks(self):
        with self.assertRaises(ValueError) as cm:
            self._run("_load_flex_susceptibilities", cross=0.3)
        self.assertIn("cross-spin blocks are nonzero", str(cm.exception))

    def test_dynamic_loader_refuses_nonzero_cross_blocks(self):
        with self.assertRaises(ValueError) as cm:
            self._run("_load_flex_susceptibilities_full", cross=0.3)
        self.assertIn("cross-spin blocks are nonzero", str(cm.exception))

    def test_green_blocks_within_roundoff_are_accepted(self):
        """The Green check uses the same relative allowance as the chi check, so
        a last-bit difference does not reject a paramagnetic producer."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, green_blocks=2)
            g = np.load(os.path.join(d, "green.npz"))["green"]
            g[1] = g[0] * (1.0 + 0.1 * sc._SPIN_DISCARD_ROUNDOFF_RATIO)
            np.savez(os.path.join(d, "green.npz"), green=g)
            sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)

    def test_two_identical_green_blocks_are_accepted(self):
        """A paramagnetic run may legitimately store two identical blocks."""
        self._run("_load_flex_susceptibilities", green_blocks=2)

    def test_unequal_green_spin_blocks_are_refused(self):
        """chi can be redundant while the Green functions are not; the pair
        bubble is built from a single propagator, so this must not pass."""
        with self.assertRaises(ValueError) as cm:
            self._run("_load_flex_susceptibilities", green_blocks=2,
                      green_down_factor=0.5)
        msg = str(cm.exception)
        self.assertIn("spin blocks that are not identical", msg)
        self.assertIn("paramagnetic", msg)


if __name__ == "__main__":
    unittest.main()
