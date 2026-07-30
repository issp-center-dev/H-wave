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
    to ``"reduced"`` for U + Hund/Ising (auto selects GENERAL for Exchange or
    PairHop since #120) and writes a 2-index chi0q, and ``hwave_sc`` then
    routes to the general S/C formulation because an inter-orbital vertex term
    is present.

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
        # Exchange/PairHop are absent: with a 2-index chi0q they are now
        # REJECTED outright (no density-diagonal vertex content, #120
        # alignment) -- pinned separately below
        for keys in (["CoulombIntra"],
                     ["CoulombIntra", "Hund"],
                     ["CoulombIntra", "CoulombInter", "Hund", "Ising"]):
            res = self._run(keys)
            for pairing, (V4, V2) in res.items():
                with self.subTest(interactions="+".join(keys), pairing=pairing):
                    np.testing.assert_allclose(
                        V2, V4, rtol=1e-9, atol=1e-11,
                        err_msg="a density-diagonal chi0q loses nothing in the "
                                "2-index reduction, so the vertices must match")

    def test_one_orbital_reduced_input_still_rejects_exchange_pairhop(self):
        """norb = 1 must NOT bypass the rejection: the one-orbital S/C
        builder gives Exchange/PairHop zero weight too, so accepting them
        silently omits the interaction (the norb shortcut used to return
        before the check -- round 7)."""
        import hwave.sc as sc

        for itype in ("Exchange", "PairHop"):
            with self.subTest(interaction=itype):
                m = np.full((1, 1, 2, 1, 1), 0.4, dtype=complex)
                with self.assertRaises(ValueError) as cm:
                    sc._warn_reduced_flex_missing_components(
                        {"CoulombIntra": m.copy(), itype: m}, 1, 2, 1, 1)
                self.assertIn(itype, str(cm.exception))

    def test_two_index_chi0q_rejects_exchange_and_pairhop(self):
        """A 2-index chi0q carries no off-density bubble at all, and
        Exchange/PairHop have no density-diagonal vertex: nothing of them
        could be dressed, and since #120 no reduced producer can generate
        such a run -- the combination is rejected, not approximated."""
        with self.assertRaises(ValueError) as cm:
            self._run(["CoulombIntra", "Exchange"])
        self.assertIn("general", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            self._run(["CoulombIntra", "PairHop"])
        self.assertIn("general", str(cm.exception))


class TestReducedFlexMissingComponentWarning(unittest.TestCase):
    """A reduced FLEX chi cannot dress the off-density S/C blocks that
    inter-orbital interactions create, so consuming the two together must say
    so instead of silently returning an approximate lambda."""

    def _inter_k(self, keys, norb=2, Nx=2, Ny=2, Nz=1):
        return {k: np.ones((norb, norb, Nx, Ny, Nz), dtype=complex)
                for k in keys}

    def test_warns_for_each_unsupported_interorbital_term(self):
        import hwave.sc as sc

        # partial types warn; the no-density-content types are rejected
        for key in ("CoulombInter", "Hund", "Ising"):
            with self.subTest(interaction=key):
                with self.assertLogs("hwave_sc", level="WARNING") as cm:
                    sc._warn_reduced_flex_missing_components(
                        self._inter_k(["CoulombIntra", key]), 2, 2, 2, 1)
                joined = "\n".join(cm.output)
                self.assertIn(key, joined)
                self.assertIn("calc_scheme='general'", joined)
        for key in ("Exchange", "PairHop"):
            with self.subTest(interaction=key):
                with self.assertRaises(ValueError) as cm:
                    sc._warn_reduced_flex_missing_components(
                        self._inter_k(["CoulombIntra", key]), 2, 2, 2, 1)
                self.assertIn(key, str(cm.exception))
                self.assertIn("general", str(cm.exception))

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

        Under the adjudicated vertex table (#113) the pair-cross slots carry
        CoulombInter (S, C) = (+1, -1) and Hund (-1, +1) per unit coupling,
        so EQUAL couplings cancel the off-density block exactly in both
        channels. (The test originally used Exchange = -PairHop, which
        cancelled in the pre-#113 case-4 placement; Exchange has since moved
        to the pair-diagonal family and no longer shares a slot with
        PairHop.) A per-term test would announce missing dressing that does
        not exist."""
        import hwave.sc as sc
        import logging

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        jp = np.zeros((norb, norb, Nx, Ny, Nz), dtype=complex)
        jp[0, 1] = 0.6
        jp[1, 0] = 0.6
        inter_k = {"CoulombIntra": np.ones((norb, norb, Nx, Ny, Nz),
                                           dtype=complex),
                   "CoulombInter": jp, "Hund": jp.copy()}
        self.assertEqual(
            sc._off_density_sc_weight(inter_k, norb, Nx, Ny, Nz), 0.0,
            "equal CoulombInter and Hund must cancel the pair-cross "
            "off-density block under the adjudicated table")

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
        # Hund alone: a PARTIAL type (density slots dressed, cross slots
        # not) with genuine off-density weight -- the warning case.
        # Exchange/PairHop are rejected outright now, pinned elsewhere.
        inter_k = {"CoulombIntra": np.ones((norb, norb, Nx, Ny, Nz),
                                           dtype=complex),
                   "Hund": jp}
        self.assertGreater(
            sc._off_density_sc_weight(inter_k, norb, Nx, Ny, Nz), 0.0)
        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            sc._warn_reduced_flex_missing_components(inter_k, norb, Nx, Ny, Nz)
        self.assertIn("Hund", "\n".join(cm.output))

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
            # Filter on nothing: the guard has several messages and an earlier
            # version of this test matched a phrase none of them contain, so the
            # list was empty however the guard behaved.
            warns = [r.getMessage() for r in records
                     if r.levelno >= logging.WARNING]
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
    """The guard lives at the loader boundary, so pin it THERE.

    The helper-level tests above cannot catch the check being dropped from
    _load_flex_susceptibilities / _load_flex_susceptibilities_full, which is
    exactly the placement that matters. Susceptibility only: the dressed Green
    function is not guarded here (tracked separately).
    """

    def _write(self, d, norb=2, nmat=4, Nx=2, Ny=2, Nz=1,
               down_factor=1.0, cross=0.0):
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
        np.savez(os.path.join(d, "green.npz"),
                 green=rc((1, nmat, nvol, norb, norb)))
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


def _bad_tail(rc, nmat, nvol, norb, nd_so):
    """Redundant at the static slice, NOT redundant elsewhere on the axis."""
    block = rc((nmat, nvol, norb, norb))
    a = np.zeros((nmat, nvol, nd_so, nd_so), dtype=complex)
    a[..., :norb, :norb] = block
    a[..., norb:, norb:] = block
    a[0, ..., norb:, norb:] *= 0.5
    return a


class TestGuardsDoNotBreakLegitimateInput(unittest.TestCase):
    """The refusals must not reject workflows that were fine before.

    Each case here is one the guard could plausibly over-reject: data that is
    absent rather than polarized, frequencies that never enter the result, and
    round-off between Green blocks.
    """

    def _write(self, d, norb=2, nmat=4, Nx=2, Ny=2, Nz=1, build=None,
               green=None, tag_up_only=False):
        nvol, nd_so = Nx * Ny * Nz, norb * 2
        rng = np.random.default_rng(11)

        def rc(shape):
            return (rng.standard_normal(shape)
                    + 1j * rng.standard_normal(shape)) * 0.01

        extra = {"chi_spin_blocks": "up_only"} if tag_up_only else {}
        for name, key in (("chiq_s", "chiq_s"), ("chiq_c", "chiq_c")):
            np.savez(os.path.join(d, name + ".npz"),
                     **{key: build(rc, nmat, nvol, norb, nd_so)},
                     chi_convention="kuroki", **extra)
        g = green(rc, nmat, nvol, norb) if green else rc(
            (1, nmat, nvol, norb, norb))
        np.savez(os.path.join(d, "green.npz"), green=g)
        return {"mode": {"param": {"Nmat": nmat}},
                "file": {"output": {"path_to_output": d}},
                "eliashberg": {"chi0q_mode": "flex"}}, norb, Nx, Ny, Nz

    @staticmethod
    def _only_up(rc, nmat, nvol, norb, nd_so):
        a = np.zeros((nmat, nvol, nd_so, nd_so), dtype=complex)
        a[..., :norb, :norb] = rc((nmat, nvol, norb, norb))
        return a

    def test_untagged_all_zero_down_block_is_refused(self):
        """An all-zero down/cross block is NOT a reliable legacy marker: a
        saturated or projected spin sector looks identical. The values alone
        cannot tell them apart, so an untagged file must not be assumed benign
        -- that was a hole in the first version of this guard."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=self._only_up)
            with self.assertRaises(ValueError) as cm:
                sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)
        self.assertIn("chi_spin_blocks='up_only'", str(cm.exception))

    def test_config_flag_accepts_an_untagged_legacy_file(self):
        """Files written before this check existed cannot carry the tag, so the
        configuration route must work for them -- otherwise the escape hatch is
        unreachable for exactly the files that need it."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=self._only_up)
            inp["eliashberg"]["accept_up_block_only"] = True
            with self.assertLogs("hwave_sc", level="WARNING") as cm:
                sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)
        msg = "\n".join(cm.output)
        # Name the real authoriser: a user override and a file assertion carry
        # different weight, and an earlier version reported the override as
        # though the file had been tagged.
        self.assertIn("accept_up_block_only", msg)
        self.assertNotIn("chi_spin_blocks", msg)

    def test_accept_up_block_only_rejects_non_boolean(self):
        """A TOML typo like the string "false" is truthy in Python; silently
        enabling an override that relaxes a correctness guard is the worst
        reading of a malformed value."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=self._only_up)
            inp["eliashberg"]["accept_up_block_only"] = "false"
            with self.assertRaises(ValueError) as cm:
                sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)
        self.assertIn("must be a boolean", str(cm.exception))

    def test_tagged_legacy_file_is_accepted(self):
        """...but a file that declares the layout is read as intended."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=self._only_up,
                                                tag_up_only=True)
            with self.assertLogs("hwave_sc", level="WARNING") as cm:
                sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)
        self.assertIn("chi_spin_blocks", "\n".join(cm.output))

    def test_static_route_also_validates_the_whole_axis(self):
        """Paramagnetism is a property of the producing run, not only of the
        slice consumed. Validating just the static slice let a file that the
        dynamic route rejects through the static one -- and with a one-block
        Green function it would then return a paramagnetic eigenvalue in
        silence."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=_bad_tail)
            with self.assertRaises(ValueError):
                sc._load_flex_susceptibilities(inp, norb, Nx, Ny, Nz)

    def test_dynamic_route_rejects_it_too(self):
        """The two routes must agree about whether a file is usable."""
        import hwave.sc as sc
        with tempfile.TemporaryDirectory() as d:
            inp, norb, Nx, Ny, Nz = self._write(d, build=_bad_tail)
            with self.assertRaises(ValueError):
                sc._load_flex_susceptibilities_full(inp, norb, Nx, Ny, Nz)

    def test_non_finite_susceptibility_is_refused(self):
        """NaN in a DISCARDED block is the dangerous case: max(0.0, nan) is 0.0
        in Python, so it would leave the redundancy maxima untouched and the
        array would be accepted as exactly redundant."""
        import hwave.sc as sc
        for bad in (np.nan, np.inf):
            for blk, sl in (("up", (slice(None), slice(None), 0, 0)),
                            ("down", (slice(None), slice(None), 2, 2)),
                            ("cross", (slice(None), slice(None), 0, 2))):
                with self.subTest(value=bad, block=blk):
                    with tempfile.TemporaryDirectory() as d:
                        inp, norb, Nx, Ny, Nz = self._write(d, build=self._para)
                        f = os.path.join(d, "chiq_s.npz")
                        z = dict(np.load(f))
                        z["chiq_s"][sl] = bad
                        np.savez(f, **z)
                        with self.assertRaises(ValueError) as cm:
                            sc._load_flex_susceptibilities(inp, norb,
                                                           Nx, Ny, Nz)
                    self.assertIn("non-finite", str(cm.exception))

    def _green_pair(self, factor):
        def g(rc, nmat, nvol, norb):
            g0 = rc((nmat, nvol, norb, norb))
            return np.stack([g0, g0 * factor])
        return g

    def _para(self, rc, nmat, nvol, norb, nd_so):
        block = rc((nmat, nvol, norb, norb))
        a = np.zeros((nmat, nvol, nd_so, nd_so), dtype=complex)
        a[..., :norb, :norb] = block
        a[..., norb:, norb:] = block
        return a


class TestChunkedCheckMatchesWholeArray(unittest.TestCase):
    """The guards accumulate per frequency to stay memory-bounded. That rewrite
    must not change any decision, including at the boundaries where it would be
    easiest to get wrong."""

    NORB = 2

    def _whole_array_verdict(self, chi):
        """The decision a whole-array evaluation would reach, written out
        independently of the implementation."""
        n = self.NORB
        up, dn = chi[..., :n, :n], chi[..., n:, n:]
        cu, cd = chi[..., :n, n:], chi[..., n:, :n]
        import hwave.sc as sc
        if float(np.max(np.abs(up - dn))) == 0.0 and \
                float(np.max(np.abs(cu))) == 0.0 and \
                float(np.max(np.abs(cd))) == 0.0:
            return False
        if float(np.max(np.abs(up))) > 0.0 and not np.any(dn) \
                and not np.any(cu) and not np.any(cd):
            return False        # legacy layout, accepted when authorised
        # PER-FREQUENCY ratios. An earlier version of this reference mirrored
        # production's single global difference over global scale, so it agreed
        # with the implementation by construction and could not detect a large
        # redundant frequency masking a small non-redundant one.
        worst = 0.0
        for w in range(chi.shape[0]):
            u, d = up[w], dn[w]
            c1, c2 = cu[w], cd[w]
            ws = float(np.max(np.abs(u)))
            wd = max(float(np.max(np.abs(u - d))),
                     float(np.max(np.abs(c1))), float(np.max(np.abs(c2))))
            worst = max(worst, wd / max(ws, np.finfo(float).tiny))
        return worst > sc._SPIN_DISCARD_ROUNDOFF_RATIO

    def _cases(self):
        import hwave.sc as sc
        n, nd_so = self.NORB, self.NORB * 2
        R = sc._SPIN_DISCARD_ROUNDOFF_RATIO
        rng = np.random.default_rng(0)
        for nfreq in (1, 3, 8):
            base = (rng.standard_normal((nfreq, 5, n, n))
                    + 1j * rng.standard_normal((nfreq, 5, n, n)))

            def mk(down=1.0, cross=0.0, zero_down=False):
                a = np.zeros((nfreq, 5, nd_so, nd_so), dtype=complex)
                a[..., :n, :n] = base
                if not zero_down:
                    a[..., n:, n:] = base * down
                if cross:
                    a[..., :n, n:] = base * cross
                    a[..., n:, :n] = base * cross
                return a

            yield "n%d paramagnetic" % nfreq, mk()
            yield "n%d polarized" % nfreq, mk(down=0.5)
            yield "n%d cross" % nfreq, mk(cross=0.3)
            yield "n%d zero-down" % nfreq, mk(zero_down=True)
            yield "n%d just-below" % nfreq, mk(down=1.0 + 0.1 * R)
            yield "n%d just-above" % nfreq, mk(down=1.0 + 10.0 * R)
            if nfreq > 1:
                c = mk()
                c[0, ..., n:, n:] *= 0.5     # only an unused frequency differs
                yield "n%d tail-only" % nfreq, c

                # A huge redundant frequency next to a small, totally
                # non-redundant one. Judged globally the mismatch looks like
                # 1e-14 and slips under the threshold; judged per frequency it
                # is 100%.
                m = np.zeros((nfreq, 5, nd_so, nd_so), dtype=complex)
                m[0, :, :n, :n] = 1e14
                m[0, :, n:, n:] = 1e14
                m[1, :, :n, :n] = 1.0        # down block left at zero
                m[1, :, n:, n:] = 0.0
                m[1, :, 0, n] = 1e-30        # keep it out of the legacy branch
                yield "n%d masked-by-scale" % nfreq, m

                # A frequency whose KEPT block is exactly zero. Redundancy
                # cannot be established there, so nonzero discarded content
                # must be refused however small it is -- max(w_scale, tiny)
                # makes the ratio enormous, which is the intended call, not an
                # accident.
                z = mk()
                z[1, :, :n, :n] = 0.0
                z[1, :, n:, n:] = 1e-300
                yield "n%d zero-kept-scale" % nfreq, z

                # ...and an all-zero frequency is fine: nothing is discarded.
                a0 = mk()
                a0[1] = 0.0
                yield "n%d all-zero-frequency" % nfreq, a0

    def test_every_case_decides_the_same_way(self):
        import hwave.sc as sc
        import logging
        # Silence the guard's warnings for this test only. Setting the level
        # without restoring it leaks into every later test that asserts on
        # warnings -- which is exactly what happened the first time.
        lg = logging.getLogger("hwave_sc")
        prev = lg.level
        lg.setLevel(logging.CRITICAL)
        self.addCleanup(lg.setLevel, prev)
        for name, chi in self._cases():
            with self.subTest(case=name):
                want_raise = self._whole_array_verdict(chi)
                try:
                    sc._check_spin_block_discarded(chi, self.NORB, "kuroki",
                                                   "x", True)
                    got_raise = False
                except ValueError:
                    got_raise = True
                self.assertEqual(got_raise, want_raise)


class TestSpinGuardDiagnostics(unittest.TestCase):
    """Acceptance is a per-frequency ratio, so the message must name the
    frequency that decided it -- quoting a global maximum could pair a small
    tail's mismatch with a huge unrelated scale from elsewhere on the axis."""

    NORB = 2

    def _chi(self, nfreq, mismatch_at, scales):
        n, nd_so = self.NORB, self.NORB * 2
        a = np.zeros((nfreq, 3, nd_so, nd_so), dtype=complex)
        for w in range(nfreq):
            a[w, :, :n, :n] = scales[w]
            a[w, :, n:, n:] = scales[w]
        for w in mismatch_at:
            a[w, :, n:, n:] = 0.0
        return a

    def _raise_message(self, chi):
        import hwave.sc as sc
        with self.assertRaises(ValueError) as cm:
            sc._check_spin_block_discarded(chi, self.NORB, "kuroki", "chi_s")
        return str(cm.exception)

    def test_sole_mismatch_at_frequency_zero_is_named(self):
        msg = self._raise_message(self._chi(3, [0], [1.0, 1e6, 1e6]))
        self.assertIn("frequency index 0", msg)

    def test_sole_mismatch_in_a_small_tail_reports_its_own_scale(self):
        """The decisive frequency has scale 1.0 while the axis maximum is 1e6;
        reporting the global scale would understate the mismatch by 1e6."""
        msg = self._raise_message(self._chi(3, [2], [1e6, 1e6, 1.0]))
        self.assertIn("frequency index 2", msg)
        self.assertIn("1.000e+00", msg)
        self.assertNotIn("1.000e+06", msg)

    def test_tied_ratios_still_report_one_of_them(self):
        msg = self._raise_message(self._chi(3, [1, 2], [1.0, 5.0, 5.0]))
        self.assertTrue("frequency index 1" in msg or "frequency index 2" in msg)


class TestMalformedShapesAreRefused(unittest.TestCase):
    """A rectangle whose trailing elements happen to number a perfect square
    used to be reshaped into a susceptibility without complaint."""

    def test_rectangular_trailing_axes_are_refused(self):
        import hwave.sc as sc
        norb = 2
        # 2*8 = 16 trailing elements; sqrt(16) = 4 = norb^2, so deriving nd_chi
        # from the element count accepted this as a 4x4.
        chi = np.arange(1 * 1 * 2 * 8, dtype=complex).reshape(1, 1, 2, 8)
        with self.assertRaises(ValueError) as cm:
            sc._expand_flex_chi(chi, norb, 1, 1, 1, "kuroki")
        self.assertIn("shape", str(cm.exception))
        with self.assertRaises(ValueError):
            sc._check_spin_block_discarded(chi, norb, "kuroki", "chi_s")

    def test_the_two_supported_layouts_still_load(self):
        """The refusal must not catch either real layout: reduced is
        (nfreq, nvol, nd, nd) and general is
        (nfreq, nvol, norb, norb, norb, norb).

        Distinct nonzero values, not zeros: the rank-6 branch flattens the pair
        axes, and an all-zero array would agree on shape however that flattening
        permuted the values."""
        import hwave.sc as sc
        for norb in (1, 2, 3):
            with self.subTest(norb=norb, layout="reduced"):
                red = np.zeros((2, 1, norb * 2, norb * 2), dtype=complex)
                out = sc._expand_flex_chi(red, norb, 1, 1, 1, "kuroki")
                self.assertEqual(out.shape[-2:],
                                 (norb * norb, norb * norb))
            with self.subTest(norb=norb, layout="general"):
                n = norb
                gen = (np.arange(2 * 1 * n ** 4, dtype=float)
                       + 1j * np.arange(2 * 1 * n ** 4, dtype=float)
                       ).reshape(2, 1, n, n, n, n)
                out = sc._expand_flex_chi(gen, norb, 1, 1, 1, "myo")
                # pins the ORDER, not just the shape: row pair (a,c), col (b,d)
                np.testing.assert_array_equal(
                    out, gen.reshape(2, 1, 1, 1, n * n, n * n))


class TestGreenTwoBlockGuard(unittest.TestCase):
    """A two-block green.npz feeds the pair bubble from block 0 only;
    differing blocks are real physics and must be rejected unless the
    user takes responsibility (mirrors the chi-side policy)."""

    def _write(self, d, blocks, chi_extra=None):
        import hwave.sc as sc

        norb, Nx, Ny, Nz, nmat = 1, 2, 1, 1, 4
        nvol = Nx * Ny * Nz
        nd_so = norb * 2
        chi = np.zeros((nmat, nvol, nd_so, nd_so), dtype=complex)
        for b in range(2):
            chi[:, :, b*norb:(b+1)*norb, b*norb:(b+1)*norb] = 0.1
        meta = dict(chi_extra or {})
        np.savez(os.path.join(d, 'chiq_s.npz'), chiq_s=chi, **meta)
        np.savez(os.path.join(d, 'chiq_c.npz'), chiq_c=chi, **meta)
        green = np.stack(blocks, axis=0)
        np.savez(os.path.join(d, 'green.npz'), green=green)
        return {"file": {"output": {"path_to_output": d}},
                "eliashberg": {}}

    def _blocks(self, second_scale=1.0, nan=False):
        norb, nvol, nmat = 1, 2, 4
        g = (np.arange(nmat * nvol * norb * norb, dtype=complex)
             .reshape(nmat, nvol, norb, norb) + 1.0)
        if nan:
            g[0, 0, 0, 0] = np.nan
        if second_scale == 1.0:
            # a true copy: multiplying by 1.0 is a COMPLEX multiply and
            # turns nan+0j into nan+nanj, which under the exact
            # per-component semantics is genuinely differing content
            return [g, g.copy()]
        return [g, g * second_scale]

    def test_identical_blocks_pass(self):
        import tempfile

        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d, self._blocks(1.0))
        green = sc._load_flex_green(inp, 1, 2, 1, 1)
        self.assertIsNotNone(green)

    def test_nonfinite_green_is_rejected_at_the_boundary(self):
        """NaN/Inf in the dressed Green function reaches the pair bubble
        directly and can surface as a saved zero eigenvalue; the loader
        rejects it up front, as the susceptibility loader does -- for
        one-block and (block-redundant) two-block files alike."""
        import tempfile

        import hwave.sc as sc

        for label, val, blocks in (
                ("nan two-block", np.nan, 2),
                ("inf two-block", np.inf, 2),
                ("nan one-block", np.nan, 1)):
            with self.subTest(case=label):
                b = self._blocks(1.0)
                b[0][0, 0, 0, 0] = val
                b[1] = b[0].copy()
                if blocks == 1:
                    b = [b[0]]
                d = tempfile.mkdtemp()
                inp = self._write(d, b)
                with self.assertRaises(ValueError) as cm:
                    sc._load_flex_green(inp, 1, 2, 1, 1)
                self.assertIn('non-finite', str(cm.exception))

    def test_partial_nan_components_are_rejected_as_nonfinite(self):
        """Entries like nan+1j (one finite component) are caught by the
        finiteness boundary before any block comparison -- the historical
        complex equal_nan hole (either-component NaN counted as equal)
        is thereby closed at the boundary."""
        import tempfile

        import hwave.sc as sc

        for kind in ("imag-differs", "real-differs"):
            with self.subTest(kind=kind):
                blocks = self._blocks(1.0)
                if kind == "imag-differs":
                    blocks[0][0, 0, 0, 0] = complex(np.nan, 1.0)
                    blocks[1][0, 0, 0, 0] = complex(np.nan, 2.0)
                else:
                    blocks[0][0, 0, 0, 0] = complex(1.0, np.nan)
                    blocks[1][0, 0, 0, 0] = complex(2.0, np.nan)
                d = tempfile.mkdtemp()
                inp = self._write(d, blocks)
                with self.assertRaises(ValueError) as cm:
                    sc._load_flex_green(inp, 1, 2, 1, 1)
                self.assertIn('non-finite', str(cm.exception))

    def test_differing_blocks_are_rejected(self):
        import tempfile

        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d, self._blocks(2.0))
        with self.assertRaises(ValueError) as cm:
            sc._load_flex_green(inp, 1, 2, 1, 1)
        self.assertIn('accept_up_block_only', str(cm.exception))

    def test_differing_blocks_pass_with_authorization(self):
        import logging
        import tempfile

        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d, self._blocks(2.0))
        inp['eliashberg']['accept_up_block_only'] = True
        with self.assertLogs('hwave_sc', level=logging.WARNING) as cm:
            green = sc._load_flex_green(inp, 1, 2, 1, 1)
        self.assertIsNotNone(green)
        self.assertTrue(any('DIFFERING' in m for m in cm.output))

    def test_comparison_is_sliced_and_short_circuits(self):
        """The guard must compare one-frequency slices (bounded memory)
        and stop at the first mismatch (early exit): instrument
        np.array_equal and check the largest array it receives and the
        call count when the FIRST frequency already differs."""
        import tempfile
        from unittest import mock

        import hwave.sc as sc

        blocks = self._blocks(1.0)
        blocks[1][0, 0, 0, 0] += 1.0     # first frequency differs
        d = tempfile.mkdtemp()
        inp = self._write(d, blocks)
        calls = []
        real_eq = np.array_equal

        def spy(a, b, equal_nan=False):
            calls.append(np.asarray(a).size)
            return real_eq(a, b, equal_nan=equal_nan)

        with mock.patch.object(np, 'array_equal', side_effect=spy):
            with self.assertRaises(ValueError):
                sc._load_flex_green(inp, 1, 2, 1, 1)
        nvol_slice = 2 * 1 * 1           # (nvol, norb, norb) elements
        self.assertTrue(calls, "the guard must compare something")
        self.assertLessEqual(max(calls), nvol_slice,
                             "comparison must be per-frequency slices")
        self.assertLessEqual(len(calls), 2,
                             "a first-frequency mismatch must short-circuit")

    def test_dynamic_route_rejects_differing_blocks(self):
        import tempfile

        import hwave.sc as sc

        d = tempfile.mkdtemp()
        inp = self._write(d, self._blocks(2.0))
        with self.assertRaises(ValueError):
            sc._load_flex_green(inp, 1, 2, 1, 1, allow_ir=True)


if __name__ == "__main__":
    unittest.main()
