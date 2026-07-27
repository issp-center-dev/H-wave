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


def _make_reduced_flex():
    """Build a spin-free ``calc_scheme="reduced"`` FLEX solver on the fixture."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex

    dirpath = os.path.join(tempfile.gettempdir(), "hwave_flex_2d_2orb_intra")
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
    return solver, dirpath


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

    def _build(self):
        import hwave.sc as sc

        flex, dirpath = _make_reduced_flex()
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
                        self._inter_k(["CoulombIntra", key]), norb=2)
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
                self._inter_k(["CoulombIntra"]), norb=2)
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
                self._inter_k(["CoulombIntra", "Hund"], norb=1), norb=1)
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
                                      pairing_type="singlet", convention="myo")
        finally:
            lg.removeHandler(handler)
        warns = [r.getMessage() for r in records
                 if r.levelno >= logging.WARNING]
        self.assertEqual(warns, [])

    def test_kuroki_convention_warns_through_compute_vertices_flex(self):
        """The warning must reach a real caller, not just the helper."""
        import hwave.sc as sc

        norb, Nx, Ny, Nz = 2, 2, 2, 1
        nd = norb * norb
        inter_k = self._inter_k(["CoulombIntra", "Hund"], norb, Nx, Ny, Nz)
        chi = np.zeros((Nx, Ny, Nz, nd, nd), dtype=complex)

        with self.assertLogs("hwave_sc", level="WARNING") as cm:
            sc._compute_vertices_flex(chi, chi, inter_k, norb, Nx, Ny, Nz,
                                      pairing_type="singlet",
                                      convention="kuroki")
        self.assertIn("Hund", "\n".join(cm.output))


if __name__ == "__main__":
    unittest.main()
