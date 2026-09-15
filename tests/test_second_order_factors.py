"""Off-site density-slot arrays and the sparse factor pack (spec 2.4, 2.5):
orientation against the ED-validated pair-space transpose, the spin
structure per type, Hermiticity, a full-period displacement, the on-site
factor matrices against Gamma, and the byte count."""
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

_IN2 = "tests/rpa/input_2orb"
_SHAPE = (4, 4, 1)


def _write_wan(path, name, norb, rows):
    rvecs = sorted({tuple(r[:3]) for r in rows})
    with open(path, "w") as fw:
        fw.write("{} in wannier90-like format for uhfk\n{}\n{}\n".format(name, norb, len(rvecs)))
        fw.write(" ".join("1" for _ in rvecs) + "\n")
        for r in rows:
            fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {: .15e} {: .15e}\n".format(*r))


def _split_for(rows_by_type, norb=2):
    """A FLEX general solver on input_2orb's geometry/transfer with the
    given off-site rows; returns (solver, split)."""
    from hwave.solver.offsite import split_locality
    import hwave.solver.flex as flex_mod
    d = tempfile.mkdtemp()
    for f in ("geom.dat", "transfer.dat"):
        shutil.copy(os.path.join(_IN2, f), d)
    idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for t, rows in rows_by_type.items():
        _write_wan(os.path.join(d, t.lower() + ".dat"), t, norb, rows)
        idict[t] = t.lower() + ".dat"
    r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
    par = {"T": 2.0, "filling": 0.5, "CellShape": list(_SHAPE), "SubShape": [1, 1, 1], "Nmat": 8,
           "IterationMax": 1, "Mix": 1.0, "EPS": 1, "flex_second_order": "takimoto"}
    s = flex_mod.FLEX(r.get_param("ham"), {}, {"mode": "FLEX", "param": par,
                                                "enable_spin_orbital": False, "calc_scheme": "general"})
    return s, split_locality(s.ham_info, s.lattice)


class TestOffsite(unittest.TestCase):

    def test_orientation_matches_the_sc_builder_on_an_asymmetric_bond(self):
        """v_ab(R) != v_ba(R): vpair[(aa),(bb)](q) must equal the S/C channel-0
        density slot C_{(aa),(bb)}(q) / 2 of the Tier-1 builder."""
        from hwave.solver.second_order import build_offsite
        from hwave.solver.offsite import sc_matrices_from_split
        rows = [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0),   # (a=0,b=1,+x) and its closure
                (0, 1, 0, 1, 1, 0.2, 0.0), (0, -1, 0, 1, 1, 0.2, 0.0)]
        s, split = _split_for({"CoulombInter": rows})
        vp = build_offsite(split, s.lattice, 2)
        nx, ny, nz = _SHAPE
        S0, C0 = sc_matrices_from_split(split, ("CoulombInter", "Hund", "Ising"), 2, nx, ny, nz)
        C0 = C0.reshape(nx * ny * nz, 4, 4)
        for a in range(2):
            for b in range(2):
                np.testing.assert_allclose(vp[0, 0, :, a, b], 0.5 * C0[:, a * 2 + a, b * 2 + b],
                                           rtol=0, atol=1e-13)
        # asymmetric: (0,1) and (1,0) differ at generic q
        self.assertGreater(np.max(np.abs(vp[0, 0, :, 0, 1] - vp[0, 0, :, 1, 0])), 1e-3)
        # Hermitian per q, all four spin pairs equal for CoulombInter
        for s1 in (0, 1):
            for s2 in (0, 1):
                np.testing.assert_allclose(vp[s1, s2], vp[0, 0], atol=1e-14)
        np.testing.assert_allclose(vp[0, 0], np.conj(np.swapaxes(vp[0, 0], -1, -2)), atol=1e-13)

    def test_spin_structure_hund_ising(self):
        from hwave.solver.second_order import build_offsite
        rows = [(1, 0, 0, 1, 1, 0.5, 0.0), (-1, 0, 0, 1, 1, 0.5, 0.0)]
        s, split = _split_for({"Hund": rows})
        vp = build_offsite(split, s.lattice, 2)
        np.testing.assert_allclose(vp[0, 0], -vp_ref(s, split, rows), atol=1e-13)
        self.assertEqual(np.abs(vp[0, 1]).max(), 0.0)
        s, split = _split_for({"Ising": rows})
        vp = build_offsite(split, s.lattice, 2)
        np.testing.assert_allclose(vp[0, 0], vp_ref(s, split, rows), atol=1e-13)
        np.testing.assert_allclose(vp[0, 1], -vp_ref(s, split, rows), atol=1e-13)

    def test_full_period_displacement_is_offsite_and_q_independent(self):
        from hwave.solver.second_order import build_offsite
        rows = [(4, 0, 0, 1, 1, 0.5, 0.0), (-4, 0, 0, 1, 1, 0.5, 0.0)]
        s, split = _split_for({"CoulombInter": rows})
        vp = build_offsite(split, s.lattice, 2)
        self.assertIsNotNone(vp)
        np.testing.assert_allclose(vp[0, 0, :, 0, 0], vp[0, 0, 0, 0, 0], atol=1e-13)
        self.assertAlmostEqual(vp[0, 0, 0, 0, 0].real, 1.0, places=12)   # 0.5 + 0.5 (both rows)

    def test_none_without_offsite_terms(self):
        from hwave.solver.second_order import build_offsite
        s, split = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0)]})
        self.assertIsNone(build_offsite(split, s.lattice, 2))


def vp_ref(s, split, rows):
    """Reference: the CoulombInter-type placement of the same rows (spin-independent)."""
    from hwave.solver.second_order import build_offsite
    from hwave.solver.offsite import split_locality
    tbl = {"CoulombInter": {((r[0], r[1], r[2]), (r[3] - 1, r[4] - 1)): r[5] for r in rows}}
    class _S:  # a minimal split with only the off-site part
        offsite_tbl = tbl
        offsite_types = ("CoulombInter",)
    return build_offsite(_S, s.lattice, 2)[0, 0]


class TestFactors(unittest.TestCase):

    def test_onsite_factor_matrices_match_gamma(self):
        from hwave.solver.second_order import build_factors, compile_onsite, gen_index
        s, split = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0), (0, 0, 0, 2, 2, 1.0, 0.0)],
                               "CoulombInter": [(0, 0, 0, 1, 2, 0.4, 0.0), (0, 0, 0, 2, 1, 0.4, 0.0)],
                               "Exchange": [(0, 0, 0, 1, 2, 0.2, 0.0), (0, 0, 0, 2, 1, 0.2, 0.0)]})
        f = build_factors(split, s.lattice, 2)
        G = compile_onsite(split.onsite_tbl, 2)
        norb, nd = 2, 4
        self.assertLessEqual(len(f.triples), 8)
        seen = set(f.triples)
        for (ss, sr, sq), A, B in zip(f.triples, f.A_on, f.B_on):
            for c, a, r, q in np.ndindex(norb, norb, norb, norb):
                self.assertAlmostEqual(
                    A[c * norb + a, r * norb + q],
                    G[gen_index(a, 0, norb), gen_index(q, sq, norb), gen_index(r, sr, norb), gen_index(c, ss, norb)])
            for rp, qp, d, b in np.ndindex(norb, norb, norb, norb):
                self.assertAlmostEqual(
                    B[rp * norb + qp, d * norb + b],
                    G[gen_index(rp, sr, norb), gen_index(d, ss, norb), gen_index(b, 0, norb), gen_index(qp, sq, norb)])
            self.assertTrue(np.abs(A).max() > 0 or np.abs(B).max() > 0)
        # every dropped triple is structurally zero
        for trip in set(np.ndindex(2, 2, 2)) - seen:
            ss, sr, sq = trip
            A = np.zeros((nd, nd), complex)
            for c, a, r, q in np.ndindex(norb, norb, norb, norb):
                A[c * norb + a, r * norb + q] = G[gen_index(a, 0, norb), gen_index(q, sq, norb),
                                                  gen_index(r, sr, norb), gen_index(c, ss, norb)]
            self.assertEqual(np.abs(A).max(), 0.0)
        self.assertIsNone(f.vpair)
        self.assertEqual(f.nbytes, 16 * (2 * len(f.triples) * nd * nd))

    def test_offsite_factors_are_packed_read_only(self):
        """A split with BOTH parts: vpair is packed into the factors,
        write-locked, and counted in nbytes."""
        from hwave.solver.second_order import build_factors, build_offsite
        rows = [(1, 0, 0, 1, 2, 0.3, 0.0), (-1, 0, 0, 2, 1, 0.3, 0.0),
                (0, 1, 0, 1, 1, 0.2, 0.0), (0, -1, 0, 1, 1, 0.2, 0.0)]
        s, split = _split_for({"CoulombIntra": [(0, 0, 0, 1, 1, 1.0, 0.0),
                                                (0, 0, 0, 2, 2, 1.0, 0.0)],
                               "CoulombInter": rows})
        f = build_factors(split, s.lattice, 2)
        self.assertIsNotNone(f.vpair)
        self.assertEqual(f.vpair.shape, (2, 2, 16, 2, 2))
        self.assertTrue(np.array_equal(f.vpair, build_offsite(split, s.lattice, 2)))
        self.assertFalse(f.vpair.flags.writeable)
        with self.assertRaises(ValueError):
            f.vpair[0, 0, 0, 0, 0] = 1.0
        self.assertEqual(f.nbytes, 16 * (2 * len(f.triples) * 16 + 4 * 16 * 4))

    def test_offsite_exchange_alone_carries_no_density_slot(self):
        """Off-site Exchange is not a density type: no vpair."""
        from hwave.solver.second_order import build_offsite
        s, split = _split_for({"Exchange": [(1, 0, 0, 1, 2, 0.3, 0.0),
                                            (-1, 0, 0, 2, 1, 0.3, 0.0)]})
        self.assertIn("Exchange", split.offsite_tbl)
        self.assertIsNone(build_offsite(split, s.lattice, 2))

    def test_factor_bytes_formula(self):
        from hwave.solver.second_order import factor_bytes
        self.assertEqual(factor_bytes(2, 16, True), 16 * (8 * 16 + 4 * 16 * 4))
        self.assertEqual(factor_bytes(3, 100, False), 16 * 8 * 81)


if __name__ == "__main__":
    unittest.main()
