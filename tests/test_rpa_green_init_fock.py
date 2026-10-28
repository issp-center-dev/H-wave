"""RPA's ``green_init`` band is UHFk's ``T + Hartree[G] + Fock[G]`` (issue #278).

``RPA._calc_trans_mod`` assembles the one-body band of a ``[file.input]
green_init`` Green function. It used to add the Hartree term only, so the
band differed from the ``trans_mod`` UHFk writes for the same G wherever a
Fock term contributes (CoulombInter, Hund, Ising, Exchange, off-site PairHop).

Oracle (independent of the RPA code under test): UHFk's own Hamiltonian
construction ``UHFk._make_ham`` (Fock on) applied to the same G -- directly
on the k-space Hamiltonian for the plain lattice, and through the files UHFk
writes (``_save_trans_mod`` / ``_save_green``, read back by
``RPA.read_init``) for the sublattice and spin-orbital set-ups. The test
G is a generic Hermitian real-space Green function with off-site content, in
three spin structures: spin-independent, spin-diagonal with unequal spins,
and spin-mixing.

UHFk projects its Hamiltonian onto the block structure it detects from T and
the interaction tables; a UHF solution has a G with that block structure, so
the projection is a no-op there. The generic test G does not, so the oracle
runs ``_make_ham`` with a single block (no projection).

Tests must be run from the repository root.
"""

import os
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as sol_rpa
import hwave.solver.uhfk as sol_uhfk
from hwave.solver import spin_orbital

LOG = {"print_level": 0, "print_step": 1}

# per row kind: (rx, ry, a, b, value) on two orbitals
ROW_KINDS = {
    "onsite_intra": [(0, 0, 0, 0, 0.7), (0, 0, 1, 1, 0.5)],
    "onsite_inter": [(0, 0, 0, 1, 0.3), (0, 0, 1, 0, 0.3)],
    "offsite_intra": [(1, 0, 0, 0, 0.3), (-1, 0, 0, 0, 0.3)],
    "offsite_inter": [(1, 0, 0, 1, 0.3), (-1, 0, 1, 0, 0.3),
                      (0, 1, 1, 1, 0.2), (0, -1, 1, 1, 0.2)],
}
TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Ising", "Exchange",
         "PairLift", "PairHop", "Coulomb")
# (type, row kind) whose Fock term (UHFk flag_fock) is nonzero for a
# spin-independent G (the issue #278 table); used to show that the
# comparison sees the Fock term. Coulomb's on-site intra-orbital rows are its
# CoulombIntra part (no Fock term).
def _fock_carrying(typ, kind):
    if typ == "Coulomb":
        return kind != "onsite_intra"
    return typ in ("CoulombInter", "Hund", "Ising", "Exchange")


# ---------------------------------------------------------------------------
# model files
# ---------------------------------------------------------------------------

def _write_table(path, title, norb, rows):
    """H-wave wannier90-like table; rows (rx, ry, a, b, value), 0-based
    orbitals. Rows on the same (R, a, b) are summed by the reader."""
    rvecs = sorted({(r[0], r[1], 0) for r in rows})
    with open(path, "w") as f:
        f.write("{}\n{}\n{}\n".format(title, norb, len(rvecs)))
        for i in range(0, len(rvecs), 15):
            f.write(" ".join("1" for _ in rvecs[i:i + 15]) + "\n")
        for (rx, ry, a, b, v) in rows:
            v = complex(v)
            f.write("{} {} 0 {} {} {!r} {!r}\n".format(
                rx, ry, a + 1, b + 1, v.real, v.imag))


def _write_geometry(path, norb):
    with open(path, "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n{}\n".format(norb))
        for a in range(norb):
            f.write("0.0 0.0 0.0\n" if a == 0 else "0.5 0.5 0.0\n")


def _transfer_rows(norb):
    rows = []
    for a in range(norb):
        for (rx, ry) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            rows.append((rx, ry, a, a, -1.0))
    if norb == 2:
        rows += [(0, 0, 0, 1, 0.3 + 0.2j), (0, 0, 1, 0, 0.3 - 0.2j),
                 (1, 0, 0, 1, 0.1 + 0.05j), (-1, 0, 1, 0, 0.1 - 0.05j)]
    return rows


def _so_transfer_rows(norb_phys):
    """Spin-orbital transfer, interleaved index 2*orb + spin, with an
    on-site spin-orbit (spin-mixing) term."""
    rows = []
    for (rx, ry, a, b, v) in _transfer_rows(norb_phys):
        for s in range(2):
            rows.append((rx, ry, 2 * a + s, 2 * b + s, v))
    soc = 0.15j
    rows += [(0, 0, 0, 3, soc), (0, 0, 3, 0, np.conj(soc)),
             (0, 0, 1, 2, -soc), (0, 0, 2, 1, np.conj(-soc))]
    return rows


def _model(d, typ_rows, *, so=False, norb=2):
    """Write geometry, transfer and the interaction tables ``typ_rows``
    ({type: rows}) into ``d``; returns the [file.input] dict."""
    _write_geometry(os.path.join(d, "geom.dat"), 2 * norb if so else norb)
    if so:
        _write_table(os.path.join(d, "transfer.dat"), "transfer", 2 * norb,
                     _so_transfer_rows(norb))
    else:
        _write_table(os.path.join(d, "transfer.dat"), "transfer", norb,
                     _transfer_rows(norb))
    inter = {"path_to_input": d, "Geometry": "geom.dat",
             "Transfer": "transfer.dat"}
    for typ, rows in typ_rows.items():
        name = typ.lower() + ".dat"
        _write_table(os.path.join(d, name), typ, norb, rows)
        inter[typ] = name
    return {"path_to_input": d, "interaction": inter}


# ---------------------------------------------------------------------------
# Green functions
# ---------------------------------------------------------------------------

def _green(shape, norb, spin, seed=7):
    """Real-space G (nvol, 2*norb, 2*norb), spin-block order, Hermitian in
    the sense G_ab(r) = conj(G_ba(-r)), with on-site and off-site content.
    ``spin``: "free" (up == down, no mixing), "diag" (unequal spins, no
    mixing) or "mixing" (general)."""
    nx, ny, nz = shape
    nvol, nd = nx * ny * nz, 2 * norb
    rng = np.random.default_rng(seed)
    a = (rng.normal(size=(nvol, nd, nd))
         + 1j * rng.normal(size=(nvol, nd, nd))) * 0.1
    a[0] += 0.4 * np.eye(nd)
    if spin != "mixing":
        a = a.reshape(nvol, 2, norb, 2, norb)
        a[:, 0, :, 1, :] = 0.0
        a[:, 1, :, 0, :] = 0.0
        if spin == "free":
            a[:, 1, :, 1, :] = a[:, 0, :, 0, :]
        a = a.reshape(nvol, nd, nd)
    # r -> -r on the periodic grid
    idx = np.arange(nvol).reshape(nx, ny, nz)
    rev = np.flip(np.roll(idx, -1, axis=(0, 1, 2)), axis=(0, 1, 2)).ravel()
    return 0.5 * (a + np.conj(np.swapaxes(a[rev], -1, -2)))


# ---------------------------------------------------------------------------
# solvers
# ---------------------------------------------------------------------------

def _uhfk(fin, cell, sub=(1, 1, 1), so=False):
    mode = {"mode": "UHFk", "flag_fock": True,
            "param": {"T": 0.1, "CellShape": list(cell), "SubShape": list(sub),
                      "Ncond": 2, "IterationMax": 1}}
    if so:
        mode["enable_spin_orbital"] = True
    u = sol_uhfk.UHFk(read_input_k.QLMSkInput(fin).get_param("ham"), LOG, mode)
    u._make_ham_trans()
    u._make_ham_inter()
    u._detect_blocks()
    u.block_info = [np.arange(u.nd)]      # no block projection (see module doc)
    return u


def _rpa(fin, cell, sub=(1, 1, 1), so=False, orbital_site=None):
    param = {"T": 0.1, "CellShape": list(cell), "SubShape": list(sub),
             "Nmat": 8, "filling": 0.5}
    if orbital_site is not None:
        param["orbital_site"] = list(orbital_site)
    mode = {"mode": "RPA", "calc_scheme": "general", "param": param}
    if so:
        mode["enable_spin_orbital"] = True
    return sol_rpa.RPA(read_input_k.QLMSkInput(fin).get_param("ham"), LOG, mode)


def _max_diff(a, b):
    return float(np.max(np.abs(np.asarray(a) - np.asarray(b))))


CELL = (4, 3, 1)
SUB_CELL = (8, 4, 1)          # folded to (4, 4, 1) by SubShape (2, 1, 1)
TOL = 1e-12


class TestGreenInitBandIsUHFkBand(unittest.TestCase):
    """Plain lattice: RPA's green_init band equals UHFk._make_ham (Fock on)
    of the same G, for every interaction type, row kind, spin structure and
    orbital_site set-up."""

    def test_each_type_row_kind_and_spin_structure(self):
        nvol, nd = CELL[0] * CELL[1], 4
        checked, fock_seen = 0, 0
        for typ in TYPES:
            for kind, rows in ROW_KINDS.items():
                with tempfile.TemporaryDirectory() as d:
                    fin = _model(d, {typ: rows})
                    try:
                        u = _uhfk(fin, CELL)
                    except (ValueError, SystemExit):
                        continue          # declaration refused by the reader
                    for spin in ("free", "diag", "mixing"):
                        G = _green(CELL, 2, spin)
                        u.Green = G.reshape(nvol, 2, 2, 2, 2).copy()
                        u.iflag_fock = True
                        u._make_ham()
                        ref = np.array(u.ham, copy=True)
                        u.iflag_fock = False
                        u._make_ham()
                        hartree_only = np.array(u.ham, copy=True)
                        for site in (None, (0, 1)):
                            r = _rpa(fin, CELL, orbital_site=site)
                            band = r._calc_trans_mod(G.copy())
                            with self.subTest(type=typ, rows=kind, spin=spin,
                                              site=site):
                                self.assertLess(_max_diff(band, ref), TOL)
                                # Hermitian band
                                self.assertLess(_max_diff(
                                    band, np.conj(np.swapaxes(band, -1, -2))), TOL)
                            checked += 1
                        if spin == "free" and _fock_carrying(typ, kind):
                            gap = _max_diff(ref, hartree_only)
                            self.assertGreater(gap, 1e-3, (typ, kind))
                            fock_seen += 1
        self.assertGreater(checked, 150)
        self.assertGreater(fock_seen, 15)

    def test_offsite_pairhop(self):
        """Off-site PairHop: RPA's vertex drops it (issue #157), but UHFk's
        mean field keeps it, and so does the green_init band (before #278
        the band carried no off-site PairHop term at all)."""
        nvol = CELL[0] * CELL[1]
        G = _green(CELL, 2, "free")
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, {"PairHop": ROW_KINDS["offsite_inter"]})
            u = _uhfk(fin, CELL)
            u.Green = G.reshape(nvol, 2, 2, 2, 2).copy()
            u._make_ham()
            r = _rpa(fin, CELL)
            band = r._calc_trans_mod(G.copy())
            bare = r._calc_trans_mod(np.zeros_like(G))
            self.assertLess(_max_diff(band, u.ham), TOL)
            self.assertGreater(_max_diff(band, bare), 1e-3)

    def test_all_types_together(self):
        nvol = CELL[0] * CELL[1]
        typ_rows = {t: ROW_KINDS["offsite_inter"] + ROW_KINDS["onsite_inter"]
                    for t in ("CoulombInter", "Hund", "Ising", "Exchange",
                              "PairLift", "PairHop")}
        typ_rows["CoulombIntra"] = ROW_KINDS["onsite_intra"]
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, typ_rows)
            u = _uhfk(fin, CELL)
            r = _rpa(fin, CELL)
            for spin in ("free", "diag", "mixing"):
                G = _green(CELL, 2, spin, seed=11)
                u.Green = G.reshape(nvol, 2, 2, 2, 2).copy()
                u._make_ham()
                with self.subTest(spin=spin):
                    self.assertLess(_max_diff(r._calc_trans_mod(G.copy()),
                                              u.ham), TOL)

    def test_inputs_not_modified(self):
        for so in (False, True):
            with self.subTest(spin_orbital=so):
                self._inputs_not_modified(so)

    def _inputs_not_modified(self, so):
        # the spin-orbital branch used to add the mean field into a view of
        # ham_trans_r, so a second call saw it twice
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, {"CoulombInter": ROW_KINDS["offsite_inter"]}, so=so)
            r = _rpa(fin, CELL, so=so)
            G = _green(CELL, 2, "free")
            g_before = G.copy()
            t_before = np.array(r.ham_info.ham_trans_r, copy=True)
            a = r._calc_trans_mod(G)
            b = r._calc_trans_mod(G)
            np.testing.assert_array_equal(G, g_before)
            np.testing.assert_array_equal(r.ham_info.ham_trans_r, t_before)
            np.testing.assert_array_equal(a, b)


def _write_uhfk_files(u, d):
    u._make_ham()
    u._save_trans_mod(os.path.join(d, "tm.npz"))
    u._save_green(os.path.join(d, "g.npz"))
    return {"path_to_input": d, "trans_mod": "tm.npz", "green_init": "g.npz"}


ALL_TYPES = {
    "CoulombIntra": ROW_KINDS["onsite_intra"],
    "CoulombInter": ROW_KINDS["offsite_inter"] + ROW_KINDS["onsite_inter"],
    "Hund": ROW_KINDS["offsite_inter"] + ROW_KINDS["onsite_inter"],
    "Ising": ROW_KINDS["offsite_intra"] + ROW_KINDS["onsite_inter"],
    "Exchange": ROW_KINDS["offsite_inter"] + ROW_KINDS["onsite_inter"],
    "PairLift": ROW_KINDS["onsite_inter"],
    "PairHop": ROW_KINDS["offsite_inter"] + ROW_KINDS["onsite_inter"],
}


class TestGreenInitEqualsTransModThroughFiles(unittest.TestCase):
    """The files UHFk writes for one G -- trans_mod and green -- read back
    by RPA define the same band (sublattice fold and spin-orbital remap
    included)."""

    def _check(self, cell, sub, so, typ_rows):
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, typ_rows, so=so)
            u = _uhfk(fin, cell, sub, so=so)
            r = _rpa(fin, cell, sub, so=so)
            norb_orig = 2
            for spin in ("free", "diag", "mixing"):
                G0 = _green(cell, norb_orig, spin, seed=3)   # original lattice
                cellvol = cell[0] * cell[1] * cell[2]
                if so:
                    Gi = spin_orbital.spinblock_to_interleaved(G0, norb_orig)
                    Gu = Gi.reshape(cellvol, 1, 2 * norb_orig, 1, 2 * norb_orig)
                else:
                    Gu = G0.reshape(cellvol, 2, norb_orig, 2, norb_orig)
                if u.has_sublattice:
                    Gu = u._reshape_green(Gu)
                u.Green = np.ascontiguousarray(Gu)
                info = r.read_init(_write_uhfk_files(u, d))
                band = r._calc_trans_mod(info["green_init"])
                with self.subTest(spin=spin):
                    self.assertLess(_max_diff(band, info["trans_mod"]), TOL)

    def test_sublattice(self):
        self._check(SUB_CELL, (2, 1, 1), False, ALL_TYPES)

    def test_sublattice_each_type(self):
        for typ, rows in ALL_TYPES.items():
            with self.subTest(type=typ):
                self._check(SUB_CELL, (2, 1, 1), False, {typ: rows})

    def test_spin_orbital(self):
        self._check(CELL, (1, 1, 1), True, ALL_TYPES)

    def test_spin_orbital_each_type(self):
        for typ, rows in ALL_TYPES.items():
            with self.subTest(type=typ):
                self._check(CELL, (1, 1, 1), True, {typ: rows})

    def test_spin_orbital_sublattice(self):
        self._check(SUB_CELL, (2, 1, 1), True, ALL_TYPES)


class TestEpsilonKUsesTheBand(unittest.TestCase):
    """The band that reaches RPA's H0(k) (``_calc_epsilon_k``) from
    green_init equals the one from the trans_mod UHFk writes for the same G."""

    def test_green_init_and_trans_mod_give_same_h0(self):
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, ALL_TYPES)
            u = _uhfk(fin, CELL)
            G = _green(CELL, 2, "free", seed=5)
            u.Green = G.reshape(-1, 2, 2, 2, 2).copy()
            files = _write_uhfk_files(u, d)
            r = _rpa(fin, CELL)
            info = r.read_init(files)
            r._calc_epsilon_k({"green_init": info["green_init"]})
            h_gi, mode_gi = np.array(r.H0_k, copy=True), r.spin_mode
            r._calc_epsilon_k({"trans_mod": info["trans_mod"]})
            self.assertEqual(mode_gi, r.spin_mode)
            self.assertLess(_max_diff(h_gi, r.H0_k), TOL)


class TestComplexBondOrientation(unittest.TestCase):
    """Orientation of the Fock term on a complex, R-asymmetric bond: a
    complex transfer t on (R, a, b) = (R, 0, 1) (and its Hermitian partner
    at (-R, 1, 0)) with a CoulombInter (or Hund) row on the same (R, a, b),
    G from a converged UHFk run. The Fock correction of the green_init band
    sits on exactly the entries the transfer occupies -- the same bond, not
    the reversed one -- and equals UHFk's mean field of the same G. A real,
    R-symmetric G could not tell the two orientations apart."""

    CELL = (5, 4, 1)

    def _case(self, t, R, V, typ):
        import hwave.qlms as qlms
        nvol, nd = self.CELL[0] * self.CELL[1], 4
        Rm = (-R[0], -R[1])
        with tempfile.TemporaryDirectory() as d:
            _write_geometry(os.path.join(d, "geom.dat"), 2)
            _write_table(os.path.join(d, "transfer.dat"), "transfer", 2,
                         [(0, 0, 0, 0, 0.1), (0, 0, 1, 1, -0.1),
                          (R[0], R[1], 0, 1, t), (Rm[0], Rm[1], 1, 0, np.conj(t))])
            _write_table(os.path.join(d, "v.dat"), typ, 2,
                         [(R[0], R[1], 0, 1, V), (Rm[0], Rm[1], 1, 0, V)])
            inter = {"path_to_input": d, "Geometry": "geom.dat",
                     "Transfer": "transfer.dat", typ: "v.dat"}
            out = os.path.join(d, "out")
            qlms.run(input_dict={
                "log": {"print_level": 0, "print_step": 1000},
                "mode": {"mode": "UHFk", "flag_fock": True,
                         "param": {"Ncond": 20, "IterationMax": 3000,
                                   "EPS": 12, "Mix": 0.5, "RndSeed": 3,
                                   "T": 0.05, "CellShape": list(self.CELL),
                                   "SubShape": [1, 1, 1]}},
                "file": {"input": {"path_to_input": "", "interaction": inter},
                         "output": {"path_to_output": out, "green": "green"}}})
            fin = {"path_to_input": out, "green_init": "green.npz",
                   "interaction": inter}
            r = _rpa({"path_to_input": d, "interaction": inter}, self.CELL)
            G = r.read_init(fin)["green_init"]
            band = r._calc_trans_mod(G)
            bare = r._calc_trans_mod(np.zeros_like(G))
            u = _uhfk({"path_to_input": d, "interaction": inter}, self.CELL)
            u.Green = G.reshape(nvol, 2, 2, 2, 2).copy()
            u._make_ham()
            self.assertLess(_max_diff(band, u.ham), TOL)

        def to_r(hk):
            return (np.fft.fftn(np.asarray(hk).reshape(*self.CELL, nd, nd),
                                axes=(0, 1, 2)) / nvol).reshape(*self.CELL, nd, nd)
        t_r, f_r = to_r(bare), to_r(band - bare)
        support = lambda x: {tuple(int(i) for i in idx) for idx in
                             zip(*np.nonzero(np.abs(x) > 1e-8))
                             if tuple(idx[:3]) != (0, 0, 0)}
        t_sup, f_sup = support(t_r), support(f_r)
        self.assertTrue(f_sup)
        self.assertLessEqual(f_sup, t_sup)       # same bond entries as t
        # the declared (R, 0, 1) entry, spin up, carries both
        key = (R[0] % self.CELL[0], R[1] % self.CELL[1], 0, 0, 1)
        self.assertIn(key, t_sup)
        self.assertIn(key, f_sup)

    def test_coulombinter_r10(self):
        self._case(0.4 + 0.3j, (1, 0), 0.5, "CoulombInter")

    def test_coulombinter_r11(self):
        self._case(0.2 - 0.45j, (1, 1), 0.7, "CoulombInter")

    def test_hund_r10(self):
        self._case(0.4 + 0.3j, (1, 0), 0.5, "Hund")


class TestRPAOffsitePairHopVertexDrop(unittest.TestCase):
    """RPA's own vertex (not the band): PairHop rows that are off-site --
    declared R != 0, or R = 0 between orbitals of different sites per
    orbital_site -- are reported and dropped from the vertex (issue #157),
    judged on the DECLARED rows."""

    def _vertex(self, rows, site=None):
        import logging
        with tempfile.TemporaryDirectory() as d:
            fin = _model(d, {"PairHop": rows})
            records = []
            handler = logging.Handler()
            handler.emit = records.append
            lg = logging.getLogger("hwave.solver.rpa")
            lg.addHandler(handler)
            try:
                r = _rpa(fin, CELL, orbital_site=site)
            finally:
                lg.removeHandler(handler)
        dropped = any("PairHop declares" in rec.getMessage() for rec in records)
        return dropped, float(np.max(np.abs(r.ham_info.ham_inter_q)))

    def test_locality(self):
        onsite = [(0, 0, 0, 1, 0.3), (0, 0, 1, 0, 0.3)]
        # on-site: kept
        dropped, w = self._vertex(onsite)
        self.assertFalse(dropped)
        self.assertGreater(w, 1e-3)
        # same rows, orbitals on different sites: dropped
        dropped, w = self._vertex(onsite, site=(0, 1))
        self.assertTrue(dropped)
        self.assertEqual(w, 0.0)
        # declared R = (Nx, 0): off-site although it folds onto R = 0
        dropped, w = self._vertex([(CELL[0], 0, 0, 1, 0.3),
                                   (-CELL[0], 0, 1, 0, 0.3)])
        self.assertTrue(dropped)
        self.assertEqual(w, 0.0)


if __name__ == "__main__":
    unittest.main()
