"""The initial-state writer and its check (``hwave.initial.writer``; design
section 9, validation R0, R1, R2).

Every solver run here is the first iteration of ``hwave`` (``IterationMax =
1``) with its ``[file.output] green`` read back. With ``Mix = 0`` the saved
Green function is the initial one (``G = (1 - Mix) G_initial + Mix
G_new``), so write -> run -> read must give back the written seed exactly.

Green functions are never compared across numpy versions bit for bit (the
eigenbasis of a degenerate level is not fixed); the cross-version reference
compares member content and keys only.
"""
import contextlib
import json
import logging
import math
import os
import shutil
import tempfile
import time
import unittest
from fractions import Fraction

import numpy as np

from tests import initial_fixtures
from tests.initial_reference_capture import dump_toml

_REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), "initial_ref",
                    "writer_reference.npz")

_CHAIN = {((1, 0, 0), 0, 0): -1.0, ((-1, 0, 0), 0, 0): -1.0}


def _spec(entries, q, subshape, *, phase=0.0, sign=1):
    from hwave.initial.patterns import PatternSpec
    return PatternSpec(q=tuple(str(Fraction(x)) for x in q), subshape=subshape,
                       phase=phase, sign=sign, modulation="cos", entries=tuple(entries))


def _e(a=0, b=None, r=(0, 0, 0), **c):
    from hwave.initial.patterns import Entry
    return Entry(a, a if b is None else b, r, **c)


def _user(spec, label=None):
    from hwave.initial import Candidate
    return Candidate.user(pattern=spec, label=label)


def _neel():
    return _user(_spec([_e(cz=1.0)], ("1/2", "1/2", "0"), (2, 2, 1)), label="Neel")


def _read_toml(path):
    import tomli
    with open(path, "rb") as f:
        return tomli.load(f)


def _edit(toml, fn):
    doc = _read_toml(toml)
    fn(doc)
    with open(toml, "w") as f:
        f.write(dump_toml(doc))


def _set_run(toml, *, initial="initial.npz", mix=0.0, green="green", energy=None):
    """``IterationMax = 1``, ``Mix``, ``[file.input] initial`` and the
    ``[file.output] green`` (and ``energy``) files."""
    def fn(doc):
        p = doc["mode"]["param"]
        p["IterationMax"] = 1
        p["Mix"] = mix
        if initial is not None:
            doc["file"]["input"]["initial"] = initial
        out = doc["file"]["output"]
        out["green"] = green
        if energy is not None:
            out["energy"] = energy
    _edit(toml, fn)


@contextlib.contextmanager
def _cwd(path):
    old = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(old)


def _run_hwave(toml):
    """``hwave input.toml`` in the TOML's directory (in process)."""
    from hwave import qlms
    logger = logging.getLogger("qlms")
    level = logger.level
    with _cwd(os.path.dirname(os.path.abspath(toml))):
        try:
            qlms.run(input_file=os.path.basename(toml))
        finally:
            logger.setLevel(level)


def _load(path):
    with np.load(path, allow_pickle=False) as data:
        return {k: data[k] for k in data.files}


def _save(path, members):
    np.savez(path, **members)


def _prov(members, key):
    return json.loads(str(members[key][()]))


def _so_square(dest, L, subshape, *, U=2.0, T=0.1, ncond=None, norb_phys=1,
               zeeman=None):
    """Spin-orbital square lattice (interleaved spin-orbitals ``2a + s``),
    hopping 1, on-site U on each physical orbital; ``zeeman``: ``{a: b}``
    adds ``-b sigma^z`` on physical orbital a."""
    t = {}
    for a in range(norb_phys):
        for s in range(2):
            i = 2 * a + s
            for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
                t[(r, i, i)] = -1.0
    for a, b in (zeeman or {}).items():
        t[((0, 0, 0), 2 * a, 2 * a)] = -b
        t[((0, 0, 0), 2 * a + 1, 2 * a + 1)] = b
    p = {"Ncond": L * L * norb_phys if ncond is None else ncond, "T": T,
         "SubShape": list(subshape)}
    return initial_fixtures.write_input(
        dest, cellshape=(L, L, 1), norb=2 * norb_phys, transfer=t,
        interactions={"CoulombIntra": {((0, 0, 0), a, a): U for a in range(norb_phys)}},
        param=p, mode={"enable_spin_orbital": True})


def _sigma_y_chain(dest, *, U=1.0, T=0.0):
    """Design R1's sigma^y fixture: one orbital on a 2 x 1 chain, Sz free,
    Fock on, an on-site ``-1/4 sigma^z`` (9-column Transfer)."""
    entries = {((0, 0, 0), 0, 0, 0, 0): -0.25, ((0, 0, 0), 0, 0, 1, 1): 0.25}
    text = initial_fixtures.spin_transfer_text(1, entries)
    return initial_fixtures.write_input(
        dest, cellshape=(2, 1, 1), transfer=_CHAIN,
        interactions={"CoulombIntra": {((0, 0, 0), 0, 0): U}},
        param={"Ncond": 2, "T": T, "SubShape": [1, 1, 1]},
        extra_files={"transfer.dat": text})


def _positions(geom_centres, subshape):
    """Site positions of the folded cell, one site per orbital (the PR 1a
    default grouping): ``R_ir + tau_a`` for ``j = a + n * ir``."""
    B = subshape
    out = []
    for iz in range(B[2]):
        for iy in range(B[1]):
            for ix in range(B[0]):
                for c in geom_centres:
                    out.append([ix + c[0], iy + c[1], iz + c[2]])
    return np.array(out, dtype=float)


def _op_of(G, *, so, norb_phys, centres, subshape, cellshape):
    """PR 1a's order parameters of a stored (folded) G at r = 0."""
    from hwave.results import order_parameters_from_density, so_density_from_green
    subvol = int(np.prod(subshape))
    if so:
        rho = so_density_from_green(G[0, 0, :, 0, :], norb_phys, subvol)
    else:
        rho = np.einsum("sata->ast", G[0])
    return order_parameters_from_density(
        rho, site_of_orbital=np.arange(norb_phys * subvol),
        site_positions=_positions(centres, subshape), cell_shape=subshape,
        lattice_shape=cellshape)


def _q_index(op, q):
    from hwave.initial.cells import same_q
    for i, qq in enumerate(op.q_frac):
        if same_q(qq, tuple(Fraction(x) for x in q)):
            return i
    raise AssertionError("q {} not among {}".format(q, op.q_frac))


class _Tmp(unittest.TestCase):
    def setUp(self):
        self.root = tempfile.mkdtemp(prefix="initial_writer_")
        self._n = 0
        self._t0 = time.perf_counter()

    def tearDown(self):
        shutil.rmtree(self.root, ignore_errors=True)
        # each fast test method stays under 5 s (tests/heavy_tests.py)
        self.assertLess(time.perf_counter() - self._t0, 5.0)

    def dest(self):
        self._n += 1
        return os.path.join(self.root, "m{}".format(self._n))

    def write(self, toml, cand, name="initial.npz", **kw):
        from hwave.initial.writer import initial_state
        out = os.path.join(os.path.dirname(toml), name)
        path = initial_state(toml, cand, out=out, **kw)
        self.assertEqual(os.path.abspath(str(path)), os.path.abspath(out))
        return out

    def square(self, L=4, subshape=(2, 2, 1), **kw):
        kw.setdefault("U", 2.0)
        kw.setdefault("T", 0.1)
        return initial_fixtures.square(self.dest(), L, subshape=subshape, **kw)


# --------------------------------------------------------------------------
# R0: write -> run -> read gives back the seed
# --------------------------------------------------------------------------

class TestRoundTripLayouts(_Tmp):
    """``test_round_trip_layouts``, one method per layout."""

    def round_trip(self, toml, cand, *, q, so=False, norb_phys=1, centres=((0, 0, 0),),
                   expect_moment=None):
        from hwave.results import read_uhfk_green, read_uhfk_order_parameters
        from hwave.initial.writer import check_initial_state
        out = self.write(toml, cand)
        _set_run(toml)
        res = check_initial_state(toml)
        self.assertTrue(res.verified)
        doc = _read_toml(toml)
        p = doc["mode"]["param"]
        B = tuple(p["SubShape"])
        L = tuple(p["CellShape"])
        written = _load(out)
        _run_hwave(toml)
        green = read_uhfk_green(toml)
        np.testing.assert_allclose(green.green, written["green"], rtol=0, atol=1e-12)
        if int(np.prod(B)) > 1:
            self.assertEqual(res.key, "green_sublattice")
            G = written["green_sublattice"]
            np.testing.assert_allclose(green.sublattice.green, G, rtol=0, atol=1e-12)
        else:
            self.assertEqual(res.key, "green")
            self.assertNotIn("green_sublattice", written)
            G = written["green"]
        # PR 1a's Fourier components at the pattern's q, output vs written G
        got = read_uhfk_order_parameters(toml)
        want = _op_of(G, so=so, norb_phys=norb_phys, centres=centres, subshape=B,
                      cellshape=L)
        i = _q_index(got, q)
        j = _q_index(want, q)
        np.testing.assert_allclose(got.spin_fourier[i], want.spin_fourier[j], rtol=0,
                                   atol=1e-10)
        np.testing.assert_allclose(got.charge_fourier[i], want.charge_fourier[j], rtol=0,
                                   atol=1e-10)
        if expect_moment is not None:
            mu = "xyz".index(expect_moment)
            self.assertGreater(abs(want.spin_fourier[j][mu]), 0.05)
        return written

    def test_one_orbital_2x2(self):
        toml = self.square()
        self.round_trip(toml, _neel(), q=("1/2", "1/2", "0"), expect_moment="z")

    def test_two_orbital_2x1x1(self):
        t = {}
        for a in range(2):
            t[((1, 0, 0), a, a)] = -1.0
            t[((-1, 0, 0), a, a)] = -1.0
        t[((0, 0, 0), 0, 1)] = -0.3
        t[((0, 0, 0), 1, 0)] = -0.3
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), norb=2, transfer=t,
            centers=[(0, 0, 0), (0.5, 0, 0)],
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0, ((0, 0, 0), 1, 1): 2.0}},
            param={"Ncond": 8, "T": 0.1, "SubShape": [2, 1, 1]})
        cand = _user(_spec([_e(0, cz=1.0), _e(1, cz=-0.5)], ("1/2", "0", "0"), (2, 1, 1)))
        self.round_trip(toml, cand, q=("1/2", "0", "0"), norb_phys=2,
                        centres=((0, 0, 0), (0.5, 0, 0)), expect_moment="z")

    def test_spin_orbital_2x2(self):
        toml = _so_square(self.dest(), 4, (2, 2, 1))
        w = self.round_trip(toml, _neel(), q=("1/2", "1/2", "0"), so=True,
                            expect_moment="z")
        self.assertEqual(w["green_sublattice"].shape, (4, 1, 8, 1, 8))
        self.assertEqual(_prov(w, "index_convention"), "interleaved")
        self.assertIs(_prov(w, "spin_orbital"), True)

    def test_subvol_one(self):
        toml = self.square(L=4, subshape=(1, 1, 1))
        ferro = _user(_spec([_e(cz=1.0)], ("0", "0", "0"), (1, 1, 1)), label="ferro z")
        self.round_trip(toml, ferro, q=("0", "0", "0"), expect_moment="z")

    def test_4x3_incommensurate_q(self):
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 3, 1), transfer=initial_fixtures.square_transfer(),
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 2.0}},
            param={"Ncond": 12, "T": 0.1, "SubShape": [4, 3, 1]})
        cand = _user(_spec([_e(cz=1.0)], ("1/4", "1/3", "0"), (4, 3, 1)))
        self.round_trip(toml, cand, q=("1/4", "1/3", "0"), expect_moment="z")

    def test_sigma_y_pattern(self):
        toml = self.square()
        cand = _user(_spec([_e(cy=1.0)], ("1/2", "1/2", "0"), (2, 2, 1)))
        w = self.round_trip(toml, cand, q=("1/2", "1/2", "0"), expect_moment="y")
        # sigma^y lives in the imaginary spin-off-diagonal slots
        self.assertGreater(np.max(np.abs(w["green_sublattice"][:, 0, :, 1, :].imag)), 0.05)

    def test_bond_order_complex_amplitude(self):
        toml = initial_fixtures.write_input(
            self.dest(), cellshape=(4, 1, 1), transfer=_CHAIN,
            interactions={"CoulombIntra": {((0, 0, 0), 0, 0): 1.0},
                          "CoulombInter": {((1, 0, 0), 0, 0): 1.5, ((-1, 0, 0), 0, 0): 1.5}},
            param={"Ncond": 4, "T": 0.1, "SubShape": [2, 1, 1]})
        cand = _user(_spec([_e(r=(1, 0, 0), c0=complex(0.6, 0.8))], ("1/2", "0", "0"),
                           (2, 1, 1)), label="bond order")
        w = self.round_trip(toml, cand, q=("1/2", "0", "0"))
        G = w["green_sublattice"]
        # the two bonds of the cell differ (a q = 1/2 bond modulation) and
        # carry an imaginary part (complex A)
        b01 = G[0, 0, 0, 0, 1]
        b10 = G[1, 0, 1, 0, 0]
        self.assertGreater(abs(b01 - b10), 1e-3)
        self.assertGreater(abs(G[..., 0, :, 0, :].imag).max(), 1e-3)


class TestMixAndSpin(_Tmp):

    def test_mix_one_keeps_pattern(self):
        from hwave.results import read_uhfk_order_parameters
        toml = self.square(U=4.0)
        out = self.write(toml, _neel())
        G = _load(out)["green_sublattice"]
        _set_run(toml, mix=1.0)
        _run_hwave(toml)
        q = ("1/2", "1/2", "0")
        got = read_uhfk_order_parameters(toml)
        want = _op_of(G, so=False, norb_phys=1, centres=((0, 0, 0),), subshape=(2, 2, 1),
                      cellshape=(4, 4, 1))
        a = np.asarray(want.spin_fourier[_q_index(want, q)])
        b = np.asarray(got.spin_fourier[_q_index(got, q)])
        overlap = float(np.real(np.vdot(a, b))) / (np.linalg.norm(a) * np.linalg.norm(b))
        self.assertGreater(overlap, 0.9)

    def test_sx_sy_from_q0_seeds(self):
        """An x- and a y-ferro seed on the sigma^y fixture give a first
        iteration with Sx > 0 and Sy > 0 in energy.dat, each equal to
        ``N_site * n_cell * uniform`` of PR 1a on the same output."""
        from hwave.results import read_energy, read_uhfk_order_parameters
        for mu, key in ((0, "cx"), (1, "cy")):
            toml = _sigma_y_chain(self.dest())
            cand = _user(_spec([_e(**{key: 1.0})], ("0", "0", "0"), (1, 1, 1)))
            self.write(toml, cand)
            _set_run(toml, mix=1.0, energy="energy.dat")
            _run_hwave(toml)
            ene = read_energy(toml)
            op = read_uhfk_order_parameters(toml)
            s = (ene.sx, ene.sy)[mu]
            self.assertGreater(s, 0.05)
            total = op.spin.shape[0] * op.n_cell * np.asarray(op.uniform)
            self.assertAlmostEqual(total[mu], s, delta=1e-10)
            self.assertAlmostEqual(total[0], ene.sx, delta=1e-10)
            self.assertAlmostEqual(total[1], ene.sy, delta=1e-10)


# --------------------------------------------------------------------------
# the check
# --------------------------------------------------------------------------

class TestCheck(_Tmp):

    def written_neel(self):
        toml = self.square()
        out = self.write(toml, _neel())
        _set_run(toml)
        return toml, out

    def assert_refused(self, toml, npz, *words):
        from hwave.initial.writer import InitialStateError, check_initial_state
        with self.assertRaises(InitialStateError) as cm:
            check_initial_state(toml, npz)
        for w in words:
            self.assertIn(w, str(cm.exception))
        return str(cm.exception)

    def test_check_refusals(self):
        from hwave.initial.writer import check_initial_state
        toml, out = self.written_neel()
        self.assertTrue(check_initial_state(toml, out).verified)
        members = _load(out)
        # missing key the solver reads: it would fall back to another member
        bad = os.path.join(self.root, "missing.npz")
        _save(bad, {k: v for k, v in members.items() if k != "green_sublattice"})
        self.assert_refused(toml, bad, "green_sublattice", "missing")
        bad2 = os.path.join(self.root, "nogreen.npz")
        _save(bad2, {k: v for k, v in members.items()
                     if k not in ("green_sublattice", "green")})
        self.assert_refused(toml, bad2, "green_sublattice", "zero")
        # wrong shape
        bad3 = os.path.join(self.root, "shape.npz")
        m = dict(members)
        m["green_sublattice"] = members["green_sublattice"][:, :, :2, :, :2]
        _save(bad3, m)
        self.assert_refused(toml, bad3, "shape", "green_sublattice")
        # SubShape 2x2 -> 4x1 (same array shape: only provenance tells)
        _edit(toml, lambda d: d["mode"]["param"].update(SubShape=[4, 1, 1]))
        self.assert_refused(toml, out, "subshape", "[2, 2, 1]", "[4, 1, 1]")
        _edit(toml, lambda d: d["mode"]["param"].update(SubShape=[2, 2, 1]))
        # Ncond and 2Sz
        _edit(toml, lambda d: d["mode"]["param"].update(Ncond=14))
        self.assert_refused(toml, out, "ncond", "16", "14")
        _edit(toml, lambda d: d["mode"]["param"].update(Ncond=16))
        _edit(toml, lambda d: d["mode"]["param"].update({"2Sz": 2}))
        self.assert_refused(toml, out, "two_sz")
        _edit(toml, lambda d: d["mode"]["param"].pop("2Sz"))
        self.assertTrue(check_initial_state(toml, out).verified)

    def test_check_refuses_spin_orbital_flag(self):
        toml = _so_square(self.dest(), 4, (2, 2, 1))
        out = self.write(toml, _neel())
        _set_run(toml)
        # the same files read as a two-orbital normal-mode input
        _edit(toml, lambda d: d["mode"].pop("enable_spin_orbital"))
        self.assert_refused(toml, out, "spin_orbital")

    def test_tampered_provenance_field_named(self):
        toml, out = self.written_neel()
        members = _load(out)
        for key, value in (("cellshape", [8, 4, 1]), ("norb_orig", 2),
                           ("index_convention", "interleaved")):
            m = dict(members)
            m[key] = np.array(json.dumps(value))
            bad = os.path.join(self.root, key + ".npz")
            _save(bad, m)
            self.assert_refused(toml, bad, key)
        m = dict(members)
        m["initial_state_schema"] = np.array("2")
        bad = os.path.join(self.root, "schema.npz")
        _save(bad, m)
        self.assert_refused(toml, bad, "initial_state_schema")
        m = dict(members)
        del m["ncond"]
        bad = os.path.join(self.root, "noncond.npz")
        _save(bad, m)
        self.assert_refused(toml, bad, "ncond", "missing")

    def test_default_path_is_the_solver_path(self):
        from hwave.initial.writer import check_initial_state
        toml, out = self.written_neel()
        self.assertTrue(check_initial_state(toml).verified)
        os.remove(out)
        msg = self.assert_refused(toml, None, "initial.npz", "path_to_input")
        self.assertIn(os.path.abspath(out), msg)
        _edit(toml, lambda d: d["file"]["input"].pop("initial"))
        self.assert_refused(toml, None, "[file.input] initial")
        _edit(toml, lambda d: d["file"]["input"].update(initial="initial.npz",
                                                        initial_uhf="x.dat"))
        self.assert_refused(toml, None, "initial_uhf")

    def test_no_provenance_unverified(self):
        from hwave.initial.writer import (
            NO_PROVENANCE_WARNING, ProvenanceMissingError, check_initial_state)
        toml = self.square()
        _set_run(toml, initial=None, mix=1.0)
        _edit(toml, lambda d: d["file"]["input"].update(initial_mode="random"))
        _run_hwave(toml)
        plain = os.path.join(os.path.dirname(toml), "output", "green.npz")
        res = check_initial_state(toml, plain)
        self.assertFalse(res.verified)
        self.assertEqual(res.key, "green_sublattice")
        self.assertIn(NO_PROVENANCE_WARNING, res.warnings)
        self.assertEqual(
            NO_PROVENANCE_WARNING,
            "no provenance: SubShape, SO convention and electron sectors cannot be "
            "confirmed beyond shape")
        with self.assertRaises(ProvenanceMissingError) as cm:
            check_initial_state(toml, plain, require_provenance=True)
        self.assertIn("no provenance", str(cm.exception))

    def test_no_provenance_still_checked(self):
        """Design section 9.3: provenance-free files are checked for shape,
        Hermiticity, the eigenvalue bounds and the trace."""
        from hwave.initial.writer import check_initial_state
        toml, out = self.written_neel()
        members = _load(out)
        G = members["green_sublattice"]
        base = {"green_sublattice": G, "green": members["green"]}
        ok = os.path.join(self.root, "ok.npz")
        _save(ok, base)
        self.assertFalse(check_initial_state(toml, ok).verified)

        def refused(name, g, *words):
            path = os.path.join(self.root, name + ".npz")
            m = dict(base)
            m["green_sublattice"] = g
            _save(path, m)
            self.assert_refused(toml, path, *words)

        g = G.copy()
        g[1, 0, 0, 0, 1] += 0.1                      # breaks G(r) = G(-r)^dagger
        refused("herm", g, "Hermitian")
        g = G.copy()
        g[0, 0, 0, 0, 0] += 0.9                      # trace kept, eigenvalue > 1
        g[0, 0, 1, 0, 1] -= 0.9
        refused("eig", g, "eigenvalue")
        refused("trace", 0.9 * G, "trace", "Ncond")
        refused("shape_norb", G[:, :, :2, :, :2], "shape")
        refused("shape_so", G.reshape(4, 1, 8, 1, 8), "shape")
        refused("shape_sub", G[:2], "shape")
        # subvol = 1: a lone 'green' of another lattice
        toml1 = self.square(L=4, subshape=(1, 1, 1))
        p = os.path.join(self.root, "g1.npz")
        _save(p, {"green": G})
        self.assert_refused(toml1, p, "shape", "green")

    def test_object_member_refused(self):
        toml, out = self.written_neel()
        members = _load(out)
        members["extra"] = np.array([{"a": 1}], dtype=object)
        bad = os.path.join(self.root, "obj.npz")
        np.savez(bad, **members)
        self.assert_refused(toml, bad, "extra", "object")

    def test_sector_and_block_checks(self):
        """2Sz mode: the spin sectors must carry (Ncond +- 2Sz)/2 and a
        spin-off-diagonal entry lies outside the solver blocks."""
        toml = self.square(two_sz=0)
        out = self.write(toml, _neel())
        _set_run(toml)
        members = _load(out)
        G = members["green_sublattice"]

        def refused(name, g, *words):
            path = os.path.join(self.root, name + ".npz")
            m = dict(members)
            m["green_sublattice"] = g
            _save(path, m)
            self.assert_refused(toml, path, *words)

        X = np.zeros_like(G)                          # every up level filled (16, 0)
        X[0, 0, :, 0, :] = np.eye(G.shape[2])
        refused("sectors", 0.99 * G + 0.01 * X, "2Sz sectors")   # 16 electrons, 8.08 up
        Y = np.zeros_like(G)                          # one electron along +x per site
        for a in range(G.shape[2]):
            Y[0, :, a, :, a] = 0.5
        refused("blocks", 0.99 * G + 0.01 * Y, "block check", "transverse")

    def test_changed_input_warns(self):
        from hwave.initial.writer import check_initial_state
        toml, out = self.written_neel()
        with open(os.path.join(os.path.dirname(toml), "transfer.dat"), "a") as f:
            f.write("\n")
        res = check_initial_state(toml, out)
        self.assertTrue(res.verified)
        self.assertTrue(any("input_sha256" in w and "Transfer" in w for w in res.warnings),
                        res.warnings)


# --------------------------------------------------------------------------
# the file contract
# --------------------------------------------------------------------------

class TestFileContract(_Tmp):

    def test_keys_and_provenance(self):
        import hwave
        from hwave.initial.writer import PROVENANCE_KEYS
        toml = self.square()
        out = self.write(toml, _neel())
        m = _load(out)
        self.assertEqual(
            set(m), {"green", "green_sublattice", "green_convention",
                     "momentum_convention"} | set(PROVENANCE_KEYS))
        self.assertEqual(str(m["green_convention"][()]), "green_slot_first")
        self.assertEqual(str(m["momentum_convention"][()]), "e_plus_ikR")
        self.assertEqual(_prov(m, "initial_state_schema"), 1)
        self.assertEqual(_prov(m, "run_id"), _neel().candidate_id + "-b2x2x1")
        self.assertEqual(_prov(m, "subshape"), [2, 2, 1])
        self.assertEqual(_prov(m, "cellshape"), [4, 4, 1])
        self.assertEqual(_prov(m, "q"), ["1/2", "1/2", "0"])
        self.assertEqual(_prov(m, "ncond"), 16)
        self.assertIsNone(_prov(m, "two_sz"))
        self.assertIsNone(_prov(m, "seed"))
        self.assertEqual(_prov(m, "seed_base"), "bare")
        self.assertEqual(_prov(m, "index_convention"), "spin_block")
        self.assertEqual(_prov(m, "hwave_version"), hwave.__version__)
        self.assertGreater(_prov(m, "h_used"), 0)
        self.assertEqual(sorted(_prov(m, "input_sha256")),
                         ["CoulombIntra", "Geometry", "Transfer"])
        # canonical JSON: sorted keys, no whitespace
        text = str(m["input_sha256"][()])
        self.assertEqual(text, json.dumps(json.loads(text), sort_keys=True,
                                          separators=(",", ":")))
        # the unfolded view is the solver's deflate of the folded array
        from hwave.solver.fold import unfold_layout
        np.testing.assert_array_equal(
            m["green"], unfold_layout(m["green_sublattice"], cellshape=(4, 4, 1),
                                      subshape=(2, 2, 1), norb_orig=1, ns=2))

    def test_provenance_pickle_free(self):
        from hwave.initial.writer import PROVENANCE_KEYS
        toml = self.square()
        out = self.write(toml, _neel())
        with np.load(out, allow_pickle=False) as data:
            for k in data.files:
                v = data[k]                              # loads without pickle
                self.assertNotEqual(v.dtype.kind, "O", k)
                if k in PROVENANCE_KEYS:
                    self.assertEqual(v.shape, ())
                    self.assertEqual(v.dtype.kind, "U", k)
                    json.loads(str(v[()]))

    def test_bytes_identical_two_runs(self):
        toml = self.square()
        a = self.write(toml, _neel(), name="a.npz")
        b = self.write(toml, _neel(), name="b.npz")
        with open(a, "rb") as fa, open(b, "rb") as fb:
            self.assertEqual(fa.read(), fb.read())
        with np.load(a, allow_pickle=False) as data:
            self.assertEqual(list(data.files), sorted(data.files))
        umask = os.umask(0)
        os.umask(umask)
        self.assertEqual(os.stat(a).st_mode & 0o777, 0o666 & ~umask)

    def test_content_identical_across_numpy(self):
        """Design R1 (plan call): the key list and the provenance members
        exactly, the arrays to 1e-12, against the committed reference. The
        reference does not fix float bytes (eigh and FFT are not bitwise
        across numpy and BLAS builds): ``h_used`` is a result of the
        field search and is compared to 1e-12 relative, ``hwave_version``
        is the running version."""
        import hwave
        from tests.initial_reference_capture import writer_reference_case
        toml, cand = writer_reference_case(self.dest())
        out = self.write(toml, cand)
        got = _load(out)
        ref = _load(_REF)
        self.assertEqual(sorted(got), sorted(ref))
        for k in sorted(ref):
            if ref[k].dtype.kind == "U":
                if k == "hwave_version":
                    self.assertEqual(_prov(got, k), hwave.__version__)
                elif k == "h_used":
                    self.assertAlmostEqual(_prov(got, k) / _prov(ref, k), 1.0, delta=1e-12)
                else:
                    self.assertEqual(str(got[k][()]), str(ref[k][()]), k)
            else:
                np.testing.assert_allclose(got[k], ref[k], rtol=0, atol=1e-12, err_msg=k)

    def test_random_same_seed_same_npz(self):
        from hwave.initial import Candidate
        from hwave.initial.patterns import PatternSpec
        toml = self.square()

        def rnd(seed):
            spec = PatternSpec(q=None, subshape=(2, 2, 1), phase=None, sign=None,
                               modulation="random", seed=seed)
            return Candidate(candidate_id="random-s{}".format(seed), family="random",
                             sources=("random",), pattern=spec, label="random")

        a = _load(self.write(toml, rnd(7), name="a.npz"))
        b = _load(self.write(toml, rnd(7), name="b.npz"))
        c = _load(self.write(toml, rnd(8), name="c.npz"))
        np.testing.assert_array_equal(a["green_sublattice"], b["green_sublattice"])
        self.assertEqual(_prov(a, "seed"), 7)
        self.assertIsNone(_prov(a, "q"))
        self.assertGreater(np.max(np.abs(a["green_sublattice"] - c["green_sublattice"])),
                           1e-3)

    def test_conflicting_run_subshape(self):
        from hwave.initial import CandidateRun
        from hwave.initial.writer import initial_state
        toml = self.square()
        run = CandidateRun(_neel(), (2, 2, 1))
        out = os.path.join(self.root, "x.npz")
        with self.assertRaises(ValueError) as cm:
            initial_state(toml, run, out=out, subshape=(4, 4, 1))
        self.assertIn("subshape", str(cm.exception))
        self.assertFalse(os.path.exists(out))
        initial_state(toml, run, out=out, subshape=(2, 2, 1))
        self.assertEqual(_prov(_load(out), "run_id"), run.run_id)
        # a pattern on another SubShape than its own is a run of that B
        para = os.path.join(self.root, "para.npz")
        from hwave.initial import Candidate
        spec = PatternSpec_para()
        cand = Candidate(candidate_id="para", family="para", sources=("para",),
                         pattern=spec, label="para")
        initial_state(toml, cand, out=para, subshape=(2, 2, 1))
        m = _load(para)
        self.assertEqual(_prov(m, "run_id"), "para-b2x2x1")
        self.assertEqual(_prov(m, "h_used"), 0.0)

    def test_refused_rows_raise(self):
        from hwave.initial.writer import InitialStateError, initial_state
        # transverse spin in 2Sz mode: outside the solver blocks
        toml = self.square(two_sz=0)
        cand = _user(_spec([_e(cx=1.0)], ("1/2", "1/2", "0"), (2, 2, 1)))
        out = os.path.join(self.root, "x.npz")
        with self.assertRaises(InitialStateError) as cm:
            initial_state(toml, cand, out=out)
        self.assertIn("block", str(cm.exception))
        self.assertFalse(os.path.exists(out))
        # a pattern q the SubShape cannot hold
        toml = self.square()
        cand = _user(_spec([_e(cz=1.0)], ("1/4", "0", "0"), (4, 1, 1)))
        with self.assertRaises(ValueError) as cm:
            initial_state(toml, cand, out=out, subshape=(2, 2, 1))
        self.assertIn("not representable", str(cm.exception))
        # para has no amplitude to refuse; para is fine. UHFr TOML as model:
        uhfr = self.uhfr_input(4)
        with self.assertRaises(ValueError) as cm:
            initial_state(uhfr, _neel(), out=out)
        self.assertIn("site_map", str(cm.exception))

    def test_so_physical_index_realisation(self):
        """Review Focus 5: a cz entry on physical orbital 1 is +1/2 / -1/2 on
        the interleaved slots 2 and 3; the written seed has S^z > 0 on
        physical orbital 1 by PR 1a's so_density_from_green, and the solver
        reads it back so."""
        from hwave.initial import facts as facts_mod
        from hwave.initial.patterns import realise
        from hwave.results import (
            order_parameters_from_density, read_uhfk_order_parameters,
            so_density_from_green)
        toml = _so_square(self.dest(), 4, (1, 1, 1), norb_phys=2)
        spec = _spec([_e(1, cz=1.0)], ("0", "0", "0"), (1, 1, 1))
        f = facts_mod.model_facts(toml)
        O = realise(spec, f, (1, 1, 1), normalise=False)
        want = np.zeros((4, 4))
        want[2, 2], want[3, 3] = 0.5, -0.5
        for k in range(O.shape[0]):
            np.testing.assert_allclose(O[k], want, atol=1e-15)
        out = self.write(toml, _user(spec))
        G = _load(out)["green"]
        rho = so_density_from_green(G[0, 0, :, 0, :], 2, 1)
        op = order_parameters_from_density(rho, site_of_orbital=[0, 1], n_cell=None)
        self.assertGreater(op.spin[1][2], 0.05)
        self.assertAlmostEqual(op.spin[0][2], 0.0, delta=1e-10)
        self.assertAlmostEqual(op.spin[1][2], 0.5 * float(np.real(G[0, 0, 2, 0, 2]
                                                                   - G[0, 0, 3, 0, 3])),
                               delta=1e-14)
        _set_run(toml)
        _run_hwave(toml)
        got = read_uhfk_order_parameters(toml)
        np.testing.assert_allclose(got.spin, op.spin, atol=1e-12)

    # UHFr ------------------------------------------------------------------

    def uhfr_input(self, L, *, U=2.0, ncond=None, two_sz=None, initial="initial.def",
                   mix=0.0, index=None):
        """A UHFr Hubbard input on an ``L x L`` periodic square lattice with
        the bonds of :func:`initial_fixtures.square` (site ``i = index(x,
        y)``, default ``x + L*y``)."""
        index = index or (lambda x, y: x + L * y)
        d = self.dest()
        os.makedirs(d)
        N = L * L
        rows = []
        for y in range(L):
            for x in range(L):
                i = index(x, y)
                for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                    j = index((x + dx) % L, (y + dy) % L)
                    for s in range(2):
                        rows.append("{} {} {} {} 1.0 0.0".format(i, s, j, s))
        with open(os.path.join(d, "trans.def"), "w") as f:
            f.write("=====\nNTransfer {}\n=====\n=====\n=====\n".format(len(rows)))
            f.write("\n".join(rows) + "\n")
        with open(os.path.join(d, "coulombintra.def"), "w") as f:
            f.write("=====\nNCoulombIntra {}\n=====\n=====\n=====\n".format(N))
            f.write("".join("{} {}\n".format(i, U) for i in range(N)))
        p = {"Nsite": N, "Ncond": N if ncond is None else ncond, "IterationMax": 1,
             "EPS": 8, "Mix": mix, "RndSeed": 1, "T": 0.1}
        if two_sz is not None:
            p["2Sz"] = two_sz
        doc = {"log": {"print_level": 0},
               "mode": {"mode": "UHFr", "param": p},
               "file": {"input": {"path_to_input": "", "Initial": initial,
                                  "interaction": {"Trans": "trans.def",
                                                  "CoulombIntra": "coulombintra.def"}},
                        "output": {"path_to_output": "output", "initial": "initial_out.def"}}}
        toml = os.path.join(d, "input.toml")
        with open(toml, "w") as f:
            f.write(dump_toml(doc))
        return toml

    def test_uhfr_site_map(self):
        from hwave.initial.writer import InitialStateError, initial_state
        L = 4
        uhfk = self.square(L=L, subshape=(2, 2, 1))
        uhfr = self.uhfr_input(L)
        site_map = [(x, y, 0, 0) for y in range(L) for x in range(L)]
        out = os.path.join(os.path.dirname(uhfr), "initial.def")
        initial_state(uhfk, _neel(), out=out, site_map=site_map, uhfr_toml=uhfr)
        # the UHFk seed of the same run
        npz = self.write(uhfk, _neel())
        Gsub = _load(npz)["green_sublattice"]
        N = L * L

        def r_index(dx, dy):                 # supercell grid 2 x 2 x 1, z fastest
            return (dy % 2) + 2 * (dx % 2)

        want = np.zeros((2 * N, 2 * N), dtype=complex)
        for i, (xi, yi, _, _) in enumerate(site_map):
            for j, (xj, yj, _, _) in enumerate(site_map):
                ai = (xi % 2) + 2 * (yi % 2)
                aj = (xj % 2) + 2 * (yj % 2)
                r = r_index(xi // 2 - xj // 2, yi // 2 - yj // 2)
                for s in range(2):
                    for t in range(2):
                        want[i + s * N, j + t * N] = Gsub[r, s, ai, t, aj]
        _run_hwave(uhfr)
        back = self.read_uhfr_initial(os.path.join(os.path.dirname(uhfr), "output",
                                                   "initial_out.def"), N)
        np.testing.assert_allclose(back, want, rtol=0, atol=1e-12)
        mine = self.read_uhfr_initial(out, N)
        np.testing.assert_allclose(mine, want, rtol=0, atol=1e-12)
        # per-site S^z equals the UHFk seed's
        sz_r = 0.5 * np.real(np.diag(back)[:N] - np.diag(back)[N:])
        for i, (x, y, _, _) in enumerate(site_map):
            a = (x % 2) + 2 * (y % 2)
            sz_k = 0.5 * np.real(Gsub[0, 0, a, 0, a] - Gsub[0, 1, a, 1, a])
            self.assertAlmostEqual(sz_r[i], sz_k, delta=1e-12)
        self.assertGreater(abs(sz_r[0]), 0.05)
        self.assertAlmostEqual(sz_r[0], -sz_r[1], delta=1e-12)
        # a map that is not a bijection, and a mismatched UHFr input
        bad = list(site_map)
        bad[1] = bad[0]
        with self.assertRaises(ValueError) as cm:
            initial_state(uhfk, _neel(), out=out, site_map=bad, uhfr_toml=uhfr)
        self.assertIn("bijection", str(cm.exception))
        with self.assertRaises(ValueError):
            initial_state(uhfk, _neel(), out=out, site_map=site_map[:-1], uhfr_toml=uhfr)
        with self.assertRaises(ValueError) as cm:
            initial_state(uhfk, _neel(), out=out, site_map=site_map,
                          uhfr_toml=self.uhfr_input(L, two_sz=2))
        self.assertIn("2Sz", str(cm.exception))
        with self.assertRaises(ValueError) as cm:
            initial_state(uhfk, _neel(), out=out, site_map=site_map,
                          uhfr_toml=self.uhfr_input(L, ncond=14))
        self.assertIn("Ncond", str(cm.exception))
        with self.assertRaises(ValueError):
            initial_state(uhfk, _neel(), out=out, site_map=site_map)
        self.assertTrue(issubclass(InitialStateError, ValueError))

    def test_uhfr_first_iteration_matches_uhfk(self):
        """The orientation of the map, checked by physics: with the same
        seed (a complex in-plane Neel, so a transposed map would reverse
        S^y) and the same Hubbard model, the first UHFr iteration (``Mix =
        1``) is the first UHFk iteration mapped by :func:`uhfr_green`. A 6 x
        6 lattice with SubShape 2 x 2 gives a 3 x 3 supercell grid, where
        ``r`` and ``-r`` differ; the UHFr sites are numbered y fastest."""
        from hwave.initial import facts as facts_mod
        from hwave.initial.writer import initial_state, uhfr_green
        from hwave.results import read_uhfk_green
        L = 6
        uhfk = self.square(L=L, subshape=(2, 2, 1))
        index = lambda x, y: y + L * x
        uhfr = self.uhfr_input(L, mix=1.0, index=index)
        site_map = sorted(((x, y, 0, 0) for y in range(L) for x in range(L)),
                          key=lambda c: index(c[0], c[1]))
        cand = _user(_spec([_e(cy=1.0)], ("1/2", "1/2", "0"), (2, 2, 1)))
        out = os.path.join(os.path.dirname(uhfr), "initial.def")
        initial_state(uhfk, cand, out=out, site_map=site_map, uhfr_toml=uhfr)
        self.write(uhfk, cand)
        _set_run(uhfk, mix=1.0)
        _run_hwave(uhfk)
        _run_hwave(uhfr)
        N = L * L
        g_k = read_uhfk_green(uhfk).sublattice.green
        p = facts_mod.model_facts(uhfk).problem
        want = uhfr_green(g_k, p, site_map)
        got = self.read_uhfr_initial(os.path.join(os.path.dirname(uhfr), "output",
                                                  "initial_out.def"), N)
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-9)
        # the seed has an imaginary S^+ part, which the first iteration keeps
        self.assertGreater(np.max(np.abs(got[:N, N:].imag)), 0.05)
        # the written file read by UHFr's own S^y formula (Green[i][j] =
        # <c^dag_i c_j>, uhfr.py _calc_phys) equals PR 1a's S^y of the UHFk
        # seed site by site (a transposed map reverses it)
        from hwave.results import order_parameters_from_density
        Gsub = _load(os.path.join(os.path.dirname(uhfk), "initial.npz"))["green_sublattice"]
        seed_spin = order_parameters_from_density(
            np.einsum("sata->ast", Gsub[0]), site_of_orbital=np.arange(4), n_cell=None).spin
        mine = self.read_uhfr_initial(out, N)
        for i, (x, y, _, _) in enumerate(site_map):
            ud, du = mine[i, i + N], mine[i + N, i]
            sy = 0.5 * float(np.real(-1j * (ud - du)))
            ir = (x % 2) + 2 * (y % 2)
            self.assertAlmostEqual(sy, seed_spin[ir][1], delta=1e-12)
        self.assertGreater(abs(seed_spin[0][1]), 0.05)

    @staticmethod
    def read_uhfr_initial(path, N):
        g = np.zeros((2 * N, 2 * N), dtype=complex)
        with open(path) as f:
            lines = f.readlines()
        count = int(lines[1].split()[1])
        rows = [ln.split() for ln in lines[5:] if ln.strip()]
        assert len(rows) == count
        for i, s, j, t, re, im in rows:
            g[int(i) + int(s) * N, int(j) + int(t) * N] = complex(float(re), float(im))
        return g

    def test_so_to_uhfr_refused(self):
        from hwave.initial.writer import InitialStateError, initial_state
        toml = _so_square(self.dest(), 4, (2, 2, 1))
        uhfr = self.uhfr_input(4)
        out = os.path.join(self.root, "x.def")
        with self.assertRaises(InitialStateError) as cm:
            initial_state(toml, _neel(), out=out,
                          site_map=[(x, y, 0, 0) for y in range(4) for x in range(4)],
                          uhfr_toml=uhfr)
        self.assertEqual(str(cm.exception),
                         "UHFr has no spin-orbital mode; use UHFk with enable_spin_orbital")
        self.assertFalse(os.path.exists(out))


def PatternSpec_para():
    from hwave.initial.patterns import PatternSpec
    return PatternSpec(q=("0", "0", "0"), subshape=(1, 1, 1), phase=0.0, sign=1,
                       modulation="none")


# --------------------------------------------------------------------------
# where the solver looks
# --------------------------------------------------------------------------

class TestSolverPath(_Tmp):

    def test_initial_key_resolves_under_path_to_input(self):
        from hwave.initial.writer import (
            InitialStateError, check_initial_state, initial_key_value, initial_state,
            solver_initial_path)
        d = self.dest()
        inputs = os.path.join(d, "inputs")
        src = initial_fixtures.square(inputs, 4, U=2.0, T=0.1, subshape=(1, 1, 1))
        toml = os.path.join(d, "input.toml")
        shutil.move(src, toml)

        def paths(doc):
            doc["file"]["input"]["path_to_input"] = "inputs"
            doc["file"]["input"]["interaction"]["path_to_input"] = "inputs"
        _edit(toml, paths)
        ferro = _user(_spec([_e(cz=1.0)], ("0", "0", "0"), (1, 1, 1)))
        out = os.path.join(d, "initial.npz")
        initial_state(toml, ferro, out=out)
        key = initial_key_value(toml, out)
        self.assertEqual(key, "../initial.npz")
        _set_run(toml, initial=key)
        self.assertTrue(os.path.samefile(solver_initial_path(toml), out))
        self.assertTrue(check_initial_state(toml).verified)
        with self.assertLogs("qlms", level="INFO") as cm:
            _run_hwave(toml)
        self.assertTrue(any("read initial from inputs/../initial.npz" in m
                            for m in cm.output))
        self.assertTrue(os.path.exists(os.path.join(d, "output", "green.npz")))
        # the naive key resolves under path_to_input
        _set_run(toml, initial="initial.npz")
        with self.assertRaises(InitialStateError) as cm:
            check_initial_state(toml)
        msg = str(cm.exception)
        self.assertIn("inputs/initial.npz", msg)
        self.assertIn("path_to_input", msg)
        # without an initial key
        _edit(toml, lambda doc: doc["file"]["input"].pop("initial"))
        self.assertIsNone(solver_initial_path(toml))
        # an absolute key is used as it is (os.path.join)
        _set_run(toml, initial=out)
        self.assertEqual(str(solver_initial_path(toml)), out)
        self.assertEqual(initial_key_value(toml, out), "../initial.npz")


if __name__ == "__main__":
    unittest.main()
