"""Spin-polarized states through the spin-orbital path (issue #88, PR 0).

A Zeeman field ``h`` (spin up shifted by ``+h``, spin down by ``-h``) can
reach a spin-orbital calculation (``[mode] enable_spin_orbital = true``)
along two routes:

* **Route 1 ("transfer")** -- the field is written into the Transfer file
  as on-site entries on the interleaved spin-orbital indices
  ``2a+1`` (``+h``) and ``2a+2`` (``-h``);
* **Route 2 ("extern")** -- an ``Extern`` file in PHYSICAL orbital indices
  ``1..norb_phys``, scaled by ``[mode.param] coeff_extern``; the spin sign
  comes from ``sigma^z``.

What this module pins:

* :class:`TestRPARoutes` -- the two routes give the same RPA chi0q and chiq
  (gate 0-E1); a spin-orbital RPA chi0 equals the non-spin-orbital
  spin-diagonal one block by block (0-E2), and so does chi for an on-site
  ``CoulombInter`` model; an Extern index outside ``1..norb_phys`` (or an
  Extern file whose header declares another orbital count) is refused
  instead of being dropped without a word.
* :class:`TestFlexRoutes` -- the spin-orbital Hartree-Fock FLEX accepts
  Extern (0-F1): one map and the converged state agree with the Transfer
  route; the keys are case-insensitive through ``hwave.qlms.run``; the
  non-spin-orbital Hartree-Fock FLEX still refuses an external field.

Every run goes through the public entry :func:`hwave.qlms.run` (or, for the
Eliashberg chi0, :func:`hwave.sc._calc_chi0q_internal`), with fixtures
written to a temporary directory.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.wan90 as wan90

_HERE = os.path.dirname(os.path.abspath(__file__))
_SPIN_DIR = os.path.join(_HERE, "equivalence_input", "spin")

#: Zeeman field of every fixture
_H = 0.1


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


# ---------------------------------------------------------------- file I/O

def _read_w90_rows(path):
    """``(header_dim, nrpts, ndegen_lines, rows)`` of a sparse W90-like file;
    each row is ``(Rx, Ry, Rz, i, j, re, im)`` with 1-based ``i, j``."""
    with open(path) as f:
        lines = [s.strip() for s in f.readlines()]
    dim, nr = int(lines[1]), int(lines[2])
    nskip = -(-nr // 15)
    deg = lines[3:3 + nskip]
    rows = []
    for ln in lines[3 + nskip:]:
        if not ln:
            continue
        v = ln.split()
        rows.append(tuple(int(x) for x in v[:5]) + (float(v[5]), float(v[6])))
    return dim, nr, deg, rows


def _write_w90_file(path, header_dim, rows, head="generated"):
    """Write ``rows`` in the sparse W90-like format (degeneracy 1 for every
    distinct R vector, in order of first appearance)."""
    rvecs = []
    for r in rows:
        if r[:3] not in rvecs:
            rvecs.append(r[:3])
    nr = max(1, len(rvecs))
    with open(path, "w") as f:
        f.write("{}\n{}\n{}\n".format(head, header_dim, nr))
        for i in range(0, nr, 15):
            f.write(" ".join(["1"] * min(15, nr - i)) + "\n")
        for r in rows:
            f.write("  %d %d %d %d %d  %.12f  %.12f\n" % r)
    return path


def _write_w90(d, header_dim, rows, name="extern.dat"):
    """An Extern-style file in ``d``; returns its path."""
    return _write_w90_file(os.path.join(d, name), header_dim, rows, head="Extern")


def _write_geom(path, centers):
    with open(path, "w") as f:
        f.write("1.0 0.0 0.0\n0.0 1.0 0.0\n0.0 0.0 1.0\n%d\n" % len(centers))
        for c in centers:
            f.write("%.12f %.12f %.12f\n" % tuple(c))


def _read_geom_centers(path):
    with open(path) as f:
        lines = [s.strip() for s in f.readlines()]
    n = int(lines[3])
    return [tuple(float(x) for x in lines[4 + i].split()) for i in range(n)]


def _interleave_rows(rows):
    """Physical-orbital rows -> the spin-diagonal spin-orbital copy on the
    interleaved 1-based file index ``2a + s + 1``."""
    out = []
    for (rx, ry, rz, i, j, re, im) in rows:
        for s in (0, 1):
            out.append((rx, ry, rz, 2 * (i - 1) + s + 1, 2 * (j - 1) + s + 1, re, im))
    return out


def _zeeman_transfer_rows(norb_phys, h):
    """Route 1: the on-site +h (up) / -h (down) rows of every orbital."""
    rows = []
    for a in range(norb_phys):
        rows.append((0, 0, 0, 2 * a + 1, 2 * a + 1, h, 0.0))
        rows.append((0, 0, 0, 2 * a + 2, 2 * a + 2, -h, 0.0))
    return rows


# ---------------------------------------------------------------- runners

def _run(d, mode, param, interaction, output, enable_spin_orbital, extra_mode=None):
    """One ``hwave.qlms.run`` call; returns the loaded output archives (and
    ``"status"``, the value ``run`` returns)."""
    import hwave.qlms as qlms
    out_dir = os.path.join(d, "output")
    mode_tbl = {"mode": mode, "enable_spin_orbital": enable_spin_orbital,
                "param": dict(param)}
    mode_tbl.update(extra_mode or {})
    inter = {"path_to_input": d}
    inter.update(interaction)
    inp = {"log": {"print_level": 0},
           "mode": mode_tbl,
           "file": {"input": {"path_to_input": d, "interaction": inter},
                    "output": dict(output, path_to_output=out_dir)}}
    with _quiet():
        status = qlms.run(input_dict=inp)
    res = {"status": status}
    for key, fname in output.items():
        with np.load(os.path.join(out_dir, fname)) as z:
            res[key] = np.array(z[key])
    return res


def _so_1orb_fixture(d, route, h=_H, extern_rows=None, extern_header_dim=1,
                     u=1.0, extern_key="Extern"):
    """The 0-E1 fixture: ``transfer_spinful.dat`` (one physical orbital,
    complex SOC 0.35+0.15i) on ``geom_so.dat``, CoulombIntra ``u``.

    ``route = "transfer"`` appends the +-h rows to the Transfer file (the
    file has no on-site diagonal rows to collide with); ``"extern"`` writes
    an Extern file in physical indices (``extern_rows`` default: ``H_11(0) =
    h``). Returns the interaction table of the run."""
    dim, _nr, _deg, rows = _read_w90_rows(os.path.join(_SPIN_DIR, "transfer_spinful.dat"))
    assert dim == 2 and not any(r[:3] == (0, 0, 0) and r[3] == r[4] for r in rows)
    if route == "transfer":
        rows = rows + _zeeman_transfer_rows(1, h)
    _write_w90_file(os.path.join(d, "transfer.dat"), dim, rows, head="Transfer")
    shutil.copy(os.path.join(_SPIN_DIR, "geom_so.dat"), os.path.join(d, "geom.dat"))
    _write_w90_file(os.path.join(d, "coulombintra.dat"), 1,
                    [(0, 0, 0, 1, 1, u, 0.0)], head="CoulombIntra")
    inter = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
             "CoulombIntra": "coulombintra.dat"}
    if route == "extern":
        if extern_rows is None:
            extern_rows = [(0, 0, 0, 1, 1, h, 0.0)]
        _write_w90(d, extern_header_dim, extern_rows)
        inter[extern_key] = "extern.dat"
    else:
        assert route == "transfer", route
    return inter


def _so_2orb_fixture(d, route, types=("CoulombIntra",)):
    """The 0-E2 fixture, in the spin-orbital encoding: a generated
    4-spin-orbital geometry (the centres of ``tests/equivalence_input/spin/
    geom.dat``, each twice) and an interleaved copy of its 2-orbital
    ``transfer.dat``; the field is ``extern_zeeman.dat`` (route "extern",
    ``H_aa(0) = 0.1``) or the equivalent Transfer rows. No SOC."""
    centers = _read_geom_centers(os.path.join(_SPIN_DIR, "geom.dat"))
    _write_geom(os.path.join(d, "geom.dat"), [c for c in centers for _s in (0, 1)])
    dim, _nr, _deg, rows = _read_w90_rows(os.path.join(_SPIN_DIR, "transfer.dat"))
    rows = _interleave_rows(rows)
    if route == "transfer":
        rows = rows + _zeeman_transfer_rows(dim, _H)
    _write_w90_file(os.path.join(d, "transfer.dat"), 2 * dim, rows, head="Transfer")
    inter = {"Geometry": "geom.dat", "Transfer": "transfer.dat"}
    if route == "extern":
        shutil.copy(os.path.join(_SPIN_DIR, "extern_zeeman.dat"), os.path.join(d, "extern.dat"))
        inter["Extern"] = "extern.dat"
    _add_2orb_interactions(d, inter, types)
    return inter


def _add_2orb_interactions(d, inter, types):
    src = {"CoulombIntra": "coulombintra.dat", "CoulombInter": "onsite_inter.dat"}
    for t in types:
        shutil.copy(os.path.join(_SPIN_DIR, src[t]), os.path.join(d, src[t]))
        inter[t] = src[t]


_E1_PARAM = {"T": 0.5, "filling": 0.4, "Nmat": 16, "coeff_extern": 1.0}
_E2_PARAM = {"T": 1.0, "filling": 0.4, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
             "Nmat": 16, "coeff_extern": 1.0}


class _TmpDirs:
    """Fresh temporary directories, removed at exit of the ``with`` block."""

    def __init__(self):
        self.dirs = []

    def new(self):
        d = tempfile.mkdtemp(prefix="so_zeeman_")
        self.dirs.append(d)
        return d

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        for d in self.dirs:
            shutil.rmtree(d, ignore_errors=True)


def _rpa_so(route, extern_rows=None, cell=(4, 1, 1), sub_shape=None,
            extern_header_dim=1, assert_no_sublattice=False):
    """Spin-orbital RPA (general, ring) on the 0-E1 fixture through
    ``hwave.qlms.run``; returns ``{"chi0q", "chiq"}``. ``sub_shape = None``
    omits ``SubShape`` (no sublattice)."""
    if assert_no_sublattice:
        assert sub_shape is None, "assert_no_sublattice needs sub_shape=None"
    param = dict(_E1_PARAM, CellShape=list(cell))
    if sub_shape is not None:
        param["SubShape"] = list(sub_shape)
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_1orb_fixture(d, route, extern_rows=extern_rows,
                                 extern_header_dim=extern_header_dim)
        return _run(d, "RPA", param, inter, {"chi0q": "chi0q.npz", "chiq": "chiq.npz"},
                    True, {"calc_scheme": "general", "calc_type": "ring"})


def _sc_chi0_general(route):
    """The spin-orbital chi0 of the Eliashberg front end on the 0-E2 fixture:
    ``hwave.sc._calc_chi0q_internal(chi0q_tensor="general")``, shape
    ``(nmat, nvol, nd, nd, nd, nd)`` in spin-block order (``nd = 4``)."""
    import hwave.sc as sc
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_2orb_fixture(d, route)
        inter["path_to_input"] = d
        inp = {"log": {"print_level": 0},
               "mode": {"mode": "RPA", "enable_spin_orbital": True, "param": dict(_E2_PARAM)},
               "file": {"input": {"path_to_input": d, "interaction": inter}}}
        with _quiet():
            return np.asarray(sc._calc_chi0q_internal(inp, chi0q_tensor="general"))


def _sc_chi_general(route, types=("CoulombIntra",)):
    """The spin-orbital RPA chiq (general, ring) on the 0-E2 fixture
    (``hwave.qlms.run``), shape ``(nmat, nvol, nd, nd, nd, nd)`` (the
    spinful archive has no block axis)."""
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_2orb_fixture(d, route, types=types)
        res = _run(d, "RPA", _E2_PARAM, inter, {"chiq": "chiq.npz"}, True,
                   {"calc_scheme": "general", "calc_type": "ring"})
    chiq = res["chiq"]
    assert chiq.ndim == 6, chiq.shape
    return chiq


def _rpa_non_so_spin_diag(types=("CoulombIntra",)):
    """The non-spin-orbital spin-diagonal RPA (general, ring) on
    ``tests/equivalence_input/spin`` (2 orbitals + ``extern_zeeman.dat``):
    ``{"chi0q": (2, nmat, nvol, n, n, n, n), "chiq": (nmat, nvol, 2n, 2n, 2n,
    2n)}`` (chiq couples the spins, in spin-block order)."""
    with _TmpDirs() as t:
        d = t.new()
        for f in ("geom.dat", "transfer.dat", "extern_zeeman.dat"):
            shutil.copy(os.path.join(_SPIN_DIR, f), os.path.join(d, f))
        inter = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "Extern": "extern_zeeman.dat"}
        _add_2orb_interactions(d, inter, types)
        return _run(d, "RPA", _E2_PARAM, inter, {"chi0q": "chi0q.npz", "chiq": "chiq.npz"},
                    False, {"calc_scheme": "general", "calc_type": "ring"})


def _ssss_slots(chi, s, n=2):
    """The all-spin-``s`` block of a spin-block 4-index tensor
    ``(..., nd, nd, nd, nd)`` with ``nd = 2 n``."""
    sl = slice(s * n, (s + 1) * n)
    return chi[..., sl, sl, sl, sl]


def _cross_spin_slots(chi, n=2):
    """Every entry of a spin-block 4-index tensor whose four spins are not
    all equal (flattened)."""
    nd = 2 * n
    spin = np.arange(nd) // n
    s1, s2, s3, s4 = np.meshgrid(spin, spin, spin, spin, indexing="ij")
    mask = ~((s1 == s2) & (s2 == s3) & (s3 == s4))
    return chi[..., mask]


_F1_PARAM = {"T": 0.5, "filling": 0.4, "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
             "Nmat": 16, "Mix": 0.5, "mixing_scheme": "linear",
             "flex_hartree_fock": True}


def _flex_so_one_map(route, key="Extern", coeff_key="coeff_extern", iteration_max=1, eps=12):
    """Spin-orbital Hartree-Fock FLEX on the 0-E1 fixture through
    ``hwave.qlms.run``: ``{"sigma", "green", "chiq", "scf_tol"}`` after
    ``iteration_max`` maps. ``key`` / ``coeff_key`` spell the Extern file key
    and ``coeff_extern`` (case-insensitivity)."""
    param = dict(_F1_PARAM, IterationMax=iteration_max, EPS=eps)
    param.pop("coeff_extern", None)
    param[coeff_key] = 1.0
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_1orb_fixture(d, route, extern_key=key)
        res = _run(d, "FLEX", param, inter,
                   {"sigma": "sigma.npz", "green": "green.npz", "chiq": "chiq.npz"},
                   True, {"calc_scheme": "general"})
    res["scf_tol"] = 10.0 ** (-eps)
    res["scf_converged"] = res["status"]["scf_converged"]
    return res


def _flex_so(route):
    """The converged spin-orbital Hartree-Fock FLEX of the 0-E1 fixture."""
    return _flex_so_one_map(route, iteration_max=200, eps=10)


def _flex_non_so_hf_with_extern():
    """Non-spin-orbital Hartree-Fock FLEX with an Extern file (refused)."""
    with _TmpDirs() as t:
        d = t.new()
        for f in ("geom.dat", "transfer.dat", "extern_zeeman.dat", "coulombintra.dat"):
            shutil.copy(os.path.join(_SPIN_DIR, f), os.path.join(d, f))
        inter = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "Extern": "extern_zeeman.dat", "CoulombIntra": "coulombintra.dat"}
        param = dict(_E2_PARAM, IterationMax=1, flex_hartree_fock=True)
        _run(d, "FLEX", param, inter, {"green": "green.npz"}, False,
             {"calc_scheme": "general"})


# ---------------------------------------------------------------- tests

class TestRPARoutes(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="so_zeeman_")
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_extern_equals_transfer_zeeman(self):                       # 0-E1
        a, b = _rpa_so(route="extern"), _rpa_so(route="transfer")
        for key in ("chi0q", "chiq"):
            with self.subTest(key=key):
                self.assertGreater(np.abs(b[key]).max(), 1e-3)          # anti-vacuity
                np.testing.assert_allclose(a[key], b[key], rtol=0,
                                           atol=1e-12 * np.abs(b[key]).max())
        # the field is load-bearing: without it chi0 differs
        c = _rpa_so(route="extern", extern_rows=[(0, 0, 0, 1, 1, 0.0, 0.0)])
        self.assertGreater(np.abs(c["chi0q"] - b["chi0q"]).max(), 1e-4)

    def test_spinful_chi0_matches_non_so_spin_diag(self):              # 0-E2
        chi0_so = _sc_chi0_general(route="extern")
        chi0_ns = _rpa_non_so_spin_diag()["chi0q"]
        self.assertEqual(chi0_so.shape[-4:], (4, 4, 4, 4))
        self.assertEqual(chi0_ns.shape[0], 2)
        # the two spin blocks differ (the field polarizes)
        self.assertGreater(np.abs(chi0_ns[0] - chi0_ns[1]).max(), 1e-4)
        for s in (0, 1):
            with self.subTest(spin=s):
                np.testing.assert_allclose(_ssss_slots(chi0_so, s), chi0_ns[s], rtol=0,
                                           atol=1e-12 * np.abs(chi0_ns).max())
        self.assertGreater(np.abs(_cross_spin_slots(chi0_so)).max(), 1e-3)

    def test_coulombinter_chi_matches_non_so(self):                     # 0-E2 (chi, CoulombInter)
        """Without SOC the spin-orbital RPA resolves to the spin-diagonal
        layout, so both archives hold the same ``(nmat, nvol, nd, nd, nd,
        nd)`` spin-block chi (``nd = 4``); they agree entry by entry, the
        all-up block included. Compared for ``CoulombInter`` only: the
        non-spin-orbital spin-diagonal solve has no on-site exchange
        crossing for the other types."""
        chi_so = _sc_chi_general(route="extern", types=("CoulombInter",))
        chi_ns = _rpa_non_so_spin_diag(types=("CoulombInter",))["chiq"]
        self.assertEqual(chi_so.shape, chi_ns.shape)
        self.assertGreater(np.abs(_ssss_slots(chi_ns, 0)).max(), 1e-3)
        np.testing.assert_allclose(_ssss_slots(chi_so, 0), _ssss_slots(chi_ns, 0),
                                   rtol=0, atol=1e-12)
        np.testing.assert_allclose(chi_so, chi_ns, rtol=0, atol=1e-12)

    def test_extern_index_out_of_range_no_sublattice(self):            # silently dropped before
        with self.assertRaisesRegex(ValueError, r"Extern file '.*extern\.dat': orbital index 2 "
                                    r"is outside 1\.\.1\."):
            _rpa_so(route="extern", extern_rows=[(0, 0, 0, 2, 2, 0.1, 0.0)], sub_shape=None,
                    assert_no_sublattice=True)

    def test_extern_index_out_of_range_sublattice(self):
        with self.assertRaisesRegex(ValueError, r"outside 1\.\.1\. In spin-orbital mode Extern "
                                    r"uses PHYSICAL"):
            _rpa_so(route="extern", extern_rows=[(0, 0, 0, 2, 2, 0.1, 0.0)], cell=(4, 1, 1),
                    sub_shape=(2, 1, 1))

    def test_extern_header_dimension_refused(self):                    # rows all in range
        for sub in (None, (2, 1, 1)):
            with self.subTest(sublattice=sub), \
                    self.assertRaisesRegex(ValueError, r"orbital index 2 is outside 1\.\.1"):
                _rpa_so(route="extern", extern_header_dim=2,
                        extern_rows=[(0, 0, 0, 1, 1, 0.1, 0.0)], cell=(4, 1, 1), sub_shape=sub)

    def test_header_dim_reader_contract(self):
        path = _write_w90(self.d, header_dim=2, rows=[(0, 0, 0, 1, 1, 0.1, 0.0)])
        self.assertEqual(wan90.read_w90_header_dim(path), 2)
        self.assertEqual(wan90.read_w90(path), {((0, 0, 0), (0, 0)): 0.1 + 0j})


class TestExternIndexChecker(unittest.TestCase):
    """:func:`hwave.solver.spin_orbital.check_extern_indices` directly."""

    def test_message_text(self):
        from hwave.solver.spin_orbital import check_extern_indices
        with self.assertRaises(ValueError) as cm:
            check_extern_indices({((0, 0, 0), (0, 2)): 0.1}, 2, "x/ext.dat")
        self.assertEqual(str(cm.exception), (
            "Extern file 'x/ext.dat': orbital index 3 is outside 1..2. In spin-orbital "
            "mode Extern uses PHYSICAL orbital indices (norb_phys = 2), while the Transfer "
            "file uses the 2·norb_phys = 4 spin-orbital indices; write the field per "
            "physical orbital (the sign per spin comes from σ^z)."))

    def test_in_range_and_header(self):
        from hwave.solver.spin_orbital import check_extern_indices
        tbl = {((0, 0, 0), (0, 1)): 0.1, ((1, 0, 0), (1, 0)): 0.2}
        check_extern_indices(tbl, 2, "e.dat")                      # no header: fine
        check_extern_indices(tbl, 2, "e.dat", declared_norb=2)     # matching header
        with self.assertRaisesRegex(ValueError, r"orbital index 4 is outside 1\.\.2\."):
            check_extern_indices(tbl, 2, "e.dat", declared_norb=4)
        with self.assertRaisesRegex(ValueError, r"orbital index 0 is outside 1\.\.2\."):
            check_extern_indices({((0, 0, 0), (-1, 0)): 0.1}, 2, "e.dat")

    def test_header_dim_reader_refuses_malformed(self):
        d = tempfile.mkdtemp(prefix="so_zeeman_")
        self.addCleanup(shutil.rmtree, d, True)
        for line2 in ("", "abc", "0", "-2", "2.0"):
            with self.subTest(line2=line2):
                p = os.path.join(d, "bad.dat")
                with open(p, "w") as f:
                    f.write("comment\n{}\n1\n1\n 0 0 0 1 1 0.1 0.0\n".format(line2))
                with self.assertRaisesRegex(ValueError, r"bad\.dat: line 2 must be the number "
                                            r"of orbitals \(a positive integer\), got "):
                    wan90.read_w90_header_dim(p)
        p = os.path.join(d, "short.dat")
        with open(p, "w") as f:
            f.write("only a comment\n")
        with self.assertRaisesRegex(ValueError, r"short\.dat: line 2 must be"):
            wan90.read_w90_header_dim(p)


class TestFlexRoutes(unittest.TestCase):
    def test_extern_equals_transfer_zeeman_one_map(self):              # 0-F1
        a, b = _flex_so_one_map("extern"), _flex_so_one_map("transfer")
        for key in ("sigma", "green", "chiq"):
            with self.subTest(key=key):
                self.assertGreater(np.abs(b[key]).max(), 1e-3)
                np.testing.assert_allclose(a[key], b[key], rtol=0,
                                           atol=1e-12 * np.abs(b[key]).max())

    def test_converged_agree(self):
        a, b = _flex_so("extern"), _flex_so("transfer")
        self.assertTrue(a["scf_converged"] and b["scf_converged"])
        np.testing.assert_allclose(a["green"], b["green"], rtol=0, atol=10 * a["scf_tol"])

    def test_mixed_case_keys(self):
        np.testing.assert_array_equal(
            _flex_so_one_map("extern", key="EXTERN", coeff_key="Coeff_Extern")["green"],
            _flex_so_one_map("extern")["green"])

    def test_non_so_hf_flex_still_refuses_extern(self):
        with self.assertRaisesRegex(ValueError, "does not support an external field"):
            _flex_non_so_hf_with_extern()


if __name__ == "__main__":
    unittest.main()
