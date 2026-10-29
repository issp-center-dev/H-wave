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
* :class:`TestScBand` -- the spin-orbital Eliashberg front end (``hwave_sc``)
  accepts Extern: its epsilon(k) equals the spinful RPA H0(k) with Extern;
  the two routes add (with a WARNING); a zero ``coeff_extern`` gives the
  Transfer band; the keys are case-insensitive; a bad Extern index or
  header is refused in calc, load and flex modes before any heavy work; the
  dynamic solver issues no false "Extern is ignored" WARNING and agrees
  between the routes.
* :class:`TestScSpectrum` -- a Zeeman field without spin-orbit coupling
  conserves S_z: the static spectrum is the union of the three pair-S_z
  sector spectra (gate 0-S1), and the two routes give the same spectrum.

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
from unittest import mock

import numpy as np

import hwave.qlmsio.wan90 as wan90
import tests.test_sc_spin_orbital as S        # the shared two-orbital model

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
            extern_header_dim=1, assert_no_sublattice=False, extern_key="Extern"):
    """Spin-orbital RPA (general, ring) on the 0-E1 fixture through
    ``hwave.qlms.run``; returns ``{"chi0q", "chiq"}``. ``sub_shape = None``
    omits ``SubShape`` (no sublattice); ``extern_header_dim`` is written
    verbatim as line 2 of the Extern file (a string writes a malformed
    header); ``extern_key`` spells the ``[file.input.interaction]`` key."""
    if assert_no_sublattice:
        assert sub_shape is None, "assert_no_sublattice needs sub_shape=None"
    param = dict(_E1_PARAM, CellShape=list(cell))
    if sub_shape is not None:
        param["SubShape"] = list(sub_shape)
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_1orb_fixture(d, route, extern_rows=extern_rows,
                                 extern_header_dim=extern_header_dim,
                                 extern_key=extern_key)
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


def _sc_chi_general(route, types=("CoulombIntra",), spinful=True):
    """The spin-orbital RPA chiq (general, ring) on the 0-E2 fixture
    (``hwave.qlms.run``), shape ``(nmat, nvol, nd, nd, nd, nd)`` in
    spin-block order (``nd = 4``).

    ``spinful = True`` forces the spinful path (``RPA._keep_spinful``, as the
    Eliashberg front end does, ``hwave.sc._calc_chi0q_internal``): without
    SOC the RPA would otherwise resolve the band to ``spin-diag`` and solve
    exactly the non-spin-orbital problem. ``spinful = False`` keeps that
    default resolution (a check of the spin-orbital ENCODING only)."""
    import hwave.solver.rpa as rpa_mod
    with _TmpDirs() as t:
        d = t.new()
        inter = _so_2orb_fixture(d, route, types=types)
        with mock.patch.object(rpa_mod.RPA, "_keep_spinful", spinful, create=True):
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
        """The SPINFUL spin-orbital chi (path forced, as hwave_sc does) per
        all-``s`` spin block equals the non-spin-orbital spin-diagonal chi.
        Compared for ``CoulombInter`` only: the non-spin-orbital
        spin-diagonal solve has no on-site exchange crossing for the other
        types (the spinful antisymmetrized vertex does)."""
        chi_so = _sc_chi_general(route="extern", types=("CoulombInter",))
        chi_ns = _rpa_non_so_spin_diag(types=("CoulombInter",))["chiq"]
        self.assertEqual(chi_so.shape, chi_ns.shape)
        # the spinful path was taken: the transverse (spin-flip) slots,
        # which the spin-diagonal layout leaves at zero, are populated
        n = 2
        sp = np.arange(2 * n) // n
        s1, s2, s3, s4 = np.meshgrid(sp, sp, sp, sp, indexing="ij")
        transverse = ~((s1 == s2) & (s3 == s4))
        self.assertEqual(np.abs(chi_ns[..., transverse]).max(), 0.0)
        self.assertGreater(np.abs(chi_so[..., transverse]).max(), 1e-3)
        for s in (0, 1):
            with self.subTest(spin=s):
                self.assertGreater(np.abs(_ssss_slots(chi_ns, s)).max(), 1e-3)
                np.testing.assert_allclose(_ssss_slots(chi_so, s), _ssss_slots(chi_ns, s),
                                           rtol=0, atol=1e-12)

    def test_so_encoding_spin_diag_equals_non_so(self):
        """Encoding check: without SOC and without forcing the spinful path,
        the spin-orbital input resolves to the spin-diagonal solve and its
        chi equals the non-spin-orbital one entry by entry."""
        chi_so = _sc_chi_general(route="extern", types=("CoulombInter",), spinful=False)
        chi_ns = _rpa_non_so_spin_diag(types=("CoulombInter",))["chiq"]
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
        for sub, key in ((None, "Extern"), ((2, 1, 1), "Extern"), (None, "EXTERN")):
            with self.subTest(sublattice=sub, key=key), \
                    self.assertRaisesRegex(ValueError, r"Extern file '.*extern\.dat': its header "
                                           r"declares 2 orbitals, but norb_phys = 1\. In "
                                           r"spin-orbital mode Extern uses PHYSICAL"):
                _rpa_so(route="extern", extern_header_dim=2,
                        extern_rows=[(0, 0, 0, 1, 1, 0.1, 0.0)], cell=(4, 1, 1), sub_shape=sub,
                        extern_key=key)

    def test_extern_malformed_header_refused(self):
        """The reader stores no dimension for a malformed line 2; the
        spin-orbital Extern check re-reads it and names the problem."""
        for key in ("Extern", "EXTERN"):
            with self.subTest(key=key), \
                    self.assertRaisesRegex(ValueError, r"extern\.dat: line 2 must be the number "
                                           r"of orbitals \(a positive integer\), got 'abc'"):
                _rpa_so(route="extern", extern_header_dim="abc", extern_key=key)

    def test_non_so_malformed_header_accepted_as_before(self):
        """Without spin-orbital mode a malformed line 2 (ignored by
        ``read_w90``) is still accepted, with the same result."""
        import hwave.qlmsio.read_input_k as read_input_k
        res = {}
        for line2 in ("2", "abc"):
            with _TmpDirs() as t:
                d = t.new()
                for f in ("geom.dat", "transfer.dat", "coulombintra.dat"):
                    shutil.copy(os.path.join(_SPIN_DIR, f), os.path.join(d, f))
                p = os.path.join(d, "coulombintra.dat")
                with open(p) as f:
                    lines = f.readlines()
                lines[1] = line2 + "\n"
                with open(p, "w") as f:
                    f.writelines(lines)
                inter = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
                         "CoulombIntra": "coulombintra.dat"}
                with _quiet():
                    reader = read_input_k.QLMSkInput(
                        {"path_to_input": d, "interaction": dict(inter, path_to_input=d)})
                dims = reader.get_param("ham").source_dimensions
                self.assertEqual(dims["coulombintra"], None if line2 == "abc" else 2)
                res[line2] = _run(d, "RPA", _E2_PARAM, inter, {"chiq": "chiq.npz"}, False,
                                  {"calc_scheme": "general", "calc_type": "ring"})["chiq"]
        np.testing.assert_array_equal(res["abc"], res["2"])

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
        with self.assertRaisesRegex(ValueError, r"^Extern file 'e\.dat': its header declares 4 "
                                    r"orbitals, but norb_phys = 2\. In spin-orbital mode "
                                    r"Extern uses PHYSICAL orbital indices \(norb_phys = 2\), "
                                    r"while the Transfer file uses the 2·norb_phys = 4 "
                                    r"spin-orbital indices; write the field per physical "
                                    r"orbital \(the sign per spin comes from σ\^z\)\.$"):
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


# ------------------------------------------------- the Eliashberg front end
#
# The two-orbital square-lattice model of tests/test_sc_spin_orbital.py
# (``S._write_model``; on-site SU(2) interaction set, optional Rashba-like
# spin-orbit hopping ``soc``), with a Zeeman field along either route.


class _Stop(Exception):
    """Ends a run once the array under test has been captured."""


@contextlib.contextmanager
def _quiet_info(level=logging.INFO):
    """Silence ``level`` and below (default INFO: WARNINGs still reach
    ``assertLogs``)."""
    prev = logging.root.manager.disable
    logging.disable(level)
    try:
        yield
    finally:
        logging.disable(prev)


def _sc_model(d, soc=0.0, transfer_zeeman=0.0, extern=None, extern_rows=None,
              extern_header_dim=S.NP, key="Extern", types=S._SU2_TYPES):
    """Write the model into ``d``; returns the interaction file map.

    ``transfer_zeeman = h`` adds ``+h`` / ``-h`` to the on-site up / down
    Transfer entries of every orbital (Route 1); ``extern = h`` writes an
    Extern file ``H_aa(0) = h`` in physical indices (Route 2), or
    ``extern_rows`` verbatim (1-based), under the key ``key`` with line 2
    ``extern_header_dim``."""
    files = S._write_model(d, so=True, soc=soc, types=types)
    if transfer_zeeman:
        path = os.path.join(d, "transfer.dat")
        dim, _nr, _deg, rows = _read_w90_rows(path)
        # 1-based interleaved file index 2a + s + 1 -> the sign of the field
        onsite = {2 * a + s + 1: (1 if s == 0 else -1)
                  for a in range(S.NP) for s in (0, 1)}
        hit = set()
        new_rows = []
        for r in rows:
            if r[:3] == (0, 0, 0) and r[3] == r[4] and r[3] in onsite:
                r = r[:5] + (r[5] + onsite[r[3]] * transfer_zeeman, r[6])
                hit.add(r[3])
            new_rows.append(r)
        assert hit == set(onsite), hit          # every on-site row exists
        S._wfile(path, "Transfer", dim, new_rows)
    if extern is not None or extern_rows is not None:
        if extern_rows is None:
            extern_rows = [(0, 0, 0, a + 1, a + 1, extern, 0.0) for a in range(S.NP)]
        _write_w90(d, extern_header_dim, extern_rows)
        files = dict(files, **{key: "extern.dat"})
    return files


def _sc_input(d, files, eli=None, coeff_extern=1.0, coeff_key="coeff_extern", **kw):
    inp = S._input(d, files, True,
                   dict({"chi0q_mode": "calc", "solver_mode": "eigenvalue"}, **(eli or {})),
                   **kw)
    if coeff_extern is not None:
        inp["mode"]["param"][coeff_key] = coeff_extern
    return inp


def _capture_sc_band(soc=0.0, transfer_zeeman=0.0, extern=None, coeff_extern=1.0,
                     key="Extern", coeff_key="coeff_extern", silence=logging.INFO):
    """The spin-block epsilon(k) ``(nd, nd, Nx, Ny, Nz)`` the static
    spin-orbital hwave_sc run diagonalizes (``chi0q_mode = "calc"``);
    log records up to ``silence`` are dropped."""
    import hwave.sc as sc
    seen = {}

    def spy_eig(eps, *a, **kw):
        seen["eps"] = np.array(eps, copy=True)
        raise _Stop()

    with _TmpDirs() as t:
        d = t.new()
        files = _sc_model(d, soc=soc, transfer_zeeman=transfer_zeeman, extern=extern, key=key)
        inp = _sc_input(d, files, coeff_extern=coeff_extern, coeff_key=coeff_key)
        with _quiet_info(silence), mock.patch.object(sc, "_calc_eigenvalues",
                                                     side_effect=spy_eig):
            try:
                sc.calc_eliashberg(inp)
            except _Stop:
                pass
    return seen["eps"]


def _to_rpa_layout(eps):
    """``(nd, nd, Nx, Ny, Nz)`` -> the RPA ``(nvol, nd, nd)`` layout."""
    nd = eps.shape[0]
    return np.moveaxis(eps.reshape(nd, nd, -1), -1, 0)


def _capture_sc_and_rpa_band(extern, soc):
    """``(epsilon_k of hwave_sc, H0(k) of the spin-orbital RPA)`` with the
    Extern route; the RPA band must be the spinful one (``soc != 0``)."""
    import hwave.qlms as qlms
    import hwave.solver.rpa as rpa_mod
    eps = _capture_sc_band(soc=soc, extern=extern)
    seen = {}
    real_eps = rpa_mod.RPA._calc_epsilon_k

    def spy_rpa(self_, green_info):
        real_eps(self_, green_info)
        seen["H0"] = np.array(self_.H0_k, copy=True)
        seen["spin_mode"] = self_.spin_mode
        raise _Stop()

    with _TmpDirs() as t:
        d = t.new()
        files = _sc_model(d, soc=soc, extern=extern)
        inp = _sc_input(d, files, out="out_rpa")
        inp["mode"]["calc_scheme"] = "general"
        inp["file"]["input"]["path_to_input"] = d
        del inp["eliashberg"]
        inp["log"] = {"print_level": 0}
        with _quiet(), mock.patch.object(rpa_mod.RPA, "_calc_epsilon_k", autospec=True,
                                         side_effect=spy_rpa):
            try:
                qlms.run(input_dict=inp)
            except _Stop:
                pass
    assert seen["spin_mode"] == "spinful", seen["spin_mode"]
    nd = eps.shape[0]
    return eps, seen["H0"].reshape(-1, nd, nd)


_HEAVY_MOCKS = ("_build_bond_green", "_calc_chi0q_internal", "_load_chi0q",
                "_load_flex_so", "_so_static_vertex", "_so_rpa_vertex",
                "_calc_eigenvalues")


def _run_sc(mode, extern_rows, extern_header_dim=S.NP, key="Extern"):
    """A static spin-orbital hwave_sc run in ``chi0q_mode = mode`` with the
    given Extern file; every heavy step raises, so a refusal must come
    first. Raises what the run raises."""
    import hwave.sc as sc
    with _TmpDirs() as t:
        d = t.new()
        files = _sc_model(d, extern_rows=extern_rows, extern_header_dim=extern_header_dim,
                          key=key)
        inp = _sc_input(d, files, eli={"chi0q_mode": mode})
        with contextlib.ExitStack() as st:
            st.enter_context(_quiet())
            for name in _HEAVY_MOCKS:
                st.enter_context(mock.patch.object(
                    sc, name, side_effect=AssertionError("heavy work reached: " + name)))
            try:
                sc.calc_eliashberg(inp)
            finally:
                out = inp["file"]["output"]["path_to_output"]
                assert not (os.path.exists(out) and os.listdir(out)), os.listdir(out)


def _run_sc_dynamic_calc(extern=None, transfer_zeeman=0.0, transfer_onsite_symmetric=True,
                         soc=0.3):
    """A full spin-orbital DYNAMIC hwave_sc run (``chi0q_mode = "calc"``,
    Nmat = 8); returns the eigenvalues it found. Logging is not silenced
    (callers check the WARNINGs). ``transfer_onsite_symmetric`` asserts the
    Transfer file has no spin-dependent on-site entry (Route 1 off)."""
    import hwave.sc as sc
    assert transfer_onsite_symmetric == (transfer_zeeman == 0.0)
    cap = {}
    real_solve = sc._so_solve_eigenvalue

    def solve(*a, **kw):
        out = real_solve(*a, **kw)
        cap["vals"] = np.array(out[0])
        return out

    with _TmpDirs() as t:
        d = t.new()
        files = _sc_model(d, soc=soc, transfer_zeeman=transfer_zeeman, extern=extern)
        inp = _sc_input(d, files, eli={"frequency": "dynamic", "num_eigenvalues": 4,
                                       "spectral_shift": "auto"}, nmat=8)
        with _quiet_info(), mock.patch.object(sc, "_so_solve_eigenvalue", side_effect=solve):
            sc.calc_eliashberg(inp)
    return cap["vals"]


class TestScBand(unittest.TestCase):
    def test_epsilon_k_equals_rpa_with_extern(self):                   # soc != 0: spinful RPA
        eps_sc, h0_rpa = _capture_sc_and_rpa_band(extern=0.1, soc=0.35)
        np.testing.assert_allclose(_to_rpa_layout(eps_sc), h0_rpa, rtol=0, atol=1e-12)
        # spin-mixing hopping present, and the field is in the band
        n = S.NP
        self.assertGreater(np.abs(h0_rpa[:, :n, n:]).max(), 0.1)
        no_field = _to_rpa_layout(_capture_sc_band(soc=0.35))
        diff = h0_rpa - no_field
        np.testing.assert_allclose(diff[:, :n, :n], np.broadcast_to(0.1 * np.eye(n), diff[:, :n, :n].shape),
                                   rtol=0, atol=1e-14)
        np.testing.assert_allclose(diff[:, n:, n:], -diff[:, :n, :n], rtol=0, atol=1e-14)

    def test_routes_add_and_warn(self):
        with self.assertLogs("hwave_sc", "WARNING") as cm:
            both = _capture_sc_band(transfer_zeeman=0.1, extern=0.1)
        self.assertTrue(any("the two fields add" in m for m in cm.output), cm.output)
        self.assertTrue(any("orbital 1" in m for m in cm.output), cm.output)
        np.testing.assert_allclose(both, _capture_sc_band(transfer_zeeman=0.2), rtol=0, atol=1e-14)

    def test_single_route_does_not_warn_and_states_sign(self):
        for kw in ({"extern": 0.1}, {"transfer_zeeman": 0.1}):
            with self.subTest(**kw), self.assertLogs("hwave_sc", "INFO") as cm:
                _capture_sc_band(silence=logging.DEBUG, **kw)
            self.assertFalse(any("the two fields add" in m for m in cm.output))
            sign = [m for m in cm.output
                    if "Zeeman sign: +h raises spin-up (favours spin-down occupation)" in m]
            self.assertEqual(len(sign), 1 if "extern" in kw else 0, cm.output)

    def test_extern_and_transfer_routes_same_band(self):
        np.testing.assert_allclose(_capture_sc_band(extern=0.1),
                                   _capture_sc_band(transfer_zeeman=0.1), rtol=0, atol=1e-14)

    def test_zero_coeff_extern_equals_transfer_band(self):
        np.testing.assert_array_equal(_capture_sc_band(extern=0.1, coeff_extern=0.0),
                                      _capture_sc_band())
        # and an absent coeff_extern defaults to 0.0
        np.testing.assert_array_equal(_capture_sc_band(extern=0.1, coeff_extern=None),
                                      _capture_sc_band())

    def test_mixed_case_keys(self):
        ref = _capture_sc_band(extern=0.1)
        self.assertGreater(np.abs(ref - _capture_sc_band()).max(), 0.05)
        np.testing.assert_array_equal(
            _capture_sc_band(extern=0.1, key="EXTERN", coeff_key="Coeff_Extern"), ref)

    def test_extern_index_out_of_range_load_and_flex(self):
        for mode in ("calc", "load", "flex"):
            for key in ("Extern", "EXTERN"):
                with self.subTest(mode=mode, key=key), self.assertRaisesRegex(
                        ValueError, r"Extern file '.*extern\.dat': orbital index 3 is outside "
                                    r"1\.\.2\. In spin-orbital mode"):
                    _run_sc(mode, extern_rows=[(0, 0, 0, 3, 3, 0.1, 0.0)], key=key)
            with self.subTest(mode=mode, header=4), self.assertRaisesRegex(
                    ValueError, r"Extern file '.*extern\.dat': its header declares 4 orbitals, "
                                r"but norb_phys = 2\. In spin-orbital mode"):
                _run_sc(mode, extern_header_dim=4, extern_rows=[(0, 0, 0, 1, 1, 0.1, 0.0)])
            with self.subTest(mode=mode, header="abc"), self.assertRaisesRegex(
                    ValueError, r"extern\.dat: line 2 must be the number of orbitals"):
                _run_sc(mode, extern_header_dim="abc", extern_rows=[(0, 0, 0, 1, 1, 0.1, 0.0)])

    def test_no_false_ignored_warning_dynamic(self):
        with self.assertNoLogs(level="WARNING"):
            lam_ext = _run_sc_dynamic_calc(extern=0.1, transfer_onsite_symmetric=True)
        lam_tr = _run_sc_dynamic_calc(transfer_zeeman=0.1, transfer_onsite_symmetric=False)
        lam_0 = _run_sc_dynamic_calc()
        self.assertLessEqual(_assignment_max(lam_ext, lam_tr), 1e-10 * np.abs(lam_tr).max())
        # the field reaches the dynamic solve
        self.assertGreater(_assignment_max(lam_ext, lam_0), 1e-6)

    def test_no_false_ignored_warning_static_flex(self):
        """The static chi0q_mode = "flex" path issues no "Extern is ignored"
        WARNING either (the spin-orbital FLEX band includes Extern); the run
        is stopped where the FLEX archives would be read."""
        import hwave.sc as sc
        with _TmpDirs() as t:
            d = t.new()
            files = _sc_model(d, soc=0.3, extern=0.1)
            inp = _sc_input(d, files, eli={"chi0q_mode": "flex"})
            with self.assertNoLogs(level="WARNING"), _quiet_info(), \
                    mock.patch.object(sc, "_load_flex_so", side_effect=_Stop()):
                with self.assertRaises(_Stop):
                    sc.calc_eliashberg(inp)


def _assignment_max(a, b):
    """Max ``|a_i - b_pi(i)|`` over the optimal one-to-one assignment."""
    from scipy.optimize import linear_sum_assignment
    a, b = np.asarray(a).ravel(), np.asarray(b).ravel()
    assert a.size == b.size, (a.size, b.size)
    cost = np.abs(a[:, None] - b[None, :])
    r, c = linear_sum_assignment(cost)
    return float(cost[r, c].max())


_SPECTRA = {}


def _sc_static_dense(h, soc=0.0, route="extern"):
    """The production static spin-orbital run (dense branch, every
    eigenvalue of the antisymmetric space) with a field ``h``:
    ``{"evs", "M", "shape"}`` -- the eigenvalues it reports and the dense
    matrix of the kernel operator it built. Cached per arguments."""
    import hwave.sc as sc
    key = (h, soc, route)
    if key in _SPECTRA:
        return _SPECTRA[key]
    cap = {}
    real_solve = sc._so_solve_eigenvalue

    def solve(A, basis, *a, **kw):
        cap["A"] = A
        out = real_solve(A, basis, *a, **kw)
        cap["evs"] = np.array(out[0])
        return out

    Nx, Ny, Nz = S.CELL
    shape = (S.ND, S.ND, Nx, Ny, Nz)
    n = int(np.prod(shape))
    with _TmpDirs() as t:
        d = t.new()
        kw = {"extern": h} if route == "extern" else {"transfer_zeeman": h}
        files = _sc_model(d, soc=soc, **kw)
        inp = _sc_input(d, files, eli={"num_eigenvalues": n}, nmat=32)
        with _quiet(), mock.patch.object(sc, "_so_solve_eigenvalue", side_effect=solve):
            sc.calc_eliashberg(inp)
    M = cap["A"].matmat(np.eye(n, dtype=complex))
    _SPECTRA[key] = {"evs": cap["evs"], "M": M, "shape": shape}
    return _SPECTRA[key]


def _sc_static_full_spectrum(h, soc=0.0):
    out = _sc_static_dense(h, soc)
    dim = np.linalg.matrix_rank(S._dense_P(out["shape"]))
    assert out["evs"].size == dim, (out["evs"].size, dim)    # every eigenvalue
    return out["evs"]


def _sector_spectrum(p, h, soc=0.0):
    """Eigenvalues of the kernel compressed onto the antisymmetric gaps of
    pair-``S_z = p`` (``p = +1``: up-up, ``0``: up-down and down-up, ``-1``:
    down-down), from an orthonormal basis of ``P Pi_p`` built here."""
    out = _sc_static_dense(h, soc)
    shape = out["shape"]
    nd = shape[0]
    spin = np.arange(nd) // (nd // 2)              # spin-block order: 0 = up
    sz = (1 - 2 * spin)[:, None] + (1 - 2 * spin)[None, :]     # 2 * (s_z + t_z)
    mask = np.broadcast_to((sz == 2 * p)[:, :, None, None, None], shape).ravel()
    P = S._dense_P(shape)
    Pi = np.diag(mask.astype(float))
    w, U = np.linalg.eigh(Pi @ P @ Pi)
    Q = U[:, w > 0.5]
    return np.linalg.eigvals(Q.conj().T @ out["M"] @ Q)


def _sc_eigs(route, h=_H, soc=0.35):
    return _sc_static_dense(h, soc, route)["evs"]


class TestScSpectrum(unittest.TestCase):                                # 0-S1
    def test_zeeman_spectrum_is_union_of_pair_sz_sectors(self):
        lam = _sc_static_full_spectrum(h=0.1, soc=0.0)
        sectors = {p: _sector_spectrum(p, h=0.1) for p in (+1, 0, -1)}
        expect = np.concatenate([sectors[p] for p in (+1, 0, -1)])
        self.assertLessEqual(_assignment_max(lam, expect), 1e-8 * max(1, np.abs(lam).max()))
        # the field splits the equal-spin sectors (without it they coincide)
        self.assertGreater(_assignment_max(sectors[+1], sectors[-1]), 1e-4)
        z = {p: _sector_spectrum(p, h=0.0) for p in (+1, -1)}
        self.assertLessEqual(_assignment_max(z[+1], z[-1]), 1e-10 * max(1, np.abs(lam).max()))

    def test_extern_and_transfer_routes_same_spectrum(self):
        a, b = _sc_eigs("extern"), _sc_eigs("transfer")
        self.assertEqual(a.size, b.size)
        self.assertLessEqual(_assignment_max(a, b), 1e-12 * max(1, np.abs(b).max()))
        self.assertGreater(_assignment_max(a, _sc_eigs("extern", h=0.0)), 1e-6)


if __name__ == "__main__":
    unittest.main()
