"""hwave_sc in spin-orbital mode with bond-resolved channels and a Green
function read from ``[eliashberg] bond_green`` (issue #270): the
configuration, the Green function built from the green.npz of a
spin-orbital bond FLEX run, the refusals of that file, the tail policy and
the S_z-leakage gates, the content and k/frequency-axis gates and the
memory admission of a borrowed Green function (with their deliberate-break
rows).

The producer files are written once (:func:`setUpModule`, converged
spin-orbital bond FLEX runs at the consumer's state point) and removed in
:func:`tearDownModule`. The rows vary the consumer (the prepared run) or a
copy of the file, and call ``sc._so_bond_green`` directly (no solve); the
on-site row, whose run record has no band, calls ``sc._so_borrowed_green``.
"""

import contextlib
import gc
import io
import logging
import os
import shutil
import tempfile
import unittest
import zipfile
from unittest import mock

import numpy as np
from requests.structures import CaseInsensitiveDict

from hwave.solver.eliashberg import static as sc_static
from hwave.solver.eliashberg import io as sc_io
from hwave.solver.eliashberg import so_solve as sc_so
from hwave.solver.eliashberg import bond as sc_bond
from hwave.solver.eliashberg import config as sc_config
from hwave.solver import fingerprint, so_bond, so_bond_archive
from tests import so_bond_archive_fixtures as AF
from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X
from tests import test_sc_spin_orbital_bond as T
from tests.so_bond_oracles import CELL
from tests.test_sc_spin_orbital_bond import (_gaps, _miss, _reference_kernel_of,
                                             _reference_vertex)
from tests.test_so_entry_rules import _PayloadReads

_W = "[mode] enable_spin_orbital = true (hwave_sc, issue #83)"
_LABEL = "hwave FLEX (spin-orbital, bond channels)"
_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond_green")
#: ZEEMAN_EXTERN with every interaction type: its band conserves S_z, its
#: interaction (PairLift, ...) does not
ZPL = F.ZEEMAN_EXTERN._replace(name="ZPL", interactions=F._ALL_TYPES)
#: the producers of the module
_PRODUCERS = {"ZEEMAN_SZ_ALL": X.ZEEMAN_SZ_ALL, "SOC_ALL_TYPES": F.SOC_ALL_TYPES, "ZPL": ZPL}
#: the producers that also write sigma.npz
_SIGMA = ("ZEEMAN_SZ_ALL", "SOC_ALL_TYPES")
#: the producer of the Nmat 8 memory row (ZEEMAN_SZ_ALL at T 0.5, Nmat 8)
_N8 = "ZEEMAN_SZ_ALL_N8"
_GREENS = {}
_ROOT = []
_DROP = object()


def setUpModule():
    """The green.npz (and, for :data:`_SIGMA`, sigma.npz) of a converged
    spin-orbital bond FLEX run of each producer, at T 0.5, Nmat 32 and the
    fixture's filling; and :data:`_N8`, ZEEMAN_SZ_ALL at Nmat 8."""
    root = tempfile.mkdtemp(prefix="sc_so_bond_green_")
    _ROOT.append(root)
    with contextlib.redirect_stdout(io.StringIO()):
        for name, fx in _PRODUCERS.items():
            d = os.path.join(root, name)
            os.makedirs(d)
            _GREENS[name] = X.write_sc_green(d, fx, sigma=name in _SIGMA)
        d = os.path.join(root, _N8)
        os.makedirs(d)
        _GREENS[_N8] = X.write_sc_green(d, X.ZEEMAN_SZ_ALL, param={"Nmat": 8})


def tearDownModule():
    for root in _ROOT:
        shutil.rmtree(root, True)
    _ROOT.clear()
    _GREENS.clear()
    _BASE.clear()
    _KERNELS.clear()


def _scratch(test):
    d = tempfile.mkdtemp(prefix="sc_so_bond_green_row_")
    test.addCleanup(shutil.rmtree, d, True)
    return d


def _copy(test, name):
    """A copy of the module green of producer ``name`` in a fresh directory."""
    path = os.path.join(_scratch(test), "green.npz")
    shutil.copyfile(_GREENS[name], path)
    return path


def _members(path):
    with np.load(path, allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def _edited(test, path, **members):
    """A copy of the archive at ``path`` with ``members`` replaced
    (``_DROP`` removes a member); nothing is re-stamped."""
    m = _members(path)
    for k, v in members.items():
        if v is _DROP:
            del m[k]
        else:
            m[k] = v
    out = os.path.join(_scratch(test), "green.npz")
    np.savez(out, **m)
    return out


def _restamped(G):
    return np.str_(fingerprint.green_digest(G, fingerprint.GREEN_LAYOUT))


def _perturbed(test, path, rel, *, restamp):
    """A copy of ``path`` whose ``green`` is ``G + rel max|G| P`` (``P`` a
    seed-5 complex pattern of unit max-norm, nonzero only S_z-off-sector),
    its ``fp_green`` re-stamped when ``restamp``."""
    G = _members(path)["green"]
    rng = np.random.default_rng(5)
    pat = rng.standard_normal(G.shape) + 1j * rng.standard_normal(G.shape)
    pat /= np.max(np.abs(pat))
    n = G.shape[-1] // 2
    P = np.zeros_like(G)
    P[..., :n, n:] = pat[..., :n, n:]
    P[..., n:, :n] = pat[..., n:, :n]
    G2 = G + rel * np.max(np.abs(G)) * P
    test.assertGreater(float(np.max(np.abs(G2 - G))), 0.0)
    extra = {"fp_green": _restamped(G2)} if restamp else {}
    return _edited(test, path, green=G2, **extra)


def _mismatch(path):
    """The ``fp_green`` mismatch text of the edited, not re-stamped file at
    ``path``."""
    m = _members(path)
    return ("'{}' (written by {}): fp_green mismatch: the Green function read from it has "
            "digest {}..., the file records {}... (both at fingerprint_version 2); its "
            "green member is not the one the run that wrote it recorded. Regenerate it "
            "with {}.".format(path, _LABEL, str(_restamped(m["green"]))[:12],
                              str(m["fp_green"])[:12], _LABEL))


def _consumer(fixture, *, param=None, eliashberg=None):
    """The prepared hwave_sc run of ``fixture`` (beta 2.0, Nmat 32;
    ``param``, ``eliashberg`` on top)."""
    return T._prepare_admission(fixture, eliashberg=eliashberg, param=param)


_BASE = {}


def _base(name):
    if name not in _BASE:
        _BASE[name] = _consumer(_PRODUCERS[name])
    return _BASE[name]


def _borrow(p, path):
    return sc_bond._so_bond_green(p.epsilon_k, p.so_run, norb_phys=p.norb_phys, filling=p.filling,
                             beta=p.beta, nmat=p.nmat, coeff_tail=p.cfg.coeff_tail,
                             bond_green=path, coeff_tail_source=p.cfg.coeff_tail_source)


@contextlib.contextmanager
def _records():
    """The INFO-and-above records of the ``hwave_sc`` and
    ``hwave.solver.so_bond`` loggers while active."""
    seen = []

    class H(logging.Handler):
        def emit(self, record):
            seen.append(record)
    h = H(logging.INFO)
    loggers = [sc_static.logger, so_bond.logger]
    levels = [lg.level for lg in loggers]
    for lg in loggers:
        lg.addHandler(h)
        lg.setLevel(logging.INFO)
    try:
        yield seen
    finally:
        for lg, lv in zip(loggers, levels):
            lg.removeHandler(h)
            lg.setLevel(lv)


def _count(records, template, *args):
    text = template % args if args else template
    return sum(1 for r in records if r.getMessage() == text)


def _warnings(records):
    return [r.getMessage() for r in records if r.levelno >= logging.WARNING]


# ------------------------------------------------------------ configuration
class TestBorrowedConfig(unittest.TestCase):
    """``bond_green`` is accepted on the spin-orbital bond path: read from
    the case-insensitive ``[eliashberg]`` table as a path; the empty string
    keeps the non-spin-orbital text; the remaining refusals keep their
    order."""

    def _parse(self, eli):
        return sc_config._read_so_bond_config(CaseInsensitiveDict(eli))

    def test_path_read(self):
        cfg = self._parse({"bond_channels": True, "Bond_Green": "g.npz"})
        self.assertEqual(cfg.green, "g.npz")
        self.assertIsNone(self._parse({"bond_channels": True}).green)
        # the band-only constructor of seven fields stays valid
        self.assertIsNone(sc_config._SOBondConfig(True, "auto", None, None, 1.0, "default", True).green)

    def test_empty_string(self):
        with self.assertRaises(ValueError) as cm:
            self._parse({"bond_channels": True, "bond_green": ""})
        self.assertEqual(str(cm.exception),
                         "[eliashberg] bond_green must be a path to a Green-function npz "
                         "file; got an empty string. Omit the key to build the bare Green "
                         "function from the transfer Hamiltonian.")
        # the non-spin-orbital reader keeps the same text
        with self.assertRaises(ValueError) as cm2:
            sc_config._read_bond_config({"bond_channels": True, "bond_green": ""})
        self.assertEqual(str(cm2.exception), str(cm.exception))

    def test_refusal_order(self):
        on = {"bond_channels": True}
        with self.assertRaises(ValueError) as cm:
            self._parse(dict(on, bond_green="g.npz", bond_diagnostics=True))
        self.assertEqual(str(cm.exception), sc_config._SO_BOND_DIAGNOSTICS_REFUSAL)


# ------------------------------------------------------------ accepted files
class TestBorrowedGreen(unittest.TestCase):
    """The Green function of a static spin-orbital bond run read from the
    green.npz of the spin-orbital bond FLEX run at its state point: the
    carrier holds the file's G (no copy, no tail), the record holds the
    file's digest, mu and state."""

    def _run(self, name):
        path = _copy(self, name)
        got = []
        real = sc_io._so_borrowed_green

        def spy(*a, **kw):
            got.append(real(*a, **kw))
            return got[-1]
        with _records() as rec, mock.patch.object(sc_bond, "_so_borrowed_green", side_effect=spy):
            so_green = _borrow(_base(name), path)
        self.assertEqual(len(got), 1)
        return path, so_green, got[0], rec

    def test_clean_conserving(self):
        path, g, borrowed, rec = self._run("ZEEMAN_SZ_ALL")
        m = _members(path)
        c = g.carrier
        self.assertEqual(c.source, "external_npz")
        self.assertEqual(g.source, "external_npz")
        self.assertEqual(g.fp_green, str(m["fp_green"]))
        self.assertIsNone(g.fp_green_source)
        self.assertEqual(g.leak, 0.0)
        self.assertEqual(c.coeff_tail_applied, 0.0)
        self.assertEqual(c.coeff_tail_requested, 1.0)
        self.assertIsNone(c.green0_tail)
        self.assertFalse(c.full_sc.flags.writeable)
        self.assertFalse(c.deflated_kw.flags.writeable)
        self.assertTrue(np.array_equal(c.deflated_kw, m["green"]))
        self.assertEqual(g.mu, float(m["state_mu"]))
        with np.load(path) as z:
            self.assertEqual(g.state, fingerprint._file_state(z))
        self.assertIsNone(g.eigenvalues)
        self.assertIs(g.conserving, True)
        self.assertIs(g.scf_converged, True)
        self.assertIs(c.deflated_kw, borrowed.used)
        self.assertFalse(borrowed.projected)
        self.assertEqual(_warnings(rec), [])
        self.assertEqual(_count(rec, sc_bond._SO_BOND_TAIL_DEFAULT_INFO), 1)
        self.assertEqual(_count(rec, sc_bond._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
                                "0.000e+00", "no"), 1)

    def test_clean_non_conserving(self):
        path, g, borrowed, rec = self._run("SOC_ALL_TYPES")
        self.assertIsNone(g.leak)
        self.assertIsNone(borrowed.leak)
        self.assertIs(g.conserving, False)
        self.assertEqual(g.fp_green, str(_members(path)["fp_green"]))
        self.assertIs(g.carrier.deflated_kw, borrowed.used)
        self.assertEqual(_warnings(rec), [])
        self.assertEqual(_count(rec, sc_bond._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
                                "not checked (the input does not conserve S_z)", "no"), 1)


class TestBorrowedRecord(unittest.TestCase):
    """The gap archive of a borrowed G records ``fp_green_source`` when its
    kernel used the projection: the writer forwards it and the member
    policy of the bond gap archive admits it then only."""

    def test_writer_and_member_policy(self):
        p = _base("ZEEMAN_SZ_ALL")
        with np.load(_GREENS["ZEEMAN_SZ_ALL"]) as z:
            state = fingerprint._file_state(z)
            fp_file = str(z["fp_green"])
        fp_used = "0" * 64
        plain = sc_so._so_writer_members("sc_so_bond_static", "test", p.so_run, fp_file, state)
        self.assertNotIn("fp_green_source", plain)
        members = sc_so._so_writer_members("sc_so_bond_static", "test", p.so_run, fp_used, state,
                                        fp_green_source=fp_file)
        self.assertEqual(str(members["fp_green_source"]), fp_file)
        self.assertEqual(str(members["fp_green"]), fp_used)
        self.assertEqual(set(members) - set(plain), {"fp_green_source"})
        for listed in (False, True):
            base = sc_so._so_bond_gap_members(listed)
            self.assertNotIn("fp_green_source", base)
            self.assertEqual(sc_so._so_bond_gap_members(listed, projected=True),
                             base | {"fp_green_source"})


# ------------------------------------------------------------ refused files
class TestBorrowedRefusals(unittest.TestCase):
    """The refusals of a ``bond_green`` file, exact texts and exception
    classes, in the order of the checks; every refusal before the payload
    read reads no payload (``_PayloadReads``)."""

    maxDiff = None

    def _refuse(self, p, path, cls=fingerprint.SOArchiveError):
        with _records() as rec, _PayloadReads() as spy:
            with self.assertRaises(Exception) as cm:
                _borrow(p, path)
        self.assertIs(type(cm.exception), cls, repr(cm.exception))
        return cm.exception, spy.reads, rec

    def _row(self, p, path, text, cls=fingerprint.SOArchiveError, read=False):
        e, reads, rec = self._refuse(p, path, cls)
        self.assertEqual(str(e), text)
        self.assertEqual(reads, [(os.path.basename(path), "green")] if read else [])
        self.assertEqual(_count(rec, sc_bond._SO_BOND_GREEN_UNCONVERGED, path), 0)
        return e

    def test_texts(self):
        self.assertEqual(sc_io.T_BG_MISSING, "[eliashberg] bond_green file not found: {path}")
        self.assertEqual(sc_io.T_BG_SHAPE,
                         "'{path}' is malformed or truncated: green shape differs from this "
                         "run ({stored} vs {expected}). Regenerate it with {label}.")
        self.assertEqual(sc_io.T_BG_NONFINITE,
                         "'{path}' is malformed or truncated: green has non-finite values at "
                         "frequency index {n}. Regenerate it with {label}.")
        self.assertEqual(sc_io.T_BG_MARKER,
                         "'{path}' is malformed or truncated: {member} is missing or older "
                         "than this build requires. Regenerate it with {label}.")
        self.assertEqual(sc_io.T_BG_KGRID,
                         "'{path}' is malformed or truncated: wavevector_index is not the k "
                         "grid {shape} of this run (CellShape / SubShape). Regenerate it with "
                         "{label}.")
        self.assertEqual(sc_io.T_BG_SCF,
                         "'{path}' is malformed or truncated: scf_converged is not a 0-d "
                         "boolean. Regenerate it with {label}.")
        self.assertEqual(sc_io._SO_BOND_GREEN_IS_FLEX_OUTPUT,
                         _W + ": [eliashberg] bond_green = '{path}' is the {kind} file of a "
                         "FLEX run, not a Green function; set bond_green to the green.npz "
                         "written by that run ([file.output] green).")
        self.assertEqual(sc_io._SO_BOND_GREEN_ONSITE_REMEDY,
                         "[eliashberg] bond_green needs the green.npz of a spin-orbital FLEX "
                         "run with [mode.param] longitudinal_bond_channels = true (the "
                         "bond-resolved run); re-run FLEX with it.")
        self.assertEqual(sc_io._SO_BOND_GREEN_IS_ARCHIVE,
                         _W + ": [eliashberg] bond_green = '{path}' is the bond archive of a "
                         "spin-orbital FLEX run, not a Green function; set bond_green to its "
                         "companion green file '{green}'.")
        self.assertEqual(sc_io._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN,
                         _W + ": [eliashberg] bond_green = '{path}' is the bond archive of a "
                         "spin-orbital FLEX run, not a Green function, and it records no "
                         "readable companion green file; set bond_green to the green.npz "
                         "written by that FLEX run ([file.output] green).")
        self.assertEqual(sc_io._SO_BOND_GREEN_NONSO_ARCHIVE,
                         _W + ": [eliashberg] bond_green = '{path}' is a bond archive, not a "
                         "Green function; set bond_green to the green.npz of a spin-orbital "
                         "FLEX run with [mode.param] longitudinal_bond_channels = true.")
        self.assertEqual(sc_bond._SO_BOND_GREEN_UNCONVERGED,
                         "[eliashberg] bond_green '%s' records scf_converged = false: the FLEX "
                         "run did not converge; its Green function is used as stored.")
        self.assertEqual(sc_bond._SO_BOND_GREEN_ACCEPTED_INFO,
                         "[eliashberg] bond_green: '%s' (written by %s) passed the archive "
                         "checks; S_z leak %s; projected: %s.")

    def test_file_kind(self):
        p = _base("ZEEMAN_SZ_ALL")
        d = _scratch(self)
        missing = os.path.join(d, "absent.npz")
        self._row(p, missing, "[eliashberg] bond_green file not found: " + missing,
                  cls=FileNotFoundError)
        npy = os.path.join(d, "green.npy")
        np.save(npy, np.zeros(3))
        self._row(p, npy, "[eliashberg] bond_green file '{}' is not an npz archive (it holds "
                  "a single .npy array); point the configuration at the .npz file the FLEX "
                  "run wrote.".format(npy), cls=ValueError)
        # a non-spin-orbital green.npz: the structure step
        path = _edited(self, _GREENS["ZEEMAN_SZ_ALL"], spin_orbital=_DROP)
        self._row(p, path, "'{}' is malformed or truncated: marker spin_orbital is missing "
                  "(absent vs present). Regenerate it with {}.".format(path, _LABEL))

    def test_archive_index(self):
        p = _base("ZEEMAN_SZ_ALL")
        d = _scratch(self)
        inputs = AF.archive_inputs(green_file=os.path.join("sub", "green.npz"))
        index = os.path.join(d, "longitudinal_bond.npz")
        so_bond_archive.write_so_bond_archive(index, output_layout="npz", log=_LOG,
                                              **inputs.kwargs)
        companion = os.path.join(d, "sub", "green.npz")
        # the companion is named only when it is a readable regular file
        self._row(p, index, sc_io._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN.format(path=index),
                  cls=ValueError)
        os.makedirs(companion)
        self._row(p, index, sc_io._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN.format(path=index),
                  cls=ValueError)
        os.rmdir(companion)
        shutil.copyfile(_GREENS["ZEEMAN_SZ_ALL"], companion)
        self._row(p, index, sc_io._SO_BOND_GREEN_IS_ARCHIVE.format(path=index, green=companion),
                  cls=ValueError)
        self.assertIn("/sub/green.npz'.", str(self._refuse(p, index, ValueError)[0]))
        no_file = _edited(self, index, green_file=_DROP)
        empty = _edited(self, index, green_file=np.str_(""))
        for path in (no_file, empty):
            with self.subTest(path=os.path.basename(path)):
                self._row(p, path, sc_io._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN.format(path=path),
                          cls=ValueError)
        nonso = os.path.join(d, "bond_archive.npz")
        np.savez(nonso, bond_archive_schema=np.int64(3))
        self._row(p, nonso, sc_io._SO_BOND_GREEN_NONSO_ARCHIVE.format(path=nonso), cls=ValueError)

    def test_flex_output_file(self):
        """sigma.npz, chi0q.npz or chiq*.npz named as bond_green: the hint
        names the file's kind instead of the structure step's "malformed"
        text; a file that also has a green member goes on to the checks."""
        p = _base("ZEEMAN_SZ_ALL")
        sigma = os.path.join(os.path.dirname(_GREENS["ZEEMAN_SZ_ALL"]), "sigma.npz")
        self.assertTrue(os.path.isfile(sigma))
        self._row(p, sigma, sc_io._SO_BOND_GREEN_IS_FLEX_OUTPUT.format(path=sigma, kind="sigma"),
                  cls=ValueError)
        d = _scratch(self)
        for kind in ("chi0q", "chiq", "chiq_s", "chiq_c"):
            with self.subTest(kind=kind):
                path = os.path.join(d, kind + ".npz")
                np.savez(path, **{kind: np.zeros(2)}, spin_orbital=True)
                self._row(p, path, sc_io._SO_BOND_GREEN_IS_FLEX_OUTPUT.format(path=path,
                                                                            kind=kind),
                          cls=ValueError)
        # a green member present: no hint (the green archive goes on to its checks)
        path = _edited(self, _GREENS["ZEEMAN_SZ_ALL"], sigma=np.zeros(2))
        self.assertEqual(_borrow(p, path).carrier.source, "external_npz")

    def test_onsite_green(self):
        ctx = T._context()
        run = T._run(ctx, bond=sc_so._so_bond_run(ctx, CELL, None))
        path = T._onsite_green_archive(self, run)
        with _PayloadReads() as spy:
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                sc_io._so_borrowed_green(path, run, norb_phys=1, nmat=T._RUN_NMAT,
                                      shape=(2, 2, 1), conserving_input=True)
        e = cm.exception
        self.assertEqual(str(e), "{} was produced by test FLEX producer (so_onsite), but this "
                         "run needs a so_bond artifact. ".format(path)
                         + sc_io._SO_BOND_GREEN_ONSITE_REMEDY)
        self.assertEqual(e.check, "capability")
        self.assertIsInstance(e.__cause__, fingerprint.SOArchiveError)
        self.assertEqual(spy.reads, [])

    def test_consumer_differs(self):
        fx = X.ZEEMAN_SZ_ALL
        path = _GREENS["ZEEMAN_SZ_ALL"]
        head = "'{}' (written by {}) does not match this run: ".format(path, _LABEL)
        again = ("Archives are reused only at the identical state point (not a warm start "
                 "across a sweep); re-run {} with this input".format(_LABEL))
        hop11 = fx._replace(hopping=tuple((R, t * 1.1) for R, t in fx.hopping))
        intra = fx._replace(interactions=tuple(
            (t, (((0, 0, 0), (0, 0), 0.9),)) if t == "CoulombIntra" else (t, r)
            for t, r in fx.interactions))
        zero2 = (((2, 0, 0), (0, 0), 0.0), ((-2, 0, 0), (0, 0), 0.0))
        declared = fx._replace(interactions=tuple(
            (t, r + zero2) if t == "CoulombInter" else (t, r) for t, r in fx.interactions))
        state = "Archives are reused only at the identical state point; "
        rerun = "re-run {} at this filling and temperature.".format(_LABEL)
        rows = (
            ("Transfer x1.1", _consumer(hop11), head + "band differs: Transfer: file 16 rows, "
             "run 16 rows; first differing row (R=(-1, -1, 0), a=1, b=1): file (0.25+0j), run "
             "(0.275+0j). " + again + " (the same Transfer/Extern/coeff_extern).", "band"),
            ("interaction", _consumer(intra), head + "interaction differs: file 53 records, "
             "run 53 records; first differing record CoulombIntra modes (0, 0, 1, 1) cells "
             "((0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)): file (0.83+0j), run (0.9+0j). "
             + again + " (the same interaction definitions).", "interaction"),
            ("declared-zero second shell", _consumer(declared), head + "declared topology "
             "differs: file 5 displacements, run 7 displacements. " + again
             + " (the same interaction definitions).", "topology"),
            ("Nmat 16", _consumer(fx, param={"Nmat": 16}), head + "state differs: nmat 32 "
             "(file) vs 16 (run). " + state + "set [mode.param] Nmat = 32 (the producer's "
             "grid).", "state"),
            ("T 0.25", _consumer(fx, param={"T": 0.25}), head + "state differs: beta 2.0 "
             "(file) vs 4.0 (run). " + state + rerun, "state"),
            # the producer's density as recorded (its last digit is round-off)
            ("filling 0.35", _consumer(fx, param={"filling": 0.35}), head + "state differs: "
             "target 0.8 (file) vs 0.7 (run) electrons per cell; density {!r} (file) vs "
             "target 0.7 (run, tol 1e-08) electrons per cell. ".format(
                 float(_members(path)["state_density"])) + state + rerun, "state"),
        )
        for label, p, text, check in rows:
            with self.subTest(row=label):
                self.assertEqual(self._row(p, path, text).check, check)
        # the residual step comes first: a file that did not reach its
        # filling, read by the Transfer x1.1 consumer
        res = _edited(self, path, state_mu_tol=np.float64(1e-6),
                      state_density=np.float64(float(_members(path)["state_target"]) + 1e-6))
        e = self._row(rows[0][1], res, "'{}' cannot be reused: the producer did not reach its "
                      "filling (residual 1.000e-06 electrons per cell > 1e-8·max(1,|t|) = "
                      "1.000e-08). Regenerate it with {} so that its density closes to that "
                      "tolerance.".format(res, _LABEL))
        self.assertEqual(e.check, "residual")
        for shells in (1, 3):
            with self.subTest(bond_max_shells=shells):
                g = _borrow(_consumer(fx, eliashberg={"Bond_Max_Shells": shells}), path)
                self.assertEqual(g.carrier.source, "external_npz")

    def test_markers(self):
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        for member, value in (("momentum_convention", _DROP), ("sc_vertex_version", _DROP),
                              ("sc_vertex_version", np.int64(2)), ("wavevector_index", _DROP)):
            with self.subTest(member=member, value=value):
                path = _edited(self, src, **{member: value})
                e = self._row(p, path, sc_io.T_BG_MARKER.format(path=path, member=member,
                                                             label=_LABEL))
                self.assertIsNone(e.check)
        index = _members(src)["wavevector_index"]
        path = _edited(self, src, wavevector_index=index[::-1].copy())
        self._row(p, path, "'{}' is malformed or truncated: wavevector_index is not the k grid "
                  "[4, 4, 1] of this run (CellShape / SubShape). Regenerate it with {}."
                  .format(path, _LABEL))

    def test_scf_converged_marker(self):
        """scf_converged must be a 0-d boolean when present (checked before
        the payload read); absent, the record has None."""
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        bad = (("string", np.str_("true")), ("int", np.int64(1)),
               ("1-d singleton", np.array([True])), ("empty", np.array([], dtype=bool)),
               ("multi-element", np.array([True, True])))
        for label, value in bad:
            with self.subTest(value=label):
                path = _edited(self, src, scf_converged=value)
                e = self._row(p, path, sc_io.T_BG_SCF.format(path=path, label=_LABEL))
                self.assertIsNone(e.check)
        path = _edited(self, src, scf_converged=_DROP)
        with _records() as rec:
            g = _borrow(p, path)
        self.assertIsNone(g.scf_converged)
        self.assertEqual(_warnings(rec), [])

    def test_unconverged(self):
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        path = _edited(self, src, scf_converged=np.bool_(False))
        with _records() as rec:
            g = _borrow(p, path)
        self.assertIs(g.scf_converged, False)
        self.assertEqual(_warnings(rec), [sc_bond._SO_BOND_GREEN_UNCONVERGED % path])
        # an edited green on top: the digest refuses it, no unconverged WARNING
        edited = _edited(self, _perturbed(self, src, 1e-6, restamp=False),
                         scf_converged=np.bool_(False))
        self._row(p, edited, _mismatch(edited), read=True)

    def test_payload(self):
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        G = _members(src)["green"]
        short = G[:, :30].copy()
        path = _edited(self, src, green=short, fp_green=_restamped(short))
        self._row(p, path, sc_io.T_BG_SHAPE.format(path=path, stored=(1, 30, 16, 2, 2),
                                                expected=(1, 32, 16, 2, 2), label=_LABEL))
        nan = G.copy()
        nan[:, 3, 0, 0, 0] = np.nan
        path = _edited(self, src, green=nan)
        with mock.patch.object(fingerprint, "green_digest",
                               wraps=fingerprint.green_digest) as digest:
            e = self._row(p, path, "'{}' is malformed or truncated: green has non-finite "
                          "values at frequency index 3. Regenerate it with {}.".format(
                              path, _LABEL), read=True)
        self.assertIsNone(e.check)
        self.assertEqual(digest.call_count, 0)
        path = _perturbed(self, src, 1e-6, restamp=False)
        e = self._row(p, path, _mismatch(path), read=True)
        self.assertEqual(e.check, "green")


# ------------------------------------------------------------ tail policy
class TestTailPolicy(unittest.TestCase):
    """No tail is applied to a Green function read from a file; the
    default bond_coeff_tail is noted at INFO, a configured non-zero one is
    a WARNING."""

    def test_tail_policy(self):
        self.assertEqual(sc_bond._SO_BOND_TAIL_NOT_APPLIED,
                         "[eliashberg] bond_coeff_tail = %r is not applied with bond_green = "
                         "'%s': the tail correction is implemented for the band Green function "
                         "only (bond_coeff_tail_applied = 0), so the bond bubble carries the "
                         "Matsubara truncation error, which decreases only about as 1/Nmat; "
                         "check the eigenvalue against Nmat.")
        self.assertEqual(sc_bond._SO_BOND_TAIL_DEFAULT_INFO,
                         "[eliashberg] bond_coeff_tail: the default 1.0 applies to the band "
                         "Green function only; the Green function read from bond_green is "
                         "used without a tail correction.")
        path = _GREENS["ZEEMAN_SZ_ALL"]
        rows = ((None, 1, None), (0.5, 0, "0.5"), (1.0, 0, "1.0"), (0.0, 0, None))
        for value, default_info, shown in rows:
            with self.subTest(bond_coeff_tail=value):
                eli = None if value is None else {"Bond_Coeff_Tail": value}
                p = _consumer(X.ZEEMAN_SZ_ALL, eliashberg=eli)
                with _records() as rec:
                    g = _borrow(p, path)
                self.assertEqual(g.carrier.coeff_tail_applied, 0.0)
                self.assertEqual(_count(rec, sc_bond._SO_BOND_TAIL_DEFAULT_INFO), default_info)
                tail = [w for w in _warnings(rec) if "bond_coeff_tail" in w]
                if shown is None:
                    self.assertEqual(tail, [])
                else:
                    self.assertEqual(tail, [
                        "[eliashberg] bond_coeff_tail = {} is not applied with bond_green = "
                        "'{}': the tail correction is implemented for the band Green "
                        "function only (bond_coeff_tail_applied = 0), so the bond bubble "
                        "carries the Matsubara truncation error, which decreases only about "
                        "as 1/Nmat; check the eigenvalue against Nmat.".format(shown, path)])


# ------------------------------------------------------------ S_z leakage (C1b-PG)
#: so_bond.project_borrowed_green's WARNING and refusal (``%s`` path, ``%.3e`` leak)
_PROJECTED_WARNING = (
    "the Green function '%s' has S_z-off-sector entries at the round-off level (leak = max "
    "off-sector |G| / max |G| = %.3e <= 1e-12: numerical leakage from a non-canonical "
    "producer); they are projected to exact zeros")
_LEAK_REFUSAL = (
    "the Green function '%s' breaks S_z (leak = max off-sector |G| / max |G| = %.3e > "
    "1e-12) although this run's input conserves S_z; it was produced from a different "
    "Hamiltonian. Regenerate it from this run's input.")
_BUBBLE_REFUSAL = ("raised ValueError: so_bond_bubble: the static bond bubble has "
                   "S_z-off-sector entries")
_GATE_TOL = 1e-12
_NON_VACUITY = 1e-3


def rel_or_abs(a, b):
    """``max|a - b| / max|b|``, or ``max|a - b|`` when ``max|b| = 0``."""
    a, b = np.asarray(a), np.asarray(b)
    num = float(np.max(np.abs(a - b)))
    den = float(np.max(np.abs(b)))
    return num / den if den != 0.0 else num


def _gate_failures(gate):
    """The failed checks of one gate: ``gate(check)`` runs the production
    call and reports each assertion as ``check(name, ok, detail="")``. A
    production exception is the failure ``"raised <Type>: <text>"``. A gate
    row asserts ``== []``, a break row that the named failure is listed."""
    failures = []

    def check(name, ok, detail=""):
        if not ok:
            failures.append(name + (" (" + detail + ")" if detail else ""))
    try:
        gate(check)
    except Exception as e:  # a production exception fails the gate
        failures.append("raised {}: {}".format(type(e).__name__, e))
    return failures


def _bounded(check, name, value, bound):
    """``value`` (a metric) finite, then ``<= bound``."""
    finite = bool(np.isfinite(value))
    check(name + " finite", finite, repr(value))
    check(name, finite and value <= bound, "%.3e > %.0e" % (value, bound))


def _has(failures, name):
    return any(f.startswith(name) for f in failures)


def _leak(G):
    """``max off-sector |G| / max |G|`` of a canonical ``(..., nd, nd)`` G,
    computed here (not by ``so_bond.green_leak``)."""
    n = G.shape[-1] // 2
    off = max(float(np.max(np.abs(G[..., :n, n:]))), float(np.max(np.abs(G[..., n:, :n]))))
    return off / float(np.max(np.abs(G)))


def _warning_records(records):
    return [(r.name, r.getMessage()) for r in records if r.levelno >= logging.WARNING]


def _clean_operator(p):
    """The production operator of the run ``p`` on the clean module green of
    ZEEMAN_SZ_ALL (the reference of gate (ii), built before any patch)."""
    g = _borrow(p, _GREENS["ZEEMAN_SZ_ALL"])
    return T._operator(p, g, T._admit(p))


def _projected_gate(p, path, *, layout, ref, seen):
    """Gate C1b-PG (ii) on the run ``p``: the re-stamped file at ``path``
    (``G + 1e-14 max|G| P``) is projected with one WARNING; the record holds
    the clean digest and the file's; the operator equals ``ref``."""
    G2 = _members(path)["green"]
    fp_clean = str(_members(_GREENS["ZEEMAN_SZ_ALL"])["fp_green"])
    fp_file = str(_members(path)["fp_green"])
    leak = _leak(G2)

    def gate(check):
        with _records() as rec:
            g = _borrow(p, path)
            adm = T._admit(p)
            static = T._operator(p, g, adm)
        check("layout " + layout, adm.layout == layout and static.layout == layout,
              adm.layout)
        warns = _warning_records(rec)
        check("one projection WARNING",
              warns == [("hwave.solver.so_bond", _PROJECTED_WARNING % (path, leak))], repr(warns))
        check("accepted INFO", _count(rec, sc_bond._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
                                      "%.3e" % leak, "yes") == 1)
        check("fp_green == clean digest", g.fp_green == fp_clean, str(g.fp_green))
        check("fp_green_source == re-stamped digest", g.fp_green_source == fp_file,
              str(g.fp_green_source))
        n = p.norb_phys
        kw = g.carrier.deflated_kw
        check("off-sector of the used G exactly zero",
              not np.any(kw[..., :n, n:]) and not np.any(kw[..., n:, :n]))
        check("digest of the used G == clean digest",
              fingerprint.green_digest(kw, fingerprint.GREEN_LAYOUT) == fp_clean)
        miss = T._miss(lambda v: T._apply(static, v), lambda v: T._apply(ref, v),
                       T._gaps(p.fixture))
        _bounded(check, "operator miss vs the clean file", miss, _GATE_TOL)
    return _gate_failures(gate)


def _refused_gate(p, path, *, seen):
    """Gate C1b-PG (iii): the re-stamped file at ``path`` (``1e-6``) is
    refused by the leak check, before any bubble."""
    leak = _leak(_members(path)["green"])
    seen["leak"] = leak

    def gate(check):
        refused = None
        with mock.patch.object(so_bond, "so_bond_bubble",
                               wraps=so_bond.so_bond_bubble) as bubble:
            try:
                g = _borrow(p, path)
                T._operator(p, g, T._admit(p))
            except ValueError as e:
                refused = e
        check("refused with the leak text", type(refused) is ValueError
              and str(refused) == _LEAK_REFUSAL % (path, leak), repr(refused))
        check("so_bond_bubble not called", bubble.call_count == 0, str(bubble.call_count))
    return _gate_failures(gate)


def _digest_gate(p, path, *, seen):
    """Gate C1b-PG (iv): the edited, not re-stamped file at ``path`` is
    refused by the digest step, before the leak is measured."""
    text = _mismatch(path)

    def gate(check):
        refused = None
        with _records() as rec, \
                mock.patch.object(so_bond, "project_borrowed_green",
                                  wraps=so_bond.project_borrowed_green) as project, \
                mock.patch.object(so_bond, "green_leak", wraps=so_bond.green_leak) as leak:
            try:
                _borrow(p, path)
            except fingerprint.SOArchiveError as e:
                refused = e
        seen["warnings"] = _warning_records(rec)
        check("refused with the fp_green mismatch text", refused is not None
              and str(refused) == text and refused.check == "green", repr(refused))
        check("project_borrowed_green not called", project.call_count == 0,
              str(project.call_count))
        check("green_leak not called", leak.call_count == 0, str(leak.call_count))
    return _gate_failures(gate)


def _conserving_rule_gate(p, path):
    """The conserving-input rule on ZPL: the band conserves S_z, the interaction does not, so
    the re-stamped file at ``path`` (``1e-6``) is used as stored: no leak
    check, no projection, no WARNING."""
    m = _members(path)

    def gate(check):
        with _records() as rec:
            g = _borrow(p, path)
        check("band conserves S_z", g.conserving is True)
        check("leak not measured", g.leak is None, repr(g.leak))
        check("fp_green == the file's", g.fp_green == str(m["fp_green"]), str(g.fp_green))
        check("no fp_green_source", g.fp_green_source is None, repr(g.fp_green_source))
        check("G used as stored", np.array_equal(g.carrier.deflated_kw, m["green"]))
        check("no WARNING", _warning_records(rec) == [], repr(_warning_records(rec)))
        check("accepted INFO", _count(rec, sc_bond._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
                                      sc_bond._SO_BOND_GREEN_LEAK_UNCHECKED, "no") == 1)
    return _gate_failures(gate)


def _identity_projection(G, *, norb_phys, source, axes=(-2, -1)):
    """Break (b1): no projection; the leak is still measured."""
    return G, so_bond.green_leak(G, norb_phys, axes=axes)


def _digest_unchecked(stored, path, *, label, green_fp, green_only=False):
    """Break (b2): the stored digest is returned without comparing."""
    return stored


_DENSE = {"SO_Bond_Layout": "dense"}


class TestBorrowedLeak(unittest.TestCase):
    """Gates C1b-PG (ii)-(iv) and the conserving-input rule: a borrowed
    Green function with S_z leakage up to 1e-12 is projected (WARNING, both
    digests recorded, the operator of the clean file), above that refused
    before any bubble; an edited file is refused by its digest before the
    leak is measured; an input that does not conserve S_z uses the file as
    stored."""

    maxDiff = None

    def test_sub_threshold_projected(self):
        fx = X.ZEEMAN_SZ_ALL
        clean = _GREENS["ZEEMAN_SZ_ALL"]
        path = _perturbed(self, clean, 1e-14, restamp=True)
        G, G2 = _members(clean)["green"], _members(path)["green"]
        # non-vacuity: the test-side dense bubble of G' has off-sector entries
        p = _base("ZEEMAN_SZ_ALL")
        nd, nvol = 2 * p.norb_phys, int(np.prod(p.shape))
        geom = so_bond.slot_geometry(p.so_run.bond.topo, p.shape)
        from tests.so_bond_fixtures import NOCAP
        plan = so_bond.plan_so_bond(nmat=p.nmat, nvol=nvol, B=geom.B, nd=nd,
                                    consumers=("pairing_static",), caps=NOCAP, layout="dense")
        ND = geom.B * nd * nd
        Xb = so_bond.so_bond_bubble(G2, None, p.beta, geom, p.shape, plan=plan, freq="static",
                                    sink=np.zeros((nvol, ND, ND), dtype=complex))
        delta = float(np.max(np.abs(G2 - G)))
        ratio = (float(np.max(np.abs(Xb[:, T._bond_off_sector_mask(nd, geom.B)])))
                 / (delta * float(np.max(np.abs(G)))))
        self.assertTrue(np.isfinite(ratio))
        self.assertGreater(ratio, _NON_VACUITY)
        for layout, eli in (("sz_packed", None), ("dense", _DENSE)):
            with self.subTest(layout=layout):
                q = _consumer(fx, eliashberg=eli)
                seen = {}
                self.assertEqual(_projected_gate(q, path, layout=layout, ref=_clean_operator(q),
                                                 seen=seen), [])
        # through the entry: the gap archive records both digests
        ref = F.run_case(fx, eliashberg={"Bond_Green": clean}, fresh=True).gap_npz
        fp_clean, fp_file = str(_members(clean)["fp_green"]), str(_members(path)["fp_green"])

        def entry(check):
            r = F.run_case(fx, eliashberg={"Bond_Green": path}, fresh=True)
            gap = r.gap_npz
            check("gap fp_green == clean digest", str(gap["fp_green"]) == fp_clean)
            check("gap fp_green_source == re-stamped digest",
                  "fp_green_source" in gap and str(gap["fp_green_source"]) == fp_file)
            check("eigenvalue.dat leak", T._dat_keys(r.eigenvalue_dat)["bond_green_leak"]
                  == "%.3e" % _leak(G2))
            _bounded(check, "eigenvalues rel_or_abs", rel_or_abs(gap["eigenvalues"],
                                                                 ref["eigenvalues"]), _GATE_TOL)
        self.assertNotIn("fp_green_source", ref)
        self.assertEqual(_gate_failures(entry), [])

    def test_above_threshold_refused(self):
        path = _perturbed(self, _GREENS["ZEEMAN_SZ_ALL"], 1e-6, restamp=True)
        seen = {}
        self.assertEqual(_refused_gate(_base("ZEEMAN_SZ_ALL"), path, seen=seen), [])
        self.assertGreater(seen["leak"], so_bond._LEAK_NUMERICAL)

    def test_digest_before_leak(self):
        path = _perturbed(self, _GREENS["ZEEMAN_SZ_ALL"], 1e-14, restamp=False)
        self.assertEqual(_digest_gate(_base("ZEEMAN_SZ_ALL"), path, seen={}), [])

    def test_conserving_input_rule(self):
        path = _perturbed(self, _GREENS["ZPL"], 1e-6, restamp=True)
        p = _base("ZPL")
        # non-vacuity: a conserving rule would refuse this file
        self.assertGreater(_leak(_members(path)["green"]), so_bond._LEAK_NUMERICAL)
        self.assertIs(sc_bond._so_conserving_input(p.so_run, p.norb_phys), False)
        self.assertEqual(_conserving_rule_gate(p, path), [])


class TestBorrowedLeakBreaks(unittest.TestCase):
    """Deliberate breaks of the S_z-leakage gates, made with call-time
    patches: each row asserts that the break reached the production call,
    then that its gate fails with the named failure."""

    def test_b1_no_projection_packed(self):
        path = _perturbed(self, _GREENS["ZEEMAN_SZ_ALL"], 1e-14, restamp=True)
        p = _base("ZEEMAN_SZ_ALL")
        ref = _clean_operator(p)
        with mock.patch.object(so_bond, "project_borrowed_green",
                               side_effect=_identity_projection) as brk:
            failures = _projected_gate(p, path, layout="sz_packed", ref=ref, seen={})
        self.assertEqual(brk.call_count, 1)
        self.assertTrue(_has(failures, _BUBBLE_REFUSAL), failures)

    def test_b1_no_projection_dense(self):
        path = _perturbed(self, _GREENS["ZEEMAN_SZ_ALL"], 1e-14, restamp=True)
        p = _consumer(X.ZEEMAN_SZ_ALL, eliashberg=_DENSE)
        ref = _clean_operator(p)
        with mock.patch.object(so_bond, "project_borrowed_green",
                               side_effect=_identity_projection) as brk:
            failures = _projected_gate(p, path, layout="dense", ref=ref, seen={})
        self.assertEqual(brk.call_count, 1)
        self.assertTrue(_has(failures, "fp_green == clean digest"), failures)
        # the G handed to the kernel is the unprojected file content
        self.assertTrue(_has(failures, "off-sector of the used G exactly zero"), failures)
        self.assertTrue(_has(failures, "digest of the used G == clean digest"), failures)

    def test_b2_digest_unchecked(self):
        path = _perturbed(self, _GREENS["ZEEMAN_SZ_ALL"], 1e-14, restamp=False)
        seen = {}
        with mock.patch.object(fingerprint, "check_green_digest",
                               side_effect=_digest_unchecked) as brk:
            failures = _digest_gate(_base("ZEEMAN_SZ_ALL"), path, seen=seen)
        self.assertEqual(brk.call_count, 1)
        self.assertTrue(_has(failures, "refused with the fp_green mismatch text"), failures)
        self.assertTrue(_has(failures, "project_borrowed_green not called"), failures)
        # accepted, with the projection WARNING
        self.assertEqual(seen["warnings"], [(
            "hwave.solver.so_bond",
            _PROJECTED_WARNING % (path, _leak(_members(path)["green"])))])

    def test_b3_conserving_rule_forced(self):
        path = _perturbed(self, _GREENS["ZPL"], 1e-6, restamp=True)
        p = _base("ZPL")
        with mock.patch.object(sc_bond, "_so_conserving_input", return_value=True) as brk:
            failures = _conserving_rule_gate(p, path)
        self.assertEqual(brk.call_count, 1)
        self.assertTrue(_has(failures, "raised ValueError: " + _LEAK_REFUSAL % (
            path, _leak(_members(path)["green"]))), failures)


# ------------------------------------------------------------ content and axes (C1b-BG)
_BREAK_MISS = 1e-3
#: the content rows: (producer, layout, [eliashberg] overrides)
_CONTENT_ROWS = (("ZEEMAN_SZ_ALL", "sz_packed", None), ("ZEEMAN_SZ_ALL", "dense", _DENSE),
                 ("SOC_ALL_TYPES", "dense", None))
_KERNELS = {}


def _sc_layout(G, shape):
    """A canonical ``(nmat, nvol, nd, nd)`` G in the ``(nd, nd, Nx, Ny, Nz,
    nmat)`` layout, reshaped here (not by ``sc._green_canonical_to_sc``)."""
    nmat, nd = G.shape[0], G.shape[-1]
    return G.reshape((nmat,) + tuple(shape) + (nd, nd)).transpose(4, 5, 1, 2, 3, 0)


def _file_kernel(name, p):
    """The reference kernel of the module green of ``name``: the validated
    array functions on the file's G (reshaped here), with no tail (a file
    carries none) and the run's ``g2_tail``. Built before any patch; one per
    producer (the reference kernel is always dense)."""
    if name not in _KERNELS:
        G = _sc_layout(_members(_GREENS[name])["green"][0], p.shape)
        v_bare, v_fluct, geom = _reference_vertex(p, G, coeff_tail=0.0)
        _KERNELS[name] = _reference_kernel_of(p, v_bare + v_fluct, G, geom,
                                              g2_tail=p.cfg.g2_tail)
    return _KERNELS[name]


def _content_gate(p, path, *, K, layout, seen):
    """Gate C1b-BG content: the production operator of the run ``p`` on the
    borrowed green at ``path`` equals the reference kernel ``K``."""
    def gate(check):
        g = _borrow(p, path)
        adm = T._admit(p)
        static = T._operator(p, g, adm)
        check("layout " + layout, static.layout == layout, static.layout)
        miss = _miss(lambda v: T._apply(static, v), K.matvec, _gaps(p.fixture))
        seen["miss"] = miss
        _bounded(check, "operator miss vs the reference kernel", miss, _GATE_TOL)
    return _gate_failures(gate)


def _axis_gate(name, p, *, seen):
    """Gate C1b-BG axis oracle: the carrier ``full_sc`` of the borrowed green
    of ``name`` equals ``G_ind(k, i w_n) = [(i w_n + mu) 1 - eps_k -
    Sigma(k, i w_n)]^-1`` with ``w_n = (2n + 1 - nmat) pi / beta``, hwave_sc's
    own ``eps_k``, the file's ``state_mu`` and the producer's sigma.npz."""
    path = _GREENS[name]
    green = _members(path)
    sigma = _members(os.path.join(os.path.dirname(path), "sigma.npz"))
    nd, nmat, shape = 2 * p.norb_phys, p.nmat, tuple(p.shape)
    nvol = int(np.prod(shape))
    n = np.arange(nmat)
    iw = 1j * (2.0 * n + 1.0 - nmat) * np.pi / p.beta
    h0 = np.moveaxis(np.asarray(p.epsilon_k).reshape(nd, nd, nvol), -1, 0)
    mu = float(green["state_mu"])
    S = sigma["sigma"]

    def gate(check):
        check("sigma.npz fp_green == green.npz fp_green",
              str(sigma["fp_green"]) == str(green["fp_green"]),
              "%s != %s" % (sigma["fp_green"], green["fp_green"]))
        check("sigma shape (1, nmat, nvol, nd, nd)", S.shape == (1, nmat, nvol, nd, nd),
              repr(S.shape))
        G_ind = np.linalg.inv((iw[:, None, None, None] + mu) * np.eye(nd)
                              - h0[None] - S[0])
        g = _borrow(p, path)
        err = rel_or_abs(g.carrier.full_sc, _sc_layout(G_ind, shape))
        seen["err"] = err
        _bounded(check, "full_sc vs G_ind", err, _GATE_TOL)
    return _gate_failures(gate)


def _scaled_bubble(real):
    """Break: the bond bubble's first argument (the Green function) x1.01."""
    def bubble(G, *args, **kw):
        return real(1.01 * G, *args, **kw)
    return bubble


def _scaled_to_sc(real):
    """Break: the carrier ``full_sc`` (the pair propagator's G) x1.01."""
    def to_sc(*args, **kw):
        return 1.01 * real(*args, **kw)
    return to_sc


class TestBorrowedContent(unittest.TestCase):
    """Gate C1b-BG content and axis oracle: the operator built on a borrowed
    Green function is the reference kernel of the file's G (and not the
    band's), and the file's G is the Dyson inverse of hwave_sc's band and the
    producer's self-energy on hwave_sc's k and frequency axes. Each break
    scales one input of the operator by 1.01 and must miss by >= 1e-3."""

    def _row(self, name, eli):
        return _consumer(_PRODUCERS[name], eliashberg=eli)

    def test_content(self):
        for name, layout, eli in _CONTENT_ROWS:
            with self.subTest(producer=name, layout=layout):
                p = self._row(name, eli)
                K = _file_kernel(name, p)
                seen = {}
                self.assertEqual(_content_gate(p, _GREENS[name], K=K, layout=layout,
                                               seen=seen), [])
                # non-vacuity: the band's operator is not the file's
                band = T._operator(p, T._green_of(p), T._admit(p))
                nv = _miss(lambda v: T._apply(band, v), K.matvec, _gaps(p.fixture))
                self.assertTrue(np.isfinite(nv))
                self.assertGreaterEqual(nv, _BREAK_MISS)

    def _break(self, target, attr, wrap):
        for name, layout, eli in _CONTENT_ROWS:
            with self.subTest(producer=name, layout=layout):
                p = self._row(name, eli)
                K = _file_kernel(name, p)
                seen = {}
                real = getattr(target, attr)
                with mock.patch.object(target, attr, side_effect=wrap(real)) as brk:
                    failures = _content_gate(p, _GREENS[name], K=K, layout=layout, seen=seen)
                self.assertEqual(brk.call_count, 1)
                self.assertTrue(_has(failures, "operator miss vs the reference kernel"),
                                failures)
                self.assertGreaterEqual(seen["miss"], _BREAK_MISS)

    def test_break_bubble_input(self):
        self._break(so_bond, "so_bond_bubble", _scaled_bubble)

    def test_break_pair_propagator(self):
        self._break(sc_bond, "_green_canonical_to_sc", _scaled_to_sc)

    def test_axis_oracle(self):
        for name in _SIGMA:
            with self.subTest(producer=name):
                seen = {}
                self.assertEqual(_axis_gate(name, _base(name), seen=seen), [])


# ------------------------------------------------------------ memory (C1b-Mem)
#: the fixed part and the member factor of the resident row of a Green
#: function read from a file (alive in every phase: the read leaves
#: reference-cycle garbage until the next cyclic collection)
_FILE_READ_FIXED = 262144
_FILE_READ_FACTOR = 4


def _member_bytes(path):
    """The bytes of the members other than ``green`` of the npz ``path``,
    each the larger of its zip record size and the bytes of the array
    ``np.load`` materialises (independent of the production header read)."""
    with zipfile.ZipFile(path) as zf:
        sizes = {i.filename[:-4]: i.file_size for i in zf.infolist()}
    with np.load(path, allow_pickle=False) as z:
        return sum(max(sizes[k], z[k].nbytes) for k in z.files if k != "green")


def _file_read_row(path):
    """The resident ``file_read`` row the admission of ``path`` carries."""
    return ("resident", "file_read", _FILE_READ_FIXED + _FILE_READ_FACTOR * _member_bytes(path))


def _fixed_only(*args, **kw):
    """Break: the ``file_read`` row without its member part (the fixed
    256 KiB alone)."""
    return (("file_read", _FILE_READ_FIXED),)


def _inflated(test, path, reps):
    """A copy of ``path`` whose provenance interaction records
    (``prov_rec_*``) are repeated ``reps`` times: a large provenance table
    that the checks still accept (the digests are not recomputed from these
    members)."""
    m = _members(path)
    names = ("prov_rec_types", "prov_rec_modes", "prov_rec_cells", "prov_rec_v")
    return _edited(test, path, **{k: np.concatenate([m[k]] * reps, axis=0) for k in names})


@contextlib.contextmanager
def _gc_held():
    """The cyclic garbage collector held off, after one full collection: the
    traced phases then keep their reference-cycle garbage alive throughout,
    the worst case for every row, independent of when a collection would
    otherwise run."""
    gc.collect()
    enabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if enabled:
            gc.enable()


def _borrowed_green_fn(p):
    """The Green function record of the run ``p`` from its own
    ``[eliashberg] bond_green``."""
    return sc_bond._so_bond_green(p.epsilon_k, p.so_run, norb_phys=p.norb_phys, filling=p.filling,
                             beta=p.beta, nmat=p.nmat, coeff_tail=p.cfg.coeff_tail,
                             bond_green=p.cfg.green,
                             coeff_tail_source=p.cfg.coeff_tail_source)


def _memory_consumer(path, nmat, eli=None):
    """The admission run of ZEEMAN_SZ_ALL at T 0.5, ``nmat``, reading
    ``path`` through ``Bond_Green`` (``eli`` on top)."""
    return T._prepare_admission(X.ZEEMAN_SZ_ALL, eliashberg=dict({"Bond_Green": path},
                                                                  **(eli or {})),
                                param={"T": 0.5, "Nmat": nmat})


def _warm(p):
    """One untraced run of the borrowed chain of ``p`` (the reader, the
    operator and the eigenvalue solve), so that no first-call import or
    cache lands inside a traced phase."""
    g = _borrowed_green_fn(p)
    static = T._operator(p, g, T._admit(p))
    T._solve_leg(p, static, "eigenvalue")


class _PhaseChecks:
    """The assertions ``T.check_phase_rows`` makes, reported to a gate's
    ``check`` instead of raised; the traced peak and the rows of each phase
    are kept in ``values``."""

    def __init__(self, check):
        self.check = check
        self.phase = None
        self.values = {}

    @contextlib.contextmanager
    def subTest(self, **kw):
        self.phase = kw.get("phase")
        try:
            yield
        finally:
            self.phase = None

    def assertEqual(self, a, b, msg=None):
        self.check("equal", a == b, "%r != %r" % (a, b))

    def assertNotEqual(self, a, b, msg=None):
        self.check("not equal", a != b, "%r == %r" % (a, b))

    def assertLessEqual(self, a, b, msg=None):
        self.values[self.phase] = (a, b)
        self.check("%s traced <= rows" % self.phase, a <= b, "%d > %d" % (a, b))


def _memory_gate(p, *, seen):
    """Gate C1b-Mem on a borrowed Green function: every traced phase of the
    run ``p`` (eigenvalue solve) within its rows, through
    ``T.check_phase_rows`` with the borrowed ``green_fn``."""
    def gate(check):
        rec = _PhaseChecks(check)
        seen["values"] = rec.values
        with _gc_held():
            seen["adm"] = T.check_phase_rows(rec, p, "eigenvalue", green_fn=_borrowed_green_fn)
    return _gate_failures(gate)


def _dropping(real, names):
    """Break: ``static_so_bond_rows`` without the ``green_static`` rows
    ``names``."""
    def rows(*args, **kw):
        return tuple(r for r in real(*args, **kw)
                     if not (r[0] == "green_static" and r[1] in names))
    return rows


def _no_resident(*args, **kw):
    """Break: ``static_so_bond_resident_rows`` without the ``file_read``
    row."""
    return ()


class TestBorrowedAdmission(unittest.TestCase):
    """Gate C1b-Mem of a borrowed Green function: the admission of a run with
    ``bond_green`` carries the resident ``file_read`` row (a band run does
    not), sized from the file's member headers, and every traced phase, the
    read included, stays within its rows (the cyclic collector held off,
    :func:`_gc_held`): packed (projected file) and dense (clean file) at
    Nmat 32, packed (projected file) at Nmat 8, where the band rows alone
    are below the read peak, and a file with a large provenance table.
    Breaks: (m1) without ``carrier_build``, ``transient`` and ``file_read``
    every row fails; (m2) without ``file_read`` the Nmat 8 row fails; (m3)
    the large-provenance row fails with the fixed 256 KiB alone. A cap
    between the fixed and the sized estimate refuses the large-provenance
    run at admission, before any member of the file is read."""

    def _rows(self, labels=None):
        """``(label, p, layout)`` of the memory rows (all, or those named in
        ``labels``), each warmed up."""
        def packed(name, nmat):
            return _memory_consumer(_perturbed(self, _GREENS[name], 1e-14, restamp=True), nmat)
        make = (("packed projected Nmat 32", lambda: packed("ZEEMAN_SZ_ALL", 32), "sz_packed"),
                ("dense clean Nmat 32",
                 lambda: _memory_consumer(_GREENS["ZEEMAN_SZ_ALL"], 32, _DENSE), "dense"),
                ("packed projected Nmat 8", lambda: packed(_N8, 8), "sz_packed"))
        rows = [(label, build(), layout) for label, build, layout in make
                if labels is None or label in labels]
        for _, p, _ in rows:
            _warm(p)
        return rows

    def test_file_read_row(self):
        for nmat in (32, 8):
            with self.subTest(nmat=nmat):
                path = _GREENS["ZEEMAN_SZ_ALL" if nmat == 32 else _N8]
                rows = T._admit(_memory_consumer(path, nmat)).plan.rows
                row = _file_read_row(path)
                self.assertGreater(row[2], _FILE_READ_FIXED)
                self.assertIn(row, rows)
                band = T._admit(T._prepare_admission(X.ZEEMAN_SZ_ALL,
                                                     param={"T": 0.5, "Nmat": nmat})).plan.rows
                self.assertEqual([r for r in band if r[1] == "file_read"], [])
                self.assertEqual(set(rows) - set(band), {row})

    def test_member_bytes_headers_only(self):
        """The admission sizes the row from the zip records and npy headers
        alone: no member is materialised (``np.load`` is never called), and
        the sum is the independent one of :func:`_member_bytes`."""
        path = _inflated(self, _GREENS["ZEEMAN_SZ_ALL"], 100)
        with mock.patch.object(np, "load", side_effect=AssertionError("np.load")) as load, \
                _PayloadReads() as spy:
            got = sc_io._so_bond_green_member_bytes(path)
        self.assertEqual(load.call_count, 0)
        self.assertEqual(spy.reads, [])
        self.assertEqual(got, _member_bytes(path))
        self.assertGreater(got, 1_000_000)
        d = _scratch(self)
        self.assertEqual(sc_io._so_bond_green_member_bytes(os.path.join(d, "absent.npz")), 0)
        npy = os.path.join(d, "green.npy")
        np.save(npy, np.zeros(3))
        self.assertEqual(sc_io._so_bond_green_member_bytes(npy), 0)

    def test_large_provenance(self):
        """A file whose provenance table is about 1 MB: its read stays
        within the sized row (m3: not within the fixed 256 KiB), and a cap
        between the two estimates refuses the run at admission, before the
        file is read."""
        path = _inflated(self, _GREENS["ZEEMAN_SZ_ALL"], 100)
        p = _memory_consumer(path, 32)
        _warm(p)
        seen = {}
        self.assertEqual(_memory_gate(p, seen=seen), [])
        self.assertIn(_file_read_row(path), seen["adm"].plan.rows)
        seen = {}
        with mock.patch.object(so_bond, "static_so_bond_resident_rows",
                               side_effect=_fixed_only) as brk:
            failures = _memory_gate(p, seen=seen)
        self.assertEqual(brk.call_count, 1)
        self.assertIn(("resident", "file_read", _FILE_READ_FIXED), seen["adm"].plan.rows)
        self.assertTrue(_has(failures, "green_static traced <= rows"), failures)
        # the cap: the sized estimate is refused, the fixed one would pass
        sized = T._admit(p).plan.peak_host
        extra = 1.25 * _FILE_READ_FACTOR * _member_bytes(path)
        cap_gb = (sized - extra / 2) / 2.0 ** 30
        small = _memory_consumer(_GREENS["ZEEMAN_SZ_ALL"], 32, {"bond_memory_cap_gb": cap_gb})
        self.assertLess(T._admit(small).plan.peak_host, sized - extra / 2)
        capped = _memory_consumer(path, 32, {"bond_memory_cap_gb": cap_gb})
        real = sc_io._so_borrowed_green
        with mock.patch.object(sc_bond, "_so_borrowed_green", side_effect=real) as read, \
                mock.patch.object(np, "load", side_effect=AssertionError("np.load")) as load:
            with self.assertRaises(ValueError) as cm:
                T._admit(capped)
        self.assertEqual((read.call_count, load.call_count), (0, 0))
        self.assertIn("but the cap is", str(cm.exception))
        self.assertIn("file_read", str(cm.exception))

    def test_phase_rows(self):
        for label, p, layout in self._rows():
            with self.subTest(row=label):
                seen = {}
                self.assertEqual(_memory_gate(p, seen=seen), [])
                self.assertEqual(seen["adm"].layout, layout)
                row = _file_read_row(p.cfg.green)
                self.assertIn(row, seen["adm"].plan.rows)
                peak, total = seen["values"]["green_static"]
                band = total - row[2]
                _LOG.info("%s: green_static traced %d, rows %d (%.3f), band rows %d (%.3f)",
                          label, peak, total, peak / total, band, peak / band)
                if label.endswith("Nmat 8"):
                    # the file_read row is needed here: the band rows alone are below the peak
                    self.assertGreater(peak, band)

    def test_m1_no_carrier_build_no_file_read(self):
        for label, p, _ in self._rows():
            with self.subTest(row=label):
                seen = {}
                real = so_bond.static_so_bond_rows
                # transient is the band-build/digest temporary the file read replaces
                drop = ("carrier_build", "transient")
                with mock.patch.object(so_bond, "static_so_bond_rows",
                                       side_effect=_dropping(real, drop)) as brk, \
                        mock.patch.object(so_bond, "static_so_bond_resident_rows",
                                          side_effect=_no_resident) as brk2:
                    failures = _memory_gate(p, seen=seen)
                self.assertEqual((brk.call_count, brk2.call_count), (1, 1))
                names = {r[1] for r in seen["adm"].plan.rows}
                self.assertNotIn("carrier_build", names)
                self.assertNotIn("transient", names)
                self.assertNotIn("file_read", names)
                self.assertTrue(_has(failures, "green_static traced <= rows"), failures)

    def test_m2_no_file_read(self):
        [(label, p, _)] = self._rows(("packed projected Nmat 8",))
        seen = {}
        with mock.patch.object(so_bond, "static_so_bond_resident_rows",
                               side_effect=_no_resident) as brk:
            failures = _memory_gate(p, seen=seen)
        self.assertEqual(brk.call_count, 1)
        self.assertNotIn("file_read", {r[1] for r in seen["adm"].plan.rows})
        # with the collector held off (_gc_held) the read's cycle garbage stays
        # alive into g2_static, so both phases fail, deterministically
        self.assertTrue(_has(failures, "green_static traced <= rows"), failures)
        self.assertTrue(_has(failures, "g2_static traced <= rows"), failures)


# ------------------------------------------------------------ public entry
_LADDER_BORROWED = (
    "[eliashberg] bond_channels=true in spin-orbital mode is a STATIC RPA-ladder bond "
    "dressing on the Green function read from bond_green = '%s' (its approximation level, "
    "e.g. FLEX-dressed, is that of the file). The bond chi-bar is built here from that Green "
    "function at RPA-ladder level; it is NOT the bond susceptibility of the run that wrote "
    "the file, so absolute lambda values are not comparable with FLEX/dynamic references.")
_FLEX_REFUSAL = (
    _W + ": [eliashberg] chi0q_mode = \"flex\" with bond_channels = true is not available: "
    "the spin-orbital bond path builds its bond bubble from a Green function (the band's, or "
    "the file named by bond_green) and reads no FLEX susceptibility. Remove chi0q_mode (or "
    "leave it unset); to use a FLEX Green function, set bond_green to the green.npz of a "
    "spin-orbital FLEX run with [mode.param] longitudinal_bond_channels = true.")
_FLEX_KEYS_UNUSED = (
    "[eliashberg] %s: not used for bond_channels=true in spin-orbital mode (the spin-orbital "
    "bond path reads no FLEX susceptibility).")
_ACCEPTED_INFO = ("[eliashberg] bond_green: '%s' (written by %s) passed the archive checks; "
                  "S_z leak %s; projected: %s.")


class TestBorrowedEntry(unittest.TestCase):
    """The public entry on a borrowed Green function: the run texts and the
    eigenvalue.dat lines."""

    def test_texts(self):
        self.assertEqual(sc_config._SO_BOND_FLEX_REFUSAL, _FLEX_REFUSAL)
        self.assertEqual(sc_config._SO_BOND_FLEX_KEYS_UNUSED, _FLEX_KEYS_UNUSED)
        self.assertEqual(sc_static._SO_BOND_RPA_LADDER_BORROWED, _LADDER_BORROWED)
        self.assertEqual(sc_bond._SO_BOND_GREEN_ACCEPTED_INFO, _ACCEPTED_INFO)

    def test_entry_outputs(self):
        path = _copy(self, "ZEEMAN_SZ_ALL")
        r = F.run_case(X.ZEEMAN_SZ_ALL, eliashberg={"Bond_Green": path}, fresh=True)
        keys = T._dat_keys(r.eigenvalue_dat)
        extra = ("bond_green", "bond_coeff_tail_applied", "bond_green_leak",
                 "bond_green_scf_converged")
        for k in T._BOND_DAT_KEYS + extra:
            self.assertIn(k, keys)
        self.assertEqual(keys["bond_green"], path)
        self.assertEqual(keys["bond_coeff_tail_applied"], "0.0")
        self.assertEqual(keys["bond_green_leak"], "0.000e+00")
        self.assertEqual(keys["bond_green_scf_converged"], "true")
        self.assertEqual(keys["bond_coeff_tail"], "1.0")
        self.assertEqual(keys["approximation"],
                         "static RPA-ladder bond dressing on a borrowed Green function")
        # the gap archive records the file's Green function and its state, not
        # the run's own (the band's mu differs: the state check would not see it)
        gap, green = r.gap_npz, _members(path)
        for k in ("fp_green", "state_beta", "state_nmat", "state_mu", "state_density",
                  "state_calc_mu", "state_target", "state_mu_tol"):
            self.assertEqual(np.asarray(gap[k]).tolist(), np.asarray(green[k]).tolist(), k)
        self.assertNotIn("fp_green_source", gap)
        band = F.run_case(X.ZEEMAN_SZ_ALL).gap_npz
        self.assertNotEqual(float(band["state_mu"]), float(green["state_mu"]))
        warns = [x.getMessage() for x in r.records if x.levelno == logging.WARNING]
        infos = [x.getMessage() for x in r.records]
        self.assertEqual(warns.count(_LADDER_BORROWED % path), 1)
        self.assertEqual(warns.count(sc_static._SO_BOND_RPA_LADDER), 0)
        accepted = _ACCEPTED_INFO % (path, _LABEL, "0.000e+00", "no")
        self.assertEqual(infos.count(accepted), 1)
        self.assertLess(infos.index(accepted), infos.index(_LADDER_BORROWED % path))

    def test_non_conserving_entry(self):
        path = _copy(self, "SOC_ALL_TYPES")
        r = F.run_case(F.SOC_ALL_TYPES, eliashberg={"Bond_Green": path}, fresh=True)
        keys = T._dat_keys(r.eigenvalue_dat)
        self.assertNotIn("bond_green_leak", keys)
        self.assertEqual(keys["bond_green"], path)
        self.assertEqual(keys["bond_green_scf_converged"], "true")
        infos = [x.getMessage() for x in r.records]
        self.assertEqual(infos.count(_ACCEPTED_INFO % (
            path, _LABEL, "not checked (the input does not conserve S_z)", "no")), 1)

    def test_internal_entry_unchanged(self):
        r = F.run_case(X.ZEEMAN_SZ_ALL, fresh=True)
        keys = T._dat_keys(r.eigenvalue_dat)
        for k in ("bond_green", "bond_green_leak", "bond_coeff_tail_applied",
                  "bond_green_scf_converged"):
            self.assertNotIn(k, keys)
        self.assertEqual(keys["approximation"], "static RPA-ladder bond dressing")
        warns = [x.getMessage() for x in r.records if x.levelno == logging.WARNING]
        self.assertEqual(warns.count(sc_static._SO_BOND_RPA_LADDER), 1)

    def test_refused_entry_logs_no_ladder(self):
        path = _copy(self, "ZEEMAN_SZ_ALL")
        seen = []

        class H(logging.Handler):
            def emit(self, record):
                seen.append(record)
        root = logging.getLogger()
        h = H(logging.INFO)
        level = root.level
        root.addHandler(h)
        root.setLevel(logging.INFO)
        try:
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                F.run_case(X.ZEEMAN_SZ_ALL, eliashberg={"Bond_Green": path},
                           param={"Nmat": 16}, fresh=True)
        finally:
            root.removeHandler(h)
            root.setLevel(level)
        self.assertIn("state differs", str(cm.exception))
        msgs = [r.getMessage() for r in seen]
        # the handler received the run's records, so the absences below hold
        self.assertIn("Calculating Green's function G(k, iwn) (spin-orbital bond)...", msgs)
        self.assertNotIn(_LADDER_BORROWED % path, msgs)
        self.assertNotIn(sc_static._SO_BOND_RPA_LADDER, msgs)


if __name__ == "__main__":
    unittest.main()
