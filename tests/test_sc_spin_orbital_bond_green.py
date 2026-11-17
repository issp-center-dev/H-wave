"""hwave_sc in spin-orbital mode with bond-resolved channels and a Green
function read from ``[eliashberg] bond_green`` (issue #270): the
configuration, the Green function built from the green.npz of a
spin-orbital bond FLEX run, the refusals of that file and the tail policy.

The producer files are written once (:func:`setUpModule`, converged
spin-orbital bond FLEX runs at the consumer's state point) and removed in
:func:`tearDownModule`. The rows vary the consumer (the prepared run) or a
copy of the file, and call ``sc._so_bond_green`` directly (no solve); the
on-site row, whose run record has no band, calls ``sc._so_borrowed_green``.
"""

import contextlib
import io
import logging
import os
import shutil
import tempfile
import unittest
from unittest import mock

import numpy as np
from requests.structures import CaseInsensitiveDict

import hwave.sc as sc
from hwave.solver import fingerprint, so_bond, so_bond_archive
from tests import so_bond_archive_fixtures as AF
from tests import so_bond_entry_fixtures as F
from tests import so_bond_flex_fixtures as X
from tests import test_sc_spin_orbital_bond as T
from tests.so_bond_oracles import CELL
from tests.test_so_entry_rules import _PayloadReads

_W = "[mode] enable_spin_orbital = true (hwave_sc, issue #83)"
_LABEL = "hwave FLEX (spin-orbital, bond channels)"
_LOG = logging.getLogger("qlms.test_sc_spin_orbital_bond_green")
#: the producers of the module (fixture, sigma.npz written)
_PRODUCERS = {"ZEEMAN_SZ_ALL": X.ZEEMAN_SZ_ALL, "SOC_ALL_TYPES": F.SOC_ALL_TYPES}
_GREENS = {}
_ROOT = []
_DROP = object()


def setUpModule():
    """The green.npz (and sigma.npz) of a converged spin-orbital bond FLEX
    run of each producer, at T 0.5, Nmat 32 and the fixture's filling."""
    root = tempfile.mkdtemp(prefix="sc_so_bond_green_")
    _ROOT.append(root)
    with contextlib.redirect_stdout(io.StringIO()):
        for name, fx in _PRODUCERS.items():
            d = os.path.join(root, name)
            os.makedirs(d)
            _GREENS[name] = X.write_sc_green(d, fx, sigma=True)


def tearDownModule():
    for root in _ROOT:
        shutil.rmtree(root, True)
    _ROOT.clear()
    _GREENS.clear()


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
    return sc._so_bond_green(p.epsilon_k, p.so_run, norb_phys=p.norb_phys, filling=p.filling,
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
    loggers = [sc.logger, so_bond.logger]
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
        return sc._read_so_bond_config(CaseInsensitiveDict(eli))

    def test_path_read(self):
        cfg = self._parse({"bond_channels": True, "Bond_Green": "g.npz"})
        self.assertEqual(cfg.green, "g.npz")
        self.assertIsNone(self._parse({"bond_channels": True}).green)
        # the band-only constructor of seven fields stays valid
        self.assertIsNone(sc._SOBondConfig(True, "auto", None, None, 1.0, "default", True).green)

    def test_empty_string(self):
        with self.assertRaises(ValueError) as cm:
            self._parse({"bond_channels": True, "bond_green": ""})
        self.assertEqual(str(cm.exception),
                         "[eliashberg] bond_green must be a path to a Green-function npz "
                         "file; got an empty string. Omit the key to build the bare Green "
                         "function from the transfer Hamiltonian.")
        # the non-spin-orbital reader keeps the same text
        with self.assertRaises(ValueError) as cm2:
            sc._read_bond_config({"bond_channels": True, "bond_green": ""})
        self.assertEqual(str(cm2.exception), str(cm.exception))

    def test_refusal_order(self):
        on = {"bond_channels": True}
        with self.assertRaises(ValueError) as cm:
            self._parse(dict(on, bond_green="g.npz", bond_diagnostics=True))
        self.assertEqual(str(cm.exception), sc._SO_BOND_DIAGNOSTICS_REFUSAL)


# ------------------------------------------------------------ accepted files
class TestBorrowedGreen(unittest.TestCase):
    """The Green function of a static spin-orbital bond run read from the
    green.npz of the spin-orbital bond FLEX run at its state point: the
    carrier holds the file's G (no copy, no tail), the record holds the
    file's digest, mu and state."""

    def _run(self, name):
        path = _copy(self, name)
        got = []
        real = sc._so_borrowed_green

        def spy(*a, **kw):
            got.append(real(*a, **kw))
            return got[-1]
        with _records() as rec, mock.patch.object(sc, "_so_borrowed_green", side_effect=spy):
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
        self.assertEqual(_count(rec, sc._SO_BOND_TAIL_DEFAULT_INFO), 1)
        self.assertEqual(_count(rec, sc._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
                                "0.000e+00", "no"), 1)

    def test_clean_non_conserving(self):
        path, g, borrowed, rec = self._run("SOC_ALL_TYPES")
        self.assertIsNone(g.leak)
        self.assertIsNone(borrowed.leak)
        self.assertIs(g.conserving, False)
        self.assertEqual(g.fp_green, str(_members(path)["fp_green"]))
        self.assertIs(g.carrier.deflated_kw, borrowed.used)
        self.assertEqual(_warnings(rec), [])
        self.assertEqual(_count(rec, sc._SO_BOND_GREEN_ACCEPTED_INFO, path, _LABEL,
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
        plain = sc._so_writer_members("sc_so_bond_static", "test", p.so_run, fp_file, state)
        self.assertNotIn("fp_green_source", plain)
        members = sc._so_writer_members("sc_so_bond_static", "test", p.so_run, fp_used, state,
                                        fp_green_source=fp_file)
        self.assertEqual(str(members["fp_green_source"]), fp_file)
        self.assertEqual(str(members["fp_green"]), fp_used)
        self.assertEqual(set(members) - set(plain), {"fp_green_source"})
        for listed in (False, True):
            base = sc._so_bond_gap_members(listed)
            self.assertNotIn("fp_green_source", base)
            self.assertEqual(sc._so_bond_gap_members(listed, green_source=True),
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
        self.assertEqual(_count(rec, sc._SO_BOND_GREEN_UNCONVERGED, path), 0)
        return e

    def test_texts(self):
        self.assertEqual(sc.T_BG_MISSING, "[eliashberg] bond_green file not found: {path}")
        self.assertEqual(sc.T_BG_SHAPE,
                         "'{path}' is malformed or truncated: green shape differs from this "
                         "run ({stored} vs {expected}). Regenerate it with {label}.")
        self.assertEqual(sc.T_BG_NONFINITE,
                         "'{path}' is malformed or truncated: green has non-finite values at "
                         "frequency index {n}. Regenerate it with {label}.")
        self.assertEqual(sc.T_BG_MARKER,
                         "'{path}' is malformed or truncated: {member} is missing or older "
                         "than this build requires. Regenerate it with {label}.")
        self.assertEqual(sc.T_BG_KGRID,
                         "'{path}' is malformed or truncated: wavevector_index is not the k "
                         "grid of CellShape {shape}. Regenerate it with {label}.")
        self.assertEqual(sc._SO_BOND_GREEN_ONSITE_REMEDY,
                         "[eliashberg] bond_green needs the green.npz of a spin-orbital FLEX "
                         "run with [mode.param] longitudinal_bond_channels = true (the "
                         "bond-resolved run); re-run FLEX with it.")
        self.assertEqual(sc._SO_BOND_GREEN_IS_ARCHIVE,
                         _W + ": [eliashberg] bond_green = '{path}' is the bond archive of a "
                         "spin-orbital FLEX run, not a Green function; set bond_green to its "
                         "companion green file '{green}'.")
        self.assertEqual(sc._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN,
                         _W + ": [eliashberg] bond_green = '{path}' is the bond archive of a "
                         "spin-orbital FLEX run, not a Green function, and it records no "
                         "readable companion green file; set bond_green to the green.npz "
                         "written by that FLEX run ([file.output] green).")
        self.assertEqual(sc._SO_BOND_GREEN_NONSO_ARCHIVE,
                         _W + ": [eliashberg] bond_green = '{path}' is a bond archive, not a "
                         "Green function; set bond_green to the green.npz of a spin-orbital "
                         "FLEX run with [mode.param] longitudinal_bond_channels = true.")
        self.assertEqual(sc._SO_BOND_GREEN_UNCONVERGED,
                         "[eliashberg] bond_green '%s' records scf_converged = false: the FLEX "
                         "run did not converge; its Green function is used as stored.")
        self.assertEqual(sc._SO_BOND_GREEN_ACCEPTED_INFO,
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
        self._row(p, index, sc._SO_BOND_GREEN_IS_ARCHIVE.format(
            path=index, green=os.path.join(d, "sub", "green.npz")), cls=ValueError)
        self.assertIn("/sub/green.npz'.", str(self._refuse(p, index, ValueError)[0]))
        no_file = _edited(self, index, green_file=_DROP)
        empty = _edited(self, index, green_file=np.str_(""))
        for path in (no_file, empty):
            with self.subTest(path=os.path.basename(path)):
                self._row(p, path, sc._SO_BOND_GREEN_IS_ARCHIVE_NO_GREEN.format(path=path),
                          cls=ValueError)
        nonso = os.path.join(d, "bond_archive.npz")
        np.savez(nonso, bond_archive_schema=np.int64(3))
        self._row(p, nonso, sc._SO_BOND_GREEN_NONSO_ARCHIVE.format(path=nonso), cls=ValueError)

    def test_onsite_green(self):
        ctx = T._context()
        run = T._run(ctx, bond=sc._so_bond_run(ctx, CELL, None))
        path = T._onsite_green_archive(self, run)
        with _PayloadReads() as spy:
            with self.assertRaises(fingerprint.SOArchiveError) as cm:
                sc._so_borrowed_green(path, run, norb_phys=1, nmat=T._RUN_NMAT,
                                      shape=(2, 2, 1), conserving_input=True)
        e = cm.exception
        self.assertEqual(str(e), "{} was produced by test FLEX producer (so_onsite), but this "
                         "run needs a so_bond artifact. ".format(path)
                         + sc._SO_BOND_GREEN_ONSITE_REMEDY)
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
                e = self._row(p, path, sc.T_BG_MARKER.format(path=path, member=member,
                                                             label=_LABEL))
                self.assertIsNone(e.check)
        index = _members(src)["wavevector_index"]
        path = _edited(self, src, wavevector_index=index[::-1].copy())
        self._row(p, path, "'{}' is malformed or truncated: wavevector_index is not the k grid "
                  "of CellShape [4, 4, 1]. Regenerate it with {}.".format(path, _LABEL))

    def test_unconverged(self):
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        path = _edited(self, src, scf_converged=np.bool_(False))
        with _records() as rec:
            g = _borrow(p, path)
        self.assertIs(g.scf_converged, False)
        self.assertEqual(_warnings(rec), [sc._SO_BOND_GREEN_UNCONVERGED % path])
        # an edited green on top: the digest refuses it, no unconverged WARNING
        edited = _edited(self, _perturbed(self, src, 1e-6, restamp=False),
                         scf_converged=np.bool_(False))
        self._row(p, edited, self._mismatch(edited), read=True)

    def _mismatch(self, path):
        m = _members(path)
        return ("'{}' (written by {}): fp_green mismatch: the Green function read from it has "
                "digest {}..., the file records {}... (both at fingerprint_version 2); its "
                "green member is not the one the run that wrote it recorded. Regenerate it "
                "with {}.".format(path, _LABEL, str(_restamped(m["green"]))[:12],
                                  str(m["fp_green"])[:12], _LABEL))

    def test_payload(self):
        p = _base("ZEEMAN_SZ_ALL")
        src = _GREENS["ZEEMAN_SZ_ALL"]
        G = _members(src)["green"]
        short = G[:, :30].copy()
        path = _edited(self, src, green=short, fp_green=_restamped(short))
        self._row(p, path, sc.T_BG_SHAPE.format(path=path, stored=(1, 30, 16, 2, 2),
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
        e = self._row(p, path, self._mismatch(path), read=True)
        self.assertEqual(e.check, "green")


# ------------------------------------------------------------ tail policy
class TestTailPolicy(unittest.TestCase):
    """No tail is applied to a Green function read from a file; the
    default bond_coeff_tail is noted at INFO, a configured non-zero one is
    a WARNING."""

    def test_tail_policy(self):
        self.assertEqual(sc._SO_BOND_TAIL_NOT_APPLIED,
                         "[eliashberg] bond_coeff_tail = %r is not applied with bond_green = "
                         "'%s': the tail correction needs the analytic asymptote of the band "
                         "Green function, which a Green function read from a file does not "
                         "provide (bond_coeff_tail_applied = 0).")
        self.assertEqual(sc._SO_BOND_TAIL_DEFAULT_INFO,
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
                self.assertEqual(_count(rec, sc._SO_BOND_TAIL_DEFAULT_INFO), default_info)
                tail = [w for w in _warnings(rec) if "bond_coeff_tail" in w]
                if shown is None:
                    self.assertEqual(tail, [])
                else:
                    self.assertEqual(tail, [
                        "[eliashberg] bond_coeff_tail = {} is not applied with bond_green = "
                        "'{}': the tail correction needs the analytic asymptote of the band "
                        "Green function, which a Green function read from a file does not "
                        "provide (bond_coeff_tail_applied = 0).".format(shown, path)])


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
        self.assertEqual(sc._SO_BOND_FLEX_REFUSAL, _FLEX_REFUSAL)
        self.assertEqual(sc._SO_BOND_FLEX_KEYS_UNUSED, _FLEX_KEYS_UNUSED)
        self.assertEqual(sc._SO_BOND_RPA_LADDER_BORROWED, _LADDER_BORROWED)
        self.assertEqual(sc._SO_BOND_GREEN_ACCEPTED_INFO, _ACCEPTED_INFO)

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
        warns = [x.getMessage() for x in r.records if x.levelno == logging.WARNING]
        infos = [x.getMessage() for x in r.records]
        self.assertEqual(warns.count(_LADDER_BORROWED % path), 1)
        self.assertEqual(warns.count(sc._SO_BOND_RPA_LADDER), 0)
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
        self.assertEqual(warns.count(sc._SO_BOND_RPA_LADDER), 1)

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
        self.assertNotIn(_LADDER_BORROWED % path, msgs)
        self.assertNotIn(sc._SO_BOND_RPA_LADDER, msgs)


if __name__ == "__main__":
    unittest.main()
