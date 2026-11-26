#!/usr/bin/env python3
"""[file.input.interaction] keys are case-insensitive.

Through the public entry ``hwave.qlms.run``: a lower-case or mixed-case
spelling of every interaction key (``geometry``, ``coulombINTER``, ...) gives
byte-identical outputs to the canonical spelling, for UHFk, UHFr, RPA and
FLEX, and two keys that differ only in case are refused. A lower-case
``geometry`` used to crash UHFk in ``_check_cellsize``.
"""

import copy
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlms as qlms

_HERE = os.path.dirname(os.path.abspath(__file__))
UHFK_DIR = os.path.join(_HERE, "uhfk")
UHFR_DIR = os.path.join(_HERE, "uhfr")
RPA_DIR = os.path.join(_HERE, "rpa", "input")


def _lower(key):
    return key.lower()


def _mixed(key):
    # 'Geometry' -> 'gEOMETRY': differs from the canonical spelling in
    # every letter
    return key.swapcase()


def _respell(cfg, spell):
    cfg = copy.deepcopy(cfg)
    inter = cfg["file"]["input"]["interaction"]
    cfg["file"]["input"]["interaction"] = {spell(k): v for k, v in inter.items()}
    return cfg


def _uhfk_config():
    # several UHFk two-body types at once: the case directories share one
    # geometry and one transfer
    inter = {
        "path_to_input": UHFK_DIR,
        "Geometry": "CoulombIntra/geom.dat",
        "Transfer": "CoulombIntra/transfer.dat",
        "CoulombIntra": "CoulombIntra/coulombintra.dat",
        "CoulombInter": "CoulombInter/coulombinter.dat",
        "Ising": "Ising/ising.dat",
        "Exchange": "Exchange/exchange.dat",
        "PairHop": "PairHop/pairhop.dat",
        "PairLift": "PairLift/pairlift.dat",
    }
    return {
        "log": {"print_level": 0},
        "mode": {"mode": "UHFk",
                 "param": {"2Sz": 0, "Ncond": 16, "IterationMax": 20,
                           "EPS": 12, "Mix": 0.5, "RndSeed": 123456789,
                           "T": 0.0, "CellShape": [4, 4, 1],
                           "SubShape": [2, 2, 1]}},
        "file": {"input": {"path_to_input": "", "interaction": inter},
                 "output": {"energy": "energy.dat", "eigen": "eigen.dat",
                            "green": "green.dat"}},
    }


def _uhfr_config():
    d = os.path.join(UHFR_DIR, "coulombintra")
    return {
        "log": {"print_level": 0},
        "mode": {"mode": "UHFr",
                 "param": {"Nsite": 16, "Ncond": 16,
                           "IterationMax": 20, "EPS": 8, "Mix": 0.5,
                           "RndSeed": 123456789, "T": 0.0}},
        "file": {"input": {"path_to_input": d,
                           "OneBodyG": "greenone.def",
                           "Initial": "output_ref/green_init.dat",
                           "interaction": {"Trans": "trans.def",
                                           "CoulombIntra": "coulombintra.def"}},
                 "output": {"energy": "energy.dat", "eigen": "eigen.dat",
                            "green": "green.dat"}},
    }


def _rpa_config(mode="RPA"):
    inter = {
        "path_to_input": RPA_DIR,
        "Geometry": "geom.dat",
        "Transfer": "transfer.dat",
        "CoulombIntra": "coulombintra.dat",
        "CoulombInter": "coulombinter.dat",
        "Extern": "extern.dat",
        # off-site bond types from the UHFk cases (same one-orbital lattice)
        "Ising": os.path.join(UHFK_DIR, "Ising", "ising.dat"),
        "Exchange": os.path.join(UHFK_DIR, "Exchange", "exchange.dat"),
        "PairHop": os.path.join(UHFK_DIR, "PairHop", "pairhop.dat"),
        "PairLift": os.path.join(UHFK_DIR, "PairLift", "pairlift.dat"),
    }
    param = {"T": 0.5, "filling": 0.5, "CellShape": [4, 4, 1],
             "SubShape": [1, 1, 1], "Nmat": 16, "coeff_extern": 0.1}
    out = {"chiq": "chiq.npz", "chi0q": "chi0q.npz"}
    if mode == "FLEX":
        # FLEX general is spin-free (no Extern field) and refuses off-site
        # Exchange / PairHop
        del inter["Extern"], param["coeff_extern"]
        del inter["Exchange"], inter["PairHop"]
        param.update({"IterationMax": 3, "Mix": 0.5, "EPS": 10})
        out = {"chiq": "chiq.npz", "sigma": "sigma.npz"}
    return {
        "log": {"print_level": 0},
        "mode": {"mode": mode, "calc_scheme": "general", "param": param},
        "file": {"input": {"path_to_input": RPA_DIR, "interaction": inter},
                 "output": out},
    }


#: numpy 1.x runs of the SAME canonical input already differ in the last
#: bits from run to run (observed with numpy 1.26.4 on macOS arm64: up to
#: 8.9e-16 in an RPA chiq of magnitude 1.3, and up to 3.7e-13 after the
#: three FLEX iterations), so byte identity of the npz arrays is only
#: asserted on numpy 2, where repeated runs are bitwise reproducible; on
#: numpy 1.x the arrays must agree within 1e-11 times their largest
#: magnitude (about 30 times the largest observed repeat-run noise). Text
#: outputs are always compared byte for byte, and the reader tests below
#: check the solver input bit for bit on every numpy.
_NUMPY1_NOISE = 1e-11
_NUMPY2 = int(np.__version__.split(".")[0]) >= 2


def _read_outputs(out_dir):
    """Every output file: text files as bytes, npz files as their arrays'
    (dtype, shape, bytes) -- the zip container itself carries timestamps."""
    result = {}
    for name in sorted(os.listdir(out_dir)):
        path = os.path.join(out_dir, name)
        if name.endswith(".npz"):
            with np.load(path, allow_pickle=False) as z:
                result[name] = {k: (z[k].dtype.str, z[k].shape, z[k].tobytes())
                                for k in sorted(z.files)}
        else:
            with open(path, "rb") as f:
                result[name] = f.read()
    return result


class _CaseMixin:
    def setUp(self):
        self._tmp = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self._tmp, ignore_errors=True)

    def _run(self, cfg, name):
        cfg = copy.deepcopy(cfg)
        out = os.path.join(self._tmp, name)
        cfg["file"]["output"]["path_to_output"] = out
        qlms.run(input_dict=cfg)
        return _read_outputs(out)

    def _check_spellings(self, cfg):
        ref = self._run(cfg, "canonical")
        self.assertTrue(ref)
        for label, spell in (("lower", _lower), ("mixed", _mixed)):
            with self.subTest(spelling=label):
                got = self._run(_respell(cfg, spell), label)
                self._assert_same_outputs(got, ref)

    def _assert_same_outputs(self, got, ref):
        self.assertEqual(sorted(got), sorted(ref))
        for name in ref:
            if isinstance(ref[name], bytes) or _NUMPY2:
                self.assertEqual(got[name], ref[name], name)
                continue
            self.assertEqual(sorted(got[name]), sorted(ref[name]), name)
            for key, (dtype, shape, raw) in ref[name].items():
                gdtype, gshape, graw = got[name][key]
                self.assertEqual((gdtype, gshape), (dtype, shape), name)
                a = np.frombuffer(raw, dtype=dtype)
                b = np.frombuffer(graw, dtype=dtype)
                if a.dtype.kind in "fc":
                    scale = max(1.0, float(np.max(np.abs(a)))) if a.size else 1.0
                    np.testing.assert_allclose(b, a, rtol=0,
                                               atol=_NUMPY1_NOISE * scale,
                                               err_msg=name + ":" + key)
                else:
                    np.testing.assert_array_equal(b, a, err_msg=name + ":" + key)

    def _check_duplicate_refused(self, cfg, key):
        cfg = copy.deepcopy(cfg)
        inter = cfg["file"]["input"]["interaction"]
        inter[key.lower()] = inter[key]
        cfg["file"]["output"]["path_to_output"] = os.path.join(self._tmp, "dup")
        with self.assertRaisesRegex(
                ValueError, r"\[file.input.interaction\] {} is given more "
                r"than once with different spellings \('{}', '{}'\)".format(
                    key, key, key.lower())):
            qlms.run(input_dict=cfg)


class TestCanonicalizeInteractionKeys(unittest.TestCase):
    """The helper itself, over every key either reader accepts."""

    def test_every_key_every_spelling(self):
        from hwave.qlmsio import keys
        for key in keys.ALL_KEYS:
            for spelling in (key, key.lower(), key.upper(), key.swapcase()):
                with self.subTest(key=key, spelling=spelling):
                    value = object()
                    out = keys.canonicalize_interaction_keys(
                        {"Unknown_Key": 1, spelling: value})
                    # order kept, value passed through, unknown key untouched
                    self.assertEqual(list(out), ["Unknown_Key", key])
                    self.assertIs(out[key], value)

    def test_every_key_duplicate_refused(self):
        from hwave.qlmsio import keys
        for key in keys.ALL_KEYS:
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, "given more than once"):
                    keys.canonicalize_interaction_keys(
                        {key: "a", key.swapcase(): "a"})

    def test_canonical_input_unchanged(self):
        from hwave.qlmsio import keys
        table = {k: k + ".dat" for k in keys.KSPACE_KEYS}
        self.assertEqual(list(keys.canonicalize_interaction_keys(table).items()),
                         list(table.items()))


class TestCallerInputNotRewritten(_CaseMixin, unittest.TestCase):
    """The public entries do not rewrite the caller's interaction table."""

    def _check(self, entry, cfg):
        inter = cfg["file"]["input"]["interaction"]
        snapshot = list(inter.items())
        try:
            entry(cfg)
        except Exception:
            pass    # only the caller's mapping matters here
        self.assertIs(cfg["file"]["input"]["interaction"], inter)
        self.assertEqual(list(inter.items()), snapshot)

    def test_qlms_run(self):
        cfg = _respell(_uhfk_config(), _mixed)
        cfg["file"]["output"]["path_to_output"] = os.path.join(self._tmp, "o")
        self._check(lambda c: qlms.run(input_dict=c), cfg)

    def test_calc_eliashberg(self):
        from hwave.sc import calc_eliashberg
        cfg = _respell(_rpa_config(), _mixed)
        cfg["file"]["output"]["path_to_output"] = os.path.join(self._tmp, "o")
        self._check(calc_eliashberg, cfg)


class TestUHFkInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_uhfk_config())

    def test_lower_case_geometry_alone(self):
        # the original report: only 'geometry' lower-cased
        cfg = _uhfk_config()
        ref = self._run(cfg, "canonical")
        inter = cfg["file"]["input"]["interaction"]
        inter["geometry"] = inter.pop("Geometry")
        self._assert_same_outputs(self._run(cfg, "geometry"), ref)

    def test_reader_tables_identical(self):
        # the solvers receive the same tables, bit for bit and in the same
        # order, whatever the spelling
        from hwave.qlmsio.read_input_k import QLMSkInput

        def tables(cfg):
            ham = QLMSkInput(cfg["file"]["input"]).get_param("ham")
            return [(k, repr(v)) for k, v in ham.items()]

        ref = tables(_uhfk_config())
        self.assertEqual([k for k, _ in ref],
                         [k for k in _uhfk_config()["file"]["input"]["interaction"]
                          if k != "path_to_input"])
        for spell in (_lower, _mixed):
            self.assertEqual(tables(_respell(_uhfk_config(), spell)), ref)

    def test_duplicate_spellings_refused(self):
        for key in ("Geometry", "Transfer", "CoulombIntra", "PairLift"):
            with self.subTest(key=key):
                self._check_duplicate_refused(_uhfk_config(), key)


class TestUHFrInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_uhfr_config())

    def test_reader_tables_identical(self):
        # the real-space reader, fed the file list qlms.run builds, gets the
        # same tables bit for bit (it reads every type, present or not)
        from unittest import mock
        import hwave.qlmsio.read_input as read_input
        seen = []
        orig = read_input.QLMSInput.__init__

        def spy(self, *a, **k):
            orig(self, *a, **k)
            seen.append([(t, repr(v)) for t, v in self.ham_param.items()])

        for cfg in (_uhfr_config(), _respell(_uhfr_config(), _lower),
                    _respell(_uhfr_config(), _mixed)):
            cfg["file"]["output"]["path_to_output"] = os.path.join(self._tmp, "o")
            with mock.patch.object(read_input.QLMSInput, "__init__", spy):
                qlms.run(input_dict=cfg)
        self.assertEqual(len(seen), 3)
        self.assertEqual(seen[1], seen[0])
        self.assertEqual(seen[2], seen[0])

    def test_duplicate_spellings_refused(self):
        self._check_duplicate_refused(_uhfr_config(), "CoulombIntra")


class TestRPAInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_rpa_config())

    def test_reader_tables_identical(self):
        from hwave.qlmsio.read_input_k import QLMSkInput

        def tables(cfg):
            ham = QLMSkInput(cfg["file"]["input"], solver_type="RPA").get_param("ham")
            return [(k, repr(v)) for k, v in ham.items()]

        ref = tables(_rpa_config())
        for spell in (_lower, _mixed):
            self.assertEqual(tables(_respell(_rpa_config(), spell)), ref)

    def test_duplicate_spellings_refused(self):
        self._check_duplicate_refused(_rpa_config(), "Extern")


class TestFLEXInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_rpa_config("FLEX"))

    def test_duplicate_spellings_refused(self):
        self._check_duplicate_refused(_rpa_config("FLEX"), "CoulombInter")


class TestOtherEntriesInteractionKeyCase(_CaseMixin, unittest.TestCase):
    """hwave_dos and hwave_sc read the same table."""

    @unittest.skipUnless(
        __import__("importlib").util.find_spec("libtetrabz") is not None,
        "libtetrabz not installed")
    def test_dos_lower_case_keys(self):
        from hwave.dos import calc_dos
        cfg = _uhfk_config()
        out = os.path.join(self._tmp, "uhfk")
        cfg["file"]["output"]["path_to_output"] = out
        qlms.run(input_dict=copy.deepcopy(cfg))
        ref = calc_dos(copy.deepcopy(cfg), ene_num=21)
        got = calc_dos(_respell(cfg, _lower), ene_num=21)
        self.assertEqual(got.ene.tobytes(), ref.ene.tobytes())
        self.assertEqual(got.dos.tobytes(), ref.dos.tobytes())

    def test_sc_spellings_give_identical_outputs(self):
        from hwave.sc import calc_eliashberg
        from tests.test_sc_trans_mod import _input, _write_model
        d = os.path.join(self._tmp, "model")
        os.makedirs(d)
        _write_model(d, 1)

        def run(name, spell):
            out = os.path.join(self._tmp, name)
            cfg = _input(d, out)
            if spell is not None:
                cfg = _respell(cfg, spell)
            calc_eliashberg(cfg)
            return _read_outputs(out)

        ref = run("canonical", None)
        self.assertTrue(ref)
        for label, spell in (("lower", _lower), ("mixed", _mixed)):
            with self.subTest(spelling=label):
                self._assert_same_outputs(run(label, spell), ref)

    def test_sc_duplicate_spellings_refused(self):
        from hwave.sc import calc_eliashberg
        cfg = _rpa_config()
        inter = cfg["file"]["input"]["interaction"]
        inter["GEOMETRY"] = inter["Geometry"]
        with self.assertRaisesRegex(
                ValueError, r"\[file.input.interaction\] Geometry is given "
                r"more than once with different spellings "
                r"\('Geometry', 'GEOMETRY'\)"):
            calc_eliashberg(cfg)


if __name__ == "__main__":
    unittest.main()
