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
                self.assertEqual(sorted(got), sorted(ref))
                for name in ref:
                    self.assertEqual(got[name], ref[name], name)

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


class TestUHFkInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_uhfk_config())

    def test_lower_case_geometry_alone(self):
        # the original report: only 'geometry' lower-cased
        cfg = _uhfk_config()
        ref = self._run(cfg, "canonical")
        inter = cfg["file"]["input"]["interaction"]
        inter["geometry"] = inter.pop("Geometry")
        self.assertEqual(self._run(cfg, "geometry"), ref)

    def test_duplicate_spellings_refused(self):
        for key in ("Geometry", "Transfer", "CoulombIntra", "PairLift"):
            with self.subTest(key=key):
                self._check_duplicate_refused(_uhfk_config(), key)


class TestUHFrInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_uhfr_config())

    def test_duplicate_spellings_refused(self):
        self._check_duplicate_refused(_uhfr_config(), "CoulombIntra")


class TestRPAInteractionKeyCase(_CaseMixin, unittest.TestCase):
    def test_spellings_give_identical_outputs(self):
        self._check_spellings(_rpa_config())

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
