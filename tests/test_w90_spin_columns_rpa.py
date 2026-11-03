#!/usr/bin/env python3
"""RPA / FLEX consuming the spin-resolved (9-column) Transfer format, issue
#323 phase 3: the acceptance cases of the request on the RPA side.

The one-body term is held as nd x nd in spin-block order and its spin
structure is judged from the matrix (spin-free / spin-diag / spinful), so a
spin-independent 9-column file takes the 7-column path, a Zeeman-like one
is spin-diag and conserves the flavour, and an SOC one is spinful and
promotes the auto scheme to general. FLEX imposes its spin-orbital
conditions on the spinful route; hwave_sc still refuses the format.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

from hwave.qlmsio import wan90
from hwave.solver import hamiltonian_k, spin_orbital
from tests.test_w90_spin_columns_uhfk import (
    _ORB, _closed, _spin_independent, _write_coulomb_intra, _write_geom,
    _write_plain_transfer, _write_spin_transfer)

_PARAM = {"T": 0.5, "Ncond": 2, "CellShape": [4, 1, 1], "Nmat": 16}

# the interaction set: CoulombIntra (always reduced-exact) plus an on-site
# inter-orbital Hund, a "conditional" type whose reduced treatment is exact
# only while the flavour is conserved -- the auto scheme consults the
# one-body structure only in the presence of such a type
_FILES = {"Geometry": "geom.dat", "Transfer": "transfer.dat",
          "CoulombIntra": "coulombintra.dat", "Hund": "hund.dat"}


def _write_hund(d, j):
    with open(os.path.join(d, "hund.dat"), "w") as f:
        f.write("Hund\n2\n1\n 1\n")
        f.write("0 0 0 1 2 {:.6f} 0.0\n0 0 0 2 1 {:.6f} 0.0\n".format(j, j))


def _ham(d, *, extern=False):
    import hwave.qlmsio.read_input_k as read_input_k
    files = dict(_FILES, path_to_input=d)
    if extern:
        files["Extern"] = "extern.dat"
    read_io = read_input_k.QLMSkInput({"path_to_input": d,
                                       "interaction": files})
    return read_io.get_param("ham"), read_io.get_param("green")


def _rpa(d, *, scheme="auto", subshape=(1, 1, 1), enable_spin_orbital=False,
         extern=False, coeff_extern=0.0):
    import hwave.solver.rpa as rpa_module
    ham, green = _ham(d, extern=extern)
    par = dict(_PARAM, SubShape=list(subshape))
    if extern:
        par["coeff_extern"] = coeff_extern
    info_mode = {"mode": "RPA", "param": par, "calc_scheme": scheme,
                 "calc_type": "ring", "enable_spin_orbital": enable_spin_orbital}
    return rpa_module.RPA(ham, {"print_level": 0, "print_step": 100},
                          info_mode), green


def _run_chi0q(d, *, scheme, enable_spin_orbital=False):
    """One hwave.qlms.run; returns the chi0q archive array."""
    import hwave.qlms as qlms
    out = os.path.join(d, "output")
    inp = {"log": {"print_level": 0},
           "mode": {"mode": "RPA", "param": dict(_PARAM, SubShape=[1, 1, 1]),
                    "calc_scheme": scheme, "calc_type": "ring",
                    "enable_spin_orbital": enable_spin_orbital},
           "file": {"input": {"path_to_input": d,
                              "interaction": dict(_FILES, path_to_input=d)},
                    "output": {"path_to_output": out, "chi0q": "chi0q.npz"}}}
    qlms.run(input_dict=inp)
    with np.load(os.path.join(out, "chi0q.npz")) as z:
        return np.array(z["chi0q"])


class _Fixture(unittest.TestCase):
    def setUp(self):
        self.dirs = []

    def tearDown(self):
        for d in self.dirs:
            shutil.rmtree(d, ignore_errors=True)

    def _dir(self, entries=None, *, plain=None, geom_norb=2, interleaved=None,
             extern_rows=None):
        d = tempfile.mkdtemp(prefix="w90spin_rpa_")
        self.dirs.append(d)
        _write_geom(d, geom_norb)
        _write_coulomb_intra(d, 2, 1.0)
        _write_hund(d, 0.3)
        if extern_rows is not None:
            with open(os.path.join(d, "extern.dat"), "w") as f:
                f.write("Extern\n2\n1\n 1\n")
                for row in extern_rows:
                    f.write(row + "\n")
        if plain is not None:
            _write_plain_transfer(d, 2, plain)
        elif interleaved is not None:
            with open(os.path.join(d, "transfer.dat"), "w") as f:
                f.write("Transfer SOI format\n{}\n1\n 1\n".format(geom_norb))
                for (r, a, b, s, t), v in _closed(interleaved).items():
                    ia = spin_orbital.interleaved_index(a - 1, s) + 1
                    ib = spin_orbital.interleaved_index(b - 1, t) + 1
                    f.write("{} {} {} {} {} {:.12f} {:.12f}\n".format(
                        *r, ia, ib, v.real, v.imag))
        else:
            _write_spin_transfer(d, 2, entries)
        return d


def _soc_entries(orb=_ORB):
    e = _spin_independent(orb)
    e[((0, 0, 0), 1, 2, 0, 1)] = 0.15 + 0.1j       # on-site SOC
    e[((1, 0, 0), 2, 1, 1, 0)] = 0.08               # Rashba-like
    return e


# orbital-diagonal transfer: the auto scheme's "flavour" is orbital AND
# spin, so the inter-orbital hopping of _ORB promotes every scheme to
# general whatever the spin structure (as it does for a 7-column file);
# the conservation cases below need a diagonal transfer
_DIAG = {k: v for k, v in _ORB.items() if k[1] == k[2]}


class TestSpinStructureFromTheMatrix(_Fixture):
    def _plain_vs_spin(self, orb, expect_scheme):
        d7 = self._dir(plain=orb)
        d9 = self._dir(_spin_independent(orb))
        s7, g7 = _rpa(d7)
        s9, g9 = _rpa(d9)
        self.assertTrue(s9.ham_info.spin_resolved_transfer)
        self.assertFalse(s7.ham_info.spinful_one_body)
        np.testing.assert_allclose(
            s9.ham_info.ham_trans_q,
            hamiltonian_k.embed_spin(s7.ham_info.ham_trans_q), atol=1e-12)
        s9._calc_epsilon_k(g9)
        self.assertEqual(s9.spin_mode, "spin-free")
        self.assertEqual(s9.preview_scheme(), s7.preview_scheme())
        self.assertEqual(s9.preview_scheme()[0], expect_scheme)
        np.testing.assert_allclose(_run_chi0q(d9, scheme="auto"),
                                   _run_chi0q(d7, scheme="auto"), atol=1e-10)

    def test_case1_spin_free_diagonal_matches_the_plain_file(self):
        self._plain_vs_spin(_DIAG, "reduced")

    def test_case1_spin_free_interorbital_matches_the_plain_file(self):
        # the inter-orbital hopping promotes both files alike
        self._plain_vs_spin(_ORB, "general")

    def test_case2_zeeman_is_spin_diag_and_conserves(self):
        e = _spin_independent(_DIAG)
        e[((0, 0, 0), 1, 1, 0, 0)] += -0.4
        e[((0, 0, 0), 1, 1, 1, 1)] += 0.4
        s, g = _rpa(self._dir(e))
        s._calc_epsilon_k(g)
        self.assertEqual(s.spin_mode, "spin-diag")
        self.assertEqual(s.preview_scheme(),
                         ("reduced", "auto:exact:diagonal_transfer"))

    def test_case3_soc_is_spinful_and_promotes_to_general(self):
        s, g = _rpa(self._dir(_soc_entries(_DIAG)))
        s._calc_epsilon_k(g)
        self.assertEqual(s.spin_mode, "spinful")
        self.assertEqual(s.preview_scheme(), ("general", "auto:mixed:transfer"))

    def test_case4_soc_matches_spin_orbital_mode(self):
        d9 = self._dir(_soc_entries())
        dso = self._dir(interleaved=_soc_entries(), geom_norb=4)
        s9, _ = _rpa(d9, scheme="general")
        sso, _ = _rpa(dso, scheme="general", enable_spin_orbital=True)
        np.testing.assert_allclose(s9.ham_info.ham_trans_q,
                                   sso.ham_info.ham_trans_q, atol=1e-12)
        np.testing.assert_allclose(
            _run_chi0q(d9, scheme="general"),
            _run_chi0q(dso, scheme="general", enable_spin_orbital=True),
            atol=1e-10)

    def test_case5_fold_keeps_the_down_block(self):
        e = _spin_independent(_ORB)
        e[((1, 0, 0), 1, 1, 1, 1)] = -0.5
        s, _ = _rpa(self._dir(e), subshape=(2, 1, 1))
        norb = s.ham_info.norb                     # 4 after the fold
        h = s.ham_info.ham_trans_q
        self.assertEqual(h.shape, (2, 2 * norb, 2 * norb))
        self.assertGreater(np.abs(h[:, :norb, :norb] - h[:, norb:, norb:]).max(), 0.1)
        np.testing.assert_allclose(h[:, :norb, norb:], 0.0, atol=1e-12)


class TestOneBodyRoutes(_Fixture):
    def test_extern_adds_sigma_z_to_the_spin_resolved_band(self):
        # the field through Extern (physical indices, h sigma^z) equals the
        # same field written into the spin-resolved rows
        h = 0.3
        d_ext = self._dir(_spin_independent(_DIAG),
                          extern_rows=["0 0 0 1 1 {:.6f} 0.0".format(h)])
        e = _spin_independent(_DIAG)
        e[((0, 0, 0), 1, 1, 0, 0)] += h
        e[((0, 0, 0), 1, 1, 1, 1)] -= h
        d_in = self._dir(e)
        s_ext, g_ext = _rpa(d_ext, extern=True, coeff_extern=1.0)
        s_in, g_in = _rpa(d_in)
        s_ext._calc_epsilon_k(g_ext)
        s_in._calc_epsilon_k(g_in)
        self.assertEqual(s_ext.spin_mode, "spin-diag")
        self.assertEqual(s_in.spin_mode, "spin-diag")
        np.testing.assert_allclose(np.sort(np.asarray(s_ext.H0_eigenvalue).ravel()),
                                   np.sort(np.asarray(s_in.H0_eigenvalue).ravel()),
                                   atol=1e-12)

    def test_trans_mod_seed_matches_the_file_band(self):
        # a trans_mod in the (nvol, nd, nd) spin-block form the run itself
        # would write reproduces the band of the spin-resolved file
        s, g = _rpa(self._dir(_soc_entries(_DIAG)), scheme="general")
        s._calc_epsilon_k(g)
        ref = np.sort(np.asarray(s.H0_eigenvalue).ravel())
        s2, g2 = _rpa(self._dir(_soc_entries(_DIAG)), scheme="general")
        g2 = dict(g2)
        g2["trans_mod"] = np.array(s2.ham_info.ham_trans_q)
        s2._calc_epsilon_k(g2)
        self.assertEqual(s2.spin_mode, "spinful")
        np.testing.assert_allclose(np.sort(np.asarray(s2.H0_eigenvalue).ravel()),
                                   ref, atol=1e-12)


class TestProvenanceMarker(_Fixture):
    """A chi0q computed from a spin-resolved Transfer is marked
    spin_resolved = True; a spin-resolved run requires the marker and every
    other run refuses it (the shapes coincide with a spin-free model of
    twice the orbitals)."""

    def _four_orbital_plain_dir(self):
        # 4-orbital spin-free model: its spin-free chi0q has the shape of
        # the spinful chi0q of a 2-orbital spin-resolved model
        d = tempfile.mkdtemp(prefix="w90spin_rpa_4orb_")
        self.dirs.append(d)
        _write_geom(d, 4)
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer\n4\n1\n 1\n")
            for a in range(1, 5):
                f.write("1 0 0 {0} {0} -1.0 0.0\n-1 0 0 {0} {0} -1.0 0.0\n".format(a))
        _write_coulomb_intra(d, 4, 1.0)
        _write_hund(d, 0.3)
        return d

    def _chi0q_file(self, d, **kw):
        _run_chi0q(d, **kw)
        return os.path.join(d, "output"), "chi0q.npz"

    def test_marker_written_and_round_trips(self):
        d = self._dir(_soc_entries(_DIAG))
        out, name = self._chi0q_file(d, scheme="general")
        with np.load(os.path.join(out, name)) as z:
            self.assertIn("spin_resolved", z.files)
            self.assertTrue(bool(z["spin_resolved"]))
            # not a spin-orbital archive
            self.assertFalse("spin_orbital" in z.files and bool(z["spin_orbital"]))
        s, _ = _rpa(d, scheme="general")
        info = s.read_init({"path_to_input": out, "chi0q_init": name})
        self.assertIn("chi0q", info)

    def test_other_run_refuses_a_spin_resolved_chi0q(self):
        d = self._dir(_soc_entries(_DIAG))
        out, name = self._chi0q_file(d, scheme="general")
        s4, _ = _rpa(self._four_orbital_plain_dir(), scheme="general")
        with self.assertRaisesRegex(ValueError, "spin_resolved = True"):
            s4.read_init({"path_to_input": out, "chi0q_init": name})

    def test_spin_resolved_run_refuses_an_unmarked_chi0q(self):
        d4 = self._four_orbital_plain_dir()
        out, name = self._chi0q_file(d4, scheme="general")
        s, _ = _rpa(self._dir(_soc_entries(_DIAG)), scheme="general")
        with self.assertRaisesRegex(ValueError, "spin_resolved"):
            s.read_init({"path_to_input": out, "chi0q_init": name})


class TestFlexGates(_Fixture):
    def _flex(self, d, scheme, extra=None):
        import hwave.solver.flex as flex_module
        ham, _ = _ham(d)
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=1, Mix=1.0)
        par.update(extra or {})
        info_mode = {"mode": "FLEX", "param": par, "calc_scheme": scheme,
                     "calc_type": "ring", "enable_spin_orbital": False}
        return flex_module.FLEX(ham, {}, info_mode)

    def test_reduced_refuses_spinful_one_body(self):
        with self.assertRaisesRegex(ValueError, "spin-resolved"):
            self._flex(self._dir(_soc_entries()), "reduced")

    def test_general_requires_the_spin_orbital_conditions(self):
        with self.assertRaisesRegex(ValueError, "flex_hartree_fock"):
            self._flex(self._dir(_soc_entries()), "general")


class TestExportRoundTrip(_Fixture):
    def test_export_writes_the_spin_resolved_format_back(self):
        d = self._dir(_soc_entries())
        s, _ = _rpa(d)
        path = os.path.join(d, "exported.dat")
        s.ham_info._export_interaction("Transfer", path)
        self.assertEqual(wan90.read_w90_header_marker(path), ("w90spin", 1))
        back = wan90.read_w90(path, allow_spin_columns=True)
        table = {k: v for k, v in s.ham_info.param_ham["Transfer"].items()
                 if abs(v) > 1.0e-12}
        self.assertEqual(set(back), set(table))
        for key, v in table.items():
            self.assertAlmostEqual(back[key], v, places=10)


if __name__ == "__main__":
    unittest.main()
