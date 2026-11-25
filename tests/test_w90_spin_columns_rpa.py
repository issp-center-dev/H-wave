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


class TestLineTwoNorbMustMatchGeometry(_Fixture):
    """Line 2 of a spin-resolved Transfer is the physical orbital count the
    reader folds with (s*N + a); RPA / FLEX decode the keys with the
    Geometry norb, so a different value is refused."""

    def _dir_header(self, header_norb, entries):
        d = self._dir(_spin_independent(_DIAG))
        with open(os.path.join(d, "transfer.dat"), "w") as f:
            f.write("Transfer w90spin/1\n{}\n1\n 1\n".format(header_norb))
            for (r, a, b, s_, t), v in _closed(entries).items():
                f.write("{} {} {} {} {} {} {} {:.12f} {:.12f}\n".format(
                    *r, a, b, s_, t, v.real, v.imag))
        return d

    def test_line_two_smaller_than_geometry_refused(self):
        # the review's case: line 2 = 1, Geometry norb = 2; without the
        # check the down spin of orbital 1 is read as the up spin of 2
        d = self._dir_header(1, {((1, 0, 0), 1, 1, 0, 0): -1.0,
                                 ((1, 0, 0), 1, 1, 1, 1): -1.0})
        for scheme in ("general", "auto"):
            with self.subTest(scheme=scheme), self.assertRaisesRegex(
                    ValueError, r"gives norb = 1, but the Geometry norb is 2 "
                                r"\(RPA/FLEX\)"):
                _rpa(d, scheme=scheme)

    def test_line_two_larger_than_geometry_refused(self):
        d = self._dir_header(4, _spin_independent(_DIAG))
        with self.assertRaisesRegex(ValueError, "gives norb = 4, but the "
                                    "Geometry norb is 2"):
            _rpa(d, scheme="general")

    def test_flex_refuses_too(self):
        import hwave.solver.flex as flex_module
        d = self._dir_header(1, {((1, 0, 0), 1, 1, 0, 0): -1.0,
                                 ((1, 0, 0), 1, 1, 1, 1): -1.0})
        ham, _ = _ham(d)
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=1, Mix=1.0,
                   flex_hartree_fock=True)
        with self.assertRaisesRegex(ValueError, "gives norb = 1"):
            flex_module.FLEX(ham, {}, {"mode": "FLEX", "param": par,
                                       "calc_scheme": "general",
                                       "calc_type": "ring",
                                       "enable_spin_orbital": False})


class TestSpinResolvedChecks(_Fixture):
    """The checks spin-orbital mode applies to its spinful one-body term
    apply to a spin-resolved Transfer as well."""

    def test_coeff_extern_is_validated(self):
        d = self._dir(_soc_entries(_DIAG),
                      extern_rows=["0 0 0 1 1 0.3 0.0"])
        with self.assertRaisesRegex(ValueError, "coeff_extern"):
            _rpa(d, scheme="general", extern=True, coeff_extern=True)

    def test_chi0q_init_needs_the_index_convention(self):
        d = self._dir(_soc_entries(_DIAG))
        _run_chi0q(d, scheme="general")
        src = os.path.join(d, "output", "chi0q.npz")
        with np.load(src) as z:
            members = {k: z[k] for k in z.files if k != "index_convention"}
        np.savez(os.path.join(d, "output", "no_conv.npz"), **members)
        s, _ = _rpa(d, scheme="general")
        with self.assertRaisesRegex(ValueError, "lacks an 'index_convention'"):
            s.read_init({"path_to_input": os.path.join(d, "output"),
                         "chi0q_init": "no_conv.npz"})

    def test_spin_orbital_archive_gets_the_specific_message(self):
        dso = self._dir(interleaved=_soc_entries(_DIAG), geom_norb=4)
        _run_chi0q(dso, scheme="general", enable_spin_orbital=True)
        s, _ = _rpa(self._dir(_soc_entries(_DIAG)), scheme="general")
        with self.assertRaisesRegex(ValueError, "produced by a spin-orbital run"):
            s.read_init({"path_to_input": os.path.join(dso, "output"),
                         "chi0q_init": "chi0q.npz"})

    def test_malformed_spin_orbital_marker_refused(self):
        d = self._dir(_soc_entries(_DIAG))
        _run_chi0q(d, scheme="general")
        with np.load(os.path.join(d, "output", "chi0q.npz")) as z:
            members = {k: z[k] for k in z.files}
        members["spin_orbital"] = np.int64(0)
        np.savez(os.path.join(d, "output", "bad_mark.npz"), **members)
        s, _ = _rpa(d, scheme="general")
        with self.assertRaisesRegex(ValueError, "malformed spin_orbital marker"):
            s.read_init({"path_to_input": os.path.join(d, "output"),
                         "chi0q_init": "bad_mark.npz"})

    def test_mismatched_or_non_finite_entries_refused(self):
        d = self._dir(_spin_independent(_DIAG))
        with open(os.path.join(d, "transfer.dat"), "a") as f:
            f.write("0 0 0 1 2 0 1 0.150000000000 0.100000000000\n")
            f.write("0 0 0 2 1 1 0 0.300000000000 -0.100000000000\n")
        with self.assertRaisesRegex(ValueError, r"not Hermitian-closed.*expected its "
                                    r"conjugate, \(0\.15-0\.1j\)"):
            _rpa(d, scheme="general")
        d = self._dir(_spin_independent(_DIAG))
        with open(os.path.join(d, "transfer.dat"), "a") as f:
            f.write("0 0 0 1 2 0 1 nan 0.0\n0 0 0 2 1 1 0 nan 0.0\n")
        with self.assertRaisesRegex(ValueError, "non-finite coefficient"):
            _rpa(d, scheme="general")

    def test_missing_spin_flip_partner_refused(self):
        d = self._dir(_spin_independent(_DIAG))
        with open(os.path.join(d, "transfer.dat"), "a") as f:
            f.write("0 0 0 1 2 0 1 0.150000000000 0.100000000000\n")
        with self.assertRaisesRegex(
                ValueError, r"NO partner entry at R=\(0, 0, 0\) \(orbital, spin\) "
                            r"= \(2, 1\), \(1, 0\)"):
            _rpa(d, scheme="general")


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

    def test_bond_channels_are_refused(self):
        # the bond-resolved spinful FLEX (#270) is built from the
        # spin-orbital entry context, which only enable_spin_orbital makes
        with self.assertRaisesRegex(
                ValueError, r"spin-resolved \(9-column\) Transfer \(FLEX\) does not "
                            r"support longitudinal_bond_channels"):
            self._flex(self._dir(_soc_entries()), "general",
                       {"flex_hartree_fock": True, "longitudinal_bond_channels": True})

    def test_offsite_refusal_does_not_point_to_the_bond_channels(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_module
        d = self._dir(_soc_entries())
        with open(os.path.join(d, "coulombinter.dat"), "w") as f:
            f.write("CoulombInter\n2\n1\n 1\n")
            f.write("1 0 0 1 1 0.2 0.0\n-1 0 0 1 1 0.2 0.0\n")
        files = dict(_FILES, path_to_input=d, CoulombInter="coulombinter.dat")
        read_io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": files})
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=1, Mix=1.0,
                   flex_hartree_fock=True)
        s = flex_module.FLEX(read_io.get_param("ham"), {},
                             {"mode": "FLEX", "param": par, "calc_scheme": "general",
                              "calc_type": "ring", "enable_spin_orbital": False})
        with self.assertRaisesRegex(
                ValueError, r"spin-resolved \(9-column\) Transfer \(FLEX\) supports "
                            r"on-site and in-cell interactions only \(.*off-site "
                            r"'CoulombInter' declaration .* is refused\.$"):
            s.solve(read_io.get_param("green"), d)


_NUMPY_BELOW_2 = int(np.__version__.split(".")[0]) < 2

_TUTORIAL_1ORB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "docs", "tutorial", "Hubbard", "FLEX", "1orb")


def _nine_column_copy(src, dst):
    """Rewrite the 7-column Transfer file ``src`` as the equivalent
    spin-resolved file (each row once per spin, s = t)."""
    with open(src) as f:
        lines = f.read().splitlines()
    nr = int(lines[2].split()[0])
    body0 = 3 + (nr + 14) // 15
    out = [lines[0] + " w90spin/1"] + lines[1:body0]
    for line in lines[body0:]:
        v = line.split()
        if not v:
            break
        for spin in (0, 1):
            out.append(" ".join(v[:5] + [str(spin), str(spin)] + v[5:7]))
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


class TestFlexSpinFreeNineColumnTutorial(unittest.TestCase):
    """The FLEX 1orb tutorial with its Transfer rewritten as 9 columns.
    The spin structure is judged from the matrix: a spin-independent file
    takes the 7-column path, so the DEFAULT calc_scheme = "auto" runs (it
    resolves as for the 7-column file) and every archive equals the
    7-column run's (bitwise with numpy >= 2, to round-off below); so does
    calc_scheme = "general" with flex_hartree_fock."""

    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="w90spin_flex_tut_")
        for name in ("geom.dat", "transfer.dat", "coulombintra.dat"):
            shutil.copy(os.path.join(_TUTORIAL_1ORB, name), self.d)
        _nine_column_copy(os.path.join(self.d, "transfer.dat"),
                          os.path.join(self.d, "transfer9.dat"))

    def tearDown(self):
        shutil.rmtree(self.d, ignore_errors=True)

    def _run(self, transfer, tag, scheme=None, hf=None):
        import tomli
        import hwave.qlms as qlms
        with open(os.path.join(_TUTORIAL_1ORB, "input.toml"), "rb") as f:
            inp = tomli.load(f)
        inp["log"] = {"print_level": 0}
        inp["mode"]["param"]["IterationMax"] = 20
        if scheme is not None:
            inp["mode"]["calc_scheme"] = scheme
        if hf is not None:
            inp["mode"]["param"]["flex_hartree_fock"] = hf
        inter = inp["file"]["input"]["interaction"]
        inter["path_to_input"] = self.d
        inter["Transfer"] = transfer
        inp["file"]["input"]["path_to_input"] = self.d
        out = os.path.join(self.d, "out_" + tag)
        inp["file"]["output"]["path_to_output"] = out
        qlms.run(input_dict=inp)
        return out

    def _assert_same_archives(self, a, b):
        for name in ("chi0q.npz", "chiq.npz", "sigma.npz", "green.npz"):
            with np.load(os.path.join(a, name)) as za, np.load(os.path.join(b, name)) as zb:
                self.assertEqual(set(za.files), set(zb.files), name)
                for key in za.files:
                    if key == "producer_version":
                        continue
                    x, y = za[key], zb[key]
                    self.assertEqual((x.dtype, x.shape), (y.dtype, y.shape), (name, key))
                    if _NUMPY_BELOW_2 and x.dtype.kind in "fc":
                        # numpy < 2: two runs of one input already differ
                        # at round-off on this platform, so only the
                        # round-off level is meaningful there
                        self.assertTrue(np.all(np.isfinite(x)) and np.all(np.isfinite(y)),
                                        (name, key))
                        np.testing.assert_allclose(x, y, rtol=1e-12, atol=1e-14,
                                                   equal_nan=False,
                                                   err_msg=str((name, key)))
                    else:
                        self.assertTrue(np.array_equal(x, y), (name, key))
                self.assertNotIn("spin_resolved", za.files, name)

    def test_default_auto_runs_and_equals_the_seven_column_run(self):
        self.assertIsNotNone(wan90.read_w90_header_marker(
            os.path.join(self.d, "transfer9.dat")))
        self._assert_same_archives(self._run("transfer9.dat", "9auto"),
                                   self._run("transfer.dat", "7auto"))

    def test_general_with_hartree_fock_equals_the_seven_column_run(self):
        self._assert_same_archives(
            self._run("transfer9.dat", "9gen", "general", True),
            self._run("transfer.dat", "7gen", "general", True))

    def test_the_caller_table_is_not_rewritten(self):
        # the 7-column view lives in FLEX's private container; the reader
        # result keeps its spin-block table for another solver
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_module
        files = {"path_to_input": self.d, "Geometry": "geom.dat",
                 "Transfer": "transfer9.dat", "CoulombIntra": "coulombintra.dat"}
        ham = read_input_k.QLMSkInput({"path_to_input": self.d,
                                       "interaction": files}).get_param("ham")
        before = dict(ham["Transfer"])
        par = {"T": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1], "Nmat": 16,
               "filling": 0.5, "IterationMax": 1}
        s = flex_module.FLEX(ham, {}, {"mode": "FLEX", "param": par})
        self.assertTrue(s.ham_info.spin_free_collapsed)
        self.assertFalse(s.ham_info.spinful_one_body)
        self.assertEqual(s.calc_scheme, "reduced")
        self.assertEqual(dict(ham["Transfer"]), before)
        self.assertEqual(ham.source_formats["transfer"], ("w90spin", 1))
        # RPA keeps the nd x nd form and judges the matrix itself
        import hwave.solver.rpa as rpa_module
        r = rpa_module.RPA(ham, {}, {"mode": "RPA", "param": par})
        self.assertTrue(r.ham_info.spin_resolved_transfer)


def _assert_same_archive_dirs(test, a, b, names):
    """Every member of each archive equal between two output directories
    (bitwise with numpy >= 2, to round-off below; NaN never passes)."""
    for name in names:
        with np.load(os.path.join(a, name)) as za, np.load(os.path.join(b, name)) as zb:
            test.assertEqual(set(za.files), set(zb.files), name)
            for key in za.files:
                if key == "producer_version":
                    continue
                x, y = za[key], zb[key]
                test.assertEqual((x.dtype, x.shape), (y.dtype, y.shape), (name, key))
                if _NUMPY_BELOW_2 and x.dtype.kind in "fc":
                    test.assertTrue(np.all(np.isfinite(x)) and np.all(np.isfinite(y)),
                                    (name, key))
                    np.testing.assert_allclose(x, y, rtol=1e-12, atol=1e-14,
                                               equal_nan=False, err_msg=str((name, key)))
                else:
                    test.assertTrue(np.array_equal(x, y), (name, key))


class TestFlexSpinFreeNineColumnWithSublattice(_Fixture):
    """A spin-independent 9-column Transfer read by FLEX as 7 columns, with
    SubShape != [1, 1, 1]: the private table container is the one the
    sublattice fold rewrites, so the reciprocal vectors (and the archives'
    wavevector_unit) equal the 7-column run's."""

    def _run(self, d, tag):
        import hwave.qlms as qlms
        out = os.path.join(d, "out_" + tag)
        par = dict(_PARAM, CellShape=[4, 1, 1], SubShape=[2, 1, 1],
                   IterationMax=3, Mix=0.5)
        inp = {"log": {"print_level": 0},
               "mode": {"mode": "FLEX", "param": par, "calc_type": "ring"},
               "file": {"input": {"path_to_input": d,
                                  "interaction": dict(_FILES, path_to_input=d)},
                        "output": {"path_to_output": out, "chi0q": "chi0q.npz",
                                   "chiq": "chiq.npz", "green": "green",
                                   "sigma": "sigma"}}}
        qlms.run(input_dict=inp)
        return out

    def test_wavevectors_and_archives_equal_the_seven_column_run(self):
        d7 = self._dir(plain=_DIAG)
        d9 = self._dir(_spin_independent(_DIAG))
        o7, o9 = self._run(d7, "7"), self._run(d9, "9")
        with np.load(os.path.join(o9, "chi0q.npz")) as z9, \
                np.load(os.path.join(o7, "chi0q.npz")) as z7:
            self.assertIn("wavevector_unit", z9.files)
            np.testing.assert_array_equal(z9["wavevector_unit"], z7["wavevector_unit"])
        _assert_same_archive_dirs(self, o9, o7, ("chi0q.npz", "chiq.npz",
                                                 "sigma.npz", "green.npz"))


class TestFlexAutoOnASpinDependentTransfer(_Fixture):
    def test_auto_resolves_to_general(self):
        # a spin flip: only the general scheme supports it, so auto goes
        # there (and then asks for the spin-orbital conditions)
        # CoulombIntra only: the interaction types alone resolve to reduced
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_module
        d = self._dir(_soc_entries(_DIAG))

        def ham():
            files = {"path_to_input": d, "Geometry": "geom.dat",
                     "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
            return read_input_k.QLMSkInput({"path_to_input": d,
                                            "interaction": files}).get_param("ham")
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=1, Mix=1.0,
                   flex_hartree_fock=True)
        s = flex_module.FLEX(ham(), {}, {"mode": "FLEX", "param": par,
                                         "calc_type": "ring"})
        self.assertEqual(s.calc_scheme, "general")
        self.assertEqual(s._scheme_resolution, "auto:mixed:transfer")
        self.assertTrue(s._so_spinful)
        self.assertFalse(s.ham_info.spin_free_collapsed)
        par.pop("flex_hartree_fock")
        with self.assertRaisesRegex(ValueError, "flex_hartree_fock"):
            flex_module.FLEX(ham(), {}, {"mode": "FLEX", "param": par,
                                         "calc_type": "ring"})

    def test_auto_resolves_to_general_in_spin_orbital_mode_too(self):
        import hwave.qlmsio.read_input_k as read_input_k
        import hwave.solver.flex as flex_module
        d = self._dir(interleaved=_soc_entries(_DIAG), geom_norb=4)
        files = {"path_to_input": d, "Geometry": "geom.dat",
                 "Transfer": "transfer.dat", "CoulombIntra": "coulombintra.dat"}
        ham = read_input_k.QLMSkInput({"path_to_input": d,
                                       "interaction": files}).get_param("ham")
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=1, Mix=1.0,
                   flex_hartree_fock=True)
        s = flex_module.FLEX(ham, {}, {"mode": "FLEX", "param": par,
                                       "calc_type": "ring",
                                       "enable_spin_orbital": True})
        self.assertEqual((s.calc_scheme, s._scheme_resolution),
                         ("general", "auto:mixed:transfer"))


class TestFlexRun(_Fixture):
    def test_run_writes_archives_without_spin_orbital_provenance(self):
        # a spin-resolved run carries no spin-orbital fingerprints (they
        # belong to the entry context of enable_spin_orbital), in FLEX as
        # in RPA; the run itself goes through the spinful route to the end
        import hwave.qlms as qlms
        d = self._dir(_soc_entries())
        out = os.path.join(d, "output")
        par = dict(_PARAM, SubShape=[1, 1, 1], IterationMax=2, Mix=0.5,
                   flex_hartree_fock=True)
        inp = {"log": {"print_level": 0},
               "mode": {"mode": "FLEX", "param": par, "calc_scheme": "general",
                        "calc_type": "ring", "enable_spin_orbital": False},
               "file": {"input": {"path_to_input": d,
                                  "interaction": dict(_FILES, path_to_input=d)},
                        "output": {"path_to_output": out, "chi0q": "chi0q.npz",
                                   "chiq": "chiq.npz", "green": "green",
                                   "sigma": "sigma"}}}
        qlms.run(input_dict=inp)
        for name in ("chi0q.npz", "chiq.npz", "green.npz", "sigma.npz"):
            with np.load(os.path.join(out, name)) as z:
                stray = [k for k in z.files if k.startswith(("fp_", "prov_"))]
                self.assertEqual(stray, [], name)
                # not a spin-orbital archive: the marker of this format is
                # spin_resolved, on chi0q / chiq as in RPA
                self.assertNotIn("spin_orbital", z.files, name)
                if name in ("chi0q.npz", "chiq.npz"):
                    self.assertIn("spin_resolved", z.files, name)
                    self.assertTrue(bool(z["spin_resolved"]), name)
        # the chi0q is admitted as chi0q_init by a spin-resolved run and
        # refused by a spin-orbital one
        s, _ = _rpa(d, scheme="general")
        info = s.read_init({"path_to_input": out, "chi0q_init": "chi0q.npz"})
        self.assertIn("chi0q", info)
        dso = self._dir(interleaved=_soc_entries(), geom_norb=4)
        sso, _ = _rpa(dso, scheme="general", enable_spin_orbital=True)
        with self.assertRaisesRegex(ValueError, "spin_resolved = True"):
            sso.read_init({"path_to_input": out, "chi0q_init": "chi0q.npz"})


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
