"""The spin-orbital FLEX: ``[mode] enable_spin_orbital = true`` with
``calc_scheme = "general"`` and ``flex_hartree_fock = true`` (issue #83).

The spin is folded into the orbital index (``nd = 2 norb``, spin-block order
``s * norb + a``); the first order is the self-consistent Hartree-Fock term
of the full spin-orbital density (``flex_hf.hf_map_spin_orbital``) and the
fluctuation part is one particle-hole channel built from the antisymmetrized
spinful vertex ``Gamma``::

    chi = [1 + chi0 Gamma]^-1 chi0
    V   = P[ Gamma (chi - chi0 / 2) Gamma ]

with ``P`` the intra-pair transpose of both pairs (the conserving
orientation of issue #251).

What this module pins:

* :class:`TestConstructionGates` -- the admitted configuration and every
  refusal, each with a message that names what to change.
* :class:`TestSpinConservingLimit` -- a spin-diagonal model written in the
  spin-orbital encoding converges to the SAME fixed point as the spin-free
  FLEX (static and fluctuation self-energy, mu, N), for the rotation-
  invariant on-site set. ``PairLift`` is excluded on purpose: it is not
  spin conserving (``Delta S_z = 2``), and the spin-free FLEX drops its
  transverse content.
* :class:`TestPhiDerivability` -- the second- and third-order terms of the
  production kernel are the functional derivatives of ``Tr (chi0 Gamma)^n``
  (direction-independent ratio on a random complex Green function), and
  the same check FAILS without ``P``.
* :class:`TestFirstOrder` -- ``hf_map_spin_orbital`` equals the test-side
  Wick functional of the exact-diagonalization gate
  (``tests/test_flex_second_order_ed_onsite.wick_hf_sigma``) for every type,
  on generic (spin-mixing) Hermitian densities.
* :class:`TestOutputs` -- the archive tags of a short run and the spin
  vector in ``physics``.

The exactness of the second order against exact diagonalization with a
spin-mixing one-body term is in ``tests/test_flex_spin_orbital_ed.py``.
"""
import contextlib
import logging
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k

#: physical orbitals of the lattice fixture
NP = 2

#: the rotation-invariant on-site set: ``J_Exchange = J_Hund - 2 J_Ising``
_J_HUND, _J_ISING = 0.15, 0.05
_ONSITE_VALUES = (("Hund", _J_HUND), ("Exchange", _J_HUND - 2 * _J_ISING),
                  ("PairHop", 0.07), ("Ising", _J_ISING), ("PairLift", 0.03))
_SU2_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Exchange", "Ising", "PairHop")


@contextlib.contextmanager
def _quiet():
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _wfile(path, head, n, rows):
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(head, n))
        for r in rows:
            f.write("  %d %d %d %d %d  %.12f  %.12f\n" % r)


def _write_model(d, so, types=_SU2_TYPES, soc=0.0, offsite=False, seed=3):
    """A two-orbital chain (``CellShape = [4, 1, 1]``) with complex hopping.

    ``so``: write it in the spin-orbital encoding (geometry of ``2 NP``
    orbitals, transfer on the interleaved file index ``2 a + s + 1``); the
    interaction files use the PHYSICAL orbital index either way. ``soc``
    adds an on-site spin-mixing block (spin-orbital encoding only).
    ``offsite`` adds a nearest-neighbour CoulombInter bond."""
    rng = np.random.default_rng(seed)
    n = 2 * NP if so else NP
    with open(os.path.join(d, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % n)
        f.write("0 0 0\n" * n)
    H0 = rng.normal(size=(NP, NP)) + 1j * rng.normal(size=(NP, NP))
    H0 = 0.3 * (H0 + H0.conj().T)
    T = 0.4 * (rng.normal(size=(NP, NP)) + 1j * rng.normal(size=(NP, NP))) + np.eye(NP)
    F = soc * (rng.normal(size=(NP, NP)) + 1j * rng.normal(size=(NP, NP)))
    rows = []
    for a in range(NP):
        for b in range(NP):
            for R, v in ((0, H0[a, b]), (1, T[a, b]), (-1, T[b, a].conjugate())):
                if so:
                    for sp in range(2):
                        rows.append((R, 0, 0, 2 * a + sp + 1, 2 * b + sp + 1, v.real, v.imag))
                else:
                    rows.append((R, 0, 0, a + 1, b + 1, v.real, v.imag))
            if so and soc:
                f_ab = F[a, b]
                rows.append((0, 0, 0, 2 * a + 1, 2 * b + 2, f_ab.real, f_ab.imag))
                rows.append((0, 0, 0, 2 * b + 2, 2 * a + 1, f_ab.real, -f_ab.imag))
    _wfile(os.path.join(d, "transfer.dat"), "Transfer", n, rows)
    offd = [(a, b) for a in range(NP) for b in range(NP) if a != b]
    spec = {"CoulombIntra": [(0, a, a, 1.0 + 0.2 * a) for a in range(NP)],
            "CoulombInter": [(0, a, b, 0.5) for a, b in offd]}
    for nm, v in _ONSITE_VALUES:
        spec[nm] = [(0, a, b, v) for a, b in offd]
    if offsite:
        spec["CoulombInter"] = spec["CoulombInter"] + [(1, 0, 0, 0.2), (-1, 0, 0, 0.2)]
    files = {}
    for nm in types:
        fn = nm.lower() + ".dat"
        _wfile(os.path.join(d, fn), nm, NP,
               [(R, 0, 0, a + 1, b + 1, v, 0.0) for R, a, b, v in spec[nm]])
        files[nm] = fn
    return files


def _build(d, files, so=True, scheme="general", nmat=16, **par_extra):
    import hwave.solver.flex as flex_mod
    inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    inter.update(files)
    io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
    par = {"T": 0.5, "filling": 0.4, "CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
           "Nmat": nmat, "IterationMax": 100, "Mix": 0.5, "EPS": 12,
           "mixing_scheme": "linear", "flex_hartree_fock": True}
    par.update(par_extra)
    with _quiet():
        s = flex_mod.FLEX(io.get_param("ham"), {}, {
            "mode": "FLEX", "param": par, "enable_spin_orbital": so,
            "calc_scheme": scheme})
    return s, io.get_param("green")


class _TmpDir(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp(prefix="flex_so_")
        self.addCleanup(shutil.rmtree, self.d, True)


class TestConstructionGates(_TmpDir):
    """What the spin-orbital FLEX admits, and the message of every refusal."""

    def test_admitted_configuration(self):
        s, _ = _build(self.d, _write_model(self.d, True, soc=0.2))
        self.assertTrue(s._flex_general)
        self.assertTrue(s._so_spinful)
        self.assertEqual(s.norb, NP)

    def test_requires_flex_hartree_fock(self):
        files = _write_model(self.d, True)
        with self.assertRaisesRegex(ValueError, r"enable_spin_orbital = true \(FLEX\) requires "
                                    r"\[mode\.param\] flex_hartree_fock = true"):
            _build(self.d, files, flex_hartree_fock=False)

    def test_refuses_longitudinal_bond_channels(self):
        files = _write_model(self.d, True)
        with self.assertRaisesRegex(ValueError, r"enable_spin_orbital = true \(FLEX\) does not "
                                    r"support longitudinal_bond_channels"):
            _build(self.d, files, longitudinal_bond_channels=True)

    def test_requires_spinful_vertex_exchange(self):
        files = _write_model(self.d, True)
        with self.assertRaisesRegex(ValueError, r"requires spinful_vertex_exchange = true"):
            _build(self.d, files, spinful_vertex_exchange=False)

    def test_reduced_scheme_still_refuses(self):
        files = _write_model(self.d, True, types=("CoulombIntra", "CoulombInter", "Hund"))
        with self.assertRaisesRegex(ValueError, r"calc_scheme='reduced' FLEX does not "
                                    r"support enable_spin_orbital"):
            _build(self.d, files, scheme="reduced", flex_hartree_fock=False)

    def test_refuses_offsite_interactions(self):
        files = _write_model(self.d, True, offsite=True)
        s, gi = _build(self.d, files)
        with self.assertRaisesRegex(ValueError, r"supports on-site interactions only .*"
                                    r"off-site 'CoulombInter' declaration"):
            with _quiet():
                s.solve(gi, self.d)

    def test_refuses_sigma_init(self):
        """Through the public route: the split ``sigma.npz`` of an earlier
        spin-orbital run, read by ``read_init``."""
        files = _write_model(self.d, True, soc=0.2)
        s, gi = _build(self.d, files, IterationMax=1)
        with _quiet():
            s.solve(gi, self.d)
            s.save_results({"path_to_output": self.d, "sigma": "sigma.npz"}, gi)
        s2, _ = _build(self.d, files)
        with _quiet():
            gi2 = s2.read_init({"path_to_input": self.d, "sigma_init": "sigma.npz"})
        with self.assertRaisesRegex(ValueError, r"does not accept sigma_init in this version"):
            with _quiet():
                s2.solve(gi2, self.d)

    def test_refuses_user_mean_field(self):
        """``trans_mod`` / ``green_init`` as ``read_init`` places them in
        the green_info dictionary."""
        for key in ("trans_mod", "green_init", "sigma_init"):
            with self.subTest(seed=key):
                s, gi = _build(self.d, _write_model(self.d, True))
                gi[key] = np.zeros(1)
                with self.assertRaisesRegex(ValueError, r"does not accept {} in this version"
                                            .format(key)):
                    with _quiet():
                        s.solve(gi, self.d)


class TestOutputPaths(_TmpDir):
    def test_single_channel_may_take_a_spin_charge_file_name(self):
        """The spin-orbital run writes ``chiq`` only, so its preflight must
        not reserve the ``chiq_s`` / ``chiq_c`` names it never writes."""
        s, gi = _build(self.d, _write_model(self.d, True), IterationMax=1)
        out = {"path_to_output": self.d, "chiq": "chiq_s.npz"}
        s.validate_output_paths(out, self.d)
        with _quiet():
            s.solve(gi, self.d)
            s.save_results(out, gi)
        with np.load(os.path.join(self.d, "chiq_s.npz")) as z:
            self.assertIn("chiq", z.files)
            self.assertEqual(str(z["chi_convention"]), "spin_orbital")
        self.assertIsNone(s._so_vertex)           # solve-scoped, released


class TestSpinConservingLimit(unittest.TestCase):
    """A spin-diagonal model in the spin-orbital encoding converges to the
    fixed point of the spin-free FLEX, block for block."""

    @classmethod
    def setUpClass(cls):
        cls.dirs = [tempfile.mkdtemp(prefix="flex_so_limit_") for _ in range(2)]
        cls.res = {}
        for so, d in zip((False, True), cls.dirs):
            s, gi = _build(d, _write_model(d, so), so=so, nmat=32,
                           flex_second_order="local")
            with _quiet():
                s.solve(gi, d)
            cls.res[so] = (s, gi)

    @classmethod
    def tearDownClass(cls):
        for d in cls.dirs:
            shutil.rmtree(d, ignore_errors=True)

    def test_converged(self):
        for so in (False, True):
            with self.subTest(spin_orbital=so):
                self.assertTrue(self.res[so][0].scf_converged)

    def test_fixed_point_equals_the_spin_free_one(self):
        (a, ga), (b, gb) = self.res[False], self.res[True]
        for key in ("sigma_static", "sigma_fluct"):
            x, y = np.asarray(ga[key]), np.asarray(gb[key])
            blk = np.zeros(x.shape[:-2] + (2 * NP, 2 * NP), dtype=complex)
            blk[..., :NP, :NP] = x
            blk[..., NP:, NP:] = x
            with self.subTest(key=key):
                self.assertGreater(np.abs(x).max(), 1e-2)        # anti-vacuity
                np.testing.assert_allclose(y, blk, rtol=0, atol=1e-10)
        self.assertAlmostEqual(a.mu, b.mu, delta=1e-10)
        self.assertAlmostEqual(ga["physics"]["NCond"], gb["physics"]["NCond"], delta=1e-10)
        for k in ("Sx", "Sy", "Sz"):
            self.assertAlmostEqual(gb["physics"][k], 0.0, delta=1e-10)


class TestPhiDerivability(unittest.TestCase):
    """The production spinful kernel is conserving order by order.

    Scaling the cached vertex ``Gamma -> x Gamma`` inside one production
    assembly (``_calc_chi0q`` -> ``_flex_compute_veff_spinful`` ->
    ``_calc_self_energy_general``), the part even in ``x`` starts with the
    second order and the odd part with the third; a two-level Richardson
    step removes the next order of each. The isolated ``Sigma_n`` must be
    proportional to ``dPhi_n / dG`` with ``Phi_n = sum Tr (chi0 Gamma)^n``,
    i.e. ``dPhi_n / Tr[Sigma_n dG]`` must not depend on the direction ``dG``
    of a random complex Green function. Without the pair transpose ``P`` it
    does (an O(1) spread)."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp(prefix="flex_so_phi_")
        types = ("CoulombIntra", "CoulombInter", "Hund", "Exchange", "PairHop", "PairLift")
        s, _ = _build(cls.d, _write_model(cls.d, True, types=types, soc=0.3), nmat=8)
        with _quiet():
            s._calc_epsilon_k({})
        cls.s = s
        cls.beta = 1.0 / s.T
        cls.nvol = s.lattice.nvol
        cls.no = 2 * NP
        cls.ndx = cls.no * cls.no
        cls.nmat = s.nmat
        cls.gam0 = np.array(s._spinful_vertex(np), copy=True)
        rng = np.random.default_rng(1)
        shape = (1, cls.nmat, cls.nvol, cls.no, cls.no)
        cls.G = 0.3 * (rng.normal(size=shape) + 1j * rng.normal(size=shape))
        cls.dG = [rng.normal(size=shape) + 1j * rng.normal(size=shape) for _ in range(4)]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def _sigma(self, x):
        s = self.s
        s._so_vertex = x * self.gam0
        try:
            c0 = s._calc_chi0q(self.G, np.zeros_like(self.G), self.beta)[0]
            _, v, _ = s._flex_compute_veff_spinful(c0)
            return np.asarray(s._calc_self_energy_general(self.G, v, self.beta))
        finally:
            s._so_vertex = None

    def _phi(self, G, n):
        c = self.s._calc_chi0q(G, np.zeros_like(G), self.beta)[0]
        m = c.reshape(self.nmat, self.nvol, self.ndx, self.ndx) @ self.gam0[None]
        p = m.copy()
        for _ in range(n - 1):
            p = p @ m
        return np.trace(p, axis1=-2, axis2=-1).sum()

    def _sigma_n(self):
        h = 0.02
        sp, sm = self._sigma(h), self._sigma(-h)
        sp2, sm2 = self._sigma(2 * h), self._sigma(-2 * h)
        even1, even2 = 0.5 * (sp + sm), 0.5 * (sp2 + sm2)
        odd1, odd2 = 0.5 * (sp - sm), 0.5 * (sp2 - sm2)
        return {2: (16.0 * even1 - even2) / (12.0 * h ** 2),        # x^4 removed
                3: (32.0 * odd1 - odd2) / (24.0 * h ** 3)}          # x^5 removed

    def _spread(self, sigma, n):
        r = []
        for dG in self.dG:
            e = 1e-6
            dphi = (self._phi(self.G + e * dG, n) - self._phi(self.G - e * dG, n)) / (2 * e)
            r.append(dphi / np.einsum("gfkmn,gfknm->", sigma, dG))
        r = np.array(r)
        return float(np.max(np.abs(r - r.mean())) / abs(r.mean()))

    def test_second_and_third_order_are_phi_derivable(self):
        sig = self._sigma_n()
        for n in (2, 3):
            with self.subTest(order=n):
                self.assertGreater(float(np.abs(sig[n]).max()), 1e-3)     # anti-vacuity
                self.assertLess(self._spread(sig[n], n), 1e-5)

    def test_untransposed_kernel_is_not(self):
        """Mutation: the same assembly with ``P`` replaced by a plain copy."""
        import hwave.solver.flex as flex_mod
        orig = flex_mod._pair_transpose_into

        def _copy(out, v, no):
            out[...] = v
        flex_mod._pair_transpose_into = _copy
        try:
            sig = self._sigma_n()
        finally:
            flex_mod._pair_transpose_into = orig
        for n in (2, 3):
            with self.subTest(order=n):
                self.assertGreater(self._spread(sig[n], n), 1e-2)


class TestFirstOrder(unittest.TestCase):
    """``flex_hf.hf_map_spin_orbital`` is the Hartree-Fock functional of the
    declared Hamiltonian on the FULL spin-orbital density, spin-off-diagonal
    blocks included: on the one-site three-orbital model of the
    exact-diagonalization gate it equals the test-side Wick functional
    ``wick_hf_sigma`` of the ED monomials (with the mirrored-row half weight
    of that module) for every type and a mixed declaration."""

    def _tables(self, d, terms):
        """The production Hartree-Fock tables of a spin-orbital FLEX solver
        built on the one-site model with the declaration ``terms``."""
        from tests import test_flex_second_order_ed_onsite as G3
        from tests.test_second_order_factors import _write_wan
        from hwave.solver import flex_hf
        import hwave.solver.flex as flex_mod
        n = G3.NORB
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(2 * n) + "0 0 0\n" * (2 * n))
        _write_wan(os.path.join(d, "transfer.dat"), "Transfer", 2 * n,
                   [(0, 0, 0, p + 1, p + 1, 0.1, 0.0) for p in range(2 * n)])
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for itype, v in terms:
            _write_wan(os.path.join(d, itype + ".dat"), itype, n, G3._rows(itype, v))
            idict[itype] = itype + ".dat"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        par = {"T": 0.5, "mu": 0.0, "CellShape": [1, 1, 1], "SubShape": [1, 1, 1],
               "Nmat": 4, "IterationMax": 1, "flex_hartree_fock": True}
        with _quiet():
            s = flex_mod.FLEX(r.get_param("ham"), {}, {
                "mode": "FLEX", "param": par, "enable_spin_orbital": True,
                "calc_scheme": "general"})
        self.assertEqual(s.norb, n)
        return flex_hf.build_flex_hf_tables(s.ham_info.param_ham, s.norb, (1, 1, 1))

    def test_equals_the_wick_functional(self):
        from tests import test_flex_second_order_ed_onsite as G3
        from hwave.solver import flex_hf
        n = G3.NORB
        cases = [[(t, 0.41)] for t in G3._TYPES]
        cases.append([("CoulombIntra", 0.8), ("Hund", 0.3), ("Exchange", 0.2),
                      ("PairHop", 0.15), ("PairLift", 0.1), ("Ising", 0.05)])
        dens = G3._random_densities(spin_block_diagonal=False)
        for terms in cases:
            d = tempfile.mkdtemp(prefix="flex_so_hf_")
            try:
                tables = self._tables(d, terms)
            finally:
                shutil.rmtree(d, ignore_errors=True)
            for k, rho in enumerate(dens):
                mine = flex_hf.hf_map_spin_orbital(rho[None], tables, (1, 1, 1), n)[0]
                ref = G3._mode_to_gen(G3.wick_hf_sigma(2 * n, G3._mode_terms(terms),
                                                       G3._gen_to_mode(rho)))
                with self.subTest(terms=[t for t, _ in terms], density=k):
                    self.assertGreater(np.abs(ref).max(), 1e-3)      # anti-vacuity
                    self.assertLess(np.abs(mine - ref).max(), 1e-12)

    def test_spin_off_diagonal_density_is_load_bearing(self):
        """PairLift and Exchange carry a spin-off-diagonal first order: zeroing
        the off-diagonal spin blocks of the density changes the result, so
        the comparison above is not blind to them."""
        from tests import test_flex_second_order_ed_onsite as G3
        from hwave.solver import flex_hf
        n = G3.NORB
        rho = G3._random_densities(spin_block_diagonal=False)[0]
        rho_bd = rho.copy()
        rho_bd[:n, n:] = 0.0
        rho_bd[n:, :n] = 0.0
        for t in ("PairLift", "Exchange"):
            d = tempfile.mkdtemp(prefix="flex_so_hf_")
            try:
                tables = self._tables(d, [(t, 0.41)])
            finally:
                shutil.rmtree(d, ignore_errors=True)
            a = flex_hf.hf_map_spin_orbital(rho[None], tables, (1, 1, 1), n)[0]
            b = flex_hf.hf_map_spin_orbital(rho_bd[None], tables, (1, 1, 1), n)[0]
            with self.subTest(itype=t):
                self.assertGreater(np.abs(a - b).max(), 1e-3)


class TestOutputs(_TmpDir):
    """The archives of a short spin-orbital run carry the single spinful
    channel ``chiq`` and its index convention; ``physics`` carries the spin
    vector."""

    def test_archive_tags(self):
        from hwave.solver import spin_orbital
        s, gi = _build(self.d, _write_model(self.d, True, soc=0.3), IterationMax=3)
        with _quiet():
            s.solve(gi, self.d)
            s.save_results({"path_to_output": self.d, "chiq": "chiq.npz", "chi0q": "chi0q.npz",
                            "sigma": "sigma.npz", "green": "green.npz",
                            "energy": "energy.dat"}, gi)
        self.assertNotIn("chiq_s", gi)
        self.assertNotIn("chiq_c", gi)
        no = 2 * NP
        with np.load(os.path.join(self.d, "chiq.npz")) as z:
            self.assertEqual(str(z["chi_convention"]), "spin_orbital")
            self.assertEqual(str(z["index_convention"]), spin_orbital.INDEX_CONVENTION)
            self.assertNotIn("chiq_s", z.files)
            self.assertNotIn("chiq_c", z.files)
            self.assertEqual(z["chiq"].shape[-4:], (no, no, no, no))
        for name in ("sigma", "green"):
            with np.load(os.path.join(self.d, name + ".npz")) as z, self.subTest(archive=name):
                self.assertEqual(str(z["index_convention"]), spin_orbital.INDEX_CONVENTION)
                self.assertEqual(z[name].shape[-2:], (no, no))
        with np.load(os.path.join(self.d, "chi0q.npz")) as z:
            self.assertEqual(str(z["index_convention"]), spin_orbital.INDEX_CONVENTION)
            self.assertEqual(z["chi0q"].shape[-4:], (no, no, no, no))
        self.assertFalse(os.path.exists(os.path.join(self.d, "chiq_s.npz")))
        phys = gi["physics"]
        for k in ("NCond", "mu", "Sx", "Sy", "Sz"):
            self.assertIn(k, phys)
        self.assertAlmostEqual(phys["NCond"], 0.4 * 2 * NP * s.lattice.nvol, delta=1e-8)
        # the energy file carries the whole spin vector, as UHFk's does
        with open(os.path.join(self.d, "energy.dat")) as fh:
            lines = dict(l.split(" = ") for l in fh.read().splitlines())
        for k in ("NCond", "Sz", "Sx", "Sy", "ChemicalPotential"):
            self.assertIn(k, lines)
        self.assertAlmostEqual(float(lines["Sx"]), phys["Sx"], places=12)
        self.assertAlmostEqual(float(lines["Sy"]), phys["Sy"], places=12)

    def test_spin_vector_from_the_density(self):
        """``_spin_vector`` is ``1/2 Tr[sigma_alpha rho]`` per cell, summed."""
        s, _ = _build(self.d, _write_model(self.d, True))
        rng = np.random.default_rng(5)
        no = 2 * NP
        A = rng.normal(size=(no, no)) + 1j * rng.normal(size=(no, no))
        rho0 = 0.1 * (A + A.conj().T)
        rho_r = np.zeros((s.lattice.nvol, no, no), complex)
        rho_r[0] = rho0
        paulis = {"Sx": np.array([[0, 1], [1, 0]]), "Sy": np.array([[0, -1j], [1j, 0]]),
                  "Sz": np.array([[1, 0], [0, -1]])}
        got = s._spin_vector(rho_r)
        m = rho0.reshape(2, NP, 2, NP)
        for k, p in paulis.items():
            # rho[s a, s' a] = <c^dag_{s a} c_{s' a}>; <S> = 1/2 sum sigma_{s s'} <c^dag_s c_s'>
            ref = 0.5 * sum(p[s1, s2] * np.trace(m[s1, :, s2, :])
                            for s1 in range(2) for s2 in range(2))
            with self.subTest(component=k):
                self.assertAlmostEqual(got[k], float(ref.real) * s.lattice.nvol, places=12)


if __name__ == "__main__":
    unittest.main()
