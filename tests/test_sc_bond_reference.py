"""The static ``hwave_sc`` bond operator against the previous implementation,
for nonzero CoulombInter (issue #239).

Issue #239 rebuilt the static bond path (``[eliashberg] bond_channels =
true``, ``frequency = "static"``) on the bond modules of FLEX and of the
dynamic bond pairing. For CoulombInter-only input at ``norb = 1`` the new
operator must equal the previous, separately implemented chain
(``resolve_interactions`` -> ``_build_bond_m0_blocks`` ->
``bare_bond_vertices`` -> ``dress_bond`` -> ``make_bond_kernel_parts``) to
round-off. The previous chain is run from the reference checkout of
``tests.test_eliashberg_bond_entry.BOND_REFERENCE_COMMIT`` in a subprocess
(the same provisioning as the bond-archive comparison there: the checkout is
pointed at by ``HWAVE_DEVELOP_CHECKOUT_BOND``; without it the test skips, and
where ``HWAVE_REQUIRE_DEVELOP_COMPARISON_BOND`` is set it fails instead).

Both revisions are fed the SAME Green function (an external carrier, so the
chemical-potential search and the tail handling of either revision play no
part) and the same interaction tables; the full operator and its bare and
fluctuation parts are applied to every basis vector, and the leading
eigenvalues of the resulting dense matrices are compared too. Cases: the
4x4 single-band fixture of ``tests/rpa/input`` (``U = 4``) with a
nearest-neighbour ``V = 0.5`` and ``1.0`` in both pairing channels, and the
``V = 1.2`` Onari point on its committed L = 16 FLEX Green function.
The interaction is symmetric under ``R -> -R``, so the bond-placement change
of issue #232 (after the reference revision) does not enter.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np

import hwave.sc as sc

NN = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))

#: The driver the reference revision runs: its own ``_build_bond_operator``
#: (previous signature) on the inputs in ``argv[1]``; the dense operator
#: parts go to ``argv[2]``.
_REFERENCE_DRIVER = r'''
import sys
import numpy as np
import hwave.sc as sc
from hwave.solver import bond_channels as bc

d = np.load(sys.argv[1])
green = d["green"]
beta = float(d["beta"])
U, V = float(d["U"]), float(d["V"])
eta = str(d["eta"])
Nx, Ny, Nz = green.shape[2:5]
kx, ky, kz = (np.linspace(0, 2 * np.pi, n, endpoint=False) for n in (Nx, Ny, Nz))
NN = ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
it = {"CoulombIntra": {((0, 0, 0), (0, 0)): U},
      "CoulombInter": {(r, (0, 0)): V for r in NN}}
inter_k = sc._build_interaction_k(kx, ky, kz, it, 1)
bond_set = bc.resolve_interactions(it["CoulombInter"], np.eye(3), 1)
carrier = sc._BondGreen(full_sc=green, deflated_kw=sc._green_sc_to_canonical(green),
                        green0_tail=None, coeff_tail_requested=0.0, coeff_tail_applied=0.0,
                        source="external_npz")
(A, n), _, attr = sc._build_bond_operator(
    bond_set, carrier, it, inter_k, {"norb": 1, "rvec": np.eye(3)}, 1, kx, ky, kz, beta, eta,
    bond_max_shells=None, bond_memory_cap_gb=8.0, g2_tail=bool(d["g2_tail"]))
eye = np.eye(n)
np.savez(sys.argv[2], **{k: np.column_stack([attr[k].matvec(eye[:, j]) for j in range(n)])
                         for k in ("full", "pp", "fl")})
'''


def _dense_parts(green, beta, U, V, eta, g2_tail):
    """The dense ``full``/``pp``/``fl`` matrices of this tree's operator."""
    Nx, Ny, Nz = green.shape[2:5]
    kx, ky, kz = (np.linspace(0, 2 * np.pi, n, endpoint=False) for n in (Nx, Ny, Nz))
    it = {"CoulombIntra": {((0, 0, 0), (0, 0)): U},
          "CoulombInter": {(r, (0, 0)): V for r in NN}}
    topo, _, split = sc._resolve_bond_topology(it, {"rvec": np.eye(3)}, 1, None)
    carrier = sc._BondGreen(full_sc=green, deflated_kw=sc._green_sc_to_canonical(green),
                            green0_tail=None, coeff_tail_requested=0.0, coeff_tail_applied=0.0,
                            source="external_npz")
    (A, n), _, attr = sc._build_bond_operator(
        topo, split, carrier, 1, kx, ky, kz, beta, eta,
        bond_max_shells=None, bond_memory_cap_gb=8.0, g2_tail=g2_tail)
    eye = np.eye(n)
    return {k: np.column_stack([attr[k].matvec(eye[:, j]) for j in range(n)])
            for k in ("full", "pp", "fl")}


def _small_green():
    """Bare Green function of the 4x4 ``tests/rpa/input`` model (T = 1,
    Nmat = 128, half filling of the band), sc.py layout."""
    geom, hr, _ = sc._read_interaction_files(
        {"file": {"input": {"interaction": {
            "path_to_input": "tests/rpa/input", "Geometry": "geom.dat",
            "Transfer": "transfer.dat"}}}})
    k = np.linspace(0, 2 * np.pi, 4, endpoint=False)
    ev, evec = sc._calc_eigenvalues(sc._build_hamiltonian_k(k, k, np.zeros(1), hr, 1))
    beta = 1.0
    mu = sc._determine_mu(ev, beta, 0.5, 1)
    return sc._build_bond_green(ev, evec, mu, beta, 128, 0.0, None).full_sc, beta


def _onari_green():
    """The committed L = 16, V = 1.2 FLEX Green function of the milestone."""
    raw = np.load(os.path.join("tests", "sc", "onari_bond", "green_L16_V1.20.npz"))["green"]
    nmat = raw.shape[1]
    return raw[0].reshape(nmat, 16, 16, 1, 1, 1).transpose(4, 5, 1, 2, 3, 0).copy(), 50.0


class TestStaticBondOperatorMatchesTheReferenceRevision(unittest.TestCase):

    #: relative bound on every matrix element (scaled by the largest one) and
    #: on the leading eigenvalues: round-off of two summation orders
    RTOL = 1e-12

    @classmethod
    def setUpClass(cls):
        from tests.test_eliashberg_bond_entry import (BOND_REFERENCE_COMMIT,
                                                      BOND_REFERENCE_PATH_ENV,
                                                      BOND_REFERENCE_REQUIRE_ENV)
        from tests.test_flex_second_order_compat import develop_checkout
        cls.dev, why = develop_checkout(expected_commit=BOND_REFERENCE_COMMIT,
                                        path_env=BOND_REFERENCE_PATH_ENV,
                                        require_env=BOND_REFERENCE_REQUIRE_ENV)
        if cls.dev is None:
            raise unittest.SkipTest(why)
        cls.work = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.work, ignore_errors=True)

    def _reference(self, green, beta, U, V, eta, g2_tail):
        from tests.test_flex_second_order_compat import reference_subprocess_env
        tag = "{}_{}_{}_{}".format(green.shape[2], V, eta, int(g2_tail))
        inp = os.path.join(self.work, "in_{}.npz".format(tag))
        out = os.path.join(self.work, "out_{}.npz".format(tag))
        np.savez(inp, green=green, beta=beta, U=U, V=V, eta=eta, g2_tail=g2_tail)
        subprocess.run([sys.executable, "-B", "-c", _REFERENCE_DRIVER, inp, out],
                       env=reference_subprocess_env(self.dev), check=True,
                       capture_output=True, cwd=self.dev)
        with np.load(out) as d:
            return {k: d[k] for k in d.files}

    def _compare(self, green, beta, U, V, eta, g2_tail):
        ref = self._reference(green, beta, U, V, eta, g2_tail)
        new = _dense_parts(green, beta, U, V, eta, g2_tail)
        for part in ("full", "pp", "fl"):
            scale = np.abs(ref[part]).max()
            self.assertGreater(scale, 0.0, part)
            err = np.abs(new[part] - ref[part]).max() / scale
            self.assertLess(err, self.RTOL, (V, eta, part, err))
        ev_ref = np.sort_complex(np.linalg.eigvals(ref["full"]))[::-1][:6]
        ev_new = np.sort_complex(np.linalg.eigvals(new["full"]))[::-1][:6]
        np.testing.assert_allclose(ev_new, ev_ref, rtol=0,
                                   atol=self.RTOL * np.abs(ev_ref).max() * 10,
                                   err_msg="{} {}".format(V, eta))

    def test_single_band_4x4(self):
        green, beta = _small_green()
        for V in (0.5, 1.0):
            for eta in ("singlet", "triplet"):
                with self.subTest(V=V, eta=eta):
                    self._compare(green, beta, 4.0, V, eta, True)

    def test_onari_l16_point(self):
        green, beta = _onari_green()
        for eta in ("triplet", "singlet"):
            with self.subTest(eta=eta):
                self._compare(green, beta, 4.0, 1.2, eta, False)


if __name__ == "__main__":
    unittest.main()
