"""Multi-orbital validation of the bond-resolved pairing by supercell
equivalence (issue #239).

A single-band model written on a two-site unit cell (``norb = 2``, half the
cells along ``x``) is the same Hamiltonian, so the multi-orbital bond pairing
must reproduce the single-band one, which is validated on its own (the
first-order oracle of ``tests/test_sc_bond_pairing_oracle.py``, the exact
diagonalization of ``tests/test_bond_pairing_ed.py``). The two-site cell
exercises every multi-orbital piece of the bond path: inter-orbital bonds in
both directions, the in-cell (``R = 0``) inter-orbital rows that a
single-band bond becomes, orbital-resolved channel blocks and the pair
crossing of the bare term, for every bond-resolved type (CoulombInter, Hund,
Ising, Exchange, PairHop).

* :class:`TestStaticKernel` -- the static bond kernel (full, fluctuation and
  bare part) of the two models agree to machine precision on gap vectors
  mapped between the two bases.
* :class:`TestDynamicKernel` -- so does the frequency-resolved dynamic
  kernel (second-order vertex of the dynamic bond bubble plus the bare term).
* :class:`TestStaticEntry` -- ``hwave_sc`` through its files: the single-band
  eigenvalues are eigenvalues of the supercell (whose spectrum also holds the
  pairs of momentum ``pi`` along ``x`` of the single-band model, folded to
  zero).
* :class:`TestDynamicEntries` -- the same for the dynamic bond pairing, after
  a FLEX bond-gate run of each model: the in-process pairing and the
  post-processing of the supercell archive.
"""
import itertools
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.sc as sc
from hwave.solver import bond_channels as bc
from tests.heavy_tests import heavy

#: single-band couplings on the nearest-neighbour x bonds of a square lattice
OFFSITE = {"CoulombInter": 0.4, "Hund": 0.15, "Ising": 0.05, "Exchange": 0.05,
           "PairHop": 0.1}
U = 3.0
#: single-band hopping rows (Rx, Ry): t
HOP = {(1, 0): -1.0, (-1, 0): -1.0, (0, 1): -1.0, (0, -1): -1.0,
       (1, 1): 0.2, (-1, -1): 0.2, (1, -1): 0.2, (-1, 1): 0.2}
#: the couplings sit on the x bonds only, the axis the supercell folds, so
#: the leading gaps live on the in-cell pairs a single-band x bond becomes
NN = ((1, 0), (-1, 0))


def _fold(rows):
    """Single-band rows ``{(Rx, Ry): v}`` -> two-site-cell rows
    ``{(RX, Ry, a, b): v}`` (site ``x = 2 X + a``)."""
    out = {}
    for (rx, ry), v in rows.items():
        for a in range(2):
            xp = a + rx
            b = xp % 2
            key = ((xp - b) // 2, ry, a, b)
            out[key] = out.get(key, 0.0) + v
    return out


def _model(norb, types=tuple(OFFSITE)):
    """``(transfer, interactions)`` as ``{(Rx, Ry, a, b): v}`` tables."""
    single = {"Transfer": dict(HOP)}
    single.update({t: {R: OFFSITE[t] for R in NN} for t in types})
    if norb == 1:
        tables = {k: {(R[0], R[1], 0, 0): v for R, v in rows.items()}
                  for k, rows in single.items()}
        tables["CoulombIntra"] = {(0, 0, 0, 0): U}
    else:
        tables = {k: _fold(rows) for k, rows in single.items()}
        tables["CoulombIntra"] = {(0, 0, a, a): U for a in range(2)}
    return tables


def _declarations(tables):
    return {t: {((R0, R1, 0), (a, b)): complex(v) for (R0, R1, a, b), v in rows.items()}
            for t, rows in tables.items() if t != "Transfer"}


# ---------------------------------------------------------------------------
# the static kernel
# ---------------------------------------------------------------------------

BETA, NMAT, MU = 2.0, 64, 0.2
NX1, NY = 8, 4


def _green(norb, nx, tables):
    kx = 2 * np.pi * np.arange(nx) / nx
    ky = 2 * np.pi * np.arange(NY) / NY
    h = np.zeros((nx, NY, norb, norb), complex)
    for (r0, r1, a, b), t in tables["Transfer"].items():
        h[:, :, a, b] += t * np.exp(1j * (kx[:, None] * r0 + ky[None, :] * r1))
    wn = (2 * np.arange(NMAT) + 1 - NMAT) * np.pi / BETA
    g = np.empty((norb, norb, nx, NY, 1, NMAT), complex)
    eye = np.eye(norb)
    for n, w in enumerate(wn):
        g[:, :, :, :, 0, n] = np.linalg.inv((1j * w + MU) * eye - h).transpose(2, 3, 0, 1)
    return g


def _kernel_parts(norb, eta, types):
    nx = NX1 // norb
    tables = _model(norb, types)
    g = _green(norb, nx, tables)
    shape = (nx, NY, 1)
    topo, view, split = sc._resolve_bond_topology(_declarations(tables), {"rvec": np.eye(3)},
                                                  norb, None)
    chi_bar = bc.bond_bubble(g, view, BETA).reshape(nx * NY, view.n_channels * norb ** 2, -1)
    S, C = sc._build_bond_vertices(topo, split, norb, shape)
    dressed = ((ch, bc.dress_batch(chi_bar[None], W, ch, l0=0, nmat=0, spatial_shape=shape)[0][0])
               for ch, W in (("spin", S), ("charge", C)))
    A, A_fl, A_pp, _ = sc._bond_kernel_parts(S, C, dressed, g, view, norb, shape, BETA, eta,
                                             topo=topo)
    return A, A_fl, A_pp


def _supercell_map():
    """``M``: a single-band gap ``Delta(k)`` (``NX1 x NY``) -> the two-site-cell
    gap ``Delta_ab(K)`` (``2 x 2 x NX1/2 x NY``): ``Delta_ab(R) = Delta(r)`` with
    ``r = 2 R + a - b`` along ``x``."""
    n2 = NX1 // 2
    k1 = 2 * np.pi * np.arange(NX1) / NX1
    K2 = 2 * np.pi * np.arange(n2) / n2
    M = np.zeros((2, 2, n2, NX1), complex)
    for a, b in itertools.product(range(2), range(2)):
        for R in range(n2):
            row = np.exp(1j * k1 * (2 * R + a - b)) / NX1
            M[a, b] += np.exp(-1j * K2[:, None] * R) * row[None, :]
    # the y axis is shared: Delta(kx, ky) -> Delta_ab(Kx, ky)
    return np.einsum("abKk,yz->abKykz", M, np.eye(NY)).reshape(4 * n2 * NY, NX1 * NY)


class TestStaticKernel(unittest.TestCase):

    def test_supercell_reproduces_the_single_band_kernel(self):
        M = _supercell_map()
        mk = (-np.arange(NX1))[:, None] % NX1, (-np.arange(NY))[None, :] % NY
        rng = np.random.default_rng(1)
        for types in [tuple(OFFSITE)] + [(t,) for t in OFFSITE]:
            for eta in ("singlet", "triplet"):
                one = _kernel_parts(1, eta, types)
                two = _kernel_parts(2, eta, types)
                for name, A1, A2 in zip(("full", "fl", "pp"), one, two):
                    with self.subTest(types=types, eta=eta, part=name):
                        v = (rng.standard_normal((NX1, NY))
                             + 1j * rng.standard_normal((NX1, NY)))
                        v = 0.5 * (v + (1 if eta == "singlet" else -1) * v[mk])
                        v = v.ravel()
                        want = M @ A1.matvec(v)
                        got = A2.matvec(M @ v)
                        self.assertLess(np.linalg.norm(got - want),
                                        1e-12 * np.linalg.norm(want) + 1e-14)


NMAT_DYN = 16


def _dynamic_kernel(norb, eta):
    """The dynamic bond pairing kernel (``BondPairKernel`` on the full
    bosonic grid) of the bare Green function, with the second-order vertex
    ``Gamma = c_s S chibar S + c_c C chibar C`` of the dynamic bond bubble and
    the bare term: the frequency-resolved legs of the kernel at ``norb``."""
    from hwave.solver import bubble
    from hwave.solver.eliashberg_bond import (COEFF, ArrayBlockSource, BondPairKernel,
                                              PairVertexUniform, instantaneous_vertex)
    from hwave.solver.eliashberg_dynamic import calc_g2_dynamic
    global NMAT
    nmat0, NMAT = NMAT, NMAT_DYN
    try:
        nx = NX1 // norb
        tables = _model(norb)
        g = _green(norb, nx, tables)
    finally:
        NMAT = nmat0
    shape = (nx, NY, 1)
    nvol, nd = nx * NY, norb * norb
    topo, view, split = sc._resolve_bond_topology(_declarations(tables), {"rvec": np.eye(3)},
                                                  norb, None)
    S, C = sc._build_bond_vertices(topo, split, norb, shape)
    S, C = (np.asarray(W).reshape(nvol, *np.shape(W)[-2:]) for W in (S, C))
    ND = view.n_channels * nd
    chib = np.zeros((NMAT_DYN, nvol, ND, ND), complex)
    for (m, mp), blk in bubble._iter_bond_dynamic(sc._green_sc_to_canonical(g), None, BETA,
                                                  view, spatial_shape=shape):
        chib[:, :, m * nd:(m + 1) * nd, mp * nd:(mp + 1) * nd] = blk.reshape(NMAT_DYN, nvol, nd, nd)
    cs, cc = COEFF[eta]
    gam = cs * S[None] @ chib @ S[None] + cc * C[None] @ chib @ C[None]
    vert = PairVertexUniform(source=ArrayBlockSource({"G": gam}, nd), slot="G",
                             B=view.n_channels, nd=nd, nmat=NMAT_DYN, nvol=nvol)
    return BondPairKernel(vert, calc_g2_dynamic(g, BETA), view, xp=np, spatial_shape=shape,
                          norb=norb, beta=BETA, nfreq=NMAT_DYN,
                          V_inst=instantaneous_vertex(S, C, nd, eta, shape, topo=topo),
                          residency="host", host_cap=10 ** 12, device_cap=10 ** 12)


class TestDynamicKernel(unittest.TestCase):
    """The frequency-resolved kernel: the supercell reproduces the
    single-band one on mapped gaps of every frequency, to machine precision."""

    def test_supercell_reproduces_the_single_band_kernel(self):
        M = _supercell_map()
        rng = np.random.default_rng(2)
        mk = ((-np.arange(NX1))[:, None, None] % NX1, (-np.arange(NY))[None, :, None] % NY,
              np.arange(NMAT_DYN)[None, None, :])
        for eta in ("singlet", "triplet"):
            with self.subTest(eta=eta):
                K1, K2 = _dynamic_kernel(1, eta), _dynamic_kernel(2, eta)
                v = (rng.standard_normal((NX1, NY, NMAT_DYN))
                     + 1j * rng.standard_normal((NX1, NY, NMAT_DYN)))
                v = 0.5 * (v + (1 if eta == "singlet" else -1) * v[mk])
                # (kx, ky, n) -> the single-band layout (1, 1, kx, ky, 1, n)
                want = K1.matvec(v.ravel()).reshape(NX1 * NY, NMAT_DYN)
                want = (M @ want).ravel()
                got = K2.matvec((M @ v.reshape(NX1 * NY, NMAT_DYN)).ravel())
                self.assertLess(np.linalg.norm(got - want), 1e-12 * np.linalg.norm(want))


# ---------------------------------------------------------------------------
# the entries, through their input files
# ---------------------------------------------------------------------------

def _write_tables(d, norb, tables):
    """Geometry and one Wannier90-like file per table; returns the
    ``[file.input.interaction]`` table."""
    with open(os.path.join(d, "geom.dat"), "w") as fw:
        fw.write("{} 0 0\n0 1 0\n0 0 1\n{}\n".format(norb, norb))
        for a in range(norb):
            fw.write("{} 0 0\n".format(a / norb))
    inter = {"path_to_input": d, "Geometry": "geom.dat"}
    for t, rows in tables.items():
        name = t.lower() + ".dat"
        with open(os.path.join(d, name), "w") as fw:
            fw.write("{} for the supercell test\n{}\n{}\n".format(t, norb, len(rows)))
            # the degeneracy list, 15 per line (Wannier90)
            for i in range(0, len(rows), 15):
                fw.write(" ".join(["1"] * min(15, len(rows) - i)) + "\n")
            for (r0, r1, a, b), v in sorted(rows.items()):
                fw.write("{:4d} {:4d} {:4d} {:4d} {:4d} {:.12f} 0.0\n".format(
                    r0, r1, 0, a + 1, b + 1, v))
        inter[t] = name
    return inter


def _eigenvalues(outdir):
    """Re(eigenvalue) of the eigenvalue-analysis rows of the channel's parity."""
    with open(os.path.join(outdir, "eigenvalue.dat")) as f:
        rows = [l.split() for l in f if l.strip() and not l.lstrip().startswith("#")]
    return np.array([float(r[1]) for r in rows if len(r) == 5 and r[4] == "1"])


class _Case(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="hwave_bond_supercell_")
        self.addCleanup(shutil.rmtree, self.tmp, ignore_errors=True)

    def _dir(self, *parts):
        d = os.path.join(self.tmp, *parts)
        os.makedirs(d, exist_ok=True)
        return d

    def assertContained(self, single, cell, rel=1e-8):
        """Every single-band eigenvalue is an eigenvalue of the supercell."""
        self.assertGreater(len(single), 0)
        for lam in np.atleast_1d(single):
            dev = np.min(np.abs(np.asarray(cell) - lam))
            self.assertLess(dev, rel * max(abs(lam), 1.0),
                            "single-band eigenvalue {} not in the supercell spectrum {}"
                            .format(lam, cell))


class TestStaticEntry(_Case):

    def _run(self, norb, eta, n_eig):
        d = self._dir("static", eta, str(norb))
        inter = _write_tables(d, norb, _model(norb))
        inp = {"mode": {"param": {"T": 1.0 / BETA, "CellShape": [NX1 // norb, NY, 1],
                                  "SubShape": [1, 1, 1], "Nmat": NMAT, "filling": 0.45}},
               "file": {"input": {"interaction": inter}, "output": {"path_to_output": d}},
               "eliashberg": {"chi0q_mode": "calc", "bond_channels": True,
                              "pairing_type": eta, "solver_mode": "eigenvalue",
                              "num_eigenvalues": n_eig}}
        sc.calc_eliashberg(inp)
        return _eigenvalues(d)

    def test_single_band_eigenvalues_are_supercell_eigenvalues(self):
        for eta in ("singlet", "triplet"):
            with self.subTest(eta=eta):
                self.assertContained(self._run(1, eta, 2), self._run(2, eta, 8))


def _flex_input(d, norb, inter, pairing):
    par = {"T": 1.0 / BETA, "filling": 0.45, "CellShape": [NX1 // norb, NY, 1],
           "SubShape": [1, 1, 1], "Nmat": NMAT, "IterationMax": 300, "Mix": 0.3, "EPS": 10,
           "flex_hartree_fock": True, "longitudinal_bond_channels": True,
           "longitudinal_bond_output_full": True, "mixing_scheme": "linear",
           "longitudinal_bond_pairing": pairing}
    return {"mode": {"mode": "FLEX", "calc_scheme": "general", "param": par},
            "file": {"input": {"path_to_input": "", "interaction": inter},
                     "output": {"path_to_output": d, "chiq_s": "chiq_s", "chiq_c": "chiq_c",
                                "green": "green", "longitudinal_bond": "longitudinal_bond"}},
            "eliashberg": {"solver_mode": "eigenvalue", "num_eigenvalues": 8}}


def _all_eigenvalues(path):
    """The complex eigenvalue of every eigenvalue-analysis row."""
    with open(path) as f:
        rows = [l.split() for l in f if l.strip() and not l.lstrip().startswith("#")]
    return np.array([float(r[1]) + 1j * float(r[2]) for r in rows if len(r) == 5])


class TestDynamicEntries(_Case):
    """A FLEX bond-gate run of each model (they converge to the same fixed
    point), then the dynamic bond pairing. The supercell kernel conserves only
    the combined (momentum x frequency) parity, so its spectrum also holds the
    odd-frequency gaps and the folded momentum-pi pairs: it is asked for
    enough eigenvalues to reach the single-band ones."""

    N_SINGLE, N_CELL = 2, 80

    def _flex(self, norb):
        import hwave.qlms as qlms
        d = self._dir("flex", str(norb))
        inter = _write_tables(d, norb, _model(norb))
        inp = _flex_input(d, norb, inter, "both")
        inp["eliashberg"]["num_eigenvalues"] = self.N_SINGLE if norb == 1 else self.N_CELL
        qlms.run(input_dict=inp)
        return d, inter

    @heavy
    def test_in_process_and_post_processing(self):
        d1, _ = self._flex(1)
        d2, inter2 = self._flex(2)
        for eta in ("singlet", "triplet"):
            with self.subTest(eta=eta, entry="in-process"):
                for d in (d1, d2):
                    with np.load(os.path.join(d, "eliashberg_bond_{}.npz".format(eta))) as z:
                        self.assertTrue(bool(z["scf_converged"]), d)
                one = _all_eigenvalues(os.path.join(d1, "eigenvalue_bond_{}.dat".format(eta)))
                two = _all_eigenvalues(os.path.join(d2, "eigenvalue_bond_{}.dat".format(eta)))
                self.assertEqual(len(one), self.N_SINGLE)
                self.assertContained(one, two, rel=1e-7)
            with self.subTest(eta=eta, entry="post-processing"):
                out = self._dir("post", eta)
                inp = {"mode": {"param": {"T": 1.0 / BETA, "CellShape": [NX1 // 2, NY, 1],
                                          "SubShape": [1, 1, 1], "Nmat": NMAT,
                                          "filling": 0.45}},
                       "file": {"input": {"path_to_flex_output": d2, "interaction": inter2},
                                "output": {"path_to_output": out}},
                       "eliashberg": {"chi0q_mode": "flex", "frequency": "dynamic",
                                      "bond_channels": True, "pairing_type": eta,
                                      "solver_mode": "eigenvalue",
                                      "num_eigenvalues": self.N_CELL}}
                sc.calc_eliashberg(inp)
                post = _all_eigenvalues(os.path.join(out, "eigenvalue.dat"))
                inproc = _all_eigenvalues(os.path.join(d2, "eigenvalue_bond_{}.dat".format(eta)))
                # both print eight significant digits
                self.assertLess(abs(post[0] - inproc[0]), 1e-7 * abs(inproc[0]))


if __name__ == "__main__":
    unittest.main()
