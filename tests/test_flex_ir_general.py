"""IR-basis (sparse-ir) FLEX for the general (full-vertex / MYO) scheme.

Mirrors tests/test_flex_ir.py but for calc_scheme='general' on a norb=2
ON-SITE fixture. Equivalence is convergence of the uniform result toward the
IR result as Nmat grows (the uniform path carries the O(beta/Nmat) artifact).
Run from the repository root.
"""
import os

import numpy as np
import pytest

try:
    import sparse_ir  # noqa: F401
    _HAVE_SPARSE_IR = True
except ImportError:
    _HAVE_SPARSE_IR = False

pytestmark = pytest.mark.skipif(not _HAVE_SPARSE_IR,
                                reason="sparse-ir not installed")

# Reuse the on-site 2-orbital fixture writer from the general-scheme tests.
from tests.test_flex_general import _write_2d_2orb_onsite_fixture


def _make_general_solver(nmat, matsubara_basis="uniform", T=2.0,
                         iteration_max=60, mix=0.3):
    """A spin-free general FLEX solver on the on-site 2-orbital fixture."""
    import tempfile
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as solver_flex
    dirpath = os.path.join(tempfile.gettempdir(), "hwave_flex_2d_2orb_onsite")
    _write_2d_2orb_onsite_fixture(dirpath)
    info_input = {'path_to_input': dirpath, 'interaction': {
        'path_to_input': dirpath, 'Geometry': 'geom.dat',
        'Transfer': 'transfer.dat', 'CoulombIntra': 'coulombintra.dat',
        'CoulombInter': 'coulombinter.dat'}}
    ham = read_input_k.QLMSkInput(info_input).get_param("ham")
    param = {'T': T, 'mu': 0.0, 'CellShape': [4, 4, 1], 'SubShape': [1, 1, 1],
             'Nmat': nmat, 'IterationMax': iteration_max, 'Mix': mix, 'EPS': 8}
    # Construct with matsubara_basis='uniform' unconditionally: flex.py still
    # has a construction-time guard (`raise ValueError` when
    # use_ir and self._flex_general), pinned by the existing
    # tests/test_flex_ir.py::test_ir_rejects_general_scheme; removing it is
    # out of scope for this task (it belongs to the SCF-dispatch task that
    # follows this one) and this test must not touch it. Instead, flip the
    # two purely-descriptive post-construction attributes that
    # `matsubara_basis='ir'` would otherwise have set. This is safe here:
    # `_init_flex_param` (the only __init__ code that reads `self.use_ir`)
    # only consults it at the guard itself and at the already-satisfied
    # `write_densified` check -- no other construction-time branch depends
    # on it, and this test never calls solve() (which is where any
    # remaining use_ir-dependent SCF dispatch lives).
    if matsubara_basis != "uniform":
        param['matsubara_basis'] = "uniform"
    info_mode = {'mode': 'FLEX', 'param': param, 'calc_scheme': 'general'}
    gi = read_input_k.QLMSkInput(info_input).get_param("green")
    solver = solver_flex.FLEX(ham, {}, info_mode)
    if matsubara_basis != "uniform":
        solver.matsubara_basis = matsubara_basis
        solver.use_ir = (matsubara_basis == "ir")
    return solver, gi


def test_grev_index_map_tau_flip_only_spatial_roll_flip():
    """§6.10: the reverse+roll used by the general IR bubble maps tau by
    flip-only (j -> nt-1-j) and each spatial axis by roll(-1)+flip
    ((n-x) mod n); block/orbital axes untouched."""
    xp = np
    g, nt, nx, ny, nz, no2 = 1, 5, 4, 4, 1, 4
    rng = np.random.default_rng(0)
    a = (rng.standard_normal((g, nt, nx, ny, nz, no2))
         + 1j * rng.standard_normal((g, nt, nx, ny, nz, no2)))
    ref = -xp.flip(xp.roll(a, -1, axis=(2, 3, 4)), axis=(1, 2, 3, 4))
    # explicit elementwise map
    out = np.empty_like(a)
    for j in range(nt):
        jm = nt - 1 - j
        for x in range(nx):
            xm = (nx - x) % nx
            for y in range(ny):
                ym = (ny - y) % ny
                for z in range(nz):
                    zm = (nz - z) % nz
                    out[:, j, x, y, z, :] = -a[:, jm, xm, ym, zm, :]
    np.testing.assert_allclose(ref, out, atol=1e-12)


def test_chi0_general_gate_ir_matches_uniform_large_nmat():
    """GATE: the general IR-native chi0 (rank-6) must agree with the uniform
    general chi0 at large Nmat, compressed onto the same bosonic nodes."""
    from hwave.solver import eliashberg_dynamic as ed
    T = 2.0
    beta = 1.0 / T
    diffs = []
    for nmat in (256, 1024):
        s, gi = _make_general_solver(nmat, "ir", T=T)
        s._calc_epsilon_k(gi)
        s._ir_setup(beta)
        axB = s._ir_axB
        g_nodes, _ = s._calc_green_ir(beta, 0.0)
        chi0_ir = s._calc_chi0q_general_ir(g_nodes, beta)[0]      # strip block

        su, giu = _make_general_solver(nmat, "uniform", T=T)
        su._calc_epsilon_k(giu)
        gu, gtail = su._calc_green(beta, 0.0)
        chi0_u = su._calc_chi0q(gu, gtail, beta)[0]               # rank-6 (l,r,a,c,b,d)
        chi0_u_nodes = ed._ir_compress(
            np.moveaxis(chi0_u, 0, -1), axB, nmat, "chi0_u_gen",
            drop_constant=True)
        chi0_u_nodes = np.moveaxis(chi0_u_nodes, -1, 0)

        scale = np.abs(chi0_u_nodes).max()
        diffs.append(np.abs(chi0_ir - chi0_u_nodes).max() / scale)
    assert diffs[-1] < 1e-2, "general chi0 gate failed: {}".format(diffs)
    assert diffs[-1] < 0.6 * diffs[0], \
        "no Nmat convergence toward IR chi0: {}".format(diffs)


def test_chi0_general_ir_orbital_covariance():
    """A norb-permutation P of the bare G permutes the general IR chi0 output
    orbital axes correspondingly (guards swapped MYO pair indices)."""
    T = 2.0
    beta = 1.0 / T
    s, gi = _make_general_solver(256, "ir", T=T)
    s._calc_epsilon_k(gi)
    s._ir_setup(beta)
    g_nodes, _ = s._calc_green_ir(beta, 0.0)          # (g, nF, r, no, no)
    chi0 = s._calc_chi0q_general_ir(g_nodes, beta)     # (g, nB, r, a,c,b,d)
    P = np.array([1, 0])                               # swap the two orbitals
    g_perm = g_nodes[:, :, :, P][:, :, :, :, P]
    chi0_perm = s._calc_chi0q_general_ir(g_perm, beta)
    chi0_from_chi0 = chi0[:, :, :, P][:, :, :, :, P][:, :, :, :, :, P] \
        [:, :, :, :, :, :, P]
    np.testing.assert_allclose(chi0_perm, chi0_from_chi0, atol=1e-10)
