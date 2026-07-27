"""P0-1 gate: dense uniform Hermiticity check for the dynamic bond kernel.

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-design.md,
section 2 (the P0-1 amendment gate). This test decides whether the
Hermitian-solver design of the static bond kernel transfers to the dynamic
(frequency-resolved) bond kernel BEFORE any production dynamic-kernel code
(Task 6+) is written. Uses ONLY the independent oracle
(tests/oracle_bond_dynamic.py, Task 2) plus the existing static
bond_channels machinery for vertex construction -- no production dynamic
kernel exists yet.

Per the brief: if this gate fails, STOP. Do not tune, do not weaken the
tolerance, do not modify the oracle. Report the measured residuals and any
structure (instantaneous vs fluctuation) and let the spec-amendment gate
decide next steps.

SKIPPED BY DEFAULT (``@pytest.mark.slow``): densifying the oracle matvec
over all 128 basis vectors, per parity, takes ~3-6 minutes total on this
grid -- the same skip-guard pattern as
``tests/test_bond_onari_milestone.py``'s slow grid-convergence test. Run
with ``HWAVE_RUN_SLOW_FIXTURES=1`` or ``pytest -m slow``.
"""
import numpy as np
import pytest

from hwave.solver import bond_channels as bc
from tests.oracle_bond_dynamic import oracle_bond_bubble, oracle_bond_matvec
from tests.test_bond_onari_milestone import _require_slow_fixtures
from tests.test_oracle_bond_dynamic import _symmetric_green

BETA = 5.0


def _dense_operator(matvec, shape):
    size = int(np.prod(shape))
    M = np.zeros((size, size), dtype=complex)
    for i in range(size):
        e = np.zeros(size, dtype=complex); e[i] = 1.0
        M[:, i] = matvec(e.reshape(shape)).ravel()
    return M


def _build_vertices(bset, U, shape):
    """S0_q/C0_q for a single-band (norb=1) on-site U + NN-bond V model,
    built the same way ``sc.py._build_bond_m0_blocks`` builds them for the
    bond machinery (CoulombIntra U on the m=0 diagonal, Hartree 2*V(q) added
    to the m=0 charge block from the SAME resolved bond set that feeds
    ``bare_bond_vertices``), then handed to
    ``bond_channels.bare_bond_vertices`` -- no new S/C convention invented
    here, only the existing construction reused for norb=1 with no declared
    on-site (Delta r = 0) V (so the Fock (ab,ab) correction is zero, matching
    tests/test_bond_channels.py's ``test_pp_vertex_singleband_fixture`` /
    ``test_bond_kernel_triplet_lambda_minus_one`` U=4, V=1 single-band
    fixtures, generalized here to the full q-grid).
    """
    Nx, Ny, Nz = shape
    norb = 1
    nd = 1
    kx = 2.0 * np.pi * np.arange(Nx) / Nx
    ky = 2.0 * np.pi * np.arange(Ny) / Ny
    kz = 2.0 * np.pi * np.arange(Nz) / Nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")

    # V(q) = sum_{m>=1} V_bond[m] e^{-i q.dr_m}; same e^{-iq.R} convention as
    # sc.py._build_bond_m0_blocks / _build_interaction_k._to_k. No on-site
    # (R=0) V is declared in this model, so v_onsite is exactly zero and the
    # m=0 Fock element gets no correction.
    V_q = np.zeros((Nx, Ny, Nz), dtype=complex)
    for m in range(1, bset.n_channels):
        Rx, Ry, Rz = bset.delta_r[m]
        phase = np.exp(-1j * (KX * Rx + KY * Ry + KZ * Rz))
        V_q += bset.v_bond[m][0, 0] * phase

    S0 = np.full((Nx, Ny, Nz, nd, nd), U, dtype=complex)
    C0 = np.full((Nx, Ny, Nz, nd, nd), U, dtype=complex)
    C0[:, :, :, 0, 0] += 2.0 * V_q

    return bc.bare_bond_vertices(bset, S0, C0, norb)


@pytest.mark.slow
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_p0_1_uniform_hermiticity_gate(pairing, request):
    _require_slow_fixtures(request)
    Nx = Ny = 4; Nz = 1; nmat = 8
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, BETA)
    # --- v1 symmetry-class preconditions (spec section 2) ---
    gmax = np.abs(green).max()
    g_flip_k = np.roll(green[:, :, ::-1, ::-1, ::-1, :], (1, 1, 1), (2, 3, 4))
    assert np.abs(g_flip_k - green).max() <= 1e-10 * gmax
    assert np.abs(green[..., ::-1] - np.conj(green)).max() <= 1e-10 * gmax
    w_complex = (green * np.roll(green[:, :, ::-1, ::-1, ::-1, ::-1],
                                 (1, 1, 1), (2, 3, 4)))[0, 0]
    # w must be real (up to round-off) BEFORE the .real truncation below
    # silently discards any imaginary part -- tol_sym = 1e-10 relative.
    assert (np.abs(w_complex.imag).max()
            <= 1e-10 * np.abs(w_complex).max())
    w = w_complex.real
    assert w.min() >= 1e-12 * w.max()
    # --- bond topology + vertices (reuse the static machinery) ---
    # dict format is {(irvec, orbvec): value} per resolve_interactions'
    # docstring / tests/test_bond_channels.py's _nn_square helper -- the
    # brief's literal {irvec: matrix} sketch does not match the actual API,
    # so this uses the real format with the same NN V=1 topology.
    inter = {((1, 0, 0), (0, 0)): 1.0, ((-1, 0, 0), (0, 0)): 1.0,
             ((0, 1, 0), (0, 0)): 1.0, ((0, -1, 0), (0, 0)): 1.0}
    bset = bc.resolve_interactions(inter, np.eye(3), 1)
    B = bset.n_channels; ND = B  # nd = 1
    assert B == 5
    # S0_q/C0_q inputs: follow the construction used in
    # tests/test_bond_channels.py for the same U=4 + NN V model (read that
    # file; the helper builds S/C with U on the m=0 block, Hartree 2V(q) in
    # C, Fock +/-V on the bond diagonal).
    S_bond, C_bond, Vpp_s, Vpp_t = _build_vertices(bset, U=4.0,
                                                   shape=(Nx, Ny, Nz))
    # --- frequency-resolved dressed F_w via plain linalg ---
    jlist = list(range(-(nmat // 2), nmat // 2))
    chib = oracle_bond_bubble(green, bset.delta_r, BETA, jlist)
    F_w = np.zeros((Nx, Ny, Nz, nmat, ND, ND), dtype=complex)
    for jt in jlist:
        j = jt + nmat // 2
        cb = chib[jt]
        for ix in np.ndindex(Nx, Ny, Nz):
            S = S_bond[ix]; C = C_bond[ix]; X = cb[ix]
            chi_s = np.linalg.solve(np.eye(ND) - X @ S, X)
            chi_c = np.linalg.solve(np.eye(ND) + X @ C, X)
            if pairing == "singlet":
                F_w[ix + (j,)] = 1.5 * S @ chi_s @ S - 0.5 * C @ chi_c @ C
            else:
                F_w[ix + (j,)] = -0.5 * S @ chi_s @ S - 0.5 * C @ chi_c @ C
    Vpp = Vpp_s if pairing == "singlet" else Vpp_t
    shape = (1, 1, Nx, Ny, Nz, nmat)
    K = _dense_operator(
        lambda p: oracle_bond_matvec(p, green, F_w, Vpp, bset.delta_r, BETA),
        shape)
    sqw = np.sqrt(w).ravel()  # gap-vector order = (k, n) raveled like phi
    Kt = (sqw[:, None] * K) / sqw[None, :]
    res = np.linalg.norm(Kt - Kt.conj().T) / np.linalg.norm(Kt)
    assert res <= 1e-10, f"P0-1 GATE FAILED: residual {res:.3e}"
