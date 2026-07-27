"""Task 6: the uniform-axis dynamic bond Eliashberg kernel.

Normative reference: docs/superpowers/specs/2026-07-27-dynamic-bond-channels-
design.md sections 3.3/3.5. The acceptance test is array-wise equality of
``bond_channels.make_bond_kernel_dynamic`` against the INDEPENDENT direct-sum
oracle ``tests/oracle_bond_dynamic.oracle_bond_matvec`` (Task 2), which may
not import any hwave transform/backend helper. The mutation battery at the
bottom demonstrates that this equality actually discriminates -- each
mutation corrupts a PRODUCTION input only and the comparison must then fail.
"""
import dataclasses

import numpy as np
import pytest
from scipy.sparse.linalg import LinearOperator

from hwave.solver import bond_channels as bc
from hwave.solver import eliashberg_dynamic as ed
from tests.oracle_bond_dynamic import oracle_bond_matvec, oracle_pair_bubble_X
from tests.test_bond_channels_dynamic import _nn_bond_set
from tests.test_bond_dynamic_hermiticity import _build_vertices
from tests.test_oracle_bond_dynamic import _symmetric_green

BETA = 5.0


def _vertices_full(bset, shape):
    """The brief's ``_vertices_full``: the existing U=4 / NN-V=1 single-band
    vertex construction from tests/test_bond_dynamic_hermiticity.py, on the
    full q-grid. Reused rather than re-derived (no new S/C convention here).
    """
    return _build_vertices(bset, U=4.0, shape=shape)


def _sandwich(V, chi_w):
    """(Nx,Ny,Nz,ND,ND) x (Nx,Ny,Nz,nmat,ND,ND) -> V chi V per (q, j)."""
    return np.einsum('xyzab,xyzjbc,xyzcd->xyzjad', V, chi_w, V)


def _oracle_F_w(pairing, S_bond, C_bond, chi_s_w, chi_c_w):
    if pairing == "singlet":
        return (1.5 * _sandwich(S_bond, chi_s_w)
                - 0.5 * _sandwich(C_bond, chi_c_w))
    return (-0.5 * _sandwich(S_bond, chi_s_w)
            - 0.5 * _sandwich(C_bond, chi_c_w))


def _fixture(nmat=8, Nx=4, Ny=4, Nz=1):
    """Green + bond topology + bare/dressed vertices for the U=4, NN-V=1
    single-band model on the (Nx,Ny,Nz) x nmat uniform grid."""
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, BETA)
    bset = _nn_bond_set()
    S_bond, C_bond, Vpp_s, Vpp_t = _vertices_full(bset, (Nx, Ny, Nz))
    chi_w = bc.bond_bubble_dynamic(green, bset, BETA)
    chi_s_w, chi_c_w, _, _ = bc.dress_bond_dynamic(chi_w, S_bond, C_bond)
    return dict(green=green, bset=bset, S_bond=S_bond, C_bond=C_bond,
                Vpp_s=Vpp_s, Vpp_t=Vpp_t,
                chi_s_w=chi_s_w, chi_c_w=chi_c_w,
                shape=(Nx, Ny, Nz), nmat=nmat)


def _random_gap(rng, Nx, Ny, Nz, nmat, norb=1):
    return (rng.normal(size=(norb, norb, Nx, Ny, Nz, nmat))
            + 1j * rng.normal(size=(norb, norb, Nx, Ny, Nz, nmat)))


# ---------------------------------------------------------------------------
# Step 1: oracle array equality -- the acceptance test
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_dynamic_bond_kernel_equals_oracle(pairing):
    Nx = Ny = 4
    Nz = 1
    nmat = 8
    fx = _fixture(nmat=nmat, Nx=Nx, Ny=Ny, Nz=Nz)
    green, bset = fx["green"], fx["bset"]

    A, vec_size = bc.make_bond_kernel_dynamic(
        fx["chi_s_w"], fx["chi_c_w"], fx["S_bond"], fx["C_bond"],
        fx["Vpp_s"], fx["Vpp_t"], green, bset, pairing, BETA)
    assert isinstance(A, LinearOperator)
    assert vec_size == 1 * 16 * nmat

    rng = np.random.default_rng(7)
    phi = _random_gap(rng, Nx, Ny, Nz, nmat)

    # Build the oracle's F_w from the SAME dressed chis: the kernel algebra
    # F = 1.5 S chi_s S - 0.5 C chi_c C is shared; what this test pins is the
    # transport (phases, q = k - k' convolution, tau-node frequency map).
    F_w = _oracle_F_w(pairing, fx["S_bond"], fx["C_bond"],
                      fx["chi_s_w"], fx["chi_c_w"])
    Vpp = fx["Vpp_s"] if pairing == "singlet" else fx["Vpp_t"]

    Y_oracle = oracle_bond_matvec(phi, green, F_w, Vpp, bset.delta_r, BETA)
    Y_prod = (A @ phi.ravel()).reshape(phi.shape)
    np.testing.assert_allclose(Y_prod, Y_oracle, rtol=1e-9, atol=1e-11)


# ---------------------------------------------------------------------------
# Step 5: B = 1 reduction and the instantaneous-term definition
# ---------------------------------------------------------------------------
def test_dynamic_bond_kernel_B1_reduces_to_scalar_kernel():
    """With a single (on-site, Delta r = 0) bond channel the enlarged index
    collapses to ND = nd = 1, the bond phases are all 1, and the fluctuation
    part must reproduce the pre-existing scalar dynamic kernel
    ``eliashberg_dynamic.eliashberg_kernel_dynamic`` exactly (spec 3.5).

    The scalar vertex is injected through the production signature by taking
    ``S_bond = 1``, ``C_bond = 0`` (so ``F = 1.5 chi_s`` for the singlet
    channel) and ``Vpp = 0`` (the scalar kernel has no instantaneous piece).
    """
    Nx = Ny = 4
    Nz = 1
    nmat = 8
    rng = np.random.default_rng(11)
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, BETA)
    bset = bc.resolve_interactions({}, np.eye(3), 1)
    assert bset.n_channels == 1 and bset.delta_r == ((0, 0, 0),)

    # Arbitrary bosonic scalar "chi" in the reality class F(q, j) =
    # conj(F(-q, -j)) -- the same class the oracle/scalar comparison in
    # tests/test_oracle_bond_dynamic.py uses.
    chi = (rng.normal(size=(Nx, Ny, Nz, nmat))
           + 1j * rng.normal(size=(Nx, Ny, Nz, nmat)))
    chi = 0.5 * (chi + np.conj(np.roll(
        chi[::-1, ::-1, ::-1, ::-1], (1, 1, 1, 1), (0, 1, 2, 3))))
    chi_s_w = chi[..., None, None]
    chi_c_w = np.zeros_like(chi_s_w)
    S_bond = np.ones((Nx, Ny, Nz, 1, 1), dtype=complex)
    C_bond = np.zeros((Nx, Ny, Nz, 1, 1), dtype=complex)
    Vpp0 = np.zeros((1, 1), dtype=complex)

    A, vec_size = bc.make_bond_kernel_dynamic(
        chi_s_w, chi_c_w, S_bond, C_bond, Vpp0, Vpp0,
        green, bset, "singlet", BETA)
    assert vec_size == Nx * Ny * Nz * nmat

    phi = _random_gap(rng, Nx, Ny, Nz, nmat)
    Y_bond = (A @ phi.ravel()).reshape(phi.shape)

    Vs_q_w = (1.5 * chi)[None, None, None, None]
    G2_w = ed.calc_g2_dynamic(green, BETA)
    Y_scalar = ed.eliashberg_kernel_dynamic(Vs_q_w, G2_w, phi, 1, BETA)

    np.testing.assert_allclose(Y_bond, Y_scalar, rtol=1e-10, atol=1e-12)


def test_instantaneous_term_equals_truncated_matsubara_sum():
    """With the fluctuation vertex switched off (chi_s = chi_c = 0 -> F_w = 0)
    the kernel output must be exactly the instantaneous (bare Cooper) term:
    flat in the external fermionic frequency, and equal to

        Y_pp(k) = -1/2 sum_{m,m'} Vpp[(m,ab),(m',cd)] e^{i k.dr_m'} A_{m;cd},
        A_{m;cd} = (T/N) sum_{k',n'} e^{-i k'.dr_m} X_cd(k', n'),

    i.e. the WINDOW (truncated) Matsubara sum of the pair bubble -- spec 3.3's
    normative definition on the uniform axis. ``X`` is taken from the oracle
    (``oracle_pair_bubble_X``), so the whole right-hand side is computed
    independently of the production transform chain.
    """
    Nx = Ny = 4
    Nz = 1
    nmat = 8
    rng = np.random.default_rng(23)
    green = _symmetric_green(1, Nx, Ny, Nz, nmat, BETA)
    bset = _nn_bond_set()
    S_bond, C_bond, Vpp_s, Vpp_t = _vertices_full(bset, (Nx, Ny, Nz))
    B = bset.n_channels
    ND = B  # norb = 1 -> nd = 1
    zero_chi = np.zeros((Nx, Ny, Nz, nmat, ND, ND), dtype=complex)

    A_op, _ = bc.make_bond_kernel_dynamic(
        zero_chi, zero_chi, S_bond, C_bond, Vpp_s, Vpp_t,
        green, bset, "singlet", BETA)
    phi = _random_gap(rng, Nx, Ny, Nz, nmat)
    Y = (A_op @ phi.ravel()).reshape(phi.shape)

    # (a) flat in the external frequency
    scale = np.abs(Y).max()
    assert scale > 0
    assert np.abs(Y - Y[..., :1]).max() <= 1e-12 * scale

    # (b) equals the independently-summed instantaneous expression
    X = oracle_pair_bubble_X(green, phi, BETA)      # (nd, Nx, Ny, Nz, nmat)
    N = Nx * Ny * Nz
    T = 1.0 / BETA
    kxg = 2.0 * np.pi * np.arange(Nx) / Nx
    kyg = 2.0 * np.pi * np.arange(Ny) / Ny
    kzg = 2.0 * np.pi * np.arange(Nz) / Nz
    KX, KY, KZ = np.meshgrid(kxg, kyg, kzg, indexing="ij")
    PH = np.empty((B, Nx, Ny, Nz), dtype=complex)
    for m in range(B):
        d = bset.delta_r[m]
        PH[m] = np.exp(1j * (KX * d[0] + KY * d[1] + KZ * d[2]))
    X_sum = X.sum(axis=-1)                          # (nd, Nx, Ny, Nz)
    A_coef = (T / N) * np.einsum('mxyz,cxyz->mc', np.conj(PH), X_sum)
    Vpp6 = np.asarray(Vpp_s).reshape(B, 1, B, 1)
    Bcoef = np.einsum('aAbB,aB->bA', Vpp6, A_coef)  # (B, nd)
    Y_ref_k = -0.5 * np.einsum('mxyz,ma->xyza', PH, Bcoef)  # (Nx,Ny,Nz,nd)
    Y_ref = np.broadcast_to(
        np.moveaxis(Y_ref_k, 3, 0).reshape(1, 1, Nx, Ny, Nz)[..., None],
        (1, 1, Nx, Ny, Nz, nmat))
    np.testing.assert_allclose(Y, Y_ref, rtol=1e-10, atol=1e-12)


# ---------------------------------------------------------------------------
# Step 6: production Hermiticity under the metric + mutation battery
# ---------------------------------------------------------------------------
def _dense_from_operator(A, vec_size):
    M = np.zeros((vec_size, vec_size), dtype=complex)
    for i in range(vec_size):
        e = np.zeros(vec_size, dtype=complex)
        e[i] = 1.0
        M[:, i] = A @ e
    return M


@pytest.mark.parametrize("pairing", ["singlet", "triplet"])
def test_production_kernel_hermitian_under_metric(pairing):
    """P0-1 on the PRODUCTION operator (the slow gate in
    tests/test_bond_dynamic_hermiticity.py established this for the oracle;
    densifying the production kernel is cheap enough to run by default).

    The metric is the pair weight w(k, n) = G(k, iw_n) G(-k, -iw_n), and the
    whitened kernel sqrt(w) K / sqrt(w) must be Hermitian.
    """
    Nx = Ny = 4
    Nz = 1
    nmat = 8
    fx = _fixture(nmat=nmat, Nx=Nx, Ny=Ny, Nz=Nz)
    green = fx["green"]
    A, vec_size = bc.make_bond_kernel_dynamic(
        fx["chi_s_w"], fx["chi_c_w"], fx["S_bond"], fx["C_bond"],
        fx["Vpp_s"], fx["Vpp_t"], green, fx["bset"], pairing, BETA)

    w_complex = (green * np.roll(green[:, :, ::-1, ::-1, ::-1, ::-1],
                                 (1, 1, 1), (2, 3, 4)))[0, 0]
    assert np.abs(w_complex.imag).max() <= 1e-10 * np.abs(w_complex).max()
    w = w_complex.real
    assert w.min() >= 1e-12 * w.max()

    K = _dense_from_operator(A, vec_size)
    sqw = np.sqrt(w).ravel()      # gap-vector order = (k, n) raveled like phi
    Kt = (sqw[:, None] * K) / sqw[None, :]
    res = np.linalg.norm(Kt - Kt.conj().T) / np.linalg.norm(Kt)
    assert res <= 1e-10, "production kernel not Hermitian: residual {:.3e}".format(res)


@pytest.mark.parametrize("mutation", ["vpp_sign", "fock_zero",
                                      "phase_conj", "transpose_mm"])
def test_mutations_break_oracle_equality(mutation):
    """Discriminating power of the Step-1 oracle equality (spec 3.5).

    Each mutation corrupts a PRODUCTION input only -- the oracle side keeps
    the pristine vertices -- and the comparison must then FAIL. A mutation
    that slipped through would mean the acceptance test does not actually
    constrain that piece of the transport. Measured relative deviations on
    this fixture: vpp_sign 2.6e-1, fock_zero 1.6e-1, phase_conj 7.8e-2,
    transpose_mm 7.8e-2 -- all ~1e8 times the 1e-9 acceptance tolerance.

    NOTE: on THIS fixture ``phase_conj`` and ``transpose_mm`` are degenerate
    (identical deviation, not a coincidence). The bond set is reversal-closed
    and the model is norb=1/centrosymmetric, so delta_r -> -delta_r merely
    permutes the channel labels by ``bond_set.reverse``, and the bubble obeys
    chi_{rev(m),rev(m')} = chi_{m',m} (its phase is e^{i k.(dr_m - dr_m')}).
    The two therefore probe the same degree of freedom here; both are kept
    because they stop being degenerate for a non-centrosymmetric or
    multi-orbital bond set.
    """
    Nx = Ny = 4
    Nz = 1
    nmat = 8
    pairing = "singlet"
    fx = _fixture(nmat=nmat, Nx=Nx, Ny=Ny, Nz=Nz)
    green, bset = fx["green"], fx["bset"]
    S_bond, C_bond = fx["S_bond"], fx["C_bond"]
    chi_s_w, chi_c_w = fx["chi_s_w"], fx["chi_c_w"]
    Vpp_s, Vpp_t = np.asarray(fx["Vpp_s"]), np.asarray(fx["Vpp_t"])
    B = bset.n_channels

    # pristine reference (oracle side)
    F_w = _oracle_F_w(pairing, S_bond, C_bond, chi_s_w, chi_c_w)
    rng = np.random.default_rng(7)
    phi = _random_gap(rng, Nx, Ny, Nz, nmat)
    Y_oracle = oracle_bond_matvec(phi, green, F_w, Vpp_s, bset.delta_r, BETA)

    # --- corrupt the production inputs -----------------------------------
    p_S, p_C = S_bond, C_bond
    p_chi_s, p_chi_c = chi_s_w, chi_c_w
    p_Vpp_s, p_Vpp_t = Vpp_s, Vpp_t
    p_bset = bset

    if mutation == "vpp_sign":
        p_Vpp_s, p_Vpp_t = -Vpp_s, -Vpp_t
    elif mutation == "fock_zero":
        # zero the Fock (m >= 1 bond-diagonal) entries of the bare vertices
        p_S, p_C = S_bond.copy(), C_bond.copy()
        for m in range(1, B):
            p_S[:, :, :, m, m] = 0.0
            p_C[:, :, :, m, m] = 0.0
    elif mutation == "phase_conj":
        # PH -> conj(PH) via delta_r -> -delta_r on the production side only
        p_bset = dataclasses.replace(
            bset,
            delta_r=tuple(tuple(-c for c in d) for d in bset.delta_r))
    elif mutation == "transpose_mm":
        # transpose the F blocks in (m, m') (norb = 1, so the enlarged index
        # IS the channel index and the block transpose is a plain transpose)
        p_chi_s = np.swapaxes(chi_s_w, -1, -2).copy()
        p_chi_c = np.swapaxes(chi_c_w, -1, -2).copy()
    else:  # pragma: no cover - parametrization is closed
        raise AssertionError(mutation)

    A, _ = bc.make_bond_kernel_dynamic(
        p_chi_s, p_chi_c, p_S, p_C, p_Vpp_s, p_Vpp_t,
        green, p_bset, pairing, BETA)
    Y_prod = (A @ phi.ravel()).reshape(phi.shape)

    np.testing.assert_raises(
        AssertionError, np.testing.assert_allclose,
        Y_prod, Y_oracle, rtol=1e-9, atol=1e-11)
