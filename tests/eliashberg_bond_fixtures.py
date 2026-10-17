"""Physical random fixtures for the bond pairing kernel tests (spec 10.1):
a k-even Hermitian hopping, the bare Green function, the bond bubble from
bubble._iter_bond_dynamic, S/C from build_sc_bond_channel and the dressed
chi from flex_bond._dress. Everything is reversal-closed by construction."""
import numpy as np

from hwave.solver import bond_channels, bubble, flex_bond, offsite
from hwave.solver.eliashberg_bond import ArrayBlockSource


def _k_even_hopping(norb, shape, rng):
    """h(k) = h(-k) Hermitian: real symmetric hopping on the first shell plus a
    random Hermitian on-site block."""
    nx, ny, nz = shape
    kx = 2 * np.pi * np.arange(nx) / nx
    ky = 2 * np.pi * np.arange(ny) / ny
    kz = 2 * np.pi * np.arange(nz) / nz
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    disp = -2.0 * (np.cos(KX) + np.cos(KY) + (np.cos(KZ) if nz > 1 else 0.0))
    onsite = rng.standard_normal((norb, norb)) + 1j * rng.standard_normal((norb, norb))
    onsite = 0.5 * (onsite + onsite.conj().T)
    t_orb = np.diag(0.5 + rng.random(norb))
    h = disp[..., None, None] * t_orb + onsite
    return h.reshape(nx * ny * nz, norb, norb)


def _green(h_k, beta, nmat, mu=0.0):
    nvol, norb, _ = h_k.shape
    w = (2 * np.arange(nmat) + 1 - nmat) * np.pi / beta
    g = np.empty((1, nmat, nvol, norb, norb), complex)
    eye = np.eye(norb)
    for n in range(nmat):
        g[0, n] = np.linalg.inv((1j * w[n] + mu) * eye[None] - h_k)
    return g


def _tables_from_pairs(inter, U, norb):
    """``{(type, R): (norb, norb) real matrix}`` (the fixture's own argument
    form) plus an on-site CoulombIntra ``U`` -> the ``read_w90`` dict-of-dicts
    ``{type: {(irvec, (a, b)): value}}`` that ``resolve_bond_topology`` and
    ``offsite.split_locality`` both speak. Returns ``(whole, onsite, offsite)``.
    """
    onsite_tbl = {"CoulombIntra": {((0, 0, 0), (a, a)): complex(U) for a in range(norb)}}
    offsite_tbl = {}
    for (t, R), mat in inter.items():
        mat = np.asarray(mat)
        tbl = offsite_tbl.setdefault(t, {})
        for a in range(norb):
            for b in range(norb):
                if mat[a, b] != 0:
                    tbl[(tuple(int(x) for x in R), (a, b))] = complex(mat[a, b])
    whole_tbl = {t: dict(v) for t, v in onsite_tbl.items()}
    for t, tbl in offsite_tbl.items():
        whole_tbl.setdefault(t, {}).update(tbl)
    return whole_tbl, onsite_tbl, offsite_tbl


def _split_from_tables(whole_tbl, onsite_tbl, offsite_tbl):
    """A :class:`hwave.solver.offsite.LocalitySplit` over already-split test
    tables (no folding), so that the fixture can call the FLEX gate's own
    ``sc_matrices_from_split`` builder. Test-only: the production split is
    always built by ``offsite.split_locality`` from a solver's ``ham_info``.
    """
    return offsite.LocalitySplit(
        whole_tbl=whole_tbl, onsite_tbl=onsite_tbl, offsite_tbl=offsite_tbl,
        offsite_prefold_tbl=offsite_tbl, offsite_types=tuple(offsite_tbl),
        has_fold=False)


def _channel0_from_gate(whole_tbl, onsite_tbl, offsite_tbl, norb, shape):
    """S0, C0 (nvol, nd, nd) through offsite.sc_matrices_from_split, the FLEX
    gate's own builder, from a minimal locality split."""
    nx, ny, nz = shape
    nd = norb * norb
    split = _split_from_tables(whole_tbl, onsite_tbl, offsite_tbl)
    S0, C0 = offsite.sc_matrices_from_split(
        split, bond_channels._LONGITUDINAL_ACTIVE_TYPES, norb, nx, ny, nz)
    return (S0.reshape(nx * ny * nz, nd, nd).copy(),
            C0.reshape(nx * ny * nz, nd, nd).copy())


def physical_fixture(norb=1, shape=(4, 4, 1), nmat=8, beta=2.0,
                     delta_r=((0, 0, 0), (1, 0, 0), (-1, 0, 0)),
                     coeffs=None, seed=0, U=1.0):
    """coeffs: {(type, R): (norb, norb) real} declared off-site couplings; default
    CoulombInter V = 0.5 on every non-zero delta_r, plus CoulombIntra U on channel 0."""
    rng = np.random.default_rng(seed)
    nx, ny, nz = shape
    nvol = nx * ny * nz
    nd = norb * norb
    h_k = _k_even_hopping(norb, shape, rng)
    green = _green(h_k, beta, nmat)
    inter = {}
    for R in delta_r:
        if tuple(R) != (0, 0, 0):
            inter[("CoulombInter", tuple(R))] = 0.5 * np.ones((norb, norb))
    if coeffs:
        inter.update(coeffs)
    whole_tbl, onsite_tbl, offsite_tbl = _tables_from_pairs(inter, U, norb)
    # the topology API of the FLEX gate takes the read_w90 dict-of-dicts
    topo = bond_channels.resolve_bond_topology(
        whole_tbl, np.eye(3), norb,
        active_types=bond_channels._LONGITUDINAL_ACTIVE_TYPES)
    view = bond_channels.BondSetView(topo)
    B = view.n_channels
    ND = B * nd
    # channel-0 Kuroki matrices for CoulombIntra U (+ the V(q) Hartree the gate
    # builder places in the charge (aa,bb) element; the fixture takes the
    # FLEX gate's own builder to stay bit-consistent with the solver)
    S0, C0 = _channel0_from_gate(whole_tbl, onsite_tbl, offsite_tbl, norb, shape)
    S = bond_channels.build_sc_bond_channel(topo, S0, "S", types=tuple(topo.coeffs))
    C = bond_channels.build_sc_bond_channel(topo, C0, "C", types=tuple(topo.coeffs))
    chibar = np.zeros((nmat, nvol, ND, ND), complex)
    for (m, mp), block in bubble._iter_bond_dynamic(green, None, beta, view,
                                                    spatial_shape=tuple(shape)):
        chibar[:, :, m * nd:(m + 1) * nd, mp * nd:(mp + 1) * nd] = block
    chi_s = np.empty_like(chibar)
    chi_c = np.empty_like(chibar)
    nb = max(1, nmat // 2)
    for l0 in range(0, nmat, nb):
        l1 = min(nmat, l0 + nb)
        chi_s[l0:l1] = flex_bond._dress(chibar[l0:l1], S, "spin", l0, nmat,
                                        tuple(shape), 1e-6, None)[0]
        chi_c[l0:l1] = flex_bond._dress(chibar[l0:l1], C, "charge", l0, nmat,
                                        tuple(shape), 1e-6, None)[0]
    store = ArrayBlockSource({"chibar": chibar, "chi_s_w": chi_s, "chi_c_w": chi_c}, nd)
    green_sc = green[0].reshape(nmat, nx, ny, nz, norb, norb).transpose(4, 5, 1, 2, 3, 0).copy()
    return dict(green=green, green_sc=green_sc, topo=topo, view=view, S=S, C=C, store=store,
                nd=nd, ND=ND, B=B, nvol=nvol, spatial_shape=tuple(shape), nmat=nmat, beta=beta,
                norb=norb, chibar=chibar, chi_s=chi_s, chi_c=chi_c, inter=inter, U=U,
                S0=S0, C0=C0, whole_tbl=whole_tbl, onsite_tbl=onsite_tbl,
                offsite_tbl=offsite_tbl)
