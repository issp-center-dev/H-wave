"""Density matrix consistency check against H-wave's _UHF_cisajs.dat.

G_{ij} = <c^dag_i c_j> in H-wave's convention is
   G^sigma_{ij} = sum_alpha (A^sigma)^*_{i, alpha} A^sigma_{j, alpha}
              = (A^sigma.conj() @ A^sigma.T)_{ij}

See docs/en/source/algorithm/uhfk_to_mvmc.rst.
"""
from __future__ import annotations

import numpy as np


class DensityMismatchError(RuntimeError):
    """Raised when bridge-built G^sigma_{ij} disagrees with H-wave's
    _UHF_cisajs.dat beyond the given tolerance."""


def density_from_amplitudes(A):
    """G_{ij} = sum_alpha conj(A_{i, alpha}) A_{j, alpha}."""
    A = np.asarray(A, dtype=np.complex128)
    return np.conj(A) @ A.T


def parse_uhf_cisajs_dat(path):
    """Parse a _UHF_cisajs.dat text file. Returns list of
    (i, s, j, t, complex value)."""
    out = []
    with open(path) as fp:
        for ln in fp:
            ln = ln.strip()
            if not ln:
                continue
            toks = ln.split()
            i, s, j, t = (int(toks[0]), int(toks[1]),
                          int(toks[2]), int(toks[3]))
            v = float(toks[4]) + 1j * float(toks[5])
            out.append((i, s, j, t, v))
    return out


def compare_against_onebodyg_uhf(
    G_up, G_down, onebodyg_uhf_path, tol=1e-10
):
    """Compare bridge-built G^up, G^down against H-wave's _UHF_cisajs.dat.

    Line by line, iterate (i, s, j, t) entries; the bridge
    side returns G^s_{ij} (s == t required, otherwise 0 for Sz-fixed).
    """
    entries = parse_uhf_cisajs_dat(onebodyg_uhf_path)
    diffs = []
    for (i, s, j, t, v_uhf) in entries:
        if s != t:
            v_bridge = 0.0 + 0.0j
        elif s == 0:
            v_bridge = complex(G_up[i, j])
        elif s == 1:
            v_bridge = complex(G_down[i, j])
        else:
            raise DensityMismatchError(
                f"unexpected spin index s={s} in {onebodyg_uhf_path}"
            )
        if abs(v_bridge - v_uhf) > tol:
            diffs.append((i, s, j, t, v_bridge, v_uhf))
    if diffs:
        head = diffs[:3]
        raise DensityMismatchError(
            f"{len(diffs)} (i, s, j, t) entries differ beyond tol={tol}; "
            f"first 3: {head}"
        )


def compare_against_green_sublattice(
    G_all: np.ndarray,
    green_path: str,
    site_positions: np.ndarray = None,
    cell_shape: np.ndarray = None,
    subshape: np.ndarray = None,
    norb_orig: int = None,
    tol: float = 1e-10,
    is_soc_sublattice_mode: bool = False,
    Ns: int = None,
    boundary_theta: tuple = (0.0, 0.0, 0.0),
) -> None:
    """Compare bridge-built physical G against H-wave's folded
    ``green_sublattice``.

    Under SOC + SubShape > [1, 1, 1], H-wave's ``greenone.dat`` is
    not used as the density-check reference. This function folds the
    bridge's physical ``G_all`` back to the supercell basis and compares
    element-wise against ``green.npz['green_sublattice']``.

    Fold formula (matches H-wave's ``_reshape_orbit_spin`` encoding):

        G_sub[R_diff]_{aa, bb} = G_all[all_i, all_j]

    where for every physical pair (r_i, r_j):
      R_diff = (r_j // B - r_i // B) mod (L / B)
      aa     = spin_i + norb_orig * within_flat(r_i % B)
      bb     = spin_j + norb_orig * within_flat(r_j % B)
      all_i  = r_i_ordinal + spin_i * Ns
      all_j  = r_j_ordinal + spin_j * Ns

    Supercell translation invariance in the physical G means every
    (r_i, r_j) pair mapping to a given (R_diff, aa, bb) class must give
    the same value; a residual > ``tol`` is a bridge bug.

    When ``is_soc_sublattice_mode=True``, the check
    switches from folding G_all back to the supercell basis to instead
    LIFTING ``green_sublattice`` to the physical basis via ``gauge_lift``
    and comparing to ``G_all`` element-wise on the full 2Ns x 2Ns block.
    This gauge-invariant lift provides a strict density gate for the
    SOC + SubShape > [1, 1, 1] path. Requires ``cell_shape``, ``subshape``,
    ``site_positions``; ``Ns`` defaults to ``site_positions.shape[0]``.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst for the gauge
    composition and physical-density convention.
    """
    G_all = np.asarray(G_all, dtype=np.complex128)

    if is_soc_sublattice_mode:
        # Apply the documented SOC + SubShape gauge lift; see
        # docs/en/source/algorithm/uhfk_to_mvmc.rst.
        assert (
            cell_shape is not None
            and subshape is not None
            and site_positions is not None
        ), (
            "SOC+SubShape density check requires cell_shape, subshape, "
            "site_positions."
        )
        cell_shape_arr = np.asarray(cell_shape, dtype=np.int64)
        subshape_arr = np.asarray(subshape, dtype=np.int64)
        site_positions_arr = np.asarray(site_positions, dtype=np.int64)
        Ns_soc = int(Ns) if Ns is not None else int(site_positions_arr.shape[0])
        assert G_all.shape == (2 * Ns_soc, 2 * Ns_soc), (
            f"G_all shape {G_all.shape} does not match (2*Ns, 2*Ns) = "
            f"({2 * Ns_soc}, {2 * Ns_soc})"
        )
        data = np.load(green_path, allow_pickle=False)
        if "green_sublattice" not in data.files:
            raise DensityMismatchError(
                f"green.npz {green_path!r} has no green_sublattice key; "
                "H-wave writes it only under has_sublattice = True"
            )
        green_sublattice = data["green_sublattice"]
        folded_cell_of = lambda r_phys: r_phys // subshape_arr
        G_ref = np.zeros((2 * Ns_soc, 2 * Ns_soc), dtype=np.complex128)
        for site_i in range(Ns_soc):
            for spin_i in (0, 1):
                for site_j in range(Ns_soc):
                    for spin_j in (0, 1):
                        all_i = site_i + spin_i * Ns_soc
                        all_j = site_j + spin_j * Ns_soc
                        G_ref[all_i, all_j] = gauge_lift(
                            green_sublattice,
                            site_i, spin_i, site_j, spin_j,
                            subshape=subshape_arr,
                            cell_shape=cell_shape_arr,
                            site_positions=site_positions_arr,
                            folded_cell_of=folded_cell_of,
                            boundary_theta=boundary_theta,
                        )
        max_diff = float(np.max(np.abs(G_ref - G_all)))
        if max_diff > tol:
            raise DensityMismatchError(
                f"SOC+SubShape gauge-lifted density mismatch: "
                f"|G_ref - G_ship|_max = {max_diff:.3e} > tol = {tol:.3e}"
            )
        return

    # The non-SOC sublattice path uses the folded-Green comparison.
    assert (
        site_positions is not None
        and cell_shape is not None
        and subshape is not None
        and norb_orig is not None
    ), (
        "Non-SOC folded-Green density check requires site_positions, "
        "cell_shape, subshape, norb_orig."
    )
    site_positions = np.asarray(site_positions, dtype=np.int64)
    cell_shape = np.asarray(cell_shape, dtype=np.int64)
    subshape = np.asarray(subshape, dtype=np.int64)
    Ns = site_positions.shape[0]
    assert G_all.shape == (2 * Ns, 2 * Ns)
    L_folded = cell_shape // subshape
    Nx, Ny, Nz = int(L_folded[0]), int(L_folded[1]), int(L_folded[2])
    Bx, By, Bz = int(subshape[0]), int(subshape[1]), int(subshape[2])
    subvol = Bx * By * Bz
    norb_folded = int(norb_orig) * subvol
    nvol_folded = Nx * Ny * Nz

    data = np.load(green_path, allow_pickle=False)
    if "green_sublattice" not in data.files:
        raise DensityMismatchError(
            f"green.npz {green_path!r} has no green_sublattice key; "
            "H-wave writes it only under has_sublattice = True"
        )
    green_sub = data["green_sublattice"]  # (nvol, ns, norb, ns, norb)
    # SOC uses ns=1; extract the (nvol, norb, norb) 2-D block.
    if green_sub.ndim != 5 or green_sub.shape[1] != 1 or green_sub.shape[3] != 1:
        raise DensityMismatchError(
            f"green.npz green_sublattice shape {green_sub.shape} does not "
            "match expected (nvol, 1, norb, 1, norb) for SOC mode"
        )
    green_sub_2d = green_sub[:, 0, :, 0, :]  # (nvol_folded, norb_folded, norb_folded)
    if green_sub_2d.shape != (nvol_folded, norb_folded, norb_folded):
        raise DensityMismatchError(
            f"green.npz green_sublattice folded shape "
            f"{green_sub_2d.shape} does not match expected "
            f"({nvol_folded}, {norb_folded}, {norb_folded}) derived from "
            f"CellShape={cell_shape.tolist()}, SubShape={subshape.tolist()}, "
            f"norb_orig={norb_orig}"
        )

    def _w_flat(w):
        return int(w[0] + Bx * (w[1] + By * w[2]))

    def _R_pack(R):
        # Matches H-wave uhfk._pack_site: R_z + Nz*(R_y + Ny*R_x).
        return int(R[2] + Nz * (R[1] + Ny * R[0]))

    diffs = []
    max_translation_residual = 0.0
    seen_classes = {}
    for i in range(Ns):
        r_i = site_positions[i]
        R_i = r_i // subshape
        w_i = r_i % subshape
        for j in range(Ns):
            r_j = site_positions[j]
            R_j = r_j // subshape
            w_j = r_j % subshape
            R_diff = (R_j - R_i) % L_folded
            idx_R = _R_pack(R_diff)
            for spin_i in (0, 1):
                for spin_j in (0, 1):
                    aa = spin_i + norb_orig * _w_flat(w_i)
                    bb = spin_j + norb_orig * _w_flat(w_j)
                    all_i = i + spin_i * Ns
                    all_j = j + spin_j * Ns
                    v_bridge = complex(G_all[all_i, all_j])
                    v_ref = complex(green_sub_2d[idx_R, aa, bb])
                    key = (idx_R, aa, bb)
                    prev = seen_classes.get(key)
                    if prev is not None:
                        d = abs(v_bridge - prev)
                        if d > max_translation_residual:
                            max_translation_residual = d
                    else:
                        seen_classes[key] = v_bridge
                    d = abs(v_bridge - v_ref)
                    if d > tol:
                        diffs.append(
                            (i, spin_i, j, spin_j, v_bridge, v_ref)
                        )
    if diffs:
        head = diffs[:3]
        raise DensityMismatchError(
            f"{len(diffs)} (i, s, j, t) entries differ beyond tol={tol} "
            f"[SOC + SubShape > 1 folded-Green check]; first 3: {head}. "
            f"Bridge G supercell translation invariance max residual = "
            f"{max_translation_residual:.3e}"
        )


def gauge_lift(green_sublattice, site_i, spin_i, site_j, spin_j,
               subshape, cell_shape, site_positions, folded_cell_of,
               boundary_theta):
    """Reconstruct G_phys[all_i, all_j] independently from green_sublattice.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst for the sublattice and
    APBC gauge composition. Convention: G[all_i, all_j] =
    sum_alpha conj(A[all_i, alpha]) * A[all_j, alpha], matching H-wave's
    _green (conj(V_a) * V_b) and the existing bridge density check.

    Independence contract: MUST NOT call build_slater_orbitals or any
    shipping-A helper. Derives purely from green_sublattice + geometry.

    Parameters
    ----------
    green_sublattice : (nvol, 1, 2*Ns_folded, 1, 2*Ns_folded) complex
        Raw H-wave 5D green.npz["green_sublattice"] under SOC mode.
    site_i, site_j : int
        Physical site indices in [0, Ns).
    spin_i, spin_j : int
        Physical spin indices in {0, 1}.
    subshape, cell_shape : (3,) int
        SubShape and CellShape from input.toml.
    site_positions : (Ns, 3) int
        Physical site coordinates from geometry_uhf.
    folded_cell_of : callable(r_phys) -> (3,) int
        Returns the folded-cell index for a physical position (r_phys //
        subshape componentwise).
    boundary_theta : array-like of length 3
        Twist ``(theta_x, theta_y, theta_z)`` in RADIANS. Each component
        in ``{0, pi}`` under Periodic / Antiperiodic. Passed WITHOUT
        conversion; callers MUST pass radians, NOT the dimensionless
        ``twist_offset = theta / (2*pi)``. The composed phase is
        ``exp(-i k_folded . dr_folded) * exp(-i theta . dr_full / L_full)``
        where ``dr_full = r_phys_j - r_phys_i`` and
        ``L_full = SubShape * L_folded``. Under PBC (all-zero theta) the
        twist factor collapses to 1 and the output is bit-identical to
        the PBC result.

    Returns
    -------
    G_phys : complex
        The physical density matrix element G[all_i, all_j] where
        all_i = site_i + spin_i * Ns, all_j = site_j + spin_j * Ns.
    """
    import numpy as np
    assert green_sublattice.ndim == 5, (
        f"green_sublattice must be H-wave 5D shape; got {green_sublattice.shape}"
    )
    gs_soc = green_sublattice[:, 0, :, 0, :]
    L_folded = cell_shape // subshape
    gs = gs_soc.reshape(
        L_folded[0], L_folded[1], L_folded[2], gs_soc.shape[1], gs_soc.shape[2],
    )
    G_k = np.fft.ifftn(gs, axes=(0, 1, 2), norm="forward")

    # H-wave folded SOC row: 2 * folded_orb + spin under enable_spin_orbital.
    # This path supports one physical orbital (norb_orig==1); folded_orb is the
    # flat sublattice index within the folded cell.
    r_phys_i = np.asarray(site_positions[site_i], dtype=np.int64)
    r_phys_j = np.asarray(site_positions[site_j], dtype=np.int64)
    fc_i = folded_cell_of(r_phys_i)
    fc_j = folded_cell_of(r_phys_j)
    so_i = r_phys_i - fc_i * subshape  # sub_offset(i)
    so_j = r_phys_j - fc_j * subshape
    folded_orb_i = int(so_i[0] + subshape[0] * (so_i[1] + subshape[1] * so_i[2]))
    folded_orb_j = int(so_j[0] + subshape[0] * (so_j[1] + subshape[1] * so_j[2]))
    aa = 2 * folded_orb_i + int(spin_i)
    bb = 2 * folded_orb_j + int(spin_j)

    # The twist gauge is k-independent per (i, j) pair, so hoist
    # it out of the k loop. Under PBC (all-zero theta) phase_twist == 1
    # and the result reduces exactly to the PBC formula. See
    # docs/en/source/algorithm/uhfk_to_mvmc.rst.
    r_diff_folded = (fc_j.astype(np.float64) + so_j.astype(np.float64)) \
                    - (fc_i.astype(np.float64) + so_i.astype(np.float64))
    r_diff_full = r_phys_j.astype(np.float64) - r_phys_i.astype(np.float64)
    L_full = (subshape * L_folded).astype(np.float64)
    theta = np.asarray(boundary_theta, dtype=np.float64)
    phase_twist = np.exp(-1j * np.dot(theta, r_diff_full / L_full))

    accum = 0.0j
    for kx in range(L_folded[0]):
        for ky in range(L_folded[1]):
            for kz in range(L_folded[2]):
                k_vec = 2.0 * np.pi * np.array([kx, ky, kz], dtype=np.float64) / L_folded
                phase_folded = np.exp(-1j * np.dot(k_vec, r_diff_folded))
                accum += G_k[kx, ky, kz, aa, bb] * phase_folded
    accum *= phase_twist
    # Divide by N_folded to match H-wave's forward-normalized fftn
    # convention; see docs/en/source/algorithm/uhfk_to_mvmc.rst.
    return accum / float(np.prod(L_folded))


def compare_against_onebodyg_uhf_general(
    G_all: np.ndarray, onebodyg_uhf_path: str, tol: float = 1e-10,
    is_soc_mode: bool = False,
) -> None:
    """Compare bridge-built 2Ns × 2Ns G against H-wave's greenone.dat
    for the General path.

    G_all[iσ, jσ'] = <c^†_{i,σ} c_{j,σ'}> in the physical basis with
    the mVMC spin-block index ``all_i = i + σ * Ns`` (site-major,
    spin-minor). greenone.dat lines are ``i s j t re im`` and are
    compared element-wise to ``G_all[i + s*Ns, j + t*Ns]``.

    Behavior by mode:

    ``is_soc_mode=False`` (Sz-diagonal path)
        The bridge produces a spin-block-diagonal ``G_all`` because
        amplitudes are Sz-fixed. Reference greenone.dat rows with
        s == t compare element-wise. Rows with s != t are also
        compared element-wise; because ``G_all`` is zero on
        off-diagonal blocks, any non-zero s != t reference value
        raises ``DensityMismatchError``. This element-wise comparison
        is the mixed-block scope-violation guard: it catches the case
        where a caller feeds SOC-tainted reference data into the General
        path.

    ``is_soc_mode=True`` (SOC path)
        Both G_all and the reference greenone.dat carry physically
        non-zero s != t entries. Every row is compared under the same
        `tol` regardless of (s, t); the scope-violation framing does
        not apply because the SOC bridge deliberately populates the
        off-diagonal spin blocks.

    Both modes share the same comparison kernel: uniform element-wise
    match under `tol`. The flag documents intent and gates any future
    mode-specific diagnostics without changing comparison behavior.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst for the physical
    density orientation and spin-block indexing convention.
    """
    G_all = np.asarray(G_all, dtype=np.complex128)
    two_ns = G_all.shape[0]
    assert two_ns % 2 == 0, "G_all must have even row count"
    Ns = two_ns // 2

    entries = parse_uhf_cisajs_dat(onebodyg_uhf_path)
    diffs = []
    for (i, s, j, t, v_uhf) in entries:
        all_i = i + s * Ns
        all_j = j + t * Ns
        v_bridge = complex(G_all[all_i, all_j])
        if abs(v_bridge - v_uhf) > tol:
            diffs.append((i, s, j, t, v_bridge, v_uhf))
    if diffs:
        head = diffs[:3]
        mode_label = "SOC" if is_soc_mode else "Sz-diagonal"
        raise DensityMismatchError(
            f"{len(diffs)} (i, s, j, t) entries differ beyond tol={tol} "
            f"[{mode_label} mode]; first 3: {head}"
        )
