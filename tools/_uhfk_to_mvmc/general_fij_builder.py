"""InOrbitalGeneral F builder (Sz-fixed 2Sz≠0 A + Sz-free non-mixed B).

See docs/en/source/algorithm/uhfk_to_mvmc.rst for the construction.
"""
from __future__ import annotations

import numpy as np


def compute_canonical_reps(
    partner_rows: np.ndarray,
    wavevector_index: np.ndarray,
) -> tuple[list[int], list[int]]:
    """Pick one canonical row per unordered {k, partner(k)} pair.

    Rules (see docs/en/source/algorithm/uhfk_to_mvmc.rst):
      - self-pair k (partner_rows[k] == k): k is canonical.
      - non-self-pair {k, partner(k)}: canonical = row whose
        wavevector_index tuple is lexicographically smaller
        (fallback tie-break: smaller row index).

    Returns
    -------
    canonical_rows : list[int]
        Sorted list of canonical row indices; length =
        (#self-pairs) + (#non-self-pairs / 2).
    self_pair_rows : list[int]
        Subset of canonical_rows whose partner is themselves.
    """
    partner_rows = np.asarray(partner_rows, dtype=np.int64)
    wavevector_index = np.asarray(wavevector_index, dtype=np.int64)
    nvol = partner_rows.shape[0]

    canonical = []
    self_pairs = []
    seen = set()
    for k_row in range(nvol):
        if k_row in seen:
            continue
        partner = int(partner_rows[k_row])
        if partner == k_row:
            canonical.append(k_row)
            self_pairs.append(k_row)
            seen.add(k_row)
            continue
        # Non-self: pick lex-smaller wavevector_index tuple.
        wv_k = tuple(int(v) for v in wavevector_index[k_row])
        wv_p = tuple(int(v) for v in wavevector_index[partner])
        if wv_k < wv_p:
            chosen = k_row
        elif wv_p < wv_k:
            chosen = partner
        else:  # exact tie → smaller row index
            chosen = min(k_row, partner)
        canonical.append(chosen)
        seen.add(k_row)
        seen.add(partner)
    canonical.sort()
    return canonical, self_pairs


def _occ_by_spin(stepped_occupation, column_spin, spin_value):
    """Return dict[k_row] = ordered list of columns occupied at that k
    with the specified spin (0 = up, 1 = down)."""
    nvol, nd = stepped_occupation.shape
    columns = [c for c in range(nd) if int(column_spin[c]) == spin_value]
    out = {}
    for k_row in range(nvol):
        out[k_row] = [c for c in columns if stepped_occupation[k_row, c] >= 0.5]
    return out


def validate_general_prerequisites(
    Ncond: int,
    stepped_occupation: np.ndarray,
    column_spin: np.ndarray,
    partner_rows: np.ndarray,
    wavevector_index: np.ndarray,
    is_soc_mode: bool = False,
) -> None:
    """Fail-fast validation before General F construction.

    Raises ValueError with a diagnostic message on:
      - Any column_spin value not in {0, 1} when ``is_soc_mode`` is False
        (mixed block; unsupported on the non-SOC path)
      - Any column_spin value not in {-1, 0, 1} when ``is_soc_mode`` is True
        (SOC path accepts mixed-block sentinel -1)
      - Ncond odd (mVMC PairProduct requires even Ne)
      - sum(stepped_occupation) != Ncond
      - Non-self canonical (k, partner): n_excess_up_k != n_excess_up_p
        OR n_excess_down_k != n_excess_down_p (same-spin excess imbalance;
        non-SOC only)
      - Self-pair canonical k: n_excess_up_k odd OR n_excess_down_k odd
        (non-SOC only)
      - When ``is_soc_mode`` is True: non-self canonical
        (k, p) with n_occ(k) != n_occ(p) OR self canonical k with odd
        n_occ(k)

    ``is_soc_mode`` (default ``False``) is the SOC dispatch flag
    threaded through by the CLI. When True, the mixed-block column_spin
    guard is relaxed to accept the -1 sentinel and SOC-specific spin-
    agnostic occupation prerequisites (n_occ balance, self-pair parity)
    are enforced as upfront preconditions before pair emission. See
    docs/en/source/algorithm/uhfk_to_mvmc.rst for partner balance.
    """
    column_spin = np.asarray(column_spin, dtype=np.int64)
    stepped_occupation = np.asarray(stepped_occupation, dtype=np.float64)
    partner_rows = np.asarray(partner_rows, dtype=np.int64)

    if not is_soc_mode:
        if np.any((column_spin != 0) & (column_spin != 1)):
            offending = np.unique(
                column_spin[(column_spin != 0) & (column_spin != 1)]
            )
            raise ValueError(
                "column_spin contains mixed-block values "
                f"{offending.tolist()}; pass is_soc_mode=True for the "
                "SOC path"
            )
    else:
        # The SOC path accepts column_spin in {-1} for mixed block, or
        # {0, 1} for degenerate SOC states.
        allowed = np.isin(column_spin, [-1, 0, 1])
        if not np.all(allowed):
            offending = np.unique(column_spin[~allowed])
            raise ValueError(
                f"SOC path: column_spin values {offending.tolist()} not "
                "recognized (expected in {-1, 0, 1})"
            )
        # Enforce spin-agnostic SOC occupation prerequisites. This mirrors
        # the pair-emission-time guards in
        # ``build_pair_list`` but surfaces them as upfront preconditions.
        for k in range(stepped_occupation.shape[0]):
            p = int(partner_rows[k])
            occ_k = int(np.round(stepped_occupation[k].sum()))
            if p == k:
                if occ_k % 2 != 0:
                    raise ValueError(
                        f"SOC self canonical block k={k}: n_occ={occ_k} "
                        "is odd"
                    )
            else:
                occ_p = int(np.round(stepped_occupation[p].sum()))
                if occ_k != occ_p:
                    raise ValueError(
                        f"SOC canonical block (k={k}, p={p}) imbalance: "
                        f"n_occ(k)={occ_k}, n_occ(p)={occ_p}"
                    )
    if Ncond % 2 != 0:
        raise ValueError(
            f"Ncond={Ncond} is odd; mVMC PairProduct requires even Ne"
        )
    total_occ = float(np.sum(stepped_occupation))
    if abs(total_occ - Ncond) > 0.5:
        raise ValueError(
            f"sum(stepped_occupation)={total_occ} != Ncond={Ncond}"
        )

    occ_up = _occ_by_spin(stepped_occupation, column_spin, 0)
    occ_dn = _occ_by_spin(stepped_occupation, column_spin, 1)
    canonical, self_pairs = compute_canonical_reps(partner_rows, wavevector_index)

    for k in canonical:
        partner = int(partner_rows[k])
        nu_k, nd_k = len(occ_up[k]), len(occ_dn[k])
        if partner == k:
            n_cross = min(nu_k, nd_k)
            excess_up = nu_k - n_cross
            excess_dn = nd_k - n_cross
            if excess_up % 2 != 0:
                raise ValueError(
                    f"self-pair canonical row {k}: excess up-spin count "
                    f"{excess_up} is odd (must be even for same-spin pairing)"
                )
            if excess_dn % 2 != 0:
                raise ValueError(
                    f"self-pair canonical row {k}: excess down-spin count "
                    f"{excess_dn} is odd (must be even)"
                )
        else:
            nu_p, nd_p = len(occ_up[partner]), len(occ_dn[partner])
            n_cross_kd = min(nu_k, nd_p)
            n_cross_dk = min(nu_p, nd_k)
            excess_up_k = nu_k - n_cross_kd
            excess_up_p = nu_p - n_cross_dk
            excess_dn_k = nd_k - n_cross_dk
            excess_dn_p = nd_p - n_cross_kd
            if excess_up_k != excess_up_p:
                raise ValueError(
                    f"canonical block (k={k}, partner={partner}): "
                    f"n_excess_up_k={excess_up_k} != n_excess_up_p={excess_up_p} "
                    f"(#up@k={nu_k}, #down@partner={nd_p}, "
                    f"#up@partner={nu_p}, #down@k={nd_k}); "
                    "spin-cross + same-spin excess pair-closure violated"
                )
            if excess_dn_k != excess_dn_p:
                raise ValueError(
                    f"canonical block (k={k}, partner={partner}): "
                    f"n_excess_down_k={excess_dn_k} != n_excess_down_p={excess_dn_p} "
                    "pair-closure violated"
                )


def build_pair_list(
    stepped_occupation: np.ndarray,
    column_spin: np.ndarray,
    canonical_rows: list[int],
    partner_rows: np.ndarray,
    is_soc_mode: bool = False,
) -> list:
    """Emit the ordered pair list for General F construction.

    General path (``is_soc_mode=False``):
      Each pair is emitted from a canonical (k, partner(k)) block as a
      dict ``{"alpha": (k_row, col, spin_label),
               "beta":  (k_row, col, spin_label)}``.
      Non-self canonical k:
        1. (up@k, down@partner) - cross_kd pairs
        2. (up@partner, down@k) - cross_dk pairs
        3. (up@k, up@partner) - same-spin up excess
        4. (down@k, down@partner) - same-spin down excess
      Self-pair canonical k:
        1. (up@k, down@k) - cross
        2. (up@k, up@k) - same-spin up excess (2-at-a-time)
        3. (down@k, down@k) - same-spin down excess (2-at-a-time)

    SOC path (``is_soc_mode=True``):
      Column-index-based pairing, spin-agnostic (column_spin is not
      consulted). Each pair is emitted as a nested tuple
      ``((k_row_alpha, col_alpha), (k_row_beta, col_beta))``. Both
      ``k_row_*`` values are folded-BZ row indices into
      ``wavevector_index`` / ``eigenvector`` so the downstream
      ``build_slater_orbitals`` SOC branch can read the eigenvector at
      the correct k (the raw column index alone is ambiguous — the
      same column index at a different k refers to a different
      single-particle eigenstate).
      Non-self canonical k (k != partner(k)):
        Pair the alpha-th occupied column at k with the alpha-th
        occupied column at partner(k). Requires
        ``n_occ(k) == n_occ(partner)``; imbalance raises ``ValueError``.
      Self canonical k (k == partner(k)):
        Pair consecutive occupied columns (0, 1), (2, 3), .... Requires
        ``n_occ(k)`` even; an odd count raises ``ValueError``.

    See docs/en/source/algorithm/uhfk_to_mvmc.rst for canonical-pair
    selection and partner balance.
    """
    partner_rows = np.asarray(partner_rows, dtype=np.int64)

    if is_soc_mode:
        stepped_occupation = np.asarray(stepped_occupation, dtype=np.float64)
        nvol, ncol = stepped_occupation.shape
        pairs = []
        for k_row in canonical_rows:
            k = int(k_row)
            partner = int(partner_rows[k])
            occ_k = [
                c for c in range(ncol)
                if stepped_occupation[k, c] >= 0.5
            ]
            if partner == k:
                if len(occ_k) % 2 != 0:
                    raise ValueError(
                        f"SOC self canonical block k={k} has odd occupied "
                        f"count {len(occ_k)}; pair-closure impossible"
                    )
                for i in range(0, len(occ_k), 2):
                    pairs.append(((k, occ_k[i]), (k, occ_k[i + 1])))
            else:
                occ_p = [
                    c for c in range(ncol)
                    if stepped_occupation[partner, c] >= 0.5
                ]
                if len(occ_k) != len(occ_p):
                    raise ValueError(
                        f"SOC non-self canonical block (k={k}, "
                        f"partner={partner}) has imbalance "
                        f"n_occ_k={len(occ_k)} != n_occ_p={len(occ_p)}"
                    )
                for alpha in range(len(occ_k)):
                    pairs.append(((k, occ_k[alpha]), (partner, occ_p[alpha])))
        return pairs

    occ_up = _occ_by_spin(stepped_occupation, column_spin, 0)
    occ_dn = _occ_by_spin(stepped_occupation, column_spin, 1)
    pairs = []
    for k in canonical_rows:
        partner = int(partner_rows[k])
        if partner == k:
            u = occ_up[k]
            d = occ_dn[k]
            n_cross = min(len(u), len(d))
            for i in range(n_cross):
                pairs.append({
                    "alpha": (k, u[i], "up"),
                    "beta":  (k, d[i], "down"),
                })
            excess_u = u[n_cross:]
            excess_d = d[n_cross:]
            for i in range(len(excess_u) // 2):
                pairs.append({
                    "alpha": (k, excess_u[2 * i], "up"),
                    "beta":  (k, excess_u[2 * i + 1], "up"),
                })
            for i in range(len(excess_d) // 2):
                pairs.append({
                    "alpha": (k, excess_d[2 * i], "down"),
                    "beta":  (k, excess_d[2 * i + 1], "down"),
                })
        else:
            uk, dk = occ_up[k], occ_dn[k]
            up_, dp = occ_up[partner], occ_dn[partner]
            n_cross_kd = min(len(uk), len(dp))
            n_cross_dk = min(len(up_), len(dk))
            for i in range(n_cross_kd):
                pairs.append({
                    "alpha": (k, uk[i], "up"),
                    "beta":  (partner, dp[i], "down"),
                })
            for i in range(n_cross_dk):
                pairs.append({
                    "alpha": (partner, up_[i], "up"),
                    "beta":  (k, dk[i], "down"),
                })
            excess_uk = uk[n_cross_kd:]
            excess_up_ = up_[n_cross_dk:]
            excess_dk = dk[n_cross_dk:]
            excess_dp = dp[n_cross_kd:]
            for i in range(len(excess_uk)):
                pairs.append({
                    "alpha": (k, excess_uk[i], "up"),
                    "beta":  (partner, excess_up_[i], "up"),
                })
            for i in range(len(excess_dk)):
                pairs.append({
                    "alpha": (k, excess_dk[i], "down"),
                    "beta":  (partner, excess_dp[i], "down"),
                })
    return pairs


def _spin_row_offset(spin_label: str, Ns_phys: int) -> int:
    """Return 0 for 'up' rows (top half of 2Ns A matrix) or Ns_phys for
    'down' rows (bottom half). Used to place the folded orbital
    amplitude at the correct spin block."""
    if spin_label == "up":
        return 0
    if spin_label == "down":
        return Ns_phys
    raise ValueError(f"unknown spin label: {spin_label!r}")


def build_slater_orbitals(
    wavevector_index: np.ndarray,
    eigenvector: np.ndarray,
    column_spin: np.ndarray,
    site_positions: np.ndarray,
    cell_shape: np.ndarray,
    subshape: np.ndarray,
    theta: np.ndarray,
    pair_list: list[dict],
    is_soc_mode: bool = False,
) -> np.ndarray:
    """Extract physical Slater orbitals for every pair member into a
    (2*Ns_phys, 2*len(pair_list)) complex matrix A.

    Column layout: for pair index p in [0, len(pair_list)):
      - A[:, 2*p]     = ψ_alpha(i, spin_alpha)  (top half if alpha is up)
      - A[:, 2*p+1]   = ψ_beta (i, spin_beta )  (top half if beta is up)

    Amplitudes follow the documented positive-Bloch convention:
      ψ_alpha(i, sigma) = v[k_alpha, row(sigma, sub_offset(r_i)), col_alpha]
                       * exp(+i k_folded · folded_cell(r_i))
                       * exp(+i theta · r_i / L_phys)
                       / sqrt(nvol_folded)

    Non-zero only when sigma == column_spin[col_alpha] (mixed block is
    unsupported on the non-SOC path).

    ``is_soc_mode`` (default ``False``) is the SOC dispatch flag
    threaded through by the CLI. When True, A is built on the unified
    mVMC spin-block row index ``all_i = r_phys + spin * Ns_phys``
    (site-major, spin-minor), and ``pair_list`` entries are
    ``(alpha_1, alpha_2)`` column-index tuples instead of spin-labelled
    dicts. H-wave's interleaved eigenvector row
    ``2 * a_folded + spin`` is permuted to spin-block order at read
    time; every (r_phys, spin) row is populated per pair column.

    Under SOC + SubShape > [1, 1, 1], the plane-wave phase carries the
    intra-cell displacement via ``folded_cell + sub_offset``. The density
    gate (``compare_against_green_sublattice(is_soc_sublattice_mode=True)``)
    validates the shipping A directly by lifting ``green_sublattice`` to
    the physical basis via ``gauge_lift``, so no dual-A / reference-A
    scaffold is required. See docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    from .sublattice_unfold import (
        decode_physical_site,
        encode_folded_orbital,
        folded_row_indices,
    )

    wavevector_index = np.asarray(wavevector_index, dtype=np.int64)
    eigenvector = np.asarray(eigenvector, dtype=np.complex128)
    column_spin = np.asarray(column_spin, dtype=np.int64)
    site_positions = np.asarray(site_positions, dtype=np.int64)
    cell_shape = np.asarray(cell_shape, dtype=np.int64)
    subshape = np.asarray(subshape, dtype=np.int64)
    theta = np.asarray(theta, dtype=np.float64)

    Ns_phys = site_positions.shape[0]
    nvol_folded = wavevector_index.shape[0]
    subvol = int(np.prod(subshape))
    # This path requires norb_orig == 1, so norb_folded equals subvol.
    norb_folded = subvol
    L_folded = (cell_shape // subshape).astype(np.float64)
    L_phys = cell_shape.astype(np.float64)

    inv_sqrt = 1.0 / np.sqrt(float(nvol_folded))
    theta_over_L = theta / L_phys

    # Pre-compute per-site (folded_cell, sub_offset, folded_orb) for row
    # lookup. sub_offset is retained separately because the SOC branch
    # plane-wave uses ``folded_cell + sub_offset``; see
    # docs/en/source/algorithm/uhfk_to_mvmc.rst.
    folded_cell = np.empty((Ns_phys, 3), dtype=np.int64)
    sub_offset_per_site = np.empty((Ns_phys, 3), dtype=np.int64)
    folded_orb_per_site = np.empty(Ns_phys, dtype=np.int64)
    for i in range(Ns_phys):
        fc, so = decode_physical_site(site_positions[i], subshape)
        folded_cell[i] = fc
        sub_offset_per_site[i] = so
        folded_orb_per_site[i] = encode_folded_orbital(0, so, 1, subshape)

    # Pre-compute per-site physical gauge factor exp(+i theta r_i / L_phys).
    phys_arg = np.einsum(
        "d,id->i", theta_over_L, site_positions.astype(np.float64)
    )
    phys_up = np.exp(+1j * phys_arg)

    if is_soc_mode:
        # On the SOC branch, A rows are indexed in mVMC
        # site-major spin-block order all_i = r_phys + spin * Ns_phys.
        # H-wave's eigenvector rows use the interleaved SOC packing
        # 2 * a_folded + spin, so we permute at read time.
        # Under the single-orbital constraint (norb_orig == 1),
        # folded_orb_per_site[i] is a_folded for site i (sublattice
        # index within its folded cell); for SubShape=[1,1,1] it is 0.
        # pair_list entries are nested tuples
        # ``((k_row_alpha, col_alpha), (k_row_beta, col_beta))``; each
        # pair member identifies a specific single-particle Bloch state
        # at that k_row (the same column index at a different k_row
        # refers to a different eigenstate, so we must NOT sum over k
        # when reading the eigenvector).
        #
        # Plane-wave sign convention: H-wave's k-space
        # Hamiltonian is FFT'd with ``ifftn(..., norm='forward')``
        # (uhfk._make_ham_trans:1147), which uses the positive-Bloch
        # convention ``c_k = (1/sqrt(N)) Σ_R exp(+i k R) c_R``. Under
        # this convention the Bloch amplitude at real-space site R is
        # ``psi_alpha(R, s) = V[k, s, alpha] * exp(-i k R) / sqrt(N)``.
        # The non-SOC branch below still uses the ``+i k R`` sign
        # (its docstring calls this "positive-Bloch convention" — that
        # name refers to the sign of the plane_wave factor in the
        # bridge, not H-wave's Fourier convention). For non-SOC every
        # shipped fixture is inversion-symmetric so ``G(r) = G(-r)``
        # in each same-spin block and the sign discrepancy is
        # numerically undetectable in the density check. Under SOC the
        # off-diagonal spin blocks satisfy ``G(r)_{s,t} = -G(-r)_{s,t}``
        # (Rashba parity anti-symmetry), so a positive-Bloch bridge
        # combined with H-wave's own convention would flip the sign of
        # every s != t density entry and trip
        # ``compare_against_onebodyg_uhf_general``. The pair emission
        # rule pairs (k, partner(k)) with partner(k) = -k mod L,
        # so k_alpha + k_beta ≡ 0 in every pair and the F matrix is
        # invariant under this sign flip; the emitted F (and therefore
        # every mVMC observable) is unchanged.
        A_soc = np.zeros(
            (2 * Ns_phys, 2 * len(pair_list)), dtype=np.complex128
        )
        k_folded_all = (
            2.0 * np.pi
            * wavevector_index.astype(np.float64) / L_folded
        )
        # Under SOC + SubShape > [1, 1, 1],
        # ``SubShape > [1, 1, 1]`` H-wave's k-space Hamiltonian folds the
        # original lattice into a supercell where each folded cell holds
        # ``subvol`` sublattice slots. Reconstructing the physical Bloch
        # amplitude from ``V[k, s + 2*ir(R), alpha]`` requires the phase
        # ``exp(-i k_folded . (folded_cell(R) + sub_offset(R)))`` — not
        # just ``k_folded . folded_cell(R)``. The extra
        # ``k_folded . sub_offset(R)`` contribution is what encodes the
        # intra-supercell displacement into the Bloch factor; without it
        # the reconstructed A columns give the right density only when
        # ``ir(r_i) == ir(r_j)`` and produce wrong off-diagonal (in ir)
        # density entries. See
        # docs/en/source/algorithm/uhfk_to_mvmc.rst.
        # Under ``SubShape = [1, 1, 1]``, ``sub_offset`` is always zero
        # so the expression reduces to the SOC-only
        # (SubShape=[1,1,1]) formula.
        kf_dot_r = np.einsum(
            "kd,id->ki", k_folded_all,
            (folded_cell + sub_offset_per_site).astype(np.float64),
        )
        # The SOC + APBC branch reads folded
        # eigenvectors under the negative-Bloch convention
        # ``psi_alpha(R) = V[k, ...] * exp(-i k R) / sqrt(N)``, so the
        # twist gauge on the *same* branch must also be
        # ``exp(-i theta r / L_phys)`` for the k -> k + theta/L shift
        # to compose cleanly (each Bloch factor picks up the same twist
        # sign as its ``exp(-i k R)`` prefactor). Under PBC ``theta = 0``
        # this multiplication is trivially 1. Under APBC using
        # ``phys_up = exp(+i theta r / L)`` here would leave every pair
        # product ``plane_wave_a * plane_wave_b`` carrying an
        # ``exp(+i 2*theta*r/L)`` factor that does NOT cancel the folded
        # ``exp(-i (k_a + k_b) . fc)`` (whose sum lands on a non-zero
        # residue in the twisted mesh); the residual would appear as a
        # per-site ``(-1)^{r_x}`` alternation in F on-site diagonal
        # blocks and trip
        # ``aggregate_general_orbital_params`` class-consistency. The
        # non-SOC (positive-Bloch) branch retains
        # ``exp(+i theta r / L)`` because its plane-wave sign is ``+i k R``.
        phys_dn = np.exp(-1j * phys_arg)
        plane_wave_soc = (
            np.exp(-1j * kf_dot_r) * inv_sqrt * phys_dn[np.newaxis, :]
        )
        for p, pair in enumerate(pair_list):
            for column_idx, member in enumerate(pair):
                k_row_member, alpha = int(member[0]), int(member[1])
                slater_col = 2 * p + column_idx
                plane_wave_member = plane_wave_soc[k_row_member]
                for spin in (0, 1):
                    # H-wave interleaved eigenvector row per site:
                    # hwave_rows[i] = 2 * a_folded_i + spin
                    hwave_rows = (
                        2 * folded_orb_per_site + spin
                    ).astype(np.int64)
                    v_vals = np.array(
                        [
                            eigenvector[k_row_member, int(hwave_rows[i]), alpha]
                            for i in range(Ns_phys)
                        ],
                        dtype=np.complex128,
                    )
                    amp = plane_wave_member * v_vals
                    # Row assignment: all_i = r_phys + spin * Ns_phys
                    # (site-major, spin-minor); see
                    # docs/en/source/algorithm/uhfk_to_mvmc.rst.
                    row_start = spin * Ns_phys
                    A_soc[row_start:row_start + Ns_phys, slater_col] = amp
        return A_soc

    A = np.zeros((2 * Ns_phys, 2 * len(pair_list)), dtype=np.complex128)

    for p, pair in enumerate(pair_list):
        for offset, member in enumerate(("alpha", "beta")):
            k_row, col, spin_label = pair[member]
            # k_folded_phys per-direction = 2 pi n_tilde / L_folded
            k_folded = (
                2.0 * np.pi
                * wavevector_index[k_row].astype(np.float64) / L_folded
            )
            # exp(+i k_folded · folded_cell_i)
            kf_dot_fc = np.einsum(
                "d,id->i", k_folded, folded_cell.astype(np.float64)
            )
            plane_wave = np.exp(+1j * kf_dot_fc) * inv_sqrt * phys_up
            # Read pure-spin nd rows.
            row_up, row_dn = np.empty(Ns_phys, dtype=np.int64), np.empty(
                Ns_phys, dtype=np.int64
            )
            for i in range(Ns_phys):
                row_up[i], row_dn[i] = folded_row_indices(
                    int(folded_orb_per_site[i]), norb_folded
                )
            if spin_label == "up":
                v_row = row_up
            else:
                v_row = row_dn
            v_vals = np.array(
                [eigenvector[k_row, int(v_row[i]), col] for i in range(Ns_phys)],
                dtype=np.complex128,
            )
            spin_offset = _spin_row_offset(spin_label, Ns_phys)
            A[spin_offset:spin_offset + Ns_phys, 2 * p + offset] = (
                plane_wave * v_vals
            )

    return A


def build_fij_general(A: np.ndarray) -> np.ndarray:
    """Assemble the antisymmetric F matrix from the 2-column-per-pair
    Slater amplitude matrix A.

    F[iσ, jσ'] = sum_{p} [ A[iσ, 2p] * A[jσ', 2p+1]
                          - A[iσ, 2p+1] * A[jσ', 2p] ]

    F.shape = (A.shape[0], A.shape[0]) and F.T = -F elementwise. See
    docs/en/source/algorithm/uhfk_to_mvmc.rst.
    """
    A = np.asarray(A, dtype=np.complex128)
    n_row = A.shape[0]
    n_col = A.shape[1]
    assert n_col % 2 == 0, "A must have even column count (pairs of orbitals)"
    F = np.zeros((n_row, n_row), dtype=np.complex128)
    for p in range(n_col // 2):
        alpha = A[:, 2 * p]
        beta = A[:, 2 * p + 1]
        F += np.outer(alpha, beta) - np.outer(beta, alpha)
    return F
