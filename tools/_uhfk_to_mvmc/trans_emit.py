"""H-wave Transfer.dat -> mVMC trans.def emitter for SOC.

Reads H-wave's Wannier90-like ``Transfer.dat`` (SOI-packed ``iWan/jWan``
indices ``2 * a_phys + spin + 1``) and emits mVMC's real-space
``trans.def`` with ``(i, s, j, t)`` rows for the full physical lattice,
including the spin-off-diagonal entries that ``vmcdry.out``'s
``FermionHubbardGC`` generator drops silently.

General complex-SOC mapping
---------------------------
ComplexUHF constructs its bare Hamiltonian as ``K = -trans``. In H-wave's
negative-Bloch physical basis, a parsed entry ``(R, s, t, v)`` maps at target
displacement ``+R`` as

``K[i, t; i+R, s] = conj(v)`` and
``trans[i, t; i+R, s] = -conj(v)``.

Thus the site endpoints remain ``i -> i+R``, while the spin endpoints are
swapped and the coefficient is conjugated and negated. On a real
spin-diagonal entry the spin swap is a no-op and the rule reduces to
``trans = -v``.

See docs/en/source/algorithm/uhfk_to_mvmc.rst for the mapping derivation,
Hermiticity convention, and boundary-wrap composition.
"""
from __future__ import annotations

import numpy as np

from hwave.solver._apbc_phase import inverse_gauge_phase


class TransEmitError(ValueError):
    """Raised on malformed Transfer.dat or trans.def emission failure."""


def _unpack_soi(iWan: int) -> tuple[int, int]:
    """Inverse of ``emit_rashba_transfer._pack_index``: ``iWan = 2 * a + s + 1``.

    Returns ``(a_phys, spin)``. With ``norb_orig = 1``,
    ``a_phys`` is always ``0`` and ``spin`` is ``{0, 1}``.
    """
    if iWan < 1:
        raise TransEmitError(
            f"invalid 1-based orbital index iWan={iWan} (must be >= 1)"
        )
    idx0 = iWan - 1
    return (idx0 // 2, idx0 % 2)


def parse_hwave_transfer(path):
    """Parse H-wave ``Transfer.dat`` (Wannier90-like format).

    File layout (see ``src/hwave/qlmsio/wan90.py::read_w90``)::

        <header line>                               # skipped
        <num_wann>                                  # skipped (unused)
        <nr>                                        # number of entries
        <ndegen block: nr ints, 15 per line>        # applied as divisor
        <data rows: rx ry rz iWan jWan re im>

    Applies the R-point degeneracy division from the Wannier90 convention
    ``H(k) = sum_R exp(i k.R) H(R) / ndegen(R)``. Mirrors the semantics of
    ``hwave.qlmsio.wan90.read_w90``: if every ndegen entry is ``1`` the
    division is a no-op (sparse H-wave format); otherwise each entry's
    coefficient is divided by ``ndegen[i]``, where ``i`` is the position of
    the entry's R-vector in file order (dense listing required for the
    non-unit case).

    Parameters
    ----------
    path : str
        Path to ``Transfer.dat``.

    Returns
    -------
    entries : list of tuples
        Each entry is ``(rx, ry, rz, iWan, jWan, val_complex)``. ``rx``,
        ``ry``, ``rz`` are signed real-space displacement components,
        ``iWan``/``jWan`` are 1-based SOI-packed orbital indices, and
        ``val_complex`` is ``(re + 1j * im) / ndegen``.

    Raises
    ------
    TransEmitError
        On file layout errors, truncated files, ndegen count mismatch, or
        a non-unit ndegen block paired with a sparse file listing (fewer
        distinct R-vectors than ``nr``).
    """
    try:
        with open(path) as fp:
            lines = fp.read().splitlines()
    except OSError as e:
        raise TransEmitError(f"cannot open Transfer.dat: {path}: {e}") from e

    if len(lines) < 3:
        raise TransEmitError(
            f"{path}: too short for Wannier90-like header "
            f"(got {len(lines)} lines, need >= 3)"
        )

    # Skip header (line 0), num_wann (line 1). Read nr (line 2).
    try:
        nr = int(lines[2].strip())
    except ValueError as e:
        raise TransEmitError(
            f"{path}: line 3 is not an integer nr (got {lines[2]!r})"
        ) from e
    if nr < 0:
        raise TransEmitError(f"{path}: negative nr = {nr}")

    # Parse ndegen block: nr integers, 15 per line -> ceil(nr / 15) lines.
    ndegen_rows = (nr + 14) // 15
    data_start = 3 + ndegen_rows
    if len(lines) < data_start:
        raise TransEmitError(
            f"{path}: truncated ndegen block (need {ndegen_rows} rows for "
            f"nr={nr}, only {len(lines) - 3} available)"
        )
    ndegen = []
    for i, line in enumerate(lines[3:data_start], start=4):
        try:
            ndegen.extend(int(x) for x in line.split())
        except ValueError as e:
            raise TransEmitError(
                f"{path}:{i}: cannot parse ndegen row: {line!r}"
            ) from e
    if len(ndegen) != nr:
        raise TransEmitError(
            f"{path}: ndegen count mismatch (declared nr={nr}, found "
            f"{len(ndegen)})"
        )
    # When every degeneracy is one, division is a no-op and the sparse
    # H-wave listing is legal (nr may be a placeholder). Only non-unit
    # ndegen requires the positional R -> ndegen[i] mapping.
    all_unit = all(d == 1 for d in ndegen)
    deg_of_r = {}  # irvec -> ndegen, keyed by first appearance in file.
    seen = set()  # (rx, ry, rz, iWan, jWan) tuples already emitted.

    entries = []
    for line_num, line in enumerate(lines[data_start:], start=data_start + 1):
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#") or stripped.startswith("!"):
            continue
        toks = stripped.split()
        if len(toks) < 7:
            raise TransEmitError(
                f"{path}:{line_num}: expected 7 columns "
                f"(rx ry rz iWan jWan re im), got {len(toks)}: {stripped!r}"
            )
        try:
            rx, ry, rz = int(toks[0]), int(toks[1]), int(toks[2])
            iWan, jWan = int(toks[3]), int(toks[4])
            re, im = float(toks[5]), float(toks[6])
        except ValueError as e:
            raise TransEmitError(
                f"{path}:{line_num}: cannot parse row: {stripped!r}"
            ) from e
        # Mirror hwave.qlmsio.wan90.read_w90's duplicate-entry rejection.
        # Two rows with the same (R, iWan, jWan) key would be silently
        # summed into trans.def by emit_trans_def below, double-counting
        # the hopping. H-wave's reader treats this as fatal; we match.
        key = (rx, ry, rz, iWan, jWan)
        if key in seen:
            raise TransEmitError(
                f"{path}:{line_num}: duplicate Transfer.dat entry "
                f"(rx, ry, rz, iWan, jWan) = {key}; H-wave's read_w90 "
                "rejects duplicate hopping entries -- matching that "
                "contract here."
            )
        seen.add(key)
        # Mirror hwave.qlmsio.wan90.read_w90 ndegen lookup.
        if all_unit:
            deg = 1
        else:
            irvec = (rx, ry, rz)
            if irvec not in deg_of_r:
                idx_r = len(deg_of_r)
                if idx_r >= len(ndegen):
                    raise TransEmitError(
                        f"{path}:{line_num}: more distinct R-points than "
                        f"declared (nr={nr})"
                    )
                deg_of_r[irvec] = ndegen[idx_r]
            deg = deg_of_r[irvec]
        entries.append((rx, ry, rz, iWan, jWan, complex(re, im) / deg))

    # Reject non-unit + sparse listing (fewer distinct Rs than declared),
    # matching hwave.qlmsio.wan90.read_w90.
    if not all_unit and len(deg_of_r) != nr:
        raise TransEmitError(
            f"{path}: non-unit ndegen requires a dense file listing all "
            f"{nr} R-points (found {len(deg_of_r)})"
        )
    return entries


def emit_trans_def(
    transfer_path, cell_shape, out_path,
    boundary_theta=None,
):
    """Emit mVMC ``trans.def`` from H-wave ``Transfer.dat``.

    For each ``(R = (rx, ry, rz), iWan, jWan, val)`` entry in
    ``Transfer.dat`` and each source site ``i_src = (ix, iy, iz)`` in
    the physical lattice, emits one row::

        i_src_flat  s_tgt  j_tgt_flat  s_src  value.real  value.imag

    where ``j_tgt = ((ix + rx) mod Lx, (iy + ry) mod Ly, (iz + rz) mod Lz)``
    (PBC unfold), ``s_src, s_tgt = unpack_soi(iWan, jWan)``, and
    ``value = -conj(val) * P`` for the mVMC-frame boundary wrap phase ``P``.
    The flat index is site-major (``i = ix + Lx * (iy + Ly * iz)``, matching
    H-wave's ``Geometry.dat`` layout). The site endpoints are unchanged while
    the spin endpoints are swapped.

    Parameters
    ----------
    transfer_path : str
        H-wave ``Transfer.dat`` path.
    cell_shape : list-like of length 3
        Physical lattice shape ``[Lx, Ly, Lz]``.
    out_path : str
        Output ``trans.def`` path.
    boundary_theta : array-like of length 3 or None, optional
        Twist ``(theta_x, theta_y, theta_z)`` in radians. ``None``
        (default) means all-PBC and no boundary phase is applied.
        Non-zero components mean twisted / antiperiodic BC in the
        corresponding direction. For a bond whose displacement takes
        the target site out of the primary cell in direction ``d``,
        the emitted row acquires the physical wrap phase
        ``exp(i * theta_d * wraps_d)`` (``= (-1)^wraps_d`` for
        ``theta_d = pi`` / APBC). Non-wrapping bonds are unchanged.
        This mirrors the sign flip StdFace bakes into vmcdry.out's
        ``trans.def`` under ``phase0 = 180`` so mVMC's ``H = -sum trans``
        recovers the physical Hamiltonian on the periodic-site frame.

    Raises
    ------
    TransEmitError
        On malformed input or unsupported SOI unpacking (e.g. a physical
        orbital index != 0 under the single-orbital constraint).

    See docs/en/source/algorithm/uhfk_to_mvmc.rst for the Transfer mapping.
    """
    entries = parse_hwave_transfer(transfer_path)
    emit_trans_def_from_entries(
        entries,
        cell_shape,
        out_path,
        boundary_theta=boundary_theta,
    )


def emit_trans_def_from_entries(
    entries, cell_shape, out_path,
    boundary_theta=None,
):
    """Emit mVMC ``trans.def`` from already parsed Transfer.dat entries.

    The entry layout and emission mapping are identical to
    :func:`emit_trans_def`; this entry point lets callers validate the parsed
    values and emit those same values without parsing the source twice.
    """
    entries = tuple(entries)
    validate_trans_def_entries(entries, cell_shape)
    Lx, Ly, Lz = (int(c) for c in cell_shape)

    # For SOC + APBC, compose the per-row wrap phase for
    # boundary-crossing bonds. ``inverse_gauge_phase(r_j_wrapped,
    # r_j_unwrapped, theta, L)`` evaluates to
    # ``exp(i * theta . (r_j_wrapped - r_j_unwrapped) / L)`` which
    # equals ``(-1)^wraps_d`` in each AP direction (theta_d = pi and
    # ``r_j_unwrapped - r_j_wrapped`` a multiple of L_d in that
    # direction, with the wrap count as the multiplier). Non-wrapping
    # bonds land on ``exp(0) = 1`` and are unchanged. See
    # docs/en/source/algorithm/uhfk_to_mvmc.rst.
    if boundary_theta is not None and np.any(
        np.abs(np.asarray(boundary_theta, dtype=np.float64)) > 1e-12
    ):
        theta_arr = np.asarray(boundary_theta, dtype=np.float64)
        L_arr = np.asarray(cell_shape, dtype=np.float64)
        apply_gauge = True
    else:
        theta_arr = None
        L_arr = None
        apply_gauge = False

    def site_index(ix, iy, iz):
        return ix + Lx * (iy + Ly * iz)

    rows = []
    for rx, ry, rz, iWan, jWan, val in entries:
        _, s_src = _unpack_soi(iWan)
        _, s_tgt = _unpack_soi(jWan)
        for iz in range(Lz):
            for iy in range(Ly):
                for ix in range(Lx):
                    i_site = site_index(ix, iy, iz)
                    jx = (ix + rx) % Lx
                    jy = (iy + ry) % Ly
                    jz = (iz + rz) % Lz
                    j_site = site_index(jx, jy, jz)
                    # Conjugate the H-wave coefficient first, then multiply
                    # by the unconjugated mVMC-frame wrap phase P:
                    # value = -conj(val) * P. Every shipped fixture has
                    # theta in {0, pi}, so P is real (+/-1) and conj(P) == P;
                    # those fixtures cannot distinguish this ordering from
                    # -conj(val * P). No fixture verifies a non-real twist.
                    v_val = -np.conjugate(val)
                    if apply_gauge:
                        r_j_wrapped = np.array(
                            [jx, jy, jz], dtype=np.float64
                        )
                        r_j_unwrapped = np.array(
                            [ix + rx, iy + ry, iz + rz], dtype=np.float64
                        )
                        wrap_phase = inverse_gauge_phase(
                            r_j_wrapped, r_j_unwrapped, theta_arr, L_arr
                        )
                        v_val = v_val * wrap_phase
                    rows.append(
                        (
                            i_site, s_tgt, j_site, s_src,
                            v_val.real, v_val.imag,
                        )
                    )

    # Mirror vmcdry.out's trans.def header layout so mVMC's parser
    # accepts the file interchangeably with the vmcdry-generated one.
    with open(out_path, "w") as fw:
        fw.write("======================== \n")
        fw.write("NTransfer      {}  \n".format(len(rows)))
        fw.write("======================== \n")
        fw.write("========i_j_s_tijs====== \n")
        fw.write("======================== \n")
        for (i, s, j, t, re, im) in rows:
            fw.write(
                "{:5d}{:6d}{:6d}{:6d}{:26.15f}{:26.15f}\n".format(
                    i, s, j, t, re, im
                )
            )


def validate_trans_def_entries(entries, cell_shape):
    """Validate constraints required to emit parsed Transfer.dat entries."""
    if len(cell_shape) != 3:
        raise TransEmitError(
            f"cell_shape must have length 3, got {list(cell_shape)}"
        )
    Lx, Ly, Lz = (int(c) for c in cell_shape)
    if Lx <= 0 or Ly <= 0 or Lz <= 0:
        raise TransEmitError(
            f"cell_shape components must be positive: {[Lx, Ly, Lz]}"
        )

    for _rx, _ry, _rz, iWan, jWan, _val in entries:
        a_src, _s_src = _unpack_soi(iWan)
        a_tgt, _s_tgt = _unpack_soi(jWan)
        # Require one physical orbital per site.
        if a_src != 0 or a_tgt != 0:
            raise TransEmitError(
                f"trans_emit requires norb_orig == 1; got "
                f"iWan={iWan} -> a_src={a_src}, jWan={jWan} -> "
                f"a_tgt={a_tgt}"
            )
