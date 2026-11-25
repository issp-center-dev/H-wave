"""Sublattice folding of geometry and R-space interaction/transfer tables.

Single definition shared by UHFk and RPA/FLEX: the two solvers previously
carried near-verbatim copies of these folds, and bug fixes (physical-orbital
stride, wrap-around canonicalization/accumulation) had to be mirrored
between them by hand.

``fold_layout`` / ``unfold_layout`` are the matching layout change of a
dense R-space array (a Green function or a Hamiltonian), in either slot
convention; ``UHFk._reshape_green`` / ``UHFk._deflate_green`` call them.
"""
import itertools
import logging

import numpy as np

logger = logging.getLogger("qlms").getChild("fold")


def reshape_geometry(geom, subshape):
    """Fold a geometry definition onto the sublattice (supercell) basis.

    Parameters
    ----------
    geom : dict
        Original geometry with 'norb', 'rvec' and 'center'.
    subshape : tuple of int
        Supercell shape (Bx, By, Bz).

    Returns
    -------
    dict
        The folded geometry (norb multiplied by the supercell volume,
        rvec scaled, orbital centers replicated per supercell cell).
    """
    Bx, By, Bz = subshape
    bvol = Bx * By * Bz

    norb = geom['norb']

    geom_new = {}
    geom_new['norb'] = geom['norb'] * bvol
    geom_new['rvec'] = np.matmul(np.diag([Bx, By, Bz]), geom['rvec'])

    sc = np.array([1.0/Bx, 1.0/By, 1.0/Bz])
    cw = [ sc * geom['center'][k] for k in range(norb) ]

    centerv = np.zeros((norb * bvol, 3), dtype=np.float64)
    k = 0
    for bz, by, bx in itertools.product(range(Bz), range(By), range(Bx)):
        for i in range(norb):
            centerv[k] = cw[i] + np.array([bx, by, bz]) * sc
            k += 1
    geom_new['center'] = centerv

    return geom_new


def reshape_interaction(ham, subshape, shape, norb_so_orig, norb_phys_orig,
                        enable_spin_orbital, *, drop_spin_block=False):
    """Fold an R-space table onto the sublattice (supercell) basis.

    Parameters
    ----------
    ham : dict
        {((rx,ry,rz), (alpha,beta)): value} in the original cell basis.
    subshape : tuple of int
        Supercell shape (Bx, By, Bz).
    shape : tuple of int
        Folded lattice shape in units of the supercell (Nx, Ny, Nz).
    norb_so_orig : int
        Geometry norb of the original cell (the spin-orbital count in
        enable_spin_orbital mode).
    norb_phys_orig : int
        Physical orbital count of the original cell (equals norb_so_orig
        when spin-orbital mode is off).
    drop_spin_block : bool, keyword-only
        True for ONE-BODY tables (Transfer, Extern) read in normal mode from
        a spin-orbital-format file: entries with an orbital index >=
        norb_phys_orig are the spin block, which the solvers ignore with a
        warning, so the fold drops them (once, with the same warning).
        False (default) for two-body tables, whose indices are physical in
        every mode: such an entry there is invalid input and raises
        ValueError instead of being dropped.
    enable_spin_orbital : bool
        True when THIS table uses spin-orbital indices (the Transfer in SO
        mode).  Two-body interaction tables always use physical orbital
        indices and pass False.

    Returns
    -------
    dict
        The folded table, with R canonicalized to [0, n) and wrap-around
        collisions accumulated.
    """
    Bx, By, Bz = subshape
    nx, ny, nz = shape

    def _reshape_orbit_(a, x):
        # Two-body interaction files use PHYSICAL orbital indices even in
        # spin-orbital mode, and the interaction tables are laid out as
        # a + norb_phys_orig * cellidx -- so fold with the physical stride.
        return a + norb_phys_orig * (x[0] + Bx * (x[1] + By * (x[2])))

    def _reshape_orbit_spin(a, x):
        a_, s_ = a % norb_so_orig, a // norb_so_orig
        return a_ + norb_so_orig * (x[0] + Bx * (x[1] + By * (x[2] + Bz * s_)))

    if enable_spin_orbital:
        _reshape_orbit = _reshape_orbit_spin
    else:
        _reshape_orbit = _reshape_orbit_

    def _round(x, n):
        # canonicalize the folded supercell index to [0, n): the consumers
        # scatter with tab_r[ir] += v, where a negative index wraps onto
        # the same slot as its positive counterpart
        return x % n

    ham_new = {}
    has_spin_dep = False
    for (irvec, orbvec), v in ham.items():
        rx, ry, rz = irvec
        alpha, beta = orbvec

        if not enable_spin_orbital and (alpha >= norb_phys_orig
                                        or beta >= norb_phys_orig):
            if not drop_spin_block:
                # two-body tables are physical-indexed in every mode, so
                # this is invalid input, not an ignorable spin term; fail
                # closed rather than fold it (the unfolded consumers
                # raise on it too)
                raise ValueError(
                    "interaction table has orbital index {} >= norb ({}) "
                    "at irvec={}: two-body interaction tables use physical "
                    "orbital indices in every mode, so this entry is "
                    "invalid input (it is not a spin-dependent term that "
                    "could be ignored).".format(
                        tuple(orbvec), norb_phys_orig, tuple(irvec)))
            # Spin-block entries of a spin-orbital-format one-body table
            # read with enable_spin_orbital=False: the solvers drop such
            # terms from H(k), and folding them with the physical stride
            # would relabel them onto genuine sublattice orbitals.
            has_spin_dep = True
            continue

        for bz, by, bx in itertools.product(range(Bz), range(By), range(Bx)):

            # original cell index of both ends
            #   x0 -> x1=x0+r
            x0, y0, z0 = bx, by, bz
            x1, y1, z1 = x0 + rx, y0 + ry, z0 + rz

            # decompose into supercell-index and cell-index within supercell
            #   x0 = 0 + x0
            #   x1 = X + xr
            xx1, xr1 = x1 // Bx, x1 % Bx
            yy1, yr1 = y1 // By, y1 % By
            zz1, zr1 = z1 // Bz, z1 % Bz

            # find orbital index within supercell
            aa = _reshape_orbit(alpha, (x0, y0, z0))
            bb = _reshape_orbit(beta, (xr1, yr1, zr1))

            # Wrap-around can map distinct original R-vectors onto the same
            # folded (ir, ov) key (e.g. when the hopping range reaches the
            # folded supercell size); accumulate so such contributions are
            # summed rather than overwritten.
            ir = (_round(xx1, nx), _round(yy1, ny), _round(zz1, nz))
            ov = (aa, bb)

            ham_new[(ir, ov)] = ham_new.get((ir, ov), 0.0) + v

    if has_spin_dep:
        # the dropped entries no longer reach the solvers' own spin-block
        # warning, so report the ignore here instead
        logger.warning(
            "table has orbital indices >= norb (spin-dependent terms) "
            "but enable_spin_orbital is False. "
            "These terms are ignored. Set enable_spin_orbital = true to use them."
        )

    return ham_new


def _layout_dims(cellshape, subshape):
    """Validate the cell/sub shapes; return (cellshape, subshape, shape)."""
    cellshape = tuple(int(n) for n in cellshape)
    subshape = tuple(int(n) for n in subshape)
    if len(cellshape) != 3 or len(subshape) != 3:
        raise ValueError(
            "cellshape and subshape must have three entries, got {} and {}"
            .format(cellshape, subshape))
    if any(b < 1 for b in subshape) or any(
            L % b != 0 for L, b in zip(cellshape, subshape)):
        raise ValueError(
            "subshape {} must divide cellshape {} on every axis".format(
                subshape, cellshape))
    shape = tuple(L // b for L, b in zip(cellshape, subshape))
    return cellshape, subshape, shape


def _within_cell_offsets(subshape, norb_orig):
    """Per folded orbital ``aa = a + norb_orig*ir``: ``a`` and the in-cell
    coordinates of ``ir = ix + Bx*(iy + By*iz)`` (x fastest)."""
    Bx, By, Bz = subshape
    orb = np.arange(norb_orig * Bx * By * Bz)
    a = orb % norb_orig
    ir = orb // norb_orig
    return a, ir % Bx, (ir // Bx) % By, (ir // (Bx * By)) % Bz


def fold_layout(green, *, cellshape, subshape, norb_orig, ns,
                hamiltonian=False):
    """Fold an R-space array from the original lattice onto the sublattice.

    ``(Lvol, ns, norb_orig, ns, norb_orig)`` -> ``(Nvol, ns, norb, ns, norb)``
    with ``norb = norb_orig * Bvol``. Site indices are C order (z fastest):
    ``r = rz + Nz*(ry + Ny*rx)``, on the original lattice ``(Lx, Ly, Lz)``
    and on the supercell lattice ``(Nx, Ny, Nz) = cellshape // subshape``
    alike. A folded orbital is ``aa = a + norb_orig*ir`` with the in-cell
    index ``ir = ix + Bx*(iy + By*iz)`` (x fastest).

    ``green_sub[R, s, aa, t, bb] = green[r, s, a, t, b]`` with
    ``r = R*B + dr`` (per axis, modulo the cell shape) and

    * ``hamiltonian=False`` (Green convention): ``dr = ri - rj``, the inverse
      of :func:`unfold_layout` with the offset on the first slot; this is
      ``UHFk._reshape_green``;
    * ``hamiltonian=True`` (Hamiltonian / transfer convention):
      ``dr = rj - ri``, the inverse of ``unfold_layout(hamiltonian=True)``
      and the same layout ``reshape_interaction`` gives a transfer table.

    The two conventions differ because a folded Green function is the
    orbital transpose of the Hamiltonian-convention object (see
    ``UHFk._deflate_green``). The result is a new C-contiguous complex128
    array; the copy is exact.

    Raises
    ------
    ValueError
        If ``green`` does not have the shape above (e.g. an already folded
        array), or ``subshape`` does not divide ``cellshape``.
    """
    cellshape, subshape, shape = _layout_dims(cellshape, subshape)
    Lx, Ly, Lz = cellshape
    Bx, By, Bz = subshape
    Nx, Ny, Nz = shape
    norb_orig, ns = int(norb_orig), int(ns)
    Lvol, Nvol = Lx * Ly * Lz, Nx * Ny * Nz

    green = np.asarray(green)
    expected = (Lvol, ns, norb_orig, ns, norb_orig)
    if green.shape != expected:
        raise ValueError(
            "fold_layout expects the original-lattice shape (Lvol, ns, "
            "norb_orig, ns, norb_orig) = {} for cellshape {}, ns = {}, "
            "norb_orig = {}; got {}".format(expected, cellshape, ns,
                                            norb_orig, green.shape))

    isite = np.arange(Nvol)
    ix0 = ((isite // (Nz * Ny)) % Nx) * Bx
    iy0 = ((isite // Nz) % Ny) * By
    iz0 = (isite % Nz) * Bz

    a, rx, ry, rz = _within_cell_offsets(subshape, norb_orig)
    if hamiltonian:
        dx = rx[np.newaxis, :] - rx[:, np.newaxis]
        dy = ry[np.newaxis, :] - ry[:, np.newaxis]
        dz = rz[np.newaxis, :] - rz[:, np.newaxis]
    else:
        dx = rx[:, np.newaxis] - rx[np.newaxis, :]
        dy = ry[:, np.newaxis] - ry[np.newaxis, :]
        dz = rz[:, np.newaxis] - rz[np.newaxis, :]

    jx = (ix0[:, None, None] + dx[None, :, :]) % Lx
    jy = (iy0[:, None, None] + dy[None, :, :]) % Ly
    jz = (iz0[:, None, None] + dz[None, :, :]) % Lz
    jsite = jz + Lz * (jy + Ly * jx)                      # (Nvol, norb, norb)

    spin = np.arange(ns)
    # out[R, s, aa, t, bb] = green[jsite[R, aa, bb], s, a[aa], t, a[bb]]
    out = green[jsite[:, None, :, None, :],
                spin[None, :, None, None, None],
                a[None, None, :, None, None],
                spin[None, None, None, :, None],
                a[None, None, None, None, :]]
    return np.ascontiguousarray(out, dtype=np.complex128)


def unfold_layout(green_sub, *, cellshape, subshape, norb_orig, ns,
                  hamiltonian=False):
    """Unfold a sublattice-basis R-space array back to the original lattice.

    ``(Nvol, ns, norb, ns, norb)`` -> ``(Lvol, ns, norb_orig, ns,
    norb_orig)``, the layouts of :func:`fold_layout`. The original site
    ``r`` decomposes per axis into the supercell site ``R`` and the in-cell
    offset ``ir``, and ``green[r, s, a, t, b] = green_sub[R, s, aa, t, bb]``
    with

    * ``hamiltonian=False`` (Green convention): ``aa = a + norb_orig*ir``,
      ``bb = b``;
    * ``hamiltonian=True`` (Hamiltonian / transfer convention): ``aa = a``,
      ``bb = b + norb_orig*ir``.

    This is ``UHFk._deflate_green``. ``unfold_layout(fold_layout(g))`` is
    ``g`` exactly, for either convention (the same flag on both calls). The
    result is a new C-contiguous complex128 array.

    Raises
    ------
    ValueError
        If ``green_sub`` does not have the shape above, or ``subshape`` does
        not divide ``cellshape``.
    """
    cellshape, subshape, shape = _layout_dims(cellshape, subshape)
    Lx, Ly, Lz = cellshape
    Bx, By, Bz = subshape
    Nx, Ny, Nz = shape
    norb_orig, ns = int(norb_orig), int(ns)
    Lvol = Lx * Ly * Lz
    Nvol = Nx * Ny * Nz
    norb = norb_orig * Bx * By * Bz

    green_sub = np.asarray(green_sub)
    expected = (Nvol, ns, norb, ns, norb)
    if green_sub.shape != expected:
        raise ValueError(
            "unfold_layout expects the sublattice shape (Nvol, ns, "
            "norb_orig*Bvol, ns, norb_orig*Bvol) = {} for cellshape {}, "
            "subshape {}, ns = {}, norb_orig = {}; got {}".format(
                expected, cellshape, subshape, ns, norb_orig,
                green_sub.shape))

    jsite = np.arange(Lvol)
    ix = (jsite // (Lz * Ly)) % Lx
    iy = (jsite // Lz) % Ly
    iz = jsite % Lz
    isite = (iz // Bz) + Nz * ((iy // By) + Ny * (ix // Bx))
    ir = (ix % Bx) + Bx * ((iy % By) + By * (iz % Bz))

    orb = np.arange(norb_orig)
    shifted = orb[None, :] + norb_orig * ir[:, None]      # (Lvol, norb_orig)
    plain = np.broadcast_to(orb[None, :], shifted.shape)
    if hamiltonian:
        aa, bb = plain, shifted
    else:
        aa, bb = shifted, plain

    spin = np.arange(ns)
    # out[r, s, a, t, b] = green_sub[isite[r], s, aa[r, a], t, bb[r, b]]
    out = green_sub[isite[:, None, None, None, None],
                    spin[None, :, None, None, None],
                    aa[:, None, :, None, None],
                    spin[None, None, None, :, None],
                    bb[:, None, None, None, :]]
    return np.ascontiguousarray(out, dtype=np.complex128)
