"""Sublattice folding of geometry and R-space interaction/transfer tables.

Single definition shared by UHFk and RPA/FLEX: the two solvers previously
carried near-verbatim copies of these folds, and bug fixes (physical-orbital
stride, wrap-around canonicalization/accumulation) had to be mirrored
between them by hand.
"""
import itertools

import numpy as np


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
                        enable_spin_orbital):
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
    for (irvec, orbvec), v in ham.items():
        rx, ry, rz = irvec
        alpha, beta = orbvec

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

    return ham_new
