"""Frozen verbatim copies of UHFk methods at revision d72f35f6.

These are oracles for the bitwise-equality tests of code lifted out of
``UHFk``: each function body is the d72f35f6 method body, unchanged, taking
the solver instance (or any object with the same attributes) as ``self``.
Do not edit them; compare new code against them in the same process.
"""
import itertools

import numpy as np


def reshape_green(self, green):
    """Convert Green function to sublattice basis.

    Parameters
    ----------
    green : ndarray
        Green function in original basis

    Returns
    -------
    ndarray
        Green function in sublattice basis
    """
    # convert green function into sublattice

    Lx,Ly,Lz = self.cellshape
    Lvol = Lx * Ly * Lz
    Bx,By,Bz = self.subshape
    Bvol = Bx * By * Bz
    Nx,Ny,Nz = self.shape
    Nvol = Nx * Ny * Nz

    norb_orig = self.norb_orig
    norb = self.norb
    ns = self.ns

    # check array size
    assert(green.shape == (Lvol,ns,norb_orig,ns,norb_orig))

    def _pack_index(x, n):
        _ix, _iy, _iz = x
        _nx, _ny, _nz = n
        return _ix + _nx * (_iy + _ny * (_iz))

    def _unpack_index(x, n):
        _nx, _ny, _nz = n
        _ix = x % _nx
        _iy = (x // _nx) % _ny
        _iz = (x // (_nx * _ny)) % _nz
        return (_ix, _iy, _iz)

    def _pack_site(x, n):
        _ix, _iy, _iz = x
        _nx, _ny, _nz = n
        return _iz + _nz * (_iy + _ny * (_ix))

    def _unpack_site(idx, n):
        _nx, _ny, _nz = n
        _iz = idx % _nz
        _iy = (idx // _nz) % _ny
        _ix = (idx // (_nz * _ny)) % _nx
        return (_ix,_iy,_iz)

    green_sub = np.zeros((Nvol,ns,norb,ns,norb), dtype=np.complex128)

    for isite in range(Nvol):
        ixx, iyy, izz = _unpack_site(isite, (Nx,Ny,Nz))
        ix0, iy0, iz0 = ixx * Bx, iyy * By, izz * Bz

        for aa, bb in itertools.product(range(norb), range(norb)):
            a, ri = aa % norb_orig, aa // norb_orig
            b, rj = bb % norb_orig, bb // norb_orig

            rix, riy, riz = _unpack_index(ri, (Bx,By,Bz))
            rjx, rjy, rjz = _unpack_index(rj, (Bx,By,Bz))

            # Green convention: the within-cell offset sits on the first
            # orbital slot (aa), matching _deflate_green's default. Folding
            # uses dr = ri - rj, the inverse of that deflate, so folding an
            # unfolded green and deflating it round-trips. (Hamiltonian/
            # transfer folding uses dr = rj - ri and lives in
            # _reshape_interaction.)
            ix = (ix0 + rix - rjx) % Lx
            iy = (iy0 + riy - rjy) % Ly
            iz = (iz0 + riz - rjz) % Lz

            jsite = _pack_site((ix,iy,iz), (Lx,Ly,Lz))

            for s, t in itertools.product(range(ns), range(ns)):
                green_sub[isite, s, aa, t, bb] = green[jsite, s, a, t, b]

    return green_sub


def deflate_green(self, green_sub, hamiltonian=False):
    """Convert Green function back to original basis.

    Parameters
    ----------
    green_sub : ndarray
        Green function in sublattice basis
    hamiltonian : bool, optional
        Select the cross-sublattice slot convention. With
        ``green[jsite,s,a,t,b] = green_sub[isite,s,aa,t,bb]``, where the
        original site ``jsite`` decomposes per axis into a supercell site
        ``isite`` and the within-cell offset ``ir``, the default (Green)
        carries ``ir`` on the FIRST orbital slot
        (``aa = a + norb_orig*ir``, ``bb = b``); ``hamiltonian=True`` carries
        it on the SECOND slot (``aa = a``, ``bb = b + norb_orig*ir``). The
        two differ because the folded Green ``self.Green`` from ``_green()``
        is the orbital transpose of the Hamiltonian-convention object (it
        stores ``G_ab(k) = sum_l conj(V_al) f_l V_bl``). Use the default for
        a Green function and ``hamiltonian=True`` for Hamiltonian / transfer
        quantities (e.g. ``trans_mod``).

    Returns
    -------
    ndarray
        Green function in original basis
    """
    # convert green function back to original lattice
    Lx,Ly,Lz = self.cellshape
    Lvol = Lx * Ly * Lz
    Bx,By,Bz = self.subshape
    Bvol = Bx * By * Bz
    Nx,Ny,Nz = self.shape
    Nvol = Nx * Ny * Nz

    norb_orig = self.norb_orig
    norb = self.norb
    ns = self.ns

    def _pack_index(x, n):
        _ix, _iy, _iz = x
        _nx, _ny, _nz = n
        return _ix + _nx * (_iy + _ny * (_iz))

    def _unpack_index(x, n):
        _nx, _ny, _nz = n
        _ix = x % _nx
        _iy = (x // _nx) % _ny
        _iz = (x // (_nx * _ny)) % _nz
        return (_ix, _iy, _iz)

    def _pack_site(x, n):
        _ix, _iy, _iz = x
        _nx, _ny, _nz = n
        return _iz + _nz * (_iy + _ny * (_ix))

    def _unpack_site(idx, n):
        _nx, _ny, _nz = n
        _iz = idx % _nz
        _iy = (idx // _nz) % _ny
        _ix = (idx // (_nz * _ny)) % _nx
        return (_ix,_iy,_iz)

    green = np.zeros((Lvol,ns,norb_orig,ns,norb_orig), dtype=np.complex128)

    for jsite in range(Lvol):
        ix, iy, iz = _unpack_site(jsite, (Lx,Ly,Lz))

        ixx, irx = ix // Bx, ix % Bx
        iyy, iry = iy // By, iy % By
        izz, irz = iz // Bz, iz % Bz

        isite = _pack_site((ixx,iyy,izz), (Nx,Ny,Nz))
        ir = _pack_index((irx,iry,irz), (Bx,By,Bz))

        for a, b in itertools.product(range(norb_orig), range(norb_orig)):
            if hamiltonian:
                aa = a
                bb = b + norb_orig * ir
            else:
                aa = a + norb_orig * ir
                bb = b
            for s, t in itertools.product(range(ns), range(ns)):
                green[jsite, s, a, t, b] = green_sub[isite, s, aa, t, bb]

    return green
