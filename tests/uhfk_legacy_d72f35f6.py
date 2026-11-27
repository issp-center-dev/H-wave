"""Frozen verbatim copies of UHFk methods at revision d72f35f6.

These are oracles for the bitwise-equality tests of code lifted out of
``UHFk``: each function body is the d72f35f6 method body, unchanged, taking
the solver instance (or any object with the same attributes) as ``self``.
Do not edit them; compare new code against them in the same process.

Revision: d72f35f6 (``src/hwave/solver/uhfk.py``: the fold methods, the
block detection, ``_make_ham``, ``_diag``, ``_green`` and
``_find_dist_group*``; and, for ``check_spin_flip_with_fixed_sz``,
``src/hwave/solver/base.py``).

The bodies call other private methods through ``self`` (``_make_ham``
calls ``_so_to_virtual_green``, ``_detect_blocks`` calls
``_check_spin_flip_with_fixed_sz``). :class:`LegacyView` wraps a solver so
that those calls reach the frozen copies too, and so that the attributes the
bodies assign (``block_info``, ``ham``, ...) land on the view, not on the
solver.
"""
import itertools
import logging

import numpy as np

from hwave.solver import hartree_fock
from hwave.solver import spin_orbital

logger = logging.getLogger("qlms").getChild("uhfk")


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


def check_spin_flip_with_fixed_sz(self, trans_abs, norb):
    """Reject one-body spin-flip terms in 2Sz-fixed mode.

    Such terms genuinely break Sz conservation, so a fixed 2Sz is
    ill-defined.  Numerical noise must not abort the run: cross-spin
    entries are compared against the overall transfer scale (Wannier
    transfers carry ~1e-14 cross-spin residues, and a k-sum amplifies
    them past any absolute epsilon), so entries below 1e-8 of the scale
    are treated as zero (and masked by the block structure downstream).

    Parameters
    ----------
    trans_abs : ndarray, shape (2*norb, 2*norb)
        Non-negative magnitudes of the one-body terms in spin-block
        order (possibly summed over k).
    norb : int
        Number of orbitals per spin sector.

    Raises
    ------
    ValueError
        If a cross-spin entry exceeds 1e-8 of the largest entry.
    """
    scale = trans_abs.max() if trans_abs.size > 0 else 0.0
    if scale <= 0.0:
        return
    cross = max(trans_abs[:norb, norb:].max(initial=0.0),
                trans_abs[norb:, :norb].max(initial=0.0))
    if cross > 1.0e-8 * scale:
        msg = ("2Sz-fixed mode is incompatible with spin-flip one-body "
               "terms: Sz is not conserved, so a fixed 2Sz is "
               "ill-defined. Use the Sz-free mode (omit 2Sz) for such "
               "systems.")
        logger.error(msg)
        raise ValueError(msg)
    if cross > 0.0:
        # nonzero but below the noise threshold: masked, not fatal --
        # say so, in case it is a deliberately weak physical term
        logger.warning(
            "2Sz-fixed mode: cross-spin one-body entries up to {:.3e} "
            "({:.1e} relative to the transfer scale) are below the "
            "noise threshold (1e-8) and are ignored. If they are "
            "physical (e.g. weak spin-orbit coupling), use the Sz-free "
            "mode.".format(cross, cross / scale))


def detect_blocks(self):
    """Detect block-diagonal structure from transfer and interaction terms.

    Analyzes the connectivity in the (spin, orbital) index space by
    combining the non-zero patterns of:
    1. Transfer Hamiltonian ham_trans (summed over k-points)
    2. Interaction terms inter_table x spin_table (expanded to nd x nd)

    The detected blocks are used by _diag() and _green() to decompose
    the problem into smaller independent sub-problems.

    Sets self.block_info: list of arrays, each containing nd-space indices
    for one independent block. Also sets self.group_nconds and
    self.block_to_group for global chemical potential determination.
    """
    logger.debug(">>> _detect_blocks")

    nd = self.nd
    norb = self.norb
    ns = self.ns

    # Build connectivity matrix in nd x nd space
    connectivity = np.zeros((nd, nd), dtype=np.float64)

    # 1. Transfer Hamiltonian: sum |ham_trans(k)| over k
    # (kept as a local: the 2Sz-fixed spin-flip check below reuses it)
    trans_conn = np.sum(np.abs(self.ham_trans), axis=0)
    connectivity += trans_conn

    # 2. Interaction terms: expand (norb_inter, norb_inter) x spin(2,2,2,2) -> (nd, nd)
    norb_inter = self.norb_phys if self.enable_spin_orbital else norb

    for type_name in self.inter_table:
        if self.inter_table[type_name] is not None:
            jab = self.inter_table[type_name]
            spin = self.spin_table[type_name]

            # Sum |J_ab(r)| over r to get orbital connectivity
            jab_sum = np.sum(np.abs(jab.reshape(-1, norb_inter, norb_inter)), axis=0)

            # For Hartree: spin_mat_h[s,t] = any(spin[s,:,:,t] != 0)
            spin_mat_h = np.any(spin != 0, axis=(1, 2)).astype(float)
            orb_connected_h = (jab_sum > 0).any(axis=1).astype(float)

            if self.enable_spin_orbital:
                # In spin-orbital mode: (s, a) -> the interleaved index
                for s in range(2):
                    for t in range(2):
                        if spin_mat_h[s, t] > 0:
                            for a in range(norb_inter):
                                if orb_connected_h[a]:
                                    connectivity[spin_orbital.interleaved_index(a, s),
                                                 spin_orbital.interleaved_index(a, t)] += 1.0
                if self.iflag_fock:
                    spin_fock = np.max(np.abs(spin), axis=(1, 3))
                    for s in range(2):
                        for t in range(2):
                            if spin_fock[s, t] > 0:
                                for a in range(norb_inter):
                                    for b in range(norb_inter):
                                        if jab_sum[a, b] > 0:
                                            connectivity[spin_orbital.interleaved_index(a, s),
                                                                 spin_orbital.interleaved_index(b, t)] += jab_sum[a, b]
            else:
                # Normal mode: (s, a) -> nd_idx = s * norb + a
                for s in range(ns):
                    for t in range(ns):
                        if spin_mat_h[s, t] > 0:
                            for a in range(norb):
                                if orb_connected_h[a]:
                                    connectivity[s * norb + a, t * norb + a] += 1.0

                if self.iflag_fock:
                    spin_mat_f = np.any(
                        spin.transpose(0, 2, 1, 3).reshape(ns * ns, ns * ns) != 0
                    ).astype(float) if ns > 1 else np.ones((1, 1))
                    spin_fock = np.max(np.abs(spin), axis=(1, 3))
                    for s in range(ns):
                        for t in range(ns):
                            if spin_fock[s, t] > 0:
                                connectivity[s * norb:(s + 1) * norb,
                                             t * norb:(t + 1) * norb] += jab_sum

    # 2Sz-fixed mode: each spin sector carries its own electron count
    # (Ncond +- 2Sz)/2, so the up and down sectors are diagonalized
    # independently and spin off-diagonal (Fock) entries are not used,
    # exactly as in the pre-block full-matrix implementation
    # (einsum 'ksasb->skab' kept only the spin-diagonal blocks).
    # Interaction terms only produce spin off-diagonal Fock entries via
    # the spin off-diagonal Green components, which vanish identically in
    # the spin-diagonal state enforced by the constraint -- so their
    # potential cross-spin connectivity is masked out here.  Without the
    # mask, such interaction patterns would merge the two sectors into
    # one mu-group and the 2Sz constraint would be silently lost.
    # One-body spin-flip terms, in contrast, genuinely break Sz
    # conservation, so a fixed 2Sz is ill-defined: reject them.
    if not self.sz_free and not self.enable_spin_orbital and ns == 2:
        self._check_spin_flip_with_fixed_sz(trans_conn, norb)
        spin_diag = np.zeros((nd, nd), dtype=bool)
        spin_diag[:norb, :norb] = True
        spin_diag[norb:, norb:] = True
        connectivity[~spin_diag] = 0.0

    # Symmetrize
    connectivity = connectivity + connectivity.T

    # Find connected components
    threshold = 1.0e-12
    adj = np.abs(connectivity) > threshold

    # BFS/label propagation
    labels = np.arange(nd)
    changed = True
    while changed:
        changed = False
        for i in range(nd):
            neighbors = np.where(adj[i])[0]
            if len(neighbors) > 0:
                min_label = min(labels[i], np.min(labels[neighbors]))
                if labels[i] != min_label:
                    labels[i] = min_label
                    changed = True
                for n in neighbors:
                    if labels[n] != min_label:
                        labels[n] = min_label
                        changed = True

    unique_labels = np.unique(labels)
    blocks = [np.where(labels == lbl)[0] for lbl in unique_labels]

    # Sort blocks by first index
    blocks.sort(key=lambda b: b[0])

    # Group blocks by spin sector for shared chemical potential.
    # Instead of distributing Ncond per block (which can violate the
    # global energy minimum), we group blocks that share the same
    # electron reservoir and find a single mu per group.
    if self.sz_free:
        # All blocks share one global mu
        # group_nconds[g] = target electron count for group g
        # block_to_group[b] = which group block b belongs to
        block_to_group = [0] * len(blocks)
        group_nconds = [self.Nconds[0]]
    else:
        # Classify blocks by spin content
        ncond_up, ncond_down = self.Nconds[0], self.Nconds[1]
        block_to_group = []
        # group 0 = up-only, group 1 = down-only, group 2 = mixed
        group_nconds_dict = {}  # group_id -> ncond
        for blk in blocks:
            if self.enable_spin_orbital:
                n_up = np.sum(spin_orbital.spin_of(np.asarray(blk)) == 0)
                n_down = np.sum(spin_orbital.spin_of(np.asarray(blk)) == 1)
            else:
                n_up = np.sum(blk < norb)
                n_down = np.sum(blk >= norb)
            if n_up > 0 and n_down == 0:
                block_to_group.append(0)
                group_nconds_dict[0] = ncond_up
            elif n_up == 0 and n_down > 0:
                block_to_group.append(1)
                group_nconds_dict[1] = ncond_down
            else:
                # Mixed block (couples up and down indices). A single mixed
                # block spanning the whole system is fine -- it just takes
                # the total ncond. But if a mixed block COEXISTS with pure
                # up/down blocks, assigning each group its own (global)
                # ncond double-counts the electrons (sum > Ncond), and the
                # per-spin split cannot be enforced on the mixed block.
                block_to_group.append(2)
                group_nconds_dict[2] = ncond_up + ncond_down
        # Guard the genuine over-count: a mixed group together with a pure
        # up/down group.
        if 2 in block_to_group and (0 in block_to_group or 1 in block_to_group):
            logger.error(
                "2Sz-fixed mode found a spin-mixed block coexisting with "
                "pure-spin blocks; the per-group electron counts would "
                "double-count (sum != Ncond). Use the Sz-free mode (global "
                "chemical potential) for such systems.")
            raise ValueError(
                "2Sz-fixed mode cannot mix a spin-coupled block with "
                "pure-spin blocks; use Sz-free mode.")
        # Build ordered list: group_nconds[g] for each unique group
        unique_groups = sorted(set(block_to_group))
        group_remap = {g: i for i, g in enumerate(unique_groups)}
        block_to_group = [group_remap[g] for g in block_to_group]
        group_nconds = [group_nconds_dict[g] for g in unique_groups]

    self.block_info = blocks
    self.block_to_group = block_to_group
    self.group_nconds = group_nconds

    # Log detected structure
    if len(blocks) == 1:
        logger.info("Block detection: single block (nd={})".format(nd))
    else:
        block_desc = ", ".join(
            ["{}".format(len(b)) for b in blocks]
        )
        logger.info("Block detection: {} blocks of sizes [{}], "
                    "{} mu-group(s) with Nconds={}".format(
                        len(blocks), block_desc,
                        len(group_nconds), group_nconds))


def so_to_virtual_green(self, G_so):
    """Spin-orbital Green function ``(nvol, 1, nd, 1, nd)`` (interleaved
    index ``2*orb + spin``) -> the virtual spin-separate form
    ``(nvol, 2, norb_phys, 2, norb_phys)``: the spin-block reorder of
    :mod:`hwave.solver.spin_orbital`, so the normal-mode interaction
    contractions apply unchanged."""
    nvol = G_so.shape[0]
    norb_phys = self.nd // 2
    G_sb = spin_orbital.interleaved_to_spinblock(
        G_so.reshape(nvol, self.nd, self.nd), norb_phys)
    return np.ascontiguousarray(G_sb).reshape(nvol, 2, norb_phys, 2, norb_phys)


def virtual_ham_to_so(self, H_virt):
    """Inverse of :meth:`_so_to_virtual_green` for a Hamiltonian
    ``(nvol, 2, norb_phys, 2, norb_phys)`` -> ``(nvol, nd, nd)`` in the
    interleaved order."""
    nvol = H_virt.shape[0]
    return np.ascontiguousarray(spin_orbital.spinblock_to_interleaved(
        H_virt.reshape(nvol, self.nd, self.nd), self.nd // 2))


def make_ham(self):
    logger.debug(">>> _make_ham")

    nx,ny,nz = self.shape
    nvol     = self.nvol
    nd       = self.nd
    ns       = self.ns

    # In spin-orbital mode, use virtual (ns=2, norb_phys) form for interactions
    # green function G_{ab,st}(r) : gab_r(r,s,a,t,b)
    if self.enable_spin_orbital:
        gab_r = self._so_to_virtual_green(self.Green)
    else:
        gab_r = self.Green

    # hamiltonian H_{ab,st}(k) : ham(k,(s,a),(t,b))
    ham = np.zeros((nvol,nd,nd), dtype=np.complex128)

    # transfer term  T_{ab}(k) (note convention)
    logger.debug("Transfer")
    ham += self.ham_trans

    # interaction term: the shared kernel adds every interaction
    # contribution IN PLACE, in the original order and with the original
    # expressions (bit-identical in normal mode; #181 Phase B). In
    # spin-orbital mode it runs in the virtual spin-block frame and the
    # sum is mapped back to the interleaved order (#83).
    if not self.enable_spin_orbital:
        hartree_fock.accumulate_hf(
            ham, gab_r, self.inter_table, self.spin_table, self.shape,
            include_fock=self.iflag_fock, debug=logger.debug)
    else:
        # accumulate onto the transfer term reordered to the virtual
        # frame: the reorder is a pure permutation, so every element sees
        # the same additions in the same order as the interleaved loop
        # this replaced (bit-identical)
        ham_virt = np.ascontiguousarray(
            spin_orbital.interleaved_to_spinblock(ham, nd // 2))
        hartree_fock.accumulate_hf(
            ham_virt, gab_r, self.inter_table, self.spin_table, self.shape,
            include_fock=self.iflag_fock, debug=logger.debug)
        ham = self._virtual_ham_to_so(
            ham_virt.reshape(nvol, 2, nd // 2, 2, nd // 2))

    # Enforce block structure: zero out cross-block entries
    if len(self.block_info) > 1:
        mask = np.zeros((nd, nd), dtype=bool)
        for blk in self.block_info:
            idx = np.array(blk)
            ix = np.ix_(idx, idx)
            mask[ix] = True
        ham[:, ~mask] = 0.0

    # store
    self.ham = ham


def diag(self):
    logger.debug(">>> _diag")

    nvol = self.nvol
    blocks = self.block_info

    # Diagonalize each block independently
    eigenvalues = []
    eigenvectors = []
    for blk in blocks:
        idx = np.array(blk)
        ix = np.ix_(idx, idx)
        mat_blk = self.ham[:, ix[0], ix[1]]  # (nvol, blk_size, blk_size)
        w, v = np.linalg.eigh(mat_blk)
        eigenvalues.append(w)
        eigenvectors.append(v)

    # Store as list of per-block arrays
    self._green_list["eigenvalue"] = eigenvalues
    self._green_list["eigenvector"] = eigenvectors

def green(self):
    logger.debug(">>> _green")

    nx,ny,nz = self.shape
    nvol     = self.nvol
    nd       = self.nd
    norb     = self.norb
    ns       = self.ns

    ws_list = self._green_list["eigenvalue"]
    vs_list = self._green_list["eigenvector"]
    blocks = self.block_info
    nblock = len(blocks)
    block_to_group = self.block_to_group
    group_nconds = self.group_nconds
    ngroup = len(group_nconds)

    # Find a single global mu per group (shared across all blocks in the group)
    self._green_list["mu"] = np.zeros(ngroup, dtype=float)

    # Group blocks by mu-group
    group_blocks = [[] for _ in range(ngroup)]
    for b in range(nblock):
        group_blocks[block_to_group[b]].append(b)

    # Find mu for each group using all eigenvalues in that group
    group_dists = [None] * nblock
    for g in range(ngroup):
        blist = group_blocks[g]
        ws_group = [ws_list[b] for b in blist]
        vs_group = [vs_list[b] for b in blist]
        ncond = group_nconds[g]

        dists, mu = self._find_dist_group(ws_group, vs_group, ncond)

        self._green_list["mu"][g] = mu
        logger.debug("mu[group {}] = {}".format(g, mu))
        for i, b in enumerate(blist):
            group_dists[b] = dists[i]

    # Build full Green's function in k-space by assembling blocks
    gab_k = np.zeros((nvol, nd, nd), dtype=np.complex128)

    for k in range(nblock):
        v = vs_list[k]
        dist = group_dists[k]

        # G_ab(k) for this block
        #   G_ab(k) = sum_l conj(V_kal) dist_kl V_kbl
        #   = Vconj @ diag(dist) @ V^T, batched over k. Explicit matmul lowers
        #   to batched GEMM (numpy einsum does not); numerically equivalent to
        #   the 'kal,kl,kbl->kab' einsum to ~1e-14 (reduction order).
        gg_blk = np.matmul(np.conjugate(v) * dist[:, None, :], v.transpose(0, 2, 1))

        # Place into full matrix
        idx = np.array(blocks[k])
        ix = np.ix_(idx, idx)
        gab_k[:, ix[0], ix[1]] = gg_blk

    # G_ab(r) = 1/V sum_k G_ab(k) e^{-ikr}
    gab_r = np.fft.fftn(gab_k.reshape(nx, ny, nz, nd, nd),
                        axes=(0, 1, 2), norm='forward')

    # store
    self.Green_prev = self.Green
    self.Green = gab_r.reshape(nvol, ns, norb, ns, norb)

def find_dist_group(self, ws_list, vs_list, ncond):
    """Find occupation distribution using a single mu across multiple blocks.

    Parameters
    ----------
    ws_list : list of ndarray
        Eigenvalues for each block in the group, shape (nvol, block_size)
    vs_list : list of ndarray
        Eigenvectors for each block, shape (nvol, block_size, block_size)
    ncond : int
        Total electron count for the entire group

    Returns
    -------
    dists : list of ndarray
        Occupation distribution for each block
    mu : float
        Chemical potential
    """
    if self.T == 0:
        return self._find_dist_group_zero_t(ws_list, vs_list, ncond)
    else:
        return self._find_dist_group_nonzero_t(ws_list, vs_list, ncond)

def find_dist_group_zero_t(self, ws_list, vs_list, ncond):
    def _ksq_table(width):
        nx,ny,nz = self.shape
        nvol = self.nvol

        kx = np.roll( (np.arange(nx) - nx//2), -nx//2) ** 2
        ky = np.roll( (np.arange(ny) - ny//2), -ny//2) ** 2
        kz = np.roll( (np.arange(nz) - nz//2), -nz//2) ** 2

        rr = np.zeros((nx,ny,nz))
        rr += np.broadcast_to(kx.reshape(nx,1,1),(nx,ny,nz))
        rr += np.broadcast_to(ky.reshape(1,ny,1),(nx,ny,nz))
        rr += np.broadcast_to(kz.reshape(1,1,nz),(nx,ny,nz))

        return np.broadcast_to(rr.reshape(nvol,1), (nvol,width))

    # Collect all eigenvalues across blocks (concatenated in block order)
    all_ww = []
    all_ksq = []
    for b, w in enumerate(ws_list):
        k_sq = _ksq_table(w.shape[1]).flatten()
        ww = w.flatten()
        all_ww.append(ww)
        all_ksq.append(k_sq)

    all_ww = np.concatenate(all_ww)
    all_ksq = np.concatenate(all_ksq)

    # Sort globally and pick lowest ncond states
    ev_idx = np.lexsort((all_ksq, all_ww))[0:ncond]

    # Distribute back to per-block arrays.
    # The concatenated index space maps to (block b, local li) by the block
    # offsets used during concatenation: a global index in [off_b, off_b+size_b)
    # belongs to block b at local index (idx - off_b). Occupied states are
    # distinct, so a single fancy-index assignment + split is exactly the
    # per-element scatter loop, vectorized.
    occ = np.zeros(all_ww.size)
    occ[ev_idx] = 1.0
    block_sizes = [w.size for w in ws_list]
    split_points = np.cumsum(block_sizes)[:-1]
    dists = [
        d.reshape(w.shape)
        for d, w in zip(np.split(occ, split_points), ws_list)
    ]

    # Chemical potential (Fermi level) at T=0: midpoint of the HOMO-LUMO
    # gap. This is the T -> 0+ limit of the finite-temperature mu equation.
    sorted_ww = np.sort(all_ww)
    ntot = sorted_ww.size
    if ncond <= 0:
        mu = sorted_ww[0]
    elif ncond >= ntot:
        mu = sorted_ww[-1]
    else:
        mu = 0.5 * (sorted_ww[ncond - 1] + sorted_ww[ncond])

    return dists, mu

def find_dist_group_nonzero_t(self, ws_list, vs_list, ncond):
    from scipy import optimize

    # Collect all eigenvalues for bracket search
    all_ev = np.sort(np.concatenate([w.flatten() for w in ws_list]))
    occupied_number = ncond

    def _fermi(t, mu, ev):
        w_ = (ev - mu) / t
        mask_ = w_ < self.ene_cutoff
        w1_ = np.where( mask_, w_, 0.0 )
        v1_ = 1.0 / (1.0 + np.exp(w1_))
        v_ = np.where( mask_, v1_, 0.0 )
        return v_

    def _calc_delta_n(mu):
        nn = 0.0
        for w, v in zip(ws_list, vs_list):
            ff = _fermi(self.T, mu, w)
            nn += np.einsum('kal,kl,kal->', np.conjugate(v), ff, v).real
        return nn - occupied_number

    # find mu s.t. <n>(mu) = N0
    is_converged = False
    if (_calc_delta_n(all_ev[0]) * _calc_delta_n(all_ev[-1])) < 0.0:
        logger.debug("+++ find mu: try bisection")
        mu, r = optimize.bisect(_calc_delta_n, all_ev[0], all_ev[-1], full_output=True, disp=False)
        is_converged = r.converged
    if not is_converged:
        logger.debug("+++ find mu: try newton")
        mu, r = optimize.newton(_calc_delta_n, all_ev[0], full_output=True)
        is_converged = r.converged
    if not is_converged:
        logger.error("find mu: not converged. abort")
        exit(1)

    dists = [_fermi(self.T, mu, w) for w in ws_list]
    return dists, mu


class LegacyView:
    """``solver`` seen through the frozen d72f35f6 methods.

    Attribute reads fall through to ``solver``; attribute writes stay on the
    view, so running a frozen method leaves ``solver`` untouched."""

    def __init__(self, solver):
        object.__setattr__(self, "_solver", solver)
        # the bodies write into the _green_list dict: give the view its own
        if hasattr(solver, "_green_list"):
            object.__setattr__(self, "_green_list", dict(solver._green_list))

    def __getattr__(self, name):
        return getattr(object.__getattribute__(self, "_solver"), name)

    _check_spin_flip_with_fixed_sz = check_spin_flip_with_fixed_sz
    _detect_blocks = detect_blocks
    _so_to_virtual_green = so_to_virtual_green
    _virtual_ham_to_so = virtual_ham_to_so
    _make_ham = make_ham
    _diag = diag
    _green = green
    _find_dist_group = find_dist_group
    _find_dist_group_zero_t = find_dist_group_zero_t
    _find_dist_group_nonzero_t = find_dist_group_nonzero_t
