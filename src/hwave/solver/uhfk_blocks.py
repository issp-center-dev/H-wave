"""Block structure and the block-projected Hartree-Fock potential of UHFk.

These are the block detection of ``UHFk._detect_blocks`` and the
interaction part of ``UHFk._make_ham``, written as functions of their
inputs so that code outside the solver (``hwave.solver.initial_facade``)
can use them without a solver instance. UHFk calls them; the arithmetic and
its order are those of the methods they came from, so the solver results
are bit-for-bit unchanged.

* :func:`detect_blocks` finds the independent blocks of the ``(nd, nd)``
  index space and the electron reservoir (mu-group) of each block.
* :func:`accumulate_hf_projected` adds the Hartree(-Fock) potential of a
  Green function onto a Hamiltonian, in place, and zeroes the entries
  between different blocks.
"""
import logging
from dataclasses import dataclass
from typing import Mapping, Tuple

import numpy as np

from . import hartree_fock
from . import spin_orbital
from .base import check_spin_flip_with_fixed_sz

# the messages below were logged by UHFk._detect_blocks; they keep its
# logger so that the run log does not change
logger = logging.getLogger("qlms").getChild("uhfk")


# eq=False: the fields hold arrays (identity comparison and hash)
@dataclass(frozen=True, eq=False)
class SolverBlocks:
    """The block structure UHFk solves with.

    Attributes
    ----------
    block_info : tuple of ndarray
        One read-only ``int`` array of ``nd``-space indices per independent
        block, ordered by first index.
    block_to_group : tuple of int
        The mu-group of each block.
    group_nconds : tuple of int
        The electron count of each mu-group.
    """
    block_info: Tuple[np.ndarray, ...]
    block_to_group: Tuple[int, ...]
    group_nconds: Tuple[int, ...]


@dataclass(frozen=True)
class HFTables:
    """The interaction tables of the Hartree-Fock potential.

    Attributes
    ----------
    inter_table, spin_table : Mapping
        Per interaction type, ``J_ab(r)`` and its spin structure, as
        ``hartree_fock.build_interaction_tables`` returns them (``None`` for
        a type that is absent).
    shape : tuple of int
        The k-grid ``(nx, ny, nz)``.
    enable_spin_orbital : bool
        The interleaved spin-orbital layout (``ns = 1``, index
        ``2*orb + spin``); the potential is then built in the spin-block
        frame and mapped back.
    norb_phys : int
        The number of physical orbitals the tables are indexed by.
    """
    inter_table: Mapping
    spin_table: Mapping
    shape: Tuple[int, int, int]
    enable_spin_orbital: bool
    norb_phys: int


def detect_blocks(ham_trans, inter_table, spin_table, *, nd, norb, ns,
                  iflag_fock, sz_free, enable_spin_orbital, nconds):
    """The independent blocks of the one-body plus mean-field Hamiltonian.

    Two indices of the ``(nd, nd)`` space are connected when the transfer
    term (summed over k) or an interaction term can couple them; the blocks
    are the connected components. With a fixed 2Sz (``sz_free`` false, normal
    mode) the spin sectors are never connected and a one-body spin-flip
    term is refused.

    Parameters
    ----------
    ham_trans : ndarray, (nvol, nd, nd)
        The transfer term ``T(k)``.
    inter_table, spin_table : Mapping
        The interaction tables (see :class:`HFTables`).
    nd, norb, ns : int
        The solver dimensions (``nd = ns * norb``).
    iflag_fock : bool
        Whether the Fock term is included.
    sz_free : bool
        Whether 2Sz is free (one mu-group) or fixed (per-spin groups).
    enable_spin_orbital : bool
        The interleaved spin-orbital layout.
    nconds : sequence of int
        The electron counts: ``(Ncond,)`` or ``(N_up, N_down)``.

    Returns
    -------
    SolverBlocks

    Raises
    ------
    ValueError
        A fixed 2Sz with a one-body spin-flip term, or with a spin-mixed
        block next to pure-spin blocks.
    """
    # Build connectivity matrix in nd x nd space
    connectivity = np.zeros((nd, nd), dtype=np.float64)

    # 1. Transfer Hamiltonian: sum |ham_trans(k)| over k
    # (kept as a local: the 2Sz-fixed spin-flip check below reuses it)
    trans_conn = np.sum(np.abs(ham_trans), axis=0)
    connectivity += trans_conn

    # 2. Interaction terms: expand (norb_inter, norb_inter) x spin(2,2,2,2) -> (nd, nd)
    # (spin-orbital mode: nd = 2 * norb_phys)
    norb_inter = nd // 2 if enable_spin_orbital else norb

    for type_name in inter_table:
        if inter_table[type_name] is not None:
            jab = inter_table[type_name]
            spin = spin_table[type_name]

            # Sum |J_ab(r)| over r to get orbital connectivity
            jab_sum = np.sum(np.abs(jab.reshape(-1, norb_inter, norb_inter)), axis=0)

            # For Hartree: spin_mat_h[s,t] = any(spin[s,:,:,t] != 0)
            spin_mat_h = np.any(spin != 0, axis=(1, 2)).astype(float)
            orb_connected_h = (jab_sum > 0).any(axis=1).astype(float)

            if enable_spin_orbital:
                # In spin-orbital mode: (s, a) -> the interleaved index
                for s in range(2):
                    for t in range(2):
                        if spin_mat_h[s, t] > 0:
                            for a in range(norb_inter):
                                if orb_connected_h[a]:
                                    connectivity[spin_orbital.interleaved_index(a, s),
                                                 spin_orbital.interleaved_index(a, t)] += 1.0
                if iflag_fock:
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

                if iflag_fock:
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
    if not sz_free and not enable_spin_orbital and ns == 2:
        check_spin_flip_with_fixed_sz(trans_conn, norb)
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
    if sz_free:
        # All blocks share one global mu
        # group_nconds[g] = target electron count for group g
        # block_to_group[b] = which group block b belongs to
        block_to_group = [0] * len(blocks)
        group_nconds = [nconds[0]]
    else:
        # Classify blocks by spin content
        ncond_up, ncond_down = nconds[0], nconds[1]
        block_to_group = []
        # group 0 = up-only, group 1 = down-only, group 2 = mixed
        group_nconds_dict = {}  # group_id -> ncond
        for blk in blocks:
            if enable_spin_orbital:
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

    for blk in blocks:
        blk.flags.writeable = False
    return SolverBlocks(tuple(blocks), tuple(block_to_group), tuple(group_nconds))


def accumulate_hf_projected(out_k, green, tables, *, include_fock, blocks,
                            debug=None):
    """Add the Hartree(-Fock) potential of ``green`` onto ``out_k`` IN
    PLACE, then zero the entries between different blocks.

    This is the interaction part of ``UHFk._make_ham``: called on
    ``0 + T(k)`` it gives the UHFk Hamiltonian bit for bit.

    Parameters
    ----------
    out_k : ndarray, complex128, (nvol, nd, nd)
        Accumulated into; in the solver's index order (interleaved in
        spin-orbital mode).
    green : ndarray, (nvol, ns, norb, ns, norb)
        The stored Green function ``G[r, s, a, t, b] = <c^dag_{sa}(r)
        c_{tb}(0)>`` (UHFk's ``Green``). Not modified.
    tables : HFTables
    include_fock : bool
        Include the Fock (exchange) terms.
    blocks : SolverBlocks or sequence of index arrays
        The block structure; a single block masks nothing.
    debug : callable, optional
        Called with the name of each interaction type (logging).

    Returns
    -------
    None
    """
    nvol = out_k.shape[0]
    nd = out_k.shape[-1]
    block_info = blocks.block_info if isinstance(blocks, SolverBlocks) else blocks

    # interaction term: the shared kernel adds every interaction
    # contribution IN PLACE, in the original order and with the original
    # expressions (bit-identical in normal mode; #181 Phase B). In
    # spin-orbital mode it runs in the virtual spin-block frame and the
    # sum is mapped back to the interleaved order (#83).
    if not tables.enable_spin_orbital:
        hartree_fock.accumulate_hf(
            out_k, green, tables.inter_table, tables.spin_table, tables.shape,
            include_fock=include_fock, debug=debug)
    else:
        norb_phys = nd // 2
        # the spin-orbital Green function (nvol, 1, nd, 1, nd) in the
        # virtual spin-separate form (nvol, 2, norb_phys, 2, norb_phys)
        g_sb = spin_orbital.interleaved_to_spinblock(
            green.reshape(nvol, nd, nd), norb_phys)
        gab_r = np.ascontiguousarray(g_sb).reshape(nvol, 2, norb_phys, 2, norb_phys)
        # accumulate onto the transfer term reordered to the virtual
        # frame: the reorder is a pure permutation, so every element sees
        # the same additions in the same order as the interleaved loop
        # this replaced (bit-identical)
        ham_virt = np.ascontiguousarray(
            spin_orbital.interleaved_to_spinblock(out_k, norb_phys))
        hartree_fock.accumulate_hf(
            ham_virt, gab_r, tables.inter_table, tables.spin_table, tables.shape,
            include_fock=include_fock, debug=debug)
        out_k[...] = spin_orbital.spinblock_to_interleaved(
            ham_virt.reshape(nvol, nd, nd), norb_phys)

    # Enforce block structure: zero out cross-block entries
    if len(block_info) > 1:
        mask = np.zeros((nd, nd), dtype=bool)
        for blk in block_info:
            idx = np.array(blk)
            ix = np.ix_(idx, idx)
            mask[ix] = True
        out_k[:, ~mask] = 0.0
