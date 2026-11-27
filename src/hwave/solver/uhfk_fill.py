"""Diagonalisation and occupation fill of UHFk, as functions.

These are ``UHFk._diag``, ``UHFk._green`` and ``UHFk._find_dist_group*``
written as functions of their inputs, so that code outside the solver
(``hwave.solver.initial_facade``) can build the solver's Green function
from a Hamiltonian without a solver instance. UHFk calls them; the
arithmetic and its order are those of the methods they came from, so the
solver results are bit-for-bit unchanged.

* :func:`diag_blocks`: the eigen-decomposition of each block of ``H(k)``.
* :func:`fill_green`: one chemical potential per mu-group, the occupations
  and the stored Green function ``G[r, s, a, t, b] = <c^dag_{sa}(r)
  c_{tb}(0)>``, with the solver's fill (``fill="solver"``) or, at T = 0, a
  uniform fill of the Fermi shell (``fill="uniform"``).
* :func:`find_dist_group`, :func:`find_dist_group_zero_t`,
  :func:`find_dist_group_nonzero_t`: the occupations of one mu-group.

At T = 0 the solver fills the ``ncond`` lowest levels of a group, exact
ties ordered by ``|k|^2`` and then by block and band order
(``np.lexsort``). A Fermi level that a degenerate shell straddles is then
filled by that order, which breaks the lattice symmetry and depends on the
eigenbasis ``eigh`` picks inside the shell. ``fill="uniform"`` fills every
level of the shell (the levels within ``shell_tol`` of the Fermi level)
by the same fraction instead, a projector onto whole degenerate subspaces.
"""
import logging
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

# the messages below were logged by the UHFk methods; they keep its logger
# so that the run log does not change
logger = logging.getLogger("qlms").getChild("uhfk")

__all__ = ["ShellRecord", "FillInfo", "MuNotConvergedError", "diag_blocks",
           "group_blocks", "assemble_green", "fill_green", "find_dist_group", "find_dist_group_zero_t",
           "find_dist_group_nonzero_t"]

_FILLS = ("solver", "uniform")


class MuNotConvergedError(RuntimeError):
    """The T > 0 chemical potential search did not converge (UHFk exits
    with status 1 on it)."""


@dataclass(frozen=True)
class ShellRecord:
    """The Fermi shell of one mu-group at T = 0.

    ``group`` is ``"all"`` (one group for both spins, Sz free or a
    spin-mixed group), ``"up"`` or ``"down"`` (2Sz fixed). ``open_shell``
    is true when the shell (the levels within ``shell_tol`` of the Fermi
    level) is only partly filled; ``shell_degeneracy`` (its number of
    levels) and ``shell_fill`` (the occupation ``(n_group - n_below) /
    g_shell`` of each of its levels under the uniform fill) are given only
    then, ``None`` otherwise. Always closed at T > 0.
    """
    group: str
    open_shell: bool
    shell_degeneracy: Optional[int]
    shell_fill: Optional[float]


@dataclass(frozen=True)
class FillInfo:
    """The chemical potential (``mu``) and the Fermi shell (``shell``) of
    each mu-group, in group order."""
    mu: Tuple[float, ...]
    shell: Tuple[ShellRecord, ...]


def diag_blocks(ham, blocks):
    """The eigenvalues and eigenvectors of each block of ``ham`` ``(nvol,
    nd, nd)``: two lists, one ``(nvol, n_b)`` and one ``(nvol, n_b, n_b)``
    array per block (``UHFk._diag``)."""
    # Diagonalize each block independently
    eigenvalues = []
    eigenvectors = []
    for blk in blocks:
        idx = np.array(blk)
        ix = np.ix_(idx, idx)
        mat_blk = ham[:, ix[0], ix[1]]  # (nvol, blk_size, blk_size)
        w, v = np.linalg.eigh(mat_blk)
        eigenvalues.append(w)
        eigenvectors.append(v)
    return eigenvalues, eigenvectors


def _ksq_table(shape, width):
    nx, ny, nz = shape
    nvol = nx * ny * nz

    kx = np.roll( (np.arange(nx) - nx//2), -nx//2) ** 2
    ky = np.roll( (np.arange(ny) - ny//2), -ny//2) ** 2
    kz = np.roll( (np.arange(nz) - nz//2), -nz//2) ** 2

    rr = np.zeros((nx,ny,nz))
    rr += np.broadcast_to(kx.reshape(nx,1,1),(nx,ny,nz))
    rr += np.broadcast_to(ky.reshape(1,ny,1),(nx,ny,nz))
    rr += np.broadcast_to(kz.reshape(1,1,nz),(nx,ny,nz))

    return np.broadcast_to(rr.reshape(nvol,1), (nvol,width))


def find_dist_group_zero_t(ws_list, vs_list, ncond, *, shape):
    """The T = 0 occupations of one mu-group (``UHFk._find_dist_group_zero_t``):
    the ``ncond`` lowest levels of the group, exact ties by ``|k|^2``, then
    by concatenation order. Returns ``(dists, mu)``, ``mu`` the midpoint of
    the HOMO-LUMO gap."""
    # Collect all eigenvalues across blocks (concatenated in block order)
    all_ww = []
    all_ksq = []
    for b, w in enumerate(ws_list):
        k_sq = _ksq_table(shape, w.shape[1]).flatten()
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
    dists = _split(occ, ws_list)

    return dists, _zero_t_mu(all_ww, ncond)


def _split(occ, ws_list):
    block_sizes = [w.size for w in ws_list]
    split_points = np.cumsum(block_sizes)[:-1]
    return [
        d.reshape(w.shape)
        for d, w in zip(np.split(occ, split_points), ws_list)
    ]


def _zero_t_mu(all_ww, ncond):
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
    return mu


def find_dist_group_nonzero_t(ws_list, vs_list, ncond, *, T, ene_cutoff):
    """The T > 0 occupations of one mu-group
    (``UHFk._find_dist_group_nonzero_t``): Fermi functions at the ``mu``
    that holds ``ncond`` electrons (bisection, then Newton). Returns
    ``(dists, mu)``; raises :class:`MuNotConvergedError` when neither
    converges."""
    from scipy import optimize

    # Collect all eigenvalues for bracket search
    all_ev = np.sort(np.concatenate([w.flatten() for w in ws_list]))
    occupied_number = ncond

    def _fermi(t, mu, ev):
        w_ = (ev - mu) / t
        mask_ = w_ < ene_cutoff
        w1_ = np.where( mask_, w_, 0.0 )
        v1_ = 1.0 / (1.0 + np.exp(w1_))
        v_ = np.where( mask_, v1_, 0.0 )
        return v_

    def _calc_delta_n(mu):
        nn = 0.0
        for w, v in zip(ws_list, vs_list):
            ff = _fermi(T, mu, w)
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
        raise MuNotConvergedError("find mu: not converged")

    dists = [_fermi(T, mu, w) for w in ws_list]
    return dists, mu


def find_dist_group(ws_list, vs_list, ncond, *, shape, T, ene_cutoff):
    """The occupations and ``mu`` of one mu-group
    (``UHFk._find_dist_group``)."""
    if T == 0:
        return find_dist_group_zero_t(ws_list, vs_list, ncond, shape=shape)
    else:
        return find_dist_group_nonzero_t(ws_list, vs_list, ncond, T=T,
                                         ene_cutoff=ene_cutoff)


def _group_label(blocks, blist, ns, norb, ngroup):
    """``"up"``/``"down"`` for a group of pure-spin blocks of a normal-mode
    solver with more than one group, else ``"all"``."""
    if ns != 2 or ngroup < 2:
        return "all"
    idx = np.concatenate([np.asarray(blocks[b]) for b in blist])
    if np.all(idx < norb):
        return "up"
    if np.all(idx >= norb):
        return "down"
    return "all"


def _shell(ws_group, ncond, shell_tol):
    """``(n_below, g_shell, in_shell mask, below mask)`` of the T = 0 Fermi
    shell of one group (the levels within ``shell_tol`` of the HOMO), on
    the concatenated eigenvalues; ``None`` when the group is empty or
    full (no Fermi level inside the spectrum)."""
    all_ww = np.concatenate([w.flatten() for w in ws_group])
    if ncond <= 0 or ncond >= all_ww.size:
        return None
    e_f = np.sort(all_ww)[ncond - 1]
    in_shell = np.abs(all_ww - e_f) <= shell_tol
    below = all_ww < e_f - shell_tol
    return int(np.count_nonzero(below)), int(np.count_nonzero(in_shell)), in_shell, below


def group_blocks(block_to_group, ngroup):
    """The block indices of each mu-group, in block order."""
    groups = [[] for _ in range(ngroup)]
    for b in range(len(block_to_group)):
        groups[block_to_group[b]].append(b)
    return groups


def assemble_green(vs_list, dists, blocks, *, shape, ns, norb):
    """The stored ``G(r)`` ``(nvol, ns, norb, ns, norb)`` of the per-block
    eigenvectors ``vs_list`` and occupations ``dists`` (the assembly and
    Fourier transform of ``UHFk._green``)."""
    nx, ny, nz = shape
    nvol = nx * ny * nz
    nd = ns * norb

    # Build full Green's function in k-space by assembling blocks
    gab_k = np.zeros((nvol, nd, nd), dtype=np.complex128)

    for k in range(len(blocks)):
        v = vs_list[k]
        dist = dists[k]

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

    return gab_r.reshape(nvol, ns, norb, ns, norb)


def fill_green(ws_list, vs_list, blocks, block_to_group, group_nconds, *, shape,
               ns, norb, T, ene_cutoff, fill="solver", shell_tol=None,
               return_info=False):
    """The stored Green function of the eigen-decomposition ``ws_list``,
    ``vs_list`` (:func:`diag_blocks`) of the blocks ``blocks``.

    With ``fill="solver"`` this is ``UHFk._green``: ``(green, mu)``,
    ``green`` ``(nvol, ns, norb, ns, norb)`` and ``mu`` one float per
    mu-group (an ``ndarray``). With ``fill="uniform"`` at T = 0 the levels
    within ``shell_tol`` of the Fermi level (the ``ncond``-th lowest level
    of the group) form the shell; the levels below it are filled and the
    shell by ``(n_group - n_below) / g_shell`` each. At T > 0 the uniform
    fill is the solver's.

    With ``return_info`` a :class:`FillInfo` is appended to the result. At
    T = 0 and ``fill="solver"`` its shell is that of ``shell_tol`` (exact
    ties when ``shell_tol`` is ``None``); the solver's fill does not change
    with it.
    """
    if fill not in _FILLS:
        raise ValueError("fill must be one of {}, got {!r}".format(_FILLS, fill))
    if shell_tol is not None and not shell_tol >= 0.0:
        raise ValueError("shell_tol must be a non-negative number, got {!r}"
                         .format(shell_tol))
    uniform = fill == "uniform" and T == 0
    if uniform and shell_tol is None:
        raise ValueError("fill='uniform' at T = 0 needs shell_tol (the width of "
                         "the Fermi shell)")
    tol = 0.0 if shell_tol is None else float(shell_tol)

    nblock = len(blocks)
    ngroup = len(group_nconds)

    # Find a single global mu per group (shared across all blocks in the group)
    mus = np.zeros(ngroup, dtype=float)

    # Group blocks by mu-group
    groups = group_blocks(block_to_group, ngroup)

    # Find mu for each group using all eigenvalues in that group
    group_dists = [None] * nblock
    shells = []
    for g in range(ngroup):
        blist = groups[g]
        ws_group = [ws_list[b] for b in blist]
        vs_group = [vs_list[b] for b in blist]
        ncond = group_nconds[g]

        dists, mu = find_dist_group(ws_group, vs_group, ncond, shape=shape, T=T,
                                    ene_cutoff=ene_cutoff)

        record = None
        if (uniform or return_info) and T == 0:
            sh = _shell(ws_group, ncond, tol)
            label = _group_label(blocks, blist, ns, norb, ngroup)
            if sh is not None and sh[0] < ncond < sh[0] + sh[1]:
                n_below, g_shell, in_shell, below = sh
                f_shell = (ncond - n_below) / g_shell
                record = ShellRecord(label, True, g_shell, f_shell)
                if uniform:
                    occ = np.zeros(in_shell.size)
                    occ[below] = 1.0
                    occ[in_shell] = f_shell
                    dists = _split(occ, ws_group)
            else:
                record = ShellRecord(label, False, None, None)
        elif return_info:
            record = ShellRecord(_group_label(blocks, blist, ns, norb, ngroup),
                                 False, None, None)
        shells.append(record)

        mus[g] = mu
        logger.debug("mu[group {}] = {}".format(g, mu))
        for i, b in enumerate(blist):
            group_dists[b] = dists[i]

    green = assemble_green(vs_list, group_dists, blocks, shape=shape, ns=ns, norb=norb)
    if return_info:
        return green, mus, FillInfo(mu=tuple(float(m) for m in mus), shell=tuple(shells))
    return green, mus
