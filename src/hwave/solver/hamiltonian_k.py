"""One-body H(k) from a Wannier90-like R-space table (refactoring plan 2.1).

Three steps, shared by UHFk, RPA/FLEX and hwave_sc -- each of which keeps
its own copy (``UHFk._make_ham_trans``, ``rpa.Interaction._make_ham_trans``,
``sc._build_hamiltonian_k``) until plan steps 2.2-2.4 switch them over; the
equivalence tests in ``tests/test_hamiltonian_k.py`` pin this module to
those three element by element.

``dense_table``
    ``{((rx, ry, rz), (a, b)): v}`` -> ``T(R)`` as ``(nx, ny, nz, p, p)``.
    Entries accumulate: a bond ``R`` and its image ``R +- L`` land on the
    same slot when the hopping range reaches the lattice size, and a
    negative ``R`` wraps exactly as numpy's negative indexing did in the
    solver copies.
``to_k``
    ``T(R)`` -> ``T(k) = sum_R e^{+ikR} T(R)`` as ``(nvol, p, p)`` on the FFT
    grid ``k = 2 pi m / N`` in C order -- the documented Wannier90-style
    sign (issue #133), i.e. an unnormalised inverse FFT.
``embed_spin``
    ``(nvol, p, p)`` -> ``(nvol, ns*p, ns*p)``, spin-major and spin-diagonal:
    the ``T ⊗ I_ns`` of a spin-independent table.

``p`` is whatever the table's indices count, and the builder never reorders
them:

* physical orbitals -- the 7-column Transfer in normal mode (``embed_spin``
  then supplies the spin axis);
* interleaved spin-orbitals ``2*orb + spin`` -- ``enable_spin_orbital``
  (UHFk keeps this order; RPA/FLEX and hwave_sc apply
  :func:`hwave.solver.spin_orbital.interleaved_to_spinblock` afterwards);
* spin-block composite indices ``spin*norb + orb`` -- the spin-resolved
  9-column Transfer as folded by the reader of issue #323, which is
  already the spin-major layout ``embed_spin`` produces, so it is used as
  is with ``p = 2*norb``.
"""
import logging

import numpy as np

logger = logging.getLogger("qlms").getChild("hamiltonian_k")


def dense_table(table, shape, dim, *, out_of_range="raise", label="Transfer"):
    """Scatter an R-space table into a dense ``(nx, ny, nz, dim, dim)`` array.

    Parameters
    ----------
    table : Mapping
        ``{((rx, ry, rz), (a, b)): value}`` with 0-based indices, as read by
        ``wan90.read_w90`` and (when a sublattice is used) folded by
        ``fold.reshape_interaction``.
    shape : tuple of int
        ``(nx, ny, nz)``: the k-grid, i.e. the cell shape in units of the
        (folded) unit cell. ``R`` is taken modulo it.
    dim : int
        Size of the index space: every ``a``, ``b`` must lie in ``[0, dim)``.
    out_of_range : {"raise", "drop"}
        What to do with an entry whose index is outside ``[0, dim)``.
        ``"drop"`` skips it and logs one warning per call (the normal-mode
        reading of a spin-orbital-format file, where the spin block of a
        7-column table is ignored); ``"raise"`` refuses it.
    label : str
        Named in the messages.

    Returns
    -------
    ndarray, complex128, shape ``(nx, ny, nz, dim, dim)``
    """
    if out_of_range not in ("raise", "drop"):
        raise ValueError("out_of_range must be 'raise' or 'drop', got {!r}"
                         .format(out_of_range))
    nx, ny, nz = (int(n) for n in shape)
    dim = int(dim)
    tab = np.zeros((nx, ny, nz, dim, dim), dtype=np.complex128)
    dropped = 0
    for (irvec, orbvec), v in table.items():
        a, b = int(orbvec[0]), int(orbvec[1])
        if not (0 <= a < dim and 0 <= b < dim):
            if out_of_range == "drop":
                dropped += 1
                continue
            raise ValueError(
                "{}: orbital index {} (1-based {}) is outside [0, {}) at "
                "R = {}".format(label, (a, b), (a + 1, b + 1), dim,
                                tuple(irvec)))
        rx, ry, rz = (int(r) for r in irvec)
        tab[rx % nx, ry % ny, rz % nz, a, b] += v
    if dropped:
        _warn_dropped(label, dropped, dim)
    return tab


def _warn_dropped(label, dropped, dim):
    # the one wording of the ignore-warning; tests count it by its
    # "These terms are ignored" phrase. Only the Transfer can be read as a
    # spin-orbital table; Extern is physical-indexed in every mode, so its
    # out-of-range rows are refused there rather than used.
    if label.lower() == "transfer":
        hint = "Set enable_spin_orbital = true to use them."
    else:
        hint = ("{} uses physical orbital indices in every mode; remove "
                "these rows.".format(label))
    logger.warning(
        "{} has {} entries with orbital indices >= {} (spin-dependent "
        "terms) but enable_spin_orbital is False. These terms are "
        "ignored. {}".format(label, dropped, dim, hint))


def drop_spin_block(table, dim, *, label="Transfer"):
    """Normalise a normal-mode one-body table once, at the boundary.

    A spin-orbital-format file read with ``enable_spin_orbital = False``
    carries its spin block as orbital indices ``>= dim`` (the physical
    orbital count); the documented behaviour is to ignore those entries
    with a warning. The solvers call this right after reading, so every
    later consumer (hermiticity check, sublattice fold, H(k) builder) sees
    an in-range table and raises on anything else instead of guarding
    individually.

    Returns
    -------
    (table, dropped)
        The input table itself and 0 when nothing is out of range (no
        copy); otherwise a new dict without the spin-block entries and
        their count, after one warning.
    """
    dim = int(dim)
    dropped = sum(1 for key in table
                  if key[1][0] >= dim or key[1][1] >= dim)
    if dropped == 0:
        return table, 0
    _warn_dropped(label, dropped, dim)
    kept = {key: v for key, v in table.items()
            if key[1][0] < dim and key[1][1] < dim}
    return kept, dropped


def to_k(tab_r):
    """``T(k) = sum_R e^{+ikR} T(R)`` on the FFT grid, as ``(nvol, p, q)``.

    ``tab_r`` is ``(nx, ny, nz, p, q)``; the k index runs in C order over
    ``(nx, ny, nz)`` with ``k_i = 2 pi m_i / n_i``. The transform is the
    unnormalised inverse FFT (``ifftn * nvol``), the sign every solver copy
    uses since issue #133.
    """
    tab_r = np.asarray(tab_r)
    nx, ny, nz, p, q = tab_r.shape
    tab_k = np.fft.ifftn(tab_r, axes=(0, 1, 2), norm="forward")
    return tab_k.reshape(nx * ny * nz, p, q)


def embed_spin(mat_k, ns=2):
    """``(nvol, p, q)`` -> ``(nvol, ns*p, ns*q)``: ``mat_k ⊗ I_ns`` with the
    spin index outer (composite index ``s*p + a``), i.e. the same
    spin-independent, spin-diagonal operator on every spin."""
    mat_k = np.asarray(mat_k)
    nvol, p, q = mat_k.shape
    return np.einsum("kab,st->ksatb", mat_k, np.eye(ns)).reshape(
        nvol, ns * p, ns * q)


def one_body_k(table, shape, dim, *, spin_embed=False, out_of_range="raise",
               label="Transfer"):
    """``dense_table`` -> ``to_k`` (-> ``embed_spin``) in one call.

    ``spin_embed=True`` is the normal-mode 7-column Transfer (``dim`` = the
    physical orbital count, result ``(nvol, 2*dim, 2*dim)``); ``False``
    returns ``(nvol, dim, dim)`` for a table whose indices already carry
    the spin (interleaved or spin-block) or for a consumer without a spin
    axis.
    """
    mat_k = to_k(dense_table(table, shape, dim, out_of_range=out_of_range,
                             label=label))
    return embed_spin(mat_k) if spin_embed else mat_k
