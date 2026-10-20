"""Spin-orbital (``[mode] enable_spin_orbital``) conventions shared by the
solvers (GitHub issue #83, step 2 of the spin-orbital roadmap).

In spin-orbital mode the spin is folded into the orbital index:

* the Geometry ``norb`` and the Transfer file count SPIN-ORBITALS, in the
  INTERLEAVED order ``i = 2*orb + spin`` (Wannier90 SOI, UHFk's native
  order);
* the two-body interaction files are written in PHYSICAL orbital indices
  ``0 <= orb < norb // 2``;
* RPA and FLEX work internally, and write their susceptibilities, in the
  SPIN-BLOCK order ``j = spin*norb_phys + orb`` (tagged
  ``index_convention = "spin_block"``); UHFk evaluates its interaction
  terms in the same spin-block ("virtual spin") frame.

Everything that converts between these orders, counts physical orbitals or
reads the switch lives here, so every solver agrees on it.
"""
import logging

import numpy as np
from requests.structures import CaseInsensitiveDict

from .backend import strict_bool

#: The ``index_convention`` tag RPA/FLEX write into chi0q/chiq files.
INDEX_CONVENTION = "spin_block"

_FLAG = "enable_spin_orbital"

logger = logging.getLogger(__name__)


def resolve_flag(info_mode, name="[mode] enable_spin_orbital", *, allow_float=False):
    """``enable_spin_orbital`` of a ``[mode]`` dict as a real boolean.

    The key lookup is case-insensitive and the value is parsed strictly
    (:func:`hwave.solver.backend.strict_bool`): a misspelled value raises
    instead of silently reading as true or false. A recognized string is
    accepted with a warning. ``allow_float`` additionally accepts an exact
    float 0.0/1.0 -- for UHFk only, whose own parser accepted them before
    the switch was shared; every other caller rejects floats."""
    if info_mode is None:
        return False
    value = CaseInsensitiveDict(info_mode).get(_FLAG, False)
    if allow_float and isinstance(value, float) and value in (0.0, 1.0):
        return bool(value)        # accepted by UHFk before the switch was shared
    if isinstance(value, str):
        logger.warning(
            "%s should be a boolean, not the string %r; write "
            "enable_spin_orbital = true/false without quotes.", name, value)
    return strict_bool(value, name)


def physical_norb(geom_norb, enable_spin_orbital, *, check_norb=None,
                  source="geom.dat"):
    """Physical orbital count from a geometry ``norb``.

    In spin-orbital mode ``geom.dat``'s ``norb`` is the spin-orbital count
    (= 2 * physical orbitals = Wannier90 num_wann), so halve it. Evenness is
    validated on ``check_norb`` (the *original*, pre-sublattice-fold value)
    so the error names the user's actual ``geom.dat`` entry, while the
    returned count is derived from ``geom_norb`` (which may be the post-fold
    value ``orig * subvol``)."""
    if not enable_spin_orbital:
        return geom_norb
    cn = check_norb if check_norb is not None else geom_norb
    if cn % 2 != 0:
        raise ValueError(
            "spin-orbital mode requires an even Geometry norb (the spin-orbital "
            "count = 2 * physical orbitals); got {} in {}".format(cn, source))
    return geom_norb // 2


def interleaved_index(orb, spin):
    """Interleaved spin-orbital index of ``(orb, spin)``."""
    return 2 * orb + spin


def spin_of(index):
    """Spin (0 = up, 1 = down) of an interleaved index (array-valued too)."""
    return index % 2


def to_spinblock_index(index, norb_phys):
    """Spin-block position of interleaved index ``index``."""
    return (index % 2) * norb_phys + index // 2


def spinblock_order(norb_phys):
    """Gather order ``g`` with ``g[j] =`` the interleaved index of spin-block
    position ``j``: ``M_spinblock = M_interleaved[g][:, g]``."""
    j = np.arange(2 * norb_phys)
    return interleaved_index(j % norb_phys, j // norb_phys)


def interleaved_to_spinblock(mat, norb_phys, axes=(-2, -1)):
    """Reorder the given axes of ``mat`` from interleaved to spin-block
    order (a copy; the identity for ``norb_phys == 1``)."""
    g = spinblock_order(norb_phys)
    out = np.asarray(mat)
    for ax in axes:
        out = np.take(out, g, axis=ax)
    return out


def spinblock_to_interleaved(mat, norb_phys, axes=(-2, -1)):
    """Inverse of :func:`interleaved_to_spinblock`."""
    inv = np.argsort(spinblock_order(norb_phys))
    out = np.asarray(mat)
    for ax in axes:
        out = np.take(out, inv, axis=ax)
    return out


def validate_chi0q_index_convention(data, enable_spin_orbital, file_name=""):
    """Reject a stored chi0q whose spin-orbital index convention is unknown.

    RPA writes ``index_convention="spin_block"`` (spin*norb+orb) into its
    chi0q.npz / chiq.npz. UHFk instead uses the interleaved (2*orb+spin)
    ordering, and chi0q files produced before the SO convention fix
    (commit 9dd9a21) carry no ``index_convention`` key at all. In spin-orbital
    mode the two orderings differ, so silently accepting such a file would mix
    conventions and corrupt the result. Outside spin-orbital mode the
    convention is irrelevant and this is a no-op.

    Parameters
    ----------
    data : Mapping
        Loaded npz (e.g. ``numpy.load(...)``); must support ``in`` and ``[]``.
    enable_spin_orbital : bool
        Whether the consuming calculation runs in spin-orbital mode.
    file_name : str, optional
        Source path, used only for the error message.

    Raises
    ------
    ValueError
        If spin-orbital mode is on and the file lacks a ``spin_block``
        ``index_convention`` marker.
    """
    if not enable_spin_orbital:
        return
    if "index_convention" not in data:
        raise ValueError(
            "chi0q file '{}' lacks an 'index_convention' marker. It predates "
            "the spin-orbital index convention fix and may use the interleaved "
            "(2*orb+spin) ordering; regenerate it with the current RPA solver "
            "(which writes index_convention='spin_block').".format(file_name)
        )
    conv = str(data["index_convention"])
    if conv != INDEX_CONVENTION:
        raise ValueError(
            "chi0q file '{}' has index_convention='{}', but spin-orbital mode "
            "requires 'spin_block' (spin*norb+orb). Regenerate it with the "
            "current RPA solver.".format(file_name, conv)
        )
