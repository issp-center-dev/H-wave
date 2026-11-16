"""The spin-orbital bond-resolved FLEX (issue #270): the admission rows of the
FLEX SCF beyond the bond plan's own (:func:`flex_so_bond_rows`, for
:func:`hwave.solver.so_bond.plan_so_bond` with the consumer ``flex_scf``).
"""
from . import flex_bond

__all__ = ["flex_so_bond_rows"]

_ENTRY = 16                     # bytes per complex128 entry
#: the per-map Green-function digest hashes at most this many bytes per copy
_DIGEST_COPY_CAP = 64 * 2 ** 20


def flex_so_bond_rows(*, nmat, nvol, nd, depth, mixing, n_types, split_seed=False):
    """The FLEX SCF rows of a spin-orbital bond plan: ``(resident_rows,
    extra_rows)`` for :func:`hwave.solver.so_bond.plan_so_bond`
    (``resident_rows=`` and ``extra_rows=``).

    With ``G = 16 nmat nvol nd^2`` (one spin-orbital Green function) and
    ``C = 16 nmat nvol nd^4`` (one channel-0 susceptibility block over the
    frequency grid):

    * resident: the persistent rows of :func:`hwave.solver.flex_bond.
      flex_scf_rows` with ``nmodes = nd`` (``state``, ``seed_envelope``,
      ``anderson_history``, ``anderson_work``, ``eigenpairs_hf``,
      ``flex_arrays``, ``hf_tables``) and ``collapses = 2 C`` (the chi0q
      and chiq buffers of the last map);
    * extra: its phase rows (``green_mu``, ``density_hf``, ``convergence``,
      ``mixing``; each row named by its phase), ``("convergence",
      "green_digest", 2 min(G, 64 MiB))`` (the per-map Green-function
      digest), ``("transport", "transport", 4 G + 4 C)`` and
      ``("post_scf", "post_scf", G + C)``.

    ``split_seed`` is False on the spin-orbital path (its seeds are
    refused), so ``seed_envelope`` is 0 by default. Every byte count is an
    ``int``."""
    nmat, nvol, nd = int(nmat), int(nvol), int(nd)
    persistent, phase = flex_bond.flex_scf_rows(nmat=nmat, nvol=nvol, nmodes=nd, depth=depth,
                                                mixing=mixing, n_types=n_types,
                                                split_seed=split_seed)
    G = _ENTRY * nmat * nvol * nd * nd
    C = _ENTRY * nmat * nvol * nd ** 4
    resident = tuple((k, int(v)) for k, v in persistent.items()) + (("collapses", int(2 * C)),)
    extra = tuple((k, k, int(v)) for k, v in phase.items()) + (
        ("convergence", "green_digest", int(2 * min(G, _DIGEST_COPY_CAP))),
        ("transport", "transport", int(4 * G + 4 * C)),
        ("post_scf", "post_scf", int(G + C)))
    return resident, extra
