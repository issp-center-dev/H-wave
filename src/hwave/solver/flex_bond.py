"""The dynamic bond pipeline of the Phase B FLEX solver (GitHub issue
#181; spec docs/superpowers/specs/2026-09-06-flex-bond-sigma-phase-b-181-design.md,
section 3): the dense bond-block store, the bubble assembly, the batched
dressing, the bond-aware self-energy transport, the memory preflight and
the output helpers.
"""
import numpy as np

from . import backend as _bk
from . import bubble as _bubble


class BondBlockStore:
    """Owner of the dense ``(nmat, nvol, ND, ND)`` bond arrays (spec 3.5).

    A context manager: construction allocates every named array (rolling
    back on a partial failure), ``release()`` runs on every exit of the
    ``with`` block. Pair access (``put_pair`` / ``get_pair``: a view
    ``(nmat, nvol, nd, nd)`` onto channel pair ``(alpha, beta)``) and
    frequency-batch access (``put_freq_batch`` / ``get_freq_batch``: a
    contiguous ``(nb, nvol, ND, ND)`` view) are the only ways in and out;
    ``detach(name)`` transfers ownership of one array to the caller.
    A future streaming store implements the same five methods.
    """

    def __init__(self, nmat, nvol, ND, nd, names):
        names = tuple(names)
        if len(set(names)) != len(names):
            raise ValueError("BondBlockStore: duplicate array names {}".format(names))
        if ND % nd != 0:
            raise ValueError("BondBlockStore: ND={} is not a multiple of nd={}".format(ND, nd))
        self.nmat, self.nvol, self.ND, self.nd = int(nmat), int(nvol), int(ND), int(nd)
        self.B = self.ND // self.nd
        self._arrays = {}
        try:
            for name in names:
                self._arrays[name] = np.zeros((self.nmat, self.nvol, self.ND, self.ND),
                                              dtype=np.complex128)
        except MemoryError:
            self._arrays.clear()
            raise
        self._released = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        self.release()
        return False

    def _arr(self, name):
        if self._released:
            raise RuntimeError("BondBlockStore: released")
        return self._arrays[name]

    def put_pair(self, name, alpha, beta, block):
        a = self._arr(name)
        nd = self.nd
        a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd] = block

    def get_pair(self, name, alpha, beta):
        a = self._arr(name)
        nd = self.nd
        return a[:, :, alpha * nd:(alpha + 1) * nd, beta * nd:(beta + 1) * nd]

    def put_freq_batch(self, name, l0, l1, batch):
        self._arr(name)[l0:l1] = batch

    def get_freq_batch(self, name, l0, l1):
        return self._arr(name)[l0:l1]

    def detach(self, name):
        a = self._arr(name)
        del self._arrays[name]
        return a

    def release(self):
        self._arrays.clear()
        self._released = True


def assemble_bubble(store, green_scf, green0_tail, beta, view, spatial_shape, workers):
    """Fill ``store['chibar']`` pair by pair from ``bubble._iter_bond_dynamic``
    (spec 3.1). ``green_scf`` is the TAIL-SUBTRACTED single-block Green
    function the general bubble consumes; ``green0_tail`` its tail."""
    for (alpha, beta_), block in _bubble._iter_bond_dynamic(
            green_scf, green0_tail, beta, view, spatial_shape=tuple(spatial_shape),
            workers=workers):
        store.put_pair("chibar", alpha, beta_, _bk.to_host(block))
        del block
