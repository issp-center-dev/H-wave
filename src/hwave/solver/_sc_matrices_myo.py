"""Spin/charge interaction matrices for the full-vertex FLEX path.

Historically this module carried a hand-maintained sibling of
``hwave.sc._build_sc_matrices_all_q`` differing in exactly one element: the
charge (ab, ab) slot, where it used the MYO (cond-mat/0407094) value ``-U'+2J``
against the Kuroki (arXiv:0902.3691) value ``-U'+J`` used by ``sc.py``.

Exact diagonalization of H-wave's documented file operators adjudicated that
difference (issue #113): per interaction type the exact charge vertex at that
slot is ``+J`` from the Hund file term and ``+J`` from the Exchange file term.
For the SU(2) Kanamori combination (Hund + Exchange at equal J) their sum
reproduces the standard literature value ``-U'+2J`` -- so the MYO matrix was
right for the COMBINATION but wrong as a per-type attribution, assigning the
whole ``2J`` to the Hund file entry and nothing to Exchange. With the
per-type values fixed in ``sc.py``, the two builders are identical, and this
module now simply re-exports the single implementation.
"""

from hwave.sc import _build_sc_matrices_all_q


def build_sc_matrices_myo(inter_k, norb, Nx, Ny, Nz):
    """One implementation: see :func:`hwave.sc._build_sc_matrices_all_q`."""
    return _build_sc_matrices_all_q(inter_k, norb, Nx, Ny, Nz)
