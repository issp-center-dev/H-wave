"""Bridge-side all-unique-classes ``orbitalidxgen.def`` emitter.

Under ``enable_spin_orbital = true`` + ``SubShape > [1, 1, 1]``,
StdFace's ``orbitalidxgen.def`` over-groups pair classes because the
FermionHubbardGC generator classifies (site_i, spin_i, site_j, spin_j)
pairs by the ORIGINAL lattice translation invariance. SOC combined with
sublattice folding keeps only supercell (folded) translation invariance,
so on-site pairs at different sublattice slots have different physical F
values — StdFace's over-grouping merges them into a single class and the
bridge's class-consistency check fires with residuals of
order 1e-1.

The emitter bypasses the over-grouping by assigning every non-redundant
``(site_i, spin_i, site_j, spin_j)`` upper-triangle pair its own unique
class idx. That trivially satisfies class-consistency (one value per
class). Downstream mVMC optimization loses class-based parameter
reduction, but the initial WF from the bridge is correct.

Format matches the mVMC ``InOrbitalGeneral`` schema (readdef.c
GetInfoOrbitalGeneral): 5 header lines + ``2*Ns^2 - Ns`` mapping rows
(upper triangle only, sign always +1) + ``NOrbitalIdx = 2*Ns^2 - Ns``
optimize-flag rows.

See docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for the emission and
class-consistency contracts.
"""
from __future__ import annotations


def emit_orbitalidx_all_unique(
    nsite: int, out_path: str, complex_type: int = 1
) -> int:
    """Emit an ``orbitalidxgen.def`` with all-unique class indices.

    Assigns each upper-triangle ``(all_i, all_j)`` with
    ``all_i = site_i + spin_i * nsite < all_j`` a distinct class idx
    ``0, 1, ..., 2*nsite*(2*nsite-1)/2 - 1``. Every row has sign = +1.

    Parameters
    ----------
    nsite : int
        Number of physical sites (Ns).
    out_path : str
        Output path for the emitted ``orbitalidxgen.def``.
    complex_type : int
        ComplexType header value (default 1 — SOC F values are complex).

    Returns
    -------
    int
        The number of classes emitted (= ``nsite * (2*nsite - 1)``).

    Raises
    ------
    ValueError
        If ``nsite < 1`` or ``complex_type not in {0, 1}``.
    """
    if int(nsite) < 1:
        raise ValueError(f"nsite must be positive; got {nsite}")
    if int(complex_type) not in (0, 1):
        raise ValueError(
            f"complex_type must be 0 or 1; got {complex_type}"
        )
    nsite = int(nsite)
    complex_type = int(complex_type)
    n_all = 2 * nsite
    n_classes = nsite * (n_all - 1)  # 2*Ns*(2*Ns-1)/2 == Ns*(2*Ns-1)

    with open(out_path, "w") as fw:
        fw.write("=============================================\n")
        fw.write(f"NOrbitalIdx {n_classes:10d}\n")
        fw.write(f"ComplexType {complex_type:10d}\n")
        fw.write("=============================================\n")
        fw.write("=============================================\n")
        # Mapping rows: (i, spn_i, j, spn_j, class_idx, sign).
        # Iterate the upper triangle in mVMC spin-block order
        # all_i = i + spin_i * nsite < all_j; assign class_idx in
        # ascending order so mVMC's InOrbitalGeneral consumer sees one
        # class per pair.
        idx = 0
        for all_i in range(n_all):
            i_site, spin_i = all_i % nsite, all_i // nsite
            for all_j in range(all_i + 1, n_all):
                j_site, spin_j = all_j % nsite, all_j // nsite
                fw.write(
                    f"    {i_site:3d}  {spin_i:d}"
                    f"     {j_site:3d}  {spin_j:d}"
                    f"     {idx:6d}  {1:2d}\n"
                )
                idx += 1
        # Optimize flags: one per class, all enabled (matches StdFace
        # default). mVMC treats these as "optimize this parameter"
        # bits; they have no effect on the initial WF, which is what
        # the bridge is producing.
        for k in range(n_classes):
            fw.write(f"    {k:3d}      1\n")

    return n_classes
