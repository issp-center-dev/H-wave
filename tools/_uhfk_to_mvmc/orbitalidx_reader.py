"""Parser for mVMC orbitalidx.def (AntiParallel) input file.

Supports:
- 3-column form: ``i j orbital_idx`` (PBC; ``NMPTrans > 0``)
- 4-column form: ``i j orbital_idx sign`` (APBC; ``NMPTrans < 0``)
Followed by ``NOrbitalIdx`` optimize-flag lines: ``idx optimize``.

See docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for path selection and
accepted formats.
"""
from __future__ import annotations

import os


class OrbitalidxFormatError(ValueError):
    """Raised on malformed orbitalidx.def content."""


def parse_orbitalidx_def(path: str) -> dict:
    """Parse an orbitalidx.def file and return a structured dict.

    Returns
    -------
    dict with keys:
        n_orbital_idx : int — total unique parameter count.
        complex_type : int — 0 (real) or 1 (complex).
        mapping : dict[(int, int)] = (idx, sign)
            (i, j) → (param idx in [0, n_orbital_idx), sign in {-1, +1}).
            For 3-column form sign defaults to +1.
        optimize_flags : dict[int] = int
            idx → optimize flag (0 or 1).
        nsite : int — inferred from max(i, j) + 1.
        has_sign_column : bool — True if the body used 4-column form.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"orbitalidx.def not found: {path}")

    with open(path) as fp:
        lines = [ln.rstrip("\n") for ln in fp]

    n_orbital_idx = None
    complex_type = None
    data_start = 0
    for idx, ln in enumerate(lines):
        stripped = ln.strip()
        if stripped.startswith("=="):
            continue
        if not stripped:
            continue
        toks = stripped.split()
        if toks[0] == "NOrbitalIdx":
            n_orbital_idx = int(toks[1])
            continue
        if toks[0] == "ComplexType":
            complex_type = int(toks[1])
            continue
        data_start = idx
        break

    if n_orbital_idx is None:
        raise OrbitalidxFormatError(f"{path}: missing NOrbitalIdx header")
    if complex_type is None:
        raise OrbitalidxFormatError(f"{path}: missing ComplexType header")

    mapping = {}
    optimize_flags = {}
    has_sign_column = None
    max_site = -1

    for ln in lines[data_start:]:
        stripped = ln.strip()
        if not stripped or stripped.startswith("=="):
            continue
        toks = stripped.split()
        if len(toks) in (3, 4):
            i, j = int(toks[0]), int(toks[1])
            idx = int(toks[2])
            sign = int(toks[3]) if len(toks) == 4 else 1
            if len(toks) == 4:
                if has_sign_column is False:
                    raise OrbitalidxFormatError(
                        f"{path}: mixed 3-column and 4-column mapping lines"
                    )
                has_sign_column = True
            else:
                if has_sign_column is True:
                    raise OrbitalidxFormatError(
                        f"{path}: mixed 4-column and 3-column mapping lines"
                    )
                has_sign_column = False
            if not (-1 <= idx < n_orbital_idx):
                raise OrbitalidxFormatError(
                    f"{path}: orbital_idx {idx} out of range "
                    f"[-1, {n_orbital_idx})"
                )
            if sign not in (-1, 1):
                raise OrbitalidxFormatError(
                    f"{path}: sign {sign} not in {{-1, +1}}"
                )
            mapping[(i, j)] = (idx, sign)
            max_site = max(max_site, i, j)
        elif len(toks) == 2:
            idx = int(toks[0])
            opt = int(toks[1])
            if not (0 <= idx < n_orbital_idx):
                raise OrbitalidxFormatError(
                    f"{path}: optimize-flag idx {idx} out of range "
                    f"[0, {n_orbital_idx})"
                )
            if opt not in (0, 1):
                raise OrbitalidxFormatError(
                    f"{path}: optimize flag {opt} not in {{0, 1}}"
                )
            optimize_flags[idx] = opt
        else:
            raise OrbitalidxFormatError(
                f"{path}: unexpected line with {len(toks)} tokens: {ln!r}"
            )

    nsite = max_site + 1
    expected_pairs = nsite * nsite
    if len(mapping) != expected_pairs:
        raise OrbitalidxFormatError(
            f"{path}: got {len(mapping)} (i, j) mapping lines but "
            f"expected nsite^2 = {expected_pairs} (nsite={nsite})"
        )
    if len(optimize_flags) != n_orbital_idx:
        raise OrbitalidxFormatError(
            f"{path}: got {len(optimize_flags)} optimize-flag lines but "
            f"expected NOrbitalIdx = {n_orbital_idx}"
        )

    return {
        "n_orbital_idx": n_orbital_idx,
        "complex_type": complex_type,
        "mapping": mapping,
        "optimize_flags": optimize_flags,
        "nsite": nsite,
        "has_sign_column": bool(has_sign_column),
    }
