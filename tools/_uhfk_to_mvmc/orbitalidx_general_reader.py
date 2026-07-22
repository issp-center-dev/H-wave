"""Parser for mVMC orbitalidx_general.def (General 6-column format).

mVMC readdef.c GetInfoOrbitalGeneral (readdef.c:2949+) reads

    <i> <spn_i> <j> <spn_j> <fij_class_idx> <fijSign>

per mapping row, followed by NOrbitalIdx optimize-flag lines. The
6-column arity distinguishes InOrbitalGeneral from the
InOrbitalAntiParallel 3/4-column format handled by
``orbitalidx_reader.parse_orbitalidx_def``.
"""
from __future__ import annotations

import os


class OrbitalidxFormatError(ValueError):
    """Raised on malformed orbitalidx_general.def content or unexpected
    mapping-line arity in detect_orbitalidx_format."""


def _iter_non_separator_lines(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"orbitalidx file not found: {path}")
    with open(path) as fp:
        for ln in fp:
            stripped = ln.strip()
            if not stripped or stripped.startswith("=="):
                continue
            yield stripped


def detect_orbitalidx_format(path: str) -> str:
    """Inspect the first mapping-line arity and return format tag.

    Returns
    -------
    "antiparallel"
        First mapping line has 3 (PBC) or 4 (APBC) whitespace-
        separated integers.
    "general"
        First mapping line has 6 integers (InOrbitalGeneral).

    Raises
    ------
    OrbitalidxFormatError
        Header missing or first mapping line arity is anything other than
        3, 4, or 6.
    """
    n_orbital_idx = None
    complex_type = None
    for stripped in _iter_non_separator_lines(path):
        toks = stripped.split()
        if toks[0] == "NOrbitalIdx":
            n_orbital_idx = int(toks[1])
            continue
        if toks[0] == "ComplexType":
            complex_type = int(toks[1])
            continue
        # First real mapping line (integers only)
        try:
            [int(t) for t in toks]
        except ValueError:
            raise OrbitalidxFormatError(
                f"orbitalidx file has non-integer mapping token: {stripped!r}"
            )
        arity = len(toks)
        if arity in (3, 4):
            return "antiparallel"
        if arity == 6:
            return "general"
        raise OrbitalidxFormatError(
            f"unexpected mapping-line arity {arity} in {path}; expected "
            f"3 or 4 (AntiParallel) or 6 (General). "
            f"Sample line: {stripped!r}"
        )
    raise OrbitalidxFormatError(
        f"orbitalidx file {path} has no mapping lines after header"
    )


def parse_orbitalidx_general_def(path: str) -> dict:
    """Parse an orbitalidx_general.def (6-column mVMC format).

    Returns
    -------
    dict with keys:
        n_orbital_idx : int — total unique parameter count from header.
        complex_type : int — must be 1 because F is complex.
        mapping : dict[(all_i, all_j)] = (idx, sign)
            all_i = i + spn_i * Nsite, all_j = j + spn_j * Nsite,
            all_i < all_j required (upper triangle).
        optimize_flags : dict[int] = int
            idx → optimize flag (0 or 1).
        nsite : int — inferred from max((i, j)) + 1.
        format : "general"
        has_sign_column : True (always for 6-column general form).

    Raises
    ------
    OrbitalidxFormatError
        Header missing/invalid; mapping row arity != 6; sign not in
        {-1, +1}; ``all_i >= all_j``; total rows != ``2 * Ns**2 - Ns``.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"orbitalidx_general.def not found: {path}")
    with open(path) as fp:
        lines = [ln.rstrip("\n") for ln in fp]

    n_orbital_idx = None
    complex_type = None
    data_start = None
    for idx, ln in enumerate(lines):
        stripped = ln.strip()
        if not stripped or stripped.startswith("=="):
            continue
        toks = stripped.split()
        if toks[0] == "NOrbitalIdx":
            n_orbital_idx = int(toks[1])
        elif toks[0] == "ComplexType":
            complex_type = int(toks[1])
        else:
            data_start = idx
            break
    if n_orbital_idx is None or complex_type is None:
        raise OrbitalidxFormatError(
            f"orbitalidx_general.def missing NOrbitalIdx or ComplexType header: {path}"
        )
    if data_start is None:
        raise OrbitalidxFormatError(
            f"orbitalidx_general.def has no mapping rows after header: {path}"
        )

    mapping = {}
    optimize_flags = {}
    all_i_max = -1
    n_mapping = 0
    for ln in lines[data_start:]:
        stripped = ln.strip()
        if not stripped or stripped.startswith("=="):
            continue
        toks = stripped.split()
        if len(toks) == 6:
            i, spn_i, j, spn_j, fij_idx, fij_sign = (int(t) for t in toks)
            if fij_sign not in (-1, 1):
                raise OrbitalidxFormatError(
                    f"sign column must be -1 or +1; got {fij_sign} in row {stripped!r}"
                )
            all_i_max = max(all_i_max, i, j)
            n_mapping += 1
            # Nsite will be determined after we know nsite; store raw for now
            mapping[(i, spn_i, j, spn_j)] = (fij_idx, fij_sign)
        elif len(toks) == 2:
            k, flag = int(toks[0]), int(toks[1])
            optimize_flags[k] = flag
        else:
            raise OrbitalidxFormatError(
                f"unexpected mapping row arity {len(toks)}: {stripped!r}"
            )

    nsite = all_i_max + 1

    # Rebuild mapping keyed by (all_i, all_j) with all_i < all_j check.
    mapping_all = {}
    for (i, spn_i, j, spn_j), (fij_idx, fij_sign) in mapping.items():
        all_i = i + spn_i * nsite
        all_j = j + spn_j * nsite
        if all_i >= all_j:
            raise OrbitalidxFormatError(
                f"all_i={all_i} >= all_j={all_j} violates upper-triangle "
                f"constraint (nsite={nsite}, row was "
                f"i={i} spn_i={spn_i} j={j} spn_j={spn_j})"
            )
        mapping_all[(all_i, all_j)] = (fij_idx, fij_sign)

    expected_rows = 2 * nsite * nsite - nsite
    if n_mapping != expected_rows:
        raise OrbitalidxFormatError(
            f"orbitalidx_general.def has {n_mapping} mapping rows; "
            f"expected 2*Ns^2 - Ns = {expected_rows} for Nsite={nsite}"
        )

    return {
        "n_orbital_idx": n_orbital_idx,
        "complex_type": complex_type,
        "mapping": mapping_all,
        "optimize_flags": optimize_flags,
        "nsite": nsite,
        "format": "general",
        "has_sign_column": True,
    }
