"""Assert Transfer.dat rows satisfy
H[s,t](R) = conj(H[t,s](-R)) for every nontrivial pair.

Called from run.sh startup and on the single-case path, and from the
manifest producer, to guard against a non-Hermitian fixture Transfer.dat.
See docs/en/source/algorithm/uhfk_to_mvmc.rst for the Transfer mapping.
"""
from __future__ import annotations

import math
import os

from tools._uhfk_to_mvmc.trans_emit import (
    parse_hwave_transfer,
    TransEmitError,
)


class TransferHermiticityError(RuntimeError):
    """Raised when Transfer.dat Hermiticity validation fails."""


def validate_transfer_entries_hermiticity(
    entries, *, source="Transfer.dat", atol=1e-12,
):
    """Validate parsed Transfer.dat entries for finite Hermitian pairs.

    ``entries`` must contain the post-``ndegen`` values returned by
    :func:`parse_hwave_transfer`. Empty entries are not themselves a
    Hermiticity violation, and declared inventory metadata is unavailable at
    this layer; callers apply their separate product or fixture policies.
    """
    # Validate every coefficient before building or comparing partner values:
    # NaN makes ordered comparisons false, and infinities are never valid
    # hopping coefficients even when paired symmetrically.
    for row_num, (dx, dy, dz, s1, s2, value) in enumerate(entries, start=1):
        if not (math.isfinite(value.real) and math.isfinite(value.imag)):
            raise TransferHermiticityError(
                f"{source}: non-finite coefficient in row {row_num} "
                f"(dR=({dx},{dy},{dz}), s1={s1}, s2={s2}, "
                f"re={value.real}, im={value.imag})"
            )

    # Index by (dx, dy, dz, s1, s2).
    index = {(r[0], r[1], r[2], r[3], r[4]): r[5] for r in entries}
    for dx, dy, dz, s1, s2, value in entries:
        re, im = value.real, value.imag
        if (dx, dy, dz) == (0, 0, 0) and s1 == s2 and abs(im) < atol:
            # A real same-spin on-site term is its own conjugate. Real
            # off-site hopping is not: it still requires the -R row.
            continue
        partner_key = (-dx, -dy, -dz, s2, s1)
        if partner_key not in index:
            raise TransferHermiticityError(
                f"{source}: missing partner for row "
                f"(dR=({dx},{dy},{dz}), s1={s1}, s2={s2}, "
                f"re={re}, im={im}); expected "
                f"(dR=({-dx},{-dy},{-dz}), s1={s2}, s2={s1})"
            )
        partner = index[partner_key]
        partner_re, partner_im = partner.real, partner.imag
        # Hermitian invariant: H[s,t](R) = conj(H[t,s](-R))
        # -> re(H[s,t](R)) = re(H[t,s](-R))
        # -> im(H[s,t](R)) = -im(H[t,s](-R))
        if abs(re - partner_re) > atol:
            raise TransferHermiticityError(
                f"{source}: wrong magnitude (real) for row "
                f"(dR=({dx},{dy},{dz}), s1={s1}, s2={s2}): "
                f"expected re={partner_re}, got {re}"
            )
        if abs(im + partner_im) > atol:
            raise TransferHermiticityError(
                f"{source}: wrong conj sign (imag) for row "
                f"(dR=({dx},{dy},{dz}), s1={s1}, s2={s2}): "
                f"expected im={-partner_im}, got {im}"
            )


def _read_declared_row_count(path):
    """Return the already parser-validated ``nr`` declaration."""
    try:
        with open(path) as fp:
            for _ in range(2):
                next(fp)
            return int(next(fp).strip())
    except (OSError, StopIteration, ValueError) as exc:
        # parse_hwave_transfer validates this header before this helper runs.
        # Reaching here means the file changed between reads, so fail closed.
        raise TransferHermiticityError(
            f"{path}: could not re-read declared nr after parsing"
        ) from exc


def check_transfer_dat_hermiticity(path, atol=1e-12):
    """Verify Transfer.dat satisfies H[s,t](R) = conj(H[t,s](-R))
    for every row except a real, same-spin on-site term, which is
    self-conjugate and therefore passes without a distinct partner.

    Raises TransferHermiticityError with a descriptive message on
    the first violation encountered, or if the file is empty,
    headerless, or malformed. Raises FileNotFoundError if path does
    not exist.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Transfer.dat not found: {path}")
    try:
        rows = parse_hwave_transfer(path)
    except TransEmitError as exc:
        raise TransferHermiticityError(str(exc)) from exc

    # Keep this inventory gate here rather than adding expected_rows to every
    # caller. The canonical parser intentionally permits sparse all-unit files
    # and returns no metadata, so the checker reads only its already-validated
    # nr declaration and still delegates all data-row parsing to that parser.
    # Deriving the expectation from the file protects every call site and
    # avoids duplicating fixture-specific counts (currently 24 and 16).
    declared_rows = _read_declared_row_count(path)
    if not rows:
        raise TransferHermiticityError(f"{path}: no data rows")

    validate_transfer_entries_hermiticity(rows, source=path, atol=atol)

    # Check the declared inventory after the pair diagnostics so an
    # asymmetric truncation keeps the more specific missing-partner error;
    # a symmetrically truncated Hermitian subset reaches and fails this gate.
    if len(rows) != declared_rows:
        raise TransferHermiticityError(
            f"{path}: declared nr={declared_rows} data rows, found {len(rows)}"
        )
