"""Energy comparison helper for the G3 gate.

Provides the canonical helper `energy_relative_delta` and the workspace
resolver `_resolve_g3_paths` used by `tests/validation/uhfk_mvmc_pairproduct/
compare.py --mode g3`.

`energy_relative_delta` computes `delta_rel = |E_mvmc - E_hwave| /
max(|E_hwave|, 1e-12)`. Both paths are explicit strings; workspace-level
resolution lives in `_resolve_g3_paths`.

`_resolve_g3_paths` reads only the canonical `${workspace}/mvmc/
zvo_out_selected.dat` produced by the run.sh normalization step
(`normalize_mvmc_output`). Raw `zvo_out_*.dat` selection is
NOT done here; that lives in run.sh so mVMC binary variants (multi-sample
outputs, `output/` vs flat layout) are handled outside G3's contract.

See docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for the validation gates.
"""
from __future__ import annotations

import math
import os
import re
from statistics import mean
from typing import List, Tuple


class EnergyCompareError(ValueError):
    """Raised for parseable but semantically invalid energy inputs."""


def _parse_hwave_energy(path: str) -> float:
    """Parse H-wave UHFk's ``energy.dat`` and return Energy_Total.

    Format (uhfk.py:_save_results writes one ``KEY = VALUE`` per line):
        Energy_Total = -3.75000000000011
        Energy_Band  = ...
    Raises FileNotFoundError if the path does not exist, EnergyCompareError
    if the file exists but no Energy_Total line is present.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"hwave energy file not found: {path}")
    pattern = re.compile(
        r"\s*Energy_Total\s*=\s*(-?\d+\.\d+(?:[eE][+-]?\d+)?)"
    )
    with open(path) as fp:
        for line in fp:
            m = pattern.match(line)
            if m:
                return float(m.group(1))
    raise EnergyCompareError(f"Energy_Total not found in {path}")


def _parse_mvmc_zvo_out(path: str) -> List[float]:
    """Return the list of per-bin ``<H>`` values from a zvo_out_*.dat file.

    Format (mVMC 1.4.0): each non-comment line contains exactly six
    numeric columns beginning with ``<H> <H^2> <doublon> <singleon>``.
    Column 0 is the total-energy expectation for that bin. Raises
    FileNotFoundError if the path does not exist, EnergyCompareError if
    any row violates the schema or no numeric rows are present.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"mvmc zvo_out file not found: {path}")
    samples: List[float] = []
    with open(path) as fp:
        for lineno, line in enumerate(fp, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            toks = line.split()
            if len(toks) != 6:
                raise EnergyCompareError(
                    f"{path}:{lineno}: expected exactly 6 columns; "
                    f"got {len(toks)}"
                )
            try:
                values = [float(token) for token in toks]
            except ValueError as exc:
                raise EnergyCompareError(
                    f"{path}:{lineno}: expected numeric columns ({exc})"
                ) from exc
            if not all(math.isfinite(value) for value in values):
                raise EnergyCompareError(
                    f"{path}:{lineno}: row contains non-finite values"
                )
            samples.append(values[0])
    if not samples:
        raise EnergyCompareError(f"no numeric rows parsed from {path}")
    return samples


def energy_relative_delta(
    hwave_energy_dat: str, mvmc_zvo_out_dat: str
) -> Tuple[float, float, float]:
    """Return ``(E_hwave, E_mvmc, delta_rel)`` from the two file paths.

    ``delta_rel = |E_mvmc - E_hwave| / max(|E_hwave|, 1e-12)`` is the
    relative energy delta used by the G3 gate. See
    docs/en/source/uhfk/tools/uhfk_to_mvmc.rst.

    `E_hwave` is `Energy_Total` from H-wave's `energy.dat`. `E_mvmc` is the
    mean of `<H>` samples from mVMC's `zvo_out_*.dat` (column 0). Both
    parsers are the canonical readers that MUST be used by G3; alternate
    parsers are out of contract.

    Parameters
    ----------
    hwave_energy_dat : str
        Path to H-wave's ``energy.dat``.
    mvmc_zvo_out_dat : str
        Path to mVMC's ``zvo_out_XXX.dat`` (or the normalized
        ``zvo_out_selected.dat``; the canonical G3 workspace resolver
        selects the latter).

    Returns
    -------
    (E_hwave, E_mvmc, delta_rel) : tuple[float, float, float]
    """
    e_hwave = _parse_hwave_energy(hwave_energy_dat)
    bins = _parse_mvmc_zvo_out(mvmc_zvo_out_dat)
    e_mvmc = mean(bins)
    delta_rel = abs(e_mvmc - e_hwave) / max(abs(e_hwave), 1e-12)
    return e_hwave, e_mvmc, delta_rel


def _resolve_g3_paths(workspace: str) -> Tuple[str, str]:
    """Return ``(hwave_energy_path, mvmc_zvo_out_selected_path)`` for G3.

    The caller (`run.sh`) MUST have
    normalized mVMC output to `${workspace}/mvmc/zvo_out_selected.dat`
    before invoking G3. This resolver never scans raw `zvo_out_*.dat`;
    the multi-file / layout-variant selection is `normalize_mvmc_output`'s
    responsibility in run.sh.

    See docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for the G3 contract.

    Raises FileNotFoundError with a documented message if either canonical
    path is missing. The message names the file so callers can distinguish
    "H-wave missed" from "normalization did not run".
    """
    hwave_energy = os.path.join(workspace, "hwave", "energy.dat")
    if not os.path.isfile(hwave_energy):
        raise FileNotFoundError(f"G3: missing {hwave_energy}")
    mvmc_selected = os.path.join(workspace, "mvmc", "zvo_out_selected.dat")
    if not os.path.isfile(mvmc_selected):
        raise FileNotFoundError(
            f"G3: missing {mvmc_selected}; run.sh must normalize the mVMC output"
        )
    return hwave_energy, mvmc_selected
