"""Readers of the UHFr and UHFk outputs.

* :func:`read_energy`: ``energy.dat`` (:class:`Energy`).
"""
from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Optional, Union

from hwave.results._frozen import _freeze
from hwave.results._io import read_key_value_text
from hwave.results.errors import UnsupportedLayoutError
from hwave.results.run import as_run

__all__ = ["Energy", "read_energy"]

#: the ``energy.dat`` keys of the scalar fields of :class:`Energy`, in lower
#: case: UHFr writes ``Energy_total`` (uhfr.py:1271-1292), UHFk
#: ``Energy_Total`` (uhfk.py:1796-1835)
_ENERGY_FIELDS = (
    ("total", "energy_total"),
    ("band", "energy_band"),
    ("free_total", "freeenergy_total"),
    ("free_band", "freeenergy_band"),
    ("ncond", "ncond"),
    ("sz", "sz"),
    ("sx", "sx"),
    ("sy", "sy"),
)

_MU = "chemicalpotential"

#: the order of :attr:`Energy.meta` ``["missing_keys"]``
_FIELD_ORDER = ("total", "band", "free_total", "free_band", "ncond", "sz", "mu",
                "sx", "sy")


@dataclasses.dataclass(frozen=True)
class Energy:
    """The content of ``energy.dat``.

    Attributes
    ----------
    values : dict[str, float]
        Every line of the file, ``{key: value}`` in file order, the keys
        stripped of their padding (e.g. ``"Energy_CoulombIntra"``).
    total, band : float or None
        ``Energy_Total``/``Energy_total`` and ``Energy_Band``/``Energy_band``.
    free_total, free_band : float or None
        The free energy ``F = E - TS`` (``FreeEnergy_*``).
    ncond, sz : float or None
        ``NCond`` and ``Sz``.
    sx, sy : float or None
        ``Sx`` and ``Sy``, the last two lines of the file; ``None`` in files
        written before they were added.
    mu : float, dict[str, float] or None
        ``ChemicalPotential``; a dict ``{suffix: value}`` when the file has
        one ``ChemicalPotential_<suffix>`` per group (UHFk mu groups ``0``,
        ``1``, ...; UHFr blocks ``spin-up``, ``spin-down``).
    solver : str
        ``"UHFr"`` or ``"UHFk"``, from the input.
    source : Path
        The file read.
    meta : dict
        ``missing_keys``: the names of the fields above whose key is absent
        from the file, in field order.

    The keys are matched case-insensitively (the solvers write them in
    different case).
    """
    values: dict
    total: Optional[float]
    band: Optional[float]
    free_total: Optional[float]
    free_band: Optional[float]
    ncond: Optional[float]
    sz: Optional[float]
    sx: Optional[float]
    sy: Optional[float]
    mu: Union[float, dict, None]
    solver: str
    source: Path
    meta: dict


def read_energy(run_or_toml, *, base_dir=None) -> Energy:
    """Read ``energy.dat`` of a run (a :class:`hwave.results.Run` or the
    path of its input TOML file).

    Raises MissingOutputError when the output is not configured or its file
    is absent, and UnsupportedLayoutError for a line that is not
    ``key = number``, for a key given twice (compared case-insensitively),
    or for a file with both ``ChemicalPotential`` and
    ``ChemicalPotential_<suffix>``."""
    run = as_run(run_or_toml, base_dir=base_dir)
    path = run.path("energy")
    values = {}
    seen = {}
    for key, text, lineno in read_key_value_text(path):
        lower = key.lower()
        if lower in seen:
            raise UnsupportedLayoutError(
                "{}: {!r}: line {} repeats the key of line {}".format(
                    path, key, lineno, seen[lower]), path=path)
        seen[lower] = lineno
        try:
            values[key] = float(text)
        except ValueError:
            raise UnsupportedLayoutError(
                "{}: {!r}: line {} has the value {!r}, not a number".format(
                    path, key, lineno, text), path=path) from None

    by_lower = {key.lower(): value for key, value in values.items()}
    fields = {name: by_lower.get(lower) for name, lower in _ENERGY_FIELDS}
    groups = {key[len(_MU) + 1:]: value for key, value in values.items()
              if key.lower().startswith(_MU + "_")}
    if _MU in by_lower and groups:
        first = next(key for key in values if key.lower().startswith(_MU + "_"))
        raise UnsupportedLayoutError(
            "{}: {!r}: the file has both ChemicalPotential and per-group "
            "ChemicalPotential_<suffix> lines".format(path, first), path=path)
    mu = by_lower.get(_MU, groups or None)
    fields["mu"] = mu
    missing = tuple(name for name in _FIELD_ORDER if fields[name] is None)
    return Energy(values=_freeze(values), mu=_freeze(mu), solver=run.mode,
                  source=Path(path), meta={"missing_keys": missing},
                  **{name: fields[name] for name, _ in _ENERGY_FIELDS})
