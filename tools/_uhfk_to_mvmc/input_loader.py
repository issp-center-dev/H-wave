"""Input loaders for the bridge CLI."""
from __future__ import annotations

import sys
import numpy as np

try:
    import tomllib  # Python 3.11+
except ImportError:
    import tomli as tomllib  # Python 3.9-3.10


def load_input_toml(path):
    """Read input.toml and return a merged param dict.

    H-wave's UHFk solver keeps two related sections in input.toml:

    - ``[mode]`` — solver-level flags such as ``flag_fock`` and
      ``enable_spin_orbital`` (see ``hwave.solver.uhfk._init_mode``,
      which reads its ``param`` argument directly from the top-level
      ``[mode]`` section, not ``[mode.param]``).
    - ``[mode.param]`` — numeric parameters such as ``CellShape``,
      ``Ncond``, ``EPS``, ``BoundaryCondition``.

    The bridge needs both. This helper returns a shallow-merged dict:
    every top-level ``[mode]`` scalar (anything that is not the
    ``param`` subtable itself) is copied in first, then the
    ``[mode.param]`` entries are layered on top. ``[mode.param]``
    values therefore take precedence if a key is present in both.

    The merge lets the CLI look up ``enable_spin_orbital`` regardless
    of whether the user placed it under ``[mode]`` or ``[mode.param]``.
    """
    with open(path, "rb") as fp:
        data = tomllib.load(fp)
    if "mode" not in data or "param" not in data["mode"]:
        raise ValueError(f"{path}: missing [mode.param] section")
    merged = {}
    for key, value in data["mode"].items():
        if key == "param":
            continue
        merged[key] = value
    merged.update(data["mode"]["param"])
    return merged


def derive_ne_per_group(toml_param):
    """Return [Ne_up, Ne_down] from Ncond and 2Sz (Sz-fixed)."""
    ncond = int(toml_param["Ncond"])
    two_sz = int(toml_param.get("2Sz", 0))
    if (ncond + two_sz) % 2 != 0 or (ncond - two_sz) % 2 != 0:
        raise ValueError(
            f"Ncond={ncond}, 2Sz={two_sz} do not give integer "
            f"Ne_up = (Ncond + 2Sz)/2 / Ne_down = (Ncond - 2Sz)/2"
        )
    return [(ncond + two_sz) // 2, (ncond - two_sz) // 2]


def load_geometry_uhf(path):
    """Parse a geometry_uhf.dat file.

    Format (matching ``src/hwave/qlmsio/wan90.py::read_geometry``):
        rows 0..2 : 3 unit_vec lines (3 floats each)
        row  3    : 1 degree line (3 floats; unused)
        rows 4..6 : 3 cell_vec lines (3 ints each; unused)
        rows 7..  : site lines, each ``R_x R_y R_z orb_idx`` (ints)

    Returns
    -------
    unit_vec : (3, 3) float
        Lattice vectors (rows = a_d).
    site_R_int : (Ns, 3) float
        Integer cell indices ``R_i`` per site, returned as float for
        downstream einsum convenience. Sites with ``orb_idx != 0`` are
        skipped because this bridge path supports one physical orbital.
    norb : int
        Number of distinct orbitals per cell (count of site lines with
        ``R == (0, 0, 0)``); used by the CLI to enforce the
        ``norb_orig == 1`` guard.

    Note
    ----
    The bridge's plane-wave construction expects ``site_R_int`` to be
    integer lattice cell indices, not Cartesian positions, because
    ``fij_builder`` dots it against ``k = 2π * wavevector_index / L``
    (both dimensionless) — see ``tools/_uhfk_to_mvmc/fij_builder.py``
    docstring.
    """
    with open(path) as fp:
        lines = [ln.rstrip("\n") for ln in fp if ln.strip()]
    if len(lines) < 7:
        raise ValueError(
            f"{path}: geometry file has fewer than 7 header lines "
            f"(got {len(lines)})"
        )
    unit_vec = np.array(
        [[float(x) for x in lines[i].split()] for i in range(3)],
        dtype=np.float64,
    )
    site_R = []
    norb = 0
    for ln in lines[7:]:
        toks = ln.split()
        if len(toks) < 4:
            continue
        R = (int(toks[0]), int(toks[1]), int(toks[2]))
        orb = int(toks[3])
        if R == (0, 0, 0):
            norb += 1
        if orb != 0:
            continue  # This bridge path supports one physical orbital.
        site_R.append(R)
    if not site_R:
        raise ValueError(f"{path}: no site lines with orb_idx==0 found")
    return unit_vec, np.array(site_R, dtype=np.float64), norb
