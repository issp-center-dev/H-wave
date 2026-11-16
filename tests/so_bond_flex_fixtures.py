"""Public-entry fixtures of the spin-orbital bond-resolved FLEX (issue #270).

The fixture tables extend :mod:`tests.so_bond_entry_fixtures` (whose
conventions apply: interleaved Transfer index, physical-orbital interaction
and Extern indices, a 4x4x1 mesh with one physical orbital):

* :data:`ZEEMAN_SZ_ALL` -- the Zeeman band of ``ZEEMAN_EXTERN`` (an Extern
  field) with every S_z-conserving interaction type, a complex off-site
  PairHop included: the S_z-packed layout;
* :data:`SOC_ALL_TYPES_CPH` -- :data:`SOC_ALL_TYPES` with a complex off-site
  PairHop: the dense layout;
* :data:`ZEEMAN_EXTERN_ONSITE` -- ``ZEEMAN_EXTERN`` with its on-site rows only
  (``B = 1``).

:func:`write_flex_case` writes the files of a fixture and returns the
``input_dict`` of a FLEX run at a fixed chemical potential;
:func:`build_solver` constructs the solver from it as ``hwave.qlms.run``
does. At module level this module imports numpy, the standard library and
the test-side fixture module only; a helper that drives production code
imports ``hwave`` inside its own body.
"""
import cmath

import numpy as np

from tests import so_bond_entry_fixtures as F
from tests.so_bond_entry_fixtures import _nn, _override

__all__ = ["V_CPH", "ZEEMAN_SZ_ALL", "SOC_ALL_TYPES_CPH", "ZEEMAN_EXTERN_ONSITE",
           "FLEX_NMAT", "FLEX_T", "onsite_rows", "interaction_tables", "fixture_mu",
           "write_flex_case", "build_solver"]


def _nn_c(v):
    """A nearest-neighbour row on orbital 0 with a complex coefficient: ``v``
    on +x and +y, ``conj(v)`` on -x and -y (the Hermitian closure of a
    complex PairHop; ``_nn(v)`` with a complex ``v`` is not Hermitian-closed
    and is refused by the compiler)."""
    c = complex(v).conjugate()
    return (((1, 0, 0), (0, 0), v), ((-1, 0, 0), (0, 0), c),
            ((0, 1, 0), (0, 0), v), ((0, -1, 0), (0, 0), c))


#: the complex PairHop coefficient of the complex-PairHop fixtures
V_CPH = 0.43 * cmath.exp(0.7j)

#: the Zeeman band of ZEEMAN_EXTERN with every S_z-conserving interaction type
#: (a complex PairHop included; no PairLift): S_z conserved, the packed layout
ZEEMAN_SZ_ALL = F.ZEEMAN_EXTERN._replace(name="ZEEMAN_SZ_ALL", interactions=(
    ("CoulombIntra", (((0, 0, 0), (0, 0), 0.83),)),
    ("CoulombInter", _nn(0.37)),
    ("Hund", _nn(0.53)),
    ("Ising", _nn(0.11)),
    ("Exchange", _nn(0.31)),
    ("PairHop", _nn_c(V_CPH)),
))

#: SOC_ALL_TYPES with its PairHop rows replaced by the complex ones
SOC_ALL_TYPES_CPH = F.SOC_ALL_TYPES._replace(name="SOC_ALL_TYPES_CPH", interactions=tuple(
    (t, _nn_c(V_CPH) if t == "PairHop" else rows) for t, rows in F.SOC_ALL_TYPES.interactions))


def onsite_rows(fixture, name):
    """``fixture`` with its R = 0 interaction rows only (B = 1)."""
    return fixture._replace(name=name, interactions=tuple(
        (t, tuple(r for r in rows if r[0] == (0, 0, 0)))
        for t, rows in fixture.interactions if any(r[0] == (0, 0, 0) for r in rows)))


#: the band and field of ZEEMAN_EXTERN with its on-site interaction rows only
ZEEMAN_EXTERN_ONSITE = onsite_rows(F.ZEEMAN_EXTERN, "ZEEMAN_EXTERN_ONSITE")


def interaction_tables(fixture):
    """The interaction tables of ``fixture`` as the compiler takes them:
    ``{type: {((Rx, Ry, Rz), (a, b)): v}}`` (physical orbitals, 0-based)."""
    return {t: {(tuple(R), tuple(ab)): v for R, ab, v in rows}
            for t, rows in fixture.interactions}


#: the Matsubara count and the temperature of :func:`write_flex_case`
FLEX_NMAT = 16
FLEX_T = 0.25


def fixture_mu(fixture):
    """The chemical potential of ``fixture`` at its filling, from the literal
    band (:func:`tests.so_bond_entry_fixtures.independent_green` at
    ``FLEX_NMAT``, ``beta = 1 / FLEX_T``)."""
    return float(F.independent_green(fixture, nmat=FLEX_NMAT, beta=1.0 / FLEX_T)[1])


def write_flex_case(tmpdir, fixture, *, gate=True, so=True, param=None, mode=None):
    """Write the input files of ``fixture`` into ``tmpdir`` and return the
    ``input_dict`` of a FLEX run.

    The files are those of :func:`tests.so_bond_entry_fixtures.write_case`
    (``so = False`` writes the non-spin-orbital files of the spin-independent
    hopping). The ``[eliashberg]`` table and ``filling`` are removed and the
    run is set to: ``[mode] mode = "FLEX"``, ``enable_spin_orbital = so``,
    ``calc_scheme = "general"``; ``[mode.param] flex_hartree_fock = true``,
    ``longitudinal_bond_channels = gate``, ``Nmat = FLEX_NMAT``, ``T =
    FLEX_T``, ``mu`` = :func:`fixture_mu` (a fixed chemical potential: no
    ``filling`` / ``Ncond``), ``coeff_tail = 0.0``, ``IterationMax = 1``,
    ``EPS = 12``, ``Mix = 0.5``. ``param`` and ``mode`` override the
    ``[mode.param]`` and ``[mode]`` tables (keys matched
    case-insensitively; a key already present keeps its spelling)."""
    if mode is not None and any(k.lower() == "param" for k in mode):
        raise ValueError("write_flex_case: pass the [mode.param] overrides as param=")
    input_dict = F.write_case(tmpdir, fixture, mode={"enable_spin_orbital": bool(so)})
    input_dict.pop("eliashberg", None)
    par = {k: v for k, v in input_dict["mode"]["param"].items() if k.lower() != "filling"}
    _override(par, {"flex_hartree_fock": True, "longitudinal_bond_channels": bool(gate),
                    "Nmat": FLEX_NMAT, "T": FLEX_T, "mu": fixture_mu(fixture),
                    "coeff_tail": 0.0, "IterationMax": 1, "EPS": 12, "Mix": 0.5},
              "[mode.param]")
    _override(par, param, "[mode.param]")
    mode_tbl = _override({"mode": "FLEX", "enable_spin_orbital": bool(so),
                          "calc_scheme": "general"}, mode, "[mode]")
    mode_tbl["param"] = par
    input_dict["mode"] = mode_tbl
    return input_dict


def build_solver(input_dict):
    """``(solver, green_info)``: the FLEX solver of ``input_dict`` and the
    reader's Green-function entries, constructed as ``hwave.qlms.run`` does
    (the ``[mode]`` table checked by ``check_mode_tables``, ``[mode.param]``
    wrapped in a ``ParamDict``)."""
    import hwave.qlmsio.read_input_k as read_input_k
    import hwave.solver.flex as flex_mod
    from hwave.qlms import check_mode_tables
    from hwave.solver.base import ParamDict
    info_mode = dict(input_dict["mode"])
    errors = check_mode_tables(info_mode)
    if errors:
        raise ValueError("build_solver: " + "; ".join(errors))
    info_mode["param"] = ParamDict(info_mode.get("param", {}) or {})
    reader = read_input_k.QLMSkInput(input_dict["file"]["input"])
    solver = flex_mod.FLEX(reader.get_param("ham"), {}, info_mode)
    return solver, reader.get_param("green")
