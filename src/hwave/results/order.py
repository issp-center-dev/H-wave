"""Order parameters of a mean-field density in one magnetic cell.

Input
-----
* ``rho``, shape ``(n_orb, 2, 2)`` complex: ``rho[o, s, t] = <c^dag_{o s}
  c_{o t}>`` for the orbital entries ``o`` of one magnetic cell, with
  ``s, t = 0 (up), 1 (down)``. This is the Green function slot as H-wave
  stores it, ``G_ab = <c^dag_a c_b>`` (uhfk.py ``UHFk._green``; uhfr.py),
  i.e. the TRANSPOSE of the operator density matrix ``rho_op`` with
  ``(rho_op)_{ts} = <c^dag_s c_t>``. A state with its spin along +y has
  ``rho[o] = [[1/2, i/2], [-i/2, 1/2]]``; reading the transpose instead
  flips the sign of S^y.
* ``site_of_orbital``, shape ``(n_orb,)``, int in ``[0, n_site)``.
* ``site_positions``, shape ``(n_site, 3)`` in fractional lattice
  coordinates, or ``None``.
* ``cell_shape`` ``B`` and ``lattice_shape`` ``L``, positive integer triples
  with ``L_a % B_a == 0``. Either may be ``None`` (no lattice): then only the
  per-site and the uniform (Q = 0) quantities are given.
* ``n_cell``, the number of magnetic cells in the run. With both shapes it
  must equal ``prod(L) / prod(B)`` and is derived from them when omitted;
  otherwise it defaults to 1.

Site defaults
-------------
When ``site_of_orbital`` is ``None`` and ``n_orb`` equals the number of cell
positions (the rows of ``site_positions``, else ``prod(cell_shape)``; one
orbital per cell position), ``site_of_orbital = arange(n_orb)``. With
several orbitals per cell position and no ``site_of_orbital``, a
``ValueError``: "site_of_orbital is required when a cell holds several
orbitals". There is no clustering and no tolerance-based grouping.

Site quantities
---------------
``rho_s = sum_{o: site(o) = s} rho[o]``, ``n_s = Tr rho_s`` and
``S_s = 1/2 sum_{t,t'} sigma_{t t'} rho_s[t, t']``, i.e.

* ``S^x_s = 1/2 Re(rho_s[0,1] + rho_s[1,0])``,
* ``S^y_s = 1/2 Re(-i (rho_s[0,1] - rho_s[1,0]))``,
* ``S^z_s = 1/2 (rho_s[0,0] - rho_s[1,1])``.

This is the formula of the solvers' ``Sx``, ``Sy``, ``Sz`` (uhfk.py
``_calc_phys``; uhfr.py ``_calc_phys``), in units of hbar as ``Sz``.

Uniform
-------
``M = N_site^-1 sum_s S_s``; ``N_site * n_cell * M`` equals the
``energy.dat`` totals ``(Sx, Sy, Sz)`` (within :func:`sum_rule_tolerance`,
since the saved Green function is the mixed one).

Fourier (exact formulas)
------------------------
* ``M(Q) = N_site^-1 sum_s S_s e^{i Q.r_s}``,
* ``dn(Q) = N_site^-1 sum_s (n_s - nbar) e^{i Q.r_s}``, with
  ``nbar = N_site^-1 sum_s n_s``,
* ``Q.r_s = 2 pi sum_a q_a r_{s,a}``, with ``r_s = site_positions[s]``
  (fractional) and ``q`` the fractional wave vector,
* ``|M(Q)| = sqrt(sum_a |M_a(Q)|^2)``.

Amplitudes are complex in general. With ``site_positions = None`` (or no
lattice), every Fourier, ordering and staggered field is ``None``.

Q set
-----
``q_a = m_a / B_a`` with integers ``m_a`` in ``(-B_a/2, B_a/2]``, stored as
``fractions.Fraction``, ``Q = 0`` included, in lexicographic order of
``(m_x, m_y, m_z)``. Charge amplitudes are reported for ``Q != 0`` (0 at
``Q = 0``).

Staggered
---------
``M_s = M(q_AF)`` with ``q_AF,a = 1/2`` on every axis where ``L_a > 1`` and
0 elsewhere. It is ``None`` unless ``B_a`` is even on all those axes (and
when no axis has ``L_a > 1``).

Ordering vector
---------------
The ``q`` with the largest amplitude. Ties within 1e-12 go to the smallest
``sum_a q_a^2`` (exact Fraction), then to the lexicographically smallest
``(q_x, q_y, q_z)``. ``None`` when the maximum is <= 1e-10.

Other relations
---------------
* ``cosine_amplitude(S, site_positions, q) = N_site^-1 sum_s cos(Q.r_s)
  S_s`` equals ``M(Q)`` only when every ``Q.r_s`` is in ``pi Z``.
* Spirals: for a single-Q coplanar spiral with ``Q != -Q`` (mod reciprocal
  lattice vectors), the local moment is ``sqrt(2) |M(Q)|``.
* ``mean_local_moment = N_site^-1 sum_s |S_s|``. It is the observable for
  patterns no commensurate ``Q`` captures, such as an antiferromagnet inside
  the cell.

Sum-rule tolerance
------------------
``sum_rule_tolerance(residual, *, mix, g_size, nvol, nd) =
(1 - mix) nvol sqrt(nd) residual g_size / 2 + nvol nd 1e-12 + 1e-12``:
the saved Green function differs from the unmixed one by
``(1 - mix) dG``, the residual is ``Rest = 2 ||dG||_F / g_size``, and
elements below 1e-12 are cut (base.py, the ``threshold``; uhfk.py and
uhfr.py after the mixing). For UHFr: ``g_size = 4 N^2``, ``nvol = 1``,
``nd = 2 N``.
"""
from __future__ import annotations

import dataclasses
import itertools
import math
import numbers
from fractions import Fraction
from pathlib import Path
from typing import Optional

import numpy as np

from hwave.results._frozen import _freeze

#: the tolerance of the Hermiticity check of ``rho``
HERMITIAN_ATOL = 1e-10
#: amplitudes within this of the maximum tie
TIE_ATOL = 1e-12
#: the largest |Im S| that cosine_amplitude drops
COSINE_IMAG_ATOL = 1e-12
#: no ordering vector when the largest amplitude is at most this
ORDER_FLOOR = 1e-10

UNITS = {"spin": "hbar (same convention as Sz in energy.dat)",
         "charge": "electrons"}


@dataclasses.dataclass(frozen=True)
class OrderParameters:
    """Order parameters of one magnetic cell (see the module docstring for
    the exact definitions).

    Attributes
    ----------
    site_of_orbital : ndarray (n_orb,) int
    site_positions : ndarray (n_site, 3) or None
        Fractional lattice coordinates.
    charge : ndarray (n_site,)
        ``n_s``.
    spin : ndarray (n_site, 3)
        ``S_s = (S^x, S^y, S^z)`` in units of hbar.
    moment : ndarray (n_site,)
        ``|S_s|``.
    uniform : ndarray (3,)
        ``M = N_site^-1 sum_s S_s``.
    uniform_norm, mean_local_moment : float
    q_frac : tuple of (Fraction, Fraction, Fraction), or None
    spin_fourier : ndarray (nq, 3) complex, or None
        ``M(Q)`` for each ``q_frac``.
    spin_amplitude : ndarray (nq,), or None
        ``|M(Q)|``.
    spin_ordering_vector : tuple of Fraction, or None
    charge_fourier : ndarray (nq,) complex, or None
        ``dn(Q)``.
    charge_amplitude : ndarray (nq,), or None
        ``|dn(Q)|``; 0 at ``Q = 0``.
    charge_ordering_vector : tuple of Fraction, or None
        Over ``Q != 0``.
    charge_order_amplitude : float or None
        The largest ``|dn(Q)|`` over ``Q != 0`` (0 when the cell has no
        ``Q != 0``).
    staggered : ndarray (3,) complex, or None
        ``M(q_AF)``.
    staggered_norm : float or None
    cell_shape, lattice_shape : tuple of int, or None
    n_cell : int
    transverse_available : bool
        False when the run cannot have spin-off-diagonal entries (UHFr with
        a fixed 2Sz); S^x and S^y are then zero by construction.
    units : dict
    source : Path or None
        ``None`` for results built by :func:`order_parameters_from_density`.
    meta : dict
        ``rho_convention``; ``fourier_unavailable`` and
        ``staggered_unavailable`` give the reason when those fields are
        ``None``.
    """
    site_of_orbital: np.ndarray
    site_positions: Optional[np.ndarray]
    charge: np.ndarray
    spin: np.ndarray
    moment: np.ndarray
    uniform: np.ndarray
    uniform_norm: float
    mean_local_moment: float
    q_frac: Optional[tuple]
    spin_fourier: Optional[np.ndarray]
    spin_amplitude: Optional[np.ndarray]
    spin_ordering_vector: Optional[tuple]
    charge_fourier: Optional[np.ndarray]
    charge_amplitude: Optional[np.ndarray]
    charge_ordering_vector: Optional[tuple]
    charge_order_amplitude: Optional[float]
    staggered: Optional[np.ndarray]
    staggered_norm: Optional[float]
    cell_shape: Optional[tuple]
    lattice_shape: Optional[tuple]
    n_cell: int
    transverse_available: bool
    units: dict
    source: Optional[Path]
    meta: dict


def _is_integer(value):
    if isinstance(value, (bool, np.bool_)):
        return False
    if isinstance(value, numbers.Integral):
        return True
    return (isinstance(value, numbers.Real) and math.isfinite(value)
            and float(value).is_integer())


def _positive_int(name, value):
    if not _is_integer(value) or int(value) < 1:
        raise ValueError("{} must be a positive integer, got {!r}".format(name, value))
    return int(value)


def _triple(name, value):
    if value is None:
        return None
    try:
        items = list(np.asarray(value, dtype=object).reshape(-1)) \
            if np.ndim(value) == 1 else None
    except (TypeError, ValueError):
        items = None
    if items is None or len(items) != 3 or not all(
            _is_integer(v) and int(v) >= 1 for v in items):
        raise ValueError("{} must be three positive integers, got {!r}"
                         .format(name, value))
    return tuple(int(v) for v in items)


def _validate_rho(rho):
    arr = np.asarray(rho)
    if not (np.issubdtype(arr.dtype, np.number) or arr.dtype == np.bool_):
        raise ValueError("rho must be a numeric array, got dtype {}".format(arr.dtype))
    if arr.ndim != 3 or arr.shape[1:] != (2, 2) or arr.shape[0] < 1:
        raise ValueError("rho must have shape (n_orb, 2, 2) with n_orb >= 1, "
                         "got {}".format(arr.shape))
    arr = np.array(arr, dtype=np.complex128, copy=True)
    if not np.all(np.isfinite(arr)):
        raise ValueError("rho has non-finite entries")
    dev = float(np.max(np.abs(arr - np.conj(np.swapaxes(arr, 1, 2)))))
    if dev > HERMITIAN_ATOL:
        raise ValueError("rho must be Hermitian in (s, t) within {:g}; "
                         "max |rho - rho^dagger| = {:.3g}".format(HERMITIAN_ATOL, dev))
    return arr


def _validate_positions(positions):
    arr = np.asarray(positions)
    if not np.issubdtype(arr.dtype, np.number) or np.iscomplexobj(arr) \
            or arr.ndim != 2 or arr.shape[1] != 3:
        raise ValueError("site_positions must be a real array of shape "
                         "(n_site, 3), got shape {}".format(arr.shape))
    arr = np.array(arr, dtype=float, copy=True)
    if not np.all(np.isfinite(arr)):
        raise ValueError("site_positions has non-finite entries")
    return arr


def _site_mapping(site_of_orbital, n_orb, positions, cell_shape):
    if site_of_orbital is None:
        if positions is not None:
            n_pos, where = len(positions), "site_positions"
        elif cell_shape is not None:
            n_pos, where = int(np.prod(cell_shape)), "cell_shape"
        else:
            raise ValueError(
                "site_of_orbital is required when neither site_positions nor "
                "cell_shape gives the number of cell positions")
        if n_orb > n_pos:
            raise ValueError(
                "site_of_orbital is required when a cell holds several "
                "orbitals ({} orbital entries, {} cell positions from {})"
                .format(n_orb, n_pos, where))
        if n_orb < n_pos:
            raise ValueError(
                "{} gives {} cell positions but rho has only {} orbital "
                "entries".format(where, n_pos, n_orb))
        return np.arange(n_orb)
    raw = np.asarray(site_of_orbital)
    if raw.ndim != 1 or len(raw) != n_orb:
        raise ValueError("site_of_orbital must have shape ({},), got {}"
                         .format(n_orb, raw.shape))
    if raw.dtype == np.bool_ or not np.issubdtype(raw.dtype, np.number) \
            or np.iscomplexobj(raw) or not all(_is_integer(v) for v in raw.tolist()):
        raise ValueError("site_of_orbital must hold integers, got {!r}"
                         .format(raw.tolist()))
    sites = raw.astype(np.int64)
    if np.any(sites < 0):
        raise ValueError("site_of_orbital values must be >= 0, got {!r}"
                         .format(sites.tolist()))
    n_site = int(sites.max()) + 1
    missing = sorted(set(range(n_site)) - set(sites.tolist()))
    if missing:
        raise ValueError("site_of_orbital must cover every site 0..{}; no "
                         "orbital is on site(s) {}".format(n_site - 1, missing))
    return sites


def _q_set(cell_shape):
    axes = [[Fraction(m, b) for m in range(-((b - 1) // 2), b // 2 + 1)]
            for b in cell_shape]
    return tuple(itertools.product(*axes))


def _ordering_vector(amplitude, q_frac, candidates):
    if not candidates:
        return None
    best = max(amplitude[i] for i in candidates)
    if best <= ORDER_FLOOR:
        return None
    tied = [q_frac[i] for i in candidates if amplitude[i] >= best - TIE_ATOL]
    return min(tied, key=lambda q: (sum(c * c for c in q), q))


def _phases(positions, q_frac):
    """``e^{i Q.r_s}``, shape ``(n_site, nq)``."""
    q = np.array([[float(c) for c in qq] for qq in q_frac], dtype=float)
    x = positions @ q.T
    x = x - np.floor(x)              # e^{2 pi i x} is periodic in x
    return np.exp(2j * np.pi * x)


def order_parameters_from_density(rho, *, site_of_orbital=None,
                                  site_positions=None, cell_shape=None,
                                  lattice_shape=None, n_cell=None,
                                  transverse_available=True) -> OrderParameters:
    """The order parameters of the one-cell density ``rho`` (the Green slot
    ``rho[o, s, t] = <c^dag_{o s} c_{o t}>``; see the module docstring for
    every definition).

    The inputs are copied; the result is frozen with read-only arrays and
    ``source = None``. Without ``site_positions``, ``cell_shape`` or
    ``lattice_shape`` the result has the per-site and uniform quantities
    only, with every Fourier, ordering and staggered field ``None``.

    Raises
    ------
    ValueError
        Naming the argument: ``rho`` not of shape ``(n, 2, 2)``, not finite
        or not Hermitian within 1e-10; ``site_of_orbital`` not of length
        ``n`` or not covering ``0..n_site-1`` (or missing when a cell holds
        several orbitals); ``site_positions`` not of shape ``(n_site, 3)``;
        ``cell_shape``/``lattice_shape`` not positive integer triples or
        ``lattice_shape`` not divisible by ``cell_shape``; ``n_cell`` not a
        positive integer, or not ``prod(lattice_shape) / prod(cell_shape)``
        when both shapes are given.
    """
    rho = _validate_rho(rho)
    n_orb = rho.shape[0]
    cell_shape = _triple("cell_shape", cell_shape)
    lattice_shape = _triple("lattice_shape", lattice_shape)
    if cell_shape is not None and lattice_shape is not None:
        bad = [a for a in range(3) if lattice_shape[a] % cell_shape[a]]
        if bad:
            raise ValueError("lattice_shape {} must be divisible by cell_shape {} "
                             "on every axis (not on axis {})"
                             .format(lattice_shape, cell_shape, bad))
    if n_cell is not None:
        n_cell = _positive_int("n_cell", n_cell)
    if cell_shape is not None and lattice_shape is not None:
        cells = int(np.prod(lattice_shape)) // int(np.prod(cell_shape))
        if n_cell is None:
            n_cell = cells
        elif n_cell != cells:
            raise ValueError(
                "n_cell {} does not match lattice_shape {} / cell_shape {} = "
                "{} magnetic cells".format(n_cell, lattice_shape, cell_shape, cells))
    elif n_cell is None:
        n_cell = 1
    positions = None if site_positions is None else _validate_positions(site_positions)
    sites = _site_mapping(site_of_orbital, n_orb, positions, cell_shape)
    n_site = int(sites.max()) + 1
    if positions is not None and positions.shape[0] != n_site:
        raise ValueError("site_positions must have shape ({}, 3) for {} sites, "
                         "got {}".format(n_site, n_site, positions.shape))

    rho_site = np.zeros((n_site, 2, 2), dtype=np.complex128)
    np.add.at(rho_site, sites, rho)
    up_down, down_up = rho_site[:, 0, 1], rho_site[:, 1, 0]
    charge = np.real(rho_site[:, 0, 0] + rho_site[:, 1, 1])
    spin = np.stack([0.5 * np.real(up_down + down_up),
                     0.5 * np.real(-1j * (up_down - down_up)),
                     0.5 * np.real(rho_site[:, 0, 0] - rho_site[:, 1, 1])], axis=1)
    moment = np.linalg.norm(spin, axis=1)
    uniform = spin.mean(axis=0)

    meta = {"rho_convention": "rho[o, s, t] = <c^dag_{o s} c_{o t}> "
                              "(the stored Green slot G_ab = <c^dag_a c_b>)"}
    fields = dict(q_frac=None, spin_fourier=None, spin_amplitude=None,
                  spin_ordering_vector=None, charge_fourier=None,
                  charge_amplitude=None, charge_ordering_vector=None,
                  charge_order_amplitude=None, staggered=None,
                  staggered_norm=None)
    absent = [name for name, value in (("site_positions", positions),
                                       ("cell_shape", cell_shape),
                                       ("lattice_shape", lattice_shape))
              if value is None]
    if absent:
        meta["fourier_unavailable"] = "no {}".format(" and ".join(absent))
    else:
        q_frac = _q_set(cell_shape)
        phase = _phases(positions, q_frac)
        spin_fourier = phase.T @ spin / n_site
        charge_fourier = phase.T @ (charge - charge.mean()) / n_site
        spin_amplitude = np.sqrt(np.sum(np.abs(spin_fourier) ** 2, axis=1))
        charge_amplitude = np.abs(charge_fourier)
        zero = q_frac.index((Fraction(0),) * 3)
        charge_amplitude[zero] = 0.0
        nonzero = [i for i in range(len(q_frac)) if i != zero]
        fields.update(
            q_frac=q_frac, spin_fourier=spin_fourier,
            spin_amplitude=spin_amplitude,
            spin_ordering_vector=_ordering_vector(
                spin_amplitude, q_frac, list(range(len(q_frac)))),
            charge_fourier=charge_fourier, charge_amplitude=charge_amplitude,
            charge_ordering_vector=_ordering_vector(
                charge_amplitude, q_frac, nonzero),
            charge_order_amplitude=float(max(
                (charge_amplitude[i] for i in nonzero), default=0.0)))
        extended = [a for a in range(3) if lattice_shape[a] > 1]
        odd = [a for a in extended if cell_shape[a] % 2]
        if not extended:
            meta["staggered_unavailable"] = "no axis with lattice_shape > 1"
        elif odd:
            meta["staggered_unavailable"] = (
                "cell_shape {} is odd on axis {} where lattice_shape {} > 1"
                .format(cell_shape, odd, lattice_shape))
        else:
            q_af = tuple(Fraction(1, 2) if a in extended else Fraction(0)
                         for a in range(3))
            i = q_frac.index(q_af)
            fields.update(staggered=spin_fourier[i],
                          staggered_norm=float(spin_amplitude[i]))

    return OrderParameters(
        site_of_orbital=_freeze(sites), site_positions=_freeze(positions),
        charge=_freeze(charge), spin=_freeze(spin), moment=_freeze(moment),
        uniform=_freeze(uniform), uniform_norm=float(np.linalg.norm(uniform)),
        mean_local_moment=float(moment.mean()),
        **{k: _freeze(v) for k, v in fields.items()},
        cell_shape=cell_shape, lattice_shape=lattice_shape, n_cell=n_cell,
        transverse_available=bool(transverse_available), units=dict(UNITS),
        source=None, meta=meta)


def cosine_amplitude(S, site_positions, q):
    """``N_site^-1 sum_s cos(Q.r_s) S_s`` with ``Q.r_s = 2 pi q.r_s``.

    ``S`` has the site as its first axis (a vector or a scalar per site);
    ``site_positions`` is ``(n_site, 3)`` fractional, ``q`` three fractional
    components. It equals ``M(Q)`` only when every ``Q.r_s`` is in ``pi Z``;
    otherwise use ``OrderParameters.spin_fourier``.

    ``S`` must be real: a complex ``S`` is accepted only when every
    imaginary part is at most ``COSINE_IMAG_ATOL`` (1e-12) in absolute
    value, and is then dropped; otherwise ``ValueError``. A bad ``q`` (not
    three finite real numbers) is a ``ValueError`` too."""
    values = np.asarray(S)
    if values.dtype == np.bool_ or not np.issubdtype(values.dtype, np.number):
        raise ValueError("S must be a real numeric array, got dtype {}"
                         .format(values.dtype))
    if np.iscomplexobj(values):
        imag = float(np.max(np.abs(values.imag), initial=0.0))
        if not imag <= COSINE_IMAG_ATOL:
            raise ValueError("S must be real (|Im S| <= {:g}); max |Im S| = {:.3g}"
                             .format(COSINE_IMAG_ATOL, imag))
        values = values.real
    values = np.array(values, dtype=float)
    positions = _validate_positions(site_positions)
    if values.ndim < 1 or values.shape[0] != positions.shape[0]:
        raise ValueError("S must have the site as its first axis ({} sites), "
                         "got shape {}".format(positions.shape[0], values.shape))
    try:
        items = list(q)
    except TypeError:
        items = None
    if items is None or len(items) != 3 or not all(
            isinstance(c, numbers.Real) and not isinstance(c, (bool, np.bool_))
            and math.isfinite(c) for c in items):
        raise ValueError("q must be three finite real numbers, got {!r}".format(q))
    qv = np.array([float(c) for c in items], dtype=float)
    cos = np.cos(2 * np.pi * (positions @ qv))
    return np.tensordot(cos, values, axes=(0, 0)) / positions.shape[0]


def sum_rule_tolerance(residual, *, mix, g_size, nvol, nd):
    """The bound on ``|sum of the saved Green function - energy.dat total|``:
    ``(1 - mix) nvol sqrt(nd) residual g_size / 2 + nvol nd 1e-12 + 1e-12``.

    ``residual`` is the last ``Rest = 2 ||dG||_F / g_size`` (from
    ``scf_status.dat``), ``mix`` the ``Mix`` parameter, ``g_size`` the
    number of elements of the Green function array, ``nvol`` the number of
    cells and ``nd`` the orbital-spin dimension of one cell. For UHFr:
    ``g_size = 4 N^2``, ``nvol = 1``, ``nd = 2 N``."""
    if isinstance(residual, bool) or not isinstance(residual, numbers.Real) \
            or not math.isfinite(residual) or residual < 0:
        raise ValueError("residual must be a finite number >= 0, got {!r}"
                         .format(residual))
    if isinstance(mix, bool) or not isinstance(mix, numbers.Real) \
            or not math.isfinite(mix):
        raise ValueError("mix must be a finite number, got {!r}".format(mix))
    g_size = _positive_int("g_size", g_size)
    nvol = _positive_int("nvol", nvol)
    nd = _positive_int("nd", nd)
    return float((1.0 - mix) * nvol * math.sqrt(nd) * residual * g_size / 2
                 + nvol * nd * 1e-12 + 1e-12)
