"""Public-entry fixtures of the spin-orbital bond path, an independent
Green function and an analytic Zeeman-sign oracle.

The fixture tables, :func:`independent_green`, :func:`analytic_dispersion`
and :func:`zeeman_sign_miss` share no code with ``hwave``: at module level
this module imports numpy and the standard library only (the test-side
:mod:`tests.so_bond_oracles` is the one further import allowed). A helper that
drives production code imports ``hwave`` inside its own body.

Conventions (stated here, not read from the solver)
---------------------------------------------------
* Every fixture is a 4x4x1 mesh with one physical orbital (``nd = 2``) unless
  stated. Hopping is listed literally as ``(R, t)`` pairs, real and spin
  independent, with the Hermitian partner listed at ``-R``.
* Transfer entries carry the interleaved spin-orbital index ``A = 2 a + s``
  (``s = 0`` up), the spin-orbital file format; interaction and Extern
  entries carry physical orbital indices; all are 0-based here and 1-based in
  the files written by :func:`write_case`.
* The one-body band is ``H_XY(k) = sum_R t_XY(R) e^{+i k.R}`` (the
  Wannier90 sign), with ``k = 2 pi (i/Nx, j/Ny, l/Nz)``, in the spin-block
  index ``X = s * norb_phys + a``. An Extern table ``H1(R)`` enters as
  ``coeff_extern * H1(k) (x) sigma^z``: ``+h`` on spin up, ``-h`` on spin
  down.
* Green functions are ``(nd, nd, Nx, Ny, Nz, nmat)``, ``G_XY(k, i w_n) =
  sum_i v_Xi conj(v_Yi) / (i w_n + mu - e_i(k))`` on the fermionic grid
  ``w_n = (2n + 1 - nmat) pi / beta``, ``n = 0 .. nmat - 1`` (that is
  ``(2m + 1) pi / beta``, ``m = -nmat/2 .. nmat/2 - 1``).
* ``filling`` is electrons per orbital per spin: the electron count is
  ``filling * nd * nvol``.
"""

import os
from typing import NamedTuple

import numpy as np


class EntryFixture(NamedTuple):
    """A public-entry fixture, hashable by content.

    ``hopping``: ``((R, t), ...)``, spin independent, on orbital 0;
    ``transfer``: further Transfer entries ``((R, (A, B), v), ...)`` in the
    interleaved index (spin-flip hopping, a Zeeman field written as Transfer);
    ``extern``: the Extern table ``((R, (a, b), v), ...)`` in physical
    orbitals (``()``: no Extern file) with ``coeff_extern``; ``interactions``:
    ``((type, ((R, (a, b), v), ...)), ...)`` in physical orbitals."""
    name: str
    shape: tuple
    norb_phys: int
    hopping: tuple
    transfer: tuple
    extern: tuple
    coeff_extern: float
    interactions: tuple
    filling: float

    @property
    def nd(self):
        return 2 * self.norb_phys

    @property
    def nvol(self):
        return int(np.prod(self.shape))


_SHAPE = (4, 4, 1)

#: nearest-neighbour t1 = -1.0 and next-nearest-neighbour t2 = 0.25
_HOPPING = (
    ((1, 0, 0), -1.0), ((-1, 0, 0), -1.0), ((0, 1, 0), -1.0), ((0, -1, 0), -1.0),
    ((1, 1, 0), 0.25), ((-1, -1, 0), 0.25), ((1, -1, 0), 0.25), ((-1, 1, 0), 0.25),
)

#: nearest-neighbour spin-flip hopping (interleaved (up, dn) = (0, 1)), 0.3 per
#: bond: the Rashba-like 0.6 (sin ky sigma_x - sin kx sigma_y) in k space,
#: Hermitian-closed
_SOC = (
    ((1, 0, 0), (0, 1), 0.3), ((1, 0, 0), (1, 0), -0.3),
    ((-1, 0, 0), (0, 1), -0.3), ((-1, 0, 0), (1, 0), 0.3),
    ((0, 1, 0), (0, 1), -0.3j), ((0, 1, 0), (1, 0), -0.3j),
    ((0, -1, 0), (0, 1), 0.3j), ((0, -1, 0), (1, 0), 0.3j),
)


def _nn(v):
    """A nearest-neighbour row ``v`` on orbital 0, on all four bonds."""
    return (((1, 0, 0), (0, 0), v), ((-1, 0, 0), (0, 0), v),
            ((0, 1, 0), (0, 0), v), ((0, -1, 0), (0, 0), v))


#: on-site CoulombIntra and the nearest-neighbour rows of every other type
_ALL_TYPES = (
    ("CoulombIntra", (((0, 0, 0), (0, 0), 0.83),)),
    ("CoulombInter", _nn(0.37)),
    ("Hund", _nn(0.53)),
    ("Ising", _nn(0.11)),
    ("Exchange", _nn(0.31)),
    ("PairHop", _nn(0.43)),
    ("PairLift", _nn(0.07)),
)

#: the S_z-conserving interactions of the Zeeman fixtures
_SZ_TYPES = (
    ("CoulombIntra", (((0, 0, 0), (0, 0), 0.83),)),
    ("CoulombInter", _nn(0.37)),
    ("Ising", _nn(0.11)),
)

#: the field of the Zeeman fixtures
ZEEMAN_H = 0.05

#: t1/t2 hopping plus nearest-neighbour spin-flip SOC; every interaction type
SOC_ALL_TYPES = EntryFixture("SOC_ALL_TYPES", _SHAPE, 1, _HOPPING, _SOC, (), 0.0,
                             _ALL_TYPES, 0.4)

#: the band of :data:`SOC_ALL_TYPES`; only the R = 0 interaction rows
SOC_ONSITE = SOC_ALL_TYPES._replace(
    name="SOC_ONSITE",
    interactions=tuple((t, tuple(r for r in rows if r[0] == (0, 0, 0)))
                       for t, rows in _ALL_TYPES
                       if any(r[0] == (0, 0, 0) for r in rows)))

#: t1/t2 hopping on both spins, the Zeeman field as an Extern file (H1(R = 0)
#: = 1 on physical orbital 1, coeff_extern = ZEEMAN_H)
ZEEMAN_EXTERN = EntryFixture("ZEEMAN_EXTERN", _SHAPE, 1, _HOPPING, (),
                             (((0, 0, 0), (0, 0), 1.0),), ZEEMAN_H, _SZ_TYPES, 0.4)

#: the same field as on-site Transfer entries: +ZEEMAN_H on spin up,
#: -ZEEMAN_H on spin down
ZEEMAN_TRANSFER = EntryFixture("ZEEMAN_TRANSFER", _SHAPE, 1, _HOPPING,
                               (((0, 0, 0), (0, 0), ZEEMAN_H),
                                ((0, 0, 0), (1, 1), -ZEEMAN_H)),
                               (), 0.0, _SZ_TYPES, 0.4)

#: the SU(2)-symmetric fixture: spin-independent t1/t2 hopping, on-site
#: CoulombIntra and the nearest-neighbour CoulombInter, Hund, Ising, Exchange
#: (= J_H - 2 J_I) and real PairHop rows
SU2_4X4 = EntryFixture("SU2_4X4", _SHAPE, 1, _HOPPING, (), (), 0.0, (
    ("CoulombIntra", (((0, 0, 0), (0, 0), 0.83),)),
    ("CoulombInter", _nn(0.37)),
    ("Hund", _nn(0.53)),
    ("Ising", _nn(0.11)),
    ("Exchange", _nn(0.31)),
    ("PairHop", _nn(0.43)),
), 0.4)


# ------------------------------------------------------------ the files
def _transfer_rows(fixture):
    """Every Transfer entry ``(R, (A, B), v)`` of ``fixture`` (interleaved),
    the hopping on both spins first."""
    rows = []
    for R, t in fixture.hopping:
        for s in (0, 1):
            rows.append((R, (s, s), t))
    rows.extend(fixture.transfer)
    return rows


def _write_table(path, head, n, rows):
    """A sparse Wannier90-style table (degeneracy 1, 1-based indices)."""
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(head, n))
        for R, (a, b), v in rows:
            v = complex(v)
            f.write("  %d %d %d %d %d  %r  %r\n" % (R[0], R[1], R[2], a + 1, b + 1,
                                                    v.real, v.imag))


def _override(section, values, where):
    """Set ``values`` into the dict ``section``, keys matched
    case-insensitively: a key already present keeps its spelling and takes
    the new value, so a section never holds two spellings of one key."""
    for key, value in (values or {}).items():
        if not isinstance(key, str):
            raise ValueError("write_case: {} keys must be strings, got {!r}".format(where, key))
        same = [k for k in section if k.lower() == key.lower()]
        section[same[0] if same else key] = value
    return section


def write_case(tmpdir, fixture, *, param=None, eliashberg=None, mode=None):
    """Write the input files of ``fixture`` into ``tmpdir`` and return the
    ``input_dict`` of an hwave_sc static bond run (beta 2.0, Nmat 32).

    ``param``, ``eliashberg`` and ``mode`` override the ``[mode.param]``,
    ``[eliashberg]`` and ``[mode]`` tables of that dict; keys are matched
    case-insensitively, a key already present keeps its spelling (the
    defaults spell ``T``, ``CellShape``, ``SubShape``, ``Nmat`` and
    ``filling`` as the run reads them) and any other key is added as given
    (``param={"coeff_tail": ...}`` lands in ``[mode.param]``). The keys the
    run reads case-insensitively are spelled in mixed case.

    ``[mode] enable_spin_orbital`` (default true) selects the files: true
    writes the spin-orbital files (``nd = 2 norb_phys`` spin-orbitals, the
    interleaved Transfer table, the Extern file); false writes the
    non-spin-orbital files of the spin-independent hopping (``norb_phys``
    orbitals), refused for a fixture with spin-dependent Transfer entries or
    an Extern table."""
    if mode is not None and any(k.lower() == "param" for k in mode):
        raise ValueError("write_case: pass the [mode.param] overrides as param=")
    mode_tbl = _override({"mode": "RPA", "Enable_Spin_Orbital": True}, mode, "[mode]")
    so = [v for k, v in mode_tbl.items() if k.lower() == "enable_spin_orbital"][0]
    if not isinstance(so, bool):
        raise ValueError("write_case: [mode] enable_spin_orbital must be a bool, got "
                         "{!r}".format(so))
    if not so and (fixture.transfer or fixture.extern):
        raise ValueError("write_case: fixture {} has spin-dependent Transfer entries or an "
                         "Extern table; it has no non-spin-orbital files".format(fixture.name))
    n = fixture.nd if so else fixture.norb_phys
    with open(os.path.join(tmpdir, "geom.dat"), "w") as f:
        f.write("1 0 0\n0 1 0\n0 0 1\n{}\n".format(n))
        f.write("0 0 0\n" * n)
    rows = (_transfer_rows(fixture) if so
            else [(R, (0, 0), t) for R, t in fixture.hopping])
    _write_table(os.path.join(tmpdir, "transfer.dat"), "Transfer", n, rows)
    inter = {"path_to_input": tmpdir, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for typ, table in fixture.interactions:
        name = typ.lower() + ".dat"
        _write_table(os.path.join(tmpdir, name), typ, fixture.norb_phys, table)
        inter[typ.upper() if typ == "Ising" else typ.lower()] = name
    par = {"T": 0.5, "CellShape": list(fixture.shape), "SubShape": [1, 1, 1],
           "Nmat": 32, "filling": fixture.filling}
    if fixture.extern:
        _write_table(os.path.join(tmpdir, "extern.dat"), "Extern", fixture.norb_phys,
                     fixture.extern)
        inter["Extern"] = "extern.dat"
        par["COEFF_EXTERN"] = fixture.coeff_extern
    _override(par, param, "[mode.param]")
    eli = _override({"Frequency": "static", "Bond_Channels": True,
                     "Solver_Mode": "eigenvalue", "Num_Eigenvalues": 2},
                    eliashberg, "[eliashberg]")
    mode_tbl["param"] = par
    return {"mode": mode_tbl,
            "file": {"input": {"path_to_input": "", "interaction": inter},
                     "output": {"path_to_output": os.path.join(tmpdir, "out")}},
            "eliashberg": eli}


# ------------------------------------------------------------ the oracles
def _momenta(shape):
    """``k`` of every mesh point, ``(Nx, Ny, Nz, 3)``."""
    axes = [2.0 * np.pi * np.arange(n) / n for n in shape]
    return np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)


def _phase(k, R):
    return np.exp(1j * (k @ np.asarray(R, dtype=float)))


def _band(fixture):
    """``H(k)`` ``(Nx, Ny, Nz, nd, nd)`` in the spin-block index, from the
    literal lists."""
    n = fixture.norb_phys
    nd = fixture.nd
    k = _momenta(fixture.shape)
    H = np.zeros(tuple(fixture.shape) + (nd, nd), dtype=complex)
    for R, (A, B), v in _transfer_rows(fixture):
        X, Y = (A % 2) * n + A // 2, (B % 2) * n + B // 2
        H[..., X, Y] += v * _phase(k, R)
    for R, (a, b), v in fixture.extern:
        h1 = fixture.coeff_extern * v * _phase(k, R)
        H[..., a, b] += h1
        H[..., n + a, n + b] -= h1
    return H


def _fermi(x, beta):
    return 0.5 * (1.0 - np.tanh(0.5 * beta * x))


def _matsubara(nmat, beta):
    return (2.0 * np.arange(nmat) + 1.0 - nmat) * np.pi / beta


def independent_green(fixture, *, nmat, beta):
    """``(G, mu)`` of ``fixture`` from its literal lists: ``H(k)`` diagonalized
    with ``numpy.linalg.eigh``, mu by bisection on the count ``filling * nd *
    nvol`` (to an interval of 1e-14), ``G`` in the ``(nd, nd, Nx, Ny, Nz,
    nmat)`` layout."""
    e, v = np.linalg.eigh(_band(fixture))
    target = fixture.filling * fixture.nd * fixture.nvol
    lo, hi = float(e.min()) - 10.0, float(e.max()) + 10.0
    while hi - lo > 1e-14:
        mid = 0.5 * (lo + hi)
        if mid in (lo, hi):
            break
        if _fermi(e - mid, beta).sum() < target:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    iw = 1j * _matsubara(nmat, beta)
    g = 1.0 / (iw[None, None, None, None, :] + mu - e[..., None])     # (x, y, z, i, n)
    G = np.einsum("xyzai,xyzbi,xyzin->abxyzn", v, v.conj(), g)
    return G, mu


def analytic_dispersion(fixture, kx, ky, kz):
    """``sum_(R, t) t cos(k.R)`` over the spin-independent hopping list (each
    listed ``R`` once, so a ``+-R`` pair gives ``2 t cos(k.R)``)."""
    kx, ky, kz = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (kx, ky, kz)))
    out = np.zeros(kx.shape)
    for R, t in fixture.hopping:
        out = out + float(np.real(t)) * np.cos(kx * R[0] + ky * R[1] + kz * R[2])
    return out


def zeeman_sign_miss(G, mu, *, fixture, h, beta, nmat):
    """``max |G - G_expected| / max |G_expected|`` with ``G_expected`` the
    spin-diagonal Green function of the analytic dispersion in a field ``h``
    on spin up: ``G_00 = 1/(i w_n + mu - e(k) - h)``, ``G_11 = 1/(i w_n + mu -
    e(k) + h)``, off-diagonal 0 (one physical orbital)."""
    if fixture.norb_phys != 1:
        raise ValueError("zeeman_sign_miss: one physical orbital only")
    k = _momenta(fixture.shape)
    eps = analytic_dispersion(fixture, k[..., 0], k[..., 1], k[..., 2])[..., None]
    iw = 1j * _matsubara(nmat, beta)
    expected = np.zeros((2, 2) + tuple(fixture.shape) + (nmat,), dtype=complex)
    expected[0, 0] = 1.0 / (iw + mu - eps - h)
    expected[1, 1] = 1.0 / (iw + mu - eps + h)
    G = np.asarray(G)
    if G.shape != expected.shape:
        raise ValueError("zeeman_sign_miss: G has shape {}, expected {}".format(
            G.shape, expected.shape))
    return float(np.abs(G - expected).max() / np.abs(expected).max())
