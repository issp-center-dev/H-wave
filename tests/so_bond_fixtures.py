"""Shared test fixtures of the spin-orbital vertex tests that need the
production readers (``tests/so_bond_oracles.py`` stays free of solver
imports).

* :func:`onsite_decl` -- the Hermitian-closed on-site declaration of a type.
* :func:`rpa_spinful_gamma` -- the antisymmetrized spinful vertex of the
  spin-orbital RPA (``rpa.spinful_bubble_vertex`` of ``ham_inter_q`` and
  ``spinful_exchange()``), read from interaction files.
* :func:`production_gamma_b` -- the production ``Gamma_B`` of a table set
  through a vertex-only admission plan.
* :func:`green_kw`, :func:`pair_kernel` -- the Green-function layout of
  the bond bubble and the bond pairing kernel of code-slot vertex blocks.
* :func:`exact` -- an anchored regex of an exact refusal message.
* :func:`chibar_static`, :func:`production_gamma_label`,
  :func:`production_first_order` -- the production spinful bond vertex in
  the label form of an exact-diagonalization fixture and its first-order
  correlator ``u^T [-chibar(q) Gamma(-q) chibar(q)] v``.
"""
import contextlib
import functools
import logging
import os
import re
import shutil
import tempfile

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
from hwave.solver import eliashberg_bond as eb
from hwave.solver import rpa
from hwave.solver.so_bond import (SLOT_CONVENTION, SOBondCaps, plan_so_bond, slot_geometry,
                                  so_bond_topology, so_bond_vertex)
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import (CELL, first_order_prediction, gamma_code_to_label,
                                   lindhard_static, mesh_index)


@contextlib.contextmanager
def quiet():
    """Silence log records up to WARNING inside the block."""
    prev = logging.root.manager.disable
    logging.disable(logging.WARNING)
    try:
        yield
    finally:
        logging.disable(prev)


def _wfile(path, head, n, rows, value_format):
    line = "  %d %d %d %d %d  " + value_format + "  " + value_format + "\n"
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(head, n))
        for r in rows:
            f.write(line % r)


def onsite_decl(t, v):
    """The Hermitian-closed on-site declaration of type ``t`` (two physical
    orbitals): the diagonal rows of CoulombIntra, the inter-orbital pair
    ``(0, 1)`` / ``(1, 0)`` of every other type."""
    if t == "CoulombIntra":
        return {t: {((0, 0, 0), (0, 0)): v, ((0, 0, 0), (1, 1)): 0.8 * v}}
    return {t: {((0, 0, 0), (0, 1)): v, ((0, 0, 0), (1, 0)): np.conj(v)}}


def rpa_spinful_gamma(decl, norb, shape, *, beta=2.0, value_format="%.17g",
                      orbital_site=None):
    """``Gamma (nvol, nd**2, nd**2)`` (``nd = 2 * norb``) of the spin-orbital
    RPA built from the declarations ``{type: {((Rx, Ry, Rz), (a, b)): v}}``
    (physical orbital indices, as the interaction files are read in
    spin-orbital mode) on the mesh ``shape``. ``value_format`` is the printf
    format of the coefficients in the files; ``orbital_site`` is passed to
    the solver parameters when given."""
    nd = 2 * norb
    d = tempfile.mkdtemp(prefix="so_vertex_reader_")
    try:
        with open(os.path.join(d, "geom.dat"), "w") as f:
            f.write("1 0 0\n0 1 0\n0 0 1\n%d\n" % nd + "0 0 0\n" * nd)
        _wfile(os.path.join(d, "transfer.dat"), "Transfer", nd,
               [(0, 0, 0, p + 1, p + 1, 0.1, 0.0) for p in range(nd)], value_format)
        inter = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat"}
        for t, tbl in decl.items():
            _wfile(os.path.join(d, t + ".dat"), t, norb,
                   [(R[0], R[1], R[2], a + 1, b + 1, complex(v).real, complex(v).imag)
                    for (R, (a, b)), v in tbl.items()], value_format)
            inter[t] = t + ".dat"
        io = read_input_k.QLMSkInput({"path_to_input": d, "interaction": inter})
        par = {"T": 1.0 / beta, "mu": 0.0, "CellShape": list(shape),
               "SubShape": [1, 1, 1], "Nmat": 4}
        if orbital_site is not None:
            par["orbital_site"] = orbital_site
        with quiet():
            s = rpa.RPA(io.get_param("ham"), {}, {
                "mode": "RPA", "param": par, "enable_spin_orbital": True,
                "calc_scheme": "general"})
        nvol = int(np.prod(shape))
        ham = np.asarray(s.ham_info.ham_inter_q).reshape(nvol, nd, nd, nd, nd)
        return rpa.spinful_bubble_vertex(ham, s.ham_info.spinful_exchange())
    finally:
        shutil.rmtree(d, ignore_errors=True)


def exact(msg):
    """The regex that matches exactly ``msg``."""
    return "^" + re.escape(msg) + "$"


def production_gamma_b(tables, norb, topo, shape):
    """The production ``Gamma_B`` ``(nvol, ND, ND)`` of the interaction
    ``tables`` (``norb`` physical orbitals) on the channels of ``topo`` (any
    object with ``delta_r`` and ``reverse``) and the mesh ``shape``, built
    through a vertex-only admission plan."""
    c = compile_interaction(tables, norb)
    plan = plan_so_bond(nmat=2, nvol=int(np.prod(shape)), B=len(topo.delta_r), nd=2 * norb,
                        consumers=("vertex",), caps=SOBondCaps(None, None, None))
    return so_bond_vertex(c, topo, shape, plan=plan)


def green_kw(g):
    """sc layout ``(nd, nd, Nx, 1, 1, nmat)`` -> the bond-bubble layout
    ``(1, nmat, Nx, nd, nd)`` (a chain along x)."""
    return np.ascontiguousarray(g[:, :, :, 0, 0, :].transpose(3, 2, 0, 1))[None]


def pair_kernel(blocks, G2, view, nd, shape, beta):
    """The bond pairing kernel of the code-slot vertex blocks ``(nfreq, nvol,
    B nd^2, B nd^2)`` on the uniform bosonic grid of the pair bubble ``G2``
    ``(nd, nd, nd, nd, Nx, Ny, Nz, nfreq)``, no instantaneous vertex (a
    frequency-independent part sits in every block)."""
    nvol = int(np.prod(shape))
    B = int(view.n_channels)
    nfreq = int(blocks.shape[0])
    vertex = eb.PairVertexUniform(eb.ArrayBlockSource({"V": blocks}, nd * nd), "V", B,
                                  nd * nd, nfreq, nvol)
    return eb.BondPairKernel(vertex, G2, view, xp=np, spatial_shape=shape, norb=nd,
                             beta=beta, nfreq=nfreq, V_inst=None, residency="host",
                             host_cap=np.inf)


@functools.lru_cache(maxsize=None)
def chibar_static(fixture):
    """The analytic static label bubble of a fixture's band (read-only
    cache by content)."""
    return lindhard_static(*fixture.band, fixture.geom)


@functools.lru_cache(maxsize=None)
def production_gamma_label(fixture, frame):
    """Production ``Gamma_B`` of ``fixture``'s tables, mapped to the label
    form by ``gamma_code_to_label`` on the production channels and embedded
    into the label geometry of ``frame`` (a fixture with the same mesh and
    band; channels matched by displacement modulo the mesh)."""
    nd, shape = fixture.nd, fixture.spatial_shape
    tables = fixture.tables_dict
    topo = so_bond_topology(compile_interaction(tables, fixture.norb_phys), CELL)
    gb = production_gamma_b(tables, fixture.norb_phys, topo, shape)
    geom_p = slot_geometry(topo, shape)
    gl_p = gamma_code_to_label(gb, geom_p, nd, SLOT_CONVENTION.slot_sign)
    geom = frame.geom
    chan = {int(mesh_index(np.asarray(d), shape)): m for m, d in enumerate(geom.delta_r)}
    npair = nd * nd
    idx = np.concatenate([chan[int(mesh_index(np.asarray(d), shape))] * npair + np.arange(npair)
                          for d in geom_p.delta_r])
    NL = len(geom.delta_r) * npair
    out = np.zeros((fixture.nvol, NL, NL), dtype=np.complex128)
    out[:, idx[:, None], idx[None, :]] = gl_p
    out.flags.writeable = False
    return out


def production_first_order(fixture, frame, sketches, qs):
    """``u^T [-chibar(q) Gamma(-q) chibar(q)] v`` ``(len(qs), len(sketches))``
    of the production vertex of ``fixture``'s tables on the labels and band
    of ``frame``."""
    return first_order_prediction(chibar_static(frame), production_gamma_label(fixture, frame),
                                  frame.spatial_shape, sketches, qs)
