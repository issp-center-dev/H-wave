"""Input builders for the ``hwave.initial`` tests.

Each builder writes a UHFk input directory (``input.toml``, ``geom.dat``,
``transfer.dat`` and one file per interaction) and returns the TOML path.
Orbital indices are 0-based here and written 1-based. A Transfer entry
``(R, a, b): v`` is the coefficient ``v`` of ``c^dag_a c_b`` on the bond
``R`` (``H0 = sum v c^dag_a c_b``), so a hopping ``t`` is written as
``-t`` (tests/test_uhf_input_regressions.py, "Trans is H = -sum t
c^dag c").
"""
import os

from tests.initial_reference_capture import dump_toml


def _w90(path, comment, norb, rows):
    """A Wannier90-like table ``{(R, a, b): complex}`` (0-based orbitals)."""
    with open(path, "w") as f:
        f.write("{}\n{}\n1\n 1\n".format(comment, norb))
        for ((rx, ry, rz), a, b), v in rows.items():
            v = complex(v)
            f.write("{} {} {} {} {} {:.15g} {:.15g}\n".format(
                rx, ry, rz, a + 1, b + 1, v.real, v.imag))


def write_input(dest, *, cellshape, transfer, norb=1, centers=None, rvec=None,
                interactions=None, param=None, mode=None, finput=None,
                extra_files=None, foutput=None):
    """Write a UHFk input under ``dest`` and return the TOML path.

    ``transfer`` and each ``interactions[name]`` are ``{(R, a, b): value}``;
    ``param`` is added to ``[mode.param]`` (``CellShape`` is set from
    ``cellshape``), ``mode`` to ``[mode]``, ``finput`` to ``[file.input]``,
    ``foutput`` to ``[file.output]``.
    ``extra_files`` maps a file name under ``dest`` to its text.
    """
    os.makedirs(dest, exist_ok=True)
    rvec = rvec if rvec is not None else ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    centers = centers if centers is not None else [(0.0, 0.0, 0.0)] * norb
    with open(os.path.join(dest, "geom.dat"), "w") as f:
        for v in rvec:
            f.write("{:.15g} {:.15g} {:.15g}\n".format(*v))
        f.write("{}\n".format(norb))
        for c in centers:
            f.write("{:.15g} {:.15g} {:.15g}\n".format(*c))
    _w90(os.path.join(dest, "transfer.dat"), "Transfer", norb, transfer)
    inter = {"path_to_input": "./", "Geometry": "geom.dat", "Transfer": "transfer.dat"}
    for name, rows in (interactions or {}).items():
        fname = name.lower() + ".dat"
        _w90(os.path.join(dest, fname), name, norb, rows)
        inter[name] = fname
    for name, text in (extra_files or {}).items():
        with open(os.path.join(dest, name), "w") as f:
            f.write(text)
    p = {"CellShape": list(cellshape), "IterationMax": 1, "EPS": 8, "Mix": 0.5,
         "RndSeed": 1}
    p.update(param or {})
    m = {"mode": "UHFk"}
    m.update(mode or {})
    m["param"] = p
    fin = {"path_to_input": ""}
    fin.update(finput or {})
    fin["interaction"] = inter
    fout = {"path_to_output": "output"}
    fout.update(foutput or {})
    doc = {"log": {"print_level": 0}, "mode": m,
           "file": {"input": fin, "output": fout}}
    toml = os.path.join(dest, "input.toml")
    with open(toml, "w") as f:
        f.write(dump_toml(doc))
    return toml


def square_transfer(t=1.0):
    """Nearest-neighbour hopping ``t`` of one orbital on the square lattice."""
    return {((1, 0, 0), 0, 0): -t, ((-1, 0, 0), 0, 0): -t,
            ((0, 1, 0), 0, 0): -t, ((0, -1, 0), 0, 0): -t}


def square(dest, L, *, U=0.0, T=0.0, ncond=None, subshape=None, fock=True,
           two_sz=None, param=None):
    """One-orbital square lattice ``L x L``, hopping 1, on-site ``U``;
    half filling unless ``ncond`` is given."""
    p = {"Ncond": L * L if ncond is None else ncond, "T": T}
    if subshape is not None:
        p["SubShape"] = list(subshape)
    if two_sz is not None:
        p["2Sz"] = two_sz
    p.update(param or {})
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U}} if U else None
    return write_input(dest, cellshape=(L, L, 1), transfer=square_transfer(),
                       interactions=inter, param=p,
                       mode=None if fock else {"flag_fock": False})


def spin_transfer_text(norb, entries):
    """A 9-column (spin-resolved) Transfer file of ``{(R, a, b, s, t): v}``
    (0-based orbitals, spins 0 = up and 1 = down), Hermitian-closed: each
    entry gets its partner ``(-R, b, a, t, s): conj(v)``. Pass it as
    ``extra_files={"transfer.dat": ...}`` to :func:`write_input`."""
    closed = {}
    for (r, a, b, s, t), v in entries.items():
        closed[(r, a, b, s, t)] = closed.get((r, a, b, s, t), 0) + complex(v)
        neg = tuple(-x for x in r)
        if (neg, b, a, t, s) != (r, a, b, s, t):
            closed[(neg, b, a, t, s)] = (closed.get((neg, b, a, t, s), 0)
                                         + complex(v).conjugate())
    lines = ["Transfer in wannier90-like format for uhfk w90spin/1",
             str(norb), "1", " 1"]
    for (r, a, b, s, t), v in closed.items():
        lines.append("{} {} {} {} {} {} {} {:.15g} {:.15g}".format(
            *r, a + 1, b + 1, s, t, v.real, v.imag))
    return "\n".join(lines) + "\n"


def spin_independent_entries(transfer):
    """The 9-column entries ``(R, a, b, s, s)`` of a spin-independent
    7-column table, one direction of each Hermitian pair
    (:func:`spin_transfer_text` adds the partners)."""
    out = {}
    for (r, a, b), v in transfer.items():
        neg = tuple(-x for x in r)
        if (neg, b, a, 0, 0) in out:
            continue
        for s in range(2):
            out[(r, a, b, s, s)] = v
    return out


#: the triangular Bravais vectors a1 = (1, 0, 0), a2 = (1/2, sqrt(3)/2, 0)
TRIANGULAR_RVEC = ((1.0, 0.0, 0.0), (0.5, 3.0 ** 0.5 / 2.0, 0.0), (0.0, 0.0, 1.0))


def triangular_transfer(t=1.0):
    """Nearest-neighbour hopping ``t`` of one orbital on the triangular
    lattice: the bonds a1, a2 and a2 - a1 (reduced (1,0), (0,1), (-1,1))
    and their negatives."""
    out = {}
    for r in ((1, 0, 0), (0, 1, 0), (-1, 1, 0)):
        out[(r, 0, 0)] = -t
        out[(tuple(-x for x in r), 0, 0)] = -t
    return out


def triangular(dest, L, *, U=0.0, T=0.0, ncond=None, subshape=None, param=None):
    """One-orbital triangular lattice ``L x L`` (:data:`TRIANGULAR_RVEC`),
    hopping 1, on-site ``U``; half filling unless ``ncond`` is given."""
    p = {"Ncond": L * L if ncond is None else ncond, "T": T}
    if subshape is not None:
        p["SubShape"] = list(subshape)
    p.update(param or {})
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U}} if U else None
    return write_input(dest, cellshape=(L, L, 1), transfer=triangular_transfer(),
                       rvec=TRIANGULAR_RVEC, interactions=inter, param=p)


def chain_transfer(t=1.0, norb=1):
    """Nearest-neighbour hopping ``t`` along x of each of ``norb``
    orbitals (orbital-diagonal)."""
    out = {}
    for a in range(norb):
        out[((1, 0, 0), a, a)] = -t
        out[((-1, 0, 0), a, a)] = -t
    return out


def chain(dest, L, *, U=0.0, T=0.0, ncond=None, subshape=None, param=None, mode=None):
    """One-orbital chain of ``L`` sites along x, hopping 1, on-site ``U``;
    half filling unless ``ncond`` is given."""
    p = {"Ncond": L if ncond is None else ncond, "T": T}
    if subshape is not None:
        p["SubShape"] = list(subshape)
    p.update(param or {})
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U}} if U else None
    return write_input(dest, cellshape=(L, 1, 1), transfer=chain_transfer(),
                       interactions=inter, param=p, mode=mode)


def two_orbital(dest, L, *, delta, t=1.0, U=(0.0, 0.0), V=0.0, T=0.0, ncond=None,
                subshape=None, param=None, mode=None):
    """Two orbitals per site on a chain of ``L`` sites: hopping ``t`` along x
    within each orbital, on-site levels ``-delta/2`` (orbital 0) and
    ``+delta/2`` (orbital 1), on-site ``U = (U_0, U_1)`` and the on-site
    inter-orbital density interaction ``V`` (CoulombInter at R = 0); one
    electron per site unless ``ncond`` is given."""
    transfer = chain_transfer(t, norb=2)
    transfer[((0, 0, 0), 0, 0)] = -0.5 * delta
    transfer[((0, 0, 0), 1, 1)] = 0.5 * delta
    inter = {}
    intra = {((0, 0, 0), a, a): u for a, u in enumerate(U) if u}
    if intra:
        inter["CoulombIntra"] = intra
    if V:
        inter["CoulombInter"] = {((0, 0, 0), 0, 1): V, ((0, 0, 0), 1, 0): V}
    p = {"Ncond": L if ncond is None else ncond, "T": T}
    if subshape is not None:
        p["SubShape"] = list(subshape)
    p.update(param or {})
    return write_input(dest, cellshape=(L, 1, 1), norb=2, transfer=transfer,
                       interactions=inter or None, param=p, mode=mode)


def so_chain(dest, L, *, U=0.0, T=0.0, soc=0.0, ncond=None, param=None):
    """One physical orbital on a chain of ``L`` sites in spin-orbital mode
    (Geometry norb = 2, spin-orbitals ``2a + s``): hopping 1 along x and a
    spin-orbit hopping ``i soc sigma^z`` on the x bond (time-reversal
    invariant), on-site ``U``; half filling unless ``ncond`` is given."""
    transfer = {}
    for s in range(2):
        transfer[((1, 0, 0), s, s)] = -1.0 + 1j * soc * (1 - 2 * s)
        transfer[((-1, 0, 0), s, s)] = -1.0 - 1j * soc * (1 - 2 * s)
    p = {"Ncond": L if ncond is None else ncond, "T": T}
    p.update(param or {})
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U}} if U else None
    return write_input(dest, cellshape=(L, 1, 1), norb=2, transfer=transfer,
                       interactions=inter, param=p, mode={"enable_spin_orbital": True},
                       centers=[(0.0, 0.0, 0.0)] * 2)


#: one coupling table per two-body type on two orbitals (0-based), bonds
#: in x and y and on-site inter-orbital entries; PairHop is complex. The
#: tables reach every Hartree, Fock and PairHop slot form of the solver.
#: One direction of each pair is listed; :func:`hermitian_closed` adds the
#: partner ``(-R, b, a): conj(v)`` the reader requires.
GOLDEN_TABLES = {
    "CoulombIntra": {((0, 0, 0), 0, 0): 2.0, ((0, 0, 0), 1, 1): 1.5},
    "CoulombInter": {((0, 0, 0), 0, 1): 0.7, ((1, 0, 0), 0, 0): 0.3,
                     ((1, 0, 0), 0, 1): 0.2, ((0, 1, 0), 1, 1): 0.25,
                     ((1, 1, 0), 1, 0): 0.15},
    "Hund": {((0, 0, 0), 0, 1): 0.3, ((1, 0, 0), 0, 1): 0.1},
    "Ising": {((0, 0, 0), 0, 1): 0.2, ((0, 1, 0), 0, 0): 0.05},
    "PairLift": {((0, 0, 0), 0, 1): 0.1},
    "Exchange": {((0, 0, 0), 0, 1): 0.25, ((1, 0, 0), 0, 0): 0.07},
    "PairHop": {((0, 0, 0), 0, 1): 0.15 + 0.05j, ((2, 0, 0), 1, 0): 0.04 - 0.03j},
}


def hermitian_closed(rows):
    """``rows`` ``{(R, a, b): v}`` with every partner ``(-R, b, a): conj(v)``
    added (an entry that is its own partner is kept once)."""
    out = dict(rows)
    for (r, a, b), v in rows.items():
        key = (tuple(-x for x in r), b, a)
        if key not in out:
            out[key] = complex(v).conjugate()
    return out


def golden_transfer(*, inter_orbital=True):
    """Two orbitals on a square lattice: hoppings 1 (orbital 0) and 0.8
    (orbital 1) along x and y, levels -0.2 and +0.2, and (with
    ``inter_orbital``) an on-site mixing 0.1 and a complex x-bond mixing
    ``0.2 + 0.1i``, Hermitian-closed."""
    out = {}
    for a, t in ((0, 1.0), (1, 0.8)):
        for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
            out[(r, a, a)] = -t
    out[((0, 0, 0), 0, 0)] = -0.2
    out[((0, 0, 0), 1, 1)] = 0.2
    if inter_orbital:
        out[((0, 0, 0), 0, 1)] = 0.1
        out[((0, 0, 0), 1, 0)] = 0.1
        out[((1, 0, 0), 0, 1)] = 0.2 + 0.1j
        out[((-1, 0, 0), 1, 0)] = 0.2 - 0.1j
    return out


def golden(dest, *, cellshape=(6, 3, 1), types=None, fock=True, two_sz=None,
           spin_orbital=False, inter_orbital=True, tables=None, T=0.1, param=None):
    """Two orbitals on a ``cellshape`` square lattice (:func:`golden_transfer`)
    with the two-body tables ``types`` of :data:`GOLDEN_TABLES` (all by
    default; ``tables`` replaces them), ``SubShape = [1, 1, 1]`` unless
    ``param`` sets it. ``spin_orbital``: the same model in
    spin-orbital mode (Geometry norb = 4, spin-orbitals ``2a + s``) with an
    extra spin-orbit x-bond hopping ``0.05i sigma^z`` between the orbitals."""
    if tables is None:
        names = tuple(GOLDEN_TABLES) if types is None else tuple(types)
        tables = {name: GOLDEN_TABLES[name] for name in names}
    tables = {name: hermitian_closed(rows) for name, rows in tables.items()}
    transfer = golden_transfer(inter_orbital=inter_orbital)
    norb, centers, mode = 2, [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)], {}
    if spin_orbital:
        so = {}
        for (r, a, b), v in transfer.items():
            for s in range(2):
                so[(r, 2 * a + s, 2 * b + s)] = v
        if inter_orbital:
            for s in range(2):
                sz = 1 - 2 * s
                so[((1, 0, 0), 2 * 0 + s, 2 * 1 + s)] += 0.05j * sz
                so[((-1, 0, 0), 2 * 1 + s, 2 * 0 + s)] -= 0.05j * sz
        transfer, norb = so, 4
        centers = [(0.0, 0.0, 0.0)] * 2 + [(0.5, 0.0, 0.0)] * 2
        mode["enable_spin_orbital"] = True
    if not fock:
        mode["flag_fock"] = False
    nvol = cellshape[0] * cellshape[1] * cellshape[2]
    p = {"Ncond": nvol, "T": T, "SubShape": [1, 1, 1]}
    if two_sz is not None:
        p["2Sz"] = two_sz
    p.update(param or {})
    return write_input(dest, cellshape=cellshape, norb=norb, transfer=transfer,
                       centers=centers, interactions=tables, param=p, mode=mode or None)


#: the Kanamori couplings of :func:`wannier_like` on its two interacting
#: orbitals (``U' = U - 2J`` is derived)
WANNIER_LIKE_TYPES = ("CoulombIntra", "CoulombInter", "Hund", "Exchange", "PairHop")


def wannier_like_transfer():
    """The one-body table of :func:`wannier_like` (0-based orbitals, physical):
    two orbitals 0 and 1 at the cell origin and a ligand 2 at (1/2, 1/2) on
    the square lattice. Orbital 0 hops 1 and orbital 1 hops 0.7 to the
    nearest neighbours, with a real nearest-neighbour mixing 0.1 and a
    complex x-bond mixing ``0.1 + 0.08i`` (Hermitian-closed); levels 0, 0.2
    and -1.5; the ligand hybridises with the four surrounding cells, 0.6
    with orbital 0 and +-0.4 (a d-like sign pattern) with orbital 1."""
    t = {}
    for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
        t[(r, 0, 0)] = -1.0
        t[(r, 1, 1)] = -0.7
        t[(r, 0, 1)] = 0.1
        t[(r, 1, 0)] = 0.1
    t[((1, 0, 0), 0, 1)] = 0.1 + 0.08j
    t[((-1, 0, 0), 1, 0)] = 0.1 - 0.08j
    t[((0, 0, 0), 1, 1)] = 0.2
    t[((0, 0, 0), 2, 2)] = -1.5
    for r, s in (((0, 0, 0), 1), ((-1, 0, 0), -1), ((0, -1, 0), -1), ((-1, -1, 0), 1)):
        neg = tuple(-x for x in r)
        t[(r, 0, 2)] = 0.6
        t[(neg, 2, 0)] = 0.6
        t[(r, 1, 2)] = 0.4 * s
        t[(neg, 2, 1)] = 0.4 * s
    return t


def wannier_like(dest, *, U=4.0, J=0.6, T=0.2, L=4, ncond=40, fock=True, types=None,
                 spin_orbital=False, subshape=(1, 1, 1), param=None):
    """A Wannier-like multi-orbital model on an ``L x L`` square lattice: the
    two interacting orbitals and the ligand of :func:`wannier_like_transfer`
    (Wannier90-format files), the on-site Kanamori couplings ``U``
    (CoulombIntra), ``U' = U - 2J`` (CoulombInter), ``J`` (Hund, Exchange
    and PairHop) on orbitals 0 and 1 (``types``: a subset of
    :data:`WANNIER_LIKE_TYPES`), ``Ncond = ncond``. ``spin_orbital``: the
    same model in spin-orbital mode (spin-orbitals ``2a + s``, interaction
    files in physical orbitals)."""
    t = wannier_like_transfer()
    Up = U - 2.0 * J
    inter = {"CoulombIntra": {((0, 0, 0), 0, 0): U, ((0, 0, 0), 1, 1): U},
             "CoulombInter": {((0, 0, 0), 0, 1): Up, ((0, 0, 0), 1, 0): Up},
             "Hund": {((0, 0, 0), 0, 1): J, ((0, 0, 0), 1, 0): J},
             "Exchange": {((0, 0, 0), 0, 1): J, ((0, 0, 0), 1, 0): J},
             "PairHop": {((0, 0, 0), 0, 1): J, ((0, 0, 0), 1, 0): J}}
    if types is not None:
        inter = {k: v for k, v in inter.items() if k in types}
    centers = [(0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.5, 0.5, 0.0)]
    norb = 3
    mode = {} if fock else {"flag_fock": False}
    if spin_orbital:
        t = {(r, 2 * a + s, 2 * b + s): v for (r, a, b), v in t.items() for s in range(2)}
        centers = [c for c in centers for _ in range(2)]
        norb = 6
        mode["enable_spin_orbital"] = True
    p = {"Ncond": ncond, "T": T, "SubShape": list(subshape), "IterationMax": 3000,
         "EPS": 10, "Mix": 0.5}
    p.update(param or {})
    return write_input(dest, cellshape=(L, L, 1), norb=norb, transfer=t, centers=centers,
                       interactions=inter, param=p, mode=mode or None)


def so_two_orbital(dest, soc, *, U=3.5, T=0.2, L=4, tp=0.3, subshape=(2, 2, 1), param=None,
                   finput=None, foutput=None):
    """Two orbitals per site on an ``L x L`` square lattice in spin-orbital
    mode (spin-orbitals ``2a + s``): nearest-neighbour hopping 1 (orbital 0)
    and 0.8 (orbital 1), next-nearest ``tp`` times that, levels -0.1 and
    +0.1, the on-site inter-orbital spin-orbit coupling ``i soc sigma^z``
    (``c^dag_0 c_1``, time-reversal invariant), on-site ``U`` on both
    orbitals, half filling (two electrons per site)."""
    t = {}
    for a, hop in ((0, 1.0), (1, 0.8)):
        for s in range(2):
            i = 2 * a + s
            for r in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)):
                t[(r, i, i)] = -hop
            if tp:
                for r in ((1, 1, 0), (-1, -1, 0), (1, -1, 0), (-1, 1, 0)):
                    t[(r, i, i)] = -tp * hop
    for s in range(2):
        sz = 1 - 2 * s
        t[((0, 0, 0), s, s)] = -0.1
        t[((0, 0, 0), 2 + s, 2 + s)] = 0.1
        t[((0, 0, 0), s, 2 + s)] = 1j * soc * sz
        t[((0, 0, 0), 2 + s, s)] = -1j * soc * sz
    p = {"Ncond": 2 * L * L, "T": T, "SubShape": list(subshape), "IterationMax": 20000,
         "EPS": 12, "Mix": 0.5}
    p.update(param or {})
    return write_input(dest, cellshape=(L, L, 1), norb=4, transfer=t,
                       centers=[(0.0, 0.0, 0.0)] * 4,
                       interactions={"CoulombIntra": {((0, 0, 0), 0, 0): U,
                                                      ((0, 0, 0), 1, 1): U}},
                       param=p, mode={"enable_spin_orbital": True}, finput=finput,
                       foutput=foutput)
