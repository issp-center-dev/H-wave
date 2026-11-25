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
                extra_files=None):
    """Write a UHFk input under ``dest`` and return the TOML path.

    ``transfer`` and each ``interactions[name]`` are ``{(R, a, b): value}``;
    ``param`` is added to ``[mode.param]`` (``CellShape`` is set from
    ``cellshape``), ``mode`` to ``[mode]``, ``finput`` to ``[file.input]``.
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
    doc = {"log": {"print_level": 0}, "mode": m,
           "file": {"input": fin, "output": {"path_to_output": "output"}}}
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
