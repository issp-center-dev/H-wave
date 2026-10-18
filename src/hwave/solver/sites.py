"""Which orbitals share a site (GitHub issue #234).

The pair-space vertex of the RPA ring, of the FLEX general path and of the
non-bond ``hwave_sc`` path gives an ``R = 0`` two-body row the on-site
(Kanamori) slot map, the exchange crossing included, while an off-site row
(``R != 0``) contributes its density (Hartree) part only. An ``R = 0`` row
between two orbitals that sit on DIFFERENT sites of the cell (for example
the two molecules of a dimer) is an inter-site interaction like any
off-site one; giving it the exchange crossing would make the result depend
on the choice of unit cell.

``[mode.param] orbital_site`` assigns each orbital of the (pre-fold, and in
spin-orbital mode physical) orbital list to a site. An ``R = 0`` row
between orbitals of different sites is then treated exactly as an off-site
row. Without the key every orbital is taken to be on one site, which is the
behaviour before the key existed.

The bond-resolved paths (``longitudinal_bond_channels``) resolve the
exchange crossing of every bond exactly, the in-cell ones included, and
keep doing so whatever the key says. The Hartree-Fock mean field (UHF and
the FLEX Hartree-Fock self-energy) keeps the exchange of every row: it is
exact at that level.
"""
import logging

logger = logging.getLogger("qlms").getChild("sites")

#: the two-body types whose ``R = 0`` inter-orbital rows the key can affect
TWO_BODY_TYPES = ("CoulombInter", "Hund", "Ising", "Exchange", "PairLift", "PairHop",
                  "Coulomb")


def parse_orbital_site(raw, norb, source="[mode.param] orbital_site"):
    """``raw`` (``None`` or a list of non-negative integers, one per
    orbital) as a tuple of ``norb`` site labels, or ``None``."""
    if raw is None:
        return None
    if isinstance(raw, (str, bytes)) or not hasattr(raw, "__len__"):
        raise ValueError("{} must be a list with one site label per orbital, got {!r}"
                         .format(source, raw))
    labels = []
    for v in raw:
        if isinstance(v, bool) or not isinstance(v, int) or v < 0:
            raise ValueError("{} must hold non-negative integers, got {!r}".format(source, raw))
        labels.append(int(v))
    if len(labels) != int(norb):
        raise ValueError("{} has {} entries but the model has {} orbitals (one label "
                         "per orbital of the original cell{})".format(
                             source, len(labels), norb,
                             "; physical orbitals in spin-orbital mode"))
    return tuple(labels)


def is_onsite(irvec, orbvec, site):
    """``True`` when the row ``(irvec, orbvec)`` acts within one site: zero
    displacement and, if ``site`` is given, both orbitals on the same
    site."""
    if tuple(int(x) for x in irvec) != (0, 0, 0):
        return False
    if site is None:
        return True
    a, b = int(orbvec[0]), int(orbvec[1])
    return site[a] == site[b]


def cross_site_rows(tables, site):
    """``{type: [(irvec, orbvec), ...]}`` of the ``R = 0`` rows between
    orbitals of different sites (empty without ``site``)."""
    out = {}
    if site is None:
        return out
    for t in TWO_BODY_TYPES:
        tbl = tables.get(t) if hasattr(tables, "get") else None
        if not tbl:
            continue
        rows = [(irvec, orbvec) for (irvec, orbvec) in tbl
                if tuple(int(x) for x in irvec) == (0, 0, 0)
                and not is_onsite(irvec, orbvec, site)]
        if rows:
            out[t] = rows
    return out


def warn_if_undeclared(tables, site, solver):
    """Warn once when ``R = 0`` inter-orbital rows are present and
    ``orbital_site`` is not declared: they are then treated as interactions
    on one site (with their exchange crossing)."""
    if site is not None:
        return
    types = []
    for t in TWO_BODY_TYPES:
        tbl = tables.get(t) if hasattr(tables, "get") else None
        if not tbl:
            continue
        if any(tuple(int(x) for x in irvec) == (0, 0, 0) and orbvec[0] != orbvec[1]
               for (irvec, orbvec) in tbl):
            types.append(t)
    if types:
        logger.warning(
            "%s: the R = 0 inter-orbital rows of %s are treated as interactions "
            "between orbitals of the same site (Kanamori form, with their exchange "
            "crossing). If some of these orbitals sit on different sites of the cell, "
            "declare [mode.param] orbital_site (one site label per orbital) so that "
            "those rows are treated like off-site rows; declare all orbitals on one "
            "site to silence this message.", solver, ", ".join(types))
