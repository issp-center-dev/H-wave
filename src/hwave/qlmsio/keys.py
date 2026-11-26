"""Canonical spelling of the ``[file.input.interaction]`` keys.

The interaction keys are case-insensitive: ``geometry``, ``GEOMETRY`` and
``Geometry`` name the same file. The readers used to accept any spelling but
store the table under the spelling the user wrote, so a solver that compared
a key with its canonical spelling (``k == "Geometry"``) mis-handled the table
(UHFk crashed in ``_check_cellsize``). The keys are therefore rewritten to
their canonical spelling ONCE, at the input boundary, by
:func:`canonicalize_interaction_keys`; everything downstream sees canonical
names only.
"""

#: Keys of the k-space reader (UHFk, RPA, FLEX, hwave_sc, hwave_dos).
KSPACE_KEYS = ("path_to_input", "Geometry", "Transfer", "Coulomb",
               "CoulombIntra", "CoulombInter", "Hund", "Ising", "PairLift",
               "Exchange", "PairHop", "Extern")

#: Keys of the real-space reader (UHFr).
REALSPACE_KEYS = ("path_to_input", "Trans", "CoulombIntra", "CoulombInter",
                  "Hund", "Ising", "PairLift", "Exchange", "PairHop",
                  "InterAll")

#: Every key either reader accepts (no two differ only in case).
ALL_KEYS = tuple(dict.fromkeys(KSPACE_KEYS + REALSPACE_KEYS))

_LOWER_TO_CANONICAL = {k.lower(): k for k in ALL_KEYS}
assert len(_LOWER_TO_CANONICAL) == len(ALL_KEYS)


def canonical_key(key):
    """Canonical spelling of an interaction key, or ``key`` itself when it
    is not a known key (the caller reports unknown keys)."""
    return _LOWER_TO_CANONICAL.get(str(key).lower(), key)


def canonicalize_interaction_keys(table,
                                  section="[file.input.interaction]"):
    """Return ``table`` as a plain dict whose keys use the canonical spelling.

    Values and key order are kept. Unknown keys are passed through unchanged
    so that each reader keeps its own "unknown keyword" message.

    Raises
    ------
    ValueError
        When two keys differ only in case (for example ``Geometry`` and
        ``geometry``): TOML accepts both, but they name the same entry and
        one of them would be silently ignored.
    """
    table = table or {}
    spellings = {}
    for key in table:
        spellings.setdefault(str(key).lower(), []).append(key)
    for group in spellings.values():
        if len(group) > 1:
            raise ValueError(
                "{} {} is given more than once with different spellings "
                "({}); keys are case-insensitive.".format(
                    section, canonical_key(group[0]),
                    ", ".join(repr(k) for k in group)))
    return {canonical_key(key): value for key, value in table.items()}
