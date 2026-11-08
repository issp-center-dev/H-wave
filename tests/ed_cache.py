"""Content-keyed sharing of exact-diagonalization results across tests.

Several ED oracle modules diagonalize the SAME Hamiltonian more than once:
the free fixture reached from two different couplings at zero amplitude,
the Hartree-Fock-only reference of two granules on one fixture, or the
case-M off-site CoulombInter ray that ``test_bond_longitudinal_ed`` and
``test_bond_vs_ed_oracle`` both need. :func:`sector` and :func:`total` hand
out a stand-in for :class:`tests.ed_oracle_util.SectorED` /
:class:`~tests.ed_oracle_util.TotalNED` whose public observables are
computed once per process for each distinct input and then shared.

What makes the sharing exact
----------------------------
* The key is the full CONTENT of the inputs, never an object identity:
  the fixture's ``L``, ``norb``, hopping table, on-site energies, ``T`` and
  ``mu``; the one-body matrix actually diagonalized (``fx.build_h1()`` when
  none is given) as raw bytes; the interaction term list in its given
  ORDER, coefficients as raw bytes (the order fixes the floating-point
  summation of the sector Hamiltonian); the engine class; the observable's
  name and its arguments (arrays as dtype, shape and bytes). Two calls
  share a result only when every one of these is bit-identical, so a
  shared result is the value the second call would have computed.
* The one normalisation: a term whose coefficient is exactly zero is
  dropped (from the key AND from what the engine is given). Such a term
  adds a signed zero to sector-Hamiltonian entries that are never ``-0``
  (they start at ``+0`` and only ever receive sums), which leaves every
  entry bit-identical, so the engine's input matrix -- and hence its
  output -- is unchanged. This is what lets the zero-coupling point of a
  coupling sweep share the free diagonalization
  (``test_ed_cache`` checks the bit identity directly).
* Results are returned READ-ONLY (``writeable = False``), so one caller
  cannot alter what the next one receives; an in-place update raises.
* The engines themselves are kept in a small LRU (they hold the sector
  eigenvectors, up to ~100 MB for the largest total-N fixture), so two
  observables of one Hamiltonian share one diagonalization without the
  process accumulating every ED it ever built; their spectra are frozen
  read-only as well.

The cache lives in this module so that every test module importing it
shares one copy when they run in the same process (one ``pytest`` or
``unittest`` invocation); separate processes simply recompute.
"""

import collections

import numpy as np

from tests import ed_oracle_util

#: observables served through the cache (the public correlators)
_OBSERVABLES = frozenset({
    "bond_correlator", "bond_correlator_transverse", "pair_correlator",
    "chi_connected", "green", "density_matrix",
})

_ENGINES = {"SectorED": ed_oracle_util.SectorED,
            "TotalNED": ed_oracle_util.TotalNED}

#: engines kept alive at once (see the module docstring)
ENGINE_LRU_SIZE = 2

_results = {}
_engines = collections.OrderedDict()
_stats = collections.Counter()


def _array_key(a):
    a = np.ascontiguousarray(a)
    return ("ndarray", a.dtype.str, a.shape, a.tobytes())


def _canon(x):
    """A hashable, exact (bitwise for floats) form of an argument."""
    if isinstance(x, np.ndarray):
        return _array_key(x)
    if isinstance(x, (list, tuple)):
        return (type(x).__name__,) + tuple(_canon(y) for y in x)
    if isinstance(x, dict):
        return ("dict",) + tuple(sorted((_canon(k), _canon(v)) for k, v in x.items()))
    if isinstance(x, (bool, int, np.integer)):
        return ("int", int(x))
    if isinstance(x, (float, complex, np.floating, np.complexfloating)):
        return ("num", np.asarray(x, dtype=complex).tobytes())
    if x is None or isinstance(x, str):
        return x
    raise TypeError("ed_cache: no exact key for argument of type {}".format(
        type(x).__name__))


def _fixture_key(fx):
    return ("EDFixture", int(fx.L), int(fx.norb), _canon(dict(fx.t)),
            _canon(tuple(fx.eps)), _canon(float(fx.T)), _canon(float(fx.mu)))


def _freeze(x):
    if isinstance(x, np.ndarray):
        x = np.array(x, copy=True)
        x.flags.writeable = False
        return x
    if isinstance(x, tuple):
        return tuple(_freeze(y) for y in x)
    return x


class CachedED:
    """Stand-in for one ``SectorED``/``TotalNED`` construction: the public
    observables are served from the process-wide cache, the engine is
    built only on a miss (or by :meth:`engine`)."""

    def __init__(self, engine_name, fx, terms=(), h1=None):
        self.fx = fx
        self._name = engine_name
        self._terms = [tuple(t) for t in terms if t[-1] != 0]
        self._h1 = fx.build_h1() if h1 is None else np.asarray(h1)
        self.key = (engine_name, _fixture_key(fx),
                    _canon(tuple(self._terms)),
                    _array_key(np.asarray(self._h1)))

    def engine(self):
        """The real engine (shared through the small LRU)."""
        eng = _engines.get(self.key)
        if eng is not None:
            _engines.move_to_end(self.key)
            _stats["engine_hit"] += 1
            return eng
        _stats["engine_miss"] += 1
        eng = _ENGINES[self._name](self.fx, terms=self._terms, h1=self._h1)
        for name in ("_ev", "_w", "_V"):
            for arr in getattr(eng, name).values():
                arr.flags.writeable = False
        _engines[self.key] = eng
        while len(_engines) > ENGINE_LRU_SIZE:
            _engines.popitem(last=False)
        return eng

    def __getattr__(self, name):
        if name not in _OBSERVABLES:
            raise AttributeError(
                "ed_cache.CachedED serves only {}; call .engine() for "
                "{!r}".format(sorted(_OBSERVABLES), name))

        def observable(*args, **kwargs):
            k = (self.key, name, _canon(args), _canon(kwargs))
            hit = _results.get(k)
            if hit is not None:
                _stats["result_hit"] += 1
                return hit
            _stats["result_miss"] += 1
            out = _freeze(getattr(self.engine(), name)(*args, **kwargs))
            _results[k] = out
            return out

        return observable


def sector(fx, terms=(), h1=None):
    """Cached stand-in for ``ed_oracle_util.SectorED(fx, terms, h1)``."""
    return CachedED("SectorED", fx, terms, h1)


def total(fx, terms=(), h1=None):
    """Cached stand-in for ``ed_oracle_util.TotalNED(fx, terms, h1)``."""
    return CachedED("TotalNED", fx, terms, h1)


def stats():
    """Hit/miss counters (a copy), for diagnostics."""
    return dict(_stats)
