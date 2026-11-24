"""Deep read-only copies for the result objects of :mod:`hwave.results`."""
from collections.abc import Mapping
from types import MappingProxyType

import numpy as np


def _freeze(obj):
    """A deep private copy of ``obj`` in which every ndarray is a fresh copy
    marked ``flags.writeable = False``.

    Mappings and lists are rebuilt (new containers: a CaseInsensitiveDict
    stays one, any other mapping becomes a dict), tuples stay tuples, and
    anything else (scalars, strings, frozen dataclasses) is returned as is.
    The input is never modified."""
    if isinstance(obj, np.ndarray):
        arr = np.array(obj, copy=True)
        if arr.dtype == object:
            for index in np.ndindex(arr.shape):
                arr[index] = _freeze(arr[index])
        arr.flags.writeable = False
        return arr
    if isinstance(obj, Mapping):
        items = {key: _freeze(value) for key, value in obj.items()}
        from requests.structures import CaseInsensitiveDict
        if isinstance(obj, CaseInsensitiveDict):
            return CaseInsensitiveDict(items)
        return items
    if isinstance(obj, list):
        return [_freeze(value) for value in obj]
    if isinstance(obj, tuple):
        values = [_freeze(value) for value in obj]
        if hasattr(obj, "_fields"):           # a namedtuple keeps its type
            return type(obj)(*values)
        return tuple(values)
    return obj


def _read_only(obj):
    """A deep read-only view of the mappings and lists of ``obj`` (a parsed
    TOML table): every mapping is wrapped in ``types.MappingProxyType`` (a
    CaseInsensitiveDict keeps its case rules behind the proxy) and every
    list becomes a tuple. Nothing reachable from the result can be changed
    in place, so a caller cannot alter what later readers see."""
    if isinstance(obj, Mapping):
        items = {key: _read_only(value) for key, value in obj.items()}
        from requests.structures import CaseInsensitiveDict
        if isinstance(obj, CaseInsensitiveDict):
            return MappingProxyType(CaseInsensitiveDict(items))
        return MappingProxyType(items)
    if isinstance(obj, (list, tuple)):
        return tuple(_read_only(value) for value in obj)
    return obj
