"""Deep read-only copies for the result objects of :mod:`hwave.results`."""
from collections.abc import Mapping

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
