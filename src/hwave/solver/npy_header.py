"""Reading an ``.npy`` array header without loading the array.

Two call sites need the shape of an array stored in an ``.npz`` member
before deciding whether to materialise it: the Eliashberg frequency probe
and the bond preflight's Green-function peek. Both used private
``numpy.lib.format`` helpers (``_check_version``/``_read_array_header``),
which numpy removed in 2.4 -- so the probe stopped working on any file,
not merely on the newer header versions the fallback existed for. The
dispatch lives here once so the two sites cannot answer the question
differently.

Only the PUBLIC ``numpy.lib.format`` surface is used. numpy exposes a
reader per header version (``read_array_header_1_0``,
``read_array_header_2_0``) but has never published one for version 3.0,
even though it writes 3.0 whenever a dtype carries field names outside
latin-1. Version 3.0's header is byte-compatible with 2.0 -- both prefix
the header dict with a 4-byte little-endian length, and they differ only
in the encoding of that dict (utf-8 vs latin-1), which numpy's own reader
decodes leniently. Reading a 3.0 header with the 2.0 reader is therefore
correct, and is what numpy's removed private dispatcher did internally.
"""

import numpy as np


class UnsupportedNpyHeaderVersion(Exception):
    """The ``.npy`` header version has no reader available here.

    Raised only for versions this module has never seen -- a future numpy
    header revision. Callers decide whether that is fatal or a reason to
    fall back to a full load.
    """


#: Header versions handled, mapped to the public reader that parses them.
#: 3.0 is read with the 2.0 reader (see the module docstring).
_READERS = {
    (1, 0): "read_array_header_1_0",
    (2, 0): "read_array_header_2_0",
    (3, 0): "read_array_header_2_0",
}


def read_npy_header_shape(fh):
    """The ``shape`` of the ``.npy`` stream ``fh``, read from its header.

    ``fh`` must be positioned at the start of the member; on return it sits
    just past the header, before the array data (which is never read).

    Parameters
    ----------
    fh : file-like
        A binary stream open on an ``.npy`` member.

    Returns
    -------
    tuple of int
        The stored array's shape.

    Raises
    ------
    UnsupportedNpyHeaderVersion
        If the header version is not one this module knows how to read.
    """
    version = np.lib.format.read_magic(fh)
    reader_name = _READERS.get(tuple(version))
    if reader_name is None:
        raise UnsupportedNpyHeaderVersion(
            "no reader available for NPY format version {!r}".format(
                tuple(version)))
    reader = getattr(np.lib.format, reader_name, None)
    if reader is None:
        # The public surface itself changed under us.
        raise UnsupportedNpyHeaderVersion(
            "numpy.lib.format.{} is unavailable (numpy {})".format(
                reader_name, np.__version__))
    shape, _fortran_order, _dtype = reader(fh)
    return shape
