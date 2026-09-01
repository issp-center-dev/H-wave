"""Reading an ``.npy`` array header without loading the array.

Two call sites need the shape of an array stored in an ``.npz`` member
before deciding whether to materialise it: the Eliashberg frequency probe
and the bond preflight's Green-function peek. Both used private
``numpy.lib.format`` helpers (``_check_version``/``_read_array_header``),
which numpy removed in 2.4 -- so the probe stopped working on any file,
not merely on the newer header versions the fallback existed for. The
dispatch lives here once so the two sites cannot answer the question
differently.

No private numpy API is used. numpy publishes a reader per header
version (``read_array_header_1_0``, ``read_array_header_2_0``) but has
never published one for version 3.0, even though it writes 3.0 whenever
a dtype carries field names outside latin-1.

Version 3.0 is NOT delegated to the 2.0 reader. The two share the
4-byte little-endian header length, but the 2.0 reader decodes the
header dict as latin-1 and applies numpy's Python-2 compatibility
filtering -- so a malformed 3.0 header carrying Python-2 long literals
(``'shape': (1L, 24L)``) is ACCEPTED by it and returns a shape, while
``np.load`` rejects the same file. Both callers here treat a readable
header as licence to skip or precede the authoritative load, so
accepting what the loader will reject turns a clear load error into a
wrong frequency count or an unrelated preflight failure. Version 3.0 is
therefore parsed here directly from the format specification: strict
utf-8, a literal dict, and a shape of plain integers.
"""

import numpy as np


class UnsupportedNpyHeaderVersion(Exception):
    """The ``.npy`` header version has no reader available here.

    Raised only for versions this module has never seen -- a future numpy
    header revision. Callers decide whether that is fatal or a reason to
    fall back to a full load.
    """


#: Header versions delegated to numpy's own published readers.
_READERS = {
    (1, 0): "read_array_header_1_0",
    (2, 0): "read_array_header_2_0",
}

#: Header versions parsed here (see the module docstring for why 3.0 is
#: not delegated).
_LOCAL_VERSIONS = ((3, 0),)


#: Upper bound on a header's DECODED length, the same limit and the same
#: quantity numpy applies (its ``max_header_size`` bounds the decoded
#: string, not the encoded bytes).
_MAX_HEADER_CHARS = 10000

#: Upper bound on the DECLARED byte length, used only to refuse an
#: absurd read before it happens. The declared length is four bytes
#: taken from the file, so without this a crafted member could make a
#: probe -- whose entire purpose is to avoid reading data -- allocate up
#: to 4 GiB and hand it to the literal parser. It is deliberately the
#: widest encoding of the character limit (utf-8 uses at most 4 bytes
#: per character) rather than the character limit itself: version 3.0
#: exists FOR field names outside latin-1, and a header of 10,000 CJK
#: characters is ~30,000 bytes yet perfectly loadable, so a byte-based
#: limit at 10,000 would reject files np.load accepts.
_MAX_HEADER_LEN = 4 * _MAX_HEADER_CHARS

#: The keys an NPY header dict must carry, exactly (numpy rejects both
#: missing and extra keys).
_REQUIRED_KEYS = frozenset(("descr", "fortran_order", "shape"))


def _read_header_3_0(fh):
    """Parse a version-3.0 header and return its ``shape``.

    Follows the NPY specification: a 4-byte little-endian length, then
    that many bytes of utf-8 holding a Python literal dict.

    The validation deliberately mirrors what ``np.load`` itself enforces
    -- exactly the three required keys, a boolean ``fortran_order``, a
    dtype descriptor numpy accepts, and a shape of non-negative plain
    integers. A probe that accepted more than the loader does would
    reintroduce the defect this parser exists to remove: reporting a
    shape for a file that cannot then be loaded, turning a clear load
    error into a wrong frequency count or an unrelated preflight
    failure.

    Every rejection raises ``ValueError`` -- including for headers whose
    keys cannot be ordered against each other -- which both callers
    already handle as "this header is unusable".
    """
    import ast

    raw_len = fh.read(4)
    if len(raw_len) != 4:
        raise ValueError("truncated NPY 3.0 header length field")
    header_len = int.from_bytes(raw_len, "little")
    if header_len > _MAX_HEADER_LEN:
        raise ValueError(
            "NPY 3.0 header declares {} bytes, above the {}-byte limit; "
            "refusing to read it".format(header_len, _MAX_HEADER_LEN))
    raw = fh.read(header_len)
    if len(raw) != header_len:
        raise ValueError("truncated NPY 3.0 header")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(
            "NPY 3.0 header is not valid utf-8: {}".format(exc))
    if len(text) > _MAX_HEADER_CHARS:
        raise ValueError(
            "NPY 3.0 header is {} characters, above the {}-character "
            "limit".format(len(text), _MAX_HEADER_CHARS))
    try:
        header = ast.literal_eval(text.strip())
    except (ValueError, SyntaxError, MemoryError, RecursionError) as exc:
        raise ValueError(
            "NPY 3.0 header is not a Python literal: {}".format(exc))
    if not isinstance(header, dict):
        raise ValueError(
            "NPY 3.0 header is {}, expected a dict".format(
                type(header).__name__))
    if _REQUIRED_KEYS != set(header):
        # Sort by repr: a header may carry keys of mixed types, and
        # plain sorted() would raise TypeError out of this function --
        # which sc's probe does not catch, so it would escape instead of
        # deferring to the authoritative loader.
        raise ValueError(
            "NPY 3.0 header keys are {}, expected exactly {}".format(
                sorted(map(repr, header)), sorted(_REQUIRED_KEYS)))
    if not isinstance(header["fortran_order"], bool):
        raise ValueError(
            "NPY 3.0 header fortran_order is not a bool: {!r}".format(
                header["fortran_order"]))
    try:
        np.lib.format.descr_to_dtype(header["descr"])
    except Exception as exc:
        raise ValueError(
            "NPY 3.0 header descr is not a valid dtype descriptor: "
            "{!r} ({})".format(header["descr"], exc))
    shape = header["shape"]
    if (not isinstance(shape, tuple)
            or not all(isinstance(n, int) and not isinstance(n, bool)
                       and n >= 0 for n in shape)):
        raise ValueError(
            "NPY 3.0 header shape is not a tuple of non-negative "
            "integers: {!r}".format(shape))
    return shape


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
    version = tuple(np.lib.format.read_magic(fh))
    if version in _LOCAL_VERSIONS:
        return _read_header_3_0(fh)
    reader_name = _READERS.get(version)
    if reader_name is None:
        raise UnsupportedNpyHeaderVersion(
            "no reader available for NPY format version {!r}".format(
                version))
    reader = getattr(np.lib.format, reader_name, None)
    if reader is None:
        # The public surface itself changed under us.
        raise UnsupportedNpyHeaderVersion(
            "numpy.lib.format.{} is unavailable (numpy {})".format(
                reader_name, np.__version__))
    shape, _fortran_order, _dtype = reader(fh)
    return shape
