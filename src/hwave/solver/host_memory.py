"""The host memory available to this process.

:func:`available_host_bytes` is the default of the spin-orbital bond memory
cap: the minimum over the REMAINING amounts of every source that gives a
value, each clamped at 0. The sources, in order (a tie goes to the first):

* ``"psutil"``: ``psutil.virtual_memory().available``;
* ``"meminfo"``: the ``MemAvailable:`` line of ``/proc/meminfo``, consulted
  only when psutil gives no value;
* ``"RLIMIT_AS"``: a finite soft ``RLIMIT_AS`` minus the process's
  ``VmSize`` (``/proc/self/status``, else psutil ``memory_info().vms``);
* ``"RLIMIT_DATA"``: a finite soft ``RLIMIT_DATA`` minus ``VmData``;
* ``"cgroup v2"`` / ``"cgroup v1"``: the remaining memory of the process's
  own cgroup and its ancestors (:func:`_cgroup_remaining`).

Every reader gives ``None`` ("no value") when its file is unreadable or
does not parse; nothing here raises. A cgroup or mountinfo file that
exists but cannot be read or parsed is named in a WARNING (logger
``qlms.solver.host_memory``). No value at all gives
``HostMemory(None, "none")``: the caller then disables its admission check
(there is no fixed fallback size).
"""

import logging
import os
import re
from typing import NamedTuple, Optional

__all__ = ["HostMemory", "available_host_bytes"]

#: the files read (module constants, so a test can point them elsewhere)
PROC_SELF_CGROUP = "/proc/self/cgroup"
PROC_SELF_MOUNTINFO = "/proc/self/mountinfo"
PROC_SELF_STATUS = "/proc/self/status"
PROC_MEMINFO = "/proc/meminfo"

#: a cgroup v1 memory limit at or above this is the kernel's "unlimited"
#: (``PAGE_COUNTER_MAX`` pages, e.g. ``9223372036854771712``)
V1_UNLIMITED_FLOOR = 2 ** 62

logger = logging.getLogger("qlms.solver.host_memory")

_OCTAL = re.compile(r"\\([0-7]{3})")
_DIGITS = re.compile(r"[0-9]+", re.ASCII)


class HostMemory(NamedTuple):
    """The available host memory: ``nbytes`` (``None`` when no source gave a
    value) and the ``source`` that set it (``"none"`` without a value)."""
    nbytes: Optional[int]
    source: str


class Mount(NamedTuple):
    """One line of ``/proc/self/mountinfo``: the mount's ``root`` within its
    filesystem, its ``mount_point``, ``fstype`` and ``super_options``
    (escapes decoded)."""
    root: str
    mount_point: str
    fstype: str
    super_options: str


def _read_text(path):
    with open(path) as f:
        return f.read()


def _reason(exc):
    """A short reason for a failed read: the OS error text, else the
    exception message."""
    if isinstance(exc, OSError):
        return exc.strerror or type(exc).__name__
    return str(exc)


def _read_cgroup_file(path):
    """The lines of a ``/proc`` cgroup file, or ``None``: a missing file
    silently, an unreadable one with a WARNING."""
    try:
        with open(path, "rb") as f:
            raw = f.read()
    except FileNotFoundError:
        return None
    except OSError as exc:
        logger.warning("%s could not be read (%s); the cgroup memory limit is not read",
                       path, _reason(exc))
        return None
    return raw.splitlines()


def _malformed_line(path, n, line):
    logger.warning("%s line %d is malformed: %r; the cgroup memory limit is not read",
                   path, n, line)


def _parse_kib_line(line, field):
    """The bytes of a ``<field>:<spaces><integer> kB`` line (``/proc/meminfo``,
    ``/proc/self/status``): ``integer * 1024``. Any other field name or unit,
    a missing unit, a non-integer or a negative value raises ``ValueError``."""
    name, sep, rest = line.partition(":")
    if not sep or name != field:
        raise ValueError("not a {} line: {!r}".format(field, line))
    parts = rest.split()
    if len(parts) != 2 or parts[1] != "kB" or not _DIGITS.fullmatch(parts[0]):
        raise ValueError("{} is not '<integer> kB': {!r}".format(field, line))
    return int(parts[0]) * 1024


def _read_kib_field(path, field):
    """The bytes of the ``field`` line of ``path``, or ``None``."""
    try:
        for line in _read_text(path).splitlines():
            if line.startswith(field + ":"):
                return _parse_kib_line(line, field)
    except (OSError, ValueError):
        return None
    return None


def _psutil():
    """``(psutil, errors)``: the module and the exceptions a read of it may
    raise (``psutil.Error``, e.g. ``AccessDenied``, when the module defines
    it), or ``(None, ())`` when psutil cannot be imported."""
    try:
        import psutil
    except ImportError:
        return None, ()
    errors = (OSError, ValueError, AttributeError, TypeError)
    base = getattr(psutil, "Error", None)
    if isinstance(base, type) and issubclass(base, BaseException):
        errors += (base,)
    return psutil, errors


def _read_psutil():
    """``psutil.virtual_memory().available`` in bytes, or ``None``."""
    psutil, errors = _psutil()
    if psutil is None:
        return None
    try:
        return int(psutil.virtual_memory().available)
    except errors:
        return None


def _read_psutil_vms():
    """The process's virtual size from psutil, or ``None``."""
    psutil, errors = _psutil()
    if psutil is None:
        return None
    try:
        return int(psutil.Process().memory_info().vms)
    except errors:
        return None


def _read_meminfo():
    """``MemAvailable`` of ``/proc/meminfo`` in bytes, or ``None``."""
    return _read_kib_field(PROC_MEMINFO, "MemAvailable")


def _read_status(field):
    """A ``kB`` field (``VmSize``, ``VmData``) of ``/proc/self/status`` in
    bytes, or ``None``."""
    return _read_kib_field(PROC_SELF_STATUS, field)


def _read_rlimit(name):
    """The finite soft limit ``resource.<name>`` in bytes, or ``None``
    (``RLIM_INFINITY``, a negative value, or no such limit)."""
    try:
        import resource
        soft, _ = resource.getrlimit(getattr(resource, name))
        soft = int(soft)
        if soft == resource.RLIM_INFINITY or soft < 0:
            return None
        return soft
    except (ImportError, AttributeError, OSError, ValueError, TypeError):
        return None


def _read_proc_self_cgroup():
    """The process's cgroup paths from ``/proc/self/cgroup``:
    ``{"v2": path, "v1 memory": path}`` (each key only when present).

    The v2 entry is hierarchy ``0`` with empty controllers (``0::/a/b``);
    the v1 memory entry is the one whose comma-separated controllers contain
    ``memory``. A line without exactly two ``:`` separators or with a path
    not starting with ``/`` makes the file malformed: ``{}``."""
    lines = _read_cgroup_file(PROC_SELF_CGROUP)
    if lines is None:
        return {}
    out = {}
    for n, raw in enumerate(lines, 1):
        line = raw.decode("utf-8", "surrogateescape")
        if not line.strip():
            continue
        parts = line.split(":")
        if len(parts) != 3 or not parts[2].startswith("/"):
            _malformed_line(PROC_SELF_CGROUP, n, line)
            return {}
        hid, controllers, path = parts
        if hid == "0" and controllers == "":
            out.setdefault("v2", path)
        elif "memory" in controllers.split(","):
            out.setdefault("v1 memory", path)
    return out


def _unescape(field):
    return _OCTAL.sub(lambda m: chr(int(m.group(1), 8)), field)


def _read_mountinfo():
    """The mounts of ``/proc/self/mountinfo``, in file order.

    Fields per ``proc(5)``: 4 ``root``, 5 ``mount point``, then optional
    fields up to the separator ``-``, then ``fstype``, ``source`` and
    ``super options``; octal escapes are decoded. A line without the
    separator or with fewer than 10 fields makes the file malformed: ``[]``
    (with a WARNING). A line that is not valid UTF-8 (a mount point with
    other bytes) is skipped with a WARNING; the other lines are kept."""
    lines = _read_cgroup_file(PROC_SELF_MOUNTINFO)
    if lines is None:
        return []
    out = []
    for n, raw in enumerate(lines, 1):
        try:
            line = raw.decode("utf-8")
        except UnicodeDecodeError:
            logger.warning("%s line %d is not valid UTF-8; that mount is ignored",
                           PROC_SELF_MOUNTINFO, n)
            continue
        if not line.strip():
            continue
        f = line.split()
        if len(f) < 10 or "-" not in f[6:] or len(f) < f.index("-", 6) + 4:
            _malformed_line(PROC_SELF_MOUNTINFO, n, line)
            return []
        sep = f.index("-", 6)
        out.append(Mount(_unescape(f[3]), _unescape(f[4]), f[sep + 1], _unescape(f[sep + 3])))
    return out


def _cgroup_dir(mount, path):
    """The directory of cgroup ``path`` through ``mount``, or ``None`` when
    the path is not visible there (a mount rooted elsewhere)."""
    root = mount.root
    if root == "/":
        rel = path
    elif path == root or path.startswith(root + "/"):
        rel = path[len(root):]
    else:
        return None
    return os.path.normpath(os.path.join(mount.mount_point, rel.lstrip("/")))


class _Malformed(Exception):
    pass


def _read_int_file(path, *, unlimited=None):
    """The integer of a one-value cgroup file; ``"max"`` (``unlimited``)
    gives ``float("inf")``; a missing file ``None``; anything else that does
    not parse, or a file that cannot be read, raises :class:`_Malformed`
    with the reason."""
    try:
        text = _read_text(path).strip()
    except FileNotFoundError:
        return None
    except (OSError, ValueError) as exc:
        raise _Malformed(_reason(exc))
    if unlimited is not None and text == unlimited:
        return float("inf")
    if not _DIGITS.fullmatch(text):
        raise _Malformed("not an integer: {!r}".format(text))
    return int(text)


def _read_level_file(path, *, unlimited=None):
    """:func:`_read_int_file`, with a failure reported as a WARNING naming
    the file and read as "no value": that level alone is ignored."""
    try:
        return _read_int_file(path, unlimited=unlimited)
    except _Malformed as exc:
        logger.warning("cgroup memory file %s could not be read (%s); that cgroup level is "
                       "ignored", path, exc)
        return None


def _inactive_file(directory, key):
    """The ``key`` value of ``<directory>/memory.stat`` in bytes; 0 when the
    file is missing, unreadable, malformed or has no such key (the
    conservative reading: the level still counts, with its whole usage). A
    file that exists but gives no value is named in one WARNING."""
    path = os.path.join(directory, "memory.stat")
    try:
        text = _read_text(path)
    except FileNotFoundError:
        return 0
    except (OSError, ValueError) as exc:
        problem = "could not be read ({})".format(_reason(exc))
    else:
        problem = "has no {} line".format(key)
        for line in text.splitlines():
            parts = line.split()
            if parts and parts[0] == key:
                if len(parts) == 2 and _DIGITS.fullmatch(parts[1]):
                    return int(parts[1])
                problem = "has a malformed {} line: {!r}".format(key, line)
                break
    logger.warning("cgroup memory file %s %s; its %s is read as 0 (that cgroup level still "
                   "counts, with its whole usage)", path, problem, key)
    return 0


_V2 = ("memory.max", "memory.current", "inactive_file")
_V1 = ("memory.limit_in_bytes", "memory.usage_in_bytes", "total_inactive_file")


def _mount_remaining(mount, path, v1):
    """The minimum of ``limit - working_set`` over the process's directory
    and its ancestors up to the mount point, or ``None`` (the path is not
    visible through the mount or no level has a limit). A level whose
    limit or usage file is missing, unreadable or malformed is skipped; the
    other levels still count."""
    directory = _cgroup_dir(mount, path)
    if directory is None:
        return None
    lim_name, use_name, inactive_key = _V1 if v1 else _V2
    top = os.path.normpath(mount.mount_point)
    best = None
    d = directory
    while True:
        limit = _read_level_file(os.path.join(d, lim_name), unlimited=None if v1 else "max")
        usage = _read_level_file(os.path.join(d, use_name))
        if limit is not None and usage is not None:
            if v1 and limit >= V1_UNLIMITED_FLOOR:
                limit = float("inf")
            if limit != float("inf"):
                working = max(0, usage - _inactive_file(d, inactive_key))
                remaining = int(limit) - working
                best = remaining if best is None else min(best, remaining)
        if d == top or os.path.dirname(d) == d:
            break
        d = os.path.dirname(d)
    return best


def _hierarchy_remaining(mounts, path, v1):
    values = [r for r in (_mount_remaining(m, path, v1) for m in mounts) if r is not None]
    return min(values) if values else None


def _cgroup_remaining():
    """``(nbytes, source)`` of the process's cgroup, or ``None``.

    The process's paths come from ``/proc/self/cgroup``, the mounts from
    ``/proc/self/mountinfo``: every ``cgroup2`` mount (v2) and every
    ``cgroup`` mount with ``memory`` in its super options (v1) is evaluated,
    and the minimum over the mounts that show the path is taken. At each
    level, from the process's directory up to the mount point, the
    remaining amount is ``limit - max(0, usage - inactive)`` (v2
    ``memory.max``, ``memory.current``, ``inactive_file``; v1
    ``memory.limit_in_bytes``, ``memory.usage_in_bytes``,
    ``total_inactive_file``); ``max`` or a v1 limit ``>= 2**62`` is no
    limit at that level. A v2 value is used whenever there is one (source
    ``"cgroup v2"``), else the v1 value (``"cgroup v1"``)."""
    paths = _read_proc_self_cgroup()
    if not paths:
        return None
    mounts = _read_mountinfo()
    if "v2" in paths:
        v2 = _hierarchy_remaining([m for m in mounts if m.fstype == "cgroup2"],
                                  paths["v2"], False)
        if v2 is not None:
            return v2, "cgroup v2"
    if "v1 memory" in paths:
        v1 = _hierarchy_remaining(
            [m for m in mounts
             if m.fstype == "cgroup" and "memory" in m.super_options.split(",")],
            paths["v1 memory"], True)
        if v1 is not None:
            return v1, "cgroup v1"
    return None


def available_host_bytes():
    """The available host memory (:class:`HostMemory`): the minimum over the
    remaining amounts of the sources in the module docstring, each clamped
    at 0; ``HostMemory(None, "none")`` when no source gives a value."""
    found = []
    ram = _read_psutil()
    if ram is not None:
        found.append((ram, "psutil"))
    else:
        ram = _read_meminfo()
        if ram is not None:
            found.append((ram, "meminfo"))
    limit = _read_rlimit("RLIMIT_AS")
    if limit is not None:
        used = _read_status("VmSize")
        if used is None:
            used = _read_psutil_vms()
        if used is not None:
            found.append((limit - used, "RLIMIT_AS"))
    limit = _read_rlimit("RLIMIT_DATA")
    if limit is not None:
        used = _read_status("VmData")
        if used is not None:
            found.append((limit - used, "RLIMIT_DATA"))
    cg = _cgroup_remaining()
    if cg is not None:
        found.append(cg)
    best = None
    for nbytes, source in found:
        nbytes = max(0, int(nbytes))
        if best is None or nbytes < best[0]:
            best = (nbytes, source)
    if best is None:
        return HostMemory(None, "none")
    return HostMemory(*best)

