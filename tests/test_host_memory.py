"""The available host memory and the memory cap of the spin-orbital bond
admission.

Every ``/proc`` and cgroup file is a file under a temporary directory: the
module constants of :mod:`hwave.solver.host_memory` are pointed at it and the
psutil and rlimit readers are patched, so no test reads the real host.
"""

import errno
import logging
import os
import sys
import tempfile
import types
import unittest
from unittest import mock

import numpy as np

from hwave.solver import host_memory, so_bond
from hwave.solver.host_memory import HostMemory, available_host_bytes
from hwave.solver.so_bond import SOBondCaps
from hwave.solver.so_records import compile_interaction
from tests.so_bond_oracles import CELL

GiB = 2 ** 30
#: the logger of the host memory module
_HM_LOG = "qlms.solver.host_memory"
MiB = 2 ** 20
#: the kernel's "unlimited" v1 memory limit (PAGE_COUNTER_MAX pages of 4 KiB)
V1_UNLIMITED = 9223372036854771712

try:
    import resource
except ImportError:          # pragma: no cover - not a POSIX host
    resource = None


def _escape(path):
    """A path as ``/proc/self/mountinfo`` writes it (octal escapes)."""
    return (path.replace("\\", "\\134").replace(" ", "\\040").replace("\t", "\\011")
            .replace("\n", "\\012"))


class _Host:
    """A fake host under a temporary directory.

    Nothing exists until a test writes it; psutil and every rlimit read
    ``None`` until set."""

    def __init__(self, case):
        tmp = tempfile.TemporaryDirectory()
        case.addCleanup(tmp.cleanup)
        self.root = tmp.name
        self.psutil = None
        self.vms = None
        self.rlimits = {}
        self.mounts = []
        for name, rel in (("PROC_SELF_CGROUP", "proc/self/cgroup"),
                          ("PROC_SELF_MOUNTINFO", "proc/self/mountinfo"),
                          ("PROC_SELF_STATUS", "proc/self/status"),
                          ("PROC_MEMINFO", "proc/meminfo")):
            p = mock.patch.object(host_memory, name, os.path.join(self.root, rel))
            p.start()
            case.addCleanup(p.stop)
        for name, fn in (("_read_psutil", lambda: self.psutil),
                         ("_read_psutil_vms", lambda: self.vms),
                         ("_read_rlimit", lambda which: self.rlimits.get(which))):
            p = mock.patch.object(host_memory, name, side_effect=fn)
            p.start()
            case.addCleanup(p.stop)

    def path(self, rel):
        return os.path.join(self.root, rel)

    def write(self, rel, text):
        p = self.path(rel)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        with open(p, "w") as f:
            f.write(text)
        return p

    def cgroup(self, *lines):
        self.write("proc/self/cgroup", "".join(line + "\n" for line in lines))

    def status(self, *lines):
        self.write("proc/self/status", "".join(line + "\n" for line in lines))

    def meminfo(self, *lines):
        self.write("proc/meminfo", "".join(line + "\n" for line in lines))

    def mount(self, rel, *, root="/", fstype="cgroup2", options=None):
        """Add a mountinfo line for the directory ``<tmp>/<rel>`` (created)."""
        mp = self.path(rel)
        os.makedirs(mp, exist_ok=True)
        if options is None:
            options = "rw,nsdelegate" if fstype == "cgroup2" else "rw,memory"
        n = len(self.mounts) + 30
        self.mounts.append("{} 25 0:{} {} {} rw,nosuid,nodev shared:{} - {} {} {}".format(
            n, n, _escape(root), _escape(mp), n, fstype, fstype, options))
        self.write("proc/self/mountinfo", "".join(m + "\n" for m in self.mounts))
        return mp

    def level(self, mount_dir, sub, *, limit=None, usage=None, stat=None, v1=False):
        """The memory files of one cgroup directory ``<mount_dir>/<sub>``."""
        d = os.path.join(mount_dir, sub.strip("/")) if sub.strip("/") else mount_dir
        os.makedirs(d, exist_ok=True)
        lim, cur = (("memory.limit_in_bytes", "memory.usage_in_bytes") if v1
                    else ("memory.max", "memory.current"))
        for name, value in ((lim, limit), (cur, usage), ("memory.stat", stat)):
            if value is not None:
                with open(os.path.join(d, name), "w") as f:
                    f.write("{}\n".format(value))
        return d


class TestAvailable(unittest.TestCase):
    """The minimum over the remaining amounts of every source."""

    def setUp(self):
        self.host = _Host(self)

    def test_minimum_of_remaining(self):
        h = self.host
        h.psutil = 10 * GiB
        h.rlimits["RLIMIT_AS"] = 8 * GiB
        h.status("Name:\tpython", "VmSize:\t 3145728 kB", "VmData:\t 1048576 kB")
        h.cgroup("0::/user.slice/job1")
        mp = h.mount("cg")
        h.level(mp, "/user.slice/job1", limit=6 * GiB, usage=2 * GiB)
        self.assertEqual(available_host_bytes(), HostMemory(4 * GiB, "cgroup v2"))
        # RLIMIT_AS 8 - VmSize 3 = 5 is below psutil once the cgroup is gone
        os.remove(h.path("proc/self/cgroup"))
        self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "RLIMIT_AS"))
        # a tie goes to the first source in the documented order
        h.psutil = 5 * GiB
        self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "psutil"))

    def test_meminfo_source(self):
        h = self.host
        h.meminfo("MemTotal:       16777216 kB", "MemFree:         1048576 kB",
                  "MemAvailable:    6291456 kB")
        self.assertEqual(available_host_bytes(), HostMemory(6 * GiB, "meminfo"))
        h.psutil = 5 * GiB
        with mock.patch.object(host_memory, "_read_meminfo",
                               wraps=host_memory._read_meminfo) as spy:
            self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "psutil"))
        spy.assert_not_called()

    def test_units(self):
        h = self.host
        h.rlimits["RLIMIT_AS"] = 8 * GiB
        h.status("VmSize:  3145728 kB")
        self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "RLIMIT_AS"))
        self.assertEqual(host_memory._parse_kib_line("VmSize:  3145728 kB", "VmSize"),
                         3 * GiB)
        for bad in ("VmSize: 3 MB", "VmSize: 3145728", "VmSize: -1 kB", "VmSize: x kB",
                    "VmSize: 3.5 kB", "VmSize: 1 kB extra"):
            with self.subTest(line=bad):
                with self.assertRaises(ValueError):
                    host_memory._parse_kib_line(bad, "VmSize")
                h.status(bad)
                self.assertEqual(available_host_bytes(), HostMemory(None, "none"))
        with self.assertRaises(ValueError):
            host_memory._parse_kib_line("VmData: 1 kB", "VmSize")

    @unittest.skipIf(resource is None, "no resource module on this host")
    def test_rlimit_infinity(self):
        # the module's own reader, under a patched resource.getrlimit
        real = _unpatched("_read_rlimit")
        with mock.patch("resource.getrlimit",
                        return_value=(resource.RLIM_INFINITY, resource.RLIM_INFINITY)):
            self.assertIsNone(real("RLIMIT_AS"))
        with mock.patch("resource.getrlimit", return_value=(-1, -1)):
            self.assertIsNone(real("RLIMIT_AS"))
        with mock.patch("resource.getrlimit", return_value=(8 * GiB, resource.RLIM_INFINITY)):
            self.assertEqual(real("RLIMIT_AS"), 8 * GiB)
        with mock.patch("resource.getrlimit", side_effect=OSError("denied")):
            self.assertIsNone(real("RLIMIT_AS"))
        self.assertIsNone(real("RLIMIT_NO_SUCH_LIMIT"))
        # through the public entry: infinity skipped, a finite limit counted
        h = self.host
        h.status("VmSize:\t 1048576 kB", "VmData:\t 1048576 kB")
        for soft, want in ((resource.RLIM_INFINITY, HostMemory(None, "none")),
                           (4 * GiB, HostMemory(3 * GiB, "RLIMIT_AS"))):
            with self.subTest(soft=soft):
                with mock.patch.object(host_memory, "_read_rlimit", side_effect=lambda w: (
                        real(w) if w == "RLIMIT_AS" else None)), \
                        mock.patch("resource.getrlimit", return_value=(soft, soft)):
                    self.assertEqual(available_host_bytes(), want)

    def test_reader_errors_are_no_value(self):
        h = self.host
        h.meminfo("MemAvailable:    6291456 kB")
        h.status("VmSize:\t 1048576 kB", "VmData:\t 1048576 kB")
        h.rlimits.update(RLIMIT_AS=8 * GiB, RLIMIT_DATA=8 * GiB)
        h.cgroup("0::/job")
        mp = h.mount("cg")
        h.level(mp, "/job", limit=6 * GiB, usage=GiB)
        self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "cgroup v2"))
        with mock.patch.object(host_memory, "open", create=True,
                               side_effect=OSError("unreadable")), \
                self.assertLogs(_HM_LOG, level="WARNING"):
            for name in ("_read_meminfo", "_cgroup_remaining", "_read_proc_self_cgroup",
                         "_read_mountinfo"):
                with self.subTest(reader=name):
                    got = getattr(host_memory, name)()
                    self.assertIn(got, (None, {}, []))
            self.assertIsNone(host_memory._read_status("VmSize"))
            self.assertEqual(available_host_bytes(), HostMemory(None, "none"))
        # a non-integer field: that source alone is skipped
        h.meminfo("MemAvailable:    x kB")
        h.status("VmSize:\t 1.5 kB", "VmData:\t 1048576 kB")
        self.assertIsNone(host_memory._read_meminfo())
        self.assertIsNone(host_memory._read_status("VmSize"))
        self.assertEqual(available_host_bytes(), HostMemory(5 * GiB, "cgroup v2"))
        # psutil: missing (ImportError) or failing (OSError, ValueError)
        real = _unpatched("_read_psutil")
        real_vms = _unpatched("_read_psutil_vms")
        with mock.patch.dict(sys.modules, {"psutil": None}):
            self.assertIsNone(real())
            self.assertIsNone(real_vms())
        for exc in (OSError("no /proc"), ValueError("parse")):
            fake = mock.Mock()
            fake.virtual_memory.side_effect = exc
            fake.Process.side_effect = exc
            with self.subTest(exc=type(exc).__name__), \
                    mock.patch.dict(sys.modules, {"psutil": fake}):
                self.assertIsNone(real())
                self.assertIsNone(real_vms())
        fake = mock.Mock()
        fake.virtual_memory.return_value = mock.Mock(available=3 * GiB)
        fake.Process.return_value.memory_info.return_value = mock.Mock(vms=2 * GiB)
        with mock.patch.dict(sys.modules, {"psutil": fake}):
            self.assertEqual(real(), 3 * GiB)
            self.assertEqual(real_vms(), 2 * GiB)

    def test_psutil_error_is_no_value(self):
        """A psutil error (``psutil.Error``, e.g. ``AccessDenied``) from
        ``virtual_memory()`` or ``Process().memory_info()`` makes that reader
        "no value"; the other sources still give the minimum."""
        h = self.host
        h.meminfo("MemAvailable:    6291456 kB")      # 6 GiB
        h.status("VmData:\t 1048576 kB")             # no VmSize: RLIMIT_AS asks psutil
        h.rlimits.update(RLIMIT_AS=4 * GiB, RLIMIT_DATA=8 * GiB)
        h.cgroup("0::/job")
        mp = h.mount("cg")
        h.level(mp, "/job", limit=6 * GiB, usage=GiB)

        class Error(Exception):
            pass

        class AccessDenied(Error):
            pass

        fake = types.ModuleType("psutil")
        fake.Error, fake.AccessDenied = Error, AccessDenied
        modules = [("module", fake, AccessDenied)]
        try:
            import psutil
        except ImportError:                  # pragma: no cover - psutil absent
            pass
        else:
            modules.append(("psutil", psutil, psutil.AccessDenied))
        for label, module, denied in modules:
            for where in ("virtual_memory", "memory_info", "both"):
                with self.subTest(module=label, where=where), \
                        mock.patch.object(host_memory, "_read_psutil",
                                          _unpatched("_read_psutil")), \
                        mock.patch.object(host_memory, "_read_psutil_vms",
                                          _unpatched("_read_psutil_vms")), \
                        mock.patch.dict(sys.modules, {"psutil": module}):
                    vm = mock.Mock(return_value=mock.Mock(available=3 * GiB))
                    info = mock.Mock(return_value=mock.Mock(vms=3 * GiB))
                    if where in ("virtual_memory", "both"):
                        vm = mock.Mock(side_effect=denied("virtual_memory"))
                    if where in ("memory_info", "both"):
                        info = mock.Mock(side_effect=denied("memory_info"))
                    with mock.patch.object(module, "virtual_memory", vm, create=True), \
                            mock.patch.object(module, "Process",
                                              mock.Mock(return_value=mock.Mock(
                                                  memory_info=info)), create=True):
                        ram = host_memory._read_psutil()
                        vms = host_memory._read_psutil_vms()
                        got = available_host_bytes()
                    self.assertEqual(ram, None if where != "memory_info" else 3 * GiB)
                    self.assertEqual(vms, None if where != "virtual_memory" else 3 * GiB)
                    # psutil RAM 3 GiB when it reads, else meminfo 6 GiB; RLIMIT_AS
                    # 4 - 3 = 1 GiB only with the vms value; RLIMIT_DATA 7 GiB;
                    # cgroup 6 - 1 = 5 GiB
                    want = {"virtual_memory": HostMemory(GiB, "RLIMIT_AS"),
                            "memory_info": HostMemory(3 * GiB, "psutil"),
                            "both": HostMemory(5 * GiB, "cgroup v2")}[where]
                    self.assertEqual(got, want)

    def test_clamped_and_data_without_vmdata_ignored(self):
        h = self.host
        h.psutil = 10 * GiB
        h.rlimits.update(RLIMIT_AS=2 * GiB, RLIMIT_DATA=GiB)
        h.status("VmSize:\t 3145728 kB")              # above RLIMIT_AS; no VmData
        self.assertEqual(available_host_bytes(), HostMemory(0, "RLIMIT_AS"))
        # without VmSize, the psutil vms value is the usage of RLIMIT_AS
        os.remove(h.path("proc/self/status"))
        h.vms = GiB + 512 * MiB
        self.assertEqual(available_host_bytes(), HostMemory(512 * MiB, "RLIMIT_AS"))
        # no usage value at all: both rlimits skipped
        h.vms = None
        self.assertEqual(available_host_bytes(), HostMemory(10 * GiB, "psutil"))
        # RLIMIT_DATA counts once VmData is there
        h.status("VmData:\t 262144 kB")
        self.assertEqual(available_host_bytes(), HostMemory(768 * MiB, "RLIMIT_DATA"))

    def test_no_measurement(self):
        self.assertEqual(available_host_bytes(), HostMemory(None, "none"))
        self.assertEqual(HostMemory._fields, ("nbytes", "source"))


#: the module's own readers, taken before any test patches them
_ORIGINAL = {name: getattr(host_memory, name)
             for name in ("_read_psutil", "_read_psutil_vms", "_read_rlimit")}


def _unpatched(name):
    """The module's own reader ``name``, bypassing the test patch."""
    return _ORIGINAL[name]


class TestCgroup(unittest.TestCase):
    """The cgroup v2 / v1 remaining memory of this process."""

    def setUp(self):
        self.host = _Host(self)

    def remaining(self):
        return host_memory._cgroup_remaining()

    def test_v2_unified(self):
        h = self.host
        h.cgroup("0::/user.slice/job1")
        mp = h.mount("cg")
        h.level(mp, "/user.slice/job1", limit=8589934592, usage=2147483648)
        self.assertEqual(self.remaining(), (6 * GiB, "cgroup v2"))
        self.assertEqual(available_host_bytes(), HostMemory(6 * GiB, "cgroup v2"))

    def test_v2_nested(self):
        h = self.host
        h.cgroup("0::/a/b")
        mp = h.mount("cg")
        h.level(mp, "/a/b", limit="max", usage=GiB)
        h.level(mp, "/a", limit=4 * GiB, usage=GiB)
        h.level(mp, "/", limit=16 * GiB, usage=15 * GiB + 512 * MiB)
        self.assertEqual(self.remaining(), (512 * MiB, "cgroup v2"))
        # every level unlimited: no value
        for sub in ("/a", "/"):
            h.level(mp, sub, limit="max")
        self.assertIsNone(self.remaining())

    def test_v1_memory(self):
        h = self.host
        h.cgroup("12:cpu,cpuacct:/slurm/uid_1/job_2", "5:memory:/slurm/uid_1/job_2",
                 "1:name=systemd:/user.slice", "0::/user.slice")
        h.mount("cpu", fstype="cgroup", options="rw,cpu,cpuacct")
        mp = h.mount("mem", fstype="cgroup", options="rw,memory")
        d = h.level(mp, "/slurm/uid_1/job_2", limit=8 * GiB, usage=3 * GiB, v1=True)
        self.assertEqual(self.remaining(), (5 * GiB, "cgroup v1"))
        with open(os.path.join(d, "memory.limit_in_bytes"), "w") as f:
            f.write("{}\n".format(V1_UNLIMITED))
        self.assertIsNone(self.remaining())
        self.assertEqual(host_memory._read_proc_self_cgroup(),
                         {"v2": "/user.slice", "v1 memory": "/slurm/uid_1/job_2"})

    def test_hybrid(self):
        h = self.host
        h.cgroup("4:memory:/job", "0::/job")
        uni = h.mount("unified")
        h.level(uni, "/job", limit="max", usage=GiB)
        h.level(uni, "/", limit="max", usage=GiB)
        mem = h.mount("memory", fstype="cgroup", options="rw,memory")
        h.level(mem, "/job", limit=6 * GiB, usage=GiB, v1=True)
        self.assertEqual(self.remaining(), (5 * GiB, "cgroup v1"))
        # a v2 limit wins even when it leaves more than the v1 one
        h.level(uni, "/job", limit=8 * GiB, usage=GiB)
        self.assertEqual(self.remaining(), (7 * GiB, "cgroup v2"))

    def test_container_root(self):
        h = self.host
        mp = h.mount("cg", root="/docker/abc")
        h.level(mp, "/", limit=4 * GiB, usage=GiB)
        h.level(mp, "/sub", limit=2 * GiB, usage=0)
        for path, want in (("/docker/abc", (3 * GiB, "cgroup v2")),
                           ("/docker/abc/sub", (2 * GiB, "cgroup v2")),
                           ("/docker/abcd", None), ("/other", None)):
            with self.subTest(path=path):
                h.cgroup("0::" + path)
                self.assertEqual(self.remaining(), want)

    def test_escaped_mount_point(self):
        h = self.host
        h.cgroup("0::/job")
        mp = h.mount("my cg")
        self.assertIn("\\040", h.mounts[0])
        h.level(mp, "/job", limit=3 * GiB, usage=GiB)
        self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
        self.assertEqual(host_memory._read_mountinfo()[0].mount_point, mp)

    def test_v2_first_mount_incompatible(self):
        h = self.host
        h.cgroup("0::/user.slice/job1")
        h.mount("other", root="/docker/other")
        mp = h.mount("cg")
        h.level(mp, "/user.slice/job1", limit=8 * GiB, usage=2 * GiB)
        self.assertEqual(self.remaining(), (6 * GiB, "cgroup v2"))

    def test_v1_first_mount_incompatible(self):
        h = self.host
        h.cgroup("3:memory:/job")
        h.mount("cpu", fstype="cgroup", options="rw,cpu")        # no memory: skipped
        h.mount("other", root="/other", fstype="cgroup", options="rw,memory")
        mp = h.mount("mem", fstype="cgroup", options="rw,memory")
        h.level(mp, "/job", limit=8 * GiB, usage=3 * GiB, v1=True)
        self.assertEqual(self.remaining(), (5 * GiB, "cgroup v1"))
        mounts = host_memory._read_mountinfo()
        self.assertEqual([m.fstype for m in mounts], ["cgroup"] * 3)
        self.assertEqual([m.root for m in mounts], ["/", "/other", "/"])

    def test_multiple_valid_roots(self):
        for order in ((("a", "/"), ("b", "/user.slice")), (("b", "/user.slice"), ("a", "/"))):
            with self.subTest(order=order):
                h = _Host(self)
                h.cgroup("0::/user.slice/job1")
                dirs = {name: h.mount(name, root=root) for name, root in order}
                h.level(dirs["a"], "/user.slice/job1", limit=8 * GiB, usage=2 * GiB)
                h.level(dirs["b"], "/job1", limit=4 * GiB, usage=GiB)
                self.assertEqual(host_memory._cgroup_remaining(), (3 * GiB, "cgroup v2"))

    def test_inactive_file_subtracted(self):
        h = self.host
        h.cgroup("0::/job")
        mp = h.mount("cg")
        h.level(mp, "/job", limit=8 * GiB, usage=6 * GiB,
                stat="anon 1024\nfile 4096\ninactive_file 2147483648\nactive_file 0")
        self.assertEqual(self.remaining(), (4 * GiB, "cgroup v2"))
        # the cache is larger than the usage: the working set clamps at 0
        h.level(mp, "/job", stat="inactive_file {}".format(7 * GiB))
        self.assertEqual(self.remaining(), (8 * GiB, "cgroup v2"))
        for stat in ("inactive_file x", "anon 1024\nactive_file 2147483648", "garbage",
                     "inactive_file", "inactive_file -1"):
            with self.subTest(stat=stat):
                h.level(mp, "/job", stat=stat)
                with self.assertLogs(_HM_LOG, level="WARNING"):
                    self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
        os.remove(os.path.join(mp, "job", "memory.stat"))
        self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
        # v1 subtracts total_inactive_file (not the local inactive_file)
        h2 = _Host(self)
        h2.cgroup("7:memory:/job")
        m1 = h2.mount("mem", fstype="cgroup", options="rw,memory")
        h2.level(m1, "/job", limit=8 * GiB, usage=6 * GiB, v1=True,
                 stat="inactive_file 1073741824\ntotal_inactive_file 2147483648")
        self.assertEqual(host_memory._cgroup_remaining(), (4 * GiB, "cgroup v1"))

    def test_memory_stat_unusable_warns(self):
        """An unreadable or malformed memory.stat, or one without the
        inactive-file key, gives one WARNING naming the file; its inactive
        file cache reads as 0 and the level still counts. A missing
        memory.stat is silent."""
        h = self.host
        h.cgroup("0::/job")
        mp = h.mount("cg")
        h.level(mp, "/job", limit=8 * GiB, usage=6 * GiB)
        stat = os.path.join(mp, "job", "memory.stat")
        with mock.patch.object(host_memory.logger, "warning") as warned:
            self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
        warned.assert_not_called()
        for text in ("inactive_file x", "inactive_file", "inactive_file -1",
                     "inactive_file 1 2", "garbage", "anon 1024\nactive_file 2147483648"):
            with self.subTest(stat=text):
                h.level(mp, "/job", stat=text)
                with self.assertLogs(_HM_LOG, level="WARNING") as cm:
                    self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
                self.assertEqual(len(cm.records), 1)
                self.assertIn(stat, cm.output[0])
                self.assertIn("inactive_file", cm.output[0])
        real_open = open

        def unreadable(path, *a, **kw):
            if path == stat:
                raise PermissionError(errno.EACCES, "Permission denied")
            return real_open(path, *a, **kw)
        h.level(mp, "/job", stat="inactive_file 1073741824")
        self.assertEqual(self.remaining(), (3 * GiB, "cgroup v2"))
        with mock.patch.object(host_memory, "open", create=True, side_effect=unreadable), \
                self.assertLogs(_HM_LOG, level="WARNING") as cm:
            self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))
        self.assertEqual(len(cm.records), 1)
        self.assertIn(stat, cm.output[0])
        self.assertIn("Permission denied", cm.output[0])
        # v1: the key is total_inactive_file
        h2 = _Host(self)
        h2.cgroup("7:memory:/job")
        m1 = h2.mount("mem", fstype="cgroup", options="rw,memory")
        h2.level(m1, "/job", limit=8 * GiB, usage=6 * GiB, v1=True,
                 stat="inactive_file 1073741824\ntotal_inactive_file x")
        with self.assertLogs(_HM_LOG, level="WARNING") as cm:
            self.assertEqual(host_memory._cgroup_remaining(), (2 * GiB, "cgroup v1"))
        self.assertEqual(len(cm.records), 1)
        self.assertIn("total_inactive_file", cm.output[0])

    def test_malformed(self):
        """A malformed /proc file gives no cgroup value (with a WARNING); a
        bad file at one level skips only that level (with a WARNING), so the
        limit of the mount point level still counts."""
        h = self.host
        h.psutil = 9 * GiB
        mp = h.mount("cg")
        good_mounts = "".join(m + "\n" for m in h.mounts)
        job = os.path.join(mp, "job")
        cg = h.path("proc/self/cgroup")
        mi = h.path("proc/self/mountinfo")

        def reset():
            h.cgroup("0::/job")
            h.write("proc/self/mountinfo", good_mounts)
            h.level(mp, "/", limit=5 * GiB, usage=GiB)
            h.level(mp, "/job", limit=3 * GiB, usage=GiB)
            self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))

        def line_text(path, n, line):
            return ("{} line {} is malformed: {!r}; the cgroup memory limit is not read"
                    .format(path, n, line))

        def level_text(name, reason):
            return ("cgroup memory file {} could not be read ({}); that cgroup level is "
                    "ignored".format(os.path.join(job, name), reason))

        short = "30 25 0:30 / {} - cgroup2".format(mp)
        no_cgroup = (None, HostMemory(9 * GiB, "psutil"))
        parent = ((4 * GiB, "cgroup v2"), HostMemory(4 * GiB, "cgroup v2"))
        for name, breakage, want, warning in (
                ("cgroup line", lambda: h.cgroup("0::/job", "garbage"), no_cgroup,
                 line_text(cg, 2, "garbage")),
                ("cgroup path", lambda: h.cgroup("0::job"), no_cgroup,
                 line_text(cg, 1, "0::job")),
                ("cgroup colons", lambda: h.cgroup("0::/job:x"), no_cgroup,
                 line_text(cg, 1, "0::/job:x")),
                ("cgroup missing", lambda: os.remove(cg), no_cgroup, None),
                ("mountinfo separator", lambda: h.write(
                    "proc/self/mountinfo", good_mounts.replace(" - ", " ")), no_cgroup,
                 line_text(mi, 1, good_mounts.replace(" - ", " ").rstrip("\n"))),
                ("mountinfo short", lambda: h.write("proc/self/mountinfo", short + "\n"),
                 no_cgroup, line_text(mi, 1, short)),
                ("mountinfo missing", lambda: os.remove(mi), no_cgroup, None),
                ("memory.max", lambda: h.level(mp, "/job", limit="12x"), parent,
                 level_text("memory.max", "not an integer: '12x'")),
                ("memory.current", lambda: h.level(mp, "/job", usage=""), parent,
                 level_text("memory.current", "not an integer: ''")),
                ("memory.current missing", lambda: os.remove(os.path.join(job, "memory.current")),
                 parent, None),
                ("memory.max missing", lambda: os.remove(os.path.join(job, "memory.max")),
                 parent, None)):
            with self.subTest(case=name):
                reset()
                breakage()
                if warning is None:
                    with self.assertNoLogs(_HM_LOG, level="WARNING"):
                        self.assertEqual(self.remaining(), want[0])
                else:
                    with self.assertLogs(_HM_LOG, level="WARNING") as cm:
                        self.assertEqual(self.remaining(), want[0])
                    self.assertEqual(cm.output, ["WARNING:{}:{}".format(_HM_LOG, warning)])
                with self.assertLogs(_HM_LOG, level="WARNING") if warning else \
                        self.assertNoLogs(_HM_LOG, level="WARNING"):
                    self.assertEqual(available_host_bytes(), want[1])
        reset()
        h.write("proc/self/mountinfo", good_mounts.replace(" - ", " "))
        with self.assertLogs(_HM_LOG, level="WARNING"):
            self.assertEqual(host_memory._read_mountinfo(), [])
        h.cgroup("0::/job", "garbage")
        with self.assertLogs(_HM_LOG, level="WARNING"):
            self.assertEqual(host_memory._read_proc_self_cgroup(), {})

    def test_bad_level_keeps_other_levels(self):
        """A malformed or unreadable leaf file does not hide a tighter parent
        limit."""
        h = self.host
        h.psutil = 64 * GiB
        h.cgroup("0::/job.slice/step")
        mp = h.mount("cg")
        h.level(mp, "/job.slice", limit=GiB, usage=0)
        leaf = h.level(mp, "/job.slice/step", limit="12x", usage=0)
        with self.assertLogs(_HM_LOG, level="WARNING") as cm:
            self.assertEqual(available_host_bytes(), HostMemory(GiB, "cgroup v2"))
        self.assertEqual(cm.output, [
            "WARNING:qlms.solver.host_memory:cgroup memory file {} could not be read "
            "(not an integer: '12x'); that cgroup level is ignored".format(
                os.path.join(leaf, "memory.max"))])
        # a level file that cannot be opened (a directory: EISDIR)
        os.remove(os.path.join(leaf, "memory.max"))
        os.mkdir(os.path.join(leaf, "memory.max"))
        with self.assertLogs(_HM_LOG, level="WARNING") as cm:
            self.assertEqual(available_host_bytes(), HostMemory(GiB, "cgroup v2"))
        self.assertEqual(cm.output, [
            "WARNING:qlms.solver.host_memory:cgroup memory file {} could not be read "
            "({}); that cgroup level is ignored".format(
                os.path.join(leaf, "memory.max"), os.strerror(errno.EISDIR))])

    def test_mountinfo_undecodable_line_skipped(self):
        h = self.host
        h.cgroup("0::/job")
        mp = h.mount("cg")
        h.level(mp, "/job", limit=3 * GiB, usage=GiB)
        good = "".join(m + "\n" for m in h.mounts).encode()
        bad = b"40 25 8:1 / /mnt/caf\xe9 rw,relatime shared:2 - ext4 /dev/sda1 rw\n"
        with open(h.path("proc/self/mountinfo"), "wb") as f:
            f.write(bad + good)
        with self.assertLogs(_HM_LOG, level="WARNING") as cm:
            mounts = host_memory._read_mountinfo()
        self.assertEqual(cm.output, [
            "WARNING:qlms.solver.host_memory:{} line 1 is not valid UTF-8; that mount is "
            "ignored".format(h.path("proc/self/mountinfo"))])
        self.assertEqual([m.mount_point for m in mounts], [mp])
        with self.assertLogs(_HM_LOG, level="WARNING"):
            self.assertEqual(self.remaining(), (2 * GiB, "cgroup v2"))


KEY = "[eliashberg] bond_memory_cap_gb"
SHELLS_KEY = "[eliashberg] bond_max_shells"
_LOG = logging.getLogger("hwave.test.host_memory")


def _chain_tables():
    """A 1-D chain (norb_phys = 1): CoulombIntra, nearest-neighbour
    CoulombInter and a weaker Hund (records of different magnitudes on
    shell 1), and the second shell declared with zero coupling only."""
    inter = {((1, 0, 0), (0, 0)): 0.9, ((-1, 0, 0), (0, 0)): 0.9,
             ((2, 0, 0), (0, 0)): 0.0, ((-2, 0, 0), (0, 0)): 0.0}
    hund = {((1, 0, 0), (0, 0)): 0.1, ((-1, 0, 0), (0, 0)): 0.1}
    return {"CoulombIntra": {((0, 0, 0), (0, 0)): 0.83}, "CoulombInter": inter, "Hund": hund}


class TestCaps(unittest.TestCase):
    """The spin-orbital memory cap values and the over-cap refusal."""

    def resolve(self, raw, avail=None):
        info = {}
        with mock.patch.object(host_memory, "available_host_bytes",
                               return_value=avail or HostMemory(None, "none")):
            with self.assertLogs(_LOG, level="INFO") as cm:
                caps = so_bond.resolve_host_cap(raw, key=KEY, log=_LOG, info_out=info)
        self.assertEqual(len(cm.records), 1)
        rec = cm.records[0]
        return caps, (rec.levelname, rec.getMessage()), info["summary"]

    def test_semantics(self):
        caps, log, summary = self.resolve(None, HostMemory(10 * GiB, "cgroup v2"))
        self.assertEqual(caps, SOBondCaps(0.8 * 10 * GiB, None, None))
        self.assertEqual(log, ("INFO", "spin-orbital bond memory cap: 8.00 GiB = 0.8 x "
                                       "available 10.00 GiB (source: cgroup v2)"))
        self.assertEqual(summary, "8.000 GiB (0.8 x available 10.000 GiB, source: cgroup v2)")

        caps, log, summary = self.resolve(None, HostMemory(None, "none"))
        self.assertEqual(caps, SOBondCaps(None, None, None))
        self.assertEqual(log, ("WARNING", "no memory limit could be determined; the admission "
                                          "check is disabled; set [eliashberg] "
                                          "bond_memory_cap_gb to enable it"))
        self.assertEqual(summary, "unlimited (no memory limit could be determined)")

        for raw, want_raw in ((2.5, 2.5), (2, 2.0), (np.float64(2.5), 2.5)):
            with self.subTest(raw=raw):
                caps, log, summary = self.resolve(raw)
                self.assertEqual(caps, SOBondCaps(want_raw * GiB, want_raw, None))
                self.assertIs(type(caps.host_raw), float)
                self.assertEqual(log, ("INFO", "spin-orbital bond memory cap: {:.2f} GiB "
                                               "([eliashberg] bond_memory_cap_gb)".format(
                                                   want_raw)))
                self.assertEqual(summary, "{:.3f} GiB ([eliashberg] bond_memory_cap_gb = {!r})"
                                 .format(want_raw, want_raw))

        for raw in (0, 0.0):
            with self.subTest(raw=raw):
                caps, log, summary = self.resolve(raw)
                self.assertEqual(caps, SOBondCaps(None, 0.0, None))
                self.assertEqual(log, ("WARNING", "[eliashberg] bond_memory_cap_gb = 0: the "
                                                  "memory admission check of the spin-orbital "
                                                  "bond path is disabled"))
                self.assertEqual(summary, "unlimited ([eliashberg] bond_memory_cap_gb = 0)")

        for raw, shown in ((-1.0, "-1.0"), (float("nan"), "nan"), (float("inf"), "inf"),
                           (True, "True"), (np.bool_(False), repr(np.bool_(False))),
                           ("8", "'8'"), (10 ** 400, repr(10 ** 400))):
            with self.subTest(raw=raw):
                info = {}
                with self.assertRaises(ValueError) as cm:
                    so_bond.resolve_host_cap(raw, key=KEY, log=_LOG, info_out=info)
                self.assertEqual(str(cm.exception),
                                 "[eliashberg] bond_memory_cap_gb must be a finite number >= 0 "
                                 "in spin-orbital mode (0 disables the check), got " + shown)
                self.assertEqual(info, {})
        self.assertEqual(so_bond.CAP_RANGE_TEXT,
                         "{key} must be a finite number >= 0 in spin-orbital mode (0 disables "
                         "the check), got {raw!r}")
        # without info_out
        with mock.patch.object(host_memory, "available_host_bytes",
                               return_value=HostMemory(GiB, "meminfo")), \
                self.assertLogs(_LOG, level="INFO"):
            self.assertEqual(so_bond.resolve_host_cap(None, key=KEY, log=_LOG),
                             SOBondCaps(0.8 * GiB, None, None))

    def test_refusal_message_tables(self):
        compiled = compile_interaction(_chain_tables(), 1)
        topo = so_bond.so_bond_topology(compiled, CELL)
        rows = so_bond.shell_rows(compiled, topo, CELL)
        off_site = [abs(complex(r["v"])) for r in compiled.records
                    if np.any(np.asarray(r["cells"]))]
        self.assertEqual([(r.shell, r.length, sorted(r.displacements)) for r in rows],
                         [(1, 1.0, [(-1, 0, 0), (1, 0, 0)]), (2, 2.0, [(-2, 0, 0), (2, 0, 0)])])
        self.assertEqual(rows[0].max_abs_v, max(off_site))
        self.assertLess(min(off_site), max(off_site))
        self.assertGreater(rows[0].max_abs_v, 0.0)
        self.assertEqual(rows[1].max_abs_v, 0.0)

        for layout in ("dense", "sz_packed"):
            with self.subTest(layout=layout):
                text = so_bond.envelope_text(layout)
                marked = [ln for ln in text.splitlines() if ln.endswith("  <- this run")]
                self.assertEqual(len(marked), 1)
                self.assertTrue(marked[0].lstrip().startswith(
                    "dense layout" if layout == "dense" else "S_z-packed layout"), marked)
        with self.assertRaisesRegex(ValueError, "^envelope_text: layout must be 'dense' or "
                                                "'sz_packed', got 'packed'$"):
            so_bond.envelope_text("packed")

        envelope = so_bond.envelope_text("dense")
        kw = dict(nmat=2, nvol=5, B=topo.delta_r.shape[0], nd=2, consumers=("vertex",))
        plan_rows = so_bond.plan_so_bond(caps=SOBondCaps(None, None, None), **kw).rows
        with self.assertRaises(ValueError) as cm:
            so_bond.plan_so_bond(caps=SOBondCaps(1024.0, 1e-6, None), cap_key=KEY,
                                 shell_rows=rows, max_shells_key=SHELLS_KEY,
                                 envelope=envelope, **kw)
        msg = str(cm.exception)
        sentence = (
            "spin-orbital bond FLEX needs an estimated peak of 0.00 GiB in phase vertex (B = 5 "
            "channels, ND = B·(2·norb_phys)² = 20, nvol = 5, Nmat = 2; one store slot = "
            "16·Nmat·nvol·ND² bytes = 0.00 GiB) but the cap is 0.00 GiB ([eliashberg] "
            "bond_memory_cap_gb). Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce Nmat "
            "or the k grid, or raise [eliashberg] bond_memory_cap_gb (0 disables the check); "
            "pruning bond shells helps only for shells the per-shell table below shows with "
            "negligible coupling.")
        row_table = "rows (phase  name  bytes):\n" + "\n".join(
            "  {}  {}  {}".format(*r) for r in plan_rows)
        shells = ("shells (shell  |dr|  max|v|):\n"
                  "  1  1  {:.6g}\n"
                  "  2  2  0".format(rows[0].max_abs_v))
        clause = ("1 outer shell(s) declared only with zero coupling; [eliashberg] "
                  "bond_max_shells = 1 removes them")
        self.assertEqual(msg, "\n".join((sentence, row_table, envelope, shells, clause)))
        self.assertEqual(row_table, "rows (phase  name  bytes):\n  resident  gamma_b  32000")

        # no zero-coupling outer shell: no clause
        with self.assertRaises(ValueError) as cm:
            so_bond.plan_so_bond(caps=SOBondCaps(1024.0, 1e-6, None), cap_key=KEY,
                                 shell_rows=rows[:1], max_shells_key=SHELLS_KEY,
                                 envelope=envelope, **kw)
        self.assertTrue(str(cm.exception).endswith(
            "shells (shell  |dr|  max|v|):\n  1  1  {:.6g}".format(rows[0].max_abs_v)))

        # without cap_key: today's text, whatever the new keywords hold
        today = (
            "spin-orbital bond FLEX needs an estimated peak of 0.00 GiB in phase vertex (B = 5 "
            "channels, ND = B·(2·norb_phys)² = 80, nvol = 5, Nmat = 2; one store slot = "
            "16·Nmat·nvol·ND² bytes = 0.00 GiB) but the cap is 0.00 GiB (caps.host_bytes). "
            "Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce Nmat or the k grid, or prune "
            "bond shells.\nrows (phase  name  bytes):\n  resident  gamma_b  512000")
        for extra in ({}, dict(shell_rows=rows, max_shells_key=SHELLS_KEY, envelope=envelope)):
            with self.subTest(extra=sorted(extra)):
                with self.assertRaises(ValueError) as cm:
                    so_bond.plan_so_bond(nmat=2, nvol=5, B=5, nd=4, consumers=("vertex",),
                                         caps=SOBondCaps(1024.0, 1e-6, None), **extra)
                self.assertEqual(str(cm.exception), today)


    def test_refusal_options(self):
        """No empty line for an empty envelope; the clause without a key and
        with no nonzero off-site shell; a topology of another interaction."""
        compiled = compile_interaction(_chain_tables(), 1)
        topo = so_bond.so_bond_topology(compiled, CELL)
        rows = so_bond.shell_rows(compiled, topo, CELL)
        zero = tuple(r._replace(max_abs_v=0.0) for r in rows)
        kw = dict(nmat=2, nvol=5, B=5, nd=2, consumers=("vertex",),
                  caps=SOBondCaps(1024.0, 1e-6, None), cap_key=KEY)
        with self.assertRaises(ValueError) as cm:
            so_bond.plan_so_bond(shell_rows=zero, **kw)
        msg = str(cm.exception)
        self.assertNotIn("\n\n", msg)
        self.assertTrue(msg.endswith(
            "gamma_b  32000\nshells (shell  |dr|  max|v|):\n  1  1  0\n  2  2  0\n"
            "2 outer shell(s) declared only with zero coupling; max_shells = 0 removes "
            "them"), msg)
        with self.assertRaises(ValueError) as cm:
            so_bond.plan_so_bond(**kw)
        self.assertTrue(str(cm.exception).endswith(
            "gamma_b  32000\nshells (shell  |dr|  max|v|):\n  (no off-site shell)"))

        nn_only = compile_interaction({"CoulombInter": {((1, 0, 0), (0, 0)): 0.9,
                                                        ((-1, 0, 0), (0, 0)): 0.9}}, 1)
        foreign = [tuple(int(x) for x in d) for d in topo.delta_r if abs(int(d[0])) == 2]
        with self.assertRaises(ValueError) as cm:
            so_bond.shell_rows(nn_only, topo, CELL)
        self.assertEqual(str(cm.exception),
                         "shell_rows: the topology has displacement(s) {} that the compiled "
                         "interaction does not declare; topo must be a topology of compiled "
                         "(so_bond_topology)".format(foreign))
        with self.assertRaisesRegex(ValueError, "^shell_rows: compiled must be a "
                                                "CompiledInteraction"):
            so_bond.shell_rows(object(), topo, CELL)


if __name__ == "__main__":
    unittest.main()
