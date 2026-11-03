"""Archive members of the S_z-sector packed layout of code-slot arrays.

* ``so_bond.sector_map_members`` writes ``block_layout`` and, for the packed
  layout, the sector map members ``sz_sector_order``, ``sz_sector_offsets``
  and ``sz_sector_index``; ``so_bond.packed_members`` names the three
  sector blocks ``<name>@0``, ``<name>@+1``, ``<name>@-1``.
* ``so_bond.check_sector_map_members`` recomputes the map members from
  ``(B, nd)`` and refuses a file whose stored members differ, naming the
  first differing position.
* ``so_bond.read_slot_member`` converts between the layouts: a packed file
  is expanded for a dense consumer (off-sector entries exactly 0); a dense
  file is packed for a packed consumer only when every off-sector entry is
  exactly 0, otherwise it is refused with the measured leak of the member
  and its classification.

The expected maps come from the independent test-side enumeration
``tests.so_bond_oracles.ph_sector_lists``.
"""
import os
import tempfile
import types
import unittest

import numpy as np

from hwave.solver import so_bond
from hwave.solver.so_bond import (PackedArrays, check_sector_map_members, packed_members,
                                  read_slot_member, sector_map_members, sz_sector_map)
from tests.so_bond_fixtures import exact
from tests.so_bond_oracles import ph_sector_lists, restrict

NORB, B = 2, 3
ND_MODES = 2 * NORB
_REV = np.array([0, 2, 1], dtype=np.int64)
SM = sz_sector_map(types.SimpleNamespace(reverse=_REV), ND_MODES)
LISTS = ph_sector_lists(B, ND_MODES)
FILE = "chiq_bond.npz"
LABEL = "H-wave 1.2 (flex)"

_NUMERICAL = "numerical leakage from a non-canonical producer"
_INCOMPATIBLE = "incompatible Hamiltonian (produced from a different, S_z-breaking input)"


def _block_diagonal(lead, seed=0):
    """A random complex code-slot array (lead + (ND, ND)) whose off-sector
    entries are exactly 0."""
    rng = np.random.default_rng(seed)
    M = np.zeros(tuple(lead) + (SM.ND, SM.ND), dtype=np.complex128)
    for l in LISTS:
        n = len(l)
        M[..., l[:, None], l] = (rng.normal(size=tuple(lead) + (n, n))
                                 + 1j * rng.normal(size=tuple(lead) + (n, n)))
    return M


def _off_sector_mask():
    sector = np.empty(SM.ND, dtype=np.int64)
    for s, l in enumerate(LISTS):
        sector[l] = s
    return sector[:, None] != sector[None, :]


def _expected_map():
    """The map members of (B, nd) from the test-side enumeration."""
    sizes = [len(l) for l in LISTS]
    return {"sz_sector_order": np.array(["0", "+1", "-1"]),
            "sz_sector_offsets": np.array([0, sizes[0], sizes[0] + sizes[1], sum(sizes)],
                                          dtype=np.int64),
            "sz_sector_index": np.concatenate(LISTS).astype(np.int64)}


def _map_msg(member, i, stored, expected, B_=B, nd=ND_MODES):
    return ("'{}' is malformed or truncated: {} differs from the map of (B = {}, nd = {}) at "
            "position {} ({} vs {}). Regenerate it with {}.".format(
                FILE, member, B_, nd, i, stored, expected, LABEL))


def _leak_msg(leak, cls):
    return ("'%s' is dense with nonzero S_z-off-sector entries (measured leak = %.3e; "
            "packing requires exactly 0) although this run conserves S_z: %s. Regenerate it "
            "with this version." % (FILE, leak, cls))


def _packed_file(M):
    data = dict(sector_map_members(SM))
    data.update(packed_members("chi_w", PackedArrays.from_dense(M, SM, source="test")))
    return data


def _dense_file(M):
    data = dict(sector_map_members(None))
    data["chi_w"] = M
    return data


class TestArchiveMembers(unittest.TestCase):
    """Construction check of the member writers: names, dtypes, shapes and
    the values of the map members."""

    def test_packed_map_members(self):
        mem = sector_map_members(SM)
        self.assertEqual(sorted(mem), ["block_layout", "sz_sector_index", "sz_sector_offsets",
                                       "sz_sector_order"])
        self.assertEqual(mem["block_layout"].dtype, np.dtype("<U16"))
        self.assertEqual(mem["block_layout"].shape, ())
        self.assertEqual(str(mem["block_layout"][()]), "sz_packed")
        ref = _expected_map()
        self.assertEqual(mem["sz_sector_order"].dtype, np.dtype("<U2"))
        for k in ("sz_sector_offsets", "sz_sector_index"):
            self.assertEqual(mem[k].dtype, np.dtype(np.int64))
        for k, v in ref.items():
            self.assertTrue(np.array_equal(mem[k], v), k)
        self.assertEqual(mem["sz_sector_index"].shape, (SM.ND,))
        self.assertEqual(mem["sz_sector_offsets"].shape, (4,))

    def test_members_are_copies(self):
        mem = sector_map_members(SM)
        mem["sz_sector_index"][0] = -7                   # writable, not the map's own array
        self.assertEqual(int(SM.index[0]), int(LISTS[0][0]))

    def test_dense_members(self):
        mem = sector_map_members(None)
        self.assertEqual(list(mem), ["block_layout"])
        self.assertEqual(str(mem["block_layout"][()]), "dense")
        self.assertEqual(mem["block_layout"].dtype, np.dtype("<U16"))

    def test_sector_map_type_refused(self):
        with self.assertRaisesRegex(ValueError, exact(
                "sector_map_members: sector_map must be a SectorMap or None, got dict")):
            sector_map_members({})

    def test_packed_members(self):
        P = PackedArrays.from_dense(_block_diagonal((2,)), SM, source="test")
        mem = packed_members("Gamma_bond", P)
        self.assertEqual(list(mem), ["Gamma_bond@0", "Gamma_bond@+1", "Gamma_bond@-1"])
        for k, b in zip(mem, P.blocks):
            self.assertIs(mem[k], b)

    def test_packed_members_refusals(self):
        with self.assertRaisesRegex(ValueError, exact(
                "packed_members: packed must be a PackedArrays, got ndarray")):
            packed_members("chi_w", np.zeros((2, 2)))
        P = PackedArrays.zeros(SM, (1,))
        for bad in ("", "chi@w", 3):
            with self.subTest(name=bad), self.assertRaisesRegex(ValueError, exact(
                    "packed_members: name must be a non-empty string without '@', got "
                    "{!r}".format(bad))):
                packed_members(bad, P)


class TestArchiveLayout(unittest.TestCase):
    """P-arch: packed write/read is bitwise; packed -> dense is exact with
    off-sector zeros; dense -> packed accepts exact zeros and refuses
    nonzero off-sector entries with the measured leak and its
    classification."""

    def test_npz_round_trip_packed(self):
        M = _block_diagonal((2, 3), seed=1)
        P = PackedArrays.from_dense(M, SM, source="test")
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, FILE)
            np.savez(path, **_packed_file(M))
            with np.load(path) as data:
                self.assertEqual(check_sector_map_members(data, B, ND_MODES, file=FILE,
                                                          producer_label=LABEL), "sz_packed")
                for s, k in enumerate(("chi_w@0", "chi_w@+1", "chi_w@-1")):
                    self.assertEqual(data[k].dtype, np.complex128)
                    self.assertEqual(data[k].tobytes(), P.blocks[s].tobytes())
                R = read_slot_member(data, "chi_w", want="sz_packed", sector_map=SM, file=FILE)
                D = read_slot_member(data, "chi_w", want="dense", sector_map=SM, file=FILE)
        self.assertIsInstance(R, PackedArrays)
        self.assertIs(R.sector_map, SM)
        for s in range(3):
            self.assertEqual(R.blocks[s].tobytes(), P.blocks[s].tobytes())
        self.assertEqual(D.tobytes(), M.tobytes())

    def test_packed_to_dense_exact(self):
        M = _block_diagonal((2, 3), seed=2)
        D = read_slot_member(_packed_file(M), "chi_w", want="dense", sector_map=SM, file=FILE)
        self.assertEqual((D.dtype, D.shape), (np.dtype(np.complex128), M.shape))
        self.assertTrue(np.array_equal(D, M))
        self.assertTrue(np.all(D[..., _off_sector_mask()] == 0))

    def test_dense_to_packed_exact_zeros(self):
        M = _block_diagonal((2, 3), seed=3)
        R = read_slot_member(_dense_file(M), "chi_w", want="sz_packed", sector_map=SM,
                             file=FILE)
        self.assertIsInstance(R, PackedArrays)
        for s, ref in enumerate(restrict(M, LISTS)):
            self.assertTrue(np.array_equal(R.blocks[s], ref))
            self.assertTrue(R.blocks[s].flags.c_contiguous)

    def test_same_layout_pass_through(self):
        M = _block_diagonal((3,), seed=4)
        D = read_slot_member(_dense_file(M), "chi_w", want="dense", sector_map=None, file=FILE)
        self.assertEqual((D.dtype, D.shape), (np.dtype(np.complex128), M.shape))
        self.assertTrue(np.array_equal(D, M))
        R = read_slot_member(_packed_file(M), "chi_w", want="sz_packed", sector_map=SM,
                             file=FILE)
        for s, ref in enumerate(restrict(M, LISTS)):
            self.assertEqual((R.blocks[s].dtype, R.blocks[s].shape),
                             (np.dtype(np.complex128), ref.shape))
            self.assertTrue(np.array_equal(R.blocks[s], ref))

    def test_unbatched_member(self):
        M = _block_diagonal((), seed=5)
        D = read_slot_member(_packed_file(M), "chi_w", want="dense", sector_map=SM, file=FILE)
        self.assertTrue(np.array_equal(D, M))
        R = read_slot_member(_dense_file(M), "chi_w", want="sz_packed", sector_map=SM,
                             file=FILE)
        for s, ref in enumerate(restrict(M, LISTS)):
            self.assertTrue(np.array_equal(R.blocks[s], ref))

    def test_dense_to_packed_refusals(self):
        # off-sector (0, +1) entry in the LAST slice: the leak is that of the whole member
        I, J = int(LISTS[0][3]), int(LISTS[1][2])
        for eps, cls in ((1e-14, _NUMERICAL), (1e-3, _INCOMPATIBLE)):
            M = _block_diagonal((2, 3), seed=6) / 10.0
            M[0, 0, int(LISTS[0][0]), int(LISTS[0][1])] = 2.0       # max |M| in slice 0
            M[1, 2, I, J] = eps
            leak = eps / 2.0
            with self.subTest(eps=eps), self.assertRaisesRegex(ValueError,
                                                                exact(_leak_msg(leak, cls))):
                read_slot_member(_dense_file(M), "chi_w", want="sz_packed", sector_map=SM,
                                 file=FILE)

    def test_dense_to_packed_non_finite(self):
        M = _block_diagonal((2,), seed=7)
        M[1, int(LISTS[2][0]), int(LISTS[0][1])] = np.nan
        with self.assertRaisesRegex(ValueError, exact(_leak_msg(float("nan"), _INCOMPATIBLE))):
            read_slot_member(_dense_file(M), "chi_w", want="sz_packed", sector_map=SM, file=FILE)
        # a finite off-sector entry with a non-finite on-sector entry in another slice
        M = _block_diagonal((2,), seed=8)
        M[0, int(LISTS[0][0]), int(LISTS[0][0])] = np.nan
        M[1, int(LISTS[1][0]), int(LISTS[0][1])] = 1e-3
        with self.assertRaisesRegex(ValueError, exact(_leak_msg(float("nan"), _INCOMPATIBLE))):
            read_slot_member(_dense_file(M), "chi_w", want="sz_packed", sector_map=SM, file=FILE)


class TestSectorMapCheck(unittest.TestCase):
    """P-arch: a stored sector map that differs from the map of (B, nd) is
    refused at its first differing position (permuted, truncated, extended,
    reordered, missing, wrong dtype or shape)."""

    def _check(self, data, B_=B, nd=ND_MODES):
        return check_sector_map_members(data, B_, nd, file=FILE, producer_label=LABEL)

    def _refused(self, data, msg, **kw):
        with self.assertRaisesRegex(ValueError, exact(msg)):
            self._check(data, **kw)

    def test_accepted(self):
        self.assertEqual(self._check(sector_map_members(SM)), "sz_packed")
        self.assertEqual(self._check(sector_map_members(None)), "dense")
        # the map does not depend on reverse
        other = sz_sector_map(types.SimpleNamespace(reverse=np.arange(B)), ND_MODES)
        self.assertEqual(self._check(sector_map_members(other)), "sz_packed")

    def test_permuted_index(self):
        data = sector_map_members(SM)
        idx = data["sz_sector_index"]
        idx[[5, 9]] = idx[[9, 5]]
        ref = _expected_map()["sz_sector_index"]
        self._refused(data, _map_msg("sz_sector_index", 5, idx[5], ref[5]))

    def test_truncated_index(self):
        data = sector_map_members(SM)
        data["sz_sector_index"] = data["sz_sector_index"][:-1]
        ref = _expected_map()["sz_sector_index"]
        self._refused(data, _map_msg("sz_sector_index", SM.ND - 1, "end of member",
                                     ref[SM.ND - 1]))

    def test_extended_index(self):
        data = sector_map_members(SM)
        data["sz_sector_index"] = np.append(data["sz_sector_index"], np.int64(3))
        self._refused(data, _map_msg("sz_sector_index", SM.ND, 3, "end of map"))

    def test_map_of_another_bond_count(self):
        # a file written for B = 3 checked against a run with B = 2
        nd = ND_MODES
        sm2 = ph_sector_lists(2, nd)
        n0 = len(sm2[0])
        self._refused(sector_map_members(SM),
                      _map_msg("sz_sector_offsets", 1, len(LISTS[0]), n0, B_=2), B_=2)

    def test_reordered_sector_order(self):
        data = sector_map_members(SM)
        data["sz_sector_order"] = np.array(["0", "-1", "+1"])
        self._refused(data, _map_msg("sz_sector_order", 1, "'-1'", "'+1'"))

    def test_offsets(self):
        data = sector_map_members(SM)
        data["sz_sector_offsets"][2] += 1
        ref = _expected_map()["sz_sector_offsets"]
        self._refused(data, _map_msg("sz_sector_offsets", 2, ref[2] + 1, ref[2]))

    def test_missing_member(self):
        data = sector_map_members(SM)
        del data["sz_sector_offsets"]
        self._refused(data, _map_msg("sz_sector_offsets", 0, "missing", 0))

    def test_dtype_and_shape(self):
        for dt in (np.float64, np.complex128, np.bool_, object):
            with self.subTest(dtype=dt):
                data = sector_map_members(SM)
                data["sz_sector_index"] = data["sz_sector_index"].astype(dt)
                self._refused(data, _map_msg("sz_sector_index", 0, "dtype {}".format(
                    np.dtype(dt)), "dtype int64"))
        data = sector_map_members(SM)
        data["sz_sector_order"] = data["sz_sector_order"].astype("S2")
        self._refused(data, _map_msg("sz_sector_order", 0, "dtype |S2", "dtype <U2"))
        data = sector_map_members(SM)
        data["sz_sector_index"] = data["sz_sector_index"].astype(np.uint64)
        data["sz_sector_index"][0] = np.iinfo(np.uint64).max
        self._refused(data, _map_msg("sz_sector_index", 0, "dtype uint64", "dtype int64"))
        data = sector_map_members(SM)
        data["sz_sector_index"] = data["sz_sector_index"].reshape(B, -1)
        self._refused(data, _map_msg("sz_sector_index", 0, "shape (3, 16)", "shape (48,)"))

    def test_value_based_dtypes(self):
        """Any integer dtype holds the integer members and any unicode width
        or byte order the sector order (``<U2`` as written, ``>U2``, ``<U5``);
        the values decide."""
        for member, dts in (("sz_sector_index", (np.int32, np.uint16, ">i8", np.int64)),
                            ("sz_sector_offsets", (np.int16, np.uint32)),
                            ("sz_sector_order", ("<U2", ">U2", "<U5"))):
            for dt in dts:
                with self.subTest(member=member, dtype=dt):
                    data = sector_map_members(SM)
                    data[member] = data[member].astype(dt)
                    self.assertEqual(data[member].dtype, np.dtype(dt))
                    self.assertEqual(self._check(data), "sz_packed")
        # a wider unicode member is compared untruncated
        data = sector_map_members(SM)
        data["sz_sector_order"] = np.array(["0", "+1x", "-1"], dtype="<U3")
        self._refused(data, _map_msg("sz_sector_order", 1, "'+1x'", "'+1'"))
        # a value difference in another integer dtype is reported by value
        data = sector_map_members(SM)
        data["sz_sector_offsets"] = data["sz_sector_offsets"].astype(np.int32)
        data["sz_sector_offsets"][2] += 1
        ref = _expected_map()["sz_sector_offsets"]
        self._refused(data, _map_msg("sz_sector_offsets", 2, ref[2] + 1, ref[2]))

    def test_block_layout(self):
        data = sector_map_members(SM)
        data["block_layout"] = np.array("packed")
        self._refused(data, _map_msg("block_layout", 0, "'packed'", "'dense' or 'sz_packed'"))
        del data["block_layout"]
        self._refused(data, _map_msg("block_layout", 0, "missing", "'dense' or 'sz_packed'"))

    def test_argument_refusals(self):
        for bad in (0, -1, 2.0, True):
            with self.subTest(B=bad), self.assertRaisesRegex(ValueError, exact(
                    "check_sector_map_members: B must be a positive integer, got "
                    "{!r}".format(bad))):
                self._check(sector_map_members(SM), B_=bad)
        with self.assertRaisesRegex(ValueError, exact(
                "sz_sector_map: nd must be a positive even integer (nd = 2*norb_phys), "
                "got 3")):
            self._check(sector_map_members(SM), nd=3)


class TestReadSlotMemberRefusals(unittest.TestCase):
    """Construction check of the reader's refusals: layouts, maps and
    missing members."""

    def test_want(self):
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: want must be 'dense' or 'sz_packed', got 'auto'")):
            read_slot_member(_dense_file(_block_diagonal((1,))), "chi_w", want="auto",
                             sector_map=SM, file=FILE)

    def test_sector_map_required(self):
        M = _block_diagonal((1,))
        for data, want in ((_dense_file(M), "sz_packed"), (_packed_file(M), "dense")):
            with self.subTest(want=want), self.assertRaisesRegex(ValueError, exact(
                    "read_slot_member: converting 'chi_w' of '{}' needs the run's sector "
                    "map, got None".format(FILE))):
                read_slot_member(data, "chi_w", want=want, sector_map=None, file=FILE)

    def test_unknown_layout(self):
        data = _dense_file(_block_diagonal((1,)))
        data["block_layout"] = np.array("packed")
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: '{}' has block_layout 'packed', expected 'dense' or "
                "'sz_packed'; check it with check_sector_map_members first".format(FILE))):
            read_slot_member(data, "chi_w", want="dense", sector_map=SM, file=FILE)

    def test_map_differs_from_run(self):
        data = _packed_file(_block_diagonal((1,)))
        data["sz_sector_index"] = data["sz_sector_index"][::-1].copy()
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: the sector map stored in '{}' differs from this run's "
                "map (B = 3, nd = 4); check it with check_sector_map_members "
                "first".format(FILE))):
            read_slot_member(data, "chi_w", want="dense", sector_map=SM, file=FILE)

    def test_missing_member(self):
        data = _packed_file(_block_diagonal((1,)))
        del data["chi_w@+1"]
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: '{}' has no member 'chi_w@+1'".format(FILE))):
            read_slot_member(data, "chi_w", want="dense", sector_map=SM, file=FILE)
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: '{}' has no member 'Gamma_bond'".format(FILE))):
            read_slot_member(_dense_file(_block_diagonal((1,))), "Gamma_bond", want="dense",
                             sector_map=None, file=FILE)

    def test_dense_shape(self):
        data = _dense_file(np.zeros((2, SM.ND, SM.ND - 1), dtype=np.complex128))
        with self.assertRaisesRegex(ValueError, exact(
                "read_slot_member: 'chi_w' of '{0}' must end in shape (ND, ND) = ({1}, {1}), "
                "got shape (2, {1}, {2})".format(FILE, SM.ND, SM.ND - 1))):
            read_slot_member(data, "chi_w", want="sz_packed", sector_map=SM, file=FILE)


if __name__ == "__main__":
    unittest.main()
