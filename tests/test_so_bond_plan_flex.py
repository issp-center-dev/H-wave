"""The FLEX consumer ``flex_scf`` of the spin-orbital bond admission plan
(``so_bond.plan_so_bond``) and the SCF rows it shares with the non-SO
bond FLEX (``flex_bond.flex_scf_rows``, ``flex_so_bond.flex_so_bond_rows``).

* ``flex_scf`` closes over ``vertex``; its rows are the resident store slot
  ``chibar`` (``W`` written in place over it by default, a second resident
  ``w_slot`` with ``w_in_place=False``), the bubble's pair block, the
  dressing batch with the static pole slice, and the W-conversion
  workspace; the plan's peak is ``1.25 x`` the largest phase total, the
  resident-only consumer phases ``vertex`` and ``flex_scf`` included.
* ``resident_rows`` adds named rows alive in every phase; a name that
  repeats a row of the plan is refused.
* ``remedy_note`` is inserted after the head sentence of the over-cap
  refusal; without it every text is the one of a92b533d.
* Every parameter set of the existing plan and estimate tests gives the
  rows, peak and reason (estimate: the whole returned dict) of a92b533d:
  the digests below were recorded at a92b533d (numpy 2.4.6 and 1.26.4
  give the same digests) by running ``existing_plan_cases`` and
  ``estimate_cases`` unchanged against that revision.
"""
import hashlib
import itertools
import unittest

import numpy as np

from hwave.solver import flex_bond, so_bond
from hwave.solver.so_bond import SOBondCaps, plan_so_bond

E = 16
NOCAP = SOBondCaps(None, None, None)
SMALL = SOBondCaps(1.0e3, 1.0e3 / 2 ** 30, None)
KEY = "[mode.param] longitudinal_bond_memory_cap_gb"


def _flex_plan(nmat, nvol, B, nd, layout="dense", caps=NOCAP, **kw):
    return plan_so_bond(nmat=nmat, nvol=nvol, B=B, nd=nd, consumers=("flex_scf",), caps=caps,
                        layout=layout, layout_reason=None if layout == "dense" else _PACKED,
                        **kw)


def _so_rows():
    from hwave.solver.flex_so_bond import flex_so_bond_rows
    return flex_so_bond_rows


# The grid generators below use only the a92b533d API; they are the ones the
# digests were recorded with.
_PACKED = "so bond layout: sz_packed (S_z conserved)"
#: the (nmat, nvol, B, nd) of every plan_so_bond call of test_so_bond_plan_static.py
#: and test_so_bond_packed_memory.py at a92b533d
_PLAN_SIZES = ((4, 5, 5, 2), (4, 15, 5, 2), (4, 5, 3, 4), (4, 15, 5, 4), (4, 8, 3, 2),
               (8, 16, 5, 2), (8, 15, 5, 4), (2, 15, 5, 4), (2, 16, 1, 2), (2, 15, 5, 2),
               (2, 8, 3, 2), (2, 15, 1, 4), (2, 15, 3, 4), (2, 72, 3, 4), (256, 1024, 5, 4))
_PLAN_CONSUMERS = (("vertex",), ("chi_static",), ("chi_dynamic",), ("pairing_static",),
                   ("chi_dynamic", "pairing_static"))
#: the static extra rows of those calls: (nmat, nvol, B, nd) -> spatial shape
_STATIC_SHAPES = {(8, 16, 5, 2): (4, 4, 1), (2, 16, 1, 2): (4, 4, 1),
                  (256, 1024, 5, 4): (32, 32, 1)}


def _canon(v):
    """A numpy-version independent canonical form (numpy scalars -> Python)."""
    if isinstance(v, dict):
        return tuple(sorted((str(k), _canon(x)) for k, x in v.items()))
    if isinstance(v, (list, tuple)):
        return tuple(_canon(x) for x in v)
    if isinstance(v, np.ndarray):
        return ("array", str(v.dtype), v.shape, _canon(v.tolist()))
    if isinstance(v, np.generic):
        return v.item()
    return v


def _digest(outcome):
    return hashlib.sha256(repr(_canon(outcome)).encode()).hexdigest()[:16]


def _plan_outcome(kw):
    try:
        p = so_bond.plan_so_bond(**kw)
    except ValueError as e:
        return ("ValueError", str(e))
    return (p.rows, float(p.peak_host), p.reason, p.nb, p.dress_capacity, p.slots, p.aliases)


def existing_plan_cases():
    """(label, plan_so_bond kwargs): the parameter sets of the existing plan tests."""
    nocap = SOBondCaps(None, None, None)
    out = []
    for size, cons, layout in itertools.product(_PLAN_SIZES, _PLAN_CONSUMERS,
                                                ("dense", "sz_packed")):
        nmat, nvol, B, nd = size
        base = dict(nmat=nmat, nvol=nvol, B=B, nd=nd, consumers=cons, caps=nocap,
                    layout=layout, layout_reason=None if layout == "dense" else _PACKED)
        dyn = "chi_dynamic" in cons
        for keep, nb in itertools.product((False, True) if dyn else (False,), (None, 2)):
            kw = dict(base, nb=nb)
            if keep:
                kw.update(keep_chi=True, vertex_batch=True)
            out.append(("{} {} {} keep={} nb={}".format(size, cons, layout, keep, nb), kw))
        if cons == ("pairing_static",) and size in _STATIC_SHAPES:
            rows = so_bond.static_so_bond_rows(
                nmat=nmat, nvol=nvol, B=B, nd=nd, layout=layout, solver_mode="eigenvalue",
                num_eigenvalues=10, max_iter=1000, spatial_shape=_STATIC_SHAPES[size])
            out.append(("{} static {}".format(size, layout), dict(base, extra_rows=rows)))
            small = SOBondCaps(1.0e3, 1.0e3 / 2 ** 30, None)
            out.append(("{} static {} over cap".format(size, layout),
                        dict(base, extra_rows=rows, caps=small)))
            out.append(("{} static {} over cap keyed".format(size, layout),
                        dict(base, extra_rows=rows, caps=small, cap_key="[eliashberg] k",
                             envelope=so_bond.envelope_text(layout),
                             max_shells_key="[eliashberg] max_shells")))
        if cons == ("chi_dynamic",) and size == (8, 16, 5, 2):
            small = SOBondCaps(1.0e3, 1.0e3 / 2 ** 30, None)
            out.append(("{} dynamic {} over cap".format(size, layout), dict(base, caps=small)))
            out.append(("{} dynamic {} over cap keyed".format(size, layout),
                        dict(base, caps=small, cap_key="[mode.param] k",
                             envelope=so_bond.envelope_text(layout))))
    return out


def _estimate_outcome(kw):
    try:
        return ("ok", flex_bond.estimate_bond_memory(**kw))
    except ValueError as e:
        return ("ValueError", str(e))


def estimate_cases():
    """(label, estimate_bond_memory kwargs): the TestEstimate grid of
    test_flex_bond_memory.py widened over every SCF-row argument."""
    base = dict(nmat=8, nvol=16, norb=2, B=3, depth=0, output_full=False, split_seed=False,
                n_types=1, freq_batch=None, cap_gb=8.0, mixing="linear")
    out = []
    for depth, mixing, of, ss, nt, hg, wip, dev in itertools.product(
            (0, 1, 4), ("linear", "anderson"), (False, True), (False, True), (0, 2),
            (False, True), (False, True), (None, int(2e8))):
        out.append(("grid", dict(base, depth=depth, mixing=mixing, output_full=of,
                                 split_seed=ss, n_types=nt, half_grid=hg, w_in_place=wip,
                                 device_available=dev)))
    for nmat, nvol, norb, B in ((4, 9, 1, 5), (16, 16, 3, 9), (32, 64, 2, 13)):
        for fb in (None, 1, 3):
            out.append(("size", dict(base, nmat=nmat, nvol=nvol, norb=norb, B=B,
                                     freq_batch=fb, cap_gb=1e3)))
    # the TestEstimate calls themselves (cap boundary and refusals)
    big = flex_bond.estimate_bond_memory(**base)
    cap = big["peak"] / 2 ** 30
    for over in (dict(cap_gb=cap * (1 + 1e-9)), dict(cap_gb=cap * (1 + 1e-9), freq_batch=2),
                 dict(cap_gb=cap * (1 + 1e-9), freq_batch=8), dict(cap_gb=1e-6),
                 dict(depth=4, mixing="anderson", output_full=True, split_seed=True,
                      n_types=2), dict(factor_bytes=12345, w2_triples=7),
                 dict(guard_freqs="static", device_available=int(2e8)),
                 dict(guard_enabled=False, device_available=int(2e8)),
                 dict(store="auto", device_available=int(2e8), pairing=True),
                 dict(store="device", device_available=int(2e9))):
        out.append(("estimate {}".format(sorted(over)), dict(base, **over)))
    return out


#: recorded at a92b533d (see the module docstring)
_PLAN_DIGESTS = (
    'b241fe91c8648923', 'd973d8aee9b0f053', '365611d74795942e', '674f9f0b12ccdfde',
    '0b2948d455fe20f9', 'bf1659c6967585c9', '141f549fcb38d237', '37d1c3c002e77127',
    '52e4a58274c2d814', 'd36db021ea309a80', '6871c439dcf2dc92', 'd94c727512440556',
    '13948083debb913b', '530236441d37f467', 'cedf01ac0e53f5db', 'b7af866e182cb47b',
    'a4c1a547d3033efd', 'f3f87ed7ef4ed1c7', '21a0dcd4147f1616', '76d92d1098f46539',
    '02f09c5a880c209f', '447f3ee767a2ef29', 'cad3855d7032a303', '4cbd4655b6ba23b4',
    'bcc0c912a494975e', '9e4948595b3bbdba', '123b76f2c4d2efbc', '271cde28f2830a48',
    '7012ba5689ade22f', '60c06d2a21845e69', '891045118377f804', '0aea6d6a566c5598',
    'd55e4356702f811d', 'ca8ebd40b6505cd5', '4b93bc81949e06a7', '7a45bd116fde0c03',
    '325c7fff3caabed6', 'da25f64c91f44c4b', '5437830acc5bc367', 'e396381b15a1d10b',
    '1149af7855de44f5', '7fb9d9050925a270', '06cf93aa3f3aaa91', 'ef5b0ea83ffe3585',
    '1e14733c3c619b45', 'dcd2ceac3ec84f27', 'b602830610add6a7', '8d1e85f998ec42ed',
    'a2ea599258bdfd1a', '92d627825e916a2a', '47e83906b373ca7b', 'd39b5e9dfc041fef',
    '7b6f841382b0d9c4', '69269da611a62046', '458997410c66c527', '4c160090e06e9f9b',
    '842324f6431c4069', '43739981dbe780fc', '9b30ef4d675a9057', '54b413af1af802c7',
    '243cc8b8737a9d47', '44686b1843766ba3', 'edb0b1fd057181f5', '5d36c3e426c0d727',
    'c6e7778ed5338361', 'd39895ff2c0978f0', '37f2a1d9f6dff317', '1170f745a3ec4e7d',
    '84bbc20d7c8c3c7c', '618e3c89674989c3', 'b5e2970b030db6e9', 'b738753a0b6d21f7',
    '5d9bc1cd6935524e', '6e594cf06b861340', '082e4df3e5d65fc4', 'c27c7d13c6aded61',
    'afa9ef0a875e3f9a', 'd57600472f1aee87', '6041f2489c6d547c', '2a9d7731cbed5d2e',
    'b52c46ba438136cb', '2e6129409298e43b', '30604f7d1e94fdb4', '79abd0e7f89bcbeb',
    '742a178c91e0cb62', '0a13f01962318f82', 'fde0df6108eecf1f', '34ae4b09101f01fb',
    '281841aad6866eed', 'd55f86b5bf649ea1', 'cc094150eff0438a', '9fbeda84fb019320',
    '0e658231d0313088', '470a0e92dcb212ca', '7d0508507f7f7743', '58c45cb9ccdac095',
    'f0586a65f978e0b4', '3af2c6a556d71141', '7674b7f343ecdbec', 'c5bba6706636f602',
    'b418e051a72940e5', 'f337e966d28d28ab', '5302f5d67f848932', 'a94d7193b31c58e0',
    'c326fea2f1d73185', '5115a660be2113bc', '1f9a4d680b342b5c', '5d0cf8d6a9ffaa2a',
    'de52698e04896526', '4abf7daa6a771fe5', 'f0808ba84182e5a0', 'f91bc6dc1e9a0e90',
    'f43725c83e57a296', '0118bdc2c47fce01', 'de812ac8019873d5', '9ea080fd0cf7fa25',
    '87d0a08c167b7fb8', '87341916e23ec1ad', 'd157660c33d942db', 'b852e5412ffb7a74',
    'a7e60f99670fc3fb', 'b2cfec2098ea1b4c', 'fb166e616a2a1583', '52686e9f6a3ff442',
    '9dd0d4b94906f52a', '745703adb11febc6', '1448d27d624ae536', 'd03bb9852828e08f',
    '872d408dc0035268', '74ae0375b258ffec', '0fa320266c31db85', '6740ddd8a3480a46',
    '980affdb40a46361', '4d84e794e6728e24', '192e3961ac556d44', 'c6c9808afae0c9e7',
    'd63b03b5a57a18cc', 'f2189df93a1a7813', '1eac51a727c4aa5a', '8d971c028e7cd359',
    'a5e22ebd377ffdff', '3d3d027606f9e32f', '02ae0ccfe1b71af8', 'a30b6dd0937f6198',
    'bae9ad9caaedf028', '7752553ee229f440', '22f0c29248a3652e', '96513a2fd72926c4',
    'ce6d23a5689d7005', '78a69101bfa8f455', '2d6f005d74374b3d', 'f0dd96647e385ec1',
    '8e2951c40d9293c0', '54d0794f0296dd9d', '6d0a30aca3472400', 'f5304cf87c10be97',
    '6bb68f43ec4ca841', '4b58e8566d2036b9', 'a2f196cace0fe9df', '82182418c5eb0f1d',
    '8dc2d3f151bf252b', '0942297d920dc826', '222a73ea1c4f89e9', 'f06f53a0ad6446e6',
    'c213207d63634824', '4bc80e5f418d0af9', 'eba1015b3b81ebfb', '8c0f83d4cb24d8c3',
    'dc4e54333ad67e8d', 'a8d76e88f8a70d3b', 'c501a95037dd942b', 'b8050497ba99a21f',
    '2ecd95d0f8487d44', 'd0a3721fb9bd2072', '4ddb93d62b20c759', '35cfc10a0808df30',
    'f8ebbb8b7fc70ab5', 'd444eb4025d70cfc', '11427f1eb1d9db9b', '2bc8b982779dfd8c',
    '5ccd9b133f55e260', '7df803a9dff429ca', '38518b1e6e2ef753', '2e25b1dfe96de070',
    'b4f09358abc22a61', 'b394b4fca1669015', '24a97966ca2a2151', '368fef77ed38198a',
    '668808aa9b68fc26', 'ace2dc67145939f2', 'd8ba282ccd8d51b5', '311cfec8a6d9bd83',
    '662cf9ce06a236a4', '9bade4aa7fdb2e7c', '98ae51760e42370c', '30400b61d9d7a1bf',
    '6cd8bf61c7ed492a', 'c43bef456cd761dd', '27e70a391b4ad35b', '83e4bf77b8884214',
    'd361951b7ee4eaa2', '55574f643532c2e8', 'a02a8891e9dba45a', '0e0b405556612dba',
    'c43db2f600bc22b0', '4cb8072acae77843', '77c978fa1f7983e5', '77c978fa1f7983e5',
    '122f73f0bbcd8809', '122f73f0bbcd8809', '741e7e201fc8e6e1', '741e7e201fc8e6e1',
    'b258d252c437bf11', 'b258d252c437bf11', '2ff74914499031a8', '2ff74914499031a8',
    'ac57cf0a45748907', 'ac57cf0a45748907', '9fb5cd6d2dd6bba5', '9fb5cd6d2dd6bba5',
    'd84e73d650c5bfdf', 'd84e73d650c5bfdf', '1dad051694afcfa6', '1dad051694afcfa6',
    'cc41775a790d3103', 'cc41775a790d3103', '80ea0afc2462665d', '80ea0afc2462665d',
    'c918f0ee413ba2f0', 'c918f0ee413ba2f0', '8b0d2811d4a257cc', '8b0d2811d4a257cc',
    'd6b44cbd34420505', 'd6b44cbd34420505', '6858a1a929da0cbe', '6858a1a929da0cbe',
    '6310d1eccbd2f2a9', '6310d1eccbd2f2a9', 'd97e14fe147d5851', 'd97e14fe147d5851',
    'ed40f17b2f867e35', 'ed40f17b2f867e35', '5511fbce201f3e72', '5511fbce201f3e72',
    'a66e4286532eb7fc', 'a66e4286532eb7fc', 'aacae35f45819a5a', 'aacae35f45819a5a',
    '42d3f071d69a3418', '42d3f071d69a3418', '57629e7d9383febd', '57629e7d9383febd',
    'cbd19d54b60cf81c', '529fd11567643cc7', 'ea5ecb466aca5d69', '88f8718f806d2923',
    '88f8718f806d2923', 'd0f5bbe30cd7aee3', '5c74a0587ab63aea', '73e7fa50aa259233',
    '115ce5687e500199', '115ce5687e500199', 'a5f7660ceeb4f6d0', 'a5f7660ceeb4f6d0',
    '5b0b93b1ed02fe4c', '5b0b93b1ed02fe4c', '2a273f8a898f532d', '2a273f8a898f532d',
    '7db8905406af5681', '7db8905406af5681', '0b8bfbb78b86bd24', '0b8bfbb78b86bd24',
    '8b5acddf19b2a9fa', '8b5acddf19b2a9fa', 'fa46a735c46c709c', 'fa46a735c46c709c',
    'b178c74a756f69e5', 'b178c74a756f69e5', 'e0a9da4b38e01fa8', 'e0a9da4b38e01fa8',
    '9477086f96ad4c04', '9477086f96ad4c04', '4ee5eb5104ca42bc', '4ee5eb5104ca42bc',
    '894877dbd27327fb', '894877dbd27327fb', '6b16de7107cdc91f', '6b16de7107cdc91f',
    '6c8bd599ba5ef82c', '6c8bd599ba5ef82c', '650b6c64afe4f542', '650b6c64afe4f542',
    'bec37a2dc6d3993b', 'bec37a2dc6d3993b', '3fb725b62ae045e1', '3fb725b62ae045e1',
    '4f2303e797d3a7d3', '4f2303e797d3a7d3', '9a5c12f2714fa7d2', '9a5c12f2714fa7d2',
    '9e32bc982abe3d9f', '9e32bc982abe3d9f', 'e6261de25a92e69e', 'e6261de25a92e69e',
    'c97b9561dd1dfc0e', 'c97b9561dd1dfc0e', '8971776b5f0fa32d', '8971776b5f0fa32d',
    '0ccb45f1843f0210', '0ccb45f1843f0210', '81221516d5873818', '81221516d5873818',
    '82cd8ccf633356d4', '82cd8ccf633356d4', 'ab5375e6fb094a6a', 'ab5375e6fb094a6a',
    '2b29ccda5aa831f4', '2b29ccda5aa831f4', '60beae384f8f469e', '60beae384f8f469e',
    '5f3379fc48996ead', '5f3379fc48996ead', '6420748216e06a53', '6420748216e06a53',
    '55502c213e7fcd90', '55502c213e7fcd90', '9d5ae334de4338f9', '9d5ae334de4338f9',
    '60b05227401a336d', '60b05227401a336d', 'eca2d9fbfe6a37b3', 'eca2d9fbfe6a37b3',
    '7ba707e2ee227399', '7ba707e2ee227399', 'c6d17af66c1a4b87', 'c6d17af66c1a4b87',
    '0f0a280b2602532e', '0f0a280b2602532e', 'bca3e7964d107a53', 'bca3e7964d107a53',
    '075382b9a969e54e', '075382b9a969e54e', '65ea25497b4e6ced', '65ea25497b4e6ced',
    '26e02d33e0324674', '26e02d33e0324674', 'b0336e315995c5e7', 'b0336e315995c5e7',
    '78b9e339f6d2aa01', '78b9e339f6d2aa01', 'a7cc3da60751d07d', 'a7cc3da60751d07d',
    'fe088097c449a7cd', 'fe088097c449a7cd', 'b5223469844ac794', 'b5223469844ac794',
    'c5cdf8e004768425', 'c5cdf8e004768425', '0ac249d78f2234f8', '0ac249d78f2234f8',
    'c164622a86cbb5f1', 'c164622a86cbb5f1', 'd5eba82bd656a329', 'd5eba82bd656a329',
    '14756a221d53df06', '14756a221d53df06', '9b28d3c9cecb3c91', '9b28d3c9cecb3c91',
    'c36b57e3037d9376', 'c36b57e3037d9376', '01df729f011db681', '01df729f011db681',
    '880ab7c02b884451', '880ab7c02b884451', '64919567d43ea291', '64919567d43ea291',
    '4f16856daaea0328', '4f16856daaea0328', 'b6f788cbc3cc9411', 'b6f788cbc3cc9411',
    'ad47b7687978b48f', 'ad47b7687978b48f', 'c05725dc213e384b', 'c05725dc213e384b',
    '01f78df5614d5f29', '01f78df5614d5f29', 'ed9600700e0d32cf', 'ed9600700e0d32cf',
    '60f1866bbaa6b286', '60f1866bbaa6b286', 'bb6acca0fd7f0e59', 'bb6acca0fd7f0e59',
    '91fbefe85f1e7cc8', '91fbefe85f1e7cc8', '893ffba44c96144e', '893ffba44c96144e',
    '4cb5149939bddd88', '4cb5149939bddd88', '6a4fc9d8b4c33255', '6a4fc9d8b4c33255',
    'e6a1b943ace449ee', 'e6a1b943ace449ee', 'ff3900525815f7b0', 'ff3900525815f7b0',
    'd1629db14296ca27', 'd1629db14296ca27', '0a7377ec155f6168', '0a7377ec155f6168',
    '5a0426c4936604eb', '5a861f947a940794', 'cc671871bc10df54', '6f20c476546e3819',
    '9b2ee35650b1ad79', '83e6f6bf79725030', '9da54a923f150605', 'e619215d90653786',
    '32aad6ec0e5b60e4', 'b014b26242790014', 'e2b09fc5f9a07aa9', '55eec8c26a5798b8',
    '47bb0cefb690e877', '9845e2f6561e3f33', '549a3a8168e9417a', '52f6608b85fbfa32',
    'c2c82600e6f6f451', 'dc7030a433ddda24', '5b2503692d42034b', 'c8e96b1e1357bb47',
    '059020e4b8545857', 'ac0d4da3af4a47eb', '24154cb03b59eb92', '0c8aa1c6c954b481',
    'e873f6f87dff7780', 'a4416e40a7c6b756', 'aec44160df5cd6a3', 'ba661f9ebff00c5c',
    'ed0820a4800919f7', '3dc40bb215acacb6', '063225b2278ccac5', '5ac11e1755e75ffc',
    'bb8201d23ca52ad6', 'b2c3f21196464fba',
)
_ESTIMATE_DIGESTS = (
    'ce0167b25dc06220', '90274f4344c11081', '06c991f66ccc4c03', '0abcb93447994131',
    'af0300f013f91255', '53aad3aa4d3c19ef', 'd5aa52ab51f913a1', '825a6de59160fd20',
    '6b6f6d5fafe66ab4', '4f1a760a00df793e', '5b71641a07b35cdd', '3cf8a76109ddb58c',
    '5beeb6d25a218cbf', '825cad65216cbe7a', 'f4b2d62a48ddc59b', 'b928d3121ef9f172',
    'a3cbe4fcb9922194', 'a6090fdd5cfa1b3d', '866cf71dd12a64ee', 'f622ce699331e87b',
    '8c065dc2d98fed54', '792dea06219be425', '0af42d8183b79d22', '0d9f26125e61fedd',
    '20b43ad48b50c7ac', '1727a17543d6bc3f', '3b91bf5c8677f040', 'a792a5fdf4c6fdf8',
    '7301ccdf80a98936', 'e572b57998539675', '982e7070c805c1e3', '1a1476bc2d040f9c',
    '3706f8a52f44e229', '50d4095f9ba77950', '823350df06622401', 'a8735a379fb96767',
    '542ca49c28614ce2', 'bd9c3d59be7edd4b', 'c39a84703adecaf4', '21ebbed502eb7302',
    '8c21d01fa37649d5', 'bcf5abf4907826bd', '088ae68fba0ac0fa', '6c336ca28caab433',
    '169936b7f401a522', 'a967ca99e4b23b48', '85e178190f594078', 'ff50fc4e0784540f',
    'a1a4d8adb35145f5', '72855e1893e1891c', 'eca6c34b7ff34375', '8cb3879c8be8ee87',
    '272e3e7aa2732839', 'd8831983cea415b0', '68928e65dd7686dc', '5ce3490fc973a493',
    '8ca927c67fd27201', '27458fd0f946de34', '514051a30b00398c', '0d15bfc5c7c8a20c',
    '660c885dd683d9ab', '229a612d25b0567c', 'ffe9352db6d0c183', '765c7b882bb5ff4b',
    '287bfc0686c4f6c8', '7d62764e351ca72a', '8a26efdeee1b7a59', 'cb67fbfee80e8c10',
    '945681303a262460', 'c0998f7aa4ec2934', '8497c2d0cde69d64', 'a6cd4849e43b9f12',
    'a7d6c1b802c96825', '9be3361e446331bf', '8142e2a5ffc31244', '6b27092de6e48511',
    '1b4d1c30dedaafee', '7ac1656d6d707699', 'b1e4cbeea8f27c13', '54e3d0af969903b1',
    'a1ea4da8d4fc1901', 'fdab61dedb7ff98d', '726a46db762bf00a', '14bff68f33976b19',
    '1e8c223c17a3db60', 'd13504944b9751bf', '880125353a8fdbe8', '8a2b3dd6fbf1efcb',
    'decbce339b3f5e88', '4d52a00e5a045c51', 'e8843fc8d9d58ceb', '466992c067190a8e',
    '7442d47e0a038141', '8223c295a55aab9a', '83b61f0da35dd415', '7097db906829c196',
    'd0476872dc08f07c', '6459f7083c9f542c', '58ccae045c52da6a', 'f01b3b6cf58b124c',
    'a911307d4d136331', '7da866ebb1060ee4', '1f14c8f3b90772fc', 'b05970669d03bfdb',
    'ae56a848dd0ad42e', 'd61ced3ab226c96f', 'bd27d51b03f74464', '6c19eda6f88371ae',
    'ae4c49486afccb94', '9e8b875179b05d0f', '11b8b25a585fa4f4', '3aecb50c03c159db',
    'b173ec7aa80f61c6', 'd0f3ea0fcc6bd00f', 'bba21ee3fc5d516e', 'efd5af68382e4b6a',
    '151048ac0a26b706', '791784a45352f56b', 'beaa95f3a428c451', '608e208659c5016f',
    'ddacecfdd7b611fe', 'a4608b9c38436c71', '927720483dc6b58d', '511b67e07c5a8738',
    '13ba950d5196fe4f', '9e76373411b4b313', '3b0c5175c27f9856', 'cfbeccbe5e70a925',
    'ce0167b25dc06220', '90274f4344c11081', '06c991f66ccc4c03', '0abcb93447994131',
    'af0300f013f91255', '53aad3aa4d3c19ef', 'd5aa52ab51f913a1', '825a6de59160fd20',
    '6b6f6d5fafe66ab4', '4f1a760a00df793e', '5b71641a07b35cdd', '3cf8a76109ddb58c',
    '5beeb6d25a218cbf', '825cad65216cbe7a', 'f4b2d62a48ddc59b', 'b928d3121ef9f172',
    'a3cbe4fcb9922194', 'a6090fdd5cfa1b3d', '866cf71dd12a64ee', 'f622ce699331e87b',
    '8c065dc2d98fed54', '792dea06219be425', '0af42d8183b79d22', '0d9f26125e61fedd',
    '20b43ad48b50c7ac', '1727a17543d6bc3f', '3b91bf5c8677f040', 'a792a5fdf4c6fdf8',
    '7301ccdf80a98936', 'e572b57998539675', '982e7070c805c1e3', '1a1476bc2d040f9c',
    '3706f8a52f44e229', '50d4095f9ba77950', '823350df06622401', 'a8735a379fb96767',
    '542ca49c28614ce2', 'bd9c3d59be7edd4b', 'c39a84703adecaf4', '21ebbed502eb7302',
    '8c21d01fa37649d5', 'bcf5abf4907826bd', '088ae68fba0ac0fa', '6c336ca28caab433',
    '169936b7f401a522', 'a967ca99e4b23b48', '85e178190f594078', 'ff50fc4e0784540f',
    'a1a4d8adb35145f5', '72855e1893e1891c', 'eca6c34b7ff34375', '8cb3879c8be8ee87',
    '272e3e7aa2732839', 'd8831983cea415b0', '68928e65dd7686dc', '5ce3490fc973a493',
    '8ca927c67fd27201', '27458fd0f946de34', '514051a30b00398c', '0d15bfc5c7c8a20c',
    '660c885dd683d9ab', '229a612d25b0567c', 'ffe9352db6d0c183', '765c7b882bb5ff4b',
    '287bfc0686c4f6c8', '7d62764e351ca72a', '8a26efdeee1b7a59', 'cb67fbfee80e8c10',
    '945681303a262460', 'c0998f7aa4ec2934', '8497c2d0cde69d64', 'a6cd4849e43b9f12',
    'a7d6c1b802c96825', '9be3361e446331bf', '8142e2a5ffc31244', '6b27092de6e48511',
    '1b4d1c30dedaafee', '7ac1656d6d707699', 'b1e4cbeea8f27c13', '54e3d0af969903b1',
    'a1ea4da8d4fc1901', 'fdab61dedb7ff98d', '726a46db762bf00a', '14bff68f33976b19',
    '1e8c223c17a3db60', 'd13504944b9751bf', '880125353a8fdbe8', '8a2b3dd6fbf1efcb',
    'decbce339b3f5e88', '4d52a00e5a045c51', 'e8843fc8d9d58ceb', '466992c067190a8e',
    '7442d47e0a038141', '8223c295a55aab9a', '83b61f0da35dd415', '7097db906829c196',
    'd0476872dc08f07c', '6459f7083c9f542c', '58ccae045c52da6a', 'f01b3b6cf58b124c',
    'a911307d4d136331', '7da866ebb1060ee4', '1f14c8f3b90772fc', 'b05970669d03bfdb',
    'ae56a848dd0ad42e', 'd61ced3ab226c96f', 'bd27d51b03f74464', '6c19eda6f88371ae',
    'ae4c49486afccb94', '9e8b875179b05d0f', '11b8b25a585fa4f4', '3aecb50c03c159db',
    'b173ec7aa80f61c6', 'd0f3ea0fcc6bd00f', 'bba21ee3fc5d516e', 'efd5af68382e4b6a',
    '151048ac0a26b706', '791784a45352f56b', 'beaa95f3a428c451', '608e208659c5016f',
    'ddacecfdd7b611fe', 'a4608b9c38436c71', '927720483dc6b58d', '511b67e07c5a8738',
    '13ba950d5196fe4f', '9e76373411b4b313', '3b0c5175c27f9856', 'cfbeccbe5e70a925',
    'ce0167b25dc06220', '90274f4344c11081', '06c991f66ccc4c03', '0abcb93447994131',
    'af0300f013f91255', '53aad3aa4d3c19ef', 'd5aa52ab51f913a1', '825a6de59160fd20',
    '6b6f6d5fafe66ab4', '4f1a760a00df793e', '5b71641a07b35cdd', '3cf8a76109ddb58c',
    '5beeb6d25a218cbf', '825cad65216cbe7a', 'f4b2d62a48ddc59b', 'b928d3121ef9f172',
    'a3cbe4fcb9922194', 'a6090fdd5cfa1b3d', '866cf71dd12a64ee', 'f622ce699331e87b',
    '8c065dc2d98fed54', '792dea06219be425', '0af42d8183b79d22', '0d9f26125e61fedd',
    '20b43ad48b50c7ac', '1727a17543d6bc3f', '3b91bf5c8677f040', 'a792a5fdf4c6fdf8',
    '7301ccdf80a98936', 'e572b57998539675', '982e7070c805c1e3', '1a1476bc2d040f9c',
    '3706f8a52f44e229', '50d4095f9ba77950', '823350df06622401', 'a8735a379fb96767',
    '542ca49c28614ce2', 'bd9c3d59be7edd4b', 'c39a84703adecaf4', '21ebbed502eb7302',
    '8c21d01fa37649d5', 'bcf5abf4907826bd', '088ae68fba0ac0fa', '6c336ca28caab433',
    '169936b7f401a522', 'a967ca99e4b23b48', '85e178190f594078', 'ff50fc4e0784540f',
    'a1a4d8adb35145f5', '72855e1893e1891c', 'eca6c34b7ff34375', '8cb3879c8be8ee87',
    '272e3e7aa2732839', 'd8831983cea415b0', '68928e65dd7686dc', '5ce3490fc973a493',
    '8ca927c67fd27201', '27458fd0f946de34', '514051a30b00398c', '0d15bfc5c7c8a20c',
    '660c885dd683d9ab', '229a612d25b0567c', 'ffe9352db6d0c183', '765c7b882bb5ff4b',
    '794d0f80a3055e1c', 'fb511aaf05dd9b42', '0ccd7b731185545f', 'e7b3fc21fd54a2c0',
    'e43c783f22a51a9b', 'aa4562f603f8da4d', '720348ae443443fb', 'a5537d74bbef74ce',
    'eefbb51741ab3c6f', '1f4ed16ddbd106c1', '2145f3d16ae91d53', 'ef7920d499c5dcbc',
    '2492b2b57eb4111b', 'c1de3096a51784f4', '5622e8f8d851b4b0', '0dcb54ba92a68f4c',
    '06ec1a77c403e7c6', 'f6cb0cbe38ecdc43', '3c607402552aa4a5', 'ae0c6c8a9935e5f7',
    '8edfa02191b4f543', '87d622dc09cea5be', 'cc0670f48d71a57a', '359ebb4ca30ae9ec',
    '4b52fa38dec61cb3', '3043a6343518ed55', '8e665d439e025496', '0560eb1c511952d9',
    'd9542bf772a5c229', 'b3ca1e0eaed01a02', 'c99dce64161907c0', 'adb644469cb92480',
    'd1859d5258fc7214', '1837594498a5c26a', '71ac148ff539b308', '79f3cd0202c04b86',
    '48561a82dd8161cf', '3abc5a8a2409b355', 'dee2fdbb9a6e58fb', '7ca8cdfd71668e63',
    'ffcdbfe744c12e62', '3496e0b7bb57dc9f', '44478192f149519c', '4e7526049e6a0e26',
    'e05268ad411fb95e', '411a7bb18adb7286', '669bc1eeb0fa878f', '036231778279b91b',
    '4e6e206f999af2af', '15dfe42c94628ff6', '148605c4f739fd7e', '09ab4ab63e7bc92f',
    '803b092e21ac837a', '3e1d71bb69cec8e9', '62b52d3a94e6685d', 'de0f47a781035e40',
    '6086fb868ae1393e', 'f8826d2e853fd1b5', '713e04ce25bde03a', '0f9a86931f7baf58',
    '71b12246a972a367', 'd1b64780776a3919', 'f7e16827de8229fd', '92425b54e75f4b5b',
    '05f6dcae914f80e7', '512a61846eb80f3c', '6f0b715fdf3d83bf', 'be6246af0feea4fe',
    'e32f69a414b86f45', '7b98ce9f1cd8375a', 'caea7c07af911040', '3dce492eba32de78',
    '2541eb1fefd8e557', '786462fa0476103d', '9495c01f5273dd40', '786462fa0476103d',
    '62e99512d3113b0c', '6086fb868ae1393e', '48b31057fdd34df6', '9850ce45eb173031',
    'c8ffcbbb843e2bd4', 'def7da684ec1713b', 'b2c3e94cab541963',
)


#: the over-cap refusal of a chi_dynamic plan (nmat 8, nvol 16, B 5, nd 2,
#: dense, cap 1000 bytes, cap_key = KEY) at a92b533d
_A92_DYNAMIC_REFUSAL = (
    "spin-orbital bond FLEX needs an estimated peak of 0.01 GiB in phase chi_dynamic (B = 5 "
    "channels, ND = B·(2·norb_phys)² = 20, nvol = 16, Nmat = 8; one store slot = "
    "16·Nmat·nvol·ND² bytes = 0.00 GiB) but the cap is 0.00 GiB ([mode.param] "
    "longitudinal_bond_memory_cap_gb). Memory grows as Nmat · nvol · B² · norb_phys⁴. Reduce "
    "Nmat or the k grid, or raise [mode.param] longitudinal_bond_memory_cap_gb (0 disables "
    "the check); pruning bond shells helps only for shells the per-shell table below shows "
    "with negligible coupling.\nrows (phase  name  bytes):\n  resident  gamma_b  102400\n"
    "  chi_dynamic  chibar  819200\n  chi_dynamic  bubble_pair_block  131072\n"
    "  chi_dynamic  dress_batch  5734400\nshells (shell  |dr|  max|v|):\n"
    "  (no off-site shell)")

#: the flex_scf rows, (nmat, nvol, B, nd, layout) -> rows (nb = min(16, nmat))
_FLEX_ROWS = {
    (16, 16, 5, 2, "dense"): (
        ("resident", "gamma_b", 102400), ("resident", "chibar", 1638400),
        ("bubble", "bubble_pair_block", 253952), ("dressing", "dress_batch", 11468800),
        ("dressing", "pole_slice", 204800), ("w_conversion", "w_batch_workspace", 6553600)),
    (16, 16, 5, 2, "sz_packed"): (
        ("resident", "gamma_b", 38400), ("resident", "chibar", 614400),
        ("bubble", "bubble_pair_block", 253952), ("dressing", "dress_batch", 4300800),
        ("dressing", "pole_slice", 76800), ("w_conversion", "w_batch_workspace", 2457600)),
    (8, 9, 5, 4, "dense"): (
        ("resident", "gamma_b", 921600), ("resident", "chibar", 7372800),
        ("bubble", "bubble_pair_block", 958464), ("dressing", "dress_batch", 51609600),
        ("dressing", "pole_slice", 1843200), ("w_conversion", "w_batch_workspace", 29491200)),
    (8, 9, 5, 4, "sz_packed"): (
        ("resident", "gamma_b", 345600), ("resident", "chibar", 2764800),
        ("bubble", "bubble_pair_block", 958464), ("dressing", "dress_batch", 19353600),
        ("dressing", "pole_slice", 691200), ("w_conversion", "w_batch_workspace", 11059200)),
}

#: flex_so_bond_rows(nmat=8, nvol=9, nd=4, depth=2, mixing="anderson", n_types=3)
_SO_RESIDENT = (("state", 20736), ("seed_envelope", 0), ("anderson_history", 147456),
                ("anderson_work", 73728), ("eigenpairs_hf", 11520), ("flex_arrays", 92160),
                ("hf_tables", 6912), ("collapses", 589824),
                ("solver_bookkeeping", 65536))
_SO_EXTRA = (("green_mu", "green_mu", 112896), ("density_hf", "density_hf", 117504),
             ("convergence", "convergence", 131328), ("mixing", "mixing", 73728),
             ("density_hf", "density_einsum_workspace", 294912),
             ("convergence", "green_digest", 36864), ("transport", "transport", 1253376),
             ("transport", "transport_einsum_workspace", 262144),
             ("post_scf", "post_scf", 313344),
             ("post_scf", "density_einsum_workspace", 294912))


def _phase_totals(rows, phases):
    resident = sum(b for ph, _, b in rows if ph == "resident")
    return {p: resident + sum(b for ph, _, b in rows if ph == p) for p in phases}


class TestRows(unittest.TestCase):

    def test_flex_rows_literal(self):
        for (nmat, nvol, B, nd, layout), want in _FLEX_ROWS.items():
            with self.subTest(nmat=nmat, nvol=nvol, B=B, nd=nd, layout=layout):
                plan = _flex_plan(nmat, nvol, B, nd, layout)
                self.assertEqual(plan.rows, want)
                self.assertEqual(plan.consumers, ("vertex", "flex_scf"))
                self.assertEqual((plan.nb, plan.dress_capacity), (min(16, nmat),) * 2)
                phases = ("vertex", "flex_scf", "bubble", "dressing", "w_conversion")
                totals = _phase_totals(want, phases)
                for p in phases:
                    self.assertEqual(plan.phase_total(p), totals[p])
                self.assertEqual(plan.peak_host, 1.25 * max(totals.values()))
                self.assertEqual(plan.peak_phase, "dressing")

    def test_aliases_and_slots(self):
        for w in (None, True):
            plan = _flex_plan(16, 16, 5, 2, w_in_place=w)
            self.assertEqual(plan.slots, ("chibar",))
            self.assertEqual(plan.aliases, (("W", "chibar"),))
            self.assertEqual(plan.rows, _FLEX_ROWS[(16, 16, 5, 2, "dense")])
        two = _flex_plan(16, 16, 5, 2, w_in_place=False)
        self.assertEqual(two.slots, ("chibar", "W"))
        self.assertEqual(two.aliases, ())
        want = list(_FLEX_ROWS[(16, 16, 5, 2, "dense")])
        want.insert(2, ("resident", "w_slot", 1638400))
        self.assertEqual(two.rows, tuple(want))
        dyn = plan_so_bond(nmat=16, nvol=16, B=5, nd=2, consumers=("chi_dynamic",),
                           caps=NOCAP)
        self.assertEqual((dyn.slots, dyn.aliases), (("chibar",), ()))
        with self.assertRaisesRegex(ValueError, r"^plan_so_bond: w_in_place must be None, "
                                                r"True or False, got 1$"):
            _flex_plan(16, 16, 5, 2, w_in_place=1)

    def test_nb_clamped_to_nmat(self):
        for cons in (("flex_scf",), ("chi_dynamic",)):
            plan = plan_so_bond(nmat=4, nvol=5, B=5, nd=2, consumers=cons, caps=NOCAP, nb=16)
            self.assertEqual((plan.nb, plan.dress_capacity), (4, 4))
            self.assertEqual(plan.rows, plan_so_bond(nmat=4, nvol=5, B=5, nd=2, consumers=cons,
                                                     caps=NOCAP, nb=4).rows)

    def test_w_in_place_default_unchanged(self):
        cases = existing_plan_cases()
        self.assertEqual(len(cases), len(_PLAN_DIGESTS))
        for (label, kw), want in zip(cases, _PLAN_DIGESTS):
            with self.subTest(label):
                self.assertEqual(_digest(_plan_outcome(kw)), want)
        self.assertEqual(sum(1 for _, kw in cases if _plan_outcome(kw)[0] == "ValueError"), 16)

    def test_w_in_place_refused_without_flex_scf(self):
        for w in (True, False):
            for cons in (("chi_dynamic",), ("vertex",), ("pairing_static",)):
                with self.subTest(w=w, cons=cons):
                    with self.assertRaises(ValueError) as cm:
                        plan_so_bond(nmat=8, nvol=4, B=5, nd=2, consumers=cons, caps=NOCAP,
                                     w_in_place=w)
                    self.assertEqual(str(cm.exception), "plan_so_bond: w_in_place applies to "
                                     "the consumer flex_scf only, got {!r}".format(w))
        plain = plan_so_bond(nmat=8, nvol=4, B=5, nd=2, consumers=("chi_dynamic",), caps=NOCAP)
        same = plan_so_bond(nmat=8, nvol=4, B=5, nd=2, consumers=("chi_dynamic",), caps=NOCAP,
                            w_in_place=None)
        self.assertEqual((plain.rows, plain.slots, plain.aliases),
                         (same.rows, same.slots, same.aliases))

    def test_flex_scf_with_chi_dynamic_refused(self):
        for cons in (("flex_scf", "chi_dynamic"), ("chi_dynamic", "flex_scf"),
                     ("flex_scf", "chi_dynamic", "pairing_static")):
            with self.subTest(cons=cons):
                with self.assertRaises(ValueError) as cm:
                    plan_so_bond(nmat=8, nvol=4, B=1, nd=2, consumers=cons, caps=NOCAP)
                self.assertEqual(str(cm.exception), "plan_so_bond: the consumer flex_scf "
                                 "cannot be combined with chi_dynamic (both own the store "
                                 "slot chibar)")

    def test_flex_scf_with_keep_chi_or_write_chibar_refused(self):
        for keep, write in ((True, False), (False, True), (True, True)):
            with self.subTest(keep_chi=keep, write_chibar=write):
                with self.assertRaises(ValueError) as cm:
                    _flex_plan(8, 4, 1, 2, keep_chi=keep, write_chibar=write)
                self.assertEqual(str(cm.exception), (
                    "plan_so_bond: the consumer flex_scf takes neither keep_chi nor "
                    "write_chibar (its store keeps no chi slot and writes no chibar archive), "
                    "got keep_chi = {!r}, write_chibar = {!r}").format(keep, write))
        # the defaults (and explicit False) are accepted
        self.assertEqual(_flex_plan(8, 4, 1, 2, keep_chi=False, write_chibar=False).rows,
                         _flex_plan(8, 4, 1, 2).rows)

    def test_resident_rows(self):
        base = _flex_plan(8, 9, 5, 4)
        extra = (("collapses", 1000), ("state", 7))
        plan = _flex_plan(8, 9, 5, 4, resident_rows=extra)
        self.assertEqual(plan.rows, base.rows + (("resident", "collapses", 1000),
                                                 ("resident", "state", 7)))
        for p in ("vertex", "flex_scf", "bubble", "dressing", "w_conversion"):
            self.assertEqual(plan.phase_total(p), base.phase_total(p) + 1007)
        # also with extra rows: a resident row lives through their phases too
        with_extra = _flex_plan(8, 9, 5, 4, resident_rows=extra,
                                extra_rows=(("transport", "transport", 5),))
        self.assertEqual(with_extra.phase_total("transport"),
                         sum(b for ph, _, b in base.rows if ph == "resident") + 1007 + 5)
        for name in ("gamma_b", "chibar", "dress_batch"):
            with self.subTest(name=name):
                with self.assertRaises(ValueError) as cm:
                    _flex_plan(8, 9, 5, 4, resident_rows=((name, 1),))
                self.assertEqual(str(cm.exception), "plan_so_bond: resident row {!r} "
                                 "duplicates a row of the plan".format(name))
        for rr, ex in (((("a", 1), ("a", 2)), ()), ((("transport", 1),),
                                                    (("transport", "transport", 5),))):
            with self.assertRaises(ValueError) as cm:
                _flex_plan(8, 9, 5, 4, resident_rows=rr, extra_rows=ex)
            self.assertIn("duplicates a row of the plan", str(cm.exception))
        for bad in ((("", 1),), (("x", -1),), (("x", 1.5),), (("x",),), ("x",)):
            with self.subTest(bad=bad):
                with self.assertRaises(ValueError) as cm:
                    _flex_plan(8, 9, 5, 4, resident_rows=bad)
                self.assertEqual(str(cm.exception), (
                    "plan_so_bond: resident_rows must hold (name, bytes) tuples with a "
                    "non-empty string name and bytes an integer >= 0, got {!r}").format(bad[0]))
        with self.assertRaises(ValueError) as cm:
            _flex_plan(8, 9, 5, 4, extra_rows=(("resident", "x", 1),))
        self.assertEqual(str(cm.exception), (
            "plan_so_bond: extra_rows must hold (phase, name, bytes) tuples with non-empty "
            "string phase and name, a phase other than 'resident' and bytes an integer >= 0, "
            "got ('resident', 'x', 1)"))

    def test_flex_so_bond_rows_literal(self):
        rows = _so_rows()
        res, ext = rows(nmat=8, nvol=9, nd=4, depth=2, mixing="anderson", n_types=3)
        self.assertEqual(res, _SO_RESIDENT)
        self.assertEqual(ext, _SO_EXTRA)
        self.assertTrue(all(type(b) is int for _, b in res))
        self.assertTrue(all(type(b) is int for _, _, b in ext))
        # the SCF rows are flex_scf_rows' with nmodes = nd
        per, phase = flex_bond.flex_scf_rows(nmat=8, nvol=9, nmodes=4, depth=2,
                                             mixing="anderson", n_types=3, split_seed=False)
        self.assertEqual(tuple(per.items()), res[:-2])
        self.assertEqual(tuple((k, k, v) for k, v in phase.items()), ext[:4])
        # the green digest's copy is capped at 64 MiB per copy
        _, big = rows(nmat=2048, nvol=4096, nd=4, depth=0, mixing="linear", n_types=1)
        self.assertIn(("convergence", "green_digest", 2 * 64 * 2 ** 20), big)
        # and the plan accepts them as they are
        plan = _flex_plan(8, 9, 5, 4, resident_rows=res, extra_rows=ext)
        self.assertEqual(plan.rows[len(_FLEX_ROWS[(8, 9, 5, 4, "dense")]):],
                         tuple(("resident",) + r for r in res) + ext)

    def test_split_seed_zero(self):
        rows = _so_rows()
        res, _ = rows(nmat=8, nvol=9, nd=4, depth=2, mixing="anderson", n_types=3)
        self.assertEqual(dict(res)["seed_envelope"], 0)
        res_t, _ = rows(nmat=8, nvol=9, nd=4, depth=2, mixing="anderson", n_types=3,
                        split_seed=True)
        self.assertEqual(dict(res_t)["seed_envelope"], 2 * 18432 + 2304)

    def test_non_so_estimate_bitwise(self):
        cases = estimate_cases()
        self.assertEqual(len(cases), len(_ESTIMATE_DIGESTS))
        for (label, kw), want in zip(cases, _ESTIMATE_DIGESTS):
            with self.subTest(label, kw=kw):
                self.assertEqual(_digest(_estimate_outcome(kw)), want)
        # the shared rows are flex_scf_rows' with nmodes = norb
        for depth, mixing, ss, nt in itertools.product((0, 3), ("linear", "anderson"),
                                                       (False, True), (0, 2)):
            est = flex_bond.estimate_bond_memory(
                nmat=8, nvol=16, norb=2, B=3, depth=depth, output_full=False, split_seed=ss,
                n_types=nt, freq_batch=None, cap_gb=8.0, mixing=mixing)
            per, phase = flex_bond.flex_scf_rows(nmat=8, nvol=16, nmodes=2, depth=depth,
                                                 mixing=mixing, n_types=nt, split_seed=ss)
            for k, v in per.items():
                self.assertEqual(est["persistent_rows"][k], v)
            for k, v in phase.items():
                self.assertEqual(est["phase_rows"][k], v)

    def test_w_batch_boundary(self):
        plan = _flex_plan(8, 9, 5, 4, nb=4)
        ND = 5 * 16
        plan.require("so_bond_w_batch", nvol=9, ND=ND, layout="dense", batch=4)
        with self.assertRaises(ValueError) as cm:
            plan.require("so_bond_w_batch", nvol=9, ND=ND, layout="dense", batch=5)
        self.assertEqual(str(cm.exception), "so_bond_w_batch: a batch of 5 slices exceeds the "
                         "plan's dress_capacity = 4 ({})".format(plan.reason))
        with self.assertRaisesRegex(ValueError, r"^so_bond_w_batch: the admission plan has "
                                                r"layout = dense but the call has layout = "
                                                r"sz_packed"):
            plan.require("so_bond_w_batch", nvol=9, ND=ND, layout="sz_packed", batch=1)
        # the bubble and the dressing of the flex plan pass their boundaries
        plan.require("so_bond_bubble", freq="all", nvol=9, B=5, nd=4, ND=ND, layout="dense",
                     nmat=8, placement="host", slots=plan.slots)
        plan.require("so_bond_dress", nvol=9, ND=ND, batch=4)
        dyn = plan_so_bond(nmat=8, nvol=9, B=5, nd=4, consumers=("chi_dynamic",), caps=NOCAP,
                           nb=4)
        with self.assertRaises(ValueError) as cm:
            dyn.require("so_bond_w_batch", nvol=9, ND=ND, layout="dense", batch=1)
        self.assertEqual(str(cm.exception), "so_bond_w_batch: the admission plan has no row "
                         "'w_batch_workspace' (closed consumers: vertex, chi_dynamic; {})"
                         .format(dyn.reason))

    def test_remedy_note(self):
        kw = dict(cap_key=KEY, envelope=so_bond.envelope_text("dense"))
        for keyed in (True, False):
            with self.subTest(keyed=keyed):
                extra = kw if keyed else {}
                with self.assertRaises(ValueError) as cm:
                    _flex_plan(8, 16, 5, 2, caps=SMALL, **extra)
                plain = str(cm.exception)
                with self.assertRaises(ValueError) as cm:
                    _flex_plan(8, 16, 5, 2, caps=SMALL, remedy_note="N.", **extra)
                noted = str(cm.exception)
                head = ("spin-orbital bond FLEX needs an estimated peak of 0.01 GiB in phase "
                        "dressing (B = 5 channels, ND = B·(2·norb_phys)² = 20, nvol = 16, "
                        "Nmat = 8; one store slot = 16·Nmat·nvol·ND² bytes = 0.00 GiB) but the "
                        "cap is 0.00 GiB ({}).".format(KEY if keyed else "caps.host_bytes"))
                first = noted.split("\n")[0]
                self.assertTrue(first.startswith(head + " N. Memory grows as "), first)
                self.assertTrue(plain.startswith(head + " Memory grows as "), plain)
                self.assertEqual(noted.replace(head + " N.", head, 1), plain)
        # without a note every existing text is byte-identical (a92b533d)
        with self.assertRaises(ValueError) as cm:
            plan_so_bond(nmat=8, nvol=16, B=5, nd=2, consumers=("chi_dynamic",), caps=SMALL,
                         cap_key=KEY)
        self.assertEqual(str(cm.exception), _A92_DYNAMIC_REFUSAL)
        with self.assertRaises(ValueError) as cm:
            plan_so_bond(nmat=8, nvol=16, B=5, nd=2, consumers=("chi_dynamic",), caps=SMALL,
                         cap_key=KEY, remedy_note=None)
        self.assertEqual(str(cm.exception), _A92_DYNAMIC_REFUSAL)
        # the solver word: Eliashberg for pairing_static (also on the static path)
        for cons in (("pairing_static",), ("chi_dynamic", "pairing_static")):
            with self.assertRaises(ValueError) as cm:
                plan_so_bond(nmat=8, nvol=16, B=5, nd=2, consumers=cons, caps=SMALL,
                             remedy_note="N.")
            self.assertTrue(str(cm.exception).startswith("spin-orbital bond Eliashberg needs"))
            self.assertIn("(caps.host_bytes). N. Memory grows as", str(cm.exception))
        rows = so_bond.static_so_bond_rows(
            nmat=8, nvol=16, B=5, nd=2, layout="dense", solver_mode="eigenvalue",
            num_eigenvalues=10, max_iter=1000, spatial_shape=(4, 4, 1))
        with self.assertRaises(ValueError) as cm:
            plan_so_bond(nmat=8, nvol=16, B=5, nd=2, consumers=("pairing_static",),
                         caps=SMALL, extra_rows=rows, remedy_note="N.", cap_key=KEY)
        msg = str(cm.exception)
        self.assertTrue(msg.startswith("spin-orbital bond Eliashberg needs"))
        self.assertIn("(" + KEY + "). N. Memory grows as Nmat · nvol · norb_phys⁴", msg)
        for bad in ("", 3, b"N."):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, r"^plan_so_bond: remedy_note must be "
                                                        r"None or a non-empty string, got "):
                    _flex_plan(8, 16, 5, 2, remedy_note=bad)


if __name__ == "__main__":
    unittest.main()
