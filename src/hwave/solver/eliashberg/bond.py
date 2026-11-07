"""The static bond-resolved pairing path of hwave_sc: topology, vertices, the bond kernel operator, the bond Green function and the diagnostics record.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 4); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import numpy as np
from scipy.sparse.linalg import LinearOperator

from hwave.solver.rpa import TAIL_ENDPOINT_CONVENTION
from hwave.solver import bond_channels
from hwave.solver import bubble
from hwave.solver import green as green_mod
from hwave.solver.eliashberg.config import _validate_bond_green_nmat
from hwave.solver.eliashberg.kernel import _green_sc_to_canonical, _green_canonical_to_sc, _BondGreen
from hwave.solver.eliashberg.eigen import _odd_harmonic_basis, _project_gap_parity
from hwave.solver.eliashberg.io import _refuse_so_green_file, _load_green_npz

logger = logging.getLogger("hwave_sc")


# Degeneracy tolerance of the opt-in bond diagnostics' eigenvalue clustering
# ([eliashberg] bond_diagnostics; spec S7.7 default).
_BOND_DIAGNOSTICS_DEG_TOL = 1.0e-3


# Number of (N_q, ND, ND)-shaped complex128 arrays that the bond-channel
# resource preflight (below) assumes are alive simultaneously. The high-water
# mark of ``_build_bond_operator`` is the moment a dressed channel is added
# to the pairing vertex: chi_bar, S, C, the pairing vertex Gamma, the dressed
# channel chi, the sandwich product V chi V and its coefficient-scaled copy
# (``PairVertexAccumulator``), plus one transient (the dressing denominator,
# the matrix-product intermediate V chi, or the kernel's transformed blocks).
# Measured at 7.1 arrays (``tests/test_sc_bond.py::TestPreflightOperatorBudget``).
# Named so the preflight's byte estimate and this list can never silently
# drift apart.
_BOND_N_Q_ARRAYS = 8


def _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat, tail_on):
    """Byte budget for the bond path, broken down by buffer family.

    Split out of :func:`_bond_resource_preflight` so the accounting can be
    asserted directly against a MEASURED peak (``tests/test_sc_bond.py``):
    a preflight that undercounts is worse than no preflight, because a run it
    waves through then OOMs anyway.

    Everything is complex128 (16 B) end to end on the bond path (the
    unified-bubble-kernel spec's "Memory contract": ``build_green`` produces
    complex128, the external-npz loader promotes to complex128 at load, so
    there is one itemsize and no input-vs-compute dtype split). ``unit =
    N_q * nmat * 16`` is one ``norb``-pair-free frequency-grid buffer, so
    ``S_in = norb**2 * unit`` is one FULL-size Green-function-shaped array
    (``(norb, norb, Nx, Ny, Nz, nmat)`` worth of complex128).

    ``q_bytes``
        ``_BOND_N_Q_ARRAYS`` q-resolved ``(N_q, ND, ND)`` arrays alive
        simultaneously in the operator assembly (see the list beside that
        constant). ``bond_bubble_static``'s output ``chi_bar`` is the first of
        these, so it is NOT counted again in ``bubble_bytes``.
    ``bubble_bytes``
        ``bond_bubble``'s simultaneous working set at its high-water mark, the
        Green carrier excluded (that is ``carrier_bytes``):
        ``BOND_BUBBLE_N4_BUFFERS`` ``norb**4``-sized plus
        ``BOND_BUBBLE_N2_BUFFERS`` ``norb**2``-sized frequency-grid buffers.
        The buffer-by-buffer list, and the ``del``/in-place discipline in
        ``bond_bubble`` that holds the real peak down to it, are documented
        beside those two constants in ``hwave.solver.bond_channels``; they are
        imported from there so the two sides cannot silently desync.
    ``carrier_bytes``
        The Green carrier's (:class:`_BondGreen`) contribution, replacing
        the pre-tail-correction ``green_bytes`` (one bare ``green_kw``
        array). Derived from TWO phases the carrier passes through (spec
        "Memory contract"), each expressed in units of ``S_in``:

        * carrier CONSTRUCTION (:func:`_build_bond_green`, transient):
          ``4 * S_in`` when ``tail_on`` (``build_green``'s canonical
          ``full_kw`` + canonical ``deflated_kw`` + ``green0_tail`` all
          alive at once, PLUS ``full_sc`` during the layout conversion,
          before the canonical ``full_kw`` is released) and ``2 * S_in``
          tail-off (``full``/``deflated`` alias one array; the sc-layout
          conversion adds the second). The external-npz branch peaks at
          ``2 * S_in`` (source + complex128-promoted copy) before the
          source is released -- covered by the same tail-off bound.

          Pass B review, round 1, measured ~5-7x ``S_in`` here on a SMALL
          fixture and (wrongly) concluded the ``4x`` bound undercounted;
          round 2 found the true cause and fixed it at the source instead
          of inflating this estimate: ``build_green``'s locals ``Vg``
          (full-size) and ``g_deflated`` were never released and stayed
          alive through the ``green0_tail``/``full_kw`` assembly (see
          ``green.build_green``'s ``del Vg, g_deflated`` -- load-bearing
          for this bound), and round 1's tiny fixture additionally paid a
          FIXED small-buffer overhead (``V_conj_t``/``VVt``/frequency-grid
          arrays, not ``S_in``-scaled) that is disproportionate relative to
          a small ``S_in`` but negligible at realistic problem sizes. With
          both effects accounted for, the measured peak is ``4 * S_in``
          plus a small FIXED residual (measured 888-952 bytes across
          several fixture sizes -- NOT ``S_in``-scaled, so it vanishes in
          relative terms as ``S_in`` grows), pinned by
          ``tests.test_sc_bond.TestPreflightMemoryBudget.
          test_bond_green_construction_peak_within_carrier_budget`` (which
          carries a small, honestly-documented absolute slack for exactly
          that residual in its own comparison -- not in this estimate).
        * carrier SETTLED / bubble phase: ``n_green * S_in`` with
          ``n_green = 3`` (``full_sc`` + canonical ``deflated_kw`` +
          ``green0_tail``) when ``tail_on`` else ``2``, PLUS one extra
          ``S_in``-sized endpoint-branch buffer while ``tail_on`` (the
          bubble's per-pair tail-endpoint correction) -- ``4 * S_in``
          tail-on, ``2 * S_in`` tail-off, both equal to their respective
          construction-phase bound.

        Both phases top out at ``4 * S_in`` (tail-on) / ``2 * S_in``
        (tail-off), so ONE closed-form bound covers the whole carrier
        lifetime: ``carrier_bytes = (4 if tail_on else 2) * norb**2 *
        unit``.
    ``chi_bar_bytes``
        one ``(N_q, ND, ND)`` array; reported separately (it is already inside
        ``q_bytes``) because it is the one q-resolved array ``bond_bubble``
        itself allocates, so ``bubble_bytes + chi_bar_bytes`` is the budget a
        measurement of ``bond_bubble`` alone must fit into.
    ``peak``
        ``q_bytes + bubble_bytes + carrier_bytes`` -- the SUM combine rule;
        what the cap is applied to.

    Parameters
    ----------
    tail_on : bool
        ``coeff_tail_applied != 0`` (NEVER the requested value -- an
        external Green function always resolves ``tail_on = False``
        regardless of what was requested, since ``coeff_tail_applied`` is
        forced to ``0.0`` on that branch).
    """
    nd = norb * norb
    B = int(bond_set.n_channels)
    ND = nd * B
    Nq = int(Nx) * int(Ny) * int(Nz)
    itemsize = 16  # complex128
    unit = Nq * int(nmat) * itemsize

    chi_bar_bytes = Nq * ND * ND * itemsize
    q_bytes = _BOND_N_Q_ARRAYS * chi_bar_bytes
    bubble_bytes = (bond_channels.BOND_BUBBLE_N4_BUFFERS * norb ** 4
                    + bond_channels.BOND_BUBBLE_N2_BUFFERS * norb ** 2) * unit
    carrier_bytes = (4 if tail_on else 2) * (norb ** 2) * unit
    return {"nd": nd, "B": B, "ND": ND, "Nq": Nq,
            "chi_bar_bytes": chi_bar_bytes,
            "q_bytes": q_bytes,
            "bubble_bytes": bubble_bytes,
            "carrier_bytes": carrier_bytes,
            "peak": q_bytes + bubble_bytes + carrier_bytes}


def _bond_resource_preflight(norb, bond_set, Nx, Ny, Nz, nmat, cap_gb, *,
                             tail_on):
    """Estimate the peak memory of the bond path and refuse to exceed the cap.

    Spec S3.2 "Resource guard": only here are ``nd``, ``N_q``, ``B`` and the
    dtype all known, so this is where the ``ND = nd*B`` blow-up is caught --
    never a silent runaway allocation. Called BEFORE the Green carrier is
    built (review fix I-2, preserved through the unified-bubble-kernel
    switch: previously the bond-channel dispatch branch built the Green
    function first, so a large-Nmat run could OOM before the cap was
    ever consulted); the estimate includes the carrier's own allocation so
    the preflight is not blind to it.

    ``tail_on`` : bool
        Whether the bond path's tail correction will be applied
        (``coeff_tail_applied != 0``); known at this call site even though
        the carrier itself is built afterwards, because an external
        ``bond_green`` always resolves ``coeff_tail_applied = 0`` (see
        :func:`_build_bond_green`) regardless of what was requested.

    See :func:`_bond_memory_estimate` for the buffer-by-buffer accounting.
    """
    est = _bond_memory_estimate(norb, bond_set, Nx, Ny, Nz, nmat, tail_on)
    B, ND, Nq = est["B"], est["ND"], est["Nq"]
    q_bytes = est["q_bytes"]
    bubble_bytes = est["bubble_bytes"]
    carrier_bytes = est["carrier_bytes"]
    peak = est["peak"]

    logger.info(
        "Bond-channel preflight: B = %d channels, ND = nd*B = %d, "
        "N_q = %d, tail_on = %s, estimated peak memory = %.3f GB "
        "(cap %.3f GB)", B, ND, Nq, tail_on, peak / 1.0e9, cap_gb)

    if peak > cap_gb * 1.0e9:
        raise ValueError(
            "[eliashberg] bond_channels: estimated peak memory {:.3f} GB "
            "exceeds bond_memory_cap_gb = {:.3f} GB. The cost is driven by "
            "B = {} bond channels (Delta r = {}) giving ND = nd*B = {} on "
            "N_q = {} q-points ({:.3f} GB of q-resolved ND x ND arrays) plus "
            "{:.3f} GB of bond-bubble work buffers plus {:.3f} GB for the "
            "Green carrier (tail_on={}). Use fewer bond channels (declare "
            "fewer off-site shells, or cut declared-zero shells with "
            "bond_max_shells), reduce the k-grid/Nmat, or raise "
            "bond_memory_cap_gb.".format(
                peak / 1.0e9, cap_gb, B, list(bond_set.delta_r), ND, Nq,
                q_bytes / 1.0e9, bubble_bytes / 1.0e9,
                carrier_bytes / 1.0e9, tail_on))
    return peak


def _bond_locality_split(interactions):
    """The pre-fold locality split of the static bond path's interactions.

    ``hwave_sc`` refuses sublattice folding (``SubShape`` must be
    ``[1, 1, 1]``), so the table read from the files is already the pre-fold
    one and :func:`hwave.solver.offsite.split_locality` runs on it directly.
    No ``orbital_site`` labels are passed: the bond channels treat the
    exchange of every bond exactly, so -- as on the FLEX bond gate -- only
    ``R = 0`` rows are on-site.
    """
    from types import SimpleNamespace
    from hwave.solver.offsite import split_locality
    return split_locality(SimpleNamespace(param_ham=interactions),
                          SimpleNamespace(subshape=(1, 1, 1)))


def _resolve_bond_topology(interactions, geom_info, norb, bond_max_shells):
    """Resolve the bond topology of the static bond path.

    The same resolver as the FLEX bond gate and the dynamic bond pairing
    (``bond_channels.resolve_bond_topology``) over every off-site type the
    bond channels carry (``bond_channels.longitudinal_bond_types``:
    CoulombInter, Hund and Ising, plus Exchange and PairHop when declared
    off-site). Channel 0 is the on-site ``R = 0`` point; the off-site shells
    follow by bond length (measured with the lattice vectors of the
    geometry file), then lexicographically. ``bond_max_shells`` is the
    resolver's ``max_shells``: it keeps that many off-site shells and
    refuses a cutoff that would drop a declared nonzero coefficient.

    Returns
    -------
    (topo, view, split)
        The :class:`~hwave.solver.bond_channels.BondTopology`, its
        :class:`~hwave.solver.bond_channels.BondSetView` and the locality
        split (:func:`_bond_locality_split`) the channel-0 block is built
        from.
    """
    split = _bond_locality_split(interactions)
    if "PairLift" in split.offsite_types:
        logger.info(
            "[eliashberg] bond_channels=true: off-site PairLift declarations "
            "carry no longitudinal (spin/charge) content and are not part of "
            "the bond channels.")
    topo = bond_channels.resolve_bond_topology(
        split.whole_tbl, geom_info["rvec"], norb, max_shells=bond_max_shells,
        active_types=bond_channels.longitudinal_bond_types(split.offsite_types))
    return topo, bond_channels.BondSetView(topo), split


def _build_bond_vertices(topo, split, norb, spatial_shape):
    """The enlarged spin/charge vertices ``S``, ``C`` on the bond basis.

    Built exactly as the FLEX bond gate builds them: the channel-0 block is
    the locality-split pair-space matrix
    (:func:`hwave.solver.offsite.sc_matrices_from_split`; the off-site
    Hartree ``V(q)`` of CoulombInter / Hund / Ising sits in its density
    slots, the Fock crossing of an ``R = 0`` row in its cross slots), and
    every off-site row of every type the topology carries enters through
    ``bond_channels.build_sc_bond_channel`` (the bond-diagonal Fock
    crossing; the PairHop bond entries).

    Returns
    -------
    S, C : ndarray, shape (nvol, ND, ND), complex
        ``ND = B * norb**2``, bond-major index ``m * norb**2 + a * norb + b``.
    """
    from hwave.solver.offsite import sc_matrices_from_split
    Nx, Ny, Nz = (int(x) for x in spatial_shape)
    nvol, nd = Nx * Ny * Nz, norb * norb
    S0, C0 = sc_matrices_from_split(
        split, bond_channels._LONGITUDINAL_ACTIVE_TYPES, norb, Nx, Ny, Nz)
    S0 = S0.reshape(nvol, nd, nd)
    C0 = C0.reshape(nvol, nd, nd)
    types = tuple(topo.coeffs)
    S = bond_channels.build_sc_bond_channel(
        topo, S0, "S", types=types, spatial_shape=(Nx, Ny, Nz))
    C = bond_channels.build_sc_bond_channel(
        topo, C0, "C", types=types, spatial_shape=(Nx, Ny, Nz))
    return S, C


def _bond_kernel_parts(S, C, dressed, green_kw, view, norb, spatial_shape,
                       beta, pairing_type, *, topo, g2_tail=False):
    """The static bond Eliashberg operator and its two parts, built from the
    bond pairing pieces of ``hwave.solver.eliashberg_bond``.

    Parameters
    ----------
    S, C : ndarray, shape (nvol, ND, ND) or (Nx, Ny, Nz, ND, ND)
        The enlarged bare spin/charge vertices (:func:`_build_bond_vertices`).
    dressed : iterable of (channel, chi)
        ``("spin", chi_s)`` and ``("charge", chi_c)``, each ``chi`` of the
        shape of ``S``: the dressed susceptibilities on the bond basis. Each
        is added to the pairing vertex as it arrives, so a generator keeps
        one channel resident at a time.
    green_kw : ndarray, shape (norb, norb, Nx, Ny, Nz, nmat)
        The full Green function (``sc.py`` layout) of the pair bubble.
    view : BondSetView
        The bond topology.
    pairing_type : {"singlet", "triplet"}
    topo : bond_channels.BondTopology or None
        The topology ``S`` and ``C`` were built on, for the singlet bare
        pairing of its off-site Exchange and PairHop rows
        (``eliashberg_bond.instantaneous_vertex``); ``None`` states that
        there are none.

    The fluctuation vertex ``Gamma = c_s S chi_s S + c_c C chi_c C``
    (``PairVertexAccumulator``) and the bare vertex (``instantaneous_vertex``:
    the channel-0 block, plus the singlet pairing of the off-site Exchange
    and PairHop rows of ``topo``) are applied by ``BondPairKernel`` on a
    one-point frequency axis, where its Matsubara transforms are the
    identity: one matvec is ``-(1/N) sum e^{-ik.R_b} Gamma_ab(k - k')
    e^{-ik'.R_a} X(k')`` with the static pair bubble ``X = G2 phi``,
    ``G2 = T sum_n G(k, iw_n) G(-k, -iw_n)``.

    The instantaneous part acts on the gap projected onto the parity sector
    of the pairing channel (``_project_gap_parity``): the bare Cooper vertex
    acts on the pair amplitude, which is even (singlet) or odd (triplet)
    under ``Delta_ab(k) -> Delta_ba(-k)``. On that sector the projection is
    the identity; on the opposite sector it keeps the bare vertex from
    acting on a pair amplitude the channel cannot carry. The fluctuation
    part is applied unprojected, as on the non-bond static path.

    Returns
    -------
    A, A_fl, A_pp : scipy.sparse.linalg.LinearOperator
        The full kernel ``K = K^fl + K^pp``, its fluctuation part and its
        instantaneous (bare) part, acting on a gap ``(norb, norb, Nx, Ny,
        Nz)`` raveled.
    vec_size : int
        ``norb**2 * Nx * Ny * Nz``.
    """
    from types import SimpleNamespace
    from hwave.solver.eliashberg import bond_dynamic as _eb
    from hwave.solver import flex_bond as _fb

    Nx, Ny, Nz = (int(x) for x in spatial_shape)
    shape = (Nx, Ny, Nz)
    nvol, nd = Nx * Ny * Nz, norb * norb
    S = np.asarray(S)
    C = np.asarray(C)
    ND = S.shape[-1]
    S = S.reshape(nvol, ND, ND)
    C = C.reshape(nvol, ND, ND)
    if pairing_type not in _eb.COEFF:
        raise ValueError(
            "unknown pairing_type '{}'. Use 'singlet' or 'triplet'.".format(
                pairing_type))

    with _fb.BondDeviceContext(np, S, C) as dev:
        acc = _eb.PairVertexAccumulator(
            dev, pairing_types=(pairing_type,), nb=1, nmat=1, nvol=nvol,
            nd=nd, spatial_shape=shape)
        for channel, chi in dressed:
            chi = np.asarray(chi).reshape(1, nvol, ND, ND)
            acc.add_channel(channel, _eb.ArrayBlockSource({"chi": chi}, nd),
                            "chi")
            del chi
        vertex = acc.finish()[pairing_type]
        del acc
    V_inst = _eb.instantaneous_vertex(S, C, nd, pairing_type, shape, topo=topo)

    # the static pair bubble on a one-point frequency axis
    G2 = bond_channels._g2_from_green(green_kw, beta, tail=g2_tail)[..., None]
    # The static path has its own resource preflight (bond_memory_cap_gb),
    # run before the Green function was built; the kernel's admission is
    # therefore not capped a second time.
    common = dict(xp=np, spatial_shape=shape, norb=norb, beta=beta, nfreq=1,
                  residency="host", host_cap=np.inf)
    kernel_fl = _eb.BondPairKernel(vertex, G2, view, **common)
    # The bare vertex lives on channel 0 alone, so its kernel needs only the
    # on-site channel and a zero fluctuation block.
    view0 = SimpleNamespace(n_channels=1, delta_r=((0, 0, 0),), reverse=(0,))
    zero = _eb.PairVertexUniform(
        source=_eb.ArrayBlockSource(
            {"zero": np.zeros((1, nvol, nd, nd), dtype=np.complex128)}, nd),
        slot="zero", B=1, nd=nd, nmat=1, nvol=nvol)
    kernel_pp = _eb.BondPairKernel(zero, G2, view0, V_inst=V_inst, **common)
    del G2

    gap_shape = (norb, norb, Nx, Ny, Nz)
    vec_size = nd * nvol

    def _fl(v):
        return kernel_fl.matvec(v)

    def _pp(v):
        phi = np.asarray(v).reshape(gap_shape)
        return kernel_pp.matvec(_project_gap_parity(phi, pairing_type).ravel())

    def _full(v):
        return _fl(v) + _pp(v)

    A, A_fl, A_pp = (LinearOperator((vec_size, vec_size), matvec=f,
                                    dtype=complex)
                     for f in (_full, _fl, _pp))
    return A, A_fl, A_pp, vec_size


def _build_bond_operator(topo, split, green_carrier, norb, kx_array, ky_array,
                         kz_array, beta, pairing_type, *, bond_max_shells,
                         bond_memory_cap_gb, precondition_opts=None,
                         green_source=None, g2_tail=False,
                         bond_coeff_tail_source="default"):
    """Build the static bond-resolved Eliashberg operator.

    The static path is the zero-bosonic-frequency member of the bond-resolved
    pairing family the FLEX bond gate and the dynamic bond pairing share, and
    it is built from the same modules:

    1. the bond bubble ``chibar(q)``: ``bubble.bond_bubble_static`` -- the
       ``Omega = 0`` slice of the per-pair pipeline that
       ``flex_bond.assemble_bubble`` streams over every bosonic frequency;
    2. the vertices ``S``, ``C`` (:func:`_build_bond_vertices`);
    3. the RPA dressing ``chi = [1 -/+ chibar W]^{-1} chibar``:
       ``bond_channels.dress_batch`` on a one-slice batch, one channel at a
       time (its conditioning guard names the channel and the q-point of a
       singular denominator);
    4. the fluctuation vertex ``Gamma = c_s S chi_s S + c_c C chi_c C``:
       ``eliashberg_bond.PairVertexAccumulator``;
    5. the bare (instantaneous) vertex: ``eliashberg_bond.instantaneous_
       vertex`` on the channel-0 block;
    6. the kernel: ``eliashberg_bond.BondPairKernel`` on a one-point
       frequency axis (:func:`_bond_kernel_parts`).

    ``topo``/``split`` come from :func:`_resolve_bond_topology` and must have
    been pre-flighted (``_bond_resource_preflight``) by the caller BEFORE
    ``green_carrier`` was built (review fix I-2: the preflight runs before the
    potentially large Green-function allocation).

    ``green_carrier`` : :class:`_BondGreen`
        The bubble reads ONLY ``deflated_kw``/``green0_tail``; the pair
        bubble ``G2`` and the pair weight read ``full_sc``.

    ``precondition_opts`` (from ``_read_bond_config``) is forwarded verbatim
    to ``bond_channels.check_hermitian_preconditions``.

    ``green_source`` is provenance only: the path the Green function was
    loaded from (``[eliashberg] bond_green``), or ``None`` when the caller
    built the bare Green function from the transfer Hamiltonian.

    ``bond_max_shells`` and ``bond_memory_cap_gb`` are likewise PROVENANCE
    ONLY here: the cutoff and the preflight already happened in the caller.

    Returns
    -------
    (A, vec_size) : tuple
        The ``scipy.sparse.linalg.LinearOperator`` kernel and its vector size,
        in the same convention as ``_make_kernel_operator``.
    provenance : dict
        Additive output-provenance record (spec S5).
    attribution : dict
        Context for the ``lambda = lambda^pp + lambda^fl`` decomposition: the
        two kernel parts (``pp``/``fl``), the full kernel and the pair weight
        (``weight``), consumed by ``bond_channels.attribute_lambda``.
    """
    Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
    shape = (Nx, Ny, Nz)
    view = bond_channels.BondSetView(topo)
    green_kw = green_carrier.full_sc

    logger.info("Computing the bond-resolved bubble chi-bar...")
    chi_bar = bubble.bond_bubble_static(
        green_carrier.deflated_kw, green_carrier.green0_tail, beta, view,
        spatial_shape=shape)

    S, C = _build_bond_vertices(topo, split, norb, shape)

    def _dressed(chi_bar, S, C):
        # one dressed channel at a time, on a one-slice batch at the zero
        # bosonic frequency (l0 = nmat = 0, so a refusal names index 0)
        for channel, W in (("spin", S), ("charge", C)):
            chi, _ = bond_channels.dress_batch(
                chi_bar[None], W, channel, l0=0, nmat=0, spatial_shape=shape)
            yield channel, chi[0]
            del chi

    logger.info("Building the bond-resolved Eliashberg operator "
                "(pairing_type=%s)...", pairing_type)
    A, A_fl, A_pp, vec_size = _bond_kernel_parts(
        S, C, _dressed(chi_bar, S, C), green_kw, view, norb, shape, beta,
        pairing_type, topo=topo, g2_tail=g2_tail)
    del chi_bar, S, C

    # v1 runtime preconditions of the Hermitian static path (spec S4.5): the
    # pair weight w = GG must be real (Hermitian) and >= 0 and the symmetrized
    # kernel K~ = -sqrt(w) Gamma sqrt(w) Hermitian, else the reported lambda is
    # not a real eigenvalue/Rayleigh quotient. Violations RAISE -- there is no
    # non-Hermitian biorthogonal fallback in v1.
    weight = bond_channels.pair_weight(green_kw, beta, g2_tail=g2_tail)
    diagnostics = bond_channels.check_hermitian_preconditions(
        A, weight, label="eliashberg bond_channels",
        parts={"pp": A_pp, "fl": A_fl},
        **(precondition_opts or {}))
    logger.info(
        "Bond kernel Hermiticity preconditions OK (%s check): "
        "||K~ - K~^dag||_F = %.3e (relative %.3e; the same residual on "
        "M = W K is %.3e), min eig w = %.3e",
        diagnostics["method"],
        diagnostics["kernel_hermiticity_residual_ktilde"],
        diagnostics["kernel_hermiticity_relative_ktilde"],
        diagnostics["kernel_hermiticity_residual"],
        diagnostics["weight_min_eigenvalue"])

    attribution = {
        "full": A,
        "pp": A_pp,
        "fl": A_fl,
        "weight": weight,
        "diagnostics": diagnostics,
    }

    provenance = {
        "bond_channels": True,
        "bond_n_channels": int(view.n_channels),
        "bond_delta_r": [list(dr) for dr in view.delta_r],
        "bond_max_shells": ("all" if bond_max_shells is None
                            else int(bond_max_shells)),
        "bond_memory_cap_gb": float(bond_memory_cap_gb),
        "g2_tail": bool(g2_tail),
        "bond_coeff_tail_requested": float(green_carrier.coeff_tail_requested),
        "bond_coeff_tail_applied": float(green_carrier.coeff_tail_applied),
        "bond_coeff_tail_source": bond_coeff_tail_source,
        "bond_endpoint_convention": TAIL_ENDPOINT_CONVENTION,
        "bond_green_source": green_carrier.source,
        "external_full_green_assumed": (green_carrier.source == "external_npz"),
        "approximation": (
            "static RPA-ladder bond dressing on "
            + ("an EXTERNAL Green function supplied via [eliashberg] "
               "bond_green (its own approximation level -- e.g. FLEX-dressed "
               "and self-consistent -- is that of the file, not of this run); "
               "the bond chi-bar is still built here from that green, so "
               "this is NOT a conserving FLEX result"
               if green_source else
               "the BARE Green function built from the transfer Hamiltonian "
               "(supply [eliashberg] bond_green to feed an external/FLEX "
               "green instead); NOT a conserving FLEX result")
            + " -- absolute lambda is not comparable to FLEX/dynamic values"),
    }
    if green_source:
        provenance["bond_green"] = str(green_source)
    if view.n_channels == 1:
        provenance["collapsed_to_pure_hubbard"] = True
    provenance["kernel_hermiticity_residual"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_residual"])
    provenance["kernel_hermiticity_residual_ktilde"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_residual_ktilde"])
    provenance["kernel_hermiticity_relative_ktilde"] = "{:.3e}".format(
        diagnostics["kernel_hermiticity_relative_ktilde"])
    provenance["kernel_hermiticity_method"] = diagnostics["method"]
    provenance["pair_weight_min_eigenvalue"] = "{:.6e}".format(
        diagnostics["weight_min_eigenvalue"])

    return (A, vec_size), provenance, attribution


def _build_bond_green(eigenvalues, eigenvectors, mu, beta, nmat,
                      coeff_tail_requested, bond_green_path,
                      spin_orbital_mode=False):
    """Build the bond path's :class:`_BondGreen` carrier.

    Internal branch (``bond_green_path is None``): flattens sc.py's
    ``(Nx, Ny, Nz, p)`` / ``(Nx, Ny, Nz, p, p)`` eigen-decomposition to the
    kernel's flattened-only ``(nvol, p)`` / ``(nvol, p, p)`` contract (C
    order -- the SAME flattening :func:`_green_sc_to_canonical` uses, so
    k-point ordering is consistent by construction), calls
    ``hwave.solver.green.build_green`` with the CALLER's already-computed
    ``mu`` (no mu search is ever re-run here), and converts the canonical
    ``full_kw`` to ``full_sc`` ONCE -- the canonical array is released
    (``del``) immediately after, so the two full-size layouts coexist only
    for the duration of that one conversion. ``coeff_tail_applied ==
    coeff_tail_requested``, ``source = "internal"``.

    External branch (``bond_green_path`` set): loads the npz via
    :func:`_load_green_npz` (promoted to complex128 at load if needed),
    ``full_sc = loaded``, ``deflated_kw = _green_sc_to_canonical(full_sc)``,
    ``green0_tail = None``, ``coeff_tail_applied = 0.0``, ``source =
    "external_npz"``. A nonzero ``coeff_tail_requested`` logs a WARNING
    (the tail correction is not applicable to an externally supplied Green
    function, whose high-frequency asymptote is not known here) rather
    than raising -- recorded as requested-but-not-applied, not an error.

    Parameters
    ----------
    eigenvalues : ndarray, shape (Nx, Ny, Nz, p)
        Only used for its VALUES on the internal branch; on the external
        branch only its SHAPE is read (to validate the loaded file's
        ``norb``/``Nx``/``Ny``/``Nz`` against the model).
    eigenvectors : ndarray, shape (Nx, Ny, Nz, p, p)
        As above.
    mu : float
        The caller's already-determined chemical potential.
    beta : float
        Inverse temperature.
    nmat : int
        Target Matsubara-frequency count for an internal build; unused on
        the external branch (the file's own frequency count wins -- the
        caller is responsible for reconciling the two, exactly as the
        pre-carrier code did).
    coeff_tail_requested : float
        The resolved ``[eliashberg] bond_coeff_tail`` value.
    bond_green_path : str or None
        ``[eliashberg] bond_green``; ``None`` selects the internal build.
    spin_orbital_mode : bool
        ``False`` (the non-spin-orbital bond path): an external file is
        checked for spin-orbital markers first and refused (the #271
        preflight does not run on the bond path). ``True`` skips the check.

    Returns
    -------
    _BondGreen
    """
    if bond_green_path is not None:
        if not spin_orbital_mode:
            _refuse_so_green_file(bond_green_path)
        Nx, Ny, Nz, norb = eigenvalues.shape
        green_sc, nfreq = _load_green_npz(
            bond_green_path, norb, Nx, Ny, Nz, label="bond_green")
        _validate_bond_green_nmat(nfreq, bond_green_path)
        if green_sc.dtype != np.complex128:
            green_sc = green_sc.astype(np.complex128)
        if coeff_tail_requested != 0.0:
            logger.warning(
                "[eliashberg] bond_coeff_tail = %s requested but "
                "bond_green = '%s' supplies an EXTERNAL Green function: "
                "the tail correction needs the analytic high-frequency "
                "asymptote of an internally built bare Green function, "
                "which is not available for a file. NOT applying it here "
                "(recorded as requested but not applied).",
                coeff_tail_requested, bond_green_path)
        return _BondGreen(
            full_sc=green_sc,
            deflated_kw=_green_sc_to_canonical(green_sc),
            green0_tail=None,
            coeff_tail_requested=float(coeff_tail_requested),
            coeff_tail_applied=0.0,
            source="external_npz")

    Nx, Ny, Nz, p = eigenvalues.shape
    nvol = Nx * Ny * Nz
    ev_flat = eigenvalues.reshape(nvol, p)
    evec_flat = eigenvectors.reshape(nvol, p, p)
    full_kw, deflated_kw, green0_tail = green_mod.build_green(
        ev_flat, evec_flat, mu, beta, nmat, coeff_tail_requested)
    full_sc = _green_canonical_to_sc(full_kw, (Nx, Ny, Nz))
    del full_kw
    return _BondGreen(
        full_sc=full_sc,
        deflated_kw=deflated_kw,
        green0_tail=green0_tail,
        coeff_tail_requested=float(coeff_tail_requested),
        coeff_tail_applied=float(coeff_tail_requested),
        source="internal")


def _bond_diagnostics_record(provenance, gap, eigenvalues, weight, norb,
                             kx_array, ky_array, kz_array,
                             deg_tol=_BOND_DIAGNOSTICS_DEG_TOL,
                             sector_eigenvalues=None, leader_available=True):
    """Opt-in character analysis of the leading bond state (``[eliashberg]
    bond_diagnostics``; spec S7.7).

    Writes into the ADDITIVE provenance record (``#`` comment lines of
    ``eigenvalue.dat``, see ``_save_results``) two things a single solve can
    honestly report about its leading state:

    1. ``bond_diagnostics_harmonics`` -- the odd-harmonic decomposition of the
       returned gap on ``_odd_harmonic_basis`` in the ``sqrt(GG)`` metric
       (``bond_channels.harmonic_decomposition``). Each number is the fraction
       of that harmonic captured by the state, in ``[0, 1]``; it is the
       V-sweep "character" quantity of spec S7.7 evaluated at ONE point.
    2. ``bond_diagnostics_eigenvalue_clusters`` -- the near-degenerate grouping
       of the computed eigenvalues (``bond_channels.cluster_eigenvalues``), so
       a user can see whether the leading state is a doublet (the odd E
       doublet of the square lattice) or a singlet branch. Only present when a
       spectrum was actually computed (``solver_mode`` includes
       ``"eigenvalue"``); the iteration path returns one vector and no
       spectrum, and a one-eigenvalue "cluster list" would be a fake.

    The ``lambda^pp`` / ``lambda^fl`` attribution sits in the same block and is
    written by the caller regardless of this flag.

    **Eigenvector orientation.** ``eigenvalues``/``gap`` come from
    ``_solve_eigenvalue``, which already returns eigenvectors as one gap array
    PER LEADING INDEX (it transposes ``scipy.sparse.linalg.eigs``'s COLUMNS);
    the row stack this function hands to ``harmonic_decomposition`` therefore
    holds one eigenvector per row, which is the convention the
    ``bond_channels`` subspace helpers (and ``track_subspace``) consume.

    Multi-point tracking across a ``V`` sweep -- following ONE invariant
    subspace through the spectrum rather than characterizing one solve -- is
    ``bond_channels.track_subspace``; a single ``calc_eliashberg`` run has only
    one point, so that helper is the documented API for sweep drivers instead.
    """
    provenance["bond_diagnostics"] = True

    basis = _odd_harmonic_basis(norb, kx_array, ky_array, kz_array)
    if basis and gap is not None:
        vec = np.asarray(gap).ravel()[np.newaxis, :]
        harmonics = bond_channels.harmonic_decomposition(vec, basis, weight)
        provenance["bond_diagnostics_harmonics"] = ", ".join(
            "{}={:.6f}".format(name, harmonics[name])
            for name in sorted(harmonics))
        provenance["bond_diagnostics_harmonics_note"] = (
            "odd-harmonic decomposition of the LEADING (returned) gap on the "
            "basis {sin kx, sin ky, sin kx (cos kx - cos ky), "
            "sin ky (cos ky - cos kx)} in the sqrt(GG) pair-weight metric; "
            "each value is the fraction of that harmonic captured by the "
            "state, in [0, 1] (spec S7.7). p_x/p_y/f_x/f_y span the SAME odd "
            "2D representation on a square lattice, so this is a character "
            "report, not a point-group separation")
        logger.info("bond diagnostics -- odd-harmonic content of the leading "
                    "state: %s", provenance["bond_diagnostics_harmonics"])
    elif not basis:
        provenance["bond_diagnostics_harmonics"] = (
            "n/a (every odd harmonic vanishes on this k-grid)")

    if not leader_available:
        # no reported leader (NaN): there is no state to describe
        provenance["bond_diagnostics_leading_cluster"] = "n/a"
        provenance["bond_diagnostics_leading_cluster_size"] = "n/a"
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "n/a (no leading eigenvalue was obtained)")
    elif sector_eigenvalues is not None and len(sector_eigenvalues) > 0:
        # the reported leader came from the channel's parity-sector solve:
        # its degeneracy is read off THAT spectrum (index 0 is the leader),
        # not off the global listing, whose rows may describe other states
        clusters = bond_channels.cluster_eigenvalues(
            np.asarray(sector_eigenvalues), deg_tol=deg_tol)
        leading = next(c for c in clusters if 0 in c)
        if eigenvalues is not None and len(eigenvalues) > 0:
            # the table's own groups, for reference (indices into its rows;
            # row 0 is the leader, whose state the global rows may list
            # again, so this grouping is not a degeneracy count)
            provenance["bond_diagnostics_eigenvalue_clusters"] = str(
                bond_channels.cluster_eigenvalues(np.asarray(eigenvalues),
                                                  deg_tol=deg_tol))
        provenance["bond_diagnostics_sector_eigenvalue_clusters"] = str(clusters)
        provenance["bond_diagnostics_leading_cluster"] = str(leading)
        provenance["bond_diagnostics_leading_cluster_size"] = len(leading)
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "channel parity-sector spectrum")
        provenance["bond_diagnostics_deg_tol"] = "{:.3e}".format(deg_tol)
        provenance["bond_diagnostics_clusters_note"] = (
            "near-degenerate groups (single-linkage on |lambda_i - lambda_j| "
            "<= deg_tol) of the eigenvalues of the requested channel's "
            "parity-sector solve, which produced the reported leader (its "
            "index 0); bond_diagnostics_leading_cluster is the group "
            "containing it, i.e. the degeneracy of the state this run "
            "returns; size 2 is the odd E doublet")
        logger.info("bond diagnostics -- sector eigenvalue clusters (deg_tol "
                    "%.1e): %s; the returned state's cluster is %s (size %d)",
                    deg_tol, clusters, leading, len(leading))
    elif eigenvalues is not None and len(eigenvalues) > 0:
        clusters = bond_channels.cluster_eigenvalues(np.asarray(eigenvalues),
                                                     deg_tol=deg_tol)
        # The cluster of the RETURNED eigenpair, i.e. of row 0 -- NOT
        # clusters[0]. cluster_eigenvalues orders by descending Re lambda over
        # every computed eigenpair, while row 0 is the eigenpair the run
        # actually reports (the channel's leader is row 0 of the table); an
        # opposite-parity mode with a larger Re lambda sits below it but is
        # still clustered here. Reporting clusters[0]'s size as "the leading
        # state's degeneracy" would then describe a different state.
        leading = next(c for c in clusters if 0 in c)
        provenance["bond_diagnostics_eigenvalue_clusters"] = str(clusters)
        provenance["bond_diagnostics_leading_cluster"] = str(leading)
        provenance["bond_diagnostics_leading_cluster_size"] = len(leading)
        provenance["bond_diagnostics_leading_cluster_source"] = (
            "eigen-analysis table")
        provenance["bond_diagnostics_deg_tol"] = "{:.3e}".format(deg_tol)
        provenance["bond_diagnostics_clusters_note"] = (
            "near-degenerate groups of the COMPUTED eigenvalues (indices into "
            "the rows below, single-linkage on |lambda_i - lambda_j| <= "
            "deg_tol, ordered by descending Re lambda over ALL computed "
            "eigenpairs -- including any opposite-parity mode below row 0). "
            "bond_diagnostics_leading_cluster is the group containing ROW 0, "
            "i.e. the degeneracy of the state this run actually returns; size "
            "2 is the odd E doublet")
        logger.info("bond diagnostics -- eigenvalue clusters (deg_tol %.1e): "
                    "%s; the returned state's cluster is %s (size %d)",
                    deg_tol, clusters, leading, len(leading))
    return provenance
