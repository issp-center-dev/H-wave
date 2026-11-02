"""The pairing vertex of the spin-orbital linearized Eliashberg equation (#83).

In spin-orbital mode (``[mode] enable_spin_orbital``) the gap is one
``2N x 2N`` matrix ``Delta_AB(k)`` over spin-orbital indices ``A = (s, a)``
(spin-block order ``s * N + a``), not a singlet and a triplet channel. Its
pairing vertex is built here from the two objects the spin-orbital RPA and
FLEX already use:

* ``Gamma`` -- the antisymmetrized spinful vertex
  (:func:`hwave.solver.rpa.spinful_bubble_vertex`: the direct tensor plus
  the on-site exchange crossing, in the bubble's pair convention), and
* ``chi = [1 + chi0 Gamma]^-1 chi0`` -- the spinful susceptibility of the
  same vertex (static here), in the bubble's ``[a, c, b, d]`` layout.

Derivation
----------
``Gamma[(x, y), (z, w)]`` multiplies ``(c+_w c_z)(c+_x c_y)``; it is
antisymmetric under the exchange crossing (``Gamma[(x,y),(z,w)] =
-Gamma[(x,z),(y,w)]``, the exchange image built in) and symmetric under the
exchange of the two particles (``Gamma[(x,y),(z,w)] = Gamma[(w,z),(y,x)]``).
Normal ordering ``c+_w c_z c+_x c_y = c+_w c+_x c_y c_z + (one-body)``
shows that, read in the pair slots, it is the antisymmetrized
particle-particle vertex ``A`` of ``H = 1/4 sum A_pq,rs c+_p c+_q c_s c_r``.
The first-order gap equation of an antisymmetrized vertex carries a 1/2
(``Delta = -1/2 A F``), so with ``Delta_ad = -sum_bc V_abcd F_bc`` (the
contraction of ``sc._make_kernel_operator``)

    V_bare[a, b, c, d] = 1/2 Gamma[(a, b), (c, d)],

written as the particle-exchange symmetric part
``1/4 (Gamma[(a,b),(c,d)] + Gamma[(d,c),(b,a)])``, equal to it for every
on-site type.

The fluctuation part is the particle-hole reducible vertex ``-Gamma chi
Gamma`` (the RPA screening of ``Gamma``; it carries the charge, the
longitudinal and the transverse spin channels and their mixing) at the
momentum transfer ``q = k - k'`` of the direct particle-hole channel:

    V_fluct[a, b, c, d](q) = -(Gamma chi(q) Gamma)[(a, b), (c, d)].

Its crossed image sits at ``k + k'``; on a gap of fermionic antisymmetry
``Delta_AB(k) = -Delta_BA(-k)`` it contributes exactly what the direct image
does, so the direct image enters with weight one and no 1/2 (the bare term,
which already contains both images, with 1/2). Reading the direct image off
the pair slots gives ``-(Gamma chi Gamma)[(d, c), (b, a)]`` at ``-q``; the
form above is the same tensor by the particle-exchange symmetry
``chi_{(x,y),(z,w)}(q) = chi_{(w,z),(y,x)}(-q)`` of the static
susceptibility (and of ``Gamma``).

With spin-independent hopping and a spin-rotation invariant interaction
this reproduces the non-spin-orbital singlet ``3/2 S chi_s S - 1/2 C chi_c C
+ 1/2 (S + C)`` and triplet ``-1/2 S chi_s S - 1/2 C chi_c C + 1/2 (C - S)``
(cross slots dropped, #262) kernels on ``phi (x) i sigma_y`` and
``phi (x) sigma_j i sigma_y``; the bare part equals an independent
first-order spin-orbital construction for every on-site type with
spin-orbit coupling (``tests/test_so_pairing_vertex.py``).

Only on-site interactions are representable: an off-site term makes
``Gamma`` depend on ``q``, and its exchange crossing then depends on
``k + k'`` (the same restriction as the spin-orbital RPA and FLEX). Such a
``Gamma`` is rejected.
"""
import numpy as np

__all__ = ["so_pairing_vertex", "so_pairing_vertex_parts",
           "so_pairing_vertex_dynamic", "so_initial_gap", "so_frequency",
           "check_so_config", "check_so_supported", "so_spin_components",
           "so_spin_decomposition", "SPIN_COMPONENTS", "SZ_COMPONENTS",
           "so_sz_decomposition", "sz_off_sector_index",
           "sz_conservation_residual", "sz_report"]

#: relative tolerance of the on-site (q-independent Gamma) check
_ONSITE_RTOL = 1.0e-10


def _validate(gam, chi, nd, spatial_shape):
    try:
        shape = tuple(int(s) for s in spatial_shape)
    except TypeError:
        raise ValueError("spatial_shape must be a 3-tuple of positive integers, got "
                         "{!r}".format(spatial_shape))
    if len(shape) != 3 or any(s <= 0 for s in shape):
        raise ValueError("spatial_shape must be a 3-tuple of positive integers, got "
                         "{!r}".format(spatial_shape))
    nvol = shape[0] * shape[1] * shape[2]
    nd = int(nd)
    want = (nvol, nd * nd, nd * nd)
    gam = np.asarray(gam)
    chi = np.asarray(chi)
    if gam.shape != want:
        raise ValueError("the spinful vertex Gamma must have shape {} (nvol, nd^2, nd^2) "
                         "for nd = {} and spatial_shape = {}, got {}".format(
                             want, nd, shape, gam.shape))
    if chi.shape != want:
        raise ValueError("the spinful susceptibility chi must have shape {} (nvol, nd^2, "
                         "nd^2), got {}".format(want, chi.shape))
    scale = float(np.abs(gam).max()) if gam.size else 0.0
    if scale > 0.0 and float(np.abs(gam - gam[:1]).max()) > _ONSITE_RTOL * scale:
        raise ValueError(
            "the spinful vertex Gamma depends on q: the spin-orbital pairing vertex "
            "supports on-site interactions only (the exchange crossing of an off-site "
            "term depends on k + k')")
    return gam, chi, nd, shape


def so_pairing_vertex_parts(gam, chi, nd, spatial_shape):
    """The bare and fluctuation parts of the spin-orbital pairing vertex.

    Parameters
    ----------
    gam : ndarray, (nvol, nd^2, nd^2)
        The antisymmetrized spinful vertex in the bubble's pair convention
        (:func:`hwave.solver.rpa.spinful_bubble_vertex`); on-site only, so
        independent of q.
    chi : ndarray, (nvol, nd^2, nd^2)
        The static spinful susceptibility ``[1 + chi0 Gamma]^-1 chi0`` in the
        bubble's ``[a, c, b, d]`` layout (zeros give the bare vertex alone).
    nd : int
        Number of spin-orbitals (spin-block order ``spin * norb + orb``).
    spatial_shape : (Nx, Ny, Nz)
        The k grid; ``nvol = Nx * Ny * Nz`` in C order.

    Returns
    -------
    (V_bare, V_fluct) : ndarrays, (nd, nd, nd, nd, Nx, Ny, Nz)
        In the layout of the general-mode ``Vs_q`` of
        ``sc._make_kernel_operator`` (``Delta_ad = -sum_bc V_abcd F_bc``):
        ``V_bare = 1/2 Gamma``, ``V_fluct = -Gamma chi Gamma`` (see the module
        docstring).
    """
    gam, chi, nd, shape = _validate(gam, chi, nd, spatial_shape)
    nvol = gam.shape[0]
    g4 = gam.reshape(nvol, nd, nd, nd, nd)
    bare = 0.25 * (g4 + g4.transpose(0, 4, 3, 2, 1))
    fluct = -(gam @ chi @ gam).reshape(nvol, nd, nd, nd, nd)

    def _layout(v):
        return np.ascontiguousarray(np.moveaxis(v, 0, -1)).reshape(
            (nd,) * 4 + shape).astype(complex, copy=False)

    return _layout(bare), _layout(fluct)


def so_pairing_vertex(gam, chi, nd, spatial_shape):
    """The spin-orbital pairing vertex ``V_bare + V_fluct``
    (:func:`so_pairing_vertex_parts`), shape ``(nd, nd, nd, nd, Nx, Ny, Nz)``."""
    bare, fluct = so_pairing_vertex_parts(gam, chi, nd, spatial_shape)
    return bare + fluct


def so_pairing_vertex_dynamic(gam, chi_w, nd, spatial_shape):
    """The frequency-resolved spin-orbital pairing vertex (issue #284).

    ``chi_w`` ``(nmat, nvol, nd^2, nd^2)`` is the spinful susceptibility at
    every bosonic frequency ``nu_m = (2 m - nmat) pi T`` of the uniform grid
    (bubble ``[a, c, b, d]`` layout, spin-block order). Returns ``V_w``
    ``(nd, nd, nd, nd, Nx, Ny, Nz, nmat)``, the layout of the dynamic kernel
    (``eliashberg_dynamic.eliashberg_kernel_dynamic``):

        V_w[a, b, c, d, q, m] = V_bare[a, b, c, d] + V_fluct[a, b, c, d](q, i nu_m),

    with the parts of :func:`so_pairing_vertex_parts`; the fluctuation part is
    read at ``+q`` and the SAME bosonic index ``m`` as ``chi_w`` (the
    orientation the exact-diagonalization gate of
    ``tests/test_so_pairing_ed.py`` fixes at nonzero frequency; the opposite
    reading ``V[d, c, b, a](-q)`` agrees only at ``nu = 0``). The bare part is
    frequency independent: on the uniform grid it enters every bosonic
    frequency, a ``delta(tau)`` in imaginary time, as in the dynamic kernel's
    other vertices. The on-site check on ``gam`` runs once."""
    chi_w = np.asarray(chi_w)
    if chi_w.ndim != 4 or chi_w.shape[0] < 1:
        raise ValueError("the dynamic spinful susceptibility must have shape "
                         "(nmat, nvol, nd^2, nd^2), got {}".format(chi_w.shape))
    gam, _, nd, shape = _validate(gam, chi_w[0], nd, spatial_shape)
    nmat = chi_w.shape[0]
    nvol = gam.shape[0]
    g4 = gam.reshape(nvol, nd, nd, nd, nd)
    bare = np.ascontiguousarray(np.moveaxis(
        0.25 * (g4 + g4.transpose(0, 4, 3, 2, 1)), 0, -1)).reshape(
            (nd,) * 4 + shape).astype(complex, copy=False)
    out = np.empty((nd,) * 4 + shape + (nmat,), dtype=complex)
    for m in range(nmat):
        fluct = -(gam @ chi_w[m] @ gam).reshape(nvol, nd, nd, nd, nd)
        out[..., m] = bare + np.moveaxis(fluct, 0, -1).reshape((nd,) * 4 + shape)
    return out


# ---------------------------------------------------------------------------
# The static hwave_sc entry in spin-orbital mode: the initial gap and the
# refusal of every combination it does not support
# ---------------------------------------------------------------------------

#: spin structure of the even (singlet-like) and odd (triplet-like) seeds,
#: ``i sigma_y`` and ``sigma_z i sigma_y``
_SEED_SPIN_EVEN = np.array([[0.0, 1.0], [-1.0, 0.0]], dtype=complex)
_SEED_SPIN_ODD = np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)

#: the odd (k -> -k) form factors of ``sc._initialize_gap``
_ODD_SEEDS = ("p_x", "p_y", "p_z", "f_x", "f_y")


def so_initial_gap(init_gap_mode, norb_phys, kx_array, ky_array, kz_array):
    """Initial spin-orbital gap ``(nd, nd, Nx, Ny, Nz)``, ``nd = 2 norb_phys``,
    spin-block order (``A = s * norb_phys + a``).

    The orbital seed of ``sc._initialize_gap`` (``1_orb (x) f(k)``) times a
    spin structure that makes it fermion-antisymmetric
    (``Delta_AB(k) = -Delta_BA(-k)``): ``i sigma_y`` (singlet) for an even
    form factor, ``sigma_z i sigma_y`` (the S_z = 0 triplet) for an odd one
    (``p_*``, ``f_*``). ``"random"`` fills the whole spin-orbital matrix (the
    caller's antisymmetric projection then keeps its admissible part). The
    result has unit norm unless the form factor vanishes on the grid."""
    from hwave import sc
    if init_gap_mode == "random":
        return sc._initialize_gap("random", 2 * int(norb_phys),
                                  kx_array, ky_array, kz_array)
    phi = sc._initialize_gap(init_gap_mode, int(norb_phys),
                             kx_array, ky_array, kz_array)
    spin = _SEED_SPIN_ODD if init_gap_mode in _ODD_SEEDS else _SEED_SPIN_EVEN
    n = int(norb_phys)
    gap = np.einsum("st,ab...->satb...", spin, phi).reshape(
        (2 * n, 2 * n) + phi.shape[2:])
    norm = np.linalg.norm(gap)
    return gap / norm if norm > 0 else gap


def _ci(mapping):
    from requests.structures import CaseInsensitiveDict
    return CaseInsensitiveDict(mapping or {})


_SO_WHAT = "[mode] enable_spin_orbital = true (hwave_sc, issue #83)"
_SO_DYN_WHAT = ("[mode] enable_spin_orbital = true with frequency = \"dynamic\" "
                "(hwave_sc, issue #284)")


def so_frequency(input_dict):
    """``[eliashberg] frequency`` of a spin-orbital run, read
    case-insensitively (as :func:`check_so_config` reads it)."""
    return _ci(input_dict.get("eliashberg", {})).get("frequency", "static")


def _check_so_dynamic_config(eli):
    """The refusals of the spin-orbital dynamic solver (issue #284) on the
    case-insensitive ``[eliashberg]`` table ``eli``, before any file is read.
    Supported: ``chi0q_mode = "flex"`` (a spin-orbital FLEX run with the full
    bosonic axis) or ``"calc"``, ``matsubara_basis = "uniform"``. Refused:
    ``chi0q_mode = "load"``, the IR basis, ``zero_chi_s`` / ``zero_chi_c`` = true
    (there is no spin/charge split of the spinful channel) and the bond-path
    keys ``parity_leakage_policy`` / ``parity_leakage_tol``. (``bond_channels``,
    ``pairing_type`` and ``eigenvalue_method`` are refused by
    :func:`check_so_config` for every frequency; the temperature sweep refuses
    spin-orbital mode in its preflight.)"""
    basis = str(eli.get("matsubara_basis", "uniform")).lower()
    if basis != "uniform":
        raise ValueError(
            "{}: [eliashberg] matsubara_basis = {!r} is not supported; the "
            "spin-orbital dynamic solver runs on the uniform Matsubara grid "
            "only (the IR basis is not available in spin-orbital mode yet). Set "
            "matsubara_basis = \"uniform\" (or remove the key).".format(
                _SO_DYN_WHAT, eli.get("matsubara_basis")))
    mode = eli.get("chi0q_mode", "load")
    if mode not in ("flex", "calc"):
        raise ValueError(
            "{}: [eliashberg] chi0q_mode = {!r} is not supported; use \"flex\" "
            "(the full-frequency chiq, chi0q and green of a spin-orbital FLEX "
            "run) or \"calc\" (the spinful chi0 is computed here at every "
            "bosonic frequency). A loaded chi0q carries no Green function for "
            "the pair bubble.".format(_SO_DYN_WHAT, mode))
    from hwave import sc
    for key in ("zero_chi_s", "zero_chi_c"):
        # strict parse, as on the bond path: a typo is refused, false is fine
        if sc._bond_bool_option(eli, key, False):
            raise ValueError(
                "{}: [eliashberg] {} has no meaning in spin-orbital mode: the "
                "pairing vertex is built from the single spinful channel, which "
                "has no spin/charge split. Remove {} (or set it to false).".format(
                    _SO_DYN_WHAT, key, key))
    for key in ("parity_leakage_policy", "parity_leakage_tol"):
        if key in eli:
            raise ValueError(
                "{}: [eliashberg] {} is a key of the bond-resolved dynamic path; "
                "the spin-orbital solver works in the fermion-antisymmetric space "
                "directly and has no channel sector. Remove {}.".format(
                    _SO_DYN_WHAT, key, key))


def check_so_config(input_dict):
    """Refuse every configuration the spin-orbital entries do not support,
    from the configuration alone (no file is read). Keys are looked up
    case-insensitively; each message names the key and what to change.

    ``frequency = "dynamic"`` (issue #284) adds the refusals of
    :func:`_check_so_dynamic_config`; everything below applies to both.

    Supported: ``frequency = "static"`` or ``"dynamic"``, ``bond_channels = false``,
    ``chi0q_mode = "calc"`` or ``"load"`` with a general (4-index) chi0q, or
    ``"flex"`` with the outputs of a spin-orbital FLEX run (``flex_chiq``,
    ``flex_green``; ``flex_chi_s`` / ``flex_chi_c`` are refused),
    ``solver_mode = "iteration" / "eigenvalue" / "both"`` with
    ``eigenvalue_method = "arnoldi"``, no ``pairing_type`` (the gap is one
    general spin-orbital matrix) and the antisymmetrized spinful vertex
    (``spinful_vertex_exchange`` not false). An ``Extern`` file (a Zeeman
    field ``coeff_extern * sigma^z H_ab(R)`` in physical orbital indices,
    issue #88) is accepted; it enters the band through
    :func:`hwave.sc._so_epsilon_k` and is checked by
    :func:`hwave.sc._read_so_extern`."""
    from hwave.solver.backend import strict_bool
    eli = _ci(input_dict.get("eliashberg", {}))
    freq = eli.get("frequency", "static")
    if freq not in ("static", "dynamic"):
        raise ValueError(
            "{}: [eliashberg] frequency = {!r} is not supported; use \"static\" "
            "or \"dynamic\" (or remove the key).".format(_SO_WHAT, freq))
    bond = eli.get("bond_channels")
    if bond is not None:
        from hwave import sc
        if sc._bond_bool_option(eli, "bond_channels", False):
            raise ValueError(
                "{}: [eliashberg] bond_channels = true is not supported; the "
                "spin-orbital pairing vertex is on-site only. Set bond_channels = "
                "false.".format(_SO_WHAT))
    if "pairing_type" in eli:
        raise ValueError(
            "{}: [eliashberg] pairing_type = {!r} has no meaning in spin-orbital "
            "mode: the gap is one general 2N x 2N spin-orbital matrix, not a "
            "singlet or a triplet channel. Remove pairing_type.".format(
                _SO_WHAT, eli["pairing_type"]))
    if freq == "dynamic":
        _check_so_dynamic_config(eli)
    method = eli.get("eigenvalue_method", "arnoldi")
    if method != "arnoldi":
        raise ValueError(
            "{}: [eliashberg] eigenvalue_method = {!r} is not supported; only "
            "\"arnoldi\" is. Set eigenvalue_method = \"arnoldi\" (or remove the "
            "key).".format(_SO_WHAT, method))
    from hwave.sc import DEFAULT_SOLVER_MODE
    solver_mode = eli.get("solver_mode", DEFAULT_SOLVER_MODE)
    if solver_mode not in ("iteration", "eigenvalue", "both"):
        raise ValueError(
            "{}: [eliashberg] solver_mode = {!r} is not supported; use "
            "\"iteration\", \"eigenvalue\" or \"both\".".format(_SO_WHAT, solver_mode))
    chi0q_mode = eli.get("chi0q_mode", "load")
    if chi0q_mode not in ("calc", "load", "flex"):
        raise ValueError(
            "{}: [eliashberg] chi0q_mode = {!r} is not supported; use \"calc\" "
            "(the spinful chi0 is computed here), \"load\" (a general chi0q of "
            "a spin-orbital RPA run) or \"flex\" (the chiq and green of a "
            "spin-orbital FLEX run).".format(_SO_WHAT, chi0q_mode))
    for key in ("flex_chi_s", "flex_chi_c"):
        if key in eli:
            raise ValueError(
                "{}: [eliashberg] {} names a spin or charge channel file, which a "
                "spin-orbital FLEX run does not write (it stores the single "
                "spinful channel as chiq). Remove {} and set flex_chiq to the "
                "chiq archive of the spin-orbital FLEX run (default "
                "\"chiq.npz\").".format(_SO_WHAT, key, key))
    tensor = eli.get("chi0q_tensor", "auto")
    if tensor not in ("auto", "general"):
        raise ValueError(
            "{}: [eliashberg] chi0q_tensor = {!r} is not supported; the "
            "spin-orbital vertex needs the general (4-index) chi0q. Set "
            "chi0q_tensor = \"general\" (or \"auto\").".format(_SO_WHAT, tensor))
    param = _ci(input_dict.get("mode", {}).get("param", {}))
    sve = param.get("spinful_vertex_exchange")
    if sve is not None and not strict_bool(sve, "[mode.param] spinful_vertex_exchange"):
        raise ValueError(
            "{}: [mode.param] spinful_vertex_exchange = false is not supported; "
            "the spin-orbital pairing vertex is built from the antisymmetrized "
            "spinful vertex. Set spinful_vertex_exchange = true (or remove the "
            "key).".format(_SO_WHAT))
    finp = _ci(input_dict.get("file", {}).get("input", {}))
    for key in ("trans_mod", "green_init"):
        if key in finp:
            raise ValueError(
                "{}: [file.input] {} is not supported: the spin-orbital "
                "Eliashberg solver does not support a mean-field one-body input "
                "yet (the susceptibility would be built from that one-body "
                "Hamiltonian while the kernel Green function and chemical "
                "potential come from the bare Transfer file; see issue #273). "
                "Remove {}.".format(_SO_WHAT, key, key))


def check_so_supported(input_dict, interactions, site, norb_phys=None):
    """:func:`check_so_config` plus the interaction tables.

    ``interactions`` are the ORIGINAL declarations ``{type: {(R, (a, b)): v}}``
    in physical orbital indices (as read, before any folding: an off-site row
    can fold onto ``R = 0``), ``site`` the ``[mode.param] orbital_site``
    labels of the physical orbitals or ``None``. Refused: any nonzero row
    with ``R != 0`` (also inside a ``Coulomb`` aggregate file) and, with
    ``orbital_site``, any ``R = 0`` row between orbitals of different sites:
    the spin-orbital pairing vertex supports on-site interactions only.
    Every orbital index must satisfy ``0 <= o < norb_phys`` (``o >= 0``
    when ``norb_phys`` is not given): an input orbital 0 (reads as -1) would
    otherwise wrap to the last orbital in the vertex, and an index past the
    physical orbitals means a file written in spin-orbital indices."""
    from hwave.solver.sites import is_onsite
    check_so_config(input_dict)
    files = _ci(input_dict.get("file", {}).get("input", {}).get("interaction", {}))
    aggregate = "Coulomb" in files
    for t, tbl in (interactions or {}).items():
        for (irvec, orbvec), v in (tbl or {}).items():
            bad = [int(o) for o in orbvec if int(o) < 0
                   or (norb_phys is not None and int(o) >= int(norb_phys))]
            if bad:
                fname = files.get("Coulomb" if (
                    aggregate and t in ("CoulombIntra", "CoulombInter")
                    and t not in files) else t, "")
                rng = ("1..{}".format(int(norb_phys)) if norb_phys is not None
                       else "starting at 1")
                raise ValueError(
                    "{}: the {} file{} has a row on orbital {}, outside the "
                    "physical orbitals {}: in spin-orbital mode the interaction "
                    "files use 1-based PHYSICAL orbital indices{} (the Geometry "
                    "norb and the Transfer file count spin-orbitals).".format(
                        _SO_WHAT, t, " '{}'".format(fname) if fname else "",
                        bad[0] + 1, rng,
                        "" if norb_phys is None else
                        ", and the model has {}".format(int(norb_phys))))
            if v == 0 or is_onsite(irvec, orbvec, site):
                continue
            R = tuple(int(x) for x in irvec)
            src = ("the Coulomb file (split into {})".format(t)
                   if aggregate and t in ("CoulombIntra", "CoulombInter")
                   else "the {} file".format(t))
            if R != (0, 0, 0):
                why = "an off-site row (R = {}, orbitals {}) in {}".format(
                    R, tuple(int(o) + 1 for o in orbvec), src)
                fix = "remove the off-site rows"
            else:
                why = ("an R = 0 row between orbitals {} of different sites "
                       "([mode.param] orbital_site) in {}".format(
                           tuple(int(o) + 1 for o in orbvec), src))
                fix = "remove the row or the orbital_site labels"
            raise ValueError(
                "{}: {} is not supported; the spin-orbital pairing vertex supports "
                "on-site interactions only (the exchange crossing of an off-site "
                "term depends on k + k'). {}.".format(_SO_WHAT, why, fix[0].upper() + fix[1:]))


# ---------------------------------------------------------------------------
# Spin structure of a spin-orbital gap (output diagnostics)
# ---------------------------------------------------------------------------

#: names of the spin components, in output order
SPIN_COMPONENTS = ("psi", "dx", "dy", "dz")

_PAULI = (np.eye(2, dtype=complex),
          np.array([[0, 1], [1, 0]], dtype=complex),
          np.array([[0, -1j], [1j, 0]], dtype=complex),
          np.array([[1, 0], [0, -1]], dtype=complex))
#: ``i sigma_y``
_I_SIGMA_Y = 1j * _PAULI[2]
#: the spin basis ``sigma_mu i sigma_y`` (mu = 0, x, y, z) of the
#: decomposition; orthogonal, ``Tr[B_mu^dagger B_nu] = 2 delta_mu,nu``
_SPIN_BASIS = tuple(m @ _I_SIGMA_Y for m in _PAULI)


def _split_spin_blocks(gap, norb_phys):
    gap = np.asarray(gap)
    n = int(norb_phys)
    if n <= 0 or gap.ndim < 2 or gap.shape[0] != 2 * n or gap.shape[1] != 2 * n:
        raise ValueError(
            "a spin-orbital gap must have shape (2 N, 2 N, ...) with N = "
            "norb_phys = {}; got {}".format(norb_phys, gap.shape))
    return gap.reshape((2, n, 2, n) + gap.shape[2:])


def so_spin_components(gap, norb_phys):
    """Singlet and triplet components of a spin-orbital gap.

    ``gap`` is ``Delta_AB(...)`` of shape ``(2N, 2N, ...)`` in spin-block
    order (``A = s * N + a``, spin 0 = up; the trailing axes are the k grid,
    or anything else). For every orbital pair ``(a, b)`` (and k) the 2 x 2
    spin block ``[Delta_ab]_st = Delta_{(s,a),(t,b)}`` is written as

        Delta_ab = (psi_ab + d_ab . sigma) i sigma_y,

    i.e. ``psi_ab = 1/2 Tr[(i sigma_y)^dagger Delta_ab]`` and
    ``d_ab,j = 1/2 Tr[(sigma_j i sigma_y)^dagger Delta_ab]``; explicitly

        psi = (D_ud - D_du) / 2,      d_x = (D_dd - D_uu) / 2,
        d_y = (D_uu + D_dd) / (2 i),  d_z = (D_ud + D_du) / 2.

    Returns ``(psi, dx, dy, dz)``, complex arrays of shape ``(N, N, ...)``.
    The four matrices ``sigma_mu i sigma_y`` are orthogonal with norm
    ``sqrt(2)``, so ``|Delta_ab|_F^2 = 2 (|psi_ab|^2 + |d_ab|^2)``."""
    g = _split_spin_blocks(gap, norb_phys)
    return tuple(0.5 * np.einsum("st,satb...->ab...", b.conj(), g)
                 for b in _SPIN_BASIS)


def so_spin_decomposition(gap, norb_phys):
    """Spin weights of a spin-orbital gap: ``{"psi", "dx", "dy", "dz"}``.

    With the components of :func:`so_spin_components`
    (``Delta_ab(k) = (psi_ab(k) + d_ab(k) . sigma) i sigma_y`` per orbital
    pair), each weight is the sum over k and the orbital pairs ``(a, b)`` of
    the component's ``|.|^2``, normalized so that the four sum to 1:

        w_psi = sum_{k,a,b} |psi_ab(k)|^2 / W,  w_dj = sum_{k,a,b} |d_ab,j(k)|^2 / W,
        W = sum_{k,a,b} (|psi_ab(k)|^2 + |d_ab(k)|^2) = |Delta|_F^2 / 2.

    ``w_psi = 1`` is a pure (spin-singlet) gap, ``w_dx + w_dy + w_dz = 1`` a
    pure spin-triplet one. The weights do not depend on the overall phase or
    normalization of the gap. A zero gap gives NaN weights."""
    comps = so_spin_components(gap, norb_phys)
    sq = [float(np.sum(np.abs(c) ** 2)) for c in comps]
    tot = sum(sq)
    if not tot > 0.0:
        return {k: float("nan") for k in SPIN_COMPONENTS}
    return {k: v / tot for k, v in zip(SPIN_COMPONENTS, sq)}


# ---------------------------------------------------------------------------
# S_z content of a spin-orbital gap and the S_z-conservation check (#88)
# ---------------------------------------------------------------------------

#: names of the S_z components, in output order
SZ_COMPONENTS = ("upup", "dndn", "singlet0", "triplet0")

#: ``sz_conserved`` holds when the residual is at most this
SZ_CONSERVED_TOL = 1.0e-12


def so_sz_decomposition(gap, norb_phys):
    """S_z weights of a spin-orbital gap:
    ``{"upup", "dndn", "singlet0", "triplet0"}``.

    ``gap`` is ``Delta_AB(...)`` of shape ``(2N, 2N, ...)`` in spin-block
    order (``A = s * N + a``, spin 0 = up along the INPUT quantization axis:
    the sigma^z basis of the Transfer spin labels and of Extern's sigma^z;
    the weights are not rotated to a field or magnetization axis). With the
    spin blocks ``D_st`` of :func:`so_spin_components`,
    ``psi = (D_ud - D_du) / 2`` and ``d_z = (D_ud + D_du) / 2``:

        W = sum |Delta|^2,  w_upup = sum |D_uu|^2 / W,  w_dndn = sum |D_dd|^2 / W,
        w_singlet0 = 2 sum |psi|^2 / W,  w_triplet0 = 2 sum |d_z|^2 / W,

    summed over every orbital pair and trailing axis; the four sum to 1
    (``|D_ud|^2 + |D_du|^2 = 2 (|psi|^2 + |d_z|^2)``). In terms of
    :func:`so_spin_decomposition`, ``w_singlet0 = w_psi``,
    ``w_triplet0 = w_dz`` and ``w_upup + w_dndn = w_dx + w_dy``. A zero gap
    gives NaN weights."""
    g = _split_spin_blocks(gap, norb_phys)

    def sq(x):
        return float(np.sum(np.abs(x) ** 2))

    ud, du = g[0, :, 1], g[1, :, 0]
    parts = (sq(g[0, :, 0]), sq(g[1, :, 1]),
             0.5 * sq(ud - du), 0.5 * sq(ud + du))
    tot = parts[0] + parts[1] + sq(ud) + sq(du)
    if not tot > 0.0:
        return {k: float("nan") for k in SZ_COMPONENTS}
    return {k: v / tot for k, v in zip(SZ_COMPONENTS, parts)}


def _sz_spin_tuples(conserving):
    """The spin labels ``(s_a, s_b, s_c, s_d)`` of the vertex
    ``V_abcd`` (``Delta_ad = -sum_bc V_abcd F_bc``) inside
    (``conserving``) or outside the S_z sector ``S(a)+S(d) = S(b)+S(c)``."""
    return tuple((sa, sb, sc_, sd)
                 for sa in (0, 1) for sb in (0, 1) for sc_ in (0, 1)
                 for sd in (0, 1)
                 if (sa + sd == sb + sc_) == conserving)


def sz_off_sector_index(norb_phys):
    """Every off-sector vertex index ``(a, b, c, d)`` (spin-block order,
    ``S(A) = A // norb_phys``): ``S(a) + S(d) != S(b) + S(c)``, i.e. ``V_abcd``
    maps a pair of one S_z onto a pair of another. 10 of the 16 spin blocks."""
    n = int(norb_phys)
    if n <= 0:
        raise ValueError("norb_phys must be positive, got {}".format(norb_phys))
    out = []
    for spins in _sz_spin_tuples(False):
        for orbs in np.ndindex(n, n, n, n):
            out.append(tuple(int(s * n + o) for s, o in zip(spins, orbs)))
    return tuple(out)


def _block_max(x, n, spins):
    """``max |x[block]|`` over the spin blocks ``spins`` of the leading
    ``len(spins[0])`` axes; every block is a basic-slicing view (no copy of
    ``x``), and the temporaries are one block at a time. A non-finite entry
    raises ``ValueError``: it would otherwise drop out of ``max`` and report
    a broken input as S_z-conserving."""
    best = 0.0
    for s in spins:
        blk = x[tuple(slice(t * n, (t + 1) * n) for t in s)]
        if blk.size:
            m = float(np.abs(blk).max())
            if not np.isfinite(m):
                raise ValueError(
                    "the S_z conservation residual needs finite input; the "
                    "spin block {} holds a non-finite entry".format(s))
            best = max(best, m)
    return best


def _leak(x, n, off, on):
    """``max |x_off-sector| / max |x|`` (max-norm; 0 if ``max |x| = 0``)."""
    x_off = _block_max(x, n, off)
    scale = max(x_off, _block_max(x, n, on))
    return x_off / scale if scale > 0.0 else 0.0


def sz_conservation_residual(V=None, green=None, norb_phys=None):
    """The S_z-conservation residual ``max(r_V, r_G)`` of a spin-orbital run.

    * ``r_V``: the vertex ``V`` ``(nd, nd, nd, nd, ...)`` (spin-block order,
      any trailing k / frequency axes), off-sector where
      ``S(a) + S(d) != S(b) + S(c)`` (:func:`sz_off_sector_index`);
    * ``r_G``: the Green function ``green`` ``(nd, nd, ...)``, off-sector in
      the spin-off-diagonal blocks ``G[:N, N:]`` and ``G[N:, :N]``;

    each ``max |X_off| / max |X|`` over all blocks and slices (0 for a zero
    array). A ``None`` argument contributes 0. The blocks are read through
    views; neither array is copied."""
    if norb_phys is None:
        raise ValueError("sz_conservation_residual needs norb_phys")
    n = int(norb_phys)
    nd = 2 * n
    r = 0.0
    if V is not None:
        V = np.asarray(V)
        if V.ndim < 4 or V.shape[:4] != (nd,) * 4:
            raise ValueError("the vertex must have shape (nd, nd, nd, nd, ...) "
                             "with nd = {}; got {}".format(nd, V.shape))
        r = max(r, _leak(V, n, _sz_spin_tuples(False), _sz_spin_tuples(True)))
    if green is not None:
        green = np.asarray(green)
        if green.ndim < 2 or green.shape[:2] != (nd, nd):
            raise ValueError("the Green function must have shape (nd, nd, ...) "
                             "with nd = {}; got {}".format(nd, green.shape))
        r = max(r, _leak(green, n, ((0, 1), (1, 0)), ((0, 0), (1, 1))))
    return float(r)


def sz_report(weights, residual, logger):
    """The S_z report of a spin-orbital solve. ``weights`` holds the
    :func:`so_sz_decomposition` of each reported gap (eigenvector ``i``);
    logs one INFO line per eigenvector and, when
    ``residual > SZ_CONSERVED_TOL``, one WARNING; returns the archive members
    ``sz_weight_<component>`` (float64 ``(n_ev,)``),
    ``sz_conservation_residual`` (float64) and ``sz_conserved`` (bool_)."""
    rows = []
    for i, w in enumerate(weights):
        rows.append([float(w[c]) for c in SZ_COMPONENTS])
        logger.info(
            "S_z decomposition of eigenvector %d: up-up = %.6f, down-down = "
            "%.6f, S_z=0 singlet = %.6f, S_z=0 triplet = %.6f (the four sum to "
            "1)", i, *rows[-1])
    residual = float(residual)
    conserved = residual <= SZ_CONSERVED_TOL
    if not conserved:
        logger.warning(
            "S_z is not conserved (residual %.3e): the S_z weights describe "
            "the composition of each eigenvector in the input quantization "
            "axis, not separate pairing channels.", residual)
    rows = np.array(rows, dtype=np.float64).reshape(-1, len(SZ_COMPONENTS))
    out = {"sz_weight_" + c: rows[:, j].copy() for j, c in enumerate(SZ_COMPONENTS)}
    out["sz_conservation_residual"] = np.float64(residual)
    out["sz_conserved"] = np.bool_(conserved)
    return out
