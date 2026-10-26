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

__all__ = ["so_pairing_vertex", "so_pairing_vertex_parts", "so_initial_gap",
           "check_so_config", "check_so_supported"]

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


def check_so_config(input_dict):
    """Refuse every configuration the spin-orbital static entry does not
    support, from the configuration alone (no file is read). Keys are looked
    up case-insensitively; each message names the key and what to change.

    Supported: ``frequency = "static"``, ``bond_channels = false``,
    ``chi0q_mode = "calc"`` or ``"load"`` with a general (4-index) chi0q, or
    ``"flex"`` with the outputs of a spin-orbital FLEX run (``flex_chiq``,
    ``flex_green``; ``flex_chi_s`` / ``flex_chi_c`` are refused),
    ``solver_mode = "iteration" / "eigenvalue" / "both"`` with
    ``eigenvalue_method = "arnoldi"``, no ``pairing_type`` (the gap is one
    general spin-orbital matrix), the antisymmetrized spinful vertex
    (``spinful_vertex_exchange`` not false) and no ``Extern`` file."""
    from hwave.solver.backend import strict_bool
    eli = _ci(input_dict.get("eliashberg", {}))
    freq = eli.get("frequency", "static")
    if freq != "static":
        raise ValueError(
            "{}: [eliashberg] frequency = {!r} is not supported; the spin-orbital "
            "Eliashberg solver is static only. Set frequency = \"static\" (or remove "
            "the key).".format(_SO_WHAT, freq))
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
    method = eli.get("eigenvalue_method", "arnoldi")
    if method != "arnoldi":
        raise ValueError(
            "{}: [eliashberg] eigenvalue_method = {!r} is not supported; only "
            "\"arnoldi\" is. Set eigenvalue_method = \"arnoldi\" (or remove the "
            "key).".format(_SO_WHAT, method))
    solver_mode = eli.get("solver_mode", "iteration")
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
    files = _ci(input_dict.get("file", {}).get("input", {}).get("interaction", {}))
    if "Extern" in files:
        raise ValueError(
            "{}: the Extern (external field) file is not used by hwave_sc. In "
            "spin-orbital mode write a Zeeman term into the Transfer file as a "
            "spin-dependent one-body term, and remove Extern.".format(_SO_WHAT))


def check_so_supported(input_dict, interactions, site, norb_phys=None):
    """:func:`check_so_config` plus the interaction tables.

    ``interactions`` are the ORIGINAL declarations ``{type: {(R, (a, b)): v}}``
    in physical orbital indices (as read, before any folding: an off-site row
    can fold onto ``R = 0``), ``site`` the ``[mode.param] orbital_site``
    labels of the physical orbitals or ``None``. Refused: any nonzero row
    with ``R != 0`` (also inside a ``Coulomb`` aggregate file) and, with
    ``orbital_site``, any ``R = 0`` row between orbitals of different sites:
    the spin-orbital pairing vertex supports on-site interactions only.
    With ``norb_phys`` given, an orbital index outside the physical orbitals
    (a file written in spin-orbital indices) is refused too."""
    from hwave.solver.sites import is_onsite
    check_so_config(input_dict)
    files = _ci(input_dict.get("file", {}).get("input", {}).get("interaction", {}))
    aggregate = "Coulomb" in files
    for t, tbl in (interactions or {}).items():
        for (irvec, orbvec), v in (tbl or {}).items():
            if norb_phys is not None and max(int(o) for o in orbvec) >= int(norb_phys):
                raise ValueError(
                    "{}: the {} file has a row on orbital {} but the model has {} "
                    "physical orbitals: in spin-orbital mode the interaction files "
                    "use PHYSICAL orbital indices 1..{} (the Geometry norb and the "
                    "Transfer file count spin-orbitals).".format(
                        _SO_WHAT, t, max(int(o) for o in orbvec) + 1, norb_phys,
                        norb_phys))
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
