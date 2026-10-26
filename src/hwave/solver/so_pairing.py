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

__all__ = ["so_pairing_vertex", "so_pairing_vertex_parts"]

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
