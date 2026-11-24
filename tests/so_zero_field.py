"""The weak-field limit h -> 0 of a spin-orbital pairing kernel without
spin-orbit coupling: the measures shared by the gates of the on-site path
(``tests/test_so_zeeman_route.py``) and of the bond path
(``tests/test_sc_spin_orbital_bond.py``).

A weak Zeeman field moves the eigenvalues of an SU(2) kernel linearly: an
SU(2) triplet (degeneracy 3m) splits into m slopes near zero (pair S_z = 0)
and m pairs of opposite slopes (up-up / down-down), a singlet does not move
at first order. The eigenvalue paths over the fields ``H`` (``H[0] = 0``,
increasing) are tracked by assignment and fitted by ``a h + b h^2``.

The fit describes a cluster only while the field does not mix it with
another cluster: besides the slope separation (the paths are tracked
without a crossing), a selected cluster ``i`` needs, for every other
cluster ``j``,

``r_ij = h_4 sqrt(||C_ij||_2 ||C_ji||_2) / |lambda_i - lambda_j| <= 0.03``

with ``C = W (K(h_1) - K(h_0)) V / h_1`` the first-order field coupling in
the biorthogonal eigenbasis of ``K(h_0) = V diag(lambda) W`` (``W = V^-1``),
``C_ij`` its block between the clusters. In a two-level model the first
term the fit neglects is of relative size ``r^2``, and condition 1 allows
``1e-3`` of the deviation, so ``r <= 0.03``. A fixture with a larger ``r``
on a selected cluster is rejected (its selected cluster is not isolated).
"""

import numpy as np

from hwave.solver import bond_channels, so_pairing

#: the leading clusters by |lambda| that conditions 1, 3 and 4 inspect
FL_LEAD = 16
#: the largest field-induced mixing ratio of a selected cluster
ISOLATION_MAX = 0.03


def track(lams):
    """Eigenvalue paths ``(n_ev, n_h)``: the eigenvalues at each field are
    assigned one to one (minimum total distance) to the paths at the previous
    field, starting from the order at the first field."""
    from scipy.optimize import linear_sum_assignment
    paths = np.empty((lams.shape[1], lams.shape[0]), dtype=complex)
    paths[:, 0] = lams[0]
    for k in range(1, lams.shape[0]):
        cost = np.abs(paths[:, k - 1][:, None] - lams[k][None, :])
        r, c = linear_sum_assignment(cost)
        paths[r, k] = lams[k][c]
    return paths


def fit(paths, H):
    """Least squares ``lambda(h) - lambda(0) = a h + b h^2`` on the nonzero
    fields: complex ``a``, ``b`` and the max absolute residual per path."""
    h = np.asarray(H[1:], float)
    A = np.stack([h, h ** 2], axis=1)
    y = (paths[:, 1:] - paths[:, :1]).T                    # (n_h - 1, n_ev)
    coef = np.linalg.lstsq(A.astype(complex), y, rcond=None)[0]
    res = np.abs(A @ coef - y).max(axis=0)
    return coef[0], coef[1], res


def classify_clusters(clusters, vecs0, gap_shape, norb_phys):
    """``"singlet"`` / ``"triplet"`` per cluster from the spin weights of its
    whole (orthonormalized) eigenspace, so the mixing of degenerate
    eigenvectors does not matter: the mean of
    :func:`so_pairing.so_spin_decomposition` over an orthonormal basis of
    the span (``w_psi >= 0.99`` singlet-like, ``w_d >= 0.99`` triplet-like);
    a mixed cluster is an ``AssertionError`` (the fixture is rejected).
    ``vecs0``: the eigenvectors at ``H[0]`` in the full flat gap (columns)."""
    kinds, weights = [], []
    for c in clusters:
        Qc, _ = np.linalg.qr(vecs0[:, list(c)])
        w = {"psi": 0.0, "d": 0.0}
        for j in range(Qc.shape[1]):
            sw = so_pairing.so_spin_decomposition(Qc[:, j].reshape(gap_shape), norb_phys)
            w["psi"] += sw["psi"] / len(c)
            w["d"] += (sw["dx"] + sw["dy"] + sw["dz"]) / len(c)
        weights.append(w)
        if w["psi"] >= 0.99:
            kinds.append("singlet")
        elif w["d"] >= 0.99:
            kinds.append("triplet")
        else:
            raise AssertionError("cluster %s is spin-mixed: %r" % (c, w))
    return kinds, weights


def cluster_mixing(K0, K1, h1, h4, lam0, clusters):
    """``r[i, j]`` (module docstring) for every pair of the ``clusters`` of
    ``lam0`` (the eigenvalues at ``h_0 = 0``), from the reduced kernels ``K0
    = K(h_0)`` and ``K1 = K(h_1)`` in one set of coordinates; the diagonal is
    0. The eigenvalues of ``K0`` are assigned one to one to ``lam0`` (they
    must agree to ``1e-8 max(1, max|lambda|)``: ``lam0`` may come from
    another solve of the same kernel, e.g. a public run's listing)."""
    from scipy.optimize import linear_sum_assignment
    lam, V = np.linalg.eig(K0)
    if lam.shape != np.shape(lam0):
        raise AssertionError("cluster_mixing: K0 has {} eigenvalues, the paths {}".format(
            lam.shape, np.shape(lam0)))
    row, col = linear_sum_assignment(np.abs(np.asarray(lam0)[:, None] - lam[None, :]))
    miss = float(np.abs(np.asarray(lam0)[row] - lam[col]).max())
    if miss > 1e-8 * max(1.0, float(np.abs(lam).max())):
        raise AssertionError("cluster_mixing: K0 does not give the eigenvalues of the paths "
                             "(miss {:.3e})".format(miss))
    of = dict(zip(row.tolist(), col.tolist()))                 # path index -> eig index
    C = np.linalg.solve(V, (np.asarray(K1) - np.asarray(K0)) @ V) / h1
    idx = [[of[i] for i in c] for c in clusters]
    n = len(clusters)
    r = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            cij = np.linalg.norm(C[np.ix_(idx[i], idx[j])], 2)
            cji = np.linalg.norm(C[np.ix_(idx[j], idx[i])], 2)
            r[i, j] = h4 * np.sqrt(cij * cji) / abs(lam0[clusters[i][0]] - lam0[clusters[j][0]])
    return r


def zero_field_measures(lams, vecs0, gap_shape, norb_phys, H, *, kernels=None):
    """The measures of the weak-field conditions on the eigenvalue paths of
    ``lams`` (rows: the fields ``H``): the paths tracked by assignment
    (:func:`track`), the fit ``lambda_i(h) - lambda_i(0) = a_i h + b_i h^2``
    on the nonzero fields (:func:`fit`); clusters at ``H[0]`` with the
    absolute ``deg_tol = 1e-8 max(1, max|lambda|)``, classified by their spin
    weights (:func:`classify_clusters`). With ``scale = max(1,
    max|lambda|)`` the conditions are: (1) ``fit <= 1``, each fit residual
    over ``max(1e-3 max_j |lambda_i(h_j) - lambda_i(0)|, 1e-11 scale)``; (2)
    ``floor = max_i |a_i| h_4 / scale >= 1e-6``; (3) ``triplet_split``: every
    triplet-like cluster has size 3m and exactly m slopes with ``|a| <= 1e-3
    max_cluster |a|``, and ``triplet_pair <= 1e-2``, the other 2m slopes
    paired with their negatives (assignment on ``|a_i + a_j|``, no self
    pairing), over ``max_cluster |a|``; (4) ``singlet <= 1e-3``, the
    singlet-like slopes over ``max_all |a|``; (5) ``shift = max_i
    |lambda_i(h_4) - lambda_i(0)| / (h_4 scale) <= 10``. Conditions 2 and 5
    are global; 1, 3 and 4 apply to the leading clusters by ``|lambda|`` (at
    most :data:`FL_LEAD`) that are isolated, farther than ``10 (max_c |a| +
    max_c' |a|) h_4`` from every other cluster, so that their paths are
    tracked without a crossing (``selected`` lists them; a triplet-like one
    must be among them).

    ``kernels = (K0, K1)``, the reduced kernels at ``H[0]`` and ``H[1]``,
    adds ``isolation``: the largest :func:`cluster_mixing` ratio of a
    selected cluster with any other cluster (``<=`` :data:`ISOLATION_MAX`,
    else the fixture is rejected); without them ``isolation`` is None."""
    from scipy.optimize import linear_sum_assignment
    scale = max(1.0, float(np.abs(lams[0]).max()))
    paths = track(lams)
    a, _, res = fit(paths, H)
    clusters = bond_channels.cluster_eigenvalues(lams[0], deg_tol=1e-8 * scale)
    kinds, _ = classify_clusters(clusters, vecs0, gap_shape, norb_phys)
    amax_all = float(np.abs(a).max())
    amax = [float(np.abs(a[list(c)]).max()) for c in clusters]
    lead = sorted(range(len(clusters)), key=lambda i: -abs(lams[0][clusters[i][0]]))[:FL_LEAD]
    selected = [i for i in lead
                if all(abs(lams[0][clusters[i][0]] - lams[0][clusters[j][0]])
                       > 10 * (amax[i] + amax[j]) * H[-1]
                       for j in range(len(clusters)) if j != i)]
    out = {"floor": amax_all * H[-1] / scale,
           "shift": float(np.abs(paths[:, -1] - paths[:, 0]).max()) / (H[-1] * scale),
           "fit": 0.0, "triplet_split": True, "triplet_zero": 0.0, "triplet_pair": 0.0,
           "singlet": 0.0, "kinds": tuple(sorted((k, len(c)) for k, c in zip(kinds, clusters))),
           "selected": tuple((kinds[i], len(clusters[i])) for i in selected),
           "clusters": len(clusters), "isolation": None}
    if kernels is not None:
        r = cluster_mixing(kernels[0], kernels[1], H[1], H[-1], lams[0], clusters)
        out["isolation"] = float(max((r[i].max() for i in selected), default=0.0))
    for i in selected:
        c, kind = list(clusters[i]), kinds[i]
        dev = np.abs(paths[c] - paths[c, :1]).max(axis=1)
        out["fit"] = max(out["fit"], float((res[c] / np.maximum(1e-3 * dev, 1e-11 * scale)).max()))
        ac = a[c]
        if kind == "triplet":
            m = len(c) // 3
            small = np.argsort(np.abs(ac))[:m]
            if len(c) % 3 or int(np.sum(np.abs(ac) <= 1e-3 * amax[i])) != m:
                out["triplet_split"] = False
                continue
            out["triplet_zero"] = max(out["triplet_zero"], float(np.abs(ac[small]).max()) / amax[i])
            rest = np.delete(ac, small)
            cost = np.abs(rest[:, None] + rest[None, :])
            np.fill_diagonal(cost, np.inf)
            rr, cc = linear_sum_assignment(cost)
            out["triplet_pair"] = max(out["triplet_pair"], float(cost[rr, cc].max()) / amax[i])
        else:
            out["singlet"] = max(out["singlet"], amax[i] / amax_all)
    return out


def assert_zero_field(test, m, *, need_isolation=True):
    """The five weak-field conditions on the measures of
    :func:`zero_field_measures`, with a triplet-like cluster among the
    selected ones and, first, the isolation of the selected clusters (a
    fixture whose selected cluster is mixed by the field is rejected;
    ``need_isolation`` asserts that it was measured)."""
    if need_isolation:
        test.assertIsNotNone(m["isolation"], m)
    if m["isolation"] is not None:
        test.assertLessEqual(m["isolation"], ISOLATION_MAX,
                             "the fixture is rejected: a selected cluster is mixed by the "
                             "field (%r)" % (m,))
    test.assertIn("triplet", {k for k, _ in m["selected"]}, m)
    test.assertLessEqual(m["fit"], 1.0, m)                                       # 1
    test.assertGreaterEqual(m["floor"], 1e-6, m)                                 # 2
    test.assertTrue(m["triplet_split"], m)                                       # 3
    test.assertLessEqual(m["triplet_pair"], 1e-2, m)
    test.assertLessEqual(m["singlet"], 1e-3, m)                                  # 4
    test.assertLessEqual(m["shift"], 10.0, m)                                    # 5
