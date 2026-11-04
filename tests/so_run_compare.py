"""Comparing the outputs of two runs that must give the same result.

Two runs of the same input are deterministic only up to the floating-point
reduction order of the linear-algebra library (an OpenBLAS build can change
the last bits from one call to the next). The comparisons here therefore
use tolerances instead of the stored bytes:

* eigenvalues agree to ``1e-12 * max(1, |lambda|)``;
* eigenvectors are compared through the subspaces they span: the leading
  gap and the stored eigenvectors carry an arbitrary phase, and a
  degenerate set an arbitrary basis, so the cosines of the principal angles
  between the two spans (the singular values of ``Qa^H Qb`` for orthonormal
  bases) must all be 1 to ``1e-10``;
* every other floating-point member agrees to ``1e-13 * max(1, max|x|)``;
  integer and boolean members agree exactly.
"""

import numpy as np

#: eigenvalue tolerance, relative to max(1, |lambda|)
EIGENVALUE_TOL = 1e-12
#: eigenvalues closer than this (relative to max(1, |lambda|)) form one
#: degenerate set
DEGENERACY_TOL = 1e-8
#: tolerance on the cosines of the principal angles between two spans
SUBSPACE_TOL = 1e-10
#: tolerance of the other floating-point members, relative to max(1, max|x|)
VALUE_TOL = 1e-13


def eigenvalues_close(a, b, tol=EIGENVALUE_TOL):
    """True when ``a`` and ``b`` agree elementwise to ``tol * max(1, |a|)``."""
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    return bool(np.all(np.abs(a - b) <= tol * np.maximum(1.0, np.abs(a))))


def degenerate_sets(vals, tol=DEGENERACY_TOL):
    """The index lists of the eigenvalues of ``vals`` that coincide to
    ``tol * max(1, |lambda|)``, in the order of ``vals``."""
    vals = np.asarray(vals).reshape(-1)
    sets, seen = [], set()
    for i in range(len(vals)):
        if i in seen:
            continue
        scale = tol * max(1.0, abs(vals[i]))
        group = [j for j in range(len(vals))
                 if j not in seen and abs(vals[j] - vals[i]) <= scale]
        seen.update(group)
        sets.append(group)
    return sets


def _orthonormal(rows):
    q, _ = np.linalg.qr(np.asarray(rows, dtype=complex).reshape(len(rows), -1).T)
    return q


def same_span(a_rows, b_rows, tol=SUBSPACE_TOL):
    """True when the vectors ``a_rows`` and ``b_rows`` (one vector per
    leading index) span the same subspace: every singular value of
    ``Qa^H Qb`` is 1 to ``tol``."""
    a_rows, b_rows = np.asarray(a_rows), np.asarray(b_rows)
    if a_rows.shape != b_rows.shape:
        return False
    s = np.linalg.svd(_orthonormal(a_rows).conj().T @ _orthonormal(b_rows),
                      compute_uv=False)
    return bool(np.all(np.abs(s - 1.0) <= tol))


def same_eigenvectors(vals, a, b, tol=SUBSPACE_TOL):
    """True when the eigenvectors ``a`` and ``b`` (``(n, ...)``, one per
    eigenvalue of ``vals``) span the same subspace in every degenerate set."""
    a, b = np.asarray(a), np.asarray(b)
    return a.shape == b.shape and all(same_span(a[g], b[g], tol)
                                      for g in degenerate_sets(vals))


def same_leading_vector(vals, a, b, vecs_a=None, vecs_b=None, tol=SUBSPACE_TOL):
    """True when the leading vectors ``a`` and ``b`` agree up to their phase
    (same norm, overlap 1). When the leading eigenvalue of ``vals`` (the
    largest real part) is degenerate, each must instead lie in the other
    run's leading span ``vecs_*`` (``(n, ...)``); without those spans only
    the norms are compared."""
    a, b = np.asarray(a, dtype=complex), np.asarray(b, dtype=complex)
    if a.shape != b.shape:
        return False
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if abs(na - nb) > VALUE_TOL * max(1.0, na):
        return False
    vals = np.asarray(vals).reshape(-1)
    lead = int(np.argmax(vals.real))
    group = next(g for g in degenerate_sets(vals) if lead in g)
    if len(group) == 1:
        return na == 0.0 or abs(abs(np.vdot(a, b)) / (na * nb) - 1.0) <= tol
    if vecs_a is None or vecs_b is None:
        return True
    for v, span in ((a, np.asarray(vecs_b)[group]), (b, np.asarray(vecs_a)[group])):
        q = _orthonormal(span)
        v = v.reshape(-1)
        if np.linalg.norm(v - q @ (q.conj().T @ v)) > tol * max(1.0, np.linalg.norm(v)):
            return False
    return True


def _close(x, y):
    if x.dtype.kind in "biu":
        return bool(np.array_equal(x, y))
    if x.size == 0:
        return True
    return bool(np.max(np.abs(x - y)) <= VALUE_TOL * max(1.0, float(np.max(np.abs(x)))))


def numeric_mismatches(a, b):
    """The keys ``(archive, member)`` of the numeric members of two runs
    (``{(archive, member): array}``) that do not agree; a key present in
    only one run is reported too. ``eigenvalues``, ``eigenvalue`` and the
    returned ``("return", "lambda")`` are eigenvalues; ``eigenvectors`` and
    the ``gap`` of a gap archive (one that records ``eigenvalues``) are
    compared through their spans."""
    keys = set(a) | set(b)
    bad = sorted(k for k in keys if k not in a or k not in b)
    for k in sorted(keys - set(bad)):
        x, y = np.asarray(a[k]), np.asarray(b[k])
        if x.dtype != y.dtype or x.shape != y.shape:
            bad.append(k)
            continue
        fn, name = k
        vals = a.get((fn, "eigenvalues"))
        if name in ("eigenvalues", "eigenvalue") or k == ("return", "lambda"):
            ok = eigenvalues_close(x, y)
        elif name == "eigenvectors" and vals is not None:
            ok = same_eigenvectors(vals, x, y)
        elif name == "gap" and vals is not None:
            ok = same_leading_vector(vals, x, y, a.get((fn, "eigenvectors")),
                                     b.get((fn, "eigenvectors")))
        else:
            ok = _close(x, y)
        if not ok:
            bad.append(k)
    return bad
