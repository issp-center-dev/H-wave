"""Parse emitted bridge output into F and project F to
one-body density G under Slater-form assumption.

Two callables, both pinned by the seven-gate contract:

- ``parse_emitted_F(workspace)`` reads
  ``${workspace}/{namelist.def, orbitalidx_general.def, zqp_orbital_uhfk.dat}``
  and reconstructs the ``(2*Nsite, 2*Nsite)`` complex antisymmetric F matrix
  that mVMC's ``InOrbitalGeneral`` reader materializes. Nsite is authoritative
  from ``namelist.def`` (via the ``ModPara`` reference) and cross-checked
  against the site range in ``orbitalidx_general.def``.

- ``pair_product_density_from_F(F, N_pairs, rank_tol=None)`` uses a
  skew-SVD projector returning the one-body density
  ``G[i, j] = <c^dag_i c_j>`` under H-wave / bridge convention
  ``G = conj(A) @ A.T``. Uses ``np.linalg.svd``; top ``2 * N_pairs`` left-
  singular vectors span the occupied subspace.

Both functions are consumed by the G0-writer-check and G2a-emitted-F gates;
they must NOT be substituted with alternate helpers. See
docs/en/source/algorithm/uhfk_to_mvmc.rst for the construction and
docs/en/source/uhfk/tools/uhfk_to_mvmc.rst for the validation gates.
"""
from __future__ import annotations

import os

import numpy as np


class NamelistFormatError(RuntimeError):
    """Raised on malformed ``namelist.def`` / ``modpara.def`` content."""


class ZqpFormatError(RuntimeError):
    """Raised on malformed ``zqp_orbital_uhfk.dat`` content."""


def _parse_namelist_def(path):
    """Return dict[key] = value_path from the mVMC namelist.def.

    mVMC namelist.def lines have the form ``<KEY>   <path>`` (whitespace-
    separated). Lines starting with ``#`` are comments and skipped.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"namelist.def not found: {path}")
    entries = {}
    with open(path) as fp:
        for ln in fp:
            stripped = ln.strip()
            if not stripped or stripped.startswith("#"):
                continue
            toks = stripped.split()
            if len(toks) < 2:
                raise NamelistFormatError(
                    f"namelist.def line does not parse as '<KEY> <path>': "
                    f"{stripped!r}"
                )
            key = toks[0]
            value = toks[1]
            entries[key] = value
    return entries


def _parse_nsite_from_modpara(path):
    """Return Nsite (int) from the mVMC modpara.def.

    modpara.def has a line ``Nsite   <N>`` after a series of header-only
    ``----`` separators; some mVMC variants pad the value with extra
    whitespace. Any line whose first token is exactly ``Nsite`` is matched.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"modpara.def not found: {path}")
    with open(path) as fp:
        for ln in fp:
            stripped = ln.strip()
            if not stripped:
                continue
            toks = stripped.split()
            if toks[0] == "Nsite":
                if len(toks) < 2:
                    raise NamelistFormatError(
                        f"modpara.def Nsite line has no value: {stripped!r}"
                    )
                try:
                    return int(toks[1])
                except ValueError as exc:
                    raise NamelistFormatError(
                        f"modpara.def Nsite value not integer: {stripped!r}"
                    ) from exc
    raise NamelistFormatError(
        f"modpara.def missing 'Nsite <N>' line: {path}"
    )


def _read_nsite(workspace):
    """Resolve Nsite via the namelist.def -> ModPara -> modpara.def chain.

    ``bridge/namelist.def`` MUST contain a ``ModPara <relative-path>`` entry;
    the resolver reads that file (relative to the workspace) and returns its
    ``Nsite`` value. This matches how mVMC's readdef.c populates ``Nsite`` for
    ``InOrbitalGeneral`` consumers.
    """
    namelist_path = os.path.join(workspace, "namelist.def")
    entries = _parse_namelist_def(namelist_path)
    if "ModPara" not in entries:
        raise NamelistFormatError(
            f"namelist.def missing 'ModPara <path>' entry: {namelist_path}"
        )
    modpara_rel = entries["ModPara"]
    if os.path.isabs(modpara_rel):
        modpara_path = modpara_rel
    else:
        modpara_path = os.path.join(workspace, modpara_rel)
    return _parse_nsite_from_modpara(modpara_path)


def _parse_zqp_orbital_uhfk(path):
    """Return dict[class_idx] = complex(re, im) from zqp_orbital_uhfk.dat.

    The file format (matching ``general_output_writer.write_zqp_orbital_general``)
    is five header lines followed by rows ``<idx> <real> <imag>``.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"zqp_orbital_uhfk.dat not found: {path}")
    params = {}
    with open(path) as fp:
        for ln in fp:
            stripped = ln.strip()
            if not stripped:
                continue
            # Skip separator lines and header keyword lines.
            if stripped.startswith("=") or stripped.startswith("N") or \
                    stripped.startswith("i"):
                continue
            toks = stripped.split()
            if len(toks) != 3:
                raise ZqpFormatError(
                    f"zqp_orbital_uhfk.dat body line has {len(toks)} tokens "
                    f"(expected 3): {stripped!r}"
                )
            try:
                idx = int(toks[0])
                re = float(toks[1])
                im = float(toks[2])
            except ValueError as exc:
                raise ZqpFormatError(
                    f"zqp_orbital_uhfk.dat body line has non-numeric tokens: "
                    f"{stripped!r}"
                ) from exc
            params[idx] = complex(re, im)
    if not params:
        raise ZqpFormatError(
            f"zqp_orbital_uhfk.dat has no body rows: {path}"
        )
    return params


def _parse_orbitalidx_general_rows(path):
    """Return list of (i, spn_i, j, spn_j, class_idx, sign) from the
    ``orbitalidx_general.def`` 6-column mapping rows.

    Unlike ``orbitalidx_general_reader.parse_orbitalidx_general_def`` (which
    returns a dict keyed by ``(all_i, all_j)`` and enforces upper-triangle +
    total-row constraints tied to nsite), this helper returns the raw rows so
    parse_emitted_F can apply the normative
    ``all_i = i + spn_i * Nsite`` mapping with an authoritative Nsite from
    the namelist chain.
    """
    if not os.path.isfile(path):
        raise FileNotFoundError(f"orbitalidx_general.def not found: {path}")
    rows = []
    saw_header_marker = False
    with open(path) as fp:
        for ln in fp:
            stripped = ln.strip()
            if not stripped or stripped.startswith("="):
                continue
            toks = stripped.split()
            if toks[0] in ("NOrbitalIdx", "ComplexType"):
                saw_header_marker = True
                continue
            if len(toks) == 6:
                try:
                    i, spn_i, j, spn_j, class_idx, sign = (
                        int(t) for t in toks
                    )
                except ValueError as exc:
                    raise NamelistFormatError(
                        f"orbitalidx_general.def mapping row has non-integer "
                        f"tokens: {stripped!r}"
                    ) from exc
                if sign not in (-1, 1):
                    raise NamelistFormatError(
                        f"orbitalidx_general.def sign column must be +/-1; "
                        f"got {sign} in {stripped!r}"
                    )
                rows.append((i, spn_i, j, spn_j, class_idx, sign))
            elif len(toks) == 2:
                # Optimize-flag lines; not needed for F reconstruction.
                continue
            else:
                raise NamelistFormatError(
                    f"orbitalidx_general.def mapping row has {len(toks)} "
                    f"tokens (expected 6): {stripped!r}"
                )
    if not saw_header_marker:
        raise NamelistFormatError(
            f"orbitalidx_general.def missing NOrbitalIdx/ComplexType header: "
            f"{path}"
        )
    if not rows:
        raise NamelistFormatError(
            f"orbitalidx_general.def has no mapping rows: {path}"
        )
    return rows


def parse_emitted_F(workspace):
    """Return the ``(2*Nsite, 2*Nsite)`` complex128 antisymmetric F matrix
    reconstructed from the mVMC-compatible bridge output under ``workspace``.

    Reads:
    - ``${workspace}/namelist.def`` (with ``ModPara <path>`` referring to a
      modpara-style file that contains ``Nsite <N>``)
    - ``${workspace}/orbitalidx_general.def`` (InOrbitalGeneral, 6-column
      mapping rows ``i spn_i j spn_j class_idx sign``)
    - ``${workspace}/zqp_orbital_uhfk.dat`` (``<idx> <re> <im>`` rows)

    Normative reconstruction (see
    docs/en/source/algorithm/uhfk_to_mvmc.rst):

    ::

        Nsite  = read from bridge/namelist.def
        idx    = parse bridge/orbitalidx_general.def
                 rows are (i, spn_i, j, spn_j, class_idx, sign)
                 index base is 0
        zqp    = parse bridge/zqp_orbital_uhfk.dat
                 rows are (class_idx, re, im)
                 params[class_idx] = complex(re, im)
        F = zeros((2*Nsite, 2*Nsite), dtype=complex128)
        for (i, spn_i, j, spn_j, class_idx, sign) in idx:
            all_i = i + spn_i * Nsite
            all_j = j + spn_j * Nsite
            F[all_i, all_j] = sign * params[class_idx]
            F[all_j, all_i] = -F[all_i, all_j]           # antisymmetry
        assert np.allclose(F, -F.T, atol=1e-12)
        return F

    Contract-pinned choices:

    - Spin-block ordering ``all_i = i + spn_i * Nsite`` (NOT site-major
      ``all_i = 2*i + spn_i``).
    - ``sign * param`` semantics (NOT ``param / sign`` or other variations).
    - Zero-based indexing.
    - Antisymmetric fill (``F[all_j, all_i] = -F[all_i, all_j]``).
    - Nsite from ``namelist.def`` (via ``ModPara``), not inferred from row
      site indices alone; the two are cross-checked and disagreement raises.
    """
    Nsite = _read_nsite(workspace)
    idx_rows = _parse_orbitalidx_general_rows(
        os.path.join(workspace, "orbitalidx_general.def")
    )
    params = _parse_zqp_orbital_uhfk(
        os.path.join(workspace, "zqp_orbital_uhfk.dat")
    )
    # Cross-check: max site index seen in rows should be < Nsite.
    site_max = max(max(r[0], r[2]) for r in idx_rows)
    if site_max + 1 > Nsite:
        raise NamelistFormatError(
            f"orbitalidx_general.def references site {site_max} but "
            f"namelist.def -> modpara.def says Nsite = {Nsite}; "
            f"the two must agree"
        )

    F = np.zeros((2 * Nsite, 2 * Nsite), dtype=np.complex128)
    for (i, spn_i, j, spn_j, class_idx, sign) in idx_rows:
        if class_idx not in params:
            raise ZqpFormatError(
                f"orbitalidx_general.def row references class_idx={class_idx}"
                f" but zqp_orbital_uhfk.dat has no such class"
            )
        all_i = i + spn_i * Nsite
        all_j = j + spn_j * Nsite
        val = sign * params[class_idx]
        F[all_i, all_j] = val
        F[all_j, all_i] = -val
    # Enforce the documented antisymmetry contract.
    if not np.allclose(F, -F.T, atol=1e-12):
        raise NamelistFormatError(
            "parse_emitted_F: reconstructed F failed antisymmetry check "
            "(F != -F.T at atol=1e-12); orbitalidx rows must not overwrite "
            "each other with inconsistent values"
        )
    return F


def pair_product_density_from_F(F, N_pairs, rank_tol=None):
    """Return the one-body density ``G[i, j] = <c^dag_i c_j>`` for the mVMC
    PairProduct wavefunction ``|Psi> = Pf(F) |vac>`` under Slater-form
    assumption.

    F has rank ``2 * N_pairs`` and decomposes as ``F = A Omega A^T`` where
    ``A`` has orthonormal columns and ``Omega = block_diag([[0,I_N],[-I_N,0]])``.

    Convention (matching H-wave / bridge ``compare_against_green_sublattice``
    which uses ``conj(A) @ A.T``):

    ::

        G = conj(A) @ A.T = conj(U_occ) @ U_occ.T

    where ``U_occ = SVD(F).U[:, :2 * N_pairs]``. NOTE: the wrong orientation
    ``U_occ @ U_occ.conj().T`` gives ``G.T`` (equivalently, ``conj(G)`` for
    Hermitian G) — that is the WRONG convention for complex SOC densities
    where off-diagonal cross-spin entries have non-zero imaginary parts.

    Algorithm (see docs/en/source/algorithm/uhfk_to_mvmc.rst):

    ::

        def pair_product_density_from_F(F, N_pairs, rank_tol=None):
            assert np.allclose(F, -F.T, atol=1e-12)
            U, s, _ = np.linalg.svd(F)
            if rank_tol is None:
                rank_tol = 1e-9  # default suitable for zero-noise inputs
            assert np.all(s[:2 * N_pairs] > 1e-10)
            tail_max = float(np.max(s[2*N_pairs:])) if s.shape[0] > 2*N_pairs else 0.0
            assert tail_max < rank_tol
            top = s[:2 * N_pairs].reshape(-1, 2)
            assert np.allclose(top[:, 0], top[:, 1], atol=max(1e-8, 10*rank_tol))
            U_occ = U[:, :2 * N_pairs]
            return np.conj(U_occ) @ U_occ.T

    ``rank_tol`` policy:

    - Zero-noise inputs (snapshot or hand-computed synthetic):
      ``rank_tol = 1e-10``.
    - Shipping inputs with ``--epsilon-noise 1e-8``: ``rank_tol = 1e-6``.
    """
    F = np.asarray(F, dtype=np.complex128)
    if not np.allclose(F, -F.T, atol=1e-12):
        raise ValueError(
            "pair_product_density_from_F: input F is not antisymmetric "
            "at atol=1e-12 (F != -F.T)"
        )
    U, s, _ = np.linalg.svd(F)
    if rank_tol is None:
        rank_tol = 1e-9  # default suitable for zero-noise inputs
    top_svs = s[: 2 * N_pairs]
    if not np.all(top_svs > 1e-10):
        raise ValueError(
            f"pair_product_density_from_F: top {2 * N_pairs} singular values "
            f"must be > 1e-10; got {top_svs}"
        )
    if s.shape[0] > 2 * N_pairs:
        tail_max = float(np.max(s[2 * N_pairs:]))
    else:
        tail_max = 0.0
    if not tail_max < rank_tol:
        raise ValueError(
            f"pair_product_density_from_F: tail singular values must be < "
            f"rank_tol = {rank_tol}; got max = {tail_max}"
        )
    top = top_svs.reshape(-1, 2)
    pair_tol = max(1e-8, 10 * rank_tol)
    if not np.allclose(top[:, 0], top[:, 1], atol=pair_tol):
        raise ValueError(
            f"pair_product_density_from_F: skew-SVD pairing check failed "
            f"(top singular values must appear in duplicate pairs at atol="
            f"{pair_tol}); got {top}"
        )
    U_occ = U[:, : 2 * N_pairs]
    return np.conj(U_occ) @ U_occ.T
