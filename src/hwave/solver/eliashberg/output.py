"""Writing the results of hwave_sc and the spectral-shift note.

Moved verbatim from hwave.sc (dev/sc-split-plan.md, step 3); the shell
re-exports every name here under its old spelling.
"""
from __future__ import annotations

import logging
import os
import numpy as np

logger = logging.getLogger("hwave_sc")


def _shifted_eigenvalue_note(solver_mode, spectral_shift, converged, n_iter,
                             iter_info, kernel_label="kernel"):
    """The ``eigenvalue.dat`` comment describing a power iteration's number.

    Split out of ``calc_eliashberg`` so the static and dynamic entry points --
    and the tests -- share one wording, and so the "this is NOT an eigenvalue"
    case can never be silently dropped from the file. ``converged`` is the
    ``eigen_converged`` verdict (``_solve_leading``'s ``info["converged"]``):
    when it is False the reported value is NaN and the note says why.
    ``spectral_shift`` None means the unshifted iteration
    (``spectral_shift = "none"``).
    """
    info = iter_info or {}
    resid = info.get("rayleigh_residual")
    resid_txt = ("{:.3e}".format(resid) if isinstance(resid, float)
                 and np.isfinite(resid) else str(resid))
    shift_txt = ("spectral_shift={!r}".format(spectral_shift)
                 if spectral_shift is not None else "no spectral_shift")
    if not converged:
        if info.get("two_cycle"):
            why = ("the iterate settled into a sign-alternating two-cycle "
                   "after {} steps -- the dominant (largest-magnitude) "
                   "eigenvalue of the iterated operator is NEGATIVE"
                   .format(n_iter))
            cyc = info.get("cycle_eigenvalue")
            if cyc is not None:
                why += " ({:.8e} on the {})".format(cyc, kernel_label)
        elif info.get("diff_converged"):
            why = ("the iterate difference converged after {} steps but the "
                   "Rayleigh check left a residual ||K v - lambda v||/||v|| = "
                   "{} (tolerance {})".format(
                       n_iter, resid_txt, info.get("rayleigh_residual_tol")))
        else:
            why = "the iteration did NOT converge after {} steps".format(n_iter)
        return (
            "solver_mode='{}' with {}: {}. eigen_converged = false: NO "
            "eigenvalue is reported (the value below is NaN, NOT an "
            "eigenvalue) and the gap written with it is the last iterate, NOT "
            "an eigenvector. Use solver_mode=\"eigenvalue\", or "
            "spectral_shift=\"auto\" (the default for the iteration), or "
            "raise max_iter.".format(solver_mode, shift_txt, why))
    if spectral_shift is None:
        return (
            "solver_mode='{}' with no spectral_shift: the value below is the "
            "SIGNED Rayleigh eigenvalue <v|K|v>/<v|v> of the converged "
            "iterate (residual {}); the iteration converged after {} steps. "
            "The unshifted iteration follows the largest-MAGNITUDE "
            "eigenvalue.".format(solver_mode, resid_txt, n_iter))
    note = (
        "solver_mode='{}' with spectral_shift={!r}: the power iteration "
        "ran on the shifted {} K + sigma*I and the value below is the "
        "SIGNED eigenvalue of K (the shift has been subtracted back), "
        "VALIDATED by the Rayleigh quotient <v|K|v>/<v|v> with residual "
        "||K v - lambda v||/||v|| = {}; the iteration converged after {} "
        "steps.".format(solver_mode, spectral_shift, kernel_label, resid_txt,
                        n_iter))
    if info.get("shift_sufficient") is False:
        note += (" WARNING: spectral_shift appears INSUFFICIENT -- the "
                 "mode found has lambda + sigma < 0, so it is the "
                 "largest-MAGNITUDE eigenvalue of K + sigma*I, not the "
                 "algebraically largest (physical) one; re-run with a "
                 "larger spectral_shift or \"auto\".")
    return note


def _save_results(output_dir, sigma, eigenvalue, eigenvalues_eig, kx_array, ky_array, kz_array,
                  gap_file="gap.dat", eigenvalue_file="eigenvalue.dat",
                  eigenvalue_match=None, provenance=None,
                  eigenvalue_note=None, eigenvalue_weights=None,
                  gap_comments=None, eigen_converged=None,
                  iteration_converged=None, analysis_note=None,
                  leading_eigenvalue=None):
    """Save gap function and eigenvalue results to files.

    Parameters
    ----------
    output_dir : str
        Output directory path.
    sigma : ndarray or None
        Gap function from iteration, shape (norb, norb, Nx, Ny, Nz).
    eigenvalue : float or None
        Leading eigenvalue from iteration.
    eigenvalues_eig : ndarray or None
        Eigenvalues from eigenvalue analysis.
    kx_array, ky_array, kz_array : ndarray
        k-point arrays.
    gap_file : str
        Output filename for gap function.
    eigenvalue_file : str
        Output filename for eigenvalues.
    provenance : dict, optional
        Extra approximation-level metadata (currently the bond-channel record,
        spec S5). Written as ``# key = value`` COMMENT lines at the top of the
        eigenvalue file, so this is a strictly ADDITIVE schema change: existing
        readers (``np.loadtxt``, ``tsweep.parse_leading_eig``) skip ``#`` lines
        and every existing key/column is untouched. When None (the default) the
        file is byte-for-byte what it always was.
    eigenvalue_note : str, optional
        Free text describing WHAT the single ``eigenvalue`` number is, written
        as ``#`` comment line(s) between the ``# Iteration eigenvalue`` header
        and the value itself. Used to state that an iteration-mode value is an
        unsigned iterate norm (and whether it converged), so it can never be
        silently confused with the signed ``lambda_rayleigh`` sitting in the
        provenance block above it. Additive in the same sense as
        ``provenance``: comment lines only, no numeric row changes.
    eigenvalue_weights : ndarray, optional
        Spin-orbital mode only: ``(n_eig, 4)`` spin weights
        ``(psi, dx, dy, dz)`` of each eigenvector
        (:func:`hwave.solver.so_pairing.so_spin_decomposition`), appended as
        four columns after the ``match`` column. None (the default) leaves the
        eigenvalue rows unchanged.
    gap_comments : list of str, optional
        Spin-orbital mode only: ``#`` comment lines written before the column
        header of the gap file. None (the default) leaves the gap file
        unchanged.
    eigen_converged : bool, optional
        Written as ``# eigen_converged = <true|false>`` -- right after
        ``# Iteration eigenvalue`` when the iteration ran, else at the top.
        ``false``: the power iteration did not converge, the iteration value
        is NaN and the gap file holds the last iterate, not an eigenvector.
        With ``iteration_converged`` given (``solver_mode = "both"``) it
        refers to the eigenvalue analysis and is written right before
        ``# Eigenvalue analysis``.
    iteration_converged : bool, optional
        ``solver_mode = "both"`` only: the power iteration's own verdict,
        written as ``# iteration_converged = <true|false>`` after
        ``# Iteration eigenvalue``.
    analysis_note : str, optional
        Where the eigenvalue analysis' reported leader comes from (the
        channel's parity-sector solve), written as ``#`` lines right after
        the ``# eigen_converged`` line.
    leading_eigenvalue : complex, optional
        The eigenvalue analysis' reported leader (NaN when none), written as
        the machine-readable line ``# leading_eigenvalue = <Re> <Im>`` right
        before ``# Eigenvalue analysis``; it equals row 0 of the table.
    """
    os.makedirs(output_dir, exist_ok=True)

    if sigma is not None:
        norb = sigma.shape[0]
        Nx, Ny, Nz = len(kx_array), len(ky_array), len(kz_array)
        filepath = os.path.join(output_dir, gap_file)
        logger.info("Saving gap function to {}".format(filepath))

        with open(filepath, "w") as fw:
            for line in (gap_comments or ()):
                fw.write("# {}\n".format(line))
            # Header
            header_parts = ["# kx", "ky", "kz"]
            for i in range(norb):
                for j in range(norb):
                    header_parts.append("Re(sigma_{}{})".format(i, j))
                    header_parts.append("Im(sigma_{}{})".format(i, j))
            fw.write(" ".join(header_parts) + "\n")

            for ix in range(Nx):
                kx = kx_array[ix]
                if kx > np.pi:
                    kx -= 2.0 * np.pi
                for iy in range(Ny):
                    ky = ky_array[iy]
                    if ky > np.pi:
                        ky -= 2.0 * np.pi
                    for iz in range(Nz):
                        kz = kz_array[iz]
                        if kz > np.pi:
                            kz -= 2.0 * np.pi
                        parts = ["{:.8f}".format(kx),
                                 "{:.8f}".format(ky),
                                 "{:.8f}".format(kz)]
                        for i in range(norb):
                            for j in range(norb):
                                val = sigma[i, j, ix, iy, iz]
                                parts.append("{:.8e}".format(val.real))
                                parts.append("{:.8e}".format(val.imag))
                        fw.write(" ".join(parts) + "\n")

    if eigenvalue is not None or eigenvalues_eig is not None:
        filepath = os.path.join(output_dir, eigenvalue_file)
        logger.info("Saving eigenvalues to {}".format(filepath))
        with open(filepath, "w") as fw:
            if provenance:
                for key in sorted(provenance):
                    fw.write("# {} = {}\n".format(key, provenance[key]))
            flag_line = (None if eigen_converged is None else
                         "# eigen_converged = {}\n".format(
                             str(bool(eigen_converged)).lower()))
            if eigenvalue is not None:
                fw.write("# Iteration eigenvalue\n")
                if iteration_converged is not None:
                    fw.write("# iteration_converged = {}\n".format(
                        str(bool(iteration_converged)).lower()))
                elif flag_line:
                    fw.write(flag_line)
                    if analysis_note:
                        fw.write("# {}\n".format(analysis_note))
                if eigenvalue_note:
                    for line in str(eigenvalue_note).splitlines():
                        fw.write("# {}\n".format(line))
                fw.write("{:.8e}\n".format(eigenvalue))
            if eigenvalues_eig is not None:
                if flag_line and (eigenvalue is None
                                  or iteration_converged is not None):
                    fw.write(flag_line)
                    if analysis_note:
                        fw.write("# {}\n".format(analysis_note))
                if leading_eigenvalue is not None:
                    lv = complex(leading_eigenvalue)
                    fw.write("# leading_eigenvalue = {:.8e} {:.8e}\n".format(
                        lv.real, lv.imag))
                fw.write("# Eigenvalue analysis\n")
                # The optional trailing `match` column is 1 when the eigenvector
                # lies (dominantly) in the channel's parity sector (even for
                # singlet, odd for triplet) and 0 when it lies in the opposite
                # sector. The leading index/Re/Im/|ev| columns are unchanged, so
                # position-based parsers (e.g. reading column 1 as Re) are
                # unaffected; only the column count grows from 4 to 5.
                if eigenvalue_weights is not None:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|  match(1=fermion-antisymmetric)  "
                             "w_psi  w_dx  w_dy  w_dz\n")
                elif eigenvalue_match is not None:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|  match(1=channel-parity)\n")
                else:
                    fw.write("# index  Re(eigenvalue)  Im(eigenvalue)  "
                             "|eigenvalue|\n")
                for i, ev in enumerate(eigenvalues_eig):
                    row = "{:4d} {:15.8e} {:15.8e} {:15.8e}".format(
                        i, ev.real, ev.imag, abs(ev))
                    if eigenvalue_match is not None:
                        row += " {:d}".format(int(bool(eigenvalue_match[i])))
                    if eigenvalue_weights is not None:
                        row += "".join(" {:15.8e}".format(w)
                                       for w in eigenvalue_weights[i])
                    fw.write(row + "\n")
