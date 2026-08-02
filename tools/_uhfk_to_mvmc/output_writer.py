"""zqp_orbital_uhfk.dat writer and Fij → orbital param aggregation.

The output text format matches mVMC's
Child_OutputOptData (mVMC-1.4.0/src/ComplexUHF/output.c:43-52) which
mVMC reads via fscanf("%d %lf %lf") for the InOrbital /
InOrbitalAntiParallel key (readdef.c:1714-1726). Aggregation follows
OutputAntiParallel_2 (output.c:245-279): param[idx] += F[i, j] * sign,
then divide by count[idx], then adds a small uniform noise to lift the
rank deficiency that mVMC's Pfaffian Slater evaluation cannot handle
when F is built from a single (k, -k) shell.

The noise step mirrors ComplexUHF's own `OutputAntiParallel_2`
(`output.c:274`), which adds `genrand_real2() * pow(10, -eps_int_slater)`
to each averaged ParamOrbital entry for exactly the same reason. The
bridge's default amplitude ``1e-8`` (see ``aggregate_orbital_params``
default and CLI ``--epsilon-noise`` help) sits in the stable plateau
``epsilon_noise <= 1e-7`` where the case_apbc Ne=4 (rank-2) E2E
reproduces UHF energy within ~1 VMC stderr at NVMCSample=10000
without significantly perturbing the physics. See
docs/en/source/uhfk/tools/uhfk_to_mvmc.rst.
"""
from __future__ import annotations

import numpy as np


def aggregate_orbital_params(
    F, mapping, n_orbital_idx, *,
    epsilon_noise=1.0e-8, complex_type=1, rng=None,
):
    """Average F[i, j] * sign per orbital idx class, then lift rank with
    uniform noise of amplitude ``epsilon_noise``.

    Parameters
    ----------
    F : (Ns, Ns) complex
        Physical-basis Fij.
    mapping : dict[(int, int)] = (idx, sign)
        From orbitalidx_reader.parse_orbitalidx_def.
    n_orbital_idx : int
        Total parameter count.
    epsilon_noise : float, optional
        Amplitude of uniform noise added to each param (default 1e-8).
        Set to 0 to disable.

        mVMC's Pfaffian Slater evaluation (slater.c) becomes
        ill-conditioned when F is exactly rank-deficient: a single
        (k, -k) shell occupied at low filling gives rank min(Ne_up, Ne_down)
        which is much smaller than ``Ns``. mVMC then samples states
        with the wrong density (e.g. <n_{0, up}> = 0 instead of
        ``Ne_up / Ns``) and the energy is wrong by integer multiples of
        ``U / 2``. ComplexUHF works around the same singularity by
        injecting noise at the same step (see
        ``mVMC-1.4.0/src/ComplexUHF/output.c:274``).

        Empirically the residual mVMC-vs-UHF bias plateaus for
        ``epsilon_noise <= 1e-7`` (case_apbc Ne=4 L=8 gives
        ``<H>`` within ~1 VMC stderr at ``NVMCSample = 10000`` of UHF),
        so the default 1e-8 sits well inside the stable region.
    complex_type : int, optional
        ``1`` (default) for complex variational parameters (noise is
        added to both real and imag parts). ``0`` for real-valued
        variational parameters (noise is added to the real part only,
        imag stays zero), matching ``orbitalidx.def``'s ``ComplexType``
        header.
    rng : numpy.random.Generator, optional
        Source of randomness for the noise. Defaults to a fresh
        ``np.random.default_rng(7919)`` so output is reproducible.

    Returns
    -------
    params : (n_orbital_idx,) complex128
        Averaged f_{idx} values plus uniform noise (with idx == -1
        entries skipped and left at 0).
    """
    F = np.asarray(F, dtype=np.complex128)
    params = np.zeros(n_orbital_idx, dtype=np.complex128)
    counts = np.zeros(n_orbital_idx, dtype=np.int64)
    for (i, j), (idx, sign) in mapping.items():
        if idx < 0:
            continue
        params[idx] += F[i, j] * sign
        counts[idx] += 1
    nonzero = counts > 0
    params[nonzero] = params[nonzero] / counts[nonzero]

    if epsilon_noise > 0:
        if rng is None:
            rng = np.random.default_rng(7919)
        noise_re = rng.uniform(-0.5, 0.5, size=n_orbital_idx) * epsilon_noise
        if complex_type == 1:
            noise_im = rng.uniform(-0.5, 0.5, size=n_orbital_idx) * epsilon_noise
            params[nonzero] = (
                params[nonzero] + noise_re[nonzero] + 1j * noise_im[nonzero]
            )
        else:
            # Real-valued variational parameters: lift only the real part.
            params[nonzero] = params[nonzero] + noise_re[nonzero]
    return params


def write_zqp_orbital(out_path, params):
    """Write the InOrbital file text in mVMC ChildOutputOptData format.

    Header (5 lines) matches the mVMC `ReadInputParameters` reader
    (readdef.c:1467-1477), which reads `NOrbitalIdx <N>` from LINE 2 via
    `sscanf(line2, "%s %d", ctmp, &idx)`.

        ======================
        NOrbitalIdx <N>
        ======================
        == i_j_OrbitalIdx ===
        ======================
        <idx> <real> <imag>
        ...
    """
    params = np.asarray(params, dtype=np.complex128)
    n = int(len(params))
    with open(out_path, "w") as fp:
        fp.write("======================\n")
        fp.write(f"NOrbitalIdx  {n}\n")
        fp.write("======================\n")
        fp.write("== i_j_OrbitalIdx ===\n")
        fp.write("======================\n")
        for i in range(n):
            fp.write(
                "{:d} {: .18e} {: .18e}\n".format(
                    i, float(params[i].real), float(params[i].imag)
                )
            )
