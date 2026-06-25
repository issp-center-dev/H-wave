#!/usr/bin/env python3
"""Generate figures for the Eliashberg solver tutorial.

Run from the sample_sc directory after completing the RPA and Eliashberg
calculations:
    $ hwave input.toml
    $ hwave_sc input.toml
    $ python plot_results.py
"""

import numpy as np

# matplotlib is imported lazily inside the plotting functions so the data-loading
# helpers (e.g. load_eigenvalues) can be imported without a matplotlib install.


def _get_plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def load_gap(filename):
    """Load gap function from output file."""
    data = np.loadtxt(filename, comments="#")
    kx = data[:, 0]
    ky = data[:, 1]
    # sigma_00 (Re, Im), sigma_01 (Re, Im), sigma_10 (Re, Im), sigma_11 (Re, Im)
    sigma = {}
    sigma[(0, 0)] = data[:, 3] + 1j * data[:, 4]
    sigma[(0, 1)] = data[:, 5] + 1j * data[:, 6]
    sigma[(1, 0)] = data[:, 7] + 1j * data[:, 8]
    sigma[(1, 1)] = data[:, 9] + 1j * data[:, 10]
    return kx, ky, sigma


def load_eigenvalues(filename):
    """Load eigenvalues from output file.

    Returns (ev_iter, ev_list, match_list). ``match_list[i]`` is True when the
    eigenvector lies (dominantly) in the channel's parity sector (even for
    singlet, odd for triplet) and False when it lies in the opposite-parity
    sector. For older output files without the ``match`` column, ``match_list``
    is ``None`` (the sector is unknown and must not be inferred).
    """
    ev_iter = None
    ev_list = []
    match_list = []
    n_rows = 0
    n_with_match = 0
    with open(filename) as f:
        in_eigenvalue_section = False
        for line in f:
            line = line.strip()
            if line.startswith("# Iteration eigenvalue"):
                continue
            elif line.startswith("# Eigenvalue analysis"):
                in_eigenvalue_section = True
                continue
            elif line.startswith("# index"):
                continue
            elif in_eigenvalue_section and line:
                parts = line.split()
                ev_list.append(float(parts[1]))
                n_rows += 1
                if len(parts) > 4:
                    n_with_match += 1
                    match_list.append(parts[4] == "1")
                else:
                    match_list.append(None)
            elif not in_eigenvalue_section and line and not line.startswith("#"):
                ev_iter = float(line.strip())
    # Trust the match column only if EVERY eigenvalue row carries it; a missing
    # or partially annotated (e.g. hand-edited) column is treated as legacy so
    # unknown-sector rows are never misclassified as opposite-parity.
    has_match = n_rows > 0 and n_with_match == n_rows
    return ev_iter, ev_list, (match_list if has_match else None)


def plot_gap_kspace(kx, ky, sigma, title, filename, N=32):
    """Plot gap function in k-space as a 2D color map."""
    # Shift to first Brillouin zone [-pi, pi]
    kx = np.where(kx > np.pi, kx - 2 * np.pi, kx)
    ky = np.where(ky > np.pi, ky - 2 * np.pi, ky)

    # Sort for proper plotting
    idx = np.lexsort((ky, kx))
    kx = kx[idx]
    ky = ky[idx]
    for key in sigma:
        sigma[key] = sigma[key][idx]

    KX = kx.reshape(N, N)
    KY = ky.reshape(N, N)

    plt = _get_plt()
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Intra-orbital: sigma_00
    val_00 = np.real(sigma[(0, 0)]).reshape(N, N)
    im0 = axes[0].pcolormesh(KX / np.pi, KY / np.pi, val_00,
                              shading="auto", cmap="RdBu_r")
    axes[0].set_xlabel(r"$k_x / \pi$")
    axes[0].set_ylabel(r"$k_y / \pi$")
    axes[0].set_title(r"$\mathrm{Re}\,\Sigma_{00}(\mathbf{k})$")
    axes[0].set_aspect("equal")
    plt.colorbar(im0, ax=axes[0])

    # Inter-orbital: sigma_01
    val_01 = np.real(sigma[(0, 1)]).reshape(N, N)
    im1 = axes[1].pcolormesh(KX / np.pi, KY / np.pi, val_01,
                              shading="auto", cmap="RdBu_r")
    axes[1].set_xlabel(r"$k_x / \pi$")
    axes[1].set_ylabel(r"$k_y / \pi$")
    axes[1].set_title(r"$\mathrm{Re}\,\Sigma_{01}(\mathbf{k})$")
    axes[1].set_aspect("equal")
    plt.colorbar(im1, ax=axes[1])

    fig.suptitle(title, fontsize=14)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {filename}")


def plot_eigenvalue_spectrum(ev_singlet, match_singlet, ev_triplet, match_triplet,
                             ev_iter_singlet, ev_iter_triplet, filename):
    """Plot the positive eigenvalue spectrum, separating the channel-parity
    sector from the opposite-parity sector.

    Only positive eigenvalues are plotted, since negative eigenvalues do not
    indicate an SC instability. The SC criterion is lambda = 1. For each channel
    the filled markers are eigenvalues whose gap lies in the channel's parity
    sector (even for singlet, odd for triplet); open markers are opposite-parity
    eigenvalues. For this centrosymmetric model parity is an exact kernel
    symmetry, so the opposite-parity modes (e.g. the even-parity modes of the
    triplet kernel that exceed 1) are unphysical for that spin channel by the
    Pauli principle.
    """
    plt = _get_plt()
    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Sector resolution is only possible when the output carries the `match`
    # column. For legacy files (match is None) fall back to the previous
    # unsplit plot so we never infer a parity sector we do not actually know.
    have_sectors = match_singlet is not None and match_triplet is not None

    def split_pos(evs, matches):
        if matches is None:
            return sorted([e for e in evs if e > 0], reverse=True), []
        phys = sorted([e for e, m in zip(evs, matches) if e > 0 and m], reverse=True)
        spur = sorted([e for e, m in zip(evs, matches) if e > 0 and not m], reverse=True)
        return phys, spur

    s_phys, s_spur = split_pos(ev_singlet, match_singlet)
    t_phys, t_spur = split_pos(ev_triplet, match_triplet)

    s_label = "Singlet (even sector)" if have_sectors else "Singlet"
    t_label = "Triplet (odd sector)" if have_sectors else "Triplet"
    ax.scatter(range(len(s_phys)), s_phys, marker="o", s=80, color="C0",
               label=s_label, zorder=3)
    ax.scatter(range(len(t_phys)), t_phys, marker="s", s=80, color="C1",
               label=t_label, zorder=3)
    # Opposite-parity sector: open markers (only when sectors are known).
    if s_spur:
        ax.scatter(range(len(s_spur)), s_spur, marker="o", s=80,
                   facecolors="none", edgecolors="C0", alpha=0.6, zorder=2)
    if t_spur:
        ax.scatter(range(len(t_spur)), t_spur, marker="s", s=80,
                   facecolors="none", edgecolors="C1", alpha=0.6, zorder=2,
                   label="Opposite-parity sector")

    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.7,
               label=r"$\lambda = 1$ (SC instability)")

    ax.set_xlabel("Eigenvalue index")
    ax.set_ylabel(r"$\lambda$")
    ax.set_title("Eigenvalue spectrum of Eliashberg equation")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {filename}")


def main():
    # Load singlet results
    kx_s, ky_s, sigma_s = load_gap("output/gap.dat")
    ev_iter_s, ev_list_s, match_s = load_eigenvalues("output/eigenvalue.dat")

    # Load triplet results
    kx_t, ky_t, sigma_t = load_gap("output/gap_triplet.dat")
    ev_iter_t, ev_list_t, match_t = load_eigenvalues("output/eigenvalue_triplet.dat")

    N = 32

    # Fig 1: Singlet gap function
    plot_gap_kspace(kx_s, ky_s, sigma_s,
                    f"Singlet gap ($\\lambda = {ev_iter_s:.4f}$)",
                    "gap_singlet.png", N)

    # Fig 2: Triplet gap function
    plot_gap_kspace(kx_t, ky_t, sigma_t,
                    f"Triplet gap ($\\lambda = {ev_iter_t:.4f}$)",
                    "gap_triplet.png", N)

    # Fig 3: Eigenvalue spectrum comparison
    plot_eigenvalue_spectrum(ev_list_s, match_s, ev_list_t, match_t,
                             ev_iter_s, ev_iter_t,
                             "eigenvalue_spectrum.png")

    # Print summary
    print(f"\nSinglet: lambda_iter = {ev_iter_s:.6f}")
    print(f"  |sigma_01|/|sigma_00| = "
          f"{np.max(np.abs(sigma_s[(0,1)]))/np.max(np.abs(sigma_s[(0,0)])):.2f}")
    print(f"Triplet: lambda_iter = {ev_iter_t:.6f}")
    print(f"  |sigma_00|/|sigma_01| = "
          f"{np.max(np.abs(sigma_t[(0,0)]))/np.max(np.abs(sigma_t[(0,1)])):.2f}")


if __name__ == "__main__":
    main()
