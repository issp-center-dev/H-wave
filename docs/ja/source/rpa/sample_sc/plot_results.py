#!/usr/bin/env python3
"""Generate figures for the Eliashberg solver tutorial.

Run from the sample_sc directory after completing the RPA and Eliashberg
calculations:
    $ hwave input.toml
    $ hwave_sc input.toml
    $ python plot_results.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    """Load eigenvalues from output file."""
    ev_iter = None
    ev_list = []
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
            elif not in_eigenvalue_section and line and not line.startswith("#"):
                ev_iter = float(line.strip())
    return ev_iter, ev_list


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


def plot_eigenvalue_spectrum(ev_singlet, ev_triplet,
                             ev_iter_singlet, ev_iter_triplet, filename):
    """Plot positive eigenvalue spectrum comparison.

    Only positive eigenvalues are plotted, since negative eigenvalues
    do not indicate SC instability. The SC instability criterion is
    lambda = 1.
    """
    fig, ax = plt.subplots(figsize=(7, 4.5))

    # Filter positive eigenvalues only
    ev_s_pos = sorted([e for e in ev_singlet if e > 0], reverse=True)
    ev_t_pos = sorted([e for e in ev_triplet if e > 0], reverse=True)

    ax.scatter(range(len(ev_s_pos)), ev_s_pos,
               marker="o", s=80, color="C0", label="Singlet", zorder=3)
    ax.scatter(range(len(ev_t_pos)), ev_t_pos,
               marker="s", s=80, color="C1", label="Triplet", zorder=3)

    ax.axhline(y=1.0, color="red", linestyle="--", alpha=0.7,
               label=r"$\lambda = 1$ (SC instability)")

    ax.set_xlabel("Eigenvalue index")
    ax.set_ylabel(r"$\lambda$")
    ax.set_title("Eigenvalue spectrum of Eliashberg equation")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {filename}")


def main():
    # Load singlet results
    kx_s, ky_s, sigma_s = load_gap("output/gap.dat")
    ev_iter_s, ev_list_s = load_eigenvalues("output/eigenvalue.dat")

    # Load triplet results
    kx_t, ky_t, sigma_t = load_gap("output/gap_triplet.dat")
    ev_iter_t, ev_list_t = load_eigenvalues("output/eigenvalue_triplet.dat")

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
    plot_eigenvalue_spectrum(ev_list_s, ev_list_t,
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
