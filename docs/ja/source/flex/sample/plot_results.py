#!/usr/bin/env python3
"""Generate figures for the FLEX solver tutorial.

Run from either the 1orb or 2orb directory after completing the FLEX
calculation:
    $ cd 1orb
    $ hwave input.toml
    $ python ../plot_results.py

Or run from the sample_flex directory to process both:
    $ python plot_results.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_npz(filename):
    """Load and return numpy npz archive."""
    return np.load(filename, allow_pickle=True)


def plot_spin_susceptibility(chi_s, kvec, wavenum_table, shape, norb,
                              title, filename):
    """Plot static spin susceptibility chi_s(q, iv_m=0) as 2D color map.

    Parameters
    ----------
    chi_s : ndarray
        Spin susceptibility, shape (nmat, nvol, nd, nd).
    kvec : ndarray
        k-point vectors.
    wavenum_table : ndarray
        Wavenum table for mapping.
    shape : tuple
        (nx, ny, nz).
    norb : int
        Number of orbitals.
    title : str
        Figure title.
    filename : str
        Output filename.
    """
    nx, ny = shape[0], shape[1]
    nmat = chi_s.shape[0]
    nd = chi_s.shape[-1]

    # Static susceptibility: ivm = 0 is at index nmat//2
    # For bosonic Matsubara, m=0 corresponds to the central index
    iv0 = nmat // 2

    # Sum over orbital diagonal components for total susceptibility
    chi_s_static = np.zeros(nx * ny, dtype=complex)
    for a in range(nd):
        chi_s_static += chi_s[iv0, :, a, a]

    chi_s_real = chi_s_static.real

    # Build q-point grid
    qx = np.array([2.0 * np.pi * i / nx for i in range(nx)])
    qy = np.array([2.0 * np.pi * j / ny for j in range(ny)])
    QX, QY = np.meshgrid(qx, qy, indexing='ij')

    chi_map = chi_s_real.reshape(nx, ny)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.pcolormesh(QX / np.pi, QY / np.pi, chi_map,
                        shading="auto", cmap="hot_r")
    ax.set_xlabel(r"$q_x / \pi$")
    ax.set_ylabel(r"$q_y / \pi$")
    ax.set_title(title)
    ax.set_aspect("equal")
    plt.colorbar(im, ax=ax, label=r"$\chi_s(\mathbf{q}, i\nu_0)$")
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {filename}")


def plot_self_energy(sigma, shape, norb, title, filename):
    """Plot imaginary part of self-energy on the Fermi surface.

    Parameters
    ----------
    sigma : ndarray
        Self-energy, shape (nblock, nmat, nvol, nd, nd).
    shape : tuple
        (nx, ny, nz).
    norb : int
        Number of orbitals.
    title : str
        Figure title.
    filename : str
        Output filename.
    """
    nx, ny = shape[0], shape[1]
    nblock = sigma.shape[0]
    nmat = sigma.shape[1]
    nd = sigma.shape[-1]

    # Lowest positive Matsubara frequency (iwn closest to 0 from above)
    iw0 = nmat // 2

    # Sum over spin blocks and orbital diagonal
    sigma_iw0 = np.zeros(nx * ny, dtype=complex)
    for g in range(nblock):
        for a in range(nd):
            sigma_iw0 += sigma[g, iw0, :, a, a]

    im_sigma = sigma_iw0.imag / nblock  # Average over blocks

    qx = np.array([2.0 * np.pi * i / nx for i in range(nx)])
    qy = np.array([2.0 * np.pi * j / ny for j in range(ny)])
    QX, QY = np.meshgrid(qx, qy, indexing='ij')

    sigma_map = im_sigma.reshape(nx, ny)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.pcolormesh(QX / np.pi, QY / np.pi, sigma_map,
                        shading="auto", cmap="RdBu_r")
    ax.set_xlabel(r"$k_x / \pi$")
    ax.set_ylabel(r"$k_y / \pi$")
    ax.set_title(title)
    ax.set_aspect("equal")
    plt.colorbar(im, ax=ax,
                  label=r"$\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_0)$")
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {filename}")


def plot_convergence(logfile, title, filename):
    """Plot convergence history from the log file.

    Parameters
    ----------
    logfile : str
        Path to log output (captured from hwave).
    title : str
        Figure title.
    filename : str
        Output filename.
    """
    iterations = []
    diffs = []
    with open(logfile) as f:
        for line in f:
            if "convergence:" in line:
                parts = line.strip().split("=")
                diff = float(parts[-1])
                iterations.append(len(iterations) + 1)
                diffs.append(diff)

    if not iterations:
        print(f"  No convergence data found in {logfile}")
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(iterations, diffs, "o-", markersize=3, color="C0")
    ax.axhline(y=1e-6, color="red", linestyle="--", alpha=0.7,
               label=r"EPS = $10^{-6}$")
    ax.set_xlabel("SCF iteration")
    ax.set_ylabel(r"$|\Delta\Sigma| / |\Sigma|$")
    ax.set_title(title)
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {filename}")


def plot_matsubara_sigma(sigma, shape, norb, title, filename):
    """Plot self-energy as a function of Matsubara frequency for selected k-points.

    Parameters
    ----------
    sigma : ndarray
        Self-energy, shape (nblock, nmat, nvol, nd, nd).
    shape : tuple
        (nx, ny, nz).
    norb : int
        Number of orbitals.
    title : str
        Figure title.
    filename : str
        Output filename.
    """
    nx, ny = shape[0], shape[1]
    nblock = sigma.shape[0]
    nmat = sigma.shape[1]
    nd = sigma.shape[-1]

    # Take first spin block, orbital (0,0)
    sigma_block = sigma[0, :, :, 0, 0]  # (nmat, nvol)

    # Select representative k-points
    # (0,0): Gamma, (nx//2, ny//2): M=(pi,pi), (nx//2, 0): X=(pi,0)
    k_gamma = 0
    k_M = (nx // 2) * ny + ny // 2
    k_X = (nx // 2) * ny

    # Positive Matsubara frequencies only
    n_pos = nmat // 2
    wn = np.arange(n_pos) * 2 + 1  # (2n+1) in units of pi*T

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    for k_idx, label, color in [(k_gamma, r"$\Gamma$", "C0"),
                                  (k_M, r"$M=(\pi,\pi)$", "C1"),
                                  (k_X, r"$X=(\pi,0)$", "C2")]:
        im_sigma = sigma_block[n_pos:, k_idx].imag
        re_sigma = sigma_block[n_pos:, k_idx].real
        ax1.plot(wn[:20], im_sigma[:20], "o-", markersize=3,
                  label=label, color=color)
        ax2.plot(wn[:20], re_sigma[:20], "o-", markersize=3,
                  label=label, color=color)

    ax1.set_xlabel(r"$\omega_n / (\pi T)$")
    ax1.set_ylabel(r"$\mathrm{Im}\,\Sigma(\mathbf{k}, i\omega_n)$")
    ax1.set_title("Imaginary part")
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel(r"$\omega_n / (\pi T)$")
    ax2.set_ylabel(r"$\mathrm{Re}\,\Sigma(\mathbf{k}, i\omega_n)$")
    ax2.set_title("Real part")
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)

    fig.suptitle(title, fontsize=12)
    plt.tight_layout()
    plt.savefig(filename, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved {filename}")


def process_directory(dirpath):
    """Process a single FLEX output directory."""
    output_dir = os.path.join(dirpath, "output")
    if not os.path.isdir(output_dir):
        print(f"Skipping {dirpath}: no output directory found")
        return

    dirname = os.path.basename(dirpath)
    print(f"\nProcessing {dirname}...")

    # Load data
    chi_s_data = load_npz(os.path.join(output_dir, "chiq_s.npz"))
    chi_s = chi_s_data["chiq_s"]

    sigma_data = load_npz(os.path.join(output_dir, "sigma.npz"))
    sigma = sigma_data["sigma"]

    # Determine shape from chi_s
    nmat = chi_s.shape[0]
    nvol = chi_s.shape[1]
    nd = chi_s.shape[-1]

    # Read shape from input.toml
    import tomli
    with open(os.path.join(dirpath, "input.toml"), "rb") as f:
        config = tomli.load(f)
    cell_shape = config["mode"]["param"]["CellShape"]
    nx, ny, nz = cell_shape
    norb = nd // 2  # nd = norb * ns (ns=2)

    label = f"{norb}-orbital" if norb > 1 else "1-orbital"
    shape = (nx, ny, nz)

    # Fig 1: Spin susceptibility
    plot_spin_susceptibility(
        chi_s, None, None, shape, norb,
        rf"Spin susceptibility $\chi_s(\mathbf{{q}})$ ({label})",
        os.path.join(dirpath, "chi_s.png")
    )

    # Fig 2: Self-energy in k-space
    plot_self_energy(
        sigma, shape, norb,
        rf"$\mathrm{{Im}}\,\Sigma(\mathbf{{k}}, i\omega_0)$ ({label})",
        os.path.join(dirpath, "sigma_kspace.png")
    )

    # Fig 3: Self-energy vs Matsubara frequency
    plot_matsubara_sigma(
        sigma, shape, norb,
        rf"Self-energy vs $\omega_n$ ({label})",
        os.path.join(dirpath, "sigma_matsubara.png")
    )

    # Print summary
    # Find peak of chi_s
    iv0 = nmat // 2
    chi_s_diag = sum(chi_s[iv0, :, a, a].real for a in range(nd))
    q_peak = np.argmax(chi_s_diag)
    qx_peak = (q_peak // ny) * 2.0 * np.pi / nx
    qy_peak = (q_peak % ny) * 2.0 * np.pi / ny
    print(f"  chi_s peak at q = ({qx_peak/np.pi:.2f}pi, {qy_peak/np.pi:.2f}pi), "
          f"value = {chi_s_diag[q_peak]:.4f}")

    # Self-energy at lowest Matsubara frequency
    nblock = sigma.shape[0]
    iw0 = sigma.shape[1] // 2
    im_sigma_avg = np.mean([sigma[g, iw0, :, 0, 0].imag for g in range(nblock)])
    print(f"  <Im Sigma(k, iw0)> = {im_sigma_avg:.6f}")


def main():
    # Determine which directories to process
    if os.path.isfile("input.toml"):
        # Running from within a sample directory
        process_directory(".")
    else:
        # Running from sample_flex directory
        for subdir in ["1orb", "2orb"]:
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), subdir)
            if os.path.isdir(path):
                process_directory(path)


if __name__ == "__main__":
    main()
