#!/usr/bin/env python3
"""Generate figures for the iron pnictide 2-orbital FLEX tutorial.

Run after completing the FLEX calculation:
    $ hwave input.toml
    $ python plot_results.py
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    nx, ny = 16, 16
    nvol = nx * ny

    # --- Load data ---
    chi_s = np.load("output/chiq_s.npz", allow_pickle=True)["chiq_s"]
    chi_c = np.load("output/chiq_c.npz", allow_pickle=True)["chiq_c"]
    sigma = np.load("output/sigma.npz", allow_pickle=True)["sigma"]

    nmat_chi = chi_s.shape[0]
    nmat_sig = sigma.shape[1]
    nd = chi_s.shape[-1]
    norb = sigma.shape[-1]  # nd_block (=norb for spin-free)

    iv0 = nmat_chi // 2  # static (bosonic omega=0)
    iw0 = nmat_sig // 2  # lowest positive Matsubara

    # q-grid
    qx = np.arange(nx) * 2.0 * np.pi / nx
    qy = np.arange(ny) * 2.0 * np.pi / ny
    QX, QY = np.meshgrid(qx, qy, indexing='ij')

    # --- Fig 1: Spin susceptibility chi_s(q) ---
    # The general (full-vertex) scheme stores the four-leg susceptibility
    # chi[iv, q, a, a', b, b'].  The physical (orbital-traced) susceptibility
    # is the sum over the orbital-diagonal legs, sum_{a,b} chi[iv, q, a,a, b,b].
    if chi_s.ndim == 6:
        chi_s_total = sum(chi_s[iv0, :, a, a, b, b].real
                          for a in range(nd) for b in range(nd))
        chi_c_total = sum(chi_c[iv0, :, a, a, b, b].real
                          for a in range(nd) for b in range(nd))
    else:
        chi_s_total = sum(chi_s[iv0, :, a, a].real for a in range(nd))
        chi_c_total = sum(chi_c[iv0, :, a, a].real for a in range(nd))

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    im0 = axes[0].pcolormesh(QX / np.pi, QY / np.pi,
                              chi_s_total.reshape(nx, ny),
                              shading="auto", cmap="hot_r")
    axes[0].set_xlabel(r"$q_x / \pi$")
    axes[0].set_ylabel(r"$q_y / \pi$")
    axes[0].set_title(r"Spin $\chi_s(\mathbf{q})$")
    axes[0].set_aspect("equal")
    plt.colorbar(im0, ax=axes[0])

    im1 = axes[1].pcolormesh(QX / np.pi, QY / np.pi,
                              chi_c_total.reshape(nx, ny),
                              shading="auto", cmap="hot_r")
    axes[1].set_xlabel(r"$q_x / \pi$")
    axes[1].set_ylabel(r"$q_y / \pi$")
    axes[1].set_title(r"Charge $\chi_c(\mathbf{q})$")
    axes[1].set_aspect("equal")
    plt.colorbar(im1, ax=axes[1])

    plt.tight_layout()
    plt.savefig("chi_spin_charge.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved chi_spin_charge.png")

    # --- Fig 2: Orbital-resolved self-energy Im Sigma(k, iw0) ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    orb_titles = [r"$\mathrm{Im}\,\Sigma_{d_{xz}}(\mathbf{k}, i\omega_0)$",
                   r"$\mathrm{Im}\,\Sigma_{d_{yz}}(\mathbf{k}, i\omega_0)$"]
    for idx in range(2):
        im_sig = sigma[0, iw0, :, idx, idx].imag.reshape(nx, ny)
        im = axes[idx].pcolormesh(QX / np.pi, QY / np.pi, im_sig,
                                   shading="auto", cmap="RdBu_r")
        axes[idx].set_xlabel(r"$k_x / \pi$")
        axes[idx].set_ylabel(r"$k_y / \pi$")
        axes[idx].set_title(orb_titles[idx])
        axes[idx].set_aspect("equal")
        plt.colorbar(im, ax=axes[idx])

    plt.tight_layout()
    plt.savefig("sigma_orbital.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved sigma_orbital.png")

    # --- Fig 3: Self-energy vs Matsubara frequency ---
    T = 0.1
    n_pos = nmat_sig // 2
    wn = (np.arange(n_pos) * 2 + 1) * np.pi * T

    # k-points: Gamma, M=(pi,0), X=(pi,pi)
    k_gamma = 0
    k_M = (nx // 2) * ny
    k_X = (nx // 2) * ny + ny // 2

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    n_show = 20

    orb_names = [r"$d_{xz}$", r"$d_{yz}$"]
    ylabels = [r"$\mathrm{Im}\,\Sigma_{d_{xz}}(i\omega_n)$",
               r"$\mathrm{Im}\,\Sigma_{d_{yz}}(i\omega_n)$"]
    for idx in range(2):
        for k_idx, k_label, color in [(k_gamma, r"$\Gamma$", "C0"),
                                        (k_M, r"$M=(\pi,0)$", "C1"),
                                        (k_X, r"$X=(\pi,\pi)$", "C2")]:
            im_sig = sigma[0, n_pos:n_pos + n_show, k_idx, idx, idx].imag
            axes[idx].plot(wn[:n_show], im_sig, "o-", markersize=3,
                           label=k_label, color=color)

        axes[idx].set_xlabel(r"$\omega_n$")
        axes[idx].set_ylabel(ylabels[idx])
        axes[idx].set_title("Orbital " + orb_names[idx])
        axes[idx].legend(fontsize=8)
        axes[idx].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig("sigma_matsubara_orbital.png", dpi=150, bbox_inches="tight")
    plt.close()
    print("Saved sigma_matsubara_orbital.png")

    # --- Summary ---
    q_peak = np.argmax(chi_s_total)
    qx_p = (q_peak // ny) * 2.0 * np.pi / nx
    qy_p = (q_peak % ny) * 2.0 * np.pi / ny
    print(f"\nSpin susceptibility peak: q=({qx_p/np.pi:.2f}pi, "
          f"{qy_p/np.pi:.2f}pi), value={chi_s_total[q_peak]:.4f}")
    print(f"  chi_s(pi,0) = {chi_s_total[(nx//2)*ny]:.4f}")
    print(f"  chi_s(0,pi) = {chi_s_total[ny//2]:.4f}")
    print(f"  chi_s(pi,pi) = {chi_s_total[(nx//2)*ny + ny//2]:.4f}")


if __name__ == "__main__":
    main()
