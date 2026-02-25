#!/usr/bin/env python3
"""Benchmark: block matrix optimization for UHFk _diag() and RPA _solve_rpa().

Measures speedup from block-diagonal decomposition by comparing:
  - Single-block (fully coupled): all indices in one block
  - Multi-block (decomposed): independent sub-blocks

Usage:
    python benchmarks/bench_block_optimization.py
"""

import time
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import hwave.solver.uhfk as uhfk_module
import hwave.solver.rpa as rpa_module


def make_uhfk_stub(norb, ns, nvol, ham_trans, block_info, Nconds):
    """Create a minimal UHFk stub for _diag() benchmarking."""
    nd = norb * ns
    stub = object.__new__(uhfk_module.UHFk)
    stub.norb = norb
    stub.ns = ns
    stub.nd = nd
    stub.nvol = nvol
    stub.ham = ham_trans.copy()
    stub.block_info = block_info
    stub.Nconds = Nconds
    stub.T = 0
    stub.shape = (nvol, 1, 1)
    stub._green_list = {}
    return stub


def make_rpa_stub(nvol):
    """Create a minimal RPA stub for _solve_rpa() benchmarking."""
    class LatticeStub:
        pass
    stub = object.__new__(rpa_module.RPA)
    stub.lattice = LatticeStub()
    stub.lattice.nvol = nvol
    return stub


def bench_uhfk_diag(norb_list, nvol, n_repeat=5):
    """Benchmark UHFk _diag() with single vs multiple blocks."""
    ns = 2
    print("=" * 72)
    print("UHFk _diag() benchmark")
    print(f"  nvol={nvol}, ns={ns}, n_repeat={n_repeat}")
    print("-" * 72)
    print(f"{'norb':>6} {'nd':>6} {'nblk':>6} {'single(ms)':>12} {'multi(ms)':>12} {'speedup':>10}")
    print("-" * 72)

    rng = np.random.RandomState(42)

    for norb in norb_list:
        nd = norb * ns

        # Build Hermitian Hamiltonian (spin-diagonal, so 2 natural blocks)
        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(ns):
            a, b = s * norb, (s + 1) * norb
            blk = rng.randn(nvol, norb, norb) + 1j * rng.randn(nvol, norb, norb)
            blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
            ham[:, a:b, a:b] = blk

        # Single block
        stub_single = make_uhfk_stub(norb, ns, nvol, ham,
                                     block_info=[list(range(nd))],
                                     Nconds=[nd])
        times_single = []
        for _ in range(n_repeat):
            stub_single.ham = ham.copy()
            t0 = time.perf_counter()
            stub_single._diag()
            times_single.append(time.perf_counter() - t0)

        # Multi block (2 spin blocks)
        stub_multi = make_uhfk_stub(norb, ns, nvol, ham,
                                    block_info=[list(range(norb)),
                                                list(range(norb, nd))],
                                    Nconds=[norb, norb])
        times_multi = []
        for _ in range(n_repeat):
            stub_multi.ham = ham.copy()
            t0 = time.perf_counter()
            stub_multi._diag()
            times_multi.append(time.perf_counter() - t0)

        t_single = np.median(times_single) * 1000
        t_multi = np.median(times_multi) * 1000
        speedup = t_single / t_multi if t_multi > 0 else float('inf')

        print(f"{norb:>6} {nd:>6} {2:>6} {t_single:>12.3f} {t_multi:>12.3f} {speedup:>10.2f}x")

    print()


def bench_rpa_solve(norb_list, nvol, nmat, n_repeat=3):
    """Benchmark RPA _solve_rpa() with single vs block-diagonal ham."""
    ns = 2
    print("=" * 72)
    print("RPA _solve_rpa() benchmark")
    print(f"  nvol={nvol}, nmat={nmat}, ns={ns}, n_repeat={n_repeat}")
    print("-" * 72)
    print(f"{'norb':>6} {'nd':>6} {'nblk':>6} {'single(ms)':>12} {'multi(ms)':>12} {'speedup':>10}")
    print("-" * 72)

    rng = np.random.RandomState(123)

    for norb in norb_list:
        nd = norb * ns

        # Build block-diagonal chi0q and ham (reduced scheme)
        chi0q = np.zeros((nmat, nvol, nd, nd), dtype=np.complex128)
        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        for s in range(ns):
            a, b = s * norb, (s + 1) * norb
            blk_chi = rng.randn(nmat, nvol, norb, norb) \
                      + 1j * rng.randn(nmat, nvol, norb, norb)
            blk_chi = (blk_chi + blk_chi.conj().transpose(0, 1, 3, 2)) / 2
            chi0q[:, :, a:b, a:b] = blk_chi

            blk_ham = rng.randn(nvol, norb, norb) \
                      + 1j * rng.randn(nvol, norb, norb)
            blk_ham = (blk_ham + blk_ham.conj().transpose(0, 2, 1)) / 2
            ham[:, a:b, a:b] = blk_ham

        # Full matrix (break block structure by adding tiny cross-block)
        ham_full = ham.copy()
        ham_full[:, 0, norb] = 1e-15
        ham_full[:, norb, 0] = 1e-15

        solver = make_rpa_stub(nvol)

        # Single block (full)
        times_single = []
        for _ in range(n_repeat):
            t0 = time.perf_counter()
            solver._solve_rpa(chi0q, ham_full)
            times_single.append(time.perf_counter() - t0)

        # Multi block (auto-detected)
        times_multi = []
        for _ in range(n_repeat):
            t0 = time.perf_counter()
            solver._solve_rpa(chi0q, ham)
            times_multi.append(time.perf_counter() - t0)

        t_single = np.median(times_single) * 1000
        t_multi = np.median(times_multi) * 1000
        speedup = t_single / t_multi if t_multi > 0 else float('inf')

        print(f"{norb:>6} {nd:>6} {2:>6} {t_single:>12.3f} {t_multi:>12.3f} {speedup:>10.2f}x")

    print()


def bench_uhfk_nblocks(nd, nvol, nblocks_list, n_repeat=5):
    """Benchmark UHFk _diag() varying number of blocks for fixed nd."""
    print("=" * 72)
    print(f"UHFk _diag() benchmark: varying number of blocks (nd={nd}, nvol={nvol})")
    print("-" * 72)
    print(f"{'nblocks':>8} {'blk_size':>10} {'time(ms)':>12} {'speedup':>10}")
    print("-" * 72)

    rng = np.random.RandomState(99)

    # Build Hermitian Hamiltonian (block-diagonal with nblocks blocks)
    t_ref = None
    for nblocks in nblocks_list:
        if nd % nblocks != 0:
            continue
        blk_size = nd // nblocks

        ham = np.zeros((nvol, nd, nd), dtype=np.complex128)
        block_info = []
        for b in range(nblocks):
            a0, a1 = b * blk_size, (b + 1) * blk_size
            blk = rng.randn(nvol, blk_size, blk_size) \
                  + 1j * rng.randn(nvol, blk_size, blk_size)
            blk = (blk + blk.conj().transpose(0, 2, 1)) / 2
            ham[:, a0:a1, a0:a1] = blk
            block_info.append(list(range(a0, a1)))

        # We need ns=1 to avoid issues with Green function construction
        stub = make_uhfk_stub(nd, 1, nvol, ham,
                              block_info=block_info,
                              Nconds=[blk_size] * nblocks)

        times = []
        for _ in range(n_repeat):
            stub.ham = ham.copy()
            t0 = time.perf_counter()
            stub._diag()
            times.append(time.perf_counter() - t0)

        t_median = np.median(times) * 1000
        if t_ref is None:
            t_ref = t_median
        speedup = t_ref / t_median if t_median > 0 else float('inf')

        print(f"{nblocks:>8} {blk_size:>10} {t_median:>12.3f} {speedup:>10.2f}x")

    print()


if __name__ == '__main__':
    print("Block Matrix Optimization Benchmark")
    print("=" * 72)
    print()

    # UHFk _diag() with varying norb
    bench_uhfk_diag(norb_list=[2, 4, 8, 16, 32], nvol=8)

    # RPA _solve_rpa() with varying norb
    bench_rpa_solve(norb_list=[2, 4, 8, 16], nvol=8, nmat=16)

    # UHFk _diag() with varying number of blocks for nd=32
    bench_uhfk_nblocks(nd=32, nvol=8,
                       nblocks_list=[1, 2, 4, 8, 16, 32])
