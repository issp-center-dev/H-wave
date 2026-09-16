"""Double-counting energies must equal -Tr(H_int rho)/2.

The unequal intra-cell/inter-cell hoppings make the two directed orbital
coherences differ. The on-site Exchange case fixes the sign and total-cell
normalization independently of the off-site regression.

Before the fix, plain UHFk gave these (energy, -Tr(H_int rho)/2) values:
Exchange off-site: (-0.000996118502, -0.540038047857), error +0.539041929355.
Ising off-site:    (+0.000927187264, +0.154288384881), error -0.153361197617.
PairHop off-site:  (-0.000001837374, -0.540038047857), error +0.540036210483
                  with Fock either enabled or disabled.
On-site Exchange and PairHop already agreed to within 1e-15.

PairHop's Hamiltonian contracts conj(G[r,v,b,u,a]) with the expectation
conj(G[r,t,b,s,a]). Relabeling its spin table gives the direct product
conj(G[r,v,b,s,a])*conj(G[r,u,b,t,a]) and the cross product
conj(G[r,u,b,s,a])*conj(G[r,v,b,t,a]). Both previously used the opposite
bond. Spin-orbital conversion only permutes basis indices, so the same
contractions apply; complex spin mixing exercises the cross product.
"""
import numpy as np
import pytest

from hwave.solver.uhfk import UHFk


@pytest.mark.parametrize("spin_orbital", [False, True], ids=["plain", "spin-orbital"])
@pytest.mark.parametrize("interaction, offset, fock", [
    ("Exchange", 0, True), ("Exchange", 1, True), ("Ising", 1, True),
    ("PairHop", 0, True), ("PairHop", 1, True), ("PairHop", 1, False),
], ids=["onsite-exchange", "offsite-exchange", "offsite-ising", "onsite-pairhop",
        "offsite-pairhop", "offsite-pairhop-no-fock"])
def test_fock_double_counting_matches_mean_field(interaction, offset, fock, spin_orbital):
    transfer = {
        ((0, 0, 0), (0, 0)): -0.4,
        ((0, 0, 0), (1, 1)): 0.6,
        ((0, 0, 0), (0, 1)): -0.3,
        ((0, 0, 0), (1, 0)): -0.3,
        ((1, 0, 0), (0, 1)): -1.0,
        ((-1, 0, 0), (1, 0)): -1.0,
    }
    norb = 2
    if spin_orbital:
        norb = 4
        transfer = {
            (r, (2*a+s, 2*b+s)): value
            for (r, (a, b)), value in transfer.items() for s in range(2)
        }
        # Complex spin mixing also exercises PairHop's cross contraction.
        transfer[((0, 0, 0), (0, 1))] = 0.12j
        transfer[((0, 0, 0), (1, 0))] = -0.12j
    ham = {
        "Geometry": {"norb": norb, "rvec": np.eye(3),
                     "center": np.zeros((norb, 3)), "degree": np.ones(norb)},
        "Transfer": transfer,
        interaction: {((offset, 0, 0), (0, 1)): 0.4,
                      ((-offset, 0, 0), (1, 0)): 0.4},
    }
    params = {"CellShape": [4, 1, 1], "SubShape": [1, 1, 1],
              "Ncond": 8, "2Sz": None, "IterationMax": 1000,
              "eps": 1e-10, "Mix": 0.5, "RndSeed": 7, "T": 0.05}
    solver = UHFk(ham, {"print_level": 0},
                  {"mode": "UHFk", "enable_spin_orbital": spin_orbital,
                   "flag_fock": fock}, params)
    solver.solve({"initial_mode": "random"}, ".")
    assert solver.physics["Rest"] < 1e-10

    # Rebuild both quantities from the same final (mixed) Green function.
    solver._make_ham()
    solver._calc_energy()
    green_k = np.fft.ifftn(
        solver.Green.reshape(*solver.shape, solver.nd, solver.nd),
        axes=(0, 1, 2), norm="forward",
    ).reshape(solver.nvol, solver.nd, solver.nd)
    # G_ij = <c_i^dag c_j>, whereas rho_ij = <c_j^dag c_i>.
    rho_k = green_k.transpose(0, 2, 1)
    h_int = solver.ham - solver.ham_trans
    expectation = np.trace(h_int @ rho_k, axis1=1, axis2=2).sum()
    assert abs(expectation.imag) < 1e-12
    expected = -expectation.real / 2
    actual = solver.physics["Ene"][interaction]
    if offset:
        green = (solver._so_to_virtual_green(solver.Green)
                 if spin_orbital else solver.Green)
        assert abs(green[1, 0, 0, 0, 1] - green[-1, 0, 0, 0, 1]) > 0.1
    print(f"{interaction} offset={offset} SO={spin_orbital} Fock={fock}: "
          f"energy={actual.real:.12f}, oracle={expected:.12f}, "
          f"error={actual.real - expected:.12g}")
    assert actual == pytest.approx(expected, rel=0, abs=1e-10)
