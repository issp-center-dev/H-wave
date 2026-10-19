"""SU(2) cross-check of the two bond-resolved RPA gates on a Heisenberg bond
(GitHub issue #181, Exchange promotion).

An off-site Ising term ``J_I`` and an off-site Exchange term ``J_E`` on the
same bond are spin-rotation invariant exactly when ``J_E = -2 J_I`` in the
input-file convention (each row declared together with its partner): the
Hamiltonian then commutes with the total ``S^-`` (pinned below by exact
diagonalization of the term builders the gates are validated with). In a
paramagnet the longitudinal spin susceptibility and the transverse one
then coincide, so the static channel-0 output of the LONGITUDINAL gate
(``longitudinal_bond_chiq_s_static``, which carries the off-site Exchange
through its bond-diagonal Fock crossing) must equal that of the
TRANSVERSE gate (``chiq_pm_static``, which carries the Exchange spin flip
on the local slots) -- two independent constructions, each pinned against
exact diagonalization on its own.

Two-sided: away from the SU(2) ratio the two differ, and with the
Exchange removed from the longitudinal side only they differ as well.
"""
import os
import shutil
import tempfile
import unittest

import numpy as np

import hwave.qlmsio.read_input_k as read_input_k
import hwave.solver.rpa as rpa_mod

_SRC = "tests/equivalence_input/orb2"
#: input coupling of the Ising rows (offsite_ising.dat: inter-orbital (1,2) bond, R = +x)
J_ISING = 0.2


def _solve(j_exchange, gate, j_exchange_longitudinal=None):
    d = tempfile.mkdtemp()
    out = tempfile.mkdtemp()
    try:
        for f in ("geom.dat", "transfer.dat", "offsite_ising.dat"):
            shutil.copy(os.path.join(_SRC, f), d)
        with open(os.path.join(_SRC, "offsite_exchange.dat")) as fr, \
                open(os.path.join(d, "offsite_exchange.dat"), "w") as fw:
            fw.write(fr.read().replace("0.200000000000", "{:.12f}".format(j_exchange)))
        idict = {"path_to_input": d, "Geometry": "geom.dat", "Transfer": "transfer.dat",
                 "Ising": "offsite_ising.dat", "Exchange": "offsite_exchange.dat"}
        par = {"T": 0.5, "filling": 0.5, "CellShape": [4, 4, 1], "SubShape": [1, 1, 1],
               "Nmat": 64, gate: True}
        calc_type = "ring"
        if gate == "transverse_bond_channels":
            # a declared zero field: spin-diag Green function with G_up == G_down
            with open(os.path.join(d, "extern.dat"), "w") as fw:
                fw.write("zero field\n2\n2\n 1 1\n 0 0 0 1 1 1.0 0.0\n 0 0 0 2 2 1.0 0.0\n")
            idict["Extern"] = "extern.dat"
            par["coeff_extern"] = 0.0
            calc_type = "ring+ladder"
        r = read_input_k.QLMSkInput({"path_to_input": d, "interaction": idict})
        info = {"mode": "RPA", "param": par, "enable_spin_orbital": False,
                "calc_scheme": "general", "calc_type": calc_type}
        s = rpa_mod.RPA(r.get_param("ham"), {}, info)
        if j_exchange_longitudinal is not None:
            tbl = s.ham_info.param_ham["Exchange"]
            for key in list(tbl):
                tbl[key] = j_exchange_longitudinal
        gi = r.get_param("green")
        s.solve(gi, out)
        key = ("chiq_pm_static" if gate == "transverse_bond_channels"
               else "longitudinal_bond_chiq_s_static")
        return np.asarray(gi[key])
    finally:
        shutil.rmtree(d, ignore_errors=True)
        shutil.rmtree(out, ignore_errors=True)


class TestSU2Ratio(unittest.TestCase):

    def test_heisenberg_ratio_commutes_with_total_spin(self):
        from tests import ed_oracle_util as eu
        from tests.test_bond_transverse_ed import _terms_exchange_offsite, _terms_ising_offsite
        fx = eu.EDFixture(3, 1, {(0, 0): -1.0}, [0.0], 0.5, 0.0)
        C = fx.annihilators()
        CD = [c.conj().T for c in C]
        Sm = sum(CD[fx.mode(j, 0, 1)] @ C[fx.mode(j, 0, 0)] for j in range(fx.L))
        for ratio, commutes in ((-2.0, True), (-1.0, False), (-4.0, False)):
            terms = (_terms_ising_offsite(fx, 0, 0, 1, J_ISING)
                     + _terms_exchange_offsite(fx, 0, 0, 1, ratio * J_ISING))
            H = eu.h_int_from_terms(fx, terms)
            dev = np.abs(H @ Sm - Sm @ H).max()
            if commutes:
                self.assertLess(dev, 1e-14, ratio)
            else:
                self.assertGreater(dev, 0.1, ratio)


class TestLongitudinalEqualsTransverse(unittest.TestCase):

    def test_heisenberg_bond(self):
        zz = _solve(-2.0 * J_ISING, "longitudinal_bond_channels")
        pm = _solve(-2.0 * J_ISING, "transverse_bond_channels")
        self.assertGreater(np.abs(zz).max(), 0.1)
        np.testing.assert_allclose(zz, pm, rtol=0, atol=1e-12)

    def test_away_from_the_su2_ratio_they_differ(self):
        zz = _solve(-1.0 * J_ISING, "longitudinal_bond_channels")
        pm = _solve(-1.0 * J_ISING, "transverse_bond_channels")
        self.assertGreater(np.abs(zz - pm).max(), 1e-3)

    def test_without_the_longitudinal_exchange_they_differ(self):
        """The longitudinal side with the Exchange coefficient zeroed (its
        bond blocks gone) misses the transverse side on the SU(2) bond."""
        zz = _solve(-2.0 * J_ISING, "longitudinal_bond_channels", j_exchange_longitudinal=0.0)
        pm = _solve(-2.0 * J_ISING, "transverse_bond_channels")
        # measured 2.5e-4 (the bond crossing is small at this temperature),
        # eight orders above the 1e-12 agreement of the full longitudinal side
        self.assertGreater(np.abs(zz - pm).max(), 1e-5)


if __name__ == "__main__":
    unittest.main()
