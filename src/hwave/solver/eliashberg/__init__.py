"""The linearized Eliashberg solver of hwave_sc, one module per layer:
config, band, vertex, kernel, eigen, io, output, bond, so_solve, static
(the entry), and the dynamic solvers dynamic, dynamic_so, bond_dynamic,
bond_io. hwave.sc is the entry shell; it re-exports only the public API
(calc_eliashberg, main, reject_spin_orbital_mode and the public
constants), and nothing under hwave.solver imports hwave.sc.
"""
