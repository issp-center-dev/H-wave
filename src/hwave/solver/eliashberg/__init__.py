"""The linearized Eliashberg solver of hwave_sc, one module per layer
(dev/sc-split-plan.md): config, band, vertex, kernel, eigen, io, output,
bond, so_solve, static (the entry), and the dynamic solvers dynamic,
dynamic_so, bond_dynamic, bond_io. hwave.sc is the shell that re-exports
what the tests and external scripts still read from it; nothing under
hwave.solver imports hwave.sc.
"""
