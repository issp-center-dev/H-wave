# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

H-wave is a Python package for Unrestricted Hartree-Fock (UHF) and Random Phase Approximation (RPA) calculations on itinerant electron systems (Hubbard and multi-orbital Hubbard models). Developed by ISSP/University of Tokyo. Licensed under GPL v3.

## Build and Install

```bash
# Install from source (uses Poetry build backend)
python3 -m pip install .

# Install pytest for testing
python3 -m pip install pytest
```

Dependencies: Python 3.7+, numpy, scipy, requests, tomli.

## Running Tests

```bash
# Run all tests
pytest -v

# Run a single test class
pytest tests/test_uhf.py::TestUHFk -v

# Run a single test
pytest tests/test_uhf.py::TestUHFk::test_CoulombIntra -v

# RPA tests
pytest tests/test_rpa_1orb.py -v
pytest tests/test_rpa_2orb.py -v
```

Tests use `unittest` with `pytest` as runner. Each test runs a solver against reference data in `tests/<solver>/<case>/output_ref/energy.dat` and compares with tolerance `atol=1.0e-8`. Tests must be run from the repository root (they use relative paths like `tests/uhfr/coulombintra`).

## CLI Entry Points

```bash
hwave <input.toml>       # Main solver
hwave_dos [args]         # Density of States
```

Defined in `pyproject.toml` under `[tool.poetry.scripts]`.

## Architecture

### Solver Pipeline (`src/hwave/qlms.py`)

The `run()` function orchestrates the full pipeline:
1. Parse TOML config → determine mode (`UHFr`, `UHFk`, or `RPA`)
2. Read input via `qlmsio` (real-space or k-space reader)
3. Instantiate solver → `solver.solve()` → `solver.save_results()`

### Solver Hierarchy (`src/hwave/solver/`)

- **`base.py`** — `solver_base`: common parameter validation, defaults, range checks
- **`uhfr.py`** — `UHFr(solver_base)`: real-space UHF for small cluster systems
- **`uhfk.py`** — `UHFk(solver_base)`: k-space UHF with FFT, periodic boundaries, band structure. Initialization has ordered phases: `_init_mode()` → `_init_param()` → `_init_lattice()` → `_init_orbit()` → `_init_wavevec()` → `_init_interaction()` → `_check_interaction()`
- **`rpa.py`** — `RPA(solver_base)`: Random Phase Approximation for susceptibilities, builds on k-space framework

### I/O Layer (`src/hwave/qlmsio/`)

- **`read_input.py`** — `QLMSInput`: reads real-space format (Trans, CoulombIntra/Inter, Hund, Exchange, Ising, PairLift, PairHop, Interall)
- **`read_input_k.py`** — `QLMSkInput`: reads k-space/Wannier90 format with FFT back-transformation
- **`wan90.py`** — Wannier90 format utilities

### Key Patterns

- **CaseInsensitiveDict** (from `requests.structures`) used throughout for flexible parameter keys
- **Logging hierarchy**: `qlms` → `qlms.solver` → `qlms.solver.uhfk` etc.
- **Config-driven**: all calculations configured via TOML input files
- **`enable_spin_orbital` mode** in UHFk: collapses spin degree of freedom (`ns=1`), currently only supports Transfer interaction term

### Input/Output Formats

- **Input**: TOML config files, Wannier90 geometry/transfer files, plain text interaction definitions
- **Output**: plain text (`energy.dat`), NumPy `.npz` files (Green's functions, susceptibilities)
