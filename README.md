# H-wave

<div align="center">
  <img src="docs/figs/Hwave_logo.png" alt="H-wave Logo" width="200">
</div>

[![Run tests](https://github.com/issp-center-dev/H-wave/actions/workflows/run_tests.yml/badge.svg)](https://github.com/issp-center-dev/H-wave/actions/workflows/run_tests.yml)
[![CI Python 3.9+](https://github.com/issp-center-dev/H-wave/actions/workflows/ci-python39.yml/badge.svg)](https://github.com/issp-center-dev/H-wave/actions/workflows/ci-python39.yml)
[![PyPI version](https://img.shields.io/pypi/v/hwave)](https://pypi.org/project/hwave/)
[![Python](https://img.shields.io/badge/python-3.9%2B-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Doc](https://img.shields.io/badge/doc-Manual-blue)](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en/doc/manual)

H-wave is a Python package for performing unrestricted Hartree-Fock (UHF) approximation and random phase approximation (RPA) for itinerant electron systems. UHF and RPA correspond to simple approximations that deal with fluctuations up to first order and enable analyses of electron correlation effects in materials at a low computational cost. The input files describing the one-body and two-body interactions are based on the Wannier90 format[1]. This allows smooth connection for the software packages that derive the effective models from first principles calculations, such as RESPACK[2], to the analyses of the effective model with H-wave.

[1] G. Pizzi et al, J. Phys.: Condens. Matter 32 165902 (2020).  
[2] K. Nakmura, Y. Yoshimoto, Y. Nomura et al., Comp. Phys. Commun. 261, 107781 (2021).

## Features

- **Unrestricted Hartree-Fock (UHF) approximation** for itinerant electron systems
- **Random Phase Approximation (RPA)** for correlation effects analysis
- **Support for various interaction types**: Coulomb intra/inter, Exchange, Hund, Ising, PairHop, PairLift
- **Wannier90 format compatibility** for seamless integration with first-principles calculations
- **Comprehensive test suite** with automated CI/CD pipeline
- **Multi-orbital Hubbard model support**

## Methods

- **Hartree-Fock approximation**: Ground state energy calculation
- **Random Phase approximation**: Susceptibility and correlation analysis

## Target Models

- Hubbard model
- Multi-orbital Hubbard model
- Extended Hubbard models with various interaction types

## Available Physical Quantities

- Ground-state energy
- Free energy
- Charge and spin susceptibilities
- Green's functions
- Eigenvalues and eigenvectors

## Requirements

- Python 3.9 or higher
- NumPy (^1.14)
- SciPy (^1.7)
- Requests (^2.28.1)
- Tomli (^2.0.1)

## Installation

### From PyPI (Recommended)

```bash
pip install hwave
```

### From Source

For development or if you need to modify the program:

```bash
git clone https://github.com/issp-center-dev/H-wave.git
cd H-wave
pip install -e .
```

### Using Poetry (Development)

```bash
git clone https://github.com/issp-center-dev/H-wave.git
cd H-wave
poetry install
```

## Quick Start

### Basic Usage

```python
import hwave.qlms

# Run UHF calculation
hwave.qlms.run(input_dict=params)
```

### Command Line Interface

```bash
# Run UHF calculation
hwave input.toml

# Calculate DOS
hwave_dos input.toml
```

## Testing

The project includes a comprehensive test suite that covers:

- **UHF calculations** (UHFr and UHFk solvers)
- **RPA calculations** with various interaction types
- **Multi-orbital systems**
- **Spin-orbital coupling**

### Running Tests

```bash
# Run all tests
python -m unittest discover tests/ -v

# Run specific test modules
python -m unittest tests.test_uhf -v
python -m unittest tests.test_rpa_1orb -v

# Using pytest (if installed)
pytest tests/ -v
```

### Continuous Integration

The project uses GitHub Actions for automated testing:

- **Multi-version testing**: Python 3.9, 3.10, 3.11, 3.12
- **Code quality checks**: flake8, black, isort, mypy
- **Automated testing**: Runs on every push and pull request to main/develop branches

## Documentation

- [User Manual](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en/doc/manual)
- [API Documentation](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en/doc/manual)
- [Tutorial Examples](https://github.com/issp-center-dev/H-wave/tree/main/docs/tutorial)

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Make your changes
4. Run tests: `python -m unittest discover tests/ -v`
5. Commit your changes: `git commit -m "Add your feature"`
6. Push to the branch: `git push origin feature/your-feature-name`
7. Create a Pull Request

### Code Quality

The project enforces code quality through:

- **Black**: Code formatting
- **isort**: Import sorting
- **flake8**: Linting
- **mypy**: Type checking

Run these tools before submitting:

```bash
black src/ tests/
isort src/ tests/
flake8 src/ tests/
mypy src/
```

## License

The distribution of the program package and the source codes for H-wave follow
GNU General Public License version 3
([GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html)).

Copyright (c) <2022-> The University of Tokyo. All rights reserved.

This software was developed with the support of
"Project for Advancement of Software Usability in Materials Science"
of The Institute for Solid State Physics, The University of Tokyo.

We would appreciate it if you cite the following article in your research with H-wave:
H-wave -- A Python package for the Hartree-Fock approximation and the random phase approximation,
Tatsumi Aoyama, Kazuyoshi Yoshimi, Kota Ido, Yuichi Motoyama, Taiki Kawamura, Takahiro Misawa, Takeo Kato, and Akito Kobayashi,
[Computer Physics Communications, 298, 109087 (2024)](https://doi.org/10.1016/j.cpc.2024.109087).

## Official page

- [H-wave project site](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en)
- [Software repository](https://github.com/issp-center-dev/H-wave)
- [User Manual](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en/doc/manual)
- [Data Repository](https://datarepo.mdcl.issp.u-tokyo.ac.jp/repo/23)

## Authors

Kazuyoshi Yoshimi (ISSP, Univ. of Tokyo),
Yuichi Motoyama (ISSP, Univ. of Tokyo),
Tatsumi Aoyama (ISSP, Univ. of Tokyo),
Kota Ido (ISSP, Univ. of Tokyo),
Takahiro Misawa (ISSP, Univ. of Tokyo),
Taiki Kawamura (Nagoya Univ.),
Akito Kobayashi (Nagoya Univ.),
Takeo Kato (ISSP, Univ. of Tokyo)
