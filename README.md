H-wave
========

H-wave is a program for performing unrestricted Hartree-Fock (UHF) approximation
for itinerant electron systems.
The UHF method approximates two-body interaction terms into one-body ones by
taking account of the fluctuations up to the first order. The wave functions
and energies are determined self-consistently by an iterative method.

H-wave treats two types of UHF approximation: one is real-space and the other is
wavenumber-space UHF methods using translational symmetry. In the wavenumber space
UHF method, the input file defining the one-body and two-body interactions is
based on the Wannier90 format, and the program can be smoothly connected to
the softwares for deriving effective models from first-principles calculations. 

### Methods
Hartree-Fock approximation

### Target models
Hubbard model, multi-orbital Hubbard model

### Available physical quantities
ground-state energy, free energy, etc.

## Requirement
Python3 with numpy, scipy, and other library packages

## Install
You can install H-wave and also get a manual for H-wave from a
[release page](https://github.com/issp-center-dev/hwave/releases).

## License
The distribution of the program package and the source codes for H-wave follow
GNU General Public License version 3
([GPL v3](https://www.gnu.org/licenses/gpl-3.0.en.html)).

Copyright (c) <2022-> The University of Tokyo. All rights reserved.

This software was developed with the support of
"Project for Advancement of Software Usability in Materials Science"
of The Institute for Solid State Physics, The University of Tokyo.

## Official page

- [H-wave project site](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en)
- [Software repository](https://github.com/issp-center-dev/H-wave)
- [User Manual](https://www.pasums.issp.u-tokyo.ac.jp/h-wave/en/doc/manual)

## Authors
Kazuyoshi Yoshimi,
Yuichi Motoyama,
Tatsumi Aoyama,
Kota Ido,
Takahiro Misawa,
Taiki Kawamura
Akito Kobayashi,
Takeo Kato
