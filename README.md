# H-wave

H-wave is a program for performing unrestricted Hartree-Fock (UHF) approximation and random phase approximation (RPA) for itinerant electron systems. UHF and RPA correspond to simple approximations that deal with fluctuations up to first order and enable analyses of electron correlation effects in materials at a low computational cost. The input files describing the one-body and two-body interactions are based on the Wannier90 format[1]. This allows smooth connection for the software packages that derive the effective models from first principles calculations, such as RESPACK[2], to the analyses of the effective model with H-wave.

[1] G. Pizzi et al, J. Phys.: Condens. Matter 32 165902 (2020).
[2] K. Nakmura, Y. Yoshimoto, Y. Nomura et al., Comp. Phys. Commun. 261, 107781 (2021).


## Methods

Hartree-Fock and Random Phase approximation

## Target models

Hubbard model, multi-orbital Hubbard model

## Available physical quantities

ground-state energy, free energy, charge and spin susceptibilities, etc.

## Requirement

Python3 with numpy, scipy, and other library packages

## Install

- From PyPI

``` bash
python3 -m pip install hwave
```

- From source (if you modify the program)

``` bash
python3 -m pip install DIRECTORY_OF_THE_REPOSITORY
```

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
