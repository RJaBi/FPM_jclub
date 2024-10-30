# jclub

A Repository to demonstrate usage of the fortran-package manager and the meson-build system for a journal club in the school of mathematics, trinity college dublin.

As I am a lattice gauge theorist, that is the examples I am using.

In particular, I have implemented the plaquette measurement and also the anisotropy and lattice spacing measurement using the Wilson flow.

Generated with cookiecutter template:
[https://github.com/SalvadorBrandolin/fortran_meson_py](https://github.com/SalvadorBrandolin/fortran_meson_py)

---
What does the code actually do?

The `main` program (run with `fpm run jclub`) reads an openQCD format gaugefield https://luscher.web.cern.ch/luscher/openQCD/ into memory and calculates the average value of the simplest closed loop called the plaquette. The plaquette is one of the simplest observables in lattice field theory. In practice this means a lot of 3x3 complex matrix multiplications.

The `w0RE` program (run with `fpm run --profile release w0RE -- PATHTOINPUTTOML.toml`) performs the R_E method to get the anisotropy (ratio of spatial to temporal lattice spacing) and the lattice spacing (by the w_0 scale) using the gradient flow. This is outlined in https://doi.org/10.48550/arXiv.1205.0781 and https://doi.org/10.1007/JHEP09%282012%29010 . In particular the gradient flow is performed using the code available in the ancillary files of the latter on arXiv at https://arxiv.org/abs/1203.4469 .

There is an example input file `G2_wflow.toml`

Also note that the bspline (https://github.com/jacobwilliams/bspline-fortran) code used seems to have an issue with the check-bounds of the default (or debug) profile for `gfortran 14.1.0-1`. It works fine in the release profile or with the `ifx` compiler


---
Install

You need a fortran compiler, the code has been tested with `gfortran 14.1.0-1` and `ifx (IFX) 2024.2.1 20240711`.

You need the fortran-package-manager https://fpm.fortran-lang.org/. This may be installed in a number of ways. The simplest will be i.e. `pip` or `conda`
