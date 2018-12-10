# ISOFT: Ice Shelf/Ocean Fluid and Thermodynamics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1422483.svg)](https://doi.org/10.5281/zenodo.1422483)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

ISOFT is a piece of software/suite of tools which I developed while
working on my PhD thesis to simulate the evolution of ice shelves
coupled to ocean plumes. It attempts to provide an object oriented,
extensible framework with which to model glacial flow using modern
Fortran. Though developed for ice shelves, it could in principle be
modified to simulated grounded ice dynamics as well. I've published
the code and documentation to GitHub in the hopes that it might be
useful to others.

As much as practical, ISOFT was kept agnostic as to whether it was
handling a 1-D or a 2-D representation of an ice shelf/plume. However,
the current implementation does explicitly assume a 1-D system in a
number of instances and would thus need to be modified to handle 2-D
problems.  1-D simulations were sufficiently fast that they could be
run in serial. Multithreading could easily be implemented in many
parts of the code where arithmetic is performed on arrays. Indeed,
most of these cases are simple enough that a compiler may be able to
parallelise them automatically. More sophisticated approaches
involving message passing would likely be necessary to make 2-D
simulations practical, but this would be far more difficult to
implement and would likely require substantial refactoring of the
code. In particular, the nonlinear solvers would likely need to be
replaced.

## Quick Start

The build-script for FORD assumes you have the libraries and header
files for [FFTW3](http://www.fftw.org/),
[HDF5](https://www.hdfgroup.org/solutions/hdf5/), an implementation of
[BLAS](http://www.netlib.org/blas/), and
[LAPACK](http://www.netlib.org/lapack/) installed on your computer. It
also requires you to have the
[virtualenv](https://virtualenv.pypa.io/en/latest/) installed. ISOFT
is known to compile with `gfortran` v6.2 and, barring regressions,
should work with all subsequent releases.  There are a number of
other, more obscure libraries upon which ISOFT depends, but the
Makefile handles the downloading and building of these itself. You can
build ISOFT on an Ubuntu-like operating system using the commands
below:
```
sudo apt-get install libfftw3-dev libhdf5-dev libopenblas-dev liblapack-dev python-virtualenv gfortran-6
git clone https://github.com/cmacmackin/isoft
cd isoft
make lib            # Builds libisoft.a
make tests          # Builds and runs unit tests
make script         # Creates a compile-script you can use to build and
                    # link your own programs with the ISOFT library
```

This builds a static library which provides derived types for
simulating the evolution of an ice shelf. You initialise these objects
in your own programs and then call their methods to start the
simulation. An example of such a program can be found in the sample
file `main.f90`. This can be built and run as follows:
```
compile.sh main.f90 isoft_example
./isoft_example
```

The `compile.sh` script was generated specifically for your build of
ISOFT and can be used, without modification from any directory on your
computer. It's call signature is
```
./compile.sh [ main-file-name [ executable-name ] ]
```
where `main-file-name` is the program file to be compiled and
`executable-name` is the name to give the resulting executable. These
default to `main.f90` and `isoft`, respectively.

## Documentation

Information on how to install and use ISOFT is available on
[GitHub Pages](https://cmacmackin.github.io/isoft). This includes
sections (taken from my thesis) on the numerical methods and code
design choices which were made. A detailed description of the API
is also provided. This documentation can be generated locally using the
[FORD tool](https://github.com/Fortran-FOSS-Programmers/ford).  To
install FORD in a virtual environment and then run it to generate the
documentation, execute
```
make doc
```
If FORD is already present on your system then you can simply run
```
ford doc.md
```

## License

ISOFT is licensed under the GNU General Public License (GPL) v3.0 or
later. The terms are provided in the file `LICENSE`. The Lesser General
Public License (LGPL) would have been used, but the Chebyshev pseudo-spectral
implementation uses the FFTW3 library, which is licensed under the GPL.

