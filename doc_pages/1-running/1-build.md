Title: Compiling the Code
Author: Chris MacMackin
Date: November 2018 

## Dependencies

ISOFT depends on a number of external pieces of software. Mostly these
are libraries which it uses, although there are also a few programs
which are used at build-time. Except for the compiler, the latter are
all Python-based and the Makefile will automatically install them in a
virtual environment. Of the libraries, those which are widely used
must be installed on your system prior to building ISOFT. The
remainder are distributed with ISOFT and will automatically be
compiled as part of the build process.  Builds have been tested with
the `gfortran` compiler and are known to work with v6.2. In principle,
subsequent releases of `gfortran` should also work, barring any
regressions.

The following libraries and programs must be installed on your operating system
prior to starting the build process. The corresponding package name in
Ubuntu is given in parentheses.

- [FFTW3](http://www.fftw.org/) (libfftw3-dev)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/) (libhdf5-dev)
- An implementation of [BLAS](http://www.netlib.org/blas/) (e.g., libopenblas-dev)
- [LAPACK](http://www.netlib.org/lapack/) (liblapack-dev)
- [virtualenv](https://virtualenv.pypa.io/en/latest/) (python-virtualenv)
- [gfortran](https://gcc.gnu.org/wiki/GFortran) (gfortran-6)

The following programs are required during the build process. They
will be installed in a virtual environment called `buildtools` in the
top of the ISOFT directory.

- [fypp](https://fypp.readthedocs.io/en/stable/fypp.html), a preprocessor for Fortran
- [FoBiS.py](https://github.com/szaghi/FoBiS/wiki), a simple build system for Fortran
- [FORD](https://github.com/Fortran-FOSS-Programmers/ford/wiki), a documentation tool for Fortran

During the build process, the following libraries will automatically be built.

- [FACE](https://szaghi.github.io/FACE/index.html)
- [factual](https://github.com/cmacmackin/factual)
- [flogging](https://cmacmackin.github.io/flogging/)
- [lapack95](http://www.netlib.org/lapack95/)
- [nitsol](https://www.osti.gov/biblio/433349)
- [PENF](http://szaghi.github.io/PENF/index.html)
- [pFUnit](http://pfunit.sourceforge.net/)

## The Makefile

A Makefile is provided which can handle the build process on Unix-like
operating systems in most cases. It features the following build
targets:

- __all__: alias for __lib__ (default)
- __all_objects__: rebuilds all object files for the library
- __clean__: deletes object files, module files, dependency files, and Emacs backup files
- __doc__: generates HTML documentation for ISOFT
- __script__: writes a script for compiling and linking programs which use the ISOFT library
- __tests__: builds and runs the unit tests
- __lib__: builds the static library of ISOFT routines and derived types

At the top of the Makefile are two definitions used for finding
libraries installed on your system. These are
```Makefile
SYSTEM_LIB := /usr/lib       # Paths containing external libraries
SYSTEM_INC := /usr/include   # Paths containing module/include files
                             # for external libraries
```
These defaults should work on Ubuntu-like operating systems. 
They may be overridden if necessary. Multiple paths may also be
provided, separated by white-space.

## Building

After having amended the Makefile as described above, building ISOFT
should be straightforward. Simply run
```
make lib
```
to compile the necessary dependencies and create the `libisoft.a` file.
Then run
```
make script
```
to generate a script which can [compile and link programs](./3-compile.html)
using ISOFT.

## Unit Tests

An extensive set of unit tests has been written, using the
[pFUnit](http://pfunit.sourceforge.net/) framework. These ensure that
methods behave as expected and, where possible, test that simulations
converge to analytically-predicted solutions. These tests can be built
and run with the command
```
make tests
```
