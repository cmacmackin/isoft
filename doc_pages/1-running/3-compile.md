Title: Linking to the Library
Author: Chris MacMackin
Date: November 2018 

As [previously discussed](./1-build.html), by running
```
make script
```
you will produce a script which can be used to compile programs that
use ISOFT. It will contain information on the location of all the
library files needed for ISOFT to run.

This script is called `compile.sh` and has the following call signature:
```
./compile.sh [ main-file-name [ executable-name ] ]
```
where `main-file-name` is the program file to be compiled and
`executable-name` is the name to give the resulting executable. These
default to `main.f90` and `isoft`, respectively.

The [example program](./2-writing.html) provided earlier can be
compiled and run as follows:
```fortran
./compile.sh
./isoft
```
The provided script assumes that your ISOFT-dependent code is all in a
single source file. It is trivial to modify the script should further
files need to be compiled and linked.

