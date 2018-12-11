Title: Writing Software Using the ISOFT Framework
Author: Chris MacMackin
Date: November 2018 

## The Basics

Before proceeding with using ISOFT, you need to initialise the logging
tool which it uses, provided by
[flogging](https://cmacmackin.github.io/flogging/). This allows you to
specify what level of information to output to the terminal and the
log file. This is done as follows:
```fortran
use logger_mod, only: error, info, logger => master_logger, logger_init
! ...
call logger_init('isoft.log', error, info, info)
logger%trivia('isoft', 'Logger initialised successfully.')
logger%info('isoft', 'Welcome to ISOFT.')
```
The first line in this code fragment imports the module providing the
[logger](https://cmacmackin.github.io/flogging/type/logger.html)
object. The parameters `error` and `info` specify different priorities
of messages which can be sent to the logger. The `master_logger` is a
global logging object which is used within the ISOFT library. It is
initialised using the
[logger_init](https://cmacmackin.github.io/flogging/proc/logger_init.html)
routine, which specifies the name of the log file and what priorities
of messages get printed to the terminal and the log file. Methods of the
`logger` object can be used to write out messages of different priority
levels: "debug", "trivia", "info", "warning", "error", and "fatal".

It is also necessary to initialise the HDF5 library, in order to allow I/O
to be performed. This is done by running
```fortran
use hdf5
! ...
integer :: hdf_err
! ...
call h5open_f(hdf_err)
if (hdf_err < 0) error stop
```
Similarly, at the end of your simulation you should deactivate HDF5 by
running
```fortran
call h5close_f(hdf_err)
if (hdf_err < 0) error stop
```

To run simulations with ISOFT, you must first build a
[[cryosphere(type)]] object using the [[cryosphere(type):initialise]]
method. This takes, as arguments, allocatable polymorphic objects of class
[[glacier(type)]] and [[basal_surface(type)]]. These objects must be
initialised with concrete types such as [[ice_shelf(type)]] and
[[plume(type)]] (likely to be the most frequently used implementations).

[[Glacier(type)]] and [[basal_surface(type)]] objects will typically
feature `initialise` methods. However, such methods are unique to each
sub-type. Thus, the concrete type of the [[glacier(type)]] or
[[basal_surface(type)]] objects must be known when they are being
initialised The simplest way to do this is to work with allocatable
variables with those concrete types and then use the intrinsic
`move_alloc` routine to transfer the object a polymorphic variable:
```fortran
type(ice_shelf), allocatable :: shelf_obj
class(glacier), allocatable  :: glacier_obj

allocate(shelf_obj)
call shelf_obj%initialise(...)
call move_alloc(shelf_obj, glacier_obj)
```

The initialisation methods allow you to choose various parameter
values, specify initial conditions, set boundary conditions, and
select which parameterisations to use in the simulation. A number of
additional objects must be created to accomplish the latter two. As
before, the arguments must be allocatable objects of polymorphic
type. However, classes representing boundary conditions and
parameterisations have constructor functions, meaning that it is not
necessary to use the `move_alloc` intrinsic. Instead, these objects
can be initialised as below:
```fortran
type(abstract_viscosity), allocatable :: viscosity
allocate(viscosity, source=newtonian_viscosity(1.0d0))
```
which creates an [[abstract_viscosity(type)]] object using the
[[newtonian_viscosity(proc)]] constructor.

A full description of all the derived types which can be used to
construct a [[cryosphere(type)]] object is beyond the scope of this
section. An overview can be found in the
[discussion of the code design](../3-codeDesign/1-representing.html)
and by consulting the [API documentation](|url|/lists/types.html).

After creating a [[cryosphere(type)]] object, you can optionally
initialise it using the output from a previous simulation. This is
done with the [[cryosphere(type):read_data]] method. The simulation is
run by integrating forward in time using the
[[cryosphere(type):integrate]] method. At any point during the
simulation, the state of the cryosphere can be written to the disk as
an HDF5 file using the [[cryosphere(type):write_data]] method. This
output can be used to initialise future simulations or be analysed
using the provided set of
[plotting scripts](../4-plotting/index.html).

## An Example

Below is an example of a simple ISOFT simulation, designed to run to a
steady state. It is provided as `main.f90` in the top directory of ISOFT.

```fortran
{!main.f90!}
```
