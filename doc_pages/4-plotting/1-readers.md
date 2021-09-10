Title: Reader Objects
Author: Chris MacMackin
Date: November 2018 

The most basic part of the plotting library is the classes it provides
for reading the HDF5 files which ISOFT produces. These can be found in
the plotting/[[readers.py]] file. The class which you will mostly use
is `ShelfPlumeCryosphere`, which loads files from simulations run with
the [[ice_shelf(type)]] and [[plume(type)]] (or [[static_plume(type)]]
or [[asym_plume(type)]]) derived types.
Readers of this type are created using the constructor
```python
ShelfPlumeCryosphere(hdf_file)
```
where `hdf_file` is the name of an HDF5 file created by ISOFT.

This reader object a number of properties which provide information on
the state of the simulation. Arrays containing simulation variables
may, in principle, be multidimensional. However, as code currently
exists only for 1-D simulations, they are all 1-D in practice (except
`uvec` and `Uvec`, which are 2-D).

- __isoft_version__: Version number of ISOFT which generated the HDF5 file
- __compilation_time__: Time and date at which the ISOFT library used
  to generate the HDF5 file was compiled
- __output_time__: Wall-clock time the HDF5 file was created
- __time__: Time within the simulation when the HDF5 file was created
- __grid__: The coordinates of the grid discretising the domain
- __glacier_type__: The type of [[glacier(type)]] object which created this data
- __h__: Array containing thickness of the glacier's ice
- __uvec__: Array containing glacier ice velocity in vector form, with
  `u = uvec[...,0]` and `v = uvec[...,1]` (if transverse present)
- __u__: Array containing longitudinal glacier velocity
- __v__: Array containing transverse glacier velocity; will raise an error if not present
- __s__: The surface elevation of the glacier, `h*(1.-1./r)`
- __b__: The basal depth of the glacier, `-h/r`
- __kappa__: Array containingTaylor coefficients in _z_ for a Lagrangian tracer field,
  with `kappa[i,...]` corresponding to the coefficient for the
  i-1<sup>th</sup> power term; raises an error if not present in
  simulation output
- __lambd__: The dimensionless melt parameter, \(\lambda\)
- __chi__: The dimensionless stretching parameter, \(\chi\)
- __basal_surface_type__: The type of [[basal_surface(type)]] object to create this data
- __D__: Array containing the thickness of the subglacial plume
- __Uvec__: Array containing plume ice velocity in vector form, with `U = Uvec[...,0]` and `V = Uvec[...,1]` (if transverse present)
- __U__: Array containing longitudinal plume velocity
- __V__: Array containing transverse plume velocity; will raise an error if not present
- __S__: Array containing plume salinity
- __T__: Array containing plume temperature
- __delta__: The buoyancy correction parameter, \(\delta\)
- __mu__: The dimensionless drag coefficient, \(\mu\)
- __nu__: The dimensionless eddy diffusivity, \(\nu\)
- __r__: The density ratio for ocean water and glacier ice, \(r\)
- __phi__: The dimensionless Coriolis parameter, \(\Phi\)

A simple example of using a reader object to make a plot is provided below:
```python
from plotting.readers import ShelfPlumeCryosphere
import matplotlib.pyplot as plt

cryo = ShelfPlumeCryosphere('isoft-0000.h5')
plt.plot(plt.grid, plt.h)
plt.xlabel('$x$')
plt.ylabel('$h$')
plt.show()
```
