Title: General Plotting Tools
Author: Chris MacMackin
Date: November 2018 

For convenience, a number of additional tools are provided for
analysing ISOFT output. These are organised into a series of
additional Python modules.

## plotting/[[calculus.py]] {: #differentiate }

Provides a `Differentiator` class which can be used to calculate
derivatives of ISOFT output using the
[Chebyshev pseudospectral](../2-numerics/2-spatial.html)
method. Rather than the previously described algorithm using the FFT,
a differentiation matrix is used here
[(Trefethen, 2000)](../6-bibliog.html#Trefethen2000). The constructor
for this object is
```python
Differentiator(size, lower=0.0, upper=1.0)
```
where `size` is the number of Chebyshev pseudo-spectral nodes in the
data to be differentiated, `lower` is the lower bound of the domain,
and `upper` is the upper bound. The differentiator object is applied
by calling it like a function, with the array to be differentiated
passed as an argument.

```python
from plotting.readers import ShelfPlumeCryosphere
from plotting.calculus import Differentiator

cryo = ShelfPlumeCryosphere('isoft-0000.h5')
diff = Differentiator(cryo.grid.size, cryo.grid[-1], cryo.grid[0])
dh_dx = diff(cryo.h)
```


## plotting/[[dimensionless_nums.py]]

This module provides routines to calculate dimensionless numbers used
in fluid dynamics. To date, only the
[Froude number](https://en.wikipedia.org/wiki/Froude_number) has been
implemented. This has routine
```python
froude(coef, U, drho, D)
```
where `coef` is a dimensionless coefficient \(1/\sqrt{\delta}\), `U =
ShelfPlumeCryosphere.Uvec` is the vector velocity, `drho` is a
density difference causing buoyant forcing, and `D` is thickness of
the fluid layer.


## plotting/[[entrainment.py]]

This provides classes which can calculate the entrainment rate for a
subglacial plume. Two such classes are provided: one for the
parameterisation of [Jenkins (1991)](../6-bibliog.html#Jenkins1991)
and another for that of
[Kochergin (1987)](../6-bibliog.html#Kochergin1987). The constructor
for the first has the form
```python
Jenkins1991Entrainment(coefficient, size, lower=0.0, upper=1.0)
```
where `coefficient` is the dimensionless entrainment coefficient
\(E_0/\delta\) (typically with value 1) and the remaining arguments
have the same meaning as those in the constructor for the
[Differentiator type](#differentiate). The constructor for the latter
has the form
```python
Kochergin1987Entrainment(coefficient, delta)
```
where `coefficient` is \(K = c_L^2 x_0/D_0\), as described in
[equation 16](../2-numerics/1-equations.html#mjx-eqn-eqent-koch-coef-nd),
and `delta` is the buoyancy correction \(\delta\).

For both of these classes, the entrainment is calculated by calling
the object like a function. For an entrainment object named `ent`, the
call signature is
```python
ent(U, D, b, rho_diff)
```
where `U` is the vector velocity of the plume, `D` is the plume
thickness, `b` is the basal depth of the ice shelf, and `rho_diff` is
the density difference between the ambient ocean and the plume. The
last argument is not used and may be omitted for the
`Jenkins1991Entrainment` class.


## plotting/[[eos.py]]

This provides a class, `LinearEos`, for calculating the plume density
according to the linearisation in
[equation 20](../2-numerics/1-equations.html#mjx-eqn-eqlin-eos). The
constructor for this class
```python
LinearEos(ref_density, beta_T, beta_S, T_ref, S_ref)
```
where the arguments are the reference density, thermal contraction
coefficient, haline contraction coefficient, reference temperature,
and reference salinity. All quantities should be scaled to be in
nondimensional units. The resulting object (called, say, `eos`) is
called like a function:
```python
eos(T, S)
```
where `T` is the plume temperature and `S` is the plume salinity.


## plotting/[[melt.py]]

This module provides classes which can calculate the melt rate of the
ice shelf and melt contributions to the plume's salinity and heat
equations. The first class is `Dallaston2015Melt` which calculates the
melt rate in the same manner as the [[dallaston2015_melt(type)]]
derived type. It has constructor
```python
Dallaston2015Melt(beta, epsilon_g, epsilon_m)
```
where `beta` is the inverse Stefan number, `epsilon_g` is the ratio of
subglacial flux to entrained flux, and `epsilon_g` is the ratio of
subshelf melt and entrained flux. The other class is `OneEquationMelt`,
which calculates the melt rate in the manner of the
[[one_equation_melt(type)]] derived type. It has constructor
```python
OneEquationMelt(coef1, coef2, fresh_sal=0., melt_temp=0.)
```
where `coef1` and `coef2` correspond to \(\zeta_1\) and \(\zeta_2\) in
[equation 18](../2-numerics/1-equations.html#mjx-eqn-eqtherm-trans-nondim),
`fresh_sal` is the salinity value which is designated as "fresh", and
`melt_temp` is the temperature at which melting occurs (scale analysis
shows that it often makes sense for these to be non-zero).

For both classes, the melt rate is calculated by calling the melt
object (here named `m`) like a function. Additionally, there are
methods `thermal_forcing` and `saline_forcing` to calculate the
contribution of melting to the plume heat and salinity equations,
respectively. These can be called as follows:
```python
m(U, p, T, S, D)
m.thermal_forcing(U, p, T, S, D)
m.saline_forcing(U, p, T, S, D)
```
In each case, `U` refers to the plume vector velocity, `p` to the
pressure at the base of the ice shelf, `T` to the plume temperature,
`S` to the plume salinity, and `D` to the plume thickness.


## plotting/[[viscosity.py]]

Finally, this module provides two classes which can calculate the
viscosity of the glacier: `NewtonianViscosity` and `GlensLaw`. The
first of these is a trivial implementation with a constructor
```python
NewtonianViscosity(value=1.0)
```
Making calls to objects of this type will return an array filled with
the value with which the viscosity object was initialised. The
second implementation is `GlensLaw`, which treats viscosity as a
power law of strain
[(equation 5)](../2-numerics/1-equations.html#mjx-eqn-eqglens-nd). It
has constructor
```python
GlensLaw(size, lower=0.0, upper=1.0, coef=1.0, index=3)
```
The first 3 arguments have the same meaning as in the constructor for
the [Differentiator class](#differentiate). The argument `coef`
corresponds to \(\xi\), as defined in
[equation 6](../2-numerics/1-equations.html#mjx-eqn-eqglens-coef-nd)
and `index` is the value of \(n\) in the equation for Glen's Law.

Both viscosity classes are called like functions and (if an instance
of them is named `visc`) have the call signature
```python
visc(uvec, temperature=-15.0, time=0)
```
where `uvec` is the vector ice velocity, `temperature` is the
temperature of the ice, and `time` is the time during the simulation
at which the calculation is being performed. The latter two arguments
are not used in either implementation of viscosity.

