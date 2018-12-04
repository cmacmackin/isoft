Title: Parameterisations and Boundary Conditions
Author: Chris MacMackin
Date: November 2018 

One of the goals of ISOFT is to allow choices of parameterisations to
easily be changed. This is achieved using the <span
style="font-variant:small-caps;">Strategy</span> pattern
[(Rouson, Xia, and Xu, 2014, Chapter 7)](../6-bibliog.html#Rouson2014),
which provides a common abstract interface to accomplish some task,
with subtypes implementing different strategies to do so. In ISOFT,
the methods in the abstract types were generally given a large number
of arguments, to ensure sufficient information is available for all
potential parameterisations. Parameter and coefficient values can be
specified for each parameterisation when initialising its object.

The only parameterisation for the ice shelf is viscosity. The general
interface is provided by the [[abstract_viscosity(type)]] type. It’s subtypes
are [[newtonian_viscosity(type)]], which returns a `uniform_field` all
cases, and [[glens_law_viscosity(type)]] which calculates the viscosity from
the ice velocity as described in [equation 5](../2-numerics/1-equations.html#mjx-eqn-eqglens-nd). Currently Glen’s law is
only implemented for the 1-D case, as attempting to implement it for
higher dimensions resulted in a compiler bug. A class diagram is provided below.

![A UML class diagram for the viscosity type. Subtypes implement all
inherited abstract methods, but this is not shown for reasons of
space.](|media|/viscosity.png)

The plume contains a few parameterisations. The subtypes of
[[abstract_entrainment(type)]] calculate an entrainment rate for the
plume. These are [[jenkins1991_entrainment(type)]] and
[[kochergin1987_entrainment(type)]], implementing
[equations 14](../2-numerics/1-equations.html#mjx-eqn-eqentrainment-jenkins-nd)
and
[15](../2-numerics/1-equations.html#mjx-eqn-eqentrainment-koch-nd),
respectively.  The [[abstract_melt_relationship(type)]] provides an
interface for calculating the melt rate of the ice, in addition to the
heat and salt fluxes resulting from melting. The one equation
approximation of
[equation 17](../2-numerics/1-equations.html#mjx-eqn-eqmelt-nondim)
was implemented as [[one_equation_melt(type)]]. A variation of this
was implemented as [[ave_one_equation_melt(type)]], implementing the
horizontally-averaged version of the one equation formulation found in
equation XXXX. The subtype [[dallaston2015_melt(type)]] provides a way
to convert from the scaling choices used by
[Dallaston, Hewitt, and Wells (2015)](../bibliog.html#Dallaston2015)
to those used in ISOFT, which was useful for writing unit
tests. Finally, the abstract type [[equation_of_state(type)]] sets out
the interface for calculating the density of water from salinity and
temperature. Subtype [[linear_eos(type)]] implements the linearised
equation of state in equation XXXX. The related
[[ave_linear_eos(type)]] provides additional methods methods for
calculating \(\overline{\rho}\) and \(\tilde{\rho}\), as defined in
equations XXX and XXX. Last, the subtype [[prescribed_eos(type)]]
calculates the density assuming no dependence on temperature and using
a prescribed salinity profile; this is also useful in unit
tests. Class diagrams for each set of derived types are provided
below.

![A UML class diagram for the entrainment type. Subtypes implement the
inherited abstract method, but this is not shown for reasons of
space.](|media|/entrainment.png)

![A UML class diagram for the melt type. Subtypes implement the
inherited abstract method, but this is not shown for reasons of
space.](|media|/melt.png)

![A UML class diagram for the equation of state type. Subtypes implement
all inherited abstract methods, but this is not shown for reasons of
space.](|media|/eos.png)

A similar approach was taken for boundary conditions and ambient ocean
properties. The types [[glacier_boundary(type)]]
and [[plume_boundary(type)]] (class diagrams below) provide interfaces
for identifying the types of boundary conditions at different locations
and determining the appropriate values. The default implementations
effectively do not specify boundary conditions and the methods must be
overridden to be useful. The interface provided by [[plume_boundary(type)]] is
quite different from that provided by [[glacier_boundary(type)]]. The latter should
ideally be refactored to be closer to the more usable interface provided
by the former. The subtypes for [[glacier_boundary(type)]] are
[[dallaston2015_glacier_boundary(type)]], which provides time-independent
boundary conditions like those described in equations XXXX and XXXX, and
[[seasonal_glacier_boundary(type)]], which modifies these conditions according
to equations  or  to allow seasonal oscillations in ice flux.

![A UML class diagram for the glacier boundary type. Subtypes implement
all inherited abstract methods, but this is not shown for reasons of
space.](|media|/icebound.png)

![A UML class diagram for the plume boundary type. Subtypes implement
all inherited abstract methods, but this is not shown for reasons of
space.](|media|/plumebound.png)

The first subtype of [[plume_boundary(type)]] is
[[simple_plume_boundary(type)]], which implements boundary conditions
of the type described in equations XXXX and XXXX.  Closely related to
this is [[dallaston2015_seasonal_boundary(type)]], which modifies the
boundary conditions according to equation XXXXX. The type which was
ultimately used in all simulations was
[[upstream_plume_boundary(type)]].  This takes a user-provided
function which specifies the inflow value of each plume variable and
then, assuming no diffusion, integrates the plume a small distance
upstream along the current basal draft of the ice shelf using
[[rksuite_90]]
[(Brankin and Gladwell, 1994)](../6-bibliog.html#Brankin1994). This
allows the plume solver itself to avoid handling narrow boundary
layers where the plume salinity and temperature change
rapidly. Outflow conditions are again defined according to
equation XXX. Ambient ocean conditions are described according to the
interface defined by the abstract type
[[ambient_conditions(type)]]. As shown in the class diagram below, at
present only one implementation is provided
([[uniform_ambient_conditions(type)]]), specifying constant ambient
salinity and temperature.

![A UML class diagram of the ambient conditions type. The subtype
implements all inherited abstract methods, but this is not shown for
reasons of space.](|media|/ambient.png)
