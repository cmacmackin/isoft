Title: Representing the Coupled Ice Shelf/Plume
Author: Chris MacMackin
Date: November 2018 

ISOFT uses a large number of *derived types* (equivalent to *classes*
in other object oriented languages) in order to model the full glacial
system. The system as a whole is contained within a [[cryosphere(type]] type,
with methods for saving/loading data to/from the disk as HDF5 files
and for integrating forward in time.  The [[cryosphere(type)]] contains
objects of abstract classes [[basal_surface(type)]] and [[glacier(type)]], the latter
representing a glacier and the former representing whatever is
underneath it. Objects of these types have their own methods for
reading and writing data, integration, and accessing useful
information. Both are general enough to allow ISOFT to model either
floating or grounded ice. Below are class diagrams showing the
inheritance/encapsulation relationships between different derived
types, as well as methods for the [[cryosphere(type)]], [[glacier(type)]], and
[[basal_surface(type)]] types.

![A UML class diagram illustrating the relationships between the key
derived types used in ISOFT. Components and methods of the types are not
shown for reasons of space.](|media|/inheritance_diagram.png)

![UML class diagram for the cryosphere derived type, providing a simplified description of the methods associated with it.](|media|/cryo.png)

![UML class diagrams for the basal surface derived type, providing a simplified description of the methods associated with it.](|media|/basal.png)

![UML class diagrams for the glacier derived type, providing a simplified description of the methods associated with it.](|media|/glacier.png)

The only concrete existing implementation of the [[glacier(type)]] class is
the [[ice_shelf(type)]] type. As its name suggests, this models the behaviour
of an ice shelf. While the implementation of the continuity equation
is agnostic towards whether the model is 1-D or 2-D, at present the
ice momentum equation is explicitly 1-D. Ideally this will be fixed in
the future. The [[ice_shelf(type)]] type may optionally feature a Lagrangian
tracer, assumed to indicate the age of the ice as would be measured
from internal reflectors.  There is stub for a grounded [[ice_sheet(type)]]
type, but its methods have not been implemented.

A few implementations of the [[basal_surface(type)]] class are available. The
most commonly used of these is the [[plume(type)]] type, modelling the 1-D
subglacial plume. In principle this can model a second
velocity component, but such a model is physically unstable. There is
also the [[static_plume(type)]] type, which does not evolve from the state with
which it is initialised or has loaded from saved data. This is useful if
a simulation is to be performed with a fixed melt rate. The
[[asym_plume(type)]] provides an implementation of a
horizontally-integrated model. Various parameters describing the
transverse profile of this plume are provided through the associated
[[plume_shape(type)]] derived type. Finally, a [[ground(type)]] type exists as a stub,
which could be fully implemented in future to provide a basal surface
with frictional information for a grounded ice sheet.

![A UML object diagram displaying the state of a typical simulation
immediately after it has been initialised. Not all components are shown
in each object for reasons of space.](|media|/object_tree.png)

All of these implementations contain *field types* (described in the
[next section](./discretisation.html)) for each variable, describing
their state. They also contain objects representing the boundary
conditions and choices of parameterisations, described in more detail
in the
[Parameterisations and Boundary Conditions](./parameterisations.md)
section. This is illustrated in the preceding figure, showing the
state of the objects at the beginning of a representative
simulation. The [[cryosphere(type)]] class is a <span
style="font-variant:small-caps;">Puppeteer</span> pattern which, as
described by
[Rouson, Xia, and Xu (2014)](../6-bibliog.html#Rouson2014),
coordinates interactions between various other classes
([[glacier(type)]] and [[basal_surface(type)]], in this case) without
them needing to directly interact with each other. Thus,
interdependencies between different modules of code are
simplified. For each time-step, the following sequence of steps
occurs, as illustrated in the UML sequence diagram which follows:

1.  The [[cryosphere(type)]] first gets the ice thickness from the [[glacier(type)]].

2.  This information is sent to the [[basal_surface(type)]] object, with which
    it can solve for its state at the current time using the QLM solver.

3.  The [[cryosphere(type)]] gets the current melt rate and/or friction
    parameters from the [[basal_surface(type)]].

4.  This information is sent to the [[glacier(type)]] object, where it is used
    to integrate its state forward in time with NITSOL.


![A UML sequence diagram illustrating the operations involved in each
time step of ISOFT. Some of the argument lists have been simplified for
reasons of space.](|media|/isoft_sequence.png)
