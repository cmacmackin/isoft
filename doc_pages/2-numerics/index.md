Title: Numerical Methods
Author: Chris MacMackin
Date: December 2018 

ISOFT simulates the evolution of ice shelves and meltwater plumes
beneath them. A cartoon diagram of such a can be found below. The ice
flow has a vertically integrated velocity \(\vec{u}\), with longitudinal and
transverse components \(u\) and \(v\), respectively. \(h\) is the thickness of the
ice shelf, while \(b\) is the depth of its lower surface below sea
level. Subglacial discharge at the grounding line, with volume flux \(Q_g\), feeds a plume of thickness \(D\) flowing underneath the ice shelf with
vertically integrated velocity \(\vec{U}\). This velocity also has
longitudinal, \(U\), and transverse, \(V\), components. The plume has a
temperature, \(T\), and salinity, \(S\), which drive melting \(m\) at the base of
the ice shelf. The plume is further fed by turbulent entrainment, \(e\),
of the ambient ocean water. This water has its own temperature and
salinity: \(T_a\) and \(S_a\), respectively.

![A cartoon diagram of an ice shelf and meltwater plume.](|media|/iceshelf.svg)

This section of the documentation describes the mathematics behind
ISOFT. First, the equations describing ice shelf and plume behaviour
are provided. The solvers used for the shelf and plume components are
described in turn. Finally, an overview of a benchmarking problem is
provided.
