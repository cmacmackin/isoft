Title: Code Design
Author: Chris MacMackin
Date: November 2018 

This section of the documentation provides an overview of the
structure of ISOFT and an explanation of the design choices made. It
is limited to a discussion of code architecture, with numerical
methods already having been described in
[the previous section](../2-numerics/). A number of design patterns
were consciously used when developing ISOFT. The names of these are
indicated in the text by small-caps.  An understanding of object
oriented programming techniques in general and object oriented
features in Fortran, in particular, will be useful when reading these
notes; see, for example,
[Rouson, Xia, and Xu (2014)](../6-bibliog.html#Rouson2014). Familiarity
with Universal Modelling Language (UML) diagrams will also be helpful
[(Rouson, Xia, and Xu, 2014)](../6-bibliog.html#Rouson2014). Note that,
in the UML diagrams in this section, names of derived types are shown
in camel case, as is the convention for class names in most object
oriented languages. However, as Fortran is case-insensitive, the
decision was made to use underscore separation within the code.
