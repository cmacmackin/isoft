Title: Provided Plotting Scripts
Author: Chris MacMackin
Date: November 2018 

The Python library described in the previous sections allows you to
produce plots with tools such as
[matplotlib](https://matplotlib.org/). This can be useful when
creating bespoke visualisations or undertaking detailed analysis of
ISOFT output. However, there are many simple types of plots which you
will want to create regularly. For these, a selection of programs have
been provided in the directory `scripts/`. Where it is necessary to
specify a scaling or parameterisation choice, the same ones are used
as in the [example program](../1-running/2-writing.html). Unless
otherwise indicated, all scripts have the following call-signature:
```
./script-name.py [ infile [ outfile ] ]
```
Here, `infile` refers to the name of the HDF5 file the data in which
is to be plotted. If not specified, it defaults to
`isoft-0000.h5`. Conversely, `outfile` is the name of a file in which
the plot should be saved, with the format taken from the file-name
extension. If this is not specified then the script will produce the
plot in an interactive graphical interface.


## scripts/[[decompose.py]]

Displays a series of plots, one for each of the equations describing
the 1-D plume
([equations 8, 9, 11, and 12](../2-numerics/1-equations.html#mjx-eqn-eqplume-cont)). In
each figure, the magnitude of the individual terms in that equation
are plotted against each other, alongside their sum (which should be
0). The script can not save plots to a file and only takes one
argument. The next plot will be produced once you have closed the
previous one.


## scripts/[[display.py]]

Produces a plot of all plume and ice shelf variables. Salinity and
temperature are shown on a seperate set of axes, as these values can
be quite different from those of the other variables.


## scripts/[[ice-continuity.py]]

Produces a plot showing the value of each term in the continuity
equation for the ice shelf
([equation 1](../2-numerics/1-equations.html#mjx-eqn-eqice-height-nd)).


## scripts/[[layers.py]]

Produces a contour plot showing the value of a Lagrangian tracer
(indicating, e.g., the age of ice). Similar to the plot shown on the
[previous page](./layers.html).


## scripts/[[plume-momentum.py]]

Produces a plot showing the value of each term in the transverse
momentum equation for the plume
([equation 9](../2-numerics/1-equations.html#mjx-eqn-eqplume-mom-x)).


## scripts/[[shock.py]]

A plot of the plume velocity, thickness, and Froude number. This can
be useful for diagnosing the presence of a hydrolic shock.
