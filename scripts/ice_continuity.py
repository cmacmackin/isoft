#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
from plotting.calculus import Differentiator
from plotting.melt import OneEquationMelt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import sys


if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0200.h5')

x = cryo.grid
diff = Differentiator(x.size, x[-1], x[0])
D = cryo.D
Uvec = cryo.Uvec
U = cryo.U
T = cryo.T
S = cryo.S
h = cryo.h
u = cryo.u
m = OneEquationMelt(0.0182, 4.86e-4, -2035.3, -54.92)
lambd = cryo.lambd

plt.figure()
u_dh = u*diff(h)
h_du = h*diff(u)
m = -cryo.lambd*m(Uvec, cryo.b, T, S, D)
dh_dt = m - u_dh - h_du
plt.plot(x, dh_dt, label=r'$h_t$')
plt.plot(x, u_dh, label=r'$uh_x$')
plt.plot(x, h_du, label=r'$u_x h$')
plt.plot(x, m, label=r'$-\lambda m$')
plt.xlabel(r'$x$')
plt.legend(loc=0)
plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
