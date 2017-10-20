#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
from plotting.melt import OneEquationMelt, Dallaston2015Melt
from plotting.entrainment import Jenkins1991Entrainment
from plotting.eos import LinearEos
from plotting.calculus import Differentiator
from plotting.dimensionless_nums import froude
import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

x = cryo.grid
diff = Differentiator(x.size, x[-1], x[0])
#m = OneEquationMelt(0.0182, 0.0238)
#e = Jenkins1991Entrainment(1.0, x.size, x[-1], x[0])
eos = LinearEos(3.05e5, 7.74e-5, 0.0271, 1., 1.)
D = cryo.D
Uvec = cryo.Uvec*1000
U = cryo.U*1000
T = cryo.T
S = cryo.S
b = cryo.b
#ent = e(Uvec, D, b)
mu = cryo.mu
nu = cryo.nu*1000
delta = cryo.delta
T_a = 1.0
S_a = 1.0
rho = eos(T, S)
rho_a = eos(T_a, S_a)

plt.plot(x, diff(D*U**2), label=r'$(DU^2)_x$')
plt.plot(x, D*(rho_a - rho)*diff(b), 
         label=r'$D(\rho_a - \rho)b_x$')
plt.plot(x, D*(rho_a - rho)*diff(-delta*D), 
         label=r'$-D(\rho_a - \rho)\delta D_x$')
plt.plot(x, nu*diff(D*diff(U)), label=r'$\nu (DU_{x})_{x}$')
plt.plot(x, -mu*np.abs(U)*U, label=r'$-\mu |U|U$')
plt.plot(x, delta/2*D**2*diff(rho), label=r'$\frac{\delta D^2}{2}\rho_x$')
#plt.plot(x, diff(D*U**2) - D*(rho_a - rho)*diff(b - delta*D) - 
#         nu*diff(D*diff(U)) + mu*abs(U)*U - delta/2*D**2*diff(rho),
#         label='Sum')
plt.legend()
plt.xlabel('$x$')
plt.tight_layout()

if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()

