#!/usr/bin/env python

import sys
from plotting.readers import ShelfPlumeCryosphere
from plotting.dimensionless_nums import froude
from plotting.eos import LinearEos
import matplotlib.pyplot as plt
from numpy import sqrt

eos = LinearEos(3.05e5, 1.409e-6, 1.336e-5, 0., 0.)
T_a = 0.0
S_a = 0.0
rho_a = eos(T_a, S_a)
U_0 = 0.196
D_0 = 43.2
rho_0 = 1030.
rho_s = 3.38e-3
g = 9.8
coef = U_0/sqrt(g*rho_s*D_0/rho_0)


if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

x = cryo.grid
rho = eos(cryo.T, cryo.S)
plt.figure(figsize=(5,5))
plt.plot(x, cryo.D, label='$D$')
plt.plot(x, cryo.U, label='$U$')
plt.plot(x, froude(coef, cryo.Uvec, rho_a - rho, cryo.D), label=r'$\rm Fr$')
plt.legend(prop={'size': 'medium'})

plt.xlabel('$x$')
plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
