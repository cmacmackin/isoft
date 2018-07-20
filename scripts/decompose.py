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

#c1 = ShelfPlumeCryosphere('nu+1.h5')
#c31 = ShelfPlumeCryosphere('nu+5.h5')
#c61 = ShelfPlumeCryosphere('nu+10.h5')
#x = c1.grid
#diff = Differentiator(x.size, x[-1], x[0])
#print diff(c1.U)
#print diff(c31.U)
#print diff(c61.U)
#plt.plot(c1.grid, c1.U, label='1')
#plt.plot(c31.grid, c31.U, label='31')
#plt.plot(c61.grid, c61.U, label='61')
#plt.legend()
#plt.show()
#quit()

if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

x = cryo.grid
diff = Differentiator(x.size, x[-1], x[0])
m = OneEquationMelt(0.0182, 4.86e-4, -2035.3, -54.92)
e = Jenkins1991Entrainment(1.0, x.size, x[-1], x[0])
eos = LinearEos(3.05e5, 1.409e-6, 1.336e-5, 0., 0.)
D = cryo.D
Uvec = cryo.Uvec
U = cryo.U
T = cryo.T
S = cryo.S
b = cryo.b
ent = e(Uvec, D, b)
mu = cryo.mu
nu = cryo.nu
delta = cryo.delta
T_a = 0.0
S_a = 0.0
rho = eos(T, S)
rho_a = eos(T_a, S_a)

plt.plot(x, D, label='D')
plt.plot(x, U, label='U')
plt.plot(x, T, label='T')
plt.plot(x, S, label='S')
plt.plot(x, cryo.b, label='b')
plt.plot(x, diff(cryo.b), label='b_x')
#plt.plot(x, np.sqrt(D*(rho_a-rho)), label='sqrt(B)')
#plt.plot(x, Udal, label='U, Dal')
plt.legend()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
plt.show()

plt.plot(x, diff(D*U), label=r'$(DU)_x$')
plt.plot(x, ent, label=r'$e$')
plt.plot(x, m(Uvec, b, T, S, D), label=r'$m$')
plt.plot(x, diff(D*U) - ent - m(Uvec, b, T, S, D), label='Sum')
#plt.plot(x, U*diff(D), label=r'$D_xU$')
#plt.plot(x, D*diff(U), label=r'$DU_x$')
plt.legend()
plt.show()

plt.plot(x, diff(D*U**2), label=r'$(DU^2)_x$')
plt.plot(x, D*(rho_a - rho)*diff(b - delta*D), 
         label=r'$D(\rho_a - \rho)(b_x - \delta D_x)$')
plt.plot(x, nu*diff(D*diff(U)), label=r'$\nu (DU_{x})_{x}$')
plt.plot(x, -mu*np.abs(U)*U, label=r'$-\mu |U|U$')
plt.plot(x, delta/2*D**2*diff(rho), label=r'$\frac{\delta D^2}{2}\rho_x$')
plt.plot(x, diff(D*U**2) - D*(rho_a - rho)*diff(b - delta*D) - 
         nu*diff(D*diff(U)) + mu*abs(U)*U - delta/2*D**2*diff(rho),
         label='Sum')
#plt.plot(x, U**2*diff(D), label=r'$D_xU^2$')
#plt.plot(x, 2*D*U*diff(U), label=r'$2DUU_x$')
plt.legend()
plt.show()

plt.plot(x, diff(D*U*T), label='$(DUT)_x$')
plt.plot(x, ent*T_a, label='$eT_a$')
plt.plot(x, nu*diff(D*diff(T)), label=r'$\nu (DT_{x})_{x}$')
plt.plot(x, -m.thermal_forcing(Uvec, b, T, S, D), label='$-\gamma_T (T-T_m)$')
plt.plot(x, diff(D*U*T) - ent*T_a - nu*diff(D*diff(T)) + 
         m.thermal_forcing(Uvec, b, T, S, D), label='Sum')
#plt.plot(x, T*U*diff(D), label=r'$D_xUT$')
#plt.plot(x, T*D*diff(U), label=r'$DU_xT$')
#plt.plot(x, D*U*diff(T), label=r'$DUT_x$')
plt.legend()
plt.show()

plt.plot(x, diff(D*U*S), label='$(DUS)_x$')
plt.plot(x, ent*S_a, label='$eS_a$')
plt.plot(x, nu*diff(D*diff(S)), label=r'$\nu (DS_{x})_{x}$')
plt.plot(x, -m.saline_forcing(Uvec, b, T, S, D), label='$-\gamma_S (S-S_m)$')
plt.plot(x, diff(D*U*S) - ent*S_a - nu*diff(D*diff(S)) + 
         m.saline_forcing(Uvec, b, T, S, D), label='Sum')
#plt.plot(x, S*U*diff(D), label=r'$D_xUS$')
#plt.plot(x, S*D*diff(U), label=r'$DU_xS$')
#plt.plot(x, D*U*diff(S), label=r'$DUS_x$')
plt.legend()
plt.show()

