#!/usr/bin/env python
from plotting.readers import ShelfPlumeCryosphere
from plotting.layers import compute_layers

import numpy as np
import numpy.ma as ma
import matplotlib
import matplotlib.pyplot as plt

import sys

if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

nz = 1200
conts = 16

x, z, k = compute_layers(cryo, nz)
#plt.plot(cryo.s[-1] - z[:,-1], k[:, -1], label=r'$k$')
#plt.plot(cryo.s[-1] - z[:,-1], 1/np.sqrt(1.05 - cryo.s[-1] + z[:,-1]) 
#         - 1/np.sqrt(1.05), label=r'$k$')
#plt.legend(loc=0)
#plt.show()

plt.figure(figsize=(8,6))
plt.plot(cryo.grid, cryo.s, 'k')
plt.plot(cryo.grid, cryo.b, 'k')
plt.contour(x, z, k, conts, linewidths=0.5, colors='k')
plt.contourf(x, z, k, conts, cmap='PuBu')
plt.xlabel('$x$')
plt.ylabel('$z$')
plt.colorbar(orientation='horizontal', label=r'$k$')
plt.tight_layout()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()
