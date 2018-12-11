#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import sys


if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

x = cryo.grid
D = cryo.D
Uvec = cryo.Uvec
U = cryo.U
T = cryo.T
S = cryo.S
b = cryo.b

r = 1.12
gamma = 4
lambd = 1e2
X = 1.0/(lambd*2.02188465743*4.3316e-4)
velocity = np.sqrt(1.0 + gamma/4.0*X - gamma/4.0*X*(1.0 - x/X)**2)
thickness = (1.0-x/X)/velocity
depth = -thickness/r

fig = plt.figure(figsize=(5,5))
gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1])
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1], sharex=ax1)

ax1.plot(x, D, label='$D$')
ax1.plot(x, U, label='$U$')
ax1.plot(x, cryo.b, label='$b$')
ax1.plot(x, depth, '--', label=r'$\hat{b}$')
ax1.plot(x, cryo.u, label=r'$u$')
ax1.get_xaxis().set_visible(False)
ax1.set_ylim(-1, 3.5)
ax1.get_yticklabels()[0].set_visible(False)
ax1.legend(prop={'size': 'medium'}, ncol=3, columnspacing=1.0, loc=2)

ax2.plot(x, T, label='$T$')
ax2.plot(x, S, label='$S$')
ax2.set_xlabel('$x$')
ax2.legend(prop={'size': 'medium'}, loc=3)

fig.tight_layout()
fig.subplots_adjust(hspace=0)

if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()

