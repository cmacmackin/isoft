#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
import matplotlib.pyplot as plt
import sys


if len(sys.argv) > 1:
    cryo = ShelfPlumeCryosphere(sys.argv[1])
else:
    cryo = ShelfPlumeCryosphere('isoft-0000.h5')

x = cryo.grid
D = cryo.D
Uvec = cryo.Uvec
U = cryo.U*1000
T = cryo.T
S = cryo.S
b = cryo.b

plt.plot(x, D, label='$D$')
plt.plot(x, U, label='$U$')
plt.plot(x, T, label='$T$')
plt.plot(x, S, label='$S$')
plt.plot(x, cryo.b, label='$b$')
plt.xlabel('$x$')
plt.tight_layout()
plt.legend()
if len(sys.argv) > 2:
    plt.savefig(sys.argv[2])
else:
    plt.show()

