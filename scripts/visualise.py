#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
from plotting.melt import Dallaston2015Melt
from plotting.entrainment import Jenkins1991Entrainment
import numpy as np
import matplotlib.pyplot as plt

cryo = ShelfPlumeCryosphere('isoft-0000.h5')
grid = cryo.grid
m = Dallaston2015Melt(0.024, 1.1e-3, 6.9e-4)
e = Jenkins1991Entrainment(1.0, grid.size, grid[-1], grid[0])

try:
    #raise Exception()
    cryo = ShelfPlumeCryosphere('isoft_termination_dump.h5')
    plt.plot(cryo.grid, cryo.h, label='h')
    plt.plot(cryo.grid, cryo.u, label='u')
    plt.plot(cryo.grid, cryo.D, label='D')
    plt.plot(cryo.grid, cryo.U, label='U')
    plt.plot(cryo.grid, cryo.T, label='T')
    plt.plot(cryo.grid, cryo.S, label='S')
    plt.legend()
    plt.show()
except:
    pass


#cryo = ShelfPlumeCryosphere('isoft_termination_dump.h5')
#plt.plot(cryo.grid, cryo.h, label='h')
#plt.plot(cryo.grid, cryo.u, label='u')
#plt.plot(cryo.grid, cryo.D, label='D')
#plt.plot(cryo.grid, cryo.U, label='U')
#plt.plot(cryo.grid, cryo.T, label='T')
#plt.plot(cryo.grid, cryo.S, label='S')
#plt.legend()
#plt.show()

for i in range(120):
    cryo = ShelfPlumeCryosphere('isoft-{0:04d}.h5'.format(i))
    plt.cla()
    plt.plot(cryo.grid, cryo.h, label='h')
    plt.plot(cryo.grid, cryo.u, label='u')
    plt.plot(cryo.grid, cryo.D, label='D')
    plt.plot(cryo.grid, cryo.U, label='U')
    plt.plot(cryo.grid, cryo.S, label='S')
    plt.plot(cryo.grid, cryo.T, label='T')
    plt.legend()
    plt.pause(0.2)

plt.cla()
plt.plot(cryo.grid, cryo.h, label='h')
plt.plot(cryo.grid, cryo.u, label='u')
plt.plot(cryo.grid, cryo.D, label='D')
plt.plot(cryo.grid, cryo.U, label='U')
plt.plot(cryo.grid, cryo.S, label='S')
plt.plot(cryo.grid, cryo.T, label='T')
plt.show()
