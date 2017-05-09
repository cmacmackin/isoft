#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
from plotting.melt import Dallaston2015Melt
from plotting.entrainment import Jenkins1991Entrainment
import numpy as np
import matplotlib.pyplot as plt

m = Dallaston2015Melt(0.024, 1.1e-3, 6.9e-4)
e = Jenkins1991Entrainment(1.0)

ss = ShelfPlumeCryosphere('steadystate50.h5')

for i in range(11):
    cryo = ShelfPlumeCryosphere('isoft-{0:04d}.h5'.format(i))
    plt.cla()
    plt.plot(cryo.grid, cryo.h - ss.h, label='h')
    # plt.plot(cryo.grid, cryo.D - ss.D, label='D')
    # plt.plot(cryo.grid, cryo.U - ss.U, label='U')
    # plt.plot(cryo.grid, cryo.S - ss.S, label='S')
    # plt.plot(cryo.grid, cryo.T - ss.T, label='T')
    plt.legend()
    plt.pause(1)

plt.cla()
plt.plot(cryo.grid, cryo.h - ss.h, label='h')
# plt.plot(cryo.grid, cryo.D - ss.D, label='D')
# plt.plot(cryo.grid, cryo.U - ss.U, label='U')
# plt.plot(cryo.grid, cryo.S - ss.S, label='S')
# plt.plot(cryo.grid, cryo.T - ss.T, label='T')
plt.legend()
plt.show()
