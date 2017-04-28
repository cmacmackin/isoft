#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
from plotting.melt import Dallaston2015Melt
from plotting.entrainment import Jenkins1991Entrainment
import numpy as np
import matplotlib.pyplot as plt

m = Dallaston2015Melt(0.024, 1.1e-3, 6.9e-4)
e = Jenkins1991Entrainment(1.0)

for i in range(6):
    cryo = ShelfPlumeCryosphere('isoft-{0:04d}.h5'.format(i))
    plt.plot(cryo.grid, cryo.D, label='D')
    plt.plot(cryo.grid, cryo.U, label='U')
    plt.plot(cryo.grid, cryo.S, label='S')
    plt.plot(cryo.grid, cryo.T, label='T')
    plt.plot(cryo.grid, e(cryo.U, cryo.D, -cryo.h/cryo.r), 
             label='e')
    plt.plot(cryo.grid, m(np.expand_dims(cryo.U, axis=1), 0, cryo.T, 
                          cryo.S, cryo.D), label='m')
    plt.legend()
    plt.show()

