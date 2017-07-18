#!/usr/bin/env python

from plotting.readers import ShelfPlumeCryosphere
import numpy as np
import matplotlib.pyplot as plt

p050 = ShelfPlumeCryosphere('050.h5')
p100 = ShelfPlumeCryosphere('100.h5')
p200 = ShelfPlumeCryosphere('200.h5')
p400 = ShelfPlumeCryosphere('400.h5')
p800 = ShelfPlumeCryosphere('800.h5')

maxpoints = p800.grid.size

plumes = [p050, p100, p200, p400]
resolution = np.array([50,100,200,400])
quad_error = np.empty(resolution.size)
max_error = np.empty(resolution.size)

for i, p in enumerate(plumes):
    stride = maxpoints // p.grid.size + 1
    assert np.all(p.grid == p800.grid[::stride])
    state = np.concatenate((p.D, p.U, p.T, p.S))
    truth = np.concatenate((p800.D[::stride], p800.U[::stride], 
                            p800.T[::stride], p800.S[::stride]))
    error = state - truth
    max_error[i] = np.max(np.abs(error))
    quad_error[i] = np.linalg.norm(error)

plt.loglog(resolution, quad_error, label = "Error (Euclidean Norm)", basex=2)
plt.loglog(resolution, max_error, label = "Error (Maximum)", basex=2)
plt.xlabel("Grid Points")
plt.legend()
plt.savefig("error.pdf")
