#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.animation as animation
import matplotlib.pyplot as plt


def cheb(size):
    '''
    Computes a Chebyshev differentiation matrix of the appropriate size.
    Returns a 1D Numpy array of Chebyshev collocation points and a 2D
    differentiation matrix.
    '''
    def c(i):
        if i%2 == 0:
            f = 1.0
        else:
            f = -1.0
        if i == 0 or i == size:
            return f*2.0
        else:
            return f*1.0
    
    if size < 1:
        raise Exception('Must have more than one Chebyshev node.')
    #~ x = 0.5*np.cos(np.linspace(0.0, 1.0, size+1) * np.pi) + 0.5
    x = np.cos(np.linspace(0.0, 1.0, size+1) * np.pi)
    x = x + 1.0
    return x

Writer = animation.writers['ffmpeg']
writer = animation.FFMpegWriter(fps=15, bitrate=1800)    
fig = plt.figure()
l1, = plt.plot([], [], label='D')
l2, = plt.plot([], [], label='U')
l3, = plt.plot([], [], label='S')
l4, = plt.plot([], [], label='T')
plt.xlim(0.0,2.0)
plt.ylim(-0.1,5.0)
plt.xlabel(r'$x$')
plt.legend()
r = 1.12

xvals = cheb(np.loadtxt('plume+0.dat').shape[0] - 1)


def init():
    l1.set_data([],[])
    l2.set_data([],[])
    l3.set_data([],[])
    l4.set_data([],[])
    return [l1, l2, l3, l4]

def update_line(num, li1, li2, li3, li4):
    data = np.loadtxt('plume+{}.dat'.format(num))
    li1.set_data(xvals, data[:,0])
    li2.set_data(xvals, data[:,1])
    li3.set_data(xvals, data[:,2])
    li4.set_data(xvals, data[:,3])
    return [li1, li2, li3, li4]

ani = animation.FuncAnimation(fig, update_line, 399,
                              init_func=init, interval=50, blit=True, 
                              fargs=(l1, l2, l3, l4))

ani.save('plume.mp4', writer=writer)

