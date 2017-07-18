#!/usr/bin/env python
#
#  calculus.py
#  This file is part of ISOFT.
#  
#  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

import numpy as np

cheb_cache = dict()

def cheb(size, lower=0.0, upper=1.0):
    '''
    Computes a Chebyshev differentiation matrix of the appropriate size.
    Returns a 1D Numpy array of Chebyshev collocation points and a 2D
    differentiation matrix.
    '''
    if (size, lower, upper) in cheb_cache:
        return cheb_cache[(size, lower, upper)]

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
    D = np.empty((size+1, size+1))
    for i in range(size + 1):
        if i != 0 and i != size: D[i,i] = -x[i]/(2.0*(1.0-x[i]**2))
        for j in range(i) + range(i+1, size+1):
            D[i,j] = c(i)/c(j)/(x[i] - x[j])
    c = (2.0*float(size**2) + 1.0)/6.0
    D[0,0] = c
    D[size,size] = -c
    # Adjust for use on [0,1] interval
    factor = (upper - lower)/2.0
    x = factor*(1 + x) + lower
    D = 2.0*D/(upper - lower)
    cheb_cache[(size, lower, upper)] = (x,D)
    return (x,D)


class Differentiator(object):
    """An object which can perform a differentiation operation on a field
    from an ISOFT run.
    
    size
        The number of Chebyshev modes in the field
    lower
        The lower boundary of the field domain
    upper
        The upper boundary of the field domain
    """
    
    def __init__(this, size, lower=0.0, upper=1.0):
        this.x, this.D = cheb(size - 1, lower, upper)

    def __call__(this, rhs):
        return this.D.dot(rhs)

