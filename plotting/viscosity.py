#!/usr/bin/env python
#
#  viscosity.py
#  This file is part of ISOFT.
#  
#  Copyright 2018 Chris MacMackin <cmacmackin@gmail.com>
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

'''Contains classes for calculating the viscosity of the ice
using different parameterisations.
'''

import numpy as np
import calculus

class NewtonianViscosity(object):
    '''A class representing Newtonian viscosity.

    value
        The viscosity value of the ice
    '''

    def __init__(this, value=1.0):
        this.val = value

    def __call__(this, uvec, temperature=-15.0, time=0):
        return this.val*np.ones_like(uvec[0,...])


class GlensLaw(object):
    '''A class representing Newtonian viscosity.

    value
        The viscosity value of the ice
    '''

    def __init__(this, size, lower=0.0, upper=1.0, coef=1.0, index=3,):
        this.diff = calculus.Differentiator(size, lower, upper)
        this.coef = coef
        this.index = float(index)

    def __call__(this, uvec, temperature=-15.0, time=0):
        if (uvec.ndim > 2):
            raise NotImplementedError('GlensLaw only implemented for 1-D '
                                      'velocity field.')
        print uvec[0,:]
        return 0.5*this.coef*abs(this.diff(uvec[0,:]))**(1./this.index - 1.)

