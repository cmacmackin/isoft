#!/usr/bin/env python
#
#  entrainment.py
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

'''Contains classes for calculating the entrainment rate of the plume
using different parameterisations.
'''

import numpy as np
import calculus

class Jenkins1991Entrainment(object):
    '''A class representing the entrainment formulation used by Jenkins et
    al. (1991).

    coefficient
        The coefficient used in the entrainment calculation
    size
        The number of Chebyshev modes in the field
    lower
        The lower boundary of the field domain
    upper
        The upper boundary of the field domain
    '''

    def __init__(this, coefficient, size, lower=0.0, upper=1.0):
        this.coef = coefficient
        this.diff = calculus.Differentiator(size, lower, upper)

    def __call__(this, U, D, b, rho_diff = None):
        return this.coef * np.linalg.norm(U, axis=-1) * np.abs(this.diff(b))


class Kochergin1987Entrainment(object):
    '''A class representing the entrainment formulation used by Kochergin 
    (1987).

    coefficient
        The coefficient $c_l^2x_0/D_0$ used in the entrainment calculation
    delta
        The ratio between the scale of the plume thickness and that of the 
        ice thickness, $D_0/h_0$
    '''

    def __init__(this, coefficient, delta):
        this.coef = coefficient
        this.delta = delta

    def __call__(this, U, D, b, rho_diff):
        Ri = this.delta*rho_diff*D/np.linalg.norm(U, axis=-1)**2
        Sm = Ri/(0.0725*(Ri + 0.186 - np.sqrt(Ri**2 - 0.316*Ri + 0.0346)))
        return this.coef*np.linalg.norm(U, axis=-1)/Sm * np.sqrt(1 + Ri/Sm)
