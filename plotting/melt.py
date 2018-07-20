#!/usr/bin/env python
#
#  melt.py
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

'''Contains classes for calculating the melt rate of the ice shelf
using different parameterisations.
'''

import numpy as np

class Dallaston2015Melt(object):
    '''A class representing the melt formulation used by Dallaston et
    al. (2015).

    beta
        Inverse Stefan number
    epsilon_g
        subglacial flux/entrained flux
    epsilon_m
        subshelf melt/entrained flux
    '''

    def __init__(self, beta, epsilon_g, epsilon_m):
        self.beta = beta
        self.epsilon_g = epsilon_g
        self.epsilon_m = epsilon_m

    def __call__(self, U, p, T, S, D):
        return np.linalg.norm(U, axis=-1)*(1 - self.epsilon_m/self.beta*T)

    def thermal_forcing(self, U, p, T, S, D):
        return self(U, p, T, S, D)*(self.beta + 1)

    def saline_forcing(self, U, p, T, S, D):
        return np.zeros(T.size)


class OneEquationMelt(object):
    '''A class representing a melt formulation similar to that used by
    Dallaston et al. (2015), but including the terms they removed
    after scale analysis.

    coef1
        The factor $Gamma_tx_0/D_0$
    coef2
        The factor $c_0T_0/L$
    fresh_sal
        The numerical salinity value given to fresh water
    melt_temp
        The numerical temperature value at which melting occurs
    '''

    def __init__(self, coef1, coef2, fresh_sal=0., melt_temp=0.):
        self.coef1 = coef1
        self.coef2 = coef2
        self.fresh_sal = fresh_sal
        self.melt_temp = melt_temp

    def forcing_value(self, U, p, T, S, D):
        return self.coef1*np.linalg.norm(U, axis=-1)*(T - self.melt_temp)        

    def __call__(self, U, p, T, S, D):
        return self.coef2*self.forcing_value(U, p, T, S, D)

    def thermal_forcing(self, U, p, T, S, D):
        return (1. - self.coef2*self.melt_temp)*self.forcing_value(U, p, T, S, D)

    def saline_forcing(self, U, p, T, S, D):
        return -self(U, p, T, S, D)*self.fresh_sal
