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

    def __init__(this, beta, epsilon_g, epsilon_m):
        this.beta = beta
        this.epsilon_g = epsilon_g
        this.epsilon_m = epsilon_m

    def __call__(this, U, p, T, S, D):
        return np.linalg.norm(U, axis=-1)*(1 - this.epsilon_m/this.beta*T)
        
