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

'''Contains classes for calculating density of water.
'''

import numpy as np
import calculus

class LinearEos(object):
    '''A class representing a linearised equation of state. It uses the
    equation:

    density = ref_density*[1 - beta_T*(T-T_ref) + beta_S*(S-S_ref)
    '''

    def __init__(this, ref_density, beta_T, beta_S, T_ref, S_ref):
        this.rd = ref_density
        this.bT = beta_T
        this.bS = beta_S
        this.Tr = T_ref
        this.Sr = S_ref

    def __call__(this, T, S):
        return this.rd*(1 - this.bT*(T - this.Tr) + this.bS*(S - this.Sr))
