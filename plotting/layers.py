#!/usr/bin/env python
#
#  layers.py
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

import numpy as np
import numpy.ma as ma

def compute_layers(shelf, vertical_resolution=300):
    """Computes the internal layers or age field for an ice shelf from
    the Taylor coefficients.

    """
    xc = shelf.grid
    zc = np.linspace(np.min(shelf.b), np.max(shelf.s), vertical_resolution)
    
    xx, zz = np.meshgrid(xc, zc)
    shelf_domain = np.logical_or(np.greater(zz, shelf.s), np.less(zz, shelf.b))
    x = ma.array(xx, mask=shelf_domain, copy=False)
    z = ma.array(zz, mask=shelf_domain, copy=False)
    
    kappa = shelf.kappa
    
    # This isn't really the most efficient way to calculate the Taylor
    # series, but array broadcasting was giving me a headache.
    k = np.zeros_like(z)
    for i in range(1, 1+kappa.shape[0]):
        k += kappa[i-1] * (shelf.s - z)**i
    return x, z, k
