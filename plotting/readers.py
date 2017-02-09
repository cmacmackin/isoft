#!/usr/bin/env python
#
#  readers.py
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

'''
Contains classes for reading simulation data from HDF5 files.
'''

import h5py

class Glacier(object):
    '''A class which reads glacier data from HDF5 output. It presents
    the data in a convenient manner for plotting.

    hdf_obj
        An HDF5 file or group.
    hdf_group
        The name of the group in `hdf_obj` which contains the ice shelf 
        data.
    '''

    def __init__(self, hdf_obj, hdf_group):
        self.groupname = hdf_group
        self.data = hdf_obj[hdf_group]
    
    @property
    def grid(self):
        return self.data['thickness'].attrs['gridpoints_dim1'][...]

    @property
    def h(self):
        return self.data['thickness'][...]

    @property
    def u():
        return self.data['velocity'][0,...]

    @property
    def v(this):
        return self.data['velocity'][1,...]

    @property
    def glacier_type(self):
        return self.attrs['glacier_type']


