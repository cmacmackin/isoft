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
        return self.data[self.data['thickness'].attrs['gridpoints_dim1']][...]

    @property
    def h(self):
        return self.data['thickness'][...]

    @property
    def u(self):
        return self.data['velocity'][0,...]

    @property
    def v(self):
        return self.data['velocity'][1,...]

    @property
    def glacier_type(self):
        return self.data.attrs['glacier_type']

    @property
    def lambd(self):
        return self.data.attrs['lambda']

    @property
    def chi(self):
        return self.data.attrs['chi']


class Plume(object):
    '''A class which reads Plume data from HDF5 output. It presents
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
        return self.data[self.data['thickness'].attrs['gridpoints_dim1']][...]

    @property
    def D(self):
        return self.data['thickness'][...]

    @property
    def U(self):
        return self.data['velocity'][0,...]

    @property
    def V(self):
        return self.data['velocity'][1,...]

    @property
    def S(self):
        return self.data['salinity'][...]

    @property
    def T(self):
        return self.data['temperature'][...]

    @property
    def basal_surface_type(self):
        return self.data.attrs['basal_type']

    @property
    def delta(self):
        return self.data.attrs['delta']

    @property
    def mu(self):
        return self.data.attrs['mu']

    @property
    def nu(self):
        return self.data.attrs['nu']

    @property
    def r(self):
        return self.data.attrs['r_val']


class ShelfPlumeCryosphere(object):
    '''An abstract class which represents a Cryosphere object containing
    at least one glacier and one basal surface.
    
    hdf_file
        The HDF5 file containing the ice shelf/plume data.
    '''
    def __init__(self, hdf_file):
        self.filename = hdf_file
        self.data = h5py.File(hdf_file, 'r')
        self.shelf = Glacier(self.data, 'glacier')
        self.plume = Plume(self.data, 'basal_surface')

    @property
    def isoft_version(self):
        return self.attrs['isoft_version']

    @property
    def compilation_time(self):
        return self.attrs['binary_compilation_time']

    @property
    def output_time(self):
        return self.attrs['data_output_time']

    @property
    def time(self):
        return self.attrs['simulation_time']

    @property
    def grid(self):
        return self.shelf.grid

    @property
    def h(self):
        return self.shelf.h

    @property
    def u(self):
        return self.shelf.u

    @property
    def v(self):
        return self.shelf.v

    @property
    def glacier_type(self):
        return self.shelf.glacier_type
    
    @property
    def lambd(self):
        return self.shelf.lambd

    @property
    def chi(self):
        return self.shelf.lambd

    @property
    def s(self):
        return self.h*(1.0 - 1.0/self.r)

    @property
    def b(self):
        return -self.h/self.r

    @property
    def D(self):
        return self.plume.D

    @property
    def U(self):
        return self.plume.U

    @property
    def V(self):
        return self.plume.V

    @property
    def S(self):
        return self.plume.S

    @property
    def T(self):
        return self.plume.T

    @property
    def basal_surface_type(self):
        return self.plume.basal_surface_type

    @property
    def delta(self):
        return self.plume.delta

    def mu(self):
        return self.plume.mu

    @property
    def nu(self):
        return self.plume.nu

    @property
    def r(self):
        return self.plume.r
