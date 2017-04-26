!
!  basal_surface.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@physics.ox.ac.uk>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!  
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
!  MA 02110-1301, USA.
!  

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module basal_surface_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract data type to model the ground or ocean below
  ! the glacier.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field
  use nitsol_mod, only: nitsol, dummy_jacv, ddot, dnrm2, iplvl, dnrm2
  use hdf5
  implicit none
  private

  character(len=10), parameter, public :: hdf_type_attr = 'basal_type'

  type, abstract, public :: basal_surface
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! An abstract data type which represents whatever lies below a [[glacier]].
    ! This could be the ground, a plume, or a fully dynamic ocean model.
    ! Methods are available to provide the coupling information between the
    ! [[glacier]] and the basal surface.
    !
  contains
    procedure(get_scalar), deferred    :: basal_melt
      !! Returns the basal melt rate.
    procedure(get_scalar), deferred    :: basal_drag_parameter
      !! Returns a value which may be needed to calculate basal drag,
      !! such as the coefficient of friction.
    procedure(get_real), deferred      :: water_density
      !! Density of the water at the basal surface.
    procedure(setter), deferred        :: update
      !! Sets the state of the basal surface
    procedure(get_i), deferred         :: data_size
      !! Returns the number of elements in the basal surface's state
      !! vector
    procedure(get_r81d), deferred      :: state_vector
      !! Returns the basal surface's state vector, a 1D array with all
      !! necessary data to describe its state.
    procedure(read_dat), deferred      :: read_data
      !! Read the basal surface data from an HDF5 file on the disc.
    procedure(write_dat), deferred     :: write_data
      !! Writes the data describing the basal surface to the disc as
      !! an HDF5 file.
    procedure(surface_solve), deferred :: solve
      !! Solves for the state of the basal surface given a particular
      !! ice shelf geometry.
  end type basal_surface

  abstract interface
    function get_scalar(this) result(property)
      import :: basal_surface
      import :: scalar_field
      class(basal_surface), intent(in) :: this
      class(scalar_field), allocatable :: property
        !! The value of whatever property of the basal surface is being
        !! returned.
    end function get_scalar

    function get_real(this) result(property)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(in) :: this
      real(r8)                         :: property
        !! The value of whatever property of the basal surface is being 
        !! returned.
    end function get_real

    function get_r81d(this) result(state_vector)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(in)          :: this
      real(r8), dimension(:), allocatable :: state_vector
        !! The state vector of the basal surface
    end function get_r81d

    subroutine setter(this, state_vector)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(inout) :: this
      real(r8), dimension(:), intent(in)  :: state_vector
        !! A real array containing the data describing the state of the
        !! basal surface.
    end subroutine setter

    subroutine time_setter(this, time)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(inout) :: this
      real(r8), intent(in)          :: time
        !! The time at which the basal surface is in the present state.
    end subroutine time_setter

    function get_i(this) result(property)
      import :: basal_surface
      class(basal_surface), intent(in) :: this
      integer                    :: property
        !! The value of whatever property of the basal surface is being
        !! returned.
    end function get_i

    subroutine read_dat(this,file_id,group_name,error)
      import :: basal_surface
      import :: hid_t
      class(basal_surface), intent(inout) :: this
      integer(hid_t), intent(in)          :: file_id
        !! The identifier for the HDF5 file/group from which the data
        !! will be read.
      character(len=*), intent(in)        :: group_name
        !! The name of the group in the HDF5 file from which to read
        !! basal surface's data.
      integer, intent(out)                :: error
        !! Flag indicating whether routine ran without error. If no
        !! error occurs then has value 0.
    end subroutine read_dat

    subroutine write_dat(this,file_id,group_name,error)
      import :: basal_surface
      import :: hid_t
      class(basal_surface), intent(in) :: this
      integer(hid_t), intent(in)       :: file_id
        !! The identifier for the HDF5 file/group in which this data is
        !! meant to be written.
      character(len=*), intent(in)     :: group_name
        !! The name to give the group in the HDF5 file storing the
        !! basal surface's data.
      integer, intent(out)             :: error
        !! Flag indicating whether routine ran without error. If no
        !! error occurs then has value 0.
    end subroutine write_dat

    subroutine surface_solve(this, ice_thickness, ice_density, &
                           ice_temperature, time, success)
      import :: basal_surface
      import :: scalar_field
      import :: r8
      class(basal_surface), intent(inout) :: this
      class(scalar_field), intent(in)     :: ice_thickness
        !! Thickness of the ice above the basal surface
      real(r8), intent(in)                :: ice_density
        !! The density of the ice above the basal surface, assumed uniform
      real(r8), intent(in)                :: ice_temperature
        !! The temperature of the ice above the basal surface, assumed uniform
      real(r8), intent(in)                :: time
        !! The time to which the basal surface should be solved
      logical, intent(out)                :: success
        !! True if the solver is successful, false otherwise
    end subroutine surface_solve
  end interface

end module basal_surface_mod
