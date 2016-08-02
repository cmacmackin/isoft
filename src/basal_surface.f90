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
  implicit none
  private

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
    procedure(get_scalar), deferred :: basal_melt
      !! Returns the basal melt rate.
    procedure(get_scalar), deferred :: basal_drag_parameter
      !! Returns a value which may be needed to calculate basal drag,
      !! such as the coefficient of friction.
    procedure(get_scalar), deferred :: water_density
      !! Density of the water at the basal surface.
    procedure(get_residual), deferred :: residual
      !! Computes the residual of the system of equations describing the
      !! basal surface.
    procedure(setter), deferred       :: update
      !! Sets the state of the basal surface
    procedure(get_i), deferred        :: data_size
      !! Returns the number of elements in the glacier's state vector
    procedure(get_r81d), deferred     :: state_vector
      !! Returns the glacier's state vector, a 1D array with all necessary 
      !! data to describe its state.
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
    
    function get_residual(this, ice_thickness, ice_density, ice_temperature) &
                                                              result(residual)
      import :: basal_surface
      import :: scalar_field
      import :: r8
      class(basal_surface), intent(in)    :: this
      class(scalar_field), intent(in)     :: ice_thickness
        !! Thickness of the ice above the basal surface
      real(r8), intent(in)                :: ice_density
        !! The density of the ice above the basal surface, assumed uniform
      real(r8), intent(in)                :: ice_temperature
        !! The temperature of the ice above the basal surface, assumed uniform
      real(r8), dimension(:), allocatable :: residual
        !! The residual of the system of equations describing the basal
        !! surface.
    end function get_residual

    pure function get_r81d(this) result(state_vector)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(in)          :: this
      real(r8), dimension(:), allocatable :: state_vector
        !! The state vector of the basal surface
    end function get_r81d

    subroutine setter(this, state_vector, time)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(inout) :: this
      real(r8), dimension(:), intent(in)  :: state_vector
        !! A real array containing the data describing the state of the
        !! basal surface.
      real(r8), intent(in), optional      :: time
        !! The time at which the basal surface is in this state. If not
        !! present then assumed to be same as previous value passed.
    end subroutine setter

    pure function get_i(this) result(property)
      import :: basal_surface
      class(basal_surface), intent(in) :: this
      integer                    :: property
        !! The value of whatever property of the basal surface is being
        !! returned.
    end function get_i
  end interface

end module basal_surface_mod
