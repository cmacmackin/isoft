!
!  glacier.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@gmail.com>
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

module glacier_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract type to represent large masses of ice, such
  ! as ice sheets and ice shelves.
  !
  use iso_fortran_env, only: r8 => real64
  use foodie, only: integrand
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private

  type, abstract, public :: glacier
    !* Author: Christopehr MacMackin
    !  Date: April 2016
    !
    ! An abstract data type which represents large masses of ice, such
    ! as ice shelves and ice sheets.
    !
  contains
    procedure(get_scalar), deferred   :: ice_thickness
      !! Returns the thickness of the ice in the glacier across the domain.
    procedure(get_r8), deferred       :: ice_density
      !! Returns the density of the ice, which is assumed to be uniform.
    procedure(get_r8), deferred       :: ice_temperature
      !! Returns the temperature of the ice, which is assumed to be uniform.
    procedure(get_residual), deferred :: residual
      !! Computes the residual of the system of equations describing the
      !! glacier.
    procedure(setter), deferred       :: update
      !! Sets the state of the glacier.
  end type glacier

  abstract interface
    function get_scalar(this) result(property)
      import :: glacier
      import :: scalar_field
      class(glacier), intent(in)       :: this
      class(scalar_field), allocatable :: property
        !! The value of whatever property of the glacier is being returned.
    end function get_scalar
    
!$    function get_vector(this) result(property)
!$      import :: glacier
!$      import :: vector_field
!$      class(glacier), intent(in)       :: this
!$      class(vector_field), allocatable :: property
!$        !! The value of whatever property of the glacier is being returned.
!$    end function get_vector
    
    function get_r8(this) result(property)
      import :: glacier
      import :: r8
      class(glacier), intent(in) :: this
      real(r8)                   :: property
        !! The value of whatever property of the glacier is being returned.
    end function get_r8

    function get_residual(this, melt_rate, basal_drag_parameter, &
                          water_density) result(residual)
      import :: glacier
      import :: scalar_field
      import :: r8
      class(glacier), intent(in)          :: this
      class(scalar_field), intent(in)     :: melt_rate
        !! Thickness of the ice above the glacier
      class(scalar_field), intent(in)     :: basal_drag_parameter
        !! A paramter, e.g. coefficient of friction, needed to calculate the
        !! drag on basal surface of the glacier.
      class(scalar_field), intent(in)     :: water_density
        !! The density of the water below the glacier
      real(r8), dimension(:), allocatable :: residual
        !! The residual of the system of equations describing the glacier
    end function get_residual

    subroutine setter(this, state_vector, time)
      import :: glacier
      import :: r8
      class(glacier), intent(inout)      :: this
      real(r8), dimension(:), intent(in) :: state_vector
        !! A real array containing the data describing the state of the
        !! glacier.
      real(r8), intent(in), optional     :: time
        !! The time at which the glacier is in this state. If not present
        !! then assumed to be same as previous value passed.
    end subroutine setter
  end interface

  abstract interface
    function thickness_func(location) result(thickness)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the [[glacier]] thickness
      ! when a concrete object is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8) :: thickness
        !! The thickness of the glacier at `location`
    end function thickness_func
      
    function velocity_func(location) result(velocity)
      !* Author: Chris MacMackin
      !  Date: April 2016
      !
      ! Abstract interface for function providing the [[glacier]] velocity
      ! when a concrete object is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8), dimension(:), allocatable :: thickness
        !! The velocity vector of the ice in the glacier at `location`
    end function velocity_func
  end interface

  public :: thickness_func, velocity_func

end module glacier_mod
