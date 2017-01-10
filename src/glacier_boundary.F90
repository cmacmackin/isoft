!
!  glacier_boundary.f90
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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module glacier_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! specify the boundary conditions for glacier types.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private
  
  type, public :: glacier_boundary
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! A type in which procedures for getting the boundary conditions
    ! of glaciers are to be specified. The descendent types can
    ! contain whatever data is needed to compute the result. It
    ! provides the routine
    ! [[glacier_boundary(type):get_boundaries_residual]] to return an
    ! array with the residuals representing deviation from satisfying
    ! the conditions. This can then be appended to a
    ! [[glacier(type)]]'s residual array.
    !
    ! When specifying boundary conditions in this way, it becomes
    ! necessary not to include the residual of the plume data at some
    ! location(s). Typically this would be adjacent to the boundaries
    ! being prescribed. In order to accomplish this, routines are
    ! provided which return arrays indicating at which boundaries data
    ! should be omitted and how much. These should be passed as
    ! arguments for the field methods to get and set the field's raw
    ! data.
    !
    ! This class effectively provides free boundary conditions. It's 
    ! type-bound procedures should be overridden to provide case-specific
    ! conditions.
    !
  contains
    procedure :: thickness_lower_bound => bound_array
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the thickness
      !! field.
    procedure :: thickness_upper_bound => bound_array
      !! Returns a 1D array which should be passed as the
      !! `exclude_upper_bound`/`provide_upper_bound` argument when
      !! getting or setting the raw representation of the thickness
      !! field.
    procedure :: velocity_lower_bound => bound_array
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the velocity
      !! field.
    procedure :: velocity_upper_bound => bound_array
      !! Returns a 1D array which should be passed as the
      !! `exclude_upper_bound`/`provide_upper_bound` argument when
      !! getting or setting the raw representation of the velocity
      !! field.
    procedure :: boundary_residuals
      !! Returns an array consisting of the difference between the
      !! required boundary values and those which actually exist. The
      !! order in which these are listed is as follows: lower
      !! thickness boundary, upper thickness boundary, lower velocity
      !! boundary, and upper velocity boundary.
    procedure :: invert_residuals
      !! Returns an estimate of the values of thickness and velocity
      !! at the boundaries from the given array of residuals. The
      !! order in which the residuals are stored in the array must be
      !! the same as in that produced by the `boundary_residuals`
      !! method.
  end type glacier_boundary

  abstract interface
  end interface

contains
 
  pure function bound_array(this)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Default implementation of the method getting lower and upper
    ! boundary information, which is then passed to the methods for
    ! getting and setting raw representations of fields. It returns a
    ! 1D array of length 2, indicating free boundaries (the raw data
    ! should represent all cells contained in the field, not excluding
    ! any near the boundaries).
    !
    class(glacier_boundary), intent(in) :: this
    integer, dimension(2) :: bound_array
    bound_array = 0
  end function bound_array

  function boundary_residuals(this, thickness, velocity, t) &
                                 result(residuals)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Default implementation of the
    ! [[glacier_boundary(type):boundary_residuals]]
    ! method. It returns a zero-length array, effectively indicating free
    ! boundaries.
    !
    class(glacier_boundary), intent(in) :: this
    class(scalar_field), intent(in)     :: thickness
      !! A field containing the thickness of the glacier
    class(vector_field), intent(in)     :: velocity
      !! A field containing the flow velocity of the glacier
    real(r8), intent(in)                :: t
      !! The time at which the boundary conditions are to be
      !! calculated.
    real(r8), allocatable, dimension(:) :: residuals
      !! An array containing the difference between the required
      !! boundary values and those which are actually present. They
      !! are stored in the order: lower thickness boundary, upper
      !! thickness boundary, lower velocity boundary, and upper
      !! velocity boundary.
    allocate(residuals(0))
    return
  end function boundary_residuals

  function invert_residuals(this, residuals, thickness, velocity, t) &
                                 result(inversion)
    class(glacier_boundary), intent(in) :: this
    real(r8), dimension(:), intent(in)  :: residuals
      !! An array containing the difference between the required
      !! boundary values and those which are actually present. The
      !! storage order must be the same as in the result of the
      !! `boundary_residuals` funciton.
    class(scalar_field), intent(in)     :: thickness
      !! A field containing the thickness of the glacier
    class(vector_field), intent(in)     :: velocity
      !! A field containing the flow velocity of the glacier
    real(r8), intent(in)                :: t
      !! The time at which the boundary conditions are to be
      !! calculated.
    real(r8), allocatable, dimension(:) :: inversion
      !! An array containing estimates of the values of the thickness
      !! at the lower boundary, thickness at the upper boundary,
      !! velocity at the lower boundary, and velocity of the upper
      !! boundary, in that order. These are calculated from the
      !! residuals.
    allocate(inversion(0))
  end function invert_residuals

end module glacier_boundary_mod
