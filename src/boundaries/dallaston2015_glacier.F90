!
!  dallaston2015_glacier.f90
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

module dallaston2015_glacier_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for the ice shelf model used by Dallaston et al. (2015).  These
  ! are Dirichlet conditions at the lower bound of the first condition
  ! as well as, for thickness, the upper bound.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use glacier_boundary_mod, only: glacier_boundary
  implicit none
  private

  type, extends(glacier_boundary), public :: dallaston2015_glacier_boundary
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! A type with procedures for getting the boundary conditions of
    ! the ice shelf model used by Dallaston et al. (2015). These are
    ! Dirichlet conditions at the lower bound of the first condition
    ! as well as, for thickness, the upper bound.
    !
    ! @TODO Consider testing the consistency conditions at the outgoing boundary?
    !
    private
    real(r8) :: thickness = 1.0_r8
      !! The thickness of the glacier at the inflowing boundary
    real(r8) :: velocity = 1.0_r8
      !! The velocity of the glacier at the inflowing boundary
  contains
    procedure :: thickness_lower_bound => dallaston2015_lower_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the thickness
      !! field.
    procedure :: thickness_upper_bound => dallaston2015_upper_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_upper_bound`/`provide_upper_bound` argument when
      !! getting or setting the raw representation of the thickness
      !! field.
    procedure :: velocity_lower_bound => dallaston2015_lower_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the velocity
      !! field.
    procedure :: boundary_residuals => dallaston2015_residuals
      !! Returns an array consisting of the difference between the
      !! required boundary values and those which actually exist. This
      !! can then be appended to a glacier's state vector. The order
      !! in which these are listed is as follows: lower thickness
      !! boundary, upper thickness boundary, lower velocity boundary,
      !! and upper velocity boundary.
    procedure :: invert_residuals => dallaston2015_invert
      !! Returns an estimate of the values of thickness and velocity
      !! at the boundaries from the given array of residuals. The
      !! order in which the residuals are stored in the array must be
      !! the same as in that produced by the `boundary_residuals`
      !! method.
  end type dallaston2015_glacier_boundary

  interface dallaston2015_glacier_boundary
    module procedure constructor
  end interface dallaston2015_glacier_boundary

contains

  pure function constructor(thickness, velocity) result(this)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Constructs a boundary condition object for an ice shelf based on
    ! the conditions used in Dallaston et al. (2015).
    !
    real(r8), intent(in) :: thickness
      !! The ice thickness at the inflowing ice shelf boundary
    real(r8), intent(in) :: velocity
      !! The longitudinal ice velocity at the inflowing ice shelf boundary
    type(dallaston2015_glacier_boundary) :: this
    this%thickness = thickness
    this%velocity = velocity
  end function constructor

  pure function dallaston2015_lower_bound(this) result(bound_array)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Indicates that one layer of cells at the lower boundary in the
    ! first dimension should be omitted. This is appropriate for
    ! Dirichlet boundary conditions.
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    integer, dimension(2) :: bound_array
    bound_array = [1,0]
  end function dallaston2015_lower_bound

  pure function dallaston2015_upper_bound(this) result(bound_array)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Indicates that one layer of cells at the upper boundary in the
    ! first dimension should be omitted for thickness. This is
    ! appropriate for Dirichlet boundary conditions.
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    integer, dimension(2) :: bound_array
    bound_array = [1,0]
  end function dallaston2015_upper_bound

  function dallaston2015_residuals(this, thickness, velocity, t) &
                                   result(residuals)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the difference between the glacier conditions of the
    ! plume and the Dirichlet conditions prescribed in the model of
    ! Dallaston et al. (2015)
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
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
    class(scalar_field), allocatable :: thickness_bound_lower, &
                                        thickness_bound_upper
    class(vector_field), allocatable :: velocity_bound
    allocate(thickness_bound_lower, &
             source=(thickness%get_boundary(-1,1) - this%thickness))
    allocate(thickness_bound_upper, &
             source=thickness%get_boundary(1,1))
    allocate(velocity_bound, &
             source=(velocity%get_boundary(-1,1) - [this%velocity]))
    residuals = [thickness_bound_lower%raw(), thickness_bound_upper%raw(), &
                 velocity_bound%raw()]
  end function dallaston2015_residuals

  function dallaston2015_invert(this, residuals, thickness, velocity, t) &
                                 result(inversion)
    class(dallaston2015_glacier_boundary), intent(in) :: this
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
    integer :: i, j
    allocate(inversion(size(residuals)))
    i = 1
    j = thickness%raw_size() - thickness%raw_size(this%thickness_lower_bound())
    inversion(i:j) = residuals(i:j) + this%thickness
    i = j + 1
    j = i + thickness%raw_size() - thickness%raw_size(this%thickness_upper_bound()) - 1
    inversion(i:j) = residuals(i:j)
    i = j + 1
    j = i + velocity%raw_size() - velocity%raw_size(this%velocity_lower_bound()) - 1
    inversion(i:j) = residuals(i:j) + this%velocity
  end function dallaston2015_invert

end module dallaston2015_glacier_boundary_mod
