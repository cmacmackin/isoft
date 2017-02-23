!
!  dallaston2015_plume.f90
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

module dallaston2015_plume_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for the plume model used by Dallaston et al. (2015). This
  ! corresponds to Dirichlet conditions at the lower bound of the
  ! first dimension.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use plume_boundary_mod, only: plume_boundary
  implicit none
  private
  
  type, extends(plume_boundary), public :: dallaston2015_plume_boundary
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! A type with procedures for getting the boundary conditions of
    ! the plume model used by Dallaston et al. (2015). These are
    ! Dirichlet conditions at lower bound of the first dimension and
    ! free conditions at the upper bound. Temperature has a free lower
    ! boundary condition as well, due to the singular nature of that
    ! boundary.
    !
    private
    real(r8) :: thickness = 0.0_r8
      !! The thickness of the plume at the inflowing boundary
    real(r8) :: velocity = 1.0_r8
      !! The velocity of the plume at the inflowing boundary
    real(r8) :: salinity = 0.0_r8
      !! The salinity of the plume at the inflowing boundary
  contains
    procedure :: thickness_lower_bound => dallaston2015_lower_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the thickness
      !! field.
    procedure :: velocity_lower_bound => dallaston2015_lower_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the velocity
      !! field.
    procedure :: salinity_lower_bound => dallaston2015_lower_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_lower_bound`/`provide_lower_bound` argument when
      !! getting or setting the raw representation of the salinity
      !! field.
    procedure :: boundary_residuals => dallaston2015_residuals
      !! Returns an array consisting of the difference between the required
      !! boundary values and those which actually exist. This can then be
      !! appended to a plume's state vector
    procedure :: set_boundaries => dallaston2015_set
  end type dallaston2015_plume_boundary

  interface dallaston2015_plume_boundary
    module procedure constructor
  end interface dallaston2015_plume_boundary

contains

  pure function constructor(thickness, velocity, salinity) result(this)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Constructs a boundary condition object for an ice shelf based on
    ! the conditions used in Dallaston et al. (2015).
    !
    real(r8), intent(in) :: thickness
      !! The water thickness at the inflowing plume boundary
    real(r8), intent(in) :: velocity
      !! The longitudinal water velocity at the inflowing plume boundary
    real(r8), intent(in) :: salinity
      !! The water salinity at the inflowing plume boundary
    type(dallaston2015_plume_boundary) :: this
    this%thickness = thickness
    this%velocity = velocity
    this%salinity = salinity
  end function constructor

  pure function dallaston2015_lower_bound(this) result(bound_array)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Indicates that one layer of cells at the lower boundary in the
    ! first dimension should be omitted. This is appropriate for
    ! Dirichlet boundary conditions.
    !
    class(dallaston2015_plume_boundary), intent(in) :: this
    integer, dimension(2) :: bound_array
    bound_array = [1,0]
  end function dallaston2015_lower_bound

  function dallaston2015_residuals(this, thickness, velocity, temperature, &
                                   salinity, t) result(residuals)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the difference between the boundary conditions of the
    ! plume and the Dirichlet conditions prescribed in the model of
    ! Dallaston et al. (2015)
    !
    class(dallaston2015_plume_boundary), intent(in)   :: this
    class(scalar_field), intent(in)     :: thickness
      !! A field containing the thickness of the plume
    class(vector_field), intent(in)     :: velocity
      !! A field containing the flow velocity of the plume
    class(scalar_field), intent(in)     :: temperature
      !! The field containing the temperature of the plume
    class(scalar_field), intent(in)     :: salinity
      !! The field containing the salinity of the plume
    real(r8), intent(in)                :: t
      !! The time at which the boundary conditions are to be calculated.
    real(r8), allocatable, dimension(:) :: residuals
      !! An array containing the difference between the required boundary
      !! values and those which are actually present.
    class(scalar_field), allocatable :: thickness_bound, salinity_bound
    class(vector_field), allocatable :: velocity_bound
    allocate(thickness_bound,   &
             source=(thickness%get_boundary(-1,1) - this%thickness))
    allocate(velocity_bound,    &
             source=(velocity%get_boundary(-1,1) - [this%velocity]))
    allocate(salinity_bound,    &
             source=(salinity%get_boundary(-1,1) - this%salinity))
    residuals = [thickness_bound%raw(), velocity_bound%raw(), &
                 salinity_bound%raw()]
  end function dallaston2015_residuals

  subroutine dallaston2015_set(this, thickness, velocity, temperature, salinity, t)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    !
    class(dallaston2015_plume_boundary), intent(inout) :: this
    class(scalar_field), intent(inout)                 :: thickness
      !! A field containing the thickness of the plume
    class(vector_field), intent(inout)                 :: velocity
      !! A field containing the flow velocity of the plume
    class(scalar_field), intent(inout)                 :: temperature
      !! The field containing the temperature of the plume
    class(scalar_field), intent(inout)                 :: salinity
      !! The field containing the salinity of the plume
    real(r8), intent(in)                               :: t
      !! The time at which the boundary conditions are to be updated.
    integer :: n
    n  = thickness%elements()
    call thickness%set_element(n,this%thickness)
    call velocity%set_element(n,1,this%velocity*this%thickness)
    call salinity%set_element(n,this%salinity*this%thickness)
  end subroutine dallaston2015_set

end module dallaston2015_plume_boundary_mod
