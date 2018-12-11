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
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
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
    real(r8) :: chi = 1.0_r8
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}$
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
    procedure :: velocity_upper_bound => dallaston2015_upper_bound
      !! Returns a 1D array which should be passed as the
      !! `exclude_upper_bound`/`provide_upper_bound` argument when
      !! getting or setting the raw representation of the velocity
      !! field.
    procedure :: thickness_lower_type => dallaston2015_lower_type
      !! Returns an array indicating what type of boundary conditions
      !! apply for thickness at the lower boundary of each
      !! dimension. The types are specified using the parameters in
      !! [boundary_types_mod].
    procedure :: velocity_lower_type => dallaston2015_lower_type
      !! Returns an array indicating what type of boundary conditions
      !! apply for velocity at the lower boundary of each
      !! dimension. The types are specified using the parameters in
      !! [boundary_types_mod].
    procedure :: velocity_upper_type => dallaston2015_upper_type
      !! Returns an array indicating what type of boundary conditions
      !! apply for velocity at the upper boundary of each
      !! dimension. The types are specified using the parameters in
      !! [boundary_types_mod].
    procedure :: boundary_residuals => dallaston2015_residuals
      !! Returns an array consisting of the difference between the
      !! required boundary values and those which actually exist. This
      !! can then be appended to a glacier's state vector. The order
      !! in which these are listed is as follows: lower thickness
      !! boundary, upper thickness boundary, lower velocity boundary,
      !! and upper velocity boundary.
  end type dallaston2015_glacier_boundary

  interface dallaston2015_glacier_boundary
    module procedure constructor
  end interface dallaston2015_glacier_boundary

contains

  pure function constructor(thickness, velocity, chi) result(this)
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
    real(r8), intent(in) :: chi
      !! The dimensionless ratio
      !! $\chi \equiv \frac{\rho_igh_0x_x}{2\eta_0u_0}$
    type(dallaston2015_glacier_boundary) :: this
    this%thickness = thickness
    this%velocity = velocity
    this%chi = chi
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

  pure function dallaston2015_lower_type(this) result(bound_type)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Specifies that the lower boundary in the first dimension has
    ! Dirichlet boundary conditions.
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    integer, dimension(:), allocatable  :: bound_type
    bound_type = [dirichlet, free_boundary]
  end function dallaston2015_lower_type

  pure function dallaston2015_upper_type(this) result(bound_type)
    !* Author: Chris MacMackin
    !  Date: January 2017
    !
    ! Specifies that the upper boundary in the first dimension has
    ! Neumann boundary conditions.
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    integer, dimension(:), allocatable  :: bound_type
    bound_type = [neumann, free_boundary]
  end function dallaston2015_upper_type

  function dallaston2015_residuals(this, thickness, velocity, viscosity, t) &
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
    class(scalar_field), intent(in)     :: viscosity
      !! A field containing the viscosity of the ice in the glacier.
    real(r8), intent(in)                :: t
      !! The time at which the boundary conditions are to be
      !! calculated.
    real(r8), allocatable, dimension(:) :: residuals
      !! An array containing the difference between the required
      !! boundary values and those which are actually present. They
      !! are stored in the order: lower thickness boundary, upper
      !! thickness boundary, lower velocity boundary, and upper
      !! velocity boundary.
    class(scalar_field), pointer :: thickness_bound,      &
                                    velocity_bound_upper, &
                                    velocity_deriv
    class(vector_field), pointer :: velocity_bound_lower
    call thickness%guard_temp(); call velocity%guard_temp(); call viscosity%guard_temp()
    call thickness%allocate_scalar_field(thickness_bound)
    thickness_bound = thickness%get_boundary(-1,1) - this%thickness
    call velocity%allocate_vector_field(velocity_bound_lower)
    velocity_bound_lower = velocity%get_boundary(-1,1) - [this%velocity]
    call velocity%allocate_scalar_field(velocity_deriv)
    velocity_deriv = velocity%component_d_dx(1,1)
    call velocity%allocate_scalar_field(velocity_bound_upper)
    velocity_bound_upper = velocity_deriv%get_boundary(1,1)  &
                         * thickness%get_boundary(1,1)       &
                         - (0.25_r8*this%chi)*thickness%get_boundary(1,1)**2 &
                         / viscosity%get_boundary(1,1)
    residuals = [thickness_bound%raw(), velocity_bound_lower%raw(), &
                 velocity_bound_upper%raw()]
    call thickness%clean_temp(); call velocity%clean_temp(); call viscosity%clean_temp()
  end function dallaston2015_residuals

end module dallaston2015_glacier_boundary_mod
