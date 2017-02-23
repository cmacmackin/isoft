!
!  simple_plume.F90
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

module simple_plume_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: February 2017
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for a typical plume model. Conditions are Dirichlet at the lower
  ! boundary and, for variables other than thickness, Neumann at the
  ! upper boundary.
  !
  use iso_fortran_env, only: r8 => real64
  use chebyshev_mod, only: differentiation_row
  use factual_mod, only: scalar_field, vector_field
  use plume_boundary_mod, only: plume_boundary
  implicit none
  private
  
  type, extends(plume_boundary), public :: simple_plume_boundary
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! Provides a derived type which specifies the boundary conditions
    ! for a typical plume model. Conditions are Dirichlet at the lower
    ! boundary and, for variables other than thickness, Neumann at the
    ! upper boundary.
    !
    private
    real(r8) :: thickness = 0.0_r8
      !! The thickness of the plume at the inflowing boundary
    real(r8) :: velocity = 1.0_r8
      !! The velocity of the plume at the inflowing boundary
    real(r8) :: salinity = 0.0_r8
      !! The salinity of the plume at the inflowing boundary
    real(r8) :: temperature = 0.0_r8
      !! The temperature of the plume at the inflowing boundary
    real(r8) :: d_velocity = 0.0_r8
      !! The first derivative of the plume velocity at the outflowing
      !! boundary
    real(r8) :: d_salinity = 0.0_r8
      !! The first derivative of the plume salinity at the outflowing
      !! boundary
    real(r8) :: d_temperature = 0.0_r8
      !! The first derivative of the plume temperature at the
      !! outflowing boundary
    real(r8), dimension(:), allocatable :: diff_operator
      !! A row of a differentiation matrix.
  contains
    procedure :: thickness_lower_bound => simple_bound
    procedure :: velocity_lower_bound => simple_bound
    procedure :: salinity_lower_bound => simple_bound
    procedure :: temperature_lower_bound => simple_bound
    procedure :: velocity_upper_bound => simple_bound
    procedure :: salinity_upper_bound => simple_bound
    procedure :: temperature_upper_bound => simple_bound
    procedure :: boundary_residuals => simple_residuals
    procedure :: set_boundaries => simple_set
  end type simple_plume_boundary

  interface simple_plume_boundary
    module procedure constructor
  end interface simple_plume_boundary

contains

  pure function constructor(thickness, velocity, salinity, temperature, &
                            d_velocity, d_salinity, d_temperature) result(this)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! Constructs a boundary condition object for an typical plume.
    !
    real(r8), intent(in), optional :: thickness
      !! The water thickness at the inflowing plume boundary
    real(r8), intent(in), optional :: velocity
      !! The longitudinal water velocity at the inflowing plume boundary
    real(r8), intent(in), optional :: salinity
      !! The water salinity at the inflowing plume boundary
    real(r8), intent(in), optional :: temperature
      !! The water temperature at the inflowing plume boundary
    real(r8), intent(in), optional :: d_velocity
      !! The first derivative of the water velocity at the outflowing
      !! boundary
    real(r8), intent(in), optional :: d_salinity
      !! The first derivative of the water salinity at the outflowing
      !! boundary
    real(r8), intent(in), optional :: d_temperature
      !! The first derivative of the water temperature at the
      !! outflowing boundary
    type(simple_plume_boundary) :: this
    if (present(thickness)) this%thickness = thickness
    if (present(velocity)) this%velocity = velocity
    if (present(salinity)) this%salinity = salinity
    if (present(temperature)) this%temperature = temperature
    if (present(d_velocity)) this%d_velocity = d_velocity
    if (present(d_salinity)) this%d_salinity = d_salinity
    if (present(d_temperature)) this%d_temperature = d_temperature
  end function constructor

  pure function simple_bound(this) result(bound_array)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    ! Indicates that one layer of cells at the boundary in the
    ! first dimension should be omitted. This is appropriate for
    ! Dirichlet or Neumann boundary conditions.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, dimension(2) :: bound_array
    bound_array = [1,0]
  end function simple_bound

  function simple_residuals(this, thickness, velocity, temperature, &
                                   salinity, t) result(residuals)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    !
    class(simple_plume_boundary), intent(in)   :: this
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
    error stop ('Not implemented')
  end function simple_residuals

  subroutine simple_set(this, thickness, velocity, temperature, salinity, t)
    !* Author: Chris MacMackin
    !  Date: February 2017
    !
    !
    class(simple_plume_boundary), intent(inout) :: this
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
    real(r8) :: e_thickness, e_val, d_e_thickness
    real(r8), dimension(2,1) :: extent
    real(r8), dimension(:), allocatable :: raw_vals
    integer  :: n
    n  = thickness%elements()
    call thickness%set_element(n,this%thickness)
    call velocity%set_element(n,1,this%velocity*this%thickness)
    call salinity%set_element(n,this%salinity*this%thickness)
    call temperature%set_element(n,this%temperature*this%thickness)
    e_thickness = thickness%get_element(1)
    if (.not. allocated(this%diff_operator)) then
      extent = thickness%domain()
      this%diff_operator = differentiation_row(n-1,1,extent(1,1),extent(2,1))
    end if
    d_e_thickness = dot_product(this%diff_operator, thickness%raw())
    raw_vals = velocity%raw()
    e_val = (e_thickness*this%d_velocity + d_e_thickness*velocity%get_element(1,1) - &
             dot_product(this%diff_operator(2:),raw_vals(2:n)))/this%diff_operator(1)
    call velocity%set_element(1,1,e_val)
    raw_vals = salinity%raw()
    e_val = (e_thickness*this%d_salinity + d_e_thickness*salinity%get_element(1) - &
             dot_product(this%diff_operator(2:),raw_vals(2:)))/this%diff_operator(1)
    call salinity%set_element(1,e_val)
    raw_vals = temperature%raw()
    e_val = (e_thickness*this%d_temperature + d_e_thickness*salinity%get_element(1) - &
             dot_product(this%diff_operator(2:),raw_vals(2:)))/this%diff_operator(1)
    call temperature%set_element(1,e_val)
  end subroutine simple_set

end module simple_plume_boundary_mod
