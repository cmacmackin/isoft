!
!  simple_plume.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3 of the License, or
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
  !  Date: March 2017
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for a 1-D plume model. Dirichlet boundary conditions are used at
  ! the grounding line, while an outflow condition is used at the end
  ! of the domain.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field, &
                         uniform_vector_field
  use plume_boundary_mod, only: plume_boundary
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  implicit none
  private

  type(uniform_scalar_field) :: dummy
  
  type, extends(plume_boundary), public :: simple_plume_boundary
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! A type with procedures for getting the boundary conditions of
    ! the plume model. Dirichlet boundary conditions are used at the
    ! grounding line. In order to approximate an outflow condition,
    ! the derivatives of velocity, temperature, and salinity are set
    ! to 0 at the end of the domain. Plume thickness is left free
    ! there, as only a single boundary condition is needed for it.
    !
    private
    real(r8) :: thickness = 0.0_r8
      !! The thickness of the plume at the inflowing boundary
    real(r8), dimension(2) :: velocity = 1.0_r8
      !! The velocity of the plume at the inflowing boundary
    real(r8) :: salinity = 0.0_r8
      !! The salinity of the plume at the inflowing boundary
    real(r8) :: temperature = 0.0_r8
      !! The tempreature of the plume at the inflowing boundary
  contains
    procedure :: thickness_bound_info => simple_thickness_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: velocity_bound_info => simple_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: temperature_bound_info => simple_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: salinity_bound_info => simple_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: thickness_bound => simple_thickness_bound
      !! Produces a field containing the boundary conditions for plume
      !! thickness at the specified location.
    procedure :: velocity_bound => simple_velocity_bound
      !! Produces a field containing the boundary conditions for plume
      !! velocity at the specified location.
    procedure :: temperature_bound => simple_temperature_bound
      !! Produces a field containing the boundary conditions for plume
      !! temperature at the specified location.
    procedure :: salinity_bound => simple_salinity_bound
      !! Produces a field containing the boundary conditions for plume
      !! salinity at the specified location.
  end type simple_plume_boundary

  interface simple_plume_boundary
    module procedure constructor
  end interface simple_plume_boundary

contains

  pure function constructor(thickness, velocity, temperature, salinity) &
                                                           result(this)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Constructs a boundary condition object for an ice shelf based on
    ! the conditions used in Dallaston et al. (2015).
    !
    real(r8), intent(in)               :: thickness
      !! The water thickness at the inflowing plume boundary
    real(r8), dimension(2), intent(in) :: velocity
      !! The longitudinal water velocity at the inflowing plume boundary
    real(r8), intent(in)               :: temperature
      !! The water temperature at the inflowing plume boundary
    real(r8), intent(in)               :: salinity
      !! The water salinity at the inflowing plume boundary
    type(simple_plume_boundary)        :: this
    this%thickness = thickness
    this%velocity = velocity
    this%temperature = temperature
    this%salinity = salinity
  end function constructor

  subroutine simple_thickness_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is free.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)               :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    integer, intent(out)              :: bound_type
      !! An integer representing what sort of boundary condition is
      !! used. The integer value corresponding to each boundary type is
      !! specified in the [[boundary_types_mod]].
    integer, intent(out)              :: bound_depth
      !! The number of layers of data-points needed to specify the
      !! boundary condition.
    select case(location)
    case(-1)
      bound_type = dirichlet
      bound_depth = 1
    case default
      bound_type = free_boundary
      bound_depth = 0
    end select
  end subroutine simple_thickness_info

  subroutine simple_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is Neumann.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)                      :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    integer, intent(out)                     :: bound_type
      !! An integer representing what sort of boundary condition is
      !! used. The integer value corresponding to each boundary type is
      !! specified in the [[boundary_types_mod]].
    integer, intent(out)                     :: bound_depth
      !! The number of layers of data-points needed to specify the
      !! boundary condition.
    select case(location)
    case(-1)
      bound_type = dirichlet
      bound_depth = 1
    case(1)
      bound_type = neumann
      bound_depth = 1
    case default
      bound_type = free_boundary
      bound_depth = 0
    end select
  end subroutine simple_info

  function simple_thickness_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the thickness boundary values for the
    ! specified boundary location.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)                      :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer            :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%thickness)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function simple_thickness_bound

  function simple_velocity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the velocity boundary values for
    ! the specified boundary location.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)                      :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(vector_field), pointer             :: bound
    call dummy%allocate_vector_field(bound)
    select case(location)
    case(-1)
      bound = uniform_vector_field(this%velocity)
    case(1)
      bound = uniform_vector_field([0._r8, 0._r8])
    case default
      bound = uniform_vector_field([0._r8, 0._r8])
    end select
  end function simple_velocity_bound

  function simple_temperature_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the temperature boundary values for
    ! the specified boundary location.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)                      :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer             :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%temperature)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function simple_temperature_bound

  function simple_salinity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the salinity boundary values for
    ! the specified boundary location.
    !
    class(simple_plume_boundary), intent(in) :: this
    integer, intent(in)                      :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer             :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%salinity)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function simple_salinity_bound

end module simple_plume_boundary_mod
