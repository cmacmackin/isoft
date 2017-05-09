!
!  dallaston2015_seasonal.f90
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

module dallaston2015_seasonal_mod
  !* Author: Christopher MacMackin
  !  Date: May 2017
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for a 1-D plume model, when subglacial discharge is oscillating
  ! over time. This corresponds to Dirichlet conditions at the
  ! grounding line and Neumann conditions (wth a gradient of 0) at the
  ! calving front.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field, &
                         uniform_vector_field
  use plume_boundary_mod, only: plume_boundary
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  implicit none
  private

  type(uniform_scalar_field) :: dummy
  
  type, extends(plume_boundary), public :: dallaston2015_seasonal_boundary
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! A type with procedures for getting the boundary conditions of
    ! the plume model. It represents the case where subglacial
    ! discharge is varying in time, altering the boundary conditions
    ! for velocity and salinity using scalings similar to those in
    ! Dallaston et al. (2015). Dirichlet boundary conditions are used
    ! at the grounding line. In order to approximate an outflow
    ! condition, the derivatives of velocity, temperature, and
    ! salinity are set to 0 at the end of the domain. Plume thickness
    ! is left free there, as only a single boundary condition is
    ! needed for it.
    !
    private
    real(r8) :: thickness = 0.1_r8
      !! The thickness of the plume at the inflowing boundary
    real(r8) :: frequency = 1.0_r8
      !! The angular frequency of the oscillations in discharge
    real(r8) :: amplitude = 1.0_r8
      !! The amplitude of the oscillations in discharge
    real(r8) :: mean = 1.0_r8
      !! The time-average of the discharge, about which it oscillates
    real(r8) :: discharge = 1.0_r8
      !! The current discharge value
    real(r8) :: temperature = 0.0_r8
      !! The tempreature of the plume at the inflowing boundary
  contains
    procedure :: thickness_bound_info => seasonal_thickness_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: velocity_bound_info => seasonal_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: temperature_bound_info => seasonal_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: salinity_bound_info => seasonal_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: thickness_bound => seasonal_thickness_bound
      !! Produces a field containing the boundary conditions for plume
      !! thickness at the specified location.
    procedure :: velocity_bound => seasonal_velocity_bound
      !! Produces a field containing the boundary conditions for plume
      !! velocity at the specified location.
    procedure :: temperature_bound => seasonal_temperature_bound
      !! Produces a field containing the boundary conditions for plume
      !! temperature at the specified location.
    procedure :: salinity_bound => seasonal_salinity_bound
      !! Produces a field containing the boundary conditions for plume
      !! salinity at the specified location.
    procedure :: set_time => seasonal_set_time
      !! Specifies the time at which to calculate the boundary
      !! conditions.
  end type dallaston2015_seasonal_boundary

  interface dallaston2015_seasonal_boundary
    module procedure constructor
  end interface dallaston2015_seasonal_boundary

contains

  pure function constructor(thickness, frequency, amplitude, mean, &
                            temperature) result(this)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Constructs a boundary condition object for an ice shelf based on
    ! the conditions used in Dallaston et al. (2015), but with
    ! seasonal variations in subglacial discharge.
    !
    real(r8), intent(in), optional        :: thickness
      !! The plume thickness at the inflowing plume boundary, defaults
      !! to 0.1
    real(r8), intent(in), optional        :: frequency
      !! The angular frequency of the oscillations in discharge,
      !! defaults to 1.0
    real(r8), intent(in), optional        :: amplitude
      !! The amplitude of the oscillations in discharge, defaults to
      !! 1.0
    real(r8), intent(in), optional        :: mean
      !! The time-average of the discharge, about which it oscillates,
      !! defaulting to 1.0
    real(r8), intent(in), optional        :: temperature
      !! The water temperature at the inflowing plume boundary,
      !! defaults to 0.0
    type(dallaston2015_seasonal_boundary) :: this
    
    if (present(thickness)) this%thickness = thickness
    if (present(frequency)) this%frequency = frequency
    if (present(amplitude)) this%amplitude = amplitude
    if (present(mean)) this%mean = mean
    if (present(temperature)) this%temperature = temperature
  end function constructor

  subroutine seasonal_thickness_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is free.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    integer, intent(out)                               :: bound_type
      !! An integer representing what sort of boundary condition is
      !! used. The integer value corresponding to each boundary type is
      !! specified in the [[boundary_types_mod]].
    integer, intent(out)                               :: bound_depth
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
  end subroutine seasonal_thickness_info

  subroutine seasonal_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is Neumann.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    integer, intent(out)                               :: bound_type
      !! An integer representing what sort of boundary condition is
      !! used. The integer value corresponding to each boundary type is
      !! specified in the [[boundary_types_mod]].
    integer, intent(out)                               :: bound_depth
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
  end subroutine seasonal_info

  function seasonal_thickness_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the thickness boundary values for the
    ! specified boundary location.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer                       :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%thickness)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function seasonal_thickness_bound

  function seasonal_velocity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the velocity boundary values for
    ! the specified boundary location.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(vector_field), pointer                       :: bound
    call dummy%allocate_vector_field(bound)
    select case(location)
    case(-1)
      bound = uniform_vector_field([this%discharge**(1._r8/3._r8),0._r8])
    case(1)
      bound = uniform_vector_field([0._r8, 0._r8])
    case default
      bound = uniform_vector_field([0._r8, 0._r8])
    end select    
  end function seasonal_velocity_bound

  function seasonal_temperature_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the temperature boundary values for
    ! the specified boundary location.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer                       :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%temperature)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function seasonal_temperature_bound

  function seasonal_salinity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the salinity boundary values for
    ! the specified boundary location.
    !
    class(dallaston2015_seasonal_boundary), intent(in) :: this
    integer, intent(in)                                :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer                       :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%discharge**(2._r8/3._r8)/this%thickness)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
  end function seasonal_salinity_bound

  subroutine seasonal_set_time(this, time)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Sets the time at which boundary conditions are to be calculated.
    !
    class(dallaston2015_seasonal_boundary), intent(inout) :: this
    real(r8), intent(in)                                  :: time
    this%discharge = this%mean + this%amplitude*sin(this%frequency*time)
  end subroutine seasonal_set_time

end module dallaston2015_seasonal_mod
