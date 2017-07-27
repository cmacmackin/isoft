!
!  upstream_plume.F90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
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

module upstream_plume_mod
  !* Author: Christopher MacMackin
  !  Date: July 2017
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for a 1-D plume model. In order to avoid boundary layer effects,
  ! an ODE solver is used to integrate the plume a little way
  ! upstream, past the boundary layer.
  !
  use iso_fortran_env, only: r8 => real64
  use logger_mod, only: logger => master_logger
  use penf, only: str
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field, &
                         uniform_vector_field
  use uniform_gradient_field_mod, only: uniform_gradient_field
  use plume_boundary_mod, only: plume_boundary
  use boundary_types_mod, only: free_boundary, dirichlet, neumann
  use rksuite_90
  implicit none
  private

  type(uniform_scalar_field) :: dummy
  
  type, extends(plume_boundary), public :: upstream_plume_boundary
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! A type with procedures for getting the boundary conditions of
    ! the plume model. For given boundary conditions, it integrates
    ! the plume upstream slightly and then returns these when asked
    ! for boundary values. This allows boundary layers, which can
    ! cause numerical difficulties, to be avoided. Inflow boundaries
    ! must be Dirichlet, while outflow boundaries for velocity,
    ! salinity, and temperature are set to have a gradient of zero.
    !
    ! @Warning The `calculate` method __must__ be called prior to use.
    !
    private
    procedure(bound_vals), nopass, pointer :: get_boundaries => null()
      !! Calculate the "actual" boundary values, used to initiate the
      !! integration, for the specified time.
    real(r8) :: distance = 0.05_r8
      !! The distance upstream which the plume should be integrated
    real(r8) :: thickness = 0.1_r8
      !! The thickness of the plume at the inflowing boundary
    real(r8), dimension(:), allocatable :: velocity
      !! The velocity of the plume at the inflowing boundary
    real(r8) :: temperature = 0.0_r8
      !! The tempreature of the plume at the inflowing boundary
    real(r8) :: salinity = 1.0_r8
      !! The salinity of the plume at the inflowing boundary
    real(r8) :: boundary_time
      !! The time at which the boundaries were most recently
      !! calculated
    real(r8), dimension(:), allocatable :: thresholds
      !! Thresholds to use when calculating error in the integration.
  contains
    procedure :: thickness_bound_info => upstream_thickness_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: velocity_bound_info => upstream_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: temperature_bound_info => upstream_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: salinity_bound_info => upstream_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: thickness_bound => upstream_thickness_bound
      !! Produces a field containing the boundary conditions for plume
      !! thickness at the specified location.
    procedure :: velocity_bound => upstream_velocity_bound
      !! Produces a field containing the boundary conditions for plume
      !! velocity at the specified location.
    procedure :: temperature_bound => upstream_temperature_bound
      !! Produces a field containing the boundary conditions for plume
      !! temperature at the specified location.
    procedure :: salinity_bound => upstream_salinity_bound
      !! Produces a field containing the boundary conditions for plume
      !! salinity at the specified location.
    procedure :: calculate => upstream_calculate
      !! Calculates the upstreamed boundary conditions for the given
      !! time and ice thickness.
  end type upstream_plume_boundary

  abstract interface
    pure subroutine bound_vals(time, D, U, T, S)
      import :: r8
      real(r8), intent(in)                :: time
        !! The time at which the boundary values are being calculated
      real(r8), intent(out)               :: D
        !! Plume thickness boundary condition
      real(r8), dimension(:), allocatable, intent(out) :: U
        !! Plume velocity boundary condition
      real(r8), intent(out)               :: T
        !! Plume temperature boundary condition
      real(r8), intent(out)               :: S
        !! Plume salinity boundary condition
    end subroutine bound_vals

    subroutine non_diff(D, U, T, S, b, DU_x, DUU_x, DUT_x, DUS_x)
      import :: scalar_field
      import :: vector_field
      class(scalar_field), intent(in)  :: D
        !! The plume thickness
      class(vector_field), intent(in)  :: U
        !! The plume velocity
      class(scalar_field), intent(in)  :: T
        !! The plume temperature
      class(scalar_field), intent(in)  :: S
        !! The plume salinity
      class(scalar_field), intent(in)  :: b
        !! The debth of the base of the ice shelf
      class(scalar_field), intent(out) :: DU_x
        !! The derivative of the product DU
      class(vector_field), intent(out) :: DUU_x
        !! The derivative of the product DUU
      class(scalar_field), intent(out) :: DUT_x
        !! The derivative of the product DUT
      class(scalar_field), intent(out) :: DUS_x
        !! The derivative of the product DUS
    end subroutine non_diff
  end interface

  interface upstream_plume_boundary
    module procedure constructor
  end interface upstream_plume_boundary

contains

  pure function constructor(bound_calculator, distance, thresholds) &
                                                       result(this)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Constructs a boundary condition object which integrates an IVP
    ! from actual boundary values to calculate the staet of the plume
    ! a little upstream. This can be used to avoid boundary layers.
    !
    procedure(bound_vals)                        :: bound_calculator
      !! Calculates the "actual" inflow boundary conditions, used to
      !! initiate the integration to find the values to use in the
      !! simulation.
    real(r8), intent(in)                         :: distance
      !! The distance upstream which the plume should be integrated.
    real(r8), dimension(:), optional, intent(in) :: thresholds
      !! The thresholds to use when evaluating the error of the
      !! integration. This is done according to the formula
      !! `abs(e) / max(magnitude_y, THRESHOLDS) <= TOLERANCE`.
    type(upstream_plume_boundary) :: this
    this%get_boundaries => bound_calculator
    this%distance = distance
    if (present(thresholds)) this%thresholds = thresholds
  end function constructor

  subroutine upstream_thickness_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is free.
    !
    class(upstream_plume_boundary), intent(in) :: this
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
  end subroutine upstream_thickness_info

  subroutine upstream_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Indicates that the lower boundary is Dirichlet and the upper
    ! boundary is Neumann.
    !
    class(upstream_plume_boundary), intent(in) :: this
    integer, intent(in)                        :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    integer, intent(out)                       :: bound_type
      !! An integer representing what sort of boundary condition is
      !! used. The integer value corresponding to each boundary type is
      !! specified in the [[boundary_types_mod]].
    integer, intent(out)                       :: bound_depth
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
  end subroutine upstream_info

  function upstream_thickness_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Returns a field containing the thickness boundary values for the
    ! specified boundary location.
    !
    class(upstream_plume_boundary), intent(in) :: this
    integer, intent(in)                        :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer               :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%thickness)
    case default
      bound = uniform_scalar_field(0._r8)
    end select    
    call bound%set_temp() ! Shouldn't need to call this, but for some
                          ! reason being set as non-temporary when
                          ! assignment subroutine returns.
  end function upstream_thickness_bound

  function upstream_velocity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Returns a field containing the velocity boundary values for
    ! the specified boundary location.
    !
    class(upstream_plume_boundary), intent(in) :: this
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
      bound = uniform_vector_field(this%velocity)
    case(1)
      bound = uniform_vector_field([0._r8, 0._r8])
    case default
      bound = uniform_vector_field([0._r8, 0._r8])
    end select
    call bound%set_temp() ! Shouldn't need to call this, but for some
                          ! reason being set as non-temporary when
                          ! assignment subroutine returns.
  end function upstream_velocity_bound

  function upstream_temperature_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Returns a field containing the temperature boundary values for
    ! the specified boundary location.
    !
    class(upstream_plume_boundary), intent(in) :: this
    integer, intent(in)                        :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer               :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%temperature)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select
    call bound%set_temp() ! Shouldn't need to call this, but for some
                          ! reason being set as non-temporary when
                          ! assignment subroutine returns.
  end function upstream_temperature_bound

  function upstream_salinity_bound(this, location) result(bound)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Returns a field containing the salinity boundary values for
    ! the specified boundary location.
    !
    class(upstream_plume_boundary), intent(in) :: this
    integer, intent(in)                        :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), pointer               :: bound
    call dummy%allocate_scalar_field(bound)
    select case(location)
    case(-1)
      bound = uniform_scalar_field(this%salinity)
    case(1)
      bound = uniform_scalar_field(0._r8)
    case default
      bound = uniform_scalar_field(0._r8)
    end select
    call bound%set_temp() ! Shouldn't need to call this, but for some
                          ! reason being set as non-temporary when
                          ! assignment subroutine returns.
  end function upstream_salinity_bound

  subroutine upstream_calculate(this, t, func, b)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Calculates the boundary values to use at the current time with
    ! the current ice thickness.
    !
    class(upstream_plume_boundary), intent(inout) :: this
    real(r8), intent(in)                          :: t
      !! The time at which to calculate the boundary values.
    procedure(non_diff)                           :: func
      !! A function which returns the non-diffusive, non-inertial
      !! components of the ODEs describing the plume.
    class(scalar_field), intent(in)               :: b
      !! The depth of the ice shelf base.

    integer :: n
    real(r8) :: D_0, U_0, S_0, T_0, DU_0, DUS_0, DUT_0
    real(r8), dimension(:), allocatable :: Uvec_0, DUvecU_0
    class(scalar_field), pointer :: b_x
    real(r8), dimension(:), allocatable :: y0, ygot, yderiv_got
    real(r8) :: xgot
    integer :: flag
    type(rk_comm_real_1d) :: comm

    call b%guard_temp()
    call this%get_boundaries(t, D_0, Uvec_0, T_0, S_0)
    U_0 = Uvec_0(1)
    DU_0 = D_0*U_0
    DUvecU_0 = DU_0*Uvec_0
    DUS_0 = DU_0*S_0
    DUT_0 = DU_0*T_0
    y0 = [DU_0, DUvecU_0, DUS_0, DUT_0]
    n = size(y0)
    allocate(ygot(n), yderiv_got(n))
    if (.not. allocated(this%thresholds)) then
      allocate(this%thresholds(n))
      this%thresholds = 1._r8
    end if
    b_x => b%d_dx(1)
    call b_x%guard_temp()

    call setup(comm, 0._r8, y0, this%distance, 1e-6_r8, this%thresholds, &
               method='h', h_start=2e-2_r8*this%distance)
    call range_integrate(comm, integrand, this%distance, xgot, ygot, &
                         yderiv_got, flag)
    if (flag /= 1) then
      call logger%error('upstream_plume_boundary%calculate',            &
                        'Could only integrate plume to x = '//str(xgot))
    end if
    this%thickness = ygot(1)**2/ygot(2)
    this%velocity = ygot(2:n-2)/ygot(1)
    this%temperature = ygot(n-1)/ygot(1)
    this%salinity = ygot(n)/ygot(1)
    !print*,'t = ', t
    !print*,this%thickness, this%velocity, this%temperature, this%salinity
    if (this%thickness < 0._r8) error stop
    call collect_garbage(comm)
    call b%clean_temp(); call b_x%clean_temp()

  contains

    function integrand(x, y) result(f)
      real(r8), intent(in)               :: x
        !! The independent variable
      real(r8), dimension(:), intent(in) :: y
        !! The dependent variable
      real(r8), dimension(size(y))       :: f
        !! The derivatives

      type(uniform_scalar_field)   :: D, T, S, DU_x, DUT_x, DUS_x
      type(uniform_vector_field)   :: Uvec, DUU_x
      type(uniform_gradient_field) :: b_loc
      D = uniform_scalar_field(y(1)**2/y(2))
      Uvec = uniform_vector_field(y(2:n-2)/y(1))
      T = uniform_scalar_field(y(n-1)/y(1))
      S = uniform_scalar_field(y(n)/y(1))
      b_loc = uniform_gradient_field(b%interpolate([x]), [b_x%interpolate([x])])
      call func(D, Uvec, T, S, b_loc, DU_x, DUU_x, DUT_x, DUS_x)
      f(1) = DU_x%get_value()
      f(2:n-2) = DUU_x%get_value()
      f(n-1) = DUT_x%get_value()
      f(n) = DUS_x%get_value()
    end function integrand
    
  end subroutine upstream_calculate

end module upstream_plume_mod
