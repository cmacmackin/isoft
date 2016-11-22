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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module glacier_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract type to represent large masses of ice, such
  ! as ice sheets and ice shelves.
  !
  use iso_fortran_env, only: r8 => real64
  !use foodie, only: integrand
  use factual_mod, only: scalar_field, vector_field
  use nitsol_mod, only: nitsol, dummy_jacv, ddot, dnrm2
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
    procedure(time_setter), deferred  :: set_time
      !! Sets the time record for this glacier.
    procedure(get_i), deferred        :: data_size
      !! Returns the number of elements in the glacier's state vector
    procedure(get_r81d), deferred     :: state_vector
      !! Returns the glacier's state vector, a 1D array with all necessary 
      !! data to describe its state.
    procedure                         :: integrate => glacier_integrate
      !! Performs a time-step of the integration, taking the state of
      !! the glacier to the specified time using the provided
      !! melt-rate data.
  end type glacier

  abstract interface
    pure function get_scalar(this) result(property)
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
    
    pure function get_r8(this) result(property)
      import :: glacier
      import :: r8
      class(glacier), intent(in) :: this
      real(r8)                   :: property
        !! The value of whatever property of the glacier is being returned.
    end function get_r8

    function get_residual(this, previous_states, melt_rate, &
                          basal_drag_parameter, water_density) result(residual)
      import :: glacier
      import :: scalar_field
      import :: r8
      class(glacier), intent(in)               :: this
      class(glacier), dimension(:), intent(in) :: previous_states
        !! The states of the glacier in the previous time steps. The
        !! first element of the array should be the most recent. The
        !! default implementation will only make use of the most
        !! recent state, but the fact that this is an array allows
        !! overriding methods to use older states for higher-order
        !! integration methods.
      class(scalar_field), intent(in)          :: melt_rate
        !! Thickness of the ice above the glacier
      class(scalar_field), intent(in)          :: basal_drag_parameter
        !! A paramter, e.g. coefficient of friction, needed to calculate the
        !! drag on basal surface of the glacier.
      real(r8), intent(in)                     :: water_density
        !! The density of the water below the glacier
      real(r8), dimension(:), allocatable      :: residual
        !! The residual of the system of equations describing the glacier
    end function get_residual

    pure function get_r81d(this) result(state_vector)
      import :: glacier
      import :: r8
      class(glacier), intent(in)          :: this
      real(r8), dimension(:), allocatable :: state_vector
        !! The state vector of the glacier
    end function get_r81d

    subroutine setter(this, state_vector)
      import :: glacier
      import :: r8
      class(glacier), intent(inout)      :: this
      real(r8), dimension(:), intent(in) :: state_vector
        !! A real array containing the data describing the state of the
        !! glacier.
    end subroutine setter

    subroutine time_setter(this, time)
      import :: glacier
      import :: r8
      class(glacier), intent(inout) :: this
      real(r8), intent(in)          :: time
        !! The time at which the glacier is in the present state.
    end subroutine time_setter

    pure function get_i(this) result(property)
      import :: glacier
      class(glacier), intent(in) :: this
      integer                    :: property
        !! The value of whatever property of the glacier is being returned.
    end function get_i
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
      !  Date: July 2016
      !
      ! Abstract interface for function providing the [[glacier]] velocity
      ! when a concrete object is being instantiated.
      !
      import :: r8
      real(r8), dimension(:), intent(in) :: location
        !! The position $\vec{x}$ at which to compute the thickness
      real(r8), dimension(:), allocatable :: velocity
        !! The velocity vector of the ice in the glacier at `location`
    end function velocity_func
  end interface

  public :: thickness_func, velocity_func

contains

  subroutine glacier_integrate(this, old_states, basal_melt, basal_drag, &
                               water_density, time)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Integrates the glacier's state forward in time by one time
    ! step. This is done using the NITSOL package of iterative Krylov
    ! solvers. If a different algorithm for the integration is
    ! desired, then this method may be overridden in the concrete
    ! implementations of the glacier type.
    !
    class(glacier), intent(inout)            :: this
    class(glacier), dimension(:), intent(in) :: old_states
      !! Previous states of the glacier, with the most recent one
      !! first.
    class(scalar_field), intent(in)          :: basal_melt
      !! The melt rate that the bottom of the glacier experiences
      !! during this time step.
    class(scalar_field), intent(in)          :: basal_drag
      !! A paramter, e.g. coefficient of friction, needed to calculate
      !! the drag on basal surface of the glacier.
    real(r8), intent(in)                     :: water_density
      !! The density of the water below the glacier.
    real(r8), intent(in)                     :: time
      !! The time to which the glacier should be integrated
    logical                                   :: first_call
    integer, save                             :: nval, kdmax = 20
    real(r8), dimension(:), allocatable       :: state
    integer, dimension(10)                    :: input
    integer, dimension(6)                     :: info
    real(r8), dimension(:), allocatable, save :: work
    real(r8), dimension(1)                    :: real_param
    integer, dimension(1)                     :: int_param
    integer                                   :: flag

    first_call = .true.
    if (.not. allocated(work)) then
      nval = this%data_size()
      allocate(work(nval*(kdmax+5) + kdmax*(kdmax+3)))
    end if
    state = this%state_vector()
    call this%set_time(time)
    input = 0
    input(4) = kdmax

    call nitsol(nval, state, nitsol_residual, dummy_jacv, 1.e-6_r8, &
                1.e-6_r8, input, info, work, real_param, int_param, &
                flag, ddot, dnrm2)

    select case(flag)
    case(0)
      write(*,*) 'Integrated glacier to time', time
    case(1)
      write(*,*) 'Reached maximum number of iterations integrating glacier'
      error stop
    case default
      write(*,*) 'NITSOL failed when integrating glacier with error code', flag
      error stop
    end select

  contains
    
    subroutine nitsol_residual(n, xcur, fcur, rpar, ipar, itrmf)
      !! A routine matching the interface expected by NITSOL which
      !! returns the residual for the glacier.
      integer, intent(in)                   :: n
        !! Dimension of the problem
      real(r8), dimension(*), intent(in)    :: xcur
        !! Array of length `n` containing the current \(x\) value
      real(r8), dimension(*), intent(out)   :: fcur
        !! Array of length `n` containing f(xcur) on output
      real(r8), dimension(*), intent(inout) :: rpar
        !! Parameter/work array
      integer, dimension(*), intent(inout)  :: ipar
        !! Parameter/work array
      integer, intent(out)                  :: itrmf
        !! Termination flag. 0 means normal termination, 1 means
        !! failure to produce f(xcur)

      ! If this is the first call of this routine then the
      ! basal_surface object will already be in the same state as
      ! reflected in xcur
      if (first_call) then
        first_call = .false.
      else
        call this%update(xcur(1:n))
      end if
      fcur(1:n) = this%residual(old_states,basal_melt,basal_drag,water_density)
      itrmf = 1
    end subroutine nitsol_residual

  end subroutine glacier_integrate

end module glacier_mod
