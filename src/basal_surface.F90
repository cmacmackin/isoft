!
!  basal_surface.f90
!  This file is part of ISOFT.
!  
!  Copyright 2016 Chris MacMackin <cmacmackin@physics.ox.ac.uk>
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

module basal_surface_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract data type to model the ground or ocean below
  ! the glacier.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field
  use nitsol_mod, only: nitsol, dummy_jacv, ddot, dnrm2
  implicit none
  private

  type, abstract, public :: basal_surface
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! An abstract data type which represents whatever lies below a [[glacier]].
    ! This could be the ground, a plume, or a fully dynamic ocean model.
    ! Methods are available to provide the coupling information between the
    ! [[glacier]] and the basal surface.
    !
  contains
    procedure(get_scalar), deferred   :: basal_melt
      !! Returns the basal melt rate.
    procedure(get_scalar), deferred   :: basal_drag_parameter
      !! Returns a value which may be needed to calculate basal drag,
      !! such as the coefficient of friction.
    procedure(get_real), deferred     :: water_density
      !! Density of the water at the basal surface.
    procedure(get_residual), deferred :: residual
      !! Computes the residual of the system of equations describing the
      !! basal surface.
    procedure(setter), deferred       :: update
      !! Sets the state of the basal surface
    procedure(time_setter), deferred  :: set_time
      !! Sets the time record for this basal surface.
    procedure(get_i), deferred        :: data_size
      !! Returns the number of elements in the glacier's state vector
    procedure(get_r81d), deferred     :: state_vector
      !! Returns the glacier's state vector, a 1D array with all necessary 
      !! data to describe its state.
    procedure                         :: solve => basal_solve
  end type basal_surface

  abstract interface
    function get_scalar(this) result(property)
      import :: basal_surface
      import :: scalar_field
      class(basal_surface), intent(in) :: this
      class(scalar_field), allocatable :: property
        !! The value of whatever property of the basal surface is being
        !! returned.
    end function get_scalar

    function get_real(this) result(property)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(in) :: this
      real(r8)                         :: property
        !! The value of whatever property of the basal surface is being 
        !! returned.
    end function get_real
    
    function get_residual(this, ice_thickness, ice_density, ice_temperature) &
                                                            result(residual)
      import :: basal_surface
      import :: scalar_field
      import :: r8
      class(basal_surface), intent(inout) :: this
      class(scalar_field), intent(in)     :: ice_thickness
        !! Thickness of the ice above the basal surface
      real(r8), intent(in)                :: ice_density
        !! The density of the ice above the basal surface, assumed uniform
      real(r8), intent(in)                :: ice_temperature
        !! The temperature of the ice above the basal surface, assumed uniform
      real(r8), dimension(:), allocatable :: residual
        !! The residual of the system of equations describing the basal
        !! surface.
    end function get_residual

    pure function get_r81d(this) result(state_vector)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(in)          :: this
      real(r8), dimension(:), allocatable :: state_vector
        !! The state vector of the basal surface
    end function get_r81d

    subroutine setter(this, state_vector)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(inout) :: this
      real(r8), dimension(:), intent(in)  :: state_vector
        !! A real array containing the data describing the state of the
        !! basal surface.
    end subroutine setter

    subroutine time_setter(this, time)
      import :: basal_surface
      import :: r8
      class(basal_surface), intent(inout) :: this
      real(r8), intent(in)          :: time
        !! The time at which the basal surface is in the present state.
    end subroutine time_setter

    pure function get_i(this) result(property)
      import :: basal_surface
      class(basal_surface), intent(in) :: this
      integer                    :: property
        !! The value of whatever property of the basal surface is being
        !! returned.
    end function get_i
  end interface

contains

  subroutine basal_solve(this, ice_thickness, ice_density, &
                         ice_temperature, time)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Solves the state of the basal surface for the specified ice
    ! properties, at the specified time. This is done using the NITSOL
    ! package of iterative Krylov solvers.
    !
    class(basal_surface), intent(inout) :: this
    class(scalar_field), intent(in)     :: ice_thickness
      !! Thickness of the ice above the basal surface
    real(r8), intent(in)                :: ice_density
      !! The density of the ice above the basal surface, assumed uniform
    real(r8), intent(in)                :: ice_temperature
      !! The temperature of the ice above the basal surface, assumed uniform
    real(r8), intent(in)                :: time
      !! The time to which the basal surface should be solved
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
      write(*,*) 'Solved for plume at time', time
    case(1)
      write(*,*) 'Reached maximum number of iterations solving for plume'
      error stop
    case default
      write(*,*) 'NITSOL failed when solving plume with error code', flag
      error stop
    end select

  contains

    subroutine nitsol_residual(n, xcur, fcur, rpar, ipar, itrmf)
      !! A routine matching the interface expected by NITSOL which
      !! returns the residual for the basal surface.
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
      fcur(1:n) = this%residual(ice_thickness,ice_density,ice_temperature)
      itrmf = 1
    end subroutine nitsol_residual
  end subroutine basal_solve

end module basal_surface_mod
