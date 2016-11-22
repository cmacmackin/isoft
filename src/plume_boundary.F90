!
!  plume_boundary.f90
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

module plume_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: September 2016
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! specify the boundary conditions for plume types. Boundaries are
  ! classified according to their location as "north", "east", "south",
  ! and "west", where west is the lower boundary in the first direction,
  ! east the upper boundary in the first direction, south the lower 
  ! boundary in the second direction, and north the upper boundary in
  ! the second direction (see diagram below).
  !
  ! ![A diagram of the arrangement of the north, south, east, and west boudnaries](|media|/boundaries.svg)
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private
  
  type, public :: plume_boundary
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! A type in which procedures for getting the boundary
    ! conditions of plumes are to be specified. The descendent types
    ! can contain whatever data is needed to compute the result. For
    ! general boundary conditions, it provides the routine 
    ! [[plume_boundary(type):get_boundaries_residual]] to return an array
    ! with the residuals representing deviation from satisfying the
    ! conditions. This can then be appended to a [[plume(type)]]'s residual 
    ! array. Given that, depending on the numerical algorithm,
    ! computational savings can be made when a Dirichlet boundary is
    ! used, routines are also provided to indicated whether this is the case
    ! and to retrieve the constant boundary values.
    !
    ! This class effectively provides free boundary conditions. It's 
    ! type-bound procedures should be overridden to provide case-specific
    ! conditions.
    ! 
    ! Boundaries are classified according to their location as
    ! "north", "east", "south", and "west", where west is the lower
    ! boundary in the first direction, east the upper boundary in the
    ! first direction, south the lower boundary in the second direction,
    ! and north the upper boundary in the second direction (see diagram
    ! below).
    !
    ! ![A diagram of the arrangement of the north, south, east, and west boudnaries](|media|/boundaries.svg)
    !
  contains
    procedure :: thickness_north_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "north" boundary for thickness. Is unused for 1D simulations.
    procedure :: thickness_east_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "east" boundary for thickness
    procedure :: thickness_south_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "south" boundary for thickness. Is unused for 1D simulations.
    procedure :: thickness_west_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for thickness
    procedure :: velocity_north_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "north" boundary for velocity. Is unused for 1D simulations.
    procedure :: velocity_east_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "east" boundary for velocity
    procedure :: velocity_south_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "south" boundary for velocity. Is unused for 1D simulations.
    procedure :: velocity_west_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for velocity
    procedure :: temperature_north_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "north" boundary for temperature. Is unused for 1D simulations.
    procedure :: temperature_east_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "east" boundary for temperature
    procedure :: temperature_south_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "south" boundary for temperature. Is unused for 1D simulations.
    procedure :: temperature_west_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for temperature
    procedure :: salinity_north_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "north" boundary for salinity. Is unused for 1D simulations.
    procedure :: salinity_east_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "east" boundary for salinity
    procedure :: salinity_south_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "south" boundary for salinity. Is unused for 1D simulations.
    procedure :: salinity_west_is_dirichlet => is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for salinity
    procedure :: thickness_north => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the north boundary.
    procedure :: thickness_east => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the east boundary.
    procedure :: thickness_south => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the south boundary.
    procedure :: thickness_west => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the west boundary.
    procedure :: velocity_north => dirichlet_none_vec
      !! Returns the value of the Dirichlet boundary condition for velocity
      !! at the north boundary.
    procedure :: velocity_east => dirichlet_none_vec
      !! Returns the value of the Dirichlet boundary condition for velocity
      !! at the east boundary.
    procedure :: velocity_south => dirichlet_none_vec
      !! Returns the value of the Dirichlet boundary condition for velocity
      !! at the south boundary.
    procedure :: velocity_west => dirichlet_none_vec
      !! Returns the value of the Dirichlet boundary condition for velocity
      !! at the west boundary.
    procedure :: temperature_north => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for temperature
      !! at the north boundary.
    procedure :: temperature_east => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for temperature
      !! at the east boundary.
    procedure :: temperature_south => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for temperature
      !! at the south boundary.
    procedure :: temperature_west => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for temperature
      !! at the west boundary.
    procedure :: salinity_north => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for salinity
      !! at the north boundary.
    procedure :: salinity_east => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for salinity
      !! at the east boundary.
    procedure :: salinity_south => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for salinity
      !! at the south boundary.
    procedure :: salinity_west => dirichlet_none
      !! Returns the value of the Dirichlet boundary condition for salinity
      !! at the west boundary.
    procedure :: boundary_residuals
      !! Returns an array consisting of the difference between the required
      !! boundary values and those which actually exist. This can then be
      !! appended to a plume's state vector
    procedure :: residual_size
      !! Provides an integer indicating the size of the array returned by
      !! [[plume_boundary(type):boundary_residuals]]. 
  end type plume_boundary

  abstract interface
  end interface

contains
 
  function is_dirichlet(this)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Defualt implementation of the `_is_dirichlet` routines, indicating
    ! that Dirichlet boundary conditions are _not_ applied.
    !
    class(plume_boundary), intent(in) :: this
    logical :: is_dirichlet
      !! Always returns `.false` in the default implementation.
    is_dirichlet = .false.
    return
  end function is_dirichlet

  function dirichlet_none(this, t)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Defualt implementation of the routines to get scalar Dirichlet
    ! boundary conditions, which returns an uninitialised variable as they 
    ! should be overriden if the result is to actually be used.
    !
    class(plume_boundary), intent(in)   :: this
    real(r8), intent(in)                :: t
      !! The time at which the boundary condition is to be calculated.
    real(r8) :: dirichlet_none
      !! The value of the boundary conditions
  end function dirichlet_none

  function dirichlet_none_vec(this, t)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Defualt implementation of the routines to get vector Dirichlet
    ! boundary conditions, which returns an uninitialised variable, as they 
    ! should be overriden if the result is to actually be used.
    !
    class(plume_boundary), intent(in)   :: this
    real(r8), intent(in)                :: t
      !! The time at which the boundary condition is to be calculated.
    real(r8), allocatable, dimension(:) :: dirichlet_none_vec
      !! The value of the boundary conditions
  end function dirichlet_none_vec

  function boundary_residuals(this, thickness, velocity, temperature, &
                              salinity, t) result(residuals)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Default implementation of the [[plume_boundary(type):boundary_residuals]]
    ! method. It returns a zero-length array, effectively indicating free
    ! boundaries.
    !
    class(plume_boundary), intent(in)   :: this
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
    allocate(residuals(0))
    return
  end function boundary_residuals

  function residual_size(this)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Default implementation of the [[plume_boundary(type):residual_size]]
    ! routine providing the size of the residual array returned by
    ! [[plume_boundary(type):boundary_residuals]]. For performance
    ! reasons, this routine does not simply use the `size()`
    ! intrinisic on the residual array and should be implemented for
    ! each, boundary condition sub-type. Preferably it will return a
    ! constant.
    !
    class(plume_boundary), intent(in) :: this
    integer :: residual_size
      !! The number of elements in the array returned by
      !! [[plume_boundary(type):boundary_residuals]]
    residual_size = 0
    return
  end function residual_size

end module plume_boundary_mod
