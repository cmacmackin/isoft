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
  ! specify the boundary conditions for plume types.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field, &
                         uniform_vector_field
  use boundary_types_mod, only: free_boundary
  implicit none
  private
  
  type, public :: plume_boundary
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! A type in which procedures for getting the boundary conditions
    ! of plumes are to be specified. The descendent types can contain
    ! whatever data is needed to compute the result.
    !
    ! This class effectively provides free boundary conditions. It's 
    ! type-bound procedures should be overridden to provide case-specific
    ! conditions.
    !
  contains
    procedure :: thickness_bound_info => bound_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: velocity_bound_info => bound_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: temperature_bound_info => bound_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: salinity_bound_info => bound_info
      !! Indicates the type and depth of the thickness boundary at
      !! different locations.
    procedure :: thickness_bound => scalar_bound
      !! Produces a field containing the boundary conditions for plume
      !! thickness at the specified location.
    procedure :: velocity_bound => vector_bound
      !! Produces a field containing the boundary conditions for plume
      !! velocity at the specified location.
    procedure :: temperature_bound => scalar_bound
      !! Produces a field containing the boundary conditions for plume
      !! temperature at the specified location.
    procedure :: salinity_bound => scalar_bound
      !! Produces a field containing the boundary conditions for plume
      !! salinity at the specified location.
    procedure :: set_time
      !! Specifies the time at which to calculate the boundary
      !! conditions.
  end type plume_boundary

contains

  subroutine bound_info(this, location, bound_type, bound_depth)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Provides information about the type of boundary, and the number
    ! of layers of data-points needed to describe it.
    !
    class(plume_boundary), intent(in) :: this
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
    bound_type = free_boundary
    bound_depth = 0
  end subroutine bound_info

  function scalar_bound(this, location)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the boundary values for the specified
    ! boundary location.
    !
    class(plume_boundary), intent(in) :: this
    integer, intent(in)               :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(scalar_field), allocatable :: scalar_bound
    allocate(uniform_scalar_field :: scalar_bound)
    scalar_bound = uniform_scalar_field(0.0_r8)
  end function scalar_bound

  function vector_bound(this, location)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns a field containing the boundary values for the specified
    ! boundary location.
    !
    class(plume_boundary), intent(in) :: this
    integer, intent(in)               :: location
      !! Which boundary information is to be provided for.  The
      !! boundary will be the one normal to dimension of number
      !! `abs(boundary)`. If the argument is negative, then the lower
      !! boundary is returned. If positive, then the upper boundary is
      !! returned.
    class(vector_field), allocatable  :: vector_bound
    allocate(uniform_vector_field :: vector_bound)
    vector_bound = uniform_vector_field([0.0_r8,0.0_r8])
  end function vector_bound

  subroutine set_time(this, time)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Sets the time at which boundary conditions are to be calculated.
    !
    class(plume_boundary), intent(inout) :: this
    real(r8), intent(in)                 :: time
  end subroutine set_time

end module plume_boundary_mod
