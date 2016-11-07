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

module dallaston2015_glacier_boundary_mod
  !* Author: Christopher MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Provides a derived type which specifies the boundary conditions
  ! for the ice shelf model used by Dallaston et al. (2015). As this
  ! is a 1-D model, only west (incoming) and east (outgoing)
  ! boundaries are specified.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use glacier_boundary_mod, only: glacier_boundary
  implicit none
  private
  
  type, extends(glacier_boundary), public :: dallaston2015_glacier_boundary
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! A type with procedures for getting the boundary conditions of
    ! the ice shelf model used by Dallaston et al. (2015). These are
    ! Dirichlet conditions at the west (incoming) boundary, a
    ! thickness of 0 at the east (outgoing) boundary, and a free
    ! velocity at the east boundary.
    !
    ! @TODO Consider testing the consistency conditions at the outgoing boundary?
    !
    private
    real(r8) :: thickness = 1.0_r8
      !! The thickness of the glacier at the inflowing boundary
    real(r8) :: velocity = 1.0_r8
      !! The velocity of the glacier at the inflowing boundary
  contains
    procedure :: thickness_west_is_dirichlet => dallaston2015_is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for thickness
    procedure :: velocity_west_is_dirichlet => dallaston2015_is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "west" boundary for velocity
    procedure :: thickness_east_is_dirichlet => dallaston2015_is_dirichlet
      !! Returns `.true.` if constant Dirichlet conditions are used at the
      !! "east" boundary for thickness
    procedure :: thickness_west => dallaston2015_thickness_west
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the west boundary.
    procedure :: velocity_west => dallaston2015_velocity
      !! Returns the value of the Dirichlet boundary condition for velocity
      !! at the west boundary.
    procedure :: thickness_east => dallaston2015_thickness_east
      !! Returns the value of the Dirichlet boundary condition for thickness
      !! at the east boundary.
  end type dallaston2015_glacier_boundary

  interface dallaston2015_glacier_boundary
    module procedure constructor
  end interface dallaston2015_glacier_boundary

contains

  pure function constructor(thickness, velocity) result(this)
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
    type(dallaston2015_glacier_boundary) :: this
    this%thickness = thickness
    this%velocity = velocity
  end function constructor
 
  function dallaston2015_is_dirichlet(this)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Specifies Dirichlet inflow (west) boundaries for ice shelf
    ! thickness, velocity, temperature, and salinity.
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    logical :: dallaston2015_is_dirichlet
      !! Returns `.true.`
    dallaston2015_is_dirichlet = .true.
    return
  end function dallaston2015_is_dirichlet

  function dallaston2015_thickness_west(this, t) result(boundary_value)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Specifies the boundary value for thickness at the west boundary.
    !
    class(dallaston2015_glacier_boundary), intent(in)   :: this
    real(r8), intent(in)                :: t
      !! The time at which the boundary condition is to be calculated.
    real(r8) :: boundary_value
      !! The value of the boundary conditions
    boundary_value = this%thickness
  end function dallaston2015_thickness_west

  function dallaston2015_thickness_east(this, t) result(boundary_value)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Specifies the boundary value for thickness at the west boundary.
    !
    class(dallaston2015_glacier_boundary), intent(in)   :: this
    real(r8), intent(in)                :: t
      !! The time at which the boundary condition is to be calculated.
    real(r8) :: boundary_value
      !! The value of the boundary conditions
    boundary_value = 0.0_r8
  end function dallaston2015_thickness_east

  function dallaston2015_velocity(this,t) result(boundary_value)
    !* Author: Chris MacMackin
    !  Date: September 2016
    !
    ! Specifies the boundary value for velocity at the west boundary
    !
    class(dallaston2015_glacier_boundary), intent(in) :: this
    real(r8), intent(in)                :: t
      !! The time at which the boundary condition is to be calculated.
    real(r8), allocatable, dimension(:) :: boundary_value
      !! The value of the boundary conditions
    allocate(boundary_value(1))
    boundary_value = this%velocity
  end function dallaston2015_velocity

end module dallaston2015_glacier_boundary_mod
