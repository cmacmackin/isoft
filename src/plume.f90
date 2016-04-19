!
!  plume.f90
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

module plume_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of the [[basal_surface]] data type,
  ! representing a buoyant plume beneath an ice shelf.
  !
  use iso_fortran_env, only: r8 => real64
  use basal_surface_mod, only: basal_surface
  use factual_mod, only: scalar_field, cheb1d_scalar_field, cheb1d_vector_field
  implicit none
  private

  type, extends(basal_surface), public :: plume
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A concrete implementation of the [[basal_surface]] abstract data type,
    ! representing the buoyant plume beneath an ice shelf.
    !
  contains
    procedure :: basal_melt => plume_melt
    procedure :: basal_drag_parameter => plume_drag_parameter
    procedure :: water_density => plume_water_density
    procedure :: residual => plume_residual
    procedure :: update => plume_update
  end type plume

contains

  function plume_melt(this) result(melt)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the melt rate at the bottom of the ice
    ! shelf due to interaction with the plume.
    !
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: melt
      !! The melt rate at the base of the ice shelf.
  end function plume_melt

  function plume_drag_parameter(this) result(drag)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns a quantity which may be necessary to determine
    ! the frictional drag the plume exerts on the bottom of the ice
    ! shelf. The plume would actually tend to exert no drag on the bottom
    ! of the ice shelf, but this method is present so that there is a
    ! consistent interface with the [[ground]] data type.
    !
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: drag
      !! The melt rate at the base of the ice sheet.
  end function plume_drag_parameter

  function plume_water_density(this) result(density)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes and returns the density of the plume water beneath the ice
    ! shelf. The density of this water would vary depending on how much 
    ! saline ambient water has been entrained into the plume versus how
    ! much fresh water has been released due to melting.
    !
    class(plume), intent(in)         :: this
    class(scalar_field), allocatable :: density
      !! The melt rate at the base of the ice sheet.
  end function plume_water_density
   
  function plume_residual(this, ice_thickness, ice_density, ice_temperature) &
                                                              result(residual)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Using the current state of the plume, this computes the residual
    ! of the system of equatiosn which is used to describe the it.
    ! The residual takes the form of a 1D array, with each element 
    ! respresenting the residual for one of the equations in the system.
    !
    class(plume), intent(in)            :: this
    class(scalar_field), intent(in)     :: ice_thickness
      !! Thickness of the ice above the plume.
    real(r8), intent(in)                :: ice_density
      !! The density of the ice above the plume, assumed uniform.
    real(r8), intent(in)                :: ice_temperature
      !! The temperature of the ice above the plume, assumed uniform.
    real(r8), dimension(:), allocatable :: residual
      !! The residual of the system of equations describing the plume.
  end function plume_residual

  subroutine plume_update(this, state_vector, time)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Updates the state of the plume from its state vector. The state
    ! vector is a real array containing the value of each of the plume's
    ! properties at each of the locations on the grid used in discretization.
    !
    class(plume), intent(inout)        :: this
    real(r8), dimension(:), intent(in) :: state_vector
      !! A real array containing the data describing the state of the
      !! plume.
    real(r8), intent(in), optional     :: time
      !! The time at which the plume is in this state. If not
      !! present then assumed to be same as previous value passed.
  end subroutine plume_update

end module plume_mod
