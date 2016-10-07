!
!  viscosity.f90
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

module viscosity_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides an abstract data type to model viscosity into a
  ! vertically integrated plume. 
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private

  type, abstract, public :: abstract_viscosity
    !* Author: Christopher MacMackin
    !  Date: October 2016
    !
    ! An abstract data type for calculating viscosity of a vertically
    ! integrated [[glacier]]. 
    !
  contains
    procedure(get_viscosity), deferred   :: ice_viscosity
      !! Returns the viscosity for the ice.
  end type abstract_viscosity

  abstract interface
    function get_viscosity(this, velocity, temperature, time) result(property)
      import :: abstract_viscosity
      import :: vector_field
      import :: scalar_field
      import :: r8
      class(abstract_viscosity), intent(in) :: this
      class(vector_field), intent(in) :: velocity
        !! The velocity field of the ice for which the velocity is
        !! being calculated
      real(r8), intent(in)            :: temperature
        !! The temperature of the ice for which viscosity is being
        !! calculated.
      real(r8), intent(in), optional  :: time
        !! The time at which the viscosity is being calculated. If not
        !! present then assumed to be same as previous value passed.
      class(scalar_field), allocatable :: property
        !! The value of the viscosity
    end function get_viscosity
  end interface

end module viscosity_mod
