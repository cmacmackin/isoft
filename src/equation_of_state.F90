!
!  equation_of_state.f90
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

module equation_of_state_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! implement an equation of state.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field
  implicit none
  private
  
  type, abstract, public :: equation_of_state
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! An abstract type with a procedure for calculating water density
    ! from its temperature and salinity.
    !
  contains
    procedure(get_property), deferred :: water_density
      !! Returns the water density for the given temperature and salinity.
 end type equation_of_state

  abstract interface
    pure function get_property(this, temperature, salinity) result(density)
      import :: r8
      import :: equation_of_state
      import :: scalar_field
      class(equation_of_state), intent(in) :: this
      class(scalar_field), intent(in)      :: temperature
        !! A field containing the temperature of the water
      class(scalar_field), intent(in)      :: salinity
        !! A field containing the salinity of the water
      class(scalar_field), allocatable     :: density
        !! A field containing the density of the water
    end function get_property
  end interface

end module equation_of_state_mod