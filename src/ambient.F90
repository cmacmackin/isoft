!
!  ambient.f90
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

module ambient_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! specify the temperature and salinity of the ambient ocean.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field
  implicit none
  private
  
  type, abstract, public :: ambient_conditions
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! An abstract type to which procedures for getting the ambient ocean
    ! conditions are to be specified. The descendent types can contain
    ! whatever data is needed to compute the result.
    !
  contains
    procedure(get_property), deferred :: ambient_temperature
      !! Returns the ambient ocean temperature
    procedure(get_property), deferred :: ambient_salinity
      !! Returns the ambient ocean temperature
  end type ambient_conditions

  abstract interface
    function get_property(this, depth, t) result(property)
      import :: r8
      import :: ambient_conditions
      import :: scalar_field
      class(ambient_conditions), intent(in) :: this
      class(scalar_field), intent(in)       :: depth
        !! A field containing the depths at which the ambient conditions
        !! are to be calculated.
      real(r8), intent(in)                  :: t
        !! The time at which the ambient conditions are to be calculated.
      class(scalar_field), pointer          :: property
        !! A field containing the ambient conditions at the depth specified
        !! for each location.
    end function get_property
  end interface

end module ambient_mod
