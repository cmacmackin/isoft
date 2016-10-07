!
!  entrainment.f90
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

module entrainment_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides an abstract data type to model entrainment into a
  ! vertically integrated plume. 
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  implicit none
  private

  type, abstract, public :: abstract_entrainment
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! An abstract data type for calculating entrainment of ambient
    ! ocean water into a vertically integrated [[plume]].
    !
  contains
    procedure(get_entrainment), deferred   :: entrainment_rate
      !! Returns the entrainment rate for ambient water into the plume.
  end type abstract_entrainment

  abstract interface
    function get_entrainment(this, velocity, thickness, depth, time) &
                                                    result(property)
      import :: abstract_entrainment
      import :: vector_field
      import :: scalar_field
      import :: r8
      class(abstract_entrainment), intent(in) :: this
      class(vector_field), intent(in)  :: velocity
        !! The velocity field of the plume into which fluid is being 
        !! entrained.
      class(scalar_field), intent(in)  :: thickness
        !! The thickness of the plume into which fluid is being
        !! entrained
      class(scalar_field), intent(in)  :: depth
        !! The depth of the upper surface of the plume into which
        !! fluid is being entrained
      real(r8), intent(in), optional   :: time
        !! The time at which the entrainment is being calculated. If not
        !! present then assumed to be same as previous value passed.
      class(scalar_field), allocatable :: property
        !! The value of the entrainment
    end function get_entrainment
  end interface

end module entrainment_mod
