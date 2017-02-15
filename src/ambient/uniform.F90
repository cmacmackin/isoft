!
!  uniform.f90
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

module uniform_ambient_mod
  !* Author: Christopher MacMackin
  !  Date: November 2016
  !  License: GPLv3
  !
  ! Provides a derived type specifying uniform ambient temperature and
  ! salinity conditions beneath an ice shelf.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, uniform_scalar_field
  use ambient_mod, only: ambient_conditions
  implicit none
  private

  
  type, extends(ambient_conditions), public :: uniform_ambient_conditions
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! An derived type with procedures for getting the ambient ocean
    ! conditions. This implementation takes these conditions to be
    ! everywhere uniform.
    !
    private
    type(uniform_scalar_field) :: temperature
    type(uniform_scalar_field) :: salinity
  contains
    procedure :: ambient_temperature => uniform_temperature
      !! Returns the ambient ocean temperature
    procedure :: ambient_salinity => uniform_salinity
      !! Returns the ambient ocean temperature
  end type uniform_ambient_conditions

  interface uniform_ambient_conditions
    module procedure constructor
  end interface uniform_ambient_conditions

contains

  function constructor(temperature, salinity) result(this)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Produces an ambient object which will return the specified
    ! salinity and temeprature values.
    !
    real(r8), intent(in), optional :: temperature
      !! The temperature of the ambient ocean. Default is 0.
    real(r8), intent(in), optional :: salinity
      !! The salinity of the ambient ocean. Default is 0.
    type(uniform_ambient_conditions) :: this
    if (present(temperature)) then
      this%temperature = uniform_scalar_field(temperature)
    else
      this%temperature = uniform_scalar_field(0.0_r8)
    end if
    if (present(salinity)) then
      this%salinity = uniform_scalar_field(salinity)
    else
      this%salinity = uniform_scalar_field(0.0_r8)
    end if
  end function constructor

  pure function uniform_temperature(this, depth, t) result(property)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the ambient ocean temperature.
    !
    class(uniform_ambient_conditions), intent(in) :: this
    class(scalar_field), intent(in)               :: depth
      !! A field containing the depths at which the ambient temperature
      !! is to be calculated.
    real(r8), intent(in)                          :: t
      !! The time at which the ambient conditions are to be calculated.
    class(scalar_field), allocatable              :: property
      !! A field containing the ambient temperature at the depth specified
      !! for each location.
    call this%temperature%allocate_scalar_field(property)
    property = this%temperature
    call property%set_temp()
  end function uniform_temperature

  pure function uniform_salinity(this, depth, t) result(property)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Returns the ambient ocean salinity.
    !
    class(uniform_ambient_conditions), intent(in) :: this
    class(scalar_field), intent(in)               :: depth
      !! A field containing the depths at which the ambient salinity
      !! is to be calculated.
    real(r8), intent(in)                          :: t
      !! The time at which the ambient conditions are to be calculated.
    class(scalar_field), allocatable              :: property
      !! A field containing the ambient salinity at the depth specified
      !! for each location.
    call this%salinity%allocate_scalar_field(property)
    property = this%salinity
    call property%set_temp()
  end function uniform_salinity

end module uniform_ambient_mod
