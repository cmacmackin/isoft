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
    procedure(get_property_dx), deferred :: water_density_derivative
      !! Returns the derivative of the water density for the given
      !! temperature and salinity.
    procedure(get_coef), deferred    :: haline_contraction
      !! Returns a (possibly approximated) haline contraction coefficient.
    procedure(get_coef), deferred    :: thermal_contraction
      !! Returns a (possibly approximated) therma contraction coefficient.
 end type equation_of_state

  abstract interface
    function get_property(this, temperature, salinity) result(density)
      import :: r8
      import :: equation_of_state
      import :: scalar_field
      class(equation_of_state), intent(in) :: this
      class(scalar_field), intent(in)      :: temperature
        !! A field containing the temperature of the water
      class(scalar_field), intent(in)      :: salinity
        !! A field containing the salinity of the water
      class(scalar_field), pointer         :: density
        !! A field containing the density of the water
    end function get_property

    function get_property_dx(this, temperature, d_temperature, salinity, &
                             d_salinity, dir) result(d_density)
      import :: r8
      import :: equation_of_state
      import :: scalar_field
      class(equation_of_state), intent(in) :: this
      class(scalar_field), intent(in)      :: temperature
        !! A field containing the temperature of the water
      class(scalar_field), intent(in)      :: d_temperature
        !! A field containing the derivative of the temperature of the
        !! water, in teh same direction as `dir`
      class(scalar_field), intent(in)      :: salinity
        !! A field containing the salinity of the water
      class(scalar_field), intent(in)      :: d_salinity
        !! A field containing the derivative of the salinity of the
        !! water, in the same direction as `dir`
      integer, intent(in)                  :: dir
        !! The direction in which to take the derivative
      class(scalar_field), pointer         :: d_density
        !! A field containing the derivative of the density of the
        !! water in direction `dir`
    end function get_property_dx

    function get_coef(this, temperature, salinity) result(coef)
      import :: r8
      import :: equation_of_state
      import :: scalar_field
      class(equation_of_state), intent(in) :: this
      class(scalar_field), intent(in)      :: temperature
      class(scalar_field), intent(in)      :: salinity
      class(scalar_field), allocatable     :: coef
    end function get_coef
  end interface

end module equation_of_state_mod
