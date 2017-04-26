!
!  glens_law.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@physics.ox.ac.uk>
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

module glens_law_mod
  !* Author: Christopher MacMackin
  !  Date: April 2017
  !  License: GPLv3
  !
  ! Provides a concrete implementation for the [[abstract_viscosity]]
  ! type using Glen's flow law.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field
  use viscosity_mod, only: abstract_viscosity
  implicit none
  private

  type, extends(abstract_viscosity), public :: glens_law_viscosity
    !* Author: Christopher MacMackin
    !  Date: April 2017
    !
    ! An implementation of Glen's flow law to describe glacier
    ! viscosity. It takes the form $$\eta = \frac{1}{2}BD^{1/n-1},$$
    ! where \(D = \sqrt{D_{ij}D_{ij}/2\) is the second invarient of
    ! the strain rate $$D_{ij} = \frac{1}{2}\left(\frac{\partial
    ! u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}
    ! \right).$$ Here, \(B\) is treated as a constant, although it may
    ! be parameterised as a function of temperature.
    !
    private
    real(r8) :: b_val = 1.0_r8
    real(r8) :: index = 3._r8
  contains
    procedure :: ice_viscosity => glens_ice_viscosity
      !! Returns the viscosity for the ice.
  end type glens_law_viscosity

  interface glens_law_viscosity
    module procedure constructor
  end interface glens_law_viscosity

contains

  pure function constructor(b_val, index) result(this)
    !* Author: Christopher MacMackin
    !  Date: April 2017
    !
    ! Instantiates an instance of a viscosity object implementing
    ! Glen's flow law.
    !
    real(r8), intent(in) :: b_val
      !! The coefficient, \(B\), in Glen's flow law.
    real(r8), intent(in) :: index
      !! The index, \(n\), in the exponent of Glen's flow law.
    type(glens_law_viscosity) :: this
      !! The viscosity object being created.
    this%b_val = b_val
    this%index = index
  end function constructor

  function glens_ice_viscosity(this, velocity, temperature, time) &
                                                    result(viscosity)
    !* Author: Christopher MacMackin
    !  Date: April 2017
    !
    ! Calculates the viscosity of ice using Glen's flow law. See the
    ! documentation of the [[glens_law_viscosity]] object for a
    ! description of this parameterisation.
    !
    class(glens_law_viscosity), intent(in) :: this
    class(vector_field), intent(in)        :: velocity
      !! The velocity field of the ice for which the velocity is
      !! being calculated
    real(r8), intent(in)                   :: temperature
      !! The temperature of the ice for which viscosity is being
      !! calculated.
    real(r8), intent(in), optional         :: time
      !! The time at which the viscosity is being calculated. If not
      !! present then assumed to be same as previous value passed.
    class(scalar_field), allocatable       :: viscosity
      !! The value of the viscosity
    call velocity%guard_temp()

    call velocity%clean_temp()
    call viscosity%set_temp()
  end function glens_ice_viscosity

end module glens_law_mod
