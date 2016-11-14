!
!  linear_eos.f90
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

module linear_eos_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! implement an equation of state.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field
  use equation_of_state_mod, only: equation_of_state
  implicit none
  private

  real(r8), parameter :: absolute_zero = -273.15_r8
  
  type, extends(equation_of_state), public :: linear_eos
    !* Author: Chris MacMackin
    !  Date: April 2016
    !
    ! A linearised implementation of the equation of state, of the
    ! form $$ \rho = \rho_0[1-\beta_T(T-T_0) + \beta_S(S-S_0)]. $$
    !
    private
    real(r8) :: ref_rho = 1.0_r8
      !! The density for the temperature and salinity about which the
      !! equation of state was linearised, \(\rho_0\).
    real(r8) :: ref_t = 0.0_r8
      !! The temperature about which the equation of state was
      !! linearised, \(T_0\).
    real(r8) :: ref_s = 0.0_r8
      !! The salinity about which the equation of state was
      !! linearised, \(S_0\).
    real(r8) :: beta_t = 0.0_r8
      !! The thermal contraction coefficient, \(\beta_T\).
    real(r8) :: beta_s = 1.0_r8
      !! The haline contraction coefficient, \(\beta_S\).
  contains
    procedure :: water_density => linear_water_density
  end type linear_eos

  interface linear_eos
    module procedure constructor
  end interface linear_eos

contains

  pure function constructor(ref_rho, ref_t, ref_s, beta_t, beta_s) result(this)
    real(r8), intent(in) :: ref_rho
      !! The density for the temperature and salinity about which the
      !! equation of state was linearised, \(\rho_0\).
    real(r8), intent(in) :: ref_t
      !! The temperature about which the equation of state was
      !! linearised, \(T_0\).
    real(r8), intent(in) :: ref_s
      !! The salinity about which the equation of state was
      !! linearised, \(S_0\).
    real(r8), intent(in) :: beta_t
      !! The thermal contraction coefficient, \(\beta_T\).
    real(r8), intent(in) :: beta_s
      !! The haline contraction coefficient, \(\beta_S\).
    type(linear_eos)     :: this
    this%ref_rho = ref_rho
    this%ref_t   = ref_t
    this%ref_s   = ref_s
    this%beta_t  = beta_t
    this%beta_s  = beta_s
  end function constructor

  pure function linear_water_density(this, temperature, salinity) result(density)
    !* Author: Chris MacMackin
    !  Date: November 2016
    !
    ! Calculates the density of the water from the temperature and
    ! salinity, using a linear equatino of state, $$ \rho =
    ! \rho_0[1-\beta_T(T-T_0) + \beta_S(S-S_0)]. $$
    class(linear_eos), intent(in)    :: this
    class(scalar_field), intent(in)  :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)  :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), allocatable :: density
      !! A field containing the density of the water
    call salinity%allocate_scalar_field(density)
    density = this%ref_rho * (1.0_r8 - this%beta_t*(temperature - this%ref_t) &
                                     + this%beta_s*(salinity - this%ref_s))
  end function linear_water_density

end module linear_eos_mod
