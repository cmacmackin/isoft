!
!  one_equation_melt.f90
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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module one_equation_melt_mod
  !* Author: Christopher MacMackin
  !  Date: May 2017
  !  License: GPLv3
  !
  ! Provides an implementation of melt similar to that used by
  ! Dallaston, Hewitt, and Wells (2015), prior to their neglecting
  ! certain terms on scaling arguments.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use melt_relationship_mod, only: abstract_melt_relationship
  implicit none
  private

  type, extends(abstract_melt_relationship), public :: one_equation_melt
    !* Author: Christopher MacMackin
    !  Date: May 2017
    !
    ! A parameterisation of melting into a plume which comes from
    ! heavily simplifying the 3 equation model. It is taken from
    ! Dallaston, Hewitt, and Wells (2015), prior to the their dropping
    ! some terms based on scaling arguments. The melt rate, as well as
    ! effect on termperature and salinity, are calculated by calling
    ! [[abstract_melt_relationship:solve_for_melt]] and then accessed
    ! using [[abstract_melt_relationship:melt_rate]],
    ! [[abstract_melt_relationship:heat_equation_terms]],
    ! [[abstract_melt_relationship:salt_equation_terms]].
    ! 
    class(scalar_field), allocatable :: forcing_values
      !! Stores the resulting forcing values.
    real(r8) :: coef1 = 0.018208_r8
      !! The unitless multiplier on the thermal forcing term,
      !! \(\Gamma_Tx_0/D_0\).
    real(r8) :: coef2 = 0.023761_r8
      !! The unitless multiplier applied to the theram forcing term to
      !! get the melt rate, \(c_oT_0/L\).
  contains
    procedure :: solve_for_melt => one_equation_solve
    procedure :: heat_equation_terms => one_equation_heat
      !! Returns the terms this melt formulation contributes to the
      !! heat equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure :: salt_equation_terms => one_equation_salt
      !! Returns the terms this melt formulation contributes to the
      !! salt equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]].
    procedure :: melt_rate => one_equation_melt_rate
      !! Returns the melt rate calculated using this formulation,
      !! after it has been solved for using 
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure :: has_heat_terms => one_equation_has_heat
      !! Whether this formulation of melting contributes any terms to
      !! a plume's heat equation.
    procedure :: has_salt_terms => one_equation_has_salt
      !! Whether this formulation of melting contributes any terms to
      !! a plume's salinity equation.
 end type one_equation_melt

  interface one_equation_melt
    module procedure constructor
  end interface one_equation_melt

contains

  pure function constructor(coef1, coef2) result(this)
    real(r8), intent(in) :: coef1
      !! The unitless multiplier on the thermal forcing term,
      !! \(\Gamma_Tx_0/D_0\).
    real(r8), intent(in) :: coef2
      !! The unitless multiplier applied to the theram forcing term to
      !! get the melt rate, \(c_oT_0/L\).
    type(one_equation_melt) :: this
      !! The newly created object representing the melt relationship.
    this%coef1 = -coef1
    this%coef2 = -coef2
  end function constructor

  subroutine one_equation_solve(this, velocity, pressure, temperature, &
                                 salinity, plume_thickness, time)
    class(one_equation_melt), intent(inout) :: this
    class(vector_field), intent(in)         :: velocity
      !! The velocity field of the plume into which fluid is melting.
    class(scalar_field), intent(in)         :: pressure
      !! The water pressure at the interface where the melting occurs.
    class(scalar_field), intent(in)         :: temperature
      !! The temperature of the plume into which fluid is melting.
    class(scalar_field), intent(in)         :: salinity
      !! The salinity of the plume into which fluid is melting.
    class(scalar_field), intent(in)         :: plume_thickness
      !! The thickness of the plume into which fluid is melting.
    real(r8), intent(in), optional          :: time
      !! The time at which the melting is being solved for. If not
      !! present then assumed to be same as previous value passed.
    call velocity%guard_temp(); call pressure%guard_temp()
    call temperature%guard_temp(); call salinity%guard_temp()
    call plume_thickness%guard_temp()
    if (.not. allocated(this%forcing_values)) then
      allocate(this%forcing_values, mold=temperature)
    end if
    this%forcing_values = this%coef1*temperature*velocity%norm()
    call velocity%clean_temp(); call pressure%clean_temp()
    call temperature%clean_temp(); call salinity%clean_temp()
    call plume_thickness%clean_temp()
  end subroutine one_equation_solve

  function one_equation_heat(this) result(heat)
    class(one_equation_melt), intent(in) :: this
    class(scalar_field), pointer         :: heat
      !! The value of the contribution made by melting/thermal
      !! transfer to the heat equation for a [[plume]]
    if (.not. allocated(this%forcing_values)) error stop ('Melt values not calculated')
    call this%forcing_values%allocate_scalar_field(heat)
    heat = this%forcing_values
  end function one_equation_heat

  function one_equation_salt(this) result(salt)
    class(one_equation_melt), intent(in) :: this
    class(scalar_field), pointer         :: salt
      !! The value of the contribution made by melting/thermal
      !! transfer to the salt equation for a [[plume]]
  end function one_equation_salt

  function one_equation_melt_rate(this) result(melt)
    class(one_equation_melt), intent(in) :: this
    class(scalar_field), pointer         :: melt
      !! The melt rate from the ice into the plume water.
    if (.not. allocated(this%forcing_values)) error stop ('Melt values not calculated')
    call this%forcing_values%allocate_scalar_field(melt)
    melt = this%coef2 * this%forcing_values
  end function one_equation_melt_rate

  pure function one_equation_has_heat(this) result(has_heat)
    class(one_equation_melt), intent(in) :: this
    logical                              :: has_heat
      !! Whether this formulation of melting contributes terms to
      !! the heat equation of the plume.
    has_heat = .true.
  end function one_equation_has_heat

  pure function one_equation_has_salt(this) result(has_salt)
    class(one_equation_melt), intent(in) :: this
    logical                              :: has_salt
      !! Whether this formulation of melting contributes terms to
      !! the salinity equation of the plume.
    has_salt = .false.
  end function one_equation_has_salt

end module one_equation_melt_mod
