!
!  dallaston2015_melt.f90
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

module dallaston2015_melt_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides an implementation of [[abstract_melt_relationship]] which
  ! mimics the simple model used by Dallaston, Hewitt, and Wells
  ! (2015) for an ice shelf melting into a vertically integrated
  ! plume.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use melt_relationship_mod, only: abstract_melt_relationship
  implicit none
  private

  type, extends(abstract_melt_relationship), public :: dallaston2015_melt
    !* Author: Christopher MacMackin
    !  Date: October 2016
    !
    ! A parameterisation of melting into a plume which comes from
    ! heavily simplifying the 3 equation model. It is taken from
    ! Dallaston, Hewitt, and Wells (2015). The melt rate, as well as
    ! effect on termperature and salinity, are calculated by calling
    ! [[abstract_melt_relationship:solve_for_melt]] and then accessed
    ! using [[abstract_melt_relationship:melt_rate]],
    ! [[abstract_melt_relationship:heat_equation_terms]],
    ! [[abstract_melt_relationship:salt_equation_terms]].
    ! 
    class(scalar_field), allocatable :: melt_values
      !! Stores the resulting melt rate
    real(r8) :: coef = 1449.29936
      !! The coefficient by which the melt rate is multiplied in order
      !! to determine the contribution to the heat equation.
    real(r8) :: melt_conversion = 6.9e-4_r8
      !! The factor to convert between the scale for melt used by
      !! Dallaston et al. (2015) and that used in ISOFT, $$
      !! \frac{m_0x_0}{D_0U_0}, $$ where \(m_0\) is the melt scale
      !! used by Dalalston et al.
  contains
    procedure :: solve_for_melt => dallaston2015_solve
    procedure :: heat_equation_terms => dallaston2015_heat
      !! Returns the terms this melt formulation contributes to the
      !! heat equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure :: salt_equation_terms => dallaston2015_salt
      !! Returns the terms this melt formulation contributes to the
      !! salt equation, after they have been solved for using
      !! [[abstract_melt_relationship:solve_for_melt]].
    procedure :: melt_rate => dallaston2015_melt_rate
      !! Returns the melt rate calculated using this formulation,
      !! after it has been solved for using 
      !! [[abstract_melt_relationship:solve_for_melt]]. 
    procedure :: has_heat_terms => dallaston2015_has_heat
      !! Whether this formulation of melting contributes any terms to
      !! a plume's heat equation.
    procedure :: has_salt_terms => dallaston2015_has_salt
      !! Whether this formulation of melting contributes any terms to
      !! a plume's salinity equation.
  end type dallaston2015_melt

  interface dallaston2015_melt
    module procedure constructor
  end interface dallaston2015_melt

contains

  pure function constructor(beta, melt_conversion) result(this)
    real(r8), intent(in) :: beta
      !! The inverse stefan number, $$ \frac{c(T_a - T_m)}{L} $$
    real(r8), intent(in) :: melt_conversion
      !! The factor to convert between the scale for melt used by
      !! Dallaston et al. (2015) and that used in ISOFT, $$
      !! \frac{m_0x_0}{D_0U_0}, $$ where \(m_0\) is the melt scale
      !! used by Dalalston et al.
    type(dallaston2015_melt) :: this
      !! The newly created object representing the melt relationship.
    this%coef = (beta + 1.0_r8)
    this%melt_conversion = melt_conversion
  end function constructor

  subroutine dallaston2015_solve(this, velocity, pressure, temperature, &
                                 salinity, plume_thickness, time)
    use factual_mod
    class(dallaston2015_melt), intent(inout) :: this
    class(vector_field), intent(in)          :: velocity
      !! The velocity field of the plume into which fluid is melting.
    class(scalar_field), intent(in)          :: pressure
      !! The water pressure at the interface where the melting occurs.
    class(scalar_field), intent(in)          :: temperature
      !! The temperature of the plume into which fluid is melting.
    class(scalar_field), intent(in)          :: salinity
      !! The salinity of the plume into which fluid is melting.
    class(scalar_field), intent(in)          :: plume_thickness
      !! The thickness of the plume into which fluid is melting.
    real(r8), intent(in), optional           :: time
      !! The time at which the melting is being solved for. If not
      !! present then assumed to be same as previous value passed.
    call velocity%guard_temp(); call pressure%guard_temp()
    call temperature%guard_temp(); call salinity%guard_temp()
    call plume_thickness%guard_temp()
    if (.not. allocated(this%melt_values)) then
      call velocity%allocate_scalar_field(this%melt_values)
    end if
    this%melt_values = velocity%norm()
    call velocity%clean_temp(); call pressure%clean_temp()
    call temperature%clean_temp(); call salinity%clean_temp()
    call plume_thickness%clean_temp()
  end subroutine dallaston2015_solve

  pure function dallaston2015_heat(this) result(heat)
    class(dallaston2015_melt), intent(in) :: this
    class(scalar_field), allocatable      :: heat
      !! The value of the contribution made by melting/thermal
      !! transfer to the heat equation for a [[plume]]
    if (.not. allocated(this%melt_values)) error stop ('Melt values not allocated')
    call this%melt_values%allocate_scalar_field(heat)
    heat = this%coef * this%melt_values
    call heat%set_temp()
  end function dallaston2015_heat

  pure function dallaston2015_salt(this) result(salt)
    class(dallaston2015_melt), intent(in) :: this
    class(scalar_field), allocatable      :: salt
      !! The value of the contribution made by melting/thermal
      !! transfer to the salt equation for a [[plume]]
  end function dallaston2015_salt

  pure function dallaston2015_melt_rate(this) result(melt)
    class(dallaston2015_melt), intent(in) :: this
    class(scalar_field), allocatable      :: melt
      !! The melt rate from the ice into the plume water.
    if (.not. allocated(this%melt_values)) error stop ('Melt values not allocated')
    call this%melt_values%allocate_scalar_field(melt)
    melt = this%melt_conversion * this%melt_values
    call melt%set_temp()
  end function dallaston2015_melt_rate

  pure function dallaston2015_has_heat(this) result(has_heat)
    class(dallaston2015_melt), intent(in) :: this
    logical                               :: has_heat
      !! Whether this formulation of melting contributes terms to
      !! the heat equation of the plume.
    has_heat = .true.
  end function dallaston2015_has_heat

  pure function dallaston2015_has_salt(this) result(has_salt)
    class(dallaston2015_melt), intent(in) :: this
    logical                               :: has_salt
      !! Whether this formulation of melting contributes terms to
      !! the salinity equation of the plume.
    has_salt = .false.
  end function dallaston2015_has_salt

end module dallaston2015_melt_mod
