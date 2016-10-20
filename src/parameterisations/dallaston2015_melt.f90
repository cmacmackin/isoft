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
    real(r8) :: beta
      !! The inverse stefan number, $$\beta = \frac{c(T_a - T_m}{L}$$
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

  function constructor(beta) result(this)
    real(r8), intent(in) :: beta
      !! The inverse stefan number, $$\beta = \frac{c(T_a - T_m}{L}$$
    type(dallaston2015_melt) :: this
      !! The newly created object representing the melt relationship.
    this%beta = beta
  end function constructor

  subroutine dallaston2015_solve(this, velocity, pressure, temperature, &
                                 salinity, plume_thickness, time)
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
    if (allocated(this%melt_values)) deallocate(this%melt_values)
    allocate(this%melt_values, source=velocity%norm())
  end subroutine dallaston2015_solve

  pure function dallaston2015_heat(this) result(heat)
    class(dallaston2015_melt), intent(in) :: this
    class(scalar_field), allocatable      :: heat
      !! The value of the contribution made by melting/thermal
      !! transfer to the heat equation for a [[plume]]
    if (.not. allocated(this%melt_values)) error stop('Melt values not allocated')
    call this%melt_values%allocate_scalar_field(heat)
    heat = (this%beta + 1.0_r8) * this%melt_values
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
    if (.not. allocated(this%melt_values)) error stop('Melt values not allocated')
    allocate(melt, source=this%melt_values)
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
