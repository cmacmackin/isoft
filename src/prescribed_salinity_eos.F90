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

#ifdef DEBUG
#define pure 
#define elemental 
#endif

module prescribed_eos_mod
  !* Author: Christopher MacMackin
  !  Date: March 2017
  !  License: GPLv3
  !
  ! Provides an equation of state where the salinity is prescribed
  ! such that \(SD = {\rm constant}\). This is useful for testing and
  ! debugging the plume model.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, uniform_scalar_field
  use equation_of_state_mod, only: equation_of_state
  implicit none
  private

  real(r8), parameter :: absolute_zero = -273.15_r8
  
  type, extends(equation_of_state), public :: prescribed_eos
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! An equation of state, depending only on salinity, where the
    ! salinity is prescribed such that \(SD = {\rm constant}\) for
    ! some specified thickness \(D\). The salinity is related to the
    ! density by the haline contraction coefficient \(\beta_S\).
    !
    private
    class(scalar_field), allocatable :: density
      !! The density calculated from the prescribed salinity
  contains
    procedure :: water_density => prescribed_water_density
    procedure, pass(rhs) :: prescribed_assign
    generic :: assignment(=) => prescribed_assign
  end type prescribed_eos

  interface prescribed_eos
    module procedure constructor
  end interface prescribed_eos

contains

  pure function constructor(const, beta_s, thickness) result(this)
    real(r8), intent(in)            :: const
      !! The constant to which \(SD\) is equal.
    real(r8), intent(in)            :: beta_s
      !! The haline contraction coefficient, \(\beta_S\), relating
      !! salinity and density.
    class(scalar_field), intent(in) :: thickness
      !! The thickness of the plume, from which the salinity is calculated.
    type(prescribed_eos)            :: this
    call thickness%allocate_scalar_field(this%density)
    this%density = const*beta_s/thickness
  end function constructor

  function prescribed_water_density(this, temperature, salinity) result(density)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Returns the density corresponding to the prescribed salinity, as
    ! calculated in the constructor.
    !
    class(prescribed_eos), intent(in) :: this
    class(scalar_field), intent(in)   :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)   :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), allocatable  :: density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    if (temperature == uniform_scalar_field(0._r8) .and. &
        salinity == uniform_scalar_field(0._r8)) then
      allocate(uniform_scalar_field :: density)
      density = uniform_scalar_field(1._r8)
    else
      allocate(density, mold=this%density)
      density = this%density
    end if
    call temperature%clean_temp(); call salinity%clean_temp()
    call density%set_temp()
  end function prescribed_water_density

  subroutine prescribed_assign(lhs, rhs)
    !* Author: Chris MacMackin
    !  Date: March 2017
    !
    ! Assigns this object to another equation of state object,
    ! allowing the definided assignment for the precalculated density
    ! field to work correctly.
    !
    class(equation_of_state), intent(out) :: lhs
    class(prescribed_eos), intent(in)     :: rhs
    select type(lhs)
    class is(prescribed_eos)
      if (allocated(lhs%density)) then
        if (.not. same_type_as(lhs%density, rhs%density)) then
          deallocate(lhs%density)
          allocate(lhs%density, mold=rhs%density)
        end if
      else
        allocate(lhs%density, mold=rhs%density)
      end if
      lhs%density = rhs%density
    class default
      error stop ("Can't assign to `equation_of_state` object of class other "// &
                  "than `prescribed_eos`.")
    end select
  end subroutine prescribed_assign

end module prescribed_eos_mod
