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
    ! density by the haline contraction coefficient \(\beta_S\). The
    ! only real use for this is testing and debugging the plume model.
    !
    private
    class(scalar_field), allocatable :: density
      !! The density calculated from the prescribed salinity
    real(r8)                         :: beta_s
      !! The haline contraction coefficient
  contains
    procedure :: water_density => prescribed_water_density
    procedure, pass(rhs) :: prescribed_assign
    generic :: assignment(=) => prescribed_assign
    procedure :: water_density_derivative => prescribed_water_deriv
    procedure :: haline_contraction => prescribed_haline_contraction
    procedure :: thermal_contraction => prescribed_thermal_contraction
  end type prescribed_eos

  interface prescribed_eos
    module procedure constructor
  end interface prescribed_eos

contains

  function constructor(const, beta_s, thickness) result(this)
    real(r8), intent(in)            :: const
      !! The constant to which \(SD\) is equal.
    real(r8), intent(in)            :: beta_s
      !! The haline contraction coefficient, \(\beta_S\), relating
      !! salinity and density.
    class(scalar_field), intent(in) :: thickness
      !! The thickness of the plume, from which the salinity is calculated.
    type(prescribed_eos)            :: this
    call thickness%guard_temp()
    allocate(this%density, mold=thickness)
    this%density = const*beta_s/thickness
    this%beta_s = beta_s
    call thickness%clean_temp()
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
    class(scalar_field), pointer  :: density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    if (temperature == uniform_scalar_field(0._r8) .and. &
        salinity == uniform_scalar_field(0._r8)) then
      ! Kludge to ensure correct ambient density is returned
      call temperature%allocate_scalar_field(density)
      density = uniform_scalar_field(0._r8)
    else
      call this%density%allocate_scalar_field(density)
      density = this%density
    end if
    call temperature%clean_temp(); call salinity%clean_temp()
    call density%set_temp()
  end function prescribed_water_density

  function prescribed_water_deriv(this, temperature, d_temperature, salinity, &
                             d_salinity, dir) result(d_density)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Calculates the derivative of the water density.
    class(prescribed_eos), intent(in) :: this
    class(scalar_field), intent(in)   :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)   :: d_temperature
      !! A field containing the derivative of the temperature of the
      !! water, in teh same direction as `dir`
    class(scalar_field), intent(in)   :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), intent(in)   :: d_salinity
      !! A field containing the derivative of the salinity of the
      !! water, in the same direction as `dir`
    integer, intent(in)               :: dir
      !! The direction in which to take the derivative
    class(scalar_field), pointer      :: d_density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    call d_temperature%guard_temp(); call d_salinity%guard_temp()
    if (temperature == uniform_scalar_field(0._r8) .and. &
        salinity == uniform_scalar_field(0._r8)) then
      ! Kludge to ensure correct ambient density is returned
      call temperature%allocate_scalar_field(d_density)
      d_density = uniform_scalar_field(0._r8)
    else
      call this%density%allocate_scalar_field(d_density)
      d_density = this%density%d_dx(1)
    end if
    call temperature%clean_temp(); call salinity%clean_temp()
    call d_temperature%clean_temp(); call d_salinity%clean_temp()
    call d_density%set_temp()
  end function prescribed_water_deriv

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


  function prescribed_haline_contraction(this, temperature, salinity) result(coef)
    !* Author: Chris MacMackin
    !  Date: June 2017
    !
    ! Returns the haline contraction coefficient.
    !
    class(prescribed_eos), intent(in)    :: this
    class(scalar_field), intent(in)  :: temperature
    class(scalar_field), intent(in)  :: salinity
    class(scalar_field), allocatable :: coef
    allocate(uniform_scalar_field :: coef)
    coef = uniform_scalar_field(this%beta_s)
  end function prescribed_haline_contraction

  function prescribed_thermal_contraction(this, temperature, salinity) result(coef)
    !* Author: Chris MacMackin
    !  Date: June 2017
    !
    ! Returns the thermal contraction coefficient.
    !
    class(prescribed_eos), intent(in)    :: this
    class(scalar_field), intent(in)  :: temperature
    class(scalar_field), intent(in)  :: salinity
    class(scalar_field), allocatable :: coef
    allocate(uniform_scalar_field :: coef)
    coef = uniform_scalar_field(0.0_r8)
  end function prescribed_thermal_contraction

end module prescribed_eos_mod
