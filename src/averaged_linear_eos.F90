!
!  averaged_linear_eos.f90
!  This file is part of ISOFT.
!  
!  Copyright 2019 Chris MacMackin <cmacmackin@gmail.com>
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

module ave_linear_eos_mod
  !* Author: Christopher MacMackin
  !  Date: August 2018
  !  License: GPLv3
  !
  ! Provides an abstract derived type which can be subtyped in order to
  ! implement an equation of state.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, uniform_scalar_field
  use equation_of_state_mod, only: equation_of_state
  implicit none
  private

  real(r8), parameter :: absolute_zero = -273.15_r8
  
  type, extends(equation_of_state), public :: ave_linear_eos
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! A linearised implementation of the equation of state which has
    ! been horizontally-integrated. The basic equation of stateis $$
    ! \rho(x,y) = \rho_0[1-\beta_T(T(x,y)-T_0) +
    ! \beta_S(S(x,y)-S_0)]. $$
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
    real(r8) :: a_DS = 1.0_r8
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\alpha_{DS} = \frac{1}{y_2 - y_1}
      !! \int^{y_2}_{y_1} f_{D}f_S dy, $$ where \(f_{D}(y)\) and
      !! \(f_S(y)\) are the shapes of the variables \(D\) and \(S\) in
      !! the transverse direction.
    real(r8) :: a_DT = 1.0_r8
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\alpha_{DT} = \frac{1}{y_2 - y_1}
      !! \int^{y_2}_{y_1} f_{D}f_T dy, $$ where \(f_{D}(y)\) and
      !! \(f_T(y)\) are the shapes of the variables \(D\) and \(T\) in
      !! the transverse direction.
    real(r8) :: a_DS_t = 1.0_r8
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\tilde{\alpha}_{DS} =
      !! \frac{1}{\alpha_{D^2}(y_2 - y_1)} \int^{y_2}_{y_1} f^2_{D}f_S
      !! dy, $$ where \(f_{D}(y)\) and \(f_S(y)\) are the shapes of
      !! the variables \(D\) and \(S\) in the transverse direction and
      !! $$\alpha_{D^2} = \frac{1}{y_2 -
      !! y_1}\int^{y_2}_{y_1}f_D^2dy.$$
    real(r8) :: a_DT_t = 1.0_r8
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\tilde{\alpha}_{DT} =
      !! \frac{1}{\alpha_{D^2}(y_2 - y_1)} \int^{y_2}_{y_1} f^2_{D}f_T
      !! dy, $$ where \(f_{D}(y)\) and \(f_T(y)\) are the shapes of
      !! the variables \(D\) and \(T\) in the transverse direction and
      !! $$\alpha_{D^2} = \frac{1}{y_2 -
      !! y_1}\int^{y_2}_{y_1}f_D^2dy.$$
  contains
    procedure :: water_density => linear_water_density
    procedure :: water_density_ave1 => linear_water_density_ave1
    procedure :: water_density_ave2 => linear_water_density_ave2
    procedure :: water_density_derivative => linear_water_deriv
    procedure :: haline_contraction => linear_haline_contraction
    procedure :: thermal_contraction => linear_thermal_contraction
  end type ave_linear_eos

  interface ave_linear_eos
    module procedure constructor
  end interface ave_linear_eos

contains

  pure function constructor(ref_rho, ref_t, ref_s, beta_t, beta_s, a_DS, &
                            a_DT, a_DS_t, a_DT_t) result(this)
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
    real(r8), intent(in), optional :: a_DS
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\alpha_{DS} = \frac{1}{y_2 - y_1}
      !! \int^{y_2}_{y_1} f_{D}f_S dy, $$ where \(f_{D}(y)\) and
      !! \(f_S(y)\) are the shapes of the variables \(D\) and \(S\) in
      !! the transverse direction. Defualt value is 1.
    real(r8), intent(in), optional :: a_DT
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\alpha_{DT} = \frac{1}{y_2 - y_1}
      !! \int^{y_2}_{y_1} f_{D}f_T dy, $$ where \(f_{D}(y)\) and
      !! \(f_T(y)\) are the shapes of the variables \(D\) and \(T\) in
      !! the transverse direction. Defualt value is 1.
    real(r8), intent(in), optional :: a_DS_t
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\tilde{\alpha}_{DS} =
      !! \frac{1}{\alpha_{D^2}(y_2 - y_1)} \int^{y_2}_{y_1} f^2_{D}f_S
      !! dy, $$ where \(f_{D}(y)\) and \(f_S(y)\) are the shapes of
      !! the variables \(D\) and \(S\) in the transverse direction and
      !! $$\alpha_{D^2} = \frac{1}{y_2 -
      !! y_1}\int^{y_2}_{y_1}f_D^2dy.$$ Defualt value is 1.
    real(r8), intent(in), optional :: a_DT_t
      !! The shape coefficient for a horizontally-integrated model. It
      !! is defined as $$\tilde{\alpha}_{DT} =
      !! \frac{1}{\alpha_{D^2}(y_2 - y_1)} \int^{y_2}_{y_1} f^2_{D}f_T
      !! dy, $$ where \(f_{D}(y)\) and \(f_T(y)\) are the shapes of
      !! the variables \(D\) and \(T\) in the transverse direction and
      !! $$\alpha_{D^2} = \frac{1}{y_2 -
      !! y_1}\int^{y_2}_{y_1}f_D^2dy.$$ Defualt value is 1.
    type(ave_linear_eos)     :: this
    this%ref_rho = ref_rho
    this%ref_t   = ref_t
    this%ref_s   = ref_s
    this%beta_t  = beta_t
    this%beta_s  = beta_s
    if (present(a_DS)) this%a_DS = a_DS
    if (present(a_DT)) this%a_DT = a_DT
    if (present(a_DS_t)) this%a_DS_t = a_DS_t
    if (present(a_DT_t)) this%a_DT_t = a_DT_t
  end function constructor

  function linear_water_density(this, temperature, salinity) result(density)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Calculates the density of the water from the temperature and
    ! salinity, using a linear equation of state, $$ \rho(x,y) =
    ! \rho_0[1-\beta_T(T(x,y)-T_0) + \beta_S(S(x,y)-S_0)]. $$
    class(ave_linear_eos), intent(in):: this
    class(scalar_field), intent(in)  :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)  :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), pointer     :: density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    call salinity%allocate_scalar_field(density)
    density = this%ref_rho * (1.0_r8 - this%beta_t*(temperature - this%ref_t) &
                                     + this%beta_s*(salinity - this%ref_s))
    call temperature%clean_temp(); call salinity%clean_temp()
    call density%set_temp()
  end function linear_water_density

  function linear_water_density_ave1(this, temperature, salinity) result(density)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Calculates one form of the horizontally-averaged density of the
    ! water from the temperature and salinity, using a linear equation
    ! of state, $$ \bar{rho}(x) =
    ! \rho_0[1-\beta_T(\alpha_{DT}T(x)-T_0) +
    ! \beta_S(\alpha_{DS}S(x)-S_0)]. $$
    class(ave_linear_eos), intent(in):: this
    class(scalar_field), intent(in)  :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)  :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), pointer     :: density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    call salinity%allocate_scalar_field(density)
    density = this%ref_rho * (1.0_r8 - this%beta_t*(this%a_DT*temperature - this%ref_t) &
                                     + this%beta_s*(this%a_DS*salinity - this%ref_s))
    call temperature%clean_temp(); call salinity%clean_temp()
    call density%set_temp()
  end function linear_water_density_ave1

  function linear_water_density_ave2(this, temperature, salinity) result(density)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Calculates another form of the horizontally-averaged density of
    ! the water from the temperature and salinity, using a linear
    ! equation of state, $$ \tilde{\rho}(x) =
    ! \rho_0[1-\beta_T(\tilde{\alpha}_{DT}T(x)-T_0) +
    ! \beta_S(\tilde{\alpha}_{DS}S(x)-S_0)]. $$
    class(ave_linear_eos), intent(in):: this
    class(scalar_field), intent(in)  :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)  :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), pointer     :: density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    call salinity%allocate_scalar_field(density)
    density = this%ref_rho * (1.0_r8 - this%beta_t*(temperature - this%ref_t) &
                                     + this%beta_s*(salinity - this%ref_s))
    call temperature%clean_temp(); call salinity%clean_temp()
    call density%set_temp()
  end function linear_water_density_ave2

  function linear_water_deriv(this, temperature, d_temperature, salinity, &
                             d_salinity, dir) result(d_density)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Calculates the derivative of the average water density from the
    ! temperature and salinity, using a linear equation of state with
    ! the second type of averaging, $$ \tilde{\rho} =
    ! \rho_0[1-\beta_T(\tilde{\alpha}_{DT}T-T_0) +
    ! \beta_S(\tilde{\alpha}_{DS}S-S_0)]. $$
    class(ave_linear_eos), intent(in):: this
    class(scalar_field), intent(in)  :: temperature
      !! A field containing the temperature of the water
    class(scalar_field), intent(in)  :: d_temperature
      !! A field containing the derivative of the temperature of the
      !! water, in teh same direction as `dir`
    class(scalar_field), intent(in)  :: salinity
      !! A field containing the salinity of the water
    class(scalar_field), intent(in)  :: d_salinity
      !! A field containing the derivative of the salinity of the
      !! water, in the same direction as `dir`
    integer, intent(in)              :: dir
      !! The direction in which to take the derivative
    class(scalar_field), pointer     :: d_density
      !! A field containing the density of the water
    call temperature%guard_temp(); call salinity%guard_temp()
    call d_temperature%guard_temp(); call d_salinity%guard_temp()
    call salinity%allocate_scalar_field(d_density)
    d_density = this%ref_rho*(this%a_DS_t*this%beta_s*d_salinity - &
                              this%a_DT_t*this%beta_t*d_temperature)
    call temperature%clean_temp(); call salinity%clean_temp()
    call d_temperature%clean_temp(); call d_salinity%clean_temp()
    call d_density%set_temp()
  end function linear_water_deriv

  function linear_haline_contraction(this, temperature, salinity) result(coef)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Returns the haline contraction coefficient.
    !
    class(ave_linear_eos), intent(in)    :: this
    class(scalar_field), intent(in)  :: temperature
    class(scalar_field), intent(in)  :: salinity
    class(scalar_field), allocatable :: coef
    allocate(uniform_scalar_field :: coef)
    coef = uniform_scalar_field(this%ref_rho*this%beta_s)
  end function linear_haline_contraction

  function linear_thermal_contraction(this, temperature, salinity) result(coef)
    !* Author: Chris MacMackin
    !  Date: August 2018
    !
    ! Returns the thermal contraction coefficient.
    !
    class(ave_linear_eos), intent(in)    :: this
    class(scalar_field), intent(in)  :: temperature
    class(scalar_field), intent(in)  :: salinity
    class(scalar_field), allocatable :: coef
    allocate(uniform_scalar_field :: coef)
    coef = uniform_scalar_field(this%ref_rho*this%beta_t)
  end function linear_thermal_contraction

end module ave_linear_eos_mod
