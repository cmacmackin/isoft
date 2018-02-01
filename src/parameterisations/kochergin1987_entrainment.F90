!
!  kochergin1987_entrainment.f90
!  This file is part of ISOFT.
!  
!  Copyright 2018 Chris MacMackin <cmacmackin@physics.ox.ac.uk>
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

module kochergin1987_entrainment_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of [[abstract_entrainment]]
  ! in the form of the parameterisation described by Kochergin (1987).
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod!, only: scalar_field, vector_field, abs
  use entrainment_mod, only: abstract_entrainment
  implicit none
  private

  type, extends(abstract_entrainment), public :: kochergin1987_entrainment
    !* Author: Christopher MacMackin
    !  Date: Feburary 2018
    !
    ! A parameterisation of entrainment (\(e\)) as described by
    ! Kochergin (1987): $$e =
    ! \frac{c_L^2}{S_m}\sqrt{|\vec{U}|^2+\frac{g'D}{S_m}}.$$ Here,
    ! \(c_L\) is an entrainment coefficient, \(\vec{U}\) is the
    ! velocity of the plume, \(g'\) is the reduced gravity, and
    ! \(S_m\) is the turbulent Schmidt number. The latter-most can be
    ! expressed as $$ S_m = \frac{\Ri}{0.0725(\Ri + 0.186 -
    ! \sqrt{\Ri^2 - 0.316\Ri + 0.0346})} $$, where \(\Ri =
    ! g'D/|\vec{U}|^2\) is the Richardson number.
    !
    private
    real(r8) :: coefficient = 1.0_r8
      !! The entrainment coefficient \(c_L^2x_0/D_0\)
    real(r8) :: delta = 0.036_r8
      !! The ratio \(D_0/h_0\)
  contains
    procedure :: entrainment_rate => kochergin1987_rate
      !! Returns the entrainment rate for ambient water into the plume.
  end type kochergin1987_entrainment

  interface kochergin1987_entrainment
    module procedure constructor
  end interface kochergin1987_entrainment

contains

  pure function constructor(coefficient, delta) result(this)
    real(r8), intent(in) :: coefficient
      !! The entrainment coefficient \(c_L^2x_0/D_0\)
    real(r8), intent(in) :: delta
      !! The ratio \(D_0/h_0\)
    type(kochergin1987_entrainment) :: this
      !! A new entrainment object
    this%coefficient = coefficient
    this%delta = delta
  end function constructor

  function kochergin1987_rate(this, velocity, thickness, depth, density_diff, time) &
                                                result(entrainment)
    !* Author: Christopher MacMackin
    !  Date: Feburary 2018
    !
    ! $$e = \frac{c_L^2}{S_m}\sqrt{|\vec{U}|^2+\frac{g'D}{S_m}}.$$
    ! Here, \(c_L\) is an entrainment coefficient, \(\vec{U}\) is the
    ! velocity of the plume, \(g'\) is the reduced gravity, and
    ! \(S_m\) is the turbulent Schmidt number. The Schmidt number is a
    ! function of the Richardson number \(\Ri = g'D/|\vec{U}|^2\):
    ! $$ S_m = \frac{\Ri}{0.0725(\Ri + 0.186 -
    ! \sqrt{\Ri^2 - 0.316\Ri + 0.0346})}. $$
    !
    class(kochergin1987_entrainment), intent(in) :: this
    class(vector_field), intent(in)            :: velocity
      !! The velocity field of the plume into which fluid is being 
      !! entrained.
    class(scalar_field), intent(in)            :: thickness
      !! The thickness of the plume into which fluid is being
      !! entrained
    class(scalar_field), intent(in)            :: depth
      !! The depth of the upper surface of the plume into which
      !! fluid is being entrained
    class(scalar_field), intent(in)            :: density_diff
      !! The difference between the ambient density and the density of
      !! the plume into which the ambient fluid is being entrained.
    real(r8), intent(in), optional             :: time
      !! The time at which the entrainment is being calculated. If not
      !! present then assumed to be same as previous value passed.
    class(scalar_field), pointer               :: entrainment
      !! The value of the entrainment
    class(scalar_field), pointer               :: Ri, Sm
    call velocity%guard_temp(); call thickness%guard_temp() 
    call depth%guard_temp(); call density_diff%guard_temp()
    call depth%allocate_scalar_field(entrainment)
    call depth%allocate_scalar_field(Ri)
    call depth%allocate_scalar_field(Sm)
    call Ri%guard_temp(); call Sm%guard_temp()
    entrainment = velocity%norm() ! Have to awkwardly split this operation to prevent ICE
    Ri = this%delta*density_diff*thickness/(entrainment**2)
    Sm = Ri/(0.0725_r8*(Ri + 0.186_r8 - sqrt(Ri**2 - 0.316_r8*Ri + 0.0346_r8)))
    entrainment = this%coefficient*entrainment/Sm * sqrt(1._r8 + Ri/Sm)
    call velocity%clean_temp(); call thickness%clean_temp()
    call depth%clean_temp(); call density_diff%clean_temp()
    call Ri%clean_temp(); call Sm%clean_temp()
    call entrainment%set_temp()
  end function kochergin1987_rate

end module kochergin1987_entrainment_mod
