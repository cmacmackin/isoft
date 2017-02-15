!
!  jenkins1991_entrainment.f90
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

module jenkins1991_entrainment_mod
  !* Author: Christopher MacMackin
  !  Date: October 2016
  !  License: GPLv3
  !
  ! Provides a concrete implementation of [[abstract_entrainment]]
  ! in the form of the parameterisation used by Jenkins (1991).
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod!, only: scalar_field, vector_field, abs
  use entrainment_mod, only: abstract_entrainment
  implicit none
  private

  type, extends(abstract_entrainment), public :: jenkins1991_entrainment
    !* Author: Christopher MacMackin
    !  Date: October 2016
    !
    ! A parameterisation of entrainment (\(e\)) as described by Jenkins
    ! (1991): $$e = E_0 |\vec{U}\sin(\theta) \simeq
    ! E_0|\vec{U}||\nabla b|.$$ Here, \(E_0\) is a coefficient typically
    ! taken to be 0.036 (the default value), \(\vec{U}\) is the velocity
    ! of the plume, \(\theta\) is the angle of slope of the ice shelf
    ! base, and \(b\) is the basal depth of the ice shelf.
    !
    private
    real(r8) :: coefficient = 1.0_r8
      !! The entrainment coefficient $E_0$
  contains
    procedure :: entrainment_rate => jenkins1991_rate
      !! Returns the entrainment rate for ambient water into the plume.
  end type jenkins1991_entrainment

  interface jenkins1991_entrainment
    module procedure constructor
  end interface jenkins1991_entrainment

contains

  pure function constructor(coefficient) result(this)
    real(r8), intent(in) :: coefficient
      !! The entrainment coefficient, $E_0$ to be used
    type(jenkins1991_entrainment) :: this
      !! A new entrainment object
    this%coefficient = coefficient
  end function constructor

  pure function jenkins1991_rate(this, velocity, thickness, depth, time) &
                                                result(entrainment)
    !* Author: Christopher MacMackin
    !  Date: October 2016
    !
    ! $$e = E_0 |\vec{U}\sin(\theta) \simeq E_0|\vec{U}||\nabla b|$$
    ! Here, \(E_0\) is a coefficient typically taken to be 0.036 (the
    ! default value), \(\vec{U}\) is the velocity of the plume, \(\theta\)
    ! is the angle of slope of the ice shelf base, and \(b\) is the
    ! basal depth of the ice shelf.
    !
    ! @Warning
    ! The calculation must be performed as
    ! ```fortran
    ! this%coefficient * depth%d_dx(1) * velocity%norm()
    ! ```
    ! with the variables in a different order than how the equation is
    ! usually formulated. If they are in the correct order then
    ! `gfortran` expects the result to be a [[vector_field(type)]]. It
    ! is not clear whether this is due to a bug in `gfortran` or in
    ! `factual`.
    !
    class(jenkins1991_entrainment), intent(in) :: this
    class(vector_field), intent(in)            :: velocity
      !! The velocity field of the plume into which fluid is being 
      !! entrained.
    class(scalar_field), intent(in)            :: thickness
      !! The thickness of the plume into which fluid is being
      !! entrained
    class(scalar_field), intent(in)            :: depth
      !! The depth of the upper surface of the plume into which
      !! fluid is being entrained
    real(r8), intent(in), optional             :: time
      !! The time at which the entrainment is being calculated. If not
      !! present then assumed to be same as previous value passed.
    class(scalar_field), allocatable           :: entrainment
      !! The value of the entrainment
    class(vector_field), allocatable           :: tmp
    call velocity%guard_temp(); call thickness%guard_temp(); call depth%guard_temp()
    call depth%allocate_scalar_field(entrainment)
    call depth%allocate_vector_field(tmp)
    tmp = .grad. depth
    entrainment = tmp%norm() ! Needed due to ICE when try to put aqll on one line. TODO: Create minimal example and submit bug report.
    entrainment = this%coefficient * entrainment * velocity%norm()
    call velocity%clean_temp(); call thickness%clean_temp(); call depth%clean_temp()
    call entrainment%set_temp()
  end function jenkins1991_rate

end module jenkins1991_entrainment_mod
