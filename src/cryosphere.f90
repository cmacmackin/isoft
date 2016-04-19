!
!  cryosphere.f90
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

module cryosphere_mod
  !* Author: Christopher MacMackin
  !  Date: April 2016
  !  License: GPLv3
  !
  ! Provides a data structure representing a system of ice sheets and/or ice
  ! shelves, as well as the ground and/or ocean which they interact with. This
  ! is the fundamental data type of the ISOFT software suite.
  !
  use iso_fortran_env, only: r8 => real64
  use foodie, only: integrand
  use basal_surface_mod, only: basal_surface
  use glacier_mod, only: glacier
  implicit none
  private

  type, extends(integrand), public :: cryosphere
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! A data structure representing glaciers, either as ice shelves or 
    ! (eventually) ice sheets. It will allow coupled systems of glaciers as
    ! well as  different basal couplings with the ocean or ground. This type
    ! is a subclass of the 
    ! [FOODIE](https://github.com/Fortran-FOSS-Programmers/FOODIE)
    ! [integrand](http://fortran-foss-programmers.github.io/FOODIE/type/integrand.html),
    ! allowing it to take advantage of that set of integration libraries for
    ! evolution in time.
    !
    private
    real(r8) :: time !! The time in the simulation
  contains
    procedure :: t => cryosphere_dt
    procedure :: local_error => cryosphere_local_error
    procedure :: integrand_multiply_integrand => cryosphere_m_cryosphere
    procedure :: integrand_multiply_real => cryosphere_m_real
    procedure, pass(rhs) :: real_multiply_integrand => real_m_cryosphere
    procedure :: add => cryosphere_add
    procedure :: sub => cryosphere_sub
    procedure :: assign_integrand => cryosphere_assign
  end type cryosphere

contains
  
  function cryosphere_dt(self,t)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Computes the derivative of the cryosphere with respect to time. As
    ! the cryosphere type is really just a pupateer for components which
    ! contain the actual physics, basically this just requires calling the
    ! equivalent functions for the object's glaciers. Note that, as it is
    ! only the ice which is being treated as being dynamic in time (with
    ! the ground being unchanging and the ocean being quasistatic), only
    ! the glacier components will actually have their derivative computed.
    !
    class(cryosphere), intent(in)  :: self
    real(r8), intent(in), optional :: t
      !! Time at which to evaluate the derivative
    class(integrand), allocatable  :: cryosphere_dt
      !! The time rate of change of the cryosphere. Has dynamic type
      !! [[cryosphere]].
  end function cryosphere_dt

  function cryosphere_local_error(lhs, rhs) result(error)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Calculates a real scalar to represent the absolute difference between
    ! two cryosphere objects. `rhs` must be a [[cryosphere]] object, or a
    ! runtime error will occur.
    !
    class(cryosphere), intent(in) :: lhs
      !! Self
    class(integrand), intent(in)  :: rhs
      !! The cryosphere object which is being compared against.
    real(r8) :: error
      !! The scalar representation of the absolute difference between these
      !! two cryospheres.
  end function cryosphere_local_error

  function cryosphere_m_cryosphere(lhs, rhs) result(product)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Multiplies one cryosphere object by another. That is to say, it 
    ! performs element-wise multiplication of the state vectors 
    ! representing the two arguments. `rhs` must be a [[cryosphere]]
    ! object, or a runtime error will occur.
    !
    class(cryosphere), intent(in) :: lhs
      !! Self
    class(integrand), intent(in)  :: rhs
      !! The cryosphere object being multiplied by.
    class(integrand), allocatable :: product
      !! The product of the two arguments. Has dynamic type [[cryosphere]].
  end function cryosphere_m_cryosphere

  function cryosphere_m_real(lhs, rhs) result(product)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Multiplies one cryosphere object by a scalar. That is to say, it 
    ! performs element-wise multiplication of the state vector 
    ! representing the cryosphere.
    !
    class(cryosphere), intent(in) :: lhs
      !! Self
    real(r8), intent(in)          :: rhs
      !! The scalar being multiplied by.
    class(integrand), allocatable :: product
      !! The product of the two arguments. Has dynamic type [[cryosphere]].
  end function cryosphere_m_real

  function real_m_cryosphere(lhs, rhs) result(product)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Multiplies one cryosphere object by a scalar. That is to say, it 
    ! performs element-wise multiplication of the state vector 
    ! representing the cryosphere.
    !
    real(r8), intent(in)          :: lhs
      !! The scalar being multiplied by.
    class(cryosphere), intent(in) :: rhs
      !! Self
    class(integrand), allocatable :: product
      !! The product of the two arguments. Has dynamic type [[cryosphere]].
  end function real_m_cryosphere

  function cryosphere_add(lhs, rhs) result(sum)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Adds one cryosphere object to another. That is to say, it performs
    ! element-wise addition of the state vectors representing the two
    ! arguments. `rhs` must be a [[cryosphere]] object, or a runtime
    ! error will occur.
    !
    class(cryosphere), intent(in) :: lhs
      !! Self
    class(integrand), intent(in)  :: rhs
      !! The cryosphere object being added.
    class(integrand), allocatable :: sum
      !! The sum of the two arguments. Has dynamic type [[cryosphere]].
  end function cryosphere_add

  function cryosphere_sub(lhs, rhs) result(difference)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Subtracts one cryosphere object from another. That is to say, it 
    ! performs element-wise addition of the state vectors representing 
    ! the two arguments. `rhs` must be a [[cryosphere]] object, or a
    ! runtime error will occur.
    !
    class(cryosphere), intent(in) :: lhs
      !! Self
    class(integrand), intent(in)  :: rhs
      !! The cryosphere object being subtracted.
    class(integrand), allocatable :: difference
      !! The difference of the two arguments. Has dynamic type [[cryosphere]].
  end function cryosphere_sub

  subroutine cryosphere_assign(lhs, rhs)
    !* Author: Christopher MacMackin
    !  Date: April 2016
    !
    ! Assigns the `rhs` cryosphere to this, `lhs`, one. All components
    ! will be the same following the assignment.
    !
    class(cryosphere), intent(inout) :: lhs
      !! Self
    class(integrand), intent(in)     :: rhs
      !! The object to be assigned. Must have dynamic type [[cryosphere]],
      !! or a runtime error will occur.
  end subroutine cryosphere_assign

end module cryosphere_mod
