!
!  uniform_gradient_field.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
!  
!  This program is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3 of the License, or
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

module uniform_gradient_field_mod
  !* Author: Chris MacMackin
  !  Date: July 2017
  !  License: GPLv3
  !
  ! Provides an extension of the uniform field type which also appears
  ! to have a uniform gradient. This was written to allow some of the
  ! same code used when solving the plume in a Runge-Kutta solver,
  ! where I pass in uniform fields rather than ones using a Chebyshev
  ! grid.
  !
  use iso_fortran_env, only: r8 => real64
  use utils_mod, only: is_nan
  use factual_mod, only: scalar_field, vector_field, uniform_scalar_field, &
                         uniform_vector_field, get_tol
  implicit none
  private
  
  type, extends(uniform_scalar_field), public :: uniform_gradient_field
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! A type of uniform field which also has a uniform gradient. Of
    ! course, this is impossible in practice, but it can be useful for
    ! tricking certain routines into working properly. Ideally a whole
    ! new derived type would be created which just holds the value and
    ! gradient at a single point, but the emphasis is on getting
    ! something quickly. Note that the gradient is not propagated
    ! across operations--the result of all overloaded operators is
    ! just a normal uniform field with no gradient.
    !
    real(r8), dimension(:), allocatable :: grad
      !! The values of the gradient in each direction.
  contains
    private
    procedure, public :: d_dx => uniform_gradient_d_dx
      !! \(\frac{\partial^n}{\partial x_i^n}({\rm field})\)
    procedure :: gradient => uniform_gradient_gradient
      !! \(\nabla {\rm field}\)
    procedure :: is_equal => uniform_gradient_is_equal
      !! Checks fields are equal within a tolerance
    procedure :: assign_field => uniform_gradient_assign
      !! \({\rm field} = {\rm field}\)
  end type uniform_gradient_field

  interface uniform_gradient_field
    module procedure constructor
 end interface uniform_gradient_field

contains

  function constructor(val, grad) result(this)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Creates a new scalar field with a uniform value across all of
    ! space but a non-zero gradient.
    !
    real(r8), intent(in)               :: val
      !! The value of the field
    real(r8), dimension(:), intent(in) :: grad
      !! An array in which the `i`th element contains the gradient in
      !! the _i_th direction. Directions corresponding to values of
      !! `i` greater than the size of the array are taken to have a
      !! gradient of zero.
    type(uniform_gradient_field)       :: this
      !! A scalar field initated based on teh arguments of this
      !! function.
    this%uniform_scalar_field = uniform_scalar_field(val)
    this%grad = grad
  end function constructor

  function uniform_gradient_d_dx(this, dir, order) result(res)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! \(\frac{\partial^{\rm order}}{\partial x_{\rm dir}^{\rm order}}{\rm field}\)
    !
    class(uniform_gradient_field), intent(in) :: this
    integer, intent(in) :: dir !! Direction in which to differentiate
    integer, optional, intent(in) :: order !! Order of the derivative, default = 1
    class(scalar_field), pointer :: res
    integer :: ord
    call this%guard_temp()
    call this%allocate_scalar_field(res)
    if (present(order)) then
      ord = order
    else
      ord = 1
    end if
    select type(res)
    class is(uniform_scalar_field)
      if (ord == 1) then
        if (dir > 0 .and. dir <= size(this%grad)) then
          res = uniform_scalar_field(this%grad(dir))
        else
          res = uniform_scalar_field(0.0_r8)
        end if
      else
        res = uniform_scalar_field(0.0_r8)
      end if
    class default
      error stop ('Non-uniform_gradient_field type allocated by '//&
                  '`allocate_scalar_field` routine.')
    end select
    call res%set_temp() ! Shouldn't need to call this, but for some
                        ! rason being set as non-temporary when
                        ! assignment subroutine returns.
    call this%clean_temp()
  end function uniform_gradient_d_dx

  function uniform_gradient_gradient(this) result(res)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! \(\nabla{\rm field}\)
    !
    class(uniform_gradient_field), intent(in) :: this
    class(vector_field), pointer :: res !! The result of this operation
    call this%guard_temp()
    call this%allocate_vector_field(res)
    select type(res)
    class is(uniform_vector_field)
      res = uniform_vector_field(this%grad)
    class default
      error stop ('Non-uniform_vector_field type allocated by '//&
                  '`allocate_vector_field` routine.')
    end select
    call this%clean_temp()
  end function uniform_gradient_gradient
  
  impure elemental subroutine uniform_gradient_assign(this,rhs)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! \({\rm field} = {\rm field}\)
    !
    class(uniform_gradient_field), intent(inout) :: this
    class(scalar_field), intent(in) :: rhs
    call rhs%guard_temp()
    select type(rhs)
    class is(uniform_gradient_field)
      this%uniform_scalar_field = rhs%uniform_scalar_field
      if (allocated(rhs%grad)) then
        this%grad = rhs%grad
      else if (allocated(this%grad)) then
        deallocate(this%grad)
      end if
      call this%unset_temp()
    class is(uniform_scalar_field)
      this%uniform_scalar_field = rhs
      if (allocated(this%grad)) deallocate(this%grad)
      call this%unset_temp()
    class default
      error stop ('Assigning incompatible type to uniform_gradient_field')
    end select
    call rhs%clean_temp()
  end subroutine uniform_gradient_assign

  logical function uniform_gradient_is_equal(this,rhs) result(iseq)
    !* Author: Chris MacMackin
    !  Date: July 2017
    !
    ! Evaluates whether two scalar fields are equal within a tolerance,
    ! specified by `set_tol`.
    !
    class(uniform_gradient_field), intent(in) :: this
    class(scalar_field), intent(in) :: rhs
    real(r8) :: normalization
    integer :: i
    call this%guard_temp(); call rhs%guard_temp()
    iseq = .true.
    select type(rhs)
    class is(uniform_gradient_field)
      iseq = (this%uniform_scalar_field == rhs%uniform_scalar_field)
      do i = 1, size(this%grad)
        if (.not. iseq) return
        normalization = abs(this%grad(i))
        if (normalization < get_tol()) normalization = 1.0_r8
        iseq = iseq .and.( ((this%grad(i)-rhs%grad(i))/normalization < &
                             get_tol()) .or. (is_nan(this%grad(i)).and. &
                                              is_nan(rhs%grad(i))) )
      end do
    class default
      iseq = (rhs == this)
    end select
    call this%clean_temp(); call rhs%clean_temp()
  end function uniform_gradient_is_equal


end module uniform_gradient_field_mod
