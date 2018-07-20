!
!  pseudospectral_block.f90
!  This file is part of ISOFT.
!  
!  Copyright 2017 Chris MacMackin <cmacmackin@gmail.com>
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

module pseudospectral_block_mod
  !* Author: Christopher MacMackin
  !  Date: march 2017
  !  License: GPLv3
  !
  ! Provides a derived type which representes a finite difference
  ! matrix/operator. This can be useful for preconditioning problems
  ! which use a spectral discretisation.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod!, only: abstract_field, scalar_field, vector_field
  use chebyshev_mod, only: collocation_points, integrate_1d
  use penf, only: str
  use logger_mod, only: logger => master_logger
  implicit none
  private

  integer, parameter :: no_extra_derivative = -1
  
  type, public :: pseudospec_block
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! A data type representing a matrix pseudospectral differentiation
    ! operator. It can be useful when preconditioning systems which
    ! use a spectral discretisation, if higher accuracy than finite
    ! difference is needed. It is inherently 1-D in its
    ! implementation. Note that multiplication of a field will simply
    ! call that field's differentiation operator, which may or may not
    ! use a pseudospectral method.
    !
    private
    real(r8), dimension(:), pointer :: xvals
      !! Coordinates of collocation points.
  contains
    private
    procedure :: pseudospec_block_solve_scalar
    procedure :: pseudospec_block_solve_vector
    generic, public :: solve_for => pseudospec_block_solve_scalar, &
                                    pseudospec_block_solve_vector
  end type pseudospec_block

  interface pseudospec_block
    module procedure constructor
  end interface pseudospec_block

contains

  function constructor(template) result(this)
    !* Author: Chris MacMackin
    !  Date: September 2017
    !
    ! Builds a Chebyshsev pseudospectral differentiation matrix block
    ! which can be used to solve the inverse problem. The result can
    ! only be used with fields having the same grid as the template.
    !
    class(abstract_field), intent(in)     :: template
      !! A scalar field with the same grid as any fields passed as
      !! arguments to the [[pseudospec_block(type):solve_for]] method.
    type(pseudospec_block)                :: this

    real(r8), dimension(:,:), allocatable :: domain
    
    domain = template%domain()
    this%xvals => collocation_points(template%elements()-1, domain(1,1), &
                                     domain(1,2))
  end function constructor

  function pseudospec_block_solve_scalar(this, rhs, bound_loc, bound_val, &
                                         good_bound) result(solution)
    !* Author: Chris MacMackin
    !  Date: September 2017
    !
    ! Solves the linear(ised) system represented by this finite
    ! difference block, for a given right hand side state vector
    ! (represented by a scalar field).
    !
    ! @Warning Currently this is only implemented for a 1-D field.
    !
    class(pseudospec_block), intent(inout) :: this
    class(scalar_field), intent(in)        :: rhs
      !! The right hand side of the linear(ised) system.
    integer, intent(in)                    :: bound_loc
      !! Which boundary is being set. The boundary will be the one
      !! normal to dimension of number `abs(boundary)`. If the
      !! argument is negative, then the lower boundary is returned. If
      !! positive, then the upper boundary is returned.
    class(scalar_field), intent(in)        :: bound_val
      !! The value of the result at the specified boundary.
    integer, intent(in), optional          :: good_bound
      !! If provided, indicates which boundary contains trusted
      !! information from which to calculate the power of the highest
      !! frequency mode. Defaults to the opposite of `bound_loc`.
    class(scalar_field), pointer           :: solution

    real(r8), dimension(:), allocatable :: sol_vector, bound_vec
    logical                             :: valid_bound
    integer                             :: bloc
 
    call rhs%guard_temp(); call bound_val%guard_temp()
    valid_bound = .true.
    select case(bound_loc)
    case(1)
      bloc = 1
    case(-1)
      bloc = size(this%xvals)
    case default
      valid_bound = .false.
      !call logger%warning('pseudospec_block%solve_for', '1-D field does not '// &
      !                    'have boundary in dimension '//str(bound_loc))
    end select
    sol_vector = rhs%raw()
    if (valid_bound) then
      bound_vec = bound_val%raw()
      call integrate_1d(sol_vector, this%xvals, bloc, bound_vec(1), good_bound)
    else
      call integrate_1d(sol_vector, this%xvals, good_bound=good_bound)
    end if
#ifdef DEBUG
    call logger%debug('pseudospec_block%solve_for', &
                      'Successfully performed Chebyshev pseudospectral integration.')
#endif
    call rhs%allocate_scalar_field(solution)
    call solution%unset_temp()
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
    call rhs%clean_temp(); call bound_val%clean_temp()
    call solution%set_temp()
  end function pseudospec_block_solve_scalar

  function pseudospec_block_solve_vector(this, rhs, bound_loc, bound_val, &
                                         good_bound) result(solution)
    !* Author: Chris MacMackin
    !  Date: September 2017
    !
    ! Solves the linear(ised) system represented by this finite
    ! difference block, for a given right hand side state vector
    ! (represented by a vector field).
    !
    ! @Warning Currently this is only implemented for a 1-D field.
    !
    ! @Bug For some reason, calls to the `vector_dimensions()` method
    ! produce a segfault when `rhs` is
    ! `class(vector_field)`. Everything works fine if it is
    ! `class(cheb1d_vector_field)`, so this is used as a workaround.
    !
    class(pseudospec_block), intent(inout) :: this
    class(cheb1d_vector_field), intent(in) :: rhs
      !! The right hand side of the linear(ised) system.
    integer, intent(in)                    :: bound_loc
      !! Which boundary is being set. The boundary will be the one
      !! normal to dimension of number `abs(boundary)`. If the
      !! argument is negative, then the lower boundary is returned. If
      !! positive, then the upper boundary is returned.
    class(vector_field), intent(in)        :: bound_val
      !! The value of the result at the specified boundary.
    integer, intent(in), optional          :: good_bound
      !! If provided, indicates which boundary contains trusted
      !! information from which to calculate the power of the highest
      !! frequency mode. Defaults to the opposite of `bound_loc`.
    class(vector_field), pointer           :: solution

    real(r8), dimension(:), allocatable :: sol_vector, bound_vec
    logical                             :: valid_bound
    integer                             :: bloc, n, i
    class(scalar_field), pointer        :: bcomponent
 
    call rhs%guard_temp(); call bound_val%guard_temp()
    valid_bound = .true.
    n = size(this%xvals)
    select case(bound_loc)
    case(1)
      bloc = 1
    case(-1)
      bloc = n
    case default
      valid_bound = .false.
      !call logger%warning('pseudospec_block%solve_for', '1-D field does not '// &
      !                    'have boundary in dimension '//str(bound_loc))
    end select
    sol_vector = rhs%raw()
    call bound_val%allocate_scalar_field(bcomponent)
    call bcomponent%guard_temp()
    do i = 1, rhs%vector_dimensions()
      bcomponent = bound_val%component(i)
      if (valid_bound) then
        bound_vec = bcomponent%raw()
        call integrate_1d(sol_vector((i-1)*n+1:i*n), this%xvals, &
                          bloc, bound_vec(1), good_bound)
      else
        call integrate_1d(sol_vector((i-1)*n+1:i*n), this%xvals, &
                          good_bound=good_bound)
      end if
    end do
    call bcomponent%clean_temp()
#ifdef DEBUG
    call logger%debug('pseudospec_block%solve_for', &
                      'Successfully performed Chebyshev pseudospectral integration.')
#endif
    call rhs%allocate_vector_field(solution)
    call solution%unset_temp()
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
    call rhs%clean_temp(); call bound_val%clean_temp()
    call solution%set_temp()
  end function pseudospec_block_solve_vector

end module pseudospectral_block_mod
