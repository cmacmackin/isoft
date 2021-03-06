!
!  fin_diff_block.f90
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

module finite_difference_block_mod
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
  use boundary_types_mod, only: dirichlet, neumann, free_boundary
  use f95_lapack, only: la_gtsvx
  use nitsol_mod, only: dnrm2
  use penf, only: str
  use logger_mod, only: logger => master_logger
  implicit none
  private

  integer, parameter :: no_extra_derivative = -1
  
  type, public :: fin_diff_block
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! A data type representing a matrix finite difference operator. It
    ! can be useful when preconditioning systems which use a spectral
    ! discretisation. It is inherently 1-D in its implementation. Note
    ! that multiplication of a field will simply call that field's
    ! differentiation operator, which may or may not use a finite
    ! difference method.
    !
    ! When constructing an instance of this type, a template field
    ! must be passed which has the same grid as any other fields which
    ! will be operated upon. Additionally, types and locations of
    ! boundary conditions must be passed.
    !
    private
    real(r8), dimension(:), allocatable :: diagonal
      !! The diagonal of the tridiagonal matrix representation of this
      !! block.
    real(r8), dimension(:), allocatable :: super_diagonal
      !! The super-diagonal of the tridiagonal matrix representation
      !! of this block.
    real(r8), dimension(:), allocatable :: sub_diagonal
      !! The sub-diagonal of the tridiagonal matrix representation of
      !! this block.
    real(r8), dimension(:), allocatable :: l_multipliers
      !! Multipliers defining the L matrix in the LU factorisation of
      !! the tridiagonal matrix representation of this block.
    real(r8), dimension(:), allocatable :: u_diagonal
      !! The diagonal of the U matrix in the LU factorisation of
      !! the tridiagonal matrix representation of this block.
    real(r8), dimension(:), allocatable :: u_superdiagonal1
      !! The first superdiagonal of the U matrix in the LU
      !! factorisation of the tridiagonal matrix representation of
      !! this block.
    real(r8), dimension(:), allocatable :: u_superdiagonal2
      !! The second superdiagonal of the U matrix in the LU
      !! factorisation of the tridiagonal matrix representation of
      !! this block.
    integer, dimension(:), allocatable  :: pivots
      !! Pivot indicies from the LU factorisation of the tridiagonal
      !! matrix representation of this block.
    integer, dimension(:), allocatable  :: boundary_locs
      !! Locations in the raw arrays which are used to specify
      !! boundary conditions.
    integer, dimension(:), allocatable  :: boundary_types
      !! The types of boundary conditions, specified using the
      !! parameters found in [[boundary_types_mod]].
    logical                             :: had_offset = .true.
      !! True if the factorisation was computed from a tridiagonal
      !! system in which an offset was added to the diagonal.
!    real(r8)                            :: magnitude
!      !! The norm of the superdiagonal of the matrix. If an iterative
!      !! method is used, this is needed to decide how to do so.
  contains
    private
    procedure :: fin_diff_block_solve_scalar
    procedure :: fin_diff_block_solve_vector
    generic, public :: solve_for => fin_diff_block_solve_scalar, &
                                    fin_diff_block_solve_vector
  end type fin_diff_block

  interface fin_diff_block
    module procedure constructor
  end interface fin_diff_block

contains

  function constructor(template, boundary_locs, boundary_types) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Build a tridiagonal matrix block for finite differences. See the
    ! end of the documentation of the [[fin_diff_block(type)]] type
    ! for a description of how boundary conditions are treated.
    !
    class(abstract_field), intent(in)           :: template
      !! A scalar field with the same grid as any fields passed as
      !! arguments to the [[fin_diff_block(type):solve_for]] method.
    integer, dimension(:), optional, intent(in) :: boundary_locs
      !! The locations in the raw representation of `rhs` for which
      !! boundary conditions are specified. Defaults to there being
      !! none.
    integer, dimension(:), optional, intent(in) :: boundary_types
      !! Integers specifying the type of boundary condition. The type
      !! of boundary condition corresponding to a given integer is
      !! specified in [[boundary_types_mod]]. Only Dirichlet and
      !! Neumann conditions are supported. Defaults to Dirichlet. The
      !! order in which they are stored must match that of
      !! `boundary_locs`.
    type(fin_diff_block)                        :: this
      !! A new finite difference operator

    integer                                   :: i, n, pos
    logical                                   :: use_cached
    real(r8), dimension(:), allocatable, save :: cached_dx_c
    real(r8), save                            :: cached_upper_bound, &
                                                 cached_lower_bound
    class(vector_field), pointer              :: grid
    class(abstract_field), allocatable, save  :: cached_field_type
    real(r8), dimension(:,:), allocatable     :: domain

    if (present(boundary_locs)) then
      this%boundary_locs = boundary_locs
    else
      allocate(this%boundary_locs(0))
    end if
    if (present(boundary_types)) then
      this%boundary_types = boundary_types
    else
      allocate(this%boundary_types(size(this%boundary_locs)))
      this%boundary_types = dirichlet
    end if

    n = template%raw_size()
    domain = template%domain()
    ! Try to use cached copy of inverse grid spacing, if available and suitable 
    if (allocated(cached_field_type)) then
      use_cached = same_type_as(template, cached_field_type) .and. &
                   abs(domain(1,1) - cached_lower_bound) < 1e-10 .and. &
                   abs(domain(1,2) - cached_upper_bound) < 1e-10 .and. &
                   n==size(cached_dx_c)
    else
      use_cached = .false.
    end if
    ! Construct an array containing the distance between grid
    ! points, if necessary
    if (.not. use_cached) then
      if (allocated(cached_field_type)) then
        deallocate(cached_field_type)
        deallocate(cached_dx_c)
      end if
#ifdef DEBUG
      call logger%debug('fin_diff_block','Calculating and caching '// &
                       'grid spacings.')
#endif
      cached_lower_bound = domain(1,1)
      cached_upper_bound = domain(1,2)
      allocate(cached_field_type, mold=template)
      allocate(cached_dx_c(n))
      call template%allocate_vector_field(grid)
      grid = template%grid_spacing()
      ! Use this array to temporarily hold the grid spacings
      cached_dx_c = grid%raw()
      cached_dx_c(2:n-1) = 0.5_r8/cached_dx_c(2:n-1)
      cached_dx_c(1) = 1._r8/cached_dx_c(1)
      cached_dx_c(n) = 1._r8/cached_dx_c(n)
    end if
    ! Create tridiagonal matrix
    allocate(this%diagonal(n))
    this%diagonal = 0._r8
    this%diagonal(1) = -cached_dx_c(1)
    this%diagonal(n) = cached_dx_c(n)
    this%super_diagonal = cached_dx_c(1:n-1)
    this%sub_diagonal = -cached_dx_c(2:n)
    ! Apply boundary conditions
    this%diagonal(1) = this%diagonal(1) + 1e-7
    do i = 1, size(this%boundary_locs)
      pos = this%boundary_locs(i)
      select case(this%boundary_types(i))
      case(dirichlet)
        this%diagonal(pos) = 1._r8
        if (pos < n) this%super_diagonal(pos) = 0._r8
        if (pos > 1) this%sub_diagonal(pos-1) = 1.e-7_r8
      case(neumann)
        if (pos == n) then
          this%diagonal(n) = cached_dx_c(n)
          this%sub_diagonal(n-1) = -cached_dx_c(n)
        else if (pos == 1) then
          this%diagonal(1) = -cached_dx_c(1)
          this%super_diagonal(1) = cached_dx_c(1)
        else
          this%super_diagonal(pos) = cached_dx_c(pos)
          this%sub_diagonal(pos-1) = cached_dx_c(pos)
          this%diagonal(pos) = 0._r8
        end if
      case(free_boundary)
        continue
      case default
        call logger%fatal('fin_diff_block','Boundary condition of '// &
                          'type other than Dirichlet or Neumann encountered.')
        error stop
      end select
    end do
    do i = 1, size(this%boundary_locs)
      pos = this%boundary_locs(i)
      select case(this%boundary_types(i))
      case(dirichlet)
        this%diagonal(pos) = 1._r8
        if (pos < n) this%super_diagonal(pos) = 0._r8
        if (pos > 1) this%sub_diagonal(pos-1) = 1.e-7_r8
      case(neumann)
        if (pos == n) then
          this%diagonal(n) = cached_dx_c(n)
          this%sub_diagonal(n-1) = -cached_dx_c(n)
        else if (pos == 1) then
          this%diagonal(1) = -cached_dx_c(1)
          this%super_diagonal(1) = cached_dx_c(1)
        else
          this%super_diagonal(pos) = cached_dx_c(pos)
          this%sub_diagonal(pos-1) = cached_dx_c(pos)
          this%diagonal(pos) = 0._r8
        end if
      case(free_boundary)
        continue
      case default
        call logger%fatal('fin_diff_block','Boundary condition of '// &
                          'type other than Dirichlet or Neumann encountered.')
        error stop
      end select
    end do
    !this%magnitude = dnrm2(n-1, this%super_diagonal, 1)
#ifdef DEBUG
    call logger%debug('fin_diff_block','Instantiated a finite difference '// &
                      'block object.')
#endif
  end function constructor

!  subroutine fin_diff_block_set_bounds(this, subdiag, diag, supdiag)
!    !* Author: Chris MacMackin
!    !  Date: June 2017
!    !
!    ! Alters the provided tridiagonal vector to incorporate the
!    ! boundary conditions of this block.
!    !
!    class(fin_diff_block), intent(in)     :: this
!    real(r8), dimension(:), intent(inout) :: subdiag
!      !! The subdiagonal
!    real(r8), dimension(:), intent(inout) :: diag
!      !! The diagonal
!    real(r8), dimension(:), intent(inout) :: supdiag
!      !! The diagonal
!    integer :: i, n, pos
!    do i = 1, size(this%boundary_locs)
!      pos = this%boundary_locs(i)
!      select case(this%boundary_types(i))
!      case(dirichlet)
!        diag(pos) = 1._r8
!        if (pos < n) supdiag(pos) = 0._r8
!        if (pos > 1) subdiag(pos-1) = 1.e-7_r8
!      case(neumann)
!        if (pos == n) then
!          diag(n) = cached_dx_c(n)
!          subdiag(n-1) = -cached_dx_c(n)
!        else if (pos == 1) then
!          diag(1) = -cached_dx_c(1)
!          supdiag(1) = cached_dx_c(1)
!        else
!          supdiag(pos) = cached_dx_c(pos)
!          subdiag(pos-1) = cached_dx_c(pos)
!          diag(pos) = 0._r8
!        end if
!      case(free_boundary)
!        continue
!      case default
!        call logger%fatal('fin_diff_block','Boundary condition of '// &
!                          'type other than Dirichlet or Neumann encountered.')
!        error stop
!      end select
!    end do
!  end subroutine fin_diff_block_set_bounds

  function fin_diff_block_solve_scalar(this, rhs, offset) result(solution)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Solves the linear(ised) system represented by this finite
    ! difference block, for a given right hand side state vector
    ! (represented by a scalar field). Optionally, the differential
    ! operator can be augmented by adding an offset, i.e. a scalar
    ! field which is added to the operator.
    !
    ! @Warning Currently this is only implemented for a 1-D field.
    !
    class(fin_diff_block), intent(inout)      :: this
    class(scalar_field), intent(in)           :: rhs
      !! The right hand side of the linear(ised) system.
    class(scalar_field), optional, intent(in) :: offset
      !! An offset to add to the differential operator
    class(scalar_field), pointer              :: solution

    real(r8), dimension(:), allocatable :: sol_vector, diag_vector
    integer                             :: flag, n, i, j
    real(r8)                            :: forward_err, &
                                           backward_err, &
                                           condition_num
    character(len=1)                    :: factor
    character(len=:), allocatable       :: msg
 
    call rhs%guard_temp()
    if (present(offset)) call offset%guard_temp()
    allocate(sol_vector(rhs%raw_size()))
    ! Allocate the arrays used to hold the factorisation of the
    ! tridiagonal matrix
    n = size(this%diagonal)
    if (.not. allocated(this%pivots)) then
      allocate(this%l_multipliers(n-1))
      allocate(this%u_diagonal(n))
      allocate(this%u_superdiagonal1(n-1))
      allocate(this%u_superdiagonal2(n-2))
      allocate(this%pivots(n))
      factor = 'N'
    else
      if (.not. present(offset) .and. .not. this%had_offset) then
        factor = 'F'
      else
        factor = 'N'
      end if
    end if

    if (present(offset)) then
#ifdef DEBUG
      if (offset%elements() /= n) then
        call logger%fatal('fin_diff_block%solve_for','Offset field has '// &
                          'different resolution than finite difference block.')
      end if
#endif
      allocate(diag_vector(n))

      where (this%diagonal == 0._r8)
        diag_vector = this%diagonal + offset%raw()
      elsewhere
        diag_vector = this%diagonal
      end where
      if (.not. any(1 == this%boundary_locs)) then
        diag_vector(1) = diag_vector(1) + offset%get_element(1)
      end if
      if (.not. any(n == this%boundary_locs)) then
        diag_vector(n) = diag_vector(n) + offset%get_element(n)
      end if
      do i=1, size(this%boundary_locs)
        j = this%boundary_locs(i)
        if ((j == 1 .or. j == n) .and. this%boundary_types(i) == free_boundary) then
          diag_vector(j) = diag_vector(j) + offset%get_element(j)
        end if
      end do
      call la_gtsvx(this%sub_diagonal, diag_vector, this%super_diagonal, &
                    rhs%raw(), sol_vector, this%l_multipliers,           &
                    this%u_diagonal, this%u_superdiagonal1,              &
                    this%u_superdiagonal2, this%pivots, factor, 'N',     &
                    forward_err, backward_err, condition_num, flag)
    else
      call la_gtsvx(this%sub_diagonal, this%diagonal, this%super_diagonal,      &
                    rhs%raw(), sol_vector, this%l_multipliers, this%u_diagonal, &
                    this%u_superdiagonal1, this%u_superdiagonal2, this%pivots,  &
                    factor, 'N', forward_err, backward_err, condition_num,      &
                    flag)
    end if
    if (flag/=0) then
      msg = 'Tridiagonal matrix solver returned with flag '//str(flag)
      call logger%error('fin_diff_block%solve_for',msg)
    else
      this%had_offset = present(offset)
#ifdef DEBUG
      msg = 'Tridiagonal matrix solver retunred with estimated condition '// &
            'number '//str(condition_num)
      call logger%debug('fin_diff_block%solve_for',msg)
#endif
    end if
    call rhs%allocate_scalar_field(solution)
    call solution%unset_temp()
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
    if (present(offset)) call offset%clean_temp()
    call rhs%clean_temp(); call solution%set_temp()
  end function fin_diff_block_solve_scalar

  function fin_diff_block_solve_vector(this, rhs, offset) result(solution)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Solves the linear(ised) system represented by this finite
    ! difference block, for a given right hand side state vector
    ! (represented by a vector field). Optionally, the differential
    ! operator can be augmented by adding an offset, i.e. a vector
    ! field which is added to the operator.
    !
    ! @Warning Currently this is only implemented for a 1-D field.
    !
    ! @Bug For some reason, calls to the `vector_dimensions()` method
    ! produce a segfault when `rhs` is
    ! `class(vector_field)`. Everything works fine if it is
    ! `class(cheb1d_vector_field)`, so this is used as a workaround.
    !
    class(fin_diff_block), intent(inout)      :: this
    class(cheb1d_vector_field), intent(in)    :: rhs
      !! The right hand side of the linear(ised) system.
    class(cheb1d_vector_field), optional, intent(in) :: offset
      !! An offset to add to the differential operator
    class(vector_field), pointer              :: solution

    real(r8), dimension(:), allocatable :: sol_vector, diag_vector
    integer                             :: flag, n, i, j, k
    class(scalar_field), pointer        :: component, ocomponent
    integer                             :: m
    real(r8)                            :: forward_err, &
                                           backward_err, &
                                           condition_num
    character(len=1)                    :: factor
    character(len=:), allocatable       :: msg
 
    call rhs%guard_temp()
    if (present(offset)) call offset%guard_temp()
    allocate(sol_vector(rhs%raw_size()))
    n = size(this%diagonal)
    ! Allocate the arrays used to hold the factorisation of the
    ! tridiagonal matrix
    if (.not. allocated(this%pivots)) then
      allocate(this%l_multipliers(n-1))
      allocate(this%u_diagonal(n))
      allocate(this%u_superdiagonal1(n-1))
      allocate(this%u_superdiagonal2(n-2))
      allocate(this%pivots(n))
      factor = 'N'
    else
      if (.not. present(offset) .and. .not. this%had_offset) then
        factor = 'F'
      else
        factor = 'N'
      end if
    end if

    call rhs%allocate_scalar_field(component)
    call component%guard_temp()
    if (present(offset)) then
#ifdef DEBUG
      if (offset%elements() /= n) then
        call logger%fatal('fin_diff_block%solve_for','Offset field has '// &
                          'different resolution than finite difference block.')
        error stop
      else if (offset%vector_dimensions() < rhs%vector_dimensions()) then
        call logger%fatal('fin_diff_block%solve_for','Offset field has '// &
                          'different number of vector components than '// &
                          'field being solved for.')
        error stop
      end if
#endif
      call offset%allocate_scalar_field(ocomponent)
      call ocomponent%guard_temp()
    end if

    do i = 1, rhs%vector_dimensions()
      if (i > 1 .and. .not. present(offset)) factor = 'F'
      component = rhs%component(i)
      if (present(offset)) then
        ocomponent = offset%component(i)
        if (i == 1) allocate(diag_vector(n))
        where (this%diagonal == 0._r8)
          diag_vector = this%diagonal + ocomponent%raw()
        elsewhere
          diag_vector = this%diagonal
        end where
        if (.not. any(1 == this%boundary_locs)) then
          diag_vector(1) = diag_vector(1) + ocomponent%get_element(1)
        end if
        if (.not. any(n == this%boundary_locs)) then
          diag_vector(n) = diag_vector(n) + ocomponent%get_element(n)
        end if
        do j=1, size(this%boundary_locs)
          k = this%boundary_locs(j)
          if ((k == 1 .or. k == n) .and. this%boundary_types(j) == free_boundary) then
            diag_vector(k) = diag_vector(k) + ocomponent%get_element(k)
          end if
        end do
        call la_gtsvx(this%sub_diagonal, diag_vector, this%super_diagonal, &
                      component%raw(), sol_vector((i-1)*n+1:i*n),          &
                      this%l_multipliers, this%u_diagonal,                 &
                      this%u_superdiagonal1, this%u_superdiagonal2,        &
                      this%pivots, factor, 'N', forward_err, backward_err, &
                      condition_num, flag)
      else
        call la_gtsvx(this%sub_diagonal, this%diagonal, this%super_diagonal, &
                      component%raw(), sol_vector((i-1)*n+1:i*n),            &
                      this%l_multipliers, this%u_diagonal,                   &
                      this%u_superdiagonal1, this%u_superdiagonal2,          &
                      this%pivots, factor, 'N', forward_err, backward_err,   &
                      condition_num, flag)
      end if
    end do
    call component%clean_temp()
    if(present(offset)) call ocomponent%clean_temp()
    if (flag/=0) then
      msg = 'Tridiagonal matrix solver returned with flag '//str(flag)
      call logger%error('fin_diff_block%solve_for',msg)
    else
      this%had_offset = present(offset)
#ifdef DEBUG
      msg = 'Tridiagonal matrix solver returned with estimated condition '// &
            'number '//str(condition_num)
      call logger%debug('fin_diff_block%solve_for',msg)
#endif
    end if
    call rhs%allocate_vector_field(solution)
    call solution%unset_temp()
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
    if (present(offset)) call offset%clean_temp()
    call rhs%clean_temp(); call solution%set_temp()
  end function fin_diff_block_solve_vector

end module finite_difference_block_mod
