!
!  jacobian_block.f90
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

module jacobian_block_mod
  !* Author: Christopher MacMackin
  !  Date: December 2016
  !  License: GPLv3
  !
  ! Provides a derived type which is useful for representing blocks in
  ! a Jacobian matrix. This follows the abstract calculus pattern,
  ! providing operators for matrix multiplication and for solving the
  ! linear(ised) system. See the documentation for the
  ! [[jacobian_block(type)]] type for more details.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field
  use f95_lapack, only: la_gtsvx
  use logger_mod, only: logger => master_logger
  implicit none
  private

  integer, parameter :: no_extra_derivative = -1
  
  type, public :: jacobian_block
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! A data type representing a submatrix of a
    ! Jacobian. Specifically, it represents the commonly occurring
    ! operation $$ \frac{\partial F}{\partial x_i} + F\Delta_i, $$
    ! where \( \Delta_i \) is the differentiation operator in the
    ! \(i\)-direction. Optionally, there can be an additional
    ! differentiation operator on the right-hand-side of this.
    !
    ! This data type is useful when constructing a preconditioner. A
    ! preconditioner, \(A\), should approximate the inverse Jacobian of
    ! a problem being solved. Rather than directly computing the
    ! inverse Jacobian, it is more efficient to approximate it. If
    ! \(d\) is the vector being preconditioned, and \(z\) is the
    ! result of applying the preconditioner, then $$ z = J^{-1}d
    ! \Rightarrow Jz = d.$$ Thus, the preconditioner can be applied
    ! by approximately solving this system for \(z\). Linearising
    ! \(J\), this system can be solved efficiently using Picard
    ! iteration.
    !
    ! Say that the Jacobian is an \(n\times n\) system of these blocks
    ! (each labeled as \(J_{j,k}\)) and that the vector being
    ! preconditioned constists of \(n\) scalar fields (\(d_j\)). Then
    ! the \(m^{th}\) estimate of the solution for the \(j^{th}\) field
    ! in the preconditioned vector (\(z^m_j\)) is the solution to $$
    ! J_{j,j}z^m_j = d_j - \sum_{\substack{k=1\ k\ne j}}^n
    ! J_{j,k}z^{m-1}_k. $$ Depending on the type of fields being used
    ! and the direction in which derivatives are being taken,
    ! \(J_{j,j}\) may be tridiaganol, meaning it can be solved
    ! efficiently. <!--Otherwise, it can be approximated that $$
    ! J_{j,j}z^m_j = \left(\frac{\partial F}{\partial x_i} +
    ! F\Delta_i\right)z^{m} \simeq \frac{\partial F}{\partial
    ! x_i}z^m_j + F\frac{\partial z^{m-1}_j}{\partial x_i}. $$ The
    ! first term in this approximation corresponds to a diagonal
    ! matrix (which can be solved trivially), while the second term is
    ! known and can be subtracted from the right-hand-side of the
    ! linear system.-->
    !
    ! For this purpose, a type-bound multiplication operator is
    ! provided for the Jacobian block type, which can be used when
    ! evaluating the right hand side of the linear system. There is
    ! also a `solve_for` method, which takes as an argument the field
    ! representing the right hand side of the system. <!--Before
    ! applying the latter operator, the `update_estimate` method must
    ! be called. This allows for the term involving derivative of the
    ! current estimate of the solution to be subtracted from the
    ! right-hand-side of the linear system, in cases where the system
    ! can not be expressed as a tridiagonal matrix.-->
    !
    private
    integer                             :: direction
      !! The direction in which any derivatives are taken.
    integer                             :: extra_derivative = no_extra_derivative
      !! The direction in which to apply a differentiation on the
      !! right-hand-side of the Jacobian block operator. Defaults
      !! none.
    class(scalar_field), allocatable    :: contents
      !!The value, \(A\), to which the Jacobian block operation is
      !!beiing applied.
    procedure(jacobian_block_bounds), pointer, nopass  :: &
                                           get_boundaries => jacobian_block_bounds
      !! A subroutine which determines the expected boundary value
      !! (and their location in the raw array) for the solution of the
      !! Jacobian block. It effectively inverts the process by which
      !! residuals are determined for boundary values.
    class(scalar_field), allocatable    :: derivative
      !! The cached derivative of `contents`
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
    real(r8), dimension(:), allocatable :: boundary_vals
      !! Expected boundary values for the solution to the Jacobian
      !! block.
    integer, dimension(:), allocatable  :: boundary_locs
      !! Locations in the raw arrays which are used to specify
      !! boundary conditions.
    real(r8)                            :: scalar_increment
      !! A scalar value which is to be added to this Jacobian block
      !! (i.e. to the diagonal).
    logical                             :: has_increment = .false.
      !! Indicates whether or not there has been a `scalar_increment`
      !! added to this block.
  contains
    private
    procedure :: jacobian_block_multiply
    procedure :: jacobian_block_add
    generic, public :: operator(*) => jacobian_block_multiply
    generic, public :: operator(+) => jacobian_block_add
    procedure, public :: solve_for => jacobian_block_solve
  end type jacobian_block

  interface jacobian_block
    module procedure constructor
  end interface jacobian_block

contains

  function constructor(source_field, direction, extra_derivative, &
                       boundaries) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Build a block in a Jacobian matrix, with the form
    ! $$\frac{\partial F}{\partial x_i} + F\Delta_i,$$ where \(F\) is
    ! a scalar field and \(\Delta_i\) is the differentiation operator
    ! in the \(i\)-direction. Additionally, a further differentiation
    ! operator may be added to the right hand side of this matrix block.
    !
    class(scalar_field), intent(in)             :: source_field
      !! A scalar field (\(F\)) making up this block of the Jacobian
    integer, intent(in)                         :: direction
      !! The direction in which field derivatives are taken.
    integer, intent(in), optional               :: extra_derivative
      !! If present, specifies the direction of a differentiation
      !! operator to be added to the right hand side of this matrix
      !! block.
    procedure(jacobian_block_bounds),  optional :: boundaries
      !! A subroutine which specifies boundary values for the solution
      !! to Jacobian system. This is only needed if you will be using
      !! the `solve_for` method.
    type(jacobian_block)            :: this
      !! A new Jacobian block
    allocate(this%contents, source=source_field)
    allocate(this%derivative, mold=this%contents)
    this%direction = direction
    this%derivative = this%contents%d_dx(this%direction)
    if (present(extra_derivative)) this%extra_derivative = extra_derivative
    if (present(boundaries)) this%get_boundaries => boundaries
  end function constructor

  function jacobian_block_multiply(this, rhs) result(product)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Provides a matrix multiplication operator between a Jacobian
    ! block and a scalar field (which corresponds to a state vector).
    !
    class(jacobian_block), intent(in) :: this
    class(scalar_field), intent(in)   :: rhs
      !! A field corresponding to a state vector being multiplied by
      !! the Jacobian block.
    class(scalar_field), allocatable  :: product
    class(scalar_field), allocatable  :: tmp
    allocate(product, mold=this%contents)
    if (this%extra_derivative==no_extra_derivative) then
      if (this%has_increment) then
        product = (this%derivative + this%scalar_increment) * rhs &
                + this%contents * rhs%d_dx(this%direction)
      else
        product = this%derivative * rhs + this%contents * rhs%d_dx(this%direction)
      end if
    else
      allocate(tmp, mold=rhs)
      tmp = rhs%d_dx(this%extra_derivative)
      if (this%has_increment) then
        product = this%derivative * tmp + this%contents * tmp%d_dx(this%direction) &
                + this%scalar_increment*rhs
      else
        product = this%derivative * tmp + this%contents * tmp%d_dx(this%direction)
      end if
    end if
  end function jacobian_block_multiply

  function jacobian_block_add(this, rhs) result(sum)
    !* Author: Chris MacMackin 
    !  Date: December 2016
    !
    ! Provides a matrix multiplication operator between a Jacobian
    ! block and a scalar field (which corresponds to a state vector).
    !
    class(jacobian_block), intent(in) :: this
    real(r8), intent(in)              :: rhs
      !! A scalar which should be added to this block
    type(jacobian_block)              :: sum
    sum = this
    sum%scalar_increment = rhs
    sum%has_increment = .true.
  end function jacobian_block_add

  function jacobian_block_solve(this, rhs) result(solution)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Solves the linear(ised) system represented by this Jacobian
    ! block, for a given right hand side state vector (represented by
    ! a scalar field).
    !
    ! @Warning Currently this is only implemented for a 1-D field. The
    ! extra derivative must be in the direction along which the field
    ! varies.
    !
    class(jacobian_block), intent(inout) :: this
    class(scalar_field), intent(in)      :: rhs
      !! The right hand side of the linear(ised) system.
    class(scalar_field), allocatable     :: solution
    real(r8), dimension(:), allocatable  :: sol_vector

    character(len=52), parameter :: error_format = '("Tridiagonal matrix '// &
         'solver returned with flag ",i5)'
    character(len=77), parameter :: success_format = '("Tridiagonal matrix '// &
         'solver returned with estimated condition number ",es8.5)'
    integer, parameter                        :: error_len = 50, &
                                                 success_len = 75
    integer                                   :: flag, n
    real(r8)                                  :: forward_err, &
                                                 backward_err, &
                                                 condition_num
    character(len=1)                          :: factor
    character(len=:), allocatable             :: msg
    class(vector_field), allocatable          :: grid
    logical                                   :: use_cached
    class(scalar_field), allocatable, save    :: cached_field_type
    real(r8), dimension(:), allocatable, save :: cached_dx_c,   &
                                                 cached_dx2_uc, &
                                                 cached_dx2_ud, &
                                                 cached_dx2_dc
    real(r8), dimension(:), allocatable       :: cont, deriv, rhs_raw
    integer                                   :: i, pos
    
    allocate(sol_vector(rhs%raw_size()))
    ! Construct tridiagonal matrix for this operation, if necessary
    if (.not. allocated(this%pivots)) then
      call logger%debug('jacobian_block_solve','Constructing tridiagonal matrix.')
      factor = 'N'
      n = this%contents%raw_size()
      ! Try to use cached copy of inverse grid spacing, if available and suitable 
      if (allocated(cached_field_type)) then
        use_cached = same_type_as(this%contents, cached_field_type) .and. &
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
          deallocate(cached_dx2_uc)
          deallocate(cached_dx2_dc)
          deallocate(cached_dx2_ud)
        end if
        allocate(cached_field_type, mold=this%contents)
        allocate(cached_dx_c(n))
        allocate(cached_dx2_uc(n-1))
        allocate(cached_dx2_dc(n-1))
        allocate(cached_dx2_ud(n))
        call this%contents%allocate_vector_field(grid)
        grid = this%contents%grid_spacing()
        ! Use this array to temporarily hold the grid spacings
        cached_dx_c = grid%raw()
        ! Work out the upwinded grid spacings
        cached_dx2_uc(1) = cached_dx_c(1)
        do i = 2, n-1
          cached_dx2_uc(i) = 2*cached_dx_c(i) - cached_dx2_uc(i-1)
        end do
        ! Work out the inverse of the downwinded grid spacings
        cached_dx2_dc(1:n-2) = 1._r8/cached_dx2_uc(1:n-2)
        cached_dx2_dc(n-1) = 1._r8/cached_dx_c(n)
        ! Work out the -2 times product of the inverse up- and
        ! down-winded grid spacings
        cached_dx2_ud(1) = cached_dx2_uc(1)**(-2)
        cached_dx2_ud(2:n-1) = -2 * cached_dx2_dc(2:n-1)/cached_dx2_uc(1:n-2)
        cached_dx2_ud(n) = cached_dx2_dc(n-1)**2
        ! Work out the inverse of product of the up-winded and
        ! centred grid spacings
        cached_dx2_uc(1) = -cached_dx2_uc(1)**(-2)
        cached_dx2_uc(2:n-1) = 1._r8/(cached_dx2_uc(2:n-1) * cached_dx_c(2:n-1))
        ! Work out the inverse of the product of the down-winded and
        ! centred grid spacings
        cached_dx2_dc(1:n-2) = cached_dx2_dc(1:n-2) / cached_dx_c(2:n-1)
        cached_dx2_dc(n-1) = -cached_dx2_dc(n-1)**2
        ! Work out half the inverse of the grid spacing
        cached_dx_c(2:n-1) = 0.5_r8/cached_dx_c(2:n-1)
        cached_dx_c(1) = 1._r8/cached_dx_c(1)
        cached_dx_c(n) = 1._r8/cached_dx_c(n)
      end if
      cont = this%contents%raw()
      if (this%extra_derivative == no_extra_derivative) then
        ! Create the tridiagonal matrix when there is no additional
        ! differentiation operator on the RHS
        this%diagonal = this%derivative%raw()
        if (this%has_increment) then
          this%diagonal = this%diagonal + this%scalar_increment
        end if
        this%diagonal(1) = this%diagonal(1) - cont(1)*cached_dx_c(1)
        this%diagonal(n) = this%diagonal(n) + cont(n)*cached_dx_c(n)
        this%super_diagonal = cont(1:n-1) * cached_dx_c(1:n-1)
        this%sub_diagonal = -cont(2:n) * cached_dx_c(2:n)
      else if (this%extra_derivative == this%direction) then
        ! Create the tridiagonal matrix when the additional
        ! differentiation operator on the RHS operates in the same
        ! direction as the derivative of the contents.
        deriv = this%derivative%raw()
        allocate(this%diagonal(n))
        this%diagonal(2:n-1) = cont(2:n-1) * cached_dx2_ud(2:n-1)
        this%diagonal(1) = cont(1) * cached_dx2_ud(1) &
                         - deriv(1) * cached_dx_c(1)
        this%diagonal(n) = cont(n) * cached_dx2_ud(n) &
                         + deriv(n) * cached_dx_c(n)
        if (this%has_increment) then
          this%diagonal = this%diagonal + this%scalar_increment
        end if
        allocate(this%super_diagonal(n-1))
        this%super_diagonal(2:n-1) = cont(2:n-1) * cached_dx2_uc(2:n-1) &
                                   + deriv(2:n-1) * cached_dx_c(2:n-1)
        this%super_diagonal(1) = -2 * cont(1) * cached_dx2_uc(1) &
                               + deriv(1) * cached_dx_c(1)
        allocate(this%sub_diagonal(n-1))
        this%sub_diagonal(1:n-2) = cont(2:n-1) * cached_dx2_dc(1:n-2) &
                                 - deriv(2:n-1) * cached_dx_c(2:n-1)
        this%sub_diagonal(n-1) = -2 * cont(n) * cached_dx2_dc(n-1) &
                               - deriv(n) * cached_dx_c(n)
      else
        error stop('Currently no support for extra derivatives other than '// &
                   'in the 1-direction')
      end if
      ! Figure out the boundary conditions for this problem
      call this%get_boundaries(rhs, this%boundary_vals, this%boundary_locs)
      do i = 1, size(this%boundary_locs)
        pos = this%boundary_locs(i)
        this%diagonal(pos) = 1._r8
        if (pos < n) this%super_diagonal(pos) = 0._r8
        if (pos > 1) this%sub_diagonal(pos-1) = 0._r8
      end do
      ! Allocate the arrays used to hold the factorisation of the
      ! tridiagonal matrix
      allocate(this%l_multipliers(n-1))
      allocate(this%u_diagonal(n))
      allocate(this%u_superdiagonal1(n-1))
      allocate(this%u_superdiagonal2(n-2))
      allocate(this%pivots(n))
    else
      factor = 'F'
    end if
    rhs_raw = rhs%raw()
    do i = 1, size(this%boundary_locs)
      pos = this%boundary_locs(i)
      rhs_raw(pos) = this%boundary_vals(i)
    end do
    call la_gtsvx(this%sub_diagonal, this%diagonal, this%super_diagonal,      &
                  rhs_raw, sol_vector, this%l_multipliers, this%u_diagonal, &
                  this%u_superdiagonal1, this%u_superdiagonal2, this%pivots,  &
                  factor, 'N', forward_err, backward_err, condition_num,      &
                  flag)
    if (flag/=0) then
      allocate(character(len=error_len) :: msg)
      write(msg,error_format) flag
      call logger%error('jacobian_block_solve',msg)
    else
      allocate(character(len=success_len) :: msg)
      write(msg,success_format) condition_num
      call logger%debug('jacobian_block_solve',msg)
    end if
    allocate(solution, mold=rhs)
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
  end function jacobian_block_solve

  subroutine jacobian_block_bounds(rhs, boundary_values, boundary_locations)
    !* Author: Chris MacMackin
    !  Date: January 2016
    !
    ! A default implementation of the `get_boundaries` procedure
    ! pointer for the `jacobian_block` type. It corresponds to free
    ! boundaries.
    !
    class(scalar_field), intent(in)                  :: rhs
      !! The scalar field representing the vector being multiplied
      !! by the inverse Jacobian (i.e. the right hand side of the
      !! Jacobian system being solved).
    real(r8), dimension(:), allocatable, intent(out) :: boundary_values
      !! The values specifying the boundary conditions.
    integer, dimension(:), allocatable, intent(out)  :: boundary_locations
      !! The locations in the raw representation of `rhs` with which
      !! each of the elements of `boundary_values` is associated.
    allocate(boundary_values(0))
    allocate(boundary_locations(0))
  end subroutine jacobian_block_bounds
  
end module jacobian_block_mod
