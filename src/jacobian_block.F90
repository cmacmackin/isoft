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
  use boundary_types_mod, only: dirichlet, neumann
  use f95_lapack, only: la_gtsvx
  use penf, only: str
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
    ! Boundary conditions are a complicated issue. When constructing
    ! an instance of this type, locations in the raw field
    ! representation can be specified for which boundary conditions
    ! are given. By default these are all Dirichlet conditions, but
    ! this behaviour may be overridden to be Neumann conditions. Note
    ! that Neumann conditions result in approximately an order of
    ! magnitude lower accuracy. When using the Jacobian block for
    ! multiplation, the type of the boundary does not matter; the
    ! value at the boundary in the result will be set to 0 or to
    ! values produced by an optional user-provided function. When
    ! using the `solve_for` method, the tridiagonal matrix will be
    ! altered in the specified locations to reflect the correct
    ! boundary conditions. This effectively mimics the sorts of
    ! Jacobian rows typically produced by boundary conditions.
    !
    private
    integer                             :: direction
      !! The direction in which any derivatives are taken.
    integer                             :: extra_derivative = no_extra_derivative
      !! The direction in which to apply a differentiation on the
      !! right-hand-side of the Jacobian block operator. Defaults
      !! none.
    class(scalar_field), allocatable    :: contents
      !! The value, \(A\), to which the Jacobian block operation is
      !! beiing applied.
    real(r8)                            :: coef = 1._r8
      !! Optional coefficient by which the the \(\partial F/\partial x
      !! \) term in the operator will be multiplied.
    procedure(jacobian_block_bounds), pointer, nopass  :: &
                                           get_boundaries
      !! A subroutine which determines the expected boundary conditions
      !! (and their location in the raw array) for the solution of the
      !! Jacobian block.
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
    integer, dimension(:), allocatable  :: boundary_types
      !! The types of boundary conditions, specified using the
      !! parameters found in [[boundary_types_mod]].
    real(r8)                            :: real_increment
      !! A scalar value which is to be added to this Jacobian block
      !! (i.e. to the diagonal).
    class(scalar_field), allocatable    :: field_increment
      !! A scalar field which is to be added to this Jacobian block
      !! (i.e. to the diagonal).
    type(jacobian_block), pointer       :: block_increment => null()
      !! Another Jacobian block which is to be added to this one
    integer                             :: has_increment = 0
      !! Indicates whether or not there has been an increment added to
      !! this block. If not, then 0. If a scalar real value has been
      !! added, then 1. If a scalar value has been added, then 2.
  contains
    private
    procedure :: jacobian_block_multiply
    procedure :: jacobian_block_add_real
    procedure :: jacobian_block_add_field
    procedure :: jacobian_block_add_block
    procedure :: jacobian_block_assign
    procedure :: get_tridiag => jacobian_block_get_tridiag
    generic, public :: operator(*) => jacobian_block_multiply
    generic, public :: operator(+) => jacobian_block_add_real,  &
                                      jacobian_block_add_field, &
                                      jacobian_block_add_block
    generic, public :: assignment(=) => jacobian_block_assign
    procedure, public :: solve_for => jacobian_block_solve
  end type jacobian_block

  interface jacobian_block
    module procedure constructor
  end interface jacobian_block

contains

  function constructor(source_field, direction, extra_derivative, &
                       boundary_locs, boundary_types,             &
                       boundary_operations, coef) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Build a block in a Jacobian matrix, with the form
    ! $$\frac{\partial F}{\partial x_i} + F\Delta_i,$$ where \(F\) is
    ! a scalar field and \(\Delta_i\) is the differentiation operator
    ! in the \(i\)-direction. Additionally, a further differentiation
    ! operator may be added to the right hand side of this matrix
    ! block.  Optional arguments allow for handling of boundary
    ! conditions. See the end of the documentation of the
    ! [[jacobian_block(type)]] type for a description of how boundary
    ! conditions are treated.
    !
    class(scalar_field), intent(in)             :: source_field
      !! A scalar field (\(F\)) making up this block of the Jacobian
    integer, intent(in)                         :: direction
      !! The direction in which field derivatives are taken.
    integer, intent(in), optional               :: extra_derivative
      !! If present, specifies the direction of a differentiation
      !! operator to be added to the right hand side of this matrix
      !! block.
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
    procedure(jacobian_block_bounds),  optional :: boundary_operations
      !! A function specifying the values to place at the boundaries
      !! of the result when using the Jacobian block for
      !! multiplication. By default, all boundaries are set to 0. The
      !! order in which the resulting values are stored should match
      !! that of `boundary_locs`.
    real(r8), optional, intent(in)              :: coef
      !! An optional coefficient by which the the \(\partial
      !! F/\partial x \) term in the operator will be
      !! multipled. Default value is 1.
    type(jacobian_block)                        :: this
      !! A new Jacobian block
    call source_field%guard_temp()
    allocate(this%contents, mold=source_field)
    allocate(this%derivative, mold=source_field)
    this%contents = source_field
    this%direction = direction
    this%derivative = this%contents%d_dx(this%direction)
    if (present(extra_derivative)) this%extra_derivative = extra_derivative
    if (present(coef)) this%coef = coef
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
    if (present(boundary_operations)) then
      this%get_boundaries => boundary_operations
    else
      this%get_boundaries => jacobian_block_bounds
    end if
    call source_field%clean_temp()
#ifdef DEBUG
    call logger%debug('jacobian_block','Instantiated a Jacobian block object.')
#endif
  end function constructor

  recursive function jacobian_block_multiply(this, rhs) result(product)
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
    class(scalar_field), pointer :: product
    class(scalar_field), pointer :: tmp
    real(r8), dimension(:), allocatable :: bounds
    integer :: i
    call rhs%guard_temp()
    call this%contents%allocate_scalar_field(product)
    call product%unset_temp()
    if (this%extra_derivative==no_extra_derivative) then
      select case(this%has_increment)
      case(0)
        product = this%coef*this%derivative * rhs + this%contents * rhs%d_dx(this%direction)
      case(1)
        product = (this%coef*this%derivative + this%real_increment) * rhs &
                + this%contents * rhs%d_dx(this%direction)
      case(2)
        product = (this%coef*this%derivative + this%field_increment) * rhs &
                + this%contents * rhs%d_dx(this%direction)
      case default
        error stop ('Invalid increment has been added.')
      end select
    else
      call rhs%allocate_scalar_field(tmp)
      call tmp%guard_temp()
      tmp = rhs%d_dx(this%extra_derivative)
      select case(this%has_increment)
      case(0)
        product = this%coef*this%derivative * tmp + this%contents * tmp%d_dx(this%direction)
      case(1)
        product = this%coef*this%derivative * tmp + this%contents * tmp%d_dx(this%direction) &
                + this%real_increment*rhs
      case(2)
        product = this%coef*this%derivative * tmp + this%contents * tmp%d_dx(this%direction) &
                + this%field_increment*rhs
      case default
        error stop ('Invalid increment has been added.')
      end select
      call tmp%clean_temp()
    end if
    if (associated(this%block_increment)) then
      product = product + this%block_increment*rhs
    end if
    call this%get_boundaries(this%contents,this%derivative,rhs,      &
                            this%boundary_locs,this%boundary_types, &
                            bounds)
    do i = 1, size(this%boundary_locs)
      call product%set_element(this%boundary_locs(i), bounds(i))
    end do
    call rhs%clean_temp()
#ifdef DEBUG
    call logger%debug('jacobian_block%multiply','Multiplied vector by '// &
                      'Jacobian block.')
#endif
    call product%set_temp()
  end function jacobian_block_multiply

  function jacobian_block_add_real(this, rhs) result(sum)
    !* Author: Chris MacMackin 
    !  Date: December 2016
    !
    ! Produces a Jacobian block which has been offset by some constant
    ! increment.
    !
    ! @Warning This operation will overwrite any previous sums which
    ! have been performed to produce `this`.
    !
    class(jacobian_block), intent(in) :: this
    real(r8), intent(in)              :: rhs
      !! A scalar which should be added to this block
    type(jacobian_block)              :: sum
    sum = this
    sum%real_increment = rhs
    sum%has_increment = sum%has_increment + 1
#ifdef DEBUG
    call logger%debug('jacobian_block%add','Added real to a Jacobian block.')
#endif
  end function jacobian_block_add_real

  function jacobian_block_add_field(this, rhs) result(sum)
    !* Author: Chris MacMackin 
    !  Date: May 2017
    !
    ! Produces a Jacobian block which has been offset by a scalar
    ! field.
    !
    ! @Warning This operation will overwrite any previous sums which
    ! have been performed to produce `this`.
    !
    class(jacobian_block), intent(in) :: this
    class(scalar_field), intent(in)   :: rhs
      !! A scalar which should be added to this block
    type(jacobian_block)              :: sum
    sum = this
    allocate(sum%field_increment, mold=rhs)
    sum%field_increment = rhs
    sum%has_increment = sum%has_increment + 2
#ifdef DEBUG
    call logger%debug('jacobian_block%add','Added field to a Jacobian block.')
#endif
  end function jacobian_block_add_field

  function jacobian_block_add_block(this, rhs) result(sum)
    !* Author: Chris MacMackin 
    !  Date: May 2017
    !
    ! Produces a Jacobian block which is the sum of two existing
    ! blocks. Boundary conditions are set by the first operand
    ! (`this`).
    !
    class(jacobian_block), intent(in)         :: this
    class(jacobian_block), target, intent(in) :: rhs
      !! A second block which should be added to this block
    type(jacobian_block)                      :: sum
    sum = this
    sum%block_increment => rhs
#ifdef DEBUG
    call logger%debug('jacobian_block%add','Added field to a Jacobian block.')
#endif
  end function jacobian_block_add_block

  subroutine jacobian_block_assign(this, rhs)
    !* Author: Chris MacMackin 
    !  Date: December 2016
    !
    ! Copies the contents of the `rhs` Jacobian block into this
    ! one. It will safely deallocate any data necessary.
    !
    class(jacobian_block), intent(out) :: this
    type(jacobian_block), intent(in)   :: rhs
    this%direction = rhs%direction
    this%extra_derivative = rhs%extra_derivative
    allocate(this%contents, mold=rhs%contents)
    allocate(this%derivative, mold=rhs%derivative)
    this%contents = rhs%contents
    this%coef = rhs%coef
    this%derivative = rhs%derivative
    this%get_boundaries => rhs%get_boundaries
    if (allocated(rhs%diagonal)) then
      this%diagonal = rhs%diagonal
      this%super_diagonal = rhs%super_diagonal
      this%sub_diagonal = rhs%sub_diagonal
      this%l_multipliers = rhs%l_multipliers
      this%u_diagonal = rhs%u_diagonal
      this%u_superdiagonal1 = rhs%u_superdiagonal1
      this%u_superdiagonal2 = rhs%u_superdiagonal2
      this%pivots = rhs%pivots
    end if
    if (allocated(rhs%boundary_locs)) this%boundary_locs = rhs%boundary_locs
    if (allocated(rhs%boundary_types)) this%boundary_types = rhs%boundary_types
    this%real_increment = rhs%real_increment
    if (allocated(rhs%field_increment)) &
         allocate(this%field_increment, source=rhs%field_increment)
    if (associated(rhs%block_increment)) this%block_increment => rhs%block_increment
    this%has_increment = rhs%has_increment
  end subroutine jacobian_block_assign

  recursive subroutine jacobian_block_get_tridiag(this, diagonal, subdiagonal, &
                                                  superdiagonal)
    !* Author: Chris MacMackin
    !  Date: May 2017
    !
    ! Computes the tridiagonal matrix used to solve for this Jacobian block.
    !
    class(jacobian_block), intent(in)                :: this
    real(r8), dimension(:), allocatable, intent(out) :: diagonal
    real(r8), dimension(:), allocatable, intent(out) :: subdiagonal
    real(r8), dimension(:), allocatable, intent(out) :: superdiagonal

    integer                                   :: n, i, pos
    class(vector_field), pointer              :: grid
    logical                                   :: use_cached
    class(scalar_field), allocatable, save    :: cached_field_type
    real(r8), dimension(:), allocatable, save :: cached_dx_c,   &
                                                 cached_dx2_uc, &
                                                 cached_dx2_ud, &
                                                 cached_dx2_dc
    real(r8), save                            :: cached_upper_bound, &
                                                 cached_lower_bound
    real(r8), dimension(:), allocatable       :: cont, deriv, diag, &
                                                 subdiag, supdiag
    real(r8), dimension(:,:), allocatable     :: domain

    n = this%contents%raw_size()
    domain = this%contents%domain()
    ! Try to use cached copy of inverse grid spacing, if available and suitable 
    if (allocated(cached_field_type)) then
      use_cached = same_type_as(this%contents, cached_field_type) .and. &
                   abs(domain(1,1) - cached_lower_bound) < 1e-10 .and.  &
                   abs(domain(1,2) - cached_upper_bound) < 1e-10 .and.  &
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
#ifdef DEBUG
      call logger%debug('jacobian_block%get_tridiag','Calculating and caching '// &
                       'grid spacings.')
#endif
      cached_lower_bound = domain(1,1)
      cached_upper_bound = domain(1,2)
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
      diagonal = this%coef*this%derivative%raw()
      select case(this%has_increment)
      case(0)
        continue
      case(1)
        diagonal = diagonal + this%real_increment
      case(2)
        diagonal = diagonal + this%field_increment%raw()
      case default
        error stop ('Invalid increment has been added.')
      end select
      diagonal(1) = diagonal(1) - cont(1)*cached_dx_c(1)
      diagonal(n) = diagonal(n) + cont(n)*cached_dx_c(n)
      superdiagonal = cont(1:n-1) * cached_dx_c(1:n-1)
      subdiagonal = -cont(2:n) * cached_dx_c(2:n)
    else if (this%extra_derivative == this%direction) then
      ! Create the tridiagonal matrix when the additional
      ! differentiation operator on the RHS operates in the same
      ! direction as the derivative of the contents.
      deriv = this%coef*this%derivative%raw()
      allocate(diagonal(n))
      diagonal(2:n-1) = cont(2:n-1) * cached_dx2_ud(2:n-1)
      diagonal(1) = cont(1) * cached_dx2_ud(1) &
                       - deriv(1) * cached_dx_c(1)
      diagonal(n) = cont(n) * cached_dx2_ud(n) &
                       + deriv(n) * cached_dx_c(n)
      !print*,diagonal(10)
      select case(this%has_increment)
      case(0)
        continue
      case(1)
        diagonal = diagonal + this%real_increment
      case(2)
        diagonal = diagonal + this%field_increment%raw()
      case default
        error stop ('Invalid increment has been added.')
      end select
      allocate(superdiagonal(n-1))
      superdiagonal(2:n-1) = cont(2:n-1) * cached_dx2_uc(2:n-1) &
                           + deriv(2:n-1) * cached_dx_c(2:n-1)
      superdiagonal(1) = -2 * cont(1) * cached_dx2_uc(1) &
                       + deriv(1) * cached_dx_c(1)
      allocate(subdiagonal(n-1))
      subdiagonal(1:n-2) = cont(2:n-1) * cached_dx2_dc(1:n-2) &
                         - deriv(2:n-1) * cached_dx_c(2:n-1)
      subdiagonal(n-1) = -2 * cont(n) * cached_dx2_dc(n-1) &
                       - deriv(n) * cached_dx_c(n)
    else
      call logger%fatal('jacobian_block%solve_for', 'Currently no support '// &
                        'for extra derivatives other than in the 1-direction')
      error stop
    end if
    if (associated(this%block_increment)) then
      call this%block_increment%get_tridiag(diag, subdiag, supdiag)
      diagonal = diagonal + diag
      subdiagonal = subdiagonal + subdiag
      superdiagonal = superdiagonal + supdiag
    end if
    ! Set the boundary conditions for this problem
    do i = 1, size(this%boundary_locs)
      pos = this%boundary_locs(i)
      if (this%boundary_types(i) == dirichlet) then
        diagonal(pos) = 1._r8
        if (pos < n) superdiagonal(pos) = 0._r8
        if (pos > 1) subdiagonal(pos-1) = 0._r8
      else if (this%boundary_types(i) == neumann) then
        if (pos == n) then
          diagonal(n) = cached_dx_c(n)
          subdiagonal(n-1) = -cached_dx_c(n)
        else if (pos == 1) then
          diagonal(1) = -cached_dx_c(1)
          superdiagonal(1) = cached_dx_c(1)
        else
          superdiagonal(pos) = cached_dx_c(pos)
          subdiagonal(pos-1) = cached_dx_c(pos)
          diagonal(pos) = 0._r8
        end if
      else
        call logger%fatal('jacobian_block%solve_for','Boundary condition of '// &
                          'type other than Dirichlet or Neumann encountered.')
        error stop
      end if
    end do
#ifdef DEBUG
    call logger%debug('jacobian_block%get_tridiag','Constructed tridiagonal matrix.')
#endif
  end subroutine jacobian_block_get_tridiag

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
    class(scalar_field), pointer         :: solution
    real(r8), dimension(:), allocatable  :: sol_vector

    integer                                   :: flag, n
    real(r8)                                  :: forward_err, &
                                                 backward_err, &
                                                 condition_num
    real(r8), dimension(:), allocatable       :: diag, subdiag, supdiag
    character(len=1)                          :: factor
    character(len=:), allocatable             :: msg

    call rhs%guard_temp()
    allocate(sol_vector(rhs%raw_size()))
    ! Construct tridiagonal matrix for this operation, if necessary
    if (.not. allocated(this%pivots)) then
      factor = 'N'
      ! Allocate the arrays used to hold the factorisation of the
      ! tridiagonal matrix
      call this%get_tridiag(diag, subdiag, supdiag)
      call move_alloc(diag, this%diagonal)
      call move_alloc(subdiag, this%sub_diagonal)
      call move_alloc(supdiag, this%super_diagonal)
      n = this%contents%raw_size()
      allocate(this%l_multipliers(n-1))
      allocate(this%u_diagonal(n))
      allocate(this%u_superdiagonal1(n-1))
      allocate(this%u_superdiagonal2(n-2))
      allocate(this%pivots(n))
    else
      factor = 'F'
    end if
!    print*,rhs%raw()
    call la_gtsvx(this%sub_diagonal, this%diagonal, this%super_diagonal,      &
                  rhs%raw(), sol_vector, this%l_multipliers, this%u_diagonal, &
                  this%u_superdiagonal1, this%u_superdiagonal2, this%pivots,  &
                  factor, 'N', forward_err, backward_err, condition_num,      &
                  flag)
    if (flag/=0) then
      msg = 'Tridiagonal matrix solver returned with flag '//trim(str(flag))
!      print*,this%sub_diagonal
!      print*,
!      print*,this%diagonal
!      print*,
!      print*,this%super_diagonal
      CALL backtrace
      call logger%error('jacobian_block%solve_for',msg)
      !print*,condition_num
!      stop
#ifdef DEBUG
    else
      msg = 'Tridiagonal matrix solver returned with estimated condition '// &
            'number '//trim(str(condition_num))
      call logger%debug('jacobian_block%solve_for',msg)
#endif
    end if
    call rhs%allocate_scalar_field(solution)
    call solution%unset_temp()
    call solution%assign_meta_data(rhs)
    call solution%set_from_raw(sol_vector)
    call rhs%clean_temp(); call solution%set_temp()
  end function jacobian_block_solve

  subroutine jacobian_block_bounds(contents, derivative, rhs,   &
                                 boundary_locs, boundary_types, &
                                 boundary_values)
    !* Author: Chris MacMackin
    !  Date: January 2016
    !
    ! A default implementation of the `get_boundaries` procedure
    ! pointer for the `jacobian_block` type. It corresponds to setting
    ! all boundaries to 0 when multiplying a field by the Jacobian
    ! block.
    !
    class(scalar_field), intent(in)                 :: contents
      !! The field used to construct the Jacobian block
    class(scalar_field), intent(in)                 :: derivative
      !! The first spatial derivative of the field used to construct
      !! the Jacobian block, in the direction specified
    class(scalar_field), intent(in)                 :: rhs
      !! The scalar field representing the vector being multiplied
      !! by Jacobian
    integer, dimension(:), allocatable, intent(in)  :: boundary_locs
      !! The locations in the raw representation of `rhs` containing
      !! the boundaries.
    integer, dimension(:), allocatable, intent(in)  :: boundary_types
      !! Integers specifying the type of boundary condition. The type
      !! of boundary condition corresponding to a given integer is
      !! specified in [[boundary_types_mod]]. Only Dirichlet and
      !! Neumann conditions are supported. The storage order must
      !! correspond to that of `boundary_locs`.
    real(r8), dimension(:), allocatable, intent(out) :: boundary_values
      !! The values to go at the boundaries when multiplying a field
      !! by the Jacobian block. The storage order must be the same as
      !! for `boundary_locs`.
    allocate(boundary_values(size(boundary_locs)))
    boundary_values = 0._r8
  end subroutine jacobian_block_bounds
  
end module jacobian_block_mod
