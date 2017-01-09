!
!  preconditioner.f90
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

module preconditioner_mod
  !* Author: Christopher MacMackin
  !  Date: December 2016
  !  License: GPLv3
  !
  ! Provides a type for preconditioning fields in an iterative
  ! solver using Picard iteration.
  !
  use iso_fortran_env, only: r8 => real64
  use factual_mod, only: scalar_field, vector_field, maxval, abs
  use logger_mod, only: logger => master_logger
  use jacobian_block_mod, only: jacobian_block
  implicit none
  private

  type, public :: preconditioner
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Uses Picard iterations to apply the inverse Jacobian of a system
    ! to a vector, to low accuracy. Rather than directly computing the
    ! inverse Jacobian, it is more efficient to approximate it. If
    ! \(d\) is the vector being preconditioned, and \(z\) is the
    ! result of applying the preconditioner, then $$ z = J^{-1}d
    ! \Rightarrow Jz = d.$$ Thus, the preconditioner can be applied by
    ! approximately solving this system for \(z\). Linearising \(J\),
    ! this system can be solved efficiently using Picard iteration.
    !
    ! Say that the Jacobian is an \(n\times n\) system of blocks of
    ! the sort implemented in the [[jacobian_block]] type (each
    ! labeled as \(J_{j,k}\)) and that the vector being preconditioned
    ! constists of \(n\) scalar fields (\(d_j\)). Then \(m^{th}\)
    ! estimate of the solution for the \(j^{th}\) field in the
    ! preconditioned vector (\(z^m_j\)) is the solution to $$
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
    private
    real(r8) :: tolerance      = 1.e-3_r8
    integer  :: max_iterations = 20
  contains
    procedure :: apply => preconditioner_apply
  end type preconditioner

  interface preconditioner
    module procedure constructor
  end interface
  
contains

  function constructor(tolerance, max_iterations) result(this)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Create a preconditioner object with the desired tolerance and
    ! maximum number of iterations.
    !
    real(r8), optional, intent(in) :: tolerance
      !! The tolerance within which to apply the inverse Jacobian.
      !! Defaults to 0.001.
    integer, optional, intent(in)  :: max_iterations
      !! The maximum number of iterations to use when applying the
      !! preconditioner. Defaults to 20.
    type(preconditioner) :: this
    if (present(tolerance)) this%tolerance = tolerance
    if (present(max_iterations)) this%max_iterations = max_iterations
  end function constructor

  subroutine preconditioner_apply(this, jacobian, vector, estimate)
    !* Author: Chris MacMackin
    !  Date: December 2016
    !
    ! Use Picard iteration to approximately multiply the state vector
    ! by the inverse Jacobian. The details for this procedure are in
    ! the documentation of the [[preconditioner(type)]] type.
    !
    class(preconditioner), intent(in)                    :: this
    class(jacobian_block), dimension(:,:), intent(inout) :: jacobian
      !! An \(n\times n\) matrix approximating the Jacobian for which
      !! the preconditioner is used.
    class(scalar_field), dimension(:), intent(in)        :: vector
      !! A vector of size \(n\) which is to be preconditioned.
    class(scalar_field), dimension(:), intent(inout)     :: estimate
      !! On entry, an initial guess for the preconditioned vector. On
      !! exit, the iteratively determined value of the preconditioned
      !! vector.

    character(len=76), parameter :: success_format = '("Picard solver '// &
         'converged with error of ",es8.5," after ",i3," iterations.")'
    character(len=76), parameter :: failure_format = '("Picard solver '// &
         'reached maximum (",i3,") iterations, with error ",es8.5,".")'
    integer, parameter                             :: msg_len = 68
    class(scalar_field), dimension(:), allocatable :: prev_estimate
    integer                                        :: i, j, k, n
    real(r8)                                       :: max_err, old_max_err
    class(scalar_field), allocatable               :: tmp_field
    logical                                        :: first
    character(len=msg_len)                         :: msg
    
    call logger%debug('preconditioner_apply','Entering function `precondition`.')
    n = size(vector)
    allocate(prev_estimate(n), mold=estimate)
    allocate(tmp_field, mold=vector(1))
#ifdef DEBUG
    if (size(jacobian,1) /= size(jacobian,2)) then
      error stop('Jacobian is not a square matrix.')
    end if
    if (size(jacobian,1) /= size(vector)) then
      error stop('Vector is of different size than Jacobian.')
    end if
    if (size(estimate) /= size(vector)) then
      error stop('Estimate is of different size than vector.')
    end if
#endif
    ! Until reached maximum number of iterations...
    do i = 1, this%max_iterations
      ! For each row of the Jacobian, solve for the diagonal block,
      ! with off-diagonal blocks applied to the previous guess and
      ! subtracted from the right-hand-side.
      max_err = 0._r8
      prev_estimate = estimate
      print*,'----------------------------------------------'
      do j = 1, n
        first = .true.
        do k = 1, n
          if (k == j) cycle
          if (first) then
            first = .false.
            tmp_field = jacobian(j,k) * estimate(k)
          else
            tmp_field = tmp_field + jacobian(j,k) * estimate(k)
          end if
        end do
        estimate(j) = jacobian(j,j)%solve_for(vector(j) - tmp_field)
        max_err = max(max_err, &
                      maxval(abs( (estimate(j)-prev_estimate(j))/(prev_estimate(j)+1e-10_r8) )))
        print*,i,j,estimate(j)%raw()
      end do
      print*,i,max_err
      ! If difference between result and previous guess is less than
      ! the tolerance, stop iterations
      if (max_err < this%tolerance) then
        write(msg,success_format) max_err, i
        call logger%debug('preconditioner_apply',msg)
        call logger%debug('preconditioner_apply','Exiting function `precondition`.')
        return
      end if
      print*,i,old_max_err,max_err,old_max_err<=max_err
      if (i > 1 .and. old_max_err <= max_err) then
        call logger%error('preconditioner_apply', &
                         'Iterations diverging. Exiting and returning previous iterate.')
        call logger%debug('preconditioner_apply','Exiting function `precondition`.')
        estimate = prev_estimate
        return
      end if
      old_max_err = max_err
    end do
    write(msg,failure_format) this%max_iterations, max_err
    call logger%error('preconditioner_apply',msg)
    call logger%debug('preconditioner_apply','Exiting function `precondition`.')
  end subroutine preconditioner_apply
  
end module preconditioner_mod
